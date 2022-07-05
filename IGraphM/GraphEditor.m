(* ::Package:: *)

(* Mathematica Package *)
(* Created by Mathematica plugin for IntelliJ IDEA *)

(* :Author: Kuba Podkalicki  *)
(* :Date: 2020-09-30 *)

Package["IGraphM`"]


PackageExport["IGGraphEditor"]

(* Dev notes
 - foo[ Dynamic @ state_ ] is no better than HoldFirst etc. I just like this convention for GUI 
 - state should not be modified during continuous actions
*)

IGGraphEditor::usage        =
    "IGGraphEditor[] typesets to an interactive graph editor.\n" <>
    "IGGraphEditor[graph] uses the given graph as the starting point.";
IGGraphEditor::multiEdge    = "Multi-edges are not supported yet. Only directed pairs are {1->2, 2->1}.";
IGGraphEditor::unknownState = "Corrupted editor state."
IGGraphEditor::oldVer       = "You need to update IGraph/M to continue working with the data stored here."

$narrowAspectRatioLimit = 10 (* plot range will be adjusted if the initial calculated as is above this limit *)
$vertexEdgeThickness = 0.5
$hoverVertexEdgeThickness = 1
$edgeThickness = 1
$activeEdgeThickness = 3
$potentialEdgeColor = Purple
$gridLinesCount = 20.


(* ::Subsection:: *)
(*IGGraphEditor*)


IGGraphEditor // Options = {
  "KeepVertexCoordinates" -> True (* bool *)
, "CreateVertexSelects"   -> True (* bool *)
, "SnapToGrid"            -> False (* bool *)
, "IndexGraph"            -> False (* bool *)
, "PerformanceLimit"      -> 450   (* _Integer *)
, VertexLabels            -> False (* | "Name" *)
, VertexSize              -> Small (* Tiny | Small | Medium | Large | ratioToDiagonal_?NumericQ*)
, DirectedEdges           -> False (* bool *)
, ImageSize               -> Automatic
}


SyntaxInformation[IGGraphEditor] = {"ArgumentsPattern" -> {_., OptionsPattern[]}};

IGGraphEditor /:
  MakeBoxes[IGGraphEditor[graph:(_|PatternSequence[]), opt:OptionsPattern[]], StandardForm] :=
  ToBoxes @ iGraphEditor[graph, opt]


iGraphEditor // Options = Options @ IGGraphEditor;


supportedGraphQ = ! MixedGraphQ[#] && SimpleGraphQ[#]&


(* ::Subsection::Closed:: *)
(*iGraphEditor*)

editorFailure[ msg_String ] := Failure[
  "GraphEditor"
, <|"Message" -> "\[WarningSign] "<> msg |>
]

iGraphEditor[graph:((_? supportedGraphQ)|PatternSequence[]), opt:OptionsPattern[]] := 
With[
  {
    packageFailure = editorFailure["Failed to load IGraphM`, make sure it is installed."]
  , performanceFailure = editorFailure["Too many vertices and edges. Change 'PerformanceLimit' if you want to try anyway."]
  },
Interpretation[
  { 
    state = GraphToEditorState[graph, opt]
  , error = False 
  , refresh
  }
, refresh[] := Module[{temp}, Catch[
    If[
      Needs@"IGraphM`" === $Failed
    , Throw[ error = packageFailure]
    ]
  
  ; temp = geStateVersionCheck @ state
  ; If[ AssociationQ @ temp
    , state = temp
    , Throw[ error = temp]
    ] 

  ; If[ 
      state["config", "vCounter"] + state["config", "eCounter"] > state["config", "PerformanceLimit"]
    , Throw[ error =  performanceFailure]  
    ]   
    
  ; error = False        
  ; geAction["UpdateEdgesShapes", Hold @ state]  
    
  ]]

; 
  PaneSelector[
  { 
    True -> Button[Dynamic @ error,  refresh[], BaseStyle->15  ]
  , False  -> Panel[
      Dynamic[Refresh[iGraphEditorPanel[Dynamic@state], None]]
    , FrameMargins->0, BaseStyle->CacheGraphics->False
    ]
  }  
  , Dynamic[ MatchQ[_Failure] @ error ]
  , ImageSize->Automatic
  
  ]
  
, GraphFromEditorState @ state
  
, Initialization :> refresh[]
]]


iGraphEditor[_Graph, OptionsPattern[]] := Failure["GraphEditor", <|"Message" -> "Currently a Graph argument needs to be Simple And Not Mixed."|>]


iGraphEditor[___] := Failure["GraphEditor", <|"Message" -> "Unknown input"|>]


iGraphEditorPanel[Dynamic@state_] := EventHandler[
    geGraphics @ Dynamic @ state
  , "MouseClicked" :> (
      geAction["MouseClicked", Dynamic @ state, CurrentValue[{"MousePosition", "Graphics"}]]
    )
  , PassEventsDown->True
  ]


(* ::Subsection:: *)
(*graphics*)


geGraphics[Dynamic @ state_ ] := DynamicModule[{range}
, DynamicWrapper[
    Graphics[
      DynamicNamespace @ {
        geHighlightsPrimitives @ Dynamic @ state    
      , Gray    
      , Dynamic @ geEdges @ Dynamic @ state      
      , Dynamic @ geVertices @ Dynamic @ state    
      }
    , PlotRange -> Dynamic @ range
    , ImageSize -> state["config", "ImageSize"]
    ]
  , range = state["config", "range"]
  ]
] (* PlotRange->Dynamic@state["config... updates at any unrelated event, vertex dragging included,
     this DynamicModule @ DynamicWrapper is here to address a bug. Why does it help? Because WRI.
   *)


(* ::Subsubsection:: *)
(*vertex*)


geVertices[Dynamic @ state_] := Table[
  geVertexShapeFunction[Dynamic@state, state["vertex", id] ]
, {id, Keys @ state["vertex"] }
]


geVertexShapeFunction[Dynamic @ state_, v_Association] :=
With[ {
  nef = $vertexEdgeThickness
, aef = $hoverVertexEdgeThickness  
, step = state["config", "snapStep"]
},
DynamicModule[
  {x = v@"pos"},
Module[
  {mouseDragged, graphics}

, graphics = { {
    EdgeForm @ AbsoluteThickness @  Dynamic[ FEPrivate`If[  FrontEnd`CurrentValue["MouseOver"], aef, nef ] ]
  , DynamicName[ 
      Disk[Dynamic[x], state["config", "realVertexSize"]]
    , v["id"]
    ]
  }
  , If[ 
      state["config", "VertexLabels"] === "Name"
    , Inset[v["name"], Offset[ {12, 12}, DynamicLocation[v["id"]]] ]
    , Nothing
    ]
  }
  
; mouseDragged = If[ state["config", "snap"]
    , "MouseDragged" :> (x = Round[CurrentValue[{"MousePosition", "Graphics"}], step])
    , "MouseDragged" :> FEPrivate`Set[x , FrontEnd`CurrentValue[{"MousePosition", "Graphics"}] ]
  ]

; EventHandler[
    graphics,
    { "MouseClicked" :> geAction["VertexClicked", Dynamic @ state, v]
    , "MouseUp"      :> geAction["UpdateVertexPosition", Dynamic @ state, v["id"], x] 
    , mouseDragged
    },
    PassEventsUp -> False
  ]
      (*I'd prefer clicked to be Queued but if I put it in an inner queued EventHandler then I can't block MouseUp from fireing*)
]]]



(* ::Subsubsection::Closed:: *)
(*edges*)


geEdges[Dynamic @ state_] := Table[
  geEdgeShapeFunction[ Dynamic@state,  state["edge", id]], {id, Keys @ state["edge"]}
]  


geEdgeShapeFunction[Dynamic @ state_, e_Association] := 
With[
  { 
    nef = $edgeThickness,
    aef = $activeEdgeThickness
  },
  EventHandler[
      { AbsoluteThickness @  Dynamic[ FEPrivate`If[  FrontEnd`CurrentValue["MouseOver"], aef, nef ] ]
      , edgeToPrimitive @ e        
      }
  , { "MouseClicked" :> (geAction["EdgeClicked", Dynamic @ state, e]) }
  , PassEventsUp -> False (* edgeclicked should not be followed by outer mouseclicked*)
  ]
]


edgeToPrimitive[e_] := If[
  e["shape"] =!= Automatic
, e["shape"]
, edgeTypeToPrimitive[e["type"]
][
    DynamicLocation[e["v1"], Automatic],
    DynamicLocation[e["v2"], Automatic]
  ]
]  


edgeTypeToPrimitive["v1"->"v2"]=Arrow[{#,#2}]&
edgeTypeToPrimitive["v2"->"v1"]=Arrow[{#2,#}]&
edgeTypeToPrimitive[UndirectedEdge["v1", "v2"]]=Line[{#,#2}]&


(* ::Subsubsection::Closed:: *)
(*selected*)


geHighlightsPrimitives[Dynamic @ state_] := With[{ selV := state["selectedVertex"] }
, { $potentialEdgeColor,
    Dynamic[
      If[ 
        stateHasSelectedVertex @ state
      , { 
          (*Disk[DynamicLocation[selV], 1.05 state["config", "realVertexSize"]]
        , *)Dashed, Line[{
            DynamicLocation[selV]
          , FrontEnd`MousePosition["Graphics",DynamicLocation[selV]]
          (*TODO: with multiple editors with selected nodes the line is shown of each one of them*)
          }]
        }
      , {}  
      ]
    ]
}    
]


(* ::Subsection::Closed:: *)
(*state version*)


$stateVersion = 1; 
(*It does not change unless the structure of state association is changed. 
  That implies new geActions etc won't be compatible with old state *)


(*ONCE $stateVersion > 1 *)
(*
  geStateVersionCheck[ state : KeyValuePattern[{"version" \[Rule] oldVersion}] ] := (*oldVersion to oldVersion+1 transition rules*)
*)


geStateVersionCheck[ state : KeyValuePattern[{"version" -> $stateVersion}] ] := state


geStateVersionCheck[ state : KeyValuePattern[{"version" -> v_}] ] :=
  Failure["GraphEditor", <|"MessageTemplate" -> IGGraphEditor::oldVer|>]


geStateVersionCheck[___] :=
  Failure["GraphEditor", <|"MessageTemplate" -> IGGraphEditor::unknownState|>]


(* ::Subsection::Closed:: *)
(*from state*)


GraphFromEditorState[state_Association] := GraphFromEditorState[state, Lookup[state, "version", 1]];


GraphFromEditorState[state_, 1] := Module[{v, e, pos, graph}
, v = stateVertexList @ state

; e = stateEdgeList @ state

; pos = If[
    state["config", "KeepVertexCoordinates"]  
  , stateGraphEmbedding @ state 
  , Automatic
  ]
  
; graph = Graph[ v, e, VertexCoordinates -> pos]  

; If[
    TrueQ @ state["config", "IndexGraph"]
  , graph = IndexGraph @ graph
  ] 

; graph  
]


(* ::Subsection:: *)
(*to state*)


optionsToConfig // Options = Options @ IGGraphEditor;

optionsToConfig[OptionsPattern[]] := Association[
  ToString[#] -> OptionValue[#] & /@ Keys @ Options[IGGraphEditor]
]


GraphToEditorState[ opt:OptionsPattern[] ] := Module[{state}
, state = <|
      "vertex"         -> <||>
    , "edge"           -> <||>
    , "selectedVertex" -> Null
    , "version"        -> $stateVersion
    , "config"         -> <|
        optionsToConfig[opt]
      , "vCounter"->0
      , "eCounter" ->0
      , "range" -> {{-1, 1}, {-1, 1}}
      , "inRangeQ" -> RegionMember[ Rectangle[{-1,-1}, {1, 1}]  ]
      |>
    |>


; state
]   


GraphToEditorState[g_Graph ? supportedGraphQ, opt:OptionsPattern[]] := Module[
  {state, v, e, pos }
, v = VertexList[g] 
; pos = GraphEmbedding @ g 
; e = EdgeList[g] 

; state = GraphToEditorState[opt]

; state["vertex"] = Association @ Map[ (#id -> #) & ] @ MapThread[createVertex, {v, pos}]
; state = stateSnapInit @ state

; state["edge"]   = toStateEdges[state, g] 

; state[ "config", "vCounter"] = Length@v
; state[ "config", "eCounter"] = Length@e
; state[ "config", "DirectedEdges"] = ! EmptyGraphQ @ g && DirectedGraphQ @ g

; geAction["UpdateRange", Dynamic @ state]
; state = stateHandleNarrowRange @ state 

; state
]


(* ::Subsection::Closed:: *)
(*state helpers*)

stateSnapInit[state_Association] := Module[{config = state["config"], result = state }

, config["snap"] = config["SnapToGrid"] =!= False 

; If[
    config["snap"]
  , config["snapStep"] = stateGetAutomaticSnapStep @ state
  ; result["vertex"] = <|#, "pos" -> Round[#pos, config["snapStep"]] |>& /@ state["vertex"]   
  ] 

; <|result, "config" -> config |>     
]



stateGetAutomaticSnapStep[state_Association] := 
  Ceiling[#, .5]*10^#2 & @@ MantissaExponent[(#2 - #)/ $gridLinesCount] & @@@ state["config", "range"]


stateHandleNarrowRange[state_Association] := Module[
  {range = state["config", "range"], lengths, aspectRatio, center, side, newState}
  
, lengths = #2-#& @@@ state["config", "range"]
; aspectRatio = #2/# & @@ Sort @ Clip[lengths, {10.^(-15), Infinity}]

; If[ aspectRatio < $narrowAspectRatioLimit
  , Return[state, Module]
  ]
  
; center = Mean @ Transpose @ range
; side = Max @ lengths / 2.
; range = Transpose[{{-1,-1},{1,1}}*side+{center,center}]

; newState = state
; newState[ "config", "range"] = range
; newState[ "config", "inRangeQ" ] = RegionMember[ Rectangle @@ Transpose@ range ]
; newState
]


toStateEdges[state_Association, graph_] := Module[{ }
, edges = EdgeList @ graph
; vertexRules = (#name -> #id) & /@ Values @ state["vertex"]
; edges = Replace[edges, vertexRules , {2}]
; Association @ Map[ (#id -> #) & ] @ MapIndexed[createEdge["e"<>ToString@First@#2, #]&, edges]

]

stateVertexList[state_Association] := Values @ state[["vertex", All, "name"]]

stateEdgeList[state_Association] := state // 
 Query["edge", Values, ((#type /. # /. state[["vertex", All, "name"]])) &]

stateGraphEmbedding[state_Association] := state // Query["vertex", Values, "pos"]

stateHasSelectedVertex[state_Association] := StringQ @ state["selectedVertex"]


(* ::Subsection:: *)
(*UI Actions*)


(* ::Subsubsection::Closed:: *)
(*$geDebug*)


(* a debug feature, enabled by default till we have a first version*)

If[ 
  TrueQ @ $geDebug
  
, geAction[args___]:=(Beep[]; Print @ Framed @ InputForm @ {args})
; Module[{$inside = False}
  , geAction /: SetDelayed[geAction[args___], rhs_] /; !TrueQ[$inside] := Block[
      {$inside = True}
    , geAction[a:PatternSequence[args]]:=(
        Print[Style[StringPadRight[{a}[[1]], 16], Bold], ":", {a}[[3;;]]]
      ; rhs
      )
    ]
  ]
]


(* ::Subsubsection:: *)
(*UpdateVertexPosition*)


geAction["UpdateVertexPosition", Dynamic @ state_, vId_String, pos: {_, _}]:=(
  state["vertex", vId, "pos"] = pos
; If[
    ! state["config", "inRangeQ"] @ pos
  , geAction["UpdateRange", Dynamic @ state]  
  ]  
; geAction["UpdateEdgesShapes", Dynamic @ state]
)


(* ::Subsubsection:: *)
(*UpdateRange*)


geAction["UpdateRange", Dynamic @ state_] := Module[
  {newBounds, vs }
, vs = vertexSizeMultiplier @ state["config", "VertexSize"] 
; newBounds = CoordinateBounds[ state//Query["vertex", All, "pos"], Scaled[ 2.5 Max[vs, 0.05 ] ] ]
; state[ "config", "range"] = newBounds
; state[ "config", "inRangeQ" ] = RegionMember[ Rectangle @@ Transpose@ newBounds ]
; state[ "config", "realVertexSize" ] = Norm[Transpose @ newBounds] * vertexSizeMultiplier @ state["config", "VertexSize"]
]

vertexSizeMultiplier[vs_?NumericQ] := vs;
vertexSizeMultiplier[vs_] := vs /. {
  Tiny -> 0.005,
  Small -> 0.015,
  Medium -> 0.05,
  Large -> 0.15 
} /. Except[_?NumericQ] -> 0.05
(* ::Subsubsection::Closed:: *)
(*MouseClicked*)


geAction["MouseClicked", Dynamic @ state_ , pos_] := Module[{newV}

, If[
    Not @ CurrentValue["AltKey"]
  , geAction["Unselect", Dynamic @ state]
  ; Return[Null, Module]
  ]
  
; newV = geAction["AddVertex", Dynamic @ state, pos ]
  
]


(* ::Subsubsection::Closed:: *)
(*VertexClicked*)


geAction["VertexClicked", Dynamic @ state_, v_Association] := Catch @ With[
  {  selectedV = state["selectedVertex"], clickedV = v["id"]}

, wasAnySelected  = stateHasSelectedVertex @ state
; wasThisSelected = selectedV === clickedV

; If[ TrueQ @ CurrentValue["AltKey"] && ! wasAnySelected
  , Throw @ geAction["RemoveVertex", Dynamic @ state, v]
  ]

; If[ wasThisSelected
  , Throw @ geAction["Unselect", Dynamic @ state]
  ]

; If[  wasAnySelected
  , Throw @ geAction["CreateEdge", Dynamic @ state, selectedV, clickedV]
  ]
  
; geAction["Select", Dynamic @ state, clickedV]
 
]


(* ::Subsubsection::Closed:: *)
(*EdgeClicked*)


geAction["EdgeClicked", Dynamic @ state_, edge_Association] := If[
  TrueQ @ CurrentValue["AltKey"]
, geAction["RemoveEdge", Dynamic@state, edge]
]  


(* ::Subsubsection::Closed:: *)
(*Select*)


geAction["Select", Dynamic @ state_, vId_String] := state["selectedVertex"] = vId;


(* ::Subsubsection::Closed:: *)
(*Unselect*)


geAction["Unselect", Dynamic @ state_] := state["selectedVertex"] = Null


(* ::Subsubsection::Closed:: *)
(*AddVertex*)


geAction["AddVertex", Dynamic @ state_, pos:{_?NumericQ, _?NumericQ}] := Module[{vertex, name}

, state["config", "vCounter"]++

; name = Check[generateUniqueVertexName @ state, CreateUUID[]]
 
; vertex = createVertex[name, pos]
; id = vertex["id"]

; If[ state["config", "snap"]
  , newV["pos"] = MapThread[Round, {vertex["pos"], state["config", "snapStep"]} ]
  ]

; state["vertex", id ] = vertex

; If[
    stateHasSelectedVertex @ state
  , geAction["CreateEdge", Dynamic @ state, state["selectedVertex"], id]
  ]
  
; If[
    state["config", "CreateVertexSelects"]
  , geAction["Select", Dynamic @ state, id] 
  ]
  
; vertex
]

generateUniqueVertexName[state_Association] := Module[{names}
  (*TODO: for big graphs just generate uuid*)
, names = stateVertexList @ state
; smallestMissingInteger @ Cases[names, _Integer]
]

smallestMissingInteger[list_List] := Block[{n = 1}
, Catch[
    Scan[
      If[# == n, n++, Throw@n] &
    , Sort @ list
    ]
  ; n
  ]
]




(* ::Subsubsection::Closed:: *)
(*RemoveVertex*)


geAction["RemoveVertex", Dynamic @ state_, v_] := With[{ id = v["id"]}
, state["vertex"]  = KeyDrop[id] @ state["vertex"]
; state["edge"]    = DeleteCases[state["edge"], KeyValuePattern[{_ -> id}] ]
; If[ state["selectedVertex"] === id, geAction["Unselect", Dynamic @ state] ] 
]


(* ::Subsubsection::Closed:: *)
(*CreateEdge*)


geAction["CreateEdge", Dynamic @ state_, selectedV_String, clickedV_String] := Module[{eId, type}
, If[ 
    Not @ newEdgeAllowedQ[state, selectedV, clickedV]
  , Beep[]
  ; Message[IGGraphEditor::multiEdge]
  ; Return[$Failed, Module]
  ]  

; eId = "e"<>ToString[++state["config", "eCounter"]]
; type =   If[state["config", "DirectedEdges"], Rule, UndirectedEdge]

; state["edge", eId ] = createEdge[eId,  type[selectedV, clickedV] ]
  
; geAction["Unselect", Dynamic @ state]  
; geAction["UpdateEdgesShapes", Dynamic @ state]

]


(* ::Subsubsection::Closed:: *)
(*RemoveEdge*)


geAction["RemoveEdge", Dynamic @ state_, edge_Association]:=(
  state["edge"] = KeyDrop[edge["id"]] @ state["edge"]
; geAction["UpdateEdgesShapes", Dynamic @ state]  
)


(* ::Subsubsection::Closed:: *)
(*ToggleEdgeType*)


geAction["ToggleEdgeType", Dynamic@state_, edge_] := With[{ type := state["edge", edge["id"], "type"] }
,  type = nextEdgeType @ type
]


(nextEdgeType[#]=#2)& @@@ Partition[{UndirectedEdge["v1", "v2"], "v1"->"v2", "v2"->"v1"}, 2, 1, {1, 1}]


(* ::Subsubsection:: *)
(*UpdateEdgesShapes*)


(*(Hold|Dynamic) is there to workaround a bug with Interpretation's Initialization which inserts evaluated Dynamic's arguments*)
geAction["UpdateEdgesShapes", (Hold|Dynamic) @ state_] := Module[{primitives, graph, vertexEncoded, vertexList}
, vertexList = state // Query["vertex", Values, "id"]

; vertexEncoded = AssociationThread[
    ArrayComponents[vertexList] -> Thread[DynamicLocation[vertexList, Automatic]]
  ]
      

; graph = Graph[
    vertexList
  , state // Query["edge", Values, Tooltip[#type /. #,#id] &]
  , VertexCoordinates->(state // Query["vertex", Values, "pos"])
  ]
  
; primitives = Normal[ ToExpression@ToBoxes[graph] ]//
  Cases[#,Tooltip[prim_, eId_, ___] :> (<|"id"->eId,"primitive"->prim|>), Infinity]&;
; ( state["edge", #id, "shape"] = ToEdgeShapeFunction[#primitive, vertexEncoded] )& /@ primitives        
]



ToEdgeShapeFunction[p:{Arrowheads[0.],Arrow[b_BezierCurve, ___]}, vertexEncoded_]:={Arrowheads[0.],Arrow[b/.(vertexEncoded)]};
ToEdgeShapeFunction[Arrow[b_BezierCurve, ___], vertexEncoded_] := Arrow[b/.(vertexEncoded)];
ToEdgeShapeFunction[p_, ___]:=(Automatic);


(* ::Section::Closed:: *)
(*helpers*)


$namePatt = _ ;


createVertex[name:$namePatt, pos:{_, _}]:= <|"name"->name, "id" -> CreateUUID[],  "pos" -> pos|>


createEdge[eId_String, (e_[v1_String,  v2_String])]:=<|
    "v1"    -> v1
  , "v2"    -> v2
  , "id"    -> eId
  , "type"  -> (e /. DirectedEdge->Rule)["v1", "v2"]
  , "shape" -> Automatic
  |>


newEdgeAllowedQ::usage = "Is supposed to test whether a new edge can be created";

newEdgeAllowedQ[state_, v1_String, v2_String] := Module[{edges}
, edges = Values @ state["edge"]
; If[ 
    state["config", "DirectedEdges"]
  , Not @ MemberQ[ edges , KeyValuePattern[{"v1" -> v1, "v2" -> v2, "type" -> ("v1"->"v2")}] ]
  , Not @ MemberQ[ edges , KeyValuePattern[{ _   -> v1,  _   -> v2, "type" -> _UndirectedEdge}] ]
  ]
  
]
