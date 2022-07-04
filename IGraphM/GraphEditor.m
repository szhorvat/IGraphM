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

supportedGraphQ = MatchQ[_Graph]


(* ::Subsection:: *)
(*IGGraphEditor*)


IGGraphEditor // Options = {
  "KeepVertexCoordinates" -> True
, "CreateVertexSelects"   -> True
(*, "QuantizeVertexPosition"-> False*)
, "IndexGraph"            -> False
, "PerformanceLimit"      -> 450
, VertexLabels            -> False
, VertexSize              -> Small
, DirectedEdges           -> False
, ImageSize               -> Automatic
}

SyntaxInformation[IGGraphEditor] = {"ArgumentsPattern" -> {_., OptionsPattern[]}};

IGGraphEditor /:
  MakeBoxes[IGGraphEditor[graph:(_|PatternSequence[]), opt:OptionsPattern[]], StandardForm] :=
  ToBoxes @ iGraphEditor[graph, opt]


iGraphEditor // Options = Options @ IGGraphEditor;





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


(* ::Subsubsection::Closed:: *)
(*vertex*)


geVertices[Dynamic @ state_] := Table[
  geVertexShapeFunction[Dynamic@state, state["vertex", id] ]
, {id, Keys @ state["vertex"] }
]


geVertexShapeFunction[Dynamic @ state_, v_Association] :=
With[ {
  nef = $vertexEdgeThickness
, aef = $hoverVertexEdgeThickness  
},
DynamicModule[
  {x = v@"pos", ef = nef, id = v["id"]},
Module[
  {mouseDragged, quantization = state["config", "QuantizeVertexPosition"]}

, mouseDragged = If[ NumericQ @ #
    , "MouseDragged" :> (x = Round[CurrentValue[{"MousePosition", "Graphics"}], #])
    , "MouseDragged" :> FEPrivate`Set[x , FrontEnd`CurrentValue[{"MousePosition", "Graphics"}] ]
  ] & @ quantization 

; EventHandler[
  { {EdgeForm @ AbsoluteThickness @ Dynamic @ ef
  , DynamicName[ 
      Disk[Dynamic[x], state["config", "realVertexSize"]] 
    , id
    ]
  }
  , Black
  , If[ 
      state["config", "VertexLabels"] === "Name"
    , Inset[v["name"], Offset[ {12, 12}, DynamicLocation[id]] ]
    , Nothing
    ]
  }
, {
    "MouseClicked" :> geAction["VertexClicked", Dynamic @ state, v]
  , "MouseEntered" :> FEPrivate`Set[ef, aef]
  , "MouseExited"  :> FEPrivate`Set[ef, nef]
  , mouseDragged
  , "MouseUp"      :> geAction["UpdateVertexPosition", Dynamic @ state, id, x]      
  }
  , PassEventsUp -> False
  ]  
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
DynamicModule[
  {ef = nef}
, EventHandler[
      { AbsoluteThickness @ Dynamic @ ef
      , edgeToPrimitive @ e        
      }
    , {
        "MouseEntered" :> FEPrivate`Set[ef, aef]
      , "MouseExited"  :> FEPrivate`Set[ef, nef]
      , "MouseClicked" :> (geAction["EdgeClicked", Dynamic @ state, e])
      }
    , PassEventsUp -> False (* edgeclicked should not be followed by outer mouseclicked*)
    ]
]]


edgeToPrimitive[e_] := Module[{shapeFunction}

, shapeFunction = If[ e["shape"] =!= Automatic
  , Return[ e["shape"], Module ]
  ]

; edgeTypeToPrimitive[ e["edge"] ][
    DynamicLocation[e[["edge", 1]], Automatic],
    DynamicLocation[e[["edge", 2]], Automatic]
  ]
]  


edgeTypeToPrimitive[(Rule|DirectedEdge)[_,_]]=Arrow[{#,#2}]&
edgeTypeToPrimitive[_[_, _]]=Line[{#,#2}]&


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


(* ::Subsection::Closed:: *)
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
  {state, v, e, pos(*, quant*)}
, v = VertexList[g] 
; pos = GraphEmbedding @ g 
(*; quant = OptionValue["QuantizeVertexPosition"]
; If[ NumericQ @ quant, pos = Round[pos, quant]]*)
; e = EdgeList[g] 

; state = GraphToEditorState[opt]

; state["vertex"] = Association @ Map[ (#id -> #) & ] @ MapThread[createVertex, {v, pos}]
; state["edge"]   = toStateEdges[state, g] 

; state[ "config", "vCounter"] = Length@v
; state[ "config", "eCounter"] = Length@e
; state[ "config", "DirectedEdges"] = ! EmptyGraphQ @ g && DirectedGraphQ @ g
; state
]


(* ::Subsection::Closed:: *)
(*state helpers*)


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
 Query["edge", Values, (#edge /. state[["vertex", All, "name"]]) &]

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
; geAction["UpdateEdgesShapes", Dynamic @ state]
)


(* ::Subsubsection:: *)
(*UpdateRange*)


geAction["UpdateRange", Dynamic @ state_, pts_List]:= Module[{newBounds, vs }

, If[
    NumericQ @ state[ "config", "realVertexSize" ] &&
      And @@ (state["config", "inRangeQ"] @ pts)
  , Return[False, Module]
  ]  
  
; vs = vertexSizeMultiplier @ state["config", "VertexSize"] 
; newBounds = CoordinateBounds[ pts, Scaled[ 2.5 Max[vs, 0.05 ] ] ]
; state[ "config", "range"] = newBounds
; state[ "config", "inRangeQ" ] = RegionMember[ Rectangle @@ Transpose@ newBounds ]
; state[ "config", "realVertexSize" ] = Norm[Transpose @ newBounds] * vertexSizeMultiplier @ state["config", "VertexSize"]
; state = stateHandleNarrowRange @ state 
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


geAction["VertexClicked", Dynamic @ state_, v_Association] := With[
  {  selectedV = state["selectedVertex"], clickedV = v["id"], hasSelectedVertex = stateHasSelectedVertex @ state }
, Which[
    TrueQ @ CurrentValue["AltKey"], If[ Not @ hasSelectedVertex, geAction["RemoveVertex", Dynamic @ state, v], Beep[]]
  , selectedV === clickedV        , geAction["Unselect", Dynamic @ state]
  , ! hasSelectedVertex           , geAction["Select", Dynamic @ state, clickedV]
  , True                          , geAction["CreateEdge", Dynamic @ state, selectedV, clickedV]                                   
  ]
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

; If[ NumericQ @ #
  , newV["pos"] = Round[vertex["pos"], #]
  ]& @ state["config", "QuantizeVertexPosition"]

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

, eId = "e"<>ToString[++state["config", "eCounter"]]
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




(* ::Subsubsection:: *)
(*UpdateEdgesShapes*)


(*(Hold|Dynamic) is there to workaround a bug with Interpretation's Initialization which inserts evaluated Dynamic's arguments*)
geAction["UpdateEdgesShapes", (Hold|Dynamic) @ state_] := Module[
  {
    primitives, graph, vertexEncoded, vertexList, coordinatesRules, graphGraphics, coordinates
  }
  (*TODO, skip if no 'curved' edges are needed*)
  
, vertexList = state // Query["vertex", Values, "id"]      

; graph = Graph[
    vertexList
  , state // Query["edge", Values, Tooltip[#edge, #id] &]
  , VertexCoordinates->(state // Query["vertex", Values, "pos"])
  ]

; graphGraphics = ToExpression @ ToBoxes[graph]

; coordinates = graphGraphics // 
    First @ Cases[#, Verbatim[GraphicsComplex][pts_, ___] :> pts, Infinity] & 


; vertexEncoded = AssociationThread[
    ArrayComponents[vertexList] -> Thread[DynamicLocation[vertexList, Automatic]]
  ]

; coordinatesRules = coordinates // AssociationThread[ArrayComponents[#, 1] -> #] &  
; coordinatesRules = <|coordinatesRules, vertexEncoded|>    

  (* There is a bug that makes both edges in Graph[{Tooltip[1 -> 2, "A"], Tooltip[1 -> 2, "B"]} ]
      labeled with a tooltip A, we need to handle this. *)
; edgesIdsCollections = Values @ state["edge"] // 
    GroupBy[Key["edge"] -> Key["id"]] // Values //
    Map[#[[1]] -> # &] // 
    Association (* <|e1 -> {e1, e2, e3}, ... |>*)

; patchEdgeId[id_]:= ((edgesIdsCollections[id] = Rest @ # ; First @ #)&@ edgesIdsCollections[id])      
; primitives = Normal[ graphGraphics ] //
  Cases[#,Tooltip[prim_, eId_, ___] :> (<|"id" -> patchEdgeId @ eId,"primitive"->prim|>), Infinity]&
  
; primitives = primitives /. Arrow[BezierCurve[a : {___}, opt___], r___] :> RuleCondition[
    Arrow[BezierCurve[a /. coordinatesRules, opt], r]
  ]
  
; geAction["UpdateRange", Dynamic @ state, coordinates]
; ( state["edge", #id, "shape"] = ToEdgeShapeFunction[#primitive, vertexEncoded] )& /@ primitives        

; {primitives, coordinates}
]



ToEdgeShapeFunction[p:{Arrowheads[0.],Arrow[b_BezierCurve, ___]}, vertexEncoded_]:={Arrowheads[0.],Arrow[b/.(vertexEncoded)]};
ToEdgeShapeFunction[Arrow[b_BezierCurve, ___], vertexEncoded_] := Arrow[b/.(vertexEncoded)];
ToEdgeShapeFunction[p_, ___]:= Automatic


(* ::Section::Closed:: *)
(*helpers*)


$namePatt = _ ;


createVertex[name:$namePatt, pos:{_, _}]:= <|"name"->name, "id" -> CreateUUID[],  "pos" -> pos|>


createEdge[eId_String, edge:(e_[v1_String,  v2_String])]:=<|
    "edge"  -> edge
  , "id"    -> eId  
  , "shape" -> Automatic
  |>


