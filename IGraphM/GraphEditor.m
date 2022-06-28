(* Mathematica Package *)
(* Created by Mathematica plugin for IntelliJ IDEA *)

(* :Author: Kuba Podkalicki  *)
(* :Date: 2020-09-30 *)

Package["IGraphM`"]


PackageExport["IGGraphEditor"]

(* Dev notes
 - foo[ Dynamic @ state_ ] is no better than HoldFirst etc. I just like this convention for GUI 
*)
IGGraphEditor::usage        =
    "IGGraphEditor[] typesets to an interactive graph editor.\n" <>
    "IGGraphEditor[graph] uses the given graph as the starting point.";
IGGraphEditor::multiEdge    = "Multi-edges are not supported yet. Only directed pairs are {1->2, 2->1}.";
IGGraphEditor::unknownState = "Corrupted editor state."
IGGraphEditor::oldVer       = "You need to update IGraph/M to continue working with the data stored here."


IGGraphEditor // Options = {
  "KeepVertexCoordinates" -> True
, "DefaultEdgeType"       -> "Directed" (* | or else undirected *)  
, "CreateVertexSelects"   -> True
(*, "QuantizeVertexPosition"-> False*)
, "IndexGraph"            -> False
}

SyntaxInformation[IGGraphEditor] = {"ArgumentsPattern" -> {OptionsPattern[]}};

IGGraphEditor /:
  MakeBoxes[IGGraphEditor[graph:(_|PatternSequence[]), opt:OptionsPattern[]], StandardForm] :=
  ToBoxes @ iGraphEditor[graph, opt]


iGraphEditor // Options = Options @ IGGraphEditor;


supportedGraphQ = ! MixedGraphQ[#] && SimpleGraphQ[#]&


iGraphEditor[graph:((_? supportedGraphQ)|PatternSequence[]), opt:OptionsPattern[]]:=Interpretation[
  { 
    state = GraphToEditorState[graph, opt]
  , error = False 
  , refresh
  }
, refresh[]:= Module[{temp}, Catch[
    If[
      Needs@"IGraphM`" === $Failed
    , Throw[ error = Failure["GraphEditor", <|"Message" -> "\[WarningSign] Failed to load IGraphM`, make sure it is installed." |>]]
    ]
  
  ; temp = geStateVersionCheck @ state
  ; If[ AssociationQ @ temp
    , state = temp
    , Throw[ error = temp]
    ]  
    
  ; error = False        
  ; geAction["UpdateEdgesShapes", Hold @ state]  
    
  ]]

; Deploy @  Panel[#, FrameMargins->0, BaseStyle->CacheGraphics->False]& @
  PaneSelector[
  { 
    True -> Button[Dynamic @ error,  refresh[]  ]
  , False  -> Dynamic[Refresh[iGraphEditorPanel[Dynamic@state], None]]
  }  
  , Dynamic[ MatchQ[_Failure] @ error ]
  , ImageSize->Automatic
  
  ]
  
, GraphFromEditorState @ state
  
, Initialization :> refresh[]
]


iGraphEditor[_Graph, OptionsPattern[]]:= Failure["GraphEditor", <|"Message" -> "Currently a Graph argument needs to be Simple And Not Mixed."|>]


iGraphEditor[___]:= Failure["GraphEditor", <|"Message" -> "Unknown input"|>]


iGraphEditorPanel[Dynamic@state_]:=EventHandler[
    geGraphics @ Dynamic @ state
  , "MouseClicked" :> (
      geAction["MouseClicked", Dynamic @ state, CurrentValue[{"MousePosition", "Graphics"}]]
    )
  , PassEventsDown->True
  ]


geGraphics[Dynamic @ state_ ]:= Graphics[
  DynamicNamespace @ {
  
    geHighlightsPrimitives @ Dynamic @ state

  , Gray

  , Dynamic @ geEdges @ Dynamic@state  

  , Dynamic @ Table[ geVertexShapeFunction[Dynamic@state, state["vertex", id] ], {id, Keys @ state["vertex"] }]

  }
, PlotRange -> state["config", "coordinateBounds"]
, ImagePadding -> 14
]


$stateVersion = 1; 
(*It does not change unless the structure of state association is changed. 
  That implies new geActions etc won't be compatible with old state *)


(*ONCE $stateVersion > 1 *)
(*
  geStateVersionCheck[ state : KeyValuePattern[{"version" \[Rule] oldVersion}] ] := (*oldVersion to oldVersion+1 transition rules*)
*)


geStateVersionCheck[ state : KeyValuePattern[{"version" -> $stateVersion}] ] := state


geStateVersionCheck[ state : KeyValuePattern[{"version" -> v_}] ]:= 
  Failure["GraphEditor", <|"MessageTemplate" -> IGGraphEditor::oldVer|>]


geStateVersionCheck[___]:= 
  Failure["GraphEditor", <|"MessageTemplate" -> IGGraphEditor::unknownState|>]


GraphFromEditorState[state_Association]:= GraphFromEditorState[state, Lookup[state, "version", 1]];


GraphFromEditorState[state_, 1]:=Module[{v,e,pos, graph}
, v = stateVertexList @ state

; e = stateEdgeList @ state

; pos = If[
    state["config", "keepVertexCoordinates"]  
  , stateGraphEmbedding @ state 
  , Automatic
  ]
  
; graph = Graph[ v, e, VertexCoordinates -> pos]  

; If[
    TrueQ @ state["config", "indexGraph"]
  , graph = IndexGraph @ graph
  ] 

; graph  
]


GraphToEditorState[ opt:OptionsPattern[]]:=<|
      "vertex"         -> <||>
    , "edge"           -> <||>
    , "selectedVertex" -> Null
    , "version"        -> $stateVersion
    , "config"         -> <|
        optionsToConfig[opt]
      , "vCounter"->0
      , "eCounter" ->0
      , "coordinateBounds" -> {{-1, 1}, {-1, 1}}|>
    |>;
   

GraphToEditorState[g_Graph ? supportedGraphQ, opt:OptionsPattern[]]:= Module[{state,v,e, pos, quant}
, v = VertexList[g] 
; pos = GraphEmbedding @ g 
(*; quant = OptionValue["QuantizeVertexPosition"]
; If[ NumericQ @ quant, pos = Round[pos, quant]]*)
; e = EdgeList[g] 

; state = GraphToEditorState[opt]

; state["vertex"] = Association @ Map[ #id -> # & ] @ MapThread[createVertex, {v, pos}] 
; state["edge"]   = toStateEdges[state, g] 

; state[ "config", "vCounter"] = Length@v
; state[ "config", "eCounter"] = Length@e
; state[ "config", "coordinateBounds"] = CoordinateBounds @ pos
; state[ "config", "defaultEdgeType"] = If[ UndirectedGraphQ @ g, "Undirected", "Directed"]      

; state
]

toStateEdges[state_Association, graph_]:=Module[{ }
, edges = EdgeList @ graph
; vertexRules = #name -> #id & /@ Values @ state["vertex"]
; edges = Replace[edges, vertexRules , {2}]

; Association @ Map[ #id -> # & ] @ MapIndexed[createEdge["e"<>ToString@First@#2, #]&, edges]

]

stateVertexList[state_Association]:= Values @ state[["vertex", All, "name"]]

stateEdgeList[state_Association] := state // 
 Query["edge", Values, ((#type /. # /. state[["vertex", All, "name"]])) &]

stateGraphEmbedding[state_Association]:= state // Query["vertex", Values, "pos"]

stateHasSelectedVertex[state_Association]:= StringQ @ state["selectedVertex"] 

$namePatt = _ ;

optionsToConfig // Options = Options @ IGGraphEditor;

optionsToConfig[OptionsPattern[]]:= Association[
  Decapitalize[#] -> OptionValue[#] & /@ Keys @ Options[IGGraphEditor]
]


geVertexShapeFunction[Dynamic @ state_, v_Association]:= 
DynamicModule[
  {x = v@"pos", ef = 1},
Module[
  {mouseDragged, quantization = state["config", "quantizeVertexPosition"]}

, mouseDragged = If[ NumericQ @ #
    , "MouseDragged" :> (x = Round[CurrentValue[{"MousePosition", "Graphics"}], #])
    , "MouseDragged" :> FEPrivate`Set[x , FrontEnd`CurrentValue[{"MousePosition", "Graphics"}] ]
  ] & @ quantization 

; EventHandler[
  { EdgeForm @ AbsoluteThickness @ Dynamic @ ef
  , DynamicName[ 
      Disk[Dynamic[x], Offset[8]]
    , v["id"]
    ]
  }
, {
    "MouseClicked" :> geAction["VertexClicked", Dynamic @ state, v]
  , "MouseEntered" :> FEPrivate`Set[ef, 3]
  , "MouseExited"  :> FEPrivate`Set[ef, 1]
  , mouseDragged
  , "MouseUp"      :> geAction["UpdateVertexPosition", Dynamic @ state, v["id"], x]      
  }
  , PassEventsUp -> False
  ]  
]]



geEdges[Dynamic @ state_]:=Table[ 
  geEdgeShapeFunction[ Dynamic@state,  state["edge", id]], {id, Keys @ state["edge"]}
]  


geEdgeShapeFunction[Dynamic @ state_, e_Association]:= DynamicModule[
  {ef = 3}
, EventHandler[
      { AbsoluteThickness @ Dynamic @ ef
      , edgeToPrimitive @ e        
      }
    , {
        "MouseEntered" :> FEPrivate`Set[ef, 5]
      , "MouseExited"  :> FEPrivate`Set[ef, 3]
      , "MouseClicked" :> (geAction["EdgeClicked", Dynamic @ state, e])
      }
    , PassEventsUp -> False (* edgeclicked should not be followed by outer mouseclicked*)
    ]
]


edgeToPrimitive[e_]:=If[ 
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


geHighlightsPrimitives[Dynamic @ state_]:= With[{ selV := state["selectedVertex"] }
, { Red, EdgeForm@Red,
    Dynamic[
      If[ 
        stateHasSelectedVertex @ state
      , { 
          Disk[DynamicLocation[selV], Offset[12]]
        , Dashed, Line[{
            DynamicLocation[selV]
          , FrontEnd`MousePosition["Graphics",DynamicLocation[selV]]
          }]
        }
      , {}  
      ]
    ]
}    
]

(* a debug feature, enabled by default till we have a first version*)

If[ 
  TrueQ @ $geDebug, 

Module[{$inside = False}
, geAction /: SetDelayed[geAction[args___], rhs_] /; !TrueQ[$inside] := Block[
    {$inside = True}
  , geAction[a:PatternSequence[args]]:=(
      Print[Style[StringPadRight[{a}[[1]], 15], Bold], ":", {a}[[3;;]]]
    ; rhs
    )
  ]
]
]

geAction["UpdateVertexPosition", Dynamic @ state_, vId_String, pos: {_, _}]:=(
  state["vertex", vId, "pos"] = pos
; geAction["UpdateEdgesShapes", Dynamic @ state]
)

geAction["MouseClicked", Dynamic @ state_ , pos_] := Module[{newV}

, If[
    Not @ CurrentValue["AltKey"]
  , geAction["Unselect", Dynamic @ state]
  ; Return[Null, Module]
  ]
  
; newV = geAction["AddVertex", Dynamic @ state, pos ]
  
]


geAction["VertexClicked", Dynamic @ state_, v_Association]:= With[
  {  selectedV = state["selectedVertex"], clickedV = v["id"], hasSelectedVertex = stateHasSelectedVertex @ state }
, Which[
    TrueQ @ CurrentValue["AltKey"], If[ Not @ hasSelectedVertex, geAction["RemoveVertex", Dynamic @ state, v], Beep[]]
  , selectedV === clickedV        , geAction["Unselect", Dynamic @ state]
  , ! hasSelectedVertex           , geAction["Select", Dynamic @ state, clickedV]
  , True                          , geAction["CreateEdge", Dynamic @ state, selectedV, clickedV]                                   
  ]
]


geAction["EdgeClicked", Dynamic @ state_, edge_Association]:=If[
  TrueQ @ CurrentValue["AltKey"]
, geAction["RemoveEdge", Dynamic@state, edge]
(*, geAction["ToggleEdgeType", Dynamic@state, edge]*)
]  


geAction["Select", Dynamic @ state_, vId_String]:= state["selectedVertex"] = vId;


geAction["Unselect", Dynamic @ state_]:= state["selectedVertex"] = Null


geAction["AddVertex", Dynamic @ state_, pos:{_?NumericQ, _?NumericQ}]:=Module[{vertex, name}

, state["config", "vCounter"]++

; name = Check[generateUniqueVertexName @ state, CreateUUID[]]
 
; vertex = createVertex[name, pos]
; id = vertex["id"]

; If[ NumericQ @ #
  , newV["pos"] = Round[vertex["pos"], #]
  ]& @ state["config", "quantizeVertexPosition"]

; state["vertex", id ] = vertex

; If[
    stateHasSelectedVertex @ state
  , geAction["CreateEdge", Dynamic @ state, state["selectedVertex"], id]
  ]
  
; If[
    state["config", "createVertexSelects"]
  , geAction["Select", Dynamic @ state, id] 
  ]
  
; vertex
]

generateUniqueVertexName[state_Association]:=Module[{names}
  (*TODO: for big graphs just generate uuid*)
, names = stateVertexList @ state
; smallestMissingInteger @ Cases[names, _Integer]
]

smallestMissingInteger[list_List]:=Block[{n = 1}
, Catch[
    Scan[
      If[# == n, n++, Throw@n] &
    , Sort @ list
    ]
  ; n
  ]
]




geAction["RemoveVertex", Dynamic @ state_, v_]:=With[{ id = v["id"]}
, state["vertex"]  = KeyDrop[id] @ state["vertex"]
; state["edge"]    = DeleteCases[state["edge"], KeyValuePattern[{_ -> id}] ]
; If[ state["selectedVertex"] === id, geAction["Unselect", Dynamic @ state] ] 
]


geAction["CreateEdge", Dynamic @ state_, selectedV_String, clickedV_String]:=Module[{eId, type}
, If[ 
    Not @ newEdgeAllowedQ[state, selectedV, clickedV]
  , Beep[]
  ; Message[IGGraphEditor::multiEdge]
  ; Return[$Failed, Module]
  ]  

; eId = "e"<>ToString[++state["config", "eCounter"]]
; type =   If[state["config", "defaultEdgeType"] === "Directed", Rule, UndirectedEdge]

; state["edge", eId ] = createEdge[eId,  type[selectedV, clickedV] ]
  
; geAction["Unselect", Dynamic @ state]  
; geAction["UpdateEdgesShapes", Dynamic @ state]

]


geAction["RemoveEdge", Dynamic @ state_, edge_Association]:=(
  state["edge"] = KeyDrop[edge["id"]] @ state["edge"]
; geAction["UpdateEdgesShapes", Dynamic @ state]  
)


geAction["ToggleEdgeType", Dynamic@state_, edge_]:= With[{ type := state["edge", edge["id"], "type"] }
,  type = nextEdgeType @ type
]


(nextEdgeType[#]=#2)& @@@ Partition[{UndirectedEdge["v1", "v2"], "v1"->"v2", "v2"->"v1"}, 2, 1, {1, 1}]


(*(Hold|Dynamic) is there to workaround a bug with Interpretation's Initialization which inserts evaluated Dynamic's arguments*)
geAction["UpdateEdgesShapes", (Hold|Dynamic) @ state_]:=Module[{primitives,graph,vertexEncoded,vertexList}
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
  Cases[#,Tooltip[prim_,eId_,___]:>(<|"id"->eId,"primitive"->prim|>),Infinity]&;
; ( state["edge", #id, "shape"] = ToEdgeShapeFunction[#primitive, vertexEncoded] )& /@ primitives        
]



ToEdgeShapeFunction[p:{Arrowheads[0.],Arrow[b_BezierCurve,___]}, vertexEncoded_]:={Arrowheads[0.],Arrow[b/.(vertexEncoded)]};
ToEdgeShapeFunction[Arrow[b_BezierCurve,___], vertexEncoded_]:=Arrow[b/.(vertexEncoded)];
ToEdgeShapeFunction[p_, ___]:=(Automatic);


geAction[args___]:=(Beep[]; Print @ Framed @ InputForm @ {args})


newEdgeAllowedQ::usage = "Is supposed to test whether a new edge can be created";

newEdgeAllowedQ[state_, v1_String, v2_String]:= Module[{newEdgeType,edges}
, newEdgeType = state["config", "defaultEdgeType"]
; edges = Values @ state["edge"]
; If[ 
    newEdgeType === "Directed"
  , Not @ MemberQ[ edges , KeyValuePattern[{"v1" -> v1, "v2" -> v2, "type" -> "v1"->"v2"}] ]  
  , Not @ MemberQ[ edges , KeyValuePattern[{   _ -> v1,  _   -> v2, "type" -> _UndirectedEdge}] ]  
  ]
  
]


createVertex[name:$namePatt, pos:{_,_}]:= <|"name"->name, "id" -> CreateUUID[],  "pos" -> pos|>


createEdge[eId_String, (e_[v1_String,  v2_String])]:=<|
    "v1"    -> v1
  , "v2"    -> v2
  , "id"    -> eId
  , "type"  -> (e /. DirectedEdge->Rule)["v1", "v2"]
  , "shape" -> Automatic
  |>
