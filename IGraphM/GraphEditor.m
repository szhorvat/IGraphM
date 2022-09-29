(* ::Package:: *)

(* Mathematica Package *)
(* Created by Mathematica plugin for IntelliJ IDEA *)

(* :Author: Kuba Podkalicki  *)
(* :Date: 2020-09-30 *)

Package["IGraphM`"]


PackageExport["IGGraphEditor"]


(*PDynamic = Dynamic*)
(* Dev notes
 - foo[ Dynamic @ state_ ] is no better than HoldFirst etc. I just like this convention for GUI
 - state should not be modified during continuous actions
*)

IGGraphEditor::usage        =
    "IGGraphEditor[] typesets to an interactive graph editor. Use " <>
        If[$OperatingSystem === "MacOSX", "\[CommandKey]", "\[AltKey]"] <> "-click to add/remove vertices/edges.\n" <>
    "IGGraphEditor[graph] uses the given graph as the starting point.";
IGGraphEditor::multiEdge    = "Multi-edges are not supported yet. Only directed pairs are {1->2, 2->1}.";
IGGraphEditor::unknownState = "Corrupted editor state.";
IGGraphEditor::oldVer       = "You need to update IGraph/M to continue working with the data stored here.";
IGGraphEditor::nofe         = "The graph editor requires a notebook interface.";

$narrowAspectRatioLimit = N @ GoldenRatio ^ 2; (* plot range will be adjusted if the initial calculated as is above this limit *)
$vertexEdgeThickness = 0.5;
$hoverVertexEdgeThickness = 2;
$edgeThickness = 1;
$activeEdgeThickness = 3;
$potentialEdgeStyle = Directive[Purple, Dashed, Thick];
$gridLinesCount = 30.;


(* ::Subsection:: *)
(*IGGraphEditor*)


IGGraphEditor // Options = {
  "KeepVertexCoordinates" -> True (* bool *)
, "CreateVertexSelects"   -> True (* bool *)
, "SnapToGrid"            -> False (* bool *)
, "IndexGraph"            -> False (* bool *)
, "PerformanceLimit"      -> 450   (* _Integer *)
, VertexLabels            -> None (* | "Name" *)
, VertexSize              -> Small (* Tiny | Small | Medium | Large | ratioToDiagonal_?NumericQ*)
, DirectedEdges           -> False (* bool *)
, ImageSize               -> 300
};


SyntaxInformation[IGGraphEditor] = {"ArgumentsPattern" -> {_., OptionsPattern[]}};

(* Show warning when used without notebooks. *)
IGGraphEditor[___] := Module[{}, If[Not@TrueQ[$Notebooks], Message[IGGraphEditor::nofe]]; Null /; False]

IGGraphEditor /:
  MakeBoxes[IGGraphEditor[graph:(_|PatternSequence[]), opt:OptionsPattern[]], StandardForm] :=
  ToBoxes @ iGraphEditor[graph, opt]


iGraphEditor // Options = Options @ IGGraphEditor;


supportedGraphQ = ! MixedGraphQ[#] && SimpleGraphQ[#]&;


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
  },
Interpretation[
  {
    state = GraphToEditorState[graph, opt]
  , error = False
  , refresh
  }
, refresh[] := Module[{temp}
  , Catch[
      error = False
    ; If[ Needs@"IGraphM`" === $Failed, Throw[ error = packageFailure]]

    ; temp = geStateVersionCheck @ state
    ; If[ AssociationQ @ temp, state = temp, Throw[ error = temp] ]

    ; iGraphEditorInitialization[state, error]  (*can throw*)  
    ]
  ]

; PaneSelector[
  {
    True -> Button[Dynamic @ error,  refresh[], BaseStyle -> 15]
  , False -> Panel[
      Dynamic[Refresh[iGraphEditorPanel[Dynamic@state], None]]
    , FrameMargins -> 0, BaseStyle -> CacheGraphics->False
    ]
  }
  , Dynamic[ MatchQ[_Failure] @ error ]
  , ImageSize -> Automatic

  ]

, GraphFromEditorState @ state

, Initialization :> refresh[]
, Deinitialization :> iGraphEditorDeinitialization[state]

]]


iGraphEditorInitialization // Attributes = {HoldAll}
iGraphEditorInitialization[state_, error_]:=Module[
  {perfFailure = editorFailure["Too many vertices and edges. Increase \"PerformanceLimit\" option to try anyway."]}
  
, If[
    state[ "vCounter"] + state[ "eCounter"] > state[ "PerformanceLimit"]
  , Throw[ error =  perfFailure ]
  ]

; ToTrackedAssociation @ state
  
; geAction["UpdateVertexSize", Hold @ state]
; geAction["UpdateEdgesShapes", Hold @ state]
    (*(Hold) is there to workaround a bug with Interpretation's Initialization
      which inserts evaluated Dynamic's arguments
    *)
    
; {}
]

iGraphEditorDeinitialization // Attributes = {HoldAll}
iGraphEditorDeinitialization[state_]:= StopTracking @ state


iGraphEditor[_Graph, OptionsPattern[]] := Failure["GraphEditor", <|"Message" -> "The input graph must be simple and must not be mixed."|>]


iGraphEditor[___] := Failure["GraphEditor", <|"Message" -> "Unknown input."|>]


iGraphEditorPanel[Dynamic@state_] := EventHandler[
    geGraphics @ Dynamic @ state
  , "MouseClicked" :> (
      geAction["MouseClicked", Dynamic @ state, CurrentValue[{"MousePosition", "Graphics"}]]
    )
  , PassEventsDown -> True
  ]



(* ::Subsection:: *)
(*State*)


(* ::Subsubsection:: *)
(*state version*)


$stateVersion = 2;
(*It does not change unless the structure of state association is changed.
  That implies new geActions etc won't be compatible with old state *)


(*ONCE $stateVersion > 1 *)
(*
  geStateVersionCheck[ state : KeyValuePattern[{"version" \[Rule] oldVersion}] ] := (*oldVersion to oldVersion+1 transition rules*)
*)


geStateVersionCheck[ state : KeyValuePattern[{"version" -> $stateVersion}] ] := state


geStateVersionCheck[ state : KeyValuePattern[{"version" -> 1}] ]:= <|
  state // KeyDrop["config"]
, state["config"] 
, "version" -> 2
|>


geStateVersionCheck[ state : KeyValuePattern[{"version" -> v_}] ] :=
  Failure["GraphEditor", <|"MessageTemplate" -> IGGraphEditor::oldVer|>]


geStateVersionCheck[___] :=
  Failure["GraphEditor", <|"MessageTemplate" -> IGGraphEditor::unknownState|>]


(* ::Subsubsection:: *)
(*from state*)


GraphFromEditorState[state_Association] := GraphFromEditorState[geStateVersionCheck @ state, $stateVersion];


GraphFromEditorState[state_, $stateVersion] := Module[{v, e, pos, graph}
, v = stateVertexList @ state

; e = stateEdgeList @ state

; pos = If[
    state[ "KeepVertexCoordinates"]
  , stateGraphEmbedding @ state
  , Automatic
  ]

; graph = Graph[ v, e, VertexCoordinates -> pos]

; If[
    TrueQ @ state[ "IndexGraph"]
  , graph = IndexGraph @ graph
  ]

; graph
]


(* ::Subsubsection::Closed:: *)
(*to state*)


standardizeOption[VertexLabels][val_] := Replace[val, Automatic -> "Name"]
standardizeOption[_][val_] := val


GraphToEditorState[ opt:OptionsPattern[] ]:=GraphToEditorState @ Association[ Options @ IGGraphEditor, opt ]

GraphToEditorState[ opts_Association ] := Module[{state}
, state = <|
      "vertex"         -> <||>
    , "edge"           -> <||>
    , "selectedVertex" -> Null
    , "version"        -> $stateVersion
    , optionsToConfig[opts]
    , "vCounter"->0
    , "eCounter" ->0
    , "range" -> {{-1, 1}, {-1, 1}}
    , "aspectRatio" -> 1
    , "inRangeQ" -> RegionMember[ Rectangle[{-1,-1}, {1, 1}]  ]
  |>

; state = stateSnapInit @ state
; state
]



optionsToConfig[options_Association] := KeyMap[ToString] @ options



GraphToEditorState[g_Graph ? supportedGraphQ, opt:OptionsPattern[]] := Module[
  {state, v, e, pos }
, v = VertexList[g]
; pos = GraphEmbedding @ g
; e = EdgeList[g]

; state = GraphToEditorState[<| Options @ IGGraphEditor, opt |>]

; state["vertex"] = Association @ Map[ (#id -> #) & ] @ MapThread[createVertex, {v, pos}]

; state["edge"]   = toStateEdges[state, g]

; state[ "vCounter"] = Length@v
; state[ "eCounter"] = Length@e
; state[ "DirectedEdges"] = Not @ UndirectedGraphQ @ g

; geAction["UpdateRange", Dynamic @ state]


; state
]


(* ::Subsubsection:: *)
(*state helpers*)


stateSnapInit[state_Association] := Module[{newState = state }

, newState["snap"] = newState["SnapToGrid"] =!= False

; If[
    newState["snap"]
  , newState["snapStep"] = stateGetAutomaticSnapStep @ newState
  ; newState["vertex"] = <|#, "pos" -> Round[#pos, newState["snapStep"]] |>& /@ newState["vertex"]
  ]

; newState
]


stateGetAutomaticSnapStep[state_Association] :=
  Ceiling[#, .5]*10^#2 & @@ MantissaExponent[(#2 - #)/ $gridLinesCount] & @@@ state["range"]


stateHandleNarrowRange[state_Association] := Module[
  {range = state[ "range"], lengths, aspectRatio, center, side, newState}

, lengths = #2-#& @@@ state[ "range"]
; aspectRatio = #2/# & @@ Sort @ Clip[lengths, {10.^(-15), Infinity}]

; If[ aspectRatio < $narrowAspectRatioLimit
  , Return[state, Module]
  ]

; center = Mean @ Transpose @ range
; side = Max @ lengths / 2.
; range = Transpose[{{-1,-1},{1,1}}*side+{center,center}]

; newState = state
; newState[ "range"] = range
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


stateHasCurvedEdges[state_Association]:= state[ "DirectedEdges"] 
(*TODO: obviously this needs to be more precise, this is a first approximation, at least until we support multigraphs*)


(* ::Subsection:: *)
(*graphics*)


(* ::Subsubsection::Closed:: *)
(*geGraphics*)


geGraphics[Dynamic @ state_ ] := DynamicModule[
  {range = state[ "range"], size = state["ImageSize"], graphicsSize}
, Deploy @ EventHandler[
    Pane[
      DynamicWrapper[
        Graphics[
          DynamicNamespace @ {
            geHighlightsPrimitives @ Dynamic @ state
          , Gray
          , geEdges @ Dynamic @ state
          , geVertices @ Dynamic @ state
          }
        , PlotRange -> Dynamic @ range
        , ImageSize -> Dynamic @ graphicsSize 
        ]
      , range = state[ "range"]
      ; graphicsSize = size = state["ImageSize"]
      , TrackedSymbols :> {state}
      ]
    , AppearanceElements -> {"ResizeArea"}  
    , ImageSize -> Dynamic[
        size
      , (size = {1, 1/state["aspectRatio"]} * #[[1]] )&
      ]
    ]
  , "MouseUp" :> (state["ImageSize"] = size)
  , PassEventsDown -> True
  ]

] (* PlotRange->Dynamic@state["config... updates at any unrelated event, vertex dragging included,
     this DynamicModule @ DynamicWrapper is here to address a bug. Why does it help? Because WRI.
   *)


(* ::Subsubsection::Closed:: *)
(*vertex*)


geVertices[Dynamic @ state_] := PDynamic[
dynamicLog["vertices"];
Table[
  geVertexShapeFunction[Dynamic@state, state[["vertex", pos ]]  ]
, {pos, state["vCounter"] }
]
]


geVertexShapeFunction[Dynamic @ state_, v_Association] :=
With[ {
  nef = $vertexEdgeThickness
, aef = $hoverVertexEdgeThickness
, step = state["snapStep"]
},
DynamicModule[
  {x = v@"pos", task},
Module[
  {graphics}

, graphics = { {
    EdgeForm @ AbsoluteThickness @  Dynamic[ FEPrivate`If[  FrontEnd`CurrentValue["MouseOver"], aef, nef ] ]
  , DynamicName[
      Disk[Dynamic[x], state["realVertexSize"]]
    , v["id"]
    ]
  }
  , If[ (*TODO, this could be a separate collection, like vertex/edges, so it could be toggled 
          with lower overhead *)
      state["VertexLabels"] === "Name"
    , Inset[v["name"], Offset[ {12, 12}, DynamicLocation[v["id"]]] ]
    , Nothing
    ]
  }


; EventHandler[
    graphics,
    { "MouseClicked" :> (TaskRemove @ task; geAction["VertexClicked", Dynamic @ state, v] )
    , "MouseUp"      :> (task = SessionSubmit @ ScheduledTask[ geAction["UpdateVertexPosition", Dynamic @ state, v["id"], x] , {0.1}])
    , If[ state[ "snap"]
      , "MouseDragged" :> (x = Round[CurrentValue[{"MousePosition", "Graphics"}], step])
      , "MouseDragged" :> FEPrivate`Set[x , FrontEnd`CurrentValue[{"MousePosition", "Graphics"}] ]
      ]
    },
    PassEventsUp -> False
  ]
      (*ScheduledTask stuff is here to prevent MouseUp firing if MouseClicked is going to happen*)
      (*I'd prefer clicked to be Queued but if I put it in an inner queued EventHandler then I can't block MouseUp from fireing*)
]]]



(* ::Subsubsection::Closed:: *)
(*edges*)


geEdges[Dynamic @ state_] := PDynamic[
dynamicLog["edges"];
Table[
  geEdgeShapeFunction[ Dynamic@state,  state[["edge"]][[ pos ]] ]
, {pos, state["eCounter"]}
]
]


geEdgeShapeFunction[Dynamic @ state_, e_Association] := EventHandler[
    edgeHoverWrapper @ edgeToPrimitive @ e      
  , { "MouseClicked" :> (geAction["EdgeClicked", Dynamic @ state, e]) }
  , PassEventsUp -> False (* edgeclicked should not be followed by outer mouseclicked*)
  ]


edgeHoverWrapper[primitive_]:=  {
  AbsoluteThickness @  Dynamic[ FEPrivate`If[  FrontEnd`CurrentValue["MouseOver"], $activeEdgeThickness, $edgeThickness ] ]
, primitive
}


If[ $VersionNumber > 12
, edgeHoverWrapper[primitive_Arrow]:=  Mouseover[
    {AbsoluteThickness @ $edgeThickness, primitive},
    {AbsoluteThickness @ $activeEdgeThickness, primitive}
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


edgeTypeToPrimitive["v1"->"v2"] = Arrow[{#,#2}]&
edgeTypeToPrimitive["v2"->"v1"] = Arrow[{#2,#}]&
edgeTypeToPrimitive[UndirectedEdge["v1", "v2"]] = Line[{#,#2}]&


(* ::Subsubsection::Closed:: *)
(*selected*)


geHighlightsPrimitives[Dynamic @ state_] := With[{ selV := state["selectedVertex"] }
, { $potentialEdgeStyle,
    PDynamic[
      dynamicLog["selection"];
      If[
        stateHasSelectedVertex @ state
      , {  
          EdgeForm @ DeleteCases[$potentialEdgeStyle, _Dashing]
        , EdgeForm @ AbsoluteThickness[ 3 * $hoverVertexEdgeThickness ]
        , Disk[DynamicLocation[selV], state[ "realVertexSize"]]
        , Line[{
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


(* ::Subsection:: *)
(*UI Actions*)


(* ::Subsubsection:: *)
(*$geDebug*)


(* a debug feature, enabled by default till we have a first version*)

$actionLevel = -1;

If[
  TrueQ @ $logDynamic
, dynamicLog[msg_]:=Print[Style[Row[{"Updating :", msg}],Red]]
]

If[
  TrueQ @ $geDebug

, geAction[args___] := (Beep[]; Print @ Framed @ InputForm @ {args})

; Module[{$inside = False}
  , geAction /: SetDelayed[geAction[args___], rhs_] /; !TrueQ[$inside] := Block[
      {$inside = True}
    , geAction[a:PatternSequence[args]]:=Internal`InheritedBlock[{ $actionLevel = $actionLevel + 1}
      , Module[{result, start = AbsoluteTime[]}
        , logAction[a]
        ; result = rhs
        ; Print[StringJoin@ConstantArray["  ", $actionLevel], "timing: ", AbsoluteTime[] - start, "[s]"]
        ; result  
        ]      
      ]
    ]
  ]
]

logAction[head_, state_, args___]:= With[
  { indent = StringJoin @ ConstantArray["- ", $actionLevel] } 
, Print[          
    Row[{ 
      indent
    , Style[head, Bold]
    , ":"
    , args
    }, BaseStyle->LineBreakWithin->False]
  ]
]


(* ::Subsubsection::Closed:: *)
(*UpdateVertexPosition*)


geAction["UpdateVertexPosition", Dynamic @ state_, vId_String, pos: {_, _}] := (
  state["vertex", vId, "pos"] = pos
; If[
    ! state["inRangeQ"] @ pos
  , geAction["UpdateRange", Dynamic @ state]
  ]
; geAction["UpdateEdgesShapes", Dynamic @ state]
)


(* ::Subsubsection::Closed:: *)
(*UpdateRange*)


geAction["UpdateRange", Dynamic @ state_] := Module[
  {newBounds, vs }

, embedding = state // Query["vertex", All, "pos"]
; If[ Length[embedding ] < 1, Return[False, Module]]

; vs = vertexSizeMultiplier @ state[ "VertexSize"]
; newBounds = handleDegeneratedRange @ CoordinateBounds[ embedding, Scaled[ 2.5 Max[vs, 0.05 ] ] ]


; state[ "range"] = newBounds
; state = stateHandleNarrowRange @ state

; state[  "aspectRatio"] = #/#2& @@ (#2-#& @@@ state[ "range"])
; state[  "ImageSize" ] = {1, 1/state[  "aspectRatio"]} * If[ListQ@#, First@#,#]& @ state[ "ImageSize"]
; state[  "inRangeQ" ] = RegionMember[ Rectangle @@ Transpose@ state[  "range"] ]
; geAction["UpdateVertexSize", Dynamic @ state]
]

handleDegeneratedRange[range : {{xmin_, xmax_}, {ymin_, ymax_}}] := Module[{}
, If[ xmin != xmax && ymin != ymax
  , Return[ range, Module ]
  ]

; If[ xmin == xmax && ymin == ymax
  , Return[{{xmin-1, xmax+1}, {ymin-1, ymax+1}}, Module]
  ]

; If[ xmin != xmax
  , {{xmin, xmax}, {ymin-#, ymax+#}} & [ (xmax - xmin)/2 ]
  , {{xmin-#, xmax+#}, {ymin, ymax}} & [ (ymax - ymin)/2 ]
  ]
]


geAction["UpdateVertexSize", _ @ state_ ] := With[{ (* _ @ for interpretation bug fix *)
  boundingBox = Transpose @ state[ "range"]
, sizeMultiplier = vertexSizeMultiplier @ state[ "VertexSize"]
}
, state[  "realVertexSize" ] = Norm[ boundingBox ] * sizeMultiplier

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

, state[ "vCounter"]++

; name = Check[generateUniqueVertexName @ state, CreateUUID[]]

; vertex = createVertex[name, pos]
; id = vertex["id"]

; If[ state["snap"]
  , newV["pos"] = MapThread[Round, {vertex["pos"], state[ "snapStep"]} ]
  ]

; state["vertex", id ] = vertex

; If[
    stateHasSelectedVertex @ state
  , geAction["CreateEdge", Dynamic @ state, state["selectedVertex"], id]
  ]

; If[
    state[ "CreateVertexSelects"]
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
; state["vCounter"]--
; state["eCounter"] = Length @ state["edge"]
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

; eId = "e"<>ToString[++state[ "eCounter"]]
; type =   If[state["DirectedEdges"], Rule, UndirectedEdge]

; state["edge", eId ] = createEdge[eId,  type[selectedV, clickedV] ]

; geAction["Unselect", Dynamic @ state]
; geAction["UpdateEdgesShapes", Dynamic @ state]

]


(* ::Subsubsection::Closed:: *)
(*RemoveEdge*)


geAction["RemoveEdge", Dynamic @ state_, edge_Association] := (
  state["edge"] = KeyDrop[edge["id"]] @ state["edge"]
; geAction["UpdateEdgesShapes", Dynamic @ state]
; state["eCounter"]--  
)


(* ::Subsubsection::Closed:: *)
(*ToggleEdgeType*)


geAction["ToggleEdgeType", Dynamic@state_, edge_] := With[{ type := state["edge", edge["id"], "type"] }
,  type = nextEdgeType @ type
]


(nextEdgeType[#]=#2)& @@@ Partition[{UndirectedEdge["v1", "v2"], "v1"->"v2", "v2"->"v1"}, 2, 1, {1, 1}]


(* ::Subsubsection:: *)
(*UpdateEdgesShapes*)


geAction["UpdateEdgesShapes", _ @ state_] := Module[{primitives,  vertexEncoded, vertexList}

, If[ ! stateHasCurvedEdges @ state, Return[False, Module]]

; primitives = extractEdgePrimitives @ state

; vertexList = state // Query["vertex", Values, "id"]
; vertexEncoded = AssociationThread[
    ArrayComponents[vertexList] -> Thread[DynamicLocation[vertexList, Automatic]]
  ]

; ( state["edge", #id, "shape"] = ToEdgeShapeFunction[#primitive, vertexEncoded] )& /@ primitives
]


extractEdgePrimitives[state_] := Module[{graph, vertexList, edgeList, embedding}

, vertexList = state // Query["vertex", Values, "id"]
; edgeList = state // Query["edge", Values, Tooltip[#type /. #,#id] &]
; embedding = state // Query["vertex", Values, "pos"]

; graph = Graph[vertexList, edgeList, VertexCoordinates -> embedding]

; Cases[
    Normal @ ToExpression @ ToBoxes @ graph , #, Infinity
  ]& /@ {
    Tooltip[prim_, eId_, ___] :> <| "id" -> eId, "primitive" -> prim |>
  , TooltipBox[prim_, eId_, ___] :> <| "id" -> ToExpression @ eId, "primitive" -> (prim /. ArrowBox -> Arrow /. BezierCurveBox -> BezierCurve) |>
  } // Flatten

]


ToEdgeShapeFunction[p : {Arrowheads[0.], Arrow[b_BezierCurve, ___]}, vertexEncoded_] := {Arrowheads[0.], Arrow[b /. vertexEncoded]};
ToEdgeShapeFunction[Arrow[b_BezierCurve, ___], vertexEncoded_] := Arrow[b /. vertexEncoded];
ToEdgeShapeFunction[p_, ___] := Automatic;


(* ::Section::Closed:: *)
(*helpers*)


$namePatt = _ ;


createVertex[name:$namePatt, pos:{_, _}]:= <|"name" -> name, "id" -> CreateUUID[],  "pos" -> pos|>


createEdge[eId_String, (e_[v1_String,  v2_String])] := <|
    "v1"    -> v1
  , "v2"    -> v2
  , "id"    -> eId
  , "type"  -> (e /. DirectedEdge -> Rule)["v1", "v2"]
  , "shape" -> Automatic
  |>


newEdgeAllowedQ::usage = "Is supposed to test whether a new edge can be created";

newEdgeAllowedQ[state_, v1_String, v2_String] := Module[{edges}
, edges = Values @ state["edge"]
; If[
    state["DirectedEdges"]
  , Not @ MemberQ[ edges , KeyValuePattern[{"v1" -> v1, "v2" -> v2, "type" -> ("v1"->"v2")}] ]
  , Not @ MemberQ[ edges , KeyValuePattern[{ _   -> v1,  _   -> v2, "type" -> _UndirectedEdge}] ]
  ]

]
