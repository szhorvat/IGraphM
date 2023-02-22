(* ::Package:: *)

(* ::Subsection::Closed:: *)
(*Package Export*)


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


(* ::Subsection::Closed:: *)
(*Helpers*)


es6Decorate // Attributes={HoldAll};

es6Decorate[f_Symbol]:=(

  f /: SetDelayed[
    f[ Verbatim[Association][spec__] ]
    , rhs_
  ]:=With[
    { pattern    = es6ExtractPattern[spec]
      , heldSymbols = es6ExtractSymbols[spec]
      , heldValues  = es6ExtractValues[spec]
    }
    , es6Decorate[f, rhs, pattern, heldSymbols, heldValues]
  ]
);

es6Decorate[foo_Symbol,rhs_,pattern_, Hold[symbols___], Hold[values___]]:=(
  SetDelayed @@ Hold[
    foo[$ES6asso:pattern],
    Block[{symbols}, Unevaluated[rhs] /. Thread[{symbols}->{values}] ]
  ]
)

es6ExtractPattern[spec__]:= KeyValuePattern[
  es6PatternToPattern /@ {spec}
]

es6ExtractSymbols[spec__]:= Apply[Join][
  es6PatternToHeldSymbol /@ {spec}
]

es6ExtractValues[spec__]:= Apply[Join][
  es6PatternToHeldValue /@ {spec}
]

With[
  { vPatt = Verbatim[Pattern]
    , vOpt  = Verbatim[Optional]
    , hPatt = HoldPattern
  }
  ,

  (* a_ a_H *)
  es6PatternToPattern[    hPatt @ vPatt[s_Symbol, b_Blank] ] := symbolName[s] -> b;
  (* "a" \[Rule] b_ *)
  es6PatternToPattern[    sym_String -> Except[_Optional] ]:= sym -> Blank[];
  es6PatternToPattern[ ___ ] = Sequence[];



  es6PatternToHeldSymbol[ hPatt @ vPatt[s_Symbol, _]      ] := Hold[s];
  es6PatternToHeldSymbol[ hPatt @ vOpt[p_, _]      ] := es6PatternToHeldSymbol @ p;
  es6PatternToHeldSymbol[ _String -> p_ ]:= es6PatternToHeldSymbol @ p;
  es6PatternToHeldSymbol[ ___ ]:=Hold[];



  (* s_ ==> ass["s"] *)
  es6PatternToHeldValue[  hPatt @ vPatt[s_Symbol, _]      ] := With[{sym = symbolName@s}
    , Hold @ $ES6asso @ sym
  ];
  (* s_:default ==> Lookup[asso, "s", default] *)
  es6PatternToHeldValue[  hPatt @ vOpt[vPatt[s_Symbol, _Blank], default_] ] := With[{sym = symbolName@s}
    , Hold @ Lookup[$ES6asso, sym, default]
  ];
  (* S \[Rule] s_ *)
  es6PatternToHeldValue[ key_ -> hPatt @ vPatt[_, _Blank] ]:= Hold @ $ES6asso @ key;

  (* S \[Rule] s_:default *)
  es6PatternToHeldValue[ key_ -> hPatt @ vOpt[vPatt[_, _Blank], default_] ]:=
    Hold @ Lookup[$ES6asso, key, default];
  es6PatternToHeldValue[  ___     ] := Hold[];
];

symbolName = Function[s, SymbolName @ Unevaluated[s], HoldFirst];




(* ::Subsection::Closed:: *)
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


supportedGraphQ = MatchQ[_Graph];


(* ::Subsection:: *)
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
  , False -> Dynamic[Refresh[iGraphEditorPanel[Dynamic@state], None]]
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


iGraphEditorPanel[Dynamic@state_] := Grid[{
  {
    iGraphGraphicsPanel @ Dynamic @ state
  , iGraphMenu @ Dynamic @ state  
  }
, { 
    iGraphModeSetter @ Dynamic @ state
  , ""   
  }
}, Alignment->{Left,Top}, Spacings->{0,0}]



iGraphModeSetter[Dynamic@state_]:=SetterBar[Dynamic@state["editorMode"], {"draw" -> "Draw", "edit" -> "Edit", "config" -> "Config"}]


iGraphMenu[Dynamic[state_]]:=PaneSelector[
{ "draw" -> Column[
    {
      menuButton["Adjust range", geAction["UpdateRange", Dynamic @ state, True], Appearance->"FramedPalette"],
      menuButton["Hide controls", Appearance->"FramedPalette", Enabled->False]
    }
  , Left
  , Spacings->0
  ]
, "edit" -> Panel[
    Pane[
      PDynamic @ If[
        AssociationQ @ state["selectedObject"]
      , PDynamic[dynamicLog["selectedObject"]; state["selectedObject"]  ]
      , "Click on object"
      ]
    , ImageSize-> ({Automatic, PDynamic[state["ImageSize"][[2]]] })
    ]
  , ImageMargins -> {0,0}
  , FrameMargins->{0,0}
  ]
, "config" -> Panel[
    Pane[
      Column[{
        "Editor config"
      , Dataset @ state[[{"version","vCounter","eCounter","range","aspectRatio", "editorMode"}]]      
      , Grid[{{"VertexLabels: ", PopupMenu[Dynamic@state["VertexLabels"], {None, "Name"}]}}]
      }, Left], ImageSize-> ({Automatic, PDynamic[state["ImageSize"][[2]]] })],ImageMargins -> {0,0}, FrameMargins->{0,0}]
}
, PDynamic @ state["editorMode"]
, ImageSize->Automatic
]

menuButton // Attributes = {HoldRest}
menuButton[args___]:=Button[args, Appearance->"FramedPalette"]


iGraphGraphicsPanel[Dynamic[state_]]:=Panel[
  geGraphics @ Dynamic @ state
, FrameMargins -> 0, BaseStyle -> CacheGraphics->False
]


(* ::Subsection::Closed:: *)
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


(* ::Subsubsection::Closed:: *)
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


(* ::Subsubsection:: *)
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
    , "editorMode" -> "draw"
    , "selectedObject" -> Null
  |>

; state["ImageSize"] = state["ImageSize"] // Replace[ { Automatic -> 300, n_?NumericQ :> {n,n} } ] 

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


(* ::Subsubsection::Closed:: *)
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
 Query["edge", Values, (#edge /. state[["vertex", All, "name"]]) &]


stateGraphEmbedding[state_Association] := state // Query["vertex", Values, "pos"]


stateHasSelectedVertex[state_Association] := StringQ @ state["selectedVertex"]


stateHasCurvedEdges[state_Association]:= True
(*TODO: obviously this needs to be more precise, this is a first approximation, at least until we support multigraphs
  
*)


(* ::Subsection:: *)
(*graphics*)


(* ::Subsubsection:: *)
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
    , ImageSize -> Dynamic[size]
    ]
  , {
      "MouseUp" :> If[ size != state["ImageSize"],  geAction["PaneResized", Dynamic @ state, size] ]      
    , "MouseClicked" :> ( 
        geAction["MouseClicked", Dynamic @ state, CurrentValue[{"MousePosition", "Graphics"}]] 
      )    
  }
  , PassEventsDown -> True
  ]

] (* PlotRange->Dynamic@state["range"] updates at any unrelated event, vertex dragging included,
     this DynamicModule @ DynamicWrapper is here to address a bug. Why does it help? Because WRI.
   *)


(* ::Subsubsection:: *)
(*vertex*)


geVertices[Dynamic @ state_] := DynamicModule[{vertexMoved}
, DynamicWrapper[
    PDynamic[
      dynamicLog["vertices"]
    ; Table[
        geVertexShapeFunction[Dynamic@state, state[["vertex", pos ]] , Dynamic @ vertexMoved ]
      , {pos, state["vCounter"] }
      ]
    ]
  , If[ ListQ @ vertexMoved
    , geAction["UpdateVertexPosition", Dynamic @ state, ##& @@ vertexMoved] 
    ; vertexMoved=Null
    ]
  , TrackedSymbols :> {vertexMoved  }
  ]
]

(* This should've been Dynamic @ Table only but I needed to add this 'listener' because 
   geAction comming from the ScheduledTask itself breaks DynamicModule variable system and
   DynamicModuleBox values wouldn't be updated resulting in the vertex possition being out of sync
   after move up... Amd we need scheduled tasks becasue of other bug :heart:
*)


geVertexShapeFunction[Dynamic @ state_, v_Association, Dynamic @ vertexMoved_] :=
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
      Disk[Dynamic[x], PDynamic@state["realVertexSize"]]
    , v["id"]
    ]
  }
  , PDynamic@If[ (*TODO, this could be a separate collection, like vertex/edges, so it could be toggled 
          with lower overhead *)
      state["VertexLabels"] === "Name"
    , Inset[v["name"], Offset[ {12, 12}, DynamicLocation[v["id"]]] ]
    , {}
    ]
  }


; EventHandler[
    graphics,
    { "MouseClicked" :> (RemoveScheduledTask @ task; geAction["VertexClicked", Dynamic @ state, v] )
    , "MouseUp"      :> (task = RunScheduledTask[ vertexMoved = {v["id"], x} , {0.1}])
    , If[ state[ "snap"]
      , "MouseDragged" :> (x = Round[CurrentValue[{"MousePosition", "Graphics"}], step])
      , "MouseDragged" :> Set[x , CurrentValue[{"MousePosition", "Graphics"}] ]
      ]
    },
    PassEventsUp -> False
  ]
      (*ScheduledTask stuff is here to prevent MouseUp firing if MouseClicked is going to happen*)
      (*I'd prefer clicked to be Queued but if I put it in an inner queued EventHandler 
        then I can't block MouseUp from firing*)
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
    edgeHoverWrapper @ edgeToPrimitive @ e /. Arrowheads[0.] -> Arrowheads[0.00001] (*patch #2748 :heart: *)    
  , { "MouseClicked" :> (geAction["EdgeClicked", Dynamic @ state, e]) }
  , PassEventsUp -> False (* edgeclicked should not be followed by outer mouseclicked*)
  ]


edgeHoverWrapper[primitive_]:=  With[
  {
    nef = $edgeThickness,
    aef = $activeEdgeThickness
  },{
  AbsoluteThickness @  Dynamic[ FEPrivate`If[  FrontEnd`CurrentValue["MouseOver"], aef, nef ] ]
, primitive
}
]


If[ 12 < $VersionNumber < 13.2
, edgeHoverWrapper[primitive_Arrow]:=  Mouseover[
    {AbsoluteThickness @ $edgeThickness, primitive},
    {AbsoluteThickness @ $activeEdgeThickness, primitive}
  ]
]


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
(*Logging*)


(* a debug feature, enabled by default till we have a first version*)

$actionLevel = -1;

If[
  TrueQ @ $logDynamic
, dynamicLog[msg_]:=Print[Style[Row[{"Updating :", msg}],Red]]
]

If[
  TrueQ @ $geDebug

, Module[{$inside = False}
  , geAction /: SetDelayed[geAction[args___], rhs_] /; !TrueQ[$inside] := Block[
      {$inside = True}
    , geAction[a:PatternSequence[args]]:=Internal`InheritedBlock[{ $actionLevel = $actionLevel + 1}
      , Module[{result, start = AbsoluteTime[]}
        , logAction[a]
        ; result = rhs
        ; If[$logTimings, Print[StringJoin@ConstantArray["  ", $actionLevel], "timing: ", AbsoluteTime[] - start, "[s]"]]
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

(* ::Subsection::Closed:: *)

(*Events*)

(* ::Subsubsection::Closed:: *)

(*PaneResized*)

geAction["PaneResized", Dynamic @ state_, size_] := (
  
    state["ImageSize"] = size
  ; geAction["UpdateRange", Dynamic @ state, True]  
)


(* ::Subsubsection::Closed:: *)
(*MouseClicked*)


geAction["MouseClicked", Dynamic @ state_ , pos_] := Module[{newV}

, If[ 
    state["editorMode"] === "edit"
  , Return @ geAction["UnselectObject", Dynamic @ state]
  ]

; If[
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
, Module[
  {wasAnySelected, wasThisSelected}  

, If[ state["editorMode"] === "edit"
  , Return @ geAction["SelectObject", Dynamic @ state, v]
  ]

; wasAnySelected  = stateHasSelectedVertex @ state
; wasThisSelected = selectedV === clickedV

; If[ TrueQ @ CurrentValue["AltKey"] && ! wasAnySelected
  , Throw @ geAction["RemoveVertex", Dynamic @ state, v]
  ]

(*; If[ wasThisSelected
  , Throw @ geAction["Unselect", Dynamic @ state]
  ]*)

; If[  wasAnySelected
  , Throw @ geAction["CreateEdge", Dynamic @ state, selectedV, clickedV]
  ]

; geAction["Select", Dynamic @ state, clickedV]

]]

geAction["SelectObject", Dynamic @ state_, v_Association]:= state["selectedObject"] = v
geAction["UnselectObject", Dynamic @ state_]:= state["selectedObject"] = Null


(* ::Subsubsection::Closed:: *)
(*EdgeClicked*)


geAction["EdgeClicked", Dynamic @ state_, edge_Association] := Module[{}
, If[ state["editorMode"] === "edit"
  , Return @ geAction["SelectObject", Dynamic @ state, edge]
  ]

; If[
    TrueQ @ CurrentValue["AltKey"]
  , geAction["RemoveEdge", Dynamic@state, edge]
  ]
]


(* ::Subsection::Closed:: *)

(*Actions*)

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


(* bounds      = {{x1,x2}, {y1, y2}} = (plot)range*)
(* boundig box = {{x1,y1}, {x2, y2}}*)
geAction["UpdateRange", Dynamic @ state_, force_:False] := Module[
  {newBounds, embedding, bounds }


, (*embedding = state // Query["vertex", All, "pos"]
; If[ Length[embedding ] < 1, Return[False, Module]]*)

bounds = Lookup[  state, "coordinateBounds", 
  state["coordinateBounds"] = CoordinateBounds @ state[["vertex", All, "pos"]] 
]

; newBounds = handleDegeneratedRange @ respectVertexSize[ bounds , state["VertexSize"] ]

; If[ ! force && AllTrue[ Transpose @ newBounds, state["inRangeQ"]],  Return[False, Module]]


; state["range"] = adjustRangeToImageSize[ newBounds, state["ImageSize"] ]

; state[  "inRangeQ" ] = RegionMember[ Rectangle @@ Transpose@ state[  "range"] ]
(* ; geAction["UpdateVertexSize", Dynamic @ state] *)
]


(* I don't like when UI's aspect ration changes unless it was done explicitly by draggin Pane's resize control*)
(* So we need to pad vertically or horizontally the minimal range in order to match ImageSize aspect ratio.*)

adjustRangeToImageSize[ bounds: {{x1_, x2_}, {y1_, y2_}}, size: {w_, h_}]:= Module[
  {rangeRatio, sizeRatio, newWidth, newHeight, centerWidth, centerHeigth, x1p, x2p, y1p, y2p}

, rangeRatio = (y2-y1)/(x2-x1) (*degenerate case should be handled before*)
; sizeRatio = h/w
; If[ 
    sizeRatio < rangeRatio (*pad horizontally*)

  , newWidth = w/h (y2-y1)
  ; centerWidth = .5 (x1+x2)
  ; x1p = centerWidth - .5 newWidth
  ; x2p = centerWidth + .5 newWidth
  ; {y1p, y2p} = {y1, y2}

  , newHeight = h/w (x2-x1)
  ; centerHeigth = .5 (y2+y1)
  ; y1p = centerHeigth - .5 newHeight
  ; y2p = centerHeigth + .5 newHeight
  ; {x1p, x2p} = {x1, x2}
  ]
  (* I guess 5y ago I could write a one line vector calculation but now...*)
; {rangeRatio, sizeRatio, newWidth, newHeight, centerWidth, centerHeigth, x1p, x2p, y1p, y2p}
  
; {{x1p, x2p}, {y1p, y2p}}    

]

respectVertexSize[bounds_, vertexSize_]:=With[
  { vs =  vertexSizeMultiplier @ vertexSize }
, CoordinateBounds[ Transpose @ bounds, Scaled[ 2.5 Max[vs, 0.05 ] ] ]
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


(* ::Subsubsection::Closed:: *)
(*UpdateVertexSize*)


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
; state["edge"]    = Select[state["edge"], FreeQ[#edge, id] & ]
; state["vCounter"]--
; state["eCounter"] = Length @ state["edge"]
; If[ state["selectedVertex"] === id, geAction["Unselect", Dynamic @ state] ]
]


(* ::Subsubsection::Closed:: *)
(*CreateEdge*)


geAction["CreateEdge", Dynamic @ state_, selectedV_String, clickedV_String] := Module[{eId, type}
, eId = CreateUUID["e-"]

; type = If[state["DirectedEdges"], Rule, UndirectedEdge]

; state["edge", eId ] = createEdge[eId,  type[selectedV, clickedV] ]
; state["eCounter"]++

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


(* ::Subsubsection::Closed:: *)
(*UpdateEdgesShapes*)


geAction["UpdateEdgesShapes", _ @ state_] := Module[{coordinates,edges,  vertexEncoded, vertexList}

, If[ ! stateHasCurvedEdges @ state, Return[False, Module]]

; {edges, coordinates} = Lookup[extractEdgePrimitives @ state, {"rawEdges", "coordinates"}]

; ( state["edge", #id, "shape" ] = #shape;  ) & /@ edges 

(* nested tracking is not supported yet so we need to manually trigger all edges update 
   we keep using eCounter for triggers because general mutations of .edges should not trigger updates *)
; IGraphM`PreciseTracking`PackagePrivate`UpdateTarget[state["eCounter"]]

; state["coordinateBounds"] =  CoordinateBounds @ coordinates
; geAction["UpdateRange", Dynamic @ state]
]


extractEdgePrimitives[state_] := Module[{ data = state }

, data = processGraphBoxes @ data

; data = patchEdgeIds @ data
; data = patchCurveNormal @ data

; data = addEdgesShapes @ data

; data
]


(* ::Subsubsubsection::Closed:: *)
(*processGraphBoxes*)


processGraphBoxes[state_]:=Module[{vertexList, edgeList, embedding, graph, graphGraphics, rawEdges}

, vertexList = state // Query["vertex", Values, "id"]
; edgeList   = state // Query["edge", Values, Tooltip[#edge, #id] &]
; embedding  = state // Query["vertex", Values, "pos"]

; graph = Graph[vertexList, edgeList, VertexCoordinates -> embedding]

; graphGraphics = ToExpression @ ToBoxes @ graph

; rawEdges = Cases[
    Normal @ graphGraphics , #, Infinity
  ]& /@ {
    Tooltip[prim_, eId_, ___] :> <| "id" -> eId, "primitive" -> prim |>
  , TooltipBox[prim_, eId_, ___] :> <| "id" -> ToExpression @ eId, "primitive" -> (prim /. ArrowBox -> Arrow /. BezierCurveBox -> BezierCurve) |>
  } // Flatten
  
; <|  state
  ,  <|vertexList, edgeList, embedding, graph, graphGraphics, rawEdges |> // ToKeyValue 
  |>

]


(* ::Subsubsubsection::Closed:: *)
(*patchCurveNormal*)


(* similar bug to this https://mathematica.stackexchange.com/q/105184/5478 *)

patchCurveNormal // es6Decorate

patchCurveNormal[ <|rawEdges_, vertexList_, graphGraphics_|> ]:= Module[

  {patchedEdges,vertexEncoded, coordinates, coordinatesRules }
  
, vertexEncoded = AssociationThread[
    ArrayComponents[vertexList] -> Thread[DynamicLocation[vertexList, Automatic]]
  ]

; coordinates = graphGraphics // 
    First @ Cases[#, Verbatim[GraphicsComplex][pts_, ___] :> pts, Infinity] & 
    
; coordinatesRules = coordinates // AssociationThread[ArrayComponents[#, 1] -> #] &  
; coordinatesRules = <|coordinatesRules, vertexEncoded|>    

; patchedEdges = rawEdges /. Arrow[BezierCurve[a : {___}, opt___], r___] :> RuleCondition[
    Arrow[BezierCurve[a /. coordinatesRules, opt], r]
  ] 
  
; <| $ES6asso
   , "rawEdges" -> patchedEdges
   , <| vertexEncoded, coordinates, coordinatesRules |> // ToKeyValue
  |>  
]


(* ::Subsubsubsection::Closed:: *)
(*patchEdgeIds*)


(* There is a bug that makes both edges in Graph[{Tooltip[1 -> 2, "A"], Tooltip[1 -> 2, "B"]} ]
      labeled with a tooltip A, we need to handle this. *)
      
patchEdgeIds // es6Decorate
patchEdgeIds[ <|edge_, rawEdges_|> ]:=Module[{edgesIdsCollections, patchEdgeId, patchedEdges}

, edgesIdsCollections = Values @ edge // 
    GroupBy[standardizeEdge@*Key["edge"] -> Key["id"]] // Values //
    Map[#[[1]] -> # &] // 
    Association (* <|e1 -> {e1, e2, e3}, ... |>*)


; patchEdgeId[id_]:= With[{idCollection = edgesIdsCollections[id]}
  , edgesIdsCollections[id] = Rest @ idCollection   
  ; First @ idCollection
  ]

; patchedEdges = MapAt[ patchEdgeId, rawEdges, {All, "id" }]

; <| $ES6asso, "rawEdges" -> patchedEdges |>  
]


standardizeEdge[e_UndirectedEdge]:=Sort @ e
standardizeEdge[directed_[v1_, v2_]]:= {v1, v2}


(* ::Subsubsubsection::Closed:: *)
(*addEdgesShapes*)


addEdgesShapes // es6Decorate
addEdgesShapes[ <|vertexEncoded_, rawEdges_|> ]:= <|
  $ES6asso
, "rawEdges" -> ( <|#, "shape" -> ToEdgeShapeFunction[#primitive, vertexEncoded] |>& /@ rawEdges )
|>


ToEdgeShapeFunction[p : {Arrowheads[0.], Arrow[b_BezierCurve, ___]}, vertexEncoded_] := {Arrowheads[0.], Arrow[b /. vertexEncoded]};
ToEdgeShapeFunction[Arrow[b_BezierCurve, ___], vertexEncoded_] := Arrow[b /. vertexEncoded];
ToEdgeShapeFunction[p_, ___] := Automatic;


(* ::Subsubsection::Closed:: *)
(*actions fallthrough*)

(* This has to be here becasue of a weird specificity of args___ vs PatternSequence *)
(* SetDelayed can't be used because it is decorated in debug mode *)
If[
  TrueQ @ $geDebug

,  DownValues[geAction] = Append[ 
     DownValues[geAction]
   , HoldPattern[geAction[args___] ] :> (Beep[]; Print @ Framed @ InputForm @ {args})
   ]

]

(* ::Section::Closed:: *)
(*helpers*)


$namePatt = _ ;


createVertex[name:$namePatt, pos:{_, _}]:= <|"name" -> name, "id" -> CreateUUID[],  "pos" -> pos|>


createEdge[eId_String, edge:(e_[v1_String,  v2_String])] := <|
    "id"    -> eId
  , "edge" -> edge    
  , "shape" -> Automatic
  |>


ToKeyValue::usage = "ToKeyValue[symbol] is a small utility that generates \"symbol\" -> symbol which shortens association assembling.";

ToKeyValue // Attributes = {HoldAll, Listable};

ToKeyValue[sym_Symbol]:= SymbolToKeyName[sym] -> sym;

ToKeyValue[Association[spec__Symbol] ]:= Association @ ToKeyValue @ {spec}


SymbolToKeyName::usage = "SymbolToKeyName[symbol_] generates symbol's symbol name and trims '$..nnn' if present";

SymbolToKeyName // Attributes = {HoldAll};

SymbolToKeyName[sym_Symbol]:= StringTrim[SymbolName @ Unevaluated @ sym, "$".. ~~ DigitCharacter..];
