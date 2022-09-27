(* ::Package:: *)

(* ::Text:: *)
(*Not everything can be covered by automatic tests so we should have this notebook to hands-on test problematic scenarios.*)


$Path=DeleteDuplicates@Flatten@{
$Path,
"E:\\Idea Projects\\LTemplate",
"E:\\Idea Projects\\IGraphM"
};


ClearAll["IGraphM`GraphEditor`*`*"]
NotebookDelete@Cells@MessagesNotebook[]
SetSelectedNotebook@MessagesNotebook[];
SetSelectedNotebook@EvaluationNotebook[];
FrontEndExecute[FrontEndToken["DeleteGeneratedCells"]]

IGraphM`GraphEditor`PackagePrivate`$geDebug=False;
IGraphM`GraphEditor`PackagePrivate`$logDynamic=True;
Get["IGraphM`"]
IGraphM`GraphEditor`PackagePrivate`$gridLinesCount=25.;


TestReport@FileNameJoin[{NotebookDirectory[], "Editor.wl"}]


SetOptions[IGGraphEditor, ImageSize->200];


(* ::Subsection:: *)
(*hands on tests*)


(* ::Subsubsection:: *)
(*Does it work at all?*)


(* ::Text:: *)
(*- Alt+Click build a graph*)
(*- Evaluate output to create a graph*)
(*- Drag outside of range should extend the range @mouseUp*)
(*- Alt+click should not trigger orange resize frame, nor should discarding a potential edge.*)


IGGraphEditor[VertexSize->Small]


(* ::Input:: *)
(*\!\( *)
(*TagBox[*)
(*DynamicModuleBox[{IGraphM`GraphEditor`PackagePrivate`state$$ = <|"vertex" -> <|"4f62303b-bf36-48ec-9ed3-79d7d141e6c5" -> <|"name" -> 1, "id" -> "4f62303b-bf36-48ec-9ed3-79d7d141e6c5", "pos" -> {-0.5349999999999999, 0.3918525382413843}|>, "39d3f441-bbfe-4cfe-8c6b-04f4e193bb67" -> <|"name" -> 4, "id" -> "39d3f441-bbfe-4cfe-8c6b-04f4e193bb67", "pos" -> {0.45500000000000007`, -0.34314746175861577`}|>, "bdbcd0d0-01d5-46b2-a071-3c742217ee58" -> <|"name" -> 6, "id" -> "bdbcd0d0-01d5-46b2-a071-3c742217ee58", "pos" -> {0.18000000000000016`, -0.5931474617586158}|>, "45c84633-138c-433b-bb0d-2475b7fdb6f8" -> <|"name" -> 2, "id" -> "45c84633-138c-433b-bb0d-2475b7fdb6f8", "pos" -> {0.09499999999999997, 0.4968525382413843}|>, "046655c3-593a-446e-a681-2e2d5ba7ef9d" -> <|"name" -> 5, "id" -> "046655c3-593a-446e-a681-2e2d5ba7ef9d", "pos" -> {1.665, 0.23185253824138424`}|>, "dd43662b-914b-47e0-94ab-a38c924f28ae" -> <|"name" -> 7, "id" -> "dd43662b-914b-47e0-94ab-a38c924f28ae", "pos" -> {-0.2699999999999999, 0.44685253824138427`}|>|>, "edge" -> <|"e3" -> <|"v1" -> "046655c3-593a-446e-a681-2e2d5ba7ef9d", "v2" -> "dd43662b-914b-47e0-94ab-a38c924f28ae", "id" -> "e3", "type" -> UndirectedEdge["v1", "v2"], "shape" -> Automatic|>|>, "selectedVertex" -> Null, "version" -> 2, "KeepVertexCoordinates" -> True, "CreateVertexSelects" -> True, "SnapToGrid" -> False, "IndexGraph" -> False, "PerformanceLimit" -> 450, "VertexLabels" -> None, "VertexSize" -> Small, "DirectedEdges" -> False, "ImageSize" -> {300, 148.63636363636365`}, "vCounter" -> 6, "eCounter" -> 1, "range" -> {{-0.8099999999999999, 1.94}, {-0.7293974617586158, 0.6331025382413843}}, "aspectRatio" -> 2.018348623853211, "inRangeQ" -> RegionMemberFunction[MeshRegion[{{-0.8099999999999999, -0.7293974617586158}, {1.94, -0.7293974617586158}, {1.94, 0.6331025382413843}, {-0.8099999999999999, 0.6331025382413843}}, {Polygon[{{1, 2, 3, 4}}]}, Method -> {"EliminateUnusedCoordinates" -> True, "DeleteDuplicateCoordinates" -> Automatic, "DeleteDuplicateCells" -> Automatic, "VertexAlias" -> Identity, "CheckOrientation" -> Automatic, "CoplanarityTolerance" -> Automatic, "CheckIntersections" -> Automatic, "BoundaryNesting" -> Automatic, "SeparateBoundaries" -> Automatic, "TJunction" -> Automatic, "PropagateMarkers" -> True, "ZeroTest" -> Automatic, "Hash" -> 6068137075315393170}], 2, Region`Mesh`CanonicalDistance[True]], "snap" -> False, "realVertexSize" -> 0.03419170071312608|>, IGraphM`GraphEditor`PackagePrivate`error$$ = False, IGraphM`GraphEditor`PackagePrivate`refresh$$}, *)
(*InterpretationBox[*)
(*PaneSelectorBox[{True->*)
(*ButtonBox[*)
(*DynamicBox[ToBoxes[IGraphM`GraphEditor`PackagePrivate`error$$, StandardForm]],*)
(*Appearance->Automatic,*)
(*BaseStyle->15,*)
(*ButtonFunction:>IGraphM`GraphEditor`PackagePrivate`refresh$$[],*)
(*Evaluator->Automatic,*)
(*Method->"Preemptive"], False->*)
(*PanelBox[*)
(*DynamicBox[ToBoxes[Refresh[IGraphM`GraphEditor`PackagePrivate`iGraphEditorPanel[Dynamic[IGraphM`GraphEditor`PackagePrivate`state$$]], None], StandardForm],*)
(*ImageSizeCache->{300., {72.25177550659812, 76.38458812976553}}],*)
(*BaseStyle->(CacheGraphics -> False),*)
(*FrameMargins->0]}, Dynamic[MatchQ[Blank[Failure]][IGraphM`GraphEditor`PackagePrivate`error$$]],*)
(*ImageSize->Automatic],*)
(*IGraphM`GraphEditor`PackagePrivate`GraphFromEditorState[IGraphM`GraphEditor`PackagePrivate`state$$]],*)
(*Deinitialization:>IGraphM`GraphEditor`PackagePrivate`iGraphEditorDeinitialization[IGraphM`GraphEditor`PackagePrivate`state$$],*)
(*DynamicModuleValues:>{{UpValues[IGraphM`GraphEditor`PackagePrivate`state$$] = {HoldPattern[Condition[IGraphM`GraphEditor`PackagePrivate`state$$[Pattern[IGraphM`PreciseTracking`PackagePrivate`key$, Blank[]], Pattern[IGraphM`PreciseTracking`PackagePrivate`rest$, BlankNullSequence[]]] = Pattern[IGraphM`PreciseTracking`PackagePrivate`val$, Blank[]], Not[TrueQ[IGraphM`PreciseTracking`PackagePrivate`$inside$7851]]]] :> Block[{IGraphM`PreciseTracking`PackagePrivate`$inside$7851 = True, IGraphM`PreciseTracking`PackagePrivate`$old = IGraphM`GraphEditor`PackagePrivate`state$$[IGraphM`PreciseTracking`PackagePrivate`key$, IGraphM`PreciseTracking`PackagePrivate`rest$]}, IGraphM`GraphEditor`PackagePrivate`state$$[IGraphM`PreciseTracking`PackagePrivate`key$, IGraphM`PreciseTracking`PackagePrivate`rest$] = IGraphM`PreciseTracking`PackagePrivate`val$; If[IGraphM`PreciseTracking`PackagePrivate`$old =!= IGraphM`PreciseTracking`PackagePrivate`val$, IGraphM`PreciseTracking`PackagePrivate`UpdateTarget[IGraphM`GraphEditor`PackagePrivate`state$$[IGraphM`PreciseTracking`PackagePrivate`key$]]]; IGraphM`PreciseTracking`PackagePrivate`val$]}}, {DownValues[IGraphM`GraphEditor`PackagePrivate`refresh$$] = {HoldPattern[IGraphM`GraphEditor`PackagePrivate`refresh$$[]] :> Module[{IGraphM`GraphEditor`PackagePrivate`temp$}, Catch[If[Needs["IGraphM`"] === $Failed, Throw[IGraphM`GraphEditor`PackagePrivate`error$$ = Failure["GraphEditor", <|"Message" -> "\[WarningSign] Failed to load IGraphM`, make sure it is installed."|>]]]; IGraphM`GraphEditor`PackagePrivate`temp$ = IGraphM`GraphEditor`PackagePrivate`geStateVersionCheck[IGraphM`GraphEditor`PackagePrivate`state$$]; If[AssociationQ[IGraphM`GraphEditor`PackagePrivate`temp$], IGraphM`GraphEditor`PackagePrivate`state$$ = IGraphM`GraphEditor`PackagePrivate`temp$, Throw[IGraphM`GraphEditor`PackagePrivate`error$$ = IGraphM`GraphEditor`PackagePrivate`temp$]]; IGraphM`GraphEditor`PackagePrivate`iGraphEditorInitialization[IGraphM`GraphEditor`PackagePrivate`state$$, IGraphM`GraphEditor`PackagePrivate`error$$]]]}}},*)
(*Initialization:>IGraphM`GraphEditor`PackagePrivate`refresh$$[]],*)
(*Setting[#, {0}]& ]\)*)


IGGraphEditor[ImageSize->100]


(* ::Text:: *)
(*- Is initial plot range adequate?*)


{ IGGraphEditor[IGGraphAtlas[123]],
  IGGraphEditor[Graph[{1, 2}, {}]],
  IGGraphEditor[Graph[{1, 2}, {1<->2}]] }


(* ::Subsubsection:: *)
(*Different input graphs*)


(* ::Text:: *)
(*Try to edit them*)


{ IGGraphEditor[IGEmptyGraph[]],
  IGGraphEditor[CompleteGraph[2, DirectedEdges -> True]],
  IGGraphEditor[CompleteGraph[2]],
  IGGraphEditor[IGEmptyGraph[1]] }


(* ::Subsubsection:: *)
(*Is it usable with larger graphs?*)


g = IGLayoutTutte@IGEdgeMap[#^1.5&, EdgeWeight]@IGDistanceWeighted@IGLayoutTutte@Graph@EdgeList@GraphData["GreatRhombicosidodecahedralGraph"];


IGGraphEditor[g]


IGGraphEditor[GraphData["GreatRhombicosidodecahedralGraph"], VertexSize->Tiny]


IGGraphEditor[ExampleData[{"NetworkGraph", "DolphinSocialNetwork"}], VertexSize->Tiny, ImageSize->555]


(* ::Subsubsection:: *)
(*Vertex size*)


IGGraphEditor[Graph@{1->2, 2->1}, VertexSize->#, ImageSize->200]& /@ {Tiny, Small, Medium, Large,0.2}


(* ::Subsubsection:: *)
(*Multi-edges*)


IGGraphEditor[
	Graph[{1->1, 1->1, 1->1, 2->2, 1->2, 1->2, 1->2, 2->3, 3->1}],
	VertexLabels -> "Name",
	ImageSize -> 800
]


(* ::Subsection:: *)
(*Notebook*)


IGGraphEditor[IGGraphAtlas[123]]


g = IGLayoutTutte@IGEdgeMap[#^1.5&, EdgeWeight]@IGDistanceWeighted@IGLayoutTutte@Graph@EdgeList@GraphData["GreatRhombicosidodecahedralGraph"]


IGGraphEditor[g]


IGGraphEditor[GraphData["GreatRhombicosidodecahedralGraph"], VertexSize->Tiny,ImageSize->700, VertexLabels->"Name"]


IGGraphEditor[Graph[{1->2,2->3,3->1}]]


Graph[{1->2,2->1,3->1}]//SimpleGraphQ


MultigraphQ


PreciseTracking`Private`$TrackedTargets
