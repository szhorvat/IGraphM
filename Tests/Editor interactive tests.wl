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
(*- Click on vertex should not trigger update vertex position*)
(*- Hover over edge/vertex should thicken it.*)


IGGraphEditor[VertexSize->Small, "SnapToGrid"->True]


IGraphM`PreciseTracking`PackagePrivate`$TrackedTargets


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


Names["*`$TrackedTargets"]


IGraphM`PreciseTracking`PackagePrivate`$TrackedTargets
