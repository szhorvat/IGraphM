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
IGraphM`GraphEditor`PackagePrivate`$geDebug=True;
IGraphM`GraphEditor`PackagePrivate`$logTimings=False;
IGraphM`GraphEditor`PackagePrivate`$logDynamic=True;

Get["IGraphM`"]
IGraphM`GraphEditor`PackagePrivate`$gridLinesCount=25.;


TestReport@FileNameJoin[{NotebookDirectory[], "Editor.wl"}]


FileNameJoin[{NotebookDirectory[], "Editor.wl"}]//NotebookOpen


SetOptions[IGGraphEditor, ImageSize->Automatic];


IGGraphEditor[IGEmptyGraph[]]//ToBoxes


IGEmptyGraph[]//InputForm


(* ::Subsection:: *)
(*hands on tests*)


(* ::Subsubsection:: *)
(*Does it work at all?*)


(* ::Text:: *)
(*- Alt+Click build a graph*)
(*- Evaluate output to create a graph*)
(*    - pay attention, is the embedding the same after you drag a vertex and valuate the editor?*)
(*- Drag outside of range should extend the range @mouseUp*)
(*- Alt+click should not trigger orange resize frame, nor should discarding a potential edge.*)
(*- Click on vertex should not trigger update vertex position*)
(*- Hover over edge/vertex should thicken it.*)
(*- Are curved edges redrawn @mouseUp?*)


IGGraphEditor[
  Graph@{1->2},
  VertexSize -> Small, 
  "DirectedEdges"->True,
  "VertexColor" -> Red,
  "ShowSidePanel" -> True  
]


(* ::Input:: *)
(*InputForm@Normal@Graph[{1, 2, 3, 4}, {Null, {{1, 2}, {2, 3}, {3, 4}}}, {FormatType -> TraditionalForm, VertexCoordinates -> {{0.17812499999999998`, 0.6773712169117063}, {0.55625, 0.42112121691170634`}, {0.7843749999999999, 0.8117462169117063}, {0.7999999999999999, 0.38362121691170636`}}}]*)


Graph


IGGraphEditor[
  Graph@{1->2},  
  Prolog -> {
  Texture[img], 
  Polygon[Scaled/@{{0,0},{0,1},{1,1},{1,0}},VertexTextureCoordinates->{{0,0},{0,1},{1,1},{1,0}}
  ]
 }
]


IGGraphEditor[Graph@{1->2,2->1},VertexSize->Small, "DirectedEdges"->True]


IGraphM`PreciseTracking`PackagePrivate`$TrackedTargets


IGGraphEditor[ImageSize->100]


(* ::Text:: *)
(*- Is initial plot range adequate?*)


{ IGGraphEditor[IGGraphAtlas[123]],
  IGGraphEditor[Graph[{1, 2}, {}]],
  IGGraphEditor[Graph[{1, 2}, {1<->2}]],
  IGGraphEditor[
	Graph[{1->1, 1->1, 1->1, 2->2, 1->2, 1->2, 1->2, 2->3, 3->1}],
	VertexLabels -> "Name"	
] }


(* ::Subsubsection:: *)
(*Different input graphs*)


(* ::Text:: *)
(*Try to edit them*)


IGGraphEditor[CompleteGraph[2, DirectedEdges -> True]]


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


(* ::Subsubsection::Closed:: *)
(*Vertex size*)


IGGraphEditor[Graph@{1->2, 2->1}, VertexSize->#, ImageSize->200]& /@ {Tiny, Small, Medium, Large,0.2}


(* ::Subsubsection:: *)
(*Multi-edges*)


IGGraphEditor[
	Graph[{1->1, 1->1, 1->1, 2->2, 1->2, 1->2, 1->2, 2->3, 3->1}],
	VertexLabels -> "Name"	
]


(* ::Subsection:: *)
(*Notebook*)


IGGraphEditor[IGGraphAtlas[123]]


g = IGLayoutTutte@IGEdgeMap[#^1.5&, EdgeWeight]@IGDistanceWeighted@IGLayoutTutte@Graph@EdgeList@GraphData["GreatRhombicosidodecahedralGraph"]


IGGraphEditor[g]


IGGraphEditor[GraphData["GreatRhombicosidodecahedralGraph"], VertexSize->Tiny,ImageSize->700, VertexLabels->"Name"]


IGGraphEditor[Graph[{1->2,2->3,3->1}]]


Graph[{1->2,2->1,3->1}]//SimpleGraphQ
