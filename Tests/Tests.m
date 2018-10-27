(* ::Package:: *)
(* MicroTest test file *)

(***** Utility functions *****)

tolEq[a_, b_, tol_ : 1*^-8 ] := Max@Abs[a-b] < tol

takeUpper[mat_?SquareMatrixQ] := Extract[mat, Subsets[Range@Length[mat], {2}]]
takeLower[mat_?SquareMatrixQ] := takeUpper@Transpose[mat]

sameGraphQ[g1_, g2_] :=
    Block[{UndirectedEdge},
      SetAttributes[UndirectedEdge, Orderless];
      Sort@VertexList[g1] === Sort@VertexList[g2] && Sort@EdgeList[g1] === Sort@EdgeList[g2]
    ]

(***** Test graphs *****)

SeedRandom[137]
ugs = RandomGraph[{10,20}];
dgs = RandomGraph[{10,30}, DirectedEdges->True];

ugi = GraphComputation`ToGraphRepresentation[RandomGraph[{12,25}], "Incidence"];
dgi = GraphComputation`ToGraphRepresentation[RandomGraph[{12,25}, DirectedEdges -> True], "Incidence"];

umulti = Graph[UndirectedEdge @@@ RandomInteger[10, {50, 2}]];
dmulti = Graph[DirectedEdge @@@ RandomInteger[10, {50, 2}]];

usweights = RandomReal[1, EdgeCount[ugs]];
wugs = Graph[ugs, EdgeWeight->usweights];

uiweights = RandomReal[1, EdgeCount[ugi]];
wugi = Graph[ugi, EdgeWeight->uiweights];

dsweights = RandomReal[1, EdgeCount[dgs]];
wdgs = Graph[dgs, EdgeWeight->dsweights];

diweights = RandomReal[1, EdgeCount[dgi]];
wdgi = Graph[dgi, EdgeWeight->diweights];

bidi = DirectedGraph[RandomGraph[{10, 20}]];

empty = Graph[{},{}];
edgeless = Graph[{1,2,3},{}];

ulist = Table[ GraphComputation`ToGraphRepresentation[RandomGraph[{n, 2n}], If[EvenQ[n], "Incidence", "Simple"]], {n, 5, 100, 5}];
dlist = Table[ GraphComputation`ToGraphRepresentation[RandomGraph[{n, 2n}, DirectedEdges -> True], If[EvenQ[n], "Incidence", "Simple"]], {n, 5, 100, 5}];

dolphin = ExampleData[{"NetworkGraph", "DolphinSocialNetwork"}];
web = ExampleData[{"NetworkGraph", "ExpandedComputationalGeometry"}];
collab = ExampleData[{"NetworkGraph", "CondensedMatterCollaborations2005"}];
football = ExampleData[{"NetworkGraph", "AmericanCollegeFootball"}];
lesmiserables = ExampleData[{"NetworkGraph", "LesMiserables"}];
terrorist  = ExampleData[{"NetworkGraph", "EastAfricaEmbassyAttacks"}];

bipartite = Graph[
  {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12},
  {3 <-> 6, 5 <-> 6, 4 <-> 7, 1 <-> 8, 2 <-> 8, 4 <-> 8, 1 <-> 9,
    2 <-> 9, 4 <-> 9, 5 <-> 10, 1 <-> 11, 4 <-> 11, 5 <-> 11, 2 <-> 12,
    5 <-> 12}
];

dbipartite = Graph[
  {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12},
  {3 -> 6, 4 -> 7, 1 -> 9, 5 -> 10, 2 -> 12, 3 -> 12, 4 -> 12, 9 -> 1,
    10 -> 1, 11 -> 1, 12 -> 2, 9 -> 3, 10 -> 3, 8 -> 5, 10 -> 5}
];

(*******************************************************************************)
MTSection["Sanity checks for Mathematica builtins"]

(* Does IsomorphicGraphQ work? *)
MT[
  IsomorphicGraphQ[Graph[{1<->2}], Graph[{"a"<->"b"}]],
  True
]

(* Does specifying edge list though vertex indices work? *)
MT[
  sameGraphQ[
    Graph[{"c","b","a"}, {"a"->"b", "b"->"c"}],
    Graph[{"c","b","a"}, {{3,2}, {2,1}}, DirectedEdges -> True]
  ],
  True
]

(* Verify GraphComputation`WeightValues *)

(* Check multigraphs *)
MT[
  GraphComputation`WeightValues@Graph[{1<->2, 1<->2}, EdgeWeight -> {3,4}],
  {3,4}
]

MT[
  GraphComputation`WeightValues[collab],
  PropertyValue[collab, EdgeWeight]
]

(* IGEdgeWeightedQ relies on the fact that this function returns a list of 1s for unweighted graphs *)
MT[
  GraphComputation`WeightValues[Graph[{1<->2, 2<->3}]],
  {1,1}
]

(* Verify that it returns 1s even when attempting to set a different default value *)
MT[
  GraphComputation`WeightValues[Graph[{1<->2, 2<->3}, Properties -> {"DefaultEdgeProperties" -> {EdgeWeight -> 3}}]],
  {1,1}
]

(* Similar checks for GraphComputation`WeightVector, which returns the vertex weight vector *)

MT[
  GraphComputation`WeightVector@Graph[{2, 1}, {1<->2, 1<->2}, EdgeWeight -> {3,4}, VertexWeight -> {5,6}],
  {5,6}
]

MT[
  GraphComputation`WeightVector[Graph[{1<->2, 2<->3}]],
  {1,1,1}
]

MT[
  GraphComputation`WeightVector[Graph[{1<->2, 2<->3}, Properties -> {"DefaultVertexProperties" -> {VertexWeight -> 3}}]],
  {1,1,1}
]

(* Tests for IncidenceMatrix *)

MT[
  Normal@IncidenceMatrix[Graph[{a, b}, {b -> a}]],
  {{1}, {-1}}
]

MT[
  Normal@IncidenceMatrix[Graph[{a, b}, {b <-> a}]],
  {{1}, {1}}
]

MT[
  Normal@IncidenceMatrix[Graph[{a, b}, {b <-> b}]],
  {{0}, {2}}
]

(* According to the documentation, the below should return {{0}, {-2}}, not {{0}, {2}}.
 * This test is to alert of any change in behaviour.
 * This is relevant for IG::fromIncidenceMatrix() in IG.h, however this function is already
 * written in a way that it should be immune to such a change.
 *)
MT[
  Normal@IncidenceMatrix[Graph[{a, b}, {b -> b}]],
  {{0}, {2}}
]

MT[
  Normal@IncidenceMatrix@Graph[Range[10], {1 <-> 8, 8 <-> 10, 6 <-> 10, 2 <-> 2, 7 <-> 1,
    8 <-> 9, 3 <-> 8, 2 <-> 7, 1 <-> 2, 3 <-> 3, 1 <-> 6, 6 <-> 5,
    3 <-> 2, 10 <-> 9, 10 <-> 9}],
  {{1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0},
    {0, 0, 0, 2, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0},
    {0, 0, 0, 0, 0, 0, 1, 0, 0, 2, 0, 0, 1, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
    {0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0},
    {0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
    {1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1},
    {0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1}}
]

MT[
  Normal@IncidenceMatrix@Graph[Range[10], {6 \[DirectedEdge] 9, 3 \[DirectedEdge] 7,
    3 \[DirectedEdge] 1, 7 \[DirectedEdge] 7, 3 \[DirectedEdge] 1,
    9 \[DirectedEdge] 2, 3 \[DirectedEdge] 4, 2 \[DirectedEdge] 9,
    10 \[DirectedEdge] 5, 10 \[DirectedEdge] 4, 6 \[DirectedEdge] 6,
    5 \[DirectedEdge] 5, 1 \[DirectedEdge] 5, 5 \[DirectedEdge] 9,
    10 \[DirectedEdge] 9}],
  {{0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0},
    {0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0},
    {0, -1, -1, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 2, 1, -1, 0},
    {-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0},
    {0, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {1, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0, 0, 1, 1},
    {0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, -1}}
]


(* The following tests are to verify that vertex properties are stored in a canonical form, regardless of in which order
   they were originally specified. This is important so that the return values of IGBlissCanonicalGraph are SameQ-comparable. *)
MT[
  AdjacencyGraph[{{0,1,0},{1,0,1},{0,1,0}}, Properties -> {1 -> {"foo" -> 1}, 2 -> {"foo" -> 2}, 3 -> {"foo" -> 3}}],
  AdjacencyGraph[{{0,1,0},{1,0,1},{0,1,0}}, Properties -> {3 -> {"foo" -> 3}, 1 -> {"foo" -> 1}, 2 -> {"foo" -> 2}}]
]

Module[{g},
  MT[
    g = AdjacencyGraph[{{0,1,0},{1,0,1},{0,1,0}}];
    PropertyValue[{g,1}, "foo"] = 1;
    PropertyValue[{g,2}, "foo"] = 2;
    PropertyValue[{g,3}, "foo"] = 3;
    g
    ,
    AdjacencyGraph[{{0,1,0},{1,0,1},{0,1,0}}, Properties -> {3 -> {"foo" -> 3}, 1 -> {"foo" -> 1}, 2 -> {"foo" -> 2}}]
  ]
]

Module[{g},
  MT[
    g = AdjacencyGraph[{{0,1,0},{1,0,1},{0,1,0}}];
    PropertyValue[{g,3}, "foo"] = 3;
    PropertyValue[{g,2}, "foo"] = 2;
    PropertyValue[{g,1}, "foo"] = 1;
    g
    ,
    AdjacencyGraph[{{0,1,0},{1,0,1},{0,1,0}}, Properties -> {3 -> {"foo" -> 3}, 1 -> {"foo" -> 1}, 2 -> {"foo" -> 2}}]
  ]
]


(* Verify that WeightedAdjacencyMatrix adds up the weights of parallel edges *)
MT[
  Normal@WeightedAdjacencyMatrix@Graph[{1 -> 2, 1 -> 2, 2 -> 3}, EdgeWeight -> {4, 5, 6}],
  {{0, 9, 0}, {0, 0, 6}, {0, 0, 0}}
]


(* Verify that even 0 weights are explicitly stored in a weighted adjacency matrix *)
MT[
  With[{sa = WeightedAdjacencyMatrix[Graph[{1 -> 2, 2 -> 1}, EdgeWeight -> {2, 0}]]},
    {sa["NonzeroPositions"], sa["NonzeroValues"], sa["Background"]}
  ],
  {{{1, 2}, {2, 1}}, {2, 0}, 0}
]


(*******************************************************************************)
MTSection["Basic"]

t = First@AbsoluteTiming[
  MT[
    <<IGraphM`,
    Null
  ]
];
Print["Package load timing: ", t]

t = First@AbsoluteTiming[
  IGraphM`LTemplate`ValidTemplateQ[IGraphM`IGraphM`PackagePrivate`template]
];
Print["Template verification timing: ", t]


IGraphM`LTemplate`UnloadTemplate[IGraphM`IGraphM`PackagePrivate`template]
t = First@AbsoluteTiming[
  IGraphM`IGraphM`PackagePrivate`LoadIGraphM[]
];
Print["Template load timing: ", t]


nameFromUsage[symname_] :=
    With[{sym = Symbol[symname]},
      First@StringCases[sym::usage, Shortest[name__] ~~ "[" :> name]]

MT[
  AllTrue[Complement[Names["IGraphM`*"], {"IGraphM", "IGMinSeparators", "$IGExportFormats"}], nameFromUsage[#] === # &],
  True
]

MT[
  MatchQ[IGVersion[], _String],
  True
]

If[$Notebooks,
  MT[
    IGDocumentation[]; NotebookClose /@ Select[Notebooks[], CurrentValue[#, WindowTitle] === "IGraph/M Documentation"&];,
    Null
  ]
]

(*******************************************************************************)
MTSection["Undirected"]

MT[
  ig = IGraphM`PackageScope`igMake[ugs],
  IGraphM`LTemplate`Classes`IG[1]
]

MT[
  ig@"directedQ"[],
  False
]

MT[
  ig@"weightedQ"[],
  False
]

MT[
  IGraphM`LTemplate`LExpressionList["IG"],
  {ig}
]

ig =.

MT[
  IGraphM`LTemplate`LExpressionList["IG"],
  {}
]


(*******************************************************************************)
MTSection["Directed"]

MT[
  ig = IGraphM`PackageScope`igMake[dgs],
  IGraphM`LTemplate`Classes`IG[2]
]

MT[
  ig@"directedQ"[],
  True
]

MT[
  ig@"weightedQ"[],
  False
]

MT[
  IGraphM`LTemplate`LExpressionList["IG"],
  {ig}
]

ig =.

MT[
  IGraphM`LTemplate`LExpressionList["IG"],
  {}
]


(*******************************************************************************)
MTSection["Weighted undirected"]

MT[
  ig = IGraphM`PackageScope`igMake[wugs],
  IGraphM`LTemplate`Classes`IG[3]
]

MT[
  ig@"directedQ"[],
  False
]

MT[
  ig@"weightedQ"[],
  True
]

MT[
  ig@"getWeights"[],
  usweights
]

MT[
  ig@"clearWeights"[];
  ig@"weightedQ"[],
  False
]

MT[
  IGraphM`LTemplate`LExpressionList["IG"],
  {ig}
]

ig =.

MT[
  IGraphM`LTemplate`LExpressionList["IG"],
  {}
]


(*******************************************************************************)
MTSection["Weighted directed"]

MT[
  ig = IGraphM`PackageScope`igMake[wdgs],
  IGraphM`LTemplate`Classes`IG[4]
]

MT[
  ig@"directedQ"[],
  True
]

MT[
  ig@"weightedQ"[],
  True
]

MT[
  ig@"getWeights"[],
  dsweights
]

MT[
  ig@"clearWeights"[];
  ig@"weightedQ"[],
  False
]

MT[
  IGraphM`LTemplate`LExpressionList["IG"],
  {ig}
]

ig =.

MT[
  IGraphM`LTemplate`LExpressionList["IG"],
  {}
]


(*******************************************************************************)
MTSection["Cycling grahps through igraph"]

compare[ig_, graph_] :=
    If[DirectedGraphQ[graph],
      IGraphM`PackageScope`igIndexVec[ig@"edgeList"[]] == List @@@ EdgeList@IndexGraph[graph]
      ,
      Sort /@ IGraphM`PackageScope`igIndexVec[ig@"edgeList"[]] == Sort /@ (List @@@ EdgeList@IndexGraph[graph])
    ]

cycleTestList = {
  ugs, dgs, ugi, dgi, umulti, dmulti,
  empty, edgeless,
  dolphin, web, football, collab,
  KaryTree[100], KaryTree[101, DirectedEdges -> True]
};

MT[
  compare[IGraphM`PackageScope`igMake[#], #],
  True
]& /@ cycleTestList

MT[
  compare[IGraphM`PackageScope`igMake[#], #],
  True
]& /@ Join[IGData[{"AllDirectedGraphs", 3}], IGData[{"AllDirectedGraphs", 4}]]

MT[
  compare[IGraphM`PackageScope`igMake[#], #],
  True
]& /@ Join[IGData[{"AllUndirectedGraphs", 3}], IGData[{"AllUndirectedGraphs", 4}]]

Print["Cycle timing: ", First@Timing[IGraphM`PackageScope`igMake /@ cycleTestList;] ]


(*******************************************************************************)
MTSection["Creation"]

(* IGGraphAtlas *)

MT[
  IGNullGraphQ@IGGraphAtlas[0],
  True
]

MT[
  IGIsomorphicQ[
    IGGraphAtlas[14],
    PathGraph@Range[4]
  ],
  True
]

(* IGKautzGraph *)

MT[
  VertexCount@IGKautzGraph[0, 0],
  1
]

MT[
  IGIsomorphicQ[
    IGKautzGraph[2, 2],
    Graph[{1 -> 5, 1 -> 6, 2 -> 7, 2 -> 8, 3 -> 9, 3 -> 10, 4 -> 11, 4 -> 12,
      5 -> 1, 5 -> 2, 6 -> 3, 6 -> 4, 7 -> 9, 7 -> 10, 8 -> 11, 8 -> 12,
      9 -> 1, 9 -> 2, 10 -> 3, 10 -> 4, 11 -> 5, 11 -> 6, 12 -> 7, 12 -> 8}]
  ],
  True
]

(* IGKaryTree *)

MT[
  IGIsomorphicQ[IGKaryTree[10, 3], KaryTree[10, 3]],
  True
]

MT[
  IGIsomorphicQ[
    IGKaryTree[9, 4, DirectedEdges -> True],
    KaryTree[9, 4, DirectedEdges -> True]
  ],
  True
]

(* IGCompleteGraph *)

MT[
  IGIsomorphicQ[IGCompleteGraph[4], CompleteGraph[4]],
  True
]

MT[
  IGIsomorphicQ[
    IGCompleteGraph[5, DirectedEdges -> True],
    CompleteGraph[5, DirectedEdges -> True]
  ],
  True
]

MT[
  IGIsomorphicQ[
    IGCompleteGraph[3, SelfLoops -> True],
    Graph[{1 <-> 1, 2 <-> 2, 3 <-> 3, 1 <-> 2, 2 <-> 3, 3 <-> 1}]
  ],
  True
]

(* IGCompleteAcyclicGraph *)

MT[
  IGIsomorphicQ[
    IGCompleteAcyclicGraph[4],
    DirectedGraph[CompleteGraph[4], "Acyclic"]
  ],
  True
]

(* IGDeBruijnGraph *)

MT[
  IGIsomorphicQ[
    IGDeBruijnGraph[3, 4],
    DeBruijnGraph[3, 4]
  ],
  True
]

(* IGChordalRing *)

MT[
  IGIsomorphicQ[
    IGChordalRing[15, {{3, 12, 4}, {7, 8, 11}}],
    Graph[{1 -> 2, 2 -> 3, 3 -> 4, 4 -> 5, 5 -> 6, 6 -> 7, 7 -> 8, 8 -> 9,
      9 -> 10, 10 -> 11, 11 -> 12, 12 -> 13, 13 -> 14, 14 -> 15, 1 -> 15,
      1 -> 4, 1 -> 8, 2 -> 14, 2 -> 10, 3 -> 7, 3 -> 14, 4 -> 7, 4 -> 11,
      2 -> 5, 5 -> 13, 6 -> 10, 2 -> 6, 7 -> 10, 7 -> 14, 5 -> 8, 1 -> 8,
      9 -> 13, 5 -> 9, 10 -> 13, 2 -> 10, 8 -> 11, 4 -> 11, 1 -> 12,
      8 -> 12, 1 -> 13, 5 -> 13, 11 -> 14, 7 -> 14, 4 -> 15, 11 -> 15}, DirectedEdges -> False]
  ],
  True
]

MT[
  IGChordalRing[15, {{3, 12, 4, 3}, {7, 8, 11, 3}}];,
  Null,
  {IGraphM::error}
]

(* IGEmptyGraph *)

MT[
  IGEmptyGraph[],
  IGEmptyGraph[0]
]

MT[
  IGNullGraphQ@IGEmptyGraph[],
  True
]

MT[
  IGIsomorphicQ[IGEmptyGraph[3], Graph[{1,2,3}, {}]],
  True
]

(* IGLCF *)

MT[
  With[{g = IGLCF[{-2, 2}, 2]},
    {EdgeCount[g], VertexCount[g]}
  ],
  {6, 4}
]

MT[
  With[{g = IGLCF[{2}, 2]},
    {EdgeCount[g], VertexCount[g]}
  ],
  {1, 2}
]

MT[
  IGIsomorphicQ[PetersenGraph[6, 2], IGLCF[{-5, 2, 4, -2, -5, 4, -4, 5, 2, -4, -2, 5}]],
  True
]

MT[
  IGNullGraphQ@IGLCF[{}],
  True
]


(* IGDegreeSequenceGame *)

MT[
  IGIsomorphicQ[
    IGDegreeSequenceGame[{0,0,0,0}, Method -> #],
    Graph[{1,2,3,4},{}]
  ],
  True
]& /{ "SimpleNoMultiple", "Simple" }

MT[
  IGIsomorphicQ[
    IGDegreeSequenceGame[{1,1,0,0}, Method -> #],
    Graph[{1,2,3,4},{1<->2}]
  ],
  True
]& /{ "SimpleNoMultiple", "Simple" }

MT[
  IGIsomorphicQ[
    IGDegreeSequenceGame[{2,2,2,2}, Method -> #],
    Graph[{1,2,3,4},{1<->2}]
  ],
  True
]& /{ "SimpleNoMultiple", "Simple", "VigerLatapy" }


(* IGKRegularGame *)

MT[
  IGIsomorphicQ[
    IGKRegularGame[0,0],
    empty
  ],
  True
]

MT[
  IGIsomorphicQ[
    IGKRegularGame[4,1],
    Graph[{1<->2, 3<->4}]
  ],
  True
]

MT[
  IGIsomorphicQ[
    IGKRegularGame[4,2],
    CycleGraph[4]
  ],
  True
]

MT[
  IGIsomorphicQ[
    IGKRegularGame[5,4],
    CompleteGraph[5]
  ],
  True
]

MT[
  IGKRegularGame[3, 3];,
  Null,
  {IGraphM::error}
]


(* IGBarabasiAlbertGame *)

MT[
  With[{g = IGBarabasiAlbertGame[100, 1]}, {VertexCount[g], EdgeCount[g]}],
  {100, 99}
]

MT[
  With[{g = IGBarabasiAlbertGame[100, 2]}, {VertexCount[g], EdgeCount[g]}],
  {100, 197}
]

(* IGGrowingGame *)

MT[
  With[{g = IGGrowingGame[100, 2]}, {VertexCount[g], EdgeCount[g]}],
  {100, 198}
]

(* IGGeometricGame *)

MT[
  VertexCount@IGGeometricGame[100, 0.1],
  100
]

(* IGMakeLattice *)

MT[
  IGIsomorphicQ[IGMakeLattice[{3,4}], GridGraph[{3,4}]],
  True
]

MT[
  PathGraphQ@IGMakeLattice[{5}],
  True
]

MT[
  IGIsomorphicQ[IGMakeLattice[{5}, "Periodic" -> True], CycleGraph[5]],
  True
]

MT[
  IGIsomorphicQ[
    IGMakeLattice[{9}, "Periodic" -> True, "Radius" -> 2],
    GraphPower[CycleGraph[9], 2]
  ],
  True
]

MT[
  IGIsomorphicQ[
    IGMakeLattice[{2, 3}, DirectedEdges -> True, "Mutual" -> True],
    DirectedGraph[GridGraph[{2, 3}]]
  ],
  True
]


(*******************************************************************************)
MTSection["Modification"]

(* IGConnectNeighborhood *)

MT[
  IGIsomorphicQ[
    IGConnectNeighborhood[Graph@{1 <-> 2, 2 <-> 3}],
    Graph[{1 <-> 2, 1 <-> 3, 2 <-> 3}]
  ],
  True
]

MT[
  IGIsomorphicQ[
    IGConnectNeighborhood[Graph@{1 <-> 2, 2 <-> 3, 1 <-> 1}],
    Graph[{1 <-> 1, 1 <-> 2, 1 <-> 3, 2 <-> 3}]
  ],
  True
]

MT[
  IGIsomorphicQ[
    IGConnectNeighborhood[Graph@{1 <-> 2, 2 <-> 3, 2 <-> 3}],
    Graph[{1 <-> 2, 1 <-> 3, 2 <-> 3, 2 <-> 3}]
  ],
  True
]

MT[
  IsomorphicGraphQ[
    IGConnectNeighborhood[CycleGraph[10], 2],
    IGMakeLattice[{10}, "Periodic" -> True, "Radius" -> 2]
  ],
  True
]

MT[
  EdgeList@IGConnectNeighborhood[empty, 2],
  {}
]

MT[
  IsomorphicGraphQ[
    IGConnectNeighborhood[ugs, 1],
    ugs
  ],
  True,
  {IGraphM::warning}
]

MT[
  sameGraphQ[
    IGConnectNeighborhood[
      Graph[{1 -> 2, 2 -> 1, 2 -> 3, 3 -> 2}],
      2
    ],
    Graph[{1 -> 2, 1 -> 3, 2 -> 1, 2 -> 3, 3 -> 1, 3 -> 2}]
  ],
  True
]

MT[
  sameGraphQ[
    IGConnectNeighborhood[
      Graph[{1 -> 2, 2 -> 3}],
      2
    ],
    Graph[{1 -> 2, 1 -> 3, 2 -> 3}]
  ],
  True
]

MT[
  sameGraphQ[
    IGConnectNeighborhood[
      Graph[{1 -> 2, 3 -> 2}],
      2
    ],
    Graph[{1 -> 2, 3 -> 2}]
  ],
  True
]


(*******************************************************************************)
MTSection["Rewiring"]

(* IGRewire *)

MT[
  VertexDegree@IGRewire[#, 100] == VertexDegree[#],
  True
] & /@ {ugs, ugi}

MT[
  With[{r = IGRewire[#, 100]},
    VertexInDegree[r] == VertexInDegree[#] &&
        VertexOutDegree[r] == VertexOutDegree[#]
  ],
  True
]& /@ {dgs, dgi}

MT[
  IsomorphicGraphQ[IGRewire[#, 0], #],
  True
]& /@ {ugs, dgs, ugi, dgi}


(* IGRewireEdges *)

MT[
  EdgeCount@IGRewireEdges[#, 0.1] == EdgeCount[#],
  True
]& /@ {ugs, dgs, ugi, dgi}

MT[
  IsomorphicGraphQ[IGRewireEdges[#, 0], #],
  True
]& /@ {ugs, dgs, ugi, dgi}

MT[
  IGIsomorphicQ[
    IGRewireEdges[empty, 0],
    empty
  ],
  True
]

MT[
  IGIsomorphicQ[
    IGRewireEdges[edgeless, 0],
    edgeless
  ],
  True
]


(*******************************************************************************)
MTSection["Acyclic graphs"]

cycGraph = Graph[{13 -> 7, 7 -> 0, 0 -> 16, 16 -> 2, 2 -> 15, 10 -> 5, 5 -> 12,
  12 -> 18, 18 -> 15, 17 -> 18, 15 -> 6, 6 -> 8, 8 -> 4, 9 -> 8,
  4 -> 19, 19 -> 11, 11 -> 1, 1 -> 20, 20 -> 3, 3 -> 4, 14 -> 19}];

MT[
  AcyclicGraphQ[EdgeDelete[cycGraph, IGFeedbackArcSet[cycGraph]]],
  True
]

MT[
  AllTrue[EdgeDelete[#, IGFeedbackArcSet[#]]& /@ ulist, AcyclicGraphQ],
  True
]

MT[
  AllTrue[EdgeDelete[#, IGFeedbackArcSet[#, Method -> "EadesLinSmyth"]]& /@ ulist, AcyclicGraphQ],
  True
]

MT[
  AllTrue[EdgeDelete[#, IGFeedbackArcSet[#]]& /@ Take[dlist, 6], AcyclicGraphQ],
  True
]

MT[
  AllTrue[EdgeDelete[#, IGFeedbackArcSet[#, Method -> "EadesLinSmyth"]]& /@ dlist, AcyclicGraphQ],
  True
]

MT[
  IGFeedbackArcSet[empty],
  {}
]

MT[
  IGFeedbackArcSet[edgeless],
  {}
]


(*******************************************************************************)
MTSection["Isomorphism"]

randomize[g_] := Graph[RandomSample@VertexList[g], EdgeList[g]]

grid = GridGraph[{3,3,3,3}];

isomorphismTests[if_] :=
    {
      MT[
        if[empty, empty],
        True
      ],
      MT[
        if[edgeless, edgeless],
        True
      ],
      MT[
        if[empty, edgeless],
        False
      ],
      MT[
        if[empty, ugs],
        False
      ],
      MT[
        if[empty, dgs],
        False
      ],
      MT[
        if[edgeless, randomize@edgeless],
        True
      ],
      MT[
        if[ugs, randomize[ugs]],
        True
      ],
      MT[
        if[ugi, randomize[ugi]],
        True
      ],
      MT[
        if[grid, randomize[grid]],
        True
      ]
    }

subisomorphismTests[if_] :=
    {
      MT[
        if[empty, empty],
        True
      ],
      MT[
        if[empty, edgeless],
        True
      ],
      MT[
        if[empty, ugs],
        True
      ],
      MT[
        if[empty, dgs],
        True
      ],
      MT[
        if[edgeless, dgs],
        True
      ],
      MT[
        if[edgeless, randomize@edgeless],
        True
      ],
      MT[
        if[ugs, randomize[ugs]],
        True
      ],
      MT[
        if[ugi, randomize[ugi]],
        True
      ],
      MT[
        if[grid, randomize[grid]],
        True
      ],
      MT[
        if[GridGraph[{3,3,3}], GridGraph[{4,4,4}]],
        True
      ],
      MT[
        if[Subgraph[ugs, RandomSample[VertexList[ugs], 5]], ugs],
        True
      ],
      MT[
        if[Subgraph[ugi, RandomSample[VertexList[ugi], 5]], ugi],
        True
      ]
    }


isomorphismTests /@ {
  IGIsomorphicQ,
  IGBlissIsomorphicQ,
  IGVF2IsomorphicQ
};

subisomorphismTests /@ {
  IGSubisomorphicQ,
  IGVF2SubisomorphicQ,
  IGLADSubisomorphicQ
};

directedIsomorphismTest[if_] := {
  MT[
    if[dgs, randomize@dgs],
    True
  ],
  MT[
    if[dgs, ugs],
    LibraryFunctionError["LIBRARY_FUNCTION_ERROR", 6],
    {IGraphM::error}
  ]
}

directedIsomorphismTest /@ { IGIsomorphicQ, IGSubisomorphicQ, IGVF2IsomorphicQ, IGVF2SubisomorphicQ  };


MT[
  Length@DeleteDuplicates[Rest@IGData[{"AllDirectedGraphs", 3}], IGIsomorphicQ],
  15
]

MT[
  Length@DeleteDuplicates[Rest@IGData[{"AllDirectedGraphs", 3}], IGVF2IsomorphicQ],
  15
]

MT[
  Length@DeleteDuplicates[Rest@IGData[{"AllDirectedGraphs", 3}], IGBlissIsomorphicQ],
  15
]

MT[
  Length@DeleteDuplicates[Rest@IGData[{"AllDirectedGraphs", 4}], IGIsomorphicQ],
  217
]

MT[
  Length@DeleteDuplicates[Rest@IGData[{"AllDirectedGraphs", 4}], IGVF2IsomorphicQ],
  217
]

MT[
  Length@DeleteDuplicates[Rest@IGData[{"AllDirectedGraphs", 4}], IGBlissIsomorphicQ],
  217
]

MT[
  IGBlissAutomorphismCount /@ IGData[{"AllDirectedGraphs", 4}],
  GroupOrder@*GraphAutomorphismGroup /@ IGData[{"AllDirectedGraphs", 4}]
]

MT[
  IGBlissAutomorphismCount /@ IGData[{"AllDirectedGraphs", 3}],
  GroupOrder@*GraphAutomorphismGroup /@ IGData[{"AllDirectedGraphs", 3}]
]

MT[
  IGBlissAutomorphismCount[ GraphData["PerkelGraph"] ],
  GraphData["PerkelGraph", "AutomorphismCount"]
]

MT[
  GroupOrder@PermutationGroup@IGBlissAutomorphismGroup[ GraphData["PerkelGraph"] ],
  GraphData["PerkelGraph", "AutomorphismCount"]
]

(* Self loop handling *)

MT[
  #[
      Graph[{1, 2, 3}, {1 <-> 2, 2 <-> 3, 3 <-> 3}],
      Graph[{1, 2, 3}, {1 <-> 2, 2 <-> 3}]
      ],
  False
]& /@ {IGIsomorphicQ, IGBlissIsomorphicQ, IGVF2IsomorphicQ, IGSubisomorphicQ, IGVF2SubisomorphicQ}

MT[
  #[
      Graph[{1, 2, 3}, {1 <-> 2, 2 <-> 3, 3 <-> 3}],
      Graph[{1, 2, 3}, {1 <-> 2, 2 <-> 3, 1 <-> 1}]
      ],
  True
]& /@ {IGIsomorphicQ, IGBlissIsomorphicQ, IGVF2IsomorphicQ, IGSubisomorphicQ, IGVF2SubisomorphicQ}

MT[
  #[
      Graph[{1, 2, 3}, {1 <-> 2, 2 <-> 3, 3 <-> 3}],
      Graph[{1, 2, 3}, {1 <-> 2, 2 <-> 3, 2 <-> 2}]
      ],
  False
]& /@ {IGIsomorphicQ, IGBlissIsomorphicQ, IGVF2IsomorphicQ, IGSubisomorphicQ, IGVF2SubisomorphicQ}


(* empty graph *)
MT[
  IGBlissCanonicalGraph[Graph[{},{}]],
  Graph[{},{}]
]

(* Ensure that directed graphs which are the same as their reversed graph are handled properly. *)
MT[
  With[{g = IGBlissCanonicalGraph[Graph[{1 -> 2, 2 -> 1}]]},
    {VertexList[g], EdgeList[g]}
  ],
  {{1, 2}, {1 \[DirectedEdge] 2, 2 \[DirectedEdge] 1}}
]

MT[
  With[{g = IGBlissCanonicalGraph[Graph[{2 -> 1, 1 -> 2}]]},
    {VertexList[g], EdgeList[g]}
  ],
  {{1, 2}, {1 \[DirectedEdge] 2, 2 \[DirectedEdge] 1}}
]

MT[
  With[{g = IGBlissCanonicalGraph[Graph[{b <-> a}]]},
    {VertexList[g], EdgeList[g]}
  ],
  {{1, 2}, {1 \[UndirectedEdge] 2}}
]


(*******************************************************************************)
MTSection["Isomorphism: coloured graphs"]

gvcol1 = {Graph[{1,2}, {1<->2}], "VertexColors" -> {1,2}};
gvcol2 = {Graph[{1,2}, {1<->2}], "VertexColors" -> {2,1}};
gvcol3 = {Graph[{1,2}, {1<->2}], "VertexColors" -> {3,4}};

MT[
  #[gvcol1, gvcol2],
  {<|1 -> 2, 2 -> 1|>}
]& /@ {IGVF2FindIsomorphisms, IGVF2FindSubisomorphisms, IGBlissGetIsomorphism, IGLADGetSubisomorphism, IGLADFindSubisomorphisms}

MT[
  #[gvcol1, gvcol2],
  1
]& /@ {IGVF2IsomorphismCount, IGVF2SubisomorphismCount, IGLADSubisomorphismCount}

MT[
  #[gvcol1, gvcol2],
  True
]& /@ {IGVF2IsomorphicQ, IGVF2SubisomorphicQ, IGBlissIsomorphicQ, IGLADSubisomorphicQ}

MT[
  #[gvcol1, gvcol3],
  False
]& /@ {IGVF2IsomorphicQ, IGVF2SubisomorphicQ, IGBlissIsomorphicQ, IGLADSubisomorphicQ}

MT[
  IGBlissAutomorphismCount[gvcol1],
  1
]

MT[
  IGBlissAutomorphismGroup[gvcol1],
  {}
]

MT[
  IGBlissCanonicalLabeling[gvcol1],
  <|1 -> 1, 2 -> 2|>
]

MT[
  IGBlissCanonicalLabeling[gvcol2],
  <|1 -> 2, 2 -> 1|>
]

MT[
  IGBlissCanonicalPermutation[gvcol1],
  {1, 2}
]

MT[
  IGBlissCanonicalPermutation[gvcol2],
  {2, 1}
]

MT[
  With[{g = IGBlissCanonicalGraph[{Graph[{b <-> a}], "VertexColors" -> <|b -> 3, a -> 1|>}]},
    {VertexList[g], EdgeList[g], IGVertexProp["Color"][g]}
  ],
  {{1, 2}, {1 \[UndirectedEdge] 2}, {1, 3}}
]

MT[
  With[{g = IGBlissCanonicalGraph[{Graph[{b <-> a}], "VertexColors" -> <|b -> 1, a -> 3|>}]},
    {VertexList[g], EdgeList[g], IGVertexProp["Color"][g]}
  ],
  {{1, 2}, {1 \[UndirectedEdge] 2}, {1, 3}}
]

MT[
  IGBlissCanonicalGraph[{Graph[{},{}], "VertexColors" -> {}}],
  Graph[{},{}]
]

MT[
  IGBlissCanonicalGraph[{Graph[{Property[1, "col" -> 5], Property[2, "col" -> 1]}, {}], "VertexColors" -> "col"}] // IGVertexProp["Color"],
  {1,5}
]

MT[
  IGBlissCanonicalGraph[{Graph[{Property[1, "col" -> 5], Property[2, "col" -> 1]}, {}], "VertexColors" -> None}] // IGVertexProp["Color"],
  {Missing["Nonexistent"], Missing["Nonexistent"]}
]



(*******************************************************************************)
MTSection["Isomorphism: multigraphs"]

gm1 = EdgeAdd[PathGraph@Range[5], 1 <-> 2];
gm2 = EdgeAdd[PathGraph@Range[5], 5 <-> 4];

MT[
  IGIsomorphicQ[gm1, gm2],
  True
]

MT[
  IGSubisomorphicQ[gm1, gm2],
  True
]

MT[
  #[gm1, gm2],
  $Failed,
  {IGraphM::blissnmg}
]& /@ {IGBlissGetIsomorphism, IGBlissIsomorphicQ}

MT[
  #[gm1],
  $Failed,
  {IGraphM::blissnmg}
]& /@ {IGBlissAutomorphismCount, IGBlissAutomorphismGroup, IGBlissCanonicalGraph, IGBlissCanonicalLabeling, IGBlissCanonicalPermutation}

MT[
  #[gm1, gm2],
  $Failed,
  {IGraphM::vf2nmg}
]& /@ {IGVF2FindIsomorphisms, IGVF2FindSubisomorphisms, IGVF2IsomorphicQ, IGVF2IsomorphismCount, IGVF2SubisomorphicQ, IGVF2SubisomorphismCount}

MT[
  #[gm1, gm2],
  LibraryFunctionError["LIBRARY_FUNCTION_ERROR", 6],
  {IGraphM::error}
]& /@ {IGLADFindSubisomorphisms, IGLADGetSubisomorphism, IGLADSubisomorphicQ}

(*******************************************************************************)
MTSection["Centralities"]

MT[
  #@IGEmptyGraph[],
  {}
]& /@ {IGBetweenness, IGEdgeBetweenness, IGCloseness, IGPageRank, IGEigenvectorCentrality, IGHubScore, IGAuthorityScore, IGConstraintScore}

MT[
  BetweennessCentrality[#] == IGBetweenness[#, Method -> "Fast"],
  True
]& /@ {ugs, ugi, dgs, dgi}

MT[
  BetweennessCentrality[#] ~tolEq~ IGBetweenness[#],
  True
]& /@ {ugs, ugi, dgs, dgi}


MT[
  EdgeBetweennessCentrality[#] == 2 IGEdgeBetweenness[#],
  True
]& /@ {ugs, ugi}

MT[
  EdgeBetweennessCentrality[#] == IGEdgeBetweenness[#],
  True
]& /@ {dgs, dgi}

MT[
  ClosenessCentrality[#] == IGCloseness[#, Normalized -> True],
  True
]& /@ {ugs, dgs, wugs, wdgs, umulti, dmulti}


MT[
  PageRankCentrality[#] ~tolEq~ IGPageRank[#],
  True
]& /@ {ugs, dgs}

(* check that edge multiplicity is taken to be equivalent with edge weights *)
MT[
  IGPageRank[Graph[{1 <-> 2, 2 <-> 3}, EdgeWeight -> {1, 2}]],
  IGPageRank[Graph[{1 <-> 2, 2 <-> 3, 2 <-> 3}]]
]

MT[
  IGPageRank[Graph[{1 <-> 2, 2 <-> 3}], Method -> "Arnoldi"] == IGPageRank[Graph[{1 <-> 2, 2 <-> 3}], Method -> "PRPACK"],
  True
]

MT[
  IGBetweenness[#],
  IGBetweennessEstimate[#, Infinity]
]& /@ {ugs, dgs}

MT[
  IGCloseness[#],
  IGClosenessEstimate[#, Infinity]
]& /@ {ugs, dgs}

MT[
  IGEdgeBetweenness[#],
  IGEdgeBetweennessEstimate[#, Infinity]
]& /@ {ugs, dgs}



(*******************************************************************************)
MTSection["Clustering coefficients"]

MT[
  IGGlobalClusteringCoefficient[#] == GlobalClusteringCoefficient[#],
  True
]& /@ ulist

MT[
  IGLocalClusteringCoefficient[#] == LocalClusteringCoefficient[#],
  True
]& /@ ulist

MT[
  IGAverageLocalClusteringCoefficient[#] == Mean@LocalClusteringCoefficient[#],
  True
]& /@ ulist

MT[
  GlobalClusteringCoefficient[umulti] == IGGlobalClusteringCoefficient[umulti],
  True
]

MT[
  LocalClusteringCoefficient[umulti] == IGLocalClusteringCoefficient[umulti],
  True
]

MT[
  Mean@LocalClusteringCoefficient[umulti] == IGAverageLocalClusteringCoefficient[umulti],
  True
]

MT[
  IGGlobalClusteringCoefficient[Graph[{1,2,3},{1<->2}]],
  0.
]

MT[
  IGGlobalClusteringCoefficient[Graph[{1,2,3},{1<->2}], "ExcludeIsolates" -> True],
  Indeterminate
]

MT[
  IGLocalClusteringCoefficient[Graph[{1,2,3,4}, {1<->2, 2<->3, 3<->1, 3<->4}]],
  N@{1, 1, 1/3, 0}
]

MT[
  IGLocalClusteringCoefficient[Graph[{1,2,3,4}, {1<->2, 2<->3, 3<->1, 3<->4}], "ExcludeIsolates" -> True],
  N@{1, 1, 1/3, Indeterminate}
]

MT[
  IGAverageLocalClusteringCoefficient[Graph[{1, 2, 3, 4}, {1 <-> 2, 2 <-> 3, 3 <-> 1, 3 <-> 4}]],
  N[7/12]
]

MT[
  IGAverageLocalClusteringCoefficient[Graph[{1, 2, 3, 4}, {1 <-> 2, 2 <-> 3, 3 <-> 1, 3 <-> 4}], "ExcludeIsolates" -> True],
  N[7/9]
]


(*******************************************************************************)
MTSection["Cliques"]

cl[g_, k_] :=
    Union[Sort /@ Values /@ IGLADFindSubisomorphisms[CompleteGraph[k], g]]

cl2[g_] :=
    Rest@Union[Sort /@ Catenate[Subsets /@ FindClique[g, Infinity, All]]]

canon[ls_] := Sort[Sort /@ ls]

rg = RandomGraph[{30,200}];

MT[
  canon[Join @@ Table[cl[rg, i], {i, 1, IGCliqueNumber[rg]}]] == canon@IGCliques[rg],
  True
]

MT[
  canon@FindClique[rg, Infinity, All] == canon@IGMaximalCliques[rg],
  True
]

MT[
  canon@IGMaximalCliques[rg, {IGCliqueNumber[rg]}] == canon@IGLargestCliques[rg],
  True
]

MT[
  Complement[canon@IGMaximalCliques[rg], canon@IGCliques[rg]],
  {}
]

MT[
  canon@IGCliques[rg] == canon@cl2[rg],
  True
]

MT[
  IGCliqueSizeCounts[rg],
  Values@KeySort@Counts[Length /@ IGCliques[rg]]
]

MT[
  IGMaximalCliqueSizeCounts[rg],
  With[{cl = FindClique[rg, Infinity, All]},
    Lookup[Counts[Length /@ cl], Range@Max[Length /@ cl], 0]
  ]
]

MT[
  IGCliques[empty],
  {}
]

MT[
  canon@IGCliques[Graph[{1,2,3},{}]],
  {{1},{2},{3}}
]

MT[
  IGMaximalCliques[empty],
  {}
]

MT[
  canon@IGMaximalCliques[Graph[{1,2,3},{}]],
  {{1},{2},{3}}
]

(* edge directions are ignored for clique finding *)
MT[
  canon@IGCliques[Graph[{1, 2}, {1 -> 2}]],
  {{1},{2},{1,2}},
  {IGraphM::warning}
]

(* edge directions are ignored for clique finding *)
MT[
  canon@IGMaximalCliques[Graph[{1, 2}, {1 -> 2}]],
  {{1,2}},
  {IGraphM::warning}
]

MT[
  IGCliqueNumber[#] == IGIndependenceNumber@GraphComplement[#],
  True
]& /@ {ugs, ugi, umulti}

MT[
  IGCliqueNumber@GraphComplement[#] == IGIndependenceNumber[#],
  True
]& /@ {ugs, ugi, umulti}

MT[
  IGCliqueNumber[#] == Length@First@FindClique[#],
  True
]& /@ {ugs, ugi, umulti}

MT[
  IGIndependenceNumber[#] == Length@First@FindIndependentVertexSet[#],
  True
]& /@ {ugs, ugi, umulti}

MT[
  IGCliqueNumber[empty],
  0
]

MT[
  IGIndependenceNumber[empty],
  0
]


(*******************************************************************************)
MTSection["Weighted cliques"]

(* weights must be integers *)
MT[
  canon@IGWeightedCliques[Graph[{1, 2, 3}, {1 <-> 2}, VertexWeight -> {1.5, 2, 3}], {0,Infinity}],
  canon[{{1}, {2}, {1, 2}, {3}}],
  {IGraphM::warning}
]

(* weights must be positive *)
MT[
  IGWeightedCliques[Graph[{1,2,3}, { 1<-> 2}, VertexWeight -> {0,1,2}], {0, Infinity}];,
  Null,
  {IGraphM::error}
]

(* edge directions are ignored *)
MT[
  canon@IGWeightedCliques[Graph[{1, 2}, {1 -> 2}, VertexWeight -> {1, 2}], {0, Infinity}],
  canon[{{1}, {2}, {1, 2}}],
  {IGraphM::warning}
]

vwg = Graph[{1 <-> 2, 1 <-> 3, 1 <-> 6, 2 <-> 3, 2 <-> 4, 2 <-> 6,
  2 <-> 9, 2 <-> 10, 3 <-> 5, 3 <-> 9, 4 <-> 6, 4 <-> 7, 4 <-> 10,
  5 <-> 9, 5 <-> 10, 6 <-> 7, 6 <-> 9, 7 <-> 8, 7 <-> 10, 8 <-> 10},
  VertexWeight -> {1, 6, 3, 1, 1, 6, 5, 8, 6, 9}];

MT[
  canon@IGWeightedCliques[vwg, {0, Infinity}],
  canon@{{6}, {4}, {6, 4}, {1}, {1, 6}, {3}, {1, 3}, {10}, {4,
    10}, {9}, {3, 9}, {6, 9}, {7}, {10, 7}, {4, 10, 7}, {4, 7}, {6, 4,
    7}, {6, 7}, {2}, {2, 9}, {2, 3, 9}, {2, 6, 9}, {2, 10}, {2, 4,
    10}, {2, 3}, {1, 2, 3}, {1, 2}, {1, 2, 6}, {2, 4}, {2, 6, 4}, {2,
    6}, {5}, {9, 5}, {3, 9, 5}, {10, 5}, {3, 5}, {8}, {7, 8}, {10, 7,
    8}, {10, 8}}
]

MT[
  canon@IGWeightedCliques[vwg, {5, Infinity}],
  canon@{{10}, {4, 10}, {9}, {3, 9}, {6, 9}, {7}, {10, 7}, {4, 10, 7}, {4,
    7}, {6, 4, 7}, {6, 7}, {2}, {2, 9}, {2, 3, 9}, {2, 6, 9}, {2,
    10}, {2, 4, 10}, {2, 3}, {1, 2, 3}, {1, 2}, {1, 2, 6}, {2, 4}, {2,
    6, 4}, {2, 6}, {5}, {9, 5}, {3, 9, 5}, {10, 5}, {3, 5}, {8}, {7,
    8}, {10, 7, 8}, {10, 8}}
]

MT[
  canon@IGWeightedCliques[vwg, {5, 6}],
  canon@{{10}, {4, 10}, {9}, {7}, {2}}
]

MT[
  IGWeightedCliques[vwg, {5, 5}],
  {{10}}
]

MT[
  canon@IGMaximalWeightedCliques[vwg, {0, Infinity}],
  canon@{{5, 10}, {1, 2, 3}, {1, 2, 6}, {2, 3, 9}, {2, 4, 6}, {2, 4, 10}, {2,
    6, 9}, {3, 5, 9}, {4, 6, 7}, {4, 7, 10}, {7, 8, 10}}
]

MT[
  canon@IGMaximalWeightedCliques[vwg, {15, Infinity}],
  canon@{{2, 3, 9}, {3, 5, 9}, {7, 8, 10}}
]

MT[
  canon@IGMaximalWeightedCliques[vwg, {15, 17}],
  canon@{{2, 3, 9}, {3, 5, 9}}
]

MT[
  canon@IGMaximalWeightedCliques[vwg, {15, 15}],
  canon@{{2, 3, 9}}
]

MT[
  IGWeightedCliqueNumber[vwg],
  20
]

MT[
  canon@IGLargestWeightedCliques[vwg],
  canon@{{10, 7, 8}}
]


(*******************************************************************************)
MTSection["Motifs and subgraph counts"]

triangleCount[am_?MatrixQ] := Tr[am.am.am]/6
triangleCount[g_?GraphQ] := triangleCount@AdjacencyMatrix[g]

MT[
  Length@IGTriangles[#],
  triangleCount[#]
] & /@ {ugs, ugi}

MT[
  IGTriangles[empty],
  {}
]

MT[
  IGTriangles[edgeless],
  {}
]

MT[
  IGMotifs[web, 3],
  {Indeterminate, Indeterminate, 42371, Indeterminate, 13453, 347,
    44442, 1250, 41, 316, 5, 3, 4, 13, 1, 0}
]

MT[
  IGMotifs[web, 4],
  {Indeterminate, Indeterminate, Indeterminate, 730640, Indeterminate,
    Indeterminate, Indeterminate, 192811, 6815, Indeterminate,
    Indeterminate, Indeterminate, 140888, 886073, 44066, Indeterminate,
    3856, 23919, 2805, 25591, 4056, 276, Indeterminate, Indeterminate,
    378328, 5261, 88, Indeterminate, Indeterminate, 12308, 137, 5404,
    111, Indeterminate, Indeterminate, 2048, 52, 365, 3, Indeterminate,
    505, 136668, 3396, 4074, 124, 2106, 243, 20333, 1245, 238, 32, 27,
    844, 13, 292, 5, 198, 3, 367, 8, 104, 3, Indeterminate, 51, 2, 0,
    2516, 208, 47, 10, 33, 0, 0, 14, 3, 0, 1351484, 37006, 818, 613, 525,
    42, 2720, 222, 20, 152, 1, 1, 0, 27, 3, 0, 15934, 47, 1, 86, 1034,
    18, 213, 4, 63, 1, 526, 696, 6, 45, 2, 3, 0, 0, 10, 0, 1, 0, 2, 11,
    0, 1, 0, 15, Indeterminate, 12, 1, 4, 0, 2, 1, 0, 2, 1, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 0, 1, 0, 0, 10, 2, 0, 2, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 6, 1, 0, 0, 0, 0, 0, 11, 0, 0, 0,
    0, 1, 0, 0, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
]

MT[
  IGMotifsTotalCount[web, 3],
  102246
]

MT[
  IGMotifsTotalCount[web, 4],
  4052036
]

MT[
  IGMotifs[dolphin, 3],
  {Indeterminate, Indeterminate, 638, 95}
]

MT[
  Length@IGTriangles[dolphin],
  95
]

MT[
  IGMotifs[dolphin, 4],
  {Indeterminate, Indeterminate, Indeterminate, Indeterminate, 709,
    Indeterminate, 2099, 768, 59, 138, 27}
]

MT[
  IGMotifs[empty, 3],
  {Indeterminate, Indeterminate, 0, 0}
]

MT[
  IGMotifs[empty, 4],
  {Indeterminate, Indeterminate, Indeterminate, Indeterminate, 0, Indeterminate, 0, 0, 0, 0, 0}
]

MT[
  IGMotifs[edgeless, 3],
  {Indeterminate, Indeterminate, 0, 0}
]

MT[
  IGMotifs[edgeless, 4],
  {Indeterminate, Indeterminate, Indeterminate, Indeterminate, 0, Indeterminate, 0, 0, 0, 0, 0}
]

MT[
  IGMotifsTotalCount[dolphin, 3],
  733
]

MT[
  IGMotifsTotalCount[dolphin, 4],
  3800
]

MT[
  IGMotifsTotalCount[empty, 3],
  0
]

MT[
  IGMotifsTotalCount[edgeless, 3],
  0
]

MT[
  IGTriadCensus[web],
  <|"003" -> 301658398, "012" -> 4601577, "102" -> 15179,
      "021D" -> 44442, "021U" -> 42371, "021C" -> 13453, "111D" -> 347,
      "111U" -> 316, "030T" -> 1250, "030C" -> 3, "201" -> 5, "120D" -> 41,
      "120U" -> 13, "120C" -> 4, "210" -> 1, "300" -> 0|>
]

MT[
  IGTriadCensus[dolphin],
  <|"003" -> 29108, "012" -> 0, "102" -> 7979, "021D" -> 0, "021U" -> 0,
      "021C" -> 0, "111D" -> 0, "111U" -> 0, "030T" -> 0, "030C" -> 0,
      "201" -> 638, "120D" -> 0, "120U" -> 0, "120C" -> 0, "210" -> 0,
      "300" -> 95|>,
  {IGraphM::warning}
]

MT[
  IGTriadCensus[empty],
  <|"003" -> 0, "012" -> 0, "102" -> 0, "021D" -> 0, "021U" -> 0,
    "021C" -> 0, "111D" -> 0, "111U" -> 0, "030T" -> 0, "030C" -> 0,
    "201" -> 0, "120D" -> 0, "120U" -> 0, "120C" -> 0, "210" -> 0,
    "300" -> 0|>,
  {IGraphM::warning}
]

MT[
  IGTriadCensus[Graph[{1,2,3},{}]],
  <|"003" -> 1, "012" -> 0, "102" -> 0, "021D" -> 0, "021U" -> 0,
    "021C" -> 0, "111D" -> 0, "111U" -> 0, "030T" -> 0, "030C" -> 0,
    "201" -> 0, "120D" -> 0, "120U" -> 0, "120C" -> 0, "210" -> 0,
    "300" -> 0|>,
  {IGraphM::warning}
]

MT[
  MapThread[SameQ, {IGTriadCensus[web] /@ Keys@IGData["MANTriadLabels"], IGMotifs[web, 3]}],
  {False, False, True, False, True, True, True, True, True, True, True, True, True, True, True, True}
]

MT[
  IGDyadCensus[empty],
  <|"Mutual" -> 0, "Asymmetric" -> 0, "Null" -> 0|>,
  {IGraphM::warning}
]

MT[
  IGDyadCensus[Graph[{1,2,3}, {}]],
  <|"Mutual" -> 0, "Asymmetric" -> 0, "Null" -> 3|>,
  {IGraphM::warning}
]

MT[
  IGDyadCensus[web],
  <|"Mutual" -> 13, "Asymmetric" -> 3927, "Null" -> 746985|>
]


(*******************************************************************************)
MTSection["Connectivity"]

MT[
  IGVertexConnectivity /@ ulist,
  VertexConnectivity /@ ulist
]

MT[
  IGVertexConnectivity /@ dlist,
  VertexConnectivity /@ dlist
]

MT[
  IGEdgeConnectivity /@ ulist,
  EdgeConnectivity /@ ulist
]

MT[
  IGEdgeConnectivity /@ dlist,
  EdgeConnectivity /@ dlist
]

MT[
  IGVertexConnectivity[ugs, 1, 5],
  VertexConnectivity[ugs, 1, 5]
]

MT[
  IGEdgeConnectivity[ugs, 1, 5],
  EdgeConnectivity[ugs, 1, 5]
]

MT[
  IGVertexConnectivity[dgs, 1, 6],
  VertexConnectivity[dgs, 1, 6]
]

MT[
  IGEdgeConnectivity[dgs, 1, 6],
  EdgeConnectivity[dgs, 1, 6]
]

MT[
  IGVertexConnectivity[empty],
  0
]

MT[
  IGVertexConnectivity[edgeless],
  0
]

MT[
  IGEdgeConnectivity[empty],
  0
]

MT[
  IGEdgeConnectivity[edgeless],
  0
]

MT[
  IGVertexConnectivity[umulti],
  VertexConnectivity[umulti]
]

MT[
  IGVertexConnectivity[dmulti],
  VertexConnectivity[dmulti]
]

MT[
  IGEdgeConnectivity[umulti],
  EdgeConnectivity[umulti]
]

MT[
  IGEdgeConnectivity[dmulti],
  EdgeConnectivity[dmulti]
]

MT[
  IGEdgeConnectivity[Graph[{2 <-> 1, 1 <-> 2, 1 <-> 3, 2 <-> 3, 3 <-> 1, 4 <-> 2}], 2, 3],
  3
]

MT[
  IGVertexConnectivity[ugs, 1, 2],
  LibraryFunctionError["LIBRARY_FUNCTION_ERROR", 6],
  {IGraphM::error}
]

MT[
  IGEdgeConnectivity[Graph[{1->2}], 1, 2],
  1
]

MT[
  IGEdgeConnectivity[Graph[{1->2}], 2, 1],
  0
]

MT[
  Sort[Sort /@ IGBiconnectedComponents[dolphin]],
  DeleteCases[Sort[Sort /@ KVertexConnectedComponents[dolphin, 2]], {_}]
]

MT[
  IGArticulationPoints[dolphin],
  {"Ripplefluke", "Scabs", "Patchback", "SN63", "Trigger", "Web", "Jet"}
]


(*******************************************************************************)
MTSection["Community detection"]

funs = {IGCommunitiesEdgeBetweenness, IGCommunitiesGreedy, IGCommunitiesInfoMAP,
  IGCommunitiesLabelPropagation, IGCommunitiesMultilevel,
  IGCommunitiesWalktrap, IGCommunitiesLeadingEigenvector, IGCommunitiesSpinGlass};

igClusterQ = Head[#] === IGClusterData&

MT[
  igClusterQ[#[dolphin]],
  True
]& /@ funs

MT[
  igClusterQ@IGCommunitiesOptimalModularity[lesmiserables],
  True
]

MT[
  c = IGCommunitiesEdgeBetweenness[dolphin];
  ToExpression@ToBoxes[c],
  c
]

MT[
  IGCompareCommunities[c, c],
  <|"VariationOfInformation" -> 0., "NormalizedMutualInformation" -> 1.,
    "SplitJoinDistance" -> 0, "UnadjustedRandIndex" -> 1.,
    "AdjustedRandIndex" -> 1.|>
]

MT[
  IGSeedRandom[42],
  Null
]

MT[
  c2 = IGCommunitiesWalktrap[dolphin];
  igClusterQ[c2],
  True
]

MT[
  IGCompareCommunities[c, c2],
  <|"VariationOfInformation" -> 0.8826360297713594`,
    "NormalizedMutualInformation" -> 0.6607128742019996`,
    "SplitJoinDistance" -> 21,
    "UnadjustedRandIndex" -> 0.8307773664727658`,
    "AdjustedRandIndex" -> 0.5847287340924315`|>
]

MT[
  IGCompareCommunities[c, c2],
  IGCompareCommunities[c2, c]
]


(*******************************************************************************)
MTSection["Shortest paths"]

(* verify handling of infinities *)
MT[
  IGDistanceMatrix[edgeless],
  { {0, Infinity, Infinity},
    {Infinity, 0, Infinity},
    {Infinity, Infinity, 0}}
]

MT[
  IGDistanceMatrix[Graph[{1, 2, 3}, {2 <-> 3}, EdgeWeight -> {1.2}]],
  {{0., Infinity, Infinity}, {Infinity, 0., 1.2}, {Infinity, 1.2, 0.}}
]

MT[
  IGDistanceMatrix[empty],
  {}
]

MT[
  IGDistanceCounts[empty],
  {}
]

MT[
  IGDistanceCounts[edgeless],
  {}
]

distanceCounts[graph_, verts_] :=
    Module[{idx, asc},
      idx = Lookup[
        AssociationThread[VertexList[graph], Range@VertexCount[graph]],
        verts
      ];
      asc = KeyDrop[{0, Infinity}]@Merge[
        {If[DirectedGraphQ[graph], 1, 1 / 2] Counts@Flatten[GraphDistanceMatrix[graph][[idx, idx]]],
          Counts@Flatten[GraphDistanceMatrix[graph][[idx, Complement[Range@VertexCount[graph], idx]]]]},
        Total
      ];
      Lookup[asc, Range@Max@Keys[asc], 0]
    ]

MT[
  IGDistanceCounts[#, All] == IGDistanceCounts[#, VertexList[#]],
  True
]& /@ {empty, edgeless, ugs, dgs, umulti, dmulti}

With[{v = RandomSample[VertexList[#], Quotient[VertexCount[#],3]]},
  MT[
    IGDistanceCounts[#, v],
    distanceCounts[#, v]
  ]
]& /@ {ugs, dgs, umulti, dmulti}

MT[
  IGDistanceHistogram[empty, 1],
  {}
]

MT[
  IGDistanceHistogram[edgeless, 1],
  {}
]

MT[
  IGAveragePathLength[empty],
  Indeterminate
]

MT[
  IGAveragePathLength[edgeless],
  Indeterminate
]

(* directed, unweighted, unconnected *)
MT[
  IGDistanceCounts[web],
  hist = Values@Rest@Most@KeySort@Counts@Flatten@GraphDistanceMatrix[web]
]

MT[
  IGAveragePathLength[web],
  N@Mean@WeightedData[Range@Length[hist], hist]
]

MT[
  {0} ~Join~ IGDistanceCounts[#] == IGDistanceHistogram[#, 1],
  True
] & /@ {dgs, dgi}

MT[
  {0} ~Join~ IGDistanceCounts[#] == IGDistanceHistogram[#, 1]/2,
  True
] & /@ {ugs, ugi}

MT[
  IGDistanceHistogram[wdgs, 0.1],
  With[{vals=takeUpper@GraphDistanceMatrix[wdgs] ~Join~ takeLower@GraphDistanceMatrix[wdgs]},
    BinCounts[vals, {0, Ceiling[Max[vals], 0.1], 0.1}]
  ]
]

MT[
  IGDistanceHistogram[wugs, 0.1],
  With[{vals=takeUpper@GraphDistanceMatrix[wugs] ~Join~ takeLower@GraphDistanceMatrix[wugs]},
    BinCounts[vals, {0, Ceiling[Max[vals], 0.1], 0.1}]
  ]
]


(* undirected, unweighted, connected *)
MT[
  IGDistanceCounts[dolphin],
  hist = Values@KeySort@Counts@takeUpper@GraphDistanceMatrix[dolphin]
]

MT[
  IGAveragePathLength[dolphin],
  N@Mean@WeightedData[Range@Length[hist], hist]
]


MT[
  IGGirth[PetersenGraph[5, 3]],
  5
]


MT[
  IGDistanceMatrix[ugs],
  GraphDistanceMatrix[ugs]
]

MT[
  IGDistanceMatrix[ugs, {1,2,4}],
  GraphDistanceMatrix[ugs][[{1,2,4}]]
]

MT[
  IGDistanceMatrix[ugs, {1,2,4}, {3,5}],
  GraphDistanceMatrix[ugs][[{1,2,4}, {3,5}]]
]

MT[
  IGDistanceMatrix[ugi],
  GraphDistanceMatrix[ugi]
]

MT[
  IGDistanceMatrix[ugi, {1,2,4}],
  GraphDistanceMatrix[ugi][[{1,2,4}]]
]

MT[
  IGDistanceMatrix[ugi, {1,2,4}, {3,5}],
  GraphDistanceMatrix[ugi][[{1,2,4}, {3,5}]]
]

MT[
  IGDistanceMatrix[dgs],
  GraphDistanceMatrix[dgs]
]

MT[
  IGDistanceMatrix[dgs, {1,2,4}],
  GraphDistanceMatrix[dgs][[{1,2,4}]]
]

MT[
  IGDistanceMatrix[dgs, {1,2,4}, {3,5}],
  GraphDistanceMatrix[dgs][[{1,2,4}, {3,5}]]
]

MT[
  IGDistanceMatrix[dgi],
  GraphDistanceMatrix[dgi]
]

MT[
  IGDistanceMatrix[dgi, {1,2,4}],
  GraphDistanceMatrix[dgi][[{1,2,4}]]
]

MT[
  IGDistanceMatrix[dgi, {1,2,4}, {3,5}],
  GraphDistanceMatrix[dgi][[{1,2,4}, {3,5}]]
]

MT[
  IGDistanceMatrix[wugs] == GraphDistanceMatrix[wugs],
  True
]

MT[
  IGDistanceMatrix[wugs, {1,2,4}] == GraphDistanceMatrix[wugs][[{1,2,4}]],
  True
]

MT[
  IGDistanceMatrix[wugs, {1,2,4}, {3,4}] == GraphDistanceMatrix[wugs][[{1,2,4}, {3,4}]],
  True
]

MT[
  IGDistanceMatrix[wugi] == GraphDistanceMatrix[wugi],
  True
]

MT[
  IGDistanceMatrix[wugi, {1,2,4}] == GraphDistanceMatrix[wugi][[{1,2,4}]],
  True
]

MT[
  IGDistanceMatrix[wugi, {1,2,4}, {3,4}] == GraphDistanceMatrix[wugi][[{1,2,4}, {3,4}]],
  True
]

MT[
  IGDistanceMatrix[wdgs] == GraphDistanceMatrix[wdgs],
  True
]

MT[
  IGDistanceMatrix[wdgs, {1,2,4}] == GraphDistanceMatrix[wdgs][[{1,2,4}]],
  True
]

MT[
  IGDistanceMatrix[wdgs, {1,2,4}, {3,4}] == GraphDistanceMatrix[wdgs][[{1,2,4}, {3,4}]],
  True
]

MT[
  IGDistanceMatrix[wdgi] == GraphDistanceMatrix[wdgi],
  True
]

MT[
  IGDistanceMatrix[wdgi, {1,2,4}] == GraphDistanceMatrix[wdgi][[{1,2,4}]],
  True
]

MT[
  IGDistanceMatrix[wdgi, {1,2,4}, {3,4}] == GraphDistanceMatrix[wdgi][[{1,2,4}, {3,4}]],
  True
]

MT[
  IGDiameter[empty],
  Infinity
]

MT[
  IGDiameter[edgeless],
  Infinity
]

MT[
  IGDiameter[#],
  GraphDiameter[#]
]& /{ugs, ugi, dgs, dgi, wugs, wugi, wdgs, wdgi, umulti, dmulti}

MT[
  IGFindDiameter[empty],
  {}
]

MT[
  IGFindDiameter[empty, "ByComponents" -> True],
  {}
]

MT[
  IGFindDiameter[edgeless],
  {}
]

MT[
  IGFindDiameter[Graph[{1, 2, 3}, {}], "ByComponents" -> True],
  {1}
]

MT[
  IGFindDiameter[PathGraph@Range[10]],
  Range[10]
]

MT[
  IGFindDiameter[Graph[{1 -> 2, 2 -> 3, 4 -> 3}]], (* "ByComponents" -> False is the default *)
  {}
]

MT[
  IGFindDiameter[Graph[{4,3,2,1}, {1 -> 2, 2 -> 3, 4 -> 3}], "ByComponents" -> True],
  {1,2,3}
]


(*******************************************************************************)
MTSection["Bipartite graphs"]

MT[
  IGBipartitePartitions[bipartite],
  {{1, 2, 3, 4, 5}, {6, 7, 8, 9, 10, 11, 12}}
]

MT[
  IGBipartitePartitions[empty],
  {{},{}}
]

(* not the only valid result, change test if implementation changes *)
MT[
  IGBipartitePartitions[Graph[{1,2,3,4}, {}]],
  {{1,2,3,4}, {}}
]

{
  MT[IGBipartiteQ[#], True],
  MT[IGBipartiteQ[#, IGBipartitePartitions[#]], True]
} & /@ Hold[
  bipartite, dbipartite,
  IGEmptyGraph[0], IGEmptyGraph[5],
  IGBipartiteGameGNM[10, 12, 33], IGBipartiteGameGNM[8, 17, 33, DirectedEdges -> True],
  IGBipartiteGameGNP[9, 6, 0.3], IGBipartiteGameGNP[21, 5, 0.6, DirectedEdges -> True]
] // ReleaseHold

MT[
  IGBipartitePartitions[dbipartite],
  {{1, 2, 3, 4, 5}, {6, 7, 8, 9, 10, 11, 12}}
]

bipartiteProjectionAM[g_, parts_] :=
    With[{im = IGBipartiteIncidenceMatrix[g, parts]},
      IGraphM`PackageScope`zeroDiagonal /@ Unitize[{im.Transpose[im], Transpose[im].im}]
    ]

MT[
  Equal[
    AdjacencyMatrix /@ IGBipartiteProjections[bipartite, IGBipartitePartitions[bipartite]],
    bipartiteProjectionAM[bipartite, IGBipartitePartitions[bipartite]]
  ],
  True
]


(*******************************************************************************)
MTSection["Chordal graphs"]

(* Chordal graphs *)

MT[
  IGChordalQ[empty],
  True
]

MT[
  IGChordalQ[edgeless],
  True
]

MT[
  IGChordalQ[CycleGraph[3]],
  True
]

MT[
  IGChordalQ[CycleGraph[4]],
  False
]


(*******************************************************************************)
MTSection["Mesh graphs"]

(* Mesh operations *)

MT[
  With[{g =
      IGMeshGraph@BoundaryDiscretizeRegion[Rectangle[], MaxCellMeasure -> Infinity]},
    {sameGraphQ[g, Graph[{1, 2, 3, 4}, {1 <-> 2, 2 <-> 3, 3 <-> 4, 4 <-> 1}]],
      GraphEmbedding[g], IGEdgeProp[EdgeWeight][g]}
  ]
  ,
  {True,
    {{0., 0.}, {1., 0.}, {1., 1.}, {0., 1.}},
    {1., 1., 1., 1.}}
]

MT[
  sameGraphQ[
    IGMeshGraph[DiscretizeRegion[Point[{0,0}]]],
    Graph[{1},{}]
  ],
  True
]

MT[
  sameGraphQ[
    IGMeshCellAdjacencyGraph[DiscretizeRegion[Point[{0,0}]], 0],
    Graph[{{0,1}},{}]
  ],
  True
]

MT[
  With[{r = DiscretizeRegion[Cuboid[], MaxCellMeasure -> 1]},
    Table[EdgeList@IGMeshCellAdjacencyGraph[r, d1, d2], {d1, 0, 3}, {d2, 0, 3}]
  ],
  {{{{0, 1} \[UndirectedEdge] {0, 2}, {0, 1} \[UndirectedEdge] {0, 3}, {0, 1} \[UndirectedEdge] {0, 4}, {0,
    1} \[UndirectedEdge] {0, 5}, {0, 1} \[UndirectedEdge] {0, 6}, {0, 1} \[UndirectedEdge] {0, 7}, {0,
    1} \[UndirectedEdge] {0, 8}, {0, 2} \[UndirectedEdge] {0, 3}, {0, 2} \[UndirectedEdge] {0, 6}, {0,
    2} \[UndirectedEdge] {0, 7}, {0, 3} \[UndirectedEdge] {0, 4}, {0, 3} \[UndirectedEdge] {0, 7}, {0,
    3} \[UndirectedEdge] {0, 8}, {0, 4} \[UndirectedEdge] {0, 8}, {0, 5} \[UndirectedEdge] {0, 6}, {0,
    5} \[UndirectedEdge] {0, 7}, {0, 5} \[UndirectedEdge] {0, 8}, {0, 6} \[UndirectedEdge] {0, 7}, {0,
    7} \[UndirectedEdge] {0, 8}}, {{0, 1} \[UndirectedEdge] {1, 1}, {0, 1} \[UndirectedEdge] {1, 3}, {0,
    1} \[UndirectedEdge] {1, 4}, {0, 1} \[UndirectedEdge] {1, 8}, {0, 1} \[UndirectedEdge] {1, 9}, {0,
    1} \[UndirectedEdge] {1, 12}, {0, 1} \[UndirectedEdge] {1, 14}, {0, 2} \[UndirectedEdge] {1, 14}, {0,
    2} \[UndirectedEdge] {1, 15}, {0, 2} \[UndirectedEdge] {1, 16}, {0, 2} \[UndirectedEdge] {1, 17}, {0,
    3} \[UndirectedEdge] {1, 7}, {0, 3} \[UndirectedEdge] {1, 9}, {0, 3} \[UndirectedEdge] {1, 10}, {0,
    3} \[UndirectedEdge] {1, 13}, {0, 3} \[UndirectedEdge] {1, 16}, {0, 4} \[UndirectedEdge] {1, 10}, {0,
    4} \[UndirectedEdge] {1, 11}, {0, 4} \[UndirectedEdge] {1, 12}, {0, 5} \[UndirectedEdge] {1, 4}, {0,
    5} \[UndirectedEdge] {1, 5}, {0, 5} \[UndirectedEdge] {1, 6}, {0, 5} \[UndirectedEdge] {1, 19}, {0,
    6} \[UndirectedEdge] {1, 2}, {0, 6} \[UndirectedEdge] {1, 3}, {0, 6} \[UndirectedEdge] {1, 6}, {0,
    6} \[UndirectedEdge] {1, 17}, {0, 7} \[UndirectedEdge] {1, 1}, {0, 7} \[UndirectedEdge] {1, 2}, {0,
    7} \[UndirectedEdge] {1, 5}, {0, 7} \[UndirectedEdge] {1, 13}, {0, 7} \[UndirectedEdge] {1, 15}, {0,
    7} \[UndirectedEdge] {1, 18}, {0, 8} \[UndirectedEdge] {1, 7}, {0, 8} \[UndirectedEdge] {1, 8}, {0,
    8} \[UndirectedEdge] {1, 11}, {0, 8} \[UndirectedEdge] {1, 18}, {0, 8} \[UndirectedEdge] {1, 19}}, {{0,
    1} \[UndirectedEdge] {2, 1}, {0, 1} \[UndirectedEdge] {2, 2}, {0, 1} \[UndirectedEdge] {2, 3}, {0,
    1} \[UndirectedEdge] {2, 5}, {0, 1} \[UndirectedEdge] {2, 7}, {0, 1} \[UndirectedEdge] {2, 8}, {0,
    1} \[UndirectedEdge] {2, 9}, {0, 1} \[UndirectedEdge] {2, 10}, {0, 1} \[UndirectedEdge] {2, 11}, {0,
    1} \[UndirectedEdge] {2, 14}, {0, 1} \[UndirectedEdge] {2, 16}, {0, 1} \[UndirectedEdge] {2, 17}, {0,
    2} \[UndirectedEdge] {2, 10}, {0, 2} \[UndirectedEdge] {2, 11}, {0, 2} \[UndirectedEdge] {2, 12}, {0,
    2} \[UndirectedEdge] {2, 13}, {0, 2} \[UndirectedEdge] {2, 14}, {0, 3} \[UndirectedEdge] {2, 5}, {0,
    3} \[UndirectedEdge] {2, 6}, {0, 3} \[UndirectedEdge] {2, 7}, {0, 3} \[UndirectedEdge] {2, 9}, {0,
    3} \[UndirectedEdge] {2, 11}, {0, 3} \[UndirectedEdge] {2, 12}, {0, 3} \[UndirectedEdge] {2, 15}, {0,
    4} \[UndirectedEdge] {2, 6}, {0, 4} \[UndirectedEdge] {2, 7}, {0, 4} \[UndirectedEdge] {2, 8}, {0,
    5} \[UndirectedEdge] {2, 2}, {0, 5} \[UndirectedEdge] {2, 3}, {0, 5} \[UndirectedEdge] {2, 4}, {0,
    5} \[UndirectedEdge] {2, 17}, {0, 5} \[UndirectedEdge] {2, 18}, {0, 6} \[UndirectedEdge] {2, 1}, {0,
    6} \[UndirectedEdge] {2, 3}, {0, 6} \[UndirectedEdge] {2, 4}, {0, 6} \[UndirectedEdge] {2, 13}, {0,
    6} \[UndirectedEdge] {2, 14}, {0, 7} \[UndirectedEdge] {2, 1}, {0, 7} \[UndirectedEdge] {2, 2}, {0,
    7} \[UndirectedEdge] {2, 4}, {0, 7} \[UndirectedEdge] {2, 9}, {0, 7} \[UndirectedEdge] {2, 10}, {0,
    7} \[UndirectedEdge] {2, 12}, {0, 7} \[UndirectedEdge] {2, 13}, {0, 7} \[UndirectedEdge] {2, 15}, {0,
    7} \[UndirectedEdge] {2, 16}, {0, 7} \[UndirectedEdge] {2, 18}, {0, 8} \[UndirectedEdge] {2, 5}, {0,
    8} \[UndirectedEdge] {2, 6}, {0, 8} \[UndirectedEdge] {2, 8}, {0, 8} \[UndirectedEdge] {2, 15}, {0,
    8} \[UndirectedEdge] {2, 16}, {0, 8} \[UndirectedEdge] {2, 17}, {0, 8} \[UndirectedEdge] {2, 18}}, {{0,
    1} \[UndirectedEdge] {3, 1}, {0, 1} \[UndirectedEdge] {3, 2}, {0, 1} \[UndirectedEdge] {3, 3}, {0,
    1} \[UndirectedEdge] {3, 4}, {0, 1} \[UndirectedEdge] {3, 5}, {0, 1} \[UndirectedEdge] {3, 6}, {0,
    2} \[UndirectedEdge] {3, 3}, {0, 2} \[UndirectedEdge] {3, 4}, {0, 3} \[UndirectedEdge] {3, 2}, {0,
    3} \[UndirectedEdge] {3, 3}, {0, 3} \[UndirectedEdge] {3, 5}, {0, 4} \[UndirectedEdge] {3, 2}, {0,
    5} \[UndirectedEdge] {3, 1}, {0, 5} \[UndirectedEdge] {3, 6}, {0, 6} \[UndirectedEdge] {3, 1}, {0,
    6} \[UndirectedEdge] {3, 4}, {0, 7} \[UndirectedEdge] {3, 1}, {0, 7} \[UndirectedEdge] {3, 3}, {0,
    7} \[UndirectedEdge] {3, 4}, {0, 7} \[UndirectedEdge] {3, 5}, {0, 7} \[UndirectedEdge] {3, 6}, {0,
    8} \[UndirectedEdge] {3, 2}, {0, 8} \[UndirectedEdge] {3, 5}, {0, 8} \[UndirectedEdge] {3, 6}}}, {{{1,
    1} \[UndirectedEdge] {0, 1}, {1, 1} \[UndirectedEdge] {0, 7}, {1, 2} \[UndirectedEdge] {0, 6}, {1,
    2} \[UndirectedEdge] {0, 7}, {1, 3} \[UndirectedEdge] {0, 1}, {1, 3} \[UndirectedEdge] {0, 6}, {1,
    4} \[UndirectedEdge] {0, 1}, {1, 4} \[UndirectedEdge] {0, 5}, {1, 5} \[UndirectedEdge] {0, 5}, {1,
    5} \[UndirectedEdge] {0, 7}, {1, 6} \[UndirectedEdge] {0, 5}, {1, 6} \[UndirectedEdge] {0, 6}, {1,
    7} \[UndirectedEdge] {0, 3}, {1, 7} \[UndirectedEdge] {0, 8}, {1, 8} \[UndirectedEdge] {0, 1}, {1,
    8} \[UndirectedEdge] {0, 8}, {1, 9} \[UndirectedEdge] {0, 1}, {1, 9} \[UndirectedEdge] {0, 3}, {1,
    10} \[UndirectedEdge] {0, 3}, {1, 10} \[UndirectedEdge] {0, 4}, {1, 11} \[UndirectedEdge] {0, 4}, {1,
    11} \[UndirectedEdge] {0, 8}, {1, 12} \[UndirectedEdge] {0, 1}, {1, 12} \[UndirectedEdge] {0, 4}, {1,
    13} \[UndirectedEdge] {0, 3}, {1, 13} \[UndirectedEdge] {0, 7}, {1, 14} \[UndirectedEdge] {0, 1}, {1,
    14} \[UndirectedEdge] {0, 2}, {1, 15} \[UndirectedEdge] {0, 2}, {1, 15} \[UndirectedEdge] {0, 7}, {1,
    16} \[UndirectedEdge] {0, 2}, {1, 16} \[UndirectedEdge] {0, 3}, {1, 17} \[UndirectedEdge] {0, 2}, {1,
    17} \[UndirectedEdge] {0, 6}, {1, 18} \[UndirectedEdge] {0, 7}, {1, 18} \[UndirectedEdge] {0, 8}, {1,
    19} \[UndirectedEdge] {0, 5}, {1, 19} \[UndirectedEdge] {0, 8}}, {{1, 1} \[UndirectedEdge] {1, 2}, {1,
    1} \[UndirectedEdge] {1, 3}, {1, 1} \[UndirectedEdge] {1, 4}, {1, 1} \[UndirectedEdge] {1, 5}, {1,
    1} \[UndirectedEdge] {1, 8}, {1, 1} \[UndirectedEdge] {1, 9}, {1, 1} \[UndirectedEdge] {1, 12}, {1,
    1} \[UndirectedEdge] {1, 13}, {1, 1} \[UndirectedEdge] {1, 14}, {1, 1} \[UndirectedEdge] {1, 15}, {1,
    1} \[UndirectedEdge] {1, 18}, {1, 2} \[UndirectedEdge] {1, 3}, {1, 2} \[UndirectedEdge] {1, 5}, {1,
    2} \[UndirectedEdge] {1, 6}, {1, 2} \[UndirectedEdge] {1, 13}, {1, 2} \[UndirectedEdge] {1, 15}, {1,
    2} \[UndirectedEdge] {1, 17}, {1, 2} \[UndirectedEdge] {1, 18}, {1, 3} \[UndirectedEdge] {1, 4}, {1,
    3} \[UndirectedEdge] {1, 6}, {1, 3} \[UndirectedEdge] {1, 8}, {1, 3} \[UndirectedEdge] {1, 9}, {1,
    3} \[UndirectedEdge] {1, 12}, {1, 3} \[UndirectedEdge] {1, 14}, {1, 3} \[UndirectedEdge] {1, 17}, {1,
    4} \[UndirectedEdge] {1, 5}, {1, 4} \[UndirectedEdge] {1, 6}, {1, 4} \[UndirectedEdge] {1, 8}, {1,
    4} \[UndirectedEdge] {1, 9}, {1, 4} \[UndirectedEdge] {1, 12}, {1, 4} \[UndirectedEdge] {1, 14}, {1,
    4} \[UndirectedEdge] {1, 19}, {1, 5} \[UndirectedEdge] {1, 6}, {1, 5} \[UndirectedEdge] {1, 13}, {1,
    5} \[UndirectedEdge] {1, 15}, {1, 5} \[UndirectedEdge] {1, 18}, {1, 5} \[UndirectedEdge] {1, 19}, {1,
    6} \[UndirectedEdge] {1, 17}, {1, 6} \[UndirectedEdge] {1, 19}, {1, 7} \[UndirectedEdge] {1, 8}, {1,
    7} \[UndirectedEdge] {1, 9}, {1, 7} \[UndirectedEdge] {1, 10}, {1, 7} \[UndirectedEdge] {1, 11}, {1,
    7} \[UndirectedEdge] {1, 13}, {1, 7} \[UndirectedEdge] {1, 16}, {1, 7} \[UndirectedEdge] {1, 18}, {1,
    7} \[UndirectedEdge] {1, 19}, {1, 8} \[UndirectedEdge] {1, 9}, {1, 8} \[UndirectedEdge] {1, 11}, {1,
    8} \[UndirectedEdge] {1, 12}, {1, 8} \[UndirectedEdge] {1, 14}, {1, 8} \[UndirectedEdge] {1, 18}, {1,
    8} \[UndirectedEdge] {1, 19}, {1, 9} \[UndirectedEdge] {1, 10}, {1, 9} \[UndirectedEdge] {1, 12}, {1,
    9} \[UndirectedEdge] {1, 13}, {1, 9} \[UndirectedEdge] {1, 14}, {1, 9} \[UndirectedEdge] {1, 16}, {1,
    10} \[UndirectedEdge] {1, 11}, {1, 10} \[UndirectedEdge] {1, 12}, {1, 10} \[UndirectedEdge] {1, 13}, {1,
    10} \[UndirectedEdge] {1, 16}, {1, 11} \[UndirectedEdge] {1, 12}, {1, 11} \[UndirectedEdge] {1, 18}, {1,
    11} \[UndirectedEdge] {1, 19}, {1, 12} \[UndirectedEdge] {1, 14}, {1, 13} \[UndirectedEdge] {1, 15}, {1,
    13} \[UndirectedEdge] {1, 16}, {1, 13} \[UndirectedEdge] {1, 18}, {1, 14} \[UndirectedEdge] {1, 15}, {1,
    14} \[UndirectedEdge] {1, 16}, {1, 14} \[UndirectedEdge] {1, 17}, {1, 15} \[UndirectedEdge] {1, 16}, {1,
    15} \[UndirectedEdge] {1, 17}, {1, 15} \[UndirectedEdge] {1, 18}, {1, 16} \[UndirectedEdge] {1, 17}, {1,
    18} \[UndirectedEdge] {1, 19}}, {{1, 1} \[UndirectedEdge] {2, 1}, {1, 1} \[UndirectedEdge] {2, 2}, {1,
    1} \[UndirectedEdge] {2, 9}, {1, 1} \[UndirectedEdge] {2, 10}, {1, 1} \[UndirectedEdge] {2, 16}, {1,
    2} \[UndirectedEdge] {2, 1}, {1, 2} \[UndirectedEdge] {2, 4}, {1, 2} \[UndirectedEdge] {2, 13}, {1,
    3} \[UndirectedEdge] {2, 1}, {1, 3} \[UndirectedEdge] {2, 3}, {1, 3} \[UndirectedEdge] {2, 14}, {1,
    4} \[UndirectedEdge] {2, 2}, {1, 4} \[UndirectedEdge] {2, 3}, {1, 4} \[UndirectedEdge] {2, 17}, {1,
    5} \[UndirectedEdge] {2, 2}, {1, 5} \[UndirectedEdge] {2, 4}, {1, 5} \[UndirectedEdge] {2, 18}, {1,
    6} \[UndirectedEdge] {2, 3}, {1, 6} \[UndirectedEdge] {2, 4}, {1, 7} \[UndirectedEdge] {2, 5}, {1,
    7} \[UndirectedEdge] {2, 6}, {1, 7} \[UndirectedEdge] {2, 15}, {1, 8} \[UndirectedEdge] {2, 5}, {1,
    8} \[UndirectedEdge] {2, 8}, {1, 8} \[UndirectedEdge] {2, 16}, {1, 8} \[UndirectedEdge] {2, 17}, {1,
    9} \[UndirectedEdge] {2, 5}, {1, 9} \[UndirectedEdge] {2, 7}, {1, 9} \[UndirectedEdge] {2, 9}, {1,
    9} \[UndirectedEdge] {2, 11}, {1, 10} \[UndirectedEdge] {2, 6}, {1, 10} \[UndirectedEdge] {2, 7}, {1,
    11} \[UndirectedEdge] {2, 6}, {1, 11} \[UndirectedEdge] {2, 8}, {1, 12} \[UndirectedEdge] {2, 7}, {1,
    12} \[UndirectedEdge] {2, 8}, {1, 13} \[UndirectedEdge] {2, 9}, {1, 13} \[UndirectedEdge] {2, 12}, {1,
    13} \[UndirectedEdge] {2, 15}, {1, 14} \[UndirectedEdge] {2, 10}, {1, 14} \[UndirectedEdge] {2, 11}, {1,
    14} \[UndirectedEdge] {2, 14}, {1, 15} \[UndirectedEdge] {2, 10}, {1, 15} \[UndirectedEdge] {2, 12}, {1,
    15} \[UndirectedEdge] {2, 13}, {1, 16} \[UndirectedEdge] {2, 11}, {1, 16} \[UndirectedEdge] {2, 12}, {1,
    17} \[UndirectedEdge] {2, 13}, {1, 17} \[UndirectedEdge] {2, 14}, {1, 18} \[UndirectedEdge] {2, 15}, {1,
    18} \[UndirectedEdge] {2, 16}, {1, 18} \[UndirectedEdge] {2, 18}, {1, 19} \[UndirectedEdge] {2, 17}, {1,
    19} \[UndirectedEdge] {2, 18}}, {{1, 1} \[UndirectedEdge] {3, 1}, {1, 1} \[UndirectedEdge] {3, 3}, {1,
    1} \[UndirectedEdge] {3, 4}, {1, 1} \[UndirectedEdge] {3, 5}, {1, 1} \[UndirectedEdge] {3, 6}, {1,
    2} \[UndirectedEdge] {3, 1}, {1, 2} \[UndirectedEdge] {3, 4}, {1, 3} \[UndirectedEdge] {3, 1}, {1,
    3} \[UndirectedEdge] {3, 4}, {1, 4} \[UndirectedEdge] {3, 1}, {1, 4} \[UndirectedEdge] {3, 6}, {1,
    5} \[UndirectedEdge] {3, 1}, {1, 5} \[UndirectedEdge] {3, 6}, {1, 6} \[UndirectedEdge] {3, 1}, {1,
    7} \[UndirectedEdge] {3, 2}, {1, 7} \[UndirectedEdge] {3, 5}, {1, 8} \[UndirectedEdge] {3, 2}, {1,
    8} \[UndirectedEdge] {3, 5}, {1, 8} \[UndirectedEdge] {3, 6}, {1, 9} \[UndirectedEdge] {3, 2}, {1,
    9} \[UndirectedEdge] {3, 3}, {1, 9} \[UndirectedEdge] {3, 5}, {1, 10} \[UndirectedEdge] {3, 2}, {1,
    11} \[UndirectedEdge] {3, 2}, {1, 12} \[UndirectedEdge] {3, 2}, {1, 13} \[UndirectedEdge] {3, 3}, {1,
    13} \[UndirectedEdge] {3, 5}, {1, 14} \[UndirectedEdge] {3, 3}, {1, 14} \[UndirectedEdge] {3, 4}, {1,
    15} \[UndirectedEdge] {3, 3}, {1, 15} \[UndirectedEdge] {3, 4}, {1, 16} \[UndirectedEdge] {3, 3}, {1,
    17} \[UndirectedEdge] {3, 4}, {1, 18} \[UndirectedEdge] {3, 5}, {1, 18} \[UndirectedEdge] {3, 6}, {1,
    19} \[UndirectedEdge] {3, 6}}}, {{{2, 1} \[UndirectedEdge] {0, 1}, {2, 1} \[UndirectedEdge] {0, 6}, {2,
    1} \[UndirectedEdge] {0, 7}, {2, 2} \[UndirectedEdge] {0, 1}, {2, 2} \[UndirectedEdge] {0, 5}, {2,
    2} \[UndirectedEdge] {0, 7}, {2, 3} \[UndirectedEdge] {0, 1}, {2, 3} \[UndirectedEdge] {0, 5}, {2,
    3} \[UndirectedEdge] {0, 6}, {2, 4} \[UndirectedEdge] {0, 5}, {2, 4} \[UndirectedEdge] {0, 6}, {2,
    4} \[UndirectedEdge] {0, 7}, {2, 5} \[UndirectedEdge] {0, 1}, {2, 5} \[UndirectedEdge] {0, 3}, {2,
    5} \[UndirectedEdge] {0, 8}, {2, 6} \[UndirectedEdge] {0, 3}, {2, 6} \[UndirectedEdge] {0, 4}, {2,
    6} \[UndirectedEdge] {0, 8}, {2, 7} \[UndirectedEdge] {0, 1}, {2, 7} \[UndirectedEdge] {0, 3}, {2,
    7} \[UndirectedEdge] {0, 4}, {2, 8} \[UndirectedEdge] {0, 1}, {2, 8} \[UndirectedEdge] {0, 4}, {2,
    8} \[UndirectedEdge] {0, 8}, {2, 9} \[UndirectedEdge] {0, 1}, {2, 9} \[UndirectedEdge] {0, 3}, {2,
    9} \[UndirectedEdge] {0, 7}, {2, 10} \[UndirectedEdge] {0, 1}, {2, 10} \[UndirectedEdge] {0, 2}, {2,
    10} \[UndirectedEdge] {0, 7}, {2, 11} \[UndirectedEdge] {0, 1}, {2, 11} \[UndirectedEdge] {0, 2}, {2,
    11} \[UndirectedEdge] {0, 3}, {2, 12} \[UndirectedEdge] {0, 2}, {2, 12} \[UndirectedEdge] {0, 3}, {2,
    12} \[UndirectedEdge] {0, 7}, {2, 13} \[UndirectedEdge] {0, 2}, {2, 13} \[UndirectedEdge] {0, 6}, {2,
    13} \[UndirectedEdge] {0, 7}, {2, 14} \[UndirectedEdge] {0, 1}, {2, 14} \[UndirectedEdge] {0, 2}, {2,
    14} \[UndirectedEdge] {0, 6}, {2, 15} \[UndirectedEdge] {0, 3}, {2, 15} \[UndirectedEdge] {0, 7}, {2,
    15} \[UndirectedEdge] {0, 8}, {2, 16} \[UndirectedEdge] {0, 1}, {2, 16} \[UndirectedEdge] {0, 7}, {2,
    16} \[UndirectedEdge] {0, 8}, {2, 17} \[UndirectedEdge] {0, 1}, {2, 17} \[UndirectedEdge] {0, 5}, {2,
    17} \[UndirectedEdge] {0, 8}, {2, 18} \[UndirectedEdge] {0, 5}, {2, 18} \[UndirectedEdge] {0, 7}, {2,
    18} \[UndirectedEdge] {0, 8}}, {{2, 1} \[UndirectedEdge] {1, 1}, {2, 1} \[UndirectedEdge] {1, 2}, {2,
    1} \[UndirectedEdge] {1, 3}, {2, 2} \[UndirectedEdge] {1, 1}, {2, 2} \[UndirectedEdge] {1, 4}, {2,
    2} \[UndirectedEdge] {1, 5}, {2, 3} \[UndirectedEdge] {1, 3}, {2, 3} \[UndirectedEdge] {1, 4}, {2,
    3} \[UndirectedEdge] {1, 6}, {2, 4} \[UndirectedEdge] {1, 2}, {2, 4} \[UndirectedEdge] {1, 5}, {2,
    4} \[UndirectedEdge] {1, 6}, {2, 5} \[UndirectedEdge] {1, 7}, {2, 5} \[UndirectedEdge] {1, 8}, {2,
    5} \[UndirectedEdge] {1, 9}, {2, 6} \[UndirectedEdge] {1, 7}, {2, 6} \[UndirectedEdge] {1, 10}, {2,
    6} \[UndirectedEdge] {1, 11}, {2, 7} \[UndirectedEdge] {1, 9}, {2, 7} \[UndirectedEdge] {1, 10}, {2,
    7} \[UndirectedEdge] {1, 12}, {2, 8} \[UndirectedEdge] {1, 8}, {2, 8} \[UndirectedEdge] {1, 11}, {2,
    8} \[UndirectedEdge] {1, 12}, {2, 9} \[UndirectedEdge] {1, 1}, {2, 9} \[UndirectedEdge] {1, 9}, {2,
    9} \[UndirectedEdge] {1, 13}, {2, 10} \[UndirectedEdge] {1, 1}, {2, 10} \[UndirectedEdge] {1, 14}, {2,
    10} \[UndirectedEdge] {1, 15}, {2, 11} \[UndirectedEdge] {1, 9}, {2, 11} \[UndirectedEdge] {1, 14}, {2,
    11} \[UndirectedEdge] {1, 16}, {2, 12} \[UndirectedEdge] {1, 13}, {2, 12} \[UndirectedEdge] {1, 15}, {2,
    12} \[UndirectedEdge] {1, 16}, {2, 13} \[UndirectedEdge] {1, 2}, {2, 13} \[UndirectedEdge] {1, 15}, {2,
    13} \[UndirectedEdge] {1, 17}, {2, 14} \[UndirectedEdge] {1, 3}, {2, 14} \[UndirectedEdge] {1, 14}, {2,
    14} \[UndirectedEdge] {1, 17}, {2, 15} \[UndirectedEdge] {1, 7}, {2, 15} \[UndirectedEdge] {1, 13}, {2,
    15} \[UndirectedEdge] {1, 18}, {2, 16} \[UndirectedEdge] {1, 1}, {2, 16} \[UndirectedEdge] {1, 8}, {2,
    16} \[UndirectedEdge] {1, 18}, {2, 17} \[UndirectedEdge] {1, 4}, {2, 17} \[UndirectedEdge] {1, 8}, {2,
    17} \[UndirectedEdge] {1, 19}, {2, 18} \[UndirectedEdge] {1, 5}, {2, 18} \[UndirectedEdge] {1, 18}, {2,
    18} \[UndirectedEdge] {1, 19}}, {{2, 1} \[UndirectedEdge] {2, 2}, {2, 1} \[UndirectedEdge] {2, 3}, {2,
    1} \[UndirectedEdge] {2, 4}, {2, 1} \[UndirectedEdge] {2, 9}, {2, 1} \[UndirectedEdge] {2, 10}, {2,
    1} \[UndirectedEdge] {2, 13}, {2, 1} \[UndirectedEdge] {2, 14}, {2, 1} \[UndirectedEdge] {2, 16}, {2,
    2} \[UndirectedEdge] {2, 3}, {2, 2} \[UndirectedEdge] {2, 4}, {2, 2} \[UndirectedEdge] {2, 9}, {2,
    2} \[UndirectedEdge] {2, 10}, {2, 2} \[UndirectedEdge] {2, 16}, {2, 2} \[UndirectedEdge] {2, 17}, {2,
    2} \[UndirectedEdge] {2, 18}, {2, 3} \[UndirectedEdge] {2, 4}, {2, 3} \[UndirectedEdge] {2, 14}, {2,
    3} \[UndirectedEdge] {2, 17}, {2, 4} \[UndirectedEdge] {2, 13}, {2, 4} \[UndirectedEdge] {2, 18}, {2,
    5} \[UndirectedEdge] {2, 6}, {2, 5} \[UndirectedEdge] {2, 7}, {2, 5} \[UndirectedEdge] {2, 8}, {2,
    5} \[UndirectedEdge] {2, 9}, {2, 5} \[UndirectedEdge] {2, 11}, {2, 5} \[UndirectedEdge] {2, 15}, {2,
    5} \[UndirectedEdge] {2, 16}, {2, 5} \[UndirectedEdge] {2, 17}, {2, 6} \[UndirectedEdge] {2, 7}, {2,
    6} \[UndirectedEdge] {2, 8}, {2, 6} \[UndirectedEdge] {2, 15}, {2, 7} \[UndirectedEdge] {2, 8}, {2,
    7} \[UndirectedEdge] {2, 9}, {2, 7} \[UndirectedEdge] {2, 11}, {2, 8} \[UndirectedEdge] {2, 16}, {2,
    8} \[UndirectedEdge] {2, 17}, {2, 9} \[UndirectedEdge] {2, 10}, {2, 9} \[UndirectedEdge] {2, 11}, {2,
    9} \[UndirectedEdge] {2, 12}, {2, 9} \[UndirectedEdge] {2, 15}, {2, 9} \[UndirectedEdge] {2, 16}, {2,
    10} \[UndirectedEdge] {2, 11}, {2, 10} \[UndirectedEdge] {2, 12}, {2, 10} \[UndirectedEdge] {2, 13}, {2,
    10} \[UndirectedEdge] {2, 14}, {2, 10} \[UndirectedEdge] {2, 16}, {2, 11} \[UndirectedEdge] {2, 12}, {2,
    11} \[UndirectedEdge] {2, 14}, {2, 12} \[UndirectedEdge] {2, 13}, {2, 12} \[UndirectedEdge] {2, 15}, {2,
    13} \[UndirectedEdge] {2, 14}, {2, 15} \[UndirectedEdge] {2, 16}, {2, 15} \[UndirectedEdge] {2, 18}, {2,
    16} \[UndirectedEdge] {2, 17}, {2, 16} \[UndirectedEdge] {2, 18}, {2, 17} \[UndirectedEdge] {2, 18}}, {{2,
    1} \[UndirectedEdge] {3, 1}, {2, 1} \[UndirectedEdge] {3, 4}, {2, 2} \[UndirectedEdge] {3, 1}, {2,
    2} \[UndirectedEdge] {3, 6}, {2, 3} \[UndirectedEdge] {3, 1}, {2, 4} \[UndirectedEdge] {3, 1}, {2,
    5} \[UndirectedEdge] {3, 2}, {2, 5} \[UndirectedEdge] {3, 5}, {2, 6} \[UndirectedEdge] {3, 2}, {2,
    7} \[UndirectedEdge] {3, 2}, {2, 8} \[UndirectedEdge] {3, 2}, {2, 9} \[UndirectedEdge] {3, 3}, {2,
    9} \[UndirectedEdge] {3, 5}, {2, 10} \[UndirectedEdge] {3, 3}, {2, 10} \[UndirectedEdge] {3, 4}, {2,
    11} \[UndirectedEdge] {3, 3}, {2, 12} \[UndirectedEdge] {3, 3}, {2, 13} \[UndirectedEdge] {3, 4}, {2,
    14} \[UndirectedEdge] {3, 4}, {2, 15} \[UndirectedEdge] {3, 5}, {2, 16} \[UndirectedEdge] {3, 5}, {2,
    16} \[UndirectedEdge] {3, 6}, {2, 17} \[UndirectedEdge] {3, 6}, {2, 18} \[UndirectedEdge] {3, 6}}}, {{{3,
    1} \[UndirectedEdge] {0, 1}, {3, 1} \[UndirectedEdge] {0, 5}, {3, 1} \[UndirectedEdge] {0, 6}, {3,
    1} \[UndirectedEdge] {0, 7}, {3, 2} \[UndirectedEdge] {0, 1}, {3, 2} \[UndirectedEdge] {0, 3}, {3,
    2} \[UndirectedEdge] {0, 4}, {3, 2} \[UndirectedEdge] {0, 8}, {3, 3} \[UndirectedEdge] {0, 1}, {3,
    3} \[UndirectedEdge] {0, 2}, {3, 3} \[UndirectedEdge] {0, 3}, {3, 3} \[UndirectedEdge] {0, 7}, {3,
    4} \[UndirectedEdge] {0, 1}, {3, 4} \[UndirectedEdge] {0, 2}, {3, 4} \[UndirectedEdge] {0, 6}, {3,
    4} \[UndirectedEdge] {0, 7}, {3, 5} \[UndirectedEdge] {0, 1}, {3, 5} \[UndirectedEdge] {0, 3}, {3,
    5} \[UndirectedEdge] {0, 7}, {3, 5} \[UndirectedEdge] {0, 8}, {3, 6} \[UndirectedEdge] {0, 1}, {3,
    6} \[UndirectedEdge] {0, 5}, {3, 6} \[UndirectedEdge] {0, 7}, {3, 6} \[UndirectedEdge] {0, 8}}, {{3,
    1} \[UndirectedEdge] {1, 1}, {3, 1} \[UndirectedEdge] {1, 2}, {3, 1} \[UndirectedEdge] {1, 3}, {3,
    1} \[UndirectedEdge] {1, 4}, {3, 1} \[UndirectedEdge] {1, 5}, {3, 1} \[UndirectedEdge] {1, 6}, {3,
    2} \[UndirectedEdge] {1, 7}, {3, 2} \[UndirectedEdge] {1, 8}, {3, 2} \[UndirectedEdge] {1, 9}, {3,
    2} \[UndirectedEdge] {1, 10}, {3, 2} \[UndirectedEdge] {1, 11}, {3, 2} \[UndirectedEdge] {1, 12}, {3,
    3} \[UndirectedEdge] {1, 1}, {3, 3} \[UndirectedEdge] {1, 9}, {3, 3} \[UndirectedEdge] {1, 13}, {3,
    3} \[UndirectedEdge] {1, 14}, {3, 3} \[UndirectedEdge] {1, 15}, {3, 3} \[UndirectedEdge] {1, 16}, {3,
    4} \[UndirectedEdge] {1, 1}, {3, 4} \[UndirectedEdge] {1, 2}, {3, 4} \[UndirectedEdge] {1, 3}, {3,
    4} \[UndirectedEdge] {1, 14}, {3, 4} \[UndirectedEdge] {1, 15}, {3, 4} \[UndirectedEdge] {1, 17}, {3,
    5} \[UndirectedEdge] {1, 1}, {3, 5} \[UndirectedEdge] {1, 7}, {3, 5} \[UndirectedEdge] {1, 8}, {3,
    5} \[UndirectedEdge] {1, 9}, {3, 5} \[UndirectedEdge] {1, 13}, {3, 5} \[UndirectedEdge] {1, 18}, {3,
    6} \[UndirectedEdge] {1, 1}, {3, 6} \[UndirectedEdge] {1, 4}, {3, 6} \[UndirectedEdge] {1, 5}, {3,
    6} \[UndirectedEdge] {1, 8}, {3, 6} \[UndirectedEdge] {1, 18}, {3, 6} \[UndirectedEdge] {1, 19}}, {{3,
    1} \[UndirectedEdge] {2, 1}, {3, 1} \[UndirectedEdge] {2, 2}, {3, 1} \[UndirectedEdge] {2, 3}, {3,
    1} \[UndirectedEdge] {2, 4}, {3, 2} \[UndirectedEdge] {2, 5}, {3, 2} \[UndirectedEdge] {2, 6}, {3,
    2} \[UndirectedEdge] {2, 7}, {3, 2} \[UndirectedEdge] {2, 8}, {3, 3} \[UndirectedEdge] {2, 9}, {3,
    3} \[UndirectedEdge] {2, 10}, {3, 3} \[UndirectedEdge] {2, 11}, {3, 3} \[UndirectedEdge] {2, 12}, {3,
    4} \[UndirectedEdge] {2, 1}, {3, 4} \[UndirectedEdge] {2, 10}, {3, 4} \[UndirectedEdge] {2, 13}, {3,
    4} \[UndirectedEdge] {2, 14}, {3, 5} \[UndirectedEdge] {2, 5}, {3, 5} \[UndirectedEdge] {2, 9}, {3,
    5} \[UndirectedEdge] {2, 15}, {3, 5} \[UndirectedEdge] {2, 16}, {3, 6} \[UndirectedEdge] {2, 2}, {3,
    6} \[UndirectedEdge] {2, 16}, {3, 6} \[UndirectedEdge] {2, 17}, {3, 6} \[UndirectedEdge] {2, 18}}, {{3,
    1} \[UndirectedEdge] {3, 4}, {3, 1} \[UndirectedEdge] {3, 6}, {3, 2} \[UndirectedEdge] {3, 5}, {3,
    3} \[UndirectedEdge] {3, 4}, {3, 3} \[UndirectedEdge] {3, 5}, {3, 5} \[UndirectedEdge] {3, 6}}}}
]


(*******************************************************************************)
MTSection["Graph colouring"]

(* Vertex colouring *)

With[{g = RandomGraph[{1000,3000}]},
  MT[
    With[{col = IGVertexColoring[g]},
      And @@ Unequal @@@ Map[col[[#]]&, EdgeList[g], {2}]
    ],
    True
  ]
]

MT[
  IGVertexColoring@IGEmptyGraph[1],
  {1}
]

MT[
  IGVertexColoring@IGEmptyGraph[],
  {}
]

MT[
  IGEdgeColoring@IGEmptyGraph[#],
  {}
]& /@ {0, 1, 2, 3}


MT[
  AllTrue[
    GraphData[10],
    Max@IGMinimumVertexColoring@GraphData[#] == GraphData[#, "ChromaticNumber"]&
  ],
  True
]


MT[
  IGChromaticNumber /@ NestList[IGMycielskian, IGEmptyGraph[], 6],
  Range[0, 6]
]


(*******************************************************************************)
MTSection["Proximity graphs"]

MT[
  #[{}],
  Graph[{},{}],
  SameTest -> IGSameGraphQ
]& /@ {IGDelaunayGraph, IGRelativeNeighborhoodGraph, IGGabrielGraph, IGLuneBetaSkeleton[#, 1.5]&, IGLuneBetaSkeleton[#, 0.5]&, IGCircleBetaSkeleton[#, 1.5]&, IGCircleBetaSkeleton[#, 0.5]&}

MT[
  #[{{0,0}}],
  Graph[{1},{}],
  SameTest -> IGSameGraphQ
]& /@ {IGDelaunayGraph, IGRelativeNeighborhoodGraph, IGGabrielGraph, IGLuneBetaSkeleton[#, 1.5]&, IGLuneBetaSkeleton[#, 0.5]&, IGCircleBetaSkeleton[#, 1.5]&, IGCircleBetaSkeleton[#, 0.5]&}

MT[
  #[{{0,0}, {1,2}}],
  CompleteGraph[2],
  SameTest -> IGSameGraphQ
]& /@ {IGDelaunayGraph, IGRelativeNeighborhoodGraph, IGGabrielGraph, IGLuneBetaSkeleton[#, 1.5]&, IGLuneBetaSkeleton[#, 0.5]&, IGCircleBetaSkeleton[#, 1.5]&, IGCircleBetaSkeleton[#, 0.5]&}

MT[
  IGDelaunayGraph[{{0,0}, {1,1}, {2,2}}],
  Graph[{1<->2, 2<->3}],
  SameTest -> IGSameGraphQ
]

MT[
  IGDelaunayGraph[{{0,0,0}, {1,1,1}, {2,2,2}}],
  Graph[{1<->2, 2<->3}],
  SameTest -> IGSameGraphQ
]

MT[
  IGDelaunayGraph[{{0,0}, {1,1}, {1,2}}],
  Graph[{1<->2, 2<->3, 1<->3}],
  SameTest -> IGSameGraphQ
]

MT[
  IGDelaunayGraph[{{0,0,0}, {1,1,1}, {1,2,2}}],
  Graph[{1<->2, 2<->3, 1<->3}],
  SameTest -> IGSameGraphQ
]

MT[
  IGDelaunayGraph[{{0, 0}, {1, 0}, {0, 1}, {0, 0}}],
  $Failed,
  {IGDelaunayGraph::dupl}
]

points2d = RandomReal[1, {100, 2}];
MT[
  IGDelaunayGraph[points2d],
  IGMeshGraph@DelaunayMesh[points2d],
  SameTest -> IGSameGraphQ
]

points3d = RandomReal[1, {100, 2}];
MT[
  IGDelaunayGraph[points3d],
  IGMeshGraph@DelaunayMesh[points3d],
  SameTest -> IGSameGraphQ
]

bsPoints = {{-0.59, -0.09}, {0.07, 0.85}, {-0.48, -0.22}, {0.56, -0.12}, {0.02, -0.35}, {-0.05, 1.}, {-0.76, 0.32}, {-0.37, 0.}, {-0.03, 0.39}, {-0.65, -0.49}, {-0.06, 0.6}, {-0.53, 0.76}, {0.91, -0.25}, {0.22, -0.72}, {0.27, -0.18}, {-0.75, -0.3}, {0.85, 0.09}, {0.63, 0.2}, {0.57, 0.16}, {-0.43, -0.51}};

MT[
  IGLuneBetaSkeleton[bsPoints, 0.8],
  Graph@{1 <-> 3, 1 <-> 7, 1 <-> 8, 1 <-> 10, 1 <-> 16, 2 <-> 6, 2 <-> 11,
    2 <-> 12, 2 <-> 18, 2 <-> 19, 3 <-> 5, 3 <-> 8, 3 <-> 10, 3 <-> 16,
    3 <-> 20, 4 <-> 9, 4 <-> 13, 4 <-> 14, 4 <-> 15, 4 <-> 17, 4 <-> 18,
    4 <-> 19, 5 <-> 8, 5 <-> 9, 5 <-> 14, 5 <-> 15, 5 <-> 20, 6 <-> 11,
    6 <-> 12, 7 <-> 8, 7 <-> 9, 7 <-> 11, 7 <-> 12, 7 <-> 16, 8 <-> 9,
    8 <-> 12, 8 <-> 15, 8 <-> 19, 9 <-> 11, 9 <-> 12, 9 <-> 15, 9 <-> 18,
    9 <-> 19, 10 <-> 16, 10 <-> 20, 11 <-> 12, 11 <-> 18, 11 <-> 19,
    13 <-> 14, 13 <-> 17, 13 <-> 18, 13 <-> 19, 14 <-> 15, 14 <-> 20,
    15 <-> 19, 16 <-> 20, 17 <-> 18, 17 <-> 19, 18 <-> 19},
  SameTest -> IGSameGraphQ
]

MT[
  IGLuneBetaSkeleton[bsPoints, 1.5],
  Graph@{10 <-> 16, 3 <-> 16, 3 <-> 8, 1 <-> 3, 1 <-> 16, 11 <-> 12, 8 <-> 9,
    1 <-> 7, 5 <-> 14, 4 <-> 13, 4 <-> 15, 9 <-> 15, 17 <-> 18,
    18 <-> 19, 3 <-> 20, 5 <-> 8, 1 <-> 8, 10 <-> 20, 9 <-> 11, 2 <-> 6,
    7 <-> 12, 5 <-> 20, 5 <-> 15, 9 <-> 19, 6 <-> 12, 2 <-> 11,
    13 <-> 17, 4 <-> 19},
  SameTest -> IGSameGraphQ
]

MT[
  IGLuneBetaSkeleton[bsPoints, 2.5],
  Graph[
    Range@Length[bsPoints],
    {10 <-> 16, 1 <-> 3, 1 <-> 16, 11 <-> 12, 8 <-> 9, 1 <-> 7, 5 <-> 14,
    4 <-> 15, 17 <-> 18, 18 <-> 19, 3 <-> 20, 10 <-> 20, 9 <-> 11,
    2 <-> 6, 7 <-> 12, 5 <-> 15, 2 <-> 11, 4 <-> 19}
  ],
  SameTest -> IGSameGraphQ
]

MT[
  IGGabrielGraph[bsPoints],
  Graph@{10 <-> 16, 3 <-> 16, 3 <-> 10, 3 <-> 8, 1 <-> 3, 1 <-> 16, 11 <-> 12,
    8 <-> 9, 1 <-> 7, 3 <-> 5, 5 <-> 14, 4 <-> 13, 4 <-> 15, 9 <-> 15,
    17 <-> 18, 18 <-> 19, 3 <-> 20, 5 <-> 8, 1 <-> 8, 10 <-> 20,
    9 <-> 11, 2 <-> 6, 7 <-> 12, 5 <-> 20, 5 <-> 15, 4 <-> 17, 9 <-> 19,
    6 <-> 12, 2 <-> 11, 13 <-> 17, 4 <-> 19},
  SameTest -> IGSameGraphQ
]

MT[
  IGRelativeNeighborhoodGraph[bsPoints],
  Graph@{10 <-> 16, 1 <-> 3, 1 <-> 16, 11 <-> 12, 8 <-> 9, 1 <-> 7, 5 <-> 14,
    4 <-> 15, 17 <-> 18, 18 <-> 19, 3 <-> 20, 1 <-> 8, 10 <-> 20,
    9 <-> 11, 2 <-> 6, 7 <-> 12, 5 <-> 20, 5 <-> 15, 9 <-> 19, 2 <-> 11,
    13 <-> 17, 4 <-> 19},
  SameTest -> IGSameGraphQ
]

MT[
  IGCircleBetaSkeleton[bsPoints, 0.9],
  Graph@{1 <-> 3, 1 <-> 7, 1 <-> 8, 1 <-> 10, 1 <-> 16, 2 <-> 6, 2 <-> 11,
    2 <-> 12, 2 <-> 18, 2 <-> 19, 3 <-> 5, 3 <-> 8, 3 <-> 10, 3 <-> 16,
    3 <-> 20, 4 <-> 9, 4 <-> 13, 4 <-> 14, 4 <-> 15, 4 <-> 17, 4 <-> 19,
    5 <-> 8, 5 <-> 9, 5 <-> 14, 5 <-> 15, 5 <-> 20, 6 <-> 11, 6 <-> 12,
    7 <-> 8, 7 <-> 9, 7 <-> 11, 7 <-> 12, 8 <-> 9, 8 <-> 12, 8 <-> 15,
    8 <-> 19, 9 <-> 11, 9 <-> 15, 9 <-> 19, 10 <-> 16, 10 <-> 20,
    11 <-> 12, 11 <-> 18, 13 <-> 14, 13 <-> 17, 13 <-> 19, 14 <-> 15,
    14 <-> 20, 15 <-> 19, 17 <-> 18, 18 <-> 19},
  SameTest -> IGSameGraphQ
]

MT[
  IGCircleBetaSkeleton[bsPoints, 1.1],
  Graph@{10 <-> 16, 1 <-> 3, 11 <-> 12, 8 <-> 9, 1 <-> 7, 5 <-> 14, 4 <-> 13,
    4 <-> 15, 17 <-> 18, 18 <-> 19, 3 <-> 20, 10 <-> 20, 9 <-> 11,
    2 <-> 6, 7 <-> 12, 5 <-> 15, 2 <-> 11, 13 <-> 17, 4 <-> 19},
  SameTest -> IGSameGraphQ
]

MT[
  IGRelativeNeighborhoodGraph@GraphEmbedding@IGTriangularLattice[5],
  IGEmptyGraph@VertexCount@IGTriangularLattice[5],
  SameTest -> IGSameGraphQ
]

MT[
  IGLuneBetaSkeleton[GraphEmbedding@IGTriangularLattice[4], 1.99],
  IGTriangularLattice[4],
  SameTest -> IGSameGraphQ
]

circlePoints[n_] := Table[{Cos[x], Sin[x]}, {x, 0, 2Pi - 2Pi/n, 2Pi/n}]

MT[
  IGGabrielGraph[circlePoints[7]],
  CycleGraph[7],
  SameTest -> IGSameGraphQ
]


(*******************************************************************************)
MTSection["Test remaining functions"]

(* IGData *)

MT[
  IGData /@ IGData[];,
  Null
]

MT[
  And @@ MapThread[IGIsomorphicQ, {IGData[{"AllDirectedGraphs", 3}], Values@IGData["MANTriadLabels"]}],
  True
]

(* IGArticulationPoints *)

connCompCount[g_] := Length@ConnectedComponents[g]

MT[
  AllTrue[
    connCompCount@VertexDelete[ugs, #] & /@ IGArticulationPoints[ugs],
    # > connCompCount[ugs] &
  ],
  True
]

MT[
  AllTrue[
    connCompCount@VertexDelete[dgs, #] & /@ IGArticulationPoints[dgs],
    # > connCompCount[dgs] &
  ],
  True
]

MT[
  AllTrue[
    connCompCount@VertexDelete[ugs, #] & /@ IGArticulationPoints[ugs],
    # > connCompCount[ugi] &
  ],
  True
]

MT[
  AllTrue[
    connCompCount@VertexDelete[dgs, #] & /@ IGArticulationPoints[dgs],
    # > connCompCount[dgi] &
  ],
  True
]


(* IGMinimumSeparators *)

MT[
  IGMinimumSeparators[#] =!= {} & /@ ulist,
  ConnectedGraphQ /@ ulist
]

(* IGAdjacentTriangleCount *)

MT[
  IGAdjacentTriangleCount[#],
  With[{am = AdjacencyMatrix[#]}, Normal@Diagonal[am.am.am]/2]
]& /@ ulist

MT[
  IGAdjacentTriangleCount[empty],
  {}
]

MT[
  IGAdjacentTriangleCount[Graph[{1,2,3},{}]],
  {0,0,0}
]

(* IGGraphicalQ *)

MT[
  IGGraphicalQ@VertexDegree[#],
  True
]& /@ ulist

MT[
  IGGraphicalQ[VertexInDegree[#], VertexOutDegree[#]],
  True
]& /@ dlist

MT[
  IGGraphicalQ[VertexOutDegree[#], VertexInDegree[#]],
  True
]& /@ dlist

MT[
  IGGraphicalQ[{1,2,1}],
  True
]

MT[
  IGGraphicalQ[{1,2,2}],
  False
]

MT[
  IGGraphicalQ[{1,0},{1,0}],
  False
]

MT[
  IGGraphicalQ[{1,0},{0,1}],
  True
]

MT[
  IGGraphicalQ[{}],
  True
]


(* Index edge list *)

MT[
  IGIndexEdgeList[IGEmptyGraph[]],
  {}
]

MT[
  IGIndexEdgeList[IGEmptyGraph[5]],
  {}
]

MT[
  IGIndexEdgeList[Graph[Range[10], {5<->6}]],
  {{5,6}}
]


MT[
  IGIndexEdgeList[#],
  List @@@ EdgeList@IndexGraph[#]
]& /@ {dgs, dgi, dmulti, bidi, dbipartite}

MT[
  Sort /@ IGIndexEdgeList[#],
  Sort /@ List @@@ EdgeList@IndexGraph[#]
]& /@ {ugs, ugi, umulti, bipartite, football, dolphin}


(*******************************************************************************)
MTSection["Utilities - General"]

(* IGUndirectedGraph *)

Outer[
  MT[
    sameGraphQ[
      IGUndirectedGraph[#1, #2],
      #1
    ],
    True
  ]&,
  {empty, edgeless, ugi, ugs, umulti},
  {"Simple", "All", "Reciprocal"}
]

g = Graph[{1,2,3,4,5}, {1->2, 1->2, 3->3, 3->4, 4->3}];

MT[
  IGUndirectedGraph[g],
  IGUndirectedGraph[g, "Simple"]
]

MT[
  sameGraphQ[
    IGUndirectedGraph[g, "Simple"],
    Graph[{1,2,3,4,5}, {1<->2, 3<->3, 3<->4}]
  ],
  True
]

MT[
  sameGraphQ[
    IGUndirectedGraph[g, "All"],
    Graph[{1, 2, 3, 4, 5}, {1 <-> 2, 1 <-> 2, 3 <-> 3, 3 <-> 4, 4 <-> 3}]
  ],
  True
]

MT[
  sameGraphQ[
    IGUndirectedGraph[g, "Reciprocal"],
    Graph[{1, 2, 3, 4, 5}, {3 <-> 3, 3 <-> 4}]
  ],
  True
]

(* IGNullGraphQ *)

MT[
  IGNullGraphQ[1], (* False for non-graph *)
  False
]

MT[
  IGNullGraphQ[empty],
  True
]

MT[
  IGNullGraphQ[edgeless],
  False
]

MT[
  IGNullGraphQ[ugi],
  False
]

MT[
  IGNullGraphQ[ugs],
  False
]

(* IGSimpleGraph *)

MT[
  sameGraphQ[
    IGSimpleGraph[#],
    #
  ],
  True
]& /@ {empty, edgeless, ugi, ugs, dgi, dgs, bidi}

ug = Graph[{7, 9, 10, 5, 6, 0, 8, 2, 4, 3, 1}, {7 <-> 9, 10 <-> 5,
  6 <-> 5, 10 <-> 0, 5 <-> 8, 9 <-> 9, 2 <-> 4, 9 <-> 3, 3 <-> 1,
  6 <-> 0, 7 <-> 3, 0 <-> 5, 7 <-> 3, 0 <-> 0, 8 <-> 6, 8 <-> 0,
  10 <-> 7, 9 <-> 0, 5 <-> 3, 3 <-> 0, 6 <-> 1, 3 <-> 10, 8 <-> 8,
  7 <-> 8, 6 <-> 8, 4 <-> 8, 6 <-> 8, 4 <-> 0, 6 <-> 3, 3 <-> 5,
  5 <-> 0, 6 <-> 0, 2 <-> 1, 5 <-> 0, 10 <-> 0, 0 <-> 0, 5 <-> 6,
  1 <-> 1, 1 <-> 2, 3 <-> 10, 0 <-> 0, 10 <-> 6, 2 <-> 1, 10 <-> 10,
  7 <-> 10, 2 <-> 10, 9 <-> 5, 1 <-> 10, 7 <-> 3, 0 <-> 10}];

MT[
  sameGraphQ[
    IGSimpleGraph[ug],
    Graph[{7, 9, 10, 5, 6, 0, 8, 2, 4, 3, 1}, {6 <-> 1, 2 <-> 1, 10 <-> 6,
      2 <-> 10, 1 <-> 10, 10 <-> 5, 6 <-> 5, 5 <-> 8, 8 <-> 6, 3 <-> 1,
      5 <-> 3, 3 <-> 10, 6 <-> 3, 10 <-> 0, 6 <-> 0, 0 <-> 5, 8 <-> 0,
      3 <-> 0, 2 <-> 4, 4 <-> 8, 4 <-> 0, 9 <-> 3, 9 <-> 0, 9 <-> 5,
      7 <-> 9, 7 <-> 3, 10 <-> 7, 7 <-> 8}]
  ],
  True
]

MT[
  sameGraphQ[
    IGSimpleGraph[ug, SelfLoops -> True],
    Graph[{7, 9, 10, 5, 6, 0, 8, 2, 4, 3, 1}, {7 <-> 9, 7 <-> 10, 7 <-> 8,
      7 <-> 3, 9 <-> 9, 9 <-> 5, 9 <-> 0, 9 <-> 3, 10 <-> 10, 10 <-> 5,
      10 <-> 6, 10 <-> 0, 10 <-> 2, 10 <-> 3, 10 <-> 1, 5 <-> 6, 5 <-> 0,
      5 <-> 8, 5 <-> 3, 6 <-> 0, 6 <-> 8, 6 <-> 3, 6 <-> 1, 0 <-> 0,
      0 <-> 8, 0 <-> 4, 0 <-> 3, 8 <-> 8, 8 <-> 4, 2 <-> 4, 2 <-> 1,
      3 <-> 1, 1 <-> 1}]
  ],
  True
]

MT[
  sameGraphQ[
    IGSimpleGraph[ug, "MultipleEdges" -> True],
    Graph[{7, 9, 10, 5, 6, 0, 8, 2, 4, 3, 1}, {7 <-> 9, 7 <-> 10, 7 <-> 10,
      7 <-> 8, 7 <-> 3, 7 <-> 3, 7 <-> 3, 9 <-> 5, 9 <-> 0, 9 <-> 3,
      10 <-> 5, 10 <-> 6, 10 <-> 0, 10 <-> 0, 10 <-> 0, 10 <-> 2, 10 <-> 3,
      10 <-> 3, 10 <-> 1, 5 <-> 6, 5 <-> 6, 5 <-> 0, 5 <-> 0, 5 <-> 0,
      5 <-> 8, 5 <-> 3, 5 <-> 3, 6 <-> 0, 6 <-> 0, 6 <-> 8, 6 <-> 8,
      6 <-> 8, 6 <-> 3, 6 <-> 1, 0 <-> 8, 0 <-> 4, 0 <-> 3, 8 <-> 4,
      2 <-> 4, 2 <-> 1, 2 <-> 1, 2 <-> 1, 3 <-> 1}]
  ],
  True
]

dg = Graph[{7, 9, 10, 5, 6, 0, 8, 2, 4, 3, 1}, {7 \[DirectedEdge] 9,
  10 \[DirectedEdge] 5, 6 \[DirectedEdge] 5, 10 \[DirectedEdge] 0,
  5 \[DirectedEdge] 8, 9 \[DirectedEdge] 9, 2 \[DirectedEdge] 4,
  9 \[DirectedEdge] 3, 3 \[DirectedEdge] 1, 6 \[DirectedEdge] 0,
  7 \[DirectedEdge] 3, 0 \[DirectedEdge] 5, 7 \[DirectedEdge] 3,
  0 \[DirectedEdge] 0, 8 \[DirectedEdge] 6, 8 \[DirectedEdge] 0,
  10 \[DirectedEdge] 7, 9 \[DirectedEdge] 0, 5 \[DirectedEdge] 3,
  3 \[DirectedEdge] 0, 6 \[DirectedEdge] 1, 3 \[DirectedEdge] 10,
  8 \[DirectedEdge] 8, 7 \[DirectedEdge] 8, 6 \[DirectedEdge] 8,
  4 \[DirectedEdge] 8, 6 \[DirectedEdge] 8, 4 \[DirectedEdge] 0,
  6 \[DirectedEdge] 3, 3 \[DirectedEdge] 5, 5 \[DirectedEdge] 0,
  6 \[DirectedEdge] 0, 2 \[DirectedEdge] 1, 5 \[DirectedEdge] 0,
  10 \[DirectedEdge] 0, 0 \[DirectedEdge] 0, 5 \[DirectedEdge] 6,
  1 \[DirectedEdge] 1, 1 \[DirectedEdge] 2, 3 \[DirectedEdge] 10,
  0 \[DirectedEdge] 0, 10 \[DirectedEdge] 6, 2 \[DirectedEdge] 1,
  10 \[DirectedEdge] 10, 7 \[DirectedEdge] 10, 2 \[DirectedEdge] 10,
  9 \[DirectedEdge] 5, 1 \[DirectedEdge] 10, 7 \[DirectedEdge] 3,
  0 \[DirectedEdge] 10}];

MT[
  sameGraphQ[
    IGSimpleGraph[dg],
    Graph[{7, 9, 10, 5, 6, 0, 8, 2, 4, 3, 1}, {6 \[DirectedEdge] 1,
      2 \[DirectedEdge] 1, 1 \[DirectedEdge] 2, 10 \[DirectedEdge] 6,
      2 \[DirectedEdge] 10, 1 \[DirectedEdge] 10, 5 \[DirectedEdge] 6,
      10 \[DirectedEdge] 5, 6 \[DirectedEdge] 5, 8 \[DirectedEdge] 6,
      5 \[DirectedEdge] 8, 6 \[DirectedEdge] 8, 3 \[DirectedEdge] 1,
      3 \[DirectedEdge] 10, 3 \[DirectedEdge] 5, 5 \[DirectedEdge] 3,
      6 \[DirectedEdge] 3, 0 \[DirectedEdge] 5, 0 \[DirectedEdge] 10,
      10 \[DirectedEdge] 0, 6 \[DirectedEdge] 0, 8 \[DirectedEdge] 0,
      3 \[DirectedEdge] 0, 5 \[DirectedEdge] 0, 4 \[DirectedEdge] 8,
      4 \[DirectedEdge] 0, 2 \[DirectedEdge] 4, 9 \[DirectedEdge] 3,
      9 \[DirectedEdge] 0, 9 \[DirectedEdge] 5, 7 \[DirectedEdge] 9,
      7 \[DirectedEdge] 3, 7 \[DirectedEdge] 8, 7 \[DirectedEdge] 10,
      10 \[DirectedEdge] 7}]
  ],
  True
]

MT[
  sameGraphQ[
    IGSimpleGraph[dg, SelfLoops -> True],
    Graph[{7, 9, 10, 5, 6, 0, 8, 2, 4, 3, 1}, {7 \[DirectedEdge] 9,
    7 \[DirectedEdge] 10, 7 \[DirectedEdge] 8, 7 \[DirectedEdge] 3,
    9 \[DirectedEdge] 9, 9 \[DirectedEdge] 5, 9 \[DirectedEdge] 0,
    9 \[DirectedEdge] 3, 10 \[DirectedEdge] 7, 10 \[DirectedEdge] 10,
    10 \[DirectedEdge] 5, 10 \[DirectedEdge] 6, 10 \[DirectedEdge] 0,
    5 \[DirectedEdge] 6, 5 \[DirectedEdge] 0, 5 \[DirectedEdge] 8,
    5 \[DirectedEdge] 3, 6 \[DirectedEdge] 5, 6 \[DirectedEdge] 0,
    6 \[DirectedEdge] 8, 6 \[DirectedEdge] 3, 6 \[DirectedEdge] 1,
    0 \[DirectedEdge] 10, 0 \[DirectedEdge] 5, 0 \[DirectedEdge] 0,
    8 \[DirectedEdge] 6, 8 \[DirectedEdge] 0, 8 \[DirectedEdge] 8,
    2 \[DirectedEdge] 10, 2 \[DirectedEdge] 4, 2 \[DirectedEdge] 1,
    4 \[DirectedEdge] 0, 4 \[DirectedEdge] 8, 3 \[DirectedEdge] 10,
    3 \[DirectedEdge] 5, 3 \[DirectedEdge] 0, 3 \[DirectedEdge] 1,
    1 \[DirectedEdge] 10, 1 \[DirectedEdge] 2, 1 \[DirectedEdge] 1}]
  ],
  True
]

MT[
  sameGraphQ[
    IGSimpleGraph[dg, "MultipleEdges" -> True],
    Graph[{7, 9, 10, 5, 6, 0, 8, 2, 4, 3, 1}, {7 \[DirectedEdge] 9,
    7 \[DirectedEdge] 10, 7 \[DirectedEdge] 8, 7 \[DirectedEdge] 3,
    7 \[DirectedEdge] 3, 7 \[DirectedEdge] 3, 9 \[DirectedEdge] 5,
    9 \[DirectedEdge] 0, 9 \[DirectedEdge] 3, 10 \[DirectedEdge] 7,
    10 \[DirectedEdge] 5, 10 \[DirectedEdge] 6, 10 \[DirectedEdge] 0,
    10 \[DirectedEdge] 0, 5 \[DirectedEdge] 6, 5 \[DirectedEdge] 0,
    5 \[DirectedEdge] 0, 5 \[DirectedEdge] 8, 5 \[DirectedEdge] 3,
    6 \[DirectedEdge] 5, 6 \[DirectedEdge] 0, 6 \[DirectedEdge] 0,
    6 \[DirectedEdge] 8, 6 \[DirectedEdge] 8, 6 \[DirectedEdge] 3,
    6 \[DirectedEdge] 1, 0 \[DirectedEdge] 10, 0 \[DirectedEdge] 5,
    8 \[DirectedEdge] 6, 8 \[DirectedEdge] 0, 2 \[DirectedEdge] 10,
    2 \[DirectedEdge] 4, 2 \[DirectedEdge] 1, 2 \[DirectedEdge] 1,
    4 \[DirectedEdge] 0, 4 \[DirectedEdge] 8, 3 \[DirectedEdge] 10,
    3 \[DirectedEdge] 10, 3 \[DirectedEdge] 5, 3 \[DirectedEdge] 0,
    3 \[DirectedEdge] 1, 1 \[DirectedEdge] 10, 1 \[DirectedEdge] 2}]
  ],
  True
]

(* IGEdgeWeightedQ, IGVertexWeightedQ *)

vwg = Graph[{1,2,3},{1<->2},VertexWeight->{1,2,3}];

MT[
  IGEdgeWeightedQ[#],
  False
]& /@ Hold[
  empty, edgeless, ugi, ugs, dgi, dgs, umulti, dmulti,
  vwg
] // ReleaseHold

MT[
  IGVertexWeightedQ[#],
  False
]& /@ Hold[
  empty, edgeless, ugi, ugs, dgi, dgs, umulti, dmulti,
  wugi, wugs, wdgi, wdgs
] // ReleaseHold

MT[
  IGEdgeWeightedQ[#],
  True
]& /@ Hold[wugi, wugs, wdgi, wdgs] // ReleaseHold

MT[
  IGVertexWeightedQ[vwg],
  True
]

MT[
  IGEdgeWeightedQ@IGEmptyGraph[],
  False
]

MT[
  IGVertexWeightedQ@IGEmptyGraph[],
  False
]

MT[
  IGVertexWeightedQ@IGEmptyGraph[],
  False
]

MT[
  IGEdgeWeightedQ@IGEmptyGraph[3, VertexWeight -> {1, 2, 3}],
  False
]

(* IGWeighedAdjacencyGraph *)

MT[
  sameGraphQ[
    IGWeightedAdjacencyGraph[{{0,1},{2,3}}],
    Graph[{1->2, 2->1, 2->2}]
  ],
  True
]

MT[
  IGEdgeProp[EdgeWeight]@IGWeightedAdjacencyGraph[{{0,1},{2,3}}],
  {1,2,3}
]

MT[
  sameGraphQ[
    IGWeightedAdjacencyGraph[{{0,1},{1,3}}],
    Graph[{1<->2, 2<->2}]
  ],
  True
]

MT[
  IGEdgeProp[EdgeWeight]@IGWeightedAdjacencyGraph[{{0,1},{1,3}}],
  {1,3}
]

MT[
  sameGraphQ[
    IGWeightedAdjacencyGraph[{"a","b"}, {{0,1},{2,3}}],
    Graph[{"a"->"b", "b"->"a", "b"->"b"}]
  ],
  True
]

(* Vertex strength *)

MT[
  IGVertexStrength[#],
  VertexDegree[#]
]& /@ {ugi, ugs, dgi, dgs, umulti, dmulti}

MT[
  IGVertexInStrength[#],
  VertexInDegree[#]
]& /@ {ugi, ugs, dgi, dgs, umulti, dmulti}

MT[
  IGVertexOutStrength[#],
  VertexOutDegree[#]
]& /@ {ugi, ugs, dgi, dgs, umulti, dmulti}

MT[
  IGVertexStrength[#],
  Function[v, IGVertexStrength[#, v]] /@ VertexList[#]
]& /@ {ugi, ugs, dgi, dgs, umulti, dmulti}

MT[
  IGVertexInStrength[#],
  Function[v, IGVertexInStrength[#, v]] /@ VertexList[#]
]& /@ {ugi, ugs, dgi, dgs, umulti, dmulti}

MT[
  IGVertexOutStrength[#],
  Function[v, IGVertexOutStrength[#, v]] /@ VertexList[#]
]& /@ {ugi, ugs, dgi, dgs, umulti, dmulti}

MT[
  IGVertexStrength[wugi],
  Function[v,
    Total[
      PropertyValue[{wugi, #}, EdgeWeight]& /@ Cases[EdgeList[wugi], (v ~ UndirectedEdge ~ _) | (_ ~ UndirectedEdge ~ v)]
    ]
  ] /@ VertexList[wugi],
  SameTest -> Equal
]

MT[
  IGVertexOutStrength[wdgi],
  Function[v,
    Total[
      PropertyValue[{wdgi, #}, EdgeWeight]& /@ Cases[EdgeList[wdgi], (v ~ DirectedEdge ~ _)]
    ]
  ] /@ VertexList[wdgi],
  SameTest -> Equal
]

MT[
  IGVertexInStrength[wdgi],
  Function[v,
    Total[
      PropertyValue[{wdgi, #}, EdgeWeight]& /@ Cases[EdgeList[wdgi], (_ ~ DirectedEdge ~ v)]
    ]
  ] /@ VertexList[wdgi],
  SameTest -> Equal
]

MT[
  IGVertexStrength[Graph[{1<->2, 1<->2, 2<->3}, EdgeWeight -> {4,5,6}]],
  {9, 15, 6}
]


(*******************************************************************************)
MTSection["Utilities - Property operations"]

(* Property map operators *)

MT[
  IGEdgeProp[EdgeWeight][Graph[{1<->2}, EdgeWeight -> {12}]],
  {12}
]

MT[
  IGEdgeProp[EdgeCapacity][Graph[{1<->2}, EdgeCapacity -> {12}]],
  {12}
]

(* Multigraphs *)
MT[
  IGEdgeProp[EdgeWeight][Graph[{1<->2, 1<->2}, EdgeWeight -> {12,5}]],
  {12,5}
]

MT[
  IGEdgeProp[EdgeCapacity][Graph[{1<->2, 1<->2}, EdgeCapacity -> {12,5}]],
  {12,5}
]

(* Custom property *)
MT[
  IGEdgeProp["Foo"][ Graph[{Property[1<->2, "Foo" -> 37]}] ],
  {37}
]

MT[
  IGEdgeProp["Foo"][ Graph[{1<->2}] ],
  {Missing["Nonexistent"]}
]

{
  MT[
    IGEdgeProp[#][IGEmptyGraph[0]],
    {}
  ],
  MT[
    IGEdgeProp[#][IGEmptyGraph[1]],
    {}
  ]
}& /@ {EdgeWeight, EdgeCapacity, "foo"}

(* Vertex properties *)

MT[
  IGVertexProp[#][IGEmptyGraph[0]],
  {}
]& /@ {EdgeWeight, EdgeCapacity, "foo"}

MT[
  IGVertexProp[VertexWeight][Graph[{1<->2}, VertexWeight -> {3,4}]],
  {3,4}
]

MT[
  IGVertexProp[VertexCapacity][Graph[{1<->2}, VertexCapacity -> {3,4}]],
  {3,4}
]

(* Custom property *)
MT[
  IGVertexProp["Foo"][ Graph[{Property[1, "Foo" -> 37]}, {1<->2}] ],
  {37, Missing["Nonexistent"]}
]
