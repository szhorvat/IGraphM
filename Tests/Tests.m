(* ::Package:: *)
(* MicroTest test file *)

tolEq[a_, b_, tol_ : 1*^-8 ] := Max@Abs[a-b] < tol

takeUpper[mat_?SquareMatrixQ] := Extract[mat, Subsets[Range@Length[mat], {2}]]
takeLower[mat_?SquareMatrixQ] := takeUpper@Transpose[mat]

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
  IGraphM`LTemplate`ValidTemplateQ[IGraphM`Private`template]
];
Print["Template verification timing: ", t]


IGraphM`LTemplate`UnloadTemplate[IGraphM`Private`template]
t = First@AbsoluteTiming[
  IGraphM`Private`LoadIGraphM[]
];
Print["Template load timing: ", t]


nameFromUsage[symname_] :=
    With[{sym = Symbol[symname]},
      First@StringCases[sym::usage, Shortest[name__] ~~ "[" :> name]]

MT[
  AllTrue[Complement[Names["IGraphM`*"], {"IGraphM"}], nameFromUsage[#] === # &],
  True
]

MT[
  MatchQ[IGVersion[], _String],
  True
]

MT[
  IGDocumentation[]; NotebookClose@First@Notebooks[];,
  Null
]

SeedRandom[137]
ug = RandomGraph[{10,20}];
dg = RandomGraph[{10,30}, DirectedEdges->True];

umulti = Graph[UndirectedEdge @@@ RandomInteger[10, {50, 2}]];
dmulti = Graph[DirectedEdge @@@ RandomInteger[10, {50, 2}]];

uweights = RandomReal[1, EdgeCount[ug]];
wug = Graph[ug, EdgeWeight->uweights];

dweights = RandomReal[1, EdgeCount[dg]];
wdg = Graph[dg, EdgeWeight->dweights];

empty = Graph[{},{}];
edgeless = Graph[{1,2,3},{}];

ulist = Table[RandomGraph[{n, 2n}], {n, 5, 100, 5}];
dlist = Table[RandomGraph[{n, 2n}, DirectedEdges -> True], {n, 10, 100, 5}];

dolphin = ExampleData[{"NetworkGraph", "DolphinSocialNetwork"}];
web = ExampleData[{"NetworkGraph", "ExpandedComputationalGeometry"}];
collab = ExampleData[{"NetworkGraph", "CondensedMatterCollaborations2005"}];
football = ExampleData[{"NetworkGraph", "AmericanCollegeFootball"}];
lesmiserables = ExampleData[{"NetworkGraph", "LesMiserables"}];

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
MTSection["Undirected"]

MT[
  ig = IGraphM`Private`igMake[ug],
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
  ig = IGraphM`Private`igMake[dg],
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
  ig = IGraphM`Private`igMake[wug],
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
  uweights
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
  ig = IGraphM`Private`igMake[wdg],
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
  dweights
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
      IGraphM`Private`igIndexVec[ig@"edgeList"[]] == List @@@ EdgeList@IndexGraph[graph]
      ,
      Sort /@ IGraphM`Private`igIndexVec[ig@"edgeList"[]] == Sort /@ (List @@@ EdgeList@IndexGraph[graph])
    ]

cycleTestList = {ug, dg, umulti, dmulti, empty, edgeless, dolphin, web, football, collab, KaryTree[100], KaryTree[101, DirectedEdges -> True]};

MT[
  compare[IGraphM`Private`igMake[#], #],
  True
]& /@ cycleTestList

Print["Cycle timing: ", First@Timing[IGraphM`Private`igMake /@ cycleTestList;] ]


(*******************************************************************************)
MTSection["Creation"]

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


(* IGMakeLattice *)

MT[
  IGIsomorphicQ[IGMakeLattice[{3,4}], GridGraph[{3,4}]],
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
        if[empty, edgeless],
        False
      ],
      MT[
        if[edgeless, randomize@edgeless],
        True
      ],
      MT[
        if[ug, randomize[ug]],
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
        if[empty, ug],
        True
      ],
      MT[
        if[edgeless, randomize@edgeless],
        True
      ],
      MT[
        if[ug, randomize[ug]],
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
        if[Subgraph[ug, RandomSample[VertexList[ug], 5]], ug],
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
    if[dg, randomize@dg],
    True
  ],
  MT[
    if[dg, ug],
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
  BetweennessCentrality[ug] == IGBetweenness[ug, Method -> "Fast"],
  True
]

MT[
  BetweennessCentrality[dg] == IGBetweenness[dg, Method -> "Fast"],
  True
]

MT[
  BetweennessCentrality[ug] ~tolEq~ IGBetweenness[ug],
  True
]

MT[
  BetweennessCentrality[dg] ~tolEq~ IGBetweenness[dg],
  True
]

MT[
  EdgeBetweennessCentrality[ug] == 2 IGEdgeBetweenness[ug],
  True
]

MT[
  EdgeBetweennessCentrality[dg] == IGEdgeBetweenness[dg],
  True
]

MT[
  ClosenessCentrality[#] == IGCloseness[#, Normalized -> True],
  True
]& /@ {ug, dg, wug, wdg, umulti, dmulti}


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
]& /@ {ug, umulti}

MT[
  IGCliqueNumber@GraphComplement[#] == IGIndependenceNumber[#],
  True
]& /@ {ug, umulti}

MT[
  IGCliqueNumber[#] == Length@First@FindClique[#],
  True
]& /@ {ug, umulti}

MT[
  IGIndependenceNumber[#] == Length@First@FindIndependentVertexSet[#],
  True
]& /@ {ug, umulti}

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
  IGWeightedCliques[Graph[{1,2,3}, { 1<-> 2}], VertexWeight -> {0,1,2}];,
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
  Length@IGTriangles[ug],
  triangleCount[ug]
]

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
  IGVertexConnectivity[ug, 1, 5],
  VertexConnectivity[ug, 1, 5]
]

MT[
  IGEdgeConnectivity[ug, 1, 5],
  EdgeConnectivity[ug, 1, 5]
]

MT[
  IGVertexConnectivity[dg, 1, 6],
  VertexConnectivity[dg, 1, 6]
]

MT[
  IGEdgeConnectivity[dg, 1, 6],
  EdgeConnectivity[dg, 1, 6]
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
  IGVertexConnectivity[ug, 1, 2],
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
  {0} ~Join~ IGDistanceCounts[dg] == IGDistanceHistogram[dg, 1],
  True
]

MT[
  {0} ~Join~ IGDistanceCounts[ug] == IGDistanceHistogram[ug, 1]/2,
  True
]

MT[
  IGDistanceHistogram[wdg, 0.1],
  With[{vals=takeUpper@GraphDistanceMatrix[wdg] ~Join~ takeLower@GraphDistanceMatrix[wdg]},
    BinCounts[vals, {0, Ceiling[Max[vals], 0.1], 0.1}]
  ]
]

MT[
  IGDistanceHistogram[wug, 0.1],
  With[{vals=takeUpper@GraphDistanceMatrix[wug] ~Join~ takeLower@GraphDistanceMatrix[wug]},
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
  IGDistanceMatrix[ug],
  GraphDistanceMatrix[ug]
]

MT[
  IGDistanceMatrix[ug, {1,2,4}],
  GraphDistanceMatrix[ug][[{1,2,4}]]
]

MT[
  IGDistanceMatrix[ug, {1,2,4}, {3,5}],
  GraphDistanceMatrix[ug][[{1,2,4}, {3,5}]]
]

MT[
  IGDistanceMatrix[dg],
  GraphDistanceMatrix[dg]
]

MT[
  IGDistanceMatrix[dg, {1,2,4}],
  GraphDistanceMatrix[dg][[{1,2,4}]]
]

MT[
  IGDistanceMatrix[dg, {1,2,4}, {3,5}],
  GraphDistanceMatrix[dg][[{1,2,4}, {3,5}]]
]

MT[
  IGDistanceMatrix[wug] == GraphDistanceMatrix[wug],
  True
]

MT[
  IGDistanceMatrix[wug, {1,2,4}] == GraphDistanceMatrix[wug][[{1,2,4}]],
  True
]

MT[
  IGDistanceMatrix[wug, {1,2,4}, {3,4}] == GraphDistanceMatrix[wug][[{1,2,4}, {3,4}]],
  True
]

MT[
  IGDistanceMatrix[wdg] == GraphDistanceMatrix[wdg],
  True
]

MT[
  IGDistanceMatrix[wdg, {1,2,4}] == GraphDistanceMatrix[wdg][[{1,2,4}]],
  True
]

MT[
  IGDistanceMatrix[wdg, {1,2,4}, {3,4}] == GraphDistanceMatrix[wdg][[{1,2,4}, {3,4}]],
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
]& /{ug, dg, wug, wdg, umulti, dmulti}

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

MT[
  IGBipartiteQ[bipartite],
  True
]

MT[
  IGBipartitePartitions[dbipartite],
  {{1, 2, 3, 4, 5}, {6, 7, 8, 9, 10, 11, 12}}
]

MT[
  IGBipartiteQ[dbipartite],
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
MTSection["Rewiring"]

(* IGRewire *)

MT[
  VertexDegree@IGRewire[ug, 100] == VertexDegree[ug],
  True
]

MT[
  With[{r = IGRewire[dg, 100]},
    VertexInDegree[r] == VertexInDegree[dg] &&
        VertexOutDegree[r] == VertexOutDegree[dg]
  ],
  True
]

MT[
  IsomorphicGraphQ[IGRewire[ug, 0], ug],
  True
]


(* IGRewireEdges *)

MT[
  EdgeCount@IGRewireEdges[ug, 0.1] == EdgeCount[ug],
  True
]

MT[
  IsomorphicGraphQ[IGRewireEdges[ug, 0], ug],
  True
]

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

(* IGConnectNeighborhood *)

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
    IGConnectNeighborhood[ug, 1],
    ug
  ],
  True,
  {IGraphM::warning}
]

MT[
  EdgeList@IGConnectNeighborhood[
    Graph[{1 -> 2, 2 -> 1, 2 -> 3, 3 -> 2}],
    2],
  {1 \[DirectedEdge] 2, 1 \[DirectedEdge] 3, 2 \[DirectedEdge] 1,
    2 \[DirectedEdge] 3, 3 \[DirectedEdge] 1, 3 \[DirectedEdge] 2}
]

MT[
  EdgeList@IGConnectNeighborhood[
    Graph[{1 -> 2, 2 -> 3}],
    2],
  {1 \[DirectedEdge] 2, 1 \[DirectedEdge] 3, 2 \[DirectedEdge] 3}
]

MT[
  EdgeList@IGConnectNeighborhood[
    Graph[{1 -> 2, 3 -> 2}],
    2],
  {1 \[DirectedEdge] 2, 3 \[DirectedEdge] 2}
]


(* IGArticulationPoints *)

connCompCount[g_] := Length@ConnectedComponents[g]

MT[
  AllTrue[
    connCompCount@VertexDelete[ug, #] & /@ IGArticulationPoints[ug],
    # > connCompCount[ug] &
  ],
  True
]

MT[
  AllTrue[
    connCompCount@VertexDelete[dg, #] & /@ IGArticulationPoints[dg],
    # > connCompCount[dg] &
  ],
  True
]

(* IGMinSeparators *)

MT[
  IGMinSeparators[#] =!= {} & /@ ulist,
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


