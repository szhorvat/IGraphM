(* ::Package:: *)
(* MicroTest test file *)

tolEq[a_, b_, tol_ : 1*^-8 ] := Max@Abs[a-b] < tol

takeUpper[mat_?SquareMatrixQ] := Extract[mat, Subsets[Range@Length[mat], {2}]]

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


(*******************************************************************************)
MTSection["Acylic graphs"]

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
  IGCliques[Graph[{},{}]],
  {}
]

MT[
  canon@IGCliques[Graph[{1,2,3},{}]],
  {{1},{2},{3}}
]

MT[
  IGMaximalCliques[Graph[{},{}]],
  {}
]

MT[
  canon@IGMaximalCliques[Graph[{1,2,3},{}]],
  {{1},{2},{3}}
]

MT[
  canon@IGCliques[Graph[{1, 2}, {1 -> 2}]],
  {{1},{2},{1,2}},
  {IGraphM::warning}
]

MT[
  canon@IGMaximalCliques[Graph[{1, 2}, {1 -> 2}]],
  {{1,2}},
  {IGraphM::warning}
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
  IGMotifsTotalCount[dolphin, 3],
  733
]

MT[
  IGMotifsTotalCount[dolphin, 4],
  3800
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
  MapThread[SameQ, {IGTriadCensus[web] /@ Keys@IGData["MANTriadLabels"], IGMotifs[web, 3]}],
  {False, False, True, False, True, True, True, True, True, True, True, True, True, True, True, True}
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

(* directed, unweighted, unconnected *)
MT[
  IGDistanceCounts[web],
  hist = Values@Rest@Most@KeySort@Counts@Flatten@GraphDistanceMatrix[web]
]

MT[
  IGAveragePathLength[web],
  N@Mean@WeightedData[Range@Length[hist], hist]
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


(*******************************************************************************)
MTSection["Shortest paths"]

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


(*******************************************************************************)
MTSection["Test remaining functions"]

MT[
  IGData /@ IGData[];,
  Null
]

MT[
  IGIsomorphicQ[IGMakeLattice[{3,4}], GridGraph[{3,4}]],
  True
]

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

MT[
  EdgeCount@IGRewireEdges[ug, 0.1] == EdgeCount[ug],
  True
]

MT[
  IsomorphicGraphQ[IGRewireEdges[ug, 0], ug],
  True
]

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

MT[
  IGMinSeparators[#] =!= {} & /@ ulist,
  ConnectedGraphQ /@ ulist
]

MT[
  And @@ MapThread[IGIsomorphicQ, {IGData[{"AllDirectedGraphs", 3}], Values@IGData["MANTriadLabels"]}],
  True
]