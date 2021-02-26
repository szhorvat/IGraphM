(* ::Package:: *)

(* ::Title:: *)
(*IGraph/M tests*)


(* ::Text:: *)
(*This is a MicroTest test file. See https://github.com/szhorvat/MicroTest*)


(* ::Section::Closed:: *)
(*Utility functions for testing*)


tolEq[a_, b_, tol_ : 1*^-8 ] := Max@Abs[a-b] < tol


takeNonDiag[mat_] := IGraphM`IGTakeUpper[mat] ~Join~ IGraphM`IGTakeLower[mat]

(* In M12.1 and later, UndirectedEdge can have 3 arguments, so we cannot canonicalize simply with Orderless. *)
If[$VersionNumber >= 12.1,
  sameGraphQ[g1_, g2_] :=
      Block[{UndirectedEdge},
        UndirectedEdge[a_, b_, rest___] /; Not@OrderedQ[{a, b}] := UndirectedEdge[b, a, rest];
        Sort@VertexList[g1] === Sort@VertexList[g2] && Sort@EdgeList[g1] === Sort@EdgeList[g2]
      ]
  ,
  sameGraphQ[g1_, g2_] :=
      Block[{UndirectedEdge},
        SetAttributes[UndirectedEdge, Orderless];
        Sort@VertexList[g1] === Sort@VertexList[g2] && Sort@EdgeList[g1] === Sort@EdgeList[g2]
      ]
]


(* RelationGraph[] is not available in 10.0 *)
relationGraph[f_, list_, opt : OptionsPattern[]]:=
    AdjacencyGraph[list, Boole@Outer[f, list, list, 1], opt]


(* ::Subsubsection::Closed:: *)
(*samePropGraphQ*)


(* ::Text:: *)
(*Check that two graphs are the same and have the same properties and property values.*)


ruleQ[_Rule | _RuleDelayed] := True
ruleQ[_] = False;


rulePropQ[{__?ruleQ}] := True;
rulePropQ[_] = False;


defRulePropQ[{_, __?ruleQ}] := True;
defRulePropQ[_] = False;


normalizeGraphOpt[opt_] :=
    MapAt[
      Which[
        rulePropQ[#], Sort[#],
        defRulePropQ[#], Prepend[Sort@Rest[#], First[#]],
        True, #
      ] &,
      opt,
      {All, 2}
    ]


samePropGraphQ[g1_, g2_] :=
    Internal`InheritedBlock[{UndirectedEdge},
      SetAttributes[UndirectedEdge, Orderless];
      VertexList[g1] === VertexList[g2] &&
          EdgeList[g1] === EdgeList[g2] &&
          normalizeGraphOpt@Options[g1] === normalizeGraphOpt@Options[g2]
    ]


(* ::Section::Closed:: *)
(*Test graphs*)


(* ::Text:: *)
(*Set a fixed seed for reproducible testing.*)


SeedRandom[137]


(* These serve to silence undefined symbol warnings in IDEA *)
a::usage = "a is a vertex name";
b::usage = "b is a vertex name";


ugs = RandomGraph[{10,20}];
dgs = RandomGraph[{10,30}, DirectedEdges->True];


ugi = GraphComputation`ToGraphRepresentation[RandomGraph[{12,25}], "Incidence"];
dgi = GraphComputation`ToGraphRepresentation[RandomGraph[{12,25}, DirectedEdges -> True], "Incidence"];


umulti = Graph[UndirectedEdge @@@ RandomInteger[10, {50, 2}]];
dmulti = Graph[DirectedEdge @@@ RandomInteger[10, {50, 2}]];


umulti2 = Graph[UndirectedEdge @@@ RandomInteger[10, {60, 2}]];
dmulti2 = Graph[DirectedEdge @@@ RandomInteger[10, {60, 2}]];
If[$VersionNumber >= 12.1,
	umulti2 = EdgeTaggedGraph[umulti2];
	dmulti2 = EdgeTaggedGraph[dmulti2];
]


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
friendship = ExampleData[{"NetworkGraph","Friendship"}];


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


(* These graphs have no non-trivial automorphisms, and therefore asymmetric. *)
(* Exclude the null and singleton graphs because they are symmetric. *)
asymmList = {
  Graph[{1 <-> 2, 2 <-> 3, 3 <-> 4, 2 <-> 4, 4 <-> 5, 5 <-> 6}],
  GraphData["FruchtGraph"],
  GraphData[{6,69}],
  GraphData[{6,95}],
  GraphData[{"Tree",{7,7}}]
};


dmultiloopy = 
	Graph[
		Range[11],
		{1->2,2->11,11->3,3->4,4->1,1->5,5->6,7->2,4->4,4->4,4->4,2->11,2->8,8->3,3->9,9->9,5->10}
	];


umultiloopy = 
    Graph[
        Range[10],
        {1 <-> 2, 2 <-> 3, 3 <-> 4, 4 <-> 1, 1 <-> 4, 1 <-> 5, 5 <-> 6, 3 <-> 7, 7 <-> 8, 5 <-> 9, 9 <-> 10, 10 <-> 7, 6 <-> 6, 6 <-> 6, 9 <-> 9, 9 <-> 3, 3 <-> 9, 9 <-> 3}
    ];


(* ::Section::Closed:: *)
(*Sanity checks for Mathematica built-ins*)


MTSection["Sanity checks for Mathematica built-ins"]


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


(* ::Subsubsection::Closed:: *)
(*Verify WeightValues and WeightVector*)


(* ::Text:: *)
(*Verify GraphComputation`WeightValues*)


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


(* ::Text:: *)
(*Similar checks for GraphComputation`WeightVector, which returns the vertex weight vector.*)


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


(* ::Subsubsection::Closed:: *)
(*Tests for IncidenceMatrix*)


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
If[$VersionNumber < 12.1, (* behaviour has been updated to match the docs in 12.1 *)

MT[
  Normal@IncidenceMatrix[Graph[{a, b}, {b -> b}]],
  {{0}, {2}}
];

MT[
  Normal@IncidenceMatrix@Graph[Range[10], {DirectedEdge[6, 9], DirectedEdge[3, 7],
    DirectedEdge[3, 1], DirectedEdge[7, 7], DirectedEdge[3, 1],
    DirectedEdge[9, 2], DirectedEdge[3, 4], DirectedEdge[2, 9],
    DirectedEdge[10, 5], DirectedEdge[10, 4], DirectedEdge[6, 6],
    DirectedEdge[5, 5], DirectedEdge[1, 5], DirectedEdge[5, 9],
    DirectedEdge[10, 9]}],
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

] (* end If *)

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


(* ::Subsubsection::Closed:: *)
(*Property canonical form*)


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


(* ::Subsubsection::Closed:: *)
(*WeightedAdjacencyMatrix*)


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


(* ::Subsubsection::Closed:: *)
(*UpperTriangularMatrixToVector*)


(* Statistics`Library`UpperTriangularMatrixToVector is used in IGTakeUpper in M \[GreaterEqual] 10.4.
   Verify that this undocumented symbol exists and that it works. *)

MT[
  Names["Statistics`Library`UpperTriangularMatrixToVector"],
  {"Statistics`Library`UpperTriangularMatrixToVector"}
]

MT[
  Statistics`Library`UpperTriangularMatrixToVector@Partition[Range[16], 4],
  {2, 3, 4, 7, 8, 12}
]


(* ::Section::Closed:: *)
(*Basic tests*)


MTSection["Basic"]


t = First@AbsoluteTiming[
  MT[
    Print@Get["IGraphM`"],
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
      First@StringCases[sym::usage, Shortest[name__] ~~ "[" :> name]
    ]

MT[
  AllTrue[Complement[Names["IGraphM`*"], {"IGraphM", "MultiEdges", "$IGExportFormats", "$IGImportFormats"}], nameFromUsage[#] === # &],
  True
]


MT[
  MatchQ[IGVersion[], _String],
  True
]


(* ::Text:: *)
(*Only test IGDocumentation when running with notebooks.*)


If[$Notebooks,
  MT[
    IGDocumentation[]; NotebookClose /@ Select[Notebooks[], CurrentValue[#, WindowTitle] === "IGraph/M Documentation"&];,
    Null
  ]
]


(* ::Text:: *)
(*Seed the igraph random number generator for reproducible testing.*)


MT[
  IGSeedRandom[137],
  Null
]


(* ::Section::Closed:: *)
(*Undirected*)


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


(* ::Section::Closed:: *)
(*Directed*)


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


(* ::Section::Closed:: *)
(*Weighted undirected*)


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


(* ::Section::Closed:: *)
(*Weighted directed*)


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


(* ::Section::Closed:: *)
(*Cycling graphs through igraph*)


MTSection["Cycling graphs through igraph"]


(* We discard the 3rd argument of edge expressions, which is present in EdgeTaggedGraphs in M12.1 *)
compare[ig_, graph_] :=
    If[DirectedGraphQ[graph],
      IGraphM`PackageScope`igIndexVec[ig@"edgeList"[]] == List @@@ (EdgeList@IndexGraph[graph])[[All, {1,2}]]
      ,
      Sort /@ IGraphM`PackageScope`igIndexVec[ig@"edgeList"[]] == Sort /@ (List @@@ EdgeList@IndexGraph[graph])[[All, {1,2}]]
    ]


cycleTestList = {
  ugs, dgs, ugi, dgi, umulti, dmulti, umulti2, dmulti2,
  empty, edgeless,
  dolphin, web, football, collab,
  KaryTree[100], KaryTree[101, DirectedEdges -> True]
};

MT[
  IGraphM`PackageScope`igMake[#], 
  #,
  SameTest -> compare
]& /@ cycleTestList


MT[
  IGraphM`PackageScope`igMake[#], 
  #,
  SameTest -> compare
]& /@ Join[IGData[{"AllDirectedGraphs", 3}], IGData[{"AllDirectedGraphs", 4}]]

MT[
  IGraphM`PackageScope`igMake[#], 
  #,
  SameTest -> compare
]& /@ Join[IGData[{"AllUndirectedGraphs", 3}], IGData[{"AllUndirectedGraphs", 4}]]


Print["Cycle timing: ", First@Timing[IGraphM`PackageScope`igMake /@ cycleTestList;] ]


(* ::Section::Closed:: *)
(*Creation: deterministic*)


MTSection["Creation: deterministic"]


(* ::Subsubsection::Closed:: *)
(*IGShorthand*)


MT[
  IGShorthand[""],
  IGEmptyGraph[],
  SameTest -> IGSameGraphQ
]

MT[
  IGShorthand["1"],
  IGEmptyGraph[1],
  SameTest -> IGSameGraphQ
]

MT[
  IGShorthand["1, 2"],
  IGEmptyGraph[2],
  SameTest -> IGSameGraphQ
]

MT[
  IGShorthand["a,b"],
  Graph[{"a", "b"}, {}],
  SameTest -> IGSameGraphQ
]

MT[
  VertexList@IGShorthand["2,1,1-2"],
  {2, 1}
]

MT[
  IGShorthand["1-2"],
  Graph[{1<->2}],
  SameTest -> IGSameGraphQ
]

MT[
  IGShorthand["1 -- 2"],
  Graph[{1<->2}],
  SameTest -> IGSameGraphQ
]

MT[
  IGShorthand["1->2"],
  Graph[{1->2}],
  SameTest -> IGSameGraphQ
]

MT[
  IGShorthand["1-->2"],
  Graph[{1->2}],
  SameTest -> IGSameGraphQ
]

MT[
  IGShorthand["1<-2"],
  Graph[{2->1}],
  SameTest -> IGSameGraphQ
]

MT[
  IGShorthand["1 <--- 2"],
  Graph[{2->1}],
  SameTest -> IGSameGraphQ
]

MT[
  IGShorthand["1:2:3:4 - 1:2:3:4"],
  CompleteGraph[4],
  SameTest -> IGSameGraphQ
]

MT[
  IGShorthand["1:2:3:4 - 5:6:7"],
  CompleteGraph[{4, 3}],
  SameTest -> IGSameGraphQ
]

MT[
  IGShorthand["1----2 -3--4 - 1"],
  CycleGraph[4],
  SameTest -> IGSameGraphQ
]

MT[
  IGShorthand["1-2-1"],
  Graph[{1 <-> 2}, VertexLabels -> "Name"]
]

MT[
  IGShorthand["1:2:3 - 3:2:1", SelfLoops -> True],
  Graph[{1<->2, 2<->3, 3<->1, 1<->1, 2<->2, 3<->3}],
  SameTest -> IGSameGraphQ
]

MT[
  IGShorthand["1:2:3 - 3:2:1", SelfLoops -> True],
  IGCompleteGraph[3, SelfLoops -> True],
  SameTest -> IGSameGraphQ
]

MT[
  IGShorthand["1:2:3 - 3:2:1", SelfLoops -> True, MultiEdges -> True],
  Graph[{1 <-> 3, 1 <-> 2, 1 <-> 1, 2 <-> 3, 2 <-> 2, 2 <-> 1, 3 <-> 3, 3 <-> 2, 3 <-> 1}],
  SameTest -> IGSameGraphQ
]

MT[
  IGShorthand["1:2:3 - 3:2:1", MultiEdges -> True],
  Graph[{1 <-> 3, 1 <-> 2, 2 <-> 3, 1 <-> 2, 2 <-> 3, 1 <-> 3}],
  SameTest -> IGSameGraphQ
]


(* ::Subsubsection::Closed:: *)
(*IGEmptyGraph*)


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


(* ::Subsubsection::Closed:: *)
(*IGLCF*)


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


(* ::Subsubsection::Closed:: *)
(*IGChordalRing*)


MT[
  IGChordalRing[15, {{3, 12, 4}, {7, 8, 11}}],
  Graph[{1 -> 2, 2 -> 3, 3 -> 4, 4 -> 5, 5 -> 6, 6 -> 7, 7 -> 8, 8 -> 9,
    9 -> 10, 10 -> 11, 11 -> 12, 12 -> 13, 13 -> 14, 14 -> 15, 1 -> 15,
    1 -> 4, 1 -> 8, 2 -> 14, 2 -> 10, 3 -> 7, 3 -> 14, 4 -> 7, 4 -> 11,
    2 -> 5, 5 -> 13, 6 -> 10, 2 -> 6, 7 -> 10, 7 -> 14, 5 -> 8, 1 -> 8,
    9 -> 13, 5 -> 9, 10 -> 13, 2 -> 10, 8 -> 11, 4 -> 11, 1 -> 12,
    8 -> 12, 1 -> 13, 5 -> 13, 11 -> 14, 7 -> 14, 4 -> 15, 11 -> 15}, DirectedEdges -> False],
  SameTest -> IGSameGraphQ
]

MT[
  IGChordalRing[15, {{3, 12, 4, 3}, {7, 8, 11, 3}}];,
  Null,
  {IGraphM::error}
]


MT[
  IGChordalRing[5, {}],
  Graph[{1, 2, 3, 4, 5}, {1 <-> 2, 1 <-> 5, 2 <-> 3, 3 <-> 4, 4 <-> 5}],
  SameTest -> IGSameGraphQ
]


MT[
  IGChordalRing[6, {2, 3, 1}],
  Graph[{1 <-> 2, 2 <-> 3, 3 <-> 4, 4 <-> 5, 5 <-> 6, 1 <-> 6, 1 <-> 3, 2 <-> 5, 3 <-> 4, 4 <-> 6, 2 <-> 5, 1 <-> 6}],
  SameTest -> IGSameGraphQ
]


(* ::Subsubsection::Closed:: *)
(*IGSquareLattice*)


MT[
  IGSquareLattice[{3,4}],
  GridGraph[{3,4}],
  SameTest -> IGSameGraphQ
]

MT[
  PathGraphQ@IGSquareLattice[{5}],
  True
]

MT[
  IGSquareLattice[{5}, "Periodic" -> True],
  CycleGraph[5],
  SameTest -> IGSameGraphQ
]

MT[
  IGSquareLattice[{9}, "Periodic" -> True, "Radius" -> 2],
  GraphPower[CycleGraph[9], 2],
  SameTest -> IGSameGraphQ
]

MT[
  IGSquareLattice[{2, 3}, DirectedEdges -> True, "Mutual" -> True],
  DirectedGraph[GridGraph[{2, 3}]],
  SameTest -> IGSameGraphQ
]

MT[
  IGSquareLattice[{2, 3, 4}, "Periodic" -> True],
  Graph[{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24}, {1 <-> 2, 1 <-> 3, 1 <-> 7, 2 <-> 4, 2 <-> 8, 3 <-> 4, 3 <-> 5, 3 <-> 9, 4 <-> 6, 4 <-> 10, 5 <-> 6, 1 <-> 5, 5 <-> 11, 2 <-> 6, 6 <-> 12, 7 <-> 8, 7 <-> 9, 7 <-> 13, 8 <-> 10, 8 <-> 14, 9 <-> 10, 9 <-> 11, 9 <-> 15, 10 <-> 12, 10 <-> 16, 11 <-> 12, 7 <-> 11, 11 <-> 17, 8 <-> 12, 12 <-> 18, 13 <-> 14, 13 <-> 15, 13 <-> 19, 14 <-> 16, 14 <-> 20, 15 <-> 16, 15 <-> 17, 15 <-> 21, 16 <-> 18, 16 <-> 22, 17 <-> 18, 13 <-> 17, 17 <-> 23, 14 <-> 18, 18 <-> 24, 19 <-> 20, 19 <-> 21, 1 <-> 19, 20 <-> 22, 2 <-> 20, 21 <-> 22, 21 <-> 23, 3 <-> 21, 22 <-> 24, 4 <-> 22, 23 <-> 24, 19 <-> 23, 5 <-> 23, 20 <-> 24, 6 <-> 24}],
  SameTest -> IGSameGraphQ
]

(* TODO is this reasonable? *)
MT[
  IGSquareLattice[{}],
  Graph[{1}, {}],
  SameTest -> IGSameGraphQ
]

MT[
  IGSquareLattice[{1, 1, 1}],
  Graph[{1}, {}],
  SameTest -> IGSameGraphQ
]

MT[
  IGSquareLattice[{0}],
  Graph[{}],
  SameTest -> IGSameGraphQ
]

MT[
  IGSquareLattice[{1, 0, 2}],
  Graph[{}],
  SameTest -> IGSameGraphQ
]


(* ::Subsubsection::Closed:: *)
(*IGTriangularLattice*)


MT[
  IGTriangularLattice[0],
  IGEmptyGraph[0]
]

MT[
  IGTriangularLattice[1],
  IGEmptyGraph[1],
  SameTest -> IGSameGraphQ
]

MT[
  IGTriangularLattice[{0,0}],
  IGEmptyGraph[0]
]

MT[
  IGTriangularLattice[{0,1}],
  IGEmptyGraph[0]
]

MT[
  IGTriangularLattice[{1,0}],
  IGEmptyGraph[0]
]

MT[
  IGTriangularLattice[{1,1}],
  IGEmptyGraph[1],
  SameTest -> IGSameGraphQ
]

MT[
  IGTriangularLattice[{1, 5}],
  PathGraph@Range[5],
  SameTest -> IGSameGraphQ
]

MT[
  IGTriangularLattice[{5,1}],
  PathGraph@Range[5],
  SameTest -> IGSameGraphQ
]

MT[
  IGTriangularLattice[4],
  Graph[{1 <-> 2, 1 <-> 3, 2 <-> 3, 2 <-> 4, 2 <-> 5, 4 <-> 5, 3 <-> 5,
   3 <-> 6, 5 <-> 6, 4 <-> 7, 4 <-> 8, 7 <-> 8, 5 <-> 8, 5 <-> 9,
   8 <-> 9, 6 <-> 9, 6 <-> 10, 9 <-> 10}],
  SameTest -> IGSameGraphQ
]

MT[
  IGTriangularLattice[{3,4}],
  Graph[{1 <-> 4, 2 <-> 5, 3 <-> 6, 4 <-> 7, 5 <-> 8, 6 <-> 9, 7 <-> 10,
   8 <-> 11, 9 <-> 12, 1 <-> 2, 2 <-> 3, 4 <-> 5, 5 <-> 6, 7 <-> 8,
   8 <-> 9, 10 <-> 11, 11 <-> 12, 2 <-> 4, 2 <-> 6, 5 <-> 7, 5 <-> 9,
   8 <-> 10, 8 <-> 12}],
  SameTest -> IGSameGraphQ
]

MT[
  IGTriangularLattice[{4,3}],
  Graph[{1 <-> 5, 2 <-> 6, 3 <-> 7, 4 <-> 8, 5 <-> 9, 6 <-> 10, 7 <-> 11,
   8 <-> 12, 1 <-> 2, 2 <-> 3, 3 <-> 4, 5 <-> 6, 6 <-> 7, 7 <-> 8,
   9 <-> 10, 10 <-> 11, 11 <-> 12, 2 <-> 5, 2 <-> 7, 4 <-> 7, 6 <-> 9,
   6 <-> 11, 8 <-> 11}],
  SameTest -> IGSameGraphQ
]

MT[
  IGTriangularLattice[4, "Periodic" -> True],
  $Failed,
  {IGTriangularLattice::pimp}
]

MT[
  IGTriangularLattice[{2, 3}, "Periodic" -> True],
  Graph[{1, 2, 3, 4, 5, 6}, {1 <-> 3, 2 <-> 4, 3 <-> 5, 4 <-> 6, 5 <-> 1, 6 <-> 2, 1 <-> 2, 2 <-> 1, 3 <-> 4, 4 <-> 3, 5 <-> 6, 6 <-> 5, 2 <-> 3, 2 <-> 3, 4 <-> 5, 4 <-> 5, 6 <-> 1, 6 <-> 1}],
  SameTest -> IGSameGraphQ
]


(* ::Subsubsection::Closed:: *)
(*IGKaryTree*)


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


(* ::Subsubsection::Closed:: *)
(*IGSymmetricTree*)


MT[
  IGSymmetricTree[{2, 3, 1}],
  Graph[{1 <-> 2, 1 <-> 3, 2 <-> 4, 2 <-> 5, 2 <-> 6, 3 <-> 7, 3 <-> 8,
   3 <-> 9, 4 <-> 10, 5 <-> 11, 6 <-> 12, 7 <-> 13, 8 <-> 14, 9 <-> 15}],
  SameTest -> IGSameGraphQ
]

MT[
  IGSymmetricTree[{}],
  IGEmptyGraph[1],
  SameTest -> IGSameGraphQ
]

MT[
  IGSymmetricTree[{5}],
  StarGraph[6],
  SameTest -> IGSameGraphQ
]

(* must stay unevaluated *)
MT[
  Head@IGSymmetricTree[{0}],
  IGSymmetricTree
]

(* must stay unevaluated *)
MT[
  Head@IGSymmetricTree[{-1}],
  IGSymmetricTree
]


(* ::Subsubsection::Closed:: *)
(*IGBetheLattice*)


MT[
  IGBetheLattice[3],
  Graph[{1 <-> 2, 1 <-> 3, 1 <-> 4, 2 <-> 5, 2 <-> 6, 3 <-> 7, 3 <-> 8, 4 <-> 9, 4 <-> 10}],
  SameTest -> IGSameGraphQ
]

MT[
  IGBetheLattice[3, 4],
  Graph[{1 <-> 2, 1 <-> 3, 1 <-> 4, 1 <-> 5, 2 <-> 6, 2 <-> 7, 2 <-> 8,
   3 <-> 9, 3 <-> 10, 3 <-> 11, 4 <-> 12, 4 <-> 13, 4 <-> 14, 5 <-> 15,
   5 <-> 16, 5 <-> 17}],
  SameTest -> IGSameGraphQ
]


(* ::Subsubsection::Closed:: *)
(*IGCompleteGraph*)


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


MT[
  IGCompleteGraph[5],
  CompleteGraph[5],
  SameTest -> IGSameGraphQ
]

MT[
  IGCompleteGraph[6, DirectedEdges -> True],
  CompleteGraph[6, DirectedEdges -> True],
  SameTest -> IGSameGraphQ
]

MT[
  IGCompleteGraph[0],
  IGEmptyGraph[],
  SameTest -> IGSameGraphQ
]

MT[
  IGCompleteGraph[{}],
  IGEmptyGraph[],
  SameTest -> IGSameGraphQ
]

MT[
  IGCompleteGraph[4, SelfLoops -> True],
  Graph[{1 <-> 1, 1 <-> 2, 1 <-> 3, 1 <-> 4, 2 <-> 2, 2 <-> 3, 2 <-> 4, 3 <-> 3, 3 <-> 4, 4 <-> 4}],
  SameTest -> IGSameGraphQ
]

MT[
  VertexList[IGCompleteGraph[{"a", "b", "c", "d"}]],
  {"a", "b", "c", "d"}
]

MT[
  IGCompleteGraph[{"a","b","c","d"}],
  CompleteGraph[4],
  SameTest -> IGIsomorphicQ
]


(* ::Subsubsection::Closed:: *)
(*IGCompleteAcyclicGraph*)


MT[
  IGIsomorphicQ[
    IGCompleteAcyclicGraph[4],
    DirectedGraph[CompleteGraph[4], "Acyclic"]
  ],
  True
]

MT[
  IGIsomorphicQ[
    IGCompleteAcyclicGraph[{"a","b","c","d"}],
    DirectedGraph[CompleteGraph[4], "Acyclic"]
  ],
  True
]

MT[
  Sort[VertexList[IGCompleteAcyclicGraph[{"a", "b", "c"}]]],
  {"a", "b", "c"}
]

MT[
  IGCompleteAcyclicGraph[{}],
  Graph[{}],
  SameTest -> IGSameGraphQ
]


(* ::Subsubsection::Closed:: *)
(*IGKautzGraph*)


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


(* ::Subsubsection::Closed:: *)
(*IGDeBruijnGraph*)


MT[
  IGDeBruijnGraph[3, 4],
  DeBruijnGraph[3, 4],
  SameTest -> IGSameGraphQ
]


(* ::Subsubsection::Closed:: *)
(*IGGraphAtlas*)


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


(* ::Subsubsection::Closed:: *)
(*IGFromNauty*)


(* ::Text:: *)
(*Graph6*)


MT[
  IGFromNauty["?"],
  Graph[{}],
  SameTest -> IGSameGraphQ
]

MT[
  IGFromNauty["@"],
  Graph[{1}, {}],
  SameTest -> IGSameGraphQ
]

MT[
  IGFromNauty["A?"],
  Graph[{1, 2}, {}],
  SameTest -> IGSameGraphQ
]

MT[
  IGFromNauty["A_"],
  Graph[{1 <-> 2}],
  SameTest -> IGSameGraphQ
]

MT[
  IGFromNauty["EEz_"],
  Graph[{1, 2, 3, 4, 5, 6}, {1 <-> 4, 2 <-> 4, 1 <-> 5, 2 <-> 5, 3 <-> 5, 1 <-> 6, 2 <-> 6, 3 <-> 6}],
  SameTest -> IGSameGraphQ
]


(* ::Text:: *)
(*Sparse6*)


MT[
  IGFromNauty[":?"],
  Graph[{}],
  SameTest -> IGSameGraphQ
]

MT[
  IGFromNauty[":@"],
  Graph[{1}, {}],
  SameTest -> IGSameGraphQ
]

MT[
  IGFromNauty[":A"],
  Graph[{1, 2}, {}],
  SameTest -> IGSameGraphQ
]

MT[
  IGFromNauty[":An"],
  Graph[{1 <-> 2}],
  SameTest -> IGSameGraphQ
]

MT[
  IGFromNauty[":Ek@_Q_Q"],
  Graph[{1, 2, 3, 4, 5, 6}, {1 <-> 4, 2 <-> 4, 1 <-> 5, 2 <-> 5, 3 <-> 5, 1 <-> 6, 2 <-> 6, 3 <-> 6}],
  SameTest -> IGSameGraphQ
]


(* ::Text:: *)
(*Digraph6*)


MT[
  IGFromNauty["&?"],
  Graph[{}],
  SameTest -> IGSameGraphQ
]

MT[
  IGFromNauty["&@?"],
  Graph[{1}, {}],
  SameTest -> IGSameGraphQ
]

MT[
  IGFromNauty["&A?"],
  Graph[{1, 2}, {}],
  SameTest -> IGSameGraphQ
]

MT[
  IGFromNauty["&AO"],
  Graph[{1 -> 2}],
  SameTest -> IGSameGraphQ
]

MT[
  IGFromNauty["&AW"],
  Graph[{1 -> 2, 2 -> 1}],
  SameTest -> IGSameGraphQ
]

MT[
  IGFromNauty["&COx_"],
  Graph[{1 -> 2, 2 -> 3, 2 -> 4, 3 -> 1, 3 -> 4, 4 -> 1}],
  SameTest -> IGSameGraphQ
]


(* ::Text:: *)
(*Bad input*)


MT[
  IGFromNauty["xx"],
  $Failed,
  {IGraphM::error}
]


(* ::Subsubsection::Closed:: *)
(*IGExpressionTree*)


(* TODO *)


MT[
  IGExpressionTree[{{1}, {2, 3}, 4}],
  Graph[{{1, 1}, {1}, {2, 1}, {2, 2}, {2}, {3}, {}}, {{1} -> {1, 1}, {} -> {1}, {2} -> {2, 1}, {2} -> {2, 2}, {} -> {2}, {} -> {3}}],
  SameTest -> IGSameGraphQ
]

MT[
  IGExpressionTree[{{1}, {2, 3}, 4}, DirectedEdges -> False],
  Graph[{{1, 1}, {1}, {2, 1}, {2, 2}, {2}, {3}, {}}, {{1} <-> {1, 1}, {} <-> {1}, {2} <-> {2, 1}, {2} <-> {2, 2}, {} <-> {2}, {} <-> {3}}],
  SameTest -> IGSameGraphQ
]


(* ::Section::Closed:: *)
(*Creation: random*)


MTSection["Creation: random"]


(* ::Subsubsection::Closed:: *)
(*IGDegreeSequenceGame*)


MT[
  IGIsomorphicQ[
    IGDegreeSequenceGame[{0,0,0,0}, Method -> #],
    Graph[{1,2,3,4},{}]
  ],
  True
]& /@ { "SimpleNoMultiple", "Simple" }

MT[
  IGIsomorphicQ[
    IGDegreeSequenceGame[{1,1,0,0}, Method -> #],
    Graph[{1,2,3,4},{1<->2}]
  ],
  True
]& /@ { "SimpleNoMultiple", "Simple" }

MT[
  IGIsomorphicQ[
    IGDegreeSequenceGame[{2,2,2,2}, Method -> #],
    CycleGraph[4]
  ],
  True
]& /@ { "SimpleNoMultiple", "VigerLatapy" }


(* ::Subsubsection::Closed:: *)
(*IGKRegularGame*)


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


(* ::Subsubsection::Closed:: *)
(*IGBarabasiAlbertGame*)


MT[
  With[{g = IGBarabasiAlbertGame[100, 1]}, {VertexCount[g], EdgeCount[g]}],
  {100, 99}
]

MT[
  With[{g = IGBarabasiAlbertGame[100, 2]}, {VertexCount[g], EdgeCount[g]}],
  {100, 197}
]


MT[
  IGTreeQ[IGBarabasiAlbertGame[10, 1], "In"],
  True
]


MT[
  DirectedGraphQ[IGBarabasiAlbertGame[10, 1, DirectedEdges -> #]],
  #
]& /@ {True, False}


MT[
  IGBarabasiAlbertGame[#,1],
  IGEmptyGraph[#]
]& /@ {0, 1}


(* ::Subsubsection::Closed:: *)
(*IGGrowingGame*)


MT[
  With[{g = IGGrowingGame[100, 2]}, {VertexCount[g], EdgeCount[g]}],
  {100, 198}
]


(* ::Subsubsection::Closed:: *)
(*IGGeometricGame*)


MT[
  VertexCount@IGGeometricGame[100, 0.1],
  100
]


(* ::Subsubsection::Closed:: *)
(*IGTreeGame*)


MT[
  IGTreeQ@IGTreeGame[12],
  True
]    


MT[
  VertexCount@IGTreeGame[67],
  67
]


MT[
  IGTreeGame[#],
  IGCompleteGraph[#],
  SameTest -> IGSameGraphQ
]& /@ Range[0,2]


MT[
  IGTreeGame[3],
  Graph[{1 <-> 2, 2 <-> 3}],
  SameTest -> IGIsomorphicQ (* all 3-vertex trees are isomorphic *)
]


MT[
  Table[IGTreeGame[5],{1000}]//DeleteDuplicatesBy[AdjacencyMatrix]//Length,
  125
]


MT[
  Table[IGTreeGame[4,Method->"PruferCode"],{500}]//DeleteDuplicatesBy[AdjacencyMatrix]//Length,
  16
]


MT[
  DirectedGraphQ@IGTreeGame[5,DirectedEdges->True],
  True
]


MT[
  IGTreeQ[IGTreeGame[5, DirectedEdges -> True]],
  True
]


MT[
  IGTreeGame[4, Method -> "PruferCode", DirectedEdges -> True],
  $Failed,
  {IGraphM::error}
]


(* ::Section::Closed:: *)
(*Similarity measures*)


MTSection["Similarity measures"]


(* ::Subsubsection::Closed:: *)
(*IGCocitationCoupling*)


cocitCoup[g_] :=
	With[{am = AdjacencyMatrix[g]},
		Normal[Transpose[am] . am] - DiagonalMatrix@VertexInDegree[g]
	]


MT[
  IGCocitationCoupling[#],
  cocitCoup[#]
]& /@ ulist

MT[
  IGCocitationCoupling[#],
  cocitCoup[#]
]& /@ dlist

MT[
  IGCocitationCoupling[dgs],
  cocitCoup[dgs]
]

MT[
  IGCocitationCoupling[dmulti],
  cocitCoup[dmulti]
]

MT[
  IGCocitationCoupling[dmulti2],
  cocitCoup[dmulti2]
]


MT[
  IGCocitationCoupling[IGEmptyGraph[]],
  {}
]

MT[
  IGCocitationCoupling[IGEmptyGraph[1]],
  {{0}}
]


(* ::Subsubsection::Closed:: *)
(*IGBibliographicCoupling*)


bibCoup[g_] :=
	With[{am = AdjacencyMatrix[g]},
		Normal[am . Transpose[am] - DiagonalMatrix@VertexOutDegree[g]]
	]


MT[
  IGBibliographicCoupling[#],
  bibCoup[#]
]& /@ ulist

MT[
  IGBibliographicCoupling[#],
  bibCoup[#]
]& /@ dlist

MT[
  IGBibliographicCoupling[dgs],
  bibCoup[dgs]
]

MT[
  IGBibliographicCoupling[dmulti],
  bibCoup[dmulti]
]

MT[
  IGBibliographicCoupling[dmulti2],
  bibCoup[dmulti2]
]


MT[
  IGBibliographicCoupling[IGEmptyGraph[]],
  {}
]

MT[
  IGBibliographicCoupling[IGEmptyGraph[1]],
  {{0}}
]


(* ::Section::Closed:: *)
(*Random walks*)


MTSection["Random walks"]


rwg=IGShorthand["1->2->3->4->2"]


MT[
  Union[IGRandomWalk[rwg, 2, 100]],
  {2, 3, 4}
]

MT[
  Union[IGRandomWalk[rwg, 1, 100]],
  {1, 2, 3, 4}
]

MT[
  Union[IGRandomEdgeIndexWalk[rwg, 2, 100]],
  {2, 3, 4}
]

MT[
  Union[IGRandomEdgeIndexWalk[rwg, 1, 100]],
  {1, 2, 3, 4}
]

MT[
  Union[IGRandomEdgeWalk[rwg, 2, 100]],
  {DirectedEdge[2, 3], DirectedEdge[3, 4], DirectedEdge[4, 2]}
]

MT[
  Union[IGRandomEdgeWalk[rwg, 1, 100]],
  {DirectedEdge[1, 2], DirectedEdge[2, 3], DirectedEdge[3, 4], DirectedEdge[4, 2]}
]

MT[
  Union[IGRandomWalk[Graph[{"a" <-> "b"}], "b", 100]],
  {"a", "b"}
]


MT[
  IGRandomWalk[Graph[{"a" <-> "b"}], "x", 100],
  $Failed,
  {VertexIndex::inv}
]


MT[
  Length[IGRandomWalk[rwg, 1, 100]],
  100
]

MT[
  Length[IGRandomEdgeWalk[rwg, 1, 100]],
  100
]

MT[
  Length[IGRandomEdgeIndexWalk[rwg, 1, 100]],
  100
]


MT[
  Union[IGRandomWalk[IGShorthand["1-2-3", EdgeWeight -> {2, 1}], 1, 100]],
  {1, 2, 3}
]

MT[
  Union[Rest[IGRandomWalk[IGShorthand["1-2-3", EdgeWeight -> {0, 1}], 1, 100]]],
  {2, 3}
]

MT[
  Union[Rest[IGRandomWalk[IGShorthand["1-2-3", EdgeWeight -> {0, 1}], 1, 100, EdgeWeight -> None]]],
  {1, 2, 3}
]


(* ::Section::Closed:: *)
(*Modification*)


MTSection["Modification"]


(* ::Subsubsection::Closed:: *)
(*IGConnectNeighborhood*)


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
    IGSquareLattice[{10}, "Periodic" -> True, "Radius" -> 2]
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


(* ::Subsubsection:: *)
(*IGVertexContract*)


(* TODO *)


(* ::Subsubsection::Closed:: *)
(*IGMycielskian*)


MT[
  IGMycielskian[CycleGraph[4]],
  Graph[{1, 2, 3, 4, 5, 6, 7, 8, 9}, {1 <-> 2, 1 <-> 4, 2 <-> 3, 3 <-> 4, 5 <-> 9, 6 <-> 9, 7 <-> 9, 8 <-> 9, 1 <-> 6, 2 <-> 5, 1 <-> 8, 4 <-> 5, 2 <-> 7, 3 <-> 6, 3 <-> 8, 4 <-> 7}],
  SameTest -> IGSameGraphQ
]


MT[
  IGMycielskian[IGEmptyGraph[]],
  Graph[{1}, {}],
  SameTest -> IGSameGraphQ
]


MT[
  IGMycielskian[IGEmptyGraph[1]],
  Graph[{1 <-> 2}],
  SameTest -> IGSameGraphQ
]


MT[
  IGMycielskian[IGMycielskian[CycleGraph[3]]],
  Graph[{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15}, {1 <-> 2, 1 <-> 3, 2 <-> 3, 4 <-> 7, 5 <-> 7, 6 <-> 7, 1 <-> 5, 2 <-> 4, 1 <-> 6, 3 <-> 4, 2 <-> 6, 3 <-> 5, 8 <-> 15, 9 <-> 15, 10 <-> 15, 11 <-> 15, 12 <-> 15, 13 <-> 15, 14 <-> 15, 1 <-> 9, 2 <-> 8, 1 <-> 10, 3 <-> 8, 2 <-> 10, 3 <-> 9, 4 <-> 14, 7 <-> 11, 5 <-> 14, 7 <-> 12, 6 <-> 14, 7 <-> 13, 1 <-> 12, 5 <-> 8, 2 <-> 11, 4 <-> 9, 1 <-> 13, 6 <-> 8, 3 <-> 11, 4 <-> 10, 2 <-> 13, 6 <-> 9, 3 <-> 12, 5 <-> 10}],
  SameTest -> IGSameGraphQ
]


(* ::Subsubsection::Closed:: *)
(*IGSmoothen*)


MT[
  IGSmoothen[IGEmptyGraph[]],
  Graph[{}],
  SameTest -> IGSameGraphQ
]

MT[
  IGSmoothen[IGEmptyGraph[4]],
  Graph[{1, 2, 3, 4}, {}],
  SameTest -> IGSameGraphQ
]


MT[
  IGSmoothen[Graph[{1 <-> 2, 1 <-> 2, 2 <-> 3}]],
  Graph[{2 <-> 2, 2 <-> 3}],
  SameTest -> IGSameGraphQ
]

MT[
  IGEdgeProp[EdgeWeight][IGSmoothen[Graph[{1 <-> 2, 1 <-> 2, 2 <-> 3}]]],
  {2., 1.}
]

MT[
  IGSmoothen[CycleGraph[5]],
  Graph[{5 <-> 5}],
  SameTest -> IGIsomorphicQ
]

MT[
  IGEdgeProp[EdgeWeight][IGSmoothen[CycleGraph[5]]],
  {5.}
]


MT[
  IGSmoothen[Graph[{1 -> 2, 2 -> 3}]],
  Graph[{1 -> 3}],
  SameTest -> IGSameGraphQ
]

MT[
  IGEdgeProp[EdgeWeight][IGSmoothen[Graph[{1 -> 2, 2 -> 3}]]],
  {2.}
]


MT[
  IGSmoothen[Graph[{1 -> 1, 1 -> 2, 2 -> 3, 3 -> 4, 5 -> 4, 5 -> 2}]],
  Graph[{1 -> 1, 1 -> 2, 2 -> 4, 5 -> 2, 5 -> 4}],
  SameTest -> IGSameGraphQ
]


MT[
  IGSmoothen[Graph[{1 -> 2, 2 <-> 3}]],
  Quiet@IGSmoothen[Graph[{1 -> 2, 2 <-> 3}]],
  {IGraphM::mixed}
]


(* ::Section::Closed:: *)
(*Rewiring*)


MTSection["Rewiring"]


(* ::Subsubsection::Closed:: *)
(*IGRewire*)


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


(* ::Subsubsection::Closed:: *)
(*IGRewireEdges*)


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


(* TODO In, Out for directed *)
(* TODO ensure that In,Out version respects SelfLoops and MultiEdges setting *)


(* ::Section::Closed:: *)
(*Acyclic graphs*)


MTSection["Acyclic graphs"]


(* ::Subsubsection::Closed:: *)
(*IGTopologicalOrdering*)


(* check that reordering the adjacency matrix makes it upper triangular *)
MT[
  Module[{g = DirectedGraph[RandomGraph[{10,20}],"Acyclic"],ord,am},
    g = IGReorderVertices[RandomSample@VertexList[g],g];
    ord = IGTopologicalOrdering[g];
    am = AdjacencyMatrix[g];
    Union@Flatten@LowerTriangularize@Normal[am[[ord,ord]]]
  ],
  {0}
]


(* check that it does not error on edgeless graph; any permutation is a valid answer *)
MT[
  Sort@IGTopologicalOrdering@IGEmptyGraph[5],
  Range[5]
]


(* behaviour with non-acylic graph *)
MT[
  IGTopologicalOrdering[CycleGraph[3, DirectedEdges -> True]],
  $Failed,
  {IGraphM::error}
]


(* ::Subsubsection::Closed:: *)
(*IGDirectedAcyclicGraph*)


MT[
  IGDirectedAcyclicGraphQ[IGEmptyGraph[#]],
  True
]& /@ Range[0,2]


MT[
  IGDirectedAcyclicGraphQ[Graph[{1 -> 2}]],
  True
]


MT[
  IGDirectedAcyclicGraphQ[CycleGraph[2, DirectedEdges -> True]],
  False
]

MT[
  IGDirectedAcyclicGraphQ[CycleGraph[3, DirectedEdges -> True]],
  False
]


(* multi-edge *)
MT[
  IGDirectedAcyclicGraphQ[Graph[{1 -> 2, 1 -> 2}]],
  True
]


(* self-loop *)
MT[
  IGDirectedAcyclicGraphQ[Graph[{1 -> 2, 1 -> 1}]],
  False
]


(* ::Subsubsection::Closed:: *)
(*IGFeedbackArcSet*)


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


(* ::Section::Closed:: *)
(*Trees*)


MTSection["Trees"]


(* ::Subsubsection::Closed:: *)
(*IGTreeQ*)


MT[
  IGTreeQ[IGEmptyGraph[]],
  False
]


MT[
  IGTreeQ@IGEmptyGraph[1],
  True
]


MT[
  IGTreeQ[Graph[{1->2}]],
  True
]


MT[
  IGTreeQ[Graph[{1<->2}]],
  True
]


MT[
  IGTreeQ[Graph@{1->2,3->2}],
  False
]

MT[
  IGTreeQ[Graph@{1->2,3->2}, "Out"],
  False
]

MT[
  IGTreeQ[Graph@{1->2,3->2}, "In"],
  True
]

MT[
  IGTreeQ[Graph@{1->2,3->2}, "All"],
  True
]


MT[
  IGTreeQ[Graph@{1->2,1->3,4->3}],
  False
]

MT[
  IGTreeQ[Graph@{1->2,1->3,4->3}, "Out"],
  False
]

MT[
  IGTreeQ[Graph@{1->2,1->3,4->3}, "In"],
  False
]

MT[
  IGTreeQ[Graph@{1->2,1->3,4->3}, "All"],
  True
]


MT[
  AllTrue[Table[IGTreeGame@RandomInteger[{1,50}], {10}], IGTreeQ],
  True
]


MT[
  AllTrue[Table[IGTreeGame[RandomInteger[{1,50}], DirectedEdges -> True], {10}], IGTreeQ],
  True
]


(* https://github.com/szhorvat/IGraphM/issues/90 *)
MT[
  IGTreeQ@Graph[Range[5], {1<->4, 1<->5, 2<->4, 2<->5}],
  False
]


(* ::Subsubsection::Closed:: *)
(*IGForestQ*)


(* ::Text:: *)
(*Undirected forests*)


(* K_0 is a forest but not a tree *)
MT[
  IGForestQ@IGEmptyGraph[#],
  True
]& /@ Range[0,3]


MT[
  IGForestQ@IGShorthand["1-2-3,4-5-6"],
  True
]


MT[
  IGForestQ@IGShorthand["1-2-3-1,4-5-6"],
  False
]


MT[
  IGForestQ[IGShorthand["1-2-3-2,4-5", MultiEdges -> True]],
  False
]


MT[
  IGForestQ[IGShorthand["1-2-3-3,4-5", MultiEdges -> True, SelfLoops -> True]],
  False
]


MT[
  IGForestQ[IGShorthand["1,2,3"]],
  True
]


MT[
  IGForestQ[IGShorthand["1-2-3-1,4"]],
  False
]


(* ::Text:: *)
(*Directed forests*)


MT[
  IGForestQ[IGShorthand["1->2->3"], "Out"],
  True
]


MT[
  IGForestQ[IGShorthand["1->2->3"], "In"],
  True
]


MT[
  IGForestQ[IGShorthand["1->2->3,2->4"], "In"],
  False
]


MT[
  IGForestQ[IGShorthand["1->2->3,2->4,5<-6<-7,6<-8"], "Out"],
  False
]


MT[
  IGForestQ[IGShorthand["1->2->3,2->4,5<-6<-7,6<-8"], "In"],
  False
]


MT[
  IGForestQ[IGShorthand["1->2->3,2->4,5<-6<-7,6<-8"], "All"],
  True
]


(* ::Subsubsection::Closed:: *)
(*IGTreelikeComponents*)


MT[
  Sort[IGTreelikeComponents[Graph[{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, {6 <-> 9, 8 <-> 9, 5 <-> 8, 3 <-> 8, 3 <-> 7, 3 <-> 4, 2 <-> 7, 1 <-> 5, 4 <-> 10}]]],
  {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}
]


MT[
  Sort[IGTreelikeComponents[CycleGraph[3]]],
  {}
]


MT[
  Sort[IGTreelikeComponents[CycleGraph[10]]],
  {}
]


MT[
  Sort[IGTreelikeComponents[PathGraph[Range[10]]]],
  {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}
]


MT[
  Sort[IGTreelikeComponents[IGEmptyGraph[]]],
  {}
]


MT[
  Sort[IGTreelikeComponents[IGShorthand["1-2-3-4-2"]]],
  {1}
]


(* directed -- directions are ignored *)
MT[
  Sort[IGTreelikeComponents[IGShorthand["1->2->3->4->2"]]],
  {1}
]


(* multigraphs *)
MT[
  IGTreelikeComponents[Graph[{1 <-> 2, 1 <-> 2, 2 <-> 3}]],
  {3}
]

MT[
  IGTreelikeComponents[Graph[{1 <-> 2, 1 <-> 2}]],
  {}
]


(* ::Subsubsection::Closed:: *)
(*IGStrahlerNumber*)


MT[
  IGStrahlerNumber[IGEmptyGraph[1]],
  {1}
]


MT[
  IGStrahlerNumber[IGEmptyGraph[]],
  $Failed,
  {IGraphM::error}
]


MT[
  IGStrahlerNumber[Graph[{1->2, 3->2}]],
  $Failed,
  {IGraphM::error}
]


MT[
  IGStrahlerNumber[Graph[{1<->2}]],
  $Failed,
  {IGraphM::error}
]


MT[
  IGStrahlerNumber@Graph[
    {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20},
    {17->9,9->8,17->15,9->14,14->3,3->19,14->11,11->7,7->4,4->16,9->1,7->5,5->6,6->10,16->18,19->12,3->2,2->13,8->20}
  ],
  {1,2,3,3,3,2,4,2,7,1,5,1,1,6,1,2,8,1,2,1}
]


(* ::Subsubsection::Closed:: *)
(*IGFromPrufer*)


MT[
  IGFromPrufer[{}],
  Graph[{1 <-> 2}],
  SameTest -> IGSameGraphQ
]


MT[
  IGFromPrufer[{0}],
  $Failed,
  {IGraphM::error}
]


MT[
  IGFromPrufer[{1}],
  Graph[{1 <-> 2, 1 <-> 3}],
  SameTest -> IGSameGraphQ
]


MT[
  IGFromPrufer[{2}],
  Graph[{1 <-> 2, 2 <-> 3}],
  SameTest -> IGSameGraphQ
]


MT[
  IGFromPrufer[{3}],
  Graph[{1, 2, 3}, {1 <-> 3, 2 <-> 3}],
  SameTest -> IGSameGraphQ
]


MT[
  IGFromPrufer[{4}],
  $Failed,
  {IGraphM::error}
]


MT[
  IGFromPrufer[{1, 5, 2}],
  Graph[{1, 2, 3, 4, 5}, {1 <-> 3, 1 <-> 5, 2 <-> 4, 2 <-> 5}],
  SameTest -> IGSameGraphQ
]


MT[
  IGFromPrufer@IGToPrufer[#],
  #,
  SameTest -> IGSameGraphQ
]& /@ Table[IGTreeGame[k], {k, 3, 10}]


MT[
  IGFromPrufer /@ Tuples[Range[4], {2}],
  Graph /@ {{1 <-> 2,1 <-> 3,1 <-> 4},{1 <-> 3,1 <-> 2,2 <-> 4},{1 <-> 2,1 <-> 3,3 <-> 4},{1 <-> 2,1 <-> 4,3 <-> 4},{2 <-> 3,1 <-> 2,1 <-> 4},{1 <-> 2,2 <-> 3,2 <-> 4},{1 <-> 2,2 <-> 3,3 <-> 4},{1 <-> 2,2 <-> 4,3 <-> 4},{2 <-> 3,1 <-> 3,1 <-> 4},{1 <-> 3,2 <-> 3,2 <-> 4},{1 <-> 3,2 <-> 3,3 <-> 4},{1 <-> 3,2 <-> 4,3 <-> 4},{2 <-> 4,1 <-> 3,1 <-> 4},{1 <-> 4,2 <-> 3,2 <-> 4},{1 <-> 4,2 <-> 3,3 <-> 4},{1 <-> 4,2 <-> 4,3 <-> 4}},
  SameTest -> (And @@ MapThread[IGSameGraphQ, {#1, #2}]&)
]


MT[
  IGFromPrufer@{5,18,14,19,16,14,9,20,8,17,17,10,14,4,4,9,16,3},
  Graph@{4 <-> 19,4 <-> 14,5 <-> 14,10 <-> 14,6 <-> 19,4 <-> 9,9 <-> 16,3 <-> 16,3 <-> 20,13 <-> 20,1 <-> 5,7 <-> 16,10 <-> 17,8 <-> 17,11 <-> 14,8 <-> 15,17 <-> 18,2 <-> 18,9 <-> 12},
  SameTest -> IGSameGraphQ
]


(* ::Subsubsection::Closed:: *)
(*IGToPrufer*)


MT[
  IGToPrufer@IGEmptyGraph[],
  $Failed,
  {IGraphM::error}
]

MT[
  IGToPrufer@IGEmptyGraph[1],
  $Failed,
  {IGraphM::error}
]


MT[
  IGToPrufer[Graph[{1 <-> 2}]],
  {}
]


MT[
  IGToPrufer[Graph[{1 <-> 2, 2 <-> 3}]],
  {2}
]


MT[
  IGToPrufer[Graph[{1 <-> 2, 1 <-> 3}]],
  {1}
]


MT[
  IGToPrufer[Graph[{1 <-> 2, 1 <-> 3, 2 <-> 3}]],
  $Failed,
  {IGraphM::error}
]


MT[
  IGToPrufer[Graph[{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, {6 <-> 9, 8 <-> 9, 5 <-> 8, 3 <-> 8, 3 <-> 7, 3 <-> 4, 2 <-> 7, 1 <-> 5, 4 <-> 10}]],
  {5, 7, 8, 9, 3, 8, 3, 4}
]


(* directed -- directions are ignored *)
MT[
  IGToPrufer[IGShorthand["1->2<-3"]],
  {2}
]


(* multigraph -- not a tree *)
MT[
  IGToPrufer[Graph[{1 <-> 2, 1 <-> 2}]],
  $Failed,
  {IGraphM::error}
]


(* ::Subsubsection::Closed:: *)
(*IGSpanningTreeCount*)


(* ::Text:: *)
(*Undirected graphs*)


MT[
  IGSpanningTreeCount[IGEmptyGraph[]],
  0
]


MT[
  IGSpanningTreeCount[IGEmptyGraph[1]],
  1
]


(* not connected, so zero *)
MT[
  IGSpanningTreeCount[IGEmptyGraph[2]],
  0
]


MT[
  IGSpanningTreeCount[CycleGraph[3]],
  3
]


(* multigraph *)
MT[
  IGSpanningTreeCount[CycleGraph[2]],
  2
]


MT[
  IGSpanningTreeCount[CompleteGraph[2]],
  1
]


MT[
  IGSpanningTreeCount[CompleteGraph[4]],
  16
]


MT[
  IGSpanningTreeCount[GridGraph[{2, 3}]],
  15
]


MT[
  IGSpanningTreeCount[Graph[{1 <-> 2, 3 <-> 4}]],
  0
]


(* ::Text:: *)
(*Directed graphs*)


MT[
  IGSpanningTreeCount[GridGraph[{2, 3}, DirectedEdges -> True]],
  4
]


MT[
  IGSpanningTreeCount[Graph[{1 -> 2}]],
  1
]


MT[
  IGSpanningTreeCount[Graph[{1 -> 2, 2 -> 1}]],
  2
]


(* ::Text:: *)
(*Directed graph with starting vertex*)


MT[
  IGSpanningTreeCount[GridGraph[{2, 3}, DirectedEdges -> True], 1],
  4
]


MT[
  IGSpanningTreeCount[GridGraph[{2, 3}, DirectedEdges -> True], 2],
  0
]


MT[
  IGSpanningTreeCount[Graph[{1 -> 2, 2 -> 1}], 2],
  1
]


(* ::Text:: *)
(*Nonexistent vertex*)


MT[
  IGSpanningTreeCount[Graph[{1 -> 2, 2 -> 1}], 3],
  $Failed,
  {VertexIndex::inv}
]


MT[
  IGSpanningTreeCount[Graph[{1 <-> 2}], 3],
  $Failed,
  {VertexIndex::inv}
]


With[{g=IGGiantComponent@RandomGraph[{9,13}]},
  MT[
    IGSpanningTreeCount[g],
    Length@DeleteDuplicatesBy[IGRandomSpanningTree[g,5000],AdjacencyMatrix]
  ]
]


(* ::Subsubsection::Closed:: *)
(*IGRandomSpanningTree*)


MT[
  IGRandomSpanningTree[IGEmptyGraph[]],
  Graph[{}],
  SameTest -> IGSameGraphQ
]


MT[
  IGRandomSpanningTree[IGEmptyGraph[1]],
  Graph[{1}, {}],
  SameTest -> IGSameGraphQ
]


MT[
  IGRandomSpanningTree[IGEmptyGraph[5]],
  Graph[{1, 2, 3, 4, 5}, {}],
  SameTest -> IGSameGraphQ
]


MT[
  IGTreeQ@IGRandomSpanningTree@IGGiantComponent@RandomGraph[{10, 20}],
  True
]


MT[
  IGForestQ@IGRandomSpanningTree@RandomGraph[{20, 20}],
  True
]


MT[
  AllTrue[IGRandomSpanningTree[IGGiantComponent@RandomGraph[{10, 20}], 10], IGTreeQ],
  True
]


MT[
  IGRandomSpanningTree[{Graph[{1, 2, 3}, {1 <-> 2}], 1}],
  Graph[{1 <-> 2}],
  SameTest -> IGSameGraphQ
]


MT[
  IGRandomSpanningTree[{Graph[{1 <-> 2, 3 <-> 4}], 3}],
  Graph[{3 <-> 4}],
  SameTest -> IGSameGraphQ
]


(* https://github.com/szhorvat/IGraphM/issues/89 *)
MT[
  IGRandomSpanningTree[{Graph[{1, 2, 3}, {1 <-> 2}], 3}],
  Graph[{3},{}],
  SameTest -> IGSameGraphQ
]


(* syntax for generating multiple trees *)
MT[
  With[{trees = IGRandomSpanningTree[IGTriangularLattice[5], 6]}, AllTrue[trees, IGTreeQ] && Length[trees] == 6],
  True
]


MT[
  IGRandomSpanningTree[CompleteGraph[5], 0],
  {}
]


MT[
  VertexCount /@ IGRandomSpanningTree[{IGDisjointUnion[{CompleteGraph[5], CycleGraph[10]}], {1, 1}}, 5],
  {5, 5, 5, 5, 5}
]


(* ::Subsubsection::Closed:: *)
(*IGSpanningTree*)


MT[
  IGSpanningTree[IGEmptyGraph[]],
  Graph[{}],
  SameTest -> IGSameGraphQ
]

MT[
  IGSpanningTree[IGEmptyGraph[1]],
  Graph[{1}, {}],
  SameTest -> IGSameGraphQ
]

MT[
  IGSpanningTree[IGEmptyGraph[2]],
  Graph[{1, 2}, {}],
  SameTest -> IGSameGraphQ
]

MT[
  IGSpanningTree[IGShorthand["1-2"]],
  Graph[{1 <-> 2}],
  SameTest -> IGSameGraphQ
]

MT[
  IGSpanningTree[IGShorthand["1->2"]],
  Graph[{1 -> 2}],
  SameTest -> IGSameGraphQ
]

MT[
  IGSpanningTree[IGShorthand["1->2<-3"]],
  Graph[{1 -> 2, 3 -> 2}],
  SameTest -> IGSameGraphQ
]

MT[
  IGSpanningTree[IGShorthand["1-2,2-3,4"]],
  Graph[{1, 2, 3, 4}, {1 <-> 2, 2 <-> 3}],
  SameTest -> IGSameGraphQ
]


MT[
  With[{t = IGTreeGame[500]}, IGSameGraphQ[t, IGSpanningTree[t]]],
  True
]


MT[
  IGEdgeProp[EdgeWeight][IGSpanningTree[Graph[{1 <-> 2}, EdgeWeight -> {34}]]],
  {34}
]


(* ::Subsubsection::Closed:: *)
(*IGUnfoldTree*)


MT[
  IGUnfoldTree[IGEmptyGraph[], {}],
  Graph[{}],
  SameTest -> IGSameGraphQ
]

MT[
  IGUnfoldTree[IGEmptyGraph[1], {1}],
  Graph[{1}, {}],
  SameTest -> IGSameGraphQ
]

MT[
  IGUnfoldTree[IGShorthand["1-2-3-1"], {2}],
  Graph[{1 <-> 2, 2 <-> 3, 1 <-> 4}],
  SameTest -> IGSameGraphQ
]

MT[
  IGUnfoldTree[IGShorthand["1-2-3-1"], {}],
  Graph[{1, 2, 3}, {}],
  SameTest -> IGSameGraphQ
]


(* ::Section::Closed:: *)
(*Isomorphism*)


MTSection["Isomorphism"]


(* ::Subsection::Closed:: *)
(*General isomorphism tests*)


grid = IGSquareLattice[{3,3,3,3}]; (* do not use GridGraph because it may crash in M11.3 *)


randomize[g_] := Graph[RandomSample@VertexList[g], EdgeList[g]]


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


(* ::Subsection::Closed:: *)
(*IGGetIsomorphism, IGGetSubisomorphism*)


MT[
  IGGetIsomorphism[IGEmptyGraph[], IGEmptyGraph[]],
  {<||>}
]


MT[
  IGGetIsomorphism[IGEmptyGraph[], IGEmptyGraph[1]],
  {}
]


MT[
  IGGetSubisomorphism[IGEmptyGraph[], IGEmptyGraph[1]],
  {<||>}
]


MT[
  IGGetSubisomorphism[IGEmptyGraph[], IGEmptyGraph[]],
  {<||>}
]


(* TODO *)


(* ::Subsection::Closed:: *)
(*Directed graphs*)


directedIsomorphismTest[if_] := {
  MT[
    if[dgs, randomize@dgs],
    True
  ],
  MT[
    if[dgs, ugs],
    $Failed,
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


(* ::Subsection::Closed:: *)
(*Automorphism groups*)


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
  GroupOrder@IGBlissAutomorphismGroup[ GraphData["PerkelGraph"] ],
  GraphData["PerkelGraph", "AutomorphismCount"]
]


MT[
  IGBlissAutomorphismCount[#],
  1
]& /@ asymmList


MT[
  IGBlissAutomorphismCount[IGEmptyGraph[#]],
  #!
]& /@ Range[0,5]

MT[
  IGBlissAutomorphismCount[IGCompleteGraph[#]],
  #!
]& /@ Range[2,5]


(* ::Subsection::Closed:: *)
(*Self-loop handling*)


(* ::Text:: *)
(*VF2 does not support self loops.*)


MT[
  #[
      Graph[{1, 2, 3}, {1 <-> 2, 2 <-> 3, 3 <-> 3}],
      Graph[{1, 2, 3}, {1 <-> 2, 2 <-> 3}]
      ],
  False
]& /@ {IGIsomorphicQ, IGBlissIsomorphicQ, IGSubisomorphicQ}

MT[
  #[
      Graph[{1, 2, 3}, {1 <-> 2, 2 <-> 3, 3 <-> 3}],
      Graph[{1, 2, 3}, {1 <-> 2, 2 <-> 3, 1 <-> 1}]
      ],
  True
]& /@ {IGIsomorphicQ, IGBlissIsomorphicQ, IGSubisomorphicQ}

MT[
  #[
      Graph[{1, 2, 3}, {1 <-> 2, 2 <-> 3, 3 <-> 3}],
      Graph[{1, 2, 3}, {1 <-> 2, 2 <-> 3, 2 <-> 2}]
      ],
  False
]& /@ {IGIsomorphicQ, IGBlissIsomorphicQ, IGSubisomorphicQ}



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
  {{1, 2}, {DirectedEdge[1, 2], DirectedEdge[2, 1]}}
]

MT[
  With[{g = IGBlissCanonicalGraph[Graph[{2 -> 1, 1 -> 2}]]},
    {VertexList[g], EdgeList[g]}
  ],
  {{1, 2}, {DirectedEdge[1, 2], DirectedEdge[2, 1]}}
]

MT[
  With[{g = IGBlissCanonicalGraph[Graph[{b <-> a}]]},
    {VertexList[g], EdgeList[g]}
  ],
  {{1, 2}, {UndirectedEdge[1, 2]}}
]


(* ::Subsection::Closed:: *)
(*LAD*)


MT[
  IGLADSubisomorphicQ[CycleGraph[15], GraphData[{"Arrangement", {5, 2}}]],
  True
]


MT[
  IGLADSubisomorphicQ[CycleGraph[15], GraphData[{"Arrangement", {5, 2}}], "Induced" -> True],
  False
]


MT[
  (IGLADSubisomorphismCount[#,dgs,"Induced"->True]&/@IGData[{"AllDirectedGraphs",3}])/IGBlissAutomorphismCount/@IGData[{"AllDirectedGraphs",3}],
  Lookup[IGTriadCensus[dgs], Keys@IGData["MANTriadLabels"]]
]


Module[{bigDirected=RandomGraph[{20,200},DirectedEdges->True], pos, res1, res2},
  MT[
    res1=IGMotifs[bigDirected,4];
    res2=(IGLADSubisomorphismCount[#,bigDirected,"Induced"->True]&/@IGData[{"AllDirectedGraphs",4}])/IGBlissAutomorphismCount/@IGData[{"AllDirectedGraphs",4}];
    pos=Position[res1,Indeterminate];
    Delete[res1,pos]==Delete[res2,pos]
    ,
    True
  ]
]


(* ::Section::Closed:: *)
(*Isomorphism: coloured graphs*)


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
  PermutationGroup[{}]
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
  {{1, 2}, {UndirectedEdge[1, 2]}, {1, 3}}
]

MT[
  With[{g = IGBlissCanonicalGraph[{Graph[{b <-> a}], "VertexColors" -> <|b -> 1, a -> 3|>}]},
    {VertexList[g], EdgeList[g], IGVertexProp["Color"][g]}
  ],
  {{1, 2}, {UndirectedEdge[1, 2]}, {1, 3}}
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


(* TODO non-trivial examples with isomorphism checking, and not with canonical labelling *)
(* TODO test for Bliss coloured graph bug that was fixed recently in the igraph core *)


(* ::Section::Closed:: *)
(*Isomorphism: multigraphs*)


MTSection["Isomorphism: multigraphs"]


(* ::Text:: *)
(*This section deals both with multigraphs and graphs that have self-loops. Several of the algorithms support neither.*)


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
  $Failed,
  {IGraphM::error}
]& /@ {IGLADFindSubisomorphisms, IGLADGetSubisomorphism, IGLADSubisomorphicQ}


MT[
  IGIsomorphicQ[
    IGShorthand["1-1-2-3-2", MultiEdges -> True, SelfLoops -> True], 
    IGShorthand["1-2-3,2-1,3-3", MultiEdges -> True, SelfLoops -> True]
  ],
  True
]


MT[
  IGIsomorphicQ[
    IGShorthand["1-1-2-3-2", MultiEdges -> True, SelfLoops -> True], 
    IGShorthand["1-2-3,2-1-1", MultiEdges -> True, SelfLoops -> True]
  ],
  False
]


MT[
  IGIsomorphicQ[
    IGShorthand["1-1-2-3-2", MultiEdges -> True, SelfLoops -> True], 
    IGShorthand["1-2-3,2-2-1", MultiEdges -> True, SelfLoops -> True]
  ],
  False
]


MT[
  IGSubisomorphicQ[Graph[{1 <-> 1}], Graph[{1 <-> 2, 2 <-> 3}]],
  False
]


MT[
  IGVF2SubisomorphicQ[Graph[{1 <-> 1}], Graph[{1 <-> 2, 2 <-> 3}]],
  $Failed,
  {IGraphM::vf2nmg}
]


MT[
  IGLADSubisomorphicQ[Graph[{1 <-> 1}], Graph[{1 <-> 2, 2 <-> 3, 2 <-> 2}]],
  True
]

MT[
  IGLADSubisomorphicQ[Graph[{1 <-> 1}], Graph[{1 <-> 2, 2 <-> 3, 3 <-> 3}], "Induced" -> True],
  True
]

MT[
  IGLADSubisomorphicQ[Graph[{1 <-> 1}], Graph[{1 <-> 2, 2 <-> 3}]],
  False
]

MT[
  IGLADSubisomorphicQ[Graph[{1 <-> 1}], Graph[{1 <-> 2, 2 <-> 3}], "Induced" -> True],
  False
]


MT[
  IGLADSubisomorphicQ[Graph[{1 <-> 1}], Graph[{1 <-> 2, 2 <-> 3, 2 <-> 3}]],
  $Failed,
  {IGraphM::error}
]


MT[
  IGSubisomorphicQ[Graph[{1 <-> 1}], Graph[{1 <-> 2, 2 <-> 3, 2 <-> 2}]],
  True
]


MT[
  IGSubisomorphicQ[Graph[{1 <-> 1}], Graph[{1 <-> 2, 2 <-> 3, 2 <-> 2, 2 <-> 2}]],
  True
]


MT[
  IGSubisomorphicQ[Graph[{1 <-> 2, 2 <-> 3, 1 <-> 2, 2 <-> 3}], Graph[{1, 2, 3}, {1 <-> 3, 3 <-> 2, 1 <-> 3}]],
  False
]


MT[
  IGSubisomorphicQ[Graph[{1 <-> 2, 2 <-> 3, 2 <-> 3}], Graph[{1, 2, 3}, {1 <-> 3, 3 <-> 2, 1 <-> 3}]],
  True
]


MT[
  IGSubisomorphicQ[Graph[{1 <-> 2, 2 <-> 3}], Graph[{1, 2, 3}, {1 <-> 3, 3 <-> 2, 1 <-> 3}]],
  True
]


(* ::Section::Closed:: *)
(*Centralities*)


MTSection["Centralities"]


MT[
  #@IGEmptyGraph[],
  {}
]& /@ {IGBetweenness, IGEdgeBetweenness, IGCloseness, IGPageRank, IGEigenvectorCentrality, IGHubScore, IGAuthorityScore, IGConstraintScore}


(* ::Subsubsection::Closed:: *)
(*Betweenness*)


(* ::Text:: *)
(*Note: BetweennessCentrality does not support multigraphs.*)


MT[
  BetweennessCentrality[#],
  IGBetweenness[#],
  SameTest -> tolEq
]& /@ {ugs, ugi, dgs, dgi}


MT[
  IGBetweenness[Graph[{1 <-> 2}], {}],
  {}
]

MT[
  IGBetweenness[IGEmptyGraph[]],
  {}
]

MT[
  IGBetweenness[IGEmptyGraph[3]],
  {0., 0., 0.}
]


MT[
  EdgeBetweennessCentrality[#],
  2 IGEdgeBetweenness[#],
  SameTest -> Equal
]& /@ {ugs, ugi}

MT[
  EdgeBetweennessCentrality[#],
  IGEdgeBetweenness[#],
  SameTest -> Equal
]& /@ {dgs, dgi}


(* test on one instance of a multigraph *)
MT[
  IGBetweenness@Graph[Range[6],{1 <-> 2,2 <-> 3,3 <-> 4,2 <-> 5,4 <-> 5,4 <-> 5,4 <-> 6}],
  {0., 4.333333333333333, 1.3333333333333333, 4.666666666666666, 2.6666666666666665, 0.},
  SameTest -> Equal
]


MT[
  IGBetweenness[umultiloopy],
  {4.666666666666667, 0.9999999999999999, 14.833333333333332, 1.9999999999999998, 10., 0., 8.75, 0., 11.25, 1.5}
]

MT[
  IGBetweenness[dmultiloopy],
  {30., 22., 31., 27., 14., 0., 0., 6.333333333333334, 0., 0., 12.666666666666668},
  SameTest -> Equal
]


(* ::Subsubsection::Closed:: *)
(*Closeness*)


MT[
  ClosenessCentrality[#],
  IGCloseness[#, Normalized -> True],
  SameTest -> Equal
]& /@ {ugs, dgs, wugs, wdgs, umulti, dmulti, umulti2, dmulti2}


MT[
  IGCloseness[Graph[{1 <-> 2}], {}],
  {}
]


MT[
  IGCloseness[IGEmptyGraph[]],
  {}
]

MT[
  IGCloseness[IGEmptyGraph[1]],
  {Indeterminate}
]

MT[
  IGCloseness[IGEmptyGraph[2]],
  {Indeterminate, Indeterminate}
]


(* must be compared using SameQ because of Indeterminate *)
Block[{Internal`$SameQTolerance = Log10[2.]*7},
  MT[
    IGCloseness[dmultiloopy],
    {0.42857142857142855, 0.2903225806451613, 0.34615384615384615, 0.34615384615384615, 1., Indeterminate, 0.24390243902439024, 0.2903225806451613, Indeterminate, Indeterminate, 0.2903225806451613}
  ]
]


MT[
	IGCloseness[dmultiloopy],
	IGCloseness@SimpleGraph[dmultiloopy]
]

MT[
	IGCloseness[umultiloopy],
	IGCloseness@SimpleGraph[umultiloopy]
]


(* ::Subsubsection::Closed:: *)
(*Range limited betweenness and closeness*)


MT[
  IGBetweenness[#],
  IGBetweennessCutoff[#, Infinity]
]& /@ {ugs, dgs}

MT[
  IGCloseness[#],
  IGClosenessCutoff[#, Infinity]
]& /@ {ugs, dgs}

MT[
  IGEdgeBetweenness[#],
  IGEdgeBetweennessCutoff[#, Infinity]
]& /@ {ugs, dgs}


(* TODO add non-Infinity example *)


(* ::Subsubsection::Closed:: *)
(*PageRank*)


MT[
  PageRankCentrality[#],
  IGPageRank[#],
  SameTest -> tolEq
]& /@ {ugs, dgs}

(* check that edge multiplicity is taken to be equivalent with edge weights *)
MT[
  IGPageRank[Graph[{1 <-> 2, 2 <-> 3}, EdgeWeight -> {1, 2}]],
  IGPageRank[Graph[{1 <-> 2, 2 <-> 3, 2 <-> 3}]]
]

MT[
  IGPageRank[Graph[{1 <-> 2, 2 <-> 3}], Method -> "Arnoldi"],
  IGPageRank[Graph[{1 <-> 2, 2 <-> 3}], Method -> "PRPACK"],
  SameTest -> Equal
]


MT[
  IGPageRank[dmultiloopy],
  {0.05456256995497844, 0.05915471396249873, 0.09522964768947766, 0.1652785569095764, 0.04262996884255929, 0.03755861336978115, 0.019440876611693456, 0.036201378901068094, 0.39942317919814296, 0.03755861336978115, 0.052961881190442726},
  SameTest -> Equal
]

MT[
  IGPageRank[umultiloopy],
  {0.11169488006788235, 0.06108975438778022, 0.15779712263544857, 0.08482491640220521, 0.07999158291232755, 0.11770088070362338, 0.10255919001820273, 0.0440584371718241, 0.17497752735059527, 0.06530570835011067},
  SameTest -> Equal
]


MT[
  IGPageRank[dmultiloopy, Method -> "Arnoldi"],
  IGPageRank[dmultiloopy, Method -> "PRPACK"],
  SameTest -> Equal
]

MT[
  IGPageRank[umultiloopy, Method -> "Arnoldi"],
  IGPageRank[umultiloopy, Method -> "PRPACK"],
  SameTest -> Equal
]


(* ::Subsubsection::Closed:: *)
(*IGEigenvectorCentrality*)


eigCent[g_] :=
    Module[{am = N@Transpose@WeightedAdjacencyMatrix[g], ev},
      If[UndirectedGraphQ[g],
        am += DiagonalMatrix@Diagonal[am]
      ];
      ev = First@Eigenvectors[am, 1];
      ev = Sign@First@MaximalBy[ev, Abs] ev; (* ensure non-negative *)
      ev / Max[ev] (* normalize such that max centrality is 1 *)
    ]


MT[
  Normalize[IGEigenvectorCentrality[#], Total],
  EigenvectorCentrality[#],
  SameTest -> tolEq
]& /@ {ugs, dgs, dolphin}


MT[
  IGEigenvectorCentrality[umulti],
  IGEigenvectorCentrality@IGWeightedSimpleGraph[umulti],
  SameTest -> Equal
]


(* weighted *)
MT[
  IGEigenvectorCentrality[lesmiserables],
  eigCent[lesmiserables],
  SameTest -> tolEq
]


MT[
  IGEigenvectorCentrality[umultiloopy],
  eigCent[umultiloopy],
  SameTest -> tolEq
]

MT[
  IGEigenvectorCentrality[dmultiloopy],
  eigCent[dmultiloopy],
  SameTest -> tolEq
]


(* ::Subsubsection:: *)
(*IGHubScore and IGAuthorityScore*)


(* TODO unweighted and weighted *)


(* ::Subsubsection:: *)
(*IGConstraintScore*)


MT[
  IGConstraintScore[friendship],
  {0.3742283950617284, 0.611111111111111, 0.611111111111111, 0.6759259259259258, 0.4598765432098765, 0.7847222222222223, 0.5, 0.5, 0.5},
  SameTest -> Equal
]


(* TODO unweighted and weighted *)


(* ::Section::Closed:: *)
(*Clustering coefficients*)


MTSection["Clustering coefficients"]


(* ::Subsection::Closed:: *)
(*Global, local, average local*)


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
  GlobalClusteringCoefficient[umulti2] == IGGlobalClusteringCoefficient[umulti2],
  True
]

MT[
  LocalClusteringCoefficient[umulti2] == IGLocalClusteringCoefficient[umulti2],
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


(* ::Subsection::Closed:: *)
(*Weighted (Barrat's)*)


MT[
  IGWeightedClusteringCoefficient[IGEmptyGraph[]],
  {},
  {IGraphM::warning}
]


MT[
  IGWeightedClusteringCoefficient[IGEmptyGraph[1]],
  {0.},
  {IGraphM::warning}
]


MT[
  IGWeightedClusteringCoefficient[Graph[{1 <-> 2}, EdgeWeight -> {2}]],
  {0., 0.}
]


MT[
  IGWeightedClusteringCoefficient[IGCompleteGraph[3, EdgeWeight -> {1, 2, 3}]],
  {1., 1., 1.}
]


MT[
  IGWeightedClusteringCoefficient[IGShorthand["1:2:3-1:2:3, 3-4, 4:5:6:7-4:5:6:7", EdgeWeight -> Range[10]]],
  {1., 1., 0.2777777777777778, 0.5454545454545454, 1., 1., 1.},
  SameTest -> Equal
]


(* ::Section::Closed:: *)
(*Cliques*)


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
]& /@ {ugs, ugi, umulti, umulti2}

MT[
  IGCliqueNumber@GraphComplement[#] == IGIndependenceNumber[#],
  True
]& /@ {ugs, ugi, umulti, umulti2}

MT[
  IGCliqueNumber[#] == Length@First@FindClique[#],
  True
]& /@ {ugs, ugi, umulti, umulti2}

MT[
  IGIndependenceNumber[#] == Length@First@FindIndependentVertexSet[#],
  True
]& /@ {ugs, ugi, umulti, umulti2}

MT[
  IGCliqueNumber[empty],
  0
]

MT[
  IGIndependenceNumber[empty],
  0
]


MT[
  canon@IGCliques[IGCompleteGraph[3], {3}],
  {{1, 2, 3}}
]


MT[
  canon@IGCliques[IGCompleteGraph[3], {2, 3}],
  canon@{{2, 3}, {1, 2}, {1, 2, 3}, {1, 3}}
]


MT[
  canon@IGCliques[IGCompleteGraph[3], 3],
  canon@{{3}, {2}, {2, 3}, {1}, {1, 2}, {1, 2, 3}, {1, 3}}
]


MT[
  canon@IGCliques[IGCompleteGraph[4], Infinity],
  canon@{{4}, {3}, {3, 4}, {2}, {2, 3}, {2, 3, 4}, {2, 4}, {1}, {1, 2}, {1, 2, 3}, {1, 2, 3, 4}, {1, 2, 4}, {1, 3}, {1, 3, 4}, {1, 4}}
]


MT[
  IGCliqueSizeCounts[IGCompleteGraph[#]],
  Table[Binomial[#,k],{k,1,#}]
]& /@ Range[0,5]


MT[
  Table[IGCliqueSizeCounts[IGCompleteGraph[#],{j}][[j]],{j,1,#}],
  Table[Binomial[#,k],{k,1,#}]
]& /@ Range[0,5]


MT[
  IGCliqueSizeCounts[GraphData["IcosahedralGraph"]],
  {12, 30, 20}
]


(* ::Section::Closed:: *)
(*Weighted cliques*)


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


(* ::Section::Closed:: *)
(*Motifs and subgraph counts*)


MTSection["Motifs and subgraph counts"]


(* ::Subsubsection::Closed:: *)
(*IGTriangles*)


triangleCount[am_?MatrixQ] := Tr[am . am . am]/6
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
  Length@IGTriangles[dolphin],
  95
]


(* ::Subsubsection::Closed:: *)
(*IGMotifs*)


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
  IGMotifs[dolphin, 3],
  {Indeterminate, Indeterminate, 638, 95}
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


(* ::Subsubsection::Closed:: *)
(*IGMotifsTotalCount*)


MT[
  IGMotifsTotalCount[web, 3],
  102246
]

MT[
  IGMotifsTotalCount[web, 4],
  4052036
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


(* ::Subsubsection::Closed:: *)
(*IGTriadCensus*)


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


(* ::Subsubsection::Closed:: *)
(*IGDyadCensus*)


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


(* ::Subsubsection::Closed:: *)
(*IGAdjacentTriangleCount*)


MT[
  IGAdjacentTriangleCount[#],
  With[{am = AdjacencyMatrix[#]}, Normal@Diagonal[am . am . am]/2]
]& /@ ulist

MT[
  IGAdjacentTriangleCount[empty],
  {}
]

MT[
  IGAdjacentTriangleCount[Graph[{1,2,3},{}]],
  {0,0,0}
]


(* ::Section::Closed:: *)
(*Matching*)


MTSection["Matching"]


MT[
  IGMatchingNumber[IGEmptyGraph[#]],
  0
]& /@ Range[0, 3]


MT[
  IGMatchingNumber[Graph[{1 <-> 2}]],
  1
]


MT[
  IGMatchingNumber[CycleGraph[6]],
  3
]

MT[
  IGMatchingNumber[CycleGraph[7]],
  3
]


MT[
  IGMaximumMatching[IGEmptyGraph[#]],
  {}
]& /@ Range[0, 3]


MT[
  IGMaximumMatching[IGShorthand["1-2-3-4-5-1-6"]],
  {UndirectedEdge[4, 5], UndirectedEdge[2, 3], UndirectedEdge[1, 6]}
]


(* ::Section::Closed:: *)
(*Connectivity*)


MTSection["Connectivity"]


(* ::Subsubsection::Closed:: *)
(*IGVertexConnectivity and IGEdgeConnectivity*)


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
  IGVertexConnectivity[umulti2],
  VertexConnectivity[umulti2]
]

MT[
  IGVertexConnectivity[dmulti2],
  VertexConnectivity[dmulti2]
]

MT[
  IGEdgeConnectivity[umulti2],
  EdgeConnectivity[umulti2]
]

MT[
  IGEdgeConnectivity[dmulti2],
  EdgeConnectivity[dmulti2]
]

MT[
  IGEdgeConnectivity[Graph[{2 <-> 1, 1 <-> 2, 1 <-> 3, 2 <-> 3, 3 <-> 1, 4 <-> 2}], 2, 3],
  3
]

MT[
  IGVertexConnectivity[ugs, 1, 2],
  $Failed,
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


(* This is a special case where the ideal result is unclear. For now, it is 0. *)
MT[
  IGEdgeConnectivity[IGEmptyGraph[1]],
  0
]

MT[
  IGEdgeConnectivity[IGEmptyGraph[0]],
  0
]


(* ::Subsubsection::Closed:: *)
(*IGBiconnectedQ*)


MT[
  IGBiconnectedQ[IGEmptyGraph[#]],
  False
]& /@ Range[0,3]


(* Special case: vertex connectivity is 1 but we consider it biconnected *)
MT[
  IGBiconnectedQ[Graph[{1 <-> 2}]],
  True
]


MT[
  IGBiconnectedQ[Graph[{1 <-> 2, 2 <-> 3}]],
  False
]


MT[
  IGBiconnectedQ[IGCompleteGraph[3]],
  True
]


MT[
  IGBiconnectedQ[IGShorthand["1-2-3-4-1,3-5-6-7-3"]],
  False
]


MT[
  IGBiconnectedQ[IGShorthand["1-2-3-4-1,3-5-6-7-3,5-4"]],
  True
]


(* multigraph *)
MT[
  IGBiconnectedQ[IGShorthand["1-2-3-4-1,3-5-6-7-3,5-3"]],
  False
]


(* edge directions are ignored *)
MT[
  IGBiconnectedQ[Graph[{1 -> 2, 2 -> 3, 1 -> 3}]],
  True
]


(* ::Subsubsection::Closed:: *)
(*IGBiconnectedComponents*)


MT[
  Sort[Sort /@ IGBiconnectedComponents[dolphin]],
  DeleteCases[Sort[Sort /@ KVertexConnectedComponents[dolphin, 2]], {_}]
]


(* TODO see articulation points detailed tests *)


(* ::Subsubsection::Closed:: *)
(*IGArticulationPoints*)


MT[
  IGArticulationPoints[dolphin],
  {"Ripplefluke", "Scabs", "Patchback", "SN63", "Trigger", "Web", "Jet"}
]


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


MT[
  IGArticulationPoints[IGEmptyGraph[#]],
  {}
]& /@ Range[0, 3]


MT[
  IGArticulationPoints[Graph[{1 <-> 2}]],
  {}
]


MT[
  IGArticulationPoints[Graph[{1 <-> 2, 2 <-> 3}]],
  {2}
]


(* ::Subsubsection::Closed:: *)
(*IGMinimumSeparators, IGMinimalSeparators, IGVertexSeparatorQ*)


MT[
  IGMinimumSeparators[#] =!= {} & /@ ulist,
  ConnectedGraphQ /@ ulist
]


gsep = IGShorthand["1-2-3-4-5-1,2-4-6-7"];

MT[
  Sort[Sort /@ IGMinimumSeparators[gsep]],
  {{4}, {6}}
]

MT[
  Sort[Sort /@ IGMinimalSeparators[gsep]],
  {{4}, {6}, {1, 4}, {2, 4}, {2, 5}}
]


MT[
  (IGVertexSeparatorQ[gsep, #1] & ) /@ IGMinimalSeparators[gsep],
  {True, True, True, True, True}
]

MT[
  (IGVertexSeparatorQ[gsep, #1] & ) /@ IGMinimumSeparators[gsep],
  {True, True}
]


MT[
  IGVertexSeparatorQ[gsep, {7}],
  False
]

MT[
  IGVertexSeparatorQ[gsep, {1}],
  False
]

MT[
  IGVertexSeparatorQ[gsep, {3, 5, 2}],
  True
]

MT[
  IGVertexSeparatorQ[gsep, {4, 7}],
  True
]

MT[
  IGVertexSeparatorQ[gsep, {}],
  False
]


MT[
  IGVertexSeparatorQ[IGEmptyGraph[], {}],
  False
]

MT[
  IGVertexSeparatorQ[IGEmptyGraph[], {1}],
  $Failed,
  {VertexIndex::inv}
]


MT[
  IGMinimalSeparators[IGEmptyGraph[#]],
  {}
]& /@ Range[0, 1]

MT[
  IGMinimumSeparators[IGEmptyGraph[#]],
  {}
]& /@ Range[0, 2]


(* ::Subsubsection::Closed:: *)
(*IGGiantComponent*)


MT[
  IGGiantComponent[Graph[{1 <-> 2}]],
  Graph[{1 <-> 2}],
  SameTest -> IGSameGraphQ
]


MT[
  IGGiantComponent[Graph[{1 <-> 2, 3 <-> 4, 4 <-> 5}]],
  Graph[{5, 4, 3}, {4 <-> 5, 3 <-> 4}],
  SameTest -> IGSameGraphQ
]


MT[
  IGGiantComponent[IGEmptyGraph[]],
  Graph[{}],
  SameTest -> IGSameGraphQ
]


MT[
  IGGiantComponent[IGEmptyGraph[1]],
  Graph[{1}, {}],
  SameTest -> IGSameGraphQ
]


MT[
  IGGiantComponent[Graph[{1 -> 2, 2 -> 3}]],
  Graph[{1 -> 2, 2 -> 3}],
  SameTest -> IGSameGraphQ
]


MT[
  IGGiantComponent[Graph[{1 -> 1}]],
  Graph[{1 -> 1}],
  SameTest -> IGSameGraphQ
]


MT[
  IGGiantComponent[Graph[{1, 2, 3}, {2 <-> 3}]],
  Graph[{2 <-> 3}],
  SameTest -> IGSameGraphQ
]


(* ::Subsubsection::Closed:: *)
(*IGConnectedQ and IGWeaklyConnectedQ*)


MT[
  IGConnectedQ /@ ulist,
  ConnectedGraphQ /@ ulist
]

MT[
  IGConnectedQ /@ dlist,
  ConnectedGraphQ /@ dlist
]

MT[
  IGWeaklyConnectedQ /@ ulist,
  ConnectedGraphQ /@ ulist
]

MT[
  IGWeaklyConnectedQ /@ dlist,
  WeaklyConnectedGraphQ /@ dlist
]


(* ::Text:: *)
(*Test a second time to verify result caching.*)


MT[
  IGConnectedQ /@ ulist,
  ConnectedGraphQ /@ ulist
]

MT[
  IGConnectedQ /@ dlist,
  ConnectedGraphQ /@ dlist
]

MT[
  IGWeaklyConnectedQ /@ ulist,
  ConnectedGraphQ /@ ulist
]

MT[
  IGWeaklyConnectedQ /@ dlist,
  WeaklyConnectedGraphQ /@ dlist
]


(* Special case: null graph *)
MT[
  IGConnectedQ[IGEmptyGraph[]],
  False
]

MT[
  IGWeaklyConnectedQ[IGEmptyGraph[]],
  False
]


MT[
  IGConnectedQ[IGEmptyGraph[1]],
  True
]

MT[
  IGWeaklyConnectedQ[IGEmptyGraph[1]],
  True
]

MT[
  IGConnectedQ[IGEmptyGraph[2]],
  False
]

MT[
  IGWeaklyConnectedQ[IGEmptyGraph[2]],
  False
]

MT[
  IGConnectedQ[IGShorthand["1->2"]],
  False
]

MT[
  IGWeaklyConnectedQ[IGShorthand["1->2"]],
  True
]


(* ::Text:: *)
(*Non-graph:*)


MT[
  IGConnectedQ[123],
  False
]

MT[
  IGWeaklyConnectedQ[123],
  False
]


(* ::Section::Closed:: *)
(*Community detection*)


MTSection["Community detection"]


funs = {IGCommunitiesEdgeBetweenness, IGCommunitiesGreedy, IGCommunitiesInfoMAP,
  IGCommunitiesLabelPropagation, IGCommunitiesMultilevel, IGCommunitiesWalktrap, 
  IGCommunitiesLeadingEigenvector, IGCommunitiesFluid[#, 4]&,
  IGCommunitiesSpinGlass[#, Method -> "Original"]&, IGCommunitiesSpinGlass[#, Method -> "Negative"]&,
  IGCommunitiesLeiden};

igClusterQ = Head[#] === IGClusterData&

MT[
  igClusterQ[#[dolphin]],
  True
]& /@ funs

MT[
  #[IGEmptyGraph[]]["Communities"],
  {}
]& /@ funs


(* ::Text:: *)
(*IGCommunitiesFluid requires specifying the number of clusters.*)


IGCommunitiesFluid[dolphin, 4]


(* ::Text:: *)
(*IGCommunitiesOptimalModularity is too slow for dolphins, so we test on lesmiserables.*)


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


(* ::Section::Closed:: *)
(*Shortest paths*)


MTSection["Shortest paths"]


(* ::Subsubsection::Closed:: *)
(*IGDistanceMatrix*)


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


(* ::Subsubsection::Closed:: *)
(*IGDistanceCounts and IGDistanceHistogram*)


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
]& /@ {empty, edgeless, ugs, dgs, umulti, dmulti, umulti2, dmulti2}

With[{v = RandomSample[VertexList[#], Quotient[VertexCount[#],3]]},
  MT[
    IGDistanceCounts[#, v],
    distanceCounts[#, v]
  ]
]& /@ {ugs, dgs, umulti, dmulti, umulti2, dmulti2}


MT[
  IGDistanceHistogram[empty, 1],
  {}
]

MT[
  IGDistanceHistogram[edgeless, 1],
  {}
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
  With[{vals = takeNonDiag@GraphDistanceMatrix[wdgs]},
    BinCounts[vals, {0, Ceiling[Max[vals], 0.1], 0.1}]
  ]
]


MT[
  IGDistanceHistogram[wugs, 0.1],
  With[{vals = takeNonDiag@GraphDistanceMatrix[wugs]},
    BinCounts[vals, {0, Ceiling[Max[vals], 0.1], 0.1}]
  ]
]


(* undirected, unweighted, connected *)
MT[
  IGDistanceCounts[dolphin],
  hist = Values@KeySort@Counts@IGTakeUpper@GraphDistanceMatrix[dolphin]
]


(* ::Subsubsection::Closed:: *)
(*IGAveragePathLength*)


MT[
  IGAveragePathLength[empty],
  Indeterminate
]

MT[
  IGAveragePathLength[edgeless],
  Indeterminate
]


MT[
  IGAveragePathLength[dolphin],
  N@Mean@WeightedData[Range@Length[hist], hist] (* hist is defined in the IGDistanceHistogram tests *)
]


(* disconnected *)
MT[
  IGAveragePathLength[IGShorthand["1-2,3-4"]],
  1.
]


(* ::Subsubsection::Closed:: *)
(*IGGirth*)


MT[
  IGGirth[PetersenGraph[5, 3]],
  5
]


(* For graphs with no cycles, we return Infinity. *)

MT[
  IGGirth[IGEmptyGraph[]],
  Infinity
]

MT[
  IGGirth[IGKaryTree[12]],
  Infinity
]


MT[
  IGGirth@CycleGraph[11],
  11
]


(* self-loops must be ignored *)
MT[
  IGGirth[IGShorthand["1-2-3-1-1", SelfLoops -> True]],
  3
]


(* ::Subsubsection::Closed:: *)
(*IGDiameter*)


MT[
  IGDiameter[empty],
  Indeterminate
]

MT[
  IGDiameter[edgeless],
  Infinity
]

MT[
  IGDiameter[#],
  GraphDiameter[#],
  SameTest -> Equal
]& /@ {ugs, ugi, dgs, dgi, wugs, wugi, wdgs, wdgi, umulti, dmulti, umulti2, dmulti2}


(* ::Subsubsection::Closed:: *)
(*IGFindDiameter*)


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


(* ::Subsubsection::Closed:: *)
(*IGLocalEfficiency and IGAverageLocalEfficiency*)


(* TODO *)


MT[
  IGLocalEfficiency[IGEmptyGraph[]],
  {}
]

MT[
  IGLocalEfficiency[IGEmptyGraph[3]],
  {0., 0., 0.}
]


MT[
  IGLocalEfficiency[dolphin],
  {0.5777777777777777, 0.46130952380952367, 0.4166666666666667, 0.6111111111111112, 0., 0.75, 0.7444444444444444, 0.5166666666666666, 0.5722222222222223, 0.7619047619047619, 0.625, 0., 0., 0.744047619047619, 0.577020202020202, 0.5555555555555555, 0.8, 0.3564814814814814, 0.7301587301587302, 0.75, 0.4699074074074072, 0.7555555555555556, 0., 0.611111111111111, 0.7555555555555556, 0.8333333333333334, 0.8333333333333334, 0.6416666666666667, 0.6000000000000001, 0.4513888888888888, 0.5583333333333333, 0., 0.3333333333333333, 0.5629629629629629, 0.5583333333333335, 0., 0.38299319727891157, 0.5575757575757578, 0.49404761904761896, 0.3333333333333333, 0.5428571428571428, 0.8, 0.6222222222222221, 0.5150793650793649, 0.4583333333333333, 0.6242424242424242, 0.3333333333333333, 0.6555555555555556, 0., 0.3333333333333333, 0.5793650793650792, 0.4185185185185186, 0.638888888888889, 0.5, 0.6349206349206349, 0.5, 0.5, 0.5263888888888889, 0., 0.6166666666666667, 0., 0.3611111111111111},
  SameTest -> Equal
]


MT[
  IGLocalEfficiency[dolphin, {}],
  {}
]


MT[
  IGAverageLocalEfficiency[dolphin],
  0.49103141908441456,
  SameTest -> Equal
]

MT[
  IGAverageLocalEfficiency[dolphin],
  Mean@IGLocalEfficiency[dolphin],
  SameTest -> Equal
]


(* ::Subsubsection::Closed:: *)
(*IGGlobalEfficiency*)


MT[
  IGGlobalEfficiency[IGEmptyGraph[0]],
  Indeterminate
]

MT[
  IGGlobalEfficiency[IGEmptyGraph[1]],
  Indeterminate
]

MT[
  IGGlobalEfficiency[IGEmptyGraph[2]],
  0.
]

MT[
  IGGlobalEfficiency[IGEmptyGraph[3]],
  0.
]


MT[
  IGGlobalEfficiency[dolphin],
  0.37921419757749875,
  SameTest -> Equal
]


(* ::Section::Closed:: *)
(*Bipartite graphs*)


MTSection["Bipartite graphs"]


(* ::Subsubsection::Closed:: *)
(*IGBipartiteQ*)


Function[graph,
  With[{g = graph},
    {MT[IGBipartiteQ[g], True],
     MT[IGBipartiteQ[g, IGBipartitePartitions[g]], True]}
  ],
  HoldAll
] /@ 
Hold[
  bipartite, dbipartite,
  IGEmptyGraph[0], IGEmptyGraph[5], CycleGraph[4],
  IGBipartiteGameGNM[10, 12, 33], IGBipartiteGameGNM[8, 17, 33, DirectedEdges -> True],
  IGBipartiteGameGNP[9, 6, 0.3], IGBipartiteGameGNP[21, 5, 0.6, DirectedEdges -> True]
] // ReleaseHold


MT[
  IGBipartiteQ@CycleGraph[3],
  False
]


(* graphs with self-loops are not bipartite *)
MT[
  IGBipartiteQ[IGShorthand["1-2-3-4-1-1", SelfLoops -> True]],
  False
]

MT[
  IGBipartiteQ[Graph[{1 <-> 1}]],
  False
]


(* ::Subsubsection::Closed:: *)
(*IGBipartitePartitions*)


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
  IGBipartitePartitions[dbipartite],
  {{1, 2, 3, 4, 5}, {6, 7, 8, 9, 10, 11, 12}}
]


MT[
  Sort /@ IGBipartitePartitions[IGShorthand["1-2-3-4-5"], 2],
  {{2, 4}, {1, 3, 5}}
]


MT[
  IGBipartitePartitions[IGShorthand["1-2-3-1"]],
  $Failed,
  {IGBipartitePartitions::nbipart}
]

MT[
  IGBipartitePartitions[IGShorthand["1-2-3-1"], 2],
  $Failed,
  {IGBipartitePartitions::nbipart}
]


(* ::Subsubsection::Closed:: *)
(*IGBipartiteProjections*)


bipartiteProjectionAM[g_, parts_] :=
    With[{im = IGBipartiteIncidenceMatrix[g, parts]},
      IGZeroDiagonal /@ Unitize[{im . Transpose[im], Transpose[im] . im}]
    ]

MT[
  Equal[
    AdjacencyMatrix /@ IGBipartiteProjections[bipartite, IGBipartitePartitions[bipartite]],
    bipartiteProjectionAM[bipartite, IGBipartitePartitions[bipartite]]
  ],
  True
]


(* ::Subsubsection::Closed:: *)
(*IGBipartiteIncidenceGraph*)


MT[
  IGBipartiteIncidenceGraph[{{1, 0, 1}, {1, 0, 0}}],
  Graph[{1, 2, 3, 4, 5}, {1 <-> 3, 1 <-> 5, 2 <-> 3}],
  SameTest -> IGSameGraphQ
]

MT[
  IGBipartiteIncidenceGraph[{{"a", "b"}, {"X", "Y", "Z"}}, {{1, 0, 1}, {1, 0, 0}}],
  Graph[{"a", "b", "X", "Y", "Z"}, {"a" <-> "X", "a" <-> "Z", "b" <-> "X"}],
  SameTest -> IGSameGraphQ
]

MT[
  Head[IGBipartiteIncidenceGraph[{{"a", "b"}, {"X", "Y"}}, {{1, 0, 1}, {1, 0, 0}}]],
  IGBipartiteIncidenceGraph,
  {IGBipartiteIncidenceGraph::bdsz}
]

MT[
  Head[IGBipartiteIncidenceGraph[{{}}]],
  IGBipartiteIncidenceGraph,
  {IGBipartiteIncidenceGraph::inv}
]


(* ::Subsubsection::Closed:: *)
(*IGBipartiteIncidenceMatrix*)


MT[
  IGBipartiteIncidenceMatrix[IGEmptyGraph[2]],
  $Failed,
  {IGBipartiteIncidenceMatrix::empty}
]

MT[
  Normal[IGBipartiteIncidenceMatrix[IGEmptyGraph[2], {{1}, {2}}]],
  {{0}}
]

MT[
  Normal[IGBipartiteIncidenceMatrix[Graph[{1, 2}, {1 <-> 2}]]],
  {{1}}
]

MT[
  Normal[IGBipartiteIncidenceMatrix[Graph[{1, 2}, {1 <-> 2}], {{"a"}, {"b"}}]],
  $Failed,
  {IGBipartiteQ::bdprt, IGBipartiteIncidenceMatrix::bdpart}
]


(* ::Section::Closed:: *)
(*Chordal graphs*)


MTSection["Chordal graphs"]


(* ::Subsubsection::Closed:: *)
(*IGChrodalQ*)


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


MT[
  IGChordalQ[IGShorthand["1-1-2-3-4-1", SelfLoops -> True, MultiEdges -> True]],
  False
]


(* ::Subsubsection::Closed:: *)
(*IGChordalCompletion*)


(* nothing to complete *)
MT[
  IGChordalCompletion[CycleGraph[3]],
  {}
]

MT[
  IGChordalCompletion[IGTriangularLattice[4]],
  {UndirectedEdge[8, 2], UndirectedEdge[9, 2], UndirectedEdge[6, 2]}
]

MT[
  IGChordalCompletion[CycleGraph[4]],
  {UndirectedEdge[3, 1]}
]

MT[
  IGChordalCompletion[IGEmptyGraph[]],
  {}
]

MT[
  IGChordalCompletion[IGEmptyGraph[1]],
  {}
]


MT[
  IGChordalCompletion[Graph[{1, 2}, {1 <-> 1}]],
  {}
]


(* ::Subsubsection::Closed:: *)
(*IGMaximumCardinalitySearch*)


MT[
  IGMaximumCardinalitySearch[Graph[{1, 2, 3, 4, 5}, {1 <-> 2, 1 <-> 5, 2 <-> 3, 2 <-> 4, 2 <-> 5, 3 <-> 4, 3 <-> 5}]],
  {5, 3, 2, 1, 4}
]

MT[
  IGMaximumCardinalitySearch[IGShorthand["A-B:C:I, B-A:C:D, C-A:B:E:H, D-B:E:F, E-C:D:F:H, F-D:E:G, G-F:H, H-C:E:G:I, I-A:H"]],
  {9, 4, 6, 8, 3, 5, 7, 2, 1}
]

MT[
  IGMaximumCardinalitySearch[GridGraph[{3, 4}]],
  {12, 5, 4, 11, 6, 3, 10, 7, 2, 9, 8, 1}
]


MT[
  IGMaximumCardinalitySearch[Graph[{1 <-> 2, 1 <-> 2}]],
  {2, 1}
]


(* ::Section::Closed:: *)
(*Mesh graphs*)


MTSection["Mesh graphs"]


(* Mesh operations *)

MT[
  With[{g = IGMeshGraph@BoundaryDiscretizeRegion[Rectangle[], MaxCellMeasure -> Infinity]},
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
  {{{{0, 1}  ~UndirectedEdge~ {0, 2}, {0, 1}  ~UndirectedEdge~ {0, 3}, {0, 1}  ~UndirectedEdge~ {0, 4}, {0,
    1}  ~UndirectedEdge~ {0, 5}, {0, 1}  ~UndirectedEdge~ {0, 6}, {0, 1}  ~UndirectedEdge~ {0, 7}, {0,
    1}  ~UndirectedEdge~ {0, 8}, {0, 2}  ~UndirectedEdge~ {0, 3}, {0, 2}  ~UndirectedEdge~ {0, 6}, {0,
    2}  ~UndirectedEdge~ {0, 7}, {0, 3}  ~UndirectedEdge~ {0, 4}, {0, 3}  ~UndirectedEdge~ {0, 7}, {0,
    3}  ~UndirectedEdge~ {0, 8}, {0, 4}  ~UndirectedEdge~ {0, 8}, {0, 5}  ~UndirectedEdge~ {0, 6}, {0,
    5}  ~UndirectedEdge~ {0, 7}, {0, 5}  ~UndirectedEdge~ {0, 8}, {0, 6}  ~UndirectedEdge~ {0, 7}, {0,
    7}  ~UndirectedEdge~ {0, 8}}, {{0, 1}  ~UndirectedEdge~ {1, 1}, {0, 1}  ~UndirectedEdge~ {1, 3}, {0,
    1}  ~UndirectedEdge~ {1, 4}, {0, 1}  ~UndirectedEdge~ {1, 8}, {0, 1}  ~UndirectedEdge~ {1, 9}, {0,
    1}  ~UndirectedEdge~ {1, 12}, {0, 1}  ~UndirectedEdge~ {1, 14}, {0, 2}  ~UndirectedEdge~ {1, 14}, {0,
    2}  ~UndirectedEdge~ {1, 15}, {0, 2}  ~UndirectedEdge~ {1, 16}, {0, 2}  ~UndirectedEdge~ {1, 17}, {0,
    3}  ~UndirectedEdge~ {1, 7}, {0, 3}  ~UndirectedEdge~ {1, 9}, {0, 3}  ~UndirectedEdge~ {1, 10}, {0,
    3}  ~UndirectedEdge~ {1, 13}, {0, 3}  ~UndirectedEdge~ {1, 16}, {0, 4}  ~UndirectedEdge~ {1, 10}, {0,
    4}  ~UndirectedEdge~ {1, 11}, {0, 4}  ~UndirectedEdge~ {1, 12}, {0, 5}  ~UndirectedEdge~ {1, 4}, {0,
    5}  ~UndirectedEdge~ {1, 5}, {0, 5}  ~UndirectedEdge~ {1, 6}, {0, 5}  ~UndirectedEdge~ {1, 19}, {0,
    6}  ~UndirectedEdge~ {1, 2}, {0, 6}  ~UndirectedEdge~ {1, 3}, {0, 6}  ~UndirectedEdge~ {1, 6}, {0,
    6}  ~UndirectedEdge~ {1, 17}, {0, 7}  ~UndirectedEdge~ {1, 1}, {0, 7}  ~UndirectedEdge~ {1, 2}, {0,
    7}  ~UndirectedEdge~ {1, 5}, {0, 7}  ~UndirectedEdge~ {1, 13}, {0, 7}  ~UndirectedEdge~ {1, 15}, {0,
    7}  ~UndirectedEdge~ {1, 18}, {0, 8}  ~UndirectedEdge~ {1, 7}, {0, 8}  ~UndirectedEdge~ {1, 8}, {0,
    8}  ~UndirectedEdge~ {1, 11}, {0, 8}  ~UndirectedEdge~ {1, 18}, {0, 8}  ~UndirectedEdge~ {1, 19}}, {{0,
    1}  ~UndirectedEdge~ {2, 1}, {0, 1}  ~UndirectedEdge~ {2, 2}, {0, 1}  ~UndirectedEdge~ {2, 3}, {0,
    1}  ~UndirectedEdge~ {2, 5}, {0, 1}  ~UndirectedEdge~ {2, 7}, {0, 1}  ~UndirectedEdge~ {2, 8}, {0,
    1}  ~UndirectedEdge~ {2, 9}, {0, 1}  ~UndirectedEdge~ {2, 10}, {0, 1}  ~UndirectedEdge~ {2, 11}, {0,
    1}  ~UndirectedEdge~ {2, 14}, {0, 1}  ~UndirectedEdge~ {2, 16}, {0, 1}  ~UndirectedEdge~ {2, 17}, {0,
    2}  ~UndirectedEdge~ {2, 10}, {0, 2}  ~UndirectedEdge~ {2, 11}, {0, 2}  ~UndirectedEdge~ {2, 12}, {0,
    2}  ~UndirectedEdge~ {2, 13}, {0, 2}  ~UndirectedEdge~ {2, 14}, {0, 3}  ~UndirectedEdge~ {2, 5}, {0,
    3}  ~UndirectedEdge~ {2, 6}, {0, 3}  ~UndirectedEdge~ {2, 7}, {0, 3}  ~UndirectedEdge~ {2, 9}, {0,
    3}  ~UndirectedEdge~ {2, 11}, {0, 3}  ~UndirectedEdge~ {2, 12}, {0, 3}  ~UndirectedEdge~ {2, 15}, {0,
    4}  ~UndirectedEdge~ {2, 6}, {0, 4}  ~UndirectedEdge~ {2, 7}, {0, 4}  ~UndirectedEdge~ {2, 8}, {0,
    5}  ~UndirectedEdge~ {2, 2}, {0, 5}  ~UndirectedEdge~ {2, 3}, {0, 5}  ~UndirectedEdge~ {2, 4}, {0,
    5}  ~UndirectedEdge~ {2, 17}, {0, 5}  ~UndirectedEdge~ {2, 18}, {0, 6}  ~UndirectedEdge~ {2, 1}, {0,
    6}  ~UndirectedEdge~ {2, 3}, {0, 6}  ~UndirectedEdge~ {2, 4}, {0, 6}  ~UndirectedEdge~ {2, 13}, {0,
    6}  ~UndirectedEdge~ {2, 14}, {0, 7}  ~UndirectedEdge~ {2, 1}, {0, 7}  ~UndirectedEdge~ {2, 2}, {0,
    7}  ~UndirectedEdge~ {2, 4}, {0, 7}  ~UndirectedEdge~ {2, 9}, {0, 7}  ~UndirectedEdge~ {2, 10}, {0,
    7}  ~UndirectedEdge~ {2, 12}, {0, 7}  ~UndirectedEdge~ {2, 13}, {0, 7}  ~UndirectedEdge~ {2, 15}, {0,
    7}  ~UndirectedEdge~ {2, 16}, {0, 7}  ~UndirectedEdge~ {2, 18}, {0, 8}  ~UndirectedEdge~ {2, 5}, {0,
    8}  ~UndirectedEdge~ {2, 6}, {0, 8}  ~UndirectedEdge~ {2, 8}, {0, 8}  ~UndirectedEdge~ {2, 15}, {0,
    8}  ~UndirectedEdge~ {2, 16}, {0, 8}  ~UndirectedEdge~ {2, 17}, {0, 8}  ~UndirectedEdge~ {2, 18}}, {{0,
    1}  ~UndirectedEdge~ {3, 1}, {0, 1}  ~UndirectedEdge~ {3, 2}, {0, 1}  ~UndirectedEdge~ {3, 3}, {0,
    1}  ~UndirectedEdge~ {3, 4}, {0, 1}  ~UndirectedEdge~ {3, 5}, {0, 1}  ~UndirectedEdge~ {3, 6}, {0,
    2}  ~UndirectedEdge~ {3, 3}, {0, 2}  ~UndirectedEdge~ {3, 4}, {0, 3}  ~UndirectedEdge~ {3, 2}, {0,
    3}  ~UndirectedEdge~ {3, 3}, {0, 3}  ~UndirectedEdge~ {3, 5}, {0, 4}  ~UndirectedEdge~ {3, 2}, {0,
    5}  ~UndirectedEdge~ {3, 1}, {0, 5}  ~UndirectedEdge~ {3, 6}, {0, 6}  ~UndirectedEdge~ {3, 1}, {0,
    6}  ~UndirectedEdge~ {3, 4}, {0, 7}  ~UndirectedEdge~ {3, 1}, {0, 7}  ~UndirectedEdge~ {3, 3}, {0,
    7}  ~UndirectedEdge~ {3, 4}, {0, 7}  ~UndirectedEdge~ {3, 5}, {0, 7}  ~UndirectedEdge~ {3, 6}, {0,
    8}  ~UndirectedEdge~ {3, 2}, {0, 8}  ~UndirectedEdge~ {3, 5}, {0, 8}  ~UndirectedEdge~ {3, 6}}}, {{{1,
    1}  ~UndirectedEdge~ {0, 1}, {1, 1}  ~UndirectedEdge~ {0, 7}, {1, 2}  ~UndirectedEdge~ {0, 6}, {1,
    2}  ~UndirectedEdge~ {0, 7}, {1, 3}  ~UndirectedEdge~ {0, 1}, {1, 3}  ~UndirectedEdge~ {0, 6}, {1,
    4}  ~UndirectedEdge~ {0, 1}, {1, 4}  ~UndirectedEdge~ {0, 5}, {1, 5}  ~UndirectedEdge~ {0, 5}, {1,
    5}  ~UndirectedEdge~ {0, 7}, {1, 6}  ~UndirectedEdge~ {0, 5}, {1, 6}  ~UndirectedEdge~ {0, 6}, {1,
    7}  ~UndirectedEdge~ {0, 3}, {1, 7}  ~UndirectedEdge~ {0, 8}, {1, 8}  ~UndirectedEdge~ {0, 1}, {1,
    8}  ~UndirectedEdge~ {0, 8}, {1, 9}  ~UndirectedEdge~ {0, 1}, {1, 9}  ~UndirectedEdge~ {0, 3}, {1,
    10}  ~UndirectedEdge~ {0, 3}, {1, 10}  ~UndirectedEdge~ {0, 4}, {1, 11}  ~UndirectedEdge~ {0, 4}, {1,
    11}  ~UndirectedEdge~ {0, 8}, {1, 12}  ~UndirectedEdge~ {0, 1}, {1, 12}  ~UndirectedEdge~ {0, 4}, {1,
    13}  ~UndirectedEdge~ {0, 3}, {1, 13}  ~UndirectedEdge~ {0, 7}, {1, 14}  ~UndirectedEdge~ {0, 1}, {1,
    14}  ~UndirectedEdge~ {0, 2}, {1, 15}  ~UndirectedEdge~ {0, 2}, {1, 15}  ~UndirectedEdge~ {0, 7}, {1,
    16}  ~UndirectedEdge~ {0, 2}, {1, 16}  ~UndirectedEdge~ {0, 3}, {1, 17}  ~UndirectedEdge~ {0, 2}, {1,
    17}  ~UndirectedEdge~ {0, 6}, {1, 18}  ~UndirectedEdge~ {0, 7}, {1, 18}  ~UndirectedEdge~ {0, 8}, {1,
    19}  ~UndirectedEdge~ {0, 5}, {1, 19}  ~UndirectedEdge~ {0, 8}}, {{1, 1}  ~UndirectedEdge~ {1, 2}, {1,
    1}  ~UndirectedEdge~ {1, 3}, {1, 1}  ~UndirectedEdge~ {1, 4}, {1, 1}  ~UndirectedEdge~ {1, 5}, {1,
    1}  ~UndirectedEdge~ {1, 8}, {1, 1}  ~UndirectedEdge~ {1, 9}, {1, 1}  ~UndirectedEdge~ {1, 12}, {1,
    1}  ~UndirectedEdge~ {1, 13}, {1, 1}  ~UndirectedEdge~ {1, 14}, {1, 1}  ~UndirectedEdge~ {1, 15}, {1,
    1}  ~UndirectedEdge~ {1, 18}, {1, 2}  ~UndirectedEdge~ {1, 3}, {1, 2}  ~UndirectedEdge~ {1, 5}, {1,
    2}  ~UndirectedEdge~ {1, 6}, {1, 2}  ~UndirectedEdge~ {1, 13}, {1, 2}  ~UndirectedEdge~ {1, 15}, {1,
    2}  ~UndirectedEdge~ {1, 17}, {1, 2}  ~UndirectedEdge~ {1, 18}, {1, 3}  ~UndirectedEdge~ {1, 4}, {1,
    3}  ~UndirectedEdge~ {1, 6}, {1, 3}  ~UndirectedEdge~ {1, 8}, {1, 3}  ~UndirectedEdge~ {1, 9}, {1,
    3}  ~UndirectedEdge~ {1, 12}, {1, 3}  ~UndirectedEdge~ {1, 14}, {1, 3}  ~UndirectedEdge~ {1, 17}, {1,
    4}  ~UndirectedEdge~ {1, 5}, {1, 4}  ~UndirectedEdge~ {1, 6}, {1, 4}  ~UndirectedEdge~ {1, 8}, {1,
    4}  ~UndirectedEdge~ {1, 9}, {1, 4}  ~UndirectedEdge~ {1, 12}, {1, 4}  ~UndirectedEdge~ {1, 14}, {1,
    4}  ~UndirectedEdge~ {1, 19}, {1, 5}  ~UndirectedEdge~ {1, 6}, {1, 5}  ~UndirectedEdge~ {1, 13}, {1,
    5}  ~UndirectedEdge~ {1, 15}, {1, 5}  ~UndirectedEdge~ {1, 18}, {1, 5}  ~UndirectedEdge~ {1, 19}, {1,
    6}  ~UndirectedEdge~ {1, 17}, {1, 6}  ~UndirectedEdge~ {1, 19}, {1, 7}  ~UndirectedEdge~ {1, 8}, {1,
    7}  ~UndirectedEdge~ {1, 9}, {1, 7}  ~UndirectedEdge~ {1, 10}, {1, 7}  ~UndirectedEdge~ {1, 11}, {1,
    7}  ~UndirectedEdge~ {1, 13}, {1, 7}  ~UndirectedEdge~ {1, 16}, {1, 7}  ~UndirectedEdge~ {1, 18}, {1,
    7}  ~UndirectedEdge~ {1, 19}, {1, 8}  ~UndirectedEdge~ {1, 9}, {1, 8}  ~UndirectedEdge~ {1, 11}, {1,
    8}  ~UndirectedEdge~ {1, 12}, {1, 8}  ~UndirectedEdge~ {1, 14}, {1, 8}  ~UndirectedEdge~ {1, 18}, {1,
    8}  ~UndirectedEdge~ {1, 19}, {1, 9}  ~UndirectedEdge~ {1, 10}, {1, 9}  ~UndirectedEdge~ {1, 12}, {1,
    9}  ~UndirectedEdge~ {1, 13}, {1, 9}  ~UndirectedEdge~ {1, 14}, {1, 9}  ~UndirectedEdge~ {1, 16}, {1,
    10}  ~UndirectedEdge~ {1, 11}, {1, 10}  ~UndirectedEdge~ {1, 12}, {1, 10}  ~UndirectedEdge~ {1, 13}, {1,
    10}  ~UndirectedEdge~ {1, 16}, {1, 11}  ~UndirectedEdge~ {1, 12}, {1, 11}  ~UndirectedEdge~ {1, 18}, {1,
    11}  ~UndirectedEdge~ {1, 19}, {1, 12}  ~UndirectedEdge~ {1, 14}, {1, 13}  ~UndirectedEdge~ {1, 15}, {1,
    13}  ~UndirectedEdge~ {1, 16}, {1, 13}  ~UndirectedEdge~ {1, 18}, {1, 14}  ~UndirectedEdge~ {1, 15}, {1,
    14}  ~UndirectedEdge~ {1, 16}, {1, 14}  ~UndirectedEdge~ {1, 17}, {1, 15}  ~UndirectedEdge~ {1, 16}, {1,
    15}  ~UndirectedEdge~ {1, 17}, {1, 15}  ~UndirectedEdge~ {1, 18}, {1, 16}  ~UndirectedEdge~ {1, 17}, {1,
    18}  ~UndirectedEdge~ {1, 19}}, {{1, 1}  ~UndirectedEdge~ {2, 1}, {1, 1}  ~UndirectedEdge~ {2, 2}, {1,
    1}  ~UndirectedEdge~ {2, 9}, {1, 1}  ~UndirectedEdge~ {2, 10}, {1, 1}  ~UndirectedEdge~ {2, 16}, {1,
    2}  ~UndirectedEdge~ {2, 1}, {1, 2}  ~UndirectedEdge~ {2, 4}, {1, 2}  ~UndirectedEdge~ {2, 13}, {1,
    3}  ~UndirectedEdge~ {2, 1}, {1, 3}  ~UndirectedEdge~ {2, 3}, {1, 3}  ~UndirectedEdge~ {2, 14}, {1,
    4}  ~UndirectedEdge~ {2, 2}, {1, 4}  ~UndirectedEdge~ {2, 3}, {1, 4}  ~UndirectedEdge~ {2, 17}, {1,
    5}  ~UndirectedEdge~ {2, 2}, {1, 5}  ~UndirectedEdge~ {2, 4}, {1, 5}  ~UndirectedEdge~ {2, 18}, {1,
    6}  ~UndirectedEdge~ {2, 3}, {1, 6}  ~UndirectedEdge~ {2, 4}, {1, 7}  ~UndirectedEdge~ {2, 5}, {1,
    7}  ~UndirectedEdge~ {2, 6}, {1, 7}  ~UndirectedEdge~ {2, 15}, {1, 8}  ~UndirectedEdge~ {2, 5}, {1,
    8}  ~UndirectedEdge~ {2, 8}, {1, 8}  ~UndirectedEdge~ {2, 16}, {1, 8}  ~UndirectedEdge~ {2, 17}, {1,
    9}  ~UndirectedEdge~ {2, 5}, {1, 9}  ~UndirectedEdge~ {2, 7}, {1, 9}  ~UndirectedEdge~ {2, 9}, {1,
    9}  ~UndirectedEdge~ {2, 11}, {1, 10}  ~UndirectedEdge~ {2, 6}, {1, 10}  ~UndirectedEdge~ {2, 7}, {1,
    11}  ~UndirectedEdge~ {2, 6}, {1, 11}  ~UndirectedEdge~ {2, 8}, {1, 12}  ~UndirectedEdge~ {2, 7}, {1,
    12}  ~UndirectedEdge~ {2, 8}, {1, 13}  ~UndirectedEdge~ {2, 9}, {1, 13}  ~UndirectedEdge~ {2, 12}, {1,
    13}  ~UndirectedEdge~ {2, 15}, {1, 14}  ~UndirectedEdge~ {2, 10}, {1, 14}  ~UndirectedEdge~ {2, 11}, {1,
    14}  ~UndirectedEdge~ {2, 14}, {1, 15}  ~UndirectedEdge~ {2, 10}, {1, 15}  ~UndirectedEdge~ {2, 12}, {1,
    15}  ~UndirectedEdge~ {2, 13}, {1, 16}  ~UndirectedEdge~ {2, 11}, {1, 16}  ~UndirectedEdge~ {2, 12}, {1,
    17}  ~UndirectedEdge~ {2, 13}, {1, 17}  ~UndirectedEdge~ {2, 14}, {1, 18}  ~UndirectedEdge~ {2, 15}, {1,
    18}  ~UndirectedEdge~ {2, 16}, {1, 18}  ~UndirectedEdge~ {2, 18}, {1, 19}  ~UndirectedEdge~ {2, 17}, {1,
    19}  ~UndirectedEdge~ {2, 18}}, {{1, 1}  ~UndirectedEdge~ {3, 1}, {1, 1}  ~UndirectedEdge~ {3, 3}, {1,
    1}  ~UndirectedEdge~ {3, 4}, {1, 1}  ~UndirectedEdge~ {3, 5}, {1, 1}  ~UndirectedEdge~ {3, 6}, {1,
    2}  ~UndirectedEdge~ {3, 1}, {1, 2}  ~UndirectedEdge~ {3, 4}, {1, 3}  ~UndirectedEdge~ {3, 1}, {1,
    3}  ~UndirectedEdge~ {3, 4}, {1, 4}  ~UndirectedEdge~ {3, 1}, {1, 4}  ~UndirectedEdge~ {3, 6}, {1,
    5}  ~UndirectedEdge~ {3, 1}, {1, 5}  ~UndirectedEdge~ {3, 6}, {1, 6}  ~UndirectedEdge~ {3, 1}, {1,
    7}  ~UndirectedEdge~ {3, 2}, {1, 7}  ~UndirectedEdge~ {3, 5}, {1, 8}  ~UndirectedEdge~ {3, 2}, {1,
    8}  ~UndirectedEdge~ {3, 5}, {1, 8}  ~UndirectedEdge~ {3, 6}, {1, 9}  ~UndirectedEdge~ {3, 2}, {1,
    9}  ~UndirectedEdge~ {3, 3}, {1, 9}  ~UndirectedEdge~ {3, 5}, {1, 10}  ~UndirectedEdge~ {3, 2}, {1,
    11}  ~UndirectedEdge~ {3, 2}, {1, 12}  ~UndirectedEdge~ {3, 2}, {1, 13}  ~UndirectedEdge~ {3, 3}, {1,
    13}  ~UndirectedEdge~ {3, 5}, {1, 14}  ~UndirectedEdge~ {3, 3}, {1, 14}  ~UndirectedEdge~ {3, 4}, {1,
    15}  ~UndirectedEdge~ {3, 3}, {1, 15}  ~UndirectedEdge~ {3, 4}, {1, 16}  ~UndirectedEdge~ {3, 3}, {1,
    17}  ~UndirectedEdge~ {3, 4}, {1, 18}  ~UndirectedEdge~ {3, 5}, {1, 18}  ~UndirectedEdge~ {3, 6}, {1,
    19}  ~UndirectedEdge~ {3, 6}}}, {{{2, 1}  ~UndirectedEdge~ {0, 1}, {2, 1}  ~UndirectedEdge~ {0, 6}, {2,
    1}  ~UndirectedEdge~ {0, 7}, {2, 2}  ~UndirectedEdge~ {0, 1}, {2, 2}  ~UndirectedEdge~ {0, 5}, {2,
    2}  ~UndirectedEdge~ {0, 7}, {2, 3}  ~UndirectedEdge~ {0, 1}, {2, 3}  ~UndirectedEdge~ {0, 5}, {2,
    3}  ~UndirectedEdge~ {0, 6}, {2, 4}  ~UndirectedEdge~ {0, 5}, {2, 4}  ~UndirectedEdge~ {0, 6}, {2,
    4}  ~UndirectedEdge~ {0, 7}, {2, 5}  ~UndirectedEdge~ {0, 1}, {2, 5}  ~UndirectedEdge~ {0, 3}, {2,
    5}  ~UndirectedEdge~ {0, 8}, {2, 6}  ~UndirectedEdge~ {0, 3}, {2, 6}  ~UndirectedEdge~ {0, 4}, {2,
    6}  ~UndirectedEdge~ {0, 8}, {2, 7}  ~UndirectedEdge~ {0, 1}, {2, 7}  ~UndirectedEdge~ {0, 3}, {2,
    7}  ~UndirectedEdge~ {0, 4}, {2, 8}  ~UndirectedEdge~ {0, 1}, {2, 8}  ~UndirectedEdge~ {0, 4}, {2,
    8}  ~UndirectedEdge~ {0, 8}, {2, 9}  ~UndirectedEdge~ {0, 1}, {2, 9}  ~UndirectedEdge~ {0, 3}, {2,
    9}  ~UndirectedEdge~ {0, 7}, {2, 10}  ~UndirectedEdge~ {0, 1}, {2, 10}  ~UndirectedEdge~ {0, 2}, {2,
    10}  ~UndirectedEdge~ {0, 7}, {2, 11}  ~UndirectedEdge~ {0, 1}, {2, 11}  ~UndirectedEdge~ {0, 2}, {2,
    11}  ~UndirectedEdge~ {0, 3}, {2, 12}  ~UndirectedEdge~ {0, 2}, {2, 12}  ~UndirectedEdge~ {0, 3}, {2,
    12}  ~UndirectedEdge~ {0, 7}, {2, 13}  ~UndirectedEdge~ {0, 2}, {2, 13}  ~UndirectedEdge~ {0, 6}, {2,
    13}  ~UndirectedEdge~ {0, 7}, {2, 14}  ~UndirectedEdge~ {0, 1}, {2, 14}  ~UndirectedEdge~ {0, 2}, {2,
    14}  ~UndirectedEdge~ {0, 6}, {2, 15}  ~UndirectedEdge~ {0, 3}, {2, 15}  ~UndirectedEdge~ {0, 7}, {2,
    15}  ~UndirectedEdge~ {0, 8}, {2, 16}  ~UndirectedEdge~ {0, 1}, {2, 16}  ~UndirectedEdge~ {0, 7}, {2,
    16}  ~UndirectedEdge~ {0, 8}, {2, 17}  ~UndirectedEdge~ {0, 1}, {2, 17}  ~UndirectedEdge~ {0, 5}, {2,
    17}  ~UndirectedEdge~ {0, 8}, {2, 18}  ~UndirectedEdge~ {0, 5}, {2, 18}  ~UndirectedEdge~ {0, 7}, {2,
    18}  ~UndirectedEdge~ {0, 8}}, {{2, 1}  ~UndirectedEdge~ {1, 1}, {2, 1}  ~UndirectedEdge~ {1, 2}, {2,
    1}  ~UndirectedEdge~ {1, 3}, {2, 2}  ~UndirectedEdge~ {1, 1}, {2, 2}  ~UndirectedEdge~ {1, 4}, {2,
    2}  ~UndirectedEdge~ {1, 5}, {2, 3}  ~UndirectedEdge~ {1, 3}, {2, 3}  ~UndirectedEdge~ {1, 4}, {2,
    3}  ~UndirectedEdge~ {1, 6}, {2, 4}  ~UndirectedEdge~ {1, 2}, {2, 4}  ~UndirectedEdge~ {1, 5}, {2,
    4}  ~UndirectedEdge~ {1, 6}, {2, 5}  ~UndirectedEdge~ {1, 7}, {2, 5}  ~UndirectedEdge~ {1, 8}, {2,
    5}  ~UndirectedEdge~ {1, 9}, {2, 6}  ~UndirectedEdge~ {1, 7}, {2, 6}  ~UndirectedEdge~ {1, 10}, {2,
    6}  ~UndirectedEdge~ {1, 11}, {2, 7}  ~UndirectedEdge~ {1, 9}, {2, 7}  ~UndirectedEdge~ {1, 10}, {2,
    7}  ~UndirectedEdge~ {1, 12}, {2, 8}  ~UndirectedEdge~ {1, 8}, {2, 8}  ~UndirectedEdge~ {1, 11}, {2,
    8}  ~UndirectedEdge~ {1, 12}, {2, 9}  ~UndirectedEdge~ {1, 1}, {2, 9}  ~UndirectedEdge~ {1, 9}, {2,
    9}  ~UndirectedEdge~ {1, 13}, {2, 10}  ~UndirectedEdge~ {1, 1}, {2, 10}  ~UndirectedEdge~ {1, 14}, {2,
    10}  ~UndirectedEdge~ {1, 15}, {2, 11}  ~UndirectedEdge~ {1, 9}, {2, 11}  ~UndirectedEdge~ {1, 14}, {2,
    11}  ~UndirectedEdge~ {1, 16}, {2, 12}  ~UndirectedEdge~ {1, 13}, {2, 12}  ~UndirectedEdge~ {1, 15}, {2,
    12}  ~UndirectedEdge~ {1, 16}, {2, 13}  ~UndirectedEdge~ {1, 2}, {2, 13}  ~UndirectedEdge~ {1, 15}, {2,
    13}  ~UndirectedEdge~ {1, 17}, {2, 14}  ~UndirectedEdge~ {1, 3}, {2, 14}  ~UndirectedEdge~ {1, 14}, {2,
    14}  ~UndirectedEdge~ {1, 17}, {2, 15}  ~UndirectedEdge~ {1, 7}, {2, 15}  ~UndirectedEdge~ {1, 13}, {2,
    15}  ~UndirectedEdge~ {1, 18}, {2, 16}  ~UndirectedEdge~ {1, 1}, {2, 16}  ~UndirectedEdge~ {1, 8}, {2,
    16}  ~UndirectedEdge~ {1, 18}, {2, 17}  ~UndirectedEdge~ {1, 4}, {2, 17}  ~UndirectedEdge~ {1, 8}, {2,
    17}  ~UndirectedEdge~ {1, 19}, {2, 18}  ~UndirectedEdge~ {1, 5}, {2, 18}  ~UndirectedEdge~ {1, 18}, {2,
    18}  ~UndirectedEdge~ {1, 19}}, {{2, 1}  ~UndirectedEdge~ {2, 2}, {2, 1}  ~UndirectedEdge~ {2, 3}, {2,
    1}  ~UndirectedEdge~ {2, 4}, {2, 1}  ~UndirectedEdge~ {2, 9}, {2, 1}  ~UndirectedEdge~ {2, 10}, {2,
    1}  ~UndirectedEdge~ {2, 13}, {2, 1}  ~UndirectedEdge~ {2, 14}, {2, 1}  ~UndirectedEdge~ {2, 16}, {2,
    2}  ~UndirectedEdge~ {2, 3}, {2, 2}  ~UndirectedEdge~ {2, 4}, {2, 2}  ~UndirectedEdge~ {2, 9}, {2,
    2}  ~UndirectedEdge~ {2, 10}, {2, 2}  ~UndirectedEdge~ {2, 16}, {2, 2}  ~UndirectedEdge~ {2, 17}, {2,
    2}  ~UndirectedEdge~ {2, 18}, {2, 3}  ~UndirectedEdge~ {2, 4}, {2, 3}  ~UndirectedEdge~ {2, 14}, {2,
    3}  ~UndirectedEdge~ {2, 17}, {2, 4}  ~UndirectedEdge~ {2, 13}, {2, 4}  ~UndirectedEdge~ {2, 18}, {2,
    5}  ~UndirectedEdge~ {2, 6}, {2, 5}  ~UndirectedEdge~ {2, 7}, {2, 5}  ~UndirectedEdge~ {2, 8}, {2,
    5}  ~UndirectedEdge~ {2, 9}, {2, 5}  ~UndirectedEdge~ {2, 11}, {2, 5}  ~UndirectedEdge~ {2, 15}, {2,
    5}  ~UndirectedEdge~ {2, 16}, {2, 5}  ~UndirectedEdge~ {2, 17}, {2, 6}  ~UndirectedEdge~ {2, 7}, {2,
    6}  ~UndirectedEdge~ {2, 8}, {2, 6}  ~UndirectedEdge~ {2, 15}, {2, 7}  ~UndirectedEdge~ {2, 8}, {2,
    7}  ~UndirectedEdge~ {2, 9}, {2, 7}  ~UndirectedEdge~ {2, 11}, {2, 8}  ~UndirectedEdge~ {2, 16}, {2,
    8}  ~UndirectedEdge~ {2, 17}, {2, 9}  ~UndirectedEdge~ {2, 10}, {2, 9}  ~UndirectedEdge~ {2, 11}, {2,
    9}  ~UndirectedEdge~ {2, 12}, {2, 9}  ~UndirectedEdge~ {2, 15}, {2, 9}  ~UndirectedEdge~ {2, 16}, {2,
    10}  ~UndirectedEdge~ {2, 11}, {2, 10}  ~UndirectedEdge~ {2, 12}, {2, 10}  ~UndirectedEdge~ {2, 13}, {2,
    10}  ~UndirectedEdge~ {2, 14}, {2, 10}  ~UndirectedEdge~ {2, 16}, {2, 11}  ~UndirectedEdge~ {2, 12}, {2,
    11}  ~UndirectedEdge~ {2, 14}, {2, 12}  ~UndirectedEdge~ {2, 13}, {2, 12}  ~UndirectedEdge~ {2, 15}, {2,
    13}  ~UndirectedEdge~ {2, 14}, {2, 15}  ~UndirectedEdge~ {2, 16}, {2, 15}  ~UndirectedEdge~ {2, 18}, {2,
    16}  ~UndirectedEdge~ {2, 17}, {2, 16}  ~UndirectedEdge~ {2, 18}, {2, 17}  ~UndirectedEdge~ {2, 18}}, {{2,
    1}  ~UndirectedEdge~ {3, 1}, {2, 1}  ~UndirectedEdge~ {3, 4}, {2, 2}  ~UndirectedEdge~ {3, 1}, {2,
    2}  ~UndirectedEdge~ {3, 6}, {2, 3}  ~UndirectedEdge~ {3, 1}, {2, 4}  ~UndirectedEdge~ {3, 1}, {2,
    5}  ~UndirectedEdge~ {3, 2}, {2, 5}  ~UndirectedEdge~ {3, 5}, {2, 6}  ~UndirectedEdge~ {3, 2}, {2,
    7}  ~UndirectedEdge~ {3, 2}, {2, 8}  ~UndirectedEdge~ {3, 2}, {2, 9}  ~UndirectedEdge~ {3, 3}, {2,
    9}  ~UndirectedEdge~ {3, 5}, {2, 10}  ~UndirectedEdge~ {3, 3}, {2, 10}  ~UndirectedEdge~ {3, 4}, {2,
    11}  ~UndirectedEdge~ {3, 3}, {2, 12}  ~UndirectedEdge~ {3, 3}, {2, 13}  ~UndirectedEdge~ {3, 4}, {2,
    14}  ~UndirectedEdge~ {3, 4}, {2, 15}  ~UndirectedEdge~ {3, 5}, {2, 16}  ~UndirectedEdge~ {3, 5}, {2,
    16}  ~UndirectedEdge~ {3, 6}, {2, 17}  ~UndirectedEdge~ {3, 6}, {2, 18}  ~UndirectedEdge~ {3, 6}}}, {{{3,
    1}  ~UndirectedEdge~ {0, 1}, {3, 1}  ~UndirectedEdge~ {0, 5}, {3, 1}  ~UndirectedEdge~ {0, 6}, {3,
    1}  ~UndirectedEdge~ {0, 7}, {3, 2}  ~UndirectedEdge~ {0, 1}, {3, 2}  ~UndirectedEdge~ {0, 3}, {3,
    2}  ~UndirectedEdge~ {0, 4}, {3, 2}  ~UndirectedEdge~ {0, 8}, {3, 3}  ~UndirectedEdge~ {0, 1}, {3,
    3}  ~UndirectedEdge~ {0, 2}, {3, 3}  ~UndirectedEdge~ {0, 3}, {3, 3}  ~UndirectedEdge~ {0, 7}, {3,
    4}  ~UndirectedEdge~ {0, 1}, {3, 4}  ~UndirectedEdge~ {0, 2}, {3, 4}  ~UndirectedEdge~ {0, 6}, {3,
    4}  ~UndirectedEdge~ {0, 7}, {3, 5}  ~UndirectedEdge~ {0, 1}, {3, 5}  ~UndirectedEdge~ {0, 3}, {3,
    5}  ~UndirectedEdge~ {0, 7}, {3, 5}  ~UndirectedEdge~ {0, 8}, {3, 6}  ~UndirectedEdge~ {0, 1}, {3,
    6}  ~UndirectedEdge~ {0, 5}, {3, 6}  ~UndirectedEdge~ {0, 7}, {3, 6}  ~UndirectedEdge~ {0, 8}}, {{3,
    1}  ~UndirectedEdge~ {1, 1}, {3, 1}  ~UndirectedEdge~ {1, 2}, {3, 1}  ~UndirectedEdge~ {1, 3}, {3,
    1}  ~UndirectedEdge~ {1, 4}, {3, 1}  ~UndirectedEdge~ {1, 5}, {3, 1}  ~UndirectedEdge~ {1, 6}, {3,
    2}  ~UndirectedEdge~ {1, 7}, {3, 2}  ~UndirectedEdge~ {1, 8}, {3, 2}  ~UndirectedEdge~ {1, 9}, {3,
    2}  ~UndirectedEdge~ {1, 10}, {3, 2}  ~UndirectedEdge~ {1, 11}, {3, 2}  ~UndirectedEdge~ {1, 12}, {3,
    3}  ~UndirectedEdge~ {1, 1}, {3, 3}  ~UndirectedEdge~ {1, 9}, {3, 3}  ~UndirectedEdge~ {1, 13}, {3,
    3}  ~UndirectedEdge~ {1, 14}, {3, 3}  ~UndirectedEdge~ {1, 15}, {3, 3}  ~UndirectedEdge~ {1, 16}, {3,
    4}  ~UndirectedEdge~ {1, 1}, {3, 4}  ~UndirectedEdge~ {1, 2}, {3, 4}  ~UndirectedEdge~ {1, 3}, {3,
    4}  ~UndirectedEdge~ {1, 14}, {3, 4}  ~UndirectedEdge~ {1, 15}, {3, 4}  ~UndirectedEdge~ {1, 17}, {3,
    5}  ~UndirectedEdge~ {1, 1}, {3, 5}  ~UndirectedEdge~ {1, 7}, {3, 5}  ~UndirectedEdge~ {1, 8}, {3,
    5}  ~UndirectedEdge~ {1, 9}, {3, 5}  ~UndirectedEdge~ {1, 13}, {3, 5}  ~UndirectedEdge~ {1, 18}, {3,
    6}  ~UndirectedEdge~ {1, 1}, {3, 6}  ~UndirectedEdge~ {1, 4}, {3, 6}  ~UndirectedEdge~ {1, 5}, {3,
    6}  ~UndirectedEdge~ {1, 8}, {3, 6}  ~UndirectedEdge~ {1, 18}, {3, 6}  ~UndirectedEdge~ {1, 19}}, {{3,
    1}  ~UndirectedEdge~ {2, 1}, {3, 1}  ~UndirectedEdge~ {2, 2}, {3, 1}  ~UndirectedEdge~ {2, 3}, {3,
    1}  ~UndirectedEdge~ {2, 4}, {3, 2}  ~UndirectedEdge~ {2, 5}, {3, 2}  ~UndirectedEdge~ {2, 6}, {3,
    2}  ~UndirectedEdge~ {2, 7}, {3, 2}  ~UndirectedEdge~ {2, 8}, {3, 3}  ~UndirectedEdge~ {2, 9}, {3,
    3}  ~UndirectedEdge~ {2, 10}, {3, 3}  ~UndirectedEdge~ {2, 11}, {3, 3}  ~UndirectedEdge~ {2, 12}, {3,
    4}  ~UndirectedEdge~ {2, 1}, {3, 4}  ~UndirectedEdge~ {2, 10}, {3, 4}  ~UndirectedEdge~ {2, 13}, {3,
    4}  ~UndirectedEdge~ {2, 14}, {3, 5}  ~UndirectedEdge~ {2, 5}, {3, 5}  ~UndirectedEdge~ {2, 9}, {3,
    5}  ~UndirectedEdge~ {2, 15}, {3, 5}  ~UndirectedEdge~ {2, 16}, {3, 6}  ~UndirectedEdge~ {2, 2}, {3,
    6}  ~UndirectedEdge~ {2, 16}, {3, 6}  ~UndirectedEdge~ {2, 17}, {3, 6}  ~UndirectedEdge~ {2, 18}}, {{3,
    1}  ~UndirectedEdge~ {3, 4}, {3, 1}  ~UndirectedEdge~ {3, 6}, {3, 2}  ~UndirectedEdge~ {3, 5}, {3,
    3}  ~UndirectedEdge~ {3, 4}, {3, 3}  ~UndirectedEdge~ {3, 5}, {3, 5}  ~UndirectedEdge~ {3, 6}}}}
]


(* ::Subsubsection:: *)
(*IGMeshCellAdjacencyMatrix*)


(* TODO many tests needed before including new implementation *)
(* TODO need to test on 1D, 2D, 3D and mixed dimensional meshes in various embedding dimensions *)


(* ::Section::Closed:: *)
(*Graph colouring*)


MTSection["Graph colouring"]


(* ::Subsubsection::Closed:: *)
(*IGVertexColoring*)


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


(* self-loops must be ignored *)
MT[
  IGVertexColoring[Graph[{1 <-> 2, 1 <-> 2, 2 <-> 2}]],
  {1, 2} (* note: {2,1} is also good *)
]


(* ::Subsubsection::Closed:: *)
(*IGEdgeColoring*)


MT[
  IGEdgeColoring@IGEmptyGraph[#],
  {}
]& /@ {0, 1, 2, 3}


MT[
  IGEdgeColoring[CycleGraph[3]],
  {1, 3, 2}
]


MT[
  IGEdgeColoring[CycleGraph[4]],
  {1, 2, 2, 1}
]


MT[
  IGEdgeColoring[IGTriangularLattice[4]],
  {4, 1, 2, 3, 1, 6, 5, 3, 4, 2, 1, 5, 3, 2, 4, 1, 2, 3}
]


(* ::Subsubsection:: *)
(*IGKVertexColoring*)


(* TODO *)


(* ::Subsubsection:: *)
(*IGKEdgeColoring*)


(* TODO *)


(* ::Subsubsection:: *)
(*IGMinimumVertexColoring*)


MT[
  AllTrue[
    GraphData[10],
    Max@IGMinimumVertexColoring@GraphData[#] == GraphData[#, "ChromaticNumber"]&
  ],
  True
]


MT[
  IGMinimumVertexColoring[IGEmptyGraph[#]],
  ConstantArray[1,#]
]& /@ Range[0, 4]


(* TODO test all heuristics *)


(* ::Subsubsection:: *)
(*IGMinimumEdgeColoring*)


MT[
  IGMinimumEdgeColoring[IGEmptyGraph[#]],
  {}
]& /@ Range[0, 4]


MT[
  IGMinimumEdgeColoring[CycleGraph[3]],
  {3, 2, 1}
]


(* Note: {1, 2, 2, 1} and {2, 1, 1, 2} are both valid. *)
MT[
  MemberQ[{{2, 1, 1, 2}, {1, 2, 2, 1}}, IGMinimumEdgeColoring[CycleGraph[4]]],
  True
]


MT[
  IGVertexColoringQ[LineGraph[IGTriangularLattice[4]], IGMinimumEdgeColoring[IGTriangularLattice[4]]],
  True
]

MT[
  Max[IGMinimumEdgeColoring[IGTriangularLattice[4]]],
  6
]


(* TODO test all heuristics *)


(* ::Subsubsection::Closed:: *)
(*IGChromaticNumber*)


MT[
  IGChromaticNumber /@ NestList[IGMycielskian, IGEmptyGraph[], 6],
  Range[0, 6]
]


MT[
  IGChromaticNumber[IGEmptyGraph[]],
  0
]


MT[
  IGChromaticNumber[IGEmptyGraph[#]],
  1
]& /@ Range[3]


(* directed edges are ignored *)
MT[
  IGChromaticNumber[Graph[{1, 2}, {1 -> 2}]],
  2
]


(* ::Subsubsection::Closed:: *)
(*IGChromaticIndex*)


MT[
  IGChromaticIndex[IGEmptyGraph[#]],
  0
]& /@ Range[0, 4]

MT[
  IGChromaticIndex[IGSquareLattice[{3, 4}]],
  4
]

MT[
  IGChromaticIndex[CycleGraph[5]],
  3
]

MT[
  IGChromaticIndex[CycleGraph[6]],
  2
]

MT[
  IGChromaticIndex[GraphData["PappusGraph"]],
  3
]


(* directed edges are ignored *)
MT[
  IGChromaticIndex[Graph[{1, 2}, {1 -> 2}]],
  1
]


(* ::Subsubsection:: *)
(*IGCliqueCover*)


(* TODO *)


(* ::Subsubsection:: *)
(*IGCliqueCoverNumber*)


(* TODO *)


(* ::Subsubsection:: *)
(*IGVertexColoringQ*)


(* TODO *)


(* ::Subsubsection:: *)
(*IGPerfectQ*)


MT[
  IGPerfectQ[IGEmptyGraph[#]],
  True
]& /@ Range[0,2]


(* this caused a crash earlier *)
MT[
  IGPerfectQ[GraphData[{"Arrangement", {5, 2}}]],
  True
]


(* TODO *)


(* ::Section::Closed:: *)
(*Proximity graphs*)


MTSection["Proximity graphs"]


(* ::Subsubsection::Closed:: *)
(*Common tests*)


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


(* ::Subsubsection::Closed:: *)
(*IGDelaunayGraph*)


MT[
  IGDelaunayGraph[{}],
  Graph[{}],
  SameTest -> IGSameGraphQ
]

MT[
  IGDelaunayGraph[{{0, 0}}],
  Graph[{1}, {}],
  SameTest -> IGSameGraphQ
]

MT[
  IGDelaunayGraph[{{0, 0, 0}}],
  Graph[{1}, {}],
  SameTest -> IGSameGraphQ
]

MT[
  IGDelaunayGraph[{{1}, {2}, {3}}],
  Graph[{1 <-> 2, 2 <-> 3}],
  SameTest -> IGSameGraphQ
]

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

points3d = RandomReal[1, {100, 3}];
MT[
  IGDelaunayGraph[points3d],
  IGMeshGraph@DelaunayMesh[points3d],
  SameTest -> IGSameGraphQ
]

MT[
  With[{g = IGDelaunayGraph[RandomReal[1, {10, 3}]]}, {VertexCount[g] == 10, ConnectedGraphQ[g]}],
  {True, True}
]


(* ::Subsubsection::Closed:: *)
(*IGLuneBetaSkeleton*)


(* ::Text:: *)
(*bsPoints will be used for testing several other proximity graph functions below.*)


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
  IGLuneBetaSkeleton[GraphEmbedding@IGTriangularLattice[4], 1.99],
  IGTriangularLattice[4],
  SameTest -> IGSameGraphQ
]


(* ::Subsubsection::Closed:: *)
(*IGGabrielGraph*)


MT[
  IGGabrielGraph[bsPoints],
  Graph@{10 <-> 16, 3 <-> 16, 3 <-> 10, 3 <-> 8, 1 <-> 3, 1 <-> 16, 11 <-> 12,
    8 <-> 9, 1 <-> 7, 3 <-> 5, 5 <-> 14, 4 <-> 13, 4 <-> 15, 9 <-> 15,
    17 <-> 18, 18 <-> 19, 3 <-> 20, 5 <-> 8, 1 <-> 8, 10 <-> 20,
    9 <-> 11, 2 <-> 6, 7 <-> 12, 5 <-> 20, 5 <-> 15, 4 <-> 17, 9 <-> 19,
    6 <-> 12, 2 <-> 11, 13 <-> 17, 4 <-> 19},
  SameTest -> IGSameGraphQ
]


circlePoints[n_] := Table[{Cos[x], Sin[x]}, {x, 0, 2Pi - 2Pi/n, 2Pi/n}]

MT[
  IGGabrielGraph[circlePoints[7]],
  CycleGraph[7],
  SameTest -> IGSameGraphQ
]


(* ::Subsubsection::Closed:: *)
(*IGRelativeNeighborhoodGraph*)


MT[
  IGRelativeNeighborhoodGraph[bsPoints],
  Graph@{10 <-> 16, 1 <-> 3, 1 <-> 16, 11 <-> 12, 8 <-> 9, 1 <-> 7, 5 <-> 14,
    4 <-> 15, 17 <-> 18, 18 <-> 19, 3 <-> 20, 1 <-> 8, 10 <-> 20,
    9 <-> 11, 2 <-> 6, 7 <-> 12, 5 <-> 20, 5 <-> 15, 9 <-> 19, 2 <-> 11,
    13 <-> 17, 4 <-> 19},
  SameTest -> IGSameGraphQ
]


MT[
  IGRelativeNeighborhoodGraph@GraphEmbedding@IGTriangularLattice[5],
  IGEmptyGraph@VertexCount@IGTriangularLattice[5],
  SameTest -> IGSameGraphQ
]


(* ::Subsubsection::Closed:: *)
(*IGCircleBetaSkeleton*)


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


(* ::Section::Closed:: *)
(*Degree sequences*)


MTSection["Degree sequences"]


(* ::Subsubsection::Closed:: *)
(*IGGraphicalQ*)


MT[
  IGGraphicalQ[{}],
  True
]

MT[
  IGGraphicalQ[{0, 0, 0}],
  True
]


MT[
  IGGraphicalQ[{}, {}],
  True
]

MT[
  IGGraphicalQ[{0}, {0}],
  True
]


(* ::Text:: *)
(*Non-matching in- and out-degree vector lengths:*)


MT[
  IGGraphicalQ[{0}, {}],
  $Failed,
  {IGraphM::error}
]

MT[
  IGGraphicalQ[{}, {0}],
  $Failed,
  {IGraphM::error}
]


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
  IGGraphicalQ[{1, 2, 3}],
  False
]

MT[
  IGGraphicalQ[{3, 2, 1}],
  False
]

MT[
  IGGraphicalQ[{3, 2, 1}, MultiEdges -> True],
  True
]

MT[
  IGGraphicalQ[{5, 2, 1}],
  False
]

MT[
  IGGraphicalQ[{5, 2, 1}, MultiEdges -> True],
  False
]

MT[
  IGGraphicalQ[{5, 2, 1}, MultiEdges -> True, SelfLoops -> True],
  True
]


MT[
  IGGraphicalQ[{4,4,4,1,1,1,1}],
  False
]


MT[
  IGGraphicalQ[{4,0,0,0,2,4,2}],
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
  IGGraphicalQ[{0}, {0}],
  True
]


MT[
  IGGraphicalQ[{1, 1, 1}],
  False
]


MT[
  IGGraphicalQ[{-1}],
  False
]

MT[
  IGGraphicalQ[{1}, {-1}],
  False
]


MT[
  IGGraphicalQ[{0}, {0, 0}],
  $Failed,
  {IGraphM::error}
]


MT[
  IGGraphicalQ[{}],
  True
]

MT[
  IGGraphicalQ[{}, {}],
  True
]


MT[
  IGGraphicalQ["not a degree sequence"],
  False
]


MT[
  IGGraphicalQ[{a,b}],
  False
]


MT[
  IGGraphicalQ[{0}],
  True
]


MT[
  IGGraphicalQ[{-1,1}],
  False
]


(* ::Subsubsection::Closed:: *)
(*IGBigraphicalQ*)


MT[
  IGBigraphicalQ[{}, {}],
  True
]

MT[
  IGBigraphicalQ[{}, {0, 0}],
  True
]

MT[
  IGBigraphicalQ[{}, {0, 0}, MultiEdges -> True],
  True
]


MT[
  IGBigraphicalQ @@ {{3, 3}, {1, 2, 3}},
  False
]

MT[
  IGBigraphicalQ @@ {{3, 1, 2}, {3, 3}},
  False
]

MT[
  IGBigraphicalQ @@ {{2, 2, 2}, {3, 3}},
  True
]

MT[
  IGBigraphicalQ @@ {{2, 2, 2}, {0, 0, 0, 3, 3}},
  True
]

MT[
  IGBigraphicalQ @@ {{2, 2, 2, 0, 0}, {3, 3}},
  True
]

MT[
  IGBigraphicalQ @@ {{1, 2, 3}, {3, 2, 1}},
  True
]

MT[
  IGBigraphicalQ @@ {{1, 2, 3}, {3, 1, 1}},
  False
]


MT[
  IGBigraphicalQ[{1, 2, 3}, {3, 1, 1}, MultiEdges -> True],
  False
]

MT[
  IGBigraphicalQ[{1, 2, 3}, {3, 3}, MultiEdges -> True],
  True
]


(* ::Subsubsection::Closed:: *)
(*IGRealizeDegreeSequence*)


(* ::Text:: *)
(*Check that degree sequences are preserved as expected.*)


Do[
  MT[
    VertexDegree@IGRealizeDegreeSequence[VertexDegree[g],Method->m],
    VertexDegree[g]
  ],
  {g,ulist}, {m,{"LargestFirst","SmallestFirst","Index"}}
]


Do[
  MT[
    Through[{VertexInDegree,VertexOutDegree}@IGRealizeDegreeSequence[VertexOutDegree[g],VertexInDegree[g],Method->m]],
    {VertexOutDegree[g],VertexInDegree[g]}
  ],
  {g,dlist}, {m,{"LargestFirst","SmallestFirst","Index"}}
]


(* ::Text:: *)
(*Non-graphical cases must fail.*)


MT[
  IGRealizeDegreeSequence[{3,2,1}, Method -> #],
  $Failed,
  {IGraphM::error}
]&/@{"SmallestFirst", "LargestFirst", "Index"}


MT[
  IGRealizeDegreeSequence[{4,4,4,1,1,1,1}, Method -> #],
  $Failed,
  {IGraphM::error}
]&/@{"SmallestFirst", "LargestFirst", "Index"}


MT[
  IGRealizeDegreeSequence[{4,0,0,0,2,4,2}, Method -> #],
  $Failed,
  {IGraphM::error}
]&/@{"SmallestFirst", "LargestFirst", "Index"}


MT[
  IGRealizeDegreeSequence[{}, Method -> #],
  IGEmptyGraph[],
  SameTest -> IGSameGraphQ
]&/@{"SmallestFirst", "LargestFirst", "Index"}


(* ::Text:: *)
(*Undirected examples.*)


MT[
  IGRealizeDegreeSequence[{3, 6, 3, 5, 5, 5, 2, 4, 4, 3}, Method -> "LargestFirst"],
  Graph[{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, {2 <-> 6, 2 <-> 5, 2 <-> 4, 2 <-> 9, 2 <-> 8, 2 <-> 10, 5 <-> 6, 4 <-> 6, 6 <-> 9, 6 <-> 8, 4 <-> 5, 3 <-> 5, 1 <-> 5, 3 <-> 4, 1 <-> 4, 8 <-> 9, 9 <-> 10, 7 <-> 8, 7 <-> 10, 1 <-> 3}],
  SameTest -> IGSameGraphQ
]

MT[
  IGRealizeDegreeSequence[{3, 6, 3, 5, 5, 5, 2, 4, 4, 3}, Method -> "SmallestFirst"],
  Graph[{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, {2 <-> 7, 4 <-> 7, 2 <-> 10, 5 <-> 10, 6 <-> 10, 2 <-> 3, 3 <-> 5, 3 <-> 6, 1 <-> 4, 1 <-> 8, 1 <-> 9, 4 <-> 6, 6 <-> 8, 6 <-> 9, 2 <-> 9, 5 <-> 9, 2 <-> 8, 5 <-> 8, 4 <-> 5, 2 <-> 4}],
  SameTest -> IGSameGraphQ
]

MT[
  IGRealizeDegreeSequence[{3, 6, 3, 5, 5, 5, 2, 4, 4, 3}, Method -> "Index"],
  Graph[{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, {1 <-> 2, 1 <-> 4, 1 <-> 5, 2 <-> 6, 2 <-> 4, 2 <-> 5, 2 <-> 8, 2 <-> 9, 3 <-> 6, 3 <-> 4, 3 <-> 5, 4 <-> 6, 4 <-> 8, 5 <-> 9, 5 <-> 10, 6 <-> 9, 6 <-> 10, 7 <-> 8, 7 <-> 9, 8 <-> 10}],
  SameTest -> IGSameGraphQ
]


MT[
  IGRealizeDegreeSequence[{2, 2, 2, 4, 2, 4, 3, 4, 2, 3}, Method -> "LargestFirst"],
  Graph[{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, {6 <-> 8, 4 <-> 8, 8 <-> 10, 7 <-> 8, 4 <-> 6, 6 <-> 10, 6 <-> 7, 4 <-> 9, 4 <-> 5, 2 <-> 3, 1 <-> 3, 1 <-> 2, 5 <-> 9, 7 <-> 10}],
  SameTest -> IGSameGraphQ
]

MT[
  IGRealizeDegreeSequence[{2, 2, 2, 4, 2, 4, 3, 4, 2, 3}, Method -> "SmallestFirst"],
  Graph[{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, {4 <-> 9, 6 <-> 9, 5 <-> 8, 4 <-> 5, 3 <-> 8, 3 <-> 6, 2 <-> 7, 2 <-> 10, 1 <-> 7, 1 <-> 10, 8 <-> 10, 6 <-> 7, 4 <-> 8, 4 <-> 6}],
  SameTest -> IGSameGraphQ
]

MT[
  IGRealizeDegreeSequence[{2, 2, 2, 4, 2, 4, 3, 4, 2, 3}, Method -> "Index"],
  Graph[{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, {1 <-> 4, 1 <-> 6, 2 <-> 8, 2 <-> 4, 3 <-> 8, 3 <-> 6, 4 <-> 7, 4 <-> 10, 5 <-> 7, 5 <-> 10, 6 <-> 8, 6 <-> 9, 7 <-> 8, 9 <-> 10}],
  SameTest -> IGSameGraphQ
]


MT[
  IGRealizeDegreeSequence[{2, 2, 2, 4, 2, 4, 3, 4, 2, 3}, Method -> "LargestFirst"],
  Graph[{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, {6 <-> 8, 4 <-> 8, 8 <-> 10, 7 <-> 8, 4 <-> 6, 6 <-> 10, 6 <-> 7, 4 <-> 9, 4 <-> 5, 2 <-> 3, 1 <-> 3, 1 <-> 2, 5 <-> 9, 7 <-> 10}],
  SameTest -> IGSameGraphQ
]

MT[
  IGRealizeDegreeSequence[{2, 2, 2, 4, 2, 4, 3, 4, 2, 3}, Method -> "SmallestFirst"],
  Graph[{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, {4 <-> 9, 6 <-> 9, 5 <-> 8, 4 <-> 5, 3 <-> 8, 3 <-> 6, 2 <-> 7, 2 <-> 10, 1 <-> 7, 1 <-> 10, 8 <-> 10, 6 <-> 7, 4 <-> 8, 4 <-> 6}],
  SameTest -> IGSameGraphQ
]

MT[
  IGRealizeDegreeSequence[{2, 2, 2, 4, 2, 4, 3, 4, 2, 3}, Method -> "Index"],
  Graph[{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, {1 <-> 4, 1 <-> 6, 2 <-> 8, 2 <-> 4, 3 <-> 8, 3 <-> 6, 4 <-> 7, 4 <-> 10, 5 <-> 7, 5 <-> 10, 6 <-> 8, 6 <-> 9, 7 <-> 8, 9 <-> 10}],
  SameTest -> IGSameGraphQ
]


(* ::Text:: *)
(*Directed examples.*)


MT[
  IGRealizeDegreeSequence[{2, 2, 1, 1, 2, 1, 1, 4, 1, 0}, {3, 1, 2, 0, 0, 2, 5, 1, 1, 0}, Method -> "LargestFirst"],
  Graph[{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, {8 -> 1, 2 -> 8, 7 -> 8, 7 -> 2, 7 -> 5, 7 -> 1, 7 -> 3, 6 -> 8, 6 -> 9, 1 -> 8, 1 -> 6, 1 -> 2, 3 -> 5, 3 -> 7, 9 -> 4}],
  SameTest -> IGSameGraphQ
]

MT[
  IGRealizeDegreeSequence[{2, 2, 1, 1, 2, 1, 1, 4, 1, 0}, {3, 1, 2, 0, 0, 2, 5, 1, 1, 0}, Method -> "SmallestFirst"],
  Graph[{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, {9 -> 8, 6 -> 8, 6 -> 1, 3 -> 8, 3 -> 2, 2 -> 5, 8 -> 7, 7 -> 1, 7 -> 8, 7 -> 5, 7 -> 2, 7 -> 3, 1 -> 6, 1 -> 9, 1 -> 4}],
  SameTest -> IGSameGraphQ
]

MT[
  IGRealizeDegreeSequence[{2, 2, 1, 1, 2, 1, 1, 4, 1, 0}, {3, 1, 2, 0, 0, 2, 5, 1, 1, 0}, Method -> "Index"],
  Graph[{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, {1 -> 8, 1 -> 2, 1 -> 5, 2 -> 8, 3 -> 8, 3 -> 1, 6 -> 7, 6 -> 8, 7 -> 9, 7 -> 6, 7 -> 1, 7 -> 3, 7 -> 2, 8 -> 5, 9 -> 4}],
  SameTest -> IGSameGraphQ
]


(* ::Subsubsection::Closed:: *)
(*IGSplitQ*)


MT[
  IGSplitQ[IGEmptyGraph[]],
  True
]


MT[
  AllTrue[GraphData /@ GraphData[;;3], IGSplitQ],
  True
]


MT[
  IGSplitQ[PathGraph[Range[4]]],
  True
]


MT[
  IGSplitQ[CycleGraph[4]],
  False
]


(* ::Subsubsection::Closed:: *)
(*IGThresholdQ*)


MT[
  IGThresholdQ[IGEmptyGraph[]],
  True
]


MT[
  AllTrue[GraphData /@ GraphData[;;3], IGThresholdQ],
  True
]


MT[
  IGThresholdQ[PathGraph[Range[4]]],
  False
]


MT[
  IGThresholdQ[CycleGraph[4]],
  False
]


(* ::Section::Closed:: *)
(*Neighbour degrees*)


MTSection["Neighbour degrees"]


(* ::Subsubsection::Closed:: *)
(*IGAverageNeighborDegree*)


MT[
  IGAverageNeighborDegree[IGEmptyGraph[]],
  {}
]

MT[
  IGAverageNeighborDegree[IGEmptyGraph[1]],
  {Indeterminate}
]

MT[
  IGAverageNeighborDegree[IGEmptyGraph[3]],
  {Indeterminate, Indeterminate, Indeterminate}
]


MT[
  IGAverageNeighborDegree[IGShorthand["1-2"]],
  {1., 1.}
]

MT[
  IGAverageNeighborDegree[IGShorthand["1,2-3"]],
  {Indeterminate, 1., 1.}
]


(* works only with simple graphs *)
MT[
  IGAverageNeighborDegree[Graph[{1 <-> 2, 1 <-> 2}]],
  $Failed,
  {IGraphM::error}
]


MT[
  IGAverageNeighborDegree[StarGraph[5]],
  {1., 4., 4., 4., 4.}
]

MT[
  IGAverageNeighborDegree[StarGraph[5], {}],
  {}
]

MT[
  IGAverageNeighborDegree[StarGraph[5], {1, 3}],
  {1., 4.}
]

MT[
  IGAverageNeighborDegree[StarGraph[5], All],
  {1., 4., 4., 4., 4.}
]


(* ::Text:: *)
(*Directed graphs:*)


MT[
  IGAverageNeighborDegree[Graph[{1 -> 2}]],
  {0., Indeterminate}
]

MT[
  IGAverageNeighborDegree[IGShorthand["1->2<->3, 1->3"]],
  {1., 1., 1.}
]

MT[
  IGAverageNeighborDegree[IGShorthand["1->2<->3, 1->3"], All, "All"],
  {3., 2.6666666666666665, 2.6666666666666665}
]

MT[
  IGAverageNeighborDegree[IGShorthand["1->2<->3, 1->3"], All, "In"],
  {Indeterminate, 1., 1.}
]

MT[
  IGAverageNeighborDegree[IGShorthand["1->2<->3, 1->3"], All, "Out"],
  {1., 1., 1.}
]

MT[
  IGAverageNeighborDegree[IGShorthand["1->2<->3, 1->3"], All, "Out", "All"],
  {1., 1.3333333333333333, 1.3333333333333333}
]

MT[
  IGAverageNeighborDegree[IGShorthand["1->2<->3, 1->3"], All, "Out", "In"],
  {Indeterminate, 1.5, 1.5}
]

MT[
  IGAverageNeighborDegree[IGShorthand["1->2<->3, 1->3"], All, "Out", Automatic],
  {1., 1., 1.}
]

MT[
  IGAverageNeighborDegree[IGShorthand["1->2<->3, 1->3"], All, "All", "Out"],
  {3., 3., 3.}
]

MT[
  IGAverageNeighborDegree[IGShorthand["1->2<->3, 1->3"], All, "All", "In"],
  {Indeterminate, 2.5, 2.5}
]

MT[
  IGAverageNeighborDegree[IGShorthand["1->2<->3, 1->3"], All, "In", "Out"],
  {2., 2., 2.}
]

MT[
  IGAverageNeighborDegree[IGShorthand["1->2<->3, 1->3"], All, "In", "All"],
  {2., 1.3333333333333333, 1.3333333333333333}
]


(* ::Subsubsection::Closed:: *)
(*IGAverageDegreeConnectivity*)


MT[
  IGAverageDegreeConnectivity[IGEmptyGraph[]],
  {}
]

MT[
  IGAverageDegreeConnectivity[IGEmptyGraph[5]],
  {}
]


MT[
  IGAverageDegreeConnectivity[StarGraph[5]],
  {4., Indeterminate, Indeterminate, 1.}
]

MT[
  IGAverageDegreeConnectivity[StarGraph[5, DirectedEdges -> True]],
  {Indeterminate, Indeterminate, Indeterminate, 0.}
]

MT[
  IGAverageDegreeConnectivity[StarGraph[5, DirectedEdges -> True], "Out"],
  {Indeterminate, Indeterminate, Indeterminate, 0.}
]

MT[
  IGAverageDegreeConnectivity[StarGraph[5, DirectedEdges -> True], "In"],
  {0.}
]

MT[
  IGAverageDegreeConnectivity[StarGraph[5, DirectedEdges -> True], "All"],
  {4., Indeterminate, Indeterminate, 1.}
]

MT[
  IGAverageDegreeConnectivity[StarGraph[5, DirectedEdges -> True], "All", "In"],
  {4.}
]

MT[
  IGAverageDegreeConnectivity[StarGraph[5, DirectedEdges -> True], "All", "Out"],
  {Indeterminate, Indeterminate, Indeterminate, 1.}
]

MT[
  IGAverageDegreeConnectivity[SimpleGraph[AdjacencyGraph[Outer[Boole @* Divisible, Range[10], Range[10]]]]],
  {0., 0.5, 0.7777777777777777}
]

MT[
  IGAverageDegreeConnectivity[IGKRegularGame[100, 3]],
  {Indeterminate, Indeterminate, 3.}
]


(* ::Section::Closed:: *)
(*Cycles*)


MTSection["Cycles"]


(* ::Subsection::Closed:: *)
(*Eulerian cycles*)


(* ::Subsubsection::Closed:: *)
(*IGEulerianQ*)


MT[
  IGEulerianQ[IGEmptyGraph[#]],
  True
]& /@ Range[0,4]

MT[
  IGEulerianQ[IGEmptyGraph[#], Closed -> True],
  True
]& /@ Range[0,4]


MT[
  IGEulerianQ[PathGraph[Range[5]]],
  True
]

MT[
  IGEulerianQ[PathGraph[Range[5]], Closed -> True],
  False
]


MT[
  IGEulerianQ[Graph[{1, 2}, {1 <-> 1}]],
  True
]

MT[
  IGEulerianQ[Graph[{1 <-> 1, 2 <-> 2}]],
  False
]

MT[
  IGEulerianQ[Graph[{1, 2}, {1 <-> 1, 1 <-> 1}]],
  True
]

MT[
  IGEulerianQ[Graph[{1 <-> 1, 2 <-> 2, 1 <-> 2}]],
  True
]

MT[
  IGEulerianQ[Graph[{1 <-> 1, 2 <-> 2, 1 <-> 2}], Closed -> True],
  False
]

MT[
  IGEulerianQ[Graph[{1 <-> 2, 1 <-> 2}], Closed -> True],
  True
]

MT[
  IGEulerianQ[Graph[{1 <-> 2, 1 <-> 2, 1 <-> 2}], Closed -> True],
  False
]

MT[
  IGEulerianQ[Graph[{1 <-> 2, 1 <-> 2, 1 <-> 2}], Closed -> False],
  True
]


MT[
  IGEulerianQ[Graph[{1 -> 2, 2 -> 3, 1 -> 3}]],
  False
]

MT[
  IGEulerianQ[Graph[{1 -> 2, 2 -> 3, 3 -> 1}]],
  True
]


(* ::Subsubsection::Closed:: *)
(*IGEulerianPath and IGEulerianPathVertices*)


MT[
  IGEulerianPath[Graph[{}]],
  {}
]

MT[
  IGEulerianPath[Graph[{1, 2}, {}]],
  {}
]

MT[
  IGEulerianPath[Graph[{1, 2}, {1 <-> 1}]],
  {UndirectedEdge[1, 1]}
]

MT[
  IGEulerianPath[CycleGraph[5]],
  {UndirectedEdge[1, 2], UndirectedEdge[2, 3], UndirectedEdge[3, 4], UndirectedEdge[4, 5], UndirectedEdge[1, 5]}
]

MT[
  IGEulerianPath[CycleGraph[5, DirectedEdges -> True]],
  {DirectedEdge[1, 2], DirectedEdge[2, 3], DirectedEdge[3, 4], DirectedEdge[4, 5], DirectedEdge[5, 1]}
]


MT[
  IGEulerianPathVertices[Graph[{}]],
  {}
]

MT[
  IGEulerianPathVertices[Graph[{1, 2, 3}, {}]],
  {}
]

MT[
  IGEulerianPathVertices[Graph[{1 <-> 2, 1 <-> 2}]],
  {1, 2, 1}
]

MT[
  IGEulerianPathVertices[Graph[{1 <-> 2, 1 <-> 2, 1 <-> 2}]],
  {1, 2, 1, 2}
]

MT[
  IGEulerianPathVertices[Graph[{1 <-> 1}]],
  {1, 1}
]

MT[
  IGEulerianPathVertices[Graph[{1 <-> 2, 1 <-> 3, 2 <-> 4, 1 <-> 4}]],
  {1, 2, 4, 1, 3}
]


(* ::Text:: *)
(*Non-Eulerian input:*)


MT[
  IGEulerianPath[Graph[{1 -> 2, 3 -> 2}]],
  $Failed,
  {IGraphM::error}
]

MT[
  IGEulerianPathVertices[Graph[{1 -> 2, 3 -> 2}]],
  $Failed,
  {IGraphM::error}
]


(* ::Section::Closed:: *)
(*Property operations*)


MTSection["Property operations"]


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


(* ::Subsection:: *)
(*Listing properties*)


proppedGraph = Graph[
	{Property[1,"Name"->"Foo"],Property[2,VertexCapacity->123],3},
	{Property[1<->2,"Bar"->44],2<->3},
	EdgeWeight->{1.23,2.34},
	Properties->{1->{"Type"->"xxx"}}
]


MT[
	Sort@IGVertexPropertyList[proppedGraph],
	{"Name","Type",VertexCapacity,VertexCoordinates,VertexShape,VertexShapeFunction,VertexSize,VertexStyle}
]


MT[
	Sort@IGEdgePropertyList[proppedGraph],
	{"Bar",EdgeShapeFunction,EdgeStyle,EdgeWeight}
]


(* ::Section::Closed:: *)
(*Matrix functions*)


MTSection["Matrix functions"]


(* ::Subsubsection::Closed:: *)
(*IGZeroDiagonal*)


mat= RandomInteger[2, {10,10}];

MT[
  Diagonal@IGZeroDiagonal[mat],
  ConstantArray[0, Length[mat]]
]

MT[
  Diagonal@IGZeroDiagonal@SparseArray[mat],
  SparseArray[{}, {Length[mat]}]
]

MT[
  Diagonal@IGZeroDiagonal[ mat[[All, 1;;5]] ],
  ConstantArray[0, 5]
]

MT[
  Diagonal@IGZeroDiagonal[ mat[[1;;5, All]] ],
  ConstantArray[0, 5]
]

(* UpperTriangularize crashes with non-rectangular sparse arrays in M \[LessEqual] 11.2 
   It is not worth working around this in IGraph/M at this moment, so we disable
   this test for affected versions. *)
If[$VersionNumber >= 11.2,
  MT[
    Diagonal@IGZeroDiagonal@SparseArray[ mat[[1;;5, All]] ],
    SparseArray[{}, {5}]
  ]
]

MT[
  IGZeroDiagonal[{{}}],
  {{}}
]


(* ::Subsubsection::Closed:: *)
(*IGTakeUpper and IGTakeLower*)


MT[
  IGTakeUpper[{{}}],
  {}
]

MT[
  IGTakeLower[{{}}],
  {}
]

MT[
  IGTakeLower[{{}, {}}],
  {}
]


MT[
  IGTakeLower[{}],
  $Failed,
  {IGTakeLower::arg}
]


MT[
  IGTakeLower[a, b],
  $Failed,
  {IGTakeLower::arg}
]

MT[
  IGTakeUpper[a, b],
  $Failed,
  {IGTakeUpper::arg}
]


MT[
  IGTakeUpper[Partition[Range[16], 4]],
  {2, 3, 4, 7, 8, 12}
]

MT[
  IGTakeLower[Partition[Range[16], 4]],
  {5, 9, 10, 13, 14, 15}
]


MT[
  IGTakeUpper[Partition[N[Range[16] + I], 4]],
  {2. + 1.*I, 3. + 1.*I, 4. + 1.*I, 7. + 1.*I, 8. + 1.*I, 12. + 1.*I}
]

MT[
  IGTakeLower[Partition[Range[16]*I, 4]],
  {5*I, 9*I, 10*I, 13*I, 14*I, 15*I}
]


MT[
  IGTakeUpper[{{1, 2, 3}, {4, 5, 6}}],
  {2, 3, 6}
]

MT[
  IGTakeLower[{{1, 2, 3}, {4, 5, 6}}],
  {4}
]


MT[
  IGTakeUpper[{{1, 2}, {3, 4}, {5, 6}, {7, 8}}],
  {2}
]

MT[
  IGTakeLower[{{1, 2}, {3, 4}, {5, 6}, {7, 8}}],
  {3, 5, 6, 7, 8}
]


MT[
  IGTakeUpper[N[{{1, 2}, {3, 4}, {5, 6}, {7, 8}}]],
  {2.}
]

MT[
  IGTakeUpper[{{1., 3., 5., 7.}, {2., 4., 6., 8.}}],
  {3., 5., 7., 6., 8.}
]


MT[
  IGTakeUpper[{{"a", "b"}, {"c", "d"}}],
  {"b"}
]

MT[
  IGTakeLower[{{"a", "b"}, {"c", "d"}}],
  {"c"}
]


(* TODO sparse arrays *)


(* ::Subsubsection::Closed:: *)
(*IGJointDegreeMatrix*)


(* null graph *)
MT[
  IGJointDegreeMatrix[IGEmptyGraph[]],
  {{}}
]


(* edgeless graph *)
MT[
  IGJointDegreeMatrix[IGEmptyGraph[10]],
  {{}}
]


(* this verifies that the diagonal is not double-counted *)
MT[
  Normal[IGJointDegreeMatrix[IGShorthand["1-2-3-4-2"]]],
  {{0, 0, 1}, {0, 1, 2}, {1, 2, 0}}
]

MT[
  Normal[IGJointDegreeMatrix[CycleGraph[12]]],
  {{0, 0}, {0, 12}}
]


MT[
  Normal[IGJointDegreeMatrix[Graph[{1 <-> 1}]]],
  {{0, 0}, {0, 1}}
]


MT[
  Normal[IGJointDegreeMatrix[Graph[{1 <-> 2, 1 <-> 2, 2 <-> 3}]]],
  {{0, 0, 1}, {0, 0, 2}, {1, 2, 0}}
]


MT[
  Normal[IGJointDegreeMatrix[Graph[{1 -> 2}]]],
  {{1}}
]


MT[
  Normal[IGJointDegreeMatrix[Graph[{1 -> 2, 1 -> 2}]]],
  {{0, 0}, {0, 2}}
]


MT[
  Normal[IGJointDegreeMatrix[Graph[{1 -> 2, 2 -> 1}]]],
  {{2}}
]


MT[
  Normal[IGJointDegreeMatrix[Graph[{1 -> 1}]]],
  {{1}}
]


MT[
  Normal[IGJointDegreeMatrix[Graph[{1 -> 2, 1 -> 2, 2 -> 2, 3 -> 2}]]],
  {{0, 0, 0, 2}, {0, 0, 0, 2}}
]


(* ::Text:: *)
(*Test normalization*)


MT[
  Normal[IGJointDegreeMatrix[Graph[{1 -> 2, 1 -> 2, 2 -> 2, 3 -> 2, 4 -> 3}], Normalized -> True]],
  {{1/5, 0, 0, 2/5}, {0, 0, 0, 2/5}}
]


MT[
  Normal[IGJointDegreeMatrix[IGShorthand["1-2-3-4-2"], Normalized -> True]],
  {{0, 0, 1/4}, {0, 1/4, 1/2}, {1/4, 1/2, 0}}
]


(* ::Text:: *)
(*Test cropping and padding*)


MT[
  Normal[IGJointDegreeMatrix[IGShorthand["1-2-3-4-2"], {2, 5}]],
  {{0, 0, 1, 0, 0}, {0, 1, 2, 0, 0}}
]


MT[
  Normal[IGJointDegreeMatrix[IGShorthand["1-2-3-4-2"], {2, 2}]],
  {{0, 0}, {0, 1}}
]


MT[
  Normal[IGJointDegreeMatrix[IGShorthand["1-2-3-4-2"], 2]],
  {{0, 0}, {0, 1}}
]


MT[
  Normal[IGJointDegreeMatrix[IGShorthand["1-2-3-4-2"], 7]],
  {{0, 0, 1, 0, 0, 0, 0}, {0, 1, 2, 0, 0, 0, 0}, {1, 2, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0}}
]


(* ::Section::Closed:: *)
(*Planar graphs*)


MTSection["Planar graphs"]


(* TODO *)


(* only works on simple graphs *)
canonicalEmbedding[emb_]:=
	Module[{cemb},
		cemb=RotateLeft[#,First@Ordering[#,1]-1]&/@KeySort[emb];
		If[OrderedQ[cemb[[1,{1,2}]]],
			cemb,
			Reverse/@cemb
		]
	]


(* ::Subsubsection::Closed:: *)
(*IGPlanarQ*)


(* ::Subsubsubsection:: *)
(*Graphs*)


(* ::Text:: *)
(*Note that PlanarGraphQ does not consider the null graph planar while IGPlanarQ does.*)


MT[
  IGPlanarQ[IGEmptyGraph[#]],
  True
]& /@ Range[0, 3]


MT[
  IGPlanarQ[IGCompleteGraph[#]],
  True
]& /@ Range[0,4]


Function[g,
  MT[
    IGPlanarQ[g],
    False
  ],
  HoldAll
] /@
Hold[
  CompleteGraph[5], CompleteGraph[7], CompleteGraph[{3,3}], 
  GraphData[{"NoncayleyTransitive",{20,5}}], GraphData[{"HamiltonLaceable",{8,5}}],
  PetersenGraph[],
  IGMeshCellAdjacencyGraph[IGLatticeMesh["PentagonType4"],1]
] // ReleaseHold


Function[g,
  MT[
    IGPlanarQ[g],
    True
  ],
  HoldAll
] /@ 
Hold[
  GraphData[{"Tadpole",{5,5}}],
  GridGraph[{3,4}],
  IGTriangularLattice[15],
  CycleGraph[7],
  IGTreeGame[123],
  GraphData[{"JohnsonSkeleton",64}],
  IGDelaunayGraph@RandomReal[{-1,1}, {100,2}],
  IGLuneBetaSkeleton[RandomReal[{-1,1}, {100,2}], 1],
  IGLuneBetaSkeleton[RandomReal[{-1,1}, {100,2}], 1.5],
  IGLuneBetaSkeleton[RandomReal[{-1,1}, {100,2}], 2],
  IGLuneBetaSkeleton[RandomReal[{-1,1}, {100,2}], 2.5],
  IGCircleBetaSkeleton[RandomReal[{-1,1}, {100,2}], 1.5],
  IGMeshGraph@IGLatticeMesh["CairoPentagonal"],
  IGMeshCellAdjacencyGraph[IGLatticeMesh["PentagonType3"],0],
  IGMeshCellAdjacencyGraph[IGLatticeMesh["PentagonType4"],2,1],
  IGMeshCellAdjacencyGraph[IGLatticeMesh["PentagonType4"],2],
  IGMeshCellAdjacencyGraph[IGLatticeMesh["Basketweave"],2],
  IGMeshCellAdjacencyGraph[IGLatticeMesh["PentagonType1"],0,2],
  GraphData[{"SierpinskiCarpet", 4}]
] // ReleaseHold


(* multigraph with self-loops *)
MT[
  IGPlanarQ[Graph[{1 <-> 2, 1 <-> 2, 2 <-> 2}]],
  True
]


(* directed graph *)
MT[
  IGPlanarQ[Graph[{1 -> 2, 2 -> 3}]],
  True
]


(* disconnected graph also containing an isolated vertex *)
MT[
  IGPlanarQ[Graph[{1, 2, 3, 4, 5, 6}, {1 -> 2, 2 -> 3, 4 -> 5}]],
  True
]


(* ::Subsubsubsection:: *)
(*Embeddings*)


MT[
  IGPlanarQ[Association[]],
  True
]


MT[
  IGPlanarQ[Association[1 -> {}]],
  True
]


MT[
  IGPlanarQ[Association[1 -> {2}]],
  False,
  {IGraphM::invemb}
]


MT[
  IGPlanarQ[Association[1 -> {2}, 2 -> {1}]],
  True
]


(* ::Text:: *)
(*A planar embedding of K4:*)


MT[
  IGPlanarQ[Association[1 -> {4, 2, 3}, 2 -> {3, 1, 4}, 3 -> {4, 1, 2}, 4 -> {2, 1, 3}]],
  True
]


(* ::Text:: *)
(*A non-planar embedding of K4:*)


MT[
  IGPlanarQ[Association[1 -> {4, 3, 2}, 2 -> {1, 4, 3}, 3 -> {4, 2, 1}, 4 -> {3, 2, 1}]],
  False
]


(* ::Subsubsection::Closed:: *)
(*IGMaximalPlanarQ*)


MT[
  IGMaximalPlanarQ[IGEmptyGraph[]],
  True
]


MT[
  IGMaximalPlanarQ[IGEmptyGraph[1]],
  True
]


MT[
  IGMaximalPlanarQ[IGEmptyGraph[2]],
  False
]


MT[
  IGMaximalPlanarQ[Graph[{1 <-> 2}]],
  True
]


MT[
  IGMaximalPlanarQ[Graph[{1 <-> 2, 2 <-> 3}]],
  False
]


MT[
  IGMaximalPlanarQ[CycleGraph[3]],
  True
]


MT[
  IGMaximalPlanarQ[CompleteGraph[4]],
  True
]


(* not planar *)
MT[
  IGMaximalPlanarQ[CompleteGraph[5]],
  False
]


MT[
  IGMaximalPlanarQ[CompleteGraph[{3, 2}]],
  False
]


(* not planar *)
MT[
  IGMaximalPlanarQ[CompleteGraph[{3, 3}]],
  False
]


(* disconnected with maximal planar components *)
MT[
  IGMaximalPlanarQ[IGDisjointUnion[{CycleGraph[3], CycleGraph[3]}]],
  False
]


(* ::Subsubsection::Closed:: *)
(*IGOuterplanarQ*)


(* ::Subsubsubsection:: *)
(*Graphs*)


MT[
  IGOuterplanarQ[IGEmptyGraph[#]],
  True
]& /@ Range[0,5]


MT[
  IGOuterplanarQ[CycleGraph[#]],
  True
]& /@ Range[1,6]


(* disconnected *)
MT[
  IGOuterplanarQ[IGDisjointUnion[{CycleGraph[3], CycleGraph[4]}]],
  True
]


(* disconnected non-outerplanar *)
MT[
  IGOuterplanarQ[IGDisjointUnion[{CompleteGraph[4], CycleGraph[5]}]],
  False
]


(* planar but not outerplanar *)
MT[
  IGOuterplanarQ[CompleteGraph[4]],
  False
]


(* not planar *)
MT[
  IGOuterplanarQ[CompleteGraph[5]],
  False
]


(* ::Subsubsubsection:: *)
(*Embeddings*)


MT[
  IGOuterplanarQ[Association[]],
  True
]


MT[
  IGOuterplanarQ[Association[1 -> {}]],
  True
]


MT[
  IGOuterplanarQ[Association[1 -> {2}]],
  False,
  {IGraphM::invemb}
]


MT[
  IGOuterplanarQ[Association[1 -> {2}, 2 -> {1}]],
  True
]


(* ::Text:: *)
(*Non-planar embedding of planar graph:*)


MT[
  IGOuterplanarQ[Association[1 -> {2, 4, 3}, 2 -> {1, 3}, 3 -> {2, 4, 1}, 4 -> {1, 3}]],
  False
]


(* ::Text:: *)
(*Planar embedding of the same:*)


MT[
  IGOuterplanarQ[Association[1 -> {2, 3, 4}, 2 -> {1, 3}, 3 -> {2, 4, 1}, 4 -> {3, 1}]],
  True
]


(* ::Subsubsection::Closed:: *)
(*IGOuterplanarEmbedding*)


MT[
  IGOuterplanarEmbedding[IGEmptyGraph[]],
  <||>
]

MT[
  IGOuterplanarEmbedding[IGEmptyGraph[2]],
  <|1 -> {}, 2 -> {}|>
]


MT[
  IGOuterplanarEmbedding[CompleteGraph[3]],
  <|1 -> {2, 3}, 2 -> {1, 3}, 3 -> {2, 1}|>
]


MT[
  IGOuterplanarEmbedding[CompleteGraph[4]],
  $Failed,
  {IGOuterplanarEmbedding::nopl}
]


MT[
  IGOuterplanarEmbedding[IGShorthand["1-2-3-4-1,1-3"]],
  <|1 -> {2, 3, 4}, 2 -> {1, 3}, 3 -> {2, 4, 1}, 4 -> {3, 1}|>
]


MT[
  IGOuterplanarEmbedding[IGDisjointUnion[{CycleGraph[3], CycleGraph[4]}]],
  <|{1, 1} -> {{1, 2}, {1, 3}}, {1, 2} -> {{1, 1}, {1, 3}}, {1, 3} -> {{1, 2}, {1, 1}}, {2, 1} -> {{2, 2}, {2, 4}}, {2, 2} -> {{2, 1}, {2, 3}}, {2, 3} -> {{2, 2}, {2, 4}}, {2, 4} -> {{2, 3}, {2, 1}}|>
]


(* ::Subsubsection::Closed:: *)
(*IGCoordinatesToEmbedding and IGEmbeddingToCoordinates*)


g   = IGDelaunayGraph@RandomReal[1, {10, 2}];
emb = IGCoordinatesToEmbedding[g];

MT[
  IGPlanarQ[emb],
  True
]

MT[
  IGOuterplanarQ[emb],
  False
]

(* Note: Graph[g, VertexCoordinates -> coords] fails in some versions (11.x), so we use SetProperty instead. *)
MT[
	canonicalEmbedding@IGCoordinatesToEmbedding@SetProperty[g, VertexCoordinates-> Thread[VertexList[g] -> IGEmbeddingToCoordinates[emb]]],
	canonicalEmbedding[emb]
]


(* ::Subsubsection::Closed:: *)
(*IGDualGraph*)


MT[
	IGDualGraph@CompleteGraph[4],
	CompleteGraph[4],
	SameTest -> IGIsomorphicQ
]

MT[
	IGDualGraph@Graph[{1<->2}],
	IGEmptyGraph[1],
	SameTest -> IGIsomorphicQ
]

MT[
	IGDualGraph@IGEmptyGraph[1],
	IGEmptyGraph[],
	SameTest -> IGIsomorphicQ
]

MT[
	IGDualGraph@IGEmptyGraph[],
	IGEmptyGraph[],
	SameTest -> IGIsomorphicQ
]

MT[
	IGDualGraph@CycleGraph[3],
	Graph[{1<->2}],
	SameTest -> IGIsomorphicQ
]


(* ::Section::Closed:: *)
(*Layout functions*)


MTSection["Layout functions"]

layoutFunctions = {
  IGLayoutBipartite, IGLayoutCircle, IGLayoutDavidsonHarel,
  If[$SystemID === "Linux-ARM", (* Raspberry Pi 1 does not have enough memory for DrL *)
    Unevaluated@Sequence[],
    Unevaluated@Sequence[IGLayoutDrL, IGLayoutDrL3D]
  ],
  IGLayoutFruchtermanReingold, IGLayoutFruchtermanReingold3D,
  IGLayoutGEM, IGLayoutGraphOpt, IGLayoutKamadaKawai, IGLayoutKamadaKawai3D,
  IGLayoutPlanar, IGLayoutRandom, IGLayoutReingoldTilford,
  IGLayoutReingoldTilfordCircular, IGLayoutSphere
};

(* Check that each layout works, and issues no errors. *)
MT[
  GraphQ[#[Graph[{1<->2}]]],
  True
]&/@ layoutFunctions
  
MT[
  ListQ@OptionValue[Options[#[Graph[{1<->2}]],VertexCoordinates],VertexCoordinates],
  True
]&/@ layoutFunctions


MT[
  GraphQ[#[IGEmptyGraph[1]]],
  True
]&/@ layoutFunctions
  
MT[
  ListQ@OptionValue[Options[#[IGEmptyGraph[1]],VertexCoordinates],VertexCoordinates],
  True
]&/@ layoutFunctions


MT[
  #[IGEmptyGraph[]],
  IGEmptyGraph[],
  SameTest -> IGSameGraphQ (* there may be options added such as GraphLayout -> {"Dimension" -> 3} *)
]&/@ layoutFunctions


(* IGLayoutTutte only works on 3-connected graphs, so we need a different test case for it. *)
MT[
  GraphQ@IGLayoutTutte[CompleteGraph[4]],
  True
]

MT[
  ListQ@OptionValue[Options[IGLayoutTutte[CompleteGraph[4]],VertexCoordinates],VertexCoordinates],
  True
]


(* TODO basic checks to ensure that layout functions work 
   and do not throw errors for reasonable edge cases *)


(* ::Section::Closed:: *)
(*Weighted graphs*)


MTSection["Weighted graphs"]


(* ::Subsubsection::Closed:: *)
(*IGEdgeWeightedQ, IGVertexWeightedQ*)


vwg = Graph[{1,2,3},{1<->2},VertexWeight->{1,2,3}];

MT[
  IGEdgeWeightedQ[#],
  False
]& /@ Hold[
  empty, edgeless, ugi, ugs, dgi, dgs, umulti, dmulti, umulti2, dmulti2,
  vwg
] // ReleaseHold

MT[
  IGVertexWeightedQ[#],
  False
]& /@ Hold[
  empty, edgeless, ugi, ugs, dgi, dgs, umulti, dmulti, umulti2, dmulti2,
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



(* ::Subsubsection::Closed:: *)
(*IGWeighedAdjacencyGraph*)


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


MT[
  IGWeightedAdjacencyGraph[{{0, 1, 2}, {1, 5, 0}, {2, 0, 1}}],
  Graph[{1 <-> 2, 1 <-> 3, 2 <-> 2, 3 <-> 3}],
  SameTest -> IGSameGraphQ
]

MT[
  IGWeightedAdjacencyGraph[{{0, 1, 2}, {1, 5, 0}, {2, 0, 1}}, 1],
  Graph[{1, 2, 3}, {1 <-> 1, 1 <-> 3, 2 <-> 2, 2 <-> 3}],
  SameTest -> IGSameGraphQ
]

MT[
  IGWeightedAdjacencyGraph[{{0, 1, 2}, {1, 5, 0}, {2, 0, 1}}, Infinity],
  Graph[{1 <-> 1, 1 <-> 2, 1 <-> 3, 2 <-> 2, 2 <-> 3, 3 <-> 3}],
  SameTest -> IGSameGraphQ
]


MT[
  IGWeightedAdjacencyGraph[{{0, 1, 0}, {2, 5, 0}, {3, 0, 1}}],
  Graph[{1 -> 2, 2 -> 1, 2 -> 2, 3 -> 1, 3 -> 3}],
  SameTest -> IGSameGraphQ
]

(* With DirectedEdges \[Rule] False, only the upper triangular part of the matrix is considered. *)
MT[
  IGWeightedAdjacencyGraph[{{0, 1, 0}, {2, 5, 0}, {3, 0, 1}}, DirectedEdges -> False],
  Graph[{1 <-> 2, 2 <-> 2, 3 <-> 3}],
  SameTest -> IGSameGraphQ
]


(* ::Subsubsection::Closed:: *)
(*IGWeigthedAdjacencyMatrix*)


MT[
  IGWeightedAdjacencyMatrix[Graph[{1 <-> 2, 2 <-> 3}]],
  SparseArray[Automatic, {3, 3}, 0, {1, {{0, 1, 3, 4}, {{2}, {1}, {3}, {2}}}, {1, 1, 1, 1}}]
]


MT[
  IGWeightedAdjacencyMatrix[Graph[{1 <-> 2, 2 <-> 3}, EdgeWeight->{xx, yy}]],
  SparseArray[Automatic, {3, 3}, 0, {1, {{0, 1, 3, 4}, {{2}, {1}, {3}, {2}}}, {xx, xx, yy, yy}}]
]


MT[
  IGWeightedAdjacencyMatrix[Graph[{1 <-> 2, 2 <-> 3}, EdgeWeight->{5., 6.}]],
  SparseArray[Automatic, {3, 3}, 0, {1, {{0, 1, 3, 4}, {{2}, {1}, {3}, {2}}}, {5., 5., 6., 6.}}]
]


MT[
  IGWeightedAdjacencyMatrix[IGEmptyGraph[1]],
  SparseArray[Automatic, {1, 1}, 0, {1, {{0, 0}, {}}, {}}]
]


(* ::Text:: *)
(*Verify that zero weights are preserved.*)


MT[
  IGWeightedAdjacencyMatrix[Graph[{1 <-> 2, 2 <-> 3}, EdgeWeight->{0., 6.}], Infinity],
  SparseArray[Automatic, {3, 3}, Infinity, {1, {{0, 1, 3, 4}, {{2}, {1}, {3}, {2}}}, {0., 0., 6., 6.}}]
]


MT[
  IGWeightedAdjacencyMatrix[Graph[{1 <-> 2, 2 <-> 3}, EdgeWeight->{0, 6}], Infinity],
  SparseArray[Automatic, {3, 3}, Infinity, {1, {{0, 1, 3, 4}, {{2}, {1}, {3}, {2}}}, {0, 0, 6, 6}}]
]


MT[
  IGWeightedAdjacencyMatrix[Graph[{1 <-> 2, 2 <-> 3}, EdgeWeight->{0, 6}], dummy],
  SparseArray[Automatic, {3, 3}, dummy, {1, {{0, 1, 3, 4}, {{2}, {1}, {3}, {2}}}, {0, 0, 6, 6}}]
]


(* ::Text:: *)
(*Verify that weights of parallel edges are added up.*)


MT[
  IGWeightedAdjacencyMatrix[Graph[{1 <-> 2, 1 <-> 2, 2 <-> 3}, EdgeWeight->{4,5,6}]],
  SparseArray[Automatic, {3, 3}, 0, {1, {{0, 1, 3, 4}, {{2}, {1}, {3}, {2}}}, {9, 9, 6, 6}}]
]


(* ::Text:: *)
(*Directed graph.*)


MT[
  IGWeightedAdjacencyMatrix[Graph[{1 -> 2, 2 -> 1, 2 -> 3}, EdgeWeight->{4,5,6}]],
  SparseArray[Automatic, {3, 3}, 0, {1, {{0, 1, 3, 3}, {{2}, {1}, {3}}}, {4, 5, 6}}]
]


(* ::Text:: *)
(*Mixed graph.*)


MT[
  IGWeightedAdjacencyMatrix[Graph[{1 -> 2, 2 -> 1, 2 <-> 3}, EdgeWeight->{4,5,6}]],
  SparseArray[Automatic, {3, 3}, 0, {1, {{0, 1, 3, 4}, {{2}, {1}, {3}, {2}}}, {4, 5, 6, 6}}]
]


(* ::Text:: *)
(*Null graph.*)


MT[
  IGWeightedAdjacencyMatrix[IGEmptyGraph[]],
  $Failed,
  {IGWeightedAdjacencyMatrix::nadj}
]


(* ::Subsubsection::Closed:: *)
(*IGUnweighted*)


(* ::Text:: *)
(*Do not use IGSameGraphQ. We must check that properties are preserved.*)


MT[
  IGUnweighted@Graph[{1 <-> 2}, EdgeWeight -> {2.5}, VertexCapacity -> {5,2}],
  Graph[{1 <-> 2}, VertexCapacity -> {5,2}],
  SameTest -> samePropGraphQ
]


MT[
  IGUnweighted@Graph[{1 <-> 2}, VertexWeight -> {3,4}],
  Graph[{1 <-> 2}, VertexWeight -> {3,4}],
  SameTest -> samePropGraphQ
]


MT[
  IGUnweighted@AdjacencyGraph[{{0,1},{1,0}}, EdgeWeight -> {2.5}, VertexCapacity -> {5,2}],
  AdjacencyGraph[{{0,1},{1,0}}, VertexCapacity -> {5,2}],
  SameTest -> samePropGraphQ
]


MT[
  IGUnweighted@IGEmptyGraph[],
  IGEmptyGraph[]
]


MT[
  IGUnweighted@IGEmptyGraph[1, VertexStyle -> {1 -> Red}],
  IGEmptyGraph[1, VertexStyle -> {1 -> Red}],
  SameTest -> samePropGraphQ
]


(* ::Subsubsection::Closed:: *)
(*Vertex strength*)


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


MT[
  #[IGEmptyGraph[]],
  {}
]& /@ {IGVertexStrength, IGVertexInStrength, IGVertexOutStrength}


(* ::Section::Closed:: *)
(*Visualization*)


MTSection["Visualization"]


(* ::Subsubsection::Closed:: *)
(*IGAdjacencyMatrixPlot*)


(* ::Text:: *)
(*Exercise function.*)


MT[
  Head[IGAdjacencyMatrixPlot[Graph[{1 <-> 2}]]],
  Graphics
]


(* ::Text:: *)
(*Test for invalid input.*)


MT[
  IGAdjacencyMatrixPlot[IGEmptyGraph[]],
  $Failed,
  {WeightedAdjacencyMatrix::nadj}
]


MT[
  IGAdjacencyMatrixPlot[Graph[{1 <-> 2}], EdgeWeight -> "foo"],
  $Failed,
  {IGAdjacencyMatrixPlot::noprop}
]


(* TODO *)


(* ::Section::Closed:: *)
(*Utility functions*)


MTSection["Utility functions"]


(* ::Subsubsection::Closed:: *)
(*IGUndirectedGraph*)


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


(* ::Text:: *)
(*Check for unchanged edge ordering with "All".*)


With[{ug = IGUndirectedGraph[#, "All"]},
  MT[
    {VertexList[ug], Sort /@ List @@@ EdgeList[ug]},
    {VertexList[#], Sort /@ List @@@ EdgeList[#]}
  ]
]& /@ Join[
  Flatten@Table[
    Graph[DirectedEdge @@@ RandomInteger[{1,k}, {m, 2}]],
    {k, 10, 100, 10},
    {m, k, 4k, k}
  ],
  {dmulti},
  dlist
] 


(* ::Subsubsection::Closed:: *)
(*IGSimpleGraph*)


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
    IGSimpleGraph[ug, MultiEdges -> True],
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

dg = Graph[{7, 9, 10, 5, 6, 0, 8, 2, 4, 3, 1}, {DirectedEdge[7, 9],
  DirectedEdge[10, 5], DirectedEdge[6, 5], DirectedEdge[10, 0],
  DirectedEdge[5, 8], DirectedEdge[9, 9], DirectedEdge[2, 4],
  DirectedEdge[9, 3], DirectedEdge[3, 1], DirectedEdge[6, 0],
  DirectedEdge[7, 3], DirectedEdge[0, 5], DirectedEdge[7, 3],
  DirectedEdge[0, 0], DirectedEdge[8, 6], DirectedEdge[8, 0],
  DirectedEdge[10, 7], DirectedEdge[9, 0], DirectedEdge[5, 3],
  DirectedEdge[3, 0], DirectedEdge[6, 1], DirectedEdge[3, 10],
  DirectedEdge[8, 8], DirectedEdge[7, 8], DirectedEdge[6, 8],
  DirectedEdge[4, 8], DirectedEdge[6, 8], DirectedEdge[4, 0],
  DirectedEdge[6, 3], DirectedEdge[3, 5], DirectedEdge[5, 0],
  DirectedEdge[6, 0], DirectedEdge[2, 1], DirectedEdge[5, 0],
  DirectedEdge[10, 0], DirectedEdge[0, 0], DirectedEdge[5, 6],
  DirectedEdge[1, 1], DirectedEdge[1, 2], DirectedEdge[3, 10],
  DirectedEdge[0, 0], DirectedEdge[10, 6], DirectedEdge[2, 1],
  DirectedEdge[10, 10], DirectedEdge[7, 10], DirectedEdge[2, 10],
  DirectedEdge[9, 5], DirectedEdge[1, 10], DirectedEdge[7, 3],
  DirectedEdge[0, 10]}];

MT[
  sameGraphQ[
    IGSimpleGraph[dg],
    Graph[{7, 9, 10, 5, 6, 0, 8, 2, 4, 3, 1}, {DirectedEdge[6, 1],
      DirectedEdge[2, 1], DirectedEdge[1, 2], DirectedEdge[10, 6],
      DirectedEdge[2, 10], DirectedEdge[1, 10], DirectedEdge[5, 6],
      DirectedEdge[10, 5], DirectedEdge[6, 5], DirectedEdge[8, 6],
      DirectedEdge[5, 8], DirectedEdge[6, 8], DirectedEdge[3, 1],
      DirectedEdge[3, 10], DirectedEdge[3, 5], DirectedEdge[5, 3],
      DirectedEdge[6, 3], DirectedEdge[0, 5], DirectedEdge[0, 10],
      DirectedEdge[10, 0], DirectedEdge[6, 0], DirectedEdge[8, 0],
      DirectedEdge[3, 0], DirectedEdge[5, 0], DirectedEdge[4, 8],
      DirectedEdge[4, 0], DirectedEdge[2, 4], DirectedEdge[9, 3],
      DirectedEdge[9, 0], DirectedEdge[9, 5], DirectedEdge[7, 9],
      DirectedEdge[7, 3], DirectedEdge[7, 8], DirectedEdge[7, 10],
      DirectedEdge[10, 7]}]
  ],
  True
]

MT[
  sameGraphQ[
    IGSimpleGraph[dg, SelfLoops -> True],
    Graph[{7, 9, 10, 5, 6, 0, 8, 2, 4, 3, 1}, {DirectedEdge[7, 9],
    DirectedEdge[7, 10], DirectedEdge[7, 8], DirectedEdge[7, 3],
    DirectedEdge[9, 9], DirectedEdge[9, 5], DirectedEdge[9, 0],
    DirectedEdge[9, 3], DirectedEdge[10, 7], DirectedEdge[10, 10],
    DirectedEdge[10, 5], DirectedEdge[10, 6], DirectedEdge[10, 0],
    DirectedEdge[5, 6], DirectedEdge[5, 0], DirectedEdge[5, 8],
    DirectedEdge[5, 3], DirectedEdge[6, 5], DirectedEdge[6, 0],
    DirectedEdge[6, 8], DirectedEdge[6, 3], DirectedEdge[6, 1],
    DirectedEdge[0, 10], DirectedEdge[0, 5], DirectedEdge[0, 0],
    DirectedEdge[8, 6], DirectedEdge[8, 0], DirectedEdge[8, 8],
    DirectedEdge[2, 10], DirectedEdge[2, 4], DirectedEdge[2, 1],
    DirectedEdge[4, 0], DirectedEdge[4, 8], DirectedEdge[3, 10],
    DirectedEdge[3, 5], DirectedEdge[3, 0], DirectedEdge[3, 1],
    DirectedEdge[1, 10], DirectedEdge[1, 2], DirectedEdge[1, 1]}]
  ],
  True
]

MT[
  sameGraphQ[
    IGSimpleGraph[dg, MultiEdges -> True],
    Graph[{7, 9, 10, 5, 6, 0, 8, 2, 4, 3, 1}, {DirectedEdge[7, 9],
    DirectedEdge[7, 10], DirectedEdge[7, 8], DirectedEdge[7, 3],
    DirectedEdge[7, 3], DirectedEdge[7, 3], DirectedEdge[9, 5],
    DirectedEdge[9, 0], DirectedEdge[9, 3], DirectedEdge[10, 7],
    DirectedEdge[10, 5], DirectedEdge[10, 6], DirectedEdge[10, 0],
    DirectedEdge[10, 0], DirectedEdge[5, 6], DirectedEdge[5, 0],
    DirectedEdge[5, 0], DirectedEdge[5, 8], DirectedEdge[5, 3],
    DirectedEdge[6, 5], DirectedEdge[6, 0], DirectedEdge[6, 0],
    DirectedEdge[6, 8], DirectedEdge[6, 8], DirectedEdge[6, 3],
    DirectedEdge[6, 1], DirectedEdge[0, 10], DirectedEdge[0, 5],
    DirectedEdge[8, 6], DirectedEdge[8, 0], DirectedEdge[2, 10],
    DirectedEdge[2, 4], DirectedEdge[2, 1], DirectedEdge[2, 1],
    DirectedEdge[4, 0], DirectedEdge[4, 8], DirectedEdge[3, 10],
    DirectedEdge[3, 10], DirectedEdge[3, 5], DirectedEdge[3, 0],
    DirectedEdge[3, 1], DirectedEdge[1, 10], DirectedEdge[1, 2]}]
  ],
  True
]



(* ::Subsubsection::Closed:: *)
(*IGReorderVertices*)


MT[
  IGReorderVertices[{}, IGEmptyGraph[]],
  IGEmptyGraph[]
]

MT[
  IGReorderVertices[{1}, IGEmptyGraph[]],
  $Failed,
  {IGReorderVertices::bdvert}
]

MT[
  IGReorderVertices[{"Awad", "Odeh", "Osama", "Abouhlaima", "Fahad",
    "Fadhil", "Ghailani", "Kherchtou", "Ali", "Hage", "Abdullah",
    "Fazul", "Salim", "Khalfan", "Owhali", "Fawwaz", "Azzam",
    "Atwah"}, terrorist]
  ,
  Graph[{"Awad", "Odeh", "Osama", "Abouhlaima", "Fahad", "Fadhil", "Ghailani",
    "Kherchtou", "Ali", "Hage", "Abdullah", "Fazul", "Salim", "Khalfan", "Owhali",
    "Fawwaz", "Azzam", "Atwah"}, {UndirectedEdge["Osama", "Salim"],
    UndirectedEdge["Osama", "Ali"], UndirectedEdge["Osama", "Abdullah"],
    UndirectedEdge["Osama", "Hage"], UndirectedEdge["Osama", "Owhali"],
    UndirectedEdge["Salim", "Abdullah"], UndirectedEdge["Salim", "Hage"],
    UndirectedEdge["Ali", "Abouhlaima"], UndirectedEdge["Ali", "Kherchtou"],
    UndirectedEdge["Ali", "Hage"], UndirectedEdge["Ali", "Fazul"],
    UndirectedEdge["Abouhlaima", "Hage"], UndirectedEdge["Kherchtou", "Fawwaz"],
    UndirectedEdge["Kherchtou", "Hage"], UndirectedEdge["Kherchtou", "Odeh"],
    UndirectedEdge["Kherchtou", "Atwah"], UndirectedEdge["Fawwaz", "Hage"],
    UndirectedEdge["Fawwaz", "Owhali"], UndirectedEdge["Fawwaz", "Azzam"],
    UndirectedEdge["Abdullah", "Odeh"], UndirectedEdge["Abdullah", "Owhali"],
    UndirectedEdge["Abdullah", "Fazul"], UndirectedEdge["Abdullah", "Atwah"],
    UndirectedEdge["Abdullah", "Fahad"], UndirectedEdge["Hage", "Odeh"],
    UndirectedEdge["Hage", "Fazul"], UndirectedEdge["Hage", "Fadhil"],
    UndirectedEdge["Odeh", "Owhali"], UndirectedEdge["Odeh", "Fazul"],
    UndirectedEdge["Odeh", "Azzam"], UndirectedEdge["Odeh", "Atwah"],
    UndirectedEdge["Odeh", "Fadhil"], UndirectedEdge["Owhali", "Fazul"],
    UndirectedEdge["Owhali", "Azzam"], UndirectedEdge["Owhali", "Atwah"],
    UndirectedEdge["Fazul", "Azzam"], UndirectedEdge["Fazul", "Atwah"],
    UndirectedEdge["Azzam", "Atwah"], UndirectedEdge["Atwah", "Fahad"],
    UndirectedEdge["Atwah", "Fadhil"], UndirectedEdge["Atwah", "Awad"],
    UndirectedEdge["Fahad", "Fadhil"], UndirectedEdge["Fahad", "Khalfan"],
    UndirectedEdge["Fahad", "Ghailani"], UndirectedEdge["Fahad", "Awad"],
    UndirectedEdge["Fadhil", "Khalfan"], UndirectedEdge["Fadhil", "Ghailani"],
    UndirectedEdge["Fadhil", "Awad"], UndirectedEdge["Khalfan", "Ghailani"],
    UndirectedEdge["Khalfan", "Awad"], UndirectedEdge["Ghailani", "Awad"]},
   {Properties -> {"Ghailani" -> {"FullName" -> "Ahmed Khalfan Ghailani",
        "Group" -> "Dar es Salaam Cell"}, "Osama" -> {"FullName" -> "Osama Bin Laden",
        "Group" -> "Planners"}, "Fazul" -> {"FullName" -> "Fazul Abdullah Mohammed",
        "Group" -> "Nairobi Cell"}, "Fadhil" -> {"FullName" -> "Mustafa Mohammed Fadhil",
        "Group" -> "Dar es Salaam Cell"},
      "Atwah" -> {"FullName" -> "Muhsin Musa Matwalli Atwah", "Group" -> "Planners"},
      "Ali" -> {"FullName" -> "Ali Mohamed", "Group" -> "Planners"},
      "Fawwaz" -> {"FullName" -> "Khalid al-Fawwaz", "Group" -> "Planners"},
      "Fahad" -> {"FullName" -> "Fahad Mohammed Ally Msalam",
        "Group" -> "Dar es Salaam Cell"}, "Odeh" -> {"FullName" -> "Mohamed Sadeek Odeh",
        "Group" -> "Nairobi Cell"}, "Khalfan" -> {"FullName" -> "Khalfan Khamis Mohamed",
        "Group" -> "Dar es Salaam Cell"}, "Hage" -> {"FullName" -> "Wadih el-Hage",
        "Group" -> "Planners"}, "Owhali" ->
       {"FullName" -> "Mohamed Rashed Daoud al-Owhali", "Group" -> "Nairobi Cell"},
      "Salim" -> {"FullName" -> "Mamdouh Salim", "Group" -> "Planners"},
      "Abdullah" -> {"FullName" -> "Abdullah Ahmed Abdullah", "Group" -> "Planners"},
      "Abouhlaima" -> {"FullName" -> "Abouhlaima", "Group" -> "Planners"},
      "Kherchtou" -> {"FullName" -> "Kherchtou", "Group" -> "Planners"},
      "Awad" -> {"FullName" -> "Hamden Khalif Allah Awad",
        "Group" -> "Dar es Salaam Cell"}, "Azzam" -> {"FullName" -> "Azzam",
        "Group" -> "Nairobi Cell"}}, EdgeStyle ->
     {UndirectedEdge["Khalfan", "Awad"] -> Thickness[0.0012],
      UndirectedEdge["Fadhil", "Awad"] -> Thickness[0.0048],
      UndirectedEdge["Osama", "Ali"] -> Thickness[0.0036],
      UndirectedEdge["Abouhlaima", "Hage"] -> Thickness[0.0016],
      UndirectedEdge["Abdullah", "Fahad"] -> Thickness[0.0048],
      UndirectedEdge["Abdullah", "Odeh"] -> Thickness[0.0048],
      UndirectedEdge["Ghailani", "Awad"] -> Thickness[0.0012],
      UndirectedEdge["Salim", "Abdullah"] -> Thickness[0.0036],
      UndirectedEdge["Salim", "Hage"] -> Thickness[0.0036],
      UndirectedEdge["Fahad", "Ghailani"] -> Thickness[0.0048],
      UndirectedEdge["Odeh", "Azzam"] -> Thickness[0.0048],
      UndirectedEdge["Abdullah", "Owhali"] -> Thickness[0.0048],
      UndirectedEdge["Fahad", "Khalfan"] -> Thickness[0.0064],
      UndirectedEdge["Osama", "Hage"] -> Thickness[0.0048],
      UndirectedEdge["Atwah", "Awad"] -> Thickness[0.0012],
      UndirectedEdge["Kherchtou", "Hage"] -> Thickness[0.0036],
      UndirectedEdge["Abdullah", "Fazul"] -> Thickness[0.0048],
      UndirectedEdge["Osama", "Abdullah"] -> Thickness[0.0048],
      UndirectedEdge["Fawwaz", "Hage"] -> Thickness[0.0048],
      UndirectedEdge["Osama", "Owhali"] -> Thickness[0.0036],
      UndirectedEdge["Fadhil", "Khalfan"] -> Thickness[0.0064],
      UndirectedEdge["Ali", "Hage"] -> Thickness[0.0016],
      UndirectedEdge["Ali", "Abouhlaima"] -> Thickness[0.0036],
      UndirectedEdge["Fazul", "Azzam"] -> Thickness[0.0048],
      UndirectedEdge["Odeh", "Atwah"] -> Thickness[0.0048],
      UndirectedEdge["Osama", "Salim"] -> Thickness[0.005200000000000001],
      UndirectedEdge["Fawwaz", "Azzam"] -> Thickness[0.0036],
      UndirectedEdge["Kherchtou", "Atwah"] -> Thickness[0.0036],
      UndirectedEdge["Atwah", "Fadhil"] -> Thickness[0.0012],
      UndirectedEdge["Azzam", "Atwah"] -> Thickness[0.0012],
      UndirectedEdge["Ali", "Fazul"] -> Thickness[0.0072],
      UndirectedEdge["Atwah", "Fahad"] -> Thickness[0.0012],
      UndirectedEdge["Fadhil", "Ghailani"] -> Thickness[0.0048],
      UndirectedEdge["Fawwaz", "Owhali"] -> Thickness[0.0036],
      UndirectedEdge["Hage", "Fadhil"] -> Thickness[0.0048],
      UndirectedEdge["Kherchtou", "Odeh"] -> Thickness[0.0036],
      UndirectedEdge["Kherchtou", "Fawwaz"] -> Thickness[0.0036],
      UndirectedEdge["Hage", "Odeh"] -> Thickness[0.0064],
      UndirectedEdge["Odeh", "Fadhil"] -> Thickness[0.0048],
      UndirectedEdge["Owhali", "Atwah"] -> Thickness[0.0012],
      UndirectedEdge["Fahad", "Awad"] -> Thickness[0.0012],
      UndirectedEdge["Owhali", "Fazul"] -> Thickness[0.0064],
      UndirectedEdge["Ali", "Kherchtou"] -> Thickness[0.0036],
      UndirectedEdge["Hage", "Fazul"] -> Thickness[0.0048],
      UndirectedEdge["Odeh", "Fazul"] -> Thickness[0.0064],
      UndirectedEdge["Fazul", "Atwah"] -> Thickness[0.0012],
      UndirectedEdge["Khalfan", "Ghailani"] -> Thickness[0.0012],
      UndirectedEdge["Abdullah", "Atwah"] -> Thickness[0.0048],
      UndirectedEdge["Owhali", "Azzam"] -> Thickness[0.0028000000000000004],
      UndirectedEdge["Fahad", "Fadhil"] -> Thickness[0.0064],
      UndirectedEdge["Odeh", "Owhali"] -> Thickness[0.0048]},
    VertexLabels -> {Placed["Name", Above], "Abdullah" -> Placed["Name", Before],
      "Azzam" -> Placed["Name", Below], "Kherchtou" -> Placed["Name", After],
      "Hage" -> Placed["Name", After], "Fazul" -> Placed["Name", Below],
      "Abouhlaima" -> Placed["Name", After], "Awad" -> Placed["Name", Below],
      "Atwah" -> Placed["Name", Below], "Fahad" -> Placed["Name", Before],
      "Ali" -> Placed["Name", Above], "Fadhil" -> Placed["Name", After],
      "Khalfan" -> Placed["Name", Before], "Odeh" -> Placed["Name", Before],
      "Osama" -> Placed["Name", Before], "Ghailani" -> Placed["Name", After],
      "Owhali" -> Placed["Name", After]}, VertexStyle ->
     {"Odeh" -> RGBColor[0.6, 0.4, 0.2], "Fazul" -> RGBColor[0.6, 0.4, 0.2],
      "Awad" -> RGBColor[0.5, 0, 0.5], "Atwah" -> GrayLevel[0], "Salim" -> GrayLevel[0],
      "Fawwaz" -> GrayLevel[0], "Osama" -> GrayLevel[0], "Ali" -> GrayLevel[0],
      "Khalfan" -> RGBColor[0.5, 0, 0.5], "Fahad" -> RGBColor[0.5, 0, 0.5],
      "Azzam" -> RGBColor[0.6, 0.4, 0.2], "Kherchtou" -> GrayLevel[0],
      "Ghailani" -> RGBColor[0.5, 0, 0.5], "Hage" -> GrayLevel[0],
      "Abdullah" -> GrayLevel[0], "Fadhil" -> RGBColor[0.5, 0, 0.5],
      "Owhali" -> RGBColor[0.6, 0.4, 0.2], "Abouhlaima" -> GrayLevel[0]},
    EdgeWeight -> {0.52, 0.36, 0.48, 0.48, 0.36, 0.36, 0.36, 0.36, 0.36, 0.16, 0.72, 0.16,
     0.36, 0.36, 0.36, 0.36, 0.48, 0.36, 0.36, 0.48, 0.48, 0.48, 0.48, 0.48, 0.64, 0.48,
     0.48, 0.48, 0.64, 0.48, 0.48, 0.48, 0.64, 0.28, 0.12, 0.48, 0.12, 0.12, 0.12, 0.12,
     0.12, 0.64, 0.64, 0.48, 0.12, 0.64, 0.48, 0.48, 0.12, 0.12, 0.12},
    VertexCoordinates -> {{0.5760000000000001, 0.059000000000000004}, {0.2285, 0.209},
     {0.2435, 0.5650000000000001}, {0.6645, 0.357}, {0.583, 0.194}, {0.704, 0.234},
     {0.7015, 0.016}, {0.5135, 0.293}, {0.41200000000000003, 0.41500000000000004},
     {0.5955, 0.424}, {0.2925, 0.461}, {0.2405, 0.049}, {0.5740000000000001, 0.513},
     {0.8170000000000001, 0.111}, {0.0115, 0.272}, {0.216, 0.33}, {0.046, 0.115}, {0.423,
     0.069}}}],
  SameTest -> samePropGraphQ
]


(* ::Subsubsection::Closed:: *)
(*IGTakeSubgraph*)


(* just edges, non-induced subgraph result expected *)
MT[
  IGTakeSubgraph[terrorist, {"Azzam" <-> "Fazul", "Fazul" <-> "Atwah"}],
  Graph[{"Azzam", "Fazul", "Atwah"}, {UndirectedEdge["Azzam", "Fazul"],
    UndirectedEdge["Atwah", "Fazul"]},
    {Properties -> {"Fazul" -> {"FullName" -> "Fazul Abdullah Mohammed",
      "Group" -> "Nairobi Cell"}, "Atwah" ->
        {"FullName" -> "Muhsin Musa Matwalli Atwah", "Group" -> "Planners"},
      "Azzam" -> {"FullName" -> "Azzam", "Group" -> "Nairobi Cell"}},
      EdgeStyle -> {UndirectedEdge["Azzam", "Fazul"] -> Thickness[0.0048],
        UndirectedEdge["Atwah", "Fazul"] -> Thickness[0.0012]},
      VertexLabels -> {Placed["Name", Above], "Fazul" -> Placed["Name", Below],
        "Azzam" -> Placed["Name", Below], "Atwah" -> Placed["Name", Below]},
      VertexStyle -> {"Atwah" -> GrayLevel[0], "Fazul" -> RGBColor[0.6, 0.4, 0.2],
        "Azzam" -> RGBColor[0.6, 0.4, 0.2]}, EdgeWeight -> {0.48, 0.12},
      VertexCoordinates -> {{0.046, 0.115}, {0.2405, 0.049}, {0.423, 0.069}}}],
  SameTest -> samePropGraphQ
]


(* subgraph *)
MT[
  IGTakeSubgraph[terrorist, Subgraph[terrorist, {"Owhali", "Fawwaz", "Fazul", "Odeh", "Kherchtou", "Ghailani"}]],
  Graph[{"Owhali", "Fawwaz", "Fazul", "Odeh", "Kherchtou", "Ghailani"},
    {UndirectedEdge["Fawwaz", "Kherchtou"], UndirectedEdge["Kherchtou", "Odeh"],
      UndirectedEdge["Fawwaz", "Owhali"], UndirectedEdge["Odeh", "Owhali"],
      UndirectedEdge["Fazul", "Odeh"], UndirectedEdge["Fazul", "Owhali"]},
    {Properties -> {"Ghailani" -> {"FullName" -> "Ahmed Khalfan Ghailani",
      "Group" -> "Dar es Salaam Cell"},
      "Fazul" -> {"FullName" -> "Fazul Abdullah Mohammed", "Group" -> "Nairobi Cell"},
      "Fawwaz" -> {"FullName" -> "Khalid al-Fawwaz", "Group" -> "Planners"},
      "Odeh" -> {"FullName" -> "Mohamed Sadeek Odeh", "Group" -> "Nairobi Cell"},
      "Owhali" -> {"FullName" -> "Mohamed Rashed Daoud al-Owhali",
        "Group" -> "Nairobi Cell"}, "Kherchtou" -> {"FullName" -> "Kherchtou",
        "Group" -> "Planners"}}, EdgeStyle ->
        {UndirectedEdge["Fazul", "Odeh"] -> Thickness[0.0064],
          UndirectedEdge["Kherchtou", "Odeh"] -> Thickness[0.0036],
          UndirectedEdge["Fawwaz", "Kherchtou"] -> Thickness[0.0036],
          UndirectedEdge["Fazul", "Owhali"] -> Thickness[0.0064],
          UndirectedEdge["Fawwaz", "Owhali"] -> Thickness[0.0036],
          UndirectedEdge["Odeh", "Owhali"] -> Thickness[0.0048]},
      VertexLabels -> {Placed["Name", Above], "Ghailani" -> Placed["Name", After],
        "Fazul" -> Placed["Name", Below], "Owhali" -> Placed["Name", After],
        "Kherchtou" -> Placed["Name", After], "Odeh" -> Placed["Name", Before]},
      VertexStyle -> {"Ghailani" -> RGBColor[0.5, 0, 0.5], "Fawwaz" -> GrayLevel[0],
        "Kherchtou" -> GrayLevel[0], "Fazul" -> RGBColor[0.6, 0.4, 0.2],
        "Odeh" -> RGBColor[0.6, 0.4, 0.2], "Owhali" -> RGBColor[0.6, 0.4, 0.2]},
      EdgeWeight -> {0.36, 0.36, 0.36, 0.48, 0.64, 0.64}, VertexCoordinates -> {{0.0115,
      0.272}, {0.216, 0.33}, {0.2405, 0.049}, {0.2285, 0.209}, {0.5135, 0.293}, {0.7015,
      0.016}}}],
  SameTest -> samePropGraphQ
]


(* check weights *)
MT[
  WeightedAdjacencyMatrix@IGTakeSubgraph[terrorist, Subgraph[terrorist, {"Owhali", "Fawwaz", "Fazul", "Odeh", "Kherchtou", "Ghailani"}]],
  WeightedAdjacencyMatrix@IGWeightedSubgraph[terrorist, {"Owhali", "Fawwaz", "Fazul", "Odeh", "Kherchtou", "Ghailani"}]
]


(* empty graph behaviour *)
MT[
  IGTakeSubgraph[terrorist, {}],
  Graph[{}, {}]
]

MT[
  IGTakeSubgraph[terrorist, IGEmptyGraph[]],
  Graph[{}, {}]
]

MT[
  IGTakeSubgraph[IGEmptyGraph[],IGEmptyGraph[]],
  IGEmptyGraph[]
]


(* check for bad input *)
MT[
  IGTakeSubgraph[Graph[{1 <-> 2, 2 <-> 3}], Graph[{3 <-> 4}]],
  $Failed,
  {IGTakeSubgraph::nsg}
]

MT[
  IGTakeSubgraph[Graph[{1 <-> 2, 2 <-> 3}], {3 <-> 3}],
  $Failed,
  {IGTakeSubgraph::nsg}
]

MT[
  IGTakeSubgraph[Graph[{1 <-> 2, 2 <-> 3}], Graph[{3, 4, 2}, {2 <-> 3}]],
  $Failed,
  {IGTakeSubgraph::nsg}
]


(* check TwoWayRule handling *)
If[$VersionNumber >= 11.2,
  MT[
    IGTakeSubgraph[CycleGraph[5], {TwoWayRule[1,2],TwoWayRule[2,3]}],
    Graph[{1<->2,2<->3}],
    SameTest->IGSameGraphQ
  ];

  MT[
    IGTakeSubgraph[CycleGraph[5], {UndirectedEdge[1,2],TwoWayRule[2,3]}],
    Graph[{1<->2,2<->3}],
    SameTest->IGSameGraphQ
  ]
]

MT[
  IGTakeSubgraph[CycleGraph[5], {UndirectedEdge[1,2],UndirectedEdge[2,3]}],
  Graph[{1<->2,2<->3}],
  SameTest->IGSameGraphQ
]


(* ::Subsubsection::Closed:: *)
(*IGDisjointUnion*)


MT[
  IGDisjointUnion[{}],
  Graph[{}],
  SameTest -> IGSameGraphQ
]


MT[
  IGDisjointUnion[Association[]],
  Graph[{}],
  SameTest -> IGSameGraphQ
]


MT[
  IGDisjointUnion[{IGShorthand["1-2-3"], IGShorthand["1-2-3"]}],
  Graph[{{1, 1} <-> {1, 2}, {1, 2} <-> {1, 3}, {2, 1} <-> {2, 2}, {2, 2} <-> {2, 3}}],
  SameTest -> IGSameGraphQ
]


MT[
  IGDisjointUnion[{IGShorthand["1-2-3"], IGShorthand["a-b-c"]}],
  Graph[{{1, 1} <-> {1, 2}, {1, 2} <-> {1, 3}, {2, "a"} <-> {2, "b"}, {2, "b"} <-> {2, "c"}}],
  SameTest -> IGSameGraphQ
]


MT[
  IGDisjointUnion[{IGShorthand["1-2-3"], IGShorthand["a-b"], CycleGraph[5]}],
  Graph[{{1, 1}, {1, 2}, {1, 3}, {2, "a"}, {2, "b"}, {3, 1}, {3, 2}, {3, 3}, {3, 4}, {3, 5}}, {{1, 1} <-> {1, 2}, {1, 2} <-> {1, 3}, {2, "a"} <-> {2, "b"}, {3, 1} <-> {3, 2}, {3, 1} <-> {3, 5}, {3, 2} <-> {3, 3}, {3, 3} <-> {3, 4}, {3, 4} <-> {3, 5}}],
  SameTest -> IGSameGraphQ
]


MT[
  IGDisjointUnion[Association["one" -> IGShorthand["1-2-3"], "two" -> IGShorthand["a-b"], "three" -> CycleGraph[5]]],
  Graph[{{"one", 1}, {"one", 2}, {"one", 3}, {"two", "a"}, {"two", "b"}, {"three", 1}, {"three", 2}, {"three", 3}, {"three", 4}, {"three", 5}}, {{"one", 1} <-> {"one", 2}, {"one", 2} <-> {"one", 3}, {"two", "a"} <-> {"two", "b"}, {"three", 1} <-> {"three", 2}, {"three", 1} <-> {"three", 5}, {"three", 2} <-> {"three", 3}, {"three", 3} <-> {"three", 4}, {"three", 4} <-> {"three", 5}}],
  SameTest -> IGSameGraphQ
]


MT[
  IGDisjointUnion[Association["one" -> IGShorthand["1-2-3"], "two" -> IGEmptyGraph[], "three" -> CycleGraph[5]]],
  Graph[{{"one", 1}, {"one", 2}, {"one", 3}, {"three", 1}, {"three", 2}, {"three", 3}, {"three", 4}, {"three", 5}}, {{"one", 1} <-> {"one", 2}, {"one", 2} <-> {"one", 3}, {"three", 1} <-> {"three", 2}, {"three", 1} <-> {"three", 5}, {"three", 2} <-> {"three", 3}, {"three", 3} <-> {"three", 4}, {"three", 4} <-> {"three", 5}}],
  SameTest -> IGSameGraphQ
]


MT[
  IGDisjointUnion[{IGEmptyGraph[], IGEmptyGraph[]}],
  Graph[{}],
  SameTest -> IGSameGraphQ
]


MT[
  IGDisjointUnion[{IGEmptyGraph[], IGEmptyGraph[2]}],
  Graph[{{2, 1}, {2, 2}}, {}],
  SameTest -> IGSameGraphQ
]


MT[
  IGDisjointUnion[{Graph[{1 -> 2}]}],
  Graph[{{1, 1} -> {1, 2}}],
  SameTest -> IGSameGraphQ
]


MT[
  IGDisjointUnion[{Graph[{1 -> 2}], Graph[{2 -> 1}]}],
  Graph[{{1, 1} -> {1, 2}, {2, 2} -> {2, 1}}],
  SameTest -> IGSameGraphQ
]


MT[
  IGDisjointUnion[{Graph[{1 -> 2}], Graph[{2 <-> 1}]}],
  $Failed,
  {IGDisjointUnion::mixed}
]


MT[
  IGDisjointUnion[{Graph[{1 -> 2, 2 <-> 3}]}],
  $Failed,
  {IGDisjointUnion::mixed}
]


MT[
  IGDisjointUnion[{Graph[{1 -> 2, 1 -> 2}], Graph[{"a" -> "b"}]}],
  Graph[{{1, 1} -> {1, 2}, {1, 1} -> {1, 2}, {2, "a"} -> {2, "b"}}],
  SameTest -> IGSameGraphQ
]


(* ::Section::Closed:: *)
(*Import/Export*)


MTSection["Import/Export"]


(* ::Subsection:: *)
(*Import*)


(* ::Subsubsection:: *)
(*Nauty*)


(* TODO *)


(* ::Subsection:: *)
(*Export*)


(* ::Subsubsection:: *)
(*GraphML*)


MT[
  IGExportString[Graph[{Property[1,"Name"->"Ana"],Property[2,"Name"->"Bob"]},{Property[1<->2,{"MyStringProp"->"foo",EdgeWeight->123.4,"MyIntProp"->37}]}], "GraphML"],
  "<?xml version='1.0' encoding='UTF-8'?>\n<!-- created by IGraph/M, http://szhorvat.net/mathematica/IGraphM -->\n<graphml xmlns='http://graphml.graphdrawing.org/xmlns'\n    xmlns:xsi='http://www.w3.org/2001/XMLSchema-instance'\n    xsi:schemaLocation='http://graphml.graphdrawing.org/xmlns http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd'>\n <key for='edge'\n     id='e_MyIntProp'\n     attr.name='MyIntProp'\n     attr.type='long' />\n <key for='edge'\n     id='e_MyStringProp'\n     attr.name='MyStringProp'\n     attr.type='string' />\n <key for='edge'\n     id='e_EdgeWeight'\n     attr.name='EdgeWeight'\n     attr.type='double' />\n <key for='node'\n     id='v_Name'\n     attr.name='Name'\n     attr.type='string' />\n <graph id='Graph'\n     edgedefault='undirected'>\n  <node id='1'>\n   <data key='v_Name'>Ana</data>\n  </node>\n  <node id='2'>\n   <data key='v_Name'>Bob</data>\n  </node>\n  <edge source='1'\n      target='2'>\n   <data key='e_MyIntProp'>37</data>\n   <data key='e_MyStringProp'>foo</data>\n   <data key='e_EdgeWeight'>123.4</data>\n  </edge>\n </graph>\n</graphml>"
]


MT[
  IGExportString[IGEmptyGraph[], "GraphML"],
  "<?xml version='1.0' encoding='UTF-8'?>\n<!-- created by IGraph/M, http://szhorvat.net/mathematica/IGraphM -->\n<graphml xmlns='http://graphml.graphdrawing.org/xmlns'\n    xmlns:xsi='http://www.w3.org/2001/XMLSchema-instance'\n    xsi:schemaLocation='http://graphml.graphdrawing.org/xmlns http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd'>\n <graph id='Graph'\n     edgedefault='undirected' />\n</graphml>"
]


MT[
  IGExportString[IGShorthand["a,b-c"], "GraphML"],
  "<?xml version='1.0' encoding='UTF-8'?>\n<!-- created by IGraph/M, http://szhorvat.net/mathematica/IGraphM -->\n<graphml xmlns='http://graphml.graphdrawing.org/xmlns'\n    xmlns:xsi='http://www.w3.org/2001/XMLSchema-instance'\n    xsi:schemaLocation='http://graphml.graphdrawing.org/xmlns http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd'>\n <graph id='Graph'\n     edgedefault='undirected'>\n  <node id='a' />\n  <node id='b' />\n  <node id='c' />\n  <edge source='b'\n      target='c' />\n </graph>\n</graphml>"
]


(* TODO check multigraph handling with and without properties *)


(* ::Section::Closed:: *)
(*Processes*)


(* ::Subsubsection::Closed:: *)
(*IGSIRProcess*)


(* ::Text:: *)
(*We check the validity of the returned TimeSeries or TemporalData by verifying that it responds to property requests.*)


MT[
  IGSIRProcess[IGEmptyGraph[1], {1, 1}, 5]["Caller"],
  TemporalData
]

MT[
  IGSIRProcess[IGEmptyGraph[1], {1, 1}]["Caller"],
  TimeSeries
]

MT[
  IGSIRProcess[CompleteGraph[5], {1, 1}]["FirstTime"],
  0.
]

MT[
  IGSIRProcess[IGEmptyGraph[3], {1, 1}]["ValueDimensions"],
  3
]


MT[
  IGSIRProcess[IGEmptyGraph[], {1, 1}],
  $Failed,
  {IGraphM::error}
]

MT[
  IGSIRProcess[IGEmptyGraph[], {1, 1}, 5],
  $Failed,
  {IGraphM::error}
]


(* ::Section::Closed:: *)
(*Q functions*)


MTSection["Q functions"]


(* ::Subsubsection::Closed:: *)
(*IGNullGraphQ*)


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
  IGNullGraphQ[IGEmptyGraph[1]],
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

MT[
  IGNullGraphQ[dgi],
  False
]

MT[
  IGNullGraphQ[dgs],
  False
]


(* ::Subsubsection::Closed:: *)
(*IGCactusQ*)


(* null graph is not cactus *)
MT[
  IGCactusQ[IGEmptyGraph[]],
  False
]

(* singleton graph is cactus *)
MT[
  IGCactusQ[IGEmptyGraph[1]],
  True
]

(* non-connected *)
MT[
  IGCactusQ[IGEmptyGraph[2]],
  False
]

(* P2 *)
MT[
  IGCactusQ[CompleteGraph[2]],
  True
]

(* triangle *)
MT[
  IGCactusQ[CompleteGraph[3]],
  True
]

MT[
  IGCactusQ[CompleteGraph[4]],
  False
]

MT[
  IGCactusQ[CycleGraph[4]],
  True
]

(* with bridge *)
MT[
  IGCactusQ[IGShorthand["1-2-3-1,1-4,4-5-6-4"]],
  True
]

(* without bridge *)
MT[
  IGCactusQ[IGShorthand["1-2-3-1,1-5-6-1"]],
  True
]


(* ::Text:: *)
(*Multigraphs*)


(* 2-edge *)
MT[
  IGCactusQ[IGShorthand["1-2-1", MultiEdges -> True]],
  True
]

(* 3-edge *)
MT[
  IGCactusQ[IGShorthand["1-2-1-2", MultiEdges -> True]],
  False
]

MT[
  IGCactusQ[IGShorthand["1-2-3-1,1-4-1,4-5-6-4", MultiEdges -> True]],
  True
]

MT[
  IGCactusQ[IGShorthand["1-2-3-1,1-4-1-4,4-5-6-4", MultiEdges -> True]],
  False
]

MT[
  IGCactusQ[IGShorthand["1-2-3-1,1-4-1,4-5-6-4-5", MultiEdges -> True]],
  False
]


(* ::Text:: *)
(*Self-loops must be ignored*)


MT[
  IGCactusQ[EdgeAdd[CycleGraph[5], UndirectedEdge[1, 1]]],
  True
]

MT[
  IGCactusQ[Graph[{1 <-> 1}]],
  True
]

MT[
  IGCactusQ[Graph[{1 <-> 1, 1 <-> 1}]],
  True
]


(* ::Text:: *)
(*Directed graphs*)


MT[
  IGCactusQ[Graph[{1 -> 2}]],
  $Failed,
  {IGCactusQ::dirg}
]


(* ::Subsubsection::Closed:: *)
(*IGCompleteQ*)


MT[
  IGCompleteQ[IGEmptyGraph[]],
  True
]


MT[
  IGCompleteQ[IGEmptyGraph[1]],
  True
]


MT[
  IGCompleteQ[IGEmptyGraph[2]],
  False
]


MT[
  IGCompleteQ[CompleteGraph[#]],
  True
]& /@ Range[2,5]


MT[
  IGCompleteQ[Graph[{1 <-> 1}]],
  True
]


MT[
  IGCompleteQ[Graph[{1 -> 1}]],
  True
]


MT[
  IGCompleteQ[Graph[{1 -> 2}]],
  False
]


MT[
  IGCompleteQ[Graph[{1 -> 2, 2 -> 1}]],
  True
]


MT[
  IGCompleteQ[Graph[{1 -> 2, 2 -> 1, 1 -> 1}]],
  True
]


MT[
  IGCompleteQ[Graph[{1 -> 2, 2 -> 1, 1 -> 1, 2 -> 1}]],
  True
]


MT[
  IGCompleteQ[Graph[{1 -> 2, 1 -> 2}]],
  False
]


MT[
  IGCompleteQ[CompleteGraph[4, DirectedEdges -> True]],
  True
]


MT[
  IGCompleteQ[Graph[{1 -> 2, 1 <-> 2}]],
  False,
  {IGraphM::mixed}
]


(* ::Subsubsection::Closed:: *)
(*IGRegularQ*)


MT[
  IGRegularQ[IGEmptyGraph[]],
  True
]


MT[
  IGRegularQ[IGEmptyGraph[1]],
  True
]


MT[
  IGRegularQ[Graph[{1 -> 2, 2 -> 1}]],
  True
]


MT[
  IGRegularQ[Graph[{1 -> 2, 1 -> 2, 2 -> 1}]],
  False
]


MT[
  IGRegularQ[Graph[{1 -> 2, 1 -> 2, 2 -> 1, 2 -> 1}]],
  True
]


MT[
  IGRegularQ[Graph[{1 -> 2, 1 <-> 2}]],
  False,
  {IGraphM::mixed}
]


MT[
  IGRegularQ[Graph[{1 <-> 2}]],
  True
]


MT[
  IGRegularQ[CompleteGraph[5]],
  True
]


MT[
  IGRegularQ[Graph[{1 <-> 2, 3 <-> 4}]],
  True
]


MT[
  IGRegularQ[IGEmptyGraph[], 1],
  False
]


MT[
  IGRegularQ[IGEmptyGraph[], 0],
  True
]


MT[
  IGRegularQ[IGEmptyGraph[5], 0],
  True
]


MT[
  IGRegularQ[Graph[{1 <-> 2, 1 <-> 2}], 1],
  False
]

MT[
  IGRegularQ[Graph[{1 <-> 2, 1 <-> 2}], 2],
  True
]


(* non-graph input *)

MT[
  IGRegularQ[1],
  False
]

MT[
  IGRegularQ[1, 2],
  False
]


(* ::Subsubsection::Closed:: *)
(*IGStronglyRegularQ and IGStronglyRegularParameters*)


(* TODO *)


MT[
  IGStronglyRegularParameters[IGEmptyGraph[]],
  {0, 0, 0, 0}
]

MT[
  IGStronglyRegularParameters[IGEmptyGraph[1]],
  {1, 0, 0, 0}
]

MT[
  IGStronglyRegularParameters[IGEmptyGraph[4]],
  {4, 0, 0, 0}
]

MT[
  IGStronglyRegularParameters[CompleteGraph[4]],
  {4, 3, 2, 0}
]


MT[
  IGStronglyRegularParameters[GraphData[{"Paley", 13}]],
  {13, 6, 2, 3}
]

MT[
  IGStronglyRegularParameters[PetersenGraph[]],
  {10, 3, 0, 1}
]

MT[
  IGStronglyRegularParameters[GraphData[{"CocktailParty", 7}]],
  {14, 12, 10, 12}
]

MT[
  IGStronglyRegularParameters[GraphData["M22Graph"]],
  {77, 16, 0, 4}
]

MT[
  IGStronglyRegularParameters[HypercubeGraph[2]],
  {4, 2, 0, 2}
]


MT[
  IGStronglyRegularParameters[HypercubeGraph[3]],
  {}
]


MT[
  IGStronglyRegularParameters[IGSymmetricTree[{3, 2}]],
  {}
]


(* ::Subsubsection::Closed:: *)
(*IGDistanceRegularQ and IGIntersectionArray*)


(* TODO *)


(* ::Subsubsection::Closed:: *)
(*IGVertexTransitiveQ*)


MT[
  IGVertexTransitiveQ[IGEmptyGraph[#]],
  True
]& /@ Range[0, 4]


MT[
  IGVertexTransitiveQ[IGCompleteGraph[#]],
  True
]& /@ Range[2, 5]


MT[
  IGVertexTransitiveQ[IGSquareLattice[{2, 5}, "Periodic" -> True]],
  True
]


MT[
  IGVertexTransitiveQ[IGSquareLattice[{8, 5}, "Periodic" -> True]],
  True
]


MT[
  IGVertexTransitiveQ[CycleGraph[5, DirectedEdges -> True]],
  True
]


MT[
  IGVertexTransitiveQ[IGShorthand["3->5->2->6->3,5->1->6->4->5,3->1->2->4->3"]],
  True
]


MT[
  IGVertexTransitiveQ /@ IGData[{"AllDirectedGraphs", 4}],
  {True, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, True, False, False, False, False, False, False, False, True, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, True, False, False, False, False, True, False, False, False, False, False, False, False, False, True}
]


MT[
  IGVertexTransitiveQ /@ IGData[{"AllUndirectedGraphs", 4}],
  {True, False, False, False, False, True, False, False, True, False, True}
]


MT[
  IGVertexTransitiveQ[IGShorthand["1-2-3"]],
  False
]


MT[
  IGVertexTransitiveQ[Graph[{1 <-> 2, 1 <-> 1}]],
  False
]

MT[
  IGVertexTransitiveQ[Graph[{1 <-> 2, 1 <-> 1, 2 <-> 2}]],
  True
]


MT[
  IGVertexTransitiveQ[Graph[{1 <-> 1, 1 <-> 2, 2 <-> 2}]],
  True
]


MT[
  IGVertexTransitiveQ[Graph[{1 <-> 2, 1 <-> 2}]],
  $Failed,
  {IGVertexTransitiveQ::nmg}
]


MT[
  IGVertexTransitiveQ[#],
  False
]& /@ asymmList


(* ::Subsubsection::Closed:: *)
(*IGEdgeTransitiveQ*)


MT[
  IGEdgeTransitiveQ[IGEmptyGraph[#]],
  True
]& /@ Range[0,4]

MT[
  IGEdgeTransitiveQ[IGCompleteGraph[#]],
  True
]& /@ Range[2,5]


MT[
  IGEdgeTransitiveQ[Graph[{1 -> 2}]],
  True
]


MT[
  IGEdgeTransitiveQ[Graph[{1, 2, 3}, {1 -> 2}]],
  True
]


MT[
  IGEdgeTransitiveQ[Graph[{1, 2, 3}, {1 <-> 2}]],
  True
]


MT[
  IGEdgeTransitiveQ[Graph[{1, 2}, {1 <-> 1}]],
  True
]


MT[
  IGEdgeTransitiveQ[Graph[{1, 2, 3}, {1 <-> 1, 3 <-> 3}]],
  True
]


MT[
  IGEdgeTransitiveQ[Graph[{1 -> 2, 2 -> 3}]],
  False
]


MT[
  IGEdgeTransitiveQ[Graph[{1 -> 2, 3 -> 2}]],
  True
]


MT[
  IGEdgeTransitiveQ[StarGraph[4, DirectedEdges -> True]],
  True
]


MT[
  IGEdgeTransitiveQ[ReverseGraph[StarGraph[4, DirectedEdges -> True]]],
  True
]


MT[
  IGEdgeTransitiveQ[IGShorthand["1->2->3,2->4"]],
  False
]


MT[
  IGEdgeTransitiveQ[Graph[{1 -> 2, 1 -> 2}]],
  $Failed,
  {IGEdgeTransitiveQ::nmg}
]


MT[
  IGEdgeTransitiveQ[Graph[{1 <-> 2, 1 <-> 2}]],
  $Failed,
  {IGEdgeTransitiveQ::nmg}
]


MT[
  IGEdgeTransitiveQ[#],
  False
]& /@ asymmList


(* This graph is not edge transitive, but its line graph is vertex transitive. *)
MT[
  IGEdgeTransitiveQ@GraphDisjointUnion[CycleGraph[3], StarGraph[4]],
  False
]


(* ::Subsubsection::Closed:: *)
(*IGSymmetricQ*)


(* TODO *)


MT[
  IGSymmetricQ[IGEmptyGraph[#]],
  True
]& /@ Range[0,3]

MT[
  IGSymmetricQ[IGCompleteGraph[#]],
  True
]& /@ Range[2,5]


MT[
  IGSymmetricQ[GraphData["DoyleGraph"]],
  True
]


MT[
  IGSymmetricQ[GraphData["ShrikhandeGraph"]],
  True
]


MT[
  IGSymmetricQ[#],
  False
]& /@ asymmList


(* ::Subsubsection::Closed:: *)
(*IGDistanceTransitiveQ*)


MT[
  IGDistanceTransitiveQ[IGEmptyGraph[#]],
  True
]& /@ Range[0, 4]

MT[
  IGDistanceTransitiveQ[IGCompleteGraph[#]],
  True
]& /@ Range[2, 4]


MT[
  IGDistanceTransitiveQ[IGShorthand["1-2,3"]],
  False
]


MT[
  IGDistanceTransitiveQ[IGShorthand["1-2,3-4"]],
  True
]


MT[
  IGDistanceTransitiveQ[IGShorthand["1->2"]],
  False
]


MT[
  IGDistanceTransitiveQ[IGShorthand["1->2->3->1"]],
  True
]


MT[
  IGDistanceTransitiveQ[IGShorthand["1<->2"]],
  True
]


MT[
  IGDistanceTransitiveQ[GraphData["ShrikhandeGraph"]],
  False
]


MT[
  IGDistanceTransitiveQ[GraphData[{"Rook", {4, 4}}]],
  True
]


MT[
  IGDistanceTransitiveQ[GraphData["IcosahedralGraph"]],
  True
]


MT[
  IGDistanceTransitiveQ[GraphData["DodecahedralGraph"]],
  True
]


MT[
  IGDistanceTransitiveQ /@ IGData[{"AllUndirectedGraphs", 4}],
  {True, False, False, False, False, True, False, False, True, False, True}
]


MT[
  IGDistanceTransitiveQ[#],
  False
]& /@ asymmList


(* oriented octahedral skeleton *)
MT[
  IGDistanceTransitiveQ[IGShorthand["3->5->2->6->3,5->1->6->4->5,3->1->2->4->3"]],
  False
]


MT[
  IGDistanceTransitiveQ /@ IGData[{"AllDirectedGraphs", 4}],
  {True, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, True, False, False, False, False, False, False, False, True, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, True, False, False, False, False, False, False, False, False, False, False, False, False, False, True}
]


(* https://doi.org/10.1016/0012-365X(80)90155-7 *)
MT[
  Module[{n=#,mods},
    mods=Rest@Union@Mod[Range[n]^2, n];
    IGDistanceTransitiveQ@relationGraph[MemberQ[mods, Mod[#1-#2, n]]&, Range[0, n-1]]
  ],
  True
]& /@ {7, 11, 19, 23}


(* with self loops *)
MT[
  IGDistanceTransitiveQ[Graph[{1 <-> 2, 1 <-> 1}]],
  False
]

MT[
  IGDistanceTransitiveQ[Graph[{1 <-> 2, 1 <-> 1, 2 <-> 2}]],
  True
]


(* multigraph *)
MT[
  IGDistanceTransitiveQ[Graph[{1 <-> 2, 1 <-> 2}]],
  $Failed,
  {IGDistanceTransitiveQ::nmg}
]


(* ::Subsubsection::Closed:: *)
(*IGSameGraphQ*)


MT[
  IGSameGraphQ[IGEmptyGraph[], IGEmptyGraph[]],
  True
]

MT[
  IGSameGraphQ[IGEmptyGraph[], IGEmptyGraph[1]],
  False
]


(* ::Subsubsection:: *)
(*IGAdjacentVerticesQ*)


(* TODO *)


MT[
  IGAdjacentVerticesQ[Graph[{1 <-> 2}], {1, 2}],
  True
]

MT[
  IGAdjacentVerticesQ[Graph[{1 <-> 2}], {2, 1}],
  True
]


(* ::Subsection::Closed:: *)
(*General tests*)


(* ::Text:: *)
(*Check that Q functions that take a certain expression type (such as a graph) will return False when given a different expression type.*)


(MT[#[1234], False]; MT[#[DUMMY], False])& /@ {
	IGBiconnectedQ, IGConnectedQ, IGWeaklyConnectedQ,
	IGBipartiteQ, 
	IGChordalQ,
	IGDirectedAcyclicGraphQ, 
	IGEdgeTransitiveQ, IGVertexTransitiveQ, IGSymmetricQ, IGDistanceTransitiveQ,
	IGRegularQ, IGStronglyRegularQ, IGDistanceRegularQ,
	IGSelfComplementaryQ,
	IGEdgeWeightedQ, IGVertexWeightedQ,
	IGForestQ, IGTreeQ,
	IGPlanarQ, IGMaximalPlanarQ, IGOuterplanarQ,
	IGNullGraphQ,
	IGTriangleFreeQ,
	IGPerfectQ,
	IGCactusQ,
	IGCompleteQ
}	


MT[IGGraphicalQ[123], False]
MT[IGGraphicalQ[123, 456], False]

MT[IGGraphicalQ[DUMMY], False]
MT[IGGraphicalQ[DUMMY1, DUMMY2], False]


MT[IGEmbeddingQ[123], False]
MT[IGEmbeddingQ[{1,2,3}], False]


MT[
  IGSameGraphQ[1, 2],
  False
]

MT[
  IGSameGraphQ[IGEmptyGraph[], 2],
  False
]


(* ::Section::Closed:: *)
(*Miscellaneous functions*)


MTSection["Miscellaneous functions"]


(* ::Subsubsection::Closed:: *)
(*IGData*)


MT[
  IGData /@ IGData[];,
  Null
]

MT[
  And @@ MapThread[IGIsomorphicQ, {IGData[{"AllDirectedGraphs", 3}], Values@IGData["MANTriadLabels"]}],
  True
]


MT[
  IGIsoclass/@IGData[{"AllDirectedGraphs", 3}],
  Range[0, 15]
]


MT[
  IGIsoclass/@IGData[{"AllDirectedGraphs", 4}],
  Range[0, 217]
]


MT[
  IGIsoclass/@IGData[{"AllUndirectedGraphs", 3}],
  Range[0, 3]
]


MT[
  IGIsoclass/@IGData[{"AllUndirectedGraphs", 4}],
  Range[0, 10]
]


(* ::Subsubsection::Closed:: *)
(*IGIndexEdgeList*)


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


(* ::Subsubsection::Closed:: *)
(*IGraphM`Developer`GetInfo*)


(* exercise GetInfo *)
MT[
  StringQ@IGraphM`Developer`GetInfo[],
  True
]


(* ::Subsubsection::Closed:: *)
(*IGraphM`Information`$Version*)


MT[
  StringQ[IGraphM`Information`$Version],
  True
]


(* ::Section::Closed:: *)
(*Finalize*)


MTSection["Finalize"]


(* No IG managed library expressions should be alive *)
MT[
  IGraphM`LTemplate`LExpressionList["IG"],
  {}
]


MT[
  IGraphM`LTemplate`LExpressionList["IGEmbedding"],
  {}
]


MT[
  IGraphM`LTemplate`LExpressionList["IGLemonGraph"],
  {}
]


MT[
  IGraphM`LTemplate`LExpressionList["IGFlann2D"],
  {}
]


(* There should be precisely one IGlobal object at any time. *)
MT[
  IGraphM`LTemplate`LExpressionList["IGlobal"],
  {IGraphM`LTemplate`Classes`IGlobal[1]}
]
