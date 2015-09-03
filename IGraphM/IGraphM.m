(* Mathematica Package  *)
(* Created by IntelliJ IDEA and wlplugin.halirutan.de *)

(* :Title: IGraph/M   *)
(* :Context: IGraphM` *)
(* :Author: szhorvat  *)
(* :Date: 2015-08-28  *)

(* :Package Version: 0.1pre *)
(* :Mathematica Version: 10.0 *)
(* :Copyright: (c) 2015 Szabolcs Horvát *)
(* :Keywords: igraph, graphs, networks, LibraryLink *)
(* :Discussion: igraph interface for Mathematica, see http://igraph.org *)

BeginPackage["IGraphM`"]

(* Privately load and configure LTemplate *)
Get["LTemplate`LTemplatePrivate`"]
ConfigureLTemplate["MessageSymbol" -> IGraphM]

(***** Usage mesages *****)

IGraphM::usage = "IGraphM is a symbol to which igraph related messages are associated.";

RecompileIGraphM::usage = "RecompileIGraphM[]";

IGData::usage =
    "IGData[] returns a list of available items.\n" <>
    "IGData[item] returns the requested item.";

IGVersion::usage = "IGVersion[] returns the version of the igraph library in use.";
IGSeedRandom::usage = "IGSeedRandom[seed] seeds the random number generator used by igraph.";

IGBetweenness::usage = "IGBetweenness[graph]";
IGEdgeBetweenness::usage = "IGEdgeBetweenness[graph]";
IGCloseness::usage = "IGCloseness[graph, options]";

IGBetweennessEstimate::usage = "IGBetweennessEstimate[graph, cutoff]";
IGEdgeBetweennessEstimate::usage = "IGEdgeBetweennessEstimate[graph, cutoff]";
IGClosenessEstimate::usage = "IGClosenessEstimate[graph, cutoff, options]";

IGRewire::usage = "IGRewire[graph, n, options]";
IGRewireEdges::usage = "IGRewireEdges[graph, p, options]";

IGDirectedAcyclicGraphQ::usage = "IGDirectedAcyclicGraphQ[graph]";
IGConnectedQ::usage = "IGConnectedQ[graph]";
IGGraphicalQ::usage =
    "IGGraphicalQ[degrees] tests if a degree sequence for an undirected simple graph is graphical.\n" <>
    "IGGraphicalQ[indegrees, outdegrees] tests if a degree sequence for a directed simple graph is graphical.";

IGIsomorphicQ::usage = "IGIsomorphicQ[graph1, graph2]";
IGSubisomorphicQ::usage = "IGSubisomorphicQ[graph, subgraph]";
IGIsoclass::usage = "IGIsoclass[graph] returns the isomorphism class of the graph. Used as the index into the vector returned by motif finding functions.";

IGBlissCanonicalPermutation::usage =
    "IGBlissCanonicalPermutation[graph, options] computes a canonical permutation of the graph vertices. " <>
    "Two graphs are isomorphic iff they have the same canonical permutation.";
IGBlissIsomorphicQ::usage = "IGBlissIsomorphicQ[graph1, graph2, options]";
IGBlissFindIsomorphism::usage = "IGBlissFindIsomorphism[graph1, graph2, options]";
IGBlissAutomorphismsCount::usage = "IGBlissAutomorphismsCount[graph]";

IGTopologicalOrdering::usage = "IGTopologicalOrdering[graph] returns a permutation that sorts the vertices in topological order. Note that the values returned are vertex indices, not vertex names.";
IGFeedbackArcSet::usage = "IGFeedbackArcSet[graph]";

IGDyadCensus::usage = "IGDyadCensus[graph]";
IGTriadCensus::usage = "IGTriadCensus[graph]";
IGMotifs::usage = "IGMotifs[graph, motifSize] returns the motif distribution of graph. See IGIsoclass for motif ordering.";
IGMotifsTotalCount::usage = "IGMotifsTotalCount[graph, motifSize]";
IGMotifsEstimateTotalCount::usage = "IGMotifsEstimate[graph, motifSize, sampleSize]";

IGDegreeSequenceGame::usage =
    "IGDegreeSequenceGame[degrees, options] generates an undirected random graph with the given degree sequence.\n" <>
    "IGDegreeSequenceGame[indegrees, outdegrees, options] generates a directed random graph with the given in- and out-degree sequences.";

IGDistanceMatrix::usage = "IGDistanceMatrix[graph]";

IGCliques::usage =
    "IGCliques[graph] returns all complete subgraphs (cliques) in graph. Note that this is different from the builtin FindCliques[], which finds maximal cliques.\n" <>
    "IGCliques[graph, {min, max}] returns all complete subgraphs between sizes min and max.\n" <>
    "IGCliques[graph, max] returns all complete subgraphs of size at most max.\n" <>
    "IGCliques[graph, {n}] returns all complete subgraphs of size n.\n";
IGMaximalCliques::usage =
    "IGMaximalCliques[graph] returns all maximal cliques in graph.\n" <>
    "IGMaximalCliques[graph, {min, max}] returns all maximal cliques between sizes min and max.\n" <>
    "IGMaximalCliques[graph, max] returns all maximal cliques of size at most max.\n" <>
    "IGMaximalCliques[graph, {n}] returns all maximal cliques of size n.\n";
IGMaximalCliquesCount::usage =
    "IGMaximalCliquesCount[graph] counts all maximal cliques in graph.\n" <>
    "IGMaximalCliquesCount[graph, {min, max}] counts all maximal cliques between sizes min and max.\n" <>
    "IGMaximalCliquesCount[graph, max] counts all maximal cliques of size at most max.\n" <>
    "IGMaximalCliquesCount[graph, {n}] counts all maximal cliques of size n.\n";
IGLargestCliques::usage = "IGLargestCliques[graph] returns the largest cliques in graph.";
IGCliqueNumber::usage = "IGCliqueNumber[graph] returns the clique number of graph. The clique number is the size of the largest clique.";

IGIndependentVertexSets::usage =
    "IGIndependentVertexSets[graphs] finds all independent vertex sets of graph.\n" <>
    "IGIndependentVertexSets[graphs, {min, max}]\n" <>
    "IGIndependentVertexSets[graphs, max]\n" <>
    "IGIndependentVertexSets[graphs, {n}]";
IGLargestIndependentVertexSets::usage = "IGLargestIndependentVertexSets[graph] finds the largest independent vertex sets of graph.";
IGMaximalIndependentVertexSets::usage = "IGMaximalIndependentVertexSets[graph] finds the maximal independent vertex sets of graph.";
IGIndependenceNumber::usage = "IGIndependenceNumber[graph] returns the independence number of graph. The independence number is the size of the largest independent vertex set.";

Begin["`Private`"]

(***** Mathematica version check *****)

(* Abort loading and leave a clean $ContextPath behind *)
packageAbort[] := (End[]; EndPackage[]; Abort[])

If[Not@OrderedQ[{10.0, 2}, {$VersionNumber, $ReleaseNumber}],
  Print["IGraphM requires Mathematica 10.0.2 or later.  Aborting."];
  packageAbort[]
]


(***** Package variables *****)

$packageVersion    = "0.1pre";
$packageDirectory  = DirectoryName[$InputFileName];
$libraryDirectory  = FileNameJoin[{$packageDirectory, "LibraryResources", $SystemID}];
$sourceDirectory   = FileNameJoin[{$packageDirectory, "LibraryResources", "Source"}];
$buildSettingsFile = FileNameJoin[{$packageDirectory, "BuildSettings.m"}];


template = LTemplate["IGraphM",
  {
    LClass["IGlobal",
      {
        LFun["init", {}, "Void"],
        LFun["seedRandom", {Integer}, "Void"],
        LFun["version", {}, "UTF8String"],

        (* Graph related functions that do not use the graph data structure *)

        LFun["graphicalQ", {{Real, 1, "Constant"} (* outdeg *), {Real, 1, "Constant"} (* indeg *)}, True|False]
      }
    ],
    LClass["IG",
      {
        (* Create *)

        LFun["fromEdgeList", {{Real, 2, "Constant"} (* edges *), Integer (* vertex count *), True|False (* directed *)}, "Void"],

        (* Weights *)

        LFun["setWeights", {{Real, 1, "Constant"}}, "Void"],
        LFun["getWeights", {}, {Real, 1}],
        LFun["clearWeights", {}, "Void"],
        LFun["weightedQ", {}, True|False],

        (* Games *)

        LFun["degreeSequenceGame", {{Real, 1, "Constant"} (* outdeg *), {Real, 1, "Constant"} (* indeg *), Integer (* method *)}, "Void"],

        (* Structure *)

        LFun["edgeCount", {}, Integer],
        LFun["vertexCount", {}, Integer],

        LFun["edgeList", {}, {Real, 2}],

        (* Testing *)

        LFun["directedQ", {}, True|False],
        LFun["dagQ", {}, True|False],
        LFun["simpleQ", {}, True|False],
        LFun["connectedQ", {}, True|False],

        (* Centrality *)

        LFun["betweenness", {}, {Real, 1}],
        LFun["edgeBetweenness", {}, {Real, 1}],
        LFun["closeness", {True|False (* normalized *)}, {Real, 1}],

        LFun["betweennessEstimate", {Real (* cutoff *)}, {Real, 1}],
        LFun["edgeBetweennessEstimate", {Real (* cutoff *)}, {Real, 1}],
        LFun["closenessEstimate", {Real (* cutoff *), True|False (* normalized *)}, {Real, 1}],

        (* Randomize *)

        LFun["rewire", {Integer (* n_trials *), True|False (* loops *)}, "Void"],
        LFun["rewireEdges", {Real (* probability *), True|False (* loops *), True|False (* multiple *)}, "Void"],

        (* Isomorphism *)

        LFun["isomorphic", {LExpressionID["IG"]}, True|False],
        LFun["subisomorphic", {LExpressionID["IG"]}, True|False],
        LFun["isoclass", {}, Integer],
        LFun["blissCanonicalPermutation", {Integer (* splitting heuristics *)}, {Real, 1}],
        LFun["blissIsomorphic", {LExpressionID["IG"], Integer (* splitting heuristics *)}, True|False],
        LFun["blissFindIsomorphism", {LExpressionID["IG"], Integer (* splitting heuristics *)}, {Real, 1}],
        LFun["blissAutomorphismsCount", LinkObject],

        (* Topological sorting and directed acylic graphs *)

        LFun["topologicalSorting", {}, {Real, 1}],
        LFun["feedbackArcSet", {True|False}, {Real, 1}],

        (* Motifs and subgraph counts *)

        LFun["dyadCensus", {}, {Integer, 1}],
        LFun["triadCensus", {}, {Real, 1}],
        LFun["motifs", {Integer (* size *), {Real, 1, "Constant"} (* cut_prob *)}, {Real, 1}],
        LFun["motifsNo", {Integer (* size *), {Real, 1, "Constant"} (* cut_prob *)}, Integer],
        LFun["motifsEstimate", {Integer (* size *), {Real, 1, "Constant"} (* cut_prob *), Integer (* sample_size *)}, Integer],

        (* Shortest paths *)

        LFun["shortestPaths", {}, {Real, 2}],

        (* Cliques *)

        LFun["cliques", LinkObject],
        LFun["maximalCliques", LinkObject],
        LFun["largestCliques", LinkObject],
        LFun["maximalCliquesCount", {Integer, Integer}, Integer],
        LFun["cliqueNumber", {}, Integer],

        (* Independent vertex sets *)

        LFun["independentVertexSets", LinkObject],
        LFun["largestIndependentVertexSets", LinkObject],
        LFun["maximalIndependentVertexSets", LinkObject],
        LFun["independenceNumber", {}, Integer]
      }
    ]
  }
];


(***** Compilation, loading and initialization *****)

$buildSettings = None;
If[FileExistsQ[$buildSettingsFile], Get[$buildSettingsFile] ]


(* Add $libraryDirectory to $LibraryPath in case package is not installed in Applications. *)
If[Not@MemberQ[$LibraryPath, $libraryDirectory],
  AppendTo[$LibraryPath, $libraryDirectory]
]


RecompileIGraphM::build = "No build settings found. Please check BuildSettings.m."

RecompileIGraphM[] :=
    Module[{},
      If[$buildSettings === None,
        Message[RecompileIGraphM::build];
        Return[$Failed]
      ];
      If[Not@DirectoryQ[$libraryDirectory],
        CreateDirectory[$libraryDirectory]
      ];
      SetDirectory[$sourceDirectory];
      Quiet@UnloadTemplate[template];
      CompileTemplate[template, {"IGlobal.cpp"},
        "ShellCommandFunction" -> Print, "ShellOutputFunction" -> Print,
        "TargetDirectory" -> $libraryDirectory,
        Sequence @@ $buildSettings
      ];
      ResetDirectory[];
      LoadIGraphM[]
    ]


igraphGlobal (* there should only be a single object of this type; it's set in LoadIGraphM[] below *)

LoadIGraphM[] :=
    If[LoadTemplate[template] === $Failed,
      $Failed
      ,
      igraphGlobal = Make["IGlobal"];
      igraphGlobal@"init"[];
    ]


(* Load library, compile if necessary. *)
If[LoadIGraphM[] === $Failed,
  Print[Style["Loading failed, trying to recompile ...", Red]];
  If[RecompileIGraphM[] === $Failed
    ,
    Print[Style["Cannot load or compile library. \[FreakedSmiley] Aborting.", Red]];
    packageAbort[]
    ,
    Print[Style["Successfully compiled and loaded the library. \[HappySmiley]", Red]];
  ]
]


(***** Helper functions *****)

nonNegIntVecQ = VectorQ[#, Internal`NonNegativeMachineIntegerQ]&

(* Zero out the diagonal of a square matrix. *)
zeroDiagonal[arg_] := UpperTriangularize[arg, 1] + LowerTriangularize[arg, -1]

(* Import compressed expressions. *)
zimport[filename_] := Uncompress@Import[filename, "String"]

(* Get an IG compatible edge list. *)
igEdgeList[g_?GraphQ] :=
    Developer`ToPackedArray@N[List @@@ EdgeList[g] /.
            Dispatch@Thread[VertexList[g] -> Range@VertexCount[g] - 1]]

(* Convert IG format vertex or edge index vector to Mathematica format. *)
igIndexVec[expr_LibraryFunctionError] := expr (* hack: allows LibraryFunctionError to fall through *)
igIndexVec[arr_] := 1 + Round[arr]

igDirectedQ[g_?GraphQ] := DirectedGraphQ[g] && Not@EmptyGraphQ[g]

(* TODO: Find out how to implement this in a more robust way. *)
igWeightedGraphQ = WeightedGraphQ[#] && PropertyValue[#, EdgeList] =!= Automatic &;

(* Create IG object from Mathematica Graph. *)
igMake[g_?GraphQ] :=
    With[{ig = Make["IG"]},
      ig@"fromEdgeList"[igEdgeList[g], VertexCount[g], igDirectedQ[g]];
      If[igWeightedGraphQ[g], ig@"setWeights"[PropertyValue[g, EdgeWeight]]];
      ig
    ]

(* Create Mathematica Graph from IG object. *)
igToGraph[ig_] :=
    Graph[
      Range[ig@"vertexCount"[]],
      igIndexVec[ig@"edgeList"[]],
      DirectedEdges -> ig@"directedQ"[]
    ]

(* Convert vertex indices to vertex names. *)
igVertexNames[graph_][indices_] := Part[VertexList[graph], indices]

igLabelValues[_, err_LibraryFunctionError] := err
igLabelValues[labels_, values_] := AssociationThread[labels, values]


(***** Public functions *****)

(*  IGData  *)

$igData = zimport@FileNameJoin[{$packageDirectory, "IGData.mz"}];
IGData[] := Keys[$igData]
IGData[item_] := Lookup[$igData, Key[item], Missing["NotAvailable"]]

(* General (global) *)

IGVersion[] := "IGraph/M " <> $packageVersion <> ", based on igraph " <> igraphGlobal@"version"[] <> ".";

IGSeedRandom[seed_?Internal`NonNegativeMachineIntegerQ] := igraphGlobal@"seedRandom"[seed]

(* Create (games) *)

Options[IGDegreeSequenceGame] = { Method -> "SimpleNoMultiple" };

igDegreeSequenceGameMethods = <| "VigerLatapy" -> 2, "SimpleNoMultiple" -> 1, "Simple" -> 0 |>

IGDegreeSequenceGame::usage = IGDegreeSequenceGame::usage <>
    StringTemplate[" Available methods: ``"][ToString@InputForm@Keys[igDegreeSequenceGameMethods]];

IGDegreeSequenceGame[degrees_?nonNegIntVecQ, opt : OptionsPattern[]] := IGDegreeSequenceGame[{}, degrees, opt]

IGDegreeSequenceGame[indegrees_?nonNegIntVecQ, outdegrees_?nonNegIntVecQ, opt : OptionsPattern[]] :=
    Block[{ig = Make["IG"]},
      Check[
        ig@"degreeSequenceGame"[outdegrees, indegrees, Lookup[igDegreeSequenceGameMethods, OptionValue[Method], -1]];
        igToGraph[ig]
        ,
        $Failed
        ,
        {LTemplate::error}
      ]
    ]

(* Testing *)

IGDirectedAcyclicGraphQ[g_?GraphQ] := Block[{ig = igMake[g]}, ig@"dagQ"[]]

IGConnectedQ[g_?GraphQ] := Block[{ig = igMake[g]}, ig@"connectedQ"[]]

IGGraphicalQ[degrees_?nonNegIntVecQ] := IGGraphicalQ[{}, degrees]
IGGraphicalQ[indeg_?nonNegIntVecQ, outdeg_?nonNegIntVecQ] := igraphGlobal@"graphicalQ"[outdeg, indeg]

(* Centrality *)

IGBetweenness[g_?GraphQ] := Block[{ig = igMake[g]}, ig@"betweenness"[]]

IGEdgeBetweenness[g_?GraphQ] := Block[{ig = igMake[g]}, ig@"edgeBetweenness"[]]

Options[IGCloseness] = { "Normalized" -> False };
IGCloseness[g_?GraphQ, opt : OptionsPattern[]] := Block[{ig = igMake[g]}, ig@"closeness"[OptionValue["Normalized"]]]

(* Centrality estimates *)

IGBetweennessEstimate[g_?GraphQ, cutoff_] := Block[{ig = igMake[g]}, ig@"betweennessEstimate"[cutoff]]

IGEdgeBetweennessEstimate[g_?GraphQ, cutoff_] := Block[{ig = igMake[g]}, ig@"edgeBetweennessEstimate"[cutoff]]

Options[IGClosenessEstimate] = { "Normalized" -> False };
IGClosenessEstimate[g_?GraphQ, cutoff_, opt : OptionsPattern[]] := Block[{ig = igMake[g]}, ig@"closenessEstimate"[cutoff, OptionValue["Normalized"]]]

(* Randomization *)

(* TODO: functions in this section should warn that edge weights will be lost *)

Options[IGRewire] = { "AllowLoops" -> False };
IGRewire[g_?GraphQ, n_Integer, opt : OptionsPattern[]] :=
    Block[{ig = igMake[g]},
      ig@"rewire"[n, OptionValue["AllowLoops"]];
      igToGraph[ig]
    ]

Options[IGRewireEdges] = { "AllowLoops" -> False, "AllowMultipleEdges" -> False };
IGRewireEdges[g_?GraphQ, p_?Internal`RealValuedNumericQ, opt : OptionsPattern[]] :=
    Block[{ig = igMake[g]},
      ig@"rewireEdges"[p, OptionValue["AllowLoops"], OptionValue["AllowMultipleEdges"]];
      igToGraph[ig]
    ]

(* Isomorphism *)

IGIsomorphicQ[g1_?GraphQ, g2_?GraphQ] := Block[{ig1 = igMake[g1], ig2 = igMake[g2]}, ig1@"isomorphic"[ManagedLibraryExpressionID@ig2]]

IGSubisomorphicQ[graph_?GraphQ, subgraph_?GraphQ] := Block[{ig1 = igMake[graph], ig2 = igMake[subgraph]}, ig1@"subisomorphic"[ManagedLibraryExpressionID@ig2]]

IGIsoclass[graph_?GraphQ] := Block[{ig = igMake[graph]}, ig@"isoclass"[]]


igBlissSplittingHeuristicsNames = {
  "First", "FirstSmallest", "FirstLargest",
  "FirstMaximallyConnected", "FirstSmallestMaximallyConnected", "FirstLargestMaximallyConnected"
};

igBlissSplittingHeuristics = AssociationThread[igBlissSplittingHeuristicsNames, Range@Length[igBlissSplittingHeuristicsNames] - 1];

IGBlissCanonicalPermutation::usage = IGBlissCanonicalPermutation::usage <>
    StringTemplate[" Available values for the \"SplittingHeuristics\" option: ``." <>
        "The permutation depends on the splitting heuristics used."][ToString@InputForm@igBlissSplittingHeuristicsNames];

Options[IGBlissCanonicalPermutation] = { "SplittingHeuristics" -> "First" };
IGBlissCanonicalPermutation[graph_?GraphQ, opt : OptionsPattern[]] :=
    Block[{ig = igMake[graph]},
      igVertexNames[graph]@igIndexVec@ig@"blissCanonicalPermutation"[Lookup[igBlissSplittingHeuristics, OptionValue["SplittingHeuristics"], -1]]
    ]

Options[IGBlissIsomorphicQ] = { "SplittingHeuristics" -> "First" };
IGBlissIsomorphicQ[graph1_?GraphQ, graph2_?GraphQ, opt : OptionsPattern[]] :=
    Block[{ig1 = igMake[graph1], ig2 = igMake[graph2]},
      ig1@"blissIsomorphic"[ManagedLibraryExpressionID[ig2], Lookup[igBlissSplittingHeuristics, OptionValue["SplittingHeuristics"], -1]]
    ]

Options[IGBlissFindIsomorphism] = { "SplittingHeuristics" -> "First" };
IGBlissFindIsomorphism[graph1_?GraphQ, graph2_?GraphQ, opt : OptionsPattern[]] :=
    Block[{ig1 = igMake[graph1], ig2 = igMake[graph2], result},
      result = igIndexVec@ig1@"blissFindIsomorphism"[ManagedLibraryExpressionID[ig2], Lookup[igBlissSplittingHeuristics, OptionValue["SplittingHeuristics"], -1]];
      If[result =!= {}, {result}, {}]
    ]

Options[IGBlissAutomorphismsCount] = { "SplittingHeuristics" -> "First" };
IGBlissAutomorphismsCount[graph_?GraphQ, opt : OptionsPattern[]] :=
    Block[{ig = igMake[graph]},
      ToExpression@ig@"blissAutomorphismsCount"[Lookup[igBlissSplittingHeuristics, OptionValue["SplittingHeuristics"], -1]]
    ]

(* Directed acylic graphs and topological ordering *)

IGTopologicalOrdering[graph_?GraphQ] := Block[{ig = igMake[graph]}, igIndexVec@ig@"topologicalSorting"[]]

Options[IGFeedbackArcSet] = { "Exact" -> True };
IGFeedbackArcSet[graph_?GraphQ, opt : OptionsPattern[]] :=
    Block[{ig = igMake[graph]},
      Part[EdgeList[graph], igIndexVec@ig@"feedbackArcSet"[OptionValue["Exact"]]]
    ]

(* Motifs and subgraph counts *)

IGDyadCensus[graph_?GraphQ] := Block[{ig = igMake[graph]}, AssociationThread[{"Mutual", "Asymmetric", "Missing"}, Round@ig@"dyadCensus"[]]]

IGTriadCensus[graph_?GraphQ] :=
    Block[{ig = igMake[graph]},
      igLabelValues[
        {"003", "012", "102", "021D", "021U", "021C", "111D", "111U", "030T", "030C", "201", "120D", "120U", "120C", "210", "300"},
        Round@ig@"triadCensus"[]
      ]
    ]

IGMotifs[graph_?GraphQ, size_?Internal`PositiveIntegerQ] :=
    Block[{ig = igMake[graph]},
      Round@Developer`FromPackedArray@ig@"motifs"[size, ConstantArray[0, size]]
    ]

IGMotifsTotalCount[graph_?GraphQ, size_?Internal`PositiveIntegerQ] :=
    Block[{ig = igMake[graph]}, ig@"motifsNo"[size, ConstantArray[0, size]] ]

IGMotifsEstimateTotalCount[graph_?GraphQ, size_?Internal`PositiveIntegerQ, sampleSize_?Internal`PositiveIntegerQ] :=
    Block[{ig = igMake[graph]}, ig@"motifsEstimate"[size, ConstantArray[0, size], sampleSize] ]

(* Shortest paths *)

IGDistanceMatrix[graph_?GraphQ] :=
    Block[{ig = igMake[graph]},
      zeroDiagonal[Transpose@Round[ig@"shortestPaths"[]] /. 0 -> Infinity] (* TODO: avoid unpacking when no infinities present *)
    ]

(* Cliques *)

IGCliques[graph_] := IGCliques[graph, Infinity]
IGCliques[graph_, max : (_Integer | Infinity)] := IGCliques[graph, {1, max}]
IGCliques[graph_, {size_}] := IGCliques[graph, {size, size}]
IGCliques[graph_?GraphQ, {min_?Internal`PositiveMachineIntegerQ, max : (_?Internal`PositiveMachineIntegerQ | Infinity)}] /; max >= min :=
    iIGCliques[graph, {min, max /. Infinity -> 0}]
iIGCliques[graph_, {min_, max_}] :=
    Block[{ig = igMake[graph]},
      igVertexNames[graph] /@ igIndexVec@ig@"cliques"[min, max]
    ]

IGMaximalCliques[graph_] := IGMaximalCliques[graph, Infinity]
IGMaximalCliques[graph_, max : (_Integer | Infinity)] := IGMaximalCliques[graph, {1, max}]
IGMaximalCliques[graph_, {size_}] := IGMaximalCliques[graph, {size, size}]
IGMaximalCliques[graph_?GraphQ, {min_?Internal`PositiveMachineIntegerQ, max : (_?Internal`PositiveMachineIntegerQ | Infinity)}] /; max >= min :=
    iIGMaximalCliques[graph, {min, max /. Infinity -> 0}]
iIGMaximalCliques[graph_, {min_, max_}] :=
    Block[{ig = igMake[graph]},
      igVertexNames[graph] /@ igIndexVec@ig@"maximalCliques"[min, max]
    ]

IGMaximalCliquesCount[graph_] := IGMaximalCliquesCount[graph, Infinity]
IGMaximalCliquesCount[graph_, max : (_Integer | Infinity)] := IGMaximalCliquesCount[graph, {1, max}]
IGMaximalCliquesCount[graph_, {size_}] := IGMaximalCliquesCount[graph, {size, size}]
IGMaximalCliquesCount[graph_?GraphQ, {min_?Internal`PositiveMachineIntegerQ, max : (_?Internal`PositiveMachineIntegerQ | Infinity)}] /; max >= min :=
    iIGMaximalCliquesCount[graph, {min, max /. Infinity -> 0}]
iIGMaximalCliquesCount[graph_, {min_, max_}] :=
    Block[{ig = igMake[graph]},
      Round@ig@"maximalCliquesCount"[min, max]
    ]

IGLargestCliques[graph_?GraphQ] :=
    Block[{ig = igMake[graph]},
      igVertexNames[graph] /@ igIndexVec@ig@"largestCliques"[]
    ]

IGCliqueNumber[graph_?GraphQ] := Block[{ig = igMake[graph]}, ig@"cliqueNumber"[]]

(* Independent vertex sets *)

IGIndependentVertexSets[graph_] := IGIndependentVertexSets[graph, Infinity]
IGIndependentVertexSets[graph_, max : (_Integer | Infinity)] := IGIndependentVertexSets[graph, {1, max}]
IGIndependentVertexSets[graph_, {size_}] := IGIndependentVertexSets[graph, {size, size}]
IGIndependentVertexSets[graph_?GraphQ, {min_?Internal`PositiveMachineIntegerQ, max : (_?Internal`PositiveMachineIntegerQ | Infinity)}] /; max >= min :=
    iIGIndependentVertexSets[graph, {min, max /. Infinity -> 0}]
iIGIndependentVertexSets[graph_, {min_, max_}] :=
    Block[{ig = igMake[graph]},
      igVertexNames[graph] /@ igIndexVec@ig@"independentVertexSets"[min, max]
    ]

IGLargestIndependentVertexSets[graph_?GraphQ] :=
    Block[{ig = igMake[graph]},
      igVertexNames[graph] /@ igIndexVec@ig@"largestIndependentVertexSets"[]
    ]

IGMaximalIndependentVertexSets[graph_?GraphQ] :=
    Block[{ig = igMake[graph]},
      igVertexNames[graph] /@ igIndexVec@ig@"maximalIndependentVertexSets"[]
    ]

IGIndependenceNumber[graph_?GraphQ] := Block[{ig = igMake[graph]}, ig@"independenceNumber"[]]

End[] (* `Private` *)

EndPackage[]