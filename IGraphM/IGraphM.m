(* Mathematica Package  *)
(* Created by Mathematica Plugin for IntelliJ IDEA, see http://wlplugin.halirutan.de/ *)

(* :Title: IGraph/M   *)
(* :Context: IGraphM` *)
(* :Author: szhorvat  *)
(* :Date: 2015-08-28  *)

(* :Package Version: 0.1.1dev *)
(* :Mathematica Version: 10.0 *)
(* :Copyright: (c) 2015 Szabolcs Horvát *)
(* :Keywords: igraph, graphs, networks, LibraryLink *)
(* :Discussion: igraph interface for Mathematica, see http://igraph.org/ *)

BeginPackage["IGraphM`"];

(* Privately load and configure LTemplate *)
Get["LTemplate`LTemplatePrivate`"];
ConfigureLTemplate["MessageSymbol" -> IGraphM];

(***** Usage mesages *****)

IGraphM::usage = "IGraphM is a symbol to which igraph related messages are associated.";

`Developer`Recompile::usage = "IGraphM`Developer`Recompile[] recompiles the IGraphM library and reloads the functions.";
PrependTo[$ContextPath, $Context <> "Developer`"];

IGDocumentation::usage = "IGDocumentation[] opens IGraph/M documentation.";

IGData::usage =
    "IGData[] returns a list of available items.\n" <>
    "IGData[item] returns the requested item.";

IGVersion::usage = "IGVersion[] returns the version of the igraph library in use.";
IGSeedRandom::usage = "IGSeedRandom[seed] seeds the random number generator used by igraph.";

IGLCF::usage = "IGLCF[shifts, repeats, vertexCount] creates a graph from LCF notation.";

IGBetweenness::usage = "IGBetweenness[graph] gives a list of betweenness centralities for the vertices of graph. Weighted graphs are supported.";
IGEdgeBetweenness::usage = "IGEdgeBetweenness[graph] gives a list of betweenness centralities for the edges of graph. Weighted graphs are supported.";
IGCloseness::usage = "IGCloseness[graph, options] gives a list of closeness centralities for the vertices of graph. Weighted graphs are supported.";

IGBetweennessEstimate::usage = "IGBetweennessEstimate[graph, cutoff] estimates vertex betweenness by consdering only paths of at most length cutoff.";
IGEdgeBetweennessEstimate::usage = "IGEdgeBetweennessEstimate[graph, cutoff] estimates edge betweenness by consdering only paths of at most length cutoff.";
IGClosenessEstimate::usage = "IGClosenessEstimate[graph, cutoff, options] estimates closeness centrality by consdering only paths of at most length cutoff.";

IGRewire::usage = "IGRewire[graph, n, options] attempts to rewire the edges of graph n times while presernving its degree sequence.";
IGRewireEdges::usage = "IGRewireEdges[graph, p, options] rewires each edge of the graph with probability p.";

IGDirectedAcyclicGraphQ::usage = "IGDirectedAcyclicGraphQ[graph] tests if graph is directed and acyclic.";
IGConnectedQ::usage = "IGConnectedQ[graph] tests if graph is strongly connected.";
IGWeaklyConnectedQ::usage = "IGWeaklyConnectedQ[graph] tests if graph is weakly connected.";
IGGraphicalQ::usage =
    "IGGraphicalQ[degrees] tests if a degree sequence for an undirected simple graph is graphical.\n" <>
    "IGGraphicalQ[indegrees, outdegrees] tests if a degree sequence for a directed simple graph is graphical.";

IGIsomorphicQ::usage = "IGIsomorphicQ[graph1, graph2] tests if graph1 and graph2 are isomorphic.";
IGSubisomorphicQ::usage = "IGSubisomorphicQ[graph, subgraph] tests if subgraph is contained within graph.";
IGIsoclass::usage = "IGIsoclass[graph] returns the isomorphism class of the graph. Used as the index into the vector returned by motif finding functions. See IGData[] to get list of graphs ordered by isoclass.";

IGBlissCanonicalPermutation::usage =
    "IGBlissCanonicalPermutation[graph, options] computes a canonical permutation of the graph vertices. " <>
    "Two graphs are isomorphic iff they have the same canonical permutation.";
IGBlissIsomorphicQ::usage = "IGBlissIsomorphicQ[graph1, graph2, options] tests if graph1 and graph2 are ismorphic using the BLISS algorithm.";
IGBlissGetIsomorphism::usage = "IGBlissGetIsomorphism[graph1, graph2, options] returns one isomorphism between graph1 and graph2, if it exists.";
IGBlissAutomorphismCount::usage = "IGBlissAutomorphismCount[graph] returns the number of automorphisms of graph.";

IGVF2IsomorphicQ::usage = "IGVF2IsomorphicQ[graph1, graph2, options] tests if graph1 and graph2 are ismorphic using the VF2 algorithm.";
IGVF2FindIsomorphisms::usage =
    "IGVF2FindIsomorphism[graph1, graph2, options] finds all isomorphisms between graph1 and graph2 using the VF2 algorithm." <>
    "IGVF2FindIsomorphism[graph1, graph2, n, options] finds at most n isomorphisms between graph1 and graph2.";
IGVF2SubisomorphicQ::usage = "IGVF2SubisomorphicQ[graph, subgraph, options] tests if subgraph is contained in graph using the VF2 algorithm.";
IGVF2FindSubisomorphisms::usage =
    "IGVF2FindSubisomorphism[graph, subgraph, options] finds all subisomorphisms from subgraph to graph using the VF2 algorithm." <>
    "IGVF2FindSubisomorphism[graph, subgraph, n, options] finds at most n subisomorphisms from subgraph to graph.";
IGVF2AutomorphismCount::usage = "IGVF2AutomorphismCount[graph] returns the number of automorphisms of graph.";
IGVF2IsomorphismCount::usage =
    "IGVF2IsomorphismCount[graph1, graph2, options] returns the number of isomorphisms between graph1 and graph2." <>
    "Note that this is not the same as simply counting the automorphisms of one graph if their vertex or edge colorings differ.";
IGVF2SubisomorphismCount::usage = "IGVF2SubisomorphismCount[subgraph, graph, options]";

IGLADSubisomorphicQ::usage = "IGLADSubisomorphicQ[subgraph, graph] tests if subgraph is contained in graph. Use the \"Induced\" -> True option to look for induced subgraphs.";
IGLADGetSubisomorphism::usage = "IGLADGetSubisomorphism[subgraph, graph] returns one isomorphism between graph1 and graph2, if it exists.";

IGTopologicalOrdering::usage =
    "IGTopologicalOrdering[graph] returns a permutation that sorts the vertices in topological order." <>
    "Note that the values returned are vertex indices, not vertex names.";
IGFeedbackArcSet::usage = "IGFeedbackArcSet[graph]";

IGDyadCensus::usage = "IGDyadCensus[graph]";
IGTriadCensus::usage = "IGTriadCensus[graph]";
IGMotifs::usage = "IGMotifs[graph, motifSize] returns the motif distribution of graph. See IGIsoclass and IGData for motif ordering.";
IGMotifsTotalCount::usage = "IGMotifsTotalCount[graph, motifSize]";
IGMotifsEstimateTotalCount::usage = "IGMotifsEstimate[graph, motifSize, sampleSize]";

IGDegreeSequenceGame::usage =
    "IGDegreeSequenceGame[degrees, options] generates an undirected random graph with the given degree sequence.\n" <>
    "IGDegreeSequenceGame[indegrees, outdegrees, options] generates a directed random graph with the given in- and out-degree sequences.";

IGKRegularGame::usage = "IGKRegularGame[n, k]";

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

IGLayoutRandom::usage = "IGLayoutRandom[graph]";
IGLayoutCircle::usage = "IGLayoutCircle[graph]";
IGLayoutSphere::usage = "IGLayoutSphere[graph]";
IGLayoutGraphOpt::usage = "IGLayoutGraphOpt[graph, options]";
IGLayoutKamadaKawai::usage = "IGLayoutKamadaKawai[graph, options]";
IGLayoutKamadaKawai3D::usage = "IGLayoutKamadaKawai3D[graph, options]";
IGLayoutFruchtermanReingold::usage = "IGLayoutFruchtermanReingold[graph, options]";
IGLayoutFruchtermanReingold3D::usage = "IGLayoutFruchtermanReingold3D[graph, options]";

IGGlobalClusteringCoefficient::usage = "IGGlobalClusteringCoefficient[graph]";
IGLocalClusteringCoefficient::usage = "IGLocalClusteringCoefficient[graph]";
IGAverageLocalClusteringCoefficient::usage = "IGAverageLocalClusteringCoefficient[graph]";
IGWeightedClusteringCoefficient::usage = "IGWeightedClusteringCoefficient[graph] computes the weighted local clustering coefficient, as defined by A. Barrat et al. (2004) http://arxiv.org/abs/cond-mat/0311416";

IGCocitationSimilarity::usage =
    "IGCocitationSimilarity[graph]\n" <>
    "IGCocitationSimilarity[graph, vertex]\n" <>
    "IGCocitationSimilarity[graph, {vertex1, vertex2, \[Ellipsis]}]";

IGBibliographicCoupling::usage =
    "IGBibliographicCoupling[graph]\n" <>
    "IGBibliographicCoupling[graph, vertex]\n" <>
    "IGBibliographicCoupling[graph, {vertex1, vertex2, \[Ellipsis]}]";

IGJaccardSimilarity::usage =
    "IGJaccardSimilarity[graph]\n" <>
    "IGJaccardSimilarity[graph, {vertex1, vertex2, \[Ellipsis]}]";

IGDiceSimilarity::usage =
    "IGDiceSimilarity[graph]\n" <>
    "IGDiceSimilarity[graph, {vertex1, vertex2, \[Ellipsis]}]";

IGInverseLogWeightedSimilarity::usage =
    "IGInverseLogWeightedSimilarity[graph]\n" <>
    "IGInverseLogWeightedSimilarity[graph, vertex]\n" <>
    "IGInverseLogWeightedSimilarity[graph, {vertex1, vertex2, \[Ellipsis]}]";

IGMaximumCardinalitySearch::usage = "IGMaximumCardinalitySearch[graph]";
IGChordalQ::usage = "IGChordalQ[graph]";
IGChordalCompletion::usage = "IGChordalCompletion[graph]";

IGMinSeparators::usage = "IGMinSeparators[graph]";

IGArticulationPoints::usage = "IGArticulationPoints[graph] finds the articulation points of graph. A vertex is an articulation point if its removal increases the number of connected components in the graph.";

IGBiconnectedComponents::usage = "IGBiconnectedComponents[graph]";

Begin["`Private`"];

(***** Mathematica version check *****)

(* Abort loading and leave a clean $ContextPath behind *)
packageAbort[] := (End[]; EndPackage[]; Abort[])

If[Not@OrderedQ[{10.0, 2}, {$VersionNumber, $ReleaseNumber}],
  Print["IGraph/M requires Mathematica 10.0.2 or later.  Aborting."];
  packageAbort[]
]


(***** Package variables *****)

$packageVersion    = "0.1.1dev";
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
        LFun["compilationDate", {}, "UTF8String"],

        (* Graph related functions that do not use the graph data structure *)

        LFun["graphicalQ", {{Real, 1, "Constant"} (* outdeg *), {Real, 1, "Constant"} (* indeg *)}, True|False]
      }
    ],
    LClass["IG",
      {
        (* Create *)

        LFun["fromEdgeList", {{Real, 2, "Constant"} (* edges *), Integer (* vertex count *), True|False (* directed *)}, "Void"],
        LFun["fromLCF", {Integer, {Real, 1, "Constant"}, Integer}, "Void"],

        (* Weights *)

        LFun["setWeights", {{Real, 1, "Constant"}}, "Void"],
        LFun["getWeights", {}, {Real, 1}],
        LFun["clearWeights", {}, "Void"],
        LFun["weightedQ", {}, True|False],

        (* Games *)

        LFun["degreeSequenceGame", {{Real, 1, "Constant"} (* outdeg *), {Real, 1, "Constant"} (* indeg *), Integer (* method *)}, "Void"],
        LFun["kRegularGame", {Integer, Integer, True|False (* directed *), True|False (* multiple *)}, "Void"],

        (* Structure *)

        LFun["edgeCount", {}, Integer],
        LFun["vertexCount", {}, Integer],

        LFun["edgeList", {}, {Real, 2}],

        (* Testing *)

        LFun["directedQ", {}, True|False],
        LFun["dagQ", {}, True|False],
        LFun["simpleQ", {}, True|False],
        LFun["connectedQ", {True|False (* strongly connected *)}, True|False],

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
        LFun["blissAutomorphismCount", LinkObject],
        LFun["vf2Isomorphic", {LExpressionID["IG"], {Integer, 1, "Constant"}, {Integer, 1, "Constant"}, {Integer, 1, "Constant"}, {Integer, 1, "Constant"}}, True|False],
        LFun["vf2FindIsomorphisms", LinkObject],
        LFun["vf2Subisomorphic", {LExpressionID["IG"], {Integer, 1, "Constant"}, {Integer, 1, "Constant"}, {Integer, 1, "Constant"}, {Integer, 1, "Constant"}}, True|False],
        LFun["vf2FindSubisomorphisms", LinkObject],
        LFun["vf2AutomorphismCount", {}, Integer],
        LFun["vf2IsomorphismCount", {LExpressionID["IG"], {Integer, 1, "Constant"}, {Integer, 1, "Constant"}, {Integer, 1, "Constant"}, {Integer, 1, "Constant"}}, Integer],
        LFun["vf2SubisomorphismCount", {LExpressionID["IG"], {Integer, 1, "Constant"}, {Integer, 1, "Constant"}, {Integer, 1, "Constant"}, {Integer, 1, "Constant"}}, Integer],
        LFun["ladSubisomorphic", {LExpressionID["IG"], True|False (* induced *)}, True|False],
        LFun["ladGetSubisomorphism", {LExpressionID["IG"], True|False (* induced *)}, {Real, 1}],

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
        LFun["independenceNumber", {}, Integer],

        (* Graph drawing (layouts) *)

        LFun["layoutRandom", {}, {Real, 2}],
        LFun["layoutCircle", {}, {Real, 2}],
        LFun["layoutSphere", {}, {Real, 2}],

        LFun["layoutGraphOpt",
          {{Real, 2, "Constant"} (* initial position *), True|False (* use initial *),
            Integer (* niter *),
            Real (* charge *), Real (* mass *), Real (* spring length *),
            Real (* spring constant *), Real (* max sa movement *)},
          {Real, 2}
        ],

        LFun["layoutKamadaKawai",
          {{Real, 2, "Constant"} (* initial position *), True|False (* use initial *),
            Integer (* maxiter *), Real (* epsilon *), Real (* kkconst *)},
          {Real, 2}
        ],

        LFun["layoutKamadaKawai3D",
          {{Real, 2, "Constant"} (* initial position *), True|False (* use initial *),
            Integer (* maxiter *), Real (* epsilon *), Real (* kkconst *)},
          {Real, 2}
        ],

        LFun["layoutFruchtermanReingold",
          {{Real, 2, "Constant"} (* initial position *), True|False (* use initial *),
            Integer (* niter *), Real (* start_temp *), Integer (* grid method *)},
          {Real, 2}
        ],

        LFun["layoutFruchtermanReingold3D",
          {{Real, 2, "Constant"} (* initial position *), True|False (* use initial *),
            Integer (* niter *), Real (* start_temp *)},
          {Real, 2}
        ],

        (* Clustering coefficient *)

        LFun["transitivityUndirected", {}, Real],
        LFun["transitivityLocalUndirected", {}, {Real, 1}],
        LFun["transitivityAverageLocalUndirected", {}, Real],
        LFun["transitivityBarrat", {}, {Real, 1}],

        (* Similarity *)

        LFun["similarityCocitation", {{Real, 1, "Constant"}}, {Real, 2}],
        LFun["similarityBibcoupling", {{Real, 1, "Constant"}}, {Real, 2}],
        LFun["similarityJaccard", {{Real, 1, "Constant"}, True|False (* self loops *)}, {Real, 2}],
        LFun["similarityDice", {{Real, 1, "Constant"}, True|False (* self loops *)}, {Real, 2}],
        LFun["similarityInverseLogWeighted", {{Real, 1, "Constant"}}, {Real, 2}],

        (* Chordal graphs *)

        LFun["maximumCardinalitySearch", {}, {Real, 1}],
        LFun["chordalQ", {}, True|False],
        LFun["chordalCompletion", {LExpressionID["IG"]}, {Real, 1}],

        (* Vertex separators *)

        LFun["minimumSizeSeparators", LinkObject],

        (* Articulation points *)

        LFun["articulationPoints", {}, {Real, 1}],
        LFun["biconnectedComponents", LinkObject]
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


Recompile::build = "No build settings found. Please check BuildSettings.m."

Recompile[] :=
    Module[{},
      If[$buildSettings === None,
        Message[Recompile::build];
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
  If[Recompile[] === $Failed
    ,
    Print[Style["Cannot load or compile library. \[FreakedSmiley] Aborting.", Red]];
    packageAbort[]
    ,
    Print[Style["Successfully compiled and loaded the library. \[HappySmiley]", Red]];
  ]
]


(***** General messages *****)

IGraphM::mixed = "Mixed graphs are not supported by IGraph/M.";
IGraphM::vf2col = "Vertex or edge color specifications for VF2 functions must be a pair of integer lists or None.";
IGraphM::lytcrd = "The graph doesn't already have existing vertex coordinates. The \"Continue\" -> True layout option will be ignored.";
IGraphM::lytdim = "The existing vertex coordinates do not have the appropriate dimension for this layout algorithm. The \"Continue\" -> True layout option will be ignored.";
IGraphM::lytcnt = "`` is not a valid value for the \"Continue\" layout option.";

(***** Helper functions *****)

nonNegIntVecQ = VectorQ[#, Internal`NonNegativeMachineIntegerQ]&
intVecQ = VectorQ[#, Developer`MachineIntegerQ]&

(* Zero out the diagonal of a square matrix. *)
zeroDiagonal[arg_] := UpperTriangularize[arg, 1] + LowerTriangularize[arg, -1]

(* Replace Infinity by 0 *)
infToZero[arg_] := Replace[arg, Infinity -> 0]

(* Return immediately from parent function if we have a LibraryFunctionError.
   Used on the immediately on return values of library functions.
   Requires a parent funtcion that uses Block. *)
check = If[MatchQ[#, _LibraryFunctionError], Return[#, Block], #]&

(* Import compressed expressions. Used in IGData. *)
zimport[filename_] := Uncompress@Import[filename, "String"]

(* Get an IG compatible edge list. *)
igEdgeList[graph_] :=
    Developer`ToPackedArray@N[List @@@ EdgeList[graph] /.
            Dispatch@Thread[VertexList[graph] -> Range@VertexCount[graph] - 1]]

(* Convert IG format vertex or edge index vector to Mathematica format. *)
igIndexVec[expr_LibraryFunctionError] := expr (* hack: allows LibraryFunctionError to fall through *)
igIndexVec[arr_] := 1 + Round[arr]

igDirectedQ[graph_] := DirectedGraphQ[graph] && Not@EmptyGraphQ[graph]

(* TODO: Find out how to implement this in a more robust way.
   We only want edge-weighted graphs, not vertex weighted ones. *)
igWeightedGraphQ = WeightedGraphQ[#] && PropertyValue[#, EdgeList] =!= Automatic &;

(* Create IG object from Mathematica Graph. *)
igMake[g_] :=
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


(* check if the argument is an igraph compatible graph *)
igGraphQ = GraphQ[#] && If[MixedGraphQ[#], Message[IGraphM::mixed]; False, True] &

(***** Public functions *****)

IGDocumentation[] := (NotebookOpen@FileNameJoin[{$packageDirectory, "Documentation", "IGDocumentation.nb"}, Saveable -> False]; Null)

(*  IGData  *)

$igData = zimport@FileNameJoin[{$packageDirectory, "IGData.mz"}];
IGData[] := Keys[$igData]
IGData[item_] := Lookup[$igData, Key[item], Missing["NotAvailable"]]

(* General (global) *)

IGVersion[] :=
    "IGraph/M " <> $packageVersion <>
    "\nigraph " <> igraphGlobal@"version"[] <> " (" <> igraphGlobal@"compilationDate"[] <> ")\n" <>
    $System;

IGSeedRandom[seed_?Internal`NonNegativeMachineIntegerQ] := igraphGlobal@"seedRandom"[seed]

(* Create *)

Options[IGLCF] = { GraphLayout -> "CircularEmbedding" };

IGLCF[shifts_?intVecQ, repeats : _?Internal`PositiveMachineIntegerQ : 1, n : (_?Internal`PositiveMachineIntegerQ | Automatic) : Automatic, opt : OptionsPattern[{IGLCF, Graph}]] :=
    Block[{ig = Make["IG"]},
      ig@"fromLCF"[Replace[n, Automatic :> Length[shifts] repeats], shifts, repeats];
      Graph[igToGraph[ig], GraphLayout -> OptionValue[GraphLayout], opt]
    ]

(* Create (games) *)

Options[IGDegreeSequenceGame] = { Method -> "SimpleNoMultiple" };

igDegreeSequenceGameMethods = <| "VigerLatapy" -> 2, "SimpleNoMultiple" -> 1, "Simple" -> 0 |>

IGDegreeSequenceGame::usage = IGDegreeSequenceGame::usage <>
    StringTemplate[" Available methods: ``"][ToString@InputForm@Keys[igDegreeSequenceGameMethods]];

IGDegreeSequenceGame[degrees_?nonNegIntVecQ, opt : OptionsPattern[]] := IGDegreeSequenceGame[{}, degrees, opt]

IGDegreeSequenceGame[indegrees_?nonNegIntVecQ, outdegrees_?nonNegIntVecQ, opt : OptionsPattern[{IGDegreeSequenceGame, Graph}]] :=
    Block[{ig = Make["IG"]},
      Check[
        ig@"degreeSequenceGame"[outdegrees, indegrees, Lookup[igDegreeSequenceGameMethods, OptionValue[Method], -1]];
        Graph[igToGraph[ig], Sequence@@FilterRules[{opt}, Options[Graph]]]
        ,
        $Failed
        ,
        {LTemplate::error}
      ]
    ]


Options[IGKRegularGame] = { "AllowMultipleEdges" -> False, DirectedEdges -> False };
IGKRegularGame[n_?Internal`PositiveMachineIntegerQ, k_?Internal`PositiveMachineIntegerQ, opt : OptionsPattern[{IGKRegularGame, Graph}]] :=
    Block[{ig = Make["IG"]},
      ig@"kRegularGame"[n, k, OptionValue[DirectedEdges], OptionValue["AllowMultipleEdges"]];
      Graph[igToGraph[ig], Sequence@@FilterRules[{opt}, Options[Graph]]]
    ]

(* Testing *)

IGDirectedAcyclicGraphQ[g_?igGraphQ] := Block[{ig = igMake[g]}, ig@"dagQ"[]]

IGConnectedQ[g_?igGraphQ] := Block[{ig = igMake[g]}, ig@"connectedQ"[True]]
IGWeaklyConnectedQ[g_?igGraphQ] := Block[{ig = igMake[g]}, ig@"connectedQ"[False]]

IGGraphicalQ[degrees_?nonNegIntVecQ] := IGGraphicalQ[{}, degrees]
IGGraphicalQ[indeg_?nonNegIntVecQ, outdeg_?nonNegIntVecQ] := igraphGlobal@"graphicalQ"[outdeg, indeg]

(* Centrality *)

IGBetweenness[g_?igGraphQ] := Block[{ig = igMake[g]}, ig@"betweenness"[]]

IGEdgeBetweenness[g_?igGraphQ] := Block[{ig = igMake[g]}, ig@"edgeBetweenness"[]]

Options[IGCloseness] = { "Normalized" -> False };
IGCloseness[g_?igGraphQ, opt : OptionsPattern[]] := Block[{ig = igMake[g]}, ig@"closeness"[OptionValue["Normalized"]]]

(* Centrality estimates *)

IGBetweennessEstimate[g_?igGraphQ, cutoff_?Positive] := Block[{ig = igMake[g]}, ig@"betweennessEstimate"@infToZero[cutoff]]

IGEdgeBetweennessEstimate[g_?igGraphQ, cutoff_?Positive] := Block[{ig = igMake[g]}, ig@"edgeBetweennessEstimate"@infToZero[cutoff]]

Options[IGClosenessEstimate] = { "Normalized" -> False };
IGClosenessEstimate[g_?igGraphQ, cutoff_?Positive, opt : OptionsPattern[]] :=
    Block[{ig = igMake[g]},
      ig@"closenessEstimate"[infToZero[cutoff], OptionValue["Normalized"]]
    ]

(* Randomization *)

(* TODO: functions in this section should warn that edge weights will be lost *)

Options[IGRewire] = { "AllowLoops" -> False };
IGRewire[g_?igGraphQ, n_?Internal`PositiveMachineIntegerQ, opt : OptionsPattern[]] :=
    Block[{ig = igMake[g]},
      ig@"rewire"[n, OptionValue["AllowLoops"]];
      igToGraph[ig]
    ]

Options[IGRewireEdges] = { "AllowLoops" -> False, "AllowMultipleEdges" -> False };
IGRewireEdges[g_?igGraphQ, p_?Internal`RealValuedNumericQ, opt : OptionsPattern[]] :=
    Block[{ig = igMake[g]},
      ig@"rewireEdges"[p, OptionValue["AllowLoops"], OptionValue["AllowMultipleEdges"]];
      igToGraph[ig]
    ]

(* Isomorphism *)

IGIsomorphicQ[g1_?igGraphQ, g2_?igGraphQ] :=
    Block[{ig1 = igMake[g1], ig2 = igMake[g2]}, ig1@"isomorphic"[ManagedLibraryExpressionID@ig2]]

IGSubisomorphicQ[subgraph_?igGraphQ, graph_?igGraphQ] :=
    Block[{ig1 = igMake[graph], ig2 = igMake[subgraph]}, ig1@"subisomorphic"[ManagedLibraryExpressionID@ig2]]

IGIsoclass[graph_?igGraphQ] := Block[{ig = igMake[graph]}, ig@"isoclass"[]]


blissSplittingHeuristicsNames = {
  "First", "FirstSmallest", "FirstLargest",
  "FirstMaximallyConnected", "FirstSmallestMaximallyConnected", "FirstLargestMaximallyConnected"
};

blissSplittingHeuristics = AssociationThread[blissSplittingHeuristicsNames, Range@Length[blissSplittingHeuristicsNames] - 1];

IGBlissCanonicalPermutation::usage = IGBlissCanonicalPermutation::usage <>
    StringTemplate[" Available values for the \"SplittingHeuristics\" option: ``." <>
        "The permutation depends on the splitting heuristics used."][ToString@InputForm@blissSplittingHeuristicsNames];

Options[IGBlissCanonicalPermutation] = { "SplittingHeuristics" -> "First" };
IGBlissCanonicalPermutation[graph_?igGraphQ, opt : OptionsPattern[]] :=
    Block[{ig = igMake[graph]},
      igVertexNames[graph]@igIndexVec@ig@"blissCanonicalPermutation"[Lookup[blissSplittingHeuristics, OptionValue["SplittingHeuristics"], -1]]
    ]

Options[IGBlissIsomorphicQ] = { "SplittingHeuristics" -> "First" };
IGBlissIsomorphicQ[graph1_?igGraphQ, graph2_?igGraphQ, opt : OptionsPattern[]] :=
    Block[{ig1 = igMake[graph1], ig2 = igMake[graph2]},
      ig1@"blissIsomorphic"[ManagedLibraryExpressionID[ig2], Lookup[blissSplittingHeuristics, OptionValue["SplittingHeuristics"], -1]]
    ]

Options[IGBlissGetIsomorphism] = { "SplittingHeuristics" -> "First" };
IGBlissGetIsomorphism[graph1_?igGraphQ, graph2_?igGraphQ, opt : OptionsPattern[]] :=
    Block[{ig1 = igMake[graph1], ig2 = igMake[graph2], result},
      result = igIndexVec@check@ig1@"blissFindIsomorphism"[ManagedLibraryExpressionID[ig2], Lookup[blissSplittingHeuristics, OptionValue["SplittingHeuristics"], -1]];
      If[result === {}, Return[{}]];
      List@AssociationThread[
        VertexList[graph1],
        igVertexNames[graph2]@result
      ]
    ]

Options[IGBlissAutomorphismCount] = { "SplittingHeuristics" -> "First" };
IGBlissAutomorphismCount[graph_?igGraphQ, opt : OptionsPattern[]] :=
    Block[{ig = igMake[graph]},
      ToExpression@check@ig@"blissAutomorphismCount"[Lookup[blissSplittingHeuristics, OptionValue["SplittingHeuristics"], -1]]
    ]


vf2ParseColors[None] := {{},{}}
vf2ParseColors[col : {_?intVecQ, _?intVecQ}] := col
vf2ParseColors[expr_] := (Message[IGraphM::vf2col]; {{},{}})

Options[IGVF2IsomorphicQ] = { "VertexColors" -> None, "EdgeColors" -> None };

IGVF2IsomorphicQ[graph1_?igGraphQ, graph2_?igGraphQ, opt : OptionsPattern[]] :=
    Block[{ig1 = igMake[graph1], ig2 = igMake[graph2], vcol1, vcol2, ecol1, ecol2},
      {vcol1, vcol2} = vf2ParseColors@OptionValue["VertexColors"];
      {ecol1, ecol2} = vf2ParseColors@OptionValue["EdgeColors"];
      ig1@"vf2Isomorphic"[ManagedLibraryExpressionID[ig2], vcol1, vcol2, ecol1, ecol2]
    ]

Options[IGVF2FindIsomorphisms] = { "VertexColors" -> None, "EdgeColors" -> None };

IGVF2FindIsomorphisms[graph1_?igGraphQ, graph2_?igGraphQ, max : (_?Internal`PositiveMachineIntegerQ | All | Infinity) : All, opt : OptionsPattern[]] :=
    Block[{ig1 = igMake[graph1], ig2 = igMake[graph2], vcol1, vcol2, ecol1, ecol2, n, result},
      n = Replace[max, All|Infinity -> -1];
      {vcol1, vcol2} = vf2ParseColors@OptionValue["VertexColors"];
      {ecol1, ecol2} = vf2ParseColors@OptionValue["EdgeColors"];
      result = igIndexVec@check@ig1@"vf2FindIsomorphisms"[ManagedLibraryExpressionID[ig2], n, vcol1, vcol2, ecol1, ecol2];
      AssociationThread[
        VertexList[graph1],
        igVertexNames[graph2][#]
      ]& /@ result
    ]

Options[IGVF2SubisomorphicQ] = { "VertexColors" -> None, "EdgeColors" -> None };

IGVF2SubisomorphicQ[subgraph_?igGraphQ, graph_?igGraphQ, opt : OptionsPattern[]] :=
    Block[{ig1 = igMake[graph], ig2 = igMake[subgraph], vcol1, vcol2, ecol1, ecol2},
      {vcol1, vcol2} = Reverse@vf2ParseColors@OptionValue["VertexColors"];
      {ecol1, ecol2} = Reverse@vf2ParseColors@OptionValue["EdgeColors"];
      ig1@"vf2Subisomorphic"[ManagedLibraryExpressionID[ig2], vcol1, vcol2, ecol1, ecol2]
    ]

Options[IGVF2FindSubisomorphisms] = { "VertexColors" -> None, "EdgeColors" -> None };

IGVF2FindSubisomorphisms[subgraph_?igGraphQ, graph_?igGraphQ, max : (_?Internal`PositiveMachineIntegerQ | All | Infinity) : All, opt : OptionsPattern[]] :=
    Block[{ig1 = igMake[graph], ig2 = igMake[subgraph], vcol1, vcol2, ecol1, ecol2, n, result},
      n = Replace[max, All|Infinity -> -1];
      {vcol1, vcol2} = Reverse@vf2ParseColors@OptionValue["VertexColors"];
      {ecol1, ecol2} = Reverse@vf2ParseColors@OptionValue["EdgeColors"];
      result = igIndexVec@check@ig1@"vf2FindSubisomorphisms"[ManagedLibraryExpressionID[ig2], n, vcol1, vcol2, ecol1, ecol2];
      AssociationThread[
        VertexList[subgraph],
        igVertexNames[graph][#]
      ]& /@ result
    ]

IGVF2AutomorphismCount[graph_?igGraphQ] :=
    Block[{ig = igMake[graph]}, ig@"vf2AutomorphismCount"[]]


Options[IGVF2IsomorphismCount] = { "VertexColors" -> None, "EdgeColors" -> None };

IGVF2IsomorphismCount[graph1_?igGraphQ, graph2_?igGraphQ, opt : OptionsPattern[]] :=
    Block[{ig1 = igMake[graph1], ig2 = igMake[graph2], vcol1, vcol2, ecol1, ecol2},
      {vcol1, vcol2} = vf2ParseColors@OptionValue["VertexColors"];
      {ecol1, ecol2} = vf2ParseColors@OptionValue["EdgeColors"];
      ig1@"vf2IsomorphismCount"[ManagedLibraryExpressionID[ig2], vcol1, vcol2, ecol1, ecol2]
    ]

Options[IGVF2SubisomorphismCount] = { "VertexColors" -> None, "EdgeColors" -> None };

IGVF2SubisomorphismCount[subgraph_?igGraphQ, graph_?igGraphQ, opt : OptionsPattern[]] :=
    Block[{ig1 = igMake[graph], ig2 = igMake[subgraph], vcol1, vcol2, ecol1, ecol2},
      {vcol1, vcol2} = Reverse@vf2ParseColors@OptionValue["VertexColors"];
      {ecol1, ecol2} = Reverse@vf2ParseColors@OptionValue["EdgeColors"];
      ig1@"vf2SubisomorphismCount"[ManagedLibraryExpressionID[ig2], vcol1, vcol2, ecol1, ecol2]
    ]


Options[IGLADSubisomorphicQ] = { "Induced" -> False };

IGLADSubisomorphicQ[subgraph_?igGraphQ, graph_?igGraphQ, opt : OptionsPattern[]] :=
    Block[{ig1 = igMake[graph], ig2 = igMake[subgraph]},
      ig1@"ladSubisomorphic"[ManagedLibraryExpressionID[ig2], OptionValue["Induced"]]
    ]

Options[IGLADGetSubisomorphism] = { "Induced" -> False };

IGLADGetSubisomorphism[subgraph_?igGraphQ, graph_?igGraphQ, opt : OptionsPattern[]] :=
    Block[{ig1 = igMake[graph], ig2 = igMake[subgraph], result},
      result = igIndexVec@check@ig1@"ladGetSubisomorphism"[ManagedLibraryExpressionID[ig2], OptionValue["Induced"]];
      If[result === {}, Return[{}]];
      List@AssociationThread[
        VertexList[subgraph],
        igVertexNames[graph][result]
      ]
    ]

(* Directed acylic graphs and topological ordering *)

IGTopologicalOrdering[graph_?igGraphQ] := Block[{ig = igMake[graph]}, igIndexVec@ig@"topologicalSorting"[]]

Options[IGFeedbackArcSet] = { "Exact" -> True };
IGFeedbackArcSet[graph_?igGraphQ, opt : OptionsPattern[]] :=
    Block[{ig = igMake[graph]},
      Part[EdgeList[graph], igIndexVec@ig@"feedbackArcSet"[OptionValue["Exact"]]]
    ]

(* Motifs and subgraph counts *)

IGDyadCensus[graph_?igGraphQ] := Block[{ig = igMake[graph]}, AssociationThread[{"Mutual", "Asymmetric", "Missing"}, Round@ig@"dyadCensus"[]]]

IGTriadCensus[graph_?igGraphQ] :=
    Block[{ig = igMake[graph]},
      AssociationThread[
        {"003", "012", "102", "021D", "021U", "021C", "111D", "111U", "030T", "030C", "201", "120D", "120U", "120C", "210", "300"},
        Round@check@ig@"triadCensus"[]
      ]
    ]

IGMotifs[graph_?igGraphQ, size_?Internal`PositiveIntegerQ] :=
    Block[{ig = igMake[graph]},
      Round@Developer`FromPackedArray@ig@"motifs"[size, ConstantArray[0, size]]
    ]

IGMotifsTotalCount[graph_?igGraphQ, size_?Internal`PositiveIntegerQ] :=
    Block[{ig = igMake[graph]}, ig@"motifsNo"[size, ConstantArray[0, size]] ]

IGMotifsEstimateTotalCount[graph_?igGraphQ, size_?Internal`PositiveIntegerQ, sampleSize_?Internal`PositiveIntegerQ] :=
    Block[{ig = igMake[graph]}, ig@"motifsEstimate"[size, ConstantArray[0, size], sampleSize] ]

(* Shortest paths *)

IGDistanceMatrix[graph_?igGraphQ] :=
    Block[{ig = igMake[graph]},
      zeroDiagonal[Round[ig@"shortestPaths"[]] /. 0 -> Infinity] (* TODO: avoid unpacking when no infinities present *)
    ]

(* Cliques *)

IGCliques[graph_] := IGCliques[graph, Infinity]
IGCliques[graph_, max : (_Integer | Infinity)] := IGCliques[graph, {1, max}]
IGCliques[graph_, {size_}] := IGCliques[graph, {size, size}]
IGCliques[graph_?igGraphQ, {min_?Internal`PositiveMachineIntegerQ, max : (_?Internal`PositiveMachineIntegerQ | Infinity)}] /; max >= min :=
    iIGCliques[graph, {min, infToZero[max]}]
iIGCliques[graph_, {min_, max_}] :=
    Block[{ig = igMake[graph]},
      igVertexNames[graph] /@ igIndexVec@ig@"cliques"[min, max]
    ]

IGMaximalCliques[graph_] := IGMaximalCliques[graph, Infinity]
IGMaximalCliques[graph_, max : (_Integer | Infinity)] := IGMaximalCliques[graph, {1, max}]
IGMaximalCliques[graph_, {size_}] := IGMaximalCliques[graph, {size, size}]
IGMaximalCliques[graph_?igGraphQ, {min_?Internal`PositiveMachineIntegerQ, max : (_?Internal`PositiveMachineIntegerQ | Infinity)}] /; max >= min :=
    iIGMaximalCliques[graph, {min, infToZero[max]}]
iIGMaximalCliques[graph_, {min_, max_}] :=
    Block[{ig = igMake[graph]},
      igVertexNames[graph] /@ igIndexVec@ig@"maximalCliques"[min, max]
    ]

IGMaximalCliquesCount[graph_] := IGMaximalCliquesCount[graph, Infinity]
IGMaximalCliquesCount[graph_, max : (_Integer | Infinity)] := IGMaximalCliquesCount[graph, {1, max}]
IGMaximalCliquesCount[graph_, {size_}] := IGMaximalCliquesCount[graph, {size, size}]
IGMaximalCliquesCount[graph_?igGraphQ, {min_?Internal`PositiveMachineIntegerQ, max : (_?Internal`PositiveMachineIntegerQ | Infinity)}] /; max >= min :=
    iIGMaximalCliquesCount[graph, {min, infToZero[max]}]
iIGMaximalCliquesCount[graph_, {min_, max_}] :=
    Block[{ig = igMake[graph]},
      Round@ig@"maximalCliquesCount"[min, max]
    ]

IGLargestCliques[graph_?igGraphQ] :=
    Block[{ig = igMake[graph]},
      igVertexNames[graph] /@ igIndexVec@ig@"largestCliques"[]
    ]

IGCliqueNumber[graph_?igGraphQ] := Block[{ig = igMake[graph]}, ig@"cliqueNumber"[]]

(* Independent vertex sets *)

IGIndependentVertexSets[graph_] := IGIndependentVertexSets[graph, Infinity]
IGIndependentVertexSets[graph_, max : (_Integer | Infinity)] := IGIndependentVertexSets[graph, {1, max}]
IGIndependentVertexSets[graph_, {size_}] := IGIndependentVertexSets[graph, {size, size}]
IGIndependentVertexSets[graph_?igGraphQ, {min_?Internal`PositiveMachineIntegerQ, max : (_?Internal`PositiveMachineIntegerQ | Infinity)}] /; max >= min :=
    iIGIndependentVertexSets[graph, {min, infToZero[max]}]
iIGIndependentVertexSets[graph_, {min_, max_}] :=
    Block[{ig = igMake[graph]},
      igVertexNames[graph] /@ igIndexVec@ig@"independentVertexSets"[min, max]
    ]

IGLargestIndependentVertexSets[graph_?igGraphQ] :=
    Block[{ig = igMake[graph]},
      igVertexNames[graph] /@ igIndexVec@ig@"largestIndependentVertexSets"[]
    ]

IGMaximalIndependentVertexSets[graph_?igGraphQ] :=
    Block[{ig = igMake[graph]},
      igVertexNames[graph] /@ igIndexVec@ig@"maximalIndependentVertexSets"[]
    ]

IGIndependenceNumber[graph_?igGraphQ] := Block[{ig = igMake[graph]}, ig@"independenceNumber"[]]

(* Graph drawing (layouts *)

getVertexCoords[graph_] :=
    With[{coords = GraphEmbedding[graph]},
      If[MatrixQ[coords],
        coords,
        Message[IGraphM::cntlayout]; {{}}
      ]
    ]

continueLayout[graph_, False, ___] := Sequence[{{}}, False]
continueLayout[graph_, True, dim_ : 2] :=
    Sequence@@Module[{coords},
      coords = getVertexCoords[graph];
      If[Not@MatchQ[Dimensions[coords], {_, dim}],
        Message[IGraphM::lytdim];
        coords = {{}}
      ];
      {coords, coords =!= {{}}}
    ]
continueLayout[graph_, cont_, ___] := ( Message[IGraphM::lytcnt, cont]; continueLayout[graph, False] )
continueLayout3D[graph_, cont_] := continueLayout[graph, cont, 3]

setVertexCoords[g_, coords_] := Graph[g, VertexCoordinates -> Thread[ VertexList[g] -> coords ]]
setVertexCoords3D[g_, coords_] := Graph3D[g, VertexCoordinates -> Thread[ VertexList[g] -> coords ]]

IGLayoutRandom[graph_?igGraphQ] :=
    Block[{ig = igMake[graph]},
      Graph[graph, VertexCoordinates -> check@ig@"layoutRandom"[]]
    ]

IGLayoutCircle[graph_?igGraphQ] :=
    Block[{ig = igMake[graph]},
      setVertexCoords[graph, check@ig@"layoutCircle"[]]
    ]

IGLayoutSphere[graph_?igGraphQ] :=
    Block[{ig = igMake[graph]},
      setVertexCoords3D[graph, check@ig@"layoutSphere"[]]
    ]


Options[IGLayoutGraphOpt] = {
  "MaxIterations" -> 500, "NodeCharge" -> 0.001, "NodeMass" -> 30, "SpringLength" -> 0,
  "SpringConstant" -> 1, "MaxStepMovement" -> 5,
  "Continue" -> False
};

IGLayoutGraphOpt[graph_?igGraphQ, opt : OptionsPattern[]] :=
    Block[{ig = igMake[graph]},
      setVertexCoords[graph,
          0.01 check@ig@"layoutGraphOpt"[continueLayout[graph, OptionValue["Continue"]],
            OptionValue["MaxIterations"], OptionValue["NodeCharge"], OptionValue["NodeMass"],
            OptionValue["SpringLength"], OptionValue["SpringConstant"], OptionValue["MaxStepMovement"]
          ]
      ]
    ]

Options[IGLayoutKamadaKawai] = {
  "MaxIterations" -> Automatic, "Epsilon" -> 0, "KamadaKawaiConstant" -> Automatic,
  "Continue" -> False
};

IGLayoutKamadaKawai[graph_?igGraphQ, opt : OptionsPattern[]] :=
    Block[{ig = igMake[graph], maxiter, kkconst},
      maxiter = Replace[OptionValue["MaxIterations"], Automatic -> 10 VertexCount[graph]];
      kkconst = Replace[OptionValue["KamadaKawaiConstant"], Automatic -> VertexCount[graph]];
      setVertexCoords[graph,
        0.5 check@ig@"layoutKamadaKawai"[continueLayout[graph, OptionValue["Continue"]],
          maxiter, OptionValue["Epsilon"], kkconst]
      ]
    ]

Options[IGLayoutKamadaKawai3D] = {
  "MaxIterations" -> Automatic, "Epsilon" -> 0, "KamadaKawaiConstant" -> Automatic,
  "Continue" -> False
};

IGLayoutKamadaKawai3D[graph_?igGraphQ, opt : OptionsPattern[]] :=
    Block[{ig = igMake[graph], maxiter, kkconst},
      maxiter = Replace[OptionValue["MaxIterations"], Automatic -> 10 VertexCount[graph]];
      kkconst = Replace[OptionValue["KamadaKawaiConstant"], Automatic -> VertexCount[graph]];
      setVertexCoords3D[graph,
        0.5 check@ig@"layoutKamadaKawai3D"[continueLayout3D[graph, OptionValue["Continue"]],
          maxiter, OptionValue["Epsilon"], kkconst]
      ]
    ]

igFruchtermanReingoldMethods = <| Automatic -> 2, False -> 1, True -> 0 |>;

Options[IGLayoutFruchtermanReingold] = {
  "MaxIterations" -> 500, "MaxMovement" -> 5, "UseGrid" -> Automatic,
  "Continue" -> False
};

IGLayoutFruchtermanReingold[graph_?igGraphQ, opt : OptionsPattern[]] :=
    Block[{ig = igMake[graph]},
      setVertexCoords[graph,
        0.25 check@ig@"layoutFruchtermanReingold"[continueLayout[graph, OptionValue["Continue"]],
          OptionValue["MaxIterations"], OptionValue["MaxMovement"], Lookup[igFruchtermanReingoldMethods, OptionValue["UseGrid"], -1]
        ]
      ]
    ]

Options[IGLayoutFruchtermanReingold3D] = {
  "MaxIterations" -> 500, "MaxMovement" -> 5,
  "Continue" -> False
};

IGLayoutFruchtermanReingold3D[graph_?igGraphQ, opt : OptionsPattern[]] :=
    Block[{ig = igMake[graph]},
      setVertexCoords[graph,
        0.25 check@ig@"layoutFruchtermanReingold3D"[continueLayout3D[graph, OptionValue["Continue"]],
          OptionValue["MaxIterations"], OptionValue["MaxMovement"]
        ]
      ]
    ]

(* Clustering coefficient *)

IGGlobalClusteringCoefficient[graph_?igGraphQ] :=
    Block[{ig = igMake[graph]}, ig@"transitivityUndirected"[]]

IGLocalClusteringCoefficient[graph_?igGraphQ] :=
    Block[{ig = igMake[graph]}, ig@"transitivityLocalUndirected"[]]

IGAverageLocalClusteringCoefficient[graph_?igGraphQ] :=
    Block[{ig = igMake[graph]}, ig@"transitivityAverageLocalUndirected"[]]

IGWeightedClusteringCoefficient[graph_?igGraphQ] :=
    Block[{ig = igMake[graph]}, ig@"transitivityBarrat"[]]


(* Similarity *)

(* for those that return a list of vectors *)
similarityFunction1[name_, post_ : Identity][graph_, All] := Block[{ig = igMake[graph]}, post@check@ig@name[{}] ]
similarityFunction1[name_, post_ : Identity][graph_, {}] := {}
similarityFunction1[name_, post_ : Identity][graph_, vs_?ListQ] :=
    Block[{ig = igMake[graph]},
      post@check@ig@name[ Check[VertexIndex[graph, #] - 1& /@ vs, Return[$Failed, Block]] ]
    ]
similarityFunction1[name_, post_ : Identity][graph_, v_] := similarityFunction1[name, First @* post][graph, {v}]

IGCocitationSimilarity[graph_?igGraphQ, vs_ : All] := similarityFunction1["similarityCocitation", Round][graph, vs]

IGBibliographicCoupling[graph_?igGraphQ, vs_ : All] := similarityFunction1["similarityBibcoupling", Round][graph, vs]

IGInverseLogWeightedSimilarity[graph_?igGraphQ, vs_ : All] := similarityFunction1["similarityInverseLogWeighted"][graph, vs]

(* for those that return a matrix *)
similarityFunction2[name_][graph_, All, loops_] := Block[{ig = igMake[graph]}, check@ig@name[{}, loops] ]
similarityFunction2[name_][graph_, {}, loops_] := {}
similarityFunction2[name_][graph_, vs_?ListQ, loops_] :=
    Block[{ig = igMake[graph]},
      check@ig@name[ Check[VertexIndex[graph, #] - 1& /@ vs, Return[$Failed, Block]], loops ]
    ]

Options[IGJaccardSimilarity] = { SelfLoops -> False };
IGJaccardSimilarity[graph_?igGraphQ, vs : (_?ListQ | All) : All, opt : OptionsPattern[]] := similarityFunction2["similarityJaccard"][graph, vs, OptionValue[SelfLoops]]

Options[IGDiceSimilarity] = { SelfLoops -> False };
IGDiceSimilarity[graph_?igGraphQ, vs : (_?ListQ | All) : All, opt : OptionsPattern[]] := similarityFunction2["similarityDice"][graph, vs, OptionValue[SelfLoops]]


(* Chordal graphs *)

IGChordalQ[graph_?igGraphQ] :=
    Block[{ig = igMake[graph]}, ig@"chordalQ"[]]

IGMaximumCardinalitySearch[graph_?igGraphQ] :=
    Block[{ig = igMake[graph]}, igVertexNames[graph]@igIndexVec@ig@"maximumCardinalitySearch"[]]

IGChordalCompletion[graph_?igGraphQ] :=
    Block[{ig = igMake[graph], new = Make["IG"], result},
      result = new@"chordalCompletion"[ManagedLibraryExpressionID[ig]];
      {Partition[igVertexNames[graph]@igIndexVec[result], 2], igToGraph[new]}
    ]

(* Vertex cuts *)

IGMinSeparators[graph_?igGraphQ] := Block[{ig = igMake[graph]}, igVertexNames[graph] /@ igIndexVec@check@ig@"minimumSizeSeparators"[]]

(* Connected components *)

IGArticulationPoints[graph_?igGraphQ] := Block[{ig = igMake[graph]}, igVertexNames[graph]@igIndexVec@check@ig@"articulationPoints"[]]

IGBiconnectedComponents[graph_?igGraphQ] := Block[{ig = igMake[graph]}, igVertexNames[graph] /@ igIndexVec@check@ig@"biconnectedComponents"[]]

End[]; (* `Private` *)

EndPackage[];