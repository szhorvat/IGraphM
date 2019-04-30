(* Mathematica Package  *)
(* Created by Mathematica Plugin for IntelliJ IDEA, see http://wlplugin.halirutan.de/ *)

(* :Title: IGraph/M   *)
(* :Context: IGraphM` *)
(* :Author: szhorvat  *)
(* :Date: 2015-08-28  *)

(* :Package Version: %%version%% *)
(* :Mathematica Version: %%mathversion%% *)
(* :Copyright: (c) 2019 Szabolcs Horvát *)
(* :Keywords: igraph, graphs, networks, LibraryLink *)
(* :Discussion: igraph interface for Mathematica, see http://szhorvat.net/mathematica/IGraphM  *)

Package["IGraphM`"]

(* igContextSetup[] is a workaround for $Context and $ContextPath not being set during loading of new-style packages
   before M11.0. See also the corresponding Block'ing of $Context and $ContextPath in Kernel/init.m.
   Every file that is part of IGraphM` must evaluate igContextSetup[igPackagePrivateSymbol] at the beginning. *)
PackageScope["igContextSetup"]
If[$VersionNumber < 11.0,
  igContextSetup::usage = "igContextSetup[localSymbol]";
  igContextSetup[packagePrivateSymbol_] :=
      ($Context = Context[packagePrivateSymbol];
       $ContextPath = {"IGraphM`PackageScope`", "IGraphM`", "System`"})
]
igContextSetup[igPackagePrivateSymbol]

PackageImport["IGraphM`LTemplate`"]

(***** Usage messages *****)

PackageExport["IGraphM"]
IGraphM::usage = "IGraphM is a symbol to which igraph related messages are associated.";

PackageExport["MultiEdges"]
MultiEdges::usage = "MultiEdges is an option for several IGraph/M functions that specifies whether multi-edges should be generated.";

IGraphM`Information`$Version::usage =
    "IGraphM`Information`$Version is a string that gives the version of the currently loaded IGraph/M package.";

IGraphM`Developer`Recompile::usage =
    "IGraphM`Developer`Recompile[] recompiles the IGraphM library and reloads the functions.";

IGraphM`Developer`GetInfo::usage =
    "IGraphM`Developer`GetInfo[] returns useful information about IGraph/M and the system it is running on, for debugging and troubleshooting purposes.";


(***** Utility functions for setting up definitions and for other load-time tasks *****)

(* These functions must be defined in the file that is loaded first, IGraphM.m, to ensure availability in all others *)

symName[s_Symbol] := SymbolName[s]
symName[s_String] := s

PackageScope["optNames"]
optNames::usage = "optNames[sym1, sym2, ...] returns the option names associated with the given symbols.";
optNames[syms___] :=
    Complement[
      symName /@ Union @@ (Options[#][[All, 1]]& /@ {syms}),
      {"MultipleEdges"} (* these option names are deprecated *)
    ]


PackageScope["amendUsage"]
amendUsage::usage = "amendUsage[symbol, stringTempl, templArg1, templArg2, ...] amends the usage message of symbol.";
amendUsage[sym_Symbol, amend_, args___] :=
    Module[{lines},
      lines = StringSplit[sym::usage, "\n"];
      lines[[1]] = lines[[1]] <> " " <> StringTemplate[amend, InsertionFunction -> (ToString[#, InputForm]&)][args];
      sym::usage = StringJoin@Riffle[lines, "\n"]
    ]


(*
	Numeric codes are for certain special types of completions. Zero means 'don't complete':

	Normal argument     0
	AbsoluteFilename    2
	RelativeFilename    3
	Color               4
	PackageName         7
	DirectoryName       8
	InterpreterType     9
*)

PackageScope["addCompletion"]
addCompletion::usage = "addCompletion[symbol, argSpec] adds FE auto-completion for symbol.";
addCompletion[fun_Symbol, argSpec_List] :=
    If[$Notebooks,
      With[{compl = SymbolName[fun] -> argSpec},
        FE`Evaluate[FEPrivate`AddSpecialArgCompletion[compl]]
      ]
    ]


(* Import compressed expressions. Used in IGData and IGLatticeMesh. *)
(* Avoid Import[] because it doesn't work during kernel initialization. *)
(* Note: Currently, the loading of all external data files is defined to be lazy. However, the data loading is still
 * triggered at package load time in the dev version by some addCompletion[] calls (see IGLatticeMesh[] and IGData[]).
 * Thus zimport[] must be defined in this file (so it is ready by the time other files load) and must not use Import[].
 *)
PackageScope["zimport"]
zimport::usage = "zimport[file] reads a Compress[]'d expression from a text file.";
zimport[filename_] :=
    Module[{stream, str},
      stream = OpenRead[filename];
      str = Read[stream, Record, RecordSeparators -> {}];
      Close[stream];
      Uncompress[str]
    ]


(***** Package variables *****)

(* $packageDirectory must be package scope because it is used in functions that read data files,
   such as IGData[], IGLatticeMesh[] *)
PackageScope["$packageDirectory"]
$packageDirectory::usage = "$packageDirectory is the directory where IGraph/M is installed.";

$packageVersion    = "%%version%% (%%date%%)";
$packageDirectory  = DirectoryName[$InputFileName];
$systemID = $SystemID;
(* On OS X libraries use libc++ ABI since M10.4 and libstdc++ ABI up to M10.3.  We need separate binaries. *)
If[$OperatingSystem === "MacOSX", $systemID = $systemID <> If[$VersionNumber <= 10.3, "-libstdc++", "-libc++"]];
$libraryDirectory  = FileNameJoin[{$packageDirectory, "LibraryResources", $systemID}];
$sourceDirectory   = FileNameJoin[{$packageDirectory, "LibraryResources", "Source"}];
$buildSettingsFile = FileNameJoin[{$packageDirectory, "BuildSettings.m"}];


IGraphM`Information`$Version = $packageVersion;


(***** Class interface definition *****)

template = LTemplate["IGraphM",
  {
    LClass["IGlobal",
      {
        LFun["init", {}, "Void"],
        LFun["seedRandom", {Integer}, "Void"],
        LFun["randomGeneratorName", {}, "UTF8String"],
        LFun["setRandomGenerator", {Integer}, "Void"],
        LFun["version", {}, "UTF8String"],
        LFun["compilationDate", {}, "UTF8String"],

        LFun["infOrNanQ", {{Real, _, "Constant"}}, True|False],

        (* Graph related functions that do not use the IG data structure *)

        LFun["erdosGallai", {{Integer, 1} (* not "Constant" because it gets modified *)}, True|False],
        LFun["graphicalQ", {{Real, 1, "Constant"} (* outdeg *), {Real, 1, "Constant"} (* indeg *)}, True|False],

        LFun["compareCommunities", {{Real, 1, "Constant"}, {Real, 1, "Constant"}, Integer (* method *)}, Real],

        LFun["incidenceToEdgeList", {{LType[SparseArray, Integer], "Constant"}, True|False}, {Integer, 2}],

        LFun["edgeListSortPairs", {{Integer, 2} (* not "Constant" because it gets modified and returned! *) }, {Integer, 2}],
        LFun["edgeListMarkWhenEitherPresent", {{Integer, 2, "Constant"}, {Integer, 1, "Constant"}}, {Integer, 1}],
        LFun["edgeListMarkWhenBothPresent", {{Integer, 2, "Constant"}, {Integer, 1, "Constant"}}, {Integer, 1}],
        LFun["edgeListReindexAfterDelete", {{Integer, 2} (* not Constant *), {Integer, 1} (* not Constant *)}, {Integer, 2}],
        LFun["edgeListReindex", {{Integer, 2} (* not Constant *), {Integer, 1, "Constant"}}, {Integer, 2}],

        LFun["symmetricTree", {{Integer, 1, "Constant"}}, {Integer, 2}],

        (* Matrix functions *)

        LFun["takeLowerInteger", {{Integer, 2, "Constant"}}, {Integer, 1}],
        LFun["takeLowerReal", {{Real, 2, "Constant"}}, {Real, 1}],
        LFun["takeLowerComplex", {{Complex, 2, "Constant"}}, {Complex, 1}],
        LFun["takeUpperInteger", {{Integer, 2, "Constant"}}, {Integer, 1}],
        LFun["takeUpperReal", {{Real, 2, "Constant"}}, {Real, 1}],
        LFun["takeUpperComplex", {{Complex, 2, "Constant"}}, {Complex, 1}],

        LFun["upperIndexPairPositions", {{Integer, 2, "Constant"}, Integer}, {Integer, 2}],
        LFun["lowerIndexPairPositions", {{Integer, 2, "Constant"}, Integer}, {Integer, 2}],

        (* used in igWeightedAdjacencyGraph *)
        LFun["upperPos", {{Integer, 2, "Constant"}, True|False (* diagonal *)}, {Integer, 1}],
        LFun["nondiagPos", {{Integer, 2, "Constant"}}, {Integer, 1}]
      }
    ],

    LClass["IG",
      {
        (* Create *)

        (* LFun["fromEdgeList", {{Real, _, "Constant"} (* edges *), Integer (* vertex count *), True|False (* directed *)}, "Void"], *)
        LFun["fromIncidenceMatrix", {{LType[SparseArray, Integer], "Constant"}, True|False (* directed *)}, "Void"],
        LFun["makeEdgeless", {Integer (* vertex count *)}, "Void"],
        (* LFun["fromEdgeListML", LinkObject], *)
        LFun["realizeDegreeSequence", {{Real, 1, "Constant"}, {Real, 1, "Constant"}, Integer}, "Void"],
        LFun["fromLCF", {Integer, {Real, 1, "Constant"}, Integer}, "Void"],
        LFun["makeLattice", {{Real, 1, "Constant"}, Integer (* nei *), True|False (* directed *), True|False (* mutual *), True|False (* periodic *)}, "Void"],
        LFun["kautz", {Integer, Integer}, "Void"],
        LFun["tree", {Integer, Integer, True|False (* directed *)}, "Void"],
        LFun["fromPrufer", {{Integer, 1, "Constant"}}, "Void"],
        LFun["completeGraph", {Integer, True|False (* directed *), True|False (* loops *)}, "Void"],
        LFun["completeCitationGraph", {Integer, True|False (* directed *)}, "Void"],
        LFun["deBruijn", {Integer, Integer}, "Void"],
        LFun["extendedChordalRing", {Integer, {Real, 2}, True|False (* directed *), True|False (* self loops *), True|False (* multiple edges *)}, "Void"],
        LFun["graphAtlas", {Integer}, "Void"],

        (* Directedness *)

        LFun["makeDirected", {}, "Void"],
        LFun["makeUndirected", {}, "Void"],

        (* Weights *)

        LFun["setWeights", {{Real, 1, "Constant"}}, "Void"],
        LFun["getWeights", {}, {Real, 1}],
        LFun["clearWeights", {}, "Void"],
        LFun["weightedQ", {}, True|False],

        (* Games *)

        LFun["treeGame", {Integer (* n *), True|False (* directed *), Integer (* method *)}, "Void"],
        LFun["degreeSequenceGame", {{Real, 1, "Constant"} (* outdeg *), {Real, 1, "Constant"} (* indeg *), Integer (* method *)}, "Void"],
        LFun["kRegularGame", {Integer, Integer, True|False (* directed *), True|False (* multiple *)}, "Void"],
        LFun["stochasticBlockModel", {{Real, 2, "Constant"}, {Integer, 1, "Constant"}, True|False (* directed *), True|False (* loops *)}, "Void"],
        LFun["forestFireGame", {Integer (* vertex count *), Real (* fwprob *), Real (* bwratio *), Integer (* nambs *), True|False (* directed *)}, "Void"],

        LFun["bipartiteGameGNM", {Integer (* n1 *), Integer (* n2 *), Integer (* m *), True|False (* directed *), True|False (* bidirectional *)}, "Void"],
        LFun["bipartiteGameGNP", {Integer (* n1 *), Integer (* n2 *), Real (* p *), True|False (* directed *), True|False (* bidirectional *)}, "Void"],

        LFun["erdosRenyiGNM", {Integer (* n *), Integer (* m *), True|False (* directed *), True|False (* loops *)}, "Void"],
        LFun["erdosRenyiGNP", {Integer (* n *), Real (* p *), True|False (* directed *), True|False (* loops *)}, "Void"],

        LFun["geometricGame", {Integer (* n *), Real (* radius *), True|False (* periodic *)}, {Real, 2} (* coordinates *)],

        LFun["barabasiAlbertGame", {Integer (* n *), Real (* power *), Real (* A *), Integer (* m *), {Real, 1, "Constant"} (* mvec *), True|False (* directed *), True|False (* totalDegree *), Integer (* method *)}, "Void"],
        LFun["barabasiAlbertGameWithStartingGraph", {Integer (* n *), Real (* power *), Real (* A *), Integer (* m *), {Real, 1, "Constant"} (* mvec *), True|False (* directed *), True|False (* totalDegree *), Integer (* method *), LExpressionID["IG"]}, "Void"],

        LFun["wattsStrogatzGame", {Integer (* dim *), Integer (* size *), Integer (* radius *), Real (* p *), True|False (* loops *), True|False (* multiple *)}, "Void"],

        LFun["staticFitnessGame", {Integer (* edges *), {Real, 1, "Constant"} (* out-fitness *), {Real, 1, "Constant"} (* in-fitness *), True|False (* loops *), True|False (* multiple *)}, "Void"],

        LFun["staticPowerLawGame", {Integer (* n *), Integer (* m *), Real (* expOut *), Real (* expIn *), True|False (* loops *), True|False (* multiple *), True|False (* finiteSizeCorrection *)}, "Void"],

        LFun["growingGame", {Integer (* n *), Integer (* m *), True|False (* directed *), True|False (* citation *)}, "Void"],

        LFun["callawayTraitsGame", {Integer (* n *), Integer (* types *), Integer (* k *), {Real, 1, "Constant"} (* type distribution *), {Real, 2, "Constant"} (* pref matrix *), True|False (* directed *)}, "Void"],

        LFun["establishmentGame", {Integer (* n *), Integer (* types *), Integer (* k *), {Real, 1, "Constant"} (* type distribution *), {Real, 2} (* pref matrix *), True|False (* directed *)}, "Void"],

        (* Modification *)

        LFun["connectNeighborhood", {Integer (* order *)}, "Void"],
        LFun["mycielski", {}, "Void"],

        (* Structure *)

        LFun["edgeCount", {}, Integer],
        LFun["vertexCount", {}, Integer],

        LFun["edgeList", {}, {Real, 2}],

        (* Testing *)

        LFun["directedQ", {}, True|False],
        LFun["dagQ", {}, True|False],
        LFun["simpleQ", {}, True|False],
        LFun["connectedQ", {True|False (* strongly connected *)}, True|False],
        LFun["treeQ", {Integer (* mode *)}, True|False],
        LFun["forestQ", {Integer (* mode *)}, True|False],
        LFun["bipartiteQ", {}, True|False],
        LFun["cactusQ", {}, True|False],
        LFun["nonSimpleCactusQ", {}, True|False],

        (* Centrality *)

        LFun["betweenness", {True|False (* nobigint *), True|False (* normalized *), {Real, 1, "Constant"} (* vertices *)}, {Real, 1}],
        LFun["edgeBetweenness", {True|False (* normalized *)}, {Real, 1}],
        LFun["closeness", {True|False (* normalized *), {Real, 1, "Constant"} (* vertices *)}, {Real, 1}],

        LFun["betweennessEstimate", {Real (* cutoff *), True|False (* nobigint *), True|False (* normalized *), {Real, 1, "Constant"} (* vertices *)}, {Real, 1}],
        LFun["edgeBetweennessEstimate", {Real (* cutoff *), True|False (* normalized *)}, {Real, 1}],
        LFun["closenessEstimate", {Real (* cutoff *), True|False (* normalized *), {Real, 1, "Constant"} (* vertices *)}, {Real, 1}],

        LFun["pageRank", {Integer (* method *), Real (* damping *), True|False (* directed *), Integer (* powerNiter *), Real (* powerEpsilon *)}, {Real, 1}],
        LFun["personalizedPageRank", {Integer (* method *), {Real, 1, "Constant"}, Real (* damping *), True|False (* directed *), Integer (* powerNiter *), Real (* powerEpsilon *)}, {Real, 1}],
        LFun["eigenvectorCentrality", {True|False (* directed *), True|False (* normalized *)}, {Real, 1}],
        LFun["hubScore", {True|False (* normalized *)}, {Real, 1}],
        LFun["authorityScore", {True|False (* normalized *)}, {Real, 1}],
        LFun["constraintScore", {}, {Real, 1}],

        (* Centralization *)

        LFun["degreeCentralization", {Integer, True|False, True|False}, Real],
        LFun["betweennessCentralization", {True|False, True|False}, Real],
        LFun["closenessCentralization", {True|False}, Real],
        LFun["eigenvectorCentralization", {True|False, True|False}, Real],
        LFun["centralization", {{Real, 1, "Constant"}, Real, True|False}, Real],

        (* Randomize *)

        LFun["rewire", {Integer (* n_trials *), True|False (* loops *)}, "Void"],
        LFun["rewireEdges", {Real (* probability *), True|False (* loops *), True|False (* multiple *)}, "Void"],
        LFun["rewireDirectedEdges", {Real (* probability *), True|False (* loops *), True|False (* outEdges *)}, "Void"],

        (* Isomorphism *)

        LFun["isomorphic", {LExpressionID["IG"]}, True|False],
        LFun["subisomorphic", {LExpressionID["IG"]}, True|False],
        LFun["getIsomorphism", {LExpressionID["IG"]}, {Real, 1}],
        LFun["getSubisomorphism", {LExpressionID["IG"]}, {Real, 1}],
        LFun["isoclass", {}, Integer],

        LFun["blissCanonicalPermutation", {Integer (* splitting heuristics *), {Integer, 1, "Constant"} (* colour *)}, {Real, 1}],
        LFun["blissIsomorphic", {LExpressionID["IG"], Integer (* splitting heuristics *), {Integer, 1, "Constant"} (* color1 *), {Integer, 1, "Constant"} (* color 2 *)}, True|False],
        LFun["blissFindIsomorphism", {LExpressionID["IG"], Integer (* splitting heuristics *), {Integer, 1, "Constant"} (* color1 *), {Integer, 1, "Constant"} (* color 2 *)}, {Real, 1}],
        LFun["blissAutomorphismCount", LinkObject],
        LFun["blissAutomorphismGroup", LinkObject],

        LFun["vf2Isomorphic", {LExpressionID["IG"], {Integer, 1, "Constant"}, {Integer, 1, "Constant"}, {Integer, 1, "Constant"}, {Integer, 1, "Constant"}}, True|False],
        LFun["vf2FindIsomorphisms", LinkObject],
        LFun["vf2Subisomorphic", {LExpressionID["IG"], {Integer, 1, "Constant"}, {Integer, 1, "Constant"}, {Integer, 1, "Constant"}, {Integer, 1, "Constant"}}, True|False],
        LFun["vf2FindSubisomorphisms", LinkObject],
        LFun["vf2IsomorphismCount", {LExpressionID["IG"], {Integer, 1, "Constant"}, {Integer, 1, "Constant"}, {Integer, 1, "Constant"}, {Integer, 1, "Constant"}}, Integer],
        LFun["vf2SubisomorphismCount", {LExpressionID["IG"], {Integer, 1, "Constant"}, {Integer, 1, "Constant"}, {Integer, 1, "Constant"}, {Integer, 1, "Constant"}}, Integer],
        LFun["vf2IsomorphicMulti", {LExpressionID["IG"]}, True|False],
        LFun["vf2SubisomorphicMulti", {LExpressionID["IG"]}, True|False],

        LFun["ladSubisomorphic", {LExpressionID["IG"], True|False (* induced *)}, True|False],
        LFun["ladSubisomorphicColored", LinkObject],
        LFun["ladGetSubisomorphism", {LExpressionID["IG"], True|False (* induced *)}, {Real, 1}],
        LFun["ladGetSubisomorphismColored", LinkObject],
        LFun["ladFindSubisomorphisms", LinkObject],
        LFun["ladCountSubisomorphisms", {LExpressionID["IG"], True|False (* induced *)}, Integer],
        LFun["ladCountSubisomorphismsColored", LinkObject],

        LFun["coloredSimpleGraph", LinkObject],

        (* Functions related to isomorphism *)

        LFun["selfComplementaryQ", {}, True|False],

        (* Topological sorting and directed acyclic graphs *)

        LFun["topologicalSorting", {}, {Real, 1}],
        LFun["feedbackArcSet", {True|False}, {Real, 1}],

        (* Motifs and subgraph counts *)

        LFun["dyadCensus", {}, {Integer, 1}],
        LFun["triadCensus", {}, {Real, 1}],
        LFun["motifs", {Integer (* size *), {Real, 1, "Constant"} (* cut_prob *)}, {Real, 1}],
        LFun["motifsNo", {Integer (* size *), {Real, 1, "Constant"} (* cut_prob *)}, Integer],
        LFun["motifsEstimate", {Integer (* size *), {Real, 1, "Constant"} (* cut_prob *), Integer (* sample_size *)}, Integer],

        LFun["triangles", {}, {Integer, 1}],
        LFun["countAdjacentTriangles", {{Real, 1, "Constant"}}, {Real, 1}],

        (* Shortest paths *)

        LFun["shortestPaths", {{Real, 1, "Constant"} (* from *), {Real, 1, "Constant"} (* to *)}, {Real, 2}],
        LFun["shortestPathCounts", {}, {Real, 1}],
        LFun["shortestPathCounts2", {{Integer, 1, "Constant"}}, {Real, 1}],
        LFun["neighborhoodSize", {{Real, 1, "Constant"}, Integer, Integer}, {Real, 1}],
        LFun["shortestPathWeightedHistogram", {Real (* bin size *), {Real, 1, "Constant"} (* from *), {Real, 1, "Constant"} (* to *), Integer (* method *)}, {Integer, 1}],
        LFun["averagePathLength", {}, Real], (* TODO; not currently in use; averagePathLengthWeighted() will call averagePathLength() in C code when needed *)
        LFun["averagePathLengthWeighted", {Integer}, Real],
        LFun["girth", {}, Real],
        LFun["radius", {}, Real],
        LFun["eccentricity", {{Real, 1, "Constant"}}, {Real, 1}],

        LFun["shortestPathsDijkstra", {{Real, 1, "Constant"} (* from *), {Real, 1, "Constant"} (* to *)}, {Real, 2}],
        LFun["shortestPathsBellmanFord", {{Real, 1, "Constant"} (* from *), {Real, 1, "Constant"} (* to *)}, {Real, 2}],
        LFun["shortestPathsJohnson", {{Real, 1, "Constant"} (* from *), {Real, 1, "Constant"} (* to *)}, {Real, 2}],

        LFun["diameter", {True|False (* by components *)}, Integer],
        LFun["findDiameter", {True|False (* by components *)}, {Real, 1}],
        LFun["diameterDijkstra", {True|False (* by components *)}, Real],
        LFun["findDiameterDijkstra", {True|False (* by components *)}, {Real, 1}],

        (* Cliques *)

        LFun["cliques", {Integer, Integer}, {Integer, 1}],
        LFun["cliqueDistribution", {Integer, Integer}, {Real, 1}],
        LFun["maximalCliques", {Integer, Integer}, {Integer, 1}],
        LFun["largestCliques", {}, {Integer, 1}],
        LFun["maximalCliquesCount", {Integer, Integer}, Integer],
        LFun["maximalCliqueDistribution", {Integer, Integer}, {Real, 1}],
        LFun["cliqueNumber", {}, Integer],
        LFun["cliquesWeighted", {Integer (* min_weight *), Integer (* max_weight *), {Real, 1, "Constant"} (* vertex_weights *), True|False (* maximal *)}, {Integer, 1}],
        LFun["largestCliquesWeighted", {{Real, 1, "Constant"} (* vertex_weights *)}, {Integer, 1}],
        LFun["cliqueNumberWeighted", {{Real, 1, "Constant"} (* vertex_weights *)}, Integer],

        (* Independent vertex sets *)

        LFun["independentVertexSets", {Integer, Integer}, {Integer, 1}],
        LFun["largestIndependentVertexSets", {}, {Integer, 1}],
        LFun["maximalIndependentVertexSets", {}, {Integer, 1}],
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

        LFun["layoutGEM",
          {{Real, 2, "Constant"} (* initial position *), True|False (* use initial *),
            Integer (* maxiter *), Real (* temp_min *), Real (* temp_max *), Real (* temp_init *)},
          {Real, 2}
        ],

        LFun["layoutDavidsonHarel",
          {{Real, 2, "Constant"} (* initial position *), True|False (* use initial *),
            Integer (* maxiter *), Integer (* fineiter *), Real (* cool_fact *),
            Real (* weight_node_dist *), Real (* weight_border *),
            Real (* weight_edge_lengths *), Real (* weight_edge_crossings *),
            Real (* weight_node_edge_dist *)},
          {Real, 2}
        ],

        LFun["layoutMDS", {{Real, 2, "Constant"}, Integer}, {Real, 2}],

        LFun["layoutReingoldTilford", {{Real, 1, "Constant"} (* roots *), True|False (* directed *)}, {Real, 2}],
        LFun["layoutReingoldTilfordCircular", {{Real, 1, "Constant"} (* roots *), True|False (* directed *)}, {Real, 2}],

        LFun["layoutDrL", {{Real, 2, "Constant"} (* initial positions *), True|False (* use initial *), Integer (* settings template *)}, {Real, 2}],
        LFun["layoutDrL3D", {{Real, 2, "Constant"} (* initial positions *), True|False (* use initial *), Integer (* settings template *)}, {Real, 2}],

        LFun["layoutBipartite", {{Integer, 1, "Constant"} (* types *), Real (* hgap *), Real (* vgap *), Integer (* maxiter *)}, {Real, 2}],

        (* Clustering coefficient *)

        LFun["transitivityUndirected", {True|False (* exclude isolates? *)}, Real],
        LFun["transitivityLocalUndirected", {True|False (* exclude isolates? *)}, {Real, 1}],
        LFun["transitivityAverageLocalUndirected", {True|False (* exclude isolates? *)}, Real],
        LFun["transitivityBarrat", {True|False (* exclude isolates? *)}, {Real, 1}],

        (* Similarity *)

        LFun["similarityCocitation", {{Real, 1, "Constant"}}, {Real, 2}],
        LFun["similarityBibcoupling", {{Real, 1, "Constant"}}, {Real, 2}],
        LFun["similarityJaccard", {{Real, 1, "Constant"}, True|False (* self loops *)}, {Real, 2}],
        LFun["similarityDice", {{Real, 1, "Constant"}, True|False (* self loops *)}, {Real, 2}],
        LFun["similarityInverseLogWeighted", {{Real, 1, "Constant"}}, {Real, 2}],

        (* Chordal graphs *)

        LFun["maximumCardinalitySearch", {}, {Real, 1}],
        LFun["chordalQ", {}, True|False],
        LFun["chordalCompletion", {}, {Real, 1}],

        (* Vertex separators *)

        LFun["minimumSizeSeparators", {}, {Integer, 1}],
        LFun["minimalSeparators", {}, {Integer, 1}],
        LFun["separatorQ", {{Real, 1, "Constant"}}, True|False],
        LFun["minSeparatorQ", {{Real, 1, "Constant"}}, True|False],

        LFun["vertexConnectivity", {}, Integer],
        LFun["edgeConnectivity", {}, Integer],
        LFun["vertexConnectivityST", {Integer, Integer}, Integer],
        LFun["edgeConnectivityST", {Integer, Integer}, Integer],
        LFun["cohesiveBlocks", LinkObject],

        (* Connected components *)

        LFun["articulationPoints", {}, {Real, 1}],
        LFun["biconnectedComponents", {}, {Integer, 1}],
        LFun["biconnectedEdgeComponents", {}, {Integer, 1}],
        LFun["biconnectedQ", {}, True|False],
        LFun["bridges", {}, {Real, 1}],
        LFun["connectedComponentSizes", {True|False (* strongly connected? *)}, {Real, 1}],

        (* Graphlets *)

        LFun["graphlets", LinkObject],
        LFun["graphletBasis", LinkObject],
        LFun["graphletProject", LinkObject],

        (* Community detection *)

        LFun["modularity", {{Real, 1, "Constant"}}, Real],
        LFun["communityEdgeBetweenness", LinkObject],
        LFun["communityWalktrap", LinkObject],
        LFun["communityFastGreedy", LinkObject],
        LFun["communityOptimalModularity", LinkObject],
        LFun["communityMultilevel", LinkObject],
        LFun["communityLabelPropagation", LinkObject],
        LFun["communityInfoMAP", LinkObject],
        LFun["communitySpinGlass", LinkObject],
        LFun["communityLeadingEigenvector", LinkObject],
        LFun["communityFluid", LinkObject],

        (* Maximum flow *)

        LFun["gomoryHuTree", {LExpressionID["IG"]}, {Real, 1}],

        LFun["dominatorTree", {LExpressionID["IG"], Integer (* root *)}, "Void"],
        LFun["immediateDominators", {Integer (* root *)}, {Real, 1}],

        (* LFun["maxFlowMatrix", {Integer (* s *), Integer (* t *)}, LType[SparseArray, Real, 2]], *)
        LFun["maxFlow", {Integer (* s *), Integer (* t *), {Real, 1, "Constant"} (* capacities *)}, {Real, 1}],
        LFun["maxFlowValue", {Integer (* s *), Integer (* t *), {Real, 1, "Constant"} (* capacities *)}, Real],

        LFun["minCutValue", {}, Real],
        LFun["minCutValueST", {Integer (* s *), Integer (* t *)}, Real],
        LFun["minCut", {}, {Real, 1}],
        LFun["minCutST", {Integer (* s *), Integer (* t *)}, {Real, 1}],
        LFun["allCutsST", {Integer (* s *), Integer (* t *)}, {Integer, 1}],
        LFun["allMinCutsST", {Integer (* s *), Integer (* t *)}, {Integer, 1}],

        (* Unfold tree *)

        LFun["unfoldTree", {LExpressionID["IG"], {Real, 1, "Constant"} (* roots *), True|False (* directed *)}, {Real, 1}],

        (* Bipartite graphs *)

        LFun["bipartitePartitions", {}, {Integer, 1}],
        LFun["bipartiteProjection", {{Integer, 1}, LExpressionID["IG"], LExpressionID["IG"]}, {Integer, 1}],

        (* Vertex contraction *)
        LFun["contractVertices", {{Real, 1, "Constant"}}, "Void"],

        (* Random walk *)
        LFun["randomWalk", {Integer, Integer}, {Real, 1}],
        LFun["randomEdgeWalk", {Integer, Integer}, {Real, 1}],

        (* Spanning tree *)
        LFun["spanningTree", {}, {Real, 1}],
        LFun["randomSpanningTree", {Integer (* 0-based vertex id *)}, {Real, 1}],

        (* Coreness *)
        LFun["coreness", {Integer (* mode *)}, {Real, 1}],

        LFun["vertexColoring", {}, {Integer, 1}],

        (* Other functions *)

        LFun["treelikeComponents", {}, {Integer, 1}],

        LFun["smoothen", {LExpressionID["IG"]}, {Integer, 1}],

        LFun["coordinatesToEmbedding", {{Real, 2, "Constant"}}, {Integer, 1}],

        LFun["perfectQ", {}, True|False],

        LFun["strahlerNumber", {}, {Integer, 1}],

        LFun["toPrufer", {}, {Integer, 1}],

        LFun["intersectionArray", LinkObject]
      }
    ],

    LClass["IGEmbedding",
      {
        LFun["set", LinkObject],
        LFun["get", LinkObject],
        LFun["validQ", {}, True|False],
        LFun["faces", {}, {Integer, 1}],
        LFun["planarQ", {Integer (* number of connected components *)}, True|False],
        LFun["dualGraph", {}, {Integer, 1}]
      }
    ],

    LClass["IGLemonGraph",
      {
        LFun["fromEdgeList", {{Integer, 2, "Constant"}, Integer}, "Void"],
        LFun["planarQ", {}, True|False],
        LFun["kuratowskiSubgraph", {}, {Integer, 1}],
        LFun["layoutPlanar", {}, {Integer, 2}],
        LFun["planarEmbedding", {}, {Integer, 1}],
        LFun["embeddingToCoordinates", {LExpressionID["IGEmbedding"]}, {Integer, 2}],

        LFun["maximumMatching", {}, {Integer, 1}],
        LFun["matchingNumber", {}, Integer]
      }
    ],

    LClass["IGFlann2D",
      {
        LFun["setPoints", {{Real, 2, "Constant"}}, "Void"],
        LFun["query", {{Real, 1, "Constant"} (* centre *), Real (* distance *)}, {Integer, 1}],
        LFun["queryMultiple", {{Real, 2, "Constant"} (* centres *), {Real, 1, "Constant"} (* distances *)}, {Integer, 1}],
        LFun["neighborCounts", {{Real, 2, "Constant"} (* centres *), {Real, 1, "Constant"}}, {Integer, 1}],
        LFun["intersectionCounts",
          {{Real, 2, "Constant"} (* centres1 *),
           {Real, 2, "Constant"} (* centres2 *),
           {Real, 1, "Constant"} (* distances *)},
          {Integer, 1}],
        LFun["unionCounts",
          {{Real, 2, "Constant"} (* centres1 *),
           {Real, 2, "Constant"} (* centres2 *),
           {Real, 1, "Constant"} (* distances *)},
          {Integer, 1}]
      }
    ]
  }
];


(***** Load build settings, if present *****)

(* $buildSettings is normally defined in the build settings files.
   When no build settings file is loaded, it must be set to None. *)
$buildSettings = None;
If[FileExistsQ[$buildSettingsFile], Get[$buildSettingsFile] ]


(***** GetInfo for troubleshooting *****)

(* GetInfo[] must be defined before the LTemplate functions are loaded because package loading is aborted
   on load failure. We want to be able to get debugging information even in this situation. *)
IGraphM`Developer`GetInfo[] :=
    Module[{res = "", igver, osver},
      res = StringJoin[res, "Mathematica version: \n", System`$Version, "; Release number: ", ToString[$ReleaseNumber], "\n\n"];
      res = StringJoin[res, "Package version: \n", $packageVersion, "\n\n"];
      res = StringJoin[res, "Package location: \n", ToString@FindFile["IGraphM`"], "\n\n"];
      res = StringJoin[res, "Library location: \n", ToString@FindLibrary["IGraphM"], "\n\n"];
      igver = Quiet@IGVersion[];
      res = StringJoin[res, "IGVersion[]: \n", If[StringQ[igver], igver, "Failed."], "\n\n"];
      res = StringJoin[res, "Build settings: \n", ToString[$buildSettings], "\n\n"];
      osver = Quiet@StringTrim@Switch[$OperatingSystem,
        "MacOSX", Import["!sw_vers", "String"],
        "Unix", Import["!uname -a", "String"] <> Import["!lsb_release -a 2>/dev/null", "String"],
        "Windows", Import["!cmd /C ver", "String"]
      ];
      res = StringJoin[res, "Operating system: \n", If[StringQ[osver], osver, "Failed."]]; (* no newline after last item *)
      res
    ]


(***** Compilation, loading and initialization *****)

(* Add $libraryDirectory to $LibraryPath in case package is not installed in Applications. *)
If[Not@MemberQ[$LibraryPath, $libraryDirectory],
  PrependTo[$LibraryPath, $libraryDirectory]
]


IGraphM`Developer`Recompile::build = "No build settings found. Please check BuildSettings.m."

IGraphM`Developer`Recompile[] :=
    Module[{},
      If[$buildSettings === None,
        Message[IGraphM`Developer`Recompile::build];
        Return[$Failed]
      ];
      If[Not@DirectoryQ[$libraryDirectory],
        CreateDirectory[$libraryDirectory]
      ];
      SetDirectory[$sourceDirectory];
      CompileTemplate[template, {"IGlobal.cpp", "IG.cpp", "IGRNG.cpp"},
        "ShellCommandFunction" -> Print, "ShellOutputFunction" -> Print,
        "TargetDirectory" -> $libraryDirectory,
        Sequence @@ $buildSettings
      ];
      ResetDirectory[];
      LoadIGraphM[]
    ]


PackageScope["igraphGlobal"]
igraphGlobal::usage =
    "igraphGlobal is the unique IGlobal object. There should only be a single object of this type; it's set in LoadIGraphM[] below.";

LoadIGraphM[] :=
    Module[{deps},
      deps = FileNameJoin[{$libraryDirectory, "dependencies.m"}];
      Check[
        If[FileExistsQ[deps], Get[deps]],
        Return[$Failed]
      ];
      If[LoadTemplate[template] === $Failed,
        Return[$Failed]
      ];
      igraphGlobal = Make["IGlobal"];
      igraphGlobal@"init"[];
    ]


(* Load library, compile if necessary. *)
If[LoadIGraphM[] === $Failed,
  Print[Style["Loading failed, trying to recompile ...", Red]];
  If[IGraphM`Developer`Recompile[] === $Failed
    ,
    Print[Style["Cannot load or compile library. \[FreakedSmiley] Aborting.", Red]];
    Abort[] (* verified that it is safe to abort while loading a new-style package *)
    ,
    Print[Style["Successfully compiled and loaded the library. \[HappySmiley]", Red]];
  ]
  ,
  If[Not@MemberQ[Stack[], BeginPackage],
    Print["IGraph/M " <> $packageVersion];
    Print@If[$Notebooks, (* The string with hyperlink would print fine in normal command line mode, but not when running in a script *)
      "Evaluate \!\(\*ButtonBox[\"IGDocumentation[]\",BaseStyle->\"Link\",ButtonData->\"paclet:IGraphM\"]\) to get started.",
      "Evaluate IGDocumentation[] to get started."
    ]
  ]
]


(***** General messages *****)

(* General::invopt is not present before Mathematica version 10.3. We set it up manually when needed. *)
If[Not@ValueQ[General::invopt],
  General::invopt = "Invalid value `1` for parameter `2`. Using default value `3`.";
]

IGraphM::mixed = "Mixed graphs are not supported by IGraph/M. Use DirectedGraph to convert undirected edges to two reciprocal directed ones.";


(***** Helper functions *****)


(* Get an IG compatible edge list. *)
(* This implementation attempts to select the fastest method based on the internal representation
   of the graph. With the "Simple" representation, IndexGraph is very fast. With "Incidence" it's
   slower than the Lookup method. With "NullGraph", performance doesn't matter.

   While GraphComputation`GraphRepresentation is an internal undocumented function, hopefully this
   is robust against changes as both branches of the If are valid ways to retrieve
   the edge list for any graph. They only differ in performance.
*)
(*
igEdgeList[graph_] :=
    Developer`ToPackedArray@If[GraphComputation`GraphRepresentation[graph] === "Simple",
      Flatten[EdgeList@IndexGraph[graph, 0], 1, If[DirectedGraphQ[graph], DirectedEdge, UndirectedEdge]]
      ,
      Lookup[
        AssociationThread[VertexList[graph], Range@VertexCount[graph] - 1],
        Flatten[EdgeList[graph], 1, If[DirectedGraphQ[graph], DirectedEdge, UndirectedEdge]]
      ]
    ]
*)
(* igEdgeList[graph_] := List @@@ EdgeList@IndexGraph[graph, 0]; *)

(* Not currently in use; was originally used in igMake and related functions.
 * See igraphGlobal@"incidenceToEdgeList" for a faster solution if need arises in the future. *)
(* Thanks to Carl Woll for the following implementation idea: http://community.wolfram.com/groups/-/m/t/1250373 *)
(*
igEdgeList[graph_?EmptyGraphQ] := {}
igEdgeList[graph_?MultigraphQ] :=
    Developer`ToPackedArray@Lookup[
      AssociationThread[VertexList[graph], Range@VertexCount[graph] - 1],
      Flatten[EdgeList[graph], 1, If[DirectedGraphQ[graph], DirectedEdge, UndirectedEdge]]
    ]
igEdgeList[graph_?UndirectedGraphQ] :=
    With[{sa = UpperTriangularize@WeightedAdjacencyMatrix[graph, EdgeWeight -> Range@EdgeCount[graph]]},
      sa["NonzeroPositions"][[Ordering @ sa["NonzeroValues"]]] - 1
    ]
igEdgeList[graph_?DirectedGraphQ] :=
    With[{sa = WeightedAdjacencyMatrix[graph, EdgeWeight -> Range@EdgeCount[graph]]},
      sa["NonzeroPositions"][[Ordering @ sa["NonzeroValues"]]] - 1
    ]
*)


(* Convert IG format vertex or edge index vector to Mathematica format. *)
PackageScope["igIndexVec"]
igIndexVec::usage = "igIndexVec[expr]";
igIndexVec[expr_LibraryFunctionError] := expr (* hack: allows LibraryFunctionError to fall through *)
igIndexVec[arr_] := 1 + Round[arr]

(* igEdgeWeightedQ: We only want edge-weighted graphs, not vertex weighted ones. *)
(* igEdgeWeightedQ = WeightedGraphQ[#] && PropertyValue[#, EdgeWeight] =!= Automatic &; *)

IGraphM::invw = "Invalid edge weight vector. Edge weights will be ignored.";

(* Create IG object from Mathematica Graph. Must be used when edge ordering matters. *)
(*
igMake[g_] :=
    With[{ig = Make["IG"]},
      ig@"fromEdgeList"[igEdgeList[g], VertexCount[g], igDirectedQ[g]];
      If[IGEdgeWeightedQ[g],
        Check[
          ig@"setWeights"[igEdgeWeights[g]],
          Message[IGraphM::invw]
        ]
      ];
      ig
    ]
*)

PackageScope["igMakeEmpty"]
igMakeEmpty::usage = "igMakeEmpty[] creates an empty IG object.";
igMakeEmpty[] := Make["IG"]

PackageScope["igMake"]
igMake::usage = "igMake[graph] creates an IG object from a Graph.";
igMake[g_] :=
    With[{ig = Make["IG"]},
      If[EmptyGraphQ[g],
        ig@"makeEdgeless"[VertexCount[g]]
        ,
        ig@"fromIncidenceMatrix"[IncidenceMatrix[g], DirectedGraphQ[g] (* empty graphs are treated as undirected in the branch above *)];
        If[IGEdgeWeightedQ[g],
          Check[
            ig@"setWeights"[igEdgeWeights[g]],
            Message[IGraphM::invw]
          ]
        ]
      ];
      ig
    ]

PackageScope["igMakeUnweighted"]
igMakeUnweighted::usage = "igMakeUnweighted[graph] is like igMake, but does not transfer weights. Use for performance when weights are not needed.";
igMakeUnweighted[g_] :=
    With[{ig = Make["IG"]},
      If[EmptyGraphQ[g],
        ig@"makeEdgeless"[VertexCount[g]]
        ,
        ig@"fromIncidenceMatrix"[IncidenceMatrix[g], DirectedGraphQ[g] (* empty graphs are treated as undirected in the branch above *)];
      ];
      ig
    ]

(*
(* Fast version. Use only for unweighted graphs and when edge ordering doesn't matter. *)
igMakeFast[g_?MultigraphQ] := igMake[g]
igMakeFast[g_?EmptyGraphQ] :=
    With[{ig = Make["IG"]},
      ig@"fromEdgeList"[{}, VertexCount[g], False];
      ig
    ]
igMakeFast[g_] :=
    With[{ig = Make["IG"]},
      If[DirectedGraphQ[g], (* empty graphs handled as undirected above *)
        ig@"fromEdgeList"[AdjacencyMatrix[g]["NonzeroPositions"] - 1, VertexCount[g], True],
        ig@"fromEdgeList"[UpperTriangularize[AdjacencyMatrix[g]]["NonzeroPositions"] - 1, VertexCount[g], False]
      ];
      ig
    ]
*)
PackageScope["igMakeFast"]
igMakeFast::usage = "igMakeFast[graph]";
igMakeFast = igMakeUnweighted; (* IncidenceMatrix-based igMake is faster than the above igMakeFast implementation *)

(*
(* Fast version. Use for graphs that may be weighted when edge ordering doesn't matter. *)
igMakeFastWeighted[g_?MultigraphQ] := igMake[g]
igMakeFastWeighted[g_?EmptyGraphQ] :=
    With[{ig = Make["IG"]},
      ig@"fromEdgeList"[{}, VertexCount[g], False];
      ig
    ]
igMakeFastWeighted[g_] :=
    With[{ig = Make["IG"]},
      If[IGEdgeWeightedQ[g],
        If[DirectedGraphQ[g], (* empty graphs handled as undirected above *)
          With[{wam = WeightedAdjacencyMatrix[g]},
            ig@"fromEdgeList"[wam["NonzeroPositions"] - 1, VertexCount[g], True];
            Check[ig@"setWeights"[wam["NonzeroValues"]], Message[IGraphM::invw]]
          ]
          ,
          With[{wam = UpperTriangularize@WeightedAdjacencyMatrix[g]},
            ig@"fromEdgeList"[wam["NonzeroPositions"] - 1, VertexCount[g], False];
            Check[ig@"setWeights"[wam["NonzeroValues"]], Message[IGraphM::invw]]
          ]
        ]
        ,
        If[DirectedGraphQ[g],
          ig@"fromEdgeList"[AdjacencyMatrix[g]["NonzeroPositions"] - 1, VertexCount[g], True],
          ig@"fromEdgeList"[UpperTriangularize[AdjacencyMatrix[g]]["NonzeroPositions"] - 1, VertexCount[g], False]
        ];
      ];
      ig
    ]
*)
PackageScope["igMakeFastWeighted"]
igMakeFastWeighted::usage = "igMakeFastWeighted[graph]";
igMakeFastWeighted = igMake; (* IncidenceMatrix-based igMake is faster than the above igMakeFast implementation *)


(* Create Mathematica Graph from IG object. *)
PackageScope["igToGraph"]
igToGraph::usage = "igToGraph[ig] converts back an IG object to a Graph.";
igToGraph[ig_, opt___] :=
    Graph[
      Range[ig@"vertexCount"[]],
      igIndexVec[ig@"edgeList"[]],
      DirectedEdges -> ig@"directedQ"[],
      opt
    ]

(* Create Mathematica Graph from IG object and assign vertex names. *)
(* This uses an undocumented syntax of Graph where the first argument is the vertex list,
   and the second argument is an edge list given as pairs of vertex indices. E.g.,
   Graph[{a,b,c}, {{1,2}, {2,3}}]
   TODO add test for this syntax
 *)
PackageScope["igToGraphWithNames"]
igToGraphWithNames::usage = "igToGraphWithNames[ig, vertexNames] converts back an IG object to a Graph and uses the give vertex names.";
igToGraphWithNames[ig_, verts_, opt___] :=
    Graph[
      verts,
      igIndexVec[ig@"edgeList"[]],
      DirectedEdges -> ig@"directedQ"[],
      opt
    ]


(* Warning: this function does not check if the graph is actually weighted! *)
PackageScope["igToWeightedGraphWithNames"]
igToWeightedGraphWithNames::usage = "igToWeightedGraphWithNames[graph]";
igToWeightedGraphWithNames[ig_, verts_] :=
    Graph[
      verts,
      igIndexVec[ig@"edgeList"[]],
      DirectedEdges -> ig@"directedQ"[],
      EdgeWeight -> ig@"getWeights"[]
    ]

(* Convert vertex indices to vertex names. *)
PackageScope["igVertexNames"]
igVertexNames::usage = "igVertexNames[graph][indices] converts vertex indices to vertex names.";
igVertexNames[graph_][indices_] := Part[VertexList[graph], indices]


(* Unpacks an index list representing vertex sets from an integer array,
   To be used in conjunction with IG::packListIntoIntTensor() *)
PackageScope["igUnpackVertexSet"]
igUnpackVertexSet::usage = "igUnpackVertexSet[graph][packed]";
igUnpackVertexSet[graph_?GraphQ][packed_] := igUnpackSetsHelper[VertexList[graph]][packed]
igUnpackVertexSet[emb_?AssociationQ][packed_] := igUnpackSetsHelper[Keys[emb]][packed]

PackageScope["igUnpackEdgeSet"]
igUnpackEdgeSet::usage = "igUnpackEdgeSet[graph][packed]";
igUnpackEdgeSet[graph_?GraphQ][packed_] := igUnpackSetsHelper[EdgeList[graph]][packed]
igUnpackSetsHelper[verts_][packed_] :=
    With[{len = First[packed]},
      partitionRagged[
        Part[verts, 1 + packed[[len + 2 ;; All]]],
        packed[[2 ;; len + 1]]
      ]
    ]

(* Convert vertex list to IG format *)
PackageScope["vss"]
vss::usage = "vss[graph][vertices]";
vss[graph_][All] := {}
vss[graph_][vs_List] := Check[VertexIndex[graph, #] - 1& /@ vs, throw[$Failed]]

PackageScope["vs"]
vs::usage = "vs[graph][vertex]";
vs[graph_][v_] := Check[VertexIndex[graph, v] - 1, throw[$Failed]]


(* Workarounds for Subgraph problems and cross-version changes. *)
PackageScope["igSubgraph"]
igSubgraph::usage = "igSubgraph[graph, spec]";
Which[
  (* In M12.0 and later Subgraph preserves properties, which makes it slow.
     We disable this when not needed, for performance. *)
  $VersionNumber >= 12.0,
  igSubgraph[args___] := Subgraph[args, Properties -> None] (* TODO verify in M12.0 final *)
  ,
  $VersionNumber >= 11.2,
  igSubgraph = Subgraph
  ,
  True, (* earlier than 11.2 *)
  (* Subgraph[Graph[{},{}], {}] does not evaluate in M11.1 and earlier *)
  igSubgraph[_, {}] := Graph[{},{}];
  igSubgraph[args___] := Subgraph[args]
];


(***** Public functions *****)

PackageExport["IGDocumentation"]
IGDocumentation::usage = "IGDocumentation[] opens the IGraph/M documentation.";

SyntaxInformation[IGDocumentation] = {"ArgumentsPattern" -> {}};
IGDocumentation[] :=
    If[$Notebooks,
      If[$VersionNumber == 11.1 && $OperatingSystem === "Unix", (* work around dysfunctional address bar in 11.1/Linux *)
        NotebookOpen@Documentation`ResolveLink["paclet:IGraphM"],
        Documentation`HelpLookupPacletURI["paclet:IGraphM"]
      ];
      ,
      Print["Built-in documentation is only available when running with a Front End.\nSee the online version at http://szhorvat.net/mathematica/IGraphM"]
    ]


PackageExport["IGVersion"]
IGVersion::usage = "IGVersion[] returns the IGraph/M version along with the version of the igraph library in use.";

SyntaxInformation[IGVersion] = {"ArgumentsPattern" -> {}};
IGVersion[] :=
    "IGraph/M " <> $packageVersion <>
    "\nigraph " <> igraphGlobal@"version"[] <> " (" <> igraphGlobal@"compilationDate"[] <> ")\n" <>
    $System;


PackageExport["IGSeedRandom"]
IGSeedRandom::usage =
    "IGSeedRandom[seed] seeds the current random number generator.\n" <>
    "IGSeedRandom[Method -> type] sets the current random number generator. Valid types are \"Mathematica\" and \"igraph\".";

Options[IGSeedRandom] = { Method -> Automatic };
SyntaxInformation[IGSeedRandom] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

igRandomMethods = <|"Mathematica" -> 0, "WolframLanguage" -> 0, "igraph" -> 1, "IGraph" -> 1|>;
IGSeedRandom::nogen =
    "`1` is not a valid random number generator type. Valid types are: " <>
    StringJoin@Riffle[Keys@DeleteDuplicates[igRandomMethods], ", "] <>
    ".";

IGSeedRandom[seed : (_?Internal`NonNegativeMachineIntegerQ | Automatic) : Automatic, opt : OptionsPattern[]] :=
    catch@Module[{},
      If[OptionValue[Method] =!= Automatic,
        check@igraphGlobal@"setRandomGenerator"[
          Lookup[igRandomMethods, OptionValue[Method],
            Message[IGSeedRandom::nogen, OptionValue[Method]]; throw[$Failed]]
        ]
      ];
      If[igraphGlobal@"randomGeneratorName"[] === "Mathematica",
        SeedRandom[seed]
        ,
        If[seed === Automatic,
          (* Notes on automatic seed generation for the igraph RNG:
              - We must get different seeds on subkernels, so we must not rely on the time only, and we must use
                the full available time resolution for the time component of the seed.
              - Reduce the dependence on the WL RNG, so that SeedRandom[fixed]; IGSeedRandom[] would not always
                have the same outcome.

             Hash[] converts AbsoluteTime[] to an integer while using the full time resolution (Round[] wouldn't).
             Hash[] seems to have a default range of 0..$MaxMachineInteger.
             BlockRandom[] is used to avoid propagating the WL random state when seeding the igraph RNG.
             igraph may only use values up to 2^31-1 on Windows (due to long int being 32-bit).
           *)
          check@igraphGlobal@"seedRandom"[ Mod[Hash@AbsoluteTime[] - BlockRandom@RandomInteger[Developer`$MaxMachineInteger], 2^31 - 1] ],
          check@igraphGlobal@"seedRandom"[seed]
        ]
      ]
    ]



(***** Backwards compatibility helpers *****)

IGraphM::depr = "`1` is deprecated. Use `2` instead.";

(* TODO: remove this eventually, along with all uses *)

PackageScope["multiEdgesOptionReplace"]
multiEdgesOptionReplace::usage =
    "multiEdgesOptionReplace@OptionValue[MultiEdges] provides backwards compatibility for the \"MultipleEdges\" option name.";
SetAttributes[multiEdgesOptionReplace, HoldAll]

multiEdgesOptionReplace[ov : OptionValue[syms_, opts_, MultiEdges]] :=
    If[KeyMemberQ[Flatten[opts], "MultipleEdges"],
      Message[IGraphM::depr, "The option name \"MultipleEdges\"", "MultiEdges"];
      OptionValue[
        syms,
        Replace[opts, (rule : Rule | RuleDelayed)["MultipleEdges", val_] :> rule[MultiEdges, val], {1}],
        MultiEdges
      ],
      ov
    ]
