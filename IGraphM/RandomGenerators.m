
(* Mathematica Package *)
(* Created by Mathematica plugin for IntelliJ IDEA *)

(* :Author: szhorvat *)
(* :Date: 2018-10-24 *)
(* :Copyright: (c) 2018-2020 Szabolcs Horvát *)

Package["IGraphM`"]
igContextSetup[igPackagePrivateSymbol]

(***********************************)
(***** Random graph generators *****)
(***********************************)


PackageExport["IGDegreeSequenceGame"]
IGDegreeSequenceGame::usage =
    "IGDegreeSequenceGame[degrees] generates an undirected random graph with the given degree sequence.\n" <>
    "IGDegreeSequenceGame[indegrees, outdegrees] generates a directed random graph with the given in- and out-degree sequences.";

Options[IGDegreeSequenceGame] = { Method -> "FastSimple" };
igDegreeSequenceGameMethods = <|"ConfigurationModel" -> 0, "ConfigurationModelSimple"->3, "FastSimple"->1, "VigerLatapy"->2|>;
amendUsage[IGDegreeSequenceGame, "Available Method options: <*Keys[igDegreeSequenceGameMethods]*>."];
igDegreeSequenceGameMethods = Join[
  igDegreeSequenceGameMethods,
  <|"Simple" -> 0, "SimpleNoMultiple" -> 1, "VigerLatapy" -> 2, "SimpleNoMultipleUniform" -> 3|> (* the old naming *)
];
SyntaxInformation[IGDegreeSequenceGame] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}, "OptionNames" -> optNames[IGDegreeSequenceGame, Graph]};

IGDegreeSequenceGame[degrees_?nonNegIntVecQ, opt : OptionsPattern[{IGDegreeSequenceGame, Graph}]] :=
    catch@applyGraphOpt[opt]@igDegreeSequenceGame[{}, degrees, OptionValue[Method]]

IGDegreeSequenceGame[indegrees_?nonNegIntVecQ, outdegrees_?nonNegIntVecQ, opt : OptionsPattern[{IGDegreeSequenceGame, Graph}]] :=
    catch@applyGraphOpt[opt]@igDegreeSequenceGame[indegrees, outdegrees, OptionValue[Method]]

igDegreeSequenceGame[indegrees_, outdegrees_, method_] :=
    (* no catch *) Block[{ig = igMakeEmpty[]},
      check@ig@"degreeSequenceGame"[outdegrees, indegrees, Lookup[igDegreeSequenceGameMethods, method, -1]];
      igToGraph[ig]
    ]


PackageExport["IGKRegularGame"]
IGKRegularGame::usage = "IGKRegularGame[n, k] generates a k-regular graph on n vertices, i.e. a graph in which all vertices have degree k.";

Options[IGKRegularGame] = { DirectedEdges -> False, MultiEdges -> False, "MultipleEdges" -> "Deprecated" };
SyntaxInformation[IGKRegularGame] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}, "OptionNames" -> optNames[IGKRegularGame, Graph]};
IGKRegularGame[n_?Internal`NonNegativeMachineIntegerQ, k_?Internal`NonNegativeMachineIntegerQ, opt : OptionsPattern[{IGKRegularGame, Graph}]] :=
    catch@Block[{ig = igMakeEmpty[]},
      check@ig@"kRegularGame"[n, k, OptionValue[DirectedEdges], multiEdgesOptionReplace@OptionValue[MultiEdges]];
      applyGraphOpt[opt]@igToGraph[ig]
    ]


PackageExport["IGStochasticBlockModelGame"]
IGStochasticBlockModelGame::usage = "IGStochasticBlockModelGame[ratesMatrix, blockSizes] samples from a stochastic block model.";

Options[IGStochasticBlockModelGame] = { SelfLoops -> False, DirectedEdges -> False };
SyntaxInformation[IGStochasticBlockModelGame] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}, "OptionNames" -> optNames[IGStochasticBlockModelGame, Graph]};
IGStochasticBlockModelGame[ratesMatrix_?SquareMatrixQ, blockSizes_?nonNegIntVecQ, opt : OptionsPattern[{IGStochasticBlockModelGame, Graph}]] :=
    catch@Block[{ig = igMakeEmpty[]},
      check@ig@"stochasticBlockModel"[Normal[ratesMatrix], Normal[blockSizes], OptionValue[DirectedEdges], OptionValue[SelfLoops]];
      applyGraphOpt[opt]@igToGraph[ig]
    ]


PackageExport["IGForestFireGame"]
IGForestFireGame::usage =
    "IGForestFireGame[n, pForward] generates a graph on n vertices from the forest fire model.\n" <>
    "IGForestFireGame[n, pForward, rBackward] specifies the backward to forward burning probability ratio (default: 1).\n" <>
    "IGForestFireGame[n, pForward, rBackward, nAmbassadors] also specifies the number of ambassador nodes in each step (default: 1).";

Options[IGForestFireGame] = { DirectedEdges -> True };
SyntaxInformation[IGForestFireGame] = {"ArgumentsPattern" -> {_, _, _., _., OptionsPattern[]}, "OptionNames" -> optNames[IGForestFireGame, Graph]};
IGForestFireGame[n_?Internal`PositiveMachineIntegerQ, fwprob_?NonNegative, bwratio : _?NonNegative : 1, nambs : _?Internal`NonNegativeIntegerQ : 1, opt : OptionsPattern[{IGForestFireGame, Graph}]] :=
    catch@Block[{ig = igMakeEmpty[]},
      check@ig@"forestFireGame"[n, fwprob, bwratio, nambs, OptionValue[DirectedEdges]];
      applyGraphOpt[opt]@igToGraph[ig]
    ]


PackageExport["IGBipartiteGameGNM"]
IGBipartiteGameGNM::usage = "IGBipartiteGameGNM[n1, n2, m] generates a bipartite random graph with n1 and n2 vertices in the two partitions and m edges.";

(* TODO: Work around Mathematica bug where bidirectional directed bipartite graphs may not be laid out. *)
Options[IGBipartiteGameGNM] = { DirectedEdges -> False, "Bidirectional" -> True, GraphLayout -> "BipartiteEmbedding" };
SyntaxInformation[IGBipartiteGameGNM] = {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}, "OptionNames" -> optNames[IGBipartiteGameGNM, Graph]};
IGBipartiteGameGNM[n1_?Internal`NonNegativeMachineIntegerQ, n2_?Internal`NonNegativeMachineIntegerQ, m_?Internal`NonNegativeMachineIntegerQ, opt : OptionsPattern[{IGBipartiteGameGNM, Graph}]] :=
    catch@Block[{ig = igMakeEmpty[]},
      check@ig@"bipartiteGameGNM"[n1, n2, m, OptionValue[DirectedEdges], OptionValue["Bidirectional"]];
      applyGraphOpt[GraphLayout -> OptionValue[GraphLayout], opt]@igToGraph[ig]
    ]


PackageExport["IGBipartiteGameGNP"]
IGBipartiteGameGNP::usage = "IGBipartiteGameGNP[n1, n2, p] generates a bipartite Bernoulli random graph with n1 and n2 vertices in the two partitions and connection probability p.";

Options[IGBipartiteGameGNP] = { DirectedEdges -> False, "Bidirectional" -> True, GraphLayout -> "BipartiteEmbedding" };
SyntaxInformation[IGBipartiteGameGNP] = {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}, "OptionNames" -> optNames[IGBipartiteGameGNP, Graph]};
IGBipartiteGameGNP[n1_?Internal`NonNegativeMachineIntegerQ, n2_?Internal`NonNegativeMachineIntegerQ, p_?NonNegative, opt : OptionsPattern[{IGBipartiteGameGNP, Graph}]] :=
    catch@Block[{ig = igMakeEmpty[]},
      check@ig@"bipartiteGameGNP"[n1, n2, p, OptionValue[DirectedEdges], OptionValue["Bidirectional"]];
      applyGraphOpt[GraphLayout -> OptionValue[GraphLayout], opt]@igToGraph[ig]
    ]


PackageExport["IGErdosRenyiGameGNM"]
IGErdosRenyiGameGNM::usage = "IGErdosRenyiGameGNM[n, m] generates a random graph with n vertices and m edges.";

Options[IGErdosRenyiGameGNM] = { DirectedEdges -> False, SelfLoops -> False };
SyntaxInformation[IGErdosRenyiGameGNM] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}, "OptionNames" -> optNames[IGErdosRenyiGameGNM, Graph]};
IGErdosRenyiGameGNM[n_?Internal`NonNegativeMachineIntegerQ, m_?Internal`NonNegativeMachineIntegerQ, opt : OptionsPattern[{IGErdosRenyiGameGNM, Graph}]] :=
    catch@Block[{ig = igMakeEmpty[]},
      check@ig@"erdosRenyiGNM"[n, m, OptionValue[DirectedEdges], OptionValue[SelfLoops]];
      applyGraphOpt[opt]@igToGraph[ig]
    ]


PackageExport["IGErdosRenyiGameGNP"]
IGErdosRenyiGameGNP::usage = "IGErdosRenyiGameGNP[n, p] generates a random graph on n vertices, in which each edge is present with probability p.";

Options[IGErdosRenyiGameGNP] = { DirectedEdges -> False, SelfLoops -> False };
SyntaxInformation[IGErdosRenyiGameGNP] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}, "OptionNames" -> optNames[IGErdosRenyiGameGNP, Graph]};
IGErdosRenyiGameGNP[n_?Internal`NonNegativeMachineIntegerQ, p_?NonNegative, opt : OptionsPattern[{IGErdosRenyiGameGNP, Graph}]] :=
    catch@Block[{ig = igMakeEmpty[]},
      check@ig@"erdosRenyiGNP"[n, p, OptionValue[DirectedEdges], OptionValue[SelfLoops]];
      applyGraphOpt[opt]@igToGraph[ig]
    ]


PackageExport["IGGeometricGame"]
IGGeometricGame::usage = "IGGeometricGame[n, radius] generates an n-vertex geometric random graph on the unit square by connecting points closer than radius.";

Options[IGGeometricGame] = {"Periodic" -> False};
SyntaxInformation[IGGeometricGame] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}, "OptionNames" -> optNames[IGGeometricGame, Graph]};
IGGeometricGame[n_?Internal`NonNegativeMachineIntegerQ, radius_?nonNegativeNumericQ, opt : OptionsPattern[{IGGeometricGame, Graph}]] :=
    catch@Block[{ig = igMakeEmpty[], coord},
      coord = check@ig@"geometricGame"[n, radius, OptionValue["Periodic"]];
      applyGraphOpt[VertexCoordinates -> coord, opt]@igToGraph[ig]
    ]


PackageExport["IGBarabasiAlbertGame"]
IGBarabasiAlbertGame::usage =
    "IGBarabasiAlbertGame[n, k] generates an n-vertex Barabási–Albert random graph by adding a new vertex with k (out-)edges in each step.\n" <>
    "IGBarabasiAlbertGame[n, {k2, k3, \[Ellipsis]}] generates an n-vertex Barabási–Albert random graph by adding a new vertex with k2, k3, \[Ellipsis] out-edges in each step.\n" <>
    "IGBarabasiAlbertGame[n, k, {\[Beta], A}] generates a Barabási–Albert random graph with preferential attachment probabilities proportional to d^\[Beta] + A where d is the vertex (in-)degree.";

Options[IGBarabasiAlbertGame] = {
  DirectedEdges -> True, "TotalDegreeAttraction" -> False,
  Method -> "PSumTree",
  "StartingGraph" -> None
};
SyntaxInformation[IGBarabasiAlbertGame] = {"ArgumentsPattern" -> {_, _, _., OptionsPattern[]}, "OptionNames" -> optNames[IGBarabasiAlbertGame, Graph]};
igBarabasiAlbertGameMethods = <|"Bag" -> 0, "PSumTree" -> 1, "PSumTreeMultiple" -> 2|>;
amendUsage[IGBarabasiAlbertGame, "Available Method options: <*Keys[igBarabasiAlbertGameMethods]*>."];

IGBarabasiAlbertGame::bdstart = "An invalid value was given for the \"StartingGraph\" option.";

IGBarabasiAlbertGame[
  n_?Internal`NonNegativeMachineIntegerQ, m : (_?Internal`PositiveMachineIntegerQ | _?nonNegIntVecQ),
  opt : OptionsPattern[{IGBarabasiAlbertGame, Graph}]] :=
    igBarabasiAlbertGame[n, m, {1,1}, OptionValue[DirectedEdges], OptionValue["TotalDegreeAttraction"], OptionValue[Method], OptionValue["StartingGraph"], opt]

IGBarabasiAlbertGame[
  n_?Internal`NonNegativeMachineIntegerQ, m : (_?Internal`PositiveMachineIntegerQ | _?nonNegIntVecQ),
  power_?NumericQ,
  opt : OptionsPattern[{IGBarabasiAlbertGame, Graph}]] :=
    igBarabasiAlbertGame[n, m, {power,1}, OptionValue[DirectedEdges], OptionValue["TotalDegreeAttraction"], OptionValue[Method], OptionValue["StartingGraph"], opt]

IGBarabasiAlbertGame[
  n_?Internal`NonNegativeMachineIntegerQ, m : (_?Internal`PositiveMachineIntegerQ | _?nonNegIntVecQ),
  {power_?NumericQ, a_?nonNegativeNumericQ},
  opt : OptionsPattern[{IGBarabasiAlbertGame, Graph}]] :=
    igBarabasiAlbertGame[n, m, {power, a}, OptionValue[DirectedEdges], OptionValue["TotalDegreeAttraction"], OptionValue[Method], OptionValue["StartingGraph"], opt]

igBarabasiAlbertGame[n_, m_, {power_, a_}, directed_, totalDegree_, method_, initial_, opt___] :=
    catch@Block[{ig = igMakeEmpty[], start},
      If[initial === None,
        check@ig@"barabasiAlbertGame"[
          n, power, a,
          If[ListQ[m], 0, m], If[ListQ[m], Prepend[m,0], {}],
          directed, totalDegree, Lookup[igBarabasiAlbertGameMethods, method, -1]
        ]
        ,
        If[Not@igGraphQ[initial],
          Message[IGBarabasiAlbertGame::bdstart];
          throw[$Failed]
        ];
        start = igMakeFast[initial];
        check@ig@"barabasiAlbertGameWithStartingGraph"[
          n, power, a,
          If[ListQ[m], 0, m], If[ListQ[m], m, {}],
          directed, totalDegree, Lookup[igBarabasiAlbertGameMethods, method, -1],
          ManagedLibraryExpressionID[start]
        ]
      ];
      applyGraphOpt[opt]@igToGraph[ig]
    ]


PackageExport["IGWattsStrogatzGame"]
IGWattsStrogatzGame::usage =
    "IGWattsStrogatzGame[n, p] generates an n-vertex Watts–Strogatz random graph using rewiring probability p.\n" <>
    "IGWattsStrogatzGame[n, p, k] rewires a lattice where each vertex is connected to its k-neighbourhood.\n" <>
    "IGWattsStrogatzGame[n, p, {dim, k}] rewires a dim dimensional lattice of n^dim vertices, where each vertex is connected to its k-neighbourhood.";

Options[IGWattsStrogatzGame] = {
  SelfLoops -> False, MultiEdges -> False, "MultipleEdges" -> "Deprecated"
};
SyntaxInformation[IGWattsStrogatzGame] = {
  "ArgumentsPattern" -> {_, _, _., OptionsPattern[]}, "OptionNames" -> optNames[IGWattsStrogatzGame, Graph]
};
IGWattsStrogatzGame[
  n_?Internal`NonNegativeMachineIntegerQ, p_?NonNegative,
  {dim_?Internal`PositiveMachineIntegerQ, k_?Internal`PositiveMachineIntegerQ},
  opt : OptionsPattern[{IGWattsStrogatzGame, Graph}]] :=
    catch@Block[{ig = igMakeEmpty[]},
      check@ig@"wattsStrogatzGame"[dim, n, k, p, OptionValue[SelfLoops], multiEdgesOptionReplace@OptionValue[MultiEdges]];
      applyGraphOpt[opt]@igToGraph[ig]
    ]
IGWattsStrogatzGame[
  n_?Internal`NonNegativeMachineIntegerQ, p_?NonNegative,
  k_?Internal`PositiveMachineIntegerQ, opt : OptionsPattern[]] :=
      IGWattsStrogatzGame[n, p, {1, k}, opt]
IGWattsStrogatzGame[
  n_?Internal`NonNegativeMachineIntegerQ, p_?NonNegative,
  opt : OptionsPattern[]] :=
    IGWattsStrogatzGame[n, p, {1, 2}, opt]


PackageExport["IGStaticFitnessGame"]
IGStaticFitnessGame::usage =
    "IGStaticFitnessGame[m, {f1, f2, \[Ellipsis]}] generates a random undirected graph with m edges where edge i <-> j is inserted with probability proportional to f_i\[Times]f_j.\n" <>
    "IGStaticFitnessGame[m, {fout1, fout2, \[Ellipsis]}, {fin1, fin2, \[Ellipsis]}] generates a random directed graph with m edges where edge i -> j is inserted with probability proportional to fout_i\[Times]fin_j.";

Options[IGStaticFitnessGame] = { SelfLoops -> False, MultiEdges -> False, "MultipleEdges" -> "Deprecated" };
SyntaxInformation[IGStaticFitnessGame] = {
  "ArgumentsPattern" -> {_, _, _., OptionsPattern[]}, "OptionNames" -> optNames[IGStaticFitnessGame, Graph]
};
IGStaticFitnessGame[
  m_?Internal`NonNegativeMachineIntegerQ,
  inFitness_?nonNegVecQ, outFitness : _?nonNegVecQ : {}, opt : OptionsPattern[{IGStaticFitnessGame, Graph}]] :=
    catch@Block[{ig = igMakeEmpty[]},
      check@ig@"staticFitnessGame"[
        m, Normal[inFitness], Normal[outFitness],
        OptionValue[SelfLoops], multiEdgesOptionReplace@OptionValue[MultiEdges]
      ];
      applyGraphOpt[opt]@igToGraph[ig]
    ]


PackageExport["IGStaticPowerLawGame"]
IGStaticPowerLawGame::usage =
    "IGStaticPowerLawGame[n, m, exp] generates a random graph with n vertices and m edges, having a power-law degree distribution with the given exponent.\n" <>
    "IGStaticPowerLawGame[n, m, expOut, expIn] generates a random directed graph with n vertices and m edges, having power-law in- and out-degree distributions with the given exponents.";

Options[IGStaticPowerLawGame] = {
  SelfLoops -> False,
  MultiEdges -> False, "MultipleEdges" -> "Deprecated",
  "FiniteSizeCorrection" -> True
};
SyntaxInformation[IGStaticPowerLawGame] = {
  "ArgumentsPattern" -> {_, _, _, _., OptionsPattern[]}, "OptionNames" -> optNames[IGStaticPowerLawGame, Graph]
};
IGStaticPowerLawGame[n_?Internal`NonNegativeMachineIntegerQ, m_?Internal`NonNegativeMachineIntegerQ, exp_?nonNegativeNumericQ, opt : OptionsPattern[{IGStaticPowerLawGame, Graph}]] :=
    catch@Block[{ig = igMakeEmpty[]},
      check@ig@"staticPowerLawGame"[n, m, exp, -1, OptionValue[SelfLoops], multiEdgesOptionReplace@OptionValue[MultiEdges], OptionValue["FiniteSizeCorrection"]];
      applyGraphOpt[opt]@igToGraph[ig]
    ]
IGStaticPowerLawGame[n_?Internal`NonNegativeMachineIntegerQ, m_?Internal`NonNegativeMachineIntegerQ, expOut_?nonNegativeNumericQ, expIn_?nonNegativeNumericQ, opt : OptionsPattern[{IGStaticPowerLawGame, Graph}]] :=
    catch@Block[{ig = igMakeEmpty[]},
      check@ig@"staticPowerLawGame"[n, m, expOut, expIn, OptionValue[SelfLoops], multiEdgesOptionReplace@OptionValue[MultiEdges], OptionValue["FiniteSizeCorrection"]];
      applyGraphOpt[opt]@igToGraph[ig]
    ]


PackageExport["IGGrowingGame"]
IGGrowingGame::usage = "IGGrowingGame[n, k] generates a growing random graph with n vertices, adding a new vertex and k new edges in each step.";

Options[IGGrowingGame] = { DirectedEdges -> False, "Citation" -> False };
SyntaxInformation[IGGrowingGame] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}, "OptionNames" -> optNames[IGGrowingGame, Graph]};
IGGrowingGame[n_?Internal`NonNegativeMachineIntegerQ, m_?Internal`NonNegativeMachineIntegerQ, opt : OptionsPattern[{IGGrowingGame, Graph}]] :=
    catch@Block[{ig = igMakeEmpty[]},
      check@ig@"growingGame"[n, m, OptionValue[DirectedEdges], OptionValue["Citation"]];
      applyGraphOpt[opt]@igToGraph[ig]
    ]


PackageExport["IGTreeGame"]
IGTreeGame::usage = "IGTreeGame[n] generates a random tree on n vertices. Sampling is uniform over the set of labelled trees.";

Options[IGTreeGame] = { Method -> "LoopErasedRandomWalk", DirectedEdges -> False };
igTreeGameMethods = <|"PruferCode" -> 0, "LoopErasedRandomWalk" -> 1|>;
amendUsage[IGTreeGame, "Available Method options: <*Keys[igTreeGameMethods]*>."];
SyntaxInformation[IGTreeGame] = {"ArgumentsPattern" -> {_, OptionsPattern[]}, "OptionNames" -> optNames[IGTreeGame, Graph]};
IGTreeGame[n_?Internal`NonNegativeIntegerQ, opt : OptionsPattern[{IGTreeGame, Graph}]] :=
    catch@Block[{ig = igMakeEmpty[]},
      check@ig@"treeGame"[n, OptionValue[DirectedEdges], Lookup[igTreeGameMethods, OptionValue[Method], -1]];
      applyGraphOpt[opt]@igToGraph[ig]
    ]


PackageExport["IGCallawayTraitsGame"]
IGCallawayTraitsGame::usage = "IGCallawayTraitsGame[n, k, typeWeights, preferenceMatrix]";

IGEstablishmentGame::prefmsym = IGCallawayTraitsGame::prefmsym = "The preference matrix must be symmetric when generating undirected graphs.";
IGEstablishmentGame::prefmdim = IGCallawayTraitsGame::prefmdim = "The preference matrix must be square and agree in size with the number of vertex types.";
IGEstablishmentGame::prefmel  = IGCallawayTraitsGame::prefmel  = "The elements of the preference matrix must be probabilities between 0 and 1.";
IGEstablishmentGame::weightnn = IGCallawayTraitsGame::weightnn = "The vertex type weights must be non-negative.";

Options[IGCallawayTraitsGame] = { DirectedEdges -> False };
SyntaxInformation[IGCallawayTraitsGame] = {
  "ArgumentsPattern" -> {_, _, _, _, OptionsPattern[]}, "OptionNames" -> optNames[IGCallawayTraitsGame, Graph]
};
IGCallawayTraitsGame[
  n_?Internal`NonNegativeMachineIntegerQ, k_?Internal`NonNegativeMachineIntegerQ,
  typeWeights_?VectorQ, prefMatrix_?MatrixQ,
  opt : OptionsPattern[{IGCallawayTraitsGame, Graph}]] :=
      catch@Block[{ig = igMakeEmpty[]},
        If[Not@nonNegVecQ[typeWeights],
          Message[IGCallawayTraitsGame::weightnn];
          Return[$Failed];
        ];
        If[Dimensions[prefMatrix] != {Length[typeWeights], Length[typeWeights]},
          Message[IGCallawayTraitsGame::prefmdim];
          Return[$Failed]
        ];
        If[Not@MatrixQ[prefMatrix, 0 <= # <= 1&],
          Message[IGCallawayTraitsGame::prefmel];
          Return[$Failed]
        ];
        If[Not@TrueQ@OptionValue[DirectedEdges] && Not@SymmetricMatrixQ[prefMatrix],
          Message[IGCallawayTraitsGame::prefmsym];
          Return[$Failed]
        ];
        check@ig@"callawayTraitsGame"[n, k, Normal[typeWeights, SparseArray], Normal[prefMatrix], OptionValue[DirectedEdges]];
        applyGraphOpt[opt]@igToGraph[ig]
      ]


PackageExport["IGEstablishmentGame"]
IGEstablishmentGame::usage = "IGEstablishmentGame[n, k, typeWeights, preferenceMatrix]";

Options[IGEstablishmentGame] = { DirectedEdges -> False };
SyntaxInformation[IGEstablishmentGame] = {
  "ArgumentsPattern" -> {_, _, _, _, OptionsPattern[]}, "OptionNames" -> optNames[IGEstablishmentGame, Graph]
};
IGEstablishmentGame[
  n_?Internal`NonNegativeMachineIntegerQ, k_?Internal`NonNegativeMachineIntegerQ,
  typeWeights_?VectorQ, prefMatrix_?MatrixQ,
  opt : OptionsPattern[{IGEstablishmentGame, Graph}]] :=
    catch@Block[{ig = igMakeEmpty[]},
      If[Not@positiveVecQ[typeWeights],
        Message[IGEstablishmentGame::weightnn];
        Return[$Failed];
      ];
      If[Dimensions[prefMatrix] != {Length[typeWeights], Length[typeWeights]},
        Message[IGEstablishmentGame::prefmdim];
        Return[$Failed]
      ];
      If[Not@MatrixQ[prefMatrix, 0 <= # <= 1&],
        Message[IGEstablishmentGame::prefmel];
        Return[$Failed]
      ];
      If[Not@TrueQ@OptionValue[DirectedEdges] && Not@SymmetricMatrixQ[prefMatrix],
        Message[IGEstablishmentGame::prefmsym];
        Return[$Failed]
      ];
      check@ig@"establishmentGame"[n, k, Normal[typeWeights, SparseArray], prefMatrix, OptionValue[DirectedEdges]];
      applyGraphOpt[opt]@igToGraph[ig]
    ]


PackageExport["IGPreferenceGame"]
IGPreferenceGame::usage = "IGPreferenceGame[n, typeWeights, preferenceMatrix]";

Options[IGPreferenceGame] = { DirectedEdges -> False, SelfLoops -> False };
SyntaxInformation[IGPreferenceGame] = {
  "ArgumentsPattern" -> {_, _, _, OptionsPattern[]}, "OptionNames" -> optNames[IGPreferenceGame, Graph]
};
IGPreferenceGame[
  n_?Internal`NonNegativeMachineIntegerQ,
  typeWeights_?VectorQ, prefMatrix_?MatrixQ,
  opt : OptionsPattern[{IGPreferenceGame, Graph}]] :=
      catch@Block[{ig = igMakeEmpty[]},
        check@ig@"preferenceGame"[n, Normal[typeWeights, SparseArray], prefMatrix, OptionValue[DirectedEdges], OptionValue[SelfLoops]];
        applyGraphOpt[opt]@igToGraph[ig]
      ]


PackageExport["IGAsymmetricPreferenceGame"]
IGAsymmetricPreferenceGame::usage = "IGAsymmetricPreferenceGame[n, typeWeightsMatrix, preferenceMatrix]";

Options[IGAsymmetricPreferenceGame] = { SelfLoops -> False };
SyntaxInformation[IGAsymmetricPreferenceGame] = {
  "ArgumentsPattern" -> {_, _, _, OptionsPattern[]}, "OptionNames" -> optNames[IGPreferenceGame, Graph]
};
IGAsymmetricPreferenceGame[
  n_?Internal`NonNegativeMachineIntegerQ,
  typeWeightsMatrix_?MatrixQ, prefMatrix_?MatrixQ,
  opt : OptionsPattern[{IGAsymmetricPreferenceGame, Graph}]] :=
      catch@Block[{ig = igMakeEmpty[]},
        check@ig@"asymmetricPreferenceGame"[n, Normal[typeWeightsMatrix, SparseArray], prefMatrix, OptionValue[SelfLoops]];
        applyGraphOpt[opt]@igToGraph[ig]
      ]
