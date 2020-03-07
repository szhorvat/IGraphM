(* Mathematica Package *)
(* Created by Mathematica plugin for IntelliJ IDEA *)

(* :Author: szhorvat *)
(* :Date: 2018-10-24 *)
(* :Copyright: (c) 2019 Szabolcs Horvát *)

Package["IGraphM`"]
igContextSetup[igPackagePrivateSymbol]

(***********************************)
(***** Clustering coefficients *****)
(***********************************)


PackageExport["IGGlobalClusteringCoefficient"]
IGGlobalClusteringCoefficient::usage = "IGGlobalClusteringCoefficient[graph] gives the global clustering coefficient of graph.";

Options[IGGlobalClusteringCoefficient] = { "ExcludeIsolates" -> False };
SyntaxInformation[IGGlobalClusteringCoefficient] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};
IGGlobalClusteringCoefficient[graph_?igGraphQ, OptionsPattern[]] :=
    Block[{ig = igMakeFast[graph]},
      sck@ig@"transitivityUndirected"[OptionValue["ExcludeIsolates"]]
    ]


PackageExport["IGLocalClusteringCoefficient"]
IGLocalClusteringCoefficient::usage = "IGLocalClusteringCoefficient[graph] gives the local clustering coefficient of each vertex.";

Options[IGLocalClusteringCoefficient] = { "ExcludeIsolates" -> False };
SyntaxInformation[IGLocalClusteringCoefficient] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};
(* See https://github.com/igraph/igraph/issues/907 for why SimpleGraph and UndirectedGraph are needed. *)
IGLocalClusteringCoefficient[graph_?igGraphQ, OptionsPattern[]] :=
    catch@Block[{ig = igMakeFast[SimpleGraph@UndirectedGraph[graph]]},
      fixInfNaN@check@ig@"transitivityLocalUndirected"[OptionValue["ExcludeIsolates"]]
    ]


PackageExport["IGAverageLocalClusteringCoefficient"]
IGAverageLocalClusteringCoefficient::usage = "IGAverageLocalClusteringCoefficient[graph] gives the average local clustering coefficient of graph.";

Options[IGAverageLocalClusteringCoefficient] = { "ExcludeIsolates" -> False };
SyntaxInformation[IGAverageLocalClusteringCoefficient] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};
IGAverageLocalClusteringCoefficient[graph_?igGraphQ, OptionsPattern[]] :=
    Block[{ig = igMakeFast[graph]},
      sck@ig@"transitivityAverageLocalUndirected"[OptionValue["ExcludeIsolates"]]
    ]


PackageExport["IGWeightedClusteringCoefficient"]
IGWeightedClusteringCoefficient::usage = "IGWeightedClusteringCoefficient[graph] gives the weighted local clustering coefficient, as defined by A. Barrat et al. (2004) http://dx.doi.org/10.1073/pnas.0400087101";

Options[IGWeightedClusteringCoefficient] = { "ExcludeIsolates" -> False };
SyntaxInformation[IGWeightedClusteringCoefficient] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};
IGWeightedClusteringCoefficient[graph_?igGraphQ, OptionsPattern[]] :=
    catch@Block[{ig = igMakeFastWeighted[graph]},
      fixInfNaN@check@ig@"transitivityBarrat"[OptionValue["ExcludeIsolates"]]
    ]
