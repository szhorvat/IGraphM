(* Mathematica Package *)
(* Created by Mathematica plugin for IntelliJ IDEA *)

(* :Author: szhorvat *)
(* :Date: 2018-10-24 *)
(* :Copyright: (c) 2018 Szabolcs Horvát *)

Package["IGraphM`"]
igContextSetup[igPackagePrivateSymbol]

(***********************************)
(***** Clustering coefficients *****)
(***********************************)


PackageExport["IGGlobalClusteringCoefficient"]
IGGlobalClusteringCoefficient::usage = "IGGlobalClusteringCoefficient[graph] returns the global clustering coefficient of graph.";

SyntaxInformation[IGGlobalClusteringCoefficient] = {"ArgumentsPattern" -> {_}};
IGGlobalClusteringCoefficient[graph_?igGraphQ] :=
    Block[{ig = igMakeFast[graph]}, ig@"transitivityUndirected"[]]


PackageExport["IGLocalClusteringCoefficient"]
IGLocalClusteringCoefficient::usage = "IGLocalClusteringCoefficient[graph] returns the local clustering coefficient of each vertex.";

SyntaxInformation[IGLocalClusteringCoefficient] = {"ArgumentsPattern" -> {_}};
IGLocalClusteringCoefficient[graph_?igGraphQ] :=
    Block[{ig = igMakeFast[SimpleGraph@UndirectedGraph[graph]]}, ig@"transitivityLocalUndirected"[]]


PackageExport["IGAverageLocalClusteringCoefficient"]
IGAverageLocalClusteringCoefficient::usage = "IGAverageLocalClusteringCoefficient[graph] returns the average local clustering coefficient of graph.";

SyntaxInformation[IGAverageLocalClusteringCoefficient] = {"ArgumentsPattern" -> {_}};
IGAverageLocalClusteringCoefficient[graph_?igGraphQ] :=
    Block[{ig = igMakeFast[graph]}, ig@"transitivityAverageLocalUndirected"[]]


PackageExport["IGWeightedClusteringCoefficient"]
IGWeightedClusteringCoefficient::usage = "IGWeightedClusteringCoefficient[graph] computes the weighted local clustering coefficient, as defined by A. Barrat et al. (2004) http://dx.doi.org/10.1073/pnas.0400087101";

SyntaxInformation[IGWeightedClusteringCoefficient] = {"ArgumentsPattern" -> {_}};
IGWeightedClusteringCoefficient[graph_?igGraphQ] :=
    Block[{ig = igMakeFastWeighted[graph]}, ig@"transitivityBarrat"[]]
