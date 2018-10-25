
(* Mathematica Package *)
(* Created by Mathematica plugin for IntelliJ IDEA *)

(* :Author: szhorvat *)
(* :Date: 2018-10-24 *)
(* :Copyright: (c) 2018 Szabolcs Horvát *)

Package["IGraphM`"]
igContextSetup[igPackagePrivateSymbol]

(***********************************************)
(***** Cliques and independent vertex sets *****)
(***********************************************)


(***** Cliques *****)

PackageExport["IGCliques"]
IGCliques::usage =
    "IGCliques[graph] returns all complete subgraphs (cliques) in graph. Note that this is different from the builtin FindCliques[], which finds maximal cliques.\n" <>
    "IGCliques[graph, {min, max}] returns all complete subgraphs between sizes min and max.\n" <>
    "IGCliques[graph, max] returns all complete subgraphs of size at most max.\n" <>
    "IGCliques[graph, {n}] returns all complete subgraphs of size n.";

SyntaxInformation[IGCliques] = {"ArgumentsPattern" -> {_, _.}};
IGCliques[graph_] := IGCliques[graph, Infinity]
IGCliques[graph_, max : (_Integer | Infinity)] := IGCliques[graph, {1, max}]
IGCliques[graph_, {size_}] := IGCliques[graph, {size, size}]
IGCliques[graph_?igGraphQ, {min_?Internal`PositiveMachineIntegerQ, max : (_?Internal`PositiveMachineIntegerQ | Infinity)}] /; max >= min :=
    igCliques[graph, {min, infToZero[max]}]
igCliques[graph_, {min_, max_}] :=
    catch@Block[{ig = igMakeFast[graph]},
      igUnpackVertexSet[graph]@check@ig@"cliques"[min, max]
    ]


PackageExport["IGCliqueSizeCounts"]
IGCliqueSizeCounts::usage =
    "IGCliqueSizeCounts[graph] computes a histogram of clique sizes in graph. The kth element of the result is the number of k-cliques.\n" <>
    "IGCliqueSizeCounts[graph, {min, max}] computes a histogram of clique sizes between min and max in graph.\n" <>
    "IGCliqueSizeCounts[graph, max] computes a histogram of clique sizes no larger than max in graph.\n" <>
    "IGCliqueSizeCounts[graph, {n}] counts cliques of size n in graph.";

SyntaxInformation[IGCliqueSizeCounts] = {"ArgumentsPattern" -> {_, _.}};
IGCliqueSizeCounts[graph_] := IGCliqueSizeCounts[graph, Infinity]
IGCliqueSizeCounts[graph_, max : (_Integer | Infinity)] := IGCliqueSizeCounts[graph, {1, max}]
IGCliqueSizeCounts[graph_, {size_}] := IGCliqueSizeCounts[graph, {size, size}]
IGCliqueSizeCounts[graph_?igGraphQ, {min_?Internal`PositiveMachineIntegerQ, max : (_?Internal`PositiveMachineIntegerQ | Infinity)}] /; max >= min :=
    igCliqueSizeCounts[graph, {min, infToZero[max]}]
igCliqueSizeCounts[graph_, {min_, max_}] :=
    catch@Block[{ig = igMakeFast[graph]},
      Round@check@ig@"cliqueDistribution"[min, max]
    ]


PackageExport["IGMaximalCliqueSizeCounts"]
IGMaximalCliqueSizeCounts::usage =
    "IGMaximalCliqueSizeCounts[graph] computes a histogram of maximal clique sizes in graph. The kth element of the result is the number of maximal k-cliques.\n" <>
    "IGMaximalCliqueSizeCounts[graph, {min, max}] computes a histogram of maximal clique sizes between min and max in graph.\n" <>
    "IGMaximalCliqueSizeCounts[graph, max] computes a histogram of maximal clique sizes no larger than max in graph.\n" <>
    "IGMaximalCliqueSizeCounts[graph, {n}] counts maximal cliques of size n in graph.";

SyntaxInformation[IGMaximalCliqueSizeCounts] = {"ArgumentsPattern" -> {_, _.}};
IGMaximalCliqueSizeCounts[graph_] := IGMaximalCliqueSizeCounts[graph, Infinity]
IGMaximalCliqueSizeCounts[graph_, max : (_Integer | Infinity)] := IGMaximalCliqueSizeCounts[graph, {1, max}]
IGMaximalCliqueSizeCounts[graph_, {size_}] := IGMaximalCliqueSizeCounts[graph, {size, size}]
IGMaximalCliqueSizeCounts[graph_?igGraphQ, {min_?Internal`PositiveMachineIntegerQ, max : (_?Internal`PositiveMachineIntegerQ | Infinity)}] /; max >= min :=
    igMaximalCliqueSizeCounts[graph, {min, infToZero[max]}]
igMaximalCliqueSizeCounts[graph_, {min_, max_}] :=
    catch@Block[{ig = igMakeFast[graph]},
      Round@check@ig@"maximalCliqueDistribution"[min, max]
    ]


PackageExport["IGMaximalCliques"]
IGMaximalCliques::usage =
    "IGMaximalCliques[graph] returns all maximal cliques in graph.\n" <>
    "IGMaximalCliques[graph, {min, max}] returns all maximal cliques between sizes min and max.\n" <>
    "IGMaximalCliques[graph, max] returns all maximal cliques of size at most max.\n" <>
    "IGMaximalCliques[graph, {n}] returns all maximal cliques of size n.";

SyntaxInformation[IGMaximalCliques] = {"ArgumentsPattern" -> {_, _.}};
IGMaximalCliques[graph_] := IGMaximalCliques[graph, Infinity]
IGMaximalCliques[graph_, max : (_Integer | Infinity)] := IGMaximalCliques[graph, {1, max}]
IGMaximalCliques[graph_, {size_}] := IGMaximalCliques[graph, {size, size}]
IGMaximalCliques[graph_?igGraphQ, {min_?Internal`PositiveMachineIntegerQ, max : (_?Internal`PositiveMachineIntegerQ | Infinity)}] /; max >= min :=
    igMaximalCliques[graph, {min, infToZero[max]}]
igMaximalCliques[graph_, {min_, max_}] :=
    catch@Block[{ig = igMakeFast[graph]},
      igUnpackVertexSet[graph]@check@ig@"maximalCliques"[min, max]
    ]


PackageExport["IGMaximalCliquesCount"]
IGMaximalCliquesCount::usage =
    "IGMaximalCliquesCount[graph] counts all maximal cliques in graph.\n" <>
    "IGMaximalCliquesCount[graph, {min, max}] counts all maximal cliques between sizes min and max.\n" <>
    "IGMaximalCliquesCount[graph, max] counts all maximal cliques of size at most max.\n" <>
    "IGMaximalCliquesCount[graph, {n}] counts all maximal cliques of size n.";

SyntaxInformation[IGMaximalCliquesCount] = {"ArgumentsPattern" -> {_, _.}};
IGMaximalCliquesCount[graph_] := IGMaximalCliquesCount[graph, Infinity]
IGMaximalCliquesCount[graph_, max : (_Integer | Infinity)] := IGMaximalCliquesCount[graph, {1, max}]
IGMaximalCliquesCount[graph_, {size_}] := IGMaximalCliquesCount[graph, {size, size}]
IGMaximalCliquesCount[graph_?igGraphQ, {min_?Internal`PositiveMachineIntegerQ, max : (_?Internal`PositiveMachineIntegerQ | Infinity)}] /; max >= min :=
    igMaximalCliquesCount[graph, {min, infToZero[max]}]
igMaximalCliquesCount[graph_, {min_, max_}] :=
    catch@Block[{ig = igMakeFast[graph]},
      Round@check@ig@"maximalCliquesCount"[min, max]
    ]


PackageExport["IGLargestCliques"]
IGLargestCliques::usage = "IGLargestCliques[graph] returns the largest cliques in graph.";

SyntaxInformation[IGLargestCliques] = {"ArgumentsPattern" -> {_}};
IGLargestCliques[graph_?igGraphQ] :=
    catch@Block[{ig = igMakeFast[graph]},
      igUnpackVertexSet[graph]@check@ig@"largestCliques"[]
    ]


PackageExport["IGCliqueNumber"]
IGCliqueNumber::usage = "IGCliqueNumber[graph] returns the clique number of graph. The clique number is the size of the largest clique.";

SyntaxInformation[IGCliqueNumber] = {"ArgumentsPattern" -> {_}};
IGCliqueNumber[graph_?igGraphQ] := Block[{ig = igMakeFast[graph]}, sck@ig@"cliqueNumber"[]]


PackageExport["IGWeightedCliques"]
IGWeightedCliques::usage = "IGWeightedCliques[graph, {min, max}] returns all complete subgraphs having total vertex weight between min and max.";

SyntaxInformation[IGWeightedCliques] = {"ArgumentsPattern" -> {_, {_, _}}};
IGWeightedCliques[graph_?igGraphQ, {min_?Internal`NonNegativeMachineIntegerQ, max : (_?Internal`PositiveMachineIntegerQ | Infinity)}] /; max >= min :=
    If[IGVertexWeightedQ[graph], igCliquesWeighted, igCliques][graph, {min, infToZero[max]}]
igCliquesWeighted[graph_, {min_, max_}] :=
    catch@Block[{ig = igMakeFast[graph]},
      igUnpackVertexSet[graph]@check@ig@"cliquesWeighted"[min, max, igVertexWeights[graph], False]
    ]


PackageExport["IGMaximalWeightedCliques"]
IGMaximalWeightedCliques::usage = "IGMaximalWeightedCliques[graph, {min, max}] returns all maximal cliques having total vertex weight between min and max.";

SyntaxInformation[IGMaximalWeightedCliques] = {"ArgumentsPattern" -> {_, {_, _}}};
IGMaximalWeightedCliques[graph_?igGraphQ, {min_?Internal`NonNegativeMachineIntegerQ, max : (_?Internal`PositiveMachineIntegerQ | Infinity)}] /; max >= min :=
    If[IGVertexWeightedQ[graph], igMaximalCliquesWeighted, igMaximalCliques][graph, {min, infToZero[max]}]
igMaximalCliquesWeighted[graph_, {min_, max_}] :=
    catch@Block[{ig = igMakeFast[graph]},
      igUnpackVertexSet[graph]@check@ig@"cliquesWeighted"[min, max, igVertexWeights[graph], True (* maximal *)]
    ]


PackageExport["IGLargestWeightedCliques"]
IGLargestWeightedCliques::usage = "IGLargestWeightedCliques[graph] returns the cliques having the largest total vertex weight in graph.";

SyntaxInformation[IGLargestWeightedCliques] = {"ArgumentsPattern" -> {_}};
IGLargestWeightedCliques[graph_?igGraphQ] :=
    If[IGVertexWeightedQ[graph], igLargestCliquesWeighted, IGLargestCliques][graph]
igLargestCliquesWeighted[graph_] :=
    catch@Block[{ig = igMakeFast[graph]},
      igUnpackVertexSet[graph]@check@ig@"largestCliquesWeighted"[igVertexWeights[graph]]
    ]


PackageExport["IGWeightedCliqueNumber"]
IGWeightedCliqueNumber::usage = "IGWeightedCliqueNumber[graph] return the maximum total vertex weight of any clique in graph.";

SyntaxInformation[IGWeightedCliqueNumber] = {"ArgumentsPattern" -> {_}};
IGWeightedCliqueNumber[graph_?igGraphQ] :=
    If[IGVertexWeightedQ[graph], igMaximumCliqueWeight, IGCliqueNumber][graph]
igMaximumCliqueWeight[graph_] :=
    Block[{ig = igMakeFast[graph]},
      sck@ig@"cliqueNumberWeighted"[igVertexWeights[graph]]
    ]



(***** Independent vertex sets *****)


PackageExport["IGIndependentVertexSets"]
IGIndependentVertexSets::usage =
    "IGIndependentVertexSets[graphs] returns all independent vertex sets of graph.\n" <>
    "IGIndependentVertexSets[graphs, {min, max}] returns all independent vertex sets of graph between sizes min and max.\n" <>
    "IGIndependentVertexSets[graphs, max] returns all independent vertex sets up to size max.\n" <>
    "IGIndependentVertexSets[graphs, {n}] returns all independent vertex sets of size n.";

SyntaxInformation[IGIndependentVertexSets] = {"ArgumentsPattern" -> {_, _.}};
IGIndependentVertexSets[graph_] := IGIndependentVertexSets[graph, Infinity]
IGIndependentVertexSets[graph_, max : (_Integer | Infinity)] := IGIndependentVertexSets[graph, {1, max}]
IGIndependentVertexSets[graph_, {size_}] := IGIndependentVertexSets[graph, {size, size}]
IGIndependentVertexSets[graph_?igGraphQ, {min_?Internal`PositiveMachineIntegerQ, max : (_?Internal`PositiveMachineIntegerQ | Infinity)}] /; max >= min :=
    igIndependentVertexSets[graph, {min, infToZero[max]}]
igIndependentVertexSets[graph_, {min_, max_}] :=
    catch@Block[{ig = igMakeFast[graph]},
      igUnpackVertexSet[graph]@check@ig@"independentVertexSets"[min, max]
    ]


PackageExport["IGLargestIndependentVertexSets"]
IGLargestIndependentVertexSets::usage = "IGLargestIndependentVertexSets[graph] finds the largest independent vertex sets of graph.";

SyntaxInformation[IGLargestIndependentVertexSets] = {"ArgumentsPattern" -> {_}};
IGLargestIndependentVertexSets[graph_?igGraphQ] :=
    catch@Block[{ig = igMakeFast[graph]},
      igUnpackVertexSet[graph]@check@ig@"largestIndependentVertexSets"[]
    ]


PackageExport["IGMaximalIndependentVertexSets"]
IGMaximalIndependentVertexSets::usage = "IGMaximalIndependentVertexSets[graph] finds the maximal independent vertex sets of graph.";

SyntaxInformation[IGMaximalIndependentVertexSets] = {"ArgumentsPattern" -> {_}};
IGMaximalIndependentVertexSets[graph_?igGraphQ] :=
    catch@Block[{ig = igMakeFast[graph]},
      igUnpackVertexSet[graph]@check@ig@"maximalIndependentVertexSets"[]
    ]


PackageExport["IGIndependenceNumber"]
IGIndependenceNumber::usage = "IGIndependenceNumber[graph] returns the independence number of graph. The independence number is the size of the largest independent vertex set.";

SyntaxInformation[IGIndependenceNumber] = {"ArgumentsPattern" -> {_}};
IGIndependenceNumber[graph_?igGraphQ] := Block[{ig = igMakeFast[graph]}, sck@ig@"independenceNumber"[]]

