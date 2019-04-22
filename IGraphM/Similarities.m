(* Mathematica Package *)
(* Created by Mathematica plugin for IntelliJ IDEA *)

(* :Author: szhorvat *)
(* :Date: 2018-10-24 *)
(* :Copyright: (c) 2019 Szabolcs Horvát *)

Package["IGraphM`"]
igContextSetup[igPackagePrivateSymbol]

(*******************************)
(***** Similarity measures *****)
(*******************************)


(* for those that return a list of vectors *)
similarityFunction1[name_, post_ : Identity][graph_, All] :=
    catch@Block[{ig = igMakeFast[graph]}, post@check@ig@name[{}] ]
similarityFunction1[name_, post_ : Identity][graph_, {}] := {}
similarityFunction1[name_, post_ : Identity][graph_, vs_List] :=
    catch@Block[{ig = igMakeFast[graph]},
      post@check@ig@name[ Check[VertexIndex[graph, #] - 1& /@ vs, Return[$Failed, Block]] ]
    ]
similarityFunction1[name_, post_ : Identity][graph_, v_] := similarityFunction1[name, First @* post][graph, {v}]


PackageExport["IGCocitationCoupling"]
IGCocitationCoupling::usage =
    "IGCocitationCoupling[graph] returns the cocitation coupling between all vertex pairs in graph. The cocitation coupling of two vertices is the number of vertices connecting to both of them (with directed edges).\n" <>
    "IGCocitationCoupling[graph, vertex] returns the cocitation coupling of vertex with all other vertices in graph.\n" <>
    "IGCocitationCoupling[graph, {vertex1, vertex2, \[Ellipsis]}] returns the cocitation coupling of vertex1, vertex2, \[Ellipsis] with all other vertices in graph.";

SyntaxInformation[IGCocitationCoupling] = {"ArgumentsPattern" -> {_, _.}};
IGCocitationCoupling[graph_?igGraphQ, vs_ : All] := similarityFunction1["similarityCocitation", Round][graph, vs]


PackageExport["IGBibliographicCoupling"]
IGBibliographicCoupling::usage =
    "IGBibliographicCoupling[graph] returns the bibliographic coupling between all vertex pairs in graph. The bibliographic coupling of two vertices is the number of vertices they both connect to (with directed edges).\n" <>
    "IGBibliographicCoupling[graph, vertex] returns the bibliographic coupling of vertex with all other vertices in graph.\n" <>
    "IGBibliographicCoupling[graph, {vertex1, vertex2, \[Ellipsis]}] returns the bibliographic coupling of vertex1, vertex2, \[Ellipsis] with all other vertices in graph.";

SyntaxInformation[IGBibliographicCoupling] = {"ArgumentsPattern" -> {_, _.}};
IGBibliographicCoupling[graph_?igGraphQ, vs_ : All] := similarityFunction1["similarityBibcoupling", Round][graph, vs]


PackageExport["IGInverseLogWeightedSimilarity"]
IGInverseLogWeightedSimilarity::usage =
    "IGInverseLogWeightedSimilarity[graph] returns the inverse log-weighted similarity between all pairs of vertices.\n" <>
    "IGInverseLogWeightedSimilarity[graph, vertex] returns the inverse log-weighted similarity of vertex to all other vertices.\n" <>
    "IGInverseLogWeightedSimilarity[graph, {vertex1, vertex2, \[Ellipsis]}] returns the inverse log-weighted similarity between the given vertices.";

SyntaxInformation[IGInverseLogWeightedSimilarity] = {"ArgumentsPattern" -> {_, _.}};
IGInverseLogWeightedSimilarity[graph_?igGraphQ, vs_ : All] := similarityFunction1["similarityInverseLogWeighted"][graph, vs]


(* for those that return a matrix *)
similarityFunction2[name_][graph_, All, loops_] :=
    catch@Block[{ig = igMakeFast[graph]}, check@ig@name[{}, loops] ]
similarityFunction2[name_][graph_, {}, loops_] := {}
similarityFunction2[name_][graph_, vs_List, loops_] :=
    catch@Block[{ig = igMakeFast[graph]},
      check@ig@name[ Check[VertexIndex[graph, #] - 1& /@ vs, Return[$Failed, Block]], loops ]
    ]


PackageExport["IGJaccardSimilarity"]
IGJaccardSimilarity::usage =
    "IGJaccardSimilarity[graph] returns the Jaccard similarity between all pairs of vertices.\n" <>
    "IGJaccardSimilarity[graph, {vertex1, vertex2, \[Ellipsis]}] returns the Jaccard similarity between the given vertices.";

Options[IGJaccardSimilarity] = { SelfLoops -> False };
SyntaxInformation[IGJaccardSimilarity] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};
IGJaccardSimilarity[graph_?igGraphQ, vs : (_List | All) : All, opt : OptionsPattern[]] := similarityFunction2["similarityJaccard"][graph, vs, OptionValue[SelfLoops]]


PackageExport["IGDiceSimilarity"]
IGDiceSimilarity::usage =
    "IGDiceSimilarity[graph] returns the Dice similarity between all pairs of vertices.\n" <>
    "IGDiceSimilarity[graph, {vertex1, vertex2, \[Ellipsis]}] returns the Dice similarity between the given vertices.";

Options[IGDiceSimilarity] = { SelfLoops -> False };
SyntaxInformation[IGDiceSimilarity] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};
IGDiceSimilarity[graph_?igGraphQ, vs : (_List | All) : All, opt : OptionsPattern[]] := similarityFunction2["similarityDice"][graph, vs, OptionValue[SelfLoops]]
