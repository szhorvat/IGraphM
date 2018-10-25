(* Mathematica Package *)
(* Created by Mathematica plugin for IntelliJ IDEA *)

(* :Author: szhorvat *)
(* :Date: 2018-10-24 *)
(* :Copyright: (c) 2018 Szabolcs Horvát *)

Package["IGraphM`"]
igContextSetup[igPackagePrivateSymbol]

(******************************************)
(***** Tree graphs and spanning trees *****)
(******************************************)


(***** Testing *****)

PackageExport["IGTreeQ"]
IGTreeQ::usage =
    "IGTreeQ[graph] tests if graph is a tree or out-tree.\n" <>
    "IGTreeQ[graph, \"Out\"] tests if graph is an out-tree (arborescence).\n" <>
    "IGTreeQ[graph, \"In\"] tests if graph is an in-tree (anti-arborescence).\n" <>
    "IGTreeQ[graph, \"All\"] ignores edge directions during the test.";

SyntaxInformation[IGTreeQ] = {"ArgumentsPattern" -> {_, _.}};
IGTreeQ[graph_?igGraphQ, mode_ : "Out"] :=
    Block[{ig = igMakeFast[graph]}, sck@ig@"treeQ"[Lookup[<|"Out" -> 1, "In" -> 2, "All" -> 3|>, mode, -1]]]
IGTreeQ[_] := False (* provide this only for the single argument form *)
addCompletion[IGTreeQ, {0, {"In", "Out", "All"}}]


(***** Spanning trees *****)

PackageExport["IGSpanningTree"]
IGSpanningTree::usage = "IGSpanningTree[graph] returns a minimum spanning tree of graph. Edge directions are ignored. Edge weights are taken into account and are preserved in the tree.";

SyntaxInformation[IGSpanningTree] = {"ArgumentsPattern" -> {_, OptionsPattern[]}, "OptionNames" -> optNames[Graph]};
IGSpanningTree[graph_?igGraphQ, opt : OptionsPattern[]] :=
    catch@Block[{ig = igMake[graph]},
      With[{indices = igIndexVec@check@ig@"spanningTree"[]},
        Graph[
          VertexList[graph],
          EdgeList[graph][[indices]],
          If[IGEdgeWeightedQ[graph], EdgeWeight -> igEdgeWeights[graph][[indices]], Unevaluated@Sequence[]],
          opt
        ]
      ]
    ]


PackageExport["IGRandomSpanningTree"]
IGRandomSpanningTree::usage =
    "IGRandomSpanningTree[graph] returns a random spanning tree of graph. All spanning trees are generated with equal probability.\n" <>
    "IGRandomSpanningTree[{graph, vertex}] returns a random spanning tree of the graph component containing vertex.\n" <>
    "IGRandomSpanningTree[spec, n] returns a list of n random spanning trees.";

(* IGRandomSpanning tree has the n argument so that we can generate multiple spanning trees without
 * having to generate a new ig object (expensive) each time. *)
igRandomSpanningTree[graph_, ig_, vid_, {opt___}] :=
    With[{indices = igIndexVec@check@ig@"randomSpanningTree"[vid]},
      Graph[
        If[vid < 0, VertexList[graph], Unevaluated@Sequence[]],
        EdgeList[graph][[indices]],
        opt
      ]
    ]

SyntaxInformation[IGRandomSpanningTree] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}, "OptionNames" -> optNames[Graph]};
IGRandomSpanningTree[graph_?igGraphQ, opt : OptionsPattern[]] :=
    catch@Block[{ig = igMake[graph]},
      igRandomSpanningTree[graph, ig, -1, {opt}]
    ]
IGRandomSpanningTree[{graph_?igGraphQ, v_}, opt : OptionsPattern[]] :=
    catch@Block[{ig = igMake[graph]},
      igRandomSpanningTree[graph, ig, vs[graph][v], {opt}]
    ]
IGRandomSpanningTree[graph_?igGraphQ, n_?Internal`NonNegativeIntegerQ, opt : OptionsPattern[]] :=
    catch@Block[{ig = igMake[graph]},
      Table[igRandomSpanningTree[graph, ig, -1, {opt}], {n}]
    ]
IGRandomSpanningTree[{graph_?igGraphQ, v_}, n_?Internal`NonNegativeIntegerQ, opt : OptionsPattern[]] :=
    catch@Block[{ig = igMake[graph]},
      Table[igRandomSpanningTree[graph, ig, vs[graph][v], {opt}], {n}]
    ]


PackageExport["IGSpanningTreeCount"]
IGSpanningTreeCount::usage =
    "IGSpanningTreeCount[graph] returns the number of spanning trees of graph.\n" <>
    "IGSpanningTreeCount[graph, vertex] returns the number of spanning trees rooted in vertex for a directed graph.";

(* Notes:
 * The null graph has no spanning trees.
 * The singleton graph has one spanning tree. *)
SyntaxInformation[IGSpanningTreeCount] = {"ArgumentsPattern" -> {_, _.}};
IGSpanningTreeCount[graph_?GraphQ /; VertexCount[graph] < 2] := VertexCount[graph]
IGSpanningTreeCount[graph_?UndirectedGraphQ] := Det@Rest@Transpose@Rest@IGKirchhoffMatrix[graph]
IGSpanningTreeCount[graph_?DirectedGraphQ] :=
    With[{km = IGKirchhoffMatrix[graph, "In"]},
      Total@Table[Det@Delete[Transpose@Delete[km, i], i], {i, VertexCount[graph]}]
    ]
IGSpanningTreeCount[graph_?GraphQ, v_] :=
    catch@Module[{i, km},
      Check[i = VertexIndex[graph, v], throw[$Failed]];
      If[VertexCount[graph] == 1,
        1,
        km = IGKirchhoffMatrix[graph, "In"];
        Det@Delete[Transpose@Delete[km, i], i]
      ]
    ]


(***** Other tree-related functions *****)

PackageExport["IGUnfoldTree"]
IGUnfoldTree::usage = "IGUnfoldTree[graph, {root1, root2, \[Ellipsis]}] performs a breadth-first search on graph starting from the given roots, and converts it to a tree or forest by replicating vertices that were found more than once. The original vertex that generated a tree node is stored in the \"OriginalVertex\" property.";

Options[IGUnfoldTree] = { DirectedEdges -> True };
SyntaxInformation[IGUnfoldTree] = {"ArgumentsPattern" -> {_, OptionsPattern[]}, "OptionNames" -> optNames[IGUnfoldTree, Graph]};
IGUnfoldTree[graph_?GraphQ, roots_List, opt : OptionsPattern[{IGUnfoldTree, Graph}]] :=
    catch@Block[{new = igMakeEmpty[], ig = igMakeFast[graph], mapping, tree},
      mapping = check@new@"unfoldTree"[ManagedLibraryExpressionID[ig], vss[graph][roots], OptionValue[DirectedEdges]];
      tree = igToGraph[new];
      applyGraphOpt[opt]@Graph[tree, Properties -> Thread[VertexList[tree] -> List /@ Thread["OriginalVertex" -> igVertexNames[graph]@igIndexVec[mapping]]]]
    ]