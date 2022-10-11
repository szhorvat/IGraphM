(* Mathematica Package *)
(* Created by Mathematica plugin for IntelliJ IDEA *)

(* :Author: szhorvat *)
(* :Date: 2018-10-24 *)
(* :Copyright: (c) 2018-2022 Szabolcs Horvát *)

Package["IGraphM`"]


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
IGTreeQ[graph_?igGraphQ, mode_String : "Out"] :=
    Block[{ig = igMakeUnweighted[graph]}, sck@ig@"treeQ"[encodeNeighborMode[mode]]]
IGTreeQ[_, _ : "Out"] := False
addCompletion[IGTreeQ, {0, {"In", "Out", "All"}}]


PackageExport["IGForestQ"]
IGForestQ::usage =
    "IGForestQ[graph] tests if graph is a forest of trees or out-trees.\n" <>
    "IGForestQ[graph, \"Out\"] tests if graph is a forest of out-trees (arborescences).\n" <>
    "IGForestQ[graph, \"In\"] tests if graph is a forest of in-trees (anti-arborescences).\n" <>
    "IGForestQ[graph, \"All\"] ignores edge directions during the test.";

SyntaxInformation[IGForestQ] = {"ArgumentsPattern" -> {_, _.}};
IGForestQ[graph_?igGraphQ, mode_String : "Out"] :=
    Block[{ig = igMakeUnweighted[graph]}, sck@ig@"forestQ"[encodeNeighborMode[mode]]]
IGForestQ[_, _ : "Out"] := False
addCompletion[IGForestQ, {0, {"In", "Out", "All"}}]


(***** Spanning trees *****)

PackageExport["IGSpanningTree"]
IGSpanningTree::usage = "IGSpanningTree[graph] gives a minimum spanning tree of graph. Edge directions are ignored. Edge weights are taken into account and are preserved in the tree.";

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
    "IGRandomSpanningTree[graph] gives a random spanning tree of graph. All spanning trees are generated with equal probability.\n" <>
    "IGRandomSpanningTree[{graph, vertex}] gives a random spanning tree of the graph component containing vertex.\n" <>
    "IGRandomSpanningTree[spec, n] gives a list of n random spanning trees.";

(* IGRandomSpanning tree has the n argument so that we can generate multiple spanning trees without
 * having to generate a new ig object (expensive) each time. *)
igRandomSpanningTree[graph_, ig_, vid_, {opt___}] :=
    With[{indices = igIndexVec@check@ig@"randomSpanningTree"[vid]},
      Graph[
        (* When a starting vertex is specified (i.e. vid >= 0), we must always include it in the vertex list.
           The C function only returns an edge list and we must handle the case when the starting vertex
           is an isolated component by itself. *)
        If[vid < 0, VertexList[graph], {VertexList[graph][[vid+1]]} ],
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
    "IGSpanningTreeCount[graph] gives the number of spanning trees of graph.\n" <>
    "IGSpanningTreeCount[graph, vertex] gives the number of spanning trees rooted in vertex for a directed graph.";

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
      i = vs1[graph][v];
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
    catch@Block[{new = igMakeEmpty[], ig = igMakeUnweighted[graph], mapping, tree},
      mapping = igIndexVec@check@new@"unfoldTree"[ManagedLibraryExpressionID[ig], vss[graph][roots], OptionValue[DirectedEdges]];
      tree = igToGraph[new];
      applyGraphOpt[opt]@Graph[tree, Properties -> Thread[VertexList[tree] -> List /@ Thread["OriginalVertex" -> igVertexNames[graph][mapping]]]]
    ]


PackageExport["IGStrahlerNumber"]
IGStrahlerNumber::usage = "IGStrahlerNumber[tree] gives the Horton–Strahler number of each vertex in a directed out-tree.";
SyntaxInformation[IGStrahlerNumber] = {"ArgumentsPattern" -> {_}};
IGStrahlerNumber[g_?igGraphQ] :=
    catch@Block[{ig = igMakeUnweighted[g]},
      check@ig@"strahlerNumber"[]
    ]
