(* Mathematica Package *)
(* Created by Mathematica Plugin for IntelliJ IDEA, see http://wlplugin.halirutan.de/ *)

(* :Author: szhorvat *)
(* :Date: 2016-06-11 *)
(* :Copyright: (c) 2018 Szabolcs Horvát *)

Package["IGraphM`"]
igContextSetup[igPackagePrivateSymbol]

(**********************************************************)
(***** Utility functions for basic graph manipulation *****)
(**********************************************************)


(***** Fast structure retrieval *****)

PackageExport["IGIndexEdgeList"]
IGIndexEdgeList::usage = "IGIndexEdgeList[graph] returns the edge list of graph in terms of vertex indices, as a packed array.";

SyntaxInformation[IGIndexEdgeList] = {"ArgumentsPattern" -> {_}};
IGIndexEdgeList[graph_?EmptyGraphQ] := {}
IGIndexEdgeList[graph_?igGraphQ] :=
    catch[1 + check@igraphGlobal@"incidenceToEdgeList"[IncidenceMatrix[graph], DirectedGraphQ[graph]]]


(***** Simple tests *****)

PackageExport["IGNullGraphQ"]
IGNullGraphQ::usage = "IGNullGraphQ[graph] tests whether graph has no vertices.";

SyntaxInformation[IGNullGraphQ] = {"ArgumentsPattern" -> {_}};
IGNullGraphQ[g_?GraphQ] := VertexCount[g] === 0
IGNullGraphQ[_] = False;


PackageExport["IGSameGraphQ"]
IGSameGraphQ::usage = "IGSameGraphQ[graph1, graph2] returns True if the given graphs have the same vertices and edges. Graph properties or edge and vertex ordering is not taken into account.";

SyntaxInformation[IGSameGraphQ] = {"ArgumentsPattern" -> {_, _}};
IGSameGraphQ[g1_?GraphQ, g2_?GraphQ] :=
    Internal`InheritedBlock[{UndirectedEdge},
      SetAttributes[UndirectedEdge, Orderless];
      Sort@VertexList[g1] === Sort@VertexList[g2] && Sort@EdgeList[g1] === Sort@EdgeList[g2]
    ]
IGSameGraphQ[_, _] := False (* return False when at least one argument is not a graph *)


(***** Adjacency list representation (also used for embeddings) *****)

PackageExport["IGAdjacencyList"]
IGAdjacencyList::usage =
    "IGAdjacencyList[graph] returns the adjacency list of graph as an association.\n" <>
    "IGAdjacencyList[graph, \"In\"] returns the adjacency list of the reverse of a directed graph.\n" <>
    "IGAdjacencyList[graph, \"All\"] considers both incoming and outgoing edges.";

igAdjacencyListSimple[vl_, am_] :=
    AssociationThread[
      vl,
      vl[[#]]& /@ am["AdjacencyLists"]
    ]

igAdjacencyListMulti[vl_, am_] :=
    GroupBy[
      Join @@ MapThread[ConstantArray, {am["NonzeroPositions"], am["NonzeroValues"]}],
      First -> Last
    ] // Map[vl[[#]]&] // KeyMap[vl[[#]]&]

SyntaxInformation[IGAdjacencyList] = {"ArgumentsPattern" -> {_, _.}};
IGAdjacencyList[graph_?GraphQ] := IGAdjacencyList[graph, "Out"]
IGAdjacencyList[graph_?EmptyGraphQ, "Out"|"In"|"All"] :=
    AssociationThread[VertexList[graph], ConstantArray[{}, VertexCount[graph]]]
(* multigraphs *)
IGAdjacencyList[graph_?MultigraphQ, "Out"] :=
    igAdjacencyListMulti[VertexList[graph], AdjacencyMatrix[graph]]
IGAdjacencyList[graph_?MultigraphQ, "In"] :=
    igAdjacencyListMulti[VertexList[graph], AdjacencyMatrix[graph]]
(* a directed graph may act like a multigraph when ignoring edge directions *)
IGAdjacencyList[graph_?DirectedGraphQ, "All"] :=
    With[{am = AdjacencyMatrix[graph]},
      igAdjacencyListMulti[VertexList[graph], am + Transpose[am]]
    ]
(* simple graphs *)
IGAdjacencyList[graph_?GraphQ, "Out"] :=
    igAdjacencyListSimple[VertexList[graph], AdjacencyMatrix[graph]]
IGAdjacencyList[graph_?GraphQ, "In"] :=
    igAdjacencyListSimple[VertexList[graph], Transpose@AdjacencyMatrix[graph]]
IGAdjacencyList[graph_?GraphQ, "All"] :=
    IGAdjacencyList[graph, "Out"] (* directed graphs were caught earlier, so this applies only to undirected ones *)


PackageExport["IGAdjacencyGraph"]
IGAdjacencyGraph::usage =
    "IGAdjacencyGraph[matrix] creates a graph from the given adjacency matrix.\n" <>
    "IGAdjacencyGraph[vertices, matrix] creates a graph with the given vertices from an adjacency matrix.\n" <>
    "IGAdjacencyGraph[adjList] creates a graph from an association representing an adjacency list.";

IGAdjacencyGraph::dir = "The adjacency list does not describe an undirected graph. Ignoring DirectedEdges -> True.";
Options[IGAdjacencyGraph] = { DirectedEdges -> Automatic };
SyntaxInformation[IGAdjacencyGraph] = {"ArgumentsPattern" -> {_, OptionsPattern[]}, "OptionNames" -> optNames[IGAdjacencyGraph, Graph]};
(* Fall back to AdjacencyGraph if the input is a matrix. *)
IGAdjacencyGraph[am_?MatrixQ, opt : OptionsPattern[]] :=
    With[{graph = AdjacencyGraph[am, opt, OptionValue[DirectedEdges]]},
      If[GraphQ[graph], graph, $Failed]
    ]
IGAdjacencyGraph[vertices_, am_?MatrixQ, opt : OptionsPattern[]] :=
    With[{graph = AdjacencyGraph[vertices, am, opt, OptionValue[DirectedEdges]]},
      If[GraphQ[graph], graph, $Failed]
    ]
(* The main purpose of this function is to handle adjacency lists stored as associations. *)
IGAdjacencyGraph[<||>, opt : OptionsPattern[]] :=
    Graph[{}, {}, opt] (* Needed because the SparseArray based implementation does not work for empty graphs. *)
expr : IGAdjacencyGraph[adjList_?AssociationQ, opt : OptionsPattern[{IGAdjacencyGraph, Graph}]] :=
    If[SubsetQ[Keys[adjList], Catenate[adjList]],
      Module[{sa, ind, symm},
        With[{sao = SystemOptions["SparseArrayOptions"]},
          Internal`WithLocalSettings[
            SetSystemOptions["SparseArrayOptions" -> "TreatRepeatedEntries" -> Total]
            ,
            ind = AssociationThread[Keys[adjList], Range@Length[adjList]];
            sa =
                SparseArray[
                  Transpose[{
                    Join @@ MapThread[ConstantArray, {Range@Length[adjList], Length /@ Values[adjList]}],
                    Lookup[ind, Catenate[adjList]]
                  }] -> 1,
                  Length[adjList] {1, 1}
                ];
            symm = Not@TrueQ@OptionValue[DirectedEdges] && SymmetricMatrixQ[sa];
            If[OptionValue[DirectedEdges] === False && Not[symm],
              Message[IGAdjacencyGraph::dir]
            ];
            If[symm,
              Graph[Keys[adjList], {Null, sa}, opt],
              Graph[Keys[adjList], {sa, Null}, opt]
            ]
            ,
            SetSystemOptions[sao]
          ]
        ]
      ]
      ,
      Message[IGAdjacencyGraph::inv, HoldForm@OutputForm[expr], adjList, "adjacency list"];
      $Failed
    ]


(***** Simple modifications *****)

PackageExport["IGReverseGraph"]
IGReverseGraph::usage = "IGReverseGraph[graph] reverses the directed edges in graph while preserving edge weights.";

SyntaxInformation[IGReverseGraph] = {"ArgumentsPattern" -> {_}};
IGReverseGraph::nmg = "Multigraphs are not currently supported.";
IGReverseGraph[g_?UndirectedGraphQ, opt : OptionsPattern[]] := Graph[g, opt]
IGReverseGraph[g_?igGraphQ, opt : OptionsPattern[]] :=
    Module[{},
      If[MultigraphQ[g],
        Message[IGReverseGraph::nmg];
        Return[$Failed]
      ];
      Graph[
        VertexList[g],
        Reverse /@ EdgeList[g],
        opt,
        Options[g, {EdgeWeight, EdgeCapacity, EdgeCost, VertexWeight, VertexCapacity}]
      ]
    ]


PackageExport["IGSimpleGraph"]
IGSimpleGraph::usage = "IGSimpleGraph[graph] converts graph to a simple graph by removing self loops and multi edges, according to the given options.";

Options[IGSimpleGraph] = { SelfLoops -> False, MultiEdges -> False, "MultipleEdges" -> "Deprecated" };
SyntaxInformation[IGSimpleGraph] = {"ArgumentsPattern" -> {_, OptionsPattern[]}, "OptionNames" -> optNames[IGSimpleGraph, Graph]};
IGSimpleGraph[g_?SimpleGraphQ, opt : OptionsPattern[{IGSimpleGraph, Graph}]] := applyGraphOpt[opt][g]
IGSimpleGraph[g_?igGraphQ, opt : OptionsPattern[{IGSimpleGraph, Graph}]] :=
    Module[{self = Not@TrueQ@OptionValue[SelfLoops], multi = Not@TrueQ@multiEdgesOptionReplace@OptionValue[MultiEdges]},
      applyGraphOpt[opt]@Which[
        self && multi, SimpleGraph[g],
        self, removeSelfLoops[g],
        multi, removeMultiEdges[g],
        True, g
      ]
    ]


PackageExport["IGUndirectedGraph"]
IGUndirectedGraph::usage = "IGUndirectedGraph[graph, conv] converts a directed graph to undirected with the given conversion method: \"Simple\" creates a single edge between connected vertices; \"All\" creates an undirected edge for each directed one and may produce a multigraph; \"Mutual\" creates a single undirected edge only between mutually connected vertices.";

(* Note: IGUndirectedGraph must ensure that:
    1. vertex names are not changed
    2. vertex ordering is not changed
 *)

SyntaxInformation[IGUndirectedGraph] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}, "OptionNames" -> optNames[Graph]};

IGUndirectedGraph[g_?UndirectedGraphQ, "Simple"|"All"|"Mutual", opt : OptionsPattern[Graph]] := Graph[g, opt]
IGUndirectedGraph[g_?igGraphQ, "Simple", opt : OptionsPattern[Graph]] :=
    Graph[VertexList[g], DeleteDuplicates[igraphGlobal@"edgeListSortPairs"@IGIndexEdgeList[g]], DirectedEdges -> False, opt]
IGUndirectedGraph[g_?igGraphQ, "All", opt : OptionsPattern[Graph]] :=
    Graph[VertexList[g], IGIndexEdgeList[g], DirectedEdges -> False, opt]
IGUndirectedGraph[g_?igGraphQ, "Mutual", opt : OptionsPattern[Graph]] :=
    With[{am = AdjacencyMatrix[g]}, (* null graph has no adjacency matrix, but this case is caught by _?UndirectedGraphQ above *)
      AdjacencyGraph[VertexList[g], Unitize[am Transpose[am]], opt]
    ]

IGUndirectedGraph[g_, "Reciprocal", opt : OptionsPattern[Graph]] := IGUndirectedGraph[g, "Mutual", opt]

IGUndirectedGraph[g_, "Each", opt : OptionsPattern[Graph]] := IGUndirectedGraph[g, "All", opt]
IGUndirectedGraph[g_, All, opt : OptionsPattern[Graph]] := IGUndirectedGraph[g, "All", opt]

IGUndirectedGraph[g_, "Collapse", opt : OptionsPattern[Graph]] := IGUndirectedGraph[g, "Simple", opt]
IGUndirectedGraph[g_, opt : OptionsPattern[Graph]] := IGUndirectedGraph[g, "Simple", opt]

addCompletion[IGUndirectedGraph, {0, {"Simple", "All", "Mutual"}}]


PackageExport["IGDirectedTree"]
IGDirectedTree::usage = "IGDirectedTree[] is deprecated. Use IGOrientTree[] instead.";
IGDirectedTree::deprec = "IGDirectedTree is deprecated and will be removed from future versions of IGraph/M. Use IGOrientTree instead.";

IGDirectedTree[tree_, root_, rest___] := IGOrientTree[tree, root, rest] (* TODO: remove eventually *)


PackageExport["IGOrientTree"]
IGOrientTree::usage =
    "IGOrientTree[tree, root] orients the edges of an undirected tree so that they point away from the root. The vertex order is not preserved: vertices will be ordered topologically.\n" <>
    "IGOrientTree[tree, root, \"In\"] orients the edges so that they point towards the root.";

SyntaxInformation[IGOrientTree] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}, "OptionNames" -> optNames[Graph]};
expr : IGOrientTree[tree_?GraphQ, root_, mode : _String : "Out", opt : OptionsPattern[Graph]] :=
    Module[{verts},
      If[Not[UndirectedGraphQ[tree] && TreeGraphQ[tree]],
        Message[IGOrientTree::inv, HoldForm@OutputForm[expr], OutputForm[tree], "undirected tree"];
        Return[$Failed]
      ];
      If[Not@VertexQ[tree, root],
        Message[IGOrientTree::inv, HoldForm@OutputForm[expr], root, "vertex"];
        Return[$Failed]
      ];
      verts = First@Last@Reap@BreadthFirstScan[tree, root, "PrevisitVertex" -> Sow];
      Switch[mode,
        "Out", Null,
        "In", verts = Reverse[verts],
        _, Message[IGOrientTree::inv, HoldForm@OutputForm[expr], mode, "mode"]; Return[$Failed]
      ];
      DirectedGraph[IGReorderVertices[verts, tree], "Acyclic", opt]
    ]
addCompletion[IGOrientTree, {0, 0, {"In", "Out"}}]


(***** Membership vectors and vertex partitions *****)

PackageExport["IGPartitionsToMembership"]
IGPartitionsToMembership::usage =
    "IGPartitionsToMembership[elements, partitions] computes a membership vector for the given partitioning of elements.\n" <>
    "IGPartitionsToMembership[graph, partitions] computes a membership vector for the given partitioning of graph's vertices.\n" <>
    "IGPartitionsToMembership[elements] is an operator that can be applied to partitions.";

IGPartitionsToMembership::invpart = "Invalid element or part specification.";
SyntaxInformation[IGPartitionsToMembership] = {"ArgumentsPattern" -> {_, _.}};
IGPartitionsToMembership[graph_?GraphQ, parts : {___List}] := IGPartitionsToMembership[VertexList[graph], parts]
IGPartitionsToMembership[elements_List, parts : {___List}] :=
    Module[{copy = parts},
      With[{union = Union @@ parts},
        If[Not[SubsetQ[elements, union] && union === Sort[Join @@ parts]],
          Message[IGPartitionsToMembership::invpart];
          Return[$Failed]
        ]
      ];
      copy[[All, All]] = Range@Length[parts];
      Lookup[AssociationThread[Flatten[parts, 1], Flatten[copy, 1]], elements, 0]
    ]
IGPartitionsToMembership[elements_][parts_] := IGPartitionsToMembership[elements, parts]


PackageExport["IGMembershipToPartitions"]
IGMembershipToPartitions::usage =
    "IGMembershipToPartitions[elements, membership] computes a partitioning of elements based on the given membership vector.\n" <>
    "IGMembershipToPartitions[graph, membership] computes a partitioning graph's vertices based on the given membership vector.\n" <>
    "IGMembershipToPartitions[elements] is an operator that can be applied to membership vectors.";

SyntaxInformation[IGMembershipToPartitions] = {"ArgumentsPattern" -> {_, _.}};
IGMembershipToPartitions[graph_?GraphQ, membership_List] := IGMembershipToPartitions[VertexList[graph], membership]
IGMembershipToPartitions[elements_List, membership_List] :=
    If[Length[elements] == Length[membership],
      Values@GroupBy[Transpose[{elements, membership}], Last -> First]
      ,
      Message[IGMembershipToPartitions::ndims, elements, membership];
      $Failed
    ]
IGMembershipToPartitions[elements_][membership_] := IGMembershipToPartitions[elements, membership]


(***** Graph combination *****)

PackageExport["IGDisjointUnion"]
IGDisjointUnion::usage = "IGDisjointUnion[{g1, g2, \[Ellipsis]}] computes a disjoint union of the graphs. Each vertex of the result will be a pair consisting of the index of the graph originally containing it and the original name of the vertex.";

IGDisjointUnion::mixed = "IGDisjointUnion does not support mixed graphs.";
SyntaxInformation[IGDisjointUnion] = {"ArgumentsPattern" -> {_, OptionsPattern[]}, "OptionNames" -> optNames[IGDisjointUnion, Graph]};
IGDisjointUnion[{} | <||>, opt : OptionsPattern[]] := Graph[{}, {}, opt]
IGDisjointUnion[glist : {__?UndirectedGraphQ}, opt : OptionsPattern[]] :=
    igDisjointUnion[Range@Length[glist], glist, False, {opt}]
IGDisjointUnion[glist : {__?DirectedGraphQ}, opt : OptionsPattern[]] :=
    igDisjointUnion[Range@Length[glist], glist, True, {opt}]
IGDisjointUnion[gasc_ /; AssociationQ[gasc] && AllTrue[gasc, UndirectedGraphQ], opt : OptionsPattern[]] :=
    igDisjointUnion[Keys[gasc], Values[gasc], False, {opt}]
IGDisjointUnion[gasc_ /; AssociationQ[gasc] && AllTrue[gasc, DirectedGraphQ], opt : OptionsPattern[]] :=
    igDisjointUnion[Keys[gasc], Values[gasc], True, {opt}]
IGDisjointUnion[glist : {__?GraphQ}, opt : OptionsPattern[]] :=
    (Message[IGDisjointUnion::mixed]; $Failed)
IGDisjointUnion[gasc_ /; AssociationQ[gasc] && AllTrue[gasc, GraphQ], opt : OptionsPattern[]] :=
    (Message[IGDisjointUnion::mixed]; $Failed)

igDisjointUnion[gnames_, glist_, directed_, {opt___}] :=
    With[{vc = VertexCount /@ glist},
      Graph[
        Join @@ MapThread[Tuples@*List, {List /@ gnames, VertexList /@ glist}],
        Join @@ ((IGIndexEdgeList /@ glist) + FoldList[Plus, 0, Most[vc]]),
        DirectedEdges -> directed,
        opt
      ]
    ]

(***** Change graph representations *****)

PackageExport["IGReorderVertices"]
IGReorderVertices::usage = "IGReorderVertices[vertices, graph] reorders the vertices of graph according to the given vertex vector. Graph properties are preserved.";

findPermutation::usage = "findPermutation[l1, l2] finds the permutation that transforms list l1 into list l2.";
findPermutation[l1_, l2_] := Ordering[l1][[Ordering@Ordering[l2]]]

IGReorderVertices::bdvert = "The provided vertex list differs from the vertex list of the graph.";

SyntaxInformation[IGReorderVertices] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}, "OptionNames" -> optNames[Graph]};
IGReorderVertices[verts_List, graph_?GraphQ, opt : OptionsPattern[]] :=
    Module[{perm, vl = VertexList[graph]},
      If[Sort[verts] =!= Sort[vl],
        Message[IGReorderVertices::bdvert];
        Return[$Failed]
      ];
      perm = findPermutation[vl, verts];
      Graph[
        verts,
        EdgeList[graph],
        opt,
        Replace[
          Options[graph],
          Verbatim[Rule][sym : VertexWeight|VertexCapacity|VertexCoordinates, val_List] :> Rule[sym, val[[perm]]],
          {1}
        ]
      ]
    ]


(***** Transfer properties to a smaller graph *****)

PackageExport["IGTake"]

IGTake::usage = "IGTake[] is deprecated. Use IGTakeSubgraph[] instead.";
IGTake::deprec = "IGTake is deprecated and will be removed from future versions of IGraph/M. Use IGTakeSubgraph instead."

IGTake[args___] := (Message[IGTake::deprec]; IGTakeSubgraph[args]) (* TODO: remove eventually *)


PackageExport["IGTakeSubgraph"]
IGTakeSubgraph::usage =
    "IGTakeSubgraph[graph, subgraph] keeps only those vertices and edges of graph which are also present in subgraph, while retaining all graph properties.\n" <>
    "IGTakeSubgraph[graph, edges] uses an edge list as the subgraph specification.";

(* Keep those elements in l1 which are also in l2. l1 may have repeating elements. *)
keepCases[l1_, l2_] := Pick[l1, Lookup[AssociationThread[l2, ConstantArray[1, Length[l2]]], l1, 0], 1]

allPropNames = {
  Properties,
  VertexWeight, VertexCapacity, VertexCoordinates,
  VertexSize, VertexShape, VertexShapeFunction, VertexStyle, VertexLabels, VertexLabelStyle,
  EdgeWeight, EdgeCapacity, EdgeCost,
  EdgeStyle, EdgeShapeFunction, EdgeLabels, EdgeLabelStyle
};

IGTakeSubgraph::nsg = "Some of the edges or vertices of the second graph are not present in the first.";

Options[IGTakeSubgraph] = Options[Graph];

SyntaxInformation[IGTakeSubgraph] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};

IGTakeSubgraph[g_?GraphQ, edges_List, opt : OptionsPattern[]] :=
    Module[{sg},
      sg = Quiet@Graph[edges, opt];
      IGTakeSubgraph[g, sg, opt] /; GraphQ[sg]
    ]

IGTakeSubgraph[g_?GraphQ, sg_?GraphQ, opt : OptionsPattern[]] :=
    Internal`InheritedBlock[{UndirectedEdge}, (* ensure that a <-> b compares equal to b <-> a *)
      SetAttributes[UndirectedEdge, Orderless];
      Module[{options, prop, vindex, eindex, sgEdgeList, vlist, elist, handleList, handleRule},
        sgEdgeList = DeleteDuplicates@EdgeList[sg];

        (* Check that sg is contained within g (ignores edge multiplicities). *)
        If[Not[SubsetQ[VertexList[g], VertexList[sg]] && SubsetQ[EdgeList[g], sgEdgeList]],
          Message[IGTakeSubgraph::nsg];
          Return[$Failed]
        ];

        (* These options contain the graph properties, but some options are for different purposes *)
        options = Association@Options[g];

        (* Used to find vertex/edge indices *)
        vindex = AssociationThread[VertexList[g], Range@VertexCount[g]];
        eindex = PositionIndex@EdgeList[g]; (* edge multiplicities must be handled *)

        vlist = VertexList[sg];
        elist = keepCases[EdgeList[g], sgEdgeList];

        (* Handle custom properties *)
        prop = If[KeyExistsQ[options, Properties],
          Properties ->
              Normal@KeyTake[
                options[Properties],
                Join[VertexList[sg], elist, {"DefaultEdgeProperties", "DefaultVertexProperties", "GraphProperties"}]
              ]
          ,
          Unevaluated@Sequence[]
        ];

        (* Handle List-type properties.  We preserve those with indices idx. *)
        handleList[idx_][propName_] :=
            If[KeyExistsQ[options, propName],
              propName -> options[propName][[idx]]
              ,
              Unevaluated@Sequence[]
            ];

        (* Handle Rule-type properties.  elems can be a list of vertices or edges. *)
        handleRule[elems_][propName_] :=
            If[KeyExistsQ[options, propName],
              Module[{rules = options[propName], default = {}},
                If[Not@ruleQ@First[rules],
                  default = {First[rules]};
                  rules = Rest[rules];
                ];
                propName -> Join[default, Normal@KeyTake[rules, elems]]
              ]
              ,
              Unevaluated@Sequence[]
            ];

        Graph[
          Graph[vlist, elist,
            prop,
            handleList[Lookup[vindex, VertexList[sg]]] /@ {VertexWeight, VertexCapacity, VertexCoordinates},
            handleList[Flatten@Lookup[eindex, elist]] /@ {EdgeWeight, EdgeCapacity, EdgeCost},
            handleRule[vlist] /@ {VertexSize, VertexShape, VertexShapeFunction, VertexStyle, VertexLabels, VertexLabelStyle},
            handleRule[elist] /@ {EdgeStyle, EdgeShapeFunction, EdgeLabels, EdgeLabelStyle},
            Normal@KeyDrop[options, allPropNames]
          ],
          opt] (* apply user options *)
      ]
    ]
