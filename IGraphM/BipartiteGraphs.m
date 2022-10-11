(* Mathematica Package *)
(* Created by Mathematica plugin for IntelliJ IDEA *)

(* :Author: szhorvat *)
(* :Date: 2018-10-24 *)
(* :Copyright: (c) 2019-2022 Szabolcs Horvát *)

Package["IGraphM`"]


(****************************)
(***** Bipartite graphs *****)
(****************************)


PackageExport["IGBipartiteQ"]
IGBipartiteQ::usage =
    "IGBipartiteQ[graph] tests if graph is bipartite.\n" <>
    "IGBipartiteQ[graph, {vertices1, vertices2}] verifies that no edges are running between the two given vertex subsets.";

SyntaxInformation[IGBipartiteQ] = {"ArgumentsPattern" -> {_, _.}};
IGBipartiteQ::bdprt = "`` are not two disjoint subsets of the graph vertices."
IGBipartiteQ[g_?igGraphQ] := Block[{ig = igMakeUnweighted[g]}, sck@ig@"bipartiteQ"[]]
IGBipartiteQ[g_?igGraphQ, {vertices1_List, vertices2_List}] :=
    With[{vertexList = VertexList[g]},
      If[Not[SubsetQ[vertexList, vertices1] && SubsetQ[vertexList, vertices2] && Intersection[vertices1, vertices2] === {}],
        Message[IGBipartiteQ::bdprt, {vertices1, vertices2}]
      ];
      EmptyGraphQ@igSubgraph[g, vertices1] && EmptyGraphQ@igSubgraph[g, vertices2]
    ]
IGBipartiteQ[_] := False (* provide this only for the single-argument form *)


PackageExport["IGBipartitePartitions"]
IGBipartitePartitions::usage =
     "IGBipartitePartitions[graph] partitions the vertices of a bipartite graph.\n" <>
     "IGBipartitePartitions[graph, vertex] ensures that the first partition which is returned contains vertex.";

SyntaxInformation[IGBipartitePartitions] = {"ArgumentsPattern" -> {_, _.}};
IGBipartitePartitions::nbipart = "The graph is not bipartite.";
IGBipartitePartitions[graph_?igGraphQ] :=
    catch@Block[{ig = igMakeUnweighted[graph], parts},
      parts = ig@"bipartitePartitions"[];
      If[MatchQ[parts, _LibraryFunctionError],
        Message[IGBipartitePartitions::nbipart];
        throw[$Failed]
      ];
      {Pick[VertexList[graph], parts, 0], Pick[VertexList[graph], parts, 1]}
    ]
IGBipartitePartitions[graph_?igGraphQ, vertex_] :=
    catch@Block[{ig = igMakeUnweighted[graph], parts, ind},
      parts = ig@"bipartitePartitions"[];
      If[MatchQ[parts, _LibraryFunctionError],
        Message[IGBipartitePartitions::nbipart];
        throw[$Failed]
      ];
      ind = vs1[graph][vertex];
      {Pick[VertexList[graph], parts, parts[[ind]] ], Pick[VertexList[graph], parts, 1 - parts[[ind]] ]}
    ]


PackageExport["IGBipartiteProjections"]
IGBipartiteProjections::usage =
    "IGBipartiteProjections[graph] gives both bipartite projections of graph. Multiplicities are returned as edge weights. Edge directions are ignored.\n" <>
    "IGBipartiteProjections[graph, {vertices1, vertices2}] returns both bipartite projections according to the specified partitioning.";

IGBipartiteProjections::bdpart = "`1` is not a valid partitioning of the vertices `2`.";
SyntaxInformation[IGBipartiteProjections] = {"ArgumentsPattern" -> {_, _.}};
IGBipartiteProjections[graph_?igGraphQ, parts : {vertices1_List, vertices2_List}] :=
    catch@Module[{ig = igMakeUnweighted[graph], ig1 = igMakeEmpty[], ig2 = igMakeEmpty[], weights},
      If[Not[Sort[Join@@parts] === Sort@VertexList[graph]],
        Message[IGBipartiteProjections::bdpart, parts, VertexList[graph]];
        throw[$Failed]
      ];
      weights = check@ig@"bipartiteProjection"[communitiesToMembership[VertexList[graph], parts], ManagedLibraryExpressionID[ig1], ManagedLibraryExpressionID[ig2]];
      With[{g1 = igToGraphWithNames[ig1, vertices1], g2 = igToGraphWithNames[ig2, vertices2]},
        {Graph[g1, EdgeWeight -> Take[weights,  EdgeCount[g1]]],
         Graph[g2, EdgeWeight -> Take[weights, -EdgeCount[g2]]]}
      ]
    ]

IGBipartiteProjections[graph_?igGraphQ] := IGBipartiteProjections[graph, IGBipartitePartitions[graph]]


(***** Bipartite incidence matrices *****)

PackageExport["IGBipartiteIncidenceMatrix"]
IGBipartiteIncidenceMatrix::usage =
    "IGBipartiteIncidenceMatrix[graph] gives the incidence matrix of a bipartite graph.\n" <>
    "IGBipartiteIncidenceMatrix[graph, {vertices1, vertices2}] uses the provided vertex partitioning.";

undirectedAdjacencyMatrix[g_?DirectedGraphQ] := With[{am = AdjacencyMatrix[g]}, am + Transpose[am]]
undirectedAdjacencyMatrix[g_] := AdjacencyMatrix[g]

IGBipartiteIncidenceMatrix::nbipart = IGBipartitePartitions::nbipart;
IGBipartiteIncidenceMatrix::bdpart  = "`` is not a valid bipartite partitioning of the graph.";
IGBipartiteIncidenceMatrix::empty   = "One of the graph partitions is empty.";

SyntaxInformation[IGBipartiteIncidenceMatrix] = {"ArgumentsPattern" -> {_, _.}};

IGBipartiteIncidenceMatrix[graph_?igGraphQ] :=
    catch@Block[{ig = igMakeUnweighted[graph], parts, posIndex},
      parts = ig@"bipartitePartitions"[];
      If[MatchQ[parts, _LibraryFunctionError],
        Message[IGBipartiteIncidenceMatrix::nbipart];
        throw[$Failed]
      ];
      posIndex = PositionIndex[parts];
      If[Sort@Keys[posIndex] =!= {0,1},
        Message[IGBipartiteIncidenceMatrix::empty];
        throw[$Failed]
      ];
      undirectedAdjacencyMatrix[graph][[posIndex[0], posIndex[1]]]
    ]

IGBipartiteIncidenceMatrix[graph_?igGraphQ, Automatic] := IGBipartiteIncidenceMatrix[graph]

IGBipartiteIncidenceMatrix[graph_?igGraphQ, parts : {vertices1_List, vertices2_List}] :=
    catch@Module[{asc},
      If[Not@Check[TrueQ@IGBipartiteQ[graph, parts], False],
        Message[IGBipartiteIncidenceMatrix::bdpart, parts];
        throw[$Failed]
      ];
      If[vertices1 === {} || vertices2 === {},
        Message[IGBipartiteIncidenceMatrix::empty];
        throw[$Failed]
      ];
      asc = AssociationThread[VertexList[graph], Range@VertexCount[graph]];
      undirectedAdjacencyMatrix[graph][[ Lookup[asc, vertices1], Lookup[asc, vertices2] ]]
    ]


PackageExport["IGBipartiteIncidenceGraph"]
IGBipartiteIncidenceGraph::usage =
    "IGBipartiteIncidenceGraph[mat] creates a bipartite graph from the given incidence matrix.\n" <>
    "IGBipartiteIncidenceGraph[{vertices1, vertices2}, mat] uses vertices1 and vertices2 as the vertex names in the two partitions.";

IGBipartiteIncidenceGraph::inv  = "`1` is not a valid bipartite incidence matrix.";
IGBipartiteIncidenceGraph::bdsz = "The vertex name lists `1` have an incompatible size with the provided incidence matrix.";
IGBipartiteIncidenceGraph::bdvl = "The vertex name lists should be disjoint. The following names are present in common: `1`.";

Options[IGBipartiteIncidenceGraph] = { DirectedEdges -> False };
SyntaxInformation[IGBipartiteIncidenceGraph] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}, "OptionNames" -> optNames[IGBipartiteIncidenceGraph, Graph]};

IGBipartiteIncidenceGraph[names : {vertices1_List, vertices2_List}, bm_?MatrixQ, opt : OptionsPattern[{IGBipartiteIncidenceGraph, Graph}]] :=
    Module[{sbm = If[emptyArrayQ[bm], {}, SparseArray[bm]], good = True},
      If[Not@MatrixQ[sbm, Internal`NonNegativeIntegerQ],
        Message[IGBipartiteIncidenceGraph::inv, bm];
        good = False
      ];
      If[good && Dimensions[sbm] =!= Length /@ names,
        Message[IGBipartiteIncidenceGraph::bdsz, names];
        good = False;
      ];
      If[good && Not@DisjointQ[vertices1, vertices2],
        Message[IGBipartiteIncidenceGraph::bdvl, Intersection @@ names];
        good = False;
      ];
      AdjacencyGraph[
        Join[vertices1, vertices2],
        ArrayFlatten[{{0, sbm},{If[TrueQ@OptionValue[DirectedEdges], SparseArray[{}, Reverse@Dimensions[sbm]], Transpose[sbm]], 0}}],
        opt
      ] /; good
    ]

IGBipartiteIncidenceGraph[bm_?MatrixQ, opt : OptionsPattern[{IGBipartiteIncidenceGraph, Graph}]] :=
    Module[{sbm = If[emptyArrayQ[bm], {}, SparseArray[bm]], good = True},
      If[Not@MatrixQ[sbm, Internal`NonNegativeIntegerQ],
        Message[IGBipartiteIncidenceGraph::inv, bm];
        good = False
      ];
      AdjacencyGraph[
        ArrayFlatten[{{0, sbm},{If[TrueQ@OptionValue[DirectedEdges], SparseArray[{}, Reverse@Dimensions[sbm]], Transpose[sbm]], 0}}],
        opt
      ] /; good
    ]
