(* Mathematica Package *)
(* Created by Mathematica plugin for IntelliJ IDEA *)

(* :Author: szhorvat *)
(* :Date: 2018-10-23 *)
(* :Copyright: (c) 2018-2020 Szabolcs Horvát *)

Package["IGraphM`"]


PackageImport["IGraphM`LTemplate`"] (* we use Make[] in the definition of lgMake[] *)

(******************************************************)
(***** Planar graphs and combinatorial embeddings *****)
(******************************************************)


(* Make an embedding *)

igIndexEmbedding[emb_?AssociationQ] :=
    With[{trans = AssociationThread[Keys[emb], Range@Length[emb]]},
      Lookup[trans, #]& /@ Values[emb]
    ]

IGraphM::invemb = "`1` is not a valid combinatorial embedding.";
embMake[emb_?AssociationQ] /; VectorQ[Values[emb], ListQ] && SubsetQ[Keys[emb], Catenate[emb]] :=
    Block[{embedding = Make["IGEmbedding"]},
      check@embedding@"set"[igIndexEmbedding[emb]-1];
      If[TrueQ@embedding@"validQ"[],
        embedding,
        Message[IGraphM::invemb, emb]; throw[$Failed]
      ]
    ]
embMake[emb_] := (Message[IGraphM::invemb, emb]; throw[$Failed])


(* Create a graph from an embedding. *)
igFromEmbedding[emb_] :=
    Graph[
      Keys[emb],
      DeleteDuplicates[Sort /@ Catenate[Function[{v, vlist}, UndirectedEdge[v, #] & /@ vlist] @@@ Normal[emb]]]
    ]



PackageExport["IGEmbeddingQ"]
IGEmbeddingQ::usage = "IGEmbeddingQ[embedding] tests if embedding represents a combinatorial embedding of a simple graph.";

SyntaxInformation[IGEmbeddingQ] = {"ArgumentsPattern" -> {_}};
IGEmbeddingQ[emb_Association?AssociationQ] :=
    Quiet[Not@MatchQ[catch@embMake[emb], _LibraryFunctionError|$Failed], IGraphM::invemb]
IGEmbeddingQ[_] := False


PackageExport["IGPlanarQ"]
IGPlanarQ::usage =
    "IGPlanarQ[graph] tests if graph is planar.\n" <>
    "IGPlanarQ[embedding] tests if a combinatorial embedding is planar.";

SyntaxInformation[IGPlanarQ] = {"ArgumentsPattern" -> {_}};
IGPlanarQ[graph_?EmptyGraphQ] := True
IGPlanarQ[graph_?igGraphQ] :=
    catch@Block[{lg = lgMake@UndirectedGraph@SimpleGraph[graph]},
      check@lg@"planarQ"[]
    ]
(* We do not use IGEmbeddingQ in the function pattern so that we can issue more useful error messages.
 * An error should be shown only if the user was likely to want to pass an embedding (which is an association),
 * but it was found not to be valid. *)
IGPlanarQ[embedding_?AssociationQ] :=
    TrueQ@catch@Block[{emb = embMake[embedding]},
      check@emb@"planarQ"[Length@ConnectedComponents@igFromEmbedding[embedding]]
    ]
IGPlanarQ[_] := False


PackageExport["IGMaximalPlanarQ"]
IGMaximalPlanarQ::usage = "IGMaximalPlanarQ[graph] tests if graph is maximal planar.";

SyntaxInformation[IGMaximalPlanarQ] = {"ArgumentsPattern" -> {_}};
IGMaximalPlanarQ[graph_?igGraphQ] /; IGPlanarQ[graph] :=
    Switch[VertexCount[graph],
      0, True, (* null graph *)
      1, True, (* singleton graph *)
      2, EdgeCount@SimpleGraph@UndirectedGraph[graph] == 1, (* below formula doesn't apply for K_2, which is maximal planar *)
      _, 3 VertexCount[graph] - 6 == EdgeCount@SimpleGraph@UndirectedGraph[graph]
    ]
IGMaximalPlanarQ[_] := False


PackageExport["IGOuterplanarQ"]
IGOuterplanarQ::usage =
    "IGOuterplanarQ[graph] tests if graph is outerplanar.\n" <>
    "IGOuterplanarQ[embedding] tests if a combinatorial embedding is outerplanar.";

outerplanarVertex::usage = "outerplanarVertex is used in the implementation of IGOuterplanarQ";

igOuterplanarAugment[graph_] :=
    EdgeAdd[
      VertexAdd[UndirectedGraph[graph], outerplanarVertex],
      UndirectedEdge[outerplanarVertex, #]& /@ VertexList[graph]
    ];

(* igConnectedOuterplanarQ expects embedding to be a planar embedding of a connected graph.
 * A planar embedding of a connected graphs is outerplanar if the largest faces contain all vertices,
 * or if it has no faces at all. *)
igConnectedOuterplanarQ[embedding_] :=
    With[{faces = IGFaces[embedding]},
      faces === {} || Sort@Keys[embedding] === Union@Extract[faces, Ordering[Length /@ faces, -1]]
    ]

SyntaxInformation[IGOuterplanarQ] = {"ArgumentsPattern" -> {_}};
IGOuterplanarQ[graph_?EmptyGraphQ] := True
IGOuterplanarQ[graph_?igGraphQ] := IGPlanarQ@igOuterplanarAugment[graph]
IGOuterplanarQ[embedding_?AssociationQ] :=
    TrueQ@catch@Module[{emb = embMake[embedding], comps},
      comps = ConnectedComponents@igFromEmbedding[embedding];
      check@emb@"planarQ"[Length[comps]] &&
           AllTrue[comps, igConnectedOuterplanarQ@KeyTake[embedding, #]&]
    ]
IGOuterplanarQ[_] := False



PackageExport["IGOuterplanarEmbedding"]
IGOuterplanarEmbedding::usage = "IGOuterplanarEmbedding[graph] gives an outerplanar combinatorial embedding of a graph.";

IGOuterplanarEmbedding::nopl = "The graph is not outer-planar.";
SyntaxInformation[IGOuterplanarEmbedding] = {"ArgumentsPattern" -> {_}};
IGOuterplanarEmbedding[graph_?igGraphQ] :=
    Module[{augmentedGraph = igOuterplanarAugment[graph], emb},
      emb = Quiet@IGPlanarEmbedding[augmentedGraph];
      If[Not@AssociationQ[emb],
        Message[IGOuterplanarEmbedding::nopl];
        Return[$Failed]
      ];
      DeleteCases[outerplanarVertex] /@ Delete[emb, Key[outerplanarVertex]]
    ]


PackageExport["IGKuratowskiEdges"]
IGKuratowskiEdges::usage = "IGKuratowskiEdges[graph] gives the edges belonging to a Kuratowski subgraph.";

SyntaxInformation[IGKuratowskiEdges] = {"ArgumentsPattern" -> {_}};
IGKuratowskiEdges[graph_?EmptyGraphQ] := {}
IGKuratowskiEdges[graph_?igGraphQ] :=
    catch@Block[{lg = lgMake[graph]},
      EdgeList[graph][[ igIndexVec@check@lg@"kuratowskiSubgraph"[] ]]
    ]


PackageExport["IGPlanarEmbedding"]
IGPlanarEmbedding::usage = "IGPlanarEmbedding[graph] gives a planar combinatorial embedding of a graph.";

(* TODO: The current implementation ignores edge multiplicities.
 * Investigate whether this can be improved. *)
SyntaxInformation[IGPlanarEmbedding] = {"ArgumentsPattern" -> {_, _.}};
IGPlanarEmbedding[graph_?igGraphQ] :=
    catch@Block[{lg = lgMake@SimpleGraph@UndirectedGraph[graph]},
      AssociationThread[VertexList[graph], igUnpackVertexSet[graph]@check@lg@"planarEmbedding"[]]
    ]


PackageExport["IGEmbeddingToCoordinates"]
IGEmbeddingToCoordinates::usage = "IGEmbeddingToCoordinates[embedding] gives the coordinates of a planar drawing based on the given combinatorial embedding.";

SyntaxInformation[IGEmbeddingToCoordinates] = {"ArgumentsPattern" -> {_}};
IGEmbeddingToCoordinates[emb_] :=
    catch@Block[{lg = lgMake@igFromEmbedding[emb], embedding = embMake[emb]},
      check@lg@"embeddingToCoordinates"[ManagedLibraryExpressionID[embedding]]
    ]


PackageExport["IGCoordinatesToEmbedding"]
IGCoordinatesToEmbedding::usage =
    "IGCoordinatesToEmbedding[graph] gives a combinatorial embedding based on the vertex coordinates of graph.\n" <>
    "IGCoordinatesToEmbedding[graph, coordinates] uses the given coordinates instead of the VertexCoordinates property.";

SyntaxInformation[IGCoordinatesToEmbedding] = {"ArgumentsPattern" -> {_, _.}};
IGCoordinatesToEmbedding[graph_?igGraphQ, coords_] :=
    catch@Block[{ig = igMakeFast[graph]},
      AssociationThread[VertexList[graph], igUnpackVertexSet[graph]@check@ig@"coordinatesToEmbedding"[coords]]
    ]
IGCoordinatesToEmbedding[graph_?igGraphQ] := IGCoordinatesToEmbedding[graph, GraphEmbedding[graph]]


PackageExport["IGFaces"]
IGFaces::usage =
    "IGFaces[graph] gives the faces of a planar graph.\n" <>
    "IGFaces[embedding] gives the faces that correspond to a combinatorial embedding.";

SyntaxInformation[IGFaces] = {"ArgumentsPattern" -> {_, _.}};
IGFaces[embedding_?AssociationQ] :=
    catch@Block[{emb = embMake[embedding]},
      igUnpackVertexSet[embedding]@check@emb@"faces"[]
    ]
IGFaces[graph_?igGraphQ] := catch@IGFaces@check@IGPlanarEmbedding[graph]


PackageExport["IGDualGraph"]
IGDualGraph::usage =
    "IGDualGraph[graph] gives the dual graph of a planar graph.\n" <>
    "IGDualGraph[embedding] gives the dual graph corresponding to a specific embedding of a graph. The embedding does not need to be planar.";

SyntaxInformation[IGDualGraph] = {"ArgumentsPattern" -> {_, OptionsPattern[]}, "OptionNames" -> optNames[Graph]};
IGDualGraph[embedding_?AssociationQ, opt : OptionsPattern[]] :=
    catch@Block[{emb = embMake[embedding]},
      With[{edges = 1 + check@emb@"dualGraph"[]},
        SimpleGraph@Graph[Range@Max[0, edges], edges, DirectedEdges -> False, opt]
      ]
    ]
IGDualGraph[graph_?igGraphQ, opt : OptionsPattern[]] := catch@IGDualGraph[check@IGPlanarEmbedding[graph], opt]
