(* Mathematica Package *)
(* Created by Mathematica plugin for IntelliJ IDEA *)

(* :Author: szhorvat *)
(* :Date: 2018-10-24 *)
(* :Copyright: (c) 2019 Szabolcs Horvát *)

Package["IGraphM`"]
igContextSetup[igPackagePrivateSymbol]

(**********************************************)
(***** Utilities for edge-weighted graphs *****)
(**********************************************)


PackageExport["IGUnweighted"]
IGUnweighted::usage = "IGUnweighted[graph] returns an unweighted version of an edge-weighted graph, while preserving other graph properties.";
SyntaxInformation[IGUnweighted] = {"ArgumentsPattern" -> {_}};
IGUnweighted[g_?IGEdgeWeightedQ] := transformGraphOptions[ FilterRules[#, Except[EdgeWeight]]& ][g]
IGUnweighted[g_?GraphQ] := g
(* A more direct implementation is

     IGUnweighted[g_?IGEdgeWeightedQ] := Graph[ VertexList[g], EdgeList[g], FilterRules[Options[g], Except[EdgeWeight]] ]

   This implementation is slightly faster if the graph has a large number of properties, such as
   ExampleData[{"NetworkGraph", "CondensedMatterCollaborations"}]. However, it is noticeable slower for weighted graphs with
   packed weight vectors and no other properties.
*)


(* TODO: consider allowing all Graph options *)
PackageExport["IGDistanceWeighted"]
IGDistanceWeighted::usage = "IGDistanceWeighted[graph] sets the weight of each edge to be the geometrical distance between its endpoints.";
Options[IGDistanceWeighted] = { DistanceFunction -> EuclideanDistance };
SyntaxInformation[IGDistanceWeighted] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};
IGDistanceWeighted[graph_?GraphQ, opt : OptionsPattern[]] :=
    With[{distanceFunction = OptionValue[DistanceFunction]},
      Switch[distanceFunction,
        EuclideanDistance, igDistanceWeightedVectorized[graph, Sqrt@Total[#^2]&],
        SquaredEuclideanDistance, igDistanceWeightedVectorized[graph, Total[#^2]&],
        ManhattanDistance, igDistanceWeightedVectorized[graph, Total@*Abs],
        _, igDistanceWeighted[graph, distanceFunction]
      ]
    ]

igDistanceWeighted[graph_, distanceFunction_] :=
    IGEdgeMap[Apply[distanceFunction], EdgeWeight -> IGEdgeVertexProp[VertexCoordinates], graph]

igDistanceWeightedVectorized[graph_, fun_] :=
    igSetEdgeProperty[
      graph,
      EdgeWeight,
      With[{pts = GraphEmbedding[graph]},
        fun@Transpose[Subtract @@ (pts[[#]] & /@ Transpose@IGIndexEdgeList[graph])]
      ]
    ]


(***** Weighted adjacency matrix handling *****)

PackageExport["IGWeightedAdjacencyGraph"]
IGWeightedAdjacencyGraph::usage =
    "IGWeightedAdjacencyGraph[matrix] creates a graph from a weighted adjacency matrix, taking 0 to mean unconnected.\n" <>
    "IGWeightedAdjacencyGraph[vertices, matrix] uses vertices as the vertex names.\n" <>
    "IGWeightedAdjacencyGraph[matrix, z] creates a graph from a weighted adjacency matrix, taking the value z to mean unconnected.\n" <>
    "IGWeightedAdjacencyGraph[vertices, matrix, z] uses vertices as the vertex names.";

(*
(* When am is a SparseArray, we need to canonicalize it and ensure that it has no explicit value that is the same as the implicit one. *)
arrayRules[am_SparseArray, u_] := ArrayRules[SparseArray[am], u]
arrayRules[am_, u_] := ArrayRules[am, u]

SyntaxInformation[IGWeightedAdjacencyGraph] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}, "OptionNames" -> optNames[WeightedAdjacencyGraph]};
IGWeightedAdjacencyGraph[wam_?SquareMatrixQ, unconnected : Except[_?OptionQ] : 0, opt : OptionsPattern[WeightedAdjacencyGraph]] :=
      WeightedAdjacencyGraph[
        SparseArray[Most@arrayRules[wam, unconnected], Dimensions[wam], Infinity],
        opt
      ]
IGWeightedAdjacencyGraph[vertices_List, wam_?SquareMatrixQ, unconnected : Except[_?OptionQ] : 0, opt : OptionsPattern[WeightedAdjacencyGraph]] :=
    WeightedAdjacencyGraph[
      vertices,
      SparseArray[Most@arrayRules[wam, unconnected], Dimensions[wam], Infinity],
      opt
    ]
*)

SyntaxInformation[IGWeightedAdjacencyGraph] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}, "OptionNames" -> optNames[Graph]};
IGWeightedAdjacencyGraph[wam_?SquareMatrixQ, unconnected : Except[_?OptionQ] : 0, opt : OptionsPattern[Graph]] :=
    IGWeightedAdjacencyGraph[Range@Length[wam], wam, unconnected, opt]
IGWeightedAdjacencyGraph[vertices_List, wam_?SquareMatrixQ, unconnected : Except[_?OptionQ] : 0, opt : OptionsPattern[Graph]] :=
    Module[{sa = SparseArray[wam, Automatic, unconnected], directedEdges = OptionValue[DirectedEdges]},
      If[Length[vertices] != Length[sa],
        Message[IGWeightedAdjacencyGraph::ndims, vertices, wam];
        Return[$Failed]
      ];
      If[directedEdges === Automatic,
        directedEdges = Not@SymmetricMatrixQ[sa]
      ];
      If[Not[directedEdges],
        sa = UpperTriangularize[sa]
      ];
      Graph[vertices, sa["NonzeroPositions"], EdgeWeight -> sa["NonzeroValues"], DirectedEdges -> directedEdges, opt]
    ]


PackageExport["IGWeightedAdjacencyMatrix"]
IGWeightedAdjacencyMatrix::usage =
    "IGWeightedAdjacencyMatrix[graph] gives the adjacency matrix of the edge weights of graph.\n" <>
        "IGWeightedAdjacencyMatrix[graph, z] gives the adjacency matrix of the edge weights of graph, using the value z to represent absent connections.";

Options[IGWeightedAdjacencyMatrix] = Options[WeightedAdjacencyMatrix];
SyntaxInformation[IGWeightedAdjacencyMatrix] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};
IGWeightedAdjacencyMatrix[graph_?IGNullGraphQ, unconnected : Except[_?OptionQ] : 0, opt : OptionsPattern[]] :=
    (Message[IGWeightedAdjacencyMatrix::nadj]; $Failed)  (* This falls back on General::nadj *)
IGWeightedAdjacencyMatrix[graph_?GraphQ, unconnected : Except[_?OptionQ] : 0, opt : OptionsPattern[]] :=
    With[{sa = WeightedAdjacencyMatrix[graph, opt]},
      SparseArray[sa["NonzeroPositions"] -> sa["NonzeroValues"], Dimensions[sa], unconnected]
    ]


(***** Test for edge and vertex weightedness separately *****)

(*
 * The original implementation was
 *     WeightedGraphQ[#] && PropertyValue[#, EdgeWeight] =!= Automatic &;
 * for testing edge-weightedness.
 * PropertyValue[g, EdgeWeight] fails on the null graph. This is why we had to test with WeightedGraphQ first.
 *
 * This solution was slow because PropertyValue takes a very long time.
 * (Benchmarked on ExampleData[{"NetworkGraph", "CondensedMatterCollaborations2005"}].)
 *
 * The alternative implementation
 *     MemberQ[PropertyList[graph], EdgeWeight]
 * does not return correct result for vertex-weighted but non-edge-weighted graphs
 *
 * Update: Version 12.0 introduces EdgeWeightedGraphQ and VertexWeightedGraphQ.
 *)

If[$VersionNumber >= 10.1, (* MinMax was added in M10.1 *)
  minMax = MinMax,
  minMax = {Min[#], Max[#]}&
];


PackageExport["IGVertexWeightedQ"]
IGVertexWeightedQ::usage = "IGVertexWeightedQ[graph] tests if graph is a vertex-weighted graph.";

SyntaxInformation[IGVertexWeightedQ] = {"ArgumentsPattern" -> {_}};
If[$VersionNumber >= 12.0,
  IGVertexWeightedQ[g_] := VertexWeightedGraphQ[g] (* new in 12.0 *)
  ,
  IGVertexWeightedQ[g_] :=
      WeightedGraphQ[g] &&
          With[{weights = Developer`ToPackedArray@GraphComputation`WeightVector[g]},
            If[First[weights] === 1 && minMax[weights] === {1, 1},
              PropertyValue[g, VertexWeight] =!= Automatic,
              True
            ]
          ]
]


PackageExport["IGEdgeWeightedQ"]
IGEdgeWeightedQ::usage = "IGEdgeWeightedQ[graph] tests if graph is an edge-weighted graph.";

SyntaxInformation[IGEdgeWeightedQ] = {"ArgumentsPattern" -> {_}};
If[$VersionNumber >= 12.0,
  IGEdgeWeightedQ[g_] := EdgeWeightedGraphQ[g] (* new in 12.0 *)
  ,
  IGEdgeWeightedQ[g_?EmptyGraphQ] := False; (* avoid error from First if graph has no edges but is vertex weighted *)
  IGEdgeWeightedQ[g_] :=
      WeightedGraphQ[g] &&
          With[{weights = Developer`ToPackedArray@GraphComputation`WeightValues[g]},
            If[First[weights] === 1 && minMax[weights] === {1, 1},
              PropertyValue[g, EdgeWeight] =!= Automatic,
              True
            ]
          ]
]


(***** Weighted graph modification *****)

(* Notes:
    - igWeightedAdjacencyGraph requires a sparse matrix input. A dense matrix cannot be used.
    - igWeightedAdjacencyGraph must ignore the background value of the sparse matrix, as this may not be 0.
      See IGWeightedSimpleGraph for an example use.
 *)

igWeightedAdjacencyGraph::dummy = "igWeightedAdjacencyGraph[vertexList, sparseAdjacencyMatrix, directedQ, loopQ, graphOptions]";

igWeightedAdjacencyGraph[vertices_, sa_SparseArray, True, True, opt___] :=
    Graph[vertices, sa["NonzeroPositions"], DirectedEdges -> True, EdgeWeight -> sa["NonzeroValues"], opt]

igWeightedAdjacencyGraph[vertices_, sa_SparseArray, True, False, opt___] :=
    With[{pos = igraphGlobal@"nondiagPos"[sa["NonzeroPositions"]]},
      Graph[vertices, sa["NonzeroPositions"][[pos]], DirectedEdges -> True, EdgeWeight -> sa["NonzeroValues"][[pos]], opt]
    ]

igWeightedAdjacencyGraph[vertices_, sa_SparseArray, False, loop_, opt___] :=
    With[{pos = igraphGlobal@"upperPos"[sa["NonzeroPositions"], loop]},
      Graph[vertices, sa["NonzeroPositions"][[pos]], DirectedEdges -> False, EdgeWeight -> sa["NonzeroValues"][[pos]], opt]
    ]

$unconnected::dummy = "$unconnected is used to denote unconnected entries in the implementation of IGWeightedSimpleGraph and IGWeightedUndirectedGraph.";

PackageExport["IGWeightedSimpleGraph"]
IGWeightedSimpleGraph::usage =
    "IGWeightedSimpleGraph[graph] combines parallel edges by adding their weights. If graph is not weighted, the resulting weights will be the edge multiplicities of graph.\n" <>
        "IGWeightedSimpleGraph[graph, comb] applies the function comb to the weights of parallel edges to compute a new weight. The default combiner is Plus.";

Options[IGWeightedSimpleGraph] = { SelfLoops -> True };
SyntaxInformation[IGWeightedSimpleGraph] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}, "OptionNames" -> optNames[IGWeightedSimpleGraph, Graph]};
IGWeightedSimpleGraph[g_?EmptyGraphQ, comb : Except[_?OptionQ] : Total, opt : OptionsPattern[{IGWeightedSimpleGraph, Graph}]] :=
    applyGraphOpt[opt][g]
IGWeightedSimpleGraph[g_?igGraphQ, comb : Except[_?OptionQ] : Total, opt : OptionsPattern[{IGWeightedSimpleGraph, Graph}]] :=
    With[{sao = SystemOptions["SparseArrayOptions"], vc = VertexCount[g]},
      Internal`WithLocalSettings[
        SetSystemOptions["SparseArrayOptions" -> "TreatRepeatedEntries" -> comb],
        igWeightedAdjacencyGraph[
          VertexList[g],
          SparseArray[IGIndexEdgeList[g] -> igEdgeWeights[g], {vc, vc}, $unconnected],
          DirectedGraphQ[g],
          TrueQ@OptionValue[SelfLoops],
          FilterRules[{opt}, Options[Graph]]
        ],
        SetSystemOptions[sao]
      ]
    ]


PackageExport["IGWeightedUndirectedGraph"]
IGWeightedUndirectedGraph::usage =
    "IGWeightedUndirectedGraph[graph] converts an edge-weighted directed graph to an undirected one. The weights of reciprocal edges added up.\n" <>
        "IGWeightedUndirectedGraph[graph, comb] applies the function comb to the weights of reciprocal edges to compute the weight of the corresponding undirected edge.\n" <>
        "IGWeightedUndirectedGraph[graph, None] converts each directed edge to an undirected one without combining their weights. The result may be a multigraph.";

IGWeightedUndirectedGraph::mg = "The input is a multigraph. Weights of parallel edges will be combined with the same combiner function as used for reciprocal edges.";
SyntaxInformation[IGWeightedUndirectedGraph] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}, "OptionNames" -> optNames[Graph]};
IGWeightedUndirectedGraph[g_?UndirectedGraphQ, comb : Except[_?OptionQ] : Total, opt : OptionsPattern[Graph]] := (* also catches empty case *)
    applyGraphOpt[opt][g]
IGWeightedUndirectedGraph[g_?igGraphQ, None, opt : OptionsPattern[Graph]] :=
    If[IGEdgeWeightedQ[g],
      (* weighted case *)
      Graph[VertexList[g], IGIndexEdgeList[g], DirectedEdges -> False, EdgeWeight -> igEdgeWeights[g]],
      (* unweighted case *)
      UndirectedGraph[g, opt]
    ]
IGWeightedUndirectedGraph[g_?igGraphQ, comb : Except[_?OptionQ] : Total, opt : OptionsPattern[Graph]] :=
    If[IGEdgeWeightedQ[g],
      (* weighted case *)
      If[MultigraphQ[g], Message[IGWeightedUndirectedGraph::mg]];
      With[{sao = SystemOptions["SparseArrayOptions"], vc = VertexCount[g]},
        Internal`WithLocalSettings[
          SetSystemOptions["SparseArrayOptions" -> "TreatRepeatedEntries" -> comb],
          igWeightedAdjacencyGraph[
            VertexList[g],
            SparseArray[(igraphGlobal@"edgeListSortPairs"@IGIndexEdgeList[g]) -> igEdgeWeights[g], {vc, vc}, $unconnected],
            False,
            True,
            opt
          ],
          SetSystemOptions[sao]
        ]
      ]
      ,
      (* unweighted case *)
      UndirectedGraph[g, opt]
    ]


PackageExport["IGWeightedVertexDelete"]
IGWeightedVertexDelete::usage =
    "IGWeightedVertexDelete[graph, vertex] deletes the given vertex while preserving edge weights.\n" <>
        "IGWeightedVertexDelete[graph, {v1, v2, \[Ellipsis]}] deletes the given set of vertices while preserving edge weights.";

(* Delete with multiple indices is slow on non-packed arrays: https://mathematica.stackexchange.com/q/187206/12 *)
fastDelete[list_, inds_] := Part[list, Delete[Range@Length[list], Transpose@Developer`ToPackedArray@List@inds]]

SyntaxInformation[IGWeightedVertexDelete] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}, "OptionNames" -> optNames[Graph]};
IGWeightedVertexDelete[g_?igGraphQ, vs_List, opt : OptionsPattern[Graph]] :=
    catch@Module[{elist, emarker, vinds},
      Check[
        vinds = VertexIndex[g, #]& /@ DeleteDuplicates[vs],
        throw[$Failed]
      ];
      elist = IGIndexEdgeList[g];
      emarker = igraphGlobal@"edgeListMarkWhenEitherPresent"[elist, vinds];
      Graph[
        fastDelete[VertexList[g], vinds],
        igraphGlobal@"edgeListReindexAfterDelete"[Pick[elist, emarker, 0], vinds],
        If[IGEdgeWeightedQ[g], EdgeWeight -> Pick[igEdgeWeights[g], emarker, 0], {}],
        DirectedEdges -> DirectedGraphQ[g],
        opt
      ]
    ]
IGWeightedVertexDelete[g_?igGraphQ, vertex_, opt : OptionsPattern[Graph]] := IGWeightedVertexDelete[g, {vertex}, opt]


PackageExport["IGWeightedSubgraph"]
IGWeightedSubgraph::usage = "IGWeightedSubgraph[graph, {v1, v2, \[Ellipsis]}] returns the subgraph induced by the given vertices while preserving edge weights.";

(* TODO support edges in subgraph spec *)
SyntaxInformation[IGWeightedSubgraph] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}, "OptionNames" -> optNames[Graph]};
IGWeightedSubgraph[g_?igGraphQ, vs_List, opt : OptionsPattern[Graph]] :=
    catch@Module[{vinds, elist, emarker},
      Check[
        vinds = VertexIndex[g, #]& /@ DeleteDuplicates[vs],
        throw[$Failed]
      ];
      elist = IGIndexEdgeList[g];
      emarker = igraphGlobal@"edgeListMarkWhenBothPresent"[elist, vinds];
      Graph[
        VertexList[g][[vinds]],
        igraphGlobal@"edgeListReindex"[Pick[elist, emarker, 1], vinds],
        If[IGEdgeWeightedQ[g], EdgeWeight -> Pick[igEdgeWeights[g], emarker, 1], {}],
        DirectedEdges -> DirectedGraphQ[g],
        opt
      ]
    ]


(***** Vertex strengths, i.e. sums of the weights of incident edges *****)

PackageExport["IGVertexStrength"]
IGVertexStrength::usage =
    "IGVertexStrength[graph] returns the sum of edge weights for edges connecting to each vertex in graph.\n" <>
        "IGVertexStrength[graph, v] returns the sum of edge weights for edges connecting to vertex v in graph.";

SyntaxInformation[IGVertexStrength] = {"ArgumentsPattern" -> {_, _.}};
IGVertexStrength[g_?IGNullGraphQ] := {}
IGVertexStrength[g_?igGraphQ] :=
    With[{am = WeightedAdjacencyMatrix[g]}, (* WeightedAdjacencyMatrix adds up weights in multigraphs. *)
      If[DirectedGraphQ[g],
        Total[am] + Total[am, {2}],
        Total[am] + Diagonal[am]
      ]
    ]
IGVertexStrength[g_?igGraphQ, v_] :=
    With[{index = VertexIndex[g, v]},
      With[{am = WeightedAdjacencyMatrix[g]},
        If[DirectedGraphQ[g],
          Total[am[[index]]] + Total[am[[;;, index]]],
          Total[am[[index]]] + am[[index, index]]
        ]
      ] /; IntegerQ[index]
    ]


PackageExport["IGVertexInStrength"]
IGVertexInStrength::usage =
    "IGVertexInStrength[graph] returns the sum of edge weights for the incoming edges of each vertex in graph.\n" <>
        "IGVertexInStrength[graph, v] returns the sum of edge weights for incoming edges of vertex v in graph.";

SyntaxInformation[IGVertexInStrength] = {"ArgumentsPattern" -> {_, _.}};
IGVertexInStrength[g_?IGNullGraphQ] := {}
IGVertexInStrength[g_?igGraphQ] :=
    With[{am = WeightedAdjacencyMatrix[g]},
      If[DirectedGraphQ[g],
        Total[am],
        Total[am] + Diagonal[am]
      ]
    ]
IGVertexInStrength[g_?igGraphQ, v_] :=
    With[{index = VertexIndex[g, v]},
      With[{am = WeightedAdjacencyMatrix[g]},
        If[DirectedGraphQ[g],
          Total[am[[;;, index]]],
          Total[am[[index]]] + am[[index, index]]
        ]
      ] /; IntegerQ[index]
    ]


PackageExport["IGVertexOutStrength"]
IGVertexOutStrength::usage =
    "IGVertexOutStrength[graph] returns the sum of edge weights for the outgoing edges of each vertex in graph.\n" <>
        "IGVertexOutStrength[graph, v] returns the sum of edge weights for outgoing edges of vertex v in graph.";

SyntaxInformation[IGVertexOutStrength] = {"ArgumentsPattern" -> {_, _.}};
IGVertexOutStrength[g_?IGNullGraphQ] := {}
IGVertexOutStrength[g_?igGraphQ] :=
    With[{am = WeightedAdjacencyMatrix[g]},
      If[DirectedGraphQ[g],
        Total[am, {2}],
        Total[am] + Diagonal[am]
      ]
    ]
IGVertexOutStrength[g_?igGraphQ, v_] :=
    With[{index = VertexIndex[g, v]},
      With[{am = WeightedAdjacencyMatrix[g]},
        If[DirectedGraphQ[g],
          Total[am[[index]]],
          Total[am[[index]]] + am[[index, index]]
        ]
      ] /; IntegerQ[index]
    ]
