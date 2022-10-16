(* Mathematica Package *)
(* Created by Mathematica plugin for IntelliJ IDEA *)

(* :Author: szhorvat *)
(* :Date: 2018-10-24 *)
(* :Copyright: (c) 2018-2022 Szabolcs Horvát *)

Package["IGraphM`"]


(*****************************************)
(***** Connectivity and maximum flow *****)
(*****************************************)


(***** Testing *****)

PackageExport["IGConnectedQ"]
IGConnectedQ::usage = "IGConnectedQ[graph] tests if graph is strongly connected.";

SyntaxInformation[IGConnectedQ] = {"ArgumentsPattern" -> {_}};
igConnectedQ[g_] := Block[{ig = igMakeUnweighted[g]}, check@ig@"connectedQ"[True]]
IGConnectedQ[g_?igGraphQ] := catch@cachedFun[igConnectedQ][g]
IGConnectedQ[_] := False


PackageExport["IGWeaklyConnectedQ"]
IGWeaklyConnectedQ::usage = "IGWeaklyConnectedQ[graph] tests if graph is weakly connected.";

SyntaxInformation[IGWeaklyConnectedQ] = {"ArgumentsPattern" -> {_}};
igWeaklyConnectedQ[g_] := Block[{ig = igMakeUnweighted[g]}, check@ig@"connectedQ"[False]]
IGWeaklyConnectedQ[g_?igGraphQ] := catch@cachedFun[igWeaklyConnectedQ][g]
IGWeaklyConnectedQ[_] := False


PackageExport["IGBiconnectedQ"]
IGBiconnectedQ::usage = "IGBiconnectedQ[graph] tests if graph is biconnected.";

SyntaxInformation[IGBiconnectedQ] = {"ArgumentsPattern" -> {_}};
IGBiconnectedQ[graph_?igGraphQ] :=
    Block[{ig = igMakeUnweighted[graph]},
      sck@ig@"biconnectedQ"[]
    ]
IGBiconnectedQ[_] := False (* for non-graphs *)


(***** Vertex cuts *****)

PackageExport["IGMinSeparators"]
IGMinSeparators::usage = "IGMinSeparators[] is deprecated. Use IGMinimumSeparators[] instead."; (* deprecated in favour of IGMinimumSeparators *)
IGMinSeparators::deprec = "IGMinSeparators is deprecated and will be removed from future versions of IGraph/M. Use IGMinimumSeparators instead.";

IGMinSeparators[graph_] := (Message[IGMinSeparators::deprec]; IGMinimumSeparators[graph]) (* TODO: remove eventually *)


PackageExport["IGMinimumSeparators"]
IGMinimumSeparators::usage = "IGMinimumSeparators[graph] gives all separator vertex sets of minimum size. A vertex set is a separator if its removal disconnects the graph. Edge directions are ignored.";

SyntaxInformation[IGMinimumSeparators] = {"ArgumentsPattern" -> {_}};
IGMinimumSeparators[graph_?igGraphQ] :=
    catch@Block[{ig = igMakeUnweighted@UndirectedGraph[graph]},
      igUnpackVertexSet[graph]@check@ig@"minimumSizeSeparators"[]
    ]


PackageExport["IGMinimalSeparators"]
IGMinimalSeparators::usage = "IGMinimalSeparators[graph] gives all minimal separator vertex sets. A vertex set is a separator if its removal disconnects the graph. Edge directions are ignored."

SyntaxInformation[IGMinimalSeparators] = {"ArgumentsPattern" -> {_}};
IGMinimalSeparators[graph_?igGraphQ] :=
    catch@Block[{ig = igMakeUnweighted[graph]},
      igUnpackVertexSet[graph]@check@ig@"minimalSeparators"[]
    ]


PackageExport["IGVertexSeparatorQ"]
IGVertexSeparatorQ::usage = "IGVertexSeparatorQ[graph, {vertex1, vertex2, \[Ellipsis]}] tests if the given set of vertices disconnects the graph. Edge directions are ignored.";

SyntaxInformation[IGVertexSeparatorQ] = {"ArgumentsPattern" -> {_, _}};
IGVertexSeparatorQ[graph_?igGraphQ, {}] := False
IGVertexSeparatorQ[graph_?igGraphQ, vs_List] :=
    catch@Block[{ig = igMakeUnweighted[graph]},
      check@ig@"separatorQ"[vss[graph][vs]]
    ]


PackageExport["IGMinimalVertexSeparatorQ"]
IGMinimalVertexSeparatorQ::usage = "IGMinimalVertexSeparatorQ[graph, {vertex1, vertex2, \[Ellipsis]}] tests if the given vertex set is a minimal separator. Edge directions are ignored.";

SyntaxInformation[IGMinimalVertexSeparatorQ] = {"ArgumentsPattern" -> {_, _}};
IGMinimalVertexSeparatorQ[graph_?igGraphQ, {}] := False
IGMinimalVertexSeparatorQ[graph_?igGraphQ, vs_List] :=
    catch@Block[{ig = igMakeUnweighted[graph]},
      check@ig@"minSeparatorQ"[vss[graph][vs]]
    ]


(***** Connected components *****)

PackageExport["IGArticulationPoints"]
IGArticulationPoints::usage = "IGArticulationPoints[graph] gives the articulation points of graph. A vertex is an articulation point if its removal increases the number of (weakly) connected components in the graph.";

SyntaxInformation[IGArticulationPoints] = {"ArgumentsPattern" -> {_}};
IGArticulationPoints[graph_?igGraphQ] :=
    catch@Block[{ig = igMakeUnweighted[graph]},
      igVertexNames[graph]@igIndexVec@check@ig@"articulationPoints"[]
    ]


PackageExport["IGBiconnectedComponents"]
IGBiconnectedComponents::usage = "IGBiconnectedComponents[graph] gives the vertices of the maximal biconnected subgraphs of graph. A graph is biconnected if the removal of any single vertex does not disconnect it. Isolated vertices are not returned.";

SyntaxInformation[IGBiconnectedComponents] = {"ArgumentsPattern" -> {_}};
IGBiconnectedComponents[graph_?igGraphQ] :=
    catch@Block[{ig = igMakeUnweighted[graph]},
      igUnpackVertexSet[graph]@check@ig@"biconnectedComponents"[]
    ]

PackageExport["IGBiconnectedEdgeComponents"]
IGBiconnectedEdgeComponents::usage = "IGBiconnectedEdgeComponents[graph] gives the edges of the maximal biconnected subgraphs of graph. A graph is biconnected if the removal of any single vertex does not disconnect it.";

SyntaxInformation[IGBiconnectedEdgeComponents] = {"ArgumentsPattern" -> {_}};
IGBiconnectedEdgeComponents[graph_?igGraphQ] :=
    catch@Block[{ig = igMakeUnweighted[graph]},
      igUnpackEdgeSet[graph]@check@ig@"biconnectedEdgeComponents"[]
    ]


PackageExport["IGBridges"]
IGBridges::usage = "IGBridges[graph] gives the bridges of graph. A bridge is an edge whose removal increases the number of (weakly) connected components in the graph.";

SyntaxInformation[IGBridges] = {"ArgumentsPattern" -> {_}};
IGBridges[graph_?igGraphQ] :=
    catch@Block[{ig = igMake[graph]},
      EdgeList[graph][[ igIndexVec@check@ig@"bridges"[] ]]
    ]


PackageExport["IGConnectedComponentSizes"]
IGConnectedComponentSizes::usage = "IGConnectedComponentSizes[graph] gives the sizes of graph's connected components in decreasing order.";

SyntaxInformation[IGConnectedComponentSizes] = {"ArgumentsPattern" -> {_}};
IGConnectedComponentSizes[graph_?igGraphQ] :=
    catch@Block[{ig = igMakeUnweighted[graph]},
      Reverse@Sort@Round@check@ig@"connectedComponentSizes"[True]
    ]


PackageExport["IGWeaklyConnectedComponentSizes"]
IGWeaklyConnectedComponentSizes::usage = "IGWeaklyConnectedComponentSizes[graph] gives the sizes of graph's weakly connected components in decreasing order.";

SyntaxInformation[IGWeaklyConnectedComponentSizes] = {"ArgumentsPattern" -> {_}};
IGWeaklyConnectedComponentSizes[graph_?igGraphQ] :=
    catch@Block[{ig = igMakeUnweighted[graph]},
      Reverse@Sort@Round@check@ig@"connectedComponentSizes"[False]
    ]


(***** Vertex and edge connectivity *****)

PackageExport["IGVertexConnectivity"]
IGVertexConnectivity::usage =
    "IGVertexConnectivity[graph] gives the smallest number of vertices whose deletion disconnects graph.\n" <>
    "IGVertexConnectivity[graph, s, t] gives the smallest number of vertices whose deletion disconnects vertices s and t in graph.";

SyntaxInformation[IGVertexConnectivity] = {"ArgumentsPattern" -> {_, _., _.}};

IGVertexConnectivity[graph_?igGraphQ] :=
    Block[{ig = igMakeUnweighted[graph]},
      sck@ig@"vertexConnectivity"[]
    ]

IGVertexConnectivity[graph_?igGraphQ, s_, t_] :=
    catch@Block[{ig = igMakeUnweighted[graph]},
      check@ig@"vertexConnectivityST"[vs[graph][s], vs[graph][t]]
    ]


PackageExport["IGEdgeConnectivity"]
IGEdgeConnectivity::usage =
    "IGEdgeConnectivity[graph] gives the smallest number of edges whose deletion disconnects graph.\n" <>
    "IGEdgeConnectivity[graph, s, t] gives the smallest number of edges whose deletion disconnects vertices s and t in graph.";

SyntaxInformation[IGEdgeConnectivity] = {"ArgumentsPattern" -> {_, _., _.}};

IGEdgeConnectivity[graph_?igGraphQ] :=
    Block[{ig = igMakeUnweighted[graph]},
      sck@ig@"edgeConnectivity"[]
    ]

IGEdgeConnectivity[graph_?igGraphQ, s_, t_] :=
    catch@Block[{ig = igMakeUnweighted[graph]},
      check@ig@"edgeConnectivityST"[vs[graph][s], vs[graph][t]]
    ]


PackageExport["IGCohesiveBlocks"]
IGCohesiveBlocks::usage = "IGCohesiveBlocks[graph] gives the cohesive block structure of a simple undirected graph.";

IGCohesiveBlocks::badarg = "The input must be a simple undirected graph.";
SyntaxInformation[IGCohesiveBlocks] = {"ArgumentsPattern" -> {_}};
IGCohesiveBlocks[graph_?igGraphQ] :=
    catch@Block[{ig = igMakeUnweighted[graph], blocks, cohesion, parents},
      If[Not@SimpleGraphQ[graph], Message[IGCohesiveBlocks::badarg]; Return[$Failed]];
      {blocks, cohesion, parents} = check@ig@"cohesiveBlocks"[];
      {igVertexNames[graph] /@ igIndexVec[blocks], Round[cohesion](*, igIndexVec[parents]*)}
    ]


(***** Minimum cuts (for weighted graphs *****)

PackageExport["IGMinimumCutValue"]
IGMinimumCutValue::usage =
    "IGMinimumCutValue[graph] gives the smallest sum of weights corresponding to an edge cut in graph.\n" <>
    "IGMinimumCutValue[graph, s, t] gives the smallest sum of weights corresponding to an s-t edge cut in graph.";

SyntaxInformation[IGMinimumCutValue] = {"ArgumentsPattern" -> {_, _., _.}};
IGMinimumCutValue[graph_?igGraphQ] :=
    Block[{ig = igMake[graph]},
      sck@ig@"minCutValue"[]
    ]
IGMinimumCutValue[graph_?igGraphQ, s_, t_] :=
    catch@Block[{ig = igMake[graph]},
      check@ig@"minCutValueST"[vs[graph][s], vs[graph][t]]
    ]


PackageExport["IGMinimumCut"]
IGMinimumCut::usage =
    "IGMinimumCut[graph] gives a minimum edge cut in a weighted graph.\n" <>
    "IGMinimumCut[graph, s, t] gives a minimum s-t edge cut in a weighted graph.";

SyntaxInformation[IGMinimumCut] = {"ArgumentsPattern" -> {_, _., _.}};
IGMinimumCut[graph_?igGraphQ] :=
    catch@Block[{ig = igMake[graph]},
      EdgeList[graph][[ igIndexVec@check@ig@"minCut"[] ]]
    ]
IGMinimumCut[graph_?igGraphQ, s_, t_] :=
    catch@Block[{ig = igMake[graph]},
      EdgeList[graph][[ igIndexVec@check@ig@"minCutST"[vs[graph][s], vs[graph][t]] ]]
    ]


PackageExport["IGGomoryHuTree"]
IGGomoryHuTree::usage = "IGGomoryHuTree[graph] gives the Gomory–Hu tree of a graph.";

(* Note: edge ordering is critical *)
SyntaxInformation[IGGomoryHuTree] = {"ArgumentsPattern" -> {_, OptionsPattern[]}, "OptionNames" -> optNames[Graph]};
IGGomoryHuTree[graph_?igGraphQ, opt : OptionsPattern[]] :=
    catch@Block[{new = igMakeEmpty[], ig = igMake[graph], flow},
      flow = check@new@"gomoryHuTree"[ManagedLibraryExpressionID[ig]];
      igToGraphWithNames[new, VertexList[graph], EdgeWeight -> flow, opt]
    ]


PackageExport["IGFindMinimumCuts"]
IGFindMinimumCuts::usage =
    "IGFindMinimumCuts[graph, s, t] gives all minimum edge cuts that disconnect s and t in a weighted graph.";

SyntaxInformation[IGFindMinimumCuts] = {"ArgumentsPattern" -> {_, _, _}};
IGFindMinimumCuts[graph_?EmptyGraphQ, s_, t_] := {}
(* Note: DirectedGraph preserves EdgeWeights even in old versions such as M10.0. *)
IGFindMinimumCuts[graph_?UndirectedGraphQ, s_, t_] :=
    Apply[UndirectedEdge, IGFindMinimumCuts[DirectedGraph[graph], s, t], {2}]
IGFindMinimumCuts[graph_?igGraphQ, s_, t_] :=
    catch@Block[{ig = igMake[graph]},
      igUnpackEdgeSet[graph]@check@ig@"allMinCutsST"[vs[graph][s], vs[graph][t]]
    ]


PackageExport["IGFindMinimalCuts"]
IGFindMinimalCuts::usage =
    "IGFindMinimalCuts[graph, s, t] gives all minimal edge cuts that disconnect s and t in graph.";

SyntaxInformation[IGFindMinimalCuts] = {"ArgumentsPattern" -> {_, _, _}};
IGFindMinimalCuts[graph_?EmptyGraphQ, s_, t_] := {}
IGFindMinimalCuts[graph_?UndirectedGraphQ, s_, t_] :=
    Apply[UndirectedEdge, IGFindMinimalCuts[DirectedGraph[graph], s, t], {2}]
IGFindMinimalCuts[graph_?igGraphQ, s_, t_] :=
    catch@Block[{ig = igMakeUnweighted[graph]},
      igUnpackEdgeSet[graph]@check@ig@"allCutsST"[vs[graph][s], vs[graph][t]]
    ]


(***** Maximum flow *****)

edgeCapacities[graph_] :=
    With[{capacities = PropertyValue[graph, EdgeCapacity]},
      If[VectorQ[capacities], capacities, {}]
    ]


PackageExport["IGMaximumFlowMatrix"]
IGMaximumFlowMatrix::usage = "IGMaximumFlowMatrix[graph, s, t] gives the flow matrix of a maximum flow from s to t.";

SyntaxInformation[IGMaximumFlowMatrix] = {"ArgumentsPattern" -> {_, _, _}};
IGMaximumFlowMatrix[graph_?igGraphQ, s_, t_] :=
    catch@Block[{ig = igMakeUnweighted[graph], flows},
      flows = check@ig@"maxFlow"[vs[graph][s], vs[graph][t], edgeCapacities[graph]];
      With[{sao = SystemOptions["SparseArrayOptions"]},
        Internal`WithLocalSettings[
          SetSystemOptions["SparseArrayOptions" -> "TreatRepeatedEntries" -> Total]
          ,
          If[DirectedGraphQ[graph],
            SparseArray[IGIndexEdgeList[graph] -> flows, VertexCount[graph] {1, 1}, 0.]
            ,
            With[{pairs = igraphGlobal@"edgeListSortPairs"[IGIndexEdgeList[graph]]},
              SparseArray[Join[pairs, Reverse[pairs, 2]] -> Join[flows, -flows], VertexCount[graph] {1, 1}, 0.]
            ]
          ]
          ,
          SetSystemOptions[sao]
        ]
      ]
    ]


PackageExport["IGMaximumFlowValue"]
IGMaximumFlowValue::usage = "IGMaximumFlowValue[graph, s, t] gives the value of the maximum flow from s to t.";

SyntaxInformation[IGMaximumFlowValue] = {"ArgumentsPattern" -> {_, _, _}};
IGMaximumFlowValue[graph_?igGraphQ, s_, t_] :=
    catch@Block[{ig = igMakeUnweighted[graph]},
      check@ig@"maxFlowValue"[vs[graph][s], vs[graph][t], edgeCapacities[graph]]
    ]


(***** Other functions *****)

PackageExport["IGGiantComponent"]
IGGiantComponent::usage = "IGGiantComponent[graph] gives the largest weakly connected component of graph.";

SyntaxInformation[IGGiantComponent] = {"ArgumentsPattern" -> {_, OptionsPattern[]}, "OptionNames" -> optNames[Graph]};
IGGiantComponent[g_?IGNullGraphQ, opt : OptionsPattern[Graph]] := Graph[g, opt]
IGGiantComponent[g_?GraphQ, opt : OptionsPattern[Graph]] :=
    (* keep Subgraph as-is for v12.0 so IGGiantComponent can preserve properties *)
    Subgraph[g, First@WeaklyConnectedComponents[g], opt]


PackageExport["IGPercolationCurve"]
IGPercolationCurve::usage =
    "IGPercolationCurve[graph] gives a percolation curve corresponding to random edge removal, as {meanDegree, largestComponentFraction} pairs.\n" <>
    "IGPercolationCurve[edges] gives the percolation curve when edges are added in the specified order.\n" <>
    "IGPercolationCurve[edges, n] assumes that there are n vertices.";

SyntaxInformation[IGPercolationCurve] = {"ArgumentsPattern" -> {_, _.}};

IGPercolationCurve::vcount = "The specified vertex count `1` is less than the actual number of vertices `2`.";

IGPercolationCurve[g_?igGraphQ] :=
    igraphGlobal@"percolationCurve"[RandomSample@IGIndexEdgeList[g] - 1, VertexCount[g]]

IGPercolationCurve[edges_List, n : (_?Developer`MachineIntegerQ | Automatic) : Automatic] :=
    Module[{vertices, ind, pairs, vc},
      pairs = edges[[All, {1, 2}]];
      If[Not@MatrixQ[pairs],
        pairs = List @@@ pairs;
      ];
      pairs = Flatten[pairs, 1];
      vertices = Union[pairs];
      ind = AssociationThread[vertices, Range@Length[vertices] - 1];
      pairs = Lookup[ind, pairs];
      vc = Replace[n, Automatic -> Length[vertices]];
      If[vc < Length[vertices],
        vc = Length[vertices];
        Message[IGPercolationCurve::vcount, vc, Length[vertices]];
      ];
      igraphGlobal@"percolationCurve"[Partition[pairs, 2], vc]
    ]


PackageExport["IGSourceVertexList"]
IGSourceVertexList::usage = "IGSourceVertexList[graph] gives the list of vertices with no incoming connections.";

SyntaxInformation[IGSourceVertexList] = {"ArgumentsPattern" -> {_}};
IGSourceVertexList[g_?GraphQ] := Pick[VertexList[g], VertexInDegree[g], 0]


PackageExport["IGSinkVertexList"]
IGSinkVertexList::usage = "IGSinkVertexList[graph] gives the list of vertices with no outgoing connections.";

SyntaxInformation[IGSinkVertexList] = {"ArgumentsPattern" -> {_}};
IGSinkVertexList[g_?GraphQ] := Pick[VertexList[g], VertexOutDegree[g], 0]


PackageExport["IGIsolatedVertexList"]
IGIsolatedVertexList::usage = "IGIsolatedVertexList[graph] gives the list of isolated vertices.";

SyntaxInformation[IGIsolatedVertexList] = {"ArgumentsPattern" -> {_}};
(* Use SimpleGraph so that vertices with only self-loops can be included. *)
IGIsolatedVertexList[g_?GraphQ] := Pick[VertexList[g], VertexDegree@igSimpleGraph[g], 0]
