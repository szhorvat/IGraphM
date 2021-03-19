(* Mathematica Package *)
(* Created by Mathematica plugin for IntelliJ IDEA *)

(* :Author: szhorvat *)
(* :Date: 2018-10-24 *)
(* :Copyright: (c) 2018-2020 Szabolcs Horvát *)

Package["IGraphM`"]


(***********************************)
(***** Miscellaneous functions *****)
(***********************************)

(***** Voronoi cells with graph distance *****)

PackageExport["IGVoronoiCells"]
IGVoronoiCells::usage = "IGVoronoiCells[graph, {v1, v2, \[Ellipsis]}] returns the sets of vertices closest to each given vertex.";

IGVoronoiCells::ivert = "The given centers `1` are not vertices of the graph.";
Options[IGVoronoiCells] = { "Tiebreaker" -> Automatic };
SyntaxInformation[IGVoronoiCells] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};
IGVoronoiCells[g_?igGraphQ, centers_List, opt : OptionsPattern[]] :=
    Module[{clist = DeleteDuplicates[centers], vlist = VertexList[g], tiebreaker = OptionValue["Tiebreaker"], idx, dmat},
      If[Not@SubsetQ[vlist, clist],
        Message[IGVoronoiCells::ivert, Complement[clist, vlist]];
        Return[$Failed]
      ];
      dmat = Transpose@IGDistanceMatrix[g, centers];
      idx = If[MatchQ[tiebreaker, Automatic|First],
        Ordering[#, 1]& /@ dmat,
        With[{min = Min[#]}, tiebreaker@Position[#, min]]& /@ dmat
      ];
      GroupBy[
        Transpose[{Extract[clist, idx], vlist}],
        First -> Last
      ]
    ]


(***** k-cores *****)

PackageExport["IGCoreness"]
IGCoreness::usage =
    "IGCoreness[graph] returns the coreness of each vertex. Coreness is the highest order of a k-core containing the vertex.\n" <>
    "IGCoreness[graph, \"In\"] considers only in-degrees in a directed graph.\n" <>
    "IGCoreness[graph, \"Out\"] considers only out-degrees in a directed graph.";

corenessModes = <|"In" -> -1, "Out" -> 1, "All" -> 0|>;
SyntaxInformation[IGCoreness] = {"ArgumentsPattern" -> {_, _.}};
expr : IGCoreness[graph_?igGraphQ, mode_String : "All"] :=
    catch@Block[{ig = igMakeFast[graph]},
      Round@check@ig@"coreness"[Lookup[corenessModes, mode, Message[IGCoreness::inv, HoldForm@OutputForm[expr], mode, "parameter"]; throw[$Failed]]]
    ]
addCompletion[IGCoreness, {0, {"In", "Out", "All"}}]


(***** Degree sequences *****)

PackageExport["IGGraphicalQ"]
IGGraphicalQ::usage =
    "IGGraphicalQ[degrees] tests if degrees is the degree sequence of any simple undirected graph.\n" <>
    "IGGraphicalQ[indegrees, outdegrees] tests if indegrees with outdegrees is the degree sequence of any simple directed graph.\n" <>
    "IGGraphicalQ[degrees, SelfLoops -> True] tests if degrees is the degree sequence of any undirected graph with at most one self-loop per vetrex.\n" <>
    "IGGraphicalQ[degrees, MultiEdges -> True] tests if degrees is the degree sequence of any undirected loop-free multigraph.";
Options[IGGraphicalQ] = { SelfLoops -> False, MultiEdges -> False };
SyntaxInformation[IGGraphicalQ] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};
IGGraphicalQ[{}, {}, opt : OptionsPattern[]] := True
IGGraphicalQ[indeg_?nonNegIntVecQ, outdeg_?nonNegIntVecQ, opt : OptionsPattern[]] :=
    sck@igraphGlobal@"graphicalQ"[outdeg, indeg, True, OptionValue[SelfLoops], OptionValue[MultiEdges]]
IGGraphicalQ[degrees_?nonNegIntVecQ, opt : OptionsPattern[]] :=
    sck@igraphGlobal@"graphicalQ"[degrees, {}, False, OptionValue[SelfLoops], OptionValue[MultiEdges]]
(*    sck@igraphGlobal@"erdosGallai"[degrees] *)(* use fast custom implementation instead of igraph *)
IGGraphicalQ[___] := False


PackageExport["IGBigraphicalQ"]
IGBigraphicalQ::usage =
    "IGBigraphicalQ[degrees1, degrees2] tests if (degrees1, degrees2) is the degree sequence of any bipartite simple graph.\n" <>
    "IGBigraphicalQ[degrees1, degrees2, MultiEdges -> True] tests if (degrees1, degrees2) is the degree sequence of any bipartite multigraph.";
Options[IGBigraphicalQ] = { MultiEdges -> False };
SyntaxInformation[IGBigraphicalQ] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};
IGBigraphicalQ[deg1_?nonNegIntVecQ, deg2_?nonNegIntVecQ, opt : OptionsPattern[]] :=
    sck@igraphGlobal@"bigraphicalQ"[deg1, deg2, OptionValue[MultiEdges]]


PackageExport["IGPotentiallyConnectedQ"]
IGPotentiallyConnectedQ::usage =
    "IGPotentiallyConnectedQ[degrees] tests if degrees is the degree sequence of any connected graph.\n" <>
    "IGPotentiallyConnectedQ[indegrees, outdegrees] tests if indegrees with outdegrees is the degree sequence of any strongly connected directed graph.";
SyntaxInformation[IGPotentiallyConnectedQ] = {"ArgumentsPattern" -> {_, _.}};
(* Undirected case:
 * (no. of edges) >= (no. of vertices) - 1
 * All degrees are at least 1.
 * The singleton graph, which has one zero-degree vertex, is an exception.
 *)
IGPotentiallyConnectedQ[{}] := False (* the null graph is disconnected *)
IGPotentiallyConnectedQ[{0}] := True (* the singleton graph is connected *)
IGPotentiallyConnectedQ[degrees_?nonNegIntVecQ] :=
    With[{tot = Total[degrees], len = Length[degrees]},
      EvenQ[tot] && tot >= 2*(len-1) && Min[degrees] > 0 (* note: do not use FreeQ[degrees, 0, {1}] because it unpacks *)
    ]
(* Directed case:
 * All in- and out-degrees must be at least 1, see e.g. https://dx.doi.org/10.1186%2Fs13660-017-1544-3, Theorem 3.1 condition (ii).
 * The singleton graph, which has one zero-degree vertex, is an exception.
 *)
IGPotentiallyConnectedQ[{}, {}] := False (* the null graph is disconnected *)
IGPotentiallyConnectedQ[{0}, {0}] := True (* the singleton graph is connected *)
IGPotentiallyConnectedQ[indeg_?nonNegIntVecQ, outdeg_?nonNegIntVecQ] :=
    Total[indeg] == Total[outdeg] && Min[indeg] > 0 && Min[outdeg] > 0
IGPotentiallyConnectedQ[___] := False


PackageExport["IGSplitQ"]
IGSplitQ::usage =
    "IGSplitQ[graph] tests if graph is a split graph.\n" <>
    "IGSplitQ[degrees] tests if degrees is the degree sequence of a split graph.";

IGSplitQ::dirg = "Directed graphs are not currently supported.";

SyntaxInformation[IGSplitQ] = {"ArgumentsPattern" -> {_}};
IGSplitQ[graph_?UndirectedGraphQ] :=
    sck@igraphGlobal@"splitQ"[ VertexDegree@SimpleGraph[graph] ]
IGSplitQ[graph_?DirectedGraphQ] := (Message[IGSplitQ::dirg]; $Failed)
IGSplitQ[graph_?MixedGraphQ] := (Message[IGraphM::mixed]; $Failed)
IGSplitQ[degrees_?nonNegIntVecQ] := sck@igraphGlobal@"splitQ"[degrees]
IGSplitQ[_] := False


PackageExport["IGThresholdQ"]
IGThresholdQ::usage =
    "IGThresholdQ[graph] tests if graph is a threshold graph.\n" <>
    "IGThresholdQ[degrees] tests if degrees form a threshold degree sequence.";

IGThresholdQ::dirg = "Directed graphs are not currently supported.";

SyntaxInformation[IGThresholdQ] = {"ArgumentsPattern" -> {_}};
IGThresholdQ[graph_?UndirectedGraphQ] :=
    sck@igraphGlobal@"thresholdQ"[ VertexDegree@SimpleGraph[graph] ]
IGThresholdQ[graph_?DirectedGraphQ] := (Message[IGThresholdQ::dirg]; $Failed)
IGThresholdQ[graph_?MixedGraphQ] := (Message[IGraphM::mixed]; $Failed)
IGThresholdQ[degrees_?nonNegIntVecQ] := sck@igraphGlobal@"thresholdQ"[degrees]


(***** Dominators *****)

PackageExport["IGDominatorTree"]
IGDominatorTree::usage = "IGDominatorTree[graph, root] returns the dominator tree of a directed graph starting from root.";
SyntaxInformation[IGDominatorTree] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}, "OptionNames" -> optNames[Graph]};
IGDominatorTree[graph_?igGraphQ, root_, opt : OptionsPattern[Graph]] :=
    catch@Block[{new = igMakeEmpty[], ig = igMakeUnweighted[graph]},
      check@new@"dominatorTree"[ManagedLibraryExpressionID[ig], vs[graph][root]];
      igToGraphWithNames[new, VertexList[graph], opt, GraphRoot -> root]
    ]


PackageExport["IGImmediateDominators"]
IGImmediateDominators::usage = "IGImmediateDominators[graph, root] returns the immediate dominator of each vertex relative to root.";
SyntaxInformation[IGImmediateDominators] = {"ArgumentsPattern" -> {_, _}};
IGImmediateDominators[graph_?igGraphQ, root_] :=
    catch@Block[{ig = igMakeUnweighted[graph], dominators, mask},
      dominators = check@ig@"immediateDominators"[vs[graph][root]];
      mask = UnitStep[dominators];
      AssociationThread[Pick[VertexList[graph], mask, 1], VertexList[graph][[Pick[igIndexVec[dominators], mask, 1]]]]
    ]


(***** Neighbour degrees and assortativity *****)

PackageExport["IGAverageNeighborDegree"]
IGAverageNeighborDegree::usage =
    "IGAverageNeighborDegree[graph] gives the average neighbour degree of the vertices of graph.\n" <>
    "IGAverageNeighborDegree[graph, {vertex1, vertex2, \[Ellipsis]}] gives the average neighbour degree of the specified vertices.\n" <>
    "IGAverageNeighborDegree[graph, All, mode] uses the given mode, \"In\", \"Out\" or \"All\", to find neighbours and degrees in directed graphs. The default is \"Out\".\n" <>
    "IGAverageNeighborDegree[graph, All, degreeMode, neighborMode] uses different modes for finding neighbours and degrees.";

SyntaxInformation[IGAverageNeighborDegree] = {"ArgumentsPattern" -> {_, _., _., _.}};
IGAverageNeighborDegree[graph_?igGraphQ, {}, degMode_String : "Out", neiMode : _String | Automatic : Automatic] := {}
IGAverageNeighborDegree[graph_?igGraphQ, vs : (_List | All) : All, degMode_String : "Out", neiMode : _String | Automatic : Automatic] :=
    catch@Block[{ig = igMakeFastWeighted[graph]},
      expectInfNaN@fixInfNaN@check@ig@"averageNeighborDegree"[vss[graph][vs], encodeNeighborMode@Replace[neiMode, Automatic -> degMode], encodeNeighborMode[degMode]]
    ]
addCompletion[IGAverageNeighborDegree, {0, 0, {"In", "Out", "All"}, {"In", "Out", "All"}}]


PackageExport["IGAverageDegreeConnectivity"]
IGAverageDegreeConnectivity::usage =
    "IGAverageDegreeConnectivity[graph] gives the average neighbour degree for vertices of degree k=1, 2, \[Ellipsis]\n" <>
    "IGAverageDegreeConnectivity[graph, mode] uses the given mode, \"In\", \"Out\" or \"All\", to find neighbours and degrees in directed graphs. The default is \"Out\".\n" <>
    "IGAverageDegreeConnectivity[graph, degreeMode, neighborMode] uses different modes for finding neighbours and degrees.";

SyntaxInformation[IGAverageDegreeConnectivity] = {"ArgumentsPattern" -> {_, _., _.}};
IGAverageDegreeConnectivity[graph_?igGraphQ, degMode_String : "Out", neiMode : _String | Automatic : Automatic] :=
    catch@Block[{ig = igMakeFastWeighted[graph]},
      expectInfNaN@fixInfNaN@check@ig@"averageDegreeConnectivity"[encodeNeighborMode@Replace[neiMode, Automatic -> degMode], encodeNeighborMode[degMode]]
    ]
addCompletion[IGAverageDegreeConnectivity, {0, {"In", "Out", "All"}, {"In", "Out", "All"}}]


(***** Eulerian paths *****)

PackageExport["IGEulerianQ"]
IGEulerianQ::usage =
    "IGEulerianQ[graph] tests if graph has a path that traverses each edge once (Eulerian path).\n" <>
    "IGEulerianQ[graph, Closed -> True] tests if graph has a cycle that traverses each edge once.";
Options[IGEulerianQ] = { Closed -> False };
SyntaxInformation[IGEulerianQ] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};
IGEulerianQ[graph_?igGraphQ, opt : OptionsPattern[]] :=
    Block[{ig = igMakeUnweighted[graph]},
      sck@ig@"eulerianQ"[OptionValue[Closed]]
    ]

PackageExport["IGEulerianPath"]
IGEulerianPath::usage =
    "IGEulerianPath[graph] returns the edges of an Eulerian path, if it exists.\n" <>
    "IGEulerianPath[graph, Closed -> True] returns an Eulerian cycle.";
Options[IGEulerianPath] = { Closed -> False };
SyntaxInformation[IGEulerianPath] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};
IGEulerianPath[graph_?igGraphQ, opt : OptionsPattern[]] :=
    catch@Block[{ig = igMakeUnweighted[graph]},
      EdgeList[graph][[ igIndexVec@check@ig@"eulerianPath"[OptionValue[Closed]] ]]
    ]

PackageExport["IGEulerianPathVertices"]
IGEulerianPathVertices::usage =
    "IGEulerianPathVertices[graph] returns the vertices of an Eulerian path, if it exists.\n" <>
    "IGEulerianPathVertices[graph, Closed -> True] returns the vertices of an Eulerian cycle.";
Options[IGEulerianPathVertices] = { Closed -> False };
SyntaxInformation[IGEulerianPathVertices] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};
IGEulerianPathVertices[graph_?igGraphQ, opt : OptionsPattern[]] :=
    catch@Block[{ig = igMakeUnweighted[graph]},
      VertexList[graph][[ igIndexVec@check@ig@"eulerianPathVertices"[OptionValue[Closed]] ]]
    ]


(***** Other functions *****)

PackageExport["IGTreelikeComponents"]
IGTreelikeComponents::usage = "IGTreelikeComponents[graph] returns the vertices that make up tree-like components.";

SyntaxInformation[IGTreelikeComponents] = {"ArgumentsPattern" -> {_}};
IGTreelikeComponents[graph_?igGraphQ] :=
    catch@Block[{ig = igMakeFast[graph]},
      igVertexNames[graph]@igIndexVec@check@ig@"treelikeComponents"[]
    ]


(* Notes on graph types:

   Directed: $Failed because some authors do define directed cacti, and it could be implemented.
   Self-loops: Ignore them.
   Multi-edges: Support.
*)
PackageExport["IGCactusQ"]
IGCactusQ::usage = "IGCactusQ[graph] tests if graph is a cactus";
SyntaxInformation[IGCactusQ] = {"ArgumentsPattern" -> {_}};
IGCactusQ[graph_?UndirectedGraphQ] :=
    catch@Block[{ig = igMakeUnweighted[graph]},
      If[SimpleGraphQ[graph],
        check@ig@"cactusQ"[],
        check@ig@"nonSimpleCactusQ"[]
      ]
    ]
(* Some authors define cactus graphs for the directed case as well. Consider implementing this in the future. *)
IGCactusQ::dirg = "IGCactusQ is not implemented for directed graphs.";
IGCactusQ[_?DirectedGraphQ] := (Message[IGCactusQ::dirg]; $Failed)
IGCactusQ[_] := False


PackageExport["IGFundamentalCycles"]
IGFundamentalCycle::usage =
    "IGFundamentalCycle[graph]\n" <>
    "IGFundamentalCycle[{graph, vertex}]\n";
Options[IGFundamentalCycles] = { "Cutoff" -> Infinity };
SyntaxInformation[IGFundamentalCycles] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};
IGFundamentalCycles[{graph_?igGraphQ, v_}, opt : OptionsPattern[]] :=
    catch@Block[{ig = igMakeUnweighted[graph]},
       igUnpackIndices@check@ig@"fundamentalCycles"[vs[graph][v], infToNeg@OptionValue["Cutoff"]]
    ]
IGFundamentalCycles[graph_?igGraphQ, opt : OptionsPattern[]] :=
    catch@Block[{ig = igMakeUnweighted[graph]},
       igUnpackIndices@check@ig@"fundamentalCycles"[-1, infToNeg@OptionValue["Cutoff"]]
    ]


PackageExport["IGMinimumCycleBasis"]
IGMinimumCycleBasis::usage = "IGMinimumCycleBasis[graph]";
Options[IGMinimumCycleBasis] = { "Cutoff" -> Infinity, "Complete" -> True };
SyntaxInformation[IGMinimumCycleBasis] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};
IGMinimumCycleBasis[graph_?igGraphQ, opt : OptionsPattern[]] :=
    catch@Block[{ig = igMakeUnweighted[graph]},
       igUnpackIndices@check@ig@"minimumCycleBasis"[infToNeg@OptionValue["Cutoff"], OptionValue["Complete"]]
    ]
