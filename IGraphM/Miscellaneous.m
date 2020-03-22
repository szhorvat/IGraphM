(* Mathematica Package *)
(* Created by Mathematica plugin for IntelliJ IDEA *)

(* :Author: szhorvat *)
(* :Date: 2018-10-24 *)
(* :Copyright: (c) 2018-2020 Szabolcs Horvát *)

Package["IGraphM`"]
igContextSetup[igPackagePrivateSymbol]

(***********************************)
(***** Miscellaneous functions *****)
(***********************************)


(***** Directed acyclic graphs *****)

PackageExport["IGDirectedAcyclicGraphQ"]
IGDirectedAcyclicGraphQ::usage = "IGDirectedAcyclicGraphQ[graph] tests if graph is directed and acyclic.";

SyntaxInformation[IGDirectedAcyclicGraphQ] = {"ArgumentsPattern" -> {_}};
(* Edgeless graphs are considered both directed and undirected by Mathematica,
   but are transferred as undirected to igraph. Therefore dagQ() would return
   False. We catch these early and return True. *)
IGDirectedAcyclicGraphQ[g_?EmptyGraphQ] := True
IGDirectedAcyclicGraphQ[g_?igGraphQ] := Block[{ig = igMakeFast[g]}, sck@ig@"dagQ"[]]
IGDirectedAcyclicGraphQ[_] := False


PackageExport["IGTopologicalOrdering"]
IGTopologicalOrdering::usage =
    "IGTopologicalOrdering[graph] returns a permutation that sorts the vertices in topological order. " <>
    "Note that the values returned are vertex indices, not vertex names.";

SyntaxInformation[IGTopologicalOrdering] = {"ArgumentsPattern" -> {_}};
(* Catch edgeless graphs early. See comment for IGDirectedAcyclicGraphQ[] *)
IGTopologicalOrdering[graph_?EmptyGraphQ] := Range@VertexCount[graph]
IGTopologicalOrdering[graph_?igGraphQ] :=
    catch@Block[{ig = igMakeFast[graph]},
      igIndexVec@check@ig@"topologicalSorting"[]
    ]


PackageExport["IGFeedbackArcSet"]
IGFeedbackArcSet::usage = "IGFeedbackArcSet[graph] computes a feedback edge set of graph. Removing these edges makes the graph acyclic.";

igFeedbackArcSetMethods = <| "IntegerProgramming" -> True, "EadesLinSmyth" -> False |>;

Options[IGFeedbackArcSet] = { Method -> "IntegerProgramming" };
SyntaxInformation[IGFeedbackArcSet] = {"ArgumentsPattern" -> {_, OptionsPattern[]} };

IGFeedbackArcSet::bdmtd =
    "Value of option Method -> `` is not one of " <>
    ToString[Keys[igFeedbackArcSetMethods], InputForm] <> ".";

amendUsage[IGFeedbackArcSet, "Available Method options: <*Keys[igFeedbackArcSetMethods]*>. \"IntegerProgramming\" is guaranteed to find a minimum feedback arc set."];

IGFeedbackArcSet[graph_?igGraphQ, opt : OptionsPattern[]] :=
    catch@Block[{ig = igMake[graph]}, (* use igMake because edge ordering matters *)
      Part[
        EdgeList[graph],
        igIndexVec@check@ig@"feedbackArcSet"[Lookup[igFeedbackArcSetMethods, OptionValue[Method], Message[IGFeedbackArcSet::bdmtd, OptionValue[Method]]; throw[$Failed]]]
      ]
    ]


(***** Voronoi cells with graph distance *****)

PackageExport["IGVoronoiCells"]
IGVoronoiCells::usage = "IGVoronoiCells[graph, {v1, v2, \[Ellipsis]}] find the sets of vertices closest to each given vertex.";

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
    "IGGraphicalQ[indegrees, outdegrees] tests if indegrees with outdegrees is the degree sequence of any simple directed graph.";

SyntaxInformation[IGGraphicalQ] = {"ArgumentsPattern" -> {_, _.}};
IGGraphicalQ[degrees_?nonNegIntVecQ] := sck@igraphGlobal@"erdosGallai"[degrees] (* use fast custom implementation instead of igraph *)
IGGraphicalQ[indeg_?nonNegIntVecQ, outdeg_?nonNegIntVecQ] := sck@igraphGlobal@"graphicalQ"[outdeg, indeg]
IGGraphicalQ[___] := False


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


PackageExport["IGAverageNeighborDegree"]
IGAverageNeighborDegree::usage =
    "IGAverageNeighborDegree[graph] gives the average neighbour degree of the vertices of graph.\n" <>
    "IGAverageNeighborDegree[graph, {vertex1, vertex2, \[Ellipsis]}] gives the average neighbour degree of the specified vertices.\n" <>
    "IGAverageNeighborDegree[graph, All, mode] uses the given mode, \"In\", \"Out\" or \"All\", to find neighbours and degrees in directed graphs. The default is \"Out\".\n" <>
    "IGAverageNeighborDegree[graph, All, degreeMode, neighborMode] uses different modes for finding neighbours and degrees.";
SyntaxInformation[IGAverageNeighborDegree] = {"ArgumentsPattern" -> {_, _., _.}};
IGAverageNeighborDegree[graph_?igGraphQ, {}, degMode_, neiMode_] := {}
IGAverageNeighborDegree[graph_?igGraphQ, vs : (_List | All) : All, degMode_ : "Out", neiMode_ : Automatic] :=
    catch@Module[{ig = igMakeFastWeighted[graph], nMode = If[neiMode === Automatic, degMode, neiMode]},
      check@ig@"averageNeighborDegree"[vss[graph][vs], encodeNeighborMode[nMode], encodeNeighborMode[degMode]]
    ]
addCompletion[IGAverageNeighborDegree, {0, 0, {"In", "Out", "All"}, {"In", "Out", "All"}}]


PackageExport["IGAverageDegreeConnectivity"]
IGAverageDegreeConnectivity::usage =
    "IGAverageDegreeConnectivity[graph] gives the average neighbour degree for vertices of degree k=1, 2, \[Ellipsis]\n" <>
    "IGAverageDegreeConnectivity[graph, mode] uses the given mode, \"In\", \"Out\" or \"All\", to find neighbours and degrees in directed graphs. The default is \"Out\".\n" <>
    "IGAverageDegreeConnectivity[graph, degreeMode, neighborMode] uses different modes for finding neighbours and degrees.";
SyntaxInformation[IGAverageDegreeConnectivity] = {"ArgumentsPattern" -> {_, _., _.}};
IGAverageDegreeConnectivity[graph_?igGraphQ, degMode_ : "Out", neiMode_ : Automatic] :=
    catch@Block[{ig = igMakeFastWeighted[graph], nMode =If[neiMode === Automatic, degMode, neiMode] },
      check@ig@"averageDegreeConnectivity"[encodeNeighborMode[nMode], encodeNeighborMode[degMode]]
    ]
addCompletion[IGAverageDegreeConnectivity, {0, {"In", "Out", "All"}, {"In", "Out", "All"}}]
