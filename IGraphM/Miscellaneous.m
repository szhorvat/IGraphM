(* Mathematica Package *)
(* Created by Mathematica plugin for IntelliJ IDEA *)

(* :Author: szhorvat *)
(* :Date: 2018-10-24 *)
(* :Copyright: (c) 2019 Szabolcs Horvát *)

Package["IGraphM`"]
igContextSetup[igPackagePrivateSymbol]

(***********************************)
(***** Miscellaneous functions *****)
(***********************************)


(***** Directed acyclic graphs *****)

PackageExport["IGDirectedAcyclicGraphQ"]
IGDirectedAcyclicGraphQ::usage = "IGDirectedAcyclicGraphQ[graph] tests if graph is directed and acyclic.";

SyntaxInformation[IGDirectedAcyclicGraphQ] = {"ArgumentsPattern" -> {_}};
IGDirectedAcyclicGraphQ[g_?igGraphQ] := Block[{ig = igMakeFast[g]}, sck@ig@"dagQ"[]]
IGDirectedAcyclicGraphQ[_] := False


PackageExport["IGTopologicalOrdering"]
IGTopologicalOrdering::usage =
    "IGTopologicalOrdering[graph] returns a permutation that sorts the vertices in topological order. " <>
    "Note that the values returned are vertex indices, not vertex names.";

SyntaxInformation[IGTopologicalOrdering] = {"ArgumentsPattern" -> {_}};
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
        Transpose[{Extract[vlist, idx], vlist}],
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
expr : IGCoreness[graph_?igGraphQ, mode_ : "All"] :=
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

PackageExport["IGCactusQ"]
IGCactusQ::usage = "IGCactusQ[graph] tests if graph is a cactus";
SyntaxInformation[IGCactusQ] = {"ArgumentsPattern" -> {_}};
IGCactusQ[graph_ /; UndirectedGraphQ[graph] && SimpleGraphQ[graph]] :=
    Block[{ig = igMakeUnweighted[graph]},
      sck@ig@"cactusQ"[]
    ]
(* Some authors define cactus graphs for the directed case as well. Consider implementing this in the future. *)
IGCactusQ::dirg = "IGCactusQ does not support directed graphs and will return False.";
IGCactusQ[_?DirectedGraphQ] := (Message[IGCactusQ::dirg]; False)
IGCactusQ[_] := False


PackageExport["IGRegularQ"]
IGRegularQ::usage =
    "IGRegularQ[graph] tests if graph is regular, i.e. all vertices have the same degree.\n" <>
    "IGRegularQ[graph, k] tests if graph is k-regular, i.e. all vertices have degree k.";
SyntaxInformation[IGRegularQ] = {"ArgumentsPattern" -> {_, _.}};
IGRegularQ[graph_?igGraphQ] :=
    If[UndirectedGraphQ[graph],
      Equal @@ VertexDegree[graph],
      Equal @@ VertexOutDegree[graph] && VertexInDegree[graph] == VertexOutDegree[graph]
    ]
IGRegularQ[graph_?EmptyGraphQ, k_Integer?NonNegative] := k == 0 (* required because First@VertexDegree[g] is not usable for the null graph *)
IGRegularQ[graph_?igGraphQ, k_Integer?NonNegative] :=
    If[UndirectedGraphQ[graph],
      With[{vd = VertexDegree[graph]},
        First[vd] == k && Equal @@ vd
      ]
      ,
      With[{vod = VertexOutDegree[graph], vid = VertexInDegree[graph]},
        First[vod] == k && Equal @@ vod && vid == vod
      ]
    ]
IGRegularQ[_] := False


PackageExport["IGCompleteQ"]
IGCompleteQ::usage = "IGCompleteQ[graph] tests if all pairs of vertices are connected in graph.";

igCompleteQ[graph_?SimpleGraphQ] :=
    With[{vc = VertexCount[graph], ec = EdgeCount[graph]},
      If[UndirectedGraphQ[graph],
        vc*(vc-1)/2 == ec,
        vc*(vc-1) == ec
      ]
    ]
igCompleteQ[graph_] := igCompleteQ@SimpleGraph[graph]

SyntaxInformation[IGCompleteQ] = {"ArgumentsPattern" -> {_}};
IGCompleteQ[graph_?igGraphQ] := igCompleteQ[graph]
IGCompleteQ[_] := False


PackageExport["IGExpressionTree"]
IGExpressionTree::usage = "IGExpressionTree[expression] constructs a tree graph from an arbitrary Mathematica expression.";

SetAttributes[exprHead, HoldFirst]
exprHead[expr_Association] := HoldForm[Association] (* Associations are atoms, but Position does descend into them *)
exprHead[expr_ /; AtomQ@Unevaluated[expr]] := HoldForm[expr]
exprHead[expr_] := Extract[Unevaluated[expr], 0, HoldForm]

Options[IGExpressionTree] = {
  VertexLabels -> "Head",
  GraphLayout -> "LayeredEmbedding",
  DirectedEdges -> True
};
SyntaxInformation[IGExpressionTree] = {"ArgumentsPattern" -> {_, OptionsPattern[]}, "OptionNames" -> optNames[IGExpressionTree, Graph]};
(* TODO: Is there an AtomQ that Position will enter into? If yes, exprHead[] will label it incorrectly. *)
IGExpressionTree[expr_, opt : OptionsPattern[{IGExpressionTree, Graph}]] :=
    Module[{vertices, edges, vertexLabels},
      (* We could reverse the vertex list to make the root the first node rather than the last one.
         But this would cause the leaves to be displayed in reverse order with LayeredEmbedding. *)
      vertices = Position[expr, _, {0, Infinity}, Heads -> False];
      edges = MapThread[Rule, {vertices[[;; -2, ;; -2]], vertices[[;; -2]]}];
      vertexLabels =
          Replace[
            OptionValue[VertexLabels],
            {
              "Head" :> Thread[vertices -> Extract[expr, vertices, exprHead]],
              Placed["Head", rest___] :> Thread[vertices -> Extract[expr, vertices, Placed[#, rest]& @* exprHead]],
              "Subexpression" :> Thread[vertices -> Extract[expr, vertices, HoldForm]],
              Placed["Subexpression", rest___] :> Thread[vertices -> Extract[expr, vertices, Placed[#, rest]& @* HoldForm]]
            }
          ];
      Graph[vertices, edges,
        VertexLabels -> vertexLabels,
        opt,
        GraphRoot -> {},
        GraphLayout -> OptionValue[GraphLayout],
        DirectedEdges -> OptionValue[DirectedEdges]
      ]
    ]
