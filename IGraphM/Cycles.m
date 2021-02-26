(* Mathematica Package *)
(* Created by Mathematica plugin for IntelliJ IDEA *)

(* :Author: szhorvat *)
(* :Date: 2020-03-23 *)
(* :Copyright: (c) 2020 Szabolcs Horvát *)

Package["IGraphM`"]


(*******************************************)
(***** Graph cycles and acyclic graphs *****)
(*******************************************)


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
