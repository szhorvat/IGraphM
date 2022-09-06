(* Mathematica Package *)
(* Created by Mathematica plugin for IntelliJ IDEA *)

(* :Author: szhorvat *)
(* :Date: 2018-10-24 *)
(* :Copyright: (c) 2018-2022 Szabolcs Horvát *)

Package["IGraphM`"]


(**************************)
(***** Chordal graphs *****)
(**************************)


PackageExport["IGChordalQ"]
IGChordalQ::usage = "IGChordalQ[graph] tests if graph is chordal.";

SyntaxInformation[IGChordalQ] = {"ArgumentsPattern" -> {_}};
IGChordalQ[graph_?igGraphQ] :=
    Block[{ig = igMakeUnweighted[graph]}, ig@"chordalQ"[]]
IGChordalQ[_] := False


PackageExport["IGMaximumCardinalitySearch"]
IGMaximumCardinalitySearch::usage = "IGMaximumCardinalitySearch[graph] assigns a rank to each vertex, from 1 to n, according to the maximum cardinality search algorithm. Visiting the vertices of the graph by decreasing rank is equivalent to always visiting the next vertex with the most already visited neighbours.";

SyntaxInformation[IGMaximumCardinalitySearch] = {"ArgumentsPattern" -> {_}};
IGMaximumCardinalitySearch[graph_?igGraphQ] :=
    catch@Block[{ig = igMakeUnweighted[graph]},
      igIndexVec@check@ig@"maximumCardinalitySearch"[]
    ]

PackageExport["IGChordalCompletion"]
IGChordalCompletion::usage = "IGChordalCompletion[graph] gives a set of edges that, when added to graph, make it chordal. The edge-set this function returns is usually not minimal.";

SyntaxInformation[IGChordalCompletion] = {"ArgumentsPattern" -> {_}};
IGChordalCompletion[graph_?igGraphQ] :=
    catch@Block[{ig = igMakeUnweighted[graph], result},
      result = check@ig@"chordalCompletion"[];
      If[DirectedGraphQ[graph], DirectedEdge, UndirectedEdge] @@@ Partition[igVertexNames[graph]@igIndexVec[result], 2]
    ]
