(* Mathematica Package *)
(* Created by Mathematica plugin for IntelliJ IDEA *)

(* :Author: szhorvat *)
(* :Date: 2018-10-24 *)
(* :Copyright: (c) 2018 Szabolcs Horvát *)

Package["IGraphM`"]
igContextSetup[igPackagePrivateSymbol]

(***************************)
(***** Graph colouring *****)
(***************************)

PackageExport["IGVertexColoring"]
IGVertexColoring::usage = "IGVertexColoring[graph] returns a vertex colouring of graph.";

SyntaxInformation[IGVertexColoring] = {"ArgumentsPattern" -> {_}};
IGVertexColoring[graph_?igGraphQ] :=
    catch@Block[{ig = igMakeFast[graph]},
      1 + check@ig@"vertexColoring"[]
    ]


PackageExport["IGEdgeColoring"]
IGEdgeColoring::usage = "IGEdgeColoring[graph] returns an edge colouring of graph.";

SyntaxInformation[IGEdgeColoring] = {"ArgumentsPattern" -> {_, _}};
IGEdgeColoring[graph_?igGraphQ] := IGVertexColoring@LineGraph@IGUndirectedGraph[graph, "All"]


PackageExport["IGKVertexColoring"]
IGKVertexColoring::usage = "IGKVertexColoring[graph, k] attempts to find a k-colouring of graph's vertices. If none exist, {} is returned.";

Options[IGKVertexColoring] = { "ForcedColoring" -> Automatic };
SyntaxInformation[IGKVertexColoring] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};
IGKVertexColoring[graph_?EmptyGraphQ, k_Integer?Positive, OptionsPattern[]] := {ConstantArray[1, VertexCount[graph]]}
IGKVertexColoring[graph_?igGraphQ, 2, OptionsPattern[]] :=
    Block[{ig = igMakeFast[graph], parts},
      parts = ig@"bipartitePartitions"[];
      If[MatchQ[parts, _LibraryFunctionError],
        {},
        {parts + 1}
      ]
    ]
IGKVertexColoring[graph_?igGraphQ, k_Integer?Positive, OptionsPattern[]] :=
    catch@igKVertexColoring[graph, k, igKColoringHeuristic[OptionValue["ForcedColoring"]][graph]]


igKVertexColoring[graph_, k_, clique_] :=
    Module[{a, n, edges, res, satExpr, fixedVars},
      If[Length[clique] > k, Return[{}]];
      n = k VertexCount[graph];
      edges = Complement[
        Sort /@ IGIndexEdgeList@SimpleGraph[graph],
        Sort /@ Subsets[clique, {2}]
      ];
      satExpr = And @@ Flatten[{
        Or @@@ Delete[Partition[a /@ Range[n], k], List /@ clique],
        MapThread[
          Or,
          {Not /@ a /@ (Range[k] + k #1),
           Not /@ a /@ (Range[k] + k #2)}
        ] & @@@ (edges - 1)
      }];
      fixedVars = Flatten@MapAt[
        Not,
        Map[Not@a[#] &, Range[k (clique - 1) + 1, k clique], {2}],
        Table[{i, i}, {i, Length[clique]}]
      ];
      satExpr = And[And @@ fixedVars, satExpr];
      res = SatisfiabilityInstances[
        satExpr,
        a /@ Range[n]
      ];
      If[res === {},
        {},
        Transpose[FirstPosition[#, True] & /@ Partition[First[res], k]]
      ]
    ]
(* TODO: The graph is always (1+Max@VertexDegree[graph])-colourable. Look up an algorithm that can manage this case, and use it. *)

(* This is the original implementation of igKVertexColoring.
 * It only adds simple symmetry breaking constraints. *)
(*
igKVertexColoring[graph_, k_, clique_] :=
    Module[{a, n, res, satExpr, fixedVars},
      If[Length[clique] > k, throw[{}]];
      n = k VertexCount[graph];
      satExpr = And @@ Flatten[{
        Or @@@ Partition[a /@ Range[n], k],
        MapThread[
          Or,
          {Not /@ a /@ (Range[k] + k #1),
            Not /@ a /@ (Range[k] + k #2)}
        ] & @@@ (IGIndexEdgeList@SimpleGraph[graph] - 1)
      }];
      fixedVars = a /@ (k (clique-1) + Range@Length[clique]);
      satExpr = And[satExpr, And @@ fixedVars];
      res = SatisfiabilityInstances[
        satExpr,
        a /@ Range[n]
      ];
      If[res === {},
        {},
        Transpose[FirstPosition[#, True] & /@ Partition[First[res], k]]
      ]
    ] *)

igKColoringHeuristic[Automatic][graph_] :=
    Module[{mdcl, lcl},
      mdcl = igKColoringHeuristic["MaxDegreeClique"][graph];
      lcl = TimeConstrained[igKColoringHeuristic["LargestClique"][graph], 0.1];
      If[ListQ[lcl] && Length[lcl] > Length[mdcl],
        lcl,
        mdcl
      ]
    ]
(* The "MaxDegreeClique" heuristic finds a clique containing a max degree node,
 * and also containing some other high-degree nodes, if possible.
 * With Mathematica's SAT solver, ordering the clique's vertices by decreasing degree
 * seems to improve performance. *)
igKColoringHeuristic["MaxDegreeClique"][graph_] :=
    Module[{g = IndexGraph@UndirectedGraph[graph], vd, v, sg, cl},
      vd = VertexDegree[g];
      v = Extract[VertexList[g], Ordering[vd, -1]];
      sg = Subgraph[g, Prepend[AdjacencyList[g, v], v]];
      cl = First@MaximalBy[
        FindClique[sg, Length /@ FindClique[sg], 512 (* limit number of cliques *)],
        Reverse@Sort@Part[vd, #] &
      ];
      Reverse@SortBy[cl, vd[[#]]&]
    ]
igKColoringHeuristic["LargestClique"][graph_] :=
    With[{g = IndexGraph@UndirectedGraph[graph]},
      First@FindClique[g]
    ]
igKColoringHeuristic[None][graph_] := {}
igKColoringHeuristic[clique_List][graph_] := Check[VertexIndex[graph, #]& /@ clique, throw[$Failed]]
igKColoringHeuristic[_][graph_] := {}


PackageExport["IGKEdgeColoring"]
IGKEdgeColoring::usage = "IGKEdgeColoring[graph, k] attempts of find a k-colouring of graph's edges. If none exist, {} is returned.";

Options[IGKEdgeColoring] = { "ForcedColoring" -> "MaxDegreeClique" };
SyntaxInformation[IGKEdgeColoring] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};
IGKEdgeColoring[graph_?igGraphQ, k_Integer?Positive, opt : OptionsPattern[]] :=
    catch@Module[{fc = OptionValue["ForcedColoring"]},
      If[ListQ[fc],
        fc = Check[EdgeIndex[graph, #]& /@ fc, throw[$Failed]]
      ];
      IGKVertexColoring[LineGraph@IGUndirectedGraph[graph, "All"], k, "ForcedColoring" -> fc]
    ]


PackageExport["IGMinimumVertexColoring"]
IGMinimumVertexColoring::usage = "IGMinimumVertexColoring[graph] finds a minimum vertex colouring of graph.";

Options[IGMinimumVertexColoring] = { "ForcedColoring" -> Automatic }
SyntaxInformation[IGMinimumVertexColoring] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};
IGMinimumVertexColoring[graph_?EmptyGraphQ, OptionsPattern[]] := ConstantArray[1, VertexCount[graph]]
IGMinimumVertexColoring[graph_?igGraphQ, opt : OptionsPattern[]] :=
    catch@Module[{k, res, clique},
      clique = igKColoringHeuristic[OptionValue["ForcedColoring"]][graph];
      k = Max[2, Length[clique]];
      While[(res = igKVertexColoring[graph, k, clique]) === {},
        k++
      ];
      First[res]
    ]


PackageExport["IGMinimumEdgeColoring"]
IGMinimumEdgeColoring::usage = "IGMinimumEdgeColoring[graph] finds a minimum edge colouring of graph.";

Options[IGMinimumEdgeColoring] = { "ForcedColoring" -> "MaxDegreeClique" };
SyntaxInformation[IGMinimumEdgeColoring] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};
IGMinimumEdgeColoring[graph_?igGraphQ, opt : OptionsPattern[]] :=
    catch@Module[{fc = OptionValue["ForcedColoring"]},
      If[ListQ[fc],
        fc = Check[EdgeIndex[graph, #]& /@ fc, throw[$Failed]]
      ];
      IGMinimumVertexColoring[LineGraph@IGUndirectedGraph[graph, "All"], "ForcedColoring" -> fc]
    ]


PackageExport["IGChromaticNumber"]
IGChromaticNumber::usage = "IGChromaticNumber[graph] returns the chromatic number of graph.";

SyntaxInformation[IGChromaticNumber] = {"ArgumentsPattern" -> {_}};
IGChromaticNumber[graph_?igGraphQ] :=
    catch@Max[check@IGMinimumVertexColoring[graph], 0]


PackageExport["IGChromaticIndex"]
IGChromaticIndex::usage = "IGChromaticIndex[graph] returns the chromatic index of graph.";

SyntaxInformation[IGChromaticIndex] = {"ArgumentsPattern" -> {_}};
IGChromaticIndex[graph_?igGraphQ] :=
    catch@Max[check@IGMinimumEdgeColoring[graph], 0]


PackageExport["IGVertexColoringQ"]
IGVertexColoringQ::usage = "IGVertexColoringQ[graph, coloring] checks whether neighbouring vertices all have differing colours.";

SyntaxInformation[IGVertexColoringQ] = {"ArgumentsPattern" -> {_, _}};
IGVertexColoringQ[graph_?igGraphQ, col_List] /; VertexCount[graph] == Length[col] :=
  And @@ UnsameQ @@@ Partition[
    col[[ Flatten@IGIndexEdgeList[graph] ]],
    2
  ]
expr : IGVertexColoringQ[graph_?igGraphQ, col_List] :=
    (Message[IGVertexColoring::inv, col, HoldForm@OutputForm[expr], "vertex coloring for this graph"];
     $Failed)
expr : IGVertexColoringQ[graph_?igGraphQ, asc_?AssociationQ] :=
    If[Sort@Keys[asc] === Sort@VertexList[graph],
      IGVertexColoringQ[graph, Lookup[asc, VertexList[graph]]]
      ,
      Message[IGVertexColoring::inv, asc, HoldForm@OutputForm[expr], "vertex coloring for this graph"];
      $Failed
    ]


PackageExport["IGCliqueCover"]
IGCliqueCover::usage = "IGCliqueCover[graph] finds a minimum clique cover of graph, i.e. a partitioning of its vertices into a smallest number of cliques.";

Options[IGCliqueCover] = { Method -> "Minimum" };
SyntaxInformation[IGCliqueCover] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};
IGCliqueCover[graph_ /; UndirectedGraphQ[graph] && IGTriangleFreeQ[graph]] :=
    With[{ies = List @@@ FindIndependentEdgeSet[graph]},
      Join[
        List /@ Complement[VertexList[graph], Join @@ ies],
        ies
      ]
    ]
IGCliqueCover[graph_?igGraphQ, OptionsPattern[]] :=
    catch@Module[{coloringFun},
      Switch[OptionValue[Method],
        Automatic | "Minimum", coloringFun = IGMinimumVertexColoring,
        "Heuristic", coloringFun = IGVertexColoring,
        _, Message[IGCliqueCover::invopt, OptionValue[Method], Method, "Minimum"]; coloringFun = IGMinimumVertexColoring
      ];
      communitiesFromMembership[graph, check@coloringFun@GraphComplement@IGUndirectedGraph[graph, "Mutual"]]
    ]


PackageExport["IGCliqueCoverNumber"]
IGCliqueCoverNumber::usage = "IGCliqueCoverNumber[graphs] finds the clique vertex cover number of graph.";

SyntaxInformation[IGCliqueCoverNumber] = {"ArgumentsPattern" -> {_}};
IGCliqueCoverNumber[graph_?igGraphQ] :=
    IGChromaticNumber@GraphComplement@IGUndirectedGraph[graph, "Mutual"]


PackageExport["IGPerfectQ"]
IGPerfectQ::usage = "IGPerfectQ[graph] tests is graph is perfect. The chromatic number of clique number is the same in every induced subgraph of a perfect graph.";

IGPerfectQ::undir = "The input graph must be undirected.";
SyntaxInformation[IGPerfectQ] = {"ArgumentsPattern" -> {_}};
IGPerfectQ[graph_?EmptyGraphQ] := True
IGPerfectQ[graph_?UndirectedGraphQ] :=
    catch@With[{g = IndexGraph@SimpleGraph[graph]}, (* IndexGraph is to work around the unreliability of Subgraph with arbitrary vertex names *)
      AllTrue[ConnectedComponents[g], check@igPerfectQ@Subgraph[g, #]&]
    ]
IGPerfectQ[graph_?GraphQ] := (Message[IGPerfectQ::undir]; False)
IGPerfectQ[_] := False

igPerfectQ[graph_] :=
    Block[{ig = igMakeFast[graph]},
      check@ig@"perfectQ"[]
    ]
