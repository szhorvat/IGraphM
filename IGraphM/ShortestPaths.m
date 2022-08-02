(* Mathematica Package *)
(* Created by Mathematica plugin for IntelliJ IDEA *)

(* :Author: szhorvat *)
(* :Date: 2018-10-24 *)
(* :Copyright: (c) 2018-2022 Szabolcs Horvát *)

Package["IGraphM`"]


(*******************************************)
(***** Shortest path related functions *****)
(*******************************************)


PackageExport["IGDistanceMatrix"]
IGDistanceMatrix::usage =
    "IGDistanceMatrix[graph] gives the shortest path length between each vertex pair in graph.\n" <>
    "IGDistanceMatrix[graph, fromVertices] gives the shortest path lengths between from the given vertices to each vertex in graph.\n" <>
    "IGDistanceMatrix[graph, fromVertices, toVertices] gives the shortest path lengths between the given vertices in graph.";

Options[IGDistanceMatrix] = {Method -> Automatic};
igDistanceMatrixMethods = <|
  "Unweighted" -> igDistanceMatrixUnweighted,
  "Dijkstra" -> igDistanceMatrixDijkstra,
  "BellmanFord" -> igDistanceMatrixBellmanFord,
  "Johnson" -> igDistanceMatrixJohnson
|>;

IGDistanceMatrix::bdmtd = "Value of option Method -> `` is not one of " <> ToString[Keys[igDistanceMatrixMethods], InputForm] <> ".";

SyntaxInformation[IGDistanceMatrix] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

amendUsage[IGDistanceMatrix, "Available Method options: <*Keys[igDistanceMatrixMethods]*>."];

IGDistanceMatrix[graph_?igGraphQ, from : (_List | All) : All, to : (_List | All) : All, opt : OptionsPattern[]] :=
    catch@Module[{method},
      method = OptionValue[Method];
      If[from === {}, Return[{}]];
      If[to === {}, Return[ConstantArray[{}, Length[from]]]];
      If[Not@MemberQ[Keys[igDistanceMatrixMethods] ~Join~ {Automatic}, method],
        Message[IGDistanceMatrix::bdmtd, method];
        Return[$Failed]
      ];
      If[method === Automatic,
        method = Which[
          Not@IGEdgeWeightedQ[graph], "Unweighted",
          TrueQ[Min@igEdgeWeights[graph] >= 0], "Dijkstra",
          True, "Johnson"
        ]
      ];
      igDistanceMatrixMethods[method][graph, vss[graph][from], vss[graph][to]]
    ]

igDistanceMatrixUnweighted[graph_, from_, to_] :=
    Block[{ig = igMakeFast[graph]},
      Round@expectInfNaN@fixInfNaN@check@ig@"shortestPaths"[from, to]
    ]

igDistanceMatrixDijkstra[graph_, from_, to_] :=
    Block[{ig = igMakeFastWeighted[graph]},
      expectInfNaN@fixInfNaN@check@ig@"shortestPathsDijkstra"[from, to]
    ]

igDistanceMatrixBellmanFord[graph_, from_, to_] :=
    Block[{ig = igMakeFastWeighted[graph]},
      expectInfNaN@fixInfNaN@check@ig@"shortestPathsBellmanFord"[from, to]
    ]

igDistanceMatrixJohnson[graph_, from_, to_] :=
    Block[{ig = igMakeFastWeighted[graph]},
      expectInfNaN@fixInfNaN@check@ig@"shortestPathsJohnson"[from, to]
    ]

(***** Shortest path tree *****)

PackageExport["IGShortestPathTree"]
IGShortestPathTree::usage = "IGShortestPathTree[graph, vertex] give the shortest path tree of graph rooted in vertex.";

Options[IGShortestPathTree] = {Method -> Automatic};
igShortestPathTreeMethods = <|
  "Unweighted" -> igSPTUnweighted,
  "Dijkstra" -> igSPTDijkstra,
  "BellmanFord" -> igSPTBellmanFord
|>;

SyntaxInformation[IGShortestPathTree] = {"ArgumentsPattern" -> {_, OptionsPattern[]}, "OptionNames" -> optNames[IGShortestPathTree, Graph]};

IGShortestPathTree::bdmtd = "Value of option Method -> `` is not one of " <> ToString[Keys[igShortestPathTreeMethods], InputForm] <> ".";

IGShortestPathTree[graph_?igGraphQ, from_, opt : OptionsPattern[{IGShortestPathTree, Graph}]] :=
    catch@Module[{method = OptionValue[Method]},
      If[Not@MemberQ[Keys[igDistanceMatrixMethods] ~Join~ {Automatic}, method],
        Message[IGShortestPathTree::bdmtd, method];
        Return[$Failed]
      ];
      If[method === Automatic,
        method = Which[
          Not@IGEdgeWeightedQ[graph], "Unweighted",
          TrueQ[Min@igEdgeWeights[graph] >= 0], "Dijkstra",
          True, "BellmanFord"
        ]
      ];
      applyGraphOpt[opt]@igShortestPathTreeMethods[method][graph, vs[graph][from]]
    ]

removeZeroes[vec_] := Delete[vec, Position[vec, 0]]

igSPTUnweighted[graph_, from_] :=
    Block[{ig = igMakeUnweighted[graph]},
      IGTakeSubgraph[
        graph,
        EdgeList[graph][[ removeZeroes@check@igIndexVec@ig@"shortestPathTreeEdges"[from] ]]
      ]
    ]

igSPTDijkstra[graph_, from_] :=
    Block[{ig = igMake[graph]},
      IGTakeSubgraph[
        graph,
        EdgeList[graph][[ removeZeroes@check@igIndexVec@ig@"shortestPathTreeEdgesDijkstra"[from] ]]
      ]
    ]

igSPTBellmanFord[graph_, from_] :=
    Block[{ig = igMake[graph]},
      IGTakeSubgraph[
        graph,
        EdgeList[graph][[ removeZeroes@check@igIndexVec@ig@"shortestPathTreeEdgesBellmanFord"[from] ]]
      ]
    ]

(***** Path length histograms and averages ****)

PackageExport["IGDistanceCounts"]
IGDistanceCounts::usage =
    "IGDistanceCounts[graph] gives a histogram of unweighted shortest path lengths between all vertex pairs. The kth element of the result is the count of shortest paths of length k. In undirected graphs, each path is counted only along one traversal direction.\n" <>
    "IGDistanceCounts[graph, fromVertices] gives a histogram of unweighted shortest path lengths from the given vertices to all others.";

SyntaxInformation[IGDistanceCounts] = {"ArgumentsPattern" -> {_, _.}};
IGDistanceCounts[graph_?igGraphQ] :=
    catch@Block[{ig = igMakeFast[graph]},
      Round@check@ig@"shortestPathCounts"[]
    ]
IGDistanceCounts[graph_?igGraphQ, {}] := {}
IGDistanceCounts[graph_?igGraphQ, vs : (_List | All)] :=
    catch@Block[{ig = igMakeFast[graph]},
      Round@check@ig@"shortestPathCounts2"[vss[graph][vs]]
    ]


PackageExport["IGNeighborhoodSize"]
IGNeighborhoodSize::usage =
    "IGNeighborhoodSize[graph, vertex] gives the number of direct neighbours of vertex, i.e. its degree.\n" <>
    "IGNeighborhoodSize[graph, All] gives the number of direct neighbours of all vertices.\n" <>
    "IGNeighborhoodSize[graph, {vertex1, vertex2, \[Ellipsis]}] gives the number of direct neighbours of the specified vertices.\n" <>
    "IGNeighborhoodSize[graph, All, max] gives the number of vertices reachable in at most max hops.\n" <>
    "IGNeighborhoodSize[graph, All, {n}] gives the number of vertices reachable in precisely n hops.\n" <>
    "IGNeighborhoodSize[graph, All, {min, max}] gives the number of vertices reachable in between min and max hops (inclusive).\n" <>
    "IGNeighborhoodSize[graph, All, {min, max}, mode] uses the given mode, \"In\", \"Out\" or \"All\", when finding neighbours in directed graphs.";

igNeighborhoodSize[graph_, vs_, {min_, max_}, mode_] :=
    Block[{ig = igMakeFast[graph]},
      Round@check@ig@"neighborhoodSize"[vss[graph][vs], min, max, encodeNeighborMode[mode]]
    ]

canonOrd[n_] := {0, n}
canonOrd[{n1_, n2_}] := {n1, n2}
canonOrd[{n_}] := {n, n}

ordQ[_?Internal`NonNegativeMachineIntegerQ |
    {_?Internal`NonNegativeMachineIntegerQ} |
    {_?Internal`NonNegativeMachineIntegerQ, _?Internal`NonNegativeMachineIntegerQ}
  ] := True
ordQ[_] := False

SyntaxInformation[IGNeighborhoodSize] = {"ArgumentsPattern" -> {_, _, _., _.}};
IGNeighborhoodSize[graph_?igGraphQ, {}, ord : _?ordQ : {1}, mode_String : "Out"] := {}
IGNeighborhoodSize[graph_?igGraphQ, vs : (_List | All), ord : _?ordQ : {1}, mode_String : "Out"] :=
    catch@igNeighborhoodSize[graph, vs, canonOrd[ord], mode]
IGNeighborhoodSize[graph_?igGraphQ, v_, ord : _?ordQ : {1}, mode_String : "Out"] :=
    catch@First@igNeighborhoodSize[graph, {v}, canonOrd[ord], mode]


PackageExport["IGDistanceHistogram"]
IGDistanceHistogram::usage =
    "IGDistanceHistogram[graph, binsize] gives a histogram of weighted all-pair shortest path lengths in graph with the given bin size. In the case of undirected graphs, path lengths are double counted.\n" <>
    "IGDistanceHistogram[graph, binsize, from] gives a histogram of weighted shortest path lengths in graph for the given starting vertices and bin size.\n" <>
    "IGDistanceHistogram[graph, binsize, from, to] gives a histogram of weighted shortest path lengths in graph for the given starting and ending vertices and bin size.";

Options[IGDistanceHistogram] = { Method -> "Dijkstra" };
SyntaxInformation[IGDistanceHistogram] = {"ArgumentsPattern" -> {_, _, _., _., OptionsPattern[]}};
igDistanceHistogramMethods = <| "Dijkstra" -> 0, "BellmanFord" -> 1 |>;
amendUsage[IGDistanceHistogram, "Available Method options: <*Keys[igDistanceHistogramMethods]*>."];
IGDistanceHistogram[graph_?igGraphQ, binsize_?positiveNumericQ, from : (_List | All) : All, to : (_List | All) : All, opt : OptionsPattern[]] :=
    catch@Block[{ig = igMakeFastWeighted[graph], fromidx, toidx},
      If[from === {} || to === {},
        Return[{}]
      ];
      If[from === All,
        fromidx = Range[0, VertexCount[graph]-1], (* Must not use {} for All. See C++ code. *)
        fromidx = vss[graph][from]
      ];
      toidx = vss[graph][to];

      check@ig@"shortestPathWeightedHistogram"[binsize, fromidx, toidx, Lookup[igDistanceHistogramMethods, OptionValue[Method], -1] ]
    ]


PackageExport["IGAveragePathLength"]
IGAveragePathLength::usage = "IGAveragePathLength[graph] returns the average of all-pair shortest path lengths of the graph. Vertex pairs between which there is no path are excluded.";

igAveragePathLengthMethods = <|
  "Unweighted" -> 0,
  "Dijkstra" -> 1,
  "BellmanFord" -> 2,
  "Johnson" -> 3
|>;

IGAveragePathLength::bdmtd = "Value of option Method -> `` is not one of " <> ToString[Keys[igAveragePathLengthMethods], InputForm] <> ".";

Options[IGAveragePathLength] = {
  Method -> Automatic,
  "ByComponents" -> True
};
SyntaxInformation[IGAveragePathLength] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};
amendUsage[IGAveragePathLength, "Available Method options: <*Keys[igAveragePathLengthMethods]*>."]
IGAveragePathLength[graph_?igGraphQ, opt : OptionsPattern[]] :=
    catch@Block[{ig = igMakeFastWeighted[graph], method = OptionValue[Method]},
      (* Automatic method selection:
          - non-weighted graph: "Unweighted"
          - weighted graph with non-negative weights: "Dijkstra"
          - otherwise (i.e. negative weights): "Johnson"
       *)
      If[method === Automatic,
        method = Which[
          Not@IGEdgeWeightedQ[graph], "Unweighted",
          TrueQ[Min@igEdgeWeights[graph] >= 0], "Dijkstra",
          True, "Johnson"
        ]
      ];
      check@ig@"averagePathLengthWeighted"[
        Lookup[igAveragePathLengthMethods, method, Message[IGAveragePathLength::bdmtd, method]; throw[$Failed]],
        OptionValue["ByComponents"]
      ]
    ]


(***** Efficiency measures *****)

PackageExport["IGGlobalEfficiency"]

IGGlobalEfficiency::usage = "IGGlobalEfficiency[graph] gives the global efficiency of graph.";
SyntaxInformation[IGGlobalEfficiency] = {"ArgumentsPattern" -> {_}};
IGGlobalEfficiency[graph_?igGraphQ] :=
    catch@Block[{ig = igMakeFastWeighted[graph]},
      check@ig@"globalEfficiency"[]
    ]


PackageExport["IGLocalEfficiency"]
IGLocalEfficiency::usage =
    "IGLocalEfficiency[graph] gives the local efficiency around each vertex of graph.\n" <>
    "IGLocalEfficiency[graph, {vertex1, vertex2, \[Ellipsis]}] gives the local efficiency around the given vertices.\n" <>
    "IGLocalEfficiency[graph, All, \"Out\"] uses outgoing edges to define the neighbourhood in a directed graph.";
Options[IGLocalEfficiency] = { DirectedEdges -> True };
SyntaxInformation[IGLocalEfficiency] = {"ArgumentsPattern" -> {_, _., _., OptionsPattern[]}};
IGLocalEfficiency[graph_?igGraphQ, {}, mode_String : "All", OptionsPattern[]] := {}
IGLocalEfficiency[graph_?igGraphQ, vs : (_List | All) : All, mode_String : "All", OptionsPattern[]] :=
    catch@Block[{ig = igMakeFastWeighted[graph]},
      (* Infinities may be returned when some of the edge weights are zero. *)
      expectInfNaN@fixInfNaN@check@ig@"localEfficiency"[vss[graph][vs], OptionValue[DirectedEdges], encodeNeighborMode[mode]]
    ]
addCompletion[IGLocalEfficiency, {0, 0, {"Out", "In", "All"}}]


PackageExport["IGAverageLocalEfficiency"]
IGAverageLocalEfficiency::usage =
    "IGAverageLocalEfficiency[graph] gives the average local efficiency of graph.\n" <>
    "IGAverageLocalEfficiency[graph, \"Out\"] uses outgoing edges to define the neighbourhood in a directed graph.";
Options[IGAverageLocalEfficiency] = { DirectedEdges -> True };
SyntaxInformation[IGAverageLocalEfficiency] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};
IGAverageLocalEfficiency[graph_?igGraphQ, mode_String : "All", OptionsPattern[]] :=
    catch@Block[{ig = igMakeFastWeighted[graph]},
      check@ig@"averageLocalEfficiency"[OptionValue[DirectedEdges], encodeNeighborMode[mode]]
    ]
addCompletion[IGAverageLocalEfficiency, {0, {"Out", "In", "All"}}]


(***** Graph diameter *****)

PackageExport["IGDiameter"]
IGDiameter::usage = "IGDiameter[graph] gives the diameter of graph.";

Options[IGDiameter] = { Method -> Automatic, "ByComponents" -> False };

igDiameterMethods = <|
  "Unweighted" -> igDiameterUnweighted,
  "Dijkstra" -> igDiameterDijkstra
|>;

SyntaxInformation[IGDiameter] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

amendUsage[IGDiameter, "Available Method options: <*Keys[igDiameterMethods]*>."];

IGDiameter::bdmtd = "Value of option Method -> `` is not one of " <> ToString[Keys[igDiameterMethods], InputForm] <> ".";

IGDiameter[graph_?igGraphQ, opt : OptionsPattern[]] :=
    Module[{method},
      method = OptionValue[Method];
      If[Not@MemberQ[Keys[igDiameterMethods] ~Join~ {Automatic}, method],
        Message[IGDiameter::bdmtd, method];
        Return[$Failed]
      ];
      If[method === Automatic,
        method = Which[
          IGEdgeWeightedQ[graph], "Dijkstra",
          True, "Unweighted"
        ]
      ];
      igDiameterMethods[method][graph, OptionValue["ByComponents"]]
    ]

igDiameterUnweighted[graph_, bycomp_] :=
    catch@Block[{ig = igMakeFast[graph]},
      Round@check@ig@"diameter"[bycomp]
    ]

igDiameterDijkstra[graph_, bycomp_] :=
    Block[{ig = igMakeFastWeighted[graph]},
      sck@ig@"diameterDijkstra"[bycomp]
    ]


PackageExport["IGFindDiameter"]
IGFindDiameter::usage = "IGFindDiameter[graph] returns a longest shortest path in graph, i.e. a shortest path with length equal to the graph diameter.";

Options[IGFindDiameter] = { Method -> Automatic, "ByComponents" -> False };

igFindDiameterMethods = <|
  "Unweighted" -> igFindDiameterUnweighted,
  "Dijkstra" -> igFindDiameterDijkstra
|>;

SyntaxInformation[IGFindDiameter] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

amendUsage[IGFindDiameter, "Available Method options: <*Keys[igFindDiameterMethods]*>."];

IGFindDiameter::bdmtd = "Value of option Method -> `` is not one of " <> ToString[Keys[igFindDiameterMethods], InputForm] <> ".";

IGFindDiameter[graph_?igGraphQ, opt : OptionsPattern[]] :=
    Module[{method},
      method = OptionValue[Method];
      If[Not@MemberQ[Keys[igFindDiameterMethods] ~Join~ {Automatic}, method],
        Message[IGFindDiameter::bdmtd, method];
        Return[$Failed]
      ];
      If[method === Automatic,
        method = Which[
          IGEdgeWeightedQ[graph], "Dijkstra",
          True, "Unweighted"
        ]
      ];
      igFindDiameterMethods[method][graph, OptionValue["ByComponents"]]
    ]

igFindDiameterUnweighted[graph_, bycomp_] :=
    catch@Block[{ig = igMakeFastWeighted[graph]},
      igVertexNames[graph]@igIndexVec@check@ig@"findDiameter"[bycomp]
    ]

igFindDiameterDijkstra[graph_, bycomp_] :=
    catch@Block[{ig = igMakeFastWeighted[graph]},
      igVertexNames[graph]@igIndexVec@check@ig@"findDiameterDijkstra"[bycomp]
    ]


(***** Pseudo-diameter *****)

PackageExport["IGPseudoDiameter"]
IGPseudoDiameter::usage = "IGPseudoDiameter[graph] gives the pseudodiameter of graph.";

Options[IGPseudoDiameter] = { Method -> Automatic, "ByComponents" -> False };

igPseudoDiameterMethods = <|
  "Unweighted" -> igPseudoDiameterUnweighted,
  "Dijkstra" -> igPseudoDiameterDijkstra
|>;

SyntaxInformation[IGPseudoDiameter] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

amendUsage[IGPseudoDiameter, "Available Method options: <*Keys[igPseudoDiameterMethods]*>."];

IGPseudoDiameter::bdmtd = "Value of option Method -> `` is not one of " <> ToString[Keys[igPseudoDiameterMethods], InputForm] <> ".";

IGPseudoDiameter[graph_?igGraphQ, v : ({}|Except[_?OptionQ]) : Automatic, opt : OptionsPattern[]] :=
    Module[{method},
      method = OptionValue[Method];
      If[Not@MemberQ[Keys[igPseudoDiameterMethods] ~Join~ {Automatic}, method],
        Message[IGPseudoDiameter::bdmtd, method];
        Return[$Failed]
      ];
      If[method === Automatic,
        method = Which[
          IGEdgeWeightedQ[graph], "Dijkstra",
          True, "Unweighted"
        ]
      ];
      igPseudoDiameterMethods[method][graph, If[v === Automatic, -1, vs[graph][v]], OptionValue["ByComponents"]]
    ]

igPseudoDiameterUnweighted[graph_, v_, bycomp_] :=
    Block[{ig = igMakeFast[graph]},
      sck@ig@"pseudoDiameter"[v, bycomp]
    ]

igPseudoDiameterDijkstra[graph_, v_, bycomp_] :=
    Block[{ig = igMakeFastWeighted[graph]},
      sck@ig@"pseudoDiameterDijkstra"[v, bycomp]
    ]


PackageExport["IGFindPseudoDiameter"]
IGFindPseudoDiameter::usage = "IGFindPseudoDiameter[graph] returns a longest shortest path in graph, i.e. a shortest path with length equal to the graph diameter.";

Options[IGFindPseudoDiameter] = { Method -> Automatic, "ByComponents" -> False };

igFindPseudoDiameterMethods = <|
  "Unweighted" -> igFindPseudoDiameterUnweighted,
  "Dijkstra" -> igFindPseudoDiameterDijkstra
|>;

SyntaxInformation[IGFindPseudoDiameter] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

amendUsage[IGFindPseudoDiameter, "Available Method options: <*Keys[igFindPseudoDiameterMethods]*>."];

IGFindPseudoDiameter::bdmtd = "Value of option Method -> `` is not one of " <> ToString[Keys[igFindPseudoDiameterMethods], InputForm] <> ".";

IGFindPseudoDiameter[graph_?igGraphQ, v : ({}|Except[_?OptionQ]) : Automatic, opt : OptionsPattern[]] :=
    Module[{method},
      method = OptionValue[Method];
      If[Not@MemberQ[Keys[igFindPseudoDiameterMethods] ~Join~ {Automatic}, method],
        Message[IGFindPseudoDiameter::bdmtd, method];
        Return[$Failed]
      ];
      If[method === Automatic,
        method = Which[
          IGEdgeWeightedQ[graph], "Dijkstra",
          True, "Unweighted"
        ]
      ];
      igFindPseudoDiameterMethods[method][graph, If[v === Automatic, -1, vs[graph][v]], OptionValue["ByComponents"]]
    ]

igFindPseudoDiameterUnweighted[graph_, v_, bycomp_] :=
    catch@Block[{ig = igMakeFastWeighted[graph]},
      igVertexNames[graph]@igIndexVec@check@ig@"findPseudoDiameter"[v, bycomp]
    ]

igFindPseudoDiameterDijkstra[graph_, v_, bycomp_] :=
    catch@Block[{ig = igMakeFastWeighted[graph]},
      igVertexNames[graph]@igIndexVec@check@ig@"findPseudoDiameterDijkstra"[v, bycomp]
    ]


(***** Graph girth *****)

PackageExport["IGGirth"]
IGGirth::usage = "IGGirth[graph] returns the length of the shortest cycle of the graph. The graph is treated as undirected, self-loops and multi-edges are ignored.";

SyntaxInformation[IGGirth] = {"ArgumentsPattern" -> {_}};
IGGirth[graph_?igGraphQ] :=
    catch@Block[{ig = igMakeFast[graph]},
      Replace[check@ig@"girth"[], 0 -> Infinity]
    ]


(***** Eccentricity *****)

PackageExport["IGEccentricity"]
IGEccentricity::usage =
    "IGEccentricity[graph] returns the eccentricity of all vertices.\n" <>
    "IGEccentricity[graph, vertex] returns the eccentricity of the given vertex.\n" <>
    "IGEccentricity[graph, {vertex1, vertex2, \[Ellipsis]}] returns the eccentricity of the given vertices.";

igEccentricity[graph_, vs_] :=
    Block[{ig = igMakeFast[graph]},
      Round@check@ig@"eccentricity"[vss[graph][vs]]
    ]

SyntaxInformation[IGEccentricity] = {"ArgumentsPattern" -> {_, _.}};
IGEccentricity[graph_?igGraphQ, {}] := {}
IGEccentricity[graph_?igGraphQ, vs : (_List | All) : All] := catch@igEccentricity[graph, vs]
IGEccentricity[graph_?igGraphQ, v_] := catch@First@igEccentricity[graph, {v}]


PackageExport["IGRadius"]
IGRadius::usage = "IGRadius[graph] returns the unweighted graph radius.";

SyntaxInformation[IGRadius] = {"ArgumentsPattern" -> {_}};
IGRadius[graph_?igGraphQ] :=
    catch@Block[{ig = igMakeFast[graph]},
      Round@check@ig@"radius"[]
    ]
