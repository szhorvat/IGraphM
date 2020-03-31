(* Mathematica Package *)
(* Created by Mathematica plugin for IntelliJ IDEA *)

(* :Author: szhorvat *)
(* :Date: 2018-10-24 *)
(* :Copyright: (c) 2018-2020 Szabolcs Horvát *)

Package["IGraphM`"]
igContextSetup[igPackagePrivateSymbol]

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
      Round@fixInfNaN@check@ig@"shortestPaths"[from, to]
    ]

igDistanceMatrixDijkstra[graph_, from_, to_] :=
    Block[{ig = igMakeFastWeighted[graph]},
      fixInfNaN@check@ig@"shortestPathsDijkstra"[from, to]
    ]

igDistanceMatrixBellmanFord[graph_, from_, to_] :=
    Block[{ig = igMakeFastWeighted[graph]},
      fixInfNaN@check@ig@"shortestPathsBellmanFord"[from, to]
    ]

igDistanceMatrixJohnson[graph_, from_, to_] :=
    Block[{ig = igMakeFastWeighted[graph]},
      fixInfNaN@check@ig@"shortestPathsJohnson"[from, to]
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

Options[IGAveragePathLength] = { Method -> Automatic };
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
        Lookup[igAveragePathLengthMethods, method, Message[IGAveragePathLength::bdmtd, method]; throw[$Failed]]
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
IGLocalEfficiency[graph_?igGraphQ, {}, mode_String, OptionsPattern[]] := {}
IGLocalEfficiency[graph_?igGraphQ, vs : (_List | All) : All, mode_String : "All", OptionsPattern[]] :=
    catch@Block[{ig = igMakeFastWeighted[graph]},
      check@ig@"localEfficiency"[vss[graph][vs], OptionValue[DirectedEdges], encodeNeighborMode[mode]]
    ]
addCompletion[IGLocalEfficiency, {0, 0, {"Out", "In", "All"}}]


PackageExport["IGAverageLocalEfficiency"]
IGAverageLocalEfficiency::usage =
    "IGAverageLocalEfficiency[graph] gives the average local efficiency of graph.\n" <>
    "IGAverageLocalEfficiency[graph, \"Out\"] uses outgoing edges to define the neighbourhood in a directed graph.";
SyntaxInformation[IGAverageLocalEfficiency] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};
Options[IGAverageLocalEfficiency] = { DirectedEdges -> True };
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
    catch@Block[{ig = igMakeFast[graph], diam},
      diam = check@ig@"diameter"[bycomp];
      (* igraph returns the number of vertices for non-connected graphs.
         We translate this to Infinity, which is the only reasonable result.
         TODO https://github.com/igraph/igraph/issues/1345 *)
      If[diam == VertexCount[graph], Infinity, diam]
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

amendUsage[IGFindDiameter, "Available Method options: <*Keys[igDiameterMethods]*>."];

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


(***** Graph girth *****)

PackageExport["IGGirth"]
IGGirth::usage = "IGGirth[graph] returns the length of the shortest cycle of the graph. The graph is treated as undirected, self-loops and multi-edges are ignored.";

SyntaxInformation[IGGirth] = {"ArgumentsPattern" -> {_}};
IGGirth[graph_?igGraphQ] :=
    catch@Block[{ig = igMakeFast[graph]},
      Round@check@ig@"girth"[]
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
