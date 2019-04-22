(* Mathematica Package *)
(* Created by Mathematica plugin for IntelliJ IDEA *)

(* :Author: szhorvat *)
(* :Date: 2018-10-24 *)
(* :Copyright: (c) 2019 Szabolcs Horvát *)

Package["IGraphM`"]
igContextSetup[igPackagePrivateSymbol]

(*******************************************)
(***** Shortest path related functions *****)
(*******************************************)


PackageExport["IGDistanceMatrix"]
IGDistanceMatrix::usage =
    "IGDistanceMatrix[graph] computes the shortest path length between each vertex pair in graph.\n" <>
    "IGDistanceMatrix[graph, fromVertices] computes the shortest path lengths between from the given vertices to each vertex in graph.\n" <>
    "IGDistanceMatrix[graph, fromVertices, toVertices] computes the shortest path lengths between the given vertices in graph.";

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
    "IGDistanceCounts[graph] computes a histogram of unweighted shortest path lengths between all vertex pairs. The kth element of the result is the count of shortest paths of length k. In undirected graphs, each path is counted only along one traversal direction.\n" <>
    "IGDistanceCounts[graph, fromVertices] computes a histogram of unweighted shortest path lengths from the given vertices to all others.";

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
    "IGNeighborhoodSize[graph, vertex] returns the number of direct neighbours of vertex.\n" <>
    "IGNeighborhoodSize[graph, {vertex1, vertex2, \[Ellipsis]}] returns the number of direct neighbours of each vertex.\n" <>
    "IGNeighborhoodSize[graph, {vertex1, vertex2, \[Ellipsis]}, max] returns the number of vertices reachable in at most max hops.\n" <>
    "IGNeighborhoodSize[graph, {vertex1, vertex2, \[Ellipsis]}, {n}] returns the number of vertices reachable in precisely n hops.\n" <>
    "IGNeighborhoodSize[graph, {vertex1, vertex2, \[Ellipsis]}, {min, max}] returns the number of vertices reachable in between min and max hops (inclusive).";

igNeighborhoodSize[graph_, vs_, {min_, max_}] :=
    Block[{ig = igMakeFast[graph]},
      Round@check@ig@"neighborhoodSize"[vss[graph][vs], min, max]
    ]

canonOrd[n_] := {0, n}
canonOrd[{n1_, n2_}] := {n1, n2}
canonOrd[{n_}] := {n, n}

ordQ[_?Internal`NonNegativeMachineIntegerQ |
    {_?Internal`NonNegativeMachineIntegerQ} |
    {_?Internal`NonNegativeMachineIntegerQ, _?Internal`NonNegativeMachineIntegerQ}
  ] := True
ordQ[_] := False

SyntaxInformation[IGNeighborhoodSize] = {"ArgumentsPattern" -> {_, _, _.}};
IGNeighborhoodSize[graph_?igGraphQ, {}, ord : _?ordQ : {1}] := {}
IGNeighborhoodSize[graph_?igGraphQ, vs : (_List | All), ord : _?ordQ : {1}] :=
    catch@igNeighborhoodSize[graph, vs, canonOrd[ord]]
IGNeighborhoodSize[graph_?igGraphQ, v_, ord : _?ordQ : {1}] :=
    catch@First@igNeighborhoodSize[graph, {v}, canonOrd[ord]]


PackageExport["IGDistanceHistogram"]
IGDistanceHistogram::usage =
    "IGDistanceHistogram[graph, binsize] computes a histogram of weighted all-pair shortest path lengths in graph with the given bin size. In the case of undirected graphs, path lengths are double counted.\n" <>
    "IGDistanceHistogram[graph, binsize, from] computes a histogram of weighted shortest path lengths in graph for the given starting vertices and bin size.\n" <>
    "IGDistanceHistogram[graph, binsize, from, to] computes a histogram of weighted shortest path lengths in graph for the given starting and ending vertices and bin size.";

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
IGAveragePathLength::usage = "IGAveragePathLength[graph] returns the average of all-pair unweighted shortest path lengths of the graph. Vertex pairs between which there is no path are excluded.";

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
      If[method === Automatic,
        method = Which[
          Not@IGEdgeWeightedQ[graph], "Unweighted",
          TrueQ[Min@igEdgeWeights[graph] >= 0], "Dijkstra",
          True, "Johnson"
        ]
      ];
      check@ig@"averagePathLengthWeighted"[ Lookup[igAveragePathLengthMethods, method, Message[IGAveragePathLength::bdmtd, method]; throw[$Failed]] ]
    ]


(***** Graph diameter *****)

PackageExport["IGDiameter"]
IGDiameter::usage = "IGDiameter[graph] computes the diameter of graph.";

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
