(* Mathematica Package *)
(* Created by Mathematica plugin for IntelliJ IDEA *)

(* :Author: szhorvat *)
(* :Date: 2018-10-24 *)
(* :Copyright: (c) 2018-2020 Szabolcs Horvát *)

Package["IGraphM`"]
igContextSetup[igPackagePrivateSymbol]

(*******************************)
(***** Centrality measures *****)
(*******************************)


PackageExport["IGBetweenness"]
IGBetweenness::usage =
    "IGBetweenness[graph] gives a list of betweenness centralities for the vertices of graph.\n" <>
    "IGBetweenness[graph, {vertex1, vertex2, \[Ellipsis]}] gives a list of betweenness centralities for the specified vertices.";

igBetweennessMethods = <| "Precise" -> False, "Fast" -> True |>;

Options[IGBetweenness] = { Method -> "Precise", Normalized -> False };

IGBetweenness::bdmtd =
    "Value of option Method -> `` is not one of " <>
        ToString[Keys[igBetweennessMethods], InputForm] <>
        ". Defaulting to " <> ToString[OptionValue[IGBetweenness, Method], InputForm] <> ".";

amendUsage[IGBetweenness, "Available Method options: <*Keys[igBetweennessMethods]*>."];

SyntaxInformation[IGBetweenness] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};

IGBetweenness[g_?igGraphQ, {}, opt : OptionsPattern[]] := {}
IGBetweenness[g_?igGraphQ, vs : (_List | All) : All,  opt : OptionsPattern[]] :=
    catch@Block[{ig = igMakeFastWeighted[g]},
      check@ig@"betweenness"[
        Lookup[igBetweennessMethods, OptionValue[Method], Message[IGBetweenness::bdmtd, OptionValue[Method]]; False],
        OptionValue[Normalized],
        vss[g][vs]
      ]
    ]


PackageExport["IGEdgeBetweenness"]
IGEdgeBetweenness::usage = "IGEdgeBetweenness[graph] gives a list of betweenness centralities for the edges of graph.";

(* Note: edge ordering is critical *)
Options[IGEdgeBetweenness] = { Normalized -> False };
SyntaxInformation[IGEdgeBetweenness] = {"ArgumentsPattern" -> {_}};
IGEdgeBetweenness[g_?igGraphQ, OptionsPattern[]] :=
    Block[{ig = igMake[g]}, sck@ig@"edgeBetweenness"[OptionValue[Normalized]]]


PackageExport["IGCloseness"]
IGCloseness::usage =
    "IGCloseness[graph] gives a list of closeness centralities for the vertices of graph.\n" <>
    "IGCloseness[graph, {vertex1, vertex2, \[Ellipsis]}] gives a list of closeness centralities for the specified vertices.";

Options[IGCloseness] = { Normalized -> False };
SyntaxInformation[IGCloseness] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};
IGCloseness[g_?igGraphQ, {}, opt : OptionsPattern[]] := {}
IGCloseness[g_?igGraphQ, vs : (_List | All) : All, opt : OptionsPattern[]] :=
    catch@Block[{ig = igMakeFastWeighted[g]},
      If[VertexCount[g] == 1, Developer`FromPackedArray, Identity] @ (* prevent {Indeterminate} packed array, which may misbehave, for single-vertex graph *)
        check@ig@"closeness"[OptionValue[Normalized], vss[g][vs]]
    ]


PackageExport["IGBetweennessEstimate"]
IGBetweennessEstimate::usage =
    "IGBetweennessEstimate[graph, cutoff] estimates vertex betweenness by considering only paths of at most length cutoff.\n" <>
    "IGBetweennessEstimate[graph, cutoff, {vertex1, vertex2, \[Ellipsis]}] estimates the betweenness of the specified vertices.";

Options[IGBetweennessEstimate] = { Method -> "Precise", Normalized -> False };
IGBetweennessEstimate::bdmtd = IGBetweenness::bdmtd;
amendUsage[IGBetweennessEstimate, "Available Method options: <*Keys[igBetweennessMethods]*>."];
SyntaxInformation[IGBetweennessEstimate] = {"ArgumentsPattern" -> {_, _, _., OptionsPattern[]}};
IGBetweennessEstimate[g_?igGraphQ, cutoff_?positiveOrInfQ, {}, opt : OptionsPattern[]] := {}
IGBetweennessEstimate[g_?igGraphQ, cutoff_?positiveOrInfQ, vs : (_List | All) : All, opt : OptionsPattern[]] :=
    catch@Block[{ig = igMakeFastWeighted[g]},
      check@ig@"betweennessEstimate"[
        infToZero[cutoff],
        Lookup[igBetweennessMethods, OptionValue[Method], Message[IGBetweennessEstimate::bdmtd, OptionValue[Method]]; False],
        OptionValue[Normalized],
        vss[g][vs]
      ]
    ]


PackageExport["IGEdgeBetweennessEstimate"]
IGEdgeBetweennessEstimate::usage = "IGEdgeBetweennessEstimate[graph, cutoff] estimates edge betweenness by considering only paths of at most length cutoff.";

(* Note: edge ordering is critical *)
Options[IGEdgeBetweennessEstimate] = { Normalized -> False };
SyntaxInformation[IGEdgeBetweennessEstimate] = {"ArgumentsPattern" -> {_, _}};
IGEdgeBetweennessEstimate[g_?igGraphQ, cutoff_?positiveOrInfQ, OptionsPattern[]] :=
    Block[{ig = igMake[g]}, sck@ig@"edgeBetweennessEstimate"[infToZero[cutoff], OptionValue[Normalized]]]


PackageExport["IGClosenessEstimate"]
IGClosenessEstimate::usage =
    "IGClosenessEstimate[graph, cutoff] estimates closeness centrality by considering only paths of at most length cutoff.\n" <>
    "IGClosenessEstimate[graph, cutoff, {vertex1, vertex2, \[Ellipsis]}] estimates the closeness centrality of the specified vertices.";

Options[IGClosenessEstimate] = { Normalized -> False };
SyntaxInformation[IGClosenessEstimate] = {"ArgumentsPattern" -> {_, _, _., OptionsPattern[]}};
IGClosenessEstimate[g_?igGraphQ, cutoff_?positiveOrInfQ, {}, opt : OptionsPattern[]] := {}
IGClosenessEstimate[g_?igGraphQ, cutoff_?positiveOrInfQ, vs : (_List | All) : All, opt : OptionsPattern[]] :=
    catch@Block[{ig = igMakeFastWeighted[g]},
      If[VertexCount[g] == 1, Developer`FromPackedArray, Identity] @ (* prevent {Indeterminate} packed array, which may misbehave, for single-vertex graph *)
        check@ig@"closenessEstimate"[infToZero[cutoff], OptionValue[Normalized], vss[g][vs]]
    ]


PackageExport["IGPageRank"]
IGPageRank::usage =
    "IGPageRank[graph] gives a list of PageRank centralities for the vertices of the graph using damping factor 0.85.\n" <>
    "IGPageRank[graph, damping] gives a list of PageRank centralities for the vertices of the graph using the given damping factor.";

igPageRankMethods = { "PowerIteration", "Arnoldi", "PRPACK" };
igPageRankMethodsAsc = AssociationThread[igPageRankMethods, Range@Length[igPageRankMethods] - 1];
igPageRankPowerOptions = <| "Epsilon" -> 0.001, "MaxIterations" -> 1000 |>;

Options[IGPageRank] = { Method -> "PRPACK", DirectedEdges -> True };
SyntaxInformation[IGPageRank] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};
amendUsage[IGPageRank, "Available Method options: <*igPageRankMethods*>."];

IGPageRank[graph_?GraphQ, damping : _?positiveNumericQ : 0.85, opt : OptionsPattern[]] :=
    catch@Block[{ig = igMakeFastWeighted[graph], method, methodOptions = {}, powerOpt},
      method = OptionValue[Method];
      If[ListQ[method],
        {method, methodOptions} = {First[method], Rest[method]};
      ];
      powerOpt = Join[igPageRankPowerOptions, Association[methodOptions]];
      check@ig@"pageRank"[
        Lookup[igPageRankMethodsAsc, method, -1], damping, OptionValue[DirectedEdges], powerOpt["MaxIterations"], powerOpt["Epsilon"]
      ]
    ]


PackageExport["IGPersonalizedPageRank"]
IGPersonalizedPageRank::usage =
    "IGPersonalizedPageRank[graph, reset] gives a list of personalized PageRank centralities for the vertices of the graph.\n" <>
    "IGPersonalizedPageRank[graph, reset, damping] gives a list of personalized PageRank centralities for the vertices of the graph using damping factor damping.";

Options[IGPersonalizedPageRank] = { Method -> "PRPACK", DirectedEdges -> True };
SyntaxInformation[IGPersonalizedPageRank] = {"ArgumentsPattern" -> {_, _, _., OptionsPattern[]}};

IGPersonalizedPageRank::invarg = "Second argument must be a vector of the same length as the vertex count of the graph.";

IGPersonalizedPageRank[graph_?GraphQ, reset_?VectorQ, damping : _?positiveNumericQ : 0.85, opt : OptionsPattern[]] :=
    catch@Block[{ig = igMakeFastWeighted[graph], method, methodOptions = {}, powerOpt},
      method = OptionValue[Method];
      If[ListQ[method],
        {method, methodOptions} = {First[method], Rest[method]};
      ];
      powerOpt = Join[igPageRankPowerOptions, Association[methodOptions]];
      If[Length[reset] != VertexCount[graph],
        Message[IGPersonalizedPageRank::invarg];
        throw[$Failed]
      ];
      check@ig@"personalizedPageRank"[
        Lookup[igPageRankMethodsAsc, method, -1], Normal[reset], damping, OptionValue[DirectedEdges], powerOpt["MaxIterations"], powerOpt["Epsilon"]
      ]
    ]


PackageExport["IGEigenvectorCentrality"]
IGEigenvectorCentrality::usage = "IGEigenvectorCentrality[graph] gives the eigenvector centrality of each vertex.";

Options[IGEigenvectorCentrality] = { DirectedEdges -> True, "Normalized" -> True };
SyntaxInformation[IGEigenvectorCentrality] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};
IGEigenvectorCentrality[graph_?GraphQ, opt : OptionsPattern[]] :=
    Block[{ig = igMakeFastWeighted[graph]},
      sck@ig@"eigenvectorCentrality"[OptionValue[DirectedEdges], OptionValue["Normalized"]]
    ]


PackageExport["IGHubScore"]
IGHubScore::usage = "IGHubScore[graph] gives Kleinberg's hub score for each vertex.";

Options[IGHubScore] = { "Normalized" -> True };
SyntaxInformation[IGHubScore] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};
IGHubScore[graph_?GraphQ, opt : OptionsPattern[]] :=
    Block[{ig = igMakeFastWeighted[graph]},
      sck@ig@"hubScore"[OptionValue["Normalized"]]
    ]


PackageExport["IGAuthorityScore"]
IGAuthorityScore::usage = "IGAuthorityScore[graph] gives Kleinberg's authority score for each vertex.";

Options[IGAuthorityScore] = { "Normalized" -> True };
SyntaxInformation[IGAuthorityScore] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};
IGAuthorityScore[graph_?GraphQ, opt : OptionsPattern[]] :=
    Block[{ig = igMakeFastWeighted[graph]},
      sck@ig@"authorityScore"[OptionValue["Normalized"]]
    ]


PackageExport["IGConstraintScore"]
IGConstraintScore::usage = "IGConstraintScore[graph] returns Burt's constraint score for each vertex.";

SyntaxInformation[IGConstraintScore] = {"ArgumentsPattern" -> {_}};
IGConstraintScore[graph_?igGraphQ] :=
    Block[{ig = igMakeFastWeighted[graph]},
      sck@ig@"constraintScore"[]
    ]


(***** Centralization *****)

PackageExport["IGDegreeCentralization"]
IGDegreeCentralization::usage =
    "IGDegreeCentralization[graph] gives the graph level centralization based on degree centralities.\n" <>
    "IGDegreeCentralization[graph, mode] uses the given mode, \"In\", \"Out\", or \"All\" to compute degrees in directed graphs.";

Options[IGDegreeCentralization] = { Normalized -> True, SelfLoops -> True };
SyntaxInformation[IGDegreeCentralization] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};
IGDegreeCentralization[graph_?igGraphQ, mode : _String : "All", opt : OptionsPattern[]] :=
    Block[{ig = igMakeFast[graph]},
      sck@ig@"degreeCentralization"[encodeNeighborMode[mode], OptionValue[SelfLoops], OptionValue[Normalized]]
    ]
addCompletion[IGDegreeCentralization, {0, {"In", "Out", "All"}}]


PackageExport["IGBetweennessCentralization"]
IGBetweennessCentralization::usage = "IGBetweennessCentralization[graph] gives the graph level centralization based on betweenness.";

IGBetweennessCentralization::bdmtd = IGBetweenness::bdmtd;
Options[IGBetweennessCentralization] = { Normalized -> True, Method -> "Precise" };
SyntaxInformation[IGBetweennessCentralization] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};
IGBetweennessCentralization[graph_?igGraphQ, opt : OptionsPattern[]] :=
    Block[{ig = igMakeFast[graph]},
      sck@ig@"betweennessCentralization"[
        Lookup[igBetweennessMethods, OptionValue[Method], Message[IGBetweennessCentralization::bdmtd, OptionValue[Method]]; False],
        OptionValue[Normalized]
      ]
    ]


PackageExport["IGClosenessCentralization"]
IGClosenessCentralization::usage = "IGClosenessCentralization[graph] gives the graph level centralization based on closeness.";

Options[IGClosenessCentralization] = { Normalized -> True };
SyntaxInformation[IGClosenessCentralization] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};
IGClosenessCentralization[graph_?igGraphQ, opt : OptionsPattern[]] :=
    Block[{ig = igMakeFast[graph]},
      sck@ig@"closenessCentralization"[OptionValue[Normalized]]
    ]


PackageExport["IGEigenvectorCentralization"]
IGEigenvectorCentralization::usage = "IGEigenvectorCentralization[graph] gives the graph level centralization based on eigenvector centralities.";

Options[IGEigenvectorCentralization] = { Normalized -> True, Scaled -> True };
SyntaxInformation[IGEigenvectorCentralization] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};
IGEigenvectorCentralization[graph_?igGraphQ, opt : OptionsPattern[]] :=
    Block[{ig = igMakeFast[graph]},
      sck@ig@"eigenvectorCentralization"[OptionValue[Scaled], OptionValue[Normalized]]
    ]
