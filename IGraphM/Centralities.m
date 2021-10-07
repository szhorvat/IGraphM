(* Mathematica Package *)
(* Created by Mathematica plugin for IntelliJ IDEA *)

(* :Author: szhorvat *)
(* :Date: 2018-10-24 *)
(* :Copyright: (c) 2018-2020 Szabolcs Horvát *)

Package["IGraphM`"]


(*******************************)
(***** Centrality measures *****)
(*******************************)


PackageExport["IGBetweenness"]
IGBetweenness::usage =
    "IGBetweenness[graph] gives a list of betweenness centralities for the vertices of graph.\n" <>
    "IGBetweenness[graph, {vertex1, vertex2, \[Ellipsis]}] gives a list of betweenness centralities for the specified vertices.";

Options[IGBetweenness] = { Normalized -> False };
SyntaxInformation[IGBetweenness] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};
IGBetweenness[g_?igGraphQ, {}, opt : OptionsPattern[]] := {}
IGBetweenness[g_?igGraphQ, vs : (_List | All) : All,  opt : OptionsPattern[]] :=
    catch@Block[{ig = igMakeFastWeighted[g]},
      check@ig@"betweenness"[
        OptionValue[Normalized],
        vss[g][vs]
      ]
    ]


PackageExport["IGSubsetBetweenness"]
IGSubsetBetweenness::usage =
    "IGSubsetBetweenness[graph] gives a list of betweenness centralities for the vertices of graph.\n" <>
    "IGSubsetBetweenness[graph, {vertex1, vertex2, \[Ellipsis]}] gives a list of betweenness centralities for the specified vertices.";

Options[IGSubsetBetweenness] = { Normalized -> False };
SyntaxInformation[IGSubsetBetweenness] = {"ArgumentsPattern" -> {_, _, _, _., OptionsPattern[]}};
(*IGSubsetBetweenness[g_?igGraphQ, {}, opt : OptionsPattern[]] := {}*)
IGSubsetBetweenness[g_?igGraphQ, s_List, t_List, vs : (_List | All) : All,  opt : OptionsPattern[]] :=
    catch@Block[{ig = igMakeFastWeighted[g]},
      check@ig@"subsetBetweenness"[
        OptionValue[Normalized],
        vss[g][vs], vss[g][s], vss[g][t]
      ]
    ]


PackageExport["IGSubsetEdgeBetweenness"]
IGSubsetEdgeBetweenness::usage =
    "IGSubsetEdgeBetweenness[graph] gives a list of betweenness centralities for the vertices of graph.\n" <>
    "IGSubsetEdgeBetweenness[graph, {vertex1, vertex2, \[Ellipsis]}] gives a list of betweenness centralities for the specified vertices.";

Options[IGSubsetEdgeBetweenness] = { Normalized -> False };
SyntaxInformation[IGSubsetEdgeBetweenness] = {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};
IGSubsetEdgeBetweenness[g_?igGraphQ, s_List, t_List, opt : OptionsPattern[]] :=
    catch@Block[{ig = igMakeFastWeighted[g]},
      check@ig@"subsetEdgeBetweenness"[
        OptionValue[Normalized],
        vss[g][s], vss[g][t]
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

Options[IGCloseness] = { Normalized -> True };
SyntaxInformation[IGCloseness] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};
IGCloseness[g_?igGraphQ, {}, opt : OptionsPattern[]] := {}
IGCloseness[g_?igGraphQ, vs : (_List | All) : All, opt : OptionsPattern[]] :=
    catch@Block[{ig = igMakeFastWeighted[g]},
      expectInfNaN@fixInfNaN@check@ig@"closeness"[OptionValue[Normalized], vss[g][vs]]
    ]


PackageExport["IGHarmonicCentrality"]
IGHarmonicCentrality::usage =
    "IGHarmonicCentrality[graph] gives the harmonic centralities for the vertices of graph.\n" <>
    "IGHarmonicCentrality[graph, {vertex1, vertex2, \[Ellipsis]}] gives the harmonic centralities for the specified vertices.";

Options[IGHarmonicCentrality] = { Normalized -> True };
SyntaxInformation[IGHarmonicCentrality] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};
IGHarmonicCentrality[g_?igGraphQ, {}, opt : OptionsPattern[]] := {}
IGHarmonicCentrality[g_?igGraphQ, vs : (_List | All) : All, opt : OptionsPattern[]] :=
    catch@Block[{ig = igMakeFastWeighted[g]},
      check@ig@"harmonicCentrality"[OptionValue[Normalized], vss[g][vs]]
    ]


PackageExport["IGBetweennessCutoff"]
IGBetweennessCutoff::usage =
    "IGBetweennessCutoff[graph, cutoff] gives the range-limited betweenness centralities by considering only paths of at most length cutoff.\n" <>
    "IGBetweennessCutoff[graph, cutoff, {vertex1, vertex2, \[Ellipsis]}] gives the range-limited betweenness centralities for the specified vertices.";

Options[IGBetweennessCutoff] = { Normalized -> False };
SyntaxInformation[IGBetweennessCutoff] = {"ArgumentsPattern" -> {_, _, _., OptionsPattern[]}};
IGBetweennessCutoff[g_?igGraphQ, cutoff_?NonNegative, {}, opt : OptionsPattern[]] := {}
IGBetweennessCutoff[g_?igGraphQ, cutoff_?NonNegative, vs : (_List | All) : All, opt : OptionsPattern[]] :=
    catch@Block[{ig = igMakeFastWeighted[g]},
      check@ig@"betweennessCutoff"[
        infToNeg[cutoff],
        OptionValue[Normalized],
        vss[g][vs]
      ]
    ]

PackageExport["IGBetweennessEstimate"]
IGBetweennessEstimate::usage = "IGBetweennessEstimate[] is deprecated. Use IGBetweennessCutoff[] instead.";
IGBetweennessEstimate::deprec = "IGBetweennessEstimate is deprecated and will be removed from future versions of IGraph/M. Use IGBetweennessCutoff instead.";
IGBetweennessEstimate[args___] := (Message[IGBetweennessEstimate::deprec]; IGBetweennessCutoff[args]) (* TODO: remove eventually *)


PackageExport["IGEdgeBetweennessCutoff"]
IGEdgeBetweennessCutoff::usage = "IGEdgeBetweennessCutoff[graph, cutoff] gives the range-limited edge betweenness centralities by considering only paths of at most length cutoff.";

(* Note: edge ordering is critical *)
Options[IGEdgeBetweennessCutoff] = { Normalized -> False };
SyntaxInformation[IGEdgeBetweennessCutoff] = {"ArgumentsPattern" -> {_, _}};
IGEdgeBetweennessCutoff[g_?igGraphQ, cutoff_?NonNegative, OptionsPattern[]] :=
    Block[{ig = igMake[g]}, sck@ig@"edgeBetweennessCutoff"[infToNeg[cutoff], OptionValue[Normalized]]]

PackageExport["IGEdgeBetweennessEstimate"]
IGEdgeBetweennessEstimate::usage = "IGEdgeBetweennessEstimate[] is deprecated. Use IGEdgeBetweennessCutoff[] instead.";
IGEdgeBetweennessEstimate::deprec = "IGEdgeBetweennessEstimate is deprecated and will be removed from future versions of IGraph/M. Use IGEdgeBetweennessCutoff instead.";
IGEdgeBetweennessEstimate[args___] := (Message[IGEdgeBetweennessEstimate::deprec]; IGEdgeBetweennessCutoff[args]) (* TODO: remove eventually *)


PackageExport["IGClosenessCutoff"]
IGClosenessCutoff::usage =
    "IGClosenessCutoff[graph, cutoff] gives the range-limited closeness centralities by considering only paths of at most length cutoff.\n" <>
    "IGClosenessCutoff[graph, cutoff, {vertex1, vertex2, \[Ellipsis]}] gives the range-limited closeness centralities for the specified vertices.";

Options[IGClosenessCutoff] = { Normalized -> True };
SyntaxInformation[IGClosenessCutoff] = {"ArgumentsPattern" -> {_, _, _., OptionsPattern[]}};
IGClosenessCutoff[g_?igGraphQ, cutoff_?NonNegative, {}, opt : OptionsPattern[]] := {}
IGClosenessCutoff[g_?igGraphQ, cutoff_?NonNegative, vs : (_List | All) : All, opt : OptionsPattern[]] :=
    catch@Block[{ig = igMakeFastWeighted[g]},
      expectInfNaN@fixInfNaN@check@ig@"closenessCutoff"[infToNeg[cutoff], OptionValue[Normalized], vss[g][vs]]
    ]

PackageExport["IGClosenessEstimate"]
IGClosenessEstimate::usage = "IGClosenessEstimate[] is deprecated. Use IGClosenessCutoff[] instead.";
IGClosenessEstimate::deprec = "IGClosenessEstimate is deprecated and will be removed from future versions of IGraph/M. Use IGClosenessCutoff instead.";
IGClosenessEstimate[args___] := (Message[IGClosenessEstimate::deprec]; IGClosenessCutoff[args]) (* TODO: remove eventually *)


PackageExport["IGNeighborhoodCloseness"]
IGNeighborhoodCloseness::usage =
    "IGNeighborhoodCloseness[graph, cutoff] gives the range-limited closeness centralities along with the number of vertices reachable within the cutoff distance.\n" <>
    "IGNeighborhoodCloseness[graph, cutoff, {vertex1, vertex2, \[Ellipsis]}] gives the range-limited closeness centralities and number of reachable vertices for the specified vertices.";

Options[IGNeighborhoodCloseness] = { Normalized -> True };
SyntaxInformation[IGNeighborhoodCloseness] = {"ArgumentsPattern" -> {_, _, _., OptionsPattern[]}};
IGNeighborhoodCloseness[g_?igGraphQ, cutoff_?NonNegative, {}, opt : OptionsPattern[]] := {}
IGNeighborhoodCloseness[g_?igGraphQ, cutoff_?NonNegative, vs : (_List | All) : All, opt : OptionsPattern[]] :=
    catch@Block[{ig = igMakeFastWeighted[g]},
      expectInfNaN@fixInfNaN@check@ig@"neighborhoodCloseness"[infToNeg[cutoff], OptionValue[Normalized], vss[g][vs]]
    ]

PackageExport["IGHarmonicCentralityCutoff"]
IGHarmonicCentralityCutoff::usage =
    "IGHarmonicCentralityCutoff[graph, cutoff] gives the range-limited harmonic centralities by considering only paths of at most length cutoff.\n" <>
    "IGHarmonicCentralityCutoff[graph, cutoff, {vertex1, vertex2, \[Ellipsis]}] gives the range-limited harmonic centralities of the specified vertices.";

Options[IGHarmonicCentralityCutoff] = { Normalized -> True };
SyntaxInformation[IGHarmonicCentralityCutoff] = {"ArgumentsPattern" -> {_, _, _., OptionsPattern[]}};
IGHarmonicCentralityCutoff[g_?igGraphQ, cutoff_?NonNegative, {}, opt : OptionsPattern[]] := {}
IGHarmonicCentralityCutoff[g_?igGraphQ, cutoff_?NonNegative, vs : (_List | All) : All, opt : OptionsPattern[]] :=
    catch@Block[{ig = igMakeFastWeighted[g]},
      If[VertexCount[g] == 1, Developer`FromPackedArray, Identity] @ (* prevent {Indeterminate} packed array, which may misbehave, for single-vertex graph *)
        check@ig@"harmonicCentralityCutoff"[infToNeg[cutoff], OptionValue[Normalized], vss[g][vs]]
    ]


PackageExport["IGPageRank"]
IGPageRank::usage =
    "IGPageRank[graph] gives a list of PageRank centralities for the vertices of the graph using damping factor 0.85.\n" <>
    "IGPageRank[graph, damping] gives a list of PageRank centralities for the vertices of the graph using the given damping factor.";

igPageRankMethods = { "PowerIteration" (* TODO: no longer supported, remove and update C++ code *), "Arnoldi", "PRPACK" };
igPageRankMethodsAsc = AssociationThread[igPageRankMethods, Range@Length[igPageRankMethods] - 1];

Options[IGPageRank] = { Method -> "PRPACK", DirectedEdges -> True };
SyntaxInformation[IGPageRank] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};
amendUsage[IGPageRank, "Available Method options: <*DeleteCases[igPageRankMethods, \"PowerIteration\"]*>."];

IGPageRank[graph_?igGraphQ, damping : _?NumericQ : 0.85, opt : OptionsPattern[]] :=
    catch@Block[{ig = igMakeFastWeighted[graph], method, methodOptions = {}},
      method = OptionValue[Method];
      If[ListQ[method],
        {method, methodOptions} = {First[method], Rest[method]};
      ];
      check@ig@"personalizedPageRank"[
        Lookup[igPageRankMethodsAsc, method, -1], (* reset *) {}, damping, OptionValue[DirectedEdges]
      ]
    ]


PackageExport["IGLinkRank"]
IGLinkRank::usage =
    "IGLinkRank[graph] gives a list of LinkRank centralities for the edges of the graph using damping factor 0.85.\n" <>
    "IGLinkRank[graph, damping] gives a list of LinkRank centralities for the edges of the graph using the given damping factor.";

Options[IGLinkRank] = { Method -> "PRPACK", DirectedEdges -> True };
SyntaxInformation[IGLinkRank] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};
amendUsage[IGLinkRank, "Available Method options: <*DeleteCases[igPageRankMethods, \"PowerIteration\"]*>."];

IGLinkRank[graph_?igGraphQ, damping : _?NumericQ : 0.85, opt : OptionsPattern[]] :=
    catch@Block[{ig = igMakeFastWeighted[graph], method, methodOptions = {}},
      method = OptionValue[Method];
      If[ListQ[method],
        {method, methodOptions} = {First[method], Rest[method]};
      ];
      check@ig@"personalizedLinkRank"[
        Lookup[igPageRankMethodsAsc, method, -1], (* reset *) {}, damping, OptionValue[DirectedEdges]
      ]
    ]


PackageExport["IGPersonalizedPageRank"]
IGPersonalizedPageRank::usage =
    "IGPersonalizedPageRank[graph, reset] gives a list of personalized PageRank centralities for the vertices of the graph with personalization vector reset.\n" <>
    "IGPersonalizedPageRank[graph, reset, damping] uses the given damping factor.\n" <>
    "IGPersonalizedPageRank[graph, \[LeftAssociation] vertex1 -> weight1, vertex2 -> weight2, \[Ellipsis] \[RightAssociation], damping] uses non-zero personalization weights only for the specified vertices.";

Options[IGPersonalizedPageRank] = { Method -> "PRPACK", DirectedEdges -> True };
SyntaxInformation[IGPersonalizedPageRank] = {"ArgumentsPattern" -> {_, _, _., OptionsPattern[]}};

IGPersonalizedPageRank::invarg = "Second argument must be a vector of the same length as the vertex count of the graph.";

IGPersonalizedPageRank[graph_?igGraphQ, reset_?VectorQ, damping : _?NumericQ : 0.85, opt : OptionsPattern[]] :=
    catch@Block[{ig = igMakeFastWeighted[graph], method, methodOptions = {}},
      method = OptionValue[Method];
      If[ListQ[method],
        {method, methodOptions} = {First[method], Rest[method]};
      ];
      If[Length[reset] != VertexCount[graph],
        Message[IGPersonalizedPageRank::invarg];
        throw[$Failed]
      ];
      check@ig@"personalizedPageRank"[
        Lookup[igPageRankMethodsAsc, method, -1], Normal[reset], damping, OptionValue[DirectedEdges]
      ]
    ]

IGPersonalizedPageRank[graph_?igGraphQ, asc_?AssociationQ, damping : _?NumericQ : 0.85, opt : OptionsPattern[]] :=
    IGPersonalizedPageRank[graph, Lookup[asc, VertexList[graph], 0], damping, opt]


PackageExport["IGPersonalizedLinkRank"]
IGPersonalizedLinkRank::usage =
    "IGPersonalizedLinkRank[graph, reset] gives a list of personalized LinkRank centralities for the edges of the graph with personalization vector reset.\n" <>
    "IGPersonalizedLinkRank[graph, reset, damping] uses the given damping factor.\n" <>
    "IGPersonalizedLinkRank[graph, \[LeftAssociation] vertex1 -> weight1, vertex2 -> weight2, \[Ellipsis] \[RightAssociation], damping] uses non-zero personalization weights only for the specified vertices.";

Options[IGPersonalizedLinkRank] = { Method -> "PRPACK", DirectedEdges -> True };
SyntaxInformation[IGPersonalizedLinkRank] = {"ArgumentsPattern" -> {_, _, _., OptionsPattern[]}};

IGPersonalizedLinkRank::invarg = IGPersonalizedPageRank::invarg;

IGPersonalizedLinkRank[graph_?igGraphQ, reset_?VectorQ, damping : _?NumericQ : 0.85, opt : OptionsPattern[]] :=
    catch@Block[{ig = igMakeFastWeighted[graph], method, methodOptions = {}},
      method = OptionValue[Method];
      If[ListQ[method],
        {method, methodOptions} = {First[method], Rest[method]};
      ];
      If[Length[reset] != VertexCount[graph],
        Message[IGPersonalizedLinkRank::invarg];
        throw[$Failed]
      ];
      check@ig@"personalizedLinkRank"[
        Lookup[igPageRankMethodsAsc, method, -1], Normal[reset], damping, OptionValue[DirectedEdges]
      ]
    ]

IGPersonalizedLinkRank[graph_?igGraphQ, asc_?AssociationQ, damping : _?NumericQ : 0.85, opt : OptionsPattern[]] :=
    IGPersonalizedLinkRank[graph, Lookup[asc, VertexList[graph], 0], damping, opt]


PackageExport["IGEigenvectorCentrality"]
IGEigenvectorCentrality::usage = "IGEigenvectorCentrality[graph] gives the eigenvector centrality of each vertex.";

Options[IGEigenvectorCentrality] = { DirectedEdges -> True, "Normalized" -> True };
SyntaxInformation[IGEigenvectorCentrality] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};
IGEigenvectorCentrality[graph_?igGraphQ, opt : OptionsPattern[]] :=
    Block[{ig = igMakeFastWeighted[graph]},
      sck@ig@"eigenvectorCentrality"[OptionValue[DirectedEdges], OptionValue["Normalized"]]
    ]


PackageExport["IGHubScore"]
IGHubScore::usage = "IGHubScore[graph] gives Kleinberg's hub score for each vertex.";

Options[IGHubScore] = { "Normalized" -> True };
SyntaxInformation[IGHubScore] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};
IGHubScore[graph_?igGraphQ, opt : OptionsPattern[]] :=
    Block[{ig = igMakeFastWeighted[graph]},
      sck@ig@"hubScore"[OptionValue["Normalized"]]
    ]


PackageExport["IGAuthorityScore"]
IGAuthorityScore::usage = "IGAuthorityScore[graph] gives Kleinberg's authority score for each vertex.";

Options[IGAuthorityScore] = { "Normalized" -> True };
SyntaxInformation[IGAuthorityScore] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};
IGAuthorityScore[graph_?igGraphQ, opt : OptionsPattern[]] :=
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

Options[IGBetweennessCentralization] = { Normalized -> True };
SyntaxInformation[IGBetweennessCentralization] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};
IGBetweennessCentralization[graph_?igGraphQ, opt : OptionsPattern[]] :=
    Block[{ig = igMakeFast[graph]},
      sck@ig@"betweennessCentralization"[
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
