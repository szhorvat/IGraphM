(* Mathematica Package *)
(* Created by Mathematica plugin for IntelliJ IDEA *)

(* :Author: szhorvat *)
(* :Date: 2018-10-24 *)
(* :Copyright: (c) 2018 Szabolcs Horvát *)

Package["IGraphM`"]
igContextSetup[igPackagePrivateSymbol]

(*******************************************)
(***** Randomization and edge rewiring *****)
(*******************************************)


(* TODO: functions in this section should warn that edge weights will be lost *)

PackageExport["IGRewire"]
IGRewire::usage = "IGRewire[graph, n] attempts to rewire the edges of graph n times while preserving its degree sequence. Weights and other graph properties are discarded.";

IGRewire::multi = "The input is a multigraph. Multi-edges are never created during the rewiring process.";
Options[IGRewire] = { SelfLoops -> False };
SyntaxInformation[IGRewire] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}, "OptionNames" -> optNames[IGRewire, Graph]};
IGRewire[g_?igGraphQ, n_?Internal`NonNegativeMachineIntegerQ, opt : OptionsPattern[{IGRewire, Graph}]] :=
    catch@Block[{ig = igMakeFast[g]},
      If[MultigraphQ[g], Message[IGRewire::multi]];
      check@ig@"rewire"[n, OptionValue[SelfLoops]];
      applyGraphOpt[opt]@igToGraphWithNames[ig, VertexList[g]]
    ]


PackageExport["IGRewireEdges"]
IGRewireEdges::usage =
    "IGRewireEdges[graph, p] rewires each edge of the graph with probability p. Weights and other graph properties are discarded.\n" <>
    "IGRewireEdges[graph, p, \"In\"] rewires the starting point of each edge with probability p. The in-degree sequence is preserved.\n" <>
    "IGRewireEdges[graph, p, \"Out\"] rewires the endpoint of each edge with probability p. The out-degree sequence is preserved.";

Options[IGRewireEdges] = { SelfLoops -> False, MultiEdges -> False };
SyntaxInformation[IGRewireEdges] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}, "OptionNames" -> optNames[IGRewireEdges, Graph]};
IGRewireEdges[g_?igGraphQ, p_?Internal`RealValuedNumericQ, mode : All|"All"|"In"|"Out" : All, opt : OptionsPattern[{IGRewireEdges, Graph}]] :=
    catch@Block[{ig = igMakeFast[g]},
      Switch[mode,
        All|"All",
        check@ig@"rewireEdges"[p, OptionValue[SelfLoops], multiEdgesOptionReplace@OptionValue[MultiEdges]]
        ,
        "In",
        check@ig@"rewireDirectedEdges"[p, OptionValue[SelfLoops], False]
        ,
        "Out",
        check@ig@"rewireDirectedEdges"[p, OptionValue[SelfLoops], True]
      ];
      applyGraphOpt[opt]@igToGraphWithNames[ig, VertexList[g]]
    ]
addCompletion[IGRewireEdges, {0, 0, {"In", "Out", "All"}}]
