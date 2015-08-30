(* Mathematica Package  *)
(* Created by IntelliJ IDEA and wlplugin.halirutan.de *)

(* :Title: IGraph/M    *)
(* :Context: IGraphM`  *)
(* :Author: szhorvat   *)
(* :Date: 2015-08-28   *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: 10.0 *)
(* :Copyright: (c) 2015 Szabolcs Horvát *)
(* :Keywords: igraph, graphs, networks, LibraryLink *)
(* :Discussion: igraph interface for Mathematica, see http://igraph.org *)

BeginPackage["IGraphM`"]

Get["LTemplate`LTemplatePrivate`"]

IGVersion::usage = "IGVersion[]";

IGBetweenness::usage = "IGBetweenness[graph]";
IGEdgeBetweenness::usage = "IGEdgeBetweenness[graph]";
IGCloseness::usage = "IGCloseness[graph, normalized]";

IGRewire::usage = "IGRewire[graph, n, options]";
IGRewireEdges::usage = "IGRewireEdges[graph, p, options]";

IGDirectedAcyclicGraphQ::usage = "IGDirectedAcyclicGraphQ[graph]";
IGConnectedQ::usage = "IGConnectedQ[graph]";

IGIsomorphic::usage = "IGIsomorphic[graph1, graph2]";
IGSubisomorphic::usage = "IGSubisomorphic[graph, subgraph]";

IGTopologicalOrdering::usage = "IGTopologicalOrdering[graph] returns a permutation that sorts the vertices in topological order.";

IGFeedbackArcSet::usage = "IGFeedbackArcSet[graph]";

IGDyadCensus::usage = "IGDyadCensus[graph]";
IGTriadCensus::usage = "IGTriadCensus[graph]";

Begin["`Private`"]

(***** Package variables *****)

$packageDirectory = DirectoryName[$InputFileName];
$libraryDirectory = FileNameJoin[{$packageDirectory, "LibraryResources", $SystemID}];
$sourceDirectory  = FileNameJoin[{$packageDirectory, "LibraryResources", "Source"}];

(* Add to $LibraryPath in case package is not installed in Applications *)
If[Not@MemberQ[$LibraryPath, $libraryDirectory],
  AppendTo[$LibraryPath, $libraryDirectory]
]


template = LTemplate["IGraphM",
  {
    LClass["IGlobal",
      {
        LFun["init", {}, "Void"],
        LFun["seedRandom", {Integer}, "Void"],
        LFun["version", {}, "UTF8String"]
      }
    ],
    LClass["IG",
      {
        (* Create *)

        LFun["fromEdgeList", {{Real, 2} (* edges *), Integer (* vertex count *), True|False (* directed *)}, "Void"],

        (* Structure *)

        LFun["edgeCount", {}, Integer],
        LFun["vertexCount", {}, Integer],

        LFun["edgeList", {}, {Real, 2}],

        (* Testing *)

        LFun["directedQ", {}, True|False],
        LFun["dagQ", {}, True|False],
        LFun["simpleQ", {}, True|False],
        LFun["connectedQ", {}, True|False],

        (* Centrality *)

        LFun["betweenness", {}, {Real, 1}],
        LFun["edgeBetweenness", {}, {Real, 1}],
        LFun["closeness", {True|False (* normalized *)}, {Real, 1}],

        (* Randomize *)

        LFun["rewire", {Integer (* n_trials *), True|False (* loops *)}, "Void"],
        LFun["rewireEdges", {Real (* probability *), True|False (* loops *), True|False (* multiple *)}, "Void"],

        (* Isomorphism *)

        LFun["isomorphic", {LExpressionID["IG"]}, True|False],
        LFun["subisomorphic", {LExpressionID["IG"]}, True|False],

        (* Topological sorting and directed acylic graphs *)

        LFun["topologicalSorting", {}, {Real, 1}],
        LFun["feedbackArcSet", {True|False}, {Real, 1}],

        (* Motifs and subgraph counts *)

        LFun["dyadCensus", {}, {Integer, 1}],
        LFun["triadCensus", {}, {Real, 1}]
      }
    ]
  }
];


(***** Compilation and loading *****)

(* TODO: separate compilation into a dedicated build file *)

packageAbort[] := (End[]; EndPackage[]; Abort[])

recompileLibrary[] :=
  Block[{$CCompiler},
    If[Not@DirectoryQ[$libraryDirectory],
      CreateDirectory[$libraryDirectory]
    ];
    $CCompiler = {
      "Compiler" -> CCompilerDriver`GenericCCompiler`GenericCCompiler,
      "CompilerName" -> "g++-mp-5",
      "CompilerInstallation" -> "/opt/local/bin",
      "SystemCompileOptions" -> {"-std=c++11", "-m64", "-fPIC", "-O2", "-framework Foundation", "-framework \"mathlink\""}
    };
    SetDirectory[$sourceDirectory];
    CompileTemplate[template, {"IGlobal.cpp"},
      "IncludeDirectories" -> {"$HOME/local/include"},
      "LibraryDirectories" -> {"$HOME/local/lib"},
      "CompileOptions" -> {"-ligraph"},
      "ShellCommandFunction" -> Print, "ShellOutputFunction" -> Print,
      "TargetDirectory" -> $libraryDirectory
    ];
    ResetDirectory[]
  ]


If[LoadTemplate[template] === $Failed,
  Print[Style["Loading failed, trying to recompile ...", Red]];
  recompileLibrary[];
  If[LoadTemplate[template] === $Failed
    ,
    Print[Style["Cannot load or compile library. \[FreakedSmiley] Aborting.", Red]];
    packageAbort[]
    ,
    Print[Style["Successfully compiled and loaded the library. \[HappySmiley]", Red]];
  ]
]


(***** Initialize library *****)

igraphGlobal = Make["IGlobal"]; (* there should only be a single object of this type *)
igraphGlobal@"init"[] (* Run initialization *)


(***** Helper functions *****)

igEdgeList[g_?GraphQ] :=
    Developer`ToPackedArray@N[List @@@ EdgeList[g] /.
            Dispatch@Thread[VertexList[g] -> Range@VertexCount[g] - 1]]

igIndexVec[vec_?ArrayQ] := 1 + Round[vec]
igIndexVec[expr_] := expr

igDirectedQ[g_?GraphQ] := DirectedGraphQ[g] && Not@EmptyGraphQ[g]

igMake[g_?GraphQ] :=
    With[{ig = Make["IG"]},
      ig@"fromEdgeList"[igEdgeList[g], VertexCount[g], igDirectedQ[g]];
      ig
    ]

igToGraph[ig_] :=
    Graph[
      Range[ig@"vertexCount"[]],
      igIndexVec[ig@"edgeList"[]],
      DirectedEdges -> ig@"directedQ"[]
    ]


(***** Public functions *****)

IGVersion[] := igraphGlobal@"version"[]

IGBetweenness[g_?GraphQ] := Module[{ig = igMake[g]}, ig@"betweenness"[]]

IGEdgeBetweenness[g_?GraphQ] := Module[{ig = igMake[g]}, ig@"edgeBetwenness"[]]

Options[IGRewire] = { "AllowLoops" -> False };
IGRewire[g_?GraphQ, n_Integer, opt : OptionsPattern[]] :=
    Module[{ig = igMake[g]},
      ig@"rewire"[n, OptionValue["AllowLoops"]];
      igToGraph[ig]
    ]

Options[IGRewireEdges] = { "AllowLoops" -> False, "AllowMultipleEdges" -> False };
IGRewireEdges[g_?GraphQ, p_?Internal`RealValuedNumericQ, opt : OptionsPattern[]] :=
    Module[{ig = igMake[g]},
      ig@"rewireEdges"[p, OptionValue["AllowLoops"], OptionValue["AllowMultipleEdges"]];
      igToGraph[ig]
    ]

IGDirectedAcyclicGraphQ[g_?GraphQ] := Module[{ig = igMake[g]}, ig@"dagQ"[]]

IGConnectedQ[g_?GraphQ] := Module[{ig = igMake[g]}, ig@"connectedQ"[]]

IGCloseness[g_?GraphQ, normalized_ : False] := Module[{ig = igMake[g]}, ig@"closeness"[normalized]]

IGIsomorphic[g1_?GraphQ, g2_?GraphQ] := Block[{ig1 = igMake[g1], ig2 = igMake[g2]}, ig1@"isomorphic"[ManagedLibraryExpressionID@ig2]]

IGSubisomorphic[graph_?GraphQ, subgraph_?GraphQ] := Block[{ig1 = igMake[graph], ig2 = igMake[subgraph]}, ig1@"subisomorphic"[ManagedLibraryExpressionID@ig2]]

IGTopologicalOrdering[graph_?GraphQ] := Block[{ig = igMake[graph]}, igIndexVec@ig@"topologicalSorting"[]]

Options[IGFeedbackArcSet] = { "Exact" -> True };
IGFeedbackArcSet[graph_?GraphQ, opt : OptionsPattern[]] :=
  Block[{ig = igMake[graph]},
    Part[EdgeList[graph], igIndexVec@ig@"feedbackArcSet"[OptionValue["Exact"]]]
  ]

IGDyadCensus[graph_?GraphQ] := Block[{ig = igMake[graph]}, ig@"dyadCensus"[]]

IGTriadCensus[graph_?GraphQ] := Block[{ig = igMake[graph]}, Round[ig@"triadCensus"[]]]

End[] (* `Private` *)

EndPackage[]