(* Mathematica Package *)
(* Created by Mathematica Plugin for IntelliJ IDEA, see http://wlplugin.halirutan.de/ *)

(* :Title: IGraphM Utilities *)
(* :Context: IGraphM`Utilities` *)
(* :Author: szhorvat *)
(* :Date: 2016-06-11 *)

(* :Copyright: (c) 2017 Szabolcs Horvát *)

BeginPackage["IGraphM`Utilities`"];

Unprotect /@ Names["IGraphM`Utilities`*"];

IGUndirectedGraph::usage = "IGUndirectedGraph[graph, conv] converts a directed graph to undedirected with the given conversion method: \"Simple\" creates a single edge between connected vertices; \"All\" creates an undirected edge for each directed one and may produce a multigraph; \"Reciprocal\" creates a single undirected edge only between reciprocally connected vertices.";

IGReverseGraph::usage = "IGReverseGraph[graph] reverses the directed edges in graph while preserving edge weights.";

Begin["`Private`"];

(* Common definitions *)
Get["IGraphM`Common`"];

SyntaxInformation[IGUndirectedGraph] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}, "OptionNames" -> optNames[Graph]};

IGUndirectedGraph[g_?igGraphQ, "Simple", opt : OptionsPattern[Graph]] := UndirectedGraph[g, opt]
IGUndirectedGraph[g_?igGraphQ, "All", opt : OptionsPattern[Graph]] := Graph[VertexList[g], UndirectedEdge @@@ EdgeList[g], opt]
IGUndirectedGraph[g_?igGraphQ, "Reciprocal", opt : OptionsPattern[Graph]] :=
    With[{am = AdjacencyMatrix[g]},
      AdjacencyGraph[VertexList[g], Unitize[am Transpose[am]], opt]
    ]

IGUndirectedGraph[g_, "Mutual", opt : OptionsPattern[Graph]] := IGUndirectedGraph[g, "Reciprocal", opt]

IGUndirectedGraph[g_, "Each", opt : OptionsPattern[Graph]] := IGUndirectedGraph[g, "All", opt]
IGUndirectedGraph[g_, All, opt : OptionsPattern[Graph]] := IGUndirectedGraph[g, "All", opt]

IGUndirectedGraph[g_, "Collapse", opt : OptionsPattern[Graph]] := IGUndirectedGraph[g, "Simple", opt]
IGUndirectedGraph[g_, opt : OptionsPattern[Graph]] := IGUndirectedGraph[g, "Simple", opt]

addCompletion[IGUndirectedGraph, {0, {"Simple", "All", "Reciprocal"}}]


IGReverseGraph::nmg = "Multigraphs are not currently supported.";

IGReverseGraph[g_?igGraphQ] :=
    Module[{},
      If[MultigraphQ[g],
        Message[IGReverseGraph::nmg];
        Return[$Failed]
      ];
      Graph[
        VertexList[g],
        Reverse /@ EdgeList[g],
        Options[g, {EdgeWeight, EdgeCapacity, EdgeCost, VertexWeight, VertexCapacity}]
      ]
    ]

(***** Finalize *****)

(* Protect all package symbols *)
With[{syms = Names["IGraphM`Utilities`*"]}, SetAttributes[syms, {Protected, ReadProtected}] ];

End[]; (* `Private` *)

EndPackage[];