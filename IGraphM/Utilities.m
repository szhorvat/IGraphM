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

IGSimpleGraph::usage = "IGSimpleGraph[graph] converts graph to a simple graph by removing self loops and multi edges, according to the given options.";

IGWeightedAdjacencyGraph::usage =
    "IGWeightedAdjacencyGraph[matrix] creates a graph from a weighted adjacency matrix, taking the 0 weight to mean unconnected.\n" <>
    "IGWeightedAdjacencyGraph[vertices, matrix] uses vertices as the vertex names.\n" <>
    "IGWeightedAdjacencyGraph[matrix, z] creates a graph from a weighted adjacency matrix, taking the weight z to mean unconnected.\n" <>
    "IGWeightedAdjacencyGraph[vertices, matrix, z] uses vertices as the vertex names.";

IGVertexWeightedQ::usage = "IGVertexWeightedQ[graph] tests if graph is a vertex weighted graph.";
IGEdgeWeightedQ::usage = "IGEdgeWeightedQ[graph] tests if graph is an edge weighted graph.";


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


SyntaxInformation[IGReverseGraph] = {"ArgumentsPattern" -> {_}};

IGReverseGraph::nmg = "Multigraphs are not currently supported.";

IGReverseGraph[g_?UndirectedGraphQ] := g

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


SyntaxInformation[IGSimpleGraph] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

Options[IGSimpleGraph] = { SelfLoops -> False, "MultipleEdges" -> False };

IGSimpleGraph[g_?SimpleGraphQ, opt : OptionsPattern[]] := g
IGSimpleGraph[g_?igGraphQ, opt : OptionsPattern[]] :=
    Module[{self = Not@TrueQ@OptionValue[SelfLoops], multi = Not@TrueQ@OptionValue["MultipleEdges"]},
      Which[
        self && multi, SimpleGraph[g],
        self, removeSelfLoops[g],
        multi, removeMultiEdges[g],
        True, g
      ]
    ]


SyntaxInformation[IGWeightedAdjacencyGraph] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}, "OptionNames" -> optNames[WeightedAdjacencyGraph]};

IGWeightedAdjacencyGraph[wam_?SquareMatrixQ, unconnected : Except[_?OptionQ] : 0, opt : OptionsPattern[WeightedAdjacencyGraph]] :=
      WeightedAdjacencyGraph[
        SparseArray[Most@ArrayRules[wam, unconnected], Dimensions[wam], Infinity],
        opt
      ]

IGWeightedAdjacencyGraph[vertices_List, wam_?SquareMatrixQ, unconnected : Except[_?OptionQ] : 0, opt : OptionsPattern[WeightedAdjacencyGraph]] :=
    WeightedAdjacencyGraph[
      vertices,
      SparseArray[Most@ArrayRules[wam, unconnected], Dimensions[wam], Infinity],
      opt
    ]


SyntaxInformation[IGVertexWeightedQ] = {"ArgumentsPattern" -> {_}};
IGVertexWeightedQ[g_] := WeightedGraphQ[g] && PropertyValue[g, VertexWeight] =!= Automatic

SyntaxInformation[IGEdgeWeightedQ] = {"ArgumentsPattern" -> {_}};
IGEdgeWeightedQ[g_] := WeightedGraphQ[g] && PropertyValue[g, EdgeWeight] =!= Automatic


(***** Finalize *****)

(* Protect all package symbols *)
With[{syms = Names["IGraphM`Utilities`*"]}, SetAttributes[syms, {Protected, ReadProtected}] ];

End[]; (* `Private` *)

EndPackage[];