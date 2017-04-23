(* Mathematica Package *)
(* Created by Mathematica Plugin for IntelliJ IDEA, see http://wlplugin.halirutan.de/ *)

(* :Title: IGraphM Utilities *)
(* :Context: IGraphM`Utilities` *)
(* :Author: szhorvat *)
(* :Date: 2016-06-11 *)

(* :Copyright: (c) 2017 Szabolcs Horvát *)

BeginPackage["IGraphM`Utilities`"];

Unprotect /@ Names["IGraphM`Utilities`*"];

IGNullGraphQ::usage = "IGNullGraphQ[graph] tests whether graph has no vertices.";

IGUndirectedGraph::usage = "IGUndirectedGraph[graph, conv] converts a directed graph to undirected with the given conversion method: \"Simple\" creates a single edge between connected vertices; \"All\" creates an undirected edge for each directed one and may produce a multigraph; \"Reciprocal\" creates a single undirected edge only between reciprocally connected vertices.";

IGReverseGraph::usage = "IGReverseGraph[graph] reverses the directed edges in graph while preserving edge weights.";

IGSimpleGraph::usage = "IGSimpleGraph[graph] converts graph to a simple graph by removing self loops and multi edges, according to the given options.";

IGUnweighted::usage = "IGUnweighted[graph] returns an unweighted version of an edge-weighted graph, while preserving other graph properties.";

IGWeightedAdjacencyGraph::usage =
    "IGWeightedAdjacencyGraph[matrix] creates a graph from a weighted adjacency matrix, taking the 0 weight to mean unconnected.\n" <>
    "IGWeightedAdjacencyGraph[vertices, matrix] uses vertices as the vertex names.\n" <>
    "IGWeightedAdjacencyGraph[matrix, z] creates a graph from a weighted adjacency matrix, taking the weight z to mean unconnected.\n" <>
    "IGWeightedAdjacencyGraph[vertices, matrix, z] uses vertices as the vertex names.";

IGVertexWeightedQ::usage = "IGVertexWeightedQ[graph] tests if graph is a vertex weighted graph.";
IGEdgeWeightedQ::usage = "IGEdgeWeightedQ[graph] tests if graph is an edge weighted graph.";

IGVertexProp::usage = "IGVertexProp[prop] is an operator that extracts the values of vertex property prop from a graph.";
IGEdgeProp::usage   = "IGEdgeProp[prop] is an operator that extracts the values of edge property prop from a graph.";

IGVertexMap::usage =
    "IGVertexMap[f, prop, graph] maps the function f to the vertex property list of property prop in graph.\n" <>
    "IGVertexMap[f, prop -> pf, graph] maps the function f to the values returned by pf[graph] and assigns the result to the vertex property prop.\n" <>
    "IGVertexMap[f, prop -> {pf1, pf2, \[Ellipsis]}, graph] threads f over {pf1[graph], pf2[graph], \[Ellipsis]} and assigns the result to the vertex property prop.\n" <>
    "IGVertexMap[f, spec] represents an operator form of IGVertexMap that can be applied to a graph.";
IGEdgeMap::usage =
    "IGEdgeMap[f, prop, graph] maps the function f to the edge property list of property prop in graph.\n" <>
    "IGEdgeMap[f, prop -> pf, graph] maps the function f to the values returned by pf[graph] and assigns the result to the edge property prop.\n" <>
    "IGEdgeMap[f, prop -> {pf1, pf2, \[Ellipsis]}, graph] threads f over {pf1[graph], pf2[graph], \[Ellipsis]} and assigns the result to the edge property prop.\n" <>
    "IGEdgeMap[f, spec] represents an operator form of IGEdgeMap that can be applied to a graph.";

Begin["`Private`"];

(* Common definitions *)
Get["IGraphM`Common`"];

IGNullGraphQ[g_?GraphQ] := VertexCount[g] === 0
IGNullGraphQ[_] = False;

SyntaxInformation[IGUndirectedGraph] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}, "OptionNames" -> optNames[Graph]};

IGUndirectedGraph[g_?UndirectedGraphQ, "Simple"|"All"|"Reciprocal", opt : OptionsPattern[Graph]] := Graph[g, opt]
IGUndirectedGraph[g_?igGraphQ, "Simple", opt : OptionsPattern[Graph]] := UndirectedGraph[g, opt]
IGUndirectedGraph[g_?igGraphQ, "All", opt : OptionsPattern[Graph]] := Graph[VertexList[g], UndirectedEdge @@@ EdgeList[g], opt]
IGUndirectedGraph[g_?igGraphQ, "Reciprocal", opt : OptionsPattern[Graph]] :=
    With[{am = AdjacencyMatrix[g]}, (* null graph has no adjacency matrix, but this case is caught by _?UndirectedGraphQ above *)
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


SyntaxInformation[IGUnweighted] = {"ArgumentsPattern" -> {_}};

IGUnweighted[g_?IGEdgeWeightedQ] := Graph[VertexList[g], EdgeList[g], FilterRules[Options[g], Except[EdgeWeight]]]
IGUnweighted[g_?GraphQ] := g


(* When am is a SparseArray, we need to canonicalize it and ensure that it has no explicit that is the same as the implicit one. *)
arrayRules[am_SparseArray, u_] := ArrayRules[SparseArray[am], u]
arrayRules[am_, u_] := ArrayRules[am, u]

SyntaxInformation[IGWeightedAdjacencyGraph] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}, "OptionNames" -> optNames[WeightedAdjacencyGraph]};

IGWeightedAdjacencyGraph[wam_?SquareMatrixQ, unconnected : Except[_?OptionQ] : 0, opt : OptionsPattern[WeightedAdjacencyGraph]] :=
      WeightedAdjacencyGraph[
        SparseArray[Most@arrayRules[wam, unconnected], Dimensions[wam], Infinity],
        opt
      ]

IGWeightedAdjacencyGraph[vertices_List, wam_?SquareMatrixQ, unconnected : Except[_?OptionQ] : 0, opt : OptionsPattern[WeightedAdjacencyGraph]] :=
    WeightedAdjacencyGraph[
      vertices,
      SparseArray[Most@arrayRules[wam, unconnected], Dimensions[wam], Infinity],
      opt
    ]


SyntaxInformation[IGVertexWeightedQ] = {"ArgumentsPattern" -> {_}};
IGVertexWeightedQ[g_] := WeightedGraphQ[g] && PropertyValue[g, VertexWeight] =!= Automatic

SyntaxInformation[IGEdgeWeightedQ] = {"ArgumentsPattern" -> {_}};
IGEdgeWeightedQ[g_] := WeightedGraphQ[g] && PropertyValue[g, EdgeWeight] =!= Automatic


gmiss = Missing["Nonexistent"];

SyntaxInformation[IGVertexProp] = {"ArgumentsPattern" -> {_}};
IGVertexProp[prop_][g_?GraphQ] := Replace[PropertyValue[{g,#}, prop]& /@ VertexList[g], $Failed -> gmiss, {1}]
IGVertexProp[prop : VertexWeight|VertexCapacity (* not VertexCoordinates! *)][g_?GraphQ] :=
    With[{values = PropertyValue[g, prop]},
      If[values === Automatic,
        ConstantArray[gmiss, VertexCount[g]],
        values
      ]
    ]

specialEdgeProps = EdgeWeight|EdgeCost|EdgeCapacity;

IGEdgeProp::nmg = "Multigraphs are only supported with the following properties: " <> ToString[List @@ specialEdgeProps] <> ".";

SyntaxInformation[IGEdgeProp] = {"ArgumentsPattern" -> {_}};
IGEdgeProp[prop : specialEdgeProps][g_?GraphQ] :=
    With[{values = PropertyValue[g, prop]},
      If[values === Automatic,
        ConstantArray[gmiss, EdgeCount[g]],
        values
      ]
    ]
IGEdgeProp[prop_][g_?GraphQ] :=
    Module[{multi = MultigraphQ[g]},
      If[multi, Message[IGEdgeProp::nmg]];
      Replace[PropertyValue[{g, #}, prop]& /@ EdgeList[g], $Failed -> gmiss, {1}] /; Not[multi]
    ]


igSetVertexProperty[g_, prop_, values_] /; Length[values] == VertexCount[g] :=
    SetProperty[g, Properties -> Thread[VertexList[g] -> List /@ Thread[prop -> values]]]
igSetVertexProperty[g_, prop : VertexWeight|VertexCapacity|VertexCoordinates, values_] /; Length[values] == VertexCount[g] :=
    Graph[VertexList[g], EdgeList[g], prop -> values, FilterRules[Options[g], Except[prop]]]
igSetVertexProperty[g_, prop_, values_] := $Failed


igSetEdgeProperty[g_, prop_, values_] /; Length[values] == EdgeCount[g] :=
    SetProperty[g, Properties -> Thread[EdgeList[g] -> List /@ Thread[prop -> values]]]
igSetEdgeProperty[g_, prop : specialEdgeProps, values_] /; Length[values] == EdgeCount[g] :=
    Graph[VertexList[g], EdgeList[g], prop -> values, FilterRules[Options[g], Except[prop]]]
igSetEdgeProperty[g_, prop_, values_] := $Failed


IGVertexMap::eprop = "Warning: `` is an edge property. It should not be used as a vertex property name.";
IGVertexMap::list = "`` did not return a list of appropriate length when applied to the graph.";
IGVertexMap::list2 = "The following functions did not return lists of approprate length when applied to the graph: ``.";

checkVertexProp[prop:
  EdgeWeight|EdgeCapacity|EdgeCost|
  EdgeStyle|EdgeShapeFunction|EdgeLabels|EdgeLabelStyle] := Message[IGVertexMap::eprop, prop]

SyntaxInformation[IGVertexMap] = {"ArgumentsPattern" -> {_, _, _.}};
IGVertexMap[fun_, prop_ -> pfun_, g_?GraphQ] :=
    Module[{values},
      checkVertexProp[prop];
      values = pfun[g];
      If[Not[ListQ[values] && Length[values] == VertexCount[g]],
        Message[IGVertexMap::list, pfun];
        Return[$Failed]
      ];
      igSetVertexProperty[g, prop, fun /@ pfun[g]]
    ]
IGVertexMap[fun_, prop_ -> pfunlist_List, g_?GraphQ] :=
    Module[{values, badpos},
      checkVertexProp[prop];
      values = Through[pfunlist[g]];
      badpos = Position[Length /@ values, Except@VertexCount[g], {1}, Heads -> False];
      If[badpos =!= {},
        Message[IGVertexMap::list2, Extract[pfunlist, badpos]];
        Return[$Failed]
      ];
      igSetVertexProperty[g, prop, MapThread[fun, Through[pfunlist[g]]]]
    ]
IGVertexMap[fun_, prop : Except[_Rule], g_?GraphQ] := IGVertexMap[fun, prop -> IGVertexProp[prop], g]
IGVertexMap[fun_, spec_][g_] := IGVertexMap[fun, spec, g]


IGEdgeMap::vprop = "Warning: `` is a vertex property. It should not be used as an edge property name.";
IGEdgeMap::list  = IGVertexMap::list;
IGEdgeMap::list2 = IGVertexMap::list2;
IGEdgeMap::nmg = IGEdgeProp::nmg;

checkEdgeProp[prop:
  VertexWeight|VertexCapacity|
  VertexSize|VertexShape|VertexShapeFunction|VertexStyle|VertexLabels|VertexLabelStyle|
  VertexCoordinates] := Message[IGEdgeMap::vprop, prop]

SyntaxInformation[IGEdgeMap] = {"ArgumentsPattern" -> {_, _, _.}};
IGEdgeMap[fun_, prop_ -> pfun_, g_?GraphQ] :=
    Module[{values},
      If[MultigraphQ[g] && Not@MatchQ[prop, specialEdgeProps],
        Message[IGEdgeMap::nmg];
        Return[$Failed]
      ];
      checkEdgeProp[prop];
      values = pfun[g];
      If[Not[ListQ[values] && Length[values] == EdgeCount[g]],
        Message[IGEdgeMap::list, pfun];
        Return[$Failed]
      ];
      igSetEdgeProperty[g, prop, fun /@ values]
    ]
IGEdgeMap[fun_, prop_ -> pfunlist_List, g_?GraphQ] :=
    Module[{values, badpos},
      If[MultigraphQ[g] && Not@MatchQ[prop, specialEdgeProps],
        Message[IGEdgeMap::nmg];
        Return[$Failed]
      ];
      checkEdgeProp[prop];
      values = Through[pfunlist[g]];
      badpos = Position[Length /@ values, Except@EdgeCount[g], {1}, Heads -> False];
      If[badpos =!= {},
        Message[IGEdgeMap::list2, Extract[pfunlist, badpos]];
        Return[$Failed]
      ];
      igSetEdgeProperty[g, prop, MapThread[fun, values]]
    ]
IGEdgeMap[fun_, prop : Except[_Rule], g_?GraphQ] := IGEdgeMap[fun, prop -> IGEdgeProp[prop], g]
IGEdgeMap[fun_, spec_][g_] := IGEdgeMap[fun, spec, g]


(***** Finalize *****)

(* Protect all package symbols *)
With[{syms = Names["IGraphM`Utilities`*"]}, SetAttributes[syms, {Protected, ReadProtected}] ];

End[]; (* `Private` *)

EndPackage[];