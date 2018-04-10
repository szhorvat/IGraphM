(* Mathematica Package *)
(* Created by Mathematica Plugin for IntelliJ IDEA, see http://wlplugin.halirutan.de/ *)

(* :Title: IGraphM Utilities *)
(* :Context: IGraphM`Utilities` *)
(* :Author: szhorvat *)
(* :Date: 2016-06-11 *)

(* :Copyright: (c) 2018 Szabolcs Horvát *)

BeginPackage["IGraphM`Utilities`"];

Unprotect /@ Names["IGraphM`Utilities`*"];

IGNullGraphQ::usage = "IGNullGraphQ[graph] tests whether graph has no vertices.";

IGUndirectedGraph::usage = "IGUndirectedGraph[graph, conv] converts a directed graph to undirected with the given conversion method: \"Simple\" creates a single edge between connected vertices; \"All\" creates an undirected edge for each directed one and may produce a multigraph; \"Reciprocal\" creates a single undirected edge only between reciprocally connected vertices.";

IGReverseGraph::usage = "IGReverseGraph[graph] reverses the directed edges in graph while preserving edge weights.";

IGSimpleGraph::usage = "IGSimpleGraph[graph] converts graph to a simple graph by removing self loops and multi edges, according to the given options.";

IGUnweighted::usage = "IGUnweighted[graph] returns an unweighted version of an edge-weighted graph, while preserving other graph properties.";

IGWeightedAdjacencyGraph::usage =
    "IGWeightedAdjacencyGraph[matrix] creates a graph from a weighted adjacency matrix, taking 0 to mean unconnected.\n" <>
    "IGWeightedAdjacencyGraph[vertices, matrix] uses vertices as the vertex names.\n" <>
    "IGWeightedAdjacencyGraph[matrix, z] creates a graph from a weighted adjacency matrix, taking the value z to mean unconnected.\n" <>
    "IGWeightedAdjacencyGraph[vertices, matrix, z] uses vertices as the vertex names.";

IGWeightedAdjacencyMatrix::usage =
    "IGWeightedAdjacencyMatrix[graph] gives the adjacency matrix of the edge weights of graph.\n" <>
    "IGWeightedAdjacencyMatrix[graph, z] gives the adjacency matrix of the edge weights of graph, using the value z to represent absent connections.";

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

IGVertexPropertyList::usage = "IGVertexPropertyList[graph] returns the list of available vertex properties in graph.";
IGEdgePropertyList::usage = "IGEdgePropertyList[graph] returns the list of available edge properties in graph.";

IGVertexStrength::usage =
    "IGVertexStrength[graph] returns the sum of edge weights for edges connecting to each vertex in graph.\n" <>
    "IGVertexStrength[graph, v] returns the sum of edge weights for edges connecting to vertex v in graph.";
IGVertexInStrength::usage =
    "IGVertexInStrength[graph] returns the sum of edge weights for the incoming edges of each vertex in graph.\n" <>
    "IGVertexInStrength[graph, v] returns the sum of edge weights for incoming edges of vertex v in graph.";
IGVertexOutStrength::usage =
    "IGVertexOutStrength[graph] returns the sum of edge weights for the outgoing edges of each vertex in graph.\n" <>
    "IGVertexOutStrength[graph, v] returns the sum of edge weights for outgoing edges of vertex v in graph.";

IGExportString::usage = "IGExportString[graph, format] generates a string corresponding to graph in the given format. See $IGExportFormats for supported formats.";

IGExport::usage =
    "IGExport[file, graph] exports graph to file in a format inferred from the file name.\n" <>
    "IGExport[file, graph, format] exports graph to file in the given format. See $IGExportFormats for supported formats.";

$IGExportFormats::usage = "$IGExportFormats is a list of export formats supported by IGExport.";

IGShorthand::usage = "IGShorthand[\"...\"] builds a graph from a shorthand notation such as \"a->b<-c\" or \"a-b,c-d\".";

IGPartitionsToMembership::usage =
    "IGPartitionsToMembership[elements, partitions] computes a membership vector for the given partitioning of elements.\n" <>
    "IGPartitionsToMembership[elements] is an operator that can be applied to partitions.";
IGMembershipToPartitions::usage =
    "IGMembershipToPartitions[elements, membership] computes a partitioning of elements based on the given membership vector." <>
    "IGMembershipToPartitions[elements] is an operator that can be applied to membership vectors.";

IGSourceVertexList::usage = "IGSourceVertexList[graph] returns the list of vertices with no incoming connections.";
IGSinkVertexList::usage = "IGSinkVertexList[graph] returns the list of vertices with no outgoing connections.";

IGReorderVertices::usage = "IGReorderVertices[vertices, graph] reorders the vertices of graph according to the given vertex vector. Graph properties are preserved.";

IGDirectedTree::usage = "IGDirectedTree[] is obsolete. Use IGOrientTree[] instead."; (* TODO: remove *)
IGOrientTree::usage =
    "IGOrientTree[tree, root] orients the edges of an undirected tree so that they point away from the root. The vertex order is not preserved: vertices will be ordered topologically.\n" <>
    "IGOrientTree[tree, root, \"In\"] orients the edges so that they point towards the root.";

IGTake::usage =
    "IGTake[graph, subgraph] keeps only those vertices and edges of graph which are also present in subgraph, while retaining all graph properties.\n" <>
    "IGTake[graph, edges] uses an edge list as the subgraph specification.";

IGZeroDiagonal::usage = "IGZeroDiagonal[matrix] replaces the diagonal of matrix with zeros.";

IGAdjacencyMatrixPlot::usage =
    "IGAdjacencyMatrixPlot[graph] plots the adjacency matrix of graph.\n" <>
    "IGAdjacencyMatrixPlot[graph, {v1, v2, \[Ellipsis]}] plots the adjacency matrix of the subgraph induced by the given vertices, using the specified vertex ordering.";

IGKirchhoffMatrix::usage =
    "IGKirchhoffMatrix[graph] returns the Kirchoff matrix, also known as Laplacian matrix of graph.\n" <>
    "IGKirchhoffMatrix[graph, \"In\"] will place the in-degrees on the diagonal instead of the out-degrees.";

IGGiantComponent::usage = "IGGiantComponent[graph] returns the largest weakly connected component of graph.";

Begin["`Private`"];

(* Common definitions *)
Get["IGraphM`Common`"];


(* Utility functions *)

If[$VersionNumber >= 10.1,
  keyValueMap = KeyValueMap,
  keyValueMap[f_, asc_] := f @@@ Normal[asc]
]


(* Main definitions *)

IGNullGraphQ[g_?GraphQ] := VertexCount[g] === 0
IGNullGraphQ[_] = False;


(* The built-in KirchhoffMatrix gives incorrect results for directed graphs. It uses the total degree
   on the diagonal. It should use the out-degree instead. *)
SyntaxInformation[IGKirchoffMatrix] = {"ArgumentsPattern" -> {_, _.}};
IGKirchhoffMatrix[graph_?GraphQ] := IGKirchhoffMatrix[graph, "Out"]
IGKirchhoffMatrix[graph_?GraphQ, "Out"] :=
    With[{am = zeroDiagonal@AdjacencyMatrix[graph]},
      DiagonalMatrix@SparseArray@Total[am, {2}] - am
    ]
IGKirchhoffMatrix[graph_?GraphQ, "In"] :=
    With[{am = zeroDiagonal@AdjacencyMatrix[graph]},
      DiagonalMatrix@SparseArray@Total[am] - am
    ]


SyntaxInformation[IGGiantComponent] = {"ArgumentsPattern" -> {_, OptionsPattern[]}, "OptionNames" -> optNames[Graph]};
IGGiantComponent[g_?IGNullGraphQ, opt : OptionsPattern[Graph]] := Graph[g, opt]
IGGiantComponent[g_?GraphQ, opt : OptionsPattern[Graph]] :=
    Subgraph[g, First@WeaklyConnectedComponents[g], opt]

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
IGReverseGraph[g_?UndirectedGraphQ, opt : OptionsPattern[]] := Graph[g, opt]
IGReverseGraph[g_?igGraphQ, opt : OptionsPattern[]] :=
    Module[{},
      If[MultigraphQ[g],
        Message[IGReverseGraph::nmg];
        Return[$Failed]
      ];
      Graph[
        VertexList[g],
        Reverse /@ EdgeList[g],
        opt,
        Options[g, {EdgeWeight, EdgeCapacity, EdgeCost, VertexWeight, VertexCapacity}]
      ]
    ]


Options[IGSimpleGraph] = { SelfLoops -> False, "MultipleEdges" -> False };
SyntaxInformation[IGSimpleGraph] = {"ArgumentsPattern" -> {_, OptionsPattern[]}, "OptionNames" -> optNames[IGSimpleGraph, Graph]};
IGSimpleGraph[g_?SimpleGraphQ, opt : OptionsPattern[]] := applyGraphOpt[opt][g]
IGSimpleGraph[g_?igGraphQ, opt : OptionsPattern[]] :=
    Module[{self = Not@TrueQ@OptionValue[SelfLoops], multi = Not@TrueQ@OptionValue["MultipleEdges"]},
      applyGraphOpt[opt]@Which[
        self && multi, SimpleGraph[g],
        self, removeSelfLoops[g],
        multi, removeMultiEdges[g],
        True, g
      ]
    ]


SyntaxInformation[IGUnweighted] = {"ArgumentsPattern" -> {_}};
IGUnweighted[g_?IGEdgeWeightedQ] := transformGraphOptions[ FilterRules[#, Except[EdgeWeight]]& ][g]
IGUnweighted[g_?GraphQ] := g
(* A more direct implementation is

     IGUnweighted[g_?IGEdgeWeightedQ] := Graph[ VertexList[g], EdgeList[g], FilterRules[Options[g], Except[EdgeWeight]] ]

   This implementation is slightly faster if the graph has a large number of properties, such as
   ExampleData[{"NetworkGraph", "CondensedMatterCollaborations"}]. However, it is noticeable slower for weighted graphs with
   packed weight vectors and no other properties.
*)


(* Weighted adjacency matrix handling *)

(*
(* When am is a SparseArray, we need to canonicalize it and ensure that it has no explicit value that is the same as the implicit one. *)
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
*)

SyntaxInformation[IGWeightedAdjacencyGraph] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}, "OptionNames" -> optNames[Graph]};
IGWeightedAdjacencyGraph[wam_?SquareMatrixQ, unconnected : Except[_?OptionQ] : 0, opt : OptionsPattern[Graph]] :=
    IGWeightedAdjacencyGraph[Range@Length[wam], wam, unconnected, opt]
IGWeightedAdjacencyGraph[vertices_List, wam_?SquareMatrixQ, unconnected : Except[_?OptionQ] : 0, opt : OptionsPattern[Graph]] :=
    Module[{sa = SparseArray[wam, Automatic, unconnected], directedEdges = OptionValue[DirectedEdges]},
      If[Length[vertices] != Length[sa],
        Message[IGWeightedAdjacencyGraph::ndims, vertices, wam];
        Return[$Failed]
      ];
      If[directedEdges === Automatic,
        directedEdges = Not@SymmetricMatrixQ[sa]
      ];
      If[Not[directedEdges],
        sa = UpperTriangularize[sa]
      ];
      Graph[vertices, sa["NonzeroPositions"], EdgeWeight -> sa["NonzeroValues"], DirectedEdges -> directedEdges, opt]
    ]


SyntaxInformation[IGWeightedAdjacencyMatrix] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};
Options[IGWeightedAdjacencyMatrix] = Options[WeightedAdjacencyMatrix];
IGWeightedAdjacencyMatrix[graph_?GraphQ, unconnected : Except[_?OptionQ] : 0, opt : OptionsPattern[]] :=
    With[{sa = WeightedAdjacencyMatrix[graph, opt]},
      SparseArray[sa["NonzeroPositions"] -> sa["NonzeroValues"], Dimensions[sa], unconnected]
    ]


(* Test for edge and vertex weightedness separately *)

(*
 * The original implementation was
 *     WeightedGraphQ[#] && PropertyValue[#, EdgeWeight] =!= Automatic &;
 * for testing edge-weightedness.
 * PropertyValue[g, EdgeWeight] fails on the null graph. This is why we had to test with WeightedGraphQ first.
 *
 * This solution was slow because PropertyValue takes a very long time.
 * (Benchmarked on ExampleData[{"NetworkGraph", "CondensedMatterCollaborations2005"}].)
 *
 * The alternative implementation
 *     MemberQ[PropertyList[graph], EdgeWeight]
 * does not return correct result for vertex-weighted but non-edge-weighted graphs
 *)

If[$VersionNumber >= 10.1, (* MinMax was added in M10.1 *)
  minMax = MinMax,
  minMax = {Min[#], Max[#]}&
];

SyntaxInformation[IGVertexWeightedQ] = {"ArgumentsPattern" -> {_}};
IGVertexWeightedQ[g_] :=
    WeightedGraphQ[g] &&
    With[{weights = Developer`ToPackedArray@GraphComputation`WeightVector[g]},
      If[First[weights] === 1 && minMax[weights] === {1, 1},
        PropertyValue[g, VertexWeight] =!= Automatic,
        True
      ]
    ]

SyntaxInformation[IGEdgeWeightedQ] = {"ArgumentsPattern" -> {_}};
IGEdgeWeightedQ[g_?EmptyGraphQ] := False (* avoid error from First if graph has no edges but is vertex weighted *)
IGEdgeWeightedQ[g_] :=
    WeightedGraphQ[g] &&
    With[{weights = Developer`ToPackedArray@GraphComputation`WeightValues[g]},
      If[First[weights] === 1 && minMax[weights] === {1, 1},
        PropertyValue[g, EdgeWeight] =!= Automatic,
        True
      ]
    ]


(* Vertex strengths, i.e. sums of the weights of incident edges *)

SyntaxInformation[IGVertexStrength] = {"ArgumentsPattern" -> {_, _.}};
IGVertexStrength[g_?IGNullGraphQ] := {}
IGVertexStrength[g_?igGraphQ] :=
    With[{am = WeightedAdjacencyMatrix[g]}, (* WeightedAdjacencyMatrix adds up weights in multigraphs. *)
      If[DirectedGraphQ[g],
        Total[am] + Total[am, {2}],
        Total[am] + Diagonal[am]
      ]
    ]
IGVertexStrength[g_?igGraphQ, v_] :=
    With[{index= VertexIndex[g, v]},
      With[{am = WeightedAdjacencyMatrix[g]},
        If[DirectedGraphQ[g],
          Total[am[[index]]] + Total[am[[;;, index]]],
          Total[am[[index]]] + am[[index, index]]
        ]
      ] /; IntegerQ[index]
    ]

SyntaxInformation[IGVertexInStrength] = {"ArgumentsPattern" -> {_, _.}};
IGVertexInStrength[g_?IGNullGraphQ] := {}
IGVertexInStrength[g_?igGraphQ] :=
    With[{am = WeightedAdjacencyMatrix[g]},
      If[DirectedGraphQ[g],
        Total[am],
        Total[am] + Diagonal[am]
      ]
    ]
IGVertexInStrength[g_?igGraphQ, v_] :=
    With[{index= VertexIndex[g, v]},
      With[{am = WeightedAdjacencyMatrix[g]},
        If[DirectedGraphQ[g],
          Total[am[[;;, index]]],
          Total[am[[index]]] + am[[index, index]]
        ]
      ] /; IntegerQ[index]
    ]

SyntaxInformation[IGVertexOutStrength] = {"ArgumentsPattern" -> {_, _.}};
IGVertexOutStrength[g_?IGNullGraphQ] := {}
IGVertexOutStrength[g_?igGraphQ] :=
    With[{am = WeightedAdjacencyMatrix[g]},
      If[DirectedGraphQ[g],
        Total[am, {2}],
        Total[am] + Diagonal[am]
      ]
    ]
IGVertexOutStrength[g_?igGraphQ, v_] :=
    With[{index= VertexIndex[g, v]},
      With[{am = WeightedAdjacencyMatrix[g]},
        If[DirectedGraphQ[g],
          Total[am[[index]]],
          Total[am[[index]]] + am[[index, index]]
        ]
      ] /; IntegerQ[index]
    ]

(***** Property transformation functions *****)

missing = Missing["Nonexistent"]; (* this will be used for non-existent property values *)

vertexWrapper::usage = "vertexWrapper is an internal symbolic wrapper used by IGVertexProp to work around retrieving properties from graphs with vertices that are lists.";

SyntaxInformation[IGVertexProp] = {"ArgumentsPattern" -> {_}};
IGVertexProp[prop_][g_?IGNullGraphQ] := {} (* some of the below may fail on null graphs, so we catch them early *)
IGVertexProp[VertexWeight][g_?GraphQ] :=
    If[IGVertexWeightedQ[g],
      GraphComputation`WeightVector[g],
      ConstantArray[missing, VertexCount[g]]
    ]
IGVertexProp[prop : (* VertexWeight| *)VertexCapacity (* not VertexCoordinates! *)][g_?GraphQ] :=
    With[{values = PropertyValue[g, prop]}, (* fails on null graph, but that is caught by the first pattern *)
      If[values === Automatic,
        ConstantArray[missing, VertexCount[g]],
        values
      ]
    ]
IGVertexProp[prop_][g_?GraphQ] :=
    (* work around PropertyValue failing when some graph vertices are lists *)
    With[{g2 = If[MemberQ[VertexList[g], _List], VertexReplace[g, v_ :> vertexWrapper[v]], g]},
      Replace[PropertyValue[{g2,#}, prop]& /@ VertexList[g2], $Failed -> missing, {1}]
    ]


specialEdgePropsPattern = EdgeWeight|EdgeCost|EdgeCapacity;

IGEdgeProp::nmg = "Multigraphs are only supported with the following properties: " <> ToString[List @@ specialEdgePropsPattern] <> ".";

SyntaxInformation[IGEdgeProp] = {"ArgumentsPattern" -> {_}};
IGEdgeProp[prop_][g_?IGNullGraphQ] := {} (* some of the below may fail on null graphs, so we catch them early *)
IGEdgeProp[EdgeWeight][g_?GraphQ] :=
    If[IGEdgeWeightedQ[g],
      GraphComputation`WeightValues[g],
      ConstantArray[missing, EdgeCount[g]]
    ]
IGEdgeProp[prop : specialEdgePropsPattern][g_?GraphQ] :=
    With[{values = PropertyValue[g, prop]}, (* fails on null graph, but that is caught by the first pattern *)
      If[values === Automatic,
        ConstantArray[missing, EdgeCount[g]],
        values
      ]
    ]
IGEdgeProp[prop_][g_?GraphQ] :=
    Module[{multi = MultigraphQ[g]},
      If[multi, Message[IGEdgeProp::nmg]];
      Replace[PropertyValue[{g, #}, prop]& /@ EdgeList[g], $Failed -> missing, {1}] /; Not[multi]
    ]


igSetVertexProperty[g_, prop_, values_] /; Length[values] == VertexCount[g] :=
    SetProperty[g, Properties -> Thread[VertexList[g] -> List /@ Thread[prop -> values]]]
igSetVertexProperty[g_, prop : VertexWeight|VertexCapacity|VertexCoordinates, values_] /; Length[values] == VertexCount[g] :=
    transformGraphOptions[ Append[FilterRules[#, Except[prop]], prop -> values]& ][g]
    (* Graph[VertexList[g], EdgeList[g], prop -> values, FilterRules[Options[g], Except[prop]]] *)
igSetVertexProperty[g_, prop_, values_] := $Failed


igSetEdgeProperty[g_, prop_, values_] /; Length[values] == EdgeCount[g] :=
    SetProperty[g, Properties -> Thread[EdgeList[g] -> List /@ Thread[prop -> values]]]
igSetEdgeProperty[g_, prop : specialEdgePropsPattern, values_] /; Length[values] == EdgeCount[g] :=
    transformGraphOptions[ Append[FilterRules[#, Except[prop]], prop -> values]& ][g]
    (* Graph[VertexList[g], EdgeList[g], prop -> values, FilterRules[Options[g], Except[prop]]] *)
igSetEdgeProperty[g_, prop_, values_] := $Failed


IGVertexMap::eprop    = "Warning: `` is an edge property. It should not be used as a vertex property name.";
IGVertexMap::propname = "Warning: Property name `` is not a symbol or string.";
IGVertexMap::list     = "`` did not return a list of appropriate length when applied to the graph.";
IGVertexMap::list2    = "The following functions did not return lists of appropriate length when applied to the graph: ``.";

checkVertexProp[prop:
  EdgeWeight|EdgeCapacity|EdgeCost|
  EdgeStyle|EdgeShapeFunction|EdgeLabels|EdgeLabelStyle] := Message[IGVertexMap::eprop, prop]
checkVertexProp[_Symbol | _String] := Null
checkVertexProp[prop_] := Message[IGVertexMap::propname, prop]

SyntaxInformation[IGVertexMap] = {"ArgumentsPattern" -> {_, _, _.}};
IGVertexMap[fun_, (Rule|RuleDelayed)[prop_, pfun_], g_?GraphQ] :=
    Module[{values},
      checkVertexProp[prop];
      values = pfun[g];
      If[Not[ListQ[values] && Length[values] == VertexCount[g]],
        Message[IGVertexMap::list, pfun];
        Return[$Failed]
      ];
      igSetVertexProperty[g, prop, fun /@ pfun[g]]
    ]
IGVertexMap[fun_, (Rule|RuleDelayed)[prop_, pfunlist_List], g_?GraphQ] :=
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
IGVertexMap[fun_, prop : Except[_Rule|_RuleDelayed], g_?GraphQ] := IGVertexMap[fun, prop -> IGVertexProp[prop], g]
IGVertexMap[fun_, spec_][g_] := IGVertexMap[fun, spec, g]


IGEdgeMap::vprop    = "Warning: `` is a vertex property. It should not be used as an edge property name.";
IGEdgeMap::propname = IGVertexMap::propname;
IGEdgeMap::list     = IGVertexMap::list;
IGEdgeMap::list2    = IGVertexMap::list2;
IGEdgeMap::nmg      = IGEdgeProp::nmg;

checkEdgeProp[prop:
  VertexWeight|VertexCapacity|
  VertexSize|VertexShape|VertexShapeFunction|VertexStyle|VertexLabels|VertexLabelStyle|
  VertexCoordinates] := Message[IGEdgeMap::vprop, prop]
checkEdgeProp[_Symbol | _String] := Null
checkEdgeProp[prop_] := Message[IGEdgeMap::propname, prop]

SyntaxInformation[IGEdgeMap] = {"ArgumentsPattern" -> {_, _, _.}};
IGEdgeMap[fun_, (Rule|RuleDelayed)[prop_, pfun_], g_?GraphQ] :=
    Module[{values},
      If[MultigraphQ[g] && Not@MatchQ[prop, specialEdgePropsPattern],
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
IGEdgeMap[fun_, (Rule|RuleDelayed)[prop_, pfunlist_List], g_?GraphQ] :=
    Module[{values, badpos},
      If[MultigraphQ[g] && Not@MatchQ[prop, specialEdgePropsPattern],
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
IGEdgeMap[fun_, prop : Except[_Rule|_RuleDelayed], g_?GraphQ] := IGEdgeMap[fun, prop -> IGEdgeProp[prop], g]
IGEdgeMap[fun_, spec_][g_] := IGEdgeMap[fun, spec, g]


(* Retrieve available edge and vertex property names *)

hasCustomProp[g_] := OptionValue[Options[g, Properties], Properties] =!= {}

standardVertexProperties = {
  VertexCoordinates,
  VertexShape, VertexShapeFunction, VertexSize, VertexStyle,
  VertexLabels, VertexLabelStyle,
  VertexWeight, VertexCapacity
};

SyntaxInformation[IGVertexPropertyList] = {"ArgumentsPattern" -> {_}};
IGVertexPropertyList[g_?IGNullGraphQ] = {};
IGVertexPropertyList[g_ /; GraphQ[g] && hasCustomProp[g]] := Sort@DeleteDuplicates[Join @@ PropertyList[{g, VertexList[g]}]]
IGVertexPropertyList[g_ /; GraphQ[g]] := Intersection[PropertyList[g], standardVertexProperties]

standardEdgeProperties = {
  EdgeStyle, EdgeShapeFunction, EdgeLabels, EdgeLabelStyle,
  EdgeWeight, EdgeCapacity, EdgeCost
};

SyntaxInformation[IGEdgePropertyList] = {"ArgumentsPattern" -> {_}};
IGEdgePropertyList[g_?EmptyGraphQ] = {};
IGEdgePropertyList[g_ /; GraphQ[g] && hasCustomProp[g]] := Sort@DeleteDuplicates[Join @@ PropertyList[{g, EdgeList[g]}]]
IGEdgePropertyList[g_ /; GraphQ[g]] := Intersection[PropertyList[g], standardEdgeProperties]


(***** Import and export *****)

$IGExportFormats = {"GraphML"};

IGExport::format   = "`` is not a recognized IGExport format.";
IGExport::infer    = "Cannot infer format of file \"``\".";
IGExport::mixed    = "Exporting mixed graphs to `` is not supported.";
IGExport::propname = "Only string or symbol property names are allowed when exporting to ``.";

igExportStringTag::usage = "igExportStringTag is a symbol used to signal use of ExportString instead of Export in IGExport.";

IGExportString[g_?GraphQ, format_]  := IGExport[igExportStringTag, g, format]
addCompletion[IGExportString, {0, $IGExportFormats}]

IGExport[file_, g_?GraphQ, format_] :=
    Switch[ToLowerCase[format], (* ToLowerCase doesn't error for non-Strings *)
      "graphml", ExportGraphML[file, g],
      _, Message[IGExport::format, format]; $Failed
    ]

IGExport[file_, g_?GraphQ] :=
    Switch[ToLowerCase@FileExtension[file],
      "graphml", IGExport[file, g, "GraphML"],
      _, Message[IGExport::infer, FileNameTake[file]]; $Failed
    ]

addCompletion[IGExport, {3, 0, $IGExportFormats}]


(* Converting to an exportable format, with properties *)

properties::usage =
    "properties[vertex, association] is a wrapped used to associate properties with a vertex when exporting GraphML." <>
    "properties[{v1, v2}, association] is used to associate properties with an edge when exporting GraphML.";

$ignoredEdgeProperties = {EdgeShapeFunction, EdgeStyle, EdgeLabelStyle};
$ignoredVertexProperties = {VertexShapeFunction, VertexShape, VertexStyle, VertexLabelStyle, VertexCoordinates, VertexSize};

(* Memoized *)
propType[propName_, g_, extractor_] := propType[propName, g, extractor] =
    With[{list = DeleteMissing[extractor[propName][g]]},
      Which[
        VectorQ[list, Developer`MachineIntegerQ], "Integer",
        VectorQ[list, Internal`RealValuedNumericQ], "Real",
        VectorQ[list, StringQ], "String",
        VectorQ[list, BooleanQ], "Boolean",
        True, "Expression"
      ]
    ]

vPropType[propName_, g_] := propType[propName, g, IGVertexProp]
ePropType[propName_, g_] := propType[propName, g, IGEdgeProp]

(* Memoized *)
vPropNames[g_] := vPropNames[g] = Complement[IGVertexPropertyList[g], $ignoredVertexProperties]
ePropNames[g_] := ePropNames[g] = Complement[IGEdgePropertyList[g], $ignoredEdgeProperties]

value[_][m_Missing] := m
value["Expression"][e_] := ToString[e, InputForm]
value["String"][e_] := e
value["Boolean"][e_] := ToString[e]
value["Integer"][e_] := ToString[e]
value["Real"][e_] := ToString@CForm@N[e]

getVertexProp[name_, g_] := value[vPropType[name, g]] /@ IGVertexProp[name][g]
getEdgeProp[name_, g_]   := value[ePropType[name, g]] /@ IGEdgeProp[name][g]

getObj[g_, getProp_, propNames_, list_] :=
    MapThread[
      properties,
      {
        list,
        With[{a = AssociationMap[getProp[#, g] &, propNames[g]]},
          If[a === <||>,
            ConstantArray[a, Length[list]],
            DeleteMissing /@ Transpose[a, AllowedHeads -> All]
          ]
        ]
      }
    ]

getVertices[g_] := getObj[g, getVertexProp, vPropNames, VertexList[g]]
getEdges[g_] := getObj[g, getEdgeProp, ePropNames, List @@@ EdgeList[g]]


(* GraphML *)

graphmlAttrType = <|"Real" -> "double", "Integer" -> "long", "String" -> "string", "Expression" -> "string", "Boolean" -> "boolean"|>;

graphmlAttrName[name_String] := graphmlAttrName[name] = name
graphmlAttrName[name_Symbol] := graphmlAttrName[name] = ToString[name]
graphmlAttrName[expr_] := (Message[IGExport::propname, "GraphML"]; Throw[$Failed, graphmlTag])

(* Memoized *)
graphmlVertexKeyID[name_] := graphmlVertexKeyID[name] = "v_" <> graphmlAttrName[name]
graphmlEdgeKeyID[name_] := graphmlEdgeKeyID[name] = "e_" <> graphmlAttrName[name]

graphmlEdgeKey[g_][name_] :=
    XMLElement["key",
      {"for" -> "edge", "id" -> graphmlEdgeKeyID[name],
        "attr.name" -> graphmlAttrName[name],
        "attr.type" -> graphmlAttrType@ePropType[name, g]}, {}
    ]

graphmlVertexKey[g_][name_] :=
    XMLElement["key",
      {"for" -> "node", "id" -> graphmlVertexKeyID[name],
        "attr.name" -> graphmlAttrName[name],
        "attr.type" -> graphmlAttrType@vPropType[name, g]}, {}
    ]

graphmlNode[properties[v_, asc_]] :=
    XMLElement["node",
      {"id" -> ToString[v]},
      graphmlData[asc, graphmlVertexKeyID]
    ]

graphmlEdge[properties[{v1_, v2_}, asc_]] :=
    XMLElement["edge",
      {"source" -> ToString[v1], "target" -> ToString[v2]},
      graphmlData[asc, graphmlEdgeKeyID]
    ]

graphmlData[asc_, keyID_] :=
    keyValueMap[
      Function[{name, value},
        XMLElement["data", {"key" -> keyID[name]}, {value}]
      ],
      asc
    ]

graphmlTemplate =
    XMLObject["Document"][
      {XMLObject["Declaration"]["Version" -> "1.0", "Encoding" -> "UTF-8"],
       XMLObject["Comment"][" created by IGraph/M, http://szhorvat.net/mathematica/IGraphM "]},
      XMLElement["graphml",
        {{"http://www.w3.org/2000/xmlns/", "xmlns"} ->
            "http://graphml.graphdrawing.org/xmlns",
          {"http://www.w3.org/2000/xmlns/", "xsi"} ->
              "http://www.w3.org/2001/XMLSchema-instance",
          {"http://www.w3.org/2001/XMLSchema-instance", "schemaLocation"} ->
              "http://graphml.graphdrawing.org/xmlns http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd"},
        {Sequence @@ #Keys,
          XMLElement["graph",
            {"id" -> "Graph", "edgedefault" -> #Directed},
            {Sequence @@ #Nodes, Sequence @@ #Edges}]}
      ], {}
    ] &;

graphmlGraph[g_?GraphQ] :=
    Internal`InheritedBlock[
      (* all memoized functions will be Block'ed *)
      {propType, vPropNames, ePropNames, graphmlEdgeKeyID, graphmlVertexKeyID, graphmlAttrName},
      If[MixedGraphQ[g],
        Message[IGExport::mixed, "GraphML"];
        Return[$Failed]
      ];
      Catch[
        graphmlTemplate[<|
          "Keys" -> Join[graphmlEdgeKey[g] /@ ePropNames[g], graphmlVertexKey[g] /@ vPropNames[g]],
          "Nodes" -> (graphmlNode /@ getVertices[g]),
          "Edges" -> (graphmlEdge /@ getEdges[g]),
          "Directed" -> If[igDirectedQ[g], "directed", "undirected"]
        |>],
        graphMLTag
      ]
    ]

ExportGraphML[file_, g_?GraphQ] :=
    With[{xml = graphmlGraph[g]},
      If[xml =!= $Failed,
        If[file === igExportStringTag,
          ExportString[xml, "XML"],
          Export[file, xml, "XML"]
        ],
        $Failed
      ]
    ]


(***** Shorthand *****)

(* Notes:

    - We cannot use removeSelfLoops[] and removeMultiEdges[] due to the need to support mixed graphs in this function.
    - The symbol sep should not be localized in IGShorthand because it is used in shSplitGroup[].
 *)

sep::usage = "sep[s] represents a separator while parsing shorthand notations in IGShorthand.";

SetAttributes[shToExpr, Listable];
shToExpr[s_String] :=
    Which[
      StringMatchQ[s, DigitCharacter ..], ToExpression[s],
      True, s
    ]

shSplitGroup[s_sep] := s
shSplitGroup[s_String] := shToExpr@StringTrim@StringSplit[s, ":"]

Options[IGShorthand] = {SelfLoops -> False, "MultipleEdges" -> False, DirectedEdges -> False};
SyntaxInformation[IGShorthand] = {"ArgumentsPattern" -> {_, OptionsPattern[]}, "OptionNames" -> optNames[IGShorthand, Graph]};
IGShorthand[s_String, opt : OptionsPattern[{IGShorthand, Graph}]] :=
    Module[{comp, seq, edgeSeq, eh, edges, vertices},
      comp = StringSplit[s, ","];
      seq = StringSplit[#, (x : "<" | "") ~~ ("-" ..) ~~ (y : "" | ">") :> sep[x <> y]] & /@ comp;
      seq = Replace[{__sep, mid___, __sep} :> {mid}] /@ seq;
      seq = Map[shSplitGroup, seq, {2}];
      edgeSeq = Select[seq, Length[#] > 1 &];
      eh = If[TrueQ@OptionValue[DirectedEdges], DirectedEdge, UndirectedEdge];
      edgeSeq = Replace[
        Catenate[Partition[#, 3, 2] & /@ edgeSeq],
        {
          {x_, sep[""], y_}   :> eh[x,y],
          {x_, sep[">"], y_}  :> DirectedEdge[x,y],
          {x_, sep["<"], y_}  :> DirectedEdge[y,x],
          {x_, sep["<>"], y_} :> Sequence[DirectedEdge[x,y], DirectedEdge[y,x]]
        },
        {1}
      ];
      edges = Catenate[Tuples /@ edgeSeq];
      Internal`InheritedBlock[{UndirectedEdge},
        SetAttributes[UndirectedEdge, Orderless];
        If[Not@TrueQ@OptionValue[SelfLoops],
          edges = DeleteCases[edges, _[x_, x_]];
        ];
        If[Not@TrueQ@OptionValue["MultipleEdges"],
          edges = DeleteDuplicates[edges];
        ]
      ];
      vertices = Catenate@DeleteCases[Level[seq, {2}], _sep];
      Graph[
        vertices,
        edges,
        Sequence@@FilterRules[{opt}, Options[Graph]],
        VertexLabels -> "Name"
      ]
    ]


(* Membership vectors *)

IGPartitionsToMembership::invpart = "Invalid element or part specification.";
SyntaxInformation[IGPartitionsToMembership] = {"ArgumentsPattern" -> {_, _.}};
IGPartitionsToMembership[elements_List, parts : {___List}] :=
    Module[{copy = parts},
      With[{union = Union @@ parts},
        If[Not[SubsetQ[elements, union] && union === Sort[Join @@ parts]],
          Message[IGPartitionsToMembership::invpart];
          Return[$Failed]
        ]
      ];
      copy[[All, All]] = Range@Length[parts];
      Lookup[AssociationThread[Flatten[parts, 1], Flatten[copy, 1]], elements, 0]
    ]
IGPartitionsToMembership[elements_][parts_] := IGPartitionsToMembership[elements, parts]

SyntaxInformation[IGMembershipToPartitions] = {"ArgumentsPattern" -> {_, _.}};
IGMembershipToPartitions[elements_List, membership_List] :=
    If[Length[elements] == Length[membership],
      Values@GroupBy[Transpose[{elements, membership}], Last -> First]
      ,
      Message[IGMembershipToPartitions::ndims, elements, membership];
      $Failed
    ]
IGMembershipToPartitions[elements_][membership_] := IGMembershipToPartitions[elements, membership]


(* Source and sink vertices *)

SyntaxInformation[IGSourceVertexList] = {"ArgumentsPattern" -> {_}};
IGSourceVertexList[g_?GraphQ] := Pick[VertexList[g], VertexInDegree[g], 0]

SyntaxInformation[IGSinkVertexList] = {"ArgumentsPattern" -> {_}};
IGSinkVertexList[g_?GraphQ] := Pick[VertexList[g], VertexOutDegree[g], 0]


(* Reorder vertices *)

findPermutation::usage = "findPermutation[l1, l2] finds the permutation that transforms list l1 into list l2.";
findPermutation[l1_, l2_] := Ordering[l1][[Ordering@Ordering[l2]]]

IGReorderVertices::bdvert = "The provided vertex list differs from the vertex list of the graph.";

SyntaxInformation[IGReorderVertices] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}, "OptionNames" -> optNames[Graph]};
IGReorderVertices[verts_List, graph_?GraphQ, opt : OptionsPattern[]] :=
    Module[{perm, vl = VertexList[graph]},
      If[Sort[verts] =!= Sort[vl],
        Message[IGReorderVertices::bdvert];
        Return[$Failed]
      ];
      perm = findPermutation[vl, verts];
      Graph[
        verts,
        EdgeList[graph],
        opt,
        Replace[
          Options[graph],
          Verbatim[Rule][sym : VertexWeight|VertexCapacity|VertexCoordinates, val_List] :> Rule[sym, val[[perm]]]
        ]
      ]
    ]


(* Direct edges of a tree *)

IGDirectedTree[tree_, root_, rest___] := IGOrientTree[tree, root, rest] (* old, obsolte naming *)

SyntaxInformation[IGOrientTree] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}, "OptionNames" -> optNames[Graph]};
expr : IGOrientTree[tree_?GraphQ, root_, mode : _String : "Out", opt : OptionsPattern[Graph]] :=
    Module[{verts},
      If[Not[UndirectedGraphQ[tree] && TreeGraphQ[tree]],
        Message[IGOrientTree::inv, HoldForm@OutputForm[expr], OutputForm[tree], "undirected tree"];
        Return[$Failed]
      ];
      If[Not@VertexQ[tree, root],
        Message[IGOrientTree::inv, HoldForm@OutputForm[expr], root, "vertex"];
        Return[$Failed]
      ];
      verts = First@Last@Reap@BreadthFirstScan[tree, root, "PrevisitVertex" -> Sow];
      Switch[mode,
        "Out", Null,
        "In", verts = Reverse[verts],
        _, Message[IGOrientTree::inv, HoldForm@OutputForm[expr], mode, "mode"]; Return[$Failed]
      ];
      DirectedGraph[IGReorderVertices[verts, tree], "Acyclic", opt]
    ]


(* Transfer properties to a smaller graph *)

(* Keep those elements in l1 which are also in l2. l1 may have repeating elements. *)
keepCases[l1_, l2_] := Pick[l1, Lookup[AssociationThread[l2, ConstantArray[1, Length[l2]]], l1, 0], 1]

allPropNames = {
  Properties,
  VertexWeight, VertexCapacity, VertexCoordinates,
  VertexSize, VertexShape, VertexShapeFunction, VertexStyle, VertexLabels, VertexLabelStyle,
  EdgeWeight, EdgeCapacity, EdgeCost,
  EdgeStyle, EdgeShapeFunction, EdgeLabels, EdgeLabelStyle
};

IGTake::nsg = "Some of the edges or vertices of the second graph are not present in the first.";

Options[IGTake] = Options[Graph];

SyntaxInformation[IGTake] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};

IGTake[g_?GraphQ, edges_List, opt : OptionsPattern[]] :=
    Module[{sg},
      sg = Quiet@Graph[edges, opt];
      IGTake[g, sg, opt] /; GraphQ[sg]
    ]

IGTake[g_?GraphQ, sg_?GraphQ, opt : OptionsPattern[]] :=
    Internal`InheritedBlock[{UndirectedEdge}, (* ensure that a <-> b compares equal to b <-> a *)
      SetAttributes[UndirectedEdge, Orderless];
      Module[{options, prop, vindex, eindex, sgEdgeList, vlist, elist, handleList, handleRule},
        sgEdgeList = DeleteDuplicates@EdgeList[sg];

        (* Check that sg is contained within g (ignores edge multiplicites). *)
        If[Not[SubsetQ[VertexList[g], VertexList[sg]] && SubsetQ[EdgeList[g], sgEdgeList]],
          Message[IGTake::nsg];
          Return[$Failed]
        ];

        (* These options contain the graph properties, but some options are for different purposes *)
        options = Association@Options[g];

        (* Used to find vertex/edge indices *)
        vindex = AssociationThread[VertexList[g], Range@VertexCount[g]];
        eindex = PositionIndex@EdgeList[g]; (* edge multiplicites must be handled *)

        vlist = VertexList[sg];
        elist = keepCases[EdgeList[g], sgEdgeList];

        (* Handle custom properties *)
        prop = If[KeyExistsQ[options, Properties],
          Properties ->
              Normal@KeyTake[
                options[Properties],
                Join[VertexList[sg], elist, {"DefaultEdgeProperties", "DefaultVertexProperties", "GraphProperties"}]
              ]
          ,
          Unevaluated@Sequence[]
        ];

        (* Handle List-type properties.  We preserve those with indices idx. *)
        handleList[idx_][propName_] :=
            If[KeyExistsQ[options, propName],
              propName -> options[propName][[idx]]
              ,
              Unevaluated@Sequence[]
            ];

        (* Handle Rule-type properties.  elems can be a list of vertices or edges. *)
        handleRule[elems_][propName_] :=
            If[KeyExistsQ[options, propName],
              Module[{rules = options[propName], default = {}},
                If[Not@ruleQ@First[rules],
                  default = {First[rules]};
                  rules = Rest[rules];
                ];
                propName -> Join[default, Normal@KeyTake[rules, elems]]
              ]
              ,
              Unevaluated@Sequence[]
            ];

        Graph[
          Graph[vlist, elist,
            prop,
            handleList[Lookup[vindex, VertexList[sg]]] /@ {VertexWeight, VertexCapacity, VertexCoordinates},
            handleList[Flatten@Lookup[eindex, elist]] /@ {EdgeWeight, EdgeCapacity, EdgeCost},
            handleRule[vlist] /@ {VertexSize, VertexShape, VertexShapeFunction, VertexStyle, VertexLabels, VertexLabelStyle},
            handleRule[elist] /@ {EdgeStyle, EdgeShapeFunction, EdgeLabels, EdgeLabelStyle},
            Normal@KeyDrop[options, allPropNames]
          ],
          opt] (* apply user options *)
      ]
    ]


(***** Matrix functions ****)

SyntaxInformation[IGZeroDiagonal] = {"ArgumentsPattern" -> {_}};
IGZeroDiagonal[mat_?MatrixQ] :=
    UpperTriangularize[mat, 1] + LowerTriangularize[mat, -1]



$unconnected::usage = "$unconnected is used to denote unconnected entries in the implementation of IGAdjacencyMatrixPlot.";

IGAdjacencyMatrixPlot::noprop = "The property `1` does not have a value for some or all edges in the graph.";
IGAdjacencyMatrixPlot::bdname = "The value of the VertexLabels option must be Automatic, \"Name\", \"Index\" or a list of rules.";

Options[IGAdjacencyMatrixPlot] = Options[MatrixPlot] ~Join~ {
  EdgeWeight -> Automatic, "UnconnectedColor" -> Automatic, VertexLabels -> Automatic, "RotateColumnLabels" -> True
};
SetOptions[IGAdjacencyMatrixPlot, Mesh -> Automatic];

SyntaxInformation[IGAdjacencyMatrixPlot] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

IGAdjacencyMatrixPlot[graph_?GraphQ, vs : (_List | All) : All, opt : OptionsPattern[]] :=
    Module[{$sizeLimit = 50, bigGraphQ, am, vind, rticks, cticks, colRules, rotFun,
            prop = OptionValue[EdgeWeight], mesh = OptionValue[Mesh], vertexLabels = OptionValue[VertexLabels], ticks = OptionValue[FrameTicks]},

      am = WeightedAdjacencyMatrix[graph, EdgeWeight -> prop];
      If[Not@MatrixQ[am],
        If[MemberQ[IGEdgeProp[prop][graph], _Missing],
          Message[IGAdjacencyMatrixPlot::noprop, prop]
        ];
        Return[$Failed]
      ];
      If[vs === All,
        vind = All,
        Check[vind = VertexIndex[graph, #]& /@ vs, Return[$Failed]];
      ];
      am = am[[vind, vind]];

      bigGraphQ = Length[am] > $sizeLimit;

      mesh = Replace[mesh, Automatic :> If[bigGraphQ, False, All]];
      rotFun = If[TrueQ@OptionValue["RotateColumnLabels"], Rotate[#, Pi/2]&, Identity];
      vertexLabels = Replace[vertexLabels, Automatic :> If[bigGraphQ, "Index", "Name"]];
      Switch[vertexLabels,
        "Index",
        {rticks, cticks} = {Automatic, Automatic};
        ,
        "Name",
        rticks = Transpose@{Range@Length[am], VertexList[graph][[vind]]};
        cticks = MapAt[rotFun, rticks, {All, 2}];
        ,
        {___?ruleQ},
        rticks = Transpose@{Range@Length[am], Replace[VertexList[graph], vertexLabels, {1}][[vind]]};
        cticks = MapAt[rotFun, rticks, {All, 2}];
        ,
        _,
        Message[IGAdjacencyMatrixPlot::bdname];
        {rticks, cticks} = {Automatic, Automatic};
      ];

      (* bring ticks to canonical form *)
      Switch[Dimensions[ticks, 2],
        {2, 2},
        Null;
        ,
        {2, _|PatternSequence[]},
        ticks = {#,#}& /@ ticks;
        ,
        _, (* catches single symbol, like Automatic, or single list *)
        ticks = {{#,#}, {#,#}}& @ ticks;
      ];
      ticks = Replace[ticks, {All|Automatic|True -> Automatic, None|False -> None}, {2}];

      ticks = MapAt[Replace[Automatic -> rticks], ticks, {1, All}];
      ticks = MapAt[Replace[Automatic -> cticks], ticks, {2, All}];

      (* construct ColorRules *)
      If[OptionValue["UnconnectedColor"] === Automatic,
        colRules = OptionValue[ColorRules]
        ,
        am = SparseArray[am["NonzeroPositions"] -> am["NonzeroValues"], Dimensions[am], $unconnected];
        colRules = {$unconnected -> OptionValue["UnconnectedColor"]};
        If[Not@MatchQ[OptionValue[ColorRules], Automatic|None],
          colRules = Flatten[{colRules, OptionValue[ColorRules]}]
        ]
      ];

      MatrixPlot[
        am,
        ColorRules -> colRules,
        FrameTicks -> ticks,
        Sequence @@ FilterRules[{opt}, FilterRules[Options[MatrixPlot], Except[ColorRules|FrameTicks]]],
        Mesh -> mesh
      ]
    ]


(***** Finalize *****)

(* Protect all package symbols *)
With[{syms = Names["IGraphM`Utilities`*"]}, SetAttributes[syms, {Protected, ReadProtected}] ];

End[]; (* `Private` *)

EndPackage[];