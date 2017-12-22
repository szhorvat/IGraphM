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
IGUnweighted[g_?IGEdgeWeightedQ] := Graph[VertexList[g], EdgeList[g], FilterRules[Options[g], Except[EdgeWeight]]]
IGUnweighted[g_?GraphQ] := g


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
        Total[am]
      ]
    ]
IGVertexStrength[g_?igGraphQ, v_] :=
    With[{index= VertexIndex[g, v]},
      With[{am = WeightedAdjacencyMatrix[g]},
        If[DirectedGraphQ[g],
          Total[am[[index]]] + Total[am[[;;, index]]],
          Total[am[[index]]]
        ]
      ] /; IntegerQ[index]
    ]

SyntaxInformation[IGVertexInStrength] = {"ArgumentsPattern" -> {_, _.}};
IGVertexInStrength[g_?IGNullGraphQ] := {}
IGVertexInStrength[g_?igGraphQ] := Total@WeightedAdjacencyMatrix[g]
IGVertexInStrength[g_?igGraphQ, v_] :=
    With[{index = VertexIndex[g, v]},
      With[{am = WeightedAdjacencyMatrix[g]},
        Total[am[[index]]]
      ] /; IntegerQ[index]
    ]

SyntaxInformation[IGVertexOutStrength] = {"ArgumentsPattern" -> {_, _.}};
IGVertexOutStrength[g_?IGNullGraphQ] := {}
IGVertexOutStrength[g_?igGraphQ] := Total[WeightedAdjacencyMatrix[g], {2}]
IGVertexOutStrength[g_?igGraphQ, v_] :=
    With[{index = VertexIndex[g, v]},
      With[{am = WeightedAdjacencyMatrix[g]},
        Total[am[[;;, index]]]
      ] /; IntegerQ[index]
    ]

(***** Property transformation functions *****)

missing = Missing["Nonexistent"]; (* this will be used for non-existent property values *)

SyntaxInformation[IGVertexProp] = {"ArgumentsPattern" -> {_}};
IGVertexProp[prop_][g_?IGNullGraphQ] := {} (* some of the below may fail on null graphs, so we catch them early *)
IGEdgeProp[VertexWeight][g_?GraphQ] :=
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
IGVertexProp[prop_][g_?GraphQ] := Replace[PropertyValue[{g,#}, prop]& /@ VertexList[g], $Failed -> missing, {1}]


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
    Graph[VertexList[g], EdgeList[g], prop -> values, FilterRules[Options[g], Except[prop]]]
igSetVertexProperty[g_, prop_, values_] := $Failed


igSetEdgeProperty[g_, prop_, values_] /; Length[values] == EdgeCount[g] :=
    SetProperty[g, Properties -> Thread[EdgeList[g] -> List /@ Thread[prop -> values]]]
igSetEdgeProperty[g_, prop : specialEdgePropsPattern, values_] /; Length[values] == EdgeCount[g] :=
    Graph[VertexList[g], EdgeList[g], prop -> values, FilterRules[Options[g], Except[prop]]]
igSetEdgeProperty[g_, prop_, values_] := $Failed


IGVertexMap::eprop    = "Warning: `` is an edge property. It should not be used as a vertex property name.";
IGVertexMap::propname = "Warning: Property name `` is not a symbol or string.";
IGVertexMap::list     = "`` did not return a list of appropriate length when applied to the graph.";
IGVertexMap::list2    = "The following functions did not return lists of approprate length when applied to the graph: ``.";

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
      If[Sort[elements] =!= Union @@ parts,
        Message[IGPartitionsToMembership::invpart];
        Return[$Failed]
      ];
      copy[[All, All]] = Range@Length[parts];
      Lookup[AssociationThread[Flatten[parts, 1], Flatten[copy, 1]], elements]
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

(***** Finalize *****)

(* Protect all package symbols *)
With[{syms = Names["IGraphM`Utilities`*"]}, SetAttributes[syms, {Protected, ReadProtected}] ];

End[]; (* `Private` *)

EndPackage[];