(* Mathematica Package  *)
(* Created by Mathematica Plugin for IntelliJ IDEA, see http://wlplugin.halirutan.de/ *)

(* :Author: szhorvat *)
(* :Date: 2018-10-23 *)
(* :Copyright: (c) 2019 Szabolcs Horvát *)

Package["IGraphM`"]
igContextSetup[igPackagePrivateSymbol]

(*********************************************)
(***** Property transformation functions *****)
(*********************************************)


(***** Edge and vertex property extractors *****)

missing = Missing["Nonexistent"]; (* this will be used for non-existent property values *)


PackageExport["IGVertexProp"]
IGVertexProp::usage = "IGVertexProp[prop] is an operator that extracts the values of vertex property prop from a graph.";

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


PackageExport["IGEdgeProp"]
IGEdgeProp::usage   = "IGEdgeProp[prop] is an operator that extracts the values of edge property prop from a graph.";

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


PackageExport["IGEdgeVertexProp"]
IGEdgeVertexProp::usage = "IGEdgeVertexProp[prop] is an operator that extracts the vertex property prop for the vertex pair corresponding to each edge.";

SyntaxInformation[IGEdgeVertexProp] = {"ArgumentsPattern" -> {_}};
IGEdgeVertexProp[prop_][g_?GraphQ] :=
    Partition[
      Part[
        IGVertexProp[prop][g],
        Flatten@IGIndexEdgeList[g]
      ],
      2
    ]


(***** Edge and vertex property mapping *****)

(* Check if parallel edges of a graph are distinguishable. *)
If[$VersionNumber >= 12.1,
  nonDistinguishableEdgesQ = MultigraphQ[#] && (Not@EdgeTaggedGraphQ[#] || Length@DeleteDuplicates@EdgeList[#] != EdgeCount[#])&
  ,
  nonDistinguishableEdgesQ = MultigraphQ
]


PackageScope["igSetVertexProperty"]
igSetVertexProperty::usage = "igSetVertexProperty[graph, prop, values]";
igSetVertexProperty[g_, prop_, values_] /; Length[values] == VertexCount[g] :=
    SetProperty[g, Properties -> Thread[VertexList[g] -> List /@ Thread[prop -> values]]]
igSetVertexProperty[g_, VertexShapeFunction, values_] /; Length[values] == VertexCount[g] :=
    SetProperty[g, VertexShapeFunction -> Thread[VertexList[g] -> values]] (* workaround for bug #87 *)
igSetVertexProperty[g_, prop : VertexWeight|VertexCapacity|VertexCoordinates, values_] /; Length[values] == VertexCount[g] :=
    transformGraphOptions[ Append[FilterRules[#, Except[prop]], prop -> values]& ][g]
(* Graph[VertexList[g], EdgeList[g], prop -> values, FilterRules[Options[g], Except[prop]]] *)
igSetVertexProperty[g_, prop_, values_] := $Failed


PackageScope["igSetEdgeProperty"]
igSetEdgeProperty::usage = "igSetEdgeProperty[graph, prop, values]";
igSetEdgeProperty[g_, prop_, values_] /; Length[values] == EdgeCount[g] :=
    SetProperty[g, Properties -> Thread[EdgeList[g] -> List /@ Thread[prop -> values]]]
igSetEdgeProperty[g_, prop : specialEdgePropsPattern, values_] /; Length[values] == EdgeCount[g] :=
    transformGraphOptions[ Append[FilterRules[#, Except[prop]], prop -> values]& ][g]
(* Graph[VertexList[g], EdgeList[g], prop -> values, FilterRules[Options[g], Except[prop]]] *)
igSetEdgeProperty[g_, prop_, values_] := $Failed



PackageExport["IGVertexMap"]
IGVertexMap::usage =
    "IGVertexMap[f, prop, graph] maps the function f to the vertex property list of property prop in graph.\n" <>
    "IGVertexMap[f, prop -> pf, graph] maps the function f to the values returned by pf[graph] and assigns the result to the vertex property prop.\n" <>
    "IGVertexMap[f, prop -> {pf1, pf2, \[Ellipsis]}, graph] threads f over {pf1[graph], pf2[graph], \[Ellipsis]} and assigns the result to the vertex property prop.\n" <>
    "IGVertexMap[f, spec] represents an operator form of IGVertexMap that can be applied to a graph.";

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


PackageExport["IGEdgeMap"]
IGEdgeMap::usage =
    "IGEdgeMap[f, prop, graph] maps the function f to the edge property list of property prop in graph.\n" <>
    "IGEdgeMap[f, prop -> pf, graph] maps the function f to the values returned by pf[graph] and assigns the result to the edge property prop.\n" <>
    "IGEdgeMap[f, prop -> {pf1, pf2, \[Ellipsis]}, graph] threads f over {pf1[graph], pf2[graph], \[Ellipsis]} and assigns the result to the edge property prop.\n" <>
    "IGEdgeMap[f, spec] represents an operator form of IGEdgeMap that can be applied to a graph.";

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
      If[nonDistinguishableEdgesQ[g] && Not@MatchQ[prop, specialEdgePropsPattern],
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
      If[nonDistinguishableEdgesQ[g] && Not@MatchQ[prop, specialEdgePropsPattern],
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


(***** Retrieve available edge and vertex property names *****)

(* In 12.1 and later, custom Graph properties are stored as AnnotationRules, not as Properties *)
If[$VersionNumber >= 12.1,
  hasCustomPropQ[g_] := OptionValue[Options[g, AnnotationRules], AnnotationRules] =!= {},
  hasCustomPropQ[g_] := OptionValue[Options[g, Properties], Properties] =!= {}
]

PackageExport["IGVertexPropertyList"]
IGVertexPropertyList::usage = "IGVertexPropertyList[graph] gives the list of available vertex properties in graph.";

standardVertexProperties = {
  VertexCoordinates,
  VertexShape, VertexShapeFunction, VertexSize, VertexStyle,
  VertexLabels, VertexLabelStyle,
  VertexWeight, VertexCapacity
};

SyntaxInformation[IGVertexPropertyList] = {"ArgumentsPattern" -> {_}};
IGVertexPropertyList[g_?IGNullGraphQ] = {};
IGVertexPropertyList[g_ /; GraphQ[g] && hasCustomPropQ[g]] := Sort@DeleteDuplicates[Join @@ PropertyList[{g, VertexList[g]}]]
IGVertexPropertyList[g_ /; GraphQ[g]] := Intersection[PropertyList[g], standardVertexProperties]

PackageExport["IGEdgePropertyList"]
IGEdgePropertyList::usage = "IGEdgePropertyList[graph] gives the list of available edge properties in graph.";

standardEdgeProperties = {
  EdgeStyle, EdgeShapeFunction, EdgeLabels, EdgeLabelStyle,
  EdgeWeight, EdgeCapacity, EdgeCost
};

SyntaxInformation[IGEdgePropertyList] = {"ArgumentsPattern" -> {_}};
IGEdgePropertyList[g_?EmptyGraphQ] = {};
IGEdgePropertyList[g_ /; GraphQ[g] && hasCustomPropQ[g]] := Sort@DeleteDuplicates[Join @@ PropertyList[{g, EdgeList[g]}]]
IGEdgePropertyList[g_ /; GraphQ[g]] := Intersection[PropertyList[g], standardEdgeProperties]
