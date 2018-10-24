(* Mathematica Package *)
(* Created by Mathematica plugin for IntelliJ IDEA *)

(* :Author: szhorvat *)
(* :Date: 2018-10-24 *)
(* :Copyright: (c) 2018 Szabolcs Horvát *)

Package["IGraphM`"]

(*********************************************)
(***** Reading and writing graph formats *****)
(*********************************************)


PackageExport["$IGExportFormats"]
$IGExportFormats::usage = "$IGExportFormats is a list of export formats supported by IGExport.";

$IGExportFormats = {"GraphML"};


IGExport::format   = "`` is not a recognized IGExport format.";
IGExport::infer    = "Cannot infer format of file \"``\".";
IGExport::mixed    = "Exporting mixed graphs to `` is not supported.";
IGExport::propname = "Only string or symbol property names are allowed when exporting to ``.";


PackageExport["IGExportString"]
IGExportString::usage = "IGExportString[graph, format] generates a string corresponding to graph in the given format. See $IGExportFormats for supported formats.";

igExportStringTag::dummy = "igExportStringTag is a symbol used to signal use of ExportString instead of Export in IGExport.";

SyntaxInformation[IGExportString] = {"ArgumentsPattern" -> {_, _}};
IGExportString[g_?GraphQ, format_]  := IGExport[igExportStringTag, g, format]
addCompletion[IGExportString, {0, $IGExportFormats}]


PackageExport["IGExport"]
IGExport::usage =
    "IGExport[file, graph] exports graph to file in a format inferred from the file name.\n" <>
    "IGExport[file, graph, format] exports graph to file in the given format. See $IGExportFormats for supported formats.";

SyntaxInformation[IGExport] = {"ArgumentsPattern" -> {_, _, _.}};
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


(***** GraphML *****)

graphmlTag::dummy = "graphmlTag is an inert symbol used as a Throw/Catch tag in the GraphML exporter.";

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
        graphmlTag
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
