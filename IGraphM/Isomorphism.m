(* Mathematica Package *)
(* Created by Mathematica plugin for IntelliJ IDEA *)

(* :Author: szhorvat *)
(* :Date: 2018-10-24 *)
(* :Copyright: (c) 2019 Szabolcs Horvát *)

Package["IGraphM`"]
igContextSetup[igPackagePrivateSymbol]


PackageExport["IGIsomorphicQ"]
IGIsomorphicQ::usage = "IGIsomorphicQ[graph1, graph2] tests if graph1 and graph2 are isomorphic.";

SyntaxInformation[IGIsomorphicQ] = {"ArgumentsPattern" -> {_, _}};
IGIsomorphicQ[g1_?igGraphQ, g2_?igGraphQ] :=
    catch@Block[{ig1 = igMakeUnweighted[g1], ig2 = igMakeUnweighted[g2]},
      If[MultigraphQ[g1] || MultigraphQ[g2],
        check@ig1@"vf2IsomorphicMulti"[ManagedLibraryExpressionID[ig2]],
        check@ig1@"isomorphic"[ManagedLibraryExpressionID[ig2]]
      ]
    ]

(*
igMultigraphVertexEdgeColors[g_] :=
    Module[{e, v},
      {v, e} = Lookup[GroupBy[igraphGlobal@"edgeListSortPairs"[IGIndexEdgeList[g]], Apply[SameQ]], {True, False}, {}];
      {
        Counts[v[[All, 1]]],
        Counts[e]
      }
    ]

igMultigraphIsomorphicQ[g1_, g2_] :=
    (* no catch, catch is in wrapping function *)
      VertexCount[g1] == VertexCount[g2] && EdgeCount[g1] == EdgeCount[g2] &&
        Module[{ig1, ig2, v1, v2, vc1, ec1, vc2, ec2},
          {vc1, ec1} = igMultigraphVertexEdgeColors[g1];
          {vc2, ec2} = igMultigraphVertexEdgeColors[g2];
          v1 = Range@VertexCount[g1];
          v2 = Range@VertexCount[g2];
          ig1 = igMake@Graph[v1, Keys[ec1]];
          ig2 = igMake@Graph[v2, Keys[ec2]];
          check@ig1@"vf2Isomorphic"[
            ManagedLibraryExpressionID[ig2],
            Lookup[vc1, v1, 0], Lookup[vc2, v2, 0],
            Values[ec1], Values[ec2]
      ]
    ]
*)


PackageExport["IGSubisomorphicQ"]
IGSubisomorphicQ::usage = "IGSubisomorphicQ[subgraph, graph] tests if subgraph is contained within graph.";

SyntaxInformation[IGSubisomorphicQ] = {"ArgumentsPattern" -> {_, _}};
IGSubisomorphicQ[subgraph_?EmptyGraphQ, graph_?igGraphQ] := VertexCount[subgraph] <= VertexCount[graph]
IGSubisomorphicQ[subgraph_?igGraphQ, graph_?igGraphQ] :=
    catch@Block[{ig1 = igMakeUnweighted[graph], ig2 = igMakeUnweighted[subgraph]},
      If[Not@SimpleGraphQ[subgraph] || Not@SimpleGraphQ[graph],
        check@ig1@"vf2SubisomorphicMulti"[ManagedLibraryExpressionID[ig2]],
        check@ig1@"subisomorphic"[ManagedLibraryExpressionID[ig2]]
      ]
    ]

(*
igMultigraphSubisomorphicQ[subgraph_?EmptyGraphQ, graph_] :=
    VertexCount[subgraph] <= VertexCount[graph]
igMultigraphSubisomorphicQ[subgraph_, graph_] :=
    (* no catch, catch is in wrapping function *)
    VertexCount[graph] >= VertexCount[subgraph] && EdgeCount[graph] >= EdgeCount[subgraph] &&
        Module[{ig, igs, v, vs, vc, ec, vcs, ecs},
          {vc, ec} = igMultigraphVertexEdgeColors[graph];
          {vcs, ecs} = igMultigraphVertexEdgeColors[subgraph];
          v = Range@VertexCount[graph];
          vs = Range@VertexCount[subgraph];
          ig = igMake@Graph[v, Keys[ec]];
          igs = igMake@Graph[vs, Keys[ecs]];
          check@ig@"vf2Subisomorphic"[
            ManagedLibraryExpressionID[igs],
            Lookup[vc, v, 0], Lookup[vcs, vs, 0],
            Values[ec], Values[ecs]
      ]
    ]
*)


PackageExport["IGGetIsomorphism"]
IGGetIsomorphism::usage = "IGGetIsomorphism[graph1, graph2] returns one isomorphism between graph1 and graph2, if it exists.";

SyntaxInformation[IGGetIsomorphism] = {"ArgumentsPattern" -> {_, _}};
IGGetIsomorphism[graph1_?IGNullGraphQ, igraph2_?IGNullGraphQ] := {<||>}
IGGetIsomorphism[graph1_?igGraphQ, graph2_?igGraphQ] :=
    catch@Block[{ig1 = igMakeFast[graph1], ig2 = igMakeFast[graph2], result},
      result = igIndexVec@check@ig1@"getIsomorphism"[ManagedLibraryExpressionID[ig2]];
      If[result === {}, Return[{}]];
      List@AssociationThread[
        VertexList[graph1],
        igVertexNames[graph2]@result
      ]
    ]


PackageExport["IGGetSubisomorphism"]
IGGetSubisomorphism::usage = "IGGetSubisomorphism[subgraph, graph] returns one subisomorphism from subgraph to graph, if it exists.";

SyntaxInformation[IGGetSubisomorphism] = {"ArgumentsPattern" -> {_, _}};
IGGetSubisomorphism[subgraph_?IGNullGraphQ, graph_?igGraphQ] := {<||>}
IGGetSubisomorphism[subgraph_?igGraphQ, graph_?igGraphQ] :=
    catch@Block[{ig1 = igMakeFast[graph], ig2 = igMakeFast[subgraph], result},
      result = igIndexVec@check@ig1@"getSubisomorphism"[ManagedLibraryExpressionID[ig2]];
      If[result === {}, Return[{}]];
      List@AssociationThread[
        VertexList[subgraph],
        igVertexNames[graph]@result
      ]
    ]


PackageExport["IGIsoclass"]
IGIsoclass::usage = "IGIsoclass[graph] returns the isomorphism class of the graph. Used as the index into the vector returned by motif finding functions. See IGData to get lists of graphs ordered by isoclass.";

SyntaxInformation[IGIsoclass] = {"ArgumentsPattern" -> {_}};
IGIsoclass[graph_?igGraphQ] := Block[{ig = igMakeFast[graph]}, sck@ig@"isoclass"[]]


(* Helper functions for handling vertex and edge colours *)

IGraphM::vcol = "The \"VertexColors\" option must be a list of integers, an association assigning integers to vertices, a vertex property name, or None.";
IGraphM::vcolm = "The vertex property `1` does not contain any values. Assuming all vertices to have the same color.";
IGraphM::vcolp = "The vertex property `1` contains values that cannot be interpreted as colors. Vertex colors must be integers.";
IGraphM::ecol = "The \"EdgeColors\" option must be a list of integers, an association assigning integers to edges, an edge property name, or None.";
IGraphM::ecolm = "The edge property `1` does not contain any values. Assuming all edges to have the same color.";
IGraphM::ecolp = "The edge property `1` contains values that cannot be interpreted as colors. Edge colors must be integers.";
IGraphM::bdecol = "Edge colors: the following edges are not in the graph: ``.";
IGraphM::bdvcol = "Vertex colors: the following vertices are not in the graph: ``.";
IGraphM::vcolcnt = "When vertex colours are specified as a list, the list length must be the same as the vertex count of the graph.";
IGraphM::ecolcnt = "When edge colours are specified as a list, the list length must be the same as the edge count of the graph.";
IGraphM::vcmm = "Only one graph is vertex coloured. Colours will be ignored.";

defaultVF2Colors = {"EdgeColors" -> None, "VertexColors" -> None};

colorCheckVertices[g_, c_] := With[{cm = Complement[Keys[c], VertexList[g]]}, If[cm =!= {}, Message[IGraphM::bdvcol, cm]]];

parseVertexColors[_][None] := {}
parseVertexColors[g_][col_?intVecQ] := (If[VertexCount[g] != Length[col], Message[IGraphM::vcolcnt]; throw[$Failed]]; col)
parseVertexColors[g_][col_?AssociationQ] := (colorCheckVertices[g, col]; Lookup[col, VertexList[g], 0])
parseVertexColors[g_][prop : (_Symbol | _String)] :=
    Module[{values = IGVertexProp[prop][g], result},
      If[MatchQ[values, {__Missing}],
        Message[IGraphM::vcolm, prop]
      ];
      result = Replace[values, _Missing -> 0, {1}];
      If[Not@intVecQ[result],
        Message[IGraphM::vcolp, prop];
        throw[$Failed]
      ];
      result
    ]
parseVertexColors[_][_] := (Message[IGraphM::vcol]; throw[$Failed])

colorCheckEdges[g_, c_] := With[{cm = Complement[Keys[c], EdgeList[g]]}, If[cm =!= {}, Message[IGraphM::bdecol, cm]]];

parseEdgeColors[_][None] := {}
parseEdgeColors[g_][col_?intVecQ] := (If[EdgeCount[g] != Length[col], Message[IGraphM::ecolcnt]; throw[$Failed]]; col)
parseEdgeColors[g_][col_?AssociationQ] :=
    Block[{TwoWayRule = UndirectedEdge},
      Internal`InheritedBlock[{UndirectedEdge},
        SetAttributes[UndirectedEdge, Orderless];
        colorCheckEdges[g, col];
        Lookup[KeyMap[Identity, col] (* allow Orderless to do its job *), EdgeList[g], 0]
      ]
    ]
parseEdgeColors[g_][prop : (_Symbol | _String)] :=
    Module[{values = IGEdgeProp[prop][g], result},
      If[MatchQ[values, {__Missing}],
        Message[IGraphM::ecolm, prop]
      ];
      result = Replace[values, _Missing -> 0, {1}];
      If[Not@intVecQ[result],
        Message[IGraphM::ecolp, prop];
        throw[$Failed]
      ];
      result
    ]
parseEdgeColors[_][_] := (Message[IGraphM::ecol]; throw[$Failed])


(***** Bliss *****)

IGraphM::blissnmg = "Bliss does not support multigraphs. Consider using IGIsomorphicQ.";
blissCheckMulti[graph_] := If[MultigraphQ[graph], Message[IGraphM::blissnmg]; throw[$Failed]]

defaultBlissColors = {"VertexColors" -> None};

blissSplittingHeuristicsNames = {
  "First", "FirstSmallest", "FirstLargest",
  "FirstMaximallyConnected", "FirstSmallestMaximallyConnected", "FirstLargestMaximallyConnected"
};

blissSplittingHeuristics = AssociationThread[blissSplittingHeuristicsNames, Range@Length[blissSplittingHeuristicsNames] - 1];



PackageExport["IGBlissCanonicalLabeling"]
IGBlissCanonicalLabeling::usage =
    "IGBlissCanonicalLabeling[graph] computes a canonical integer labelling of the graph vertices. Using this labelling brings representations of isomorphic graphs to the same form.\n" <>
    "IGBlissCanonicalLabeling[{graph, colorSpec}] computes a canonical integer labelling for the vertices of a vertex coloured graph.";

amendUsage[IGBlissCanonicalLabeling,
  " Available values for the \"SplittingHeuristics\" option: ``. The labelling depends on the splitting heuristics used.",
  blissSplittingHeuristicsNames
];


Options[IGBlissCanonicalLabeling] = { "SplittingHeuristics" -> "First" };
SyntaxInformation[IGBlissCanonicalLabeling] = {"ArgumentsPattern" -> {{__}, OptionsPattern[]}};
IGBlissCanonicalLabeling[graph_?igGraphQ, opt : OptionsPattern[]] :=
    catch@Block[{ig = igMakeFast[graph]},
      blissCheckMulti[graph];
      AssociationThread[
        VertexList[graph],
        igIndexVec@check@ig@"blissCanonicalPermutation"[Lookup[blissSplittingHeuristics, OptionValue["SplittingHeuristics"], -1], {}]
      ]
    ]
IGBlissCanonicalLabeling[{graph_?igGraphQ, col : OptionsPattern[]}, opt : OptionsPattern[]] :=
    catch@Block[{ig = igMakeFast[graph], vcol},
      blissCheckMulti[graph];
      vcol = parseVertexColors[graph]@OptionValue[defaultBlissColors, {col}, "VertexColors"];
      AssociationThread[
        VertexList[graph],
        igIndexVec@check@ig@"blissCanonicalPermutation"[Lookup[blissSplittingHeuristics, OptionValue["SplittingHeuristics"], -1], vcol]
      ]
    ]



PackageExport["IGBlissCanonicalPermutation"]
IGBlissCanonicalPermutation::usage =
    "IGBlissCanonicalPermutation[graph] returns a permutation that, when applied to the adjacency matrices of isomorphic graphs, brings them to the same form.\n" <>
    "IGBlissCanonicalPermutation[{graph, colorSpec}] returns the canonical vertex permutation of a vertex coloured graph.";

Options[IGBlissCanonicalPermutation] = { "SplittingHeuristics" -> "First" };
SyntaxInformation[IGBlissCanonicalPermutation] = {"ArgumentsPattern" -> {{__}, OptionsPattern[]}};
IGBlissCanonicalPermutation[graph_?igGraphQ, opt : OptionsPattern[]] :=
    catch@Block[{ig = igMakeFast[graph]},
      blissCheckMulti[graph];
      InversePermutation@igIndexVec@check@ig@"blissCanonicalPermutation"[Lookup[blissSplittingHeuristics, OptionValue["SplittingHeuristics"], -1], {}]
    ]
IGBlissCanonicalPermutation[{graph_?igGraphQ, col : OptionsPattern[]}, opt : OptionsPattern[]] :=
    catch@Block[{ig = igMakeFast[graph], vcol},
      blissCheckMulti[graph];
      vcol = parseVertexColors[graph]@OptionValue[defaultBlissColors, {col}, "VertexColors"];
      InversePermutation@igIndexVec@check@ig@"blissCanonicalPermutation"[Lookup[blissSplittingHeuristics, OptionValue["SplittingHeuristics"], -1], vcol]
    ]


PackageExport["IGBlissCanonicalGraph"]
IGBlissCanonicalGraph::usage =
    "IGBlissCanonicalGraph[graph] returns a canonical graph of graph, based on the canonical integer labelling.\n" <>
    "IGBlissCanonicalGraph[{graph, colorSpec}] returns a canonical graph of a vertex coloured graph, based on the canonical integer labelling. Vertex colours will be stored in the \"Color\" vertex property.";

Options[IGBlissCanonicalGraph] = { "SplittingHeuristics" -> "First" };
SyntaxInformation[IGBlissCanonicalGraph] = {"ArgumentsPattern" -> {{__}, OptionsPattern[]}};
IGBlissCanonicalGraph[graph_?IGNullGraphQ, opt : OptionsPattern[]] := Graph[{},{}] (* the empty graph has no adjacency matrix *)
IGBlissCanonicalGraph[graph_?igGraphQ, opt : OptionsPattern[]] :=
    catch@With[{perm = check@IGBlissCanonicalPermutation[graph, opt], am = AdjacencyMatrix[graph]},
      AdjacencyGraph[ am[[perm, perm]], DirectedEdges -> DirectedGraphQ[graph] ]
    ]
IGBlissCanonicalGraph[spec : {graph_?igGraphQ, col : OptionsPattern[]}, opt : OptionsPattern[]] :=
      catch@Module[{perm, am, vcol},
        vcol = parseVertexColors[graph]@OptionValue[defaultBlissColors, {col}, "VertexColors"];
        If[vcol === {}, (* vertex colours set to None or the graph has no vertices *)
          IGBlissCanonicalGraph[graph]
          ,
          perm = check@IGBlissCanonicalPermutation[spec, opt];
          am = AdjacencyMatrix[graph];
          AdjacencyGraph[
            am[[perm, perm]],
            DirectedEdges -> DirectedGraphQ[graph],
            Properties -> Thread[ Range@VertexCount[graph] -> List /@ Thread[ "Color" -> vcol[[perm]] ] ]
          ]
        ]
      ]


PackageExport["IGBlissIsomorphicQ"]
IGBlissIsomorphicQ::usage =
    "IGBlissIsomorphicQ[graph1, graph2] tests if graph1 and graph2 are isomorphic using the Bliss algorithm.\n" <>
    "IGBlissIsomorphicQ[{graph1, colorSpec}, {graph2, colorSpec}] tests if two vertex coloured graphs are isomorphic using the Bliss algorithm.";

Options[IGBlissIsomorphicQ] = { "SplittingHeuristics" -> "First" };
SyntaxInformation[IGBlissIsomorphicQ] = {"ArgumentsPattern" -> {{__}, {__}, OptionsPattern[]}};
IGBlissIsomorphicQ[graph1_?igGraphQ, graph2_?igGraphQ, opt : OptionsPattern[]] :=
    catch@Block[{ig1 = igMakeFast[graph1], ig2 = igMakeFast[graph2]},
      blissCheckMulti /@ {graph1, graph2};
      check@ig1@"blissIsomorphic"[ManagedLibraryExpressionID[ig2], Lookup[blissSplittingHeuristics, OptionValue["SplittingHeuristics"], -1], {}, {}]
    ]
IGBlissIsomorphicQ[{graph1_?igGraphQ, col1 : OptionsPattern[]}, {graph2_?igGraphQ, col2 : OptionsPattern[]}, opt : OptionsPattern[]] :=
    catch@Block[{ig1 = igMakeFast[graph1], ig2 = igMakeFast[graph2], vcol1, vcol2},
      blissCheckMulti /@ {graph1, graph2};
      vcol1 = parseVertexColors[graph1]@OptionValue[defaultBlissColors, {col1}, "VertexColors"];
      vcol2 = parseVertexColors[graph2]@OptionValue[defaultBlissColors, {col2}, "VertexColors"];
      check@ig1@"blissIsomorphic"[ManagedLibraryExpressionID[ig2], Lookup[blissSplittingHeuristics, OptionValue["SplittingHeuristics"], -1], vcol1, vcol2]
    ]


PackageExport["IGBlissGetIsomorphism"]
IGBlissGetIsomorphism::usage =
    "IGBlissGetIsomorphism[graph1, graph2] returns one isomorphism between graph1 and graph2, if it exists.\n" <>
    "IGBlissGetIsomorphism[{graph1, colorSpec}, {graph2, colorSpec}] returns one isomorphism between two vertex colored graphs, if it exists.";

Options[IGBlissGetIsomorphism] = { "SplittingHeuristics" -> "First" };
SyntaxInformation[IGBlissGetIsomorphism] = {"ArgumentsPattern" -> {{__}, {__}, OptionsPattern[]}};
IGBlissGetIsomorphism[graph1_?IGNullGraphQ, graph2_?IGNullGraphQ] := {<||>}
IGBlissGetIsomorphism[graph1_?igGraphQ, graph2_?igGraphQ, opt : OptionsPattern[]] :=
    catch@Block[{ig1 = igMakeFast[graph1], ig2 = igMakeFast[graph2], result},
      blissCheckMulti /@ {graph1, graph2};
      result = igIndexVec@check@ig1@"blissFindIsomorphism"[ManagedLibraryExpressionID[ig2], Lookup[blissSplittingHeuristics, OptionValue["SplittingHeuristics"], -1], {}, {}];
      If[result === {}, Return[{}]];
      List@AssociationThread[
        VertexList[graph1],
        igVertexNames[graph2]@result
      ]
    ]
IGBlissGetIsomorphism[{graph1_?igGraphQ, col1 : OptionsPattern[]}, {graph2_?igGraphQ, col2 : OptionsPattern[]}, opt : OptionsPattern[]] :=
    catch@Block[{ig1 = igMakeFast[graph1], ig2 = igMakeFast[graph2], result, vcol1, vcol2},
      blissCheckMulti /@ {graph1, graph2};
      vcol1 = parseVertexColors[graph1]@OptionValue[defaultBlissColors, {col1}, "VertexColors"];
      vcol2 = parseVertexColors[graph2]@OptionValue[defaultBlissColors, {col2}, "VertexColors"];
      If[IGNullGraphQ[graph1] && IGNullGraphQ[graph2],
        Return[{<||>}]
      ];
      result = igIndexVec@check@ig1@"blissFindIsomorphism"[ManagedLibraryExpressionID[ig2], Lookup[blissSplittingHeuristics, OptionValue["SplittingHeuristics"], -1], vcol1, vcol2];
      If[result === {}, Return[{}]];
      List@AssociationThread[
        VertexList[graph1],
        igVertexNames[graph2]@result
      ]
    ]



PackageExport["IGBlissAutomorphismCount"]
IGBlissAutomorphismCount::usage =
    "IGBlissAutomorphismCount[graph] returns the number of automorphisms of graph.\n" <>
    "IGBlissAutomorphismCount[{graph, colorSpec}] returns the number of automorphisms of a vertex coloured graph.";

Options[IGBlissAutomorphismCount] = { "SplittingHeuristics" -> "First" };
SyntaxInformation[IGBlissAutomorphismCount] = {"ArgumentsPattern" -> {{__}, OptionsPattern[]}};
IGBlissAutomorphismCount[graph_?igGraphQ, opt : OptionsPattern[]] :=
    catch@Block[{ig = igMakeFast[graph]},
      blissCheckMulti[graph];
      ToExpression@check@ig@"blissAutomorphismCount"[Lookup[blissSplittingHeuristics, OptionValue["SplittingHeuristics"], -1], {}]
    ]

IGBlissAutomorphismCount[{graph_?igGraphQ, col : OptionsPattern[]}, opt : OptionsPattern[]] :=
    catch@Block[{ig = igMakeFast[graph], vcol},
      blissCheckMulti[graph];
      vcol = parseVertexColors[graph]@OptionValue[defaultBlissColors, {col}, "VertexColors"];
      ToExpression@check@ig@"blissAutomorphismCount"[Lookup[blissSplittingHeuristics, OptionValue["SplittingHeuristics"], -1], vcol]
    ]



PackageExport["IGBlissAutomorphismGroup"]
IGBlissAutomorphismGroup::usage =
    "IGBlissAutomorphismGroup[graph] returns a set of generators for the automorphism group of graph. It is not guaranteed to be minimal.\n" <>
    "IGBlissAutomorphismGroup[{graph, colorSpec}] returns a set of generators for the automorphism group of a vertex coloured graph.";

Options[IGBlissAutomorphismGroup] = { "SplittingHeuristics" -> "First" };
SyntaxInformation[IGBlissAutomorphismGroup] = {"ArgumentsPattern" -> {{__}, OptionsPattern[]}};
IGBlissAutomorphismGroup[graph_?GraphQ, opt : OptionsPattern[]] :=
    catch@Block[{ig = igMakeFast[graph]},
      blissCheckMulti[graph];
      igIndexVec@check@ig@"blissAutomorphismGroup"[Lookup[blissSplittingHeuristics, OptionValue["SplittingHeuristics"], -1], {}]
    ]

IGBlissAutomorphismGroup[{graph_?GraphQ, col : OptionsPattern[]}, opt : OptionsPattern[]] :=
    catch@Block[{ig = igMakeFast[graph], vcol},
      blissCheckMulti[graph];
      vcol = parseVertexColors[graph]@OptionValue[defaultBlissColors, {col}, "VertexColors"];
      igIndexVec@check@ig@"blissAutomorphismGroup"[Lookup[blissSplittingHeuristics, OptionValue["SplittingHeuristics"], -1], vcol]
    ]


(***** VF2 *****)

(* This function has been updated to filter both multi-graphs and graphs with self-loops.
   The name stays the same, but keep in mind that self-loops are disallowed too. *)
IGraphM::vf2nmg = "VF2 does not support non-simple graphs. Consider using IGIsomorphicQ.";
vf2CheckMulti[graph_] := If[Not@SimpleGraphQ[graph], Message[IGraphM::vf2nmg]; throw[$Failed]]


PackageExport["IGVF2IsomorphicQ"]
IGVF2IsomorphicQ::usage =
    "IGVF2IsomorphicQ[graph1, graph2] tests if graph1 and graph2 are isomorphic using the VF2 algorithm.\n" <>
    "IGVF2IsomorphicQ[{graph1, colorSpec}, {graph2, colorSpec}] tests if vertex or edge coloured graphs graph1 and graph2 are isomorphic.";

SyntaxInformation[IGVF2IsomorphicQ] = {"ArgumentsPattern" -> {{__}, {__}}};

IGVF2IsomorphicQ[{graph1_?igGraphQ, opt1 : OptionsPattern[]}, {graph2_?igGraphQ, opt2 : OptionsPattern[]}] :=
    catch@Block[{ig1 = igMake[graph1], ig2 = igMake[graph2], vcol1, vcol2, ecol1, ecol2},
      vf2CheckMulti /@ {graph1, graph2};
      vcol1 = parseVertexColors[graph1]@OptionValue[defaultVF2Colors, {opt1}, "VertexColors"];
      vcol2 = parseVertexColors[graph2]@OptionValue[defaultVF2Colors, {opt2}, "VertexColors"];
      ecol1 = parseEdgeColors[graph1]@OptionValue[defaultVF2Colors, {opt1}, "EdgeColors"];
      ecol2 = parseEdgeColors[graph2]@OptionValue[defaultVF2Colors, {opt2}, "EdgeColors"];
      check@ig1@"vf2Isomorphic"[ManagedLibraryExpressionID[ig2], vcol1, vcol2, ecol1, ecol2]
    ]

IGVF2IsomorphicQ[graph1_?igGraphQ, graph2_?igGraphQ] :=
    catch@Block[{ig1 = igMake[graph1], ig2 = igMake[graph2]},
      vf2CheckMulti /@ {graph1, graph2};
      check@ig1@"vf2Isomorphic"[ManagedLibraryExpressionID[ig2], {}, {}, {}, {}]
    ]


PackageExport["IGVF2FindIsomorphisms"]
IGVF2FindIsomorphisms::usage =
    "IGVF2FindIsomorphisms[graph1, graph2] finds all isomorphisms between graph1 and graph2 using the VF2 algorithm.\n" <>
    "IGVF2FindIsomorphisms[graph1, graph2, n] finds at most n isomorphisms between graph1 and graph2.\n" <>
    "IGVF2FindIsomorphisms[{graph1, colorSpec}, {graph2, colorSpec}] finds all isomorphisms between vertex or edge coloured graphs graph1 and graph2.\n" <>
    "IGVF2FindIsomorphisms[{graph1, colorSpec}, {graph2, colorSpec}, n] finds at most n isomorphisms between vertex or edge coloured graphs graph1 and graph2.";

SyntaxInformation[IGVF2FindIsomorphisms] = {"ArgumentsPattern" -> {{__}, {__}, _.}};

IGVF2FindIsomorphisms[{graph1_?igGraphQ, opt1 : OptionsPattern[]}, {graph2_?igGraphQ, opt2 : OptionsPattern[]}, max : (_?Internal`PositiveMachineIntegerQ | All | Infinity) : All] :=
    catch@Block[{ig1 = igMake[graph1], ig2 = igMake[graph2], vcol1, vcol2, ecol1, ecol2, n, result},
      vf2CheckMulti /@ {graph1, graph2};
      n = Replace[max, All|Infinity -> -1];
      vcol1 = parseVertexColors[graph1]@OptionValue[defaultVF2Colors, {opt1}, "VertexColors"];
      vcol2 = parseVertexColors[graph2]@OptionValue[defaultVF2Colors, {opt2}, "VertexColors"];
      ecol1 = parseEdgeColors[graph1]@OptionValue[defaultVF2Colors, {opt1}, "EdgeColors"];
      ecol2 = parseEdgeColors[graph2]@OptionValue[defaultVF2Colors, {opt2}, "EdgeColors"];
      result = igIndexVec@check@ig1@"vf2FindIsomorphisms"[ManagedLibraryExpressionID[ig2], n, vcol1, vcol2, ecol1, ecol2];
      AssociationThread[
        VertexList[graph1],
        igVertexNames[graph2][#]
      ]& /@ result
    ]

IGVF2FindIsomorphisms[graph1_?igGraphQ, graph2_?igGraphQ, max : (_?Internal`PositiveMachineIntegerQ | All | Infinity) : All] :=
    catch@Block[{ig1 = igMake[graph1], ig2 = igMake[graph2], n, result},
      vf2CheckMulti /@ {graph1, graph2};
      n = Replace[max, All|Infinity -> -1];
      result = igIndexVec@check@ig1@"vf2FindIsomorphisms"[ManagedLibraryExpressionID[ig2], n, {}, {}, {}, {}];
      AssociationThread[
        VertexList[graph1],
        igVertexNames[graph2][#]
      ]& /@ result
    ]


PackageExport["IGVF2GetIsomorphism"]
IGVF2GetIsomorphism::usage = "IGVF2GetIsomorphism[graph1, graph2] returns one isomorphism between graph1 and graph2, if it exists.\n" <>
    "IGVF2GetIsomorphism[{graph1, colorSpec}, {graph2, colorSpec}] returns one isomorphism between two vertex or edge colored graphs, if it exists.";

SyntaxInformation[IGVF2GetIsomorphism] = {"ArgumentsPattern" -> {{__}, {__}}};

IGVF2GetIsomorphism[graph1_?igGraphQ, graph2_?igGraphQ] :=
    IGVF2FindIsomorphisms[graph1, graph2, 1]
IGVF2GetIsomorphism[cg1: {graph1_?igGraphQ, opt1 : OptionsPattern[]}, cg2: {graph2_?igGraphQ, opt2 : OptionsPattern[]}] :=
    IGVF2FindIsomorphisms[cg1, cg2, 1]


PackageExport["IGVF2SubisomorphicQ"]
IGVF2SubisomorphicQ::usage =
    "IGVF2SubisomorphicQ[subgraph, graph] tests if subgraph is contained in graph using the VF2 algorithm.\n" <>
    "IGVF2SubisomorphicQ[{subgraph, colorSpec}, {graph, colorSpec}] tests if vertex or edge coloured subgraph is contained in graph.";

SyntaxInformation[IGVF2SubisomorphicQ] = {"ArgumentsPattern" -> {{__}, {__}}};

IGVF2SubisomorphicQ[{subgraph_?igGraphQ, optsub : OptionsPattern[]}, {graph_?igGraphQ, opt : OptionsPattern[]}] :=
    catch@Block[{ig1 = igMake[graph], ig2 = igMake[subgraph], vcol, vcolsub, ecol, ecolsub},
      vf2CheckMulti /@ {subgraph, graph};
      vcol    = parseVertexColors[graph]@OptionValue[defaultVF2Colors, {opt}, "VertexColors"];
      vcolsub = parseVertexColors[subgraph]@OptionValue[defaultVF2Colors, {optsub}, "VertexColors"];
      ecol    = parseEdgeColors[graph]@OptionValue[defaultVF2Colors, {opt}, "EdgeColors"];
      ecolsub = parseEdgeColors[subgraph]@OptionValue[defaultVF2Colors, {optsub}, "EdgeColors"];
      check@ig1@"vf2Subisomorphic"[ManagedLibraryExpressionID[ig2], vcol, vcolsub, ecol, ecolsub]
    ]

IGVF2SubisomorphicQ[subgraph_?igGraphQ, graph_?igGraphQ] :=
    catch@Block[{ig1 = igMake[graph], ig2 = igMake[subgraph]},
      vf2CheckMulti /@ {subgraph, graph};
      check@ig1@"vf2Subisomorphic"[ManagedLibraryExpressionID[ig2], {}, {}, {}, {}]
    ]


PackageExport["IGVF2FindSubisomorphisms"]
IGVF2FindSubisomorphisms::usage =
    "IGVF2FindSubisomorphisms[subgraph, graph] finds all subisomorphisms from subgraph to graph using the VF2 algorithm.\n" <>
    "IGVF2FindSubisomorphisms[subgraph, graph, n] finds at most n subisomorphisms from subgraph to graph.\n" <>
    "IGVF2FindSubisomorphisms[{subgraph, colorSpec}, {graph, colorSpec}] finds all subisomorphisms from vertex or edge coloured subgraph to graph.\n" <>
    "IGVF2FindSubisomorphisms[{subgraph, colorSpec}, {graph, colorSpec}, n] finds at most n subisomorphisms from vertex or edge coloured subgraph to graph.";

SyntaxInformation[IGVF2FindSubisomorphisms] = {"ArgumentsPattern" -> {{__}, {__}, _.}};

IGVF2FindSubisomorphisms[{subgraph_?igGraphQ, optsub : OptionsPattern[]}, {graph_?igGraphQ, opt : OptionsPattern[]}, max : (_?Internal`PositiveMachineIntegerQ | All | Infinity) : All] :=
    catch@Block[{ig1 = igMake[graph], ig2 = igMake[subgraph], vcol, vcolsub, ecol, ecolsub, n, result},
      vf2CheckMulti /@ {subgraph, graph};
      n = Replace[max, All|Infinity -> -1];
      vcol    = parseVertexColors[graph]@OptionValue[defaultVF2Colors, {opt}, "VertexColors"];
      vcolsub = parseVertexColors[subgraph]@OptionValue[defaultVF2Colors, {optsub}, "VertexColors"];
      ecol    = parseEdgeColors[graph]@OptionValue[defaultVF2Colors, {opt}, "EdgeColors"];
      ecolsub = parseEdgeColors[subgraph]@OptionValue[defaultVF2Colors, {optsub}, "EdgeColors"];
      result = igIndexVec@check@ig1@"vf2FindSubisomorphisms"[ManagedLibraryExpressionID[ig2], n, vcol, vcolsub, ecol, ecolsub];
      AssociationThread[
        VertexList[subgraph],
        igVertexNames[graph][#]
      ]& /@ result
    ]

IGVF2FindSubisomorphisms[subgraph_?igGraphQ, graph_?igGraphQ, max : (_?Internal`PositiveMachineIntegerQ | All | Infinity) : All] :=
    catch@Block[{ig1 = igMake[graph], ig2 = igMake[subgraph], n, result},
      vf2CheckMulti /@ {subgraph, graph};
      n = Replace[max, All|Infinity -> -1];
      result = igIndexVec@check@ig1@"vf2FindSubisomorphisms"[ManagedLibraryExpressionID[ig2], n, {}, {}, {}, {}];
      AssociationThread[
        VertexList[subgraph],
        igVertexNames[graph][#]
      ]& /@ result
    ]


PackageExport["IGVF2GetSubisomorphism"]
IGVF2GetSubisomorphism::usage = "IGVF2GetSubisomorphism[subgraph, graph] returns one subisomorphism from subgraph to graph, if it exists.\n" <>
    "IGVF2GetSubisomorphism[{subgraph, colorSpec}, {graph, colorSpec}] returns one subisomorphism from a vertex or edge coloured subgraph to graph, if it exists.";

SyntaxInformation[IGVF2GetSubisomorphism] = {"ArgumentsPattern" -> {{__}, {__}}};

IGVF2GetSubisomorphism[subgraph_?igGraphQ, graph_?igGraphQ] :=
    IGVF2FindSubisomorphisms[subgraph, graph, 1]
IGVF2GetSubisomorphism[cg1: {subgraph_?igGraphQ, opt1 : OptionsPattern[]}, cg2: {graph_?igGraphQ, opt2 : OptionsPattern[]}] :=
    IGVF2FindSubisomorphisms[cg1, cg2, 1]


PackageExport["IGVF2IsomorphismCount"]
IGVF2IsomorphismCount::usage =
    "IGVF2IsomorphismCount[graph1, graph2] returns the number of isomorphisms between graph1 and graph2.\n" <>
    "IGVF2IsomorphismCount[{graph1, colorSpec}, {graph2, colorSpec}] returns the number of isomorphisms between vertex or edge coloured graphs graph1 and graph2. Note that this is not the same as simply counting the automorphisms of one graph if their colourings differ.";

SyntaxInformation[IGVF2IsomorphismCount] = {"ArgumentsPattern" -> {{__}, {__}}};

IGVF2IsomorphismCount[{graph1_?igGraphQ, opt1 : OptionsPattern[]}, {graph2_?igGraphQ, opt2 : OptionsPattern[]}] :=
    catch@Block[{ig1 = igMake[graph1], ig2 = igMake[graph2], vcol1, vcol2, ecol1, ecol2},
      vf2CheckMulti /@ {graph1, graph2};
      vcol1 = parseVertexColors[graph1]@OptionValue[defaultVF2Colors, {opt1}, "VertexColors"];
      vcol2 = parseVertexColors[graph2]@OptionValue[defaultVF2Colors, {opt2}, "VertexColors"];
      ecol1 = parseEdgeColors[graph1]@OptionValue[defaultVF2Colors, {opt1}, "EdgeColors"];
      ecol2 = parseEdgeColors[graph2]@OptionValue[defaultVF2Colors, {opt2}, "EdgeColors"];
      check@ig1@"vf2IsomorphismCount"[ManagedLibraryExpressionID[ig2], vcol1, vcol2, ecol1, ecol2]
    ]

IGVF2IsomorphismCount[graph1_?igGraphQ, graph2_?igGraphQ] :=
    catch@Block[{ig1 = igMake[graph1], ig2 = igMake[graph2]},
      vf2CheckMulti /@ {graph1, graph2};
      check@ig1@"vf2IsomorphismCount"[ManagedLibraryExpressionID[ig2], {}, {} ,{}, {}]
    ]


PackageExport["IGVF2SubisomorphismCount"]
IGVF2SubisomorphismCount::usage =
    "IGVF2SubisomorphismCount[subgraph, graph] returns the number of mappings from subgraph to graph.\n" <>
    "IGVF2SubisomorphismCount[{subgraph, colorSpec}, {graph, colorSpec}] returns the number of mappings from vertex or edge coloured subgraph to graph.";

SyntaxInformation[IGVF2SubisomorphismCount] = {"ArgumentsPattern" -> {{__}, {__}}};

IGVF2SubisomorphismCount[{subgraph_?igGraphQ, optsub : OptionsPattern[]}, {graph_?igGraphQ, opt : OptionsPattern[]}] :=
    catch@Block[{ig1 = igMake[graph], ig2 = igMake[subgraph], vcol, vcolsub, ecol, ecolsub},
      vf2CheckMulti /@ {subgraph, graph};
      vcol    = parseVertexColors[graph]@OptionValue[defaultVF2Colors, {opt}, "VertexColors"];
      vcolsub = parseVertexColors[subgraph]@OptionValue[defaultVF2Colors, {optsub}, "VertexColors"];
      ecol    = parseEdgeColors[graph]@OptionValue[defaultVF2Colors, {opt}, "EdgeColors"];
      ecolsub = parseEdgeColors[subgraph]@OptionValue[defaultVF2Colors, {optsub}, "EdgeColors"];
      check@ig1@"vf2SubisomorphismCount"[ManagedLibraryExpressionID[ig2], vcol, vcolsub, ecol, ecolsub]
    ]

IGVF2SubisomorphismCount[subgraph_?igGraphQ, graph_?igGraphQ] :=
    catch@Block[{ig1 = igMake[graph], ig2 = igMake[subgraph]},
      vf2CheckMulti /@ {subgraph, graph};
      check@ig1@"vf2SubisomorphismCount"[ManagedLibraryExpressionID[ig2], {}, {}, {}, {}]
    ]


(***** LAD *****)

defaultLADColors = {"VertexColors" -> None};


PackageExport["IGLADSubisomorphicQ"]
IGLADSubisomorphicQ::usage =
    "IGLADSubisomorphicQ[subgraph, graph] tests if subgraph is contained in graph. Use the \"Induced\" -> True option to look for induced subgraphs.\n" <>
    "IGLADSubisomorphicQ[{subgraph, colorSpec}, {graph, colorSpec}] tests if a vertex coloured subgraph is contained in graph.";

Options[IGLADSubisomorphicQ] = { "Induced" -> False };
SyntaxInformation[IGLADSubisomorphicQ] = {"ArgumentsPattern" -> {{__}, {__}, OptionsPattern[]}};

IGLADSubisomorphicQ[subgraph_?igGraphQ, graph_?igGraphQ, opt : OptionsPattern[]] :=
    Block[{ig1 = igMakeFast[graph], ig2 = igMakeFast[subgraph]},
      sck@ig1@"ladSubisomorphic"[ManagedLibraryExpressionID[ig2], OptionValue["Induced"]]
    ]

IGLADSubisomorphicQ[{subgraph_?igGraphQ, colsub : OptionsPattern[]}, {graph_?igGraphQ, col : OptionsPattern[]}, opt : OptionsPattern[]] :=
    catch@Block[{vcol, vcolsub},
      vcol    = parseVertexColors[graph]@OptionValue[defaultLADColors, {col}, "VertexColors"];
      vcolsub = parseVertexColors[subgraph]@OptionValue[defaultLADColors, {colsub}, "VertexColors"];
      If[vcol === {} || vcolsub === {},
        If[vcol =!= vcolsub, Message[IGraphM::vcmm]];
        IGLADSubisomorphicQ[subgraph, graph, opt]
        ,
        Block[{ig1 = igMakeFast[graph], ig2 = igMakeFast[subgraph]},
          check@ig1@"ladSubisomorphicColored"[
            ManagedLibraryExpressionID[ig2], Boole@TrueQ@OptionValue["Induced"],
            Flatten@Position[vcol, #, {1}] - 1& /@ vcolsub
          ]
        ]
      ]
    ]


PackageExport["IGLADGetSubisomorphism"]
IGLADGetSubisomorphism::usage =
    "IGLADGetSubisomorphism[subgraph, graph] returns one subisomorphism from subgraph to graph, if it exists.\n" <>
    "IGLADGetSubisomorphism[{subgraph, colorSpec}, {graph, colorSpec}] returns one subisomorphism from a vertex coloured subgraph to graph.";

Options[IGLADGetSubisomorphism] = { "Induced" -> False };
SyntaxInformation[IGLADGetSubisomorphism] = {"ArgumentsPattern" -> {{__}, {__}, OptionsPattern[]}};

IGLADGetSubisomorphism[subgraph_?igGraphQ, graph_?igGraphQ, opt : OptionsPattern[]] :=
    catch@Block[{ig1 = igMakeFast[graph], ig2 = igMakeFast[subgraph], result},
      result = igIndexVec@check@ig1@"ladGetSubisomorphism"[ManagedLibraryExpressionID[ig2], OptionValue["Induced"]];
      If[result === {}, Return@If[IGNullGraphQ[subgraph], {<||>}, {}]];
      List@AssociationThread[
        VertexList[subgraph],
        igVertexNames[graph][result]
      ]
    ]

IGLADGetSubisomorphism[{subgraph_?igGraphQ, colsub : OptionsPattern[]}, {graph_?igGraphQ, col : OptionsPattern[]}, opt : OptionsPattern[]] :=
    catch@Block[{vcol, vcolsub},
      vcol    = parseVertexColors[graph]@OptionValue[defaultLADColors, {col}, "VertexColors"];
      vcolsub = parseVertexColors[subgraph]@OptionValue[defaultLADColors, {colsub}, "VertexColors"];
      If[vcol === {} || vcolsub === {},
        If[vcol =!= vcolsub, Message[IGraphM::vcmm]];
        IGLADGetSubisomorphism[subgraph, graph, opt]
        ,
        Block[{ig1 = igMakeFast[graph], ig2 = igMakeFast[subgraph], result},
          result = igIndexVec@check@ig1@"ladGetSubisomorphismColored"[
            ManagedLibraryExpressionID[ig2], Boole@TrueQ@OptionValue["Induced"],
            Flatten@Position[vcol, #, {1}] - 1& /@ vcolsub
          ];
          If[result === {}, Return@If[IGNullGraphQ[subgraph], {<||>}, {}]];
          List@AssociationThread[
            VertexList[subgraph],
            igVertexNames[graph][result]
          ]
        ]
      ]
    ]


PackageExport["IGLADFindSubisomorphisms"]
IGLADFindSubisomorphisms::usage =
    "IGLADFindSubisomorphisms[subgraph, graph] finds all subisomorphisms from subgraph to graph.\n" <>
    "IGLADFindSubisomorphisms[{subgraph, colorSpec}, {graph, colorSpec}] finds all subisomorphisms from a vertex coloured subgraph to graph.";

Options[IGLADFindSubisomorphisms] = { "Induced" -> False };
SyntaxInformation[IGLADFindSubisomorphisms] = {"ArgumentsPattern" -> {{__}, {__}, OptionsPattern[]}};

IGLADFindSubisomorphisms[subgraph_?IGNullGraphQ, graph_?igGraphQ, opt : OptionsPattern[]] := {<||>} 
IGLADFindSubisomorphisms[subgraph_?igGraphQ, graph_?igGraphQ, opt : OptionsPattern[]] :=
    catch@Block[{ig1 = igMakeFast[graph], ig2 = igMakeFast[subgraph], result},
      result = igIndexVec@check@ig1@"ladFindSubisomorphisms"[ManagedLibraryExpressionID[ig2], Boole@TrueQ@OptionValue["Induced"], {}];
      AssociationThread[
        VertexList[subgraph],
        igVertexNames[graph][#]
      ]& /@ result
    ]

IGLADFindSubisomorphisms[{subgraph_?igGraphQ, colsub : OptionsPattern[]}, {graph_?igGraphQ, col : OptionsPattern[]}, opt : OptionsPattern[]] :=
    catch@Block[{ig1 = igMakeFast[graph], ig2 = igMakeFast[subgraph], result, vcol, vcolsub, domain},
      vcol    = parseVertexColors[graph]@OptionValue[defaultLADColors, {col}, "VertexColors"];
      vcolsub = parseVertexColors[subgraph]@OptionValue[defaultLADColors, {colsub}, "VertexColors"];
      If[IGNullGraphQ[subgraph], Return[{<||>}]]; (* special case: one match for the null pattern for consistency *)
      If[vcol === {} || vcolsub === {},
        If[vcol =!= vcolsub, Message[IGraphM::vcmm]];
        domain = {}
        ,
        domain = Flatten@Position[vcol, #, {1}] - 1& /@ vcolsub;
      ];
      result = igIndexVec@check@ig1@"ladFindSubisomorphisms"[ManagedLibraryExpressionID[ig2], Boole@TrueQ@OptionValue["Induced"], domain];
      AssociationThread[
        VertexList[subgraph],
        igVertexNames[graph][#]
      ]& /@ result
    ]


PackageExport["IGLADSubisomorphismCount"]
IGLADSubisomorphismCount::usage =
    "IGLADSubisomorphismCount[subgraph, graph] counts subisomorphisms from subgraph to graph.\n" <>
    "IGLADSubisomorphismCount[{subgraph, colorSpec}, {graph, colorSpec}] counts subisomorphisms from a vertex coloured subgraph to graph.";

Options[IGLADSubisomorphismCount] = { "Induced" -> False };
SyntaxInformation[IGLADSubisomorphismCount] = {"ArgumentsPattern" -> {{__}, {__}, OptionsPattern[]}};

IGLADSubisomorphismCount[subgraph_?igGraphQ, graph_?igGraphQ, opt : OptionsPattern[]] :=
    Block[{ig1 = igMakeFast[graph], ig2 = igMakeFast[subgraph], result},
      sck@ig1@"ladCountSubisomorphisms"[ManagedLibraryExpressionID[ig2], OptionValue["Induced"]]
    ]

IGLADSubisomorphismCount[{subgraph_?igGraphQ, colsub : OptionsPattern[]}, {graph_?igGraphQ, col : OptionsPattern[]}, opt : OptionsPattern[]] :=
    catch@Block[{ig1 = igMakeFast[graph], ig2 = igMakeFast[subgraph], result, vcol, vcolsub, domain},
      vcol    = parseVertexColors[graph]@OptionValue[defaultLADColors, {col}, "VertexColors"];
      vcolsub = parseVertexColors[subgraph]@OptionValue[defaultLADColors, {colsub}, "VertexColors"];
      If[vcol === {} || vcolsub === {},
        If[vcol =!= vcolsub, Message[IGraphM::vcmm]];
        domain = {}
        ,
        domain = Flatten@Position[vcol, #, {1}] - 1& /@ vcolsub;
      ];
      check@ig1@"ladCountSubisomorphismsColored"[ManagedLibraryExpressionID[ig2], Boole@TrueQ@OptionValue["Induced"], domain]
    ]


(***** Vertex and edge transitivity *****)

PackageExport["IGVertexTransitiveQ"]
IGVertexTransitiveQ::usage = "IGVertexTransitiveQ[graph] tests if graph is vertex transitive.";

IGVertexTransitiveQ::nmg = "Multigraphs are not supported.";
SyntaxInformation[IGVertexTransitiveQ] = {"ArgumentsPattern" -> {_}};
IGVertexTransitiveQ[graph_?EmptyGraphQ] = True;
IGVertexTransitiveQ[graph_?igGraphQ] :=
    If[MultigraphQ[graph],
      Message[IGVertexTransitiveQ::nmg];
      $Failed
      ,
      With[{elems = Range@VertexCount[graph]},
        GroupOrbits[PermutationGroup@IGBlissAutomorphismGroup[graph], elems] === {elems}
      ]
    ]
IGVertexTransitiveQ[_] = False;


PackageExport["IGEdgeTransitiveQ"]
IGEdgeTransitiveQ::usage = "IGEdgeTransitiveQ[graph] tests if graph is edge transitive.";

IGEdgeTransitiveQ::nmg = IGVertexTransitiveQ::nmg;
SyntaxInformation[IGEdgeTransitiveQ] = {"ArgumentsPattern" -> {_}};
IGEdgeTransitiveQ[graph_?igGraphQ] :=
    If[MultigraphQ[graph],
      Message[IGEdgeTransitiveQ::nmg];
      $Failed
      ,
      IGVertexTransitiveQ@LineGraph[graph]
    ]
IGEdgeTransitiveQ[_] = False;


PackageExport["IGSymmetricQ"]
IGSymmetricQ::usage = "IGSymmetricQ[graph] tests if graph is symmetric, i.e. it is both vertex transitive and edge transitive.";

IGSymmetricQ::nmg = IGVertexTransitiveQ::nmg;
SyntaxInformation[IGSymmetricQ] = {"ArgumentsPattern" -> {_}};
IGSymmetricQ[graph_?igGraphQ] :=
    If[MultigraphQ[graph],
      Message[IGSymmetricQ::nmg];
      $Failed
      ,
      IGVertexTransitiveQ[graph] && IGEdgeTransitiveQ[graph]
    ]
IGSymmetricQ[_] = False;


(***** Homeomorphism *****)

PackageExport["IGHomeomorphicQ"]
IGHomeomorphicQ::usage = "IGHomeomorphicQ[graph1, graph2] tests if graph1 and graph2 are homeomorphic. Edge directions are ignored.";

SyntaxInformation[IGHomeomorphicQ] = {"ArgumentsPattern" -> {_, _}};
IGHomeomorphicQ[g1_?igGraphQ, g2_?igGraphQ] :=
    catch@IGIsomorphicQ[check@IGSmoothen[g1], check@IGSmoothen[g2]]


(***** Other functions related to isomorphism *****)

PackageExport["IGSelfComplementaryQ"]
IGSelfComplementaryQ::usage = "IGSelfComplementaryQ[graph] tests if graph is self-complementary.";

IGSelfComplementaryQ::nmg = "`1` is not a simple graph.";
SyntaxInformation[IGSelfComplementaryQ] = {"ArgumentsPattern" -> {_}};
IGSelfComplementaryQ[graph_?igGraphQ] :=
    catch@Block[{ig = igMakeFast[graph]},
      If[Not@SimpleGraphQ[graph],
        Message[IGSelfComplementaryQ::nmg, OutputForm[graph]];
        throw[$Failed]
      ];
      If[UndirectedGraphQ[graph],
        If[
          Sort@VertexDegree[graph] == Sort[VertexCount[graph] - 1 - VertexDegree[graph]]
          ,
          check@ig@"selfComplementaryQ"[],
          False
        ]
        , (* directed case *)
        If[
          Sort@VertexInDegree[graph] == Sort[VertexCount[graph] - 1 - VertexInDegree[graph]] &&
          Sort@VertexOutDegree[graph] == Sort[VertexCount[graph] - 1 - VertexOutDegree[graph]]
          ,
          check@ig@"selfComplementaryQ"[],
          False
        ]
      ]
    ]
IGSelfComplementaryQ[_] := False


PackageExport["IGColoredSimpleGraph"]
IGColoredSimpleGraph::usage = "IGColoredSimpleGraph[graph] encodes a non-simple graph as an edge- and vertex-colored simple graph, returned as {simpleGraph, \"VertexColors\" -> vcol, \"EdgeColors\" -> ecol} where vertex colors represent self-loop multiplicities and edge colors represent edge multiplicities. The output is suitable for use by isomorphism functions.";
SyntaxInformation[IGColoredSimpleGraph] = {"ArgumentsPattern" -> {_}};
IGColoredSimpleGraph[graph_?igGraphQ] :=
    catch@Block[{ig = igMakeUnweighted[graph], new = igMakeEmpty[], vc, ec},
      {vc, ec} = check@new@"coloredSimpleGraph"[ManagedLibraryExpressionID[ig]];
      (*igSetVertexProperty[igSetEdgeProperty[igToGraph[new], "EdgeColor", ec], "VertexColor", vc]*)
      {igToGraph[new], "VertexColors" -> vc, "EdgeColors" -> ec}
    ]
