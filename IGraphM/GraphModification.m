
(* Mathematica Package *)
(* Created by Mathematica plugin for IntelliJ IDEA *)

(* :Author: szhorvat *)
(* :Date: 2018-10-24 *)
(* :Copyright: (c) 2018 Szabolcs Horvát *)

Package["IGraphM`"]
igContextSetup[igPackagePrivateSymbol]

(****************************************)
(***** Functions that modify graphs *****)
(****************************************)


(* Warning: this function doesn't preserve the edge ordering or any graph properties *)
(* vertexRename[names_][graph_] :=
    If[VertexCount[graph] == 0,
      graph,
      AdjacencyGraph[names, AdjacencyMatrix[graph], DirectedEdges -> DirectedGraphQ[graph]]
    ] *)

PackageExport["IGConnectNeighborhood"]
IGConnectNeighborhood::usage =
    "IGConnectNeighborhood[graph] connects each vertex in graph to its 2nd order neighbourhood.\n" <>
    "IGConnectNeighborhood[graph, k] connects each vertex in graph to its order k neighbourhood. Weights and other graph properties are discarded.";

SyntaxInformation[IGConnectNeighborhood] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}, "OptionNames" -> optNames[Graph]};
IGConnectNeighborhood[graph_?igGraphQ, k : _?Internal`NonNegativeMachineIntegerQ : 2, opt : OptionsPattern[Graph]] :=
    catch@Block[{ig = igMakeFast[graph]},
      check@ig@"connectNeighborhood"[k];
      applyGraphOpt[opt][igToGraphWithNames[ig, VertexList[graph]]]
    ]


PackageExport["IGMycielskian"]
IGMycielskian::usage = "IGMycielskian[graph] returns the Mycielskian of graph.";

SyntaxInformation[IGMycielskian] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};
IGMycielskian[graph_?igGraphQ, opt : OptionsPattern[Graph]] :=
    catch@Block[{ig = igMakeFast[graph]},
      check@ig@"mycielski"[];
      applyGraphOpt[opt]@igToGraph[ig]
    ]


PackageExport["IGVertexContract"]
IGVertexContract::usage = "IGVertexContract[g, {{v1, v2, \[Ellipsis]}, \[Ellipsis]}] returns a graph in which the specified vertex sets are contracted into single vertices.";

IGVertexContract::inv = "The vertices `` are not present in the graph.";
IGVertexContract::vset = "`` must be a list of disjoint vertex sets.";

Options[IGVertexContract] = { SelfLoops -> False, MultiEdges -> False, "MultipleEdges" -> "Deprecated" };
SyntaxInformation[IGVertexContract] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}, "OptionNames" -> optNames[IGVertexContract, Graph]};
IGVertexContract[graph_?igGraphQ, sets : {___List}, opt : OptionsPattern[{IGVertexContract, Graph}]] :=
    catch@Module[{ig = igMakeFast[graph], allElements = Join @@ sets, fullSets, g, self, multi},
      If[Not@DuplicateFreeQ[allElements],
        Message[IGVertexContract::vset, sets];
        throw[$Failed]
      ];
      If[Not@SubsetQ[VertexList[graph], allElements],
        Message[IGVertexContract::inv, Complement[allElements, VertexList[graph]]];
        throw[$Failed]
      ];
      fullSets = Join[sets, List /@ Complement[VertexList[graph], allElements]];
      check@ig@"contractVertices"@communitiesToMembership[
        VertexList[graph],
        fullSets
      ];
      self = Not@TrueQ@OptionValue[SelfLoops];
      multi = Not@TrueQ@multiEdgesOptionReplace@OptionValue[MultiEdges];
      g = igToGraphWithNames[ig, fullSets[[All,1]] ];
      applyGraphOpt[opt]@Which[
        self && multi, SimpleGraph[g],
        self, removeSelfLoops[g],
        multi, removeMultiEdges[g],
        True, g
      ]
    ]

IGVertexContract[graph_?igGraphQ, arg_, opt : OptionsPattern[]] := Null /; Message[IGVertexContract::vset, arg]


PackageExport["IGSmoothen"]
IGSmoothen::usage = "IGSmoothen[graph] suppresses degree-2 vertices, thus obtaining the smallest topologically equivalent graph. Edge directions are discarded. The weights of merged edges are added up.";

SyntaxInformation[IGSmoothen] = {"ArgumentsPattern" -> {_, OptionsPattern[]}, "OptionNames" -> optNames[Graph]};
IGSmoothen[graph_?igGraphQ, opt : OptionsPattern[Graph]] :=
    catch@Module[{ig = igMakeEmpty[], ig2 = igMake[graph], graph2, deletedIndices},
      deletedIndices = igIndexVec@check@ig@"smoothen"[ManagedLibraryExpressionID[ig2]];
      graph2 = igToWeightedGraphWithNames[ig, VertexList[graph]];
      IGWeightedVertexDelete[graph2, VertexList[graph][[ deletedIndices ]], opt]
    ]

