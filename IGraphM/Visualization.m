(* Mathematica Package *)
(* Created by Mathematica plugin for IntelliJ IDEA *)

(* :Author: szhorvat *)
(* :Date: 2018-10-24 *)
(* :Copyright: (c) 2018-2019 Szabolcs Horvát *)

Package["IGraphM`"]


(*****************************************************)
(***** Utility functions for graph visualization *****)
(*****************************************************)


PackageExport["IGAdjacencyMatrixPlot"]
IGAdjacencyMatrixPlot::usage =
    "IGAdjacencyMatrixPlot[graph] plots the adjacency matrix of graph.\n" <>
    "IGAdjacencyMatrixPlot[graph, {v1, v2, \[Ellipsis]}] plots the adjacency matrix of the subgraph induced by the given vertices, using the specified vertex ordering.";

$unconnected::dummy = "$unconnected is used to denote unconnected entries in the implementation of IGAdjacencyMatrixPlot.";

IGAdjacencyMatrixPlot::noprop = "The property `1` does not have a value for some or all edges in the graph.";
IGAdjacencyMatrixPlot::bdname = "The value of the VertexLabels option must be Automatic, \"Name\", \"Index\" or a list of rules.";

Options[IGAdjacencyMatrixPlot] = Options[MatrixPlot] ~Join~ {
  EdgeWeight -> Automatic, "UnconnectedColor" -> Automatic, VertexLabels -> Automatic, "RotateColumnLabels" -> True
};
SetOptions[IGAdjacencyMatrixPlot, Mesh -> Automatic];

SyntaxInformation[IGAdjacencyMatrixPlot] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

IGAdjacencyMatrixPlot[graph_?GraphQ, vs : (_List | All) : All, opt : OptionsPattern[]] :=
    catch@Module[
        {
          $sizeLimit = 50, bigGraphQ, am, vind, rticks, cticks, colRules, rotFun,
          prop = OptionValue[EdgeWeight], mesh = OptionValue[Mesh],
          vertexLabels = OptionValue[VertexLabels], ticks = OptionValue[FrameTicks],
          unconnCol = OptionValue["UnconnectedColor"]
        },

      am = WeightedAdjacencyMatrix[graph, EdgeWeight -> prop];
      If[Not@MatrixQ[am],
        If[MemberQ[IGEdgeProp[prop][graph], _Missing],
          Message[IGAdjacencyMatrixPlot::noprop, prop]
        ];
        throw[$Failed]
      ];
      If[vs === All,
        vind = All,
        Check[vind = VertexIndex[graph, #]& /@ vs, throw[$Failed]];
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
      If[unconnCol =!= Automatic && unconnCol =!= None && Not@ColorQ[unconnCol],
        Message[IGAdjacencyMatrixPlot::invopt, unconnCol, "UnconnectedColor", Automatic];
        unconnCol = Automatic;
      ];
      If[unconnCol === Automatic,
        colRules = OptionValue[ColorRules]
        ,
        am = SparseArray[am["NonzeroPositions"] -> am["NonzeroValues"], Dimensions[am], $unconnected];
        colRules = {$unconnected -> unconnCol};
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
