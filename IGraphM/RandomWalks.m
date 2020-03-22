(* Mathematica Package *)
(* Created by Mathematica plugin for IntelliJ IDEA *)

(* :Author: szhorvat *)
(* :Date: 2018-10-24 *)
(* :Copyright: (c) 2018-2020 Szabolcs Horvát *)

Package["IGraphM`"]
igContextSetup[igPackagePrivateSymbol]

(************************)
(***** Random walks *****)
(************************)


IGraphM::invws = "`1` is not a valid edge weight specification for `2`";

igMakeWithWeightSpec[graph_, edgeWeight_ ] :=
    Switch[edgeWeight,
      Automatic,
      igMake[graph]
      ,
      None,
      igMakeUnweighted[graph]
      ,
      _,
      (* Before setting new weights, remove the old ones with IGUnweighted to work around a but in SetProperty. *)
      With[{newGraph = SetProperty[IGUnweighted[graph], EdgeWeight -> edgeWeight]},
        If[GraphQ[newGraph],
          igMake[newGraph],
          Message[IGraphM::invws, edgeWeight, OutputForm[graph]]; throw[$Failed]
        ]
      ]
    ];


PackageExport["IGRandomWalk"]
IGRandomWalk::usage = "IGRandomWalk[graph, start, steps] takes a random walk of length steps on graph, starting at vertex 'start'. The list of traversed vertices is returned.";

Options[IGRandomWalk] = { EdgeWeight -> Automatic };
SyntaxInformation[IGRandomWalk] = {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};
IGRandomWalk[graph_?igGraphQ, start_, steps_?Internal`NonNegativeMachineIntegerQ, OptionsPattern[]] :=
    catch@Block[{ig = igMakeWithWeightSpec[graph, OptionValue[EdgeWeight]]},
      Part[
        VertexList[graph],
        igIndexVec@check@ig@"randomWalk"[vs[graph][start], steps]
      ]
    ]


PackageExport["IGRandomEdgeIndexWalk"]
IGRandomEdgeIndexWalk::usage = "IGRandomEdgeIndexWalk[graph, start, steps] takes a random walk of length steps on graph, starting at vertex 'start'. The list of indices for traversed edges is returned.";

Options[IGRandomEdgeIndexWalk] = { EdgeWeight -> Automatic };
SyntaxInformation[IGRandomEdgeIndexWalk] = {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};
IGRandomEdgeIndexWalk[graph_?igGraphQ, start_, steps_?Internal`NonNegativeMachineIntegerQ, OptionsPattern[]] :=
    catch@Block[{ig = igMakeWithWeightSpec[graph, OptionValue[EdgeWeight]]},
      igIndexVec@check@ig@"randomEdgeWalk"[vs[graph][start], steps]
    ]


PackageExport["IGRandomEdgeWalk"]
IGRandomEdgeWalk::usage = "IGRandomEdgeWalk[graph, start, steps] takes a random walk of length steps on graph, starting at vertex 'start'. The list of traversed edges is returned.";

Options[IGRandomEdgeWalk] = { EdgeWeight -> Automatic };
SyntaxInformation[IGRandomEdgeWalk] = {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};
IGRandomEdgeWalk[graph_?igGraphQ, start_, steps_?Internal`NonNegativeMachineIntegerQ, OptionsPattern[]] :=
    catch@Block[{ig = igMakeWithWeightSpec[graph, OptionValue[EdgeWeight]]},
      Part[
        EdgeList[graph],
        igIndexVec@check@ig@"randomEdgeWalk"[vs[graph][start], steps]
      ]
    ]
