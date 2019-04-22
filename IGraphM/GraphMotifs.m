(* Mathematica Package *)
(* Created by Mathematica plugin for IntelliJ IDEA *)

(* :Author: szhorvat *)
(* :Date: 2018-10-24 *)
(* :Copyright: (c) 2019 Szabolcs Horvát *)

Package["IGraphM`"]
igContextSetup[igPackagePrivateSymbol]

(**********************************************)
(***** Graph motifs, subgraphs, triangles *****)
(**********************************************)


PackageExport["IGDyadCensus"]
IGDyadCensus::usage = "IGDyadCensus[graph] classifies dyad in the graph into mutual, asymmetric or null states.";

SyntaxInformation[IGDyadCensus] = {"ArgumentsPattern" -> {_}};
IGDyadCensus[graph_?igGraphQ] := Block[{ig = igMakeFast[graph]}, AssociationThread[{"Mutual", "Asymmetric", "Null"}, Round@ig@"dyadCensus"[]]]


PackageExport["IGTriadCensus"]
IGTriadCensus::usage = "IGTriadCensus[graph] classifies triads in the graph into 16 possible states, labelled using MAN (mutual, asymmetric, null) notation.";

SyntaxInformation[IGTriadCensus] = {"ArgumentsPattern" -> {_}};
IGTriadCensus[graph_?igGraphQ] :=
    catch@Block[{ig = igMakeFast[graph]},
      AssociationThread[
        {"003", "012", "102", "021D", "021U", "021C", "111D", "111U", "030T", "030C", "201", "120D", "120U", "120C", "210", "300"},
        Round@check@ig@"triadCensus"[]
      ]
    ]


PackageExport["IGMotifs"]
IGMotifs::usage = "IGMotifs[graph, motifSize] returns the motif distribution of graph. See IGIsoclass and IGData for motif ordering.";

Options[IGMotifs] = { DirectedEdges -> Automatic };
SyntaxInformation[IGMotifs] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};
IGMotifs[graph_?igGraphQ, size_?Internal`PositiveIntegerQ, opt : OptionsPattern[]] :=
    catch@Block[{ig = igMakeFast[graph]},
      Switch[OptionValue[DirectedEdges],
        True, ig@"makeDirected"[],
        False, ig@"makeUndirected"[]
      ];
      Round@Developer`FromPackedArray@check@ig@"motifs"[size, ConstantArray[0, size]
      ]
    ]


PackageExport["IGMotifsTotalCount"]
IGMotifsTotalCount::usage = "IGMotifsTotalCount[graph, motifSize] returns the total count of motifs (connected subgraphs) of the given size in the graph.";

Options[IGMotifsTotalCount] = { DirectedEdges -> Automatic };
SyntaxInformation[IGMotifsTotalCount] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};
IGMotifsTotalCount[graph_?igGraphQ, size_?Internal`PositiveIntegerQ, opt : OptionsPattern[]] :=
    Block[{ig = igMakeFast[graph]},
      Switch[OptionValue[DirectedEdges],
        True, ig@"makeDirected"[],
        False, ig@"makeUndirected"[]
      ];
      sck@ig@"motifsNo"[size, ConstantArray[0, size]]
    ]


PackageExport["IGMotifsEstimateTotalCount"]
IGMotifsEstimateTotalCount::usage = "IGMotifsEstimateTotalCount[graph, motifSize, sampleSize] estimates the total count of motifs (connected subgraphs) of the given size in graph, based on a sample of the given size.";

Options[IGMotifsEstimateTotalCount] = { DirectedEdges -> Automatic };
SyntaxInformation[IGMotifsEstimateTotalCount] = {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};
IGMotifsEstimateTotalCount[graph_?igGraphQ, size_?Internal`PositiveIntegerQ, sampleSize_?Internal`PositiveIntegerQ, opt : OptionsPattern[]] :=
    Block[{ig = igMakeFast[graph]},
      Switch[OptionValue[DirectedEdges],
        True, ig@"makeDirected"[],
        False, ig@"makeUndirected"[]
      ];
      sck@ig@"motifsEstimate"[size, ConstantArray[0, size], sampleSize]
    ]


PackageExport["IGTriangles"]
IGTriangles::usage = "IGTriangles[graph] lists all triangles in the graph. Edge directions are ignored.";

SyntaxInformation[IGTriangles] = {"ArgumentsPattern" -> {_}};
IGTriangles[graph_?igGraphQ] :=
    catch@Block[{ig = igMakeFast[graph]},
      Partition[igVertexNames[graph]@igIndexVec@check@ig@"triangles"[], 3]
    ]


PackageExport["IGAdjacentTriangleCount"]
IGAdjacentTriangleCount::usage =
    "IGAdjacentTriangleCount[graph] counts the triangles each vertex participates in. Edge directions are ignored.\n" <>
    "IGAdjacentTriangleCount[graph, vertex] counts the triangles vertex participates in.\n" <>
    "IGAdjacentTriangleCount[graph, {vertex1, vertex2, \[Ellipsis]}] counts the triangles the specified vertices participate in.";

igAdjacentTriangleCount[graph_, vs_] :=
    (* no catch *) Block[{ig = igMakeFast[graph]},
      Round@check@ig@"countAdjacentTriangles"[vss[graph][vs]]
    ]

SyntaxInformation[IGAdjacentTriangleCount] = {"ArgumentsPattern" -> {_, _.}};
IGAdjacentTriangleCount[graph_?igGraphQ, {}] := {}
IGAdjacentTriangleCount[graph_?igGraphQ, vs : (_List | All) : All] := catch@igAdjacentTriangleCount[graph, vs]
IGAdjacentTriangleCount[graph_?igGraphQ, v_] := catch@First@igAdjacentTriangleCount[graph, {v}]


PackageExport["IGTriangleFreeQ"]
IGTriangleFreeQ::usage = "IGTriangleFreeQ[graph] checks if graph is triangle-free.";

SyntaxInformation[IGTriangleFreeQ] = {"ArgumentsPattern" -> {_}};
IGTriangleFreeQ[graph_?igGraphQ] :=
    catch@Block[{ig = igMakeFast[graph]},
      Total[check@ig@"countAdjacentTriangles"[{}]] == 0
    ]
IGTriangleFreeQ[_] := False
