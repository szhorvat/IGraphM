(* Mathematica Package *)
(* Created by Mathematica plugin for IntelliJ IDEA *)

(* :Author: szhorvat *)
(* :Date: 2018-10-24 *)
(* :Copyright: (c) 2018-2020 Szabolcs Horvát *)

Package["IGraphM`"]
igContextSetup[igPackagePrivateSymbol]

PackageImport["IGraphM`LTemplate`"]

(******************************************)
(***** Deterministic graph generators *****)
(******************************************)


PackageExport["IGEmptyGraph"]
IGEmptyGraph::usage =
    "IGEmptyGraph[] gives a graph with no edges or vertices.\n" <>
    "IGEmptyGraph[n] gives a graph with no edges and n vertices.";

SyntaxInformation[IGEmptyGraph] = {"ArgumentsPattern" -> {_., OptionsPattern[]}, "OptionNames" -> optNames[Graph]};
IGEmptyGraph[n : _?Internal`NonNegativeIntegerQ : 0, opt : OptionsPattern[Graph]] := Graph[Range[n], {}, opt]


PackageExport["IGLCF"]
IGLCF::usage =
    "IGLCF[shifts, repeats] creates a graph from LCF notation.\n" <>
    "IGLCF[shifts, repeats, vertexCount] creates a graph from LCF notation with the number of vertices specified.";
Options[IGLCF] = { GraphLayout -> circularEmbedding };
SyntaxInformation[IGLCF] = {"ArgumentsPattern" -> {_, _., _., OptionsPattern[]}, "OptionNames" -> optNames[IGLCF, Graph]};
IGLCF[shifts_?intVecQ, repeats : _?Internal`PositiveMachineIntegerQ : 1, n : (_?Internal`PositiveMachineIntegerQ | Automatic) : Automatic, opt : OptionsPattern[{IGLCF, Graph}]] :=
    catch@Block[{ig = igMakeEmpty[], numVertices},
      numVertices = Replace[n, Automatic :> Length[shifts] repeats];
      If[numVertices < 2,
        IGEmptyGraph[numVertices]
        ,
        check@ig@"fromLCF"[numVertices, shifts, repeats];
        applyGraphOpt[GraphLayout -> OptionValue[GraphLayout], opt]@igToGraph[ig]
      ]
    ]


PackageExport["IGRealizeDegreeSequence"]
IGRealizeDegreeSequence::usage =
    "IGRealizeDegreeSequence[degseq] gives an undirected graph having the given degree sequence.\n" <>
    "IGRealizeDegreeSequence[outdegseq, indegseq] gives a directed graph having the given out- and in-degree sequences.\n";

Options[IGRealizeDegreeSequence] = { Method -> "SmallestFirst" };
igRealizeDegreeSequenceMethods = <|"SmallestFirst" -> 0, "LargestFirst" -> 1, "Index" -> 2|>;
amendUsage[IGRealizeDegreeSequence, "Available Method options: <*Keys[igRealizeDegreeSequenceMethods]*>."];
SyntaxInformation[IGRealizeDegreeSequence] = {"ArgumentsPattern" -> {_, OptionsPattern[]}, "OptionNames" -> optNames[IGRealizeDegreeSequence, Graph]};
IGRealizeDegreeSequence[degrees_?intVecQ, opt : OptionsPattern[{IGRealizeDegreeSequence, Graph}]] :=
    catch@Block[{ig = igMakeEmpty[]},
      check@ig@"realizeDegreeSequence"[degrees, {}, Lookup[igRealizeDegreeSequenceMethods, OptionValue[Method], -1]];
      applyGraphOpt[opt]@igToGraph[ig]
    ]
IGRealizeDegreeSequence[outdeg_?intVecQ, indeg_?intVecQ, opt : OptionsPattern[{IGRealizeDegreeSequence, Graph}]] :=
    catch@Block[{ig = igMakeEmpty[]},
      check@ig@"realizeDegreeSequence"[outdeg, indeg, Lookup[igRealizeDegreeSequenceMethods, OptionValue[Method], -1]];
      applyGraphOpt[opt]@igToGraph[ig]
    ]


PackageExport["IGKaryTree"]
IGKaryTree::usage =
    "IGKaryTree[n] gives a binary tree with n vertices.\n" <>
    "IGKaryTree[n, k] gives a k-ary tree with n vertices.";

Options[IGKaryTree] = {
  DirectedEdges -> False
};
SyntaxInformation[IGKaryTree] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}, "OptionNames" -> optNames[IGKaryTree, Graph]};
IGKaryTree[m_?Internal`NonNegativeMachineIntegerQ, n : _?Internal`PositiveMachineIntegerQ : 2, opt : OptionsPattern[{IGKaryTree, Graph}]] :=
    catch@Block[{ig = igMakeEmpty[]},
      check@ig@"tree"[m, n, OptionValue[DirectedEdges]];
      applyGraphOpt[opt]@igToGraph[ig]
    ]


PackageExport["IGSymmetricTree"]
IGSymmetricTree::usage = "IGSymmetricTree[{k1, k2, \[Ellipsis]}] gives a tree where vertices in the (i+1)st layer have k_i children.";

Options[IGSymmetricTree] = { GraphLayout -> "RadialEmbedding" };
SyntaxInformation[IGSymmetricTree] = {"ArgumentsPattern" -> {_, OptionsPattern[]}, "OptionNames" -> optNames[IGSymmetricTree, Graph]};
IGSymmetricTree[splits_?posIntVecQ, opt : OptionsPattern[{IGSymmetricTree, Graph}]] :=
    With[{edges = igraphGlobal@"symmetricTree"[splits]},
      Graph[
        Range[Length[edges] + 1], edges + 1,
        GraphLayout -> OptionValue[GraphLayout],
        opt
      ]
    ]


PackageExport["IGBetheLattice"]
IGBetheLattice::usage =
    "IGBetheLattice[n] gives the first n layers of a Bethe lattice with coordination number 3.\n" <>
    "IGBetheLattice[n, k] gives the first n layers of a Bethe lattice with coordination number k.";

SyntaxInformation[IGBetheLattice] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}, "OptionNames" -> optNames[IGBetheLattice, Graph]};
IGBetheLattice[n_?Internal`PositiveMachineIntegerQ, k : (i_Integer /; i > 1) : 3, opt : OptionsPattern[]] :=
    IGSymmetricTree[ReplacePart[ConstantArray[k-1, n-1], 1 -> k], opt]


PackageExport["IGFromPrufer"]
IGFromPrufer::usage = "IGFromPrufer[sequence] constructs a tree from a Prüfer sequence.";

SyntaxInformation[IGFromPrufer] = {"ArgumentsPattern" -> {_, OptionsPattern[]}, "OptionNames" -> optNames[IGFromPrufer, Graph]};
IGFromPrufer[vec_?intVecQ, opt : OptionsPattern[Graph]] :=
    catch@Block[{ig = igMakeEmpty[]},
      check@ig@"fromPrufer"[vec - 1];
      applyGraphOpt[opt]@igToGraph[ig]
    ]


PackageExport["IGToPrufer"]
IGToPrufer::usage = "IGToPrufer[tree] gives the Prüfer sequence of a tree.";

SyntaxInformation[IGToPrufer] = {"ArgumentsPattern" -> {_}};
IGToPrufer[graph_?igGraphQ] :=
    sck@Block[{ig = igMakeUnweighted[graph]},
      igIndexVec@ig@"toPrufer"[]
    ]


PackageExport["IGCompleteGraph"]
IGCompleteGraph::usage =
    "IGCompleteGraph[n] gives a complete graph on n vertices.\n" <>
    "IGCompleteGraph[vertices] gives a complete graph on the given vertices.";

IGCompleteGraph::dupl = "The given vertex list must not contain duplicates.";

Options[IGCompleteGraph] = {
  DirectedEdges -> False, SelfLoops -> False,
  GraphLayout -> "CircularEmbedding"
};
SyntaxInformation[IGCompleteGraph] = {"ArgumentsPattern" -> {_, OptionsPattern[]}, "OptionNames" -> optNames[IGCompleteGraph, Graph]};
IGCompleteGraph[m_?Internal`NonNegativeMachineIntegerQ, opt : OptionsPattern[{IGCompleteGraph, Graph}]] :=
    catch@Block[{ig = igMakeEmpty[]},
      check@ig@"completeGraph"[m, OptionValue[DirectedEdges], OptionValue[SelfLoops]];
      applyGraphOpt[GraphLayout -> OptionValue[GraphLayout], opt]@igToGraph[ig]
    ]
IGCompleteGraph[verts_List, opt : OptionsPattern[{IGCompleteGraph, Graph}]] :=
    catch@Block[{ig = igMakeEmpty[]},
      If[Not@DuplicateFreeQ[verts],
        Message[IGCompleteGraph::dupl];
        throw[$Failed]
      ];
      check@ig@"completeGraph"[Length[verts], OptionValue[DirectedEdges], OptionValue[SelfLoops]];
      applyGraphOpt[GraphLayout -> OptionValue[GraphLayout], opt]@igToGraphWithNames[ig, verts]
    ]


PackageExport["IGCompleteAcyclicGraph"]
IGCompleteAcyclicGraph::usage =
    "IGCompleteAcyclicGraph[n] gives a complete acyclic directed graph on n vertices.\n" <>
    "IGCompleteAcyclicGraph[vertices] gives a complete acyclic directed graph on the given vertices.";

IGCompleteAcyclicGraph::dupl = IGCompleteGraph::dupl;

SyntaxInformation[IGCompleteAcyclicGraph] = {"ArgumentsPattern" -> {_, OptionsPattern[]}, "OptionNames" -> optNames[Graph]};
IGCompleteAcyclicGraph[m_?Internal`NonNegativeMachineIntegerQ, opt : OptionsPattern[Graph]] :=
    catch@Block[{ig = igMakeEmpty[]},
      check@ig@"completeCitationGraph"[m, True];
      applyGraphOpt[opt]@igToGraph[ig]
    ]
IGCompleteAcyclicGraph[verts_List, opt : OptionsPattern[Graph]] :=
    catch@Block[{ig = igMakeEmpty[]},
      If[Not@DuplicateFreeQ[verts],
        Message[IGCompleteAcyclicGraph::dupl];
        throw[$Failed]
      ];
      check@ig@"completeCitationGraph"[Length[verts], True];
      applyGraphOpt[opt]@igToGraphWithNames[ig, Reverse[verts]]
    ]


PackageExport["IGKautzGraph"]
IGKautzGraph::usage = "IGKautzGraph[m, n] gives a Kautz graph on m+1 characters with string length n+1.";

SyntaxInformation[IGKautzGraph] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}, "OptionNames" -> optNames[Graph]};
IGKautzGraph[m_?Internal`NonNegativeMachineIntegerQ, n_?Internal`NonNegativeMachineIntegerQ, opt : OptionsPattern[Graph]] :=
    catch@Block[{ig = igMakeEmpty[]},
      check@ig@"kautz"[m, n];
      applyGraphOpt[opt]@igToGraph[ig]
    ]


PackageExport["IGDeBruijnGraph"]
IGDeBruijnGraph::usage = "IGDeBruijnGraph[m, n] gives a De Bruijn graph on m characters and string length n.";

SyntaxInformation[IGDeBruijnGraph] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}, "OptionNames" -> optNames[Graph]};
IGDeBruijnGraph[m_?Internal`NonNegativeMachineIntegerQ, n_?Internal`NonNegativeMachineIntegerQ, opt : OptionsPattern[Graph]] :=
    catch@Block[{ig = igMakeEmpty[]},
      check@ig@"deBruijn"[m, n];
      applyGraphOpt[opt]@igToGraph[ig]
    ]


PackageExport["IGChordalRing"]
IGChordalRing::usage = "IGChordalRing[n, w] gives an extended chordal ring on n vertices, based on the vector or matrix w.";

Options[IGChordalRing] = {
  GraphLayout -> circularEmbedding,
  DirectedEdges -> False,
  SelfLoops -> True,
  MultiEdges -> True, "MultipleEdges" -> "Deprecated"
};
SyntaxInformation[IGChordalRing] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}, "OptionNames" -> optNames[IGChordalRing, Graph]};
IGChordalRing::nv = "A chordal ring must have at least 3 vertices.";
IGChordalRing[n_?Developer`MachineIntegerQ /; n < 3, _, OptionsPattern[]] :=
    (Message[IGChordalRing::nv]; $Failed)
IGChordalRing[n_?Developer`MachineIntegerQ, {}, opt : OptionsPattern[{IGChordalRing, Graph}]] :=
    CycleGraph[n, opt]
IGChordalRing[n_?Developer`MachineIntegerQ, w_?intVecQ, opt : OptionsPattern[{IGChordalRing, Graph}]] :=
    IGChordalRing[n, {w}, opt]
IGChordalRing[n_?Developer`MachineIntegerQ, w_?intMatQ, opt : OptionsPattern[{IGChordalRing, Graph}]] :=
    catch@Block[{ig = igMakeEmpty[]},
      check@ig@"extendedChordalRing"[n, w, OptionValue[DirectedEdges], OptionValue[SelfLoops], multiEdgesOptionReplace@OptionValue[MultiEdges]];
      applyGraphOpt[GraphLayout -> OptionValue[GraphLayout], opt]@igToGraph[ig]
    ]


PackageExport["IGGraphAtlas"]
IGGraphAtlas::usage =
    "IGGraphAtlas[n] gives graph number n from An Atlas of Graphs by Ronald C. Read and Robin J. Wilson, Oxford University Press, 1998. " <>
    "This function is provided for convenience; if you are looking for a specific named graph, use the builtin GraphData function.";

SyntaxInformation[IGGraphAtlas] = {"ArgumentsPattern" -> {_, OptionsPattern[]}, "OptionNames" -> optNames[IGGraphAtlas, Graph]};
IGGraphAtlas[n_?Internal`NonNegativeMachineIntegerQ, opt : OptionsPattern[Graph]] :=
    catch@Block[{ig = igMakeEmpty[]},
      check@ig@"graphAtlas"[n];
      applyGraphOpt[opt]@igToGraph[ig]
    ]


PackageExport["IGTriangularLattice"]
IGTriangularLattice::usage =
    "IGTriangularLattice[n] generates a triangular lattice graph on a size n equilateral triangle using n(n+1)/2 vertices.\n" <>
    "IGTriangularLattice[{m, n}] generates a triangular lattice graph on an m by n rectangle.";

IGTriangularLattice::pimp = "Periodic triangular lattices are not implemented for the equilateral triangle case.";
Options[IGTriangularLattice] = { DirectedEdges -> False, "Periodic" -> False };
SyntaxInformation[IGTriangularLattice] = {"ArgumentsPattern" -> {_, OptionsPattern[]}, "OptionNames" -> optNames[IGTriangularLattice, Graph]};
IGTriangularLattice[{n_?Internal`NonNegativeIntegerQ, m_?Internal`NonNegativeIntegerQ}, opt : OptionsPattern[{IGTriangularLattice, Graph}]] :=
    With[{nm = n (#1 - 1) + #2 &, per = TrueQ@OptionValue["Periodic"]},
      Graph[
        Range[n m],
        Flatten[
          Join[
            Table[nm[i, j] -> nm[Mod[i + 1, m, 1], j], {i, If[per, m, m - 1]}, {j, n}],
            Table[nm[i, j] -> nm[i, Mod[j + 1, n, 1]], {i, m}, {j, If[per, n, n - 1]}],
            Table[nm[i, Mod[j + Boole@OddQ[j], n, 1]] -> nm[Mod[i + 1, m, 1], Mod[j + Boole@EvenQ[j], n, 1]], {i, If[per, m, m - 1]}, {j, If[per, n, n - 1]}]
          ],
          3
        ],
        DirectedEdges -> OptionValue[DirectedEdges],
        FilterRules[{opt}, Options[Graph]],
        If[per,
          {},
          VertexCoordinates -> Join@@Table[{x + 1/4 (-1)^y, Sqrt[3]/2 y}, {x, 1, m}, {y, 1, n}]
        ]
      ]
    ]
IGTriangularLattice[n_?Internal`NonNegativeIntegerQ, opt : OptionsPattern[{IGTriangularLattice, Graph}]] :=
    If[TrueQ@OptionValue["Periodic"],
      Message[IGTriangularLattice::pimp]; $Failed,
      Graph[
        Range[n (n + 1) / 2],
        Flatten@Table[
          With[{a = i (i - 1) / 2 + j},
            {a -> a + i, a -> a + i + 1, a + i -> a + i + 1}
          ],
          {i, n - 1}, {j, i}
        ],
        DirectedEdges -> OptionValue[DirectedEdges],
        FilterRules[{opt}, Options[Graph]],
        VertexCoordinates -> Join @@ Table[{j - i / 2, (n - i) Sqrt[3] / 2}, {i, n}, {j, i}]
      ]
    ]


PackageExport["IGSquareLattice"]
IGSquareLattice::usage = "IGSquareLattice[{d1, d2, \[Ellipsis]}] generates a square grid graph of the given dimensions.";

igSquareLattice[dims_, directed_, mutual_, radius_, periodic_, {opt___}] :=
    catch@Block[{ig = igMakeEmpty[]},
      check@ig@"makeLattice"[dims, radius, directed, mutual, periodic];
      If[Length[dims] === 2 && Not@TrueQ[periodic],
        applyGraphOpt[opt, GraphLayout -> {"GridEmbedding", "Dimension" -> dims}]@igToGraph[ig],
        applyGraphOpt[opt]@igToGraph[ig]
      ]
    ]

Options[IGSquareLattice] = { "Radius" -> 1, DirectedEdges -> False, "Mutual" -> False, "Periodic" -> False };
SyntaxInformation[IGSquareLattice] = {"ArgumentsPattern" -> {_, OptionsPattern[]}, "OptionNames" -> optNames[IGSquareLattice, Graph]};
IGSquareLattice[dims_?nonNegIntVecQ, opt : OptionsPattern[{IGSquareLattice, Graph}]] :=
    igSquareLattice[dims, OptionValue[DirectedEdges], OptionValue["Mutual"], OptionValue["Radius"], OptionValue["Periodic"], {opt}]


(* TODO: remove eventually *)
PackageExport["IGMakeLattice"]
IGMakeLattice::usage = "IGMakeLattice[{d1, d2, \[Ellipsis]}] generates a square grid graph of the given dimensions. It is a synonym for IGSquareLattice";

(* IGMakeLattice is now a synonym for IGSquareLattice, however, it needs to have its separate set of Options *)
Options[IGMakeLattice] = { "Radius" -> 1, DirectedEdges -> False, "Mutual" -> False, "Periodic" -> False };
SyntaxInformation[IGMakeLattice] = {"ArgumentsPattern" -> {_, OptionsPattern[]}, "OptionNames" -> optNames[IGMakeLattice, Graph]};
IGMakeLattice[dims_?nonNegIntVecQ, opt : OptionsPattern[{IGMakeLattice, Graph}]] :=
    igSquareLattice[dims, OptionValue[DirectedEdges], OptionValue["Mutual"], OptionValue["Radius"], OptionValue["Periodic"], {opt}]


(***** Shorthand notation *****)

PackageExport["IGShorthand"]
IGShorthand::usage = "IGShorthand[\"...\"] builds a graph from a shorthand notation such as \"a->b<-c\" or \"a-b,c-d\".";

(* Notes:

    - We cannot use removeSelfLoops[] and removeMultiEdges[] due to the need to support mixed graphs in this function.
    - The symbol sep should not be localized in IGShorthand because it is used in shSplitGroup[].
 *)

sep::dummy = "sep[s] represents a separator while parsing shorthand notations in IGShorthand.";

SetAttributes[shToExpr, Listable];
shToExpr[s_String] :=
    Which[
      StringMatchQ[s, DigitCharacter ..], ToExpression[s],
      True, s
    ]

shSplitGroup[s_sep] := s
shSplitGroup[s_String] := shToExpr@StringTrim@StringSplit[s, ":"]

Options[IGShorthand] = {
  SelfLoops -> False,
  MultiEdges -> False, "MultipleEdges" -> "Deprecated",
  DirectedEdges -> False,
  VertexLabels -> "Name"
};
SyntaxInformation[IGShorthand] = {"ArgumentsPattern" -> {_, OptionsPattern[]}, "OptionNames" -> optNames[IGShorthand, Graph]};
IGShorthand[s_String, opt : OptionsPattern[{IGShorthand, Graph}]] :=
    Module[{comp, seq, edgeSeq, eh, edges, vertices},
      comp = StringSplit[s, ","|"\n"];
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
      canonicalEdgeBlock[
        If[Not@TrueQ@OptionValue[SelfLoops],
          edges = DeleteCases[edges, _[x_, x_]];
        ];
        If[Not@TrueQ@multiEdgesOptionReplace@OptionValue[MultiEdges],
          edges = DeleteDuplicates[edges];
        ]
      ];
      vertices = Catenate@DeleteCases[Level[seq, {2}], _sep];
      Graph[
        vertices,
        edges,
        Sequence @@ FilterRules[{opt}, FilterRules[Options[Graph], Except[VertexLabels|DirectedEdges]]],
        VertexLabels -> OptionValue[VertexLabels]
      ]
    ]


PackageExport["IGExpressionTree"]
IGExpressionTree::usage = "IGExpressionTree[expression] constructs a tree graph from an arbitrary Mathematica expression.";

SetAttributes[exprHead, HoldFirst]
exprHead[expr_Association] := HoldForm[Association] (* Associations are atoms, but Position does descend into them *)
exprHead[expr_ /; AtomQ@Unevaluated[expr]] := HoldForm[expr]
exprHead[expr_] := Extract[Unevaluated[expr], 0, HoldForm]

Options[IGExpressionTree] = {
  VertexLabels -> "Head",
  GraphLayout -> "LayeredEmbedding",
  DirectedEdges -> True
};
SyntaxInformation[IGExpressionTree] = {"ArgumentsPattern" -> {_, OptionsPattern[]}, "OptionNames" -> optNames[IGExpressionTree, Graph]};
(* TODO: Is there an AtomQ that Position will enter into? If yes, exprHead[] will label it incorrectly. *)
IGExpressionTree[expr_, opt : OptionsPattern[{IGExpressionTree, Graph}]] :=
    Module[{vertices, edges, vertexLabels},
      (* We could reverse the vertex list to make the root the first node rather than the last one.
         But this would cause the leaves to be displayed in reverse order with LayeredEmbedding. *)
      vertices = Position[expr, _, {0, Infinity}, Heads -> False];
      edges = MapThread[Rule, {vertices[[;; -2, ;; -2]], vertices[[;; -2]]}];
      vertexLabels =
          Replace[
            OptionValue[VertexLabels],
            {
              "Head" :> Thread[vertices -> Extract[expr, vertices, exprHead]],
              Placed["Head", rest___] :> Thread[vertices -> Extract[expr, vertices, Placed[#, rest]& @* exprHead]],
              "Subexpression" :> Thread[vertices -> Extract[expr, vertices, HoldForm]],
              Placed["Subexpression", rest___] :> Thread[vertices -> Extract[expr, vertices, Placed[#, rest]& @* HoldForm]]
            }
          ];
      Graph[vertices, edges,
        VertexLabels -> vertexLabels,
        opt,
        GraphRoot -> {},
        GraphLayout -> OptionValue[GraphLayout],
        DirectedEdges -> OptionValue[DirectedEdges]
      ]
    ]
