(* Mathematica Package *)
(* Created by Mathematica Plugin for IntelliJ IDEA, see http://wlplugin.halirutan.de/ *)

(* :Author: szhorvat *)
(* :Date: 2016-06-12 *)

(* :Copyright: (c) 2018 Szabolcs Horvát *)


(* The following definitions are used in multiple, independently loadable packages. *)

(* General::invopt is not present before Mathematica version 10.3. We set it up manually when needed. *)
If[Not@ValueQ[General::invopt],
  General::invopt = "Invalid value `1` for parameter `2`. Using default value `3`.";
]

IGraphM`IGraphM::mixed = "Mixed graphs are not supported by IGraph/M. Use DirectedGraph to convert undirected edges to two reciprocal directed ones.";

igGraphQ::usage = "igGraphQ[g] checks if g is an igraph-compatible graph.";
igGraphQ = GraphQ[#] && If[MixedGraphQ[#], Message[IGraphM`IGraphM::mixed]; False, True] &;

igDirectedQ::usage = "igDirectedQ[g] checks if g is a directed graph. Empty graphs are considered undirected.";
igDirectedQ[graph_] := DirectedGraphQ[graph] && Not@EmptyGraphQ[graph]

amendUsage::usage = "amendUsage[symbol, stringTempl, templArg1, templArg2, ...] amends the usage message of symbol.";
amendUsage[sym_Symbol, amend_, args___] :=
    Module[{lines},
      lines = StringSplit[sym::usage, "\n"];
      lines[[1]] = lines[[1]] <> " " <> StringTemplate[amend, InsertionFunction -> (ToString[#, InputForm]&)][args];
      sym::usage = StringJoin@Riffle[lines, "\n"]
    ]

optNames::usage = "optNames[sym1, sym2, ...] returns the option names associated with the given symbols.";
optNames[syms___] := Union @@ (Options[#][[All, 1]]& /@ {syms})


applyGraphOpt::usage = "applyGraphOpt[options][graph] applies the given options to graph.";
applyGraphOpt[opt___][graph_] := Graph[graph, Sequence@@FilterRules[{opt}, Options[Graph]]]

applyGraphOpt3D::usage = "applyGraphOpt3D[options][graph] applies the given options to graph using Graph3D.";
applyGraphOpt3D[opt___][graph_] := Graph3D[graph, Sequence@@FilterRules[{opt}, Options[Graph3D]]]


zeroDiagonal::usage = "zeroDiagonal[mat] replaces the diagonal of a matrix with zeros.";
zeroDiagonal[mat_] := UpperTriangularize[mat, 1] + LowerTriangularize[mat, -1]


adjacencyGraph::usage = "adjacencyGraph[vertices, sparseAM, directed]";
adjacencyGraph[vs_, sa_, True] := Graph[vs, {sa, Null}]
adjacencyGraph[vs_, sa_, False] := Graph[vs, {Null, sa}]

removeSelfLoops::usage = "removeSelfLoops[graph] removes any self loops from graph.";
removeSelfLoops[g_?LoopFreeGraphQ] := g (* catches empty case *)
removeSelfLoops[g_] := adjacencyGraph[VertexList[g], zeroDiagonal@AdjacencyMatrix[g], DirectedGraphQ[g]]

removeMultiEdges::usage = "removeMultiEdges[graph] removes any multi-edges from graph.";
removeMultiEdges[g_?MultigraphQ] := adjacencyGraph[VertexList[g], Unitize@AdjacencyMatrix[g], DirectedGraphQ[g]]
removeMultiEdges[g_] := g


$graphLink::usage = "$graphLink is a loopback link used to convert atomic graphs to a compound form.";
transformGraphOptions::usage = "transformGraphOptions[fun][graph] applies fun to the list of options stored in graph.";
transformGraphOptions[fun_][g_?GraphQ] :=
    (
      If[Not@MemberQ[Links[], $graphLink],
        $graphLink = LinkCreate[LinkMode -> Loopback];
      ];
      With[
        {
          expr = AbortProtect[
            LinkWrite[$graphLink, g];
            LinkRead[$graphLink, Hold]
          ]
        },
        Replace[expr, Hold@Graph[v_, e_, opt : _ : {}, rest___] :> Graph[v, e, fun[opt], rest]]
      ]
    )


ruleQ::usage = "ruleQ[expr] gives True if expr is a rule, False otherwise.";
ruleQ[_Rule | _RuleDelayed] = True;
ruleQ[_] = False;


partitionRagged::usage = "partitionRagged[list, lengths] partitions list into parts of the given lengths.";
If[$VersionNumber >= 11.2,
  partitionRagged = TakeList,
  partitionRagged[v_List, l_?VectorQ] := MapThread[Take[v, {#1, #2}] &, Module[{a = Accumulate[l]}, {a - l + 1, a}]]
]


(*
	Numeric codes are for certain special types of completions. Zero means 'don't complete':

	Normal argument     0
	AbsoluteFilename    2
	RelativeFilename    3
	Color               4
	PackageName         7
	DirectoryName       8
	InterpreterType     9
*)

addCompletion::usage = "addCompletion[symbol, argSpec] adds FE auto-completion for symbol.";
addCompletion[fun_Symbol, argSpec_List] :=
    If[$Notebooks,
      With[{compl = SymbolName[fun] -> argSpec},
        FE`Evaluate[FEPrivate`AddSpecialArgCompletion[compl]]
      ]
    ]
