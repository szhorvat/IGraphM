(* Mathematica Package *)
(* Created by Mathematica Plugin for IntelliJ IDEA, see http://wlplugin.halirutan.de/ *)

(* :Author: szhorvat *)
(* :Date: 2016-06-12 *)

(* :Copyright: (c) 2017 Szabolcs Horvát *)


(* The following definitions are used in multiple, independently loadable packages. *)

IGraphM::mixed = "Mixed graphs are not supported by IGraph/M.";

(* check if the argument is an igraph compatible graph *)
igGraphQ = GraphQ[#] && If[MixedGraphQ[#], Message[IGraphM::mixed]; False, True] &;

(* Check if graph is directed. Empty graphs are considered undrected. *)
igDirectedQ[graph_] := DirectedGraphQ[graph] && Not@EmptyGraphQ[graph]

amendUsage[sym_Symbol, amend_, args___] :=
    Module[{lines},
      lines = StringSplit[sym::usage, "\n"];
      lines[[1]] = lines[[1]] <> " " <> StringTemplate[amend, InsertionFunction -> (ToString[#, InputForm]&)][args];
      sym::usage = StringJoin@Riffle[lines, "\n"]
    ]

optNames[syms___] := Union @@ (Options[#][[All, 1]]& /@ {syms})


applyGraphOpt[opt___][graph_] := Graph[graph, Sequence@@FilterRules[{opt}, Options[Graph]]]
applyGraphOpt3D[opt___][graph_] := Graph3D[graph, Sequence@@FilterRules[{opt}, Options[Graph3D]]]


(* Zero out the diagonal of a square matrix. *)
zeroDiagonal[arg_] := UpperTriangularize[arg, 1] + LowerTriangularize[arg, -1]


removeSelfLoops[g_?LoopFreeGraphQ] := g
removeSelfLoops[g_?igGraphQ] := AdjacencyGraph[VertexList[g], zeroDiagonal@AdjacencyMatrix[g], DirectedEdges -> DirectedGraphQ[g]]

removeMultiEdges[g_ /; igGraphQ[g] && MultigraphQ[g]] := AdjacencyGraph[VertexList[g], Unitize@AdjacencyMatrix[g], DirectedEdges -> DirectedGraphQ[g]]
removeMultiEdges[g_] := g

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

addCompletion[fun_Symbol, argSpec_List] :=
    If[$Notebooks,
      With[{compl = SymbolName[fun] -> argSpec},
        FE`Evaluate[FEPrivate`AddSpecialArgCompletion[compl]]
      ]
    ]
