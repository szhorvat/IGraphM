(* Mathematica Package *)
(* Created by Mathematica plugin for IntelliJ IDEA *)

(* :Author: szhorvat *)
(* :Date: 2018-10-24 *)
(* :Copyright: (c) 2018 Szabolcs Horvát *)

(* LOADING IS CURRENTLY DISABLED
   due to bug igraph core https://github.com/igraph/igraph/issues/869 *)
(*
Package["IGraphM`"]
igContextSetup[igPackagePrivateSymbol]
*)

(*********************)
(***** Graphlets *****)
(*********************)

PackageExport["IGGraphlets"]
IGGraphlets::usage =
    "IGGraphlets[graph]\n" <>
    "IGGraphlets[graph, nIterations]";

SyntaxInformation[IGGraphlets] = {"ArgumentsPattern" -> {_, _.}};
IGGraphlets[graph_?igGraphQ, niter : _?Internal`PositiveMachineIntegerQ : 1000] :=
    catch@Block[{ig = igMake[graph], basis, mu},
      {basis, mu} = check@ig@"graphlets"[niter];
      {igVertexNames[graph] /@ igIndexVec[basis], mu}
    ]


PackageExport["IGGraphletBasis"]
IGGraphletBasis::usage = "IGGraphletBasis[graph]";

SyntaxInformation[IGGraphletBasis] = {"ArgumentsPattern" -> {_}};
IGGraphletBasis[graph_?igGraphQ] :=
    catch@Block[{ig = igMake[graph], basis, thresholds},
      {basis, thresholds} = check@ig@"graphletBasis"[];
      {igVertexNames[graph] /@ igIndexVec[basis], thresholds}
    ]


PackageExport["IGGraphletProject"]
IGGraphletProject::usage =
    "IGGraphletProject[graph, cliques]\n" <>
    "IGGraphletProject[graph, cliques, nIterations]";

SyntaxInformation[IGGraphletProject] = {"ArgumentsPattern" -> {_, _, _.}};
IGGraphletProject[graph_?igGraphQ, cliques : {__List}, niter : _?Internal`PositiveMachineIntegerQ : 1000] :=
    catch@Block[{ig = igMake[graph], clq},
      check@ig@"graphletProject"[Map[vss[graph], cliques, {2}], niter]
    ]
