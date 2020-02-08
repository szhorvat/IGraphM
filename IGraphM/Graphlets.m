(* Mathematica Package *)
(* Created by Mathematica plugin for IntelliJ IDEA *)

(* :Author: szhorvat *)
(* :Date: 2018-10-24 *)
(* :Copyright: (c) 2019 Szabolcs Horvát *)

Package["IGraphM`"]
igContextSetup[igPackagePrivateSymbol]


(*********************)
(***** Graphlets *****)
(*********************)

PackageExport["IGGraphlets"]
IGGraphlets::usage = "IGGraphlets[graph] decomposes a weighted graph into a sum of cliques."; (* TODO *)

Options[IGGraphlets] = { MaxIterations -> 1000 };
SyntaxInformation[IGGraphlets] = {"ArgumentsPattern" -> {_, _.}};
IGGraphlets[graph_?igGraphQ, OptionsPattern[]] :=
    catch@Block[{ig = igMake[graph], basis, mu},
      {basis, mu} = check@ig@"graphlets"[OptionValue[MaxIterations]];
      AssociationThread[
        igVertexNames[graph] /@ igIndexVec[basis],
        mu
      ]
    ]


PackageExport["IGGraphletBasis"]
IGGraphletBasis::usage = "IGGraphletBasis[graph] computes a candidate clique basis."; (* TODO *)

SyntaxInformation[IGGraphletBasis] = {"ArgumentsPattern" -> {_}};
IGGraphletBasis[graph_?igGraphQ] :=
    catch@Block[{ig = igMake[graph], basis, thresholds},
      {basis, thresholds} = check@ig@"graphletBasis"[];
      AssociationThread[
        igVertexNames[graph] /@ igIndexVec[basis],
        thresholds
      ]
    ]


PackageExport["IGGraphletProject"]
IGGraphletProject::usage = "IGGraphletProject[graph, cliques] projects a weighted graph onto the given clique basis."; (* TODO *)

Options[IGGraphletProject] = { MaxIterations -> 1000 };
SyntaxInformation[IGGraphletProject] = {"ArgumentsPattern" -> {_, _, _.}};
IGGraphletProject[graph_?igGraphQ, cliques : {__List}, OptionsPattern[]] :=
    catch@Block[{ig = igMake[graph]},
      AssociationThread[
        cliques,
        check@ig@"graphletProject"[vss[graph] /@ cliques, OptionValue[MaxIterations]]
      ]
    ]
