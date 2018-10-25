(* Mathematica Package *)
(* Created by Mathematica plugin for IntelliJ IDEA *)

(* :Author: szhorvat *)
(* :Date: 2018-10-24 *)
(* :Copyright: (c) 2018 Szabolcs Horvát *)

Package["IGraphM`"]
igContextSetup[igPackagePrivateSymbol]

(***********************************************)
(***** Matrix utilities and graph matrices *****)
(***********************************************)

(* TODO: Add TakeUpper, TakeLower, TakeNondiagonal *)
(* TODO: Replace any uses of TakeUpper, etc. in other files with new efficient versions *)

PackageExport["IGZeroDiagonal"]
IGZeroDiagonal::usage = "IGZeroDiagonal[matrix] replaces the diagonal of matrix with zeros.";

SyntaxInformation[IGZeroDiagonal] = {"ArgumentsPattern" -> {_}};
IGZeroDiagonal[mat_?MatrixQ] :=
    UpperTriangularize[mat, 1] + LowerTriangularize[mat, -1]


PackageExport["IGKirchhoffMatrix"]
IGKirchhoffMatrix::usage =
    "IGKirchhoffMatrix[graph] returns the Kirchhoff matrix, also known as Laplacian matrix of graph.\n" <>
    "IGKirchhoffMatrix[graph, \"In\"] will place the in-degrees on the diagonal instead of the out-degrees.";

(* The built-in KirchhoffMatrix gives incorrect results for directed graphs. It uses the total degree
   on the diagonal. It should use the out-degree instead.
   KirchhoffMatrix also ignores multi-edges. IGKirchhoffMatrix takes them into account.
 *)
SyntaxInformation[IGKirchhoffMatrix] = {"ArgumentsPattern" -> {_, _.}};
IGKirchhoffMatrix::nkir = "The null graph does not have a Kirchhoff matrix.";
IGKirchhoffMatrix[graph_?IGNullGraphQ, _] := (Message[IGKirchhoffMatrix::nkir]; $Failed)
IGKirchhoffMatrix[graph_?GraphQ] := IGKirchhoffMatrix[graph, "Out"]
IGKirchhoffMatrix[graph_?GraphQ, "Out"] :=
    With[{am = zeroDiagonal@AdjacencyMatrix[graph]},
      DiagonalMatrix@SparseArray@Total[am, {2}] - am
    ]
IGKirchhoffMatrix[graph_?GraphQ, "In"] :=
    With[{am = zeroDiagonal@AdjacencyMatrix[graph]},
      DiagonalMatrix@SparseArray@Total[am] - am
    ]
addCompletion[IGKirchhoffMatrix, {0, {"In", "Out"}}]
