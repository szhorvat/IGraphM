(* Mathematica Package *)
(* Created by Mathematica plugin for IntelliJ IDEA *)

(* :Author: szhorvat *)
(* :Date: 2018-10-24 *)
(* :Copyright: (c) 2018-2020 Szabolcs Horvát *)

Package["IGraphM`"]


(***********************************************)
(***** Matrix utilities and graph matrices *****)
(***********************************************)


PackageExport["IGZeroDiagonal"]
IGZeroDiagonal::usage = "IGZeroDiagonal[matrix] replaces the diagonal of matrix with zeros.";

SyntaxInformation[IGZeroDiagonal] = {"ArgumentsPattern" -> {_}};
IGZeroDiagonal[mat_?MatrixQ] :=
    UpperTriangularize[mat, 1] + LowerTriangularize[mat, -1]


PackageExport["IGKirchhoffMatrix"]
IGKirchhoffMatrix::usage =
    "IGKirchhoffMatrix[graph] gives the Kirchhoff matrix, also known as Laplacian matrix of graph.\n" <>
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


PackageExport["IGJointDegreeMatrix"]
IGJointDegreeMatrix::usage =
    "IGJointDegreeMatrix[graph] gives the joint degree matrix of graph. Element i,j of the matrix contains the number of edges connecting degree-i and degree-j vertices.\n" <>
    "IGJointDegreeMatrix[graph, d] gives the d by d joint degree matrix of graph, up to degree d.\n" <>
    "IGJointDegreeMatrix[graph, {dOut, dIn}] gives the dOut by dIn joint degree matrix of graph."

Options[IGJointDegreeMatrix] = { Normalized -> False };
SyntaxInformation[IGJointDegreeMatrix] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};
(* Justification for returning {{}} for empty graph:
 * {{}} is MatrixQ, and it is more general not to fail for this case.
 * However, it is not SquareMatrixQ, and transposing yields {}. *)
IGJointDegreeMatrix[graph_?EmptyGraphQ, opt : OptionsPattern[]] := {{}}
IGJointDegreeMatrix[graph_?igGraphQ, opt : OptionsPattern[]] :=
    With[{sao = SystemOptions["SparseArrayOptions"]},
      Internal`WithLocalSettings[
        SetSystemOptions["SparseArrayOptions" -> "TreatRepeatedEntries" -> Total]
        ,
        Module[{a, b, pairs, res},
          {a, b} = Transpose@IGIndexEdgeList[graph];
          If[UndirectedGraphQ[graph],
            {a, b} = {Join[a, b], Join[b, a]};
          ];
          pairs = Transpose@{VertexOutDegree[graph][[a]], VertexInDegree[graph][[b]]};
          res = SparseArray[pairs -> ConstantArray[1, Length[pairs]]];
          If[TrueQ@OptionValue[Normalized],
            res = If[UndirectedGraphQ[graph], 2, 1] res / Total[res, 2]
          ];
          (* correct the diagonal entries for the undirected case *)
          If[UndirectedGraphQ[graph],
            res = res - DiagonalMatrix[Diagonal[res]/2]
          ];
          res
        ]
        ,
        SetSystemOptions[sao]
      ]
    ]
IGJointDegreeMatrix[graph_?igGraphQ, maxDeg : (_?Internal`PositiveIntegerQ | {_?Internal`PositiveIntegerQ, _?Internal`PositiveIntegerQ}), opt : OptionsPattern[]] :=
    Module[{sa = IGJointDegreeMatrix[graph, opt], max, mdOut, mdIn},
      If[ListQ[maxDeg],
        {mdOut, mdIn} = maxDeg,
        {mdOut, mdIn} = {maxDeg, maxDeg}
      ];
      max = Max[maxDeg, Dimensions[sa]];
      SparseArray[sa, {max, max}][[1;;mdOut, 1;;mdIn]]
    ]


(* Extracting the upper and lower triangular parts *)

upperSparse[sa_SparseArray] :=
    Module[{idx, rows, cols},
      {rows, cols} = Dimensions[sa];
      idx = igraphGlobal@"upperIndexPairPositions"[sa["NonzeroPositions"], cols];
      SparseArray[
        idx[[All, 2]] -> sa["NonzeroValues"][[ idx[[All, 1]] ]],
        List@If[cols > rows, rows * (rows - 1) / 2 + rows * (cols - rows), cols * (cols - 1) / 2],
        sa["Background"]
      ]
    ]

lowerSparse[sa_SparseArray] :=
    Module[{idx, rows, cols},
      {rows, cols} = Dimensions[sa];
      idx = igraphGlobal@"lowerIndexPairPositions"[sa["NonzeroPositions"], cols];
      SparseArray[
        idx[[All, 2]] -> sa["NonzeroValues"][[ idx[[All, 1]] ]],
        List@If[rows > cols, cols * (cols - 1) / 2 + cols * (rows - cols), rows * (rows - 1) / 2],
        sa["Background"]
      ]
    ]

upperDense[mat_?SquareMatrixQ] := Statistics`Library`UpperTriangularMatrixToVector[mat]
upperDense[mat_] := Join @@ Pick[mat, UpperTriangularize[ConstantArray[1, Dimensions[mat]], 1], 1]

lowerDense[mat_] := Join @@ Pick[mat, LowerTriangularize[ConstantArray[1, Dimensions[mat]], -1], 1]


PackageExport["IGTakeUpper"]

IGTakeUpper::arg = IGTakeLower::arg = "A single matrix argument is expected.";

IGTakeUpper::usage = "IGTakeUpper[matrix] extracts the elements of a matrix that are above the diagonal.";
SyntaxInformation[IGTakeUpper] = {"ArgumentsPattern" -> {_}};
IGTakeUpper[matrix_?MatrixQ] :=
    With[{mat = Developer`ToPackedArray[matrix]},
      Switch[mat,
        m_ /; Developer`PackedArrayQ[m, Integer, 2], igraphGlobal@"takeUpperInteger"[mat],
        m_ /; Developer`PackedArrayQ[m, Real, 2], igraphGlobal@"takeUpperReal"[mat],
        m_ /; Developer`PackedArrayQ[m, Complex, 2], igraphGlobal@"takeUpperComplex"[mat],
        _SparseArray, upperSparse[mat],
        _, upperDense[mat]
      ]
    ]
IGTakeUpper[___] := (Message[IGTakeUpper::arg]; $Failed)


PackageExport["IGTakeLower"]

IGTakeLower::usage = "IGTakeLower[matrix] extracts the elements of a matrix that are below the diagonal.";
SyntaxInformation[IGTakeLower] = {"ArgumentsPattern" -> {_}};
IGTakeLower[matrix_?MatrixQ] :=
    With[{mat = Developer`ToPackedArray[matrix]},
      Switch[mat,
        m_ /; Developer`PackedArrayQ[m, Integer, 2], igraphGlobal@"takeLowerInteger"[mat],
        m_ /; Developer`PackedArrayQ[m, Real, 2], igraphGlobal@"takeLowerReal"[mat],
        m_ /; Developer`PackedArrayQ[m, Complex, 2], igraphGlobal@"takeLowerComplex"[mat],
        _SparseArray, lowerSparse[mat],
        _, lowerDense[mat]
      ]
    ]
IGTakeLower[___] := (Message[IGTakeLower::arg]; $Failed)
