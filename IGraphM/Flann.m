(* Mathematica Package *)
(* Created by Mathematica plugin for IntelliJ IDEA *)

(* :Author: szhorvat *)
(* :Date: 2018-10-24 *)
(* :Copyright: (c) 2018-2020 Szabolcs Horvát *)

Package["IGraphM`"]
igContextSetup[igPackagePrivateSymbol]

PackageImport["IGraphM`LTemplate`"] (* we use Make in makeFlann *)

(********************************************************)
(***** Helper functions for the nanoflann interface *****)
(********************************************************)

unpack[packed_] :=
    With[{len = First[packed]},
      partitionRagged[packed[[len + 2 ;; All]], packed[[2 ;; len + 1]]]
    ]


PackageScope["makeFlann"]
makeFlann::usage = "makeFlann[]";
makeFlann[pts_] :=
    With[
      {
        flann = Switch[
          Dimensions[pts],
          {_, 2}, Make["IGFlann2D"],
          {_, 3}, Make["IGFlann3D"],
          _, Message[IGraphM::flanndim]; throw[$Failed]
        ]
      },
      flann@"setPoints"[pts];
      flann
    ]


PackageScope["flannQueryMultiple"]
flannQueryMultiple::usage = "flannQueryMultiple[]";
flannQueryMultiple[flann_, centres_, dists_] := unpack[flann@"queryMultiple"[centres, dists]]


IGraphM::flanndim = "Geometric nearest neighbor computations are only supported in 2 or 3 dimensions.";