(* Mathematica Package *)
(* Created by Mathematica plugin for IntelliJ IDEA *)

(* :Author: szhorvat *)
(* :Date: 2018-10-24 *)
(* :Copyright: (c) 2019 Szabolcs Horvát *)

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
    With[{flann = Make["IGFlann2D"]},
      flann@"setPoints"[pts];
      flann
    ]


PackageScope["flannQueryMultiple"]
flannQueryMultiple::usage = "flannQueryMultiple[]";
flannQueryMultiple[flann_, centres_, dists_] := unpack[flann@"queryMultiple"[centres, dists]]
