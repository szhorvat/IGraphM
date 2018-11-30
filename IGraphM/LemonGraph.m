(* Mathematica Package *)
(* Created by Mathematica plugin for IntelliJ IDEA *)

(* :Author: szhorvat *)
(* :Date: 2018-11-30 *)
(* :Copyright: (c) 2018 Szabolcs Horvát *)

Package["IGraphM`"]
igContextSetup[igPackagePrivateSymbol]

PackageImport["IGraphM`LTemplate`"] (* we use Make[] in the definition of lgMake[] *)

(**********************************************************)
(***** Utility functions for communicating with LEMON *****)
(**********************************************************)


(* Make a LemonGraph object *)

(* lgMake is package scope because it is also used in IGLayoutPlanar, defined in GraphLayouts.m *)
PackageScope["lgMake"]
lgMake::usage = "lgMake[graph] converts graph to a LemonGraph object.";
lgMake[g_] :=
    With[{lg = Make["LemonGraph"]},
      lg@"fromEdgeList"[IGIndexEdgeList[g]-1, VertexCount[g]];
      lg
    ]
