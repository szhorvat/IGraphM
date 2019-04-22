(* Mathematica Package *)
(* Created by Mathematica plugin for IntelliJ IDEA *)

(* :Author: szhorvat *)
(* :Date: 2018-11-30 *)
(* :Copyright: (c) 2019 Szabolcs Horvát *)

Package["IGraphM`"]
igContextSetup[igPackagePrivateSymbol]

(*******************************)
(***** Matching algorithms *****)
(*******************************)


PackageExport["IGMaximumMatching"]
IGMaximumMatching::usage = "IGMaximumMatching[graph] finds a maximum matching of graph. Edge weights are ignored.";
SyntaxInformation[IGMaximumMatching] = {"ArgumentsPattern" -> {_}};
IGMaximumMatching[graph_?igGraphQ] :=
    catch@Block[{lg = lgMake[graph]},
      EdgeList[graph][[ igIndexVec@check@lg@"maximumMatching"[] ]]
    ]


PackageExport["IGMatchingNumber"]
IGMatchingNumber::usage = "IGMatchingNumber[graph] returns the matching number of graph.";
SyntaxInformation[IGMatchingNumber] = {"ArgumentsPattern" -> {_}};
IGMatchingNumber[graph_?igGraphQ] :=
    Block[{lg = lgMake[graph]},
      sck@lg@"matchingNumber"[]
    ]
