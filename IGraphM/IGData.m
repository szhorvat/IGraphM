(* Mathematica Package *)
(* Created by Mathematica plugin for IntelliJ IDEA *)

(* :Author: szhorvat *)
(* :Date: 2018-10-23 *)
(* :Copyright: (c) 2018 Szabolcs Horvát *)

Package["IGraphM`"]
igContextSetup[igPackagePrivateSymbol]

PackageExport["IGData"]
IGData::usage =
    "IGData[] returns a list of available items.\n" <>
    "IGData[item] returns the requested item.";


(*  IGData  *)

$igData := $igData = zimport@FileNameJoin[{$packageDirectory, "IGData.mz"}];
$igDataCategories := $igDataCategories = GroupBy[Select[Keys[$igData], ListQ], First];
$igDataAll := $igDataAll = Join[$igData, $igDataCategories];

SyntaxInformation[IGData] = {"ArgumentsPattern" -> {_.}};
IGData[] := Keys[$igData]
IGData[item_] := Lookup[$igDataAll, Key[item], Missing["NotAvailable"]]

addCompletion[IGData, {Join[Keys[$igDataCategories], Select[Keys[$igData], StringQ]]}]; (* REPLACE-TAG *)
