(* Mathematica Package *)
(* Created by Mathematica plugin for IntelliJ IDEA *)

(* :Author: szhorvat *)
(* :Date: 2018-11-28 *)

(* :Mathematica Version: 10.0 *)
(* :Copyright: (c) 2018 szhorvat *)

(* This is the official installer script for IGraph/M, located at https://raw.githubusercontent.com/szhorvat/IGraphM/master/IGInstaller.m *)

If[{$VersionNumber, $ReleaseNumber} < {10.0, 2},
  (* Note: This is a bit sloppy as the $Version of 10.0.1 is also shown as 10.0, with no release number. *)
  Print["IGraph/M requires Mathematica 10.0.2 or later. You are current running Mathematica ", $Version];
  Abort[]
]

(* We set up a package context to avoid Global` pollution and conflicts with any existing Global` symbols. *)
BeginPackage["IGInstaller`"]
Begin["`Private`"]

check[$Failed] := Throw[$Failed]
check[val_] := val

updateIGraphM[] :=
  Catch@Module[{json, download, target, msg},
    Check[
      json = check@Import["https://api.github.com/repos/szhorvat/IGraphM/releases", "JSON"];
      download = Lookup[First@Lookup[First[json], "assets"], "browser_download_url"];
      msg = "Downloading IGraph/M " <> Lookup[First[json], "tag_name"] <> " ...";
      target = FileNameJoin[{CreateDirectory[], "IGraphM.paclet"}];
      If[$Notebooks,
        PrintTemporary@Labeled[ProgressIndicator[Appearance -> "Necklace"], msg, Right],
        Print[msg]
      ];
      URLSave[download, target]
      ,
      Throw[$Failed]
    ];
    If[FileExistsQ[target], PacletManager`PacletInstall[target], $Failed]
  ]

Print["The currently installed versions of IGraph/M are: ", PacletManager`PacletFind["IGraphM"]]
Module[{res},
  res = Check[updateIGraphM[], $latest, {PacletManager`PacletInstall::samevers}];
  Switch[res,
    $Failed,
    Print["Installation failed. Please install IGraph/M manually. ", Hyperlink["https://github.com/szhorvat/IGraphM#installation"]]
    ,
    $latest,
    Print["The latest version is already installed."],
    _,
    Print["Installing IGraph/M is complete: ", res, "."];
    Print["It can now be loaded using the command Get[\"IGraphM`\"]."]
  ]
]

End[] (* `Private` *)
EndPackage[]