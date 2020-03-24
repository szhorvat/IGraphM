(* Mathematica Package *)
(* Created by Mathematica plugin for IntelliJ IDEA *)

(* :Author: szhorvat *)
(* :Date: 2018-11-28 *)

(* :Mathematica Version: 10.0 *)
(* :Copyright: (c) 2018-2020 Szabolcs Horvát *)

(* This is the official installer script for IGraph/M, located at https://raw.githubusercontent.com/szhorvat/IGraphM/master/IGInstaller.m *)

If[Not@OrderedQ[{10.0, 2}, {$VersionNumber, $ReleaseNumber}],
  (* Note: This is a bit sloppy as the $Version of 10.0.1 is also shown as 10.0, with no release number. *)
  Print["IGraph/M requires Mathematica 10.0.2 or later. You are current running Mathematica ", $Version];
  Abort[]
]

(* On the Raspberry Pi, we generally support only the latest version of the Wolfram Engine. *)
If[$SystemID === "Linux-ARM",
  If[Not@OrderedQ[{12.0, 0}, {$VersionNumber, $ReleaseNumber}],
    Print["On the Raspberry Pi, IGraph/M requires the Wolfram Engine 12.0 or later. You are currently running the Wolfram Engine ", $Version];
    Abort[]
  ]
]

(* We set up a package context to avoid Global` pollution and conflicts with any existing Global` symbols. *)
BeginPackage["IGInstaller`"]

(* In M12.1, several PacletManager symbols have moved into the System context. We can no longer refer to these symbols
   by their full context without breaking compatibility, therefore we load PacletManager` here. *)
Needs["PacletManager`"]

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
    If[FileExistsQ[target], PacletInstall[target], $Failed]
  ]

Print["The currently installed versions of IGraph/M are: ", #["Version"]& /@ PacletFind["IGraphM"]]
Module[{res, loadCommand},
  res = Check[updateIGraphM[], $latest, {PacletInstall::samevers}];
  Switch[res,
    $Failed,
    Print[
      "Installation failed. Please install IGraph/M manually. ",
      If[$Notebooks, Hyperlink, Identity]["https://github.com/szhorvat/IGraphM#installation"]
    ]
    ,
    $latest,
    Print["The latest version is already installed."],
    _,
    Print["Installing IGraph/M is complete: ", res, "."];
    loadCommand =
        If[$Notebooks,
          Button[
            Defer[Get["IGraphM`"]],
            FrontEndExecute[{
              FrontEnd`SelectionMove[FrontEnd`EvaluationCell[], After, CellGroup],
              FrontEnd`NotebookWrite[FrontEnd`EvaluationNotebook[], #, All],
              FrontEnd`SelectionEvaluateCreateCell[FrontEnd`EvaluationNotebook[]]
            }] &,
            Evaluator -> None,
            Appearance -> None,
            BaseStyle -> "Link"
          ]
          ,
          "<< IGraphM`"
        ];
    Print[
      "It can now be loaded using the command ",
      loadCommand,
      (* On some systems, there may be problems unloading the library. *)
      If[MemberQ[$Packages, "IGraphM`"],
        Unevaluated@Sequence[
          Style["\nWarning: ", Bold],
          "An older version of IGraph/M was already loaded.\nIt is recommended to restart the Mathematica kernel before loading the new version."
        ],
        Unevaluated@Sequence[]
      ]
    ];
  ]
]

End[] (* `Private` *)
EndPackage[]