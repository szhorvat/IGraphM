(* ::Package:: *)

(* Mathematica Source File  *)
(* Created by Mathematica Plugin for IntelliJ IDEA *)
(* :Author: szhorvat *)
(* :Date: 2016-12-02 *)

(* This script requires the following versions of Mathematica:
 * 10.0 for writing 10.0-compatible notebooks
 * 10.4 for evaluating the documentation notebook and obtaining a result that renders correctly in 10.x
 * 11.1 for creating old-style documentation indexes
 * 11.2 for creating new-style documentation indexes
 *)

$appName = "IGraphM";
$LTemplateRepo = "~/Repos/LTemplate";

Print["\nScript running in Mathematica ", $Version];

(* Prints an error message and aborts the script *)
printAbort[str_] := (Print["ABORTING: ", Style[str, Red, Bold]]; Quit[])
If[$VersionNumber < 10.0, printAbort["Mathematica 10.0 or later required."]]


(* Determine script directory, whether the script is run in the terminal,
   is read with Get[], or is evaluated within the front end. *)
$dir =
    Which[
      $InputFileName =!= "", DirectoryName[$InputFileName],
      $Notebooks, NotebookDirectory[],
      True, printAbort["Cannot determine script directory."]
    ]


(* Verify that we have git *)
rg = RunProcess[{"git", "--version"}];
If[rg === $Failed || rg["ExitCode"] != 0,
  printAbort["git is not available."]
]


SetDirectory[$dir]


$buildDir = FileNameJoin[{$dir, "release"}]


If[DirectoryQ[$buildDir],
  Print["Deleting release directory."];
  DeleteDirectory[$buildDir, DeleteContents->True]
]


$appSource = FileNameJoin[{$dir, $appName}]
$appTarget = FileNameJoin[{$buildDir, $appName}]


CreateDirectory[$appTarget, CreateIntermediateDirectories -> True]


(* Copy tracked files into the release directory *)
Print["git-archive IGraphM"]

$gitArch = FileNameJoin[{$buildDir, $appName<>".tar"}]

RunProcess[{"git", "archive", "--format", "tar", "-o", $gitArch, "HEAD:IGraphM"}]

ExtractArchive[$gitArch, $appTarget]
DeleteFile[$gitArch]


SetDirectory[$appTarget]


Print["Loading IGraphM paclet..."]
PacletDirectoryAdd[$dir]
Needs[$appName <> "`"]


source = Import["IGraphM.m", "String"];


igDataCompletionRepl = StringTemplate["addCompletion[IGData, ``];"]@
    ToString[{Join[Keys[IGraphM`Private`$igDataCategories],
      Select[Keys[IGraphM`Private`$igData], StringQ]]},
      InputForm
    ]

repl = {
  "Get[\"LTemplate`LTemplatePrivate`\"]" -> "Get[\"IGraphM`LTemplate`LTemplatePrivate`\"]",
  "\"LazyLoading\" -> False" -> "\"LazyLoading\" -> True",
  "addCompletion[IGData, {Join[Keys[$igDataCategories], Select[Keys[$igData], StringQ]]}];" -> igDataCompletionRepl
};


source = StringReplace[source, repl, Length[repl]];


template = StringTemplate[source, Delimiters -> {"%%", "%%", "<%", "%>"}];


versionData = <|
  "version" -> Lookup[PacletInformation[$appName], "Version"],
  "mathversion" -> Lookup[PacletInformation[$appName], "WolframVersion"],
  "date" -> DateString[{"MonthName", " ", "DayShort", ", ", "Year"}]
|>


source = template[versionData];


Export[FileNameJoin[{$appTarget, "IGraphM.m"}], source, "String"]


(* Replace unicode characters with their Mathematica FullForm, to ensure that
   source file work identically regardless of the value of $CharacterEncoding *)
Print["\nRe-encoding source files as ASCII"]
AddPath["PackageTools"]
Needs["PackageTools`"]
With[{$appTarget = $appTarget},
  MRun[
    MCode[
      UsingFrontEnd[
        NotebookSave@NotebookOpen[FileNameJoin[{$appTarget, #}]]& /@
            {"IGraphM.m", "PacletInfo.m", "Common.m", "Utilities.m"}
      ]
    ]
  ]
]

DeleteFile["BuildSettings.m"]


(***** Copy LTemplate *****)

Print["\ngit-archive LTemplate"]

SetDirectory@CreateDirectory["LTemplate"]

RunProcess[{"git", "archive", "--format", "tar", "-o", "LTemplate.tar", "HEAD:LTemplate", "--remote=" <> $LTemplateRepo}]

ExtractArchive["LTemplate.tar"]

DeleteDirectory["Documentation", DeleteContents -> True]
DeleteFile[{"IncludeFiles/Doxyfile", "LTemplate.tar"}]

ResetDirectory[] (* LTemplate *)


(***** Clean up documentation *****)

Print["\nProcessing documentation."]
SetDirectory@FileNameJoin[{"Documentation", "English", "Tutorials"}]

Print["Evaluating..."]
With[{$buildDir = $buildDir},
  MRun[
    MCode[
      PacletManager`PacletDirectoryAdd[$buildDir];
      SetOptions[First[$Output], FormatType -> StandardForm];
      Internal`$MessageMenu = False; (* get 10.0-compatible message formatting from 11.0 *)
      RewriteNotebook[
        NBFEProcess[NotebookEvaluate[#, InsertResults -> True]&]
      ]["IGDocumentation.nb"]
    ],
    (* Note 1: 10.0 can't evaluate from a MathLink-controlled kernel due to a bug, use newer version *)
    (* Note 2: Use at least 10.4, which supports PlotTheme. *)
    (* Note 3: Use at least 11.0, which supports the Region property of PolyhedronData. *)
    "11.0"
  ]
]

Print["Rewriting..."]

taggingRules = {
  "ModificationHighlight" -> False,
  "ColorType" -> "TutorialColor",
  "Metadata" -> {
    "built" -> ToString@DateList[],
    "history" -> {versionData["version"], "", "", ""},
    "context" -> $appName <> "`",
    "keywords" -> {"igraph", "IGraph/M", "IGraphM"},
    "specialkeywords" -> {},
    "tutorialcollectionlinks" -> {},
    "index" -> True,
    "label" -> "IGraph/M Documentation",
    "language" -> "en",
    "paclet" -> $appName,
    "status" -> "None",
    "summary" -> "IGraph/M is the igraph interface for Mathematica.",
    "synonyms" -> {},
    "tabletags" -> {},
    "title" -> "IGraph/M",
    "titlemodifier" -> "",
    "windowtitle" -> "IGraph/M Documentation",
    "type" -> "Tutorial",
    "uri" -> "IGraphM/tutorial/IGDocumentation"}
};

With[{taggingRules = taggingRules},
  MRun[
    MCode[
      RewriteNotebook[
        NBSetOptions[TaggingRules -> taggingRules] /*
            NBSetOptions[Saveable -> False] /*
            NBRemoveChangeTimes /*
            NBResetWindow /*
            NBDisableSpellCheck /*
            NBDeleteOutputByTag /*
            NBDeleteCellTags["DeleteOutput"] /*
            NBRemoveOptions[{PrivateNotebookOptions, Visible, ShowCellTags}] /*
            NBSetOptions[StyleDefinitions -> NBImport["Stylesheet.nb"]]
      ]["IGDocumentation.nb"];
    ],
    "10.0" (* write notebooks with 10.0 to avoid InsufficientVersionWarning *)
  ]
]

DeleteFile["Stylesheet.nb"]

SetDirectory[".."]
Print["Indexing for 11.1- ..."]
MRun[
  MCode[
    Needs["DocumentationSearch`"];
    indexDir = CreateDirectory["Index"];
    ind = DocumentationSearch`NewDocumentationNotebookIndexer[indexDir];
    DocumentationSearch`AddDocumentationDirectory[ind, "Tutorials"];
    DocumentationSearch`CloseDocumentationNotebookIndexer[ind];
    indexSpellDir = CreateDirectory["SpellIndex"];
    DocumentationSearch`CreateSpellIndex[indexDir, indexSpellDir];
  ],
  "11.1"
]

Print["Indexing for 11.2+ ..."]
MRun[
  MCode[
    Needs["DocumentationSearch`"];
    DocumentationSearch`CreateDocumentationIndex[Directory[], Directory[], "TextSearchIndex", "UseWolframLanguageData" -> False];
  ],
  "11.2"
]
ResetDirectory[] (* .. *)


ResetDirectory[] (* Documentation *)


ResetDirectory[] (* $appTarget *)

Print["\nPacking paclet..."]
PackPaclet[$appTarget]


Print["\nFinished.\n"]


ResetDirectory[] (* $dir *)
