(* ::Package:: *)

(* Mathematica Source File  *)
(* Created by Mathematica Plugin for IntelliJ IDEA *)
(* :Author: szhorvat *)
(* :Date: 2016-12-02 *)

(* This script uses https://github.com/szhorvat/PackageTools *)

(* This script requires the following versions of Mathematica:
 * 11.0 for evaluating the documentation notebook and writing it in a 11.0-compatible form
 * 11.1 for creating old-style documentation indexes
 * 11.2 for creating new-style documentation indexes
 *)

$appName = "IGraphM";
$LTemplateRepo = "~/Repos/LTemplate"; (* location of the LTemplate repository *)

Print["\nScript running in Mathematica ", $Version, ".\n"];

(* Prints an error message and aborts the script *)
printAbort[str_] := (Print["ABORTING: ", Style[str, Red, Bold]]; Quit[])
If[$VersionNumber < 11.0, printAbort["Mathematica 11.0 or later required."]]


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

(* Include license file. *)
CopyFile["LICENSE.txt", FileNameJoin[{$appTarget, "LICENSE.txt"}]]

SetDirectory[$appTarget]


Print["Loading IGraphM paclet..."]
PacletDirectoryAdd[$dir]
Needs[$appName <> "`"]


Print["Making source replacements..."]

replaceLines[repl_][file_] :=
    Module[{source},
      source = Import[file, "String"];
      source = StringReplace[source, repl];
      Export[file, source, "String"]
    ]


replaceLines[{
    "Get[\"LTemplate`LTemplatePrivate`\"]; (* REPLACE-TAG *)" -> "Get[\"IGraphM`LTemplate`LTemplatePrivateNoCompile`\"]",
    "ConfigureLTemplate[\"MessageSymbol\" -> IGraphM, \"LazyLoading\" -> False] (* REPLACE-TAG *)"
        -> "ConfigureLTemplate[\"MessageSymbol\" -> IGraphM, \"LazyLoading\" -> True]"
  }]@FileNameJoin[{$appTarget, "LTemplate.m"}]


igDataCompletionRepl = StringTemplate["addCompletion[IGData, ``];"]@
    ToString[{Join[Keys[IGraphM`IGData`PackagePrivate`$igDataCategories],
      Select[Keys[IGraphM`IGData`PackagePrivate`$igData], StringQ]]},
      InputForm
    ]

replaceLines[{
  "addCompletion[IGData, {Join[Keys[$igDataCategories], Select[Keys[$igData], StringQ]]}]; (* REPLACE-TAG *)" -> igDataCompletionRepl
  }]@FileNameJoin[{$appTarget, "IGData.m"}]


igLatticeMeshCompletionRepl = StringTemplate["addCompletion[IGLatticeMesh, ``];"]@
    ToString[{IGLatticeMesh[]},
      InputForm
    ]

replaceLines[{
  "addCompletion[IGLatticeMesh, {IGLatticeMesh[]}]; (* REPLACE-TAG *)" -> igLatticeMeshCompletionRepl
  }]@FileNameJoin[{$appTarget, "Meshes.m"}]


applyTemplate[templateData_][file_] :=
    Module[{source, template},
      source = Import[file, "String"];
      template = StringTemplate[source, Delimiters -> {"%%", "%%", "<%", "%>"}];
      source = template[templateData];
      Export[file, source, "String"];
    ]

versionData = <|
  "version" -> Lookup[PacletInformation[$appName], "Version"],
  "mathversion" -> Lookup[PacletInformation[$appName], "WolframVersion"],
  "date" -> DateString[{"MonthName", " ", "DayShort", ", ", "Year"}]
|>

applyTemplate[versionData] /@ FileNames["*.m", $appTarget]


(* Replace unicode characters with their Mathematica FullForm, to ensure that
   source file work identically regardless of the value of $CharacterEncoding *)
Print["\nRe-encoding source files as ASCII"]
AddPath["PackageTools"]
Needs["PackageTools`"]
With[{$appTarget = $appTarget, files = AbsoluteFileName /@ FileNames["*.m", $appTarget]},
  MRun[
    MCode[
      UsingFrontEnd[
        NotebookSave@NotebookOpen[#]& /@ files
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

Print["Evaluating documentation notebook..."]
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

Print["Rewriting documentation notebook..."]

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
        NBDeleteCellByTag /*
        NBDeleteOutputByTag /*
        NBDeleteCellTags["DeleteOutput"] /*
        NBRemoveOptions[{PrivateNotebookOptions, Visible, ShowCellTags}] /*
        NBSetOptions[StyleDefinitions -> NBImport["Stylesheet.nb"]]
      ]["IGDocumentation.nb"];
    ],
    "11.0" (* write notebooks with 11.0 to avoid InsufficientVersionWarning *)
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

(* Warning: Run this in 11.2 only. This puts the index where 11.2 looks for it.
   11.3+ produces incompatible indexes that break search in 11.2. *)
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
