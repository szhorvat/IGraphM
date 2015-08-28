(* Mathematica Package  *)
(* Created by IntelliJ IDEA and wlplugin.halirutan.de *)

(* :Title: IGraph/M    *)
(* :Context: IGraphM`  *)
(* :Author: szhorvat   *)
(* :Date: 2015-08-28   *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: 10.0 *)
(* :Copyright: (c) 2015 Szabolcs Horvát *)
(* :Keywords: igraph, graphs, networks, LibraryLink *)
(* :Discussion: igraph interface for Mathematica, see http://igraph.org *)

BeginPackage["IGraphM`"]

Get["LTemplate`LTemplatePrivate`"]

IGVersion::usage = "";


Begin["`Private`"]

$packageDirectory = DirectoryName[$InputFileName];
$libraryDirectory = FileNameJoin[{$packageDirectory, "LibraryResources", $SystemID}];
$sourceDirectory  = FileNameJoin[{$packageDirectory, "LibraryResources", "Source"}];


If[Not@MemberQ[$LibraryPath, $libraryDirectory],
  AppendTo[$LibraryPath, $libraryDirectory]
]


template = LTemplate["IGraphM",
  {
    LClass["IGlobal",
      {
        LFun["init", {}, "Void"],
        LFun["version", {}, "UTF8String"]
      }
    ],

    LClass["IG",
      {
        LFun["fromEdgeList", {{Real, 2}, Integer, True | False}, "Void"],
        LFun["edgeCount", {}, Integer],
        LFun["vertexCount", {}, Integer],
        LFun["edgeList", {}, {Real, 2}],
        LFun["directedQ", {}, True | False],
        LFun["dagQ", {}, True | False],
        LFun["simpleQ", {}, True | False],
        LFun["connectedQ", {}, True | False]
      }
    ]
  }
];


recompileLibrary[] :=
  Block[{$CCompiler},
    If[Not@DirectoryQ[$libraryDirectory],
      CreateDirectory[$libraryDirectory]
    ];
    $CCompiler = {
      "Compiler" -> CCompilerDriver`GenericCCompiler`GenericCCompiler,
      "CompilerName" -> "g++-mp-5",
      "CompilerInstallation" -> "/opt/local/bin",
      "SystemCompileOptions" -> {"-std=c++11", "-m64", "-fPIC", "-O2", "-framework Foundation", "-framework \"mathlink\""}
    };
    SetDirectory[$sourceDirectory];
    CompileTemplate[template,
      "IncludeDirectories" -> {"/opt/local/include"},
      "LibraryDirectories" -> {"/opt/local/lib"},
      "CompileOptions" -> {"-ligraph"},
      "ShellCommandFunction" -> Print, "ShellOutputFunction" -> Print,
      "TargetDirectory" -> $libraryDirectory
    ];
    ResetDirectory[]
  ]


If[LoadTemplate[template] === $Failed,
  Print["Loading failed, trying to recompile ..."];
  recompileLibrary[];
  If[LoadTemplate[template] === $Failed,
    Print[Style["Cannot load or compile library.  Aborting.", Red]];
    Abort[]
  ]
]


igraphGlobal = Make["IGlobal"]; (* there should only be a single object of this type *)
igraphGlobal@"init"[] (* Run initialization *)


IGVersion[] := igraphGlobal@"version"[]


End[] (* `Private` *)

EndPackage[]