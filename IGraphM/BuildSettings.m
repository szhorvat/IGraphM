
(*
 * This file contains settings used to compile IGraphM from source.
 * The values below are for my own setup. If you wish to compile from source,
 * you will need to adjust them to match your own system.
 *
 * Note that the igraph library (e.g., -ligraph or libigraph.a) must be specified within this file.
 *)

Switch[$OperatingSystem,

  "MacOSX", (* Compilation settings for OS X *)
  $buildSettings = {
    (*"Compiler" -> CCompilerDriver`GenericCCompiler`GenericCCompiler,
    "CompilerInstallation" -> "/opt/local/bin",
    "CompilerName" -> "clang++-mp-7.0",
    "SystemCompileOptions" -> "-O3 -m64 -fPIC -framework Foundation -framework mathlink",*)

    (* IGraphM requires C++11. With OS X's default compiler this is supported 10.9 and up only,
       thus we need to override the default -mmacosx-version-min=10.6 option. *)
    "CompileOptions" -> {"-mmacosx-version-min=10.9", "-fsanitize=address", " -fsanitize-recover=address", "-fno-omit-frame-pointer"},

    (* Statically link the igraph library *)
    "ExtraObjectFiles" -> {"$HOME/local/lib/libigraph.a", "$HOME/local/lib/libgmp.a", "$HOME/local/lib/libglpk.a", "$HOME/local/lib/libemon.a"},

    (* Set igraph location *)
    "IncludeDirectories" -> {"$HOME/local/include"},
    "LibraryDirectories" -> {"$HOME/local/lib"}
  },

  "Unix", (* Compilation settings for Linux *)
  $buildSettings = {
    "CompileOptions" -> {
      (* Compile with -static-libgcc on non-RPi Linux for better compatibility with various distros *)
      If[$SystemID =!= "Linux-ARM",
        Unevaluated@Sequence["-static-libgcc", "-D_GLIBCXX_USE_CXX11_ABI=0", "-flto"],
        Unevaluated@Sequence[]
      ]
    },

    (* Statically link the igraph library *)
    "ExtraObjectFiles" -> {"$HOME/local/lib/libigraph.a", "$HOME/local/lib/libgmp.a", "$HOME/local/lib/libglpk.a", "$HOME/local/lib/libemon.a"},

    (* Set igraph location *)
    "IncludeDirectories" -> {"$HOME/local/include"},
    "LibraryDirectories" -> {"$HOME/local/lib"}
  },

  "Windows",
  $buildSettings = { 
    "CompileOptions" -> {"/EHsc", "/wd4244", "/DNOMINMAX"},
    "IncludeDirectories" -> {"C:\\msys64\\home\\%USERNAME%\\local\\include", "C:\\msys64\\home\\%USERNAME%\\lemon\\include"},
    "ExtraObjectFiles" -> {"C:\\msys64\\home\\%USERNAME%\\local\\lib\\libigraph.dll.a", "C:\\msys64\\home\\%USERNAME%\\lemon\\lib\\lemon.lib"}
  }
]
