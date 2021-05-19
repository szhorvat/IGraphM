
(*
 * This file contains settings used to compile IGraphM from source.
 * The values below are for my own setup. If you wish to compile from source,
 * you will need to adjust them to match your own system.
 *
 * Note that the igraph library (e.g., -ligraph or libigraph.a) must be specified within this file.
 *)

$useAddressSanitizer = False;

Switch[$OperatingSystem,

  "MacOSX", (* Compilation settings for OS X *)
  $buildSettings = {
    "CompileOptions" -> {
      If[$useAddressSanitizer,
        Unevaluated@Sequence["-g", "-Og", "-fno-omit-frame-pointer", "-fsanitize=address", "-fsanitize=undefined"],
        Unevaluated@Sequence["-flto=thin", "-O3"]
      ],
      "-fvisibility=hidden",
      "-framework Accelerate", (* for BLAS and LAPACK *)
      "-mmacosx-version-min=10.9", (* earliest supported macOS version---required for C++11 *)
      With[{res = Quiet@RunProcess[{"xcrun", "--sdk", "macosx", "--show-sdk-path"}]},
        (* If the SDK version can be determined, use it. igraph will be compiled with this SDK by default,
           and if we don't match it for IGraph/M, the linker will give errors. *)
        If[Not@FailureQ[res] && res["ExitCode"] == 0,
          "-isysroot " <> StringTrim[res["StandardOutput"]],
          Unevaluated@Sequence[]
        ]
      ]
    },

    (* Statically link the igraph library *)
    "ExtraObjectFiles" -> {"$HOME/local/lib/libigraph.a", "$HOME/local/lib/libglpk.a", "$HOME/local/lib/libemon.a"},

    (* Set igraph location *)
    "IncludeDirectories" -> {"$HOME/local/include"},
    "LibraryDirectories" -> {"$HOME/local/lib"}
  },

  "Unix", (* Compilation settings for Linux *)
  $buildSettings = {
    "CompileOptions" -> {
      If[$SystemID =!= "Linux-ARM",
        (* Compile with -static-libgcc on non-RPi Linux for better compatibility with older distros *)
        (* Do not use -flto at this point when compiling on Ubuntu 16.04 as it leads to crashes when igraph returns an error *)
        Unevaluated@Sequence["-static-libgcc", "-D_GLIBCXX_USE_CXX11_ABI=0"(*, "-flto"*)],
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
  With[{$depDir = "%USERPROFILE%\\Documents\\Packages"},
    $buildSettings = { 
      "CompileOptions" -> {"/EHsc", "/GL", "/wd4244", "/DNOMINMAX", "/DIGRAPH_STATIC"},
      "IncludeDirectories" -> {$depDir <> "\\igraph\\include", $depDir <> "\\lemon\\include"},
      "ExtraObjectFiles" -> {$depDir <> "\\igraph\\lib\\igraph.lib", $depDir <> "\\lemon\\lib\\lemon.lib", $depDir <> "\\glpk-5.0\\w64\\glpk.lib"}
    }
  ]
]
