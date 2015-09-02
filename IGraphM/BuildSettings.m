
(*
 * This file contains settings used to compile IGraphM from source.
 * The values below are for my own setup. If you wish to compile from source,
 * you will need to adjust them to match your own system.
 *)

Switch[$OperatingSystem,

  "MacOSX", (* Compilation settings for OS X *)
  $buildSettings = {
    (* IGraphM requires C++11. With OS X's default compiler this is supported 10.9 and up only,
       thus we need to override the default -mmacosx-version-min=10.6 option. *)
    "CompileOptions" -> {"-std=c++11", "-mmacosx-version-min=10.9"},

    (* Set igraph location *)
    "IncludeDirectories" -> {"$HOME/local/include"},
    "LibraryDirectories" -> {"$HOME/local/lib"}
  },

  "Unix", (* Compilations settings for Linux *)
  $buildSettings = {
    "CompileOptions" -> {"-std=c++11"}
  }

]
