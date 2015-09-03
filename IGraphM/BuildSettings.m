
(*
 * This file contains settings used to compile IGraphM from source.
 * The values below are for my own setup. If you wish to compile from source,
 * you will need to adjust them to match your own system.
 *
 * Note that the igraph (e.g., -ligraph) library must be specified in this file.
 *)

Switch[$OperatingSystem,

  "MacOSX", (* Compilation settings for OS X *)
  $buildSettings = {
    (* IGraphM requires C++11. With OS X's default compiler this is supported 10.9 and up only,
       thus we need to override the default -mmacosx-version-min=10.6 option. *)
    "CompileOptions" -> {"-std=c++11", "-O3", "-flto", "-mmacosx-version-min=10.9", "$HOME/local/lib/libigraph.a  $HOME/local/lib/libgmp.a"},

    (* Set igraph location *)
    "IncludeDirectories" -> {"$HOME/local/include"},
    "LibraryDirectories" -> {"$HOME/local/lib"}

    (*
       If you use MacPorts, the easy way is to install igraph,

       sudo port install igraph

       then use

       "CompileOptions" -> {"-std=c++11", "-mmacosx-version-min=10.9"},
       "Libraries" -> {"igraph"},
       "IncludeDirectories" -> {"/opt/local/include"},
       "LibraryDirectories" -> {"/opt/local/lib"}

     *)
  },

  "Unix", (* Compilation settings for Linux *)
  $buildSettings = {
    "CompileOptions" -> {"-std=c++11", "-O3", "$HOME/local/lib/libigraph.a  $HOME/local/lib/libgmp.a"},

    (* Set igraph location *)
    "IncludeDirectories" -> {"$HOME/local/include"},
    "LibraryDirectories" -> {"$HOME/local/lib"}
  }

]
