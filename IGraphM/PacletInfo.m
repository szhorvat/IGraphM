(* ::Package:: *)

(* Paclet Info File *)

Paclet[
    Name -> "IGraphM",
    Version -> "0.5.90",
    MathematicaVersion -> "11.0+",
    Description -> "IGraph/M \[Dash] the igraph interface for Mathematica.",
    Creator -> "Szabolcs Horvát <szhorvat@gmail.com>",
    URL -> "http://szhorvat.net/mathematica/IGraphM",
    Thumbnail -> "Logo.png",
    "Icon" -> "Logo.png",
    "Keywords" -> {"igraph", "graph theory", "network analysis"},
    SystemID -> {"MacOSX-x86-64", "MacOSX-ARM4", "Windows-x86-64", "Linux-x86-64"(*, "Linux-ARM"()*)},
    Extensions -> 
        {    
            {"Kernel", Root -> ".", Context -> "IGraphM`"},
            {"LibraryLink"},
            {"Documentation", MainPage -> "Tutorials/IGDocumentation"}
        }
]
