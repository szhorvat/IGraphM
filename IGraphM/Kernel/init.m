
(* Mathematica version check *)
If[Not@OrderedQ[{10.0, 2}, {$VersionNumber, $ReleaseNumber}],
  Print["IGraph/M requires Mathematica 10.0.2 or later.  Aborting."];
  Abort[]
]

Unprotect["IGraphM`*", "IGraphM`Developer`*", "IGraphM`Information`*"];

Get["IGraphM`IGraphM`"]

With[{syms = Join @@ Names /@ {"IGraphM`*", "IGraphM`Information`*", "IGraphM`Developer`*"}},
  SetAttributes[syms, {Protected, ReadProtected}]
]
