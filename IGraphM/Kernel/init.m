
(* Mathematica version check *)
If[Not@OrderedQ[{10.0, 2}, {$VersionNumber, $ReleaseNumber}],
  Print["IGraph/M requires Mathematica 10.0.2 or later.  Aborting."];
  Abort[]
]


Unprotect["IGraphM`*", "IGraphM`Developer`*", "IGraphM`Information`*"];

If[$VersionNumber < 11.0,
  (* Workaround for the $Context and $ContextPath not being set up during loading in new style packages
     in Mathematica versions earlier than 11.0.
     See also igContextSetup[], defined in IGraphM.m *)
  Block[{$Context, $ContextPath},
    Get["IGraphM`IGraphM`"];
  ];
  PrependTo[$ContextPath, "IGraphM`"]
  ,
  Get["IGraphM`IGraphM`"];
]

With[{IGraphM`IGraphM`PackagePrivate`syms = Join @@ Names /@ {"IGraphM`*", "IGraphM`Information`*", "IGraphM`Developer`*"}},
  SetAttributes[IGraphM`IGraphM`PackagePrivate`syms, {Protected, ReadProtected}]
];
