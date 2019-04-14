(* Kernel/init.m *)

(* Reminder: Avoid mentioning any non-System` symbols in this file,
   otherwise they will be created in Global` when the package is loaded. *)

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
  PrependTo[$ContextPath, "IGraphM`"];
  ,
  (* In M11.0 and later, load as normal *)
  Get["IGraphM`IGraphM`"];
]

(* Protect all package symbols *)
SetAttributes[
  Evaluate@Flatten[Names /@ {"IGraphM`*", "IGraphM`Information`*", "IGraphM`Developer`*"}],
  {Protected, ReadProtected}
]
