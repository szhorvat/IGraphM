(* Kernel/init.m *)

(* Reminder: Avoid mentioning any non-System` symbols in this file,
   otherwise they will be created in Global` when the package is loaded. *)

If[$CloudEvaluation,
  Print["IGraph/M cannot run in the Wolfram Cloud because the cloud does not support LibraryLink.  Aborting."];
  Abort[]
]

(* TODO!
   In M12.1.0, the FE will hang if Needs["IGraphM`"] is evaluated within a certain amount of time
   after kernel startup. The hung FE can be recovered by killing the kernel.

   Here we use the workaround of delaying package loading with a Pause[] if kernel startup
   happened recently.

   What we know so far about the hang:
     - It is not caused by addCompletion[]
     - Disabling SyntaxInformation using (Unprotect[SyntaxInformation]; SyntaxInformation = dummy) avoids the hang.
 *)
If[$VersionNumber == 12.1 && SessionTime[] < 3,
  Pause[1.5]
]

(* Mathematica version check *)
If[Not@OrderedQ[{10.0, 2}, {$VersionNumber, $ReleaseNumber}],
  Print["IGraph/M requires Mathematica 10.0.2 or later.  Aborting."];
  Abort[]
]

(* Special version check on the Raspberry Pi. *)
(* The Wolfram Engine is free on the Raspberry Pi. We only support the latest release at any time. *)
If[$SystemID === "Linux-ARM" && $VersionNumber < 12.0,
  Print["On the Raspberry Pi, IGraph/M only supports the Wolfram Engine version 12.0 or later.  Aborting."];
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

(* Show welcome message by returning it from the Needs/Get call *)
Column@{
  "IGraph/M " <> IGraphM`IGraphM`PackagePrivate`$packageVersion,
  If[$Notebooks,
    "Evaluate \!\(\*ButtonBox[\"IGDocumentation[]\",BaseStyle->\"Link\",ButtonData->\"paclet:IGraphM\"]\) to get started.",
    "Evaluate IGDocumentation[] to get started."
  ]
}
