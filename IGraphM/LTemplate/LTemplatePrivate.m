(* Mathematica Package *)

(* :Copyright: (c) 2016 Szabolcs Horvát *)
(* :License: MIT license, see LICENSE.txt *)

(*
 * To include LTemplate privately in another package, load LTemplatePrivate.m using Get[],
 * then immediately call ConfigureLTemplate[].
 *)

BeginPackage["`LTemplate`", {"SymbolicC`", "CCodeGenerator`", "CCompilerDriver`"}]

`Private`$private = True;
Quiet[
  Get@FileNameJoin[{DirectoryName[$InputFileName], "LTemplateInner.m"}],
  General::shdw (* suppress false shadowing warnings if public LTemplate was loaded first *)
]

EndPackage[]
