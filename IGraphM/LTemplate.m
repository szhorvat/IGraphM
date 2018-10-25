
(* This is a traditional package that uses BeginPackage instead of Package.
   It loads LTemplate into the correct context. *)

BeginPackage["IGraphM`"]

Get["LTemplate`LTemplatePrivate`"]; (* REPLACE-TAG *)
ConfigureLTemplate["MessageSymbol" -> IGraphM, "LazyLoading" -> False] (* REPLACE-TAG *)

EndPackage[]