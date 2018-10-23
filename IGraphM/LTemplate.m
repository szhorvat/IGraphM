
(* This is a traditional package that uses BeginPackage instead of Package.
   It loads LTemplate into the correct context. *)

BeginPackage["IGraphM`"]

(* NOTE: replaced in build script. Remember to update build script if editing these two lines. *)
Get["LTemplate`LTemplatePrivate`"];
ConfigureLTemplate["MessageSymbol" -> IGraphM, "LazyLoading" -> False]

EndPackage[]