(* Mathematica Package *)

(* :Copyright: (c) 2016 Szabolcs Horvát *)
(* :License: MIT license, see LICENSE.txt *)

(* This file is read directly with Get in LTemplate.m or LTemplatePrivate.m *)

LTemplate::usage = "LTemplate[name, {LClass[\[Ellipsis]], LClass[\[Ellipsis]], \[Ellipsis]}] represents a library template.";
LClass::usage = "LClass[name, {fun1, fun2, \[Ellipsis]}] represents a class within a template.";
LFun::usage =
    "LFun[name, {arg1, arg2, \[Ellipsis]}, ret] represents a class member function with the given name, argument types and return type.\n" <>
    "LFun[name, LinkObject, LinkObject] represents a function that uses MathLink/WSTP based passing. The shorthand LFun[name, LinkObject] can also be used.";

TranslateTemplate::usage = "TranslateTemplate[template] translates the template into C++ code.";

LoadTemplate::usage = "LoadTemplate[template] loads the library defined by the template. The library must already be compiled.";
UnloadTemplate::usage = "UnloadTemplate[template] attempts to unload the library defined by the template.";

CompileTemplate::usage =
    "CompileTemplate[template] compiles the library defined by the template. Required source files must be present in the current directory.\n" <>
    "CompileTemplate[template, {file1, \[Ellipsis]}] includes additional source files in the compilation.";

FormatTemplate::usage = "FormatTemplate[template] formats the template in an easy to read way.";

ValidTemplateQ::usage = "ValidTemplateQ[template] returns True if the template syntax is valid.";

Make::usage = "Make[class] creates an instance of class.";

LExpressionList::usage = "LExpressionList[class] returns all existing instances of class.";

LClassContext::usage = "LClassContext[] returns the context where class symbols are created.";

LExpressionID::usage = "LExpressionID[name] represents the data type corresponding to LClass[name, \[Ellipsis]] in templates.";

ConfigureLTemplate::usage = "ConfigureLTemplate[options] must be called after loading the LTemplate package privately.";

Begin["`Private`"] (* Begin Private Context *)

(* Private for now, use LFun[name, LinkObject, LinkObject] instead. *)
LOFun::usage =
    "LOFun[name] represents a class member funtion that uses LinkObject for passing and returning arguments.\n" <>
    "It is equivalent to LFun[name, LinkObject, LinkObject].";

(* Mathematica version checks *)

packageAbort[] := (End[]; EndPackage[]; Abort[]) (* Avoid polluting the context path when aborting early. *)

minVersion = {10.0, 0}; (* oldest supported Mathematica version *)
maxVersion = {10.4, 0}; (* latest Mathematica version the package was tested with *)
version    = {$VersionNumber, $ReleaseNumber}
versionString[{major_, release_}] := StringJoin[ToString /@ {NumberForm[major, {Infinity, 1}], ".", release}]

If[Not@OrderedQ[{minVersion, version}],
  Print["LTemplate requires at least Mathematica version " <> versionString[minVersion] <> ".  Aborting."];
  packageAbort[];
]

(* We need to rely on implementation details of SymbolicC, so warn users of yet untested new Mathematica versions. *)
If[Not@OrderedQ[{version, maxVersion}] && Not[$private],
  Print[
    StringTemplate[
      "WARNING: LTemplate has not yet been tested with Mathematica ``.\n" <>
      "The latest supported Mathematica version is ``.\n" <>
      "Please report any issues you find to szhorvat at gmail.com."
    ][versionString[version], versionString[maxVersion]]
  ]
]

(*************** Package configuration ****************)

(* Set up package global variables *)

$packageDirectory = DirectoryName[$InputFileName];
$includeDirectory = FileNameJoin[{$packageDirectory, "IncludeFiles"}];

$messageSymbol (* set by ConfigureLTemplate[] *)


LibraryFunction::noinst = "Managed library expression instance does not exist.";

LTemplate::nofun = "Function `` does not exist.";


Options[ConfigureLTemplate] = { "MessageSymbol" -> LTemplate }

ConfigureLTemplate[opt : OptionsPattern[]] :=
    With[{sym = OptionValue["MessageSymbol"]},
      $messageSymbol = sym;
      sym::info    = "``";
      sym::warning = "``";
      sym::error   = "``";
      sym::assert  = "Assertion failed: ``.";
    ]

LClassContext[] = Context[LTemplate] <> "Classes`";


(***************** SymbolicC extensions *******************)

(* CDeclareAssign[type, var, value] represents
     type var = value;
*)

SymbolicC`Private`IsCExpression[ _CDeclareAssign ] := True

GenerateCode[CDeclareAssign[typeArg_, idArg_, rhs_], opts : OptionsPattern[]] :=
    Module[{type, id},
      type = Flatten[{typeArg}];
      id = Flatten[{idArg}];
      type = Riffle[ Map[ GenerateCode[#, opts] &, type], " "];
      id = Riffle[ Map[ GenerateCode[#, opts] &, id], ", "];
      GenerateCode[CAssign[type <> " " <> id, rhs], opts]
    ]

(* CInlineCode["some code"] will prevent semicolons from being added at the end of "some code" when used in a list *)

GenerateCode[CInlineCode[arg_], opts : OptionsPattern[]] := GenerateCode[arg, opts]

(* CTryCatch[tryCode, catchArg, catchCode] represents
     try { tryCode } catch (catchArg) { catchCode }
*)

GenerateCode[CTryCatch[try_, arg_, catch_], opts: OptionsPattern[]] :=
    Module[{},
      "try " <> GenerateCode[CBlock[try], opts] <> "\n" <>
          "catch (" <> SymbolicC`Private`formatArgument[arg, opts] <> ")\n" <>
          GenerateCode[CBlock[catch], opts]
    ]


(****************** Generic template processing ****************)

(*
 Normalizing a template will:
  - Wrap a bare LClass with LTemplate.  This way a bare LClass can be used as a shorter notation for a single-class template.
  - Convert type names to a canonical form
  - Convert LFun[name, LinkObject, LinkObject] to LOFun[name]
*)

normalizeTypesRules = Dispatch@{
  Verbatim[_Integer] -> Integer,
  Verbatim[_Real] -> Real,
  Verbatim[_Complex] -> Complex,
  Verbatim[True|False] -> "Boolean",
  Verbatim[False|True] -> "Boolean"
};

normalizeFunctionsRules = Dispatch@{
  LFun[name_, LinkObject, LinkObject] :> LOFun[name],
  LFun[name_, LinkObject] :> LOFun[name]
};

normalizeTemplate[c : LClass[name_, funs_]] := normalizeTemplate[LTemplate[name, {c}]]
normalizeTemplate[t : LTemplate[name_, classes_]] := t /. normalizeTypesRules /. normalizeFunctionsRules
normalizeTemplate[t_] := t


ValidTemplateQ::template = "`` is not a valid template. Templates must follow the syntax LTemplate[name, {class1, class2, \[Ellipsis]}].";
ValidTemplateQ::class    = "In ``: `` is not a valid class. Classes must follow the syntax LClass[name, {fun1, fun2, \[Ellipsis]}.";
ValidTemplateQ::fun      = "In ``: `` is not a valid function. Functions must follow the syntax LFun[name, {arg1, arg2, \[Ellipsis]}, ret].";
ValidTemplateQ::string   = "In ``: String expected instead of ``";
ValidTemplateQ::name     = "In ``: `` is not a valid name. Names must start with a letter and may only contain letters and digits.";
ValidTemplateQ::type     = "In ``: `` is not a valid type.";
ValidTemplateQ::rettype  = "In ``: `` is not a valid return type.";
ValidTemplateQ::dupclass = "In ``: Class `` appears more than once.";
ValidTemplateQ::dupfun   = "In ``: Function `` appears more than once.";

ValidTemplateQ[tem_] := validateTemplate@normalizeTemplate[tem]


validateTemplate[tem_] := (Message[ValidTemplateQ::template, tem]; False)
validateTemplate[LTemplate[name_String, classes_List]] :=
    Block[{classlist = {}, location = "template"},
      validateName[name] && (And @@ validateClass /@ classes)
    ]

(* must be called within validateTemplate, uses location, classlist *)
validateClass[class_] := (Message[ValidTemplateQ::class, location, class]; False)
validateClass[LClass[name_, funs_List]] :=
    Block[{funlist = {}, inclass, nameValid},
      nameValid = validateName[name];
      If[MemberQ[classlist, name], Message[ValidTemplateQ::dupclass, location, name]; Return[False]];
      AppendTo[classlist, name];
      inclass = name;
      Block[{location = StringTemplate["class ``"][inclass]},
        nameValid && (And @@ validateFun /@ funs)
      ]
    ]

(* must be called within validateClass, uses location, funlist *)
validateFun[fun_] := (Message[ValidTemplateQ::fun, location, fun]; False)
validateFun[LFun[name_, args_List, ret_]] :=
    Block[{nameValid},
      nameValid = validateName[name];
      If[MemberQ[funlist, name], Message[ValidTemplateQ::dupfun, location, name]; Return[False]];
      AppendTo[funlist, name];
      Block[{location = StringTemplate["class ``, function ``"][inclass, name]},
        nameValid && (And @@ validateType /@ args) && validateReturnType[ret]
      ]
    ]
validateFun[LOFun[name_]] := validateName[name]

(* TODO: Handle other types such as images, sparse arrays, LibraryDataType, etc. *)

numericPattern = Integer|Real|Complex;
sparseArrayPattern = LibraryDataType[SparseArray, numericPattern, PatternSequence[] | _Integer?Positive]; (* disallow SparseArray without explicit element type specification *)
passingMethodPattern = PatternSequence[]|"Shared"|"Manual"|"Constant"|Automatic;

(* must be called within validateTemplate, uses location *)
validateType[numericPattern|"Boolean"|"UTF8String"|sparseArrayPattern|LExpressionID[_String]] := True
validateType[{numericPattern, (_Integer?Positive) | Verbatim[_], passingMethodPattern}] := True
validateType[{sparseArrayPattern, passingMethodPattern}] := True
validateType[type_] := (Message[ValidTemplateQ::type, location, type]; False)

validateReturnType["Void"] := True
validateReturnType[type : LExpressionID[___]] := (Message[ValidTemplateQ::rettype, location, type]; False)
validateReturnType[type_] := validateType[type]

(* must be called within validateTemplate, uses location *)
validateName[name_] := (Message[ValidTemplateQ::string, location, name]; False)
validateName[name_String] :=
    If[StringMatchQ[name, RegularExpression["[a-zA-Z][a-zA-Z0-9]*"]],
      True,
      Message[ValidTemplateQ::name, location, name]; False
    ]



(***********  Translate template to library code  **********)

TranslateTemplate[tem_] :=
    With[{t = normalizeTemplate[tem]},
      If[validateTemplate[t],
        ToCCodeString[transTemplate[t], "Indent" -> 1],
        $Failed
      ]
    ]


libFunArgs = {{"WolframLibraryData", "libData"}, {"mint", "Argc"}, {"MArgument *", "Args"}, {"MArgument", "Res"}};
linkFunArgs = {{"WolframLibraryData", "libData"}, {"MLINK", "mlp"}};
libFunRet  = "extern \"C\" DLLEXPORT int";

excType = "const mma::LibraryError &";
excName = "libErr";

varPrefix = "var";
var[k_] := varPrefix <> IntegerString[k]

includeName[classname_String] := classname <> ".h"

collectionName[classname_String] := classname <> "_collection"

collectionType[classname_String] := "std::map<mint, " <> classname <> " *>"

managerName[classname_String] := classname <> "_manager_fun"

fullyQualifiedSymbolName[sym_Symbol] := Context[sym] <> SymbolName[sym]


setupCollection[classname_String] := {
  CDeclare[collectionType[classname], collectionName[classname]],
  "",
  CFunction["DLLEXPORT void", managerName[classname], {"WolframLibraryData libData", "mbool mode", "mint id"},
    CInlineCode@StringTemplate[ (* TODO: Check if id exists, use assert *)
      "\
if (mode == 0) { // create
  `collection`[id] = new `class`();
} else {  // destroy
  if (`collection`.find(id) == `collection`.end()) {
    libData->Message(\"noinst\");
    return;
  }
  delete `collection`[id];
  `collection`.erase(id);
}\
"][<|"collection" -> collectionName[classname], "class" -> classname|>]
  ],
  "",
  CFunction[libFunRet, classname <> "_get_collection", libFunArgs,
    {
    (* Attention: make sure stuff called here won't throw LibraryError *)
      transRet[
        {Integer, 1},
        CCall["mma::detail::get_collection", collectionName[classname]]
      ],
      CReturn["LIBRARY_NO_ERROR"]
    }
  ],
  "",""
}


registerClassManager[classname_String] :=
    CBlock[{
      "int err",
      StringTemplate[
        "err = (*libData->registerLibraryExpressionManager)(\"`class`\", `manager`)"
      ][<|"class" -> classname, "manager" -> managerName[classname]|>],
      "if (err != LIBRARY_NO_ERROR) return err"
    }]


unregisterClassManager[classname_String] :=
    StringTemplate["(*libData->unregisterLibraryExpressionManager)(\"``\")"][classname]


transTemplate[LTemplate[libname_String, classes_]] :=
    Block[{classlist = {}, classTranslations},
      classTranslations = transClass /@ classes;
      {
        "",
        CInclude["LTemplate.h"],
        CInclude["LTemplateHelpers.h"],
        CInclude /@ includeName /@ classlist,
        "","",

        CInlineCode@"namespace mma {",
        "",
        CDeclare["WolframLibraryData", "libData"],
        "",
        CDefine["LTEMPLATE_MESSAGE_SYMBOL", CString[fullyQualifiedSymbolName[$messageSymbol]]],
        "",
        CInclude["LTemplate.inc"],
        "",
        CInlineCode@"} // namespace mma",

        "","",

        setupCollection /@ classlist,

        CFunction["extern \"C\" DLLEXPORT mint",
          "WolframLibrary_getVersion", {},
          "return WolframLibraryVersion"
        ],
        "",
        CFunction["extern \"C\" DLLEXPORT int",
          "WolframLibrary_initialize", {"WolframLibraryData libData"},
          {
            CAssign["mma::libData", "libData"],
            registerClassManager /@ classlist,
            "return LIBRARY_NO_ERROR"
          }
        ],
        "",
        CFunction["extern \"C\" DLLEXPORT void",
          "WolframLibrary_uninitialize", {"WolframLibraryData libData"},
          {
            unregisterClassManager /@ classlist,
            "return"
          }
        ],
        "","",
        classTranslations
      }
    ]


(* must be called within transTemplate *)
transClass[LClass[classname_String, funs_]] :=
    Block[{},
      AppendTo[classlist, classname];
      transFun[classname] /@ funs
    ]


funName[classname_][name_] := classname <> "_" <> name

transFun[classname_][LFun[name_String, args_List, ret_]] :=
    Block[{index = 0},
      {
        CFunction[libFunRet, funName[classname][name], libFunArgs,
          {
          (* TODO: check Argc is correct, use assert *)
            "const mint id = MArgument_getInteger(Args[0])",
            CInlineCode@StringTemplate[
              "if (`1`.find(id) == `1`.end()) { libData->Message(\"noinst\"); return LIBRARY_FUNCTION_ERROR; }"
            ][collectionName[classname]],
            "",
            CTryCatch[
            (* try *) {
              transArg /@ args,
              "",
              transRet[
                ret,
                CPointerMember[CArray[collectionName[classname], "id"], CCall[name, var /@ Range@Length[args]]]
              ]
            },
            (* catch *) {excType, excName},
              {
                CMember[excName, "report()"],
                CReturn[CMember[excName, "error_code()"]]
              }
            ],
            "",
            CReturn["LIBRARY_NO_ERROR"]
          }
        ],
        "", ""
      }
    ]

transFun[classname_][LOFun[name_String]] :=
    {
      CFunction[libFunRet, funName[classname][name], linkFunArgs,
        {
          CTryCatch[
            (* try *) {
              CInlineCode@StringTemplate[
"
int id;
int args = 2;

if (! MLTestHeadWithArgCount(mlp, \"List\", &args))
  return LIBRARY_FUNCTION_ERROR;
if (! MLGetInteger(mlp, &id))
  return LIBRARY_FUNCTION_ERROR;
if (`collection`.find(id) == `collection`.end()) {
  libData->Message(\"noinst\");
  return LIBRARY_FUNCTION_ERROR;
}
`collection`[id]->`funname`(mlp);
"
              ][<| "collection" -> collectionName[classname], "funname" -> name |>]
            },
            (* catch *) {excType, excName},
            {
              (* TODO: clean up link *)
              CMember[excName, "report()"],
              CReturn[CMember[excName, "error_code()"]]
            }
          ],
          "",
          CReturn["LIBRARY_NO_ERROR"]
        }
      ]
    }

transArg[type_] :=
    Module[{name, cpptype, getfun, setfun},
      index++;
      name = var[index];
      {cpptype, getfun, setfun} = type /. types;
      {
        CDeclareAssign[cpptype, name, StringTemplate["`1`(Args[`2`])"][getfun, index]]
      }
    ]

transRet[type_, value_] :=
    Module[{name = "res", cpptype, getfun, setfun},
      {cpptype, getfun, setfun} = type /. types;
      {
        CDeclareAssign[cpptype, name, value],
        CCall[setfun, {"Res", name}]
      }
    ]

transRet["Void", value_] := value

types = Dispatch@{
  Integer -> {"mint", "MArgument_getInteger", "MArgument_setInteger"},
  Real -> {"double", "MArgument_getReal", "MArgument_setReal"},
  Complex -> {"std::complex<double>", "mma::detail::getComplex", "mma::detail::setComplex"},
  "Boolean" -> {"bool", "MArgument_getBoolean", "MArgument_setBoolean"},
  "UTF8String" -> {"const char *", "mma::detail::getString", "mma::detail::setString"},
  {Integer, __} -> {"mma::IntTensorRef", "mma::detail::getTensor<mint>", "mma::detail::setTensor<mint>"},
  {Real, __} -> {"mma::RealTensorRef", "mma::detail::getTensor<double>", "mma::detail::setTensor<double>"},
  {Complex, __} -> {"mma::ComplexTensorRef", "mma::detail::getTensor< mma::complex_t >", "mma::detail::setTensor< mma::complex_t >"},
  LibraryDataType[SparseArray, Integer, ___] | {LibraryDataType[SparseArray, Integer, ___], ___} -> {"mma::SparseArrayRef<mint>", "mma::detail::getSparseArray<mint>", "mma::detail::setSparseArray<mint>"},
  LibraryDataType[SparseArray, Real, ___] | {LibraryDataType[SparseArray, Real, ___], _} -> {"mma::SparseArrayRef<double>", "mma::detail::getSparseArray<double>", "mma::detail::setSparseArray<double>"},
  LibraryDataType[SparseArray, Complex, ___] | {LibraryDataType[SparseArray, Complex, ___], ___} -> {"mma::SparseArrayRef< mma::complex_t >", "mma::detail::getSparseArray< mma::complex_t >", "mma::detail::setSparseArray< mma::complex_t >"},

  (* This is a special type that translates integer managed expression IDs on the Mathematica side
     into a class reference on the C++ side. It cannot be returned. *)
  LExpressionID[classname_String] :> {classname <> " &", "mma::detail::getObject<" <> classname <> ">(" <> collectionName[classname]  <> ")", ""}
};



(**************** Load library ***************)

(* TODO: Break out loading and compilation into separate files
   This is to make it easy to include them in other projects *)

getCollection (* underlies LClassInstances, the get_collection library function is associated with it in loadClass *)

symName[classname_String] := LClassContext[] <> classname


LoadTemplate[tem_] :=
    With[{t = normalizeTemplate[tem]},
      If[validateTemplate[t],
        Check[loadTemplate[t], $Failed],
        $Failed
      ]
    ]

loadTemplate[tem : LTemplate[libname_String, classes_]] := (
  Quiet@unloadTemplate[tem];
  (* Use FindLibrary to fail early when the library is not found
     Warning: using LibraryLoad instead of FindLibrary here would
     prevent unloading from working at least on OS X.
  *)
  If[FindLibrary[libname] =!= $Failed,
    loadClass[libname] /@ classes,
    Message[LibraryFunction::notfound, libname]
  ];
)

loadClass[libname_][tem : LClass[classname_String, funs_]] := (
  ClearAll[#]& @ symName[classname];
  loadFun[libname, classname] /@ funs;
  With[{sym = Symbol@symName[classname]},
    MessageName[sym, "usage"] = formatTemplate[tem];
    sym[id_Integer][(f_String)[___]] /; (Message[LTemplate::nofun, StringTemplate["``::``"][sym, f]]; False) := $Failed;
    getCollection[sym] = LibraryFunctionLoad[libname, funName[classname]["get_collection"], {}, {Integer, 1}];
  ];
)

loadFun[libname_, classname_][LFun[name_String, args_List, ret_]] :=
    With[{classsym = Symbol@symName[classname]},
      With[{lfun = LibraryFunctionLoad[libname, funName[classname][name], Prepend[args /. loadingTypes, Integer], ret]},
        classsym[id_Integer]@name[arguments___] := lfun[id, arguments]
      ]
    ]

loadFun[libname_, classname_][LOFun[name_String]] :=
    With[{classsym = Symbol@symName[classname]},
      With[{lfun = LibraryFunctionLoad[libname, funName[classname][name], LinkObject, LinkObject]},
        classsym[id_Integer]@name[arguments___] := lfun[id, {arguments}]
      ]
    ]

(* For types that need to be translated to LibraryFunctionLoad compatible forms before loading. *)
loadingTypes = Dispatch@{ LExpressionID[_] -> Integer };


UnloadTemplate[tem_] :=
    With[{t = normalizeTemplate[tem]},
      If[validateTemplate[t],
        unloadTemplate[t],
        $Failed
      ]
    ]

unloadTemplate[LTemplate[libname_String, classes_]] :=
    Module[{res},
      res = LibraryUnload[libname];
      With[{syms = Symbol /@ symName /@ Cases[classes, LClass[name_, __] :> name]},
        ClearAll /@ syms;
        Unset[getCollection[#]]& /@ syms;
      ];
      res
    ]


(* TODO: verify class exists for Make and LExpressionList *)

Make[class_Symbol] := Make@SymbolName[class]
Make[classname_String] := CreateManagedLibraryExpression[classname, Symbol@symName[classname]]


LExpressionList[class_Symbol] := class /@ getCollection[class][]
LExpressionList[classname_String] := LExpressionList@Symbol@symName[classname]


(********************* Compile template ********************)

CompileTemplate[tem_, sources_List, opt : OptionsPattern[CreateLibrary]] :=
    With[{t = normalizeTemplate[tem]},
      If[validateTemplate[t],
        compileTemplate[t, sources, opt],
        $Failed
      ]
    ]

CompileTemplate[tem_, opt : OptionsPattern[CreateLibrary]] := CompileTemplate[tem, {}, opt]

compileTemplate[tem: LTemplate[libname_String, classes_], sources_, opt : OptionsPattern[CreateLibrary]] :=
    Catch[
      Module[{sourcefile, code, includeDirs, classlist, print},
        print[args__] := Apply[Print, Style[#, Darker@Blue]& /@ {args}];

        print["Current directory is: ", Directory[]];
        classlist = Cases[classes, LClass[s_String, __] :> s];
        sourcefile = "LTemplate-" <> libname <> ".cpp";
        If[Not@FileExistsQ[#],
          print["File ", #, " does not exist.  Aborting."]; Throw[$Failed, compileTemplate]
        ]& /@ (# <> ".h"&) /@ classlist;
        print["Unloading library ", libname, " ..."];
        Quiet@LibraryUnload[libname];
        print["Generating library code ..."];
        code = TranslateTemplate[tem];
        If[FileExistsQ[sourcefile], print[sourcefile, " already exists and will be overwritten."]];
        Export[sourcefile, code, "String"];
        print["Compiling library code ..."];
        includeDirs = Flatten[{OptionValue["IncludeDirectories"], $includeDirectory}];
        CreateLibrary[
          AbsoluteFileName /@ Flatten[{sourcefile, sources}], libname,
          "IncludeDirectories" -> includeDirs,
          Sequence@@FilterRules[{opt}, Except["IncludeDirectories"]]]
      ],
      compileTemplate
    ]


(****************** Pretty print a template ********************)

FormatTemplate[template_] :=
    Block[{LFun, LClass, LTemplate},
    (* If the template is invalid, we report errors but we do not abort.
       Pretty-printing is still useful for invalid templates to facilitate finding mistakes.
     *)
      ValidTemplateQ[template];
      formatTemplate@normalizeTemplate[template]
    ]

(* Used in loadClass[] without FormatTemplate wrapper to set usage messages.
   TODO: fix this redundancy.
 *)
formatTemplate[template_] :=
    Block[{LFun, LOFun, LClass, LTemplate},
      With[{tem = template /. normalizeTypesRules /. normalizeFunctionsRules /. {LExpressionID -> "LExpressionID"}},
        LFun[name_, args_, ret_] := StringTemplate["`` ``(``)"][ToString[ret], name, StringJoin@Riffle[ToString /@ args, ", "]];
        LOFun[name_] := StringTemplate["LinkObject ``(LinkObject)"][name];
        LClass[name_, funs_] := StringTemplate["class ``:\n``"][name, StringJoin@Riffle["    " <> ToString[#] & /@ funs, "\n"]];
        LTemplate[name_, classes_] := StringTemplate["template ``\n\n"][name] <> Riffle[ToString /@ classes, "\n\n"];
        tem
      ]
    ]


End[] (* End Private Context *)
