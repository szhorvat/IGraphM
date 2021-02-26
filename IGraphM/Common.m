(* Mathematica Package *)
(* Created by Mathematica Plugin for IntelliJ IDEA, see http://wlplugin.halirutan.de/ *)

(* :Author: szhorvat *)
(* :Date: 2016-06-12 *)
(* :Copyright: (c) 2016-2020 Szabolcs Horvát *)


Package["IGraphM`"]


(********************************************)
(***** Common (package scope) functions *****)
(********************************************)


PackageScope["igGraphQ"]
igGraphQ::usage = "igGraphQ[g] checks if g is an igraph-compatible graph.";
igGraphQ = GraphQ[#] && If[MixedGraphQ[#], Message[IGraphM::mixed]; False, True] &;


PackageScope["igDirectedQ"]
igDirectedQ::usage = "igDirectedQ[g] checks if g is a directed graph. Empty graphs are considered undirected.";
igDirectedQ[graph_] := DirectedGraphQ[graph] && Not@EmptyGraphQ[graph]


PackageScope["applyGraphOpt"]
applyGraphOpt::usage = "applyGraphOpt[options][graph] applies the given options to graph.";
applyGraphOpt[opt___][graph_] := Graph[graph, Sequence@@FilterRules[{opt}, Options[Graph]]]


PackageScope["applyGraphOpt3D"]
applyGraphOpt3D::usage = "applyGraphOpt3D[options][graph] applies the given options to graph using Graph3D.";
applyGraphOpt3D[opt___][graph_] := Graph3D[graph, Sequence@@FilterRules[{opt}, Options[Graph3D]]]


PackageScope["zeroDiagonal"]
zeroDiagonal::usage = "zeroDiagonal[mat] replaces the diagonal of a matrix with zeros.";
zeroDiagonal[mat_] := UpperTriangularize[mat, 1] + LowerTriangularize[mat, -1]


PackageScope["adjacencyGraph"]
(* Fast version of directed or undirected adjacency matrix -> graph conversion.
   The matrix is not checked to be symmetric in the undirected case. *)
adjacencyGraph::usage = "adjacencyGraph[vertices, sparseAM, directed]";
adjacencyGraph[vs_, sa_, True] := Graph[vs, {sa, Null}]
adjacencyGraph[vs_, sa_, False] := Graph[vs, {Null, sa}]


PackageScope["removeSelfLoops"]
removeSelfLoops::usage = "removeSelfLoops[graph] removes any self loops from graph.";
removeSelfLoops[g_?LoopFreeGraphQ] := g (* catches empty case *)
removeSelfLoops[g_] := adjacencyGraph[VertexList[g], zeroDiagonal@AdjacencyMatrix[g], DirectedGraphQ[g]]


PackageScope["removeMultiEdges"]
removeMultiEdges::usage = "removeMultiEdges[graph] removes any multi-edges from graph.";
removeMultiEdges[g_?MultigraphQ] := adjacencyGraph[VertexList[g], Unitize@AdjacencyMatrix[g], DirectedGraphQ[g]]
removeMultiEdges[g_] := g


PackageScope["transformGraphOptions"]

(* Use 'dummy' instead of 'usage' to prevent IntelliJ from making the symbol known across files *)
$graphLink::dummy = "$graphLink is a loopback link used to convert atomic graphs to a compound form.";

transformGraphOptions::usage = "transformGraphOptions[fun][graph] applies fun to the list of options stored in graph.";
transformGraphOptions[fun_][g_?GraphQ] :=
    (
      If[Not@MemberQ[Links[], $graphLink],
        $graphLink = LinkCreate[LinkMode -> Loopback];
      ];
      With[
        {
          expr = AbortProtect[
            LinkWrite[$graphLink, g];
            LinkRead[$graphLink, Hold]
          ]
        },
        Replace[expr, Hold@Graph[v_, e_, opt : _ : {}, rest___] :> Graph[v, e, fun[opt], rest]]
      ]
    )


PackageScope["ruleQ"]
ruleQ::usage = "ruleQ[expr] gives True if expr is a Rule or RuleDelayed, False otherwise.";
ruleQ[_Rule | _RuleDelayed] = True;
ruleQ[_] = False;


PackageScope["partitionRagged"]
partitionRagged::usage = "partitionRagged[list, lengths] partitions list into parts of the given lengths.";
If[$VersionNumber >= 11.2,
  partitionRagged = TakeList,
  partitionRagged[v_List, l_?VectorQ] := MapThread[Take[v, {#1, #2}] &, With[{a = Accumulate[l]}, {a - l + 1, a}]]
]

(* Retrieving edge or vertex weights this way is much faster than using PropertyValue *)
PackageScope["igEdgeWeights"]
igEdgeWeights::usage = "igEdgeWeights[graph] returns the edge weight vector of graph.";
igEdgeWeights = GraphComputation`WeightValues;

PackageScope["igVertexWeights"]
igVertexWeights ::usage = "igVertexWeights[graph] returns the vertex weight vector of graph.";
igVertexWeights = GraphComputation`WeightVector;


(***** Tools for error handling *****)

(* TODO: we need a better error handling framework.
 * The message tag should be encapsulated in the thrown object.
 * catch[] should have a version that takes a head and reports the message under that head (head::tag)
 *)

(* Use 'dummy' instead of 'usage' to prevent IntelliJ from making the symbol known across files *)
igTag::dummy = "igTag is a private tag for Throw/Catch within IGraphM.";

PackageScope["throw"]
throw::usage = "throw[val]";
throw[val_] := Throw[val, igTag]

PackageScope["catch"]
catch::usage = "catch[expr]";
SetAttributes[catch, HoldFirst]
catch[expr_] := Catch[expr, igTag]

PackageScope["check"]
check::usage = "check[val]";
check[val_LibraryFunctionError] := throw[$Failed] (* this was originally throw[val] *)
check[$Failed] := throw[$Failed]
check[HoldPattern[LibraryFunction[___][___]]] := throw[$Failed]
check[val_] := val

PackageScope["sck"]
sck::usage = "sck[val]";
sck[HoldPattern[LibraryFunction[___][___]]] := $Failed
sck[val_LibraryFunctionError] := $Failed (* this was originally val *)
sck[val_] := val


(***** Functions for argument checking  *****)

(* Note: VectorQ, MatrixQ, etc. are optimized with certain second arguments. For example,
 * VectorQ[#, NumericQ]& will immediately return True for a packed array without evaluating
 * NumericQ for each array argument separately. Below it is noted which of these second arguments
 * the check is fast for. *)

PackageScope["nonNegIntVecQ"]
nonNegIntVecQ::usage = "nonNegIntVecQ[vec]";
nonNegIntVecQ = VectorQ[#, Internal`NonNegativeMachineIntegerQ]&; (* verified fast M11.1+ *)

PackageScope["posIntVecQ"]
posIntVecQ::usage = "posIntVecQ[vec]";
posIntVecQ = VectorQ[#, Internal`PositiveMachineIntegerQ]&; (* verified fast M11.1+ *)

PackageScope["intVecQ"]
intVecQ::usage = "intVecQ[]";
intVecQ =
    If[$VersionNumber < 11.0,
      VectorQ[#, IntegerQ]&, (* In M10.4 and earlier VectorQ[{}, Developer`MachineIntegerQ] returns False. M11.0+ is fine. *)
      VectorQ[#, Developer`MachineIntegerQ]& (* verified fast *)
    ];

PackageScope["intMatQ"]
intMatQ::usage = "intMatQ[mat]";
intMatQ =
    If[$VersionNumber < 11.0,
      MatrixQ[#, IntegerQ]&, (* In M10.4 and earlier MatrixQ[{{}}, Developer`MachineIntegerQ] returns False. M11.0+ is fine. *)
      MatrixQ[#, Developer`MachineIntegerQ]& (* verified fast *)
    ];

PackageScope["positiveNumericQ"]
positiveNumericQ::usage = "positiveNumericQ[num]";
positiveNumericQ = NumericQ[#] && TrueQ@Positive[#]&;

PackageScope["nonNegativeNumericQ"]
nonNegativeNumericQ::usage = "nonNegativeNumericQ[num]";
nonNegativeNumericQ = NumericQ[#] && TrueQ@NonNegative[#]&;

PackageScope["positiveVecQ"]
positiveVecQ::usage = "positiveVecQ[vec]";
positiveVecQ = VectorQ[#, Positive]&; (* NOT fast *)

PackageScope["nonNegVecQ"]
nonNegVecQ::usage = "nonNegVecQ[vec]";
nonNegVecQ = VectorQ[#, NonNegative]&; (* NOT fast *)

PackageScope["emptyArrayQ"]
emptyArrayQ::usage = "emptyArrayQ[arr]";
emptyArrayQ[arr_] := MemberQ[Dimensions[arr], 0]

PackageScope["positiveOrInfQ"]
positiveOrInfQ::usage = "positiveOrInfQ[val]";
positiveOrInfQ = TrueQ@Positive[#]&;

(* Replace Infinity by 0 *)
PackageScope["infToZero"]
infToZero::usage = "infToZero[arg] returns 0 if arg === Infinity.";
infToZero[arg_] := Replace[arg, Infinity -> 0]

(* Replace Infinity by -1 *)
PackageScope["infToNeg"]
infToNeg::usage = "infToNeg[arg] returns 0 if arg === Infinity.";
infToNeg[arg_] := Replace[arg, Infinity -> -1]

(* Temporarily disable floating point exception checking in LibraryLink *)
PackageScope["expectInfNaN"]
expectInfNaN::usage = "expectInfNaN[expr] evaluates expr with settings that allow LibraryLink functions to return Infinity or Indeterminate.";
If[$VersionNumber >= 12.2,
  SetAttributes[expectInfNaN, HoldAll];
  expectInfNaN[expr_] :=
      With[{llo = SystemOptions["LibraryLinkOptions" -> "TestFloatingPointExceptions"]},
        Internal`WithLocalSettings[
          SetSystemOptions["LibraryLinkOptions" -> "TestFloatingPointExceptions" -> False],
          expr,
          SetSystemOptions[llo]
        ]
      ]
  ,
  expectInfNaN[expr_] := expr
]

(* Unpack array containing infinities or indeterminates *)
(* TODO: Test on all platforms that unpacking such arrays produces usable Infinity and Indeterminate *)
PackageScope["fixInfNaN"]
fixInfNaN::usage = "fixInfNaN[array] unpacks array if it contains Inf or NaN.";
fixInfNaN[arr_?Developer`PackedArrayQ] := If[igraphGlobal@"infOrNanQ"[arr], Developer`FromPackedArray[arr], arr];
fixInfNaN[arr_] := arr

(***** Workarounds for old versions missing some functions *****)

PackageScope["keyValueMap"]
keyValueMap::usage = "keyValueMap[fun, asc] is a v10.0-compatible replacement for KeyValueMap.";
If[$VersionNumber >= 10.1,
  keyValueMap = KeyValueMap,
  keyValueMap[f_, asc_] := f @@@ Normal[asc]
]


PackageScope["connectedGraphComponents"]
connectedGraphComponents::usage = "connectedGraphComponents[graph]";
If[$VersionNumber >= 10.4,
  connectedGraphComponents = WeaklyConnectedGraphComponents,
  connectedGraphComponents[g_] := Subgraph[g, #]& /@ WeaklyConnectedComponents[g]
]


PackageScope["circularEmbedding"]
circularEmbedding::usage = "circularEmbedding is a GraphLayout.";
If[$VersionNumber < 12.0,
  circularEmbedding = "CircularEmbedding",
  circularEmbedding = {"CircularEmbedding", "OptimalOrder" -> False}
]


PackageScope["canonicalEdgeBlock"]
canonicalEdgeBlock::usage = "canonicalEdgeBlock[expr] evaluates expression while making sure that all UndirectedEdge expressions inside are ordered canonically.";
SetAttributes[canonicalEdgeBlock, HoldAll]

(* In M12.1 and later, UndirectedEdge can have 3 arguments, so we cannot canonicalize simply with Orderless. *)
(* TODO  The workaround /; Not@OrderedQ[{a, b}] is 10x slower than Orderless! *)
If[$VersionNumber >= 12.1,
  canonicalEdgeBlock[expr_] :=
      Internal`InheritedBlock[{UndirectedEdge},
        Unprotect[UndirectedEdge];
        UndirectedEdge[a_, b_, rest___] /; Not@OrderedQ[{a, b}] := UndirectedEdge[b, a, rest];
        expr
      ]
  ,
  canonicalEdgeBlock[expr_] :=
      Internal`InheritedBlock[{UndirectedEdge},
        SetAttributes[UndirectedEdge, Orderless];
        expr
      ]
]


(* In M12.1 and later, custom graph properties are stored in the AnnotationRules option, not in the Properties option.
   PropOptName is defined  to be the correct property-holding option name in each version. *)
PackageScope["PropOptName"]
PropOptName::usage = "PropOptName is AnnotationRules in M >= 12.1 and Properties in M < 12.0.";
If[$VersionNumber >= 12.1,
  PropOptName = AnnotationRules,
  PropOptName = Properties
]


PackageScope["encodeNeighborMode"]
encodeNeighborMode::usage = "encodeNeighborMode[mode] translates Out, In, All to 1, 2, 3. To be used with the igNeighborMode() C++ function.";
encodeNeighborMode[mode_] := Lookup[<|"Out" -> 1, "In" ->2, "All" -> 3 (*, All -> 3 *)|>, mode, -1]
