(* ::Package:: *)

toState = IGraphM`GraphEditor`PackagePrivate`GraphToEditorState;


VerificationTest[
  {
    SetOptions[IGGraphEditor,ImageSize->200]
  , toState[][ "ImageSize"]
  , SetOptions[IGGraphEditor,ImageSize->Automatic]
  } [[2]]
, 200
, TestID -> "SetOptions"
]


VerificationTest[
  state = toState @ Graph[{1<->2}]
; state[ "DirectedEdges"]
, False
, TestID -> "Two way !directed"
]


VerificationTest[
  state = toState @ Graph[{1, 2}, {}]
; state[ "DirectedEdges"]
, False
, TestID -> "Empty !directed"
]


VerificationTest[
  state = toState @ Graph[{1->2, 2->1}]
; state[ "DirectedEdges"]
, True
, TestID -> "directed graph directed"
]


VerificationTest[
  config = toState[PathGraph[{1,2}, VertexLabels->"Index"],DirectedEdges->True, VertexLabels->"Name"]
; config["DirectedEdges"]
, False
, TestID -> "Input graph DirectedEdges > options DirectedEdges"
]


VerificationTest[
  config["VertexLabels"]
, "Name"
, TestID -> "Input graph VertexLabels < options VertexLabels"
]


VerificationTest[
  #2-#&@@@config["range"] // Apply[#2-#&] // # < 10.^-15&
, True
, TestID -> "Path graph narrow range"
]


VerificationTest[
  MatchQ[
  IGraphM`GraphEditor`PackagePrivate`extractEdgePrimitives@toState@Graph@{1->2,2->1}
, {KeyValuePattern[{"id" -> _String, "primitive" -> Verbatim[Arrow][_BezierCurve,_]}]..}  
]
, True
, TestID -> "extractEdgePrimitives"
]


v1boxes = "1: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";


VerificationTest[
  ToExpression@Uncompress@v1boxes
, Graph[{1, 2}, {}, {VertexCoordinates -> {{0.00001, -0.000030000000000000004}, {0.00001, -0.00001}}}]
, TestID -> "Handling legacy typeset"
]
