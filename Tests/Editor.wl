(* ::Package:: *)

toState = IGraphM`GraphEditor`PackagePrivate`GraphToEditorState;


VerificationTest[
  {
    SetOptions[IGGraphEditor,ImageSize->200]
  , toState[]["config", "ImageSize"]
  , SetOptions[IGGraphEditor,ImageSize->Automatic]
  } [[2]]
, 200
, TestID -> "SetOptions"
]


VerificationTest[
  state = toState @ Graph[{1<->2}]
; state["config", "DirectedEdges"]
, False
, TestID -> "Two way !directed"
]


VerificationTest[
  state = toState @ Graph[{1, 2}, {}]
; state["config", "DirectedEdges"]
, False
, TestID -> "Empty !directed"
]


VerificationTest[
  state = toState @ Graph[{1->2, 2->1}]
; state["config", "DirectedEdges"]
, True
, TestID -> "directed graph directed"
]


VerificationTest[
  config = toState[PathGraph[{1,2}, VertexLabels->"Index"],DirectedEdges->True, VertexLabels->"Name"]["config"]
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
