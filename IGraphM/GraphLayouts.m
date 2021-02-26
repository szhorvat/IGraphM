(* Mathematica Package *)
(* Created by Mathematica plugin for IntelliJ IDEA *)

(* :Author: szhorvat *)
(* :Date: 2018-10-24 *)
(* :Copyright: (c) 2018-2020 Szabolcs Horvát *)

Package["IGraphM`"]


(*************************)
(***** Graph layouts *****)
(*************************)


(* Common messages for layout functions. *)

IGraphM::lytcrd = "The graph doesn't already have existing vertex coordinates. The \"Continue\" -> True layout option will be ignored.";
IGraphM::lytdim = "The existing vertex coordinates do not have the appropriate dimension for this layout algorithm. The \"Continue\" -> True layout option will be ignored.";
IGraphM::lytcnt = "`` is not a valid value for the \"Continue\" layout option.";
IGraphM::lytaln = "`` is not a valid value for the \"Align\" layout option."


(***** Helper functions for graph layouts *****)

getVertexCoords[graph_] :=
    With[{coords = GraphEmbedding[graph]},
      If[MatrixQ[coords],
        coords,
        Message[IGraphM::lytcrd]; {{}}
      ]
    ]

continueLayout[graph_, False, ___] := Sequence[{{}}, False]
continueLayout[graph_, True, scale_ : 1, dim_ : 2] :=
    Sequence@@Module[{coords},
      coords = getVertexCoords[graph];
      If[coords =!= {{}} && Not@MatchQ[Dimensions[coords], {_, dim}],
        Message[IGraphM::lytdim];
        coords = {{}}
      ];
      {coords / scale, coords =!= {{}}}
    ]
continueLayout[graph_, cont_, ___] := ( Message[IGraphM::lytcnt, cont]; continueLayout[graph, False] )
continueLayout3D[graph_, cont_, scale_ : 1] := continueLayout[graph, cont, scale, 3]

(* These should simply use Graph and Graph3D but due to a Mathematica bug,
   that won't work on some graphs, such as those returns by KaryTree.
   Thus we use SetProperty with GraphLayout instead.
*)
setVertexCoords[g_, coords_] := SetProperty[g, {VertexCoordinates -> Thread[ VertexList[g] -> coords ], PlotRange -> All}]
setVertexCoords3D[g_, coords_] := SetProperty[g, {GraphLayout -> {"Dimension" -> 3}, VertexCoordinates -> Thread[ VertexList[g] -> coords ], PlotRange -> All}]


igAlign[{}] := {}
igAlign[pts_] := PrincipalComponents[pts]

align[True] = igAlign;
align[False] = Identity;
align[val_] := (Message[IGraphM::lytaln, val]; Identity)


(***** Main layout functions *****)

PackageExport["IGLayoutRandom"]
IGLayoutRandom::usage = "IGLayoutRandom[graph] lays out vertices randomly in the unit square.";

SyntaxInformation[IGLayoutRandom] = {"ArgumentsPattern" -> {_, OptionsPattern[]}, "OptionNames" -> optNames[Graph]};
IGLayoutRandom[graph_?igGraphQ, opt : OptionsPattern[Graph]] :=
    catch@Block[{ig = igMakeFast[graph]},
      applyGraphOpt[opt]@setVertexCoords[graph, check@ig@"layoutRandom"[]]
    ]


PackageExport["IGLayoutCircle"]
IGLayoutCircle::usage = "IGLayoutCircle[graph] lays out vertices on a circle.";

Options[IGLayoutCircle] = { "Rotation" -> 0, Reverse -> False };
SyntaxInformation[IGLayoutCircle] = {"ArgumentsPattern" -> {_, OptionsPattern[]}, "OptionNames" -> optNames[IGLayoutCircle, Graph]};
IGLayoutCircle[graph_?igGraphQ, opt : OptionsPattern[{IGLayoutCircle, Graph}]] :=
    catch@Block[{ig = igMakeFast[graph]},
      applyGraphOpt[opt]@setVertexCoords[
        graph,
        Composition[
          RotationTransform[OptionValue["Rotation"]],
          If[TrueQ@OptionValue[Reverse], ScalingTransform[{1, -1}], Identity]
        ] @ check@ig@"layoutCircle"[]
      ]
    ]


PackageExport["IGLayoutSphere"]
IGLayoutSphere::usage = "IGLayoutSphere[graph] lays out vertices approximately uniformly distributed on a sphere.";

SyntaxInformation[IGLayoutSphere] = {"ArgumentsPattern" -> {_, OptionsPattern[]}, "OptionNames" -> optNames[Graph3D]};
IGLayoutSphere[graph_?igGraphQ, opt : OptionsPattern[Graph3D]] :=
    catch@Block[{ig = igMakeFast[graph]},
      applyGraphOpt3D[opt]@setVertexCoords3D[graph, check@ig@"layoutSphere"[]]
    ]


PackageExport["IGLayoutGraphOpt"]
IGLayoutGraphOpt::usage = "IGLayoutGraphOpt[graph, options] lays out the graph using the GraphOpt algorithm.";

Options[IGLayoutGraphOpt] = {
  "MaxIterations" -> 500, "NodeCharge" -> 0.001, "NodeMass" -> 30, "SpringLength" -> 0,
  "SpringConstant" -> 1, "MaxStepMovement" -> 5,
  "Continue" -> False, "Align" -> True
};

SyntaxInformation[IGLayoutGraphOpt] = {"ArgumentsPattern" -> {_, OptionsPattern[]}, "OptionNames" -> optNames[IGLayoutGraphOpt, Graph]};

IGLayoutGraphOpt[graph_?igGraphQ, opt : OptionsPattern[{IGLayoutGraphOpt,Graph}]] :=
    catch@Block[{ig = igMakeFastWeighted[graph], scale = 0.01},
      applyGraphOpt[opt]@setVertexCoords[graph,
          scale align[OptionValue["Align"]]@check@ig@"layoutGraphOpt"[continueLayout[graph, OptionValue["Continue"], scale],
            OptionValue["MaxIterations"], OptionValue["NodeCharge"], OptionValue["NodeMass"],
            OptionValue["SpringLength"], OptionValue["SpringConstant"], OptionValue["MaxStepMovement"]
          ]
      ]
    ]


PackageExport["IGLayoutKamadaKawai"]
IGLayoutKamadaKawai::usage = "IGLayoutKamadaKawai[graph, options] lays out the graph using the Kamada–Kawai algorithm (similar to \"SpringEmbedding\").";

Options[IGLayoutKamadaKawai] = {
  "MaxIterations" -> Automatic, "Epsilon" -> 0, "KamadaKawaiConstant" -> Automatic,
  "Continue" -> False, "Align" -> True
};

SyntaxInformation[IGLayoutKamadaKawai] = {"ArgumentsPattern" -> {_, OptionsPattern[]}, "OptionNames" -> optNames[IGLayoutKamadaKawai, Graph]};

IGLayoutKamadaKawai[graph_?igGraphQ, opt : OptionsPattern[{IGLayoutKamadaKawai,Graph}]] :=
    catch@Block[{ig = igMakeFastWeighted[graph], maxiter, kkconst, scale = 0.5},
      maxiter = Replace[OptionValue["MaxIterations"], Automatic :> 10 VertexCount[graph]];
      kkconst = Replace[OptionValue["KamadaKawaiConstant"], Automatic :> Max[1, VertexCount[graph]]];
      applyGraphOpt[opt]@setVertexCoords[graph,
        scale align[OptionValue["Align"]]@check@ig@"layoutKamadaKawai"[continueLayout[graph, OptionValue["Continue"], scale],
          maxiter, OptionValue["Epsilon"], kkconst]
      ]
    ]


PackageExport["IGLayoutKamadaKawai3D"]
IGLayoutKamadaKawai3D::usage = "IGLayoutKamadaKawai3D[graph, options] lays out the graph in 3D using the Kamada–Kawai algorithm (similar to \"SpringEmbedding\").";

Options[IGLayoutKamadaKawai3D] = {
  "MaxIterations" -> Automatic, "Epsilon" -> 0, "KamadaKawaiConstant" -> Automatic,
  "Continue" -> False, "Align" -> True
};

SyntaxInformation[IGLayoutKamadaKawai3D] = {"ArgumentsPattern" -> {_, OptionsPattern[]}, "OptionNames" -> optNames[IGLayoutKamadaKawai3D, Graph3D]};

IGLayoutKamadaKawai3D[graph_?igGraphQ, opt : OptionsPattern[{IGLayoutKamadaKawai3D,Graph3D}]] :=
    catch@Block[{ig = igMakeFastWeighted[graph], maxiter, kkconst, scale = 0.5},
      maxiter = Replace[OptionValue["MaxIterations"], Automatic :> 10 VertexCount[graph]];
      kkconst = Replace[OptionValue["KamadaKawaiConstant"], Automatic :> Max[1, VertexCount[graph]]];
      applyGraphOpt3D[opt]@setVertexCoords3D[graph,
        scale align[OptionValue["Align"]]@check@ig@"layoutKamadaKawai3D"[continueLayout3D[graph, OptionValue["Continue"], scale],
          maxiter, OptionValue["Epsilon"], kkconst]
      ]
    ]


toBounds[v : {_?NumericQ, _?NumericQ}] := {v, v}
toBounds[m : {{_?NumericQ, _?NumericQ}, {_?NumericQ, _?NumericQ}}] := m
toBounds[_] := throw[$Failed]

makeCons[g_, asc_?AssociationQ, scale] :=
    With[{size = (Max[Abs[asc]] + 2.0 Sqrt@VertexCount[g])},
      Join @@ Transpose[
        Lookup[
          (toBounds /@ N[asc]) / scale, VertexList[g],
          {{-1, -1}, {1, 1}} size
        ],
        {3, 1, 2}
      ]
    ]

PackageExport["IGLayoutFruchtermanReingold"]
IGLayoutFruchtermanReingold::usage = "IGLayoutFruchtermanReingold[graph, options] lays out the graph using the Fruchterman–Reingold algorithm (similar to \"SpringElectricalEmbedding\").";

igFruchtermanReingoldMethods = <| Automatic -> 2, False -> 1, True -> 0 |>;

Options[IGLayoutFruchtermanReingold] = {
  "MaxIterations" -> 500, "MaxMovement" -> 5, "UseGrid" -> Automatic,
  "Continue" -> False, "Align" -> Automatic,
  "Constraints" -> None
};

SyntaxInformation[IGLayoutFruchtermanReingold] = {"ArgumentsPattern" -> {_, OptionsPattern[]}, "OptionNames" -> optNames[IGLayoutFruchtermanReingold, Graph]};

IGLayoutFruchtermanReingold[graph_?igGraphQ, opt : OptionsPattern[{IGLayoutFruchtermanReingold,Graph}]] :=
    catch@Module[{ig = igMakeFastWeighted[graph], scale = 0.25, constraints, cons, al},
      constraints = Replace[OptionValue["Constraints"], None -> <||>];
      cons = Length[constraints] > 0;
      al = Replace[OptionValue["Align"], Automatic -> Not[cons]];
      constraints = If[cons, makeCons[graph, constraints, scale], {{},{},{},{}}];
      applyGraphOpt[opt]@setVertexCoords[graph,
        scale align[al]@check@ig@"layoutFruchtermanReingold"[continueLayout[graph, OptionValue["Continue"], scale],
          OptionValue["MaxIterations"], OptionValue["MaxMovement"], Lookup[igFruchtermanReingoldMethods, OptionValue["UseGrid"], -1],
          cons, Sequence @@ constraints
        ]
      ]
    ]


PackageExport["IGLayoutFruchtermanReingold3D"]
IGLayoutFruchtermanReingold3D::usage = "IGLayoutFruchtermanReingold3D[graph, options] lays out the graph using the Fruchterman–Reingold algorithm (similar to \"SpringElectricalEmbedding\").";

Options[IGLayoutFruchtermanReingold3D] = {
  "MaxIterations" -> 500, "MaxMovement" -> 5,
  "Continue" -> False, "Align" -> True
};

SyntaxInformation[IGLayoutFruchtermanReingold3D] = {"ArgumentsPattern" -> {_, OptionsPattern[]}, "OptionNames" -> optNames[IGLayoutFruchtermanReingold3D, Graph3D]};

IGLayoutFruchtermanReingold3D[graph_?igGraphQ, opt : OptionsPattern[{IGLayoutFruchtermanReingold3D,Graph3D}]] :=
    catch@Block[{ig = igMakeFastWeighted[graph], scale = 0.25},
      applyGraphOpt3D[opt]@setVertexCoords3D[graph,
        scale align[OptionValue["Align"]]@check@ig@"layoutFruchtermanReingold3D"[continueLayout3D[graph, OptionValue["Continue"], scale],
          OptionValue["MaxIterations"], OptionValue["MaxMovement"]
        ]
      ]
    ]


PackageExport["IGLayoutGEM"]
IGLayoutGEM::usage = "IGLayoutGEM[graph, options] lays out the graph using the GEM algorithm.";

Options[IGLayoutGEM] = {
  "MaxIterations" -> Automatic, "Continue" -> False, "Align" -> True,
  "MaxTemperature" -> Automatic, "MinTemperature" -> 1/10, "InitTemperature" -> Automatic
};

SyntaxInformation[IGLayoutGEM] = {"ArgumentsPattern" -> {_, OptionsPattern[]}, "OptionNames" -> optNames[IGLayoutGEM, Graph]};

IGLayoutGEM[graph_?igGraphQ, opt : OptionsPattern[{IGLayoutGEM,Graph}]] :=
    catch@Block[{ig = igMakeFastWeighted[graph], maxiter, maxtemp, inittemp, scale = 3*^-3},
      maxiter = Replace[OptionValue["MaxIterations"], Automatic :> 40 VertexCount[graph]^2];
      maxtemp = Replace[OptionValue["MaxTemperature"], Automatic :> Max[1, VertexCount[graph]]];
      inittemp = Replace[OptionValue["InitTemperature"], Automatic :> Max[1, Sqrt@VertexCount[graph]]];
      applyGraphOpt[opt]@setVertexCoords[graph,
        scale align[OptionValue["Align"]]@check@ig@"layoutGEM"[continueLayout[graph, OptionValue["Continue"], scale],
          maxiter, maxtemp, OptionValue["MinTemperature"], inittemp
        ]
      ]
    ]


PackageExport["IGLayoutDavidsonHarel"]
IGLayoutDavidsonHarel::usage = "IGLayoutDavidsonHarel[graph, options] lays out the graph using the Davidson–Harel algorithm, based on simulated annealing.";

Options[IGLayoutDavidsonHarel] = {
  "MaxIterations" -> 10, "Continue" -> False, "Align" -> True,
  "FineTuningIterations" -> Automatic, "CoolingFactor" -> 0.75,
  "NodeDistanceWeight" -> 1.0, "BorderDistanceWeight" -> 0.0, "EdgeLengthWeight" -> Automatic,
  "EdgeCrossingWeight" -> Automatic, "EdgeDistanceWeight" -> Automatic
};

SyntaxInformation[IGLayoutDavidsonHarel] = {"ArgumentsPattern" -> {_, OptionsPattern[]}, "OptionNames" -> optNames[IGLayoutDavidsonHarel, Graph]};

IGLayoutDavidsonHarel[graph_?igGraphQ, opt : OptionsPattern[{IGLayoutDavidsonHarel,Graph}]] :=
    catch@Block[{ig = igMakeFastWeighted[graph], tuneiter, edgelenw, edgecrossw, edgedistw, dens, scale = 0.1},
      dens = If[VertexCount[graph] <= 1, 0, GraphDensity[graph]];
      tuneiter = Replace[OptionValue["FineTuningIterations"], Automatic :> Max[10, Round@Log[2, VertexCount[graph]]]];
      edgelenw = Replace[OptionValue["EdgeLengthWeight"], Automatic :> dens/10];
      edgecrossw = Replace[OptionValue["EdgeCrossingWeight"], Automatic :> 1 - dens];
      edgedistw = Replace[OptionValue["EdgeDistanceWeight"], Automatic :> 1 - dens / 5];
      applyGraphOpt[opt]@setVertexCoords[graph,
        scale align[OptionValue["Align"]]@check@ig@"layoutDavidsonHarel"[continueLayout[graph, OptionValue["Continue"], scale],
          OptionValue["MaxIterations"],
          tuneiter, OptionValue["CoolingFactor"], OptionValue["NodeDistanceWeight"],
          OptionValue["BorderDistanceWeight"], edgelenw, edgecrossw, edgedistw
        ]
      ]
    ]

(*

PackageScope["IGLayoutMDS"]
IGLayoutMDS::usage = "IGLayoutMDS[graph]";

IGLayoutMDS[graph_?igGraphQ, dim : (2|3) : 2, Optional[distMatrix_?SquareMatrixQ, Automatic]] :=
    catch@Block[{ig = igMake[graph]},
      setVertexCoords[graph,
        check@ig@"layoutMDS"[
          Replace[distMatrix, Automatic -> {{}}],
          dim
        ]
      ]
    ]
*)

PackageExport["IGLayoutReingoldTilford"]
IGLayoutReingoldTilford::usage = "IGLayoutReingoldTilford[graph, options] lays out a tree using the Reingold–Tilford algorithm.";

(* TODO: Do this in C eventually as a workaround for connectedGraphComponents/Subgraph unreliability *)
chooseRoots[graph_?UndirectedGraphQ] := First@GraphCenter[#]& /@ WeaklyConnectedGraphComponents[graph]
chooseRoots[graph_? (IGForestQ[#, "Out"]&)] := First@TopologicalSort[#]& /@ WeaklyConnectedGraphComponents[graph]
chooseRoots[graph_? (IGForestQ[#, "In"]&)] := Last@TopologicalSort[#]& /@ WeaklyConnectedGraphComponents[graph]
chooseRoots[graph_] := First@GraphCenter[#]& /@ WeaklyConnectedGraphComponents@UndirectedGraph[graph]

Options[IGLayoutReingoldTilford] = {
  "RootVertices" -> Automatic, "Rotation" -> 0,
  "LayerHeight" -> 1, "LeafDistance" -> 1
};

SyntaxInformation[IGLayoutReingoldTilford] = {"ArgumentsPattern" -> {_, OptionsPattern[]}, "OptionNames" -> optNames[IGLayoutReingoldTilford, Graph]};

IGLayoutReingoldTilford[graph_?igGraphQ, opt : OptionsPattern[{IGLayoutReingoldTilford,Graph}]] :=
    catch@Block[{ig = igMakeFast[graph], roots},
      roots = vss[graph]@Replace[OptionValue["RootVertices"], Automatic :> chooseRoots[graph]];
      applyGraphOpt[opt]@setVertexCoords[graph,
        Composition[
          RotationTransform[OptionValue["Rotation"]],
          ScalingTransform[{OptionValue["LeafDistance"], -OptionValue["LayerHeight"]}]
        ] @ check@ig@"layoutReingoldTilford"[roots, False]
      ]
    ]


PackageExport["IGLayoutReingoldTilfordCircular"]
IGLayoutReingoldTilfordCircular::usage = "IGLayoutReingoldTilfordCircular[graph, options] lays out a tree radially using the Reingold–Tilford algorithm.";

Options[IGLayoutReingoldTilfordCircular] = { "RootVertices" -> Automatic, "Rotation" -> 0 };

SyntaxInformation[IGLayoutReingoldTilfordCircular] = {"ArgumentsPattern" -> {_, OptionsPattern[]}, "OptionNames" -> optNames[IGLayoutReingoldTilfordCircular, Graph]};

IGLayoutReingoldTilfordCircular[graph_?igGraphQ, opt : OptionsPattern[{IGLayoutReingoldTilfordCircular,Graph}]] :=
    catch@Block[{ig = igMakeFast[graph], roots},
      roots = vss[graph]@Replace[OptionValue["RootVertices"], Automatic :> chooseRoots[graph]];
      applyGraphOpt[opt]@setVertexCoords[graph,
        RotationTransform[OptionValue["Rotation"]] @ check@ig@"layoutReingoldTilfordCircular"[roots, False]
      ]
    ]


PackageExport["IGLayoutDrL"]
IGLayoutDrL::usage = "IGLayoutDrL[graph, options] lays out the graph using the DrL layout generator.";
PackageExport["IGLayoutDrL3D"]
IGLayoutDrL3D::usage = "IGLayoutDrL3D[graph, options] lays out the graph in 3D using the DrL layout generator.";

Options[IGLayoutDrL] = { "Settings" -> "Default", "Continue" -> False, "Align" -> True };
Options[IGLayoutDrL3D] = { "Settings" -> "Default", "Continue" -> False, "Align" -> True };

SyntaxInformation[IGLayoutDrL] = {"ArgumentsPattern" -> {_, OptionsPattern[]}, "OptionNames" -> optNames[IGLayoutDrL, Graph]};
SyntaxInformation[IGLayoutDrL3D] = {"ArgumentsPattern" -> {_, OptionsPattern[]}, "OptionNames" -> optNames[IGLayoutDrL3D, Graph3D]};

igLayoutDrLSettings = {"Default", "Coarsen", "Coarsest", "Refine", "Final"};
igLayoutDrLSettingsAsc = AssociationThread[igLayoutDrLSettings, Range@Length[igLayoutDrLSettings]];


amendUsage[#, "Possible values for the \"Settings\" option are <*igLayoutDrLSettings*>."]& /@ {IGLayoutDrL, IGLayoutDrL3D};

IGLayoutDrL::conn = "IGLayoutDrL may fail on disconnected graphs. Use on connected graphs only.";

IGLayoutDrL[graph_?igGraphQ, opt : OptionsPattern[{IGLayoutDrL,Graph}]] :=
    catch@Block[{ig = igMakeFastWeighted[graph], scale = 0.25},
      If[Not@WeaklyConnectedGraphQ[graph], Message[IGLayoutDrL::conn]];
      applyGraphOpt[opt]@setVertexCoords[graph,
        scale align[OptionValue["Align"]]@check@ig@"layoutDrL"[continueLayout[graph, OptionValue["Continue"], scale],
          Lookup[igLayoutDrLSettingsAsc, OptionValue["Settings"], -1]
        ]
      ]
    ]

IGLayoutDrL3D::conn = IGLayoutDrL::conn;

IGLayoutDrL3D[graph_?igGraphQ, opt : OptionsPattern[{IGLayoutDrL3D,Graph3D}]] :=
    catch@Block[{ig = igMakeFastWeighted[graph], scale = 0.25},
      If[Not@WeaklyConnectedGraphQ[graph], Message[IGLayoutDrL3D::conn]];
      applyGraphOpt[opt]@setVertexCoords3D[graph,
        scale align[OptionValue["Align"]]@check@ig@"layoutDrL3D"[continueLayout3D[graph, OptionValue["Continue"], scale],
          Lookup[igLayoutDrLSettingsAsc, OptionValue["Settings"], -1]
        ]
      ]
    ]


PackageExport["IGLayoutTutte"]
IGLayoutTutte::usage = "IGLayoutTutte[graph, options] lays out a 3-vertex-connected planar graph using the Tutte embedding.";

Options[IGLayoutTutte] = { "OuterFace" -> Automatic };
SyntaxInformation[IGLayoutTutte] = {"ArgumentsPattern" -> {_, OptionsPattern[]}, "OptionNames" -> optNames[IGLayoutTutte, Graph]};
IGLayoutTutte::n3c = "The graph is not 3-vertex-connected and a Tutte embedding cannot be computed. Vertex coordinates will not be set.";
IGLayoutTutte::npl = "The graph is not planar and a Tutte embedding cannot be computed. Vertex coordinates will not be set.";
IGLayoutTutte::bdface = "`1` is not a face of the graph.";
IGLayoutTutte[graph_?igGraphQ, opt : OptionsPattern[{IGLayoutTutte, Graph}]] :=
    catch@Module[{mat, faces, facePos, face, nonFace, pts, circPts, restPts, compl},
      If[Not@KVertexConnectedGraphQ[UndirectedGraph[graph], 3],
        Message[IGLayoutTutte::n3c];
        Return[graph]
      ];
      mat =
          zeroDiagonal@N@If[IGEdgeWeightedQ[graph] && VectorQ[IGEdgeProp[EdgeWeight][graph], NumericQ],
            WeightedAdjacencyMatrix@Quiet[IGWeightedUndirectedGraph[graph], IGWeightedUndirectedGraph::mg],
            AdjacencyMatrix@SimpleGraph[graph]
          ];
      mat = mat/Total[mat];
      faces = IGFaces@IndexGraph[graph];
      If[Not@ListQ[faces],
        Message[IGLayoutTutte::npl];
        Return[graph]
      ];
      face = OptionValue["OuterFace"];
      If[face === Automatic,
        face = faces[[ First@Ordering[Length /@ faces, -1] ]]
        ,
        If[Not@ListQ[face],
          Message[IGLayoutTutte::bdface, OptionValue["OuterFace"]];
          throw[$Failed]
        ];
        face = 1 + vss[graph][face];
        facePos = Position[Sort /@ faces, Sort[face], {1}];
        If[facePos === {},
          Message[IGLayoutTutte::bdface, OptionValue["OuterFace"]];
          throw[$Failed]
        ];
        face = Extract[faces, First[facePos]];
      ];
      nonFace = Complement[Range@VertexCount[graph], face];
      circPts = N@CirclePoints@Length[face];
      compl = mat[[nonFace, face]].circPts;
      restPts =
          LinearSolve[
            mat[[nonFace, nonFace]] - IdentityMatrix[Length[nonFace], SparseArray],
            -compl
          ];
      pts = ConstantArray[0., {VertexCount[graph], 2}];
      pts[[face]] = circPts;
      pts[[nonFace]] = restPts;
      applyGraphOpt[opt]@setVertexCoords[graph, pts]
    ]


PackageExport["IGLayoutBipartite"]
IGLayoutBipartite::usage = "IGLayoutBipartite[graph, options] lays out a bipartite graph, minimizing the number of edge crossings. Partitions can be specified manually using the \"BipartitePartitions\" option.";

blockLayout[n_, gap_, height_, mult_, offset_] :=
    Module[{nrow = Ceiling[1.01 height/gap]},
      Table[gap QuotientRemainder[i, nrow] mult + offset, {i, 0, n - 1}]
    ]

(* Deletes elements in list1 that are also present in list2.
   Assumes that both list1 and list2 are free of duplicates.
   Thanks to Mr. Wizard: https://mathematica.stackexchange.com/a/1295 *)
unsortedComplement[list1_, {}] := list1
unsortedComplement[list1_, list2_] := Drop[DeleteDuplicates@Join[list2, list1], Length[list2]]

Options[IGLayoutBipartite] = {
  "VertexGap" -> 0.1,
  "PartitionGap" -> 1,
  "Orientation" -> Vertical,
  MaxIterations -> 100,
  "BipartitePartitions" -> Automatic
};
SyntaxInformation[IGLayoutBipartite] = {"ArgumentsPattern" -> {_, OptionsPattern[]}, "OptionNames" -> optNames[IGLayoutBipartite, Graph]};
IGLayoutBipartite::invopt = "The option value `` is not valid.";
IGLayoutBipartite::notbp = "Cannot determine partitions because the graph is not bipartite. Specify partitions explicitly using the \"BipartitePartitions\" option.";
IGLayoutBipartite::bdprt = "`` is not a valid partitioning vertices.";
IGLayoutBipartite[graph_ /; VertexCount[graph] == 0, opt : OptionsPattern[{IGLayoutBipartite, Graph}]] :=
    applyGraphOpt[opt][graph]
IGLayoutBipartite[graph_?igGraphQ, opt : OptionsPattern[{IGLayoutBipartite, Graph}]] :=
    catch@Module[{sg, ig, parts, isolated, connected, vertical, types, coord, coordAsc, min, max},
      parts = OptionValue["BipartitePartitions"];
      If[parts === Automatic,
        isolated = Pick[VertexList[graph], VertexDegree[graph], 0];
        connected = unsortedComplement[VertexList[graph], isolated];
        sg = igSubgraph[graph, connected];
        ig = igMakeFast[sg];
        If[Not@ig@"bipartiteQ"[],
          Message[IGLayoutBipartite::notbp];
          Return[$Failed]
        ];
        types = ig@"bipartitePartitions"[]
        ,
        If[Not[ MatchQ[parts, {_List, _List}] && SubsetQ[VertexList[graph], Join@@parts] && Intersection@@parts === {} ],
          Message[IGLayoutBipartite::bdprt, "BipartitePartitions" -> parts];
          Return[$Failed]
        ];
        connected = Join @@ parts;
        isolated = Complement[VertexList[graph], connected];
        sg = igSubgraph[graph, connected];
        ig = igMakeFast[sg];
        types = communitiesToMembership[VertexList[sg], parts];
      ];

      Switch[OptionValue["Orientation"],
        Horizontal|"Horizontal", vertical = False,
        Vertical|"Vertical", vertical = True,
        _, Message[IGLayoutBipartite::invopt, "Orientation" -> OptionValue["Orientation"]]; Return[$Failed]
      ];

      If[Not@positiveNumericQ@OptionValue["VertexGap"],
        Message[IGLayoutBipartite::invopt, "VertexGap" -> OptionValue["VertexGap"]]; Return[$Failed]
      ];

      If[Not@positiveNumericQ@OptionValue["PartitionGap"],
        Message[IGLayoutBipartite::invopt, "PartitionGap" -> OptionValue["PartitionGap"]]; Return[$Failed]
      ];

      (* use below instead of Min/Max to avoid obtaining Infinity for empty lists *)
      min[{}] = 0; min[x_] := Min[x];
      max[{}] = 0; max[x_] := Max[x];

      coord = check@ig@"layoutBipartite"[types, OptionValue["VertexGap"], OptionValue["PartitionGap"], OptionValue[MaxIterations]];
      coordAsc = Join[
        AssociationThread[connected -> coord],
        AssociationThread[
          isolated -> blockLayout[
            Length[isolated],
            OptionValue["VertexGap"], OptionValue["PartitionGap"],
            If[vertical, {-1, -1}, {1, -1}],
            If[vertical,
              {min[coord[[All, 1]]] - 1.5 OptionValue["VertexGap"], OptionValue["PartitionGap"]}
              ,
              {max[coord[[All, 1]]] + 1.5 OptionValue["VertexGap"], OptionValue["PartitionGap"]}
            ]
          ]
        ]
      ];

      applyGraphOpt[opt]@setVertexCoords[graph,
        If[vertical, RotationTransform[Pi/2], Identity]@Lookup[coordAsc, VertexList[graph]]
      ]
    ]


PackageExport["IGLayoutPlanar"]
IGLayoutPlanar::usage = "IGLayoutPlanar[graph, options] lays out a planar graph using Schnyder's algorithm.";

SyntaxInformation[IGLayoutPlanar] = {"ArgumentsPattern" -> {_, OptionsPattern[]}, "OptionNames" -> optNames[IGLayoutPlanar, Graph]};
IGLayoutPlanar[graph_?igGraphQ, opt : OptionsPattern[Graph]] :=
    catch@Block[{lg = lgMake@UndirectedGraph@SimpleGraph[graph]},
      applyGraphOpt[opt]@setVertexCoords[graph, check@lg@"layoutPlanar"[]]
    ]