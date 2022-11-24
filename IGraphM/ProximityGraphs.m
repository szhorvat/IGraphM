(* Mathematica Package *)
(* Created by Mathematica plugin for IntelliJ IDEA *)

(* :Author: szhorvat *)
(* :Date: 2018-10-23 *)
(* :Copyright: (c) 2018-2022 Szabolcs Horvát *)

Package["IGraphM`"]


PackageImport["TriangleLink`"]
PackageImport["TetGenLink`"]

(****************************)
(***** Proximity graphs *****)
(****************************)


PackageExport["IGDelaunayGraph"]
IGDelaunayGraph::usage = "IGDelaunayGraph[points] gives the Delaunay graph of the given points.";

IGDelaunayGraph::dim  = "Delaunay graph computation is currently only supported in 2 and 3 dimensions.";
IGDelaunayGraph::dupl = "Remove any duplicate points before the Delaunay graph computation.";
IGDelaunayGraph::fail = "Could not compute Delaunay triangulation."; (* This is triggered if *)

(* Replacement for TriangleDelaunay[]; Removes duplicate points from result. *)
triangleDelaunay[points_] :=
    Module[{in, out, pts, elements},
      in = TriangleCreate[];
      TriangleSetPoints[in, points];
      out = TriangleTriangulate[in, "-jQ "];
      If[Not@ManagedLibraryExpressionQ[out], Return[$Failed]];
      pts = TriangleGetPoints[out];
      elements = TriangleGetElements[out];
      If[elements === {}, Return[$Failed]];
      {pts, elements}
    ]

delaunayEdges1D[points_] :=
    With[{pts = N[points]},
      If[DuplicateFreeQ[pts], (* TODO: Eventually replace with a DuplicateFreeQ equivalent that uses tolerances *)
        Partition[Ordering[pts], 2, 1],
        Message[IGDelaunayGraph::dupl]; throw[$Failed]
      ]
    ]

delaunayEdges2D[points_] :=
    Switch[Length[points],
      0 | 1, {},

      2,
      If[Unequal@@N[points],
        {{1, 2}},
        Message[IGDelaunayGraph::dupl]; throw[$Failed]
      ],

      _,
      (* Workaround for: TriangleLink crashes when given all-identical identical points.
         We only check the first two points in the list to avoid unpacking. *)
      If[N[points[[1]]] == N[points[[2]]],
        Message[IGDelaunayGraph::dupl]; throw[$Failed]
      ];
      Module[{res, pts, triangles, v1, v2},
        res = triangleDelaunay[points];
        If[res === $Failed,
          (* triangleDelaunay[] failed: check if points are collinear and if yes, fall back to 1D Delaunay *)
          pts = PrincipalComponents@N[points];
          {v1, v2} = Variance[pts];
          If[v2/v1 < 10^Internal`$EqualTolerance $MachineEpsilon,
            delaunayEdges1D[ pts[[All,1]] ],
            Message[IGDelaunayGraph::fail]; throw[$Failed]
          ]
          ,
          (* triangleDelaunay[] succeeded: proceed as usual *)
          {pts, triangles} = res;
          If[Length[pts] == Length[points],
            DeleteDuplicates@igraphGlobal@"edgeListSortPairs"[
              Join @@ Transpose /@ Subsets[Transpose[triangles], {2}]
            ],
            Message[IGDelaunayGraph::dupl]; throw[$Failed]
          ]
        ]
      ]
    ]

delaunayEdges3D[points_] :=
    Switch[Length[points],
      0 | 1, {},

      2,
      If[Unequal@@N[points],
        {{1, 2}},
        Message[IGDelaunayGraph::dupl]; throw[$Failed]
      ],

      (* 3: TetGenDelaunay fails gracefully for 3 points, thus we do not need an extra case for that. *)

      _,
      Module[{res, pts, tetrahedra, v1, v2, v3},
        res = Quiet[TetGenDelaunay[points], TetGenDelaunay::tetfc];
        If[res === $Failed,
          (* TetGenDelaunay failed: check if points are collinear and if yes, fall back to 2D Delaunay *)
          pts = PrincipalComponents@N[points];
          {v1, v2, v3} = Variance[pts];
          If[v3/v1 < 10^Internal`$EqualTolerance $MachineEpsilon,
            delaunayEdges2D[ pts[[All,{1,2}]] ],
            Message[IGDelaunayGraph::fail]; throw[$Failed]
          ]
          ,
          {pts, tetrahedra} = res;
          If[Length[pts] == Length[points],
            DeleteDuplicates@igraphGlobal@"edgeListSortPairs"[
              Join @@ Transpose /@ Subsets[Transpose[tetrahedra], {2}]
            ],
            Message[IGDelaunayGraph::dupl]; throw[$Failed]
          ]
        ]
      ]
    ]


SyntaxInformation[IGDelaunayGraph] = {"ArgumentsPattern" -> {_, OptionsPattern[]}, "OptionNames" -> optNames[Graph, Graph3D]};
IGDelaunayGraph[{}, opt : OptionsPattern[Graph]] := IGEmptyGraph[0, opt]
IGDelaunayGraph[points_?(MatrixQ[#, NumericQ]&), opt : OptionsPattern[{Graph, Graph3D}]] :=
    catch@Switch[Last@Dimensions[points], (* points is known to be a matrix, so we know that Length@Dimensions[points] == 2 *)
      1,   Graph[Range@Length[points], delaunayEdges1D[points], DirectedEdges -> False, opt, VertexCoordinates -> ArrayPad[points, {{0, 0}, {0, 1}}]],
      2,   Graph[Range@Length[points], delaunayEdges2D[points], DirectedEdges -> False, opt, VertexCoordinates -> points],
      3, Graph3D[Range@Length[points], delaunayEdges3D[points], DirectedEdges -> False, opt, VertexCoordinates -> points],
      _, Message[IGDelaunayGraph::dim]; throw[$Failed]
    ]


(* The following beta-skeleton computation is based on the implementation of Henrik Schumacher
   https://mathematica.stackexchange.com/a/183391/12 *)

(* Notes:

   - The circle-based and lune-based beta skeletons are distinct.
   - For beta <= 1, the two definitions coincide.
   - For beta >= 1, the circle-based beta skeleton is a subgraph of the lune-based beta skeleton,
      which is a subgraph of the Gabriel graph, which is a subgraph of the Delaunay graph.
   - For beta = 1, the beta skeleton coincides with the Gabriel graph (either definition)
   - For beta = 2, the lune-based beta skeleton coincides with the relative neighborhood graph
 *)

IGraphM::bsdim2 = "Beta skeleton computation is only supported in 2 dimensions for beta < 1 or circle-based beta skeletons.";
IGraphM::bsdim3 = "Beta skeleton computation is only supported in 2 or 3 dimensions.";
IGraphM::bsdupl = "Duplicate points must be removed before beta skeleton computations.";

(*
betaSkeletonEdgeSuperset[pts_, beta_ /; beta >= 1] :=
    Module[{mesh, edges, edgeLengths, p, q},
      If[Not@MatchQ[Dimensions[pts], {_,2}],
        Message[IGraphM::bsdim];
        throw[$Failed]
      ];
      mesh = DelaunayMesh[pts];
      If[Head[mesh] === MeshRegion
        ,
        If[MeshCellCount[mesh, 0] =!= Length[pts],
          Message[IGraphM::bsdupl];
          throw[$Failed]
        ];
        edges = MeshCells[mesh, 1, "Multicells" -> True][[1, 1]];
        edgeLengths = PropertyValue[{mesh, 1}, MeshCellMeasure];
        {p, q} = pts[[#]]& /@ Transpose[edges];
        ,
        edges = Subsets[Range@Length[pts], {2}];
        {p, q} = pts[[#]]& /@ Transpose[edges];
        edgeLengths = Sqrt[Dot[Subtract[p, q]^2, ConstantArray[1., 2]]];
      ];
      {edges, edgeLengths, p, q}
    ]
*)

betaSkeletonEdgeSuperset[pts_, beta_ /; beta >= 1] :=
    Module[{edges, edgeLengths, p, q, dim},
      Switch[Dimensions[pts],
        {_, 2},
        edges = check@delaunayEdges2D[pts];
        dim = 2;
        ,
        {_, 3},
        edges = check@delaunayEdges3D[pts];
        dim = 3;
        ,
        _,
        Message[IGraphM::bsdim3];
        throw[$Failed]
      ];
      {p, q} = pts[[#]]& /@ Transpose[edges];
      edgeLengths = Sqrt@Dot[Subtract[p, q]^2, ConstantArray[1., dim]];
      {edges, edgeLengths, p, q}
    ]

betaSkeletonEdgeSuperset[pts_, beta_ /; beta < 1] :=
    Module[{edges, edgeLengths, p, q},
      If[Not@MatchQ[Dimensions[pts], {_,2}],
        Message[IGraphM::bsdim2];
        throw[$Failed]
      ];
      edges = Subsets[Range@Length[pts], {2}];
      {p, q} = pts[[#]]& /@ Transpose[edges];
      edgeLengths = Sqrt[Dot[Subtract[p, q]^2, ConstantArray[1., 2]]];
      {edges, edgeLengths, p, q}
    ]

(* beta >= 1, lune-based *)
(*
igLuneBetaSkeletonEdges[pts_, beta_] :=
    Module[{nf, edges, edgeLengths, p, q, r, centres1, centres2},
      nf = Nearest[pts -> Automatic];

      {edges, edgeLengths, p, q} = betaSkeletonEdgeSuperset[pts, beta];

      r = 0.5 beta;

      centres1 = p + (r-1) (p-q);
      centres2 = q + (r-1) (q-p);

      Pick[
        edges,
        MapThread[
          Function[{c1, c2, d}, Length@Intersection[nf[c1, {Infinity, d}], nf[c2, {Infinity, d}]]],
          {centres1, centres2, r edgeLengths (1 + 10^Internal`$EqualTolerance $MachineEpsilon)}
        ],
        2
      ]
    ]
*)
igLuneBetaSkeletonEdges[pts_, beta_] :=
    Module[{flann, edges, edgeLengths, p, q, r, centres1, centres2, dists},
      {edges, edgeLengths, p, q} = betaSkeletonEdgeSuperset[pts, beta];

      r = 0.5 beta;

      centres1 = p + (r-1) (p-q);
      centres2 = q + (r-1) (q-p);

      (* Increase distances by the standard relative tolerance to include boundaries. *)
      dists = r edgeLengths (1 + 10^Internal`$EqualTolerance $MachineEpsilon);

      (* intersectionCounts[] excludes the edge endpoints, thus we only need to check
         that the intersection is empty. *)
      flann = makeFlann[pts];
      Pick[
        edges,
        Unitize[flann@"intersectionCounts"[centres1, centres2, dists, edges, True]],
        0
      ]
    ]


(* The relative neighbourhood graph is defined in terms of an open exclusion region,
 * i.e. points on the region boundaries are not considered. *)
igRelativeNeighborhoodGraphEdges[pts_] :=
    Module[{flann, edges, edgeLengths, p, q, dists},
      {edges, edgeLengths, p, q} = betaSkeletonEdgeSuperset[pts, 2];

      (* Decrease distances by the standard relative tolerance to include boundaries. *)
      dists = edgeLengths (1 - 10^Internal`$EqualTolerance $MachineEpsilon);

      (* intersectionCounts[] excludes the edge endpoints, thus we only need to check
         that the intersection is empty. *)
      flann = makeFlann[pts];
      Pick[
        edges,
        Unitize[flann@"intersectionCounts"[p, q, dists, edges, True]],
        0
      ]
    ]


(* beta >= 1, circle-based *)
(*
igCircleBetaSkeletonEdges[pts_, beta_] :=
    Module[{nf, edges, edgeLengths, p, q, r, centres1, centres2},
      nf = Nearest[pts -> Automatic];

      {edges, edgeLengths, p, q} = betaSkeletonEdgeSuperset[pts, beta];

      r = 0.5 beta;

      With[{mid = 0.5 (p+q), perp = Sqrt[r^2 - 0.25] RotationTransform[Pi/2][p-q]},
        centres1 = mid + perp;
        centres2 = mid - perp;
      ];

      Pick[
        edges,
        MapThread[
          Function[{c1, c2, d}, Length@Union[nf[c1, {Infinity, d}], nf[c2, {Infinity, d}]]],
          {centres1, centres2, r edgeLengths (1 + 10^Internal`$EqualTolerance $MachineEpsilon)}
        ],
        2
      ]
    ]
*)
igCircleBetaSkeletonEdges[pts_, beta_] :=
    Module[{flann, edges, edgeLengths, p, q, r, centres1, centres2, dists},
      If[Not@MatchQ[Dimensions[pts], {_, 2}],
        Message[IGraphM::bsdim2];
        throw[$Failed]
      ];

      {edges, edgeLengths, p, q} = betaSkeletonEdgeSuperset[pts, beta];

      r = 0.5 beta;

      With[{mid = 0.5 (p+q), perp = Sqrt[r^2 - 0.25] RotationTransform[Pi/2][p-q]},
        centres1 = mid + perp;
        centres2 = mid - perp;
      ];

      (* Increase distances by the standard relative tolerance to include boundaries. *)
      dists = r edgeLengths (1 + 10^Internal`$EqualTolerance $MachineEpsilon);

      flann = makeFlann[pts];
      Pick[
        edges,
        Unitize[flann@"unionCounts"[centres1, centres2, dists, edges, True]],
        0
      ]
    ]

(* beta = 1, both lune and circle *)
(*
igGabrielGraphEdges[pts_] :=
    Module[{nf, edges, edgeLengths, p, q},
      nf = Nearest[pts -> Automatic];

      {edges, edgeLengths, p, q} = betaSkeletonEdgeSuperset[pts, 1];

      Pick[
        edges,
        MapThread[
          Function[{c, d}, Length@nf[c, {Infinity, d}]],
          {(p+q)/2, 0.5 edgeLengths (1 + 10^Internal`$EqualTolerance $MachineEpsilon)}
        ],
        2
      ]
    ]
*)
igGabrielGraphEdges[pts_] :=
    Module[{flann, edges, edgeLengths, p, q, dists},
      {edges, edgeLengths, p, q} = betaSkeletonEdgeSuperset[pts, 1];

      (* Increase distances by the standard relative tolerance to include boundaries. *)
      dists = 0.5 edgeLengths (1 + 10^Internal`$EqualTolerance $MachineEpsilon);

      flann = makeFlann[pts];
      Pick[
        edges,
        Unitize[flann@"neighborCounts"[(p+q)/2, dists, edges, True]],
        0
      ]
    ]

(* 0 < beta < 1, both lune and circle *)
igBetaSkeletonEdges0[pts_, beta_] :=
    Module[{flann, edges, edgeLengths, p, q, r, centres1, centres2, dists},
      If[Not@MatchQ[Dimensions[pts], {_, 2}],
        Message[IGraphM::bsdim2];
        throw[$Failed]
      ];

      {edges, edgeLengths, p, q} = betaSkeletonEdgeSuperset[pts, beta];

      r = 0.5 / beta;

      With[{mid = 0.5 (p+q), perp = Sqrt[r^2 - 0.25] RotationTransform[Pi/2][p-q]},
        centres1 = mid + perp;
        centres2 = mid - perp;
      ];

      (* Increase distances by the standard relative tolerance to include boundaries. *)
      dists = r edgeLengths (1 + 10^Internal`$EqualTolerance $MachineEpsilon);

      flann = makeFlann[pts];
      Pick[
        edges,
        Unitize[flann@"intersectionCounts"[centres1, centres2, dists, edges, True]],
        0
      ]
    ]


PackageExport["IGLuneBetaSkeleton"]
IGLuneBetaSkeleton::usage = "IGLuneBetaSkeleton[points, beta] gives the lune-based beta skeleton of the given points.";

(* Note: pts is numericized with N[] to avoid crash in M10.0.2.  M10.3 does not crash. *)
igLuneBetaSkeleton[pts_, beta_, opt___] :=
    catch@If[Length[pts] < 2, IGEmptyGraph[Length[pts], opt],
      With[{
        edges = Which[
          beta  > 1, igLuneBetaSkeletonEdges[N[pts], beta],
          beta == 1, igGabrielGraphEdges[N[pts]],
          beta  < 1, igBetaSkeletonEdges0[N[pts], beta]
        ]
      },
        Graph[Range@Length[pts], edges, DirectedEdges -> False, opt, VertexCoordinates -> pts]
      ]
    ]

SyntaxInformation[IGLuneBetaSkeleton] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}, "OptionNames" -> optNames[Graph]};
IGLuneBetaSkeleton[pts : {} | _?(MatrixQ[#, NumericQ]&), beta_?positiveNumericQ, opt : OptionsPattern[Graph]] :=
    igLuneBetaSkeleton[pts, beta, opt]


PackageExport["IGCircleBetaSkeleton"]
IGCircleBetaSkeleton::usage = "IGCircleBetaSkeleton[points, beta] gives the circle-based beta skeleton of the given points.";

(* Note: pts is numericized with N[] to avoid crash in M10.0.2.  M10.3 does not crash. *)
igCircleBetaSkeleton[pts_, beta_, opt___] :=
    catch@If[Length[pts] < 2, IGEmptyGraph[Length[pts], opt],
      With[{
        edges = Which[
          beta >  1, igCircleBetaSkeletonEdges[N[pts], beta],
          beta == 1, igGabrielGraphEdges[N[pts]],
          beta <  1, igBetaSkeletonEdges0[N[pts], beta]
        ]
      },
        Graph[Range@Length[pts], edges, DirectedEdges -> False, opt, VertexCoordinates -> pts]
      ]
    ]

SyntaxInformation[IGCircleBetaSkeleton] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}, "OptionNames" -> optNames[Graph]};
IGCircleBetaSkeleton[pts : {} | _?(MatrixQ[#, NumericQ]&), beta_?positiveNumericQ, opt : OptionsPattern[Graph]] :=
    igCircleBetaSkeleton[pts, beta, opt]


PackageExport["IGRelativeNeighborhoodGraph"]
IGRelativeNeighborhoodGraph::usage = "IGRelativeNeighborhoodGraph[points] gives the relative neighbourhood graph of the given points.";

SyntaxInformation[IGRelativeNeighborhoodGraph] = {"ArgumentsPattern" -> {_, OptionsPattern[]}, "OptionNames" -> optNames[Graph]};
IGRelativeNeighborhoodGraph[pts : {} | _?(MatrixQ[#, NumericQ]&), opt : OptionsPattern[Graph]] :=
    catch@If[Length[pts] < 2, IGEmptyGraph[Length[pts], opt],
      With[{ edges = igRelativeNeighborhoodGraphEdges[N[pts]] },
        Graph[Range@Length[pts], edges, DirectedEdges -> False, opt, VertexCoordinates -> pts]
      ]
    ]


PackageExport["IGGabrielGraph"]
IGGabrielGraph::usage = "IGGabrielGraph[points] gives the Gabriel graph of the given points.";

SyntaxInformation[IGGabrielGraph] = {"ArgumentsPattern" -> {_, OptionsPattern[]}, "OptionNames" -> optNames[Graph]};
IGGabrielGraph[pts : {} | _?(MatrixQ[#, NumericQ]&), opt : OptionsPattern[Graph]] :=
    igLuneBetaSkeleton[pts, 1, opt]


PackageExport["IGBetaWeightedGabrielGraph"]
IGBetaWeightedGabrielGraph::usage = "IGBetaWeightedGabrielGraph[points]";

SyntaxInformation[IGBetaWeightedGabrielGraph] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}, "OptionNames" -> optNames[Graph]};
IGBetaWeightedGabrielGraph[pts : {} | _?(MatrixQ[#, NumericQ]&), beta : _?Positive : Infinity, opt : OptionsPattern[Graph]] :=
    catch@Module[{edges, flann, betas, mask},
      Switch[Dimensions[pts],
        {_, 2},
        edges = check@delaunayEdges2D[pts];
        dim = 2;
        ,
        {_, 3},
        edges = check@delaunayEdges3D[pts];
        dim = 3;
        ,
        _,
        Message[IGraphM::bsdim3];
        throw[$Failed]
      ];

      flann = makeFlann[pts];

      betas = fixInfNaN@expectInfNaN@check@flann@"edgeBetas"[edges, infToNeg[beta]];
      mask = Unitize[betas];

      Graph[
        Range@Length[pts],
        Pick[edges, mask, 1],
        DirectedEdges -> False,
        opt,
        EdgeWeight -> Pick[betas, mask, 1],
        VertexCoordinates -> pts
      ]
    ]
