(* Mathematica Package *)
(* Created by Mathematica plugin for IntelliJ IDEA *)

(* :Author: szhorvat *)
(* :Date: 2018-10-24 *)
(* :Copyright: (c) 2018 Szabolcs Horvát *)

Package["IGraphM`"]

(******************************)
(***** Geometrical meshes *****)
(******************************)


(***** Converting meshes to graphs or matrices *****)

meshQ[_?MeshRegionQ] := True
meshQ[_?BoundaryMeshRegionQ] := True
meshQ[_] := False


PackageExport["IGMeshGraph"]
IGMeshGraph::usage = "IGMeshGraph[mesh] converts the edges and vertices of a geometrical mesh to a weighted graph.";

IGMeshGraph::noprop = "The edge property `1` is not present in the mesh.";

Options[IGMeshGraph] = { EdgeWeight -> MeshCellMeasure };
SyntaxInformation[IGMeshGraph] = {"ArgumentsPattern" -> {_, OptionsPattern[]}, "OptionNames" -> optNames[IGMeshGraph, Graph]};
IGMeshGraph[mesh_?meshQ, opt : OptionsPattern[{IGMeshGraph, Graph}]] :=
    Module[{edgeWeightRule, ew = OptionValue[EdgeWeight], pv},
      edgeWeightRule = Switch[ew,
        None,
        Unevaluated@Sequence[]
        ,
        _String | _Symbol,
        pv = PropertyValue[{mesh, 1}, ew];
        If[pv === $Failed,
          Message[IGMeshGraph::noprop, ew];
          Return[$Failed]
        ];
        EdgeWeight -> pv
        ,
        _List,
        EdgeWeight -> ew
        ,
        _,
        Message[IGMeshGraph::invopt, ew, EdgeWeight, OptionValue[IGMeshGraph, EdgeWeight]];
        EdgeWeight -> PropertyValue[{mesh, 1}, MeshCellMeasure]
      ];
      Graph[
        (* We return an index graph. *)
        Range@MeshCellCount[mesh, 0],
        (* The undocumented option "Multicells" -> True causes the result to come in the form
         * {Line[{{1,2},{2,3}}]} instead of {Line[{1,2}], Line[{2,3}]} *)
        Join @@ (MeshCells[mesh, 1, "Multicells" -> True][[All, 1]]),
        edgeWeightRule,
        Sequence @@ FilterRules[FilterRules[{opt}, Options[Graph]], Except@Options[IGMeshGraph]],
        VertexCoordinates -> MeshCoordinates[mesh]
      ]
    ]



(* Thanks to Henrik Schumacher for the following set of cell adjacency functions!
   https://mathematica.stackexchange.com/a/160457/12 *)

igMeshCellAdjacencyMatrix[mesh_, 0, 0] :=
    With[{edges = Developer`ToPackedArray[MeshCells[mesh, 1][[All, 1]]]},
      SparseArray[
        Rule[
          Join[edges, Reverse[edges, {2}]],
          ConstantArray[1, 2 Length[edges]]
        ],
        {MeshCellCount[mesh, 0], MeshCellCount[mesh, 0]}
      ]
    ]

igMeshCellAdjacencyMatrix[mesh_, d_, 0] :=
    Module[{pts, cells, A, lens, nn},
      pts = MeshCoordinates[mesh];
      (* The Check[] below is meant to catch errors such as those resulting from
         trying to get 3D mesh cells of a BoundaryDiscretizeRegion[Ball[]]. *)
      cells = Developer`ToPackedArray[ Check[MeshCells[mesh, d], throw[$Failed]][[All, 1]] ];
      lens = Length /@ cells;
      nn = Total[lens];
      A = SparseArray @@ {Automatic, {Length[cells], Length[pts]},
        0, {1, {Developer`ToPackedArray[Join[{0}, Accumulate[lens]]],
          ArrayReshape[Flatten[Sort /@ cells], {nn, 1}]},
          ConstantArray[1, nn]}}
    ]

igMeshCellAdjacencyMatrix[mesh_, 0, d_] :=
    Transpose[igMeshCellAdjacencyMatrix[mesh, d, 0]]

igMeshCellAdjacencyMatrix[mesh_, d1_, d2_] :=
    With[{B = igMeshCellAdjacencyMatrix[mesh, d1, 0].igMeshCellAdjacencyMatrix[mesh, 0, d2]},
      SparseArray[
        If[d1 == d2,
          UnitStep[B - DiagonalMatrix[Diagonal[B]] - d1],
          UnitStep[B - (Min[d1, d2] + 1)]
        ]
      ]
    ]

checkDimension[dim_, d_, sym_] :=
    If[d > dim,
      Message[sym::bddim, d, dim];
      throw[$Failed]
    ]


PackageExport["IGMeshCellAdjacencyGraph"]
IGMeshCellAdjacencyGraph::usage =
    "IGMeshCellAdjacencyGraph[mesh, d] returns the connectivity structure of d-dimensional cells in mesh as a graph.\n" <>
    "IGMeshCellAdjacencyGraph[mesh, d1, d2] returns the connectivity structure of d1 and d2 dimensional cells in mesh as a bipartite graph.";

IGMeshCellAdjacencyGraph::bddim = "The requested dimension, `1`, is greater than the dimension of the mesh, `2`.";

Options[IGMeshCellAdjacencyGraph] = { VertexCoordinates -> None };
SyntaxInformation[IGMeshCellAdjacencyGraph] = {"ArgumentsPattern" -> {_, _, _., OptionsPattern[]}, "OptionNames" -> optNames[IGMeshCellAdjacencyGraph, Graph]};
IGMeshCellAdjacencyGraph[mesh_?meshQ, d1_?Internal`NonNegativeIntegerQ, d2_?Internal`NonNegativeIntegerQ, opt : OptionsPattern[{IGMeshCellAdjacencyGraph, Graph}]] :=
    catch@Module[{coord, dim = RegionDimension[mesh]},
      With[{dim = dim},
        Scan[checkDimension[dim, #, IGMeshCellAdjacencyGraph]&, {d1, d2}]
      ];
      Switch[OptionValue[VertexCoordinates],
        None,
        coord = Automatic
        ,
        Automatic,
        coord =
            If[dim < 2 || dim > 3,
              Automatic, (* can't use coordinates from less than 2 or more than 3 dimensions *)
              Check[
                If[d1 == d2,
                  PropertyValue[{mesh, {d1, All}}, MeshCellCentroid],
                  Join[
                    PropertyValue[{mesh, {d1, All}}, MeshCellCentroid],
                    PropertyValue[{mesh, {d2, All}}, MeshCellCentroid]
                  ]
                ],
                Automatic (* if the calculation failed, use automatic coordinates *)
              ]
            ]
        ,
        _,
        coord = OptionValue[VertexCoordinates]
      ];
      If[d1 == d2,
        AdjacencyGraph[
          MeshCellIndex[mesh, d1], igMeshCellAdjacencyMatrix[mesh, d1, d2],
          VertexCoordinates -> coord, FilterRules[{opt}, Except[VertexCoordinates]]
        ],
        IGBipartiteIncidenceGraph[
          {MeshCellIndex[mesh, d1], MeshCellIndex[mesh, d2]}, igMeshCellAdjacencyMatrix[mesh, d1, d2],
          VertexCoordinates -> coord, FilterRules[{opt}, Except[VertexCoordinates]]
        ]
      ]
    ]
IGMeshCellAdjacencyGraph[mesh_?meshQ, d_?Internal`NonNegativeIntegerQ, opt : OptionsPattern[Graph]] :=
    IGMeshCellAdjacencyGraph[mesh, d, d, opt]


PackageExport["IGMeshCellAdjacencyMatrix"]
IGMeshCellAdjacencyMatrix::usage =
    "IGMeshCellAdjacencyMatrix[mesh, d] returns the adjacency matrix of d-dimensional cells in mesh.\n" <>
    "IGMeshCellAdjacencyMatrix[mesh, d1, d2] returns the incidence matrix of d1- and d2-dimensional cells in mesh.";

IGMeshCellAdjacencyMatrix::bddim = IGMeshCellAdjacencyGraph::bddim;
SyntaxInformation[IGMeshCellAdjacencyMatrix] = {"ArgumentsPattern" -> {_, _, _.}};
IGMeshCellAdjacencyMatrix[mesh_?meshQ, d1_?Internal`NonNegativeIntegerQ, d2_?Internal`NonNegativeIntegerQ] :=
    catch@With[{dim = RegionDimension[mesh]},
      Scan[checkDimension[dim, #, IGMeshCellAdjacencyMatrix]&, {d1, d2}];
      igMeshCellAdjacencyMatrix[mesh, d1, d2]
    ]
IGMeshCellAdjacencyMatrix[mesh_?meshQ, d_?Internal`NonNegativeIntegerQ] :=
    IGMeshCellAdjacencyMatrix[mesh, d, d]


(***** Mesh generation *****)

PackageExport["IGLatticeMesh"]
IGLatticeMesh::usage =
    "IGLatticeMesh[type] creates a mesh of the lattice of the specified type.\n" <>
    "IGLatticeMesh[type, {m, n}] creates a lattice of n by m unit cells.\n" <>
    "IGLatticeMesh[type, region] creates a lattice from the points that fall within region.\n" <>
    "IGLatticeMesh[] returns a list of available lattice types.";

$igLatticeData := $igLatticeData = zimport@FileNameJoin[{$packageDirectory, "IGLatticeData.mz"}];
$igLatticeUnits := $igLatticeUnits = $igLatticeData["UnitCells"];
$igLatticeVectors := $igLatticeVectors = $igLatticeData["TranslationVectors"];

IGLatticeMesh::noval = "`1` is not a know lattice type. Evaluate IGLatticeMesh[] to see a list of valid lattice types.";
IGLatticeMesh::regcst = "The second argument must be a constant (parameter free) region.";
IGLatticeMesh::regdim = "The second argument must be a 2-dimensional region.";
IGLatticeMesh::regbnd = "The second argument must be a bounded region.";
IGLatticeMesh::regemp = "The given region does not contain any lattice points.";
SyntaxInformation[IGLatticeMesh] = {"ArgumentsPattern" -> {_., _., OptionsPattern[]}, "OptionNames" -> optNames[MeshRegion]};
IGLatticeMesh[] := Keys[$igLatticeUnits]
IGLatticeMesh[name_String, dims : {_?Internal`PositiveIntegerQ, _?Internal`PositiveIntegerQ} : {7, 7}, opt : OptionsPattern[]] :=
      catch@Module[{m, n, grid, polys, pts, newpts, nf, ratio = 0},
        If[Not@KeyExistsQ[$igLatticeUnits, name],
          Message[IGLatticeMesh::noval, name];
          throw[$Failed]
        ];
        {m, n} = dims;
        grid = Join @@ Table[{i - Floor[ratio j], j}, {i, 0, m-1}, {j, 0, n-1}];
        igLatticeMesh[$igLatticeUnits[name], $igLatticeVectors[name], grid.$igLatticeVectors[name], {opt}]
      ]
IGLatticeMesh[name_String, reg_?RegionQ, opt : OptionsPattern[]] :=
    catch@Module[{boundPoints, xy, grid, trpts},
      If[Not@KeyExistsQ[$igLatticeUnits, name],
        Message[IGLatticeMesh::noval, name];
        throw[$Failed]
      ];
      If[Not@ConstantRegionQ[reg],
        Message[IGLatticeMesh::regcst];
        throw[$Failed]
      ];
      If[RegionEmbeddingDimension[reg] != 2,
        Message[IGLatticeMesh::regdim];
        throw[$Failed]
      ];
      If[Not@BoundedRegionQ[reg],
        Message[IGLatticeMesh::regbnd];
        throw[$Failed]
      ];
      If[RegionDimension[reg] == -Infinity, (* empty region *)
        Message[IGLatticeMesh::regemp];
        throw[$Failed]
      ];
      boundPoints = Tuples@Check[RegionBounds[reg, "Sufficient"], throw[$Failed]];
      xy = Transpose[boundPoints.Inverse[$igLatticeVectors[name]]];
      grid = Tuples[Range @@@ Transpose@{Floor[Min /@ xy], Ceiling[Max /@ xy]}];
      trpts = Select[grid.$igLatticeVectors[name], RegionMember[reg]];
      If[trpts === {},
        Message[IGLatticeMesh::regemp];
        throw[$Failed]
      ];
      igLatticeMesh[$igLatticeUnits[name], $igLatticeVectors[name], trpts, {opt}]
    ]
igLatticeMesh[unit_, vec_, trpts_, {opt___}] :=
    Module[{polys, pts, newpts, nf},
      polys = Flatten[
        Function[tr, Replace[unit, pt_ :> pt + tr, {2}]] /@ trpts,
        1
      ];
      pts = Join @@ polys;
      newpts = DeleteDuplicatesBy[pts, Round[N[#], 1*^-6] &];
      newpts = newpts[[ Ordering[ newpts.RotationTransform[-Pi/2 + 0.01][ vec[[2]] ] ] ]];
      nf = Nearest[N[newpts] -> Range@Length[newpts]];
      MeshRegion[
        newpts,
        Polygon@Flatten[Map[nf, N[polys], {2}], {{1}, {2, 3}}],
        opt
      ]
    ]
addCompletion[IGLatticeMesh, {IGLatticeMesh[]}];
