## Revision history for [IGraph/M](README.md)

#### Unreleased

Added:

  - `IGHarmonicCentrality` and `IGHarmonicCentralityCutoff` compute the harmonic centrality and range-limited harmonic centrality.
  - `IGLinkRank` and `IGPersonalizedLinkRank` compute the equivalent of PageRank for edges.
  - `IGNeighborhoodCloseness` computes the range-limited closeness centrality, as well as the number of vertices reachable within the given range.
  - `IGFamousGraph` exposes the igraph C library's built-in graph database.
  - `IGPreferenceGame` and `IGAsymmetricPreferenceGame` create non-growing random graphs based on vertex types.
  - `IGReingoldTilford` and `IGReingoldTilfordCircular` now support the `DirectedEdges` option.
  - `IGFruchtermanReingold` now supports constraining the coordinates of a subset of vertices.
  - `IGPercolationCurve` for efficiently computing the size of the largest component as a function of mean degree while removing edges.
  - `IGShortestPathTree` for computing a shortest path tree rooted in a given vertex.
  - `IGGraphEditor` is an experimental interactive graph editor.
  - Experimental progress reporting functionality through functions in the ``IGraphM`Progress` `` context.

Changed:

 - `IGConnectedQ` and `IGWeaklyConnectedQ` now consider the null graph to be disconnected; this is consistent with other functions such as `IGTreeQ`.
 - `IGAveragePathLength` now has a `"ByComponents"` option, controlling the handling of disconnected graphs.
 - Centrality functions:
    * `IGCloseness` now computes the normalized closeness, i.e. the inverse of the mean distance to other vertices, by default. Use `Normalized -> False` to get the previous behaviour.
    * `IGCloseness` now uses the distances only to reachable vertices when computing the closeness. In undirected disconnected graphs, it effectively computes the closeness per component. For isolated vertices (or sinks in directed graphs) it now returns `Indeterminate`.
    * `IGBetweenness` and `IGBetweennessCentralization` no longer uses the `Method` option. The calculations are always fast and precise. The precision has been improved.
    * `IGBetweennessEstimate`, `IGEdgeBetweennessEstimate` and `IGClosenessEstimate` have been renamed to `IGBetweennessCutoff`, `IGEdgeBetweennessCutoff` and `IGClosenessCutoff`.
 - `IGRelativeNeighborhoodGraph` now assumes the `β -> 2`, `β < 2` limit instead of `β = 2`.
 - `IGGirth` now returns `Infinity` for the null graph.
 - `IGDiameter` now returns `Indeterminate` for the null graph.
 - `IGChordalQ`, `IGChordalCompletion` and `IGMaximumCardinalitySearch` now support non-simple graphs.
 - `IGReingoldTilford` and `IGReingoldTilfordCircular` use a new automatic root selection algorithm. The root selection heuristic may change in the future without notice. Specify roots manually for a consistent result.
 - `IGPotentiallyConnectedQ` no longer supports directed sequences. This feature was flawed in 0.5. It may be re-added in a future version.
 - `IGLayoutKamadaKawai` and `IGLayoutKamadaKawai3D` perform more iterations by default, and produce more pleasing layouts.
 - LAD isomorphism functions now support self-loops.
 - Motif finder functions now support size 5 and 6 undirected motifs.
 - The behaviour of the random number generator is now consistent between platforms, meaning that with a given seed, a randomized IGraph/M function will return the same result on all platforms. However, result will now be different from version 0.5 for the same seed.

Fixed:

 - `IGPageRank` and `IGPersonalizedPageRank` will now warn if the calculation did not converge with the `"Arnoldi"` method. This happens only in rare cases.
 - `IGPersonalizedPageRank`: the default `"PRPACK"` method returned an incorrect result when the graph was not connected and the personalization vector was not uniform.
 - Several community detection functions did not handle zero or one-vertex graphs correctly.
 - `IGVertexMap` would evaluate the mapped functions twice instead of once.
 - `IGMaximumCardinalitySearch` returned incorrect ranks for graphs whose vertex names differed from their vertex indices.
 - `IGDistanceWeighted` no longer fails on edgeless graphs.
 - `IGMotifsVertexParticipation` would fail with Mathematica 12.2 or later.
 - `IGCallawayTraisGame` no longer rejects zeros in the preference matrix.
 - An error would be triggered when some functions returned a zero-size matrix.
 - Fixed a memory leak in the Nauty format reader.
 - Fixed a conflict with Mathematica 13.0's built-in isomorphism functions which could lead to a crash on Linux.

Other:

 - IGraph/M now requires Mathematica 11.0 or later; on the Raspberry Pi it requires Wolfram Engine 12.2 or later.
 - More robust error handling: when certain serious errors occur in the igraph C library, the Mathematica kernel is no longer forced to shut down.
 - IGraph/M got leaner: smaller binary sizes.
 - Documentation improvements.

#### v0.5.1

Changes to existing functions:

 - `IGReverseGraph` now supports non-simple graphs.
 - `IGLuneBetaSkeleton` is now interruptible for β < 1.

Fixes:

 - Several functions failed to check if the graph was mixed (having both directed and undirected edges). This is now corrected.
 - `IGLayoutBipartite` now uses partitions which are consistent with `IGBipartitePartitions`.

Other changes:

 - Many corrections and improvements to the documentation.

#### v0.5

New functions:

 - `IGSplitQ` recognizes split graphs and their degree sequences.
 - `IGThresholdQ` recognizes threshold graphs and their degree sequences.
 - `IGPotentiallyConnectedQ` recognizes the degree sequences of connected graphs.
 - `IGBigraphicalQ` recognizes the degree sequence pairs of bipartite graphs.
 - `IGEulerianQ` tests if a graph has an Eulerian path; `IGEulerianPath` and `IGEulerianPathVertices` find it.

Changes to existing functions:

 - `IGGraphicalQ` and `IGRealizeDegreeSequence` now support some types of non-simple graphs.
 - `IGRealizeDegreeSequence` now takes arguments in the order `indegrees`, `outdegrees` for consistency with other functions (previously it was `outdegrees`, `indegrees`).
 - `IGEigenvectorCentrality` sometimes returned incorrect values for isolated vertices in weighted graphs.
 - `IGColoredSimpleGraph` no longer discards vertex names.
 - `IGModularity` now supports directed graphs.
 - `IGModularity` and `IGCommunitiesMultilevel` now have a resolution parameter.
 - `IGAdjacencyMatrixPlot` now allows `None` to be specified as the colour representing non-existing edges.
 - `IGEigenvectorCentrality` assumes the adjacency matrix of undirected graphs to have twice the number of self-loops for each vertex on the diagonal. This makes the results consistent between an undirected graph and its directed equivalent when each edge is replaced by a mutual edge pair.

Fixes:

 - `IGLayoutReingoldTilford` no longer flips the layout.
 - `IGLayoutReingoldTilford` no longer draws overlapping tree branches.
 - `IGBarabasiAlbertGame` now allows negative values for β.
 - `IGBetweennessEstimate` sometimes returned incorrect results with finite cutoffs. This is now corrected.
 - `IGCommunitiesLeiden`: fix incorrect results when self-loops are present.
 - `IGEigenvectorCentrality`: fix incorrect results for isolated vertices and for vertices with self-loops.
 - `IGGraphicalQ` would return incorrect results when the second argument was `{}`. This is now corrected.

Other changes:

 - `IGDegreeSequenceGame`'s `"ConfigurationModelSimple"` method is now much faster.
 - More robust error handling.

Notes:

 - IGraph/M 0.5 has been tested with Mathematica 10.3 and later only. Some effort has been made to allow it to work with Mathematica 10.0.2, but it has not been tested and compatibility is not guaranteed.

#### v0.4

New functions:

 - Deterministic graph generators: `IGKautzGraph`, `IGCompleteGraph`, `IGCompleteAcyclicGraph`, `IGDeBruijnGraph`, `IGChordalRing`, `IGEmptyGraph`, `IGRealizeDegreeSequence`, `IGFromPrufer`, `IGToPrufer`, `IGKaryTree`, `IGSymmetricTree`, `IGBetheLattice`, `IGTriangularLattice`, `IGMycielskian`, `IGExpressionTree`, `IGShorthand`, `IGFromNauty`.
 - Random graph generators: `IGWattsStrogatzGame`, `IGCallawayTraitsGame`, `IGEstablishmentGame`, `IGTreeGame`, `IGErdosRenyiGameGNM`, `IGErdosRenyiGameGNP`.
 - Weighted graph functions: `IGWeightedSimpleGraph`, `IGWeightedUndirectedGraph`, `IGWeightedVertexDelete`, `IGWeightedSubgraph`, `IGUnweighted`, `IGDistanceWeighted`, `IGWeightedAdjacencyGraph`, `IGVertexWeightedQ`, `IGEdgeWeightedQ`, `IGVertexStrength`, `IGVertexInStrength`, `IGVertexOutStrength`.
 - Community detection: `IGCommunitesFluid`, `IGCommunitiesLeiden`.
 - Graph colouring functions: `IGVertexColoring`, `IGEdgeColoring`, `IGKVertexColoring`, `IGKEdgeColoring`, `IGMinimumVertexColoring`, `IGMinimumEdgeColoring`, `IGChromaticNumber`, `IGChromaticIndex`, `IGVertexColoringQ`.
 - Clique cover: `IGCliqueCover`, `IGCliqueCoverNumber`.
 - Mesh/graph conversion: `IGMeshGraph`, `IGMeshCellAdjacencyMatrix`, `IGMeshCellAdjacencyGraph`.
 - Lattice generation: `IGLatticeMesh`, `IGTriangularLattice`.
 - Proximity graphs: `IGDelaunayGraph`, `IGGabrielGraph`, `IGRelativeNeighborhoodGraph`, `IGLuneBetaSkeleton`, `IGCircleBetaSkeleton`.
 - Centralization: `IGDegreeCentralization`, `IGBetweennessCentralization`, `IGClosenessCentralization`, `IGEigenvectorCentralization`.
 - Planar graphs: `IGPlanarQ`, `IGMaximalPlanarQ`, `IGOuterplanarQ`, `IGKuratowskiEdges`, `IGFaces`, `IGDualGraph`, `IGEmbeddingQ`, `IGPlanarEmbedding`, `IGOuterplanarEmbedding`, `IGCoordinatesToEmbedding`, `IGEmbeddingToCoordinates`, `IGLayoutPlanar`, `IGLayoutTutte`.
 - Adjacency lists and embeddings: `IGAdjacencyList`, `IGAdjacencyGraph`.
 - Spanning trees and other tree-related functionality: `IGSpanningTree`, `IGRandomSpanningTree`, `IGSpanningTreeCount`, `IGUnfoldTree`, `IGTreeQ`, `IGForestQ`, `IGTreelikeComponents`, `IGTreeGame`, `IGStrahlerNumber`, `IGOrientTree`.
 - Matching functions: `IGMaximumMatching`, `IGMatchingNumber`.
 - Dominance: `IGDominatorTree`, `IGImmediateDominators`
 - Maximum flow: `IGMaximumFlowValue`, `IGMaximumFlowMatrix`.
 - Isomorphism: `IGGetIsomorphism` and `IGGetSubisomorphism` (they work with multigraphs), `IGColoredSimpleGraph` (for transforming multigraph isomorphism to coloured graph isomorphism).
 - Transitivity: `IGVertexTransitiveQ`, `IGEdgeTransitiveQ`, `IGSymmetricQ`, `IGDistanceTransitiveQ`.
 - Regular graphs: `IGRegularQ`, `IGStronglyRegularQ`, `IGStronglyRegularParameters`, `IGDistanceRegularQ`, `IGIntersectonArray`.
 - A framework for easy property transformations and graph styling: `IGVertexProp`, `IGEdgeProp`, `IGEdgeVertexProp`, `IGVertexMap`, `IGEdgeMap`, `IGVertexPropertyList`, `IGEdgePropertyList`.
 - Import functions: `IGImport`, `IGImportString`, `$IGImportFormats`; support for importing Graph6, Digraph6 and Sparse6.
 - Export functions: `IGExport`, `IGExportString`, `$IGExportFormats`; support for exporting standards-compliant GraphML that can be read by other igraph interfaces (R, Python).
 - Matrix functions: `IGZeroDiagonal`, `IGTakeUpper`, `IGTakeLower`, `IGAdjacencyMatrixPlot`, `IGKirchhoffMatrix`, `IGJointDegreeMatrix`.
 - Bipartite graphs: `IGBipartiteIncidenceMatrix`, `IGBipartiteIncidenceGraph`, `IGBipartiteProjections`.
 - Random walks:  `IGRandomEdgeWalk`, `IGRandomEdgeIndexWalk`
 - Connectivity: `IGGiantComponent`, `IGMinimalSeparators`, `IGFindMinimumCuts`, `IGFindMinimalCuts`, `IGBridges`, `IGConnectedComponentSizes`, `IGWeaklyConnectedComponentSizes`, `IGGomoryHuTree`.
 - Efficiency measures: `IGGlobalEfficiency`, `IGLocalEfficiency`, `IGAverageLocalEfficiency`.
 - Neighbour degrees: `IGAverageNeighborDegree`, `IGAverageDegreeConnectivity`
 - Other testing functions: `IGTriangleFreeQ`, `IGSelfComplementaryQ`,  `IGCactusQ`,`IGNullGraphQ`, `IGHomeomorphicQ`, `IGAdjacentVerticesQ`.
 - Other new functions: `IGNeighborhoodSize`, `IGCoreness`, `IGVoronoiCells`, `IGSmoothen`, `IGSimpleGraph`, `IGPartitionsToMembership`, `IGMembershipToPartitions`, `IGSinkVertexList`, `IGSourceVertexList`, `IGIsolatedVertexList`, `IGReorderVertices`, `IGTakeSubgraph`, `IGDisjointUnion`, `IGTryUntil`, `IGVertexAssociate`.
 - Added `IGIndexEdgeList` for retrieving the edge list of a graph in terms of vertex indices. This function is very fast and returns a packed array. It facilitates the efficient implementation of graph processing functions in pure _Mathematica_ code, or interfacing with C libraries.

Updates to existing functions:

 - Community detection: Several functions support the `"ClusterCount"` option now.
 - `IGRewireEdges` now supports rewiring only the start or endpoint of directed edges (instead of both).
 - `IGBipartiteQ` now supports checking that a given partitioning is valid for a bipartite graph.
 - `IGBipartitePartitions` now provides control over the ordering of partitions using its second argument.
 - Isomorphism functions now ignore the directedness of empty graphs.
 - Isomorphism functions can now take vertex or edge colours from graph attributes.
 - `IGBlissCanonicalGraph` will now include include vertex colours into its output when appropriate, encoded as the `"Color"` vertex property.
 - `IGDistanceCounts` now optionally takes a list of starting vertices.
 - `IGBetweenness(Estimate)` and `IGCloseness(Estimate)` now optionally take a set of vertices to do the calculation on.
 - `IGBetweenness(Estimate)` and `IGEdgeBetweenness(Estimate)` now have the `Normalized` option.
 - `IGRandomWalk` now supports edge weights and the `EdgeWeights` option. Use `EdgeWeights -> None` to ignore weights, and thus restore the previous behaviour.
 - Clustering coefficient functions now support the `"ExcludeIsolates"` option.

Incompatible changes from IGraph/M 0.3:

 - A flat namespace structure is used. Functions from ``IGraphM`Utilities` `` have been moved to ``IGraphM` ``.
 - Renamed `"MultipleEdges"` option to `MultiEdges` for convenient typing and auto-completion.
 - Renamed `IGMinSeparators` to `IGMinimumSeparators`.
 - Renamed `IGMakeLattice` to `IGSquareLattice`. The name `IGMakeLattice` works, but it is deprecated.

Other changes:

 - Improved compatibility with _Mathematica_ versions 11.2–12.1; handling of `TwoWayRule` as an edge specification and support for edge tagged graphs.
 - Bug fixes, performance improvements, documentation updates, and general polish.

Notes:

 - IGraph/M 0.4 will be the last version of IGraph/M to support _Mathematica_ 10.0. Starting with IGraph/M 0.5, it will only work [with versions of _Mathematica_ that are still supported by Wolfram Research](http://support.wolfram.com/kb/22107).

#### v0.3.0

 - Compatibility with *Mathematica* 11.1.
 - Additional random graph generators (`IGBarabasiAlbertGame`, `IGStaticFitnessGame`, `IGStaticPowerLawGame`, `IGGeometricGame`, `IGGrowingGame`, `IGForestFireGame`).
 - Bipartite graph layout, `IGLayoutBipartite`.
 - Performance improvements in `IGMaximalCliqueSizeCounts`.
 - `IGEccentricity` and `IGRadius`.
 - `IGVertexContract`, to contract multiple vertex sets simultaneously.
 - `IGRandomWalk`, to take a random walk on a graph.
 - ``IGraphM`Utilities` `` package.
 - Significantly improved package loading performance.
 - Bug fixes and polish.

#### v0.2.2

This is a bugfix release with several minor fixes.  Important changes to be aware of:

 - VF2 isomorphism functions now work when *both* vertex and edge colours are specified at the same time.
 - `IGBlissCanonicalGraph` is now suitable for filtering isomorphic duplicates. Its output may be different from previous releases.

#### v0.2.1

This is a bugfix release.  The following changes require special mention:

 - `IGFeedbackArcSet` could return wrong results for some graphs. This is now fixed.
 - The `"RemovedEdges"` property returned by `IGCommunitiesEdgeBetweenness` could be incorrect for some graphs.  This is now fixed.  Other properties were not affected.
 - The Hierarchical Clustering package is no longer auto-loaded by IGraph/M to avoid shadowing `DistanceMatrix`, a new built-in function added in Mathematica 10.3.  Load this package manually (``<<HierarchicalClustering` ``) to work with the `"HierarchicalClusters"` property of `IGClusterData` objects.

A number of other small bugs were also fixed.

#### v0.2.0

- Significant performance improvements for many functions
- New and extended functions for shortest path calculations (extended `IGDistanceMatrix`, `IGDistanceCounts`, `IGDistanceHistogram`, `IGDiameter`, `IGFindDiameter`)
- Support for weighted clique calculations (`IGWeightedCliques`, `IGMaximalWeightedCliques`, `IGLargestWeightedCliques`, `IGWeightedCliqueNumber`)
- Additional new functions
- Syntax highlighting for functions
- Raspberry Pi support
- Compatibility with _Mathematica_ 10.4 on OS X.
- Bug fixes
