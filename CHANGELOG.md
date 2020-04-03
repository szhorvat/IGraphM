## Revision history for [IGraph/M](README.md)

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
