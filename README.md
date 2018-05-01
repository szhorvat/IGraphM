[![GitHub (pre-)release](https://img.shields.io/github/release/szhorvat/IGraphM/all.svg)](https://github.com/szhorvat/IGraphM/releases)
[![Join the chat at https://gitter.im/IGraphM/Lobby](https://badges.gitter.im/IGraphM/Lobby.svg)](https://gitter.im/IGraphM/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

# [IGraph/M – igraph for Mathematica][main]

The IGraph/M package provides a [_Mathematica_](http://www.wolfram.com/mathematica/) interface to the popular [igraph](http://igraph.org/) network analysis and graph theory package, as well as many other functions for working with graphs in _Mathematica_.  Check out the [blog post][main] for an overview.

## Installation

First, note that the system requirements are _Mathematica_ 10.0 or later, 64-bit Windows/macOS/Linux, or Raspberry Pi.

The simplest way to install automatically from the internet is to evaluate the following commands in Mathematica:

```mathematica
Needs["PacletManager`"];
PacletInstall["IGraphM", "Site" -> "http://raw.githubusercontent.com/paclets/PacletServer/master"]
```

To update from an older version, use:

```mathematica
Needs["PacletManager`"];
PacletUpdate["IGraphM", "Site" -> "http://raw.githubusercontent.com/paclets/PacletServer/master"]
```

IGraph/M can also be installed manually in the same way as any _Mathematica_ application distributed as a paclet.

Download the `.paclet` file from [the GitHub releases page](https://github.com/szhorvat/MaTeX/releases), and [install it using the `PacletInstall` function in Mathematica](http://mathematica.stackexchange.com/q/141887/12).  For example, assuming that the file `IGraphM-0.3.99.1.paclet` was downloaded into the directory `~/Downloads`, evaluate

```mathematica
Needs["PacletManager`"]
PacletInstall["~/Downloads/IGraphM-0.3.99.1.paclet"]
```

IGraph/M requires Mathematica 10.0.2 or later.  Binaries are included for Windows 64-bit, OS X 10.9 or later, Linux x86_64 and Raspbian (Linux ARM on Raspberry Pi).  For other operating systems the package must be compiled from source (see [Development.md](Development.md) for guidance).

After installation, the package can now be loaded with

    << IGraphM`

Check that it works by evaluating `IGVersion[]`, then continue to the documentation with `IGDocumentation[]`.

## Documentation

To open the documentation notebook, evaluate

    Needs["IGraphM`"]
    IGDocumentation[]

or search for "igraphm" in Mathematica's Documentation Center.

The documentation is not yet complete and contributions are very welcome.  If you would like to help out with expanding the documentation, send me an email.

For additional details about functions, or for paper references for the methods used, also check [the igraph documentation pages](http://igraph.org/c/doc/).

## Contributions

Contributions to IGraph/M are very welcome!  igraph is a large graph library with diverse functionality.  I primarily focused on providing an interface to functions that I need myself, and I do not have time to cover all of igraph.  However, the main framework is there, and adding new functions is relatively quick and easy.

If you are interested in extending IGraph/M, send me an email to get technical guidance.  IGraph/M uses the [LTemplate package][1] to simplify writing LibraryLink code, and acts as a driver for LTemplate development.  I recommend starting by reading the LTemplate tutorial.

Please see [Development.md](Development.md) for additional information.

If you do not program in C++, there are many other ways you can still help out:

 - Implement functions in pure Mathematica code
 - Expand the (currently incomplete) documentation
 - Write unit tests
 - Just test the package and try to find problems

You can contact me by email.  Evaluate the following in Mathematica to get my email address:

    Uncompress["1:eJxTTMoPChZiYGAorsrILypLLHFIz03MzNFLzs8FAG/xCOI="]

## Bugs and troubleshooting

IGraph/M is currently under development, and a few bugs are to be expected.  However, I try not to release a new version until most problems I know of are fixed.  If you do find a problem, please [open an issue on GitHub](https://github.com/szhorvat/IGraphM/issues) or write [in the chatroom](https://gitter.im/IGraphM/Lobby). **In the bug report include the output of ``IGraphM`Developer`GetInfo[]``.**

### Troubleshooting

  * **"Cannot open ``IGraphM` ``"**

    This message will be shown in the following situations:

    - IGraph/M is not installed. Please follow the installation instructions above carefully.

    - IGraph/M is not compatible with your system. Please review the requirements in the Installation section above.  Additional symptoms will be that `PacletFind["IGraphM"]` returns `{}` but `PacletFind["IGraphM", "MathematicaVersion" -> All, "SystemID" -> All]` returns a non-empty list.

  * **"Cannot open ``LTemplate`LTemplatePrivate` ``"**

    Loading fails with this error if one tries to use the `master` branch (development branch) without the necessary dependencies.

    Please download IGraph/M from the [GitHub releases page instead](https://github.com/szhorvat/IGraphM/releases).  It includes everything needed to run the package.

    Do not clone the git repository and do not use the `master` branch unless you want to develop IGraph/M.

  * **"No build settings found. Please check `BuildSettings.m`"**

    This error may be shown when trying to load IGraph/M on an incompatible platform.  Currently, the following platforms are supported: Windows 64-bit, OS X 10.9 or later, Linux x86_64 with gcc 5 or later, Raspberry Pi with Raspbian Stretch.

    It may be possible to run IGraph/M on other platforms, however, it will be necessary to compile it from source.

    If you see this error on a supported platform, please file a bug report and include the output of ``IGraphM`Developer`GetInfo[]``.


### Known issues and workarounds

   * `Graph[Graph[...], ...]` is returned

     Sometimes layout functions may return an expression which looks like

         Graph[ Graph[...], VertexCoordinates -> {...} ]

     or similar. A property does not get correctly applied to the graph.  This is due to a bug in Mathematica. I believe I have worked around most of these issues, but if you encounter them, one possible workaround is to cycle the graph `g` through some other representation, e.g. `g = Uncompress@Compress[g]`.

   * Graphlet decomposition functions may crash the kernel. This is due to [a bug in the igraph C core](https://github.com/igraph/igraph/issues/869).

   * The graphs returned by `IGBipartiteGameGNM` and `IGBipartiteGameGNP` may not render when using the `DirectedEdges -> True` and `"Bidirectional" -> True` options.  This is due to a bug in Mathematica's  `"BipartiteEmbdding"` graph layout and can be corrected by passing `GraphLayout -> Automatic` to these functions.

   * LAD functions may crash with certain inputs.  [This is a bug in the igraph C core.](https://github.com/szhorvat/IGraphM/issues/1)

   * `IGDocumentation[]` may not work with Mathematica 11.1 on Linux.  Please enter `IGraphM/IGDocumentation` in the Documentation Center address bar to access it (or simply search for `igraph` in the Documentation Center).

   * See also https://github.com/szhorvat/IGraphM/issues

## Revision history

##### v0.4.0dev (work in progress)

 - New deterministic graph generators: `IGKautzGraph`, `IGCompleteGraph`, `IGCompleteAcyclicGraph`, `IGDeBruijnGraph`, `IGChordalRing`, `IGEmptyGraph`, `IGRealizeDegreeSequence`, `IGFromPrufer`, `IGKaryTree`, `IGSymmetricTree`, `IGBetheLattice`, `IGTriangularLattice`.
 - New random graph generators: `IGWattsStrogatzGame`, `IGCallawayTraitsGame`, `IGEstablishmentGame`, `IGTreeGame`, `IGErdosRenyiGameGNM`, `IGErdosRenyiGameGNP`.
 - New weighted graph functions: `IGWeightedSimpleGraph`, `IGWeightedUndirectedGraph`, `IGWeightedVertexDelete`, `IGWeightedSubgraph`.
 - New graph colouring functions: `IGVertexColoring`, `IGEdgeColoring`, `IGKVertexColoring`, `IGKEdgeColoring`, `IGMinimumVertexColoring`, `IGMinimumEdgeColoring`, `IGChromaticNumber`, `IGChromaticIndex`.
 - New functions for mesh/graph conversion: `IGMeshGraph`, `IGMeshCellAdjacencyMatrix`, `IGMeshCellAdjacencyGraph`.
 - New functions for lattice generation: `IGLatticeMesh`, `IGTriangularLattice`.
 - New functions for centralization: `IGDegreeCentralization`, `IGBetweennessCentralization`, `IGClosenessCentralization`, `IGEigenvectorCentralization`.
 - Community detection: several functions support the `"ClusterCount"` option now; added `IGCommunitiesFluid`.
 - Added `IGIndexEdgeList` for retrieving the edge list of a graph in terms of vertex indices. This function is very fast and returns a packed array. It facilitates the efficient implementation of graph processing functions in pure Mathematica code, or interfacing with C libraries.
 - Other new functions: `IGVertexTransitiveQ`, `IGEdgeTransitiveQ`, `IGSymmetricQ`, `IGTriangleFreeQ`, `IGSelfComplementaryQ`, `IGSpanningTree`, `IGRandomSpanningTree`, `IGSpanningTreeCount`, `IGUnfoldTree`, `IGTreeQ`, `IGTreelikeComponents`, `IGRandomEdgeWalk`, `IGRandomEdgeIndexWalk`, `IGBipartiteIncidenceMatrix`, `IGBipartiteIncidenceGraph`, `IGBipartiteProjections`, `IGNeighborhoodSize`, `IGCoreness`, `IGVoronoiCells`, `IGMinimalSeparators`, `IGBridges`, `IGConnectedComponentSizes`, `IGWeaklyConnectedComponentSizes`, `IGJointDegreeMatrix`, `IGMycielskian`, `IGSmoothen`, `IGHomeomorphicQ`.
 - Updates:
    * `IGRewireEdges` now supports rewiring only the start or endpoint of directed edges (instead of both).
    * `IGBipartiteQ` now supports checking that a given partitioning is valid for a bipartite graph.
    * `IGBipartitePartitions` now provides control over the ordering of partitions using its second argument.
    * Isomorphism functions now ignore the directedness of empty graphs.
    * `IGDistanceCounts` now optionally takes a list of starting vertices.
    * Isomorphism functions can now take vertex or edge colours from graph attributes.
    * `IGBetweenness(Estimate)` and `IGCloseness(Estimate)` now optionally take a set of vertices to do the calculation on.
    * `IGBetweenness(Estimate)` and `IGEdgeBetweenness(Estimate)` now have the `Normalized` option.
    * `IGBlissCanonicalGraph` will now include include vertex colours into its output when appropriate, encoded as the `"Color"` vertex property.
    * `IGRandomWalk` now supports edge weights and the `EdgeWeights` option. Use `EdgeWeights -> None` to ignore weights, and thus restore the previous behaviour.
 - Renamed `IGMinSeparators` to `IGMinimumSeparators`.
 - New utility functions:
    * Weighted graphs: `IGUnweighted`, `IGWeightedAdjacencyGraph`, `IGVertexWeightedQ`, `IGEdgeWeightedQ`, `IGVertexStrength`, `IGVertexInStrength`, `IGVertexOutStrength`.
    * Easier property handling and graph styling:  `IGVertexProp`, `IGEdgeProp`, `IGVertexMap`, `IGEdgeMap`, `IGVertexPropertyList`, `IGEdgePropertyList`.
    * Export functions: `IGExport`, `IGExportString`, `$IGExportFormats`; support for exporting standards-compliant GraphML that can be read by other igraph interfaces (R, Python).
    * Other: `IGNullGraphQ`, `IGSimpleGraph`, `IGShorthand`, `IGPartitionsToMembership`, `IGMembershipToPartitions`, `IGSinkVertexList`, `IGSourceVertexList`, `IGReorderVertices`, `IGOrientTree`, `IGTake`, `IGZeroDiagonal`, `IGAdjacencyMatrixPlot`, `IGKirchhoffMatrix`, `IGGiantComponent`, `IGDisjointUnion`.
 - Improved compatibility with Mathematica 11.2, handling of `TwoWayRule` as an edge specification.
 - Bug fixes, performance improvements, documentation updates, and polish.

##### v0.3.0

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

##### v0.2.2

This is a bugfix release with several minor fixes.  Important changes to be aware of:

 - VF2 isomorphism functions now work when *both* vertex and edge colours are specified at the same time.
 - `IGBlissCanonicalGraph` is now suitable for filtering isomorphic duplicates. Its output may be different form previous releases.

##### v0.2.1

This is a bugfix release.  The following changes require special mention:

 - `IGFeedbackArcSet` could return wrong results for some graphs. This is now fixed.
 - The `"RemovedEdges"` property returned by `IGCommunitiesEdgeBetweenness` could be incorrect for some graphs.  This is now fixed.  Other properties were not affected.
 - The Hierarchical Clustering package is no longer auto-loaded by IGraph/M to avoid shadowing `DistanceMatrix`, a new built-in function added in Mathematica 10.3.  Load this package manually (``<<HierarchicalClustering` ``) to work with the `"HierarchicalClusters"` property of `IGClusterData` objects.

A number of other small bugs were also fixed.

##### v0.2.0

- Significant performance improvements for many functions
- New and extended functions for shortest path calculations (extended `IGDistanceMatrix`, `IGDistanceCounts`, `IGDistanceHistogram`, `IGDiameter`, `IGFindDiameter`)
- Support for weighted clique calculations (`IGWeightedCliques`, `IGMaximalWeightedCliques`, `IGLargestWeightedCliques`, `IGWeightedCliqueNumber`)
- Additional new functions
- Syntax highlighting for functions
- Raspberry Pi support
- Compatibility with Mathematica 10.4 on OS X.
- Bug fixes

## License

The IGraph/M source code is released under the MIT license.

igraph (and consequently the IGraph/M binary packages) can be distributed under the terms of the [GPLv2](http://opensource.org/licenses/GPL-2.0).

See [LICENSE.md](IGraphM/LICENSE.md) for details.

 [1]: https://github.com/szhorvat/LTemplate/
 [main]: http://szhorvat.net/mathematica/IGraphM
 [chat]: https://gitter.im/IGraphM/Lobby
