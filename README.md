[![GitHub (pre-)release](https://img.shields.io/github/release/szhorvat/IGraphM/all.svg)](https://github.com/szhorvat/IGraphM/releases)
[![Join the chat at https://gitter.im/IGraphM/Lobby](https://badges.gitter.im/IGraphM/Lobby.svg)](https://gitter.im/IGraphM/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

# [IGraph/M – igraph for Mathematica](http://szhorvat.net/mathematica/IGraphM)

**For the impatient:**

 - Download the `.paclet` file from [the releases page](https://github.com/szhorvat/IGraphM/releases).
 - Install it using the `PacletInstall` function, as explained [here](http://mathematica.stackexchange.com/q/141887/12).
 - Evaluate ``<<IGraphM` `` to load the package.

Requirements: _Mathematica_ 10.0 or later, 64-bit Windows/macOS/Linux, or Raspberry Pi.

Do _not_ use the `master` branch of the GitHub repository. It contains only the source code of the package, and cannot be used without building it first.

## Introduction

#### What is IGraph/M?

IGraph/M is a *Mathematica* interface to the [igraph](http://igraph.org/) graph manipulation and analysis library.  It works directly with *Mathematica*'s builtin `Graph` expressions and focuses on functionality that is either not already present in *Mathematica* (such as motifs, minimum feedback arc set, random rewiring, coloured graph isomorphism, etc.), or for which igraph uses different implementations (isomorphism, sampling random graphs with given degree sequence, etc.)

#### What IGraph/M is not

IGraph/M is *not a replacement* for Mathematica's graph manipulation functionality.  Instead, it is meant to complement it.  Thus it works directly with standard `Graph` expressions instead of introducing its own graph type.  Functions for trivial tasks such as adding or removing vertices or edges, returning the vertex or edge count, creating standard graphs like cycle graphs, complete graphs, etc. are not provided.

#### Why create IGraph/M?

igraph is one of the most complete open source graph manipulation libraries available.  Many of its functions are of use to *Mathematica* users, either because equivalents don't already exist in *Mathematica*, or because they can be used to verify *Mathematica*'s own results.  IGraph/M also includes several utility functions that are not based on the igraph.

#### Functionality highlights

 - Interruption support: using Evaluate → Abort Evaluation in Mathematica works with most IGraph/M functions.
 - Centrality measures for weighted graphs.
 - Fast estimates of vertex betweenness, edge betweenness and closeness centrality; for large graphs.
 - Community detection algorithms.
 - Minimum feedback arc set for weighted and unweighted graphs.
 - Vertex and edge colouring.
 - Find all cliques (not just maximal ones); count cliques of different sizes without storing them; work with weighted cliques.
 - Count 3- and 4-motifs; find triangles.
 - Histogram for shortest path lengths (weighted and unweighted).
 - Rewire edges, keeping either the density or the degree sequence.
 - Alternative fast algorithms for isomorphism testing: Bliss, VF2, LAD.
 - Subgraph isomorphism (including induced subgraphs with LAD).
 - Isomorphism for edge or vertex coloured graphs; multigraph isomorphism implementation based on edge colouring.
 - Alternative algorithms for generating random graphs with given degree sequence; test for degree sequence graphicality.
 - Additional layout algorithms: most work with weighted graphs and can continue the layout optimization starting from a given set of vertex positions.
 - Biconnected components, articulation points, find all minimum vertex cuts.
 - Functions for easy graph property transformations and graph styling.
 - Functions for converting geometric meshes to graphs.
 - Several other specialized functions not mentioned here ...

The documentation contains many examples and can be accessed using the `IGDocumentation[]` command.

## Installation

IGraph/M can be installed like any other Mathematica application distributed as a paclet.

Download the `.paclet` file from [the GitHub releases page](https://github.com/szhorvat/MaTeX/releases), and [install it using the `PacletInstall` function in Mathematica](http://mathematica.stackexchange.com/q/141887/12).  For example, assuming that the file `IGraphM-0.3.94.paclet` was downloaded into the directory `~/Downloads`, evaluate

        Needs["PacletManager`"]
        PacletInstall["~/Downloads/IGraphM-0.3.94.paclet"]

IGraph/M requires Mathematica 10.0 or later.  Binaries are included for Windows 64-bit, OS X 10.9 or later, Linux x86_64 and Raspbian (Linux ARM on Raspberry Pi).  For other operating systems the package must be compiled from source (see [Development.md](Development.md) for guidance).

The package can now be loaded with

    << IGraphM`

Check the available functions with

    ?IGraphM`*

and test that it works by evaluating `IGVersion[]`.

## Documentation

Use

    <<IGraphM`
    IGDocumentation[]

to open the documentation notebook, or search for "igraphm" in Mathematica's Documentation Centre.

The documentation is not yet complete and contributions are very welcome.  If you would like to help out with the documentation, send me an email.

For additional details about functions, or for references, check also [the igraph documentation pages](http://igraph.org/c/doc/).

## Contributions

Contributions to IGraph/M are very welcome!  igraph is a large graph library with diverse functionality.  I primarily focused on providing an interface to functions that I need myself, and I do not have time to cover all of igraph.  However, the main framework is there, and adding new functions is relatively quick and easy.

If you are interested in extending IGraph/M, send me an email to get technical guidance.  IGraph/M uses the [LTemplate package][1] to simplify writing LibraryLink code, and acts as a driver for LTemplate development.  I recommend starting by reading the LTemplate tutorial.

Please see [Development.md](Development.md) for additional information.

#### Desired but not yet completed functionality

 - hierarchical random graphs
 - spectral coarse graining
 - maximum flows, minimum cuts
 - additional structural properties
 - update LAD isomorphism finder

Remember, if you need to use any of these from *Mathematica* today, there is always [IGraphR][2].

#### Other things you can help with

 - **compile Windows binaries, compile igraph for Windows**
 - write usage messages
 - write documentation
 - write unit tests

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

    This error may be shown when trying to load IGraph/M on an incompatible platform.  Currently, the following platforms are supported: Windows 64-bit, OS X 10.9 or later, Linux x86_64 with gcc 4.8 or later, Raspberry Pi with Raspbian Jessie.

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

 - New deterministic graph generators: `IGKautzGraph`, `IGKaryTree`, `IGCompleteGraph`, `IGCompleteAcyclicGraph`, `IGDeBruijnGraph`, `IGChordalRing`, `IGEmptyGraph`.
 - New random graph generators: `IGWattsStrogatzGame`, `IGCallawayTraitsGame`, `IGEstablishmentGame`.
 - Community detection: several functions support the `"ClusterCount"` option now; added `IGCommunitiesFluid`.
 - Other new functions: `IGVertexTransitiveQ`, `IGEdgeTransitiveQ`, `IGSymmetricQ`, `IGMinimalSeparators`, `IGSpanningTree`, `IGRandomEdgeWalk`, `IGRandomEdgeIndexWalk`, `IGVertexColoring`, `IGEdgeColoring`, `IGBipartiteIncidenceMatrix`, `IGBipartiteIncidenceGraph`, `IGMeshGraph`, `IGMeshCellAdjacencyMatrix`, `IGMeshCellAdjacencyGraph`, `IGNeighborhoodSize`, `IGBipartiteProjections`.
 - Updates:
    * `IGRewireEdges` now supports rewiring only the start or endpoint of directed edges (instead of both).
    * `IGBipartiteQ` now supports checking that a given partitioning is valid for a bipartite graph.
    * `IGBipartitePartitions` now provides control over the ordering of partitions using its second argument.
    * Isomorphism functions now ignore the directedness of empty graphs.
    * `IGDistanceCounts` now optionally takes a list of starting vertices.
    * Isomorphism functions can now take vertex or edge colours from graph attributes.
 - Renamed `IGMinSeparators` to `IGMinimumSeparators`.
 - New utility functions:
    * Weighted graphs: `IGUnweighted`, `IGWeightedAdjacencyGraph`, `IGVertexWeightedQ`, `IGEdgeWeightedQ`, `IGVertexStrength`, `IGVertexInStrength`, `IGVertexOutStrength`.
    * Easier property handling and graph styling:  `IGVertexProp`, `IGEdgeProp`, `IGVertexMap`, `IGEdgeMap`, `IGVertexPropertyList`, `IGEdgePropertyList`.
    * Export functions: `IGExport`, `IGExportString`, `$IGExportFormats`; support for exporting standards-compliant GraphML that can be read by other igraph interfaces (R, Python).
    * Other: `IGNullGraphQ`, `IGSimpleGraph`, `IGShorthand`.
 - Improved compatibility with Mathematica 11.2, handling of `TwoWayRule` as an edge specification.
 - Bug fixes, performance improvements, and polish.

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

The IGraph/M source code is released under the MIT license, see [LICENSE.txt](IGraphM/LICENSE.txt).

igraph (and consequently the IGraph/M binary packages) can be distributed under the terms of the [GPLv2](http://opensource.org/licenses/GPL-2.0).

 [1]: https://github.com/szhorvat/LTemplate/
 [2]: http://szhorvat.net/pelican/using-igraph-from-mathematica.html
 [chat]: https://gitter.im/IGraphM/Lobby
