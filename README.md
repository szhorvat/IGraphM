# IGraph/M – igraph for Mathematica

## What is IGraph/M?

IGraph/M is a *Mathematica* interface to the [igraph](http://igraph.org/) graph manipulation and analysis library.  It works directly with *Mathematica*'s builtin `Graph` datatype and focuses on functionality that is either not already present in *Mathematica* (such as motifs, minimum feedback arc set, random rewiring, coloured graph isomorphism, etc.), or for which igraph uses different implementations (isomorphism, sampling random graphs with given degree sequence, etc.)

## What IGraph/M is not

IGraph/M is *not a replacement* for Mathematica's graph manipulation functionality.  Instead it is meant to complement it.  Thus it works directly with standard `Graph` objects instead of introducing its own graph type.  Functions for trivial tasks such as adding or removing vertices or edges, returning the vertex or edge count, creating standard graphs like cycle graphs, complete graphs, etc. are not provided.

## Why create IGraph/M?

igraph is one of the most complete open source graph manipulation packages available.  Many of its functions are of use to Mathematica users, either because equivalents don't already exist in Mathematica, or because they can be used to verify Mathematica's own results.

### Functionality highlights

 - Interruption support: using Evaluate → Abort Evaluation in Mathematica works with most IGraph/M functions.
 - Centrality measures for weighted graphs.
 - Fast estimates of vertex betweenness, edge betweenness and closeness centrality; for large graphs.
 - Community detection algorithms.
 - Minimum feedback arc set for weighted and unweighted graphs.
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
 - Several other specialized functions not mentioned here ...

The documentation contains many examples and can be accessed using the `IGDocumentation[]` command.

## Installation

IGraph/M can be installed like any other Mathematica application.

 - [Download the zip archive from GitHub's releases page](https://github.com/szhorvat/IGraphM/releases). Use the "Source code (zip)" link.
 - Open Mathematica's "Applications" directory by evaluating `SystemOpen@FileNameJoin[{$UserBaseDirectory, "Applications"}]`.
 - Unzip the downloaded archive.  It will contain this `README.md` file, and a directory called `IGraphM`.  Move the `IGraphM` directory into Mathematica's "Applications" directory.  If earlier versions of the package were installed, they must be fully removed first.

 IGraph/M requires Mathematica 10.0 or later.  Binaries are included for Windows 64-bit, OS X 10.9 or later, Linux x86_64 and Raspbian (Linux ARM on Raspberry Pi).  For other operating systems the package must be compiled from source (see [Development.md](Development.md) for guidance).

The package can be loaded with

    << IGraphM`

Check the available functions with

    ?IGraphM`*

and test that it works by evaluating `IGVersion[]`.

## Documentation

Currently IGraph/M is still incomplete and under development.  Use

    <<IGraphM`
    IGDocumentation[]

to open the documentation notebook.

The documentation is not yet complete and contributions are most welcome.  If you would like to help out with the documentation, please see the Contributions section below.

For additional details about functions or for references, check also [the igraph documentation pages](http://igraph.org/c/doc/).

## Contributions

Contributions to IGraph/M are most welcome!  igraph is a large graph library with diverse functionality.  I primarily focused on providing an interface to functions that I need myself, and I do not have time to cover all of igraph.  However, the main framework is there, and adding new functions is relatively quick and easy.

If you are interested in extending IGraph/M, send me an email to get technical guidance.  IGraph/M uses the [LTemplate package][1] to simplify writing LibraryLink code, and acts as a driver for LTemplate development.  I recommend starting by reading the LTemplate tutorial.

Please see [Development.md](Development.md) for additional information.

#### Desired but not yet completed functionality

 - hierarchical random graphs
 - spectral coarse graining
 - maximum flows, minimum cuts
 - additional structural properties

Remember, if you need to use any of these from *Mathematica* today, there is always [IGraphR][2].

#### Other things you can help with

 - **compile Windows binaries, compile igraph for Windows**
 - write usage messages
 - write documentation
 - write unit tests

You can contact me in email.  Evaluate the following in Mathematica to get my email address:

    Uncompress["1:eJxTTMoPChZiYGAorsrILypLLHFIz03MzNFLzs8FAG/xCOI="]

## Bugs and troubleshooting

IGraph/M is currently under development, and a few bugs are to be expected.  However, I try not to release a new version until most problems I know of are fixed.  If you do find a problem, please [open an issue on GitHub](https://github.com/szhorvat/IGraphM/issues). **In the bug report include the output of ``IGraphM`Developer`GetInfo[]``.**

### Troubleshooting

  * **"Cannot open ``LTemplate`LTemplatePrivate` ``"**

    Loading fails with this error if one tries to use the `master` branch (development branch) without the necessary dependencies.

    Please download IGraph/M from the releases page instead: https://github.com/szhorvat/IGraphM/releases  It includes everything needed to run the package.  

    Do not clone the git repository and do not use the `master` branch unless you want to develop IGraph/M.

  * **"No build settings found. Please check `BuildSettings.m`"**

    This error is shown when trying to load IGraph/M on an incompatible platform.  Currently the following platforms are supported: Windows 64-bit, OS X 10.9 or later, Linux x86_64 with gcc 4.8 or later, Raspberry Pi with Raspbian Jessie.

    It may be possible to run IGraph/M on other platforms, however, it will be necessary to compile it form source.

    If you see this error on a supported platform, please file a bug report and include the output of ``IGraphM`Developer`GetInfo[]``.

  * **``<<IGraphM` `` hangs on Raspberry Pi**

    On a Raspberry Pi 1 computer, it takes about 30-40 seconds to load the package.  Please wait patiently.

### Known issues and workarounds

   * `Graph[Graph[...], ...]` is returned

     Sometimes layout functions may return an expression which looks like

        Graph[ Graph[...], VertexCoordinates -> {...} ]

     or similar. A property does not get correctly applied to the graph.  This is due to a bug in Mathematica. I believe I have worked around most of these issues, but if you encounter them, one possible workaround is to cycle the graph `g` through some other representation, e.g. `g = Uncompress@Compress[g]`.

   * Graphlet decomposition functions may crash the kernel. This is due to [a bug in the igraph C core](https://github.com/igraph/igraph/issues/869).

   * The graphs returned by `IGBipartiteGameGNM` and `IGBipartiteGameGNP` may not render when using the `DirectedEdges -> True` and `"Bidirectional" -> True` options.  This is due to a bug in Mathematica's  `"BipartiteEmbdding"` graph layout and can be corrected by passing `GraphLayout -> Automatic` to these functions.

   * LAD functions may crash with certain inputs.  [This is a bug in the igraph C core.](https://github.com/szhorvat/IGraphM/issues/1)

   * See also https://github.com/szhorvat/IGraphM/issues

## Revision history

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

 [1]: https://bitbucket.org/szhorvat/ltemplate
 [2]: http://szhorvat.net/pelican/using-igraph-from-mathematica.html
