# IGraph/M – igraph for Mathematica

## What is IGraph/M?

IGraph/M is a *Mathematica* interface to the [igraph](http://igraph.org/) graph manipulation and analysis library.  It works directly with *Mathematica*'s builtin `Graph` datatype and focuses on functionality that is either not already present in *Mathematica* (such as motifs, minimum feedback arc set, random rewiring, etc.), or for which igraph uses different implementations (isomorphism, sampling random graphs with given degree sequence, clique counts, etc.)

## What IGraph/M is not

IGraph/M is *not a replacement* for Mathematica's graph manipulation functionality.  Instead it is meant to complement it.  Thus it works directly with standard `Graph` objects instead of introducing its own graph type.  Functions for trivial tasks such as adding or removing vertices or edges, returning the vertex or edge count, creating standard graphs like cycle graphs, complete graphs, etc. are not provided.

## Why create IGraph/M?

igraph is one of the broadest open source graph manipulation packages available.  Many of its functions are of use to Mathematica users, either because equivalents don't already exist in Mathematica, or because they can be used to verify Mathematica's own results.  My previous package, [IGraphR][2], already provides relatively convenient access to igraph's R interface from Mathematica, but unfortunately its underlying technology, [RLink](http://reference.wolfram.com/language/RLink/guide/RLink.html), suffers from reliability and performance problems.  IGraph/M uses igraph's C interface through [LibraryLink](http://reference.wolfram.com/language/LibraryLink/tutorial/Overview.html), which makes it much faster, more robust and easier to use with Mathematica's parallel tools.  My main motivation for starting IGraph/M was to get better performance and reliable parallelization.

### Functionality highlights

 - Centrality measures for weighted graphs
 - Estimates of vertex betweenness, edge betweenness and closeness centrality; for large graphs
 - Community detection algorithms
 - Minimum feedback arc set for weighted and unweighted graphs
 - Find all cliques (not just maximal ones), count cliques without listing them
 - Count 3- and 4-motifs, list triangles
 - Rewire edges, keeping either the density or the degree sequence
 - Alternative algorithms for isomorphism testing: BLISS, VF2, LAD
 - Subgraph isomorphism (including induced subgraphs with LAD)
 - Isomorphism for coloured graphs
 - Test if a degree sequence is graphical
 - Alternative algorithms for generating random graphs with given degree sequence
 - Additional layout algorithms, layouts for weighted graphs
 - Layout algorithms can continue from existing vertex coordinates
 - Biconnected components, articulation points, find all minimum vertex cuts
 - Several other specialized functions not mentioned here ...

## Installation

IGraph/M can be installed like any other Mathematica application.

 - [Download the zip archive from GitHub](https://github.com/szhorvat/IGraphM/releases)
 - Open Mathematica's "Applications" directory by evaluating `SystemOpen@FileNameJoin[{$UserBaseDirectory, "Applications"}]`
 - Unzip the archive, find the `IGraphM` directory, and move it to Mathematica's "Applications" directory.  If earlier versions of the package were installed, they must be fully removed first.

 IGraph/M requires Mathematica 10 or later.  Binaries are included for Windows 64 bit, OS X 10.9 or later, and Linux x86_64.  For other operating systems the package must be compiled from source (see [Development.md](Development.md) for guidance).

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

The documentation is not yet ready and contributions are most welcome.  If you would like to help out with the documentation, please see below.

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

IGraph/M is currently under development, and a few bugs are to be expected.  However, I try not to release a new version until most problems I know of are fixed.  If you do find a problem, please [open an issue on GitHub](https://github.com/szhorvat/IGraphM/issues).

During the version 0.1.x series IGraph/M should be considered unstable.  Function names, options, and the structure of return values may change without notice.  Thus I do encourage you to use it interactively, but at this point I do not recommend basing other packages on IGraph/M.  Starting from version 0.2.0 I will try to avoid making any breaking changes.

### Known issues and workarounds

 * `Graph[Graph[...], ...]` returned

   Sometimes layout functions may return an expression which looks like

        Graph[ Graph[...], VertexCoordinates -> {...} ]

   or similar. A property does not get correctly applied to the graph.  This is due to a bug in Mathematica. I believe I have worked around most of these issues, but if you encounter them, one possible workaround is to cycle the graph `g` through some other representation, e.g. `g = Uncompress@Compress[g]`.

 * Graphlet decomposition functions may crash the kernel. This is due to an igraph bug.

## License

The IGraph/M source code is released under the MIT license, see LICENSE.txt.

igraph (and consequently the IGraph/M binary packages) can be distributed under the terms of the [GPLv2](http://opensource.org/licenses/GPL-2.0).

 [1]: https://bitbucket.org/szhorvat/ltemplate
 [2]: http://szhorvat.net/pelican/using-igraph-from-mathematica.html
