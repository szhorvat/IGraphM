# IGraph/M – igraph for Mathematica

## What is IGraph/M?

IGraph/M is a *Mathematica* interface to the [igraph](http://igraph.org/) graph manipulation and analysis library.  It works directly with *Mathematica*'s builtin `Graph` datatype and focuses on functionality that is either not already present in *Mathematica* (such as motifs, minimum feedback arc set, random rewiring, etc.), or for which igraph uses different implementations (isomorphism, sampling random graphs with given degree sequence, clique counts, etc.)

## What IGraph/M is not

IGraph/M is *not a replacement* for Mathematica's graph manipulation functionality.  Instead it is meant to complement it.  Thus it works directly with standard `Graph` objects instead of introducing its own graph type.  Functions for trivial tasks such as adding or removing vertices or edges, returning the vertex or edge count, creating standard graphs like cycle graphs, complete graphs, etc. are not provided.

## Why create IGraph/M?

igraph is one of the broadest open source graph manipulation packages available.  Many of its functions are of use to Mathematica users, either because equivalents don't already exist in Mathematica, or because they can be used to verify Mathematica's own results.  My previous package, [IGraphR][2], already provides relatively convenient access to igraph's R interface from Mathematica, but unfortunately its underlying technology, [RLink](http://reference.wolfram.com/language/RLink/guide/RLink.html), suffers from reliability and performance problems.  IGraph/M uses igraph's C interface through [LibraryLink](http://reference.wolfram.com/language/LibraryLink/tutorial/Overview.html), which makes it much faster, more robust and easier to use with Mathematica's parallel tools.  My main motivation for starting IGraph/M was to get better performance and reliable parallelization.

### Functionality provided that is not present in Mathematica

 - Vertex betweenness centrality for weighted graphs
 - Estimates of vertex betweenness, edge betweenness and closeness centrality; for large graphs
 - Minimum feedback arc set for weighted and unweighted graphs
 - Find all cliques (not just maximal ones), count cliques without listing them
 - Count 3- and 4-motifs
 - Rewire edges, keeping either the density or the degree sequence
 - Alternative algorithms for isomorphism testing: BLISS, VF2, LAD
 - Subgraph isomorphism (including induced subgraphs with LAD)
 - Isomorphism for coloured graphs
 - Test if a degree sequence is graphical
 - Alternative algorithms for generating random graphs with given degree sequence
 - Layout algorithms that take weights into account

## Installation

IGraph/M can be installed like any other Mathematica application.

 - [Download the zip archive from GitHub](https://github.com/szhorvat/IGraphM/releases)
 - Open the Mathematica's "Applications" directory by evaluating `SystemOpen@FileNameJoin[{$UserBaseDirectory, "Applications"}]`
 - Unzip the archive, find the `IGraphM` directory, and move it to Mathematica's "Applications" directory.

**Currently IGraph/M comes only with binaries for OS X 10.9 or later and Linux.**  A Windows version is not available at the moment because I did not manage to compile igraph on Windows.  If you can help with creating a Windows binary, please see the Contributions section below.

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

Contributions to IGraph/M are most welcome!  igraph is a large graph library with diverse functionality.  I primarily focused on providing an interface to functions that I need myself, and I do not have time to cover all igraph functions.  However, the main framework is there, and adding new functions is relatively quick and easy.

If you are interested in extending IGraph/M, send me an email to get technical guidance.  IGraph/M uses the [LTemplate package][1] to simplify writing LibraryLink code, and acts as a driver for LTemplate development.  I recommend starting by reading the LTemplate tutorial.

####Compiling from source

To compile IGraph/M, you need will:

 - A C++ compiler with C++11 support.  I used GCC 4.9 and clang 3.6.
 - A recent development version of [igraph](https://github.com/igraph/igraph). You will need to compile it yourself. igraph 0.7.1 is not compatible.
 - The [LTemplate Mathematica package][1].  Please download and install it.
 - git for cloning the repository.

Then follow these steps:

 1. Clone the IGraphM repository and check out the master branch (do not use the release branch).  
 2. Edit `BuildSettings.m` and set the necessary paths to your igraph installation.  The available options are the same as for [CreateLibrary](http://reference.wolfram.com/language/CCompilerDriver/ref/CreateLibrary.html).
 3. Append the repository's root directory (i.e. the same directory where this readme file is found) to Mathematica's `$Path`.
 4. Load the package with ``<<IGraphM` ``.  It should automatically be compiled. It can also be recompiled using ``IGraphM`Developer`Recompile[]``.


#### Desired but not yet completed functionality

 - hierarchical random graphs
 - spectral coarse graining
 - layout algorithms that take edge weights into account
 - community detection
 - graphlets

Remember, if you need to use any of these from *Mathematica* today, there is always [IGraphR][2].

#### Other things you can help with

 - **compile Windows binaries, compile igraph for Windows**
 - write usage messages
 - write documentation
 - write unit tests

You can contact me in email.  Evaluate the following in Mathematica to get my email address:

    Uncompress["1:eJxTTMoPChZiYGAorsrILypLLHFIz03MzNFLzs8FAG/xCOI="]

## License

The IGraph/M source code is released under the MIT license, see LICENSE.txt.

igraph (and consequently the IGraph/M binary packages) can be distributed under the terms of the GPLv2.

 [1]: https://bitbucket.org/szhorvat/ltemplate
 [2]: http://szhorvat.net/pelican/using-igraph-from-mathematica.html
