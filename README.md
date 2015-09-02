#IGraph/M – igraph for Mathematica

##What is IGraph/M?

IGraph/M is a *Mathematica* interface to the [igraph](http://igraph.org/) graph manipulation and analysis library.  It works directly with *Mathematica*'s builtin `Graph` datatype and focuses on functionality that is either not already present in *Mathematica* (such as motifs, minimum feedback arc set, random rewiring, etc.), or for which igraph uses different implementations (isomorphism, sampling random graphs with given degree sequence, clique counts, etc.)

##Why make IGraph/M?

igraph is one of the broadest graph manipulation packages available.  Many of its functions are of use to Mathematica users, either because equivalents don't already exist in Mathematica, or because they can be used to verify Mathematica's own results.  My previous package, [IGraphR][2], already provides relatively convenient access to all igraph's R interface from Mathematica, but unfortunately its underlying technology, [RLink](http://reference.wolfram.com/language/RLink/guide/RLink.html), suffers from reliability and performance problems.  IGraph/M uses igraph's C interface through [LibraryLink](http://reference.wolfram.com/language/LibraryLink/tutorial/Overview.html), which makes it much faster, more robust and easier to use with Mathematica's parallel tools.

##Installation

Presently IGraph/M does not come with precompiled binary components, so you need to compile it yourself.  It also depends on the [LTemplate package][1], which needs to be installed (in the future it will embed LTemplate).  The installation steps are as follows:

 1. Place the `IGraphM` directory in `$UserBaseDirectory/Applications`
 2. Install `LTemplate` by placing it in the same location, `$UserBaseDirectory/Applications`
 3. Make sure the igraph headers and library are installed on your system.
 4. Edit `BuildSettings.m` within the `IGraphM` directory and set the location of the igraph library on your system.  Set up the compiler and any special compiler options to use.  Note that IGraph/M requires a C++ compiler with C++11 support.  See the [CCompilerDriver User Guide](http://reference.wolfram.com/language/CCompilerDriver/tutorial/Overview.html) for details.
 5. Start Mathematica, load the package using ``<<IGraphM` ``, and evaluate `RecompileIGraphM[]`.  If there were no errors, IGraph/M is ready to use.  Check ``?IGraphM`*`` to see what functions are available.
 
Currently IGraph/M is still incomplete and under development.  The documentation is not yet ready, but the usage messages should be descriptive.  For details on how the different functions work, refer to [the C/igraph documentation](http://igraph.org/c/doc/).

##Contributions

Contributions to IGraph/M are most welcome.  igraph is a large graph library with diverse functionality.  I primarily focused on providing an interface to functions that I need myself, and I do not have time to cover all igraph functions.  If you are interested in extending IGraph/M, send me an email to get technical guidance.  IGraph/M is uses the [LTemplate package][1], and acts as a driver for LTemplate development.  I recommend starting by reading the LTemplate tutorial.

Desired but not yet completed functionality:

 - hierarchical random graphs
 - spectral coarse graining
 - layout algorithms
 - community detection
 - graphlets

Remember, if you need to use any of these from *Mathematica* today, there is always [IGraphR][2].

##License

The IGraph/M source code is released under the MIT license, see LICENSE.txt.

 [1]: https://bitbucket.org/szhorvat/ltemplate
 [2]: http://szhorvat.net/pelican/using-igraph-from-mathematica.html
