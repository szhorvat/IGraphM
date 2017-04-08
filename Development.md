# Development guide for IGraph/M

Contributions to IGraph/M are most welcome!  igraph is a large graph library with diverse functionality.  I primarily focused on providing an interface to functions that I need myself, and I do not have time to cover all igraph functions.  However, the main framework is there, and adding new functions is relatively quick and easy.

If you are interested in extending IGraph/M, send me an email to get technical guidance.  IGraph/M uses the [LTemplate package][1] to simplify writing LibraryLink code, and acts as a driver for LTemplate development.  I recommend starting by reading the LTemplate tutorial.

## Compiling igraph

The following tools are required for building: a C++11-capable C++ compiler (gcc 4.8 works), `autoconf`, `automake`, `libtool`, `flex`, `bison`.

### OS X and Linux

In this guide, I will use `$HOME/local` as the installation directory for the various needed libraries.  You may want to use a different location for your system.

On OS X, set the following environment variable before compiling the libraries:

    export MACOSX_DEPLOYMENT_TARGET=10.9

##### GMP

[Download GMP](https://gmplib.org/), extract it to a directory with no spaces in its path (!) and compile with

    ./configure --prefix=$HOME/local --with-pic
    make
    make check

If the tests have passed, install it with `make install`.

To maximize compatibility with different types of CPUs, consider using the `--host=...` option to the `configure` script. Use the output of `configfsf.guess` as the target host.

##### igraph

Clone [this fork of the igraph](https://github.com/szhorvat/igraph) and check out the `IGraphM` branch. This fork is identical to the main igraph repository, except for a few small temporary patches that the latest version of IGraph/M may depend on. Compile as follows:

    export CPPFLAGS=-I$HOME/local/include LDFLAGS=-L$HOME/local/lib
    ./bootstrap.sh
    ./configure --prefix=$HOME/local --with-pic  --disable-graphml
    make
    make check

If the tests have passed, install it with `make install`.

### Windows

One option for compiling igraph on Windows is to use an MSYS shell to run the configure script.  Instructions for installing MSYS2 and the MinGW-w64 compiler are found at https://wiki.qt.io/MSYS2.  Install them in a directory with no spaces in its path.  The following instructions assume that libraries will be installed in `$HOME/local`.

##### GMP

Once the toolchain is set up, we can compile GMP.  [Download](https://gmplib.org/) and extract it.  Compile using

    ./configure --prefix=$HOME/local
    make
    make check

To maximize compatibility with different types of CPUs, consider using the `--host=...` option to the `configure` script. Use the output of `configfsf.guess` as the target host.

If the tests have passed, install it with `make install`

##### igraph

Clone [this fork of the igraph](https://github.com/szhorvat/igraph) and check out the `IGraphM` branch. This fork is identical to the main igraph repository, except for a few small temporary patches that the latest version of IGraph/M may depend on.  Compile and install as follows:

    export CPPFLAGS="-I$HOME/local/include -DMSDOS" LDFLAGS=-L$HOME/local/lib
    ./bootstrap.sh
    ./configure --prefix=$HOME/local --disable-graphml
    make
    make install

This will produce a DLL named `libigraph-0.dll` in `$HOME/local/bin`.  It must be copied into `IGraphM/LibraryResources/Windows-x86-64`.  When using the above version of MinGW-w64, it is also necessary to copy the dependencies `libgcc_s_seh-1.dll`, `libstdc++-6.dll` and `libwinpthread-1.dll` to the same directory.

IGraph/M needs to be told about what dependencies it has to load.  This is done by creating a file named `dependencies.m` in the same directory and adding `LibraryLoad` calls to it in the appropriate order.  For an example see `dependencies.m` on the `release` branch of the IGraph/M GitHub repo.

## Compiling IGraph/M

To compile IGraph/M, you will need:

 - A C++ compiler with C++11 support.  I used GCC 4.8.5 on Linux and clang 3.7 on OS X.
 - The [LTemplate Mathematica package][1].  Please download and install it.
 - git for cloning the repository.

Then follow these steps:

 1. Clone the IGraphM repository and check out the master branch (do not use the release branch).
 2. Edit `BuildSettings.m` and set the path to your igraph installation, where necessary.  The available options are the same as for [CreateLibrary](http://reference.wolfram.com/language/CCompilerDriver/ref/CreateLibrary.html).
 3. Append the repository's root directory (i.e. the same directory where this readme file is found) to Mathematica's `$Path`.
 4. Load the package with ``<<IGraphM` ``.  It should automatically be compiled. It can also be recompiled using ``IGraphM`Developer`Recompile[]``.

## Adding new functions to IGraph/M

If you decide to join IGraph/M development, please contact me for guidance.  The basics are covered below.

IGraph/M uses the LTemplate package to interface with C++. First please go through the [LTemplate tutorial][1].

#### Simple functions

To see the basic structure of a function, let us look at how a simple one such as `IGConnectedQ` is implemented.  On the C++ side we add a new member function to the `IG` class:

    bool connectedQ(bool strong) const {
        igraph_bool_t res;
        igCheck(igraph_is_connected(&graph, &res, strong ? IGRAPH_STRONG : IGRAPH_WEAK));
        return res;
    }

The `igCheck` function must wrap any igraph function that returns error codes.  It ensures that the error is passed on to Mathematica and reported appropriately.

On the Mathematica side we add the function to the library template as

    LFun["connectedQ", {True|False (* strongly connected *)}, True|False]

And finally we write a Mathematica function that calls this member function on an `IG` object.

    IGConnectedQ[g_?igGraphQ] :=
        catch@Block[{ig = igMake[g]},
            check@ig@"connectedQ"[True]
        ]

The pattern `g_?igGraphQ` will only match `Graph` objects compatible with igraph.  `igMake` will convert a Mathematica graph to an `IG` object.  Calls to `IG` member functions must be wrapped in `check`, which will catch any library errors and `Throw` them as exceptions.  `catch` is always used at the outermost level to catch and handle these exceptions.

`igMake` is able to handle and convert to igraph format any kind of Mathematica graph for which `igGraphQ` returns `True`.  It also transfers edge weights to igraph.  However, it is not very fast.  Two special alternatives are provided:

 * `igMakeFast` is very fast, but it only works with unweighted graphs and it does not preserve edge ordering.  We can use it in `IGConnectedQ` because this function does not make use of weights and its result does not depend on the ordering of edges.

 * `igMakeFastWeighted` is a bit slower, it can handle edge weights, but it also does not preserve edge ordering.  Use it with functions such as `IGBetweenness` which make use of weights, but their result is unaffected by edge ordering.

  * `igMake` handles weights and preserves edge ordering, but it is the slowest.  Use it with functions such as `IGEdgeBetweeness` where the result must be ordered accordingly with `EdgeList`, or which are affected by edge ordering in other ways.

`igMakeFast` could be used for `IGConnectedQ` to improve performance.

#### Working with edge weights

`igMake` will automatically transfer edge weights to the `IG` object.  To pass these weights to an igraph function as an `igraph_vector_t *`, use `passWeights()`.  See the `betweenness` member function for a simple example.  `passWeights` returns `NULL` for unweighted graphs.

#### Creating graphs

Let us look at a simple function that creates graphs, such as `IGGraphAtlas`.

    void graphAtlas(mint n) {
        destroy();
        igConstructorCheck(igraph_atlas(&graph, n));
    }

On the C++ side, the graph creator function must first `destroy()` the igraph graph that is currently stored in the `IG` object.  `IG` objects *always* store some graph, and are initialized with an empty one.

Now we are ready to call the igraph function that creates the graph (`igraph_atlas`).

Instead of `igCheck` we must use `igConstructorCheck` to handle igraph error codes.  This function ensures that the graph is not left in an inconsistent state.

On the Mathematica side, we write

    IGGraphAtlas[n_?Internal`NonNegativeMachineIntegerQ, opt : OptionsPattern[Graph]] :=
        catch@Block[{ig = Make["IG"]},
          check@ig@"graphAtlas"[n];
          applyGraphOpt[opt]@igToGraph[ig]
        ]

Since the `IG` object is not created from an existing graph, we use `Make["IG"]` to make an empty one. `igToGraph` converts an `IG` object back to a Mathematica graph.

Most IGraph/M functions that create graphs take all standard `Graph` options such as `GraphLayout`, `VertexStyle`, etc.  These can be applied using `applyGraphOpt`, as in the example above.



  [1]: https://github.com/szhorvat/LTemplate/
