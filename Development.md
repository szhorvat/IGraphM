# Development guide for IGraph/M

Contributions to IGraph/M are most welcome!  igraph is a large graph library with diverse functionality.  I primarily focused on providing an interface to functions that I need myself, and I do not have time to cover all igraph functions.  However, the main framework is there, and adding new functions is relatively quick and easy.

If you are interested in extending IGraph/M, send me an email to get technical guidance.  IGraph/M uses the [LTemplate package][1] to simplify writing LibraryLink code, and acts as a driver for LTemplate development.  I recommend starting by reading the LTemplate tutorial.

Help is also welcome with writing or editing documentation, adding examples to the documentation, and writing unit tests.


## Compiling igraph

The following tools are required for building the igraph C library: a C++ compiler (gcc 4.8 works), `autoconf`, `automake`, `libtool`, `flex`, `bison`.

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

##### GLPK

igraph already includes GLPK, but an external GLPK can be used for improved performance in `IGCommunitiesOptimalModularity` and `IGFeedbackArcSet`.

If desired, [download GLPK](https://www.gnu.org/software/glpk/), and compile it the same way as GMP.

    ./configure --prefix=$HOME/local --with-pic
    make
    make check

If the tests have passed, install it with `make install`.

When compiling igraph, pass `--with-external-glpk` to the `configure` script.

##### igraph

Clone [this fork of igraph](https://github.com/szhorvat/igraph) and check out the `IGraphM-040` branch. This fork is identical to the main igraph repository, except for a few small temporary patches that the latest version of IGraph/M may depend on. Compile as follows:

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

If the tests have passed, install it with `make install`.

##### GLPK

igraph already includes GLPK, but an external GLPK can be used for improved performance in `IGCommunitiesOptimalModularity` and `IGFeedbackArcSet`.

If desired, [download GLPK](https://www.gnu.org/software/glpk/), and compile it the same way as GMP.

    ./configure --prefix=$HOME/local --disable-reentrant
    make
    make check

If the tests have passed, install it with `make install`.

When compiling igraph, pass `--with-external-glpk` to the `configure` script.

##### igraph

Clone [this fork of igraph](https://github.com/szhorvat/igraph) and check out the `IGraphM` branch. This fork is identical to the main igraph repository, except for a few small temporary patches that the latest version of IGraph/M may depend on.  Compile and install as follows:

    export CPPFLAGS="-I$HOME/local/include -DMSDOS" LDFLAGS=-L$HOME/local/lib
    ./bootstrap.sh
    ./configure --prefix=$HOME/local --disable-graphml
    make
    make install

This will produce a DLL named `libigraph-0.dll` in `$HOME/local/bin`.  It must be copied into `IGraphM/LibraryResources/Windows-x86-64`.  When using the above version of MinGW-w64, it is also necessary to copy the dependencies `libgcc_s_seh-1.dll`, `libstdc++-6.dll` and `libwinpthread-1.dll` to the same directory.

IGraph/M needs to be told about what dependencies it has to load.  This is done by creating a file named `dependencies.m` in the same directory and adding `LibraryLoad` calls to it in the appropriate order.  For an example see `dependencies.m` on the `release` branch of the IGraph/M GitHub repo.


## Compiling LEMON

Follow the installation guide of LEMON: http://lemon.cs.elte.hu/trac/lemon/wiki/InstallGuide

When running `cmake`, use the option `-DCMAKE_INSTALL_PREFIX:PATH=$HOME/local` to set the correct installation location.

When finished compiling, install it with `make install`.


## Compiling IGraph/M

To compile IGraph/M, you will need:

 - A C++ compiler with C++11 support.  On Linux, use GCC 4.8 or later (GCC 4.8.5 is known to work). On macOS, use the command line tools of Xcode 6 or later.
 - The [LTemplate Mathematica package][1].  Please download and install it.
 - git, for cloning the repository.

Then follow these steps:

 1. Clone the IGraphM repository and check out the `master` branch (do not use any other branch).
 2. Edit `BuildSettings.m` and set the path to your igraph installation, where necessary.  The available options are the same as for [CreateLibrary](http://reference.wolfram.com/language/CCompilerDriver/ref/CreateLibrary.html).
 3. Append the repository's root directory (i.e. the same directory where this `README.md` file is found) to Mathematica's `$Path`. At the same time, ensure that you do not have the paclet version installed. If you do, remove it with `PacletUninstall["IGraphM"]` before proceeding.
 4. Load the package with ``<<IGraphM` ``.  It should automatically be compiled. It can also be recompiled using ``IGraphM`Developer`Recompile[]``.

## Adding new functions to IGraph/M

If you decide to join IGraph/M development, please contact me for guidance.  The basics are covered below.

IGraph/M uses the LTemplate package to interface with C++. First please go through the [LTemplate tutorial][1].

#### Simple functions

To see the basic structure of a function, let us look at how a simple one such as `IGConnectedQ` is implemented.  On the C++ side we add a new member function to the `IG` class:

```c++
bool connectedQ(bool strong) const {
    igraph_bool_t res;
    igCheck(igraph_is_connected(&graph, &res, strong ? IGRAPH_STRONG : IGRAPH_WEAK));
    return res;
}
```

The `igCheck` function must wrap any igraph function that returns error codes.  It ensures that the error is passed on to Mathematica and reported appropriately.

On the Mathematica side we add the function to the library template in `IGraphM.m` as

```mma
LFun["connectedQ", {True|False (* strongly connected *)}, True|False]
```

And finally we write a Mathematica function that calls this member function on an `IG` object. It is found in `Connectivity.m`.

```mathematica
IGConnectedQ[g_?igGraphQ] :=
    catch@Block[{ig = igMake[g]},
        check@ig@"connectedQ"[True]
    ]
```

The pattern `g_?igGraphQ` will only match `Graph` objects compatible with igraph.  `igMake` will convert a Mathematica graph to an `IG` object.  Calls to `IG` member functions must be wrapped in `check`, which will catch any library errors and `Throw` them as exceptions.  `catch` is always used at the outermost level to catch and handle these exceptions.

`igMake` is able to handle and convert to igraph format any kind of Mathematica graph for which `igGraphQ` returns `True`.  It also transfers edge weights to igraph. 
 
When weights are not required, use `igMakeUnweighted` for better performance. `igMakeUnweighted` could be used for `IGConnectedQ`.

**Note:** In the past, `igMake` used a much less efficient method to convert graphs. Special alternative functions were used for those cases when edge ordering did not need to be preserved. Thus you might come access `igMakeFast` and `igMakeFastWeighted` in the code base. Consider these equivalent to `igMakeUnweighted` and `igMake`, respectively.  

#### Working with edge weights

`igMake` will automatically transfer edge weights to the `IG` object.  To pass these weights to an igraph function as an `igraph_vector_t *`, use `passWeights()`.  See the `betweenness` member function for a simple example.  `passWeights` returns `NULL` for unweighted graphs.

#### Creating graphs

Let us look at a simple function that creates graphs, such as `IGGraphAtlas`, defined in `DeterministicGenerators.m`. The C++ member function is

```c++
void graphAtlas(mint n) {
    destroy();
    igConstructorCheck(igraph_atlas(&graph, n));
}
```

On the C++ side, the graph creator function must first `destroy()` the igraph graph that is currently stored in the `IG` object.  `IG` objects _always_ store some graph, and are initialized with an empty one.

Now we are ready to call the igraph function that creates the graph (`igraph_atlas`).

Instead of `igCheck` we must use `igConstructorCheck` to handle igraph error codes.  This function ensures that the graph is not left in an inconsistent state.

On the Mathematica side, we write

```mathematica
IGGraphAtlas[n_?Internal`NonNegativeMachineIntegerQ, opt : OptionsPattern[Graph]] :=
    catch@Block[{ig = igMakeEmpty[]},
      check@ig@"graphAtlas"[n];
      applyGraphOpt[opt]@igToGraph[ig]
    ]
```

Since the `IG` object is not created from an existing graph, we use `igMakeEmpty[]` to make an empty one. `igMakeEmpty[]` is equivalent to `Make["IG"]`, a construct you will be familiar with from [the LTemplate tutorial][1].

`igToGraph` converts an `IG` object back to a Mathematica graph. 

Most IGraph/M functions that create graphs take all standard `Graph` options such as `GraphLayout`, `VertexStyle`, etc.  These can be applied using `applyGraphOpt`, as in the example above.

#### Functions that do not operate on the `igraph_t` data structure

Functions that would otherwise be free-standing (not member functions) should go in the `IGlobal` class.


## Running and writing unit tests

IGraph/M's tests rely on the [MicroTest][2] utility package. Please install this first.

To run the current tests, evaluate the cells in `Tests/RunTests.nb` one by one. If no red text is printed in the notebook then the tests have passed.

To add more tests, edit `Tests.m`. Unlike the rest of the IGraph/M sources, it is recommended to edit this file with the Mathematica Front End.

The basic format of a test is `MT[expression, expected]`. For example,

```mathematica
MT[
  1+1,
  2
]
```

If evaluation generates messages, the test will fail. To specify expected messages, use

```mathematica
MT[
  1/0,
  ComplexInfinity,
  {Power::infy}
]
```

An `MT` expression is normally inert, and does not evaluate. To evaluate a test, wrap it in `MTRun`. To evaluate all tests in `Tests.m`, `MTRun@Get["Tests.m"]` is used.


  [1]: https://github.com/szhorvat/LTemplate/
  [2]: https://github.com/szhorvat/MicroTest
