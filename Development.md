# Development guide for IGraph/M

Contributions to IGraph/M are most welcome!  

IGraph/M serves both to expose the [igraph C library][3] to the Wolfram Language, as well as to provide independent graph theory and network analysis functionality.

igraph is a large graph library with diverse functionality.  I primarily focused on providing an interface to functions that I need myself, and I do not have time to cover all igraph functions.  However, the main framework is there, and adding new functions is relatively quick and easy.

If you are interested in extending IGraph/M, send me an email to get technical guidance.  IGraph/M uses the [LTemplate package][1] to simplify writing LibraryLink code, and acts as a driver for LTemplate development.  I recommend starting by reading the LTemplate tutorial.

Help is also very welcome with writing or editing documentation, adding examples to the documentation, and writing unit tests.


## Compiling igraph

The following tools are required for building the igraph C library: a C++ compiler with C++11 and C99 support, `cmake`, `flex`, `bison`. See the [installation instructions](https://igraph.org/c/html/latest/igraph-Installation.html) for the basic steps.

On Windows, we use MSVC.

In this guide, I will use `$HOME/local` as the installation directory for the various needed libraries.  You may want to use a different location for your system.

#### Dependencies

All igraph dependencies that are required for features that IGraph/M uses are already included with igraph.

On OS X, set the following environment variable before compiling the libraries:

    export MACOSX_DEPLOYMENT_TARGET=10.9

#### Running CMake

Additional options to pass to CMake:

 - It is recommended to use link-time optimization to achieve good performance, thus use `-DIGRAPH_ENABLE_LTO=ON`.
 - We do not need GraphML support and want to avoid a dependency on libxml2, so use `-DIGRAPH_GRAPHML_SUPPORT=OFF`.
 - On Linux and macOS, use `-DCMAKE_POSITION_INDEPENDENT_CODE=ON`, as we will eventually be creating a shared library.
 - On Windows, use `-DCMAKE_MSVC_RUNTIME_LIBRARY=MultiThreaded`, as this uses the `/MT` option, just like Mathematica's `CreateLibrary`.

#### Compiling for release

When compiling for release, target the appropriate macOS version:

    export MACOSX_DEPLOYMENT_TARGET=10.9

On Linux, make sure that `_GLIBCXX_USE_CXX11_ABI=0` is defined for broader compatibility. All dependencies, such as LEMON, should also be compiled with the same setting.

For CMake, set `-DIGRAPH_OPENMP_SUPPORT=OFF` to disable OpenMP support, as it only provides a modest improvement for PageRank.

## Compiling LEMON

Follow the installation guide of LEMON: http://lemon.cs.elte.hu/trac/lemon/wiki/InstallGuide

When running `cmake`, use the option `-DCMAKE_INSTALL_PREFIX:PATH=$HOME/local` to set the correct installation location.

An example command that also disables unneeded dependencies is

    cmake -DCMAKE_INSTALL_PREFIX=$HOME/local -DCMAKE_BUILD_TYPE=Release -DCMAKE_POSITION_INDEPENDENT_CODE=True -DLEMON_ENABLE_GLPK=False -DLEMON_ENABLE_COIN=False -DLEMON_ENABLE_ILOG=False -DLEMON_ENABLE_SOPLEX=False ..

When finished compiling, install it with `make install`.

When compiling for release, edit `CMakeLists.txt` and set: `CMAKE_MINIMUM_REQUIRED(VERSION 3.15)`. This will allow `DCMAKE_MSVC_RUNTIME_LIBRARY=MultiThreaded` to take effect on Windows. Optionally, set `-DCMAKE_INTERPROCEDURAL_OPTIMIZATION=ON` for better performance.


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

```mma
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

```mma
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

To run the current tests, evaluate the cells in `Tests/RunTests.nb` one by one. A summary is given at the end indicating the number of failures. Failures are also reported in red text.

To add more tests, edit `Tests.m`. Unlike the rest of the IGraph/M sources, it is recommended to edit this file with the Mathematica Front End.

The basic format of a test is `MT[expression, expected]`. For example,

```mma
MT[
  1+1,
  2
]
```

If evaluation generates messages, the test will fail. To specify expected messages, use

```mma
MT[
  1/0,
  ComplexInfinity,
  {Power::infy}
]
```

An `MT` expression is normally inert, and does not evaluate. To evaluate a test, wrap it in `MTRun`. To evaluate all tests in `Tests.m`, `MTRun@Get["Tests.m"]` is used.

The `Tests/MakeTests.nb` notebook provides an easy way to generate tests that can be pasted into `Tests.m`.


  [1]: https://github.com/szhorvat/LTemplate/
  [2]: https://github.com/szhorvat/MicroTest
  [3]: https://igraph.org
