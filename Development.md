# Development guide for IGraph/M

Contributions to IGraph/M are most welcome!  igraph is a large graph library with diverse functionality.  I primarily focused on providing an interface to functions that I need myself, and I do not have time to cover all igraph functions.  However, the main framework is there, and adding new functions is relatively quick and easy.

If you are interested in extending IGraph/M, send me an email to get technical guidance.  IGraph/M uses the [LTemplate package][1] to simplify writing LibraryLink code, and acts as a driver for LTemplate development.  I recommend starting by reading the LTemplate tutorial.

## Compiling igraph

### OS X and Linux

In this guide I will use `$HOME/local` as the installation directory for the various needed libraries.  You may want to use a different location for your system.

On OS X, set the following environment variable before compiling the libraries:

    export MACOSX_DEPLOYMENT_TARGET=10.9

##### GMP

[Download GMP](https://gmplib.org/), extract it to a directory with no spaces in its path (!) and compile with

    ./configure --prefix=$HOME/local --with-pic
    make
    make check

If the tests have passed, install it with `make install`

##### igraph

Clone [this fork of the igraph](https://github.com/szhorvat/igraph) and check out the `IGraphM` branch. This fork is identical to the main igraph repository, except for a few small temporary patches that the latest version of IGraph/M may depend on. Compile as follows:

    export CPPFLAGS=-I$HOME/local LDFLAGS=-L$HOME/local
    ./bootstrap.sh
    ./configure --prefix=$HOME/local --with-pic  --disable-graphml
    make
    make check

If the tests have passed, install it with `make install`.

### Windows

One option for compiling igraph on Windows is to use an MSYS shell to run the configure script.  Instructions for installing MSYS2 and the MinGW-w64 compiler are found at https://wiki.qt.io/MSYS2.  Install them in a directory with no spaces in its path.  The following instructions assume that libraries will be installed in `$HOME/local`.

Once it is set up, we need to compile GMP.  [Download]((https://gmplib.org/)) and extract it.  Compile using

    ./configure --prefix=$HOME/local --with-pic
    make
    make check

If the tests have passed, install it with `make install`

Clone [this fork of the igraph](https://github.com/szhorvat/igraph) and check out the `IGraphM` branch. This fork is identical to the main igraph repository, except for a few small temporary patches that the latest version of IGraph/M may depend on.  Compile and install as follows:

    export CPPFLAGS="-I$HOME/local -DMSDOS" LDFLAGS=-L$HOME/local
    ./bootstrap.sh
    ./configure --prefix=$HOME/local --disable-graphml
    make
    make install

This will produce a DLL named `libigraph-0.dll` in `$HOME/bin`.  It must be copied into `IGraphM/LibraryResources/Windows-x86-64`.  When using the above version of MinGW-w64, it is also necessary to copy the dependencies `libgcc_s_seh-1.dll`, `libstdc++-6.dll` and `libwinpthread-1.dll` to the same directory.

## Compiling IGraph/M

To compile IGraph/M, you will need:

 - A C++ compiler with C++11 support.  I used GCC 4.9 and clang 3.6.
 - The [LTemplate Mathematica package][1].  Please download and install it.
 - git for cloning the repository.

Then follow these steps:

 1. Clone the IGraphM repository and check out the master branch (do not use the release branch).  
 2. Edit `BuildSettings.m` and set the path to your igraph installation, where necessary.  The available options are the same as for [CreateLibrary](http://reference.wolfram.com/language/CCompilerDriver/ref/CreateLibrary.html).
 3. Append the repository's root directory (i.e. the same directory where this readme file is found) to Mathematica's `$Path`.
 4. Load the package with ``<<IGraphM` ``.  It should automatically be compiled. It can also be recompiled using ``IGraphM`Developer`Recompile[]``.

  [1]: https://bitbucket.org/szhorvat/ltemplate
