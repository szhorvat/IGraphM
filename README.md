[![Discourse topics](https://img.shields.io/discourse/topics?color=limegreen&server=https%3A%2F%2Figraph.discourse.group)](https://igraph.discourse.group)
[![GitHub (pre-)release](https://img.shields.io/github/release/szhorvat/IGraphM/all.svg)](https://github.com/szhorvat/IGraphM/releases)
[![GitHub All Releases](https://img.shields.io/github/downloads/szhorvat/IGraphM/total.svg)](https://github.com/szhorvat/IGraphM/releases)
[![Contributions welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg)](https://github.com/szhorvat/IGraphM#contributions)
[![DOI](https://zenodo.org/badge/41793262.svg)](https://zenodo.org/badge/latestdoi/41793262)

# [IGraph/M – igraph for _Mathematica_][main]

The IGraph/M package provides a [Wolfram Language](https://www.wolfram.com/) interface to the popular [igraph](https://igraph.org/) network analysis and graph theory library, as well as many other useful functions for working with graphs in _Mathematica_.  Check out the [blog post][main] for an overview.


## 📦 Installation

The system requirements are _Mathematica_ 10.0.2 or later, 64-bit Windows/macOS/Linux, or Raspberry Pi.

To install the package automatically, simply evaluate the following in _Mathematica_:

```mathematica
Get["https://raw.githubusercontent.com/szhorvat/IGraphM/master/IGInstaller.m"]
```

IGraph/M can also be installed manually in the same way as any _Mathematica_ application distributed as a paclet.

Download the `.paclet` file from [the GitHub releases page](https://github.com/szhorvat/IGraphM/releases), and [install it using the `PacletInstall` function in Mathematica](http://mathematica.stackexchange.com/q/141887/12).  For example, assuming that the file `IGraphM-0.5.1.paclet` was downloaded into the directory `~/Downloads`, evaluate:

```mathematica
Needs["PacletManager`"]
PacletInstall["~/Downloads/IGraphM-0.5.1.paclet"]
```

After installation, the package can now be loaded with:

```mathematica
Needs["IGraphM`"]
```

Verify that it works by evaluating `IGVersion[]`, then continue to the documentation with  the `IGDocumentation[]` command.

To uninstall all currently installed versions of IGraph/M, evaluate `PacletUninstall["IGraphM"]`. This will remove all traces of IGraph/M from your system.


## 📖 Documentation

To open the documentation notebook, evaluate:

```mathematica
Needs["IGraphM`"]
IGDocumentation[]
```

Alternatively, search for "igraphm" in _Mathematica_'s Documentation Center.

[A web-based preview](http://szhorvat.net/mathematica/IGDocumentation/) is also available.

The documentation is not yet complete and contributions are very welcome.  If you would like to help out with expanding the documentation, send me an email.

For additional details about functions, also check [the igraph documentation pages](http://igraph.org/c/doc/).


## 🛠️ Contributions

**Help wanted with editing documentation and writing unit tests! Only basic _Mathematica_ knowledge is required for this.**

Contributions to IGraph/M are very welcome!  _Mathematica_ programmers of all levels of expertise can contribute.

In order of increasing difficulty, help is needed with the following tasks:

 - Just play with the package, read the documentation, and try to find problems.
 - Create examples for the documentation or edit the documentation.
 - Write formal unit tests.
 - Implement new functions in pure Wolfram Language code.
 - Expose more igraph library functions to IGraph/M (C++ knowledge required)
 - Implement entirely new functions in C++.

If you are interested in contributing, send me an email for guidance. Evaluate the following in _Mathematica_ to get my email address:

    Uncompress["1:eJxTTMoPChZiYGAorsrILypLLHFIz03MzNFLzs8FAG/xCOI="]

C++ programmers should look at [Development.md](Development.md) for additional information.


## 💡Additional information

 - [Support and Troubleshooting](SUPPORT.md)
 - [Revision history](CHANGELOG.md)

## ⚖️ License

IGraph/M is licensed under the [GNU GPLv3](https://opensource.org/licenses/gpl-3.0.html). See [LICENSE.txt](LICENSE.txt) for details.

 [ltemplate]: https://github.com/szhorvat/LTemplate/
 [main]: http://szhorvat.net/mathematica/IGraphM

