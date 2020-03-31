## Support and Troubleshooting for [IGraph/M](README.md)

For support, use the [official igraph forum](https://igraph.discourse.group/) or [Mathematica StackExchange](https://mathematica.stackexchange.com/).

IGraph/M is currently under development, and a few bugs are to be expected.  However, I try not to release a new version until most problems I know of are fixed.  If you do find a problem, please [open an issue on GitHub](https://github.com/szhorvat/IGraphM/issues). **In the bug report include the output of ``IGraphM`Developer`GetInfo[]``.**


#### Common problems

  * **"Cannot open ``IGraphM` ``"**

    This message will be shown in the following situations:

    - IGraph/M is not installed. Please follow the installation instructions above carefully.

    - IGraph/M is not compatible with your system. Please review the requirements in the Installation section above.  Additional symptoms will be that `PacletFind["IGraphM"]` returns `{}` but `PacletFind["IGraphM", "MathematicaVersion" -> All, "SystemID" -> All]` returns a non-empty list.

  * **"Cannot open ``LTemplate`LTemplatePrivate` ``"**

    Loading fails with this error if one tries to use the `master` branch (development branch) without the necessary dependencies.

    Please download IGraph/M from the [GitHub releases page instead](https://github.com/szhorvat/IGraphM/releases).  It includes everything needed to run the package.

    Do not clone the git repository and do not use the `master` branch unless you want to develop IGraph/M.

  * **"No build settings found. Please check `BuildSettings.m`"**

    This error may be shown when trying to load IGraph/M on an incompatible platform.  Currently, the following platforms are supported: Windows 64-bit, OS X 10.9 or later, Linux x86_64 with glibc 2.14 or later, Raspberry Pi with Raspbian Stretch.

    It may be possible to run IGraph/M on other platforms, however, it will be necessary to compile it from source.

    If you see this error on a supported platform, please file a bug report and include the output of ``IGraphM`Developer`GetInfo[]``.

  * **`The function IGlobal_get_collection was not loaded from the file ...` or similar error appears on Linux**

    Check that you are using a Linux distribution with glibc 2.14 or later. This is the minimum requirement. If your Linux distribution does have the required glibc version, then please open a new issue.

#### Known issues and workarounds

   * `Graph[Graph[...], ...]` is returned

     Sometimes layout functions may return an expression which looks like

         Graph[ Graph[...], VertexCoordinates -> {...} ]

     or similar. A property does not get correctly applied to the graph.  This is due to a bug in Mathematica. I believe I have worked around most of these issues, but if you encounter them, one possible workaround is to cycle the graph `g` through some other representation, e.g. `g = Uncompress@Compress[g]`.

   * The graphs returned by `IGBipartiteGameGNM` and `IGBipartiteGameGNP` may not render when using the `DirectedEdges -> True` and `"Bidirectional" -> True` options.  This is due to a bug in _Mathematica_'s `"BipartiteEmbedding"` graph layout and can be corrected by passing `GraphLayout -> Automatic` to these functions. Alternatively, they may be rendered using `IGLayoutBipartite`.

   * `IGDocumentation[]` may not work with _Mathematica_ 11.1 on Linux.  Please enter `IGraphM/IGDocumentation` in the Documentation Center address bar to access it (or simply search for `igraph` in the Documentation Center).

   * `IGZeroDiagonal` will crash with non-square sparse matrices in _Mathematica_ 11.1 and earlier. This is a bug in those versions of _Mathematica_.

   * See also https://github.com/szhorvat/IGraphM/issues
