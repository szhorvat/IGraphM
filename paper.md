---
title: 'IGraph/M: graph theory and network analysis for Mathematica'
tags:
  - Mathematica
  - Wolfram Language
  - graph theory
  - network analysis
authors:
  - name: Szabolcs Horvát
    orcid: 0000-0002-3100-523X
    corresponding: true
    affiliation: "1, 2"
  - name: Jakub Podkalicki
    affiliation: 3
  - name: Gábor Csárdi
    orcid: 0000-0001-7098-9676
    affiliation: 4
  - name: Tamás Nepusz
    orcid: 0000-0002-1451-338X
    affiliation: 5
  - name: Vincent Traag
    orcid: 0000-0003-3170-3879
    affiliation: 6
  - name: Fabio Zanini
    orcid: 0000-0001-7097-8539
    affiliation: 7
  - name: Daniel Noom
    affiliation: 8
affiliations:
 - name: Max Planck Institute for Cell Biology and Genetics, Dresden, Germany \newline
   index: 1
 - name: Center for Systems Biology Dresden, Dresden, Germany \newline
   index: 2
 - name: Independent Developer, Poland \newline
   index: 3
 - name: RStudio \newline
   index: 4
 - name: Independent Developer, Hungary \newline
   index: 5
 - name: Centre for Science and Technology Studies, Leiden University, Leiden, Netherlands \newline
   index: 6
 - name: Lowy Cancer Research Centre, University of New South Wales, Kensington, NSW, Australia \newline
   index: 7
 - name: Independent Developer, Netherlands \newline
   index: 8

date: 29 August 2022
bibliography: paper.bib

---

# Summary

IGraph/M[^1] is an efficient general purpose graph theory and network analysis package for Mathematica [@mathematica]. IGraph/M serves as the Wolfram Language interfaces to the igraph C library [@igraph; @Csardi2006], and also provides several unique pieces of functionality not yet present in igraph, but made possible by combining its capabilities with Mathematica's. The package is designed to support both graph theoretical research as well as the analysis of large-scale empirical networks.


# Statement of need

Mathematica contains [extensive built-in functionality for working with graphs](https://reference.wolfram.com/language/guide/GraphsAndNetworks.html). IGraph/M extends this graph framework with many new functions not otherwise available in the Wolfram Language, and also provides alternative and more featureful open-source implementations of many of Mathematica's existing built-ins. This makes it possible for Wolfram Language users to easily double-check results, just as Python and R users can already do thanks to the multiple different graph packages available in those languages. This is particularly useful in graph theory where many results are just as difficult to verify as to compute.

The only other independent graph theory package for Mathematica was [Combinatorica](https://reference.wolfram.com/language/Combinatorica/guide/CombinatoricaPackage) [@combinatorica], which has been mostly unmaintained since the introduction of the `Graph` expression type with the release of Mathematica 8.0 in 2011. Despite this, not all of Combinatorica's functions have built-in equivalents in Mathematica. IGraph/M provides replacements for almost all of these old Combinatorica functions while offering much better performance thanks to being implemented in a mix of C and C++ instead of the Wolfram Language.


# Design goals and features

One of the major appeals of Mathematica is its tightly integrated nature: different functionality areas of the system can smoothly and seamlessly interoperate with each other. In order to preserve this benefit, a major design goal for IGraph/M was to integrate well into the rest of the system. This is achieved by working directly with Mathematica's native `Graph` data type, which is transparently converted to igraph's internal format as needed. This makes IGraph/M different from igraph's other high-level interfaces: igraph's internal graph data structure is not exposed to users and the package does not provide trivial operations which are already present in Mathematica, such as adding or removing vertices and edges. Instead, priority is given to functionality that delivers a true benefit over Mathematica's own built-ins. While some of IGraph/M's functions may appear to simply duplicate built-in functions, almost all of them provide additional features not otherwise available in Mathematica. For example, unlike the built-in `BetweenessCentrality[]` function, `IGBetweenness[]` supports weighted graphs; unlike the built-in `PageRankCentrality[]`, `IGPageRank[]` takes into consideration self-loops and parallel edges; in contrast to `IsomorphicGraphQ`, `IGIsomorphicQ` supports non-simple graphs as well as vertex and edge colours; etc. In order to make it easy to distinguish built-in symbols from those of IGraph/M, all functions in this package have names starting with `IG`. The naming and interface of functions is chosen so as to be familiar both to Mathematica users and users of igraph's Python and R interfaces.

IGraph/M fully leverages the igraph C library's capabilities to integrate into high-level host languages: Almost all computations are interruptible in the usual manner. This feature is particularly important for a graph theory package that provides multiple algorithms with exponential or slower time complexity. Despite being implemented in a compiled language, IGraph/M uses Mathematica's built-in random number generator by default. Therefore, all its stochastic algorithms respect seeds set with `SeedRandom[]` or `BlockRandom[]` the same way as built-in functions would.

IGraph/M aims to exploit the interactive features of Mathematica notebooks to improve user productivity. In this spirit, it includes an interactive graph editor, `IGGraphEditor[]`, and supports dynamically displaying the progress of many functions.


# Use cases and unique features

IGraph/M provides multiple unique features that are not present in the core igraph library. Examples include exact graph colouring [@VanGelder2008], functions for working with planar graphs and combinatorial embeddings, proximity graph functions [@Toussaint1980; @Kirkpatrick1985] and aids for working with spatial networks, functions for performing tests related to the graph automorphism group, and several others. Some features, such as `IGRealizeDegreeSequence`, are based on original research of the authors [@Horvat2021]. Additionally, there are many convenience functions that are helpful when working with Mathematica's `Graph`, including a simplified system for working with edge and vertex attributes based on the concept of mapping functions over attribute values (`IGEdgeMap`, `IGVertexMap`).

The following examples show a few of these features while also demonstrating the concise idioms made possible by this package. The output is shown in \autoref{fig:fig1}. Many more examples are found [in the package's documentation](http://szhorvat.net/mathematica/IGDocumentation/).

Load the package:

```
In[1]:= << IGraphM`
Out[1]= IGraph/M 0.6.2 (July 25, 2022)
        Evaluate IGDocumentation[] to get started.
```

Create a random maze as a random spanning tree of a regular lattice confined to a hexagon, and colour it according to betweenness centrality:

```
In[2]:=
g = IGMeshGraph@IGLatticeMesh["CairoPentagonal", Polygon@CirclePoints[6, 6]];
t = IGRandomSpanningTree[g,
        VertexCoordinates -> GraphEmbedding[g], GraphStyle -> "ThickEdge"];
IGEdgeMap[ColorData["Rainbow"], EdgeStyle -> IGEdgeBetweenness /* Rescale, t]
```

Compute and visualize a minimum vertex and edge colouring of the Gabriel graph of a set of spatial points:

```
In[3]:=
IGGabrielGraph[RandomPoint[Disk[], 30],
    GraphStyle -> "Monochrome", VertexSize -> 1, EdgeStyle -> Thickness[1/40]
] //
  IGVertexMap[ColorData[106], VertexStyle -> IGMinimumVertexColoring] //
  IGEdgeMap[ColorData[106], EdgeStyle -> IGMinimumEdgeColoring]
```

![Output of the above example code. `Out[2]` is on the left, `Out[3]` on the right.\label{fig:fig1}](fig1.pdf){ width=90% }


# Implementation notes

IGraph/M is built using LTemplate [@LTemplate], an open-source system that makes it easier to extend Mathematica using C++ code. IGraph/M also serves as the primary driver of LTemplate development. Planar graph-related functionality is implemented using the LEMON graph library [@lemon; @Dezso2011]. Spatial graph functions make use of the nanoflann nearest neighbour search library [@nanoflann].


# Acknowledgements

This project has been made possible in part by grant number 2021-237461 from the Chan Zuckerberg Initiative DAF, an advised fund of Silicon Valley Community Foundation. We thank Juho Lauri for feedback and ideas, as well as for guidance with the implementation of the SAT-based exact graph colouring functionality. Henrik Schumacher provided the basis for the implementation of the mesh/graph conversion functions and some of the proximity graph functions. Patrick Scheibe made available his [Wolfram Language plugin for IntelliJ IDEA](http://wlplugin.halirutan.de/), without which this project would have been much more difficult to manage. Sz.\ H.\ is grateful to Carl Modes for ongoing support and feedback regarding the paper, and to the [Wissenchaftskolleg zu Berlin](https://www.wiko-berlin.de/) for time to work on spatial networks functionality in this package.


# References

 [^1]: Available at [szhorvat.net/mathematica/IGraphM](http://szhorvat.net/mathematica/IGraphM) and [github.com/szhorvat/IGraphM](https://github.com/szhorvat/IGraphM)
