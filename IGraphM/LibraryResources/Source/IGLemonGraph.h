/*
 * Copyright (c) 2018-2020 Szabolcs Horvát.
 *
 * See the file LICENSE.txt for copying permission.
 */

#ifndef IG_LEMON_GRAPH_H
#define IG_LEMON_GRAPH_H

#include "IGEmbedding.h"

#include "LTemplate.h"
#include "mlstream.h"

// WolframLibrary.h defines True and False, which conflicts with LEMON.
// Thus we undefine them here.
#undef True
#undef False
#include <lemon/static_graph.h>
#include <lemon/planarity.h>
#include <lemon/connectivity.h>
#include <lemon/matching.h>
#include <lemon/adaptors.h>

#include <vector>
#include <utility>
#include <algorithm>
#include <functional>
#include <numeric>

class IGLemonGraph {
    typedef lemon::StaticDigraph Digraph;
    typedef lemon::Undirector<const Digraph> Graph;

    Digraph digraph; // underlying directed graph
    Graph graph; // undirected graph adaptor

    Digraph::ArcMap<mint> edgeIndex; // mapping back to the original edge indices

public:
    IGLemonGraph() :
        graph(digraph),
        edgeIndex(digraph)
    { }

    void fromEdgeList(mma::IntMatrixRef edgeList, mint vertexCount) {

        std::vector<mint> indices(edgeList.rows());
        std::iota(indices.begin(), indices.end(), 0);

        std::sort(
            indices.begin(), indices.end(),
            [&edgeList](mint i, mint j) { return edgeList(i,0) < edgeList(j,0); }
        );

        std::vector<std::pair<int, int>> edges;
        edges.reserve(indices.size());
        for (const auto &i : indices)
            edges.push_back(std::pair<int,int>(edgeList(i,0), edgeList(i,1)));

        digraph.build(vertexCount, edges.begin(), edges.end());

        for (int i=0; i < indices.size(); ++i)
            edgeIndex[digraph.arc(i)] = indices[i];
    }


    /* Planarity */

    bool planarQ() const {
        return lemon::checkPlanarity(graph);
    }

    mma::IntTensorRef kuratowskiSubgraph() const {
        lemon::PlanarEmbedding<Graph> embedding(graph);

        if (embedding.run(true))
            return mma::makeVector<mint>(0);

        std::vector<mint> kur;
        for (Graph::EdgeIt e(graph); e != lemon::INVALID; ++e)
            if (embedding.kuratowski(e))
                kur.push_back(edgeIndex[e]);

        return mma::makeVector<mint>(kur.size(), kur.data());
    }

    mma::IntTensorRef layoutPlanar() const {        

        // work around crash in LEMON under the following conditions
        if (digraph.arcNum() == 0 && digraph.nodeNum() < 3) {
            auto coord = mma::makeMatrix<mint>(digraph.nodeNum(), 2);
            for (int i=0; i < digraph.nodeNum(); ++i) {
                coord(i,0) = 0;
                coord(i,1) = i;
            }
            return coord;
        }

        lemon::PlanarEmbedding<Graph> embedding(graph);

        if (! embedding.run(false))
            throw mma::LibraryError("layoutPlanar: The graph is not planar.");

        lemon::PlanarDrawing<Graph> drawing(graph);
        drawing.run(embedding.embeddingMap());

        auto coord = mma::makeMatrix<mint>(digraph.nodeNum(), 2);

        for (Graph::NodeIt n(graph); n != lemon::INVALID; ++n) {
            coord(digraph.index(n), 0) = drawing[n].x;
            coord(digraph.index(n), 1) = drawing[n].y;
        }

        return coord;
    }

    mma::IntTensorRef planarEmbedding() const {
        lemon::PlanarEmbedding<Graph> embedding(graph);

        if (! embedding.run())
            throw mma::LibraryError("planarEmbedding: The graph is not planar.");

        std::vector<std::vector<mint>> orderings(digraph.nodeNum());

        for (int i=0; i < digraph.nodeNum(); ++i) {
            Graph::Arc a;
            graph.firstOut(a, digraph.node(i));
            if (a == lemon::INVALID)
                continue;

            Graph::Arc arc = a;
            do {
                orderings[i].push_back(digraph.index(graph.target(arc)));
                arc = embedding.next(arc);
            } while (arc != a);
        }

        mint total_length = 1;
        for (const auto &el : orderings)
            total_length += el.size() + 1;

        auto result = mma::makeVector<mint>(total_length);
        result[0] = orderings.size();
        mint k=1;
        for (const auto &el : orderings)
            result[k++] = el.size();
        for (const auto &el : orderings) {
            std::copy(el.begin(), el.end(), result.data() + k);
            k += el.size();
        }

        return result;
    }

    // embedding must be a valid combinatorial embedding for the graph
    mma::IntTensorRef embeddingToCoordinates(const IGEmbedding &embedding) const {

        // work around crash in LEMON under the following conditions
        if (digraph.arcNum() == 0 && digraph.nodeNum() < 3) {
            auto coord = mma::makeMatrix<mint>(digraph.nodeNum(), 2);
            for (int i=0; i < digraph.nodeNum(); ++i) {
                coord(i,0) = 0;
                coord(i,1) = i;
            }
            return coord;
        }

        if (! embedding.planarQ(lemon::countConnectedComponents(graph)))
            throw mma::LibraryError("embeddingToCoordinates: The embedding is not planar.");

        const auto &emb = embedding.embedding;
        Graph::ArcMap<Graph::Arc> arcmap(graph);

        for (mint v1=0; v1 < emb.size(); ++v1) {
            for (mint i=0; i < emb[v1].size(); ++i) {
                Graph::Arc a1 = lemon::findArc(graph, digraph.node(v1), digraph.node(emb[v1][i]));
                Graph::Arc a2 = lemon::findArc(graph, digraph.node(v1), digraph.node(emb[v1][(i+1) % emb[v1].size()]));
                if (a1 == lemon::INVALID || a2 == lemon::INVALID)
                    throw mma::LibraryError("embeddingToCoordinates: invalid graph");
                arcmap[a1] = a2;
            }
        }

        lemon::PlanarDrawing<Graph> drawing(graph);

        drawing.run(arcmap);

        auto coord = mma::makeMatrix<mint>(digraph.nodeNum(), 2);
        for (Graph::NodeIt n(graph); n != lemon::INVALID; ++n) {
            coord[2*digraph.index(n)  ] = drawing[n].x;
            coord[2*digraph.index(n)+1] = drawing[n].y;
        }

        return coord;
    }


    /* Matching */

    mma::IntTensorRef maximumMatching() const {

        lemon::MaxMatching<Graph> matching(graph);
        matching.run();

        std::vector<mint> mat;
        for (Graph::EdgeIt e(graph); e != lemon::INVALID; ++e)
            if (matching.matching(e))
                mat.push_back(edgeIndex[e]);

        return mma::makeVector<mint>(mat.size(), mat.data());
    }

    mint matchingNumber() const {
        lemon::MaxMatching<Graph> matching(graph);
        matching.run();
        return matching.matchingSize();
    }
};

#endif // IG_LEMON_GRAPH_H
