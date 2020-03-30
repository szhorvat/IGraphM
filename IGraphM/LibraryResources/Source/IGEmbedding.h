/*
 * Copyright (c) 2018-2020 Szabolcs Horvát.
 *
 * See the file LICENSE.txt for copying permission.
 */

#ifndef IG_EMBEDDING_H
#define IG_EMBEDDING_H

#include "LTemplate.h"
#include "mlstream.h"

#include <vector>
#include <set>
#include <map>
#include <algorithm>

/*
 * This class represents a combinatorial embedding, i.e. a cyclic ordering of the neighbours
 * of each vertex. Vertices are represented by non-negative integers (0, 1, ...), i.e. indexing
 * is 0-based.
 */
struct IGEmbedding {
    std::vector<std::vector<mint>> embedding;

    void set(MLINK link) {
        mlStream ml{link, "embeddingSet"};

        ml >> mlCheckArgs(1) >> embedding;

        ml.newPacket();
        ml << mlSymbol("Null");
    }

    void get(MLINK link) {
        mlStream ml{link, "embeddingGet"};

        ml >> mlCheckArgs(0);

        ml.newPacket();
        ml << embedding;
    }

    /* Verify that the combinatorial embedding is valid for a simple graph, i.e.:
     *  - there are no repeated elements in any of the sublists
     *  - there are no self-loops
     *  - the reverse of every directed arc is also present
     */
    bool validQ() const {
        typedef std::pair<mint, mint> arc;

        std::set<arc> arcs;
        for (mint v1=0; v1 < embedding.size(); ++v1) {
            for (const auto &v2 : embedding[v1]) {
                if (v1 == v2)
                    return false;
                arc a = {v1, v2};
                if (arcs.find(a) != arcs.end())
                    return false;
                arcs.insert(a);
            }
        }

        for (const auto &a : arcs)
            if (arcs.find(arc{a.second, a.first}) == arcs.end())
                return false;

        return true;
    }    

    // Note that redundant faces are returned for non-connected graphs.
    std::vector<std::vector<mint>> findFaces() const {
        std::vector<std::vector<mint>> faces;

        std::vector<std::vector<bool>> visited;
        visited.reserve(embedding.size());
        for (const auto &v : embedding)
            visited.emplace_back(v.size(), false);

        for (mint i=0; i < embedding.size(); ++i)
            for (mint j=0; j < embedding[i].size(); ++j) {
                if (visited[i][j])
                    continue;

                faces.emplace_back();
                auto &face = faces.back();

                mint curr = i;
                mint nexti = j;
                mint prev;
                do {
                    face.push_back(curr);
                    visited[curr][nexti] = true;

                    prev = curr;
                    curr = embedding[curr][nexti];
                    auto it = std::find(embedding[curr].begin(), embedding[curr].end(), prev);

                    if (it == embedding[curr].end())
                        mma::LibraryError("findFaces: invalid embedding.");

                    nexti = (it - embedding[curr].begin() + embedding[curr].size() - 1) % embedding[curr].size();
                } while (! visited[curr][nexti]);
            }

        return faces;
    }

    mma::IntTensorRef faces() const {
        auto faces = findFaces();

        mint total_length = 1;
        for (const auto &el : faces)
            total_length += el.size() + 1;

        auto result = mma::makeVector<mint>(total_length);
        result[0] = faces.size();
        mint k=1;
        for (const auto &el : faces)
            result[k++] = el.size();
        for (const auto &el : faces) {
            std::copy(el.begin(), el.end(), result.data() + k);
            k += el.size();
        }

        return result;
    }

    bool planarQ(mint componentCount) const {
        mint vertexCount = embedding.size();

        mint edgeCount = 0;
        for (const auto &v : embedding)
            edgeCount += v.size();
        edgeCount /= 2;

        if (edgeCount == 0)
            return true;

        mint faceCount = findFaces().size();

        return (vertexCount - edgeCount + faceCount == 2*componentCount);
    }

    mma::IntMatrixRef dualGraph() const {
        auto faces = findFaces();

        // an edge is a pair of vertices
        typedef std::pair<mint, mint> edge;

        // which two faces is an edge incident to?
        // in other words, which edge of the primal transforms into which edge of the dual?
        // the mapping we construct results in a multigraph
        std::map<edge, edge> edgeFace;

        for (mint face=0; face < faces.size(); ++face) {
            const auto &faceVertices = faces[face];
            for (mint j=0; j < faceVertices.size(); ++j) {
                mint s = faceVertices[j];
                mint t = faceVertices[(j+1) % faceVertices.size()];

                if (s < t)
                    edgeFace[edge{s,t}].first = face;
                else
                    edgeFace[edge{t,s}].second = face;
            }
        }

        auto result = mma::makeMatrix<mint>(edgeFace.size(), 2);
        mint i=0;
        for (const auto &el : edgeFace) {
            result(i,0) = el.second.first;
            result(i,1) = el.second.second;
            i++;
        }

        return result;
    }

};

#endif // IG_EMBEDDING_H
