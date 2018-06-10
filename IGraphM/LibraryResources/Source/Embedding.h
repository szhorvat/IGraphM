#ifndef EMBEDDING_H
#define EMBEDDING_H

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
struct Embedding {
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

    mma::IntTensorRef dualGraph() const {
        auto faceNode = findFaces();

        std::vector<std::vector<mint>> nodeFace;

        for (mint i=0; i < faceNode.size(); ++i)
            for (const auto &node : faceNode[i]) {
                if (node >= nodeFace.size())
                    nodeFace.resize(node+1);
                nodeFace[node].push_back(i);
            }

        typedef std::pair<mint, mint> arc;

        std::map<arc, mint> conn;

        for (mint i=0; i < faceNode.size(); ++i)
            for (const auto &node : faceNode[i])
                for (const auto &face : nodeFace[node])
                    if (i != face)
                        conn[arc(std::min(i, face), std::max(i, face))] += 1;

        auto edges =
                mma::makeVector<mint>(
                    1 + 2*std::count_if(conn.begin(), conn.end(), [](const std::pair<arc, mint> &p) { return p.second >= 4; })
                );

        edges[0] = faceNode.size();
        mint i=0;
        for (const auto &p : conn)
            if (p.second >= 4) {
                edges[1 + 2*i    ] = p.first.first;
                edges[1 + 2*i + 1] = p.first.second;
                i++;
            }

        return edges;
    }
};

#endif // EMBEDDING_H
