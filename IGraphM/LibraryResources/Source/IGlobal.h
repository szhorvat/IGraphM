/*
 * Copyright (c) 2017 Szabolcs Horvát.
 *
 * See the file LICENSE.txt for copying permission.
 */

#ifndef IGLOBAL_H
#define IGLOBAL_H

#include "IGCommon.h"

#include <random>
#include <set>
#include <map>
#include <functional>
#include <cmath>


void igWarningHandler(const char *reason, const char *file, int line, int igraph_errno);
void igErrorHandler(const char *reason, const char *file, int line, int igraph_errno);
int  igInterruptionHandler(void *);


class IGlobal {
public:
    void init() {
        igraph_set_error_handler(igErrorHandler);
        igraph_set_warning_handler(igWarningHandler);
        igraph_set_interruption_handler(igInterruptionHandler);
        igCheck(igraph_rng_seed(igraph_rng_default(), std::random_device()()));
    }

    ~IGlobal() { }

    const char *version() {
        const char *ver;
        igCheck(igraph_version(&ver, nullptr, nullptr, nullptr));
        return ver;
    }

    const char *compilationDate() { return __DATE__; }

    void seedRandom(mint s) {
        igCheck(igraph_rng_seed(igraph_rng_default(), s));
    }

    // Graph related functions that do not use the graph data structure

    // Fast implementation based on Z. Kiraly: Recognizing graphic degree sequences and generating all realizations.
    bool erdosGallai(mma::IntTensorRef deg) {
        if (deg.length() == 0)
            return true;

        mint n, w, b, s, c;
        n = deg.length();

        mint sum = 0;
        for (const auto &el : deg)
            sum += el;
        if (sum % 2 == 1)
            return false;

        std::sort(deg.begin(), deg.end(), std::greater<mint>());

        w = n; b = 0; s = 0; c = 0;
        for (mint k=1; k <= n; ++k) {
            b += deg[k-1];
            c += w-1;
            while (w > k && deg[w-1] <= k) {
                s += deg[w-1];
                c -= k;
                w--;
            }
            if (b > c+s)
                return false;
            else if (w == k)
                return true;
        }
        massert(false);
    }

    bool graphicalQ(mma::RealTensorRef outdeg, mma::RealTensorRef indeg) {
        igraph_vector_t ig_outdeg = igVectorView(outdeg);
        igraph_vector_t ig_indeg  = igVectorView(indeg);
        igraph_bool_t res;
        if (indeg.length() == 0)
            igCheck(igraph_is_graphical_degree_sequence(&ig_outdeg, nullptr, &res));
        else
            igCheck(igraph_is_graphical_degree_sequence(&ig_outdeg, &ig_indeg, &res));
        return res;
    }

    bool infOrNanQ(mma::RealTensorRef t) {
        for (double *x = t.begin(); x != t.end(); ++x)
            if (std::isnan(*x) || std::isinf(*x))
                return true;
        return false;
    }

    mma::IntTensorRef incidenceToEdgeList(mma::SparseMatrixRef<mint> im, bool directed) {
        auto edgeList = mma::makeMatrix<mint>(im.cols(), 2);
        if (directed) {
            for (auto it = im.begin(); it != im.end(); ++it) {
                switch (*it) {
                case -1:
                    edgeList[2*it.col()] = it.row();
                    break;
                case  1:
                    edgeList[2*it.col() + 1] = it.row();
                    break;
                case  2:
                case -2:
                    edgeList[2*it.col()] = it.row();
                    edgeList[2*it.col() + 1] = it.row();
                    break;
                default:
                    throw mma::LibraryError("Invalid incidence matrix.");
                }
            }
        } else {
            for (auto &el : edgeList)
                el = -1;
            for (auto it = im.begin(); it != im.end(); ++it) {
                switch (*it) {
                case  1:
                    if (edgeList[2*it.col()] == -1)
                        edgeList[2*it.col()] = it.row();
                    else
                        edgeList[2*it.col() + 1] = it.row();
                    break;
                case  2:
                    edgeList[2*it.col()] = it.row();
                    edgeList[2*it.col() + 1] = it.row();
                    break;
                default:
                    throw mma::LibraryError("Invalid incidence matrix.");
                }
            }
        }
        return edgeList;
    }

    // Sorts pairs of integers; used for canonicalizing an undirected edge list, or for directed->unirected conversion
    mma::IntTensorRef edgeListSortPairs(mma::IntMatrixRef pairs) {
        if (pairs.cols() != 2)
            throw mma::LibraryError("sortPairs: n-by-2 matrix expected.");
        for (int i=0; i < pairs.rows(); ++i)
            if (pairs(i,0) > pairs(i,1))
                std::swap(pairs(i,0), pairs(i,1));
        return pairs;
    }   

    mma::IntMatrixRef removeSelfLoops(mma::IntMatrixRef pairs) {
        if (pairs.cols() != 2)
            throw mma::LibraryError("removeSelfLoops: n-by-2 matrix expected.");

        std::vector<mint> result;
        for (int i=0; i < pairs.rows(); ++i)
            if (pairs(i,0) != pairs(0,1)) {
                result.push_back(pairs(i,0));
                result.push_back(pairs(i,1));
            }
        return mma::makeMatrix<mint>(result.size() / 2, 2, result.data());
    }

    mma::IntTensorRef edgeListMarkVertices1(mma::IntMatrixRef pairs, mma::IntTensorRef vertices) {
        std::set<mint> verts(vertices.begin(), vertices.end());
        auto markers = mma::makeVector<mint>(pairs.rows());

        for (int i=0; i < pairs.rows(); ++i)
            markers[i] = static_cast<mint>(verts.find(pairs(i,0)) != verts.end() || verts.find(pairs(i,1)) != verts.end());
        return markers;
    }

    mma::IntTensorRef edgeListMarkVertices2(mma::IntMatrixRef pairs, mma::IntTensorRef vertices) {
        std::set<mint> verts(vertices.begin(), vertices.end());
        auto markers = mma::makeVector<mint>(pairs.rows());

        for (int i=0; i < pairs.rows(); ++i)
            markers[i] = static_cast<mint>(verts.find(pairs(i,0)) != verts.end() && verts.find(pairs(i,1)) != verts.end());
        return markers;
    }

    mma::IntMatrixRef edgeListDecVertices(mma::IntMatrixRef pairs, mma::IntTensorRef vertices) {
        std::sort(vertices.begin(), vertices.end());
        for (auto &el : pairs) {
            auto i = std::upper_bound(vertices.begin(), vertices.end(), el);
            el -= (i - vertices.begin());
        }
        return pairs;
    }

    mma::IntMatrixRef edgeListReindex(mma::IntMatrixRef pairs, mma::IntTensorRef vertices) {
        std::map<mint, mint> posIndex;
        for (mint i=0; i < vertices.size(); ++i)
            posIndex.insert({vertices[i], i+1});
        for (auto &el : pairs)
            el = posIndex[el];
        return pairs;
    }

    mma::IntTensorRef symmetricTree(mma::IntTensorRef splits) {
        mint vcount = 1;
        mint c = 1;
        for (const auto &s : splits) {
            c *= s;
            vcount += c;
        }
        auto edges = mma::makeMatrix<mint>(vcount-1, 2);

        mint j1=0, j2=1;
        mint p=1;
        mint ec=0;
        for (const auto &s : splits) {
            for (mint i=0; i < p; ++i) {
                for (mint j=0; j < s; ++j) {
                    edges(ec, 0) = j1;
                    edges(ec, 1) = j2;
                    ec++;
                    j2++;
                }
                j1++;
                mma::check_abort();
            }
            p *= s;
        }

        return edges;
    }
};

#endif // IGLOBAL_H

