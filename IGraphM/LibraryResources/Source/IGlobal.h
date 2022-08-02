/*
 * Copyright (c) 2016-2022 Szabolcs Horvát.
 *
 * See the file LICENSE.txt for copying permission.
 */

#ifndef IGLOBAL_H
#define IGLOBAL_H

#include "IGCommon.h"
#include "IGRNG.h"

#include <set>
#include <map>
#include <functional>
#include <cmath>


igraph_warning_handler_t igWarningHandler;
igraph_error_handler_t igErrorHandler;
igraph_fatal_handler_t igFatalHandler;
igraph_interruption_handler_t igInterruptionHandler;
igraph_progress_handler_t  igProgressHandler;

extern double progress_reporting_granularity;

class IGlobal {
public:
    void init() {
        igraph_set_error_handler(igErrorHandler);
        igraph_set_warning_handler(igWarningHandler);
        igraph_set_interruption_handler(igInterruptionHandler);
        igraph_set_fatal_handler(igFatalHandler);

        rngInit();
    }

    ~IGlobal() { }

    const char *version() {
        const char *ver;
        igraph_version(&ver, nullptr, nullptr, nullptr);
        return ver;
    }

    const char *compilationDate() { return __DATE__; }

    /* Progress handler */

    void setProgressReporting(bool enabled, double granularity) {
        if (granularity < 0) {
            throw mma::LibraryError("Progress reporting granularity must not be negative.");
        }
        progress_reporting_granularity = granularity;
        if (enabled) {
            igraph_set_progress_handler(igProgressHandler);
        } else {
            igraph_set_progress_handler(nullptr);
            igProgressHandler("", 0.0, nullptr);
        }
    }

    /* Random number generators */

    void setRandomGenerator(mint id) { rngSet(id); }

    void seedRandom(mint s) { rngSeed(s); }

    const char *randomGeneratorName() { return rngName(); }

    /* Functions useful in argument and result checking */

    // Does the array contain any Inf or NaN?
    bool infOrNanQ(mma::RealTensorRef t) {
        for (const auto &el : t)
            if (std::isnan(el) || std::isinf(el))
                return true;
        return false;
    }

    // Does the array contain positive numbers only?
    bool posArrQ(mma::RealTensorRef t) {
        for (const auto &el : t)
            if (el <= 0)
                return false;
        return true;
    }

    // Does the array contain negative numbers only?
    bool negArrQ(mma::RealTensorRef t) {
        for (const auto &el : t)
            if (el >= 0)
                return false;
        return true;
    }

    // Does the array contain non-positive numbers only?
    bool nonPosArrQ(mma::RealTensorRef t) {
        for (const auto &el : t)
            if (el > 0)
                return false;
        return true;
    }

    // Does the array contain non-negative numbers only?
    bool nonNegArrQ(mma::RealTensorRef t) {
        for (const auto &el : t)
            if (el < 0)
                return false;
        return true;
    }


    /* Graph related functions that do not use the IG data structure */

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

    bool graphicalQ(mma::RealTensorRef outdeg, mma::RealTensorRef indeg, bool directed, bool loops, bool multi) {
        igraph_vector_t ig_outdeg = igVectorView(outdeg);
        igraph_vector_t ig_indeg  = igVectorView(indeg);
        igraph_bool_t res;
        igraph_edge_type_sw_t et = (loops ? IGRAPH_LOOPS_SW : IGRAPH_SIMPLE_SW) | (multi ? IGRAPH_MULTI_SW : IGRAPH_SIMPLE_SW);
        if (! directed)
            igCheck(igraph_is_graphical(&ig_outdeg, nullptr, et, &res));
        else
            igCheck(igraph_is_graphical(&ig_outdeg, &ig_indeg, et, &res));
        return res;
    }

    bool bigraphicalQ(mma::RealTensorRef deg1, mma::RealTensorRef deg2, bool multi) {
        igraph_vector_t ig_deg1 = igVectorView(deg1);
        igraph_vector_t ig_deg2 = igVectorView(deg2);
        igraph_bool_t res;
        igCheck(igraph_is_bigraphical(&ig_deg1, &ig_deg2, multi ? IGRAPH_MULTI_SW : IGRAPH_SIMPLE_SW, &res));
        return res;
    }

    // Recognize the degree sequence of a split graph.
    // https://en.wikipedia.org/wiki/Split_graph#Degree_sequences
    // Assumption: d[i] >= 0 for all i.
    bool splitQ(mma::IntTensorRef d) {
        std::sort(d.begin(), d.end(), std::greater<mint>());

        mint dsum1 = 0;
        mint m = 0;
        for (; m < d.size(); ++m) {
            if (d[m] < m)
                break;
            dsum1 += d[m];
        }

        mint dsum2 = 0;
        for (mint i=m; i < d.size(); ++i) {
            dsum2 += d[i];
        }

        return dsum1 == m*(m-1) + dsum2;
    }

    // Recognize the degree sequence of a threshold graph
    // Algorithm: Repeatedly remove isolated or dominating vertices.
    // If nothing remains, the graph is a threshold graph.
    bool thresholdQ(mma::IntTensorRef d) {
        std::sort(d.begin(), d.end(), std::greater<mint>());

        mint n = d.size();
        mint l = 0, h = n-1;
        mint offset = 0; // after some vertex removals, actual degrees are 'offset' lower than stored degrees

        if (n == 0)
            return true;

        if (d[l] > n-1 || d[h] < 0)
            return false;

        while (n > 0) {
            // Remvove isolated vertices
            while (d[h] == offset && n > 0) {
                h--;
                n--;
            }

            if (n == 0)
                break;

            // There must be a dominating vertex now
            if (d[l] == n-1 + offset) {
                l++;
                n--;
                offset += 1;
            } else {
                return false;
            }
        }

        return true;
    }

    double compareCommunities(mma::RealTensorRef c1, mma::RealTensorRef c2, mint m) const {
        igraph_community_comparison_t method;
        switch (m) {
        case 0: method = IGRAPH_COMMCMP_VI; break;
        case 1: method = IGRAPH_COMMCMP_NMI; break;
        case 2: method = IGRAPH_COMMCMP_SPLIT_JOIN; break;
        case 3: method = IGRAPH_COMMCMP_RAND; break;
        case 4: method = IGRAPH_COMMCMP_ADJUSTED_RAND; break;
        default: throw mma::LibraryError("Invalid community comparison method.");
        }

        igraph_vector_t comm1 = igVectorView(c1);
        igraph_vector_t comm2 = igVectorView(c2);
        double res;
        igCheck(igraph_compare_communities(&comm1, &comm2, &res, method));
        return res;
    }

    // Compute an edge list from a sparse incidence matrix returned by IncidenceMatrix[]
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

    // Sorts pairs of integers (an edge list)
    // Used for canonicalizing an undirected edge list, or for directed->unirected conversion
    mma::IntTensorRef edgeListSortPairs(mma::IntMatrixRef pairs) {
        if (pairs.cols() != 2)
            throw mma::LibraryError("sortPairs: n-by-2 matrix expected.");
        for (mint i=0; i < pairs.rows(); ++i)
            if (pairs(i,0) > pairs(i,1))
                std::swap(pairs(i,0), pairs(i,1));
        return pairs;
    }

    /*
    // Remove self-loops from an edge list
    mma::IntMatrixRef edgeListRemoveLoops(mma::IntMatrixRef pairs) {
        if (pairs.cols() != 2)
            throw mma::LibraryError("edgeListRemoveLoops: n-by-2 matrix expected.");

        std::vector<mint> result;
        for (mint i=0; i < pairs.rows(); ++i)
            if (pairs(i,0) != pairs(0,1)) {
                result.push_back(pairs(i,0));
                result.push_back(pairs(i,1));
            }
        return mma::makeMatrix<mint>(result.size() / 2, 2, result.data());
    }
    */

    /*
    mma::IntTensorRef edgeListNonLoopPositions(mma::IntMatrixRef pairs) {
        mint m = pairs.rows();

        std::vector<mint> vec;
        vec.reserve(m);
        for (mint i=0; i < m; ++i)
            if (pairs(0,i) != pairs(i,1))
                vec.push_back(i);

        return mma::makeVector<mint>(vec.size(), vec.data());
    }
    */

    // Mark edges when either endpoint is present in 'vertices'
    mma::IntTensorRef edgeListMarkWhenEitherPresent(mma::IntMatrixRef pairs, mma::IntTensorRef vertices) {
        std::set<mint> verts(vertices.begin(), vertices.end());
        auto markers = mma::makeVector<mint>(pairs.rows());

        for (mint i=0; i < pairs.rows(); ++i)
            markers[i] = static_cast<mint>(verts.find(pairs(i,0)) != verts.end() || verts.find(pairs(i,1)) != verts.end());
        return markers;
    }

    // Mark edges when both endpoints are present in 'vertices'
    mma::IntTensorRef edgeListMarkWhenBothPresent(mma::IntMatrixRef pairs, mma::IntTensorRef vertices) {
        std::set<mint> verts(vertices.begin(), vertices.end());
        auto markers = mma::makeVector<mint>(pairs.rows());

        for (mint i=0; i < pairs.rows(); ++i)
            markers[i] = static_cast<mint>(verts.find(pairs(i,0)) != verts.end() && verts.find(pairs(i,1)) != verts.end());
        return markers;
    }

    // Recompute indices in an edge list ('pairs') after 'vertices' have just been deleted
    // Precondition: none of the 'vertices' must be present in 'pairs'
    // Used in IGWeightedVertexDelete
    mma::IntMatrixRef edgeListReindexAfterDelete(mma::IntMatrixRef pairs, mma::IntTensorRef vertices) {
        std::sort(vertices.begin(), vertices.end());
        for (auto &el : pairs) {
            auto i = std::upper_bound(vertices.begin(), vertices.end(), el);
            el -= (i - vertices.begin());
        }
        return pairs;
    }

    // Replaces the vertex numbers in 'pairs' by their 1-based index in the 'vertices' array
    // Precondition: all vertices in 'pairs' are also present in 'vertices'
    // Used in IGWeightedSubgraph
    mma::IntMatrixRef edgeListReindex(mma::IntMatrixRef pairs, mma::IntTensorRef vertices) {
        std::map<mint, mint> posIndex;
        for (mint i=0; i < vertices.size(); ++i)
            posIndex.insert({vertices[i], i+1});
        for (auto &el : pairs)
            el = posIndex[el];
        return pairs;
    }

    // Compute the edge list of a symmetric tree where nodes at level k have splits[k] children
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

    /* Matrix functions */

    template<typename T>
    mma::TensorRef<T> takeLower(mma::MatrixRef<T> mat) const {
        const mint nrow = mat.rows();
        const mint ncol = mat.cols();

        mint len;
        if (nrow > ncol)
            len = ncol*(ncol-1)/2 + ncol*(nrow-ncol);
        else
            len = nrow*(nrow-1)/2;

        auto result = mma::makeVector<T>(len);
        T *r = result.begin();
        for (mint i=0; i < nrow; ++i)
            for (mint j=0; j < i && j < ncol; j++) {
                *r = mat(i,j);
                r++;
            }

        return result;
    }

    mma::IntTensorRef takeLowerInteger(mma::IntMatrixRef mat) { return takeLower(mat); }
    mma::RealTensorRef takeLowerReal(mma::RealMatrixRef mat) { return takeLower(mat); }
    mma::ComplexTensorRef takeLowerComplex(mma::ComplexMatrixRef mat) { return takeLower(mat); }

    template<typename T>
    mma::TensorRef<T> takeUpper(mma::MatrixRef<T> mat) const {
        const mint nrow = mat.rows();
        const mint ncol = mat.cols();

        mint len;
        if (ncol > nrow)
            len = nrow*(nrow-1)/2 + nrow*(ncol-nrow);
        else
            len = ncol*(ncol-1)/2;

        auto result = mma::makeVector<T>(len);
        T *r = result.begin();
        for (mint i=0; i < nrow; ++i)
            for (mint j=i+1; j < ncol; j++) {
                *r = mat(i,j);
                r++;
            }

        return result;
    }

    mma::IntTensorRef takeUpperInteger(mma::IntMatrixRef mat) { return takeUpper(mat); }
    mma::RealTensorRef takeUpperReal(mma::RealMatrixRef mat) { return takeUpper(mat); }
    mma::ComplexTensorRef takeUpperComplex(mma::ComplexMatrixRef mat) { return takeUpper(mat); }

    // Input: index-pair list from a sparse matrix; no. of columns of matrix
    // Output: index of elements that are below the diagonal; index of same elements in result of IGTakeLower
    mma::IntMatrixRef lowerIndexPairPositions(mma::IntMatrixRef pairs, mint cols) {
        std::vector<mint> result;
        result.reserve(2*pairs.rows());
        for (mint i=0; i < pairs.rows(); ++i) {
            mint r = pairs(i,0);
            mint c = pairs(i,1);
            if (r > c) {
                result.push_back(i+1); // index in value list
                result.push_back( ( r > cols ? (cols-1)*cols/2 + cols*(r-cols-1) : (r-1)*(r-2)/2 ) + c  ); // index in final result
            }
        }
        return mma::makeMatrix<mint>(result.size()/2, 2, result.data());
    }

    // See comment for lowerIndexPairPositions(); works identically but used for extracting above-diagonal elements
    mma::IntMatrixRef upperIndexPairPositions(mma::IntMatrixRef pairs, mint cols) {
        std::vector<mint> result;
        result.reserve(2*pairs.rows());
        for (mint i=0; i < pairs.rows(); ++i) {
            mint r = pairs(i,0);
            mint c = pairs(i,1);
            if (r < c) {
                result.push_back(i+1); // index in value list
                result.push_back( (r-1)*cols - (r+1)*r/2 + c ); // index in final result
            }
        }
        return mma::makeMatrix<mint>(result.size()/2, 2, result.data());
    }

    // Input: index-pair list from a sparse matrix
    // Output: index of pairs above diagonal; if 'diag', on-diagonal pairs are included too
    mma::IntTensorRef upperPos(mma::IntMatrixRef pairs, bool diag) {
        std::vector<mint> result;
        result.reserve(pairs.rows());
        for (mint i=0; i < pairs.rows(); ++i) {
            mint r = pairs(i,0);
            mint c = pairs(i,1);
            if (diag) {
                if (r <= c)
                    result.push_back(i+1); // index in value list
            } else {
                if (r < c)
                    result.push_back(i+1); // index in value list
            }
        }
        return mma::makeVector<mint>(result.size(), result.data());
    }

    // Input: index-pair list from a sparse matrix
    // Output: index of pairs not on the diagonal
    mma::IntTensorRef nondiagPos(mma::IntMatrixRef pairs) {
        std::vector<mint> result;
        result.reserve(pairs.rows());
        for (mint i=0; i < pairs.rows(); ++i) {
            mint r = pairs(i,0);
            mint c = pairs(i,1);

            if (r != c)
                result.push_back(i+1); // index in value list
        }
        return mma::makeVector<mint>(result.size(), result.data());
    }


    mma::IntTensorRef fromNauty(const char *str);

    mma::RealTensorRef percolationCurve(mma::IntMatrixRef edges, mint n);
};

#endif // IGLOBAL_H

