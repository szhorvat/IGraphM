/*
 * Copyright (c) 2019 Szabolcs Horvát.
 *
 * See the file LICENSE.txt for copying permission.
 */

#ifndef IG_FLANN3D_H
#define IG_FLANN3D_H

#include "nanoflann.hpp"
#include "IGFlannCommon.h"
#include <LTemplate.h>
#include <vector>

using namespace nanoflann;


class IGFlann3D {

    typedef KDTreeSingleIndexAdaptor<
            L2_Simple_Adaptor<double, PointSet>,
            PointSet,
            3 /* dim */
        > kdtree_t;

    PointSet *ps = nullptr;
    kdtree_t *kdtree = nullptr;

    void clear() {
        delete kdtree;
        delete ps;
    }

public:

    ~IGFlann3D() { clear(); }

    // must be called immediately after creating object
    void setPoints(mma::RealMatrixRef pts) {
        if (pts.cols() != 3)
            throw mma::LibraryError("setPoints: input must be three-dimensional");

        clear();

        ps = new PointSet(pts); // this will clone pts, so it is safe to pass pts as "Constant"
        kdtree = new kdtree_t(3, *ps, KDTreeSingleIndexAdaptorParams(10));

        kdtree->buildIndex();
    }


    // returns indices of points within distance d of pt
    mma::IntTensorRef query(mma::RealTensorRef pt, double d) {
        if (pt.size() != 3)
            throw mma::LibraryError("query: query point must be three-dimensional");
        if (d <= 0)
            throw mma::LibraryError("query: query distance must be positive");

        std::vector<std::pair<size_t, double>> ret_matches;

        SearchParams params;

        kdtree->radiusSearch(pt.data(), d*d, ret_matches, params);

        auto res = mma::makeVector<mint>(ret_matches.size());
        for (mint i=0; size_t(i) < ret_matches.size(); ++i)
            res[i] = ret_matches[i].first + 1;
        return res;
    }


    // like query(), but for multiple query points and distances
    mma::IntTensorRef queryMultiple(mma::RealMatrixRef pts, mma::RealTensorRef dists) {
        if (pts.cols() != 3)
            throw mma::LibraryError("queryMultiple: query points must be three-dimensional");
        if (pts.rows() != dists.size())
            throw mma::LibraryError("queryMultiple: there must be the same number of query distances as query points");

        std::vector<mint> results;
        std::vector<mint> sizes;

        SearchParams params;
        std::vector<std::pair<size_t, double>> ret_matches;

        for (mint i=0; i < pts.rows(); ++i) {
            double d = dists[i];
            if (d <= 0)
                throw mma::LibraryError("queryMultiple: query distances must be positive");

            kdtree->radiusSearch(&pts(i,0), d*d, ret_matches, params);

            sizes.push_back(ret_matches.size());
            for (const auto &el : ret_matches)
                results.push_back(el.first + 1);
        }

        auto res = mma::makeVector<mint>(1 + sizes.size() + results.size());
        mint i=0;
        res[i++] = sizes.size();
        for (const auto &el : sizes)
            res[i++] = el;
        for (const auto &el : results)
            res[i++] = el;
        return res;
    }


    // how many neighbours does each point have within the given distance?
    mma::IntTensorRef neighborCounts(mma::RealMatrixRef pts, mma::RealTensorRef dists) {
        mint n = dists.size();

        if (pts.cols() != 3)
            throw mma::LibraryError("neighborCounts: query points must be three-dimensional");
        if (pts.rows() != n)
            throw mma::LibraryError("neighborCounts: there must be the same number of query distances as query points");

        SearchParams params;
        std::vector<std::pair<size_t, double>> ret_matches;

        auto res = mma::makeVector<mint>(n);

        for (mint i=0; i < n; ++i) {
            double d = dists[i];
            if (d <= 0) {
                res.free(); // result vector will not be returned, so free it
                throw mma::LibraryError("neighborCounts: query distances must be positive");
            }

            kdtree->radiusSearch(&pts(i,0), d*d, ret_matches, params);

            res[i] = ret_matches.size();
        }

        return res;
    }


    // how many points are within the given radius of both centres?
    mma::IntTensorRef intersectionCounts(mma::RealMatrixRef centres1, mma::RealMatrixRef centres2, mma::RealTensorRef dists) {
        mint n = dists.size();

        if (centres1.cols() != 3)
            throw mma::LibraryError("intersectionCounts: query points must be three-dimensional");
        if (centres2.cols() != 3)
            throw mma::LibraryError("intersectionCounts: query points must be three-dimensional");
        if (centres1.rows() != n || centres2.rows() != n)
            throw mma::LibraryError("intersectionCounts: there must be the same number of query distances as query points");

        SearchParams params;
        std::vector<std::pair<size_t, double>> ret_matches;

        auto res = mma::makeVector<mint>(n);

        for (mint i=0; i < n; ++i) {
            double d = dists[i];
            if (d <= 0) {
                res.free(); // result vector will not be returned, so free it
                throw mma::LibraryError("intersectionCounts: query distances must be positive");
            }

            kdtree->radiusSearch(&centres1(i,0), d*d, ret_matches, params);

            mint count = 0;
            for (const auto &el : ret_matches) {
                double pd2 =
                      sqr(ps->kdtree_get_pt(el.first, 0) - centres2(i, 0)) +
                      sqr(ps->kdtree_get_pt(el.first, 1) - centres2(i, 1)) +
                      sqr(ps->kdtree_get_pt(el.first, 2) - centres2(i, 2));
                if (pd2 < d*d)
                    count++;
            }

            res[i] = count;
        }

        return res;
    }


    // how many points are within the given radius of either centre?
    mma::IntTensorRef unionCounts(mma::RealMatrixRef centres1, mma::RealMatrixRef centres2, mma::RealTensorRef dists) {
        mint n = dists.size();

        if (centres1.cols() != 3)
            throw mma::LibraryError("unionCounts: query points must be three-dimensional");
        if (centres2.cols() != 3)
            throw mma::LibraryError("unionCounts: query points must be three-dimensional");
        if (centres1.rows() != n || centres2.rows() != n)
            throw mma::LibraryError("unionCounts: there must be the same number of query distances as query points");

        SearchParams params;
        std::vector<std::pair<size_t, double>> ret_matches;

        auto res = mma::makeVector<mint>(n);

        for (mint i=0; i < n; ++i) {
            double d = dists[i];
            if (d <= 0) {
                res.free(); // result vector will not be returned, so free it
                throw mma::LibraryError("unionCounts: query distances must be positive");
            }

            // the number of points within |A union B| = |A| + |B| - |A intersect B|

            kdtree->radiusSearch(&centres2(i,0), d*d, ret_matches, params);
            mint count = ret_matches.size();

            kdtree->radiusSearch(&centres1(i,0), d*d, ret_matches, params);
            count += ret_matches.size();

            for (const auto &el : ret_matches) {
                double pd2 =
                      sqr(ps->kdtree_get_pt(el.first, 0) - centres2(i, 0)) +
                      sqr(ps->kdtree_get_pt(el.first, 1) - centres2(i, 1)) +
                      sqr(ps->kdtree_get_pt(el.first, 2) - centres2(i, 2));
                if (pd2 < d*d)
                    count--;
            }

            res[i] = count;
        }

        return res;
    }

};

#endif // IG_FLANN3D_H
