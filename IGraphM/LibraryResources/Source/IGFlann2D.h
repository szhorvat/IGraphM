/*
 * Copyright (c) 2018 Szabolcs Horvát.
 *
 * See the file LICENSE.txt for copying permission.
 */

#include "nanoflann.hpp"
#include <LTemplate.h>
#include <vector>

using namespace nanoflann;


template<typename T> inline T sqr(T x) { return x*x; }


class PointSet2D {
    mma::RealMatrixRef pts;

public:
    PointSet2D(mma::RealMatrixRef pts) : pts(pts.clone()) { }

    ~PointSet2D() { pts.free(); }

    size_t kdtree_get_point_count() const { return pts.rows(); }

    double kdtree_get_pt(const size_t idx, const size_t dim) const {
        if (dim == 0)
            return pts(idx,0);
        else
            return pts(idx,1);
    }

    template <class BBOX>
    bool kdtree_get_bbox(BBOX& /* bb */) const { return false; }
};


class IGFlann2D {

    typedef KDTreeSingleIndexAdaptor<
            L2_Simple_Adaptor<double, PointSet2D>,
            PointSet2D,
            2 /* dim */
        > kdtree_t;

    PointSet2D *ps = nullptr;
    kdtree_t *kdtree = nullptr;

    void clear() {
        delete kdtree;
        delete ps;
    }

public:

    ~IGFlann2D() { clear(); }

    // must be called immediately after creating object
    void setPoints(mma::RealMatrixRef pts) {
        if (pts.cols() != 2)
            throw mma::LibraryError("setPoints: input must be two-dimensional");

        clear();

        ps = new PointSet2D(pts); // this will clone pts, so it is safe to pass pts as "Constant"
        kdtree = new kdtree_t(2, *ps, KDTreeSingleIndexAdaptorParams(10));

        kdtree->buildIndex();
    }


    // returns indices of points within distance d of pt
    mma::IntTensorRef query(mma::RealTensorRef pt, double d) {
        if (pt.size() != 2)
            throw mma::LibraryError("query: query point must be two-dimensional");
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
        if (pts.cols() != 2)
            throw mma::LibraryError("queryMultiple: query points must be two-dimensional");
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

        if (pts.cols() != 2)
            throw mma::LibraryError("neighborCounts: query points must be two-dimensional");
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

        if (centres1.cols() != 2)
            throw mma::LibraryError("intersectionCounts: query points must be two-dimensional");
        if (centres2.cols() != 2)
            throw mma::LibraryError("intersectionCounts: query points must be two-dimensional");
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
                      sqr(ps->kdtree_get_pt(el.first, 1) - centres2(i, 1));
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

        if (centres1.cols() != 2)
            throw mma::LibraryError("unionCounts: query points must be two-dimensional");
        if (centres2.cols() != 2)
            throw mma::LibraryError("unionCounts: query points must be two-dimensional");
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
                      sqr(ps->kdtree_get_pt(el.first, 1) - centres2(i, 1));
                if (pd2 < d*d)
                    count--;
            }

            res[i] = count;
        }

        return res;
    }

};
