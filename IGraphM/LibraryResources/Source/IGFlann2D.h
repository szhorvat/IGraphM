/*
 * Copyright (c) 2018-2022 Szabolcs Horvát.
 *
 * See the file LICENSE.txt for copying permission.
 */

#ifndef IG_FLANN2D_H
#define IG_FLANN2D_H

#include "nanoflann.hpp"
#include "IGFlannCommon.h"
#include "IGCommon.h"
#include <vector>

using namespace nanoflann;


class IGFlann2D {

    typedef KDTreeSingleIndexAdaptor<
            L2_Simple_Adaptor<double, PointSet>,
            PointSet,
            2, /* dim */
            mint /* AccessorType */
        > kdtree_t;

    PointSet *ps = nullptr;
    kdtree_t *kdtree = nullptr;

    void clear() {
        delete kdtree;
        delete ps;
    }


    // Helper result class to count points within a neighbourhood of a given radius.
    // Set short_circuit = true to quickly determine whether the neighbourhood
    // contains at least one point. In this case size() returns 0 if no point
    // was found, and 1 if at least one point was found.
    class NeighborCounts {
        const double radius; // L2 search radius
        const bool short_circuit; // whether to stop after one point has been found
        const mint v1, v2; // edge endpoits; excluded from the count
        const PointSet *ps;

        size_t count;

    public:
        NeighborCounts(double radius_, bool short_circuit_,
                       mint v1_, mint v2_, const PointSet *ps_)
            : radius(radius_), short_circuit(short_circuit_), v1(v1_), v2(v2_), ps(ps_)
        {
            init();
        }

        void init() { clear(); }
        void clear() { count = 0; }

        size_t size() const { return count; }

        bool full() const { return true; }

        bool addPoint(double dist, mint index) {
            if (dist < radius && index != v1 && index != v2) {
                count++;
                if (short_circuit)
                    return false;
            }
            return true;
        }

        double worstDist() const { return radius; }
    };


    // Helper result class for listing the points in the lune's bounding disk
    // and counting how many are within the lune.
    class IntersectionCounts {
        const double radius; // L2 search radius
        const double beta_radius;
        const bool short_circuit; // whether to stop after one point has been found
        const mint v1, v2; // edge endpoits; excluded from the count
        const double *centre1, *centre2; // disc centres
        const PointSet *ps;

        size_t count;

    public:
        IntersectionCounts(
                    double radius_, double beta_radius_,
                    bool short_circuit_,
                    mint v1_, mint v2_, const double *centre1_, const double *centre2_, const PointSet *ps_)
            : radius(radius_), beta_radius(beta_radius_),
              short_circuit(short_circuit_),
              v1(v1_), v2(v2_),
              centre1(centre1_), centre2(centre2_), ps(ps_)
        {
            init();
        }

        void init() { clear(); }
        void clear() { count = 0; }

        size_t size() const { return count; }

        bool full() const { return true; }

        bool addPoint(double dist, mint index) {
            if (dist < radius && index != v1 && index != v2) {
                double pd1 = sqdist2(&ps->pts(index,0), centre1);
                double pd2 = sqdist2(&ps->pts(index,0), centre2);
                if (pd1 < beta_radius && pd2 < beta_radius) {
                    count++;
                    if (short_circuit)
                        return false;
                }
            }
            return true;
        }

        double worstDist() const { return radius; }
    };

public:

    ~IGFlann2D() { clear(); }

    // must be called immediately after creating object
    void setPoints(mma::RealMatrixRef pts) {
        if (pts.cols() != 2)
            throw mma::LibraryError("setPoints: input must be two-dimensional");

        clear();

        ps = new PointSet(pts); // this will clone pts, so it is safe to pass pts as "Constant"
        kdtree = new kdtree_t(2, *ps, KDTreeSingleIndexAdaptorParams(10));

        kdtree->buildIndex();
    }


    // returns indices of points within distance d of pt
    mma::IntTensorRef query(mma::RealTensorRef pt, double d) {
        if (pt.size() != 2)
            throw mma::LibraryError("query: query point must be two-dimensional");
        if (d <= 0)
            throw mma::LibraryError("query: query distance must be positive");

        std::vector<std::pair<mint, double>> ret_matches;

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
        std::vector<std::pair<mint, double>> ret_matches;

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
    // if short_circuit == true then 1 will be returned whenever at least one point exists
    // edge endpoints are excluded
    mma::IntTensorRef neighborCounts(
            mma::RealMatrixRef pts,
            mma::RealTensorRef dists,
            mma::IntMatrixRef edges,
            bool short_circuit) {

        mint n = dists.size();

        if (pts.cols() != 2)
            throw mma::LibraryError("neighborCounts: query points must be two-dimensional");
        if (pts.rows() != n)
            throw mma::LibraryError("neighborCounts: there must be the same number of query distances as query points");
        if (edges.cols() != 2)
            throw mma::LibraryError("neighborCounts: edge matrix must have two columns");
        if (edges.rows() != n)
            throw mma::LibraryError("neighborCounts: there must be the same number of edges points as query points");

        SearchParams params;

        auto res = mma::makeVector<mint>(n);
        LTGuard<mma::IntTensorRef> res_guard(res);

        for (mint i=0; i < n; ++i) {
            double d = dists[i];
            if (d <= 0) {
                throw mma::LibraryError("neighborCounts: query distances must be positive");
            }

            NeighborCounts nc(d*d, short_circuit, edges(i,0)-1, edges(i,1)-1, ps);
            res[i] = kdtree->radiusSearchCustomCallback(&pts(i,0), nc, params);
        }

        res_guard.deactivate();
        return res;
    }


    // how many points are within the given radius of both centres?
    // if short_circuit == true then 1 will be returned whenever at least one point exists
    // edge endpoints are excluded
    mma::IntTensorRef intersectionCounts(
            mma::RealMatrixRef centres1, mma::RealMatrixRef centres2,
            mma::RealTensorRef dists,
            mma::IntMatrixRef edges,
            bool short_circuit) {

        mint n = dists.size();

        if (centres1.cols() != 2)
            throw mma::LibraryError("intersectionCounts: query points must be two-dimensional");
        if (centres2.cols() != 2)
            throw mma::LibraryError("intersectionCounts: query points must be two-dimensional");
        if (centres1.rows() != n || centres2.rows() != n)
            throw mma::LibraryError("intersectionCounts: there must be the same number of query distances as query points");
        if (edges.cols() != 2)
            throw mma::LibraryError("intersectionCounts: edge matrix must have two columns");
        if (edges.rows() != n)
            throw mma::LibraryError("intersectionCounts: there must be the same number of edges points as query points");

        SearchParams params;

        auto res = mma::makeVector<mint>(n);
        LTGuard<mma::IntTensorRef> res_guard(res); // automatically free 'mat' upon premature exit from the function

        for (mint i=0; i < n; ++i) {
            double d = dists[i];
            if (d <= 0) {
                throw mma::LibraryError("intersectionCounts: query distances must be positive");
            }

            // midpoint between the two centres
            double c[2] = { 0.5*(centres1(i,0) + centres2(i,0)),
                            0.5*(centres1(i,1) + centres2(i,1)) };

            // squared radius of the lune's bounding disk
            // this is the squared half height of the lune
            double r2 = d*d - sqdist2(c, &centres1(i,0));

            IntersectionCounts ic(r2, d*d, short_circuit, edges(i,0)-1, edges(i,1)-1, &centres1(i, 0), &centres2(i, 0), ps);
            res[i] = kdtree->radiusSearchCustomCallback(c, ic, params);

            mma::check_abort();
        }

        res_guard.deactivate();
        return res;
    }


    // how many points are within the given radius of either centre?
    mma::IntTensorRef unionCounts(
            mma::RealMatrixRef centres1, mma::RealMatrixRef centres2,
            mma::RealTensorRef dists,
            mma::IntMatrixRef edges,
            bool short_circuit) {

        if (! short_circuit)
            throw mma::LibraryError("unionCounts: short_circuit=false is not yet implemented");

        mint n = dists.size();

        if (centres1.cols() != 2)
            throw mma::LibraryError("unionCounts: query points must be two-dimensional");
        if (centres2.cols() != 2)
            throw mma::LibraryError("unionCounts: query points must be two-dimensional");
        if (centres1.rows() != n || centres2.rows() != n)
            throw mma::LibraryError("unionCounts: there must be the same number of query distances as query points");
        if (edges.cols() != 2)
            throw mma::LibraryError("unionCounts: edge matrix must have two columns");
        if (edges.rows() != n)
            throw mma::LibraryError("unionCounts: there must be the same number of edges points as query points");

        SearchParams params;
        params.sorted = false;

        auto res = mma::makeVector<mint>(n);
        LTGuard<mma::IntTensorRef> res_guard(res);

        for (mint i=0; i < n; ++i) {
            double d = dists[i];
            if (d <= 0) {
                throw mma::LibraryError("unionCounts: query distances must be positive");
            }

            // full counting not currently implemented; an implementation can be based on:
            // the number of points within |A union B| = |A| + |B| - |A intersect B|

            NeighborCounts nc(d*d, /* short_circuit */ true, edges(i,0)-1, edges(i,1)-1, ps);
            res[i] = kdtree->radiusSearchCustomCallback(&centres1(i,0), nc, params);

            if (! res[i]) {
                nc.clear();

                res[i] = kdtree->radiusSearchCustomCallback(&centres2(i,0), nc, params);
            }
        }

        res_guard.deactivate();
        return res;
    }

};

#endif // IG_FLANN2D_H
