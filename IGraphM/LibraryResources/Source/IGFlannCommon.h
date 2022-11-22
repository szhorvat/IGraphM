/*
 * Copyright (c) 2019-2022 Szabolcs Horvát.
 *
 * See the file LICENSE.txt for copying permission.
 */

#ifndef IGFLANNCOMMON_H
#define IGFLANNCOMMON_H

#include "nanoflann.hpp"
#include <LTemplate.h>

// Square a number
template<typename T>
inline T sqr(T x) { return x*x; }

// Eculidean distance in 2D
template<typename T>
inline T sqdist2(const T *pt1, const T *pt2) {
    return sqr(pt1[0] - pt2[0]) + sqr(pt1[1] - pt2[1]);
}

// Eculidean distance in 3D
template<typename T>
inline T sqdist3(const T *pt1, const T *pt2) {
    return sqr(pt1[0] - pt2[0]) + sqr(pt1[1] - pt2[1]) + sqr(pt1[2] - pt2[2]);
}


// Point set to be used with nanoflann
// Supports both 2D and 3D point sets
struct PointSet {
    mma::RealMatrixRef pts;

    PointSet(mma::RealMatrixRef pts) : pts(pts.clone()) { }

    ~PointSet() { pts.free(); }

    size_t kdtree_get_point_count() const { return pts.rows(); }

    double kdtree_get_pt(const size_t idx, const size_t dim) const {
        return pts(idx, dim);
    }

    template <class BBOX>
    bool kdtree_get_bbox(BBOX& /* bb */) const { return false; }
};

#endif // IGFLANNCOMMON_H
