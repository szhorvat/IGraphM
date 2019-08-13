#ifndef IGFLANNCOMMON_H
#define IGFLANNCOMMON_H

#include "nanoflann.hpp"
#include <LTemplate.h>

// Square a number
template<typename T>
inline T sqr(T x) { return x*x; }


// Point set to be used with nanoflann
// Supports both 2D and 3D point sets
class PointSet {
    mma::RealMatrixRef pts;

public:
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
