#ifndef IGCOMMON_H
#define IGCOMMON_H

extern "C" { // workaround for igraph_version() C++ compatibility bug in igraph <= 0.7.1
#include <igraph/igraph.h>
}

#include "LTemplate.h"

#include <algorithm>



inline igraph_vector_t igVectorView(mma::RealTensorRef t) {
    static double dummy = 0.0; // work around igraph not liking zero-length vectors will NULL pointers
    igraph_vector_t vec;
    mint len = t.length();
    igraph_vector_view(&vec, len == 0 ? &dummy : t.data(), len);
    return vec;
}


struct igVector {
    igraph_vector_t vec;

    igVector() { igraph_vector_init(&vec, 0); }
    ~igVector() { igraph_vector_destroy(&vec); }

    long length() const { return igraph_vector_size(&vec); }

    igraph_real_t *begin() { return &(VECTOR(vec)[0]); }
    igraph_real_t *end() { return begin() + length(); }

    const igraph_real_t *begin() const { return &(VECTOR(vec)[0]); }
    const igraph_real_t *end() const { return begin() + length(); }

    void clear() { igraph_vector_clear(&vec); }

    void copyFromMTensor(mma::RealTensorRef t) {
        igraph_vector_t from = igVectorView(t);
        igraph_vector_update(&vec, &from);
    }

    mma::RealTensorRef makeMTensor() const {
        mma::RealTensorRef res = mma::makeVector<double>(length());
        std::copy(begin(), end(), res.begin());
        return res;
    }
};


struct igMatrix {
    igraph_matrix_t mat;

    igMatrix() { igraph_matrix_init(&mat, 0, 0); }
    ~igMatrix() { igraph_matrix_destroy(&mat); }

    long length() const { return mat.data.end - mat.data.stor_begin; }
    long nrow() const { return mat.nrow; }
    long ncol() const { return mat.ncol; }

    igraph_real_t *begin() { return mat.data.stor_begin; }
    igraph_real_t *end() { return mat.data.end; }

    const igraph_real_t *begin() const { return mat.data.stor_begin; }
    const igraph_real_t *end() const { return mat.data.end; }

    mma::RealMatrixRef makeMTensor() const {
        mma::RealMatrixRef res = mma::makeMatrix<double>(mat.nrow, mat.ncol);
        std::copy(begin(), end(), res.begin());
        return res;
    }
};


inline void igCheck(int err) {
    if (! err) return;
    throw mma::LibraryError(igraph_strerror(err));
}


#endif // IGCOMMON_H

