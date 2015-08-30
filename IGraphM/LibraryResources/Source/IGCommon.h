#ifndef IGCOMMON_H
#define IGCOMMON_H

extern "C" { // workaround for igraph_version() C++ compatibility bug in igraph <= 0.7.1
#include <igraph/igraph.h>
}

#include "LTemplate.h"

#include <algorithm>


struct igVector {
    igraph_vector_t vec;

    igVector() { igraph_vector_init(&vec, 0); }
    ~igVector() { igraph_vector_destroy(&vec); }

    long length() const { return igraph_vector_size(&vec); }

    igraph_real_t *begin() { return &(VECTOR(vec)[0]); }
    igraph_real_t *end() { return begin() + length(); }

    const igraph_real_t *begin() const { return &(VECTOR(vec)[0]); }
    const igraph_real_t *end() const { return begin() + length(); }

    mma::RealTensorRef makeMTensor() const {
        mma::RealTensorRef res = mma::makeVector<double>(length());
        std::copy(begin(), end(), res.begin());
        return res;
    }
};


inline igraph_vector_t igVectorView(mma::RealTensorRef &t) {
    igraph_vector_t vec;
    igraph_vector_view(&vec, t.data(), t.length());
    return vec;
}


inline void igCheck(int err) {
    if (! err) return;
    throw mma::LibraryError(igraph_strerror(err));
}


#endif // IGCOMMON_H

