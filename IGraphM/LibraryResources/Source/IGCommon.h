#ifndef IGCOMMON_H
#define IGCOMMON_H

extern "C" { // workaround for igraph_version() C++ compatibility bug in igraph <= 0.7.1
#include <igraph/igraph.h>
}

#include "mlstream.h"
#include "LTemplate.h"

#include <algorithm>
#include <string>
#include <sstream>
#include <type_traits>


static_assert(std::is_same<double, igraph_real_t>::value, "IGraphM assumes igraph_real_t to be double.");


inline igraph_vector_t igVectorView(mma::RealTensorRef t) {
    static double dummy = 0.0; // work around igraph not liking zero-length vectors will NULL pointers
    igraph_vector_t vec;
    mint len = t.length();
    igraph_vector_view(&vec, len == 0 ? &dummy : t.data(), len);
    return vec;
}


// RAII for igraph_vector_t
class igVector {
    bool moved;

    igVector(const igVector &source); // avoid accidental expensive implicit copy

public:
    igraph_vector_t vec;

    igVector() : moved(false) { igraph_vector_init(&vec, 0); }
    igVector(igVector &&source) : moved(false) { vec = source.vec; source.moved = true; }
    igVector(const igraph_vector_t *source) : moved(false) { igraph_vector_copy(&vec, source); }
    ~igVector() { if (!moved) igraph_vector_destroy(&vec); }

    long length() const { return vec.end - vec.stor_begin; }

    igraph_real_t *begin() { return vec.stor_begin; }
    igraph_real_t *end() { return vec.end; }

    const igraph_real_t *begin() const { return vec.stor_begin; }
    const igraph_real_t *end() const { return vec.end; }

    void clear() { igraph_vector_clear(&vec); }
    void resize(long newsize) { igraph_vector_resize(&vec, newsize); }

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


// RAII for igraph_vector_t
// note that igraph_integer_t and mint may not be the same type
struct igIntVector {
    igraph_vector_int_t vec;

    igIntVector() { igraph_vector_int_init(&vec, 0); }
    ~igIntVector() { igraph_vector_int_destroy(&vec); }

    long length() const { return vec.end - vec.stor_begin; }

    igraph_integer_t *begin() { return vec.stor_begin; }
    igraph_integer_t *end() { return vec.end; }

    const igraph_integer_t *begin() const { return vec.stor_begin; }
    const igraph_integer_t *end() const { return vec.end; }

    void clear() { igraph_vector_int_clear(&vec); }
    void resize(long newsize) { igraph_vector_int_resize(&vec, newsize); }

    void copyFromMTensor(mma::IntTensorRef t) {
        resize(t.length());
        std::copy(t.begin(), t.end(), begin());
    }

    mma::IntTensorRef makeMTensor() const {
        mma::IntTensorRef res = mma::makeVector<mint>(length());
        std::copy(begin(), end(), res.begin());
        return res;
    }
};


// RAII for igraph_maxtrix_t
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

    void copyFromMTensor(mma::RealMatrixRef t) {
        igraph_vector_t from = igVectorView(t);
        igraph_vector_update(&mat.data, &from);
        mat.nrow = t.cols();
        mat.ncol = t.rows();
        igraph_matrix_transpose(&mat);
    }

    mma::RealMatrixRef makeMTensor() const {
        mma::RealMatrixRef res = mma::makeMatrix<double>(mat.ncol, mat.nrow);
        std::copy(begin(), end(), res.begin());
        return res;
    }
};


// RAII for igraph_vector_ptr_t
struct igList {
    igraph_vector_ptr_t list;

    igList() {
        igraph_vector_ptr_init(&list, 0);
        IGRAPH_VECTOR_PTR_SET_ITEM_DESTRUCTOR(&list, igraph_vector_destroy);
    }
    ~igList() { igraph_vector_ptr_destroy_all(&list); }

    long length() const { return igraph_vector_ptr_size(&list); }

    void push(igraph_vector_t *vec) { igraph_vector_ptr_push_back(&list, vec); }
};


// extend mlstream with igraph-specific types

inline mlStream & operator << (mlStream &ml, const igraph_vector_t &vec) {
    if (! MLPutReal64List(ml.link(), vec.stor_begin, vec.end - vec.stor_begin))
        ml.error("cannot return vector");
    return ml;
}


inline mlStream & operator << (mlStream &ml, const igVector &vec) { return ml << vec.vec; }


inline mlStream & operator << (mlStream &ml, const igList &list) {
    long len = list.length();
    if (! MLPutFunction(ml.link(), "List", len))
        ml.error("cannot return vector list");
    for (int i=0; i < len; ++i)
        ml << *static_cast<igraph_vector_t *>(VECTOR(list.list)[i]);
    return ml;
}


inline mlStream & operator >> (mlStream &ml, igIntVector &vec) {
    int *data;
    int length;
    // igraph_integer_t is an int, so we try to use the corresponding MathLink type: Integer32
    if (! MLGetInteger32List(ml.link(), &data, &length))
        ml.error("Integer32List expected");
    vec.resize(length);
    std::copy(data, data+length, vec.begin());
    MLReleaseInteger32List(ml.link(), data, length);
    return ml;
}


// check igraph's error codes and abort if necessary
inline void igCheck(int err) {
    if (! err) return;
    std::ostringstream msg;
    msg << "igraph returned with error: " << igraph_strerror(err);
    throw mma::LibraryError(msg.str());
}


#endif // IGCOMMON_H

