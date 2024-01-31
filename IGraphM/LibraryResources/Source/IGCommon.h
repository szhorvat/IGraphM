/*
 * Copyright (c) 2016-2022 Szabolcs Horvát.
 *
 * See the file LICENSE.txt for copying permission.
 */

#ifndef IG_COMMON_H
#define IG_COMMON_H

#include <igraph/igraph.h>

#include "mlstream.h"
#include "LTemplate.h"

#include <algorithm>
#include <string>
#include <sstream>
#include <type_traits>


static_assert(std::is_same<double, igraph_real_t>::value, "IGraphM assumes igraph_real_t to be double.");
static_assert(sizeof(igraph_integer_t) == sizeof(mint), "IGraphM assumes igraph_integer_t to have the same size as mint.");
static_assert(std::is_same<bool, igraph_bool_t>::value, "IGraphM assumes igraph_bool_t to be bool.");
// See mlstream extractors, which are defined for integer types of particular widths.

// 'mlint' is the MathLink integer type that has the same size as igraph_integer_t
// Note that this may not be identical to igraph_integer_t. For example,
// igraph_integer_t may be 'long long' while mlint may be mlint64 = long,
// both 64-bit, yet considered distinct by the compiler.
#if IGRAPH_INTEGER_SIZE == 32
typedef int mlint;
#define MLGetIGInteger MLGetInteger32
#define MLGetIGIntegerList MLGetInteger32List
#define MLReleaseIGIntegerList MLReleaseInteger32List
#define MLPutIGIntegerList MLPutInteger32List
#define MLGetIGIntegerArray MLGetInteger32Array
#define MLReleaseIGIntegerArray MLReleaseInteger32Array
#define MLPutIGIntegerArray MLPutInteger32Array
#elif IGRAPH_INTEGER_SIZE == 64
typedef mlint64 mlint;
#define MLGetIGInteger MLGetInteger64
#define MLGetIGIntegerList MLGetInteger64List
#define MLReleaseIGIntegerList MLReleaseInteger64List
#define MLPutIGIntegerList MLPutInteger64List
#define MLGetIGIntegerArray MLGetInteger64Array
#define MLReleaseIGIntegerArray MLReleaseInteger64Array
#define MLPutIGIntegerArray MLPutInteger64Array
#else
#error "IGraph/M expected IGRAPH_INTEGER_SIZE to be either 32 or 64."
#endif

static_assert(sizeof(mlint) == sizeof(igraph_integer_t), "Mismatch in 'mlint' and 'igraph_integer_t' sizes.");


/************************
 **** Error checking ****
 ************************/

// check igraph's error codes and abort if necessary
inline void igCheck(igraph_error_t err) {
    if (err == IGRAPH_SUCCESS) return;
    std::ostringstream msg;
    msg << "igraph returned with error: " << igraph_strerror(err) << ".";
    throw mma::LibraryError(msg.str());
}


/*******************************
 **** RAII for igraph types ****
 *******************************/

inline igraph_vector_t igVectorView(mma::RealTensorRef t) {
    static double dummy = 0.0; // work around igraph not liking zero-length vectors with NULL pointers
    igraph_vector_t vec;
    mint len = t.length();
    igraph_vector_view(&vec, len == 0 ? &dummy : t.data(), len);
    return vec;
}

inline igraph_vector_int_t igIntVectorView(mma::IntTensorRef t) {
    static igraph_integer_t dummy = 0; // work around igraph not liking zero-length vectors with NULL pointers
    igraph_vector_int_t vec;
    mint len = t.length();
    igraph_vector_int_view(&vec, len == 0 ? &dummy : t.data(), len);
    return vec;
}


// RAII for igraph_vector_t
class igVector {
public:
    igraph_vector_t vec;

    igVector() { igraph_vector_init(&vec, 0); }

    igVector(igVector &&source) noexcept {
        vec = source.vec;
        source.vec.stor_begin = nullptr;
    }

    igVector(const igraph_vector_t *source) { igraph_vector_init_copy(&vec, source); }

    explicit igVector(igraph_integer_t len) { igraph_vector_init(&vec, len); }

    igVector(const igVector &igv) : igVector() { igraph_vector_init_copy(&vec, &igv.vec); }

    igVector & operator = (const igVector &igv) {
        igraph_vector_update(&vec, &igv.vec);
        return *this;
    }

    // it is safe to call igraph_vector_destroy on a vector where vec.stor_begin == NULL
    ~igVector() { igraph_vector_destroy(&vec); }

    igraph_integer_t length() const { return vec.end - vec.stor_begin; }
    igraph_integer_t size() const { return length(); }

    igraph_real_t *begin() { return vec.stor_begin; }
    igraph_real_t *end() { return vec.end; }

    const igraph_real_t *begin() const { return vec.stor_begin; }
    const igraph_real_t *end() const { return vec.end; }

    igraph_real_t & operator [] (size_t i) { return begin()[i]; }
    const igraph_real_t & operator [] (size_t i) const { return begin()[i]; }

    void clear() { igraph_vector_clear(&vec); }
    void resize(igraph_integer_t newsize) { igCheck(igraph_vector_resize(&vec, newsize)); }
    void reserve(igraph_integer_t newsize) { igCheck(igraph_vector_reserve(&vec, newsize)); }

    void push_back(igraph_real_t el) { igCheck(igraph_vector_push_back(&vec, el)); }

    void copyFromMTensor(mma::RealTensorRef t) {
        igraph_vector_t from = igVectorView(t);
        igCheck(igraph_vector_update(&vec, &from));
    }

    mma::RealTensorRef makeMTensor() const { return mma::makeVector<double>(length(), begin()); }
};


// RAII for igraph_vector_int_t
// note that igraph_integer_t and mint may not be the same type
/* TODO move constructor like igVector */
class igIntVector {

    // avoid accidental implicit copy
    igIntVector(const igIntVector &) = delete;

public:
    igraph_vector_int_t vec;

    igIntVector() { igraph_vector_int_init(&vec, 0); }
    ~igIntVector() { igraph_vector_int_destroy(&vec); }

    igIntVector(igIntVector &&source) noexcept {
        vec = source.vec;
        source.vec.stor_begin = nullptr;
    }

    igIntVector(const igraph_vector_int_t *source) { igraph_vector_int_init_copy(&vec, source); }

    igIntVector & operator = (const igIntVector &igv) {
        igraph_vector_int_update(&vec, &igv.vec);
        return *this;
    }

    explicit igIntVector(igraph_integer_t len) { igraph_vector_int_init(&vec, len); }

    igraph_integer_t length() const { return vec.end - vec.stor_begin; }
    igraph_integer_t size() const { return length(); }

    igraph_integer_t *begin() { return vec.stor_begin; }
    igraph_integer_t *end() { return vec.end; }

    const igraph_integer_t *begin() const { return vec.stor_begin; }
    const igraph_integer_t *end() const { return vec.end; }

    igraph_integer_t & operator [] (size_t i) { return begin()[i]; }
    const igraph_integer_t & operator [] (size_t i) const { return begin()[i]; }

    void clear() { igraph_vector_int_clear(&vec); }    
    void resize(igraph_integer_t newsize) { igCheck(igraph_vector_int_resize(&vec, newsize)); }
    void reserve(igraph_integer_t newsize) { igCheck(igraph_vector_int_reserve(&vec, newsize)); }

    void push_back(igraph_real_t el) { igCheck(igraph_vector_int_push_back(&vec, el)); }

    void copyFromMTensor(mma::IntTensorRef t) {
        resize(t.length());
        std::copy(t.begin(), t.end(), begin());
    }

    mma::IntTensorRef makeMTensor() const { return mma::makeVector<mint>(length(), begin()); }
};


// RAII for igraph_vector_bool_t
class igBoolVector {

    // avoid accidental implicit copy
    igBoolVector(const igBoolVector &) = delete;
    igBoolVector & operator = (const igBoolVector &) = delete;

public:

    igraph_vector_bool_t vec;

    explicit igBoolVector(igraph_integer_t len) { igraph_vector_bool_init(&vec, len); }
    igBoolVector() { igraph_vector_bool_init(&vec, 0); }
    ~igBoolVector() { igraph_vector_bool_destroy(&vec); }

    igraph_integer_t length() const { return vec.end - vec.stor_begin; }
    igraph_integer_t size() const { return length(); }

    igraph_bool_t *begin() { return vec.stor_begin; }
    igraph_bool_t *end() { return vec.end; }

    const igraph_bool_t *begin() const { return vec.stor_begin; }
    const igraph_bool_t *end() const { return vec.end; }

    igraph_bool_t & operator [] (size_t i) { return begin()[i]; }
    const igraph_bool_t & operator [] (size_t i) const { return begin()[i]; }

    void clear() { igraph_vector_bool_clear(&vec); }
    void resize(igraph_integer_t newsize) { igCheck(igraph_vector_bool_resize(&vec, newsize)); }
    void reserve(igraph_integer_t newsize) { igCheck(igraph_vector_bool_reserve(&vec, newsize)); }

    void push_back(igraph_real_t el) { igCheck(igraph_vector_bool_push_back(&vec, el)); }

    mma::IntTensorRef makeMTensor() const { return mma::makeVector<mint>(length(), begin()); }
};


// RAII for igraph_maxtrix_t
class igMatrix {

    // avoid accidental implicit copy
    igMatrix(const igMatrix &) = delete;
    igMatrix & operator = (const igMatrix &) = delete;

public:

    igraph_matrix_t mat;

    igMatrix() { igraph_matrix_init(&mat, 0, 0); }
    ~igMatrix() { igraph_matrix_destroy(&mat); }

    igraph_integer_t length() const { return mat.data.end - mat.data.stor_begin; }
    igraph_integer_t nrow() const { return mat.nrow; }
    igraph_integer_t ncol() const { return mat.ncol; }

    igraph_real_t *begin() { return mat.data.stor_begin; }
    igraph_real_t *end() { return mat.data.end; }

    const igraph_real_t *begin() const { return mat.data.stor_begin; }
    const igraph_real_t *end() const { return mat.data.end; }

    igraph_real_t & operator () (size_t i, size_t j) { return MATRIX(mat, i, j); }
    const igraph_real_t & operator () (size_t i, size_t j) const { return MATRIX(mat, i, j); }

    void copyFromMTensor(mma::RealMatrixRef t) {
        igraph_vector_t from = igVectorView(t);
        igCheck(igraph_vector_update(&mat.data, &from));
        // Mathematica uses row-major storage, igraph uses column-major storage
        // thus we need to reverse the row/column counts, then transpose
        mat.nrow = t.cols();
        mat.ncol = t.rows();
        igCheck(igraph_matrix_transpose(&mat));
    }

    mma::RealMatrixRef makeMTensor() const { return mma::makeMatrixTransposed<double>(mat.nrow, mat.ncol, begin()); }
};


// RAII for igraph_maxtrix_t
class igIntMatrix {

    // avoid accidental implicit copy
    igIntMatrix(const igIntMatrix &) = delete;
    igIntMatrix & operator = (const igIntMatrix &) = delete;

public:

    igraph_matrix_int_t mat;

    igIntMatrix() { igraph_matrix_int_init(&mat, 0, 0); }
    ~igIntMatrix() { igraph_matrix_int_destroy(&mat); }

    igraph_integer_t length() const { return mat.data.end - mat.data.stor_begin; }
    igraph_integer_t nrow() const { return mat.nrow; }
    igraph_integer_t ncol() const { return mat.ncol; }

    igraph_integer_t *begin() { return mat.data.stor_begin; }
    igraph_integer_t *end() { return mat.data.end; }

    const igraph_integer_t *begin() const { return mat.data.stor_begin; }
    const igraph_integer_t *end() const { return mat.data.end; }

    igraph_integer_t & operator () (size_t i, size_t j) { return MATRIX(mat, i, j); }
    const igraph_integer_t & operator () (size_t i, size_t j) const { return MATRIX(mat, i, j); }

    void copyFromMTensor(mma::IntMatrixRef t) {
        igraph_vector_int_t from = igIntVectorView(t);
        igCheck(igraph_vector_int_update(&mat.data, &from));
        // Mathematica uses row-major storage, igraph uses column-major storage
        // thus we need to reverse the row/column counts, then transpose
        mat.nrow = t.cols();
        mat.ncol = t.rows();
        igCheck(igraph_matrix_int_transpose(&mat));
    }

    mma::RealMatrixRef makeMTensor() const { return mma::makeMatrixTransposed<double>(mat.nrow, mat.ncol, begin()); }
};


// RAII for igraph_vector_ptr_t containing igraph_vector_t
template<typename ElemType, void DestroyElem(ElemType *)>
class igPtrVector {

    void destroy_items() {
        for (void **ptr = list.stor_begin; ptr < list.end; ++ptr)
            if (*ptr)
                DestroyElem(reinterpret_cast<ElemType *>(*ptr));
    }

public:
    igraph_vector_ptr_t list;

    igPtrVector() {
        igraph_vector_ptr_init(&list, 0);
    }
    ~igPtrVector() {
        // we destroy items manually ...
        destroy_items();

        // ... and avoid calling any items destructors that may have been set
        igraph_vector_ptr_set_item_destructor(&list, nullptr);
        igraph_vector_ptr_free_all(&list);
        igraph_vector_ptr_destroy(&list);
    }

    void clear() {
        // this mirrors igraph_vector_ptr_clear(), but does not call
        // any item desctructors automatically to avoid double-free
        destroy_items();
        list.end = list.stor_begin;
    }

    igraph_integer_t length() const { return list.end - list.stor_begin; }
    igraph_integer_t size() const { return length(); }

    void push(igraph_vector_t *vec) { igraph_vector_ptr_push_back(&list, vec); }

    ElemType **begin() { return reinterpret_cast<ElemType **>(list.stor_begin); }
    ElemType **end() { return reinterpret_cast<ElemType **>(list.end); }

    const ElemType **begin() const { return reinterpret_cast<const ElemType **>(list.stor_begin); }
    const ElemType **end() const { return reinterpret_cast<const ElemType **>(list.end); }

    const ElemType *operator [] (igraph_integer_t i) const { return static_cast<const ElemType *>(list.stor_begin[i]); }
};


class igGraphList {
public:
    igraph_graph_list_t list;

    igGraphList() { igraph_graph_list_init(&list, 0); }
    ~igGraphList() { igraph_graph_list_destroy(&list); }

    void clear() {
        igraph_graph_list_clear(&list);
    }

    igraph_integer_t length() const { return list.end - list.stor_begin; }
    igraph_integer_t size() const { return length(); }

    const igraph_t *operator [] (igraph_integer_t i) const {
        /* TODO safety! use API? */
        return &list.stor_begin[i];
    }
};


/* TODO rename to igIntVectorList */
class igList {
public:
    igraph_vector_int_list_t list;

    igList() { igraph_vector_int_list_init(&list, 0); }
    ~igList() { igraph_vector_int_list_destroy(&list); }

    void clear() {
        igraph_vector_int_list_clear(&list);
    }

    void reserve(igraph_integer_t capacity) {
        /* TODO check error */
        igraph_vector_int_list_reserve(&list, capacity);
    }

    igraph_integer_t length() const { return list.end - list.stor_begin; }
    igraph_integer_t size() const { return length(); }

    const igraph_vector_int_t *operator [] (igraph_integer_t i) const {
        /* TODO safety! use API? */
        return &list.stor_begin[i];
    }
};

class igVectorList {
public:
    igraph_vector_list_t list;

    igVectorList() { igraph_vector_list_init(&list, 0); }
    ~igVectorList() { igraph_vector_list_destroy(&list); }

    void clear() {
        igraph_vector_list_clear(&list);
    }

    igraph_integer_t length() const { return list.end - list.stor_begin; }
    igraph_integer_t size() const { return length(); }

    const igraph_vector_t *operator [] (igraph_integer_t i) const {
        /* TODO safety! use API? */
        return &list.stor_begin[i];
    }
};


/* Generic wrapper for existing instances of igraph vector types */

template<typename Vec>
class igWrapper {
    Vec &vec;

public:
    constexpr igWrapper(Vec &vec) : vec(vec) { }

    decltype (Vec::stor_begin) begin() { return vec.stor_begin; }
    decltype (Vec::end) end() { return vec.end; }

    constexpr const decltype (Vec::stor_begin) begin() const { return vec.stor_begin; }
    constexpr const decltype (Vec::end) end() const { return vec.end; }

    constexpr const decltype (Vec::stor_begin) cbegin() const { return vec.stor_begin; }
    constexpr const decltype (Vec::end) cend() const { return vec.end; }

    constexpr size_t size() const { return end() - begin(); }
};

// returns an igWrapper object, which works with range-based for loops
template<typename Vec>
constexpr igWrapper<Vec> igWrap(Vec &vec) { return igWrapper<Vec>(vec); }


/********************************************
 **** RAII for LTemplate data structures ****
 ********************************************/

// A guard class to automatically free an LTemplate / LibraryLink data structure
template<typename Tensor>
class LTGuard {
    Tensor &t;

    bool active = true;

public:
    LTGuard(Tensor &t) : t(t) { }
    ~LTGuard() { if (active) t.free(); }

    Tensor &tensor() { return t; }

    void activate() { active = true; }
    void deactivate() { active = false; }
};


/****************************************************
 **** Extend mlstream with igraph-specific types ****
 ****************************************************/

inline mlStream & operator << (mlStream &ml, const igraph_vector_t &vec) {
    if (! MLPutReal64List(ml.link(), vec.stor_begin, vec.end - vec.stor_begin))
        ml.error("cannot return vector");
    return ml;
}


inline mlStream & operator << (mlStream &ml, const igVector &vec) { return ml << vec.vec; }


inline mlStream & operator << (mlStream &ml, const igraph_vector_int_t &vec) {
    if (! MLPutIGIntegerList(ml.link(), reinterpret_cast<mlint64 *>(vec.stor_begin), vec.end - vec.stor_begin))
        ml.error("cannot return integer vector");
    return ml;
}


inline mlStream & operator << (mlStream &ml, const igIntVector &vec) { return ml << vec.vec; }


inline mlStream & operator << (mlStream &ml, const igList &list) {
    igraph_integer_t len = list.length();
    if (! MLPutFunction(ml.link(), "List", len))
        ml.error("cannot return vector list");
    for (igraph_integer_t i=0; i < len; ++i)
        ml << VECTOR(list.list)[i]; /* TODO safe API */
    return ml;
}

inline mlStream & operator << (mlStream &ml, const igVectorList &list) {
    igraph_integer_t len = list.length();
    if (! MLPutFunction(ml.link(), "List", len))
        ml.error("cannot return vector list");
    for (igraph_integer_t i=0; i < len; ++i)
        ml << VECTOR(list.list)[i]; /* TODO safe API */
    return ml;
}

/* TODO protect against larger-than-int array sizes */
inline mlStream & operator << (mlStream &ml, const igMatrix &mat) {
    int ok;

    int dims[2];
    dims[0] = mat.ncol();
    dims[1] = mat.nrow();
    if (mat.nrow() == 0) {
        // 0-column matrices are represented as {{}, ... {}} in Mathematica.
        // However, 0-row matrices can only be represented as {}, which
        // loses the column-count information, and cannot be Transpose[]d
        // to a 0-column matrix. Therefore, we do the transposition
        // manually by swapping dimensions (which is safe since there are no
        // elements in the matrix).
        std::swap(dims[0], dims[1]);
        ok = MLPutReal64Array(ml.link(), mat.begin(), dims, nullptr, 2);
    } else {
        ok = MLPutFunction(ml.link(), "Transpose", 1) &&
             MLPutReal64Array(ml.link(), mat.begin(), dims, nullptr, 2);
    }
    if (! ok)
        ml.error("cannot return matrix");
    return ml;
}

/* TODO protect against larger-than-int array sizes */
inline mlStream & operator << (mlStream &ml, const igIntMatrix &mat) {
    int ok;

    int dims[2];
    dims[0] = mat.ncol();
    dims[1] = mat.nrow();
    if (mat.nrow() == 0) {
        // 0-column matrices are represented as {{}, ... {}} in Mathematica.
        // However, 0-row matrices can only be represented as {}, which
        // loses the column-count information, and cannot be Transpose[]d
        // to a 0-column matrix. Therefore, we do the transposition
        // manually by swapping dimensions (which is safe since there are no
        // elements in the matrix).
        std::swap(dims[0], dims[1]);
        ok = MLPutIGIntegerArray(ml.link(), reinterpret_cast<const mlint64 *>(mat.begin()), dims, nullptr, 2);
    } else {
        ok = MLPutFunction(ml.link(), "Transpose", 1) &&
             MLPutIGIntegerArray(ml.link(), reinterpret_cast<const mlint64 *>(mat.begin()), dims, nullptr, 2);
    }
    if (! ok)
        ml.error("cannot return matrix");
    return ml;
}


inline mlStream & operator >> (mlStream &ml, igList &list) {
    int len;
    if (! MLTestHead(ml.link(), "List", &len))
        ml.error("List of lists expected");
    list.clear();
    list.reserve(len);
    for (int i=0; i < len; ++i) {
        igraph_vector_int_t vec;
        mlint *data;
        int listlen;
        if (! MLGetIGIntegerList(ml.link(), &data, &listlen)) {
            // raise error
            ml.error("Real64List expected in list of lists");
        }
        igraph_vector_int_init_array(&vec, reinterpret_cast<igraph_integer_t *>(data), listlen); // TODO check success
        igraph_vector_int_list_push_back(&list.list, &vec);
        MLReleaseIGIntegerList(ml.link(), data, listlen);
        VECTOR(list.list)[i] = vec;
    }
    return ml;
}


inline mlStream & operator >> (mlStream &ml, igVector &vec) {
    double *data;
    int length;
    if (! MLGetReal64List(ml.link(), &data, &length))
        ml.error("Real64List expected");
    vec.resize(length);
    std::copy(data, data+length, vec.begin());
    MLReleaseReal64List(ml.link(), data, length);
    return ml;
}


inline mlStream & operator >> (mlStream &ml, igMatrix &mat) {
    double *data;
    int *dims;
    char **heads;
    int depth;
    if (! MLGetReal64Array(ml.link(), &data, &dims, &heads, &depth))
        ml.error("Real64 matrix expected");
    if (depth != 2)
        ml.error("Real64 matrix expected, depth doesn't match");

    int length = 1;
    for (int i=0; i < depth; ++i)
        length *= dims[i];
    igraph_vector_resize(&mat.mat.data, length);
    std::copy(data, data+length, mat.begin());
    mat.mat.nrow = dims[1];
    mat.mat.ncol = dims[0];
    MLReleaseReal64Array(ml.link(), data, dims, heads, depth);
    return ml;
}


inline mlStream & operator >> (mlStream &ml, igIntMatrix &mat) {
    mlint *data;
    int *dims;
    char **heads;
    int depth;
    if (! MLGetIGIntegerArray(ml.link(), &data, &dims, &heads, &depth))
        ml.error("Integer matrix expected");
    if (depth != 2)
        ml.error("Integer matrix expected, depth doesn't match");

    int length = 1;
    for (int i=0; i < depth; ++i)
        length *= dims[i];
    igraph_vector_int_resize(&mat.mat.data, length);
    std::copy(data, data+length, mat.begin());
    mat.mat.nrow = dims[1];
    mat.mat.ncol = dims[0];
    MLReleaseIGIntegerArray(ml.link(), data, dims, heads, depth);
    return ml;
}


inline mlStream & operator >> (mlStream &ml, igIntVector &vec) {
    mlint *data;
    int length;
    if (! MLGetIGIntegerList(ml.link(), &data, &length))
        ml.error("Integer32List expected");
    vec.resize(length);
    std::copy(data, data+length, vec.begin());
    MLReleaseIGIntegerList(ml.link(), data, length);
    return ml;
}


inline mlStream & operator >> (mlStream &ml, igBoolVector &vec) {
    unsigned char *data;
    int length;
    // igraph_bool_t is a bool, so we use the smallest MathLink type: Integer8
    if (! MLGetInteger8List(ml.link(), &data, &length))
        ml.error("Integer32List expected");
    vec.resize(length);
    std::copy(data, data+length, vec.begin());
    MLReleaseInteger8List(ml.link(), data, length);
    return ml;
}


inline mlStream & operator >> (mlStream &ml, igraph_bool_t &b) {
    unsigned char c;
    MLGetInteger8(ml.link(), &c);
    b = c;
    return ml;
}


/*********************************
 **** Other utility functions ****
 *********************************/

// packs an igList (usually representing vertex or edge sets) into
// a single IntTensor for fast transfer
inline mma::IntTensorRef packListIntoIntTensor(const igList &list) {
    std::vector<igraph_integer_t> lengths;
    igraph_integer_t list_length = list.length();
    mint total_length = 0;
    for (igraph_integer_t i=0; i < list_length; ++i) {
        igraph_integer_t len = igraph_vector_int_size(&VECTOR(list.list)[i]); /* TODO proper API */
        total_length += len;
        total_length += 1;
        lengths.push_back(len);
    }
    total_length += 1;
    mma::IntTensorRef t = mma::makeVector<mint>(total_length);
    t[0] = list_length;
    std::copy(lengths.begin(), lengths.end(), t.begin() + 1);
    mint *ptr = t.begin() + 1 + list_length;
    for (igraph_integer_t i=0; i < list_length; ++i) {
        igraph_integer_t *b = &VECTOR(VECTOR(list.list)[i])[0]; /* TODO proper API */
        std::copy(b, b+lengths[i], ptr);
        ptr += lengths[i];
    }
    return t;
}


inline igraph_neimode_t igNeighborMode(mint mode, const char *context) {
    switch (mode) {
    case 1: return IGRAPH_OUT;
    case 2: return IGRAPH_IN;
    case 3: return IGRAPH_ALL;
    default: throw mma::LibraryError(std::string("Invalid neighbour mode for ") + context + ".");
    }
}


#endif // IG_COMMON_H
