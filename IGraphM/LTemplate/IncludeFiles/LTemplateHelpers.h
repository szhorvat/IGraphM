#ifndef LTEMPLATE_HELPERS_H
#define LTEMPLATE_HELPERS_H

#include "LTemplate.h"
#include <algorithm>

namespace mma {
namespace detail {

// Functions for getting and setting arguments and return values

template<typename T>
inline TensorRef<T> getTensor(MArgument marg) { return TensorRef<T>(MArgument_getMTensor(marg)); }

template<typename T>
inline void setTensor(MArgument marg, TensorRef<T> &val) { MArgument_setMTensor(marg, val.tensor()); }

template<typename T>
inline SparseArrayRef<T> getSparseArray(MArgument marg) { return MArgument_getMSparseArray(marg); }

template<typename T>
inline void setSparseArray(MArgument marg, SparseArrayRef<T> &val) { MArgument_setMSparseArray(marg, val.sparseArray()); }

inline complex_t getComplex(MArgument marg) {
    mcomplex c = MArgument_getComplex(marg);
    return complex_t(c.ri[0], c.ri[1]);
}

inline void setComplex(MArgument marg, complex_t val) {
    mcomplex *c = reinterpret_cast<mcomplex *>(&val);
    MArgument_setComplex(marg, *c);
}


inline const char *getString(MArgument marg) {
    return const_cast<const char *>(MArgument_getUTF8String(marg));
}

inline void setString(MArgument marg, const char *val) {
    MArgument_setUTF8String(marg, const_cast<char *>(val));
}


template<typename Collection>
inline IntTensorRef get_collection(const Collection &collection) {
    IntTensorRef ids = makeVector<mint>(collection.size());

    typename Collection::const_iterator i = collection.begin();
    mint *j = ids.begin();
    for (; i != collection.end(); ++i, ++j)
        *j = i->first;

    return ids;
}


template<typename T>
class getObject {
    std::map<mint, T *> &collection;
public:
    explicit getObject(std::map<mint, T *> &coll) : collection(coll) { }
    T & operator () (MArgument marg) { return *(collection[MArgument_getInteger(marg)]); }
};


} // namespace detail
} // namespace mma

#endif // LTEMPLATE_HELPERS_H
