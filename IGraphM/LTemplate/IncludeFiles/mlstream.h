#ifndef MLSTREAM_H
#define MLSTREAM_H

/** \file mlstream.h
 * mlstream.h is an auxiliary header for LTemplate to ease reading arguments and returning values through MathLink.
 *
 * LTemplate itself does not depend on mlstream.h, so if you don't use this header,
 * feel free to remove it from your project.  mlstream.h is not meant as a general
 * MathLink interface.  It is specifically designed for handling arguments and return
 * values in conjunction with LTemplate and LinkObject based functions.
 *
 * Example usage:
 * \code
 * void addMult(MLINK link) {
 *   mlStream ml(link, "addMult"); // any errors messages will mention the context "addMult"
 *
 *   int i, j;
 *   ml >> mlCheckArgs(2) // read off the head List and check argument count
 *      >> i >> j;        // read two integer arguments
 *
 *   // compute the result
 *   int sum = i+j;
 *   int prod = i*j;
 *
 *   // alias for MLNewPacket, must be used before returning the result
 *   ml.newPacket();
 *
 *   // we return two results in a list
 *   ml << mlHead("List", 2)
 *      << sum << prod;
 * }
 * \endcode
 */

#include "LTemplate.h"

#include <vector>
#include <list>
#include <string>
#include <sstream>

/// Wrapper for MLINK to allow using extractors and inserters
class mlStream {
    MLINK lp;
    std::string context;

public:
    explicit mlStream(MLINK lp_) : lp(lp_) { }
    mlStream(MLINK lp_, std::string context_) : lp(lp_), context(context_) { }

    MLINK link() { return lp; }

    void error(const std::string err) {
        std::ostringstream msg;
        if (! context.empty())
            msg << context << ": ";
        msg << err << ".";
        throw mma::LibraryError(msg.str());
    }

    void newPacket() {
        MLNewPacket(lp);
    }
};


// Special

/// Must be the first item extracted from mlStream, checks number of arguments.
struct mlCheckArgs {
    int argc;

    explicit mlCheckArgs(int argc_) : argc(argc_) { }
};

inline mlStream & operator >> (mlStream &ml, const mlCheckArgs &ca) {
    int count;

    if (! MLTestHead(ml.link(), "List", &count))
        ml.error("argument check: head \"List\" expected");

    if (count != ca.argc){
        std::ostringstream msg;
        msg << ca.argc << " arguments expected, " << count << " received";
        ml.error(msg.str());
    }

    return ml;
}


/// Used for inserting a head with the given argument count into mlStream.
/** Tyically used with the head List when returning multiple results. */
struct mlHead {
    const char *head;
    int argc;

    mlHead(const char *head_, int argc_) : head(head_), argc(argc_) { }
};

inline mlStream & operator << (mlStream &ml, const mlHead &head) {
    if (! MLPutFunction(ml.link(), head.head, head.argc)) {
        std::ostringstream msg;
        msg << "Cannot put head " << head.head << " with " << head.argc << " arguments";
        ml.error(msg.str());
    }
    return ml;
}


/// Used for inserting a symbol into mlStream
struct mlSymbol {
    const char *symbol;

    explicit mlSymbol(const char *symbol_) : symbol(symbol_) { }
};

inline mlStream & operator << (mlStream &ml, const mlSymbol &symbol) {
    if (! MLPutSymbol(ml.link(), symbol.symbol)) {
        std::ostringstream msg;
        msg << "Cannot put symbol " << symbol.symbol;
        ml.error(msg.str());
    }
    return ml;
}


struct mlDiscard {
    const int count;
    explicit mlDiscard(int count_ = 1) : count(count_) {}
};

inline mlStream & operator >> (mlStream &ml, const mlDiscard &drop) {
    for (int i=0; i < drop.count; ++i)
        if (! MLTransferExpression(NULL, ml.link()))
            ml.error("Cannot discard expression");
    return ml;
}


// Basic types (integer and floating point)

#define MLSTREAM_DEF_BASIC_GET(MTYPE, CTYPE) \
    inline mlStream & operator >> (mlStream &ml, CTYPE &x) { \
        if (! MLGet ## MTYPE(ml.link(), &x)) \
            ml.error(#MTYPE " expected"); \
        return ml; \
    }

MLSTREAM_DEF_BASIC_GET(Integer16, short)
MLSTREAM_DEF_BASIC_GET(Integer32, int)
MLSTREAM_DEF_BASIC_GET(Integer64, mlint64)
MLSTREAM_DEF_BASIC_GET(Real32, float)
MLSTREAM_DEF_BASIC_GET(Real64, double)
MLSTREAM_DEF_BASIC_GET(Real128, mlextended_double)

#ifdef MLSTREAM_32BIT_INT_AND_LONG
inline mlStream & operator >> (mlStream &ml, long &x) { ml >> reinterpret_cast<int &>(x); }
#endif


#define MLSTREAM_DEF_BASIC_PUT(MTYPE, CTYPE) \
    inline mlStream & operator << (mlStream &ml, CTYPE x) { \
        if (! MLPut ## MTYPE(ml.link(), x)) \
            ml.error("Cannot return " #MTYPE); \
        return ml; \
    }

MLSTREAM_DEF_BASIC_PUT(Integer16, short)
MLSTREAM_DEF_BASIC_PUT(Integer32, int)
MLSTREAM_DEF_BASIC_PUT(Integer64, mlint64)
MLSTREAM_DEF_BASIC_PUT(Real32, float)
MLSTREAM_DEF_BASIC_PUT(Real64, double)
MLSTREAM_DEF_BASIC_PUT(Real128, mlextended_double)

#ifdef MLSTREAM_32BIT_INT_AND_LONG
inline mlStream & operator << (mlStream &ml, long x) { ml << static_cast<int>(x); }
#endif


// Strings

inline mlStream & operator >> (mlStream &ml, std::string &s) {
    const unsigned char *sp;
    int bytes, chars;
    if (! MLGetUTF8String(ml.link(), &sp, &bytes, &chars))
        ml.error("String expected");
    s.assign(reinterpret_cast<const char *>(sp), bytes);
    return ml;
}


inline mlStream & operator << (mlStream &ml, const std::string &s) {
    if (! MLPutUTF8String(ml.link(), reinterpret_cast<const unsigned char *>(s.c_str()), s.size()))
        ml.error("Cannot return UTF8 string");
    return ml;
}

inline mlStream & operator << (mlStream &ml, const char *s) {
    if (! MLPutString(ml.link(), s))
        ml.error("Cannot return string");
    return ml;
}

// TensorRef

#define MLSTREAM_DEF_TENSOR_PUT(MTYPE, CTYPE) \
    inline mlStream & operator << (mlStream &ml, mma::TensorRef<CTYPE> t) { \
        const int maxrank  = 16; \
        const int rank = t.rank(); \
        const mint *mdims = t.dimensions(); \
        int dims[maxrank]; \
        massert(rank <= maxrank); \
        std::copy(mdims, mdims + rank, dims); \
        if (! MLPut ## MTYPE ## Array(ml.link(), t.data(), dims, NULL, rank)) \
            ml.error("Cannot return " #CTYPE " tensor"); \
        return ml; \
    }

// we need to define this for both Integer32 and Integer64 as mint could be either
MLSTREAM_DEF_TENSOR_PUT(Integer32, int)
MLSTREAM_DEF_TENSOR_PUT(Integer64, mlint64)
MLSTREAM_DEF_TENSOR_PUT(Real64, double)

// TODO support complex

// Standard containers -- list

template<typename T>
inline mlStream & operator << (mlStream &ml, const std::list<T> &ls) {
    ml << mlHead("List", ls.size());
    for (typename std::list<T>::const_iterator i = ls.begin(); i != ls.end(); ++i)
        ml << *i;
    return ml;
}


// Standard containers -- vector

template<typename T>
inline mlStream & operator << (mlStream &ml, const std::vector<T> &vec) {
    ml << mlHead("List", vec.size());
    for (typename std::vector<T>::const_iterator i = vec.begin(); i != vec.end(); ++i)
        ml << *i;
    return ml;
}

#define MLSTREAM_DEF_VEC_PUT(MTYPE, CTYPE) \
    inline mlStream & operator << (mlStream &ml, const std::vector<CTYPE> &vec) { \
        if (! MLPut ## MTYPE ## List(ml.link(), &vec[0], vec.size())) \
            ml.error("Cannot return vector of " #MTYPE); \
        return ml; \
    }

MLSTREAM_DEF_VEC_PUT(Integer16, short)
MLSTREAM_DEF_VEC_PUT(Integer32, int)
MLSTREAM_DEF_VEC_PUT(Integer64, mlint64)
MLSTREAM_DEF_VEC_PUT(Real32, float)
MLSTREAM_DEF_VEC_PUT(Real64, double)
MLSTREAM_DEF_VEC_PUT(Real128, mlextended_double)

#ifdef MLSTREAM_32BIT_INT_AND_LONG
inline mlStream & operator << (mlStream &ml, const std::vector<long> &vec) {
    ml << reinterpret_cast<const std::vector<int> &>(vec);
}
#endif

#define MLSTREAM_DEF_VEC_GET(MTYPE, CTYPE) \
    inline mlStream & operator >> (mlStream &ml, std::vector<CTYPE> &vec) { \
        CTYPE *data; \
        int count; \
        if (! MLGet ## MTYPE ## List(ml.link(), &data, &count)) \
            ml.error(#MTYPE " list expected"); \
        vec.resize(count); \
        std::copy(data, data+count, vec.begin()); \
        MLRelease ## MTYPE ## List(ml.link(), data, count); \
        return ml; \
    }

MLSTREAM_DEF_VEC_GET(Integer16, short)
MLSTREAM_DEF_VEC_GET(Integer32, int)
MLSTREAM_DEF_VEC_GET(Integer64, mlint64)
MLSTREAM_DEF_VEC_GET(Real32, float)
MLSTREAM_DEF_VEC_GET(Real64, double)
MLSTREAM_DEF_VEC_GET(Real128, mlextended_double)

#ifdef MLSTREAM_32BIT_INT_AND_LONG
inline mlStream & operator >> (mlStream &ml, std::vector<long> &vec) {
  ml >> reinterpret_cast<std::vector<int> &>(vec);
}
#endif


#endif // MLSTREAM_H

