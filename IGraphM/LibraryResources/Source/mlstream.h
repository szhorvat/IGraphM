#ifndef MLSTREAM_H
#define MLSTREAM_H

#include "LTemplate.h"

#include <list>
#include <string>
#include <sstream>

class mlStream {
    MLINK lp;
    std::string context;

public:
    mlStream(MLINK lp_) : lp(lp_) { }
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


struct mlCheckArgs {
    int argc;

    mlCheckArgs(int argc_) : argc(argc_) { }
};

struct mlHead {
    const char *head;
    int argc;

    mlHead(const char *head_, int argc_) : head(head_), argc(argc_) { }
};

inline mlStream & operator >> (mlStream &ml, short &x) {
    if (! MLGetInteger16(ml.link(), &x))
        ml.error("Integer32 expected");
    return ml;
}

inline mlStream & operator >> (mlStream &ml, int &x) {
    if (! MLGetInteger32(ml.link(), &x))
        ml.error("Integer32 expected");
    return ml;
}

inline mlStream & operator >> (mlStream &ml, mlint64 &x) {
    if (! MLGetInteger64(ml.link(), &x))
        ml.error("Integer32 expected");
    return ml;
}

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
        ml.error("Cannot return string.");
    return ml;
}

inline mlStream & operator << (mlStream &ml, short x) {
    if (! MLPutInteger16(ml.link(), x))
        ml.error("Cannot return Integer16.");
    return ml;
}

inline mlStream & operator << (mlStream &ml, int x) {
    if (! MLPutInteger32(ml.link(), x))
        ml.error("Cannot return Integer32.");
    return ml;
}

inline mlStream & operator << (mlStream &ml, long x) {
    if (! MLPutInteger64(ml.link(), x))
        ml.error("Cannot return Integer64.");
    return ml;
}

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

inline mlStream & operator << (mlStream &ml, const mlHead &head) {
    if (! MLPutFunction(ml.link(), head.head, head.argc)) {
        std::ostringstream msg;
        msg << "Cannot put head " << head.head << " with " << head.argc << " arguments";
        ml.error(msg.str());
    }
    return ml;
}

template<typename T>
inline mlStream & operator << (mlStream &ml, const std::list<T> &ls) {
    ml << mlHead("List", ls.size());
    for (typename std::list<T>::const_iterator i = ls.begin(); i != ls.end(); ++i)
        ml << *i;
    return ml;
}


#endif // MLSTREAM_H

