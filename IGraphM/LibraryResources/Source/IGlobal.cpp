/*
 * Copyright (c) 2019 Szabolcs Horvát.
 *
 * See the file LICENSE.txt for copying permission.
 */

#include "IGlobal.h"
#include <cstring>


int igInterruptionHandler(void *) {
    if (mma::libData->AbortQ()) {
        IGRAPH_FINALLY_FREE();
        return IGRAPH_INTERRUPTED;
    }
    return IGRAPH_SUCCESS;
}


void igWarningHandler(const char *reason, const char *file, int line, int /* igraph_errno */) {
    std::ostringstream msg;
    msg << file << ":" << line << " - " << reason;
    mma::message(msg.str(), mma::M_WARNING);
}


void igErrorHandler(const char *reason, const char *file, int line, int /* igraph_errno */) {
    // avoid printing empty messages
    if (strlen(reason) != 0) {
        std::ostringstream msg;
        msg << file << ":" << line << " - " << reason;
        mma::message(msg.str(), mma::M_ERROR);
    }
    IGRAPH_FINALLY_FREE();
}


/***** Read Graph6, Digraph6 and Sparse6 *****/

// MAXBYTE conflicts with winnt.h on Windows
#undef MAXBYTE

constexpr int BIAS6 = 63;
constexpr int SMALLN = 62;
constexpr int MAXBYTE = 126;
constexpr int SMALLISHN = 258047;

/* Get size of graph out of graph6, digraph6 or sparse6 string. */
int graphsize(const char *s) {
    const char *p;
    int n;

    if (s[0] == ':' || s[0] == '&') p = s+1;
    else                            p = s;
    n = *p++ - BIAS6;

    if (n > SMALLN)
    {
        n = *p++ - BIAS6;
        if (n > SMALLN)
        {
            n = *p++ - BIAS6;
            n = (n << 6) | (*p++ - BIAS6);
            n = (n << 6) | (*p++ - BIAS6);
            n = (n << 6) | (*p++ - BIAS6);
            n = (n << 6) | (*p++ - BIAS6);
            n = (n << 6) | (*p++ - BIAS6);
        }
        else
        {
            n = (n << 6) | (*p++ - BIAS6);
            n = (n << 6) | (*p++ - BIAS6);
        }
    }
    return n;
}


#define SIZELEN(n) ((n)<=SMALLN?1:((n)<=SMALLISHN?4:8))
        /* length of size code in bytes */
#define G6BODYLEN(n) \
   (((size_t)(n)/12)*((size_t)(n)-1) + (((size_t)(n)%12)*((size_t)(n)-1)+11)/12)
#define G6LEN(n) (SIZELEN(n) + G6BODYLEN(n))
  /* exact graph6 string length excluding \n\0
     This twisted expression works up to n=227023 in 32-bit arithmetic
     and for larger n if size_t has 64 bits.  */
#define D6BODYLEN(n) \
   ((n)*(size_t)((n)/6) + (((n)*(size_t)((n)%6)+5)/6))
#define D6LEN(n) (1 + SIZELEN(n) + D6BODYLEN(n))
  /* exact digraph6 string length excluding \n\0
     This twisted expression works up to n=160529 in 32-bit arithmetic
     and for larger n if size_t has 64 bits.  */


mma::IntTensorRef IGlobal::fromNauty(const char *str) {
    size_t len = strlen(str);

    // The result is returned as an integer vector with the following elements:
    // 1st value: directed?
    // 2nd value: vertex count
    // Rest: edges as vertex pairs
    std::vector<mint> result;

    if (*str == ':')
    {
        throw mma::LibraryError("Sparse6 is not yet implemented.");

        str++; // skip ':' character
        int n = graphsize(str);

        result.push_back(0); // undirected
        result.push_back(n); // vertex count

        int k = std::ceil(std::log2(n-1));

        mint vertex = 0;

        const char *p = str + SIZELEN(n);

        char val = *p - BIAS6;
        char offset = 6; // bits still unread in the current value

        bool b = val & (1 << (offset-1));

        offset -= 1;
        if (offset == 0) {
            val = *p++ - BIAS6;
            offset = 6;
        }

        int x = 0;
        int remaining = k; // bits to still read
        while (remaining > 0) {
            if (remaining > offset) {
                x = x << offset;
                x += val & ((1 << offset) - 1);
                remaining -= offset;
                val = *p++ - BIAS6;
                offset = 6;
            } else {
                x = x << remaining;
                x += (val >> (offset - remaining)) & ((1 << remaining) - 1);
                offset -= remaining;
            }
        }

        if (b) vertex++;
        if (x > vertex) {
            vertex = x;
        } else {
            result.push_back(x);
            result.push_back(vertex);
        }
    }
    else if (*str == ';')
    {
        throw mma::LibraryError("Incremental Sparse6 is not implemented.");
    }
    else if (*str == '&')
    {
        str++; // skip '&' character
        int n = graphsize(str);

        result.push_back(1); // directed
        result.push_back(n); // vertex count

        const int d6len = D6LEN(n);
        if (len < d6len)
            throw mma::LibraryError("Truncated Digraph6 line.");
        else if (len > d6len)
            throw mma::LibraryError("Unexpected characters at the end of Digraph6 line.");

        const char *p = str + SIZELEN(n);

        int i=0, j=0;
        while (true) {
            char byte = *p++ - BIAS6;
            for (int k=0; k < 6; ++k) {
                if (byte & 32) {
                    result.push_back(i);
                    result.push_back(j);
                }
                j++;
                if (j==n) {
                    j=0;
                    i++;
                }
                if (i == n)
                    goto d6done;
                byte <<= 1;
            }
        }
    d6done:

        return mma::makeVector<mint>(result.size(), result.data());
    }
    else if (*str < BIAS6 || *str > MAXBYTE) {
        throw mma::LibraryError("Invalid Graph6 data.");
    }
    else
    {
        int n = graphsize(str);

        result.push_back(0); // undirected
        result.push_back(n); // vertex count

        const int g6len = G6LEN(n);
        if (len < g6len)
            throw mma::LibraryError("Truncated Graph6 line.");
        else if (len > g6len)
            throw mma::LibraryError("Unexpected characters at the end of Graph6 line.");

        const char *p = str + SIZELEN(n);

        int i=0, j=1;

        while (true) {
            char byte = *p++ - BIAS6;
            for (int k=0; k < 6; ++k) {
                if (byte & 32) {
                    result.push_back(i);
                    result.push_back(j);
                }
                i++;
                if (i==j) {
                    i=0;
                    j++;
                }
                if (j == n)
                    goto g6done;
                byte <<= 1;
            }
        }
    g6done:

        return mma::makeVector<mint>(result.size(), result.data());
    }
}
