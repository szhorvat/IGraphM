/*
 * Copyright (c) 2016-2020 Szabolcs Horvát.
 *
 * See the file LICENSE.txt for copying permission.
 */

#include "IGlobal.h"
#include <cstring>
#include <csetjmp>

/* Detecting ASan with GCC:
 *   https://gcc.gnu.org/onlinedocs/cpp/Common-Predefined-Macros.html
 * Detecting ASan with Clang:
 *   https://clang.llvm.org/docs/AddressSanitizer.html#conditional-compilation-with-has-feature-address-sanitizer
 */
#if defined(__SANITIZE_ADDRESS__)
#  define IG_SANITIZER_AVAILABLE 1
#elif defined(__has_feature)
#  if __has_feature(address_sanitizer)
#    define IG_SANITIZER_AVAILABLE 1
#  endif
#endif

#ifdef IG_SANITIZER_AVAILABLE
#include <cstdio>
#include <sanitizer/asan_interface.h>
#endif

static inline bool is_punctuated(const char *str) {
    const size_t len = strlen(str);
    if (len == 0)
        return true;
    if (str[len-1] == '.' || str[len-1] == '!' || str[len-1] == '?')
        return true;
    return false;
}

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
    // If the message is not puncutated, add a fullstop.
    if (! is_punctuated(reason))
        msg << '.';
    mma::message(msg.str(), mma::M_WARNING);
}


void igErrorHandler(const char *reason, const char *file, int line, int /* igraph_errno */) {
    // avoid printing empty messages
    if (strlen(reason) != 0) {
        std::ostringstream msg;
        msg << file << ":" << line << " - " << reason;
        // If the message is not puncutated, add a fullstop.
        if (! is_punctuated(reason))
            msg << '.';
        mma::message(msg.str(), mma::M_ERROR);
    }
    IGRAPH_FINALLY_FREE();
}


void igFatalHandler(const char *reason, const char *file, int line) {

#ifdef IG_SANITIZER_AVAILABLE
    // When using AddressSanitizer, also print the fatal error to the terminal,
    // and show a stack trace.
    fprintf(stderr, "\nFatal igraph error at %s:%i : %s\n", file, line, reason);
    fprintf(stderr, "Stack trace:\n");
    __sanitizer_print_stack_trace();
#endif

    // Use a separate block to ensure that the destructor of msg is run before the mma::fatal_error()
    // mma::fatal_error() does a longjmp(), which may prevent the destructor from being called.
    {
        std::ostringstream msg;
        msg << "Unexpected error in igraph. Please report this as a bug, along with the steps to reproduce it.\n";
        msg << file << ":" << line << " - " << reason;
        // If the message is not puncutated, add a fullstop.
        if (! is_punctuated(reason))
            msg << '.';
        mma::message(msg.str(), mma::M_ERROR);
    }

    IGRAPH_FINALLY_CLEAN(IGRAPH_FINALLY_STACK_SIZE());
    mma::fatal_error();
}


double progress_reporting_granularity = 1.0;

int igProgressHandler(const char *message, igraph_real_t percent, void *data) {
    static igraph_real_t last_percent = 100;

    // Only report progress if it has increased by at least the given granularity,
    // if it has reached 100%, or if it has decreased since the last report.
    if (percent < last_percent || percent >= 100 || percent >= last_percent + progress_reporting_granularity) {
        MLINK link = mma::libData->getMathLink(mma::libData);
        MLPutFunction(link, "EvaluatePacket", 1);
            MLPutFunction(link, "Set", 2);
                MLPutFunction(link, "List", 2);
                    MLPutSymbol(link, "IGraphM`Progress`Message");
                    MLPutSymbol(link, "IGraphM`Progress`Percent");
                MLPutFunction(link, "List", 2);
                    MLPutString(link, message);
                    MLPutReal64(link, percent);
        mma::libData->processMathLink(link);
        int pkt = MLNextPacket(link);
        if (pkt == RETURNPKT)
            MLNewPacket(link);

        last_percent = percent;
    }

    return IGRAPH_SUCCESS;
}


/***** Read Graph6, Digraph6 and Sparse6 *****/

// MAXBYTE conflicts with winnt.h on Windows
#undef MAXBYTE

constexpr int BIAS6 = 63;
constexpr int SMALLN = 62;
constexpr int MAXBYTE = 126;
constexpr int SMALLISHN = 258047;

/* From gtools.c in nauty: */
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
    size_t len;
    {
        // scan for illegal characters and determine string length
        const char *p = str;
        if (*p == ':' || *p == ';' || *p == '&')
            p++;
        while (*p >= BIAS6 && *p <= MAXBYTE)
            p++;
        if (*p != '\0')
            throw mma::LibraryError("Illegal character in Graph6, Digraph6 or Sparse6 line.");
        len = p-str;
    }

    // The result is returned as an integer vector with the following elements:
    // 1st value: directed?
    // 2nd value: vertex count
    // Rest: edges as vertex pairs
    std::vector<mint> result;

    if (*str == ':') // non-incremental sparse6:
    {
        str++; // skip ':' character
        if (*str == '\0')
            throw mma::LibraryError("Truncated Sparse6 line.");
        int n = graphsize(str);

        result.push_back(0); // undirected
        result.push_back(n); // vertex count       

        int k; // no. of bits needed to represent n-1
        if (n > 1) {
            k  = sizeof(int);
            while (!( (n-1) & (1 << (k-1)) ))
                k--;
        } else {
            k = 1;
        }

        mint vertex = 0;

        const char *p = str + SIZELEN(n);

        char val;    // current value
        char offset; // bits still unread in the current value

        if (n == 0) goto s6done;

#define S6_GETCHAR() \
    if (*p == '\0') goto s6done; \
    val = *p++ - BIAS6; \
    offset = 6;

        offset = 0;
        while (true) {
            if (offset == 0) {
                S6_GETCHAR();
            }

            bool b = val & (1 << (offset-1));

            offset -= 1;
            if (offset == 0) {
                S6_GETCHAR();
            }

            int x = 0;
            int remaining = k; // bits to still read to complete x
            while (remaining > 0) {
                if (remaining > offset) {
                    x = x << offset;
                    x += val & ((1 << offset) - 1);
                    remaining -= offset;
                    S6_GETCHAR();
                } else {
                    x = x << remaining;
                    x += (val >> (offset - remaining)) & ((1 << remaining) - 1);
                    offset -= remaining;
                    remaining = 0;
                }
            }

            if (b) vertex++;
            if (x > vertex) {
                vertex = x;
            } else if (vertex < n) {
                result.push_back(x);
                result.push_back(vertex);
            }

        }

#undef S6_GETCHAR

    s6done:
        return mma::makeVector<mint>(result.size(), result.data());
    }
    else if (*str == ';') // incremental sparse6:
    {
        throw mma::LibraryError("Incremental Sparse6 is not implemented.");
    }
    else if (*str == '&') // digraph6:
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

        if (n <= 1) goto d6done;

        while (true) {
            unsigned char byte = *p++ - BIAS6;
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
    else if (*str < BIAS6 || *str > MAXBYTE) { // invalid character for graph6 data
        throw mma::LibraryError("Invalid Graph6 data.");
    }
    else // graph6:
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

        if (n <= 1) goto g6done;

        while (true) {
            unsigned char byte = *p++ - BIAS6;
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
