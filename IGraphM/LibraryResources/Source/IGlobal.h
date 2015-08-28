#ifndef IGLOBAL_H
#define IGLOBAL_H

extern "C" {
#include <igraph/igraph.h>
}

#include "LTemplate.h"

#include <random>

class IGlobal {
public:
    void init() {
        igraph_set_error_handler(igraph_error_handler_printignore);
        igraph_rng_seed(igraph_rng_default(), std::random_device()());
    }

    ~IGlobal() { }

    const char *version() {
        const char *ver;
        igraph_version(&ver, NULL, NULL, NULL);
        return ver;
    }

    void seedRandom(mint s) {
        igraph_rng_seed(igraph_rng_default(), s);
    }
};

#endif // IGLOBAL_H

