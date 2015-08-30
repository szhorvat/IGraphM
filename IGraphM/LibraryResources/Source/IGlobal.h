#ifndef IGLOBAL_H
#define IGLOBAL_H

#include "IGCommon.h"

#include <random>


void igWarningHandler(const char *reason, const char *file, int line, int igraph_errno);
void igErrorHandler(const char *reason, const char *file, int line, int igraph_errno);
int  igInterruptionHandler(void *);


class IGlobal {
public:
    void init() {
        igraph_set_error_handler(igErrorHandler);
        igraph_set_warning_handler(igWarningHandler);
        igraph_set_interruption_handler(igInterruptionHandler);
        igCheck(igraph_rng_seed(igraph_rng_default(), std::random_device()()));
    }

    ~IGlobal() { }

    const char *version() {
        const char *ver;
        igCheck(igraph_version(&ver, NULL, NULL, NULL));
        return ver;
    }

    void seedRandom(mint s) {
        igraph_rng_seed(igraph_rng_default(), s);
    }
};

#endif // IGLOBAL_H

