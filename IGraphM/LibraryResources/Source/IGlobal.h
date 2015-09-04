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

    const char *compilationDate() { return __DATE__; }

    void seedRandom(mint s) {
        igCheck(igraph_rng_seed(igraph_rng_default(), s));
    }

    // Graph related functions that do not use the graph data structure

    bool graphicalQ(mma::RealTensorRef outdeg, mma::RealTensorRef indeg) {
        igraph_vector_t ig_outdeg = igVectorView(outdeg);
        igraph_vector_t ig_indeg  = igVectorView(indeg);
        igraph_bool_t res;
        if (indeg.length() == 0)
            igCheck(igraph_is_graphical_degree_sequence(&ig_outdeg, NULL, &res));
        else
            igCheck(igraph_is_graphical_degree_sequence(&ig_outdeg, &ig_indeg, &res));
        return res;
    }
};

#endif // IGLOBAL_H

