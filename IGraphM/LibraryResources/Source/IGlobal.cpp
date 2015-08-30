
#include "IGlobal.h"


int igInterruptionHandler(void *) {
    if (mma::libData->AbortQ()) {
        IGRAPH_FINALLY_FREE();
        return IGRAPH_INTERRUPTED;
    }
    return IGRAPH_SUCCESS;
}


void igWarningHandler(const char *reason, const char *file, int line, int igraph_errno) {
    mma::message(reason, mma::M_WARNING);
}


void igErrorHandler(const char *reason, const char *file, int line, int igraph_errno) {
    mma::message(reason, mma::M_ERROR);
    IGRAPH_FINALLY_FREE();
}
