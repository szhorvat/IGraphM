
#include "IGlobal.h"
#include <cstring>


int igInterruptionHandler(void *) {
    if (mma::libData->AbortQ()) {
        IGRAPH_FINALLY_FREE();
        return IGRAPH_INTERRUPTED;
    }
    return IGRAPH_SUCCESS;
}


void igWarningHandler(const char *reason, const char *file, int line, int igraph_errno) {
    std::ostringstream msg;
    msg << file << ":" << line << " - " << reason;
    mma::message(msg.str(), mma::M_WARNING);
}


void igErrorHandler(const char *reason, const char *file, int line, int igraph_errno) {
    // avoid printing empty messages
    if (strlen(reason) != 0) {
        std::ostringstream msg;
        msg << file << ":" << line << " - " << reason;
        mma::message(msg.str(), mma::M_ERROR);
    }
    IGRAPH_FINALLY_FREE();
}
