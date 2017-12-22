/*
 * Copyright (c) 2017 Szabolcs Horvát.
 *
 * See the file LICENSE.txt for copying permission.
 */

#ifndef IGLOBAL_H
#define IGLOBAL_H

#include "IGCommon.h"

#include <random>
#include <cmath>


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

    bool infOrNanQ(mma::RealTensorRef t) {
        for (double *x = t.begin(); x != t.end(); ++x)
            if (std::isnan(*x) || std::isinf(*x))
                return true;
        return false;
    }

    mma::IntTensorRef incidenceToEdgeList(mma::SparseMatrixRef<mint> im, bool directed) {
        auto edgeList = mma::makeMatrix<mint>(im.cols(), 2);
        if (directed) {
            for (auto it = im.begin(); it != im.end(); ++it) {
                switch (*it) {
                case -1:
                    edgeList[2*it.col()] = it.row();
                    break;
                case  1:
                    edgeList[2*it.col() + 1] = it.row();
                    break;
                case  2:
                case -2:
                    edgeList[2*it.col()] = it.row();
                    edgeList[2*it.col() + 1] = it.row();
                    break;
                default:
                    throw mma::LibraryError("Invalid incidence matrix.");
                }
            }
        } else {
            for (auto &el : edgeList)
                el = -1;
            for (auto it = im.begin(); it != im.end(); ++it) {
                switch (*it) {
                case  1:
                    if (edgeList[2*it.col()] == -1)
                        edgeList[2*it.col()] = it.row();
                    else
                        edgeList[2*it.col() + 1] = it.row();
                    break;
                case  2:
                    edgeList[2*it.col()] = it.row();
                    edgeList[2*it.col() + 1] = it.row();
                    break;
                default:
                    throw mma::LibraryError("Invalid incidence matrix.");
                }
            }
        }
        return edgeList;
    }
};

#endif // IGLOBAL_H

