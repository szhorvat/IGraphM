#ifndef IG_H
#define IG_H


extern "C" { // workaround for igraph_version() C++ compatibility bug
#include <igraph/igraph.h>
}

#include "LTemplate.h"
#include <algorithm>


struct igVector {
    igraph_vector_t vec;

    igVector() { igraph_vector_init(&vec, 0); }
    ~igVector() { igraph_vector_destroy(&vec); }

    long length() const { return igraph_vector_size(&vec); }

    igraph_real_t *begin() { return &(VECTOR(vec)[0]); }
    igraph_real_t *end() { return begin() + length(); }
};


inline igraph_vector_t ig_view(mma::RealTensorRef &t) {
    igraph_vector_t vec;
    igraph_vector_view(&vec, t.data(), t.length());
    return vec;
}


inline void igCheck(int err) {
    if (err)
        throw mma::LibraryError(igraph_strerror(err));
}


class IG {
    igraph_t graph;

public:
    IG() {
        igraph_empty(&graph, 0, false);
    }

    ~IG() {
        igraph_destroy(&graph);
    }

    // Create

    void fromEdgeList(mma::RealTensorRef v, mint n, bool directed) {
        igraph_destroy(&graph);
        igraph_vector_t edgelist = ig_view(v);
        igraph_create(&graph, &edgelist, n, directed); // check for error manually
    }

    // Structure

    mma::RealTensorRef edgeList() const {
        igVector vec;
        igCheck(igraph_get_edgelist(&graph, &vec.vec, false));
        mma::RealTensorRef res = mma::makeMatrix<double>(vec.length() / 2, 2);
        std::copy(vec.begin(), vec.end(), res.begin());
        return res;
    }

    mint edgeCount() const { return igraph_ecount(&graph); }

    mint vertexCount() const { return igraph_vcount(&graph); }

    // Testing

    bool directedQ() const { return igraph_is_directed(&graph); }

    bool dagQ() const {
        igraph_bool_t res;
        igCheck(igraph_is_dag(&graph, &res));
        return res;
    }

    bool simpleQ() const {
        igraph_bool_t res;
        igCheck(igraph_is_simple(&graph, &res));
        return res;
    }

    // TODO handle strong/weak for directed
    bool connectedQ() const {
        igraph_bool_t res;
        igCheck(igraph_is_connected(&graph, &res, IGRAPH_STRONG));
        return res;
    }

    // Centrality measures

    mma::RealTensorRef betweenness() const {
        igVector vec;
        igCheck(igraph_betweenness(&graph, &vec.vec, igraph_vss_all(), true, NULL, true));
        mma::RealTensorRef res = mma::makeVector<double>(vec.length());
        std::copy(vec.begin(), vec.end(), res.begin());
        return res;
    }

    mma::RealTensorRef closeness(bool normalized) const {
        igVector vec;
        igCheck(igraph_closeness(&graph, &vec.vec, igraph_vss_all(), IGRAPH_OUT, NULL, normalized));
        mma::RealTensorRef res = mma::makeVector<double>(vec.length());
        std::copy(vec.begin(), vec.end(), res.begin());
        return res;
    }

    mma::RealTensorRef edgeBetweenness() const {
        igVector vec;
        igCheck(igraph_edge_betweenness(&graph, &vec.vec, true, NULL));
        mma::RealTensorRef res = mma::makeVector<double>(vec.length());
        std::copy(vec.begin(), vec.end(), res.begin());
        return res;
    }

    // Randomize

    void rewire(mint n, bool loops) {
        if (n > std::numeric_limits<igraph_integer_t>::max())
            throw mma::LibraryError("igraph rewire: Requested number rewiring trials too large.");
        igCheck(igraph_rewire(&graph, n, loops ? IGRAPH_REWIRING_SIMPLE_LOOPS : IGRAPH_REWIRING_SIMPLE));
    }

    // Isomorphism

    bool isomorphic(IG &ig) {
        igraph_bool_t res;
        igraph_isomorphic(&graph, &ig.graph, &res);
        return res;
    }

    bool subisomorphic(IG &ig) {
        igraph_bool_t res;
        igraph_subisomorphic(&graph, &ig.graph, &res);
        return res;
    }
};

#endif // IG_H

