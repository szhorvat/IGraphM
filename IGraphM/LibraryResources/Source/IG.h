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

    const igraph_real_t *begin() const { return &(VECTOR(vec)[0]); }
    const igraph_real_t *end() const { return begin() + length(); }

    mma::RealTensorRef makeMTensor() const {
        mma::RealTensorRef res = mma::makeVector<double>(length());
        std::copy(begin(), end(), res.begin());
        return res;
    }
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
        igraph_create(&graph, &edgelist, n, directed); // TODO check for error manually
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
        return vec.makeMTensor();
    }

    mma::RealTensorRef closeness(bool normalized) const {
        igVector vec;
        igCheck(igraph_closeness(&graph, &vec.vec, igraph_vss_all(), IGRAPH_OUT, NULL, normalized));
        return vec.makeMTensor();
    }

    mma::RealTensorRef edgeBetweenness() const {
        igVector vec;
        igCheck(igraph_edge_betweenness(&graph, &vec.vec, true, NULL));
        return vec.makeMTensor();
    }

    // Randomize

    void rewire(mint n, bool loops) {
        if (n > std::numeric_limits<igraph_integer_t>::max())
            throw mma::LibraryError("igraph rewire: Requested number of rewiring trials too large.");
        igCheck(igraph_rewire(&graph, n, loops ? IGRAPH_REWIRING_SIMPLE_LOOPS : IGRAPH_REWIRING_SIMPLE));
    }

    void rewireEdges(double prob, bool loops, bool multiple) {
        igCheck(igraph_rewire_edges(&graph, prob, loops, multiple));
    }

    // Isomorphism

    bool isomorphic(IG &ig) {
        igraph_bool_t res;
        igCheck(igraph_isomorphic(&graph, &ig.graph, &res));
        return res;
    }

    bool subisomorphic(IG &ig) {
        igraph_bool_t res;
        igCheck(igraph_subisomorphic(&graph, &ig.graph, &res));
        return res;
    }

    // Topological sorting, directed acylic graphs

    // see also dagQ() under Testing

    mma::RealTensorRef topologicalSorting() const {
        igVector vec;
        igCheck(igraph_topological_sorting(&graph, &vec.vec, IGRAPH_OUT));
        return vec.makeMTensor();
    }

    mma::RealTensorRef feedbackArcSet(bool exact) const {
        igVector vec;
        igCheck(igraph_feedback_arc_set(&graph, &vec.vec, NULL, exact ? IGRAPH_FAS_EXACT_IP : IGRAPH_FAS_APPROX_EADES));
        return vec.makeMTensor();
    }

    // Motifs and subgraph counts

    mma::IntTensorRef dyadCensus() const {
        igraph_integer_t mut, asym, none;
        igCheck(igraph_dyad_census(&graph, &mut, &asym, &none));
        mma::IntTensorRef res = mma::makeVector<mint>(3);
        res[0] = mut;
        res[1] = asym;
        res[2] = none;
        return res;
    }

    mma::RealTensorRef triadCensus() const {
        igVector vec;
        igCheck(igraph_triad_census(&graph, &vec.vec));
        return vec.makeMTensor();
    }
};

#endif // IG_H

