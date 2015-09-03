#ifndef IG_H
#define IG_H


#include "IGCommon.h"


class IG {    
    igraph_t graph;
    igVector weights;
    bool weighted;

    void empty() { igraph_empty(&graph, 0, false); }

    void igConstructorCheck(int err) {
        if (! err) return;
        empty(); // make sure 'graph' is not left uninitialized
        std::ostringstream msg;
        msg << "igraph returned with error: " << igraph_strerror(err);
        throw mma::LibraryError(msg.str());
    }

    void destroy() {
        igraph_destroy(&graph);
        clearWeights();
    }

    igraph_bliss_sh_t blissIntToSplitting(mint sh) const {
        switch (sh) {
        case 0: return IGRAPH_BLISS_F;
        case 1: return IGRAPH_BLISS_FL;
        case 2: return IGRAPH_BLISS_FLM;
        case 3: return IGRAPH_BLISS_FM;
        case 4: return IGRAPH_BLISS_FS;
        case 5: return IGRAPH_BLISS_FSM;
        default: throw mma::LibraryError("bliss: Unknown splitting heuristic.");
        }
    }

public:
    IG() : weighted{false} { empty(); }

    ~IG() {
        igraph_destroy(&graph);
    }

    // Create (basic)

    void fromEdgeList(mma::RealTensorRef v, mint n, bool directed) {
        destroy();
        igraph_vector_t edgelist = igVectorView(v);
        igConstructorCheck(igraph_create(&graph, &edgelist, n, directed));
    }

    // Weights

    void setWeights(mma::RealTensorRef w) {
        weighted = true;
        weights.copyFromMTensor(w);
    }

    void clearWeights() {
        weighted = false;
        weights.clear();
    }

    mma::RealTensorRef getWeights() const {
        if (! weighted)
            mma::message("Graph is not weighted. Returning empty weight vector.", mma::M_WARNING);
        return weights.makeMTensor();
    }

    bool weightedQ() const { return weighted; }

    // Create (games)

    void degreeSequenceGame(mma::RealTensorRef outdeg, mma::RealTensorRef indeg, mint method) {
        igraph_vector_t ig_indeg = igVectorView(indeg);
        igraph_vector_t ig_outdeg = igVectorView(outdeg);
        igraph_degseq_t ig_method;
        switch (method) {
        case 0: ig_method = IGRAPH_DEGSEQ_SIMPLE; break;
        case 1: ig_method = IGRAPH_DEGSEQ_SIMPLE_NO_MULTIPLE; break;
        case 2: ig_method = IGRAPH_DEGSEQ_VL; break;
        default: throw mma::LibraryError("degreeSequenceGame: unknown method option.");
        }

        destroy();
        int err;
        if (indeg.length() == 0)
            err = igraph_degree_sequence_game(&graph, &ig_outdeg, NULL, ig_method);
        else
            err = igraph_degree_sequence_game(&graph, &ig_outdeg, &ig_indeg, ig_method);
        igConstructorCheck(err);
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
        const igraph_vector_t *weightList = weighted ? &weights.vec : NULL;
        igCheck(igraph_betweenness(&graph, &vec.vec, igraph_vss_all(), true, weightList, true));
        return vec.makeMTensor();
    }

    mma::RealTensorRef edgeBetweenness() const {
        igVector vec;
        const igraph_vector_t *weightList = weighted ? &weights.vec : NULL;
        igCheck(igraph_edge_betweenness(&graph, &vec.vec, true, weightList));
        return vec.makeMTensor();
    }

    mma::RealTensorRef closeness(bool normalized) const {
        igVector vec;
        const igraph_vector_t *weightList = weighted ? &weights.vec : NULL;
        igCheck(igraph_closeness(&graph, &vec.vec, igraph_vss_all(), IGRAPH_OUT, weightList, normalized));
        return vec.makeMTensor();
    }

    // Centrality measures (estimates)

    mma::RealTensorRef betweennessEstimate(double cutoff) const {
        igVector vec;
        const igraph_vector_t *weightList = weighted ? &weights.vec : NULL;
        igCheck(igraph_betweenness_estimate(&graph, &vec.vec, igraph_vss_all(), true, cutoff, weightList, true));
        return vec.makeMTensor();
    }

    mma::RealTensorRef edgeBetweennessEstimate(double cutoff) const {
        igVector vec;
        const igraph_vector_t *weightList = weighted ? &weights.vec : NULL;
        igCheck(igraph_edge_betweenness_estimate(&graph, &vec.vec, true, cutoff, weightList));
        return vec.makeMTensor();
    }

    mma::RealTensorRef closenessEstimate(double cutoff, bool normalized) const {
        igVector vec;
        const igraph_vector_t *weightList = weighted ? &weights.vec : NULL;
        igCheck(igraph_closeness_estimate(&graph, &vec.vec, igraph_vss_all(), IGRAPH_OUT, cutoff, weightList, normalized));
        return vec.makeMTensor();
    }


    // Randomize

    void rewire(mint n, bool loops) {
        if (n > std::numeric_limits<igraph_integer_t>::max())
            throw mma::LibraryError("rewire: Requested number of rewiring trials too large.");
        igCheck(igraph_rewire(&graph, n, loops ? IGRAPH_REWIRING_SIMPLE_LOOPS : IGRAPH_REWIRING_SIMPLE));
    }

    void rewireEdges(double prob, bool loops, bool multiple) {
        igCheck(igraph_rewire_edges(&graph, prob, loops, multiple));
    }

    // Isomorphism (general)

    bool isomorphic(const IG &ig) const {
        igraph_bool_t res;
        igCheck(igraph_isomorphic(&graph, &ig.graph, &res));
        return res;
    }

    bool subisomorphic(const IG &ig) const {
        igraph_bool_t res;
        igCheck(igraph_subisomorphic(&graph, &ig.graph, &res));
        return res;
    }

    mint isoclass() const {
        igraph_integer_t res;
        igCheck(igraph_isoclass(&graph, &res));
        return res;
    }

    // Isomorphism (bliss)

    mma::RealTensorRef blissCanonicalPermutation(mint splitting) {
        igVector vec;
        igCheck(igraph_canonical_permutation(&graph, &vec.vec, blissIntToSplitting(splitting), NULL));
        return vec.makeMTensor();
    }

    bool blissIsomorphic(const IG &ig, mint splitting) {
        igraph_bool_t res;
        igCheck(igraph_isomorphic_bliss(&graph, &ig.graph, &res, NULL, NULL, blissIntToSplitting(splitting), blissIntToSplitting(splitting), NULL, NULL));
        return res;
    }

    mma::RealTensorRef blissFindIsomorphism(const IG &ig, mint splitting) {
        igraph_bool_t res;
        igVector map;
        igCheck(igraph_isomorphic_bliss(&graph, &ig.graph, &res, &map.vec, NULL, blissIntToSplitting(splitting), blissIntToSplitting(splitting), NULL, NULL));
        if (res)
            return map.makeMTensor();
        else
            return mma::makeVector<double>(0);
    }

    void blissAutomorphismsCount(MLINK link) {
        igraph_bliss_info_t info;
        mlStream ml{link, "blissAutomorphismsCount"};
        int splitting;

        ml >> mlCheckArgs(1) >> splitting;

        igCheck(igraph_automorphisms(&graph, blissIntToSplitting(splitting), &info));

        ml.newPacket();
        ml << info.group_size;
        std::free(info.group_size);
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
        igCheck(igraph_feedback_arc_set(&graph, &vec.vec, weighted ? &weights.vec : NULL, exact ? IGRAPH_FAS_EXACT_IP : IGRAPH_FAS_APPROX_EADES));
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

    mma::RealTensorRef motifs(mint size, mma::RealTensorRef cut_prob) const {
        igVector vec;
        igraph_vector_t ig_cut_prob = igVectorView(cut_prob);
        igCheck(igraph_motifs_randesu(&graph, &vec.vec, size, &ig_cut_prob));
        return vec.makeMTensor();
    }

    mint motifsNo(mint size, mma::RealTensorRef cut_prob) const {
        igraph_integer_t res;
        igraph_vector_t ig_cut_prob = igVectorView(cut_prob);
        igCheck(igraph_motifs_randesu_no(&graph, &res, size, &ig_cut_prob));
        return res;
    }

    mint motifsEstimate(mint size, mma::RealTensorRef cut_prob, mint sample_size) const {
        igraph_integer_t res;
        igraph_vector_t ig_cut_prob = igVectorView(cut_prob);
        igCheck(igraph_motifs_randesu_estimate(&graph, &res, size, &ig_cut_prob, sample_size, NULL));
        return res;
    }

    // Shortest paths

    mma::RealMatrixRef shortestPaths() const {
        igMatrix res;
        igCheck(igraph_shortest_paths(&graph, &res.mat, igraph_vss_all(), igraph_vss_all(), IGRAPH_OUT));
        return res.makeMTensor();
    }

    // Cliques

    void cliques(MLINK link) const {
        mlStream ml{link, "cliques"};
        int min, max;
        ml >> mlCheckArgs(2) >> min >> max;

        igList list;
        igCheck(igraph_cliques(&graph, &list.list, min, max));

        ml.newPacket();
        ml << list;
    }

    void maximalCliques(MLINK link) const {
        mlStream ml{link, "maximalCliques"};
        int min, max;
        ml >> mlCheckArgs(2) >> min >> max;

        igList list;
        igCheck(igraph_maximal_cliques(&graph, &list.list, min, max));

        ml.newPacket();
        ml << list;
    }

    mint maximalCliquesCount(int min, int max) const {
        igraph_integer_t res;
        igCheck(igraph_maximal_cliques_count(&graph, &res, min, max));
        return res;
    }

    void largestCliques(MLINK link) const {
        mlStream ml{link, "largestCliques"};
        ml >> mlCheckArgs(0);

        igList list;
        igCheck(igraph_largest_cliques(&graph, &list.list));

        ml.newPacket();
        ml << list;
    }

    mint cliqueNumber() const {
        igraph_integer_t res;
        igCheck(igraph_clique_number(&graph, &res));
        return res;
    }

    // Independent vertex sets

    void independentVertexSets(MLINK link) const {
        mlStream ml{link, "independentVertexSets"};
        int min, max;
        ml >> mlCheckArgs(2) >> min >> max;

        igList list;
        igCheck(igraph_independent_vertex_sets(&graph, &list.list, min, max));

        ml.newPacket();
        ml << list;
    }

    void largestIndependentVertexSets(MLINK link) const {
        mlStream ml{link, "largestIndependentVertexSets"};
        ml >> mlCheckArgs(0);

        igList list;
        igCheck(igraph_largest_independent_vertex_sets(&graph, &list.list));

        ml.newPacket();
        ml << list;
    }

    void maximalIndependentVertexSets(MLINK link) const {
        mlStream ml{link, "maximalIndependentVertexSets"};
        ml >> mlCheckArgs(0);

        igList list;
        igCheck(igraph_maximal_independent_vertex_sets(&graph, &list.list));

        ml.newPacket();
        ml << list;
    }

    mint independenceNumber() const {
        igraph_integer_t res;
        igCheck(igraph_independence_number(&graph, &res));
        return res;
    }
};

#endif // IG_H

