#ifndef IG_H
#define IG_H

#include "IGCommon.h"

#include <list>

class IG;

extern std::map<mint, IG *> IG_collection; // TODO this is a hack pending proper implementation in LTemplate

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

    // return the weights if weighted, return NULL otherwise
    // use this to pass weights to igraph functions
    const igraph_vector_t *passWeights() const { return weighted ? &weights.vec : NULL; }

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

    void fromLCF(mint n, mma::RealTensorRef v, mint repeats) {
        destroy();
        igraph_vector_t shifts = igVectorView(v);
        igConstructorCheck(igraph_lcf_vector(&graph, n, &shifts, repeats));
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

    void kRegularGame(mint n, mint k, bool directed, bool multiple) {
        destroy();
        igConstructorCheck(igraph_k_regular_game(&graph, n, k, directed, multiple));
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

    bool connectedQ(bool strong) const {
        igraph_bool_t res;
        igCheck(igraph_is_connected(&graph, &res, strong ? IGRAPH_STRONG : IGRAPH_WEAK));
        return res;
    }

    // Centrality measures

    mma::RealTensorRef betweenness() const {
        igVector vec;
        igCheck(igraph_betweenness(&graph, &vec.vec, igraph_vss_all(), true, passWeights(), true));
        return vec.makeMTensor();
    }

    mma::RealTensorRef edgeBetweenness() const {
        igVector vec;
        igCheck(igraph_edge_betweenness(&graph, &vec.vec, true, passWeights()));
        return vec.makeMTensor();
    }

    mma::RealTensorRef closeness(bool normalized) const {
        igVector vec;
        igCheck(igraph_closeness(&graph, &vec.vec, igraph_vss_all(), IGRAPH_OUT, passWeights(), normalized));
        return vec.makeMTensor();
    }

    // Centrality measures (estimates)

    mma::RealTensorRef betweennessEstimate(double cutoff) const {
        igVector vec;
        igCheck(igraph_betweenness_estimate(&graph, &vec.vec, igraph_vss_all(), true, cutoff, passWeights(), true));
        return vec.makeMTensor();
    }

    mma::RealTensorRef edgeBetweennessEstimate(double cutoff) const {
        igVector vec;
        igCheck(igraph_edge_betweenness_estimate(&graph, &vec.vec, true, cutoff, passWeights()));
        return vec.makeMTensor();
    }

    mma::RealTensorRef closenessEstimate(double cutoff, bool normalized) const {
        igVector vec;
        igCheck(igraph_closeness_estimate(&graph, &vec.vec, igraph_vss_all(), IGRAPH_OUT, cutoff, passWeights(), normalized));
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

    void blissAutomorphismCount(MLINK link) {
        igraph_bliss_info_t info;
        mlStream ml{link, "blissAutomorphismsCount"};
        int splitting;

        ml >> mlCheckArgs(1) >> splitting;

        igCheck(igraph_automorphisms(&graph, blissIntToSplitting(splitting), &info));

        ml.newPacket();
        ml << info.group_size;
        std::free(info.group_size);
    }

    // Isomorphism (VF2)

    bool vf2Isomorphic(
            const IG &ig, mma::IntTensorRef vcol1, mma::IntTensorRef vcol2,
                          mma::IntTensorRef ecol1, mma::IntTensorRef ecol2) const
    {
        igIntVector vc1; vc1.copyFromMTensor(vcol1);
        igIntVector vc2; vc2.copyFromMTensor(vcol2);
        igIntVector ec1; ec1.copyFromMTensor(ecol1);
        igIntVector ec2; ec2.copyFromMTensor(ecol2);

        igraph_bool_t res;
        igCheck(igraph_isomorphic_vf2(
                    &graph, &ig.graph,
                    vcol1.length() == 0 ? NULL : &vc1.vec, vcol2.length() == 0 ? NULL : &vc2.vec,
                    ecol1.length() == 0 ? NULL : &ec1.vec, ecol2.length() == 0 ? NULL : &ec2.vec,
                    &res, NULL, NULL, NULL, NULL, NULL));
        return res;
    }

    void vf2FindIsomorphisms(MLINK link) const {
        mlStream ml(link, "vf2Isomorphism");

        mint id; // expression ID
        igIntVector vc1, vc2, ec1, ec2;

        struct VF2data {
            std::list<igVector> list;
            long remaining; // remaining number of isomorphisms to find, negative value will run until all are found
        } vf2data;

        ml >> mlCheckArgs(6) >> id >> vf2data.remaining >> vc1 >> vc2 >> ec1 >> ec2;

        struct {
            static igraph_bool_t handle(const igraph_vector_t *map12,  const igraph_vector_t *map21, void *arg) {
                VF2data &data = *static_cast<VF2data *>(arg);
                data.list.push_back(map12);
                data.remaining--;
                return data.remaining != 0; // negative will run until all are found
            }
        } isohandler;

        igCheck(igraph_isomorphic_function_vf2(
                    &graph, &IG_collection[id]->graph,
                    vc1.length() == 0 ? NULL : &vc1.vec, vc2.length() == 0 ? NULL : &vc2.vec,
                    ec1.length() == 0 ? NULL : &ec1.vec, ec2.length() == 0 ? NULL : &ec2.vec,
                    NULL, NULL, &isohandler.handle, NULL, NULL, &vf2data));

        ml.newPacket();
        ml << vf2data.list;
    }

    bool vf2Subisomorphic(
            const IG &ig, mma::IntTensorRef vcol1, mma::IntTensorRef vcol2,
                          mma::IntTensorRef ecol1, mma::IntTensorRef ecol2) const
    {
        igIntVector vc1; vc1.copyFromMTensor(vcol1);
        igIntVector vc2; vc2.copyFromMTensor(vcol2);
        igIntVector ec1; ec1.copyFromMTensor(ecol1);
        igIntVector ec2; ec2.copyFromMTensor(ecol2);

        igraph_bool_t res;
        igCheck(igraph_subisomorphic_vf2(
                    &graph, &ig.graph,
                    vcol1.length() == 0 ? NULL : &vc1.vec, vcol2.length() == 0 ? NULL : &vc2.vec,
                    ecol1.length() == 0 ? NULL : &ec1.vec, ecol2.length() == 0 ? NULL : &ec2.vec,
                    &res, NULL, NULL, NULL, NULL, NULL));
        return res;
    }

    void vf2FindSubisomorphisms(MLINK link) const {
        mlStream ml(link, "vf2Isomorphism");

        mint id; // expression ID
        igIntVector vc1, vc2, ec1, ec2;

        struct VF2data {
            std::list<igVector> list;
            long remaining; // remaining number of isomorphisms to find, negative value will run until all are found
        } vf2data;

        ml >> mlCheckArgs(6) >> id >> vf2data.remaining >> vc1 >> vc2 >> ec1 >> ec2;

        struct {
            static igraph_bool_t handle(const igraph_vector_t *map12,  const igraph_vector_t *map21, void *arg) {
                VF2data &data = *static_cast<VF2data *>(arg);
                data.list.push_back(map21);
                data.remaining--;
                return data.remaining != 0; // negative will run until all are found
            }
        } isohandler;

        igCheck(igraph_subisomorphic_function_vf2(
                    &graph, &IG_collection[id]->graph,
                    vc1.length() == 0 ? NULL : &vc1.vec, vc2.length() == 0 ? NULL : &vc2.vec,
                    ec1.length() == 0 ? NULL : &ec1.vec, ec2.length() == 0 ? NULL : &ec2.vec,
                    NULL, NULL, &isohandler.handle, NULL, NULL, &vf2data));

        ml.newPacket();
        ml << vf2data.list;
    }

    mint vf2AutomorphismCount() const {
        igraph_integer_t res;
        igCheck(igraph_count_isomorphisms_vf2(&graph, &graph, NULL, NULL, NULL, NULL, &res, NULL, NULL, NULL));
        return res;
    }

    mint vf2IsomorphismCount(
            const IG &ig, mma::IntTensorRef vcol1, mma::IntTensorRef vcol2,
                          mma::IntTensorRef ecol1, mma::IntTensorRef ecol2) const
    {
        igIntVector vc1; vc1.copyFromMTensor(vcol1);
        igIntVector vc2; vc2.copyFromMTensor(vcol2);
        igIntVector ec1; ec1.copyFromMTensor(ecol1);
        igIntVector ec2; ec2.copyFromMTensor(ecol2);

        igraph_integer_t res;

        igCheck(igraph_count_isomorphisms_vf2(
                    &graph, &ig.graph,
                    vcol1.length() == 0 ? NULL : &vc1.vec, vcol2.length() == 0 ? NULL : &vc2.vec,
                    ecol1.length() == 0 ? NULL : &ec1.vec, ecol2.length() == 0 ? NULL : &ec2.vec,
                    &res, NULL, NULL, NULL));
        return res;
    }

    mint vf2SubisomorphismCount(
            const IG &ig, mma::IntTensorRef vcol1, mma::IntTensorRef vcol2,
                          mma::IntTensorRef ecol1, mma::IntTensorRef ecol2) const
    {
        igIntVector vc1; vc1.copyFromMTensor(vcol1);
        igIntVector vc2; vc2.copyFromMTensor(vcol2);
        igIntVector ec1; ec1.copyFromMTensor(ecol1);
        igIntVector ec2; ec2.copyFromMTensor(ecol2);

        igraph_integer_t res;
        igCheck(igraph_count_subisomorphisms_vf2(
                    &graph, &ig.graph,
                    vcol1.length() == 0 ? NULL : &vc1.vec, vcol2.length() == 0 ? NULL : &vc2.vec,
                    ecol1.length() == 0 ? NULL : &ec1.vec, ecol2.length() == 0 ? NULL : &ec2.vec,
                    &res, NULL, NULL, NULL));
        return res;
    }

    // Isomorphism (LAD)

    bool ladSubisomorphic(const IG &ig, bool induced) const {
        igraph_bool_t res;
        igCheck(igraph_subisomorphic_lad(&ig.graph, &graph, NULL, &res, NULL, NULL, induced, 0));
        return res;
    }

    mma::RealTensorRef ladGetSubisomorphism(const IG &ig, bool induced) const {
        igraph_bool_t iso;
        igVector map;
        igCheck(igraph_subisomorphic_lad(&ig.graph, &graph, NULL, &iso, &map.vec, NULL, induced, 0));
        return map.makeMTensor();
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
        igCheck(igraph_feedback_arc_set(&graph, &vec.vec, passWeights(), exact ? IGRAPH_FAS_EXACT_IP : IGRAPH_FAS_APPROX_EADES));
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

    // Graph drawing (layouts)

    mma::RealTensorRef layoutRandom() const {
        igMatrix mat;
        igCheck(igraph_layout_random(&graph, &mat.mat));
        return mat.makeMTensor();
    }

    mma::RealTensorRef layoutCircle() const {
        igMatrix mat;
        igCheck(igraph_layout_circle(&graph, &mat.mat, igraph_vss_all()));
        return mat.makeMTensor();
    }

    mma::RealTensorRef layoutSphere() const  {
        igMatrix mat;
        igCheck(igraph_layout_sphere(&graph, &mat.mat));
        return mat.makeMTensor();
    }

    mma::RealTensorRef layoutGraphOpt(
            mma::RealMatrixRef initial, bool use_seed,
            mint niter,
            double node_charge, double node_mass, double spring_length,
            double spring_constant, double max_sa_movement) const
    {
        igMatrix mat;
        mat.copyFromMTensor(initial);
        igCheck(igraph_layout_graphopt(&graph, &mat.mat, niter, node_charge, node_mass, spring_length, spring_constant, max_sa_movement, use_seed));
        return mat.makeMTensor();
    }

    mma::RealTensorRef layoutKamadaKawai(
            mma::RealMatrixRef initial, bool use_seed,
            mint maxiter, double epsilon, double kkconst) const
    {
        igMatrix mat;
        mat.copyFromMTensor(initial);
        igCheck(igraph_layout_kamada_kawai(
                    &graph, &mat.mat, use_seed, maxiter, epsilon, kkconst, passWeights(),
                    NULL, NULL, NULL, NULL));
        return mat.makeMTensor();
    }

    mma::RealTensorRef layoutKamadaKawai3D(
            mma::RealMatrixRef initial, bool use_seed,
            mint maxiter, double epsilon, double kkconst) const
    {
        igMatrix mat;
        mat.copyFromMTensor(initial);
        igCheck(igraph_layout_kamada_kawai_3d(
                    &graph, &mat.mat, use_seed, maxiter, epsilon, kkconst, passWeights(),
                    NULL, NULL, NULL, NULL, NULL, NULL));
        return mat.makeMTensor();
    }

    mma::RealTensorRef layoutFruchtermanReingold(
            mma::RealMatrixRef initial, bool use_seed,
            mint niter, double start_temp, mint grid_method) const
    {
        igMatrix mat;
        mat.copyFromMTensor(initial);
        igraph_layout_grid_t grid;
        switch (grid_method) {
        case 0: grid = IGRAPH_LAYOUT_GRID; break;
        case 1: grid = IGRAPH_LAYOUT_NOGRID; break;
        case 2: grid = IGRAPH_LAYOUT_AUTOGRID; break;
        default: throw mma::LibraryError("layoutFruchtermanReingold: unknown method.");
        }

        igCheck(igraph_layout_fruchterman_reingold(
                    &graph, &mat.mat, use_seed, niter, start_temp, grid, passWeights(),
                    NULL, NULL, NULL, NULL));
        return mat.makeMTensor();
    }

    mma::RealTensorRef layoutFruchtermanReingold3D(
            mma::RealMatrixRef initial, bool use_seed,
            mint niter, double start_temp) const
    {
        igMatrix mat;
        mat.copyFromMTensor(initial);
        igCheck(igraph_layout_fruchterman_reingold_3d(
                    &graph, &mat.mat, use_seed, niter, start_temp, passWeights(),
                    NULL, NULL, NULL, NULL, NULL, NULL));
        return mat.makeMTensor();
    }

    // Clusterig coefficient

    double transitivityUndirected() const {
        double res;
        igCheck(igraph_transitivity_undirected(&graph, &res, IGRAPH_TRANSITIVITY_ZERO));
        return res;
    }

    mma::RealTensorRef transitivityLocalUndirected() const {
        igVector vec;
        igCheck(igraph_transitivity_local_undirected(&graph, &vec.vec, igraph_vss_all(), IGRAPH_TRANSITIVITY_ZERO));
        return vec.makeMTensor();
    }

    double transitivityAverageLocalUndirected() const {
        double res;
        igCheck(igraph_transitivity_avglocal_undirected(&graph, &res, IGRAPH_TRANSITIVITY_ZERO));
        return res;
    }

    mma::RealTensorRef transitivityBarrat() const {
        igVector vec;
        igCheck(igraph_transitivity_barrat(&graph, &vec.vec, igraph_vss_all(), passWeights(), IGRAPH_TRANSITIVITY_ZERO));
        return vec.makeMTensor();
    }

    // Similarity measures

    mma::RealTensorRef similarityCocitation(mma::RealTensorRef vs) const {
        igMatrix mat;
        igraph_vector_t vsvec = igVectorView(vs);
        igCheck(igraph_cocitation(&graph, &mat.mat, vs.length() == 0 ? igraph_vss_all() : igraph_vss_vector(&vsvec)));
        return mat.makeMTensor();
    }

    mma::RealTensorRef similarityBibcoupling(mma::RealTensorRef vs) const {
        igMatrix mat;
        igraph_vector_t vsvec = igVectorView(vs);
        igCheck(igraph_bibcoupling(&graph, &mat.mat, vs.length() == 0 ? igraph_vss_all() : igraph_vss_vector(&vsvec)));
        return mat.makeMTensor();
    }

    mma::RealTensorRef similarityJaccard(mma::RealTensorRef vs, bool loops) const {
        igMatrix mat;
        igraph_vector_t vsvec = igVectorView(vs);
        igCheck(igraph_similarity_jaccard(&graph, &mat.mat, vs.length() == 0 ? igraph_vss_all() : igraph_vss_vector(&vsvec), IGRAPH_OUT, loops));
        return mat.makeMTensor();
    }

    mma::RealTensorRef similarityDice(mma::RealTensorRef vs, bool loops) const {
        igMatrix mat;
        igraph_vector_t vsvec = igVectorView(vs);
        igCheck(igraph_similarity_dice(&graph, &mat.mat, vs.length() == 0 ? igraph_vss_all() : igraph_vss_vector(&vsvec), IGRAPH_OUT, loops));
        return mat.makeMTensor();
    }

    mma::RealTensorRef similarityInverseLogWeighted(mma::RealTensorRef vs) const {
        igMatrix mat;
        igraph_vector_t vsvec = igVectorView(vs);
        igCheck(igraph_similarity_inverse_log_weighted(&graph, &mat.mat, vs.length() == 0 ? igraph_vss_all() : igraph_vss_vector(&vsvec), IGRAPH_OUT));
        return mat.makeMTensor();
    }

    // Chordal graphs

    mma::RealTensorRef maximumCardinalitySearch() const {
        igVector vec;
        igCheck(igraph_maximum_cardinality_search(&graph, &vec.vec, NULL));
        return vec.makeMTensor();
    }

    bool chordalQ() const {
        igraph_bool_t res;
        igCheck(igraph_is_chordal(&graph, NULL, NULL, &res, NULL, NULL));
        return res;
    }

    // Note that this is a constructor!
    mma::RealTensorRef chordalCompletion(const IG &source) {
        igraph_bool_t chordal;
        igVector fillin;
        destroy();
        igConstructorCheck(igraph_is_chordal(&source.graph, NULL, NULL, &chordal, &fillin.vec, &graph));
        return fillin.makeMTensor();
    }

    // Vertex separators

    void minimumSizeSeparators(MLINK link) const {
        igList list;
        mlStream ml(link, "minimumSizeSeparators");
        ml >> mlCheckArgs(0);

        igCheck(igraph_minimum_size_separators(&graph, &list.list));

        ml.newPacket();
        ml << list;
    }

    // Connected components

    mma::RealTensorRef articulationPoints() const {
        igVector vec;
        igCheck(igraph_articulation_points(&graph, &vec.vec));
        return vec.makeMTensor();
    }

    void biconnectedComponents(MLINK link) const {
        mlStream ml(link, "biconnectedComponents");
        ml >> mlCheckArgs(0);

        igList list;
        igraph_integer_t count;
        igCheck(igraph_biconnected_components(&graph, &count, NULL, NULL, &list.list, NULL));

        ml.newPacket();
        ml << list;
    }
};

#endif // IG_H

