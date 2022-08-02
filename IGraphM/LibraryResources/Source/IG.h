/*
 * Copyright (c) 2016-2022 Szabolcs Horvát.
 *
 * See the file LICENSE.txt for copying permission.
 */

#ifndef IG_H
#define IG_H

#include "IGCommon.h"

#include <list>
#include <map>
#include <set>
#include <tuple>
#include <iomanip>
#include <cmath>

class IG;

class IG {
    igraph_t graph;
    igVector weights;
    bool weighted;

    void empty() { igraph_empty(&graph, 0, false); }   

    void igConstructorCheck(int err) {
        if (! err) return;
        empty(); // make sure 'graph' is not left uninitialized
        std::ostringstream msg;
        msg << "igraph returned with error: " << igraph_strerror(err) << ".";
        throw mma::LibraryError(msg.str());
    }

    void destroy() {
        igraph_destroy(&graph);
        clearWeights();
    }

    // return the weights if weighted, return NULL otherwise
    // use this to pass weights to igraph functions
    const igraph_vector_t *passWeights() const { return weighted ? &weights.vec : nullptr; }

    /* Community detection helpers */

    void computeMembership(const igraph_integer_t n_communities, const igMatrix &merges, igVector &membership) const {
        igraph_integer_t vc = igraph_vcount(&graph);
        igCheck(igraph_community_to_membership(
                    &merges.mat, vc, vc - n_communities /* steps */,
                    &membership.vec, nullptr /* csize */
                    ));
    }

    /* Isomorphism helpers */

    // translate integer to Bliss splitting heuristic enum value
    static igraph_bliss_sh_t blissIntToSplitting(mint sh) {
        switch (sh) {
        case 0: return IGRAPH_BLISS_F;
        case 1: return IGRAPH_BLISS_FS;
        case 2: return IGRAPH_BLISS_FL;
        case 3: return IGRAPH_BLISS_FM;
        case 4: return IGRAPH_BLISS_FSM;
        case 5: return IGRAPH_BLISS_FLM;
        default: throw mma::LibraryError("bliss: Unknown splitting heuristic.");
        }
    }

    // changes ig1 to have the directedness of ig2
    static void matchDirectedness(IG &ig1, const IG &ig2) {
        if (ig2.directedQ())
            igCheck(igraph_to_directed(&ig1.graph, IGRAPH_TO_DIRECTED_ARBITRARY));
        else
            igCheck(igraph_to_undirected(&ig1.graph, IGRAPH_TO_UNDIRECTED_EACH, nullptr));
    }

    // this function is called in isomorphism functions, to ensure that if one
    // graph has no edges, it will be converted to the directedness of the other graph.
    void emptyMatchDirectedness(IG &ig) {
        if (igraph_ecount(&graph) == 0) {
            matchDirectedness(*this, ig);
        }
        else if (igraph_ecount(&ig.graph) == 0) {
            matchDirectedness(ig, *this);
        }
    }

public:
    IG() : weighted{false} { empty(); }

    ~IG() {
        igraph_destroy(&graph);
    }

    IG(const IG &) = delete;
    IG &operator = (const IG &) = delete;

    // Create (basic)

    void fromEdgeList(mma::RealTensorRef v, mint n, bool directed) {
        destroy();
        igraph_vector_t edgelist = igVectorView(v);
        igConstructorCheck(igraph_create(&graph, &edgelist, n, directed));
    }

    // this overload is used in the implementation of other functions
    void fromEdgeList(const igVector &edges, mint n, bool directed) {
        destroy();
        igConstructorCheck(igraph_create(&graph, &edges.vec, n, directed));
    }

    void fromIncidenceMatrix(mma::SparseMatrixRef<mint> im, bool directed);

    /* Creates an undirected graph with n vertices and no edges. */
    void makeEdgeless(mint n) {
        destroy();
        igConstructorCheck(igraph_empty(&graph, n /* vertex count */, false /* undirected */));
    }

    void fromEdgeListML(MLINK link);

    void realizeDegreeSequence(mma::RealTensorRef outdegs, mma::RealTensorRef indegs, bool loops, bool multi, mint method) {
        destroy();
        igraph_vector_t out = igVectorView(outdegs);
        igraph_vector_t in = igVectorView(indegs);        

        igraph_edge_type_sw_t et = (loops ? IGRAPH_LOOPS_SW : IGRAPH_SIMPLE_SW) | (multi ? IGRAPH_MULTI_SW : IGRAPH_SIMPLE_SW);

        igraph_realize_degseq_t ig_method;
        switch (method) {
        case 0: ig_method = IGRAPH_REALIZE_DEGSEQ_SMALLEST; break;
        case 1: ig_method = IGRAPH_REALIZE_DEGSEQ_LARGEST; break;
        case 2: ig_method = IGRAPH_REALIZE_DEGSEQ_INDEX; break;
        default: throw mma::LibraryError("Unknown degree sequence realization method.");
        }

        igConstructorCheck(igraph_realize_degree_sequence(&graph, &out, indegs.size() == 0 ? nullptr : &in, et, ig_method));
    }

    void fromLCF(mint n, mma::RealTensorRef v, mint repeats) {
        destroy();
        igraph_vector_t shifts = igVectorView(v);
        igConstructorCheck(igraph_lcf_vector(&graph, n, &shifts, repeats));
    }

    void makeLattice(mma::RealTensorRef dims, mint nei, bool directed, bool mutual, bool periodic) {
        destroy();
        igraph_vector_t igdims = igVectorView(dims);
        igConstructorCheck(igraph_lattice(&graph, &igdims, nei, directed, mutual, periodic));
    }

    void graphAtlas(mint n) {
        destroy();
        igConstructorCheck(igraph_atlas(&graph, n));
    }

    void famous(const char *name) {
        destroy();
        igConstructorCheck(igraph_famous(&graph, name));
        mma::disownString(name);
    }

    void kautz(mint m, mint n) {
        destroy();
        igConstructorCheck(igraph_kautz(&graph, m, n));
    }

    void tree(mint n, mint k, bool directed) {
        destroy();
        igConstructorCheck(igraph_tree(&graph, n, k, directed ? IGRAPH_TREE_OUT : IGRAPH_TREE_UNDIRECTED));    
    }   

    void fromPrufer(mma::IntTensorRef prufer) {
        destroy();
        igIntVector vec;
        vec.copyFromMTensor(prufer);
        igConstructorCheck(igraph_from_prufer(&graph, &vec.vec));
    }

    mma::IntTensorRef toPrufer() const {
        igIntVector vec;
        igCheck(igraph_to_prufer(&graph, &vec.vec));
        return vec.makeMTensor();
    }

    void completeGraph(mint n, bool directed, bool loops) {
        destroy();
        igConstructorCheck(igraph_full(&graph, n, directed, loops));
    }

    void completeCitationGraph(mint n, bool directed) {
        destroy();
        igConstructorCheck(igraph_full_citation(&graph, n, directed));
    }

    void deBruijn(mint m, mint n) {
        destroy();
        igConstructorCheck(igraph_de_bruijn(&graph, m, n));
    }

    void extendedChordalRing(mint n, mma::RealMatrixRef mat, bool directed, bool selfLoops, bool multiEdges) {
        igMatrix w;

        destroy();
        w.copyFromMTensor(mat);
        igConstructorCheck(igraph_extended_chordal_ring(&graph, n, &w.mat, directed));
        if (! (multiEdges && selfLoops))
            igCheck(igraph_simplify(&graph, !multiEdges, !selfLoops, nullptr));
    }

    void coloredSimpleGraph(MLINK link) {
        mlStream ml{link, "coloredSimpleGraph"};

        mint id;
        ml >> mlCheckArgs(1) >> id;

        IG &ig = mma::getInstance<IG>(id);

        igIntVector vc, ec;

        destroy();
        igCheck(igraph_simplify_and_colorize(&ig.graph, &graph, &vc.vec, &ec.vec));

        ml.newPacket();
        ml << mlHead("List", 2)
           << vc << ec;
    }

    // Weights

    void setWeights(mma::RealTensorRef w) {
        weighted = true;
        weights.copyFromMTensor(w);
    }

    void clearWeights() {
        if (weighted) {
            weighted = false;
            weights.clear();
        }
    }

    mma::RealTensorRef getWeights() const {
        if (! weighted)
            mma::message("Graph is not weighted. Returning empty weight vector.", mma::M_WARNING);
        return weights.makeMTensor();
    }

    bool weightedQ() const { return weighted; }

    // Change directedness

    void makeDirected() {
        igraph_to_directed(&graph, IGRAPH_TO_DIRECTED_MUTUAL);
    }

    void makeUndirected() {
        igraph_to_undirected(&graph, IGRAPH_TO_UNDIRECTED_COLLAPSE, nullptr);
    }

    // Create (games)

    void treeGame(mint n, bool directed, mint method) {
        destroy();
        igraph_random_tree_t ig_method;
        switch (method) {
        case 0:
            ig_method = IGRAPH_RANDOM_TREE_PRUFER;
            break;
        case 1:
            ig_method = IGRAPH_RANDOM_TREE_LERW;
            break;
        default:
            throw mma::LibraryError("Unknown random tree generation method.");
        }
        igConstructorCheck(igraph_tree_game(&graph, n, directed, ig_method));
    }

    void degreeSequenceGame(mma::RealTensorRef outdeg, mma::RealTensorRef indeg, mint method) {
        igraph_vector_t ig_indeg = igVectorView(indeg);
        igraph_vector_t ig_outdeg = igVectorView(outdeg);
        igraph_degseq_t ig_method;
        switch (method) {
        case 0: ig_method = IGRAPH_DEGSEQ_SIMPLE; break;
        case 1: ig_method = IGRAPH_DEGSEQ_SIMPLE_NO_MULTIPLE; break;
        case 2: ig_method = IGRAPH_DEGSEQ_VL; break;
        case 3: ig_method = IGRAPH_DEGSEQ_SIMPLE_NO_MULTIPLE_UNIFORM; break;
        default: throw mma::LibraryError("degreeSequenceGame: unknown method option.");
        }

        destroy();
        int err;
        if (indeg.length() == 0)
            err = igraph_degree_sequence_game(&graph, &ig_outdeg, nullptr, ig_method);
        else
            err = igraph_degree_sequence_game(&graph, &ig_outdeg, &ig_indeg, ig_method);
        igConstructorCheck(err);
    }

    void kRegularGame(mint n, mint k, bool directed, bool multiple) {
        destroy();
        igConstructorCheck(igraph_k_regular_game(&graph, n, k, directed, multiple));
    }

    void stochasticBlockModel(mma::RealMatrixRef tmat, mma::IntTensorRef tsizes, bool directed, bool loops) {
        igIntVector sizes;
        sizes.copyFromMTensor(tsizes);
        igMatrix mat;
        mat.copyFromMTensor(tmat);
        igraph_integer_t n = 0;
        for (const auto &el : sizes)
            n += el;
        destroy();
        igConstructorCheck(igraph_sbm_game(&graph, n, &mat.mat, &sizes.vec, directed, loops));
    }

    void forestFireGame(mint n, double fwprob, double bwratio, mint nambs, bool directed) {
        destroy();
        igConstructorCheck(igraph_forest_fire_game(&graph, n, fwprob, bwratio, nambs, directed));
    }

    void bipartiteGameGNM(mint n1, mint n2, mint m, bool directed, bool bidirectional) {
        destroy();
        igConstructorCheck(igraph_bipartite_game(
                               &graph, nullptr, IGRAPH_ERDOS_RENYI_GNM,
                               n1, n2, 0, m, directed, bidirectional ? IGRAPH_ALL : IGRAPH_OUT));
    }

    void bipartiteGameGNP(mint n1, mint n2, double p, bool directed, bool bidirectional) {
        destroy();
        igConstructorCheck(igraph_bipartite_game(
                               &graph, nullptr, IGRAPH_ERDOS_RENYI_GNP,
                               n1, n2, p, 0, directed, bidirectional ? IGRAPH_ALL : IGRAPH_OUT));
    }

    void erdosRenyiGNM(mint n, mint m, bool directed, bool loops) {
        destroy();
        igConstructorCheck(igraph_erdos_renyi_game_gnm(&graph, n, m, directed, loops));
    }

    void erdosRenyiGNP(mint n, double p, bool directed, bool loops) {
        destroy();
        igConstructorCheck(igraph_erdos_renyi_game_gnp(&graph, n, p, directed, loops));
    }

    mma::RealTensorRef geometricGame(mint n, double radius, bool periodic) {        
        destroy();
        igVector x, y;
        igConstructorCheck(igraph_grg_game(
                               &graph, n, radius, periodic,
                               &x.vec, &y.vec));

        const mint len = x.length();
        mma::RealMatrixRef coord = mma::makeMatrix<double>(len, 2);
        auto xp = x.begin();
        auto yp = y.begin();
        for (mint i=0; i < len; ++i) {
            coord(i, 0) = *xp++;
            coord(i, 1) = *yp++;
        }
        return coord;
    }

    void barabasiAlbertGame(mint n, double power, double A, mint m, mma::RealTensorRef mtens, bool directed, bool totalDegree, mint method) {
        destroy();
        igraph_vector_t mvec = igVectorView(mtens);
        igraph_barabasi_algorithm_t algo;
        switch (method) {
        case 0: algo = IGRAPH_BARABASI_BAG; break;
        case 1: algo = IGRAPH_BARABASI_PSUMTREE; break;
        case 2: algo = IGRAPH_BARABASI_PSUMTREE_MULTIPLE; break;
        default:
            empty();
            throw mma::LibraryError("Unknown method for Barabasi-Albert game.");
        }

        igConstructorCheck(igraph_barabasi_game(&graph, n, power, m, &mvec, totalDegree, A, directed, algo, nullptr));
    }

    void barabasiAlbertGameWithStartingGraph(mint n, double power, double A, mint m, mma::RealTensorRef mtens, bool directed, bool totalDegree, mint method, IG &start) {
        destroy();
        igraph_vector_t mvec = igVectorView(mtens);
        igraph_barabasi_algorithm_t algo;
        switch (method) {
        case 0: algo = IGRAPH_BARABASI_BAG; break;
        case 1: algo = IGRAPH_BARABASI_PSUMTREE; break;
        case 2: algo = IGRAPH_BARABASI_PSUMTREE_MULTIPLE; break;
        default:
            empty();
            throw mma::LibraryError("Unknown method for Barabasi-Albert game.");
        }
        // inherit directnedness from starting graph if it is not empty
        if (start.edgeCount() > 0)
            directed = start.directedQ();
        else if (directed)
            igraph_to_directed(&start.graph, IGRAPH_TO_DIRECTED_ARBITRARY);
        igConstructorCheck(igraph_barabasi_game(&graph, n, power, m, &mvec, totalDegree, A, directed, algo, &start.graph));
    }

    void wattsStrogatzGame(mint dim, mint size, mint radius, double p, bool loops, bool multiple) {
        destroy();
        igConstructorCheck(igraph_watts_strogatz_game(&graph, dim, size, radius, p, loops, multiple));
    }

    void staticFitnessGame(mint m /* edges */, mma::RealTensorRef fit_in_ten, mma::RealTensorRef fit_out_ten, bool loops, bool multiple) {
        destroy();
        igraph_vector_t fit_in = igVectorView(fit_in_ten);
        igraph_vector_t fit_out = igVectorView(fit_out_ten);
        igConstructorCheck(igraph_static_fitness_game(
                               &graph, m,
                               &fit_in, fit_out_ten.length() == 0 ? nullptr : &fit_out,
                               loops, multiple));
    }

    void staticPowerLawGame(mint n, mint m, double exp_out, double exp_in, bool loops, bool multiple, bool finite_size_correction) {
        destroy();
        igConstructorCheck(igraph_static_power_law_game(&graph, n, m, exp_out, exp_in, loops, multiple, finite_size_correction));
    }

    void growingGame(mint n, mint m, bool directed, bool citation) {
        destroy();
        igConstructorCheck(igraph_growing_random_game(&graph, n, m, directed, citation));
    }

    void callawayTraitsGame(mint n, mint edge_per_step, mma::RealTensorRef type_distr_ten, mma::RealMatrixRef pref_matrix_ten, bool directed) {
        destroy();
        igraph_vector_t type_distr = igVectorView(type_distr_ten);
        igMatrix pref_matrix;
        pref_matrix.copyFromMTensor(pref_matrix_ten);
        igConstructorCheck(igraph_callaway_traits_game(&graph, n, type_distr_ten.size(), edge_per_step, &type_distr, &pref_matrix.mat, directed, nullptr));
    }

    void establishmentGame(mint n, mint k, mma::RealTensorRef type_distr_ten, mma::RealMatrixRef pref_matrix_ten, bool directed) {
        destroy();
        igraph_vector_t type_distr = igVectorView(type_distr_ten);
        igMatrix pref_matrix;
        pref_matrix.copyFromMTensor(pref_matrix_ten);
        igConstructorCheck(igraph_establishment_game(&graph, n, type_distr_ten.size(), k, &type_distr, &pref_matrix.mat, directed, nullptr));
    }

    void preferenceGame(mint n, mma::RealTensorRef type_distr_ten, mma::RealMatrixRef pref_matrix_ten, bool directed, bool loops) {
        destroy();
        igraph_vector_t type_distr = igVectorView(type_distr_ten);
        igMatrix pref_matrix;
        pref_matrix.copyFromMTensor(pref_matrix_ten);
        igConstructorCheck(igraph_preference_game(&graph, n, type_distr_ten.size(), &type_distr, false, &pref_matrix.mat, nullptr, directed, loops));
    }

    void asymmetricPreferenceGame(mint n, mma::RealMatrixRef type_distr_ten, mma::RealMatrixRef pref_matrix_ten, bool loops) {
        destroy();
        igMatrix type_distr;
        type_distr.copyFromMTensor(type_distr_ten);
        igMatrix pref_matrix;
        pref_matrix.copyFromMTensor(pref_matrix_ten);
        igConstructorCheck(igraph_asymmetric_preference_game(&graph, n, type_distr.nrow(), type_distr.ncol(), &type_distr.mat, &pref_matrix.mat, nullptr, nullptr, loops));
    }

    // Modification

    void connectNeighborhood(mint order) {
        igCheck(igraph_connect_neighborhood(&graph, order, IGRAPH_OUT));
    }

    void mycielski();

    // Structure

    mma::RealTensorRef edgeList() const {
        igVector vec;
        igCheck(igraph_get_edgelist(&graph, &vec.vec, false));
        mma::RealMatrixRef res = mma::makeMatrix<double>(vec.length() / 2, 2);
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

    bool multiQ() const {
        igraph_bool_t res;
        igCheck(igraph_has_multiple(&graph, &res));
        return res;
    }

    bool connectedQ(bool strong) const {
        igraph_bool_t res;
        igCheck(igraph_is_connected(&graph, &res, strong ? IGRAPH_STRONG : IGRAPH_WEAK));
        return res;
    }

    bool treeQ(mint mode) const {
        igraph_bool_t res;
        igCheck(igraph_is_tree(&graph, &res, nullptr, igNeighborMode(mode, "treeQ")));
        return res;
    }

    bool forestQ(mint mode) const;

    bool bipartiteQ() const {
        igraph_bool_t res;
        igCheck(igraph_is_bipartite(&graph, &res, nullptr));
        return res;
    }

    // Centrality measures

    mma::RealTensorRef betweenness(bool normalized, mma::RealTensorRef vs) const {
        igVector res;
        igraph_vector_t vsvec = igVectorView(vs);
        igCheck(igraph_betweenness(
                &graph, &res.vec, vs.length() == 0 ? igraph_vss_all() : igraph_vss_vector(&vsvec),
                true, passWeights()));
        auto t = res.makeMTensor();
        if (normalized) {
            double vcount = vertexCount();
            double norm = 1.0;
            if (vcount > 2) {
                if (! directedQ())
                    norm = 2.0;
                norm = norm / ((vcount-1)*(vcount-2));
            }

            for (auto &el : t)
                el *= norm;
        }
        return t;
    }

    mma::RealTensorRef subsetBetweenness(bool normalized, mma::RealTensorRef vs, mma::RealTensorRef sources, mma::RealTensorRef targets) const {
        igVector res;
        igraph_vector_t vsvec = igVectorView(vs);
        igraph_vector_t s = igVectorView(sources);
        igraph_vector_t t = igVectorView(targets);
        igCheck(
            igraph_betweenness_subset(
                        &graph, &res.vec,
                        vs.length() == 0 ? igraph_vss_all() : igraph_vss_vector(&vsvec),
                        true,
                        igraph_vss_vector(&s), igraph_vss_vector(&t),
                        passWeights())
                    );
        auto result = res.makeMTensor();
        if (normalized) {
            double scount = sources.size();
            double tcount = targets.size();
            double norm = 1.0;
            if (vcount > 2) {
                if (! directedQ())
                    norm = 2.0;
                norm = norm / ((vcount-1)*(vcount-2));
            }

            for (auto &el : result)
                el *= norm;
        }
        return result;
    }

    mma::RealTensorRef edgeBetweenness(bool normalized) const {
        igVector res;
        igCheck(igraph_edge_betweenness(&graph, &res.vec, true, passWeights()));
        auto t = res.makeMTensor();
        if (normalized) {
            double vcount = vertexCount();
            double norm = 1.0;
            if (vcount > 2) {
                if (! directedQ())
                    norm = 2.0;
                norm = norm / (vcount*(vcount-1));
            }

            for (auto &el : t)
                el *= norm;
        }
        return t;
    }

    mma::RealTensorRef subsetEdgeBetweenness(bool normalized, mma::RealTensorRef sources, mma::RealTensorRef targets) const {
        igVector res;
        igraph_vector_t s = igVectorView(sources);
        igraph_vector_t t = igVectorView(targets);
        igCheck(
            igraph_edge_betweenness_subset(
                        &graph, &res.vec,
                        igraph_ess_all(IGRAPH_EDGEORDER_ID),
                        true,
                        igraph_vss_vector(&s), igraph_vss_vector(&t),
                        passWeights())
                    );
        return res.makeMTensor();
    }

    mma::RealTensorRef closeness(bool normalized, mma::RealTensorRef vs) const {
        igVector res;
        igraph_vector_t vsvec = igVectorView(vs);
        igCheck(igraph_closeness(
                &graph, &res.vec, nullptr, nullptr, vs.length() == 0 ? igraph_vss_all() : igraph_vss_vector(&vsvec),
                IGRAPH_OUT, passWeights(), normalized));
        return res.makeMTensor();
    }

    mma::RealTensorRef harmonicCentrality(bool normalized, mma::RealTensorRef vs) const {
        igVector res;
        igraph_vector_t vsvec = igVectorView(vs);
        igCheck(igraph_harmonic_centrality(
                &graph, &res.vec, vs.length() == 0 ? igraph_vss_all() : igraph_vss_vector(&vsvec),
                IGRAPH_OUT, passWeights(), normalized));
        return res.makeMTensor();
    }

    mma::RealTensorRef personalizedPageRank(mint method, mma::RealTensorRef treset, double damping, bool directed) const {
        igraph_pagerank_algo_t algo;
        igraph_arpack_options_t arpack_options;
        switch (method) {
        case 0:
            throw mma::LibraryError("Power iteration is no longer supported for PageRank calculation.");
        case 1:
            algo = IGRAPH_PAGERANK_ALGO_ARPACK;
            igraph_arpack_options_init(&arpack_options);
            break;
        case 2:
            algo = IGRAPH_PAGERANK_ALGO_PRPACK;
            break;
        default: throw mma::LibraryError("Unknown PageRank method.");
        }

        igraph_vector_t reset = igVectorView(treset);

        igVector vector;
        double value;

        igCheck(igraph_personalized_pagerank(
                    &graph, algo, &vector.vec, &value, igraph_vss_all(),
                    directed, damping, treset.size() == 0 ? nullptr : &reset, passWeights(), &arpack_options));

        // 7 binary digits of tolerance, i.e. a factor of 128, if what Mathematica's Equal uses
        // With igraph's vendored ARPACK, 5 digits were sufficient in a few tests,
        // but with an external ARPACK 3.8.0, 7 digits were needed to avoid false positives.
        if (std::fabs(value - 1.0) > 128*std::numeric_limits<double>::epsilon()) {
            std::ostringstream msg;
            msg << std::setprecision(std::numeric_limits<double>::digits10) <<
                   "The PageRank eigenvalue should be 1, actual eigenvalue is " << value << ". Possible convergence problem.";
            mma::message(msg.str(), mma::M_WARNING);
        }

        return vector.makeMTensor();
    }

    mma::RealTensorRef personalizedLinkRank(mint method, mma::RealTensorRef treset, double damping, bool directed) const {
        igraph_pagerank_algo_t algo;
        igraph_arpack_options_t arpack_options;
        switch (method) {
        case 0:
            throw mma::LibraryError("Power iteration is no longer supported for LinkRank calculation.");
        case 1:
            algo = IGRAPH_PAGERANK_ALGO_ARPACK;
            igraph_arpack_options_init(&arpack_options);
            break;
        case 2:
            algo = IGRAPH_PAGERANK_ALGO_PRPACK;
            break;
        default: throw mma::LibraryError("Unknown LinkRank method.");
        }

        igraph_vector_t reset = igVectorView(treset);

        igVector page_rank;
        double value;

        igCheck(igraph_personalized_pagerank(
                    &graph, algo, &page_rank.vec, &value, igraph_vss_all(),
                    directed, damping, treset.size() == 0 ? nullptr : &reset, passWeights(), &arpack_options));

        // 7 binary digits of tolerance, i.e. a factor of 128, if what Mathematica's Equal uses
        // With igraph's vendored ARPACK, 5 digits were sufficient in a few tests,
        // but with an external ARPACK 3.8.0, 7 digits were needed to avoid false positives.
        if (std::fabs(value - 1.0) > 128*std::numeric_limits<double>::epsilon()) {
            std::ostringstream msg;
            msg << std::setprecision(std::numeric_limits<double>::digits10) <<
                   "The LinkRank eigenvalue should be 1, actual eigenvalue is " << value << ". Possible convergence problem.";
            mma::message(msg.str(), mma::M_WARNING);
        }

        auto link_rank = mma::makeVector<double>(edgeCount());
        std::fill(link_rank.begin(), link_rank.end(), 0.0);

        igraph_neimode_t mode = directed ? IGRAPH_OUT : IGRAPH_ALL;

        igVector strength;
        igCheck(igraph_strength(&graph, &strength.vec, igraph_vss_all(), mode, IGRAPH_LOOPS, passWeights()));

        igVector out_edges;
        mint vcount = vertexCount();
        for (mint v=0; v < vcount; ++v) {
            igCheck(igraph_incident(&graph, &out_edges.vec, v, mode));

            for (const auto &e : out_edges) {
                double weight = weightedQ() ? weights[e] : 1.0;
                link_rank[e] += strength[v] == 0 ? 0.0 : page_rank[v] * weight / strength[v];
            }
        }

        return link_rank;
    }

    mma::RealTensorRef eigenvectorCentrality(bool directed, bool normalized) const {
        igVector vector;
        double value;
        igraph_arpack_options_t options;
        igraph_arpack_options_init(&options);
        igCheck(igraph_eigenvector_centrality(&graph, &vector.vec, &value, directed, normalized, passWeights(), &options));
        return vector.makeMTensor();
    }

    mma::RealTensorRef hubScore(bool normalized) const {
        igVector vector;
        double value;
        igraph_arpack_options_t options;
        igraph_arpack_options_init(&options);
        igCheck(igraph_hub_score(&graph, &vector.vec, &value, normalized, passWeights(), &options));
        if (value < -std::pow( std::numeric_limits<double>::epsilon(), 2.0/3 ) ) {
            std::ostringstream msg;
            msg << "Hub score eigenvalue should be positive, actual eigenvalue is " << value << ". Possible convergence problem.";
            mma::message(msg.str(), mma::M_WARNING);
        }
        return vector.makeMTensor();
    }

    mma::RealTensorRef authorityScore(bool normalized) const {
        igVector vector;
        double value;
        igraph_arpack_options_t options;
        igraph_arpack_options_init(&options);
        igCheck(igraph_authority_score(&graph, &vector.vec, &value, normalized, passWeights(), &options));
        if (value < -std::pow( std::numeric_limits<double>::epsilon(), 2.0/3 ) ) {
            std::ostringstream msg;
            msg << "Authority score eigenvalue should be positive, actual eigenvalue is " << value << ". Possible convergence problem.";
            mma::message(msg.str(), mma::M_WARNING);
        }
        return vector.makeMTensor();
    }

    mma::RealTensorRef constraintScore() const {
        igVector vec;
        igCheck(igraph_constraint(&graph, &vec.vec, igraph_vss_all(), passWeights()));
        return vec.makeMTensor();
    }

    // Range-limited centrality measures

    mma::RealTensorRef betweennessCutoff(double cutoff, bool normalized, mma::RealTensorRef vs) const {
        igVector res;
        igraph_vector_t vsvec = igVectorView(vs);
        igCheck(igraph_betweenness_cutoff(&graph, &res.vec, vs.length() == 0 ? igraph_vss_all() : igraph_vss_vector(&vsvec), true, passWeights(), cutoff));
        auto t = res.makeMTensor();
        if (normalized) {
            double vcount = vertexCount();
            double norm = 1.0;
            if (vcount > 2) {
                if (! directedQ())
                    norm = 2.0;
                norm = norm / ((vcount-1)*(vcount-2));
            }

            for (auto &el : t)
                el *= norm;
        }
        return t;
    }

    mma::RealTensorRef edgeBetweennessCutoff(double cutoff, bool normalized) const {
        igVector res;
        igCheck(igraph_edge_betweenness_cutoff(&graph, &res.vec, true, passWeights(), cutoff));
        auto t = res.makeMTensor();
        if (normalized) {
            double vcount = vertexCount();
            double norm = 1.0;
            if (vcount > 2) {
                if (! directedQ())
                    norm = 2.0;
                norm = norm / (vcount*(vcount-1));
            }

            for (auto &el : t)
                el *= norm;
        }
        return t;
    }

    mma::RealTensorRef closenessCutoff(double cutoff, bool normalized, mma::RealTensorRef vs) const {
        igVector res;
        igraph_vector_t vsvec = igVectorView(vs);
        igCheck(igraph_closeness_cutoff(&graph, &res.vec, nullptr, nullptr, vs.length() == 0 ? igraph_vss_all() : igraph_vss_vector(&vsvec), IGRAPH_OUT, passWeights(), normalized, cutoff));
        return res.makeMTensor();
    }

    mma::RealTensorRef neighborhoodCloseness(double cutoff, bool normalized, mma::RealTensorRef vs) const {
        igVector res, reachable;
        igraph_vector_t vsvec = igVectorView(vs);
        igCheck(igraph_closeness_cutoff(&graph, &res.vec, &reachable.vec, nullptr, vs.length() == 0 ? igraph_vss_all() : igraph_vss_vector(&vsvec), IGRAPH_OUT, passWeights(), normalized, cutoff));
        auto ret = mma::makeMatrix<double>(2, res.size());
        std::copy(res.begin(), res.end(), ret.data());
        std::copy(reachable.begin(), reachable.end(), ret.data() + res.size());
        return ret;
    }

    mma::RealTensorRef harmonicCentralityCutoff(double cutoff, bool normalized, mma::RealTensorRef vs) const {
        igVector res;
        igraph_vector_t vsvec = igVectorView(vs);
        igCheck(igraph_harmonic_centrality_cutoff(&graph, &res.vec, vs.length() == 0 ? igraph_vss_all() : igraph_vss_vector(&vsvec), IGRAPH_OUT, passWeights(), normalized, cutoff));
        return res.makeMTensor();
    }


    // Centralization

    double degreeCentralization(mint mode, bool loops, bool normalized) const {
        igraph_real_t result;
        igCheck(igraph_centralization_degree(&graph, nullptr, igNeighborMode(mode, "degree centralization"), loops, &result, nullptr, normalized));
        return result;
    }

    double betweennessCentralization(bool normalized) const {
        igraph_real_t result;
        igCheck(igraph_centralization_betweenness(&graph, nullptr, true, &result, nullptr, normalized));
        return result;
    }

    double closenessCentralization(bool normalized) const {
        igraph_real_t result;
        igCheck(igraph_centralization_closeness(&graph, nullptr, IGRAPH_OUT, &result, nullptr, normalized));
        return result;
    }

    double eigenvectorCentralization(bool scale, bool normalized) const {
        igraph_real_t result;
        igraph_arpack_options_t options;
        igraph_arpack_options_init(&options);
        igCheck(igraph_centralization_eigenvector_centrality(&graph, nullptr, nullptr, true, scale, &options, &result, nullptr, normalized));
        return result;
    }

    static double centralization(mma::RealTensorRef vec, double tmax, bool normalized) {
        igraph_vector_t scores = igVectorView(vec);
        return igraph_centralization(&scores, tmax, normalized);
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

    void rewireDirectedEdges(double prob, bool loops, bool out) {
        igCheck(igraph_rewire_directed_edges(&graph, prob, loops, out ? IGRAPH_OUT : IGRAPH_IN));
    }

    // Isomorphism (general)

    bool isomorphic(IG &ig) {
        igraph_bool_t res;
        emptyMatchDirectedness(ig);
        igCheck(igraph_isomorphic(&graph, &ig.graph, &res));
        return res;
    }

    bool subisomorphic(IG &ig) {
        igraph_bool_t res;
        emptyMatchDirectedness(ig);
        igCheck(igraph_subisomorphic(&graph, &ig.graph, &res));
        return res;
    }

    mma::RealTensorRef getIsomorphism(IG &ig);

    mma::RealTensorRef getSubisomorphism(IG &ig);

    mint isoclass() const {
        igraph_integer_t res;
        igCheck(igraph_isoclass(&graph, &res));
        return res;
    }

    // Isomorphism (bliss)

    mma::RealTensorRef blissCanonicalPermutation(mint splitting, mma::IntTensorRef col) const {
        igVector vec;
        igIntVector colvec;
        colvec.copyFromMTensor(col);
        igCheck(igraph_canonical_permutation(&graph, col.length() == 0 ? nullptr : &colvec.vec, &vec.vec, blissIntToSplitting(splitting), nullptr));
        return vec.makeMTensor();
    }

    bool blissIsomorphic(IG &ig, mint splitting, mma::IntTensorRef col1, mma::IntTensorRef col2) {
        igraph_bool_t res;
        igIntVector colvec1, colvec2;
        colvec1.copyFromMTensor(col1);
        colvec2.copyFromMTensor(col2);
        emptyMatchDirectedness(ig);
        igCheck(igraph_isomorphic_bliss(
                    &graph, &ig.graph, col1.length() == 0 ? nullptr : &colvec1.vec, col2.length() == 0 ? nullptr : &colvec2.vec,
                    &res, nullptr, nullptr, blissIntToSplitting(splitting), nullptr, nullptr));
        return res;
    }

    mma::RealTensorRef blissFindIsomorphism(IG &ig, mint splitting, mma::IntTensorRef col1, mma::IntTensorRef col2) {
        igraph_bool_t res;
        igIntVector colvec1, colvec2;
        colvec1.copyFromMTensor(col1);
        colvec2.copyFromMTensor(col2);
        emptyMatchDirectedness(ig);
        igVector map;
        igCheck(igraph_isomorphic_bliss(
                    &graph, &ig.graph, col1.length() == 0 ? nullptr : &colvec1.vec, col2.length() == 0 ? nullptr : &colvec2.vec,
                    &res, &map.vec, nullptr, blissIntToSplitting(splitting), nullptr, nullptr));
        if (res)
            return map.makeMTensor();
        else
            return mma::makeVector<double>(0);
    }

    void blissAutomorphismCount(MLINK link) const {
        igraph_bliss_info_t info;
        mlStream ml{link, "blissAutomorphismsCount"};
        int splitting;
        igIntVector colors;
        ml >> mlCheckArgs(2) >> splitting >> colors;

        igCheck(igraph_automorphisms(&graph, colors.length() == 0 ? nullptr : &colors.vec, blissIntToSplitting(splitting), &info));

        ml.newPacket();
        ml << info.group_size;
        igraph_free(info.group_size);
    }

    // returns the generators of the group
    void blissAutomorphismGroup(MLINK link) const {
        mlStream ml{link, "blissAutomorphismGroup"};
        int splitting;
        igIntVector colors;
        ml >> mlCheckArgs(2) >> splitting >> colors;

        igList list;
        igCheck(igraph_automorphism_group(&graph, colors.length() == 0 ? nullptr : &colors.vec, &list.list, blissIntToSplitting(splitting), nullptr));

        ml.newPacket();
        ml << list;
    }

    // Isomorphism (VF2)

    bool vf2Isomorphic(
            IG &ig, mma::IntTensorRef vcol1, mma::IntTensorRef vcol2,
                    mma::IntTensorRef ecol1, mma::IntTensorRef ecol2)
    {
        igIntVector vc1; vc1.copyFromMTensor(vcol1);
        igIntVector vc2; vc2.copyFromMTensor(vcol2);
        igIntVector ec1; ec1.copyFromMTensor(ecol1);
        igIntVector ec2; ec2.copyFromMTensor(ecol2);

        emptyMatchDirectedness(ig);

        igraph_bool_t res;
        igCheck(igraph_isomorphic_vf2(
                    &graph, &ig.graph,
                    vcol1.length() == 0 ? nullptr : &vc1.vec, vcol2.length() == 0 ? nullptr : &vc2.vec,
                    ecol1.length() == 0 ? nullptr : &ec1.vec, ecol2.length() == 0 ? nullptr : &ec2.vec,
                    &res, nullptr, nullptr, nullptr, nullptr, nullptr));
        return res;
    }

    void vf2FindIsomorphisms(MLINK link);

    bool vf2Subisomorphic(
            IG &ig, mma::IntTensorRef vcol1, mma::IntTensorRef vcol2,
                    mma::IntTensorRef ecol1, mma::IntTensorRef ecol2)
    {
        igIntVector vc1; vc1.copyFromMTensor(vcol1);
        igIntVector vc2; vc2.copyFromMTensor(vcol2);
        igIntVector ec1; ec1.copyFromMTensor(ecol1);
        igIntVector ec2; ec2.copyFromMTensor(ecol2);

        emptyMatchDirectedness(ig);

        igraph_bool_t res;
        igCheck(igraph_subisomorphic_vf2(
                    &graph, &ig.graph,
                    vcol1.length() == 0 ? nullptr : &vc1.vec, vcol2.length() == 0 ? nullptr : &vc2.vec,
                    ecol1.length() == 0 ? nullptr : &ec1.vec, ecol2.length() == 0 ? nullptr : &ec2.vec,
                    &res, nullptr, nullptr, nullptr, nullptr, nullptr));
        return res;
    }

    void vf2FindSubisomorphisms(MLINK link);

    mint vf2IsomorphismCount(
            IG &ig, mma::IntTensorRef vcol1, mma::IntTensorRef vcol2,
                    mma::IntTensorRef ecol1, mma::IntTensorRef ecol2)
    {
        igIntVector vc1; vc1.copyFromMTensor(vcol1);
        igIntVector vc2; vc2.copyFromMTensor(vcol2);
        igIntVector ec1; ec1.copyFromMTensor(ecol1);
        igIntVector ec2; ec2.copyFromMTensor(ecol2);

        emptyMatchDirectedness(ig);

        igraph_integer_t res;
        igCheck(igraph_count_isomorphisms_vf2(
                    &graph, &ig.graph,
                    vcol1.length() == 0 ? nullptr : &vc1.vec, vcol2.length() == 0 ? nullptr : &vc2.vec,
                    ecol1.length() == 0 ? nullptr : &ec1.vec, ecol2.length() == 0 ? nullptr : &ec2.vec,
                    &res, nullptr, nullptr, nullptr));
        return res;
    }

    mint vf2SubisomorphismCount(
            IG &ig, mma::IntTensorRef vcol1, mma::IntTensorRef vcol2,
                    mma::IntTensorRef ecol1, mma::IntTensorRef ecol2)
    {
        igIntVector vc1; vc1.copyFromMTensor(vcol1);
        igIntVector vc2; vc2.copyFromMTensor(vcol2);
        igIntVector ec1; ec1.copyFromMTensor(ecol1);
        igIntVector ec2; ec2.copyFromMTensor(ecol2);

        emptyMatchDirectedness(ig);

        igraph_integer_t res;
        igCheck(igraph_count_subisomorphisms_vf2(
                    &graph, &ig.graph,
                    vcol1.length() == 0 ? nullptr : &vc1.vec, vcol2.length() == 0 ? nullptr : &vc2.vec,
                    ecol1.length() == 0 ? nullptr : &ec1.vec, ecol2.length() == 0 ? nullptr : &ec2.vec,
                    &res, nullptr, nullptr, nullptr));
        return res;
    }

    // this function is used in the implementation of others
    // vc and ec are output arguments
    // this is a constructor!
    void createColoredSimpleGraph(const IG &ig, igIntVector &vc, igIntVector &ec) {
        destroy();
        igConstructorCheck(igraph_simplify_and_colorize(&ig.graph, &graph, &vc.vec, &ec.vec));
    }

    bool vf2IsomorphicMulti(IG &ig);

    bool vf2SubisomorphicMulti(IG &ig);

    // Isomorphism (LAD)

    bool ladSubisomorphic(IG &ig, bool induced) {
        igraph_bool_t res;
        emptyMatchDirectedness(ig);
        igCheck(igraph_subisomorphic_lad(&ig.graph, &graph, nullptr, &res, nullptr, nullptr, induced, 0));
        return res;
    }

    void ladSubisomorphicColored(MLINK link) {
        mlStream ml{link, "ladSubisomorphicColored"};
        mint id;
        igraph_bool_t induced;
        igList domain;
        ml >> mlCheckArgs(3)
           >> id >> induced >> domain;

        IG &ig = mma::getInstance<IG>(id);
        emptyMatchDirectedness(ig);

        igraph_bool_t res;
        igCheck(igraph_subisomorphic_lad(&ig.graph, &graph, &domain.list, &res, nullptr, nullptr, induced, 0));

        ml.newPacket();
        if (res)
            ml << mlSymbol("True");
        else
            ml << mlSymbol("False");
    }

    mma::RealTensorRef ladGetSubisomorphism(IG &ig, bool induced) {
        igraph_bool_t iso;
        igVector map;
        emptyMatchDirectedness(ig);
        igCheck(igraph_subisomorphic_lad(&ig.graph, &graph, nullptr, &iso, &map.vec, nullptr, induced, 0));
        return map.makeMTensor();
    }

    void ladGetSubisomorphismColored(MLINK link) {
        mlStream ml{link, "ladGetSubisomorphismColored"};
        mint id;
        igraph_bool_t induced;
        igList domain;
        ml >> mlCheckArgs(3)
           >> id >> induced >> domain;

        IG &ig = mma::getInstance<IG>(id);
        emptyMatchDirectedness(ig);

        igraph_bool_t iso;
        igVector map;
        igCheck(igraph_subisomorphic_lad(&ig.graph, &graph, &domain.list, &iso, &map.vec, nullptr, induced, 0));

        ml.newPacket();
        ml << map;
    }

    void ladFindSubisomorphisms(MLINK link) {
        mlStream ml{link, "ladFindSubisomorphism"};
        mint id;
        igraph_bool_t induced;
        igList domain;
        ml >> mlCheckArgs(3) >> id >> induced >> domain;

        IG &ig = mma::getInstance<IG>(id);
        emptyMatchDirectedness(ig);

        igList list;
        igraph_bool_t iso;
        igCheck(igraph_subisomorphic_lad(&ig.graph, &graph, domain.length() == 0 ? nullptr : &domain.list, &iso, nullptr, &list.list, induced, 0));

        ml.newPacket();
        ml << list;
    }

    mint ladCountSubisomorphisms(IG &ig, bool induced) {
        // special case: one match for the null pattern for consistency
        if (ig.vertexCount() == 0)
            return 1;

        igraph_bool_t iso;
        igList list;
        emptyMatchDirectedness(ig);
        igCheck(igraph_subisomorphic_lad(&ig.graph, &graph, nullptr, &iso, nullptr, &list.list, induced, 0));
        return list.length();
    }

    void ladCountSubisomorphismsColored(MLINK link) {
        mlStream ml{link, "ladCountSubisomorphismsColored"};
        mint id;
        igraph_bool_t induced;
        igList domain;
        ml >> mlCheckArgs(3) >> id >> induced >> domain;

        IG &ig = mma::getInstance<IG>(id);

        // special case: one match for the null pattern for consistency
        if (ig.vertexCount() == 0) {
            ml.newPacket();
            ml << 1;
            return;
        }

        emptyMatchDirectedness(ig);

        igList list;
        igraph_bool_t iso;
        igCheck(igraph_subisomorphic_lad(&ig.graph, &graph, domain.length() == 0 ? nullptr : &domain.list, &iso, nullptr, &list.list, induced, 0));

        ml.newPacket();
        ml << list.length();
    }

    // Isomorphism-related tests

    // the input is expected to be a simple graph
    bool selfComplementaryQ() const {
        IG comp;
        comp.complement(*this);

        igraph_bool_t iso;
        igCheck(igraph_isomorphic(&graph, &comp.graph, &iso));

        return iso;
    }

    bool vertexTransitiveQ(mint splitting) const;
    bool edgeTransitiveQ(mint splitting) const;
    bool symmetricQ(mint splitting) const;
    bool distanceTransitiveQ(mint splitting) const;

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

    mma::RealMatrixRef motifsParticipation(mint size, mma::RealTensorRef cut_prob) const {
        igraph_vector_t ig_cut_prob = igVectorView(cut_prob);

        struct MotifCounter {
            static igraph_bool_t handler(const igraph_t *, igraph_vector_t *vids, int isoclass, void *data) {
                auto *mat = static_cast<mma::RealMatrixRef *>(data);
                long len = igraph_vector_size(vids);
                for (int i=0; i < len; ++i)
                    (*mat)(VECTOR(*vids)[i], isoclass) += 1;
                return false;
            }

            static int isoclass_count(bool directed, int vcount) {
                switch (vcount) {
                case 3: return directed ? 16 : 4;
                case 4: return directed ? 218 : 11;
                case 5:
                    if (! directed) return 34;
                case 6:
                    if (! directed) return 156;
                default:
                    throw mma::LibraryError("Unsupported motif size.");
                }
            }
        };

        bool dir = directedQ();

        mma::RealMatrixRef mat = mma::makeMatrix<double>(vertexCount(), MotifCounter::isoclass_count(dir, size));
        for (auto &el : mat) el = 0;

        LTGuard<mma::RealTensorRef> guard(mat); // automatically free 'mat' upon premature exit from the function

        igCheck(igraph_motifs_randesu_callback(&graph, size, &ig_cut_prob, &(MotifCounter::handler), &mat));

        for (mint i=0; i < mat.rows(); ++i)
            switch (size) {
            case 3:
                if (dir) {
                    mat(i,0) = IGRAPH_NAN;
                    mat(i,1) = IGRAPH_NAN;
                    mat(i,3) = IGRAPH_NAN;
                } else {
                    mat(i,0) = IGRAPH_NAN;
                    mat(i,1) = IGRAPH_NAN;
                }
                break;
            case 4:
                if (dir) {
                    mat(i,0) = IGRAPH_NAN;
                    mat(i,1) = IGRAPH_NAN;
                    mat(i,2) = IGRAPH_NAN;
                    mat(i,4) = IGRAPH_NAN;
                    mat(i,5) = IGRAPH_NAN;
                    mat(i,6) = IGRAPH_NAN;
                    mat(i,9) = IGRAPH_NAN;
                    mat(i,10) = IGRAPH_NAN;
                    mat(i,11) = IGRAPH_NAN;
                    mat(i,15) = IGRAPH_NAN;
                    mat(i,22) = IGRAPH_NAN;
                    mat(i,23) = IGRAPH_NAN;
                    mat(i,27) = IGRAPH_NAN;
                    mat(i,28) = IGRAPH_NAN;
                    mat(i,33) = IGRAPH_NAN;
                    mat(i,34) = IGRAPH_NAN;
                    mat(i,39) = IGRAPH_NAN;
                    mat(i,62) = IGRAPH_NAN;
                    mat(i,120) = IGRAPH_NAN;
                } else {
                    mat(i,0) = IGRAPH_NAN;
                    mat(i,1) = IGRAPH_NAN;
                    mat(i,2) = IGRAPH_NAN;
                    mat(i,3) = IGRAPH_NAN;
                    mat(i,5) = IGRAPH_NAN;
                }
                break;
            case 5:
                if (! dir) {
                    mat(i,0) = IGRAPH_NAN;
                    mat(i,1) = IGRAPH_NAN;
                    mat(i,2) = IGRAPH_NAN;
                    mat(i,3) = IGRAPH_NAN;
                    mat(i,4) = IGRAPH_NAN;
                    mat(i,5) = IGRAPH_NAN;
                    mat(i,6) = IGRAPH_NAN;
                    mat(i,7) = IGRAPH_NAN;
                    mat(i,8) = IGRAPH_NAN;
                    mat(i,9) = IGRAPH_NAN;
                    mat(i,10) = IGRAPH_NAN;
                    mat(i,12) = IGRAPH_NAN;
                    mat(i,19) = IGRAPH_NAN;
                    break;
                }
                /* otherwise fall through to default */
            case 6:
                if (! dir) {
                    mat(i,0) = IGRAPH_NAN;
                    mat(i,1) = IGRAPH_NAN;
                    mat(i,2) = IGRAPH_NAN;
                    mat(i,3) = IGRAPH_NAN;
                    mat(i,4) = IGRAPH_NAN;
                    mat(i,5) = IGRAPH_NAN;
                    mat(i,6) = IGRAPH_NAN;
                    mat(i,7) = IGRAPH_NAN;
                    mat(i,8) = IGRAPH_NAN;
                    mat(i,9) = IGRAPH_NAN;
                    mat(i,10) = IGRAPH_NAN;
                    mat(i,11) = IGRAPH_NAN;
                    mat(i,12) = IGRAPH_NAN;
                    mat(i,13) = IGRAPH_NAN;
                    mat(i,14) = IGRAPH_NAN;
                    mat(i,15) = IGRAPH_NAN;
                    mat(i,16) = IGRAPH_NAN;
                    mat(i,17) = IGRAPH_NAN;
                    mat(i,18) = IGRAPH_NAN;
                    mat(i,19) = IGRAPH_NAN;
                    mat(i,20) = IGRAPH_NAN;
                    mat(i,21) = IGRAPH_NAN;
                    mat(i,22) = IGRAPH_NAN;
                    mat(i,23) = IGRAPH_NAN;
                    mat(i,24) = IGRAPH_NAN;
                    mat(i,25) = IGRAPH_NAN;
                    mat(i,26) = IGRAPH_NAN;
                    mat(i,27) = IGRAPH_NAN;
                    mat(i,28) = IGRAPH_NAN;
                    mat(i,29) = IGRAPH_NAN;
                    mat(i,30) = IGRAPH_NAN;
                    mat(i,31) = IGRAPH_NAN;
                    mat(i,32) = IGRAPH_NAN;
                    mat(i,33) = IGRAPH_NAN;
                    mat(i,35) = IGRAPH_NAN;
                    mat(i,38) = IGRAPH_NAN;
                    mat(i,44) = IGRAPH_NAN;
                    mat(i,50) = IGRAPH_NAN;
                    mat(i,51) = IGRAPH_NAN;
                    mat(i,54) = IGRAPH_NAN;
                    mat(i,74) = IGRAPH_NAN;
                    mat(i,77) = IGRAPH_NAN;
                    mat(i,89) = IGRAPH_NAN;
                    mat(i,120) = IGRAPH_NAN;
                    break;
                }
                /* otherwise fall through to default */
            default:
                massert("motifsParticipation: this line should never be reached");
            }

        guard.deactivate();
        return mat;
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
        igCheck(igraph_motifs_randesu_estimate(&graph, &res, size, &ig_cut_prob, sample_size, nullptr));
        return res;
    }

    mint motifsEstimateVerts(mint size, mma::RealTensorRef cut_prob, mma::RealTensorRef vs) const {
        igraph_integer_t res;
        igraph_vector_t ig_cut_prob = igVectorView(cut_prob);
        igraph_vector_t ig_vs = igVectorView(vs);
        igCheck(igraph_motifs_randesu_estimate(&graph, &res, size, &ig_cut_prob, 0, &ig_vs));
        return res;
    }

    mma::IntTensorRef triangles() const {
        igIntVector vec;
        igCheck(igraph_list_triangles(&graph, &vec.vec));
        return vec.makeMTensor();
    }

    mma::RealTensorRef countAdjacentTriangles(mma::RealTensorRef vs) const {
        igraph_vector_t vsvec = igVectorView(vs);
        igVector res;
        igCheck(igraph_adjacent_triangles(&graph, &res.vec, vs.length() == 0 ? igraph_vss_all() : igraph_vss_vector(&vsvec)));
        return res.makeMTensor();
    }

    // Shortest paths

    mma::RealMatrixRef shortestPaths(mma::RealTensorRef from, mma::RealTensorRef to) const {
        igMatrix res;
        igraph_vector_t fromv = igVectorView(from);
        igraph_vector_t tov   = igVectorView(to);
        igCheck(igraph_shortest_paths(
                    &graph, &res.mat,
                    from.length() == 0 ? igraph_vss_all() : igraph_vss_vector(&fromv),
                    to.length() == 0 ? igraph_vss_all() : igraph_vss_vector(&tov),
                    IGRAPH_OUT));
        return res.makeMTensor();
    }

    mma::RealTensorRef shortestPathsDijkstra(mma::RealTensorRef from, mma::RealTensorRef to) const {
        igMatrix res;
        igraph_vector_t fromv = igVectorView(from);
        igraph_vector_t tov   = igVectorView(to);
        igCheck(igraph_shortest_paths_dijkstra(
                    &graph, &res.mat,
                    from.length() == 0 ? igraph_vss_all() : igraph_vss_vector(&fromv),
                    to.length() == 0 ? igraph_vss_all() : igraph_vss_vector(&tov),
                    passWeights(), IGRAPH_OUT));
        return res.makeMTensor();
    }

    mma::RealTensorRef shortestPathsBellmanFord(mma::RealTensorRef from, mma::RealTensorRef to) const {
        igMatrix res;
        igraph_vector_t fromv = igVectorView(from);
        igraph_vector_t tov   = igVectorView(to);
        igCheck(igraph_shortest_paths_bellman_ford(
                    &graph, &res.mat,
                    from.length() == 0 ? igraph_vss_all() : igraph_vss_vector(&fromv),
                    to.length() == 0 ? igraph_vss_all() : igraph_vss_vector(&tov),
                    passWeights(), IGRAPH_OUT));
        return res.makeMTensor();
    }

    mma::RealTensorRef shortestPathsJohnson(mma::RealTensorRef from, mma::RealTensorRef to) const {
        igMatrix res;
        igraph_vector_t fromv = igVectorView(from);
        igraph_vector_t tov   = igVectorView(to);
        igCheck(igraph_shortest_paths_johnson(
                    &graph, &res.mat,
                    from.length() == 0 ? igraph_vss_all() : igraph_vss_vector(&fromv),
                    to.length() == 0 ? igraph_vss_all() : igraph_vss_vector(&tov),
                    passWeights()));
        return res.makeMTensor();
    }

    // Shortest path tree

    mma::IntTensorRef shortestPathTreeEdges(mint from) const {
        igLongVector pred_e;

        igCheck(igraph_get_shortest_paths(
                    &graph, nullptr, nullptr,
                    from, igraph_vss_all(),
                    IGRAPH_OUT,
                    nullptr, &pred_e.vec));

        return pred_e.makeMTensor();
    }

    mma::IntTensorRef shortestPathTreeEdgesDijkstra(mint from) const {
        igLongVector pred_e;

        igCheck(igraph_get_shortest_paths_dijkstra(
                    &graph, nullptr, nullptr,
                    from, igraph_vss_all(),
                    passWeights(), IGRAPH_OUT,
                    nullptr, &pred_e.vec));

        return pred_e.makeMTensor();
    }

    mma::IntTensorRef shortestPathTreeEdgesBellmanFord(mint from) const {
        igLongVector pred_e;

        igCheck(igraph_get_shortest_paths_bellman_ford(
                    &graph, nullptr, nullptr,
                    from, igraph_vss_all(),
                    passWeights(), IGRAPH_OUT,
                    nullptr, &pred_e.vec));

        return pred_e.makeMTensor();
    }

    // note that this is a constructor!
    void shortestPathTree(const IG &ig, mint from) {
        igLongVector pred_v;

        igCheck(igraph_get_shortest_paths_dijkstra(
                    &ig.graph, nullptr, nullptr,
                    from, igraph_vss_all(),
                    passWeights(), IGRAPH_OUT,
                    &pred_v.vec, nullptr));

        igVector edges;
        edges.reserve(2*ig.vertexCount());

        for (mint v=0; v < pred_v.size(); ++v) {
            mint u = pred_v[v];
            if (u >= 0) {
                edges.push_back(u);
                edges.push_back(v);
            }
        }

        fromEdgeList(edges, ig.vertexCount(), true);
    }

    // note that this is a constructor!
    void shortestPathTreeDijkstra(const IG &ig, mint from) {
        igLongVector pred_v, pred_e;

        igCheck(igraph_get_shortest_paths_dijkstra(
                    &ig.graph, nullptr, nullptr,
                    from, igraph_vss_all(),
                    passWeights(), IGRAPH_OUT,
                    &pred_v.vec, &pred_e.vec));

        igVector edges;
        edges.reserve(2*ig.vertexCount());

        weighted = true;
        weights.reserve(2*ig.vertexCount());

        for (mint v=0; v < pred_v.size(); ++v) {
            mint u = pred_v[v];
            if (u >= 0) {
                edges.push_back(u);
                edges.push_back(v);
                weights.push_back(ig.weights[pred_e[v]]);
            }
        }

        fromEdgeList(edges, ig.vertexCount(), true);
    }

    // note that this is a constructor!
    void shortestPathTreeBellmanFord(const IG &ig, mint from) {
        igLongVector pred_v, pred_e;

        igCheck(igraph_get_shortest_paths_dijkstra(
                    &ig.graph, nullptr, nullptr,
                    from, igraph_vss_all(),
                    passWeights(), IGRAPH_OUT,
                    &pred_v.vec, &pred_e.vec));

        igVector edges;
        edges.reserve(2*ig.vertexCount());

        weighted = true;
        weights.reserve(2*ig.vertexCount());

        for (mint v=0; v < pred_v.size(); ++v) {
            mint u = pred_v[v];
            if (u >= 0) {
                edges.push_back(u);
                edges.push_back(v);
                weights.push_back(ig.weights[pred_e[v]]);
            }
        }

        fromEdgeList(edges, ig.vertexCount(), true);
    }

    // Neighbourhoods

    mma::RealTensorRef neighborhoodSize(mma::RealTensorRef vs, mint mindist, mint maxdist, mint mode) const {
        igVector res;
        igraph_vector_t vsvec = igVectorView(vs);
        igCheck(igraph_neighborhood_size(&graph, &res.vec, vs.length() == 0 ? igraph_vss_all() : igraph_vss_vector(&vsvec), maxdist, igNeighborMode(mode, "neighborhoodSize"), mindist));
        return res.makeMTensor();
    }

    mma::RealTensorRef averageNeighborDegree(mma::RealTensorRef vs, mint neiMode, mint degMode) const {
        igVector res;
        igraph_vector_t vsvec = igVectorView(vs);
        const char *context = "average neighbor degree";

        igCheck(igraph_avg_nearest_neighbor_degree(&graph, vs.length() == 0 ? igraph_vss_all() : igraph_vss_vector(&vsvec), igNeighborMode(neiMode, context), igNeighborMode(degMode, context), &res.vec, nullptr, passWeights()));

        return res.makeMTensor();
    }

    mma::RealTensorRef averageDegreeConnectivity(mint neiMode, mint degMode) const {
        igVector res;
        const char *context = "average neighbor connectivity";

        igCheck(igraph_avg_nearest_neighbor_degree(&graph, igraph_vss_all(), igNeighborMode(neiMode, context), igNeighborMode(degMode, context), nullptr, &res.vec, passWeights()));

        return res.makeMTensor();
    }

    // Shortest path length histograms

    mma::RealTensorRef shortestPathCounts() const {
        igVector res;
        double unconnected;
        igCheck(igraph_path_length_hist(&graph, &res.vec, &unconnected, true));
        return res.makeMTensor();
    }

    mma::RealTensorRef shortestPathCounts2(mma::IntTensorRef vs) const {
        if (vs.length() == 0)
            return shortestPathCounts();

        std::vector<double> counts;
        igMatrix mat;
        igraph_integer_t vc = igraph_vcount(&graph);

        bool directed = directedQ();

        for (auto vertex : vs) {
            mma::check_abort();

            igCheck(igraph_shortest_paths(&graph, &mat.mat, igraph_vss_1(igraph_integer_t(vertex)), igraph_vss_all(), IGRAPH_OUT));

            for (igraph_integer_t i=0; i < vc; ++i) {
                igraph_real_t l = VECTOR(mat.mat.data)[i];
                if (l != 0 && l != IGRAPH_INFINITY) {
                    size_t k = size_t(l);
                    if (k > counts.size())
                        counts.resize(k, 0.0);
                    counts[k-1] += 1;
                }
            }
            if (! directed) {
                for (auto i : vs) {
                    igraph_real_t l = VECTOR(mat.mat.data)[i];
                    if (l != 0 && l != IGRAPH_INFINITY)
                        counts[size_t(l)-1] -= 0.5;
                }
            }
        }

        return mma::makeVector<double>(counts.size(), counts.data());
    }

    mma::IntTensorRef shortestPathWeightedHistogram(double binsize, mma::RealTensorRef from, mma::RealTensorRef to, mint method) const {
        if (weighted && igraph_vector_min(&weights.vec) < 0)
            throw mma::LibraryError("shortestPathWeightedHistogram: Negative edge weights are not supported.");

        std::vector<mint> hist;
        igMatrix mat;

        const bool toall = to.length() == 0;
        igraph_vector_t tov = igVectorView(to);
        const mint tolen = toall ? vertexCount() : to.length();

        for (mint i=0; i < from.length(); ++i) {
            mma::check_abort();

            switch (method) {
            case 0:
                igCheck(igraph_shortest_paths_dijkstra(&graph, &mat.mat, igraph_vss_1(igraph_integer_t(from[i])), toall ? igraph_vss_all() : igraph_vss_vector(&tov), passWeights(), IGRAPH_OUT));
                break;
            case 1:
                igCheck(igraph_shortest_paths_bellman_ford(&graph, &mat.mat, igraph_vss_1(igraph_integer_t(from[i])), toall ? igraph_vss_all() : igraph_vss_vector(&tov), passWeights(), IGRAPH_OUT));
                break;
            default:
                throw mma::LibraryError("shortestPathWeightedHistogram: Unknown method.");
            }

            for (igraph_integer_t j=0; j < tolen; ++j) {
                if (toall) {
                    if (from[i] == j)
                        continue;
                } else {
                    if (from[i] == to[j])
                        continue;
                }
                double length = VECTOR(mat.mat.data)[j];
                if (igraph_is_inf(length))
                    continue;
                mint idx = std::floor(length / binsize);
                if (idx >= hist.size()) {
                    hist.reserve(std::ceil(1.5*(idx + 1)));
                    hist.resize(idx+1, 0);
                }
                hist[idx] += 1;
            }
        }
        mma::IntTensorRef res = mma::makeVector<mint>(hist.size(), hist.data());
        return res;
    }

    // Diameter

    double diameter(bool components) const {
        double diam;
        igCheck(igraph_diameter(&graph, &diam, nullptr, nullptr, nullptr, nullptr, true, components));
        return diam;
    }

    mma::RealTensorRef findDiameter(bool components) const {
        igVector vertex_path;
        igCheck(igraph_diameter(&graph, nullptr, nullptr, nullptr, &vertex_path.vec, nullptr, true, components));
        return vertex_path.makeMTensor();
    }

    /* TODO */
    mma::RealTensorRef findDiameterEdges(bool components) const {
        igVector edge_path;
        igCheck(igraph_diameter(&graph, nullptr, nullptr, nullptr, nullptr, &edge_path.vec, true, components));
        return edge_path.makeMTensor();
    }

    double diameterDijkstra(bool components) const {
        double diam;
        igCheck(igraph_diameter_dijkstra(&graph, passWeights(), &diam, nullptr, nullptr, nullptr, nullptr, true, components));
        return diam;
    }

    mma::RealTensorRef findDiameterDijkstra(bool components) const {
        igVector vertex_path;
        igCheck(igraph_diameter_dijkstra(&graph, passWeights(), nullptr, nullptr, nullptr, &vertex_path.vec, nullptr, true, components));
        return vertex_path.makeMTensor();
    }

    /* TODO */
    mma::RealTensorRef findDiameterDijkstraEdges(bool components) const {
        igVector edge_path;
        igCheck(igraph_diameter_dijkstra(&graph, passWeights(), nullptr, nullptr, nullptr, nullptr, &edge_path.vec, true, components));
        return edge_path.makeMTensor();
    }

    double pseudoDiameter(mint start, bool components) const {
        double diam;

        igraph_integer_t f, t;
        igCheck(igraph_pseudo_diameter(&graph, &diam, start, &f, &t, true, components));

        return diam;
    }

    mma::RealTensorRef findPseudoDiameter(mint start, bool components) const {
        igraph_integer_t from, to;
        igVector vertex_path;

        igCheck(igraph_pseudo_diameter(&graph, nullptr, start, &from, &to, true, components));

        if (from < 0)
            return vertex_path.makeMTensor();

        igCheck(igraph_get_shortest_path(&graph, &vertex_path.vec, nullptr, from, to, IGRAPH_OUT));

        return vertex_path.makeMTensor();
    }

    double pseudoDiameterDijkstra(mint start, bool components) const {
        double diam;
        igCheck(igraph_pseudo_diameter_dijkstra(&graph, passWeights(), &diam, start, nullptr, nullptr, true, components));
        return diam;
    }

    mma::RealTensorRef findPseudoDiameterDijkstra(mint start, bool components) const {
        igraph_integer_t from, to;
        igVector vertex_path;

        igCheck(igraph_pseudo_diameter_dijkstra(&graph, passWeights(), nullptr, start, &from, &to, true, components));

        if (from < 0)
            return vertex_path.makeMTensor();

        igCheck(igraph_get_shortest_path_dijkstra(&graph, &vertex_path.vec, nullptr, from, to, passWeights(), IGRAPH_OUT));

        return vertex_path.makeMTensor();
    }

    double averagePathLength(bool components) const {
        double res;
        igCheck(igraph_average_path_length(&graph, &res, nullptr, /* directed */ true, components));
        return res;
    }

    double averagePathLengthDijkstra(bool components) const {
        double res;
        igCheck(igraph_average_path_length_dijkstra(&graph, &res, nullptr, passWeights(), /* directed */ true, components));
        return res;
    }

    double averagePathLengthWeighted(mint method, bool components) const {
        /* method values:
         * 0 - unweighted
         * 1 - Dijkstra
         * 2 - Bellman-Ford
         * 3 - Johnson
         */
        if (method == 0 || ! weightedQ())
            return averagePathLength(components);
        if (method == 1)
            return averagePathLengthDijkstra(components);
        if (method < 0 || method > 3)
            throw mma::LibraryError("averagePathLengthWeighted: Unknown method.");

        igMatrix mat;
        int n = vertexCount();

        long cnt = 0;
        double sum = 0.0;
        for (int i=0; i < n; ++i) {
            mma::check_abort();

            switch (method) {
            /*
            case 1:
                igCheck(igraph_shortest_paths_dijkstra(&graph, &mat.mat, igraph_vss_1(i), igraph_vss_all(), passWeights(), IGRAPH_OUT));
                break;
            */
            case 2:
                igCheck(igraph_shortest_paths_bellman_ford(&graph, &mat.mat, igraph_vss_1(i), igraph_vss_all(), passWeights(), IGRAPH_OUT));
                break;
            case 3:
                igCheck(igraph_shortest_paths_johnson(&graph, &mat.mat, igraph_vss_1(i), igraph_vss_all(), passWeights()));
                break;
            }

            for (int j=0; j < mat.ncol(); ++j)
                if (i != j) {
                    if (VECTOR(mat.mat.data)[j] != IGRAPH_INFINITY) {
                        cnt += 1;
                        sum += VECTOR(mat.mat.data)[j];
                    } else {
                        if (! components)
                            return IGRAPH_INFINITY;
                    }
                }
        }

        return sum/cnt;
    }

    double globalEfficiency() const {
        igraph_real_t res;
        igCheck(igraph_global_efficiency(&graph, &res, passWeights(), true));
        return res;
    }

    mma::RealTensorRef localEfficiency(mma::RealTensorRef vs, bool directed, mint mode) const {
        igVector vec;
        igraph_vector_t vsvec = igVectorView(vs);
        igCheck(igraph_local_efficiency(&graph, &vec.vec, vs.length() == 0 ? igraph_vss_all() : igraph_vss_vector(&vsvec), passWeights(), directed, igNeighborMode(mode, "local efficiency")));
        return vec.makeMTensor();
    }

    double averageLocalEfficiency(bool directed, mint mode) const {
        igraph_real_t res;
        igCheck(igraph_average_local_efficiency(&graph, &res, passWeights(), directed, igNeighborMode(mode, "average local efficiency")));
        return res;
    }

    mint girth() const {
        igraph_integer_t res;
        igCheck(igraph_girth(&graph, &res, nullptr));
        return res;
    }

    mma::RealTensorRef eccentricity(mma::RealTensorRef vs) const {
        igraph_vector_t vsvec = igVectorView(vs);
        igVector res;
        igCheck(igraph_eccentricity(&graph, &res.vec, vs.length() == 0 ? igraph_vss_all() : igraph_vss_vector(&vsvec), IGRAPH_OUT));
        return res.makeMTensor();
    }

    double radius() const {
        double res;
        igCheck(igraph_radius(&graph, &res, IGRAPH_OUT));
        return res;
    }

    // Cliques

    mma::IntTensorRef cliques(mint min, mint max) const {
        igList list;
        igCheck(igraph_cliques(&graph, &list.list, min, max));
        return packListIntoIntTensor(list);
    }

    mma::RealTensorRef cliqueDistribution(mint min, mint max) const {
        igVector hist;
        igCheck(igraph_clique_size_hist(&graph, &hist.vec, min, max));
        return hist.makeMTensor();
    }

    mma::IntTensorRef maximalCliques(mint min, mint max) const {
        igList list;
        igCheck(igraph_maximal_cliques(&graph, &list.list, min, max));
        return packListIntoIntTensor(list);
    }

    mint maximalCliquesCount(mint min, mint max) const {
        igraph_integer_t res;
        igCheck(igraph_maximal_cliques_count(&graph, &res, min, max));
        return res;
    }

    // mma::IntTensorRef
    mma::RealTensorRef
    maximalCliqueDistribution(mint min, mint max) const {
        /*
        struct clique_data {
            std::vector<mint> hist;
        } cd;

        cd.hist.reserve(50);

        // TODO interruptability
        struct {
            static igraph_bool_t handle(igraph_vector_t *clique, void *data) {
                clique_data *cd = static_cast<clique_data *>(data);
                long clique_size = igraph_vector_size(clique);
                if (cd->hist.size() < clique_size) {
                    if (cd->hist.capacity() < clique_size)
                        cd->hist.reserve(2*clique_size);
                    cd->hist.resize(clique_size, 0);
                }
                cd->hist[clique_size-1] += 1;
                igraph_vector_destroy(clique);
                igraph_free(clique);
                return true;
            }
        } clique_counter;

        igCheck(igraph_maximal_cliques_callback(&graph, &clique_counter.handle, &cd, min, max));

        return mma::makeVector<mint>(cd.hist.size(), cd.hist.data());
        */
        igVector hist;
        igCheck(igraph_maximal_cliques_hist(&graph, &hist.vec, min, max));
        return hist.makeMTensor();
    }

    mma::IntTensorRef largestCliques() const {
        igList list;
        igCheck(igraph_largest_cliques(&graph, &list.list));
        return packListIntoIntTensor(list);
    }

    mint cliqueNumber() const {
        igraph_integer_t res;
        igCheck(igraph_clique_number(&graph, &res));
        return res;
    }

    // Weighted cliques

    mma::IntTensorRef cliquesWeighted(mint wmin, mint wmax, mma::RealTensorRef vertex_weights, bool maximal) const {
        igraph_vector_t weights = igVectorView(vertex_weights);
        igList list;
        igCheck(igraph_weighted_cliques(&graph, &weights, &list.list, wmin, wmax, maximal));
        return packListIntoIntTensor(list);
    }

    mma::IntTensorRef largestCliquesWeighted(mma::RealTensorRef vertex_weights) const {
        igraph_vector_t weights = igVectorView(vertex_weights);
        igList list;
        igCheck(igraph_largest_weighted_cliques(&graph, &weights, &list.list));
        return packListIntoIntTensor(list);
    }

    mint cliqueNumberWeighted(mma::RealTensorRef vertex_weights) const {
        igraph_vector_t weights = igVectorView(vertex_weights);
        double res;
        igCheck(igraph_weighted_clique_number(&graph, &weights, &res));
        return mint(res);
    }

    // Independent vertex sets

    mma::IntTensorRef independentVertexSets(mint min, mint max) const {
        igList list;
        igCheck(igraph_independent_vertex_sets(&graph, &list.list, min, max));
        return packListIntoIntTensor(list);
    }

    mma::IntTensorRef largestIndependentVertexSets() const {
        igList list;
        igCheck(igraph_largest_independent_vertex_sets(&graph, &list.list));
        return packListIntoIntTensor(list);
    }

    mma::IntTensorRef maximalIndependentVertexSets() const {
        igList list;
        igCheck(igraph_maximal_independent_vertex_sets(&graph, &list.list));
        return packListIntoIntTensor(list);
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
                    nullptr, nullptr, nullptr, nullptr));
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
                    nullptr, nullptr, nullptr, nullptr, nullptr, nullptr));
        return mat.makeMTensor();
    }

    mma::RealTensorRef layoutFruchtermanReingold(
            mma::RealMatrixRef initial, bool use_seed,
            mint niter, double start_temp, mint grid_method,
            bool constrain,
            mma::RealTensorRef minx, mma::RealTensorRef miny,
            mma::RealTensorRef maxx, mma::RealTensorRef maxy) const
    {
        igMatrix mat;
        mat.copyFromMTensor(initial);
        igraph_layout_grid_t grid;
        switch (grid_method) {
        case 0: grid = IGRAPH_LAYOUT_GRID; break;
        case 1: grid = IGRAPH_LAYOUT_NOGRID; break;
        case 2: grid = IGRAPH_LAYOUT_AUTOGRID; break;
        default: throw mma::LibraryError("layoutFruchtermanReingold: Unknown method option.");
        }

        igraph_vector_t iminx, imaxx, iminy, imaxy;
        if (constrain) {
            iminx = igVectorView(minx);
            imaxx = igVectorView(maxx);
            iminy = igVectorView(miny);
            imaxy = igVectorView(maxy);
        }

        igCheck(igraph_layout_fruchterman_reingold(
                    &graph, &mat.mat, use_seed, niter, start_temp, grid, passWeights(),
                    constrain ? &iminx : nullptr, constrain ? &imaxx : nullptr,
                    constrain ? &iminy : nullptr, constrain ? &imaxy : nullptr));
        return mat.makeMTensor();
    }

    mma::RealTensorRef layoutFruchtermanReingold3D(
            mma::RealMatrixRef initial, bool use_seed,
            mint niter, double start_temp,
            bool constrain,
            mma::RealTensorRef minx, mma::RealTensorRef miny, mma::RealTensorRef minz,
            mma::RealTensorRef maxx, mma::RealTensorRef maxy, mma::RealTensorRef maxz) const
    {
        igMatrix mat;
        mat.copyFromMTensor(initial);

        igraph_vector_t iminx, imaxx, iminy, imaxy, iminz, imaxz;
        if (constrain) {
            iminx = igVectorView(minx);
            imaxx = igVectorView(maxx);
            iminy = igVectorView(miny);
            imaxy = igVectorView(maxy);
            iminz = igVectorView(minz);
            imaxz = igVectorView(maxz);
        }

        igCheck(igraph_layout_fruchterman_reingold_3d(
                    &graph, &mat.mat, use_seed, niter, start_temp, passWeights(),
                    constrain ? &iminx : nullptr, constrain ? &imaxx : nullptr,
                    constrain ? &iminy : nullptr, constrain ? &imaxy : nullptr,
                    constrain ? &iminz : nullptr, constrain ? &imaxz : nullptr));
        return mat.makeMTensor();
    }

    mma::RealTensorRef layoutGEM(
            mma::RealTensorRef initial, bool use_seed,
            mint maxiter, double temp_max, double temp_min, double temp_init) const
    {
        // workaround for 0- and 1-vertex graph
        if (vertexCount() == 0)
            return mma::makeMatrix<double>(0, 2);
        if (vertexCount() == 1) {
            if (use_seed)
                return initial.clone();
            else {
                auto res = mma::makeMatrix<double>(1, 2);
                res(0,0) = res(0,1) = 0.0;
                return res;
            }
        }

        igMatrix mat;
        mat.copyFromMTensor(initial);
        igCheck(igraph_layout_gem(&graph, &mat.mat, use_seed, maxiter, temp_max, temp_min, temp_init));
        return mat.makeMTensor();
    }

    mma::RealTensorRef layoutDavidsonHarel(
            mma::RealTensorRef initial, bool use_seed,
            mint maxiter, mint fineiter, double cool_fact,
            double weight_node_dist, double weight_border,
            double weight_edge_lengths, double weight_edge_crossings,
            double weight_node_edge_dist) const
    {
        igMatrix mat;
        mat.copyFromMTensor(initial);
        igCheck(igraph_layout_davidson_harel(
                    &graph, &mat.mat, use_seed,
                    maxiter, fineiter, cool_fact,
                    weight_node_dist, weight_border, weight_edge_lengths, weight_edge_crossings,weight_node_edge_dist));
        return mat.makeMTensor();
    }

    /*
    mma::RealTensorRef layoutMDS(mma::RealMatrixRef dist, mint dim) const {
        igMatrix dmat, mat;
        dmat.copyFromMTensor(dist);
        igCheck(igraph_layout_mds(&graph, &mat.mat, dist.cols() == 0 ? nullptr : &dmat.mat, dim, nullptr));
        return mat.makeMTensor();
    }
     */

    mma::RealTensorRef layoutReingoldTilford(mma::RealTensorRef roots, bool directed) const {
        igMatrix mat;
        igraph_vector_t rootvec = igVectorView(roots);
        igCheck(igraph_layout_reingold_tilford(
                &graph, &mat.mat, directed ? IGRAPH_OUT : IGRAPH_ALL,
                roots.length() == 0 ? nullptr : &rootvec, nullptr));
        return mat.makeMTensor();
    }

    mma::RealTensorRef layoutReingoldTilfordCircular(mma::RealTensorRef roots, bool directed) const {
        igMatrix mat;
        igraph_vector_t rootvec = igVectorView(roots);
        igCheck(igraph_layout_reingold_tilford_circular(
                &graph, &mat.mat, directed ? IGRAPH_OUT : IGRAPH_ALL,
                roots.length() == 0 ? nullptr : &rootvec, nullptr));
        return mat.makeMTensor();
    }

    mma::RealTensorRef layoutDrL(mma::RealTensorRef initial, bool use_seed, mint opt) const {
        igMatrix mat;
        mat.copyFromMTensor(initial);

        igraph_layout_drl_default_t templ;
        switch (opt) {
        case 1: templ = IGRAPH_LAYOUT_DRL_DEFAULT; break;
        case 2: templ = IGRAPH_LAYOUT_DRL_COARSEN; break;
        case 3: templ = IGRAPH_LAYOUT_DRL_COARSEST; break;
        case 4: templ = IGRAPH_LAYOUT_DRL_REFINE; break;
        case 5: templ = IGRAPH_LAYOUT_DRL_FINAL; break;
        default: throw mma::LibraryError("Invalid settings template for DrL layout.");
        }

        igraph_layout_drl_options_t options;
        igCheck(igraph_layout_drl_options_init(&options, templ));

        igCheck(igraph_layout_drl(&graph, &mat.mat, use_seed, &options, passWeights(), nullptr));

        return mat.makeMTensor();
    }

    mma::RealTensorRef layoutDrL3D(mma::RealTensorRef initial, bool use_seed, mint opt) const {
        // workaround for 0- and 1-vertex graph
        if (vertexCount() == 0)
            return mma::makeMatrix<double>(0, 3);
        if (vertexCount() == 1) {
            if (use_seed)
                return initial.clone();
            else {
                auto res = mma::makeMatrix<double>(1, 3);
                res(0,0) = res(0,1) = res(0,2) = 0.0;
                return res;
            }
        }

        igMatrix mat;
        mat.copyFromMTensor(initial);

        igraph_layout_drl_default_t templ;
        switch (opt) {
        case 1: templ = IGRAPH_LAYOUT_DRL_DEFAULT; break;
        case 2: templ = IGRAPH_LAYOUT_DRL_COARSEN; break;
        case 3: templ = IGRAPH_LAYOUT_DRL_COARSEST; break;
        case 4: templ = IGRAPH_LAYOUT_DRL_REFINE; break;
        case 5: templ = IGRAPH_LAYOUT_DRL_FINAL; break;
        default: throw mma::LibraryError("Invalid settings template for 3D DrL layout.");
        }

        igraph_layout_drl_options_t options;
        igCheck(igraph_layout_drl_options_init(&options, templ));

        igCheck(igraph_layout_drl_3d(&graph, &mat.mat, use_seed, &options, passWeights(), nullptr));

        return mat.makeMTensor();
    }

    mma::RealTensorRef layoutBipartite(mma::IntTensorRef typevec, double hgap, double vgap, mint maxiter) const {
        igMatrix mat;

        igBoolVector types(typevec.length());
        std::copy(typevec.begin(), typevec.end(), types.begin());

        igCheck(igraph_layout_bipartite(&graph, &types.vec, &mat.mat, hgap, vgap, maxiter));
        return mat.makeMTensor();
    }

    // Clusterig coefficient

    double transitivityUndirected(bool exclude) const {
        if (directedQ())
            mma::message("Edge directions are ignored for clustering coefficient calculations.", mma::M_WARNING);

        double res;
        igCheck(igraph_transitivity_undirected(&graph, &res, exclude ? IGRAPH_TRANSITIVITY_NAN : IGRAPH_TRANSITIVITY_ZERO));
        return res;
    }

    mma::RealTensorRef transitivityLocalUndirected(bool exclude) const {
        if (directedQ())
            mma::message("Edge directions are ignored for clustering coefficient calculations.", mma::M_WARNING);

        igVector vec;
        igCheck(igraph_transitivity_local_undirected(&graph, &vec.vec, igraph_vss_all(), exclude ? IGRAPH_TRANSITIVITY_NAN : IGRAPH_TRANSITIVITY_ZERO));
        return vec.makeMTensor();
    }

    double transitivityAverageLocalUndirected(bool exclude) const {
        if (directedQ())
            mma::message("Edge directions are ignored for clustering coefficient calculations.", mma::M_WARNING);

        double res;
        igCheck(igraph_transitivity_avglocal_undirected(&graph, &res, exclude ? IGRAPH_TRANSITIVITY_NAN : IGRAPH_TRANSITIVITY_ZERO));
        return res;
    }

    mma::RealTensorRef transitivityBarrat(bool exclude) const {
        if (directedQ())
            mma::message("Edge directions are ignored for clustering coefficient calculations.", mma::M_WARNING);

        igVector vec;
        igCheck(igraph_transitivity_barrat(&graph, &vec.vec, igraph_vss_all(), passWeights(), exclude ? IGRAPH_TRANSITIVITY_NAN : IGRAPH_TRANSITIVITY_ZERO));
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
        igCheck(igraph_maximum_cardinality_search(&graph, &vec.vec, nullptr));
        return vec.makeMTensor();
    }

    bool chordalQ() const {
        igraph_bool_t res;
        igCheck(igraph_is_chordal(&graph, nullptr, nullptr, &res, nullptr, nullptr));
        return res;
    }

    mma::RealTensorRef chordalCompletion() const {
        igraph_bool_t chordal;
        igVector fillin;
        igCheck(igraph_is_chordal(&graph, nullptr, nullptr, &chordal, &fillin.vec, nullptr));
        return fillin.makeMTensor();
    }

    // Vertex separators

    mma::IntTensorRef minimumSizeSeparators() const {
        igList list;
        igCheck(igraph_minimum_size_separators(&graph, &list.list));
        return packListIntoIntTensor(list);
    }

    mma::IntTensorRef minimalSeparators() const {
        igList list;
        igCheck(igraph_all_minimal_st_separators(&graph, &list.list));
        return packListIntoIntTensor(list);
    }

    bool separatorQ(mma::RealTensorRef vs) const {
        igraph_bool_t res;
        igraph_vector_t vsvec = igVectorView(vs);

        igCheck(igraph_is_separator(&graph, igraph_vss_vector(&vsvec), &res));
        return res;
    }

    bool minSeparatorQ(mma::RealTensorRef vs) const {
        igraph_bool_t res;
        igraph_vector_t vsvec = igVectorView(vs);

        igCheck(igraph_is_minimal_separator(&graph, igraph_vss_vector(&vsvec), &res));
        return res;
    }

    // Maximum flow

    // note that this is a constructor!
    mma::RealTensorRef gomoryHuTree(const IG &ig) {
        igVector flows;
        destroy();
        igConstructorCheck(igraph_gomory_hu_tree(&ig.graph, &graph, &flows.vec, ig.passWeights()));
        return flows.makeMTensor();
    }

    // note that this is a constructor!
    void dominatorTree(const IG &ig, mint root) {
        destroy();
        igConstructorCheck(igraph_dominator_tree(&ig.graph, root, nullptr, &graph, nullptr, IGRAPH_OUT));
    }

    mma::RealTensorRef immediateDominators(mint root) const {
        igVector dom;
        igCheck(igraph_dominator_tree(&graph, root, &dom.vec, nullptr, nullptr, IGRAPH_OUT));
        for (auto &el : dom) {
            if (! (el >= 0))
                el = -1;
        }
        return dom.makeMTensor();
    }

    /*
    mma::SparseArrayRef<double> maxFlowMatrix(mint s, mint t) const {
        igVector flow;
        igCheck(igraph_maxflow(&graph, nullptr, &flow.vec, nullptr, nullptr, nullptr, s, t, passWeights(), nullptr));

        mint ec = edgeCount();

        massert(ec == flow.length());

        auto dims = mma::makeVector<mint>(2);
        dims[0] = dims[1] = vertexCount();

        if (directedQ()) {
            auto pos = mma::makeMatrix<mint>(ec, 2);
            auto val = mma::makeVector<double>(ec);

            for (mint i=0; i < ec; ++i) {
                long s = IGRAPH_FROM(&graph, i), t = IGRAPH_TO(&graph, i);
                pos(i, 0) = s+1;
                pos(i, 1) = t+1;
                val[i] = flow[i];
            }

            return mma::makeSparseArray<double>(pos, val, dims);
        } else {
            auto pos = mma::makeMatrix<mint>(2*ec, 2);
            auto val = mma::makeVector<double>(2*ec);

            for (mint i=0; i < ec; ++i) {
                long s = IGRAPH_FROM(&graph, i), t = IGRAPH_TO(&graph, i);
                if (s > t)
                    std::swap(s, t);
                pos(2*i, 0) = s+1;
                pos(2*i, 1) = t+1;
                pos(2*i+1, 0) = t+1;
                pos(2*i+1, 1) = s+1;
                val[2*i] = flow[i];
                val[2*i+1] = -flow[i];
            }

            return mma::makeSparseArray<double>(pos, val, dims);
        }
    }
    */

    mma::RealTensorRef maxFlow(mint s, mint t, mma::RealTensorRef tcapacity) const {
        igraph_vector_t capacity = igVectorView(tcapacity);
        igVector flows;
        igCheck(igraph_maxflow(&graph, nullptr, &flows.vec, nullptr, nullptr, nullptr, s, t, tcapacity.length() == 0 ? nullptr : &capacity, nullptr));
        return flows.makeMTensor();
    }

    double maxFlowValue(mint s, mint t, mma::RealTensorRef tcapacity) const {
        igraph_vector_t capacity = igVectorView(tcapacity);
        double value;
        igCheck(igraph_maxflow_value(&graph, &value, s, t, tcapacity.length() == 0 ? nullptr : &capacity, nullptr));
        return value;
    }

    double minCutValue() const {
        double value;
        igCheck(igraph_mincut_value(&graph, &value, passWeights()));
        return value;
    }

    double minCutValueST(mint s, mint t) const {
        double value;
        igCheck(igraph_st_mincut_value(&graph, &value, s, t, passWeights()));
        return value;
    }

    mma::RealTensorRef minCut() const {
        igVector cut;
        double value;
        igCheck(igraph_mincut(&graph, &value, nullptr, nullptr, &cut.vec, passWeights()));
        return cut.makeMTensor();
    }

    mma::RealTensorRef minCutST(mint s, mint t) const {
        igVector cut;
        double value;
        igCheck(igraph_st_mincut(&graph, &value, &cut.vec, nullptr, nullptr, s, t, passWeights()));
        return cut.makeMTensor();
    }

    mma::IntTensorRef allCutsST(mint s, mint t) const {
        igList list;
        igList partitions;
        igCheck(igraph_all_st_cuts(&graph, &list.list, &partitions.list /* igraph doesn't allow a nullptr here */, s, t));
        return packListIntoIntTensor(list);
    }

    mma::IntTensorRef allMinCutsST(mint s, mint t) const {
        igList list;
        igCheck(igraph_all_st_mincuts(&graph, nullptr, &list.list, nullptr, s, t, passWeights()));
        return packListIntoIntTensor(list);
    }


    // Connected components

    mma::RealTensorRef articulationPoints() const {
        igVector vec;
        igCheck(igraph_articulation_points(&graph, &vec.vec));
        return vec.makeMTensor();
    }

    mma::IntTensorRef biconnectedComponents() const {
        igList list;
        igraph_integer_t count;
        igCheck(igraph_biconnected_components(&graph, &count, nullptr, nullptr, &list.list, nullptr));
        return packListIntoIntTensor(list);
    }

    mma::IntTensorRef biconnectedEdgeComponents() const {
        igList list;
        igraph_integer_t count;
        igCheck(igraph_biconnected_components(&graph, &count, nullptr, &list.list, nullptr, nullptr));
        return packListIntoIntTensor(list);
    }

    bool biconnectedQ() const {
        // igraph_biconnected_components does not return isolated vertices, nor does it include them
        // in the number of components. Thus to check if a graph is biconnected, it is not sufficient
        // to check if the returned number of biconnected components is 1. We also need to check for
        // isolated vertices. This is done most simply by checking if the graph is weakly connected.
        if (! connectedQ(false))
            return false;
        igraph_integer_t count;
        igCheck(igraph_biconnected_components(&graph, &count, nullptr, nullptr, nullptr, nullptr));
        return count == 1;
    }

    mma::RealTensorRef bridges() const {
        igVector res;
        igCheck(igraph_bridges(&graph, &res.vec));
        return res.makeMTensor();
    }

    mma::RealTensorRef connectedComponentSizes(bool mode) const {
        igVector csize;
        igCheck(igraph_clusters(&graph, nullptr, &csize.vec, nullptr, mode ? IGRAPH_STRONG : IGRAPH_WEAK));
        return csize.makeMTensor();
    }

    // Connectivity

    mint vertexConnectivity() const {
        igraph_integer_t res;
        igCheck(igraph_vertex_connectivity(&graph, &res, true));
        return res;
    }

    mint edgeConnectivity() const {
        igraph_integer_t res;
        igCheck(igraph_edge_connectivity(&graph, &res, true));
        return res;
    }

    mint vertexConnectivityST(mint s, mint t) const {
        igraph_integer_t res;
        igCheck(igraph_st_vertex_connectivity(&graph, &res, s, t, IGRAPH_VCONN_NEI_ERROR));
        return res;
    }

    mint edgeConnectivityST(mint s, mint t) const {
        igraph_integer_t res;
        igCheck(igraph_st_edge_connectivity(&graph, &res, s, t));
        return res;
    }

    void cohesiveBlocks(MLINK link) const {
        mlStream ml{link, "cohesiveBlocks"};
        ml >> mlCheckArgs(0);

        igList blocks;
        igVector cohesion;
        igVector parents;

        igCheck(igraph_cohesive_blocks(&graph, &blocks.list, &cohesion.vec, &parents.vec, nullptr));

        ml.newPacket();
        ml << mlHead("List", 3) << blocks << cohesion << parents;
    }

    // Graphlets

    void graphletBasis(MLINK link) const {
        mlStream ml{link, "graphletBasis"};
        ml >> mlCheckArgs(0);

        igList cliques;
        igVector thresholds;
        igCheck(igraph_graphlets_candidate_basis(&graph, passWeights(), &cliques.list, &thresholds.vec));

        ml.newPacket();
        ml << mlHead("List", 2) << cliques << thresholds;
    }

    void graphletProject(MLINK link) const  {
        mlStream ml{link, "graphletProject"};

        igList cliques;
        int niter;

        ml >> mlCheckArgs(2) >> cliques >> niter;

        igVector mu;
        igCheck(igraph_graphlets_project(&graph, passWeights(), &cliques.list, &mu.vec, false, niter));

        ml.newPacket();
        ml << mu;
    }

    void graphlets(MLINK link) const {
        mlStream ml{link, "graphlets"};
        int niter;
        ml >> mlCheckArgs(1) >> niter;

        igVector mu;
        igList cliques;
        igCheck(igraph_graphlets(&graph, passWeights(), &cliques.list, &mu.vec, niter));

        ml.newPacket();
        ml << mlHead("List", 2) << cliques << mu;
    }

    // Community detection

    double modularity(mma::RealTensorRef t, double resolution, bool directed) const {
        igraph_vector_t membership = igVectorView(t);
        double Q;
        igCheck(igraph_modularity(&graph, &membership, passWeights(), resolution, directed, &Q));
        return Q;
    }

    /*
    mma::IntTensorRef splitJoinDistance(mma::RealTensorRef c1, mma::RealTensorRef c2) const {
        igraph_integer_t d1, d2;

        igraph_vector_t comm1 = igVectorView(c1);
        igraph_vector_t comm2 = igVectorView(c2);

        igCheck(igraph_split_join_distance(&comm1, &comm2, &d1, &d2));

        mma::IntTensorRef res = mma::makeVector<mint>(2);
        res[0] = d1;
        res[1] = d2;
        return res;
    }*/

    void communityEdgeBetweenness(MLINK link) const {
        mlStream ml{link, "communityEdgeBetweenness"};

        igraph_integer_t n_communities;
        ml >> mlCheckArgs(1) >> n_communities;

        igVector result, betweenness, bridges, membership, modularity;
        igMatrix merges;

        if (n_communities == 0) { // automatic communities based on max modularity
            igCheck(igraph_community_edge_betweenness(
                        &graph, &result.vec, &betweenness.vec, &merges.mat, &bridges.vec,
                        &modularity.vec, &membership.vec,
                        true /* directed */, passWeights()
                        ));
        } else { // communities based on cluster count
            igCheck(igraph_community_edge_betweenness(
                        &graph, &result.vec, &betweenness.vec, &merges.mat, &bridges.vec,
                        nullptr /* modularity */, nullptr /* membership */,
                        true /* directed */, passWeights()
                        ));

            computeMembership(n_communities, merges, membership);
        }

        ml.newPacket();
        ml << mlHead("List", 6)
           << result << merges << betweenness << bridges << membership;

        if (n_communities == 0)
            ml << modularity;
        else
            ml << mlSymbol("None");
    }

    void communityWalktrap(MLINK link) const {
        mlStream ml{link, "communityWalktrap"};
        igraph_integer_t steps, n_communities;
        ml >> mlCheckArgs(2) >> steps >> n_communities;

        igVector modularity, membership;
        igMatrix merges;

        if (n_communities == 0) {
            igCheck(igraph_community_walktrap(&graph, passWeights(), steps, &merges.mat, &modularity.vec, &membership.vec));
        } else {
            igCheck(igraph_community_walktrap(&graph, passWeights(), steps, &merges.mat, nullptr, nullptr));
            computeMembership(n_communities, merges, membership);
            modularity.resize(1);
            igraph_modularity(&graph, &membership.vec, passWeights(), /* resolution= */ 1, /* directed= */ false, modularity.begin() /* ptr to first vec elem */);
        }

        ml.newPacket();
        ml << mlHead("List", 3)
           << merges << membership << modularity;
    }

    void communityFastGreedy(MLINK link) const {
        mlStream ml{link, "communityFastGreedy"};
        ml >> mlCheckArgs(0);

        igVector modularity, membership;
        igMatrix merges;

        igCheck(igraph_community_fastgreedy(&graph, passWeights(), &merges.mat, &modularity.vec, &membership.vec));

        ml.newPacket();
        ml << mlHead("List", 3)
           << merges << modularity << membership;
    }

    void communityMultilevel(MLINK link) const {
        mlStream ml{link, "communityMultilevel"};
        double resolution;
        ml >> mlCheckArgs(1) >> resolution;

        igVector modularity, membership;
        igMatrix memberships;

        igCheck(igraph_community_multilevel(&graph, passWeights(), resolution, &membership.vec, &memberships.mat, &modularity.vec));

        ml.newPacket();
        ml << mlHead("List", 3) << modularity << membership << memberships;
    }

    void communityLeiden(MLINK link) const {
        mlStream ml{link, "communityLeiden"};

        igraph_real_t resolution, beta;
        mint method;

        igVector vertexWeight;
        igraph_vector_t *vertex_weight_ptr = &vertexWeight.vec;

        ml >> mlCheckArgs(4) >> resolution >> beta >> method >> vertexWeight;

        switch (method) {
        case 1: /* degree/strength */
            if (weightedQ()) {
                igraph_strength(&graph, vertex_weight_ptr, igraph_vss_all(), IGRAPH_ALL, /* loops = */ true, passWeights());
                double total_edge_weight = 0.0;
                for (const auto &w : weights)
                    total_edge_weight += w;
                resolution = 0.5*resolution / total_edge_weight;
            }
            else
            {
                igraph_degree(&graph, vertex_weight_ptr, igraph_vss_all(), IGRAPH_ALL, /* loops = */ true);
                resolution = 0.5*resolution / edgeCount();
            }
            break;

        case 2: /* constant */
            vertex_weight_ptr = nullptr;
            break;

        case 3:
            if (vertexWeight.size() != vertexCount())
                throw mma::LibraryError("The length of the vertex weight vector must be the same as the number of vertices.");
            break;

        default:
            throw mma::LibraryError("communityLeiden: Invalid vertex weight computation method.");
        }

        igVector membership;
        igraph_integer_t nb_clusters;
        igraph_real_t quality;

        igCheck(igraph_community_leiden(
                    &graph, passWeights(), vertex_weight_ptr,
                    resolution, beta,
                    false, &membership.vec,
                    &nb_clusters, &quality));

        ml.newPacket();
        ml << mlHead("List", 2) << membership << quality;
    }

    void communityLabelPropagation(MLINK link) const {
        mlStream ml{link, "communityLabelPropagation"};

        igVector initial;
        igBoolVector fixed;
        ml >> mlCheckArgs(2) >> initial >> fixed;

        igVector membership;
        double modularity;

        igCheck(igraph_community_label_propagation(
                    &graph, &membership.vec, passWeights(),
                    initial.length() == 0 ? nullptr : &initial.vec,
                    fixed.length() == 0 ? nullptr : &fixed.vec,
                    &modularity));

        ml.newPacket();
        ml << mlHead("List", 2) << membership << modularity;
    }

    void communityInfoMAP(MLINK link) const {
        mlStream ml{link, "communityInfoMAP"};

        int trials;
        igVector v_weights;
        ml >> mlCheckArgs(2) >> trials >> v_weights;

        igVector membership;
        double codelen;

        igCheck(igraph_community_infomap(
                    &graph, passWeights(), v_weights.length() == 0 ? nullptr : &v_weights.vec,
                    trials, &membership.vec, &codelen));

        ml.newPacket();
        ml << mlHead("List", 2) << membership << codelen;
    }

    void communityOptimalModularity(MLINK link) const {
        mlStream ml{link, "communityOptimalModularity"};
        ml >> mlCheckArgs(0);

        igVector membership;
        double modularity;

        igCheck(igraph_community_optimal_modularity(&graph, &modularity, &membership.vec, passWeights()));

        ml.newPacket();
        ml << mlHead("List", 2) << modularity << membership;
    }

    void communitySpinGlass(MLINK link) const {
        mlStream ml{link, "communitySpinGlass"};

        igraph_integer_t spins, update_rule, method;
        igraph_bool_t parupdate;
        double starttemp, stoptemp, coolfact, gamma, gamma_minus;
        ml >> mlCheckArgs(9)
            >> spins // number of spins, i.e. max number of clusters
            >> parupdate // parallel update for IGRAPH_SPINCOMM_INP_NEG
            >> starttemp >> stoptemp >> coolfact // start, stop temperature, cooling factor
            >> update_rule
            >> gamma
            >> method
            >> gamma_minus;

        igraph_spinglass_implementation_t method_e;
        switch (method) {
        case 0: method_e = IGRAPH_SPINCOMM_IMP_ORIG; break;
        case 1: method_e = IGRAPH_SPINCOMM_IMP_NEG; break;
        default:
            throw mma::LibraryError("communitySpinGlass: Invalid method option.");
        }

        igraph_spincomm_update_t update_rule_e;
        switch (update_rule) {
        case 0: update_rule_e = IGRAPH_SPINCOMM_UPDATE_SIMPLE; break;
        case 1: update_rule_e = IGRAPH_SPINCOMM_UPDATE_CONFIG; break;
        default:
            throw mma::LibraryError("communitySpinGlass: Invalid update rule option.");
        }

        double modularity, temperature;
        igVector membership;

        igCheck(igraph_community_spinglass(
                    &graph, passWeights(),
                    &modularity, &temperature, &membership.vec, nullptr,
                    spins, parupdate, starttemp, stoptemp, coolfact,
                    update_rule_e, gamma, method_e, gamma_minus
                    ));

        ml.newPacket();
        ml << mlHead("List", 3) << membership << modularity << temperature;
    }

    void communityLeadingEigenvector(MLINK link) const {
        mlStream ml{link, "communityLeadingEigenvector"};

        igraph_integer_t steps, n_communities;
        ml >> mlCheckArgs(2) >> steps >> n_communities;

        igraph_arpack_options_t options;
        igraph_arpack_options_init(&options);

        igVector membership, finalMembership, eigenvalues;
        igMatrix merges;
        igList eigenvectors;
        igraph_real_t modularity;

        if (n_communities == 0) {
            igCheck(igraph_community_leading_eigenvector(
                        &graph, passWeights(),
                        &merges.mat, &membership.vec,
                        steps, &options, &modularity, false,
                        &eigenvalues.vec, &eigenvectors.list, // eigenvalues, eigenvectors
                        nullptr, // history
                        nullptr, nullptr // callback
                        ));
            finalMembership = membership;
        } else {
            igCheck(igraph_community_leading_eigenvector(
                        &graph, passWeights(),
                        &merges.mat, &membership.vec /* membership */,
                        steps, &options, nullptr /* modularity */, false,
                        &eigenvalues.vec, &eigenvectors.list, // eigenvalues, eigenvectors
                        nullptr, // history
                        nullptr, nullptr // callback
                        ));
            finalMembership = membership;
            igraph_integer_t cc = 1 + static_cast<igraph_integer_t>( *std::max_element(membership.begin(), membership.end()) );
            igCheck(igraph_le_community_to_membership(&merges.mat, cc - n_communities, &membership.vec, nullptr));
            igCheck(igraph_modularity(&graph, &membership.vec, passWeights(), /* resolution= */ 1, /* directed= */ false, &modularity));
        }

        ml.newPacket();
        ml << mlHead("List", 6) << membership << finalMembership << merges << eigenvalues << eigenvectors << modularity;
    }

    void communityFluid(MLINK link) const {
        mlStream ml{link, "communityFluid"};
        igraph_integer_t nc;
        ml >> mlCheckArgs(1) >> nc;

        igVector membership;
        igraph_real_t modularity;
        igCheck(igraph_community_fluid_communities(&graph, nc, &membership.vec, &modularity));

        ml.newPacket();
        ml << mlHead("List", 2) << membership << modularity;
    }

    // Unfold tree

    // note this is a constructor!
    mma::RealTensorRef unfoldTree(const IG &source, mma::RealTensorRef troots, bool directed) {
        igraph_vector_t roots = igVectorView(troots);
        igVector mapping;
        destroy();
        igConstructorCheck(igraph_unfold_tree(&source.graph, &graph, directed ? IGRAPH_OUT : IGRAPH_ALL, &roots, &mapping.vec));
        return mapping.makeMTensor();
    }

    // Bipartite partitions

    mma::IntTensorRef bipartitePartitions() const {
        igraph_bool_t res;
        igBoolVector map;

        igCheck(igraph_is_bipartite(&graph, &res, &map.vec));
        if (! res)
            throw mma::LibraryError();

        return map.makeMTensor();
    }

    mma::IntTensorRef bipartiteProjection(mma::IntTensorRef typevec, IG &p1, IG &p2) const {
        igBoolVector types(typevec.length());
        std::copy(typevec.begin(), typevec.end(), types.begin());

        igVector mult1, mult2;

        igCheck(igraph_bipartite_projection(&graph, &types.vec, &p1.graph, &p2.graph, &mult1.vec, &mult2.vec, -1));

        auto res = mma::makeVector<mint>(mult1.length() + mult2.length());
        std::copy(mult1.begin(), mult1.end(), res.begin());
        std::copy(mult2.begin(), mult2.end(), res.begin() + mult1.length());

        return res;
    }

    // Contract vertices

    // note this is a constructor!
    void contractVertices(mma::RealTensorRef t) {
        igraph_vector_t mapping = igVectorView(t);
        igCheck(igraph_contract_vertices(&graph, &mapping, nullptr));
    }

    // Random walks

    mma::RealTensorRef randomWalk(mint start, mint steps) const {
        igVector walk;

        /* Use random_edge_walk for the weighted case until random_wakl gets support for weights.
         * This is fine, as random_edge_walk returns the same result as random_walk if the weight
         * vector contains all 1s and the igraph RNG is seeded with the same value.
         */
        if (weightedQ()) {
            if (steps == 0)
                return mma::makeVector<double>(0);

            igCheck(igraph_random_edge_walk(&graph, passWeights(), &walk.vec, start, IGRAPH_OUT, steps-1, IGRAPH_RANDOM_WALK_STUCK_RETURN));

            auto result = mma::makeVector<double>(walk.length() + 1);
            long last = start;
            result[0] = last;
            for (mint i=1; i < result.size(); ++i) {
                last = IGRAPH_OTHER(&graph, walk[i-1], last);
                result[i] = last;
            }
            return result;
        } else {
            igCheck(igraph_random_walk(&graph, &walk.vec, start, IGRAPH_OUT, steps, IGRAPH_RANDOM_WALK_STUCK_RETURN));
            return walk.makeMTensor();
        }
    }

    mma::RealTensorRef randomEdgeWalk(mint start, mint steps) const {
        igVector walk;
        igCheck(igraph_random_edge_walk(&graph, passWeights(), &walk.vec, start, IGRAPH_OUT, steps, IGRAPH_RANDOM_WALK_STUCK_RETURN));
        return walk.makeMTensor();
    }

    // Spanning tree

    mma::RealTensorRef spanningTree() const {
        igVector vector;
        igCheck(igraph_minimum_spanning_tree(&graph, &vector.vec, passWeights()));
        return vector.makeMTensor();
    }

    mma::RealTensorRef randomSpanningTree(mint vid) const {
        igVector vector;
        igCheck(igraph_random_spanning_tree(&graph, &vector.vec, vid));
        return vector.makeMTensor();
    }

    // Vertex colouring

    mma::IntTensorRef vertexColoring() const {
        igIntVector igcolor;
        igCheck(igraph_vertex_coloring_greedy(&graph, &igcolor.vec, IGRAPH_COLORING_GREEDY_COLORED_NEIGHBORS));
        return igcolor.makeMTensor();
    }

    // Coreness

    mma::RealTensorRef coreness(mint mode) const {
        igVector cores;
        igraph_neimode_t m;
        switch (mode) {
        case 1:
            m = IGRAPH_OUT; break;
        case -1:
            m = IGRAPH_IN; break;
        case 0:
            m = IGRAPH_ALL; break;
        default:
            throw mma::LibraryError("coreness: Invalid mode parameter.");
        }
        igCheck(igraph_coreness(&graph, &cores.vec, m));
        return cores.makeMTensor();
    }

    // Eulerian paths and cycles

    bool eulerianQ(bool closed) const {
        igraph_bool_t has_path, has_cycle;
        igCheck(igraph_is_eulerian(&graph, &has_path, &has_cycle));
        return closed ? has_cycle : has_path;
    }

    mma::RealTensorRef eulerianPath(bool closed) const {
        igVector edges;

        if (closed)
            igCheck(igraph_eulerian_cycle(&graph, &edges.vec, nullptr));
        else
            igCheck(igraph_eulerian_path(&graph, &edges.vec, nullptr));

        return edges.makeMTensor();
    }

    mma::RealTensorRef eulerianPathVertices(bool closed) const {
        igVector vertices;

        if (closed)
            igCheck(igraph_eulerian_cycle(&graph, nullptr, &vertices.vec));
        else
            igCheck(igraph_eulerian_path(&graph, nullptr, &vertices.vec));

        return vertices.makeMTensor();
    }

    // Other

    mma::IntTensorRef treelikeComponents() const {
        std::vector<mint> comps;

        igraph_adjlist_t al;
        // This is the only igraph function used in treelikeComponents() that can fail.
        // Thus a simple igCheck() will do, no need to worry about freeing other igraph-specific data structures.
        igCheck(igraph_adjlist_init(&graph, &al, IGRAPH_ALL, IGRAPH_NO_LOOPS, IGRAPH_MULTIPLE));

        mint vcount = vertexCount();
        std::vector<mint> degrees;
        degrees.reserve(vcount);
        for (mint i=0; i < vcount; ++i) {
            mint deg = igraph_vector_int_size(igraph_adjlist_get(&al, i));
            degrees.push_back(deg);
            if (deg == 0)
                comps.push_back(i);
        }

        // finddegree-1 nodes, and update the degree vector as if they were deleted
        // repeat until nothing changes
        bool changed;
        do {
            changed = false;
            for (mint i=0; i < vcount; ++i) {
                if (degrees[i] == 1) {
                    changed = true;
                    degrees[i]--;
                    comps.push_back(i);
                    igraph_vector_int_t *neighbors = igraph_adjlist_get(&al, i);
                    mint ncount = igraph_vector_int_size(neighbors);
                    for (long j=0; j < ncount; ++j) {
                        igraph_integer_t v = VECTOR(*neighbors)[j];
                        if (degrees[v] > 0) {
                            degrees[v]--;
                            if (degrees[v] == 0)
                                comps.push_back(v);
                        }
                    }
                }
            }
        } while (changed);

        igraph_adjlist_destroy(&al);

        return mma::makeVector<mint>(comps.size(), comps.data());
    }


    // note that this is a constructor!
    mma::IntTensorRef smoothen(const IG &ig) {
        struct edge {
            igraph_integer_t v1, v2;
            double weight;
            int key; // the key serves to distinguish parallel eges

            bool operator < (const edge &e) const {
                return std::tie(v1, v2, weight, key) < std::tie(e.v1, e.v2, e.weight, e.key);
            }
        };

        mint vcount = ig.vertexCount();
        mint ecount = ig.edgeCount();

        std::vector<std::set<edge>> inclist(vcount);

        // build a vertex-edge incidence list
        for (int i=0; i < ecount; ++i) {
            edge e = {
                IGRAPH_FROM(&ig.graph, i),
                IGRAPH_TO(&ig.graph, i),
                ig.weightedQ() ? ig.weights[i] : 1.0,
                i
            };

            inclist[e.v1].insert(e);
            inclist[e.v2].insert(e);
        }

        // This vector records the vertices that are smoothed out, so that
        // later they can be deleted in pure Mathematica code.
        std::vector<mint> deletedVertices;

        int key = ecount;
        for (igraph_integer_t v=0; v < inclist.size(); ++v) {
            auto &es = inclist[v];
            // smooth degree-2 vertices
            if (es.size() == 2) {
                auto it = es.begin();
                edge e1 = *it++;
                edge e2 = *it++;
                if (e1.v1 == e1.v2 || e2.v1 == e2.v2)
                    continue; // skip self loops

                es.clear();

                igraph_integer_t v1 = e1.v1 != v ? e1.v1 : e1.v2;
                igraph_integer_t v2 = e2.v1 != v ? e2.v1 : e2.v2;

                edge combined = {v1, v2, e1.weight + e2.weight, key++};

                inclist[v1].erase(e1);
                inclist[v1].insert(combined);

                inclist[v2].erase(e2);
                inclist[v2].insert(combined);

                deletedVertices.push_back(v);
            }
        }

        std::set<edge> result;
        for (const auto &s : inclist)
            result.insert(s.begin(), s.end());

        igVector igedges(2*result.size());
        {
            mint i=0;
            for (const auto &e : result) {
                igedges[i++] = e.v1;
                igedges[i++] = e.v2;
            }
        }

        fromEdgeList(igedges, vcount, false);

        weighted = true;
        weights.resize(result.size());
        {
            mint i=0;
            for (const auto &e : result)
                weights[i++] = e.weight;
        }

        return mma::makeVector<mint>(deletedVertices.size(), deletedVertices.data());
    }


    mma::IntTensorRef smoothenDirected(const IG &ig) {
        struct edge {
            igraph_integer_t v1, v2;
            double weight;
            int key; // the key serves to distinguish parallel eges

            bool operator < (const edge &e) const {
                return std::tie(v1, v2, weight, key) < std::tie(e.v1, e.v2, e.weight, e.key);
            }
        };

        mint vcount = ig.vertexCount();
        mint ecount = ig.edgeCount();

        std::vector<std::set<edge>> inclist_out(vcount);
        std::vector<std::set<edge>> inclist_in(vcount);

        // build a vertex-edge incidence list
        for (int i=0; i < ecount; ++i) {
            edge e = {
                IGRAPH_FROM(&ig.graph, i),
                IGRAPH_TO(&ig.graph, i),
                ig.weightedQ() ? ig.weights[i] : 1.0,
                i
            };

            inclist_out[e.v1].insert(e);
            inclist_in[e.v2].insert(e);
        }

        // This vector records the vertices that are smoothed out, so that
        // later they can be deleted in pure Mathematica code.
        std::vector<mint> deletedVertices;

        int key = ecount;
        for (igraph_integer_t v=0; v < inclist_out.size(); ++v) {
            auto &es_out = inclist_out[v];
            auto &es_in = inclist_in[v];
            // smooth degree-2 vertices
            if (es_out.size() == 1 && es_in.size() == 1) {
                edge e_out = *es_out.begin();
                edge e_in = *es_in.begin();

                if (e_in.v1 == e_in.v2 || e_out.v1 == e_out.v2)
                    continue; // skip self-loops

                es_in.clear();
                es_out.clear();

                igraph_integer_t v_out = e_out.v2;
                igraph_integer_t v_in = e_in.v1;

                edge combined = {v_in, v_out, e_in.weight + e_out.weight, key++};

                inclist_out[v_in].erase(e_in);
                inclist_out[v_in].insert(combined);

                inclist_in[v_out].erase(e_out);
                inclist_in[v_out].insert(combined);

                deletedVertices.push_back(v);
            }
        }

        std::set<edge> result;
        for (const auto &s : inclist_out)
            result.insert(s.begin(), s.end());
        for (const auto &s : inclist_in)
            result.insert(s.begin(), s.end());

        igVector igedges(2*result.size());
        {
            mint i=0;
            for (const auto &e : result) {
                igedges[i++] = e.v1;
                igedges[i++] = e.v2;
            }
        }

        fromEdgeList(igedges, vcount, true);

        weighted = true;
        weights.resize(result.size());
        {
            mint i=0;
            for (const auto &e : result)
                weights[i++] = e.weight;
        }

        return mma::makeVector<mint>(deletedVertices.size(), deletedVertices.data());
    }


    // warning: this function modifies the graph (converts to a double-connected directed representation).
    mma::IntTensorRef coordinatesToEmbedding(mma::RealMatrixRef coord) {
        igCheck(igraph_to_undirected(&graph, IGRAPH_TO_UNDIRECTED_COLLAPSE, nullptr));
        igCheck(igraph_simplify(&graph, true, true, nullptr));
        igCheck(igraph_to_directed(&graph, IGRAPH_TO_DIRECTED_MUTUAL));

        long ecount = edgeCount();
        long vcount = vertexCount();

        if (coord.cols() != 2)
            throw mma::LibraryError("coordinatesToEmbedding: Two-dimensional coordinates expected.");

        if (coord.rows() != vcount)
            throw mma::LibraryError("coordinatesToEmbedding: The number of coordinate-pairs should be the same as the vertex count.");

        igraph_inclist_t inclist;
        igraph_inclist_init(&graph, &inclist, IGRAPH_OUT, IGRAPH_NO_LOOPS);

        std::vector<double> angles(ecount);

        for (long i=0; i < ecount; ++i) {
            igraph_integer_t from = IGRAPH_FROM(&graph, i), to = IGRAPH_TO(&graph, i);

            angles[i] = std::atan2(coord(to,1) - coord(from,1), coord(to,0) - coord(from,0));
        }

        for (long i=0; i < vcount; ++i) {
            igraph_vector_int_t *edges = igraph_inclist_get(&inclist, i);

            std::sort(VECTOR(*edges), VECTOR(*edges) + igraph_vector_int_size(edges),
                      [&] (igraph_integer_t e1, igraph_integer_t e2) { return angles[e1] < angles[e2]; } );
        }

        std::vector<std::vector<mint>> emb(vcount);
        for (long i=0; i < vcount; ++i) {
            igraph_vector_int_t *edges = igraph_inclist_get(&inclist, i);

            long ncount = igraph_vector_int_size(edges);
            emb[i].reserve(ncount);
            for (long j=0; j < ncount; ++j)
                emb[i].push_back(IGRAPH_TO(&graph, VECTOR(*edges)[j]));
        }

        igraph_inclist_destroy(&inclist);

        mint total_length = 1;
        for (const auto &el : emb)
            total_length += el.size() + 1;

        auto result = mma::makeVector<mint>(total_length);
        result[0] = emb.size();
        mint k=1;
        for (const auto &el : emb)
            result[k++] = el.size();
        for (const auto &el : emb) {
            std::copy(el.begin(), el.end(), result.data() + k);
            k += el.size();
        }

        return result;
    }


    // this function is used in the implementation of others
    // this is a constructor!
    void complement(const IG &ig) {
        destroy();
        igConstructorCheck(igraph_complementer(&graph, &ig.graph, /* loops = */ false));
    }


    // check if a graph is perfect
    // this function assumes the input to be undirected, simple, connected
    bool perfectQ() { /* TODO: re-add const once ladSubisomoprhic could be made const */
        // bipartite graphs and chordal graphs are perfect
        if (bipartiteQ() || chordalQ())
            return true;

        // possibly more optimizations found here: http://www.or.uni-bonn.de/~hougardy/paper/ClassesOfPerfectGraphs.pdf

        IG comp;
        comp.complement(*this);

        // weak perfect graph theorem:
        // a graph is perfect iff its complement is perfect
        if (comp.bipartiteQ() || comp.chordalQ())
            return true;

        mint g = girth();
        if (g > 3 && g % 2 == 1)
            return false;
        mint cg = comp.girth();
        if (cg > 3 && cg % 2 == 1)
            return false;

        // since bipartiteQ() also catches trees, at this point g and cg are both at least 3.
        // for trees, their value would have been 0

        // strong perfect graph theorem:
        // a graph is perfect iff neither it or its complement contains an induced odd cycle of length >= 5
        mint vcount = vertexCount();
        mint start = std::min(g, cg);
        start = start % 2 == 0 ? start+1 : start+2;
        IG cycle;
        for (mint s = start; s <= vcount; s += 2) {
            igVector edges(2*s);
            for (int i=0; i < s; ++i) {
                edges[2*i] = i;
                edges[2*i + 1] = i+1 < s ? i+1 : 0;
            }

            cycle.fromEdgeList(edges, s, false);

            if (s > g && ladSubisomorphic(cycle, /* induced = */ true))
                return false;

            if (s > cg && comp.ladSubisomorphic(cycle, /* induced = */ true))
                return false;

            mma::check_abort();
        }

        return true;
    }

    // The input is expected to be a simple undirected graph
    bool cactusQ() const {
        // The below edge and vertex count checks are valid for simple graphs

        // An upper bound on the number of edges E of a cactus is
        // E <= 3*(V-1)/2 where V is the number of vertices.
        if (edgeCount() > 3*(vertexCount() - 1)/2)
            return false;

        // A lower bound is E >= V-1
        // We do not consider the null graph a cactus
        if (edgeCount() < vertexCount() - 1)
            return  false;

        // A cactus is connected
        if (! connectedQ(/* strong= */ false))
            return false;

        igList edge_comps, vertex_comps;
        igCheck(igraph_biconnected_components(&graph, nullptr, nullptr, &edge_comps.list, &vertex_comps.list, nullptr));

        // We verify that each bi-connected component is either K_2 or an induced cycle.
        // Knowing that the component is biconnected, it is sufficient to check that E <= V.
        long no_comps = vertex_comps.size();
        for (long i=0; i < no_comps; ++i) {
            igraph_integer_t vcount = igraph_vector_size(vertex_comps[i]);
            igraph_integer_t ecount = igraph_vector_size(edge_comps[i]);

            if (ecount > vcount)
                return false;
        }

        return true;
    }

    bool nonSimpleCactusQ() {
        igraph_simplify(&graph, /* multiple= */ false, /* loops= */ true, nullptr);

        // An upper bound on the number of edges when multi-edges are allowed is E <= 2*(V-1)
        // The cactus with the most edges is like an a tree with each edge having multiplicity 2.
        if (edgeCount() > 2*(vertexCount() - 1))
            return false;

        // A lower bound is E >= V-1
        // We do not consider the null graph a cactus
        if (edgeCount() < vertexCount() - 1)
            return  false;

        // A cactus is connected
        if (! connectedQ(/* strong= */ false))
            return false;

        igList edge_comps, vertex_comps;
        igCheck(igraph_biconnected_components(&graph, nullptr, nullptr, &edge_comps.list, &vertex_comps.list, nullptr));

        // We verify that each bi-connected component is either K_2 or an induced cycle.
        // Knowing that the component is biconnected, it is sufficient to check that E <= V.
        long no_comps = vertex_comps.size();
        for (long i=0; i < no_comps; ++i) {
            igraph_integer_t vcount = igraph_vector_size(vertex_comps[i]);
            igraph_integer_t ecount = igraph_vector_size(edge_comps[i]);

            if (ecount > vcount)
                return false;
        }

        return true;
    }

    mma::IntTensorRef strahlerNumber() const {

        // special case for 1-vertex graph, which is not directedQ()
        if (vertexCount() == 1) {
            auto res = mma::makeVector<mint>(1);
            res[0] = 1;
            return res;
        }

        if (! directedQ() || ! treeQ(IGRAPH_OUT))
            throw mma::LibraryError("strahlerNumber: the graph is not a directed out-tree.");

        igVector degree;
        igCheck(igraph_degree(&graph, &degree.vec, igraph_vss_all(), IGRAPH_IN, false));

        // Find the root---this will not go out of bounds since we already verified that the graph is an out-tree.
        long root;
        for (root=0; root < degree.length(); ++root)
            if (degree[root] == 0)
                break;

        igVector postorder, parent;
        igCheck(igraph_dfs(&graph, root, IGRAPH_OUT, false, nullptr, &postorder.vec, &parent.vec, nullptr, nullptr, nullptr, nullptr));

        auto res = mma::makeVector<mint>(postorder.length());
        std::fill(res.begin(), res.end(), 1);

        for (long c : postorder) {
            long p = parent[c];
            if (res[p] <= res[c])
                res[p] = res[c] + 1;
        }

        return res;
    }

    // The input is expected to be an undirected simple graph
    // Computes the intersection array {b, c} of an undirected distance regular graph.
    // If the graph is not distance regular, {} is returned.
    // Non-connected graphs are not excluded.
    void intersectionArray(MLINK link) const;

    void sirProcess(MLINK link) const {
        mlStream ml{link, "sirProcess"};

        double beta, gamma;
        igraph_integer_t n;

        ml >> mlCheckArgs(3) >> beta >> gamma >> n;

        igPtrVector<igraph_sir_t, igraph_sir_destroy> sirList;

        igCheck(igraph_sir(&graph, beta, gamma, n, &sirList.list));

        ml.newPacket();

        ml << mlHead("List", sirList.size());
        for (const auto &sir : sirList) {
            ml << mlHead("List", 4);
            ml << sir->times << sir->no_i << sir->no_r << sir->no_s;
        }
    }

    // Let k be the number of edges within the neighbourhood of a vertex
    // and m be the number of edges of its neighbourhood which connect to
    // outside of this neighbourhood. Then the local density is k / (m + k).
    mma::RealTensorRef localDensity() const {
        // TODO self-loop and multi-edge handling.

        igraph_adjlist_t al;
        igraph_adjlist_init(&graph, &al, IGRAPH_ALL, IGRAPH_NO_LOOPS, IGRAPH_MULTIPLE);

        mint n = vertexCount();

        auto res = mma::makeVector<double>(n);

        std::vector<mint> nei_flags(n);

        for (mint i=0; i < n; ++i) {
            mint int_count = 0, ext_count = 0;

            igraph_vector_int_t *i_neis = igraph_adjlist_get(&al, i);
            mint di = igraph_vector_int_size(i_neis);

            // mark neighbours of i
            for (mint j=0; j < di; ++j) {
                nei_flags[ VECTOR(*i_neis)[j] ] = i+1;
            }
            for (mint j=0; j < di; ++j) {
                auto v = VECTOR(*i_neis)[j];
                igraph_vector_int_t *v_neis = igraph_adjlist_get(&al, v);
                mint dv = igraph_vector_int_size(v_neis);

                for (mint k=0; k < dv; ++k) {
                    mint u = VECTOR(*v_neis)[k];

                    if (nei_flags[u] == i+1) {
                        int_count += 1;
                    } else if (u == i) {
                        int_count += 2;
                    } else {
                        ext_count++;
                    }
                }
            }

            massert(int_count % 2 == 0);
            int_count /= 2;

            res[i] = ext_count == 0 ? 0.0 : double(int_count) / double(int_count + ext_count);
        }

        return res;
    }

    double assortativityNominal(mma::RealTensorRef types, bool directed) const {
        igraph_vector_t types_vec = igVectorView(types);
        igraph_real_t res;

        igCheck(igraph_assortativity_nominal(&graph, &types_vec, &res, directed));

        return res;
    }

    double assortativityValued(mma::RealTensorRef values, mma::RealTensorRef values_in, bool directed) const {
        igraph_vector_t values_vec = igVectorView(values);
        igraph_vector_t values_in_vec = igVectorView(values_in);
        igraph_real_t res;

        igCheck(igraph_assortativity(&graph, &values_vec, directed ? &values_in_vec : nullptr, &res, directed));

        return res;
    }

    double assortativityDegree(bool directed) const {
        igraph_real_t res;

        igCheck(igraph_assortativity_degree(&graph, &res, directed));

        return res;
    }
};

#endif // IG_H
