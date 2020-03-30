/*
 * Copyright (c) 2016-2020 Szabolcs Horvát.
 *
 * See the file LICENSE.txt for copying permission.
 */

#include "IG.h"
#include <algorithm>
#include <unordered_map>

/**** Create (basic) ****/

void IG::fromIncidenceMatrix(mma::SparseMatrixRef<mint> im, bool directed) {
    igVector edgeList(2*im.cols());
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
                throw mma::LibraryError("fromIncidenceMatrix: Invalid incidence matrix.");
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
                throw mma::LibraryError("fromIncidenceMatrix: Invalid incidence matrix.");
            }
        }
    }

    destroy();
    igConstructorCheck(igraph_create(&graph, &edgeList.vec, im.rows() /* vertex count */, directed));
}


void IG::fromEdgeListML(MLINK link) {
    mlStream ml{link, "fromEdgeListML"};
    igMatrix mat;
    igraph_bool_t directed;
    igraph_integer_t n;
    ml >> mlCheckArgs(3) >> n >> directed;
    int argc;
    if (! MLTestHead(link, "Graph", &argc))
        ml.error("Head Graph expected");
    ml >> mlDiscard(1);
    if (! MLTestHead(link, "List", &argc))
        ml.error("Head List expected");
    if (! directed) {
        ml >> mlDiscard(1);
    }
    ml >> mat;

    for (double *v = mat.begin(); v != mat.end(); ++v) {
        (*v) -= 1;
    }

    destroy();
    igConstructorCheck(igraph_create(&graph, &mat.mat.data, n, directed));

    ml.newPacket();
    ml << mlSymbol("Null");
}


/**** Modification ****/

void IG::mycielski() {
    igraph_integer_t vcount = vertexCount();
    igraph_integer_t ecount = edgeCount();

    // special case: for the empty graph, just add a vertex
    if (vcount == 0) {
        igCheck(igraph_add_vertices(&graph, 1, nullptr));
        return;
    }

    // special case: for the singleton graph, just add a vertex and connect it
    if (vcount == 1) {
        igCheck(igraph_add_vertices(&graph, 1, nullptr));
        igCheck(igraph_add_edge(&graph, 0, 1));
        return;
    }

    igCheck(igraph_add_vertices(&graph, vcount+1, nullptr));
    igraph_integer_t w = 2*vcount;

    igVector edges(2*(vcount + 2*ecount));
    igraph_integer_t ec = 0;
    for (igraph_integer_t u=vcount; u < w; ++u) {
        edges[ec++] = w;
        edges[ec++] = u;
    }
    for (igraph_integer_t i=0; i < ecount; ++i) {
        igraph_integer_t v1 = IGRAPH_FROM(&graph, i);
        igraph_integer_t v2 = IGRAPH_TO(&graph, i);
        igraph_integer_t u1 = v1 + vcount;
        igraph_integer_t u2 = v2 + vcount;
        edges[ec++] = u1;
        edges[ec++] = v2;
        edges[ec++] = u2;
        edges[ec++] = v1;
    }
    igCheck(igraph_add_edges(&graph, &edges.vec, nullptr));
}


/**** Testing properties ****/

bool IG::forestQ(mint mode) const {
    igraph_neimode_t imode;
    switch (mode) {
    case 1:
        imode = IGRAPH_OUT; break;
    case 2:
        imode = IGRAPH_IN; break;
    case 3:
        imode = IGRAPH_ALL; break;
    default:
        throw mma::LibraryError("forestQ: invalid mode");
    }

    igGraphList components;

    // connected components with less than 3 vertices are always trees, even if directed
    igraph_decompose(&graph, &components.list, IGRAPH_WEAK, -1, 3);

    bool res = true;

    long n = components.length();
    for (long i=0; i < n; ++i) {
        igraph_bool_t is_tree;
        igraph_is_tree(components[i], &is_tree, nullptr, imode);
        res = res && is_tree;
        if (! res)
            break;
    }

    return res;
}


/**** Isomorphism ****/

mma::RealTensorRef IG::getIsomorphism(IG &ig) {
    emptyMatchDirectedness(ig);

    igraph_bool_t iso;
    igVector map;

    if (multiQ() || ig.multiQ()) {
        IG g1, g2;
        igIntVector vc1, vc2, ec1, ec2;
        g1.createColoredSimpleGraph(*this, vc1, ec1);
        g2.createColoredSimpleGraph(ig, vc2, ec2);

        igCheck(igraph_isomorphic_vf2(
                    &g1.graph, &g2.graph,
                    &vc1.vec, &vc2.vec, &ec1.vec, &ec2.vec,
                    &iso, &map.vec, nullptr, nullptr, nullptr, nullptr));
    } else {
        igCheck(igraph_isomorphic_bliss(&graph, &ig.graph,
                                        nullptr, nullptr, &iso,
                                        &map.vec, nullptr, IGRAPH_BLISS_F, nullptr, nullptr));
    }

    if (iso)
        return map.makeMTensor();
    else
        return mma::makeVector<double>(0);
}


struct MultigraphColors {
    igIntVector vc1, vc2, ec1, ec2;
};


struct ColorReducedMultigraphIsoCompat {
    static igraph_bool_t vertexCompat(
                const igraph_t *, const igraph_t *,
                const igraph_integer_t v1, const igraph_integer_t v2,
                void *arg)
    {
        MultigraphColors *colors = static_cast<MultigraphColors *>(arg);
        return colors->vc1[v1] >= colors->vc2[v2];
    }

    static igraph_bool_t edgeCompat(
                const igraph_t *, const igraph_t *,
                const igraph_integer_t e1, const igraph_integer_t e2,
                void *arg)
    {
        MultigraphColors *colors = static_cast<MultigraphColors *>(arg);
        return colors->ec1[e1] >= colors->ec2[e2];
    }
};


mma::RealTensorRef IG::getSubisomorphism(IG &ig) {
    emptyMatchDirectedness(ig);

    igraph_bool_t iso;
    igVector map;

    if (! (simpleQ() && ig.simpleQ())) {
        IG g1, g2;

        MultigraphColors graph_colors;

        g1.createColoredSimpleGraph(*this, graph_colors.vc1, graph_colors.ec1);
        g2.createColoredSimpleGraph(ig, graph_colors.vc2, graph_colors.ec2);

        igCheck(igraph_subisomorphic_vf2(
                    &g1.graph, &g2.graph,
                    nullptr, nullptr, nullptr, nullptr,
                    &iso, nullptr, &map.vec,
                    &ColorReducedMultigraphIsoCompat::vertexCompat, &ColorReducedMultigraphIsoCompat::edgeCompat, &graph_colors));
    } else {
        igCheck(igraph_subisomorphic_vf2(
                    &graph, &ig.graph, nullptr, nullptr, nullptr, nullptr,
                    &iso, nullptr, &map.vec, nullptr, nullptr, nullptr));
    }

    if (iso)
        return map.makeMTensor();
    else
        return mma::makeVector<double>(0);
}


void IG::vf2FindIsomorphisms(MLINK link) {
    mlStream ml{link, "vf2Isomorphism"};

    mint id; // expression ID
    igIntVector vc1, vc2, ec1, ec2;

    struct VF2data {
        std::list<igVector> list;
        mlint64 remaining; // remaining number of isomorphisms to find, negative value will run until all are found
    } vf2data;

    ml >> mlCheckArgs(6) >> id >> vf2data.remaining >> vc1 >> vc2 >> ec1 >> ec2;

    IG &ig = mma::getInstance<IG>(id);
    emptyMatchDirectedness(ig);

    struct {
        static igraph_bool_t handle(const igraph_vector_t *map12,  const igraph_vector_t * /* map21 */, void *arg) {
            VF2data &data = *static_cast<VF2data *>(arg);
            data.list.push_back(map12);
            data.remaining--;
            return data.remaining != 0; // negative will run until all are found
        }
    } isohandler;

    igCheck(igraph_isomorphic_function_vf2(
                &graph, &ig.graph,
                vc1.length() == 0 ? nullptr : &vc1.vec, vc2.length() == 0 ? nullptr : &vc2.vec,
                ec1.length() == 0 ? nullptr : &ec1.vec, ec2.length() == 0 ? nullptr : &ec2.vec,
                nullptr, nullptr, &isohandler.handle, nullptr, nullptr, &vf2data));

    ml.newPacket();
    ml << vf2data.list;
}


void IG::vf2FindSubisomorphisms(MLINK link) {
    mlStream ml{link, "vf2Isomorphism"};

    mint id; // expression ID
    igIntVector vc1, vc2, ec1, ec2;

    struct VF2data {
        std::list<igVector> list;
        mlint64 remaining; // remaining number of isomorphisms to find, negative value will run until all are found
    } vf2data;

    ml >> mlCheckArgs(6) >> id >> vf2data.remaining >> vc1 >> vc2 >> ec1 >> ec2;

    IG &ig = mma::getInstance<IG>(id);
    emptyMatchDirectedness(ig);

    struct IsoHandler {
        static igraph_bool_t handle(const igraph_vector_t * /* map12 */,  const igraph_vector_t *map21, void *arg) {
            VF2data &data = *static_cast<VF2data *>(arg);
            data.list.push_back(map21);
            data.remaining--;
            return data.remaining != 0; // negative will run until all are found
        }
    };

    igCheck(igraph_subisomorphic_function_vf2(
                &graph, &ig.graph,
                vc1.length() == 0 ? nullptr : &vc1.vec, vc2.length() == 0 ? nullptr : &vc2.vec,
                ec1.length() == 0 ? nullptr : &ec1.vec, ec2.length() == 0 ? nullptr : &ec2.vec,
                nullptr, nullptr, &IsoHandler::handle, nullptr, nullptr, &vf2data));

    ml.newPacket();
    ml << vf2data.list;
}


bool IG::vf2IsomorphicMulti(IG &ig) {
    emptyMatchDirectedness(ig);

    if (directedQ() != ig.directedQ())
        throw mma::LibraryError("Cannot compare directed and undirected graphs.");

    if (vertexCount() != ig.vertexCount() || edgeCount() != ig.edgeCount())
        return false;

    IG g1, g2;
    igIntVector vc1, vc2, ec1, ec2;
    g1.createColoredSimpleGraph(*this, vc1, ec1);
    g2.createColoredSimpleGraph(ig, vc2, ec2);

    igraph_bool_t iso;
    igCheck(igraph_isomorphic_vf2(
                &g1.graph, &g2.graph,
                &vc1.vec, &vc2.vec, &ec1.vec, &ec2.vec,
                &iso, nullptr, nullptr, nullptr, nullptr, nullptr));

    return iso;
}


bool IG::vf2SubisomorphicMulti(IG &ig) {
    emptyMatchDirectedness(ig);

    if (directedQ() != ig.directedQ())
        throw mma::LibraryError("Cannot compare directed and undirected graphs.");

    if (vertexCount() < ig.vertexCount() || edgeCount() < ig.edgeCount())
        return false;

    IG g1, g2;

    MultigraphColors graph_colors;

    g1.createColoredSimpleGraph(*this, graph_colors.vc1, graph_colors.ec1);
    g2.createColoredSimpleGraph(ig, graph_colors.vc2, graph_colors.ec2);

    igraph_bool_t iso;
    igCheck(igraph_subisomorphic_vf2(
                &g1.graph, &g2.graph, nullptr, nullptr, nullptr, nullptr,
                &iso, nullptr, nullptr,
                &ColorReducedMultigraphIsoCompat::vertexCompat, &ColorReducedMultigraphIsoCompat::edgeCompat, &graph_colors));

    return iso;
}


/**** Regularity properties ****/

// The input is expected to be an undirected simple graph
// Computes the intersection array {b, c} of an undirected distance regular graph.
// If the graph is not distance regular, {} is returned.
// Non-connected graphs are not excluded.
void IG::intersectionArray(MLINK link) const {
    mlStream ml{link, "intersectionArray"};
    ml >> mlCheckArgs(0);

    bool distanceRegular = true;
    size_t vcount = vertexCount();

    // special case: null graph
    if (vcount == 0) {
        ml.newPacket();
        ml << mlHead("List", 2)
           << mlHead("List", 0) << mlHead("List", 0);
        return;
    }

    // get shortest path matrix
    // TODO refactor code so that the entire matrix does not need to be kept in memory
    igMatrix dm;
    igCheck(igraph_shortest_paths(&graph, &dm.mat, igraph_vss_all(), igraph_vss_all(), IGRAPH_OUT));

    // compute graph diameter
    igraph_real_t diam = 0;
    for (const auto &el : dm)
        if (el != IGRAPH_INFINITY && el > diam)
            diam = el;

    std::vector<mint> bvec(diam+1, -1), cvec(diam+1, -1);

    igraph_adjlist_t al;
    igCheck(igraph_adjlist_init(&graph, &al, IGRAPH_ALL));

    for (size_t u=0; u < vcount; ++u) {
        for (size_t v=0; v < vcount; ++v) {
            igraph_real_t i = dm(v, u);
            if (i == IGRAPH_INFINITY)
                continue;

            mint b=0, c=0;

            igraph_vector_int_t *neighbors = igraph_adjlist_get(&al, v);
            for (const auto &n : igWrap(*neighbors)) {
                if (dm(n, u) == i-1)
                    c++;
                if (dm(n, u) == i+1)
                    b++;
            }

            if (bvec[i] != b) {
                if (bvec[i] == -1) {
                    bvec[i] = b;
                } else {
                    distanceRegular = false;
                    goto end;
                }
            }

            if (cvec[i] != c) {
                if (cvec[i] == -1) {
                    cvec[i] = c;
                } else {
                    distanceRegular = false;
                    goto end;
                }
            }
        }
    }

end:
    igraph_adjlist_destroy(&al);

    ml.newPacket();
    if (distanceRegular) {
        bvec.pop_back(); // remove last element
        cvec.erase(cvec.begin()); // remove first element
        ml << mlHead("List", 2)
           << bvec << cvec;
    } else {
        ml << mlHead("List", 0);
    }
}


/**** Transitivity properties ****/

// Used for equivalence class (orbit) computation.
// Each orbit element points to another one that it is equivalent to, or to itself.
template<typename T>
class OrbitElement {
    T val;
    OrbitElement *equiv;

    OrbitElement(const OrbitElement &) = delete;
    OrbitElement & operator = (const OrbitElement &) = delete;

public:

    OrbitElement() { equiv = this; }

    void setValue(const T &newVal) { val = newVal; }
    T value() const { return val; }

    OrbitElement *getClassElem() {
        OrbitElement *final = equiv;
        while (final != final->equiv)
            final = final->equiv;

        // If updating is needed, do a full second pass and update each node of the chain.
        if (final != equiv) {
            OrbitElement *e1 = this;
            while (e1 != e1->equiv) {
                OrbitElement *e2 = e1->equiv;
                e1->equiv = final;
                e1 = e2;
            }
        }

        return equiv;
    }

    void updateClass(OrbitElement *elem) {
        OrbitElement *newClass = elem->getClassElem();
        getClassElem()->equiv = newClass;
        equiv = newClass;
    }

};


// Orbits of vertices. Used for vertex transitivity.
class GroupOrbits {

    typedef OrbitElement<mint> Element;

    const mint n;
    Element *elems;

public:

    GroupOrbits(mint n) : n(n) {
        elems = new Element[n];
        for (mint i=0; i < n; ++i)
            elems[i].setValue(i);
    }

    ~GroupOrbits() { delete [] elems; }

    template<typename GeneratorArray>
    void computeOrbits(const GeneratorArray &generators) {
        for (mint i=0; i < n; ++i) {
            auto &el = elems[i];
            for (auto &gen : generators)
                el.updateClass(&elems[ gen[el.value()] ]);
        }
    }

    // How many orbits? Vertex transitive graphs will have one.
    mint orbitCount() {
        std::set<mint> repr;
        for (mint i=0; i < n; ++i) {
            auto &el = elems[i];
            repr.insert(el.getClassElem()->value());
        }
        return repr.size();
    }
};


class GroupEdgeOrbits {

    typedef std::pair<mint, mint> VertexPair;
    typedef OrbitElement<VertexPair> Element;

    const bool directed;
    const mint vcount;
    const mint n;
    Element *elems;

    // Using std::map<VertexPair, mint> here is slower.
    // We could use std::unordered_map<VertexPair, mint, CustomHash> with a naive hash combiner,
    // but that does not gain any noticeable performance.
    // Instead, we convert each ordered pair into its index in the flattened adjacency matrix
    // (see pairToInt()) and use that integer for lookup.
    std::unordered_map<mint, mint> index;

    mint pairToInt(const VertexPair &p) const { return p.first * vcount + p.second; }

    // return the index of a pair in the elems array
    mint pairIndex(const VertexPair &p) const {
        return index.at(pairToInt(p));
    }

public:

    GroupEdgeOrbits(const igraph_t *graph) :
        directed(igraph_is_directed(graph)),
        vcount(igraph_vcount(graph)),
        n(directed ? igraph_ecount(graph) : 2*igraph_ecount(graph))
    {
        elems = new Element[n];

        igraph_eit_t eit;
        igCheck(igraph_eit_create(graph, igraph_ess_all(IGRAPH_EDGEORDER_ID), &eit));

        mint i=0;

        for (IGRAPH_EIT_RESET(eit); ! IGRAPH_EIT_END(eit); IGRAPH_EIT_NEXT(eit)) {
            igraph_integer_t e = IGRAPH_EIT_GET(eit);
            elems[i].setValue(VertexPair{IGRAPH_FROM(graph, e), IGRAPH_TO(graph, e)});
            index.insert({pairToInt(elems[i].value()), i});
            i++;
        }

        // for directed, we are done
        if (directed)
            goto end;

        // for undirected, also set reverse edges
        for (IGRAPH_EIT_RESET(eit); ! IGRAPH_EIT_END(eit); IGRAPH_EIT_NEXT(eit)) {
            igraph_integer_t e = IGRAPH_EIT_GET(eit);
            elems[i].setValue(VertexPair{IGRAPH_TO(graph, e), IGRAPH_FROM(graph, e)});
            index.insert({pairToInt(elems[i].value()), i});
            i++;
        }

    end:
        igraph_eit_destroy(&eit);
    }

    ~GroupEdgeOrbits() { delete [] elems; }

    template<typename GeneratorArray>
    void computeOrbits(const GeneratorArray &generators) {
        for (mint i=0; i < n; ++i) {
            mma::check_abort();
            auto &el = elems[i];
            for (auto &gen : generators) {
                VertexPair p = {gen[el.value().first], gen[el.value().second]};
                el.updateClass(&elems[ pairIndex(p) ]);
            }
        }
    }

    std::set<VertexPair> orbitRepresentatives() {
        std::set<VertexPair> repr;
        for (mint i=0; i < n; ++i) {
            auto &el = elems[i];
            repr.insert(el.getClassElem()->value());
        }
        return repr;
    }

    VertexPair getRepresentative(const VertexPair &p) {
        return elems[pairIndex(p)].getClassElem()->value();
    }
};


// Orbits of vertex pairs. Used for distance transitivity.
class GroupPairOrbits {

    typedef std::pair<mint, mint> VertexPair;
    typedef OrbitElement<VertexPair> Element;

    const mint n;
    Element *elems;

    // return the index of a pair in the elems array
    mint pairIndex(const VertexPair &p) const {
        return n*p.first + p.second;
    }

public:

    GroupPairOrbits(mint n) : n(n) {
        elems = new Element[n*n];
        for (mint i=0; i < n; ++i)
            for (mint j=0; j < n; ++j) {
                auto p = VertexPair{i,j};
                elems[pairIndex(p)].setValue(p);
            }
    }

    ~GroupPairOrbits() { delete [] elems; }

    template<typename GeneratorArray>
    void computeOrbits(const GeneratorArray &generators) {
        for (mint i=0; i < n*n; ++i) {
            mma::check_abort();
            auto &el = elems[i];
            for (auto &gen : generators) {
                VertexPair p = {gen[el.value().first], gen[el.value().second]};
                el.updateClass(&elems[ pairIndex(p) ]);
            }
        }
    }

    std::set<VertexPair> orbitRepresentatives() {
        std::set<VertexPair> repr;
        for (mint i=0; i < n*n; ++i) {
            auto &el = elems[i];
            repr.insert(el.getClassElem()->value());
        }
        return repr;
    }
};


// The input is expected not to have multi-edges. Self-loops are permissible.
bool IG::vertexTransitiveQ(mint splitting) const {

    // Handle trivial edge cases
    // vcount == 2 is not necessarily transitive in the directed case
    mint vcount = vertexCount();
    if (vcount <= 1)
        return true;

    // List of automorphism group generators
    igList list;
    igCheck(igraph_automorphism_group(&graph, nullptr, &list.list, blissIntToSplitting(splitting), nullptr));

    mint gcount = list.size();
    // If there are no non-trivial automorphisms,
    // and the graph has more than one vertex,
    // then it is not vertex transitive.
    if (gcount == 0)
        return false;

    // copy generators to std::vectors
    std::vector<std::vector<mint>> generators(gcount);
    for (mint i=0; i < gcount; ++i)
        generators[i].assign(igWrap(*list[i]).cbegin(), igWrap(*list[i]).cend());

    GroupOrbits orbits(vcount);
    orbits.computeOrbits(generators);
    if (orbits.orbitCount() != 1)
        return false;

    return true;
}


// The input is expected not to have multi-edges. Self-loops are permissible.
bool IG::edgeTransitiveQ(mint splitting) const {

    // Handle trivial edge cases
    // A graph with no edges or a single edge is always edge transitive
    if (edgeCount() <= 1)
        return true;

    // List of automorphism group generators
    igList list;
    igCheck(igraph_automorphism_group(&graph, nullptr, &list.list, blissIntToSplitting(splitting), nullptr));

    mint gcount = list.size();
    // If there are no non-trivial automorphisms,
    // and the graph has more than one edge,
    // then it cannot be edge transitive.
    if (gcount == 0)
        return false;

    // copy generators to std::vectors
    std::vector<std::vector<mint>> generators(gcount);
    for (mint i=0; i < gcount; ++i)
        generators[i].assign(igWrap(*list[i]).cbegin(), igWrap(*list[i]).cend());

    GroupEdgeOrbits orbits(&graph);
    orbits.computeOrbits(generators);
    auto repr = orbits.orbitRepresentatives();

    // More than 2 classes => not edge transitive
    if (repr.size() > 2)
        return false;

    // Precisely one class <=> arc transitive
    if (repr.size() == 1)
        return true;

    // If we reach here, then repr.size() == 2
    massert(repr.size() == 2);

    // If graph is not arc transitive but it is directed, then it is not edge transitive
    if (directedQ())
        return false;

    // In the undirected case, there may be two classes, the first containing all
    // connected pairs, and the second containing precisely the reverse of these pairs.
    // We take the representative of the first class, and check if its reverse is
    // in the other class. Then the graph is edge transitive.
    auto edge = *repr.begin();
    if (orbits.getRepresentative({edge.second, edge.first}) != edge)
        return true;
    else
        return false;
}


// "Symmetric" here means both vertex transitive and edge transitive.
// This function computes the automorphism group only once.
// The input is expected not to have multi-edges. Self-loops are permissible.
bool IG::symmetricQ(mint splitting) const {

    // Handle trivial edge cases
    // The null graph, singleton graph and single vertex with self-loop are symmetric
    mint vcount = vertexCount();
    if (vcount <= 1)
        return true;

    if (edgeCount() == 0)
        return true;

    // List of automorphism group generators
    igList list;
    igCheck(igraph_automorphism_group(&graph, nullptr, &list.list, blissIntToSplitting(splitting), nullptr));

    mint gcount = list.size();
    // If there are no non-trivial automorphisms,
    // and the graph has more than one vertex,
    // then it cannot be vertex transitive
    if (gcount == 0)
        return false;

    // copy generators to std::vectors
    std::vector<std::vector<mint>> generators(gcount);
    for (mint i=0; i < gcount; ++i)
        generators[i].assign(igWrap(*list[i]).cbegin(), igWrap(*list[i]).cend());

    // Check vertex transitivity.
    {
        GroupOrbits orbits(vcount);
        orbits.computeOrbits(generators);
        if (orbits.orbitCount() != 1)
            return false;
    }

    // If the graph is vertex transitive, also check edge transitivity.
    // Warning: Even if the graph is vertex and edge transitive,
    // connected vertex pairs may have two orbits (Doyle graph)
    {
        GroupEdgeOrbits orbits(&graph);
        orbits.computeOrbits(generators);
        auto repr = orbits.orbitRepresentatives();

        // More than 2 classes => not edge transitive
        if (repr.size() > 2)
            return false;

        // Precisely one class <=> arc transitive
        if (repr.size() == 1)
            return true;

        // If we reach here, then repr.size() == 2
        massert(repr.size() == 2);

        // If graph is not arc transitive but it is directed, then it is not edge transitive
        if (directedQ())
            return false;

        // In the undirected case, there may be two classes, the first containing all
        // connected pairs, and the second containing precisely the reverse of these pairs.
        // We take the representative of the first class, and check if its reverse is
        // in the other class. Then the graph is edge transitive.
        auto edge = *repr.begin();
        if (orbits.getRepresentative({edge.second, edge.first}) != edge)
            return true;
        else
            return false;
    }
}


// The input is expected not to have multi-edges. Self-loops are permissible.
bool IG::distanceTransitiveQ(mint splitting) const {

    // Handle trivial edge cases
    // vcount == 2 is not necessarily transitive in the directed case
    mint vcount = vertexCount();
    if (vcount <= 1)
        return true;

    // List of automorphism group generators
    igList list;
    igCheck(igraph_automorphism_group(&graph, nullptr, &list.list, blissIntToSplitting(splitting), nullptr));

    // If there are no non-trivial automorphisms,
    // and the graph has more than one vertex,
    // then it is not vertex or distance transitive.
    mint gcount = list.size();
    if (gcount == 0) // no non-trivial automorphisms
        return false;

    // massert(vcount == igraph_vector_size(*list.begin()));

    // Copy generators to std::vectors
    std::vector<std::vector<mint>> generators(gcount);
    for (mint i=0; i < gcount; ++i)
        generators[i].assign(igWrap(*list[i]).cbegin(), igWrap(*list[i]).cend());

    {
        GroupOrbits orbits(vcount);
        orbits.computeOrbits(generators);
        if (orbits.orbitCount() != 1)
            return false;
    }

    std::set<std::pair<mint,mint>> repr;

    {
        // Pair orbits calculation requires vcount^2 memory
        // Keep it in its own block so the memory is released as soon as it is not needed
        // The distance matrix calculation below also requires vcount^2 memory
        GroupPairOrbits orbits(vcount);
        orbits.computeOrbits(generators);
        repr = orbits.orbitRepresentatives();
    }

    // TODO: Do not compute *all*-pairs shortest paths as it is not really needed.
    // Not a high priority because the orbit calculation already takes vcount^2 memory,
    // and the automorphism group computation usually takes longer than the distance matrix computation.
    igMatrix dm;
    igCheck(igraph_shortest_paths(&graph, &dm.mat, igraph_vss_all(), igraph_vss_all(), IGRAPH_OUT));

    std::set<igraph_real_t> ds;
    for (const auto &el : repr) {
        igraph_real_t dist = dm(el.first, el.second);
        if (ds.find(dist) != ds.end())
            return false;
        ds.insert(dist);
    }

    return true;
}
