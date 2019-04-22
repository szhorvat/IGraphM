/*
 * Copyright (c) 2019 Szabolcs Horvát.
 *
 * See the file LICENSE.txt for copying permission.
 */

#include "IG.h"

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


/**** Testing ****/

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
        g1.coloredSimpleGraph(*this, vc1, ec1);
        g2.coloredSimpleGraph(ig, vc2, ec2);

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

        g1.coloredSimpleGraph(*this, graph_colors.vc1, graph_colors.ec1);
        g2.coloredSimpleGraph(ig, graph_colors.vc2, graph_colors.ec2);

        igCheck(igraph_subisomorphic_vf2(
                    &g1.graph, &g2.graph,
                    nullptr, nullptr, nullptr, nullptr,
                    &iso, &map.vec, nullptr,
                    &ColorReducedMultigraphIsoCompat::vertexCompat, &ColorReducedMultigraphIsoCompat::edgeCompat, &graph_colors));
    } else {
        igCheck(igraph_subisomorphic_vf2(
                    &graph, &ig.graph, nullptr, nullptr, nullptr, nullptr,
                    &iso, &map.vec, nullptr, nullptr, nullptr, nullptr));
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
    g1.coloredSimpleGraph(*this, vc1, ec1);
    g2.coloredSimpleGraph(ig, vc2, ec2);

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

    g1.coloredSimpleGraph(*this, graph_colors.vc1, graph_colors.ec1);
    g2.coloredSimpleGraph(ig, graph_colors.vc2, graph_colors.ec2);

    igraph_bool_t iso;
    igCheck(igraph_subisomorphic_vf2(
                &g1.graph, &g2.graph, nullptr, nullptr, nullptr, nullptr,
                &iso, nullptr, nullptr,
                &ColorReducedMultigraphIsoCompat::vertexCompat, &ColorReducedMultigraphIsoCompat::edgeCompat, &graph_colors));

    return iso;
}
