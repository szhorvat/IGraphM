(* Mathematica Source File *)
(* Created by Mathematica Plugin for IntelliJ IDEA *)

(* :Author: szhorvat *)
(* :Date: 2019-05-01 *)
(* :Copyright: (c) 2019 Szabolcs Horvát *)

Package["IGraphM`"]
igContextSetup[igPackagePrivateSymbol]

(**********************)
(***** Regularity *****)
(**********************)

PackageExport["IGRegularQ"]
IGRegularQ::usage =
    "IGRegularQ[graph] tests if graph is regular, i.e. all vertices have the same degree.\n" <>
    "IGRegularQ[graph, k] tests if graph is k-regular, i.e. all vertices have degree k.";
SyntaxInformation[IGRegularQ] = {"ArgumentsPattern" -> {_, _.}};
IGRegularQ[graph_?igGraphQ] :=
    If[UndirectedGraphQ[graph],
      Equal @@ VertexDegree[graph],
      Equal @@ VertexOutDegree[graph] && VertexInDegree[graph] == VertexOutDegree[graph]
    ]
IGRegularQ[graph_?EmptyGraphQ, k_Integer?NonNegative] := k == 0 (* required because First@VertexDegree[g] is not usable for the null graph *)
IGRegularQ[graph_?igGraphQ, k_Integer?NonNegative] :=
    If[UndirectedGraphQ[graph],
      With[{vd = VertexDegree[graph]},
        First[vd] == k && Equal @@ vd
      ]
      ,
      With[{vod = VertexOutDegree[graph], vid = VertexInDegree[graph]},
        First[vod] == k && Equal @@ vod && vid == vod
      ]
    ]
IGRegularQ[_] := False


PackageExport["IGStronglyRegularQ"]
IGStronglyRegularQ::usage = "IGStronglyRegularQ[graph] tests if graph is strongly regular.";

igStronglyRegularQ[g_?EmptyGraphQ] := True
igStronglyRegularQ[g_ /; UndirectedGraphQ[g] && SimpleGraphQ[g]] :=
    IGRegularQ[g] && (IGCompleteQ[g] ||
    Module[{vc = VertexCount[g], am = AdjacencyMatrix[g], k = First@VertexDegree[g], lambda, mu},
      lambda = Dot @@ am[[ First[am["NonzeroPositions"]] ]];
      mu = k*(k - lambda - 1) / (vc - k - 1);
      If[
        IntegerQ[mu]
        (* The optimization below seems to have minimal effect, thus for the moment it is excluded.
           Source of the formula: "Strongly regular graphs" by Peter J. Cameron *)
        (* && IntegerQ[1/2 (vc - 1 + ((vc - 1) (mu - lambda) - 2 k)/ Sqrt[(mu - lambda)^2 + 4 (k - mu)])] *)
        ,
        am.am == k IdentityMatrix[vc, SparseArray] + lambda am + mu (ConstantArray[1, {vc, vc}, SparseArray] - IdentityMatrix[vc, SparseArray] - am)
        ,
        False
      ]
    ])

(* TODO implement for directed graphs including parameters
   https://homepages.cwi.nl/~aeb/math/dsrg/dsrg.html *)
IGStronglyRegularQ::dirg = "Directed graphs are not supported.";
igStronglyRegularQ[g_?DirectedGraphQ] := (Message[IGStronglyRegularQ::dirg]; $Failed)

IGStronglyRegularQ::nsg  = "Non-simple graphs are not supported.";
igStronglyRegularQ[g_] := (Message[IGStronglyRegularQ::nsg]; $Failed)

SyntaxInformation[IGStronglyRegularQ] = {"ArgumentsPattern" -> {_}};
IGStronglyRegularQ[graph_?igGraphQ] := igStronglyRegularQ[graph]
IGStronglyRegularQ[_] := False


PackageExport["IGStronglyRegularParameters"]
IGStronglyRegularParameters::usage = "IGStronglyRegularParameters[graph] returns the parameters {v, k, \[Lambda], \[Mu]} of a strongly regular graph. For non-strongly-regular graphs {} is returned.";
SyntaxInformation[IGStronglyRegularParameters] = {"ArgumentsPattern" -> {_}};
IGStronglyRegularParameters[g_?EmptyGraphQ] := {VertexCount[g], 0, 0, 0}
IGStronglyRegularParameters[g_ /; SimpleGraphQ[g] && IGCompleteQ[g]] := {VertexCount[g], VertexCount[g]-1, VertexCount[g]-2, 0}
IGStronglyRegularParameters[g_?IGStronglyRegularQ] :=
    Module[{vc = VertexCount[g], am = AdjacencyMatrix[g], k = First@VertexDegree[g], lambda, mu},
      lambda = Dot @@ am[[ First[am["NonzeroPositions"]] ]];
      mu = k * (k - lambda - 1) / (vc - k - 1);
      {vc, k, lambda, mu}
    ]
IGStronglyRegularParameters[g_?GraphQ] := {}


PackageExport["IGIntersectionArray"]
IGIntersectionArray::usage = "IGIntersectionArray[graph] computes the intersection array {b, c} of a distance regular graph. For non-distance-regular graphs {} is returned.";
IGIntersectionArray::dirg = "Directed graphs are not supported.";
IGIntersectionArray::nsg  = "Non-simple graphs are not supported.";
SyntaxInformation[IGIntersectionArray] = {"ArgumentsPattern" -> {_}};
IGIntersectionArray[graph_?igGraphQ] :=
    Which[
      UndirectedGraphQ[graph] && SimpleGraphQ[graph],
      If[IGRegularQ[graph],
        catch@Block[{ig = igMakeUnweighted[graph]},
          check@ig@"intersectionArray"[]
        ],
        {}
      ]
      ,
      DirectedGraphQ[graph],
      Message[IGIntersectionArray::dirg]; $Failed
      ,
      True, (* non-simple undirected graph *)
      Message[IGIntersectionArray::nsg]; $Failed
    ]


PackageExport["IGDistanceRegularQ"]
IGDistanceRegularQ::usage = "IGDistanceRegularQ[graph] tests if graph is distance regular.";
SyntaxInformation[IGDistanceRegularQ] = {"ArgumentsPattern" -> {_}};
IGDistanceRegularQ[graph_?igGraphQ] := catch[check@IGIntersectionArray[graph] =!= {}]
IGDistanceRegularQ[_] := False


(************************)
(***** Transitivity *****)
(************************)

PackageExport["IGVertexTransitiveQ"]
IGVertexTransitiveQ::usage = "IGVertexTransitiveQ[graph] tests if graph is vertex transitive.";

IGVertexTransitiveQ::nmg = "Multigraphs are not supported.";
SyntaxInformation[IGVertexTransitiveQ] = {"ArgumentsPattern" -> {_}};
IGVertexTransitiveQ[graph_?EmptyGraphQ] = True;
IGVertexTransitiveQ[graph_?igGraphQ] :=
    If[MultigraphQ[graph],
      Message[IGVertexTransitiveQ::nmg];
      $Failed
      ,
      IGRegularQ[graph] &&
      With[{elems = Range@VertexCount[graph]},
        GroupOrbits[IGBlissAutomorphismGroup[graph], elems] === {elems}
      ]
    ]
IGVertexTransitiveQ[_] = False;


PackageExport["IGEdgeTransitiveQ"]
IGEdgeTransitiveQ::usage = "IGEdgeTransitiveQ[graph] tests if graph is edge transitive.";

IGEdgeTransitiveQ::nmg = IGVertexTransitiveQ::nmg;
SyntaxInformation[IGEdgeTransitiveQ] = {"ArgumentsPattern" -> {_}};
IGEdgeTransitiveQ[graph_?igGraphQ] :=
    If[MultigraphQ[graph],
      Message[IGEdgeTransitiveQ::nmg];
      $Failed
      ,
      IGVertexTransitiveQ@LineGraph[graph]
    ]
IGEdgeTransitiveQ[_] = False;


PackageExport["IGSymmetricQ"]
IGSymmetricQ::usage = "IGSymmetricQ[graph] tests if graph is symmetric, i.e. it is both vertex transitive and edge transitive.";

IGSymmetricQ::nmg = IGVertexTransitiveQ::nmg;
SyntaxInformation[IGSymmetricQ] = {"ArgumentsPattern" -> {_}};
IGSymmetricQ[graph_?igGraphQ] :=
    If[MultigraphQ[graph],
      Message[IGSymmetricQ::nmg];
      $Failed
      ,
      IGVertexTransitiveQ[graph] && IGEdgeTransitiveQ[graph]
    ]
IGSymmetricQ[_] = False;


PackageExport["IGDistanceTransitiveQ"]
IGDistanceTransitiveQ::usage = "IGDistanceTransitiveQ[graph] tests if graph is distance-transitive.";

IGDistanceTransitiveQ::nmg = IGVertexTransitiveQ::nmg;
SyntaxInformation[IGDistanceTransitiveQ] = {"ArgumentsPattern" -> {_}};
(* TODO: Do not compute entire distance matrix unless needed. *)
(* TODO: Implement group orbits in C++ *)
IGDistanceTransitiveQ[graph_?igGraphQ] :=
    If[MultigraphQ[graph],
      Message[IGDistanceTransitiveQ::nmg];
      $Failed
      ,
      IGRegularQ[graph] && (* exclude non-regular graphs early, for performance *)
      Block[{ig = igMakeUnweighted[graph], group, elems}, (* work on a single IG object to avoid translating the graph twice *)
        (* Warning: PermutationGroup should take Cycles expressions, but GroupOrbits does work with permutation lists *)
        group = PermutationGroup@igIndexVec@check@ig@"blissAutomorphismGroup"[0, {}];
        elems = Range@VertexCount[graph];
        GroupOrbits[group, elems] === {elems} && (* exclude non-vertex-transitive graphs early, for performance *)
        With[{dm = Round@fixInfNaN@check@ig@"shortestPaths"[{}, {}]}, (* only compute distance matrix if needed *)
        DuplicateFreeQ[
            Extract[dm, #] & /@ GroupOrbits[group, Tuples[elems, 2]][[All, 1]]
          ]
        ]
      ]
    ]
IGDistanceTransitiveQ[_] = False;