(* Mathematica Package *)
(* Created by Mathematica plugin for IntelliJ IDEA *)

(* :Author: szhorvat *)
(* :Date: 2018-10-24 *)
(* :Copyright: (c) 2019-2020 Szabolcs Horvát *)

Package["IGraphM`"]


PackageImport["HierarchicalClustering`"]

(*******************************)
(***** Community detection *****)
(*******************************)

(* IGClusterData data structure *)

clusterAscQ[asc_?AssociationQ] := AllTrue[{"Elements", "Communities", "Algorithm"}, KeyExistsQ[asc, #]&]
clusterAscQ[_] := False

igClusterDataQ[IGClusterData[asc_]] := clusterAscQ[asc]
igClusterDataQ[_] := False

hierarchicalQ[asc_] := KeyExistsQ[asc, "Merges"]

PackageExport["IGClusterData"]
IGClusterData::usage = "IGClusterData[association] represents the output of community detection functions. Properties can be queried using IGClusterData[\[Ellipsis]][\"property\"].";

IGClusterData::hier   = "The provided clustering is not hierarchical.";
IGClusterData::noprop = "There is no property `` for IGClusterData objects.";
IGClusterData::dnconn = "The dendrogram is not connected. Unconnected dendrograms are not currently supported. Use the \"Merges\" data directly."; (* TODO handle this case *)

IGClusterData[asc_?AssociationQ]["Properties"] :=
    Sort@Join[
      Keys[asc],
      {
        "Properties", "ElementCount",
        If[hierarchicalQ[asc],
          Unevaluated@Sequence["HierarchicalClusters", "Tree"],
          Unevaluated@Sequence[]
        ]
      }
    ]

mergesToHierarchy[asc_] :=
    Module[{cl, elems, merge, cnt, n, s},
      elems = Switch[asc@"Algorithm",
        "LeadingEigenvector", asc@"FinalCommunities",
        _,                    asc@"Elements"
      ];
      n = Length[elems];
      If[Sort@Flatten@asc@"Merges" =!= Range[n + Length@asc@"Merges" - 1],
        Message[IGClusterData::dnconn];
        Return[$Failed]
      ];
      s = Switch[asc@"Algorithm",
        "EdgeBetweenness", 1 + Length@asc@"RemovedEdges" - asc@"Bridges",
        _,                 Range[n]
      ];
      MapThread[Set[cl[#2], #1]&,
        {elems, Range[n]}
      ];
      cnt[Cluster[_, _, _, n1_, n2_]] := n1 + n2;
      cnt[_] := 1;
      merge[{j_, k_}, {i_}] := cl[i + n] = Cluster[cl[j], cl[k], s[[i]], cnt@cl[j], cnt@cl[k]];
      Last@MapIndexed[merge, asc["Merges"]]
    ]

IGClusterData[asc_?AssociationQ]["HierarchicalClusters"] :=
    If[hierarchicalQ[asc],
      mergesToHierarchy[asc],
      Message[IGClusterData::hier]; $Failed
    ]

mergesToTree[asc_] :=
    Module[{mc, root, elems, ec, merges = asc["Merges"], leafIndices, s},
      elems = Switch[asc@"Algorithm",
        "LeadingEigenvector", asc@"FinalCommunities",
        _,                    asc@"Elements"
      ];
      mc = Length[merges];
      ec = Length[elems];
      If[Sort@Flatten[merges] =!= Range[ec + mc - 1],
        Message[IGClusterData::dnconn];
        Return[$Failed]
      ];
      root = ec + mc;
      leafIndices = Intersection[Range[ec], Flatten[merges]];
      s = Switch[asc@"Algorithm",
        "EdgeBetweenness", 1 + Length@asc@"RemovedEdges" - asc@"Bridges",
        _,                 Range[mc]
      ];
      TreeGraph[
        Append[Flatten[merges], root],
        Append[ec + Range[mc] /. n_Integer :> Sequence[n, n], root],
        GraphLayout -> {"LayeredEmbedding", "RootVertex" -> root},
        DirectedEdges -> True,
        EdgeShapeFunction -> "Line",
        EdgeStyle -> Gray,
        VertexShapeFunction -> ({}&), (* Using VertexShapeFunction -> None does not work in M10.0-10.3 due to a bug. *)
        VertexLabels -> Thread[ leafIndices -> (Placed[Short[#], Below]&) /@ elems[[ leafIndices ]] ],
        VertexWeight -> Join[Thread[ec + Range[mc] -> s], Thread[leafIndices -> 0]],
        Properties -> {"GraphProperties" -> {"LeafLabels" -> AssociationThread[leafIndices, elems[[ leafIndices ]] ]}}
      ]
    ]

IGClusterData[asc_?AssociationQ]["Tree"] :=
    If[hierarchicalQ[asc],
      mergesToTree[asc],
      Message[IGClusterData::hier]; $Failed
    ]

IGClusterData[asc_?AssociationQ]["ElementCount"] := Length[asc["Elements"]]

IGClusterData[asc_?AssociationQ][key_String] := Lookup[asc, key, Message[IGClusterData::noprop, key]; $Failed]

IGClusterData[asc_?AssociationQ][keys_List] := IGClusterData[asc] /@ keys

igClusterDataInformation[IGClusterData[asc_?AssociationQ]] :=
    <|
      "Algorithm" -> asc["Algorithm"],
      "ElementCount" -> Length@asc["Elements"],
      "ClusterCount" -> Length@asc["Communities"],
      "Hierarchical" -> hierarchicalQ[asc],
      "Modularity" -> If[KeyExistsQ[asc, "Modularity"], Max@asc["Modularity"], None],
      KeyTake[asc, {"Quality", "CodeLength", "FinalTemperature"}]
    |>

(* Switch to using built-in summary boxes. Since the construction syntax hasn't changed between 10.0-11.2,
   it seems to be relatively safe to use this undocumented functionality *)
IGClusterData /: MakeBoxes[c : IGClusterData[asc_?clusterAscQ], form : (StandardForm|TraditionalForm)] :=
    With[{info = igClusterDataInformation[c]},
      BoxForm`ArrangeSummaryBox[
        IGClusterData,
        c,
        $igClusterDataIcon,
        {
          BoxForm`SummaryItem[{"Elements: ", info["ElementCount"]}],
          BoxForm`SummaryItem[{"Communities: ", info["ClusterCount"]}]
        },
        {
          BoxForm`SummaryItem[{"Community sizes: ", Short[Length /@ asc@"Communities", 0.35]}],
          BoxForm`SummaryItem[{"Modularity: ", Replace[info["Modularity"], None -> "unavailable"]}],
          If[KeyExistsQ[info, "Quality"],
            BoxForm`SummaryItem[{"Quality: ", info@"Quality"}],
            Unevaluated@Sequence[]
          ],
          If[KeyExistsQ[info, "CodeLength"],
            BoxForm`SummaryItem[{"Code length: ", info@"CodeLength"}],
            Unevaluated@Sequence[]
          ],
          BoxForm`SummaryItem[{"Hierarchical: ", info["Hierarchical"]}],
          BoxForm`SummaryItem[{"Algorithm: ", info["Algorithm"]}]
        },
        form
      ]
    ]

If[$VersionNumber >= 12.0,
  IGClusterData /: Information[c : IGClusterData[asc_?clusterAscQ]] :=
      InformationData[Prepend[igClusterDataInformation[c], "ObjectType" -> "ClusterData"]];

  IGClusterData /: Information[c : IGClusterData[asc_?clusterAscQ], All] := Information[c];

  IGClusterData /: Information[c : IGClusterData[asc_?clusterAscQ], key_] :=
      With[{info = igClusterDataInformation[c]},
        If[KeyExistsQ[info, key], info[key], c[key]]
      ];
]

(* Provide short formatting for text mode as well. *)
Format[c : IGClusterData[asc_?clusterAscQ], OutputForm] :=
    StringTemplate["IGClusterData[\"Elements\" -> <``>, \"Communities\" -> <``>]"][c["ElementCount"], Length@c["Communities"]]


igClusterData[graph_][asc_] := IGClusterData@Join[<|"Elements" -> VertexList[graph]|>, asc]


(***** Helper functions *****)

PackageScope["communitiesFromMembership"]
communitiesFromMembership::usage = "communitiesFromMembership[graph, membership] converts from membership vector to vertex partitions.";
communitiesFromMembership[graph_, membership_] := Values@GroupBy[Transpose[{VertexList[graph], Round[membership]}], Last -> First]

(* guarantees community ordering, used in LeadingEigenvector *)
communitiesFromMembershipLE[graph_, membership_] := Values@KeySort@GroupBy[Transpose[{VertexList[graph], Round[membership]}], Last -> First]

PackageScope["communitiesToMembership"]
communitiesToMembership::usage = "communitiesToMembership[elements, partitions] converts from vertex partitions to a membership vector.";
communitiesToMembership[elems_, communities_] :=
    Module[{copy = communities},
      copy[[All, All]] = Range@Length[communities] - 1;
      Lookup[AssociationThread[Flatten[communities, 1], Flatten[copy, 1]], elems]
    ]

IGraphM::invcomm = "Invalid community specification for the given vertices.";

communitiesToMembershipChecked[elems_, communities_] :=
    If[Sort[elems] === Union@@communities && Length[elems] === Length@Flatten[communities, 1],
      communitiesToMembership[elems, communities],
      Message[IGraphM::invcomm]; throw[$Failed]
    ]


(***** Modularity computations *****)

PackageExport["IGModularity"]
IGModularity::usage =
    "IGModularity[graph, {{v11, v12, \[Ellipsis]}, {v21, v22, \[Ellipsis]}, \[Ellipsis]}] gives the modularity the specified partitioning of graph's vertices into communities. Edge directions are ignored.\n" <>
    "IGModularity[graph, clusterdata] uses the partitioning specified by an IGClusterData object.";

Options[IGModularity] = { "Resolution" -> 1, DirectedEdges -> False };
SyntaxInformation[IGModularity] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};
IGModularity[graph_?igGraphQ, clusters_?igClusterDataQ, opt : OptionsPattern[]] := IGModularity[graph, clusters["Communities"], opt]
IGModularity[graph_?igGraphQ, communities : {__List}, opt : OptionsPattern[]] :=
    catch@Block[{ig = igMakeFastWeighted[graph]},
      check@ig@"modularity"[
        communitiesToMembershipChecked[VertexList[graph], communities],
        OptionValue["Resolution"],
        Replace[OptionValue[DirectedEdges], Automatic -> UndirectedGraphQ[graph]]
      ]
    ]


PackageExport["IGCompareCommunities"]
IGCompareCommunities::usage =
    "IGCompareCommunities[clusterdata1, clusterdata2] compares two community structures given as IGClusterData objects using all available methods.\n" <>
    "IGCompareCommunities[clusterdata1, clusterdata2, method] compares two community structures using method.\n" <>
    "IGCompareCommunities[clusterdata1, clusterdata2, {method1, \[Ellipsis]}] compares two community structures using each given method.\n" <>
    "IGCompareCommunities[graph, communities1, communities2] compares two partitionings of the graph vertices into communities using all available methods.\n" <>
    "IGCompareCommunities[graph, communities1, communities2, method] compares two community structures using method.\n" <>
    "IGCompareCommunities[graph, communities1, communities2, {method1, \[Ellipsis]}] compares two community structures using each given method.";

igCompareCommunitiesMethods = {
  "VariationOfInformation",
  "NormalizedMutualInformation",
  "SplitJoinDistance",
  "UnadjustedRandIndex",
  "AdjustedRandIndex"
};

igCompareCommunitiesMethodsAsc = AssociationThread[igCompareCommunitiesMethods, Range@Length[igCompareCommunitiesMethods] - 1]

igCompareCommunities[elems_, c1_, c2_, method_] :=
    Block[{res},
      res = check@igraphGlobal@"compareCommunities"[communitiesToMembershipChecked[elems, c1], communitiesToMembershipChecked[elems, c2],
        Lookup[igCompareCommunitiesMethodsAsc, method, -1]
      ];
      If[method === "SplitJoinDistance", res = Round[res]];
      <| method -> res |>
    ]

amendUsage[IGCompareCommunities, "Available methods: <*igCompareCommunitiesMethods*>."];

SyntaxInformation[IGCompareCommunities] = {"ArgumentsPattern" -> {_, _, _., _.}};

IGCompareCommunities[graph_?igGraphQ, comm1 : {__List}, comm2 : {__List}, methods : {__String} : igCompareCommunitiesMethods] :=
    catch[
      Join @@ (igCompareCommunities[VertexList[graph], comm1, comm2, #]& /@ methods)
    ]

IGCompareCommunities[graph_?igGraphQ, comm1 : {__List}, comm2 : {__List}, method_String] :=
    catch@igCompareCommunities[VertexList[graph], comm1, comm2, method]

IGCompareCommunities::diff = "The compared cluster objects must contain exactly the same elements"

IGCompareCommunities[c1_?igClusterDataQ, c2_?igClusterDataQ, methods : {__String} : igCompareCommunitiesMethods] :=
    catch[
      If[c1@"Elements" =!= c2@"Elements", Message[IGCompareCommunities::diff]; throw[$Failed]];
      Join @@ (igCompareCommunities[c1@"Elements", c1@"Communities", c2@"Communities", #]& /@ methods)
    ]

IGCompareCommunities[c1_?igClusterDataQ, c2_?igClusterDataQ, method_String] := IGCompareCommunities[c1, c2, {method}]


PackageExport["IGCommunitiesEdgeBetweenness"]
IGCommunitiesEdgeBetweenness::usage = "IGCommunitiesEdgeBetweenness[graph] finds communities using the Girvan–Newman algorithm.";

(* Note: edge ordering matters, use igMake instead of igMakeFastWeighted. *)
Options[IGCommunitiesEdgeBetweenness] = { "ClusterCount" -> Automatic };
SyntaxInformation[IGCommunitiesEdgeBetweenness] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};
IGCommunitiesEdgeBetweenness[graph_?igGraphQ, OptionsPattern[]] :=
    catch@Module[{ig = igMake[graph], result, merges, betweenness, bridges, modularity, membership, removed, clusterCount},
      clusterCount = OptionValue["ClusterCount"];
      If[clusterCount === Automatic,
        clusterCount = 0
        ,
        If[Not@Internal`PositiveMachineIntegerQ[clusterCount],
          Message[IGCommunitiesEdgeBetweenness::invopt, clusterCount, "ClusterCount", Automatic];
          clusterCount = 0;
        ]
      ];
      {result, merges, betweenness, bridges, membership, modularity} = check@ig@"communityEdgeBetweenness"[clusterCount];
      removed = Part[EdgeList[graph], igIndexVec[result]];
      igClusterData[graph]@<|
        "Communities" -> communitiesFromMembership[graph, membership],
        If[modularity === None, Unevaluated@Sequence[], "Modularity" -> modularity],
        "Merges" -> igIndexVec[merges],
        "RemovedEdges" -> removed,
        "EdgeBetweenness" -> betweenness,
        "Bridges" -> Round[bridges],
        "Algorithm" -> "EdgeBetweenness"
      |>
    ]


PackageExport["IGCommunitiesWalktrap"]
IGCommunitiesWalktrap::usage =
    "IGCommunitiesWalktrap[graph] finds communities via short random walks (of length 4 by default).\n" <>
    "IGCommunitiesWalktrap[graph, steps] finds communities via random walks of length steps.";

Options[IGCommunitiesWalktrap] = { "ClusterCount" -> Automatic };
SyntaxInformation[IGCommunitiesWalktrap] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};
IGCommunitiesWalktrap[graph_?igGraphQ, steps : _?Internal`PositiveMachineIntegerQ : 4, opt : OptionsPattern[]] :=
    catch@Module[{ig = igMakeFastWeighted[graph], merges, modularity, membership, clusterCount},
      clusterCount = OptionValue["ClusterCount"];
      If[clusterCount === Automatic,
        clusterCount = 0
        ,
        If[Not@Internal`PositiveMachineIntegerQ[clusterCount],
          Message[IGCommunitiesWalktrap::invopt, clusterCount, "ClusterCount", Automatic];
          clusterCount = 0;
        ]
      ];
      {merges, membership, modularity} = check@ig@"communityWalktrap"[steps, clusterCount];
      igClusterData[graph]@<|
        "Communities" -> communitiesFromMembership[graph, membership],
        "Modularity" -> modularity,
        "Merges" -> igIndexVec[merges],
        "Algorithm" -> "Walktrap"
      |>
    ]


PackageExport["IGCommunitiesGreedy"]
IGCommunitiesGreedy::usage = "IGCommunitiesGreedy[graph] finds communities using greedy optimization of modularity.";

SyntaxInformation[IGCommunitiesGreedy] = {"ArgumentsPattern" -> {_}};
IGCommunitiesGreedy[graph_?igGraphQ] :=
    catch@Module[{ig = igMakeFastWeighted[graph], merges, modularity, membership},
      {merges, modularity, membership} = check@ig@"communityFastGreedy"[];
      igClusterData[graph]@<|
        "Communities" -> communitiesFromMembership[graph, membership],
        "Modularity" -> modularity,
        "Merges" -> igIndexVec[merges],
        "Algorithm" -> "Greedy"
      |>
    ]


PackageExport["IGCommunitiesOptimalModularity"]
IGCommunitiesOptimalModularity::usage = "IGCommunitiesOptimalModularity[graph] finds communities by maximizing the modularity through integer programming.";

SyntaxInformation[IGCommunitiesOptimalModularity] = {"ArgumentsPattern" -> {_}};
IGCommunitiesOptimalModularity[graph_?igGraphQ] :=
    catch@Module[{ig = igMakeFastWeighted[graph], modularity, membership},
      {modularity, membership} = check@ig@"communityOptimalModularity"[];
      igClusterData[graph]@<|
        "Communities" -> communitiesFromMembership[graph, membership],
        "Modularity" -> {modularity},
        "Algorithm" -> "OptimalModularity"
      |>
    ]


PackageExport["IGCommunitiesMultilevel"]
IGCommunitiesMultilevel::usage = "IGCommunitiesMultilevel[graph] finds communities using the Louvain method.";

Options[IGCommunitiesMultilevel] = { "Resolution" -> 1 };
SyntaxInformation[IGCommunitiesMultilevel] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};
IGCommunitiesMultilevel[graph_?igGraphQ, opt : OptionsPattern[]] :=
    catch@Module[{ig = igMakeFastWeighted[graph], modularity, membership, memberships},
      {modularity, membership, memberships} = check@ig@"communityMultilevel"[N@OptionValue["Resolution"]];
      igClusterData[graph]@<|
        "Communities" -> communitiesFromMembership[graph, membership],
        "Modularity" -> modularity,
        "MultilevelCommunities" -> (communitiesFromMembership[graph, #]&) /@ memberships,
        "Algorithm" -> "Multilevel"
      |>
    ]


PackageExport["IGCommunitiesLeiden"]
IGCommunitiesLeiden::usage = "IGCommunitiesLeiden[graph] finds communities using the Leiden method.";

Options[IGCommunitiesLeiden] = {
  "Resolution" -> 1,
  "Beta" -> 0.01,
  VertexWeight -> "NormalizedStrength"
};

leidenMethods = <| "NormalizedStrength" -> 1, "Constant" -> 2, "VertexWeight" -> 3, VertexWeight -> 3 |>;

SyntaxInformation[IGCommunitiesLeiden] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};
IGCommunitiesLeiden[graph_?igGraphQ, OptionsPattern[]] :=
    catch@Module[{ig = igMakeFastWeighted[graph], membership, quality},
      {membership, quality} = check@ig@"communityLeiden"[
        N@OptionValue["Resolution"], N@OptionValue["Beta"],
        Lookup[leidenMethods, OptionValue[VertexWeight], 0],
        If[OptionValue[VertexWeight] === "VertexWeight", igVertexWeights[graph], {}]
      ];
      igClusterData[graph]@<|
        "Communities" -> communitiesFromMembership[graph, membership],
        "Quality" -> quality,
        "Algorithm" -> "Leiden"
      |>
    ]


PackageExport["IGCommunitiesLabelPropagation"]
IGCommunitiesLabelPropagation::usage = "IGCommunitiesLabelPropagation[graph] finds communities by assigning labels to each vertex and then updating them by majority voting in the neighbourhood of the vertex.";

Options[IGCommunitiesLabelPropagation] = { "Initial" -> None, "Fixed" -> None };
SyntaxInformation[IGCommunitiesLabelPropagation] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

IGCommunitiesLabelPropagation::invopt = "The option value `` is not valid.";

IGCommunitiesLabelPropagation[graph_?igGraphQ, opt : OptionsPattern[]] :=
    catch@Module[{ig = igMakeFastWeighted[graph], membership, modularity, initial, fixed, vl},
      initial = OptionValue["Initial"];
      fixed = OptionValue["Fixed"];
      vl = VertexList[graph];
      If[initial =!= None,
        If[ Not@MatchQ[initial, {__List}] || Not@SubsetQ[vl, Flatten[initial, 1]],
          Message[IGCommunitiesLabelPropagation::invopt, "\"Initial\""];
          throw[$Failed];
        ];
        initial = Prepend[initial, Complement[vl, Flatten[initial, 1]]];
        initial = communitiesToMembership[vl, initial] - 1;
        ,
        initial = {};
      ];
      If[fixed =!= None,
        If[Not@SubsetQ[vl, fixed],
          Message[IGCommunitiesLabelPropagation::invopt, "\"Fixed\""];
          throw[$Failed];
        ];
        fixed = {Complement[VertexList[graph], fixed], fixed};
        fixed = communitiesToMembership[vl, fixed];
        ,
        fixed = {};
      ];
      {membership, modularity} = check@ig@"communityLabelPropagation"[initial, fixed];
      igClusterData[graph]@<|
        "Communities" -> communitiesFromMembership[graph, membership],
        "Modularity" -> {modularity},
        "Algorithm" -> "LabelPropagation"
      |>
    ]


PackageExport["IGCommunitiesInfoMAP"]
IGCommunitiesInfoMAP::usage =
    "IGCommunitiesInfoMAP[graph] finds communities using the InfoMAP algorithm. The default number of trials is 10.\n" <>
    "IGCommunitiesInfoMAP[graph, trials]";

SyntaxInformation[IGCommunitiesInfoMAP] = {"ArgumentsPattern" -> {_, _.}};
IGCommunitiesInfoMAP[graph_?igGraphQ, trials_ : 10] :=
    catch@Module[{ig = igMakeFastWeighted[graph], membership, codeLength, vertexWeights},
      vertexWeights = If[IGVertexWeightedQ[graph], igVertexWeights[graph], {}];
      {membership, codeLength} = check@ig@"communityInfoMAP"[trials, vertexWeights];
      igClusterData[graph]@<|
        "Communities" -> communitiesFromMembership[graph, membership],
        "CodeLength" -> codeLength,
        "Algorithm" -> "InfoMAP"
      |>
    ]


PackageExport["IGCommunitiesSpinGlass"]
IGCommunitiesSpinGlass::usage = "IGCommunitiesSpinGlass[graph] finds communities using a spin glass model and simulated annealing.";

Options[IGCommunitiesSpinGlass] = {
  "UpdateRule" -> "Configuration",
  "ParallelUpdating" -> False,
  "SpinCount" -> 25,
  "StartingTemperature" -> 1.,
  "StoppingTemperature" -> 0.01,
  "CoolingFactor" -> 0.99,
  "Gamma" -> 1,
  "GammaMinus" -> 1,
  Method -> Automatic
};

SyntaxInformation[IGCommunitiesSpinGlass] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

igSpinGlassUpdateRules = {"Simple", "Configuration"};
igSpinGlassUpdateRulesAsc = AssociationThread[igSpinGlassUpdateRules, Range@Length[igSpinGlassUpdateRules] - 1];
igSpinGlassMethods = {"Original", "Negative"};
igSpinGlassMethodsAsc = AssociationThread[igSpinGlassMethods, Range@Length[igSpinGlassMethods] - 1];

amendUsage[IGCommunitiesSpinGlass,
  "Available \"UpdateRule\" option values: <*igSpinGlassUpdateRules*>. Available Method options: <*igSpinGlassMethods*>."
];

IGCommunitiesSpinGlass[graph_?igGraphQ, opt : OptionsPattern[]] :=
    catch@Module[{ig = igMakeFastWeighted[graph], modularity, membership, temp, method},
      method = OptionValue[Method];
      If[method === Automatic,
        method = If[
          IGEdgeWeightedQ[graph] && TrueQ@NonPositive@Min@igEdgeWeights[graph],
          "Negative",
          "Original"
        ]
      ];
      {membership, modularity, temp} = check@ig@"communitySpinGlass"[
        OptionValue["SpinCount"], Boole@OptionValue["ParallelUpdating"],
        OptionValue["StartingTemperature"], OptionValue["StoppingTemperature"], OptionValue["CoolingFactor"],
        Lookup[igSpinGlassUpdateRulesAsc, OptionValue["UpdateRule"], -1],
        OptionValue["Gamma"],
        Lookup[igSpinGlassMethodsAsc, method, -1],
        OptionValue["GammaMinus"]
      ];
      igClusterData[graph]@<|
        "Communities" -> communitiesFromMembership[graph, membership],
        "Modularity" -> {modularity},
        "FinalTemperature" -> temp,
        "Algorithm" -> "SpinGlass"
      |>
    ]


PackageExport["IGCommunitiesLeadingEigenvector"]
IGCommunitiesLeadingEigenvector::usage = "IGCommunitiesLeadingEigenvector[graph] finds communities based on the leading eigenvector of the modularity matrix.";

Options[IGCommunitiesLeadingEigenvector] = { "Steps" -> Automatic, "ClusterCount" -> Automatic };
SyntaxInformation[IGCommunitiesLeadingEigenvector] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

IGCommunitiesLeadingEigenvector[graph_?igGraphQ, opt : OptionsPattern[]] :=
    catch@Module[{ig = igMakeFastWeighted[graph], modularity, membership, finalMembership, merges, eval, evec, clusterCount},
      clusterCount = OptionValue["ClusterCount"];
      If[clusterCount === Automatic,
        clusterCount = 0
        ,
        If[Not@Internal`PositiveMachineIntegerQ[clusterCount],
          Message[IGCommunitiesLeadingEigenvector::invopt, clusterCount, "ClusterCount", Automatic];
          clusterCount = 0;
        ]
      ];
      {membership, finalMembership, merges, eval, evec, modularity} = check@ig@"communityLeadingEigenvector"[
        Replace[OptionValue["Steps"], Automatic :> VertexCount[graph]],
        clusterCount
      ];
      igClusterData[graph]@<|
        "Communities" -> communitiesFromMembershipLE[graph, membership],
        "FinalCommunities" -> communitiesFromMembershipLE[graph, finalMembership],
        "Merges" -> igIndexVec[merges], (* TODO handle partial dendrogram in IGClusterData methods *)
        "Modularity" -> {modularity},
        "Eigenvalues" -> eval,
        "Eigenvectors" -> evec,
        "Algorithm" -> "LeadingEigenvector"
      |>
    ]


PackageExport["IGCommunitiesFluid"]
IGCommunitiesFluid::usage = "IGCommunitiesFluid[graph, clusterCount] finds communities using the fluid communities algorithm.";

SyntaxInformation[IGCommunitiesFluid] = {"ArgumentsPattern" -> {_, _}};
IGCommunitiesFluid[graph_?igGraphQ, clusterCount_] :=
    catch@Module[{ig = igMakeFast[graph], membership, modularity},
      {membership, modularity} = check@ig@"communityFluid"[clusterCount];
      igClusterData[graph]@<|
        "Communities" -> communitiesFromMembership[graph, membership],
        "Modularity" -> {modularity},
        "Algorithm" -> "FluidCommunities"
      |>
    ]
