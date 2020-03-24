(* Mathematica Package *)
(* Created by Mathematica plugin for IntelliJ IDEA *)

(* :Author: szhorvat *)
(* :Date: 2020-01-25 *)
(* :Copyright: (c) 2020 Szabolcs Horvát *)

Package["IGraphM`"]
igContextSetup[igPackagePrivateSymbol]

(*******************************)
(***** Processes on graphs *****)
(*******************************)

bundleSIR[{time_, i_, r_, s_}] :=
    TimeSeries[
      N@Transpose[{s,i,r}], {time}, (* convert to Real with N[] for better performance *)
      ValueDimensions -> 3,
      (* 0th order interpolation, uses the value at the last available data point *)
      ResamplingMethod -> {"Interpolation", "HoldFrom" -> Left}
    ]

PackageExport["IGSIRProcess"]
IGSIRProcess::usage =
    "IGSIRProcess[graph, {\[Beta], \[Gamma]}] runs a stochastic epidemic SIR model on graph with infection rate \[Beta] and recovery rate \[Gamma], and returns a time series of {S, I, R} values.\n" <>
    "IGSIRProcess[graph, {\[Beta], \[Gamma]}, n] performs n SIR model runs.";
SyntaxInformation[IGSIRProcess] = {"ArgumentsPattern" -> {_, _, _.}};
IGSIRProcess[graph_?igGraphQ, {beta_?NonNegative, gamma_?NonNegative}, n_?Internal`PositiveMachineIntegerQ] :=
    catch@Module[{ig = igMakeUnweighted[graph], result},
      result = check@ig@"sirProcess"[beta, gamma, n];
      TemporalData[
        bundleSIR /@ result,
        (* Descriptive names for the three components: *)
        MetaInformation -> {"ComponentNames" -> {"S", "I", "R"}}
      ]
    ]
IGSIRProcess[graph_?igGraphQ, {beta_?NonNegative, gamma_?NonNegative}] :=
    catch@TimeSeries[
      check@IGSIRProcess[graph, {beta, gamma}, 1],
      ResamplingMethod -> {"Interpolation", "HoldFrom" -> Left}
    ]
