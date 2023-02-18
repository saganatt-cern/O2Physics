// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file taskFlow.cxx
/// \author Katarina Krizkova Gajdosova <katarina.gajdosova@cern.ch>, CERN
/// \author Maja Kabus <maja.kabus@cern.ch>, CERN

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "CCDB/BasicCCDBManager.h"
#include "Framework/GroupSlicer.h"
#include "Framework/StepTHn.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "CommonConstants/MathConstants.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Multiplicity.h"
#include "PWGCF/Core/CorrelationContainer.h"
#include "PWGCF/Core/PairCuts.h"
#include "DataFormatsParameters/GRPObject.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"

#include <TH1F.h>
#include <cmath>
#include <TDirectory.h>
#include <THn.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace constants::math;
using namespace o2::aod::hf_cand;
using namespace o2::aod::hf_cand_2prong;
using namespace o2::analysis::hf_cuts_d0_to_pi_k;

struct HfTaskFlow {
  Service<O2DatabasePDG> pdg;

  //  configurables for processing options
  Configurable<bool> processRun2{"processRun2", false, "Flag to run on Run 2 data"};
  Configurable<bool> processRun3{"processRun3", true, "Flag to run on Run 3 data"};
  Configurable<bool> processMc{"processMc", false, "Flag to run on MC"};
  Configurable<int> nMixedEvents{"nMixedEvents", 5, "Number of mixed events per event"};
  //  configurables for collisions
  Configurable<float> zVertexMax{"zVertexMax", 7.0f, "Accepted z-vertex range"};
  //  configurables for associated particles
  Configurable<float> etaTrackAssocMax{"etaTrackAssocMax", 0.8f, "max. eta of associated tracks"};
  Configurable<float> ptTrackAssocMin{"ptTrackAssocMin", 0.5f, "min. pT of associated tracks"};
  //  configurables for HF candidates
  Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Selection Flag for D0bar"};
  Configurable<double> yCandMax{"yCandMax", -1., "max. cand. rapidity"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_d0_to_pi_k::vecBinsPt}, "pT bin limits"};

  //  Collision filters
  //  FIXME: The filter is applied also on the candidates! Beware!
  Filter collisionVtxZFilter = nabs(aod::collision::posZ) < zVertexMax;
  using aodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults>>;
  using aodLabelCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::Mults>>;

  //  Charged track filters
  Filter trackFilter = (nabs(aod::track::eta) < etaTrackAssocMax) &&
                       (aod::track::pt > ptTrackAssocMin) &&
                       requireGlobalTrackWoPtEtaInFilter();
  using aodTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksDCA, aod::TrackSelection>>;

  //  MC collision filters
  Filter mcCollisionFilter = nabs(aod::mccollision::posZ) < zVertexMax;
  using aodMcCollisions = soa::Filtered<aod::McCollisions>;

  //  HF candidate filter
  //  TODO: use Partition instead of filter
  Filter candidateFilter = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar;
  using hfCandidates = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0>>;

  //  configurables for containers
  ConfigurableAxis axisVertex{"axisVertex", {14, -7, 7}, "vertex axis for histograms"};
  ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {72, -PIHalf, PIHalf * 3}, "delta phi axis for histograms"};
  ConfigurableAxis axisDeltaEta{"axisDeltaEta", {36, -1.8, 1.8}, "delta eta axis for histograms"};
  ConfigurableAxis axisDeltaEtaMFT{"axisDeltaEtaMFT", {40, 1.0, 5.0}, "delta eta axis for histograms with MFT"};
  ConfigurableAxis axisPtTrigger{"axisPtTrigger", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0}, "pt trigger axis for histograms"};
  ConfigurableAxis axisPtAssoc{"axisPtAssoc", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0}, "pt associated axis for histograms"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 100.1}, "multiplicity axis for histograms"};
  ConfigurableAxis axisVertexEfficiency{"axisVertexEfficiency", {10, -10, 10}, "vertex axis for efficiency histograms"};
  ConfigurableAxis axisEtaEfficiency{"axisEtaEfficiency", {20, -1.0, 1.0}, "eta axis for efficiency histograms"};
  ConfigurableAxis axisPtEfficiency{"axisPtEfficiency", {VARIABLE_WIDTH, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0}, "pt axis for efficiency histograms"};
  //  TODO: flow of HF will need to be done vs. invariant mass, in the signal and side-band regions
  //        either 1) add invariant mass axis or 2) define several containers for different inv. mass regions
  //        Note: don't forget to check inv. mass separately for D0 and D0bar candidate
  ConfigurableAxis axisMass{"axisMass", {30, 1.7, 2.0}, "axis of invariant mass of HF candidates"};

  HistogramRegistry registry{"registry"};

  OutputObj<CorrelationContainer> sameTpcTpcHH{"sameEventTpcTpcHH"};
  OutputObj<CorrelationContainer> sameTpcMftHH{"sameEventTpcMftHH"};
  OutputObj<CorrelationContainer> sameTpcTpcHfH{"sameEventTpcTpcHfH"};
  OutputObj<CorrelationContainer> mixedTpcTpcHH{"mixedEventTpcTpcHH"};
  OutputObj<CorrelationContainer> mixedTpcMftHH{"mixedEventTpcMftHH"};
  OutputObj<CorrelationContainer> mixedTpcTpcHfH{"mixedEventTpcTpcHfH"};

  //  =========================
  //      init()
  //  =========================
  void init(o2::framework::InitContext&)
  {
    //  EVENT HISTOGRAMS
    registry.add("hEventCounter", "hEventCounter", {HistType::kTH1F, {{3, 0.5, 3.5}}});
    //  set axes of the event counter histogram
    const int nBins = 3;
    std::string labels[nBins];
    labels[0] = "all";
    labels[1] = "after trigger selection (Run 2)";
    labels[2] = "after Physics selection";
    for (int iBin = 0; iBin < nBins; iBin++) {
      registry.get<TH1>(HIST("hEventCounter"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
    }
    registry.add("sameTpcTpcHH/hMultiplicity", "hMultiplicity", {HistType::kTH1F, {{500, 0, 500}}});
    registry.add("sameTpcTpcHH/hVtxZ", "hVtxZ", {HistType::kTH1F, {{400, -50, 50}}});
    registry.add("sameTpcTpcHH/hNtracks", "hNtracks", {HistType::kTH1F, {{500, 0, 500}}});
    registry.add("sameTpcTpcHfH/hMultiplicity", "hMultiplicity", {HistType::kTH1F, {{500, 0, 500}}});
    registry.add("sameTpcTpcHfH/hVtxZ", "hVtxZ", {HistType::kTH1F, {{400, -50, 50}}});
    registry.add("sameTpcMftHH/hMultiplicity", "hMultiplicity", {HistType::kTH1F, {{500, 0, 500}}});
    registry.add("sameTpcMftHH/hMCMultiplicity", "hMultiplicity", {HistType::kTH1F, {{500, 0, 500}}});
    registry.add("sameTpcMftHH/hVtxZ", "hVtxZ", {HistType::kTH1F, {{400, -50, 50}}});
    registry.add("sameTpcMftHH/hNtracks", "hNtracks", {HistType::kTH1F, {{500, 0, 500}}});

    //  histograms for event mixing
    const int maxMixBin = axisMultiplicity->size() * 14; // 14 bins for z-vertex
    registry.add("hEventCountMixing", "bin", {HistType::kTH1F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}}});
    registry.add("hEventCountMFTMixing", "bin", {HistType::kTH1F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}}});
    registry.add("hEventCountHFMixing", "bin", {HistType::kTH1F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}}});
    registry.add("hEventCountSame", "bin", {HistType::kTH1F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}}});
    registry.add("mixedTpcTpcHH/hMultiplicity", "hMultiplicity", {HistType::kTH1F, {{500, 0, 500}}});
    registry.add("mixedTpcTpcHH/hVtxZ", "hVtxZ", {HistType::kTH1F, {{100, -10, 10}}});
    registry.add("mixedTpcTpcHH/hNtracks", "hNtracks", {HistType::kTH1F, {{500, 0, 500}}});
    registry.add("mixedTpcTpcHfH/hMultiplicity", "hMultiplicity", {HistType::kTH1F, {{500, 0, 500}}});
    registry.add("mixedTpcTpcHfH/hVtxZ", "hVtxZ", {HistType::kTH1F, {{100, -10, 10}}});
    registry.add("mixedTpcTpcHfH/hNtracks", "hNtracks", {HistType::kTH1F, {{500, 0, 500}}});
    registry.add("mixedTpcMftHH/hMultiplicity", "hMultiplicity", {HistType::kTH1F, {{500, 0, 500}}});
    registry.add("mixedTpcMftHH/hVtxZ", "hVtxZ", {HistType::kTH1F, {{100, -10, 10}}});
    registry.add("mixedTpcMftHH/hNtracks", "hNtracks", {HistType::kTH1F, {{500, 0, 500}}});

    //  TRACK HISTOGRAMS
    //  histograms for associated particles
    registry.add("sameTpcTpcHH/hYields", "multiplicity vs pT vs eta", {HistType::kTH3F, {{200, 0, 200, "multiplicity"}, {40, 0, 20, "p_{T}"}, {100, -2, 2, "#eta"}}});
    registry.add("sameTpcTpcHH/hEtaPhi", "multiplicity vs eta vs phi", {HistType::kTH3F, {{200, 0, 200, "multiplicity"}, {100, -2, 2, "#eta"}, {200, 0, 2 * PI, "#varphi"}}});
    registry.add("sameTpcTpcHH/hPt", "pT", {HistType::kTH1F, {{100, 0, 10, "p_{T}"}}});
    registry.add("sameTpcTpcHH/hEta", "eta", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("sameTpcTpcHH/hPhi", "phi", {HistType::kTH1F, {{100, 0, 2 * PI, "#varphi"}}});

    //  histograms for particles in event mixing
    registry.add("mixedTpcTpcHH/hPt", "pT", {HistType::kTH1F, {{100, 0, 10, "p_{T}"}}});
    registry.add("mixedTpcTpcHH/hEta", "eta", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("mixedTpcTpcHH/hPhi", "phi", {HistType::kTH1F, {{100, 0, 2 * PI, "#varphi"}}});

    //  histograms for MFT tracks
    registry.add("sameTpcMftHH/hEtaPhi", "multiplicity vs eta vs phi in MFT", {HistType::kTH3F, {{200, 0, 200, "multiplicity"}, {100, -2, 2, "#eta"}, {200, 0, 2 * PI, "#varphi"}}});
    registry.add("sameTpcMftHH/hEta", "eta", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("sameTpcMftHH/hPhi", "phi", {HistType::kTH1F, {{100, 0, 2 * PI, "#varphi"}}});
    registry.add("sameTpcMftHH/hPt", "pt", {HistType::kTH1F, {{100, 0, 10, "p_{T}"}}});

    //  histograms for particles in event mixing
    registry.add("mixedTpcMftHH/hPt", "pT", {HistType::kTH1F, {{100, 0, 10, "p_{T}"}}});
    registry.add("mixedTpcMftHH/hEta", "eta", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("mixedTpcMftHH/hPhi", "phi", {HistType::kTH1F, {{100, 0, 2 * PI, "#varphi"}}});

    //  histograms for candidates
    auto vbins = (std::vector<double>)binsPt;

    registry.add("hPtCand", "2-prong candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0, 10.}}});
    registry.add("hPtProng0", "2-prong candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0, 10.}}});
    registry.add("hPtProng1", "2-prong candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0, 10.}}});
    registry.add("hMass", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{500, 0., 5.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLength", "2-prong candidates;decay length (cm);entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthXY", "2-prong candidates;decay length xy (cm);entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0", "2-prong candidates;prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1", "2-prong candidates;prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0d0", "2-prong candidates;product of DCAxy to prim. vertex (cm^{2});entries", {HistType::kTH2F, {{500, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCTS", "2-prong candidates;cos #it{#theta}* (D^{0});entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCt", "2-prong candidates;proper lifetime (D^{0}) * #it{c} (cm);entries", {HistType::kTH2F, {{120, -20., 100.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPA", "2-prong candidates;cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEtaCand", "2-prong candidates;candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hSelectionStatus", "2-prong candidates;selection status;entries", {HistType::kTH2F, {{5, -0.5, 4.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hImpParErr", "2-prong candidates;impact parameter error (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLenErr", "2-prong candidates;decay length error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLenXYErr", "2-prong candidates;decay length xy error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    //  histograms for candidates in event mixing
    registry.add("mixedTpcTpcHfH/hPt", "pT", {HistType::kTH1F, {{100, 0, 10, "p_{T}"}}});
    registry.add("mixedTpcTpcHfH/hEta", "eta", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("mixedTpcTpcHfH/hPhi", "phi", {HistType::kTH1F, {{100, 0, 2 * PI, "#varphi"}}});

    //  MONTE CARLO
    registry.add("sameMcTpcMftHH/hMultiplicity", "hMultiplicity", {HistType::kTH1F, {{500, 0, 500}}});
    registry.add("sameMcTpcMftHH/hVtxZ", "hVtxZ", {HistType::kTH1F, {{400, -50, 50}}});
    registry.add("sameMcTpcMftHH/hEta", "eta", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("sameMcTpcMftHH/hPhi", "phi", {HistType::kTH1F, {{100, 0, 2 * PI, "#varphi"}}});
    registry.add("sameMcTpcMftHH/hPt", "pt", {HistType::kTH1F, {{100, 0, 10, "p_{T}"}}});
    registry.add("sameMcTpcMftHH/hMultiplicityMft", "hMultiplicity", {HistType::kTH1F, {{500, 0, 500}}});
    registry.add("sameMcTpcMftHH/hEtaMft", "eta", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("sameMcTpcMftHH/hPhiMft", "phi", {HistType::kTH1F, {{100, 0, 2 * PI, "#varphi"}}});
    registry.add("sameMcTpcMftHH/hPtMft", "pt", {HistType::kTH1F, {{100, 0, 10, "p_{T}"}}});
    //  event mixing
    registry.add("hEventCountMcMixing", "bin", {HistType::kTH1F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}}});
    registry.add("mixedMcTpcMftHH/hMultiplicity", "hMultiplicity", {HistType::kTH1F, {{500, 0, 500}}});
    registry.add("mixedMcTpcMftHH/hVtxZ", "hVtxZ", {HistType::kTH1F, {{100, -10, 10}}});
    registry.add("mixedMcTpcMftHH/hEta", "eta", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("mixedMcTpcMftHH/hPhi", "phi", {HistType::kTH1F, {{100, 0, 2 * PI, "#varphi"}}});
    registry.add("mixedMcTpcMftHH/hPt", "pt", {HistType::kTH1F, {{100, 0, 10, "p_{T}"}}});
    registry.add("mixedMcTpcMftHH/hEtaMft", "eta", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("mixedMcTpcMftHH/hPhiMft", "phi", {HistType::kTH1F, {{100, 0, 2 * PI, "#varphi"}}});
    registry.add("mixedMcTpcMftHH/hPtMft", "pt", {HistType::kTH1F, {{100, 0, 10, "p_{T}"}}});
    registry.add("mixedMcTpcMftHH/hMultiplicityAllMcParticles", "hMultiplicity", {HistType::kTH1F, {{1000, 0, 1000}}});

    //  set axes of the correlation container
    std::vector<AxisSpec> corrAxis = {{axisDeltaEta, "#Delta#eta"},
                                      {axisPtAssoc, "p_{T} (GeV/c)"},
                                      {axisPtTrigger, "p_{T} (GeV/c)"},
                                      {axisMultiplicity, "multiplicity"},
                                      {axisDeltaPhi, "#Delta#varphi (rad)"},
                                      {axisVertex, "z-vtx (cm)"}};
    std::vector<AxisSpec> effAxis = {{axisEtaEfficiency, "#eta"},
                                     {axisPtEfficiency, "p_{T} (GeV/c)"},
                                     {axisVertexEfficiency, "z-vtx (cm)"}};
    std::vector<AxisSpec> userAxis = {{axisMass, "m_{inv} (GeV/c^{2})"}};
    std::vector<AxisSpec> corrAxisMFT = {{axisDeltaEtaMFT, "#Delta#eta"},
                                        {axisPtAssoc, "p_{T} (GeV/c)"},
                                        {axisPtTrigger, "p_{T} (GeV/c)"},
                                        {axisMultiplicity, "multiplicity"},
                                        {axisDeltaPhi, "#Delta#varphi (rad)"},
                                        {axisVertex, "z-vtx (cm)"}};

    sameTpcTpcHH.setObject(new CorrelationContainer("sameEventTpcTpcHH", "sameEventTpcTpcHH", corrAxis, effAxis, {}));
    sameTpcMftHH.setObject(new CorrelationContainer("sameEventTpcMftHH", "sameEventTpcMftHH", corrAxisMFT, effAxis, {}));
    sameTpcTpcHfH.setObject(new CorrelationContainer("sameEventTpcTpcHfH", "sameEventTpcTpcHfH", corrAxis, effAxis, userAxis));
    mixedTpcTpcHH.setObject(new CorrelationContainer("mixedEventTpcTpcHH", "mixedEventTpcTpcHH", corrAxis, effAxis, {}));
    mixedTpcMftHH.setObject(new CorrelationContainer("mixedEventTpcMftHH", "mixedEventTpcMftHH", corrAxisMFT, effAxis, {}));
    mixedTpcTpcHfH.setObject(new CorrelationContainer("mixedEventTpcTpcHfH", "mixedEventTpcTpcHfH", corrAxis, effAxis, userAxis));
  }

  //  ---------------
  //    templates
  //  FIXME: Some collisions are rejected here, what causes (part of) differences with the D0 task
  //  ---------------
  template <typename TCollision>
  bool isCollisionSelected(TCollision collision, bool fillHistograms = false)
  {
    if (processRun2 == true) {
      //  Run 2: trigger selection for data case
      if (fillHistograms)
        registry.fill(HIST("hEventCounter"), 1);
      if (!processMc) {
        if (!collision.alias()[kINT7]) {
          return false;
        }
      }
      //  Run 2: further offline selection
      if (fillHistograms)
        registry.fill(HIST("hEventCounter"), 2);
      if (!collision.sel7()) {
        return false;
      }
      if (fillHistograms)
        registry.fill(HIST("hEventCounter"), 3);
    } else {
      //  Run 3: selection
      if (fillHistograms)
        registry.fill(HIST("hEventCounter"), 1);
      if (!collision.sel8()) {
        return false;
      }
      if (fillHistograms)
        registry.fill(HIST("hEventCounter"), 3);
    }
    return true;
  }

  template <typename TTracks>
  void fillQA(float multiplicity, TTracks tracks)
  {
    int Ntracks = 0;
    for (auto& track1 : tracks) {
      Ntracks++;
      registry.fill(HIST("sameTpcTpcHH/hPt"), track1.pt());
      registry.fill(HIST("sameTpcTpcHH/hEta"), track1.eta());
      registry.fill(HIST("sameTpcTpcHH/hPhi"), track1.phi());
      registry.fill(HIST("sameTpcTpcHH/hYields"), multiplicity, track1.pt(), track1.eta());
      registry.fill(HIST("sameTpcTpcHH/hEtaPhi"), multiplicity, track1.eta(), track1.phi());
    }
    registry.fill(HIST("sameTpcTpcHH/hNtracks"), Ntracks);
  }

  template <typename TTracks>
  void fillMixingQA(float multiplicity, float vz, TTracks tracks)
  {
    registry.fill(HIST("mixedTpcTpcHH/hMultiplicity"), multiplicity);
    registry.fill(HIST("mixedTpcTpcHH/hVtxZ"), vz);

    int Ntracks = 0;
    for (auto& track1 : tracks) {
      Ntracks++;
      registry.fill(HIST("mixedTpcTpcHH/hPt"), track1.pt());
      registry.fill(HIST("mixedTpcTpcHH/hEta"), track1.eta());
      registry.fill(HIST("mixedTpcTpcHH/hPhi"), track1.phi());
    }
    registry.fill(HIST("mixedTpcTpcHH/hNtracks"), Ntracks);
  }

  template <typename TTracks>
  void fillHFMixingQA(float multiplicity, float vz, TTracks tracks)
  {
    registry.fill(HIST("mixedTpcTpcHfH/hMultiplicity"), multiplicity);
    registry.fill(HIST("mixedTpcTpcHfH/hVtxZ"), vz);

    int Ntracks = 0;
    for (auto& track1 : tracks) {
      Ntracks++;
      registry.fill(HIST("mixedTpcTpcHfH/hPt"), track1.pt());
      registry.fill(HIST("mixedTpcTpcHfH/hEta"), track1.eta());
      registry.fill(HIST("mixedTpcTpcHfH/hPhi"), track1.phi());
    }
    registry.fill(HIST("mixedTpcTpcHfH/hNtracks"), Ntracks);
  }

  template <typename TTracks>
  void fillMFTQA(float multiplicity, TTracks tracks)
  {
    int Ntracks = 0;
    for (auto& track1 : tracks) {
      Ntracks++;
      registry.fill(HIST("sameTpcMftHH/hEta"), track1.eta());
      float phi = track1.phi();
      o2::math_utils::bringTo02Pi(phi);
      registry.fill(HIST("sameTpcMftHH/hPhi"), phi);
      registry.fill(HIST("sameTpcMftHH/hPt"), track1.pt());
      registry.fill(HIST("sameTpcMftHH/hEtaPhi"), multiplicity, track1.eta(), phi);
    }
    registry.fill(HIST("sameTpcMftHH/hNtracks"), Ntracks);
  }

  template <typename TTracks>
  void fillMFTMixingQA(float multiplicity, float vz, TTracks tracks)
  {
    registry.fill(HIST("mixedTpcMftHH/hMultiplicity"), multiplicity);
    registry.fill(HIST("mixedTpcMftHH/hVtxZ"), vz);

    int Ntracks = 0;
    for (auto& track1 : tracks) {
      Ntracks++;
      float phi = track1.phi();
      o2::math_utils::bringTo02Pi(phi);
      registry.fill(HIST("mixedTpcMftHH/hPt"), track1.pt());
      registry.fill(HIST("mixedTpcMftHH/hEta"), track1.eta());
      registry.fill(HIST("mixedTpcMftHH/hPhi"), phi);
    }
    registry.fill(HIST("mixedTpcMftHH/hNtracks"), Ntracks);
  }

  //  TODO: Check how to put this into a Filter
  template <typename TTrack>
  bool isAcceptedCandidate(TTrack candidate)
  {
    if (!(candidate.hfflag() & 1 << DecayType::D0ToPiK)) {
      return false;
    }
    if (yCandMax >= 0. && std::abs(yD0(candidate)) > yCandMax) {
      return false;
    }
    return true;
  }

  //  TODO: Note: we do not need all these plots since they are in D0 and Lc task -> remove it after we are sure this works
  template <typename TTracks>
  void fillCandidateQA(TTracks candidates)
  {
    for (auto& candidate : candidates) {
      if (!isAcceptedCandidate(candidate)) {
        continue;
      }

      if (candidate.isSelD0() >= selectionFlagD0) {
        registry.fill(HIST("hMass"), invMassD0ToPiK(candidate), candidate.pt());
      }
      if (candidate.isSelD0bar() >= selectionFlagD0bar) {
        registry.fill(HIST("hMass"), invMassD0barToKPi(candidate), candidate.pt());
      }

      registry.fill(HIST("hPtCand"), candidate.pt());
      registry.fill(HIST("hPtProng0"), candidate.ptProng0());
      registry.fill(HIST("hPtProng1"), candidate.ptProng1());
      registry.fill(HIST("hDecLength"), candidate.decayLength(), candidate.pt());
      registry.fill(HIST("hDecLengthXY"), candidate.decayLengthXY(), candidate.pt());
      registry.fill(HIST("hd0Prong0"), candidate.impactParameter0(), candidate.pt());
      registry.fill(HIST("hd0Prong1"), candidate.impactParameter1(), candidate.pt());
      registry.fill(HIST("hd0d0"), candidate.impactParameterProduct(), candidate.pt());
      registry.fill(HIST("hCTS"), cosThetaStarD0(candidate), candidate.pt());
      registry.fill(HIST("hCt"), ctD0(candidate), candidate.pt());
      registry.fill(HIST("hCPA"), candidate.cpa(), candidate.pt());
      registry.fill(HIST("hEtaCand"), candidate.eta(), candidate.pt());
      registry.fill(HIST("hSelectionStatus"), candidate.isSelD0() + (candidate.isSelD0bar() * 2), candidate.pt());
      registry.fill(HIST("hImpParErr"), candidate.errorImpactParameter0(), candidate.pt());
      registry.fill(HIST("hImpParErr"), candidate.errorImpactParameter1(), candidate.pt());
      registry.fill(HIST("hDecLenErr"), candidate.errorDecayLength(), candidate.pt());
      registry.fill(HIST("hDecLenXYErr"), candidate.errorDecayLengthXY(), candidate.pt());
    }
  }

  template <typename TTrack>
  bool isMcParticleSelected(TTrack& track, bool tpcacceptance)
  {
    //  remove MC particles with charge = 0
    int8_t sign = 0;
    TParticlePDG* pdgparticle = pdg->GetParticle(track.pdgCode());
    if (pdgparticle != nullptr) {
      sign = (pdgparticle->Charge() > 0) ? 1.0 : ((pdgparticle->Charge() < 0) ? -1.0 : 0.0);
    }
    if (sign == 0) {
      return false;
    }

    //  MC particle has to be primary
    if (!track.isPhysicalPrimary()) {
      return false;
    }

    //  MC particle pt
    if (track.pt() < ptTrackAssocMin) {
      return false;
    }

    //  MC particle eta, either within TPC or MFT acceptance
    if (tpcacceptance) {
      if (track.eta() < -etaTrackAssocMax || track.eta() > etaTrackAssocMax) {
        return false;
      }
    } else {
      if (track.eta() < -4.0 || track.eta() > -2.0) {
        return false;
      }
    }

    return true;
  }

  template <typename TTracks>
  int fillMcParticleQA(TTracks mcparticles, bool fillHist)
  {
    int Nparticles = 0;
    for (auto& mcparticle : mcparticles) {
      if (!isMcParticleSelected(mcparticle, true)) {
        continue;
      }
      Nparticles++;
      if (fillHist) {
        registry.fill(HIST("sameMcTpcMftHH/hPt"), mcparticle.pt());
        registry.fill(HIST("sameMcTpcMftHH/hEta"), mcparticle.eta());
        registry.fill(HIST("sameMcTpcMftHH/hPhi"), mcparticle.phi());
      }
    }
    if (fillHist) {
      registry.fill(HIST("sameMcTpcMftHH/hMultiplicity"), Nparticles);
    }
    return Nparticles;
  }

  template <typename TTracks>
  void fillMcParticleMftQA(TTracks mcparticles, bool fillHist)
  {
    int Nparticles = 0;
    for (auto& mcparticle : mcparticles) {
      if (!isMcParticleSelected(mcparticle, false)) {
        continue;
      }
      Nparticles++;
      if (fillHist) {
        registry.fill(HIST("sameMcTpcMftHH/hPtMft"), mcparticle.pt());
        registry.fill(HIST("sameMcTpcMftHH/hEtaMft"), mcparticle.eta());
        registry.fill(HIST("sameMcTpcMftHH/hPhiMft"), mcparticle.phi());
      }
    }
    if (fillHist) {
      registry.fill(HIST("sameMcTpcMftHH/hMultiplicityMft"), Nparticles);
    }
  }

  template <typename TTracks>
  void fillMcMixingQA(float multiplicity, float vz, TTracks mcparticles)
  {
    registry.fill(HIST("mixedMcTpcMftHH/hMultiplicity"), multiplicity);
    registry.fill(HIST("mixedMcTpcMftHH/hVtxZ"), vz);

    for (auto& mcparticle : mcparticles) {
      if (!isMcParticleSelected(mcparticle, true)) {
        continue;
      }
      registry.fill(HIST("mixedMcTpcMftHH/hPt"), mcparticle.pt());
      registry.fill(HIST("mixedMcTpcMftHH/hEta"), mcparticle.eta());
      registry.fill(HIST("mixedMcTpcMftHH/hPhi"), mcparticle.phi());
    }
  }

  template <typename TTracks>
  void fillMcMixingMftQA(TTracks mcparticles)
  {
    for (auto& mcparticle : mcparticles) {
      if (!isMcParticleSelected(mcparticle, false)) {
        continue;
      }
      registry.fill(HIST("mixedMcTpcMftHH/hPtMft"), mcparticle.pt());
      registry.fill(HIST("mixedMcTpcMftHH/hEtaMft"), mcparticle.eta());
      registry.fill(HIST("mixedMcTpcMftHH/hPhiMft"), mcparticle.phi());
    }
  }

  template <CorrelationContainer::CFStep step, typename TTarget, typename TTracksTrig, typename TTracksAssoc>
  void fillCorrelations(TTarget target, TTracksTrig tracks1, TTracksAssoc tracks2, float multiplicity, float posZ)
  {
    auto triggerWeight = 1;
    auto associatedWeight = 1;

    for (auto& track1 : tracks1) {

      //  additional particle selection on MCparticles
      if constexpr (step <= CorrelationContainer::kCFStepTracked) {
        if (!isMcParticleSelected(track1, true)) {//  true to have particles within TPC acceptance
          continue;
        }
      }

      float eta1 = track1.eta();
      float pt1 = track1.pt();
      float phi1 = track1.phi();
      o2::math_utils::bringTo02Pi(phi1);

      //  TODO: add getter for NUE trigger efficiency here

      //  calculating inv. mass to be filled into the container below
      //  Note: this is needed only in case of HF-hadron correlations
      bool fillingHFcontainer = false;
      double invmass = 0;
      if constexpr (std::is_same_v<hfCandidates, TTracksTrig>) {
        //  TODO: Check how to put this into a Filter
        if (!isAcceptedCandidate(track1)) {
          continue;
        }
        fillingHFcontainer = true;
        invmass = invMassD0ToPiK(track1);
      }

      //  fill single-track distributions
      if (!fillingHFcontainer) {
        target->getTriggerHist()->Fill(step, pt1, multiplicity, posZ, triggerWeight);
      } else {
        target->getTriggerHist()->Fill(step, pt1, multiplicity, posZ, invmass, triggerWeight);
      }

      for (auto& track2 : tracks2) {

        //  case of h-h correlations where the two types of tracks are the same
        //  this avoids autocorrelations and double counting of particle pairs
        if constexpr (std::is_same_v<TTracksAssoc, TTracksTrig>) {
          if (track1.index() <= track2.index()) {
            continue;
          }
        }

        //  in case of HF-h correlations, remove candidate daughters from the pool of associated hadrons
        //  with which the candidate is being correlated
        if constexpr (std::is_same_v<hfCandidates, TTracksTrig>) {
          if ((track1.prong0Id() == track2.globalIndex()) || (track1.prong1Id() == track2.globalIndex())) {
            continue;
          }
        }

        //  additional particle selection on MCparticles
        if constexpr (step <= CorrelationContainer::kCFStepTracked) {
          if (!isMcParticleSelected(track2, false)) {// false for particles within MFT acceptance
            continue;
          }
        }

        float eta2 = track2.eta();
        float pt2 = track2.pt();
        float phi2 = track2.phi();
        o2::math_utils::bringTo02Pi(phi2);

        //  TODO: add getter for NUE associated efficiency here

        //  TODO: add pair cuts on phi*

        float deltaPhi = phi1 - phi2;
        //  set range of delta phi in (-pi/2 , 3/2*pi)
        deltaPhi = RecoDecay::constrainAngle(deltaPhi, -0.5 * PI);

        if (!fillingHFcontainer) {
          //  fill pair correlations
          target->getPairHist()->Fill(step,
                                      eta1 - eta2, pt2, pt1, multiplicity, deltaPhi, posZ,
                                      triggerWeight * associatedWeight);
        } else {
          target->getPairHist()->Fill(step,
                                      eta1 - eta2, pt2, pt1, multiplicity, deltaPhi, posZ, invmass,
                                      triggerWeight * associatedWeight);
        }
      }
    }
  }

  template <typename TTracksTrig, typename TTracksAssoc, typename TLambda>
  void mixCollisions(aodCollisions& collisions, TTracksTrig& tracks1, TTracksAssoc& tracks2, TLambda getPartsSize, OutputObj<CorrelationContainer>& corrContainer)
  {
    using BinningType = FlexibleBinningPolicy<std::tuple<decltype(getPartsSize)>, aod::collision::PosZ, decltype(getPartsSize)>;
    BinningType binningWithTracksSize{{getPartsSize}, {axisVertex, axisMultiplicity}, true};

    auto tracksTuple = std::make_tuple(tracks1, tracks2);
    Pair<aodCollisions, TTracksTrig, TTracksAssoc, BinningType> pair{binningWithTracksSize, nMixedEvents, -1, collisions, tracksTuple};

    for (auto& [collision1, tracks1, collision2, tracks2] : pair) {

      if (!(isCollisionSelected(collision1, false))) {
        continue;
      }
      if (!(isCollisionSelected(collision2, false))) {
        continue;
      }

      auto binningValues = binningWithTracksSize.getBinningValues(collision1, collisions);
      int bin = binningWithTracksSize.getBin(binningValues);

      // get multiplicity of charged hadrons, which is used for slicing in mixing
      const auto multiplicity = getPartsSize(collision1);
      const auto vz = collision1.posZ();

      if constexpr (std::is_same_v<hfCandidates, TTracksTrig>) {
        registry.fill(HIST("hEventCountHFMixing"), bin);
        fillHFMixingQA(multiplicity, vz, tracks1);
      } 
      if constexpr (std::is_same_v<aodTracks, TTracksTrig>) {
        registry.fill(HIST("hEventCountMixing"), bin);
        fillMixingQA(multiplicity, vz, tracks1);
      }
      if constexpr (std::is_same_v<mftTracks, TTracksAssoc>) {
        registry.fill(HIST("hEventCountMFTMixing"), bin);
        fillMFTMixingQA(multiplicity, vz, tracks2);
      }

      corrContainer->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
      fillCorrelations<CorrelationContainer::kCFStepReconstructed>(corrContainer, tracks1, tracks2, multiplicity, collision1.posZ());
    }
  }

  template <typename TTracksTrig, typename TTracksAssoc, typename TLambda>
  void mixMcCollisions(aodMcCollisions& mcCollisions, TTracksTrig& tracks1, TTracksAssoc& tracks2, TLambda getPartsSize, OutputObj<CorrelationContainer>& corrContainer, aodCollisions& collisions)
  {
    using BinningType = FlexibleBinningPolicy<std::tuple<decltype(getPartsSize)>, aod::McCollision::PosZ, decltype(getPartsSize)>;
    BinningType binningWithTracksSize{{getPartsSize}, {axisVertex, axisMultiplicity}, true};

    auto tracksTuple = std::make_tuple(tracks1, tracks2);
    Pair<aodMcCollisions, TTracksTrig, TTracksAssoc, BinningType> pair{binningWithTracksSize, nMixedEvents, -1, mcCollisions, tracksTuple};

    for (auto& [collision1, tracks1, collision2, tracks2] : pair) {

      //  TODO: remove mcCollisions which don't have any reconstructed collision
      registry.fill(HIST("mixedMcTpcMftHH/hMultiplicityAllMcParticles"), getPartsSize(collision1));

      auto binningValues = binningWithTracksSize.getBinningValues(collision1, mcCollisions);
      int bin = binningWithTracksSize.getBin(binningValues);

      // get multiplicity of charged hadrons, which is used for slicing in mixing
      const auto multiplicity = fillMcParticleQA(tracks1, false);
      const auto vz = collision1.posZ();

      registry.fill(HIST("hEventCountMcMixing"), bin);
      fillMcMixingQA(multiplicity, vz, tracks1);
      fillMcMixingMftQA(tracks1);

      corrContainer->fillEvent(multiplicity, CorrelationContainer::kCFStepVertex);
      fillCorrelations<CorrelationContainer::kCFStepVertex>(corrContainer, tracks1, tracks2, multiplicity, collision1.posZ());
    }
  }

  // =====================================
  //    process same event correlations: h-h case
  // =====================================
  void processSameTpcTpcHH(aodCollisions::iterator const& collision,
                           aodTracks const& tracks)
  {
    if (!(isCollisionSelected(collision, true))) {
      return;
    }

    const auto multiplicity = tracks.size();
    registry.fill(HIST("sameTpcTpcHH/hMultiplicity"), multiplicity);
    registry.fill(HIST("sameTpcTpcHH/hVtxZ"), collision.posZ());

    BinningPolicyBase<2> baseBinning{{axisVertex, axisMultiplicity}, true};
    int bin = baseBinning.getBin(std::make_tuple(collision.posZ(), multiplicity));
    registry.fill(HIST("hEventCountSame"), bin);

    sameTpcTpcHH->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
    fillQA(multiplicity, tracks);
    fillCorrelations<CorrelationContainer::kCFStepReconstructed>(sameTpcTpcHH, tracks, tracks, multiplicity, collision.posZ());
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcTpcHH, "Process same-event correlations for h-h case", true);

  // =====================================
  //    process same event correlations: HF-h case
  // =====================================
  void processSameTpcTpcHfH(aodCollisions::iterator const& collision,
                            aodTracks const& tracks,
                            hfCandidates const& candidates)
  {
    if (!(isCollisionSelected(collision, false))) {
      return;
    }

    const auto multiplicity = tracks.size();
    registry.fill(HIST("sameTpcTpcHfH/hMultiplicity"), multiplicity);
    registry.fill(HIST("sameTpcTpcHfH/hVtxZ"), collision.posZ());

    sameTpcTpcHfH->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
    fillCandidateQA(candidates);
    fillCorrelations<CorrelationContainer::kCFStepReconstructed>(sameTpcTpcHfH, candidates, tracks, multiplicity, collision.posZ());
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcTpcHfH, "Process same-event correlations for HF-h case", true);

  // =====================================
  //    process same event correlations: h-MFT case
  // =====================================
  Filter mfttrackFilter = (aod::fwdtrack::eta > -3.9f) && (aod::fwdtrack::eta < -2.0f);
  using mftTracks = soa::Filtered<aod::MFTTracks>;

  void processSameTpcMftHH(aodCollisions::iterator const& collision,
                           aodTracks const& tracks,
                           mftTracks const& mfttracks)
  {
    if (!(isCollisionSelected(collision, false))) {
      return;
    }

    const auto multiplicity = tracks.size();
    registry.fill(HIST("sameTpcMftHH/hMultiplicity"), multiplicity);
    registry.fill(HIST("sameTpcMftHH/hVtxZ"), collision.posZ());

    sameTpcMftHH->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
    fillMFTQA(multiplicity, mfttracks);
    fillCorrelations<CorrelationContainer::kCFStepReconstructed>(sameTpcMftHH, tracks, mfttracks, multiplicity, collision.posZ());
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcMftHH, "Process same-event correlations for h-MFT case", true);

  // =====================================
  //    process mixed event correlations: h-h case
  // =====================================
  void processMixedTpcTpcHH(aodCollisions& collisions,
                            aodTracks& tracks)
  {
    //  we want to group collisions based on charged-track multiplicity
    auto getTracksSize = [&tracks](aodCollisions::iterator const& col) {
      auto associatedTracks = tracks.sliceByCached(o2::aod::track::collisionId, col.globalIndex()); // it's cached, so slicing/grouping happens only once
      auto size = associatedTracks.size();
      return size;
    };

    mixCollisions(collisions, tracks, tracks, getTracksSize, mixedTpcTpcHH);
  }
  PROCESS_SWITCH(HfTaskFlow, processMixedTpcTpcHH, "Process mixed-event correlations for h-h case", true);

  // =====================================
  //    process mixed event correlations: h-h case
  // =====================================
  void processMixedTpcTpcHfH(aodCollisions& collisions,
                             aodTracks& tracks,
                             hfCandidates& candidates)
  {
    //  we want to group collisions based on charged-track multiplicity
    auto getTracksSize = [&tracks](aodCollisions::iterator const& col) {
      auto associatedTracks = tracks.sliceByCached(o2::aod::track::collisionId, col.globalIndex());
      auto size = associatedTracks.size();
      return size;
    };

    mixCollisions(collisions, candidates, tracks, getTracksSize, mixedTpcTpcHfH);
  }
  PROCESS_SWITCH(HfTaskFlow, processMixedTpcTpcHfH, "Process mixed-event correlations for HF-h case", true);

  // =====================================
  //    process mixed event correlations: h-MFT case
  // =====================================
  void processMixedTpcMftHH(aodCollisions& collisions,
                            aodTracks& tracks,
                            mftTracks& mfttracks)
  {
    //  we want to group collisions based on charged-track multiplicity
    auto getTracksSize = [&tracks](aodCollisions::iterator const& col) {
      auto associatedTracks = tracks.sliceByCached(o2::aod::track::collisionId, col.globalIndex()); // it's cached, so slicing/grouping happens only once
      auto size = associatedTracks.size();
      return size;
    };

    mixCollisions(collisions, tracks, mfttracks, getTracksSize, mixedTpcMftHH);
  }
  PROCESS_SWITCH(HfTaskFlow, processMixedTpcMftHH, "Process mixed-event correlations for h-MFT case", true);
/*
  // =====================================
  //    Monte Carlo reconstructed: process same event correlations: h-MFT case
  // =====================================
  Preslice<aod::McParticles> perMcCollision = aod::mcparticle::mcCollisionId;
  void processMcRecoSameTpcMftHH(aodLabelCollisions::iterator const& collision,
                                aodTracks const& tracks,
                                aod::MFTTracks const& mfttracks,
                                aod::McParticles const& mcParticles)
  {
    //  calculate MCgen multiplicity, which will be used as Nch axis in the containers
    if (collision.mcCollisionId() < 0) {
      return;
    }
    auto groupedParticles = mcParticles.sliceBy(perMcCollision, collision.mcCollisionId());
    const auto mcMultiplicity = fillMcParticleQA(groupedParticles, false);
    registry.fill(HIST("sameTpcMftHH/hMCMultiplicity"), mcMultiplicity);

    if (!(isCollisionSelected(collision, false))) {
     return;
    }
  
    const auto multiplicity = tracks.size();
    registry.fill(HIST("sameTpcMftHH/hMultiplicity"), multiplicity);
    registry.fill(HIST("sameTpcMftHH/hVtxZ"), collision.posZ());

    sameTpcMftHH->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
    fillMFTQA(multiplicity, mfttracks);
    fillCorrelations<CorrelationContainer::kCFStepReconstructed>(sameTpcMftHH, tracks, mfttracks, multiplicity, collision.posZ());
  }
  PROCESS_SWITCH(HfTaskFlow, processMcRecoSameTpcMftHH, "Process MCreconstructed with labels same-event correlations for h-MFT case", true);
*/
  // =====================================
  //    Monte Carlo: process same event correlations: h-MFT case
  // =====================================
  //  TODO: do partitions to select particles in TPC or MFT acceptance
  void processMcSameTpcMftHH(aodMcCollisions::iterator const& mcCollision,
                             aod::McParticles const& mcParticles,
                             aodCollisions const& collisions)
  {
    if (collisions.size() == 0) {
      return;
    }

    registry.fill(HIST("sameMcTpcMftHH/hVtxZ"), mcCollision.posZ());

    const auto multiplicity = fillMcParticleQA(mcParticles, true);
    fillMcParticleMftQA(mcParticles, true);
    sameTpcMftHH->fillEvent(multiplicity, CorrelationContainer::kCFStepVertex);
    fillCorrelations<CorrelationContainer::kCFStepVertex>(sameTpcMftHH, mcParticles, mcParticles, multiplicity, mcCollision.posZ());
  }
  PROCESS_SWITCH(HfTaskFlow, processMcSameTpcMftHH, "Process MC same-event correlations for h-MFT case", true);

  // =====================================
  //    Monte Carlo: process mixed event correlations: h-MFT case
  // =====================================
  void processMcMixedTpcMftHH(aodMcCollisions& mcCollisions,
                              aod::McParticles& mcParticles,
                              aodCollisions& collisions)
  {
    //  we want to group collisions based on charged-track multiplicity
    auto getTracksSize = [&mcParticles](aodMcCollisions::iterator const& col) {
      auto associatedTracks = mcParticles.sliceByCached(o2::aod::mcparticle::mcCollisionId, col.globalIndex()); // it's cached, so slicing/grouping happens only once
      int nParticles = 0;
      for (auto mcParticle : associatedTracks) {
        //  MC particle has to be primary
        if (!mcParticle.isPhysicalPrimary()) {
          continue;
        }
        //  MC particle pt
        if (mcParticle.pt() < 0.2) {
          continue;
        }
        //  MC particle eta, either within TPC or MFT acceptance
        if (mcParticle.eta() < -0.8 || mcParticle.eta() > 0.8) {
          continue;
        }
        nParticles++;
      }
      //auto size = associatedTracks.size();
      auto size = nParticles;
      return size;
    };

    mixMcCollisions(mcCollisions, mcParticles, mcParticles, getTracksSize, mixedTpcMftHH, collisions);
  }
  PROCESS_SWITCH(HfTaskFlow, processMcMixedTpcMftHH, "Process MC mixed-event correlations for h-MFT case", true);

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfTaskFlow>(cfgc)};
}
