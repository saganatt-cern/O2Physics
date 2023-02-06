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
#include <ctime>
#include <iomanip>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace constants::math;
using namespace o2::aod::hf_cand;
using namespace o2::aod::hf_cand_2prong;
using namespace o2::analysis::hf_cuts_d0_to_pi_k;

struct HfTaskFlow {
  Service<O2DatabasePDG> pdg;

  //  configurables for CCDB
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> path{"ccdbPath", "Users/k/kgajdoso/efficiency", "base path to the ccdb object"};
  Configurable<std::string> url{"ccdbUrl", "http://ccdb-test.cern.ch", "url of the ccdb repository"};
  //Configurable<std::string> date{"ccdbDate", "20221025", "date of the ccdb file"};

  //  EFFICIENCY
  TList* listEfficiency = nullptr;
  THn* hEfficiencyTrig = nullptr;
  TH1* hEfficiencyTrigHF = nullptr;
  THn* hEfficiencyAssoc = nullptr;
  bool isEfficiencyLoaded = false;

  //  PHYSICS PLOTS
  //  configurables for processing options
  Configurable<bool> processRun2{"processRun2", "false", "Flag to run on Run 2 data"};
  Configurable<bool> processRun3{"processRun3", "true", "Flag to run on Run 3 data"};
  Configurable<bool> processMc{"processMc", "false", "Flag to run on MC"};
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
  //Configurable<int> selectionFlagHfCand{"selectionFlagHfCand", 1, "Selection Flag for HF flagged candidates"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_d0_to_pi_k::vecBinsPt}, "pT bin limits"};

  //  Collision filters
  //  FIXME: The filter is applied also on the candidates! Beware!
  Filter collisionVtxZFilter = nabs(aod::collision::posZ) < zVertexMax;
  using aodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::McCollisionLabels>>;

  //  Charged track filters
  Filter trackFilter = (nabs(aod::track::eta) < etaTrackAssocMax) &&
                       (aod::track::pt > ptTrackAssocMin) &&
                       requireGlobalTrackWoPtEtaInFilter();
  using aodTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksDCA, aod::TrackSelection, aod::McTrackLabels>>;

  //  MC collision filters
  Filter mcCollisionFilter = nabs(aod::mccollision::posZ) < zVertexMax;
  using aodMcCollisions = soa::Filtered<aod::McCollisions>;

  //  MC particle filters
  Filter mcParticlesFilter = (nabs(aod::mcparticle::eta) < etaTrackAssocMax) &&
                             (aod::mcparticle::pt > ptTrackAssocMin); //&&
                             //(aod::mcparticle::sign != 0)
  using aodMcParticles = soa::Filtered<aod::McParticles>;

  //  HF candidate filter
  //  TODO: use Partition instead of filter
  Filter candidateFilter = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar;
  using hfCandidates = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0>>;
  //  HF MCrec candidate filter
  //Filter candidateMcFilter = aod::hf_sel_candidate_d0::isRecoHfFlag >= selectionFlagHfCand; 
  //using hfMcRecCandidates = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec>>;
  //  HF MCgen candidates
  using hfMcGenCandidates = soa::Join<aod::McParticles, aod::HfCand2ProngMcGen>;

  //  configurables for containers
  ConfigurableAxis axisVertex{"axisVertex", {14, -7, 7}, "vertex axis for histograms"};
  ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {72, -PIHalf, PIHalf * 3}, "delta phi axis for histograms"};
  ConfigurableAxis axisDeltaEta{"axisDeltaEta", {48, -2.4, 2.4}, "delta eta axis for histograms"};
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
  OutputObj<CorrelationContainer> mixedTpcTpcHfH{"mixedEventTpcTpcHfH"};

  //  =========================
  //      init()
  //  =========================
  void init(o2::framework::InitContext&)
  {
    //  CCDB settings
    ccdb->setURL(url.value);
    // Enabling object caching, otherwise each call goes to the CCDB server
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    //  No later than now, will be replaced by the value of the train creation
    long now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);
    //  put the date (date of creation of the efficiency file in CCDB) in a correct timestamp format
    //std::tm cfgtm = {};
    //std::stringstream ss(date.value); ss >> std::get_time(&cfgtm, "%Y%m%d");
    //cfgtm.tm_hour = 12;
    //long timestamp = std::mktime(&cfgtm) * 1000;
    
    //  get efficiency from CCDB
    //listEfficiency = ccdb->getForTimeStamp<TList>(path.value, timestamp);

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
    registry.add("hMultiplicity", "hMultiplicity", {HistType::kTH1F, {{500, 0, 500}}});
    registry.add("hVtxZ", "hVtxZ", {HistType::kTH1F, {{400, -50, 50}}});
    registry.add("hNtracks", "hNtracks", {HistType::kTH1F, {{500, 0, 500}}});

    //  histograms for event mixing
    const int maxMixBin = axisMultiplicity->size() * 14; // 14 bins for z-vertex
    registry.add("hEventCountMixing", "bin", {HistType::kTH1F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}}});
    registry.add("hEventCountHFMixing", "bin", {HistType::kTH1F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}}});
    registry.add("hEventCountSame", "bin", {HistType::kTH1F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}}});
    registry.add("hMultiplicityMixing", "hMultiplicityMixing", {HistType::kTH1F, {{500, 0, 500}}});
    registry.add("hVtxZMixing", "hVtxZMixing", {HistType::kTH1F, {{100, -10, 10}}});
    registry.add("hNtracksMixing", "hNtracksMixing", {HistType::kTH1F, {{500, 0, 500}}});
    registry.add("hMultiplicityHFMixing", "hMultiplicityHFMixing", {HistType::kTH1F, {{500, 0, 500}}});
    registry.add("hVtxZHFMixing", "hVtxZHFMixing", {HistType::kTH1F, {{100, -10, 10}}});
    registry.add("hNtracksHFMixing", "hNtracksHFMixing", {HistType::kTH1F, {{500, 0, 500}}});

    //  TRACK HISTOGRAMS
    //  histograms for associated particles
    registry.add("hYields", "multiplicity vs pT vs eta", {HistType::kTH3F, {{200, 0, 200, "multiplicity"}, {40, 0, 20, "p_{T}"}, {100, -2, 2, "#eta"}}});
    registry.add("hEtaPhi", "multiplicity vs eta vs phi", {HistType::kTH3F, {{200, 0, 200, "multiplicity"}, {100, -2, 2, "#eta"}, {200, 0, 2 * PI, "#varphi"}}});
    registry.add("hPt", "pT", {HistType::kTH1F, {{100, 0, 10, "p_{T}"}}});
    registry.add("hEta", "eta", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("hPhi", "phi", {HistType::kTH1F, {{100, 0, 2 * PI, "#varphi"}}});
    registry.add("hVzEta", "eta vs. Vz", {HistType::kTH2F, {{100, -4, 4, "#eta"}, {20, -10, 10, "Vz"}}});

    //  histograms for particles in event mixing
    registry.add("hPtMixing", "pT", {HistType::kTH1F, {{100, 0, 10, "p_{T}"}}});
    registry.add("hEtaMixing", "eta", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("hPhiMixing", "phi", {HistType::kTH1F, {{100, 0, 2 * PI, "#varphi"}}});

    //  histograms for MFT tracks
    registry.add("hEtaPhiMFT", "multiplicity vs eta vs phi in MFT", {HistType::kTH3F, {{200, 0, 200, "multiplicity"}, {100, -2, 2, "#eta"}, {200, 0, 2 * PI, "#varphi"}}});
    registry.add("hEtaMFT", "etaMFT", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("hPhiMFT", "phiMFT", {HistType::kTH1F, {{100, 0, 2 * PI, "#varphi"}}});

    //  histograms for candidates
    auto vbins = (std::vector<double>)binsPt;

    registry.add("hPtCand", "2-prong candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0, 10., "p_{T}"}}});
    registry.add("hEtaCand", "2-prong candidates;candidate #eta;entries", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("hPhiCand", "2-prong candidates;candidate #varphi;entries", {HistType::kTH1F, {{100, 0, 2 * PI, "#varphi"}}});
    registry.add("hMass", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{500, 0., 5.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    //  histograms for candidates in event mixing
    registry.add("hPtHFMixing", "pT", {HistType::kTH1F, {{100, 0, 10, "p_{T}"}}});
    registry.add("hEtaHFMixing", "eta", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("hPhiHFMixing", "phi", {HistType::kTH1F, {{100, 0, 2 * PI, "#varphi"}}});

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

    //  EFFICIENCY histograms
    registry.add("hEffTrigCheck", "efficiency; pT", {HistType::kTProfile, {{axisPtEfficiency}}});
    registry.add("hEffTrigHFCheck", "efficiency; pT", {HistType::kTProfile, {{axisPtEfficiency}}});
    registry.add("hEffAssocCheck", "efficiency; pT", {HistType::kTProfile, {{axisPtEfficiency}}});

    //  MC histograms
    registry.add("hMcVtxZ", "hMcVtxZ", {HistType::kTH1F, {{400, -50, 50}}});
    registry.add("hMcMultiplicity", "hMcMultiplicity", {HistType::kTH1F, {{500, 0, 500}}});
    registry.add("hMcMultiplicityPrimary", "hMcMultiplicityPrimary", {HistType::kTH1F, {{500, 0, 500}}});
    registry.add("hMcPt", "pT", {HistType::kTH1F, {{100, 0, 10, "p_{T}"}}});
    registry.add("hMcEta", "eta", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("hMcPhi", "phi", {HistType::kTH1F, {{100, 0, 2 * PI, "#varphi"}}});

    registry.add("hMcHfPt", "pT", {HistType::kTH1F, {{100, 0, 10, "p_{T}"}}});
    registry.add("hMcHfEta", "eta", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("hMcHfPhi", "phi", {HistType::kTH1F, {{100, 0, 2 * PI, "#varphi"}}});
    registry.add("hMcMultiplicity_HFtest", "hMcMultiplicity_HFtest", {HistType::kTH1F, {{500, 0, 500}}});
    registry.add("hMcMultiplicityPrimary_HFtest", "hMcMultiplicityPrimary_HFtest", {HistType::kTH1F, {{500, 0, 500}}});

    registry.add("hCheckPt", "pT", {HistType::kTH1F, {{100, 0, 10, "p_{T}"}}});
    registry.add("hCheckPtMC", "pT", {HistType::kTH1F, {{100, 0, 10, "p_{T}"}}});
    registry.add("hCheckPtCorr", "pT", {HistType::kTH1F, {{100, 0, 10, "p_{T}"}}});

    //  CORRELATION CONTAINERS
    sameTpcTpcHH.setObject(new CorrelationContainer("sameEventTpcTpcHH", "sameEventTpcTpcHH", corrAxis, effAxis, {}));
    sameTpcMftHH.setObject(new CorrelationContainer("sameEventTpcMftHH", "sameEventTpcMftHH", corrAxis, effAxis, {}));
    sameTpcTpcHfH.setObject(new CorrelationContainer("sameEventTpcTpcHfH", "sameEventTpcTpcHfH", corrAxis, effAxis, userAxis));
    mixedTpcTpcHH.setObject(new CorrelationContainer("mixedEventTpcTpcHH", "mixedEventTpcTpcHH", corrAxis, effAxis, {}));
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
  void fillQA(float multiplicity, float vz, TTracks tracks)
  {
    int Ntracks = 0;
    for (auto& track1 : tracks) {
      Ntracks++;
      registry.fill(HIST("hPt"), track1.pt());
      registry.fill(HIST("hEta"), track1.eta());
      registry.fill(HIST("hPhi"), track1.phi());
      registry.fill(HIST("hVzEta"), track1.eta(), vz);
      registry.fill(HIST("hYields"), multiplicity, track1.pt(), track1.eta());
      registry.fill(HIST("hEtaPhi"), multiplicity, track1.eta(), track1.phi());
    }
    registry.fill(HIST("hNtracks"), Ntracks);
  }

  template <typename TTracks>
  void fillMixingQA(float multiplicity, float vz, TTracks tracks)
  {
    registry.fill(HIST("hMultiplicityMixing"), multiplicity);
    registry.fill(HIST("hVtxZMixing"), vz);

    int Ntracks = 0;
    for (auto& track1 : tracks) {
      Ntracks++;
      registry.fill(HIST("hPtMixing"), track1.pt());
      registry.fill(HIST("hEtaMixing"), track1.eta());
      registry.fill(HIST("hPhiMixing"), track1.phi());
    }
    registry.fill(HIST("hNtracksMixing"), Ntracks);
  }

  template <typename TTracks>
  void fillHFMixingQA(float multiplicity, float vz, TTracks tracks)
  {
    registry.fill(HIST("hMultiplicityHFMixing"), multiplicity);
    registry.fill(HIST("hVtxZHFMixing"), vz);

    int Ntracks = 0;
    for (auto& track1 : tracks) {
      Ntracks++;
      registry.fill(HIST("hPtHFMixing"), track1.pt());
      registry.fill(HIST("hEtaHFMixing"), track1.eta());
      registry.fill(HIST("hPhiHFMixing"), track1.phi());
    }
    registry.fill(HIST("hNtracksHFMixing"), Ntracks);
  }

  template <typename TTracks>
  void fillMFTQA(float multiplicity, TTracks tracks)
  {
    for (auto& track1 : tracks) {
      registry.fill(HIST("hEtaMFT"), track1.eta());
      float phi = track1.phi();
      o2::math_utils::bringTo02Pi(phi);
      registry.fill(HIST("hPhiMFT"), phi);
      registry.fill(HIST("hEtaPhiMFT"), multiplicity, track1.eta(), phi);
    }
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
      registry.fill(HIST("hEtaCand"), candidate.pt());
      registry.fill(HIST("hPhiCand"), candidate.pt());
    }
  }

  template <CorrelationContainer::CFStep step, typename TTrack>
  bool isMcParticleSelected(TTrack& track)
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
    if constexpr (step <= CorrelationContainer::kCFStepAnaTopology) {
      return track.isPhysicalPrimary();
    }
    return true;
  }

  template <typename TTrack>
  bool isMcHfParticleSelected(TTrack& track)
  {
    //  check if we are selecting HF particle that we desire
    if (!(std::abs(track.flagMcMatchGen()) == 1 << DecayType::D0ToPiK)) {
      return false;
    }
    //  check rapidity of the HF particle
    if (yCandMax >= 0. && std::abs(RecoDecay::y(array{track.px(), track.py(), track.pz()}, RecoDecay::getMassPDG(track.pdgCode()))) > yCandMax) {
      return false;
    }
    return true;
  }

  template <CorrelationContainer::CFStep step, typename TTracks>
  int fillMcParticleQA(TTracks mcparticles, bool fillHist)
  {
    int Nparticles = 0;
    for (auto& mcparticle : mcparticles) {
      if (!isMcParticleSelected<step>(mcparticle)) {
        continue;
      }
      Nparticles++;
      if (fillHist) {
        registry.fill(HIST("hMcPt"), mcparticle.pt());
        registry.fill(HIST("hMcEta"), mcparticle.eta());
        registry.fill(HIST("hMcPhi"), mcparticle.phi());
      }
    }
    if (fillHist) {
      registry.fill(HIST("hMcMultiplicityPrimary"), Nparticles);
    }
    return Nparticles;
  }

  template <typename TTracks>
  int fillMcHfCandidateQA(TTracks mcparticles)
  {
    int Nparticles = 0;
    for (auto& mcparticle : mcparticles) {
      if (!isMcHfParticleSelected(mcparticle)) {
        continue;
      }
      Nparticles++;
      registry.fill(HIST("hMcHfPt"), mcparticle.pt());
      registry.fill(HIST("hMcHfEta"), mcparticle.eta());
      registry.fill(HIST("hMcHfPhi"), mcparticle.phi());
    }
    return Nparticles;
  }

  template <CorrelationContainer::CFStep step, typename TTarget, typename TTracksTrig, typename TTracksAssoc>
  void fillCorrelations(TTarget target, TTracksTrig tracks1, TTracksAssoc tracks2, float multiplicity, float posZ)
  {
    float triggerWeight = 1;
    float associatedWeight = 1;

    //  cache efficiency, because every time we get efficiency from the multi-dimensional histogram,
    //  we have to find a corresponding bin of pT, eta, etc., and this would be done within the nested loops
    float* associatedEfficiency = nullptr;
    if constexpr (step == CorrelationContainer::kCFStepCorrected) {
      if (hEfficiencyAssoc) {
        associatedEfficiency = new float[tracks2.size()];
        int i = 0;
        for (auto& track : tracks2) {
          associatedEfficiency[i++] = getEfficiency(hEfficiencyAssoc, track.eta(), track.pt(), multiplicity, posZ);
        }
      }
    }

    for (auto& track1 : tracks1) {

      //  in case of MC-generated, do additional selection on MCparticles : charge and isPhysicalPrimary
      if constexpr (step <= CorrelationContainer::kCFStepTracked) {
        if constexpr (std::is_same_v<hfCandidates, TTracksTrig>) {
          if (!isMcHfParticleSelected(track1)) {
            continue;
          }
        } else {
          if (!isMcParticleSelected<step>(track1)) {
            continue;
          }
        }
      }

      float eta1 = track1.eta();
      float pt1 = track1.pt();
      float phi1 = track1.phi();
      o2::math_utils::bringTo02Pi(phi1);

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

      //  TODO: add getter for trigger efficiency here
      if constexpr (step == CorrelationContainer::kCFStepCorrected) {
        if (fillingHFcontainer) {
          if (hEfficiencyTrigHF) {
            float trigEff = getEfficiencyHF(hEfficiencyTrigHF, track1.pt());
            if (trigEff == 0) {
              printf("I have 0 efficiency \n");
              trigEff = 1;
            }
            triggerWeight = 1./trigEff;
            registry.fill(HIST("hEffTrigHFCheck"), track1.pt(), 1./triggerWeight); // TODO: remove this, it is just a crosscheck to see if I get efficiency back
          }
        } else {
          if (hEfficiencyTrig) {
            float trigEff = getEfficiency(hEfficiencyTrig, track1.eta(), track1.pt(), multiplicity, posZ);
            if (trigEff == 0) {
              printf("I have 0 efficiency \n");
              trigEff = 1;
            }
            triggerWeight = 1./trigEff;
            registry.fill(HIST("hEffTrigCheck"), track1.pt(), 1./triggerWeight); // TODO: remove this, it is just a crosscheck to see if I get efficiency back
          }
        }
      }

      //  TODO: temporary checks. Careful! It is filled twice if I run both LF and HF case
      if constexpr (step == CorrelationContainer::kCFStepVertex) {  
        registry.fill(HIST("hCheckPtMC"), track1.pt());
      }
      if constexpr (step == CorrelationContainer::kCFStepReconstructed) {  
        registry.fill(HIST("hCheckPt"), track1.pt());
      }
      if constexpr (step == CorrelationContainer::kCFStepCorrected) {  
        registry.fill(HIST("hCheckPtCorr"), track1.pt(), triggerWeight);
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

        //  in case of MC-generated, do additional selection on MCparticles : charge and isPhysicalPrimary
        if constexpr (step <= CorrelationContainer::kCFStepTracked) {
          if (!isMcParticleSelected<step>(track2)) {
            continue;
          }
        }

        float eta2 = track2.eta();
        float pt2 = track2.pt();
        float phi2 = track2.phi();
        o2::math_utils::bringTo02Pi(phi2);

        //  TODO: add pair cuts on phi*

        //  TODO: add getter for associated efficiency here
        if constexpr (step == CorrelationContainer::kCFStepCorrected) {
          if (hEfficiencyAssoc) {
            associatedWeight = 1./associatedEfficiency[track2.filteredIndex()];
            registry.fill(HIST("hEffAssocCheck"), track2.pt(), 1./associatedWeight); // TODO: remove this, it is just a crosscheck to see if I get efficiency back
          }
        }

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

      const auto multiplicity = tracks2.size(); // get multiplicity of charged hadrons, which is used for slicing in mixing
      const auto vz = collision1.posZ();

      if constexpr (std::is_same_v<hfCandidates, TTracksTrig>) {
        registry.fill(HIST("hEventCountHFMixing"), bin);
        fillHFMixingQA(multiplicity, vz, tracks1);
      } else {
        registry.fill(HIST("hEventCountMixing"), bin);
        fillMixingQA(multiplicity, vz, tracks1);
      }

      //  fill raw correlations
      corrContainer->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
      fillCorrelations<CorrelationContainer::kCFStepReconstructed>(corrContainer, tracks1, tracks2, multiplicity, collision1.posZ());
      //  TODO fill corrected correlations
    }
  }

  // =====================================
  //    get efficiency
  // =====================================
  void loadEfficiency(uint64_t timestamp, std::string prefixTrig, std::string prefixAssoc)
  {
    if (isEfficiencyLoaded) { // we don't need to load efficiency every time a process() is called
      return;
    }

    if (path.value.empty() == false) { // ccdb path to efficiency histograms

      listEfficiency = ccdb->getForTimeStamp<TList>(path.value, timestamp);

      hEfficiencyTrig = (THn*)listEfficiency->FindObject(Form("%sefficiency", prefixTrig.c_str()));
      if (hEfficiencyTrig == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", path.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram for trigger particles from %s", path.value.c_str());

      hEfficiencyTrigHF = (TH1*)listEfficiency->FindObject("HFefficiency");
      if (hEfficiencyTrigHF == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for trigger HF particles from %s", path.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram for trigger HF particles from %s", path.value.c_str());

      hEfficiencyAssoc = (THn*)listEfficiency->FindObject(Form("%sefficiency", prefixAssoc.c_str()));
      if (hEfficiencyAssoc == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for associated particles from %s", path.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram for associated particles from %s", path.value.c_str());
    }
    isEfficiencyLoaded = true;
  }

  double getEfficiency(THn* hEfficiency, float eta, float pt, float multiplicity, float posZ)
  {
    double efficiency;
    int effVariables[4];
    effVariables[0] = hEfficiency->GetAxis(0)->FindBin(eta);
    effVariables[1] = hEfficiency->GetAxis(1)->FindBin(pt);
    effVariables[2] = hEfficiency->GetAxis(2)->FindBin(multiplicity);
    effVariables[3] = hEfficiency->GetAxis(3)->FindBin(posZ);
    efficiency = hEfficiency->GetBinContent(effVariables);
    return efficiency;
  }

  double getEfficiencyHF(TH1* hEfficiency, float pt)
  {
    double efficiency;
    efficiency = hEfficiency->GetBinContent(hEfficiency->FindBin(pt));
    return efficiency;
  }

  // =====================================
  //    process same event correlations: h-h case
  // =====================================
  void processSameTpcTpcHH(aodCollisions::iterator const& collision,
                            aod::BCsWithTimestamps const&,
                           aodTracks const& tracks)
  {
    const auto multiplicity = tracks.size();
    sameTpcTpcHH->fillEvent(multiplicity, CorrelationContainer::kCFStepBiasStudy);
    fillCorrelations<CorrelationContainer::kCFStepBiasStudy>(sameTpcTpcHH, tracks, tracks, multiplicity, collision.posZ());

    if (!(isCollisionSelected(collision, true))) {
      return;
    }

    //  the event histograms below are only filled for h-h case
    //  because there is a possibility of double-filling if more correlation
    //  options are ran at the same time
    //  temporary solution, since other correlation options always have to be ran with h-h, too
    //  TODO: rewrite it in a more intelligent way
    //const auto multiplicity = tracks.size();
    registry.fill(HIST("hMultiplicity"), multiplicity);
    registry.fill(HIST("hVtxZ"), collision.posZ());

    BinningPolicyBase<2> baseBinning{{axisVertex, axisMultiplicity}, true};
    int bin = baseBinning.getBin(std::make_tuple(collision.posZ(), multiplicity));
    registry.fill(HIST("hEventCountSame"), bin);

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    loadEfficiency(bc.timestamp(), "", "");

    fillQA(multiplicity, collision.posZ(), tracks);
    //  fill raw correlations
    sameTpcTpcHH->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
    fillCorrelations<CorrelationContainer::kCFStepReconstructed>(sameTpcTpcHH, tracks, tracks, multiplicity, collision.posZ());
    //  fill corrected correlations
    if (hEfficiencyTrig || hEfficiencyAssoc) {
      sameTpcTpcHH->fillEvent(multiplicity, CorrelationContainer::kCFStepCorrected);
      fillCorrelations<CorrelationContainer::kCFStepCorrected>(sameTpcTpcHH, tracks, tracks, multiplicity, collision.posZ());
    }
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcTpcHH, "Process same-event correlations for h-h case", true);

  // =====================================
  //    process same event correlations: HF-h case
  // =====================================
  void processSameTpcTpcHfH(aodCollisions::iterator const& collision,
                            aod::BCsWithTimestamps const&,
                            aodTracks const& tracks,
                            hfCandidates const& candidates)
  {
    if (!(isCollisionSelected(collision, false))) {
      return;
    }

    const auto multiplicity = tracks.size();

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    loadEfficiency(bc.timestamp(), "", "");

    fillCandidateQA(candidates);
    //  fill raw correlations
    sameTpcTpcHfH->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
    fillCorrelations<CorrelationContainer::kCFStepReconstructed>(sameTpcTpcHfH, candidates, tracks, multiplicity, collision.posZ());
    //  fill corrected correlations
    if (hEfficiencyTrig || hEfficiencyAssoc) {
      sameTpcTpcHfH->fillEvent(multiplicity, CorrelationContainer::kCFStepCorrected);
      fillCorrelations<CorrelationContainer::kCFStepCorrected>(sameTpcTpcHfH, candidates, tracks, multiplicity, collision.posZ());
    }
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcTpcHfH, "Process same-event correlations for HF-h case", true);

  // =====================================
  //    process same event correlations: h-MFT case
  // =====================================
  void processSameTpcMftHH(aodCollisions::iterator const& collision,
                           aodTracks const& tracks,
                           aod::MFTTracks const& mfttracks)
  {
    if (!(isCollisionSelected(collision, false))) {
      return;
    }

    const auto multiplicity = tracks.size();

    sameTpcMftHH->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
    fillMFTQA(multiplicity, mfttracks);
    fillCorrelations<CorrelationContainer::kCFStepReconstructed>(sameTpcMftHH, tracks, mfttracks, multiplicity, collision.posZ());
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcMftHH, "Process same-event correlations for h-MFT case", true);

  //  TODO: add also MFT option
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
  //    process mixed event correlations: HF-h case
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
  //    process same event MC correlations: h-h case
  // =====================================
  void processMcSameTpcTpcHH(aodMcCollisions::iterator const& mcCollision,
                             aodMcParticles const& mcParticles,
                             aodCollisions const& collisions)
  {
    const auto multiplicity = mcParticles.size(); // Note: these are all MC particles after selection (not only primary)
    registry.fill(HIST("hMcVtxZ"), mcCollision.posZ());
    registry.fill(HIST("hMcMultiplicity"), multiplicity);

    //  fill correlations for all MC collisions
    const auto multPrimaryCharge0 = fillMcParticleQA<CorrelationContainer::kCFStepVertex>(mcParticles, false);
    sameTpcTpcHH->fillEvent(multPrimaryCharge0, CorrelationContainer::kCFStepAll);
    fillCorrelations<CorrelationContainer::kCFStepAll>(sameTpcTpcHH, mcParticles, mcParticles, multPrimaryCharge0, mcCollision.posZ());

    if (collisions.size() == 0) {
      return;
    }

    //  fill correlations for MC collisions that have a reconstructed collision
    const auto multPrimaryCharge0 = fillMcParticleQA<CorrelationContainer::kCFStepVertex>(mcParticles, true);
    sameTpcTpcHH->fillEvent(multPrimaryCharge0, CorrelationContainer::kCFStepVertex);
    fillCorrelations<CorrelationContainer::kCFStepVertex>(sameTpcTpcHH, mcParticles, mcParticles, multPrimaryCharge0, mcCollision.posZ());
  }
  PROCESS_SWITCH(HfTaskFlow, processMcSameTpcTpcHH, "Process MC-gen level same-event correlations for h-h case", true);

  // =====================================
  //    process same event MC correlations: HF-h case
  // =====================================
  void processMcSameTpcTpcHfH(aodMcCollisions::iterator const& mcCollision,
                             aodMcParticles const& mcParticles,
                             hfMcGenCandidates const& mcCandidates,
                             aodCollisions const& collisions)
  {
    const auto multiplicity = mcParticles.size(); // Note: these are all MC particles after selection (not only primary)
    registry.fill(HIST("hMcMultiplicity_HFtest"), multiplicity);

    //  fill correlations for all MC collisions
    sameTpcTpcHfH->fillEvent(multiplicity, CorrelationContainer::kCFStepAll);
    fillCorrelations<CorrelationContainer::kCFStepAll>(sameTpcTpcHfH, mcCandidates, mcParticles, multiplicity, mcCollision.posZ());

    if (collisions.size() == 0) {
      return;
    }

    //  fill correlations for MC collisions that have a reconstructed collision
    const auto multPrimaryCharge0 = fillMcParticleQA<CorrelationContainer::kCFStepVertex>(mcParticles, false);
    registry.fill(HIST("hMcMultiplicityPrimary_HFtest"), multPrimaryCharge0);
    fillMcHfCandidateQA(mcCandidates);
    sameTpcTpcHfH->fillEvent(multPrimaryCharge0, CorrelationContainer::kCFStepVertex);
    fillCorrelations<CorrelationContainer::kCFStepVertex>(sameTpcTpcHfH, mcCandidates, mcParticles, multPrimaryCharge0, mcCollision.posZ());
  }
  PROCESS_SWITCH(HfTaskFlow, processMcSameTpcTpcHfH, "Process MC-gen level same-event correlations for HF-h case", true);

  // =====================================
  //    process efficiency
  // =====================================
  int GetSpecies(int pdgCode)
  {
    switch (pdgCode) {
      case 211: // pion
      case -211:
        return 0;
      case 321: // Kaon
      case -321:
        return 1;
      case 2212: // proton
      case -2212:
        return 2;
      default:
        return 3;
    }
  }

  Preslice<aod::Tracks> perCollision = aod::track::collisionId;
  void processMCEfficiency(aodMcCollisions::iterator const& mcCollision, aodMcParticles const& mcParticles, aodCollisions const& collisions, aodTracks const& tracks)
  {
    auto multiplicity = 0;
    for (auto& mcParticle : mcParticles) {
      if(isMcParticleSelected<CorrelationContainer::kCFStepAll>(mcParticle)) {
        multiplicity++;
      }
    }
    //auto multiplicity = mcCollision.multiplicity();
    for (auto& mcParticle : mcParticles) {
      if (isMcParticleSelected<CorrelationContainer::kCFStepAll>(mcParticle)) {
        sameTpcTpcHH->getTrackHistEfficiency()->Fill(CorrelationContainer::MC, mcParticle.eta(), mcParticle.pt(), GetSpecies(mcParticle.pdgCode()), multiplicity, mcCollision.posZ());
      }
    }

    for (auto& collision : collisions) {
      auto groupedTracks = tracks.sliceBy(perCollision, collision.globalIndex());

      for (auto& track : groupedTracks) {
        if (track.has_mcParticle()) {
          const auto& mcParticle = track.mcParticle();
          if (mcParticle.isPhysicalPrimary()) {
            sameTpcTpcHH->getTrackHistEfficiency()->Fill(CorrelationContainer::RecoPrimaries, mcParticle.eta(), mcParticle.pt(), GetSpecies(mcParticle.pdgCode()), multiplicity, mcCollision.posZ());
          }
          sameTpcTpcHH->getTrackHistEfficiency()->Fill(CorrelationContainer::RecoAll, mcParticle.eta(), mcParticle.pt(), GetSpecies(mcParticle.pdgCode()), multiplicity, mcCollision.posZ());
          // LOGF(info, "Filled track %d", track.globalIndex());
        } else {
          // fake track
          sameTpcTpcHH->getTrackHistEfficiency()->Fill(CorrelationContainer::Fake, track.eta(), track.pt(), 0, multiplicity, mcCollision.posZ());
        }
      }
    }
  }
  PROCESS_SWITCH(HfTaskFlow, processMCEfficiency, "MC: Extract efficiencies", false);

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfTaskFlow>(cfgc)};
}
