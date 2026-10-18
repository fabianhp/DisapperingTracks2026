//////////////////////////////////////////////////////////////
// ATLAS Disappearing Track Search (13 TeV, Run 2)
// Flujo de cortes para eficiencia:
//  (1) ALL → (2) Matching (MLM) → (3) Trigger (genMET>140)
//  (4) Lepton veto → (5) Jet pT / Δφ (reco-MET)
//  (6) Tracklet reco (eff) → (7) Signal Region (pT>100 GeV)
// ---------------------------------------------------------
// Fabián Hernández Pinto, 2025
//////////////////////////////////////////////////////////////

#include "Pythia8/Pythia.h"
#include "ToyDetector-ATLAS.h"
#include "Matching.h"
#include "trackletEfficiency.h"

#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <cmath>
#include <iostream>
#include <vector>

using namespace Pythia8;
//4 INPUT PARAMETERS TO IDENTIFY PATH TO LHE.GZ FILE
int main(int argc, char *argv[]) {
  string a = argv[1]; // ENTER NEUTRAL VECTOR MASS VALUE
  string b = argv[2]; // ENTER MASS SPLIT VALUE
  string c = argv[3]; // ENTER LAMBDA_HV /HIGGS PORTAL VALUE
  string d = argv[4]; // FOLDER (RELATED WITH MASS VALUE)----> TO CORRECT
  #ifdef GZIPSUPPORT
    bool useGzip = true;
  #else
    bool useGzip = false;
  #endif
    cout << " useGzip = " << useGzip << endl;
  Pythia pythia;
  Event& event = pythia.event;
  pythia.readString("Beams:frameType = 4");
  pythia.readString("Beams:LHEF =/user/f/fhpinto/projects/MG5_aMC_v2_9_21/vquint-"+d+"/Events/run-"+a+"-"+b+"-"+c+"/unweighted_events.lhe.gz");
  pythia.init();
    //DATA FILES
  ofstream effCutFlow;
  effCutFlow.open("RPV_GRID/effCutFlow_GRID-"+a+"-"+b+"-"+c+".dat"); 
  ofstream decays;
  decays.open("RPV_GRID/decays_GRID-"+a+"-"+b+"-"+c+".dat"); 
  ofstream benchmarkParamsEff;
  benchmarkParamsEff.open("RPV_GRID/effParams_GRID-"+a+"-"+b+"-"+c+".dat"); 

  ofstream promptElec;
  promptElec.open("RPV_GRID/promptElec_GRID-"+a+"-"+b+"-"+c+".dat"); 
  ToyDetector detector;

  // -----------------------------
  // Contadores
  // -----------------------------
  int nTotalEvents     = 0;  // ALL
  int nPassedMatching  = 0;  // MLM matching
  int nPassedTrigger   = 0;  // genMET trigger prefilter
  int nPassedLeptonVeto= 0;
  int nPassedJetCuts   = 0;
  int nTrackletReco    = 0;
  int nPassedSR        = 0;

  // -----------------------------
  // Parámetros
  // -----------------------------
  const double trigger_MET_thresh = 140.0;
  const double prefilter_leadJet_thresh = 140.0;
  const double Rmatch = 0.5;
  const double pTmin  = 30.0;
  const double leadJetCut = 140.0;
  const double dphiCut = 1.0;
  const double trkletCut = 100.0;
  int nAbort = 50, iAbort = 0;

  // -----------------------------
  // Loop de eventos
  // -----------------------------
  for (int iEvent = 0; ; ++iEvent) {
    if (!pythia.next()) {
      if (pythia.info.atEndOfFile()) break;
      if (++iAbort < nAbort) continue;
      break;
    }
    ++nTotalEvents;

    // =====================================================
    // (2) Matching ME ↔ jets (MLM)
    // =====================================================
    std::vector<Vec4> MEpartons;
    for (int i = 0; i < event.size(); ++i) {
      int st = event[i].status();
      if (st == 23 || st == -23) {
        int pid = std::abs(event[i].id());
        if (pid == 21 || (pid >= 1 && pid <= 6))
          MEpartons.push_back(event[i].p());
      }
    }

    std::vector<fastjet::PseudoJet> particles;
    for (int i = 0; i < event.size(); ++i) {
      if (!event[i].isFinal()) continue;
      if (!event[i].isVisible()) continue;
      particles.emplace_back(event[i].px(), event[i].py(), event[i].pz(), event[i].e());
    }

    fastjet::ClusterSequence clustSeq(particles, fastjet::JetDefinition(fastjet::antikt_algorithm, Rmatch));
    std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(clustSeq.inclusive_jets(pTmin));

    bool passMatching = true;
    for (auto &p : MEpartons) {
      bool matched = false;
      for (auto &j : jets) {
        double deta = p.eta() - j.eta();
        double dphi = std::fabs(p.phi() - j.phi());
        if (dphi > M_PI) dphi = 2 * M_PI - dphi;
        double dR = std::sqrt(deta*deta + dphi*dphi);
        if (dR < Rmatch && j.pt() > pTmin) { matched = true; break; }
      }
      if (!matched) { passMatching = false; break; }
    }

    if (!passMatching) continue;
    ++nPassedMatching;

    // =====================================================
    // generator MET
    // =====================================================
    Vec4 genMETvec(0.,0.,0.,0.);
    for (int i = 0; i < event.size(); ++i) {
      const Particle& p = event[i];
      if (!p.isFinal()) continue;
      if (!p.isVisible()) genMETvec += p.p();
    }
    double genMET = genMETvec.pT();

    // quick leading jet estimate (gen-level)
    std::vector<fastjet::PseudoJet> quickInputs;
    for (int i = 0; i < event.size(); ++i) {
      const Particle& p = event[i];
      if (!p.isFinal() || !p.isVisible()) continue;
      if (p.pT() < 0.5) continue;
      quickInputs.emplace_back(p.px(), p.py(), p.pz(), p.e());
    }
    double quick_leadJet_pt = 0.0;
    if (!quickInputs.empty()) {
      fastjet::ClusterSequence quickSeq(quickInputs, fastjet::JetDefinition(fastjet::antikt_algorithm, 0.4));
      auto quickJets = fastjet::sorted_by_pt(quickSeq.inclusive_jets(20.0));
      if (!quickJets.empty()) quick_leadJet_pt = quickJets[0].pt();
    }

    bool passTrigger = (genMET >= trigger_MET_thresh || quick_leadJet_pt >= prefilter_leadJet_thresh);
    if (!passTrigger) continue;
    ++nPassedTrigger;

    // =====================================================
    // (4) Detector reconstruction
    // =====================================================
    if (!detector.getObjects(event)) continue;

    // Lepton veto antes de tracklets
    bool hasLepton = (!detector.electrons.empty() || !detector.muons.empty());
    if (hasLepton) continue;
    ++nPassedLeptonVeto;

    // =====================================================
    // (5) Jet cuts (pT, Δφ usando MET reconstruida)
    // =====================================================
    if (detector.jets.empty()) continue;
    double leadJetPT = detector.jets.at(0).pT();

    double dphiMin = 10.0;
    for (auto &j : detector.jets) {
      double dphi = std::fabs(j.phi() - detector.MET.phi()); // <-- MET reconstruida
      if (dphi > M_PI) dphi = 2 * M_PI - dphi;
      if (dphi < dphiMin) dphiMin = dphi;
    }

    if (leadJetPT < leadJetCut) continue;
    if (dphiMin < dphiCut) continue;
    ++nPassedJetCuts;

    // =====================================================
    // (6) Tracklet reconstruction (efficiency)
    // =====================================================
    if (!detector.getTracklets(event)) continue;
    ++nTrackletReco;

    // =====================================================
    // (7) Signal Region: pT(tracklet) > 100 GeV
    // =====================================================
    double leadTrkPT = detector.tracklets.at(0).pT();
    if (leadTrkPT < trkletCut) continue;
    ++nPassedSR;
  }

  // =====================================================
  // Resultados finales
  // =====================================================
  double sigma = pythia.info.sigmaGen(); // in mb
  double width_p;
  double width_pp;  
  if((pythia.particleData.tau0(9000006)>0)){ 
    width_pp = 3.*1e11*6.581e-25/pythia.particleData.tau0(9000006);
  }  
  cout << "####################################################################" << endl;
  cout << " Cross section is "  << sigma * 1e12 << " [fb] "<<endl;
  cout << "####################################################################" << endl;
  benchmarkParamsEff<<setprecision(9)<<" "<<pythia.particleData.tau0(9000006)<<" "<<width_p<<" "<<width_pp<<setprecision(9)<<fixed<<" "<<nPassedSR/10000<<" "<<nPassedSR<<" "<<a<<" "<<b<<" "<<c<<" "<<sigma*1e12<<endl;;
    
  std::cout << "\n=============================================\n";
  std::cout << "     ATLAS DT — Efficiency Breakdown\n";
  std::cout << "---------------------------------------------\n";
  std::cout << " (1) All events            : " << nTotalEvents << "\n";
  std::cout << " (2) Passed matching (MLM) : " << nPassedMatching
            << " (" << 100.0*nPassedMatching/nTotalEvents << "%)\n";
  std::cout << " (3) Passed trigger (genMET>140) : " << nPassedTrigger
            << " (" << 100.0*nPassedTrigger/nPassedMatching << "%)\n";
  std::cout << " (4) Lepton veto           : " << nPassedLeptonVeto
            << " (" << 100.0*nPassedLeptonVeto/nPassedTrigger << "%)\n";
  std::cout << " (5) Jet pT/Δφ (reco-MET)  : " << nPassedJetCuts
            << " (" << 100.0*nPassedJetCuts/nPassedLeptonVeto << "%)\n";
  std::cout << " (6) Tracklet reco (eff)   : " << nTrackletReco
            << " (" << 100.0*nTrackletReco/nPassedJetCuts << "%)\n";
  std::cout << " (7) SR (pT>100 GeV)       : " << nPassedSR
            << " (" << 100.0*nPassedSR/nPassedJetCuts << "%)\n";
  std::cout << "=============================================\n";
    
  pythia.stat();
  effCutFlow.close();  
  benchmarkParamsEff.close();
  decays.close();
  promptElec.close();
    //massNtrk.close();
  return 0;
}
