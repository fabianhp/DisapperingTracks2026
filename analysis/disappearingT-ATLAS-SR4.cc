//////////////////////////////////////////////////////////////
// ATLAS Disappearing Track Search, 137 fb-1 (arXiv:2603.08315)
// Recast de las regiones de cuatro capas SR4High y SR4Mid.
//
// Una sola pasada por los eventos de una masa entrega la
// eficiencia para toda una malla de ctau: la probabilidad de
// tracklet se integra analiticamente (trackletProb.h) en vez
// de usar el vertice de decaimiento sorteado.
//
// Objetos (jets, leptones aislados, MET): ToyDetector-ATLAS.
// Seleccion de evento (Tabla 1 de ATLAS):
//   - cero leptones aislados (e, mu con pT > 10 GeV, |eta| < 2.5)
//   - >= 1 jet con pT > 20 GeV y |eta| < 2.8; jet lider pT > 100 GeV
//   - |dphi(jet lider, MET)| > 1.0
//   - MET > 300 GeV (SR4High)  o  150 < MET < 300 GeV (SR4Mid)
// Tracklet (particula cargada de vida larga a nivel de verdad):
//   - pT > 60 GeV, 0.1 < |eta| < 1.9
//   - dR > 0.4 respecto de jets; aislamiento de trazas y calorimetro
//   - mapa de eficiencia en (|eta|, radio de decaimiento)
//
// Uso:
//   ./dtSR4 <eventos.lhe.gz> <mapaEficiencia.dat> <salida.dat>
//           [idLLP=1000024] [nJetMax=2] [qCut=45] [estable=0]
//
// Salida: tabla  ctau[mm]  effHigh  errHigh  effMid  errMid
// con la seccion eficaz tras el matching en el encabezado.
// N_sig = sigma * L * eff se calcula despues.
//
// Fabian Hernandez Pinto, 2026
//////////////////////////////////////////////////////////////

#include "Pythia8/Pythia.h"
// El header se llama JetMatchingHooks.h desde Pythia 8.31x.
#if __has_include("Pythia8Plugins/JetMatchingHooks.h")
#include "Pythia8Plugins/JetMatchingHooks.h"
#else
#include "Pythia8Plugins/JetMatching.h"
#endif
#include "trackletProb.h"
#include "ToyDetector-ATLAS.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

using namespace Pythia8;

// -----------------------------
// Parametros de la seleccion
// -----------------------------
// Jets, leptones y MET los construye ToyDetector (anti-kT R = 0.4,
// jets con pT > 20 GeV; leptones con pT > 10 GeV, |eta| < 2.5, aislados).
const double etaMaxVis     = 4.5;    // particulas usadas en el aislamiento
const double jetEtaMax     = 2.8;
const double leadJetCut    = 100.0;
const double dphiCut       = 1.0;
const double metHighCut    = 300.0;
const double metMidLo      = 150.0;

const double trkPtMin      = 60.0;   // pT de verdad: sin smearing
const double trkEtaMin     = 0.1;
const double trkEtaMax     = 1.9;
const double trkJetDR      = 0.4;
const double trkIsoCone    = 0.4;
const double trkIsoTrkMax  = 0.04;   // pTcone40 / pT
const double trkIsoCaloMax = 5.0;    // ET en el cono [GeV]
const double isoTrackPtMin = 1.0;

// Malla de ctau [mm]
const double ctauLo = 1.0;
const double ctauHi = 1.0e4;
const int    nCtau  = 61;

bool isNeutrino(int idAbs) { return idAbs == 12 || idAbs == 14 || idAbs == 16; }

// Las particulas nuevas que Pythia crea desde los bloques QNUMBERS del LHE
// quedan marcadas como visibles. Una particula final nueva, neutra y sin
// color (p.ej. V0 del quintuplete) no deja senal en el detector: debe
// contar en la MET y no entrar en los jets.
bool isInvisibleBSM(const Pythia8::Particle& p) {
  return p.idAbs() >= 1000000 && p.isFinal() && !p.isCharged() && p.colType() == 0;
}

// Distancia en (eta, phi) entre dos cuadrimomentos.
double dREtaPhi(const Pythia8::Vec4& a, const Pythia8::Vec4& b) {
  double dphi = std::fabs(a.phi() - b.phi());
  if (dphi > M_PI) dphi = 2.0 * M_PI - dphi;
  const double deta = a.eta() - b.eta();
  return std::sqrt(deta * deta + dphi * dphi);
}

int main(int argc, char* argv[]) {

  if (argc < 4) {
    std::cerr << "Uso: " << argv[0]
              << " <eventos.lhe.gz> <mapaEficiencia.dat> <salida.dat>"
              << " [idLLP=1000024] [nJetMax=2] [qCut=45] [estable=0]\n";
    return 1;
  }
  const std::string lheFile = argv[1];
  const std::string mapFile = argv[2];
  const std::string outFile = argv[3];
  const int    idLLP     = (argc > 4) ? std::atoi(argv[4]) : 1000024;
  const int    nJetMax   = (argc > 5) ? std::atoi(argv[5]) : 2;
  const double qCut      = (argc > 6) ? std::atof(argv[6]) : 45.0;
  const bool   stableLLP = (argc > 7) ? (std::atoi(argv[7]) != 0) : false;

  EffMap effMap;
  if (!effMap.load(mapFile)) {
    std::cerr << "ERROR: no se pudo leer el mapa de eficiencia " << mapFile
              << " (o esta vacio).\n";
    return 1;
  }
  const std::vector<double> ctaus = logGrid(ctauLo, ctauHi, nCtau);

  // -----------------------------
  // Pythia + matching MLM
  // -----------------------------
  Pythia pythia;
  Event& event = pythia.event;
  pythia.readString("Beams:frameType = 4");
  pythia.readString("Beams:LHEF = " + lheFile);

  // Los ajustes JetMatching:* solo tienen efecto si se registra el
  // hook. Sin setUserHooksPtr no se aplica ningun veto de matching.
  pythia.readString("JetMatching:merge = on");
  pythia.readString("JetMatching:scheme = 1");
  pythia.readString("JetMatching:setMad = off");
  pythia.readString("JetMatching:qCut = " + std::to_string(qCut));
  pythia.readString("JetMatching:nQmatch = 4");
  pythia.readString("JetMatching:nJetMax = " + std::to_string(nJetMax));
  pythia.readString("JetMatching:jetAlgorithm = 2");
  pythia.readString("JetMatching:slowJetPower = 1");
  pythia.readString("JetMatching:coneRadius = 1.0");
  pythia.readString("JetMatching:etaJetMax = 5.0");
  pythia.readString("JetMatching:doShowerKt = off");
  std::shared_ptr<JetMatchingMadgraph> matching
      = std::make_shared<JetMatchingMadgraph>();
  pythia.setUserHooksPtr(matching);   // Pythia 8.3; en 8.2 usar puntero crudo

  if (stableLLP)
    pythia.readString(std::to_string(idLLP) + ":mayDecay = off");

  if (!pythia.init()) {
    std::cerr << "ERROR: fallo pythia.init()\n";
    return 1;
  }

  // Las particulas nuevas que Pythia crea desde los bloques QNUMBERS del
  // LHE quedan marcadas como visibles. Una neutra, sin color y estable
  // (p.ej. V0 del quintuplete) debe ser invisible para ToyDetector: si no,
  // entraria en los jets y la MET saldria nula.
  for (int id = 1; id != 0; id = pythia.particleData.nextId(id)) {
    if (id < 1000000) continue;
    if (pythia.particleData.chargeType(id) != 0) continue;
    if (pythia.particleData.colType(id) != 0) continue;
    if (pythia.particleData.mayDecay(id)) continue;
    if (!pythia.particleData.isVisible(id)) continue;
    pythia.particleData.isVisible(id, false);
    if (id < 9000000) continue;          // solo se anuncian las del modelo
    std::cout << " AVISO: particula " << id << " ("
              << pythia.particleData.name(id)
              << ") neutra y estable: se trata como invisible.\n";
  }

  // Detector
  ToyDetector detector;
  detector.idLLP = idLLP;

  // -----------------------------
  // Acumuladores (indice 0 = SR4High, 1 = SR4Mid)
  // -----------------------------
  double sumWgen = 0.0;                         // eventos tras el matching
  double nLepVeto = 0.0, nJetCuts = 0.0;
  double nMET[2] = {0.0, 0.0};
  double nCand[2] = {0.0, 0.0};                 // con >= 1 candidato a tracklet
  std::vector<double> sumW[2], sumW2[2];
  for (int r = 0; r < 2; ++r) {
    sumW[r].assign(ctaus.size(), 0.0);
    sumW2[r].assign(ctaus.size(), 0.0);
  }
  // Control cruzado con el vertice sorteado por el generador/Pythia
  double sumWmc[2] = {0.0, 0.0}, sumW2mc[2] = {0.0, 0.0};
  double nDecayedLLP = 0.0;

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
    // Los eventos vetados por el matching no llegan aqui: next() los
    // salta internamente y sigmaGen() ya los descuenta.
    const double wEvt = pythia.info.weight();
    sumWgen += wEvt;

    // ---------------------------------------------------
    // (1) Particulas cargadas de vida larga: ultima copia
    // ---------------------------------------------------
    std::vector<int> iLLP;
    for (int i = 0; i < event.size(); ++i) {
      if (event[i].idAbs() != idLLP) continue;
      if (event[i].iBotCopyId() != i) continue;
      iLLP.push_back(i);
    }

    // ---------------------------------------------------
    // (2) Objetos: ToyDetector-ATLAS
    //     Jets anti-kT, leptones aislados y MET = -(jets + leptones).
    //     La LLP no entra ni en jets ni en la MET, igual que el
    //     tracklet en la MET de ATLAS.
    // ---------------------------------------------------
    if (!detector.getObjects(event)) continue;
    const int nLeptons = int(detector.leptons.size());
    const double met    = detector.MET.pT();
    const double metPhi = detector.MET.phi();

    std::vector<Vec4> jets;                // jets del analisis: |eta| < 2.8
    for (std::size_t j = 0; j < detector.jets.size(); ++j)
      if (std::fabs(detector.jets[j].eta()) < jetEtaMax)
        jets.push_back(detector.jets[j]);

    // Particulas visibles finales, para el aislamiento del tracklet.
    std::vector<int> iVis;
    for (int i = 0; i < event.size(); ++i) {
      const Particle& p = event[i];
      if (!p.isFinal() || !p.isVisible()) continue;
      if (p.idAbs() == idLLP) continue;
      if (isNeutrino(p.idAbs())) continue;
      if (isInvisibleBSM(p)) continue;
      if (std::fabs(p.eta()) > etaMaxVis) continue;
      iVis.push_back(i);
    }

    // ---------------------------------------------------
    // (3) Seleccion de evento
    // ---------------------------------------------------
    if (nLeptons > 0) continue;
    nLepVeto += wEvt;

    if (jets.empty()) continue;
    if (jets[0].pT() < leadJetCut) continue;
    double dphi = std::fabs(jets[0].phi() - metPhi);
    if (dphi > M_PI) dphi = 2.0 * M_PI - dphi;
    if (dphi < dphiCut) continue;
    nJetCuts += wEvt;

    int region = -1;
    if (met > metHighCut) region = 0;
    else if (met > metMidLo) region = 1;
    if (region < 0) continue;
    nMET[region] += wEvt;

    // ---------------------------------------------------
    // (4) Candidatos a tracklet
    // ---------------------------------------------------
    std::vector<int> iCand;
    for (std::size_t n = 0; n < iLLP.size(); ++n) {
      const int i = iLLP[n];
      const Particle& c = event[i];
      if (c.pT() < trkPtMin) continue;
      const double absEta = std::fabs(c.eta());
      if (absEta < trkEtaMin || absEta > trkEtaMax) continue;

      bool nearJet = false;
      for (std::size_t j = 0; j < jets.size(); ++j)
        if (dREtaPhi(jets[j], c.p()) < trkJetDR) { nearJet = true; break; }
      if (nearJet) continue;

      // Aislamiento; se excluyen los productos de decaimiento de la LLP.
      double isoTrk = 0.0, isoCalo = 0.0;
      for (std::size_t j = 0; j < iVis.size(); ++j) {
        const Particle& q = event[iVis[j]];
        if (q.mother1() == i) continue;
        if (RRapPhi(c.p(), q.p()) > trkIsoCone) continue;
        if (q.isCharged() && q.pT() > isoTrackPtMin) isoTrk += q.pT();
        if (q.idAbs() != 13) isoCalo += q.eT();
      }
      if (isoTrk / c.pT() > trkIsoTrkMax) continue;
      if (isoCalo > trkIsoCaloMax) continue;

      iCand.push_back(i);
    }
    if (iCand.empty()) continue;
    nCand[region] += wEvt;

    // ---------------------------------------------------
    // (5) Peso del evento para cada ctau:
    //     w = 1 - prod_i (1 - P_i)   (al menos un tracklet)
    // ---------------------------------------------------
    for (std::size_t k = 0; k < ctaus.size(); ++k) {
      double q = 1.0;
      for (std::size_t n = 0; n < iCand.size(); ++n) {
        const Particle& c = event[iCand[n]];
        q *= 1.0 - effMap.prob(c.pT(), c.m(), std::fabs(c.eta()), ctaus[k]);
      }
      const double w = wEvt * (1.0 - q);
      sumW[region][k]  += w;
      sumW2[region][k] += w * w;
    }

    // ---------------------------------------------------
    // (6) Control cruzado: vertice sorteado, si la LLP decayo
    //     en el registro de eventos. Vale solo para el ctau
    //     con que se genero la muestra.
    // ---------------------------------------------------
    double qMC = 1.0;
    bool anyDecayed = false;
    for (std::size_t n = 0; n < iCand.size(); ++n) {
      const Particle& c = event[iCand[n]];
      if (c.isFinal()) continue;
      anyDecayed = true;
      const double rDec = std::sqrt(c.xDec() * c.xDec() + c.yDec() * c.yDec());
      qMC *= 1.0 - effMap.effAt(rDec, std::fabs(c.eta()));
    }
    if (anyDecayed) {
      nDecayedLLP += wEvt;
      const double w = wEvt * (1.0 - qMC);
      sumWmc[region]  += w;
      sumW2mc[region] += w * w;
    }
  }

  // =====================================================
  // Resultados
  // =====================================================
  pythia.stat();
  const double sigmaFb = pythia.info.sigmaGen() * 1e12;   // mb -> fb
  const double ctau0   = pythia.particleData.tau0(idLLP); // mm
  const double mLLP    = pythia.particleData.m0(idLLP);

  if (sumWgen <= 0.0) {
    std::cerr << "ERROR: ningun evento aceptado.\n";
    return 1;
  }

  std::ofstream out(outFile.c_str());
  out << std::setprecision(8) << std::scientific;
  out << "# ATLAS DT 137 fb-1 (arXiv:2603.08315), recast SR4High / SR4Mid\n";
  out << "# idLLP                 = " << idLLP << "\n";
  out << "# masa LLP [GeV]        = " << mLLP << "\n";
  out << "# sigma matching [fb]   = " << sigmaFb << "\n";
  out << "# eventos leidos        = " << pythia.info.nTried() << "\n";
  out << "# eventos tras matching = " << sumWgen << "\n";
  out << "# tras veto de leptones = " << nLepVeto << "\n";
  out << "# tras jet lider y dphi = " << nJetCuts << "\n";
  out << "# MET > 300 (SR4High)   = " << nMET[0] << "\n";
  out << "# 150 < MET < 300 (Mid) = " << nMET[1] << "\n";
  out << "# con candidato (High)  = " << nCand[0] << "\n";
  out << "# con candidato (Mid)   = " << nCand[1] << "\n";
  if (nDecayedLLP > 0.0) {
    out << "# control vertice sorteado, ctau0 [mm] = " << ctau0 << "\n";
    for (int r = 0; r < 2; ++r)
      out << "#   eff" << (r == 0 ? "High" : "Mid ") << " = "
          << sumWmc[r] / sumWgen << " +- "
          << std::sqrt(sumW2mc[r]) / sumWgen << "\n";
  }
  out << "# ctau[mm]  effHigh  errHigh  effMid  errMid\n";
  for (std::size_t k = 0; k < ctaus.size(); ++k) {
    out << ctaus[k]
        << "  " << sumW[0][k] / sumWgen
        << "  " << std::sqrt(sumW2[0][k]) / sumWgen
        << "  " << sumW[1][k] / sumWgen
        << "  " << std::sqrt(sumW2[1][k]) / sumWgen << "\n";
  }
  out.close();

  std::cout << "\n=============================================\n";
  std::cout << "  ATLAS DT 137 fb-1: flujo de cortes\n";
  std::cout << "---------------------------------------------\n";
  std::cout << " sigma tras matching [fb] : " << sigmaFb << "\n";
  std::cout << " eventos leidos           : " << pythia.info.nTried() << "\n";
  std::cout << " tras matching            : " << sumWgen << "\n";
  std::cout << " veto de leptones         : " << nLepVeto << "\n";
  std::cout << " jet lider y dphi         : " << nJetCuts << "\n";
  std::cout << " MET > 300 (SR4High)      : " << nMET[0] << "\n";
  std::cout << " 150 < MET < 300 (SR4Mid) : " << nMET[1] << "\n";
  std::cout << " con candidato (High/Mid) : " << nCand[0] << " / " << nCand[1] << "\n";
  std::cout << " tabla escrita en         : " << outFile << "\n";
  std::cout << "=============================================\n";

  return 0;
}
