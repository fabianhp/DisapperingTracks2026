#include "ToyDetector-ATLAS.h"
#include <cmath>
#include <algorithm>
#include <iostream>

namespace Pythia8 {

// ======================================================
// ============   CONFIGURACIÓN FASTJET   ===============
// ======================================================
double ToyDetector::Rparam = 0.4;
fastjet::Strategy ToyDetector::strategy = fastjet::Best;
fastjet::RecombinationScheme ToyDetector::recombScheme = fastjet::E_scheme;
fastjet::JetDefinition* ToyDetector::jetDef =
    new fastjet::JetDefinition(fastjet::antikt_algorithm,
                               ToyDetector::Rparam,
                               ToyDetector::recombScheme,
                               ToyDetector::strategy);

double ToyDetector::REtaPhi(const Vec4& a, const Vec4& b) const {
  double deta = a.eta() - b.eta();
  double dphi = fabs(a.phi() - b.phi());
  if (dphi > M_PI) dphi = 2 * M_PI - dphi;
  return sqrt(deta * deta + dphi * dphi);
}

// ======================================================
// ============   OBJETOS ESTÁNDAR (JETS, LEPTONS) ======
// ======================================================
bool ToyDetector::getObjects(const Event& event) {

  if (event.size() < 1) return false;

  jets.clear();
  electrons.clear();
  muons.clear();
  leptons.clear();
  tracklets.clear();
  MET = Vec4(0.0, 0.0, 0.0, 0.0);
  METnaive = Vec4(0.0, 0.0, 0.0, 0.0);
  fjInputs.clear();

  for (int i = 0; i < event.size(); ++i) {
    const Particle& p = event[i];
    if (!p.isFinal()) continue;

    // [cambio] La particula cargada de vida larga no entra en los jets
    // ni en la MET (solo es final si se la dejo estable).
    if (idLLP != 0 && p.idAbs() == idLLP) continue;

    // Invisibles que salen del detector
    if (!p.isVisible()) {
      double prodR = sqrt(pow2(p.xProd()) + pow2(p.yProd()));
      double prodZ = fabs(p.zProd());
      if (prodR < 1100.0 && prodZ < 2800.0)
        METnaive += p.p();
      continue;
    }

    // Fiducial region
    double prodR = sqrt(pow2(p.xProd()) + pow2(p.yProd()));
    double prodZ = fabs(p.zProd());
    if (prodR > 1100.0 || prodZ > 2800.0) continue;
    if (fabs(p.eta()) > 4.5) continue;

    // Leptones aislados
    if (std::abs(p.id()) == 11 || std::abs(p.id()) == 13) {
      if (fabs(p.eta()) > 2.5 || p.pT() < 10.0) continue;

      double ptSum = 0.0;
      for (int j = 0; j < event.size(); ++j) {
        if (i == j) continue;
        const Particle& q = event[j];
        if (!q.isFinal() || !q.isCharged()) continue;
        if (REtaPhi(p.p(), q.p()) < 0.2) ptSum += q.pT();
      }

      if (std::abs(p.id()) == 11 && ptSum / p.pT() > 0.1) continue;
      if (std::abs(p.id()) == 13 && ptSum > 1.8) continue;

      // [cambio] antes todos los leptones iban a "electrons"
      if (std::abs(p.id()) == 11) electrons.push_back(p.p());
      else                        muons.push_back(p.p());
      leptons.push_back(p.p());
      MET -= p.p();
      continue;
    }

    // Hadrones → jets
    fjInputs.emplace_back(p.px(), p.py(), p.pz(), p.e());
  }

  // Clusterización FastJet
  fastjet::ClusterSequence clustSeq(fjInputs, *jetDef);
  std::vector<fastjet::PseudoJet> sortedJets =
      fastjet::sorted_by_pt(clustSeq.inclusive_jets(20.0));

  for (auto& j : sortedJets) {
    Vec4 pj(j.px(), j.py(), j.pz(), j.e());
    jets.push_back(pj);
    MET -= pj;
  }

  return true;
}

// ======================================================
// ============   DISAPPEARING TRACKLETS   ==============
// ======================================================
bool ToyDetector::getTracklets(const Event& event) {
  tracklets.clear();

  for (int i = 0; i < event.size(); ++i) {
    const Particle& p = event[i];

    // Solo partículas cargadas
    if (p.charge() == 0) continue;
    // [cambio] solo la particula de vida larga: sin esto, hiperones de
    // alto pT que decaen en vuelo dentro de un jet contaban como tracklets
    if (idLLP != 0 && p.idAbs() != idLLP) continue;
    // [cambio] ultima copia, para no contar dos veces la misma particula
    if (p.iBotCopyId() != i) continue;
    // Debe haber decaído (tiene hijas)
    if (p.daughter1() <= 0) continue;

    // Radio de decaimiento
    double vx = p.xDec() - p.xProd();
    double vy = p.yDec() - p.yProd();
    double rDec = sqrt(vx * vx + vy * vy);

    // Región fiducial del detector interno
    // [cambio] la ventana radial la define el mapa de eficiencia
    if (fabs(p.eta()) > trkEtaMax || fabs(p.eta()) < trkEtaMin) continue;
    if (p.pT() < trkPtMin) continue;

    // Eficiencia ATLAS (mapa en radio de decaimiento y |eta|)
    if (effMap == nullptr) continue;
    double eff = effMap->effAt(rDec, fabs(p.eta()));
    if (eff <= 0.0) continue;

    double rnd = rndm.flat();
    if (rnd > eff) continue;

    tracklets.push_back(p.p());

#ifdef DEBUG_TRACKLETS
    std::cout << "[TRACKLET] id=" << p.id()
              << " pT=" << p.pT()
              << " eta=" << p.eta()
              << " rDec=" << rDec
              << " eff=" << eff << std::endl;
#endif
  }

  // =========================================================
  // Añadir MET proveniente del V0 invisible que escapa
  // =========================================================
  for (int i = 0; i < event.size(); ++i) {
    const Particle& p = event[i];
    if (p.charge() == 0) continue;     // solo V+
    if (p.daughter1() <= 0) continue;  // que haya decaído

    double vx = p.xDec() - p.xProd();
    double vy = p.yDec() - p.yProd();
    double rDec = sqrt(vx * vx + vy * vy);
    double zDec = fabs(p.zDec());
    if (rDec > 1100.0 || zDec > 2800.0) continue; // fuera del detector

    // Hijas invisibles que escapan
    for (int j = p.daughter1(); j <= p.daughter2(); ++j) {
      if (j <= 0 || j >= event.size()) continue;
      const Particle& child = event[j];
      if (!child.isVisible()) {
        double rDecChild = sqrt(pow2(child.xProd()) + pow2(child.yProd()));
        double zDecChild = fabs(child.zProd());
        if (rDecChild > 1100.0 || zDecChild > 2800.0) {
          METnaive += child.p();
#ifdef DEBUG_MET
          std::cout << "[MET_ADD] Invisible V0 id=" << child.id()
                    << " pT=" << child.pT()
                    << " eta=" << child.eta()
                    << " added to MET" << std::endl;
#endif
        }
      }
    }
  }

  return (!tracklets.empty());
}

// ======================================================
// ============   IMPRESIÓN DE OBJETOS   ================
// ======================================================
void ToyDetector::printObjects() {
  if (!jets.empty()) {
    std::cout << "Jets (pT, eta, phi):\n";
    for (auto& j : jets)
      std::cout << j.pT() << "  " << j.eta() << "  " << j.phi() << "\n";
  }
  if (!electrons.empty()) {
    std::cout << "Electrons (pT, eta, phi):\n";
    for (auto& e : electrons)
      std::cout << e.pT() << "  " << e.eta() << "  " << e.phi() << "\n";
  }
  if (!muons.empty()) {
    std::cout << "Muons (pT, eta, phi):\n";
    for (auto& m : muons)
      std::cout << m.pT() << "  " << m.eta() << "  " << m.phi() << "\n";
  }
  if (!tracklets.empty()) {
    std::cout << "Tracklets (pT, eta, phi):\n";
    for (auto& t : tracklets)
      std::cout << t.pT() << "  " << t.eta() << "  " << t.phi() << "\n";
  }
  std::cout << "MET pT: " << MET.pT()
            << "  phi: " << MET.phi() << std::endl;
  std::cout << "------------------------------------\n";
}

} // namespace Pythia8
