//////////////////////////////////////////////////////////////
// ToyDetector-ATLAS.h
// Emulacion simple del detector ATLAS para el analisis de
// disappearing tracks, basada en ToyDetector (N. Desai).
//
// Este header se reconstruyo a partir de ToyDetector-ATLAS.cc,
// porque el original no venia con el codigo. Declara los
// miembros que usa ese archivo y agrega tres ajustes:
//   - idLLP   : codigo PDG de la particula cargada de vida larga.
//   - effMap  : mapa de eficiencia de tracklets (trackletProb.h),
//               en lugar de trackletEfficiency.h.
//   - trkPtMin, trkEtaMin, trkEtaMax : cortes del tracklet.
//////////////////////////////////////////////////////////////
#ifndef TOYDETECTOR_ATLAS_H
#define TOYDETECTOR_ATLAS_H

#include "Pythia8/Pythia.h"
#include "trackletProb.h"

#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>

#include <vector>

namespace Pythia8 {

class ToyDetector {

public:

  ToyDetector() : idLLP(0), effMap(nullptr), trkPtMin(100.0),
                  trkEtaMin(0.0), trkEtaMax(1.9) {}

  // Construye jets, leptones aislados y MET a partir del evento.
  bool getObjects(const Event& event);

  // Tracklets con el vertice de decaimiento sorteado por Pythia y un
  // sorteo de la eficiencia. Requiere effMap. Devuelve true si hay
  // al menos uno.
  bool getTracklets(const Event& event);

  void printObjects();

  // Semilla del generador usado en el sorteo de la eficiencia.
  void seed(int s) { rndm.init(s); }

  // Objetos reconstruidos
  std::vector<Vec4> jets, electrons, muons, leptons, tracklets;
  Vec4 MET;        // -(suma de jets y leptones)
  Vec4 METnaive;   // suma de las particulas invisibles

  // Configuracion
  int           idLLP;       // particula cargada de vida larga (0 = ninguna)
  const EffMap* effMap;      // mapa de eficiencia de tracklets
  double        trkPtMin, trkEtaMin, trkEtaMax;

private:

  double REtaPhi(const Vec4& a, const Vec4& b) const;

  std::vector<fastjet::PseudoJet> fjInputs;
  Rndm rndm;

  static double Rparam;
  static fastjet::Strategy strategy;
  static fastjet::RecombinationScheme recombScheme;
  static fastjet::JetDefinition* jetDef;
};

} // namespace Pythia8

#endif
