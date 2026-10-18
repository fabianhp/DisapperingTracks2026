//////////////////////////////////////////////////////////////
// trackletProb.h
// Probabilidad analitica de que una particula cargada de vida
// larga produzca un tracklet, dado un mapa de eficiencia en
// (|eta|, radio transversal de decaimiento).
//
// Formato del archivo de mapa (una linea por bin, '#' = comentario):
//
//     etaLo  etaHi  rLo[mm]  rHi[mm]  eficiencia
//
// La eficiencia se toma constante dentro de cada bin. Los bins
// no deben solaparse. Fuera de los bins la eficiencia es cero.
//
// Este archivo NO trae ningun mapa: hay que llenarlo con los
// valores de ATLAS (HEPData, o Fig. 2(a) de arXiv:2603.08315
// digitalizada para tracklets de cuatro capas).
//////////////////////////////////////////////////////////////
#ifndef TRACKLETPROB_H
#define TRACKLETPROB_H

#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

struct EffBin {
  double etaLo, etaHi, rLo, rHi, eff;
};

class EffMap {
public:
  std::vector<EffBin> bins;

  // Lee el mapa. Devuelve false si no se pudo abrir o quedo vacio.
  bool load(const std::string& fileName) {
    bins.clear();
    std::ifstream in(fileName.c_str());
    if (!in.good()) return false;
    std::string line;
    while (std::getline(in, line)) {
      std::size_t hash = line.find('#');
      if (hash != std::string::npos) line = line.substr(0, hash);
      std::istringstream ss(line);
      EffBin b;
      if (!(ss >> b.etaLo >> b.etaHi >> b.rLo >> b.rHi >> b.eff)) continue;
      if (b.rHi <= b.rLo || b.etaHi <= b.etaLo) continue;
      if (b.eff < 0.0 || b.eff > 1.0) continue;
      bins.push_back(b);
    }
    return !bins.empty();
  }

  // Eficiencia para un decaimiento en el radio transversal r [mm].
  // Se usa solo para el control cruzado con vertices sorteados.
  double effAt(double r, double absEta) const {
    for (std::size_t k = 0; k < bins.size(); ++k) {
      const EffBin& b = bins[k];
      if (absEta >= b.etaLo && absEta < b.etaHi && r >= b.rLo && r < b.rHi)
        return b.eff;
    }
    return 0.0;
  }

  // Probabilidad de tracklet para una particula con (pT, m, |eta|)
  // y longitud de decaimiento propia ctau [mm]:
  //
  //   lambda_T = (pT / m) * ctau
  //   P = sum_bins eff * [ exp(-rLo/lambda_T) - exp(-rHi/lambda_T) ]
  //
  // Es la integral exacta de la exponencial de decaimiento sobre
  // el mapa, sin sortear el vertice.
  double prob(double pT, double m, double absEta, double ctau) const {
    if (pT <= 0.0 || m <= 0.0 || ctau <= 0.0) return 0.0;
    const double lambdaT = pT / m * ctau;
    double P = 0.0;
    for (std::size_t k = 0; k < bins.size(); ++k) {
      const EffBin& b = bins[k];
      if (absEta < b.etaLo || absEta >= b.etaHi) continue;
      P += b.eff * (std::exp(-b.rLo / lambdaT) - std::exp(-b.rHi / lambdaT));
    }
    if (P < 0.0) P = 0.0;
    if (P > 1.0) P = 1.0;
    return P;
  }
};

// Malla logaritmica de ctau [mm].
inline std::vector<double> logGrid(double lo, double hi, int n) {
  std::vector<double> v;
  if (n < 2 || lo <= 0.0 || hi <= lo) return v;
  for (int i = 0; i < n; ++i)
    v.push_back(lo * std::pow(hi / lo, double(i) / double(n - 1)));
  return v;
}

#endif
