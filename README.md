# Disappearing tracks from electroweak vector dark matter at the LHC

Code, model files and inputs for the recast of the ATLAS full Run-2
disappearing-track search used in

> P. Escalona, F. Hernández-Pinto and A. R. Zerwekh,
> *Disappearing Tracks from Electroweak Vector Dark Matter at the LHC*,
> [arXiv:2610.14182](https://arxiv.org/abs/2610.14182) [hep-ph].

The signal is a spin-1 dark matter candidate in the real five-dimensional
representation of SU(2)<sub>L</sub> with Y = 0 (vector quintuplet). The singly
charged state V<sup>±</sup> decays through V<sup>±</sup> → V<sup>0</sup> π<sup>±</sup>
after a macroscopic distance and is reinterpreted with the four-layer
signal region SR4High of the ATLAS search
(JHEP 07 (2026) 152, [arXiv:2603.08315](https://arxiv.org/abs/2603.08315)).
The tracklet reconstruction probability is integrated analytically over the
decay radius using the efficiency published by ATLAS, so one event sample per
mass covers the full lifetime range.

If you use this code, please cite the paper above (see `CITATION.cff`).

## Contents

```
model/
  VQUINT_UFO/            UFO model of the vector quintuplet (MadGraph5_aMC@NLO)
  feynrules/             FeynRules source (vquint.fr), corrected Lagrangian and
                         Mathematica notebook used to generate the UFO
generation/
  wino/                  MadGraph grid of pure winos (MSSM_SLHA2), for validation
  quintuplet/            MadGraph grid of the vector quintuplet
analysis/
  disappearingT-ATLAS-SR4.cc   event selection of SR4High / SR4Mid (Pythia 8)
  ToyDetector-ATLAS.{h,cc}     simplified ATLAS detector response
  trackletProb.h               analytic tracklet probability vs. ctau
  compilar.sh                  build script  ->  ./dtSR4
  run_analysis.sh              runs dtSR4 over every mass of a grid
  limites.py                   N_sig and mass limits; wino validation plot
  limites_quint.py             quintuplet limits in (M_V, ctau) and (M_V, Delta m)
  *.dat                        inputs (see "Input data" below)
legacy/                  earlier exploratory code (not used in the paper)
```

Comments and console messages in the scripts are in Spanish.

## Requirements

- MadGraph5_aMC@NLO 2.9.21, run with Python 3.7–3.12 (the grid scripts check
  this; set `MG5PYTHON=python3.11` or similar if needed)
- Pythia 8.3 (tested with 8.317 and 8.318), compiled with gzip support
- FastJet 3.4
- Python 3 (standard library); Matplotlib and NumPy only for the plots
- A C++17 compiler

## Workflow

All scripts read their settings from environment variables with sensible
defaults; the header of each script lists them.

### 1. Validation with pure winos

```bash
cd generation/wino
MG5=/path/to/MG5_aMC_v2_9_21 MASSES="400 600 800 1000" ./run_grid.sh
#  -> grid_summary.dat  (masses, cross sections and LHE files)

cd ../../analysis
PYTHIA8=/path/to/pythia8317 FASTJET=/path/to/fastjet ./compilar.sh
SUMMARY=../generation/wino/grid_summary.dat ./run_analysis.sh
#  -> resultados/out_m<M>.dat  (SR4High/SR4Mid efficiency vs ctau)

python3 limites.py --summary ../generation/wino/grid_summary.dat \
                   --xsec xsec_wino_13TeV.dat --lumi 137
#  -> nsig.dat, limites.dat and the comparison with the ATLAS wino limit
```

### 2. Vector quintuplet

```bash
cd generation/quintuplet
MG5=/path/to/MG5_aMC_v2_9_21 ./run_grid_quint.sh
#  default: m(V0) = 100-1000 GeV, Delta m = 166 MeV, kappa = 0, 1/xi = 1/5,
#  production of V+V-, V+-V0 and the channels with V++ (VPP=1)

cd ../../analysis
SUMMARY=../generation/quintuplet/grid_summary.dat IDLLP=9000006 \
OUTDIR=resultados_quint ./run_analysis.sh

python3 limites.py --summary ../generation/quintuplet/grid_summary.dat \
                   --resultados resultados_quint --xsec LO
python3 limites_quint.py \
    --dm-teoria dm_teoria_lambda0.01_Q1.dat dm_teoria_Q10_lambda0.01.dat
#  -> limites_ctau.dat, limites_dm.dat, ctau_limit.pdf, dm_mv_all.pdf
```

`IDLLP=9000006` is the PDG code of V<sup>+</sup> in the UFO. The quintuplet
uses the leading-order cross section from MadGraph (`--xsec LO`).

### Statistical treatment

95% CL counting in SR4High with the ATLAS model-independent upper limits:
N<sub>sig</sub> > 6.0 (observed) and 3.7 (expected) at 137 fb<sup>-1</sup>.
The 3 ab<sup>-1</sup> projections use N<sub>sig</sub> ≥ 3 (no background) and
N<sub>sig</sub> ≥ 9.5 (CL<sub>s</sub> with the SR4High background scaled
linearly, B ≈ 15).

## Input data

| File | Content | Source |
|---|---|---|
| `mapa_ATLAS137fb_Fig2a_4capas.dat` | four-layer tracklet efficiency vs. decay radius | ATLAS, JHEP 07 (2026) 152, Fig. 2(a) |
| `atlas_obs_JHEP07_2026_152_fig12a.dat`, `atlas_exp_...` | observed/expected wino limit in (mass, ctau) | same, Fig. 12(a) |
| `xsec_wino_13TeV.dat` | NLO+NLL wino cross sections at 13 TeV | LHC SUSY Cross Section Working Group |
| `dm_teoria_Q10_lambda*.dat` | one-loop mass splitting, μ = 10 M<sub>V</sub>, λ<sub>HV</sub> = 0.01, 0.1, 1 | this work |
| `dm_teoria_lambda0.01_Q1.dat` | one-loop mass splitting, μ = M<sub>V</sub>, λ<sub>HV</sub> = 0.01 | this work |

The ATLAS curves were extracted from the vector graphics of the published
figures; see the header of each file.

## Credits

- `ToyDetector-ATLAS` is adapted from
  [ToyDetector](https://github.com/NishitaDesai/ToyDetector) by N. Desai.
- The vector quintuplet model was introduced in P. Escalona et al.,
  *Vector dark matter from the 5-dimensional representation of SU(2)<sub>L</sub>*,
  JHEP 11 (2024) 162, [arXiv:2407.14010](https://arxiv.org/abs/2407.14010).
- MadGraph5_aMC@NLO, Pythia 8, FastJet and FeynRules are used under their own
  licenses and are not distributed here.

## License

The code in this repository is released under the MIT License (see `LICENSE`),
except for third-party material noted above.
