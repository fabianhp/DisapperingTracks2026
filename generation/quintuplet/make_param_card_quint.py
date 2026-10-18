#!/usr/bin/env python3
"""Genera una param_card del modelo VQUINT_UFO para una masa del quintuplete.

Cambia en la plantilla (la param_card por defecto del modelo):
  MASS   : m(V0) = M, m(V+) = M + dm, m(V++) = M + 4*dm   (dm2 = 3*dm)
           u, d, s, c sin masa (esquema de 4 sabores, necesario para el matching)
  YUKAWA : u, d, s, c a cero, en consistencia con las masas
  FRBLOCK: lambdav, kappa, xiinv, alpha3 solo si se piden por linea de comandos
  DECAY  : V0 estable;  V+ -> V0 pi+ ;  V++ -> V+ pi+   (100% cada uno)

El pion no existe en el UFO: los decaimientos van en la tabla DECAY con el
codigo PDG 211 y los ejecuta Pythia, igual que en la validacion con winos.

Anchos por defecto (limite de multiplete pesado, f_pi = 130.4 MeV):
  Gamma(V+  -> V0 pi+) = 6 GF^2 fpi^2 Vud^2 dm^3  sqrt(1 - mpi^2/dm^2 ) / pi
  Gamma(V++ -> V+ pi+) = 4 GF^2 fpi^2 Vud^2 dm2^3 sqrt(1 - mpi^2/dm2^2) / pi
El de V+ se puede fijar a mano con --ctau [mm]; el analisis repondera en ctau,
asi que su valor solo importa para el control con vertices sorteados.

Uso:
  make_param_card_quint.py plantilla.dat M[GeV] salida.dat [--dm 0.166] [--ctau mm]
                           [--lambdav x] [--kappa x] [--xiinv x] [--alpha3 x]
"""
import argparse
import math
import re

HBARC_GEV_MM = 1.973269804e-13
GF, FPI, VUD, MPI = 1.16637e-5, 0.1304, 0.9743, 0.13957
ID_VPP, ID_VP, ID_V0 = 9000005, 9000006, 9000007


def fmt(x):
    return "%.7e" % x


def width_pion(coef, dm):
    if dm <= MPI:
        return 0.0
    return coef * GF**2 * FPI**2 * VUD**2 * dm**3 * math.sqrt(1.0 - (MPI / dm)**2) / math.pi


def edit_block(lines, block, values):
    out, inside, seen = [], False, set()
    for line in lines:
        s = line.strip()
        low = s.lower()
        if low.startswith("block"):
            inside = (low.split()[1] == block.lower())
        elif low.startswith("decay"):
            inside = False
        elif inside and s and not s.startswith("#"):
            try:
                key = int(s.split()[0])
            except ValueError:
                key = None
            if key in values:
                comment = line[line.index("#"):].rstrip() if "#" in line else ""
                line = "  %d %s %s\n" % (key, fmt(values[key]), comment)
                seen.add(key)
        out.append(line)
    missing = set(values) - seen
    if missing:
        raise SystemExit("ERROR: en Block %s faltan las entradas %s" % (block, sorted(missing)))
    return out


def edit_decay(lines, pdg, width, channels):
    out, done = [], False
    for line in lines:
        if re.match(r"(?i)^\s*DECAY\s+%d\b" % pdg, line):
            out.append("DECAY %d %s\n" % (pdg, fmt(width)))
            for br, ids in channels:
                out.append("   %s   %d   %s\n" % (fmt(br), len(ids), "  ".join(str(i) for i in ids)))
            done = True
            continue
        out.append(line)
    if not done:
        raise SystemExit("ERROR: la plantilla no tiene DECAY %d" % pdg)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("template")
    ap.add_argument("mass", type=float, help="masa de V0 [GeV]")
    ap.add_argument("output")
    ap.add_argument("--dm", type=float, default=0.166, help="m(V+) - m(V0) [GeV]")
    ap.add_argument("--ctau", type=float, default=None, help="ctau de V+ [mm]")
    for p in ("lambdav", "kappa", "xiinv", "alpha3"):
        ap.add_argument("--" + p, type=float, default=None)
    a = ap.parse_args()
    if a.dm <= MPI:
        raise SystemExit("ERROR: dm debe superar la masa del pion (%.5f GeV)" % MPI)

    with open(a.template) as f:
        lines = f.readlines()

    dm2 = 3.0 * a.dm
    lines = edit_block(lines, "mass", {1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0,
                                       ID_V0: a.mass, ID_VP: a.mass + a.dm,
                                       ID_VPP: a.mass + a.dm + dm2})
    lines = edit_block(lines, "yukawa", {1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0})
    fr = {}
    for code, name in ((1, "alpha3"), (2, "kappa"), (3, "xiinv"), (4, "lambdav")):
        if getattr(a, name) is not None:
            fr[code] = getattr(a, name)
    if fr:
        lines = edit_block(lines, "frblock", fr)

    w_vp = HBARC_GEV_MM / a.ctau if a.ctau else width_pion(6.0, a.dm)
    w_vpp = width_pion(4.0, dm2)
    lines = edit_decay(lines, ID_V0, 0.0, [])
    lines = edit_decay(lines, ID_VP, w_vp, [(1.0, (ID_V0, 211))])
    lines = edit_decay(lines, ID_VPP, w_vpp, [(1.0, (ID_VP, 211))])

    header = ("## Quintuplete vectorial: m(V0) = %.3f, m(V+) = %.3f, m(V++) = %.3f GeV\n"
              "## ctau(V+) = %.3f mm, ctau(V++) = %.4f mm\n"
              % (a.mass, a.mass + a.dm, a.mass + a.dm + dm2,
                 HBARC_GEV_MM / w_vp, HBARC_GEV_MM / w_vpp))
    with open(a.output, "w") as f:
        f.write(header)
        f.writelines(lines)


if __name__ == "__main__":
    main()
