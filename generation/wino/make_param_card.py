#!/usr/bin/env python3
"""Genera una param_card de MSSM_SLHA2 para un wino puro (tipo AMSB).

A partir de una tarjeta plantilla cambia solo lo necesario:
  MASS   : m(n1) = M,  m(x1+) = M + dm
  NMIX   : n1 = wino puro, n2 = bino, n3/n4 = higgsinos
  UMIX, VMIX : identidad (x1+ = wino cargado puro)
  MSOFT  : M1 = pesado, M2 = M ;  HMIX : mu = pesado
  DECAY 1000024 : x1+ -> n1 pi+ (dominante), n1 e+ nu, n1 mu+ nu,
                  con el ancho fijado por ctau.

Uso:
  make_param_card.py plantilla.dat M[GeV] salida.dat [--dm 0.160] [--ctau 59.96]
"""
import argparse
import math
import re

HBARC_GEV_MM = 1.973269804e-13   # hbar*c en GeV*mm
HEAVY = 1.0e5                    # escala de desacoplo [GeV]
R2 = 1.0 / math.sqrt(2.0)

NMIX = {(1, 1): 0.0, (1, 2): 1.0, (1, 3): 0.0, (1, 4): 0.0,   # n1 = wino
        (2, 1): 1.0, (2, 2): 0.0, (2, 3): 0.0, (2, 4): 0.0,   # n2 = bino
        (3, 1): 0.0, (3, 2): 0.0, (3, 3): R2,  (3, 4): R2,    # higgsinos
        (4, 1): 0.0, (4, 2): 0.0, (4, 3): R2,  (4, 4): -R2}
IDENT2 = {(1, 1): 1.0, (1, 2): 0.0, (2, 1): 0.0, (2, 2): 1.0}

# Fracciones de decaimiento del wino cargado con dm ~ 160 MeV.
BR_PI, BR_E, BR_MU = 0.977, 0.020, 0.003


def fmt(x):
    return "%.7e" % x


def edit_block(lines, block, values, nkeys):
    """Reemplaza los valores de 'block'. values: {clave(s): valor}."""
    out, inside, seen = [], False, set()
    for line in lines:
        s = line.strip()
        if s.upper().startswith("BLOCK"):
            inside = (s.split()[1].upper() == block.upper()
                      and not s.upper().startswith("BLOCK QNUMBERS"))
        elif s.upper().startswith("DECAY"):
            inside = False
        elif inside and s and not s.startswith("#"):
            parts = s.split()
            try:
                key = tuple(int(p) for p in parts[:nkeys])
            except ValueError:
                key = None
            if key in values:
                comment = line[line.index("#"):].rstrip() if "#" in line else ""
                line = "      %s %s %s\n" % (" ".join(str(k) for k in key),
                                             fmt(values[key]), comment)
                seen.add(key)
        out.append(line)
    missing = set(values) - seen
    if missing:
        raise SystemExit("ERROR: en BLOCK %s faltan las entradas %s"
                         % (block, sorted(missing)))
    return out


def edit_decay(lines, pdg, width, channels):
    """Reemplaza la tabla DECAY de 'pdg' (ancho y modos)."""
    out, skipping, done = [], False, False
    for line in lines:
        s = line.strip()
        if re.match(r"(?i)^DECAY\s+%d\b" % pdg, s):
            out.append("DECAY  %d   %s\n" % (pdg, fmt(width)))
            out.append("#  BR             NDA  ID1    ID2   ...\n")
            for br, ids in channels:
                out.append("   %s   %d    %s\n"
                           % (fmt(br), len(ids), "  ".join(str(i) for i in ids)))
            skipping, done = True, True
            continue
        if skipping:
            # Saltar los modos antiguos (lineas que empiezan con un numero).
            if s and not s.startswith("#") and not s.upper().startswith(("DECAY", "BLOCK")):
                continue
            if s.startswith("#  BR"):
                continue
            skipping = False
        out.append(line)
    if not done:
        raise SystemExit("ERROR: la plantilla no tiene DECAY %d" % pdg)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("template")
    ap.add_argument("mass", type=float, help="masa del neutralino n1 [GeV]")
    ap.add_argument("output")
    ap.add_argument("--dm", type=float, default=0.160,
                    help="m(x1+) - m(n1) [GeV] (por defecto 0.160)")
    ap.add_argument("--ctau", type=float, default=59.96,
                    help="ctau de x1+ [mm] (por defecto 59.96 = 0.2 ns)")
    a = ap.parse_args()
    if a.dm <= 0.13957:
        raise SystemExit("ERROR: dm debe superar la masa del pion (0.13957 GeV)")

    with open(a.template) as f:
        lines = f.readlines()

    lines = edit_block(lines, "MASS",
                       {(1000022,): a.mass, (1000024,): a.mass + a.dm}, 1)
    lines = edit_block(lines, "NMIX", NMIX, 2)
    lines = edit_block(lines, "UMIX", IDENT2, 2)
    lines = edit_block(lines, "VMIX", IDENT2, 2)
    lines = edit_block(lines, "MSOFT", {(1,): HEAVY, (2,): a.mass}, 1)
    lines = edit_block(lines, "HMIX", {(1,): HEAVY}, 1)

    width = HBARC_GEV_MM / a.ctau
    lines = edit_decay(lines, 1000024, width,
                       [(BR_PI, (1000022, 211)),
                        (BR_E, (1000022, -11, 12)),
                        (BR_MU, (1000022, -13, 14))])

    header = ("## Wino puro: m(n1) = %.3f GeV, m(x1+) = %.3f GeV, "
              "ctau(x1+) = %.2f mm\n" % (a.mass, a.mass + a.dm, a.ctau))
    with open(a.output, "w") as f:
        f.write(header)
        f.writelines(lines)


if __name__ == "__main__":
    main()
