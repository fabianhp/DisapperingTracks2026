#!/usr/bin/env python3
"""Convierte las eficiencias de dtSR4 en senal esperada y limites de masa.

Lee:
  - grid_summary.dat            (masas simuladas; lo escribe run_grid.sh)
  - resultados/out_m<M>.dat     (eficiencias por ctau; los escribe dtSR4)
  - xsec_wino_13TeV.dat         (secciones eficaces NLO+NLL, columna 'total')

Calcula, para cada masa y cada ctau:
    N_sig = sigma * L * eficiencia
y, para cada ctau, la masa hasta la cual N_sig supera el limite S95 de la
region (interpolando log N_sig linealmente en la masa).

Escribe:
  - nsig.dat      N_sig en SR4High y SR4Mid para cada (masa, ctau)
  - limites.dat   limite de masa observado y esperado para cada ctau
  - limites.pdf / limites.png   figura en el plano (masa, ctau), con el estilo
                  de la Fig. 3 del borrador (solo si matplotlib esta instalado)

Uso:
  python3 limites.py
  python3 limites.py --summary grid_summary.dat --resultados resultados \
                     --xsec xsec_wino_13TeV.dat --lumi 137
"""
import argparse
import math
import os
import sys

C_MM_PER_NS = 299.792458

# Limites superiores al 95% CL sobre el numero de eventos de senal.
# ATLAS, JHEP 07 (2026) 152, Tabla 6.
S95 = {"High": {"obs": 6.0, "exp": 3.7},
       "Mid": {"obs": 12.4, "exp": 12.3}}

# Limites de masa del wino citados en el texto de ATLAS (ajuste completo
# con todas las regiones), para comparar: (tau [ns], observado, esperado).
ATLAS_REF = [(0.2, 670.0, 770.0), (1.0, 880.0, 1020.0)]


def read_table(path, ncol):
    rows = []
    with open(path) as f:
        for line in f:
            line = line.split("#")[0].strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) >= ncol:
                rows.append(parts)
    return rows


def read_xsec(path):
    """Devuelve lista ordenada de (masa, sigma_total_fb)."""
    out = [(float(r[0]), float(r[-1])) for r in read_table(path, 2)]
    return sorted(out)


def xsec_at(table, m):
    """Interpolacion log-lineal en la masa."""
    if m < table[0][0] or m > table[-1][0]:
        raise SystemExit("ERROR: masa %g GeV fuera de la tabla de secciones eficaces" % m)
    for (m0, x0), (m1, x1) in zip(table, table[1:]):
        if m0 <= m <= m1:
            if m == m0:
                return x0
            f = (m - m0) / (m1 - m0)
            return math.exp(math.log(x0) + f * (math.log(x1) - math.log(x0)))
    return table[-1][1]


def mass_limit(masses, nsig, s95):
    """Mayor masa con N_sig = s95, interpolando log(N) en la masa.

    Devuelve (valor, estado): estado es 'ok', 'debajo' (ninguna masa
    excluida) o 'encima' (todas las masas simuladas estan excluidas).
    """
    excl = [n > s95 for n in nsig]
    if not any(excl):
        return float("nan"), "debajo"
    if excl[-1]:
        return float("nan"), "encima"
    for i in range(len(masses) - 1, 0, -1):
        if excl[i - 1] and not excl[i]:
            n0, n1 = nsig[i - 1], nsig[i]
            if n1 <= 0.0:
                return masses[i - 1], "ok"
            f = (math.log(n0) - math.log(s95)) / (math.log(n0) - math.log(n1))
            return masses[i - 1] + f * (masses[i] - masses[i - 1]), "ok"
    return float("nan"), "debajo"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", default="grid_summary.dat")
    ap.add_argument("--resultados", default="resultados")
    ap.add_argument("--xsec", default="xsec_wino_13TeV.dat",
                    help="tabla de secciones eficaces, o LO para usar la sigma tras el "
                         "matching que dtSR4 escribe en cada out_m*.dat")
    ap.add_argument("--kfactor", type=float, default=1.0,
                    help="factor K que multiplica la seccion eficaz (por defecto 1)")
    ap.add_argument("--atlas", default="atlas_obs_JHEP07_2026_152_fig12a.dat",
                    help="limite observado de ATLAS: archivo con columnas masa[GeV] ctau[mm]")
    ap.add_argument("--atlas-exp", dest="atlas_exp",
                    default="atlas_exp_JHEP07_2026_152_fig12a.dat",
                    help="limite esperado de ATLAS, mismo formato")
    ap.add_argument("--xlabel", default=None, help="rotulo del eje x (acepta mathtext)")
    ap.add_argument("--figura", default="limites", help="nombre de la figura, sin extension")
    ap.add_argument("--lumi", type=float, default=137.0, help="fb^-1")
    ap.add_argument("--salida", default=".", help="carpeta de salida")
    a = ap.parse_args()

    use_lo = (a.xsec.upper() == "LO")
    xs_table = None if use_lo else read_xsec(a.xsec)
    sig_lo = {}
    masses = sorted({float(r[0]) for r in read_table(a.summary, 9)})

    data = {}       # masa -> lista de (ctau, effH, errH, effM, errM)
    for m in masses:
        path = os.path.join(a.resultados, "out_m%g.dat" % m)
        if not os.path.isfile(path):
            print("AVISO: falta %s, se omite m = %g GeV" % (path, m))
            continue
        data[m] = [tuple(float(x) for x in r[:5]) for r in read_table(path, 5)]
        with open(path) as f:
            for line in f:
                if line.startswith("# sigma matching [fb]"):
                    sig_lo[m] = float(line.split("=")[1])
    masses = sorted(data)
    if len(masses) < 2:
        raise SystemExit("ERROR: se necesitan al menos dos masas analizadas")
    ctaus = [row[0] for row in data[masses[0]]]
    for m in masses:
        if len(data[m]) != len(ctaus):
            raise SystemExit("ERROR: la malla de ctau difiere entre masas")

    if use_lo:
        xs = {m: a.kfactor * sig_lo[m] for m in masses}
    else:
        xs = {m: a.kfactor * xsec_at(xs_table, m) for m in masses}

    # ---- N_sig por (masa, ctau)
    nsig = {"High": {}, "Mid": {}}
    with open(os.path.join(a.salida, "nsig.dat"), "w") as f:
        f.write("# N_sig = sigma * L * eficiencia,  L = %g fb^-1\n" % a.lumi)
        f.write("# masa[GeV]  sigma[fb]  tau[ns]  ctau[mm]  N_SR4High  err  N_SR4Mid  err\n")
        for m in masses:
            for k, (ct, eh, dh, em, dm) in enumerate(data[m]):
                norm = xs[m] * a.lumi
                nsig["High"].setdefault(k, []).append(norm * eh)
                nsig["Mid"].setdefault(k, []).append(norm * em)
                f.write("%g  %.5g  %.5g  %.5g  %.5g  %.3g  %.5g  %.3g\n"
                        % (m, xs[m], ct / C_MM_PER_NS, ct,
                           norm * eh, norm * dh, norm * em, norm * dm))

    # ---- Limites de masa por ctau
    lim = []
    with open(os.path.join(a.salida, "limites.dat"), "w") as f:
        f.write("# Limite de masa [GeV] al 95% CL por conteo en una region.\n")
        f.write("# nan = sin limite dentro de las masas simuladas (%g-%g GeV):\n"
                % (masses[0], masses[-1]))
        f.write("#   '<' ninguna masa excluida, '>' todas excluidas.\n")
        f.write("# 'mejor' = observado de la region con mejor limite esperado.\n")
        f.write("# tau[ns]  ctau[mm]  High_obs  High_exp  Mid_obs  Mid_exp  mejor_obs  region\n")
        for k, ct in enumerate(ctaus):
            row = {}
            for reg in ("High", "Mid"):
                for kind in ("obs", "exp"):
                    row[reg, kind] = mass_limit(masses, nsig[reg][k], S95[reg][kind])
            best = None
            for reg in ("High", "Mid"):
                v, st = row[reg, "exp"]
                score = masses[-1] + 1 if st == "encima" else (v if st == "ok" else -1)
                if best is None or score > best[0]:
                    best = (score, reg)
            breg = best[1] if best[0] > 0 else "-"
            bobs = row[breg, "obs"] if breg != "-" else (float("nan"), "debajo")

            def s(x):
                v, st = x
                return "%.1f" % v if st == "ok" else ("nan>" if st == "encima" else "nan<")
            f.write("%.5g  %.5g  %s  %s  %s  %s  %s  %s\n"
                    % (ct / C_MM_PER_NS, ct, s(row["High", "obs"]), s(row["High", "exp"]),
                       s(row["Mid", "obs"]), s(row["Mid", "exp"]), s(bobs), breg))
            lim.append((ct / C_MM_PER_NS, row))

    # ---- Resumen en pantalla
    print("Masas analizadas [GeV]:", " ".join("%g" % m for m in masses))
    refs = [] if use_lo else ATLAS_REF
    if use_lo:
        print("Seccion eficaz: LO tras el matching, K = %g" % a.kfactor)
        print("\nLimites de masa en SR4High (conteo en un bin):")
        print("  tau[ns]   ctau[mm]   observado   esperado")
        for tau in (0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 1.0):
            k = min(range(len(ctaus)), key=lambda i: abs(math.log(ctaus[i] / (tau * C_MM_PER_NS))))
            row = lim[k][1]

            def s2(x):
                v, st = x
                return "%7.0f" % v if st == "ok" else ("  > %g" % masses[-1] if st == "encima"
                                                       else "  < %g" % masses[0])
            print("  %6.3f   %7.1f    %s     %s"
                  % (lim[k][0], ctaus[k], s2(row["High", "obs"]), s2(row["High", "exp"])))
    else:
        print("\nComparacion con ATLAS (SR4High, conteo en un bin):")
        print("  tau[ns]   recast obs   recast exp   ATLAS obs   ATLAS exp")
    for tau, aobs, aexp in refs:
        k = min(range(len(ctaus)), key=lambda i: abs(math.log(ctaus[i] / (tau * C_MM_PER_NS))))
        row = lim[k][1]

        def s(x):
            v, st = x
            return "%7.0f" % v if st == "ok" else ("   > %g" % masses[-1] if st == "encima"
                                                   else "   < %g" % masses[0])
        print("  %5.2f    %s      %s      %7.0f     %7.0f"
              % (lim[k][0], s(row["High", "obs"]), s(row["High", "exp"]), aobs, aexp))
    print("\nEscritos: nsig.dat, limites.dat")

    # ---- Figura
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib no esta instalado: se omite la figura.")
        return

    # Estilo de la Fig. 3 del borrador: matplotlib por defecto con tipografia
    # Computer Modern, marco completo, rejilla gris, leyenda arriba a la izquierda.
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["cmr10", "Computer Modern Roman", "DejaVu Serif"],
        "mathtext.fontset": "cm",
        "axes.formatter.use_mathtext": True,
        "axes.unicode_minus": False,
        "font.size": 12,
        "axes.labelsize": 16,
        "legend.fontsize": 12,
    })
    fig, ax = plt.subplots(figsize=(6.4, 4.8))

    ctau_mm = [t * C_MM_PER_NS for t, r in lim]

    def curve(kind):
        # Fuera del rango simulado el valor es nan y la curva se corta.
        return [r["High", kind][0] if r["High", kind][1] == "ok" else float("nan")
                for t, r in lim]

    # ATLAS en azul, recast en naranja; observado con trazo lleno o a
    # trazos, esperado punteado.
    if not use_lo:
        for path, ls, label in ((a.atlas, "--", r"ATLAS observed (137 fb$^{-1}$)"),
                                (a.atlas_exp, ":", "ATLAS expected")):
            if path and os.path.isfile(path):
                pts = [(float(r[0]), float(r[1])) for r in read_table(path, 2)]
                ax.plot([q[0] for q in pts], [q[1] for q in pts], ls, color="C0", lw=1.5,
                        label=label)
    ax.plot(curve("obs"), ctau_mm, "-", color="C1", lw=1.8,
            label="Our recast, observed" if not use_lo else "Observed")
    ax.plot(curve("exp"), ctau_mm, ":", color="C1", lw=1.8,
            label="Our recast, expected" if not use_lo else "Expected")

    ax.set_yscale("log")
    ax.set_ylim(1.0, 3.0e3)
    ax.set_xlim(masses[0], max(masses[-1], 1000.0))
    ax.set_ylabel(r"$c\tau$ [mm]")
    ax.set_xlabel(a.xlabel or (r"$M_{V}$ [GeV]" if use_lo
                               else r"$m_{\tilde{\chi}^{\pm}}$ [GeV]"))
    ax.grid(True, which="major", color="#b0b0b0", lw=0.8)
    ax.legend(loc="upper left", fontsize=11)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        out = os.path.join(a.salida, a.figura + "." + ext)
        fig.savefig(out, dpi=200)
        print("Escrita la figura:", out)


if __name__ == "__main__":
    sys.exit(main())
