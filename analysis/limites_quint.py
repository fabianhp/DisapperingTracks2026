#!/usr/bin/env python3
"""Limites del quintuplete vectorial en los planos (M_V, ctau) y (M_V, Delta m).

Se corre DESPUES de   python3 limites.py --xsec LO   y lee su nsig.dat
(N_sig en SR4High para cada masa y ctau, a 137 fb^-1).

Curvas que calcula (95% CL, conteo en SR4High):
  - Run 2 observado  : N_sig > 6.0   (S95 observado de ATLAS)
  - Run 2 esperado   : N_sig > 3.7   (S95 esperado de ATLAS)
  - Proyeccion B = 0 : N_sig(L) > 3.0
  - Proyeccion B = 0.68*L/137 : umbral CLs con n_obs = B   (9.5 para 3 ab^-1)

Plano (M_V, Delta m): para cada masa y cada Delta m se calcula
    ctau = hbar c / Gamma(V+- -> V0 pi+-)        (Ec. de la anchura del paper)
y se evalua N_sig en ese ctau interpolando en la malla de ctau.

Escribe:
  limites_ctau.dat, limites_dm.dat
  ctau_limit.pdf/.png      plano (M_V, ctau)
  dm_mv_all.pdf/.png       plano (M_V, Delta m)

Uso:
  python3 limites_quint.py
  python3 limites_quint.py --dm-teoria dm_teoria.dat      # columnas: M_V[GeV]  Delta m[GeV]
  python3 limites_quint.py --lumi-proy 3000 --dm-min 0.14 --dm-max 0.30
"""
import argparse
import math
import os
import sys

HBARC_GEV_MM = 1.973269804e-13
GF = 1.1663787e-5        # GeV^-2
FPI = 0.1302             # GeV
VUD = 0.97435
MPI = 0.13957039         # GeV

LUMI_RUN2 = 137.0
S95_OBS, S95_EXP = 6.0, 3.7      # ATLAS, SR4High
BKG_RUN2 = 0.68                  # ATLAS, SR4High


def width_pi(mv0, dm):
    """Gamma(V+- -> V0 pi+-) [GeV]; mv0 = masa de V0, dm = M(V+-) - M(V0)."""
    if dm <= MPI:
        return 0.0
    mp = mv0 + dm
    g4_over_mw4 = 32.0 * GF * GF
    s = 2.0 * mp - dm
    pref = (g4_over_mw4 * FPI ** 2 * VUD ** 2 * dm ** 2 * s ** 2
            * math.sqrt((dm ** 2 - MPI ** 2) * (s ** 2 - MPI ** 2))
            / (512.0 * math.pi * mp ** 5 * (mp - dm) ** 2))
    brk = (2.0 * (mp - dm) ** 2 * (5.0 * mp ** 2 - MPI ** 2)
           + (mp ** 2 - MPI ** 2) ** 2 + (mp - dm) ** 4)
    return pref * brk


MMU = 0.1056584          # GeV


def width_lep(dm):
    """Gamma(V+- -> V0 e nu) + Gamma(V+- -> V0 mu nu) [GeV], limite M_V >> dm.
    Comprobado contra MadGraph con el UFO (dm = 166 MeV: 2.47e-16 GeV)."""
    if dm <= 0.0:
        return 0.0
    ge = 6.0 * GF * GF * dm ** 5 / (15.0 * math.pi ** 3)
    gm = 0.0
    x = MMU / dm
    if x < 1.0:
        r = math.sqrt(1.0 - x * x)
        gm = ge * r * (1.0 - 4.5 * x ** 2 - 4.0 * x ** 4
                       + 15.0 * x ** 4 / (2.0 * r) * math.atanh(r))
    return ge + gm


LEPTONS = True


def width_tot(mv0, dm):
    return width_pi(mv0, dm) + (width_lep(dm) if LEPTONS else 0.0)


def ctau_mm(mv0, dm):
    g = width_tot(mv0, dm)
    return float("inf") if g <= 0.0 else HBARC_GEV_MM / g


def pois_cdf(n, mu):
    t, s = math.exp(-mu), 0.0
    for k in range(n + 1):
        s += t
        t *= mu / (k + 1)
    return s


def cls_threshold(b):
    """Senal excluida al 95% CL (CLs) si se observa n = b (redondeado)."""
    if b <= 0.0:
        return 3.0
    n = int(round(b))
    lo, hi = 0.0, 50.0 + 5.0 * b
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if pois_cdf(n, mid + b) / pois_cdf(n, b) > 0.05:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def read_nsig(path):
    data = {}
    with open(path) as f:
        for line in f:
            line = line.split("#")[0].split()
            if len(line) >= 5:
                data.setdefault(float(line[0]), []).append((float(line[3]), float(line[4])))
    masses = sorted(data)
    ctaus = [c for c, _ in data[masses[0]]]
    n = {m: [x for _, x in data[m]] for m in masses}
    return masses, ctaus, n


def n_at_ctau(ctaus, nrow, ct):
    """N_sig en un ctau cualquiera: interpolacion log-log en la malla.
    Por encima de la malla se extrapola como 1/ctau; por debajo vale 0."""
    if not math.isfinite(ct):
        return 0.0
    if ct < ctaus[0]:
        return 0.0
    if ct >= ctaus[-1]:
        return nrow[-1] * ctaus[-1] / ct
    for i in range(len(ctaus) - 1):
        if ctaus[i] <= ct <= ctaus[i + 1]:
            a, b = nrow[i], nrow[i + 1]
            if a <= 0.0 or b <= 0.0:
                f = (ct - ctaus[i]) / (ctaus[i + 1] - ctaus[i])
                return a + f * (b - a)
            f = math.log(ct / ctaus[i]) / math.log(ctaus[i + 1] / ctaus[i])
            return math.exp(math.log(a) + f * (math.log(b) - math.log(a)))
    return 0.0


def mass_limit(masses, nsig, thr):
    """Mayor masa con N_sig = thr (interpolando log N en la masa).
    nan si no se excluye ninguna masa; +inf si se excluyen todas las simuladas."""
    excl = [x > thr for x in nsig]
    if not any(excl):
        return float("nan")
    if excl[-1]:
        return float("inf")
    for i in range(len(masses) - 1, 0, -1):
        if excl[i - 1] and not excl[i]:
            n0, n1 = nsig[i - 1], nsig[i]
            if n1 <= 0.0:
                return masses[i - 1]
            f = math.log(n0 / thr) / math.log(n0 / n1)
            return masses[i - 1] + f * (masses[i] - masses[i - 1])
    return float("nan")


def n_at(masses, ctaus, n, m, ct):
    """N_sig en una masa y ctau cualesquiera (log-lineal en la masa)."""
    if m <= masses[0]:
        return n_at_ctau(ctaus, n[masses[0]], ct)
    if m >= masses[-1]:
        return n_at_ctau(ctaus, n[masses[-1]], ct)
    for m0, m1 in zip(masses, masses[1:]):
        if m0 <= m <= m1:
            a, b = n_at_ctau(ctaus, n[m0], ct), n_at_ctau(ctaus, n[m1], ct)
            f = (m - m0) / (m1 - m0)
            if a <= 0.0 or b <= 0.0:
                return a + f * (b - a)
            return math.exp(math.log(a) + f * (math.log(b) - math.log(a)))
    return 0.0


def fmt(v, masses):
    if math.isnan(v):
        return "nan<"
    if math.isinf(v):
        return "nan>"
    return "%.1f" % v


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nsig", default="nsig.dat")
    ap.add_argument("--lumi-proy", dest="lproy", type=float, default=3000.0, help="fb^-1")
    ap.add_argument("--dm-min", dest="dmmin", type=float, default=0.100)
    ap.add_argument("--dm-max", dest="dmmax", type=float, default=0.300)
    ap.add_argument("--dm-teoria", dest="teoria", nargs="+", default=None,
                    help="uno o dos archivos con el split radiativo (M_V[GeV]  Delta m[GeV]); "
                         "con dos (p. ej. Q = 1 y Q = 10) se dibuja la banda entre ambos")
    ap.add_argument("--etiquetas", nargs="+", default=None,
                    help="rotulo de cada archivo de --dm-teoria (mathtext); si se dan, "
                         "cada archivo se dibuja como una curva en vez de una banda")
    ap.add_argument("--etiqueta-banda", dest="etiqueta_banda",
                    default=r"$\Delta m(M_V)$, $\lambda_{HV}=0.01$, $\mu = M_V$ to $10\,M_V$",
                    help="rotulo de la banda cuando se dan dos archivos en --dm-teoria")
    ap.add_argument("--dm-logx", dest="dmlogx", action="store_true",
                    help="eje de masas logaritmico en el plano (M_V, Delta m)")
    ap.add_argument("--sin-leptones", dest="nolep", action="store_true",
                    help="usar solo la anchura a pion (sin los modos leptonicos)")
    ap.add_argument("--salida", default=".")
    a = ap.parse_args()

    global LEPTONS
    LEPTONS = not a.nolep
    masses, ctaus, n = read_nsig(a.nsig)
    scale = a.lproy / LUMI_RUN2
    bproy = BKG_RUN2 * scale
    thr_b = cls_threshold(bproy)
    # (etiqueta, factor de luminosidad, umbral)
    curves = [
        ("run2_obs", 1.0, S95_OBS),
        ("run2_exp", 1.0, S95_EXP),
        ("proy_B0", scale, 3.0),
        ("proy_B", scale, thr_b),
    ]
    print("Masas [GeV]:", " ".join("%g" % m for m in masses))
    print("Proyeccion: L = %g fb^-1, B = %.1f  ->  N_sig > %.1f (CLs);  B = 0 -> N_sig > 3.0"
          % (a.lproy, bproy, thr_b))

    # ------------------------------------------------ plano (M_V, ctau)
    lim_ct = {c[0]: [] for c in curves}
    for k, ct in enumerate(ctaus):
        col = [n[m][k] for m in masses]
        for name, sc, thr in curves:
            lim_ct[name].append(mass_limit(masses, [sc * x for x in col], thr))
    with open(os.path.join(a.salida, "limites_ctau.dat"), "w") as f:
        f.write("# Limite de masa M_V [GeV] al 95% CL en SR4High, para cada ctau.\n")
        f.write("# nan< ninguna masa excluida; nan> todas las simuladas excluidas.\n")
        f.write("# ctau[mm]  run2_obs  run2_exp  proy_B0  proy_B%.0f\n" % bproy)
        for k, ct in enumerate(ctaus):
            f.write("%.5g  %s\n" % (ct, "  ".join(fmt(lim_ct[c[0]][k], masses) for c in curves)))

    # ------------------------------------------------ plano (M_V, Delta m)
    dms = []
    x = a.dmmin if LEPTONS else max(a.dmmin, MPI + 1.0e-5)
    while x <= a.dmmax + 1e-9:
        dms.append(x)
        # paso fino justo sobre el umbral del pion, donde la anchura cambia rapido
        x += 0.0001 if MPI - 0.0002 < x < MPI + 0.004 else 0.001
    lim_dm = {c[0]: [] for c in curves}
    for dm in dms:
        col = [n_at_ctau(ctaus, n[m], ctau_mm(m, dm)) for m in masses]
        for name, sc, thr in curves:
            lim_dm[name].append(mass_limit(masses, [sc * x for x in col], thr))
    mref = masses[len(masses) // 2]
    with open(os.path.join(a.salida, "limites_dm.dat"), "w") as f:
        f.write("# Limite de masa M_V [GeV] al 95% CL en SR4High, para cada Delta m.\n")
        f.write("# ctau calculado con Gamma(pi)%s; la columna ctau es para M_V = %g GeV.\n"
                % (" + Gamma(e nu) + Gamma(mu nu)" if LEPTONS else "", mref))
        f.write("# dm[GeV]  ctau[mm]  run2_obs  run2_exp  proy_B0  proy_B%.0f\n" % bproy)
        for i, dm in enumerate(dms):
            f.write("%.6f  %.5g  %s\n" % (dm, ctau_mm(mref, dm),
                                           "  ".join(fmt(lim_dm[c[0]][i], masses) for c in curves)))

    # ------------------------------------------------ numeros para el texto
    def best(vals, xs):
        ok = [(v, x) for v, x in zip(vals, xs) if math.isfinite(v)]
        if any(math.isinf(v) for v in vals):
            return None
        return max(ok) if ok else (float("nan"), float("nan"))

    def at(vals, xs, x0):
        i = min(range(len(xs)), key=lambda j: abs(math.log(xs[j] / x0)))
        return vals[i]

    names = {"run2_obs": "Run 2 observado", "run2_exp": "Run 2 esperado",
             "proy_B0": "Proyeccion, B = 0", "proy_B": "Proyeccion, B = %.0f" % bproy}
    print("\nPlano (M_V, ctau):")
    print("  %-22s  M_V max [GeV]  en ctau [mm]   M_V(ctau = 10 mm)" % "")
    for name, _, _ in curves:
        b = best(lim_ct[name], ctaus)
        v10 = at(lim_ct[name], ctaus, 10.0)
        if b is None:
            print("  %-22s  > %g (fuera de la malla de masas)   %s"
                  % (names[name], masses[-1], fmt(v10, masses)))
        else:
            print("  %-22s  %10.0f   %10.0f        %s" % (names[name], b[0], b[1], fmt(v10, masses)))
    print("\nPlano (M_V, Delta m):")
    print("  %-22s  M_V max [GeV]  en dm [MeV]   M_V(dm = 166 MeV)" % "")
    for name, _, _ in curves:
        b = best(lim_dm[name], dms)
        v166 = at(lim_dm[name], dms, 0.166)
        if b is None:
            print("  %-22s  > %g (fuera de la malla de masas)   %s"
                  % (names[name], masses[-1], fmt(v166, masses)))
        else:
            print("  %-22s  %10.0f   %10.1f        %s"
                  % (names[name], b[0], 1e3 * b[1], fmt(v166, masses)))

    # ------------------------------------------------ cruce con el split radiativo
    teos = []
    for path in (a.teoria or []):
        teo = []
        with open(path) as f:
            for line in f:
                p = line.split("#")[0].replace(",", " ").split()
                if len(p) >= 2:
                    teo.append((float(p[0]), float(p[1])))
        teo.sort()
        teos.append(teo)
        print("\nSplit radiativo (%s):" % path)
        fine = []
        for (m0, d0), (m1, d1) in zip(teo, teo[1:]):     # interpolacion log-log
            for j in range(40):
                t = j / 40.0
                fine.append((math.exp(math.log(m0) + t * math.log(m1 / m0)),
                             math.exp(math.log(d0) + t * math.log(d1 / d0))))
        fine.append(teo[-1])
        fine = [(m, d) for m, d in fine if masses[0] <= m <= masses[-1]]
        for m, d in fine:
            if d <= MPI:
                print("  Delta m cruza m_pi en M_V = %.0f GeV (ctau = %.2g m alli)"
                      % (m, 1e-3 * ctau_mm(m, d)))
                break
        for name, sc, thr in curves:
            ex = [sc * n_at(masses, ctaus, n, m, ctau_mm(m, d)) > thr for m, d in fine]
            if not any(ex):
                print("  %-22s  no excluida en ninguna masa" % names[name])
                continue
            segs, start = [], None
            for (m, d), e in zip(fine, ex):
                if e and start is None:
                    start = m
                if not e and start is not None:
                    segs.append((start, m))
                    start = None
            if start is not None:
                segs.append((start, fine[-1][0]))
            print("  %-22s  excluida para M_V en %s GeV"
                  % (names[name], ", ".join("[%.0f, %.0f]" % sg for sg in segs)))

    print("\nEscritos: limites_ctau.dat, limites_dm.dat")

    # ------------------------------------------------ figuras
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib no esta instalado: se omiten las figuras.")
        return
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["cmr10", "Computer Modern Roman", "DejaVu Serif"],
        "mathtext.fontset": "cm",
        "axes.formatter.use_mathtext": True,
        "axes.unicode_minus": False,
        "font.size": 12,
        "axes.labelsize": 16,
        "legend.fontsize": 11,
    })
    lab_l = ("%g ab$^{-1}$" % (a.lproy / 1000.0)) if a.lproy >= 1000 else ("%g fb$^{-1}$" % a.lproy)
    style = {
        "run2_obs": ("-", "C1", 1.8, r"Run 2 observed (137 fb$^{-1}$)"),
        "run2_exp": (":", "C1", 1.8, "Run 2 expected"),
        "proy_B0": ("--", "C0", 1.6, r"%s, $B=0$ ($N_{\rm sig}\geq 3$)" % lab_l),
        "proy_B": ("-.", "C0", 1.6, r"%s, $B\approx %.0f$ ($N_{\rm sig}\geq %.1f$)"
                   % (lab_l, bproy, thr_b)),
    }

    def clean(vals):
        return [v if math.isfinite(v) else float("nan") for v in vals]

    def draw(ax, lims, ys):
        for name, _, _ in curves:
            ls, col, lw, lab = style[name]
            xs = clean(lims[name])
            ax.plot(xs, ys, ls, color=col, lw=lw, label=lab)
            if name == "run2_obs":
                # 'todas las masas excluidas' se rellena hasta la ultima masa simulada
                edge = [masses[-1] if math.isinf(v) else (masses[0] if math.isnan(v) else v)
                        for v in lims[name]]
                ax.fill_betweenx(ys, masses[0], edge, color=col, alpha=0.15, lw=0)
        ax.set_xlim(masses[0], max(masses[-1], 1000.0))
        ax.set_xlabel(r"$M_{V}$ [GeV]")
        ax.grid(True, which="major", color="#b0b0b0", lw=0.8)

    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    draw(ax, lim_ct, ctaus)
    ax.set_yscale("log")
    ax.set_ylim(1.0, 3.0e3)
    ax.set_ylabel(r"$c\tau$ [mm]")
    ax.legend(loc="upper left")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        out = os.path.join(a.salida, "ctau_limit." + ext)
        fig.savefig(out, dpi=200)
        print("Escrita la figura:", out)

    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    draw(ax, lim_dm, dms)
    if teos and a.etiquetas:
        for teo, lab, ls in zip(teos, a.etiquetas, ["--", "-.", ":"] * 3):
            pts = [(m, d) for m, d in teo if masses[0] <= m <= max(masses[-1], 1000.0)]
            ax.plot([m for m, d in pts], [d for m, d in pts], ls, color="C3", lw=1.6, label=lab)
    elif len(teos) == 2:
        import numpy as np
        xs = np.linspace(masses[0], masses[-1], 400)

        def at_m(teo, x):
            lm = [math.log(m) for m, d in teo]
            ld = [math.log(d) for m, d in teo]
            return np.exp(np.interp(np.log(x), lm, ld))
        ax.fill_between(xs, at_m(teos[0], xs), at_m(teos[1], xs), color="C3", alpha=0.25, lw=0,
                        label=a.etiqueta_banda)
        for teo in teos:
            ax.plot(xs, at_m(teo, xs), "--", color="C3", lw=1.2)
    elif teos:
        ax.plot([m for m, d in teos[0]], [d for m, d in teos[0]], "--", color="C3", lw=1.6,
                label=r"One-loop $\Delta m(M_V)$")
    ax.axhline(MPI, color="0.4", lw=0.8, ls=":")
    if a.dmlogx:
        from matplotlib.ticker import FixedLocator, FixedFormatter, NullFormatter
        ax.set_xscale("log")
        ticks = [t for t in (100, 200, 300, 500, 700, 1000, 1500, 2000) if masses[0] <= t <= masses[-1]]
        ax.xaxis.set_major_locator(FixedLocator(ticks))
        ax.xaxis.set_major_formatter(FixedFormatter(["%d" % t for t in ticks]))
        ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_ylim(a.dmmin, a.dmmax)
    ax.set_ylabel(r"$\Delta m$ [GeV]")
    ax.legend(loc="upper right")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        out = os.path.join(a.salida, "dm_mv_all." + ext)
        fig.savefig(out, dpi=200)
        print("Escrita la figura:", out)


if __name__ == "__main__":
    sys.exit(main())
