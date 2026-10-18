#!/bin/bash
# ------------------------------------------------------------------
# Genera con MadGraph5_aMC@NLO la malla de masas del wino puro
# (validacion AMSB del recast de ATLAS disappearing track).
#
# Procesos: p p > x1+ x1-, x1+ n1, x1- n1  (+ 0..NJETMAX jets, MLM)
# Una sola carpeta de proceso; por cada masa se genera la param_card
# con make_param_card.py y se lanza una corrida.
#
# Uso:
#   MG5=/ruta/a/MG5_aMC_v2_9_21 ./run_grid.sh
#   MASSES="400 600" NEVENTS=20000 ./run_grid.sh
#
# Todas las variables de abajo se pueden sobreescribir desde el entorno.
# ------------------------------------------------------------------
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"

MG5="${MG5:-$HOME/MG5_aMC_v2_9_21}"                  # carpeta de MadGraph
MG5PYTHON="${MG5PYTHON:-}"   # interprete para MadGraph (p.ej. python3.11); vacio = el de mg5_aMC
TEMPLATE="${TEMPLATE:-$HERE/param_card_template.dat}"
PROCDIR="${PROCDIR:-$HERE/wino_dt}"                  # carpeta de proceso
CARDS="${CARDS:-$HERE/cards}"                        # param cards generadas
SUMMARY="${SUMMARY:-$HERE/grid_summary.dat}"

MASSES="${MASSES:-100 200 300 400 500 600 700 800 900 1000 1100}"   # m(n1) [GeV]
DM="${DM:-0.160}"            # m(x1+) - m(n1) [GeV]
CTAU="${CTAU:-59.96}"        # ctau de x1+ [mm]  (0.2 ns)
NEVENTS="${NEVENTS:-50000}"
EBEAM="${EBEAM:-6500}"       # GeV por haz
NJETMAX="${NJETMAX:-2}"      # jets extra en el elemento de matriz (0, 1 o 2)
NCORES="${NCORES:-}"         # vacio = todos los nucleos

# Escala de matching. "auto": qCut = M/4 y xqcut = M/6.
# Un numero fija xqcut [GeV] para todas las masas, y qCut = 1.5 * xqcut
# salvo que se de QCUT.
XQCUT="${XQCUT:-auto}"
QCUT="${QCUT:-}"

[ -x "$MG5/bin/mg5_aMC" ] || { echo "ERROR: no encuentro $MG5/bin/mg5_aMC (definir MG5=...)"; exit 1; }
[ -f "$TEMPLATE" ] || { echo "ERROR: falta la plantilla $TEMPLATE"; exit 1; }
mkdir -p "$CARDS"

# MadGraph 2.9.x no funciona con Python 3.13 o posterior (falla al importar
# modelos o al escribir la carpeta de proceso). Se comprueba antes de empezar.
PYMG="${MG5PYTHON:-python3}"
command -v "$PYMG" >/dev/null 2>&1 || { echo "ERROR: no encuentro el interprete '$PYMG' (definir MG5PYTHON=...)"; exit 1; }
PYVER=$("$PYMG" -c 'import sys; print("%d.%d" % sys.version_info[:2])')
if ! "$PYMG" -c 'import sys; sys.exit(0 if (3,7) <= sys.version_info[:2] <= (3,12) else 1)'; then
  echo "ERROR: '$PYMG' es Python $PYVER y MadGraph 2.9.x necesita Python entre 3.7 y 3.12."
  echo "       Interpretes disponibles en este equipo:"
  for c in python3.12 python3.11 python3.10 python3.9 python3.8; do
    command -v "$c" >/dev/null 2>&1 && echo "         $c  ->  $(command -v "$c")"
  done
  echo "       Volver a correr anteponiendo, por ejemplo:  MG5PYTHON=python3.11"
  exit 1
fi
"$PYMG" -c 'import six' 2>/dev/null || { echo "ERROR: a '$PYMG' le falta el modulo six (instalar con: $PYMG -m pip install --user six)"; exit 1; }
MG5PYTHON="$PYMG"

mg5() { if [ -n "$MG5PYTHON" ]; then "$MG5PYTHON" "$MG5/bin/mg5_aMC" "$@"; else "$MG5/bin/mg5_aMC" "$@"; fi; }

# ---------------- 1. Carpeta de proceso (una sola vez) ----------------
if [ ! -d "$PROCDIR" ]; then
  PROC="$HERE/proc_card_wino.dat"
  {
    echo "import model MSSM_SLHA2"
    echo "define p = g u c d s u~ c~ d~ s~"
    echo "define j = g u c d s u~ c~ d~ s~"
    echo "generate p p > x1+ x1- @0"
    echo "add process p p > x1+ n1 @0"
    echo "add process p p > x1- n1 @0"
    JETS=""
    for n in $(seq 1 "$NJETMAX"); do
      JETS="$JETS j"
      echo "add process p p > x1+ x1-$JETS @$n"
      echo "add process p p > x1+ n1$JETS @$n"
      echo "add process p p > x1- n1$JETS @$n"
    done
    echo "output $PROCDIR"
  } > "$PROC"
  mg5 "$PROC" 2>&1 | tee "$HERE/generacion_proceso.log" || true
  # MadGraph devuelve 0 aunque falle: se comprueba que la carpeta quedo completa.
  if [ ! -f "$PROCDIR/SubProcesses/procdef_mg5.dat" ] || [ ! -f "$PROCDIR/bin/generate_events" ]; then
    echo
    echo "ERROR: MadGraph no termino de crear $PROCDIR. Mensajes de error de la generacion:"
    grep -n -i -A8 "error\|interrupted\|not executed\|InvalidCmd" "$HERE/generacion_proceso.log" | head -60
    echo "(registro completo en $HERE/generacion_proceso.log)"
    mv "$PROCDIR" "${PROCDIR}_incompleta_$(date +%H%M%S)" 2>/dev/null || true
    exit 1
  fi
  echo "$NJETMAX" > "$PROCDIR/NJETMAX"
else
  if [ ! -f "$PROCDIR/SubProcesses/procdef_mg5.dat" ]; then
    echo "ERROR: $PROCDIR existe pero esta incompleta (una generacion anterior fallo)."
    echo "       Borrarla (rm -rf $PROCDIR) y volver a correr este script."
    exit 1
  fi
  OLD="$(cat "$PROCDIR/NJETMAX" 2>/dev/null || echo '?')"
  if [ "$OLD" != "$NJETMAX" ]; then
    echo "ERROR: $PROCDIR se creo con NJETMAX=$OLD y ahora se pide $NJETMAX."
    echo "       Borrar la carpeta o usar otro PROCDIR."
    exit 1
  fi
fi

# ---------------- 2. Una corrida por masa ----------------
if [ ! -f "$SUMMARY" ]; then
  echo "# m_n1[GeV]  m_x1[GeV]  ctau[mm]  nJetMax  xqcut[GeV]  qCut[GeV]  sigma_LO[pb]  nEventos  LHE" > "$SUMMARY"
fi

for M in $MASSES; do
  RUN="wino_m${M}"
  LHE="$PROCDIR/Events/$RUN/unweighted_events.lhe.gz"
  if [ -f "$LHE" ]; then
    echo ">>> $RUN ya existe, se salta."
    continue
  fi
  rm -rf "$PROCDIR/Events/$RUN"          # restos de una corrida incompleta

  if [ "$XQCUT" = "auto" ]; then
    XQ=$(python3 -c "print('%.1f' % max(10.0, $M/6.0))")
    QC=$(python3 -c "print('%.1f' % max(15.0, $M/4.0))")
  else
    XQ="$XQCUT"
    QC="${QCUT:-$(python3 -c "print('%.1f' % (1.5*$XQCUT))")}"
  fi
  if [ "$NJETMAX" -eq 0 ]; then ICKKW=0; XQ=0; QC=0; else ICKKW=1; fi

  CARD="$CARDS/param_card_${RUN}.dat"
  python3 "$HERE/make_param_card.py" "$TEMPLATE" "$M" "$CARD" --dm "$DM" --ctau "$CTAU"

  CMD="$CARDS/launch_${RUN}.txt"
  {
    if [ -n "$NCORES" ]; then
      echo "set run_mode 2"
      echo "set nb_core $NCORES"
    fi
    echo "launch $PROCDIR --name=$RUN"
    echo "shower=OFF"
    echo "detector=OFF"
    echo "analysis=OFF"
    echo "madspin=OFF"
    echo "reweight=OFF"
    echo "done"
    echo "$CARD"
    echo "set run_card nevents $NEVENTS"
    echo "set run_card ebeam1 $EBEAM"
    echo "set run_card ebeam2 $EBEAM"
    echo "set run_card ickkw $ICKKW"
    echo "set run_card xqcut $XQ"
    echo "set run_card maxjetflavor 4"
    echo "set run_card asrwgtflavor 4"
    echo "set run_card time_of_flight 0"
    echo "set run_card use_syst False"
    echo "set run_card iseed $M"
    echo "done"
  } > "$CMD"

  echo ">>> Corriendo $RUN (xqcut = $XQ GeV, qCut = $QC GeV)"
  mg5 "$CMD"

  [ -f "$LHE" ] || { echo "ERROR: $RUN no produjo $LHE"; exit 1; }
  SIGMA=$(zcat "$LHE" | awk '/Integrated weight \(pb\)/ {s=$NF} END {print s}')
  NEV=$(zcat "$LHE" | grep -c '<event>' || true)
  MX=$(python3 -c "print('%.3f' % ($M + $DM))")
  echo "$M  $MX  $CTAU  $NJETMAX  $XQ  $QC  $SIGMA  $NEV  $LHE" >> "$SUMMARY"
done

echo
echo "Listo. Resumen en $SUMMARY:"
column -t "$SUMMARY" 2>/dev/null || cat "$SUMMARY"
