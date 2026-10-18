#!/bin/bash
# Corre dtSR4 sobre todas las masas listadas en grid_summary.dat.
# Uso:  ./run_analysis.sh            (variables opcionales abajo)
set -euo pipefail
SUMMARY="${SUMMARY:-grid_summary.dat}"                 # lo escribe run_grid.sh
MAP="${MAP:-mapa_ATLAS137fb_Fig2a_4capas.dat}"         # mapa de eficiencia
OUTDIR="${OUTDIR:-resultados}"
IDLLP="${IDLLP:-1000024}"                              # particula cargada de vida larga

[ -x ./dtSR4 ]   || { echo "ERROR: falta ./dtSR4 (correr ./compilar.sh)"; exit 1; }
[ -f "$SUMMARY" ] || { echo "ERROR: falta $SUMMARY"; exit 1; }
[ -f "$MAP" ]     || { echo "ERROR: falta $MAP"; exit 1; }
mkdir -p "$OUTDIR"

grep -v '^#' "$SUMMARY" | while read -r M MX CTAU NJ XQ QC SIG NEV LHE; do
  OUT="$OUTDIR/out_m${M}.dat"
  if [ -s "$OUT" ]; then echo ">>> m = $M GeV ya analizado, se salta."; continue; fi
  [ -f "$LHE" ] || { echo "ERROR: no existe $LHE"; exit 1; }
  echo ">>> Analizando m = $M GeV (nJetMax = $NJ, qCut = $QC GeV)"
  nice -n 19 ./dtSR4 "$LHE" "$MAP" "$OUT" "$IDLLP" "$NJ" "$QC" > "$OUTDIR/log_m${M}.txt" 2>&1 < /dev/null \
    || { echo "ERROR: fallo dtSR4 para m = $M (ver $OUTDIR/log_m${M}.txt)"; exit 1; }
done
echo "Listo. Tablas de eficiencia en $OUTDIR/"
