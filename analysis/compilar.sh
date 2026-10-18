#!/bin/bash
# Compila el analisis. Ajustar las dos rutas (o pasarlas por el entorno):
#   PYTHIA8 = carpeta de Pythia 8.3 ya compilada (contiene include/ y lib/)
#   FASTJET = prefijo de instalacion de FastJet (contiene include/ y lib/)
PYTHIA8=${PYTHIA8:-$HOME/hep/pythia8317}
FASTJET=${FASTJET:-$HOME/hep/fastjet}
g++ -O2 -std=c++17 -W -Wall -Wshadow disappearingT-ATLAS-SR4.cc ToyDetector-ATLAS.cc -o dtSR4 \
  -I$PYTHIA8/include -I$FASTJET/include -L$PYTHIA8/lib -L$FASTJET/lib \
  -Wl,-rpath,$PYTHIA8/lib -Wl,-rpath,$FASTJET/lib -lpythia8 -lfastjet -lz -ldl \
  && echo "OK: ejecutable dtSR4 creado"
