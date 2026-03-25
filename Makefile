# ================================================
# Makefile for Pythia 8.310 + FastJet + Matching + DT
# ================================================

PYTHIA8   = /user/f/fhpinto/RPVLLP/HEPsoftware/pythia8310
FASTJET3  = /user/f/fhpinto/RPVLLP/HEPsoftware/fastjet-install

PYTHIALIB = $(PYTHIA8)/lib
PYTHIAINC = $(PYTHIA8)/include
FJLIB     = $(FASTJET3)/lib
FJINC     = $(FASTJET3)/include

OBJECTS = disappearingT-ATLAS.o ToyDetector-ATLAS.o Matching.o trackletEfficiency.o
INCLUDE = -I$(PYTHIAINC) -I$(PYTHIAINC)/Pythia8 -I$(PYTHIAINC)/Pythia8Plugins -I$(FJINC) -I./

CXXFLAGS = -O2 -std=c++11 $(INCLUDE)
LDFLAGS  = -L$(PYTHIALIB) -Wl,-rpath,$(PYTHIALIB) -L$(FJLIB)
LDLIBS   = -lpythia8 -lfastjet

# Target principal
disappearingT-ATLAS: $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) $(LDFLAGS) $(LDLIBS)

# Regla genérica para compilar .cc -> .o
%.o: %.cc
	$(CXX) -c $< $(CXXFLAGS) -o $@

# Limpieza
clean:
	rm -f $(OBJECTS) disappearingT-ATLAS
