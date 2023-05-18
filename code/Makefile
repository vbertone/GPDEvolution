CXX = clang++
CXXFLAGS += -O3 -fPIC -std=c++17

# LHAPDF
LHAPDFINCS = $(shell lhapdf-config --cppflags)
LHAPDFLIBS = $(shell lhapdf-config --ldflags)

# APFEL++
APFELPPINCS = $(shell apfelxx-config --cppflags)
APFELPPLIBS = $(shell apfelxx-config --ldflags)

# Boost
BOOSTINCS = -I/usr/local/include/
BOOSTLIBS = -L/usr/local/include/ -lboost_system

# Now set up the compiler and link flags and libs
CXXFLAGS += $(LHAPDFINCS) $(APFELPPINCS) $(BOOSTINCS)
LDFLAGS  += $(LHAPDFINCS) $(APFELPPINCS) $(BOOSTINCS)

CLIBS += $(LHAPDFLIBS) $(APFELPPLIBS) $(BOOSTLIBS)

install : all
all : Convolutions GPDEvolution GPDEvolutionPol GPDEvolutionTrans GPDEvolutionGhost GPDEvolutionERBL GPDEvolutionPion GPDEvolutionDer GPDEvolutionDterm GPDEvolutionOperator

Convolutions: Convolutions.o 
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS)

GPDEvolution: GPDEvolution.o 
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS)

GPDEvolutionPol: GPDEvolutionPol.o 
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS)

GPDEvolutionTrans: GPDEvolutionTrans.o 
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS)

GPDEvolutionGhost: GPDEvolutionGhost.o 
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS)

GPDEvolutionERBL: GPDEvolutionERBL.o 
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS)

GPDEvolutionPion: GPDEvolutionPion.o 
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS)

GPDEvolutionDer: GPDEvolutionDer.o 
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS)

GPDEvolutionDterm: GPDEvolutionDterm.o 
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS)

GPDEvolutionOperator: GPDEvolutionOperator.o 
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS)

.SUFFIXES : .cc .o .f .c

.cxx.o:	 
	$(CXX) $(CXXFLAGS) -c $< 

.f.o:	 
	$(F77)  -c $< 

clean:
	rm -rf *.lo *.o *.la Convolutions GPDEvolution GPDEvolutionPol GPDEvolutionTrans GPDEvolutionGhost GPDEvolutionERBL GPDEvolutionPion GPDEvolutionDer GPDEvolutionDterm GPDEvolutionOperator *~
