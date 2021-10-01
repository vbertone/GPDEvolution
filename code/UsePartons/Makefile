CXX = clang++

CXXFLAGS += -O3 -fPIC# -std=c++11

# LHAPDF
LHAPDFINCS = $(shell lhapdf-config --cppflags)
LHAPDFLIBS = $(shell lhapdf-config --ldflags)

# APFEL++
APFELPPINCS = $(shell apfelxx-config --cppflags)
APFELPPLIBS = $(shell apfelxx-config --ldflags)

# Boost
BOOSTINCS = -I/usr/local/include/
BOOSTLIBS = -L/usr/local/lib/ -lboost_system

# Partons
PARTONSINCS = -I/usr/local/include/
PARTONSLIBS = -L/usr/local/lib/ -lElementaryUtils -lNumA++ -lPARTONS

# Now set up the compiler and link flags and libs
CXXFLAGS += $(LHAPDFINCS) $(APFELPPINCS) $(BOOSTINCS) $(PARTONSINCS)
LDFLAGS  += $(LHAPDFINCS) $(APFELPPINCS) $(BOOSTINCS) $(PARTONSINCS)

CLIBS += $(LHAPDFLIBS) $(APFELPPLIBS) $(BOOSTLIBS) $(PARTONSLIBS)

install : all
all : GPDEvolution

GPDEvolution: GPDEvolution.o 
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS)

.SUFFIXES : .cc .o .f .c

.cxx.o:	 
	$(CXX) $(CXXFLAGS) -c $< 

.f.o:	 
	$(F77)  -c $< 

clean:
	rm -rf *.lo *.o *.la GPDEvolution *~
