CXX = clang++

CXXFLAGS += -O3 -fPIC -std=c++11

# LHAPDF
LHAPDFINCS = #$(shell lhapdf-config --cppflags)
LHAPDFLIBS = #$(shell lhapdf-config --ldflags)

# APFEL++
APFELPPINCS = #$(shell apfelxx-config --cppflags)
APFELPPLIBS = #$(shell apfelxx-config --ldflags)

# ElementaryUtile
PARTONSINCS = -I/usr/local/include
PARTONSLIBS = -L/usr/local/lib/ -lElementaryUtils -L/usr/local/lib/ -lPARTONS

# Now set up the compiler and link flags and libs
CXXFLAGS += $(LHAPDFINCS) $(APFELPPINCS) $(PARTONSINCS)
LDFLAGS  += $(LHAPDFINCS) $(APFELPPINCS) $(PARTONSINCS)

CLIBS += $(LHAPDFLIBS) $(APFELPPLIBS) $(PARTONSLIBS)

install : all
all : Evolution

Evolution: Evolution.o 
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS)

.SUFFIXES : .cc .o .f .c .cpp

.cxx.o:	 
	$(CXX) $(CXXFLAGS) -c $< 

.f.o:	 
	$(F77)  -c $< 

clean:
	rm -rf *.lo *.o *.la Evolution *~
