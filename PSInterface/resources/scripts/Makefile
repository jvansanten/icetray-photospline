##################################################################
# Makefile
# Builds the Photon Simulation Interface sample scripts
#
# Note, set PHOTONICSINC PHOTONICSLIB and PTDLIB PTDINC
# to library and include paths for Photonics and ptd
# before running make!
#
# Thomas Burgess
# (c) The IceCube Collaboration
#
# $Revision: $
# $Author: $
# $Date:$
# $Id:$
##################################################################

GXX=g++ 
CXXFLAGS+=`root-config --cflags` 
CXXFLAGS+=-O3 -g -DENABLE_PSI_LOGGING -DPSI_LOGGING_LEVEL_INFO -I../../public
CXXFLAGS+=-I${PTDINC} -I${PHOTONICSINC}
LDFLAGS+=`root-config --libs` 
LDFLAGS+= -L../.. -lPSInterface
LDFLAGS+=-L${PTDLIB} -lptd
LDFLAGS+=-L${PHOTONICSLIB} -lphotonics -llevel2amasim -lphotoamasim

all: simple directhits cowboyTest simplePh

directhits: directhits.o psimaker.o
	${CXX} ${LDFLAGS} $^ -o $@

simple: simple.o psimaker.o
	${CXX} ${LDFLAGS} $^ -o $@

cowboyTest: cowboyTest.o psimaker.o
	${CXX} ${LDFLAGS} $^ -o $@

simplePh: simplePh.o  psimaker.o
	${CXX} ${LDFLAGS} $^ -o $@

# Compile a C++ file
%.o:%.C
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: clean

clean:
	rm -rf simple cowboyTest directhits \
	       cowboyTest.o simple.o directhits.o psimaker.o
