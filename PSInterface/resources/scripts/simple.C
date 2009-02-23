/**
 *@file
 *@brief Simple example for PSInterface
 *
 * Thomas Burgess
 * (c) the IceCube Collaboration
 *
 * $Revision: $
 * $Author$
 * $Date$
 * $Id$
 */

#ifndef __CINT__
#include <iostream>
using namespace std;
#include "TRandom.h"
#include "TBenchmark.h"
#include "PSInterface/PSInterface.h"
#include "PSInterface/PSI_Coordinate.h"
#include "psimaker.h"
#endif

/**
 * PSInterface simple test
 *
 * Run with ".x simple.C" in root or  "root -q -l -b simple.C"
 *  from commandline
 * Remeber to set LD_LIBRARY_PATH so root can find
 * - libPSInterface.so
 * - liblevel2amasim.so
 * - libphotoamasim.so
 * - libphotonics.so
 * (or copy the libraries(or softlink) to current directory)
 *
 * This macro can also be compiled if needed
 *<pre>
 * g++ simple.C -I../../public/	\
 * `root-config --libs --cflags` -L../.. -lPSInterface \
 * -L${PTDLIB} -lptd -I${PTDINC} -I{PHOTONICSINC} \ 
 * -L${PHOTONICSLIB} -llevel2amasim -lphotoamasim -lphotonics \
 * -DENABLE_PSI_LOGGING -o simple
 *</pre>
 *
 *@author Thomas Burgess
 */
int simple()
{
#ifdef __CINT__
    //Load the psi library if using root interactive
    if ( gSystem->Load("libPSInterface.so") < 0 ) {
	return -1;
    }
    gROOT->ProcessLine(".L psimaker.C");
#endif

    PSInterface* psi = 0;
    if ( ( psi = psimaker(
	       "PTD", 
	       "level1_table.list",
	       "level2_table.list",
	       "ptd-tables/tables01_bulk_mam/31ce5.9s21h2.pb.125.3302_011.nta",
	       500 ) ) == 0 ) {
	cerr << "Somehow no PSI was created, exiting" << endl;
	return 0;
    }

    //Downlooking optical module at (0,0,0)
    double OM_x=0, OM_y=0, OM_z=0, OM_o=-1;

    //100 GeV muon track at (0,0,-10) with angles (45,180) for 10m
    double TR_x=0, TR_y=0, TR_z=-10, TR_th=45, TR_ph=180,
	TR_length=10, TR_energy=100;
    int TR_type = 0;

    //Create a coordinate for the track and om
    PSI_Coordinate* coord = 0 ;
    if ( ( coord = psi->MakeCoordinate(
	       OM_x,OM_y,OM_z,OM_o,
	       TR_x,TR_y,TR_z,TR_th,TR_ph,
	       TR_length,TR_energy,TR_type))==0) {
	cout << "Failed to create coordinate" << endl;
	return -1;
    }

    //Get the mean amplitude
    double amplitude=-1;
    if (!psi->GetMeanAmplitude(amplitude,coord)) {
	cerr << "Failed to get mean ampitude" << endl; 
	return -1;
    } 
    cout << "Amplitude=" << amplitude <<endl;

    //Use the amplitude to simulate the number of hits
    cout << "Number of hits = " 
	 << psi->PoissonRandomNumber(amplitude,gRandom->Rndm()) << endl;

    //Get time delay for one hit
    double dt=-1;
    double random = 0.24; // or use gRandom->Rndm()
    if (!psi->GetTimeDelay(dt,random,coord)) {
	cerr << "Failed to get time delay" << endl;
	return -1;
    }
    cout << "Time delay=" << dt <<endl;

    //Get hit probability from a random time delay (0-100ns)
    double probability=-1;
    double timeDelay = 54.4131; //Or use gRandom->Rndm()*100
    if (!psi->GetProbability(probability,timeDelay,coord)) {
	cerr << "Failed to get probability" << endl;
	return -1;
    }
    cout << "Probability=" << probability <<endl; 

    delete coord;
    delete psi;

    return 0;
}

#ifndef __CINT__
int main() { return simple(); }
#endif
