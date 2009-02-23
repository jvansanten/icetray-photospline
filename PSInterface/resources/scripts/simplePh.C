/**
 *@file
 *@brief Simple compilable timing example for PSI_PhotonicsPhotorec
 *
 *@author Johan Lundberg
 * (c) the IceCube Collaboration
 */


#ifndef __CINT__
#include <iostream>
using namespace std;
#include "TRandom.h"
#include "TBenchmark.h"

#include "PSInterface/PSInterface.h"
#include "PSInterface/PSI_Photonics.h"
#include "PSInterface/PSI_Coordinate.h"
#include "PSInterface/PSI_Coordinate_Photonics.h"
#include "PSInterface/PSI_PhotonicsPhotorec.h"
#include "PSInterface/PSInterface.h"
#include "PSInterface/PSI_Coordinate.h"
#include "psimaker.h"
#endif
 

/**
 * PSI_PhotonicsPhotorec simple test
 *
 * Run with ".x simplePhotorec.C" in root or  "root -q -l -b simplePhotorec.C"
 *  from commandline
 *
 *@author Thomas Burgess
 */
void simplePhotorec() 
{
#ifdef __CINT__
    //Load the psi library if using root interactive
    if ( gSystem->Load("libPSInterface.so") < 0 ) {
	return -1;
    }   
    gROOT->ProcessLine(".L psimaker.C");
#endif

    //Initialize
    PSI_PhotonicsPhotorec rec;
    if ( !rec.LoadTables(2,"level2_muon_photorec.list") ) {
	cerr << "Failed to load photonics tables, exiting" << endl;
	return;
    }
    
    //Create a photonics coordinate for a 
    //downlooking optical module at (0,0,0) and a
    //100 GeV muon track at (0,0,-10) with angles (45,180) for 10m
    double OM_x=0, OM_y=0, OM_z=0, OM_o=-1;
    double TR_x=0, TR_y=0, TR_z=-10, TR_th=45, TR_ph=180,
	TR_length=10, TR_energy=100;
    int TR_type = 0;

    cout << "performing photorec requests "<< endl;

    long ii=0;
    long jj=0;

    double amplitude = -1;
    double probability = -1;

    double outdum = 0;

    PSI_Coordinate_Photonics *coordinate = 
	new PSI_Coordinate_Photonics(
	    OM_x,OM_y,OM_z,OM_o,
	    1,1,TR_z,TR_th,TR_ph,
	    TR_length,TR_energy,TR_type);

    //Check coordinate
    if ( coordinate == 0 ) {
	cerr << "Failed to make coordinate" << endl;
	return;
    }

    double outdum1=0;
    double outdum2=0;
    double outdum3=0;
    double outdum4=0;
    double outdum5=0;
    double outdum6=0;
    double outdum7=0;
    double outdum8=0;
    
    for (double q=1;q<600;q++){
	for (double j=1;j<200;j++){
	    for (double timedelay=-20;timedelay<2000;timedelay+=30){
		ii++;

		//Define amplitude and probability to fetch
		if (amplitude<0)
		    jj++;
    
		//Make photorec query and print result
		if ( !rec.Get(coordinate,timedelay,amplitude,probability) ) {
		    cerr << "Failed to get photorec result" << endl;
		    return;
		}

		outdum+=amplitude+probability;
		outdum1=(amplitude>outdum1)?amplitude:outdum1;
		outdum2=(probability>outdum2)?probability:outdum2;
		outdum3=
		    (amplitude*probability>outdum3)?
		    amplitude*probability:outdum3;
		outdum4+=outdum1;
		outdum5+=outdum2;
		outdum6+=outdum3;
		outdum7+=amplitude;
		outdum8+=probability;

	    }
	}
    }
    delete coordinate;

    cout << "number of request: " << ii
	 << "\ndummy values:" 
	 << outdum << " " 
	 << outdum1<< " " 
	 << outdum2<< " " 
	 << outdum3<< " " 
	 << outdum4<< " " 
	 << outdum5<< " " 
	 << outdum6<< " " 
	 << outdum7<< " " 
	 << outdum8	 
	 << "\nLast Amplitude="<<amplitude<<" pdf="<<probability
	 << "\nresult table result fraction, " <<  (ii-jj)*100.0/ii 
	 << "% " << endl;
}



//#ifndef CINT
int main() { simplePhotorec(); return 0; }
//#endif
