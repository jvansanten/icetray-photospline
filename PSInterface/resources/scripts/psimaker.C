/**
 *@file
 *@brief Initalizes a new PSInterface instance
 *
 *@author Thomas Burgess
 *
 * (c) the IceCube Collaboration
 *
 * $Revision:$
 * $Author:$
 * $Date:$
 * $Id:$
 */

#ifndef __CINT__
#include <iostream>
using namespace std;
#include "TBenchmark.h"
#include "PSInterface/PSInterface.h"
#include "PSInterface/PSI_Dummy.h"
#include "PSInterface/PSI_PTD.h"
#include "PSInterface/PSI_Photonics.h"
#include "PSInterface/PSI_Coordinate.h"
#include "PSInterface/PSI_Coordinate_PTD.h"
#include "PSInterface/PSI_Coordinate_Photonics.h"
#include "psimaker.h"
#endif

PSInterface* psimaker(
    string interface, 
    string photonicsLevel1TableList, 
    string photonicsLevel2TableList,
    string PTDTable="",
    int PTDMode=0) 
{
    TBenchmark bm;
    bm.Start("psimaker");

    //Pointer to psinterface to create
    PSInterface* psi;  

    if (interface=="Photonics" ) {
	//To use photonics instead
	PSI_Photonics* psi_photonics = new PSI_Photonics();
	if (photonicsLevel1TableList.empty()&&
	    photonicsLevel2TableList.empty()) {
	    cerr << "No tables to load" << endl;
	    return 0;
	}
	if ((!photonicsLevel1TableList.empty()) &&
	    (!psi_photonics->LoadTables(1,photonicsLevel1TableList,1)) ) {
	    cerr << "Failed to load photonics level 1 tables" << endl;
	    return 0;
	}
	if ((!photonicsLevel2TableList.empty()) &&
	    (!psi_photonics->LoadTables(2,photonicsLevel2TableList,1)) ) {
	    cerr << "Failed to load photonics level 2 tables" << endl;
	    return 0;
	}
	psi_photonics->SetInterpolationMode(7);
	psi_photonics->SetPhotonicsVerbosity(1,1);
	psi_photonics->SetPhotonicsVerbosity(2,1);
	psi = psi_photonics;
    } else if ( interface=="PTD" ) {
	PSI_PTD* psi_ptd = new PSI_PTD();
	if ((!PTDTable.empty())&&(!psi_ptd->LoadTable(PTDTable))) {
	    cerr << "Failed to load PTD tables" << endl;
	    return 0;
	}
	if (!psi_ptd->SetPTDMode(PTDMode)) {
	    cerr << "Failed to set PTD mode" << endl;
	    return 0;
	}
	psi = psi_ptd;
    } else if ( interface=="Dummy" ) {
        psi = new PSI_Dummy();
    } else {
	cerr << "Unknown PSInterface " << interface << endl;
	return 0;
    }
    
    bm.Stop("psimaker");
    bm.Show("psimaker");

    return psi;
}
