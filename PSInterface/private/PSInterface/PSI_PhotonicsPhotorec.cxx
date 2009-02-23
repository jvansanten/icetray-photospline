/**
 *@file
 *@brief Implementation of simple photonics Photorec PSInterface 
 *
 *@author Thomas Burgess
 *
 * (c) the IceCube Collaboration
 *
 * $Revision: $
 * $Author: $
 * $Date:  $
 * $Id:  $
 */


//Standard C/C++ includes
#include <iostream>
using std::ostream;

//Photonics includes
extern "C" {
#include "photoamasim.h"
#include "level2.h"    
}

//Local includes
#include "PSInterface/PSI_Logging.h"
#include "PSInterface/PSI_Photonics.h"
#include "PSInterface/PSI_PhotonicsPhotorec.h"
#include "PSInterface/PSI_Coordinate_Photonics.h"

#ifndef PSI_DISABLE_ROOT
//Root includes
#include "TObject.h" //Needed for ClassImp
///Root class implementation macro
ClassImp(PSI_PhotonicsPhotorec);
#endif

//Default constructor
PSI_PhotonicsPhotorec::PSI_PhotonicsPhotorec() :
    PSI_Photonics()
{
    //Set photorec friendly interpolation mode by default (full interpolation)
    SetInterpolationMode(255);
}

//Destructor
PSI_PhotonicsPhotorec::~PSI_PhotonicsPhotorec()
{
}

//Print instance to stream
void PSI_PhotonicsPhotorec::Print(
    ostream &o ) const
{
    PSI_Photonics::Print(o);
    o << "PSI_Photonics";
}
    
//Get photorec request
bool PSI_PhotonicsPhotorec::Get( 
    PSI_Coordinate_Photonics* coordinate,	      
    const double &delay,
    double &amplitude,
    double &probability)
{
    //Set unphysical default values
    amplitude = -1;
    probability = -1;
    
    //Check coordinate
    if ( coordinate == 0 ) {
	log_error("Get failed, no coordinate supplied");
	return false;
    }
    
    if ( !IsTablesLoaded(2) ) {
	log_error("Get failed, no level2 tables loaded");
	return false;
    }
    
    coordinate->CalculateCoordinate();
    bool retval;
    if ( coordinate->IsFinite() ) {   
	retval = photonics_cppio_obj_.get_level2_photorec_finite(      
	    coordinate->GetCoordinateTheta(),
	    coordinate->GetCoordinateRho(),
	    coordinate->GetCoordinatePhi(),
	    coordinate->GetCoordinateL(),
	    coordinate->GetCoordinateStopL(),
	    coordinate->GetCoordinateZSrc(),
	    coordinate->GetCoordinateStopZSrc(),
	    delay,
	    &amplitude,
	    &probability,	    
	    GetInterpolationMode(),
	    0);
    } else {
	retval = photonics_cppio_obj_.get_level2_photorec(
	    coordinate->GetCoordinateTheta(),
	    coordinate->GetCoordinateRho(),
	    coordinate->GetCoordinatePhi(),
	    coordinate->GetCoordinateL(),
	    coordinate->GetCoordinateZSrc(),
	    delay,
	    &amplitude,
	    &probability,
	    GetInterpolationMode(),
	    0);
    }

    //Get energy
    double energy = coordinate->GetTrackEnergy();

    if (retval==false) return retval;

    //Don't scale for 0 energy
    if ( energy == 0 ) { 
	return retval;
    }

    //Never scale with energies below 1
    if ( energy <= 1 ) { 
	energy = 1; 
    }

    //Scale amplitude with  light factor
    amplitude *= light(coordinate->GetTrackType(), energy );
    
    return retval;

}

bool PSI_PhotonicsPhotorec::MeanAmplitude(
    double &amplitude,
    PSI_Coordinate* coordinate)
{
    log_error("MeanAMplitude Not Available in photorec");
    return false;
}

bool PSI_PhotonicsPhotorec::TimeDelay(
    double &timeDelay,
    const double &random,
    PSI_Coordinate* coordinate)
{
    log_error("TimeDelay Not Available in photorec");
    return false;
}


bool PSI_PhotonicsPhotorec::Probability(
    double &probability,
    const double &timeDelay,
    PSI_Coordinate* coordinate)
{ 
    log_error("Probability Not Available in photorec");
    return false;
}
