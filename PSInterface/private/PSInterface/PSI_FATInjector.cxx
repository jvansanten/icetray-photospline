/**
 *@file
 *@brief Implementation of Minimal class implementing PSInterface
 *
 *@author Thomas Burgess
 *
 * (c) the IceCube Collaboration
 *
 * $Revision: 1.10 $
 * $Author: burgess $
 * $Date: 2007-05-07 04:55:00 -0400 (Mon, 07 May 2007) $
 * $Id: PSI_FATInjector.cxx 32147 2007-05-07 08:55:00Z burgess $
 */

//C/C++ includes
#include <iostream>
#include <ostream>
#include <cstdlib>
#include <cmath>

#ifndef PSI_DISABLE_ICE3
// GNU Scientific Library includes
#include <gsl/gsl_sf_gamma.h>  // gamma function from special-function package
#endif

//Local includes
#include "PSInterface/PSInterface.h"
#include "PSInterface/PSI_Logging.h"
#include "PSInterface/PSI_Coordinate.h"
#include "PSInterface/PSI_FATInjector.h"

#ifndef PSI_DISABLE_ICE3

#include <dataclasses/I3Units.h>

const double PHOTON_DENSITY = 1e6;
const double START_TIME = -1.5 * I3Units::ns;
const double END_TIME = 1.5 * I3Units::ns;

#else

const double PHOTON_DENSITY = 1e6;
const double START_TIME = -1.5 ;
const double END_TIME = 1.5 ;

#endif


#ifndef PSI_DISABLE_ROOT
//Root includes
#include "TObject.h" //Needed for ClassImp
//PSI includes
ClassImp(PSI_FATInjector);
#endif

//Default constructor
PSI_FATInjector::PSI_FATInjector() : PSInterface()
{
    log_debug("PSI_FATInjector constructor");
}

// Destructor
PSI_FATInjector::~PSI_FATInjector()
{
    log_debug("PSI_FATInjector destructor");
}

// Copy constructor 
PSI_FATInjector::PSI_FATInjector ( const PSI_FATInjector& psinterface ) 
    : PSInterface( psinterface )
{
    log_debug("PSI_FATInjector copy constructor");
}

// Assignment operator 
const PSI_FATInjector& PSI_FATInjector::operator=( 
    const PSI_FATInjector& psinterface )
{
    log_debug("PSI_FATInjector assignment operator");
    
    //Check for self assignment
    if ( this == &psinterface ) {
        return *this;
    }
    
    //Call base class copy constructor
    PSInterface::operator=( psinterface );

    //Return assigned coordinate
    return *this;
}

// Print instance to stream
void PSI_FATInjector::Print(
    ostream &o ) const
{
    log_debug("PSI_FATInjector Print");
}

// Creates a PSI_Coordinate
PSI_Coordinate* PSI_FATInjector::Coordinate( 
    const double &opticalModuleX,
    const double &opticalModuleY,
    const double &opticalModuleZ,
    const double &opticalModuleOrientation,
    const double &trackX,
    const double &trackY,
    const double &trackZ,
    const double &trackTheta,
    const double &trackPhi,
    const double &trackLength,
    const double &trackEnergy,
    const int    &trackType)
{
    log_debug("PSI_FATInjector Coordinate");

    PSI_Coordinate* coord = new PSI_Coordinate(
	opticalModuleX,
	opticalModuleY,
	opticalModuleZ,
	opticalModuleOrientation,
	trackX,
	trackY,
	trackZ,
	trackTheta,
	trackPhi,
	trackLength,
	trackEnergy,
	trackType);
    return coord;
}

// Gets a DUMMY mean amplitude
bool PSI_FATInjector::MeanAmplitude(
    double &amplitude,
    PSI_Coordinate* coordinate)
{
    log_debug("PSI_FATInjector MeanAmplitude");

    //Check coordinate
    if ( coordinate == 0 ) {
	log_error("MeanAmplitude failed, no coordinate supplied");
	return false;
    }

    if(coordinate->GetTrackType()==1){
      //cascades
      amplitude = PHOTON_DENSITY * coordinate->GetTrackEnergy()/ 1000;
    }

    return true;
}

// Gets a DUMMY Time Delay
bool PSI_FATInjector::TimeDelay(
    double &timeDelay,
    const double &random,
    PSI_Coordinate* coordinate)
{
    log_debug("PSI_FATInjector TimeDelay");

    //Check coordinate
    if ( coordinate == 0 ) {
	log_error("TimeDelay failed, no coordinate supplied");
	return false;
    }

    timeDelay = random * (END_TIME - START_TIME) + START_TIME;

    return true;
}

//Gets a dummy Probability
bool PSI_FATInjector::Probability(double &probability,
                            const double &timeDelay,
                            PSI_Coordinate* coordinate)
{
    log_debug("PSI_FATInjector Probability");

    //Check coordinate
    if ( coordinate == 0 ) {
      log_error("Probability failed, no coordinate supplied");
      return false;
    }

    if(timeDelay <= END_TIME && 
       timeDelay >= START_TIME) probability = 1;
    else probability = 0;

    return true;
}

