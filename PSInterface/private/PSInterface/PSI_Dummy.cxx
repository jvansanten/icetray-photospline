/**
 *@file
 *@brief Implementation of Minimal class implementing PSInterface
 *
 *@author Thomas Burgess
 *
 * (c) the IceCube Collaboration
 *
 * $Revision: 1.10 $
 * $Author$
 * $Date$
 * $Id$
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
#include "PSInterface/PSI_Dummy.h"


#ifndef PSI_DISABLE_ROOT
//Root includes
#include "TObject.h" //Needed for ClassImp
//PSI includes
ClassImp(PSI_Dummy);
#endif

const double BETA = 1./3.2; //Decay length of PSI_Dummy in units of meters.

//Default constructor
PSI_Dummy::PSI_Dummy() : PSInterface()
{
    log_debug("PSI_Dummy constructor");
}

// Destructor
PSI_Dummy::~PSI_Dummy()
{
    log_debug("PSI_Dummy destructor");
}

// Copy constructor 
PSI_Dummy::PSI_Dummy ( const PSI_Dummy& psinterface ) 
    : PSInterface( psinterface )
{
    log_debug("PSI_Dummy copy constructor");
}

// Assignment operator 
const PSI_Dummy& PSI_Dummy::operator=( 
    const PSI_Dummy& psinterface )
{
    log_debug("PSI_Dummy assignment operator");
    
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
void PSI_Dummy::Print(
    ostream &o ) const
{
    log_debug("PSI_Dummy Print");
}

// Creates a PSI_Coordinate
PSI_Coordinate* PSI_Dummy::Coordinate( 
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
    log_debug("PSI_Dummy Coordinate");

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
bool PSI_Dummy::MeanAmplitude(
    double &amplitude,
    PSI_Coordinate* coordinate)
{
    log_debug("PSI_Dummy MeanAmplitude");

    //Check coordinate
    if ( coordinate == 0 ) {
	log_error("MeanAmplitude failed, no coordinate supplied");
	return false;
    }

    if(coordinate->GetTrackType()==0){
      //tracks
      if((coordinate->CalcDistanceTrackToEmission()>0) &&
	 coordinate->CalcDistanceTrackToEmission() < 
	 coordinate->GetTrackLength() ){
	amplitude = 
	  coordinate->GetTrackEnergy()*exp(-BETA*fabs(coordinate->CalcRho()));
      }else{
	amplitude = 0;
      }
    }else{
      //cascades
      amplitude = 
	coordinate->GetTrackEnergy()*exp(-BETA*coordinate->CalcDistance());
    }

    return true;
}

// Gets a DUMMY Time Delay
bool PSI_Dummy::TimeDelay(
    double &timeDelay,
    const double &random,
    PSI_Coordinate* coordinate)
{
    log_debug("PSI_Dummy TimeDelay");

    //Check coordinate
    if ( coordinate == 0 ) {
	log_error("TimeDelay failed, no coordinate supplied");
	return false;
    }

    double x1=random;
    double x2=random;
    double z=sqrt(-2*log(x1))*sin(2*M_PI*x2);

    if(coordinate->GetTrackType()==0){
      //tracks
      timeDelay=coordinate->CalcRho()*z+coordinate->CalcOnTrackDistance();
    }else{
      //cascades
      timeDelay=coordinate->CalcDistance()*z;
    }

    return true;
}

//Gets a dummy Probability
bool PSI_Dummy::Probability(double &probability,
                            const double &timeDelay,
                            PSI_Coordinate* coordinate)
{
    log_debug("PSI_Dummy Probability");

    //Check coordinate
    if ( coordinate == 0 ) {
      log_error("Probability failed, no coordinate supplied");
      return false;
    }

    // use the pandel function for probability estimate
    // parameters are a mixture between shower and muons
    double lambda = 40.;  // meter
    double tau = 500.; // ns
    double x0 = 100.;  // meter
    double speed = 2.9e-1;  // meter/ns

    // the function is taken from astro-ph/0407044
    // with some mathmatics it's brought to the form
    // (a * b * c) / norm

    // distance in radiation length
    double distance;
    if (coordinate->GetTrackType() == 0) {
      // for muons use calc rho
      distance = coordinate->CalcRho()/lambda;
    } else {
      // for showers, use distance to vertex
      distance = coordinate->CalcDistance()/lambda;
    }

    // inverse time for light travelling one radation length
    double rho = (1.0/tau) + (speed/x0);
    
    double a = pow(rho, distance);
    double b = pow(timeDelay, distance - 1.0);
    double c = exp(-rho * timeDelay);

    double norm=1;
#ifndef PSI_DISABLE_ICE3
    if (distance <= 0 || distance > GSL_SF_GAMMA_XMAX) {
	log_warn("Parameter to gamma function out of range (distance=%f)", distance);
	norm = 1.0;
    } else {
	norm = gsl_sf_gamma(distance);
    }
#endif
    probability = (a * b * c)/norm;
    return true;
}

