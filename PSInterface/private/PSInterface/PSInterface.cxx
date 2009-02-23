/**
 *@file
 *@brief Implementation of PSInterface - Photon Simulation Interface
 *
 * Thomas Burgess
 * (c) the IceCube Collaboration
 *
 * $Revision: 1.31 $
 * $Author$
 * $Date$
 * $Id$
 */

//C/C++ includes
#include <iostream>
#include <ostream>
using std::ostream;
#include <cmath>

//Local includes
#include "PSInterface/PSInterface.h"
#include "PSInterface/PSI_Logging.h"

#ifndef PSI_DISABLE_ROOT
//Root includes
#include "TObject.h" //Needed for ClassImp
///Root class implementation macro
ClassImp(PSInterface);
#endif

//Default constructor
PSInterface::PSInterface()
{
}

//Destructor
PSInterface::~PSInterface()
{
}

//Copy constructor 
PSInterface::PSInterface ( 
    const PSInterface& psinterface )
{
}

const PSInterface& PSInterface::operator=( 
    const PSInterface& psinterface )
{
    //Check for self assignment
    if ( this == &psinterface ) {
	return *this;
    }
    
    //Return assigned psinterface
    return *this;
}
    
//Printing operator for this and all derived classes
ostream& operator<< (
    ostream &o,
    const PSInterface &psi )
{
    psi.Print(o);
    return o;    
}

//Printing operator for this and all derived classes
ostream& operator<< (
    ostream &o,
    const PSInterface *psi )
{    
    psi->Print(o);
    return o;
}

//Create a coordinate
PSI_Coordinate* PSInterface::MakeCoordinate( 
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
    PSI_Coordinate* coordinate = Coordinate(
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

    if ( coordinate == 0) {
	log_fatal("MakeCoordinate failed to create coordinate");
    }

    return coordinate;
}

//Get a mean amplitude
bool PSInterface::GetMeanAmplitude(
    double &amplitude,
    PSI_Coordinate* coordinate)
{
    if ( coordinate == 0 ) {
	log_fatal("GetMeanAmplitude failed becase it was called with"
		  " null coordinate");
    }
    
    if (! MeanAmplitude(amplitude, coordinate) ) {
	log_debug("MeanAmplitude failed.");
	return false;
    }
    return true;
}

//Get hit time delay
bool PSInterface::GetTimeDelay(
    double &timeDelay,
    const double &random,
    PSI_Coordinate* coordinate)
{
    if ( coordinate == 0 ) {
	log_fatal("GetTimeDelay failed becase it was called with"
		  " null coordinate");
    }

    if ( ( random < 0 ) || ( random > 1) ) {
	log_error("Invalid random number %f in GetTimeDelay. Should be [0,1]", 
		 random);
	return false;
    }

    if (! TimeDelay(timeDelay, random, coordinate) ) {
	log_error("TimeDelay failed. (timeDelay,random)=%f,%f", 
		 timeDelay, random);
	return false;
    }

    return true;
}

//Get a hit probability
bool PSInterface::GetProbability(
    double &probability,
    const double &timeDelay,
    PSI_Coordinate* coordinate)
{
    if ( coordinate == 0 ) {
	log_fatal("GetProbability failed becase it was called with"
		  " null coordinate");
    }

    if ( timeDelay < 0 ){
	log_error("Invalid time delay %f<0 in GetProbability.", 
		 timeDelay);
	return false;
    }

    if ( ! Probability(probability, timeDelay, coordinate) ) {
	log_error("GetProbability failed. (probability,timeDelay)=%f,%f", 
		 probability,timeDelay);
	return false;
    }
    
    return true;
}

//Return a Poisson Random Mumber
int PSInterface::PoissonRandomNumber( 
    const double &lambda,
    const double &random ) const
{
    //Check random number
    if ( ( random < 0 ) || ( random >= 1 ) ) {
        log_error("PoissonRandomNumber failed, bad random number\n");
	return -1;
    }
    
    //Initialize poisson
    double P = exp( -lambda );          //probability Po(X=1)
    double sum = P;                     //cumulant
    
    //Check if P allready exceeds random
    if ( sum >= random ) {
        return 0;
    }
    
    //Poisson outcome
    unsigned int k = 0;
    
    //Loop over all outcomes
    for ( k=1; k < maxPoisson_; ++k ) {
        P *= lambda / (double) k;       //Find next probability
        sum += P;                       //Increase cumulant
        if ( sum >= random ) break;     //Leave loop when done
    }

    //return outcome
    return k;
}

//Print instance to stream
void PSInterface::Print(
    ostream &o ) const
{
    o << "PSInterface";
}

bool PSInterface::SetAngularSelection(int const& level, float const& low, float const& high)
{
    log_debug("Setting level %d angular selection to %f<angle<%f",
	      level,low,high);
    // store angular selection so it can be queried by modules using service
    angularSelectLow_ = low;
    angularSelectHigh_ = high;

    return true;    
}
