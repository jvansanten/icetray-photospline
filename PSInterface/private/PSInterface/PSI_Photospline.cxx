/**
 *@file
 *@brief Implementation of simple photonics PSInterface
 *
 *@author Thomas Burgess
 *
 * (c) the IceCube Collaboration
 *
 * $Revision: 1.36 $
 * $Author: nwhitehorn $
 * $Date: 2009-02-23 10:54:11 -0600 (Mon, 23 Feb 2009) $
 * $Id: PSI_Photospline.cxx 52769 2009-02-23 16:54:11Z nwhitehorn $
 */


//Standard C/C++ includes
#include <iostream>
#include <cmath>
using std::ostream;

//Spline includes
extern "C" {
	#include "bspline.h"
	#include "splinetable.h"
}

//Local includes
#include "PSInterface/PSInterface.h"
#include "PSInterface/PSI_Logging.h"
#include "PSInterface/PSI_Photospline.h"
#include "PSInterface/PSI_Coordinate_Photonics.h" 
#include "PSInterface/PSI_Coordinate.h" 

#ifndef PSI_DISABLE_ROOT
//Root includes
#include "TObject.h" //Needed for ClassImp
///Root class implementation macro
ClassImp(PSI_Photospline);
#endif


//Default constructor
PSI_Photospline::PSI_Photospline() :
    PSInterface()
{
    splinetable.coefficients = NULL;
}

//Destructor
PSI_Photospline::~PSI_Photospline()
{
#warning Leaks memory. Should free() splinetable on exit
}

//Load photonics tables
bool PSI_Photospline::LoadTables( 
    string fileName)
{
    return (readsplinefitstable(fileName.c_str(), &splinetable) == 0);
}

//Clear photonics tables 
bool PSI_Photospline::ClearTables()
{
#warning BROKEN
return false;
}

#if 0
//Get the photonics residual time convention for loadade tables
bool PSI_Photospline::GetRefRefractive(
    const int& level,
    double &ref_ng,
    double &ref_np) 
{
#warning Missing
return false;
}

//Set the photonics residual time convention (ngroup) for loadade tables
bool PSI_Photospline::SetRefRefractive(
    const int& level, 
    const double& ref_ng){
#warning Missing
//Check if photonics tables are loaded
}

bool PSI_Photospline::IsTablesLoaded( 
    const int &level ) const
{
    return (splinetable.coefficients != NULL);
}
#endif

//Print instance to stream
void PSI_Photospline::Print(
    ostream &o ) const
{
    PSInterface::Print(o);
    o << "PSI_Photospline";
}
    
//Creates a PSI_Coordinate
PSI_Coordinate* PSI_Photospline::Coordinate( 
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
    PSI_Coordinate_Photonics* coordinate = new PSI_Coordinate_Photonics(
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
    return coordinate;
}

//Gets a  mean amplitude
bool PSI_Photospline::MeanAmplitude(
    double &amplitude,
    PSI_Coordinate* coordinate)
{
    //Initialize amplitude
    amplitude = -1;

    //Get and check coordinate
    PSI_Coordinate_Photonics* coord = 
	dynamic_cast<PSI_Coordinate_Photonics*>(coordinate);
    if ( coord == 0 ) {
	log_trace("MeanAmplitude failed, "
		  "did not recieve a Photonics coordinate");
	return false;
    }
    
    //Calculate the photonics coordinate
    coord->CalculateCoordinate();
    
    //Get the amplitude

    {
	double x[3];
	int centers[3];

	x[0] = coord->photonics_zsrc_;
	x[1] = coord->photonics_theta_;
	x[2] = coord->photonics_l_;
	
	if (tablesearchcenters(&splinetable, x, centers) != 0)
		amplitude = 0.0;
	else
		amplitude = exp(ndsplineeval(&splinetable, x, centers));
    }
	

#if 0
	amplitude = get_photon_density(
	    coord->photonics_zsrc_,
	    coord->photonics_theta_,
	    type,
	    coord->photonics_l_,
	    coord->photonics_rho_,
	    coord->photonics_phi_*180/M_PI,
	    interpolationMode_,
	    &coord->tableId_,
	    &coord->lBin_,
	    &coord->rhoBin_,
	    &coord->phiBin_);
#endif
    
    if ( coord->trackType_ < 3  && 
	 coord->trackEnergy_ >= 0 ) {
	//Scale amplitude with energy
	double energy = coord->trackEnergy_;

	//If energy is too low, use 1 GeV so light is well behaved
	if ( energy < 1 ) { 
	    log_trace(
		"Energy %f<1 GeV, "
		"using 1 GeV (lowest possible) for light factor "
		"calculation", 
		energy);
	    energy = 1;
	}

	double lightFactor = 1.0;
//	    photonics_cppio_obj_.light(coord->trackType_, energy );

	log_debug("Scaling amplitude %f amplitude with light factor %f",
		  amplitude, lightFactor);
	amplitude = amplitude*lightFactor;
    } else {
	log_trace(
	    "MeanAmplitude, will not scale amplitude with light factor."
	    "Reason: Either track energy is <0 or track type unknown."
	    "track type = %d, track energy = %f",
	    coord->trackType_, coord->trackEnergy_);
    }
    
    log_debug("Got MeanAmplitude=%f from Photonics",amplitude);
    
    return true;
}

//Gets a  Time Delay
bool PSI_Photospline::TimeDelay(
    double &timeDelay,
    const double &random,
    PSI_Coordinate* coordinate)
{
    #warning Missing
    return false;
}

//Gets a  Probability
bool PSI_Photospline::Probability(
    double &probability,
    const double &timeDelay,
    PSI_Coordinate* coordinate)
{
    //Initialize probability
    probability = -1;

    #warning Missing
    return false;
}

