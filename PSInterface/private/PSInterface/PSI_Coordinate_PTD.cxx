/**
 *@file
 *@brief Implementation of PTD track-OM coordinate class
 *
 *@author Thomas Burgess
 *
 * (c) the IceCube Collaboration
 *
 * $Revision: 1.4 $
 * $Author$
 * $Date$
 * $Id$
 */

//
// Inludes, namespaces, predeclarations
//

//Standard C/C++ includes
#include <ostream>
using std::ostream;
#include <string>
using std::string;
#include <cmath>

//Local includes
#include "PSInterface/PSI_Logging.h"
#include "PSInterface/PSI_Coordinate.h"
#include "PSInterface/PSI_Coordinate_PTD.h"
//PTD includes
extern "C" {
#include "ptdread.h"
}

#ifndef PSI_DISABLE_ROOT
//Root includes
#include "TObject.h" //Needed for ClassImp
///ROOT class implementation macro
ClassImp( PSI_Coordinate_PTD );
#endif


//Default Constructor
PSI_Coordinate_PTD::PSI_Coordinate_PTD(	
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
    const int    &trackType,
    const double &zOrigin) 
    :
    PSI_Coordinate(
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
	trackType),
    zOrigin_(zOrigin)
{
    //Reset complete photonVals struct
    photonVals_.rho=0;
    photonVals_.zet=0;
    photonVals_.theta=0;
    photonVals_.phi=0;
    photonVals_.time=0;
    photonVals_.amp=0;
    photonVals_.random=0;
    photonVals_.dummy=0;
    photonVals_.iRho=0;
    photonVals_.iZ=0;
    photonVals_.iTheta=0;
    photonVals_.iPhi=0;
    photonVals_.iT=0;
    photonVals_.dummy_zsize=0;
    photonVals_.dummy_zlow=0;    
}

//Copy constructor
PSI_Coordinate_PTD::PSI_Coordinate_PTD( 
    const PSI_Coordinate_PTD& coordinate ) :
    PSI_Coordinate( coordinate ),
    photonVals_(coordinate.photonVals_),
    zOrigin_(coordinate.zOrigin_)
{

}

//Virtual Destructor
PSI_Coordinate_PTD::~PSI_Coordinate_PTD()
{

}

//Assignment operator
const PSI_Coordinate_PTD& PSI_Coordinate_PTD::operator=( 
    const PSI_Coordinate_PTD& coordinate )
{
    //Check for self assignment
    if ( this == &coordinate ) {
        return *this;
    }
    
    //Call base class copy constructor
    PSI_Coordinate::operator=( coordinate );

    //Copy private member variables
    photonVals_=coordinate.photonVals_;
    zOrigin_=coordinate.zOrigin_;

    //Return assigned coordinate
    return *this;
}

bool PSI_Coordinate_PTD::CalculateCoordinate( 
    const simParams &simparams )
{
    //Get and check rho
    const double rho = CalcRho();
    if ( ( rho >  simparams.rh_max ) || ( rho <  simparams.rh_min ) ) {
	log_warn(
	    "Rho %f was outside range allowable by table [%f,%f]. ", 
	    rho, simparams.rh_max, simparams.rh_min);
	return false;
    }

    //Convert on track distance to z_vertex and check result
    double z_vertex = CalcOnTrackDistance();
    //Check if z is allowed by table
    if ( ( z_vertex > simparams.ra_max ) || 
	 ( z_vertex < simparams.ra_min ) ) {
	log_warn(
	    "z_vertex %f was outside range allowable by table [%f,%f]. ", 
	    z_vertex, simparams.ra_max, simparams.ra_min);
	return false;
    }

    //Length from vertex to cherenkov emission point
    double d_vtx_emit= CalcDistanceTrackToEmission();

    //Get track direction
    double dirx, diry, dirz;
    CalcTrackDirection(dirx,diry,dirz);

    //Get track to optical module vector
    double ttomx, ttomy, ttomz;
    CalcTrackToOpticalModuleVector(ttomx,ttomy,ttomz);

    //Calculate vector OM to cherenkov emission point
    double dx_pmtemit[3] = {
 	-d_vtx_emit*dirx-ttomx,
 	-d_vtx_emit*diry-ttomy,
	-d_vtx_emit*dirz-ttomz};
    
    //Length of dx_pmtemit
    double d_pmtemit= rho/sin_cherenkov_;

    //Calculate cos theta     
    double costheta = 0.0;     
    if ( d_pmtemit != 0.0 ) {
	costheta = 
	    dx_pmtemit[2]*GetOpticalModuleOrientation()/d_pmtemit; 
    }

    //Calculate theta
    double theta = 0.0;
    if (costheta >= 1.0) {
	theta =0;
    } else if (costheta <= -1.0) {
	theta = acos(-1.0);
    } else {
	theta = acos(costheta);
    }

    //Check if theta is allowed by table
    if ( ( theta >  simparams.th_max ) || ( theta <  simparams.th_min ) ) {
	log_warn(
	    "Theta %f was outside range allowable by table [%f,%f]. ", 
	    theta, simparams.th_max, simparams.th_min);
	return false;
    }
    
    //phivec_pmt is the x-product (dx_pmtemit)x(pmtdirection)
    double phivec_pmt[3] = {	 	
	dx_pmtemit[1]*GetOpticalModuleOrientation(),
	-dx_pmtemit[0]*GetOpticalModuleOrientation(),
	0};
    double d_phi_pmt = sqrt( phivec_pmt[0]*phivec_pmt[0] + 
			     phivec_pmt[1]*phivec_pmt[1]);

    //phivec_mu is the x-product (dx_pmtemit) x (track direction). 
    //It has length (d_pmtemit*sinth_cherenkov).  
    //Do not calculate phivec_mu(3) since it isn't used
    //due to some mix up of coordinate system between detector and amanda
    //we have to throw in a -1. for all mu_dir exept z           
    double phivec_mu[3] = { 
	+dx_pmtemit[1]*dirz-dx_pmtemit[2]*diry,
	+dx_pmtemit[2]*dirx-dx_pmtemit[0]*dirz,
	0};

    double d_phi_mu = d_pmtemit*sin_cherenkov_;

    double cosphi = 0;
    if ( ( d_phi_pmt>1E-8) && (d_phi_mu> 1E-8) ) {
	cosphi=(phivec_pmt[0]*phivec_mu[0]+
		phivec_pmt[1]*phivec_mu[1] ) / (d_phi_pmt*d_phi_mu);
    }
    
    double phi = 0;
    if ( cosphi>= 1. ) {
	phi= acos( 1.);
    } else if ( cosphi <= -1. ) {
	phi= acos( -1. );
    } else {
	phi= acos( cosphi );
    }

    //Correct phi so it fits the tables
    if ( ( phi > simparams.ph_max ) && 
	 ( simparams.ph_max/M_PI < 1.5 ) ){
	phi=2*M_PI-phi;
    }
    
    //Check if phi is allowed by table
    if ( ( phi >  simparams.ph_max ) || ( phi <  simparams.ph_min ) ) {
	log_warn(
	    "Phi %f was outside range allowable by table [%f,%f]. ", 
	    phi, simparams.ph_max, simparams.ph_min);
	return false;
    }

    //Assign values
    photonVals_.rho = rho;
    photonVals_.zet = z_vertex;
    photonVals_.theta = theta;
    photonVals_.phi = M_PI-phi;
    
    return true;
}

//Print coordinate to stream
void PSI_Coordinate_PTD::Print( 
    ostream &o ) const    
{    
    PSI_Coordinate::Print(o);
    o << "\n\tPSI_Coordinate_PTD"
	///Print everything in photonVals_ struct
      << "\n\t\trho=" << photonVals_.rho
      << "\n\t\tzet=" << photonVals_.zet
      << "\n\t\ttheta= " << photonVals_.theta
      << "\n\t\tphi=" << photonVals_.phi
      << "\n\t\ttime=" << photonVals_.time
      << "\n\t\tamp=" << photonVals_.amp
      << "\n\t\trandom=" << photonVals_.random
      << "\n\t\tdummy=" << photonVals_.dummy
      << "\n\t\tiRho=" << photonVals_.iRho
      << "\n\t\tiZ=" << photonVals_.iZ
      << "\n\t\tiTheta=" << photonVals_.iTheta
      << "\n\t\tiPhi=" << photonVals_.iPhi
      << "\n\t\tiT=" << photonVals_.iT
      << "\n\t\tidummy_zsize=" << photonVals_.dummy_zsize
      << "\n\t\tdummy_zlow=" << photonVals_.dummy_zlow;
}
