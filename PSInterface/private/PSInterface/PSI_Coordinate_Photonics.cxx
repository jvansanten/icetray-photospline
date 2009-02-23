/**
 *@file
 *@brief Implementation of Photonics track-OM coordinate class
 *
 * Thomas Burgess & Daan Hubert
 * (c) the IceCube Collaboration
 *
 * $Revision $
 * $Author $
 * $Date $
 * $Id $
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
#include "PSInterface/PSI_Coordinate_Photonics.h"

#ifndef PSI_DISABLE_ROOT
//Root includes
#include "TObject.h" //Needed for ClassImp
///ROOT class implementation macro
ClassImp( PSI_Coordinate_Photonics );
#endif

PSI_Coordinate_Photonics::PSI_Coordinate_Photonics(
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
    tableSetId_(-1),
    tableId_(-1),
    lBin_(-1),
    phiBin_(-1),
    rhoBin_(-1),
    stopTableId_(-1),
    stopLBin_(-1),
    stopRhoBin_(-1),
    stopPhiBin_(-1),
    photonics_zsrc_(-1),
    photonics_theta_(-1),
    photonics_l_(-1),
    photonics_rho_(-1),
    photonics_phi_(-1),    
    photonics_dist_(-1),
    photonics_stop_zsrc_(-1),
    photonics_stop_l_(-1),
    photonics_stop_dist_(-1),
    zOrigin_(0),
    finite_(false),
    calculated_(false)
{
}

PSI_Coordinate_Photonics::PSI_Coordinate_Photonics(
    const PSI_Coordinate_Photonics& coordinate ) 
    : 
    PSI_Coordinate( coordinate ),
    tableSetId_(coordinate.tableSetId_),
    tableId_(coordinate.tableId_),
    lBin_(coordinate.lBin_),
    phiBin_(coordinate.phiBin_),
    rhoBin_(coordinate.rhoBin_),
    stopTableId_(coordinate.stopTableId_),
    stopLBin_(coordinate.stopLBin_),
    stopRhoBin_(coordinate.stopRhoBin_),
    stopPhiBin_(coordinate.stopPhiBin_),
    photonics_zsrc_(coordinate.photonics_zsrc_),
    photonics_theta_(coordinate.photonics_theta_),
    photonics_l_(coordinate.photonics_l_),
    photonics_rho_(coordinate.photonics_rho_),
    photonics_phi_(coordinate.photonics_phi_),
    photonics_dist_(coordinate.photonics_dist_),
    photonics_stop_zsrc_(coordinate.photonics_stop_zsrc_),
    photonics_stop_l_(coordinate.photonics_stop_l_),
    photonics_stop_dist_(coordinate.photonics_stop_dist_),
    zOrigin_(coordinate.zOrigin_),
    finite_(coordinate.finite_),
    calculated_(coordinate.calculated_)
{

}

PSI_Coordinate_Photonics::~PSI_Coordinate_Photonics()
{    

}

const PSI_Coordinate_Photonics& PSI_Coordinate_Photonics::operator=( 
    const PSI_Coordinate_Photonics& coordinate ) 
{
    //Check for self assignment
    if ( this == &coordinate ) {
        return *this;
    }
	
    //Call base class assignment operator
    PSI_Coordinate::operator=( coordinate );

    //Copy private memebers
    tableSetId_ = coordinate.tableSetId_;
    tableId_ = coordinate.tableId_;
    lBin_ = coordinate.lBin_;
    phiBin_ = coordinate.phiBin_;
    rhoBin_ = coordinate.rhoBin_;
    stopTableId_ = coordinate.stopTableId_;
    stopLBin_ = coordinate.stopLBin_;
    stopRhoBin_ = coordinate.stopRhoBin_;
    stopPhiBin_ = coordinate.stopPhiBin_;
    photonics_zsrc_ = coordinate.photonics_zsrc_;
    photonics_theta_ = coordinate.photonics_theta_;
    photonics_l_ = coordinate.photonics_l_;
    photonics_rho_ = coordinate.photonics_rho_;
    photonics_phi_ = coordinate.photonics_phi_;
    photonics_dist_ = coordinate.photonics_dist_;
    photonics_stop_zsrc_= coordinate.photonics_stop_zsrc_;
    photonics_stop_l_ = coordinate.photonics_stop_l_;
    photonics_stop_dist_ = coordinate.photonics_stop_dist_;
    zOrigin_ = coordinate.zOrigin_;
    finite_ = coordinate.finite_;
    calculated_ = coordinate.calculated_;

    //Return assigned coordinate
    return *this;
}

void PSI_Coordinate_Photonics::Clear()
{
    //Calll base class clear first
    PSI_Coordinate::Clear();
    tableSetId_ = -1;
    tableId_ = -1;
    lBin_ = -1;
    phiBin_ = -1;
    rhoBin_ = -1;
    stopTableId_ = -1;
    stopLBin_ = -1;
    stopRhoBin_ = -1;
    stopPhiBin_ = -1;
    photonics_zsrc_ = -1;
    photonics_theta_ = -1;
    photonics_l_ = -1;
    photonics_rho_ = -1;
    photonics_phi_ = -1;    
    photonics_dist_ = -1;
    photonics_stop_zsrc_ = -1;
    photonics_stop_l_ = -1;
    photonics_stop_dist_ = -1;
    zOrigin_ = 0;
    finite_ = false;
    calculated_ = false;
}

void PSI_Coordinate_Photonics::CalculatePPerp( 
    double & x, double &y, double &z )
{
  x = -cosPhiRadians_*cosThetaRadians_;
  y = -sinPhiRadians_*cosThetaRadians_;
  z = sinThetaRadians_ ;
}

void PSI_Coordinate_Photonics::CalculateCoordinate( )
{
    if ( calculated_ ) {
	return;
    }
    calculated_=true;

    //Set finite state
    if ( ( trackType_== 0 ) && ( trackLength_ > 0 ) ) {
	//Finite muon
	finite_=true;
    } else {
	//Shower or infinite muon
	finite_=false;
    }

    //Calculate x-axis in track coordinate system
    double pperpx,pperpy,pperpz;
    CalculatePPerp(pperpx,pperpy,pperpz);

    //Get impact parameter (rho)
    double rho = CalcRho();

    //Calculate cosine of azimuth angle (phi in photonics) 
    //(angle between pperp and rho)
    float cosazi = 1.0;
    if ( rho > 0 ) {
	//Get impact vector (rho-vector)
	double rhox=0, rhoy=0, rhoz=0;
	CalcRhoVector(rhox,rhoy,rhoz);
	cosazi=(rhox*pperpx+rhoy*pperpy+rhoz*pperpz)/rho;
	//Restrict cosine to -1,1
	if ( cosazi > 1.0 ) { 
	    cosazi =  1.0;
	} else if ( cosazi < -1.0 ) {
	    cosazi = -1.0;
	}
    }
    photonics_zsrc_ = zOrigin_ + trackZ_;
    photonics_theta_ = trackTheta_;
    photonics_l_ = CalcOnTrackDistance();
    photonics_rho_ = rho;
    photonics_phi_ = acos(cosazi);
    photonics_dist_ = CalcDistance();
    
    //For finte coordinates also calculation stop information
    if ( finite_ ) {
	//track direction
	double dirx=0, diry=0, dirz=0;
	CalcTrackDirection(dirx,diry,dirz);
	
	//Calculate track end coordinate
 	float tendx=trackX_+dirx*trackLength_;
 	float tendy=trackY_+diry*trackLength_;
 	float tendz=trackZ_+dirz*trackLength_;
	//track end to OM vector
	double ttomendx=opticalModuleX_-tendx;
	double ttomendy=opticalModuleY_-tendy;
	double ttomendz=opticalModuleZ_-tendz;
	photonics_stop_zsrc_=zOrigin_+tendz;

	photonics_stop_l_=ttomendx*dirx+ttomendy*diry+ttomendz*dirz;
	photonics_stop_dist_= 
	    sqrt(ttomendx*ttomendx+ttomendy*ttomendy+ttomendz*ttomendz);
    }
}


//Print coordinate to stream
void PSI_Coordinate_Photonics::Print( 
    ostream &o ) const    
{    
    PSI_Coordinate::Print(o);
    o << "\n\tPSI_Coordinate_Photonics:"
      << "\n\t\tz origin=" << zOrigin_
      << "\n\t\tCoordinate:"
      << "\n\t\t\tzSrc=" << photonics_zsrc_
      << "\n\t\t\ttheta=" << photonics_theta_ 
      << "\n\t\t\tl=" << photonics_l_ 
      << "\n\t\t\trho=" << photonics_rho_ 
      << "\n\t\t\tphi=" << photonics_phi_ 
      << "\n\t\t\tdist=" << photonics_dist_
      << "\n\t\t\tstop zSrc=" << photonics_stop_zsrc_
      << "\n\t\t\tstop l=" << photonics_stop_l_ 
      << "\n\t\t\tstop dist=" << photonics_stop_dist_
      << "\n\t\tBin information:"
      << "\n\t\t\tTable set id=" << tableSetId_
      << "\n\t\t\tTable id=" << tableId_ 
      << "\n\t\t\tl bin=" << lBin_
      << "\n\t\t\tphi bin=" << phiBin_
      << "\n\t\t\trho bin=" << rhoBin_
      << "\n\t\t\tstop table id=" << stopTableId_ 
      << "\n\t\t\tstop l bin=" << stopLBin_ 
      << "\n\t\t\tstop rho bin=" << stopRhoBin_
      << "\n\t\t\tstop phi bin=" << stopPhiBin_;    
    if ( finite_ ) o << "\n\tCoordinate is finite";
    if ( calculated_ ) o << "\n\tCoordinate is calculated";
}
