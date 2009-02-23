/**
 *@file
 *@brief Implementation of basic track-Optical Module coordinate class
 *
 * Thomas Burgess
 * (c) the IceCube Collaboration
 *
 * $Revision: 1.5 $
 * $Author$
 * $Date$
 * $Id$
 */

//Standard C/C++ includes
#include <ostream>
using std::ostream;
#include <cmath>

//Local includes
#include "PSInterface/PSI_Coordinate.h"
#include "PSInterface/PSI_Logging.h"

#ifndef PSI_DISABLE_ROOT
//Root includes
#include "TObject.h" //Needed for ClassImp
//ROOT class implementation macro
ClassImp(PSI_Coordinate);
#endif

//Default constructor
PSI_Coordinate::PSI_Coordinate(
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
#ifndef PSI_DISABLE_ROOT
    TObject(),
#endif
    opticalModuleX_(opticalModuleX),
    opticalModuleY_(opticalModuleY),
    opticalModuleZ_(opticalModuleZ),
    opticalModuleOrientation_(opticalModuleOrientation),
    trackX_(trackX),
    trackY_(trackY),
    trackZ_(trackZ),
    trackTheta_(trackTheta),
    trackPhi_(trackPhi),
    trackLength_(trackLength),
    trackEnergy_(trackEnergy),
    trackType_(trackType),
    thetaRadians_(0),
    phiRadians_(0),
    sinThetaRadians_(0),
    cosThetaRadians_(1),
    sinPhiRadians_(0),
    cosPhiRadians_(1),    
    trackDirectionCalculated_(false),
    trackDirectionX_(0.0),
    trackDirectionY_(0.0),
    trackDirectionZ_(0.0),
    trackToOpticalModuleVectorCalculated_(false),
    trackToOpticalModuleVectorX_(0.0),
    trackToOpticalModuleVectorY_(0.0),
    trackToOpticalModuleVectorZ_(0.0),
    onTrackDistanceCalculated_(false),
    onTrackDistance_(0.0),
    distanceCalculated_(false),
    distance_(0.0),
    rhoCalculated_(false),
    rho_(0.0),
    rhoVectorCalculated_(false),
    rhoVectorX_(0.0),
    rhoVectorY_(0.0),
    rhoVectorZ_(0.0),
    distanceTrackToEmissionCalculated_(false),
    distanceTrackToEmission_(0),
    distanceEmissionToOpticalModuleCalculated_(false),
    distanceEmissionToOpticalModule_(0),
    tGeoCalculated_(false),
    tGeo_(0),
    c_vac_m_ns_(0.2997924),
    n_phase_(1.3195),
    n_group_(1.35634)
{
    //Calculate cherenkov properites
    const double n_phase2 = n_phase_*n_phase_;
    tan_cherenkov_ = sqrt(n_phase2-1.0);
    sin_cherenkov_ = tan_cherenkov_/n_phase_;
    cos_cherenkov_ = 1/n_phase_;
    cherenkov_ = acos(cos_cherenkov_);
    
    //Calculate track angle properties
    thetaRadians_=fabs(trackTheta_*M_PI/180.0-M_PI);
    phiRadians_=trackPhi_*M_PI/180.0-M_PI;
    sinThetaRadians_=sin(thetaRadians_);
    cosThetaRadians_=cos(thetaRadians_);
    sinPhiRadians_=sin(phiRadians_);
    cosPhiRadians_=cos(phiRadians_);
}

//Copy constructor
PSI_Coordinate::PSI_Coordinate( const PSI_Coordinate& coordinate ) 
    :
#ifndef PSI_DISABLE_ROOT
    TObject( coordinate ),
#endif
    opticalModuleX_(coordinate.opticalModuleX_),
    opticalModuleY_(coordinate.opticalModuleY_),
    opticalModuleZ_(coordinate.opticalModuleZ_),
    opticalModuleOrientation_(
	coordinate.opticalModuleOrientation_),
    trackX_(coordinate.trackX_),
    trackY_(coordinate.trackY_),
    trackZ_(coordinate.trackZ_),
    trackTheta_(coordinate.trackTheta_),
    trackPhi_(coordinate.trackPhi_),
    trackLength_(coordinate.trackLength_),
    trackEnergy_(coordinate.trackEnergy_),
    trackType_(coordinate.trackType_),
    thetaRadians_(coordinate.thetaRadians_),
    phiRadians_(coordinate.phiRadians_),
    sinThetaRadians_(coordinate.sinThetaRadians_),
    cosThetaRadians_(coordinate.cosThetaRadians_),
    sinPhiRadians_(coordinate.sinPhiRadians_),
    cosPhiRadians_(coordinate.cosPhiRadians_),
    trackDirectionCalculated_(coordinate.trackDirectionCalculated_),
    trackDirectionX_(coordinate.trackDirectionX_),
    trackDirectionY_(coordinate.trackDirectionY_),
    trackDirectionZ_(coordinate.trackDirectionZ_),    
    trackToOpticalModuleVectorCalculated_(
	coordinate.trackToOpticalModuleVectorCalculated_),
    trackToOpticalModuleVectorX_(coordinate.trackToOpticalModuleVectorX_),
    trackToOpticalModuleVectorY_(coordinate.trackToOpticalModuleVectorY_),
    trackToOpticalModuleVectorZ_(coordinate.trackToOpticalModuleVectorZ_),
    onTrackDistanceCalculated_(coordinate.onTrackDistanceCalculated_),
    onTrackDistance_(coordinate.onTrackDistance_),
    distanceCalculated_(coordinate.distanceCalculated_),
    distance_(coordinate.distance_),
    rhoCalculated_(coordinate.rhoCalculated_),
    rho_(coordinate.rho_),
    rhoVectorCalculated_(coordinate.rhoVectorCalculated_),
    rhoVectorX_(coordinate.rhoVectorX_),
    rhoVectorY_(coordinate.rhoVectorY_),
    rhoVectorZ_(coordinate.rhoVectorZ_),
    distanceTrackToEmissionCalculated_(
	coordinate.distanceTrackToEmissionCalculated_),
    distanceTrackToEmission_(
	coordinate.distanceTrackToEmission_),
    distanceEmissionToOpticalModuleCalculated_(
	coordinate.distanceEmissionToOpticalModuleCalculated_),
    distanceEmissionToOpticalModule_(
	coordinate.distanceEmissionToOpticalModule_),
    tGeoCalculated_(coordinate.tGeoCalculated_),
    tGeo_(coordinate.tGeo_),
    c_vac_m_ns_(coordinate.c_vac_m_ns_),
    n_phase_(coordinate.n_phase_),
    n_group_(coordinate.n_group_),
    tan_cherenkov_(coordinate.tan_cherenkov_),
    sin_cherenkov_(coordinate.sin_cherenkov_),
    cos_cherenkov_(coordinate.cos_cherenkov_),
    cherenkov_(coordinate.cherenkov_)
{
}

//Destructor
PSI_Coordinate::~PSI_Coordinate()
{
}

//Assignment operator
const PSI_Coordinate& PSI_Coordinate::operator=( 
    const PSI_Coordinate& coordinate ) 
{
    //Check for self assignment
    if ( this == &coordinate ) {
        return *this;
    }
    
#ifndef PSI_DISABLE_ROOT
    //Call base class copy constructor
    TObject::operator=( coordinate );
#endif

    //Copy private member variables
    opticalModuleX_=coordinate.opticalModuleX_;
    opticalModuleY_=coordinate.opticalModuleY_;
    opticalModuleZ_=coordinate.opticalModuleZ_;
    opticalModuleOrientation_=coordinate.opticalModuleOrientation_;
    trackX_=coordinate.trackX_;
    trackY_=coordinate.trackY_;
    trackZ_=coordinate.trackZ_;
    trackTheta_=coordinate.trackTheta_;
    trackPhi_=coordinate.trackPhi_;
    trackLength_=coordinate.trackLength_;
    trackEnergy_=coordinate.trackEnergy_;
    trackType_=coordinate.trackType_;
    opticalModuleX_ = coordinate.opticalModuleX_;
    opticalModuleY_ = coordinate.opticalModuleY_;
    opticalModuleZ_ = coordinate.opticalModuleZ_;
    opticalModuleOrientation_ = coordinate.opticalModuleOrientation_;
    trackX_ = coordinate.trackX_;
    trackY_ = coordinate.trackY_;
    trackZ_ = coordinate.trackZ_;
    trackTheta_ = coordinate.trackTheta_;
    trackPhi_ = coordinate.trackPhi_;
    trackLength_ = coordinate.trackLength_;
    trackEnergy_ = coordinate.trackEnergy_;
    trackType_ = coordinate.trackType_;

    thetaRadians_=coordinate.thetaRadians_;
    phiRadians_=coordinate.phiRadians_;
    sinThetaRadians_=coordinate.sinThetaRadians_;
    cosThetaRadians_=coordinate.cosThetaRadians_;
    sinPhiRadians_=coordinate.sinPhiRadians_;
    cosPhiRadians_=coordinate.cosPhiRadians_;

    trackDirectionCalculated_ = coordinate.trackDirectionCalculated_;
    trackDirectionX_ = coordinate.trackDirectionX_;
    trackDirectionY_ = coordinate.trackDirectionY_;
    trackDirectionZ_ = coordinate.trackDirectionZ_;

    trackToOpticalModuleVectorCalculated_ = 
	coordinate.trackToOpticalModuleVectorCalculated_;
    trackToOpticalModuleVectorX_ = coordinate.trackToOpticalModuleVectorX_;
    trackToOpticalModuleVectorY_ = coordinate.trackToOpticalModuleVectorY_;
    trackToOpticalModuleVectorZ_ = coordinate.trackToOpticalModuleVectorZ_;

    onTrackDistanceCalculated_ = coordinate.onTrackDistanceCalculated_;
    onTrackDistance_ = coordinate.onTrackDistance_;

    distanceCalculated_ = coordinate.distanceCalculated_;
    distance_ = coordinate.distance_;

    rhoCalculated_ = coordinate.rhoCalculated_;
    rho_ = coordinate.rho_;

    rhoVectorCalculated_ = coordinate.rhoVectorCalculated_;
    rhoVectorX_ = coordinate.rhoVectorX_;
    rhoVectorY_ = coordinate.rhoVectorY_;
    rhoVectorZ_ = coordinate.rhoVectorZ_;

    distanceTrackToEmissionCalculated_=
	coordinate.distanceTrackToEmissionCalculated_;
    distanceTrackToEmission_=
	coordinate.distanceTrackToEmission_;

    distanceEmissionToOpticalModuleCalculated_=
	coordinate.distanceEmissionToOpticalModuleCalculated_;
    distanceEmissionToOpticalModule_=
	coordinate.distanceEmissionToOpticalModule_;

    tGeoCalculated_=coordinate.tGeoCalculated_;
    tGeo_=coordinate.tGeo_;

    c_vac_m_ns_=coordinate.c_vac_m_ns_;
    n_phase_=coordinate.n_phase_;
    n_group_=coordinate.n_group_;
    tan_cherenkov_=coordinate.tan_cherenkov_;
    sin_cherenkov_=coordinate.sin_cherenkov_;
    cos_cherenkov_=coordinate.cos_cherenkov_;
    cherenkov_=coordinate.cherenkov_;
    opticalModuleX_=coordinate.opticalModuleX_;
    opticalModuleY_=coordinate.opticalModuleY_;
    opticalModuleZ_=coordinate.opticalModuleZ_;
    opticalModuleOrientation_=coordinate.opticalModuleOrientation_;
    trackX_=coordinate.trackX_;
    trackY_=coordinate.trackY_;
    trackZ_=coordinate.trackZ_;
    trackTheta_=coordinate.trackTheta_;
    trackPhi_=coordinate.trackPhi_;
    trackLength_=coordinate.trackLength_;
    trackEnergy_=coordinate.trackEnergy_;
    trackType_=coordinate.trackType_;
    opticalModuleX_ = coordinate.opticalModuleX_;
    opticalModuleY_ = coordinate.opticalModuleY_;
    opticalModuleZ_ = coordinate.opticalModuleZ_;
    opticalModuleOrientation_ = coordinate.opticalModuleOrientation_;
    trackX_ = coordinate.trackX_;
    trackY_ = coordinate.trackY_;
    trackZ_ = coordinate.trackZ_;
    trackTheta_ = coordinate.trackTheta_;
    trackPhi_ = coordinate.trackPhi_;
    trackLength_ = coordinate.trackLength_;
    trackEnergy_ = coordinate.trackEnergy_;
    trackType_ = coordinate.trackType_;

    thetaRadians_=coordinate.thetaRadians_;
    phiRadians_=coordinate.phiRadians_;
    sinThetaRadians_=coordinate.sinThetaRadians_;
    cosThetaRadians_=coordinate.cosThetaRadians_;
    sinPhiRadians_=coordinate.sinPhiRadians_;
    cosPhiRadians_=coordinate.cosPhiRadians_;

    trackDirectionCalculated_ = coordinate.trackDirectionCalculated_;
    trackDirectionX_ = coordinate.trackDirectionX_;
    trackDirectionY_ = coordinate.trackDirectionY_;
    trackDirectionZ_ = coordinate.trackDirectionZ_;

    trackToOpticalModuleVectorCalculated_ = 
	coordinate.trackToOpticalModuleVectorCalculated_;
    trackToOpticalModuleVectorX_ = coordinate.trackToOpticalModuleVectorX_;
    trackToOpticalModuleVectorY_ = coordinate.trackToOpticalModuleVectorY_;
    trackToOpticalModuleVectorZ_ = coordinate.trackToOpticalModuleVectorZ_;

    onTrackDistanceCalculated_ = coordinate.onTrackDistanceCalculated_;
    onTrackDistance_ = coordinate.onTrackDistance_;

    distanceCalculated_ = coordinate.distanceCalculated_;
    distance_ = coordinate.distance_;

    rhoCalculated_ = coordinate.rhoCalculated_;
    rho_ = coordinate.rho_;

    rhoVectorCalculated_ = coordinate.rhoVectorCalculated_;
    rhoVectorX_ = coordinate.rhoVectorX_;
    rhoVectorY_ = coordinate.rhoVectorY_;
    rhoVectorZ_ = coordinate.rhoVectorZ_;

    distanceTrackToEmissionCalculated_=
	coordinate.distanceTrackToEmissionCalculated_;
    distanceTrackToEmission_=
	coordinate.distanceTrackToEmission_;

    distanceEmissionToOpticalModuleCalculated_=
	coordinate.distanceEmissionToOpticalModuleCalculated_;
    distanceEmissionToOpticalModule_=
	coordinate.distanceEmissionToOpticalModule_;

    tGeoCalculated_=coordinate.tGeoCalculated_;
    tGeo_=coordinate.tGeo_;

    c_vac_m_ns_=coordinate.c_vac_m_ns_;
    n_phase_=coordinate.n_phase_;
    n_group_=coordinate.n_group_;
    tan_cherenkov_=coordinate.tan_cherenkov_;
    sin_cherenkov_=coordinate.sin_cherenkov_;
    cos_cherenkov_=coordinate.cos_cherenkov_;
    cherenkov_=coordinate.cherenkov_;

    //Return assigned coordinate
    return *this;
}

void PSI_Coordinate::Clear() {
    opticalModuleX_=0;
    opticalModuleY_=0;
    opticalModuleZ_=0;
    opticalModuleOrientation_=0;
    trackX_=0;
    trackY_=0;
    trackZ_=0;
    trackTheta_=0;
    trackPhi_=0;
    trackLength_=0;
    trackEnergy_=0;
    trackType_=0;
    opticalModuleX_ = 0;
    opticalModuleY_ = 0;
    opticalModuleZ_ = 0;
    opticalModuleOrientation_ = 0;
    trackX_ = 0;
    trackY_ = 0;
    trackZ_ = 0;
    trackTheta_ = 0;
    trackPhi_ = 0;
    trackLength_ = 0;
    trackEnergy_ = 0;
    trackType_ = 0;
    thetaRadians_=0;
    phiRadians_=0;
    sinThetaRadians_=0;
    cosThetaRadians_=1;
    sinPhiRadians_=0;
    cosPhiRadians_=1;
    trackDirectionCalculated_ = false;
    trackDirectionX_ = 0;
    trackDirectionY_ = 0;
    trackDirectionZ_ = 0;
    trackToOpticalModuleVectorCalculated_ = false;
    trackToOpticalModuleVectorX_ = 0;
    trackToOpticalModuleVectorY_ = 0;
    trackToOpticalModuleVectorZ_ = 0;
    onTrackDistanceCalculated_ = false;
    onTrackDistance_ = 0;
    distanceCalculated_ = false;
    distance_ = 0;
    rhoCalculated_ =false;
    rho_ = 0;
    rhoVectorCalculated_=false;
    rhoVectorX_=0;
    rhoVectorY_=0;
    rhoVectorZ_=0;
    distanceTrackToEmissionCalculated_=false;
    distanceTrackToEmission_=0;
    distanceEmissionToOpticalModuleCalculated_= false;
    distanceEmissionToOpticalModule_= 0;
    tGeoCalculated_=false;
    tGeo_=0;
}

void PSI_Coordinate::Set(
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
    //Clear old data
    Clear();

    //Set new data
    opticalModuleX_=opticalModuleX;
    opticalModuleY_=opticalModuleY;
    opticalModuleZ_=opticalModuleZ;
    opticalModuleOrientation_=opticalModuleOrientation;
    trackX_=trackX;
    trackY_=trackY;
    trackZ_=trackZ;
    trackTheta_=trackTheta;
    trackPhi_=trackPhi;
    trackLength_=trackLength;
    trackEnergy_=trackEnergy;
    trackType_=trackType;

    //Calculate new track angle properties
    thetaRadians_=fabs(trackTheta_*M_PI/180.0-M_PI);
    phiRadians_=trackPhi_*M_PI/180.0-M_PI;
    sinThetaRadians_=sin(thetaRadians_);
    cosThetaRadians_=cos(thetaRadians_);
    sinPhiRadians_=sin(phiRadians_);
    cosPhiRadians_=cos(phiRadians_);
}

double PSI_Coordinate::CalcDistance() {
    if ( !distanceCalculated_ ) {
	double ttomx,ttomy,ttomz;
	CalcTrackToOpticalModuleVector(ttomx,ttomy,ttomz);    
	distanceCalculated_ = true;
	distance_ = sqrt(ttomx*ttomx+ttomy*ttomy+ttomz*ttomz);
    }
    return distance_;
}	

double PSI_Coordinate::CalcOnTrackDistance() {    
    if ( !onTrackDistanceCalculated_ ) {
	double ttomx,ttomy,ttomz;
	CalcTrackToOpticalModuleVector(ttomx,ttomy,ttomz);
	double dirx,diry,dirz;
	CalcTrackDirection(dirx,diry,dirz);
	onTrackDistanceCalculated_ = true;
	onTrackDistance_ = ttomx*dirx+ttomy*diry+ttomz*dirz;
    }
    return onTrackDistance_;
}

void PSI_Coordinate::CalcTrackDirection(    
    double &x,
    double &y,
    double &z) 
{
    if ( !trackDirectionCalculated_ ) {
 	//Transform simulation angles to amanda angles
	trackDirectionCalculated_ = true;
	trackDirectionX_=sinThetaRadians_*cosPhiRadians_;
	trackDirectionY_=sinThetaRadians_*sinPhiRadians_;
	trackDirectionZ_=cosThetaRadians_;
    }
    x = trackDirectionX_;
    y = trackDirectionY_;
    z = trackDirectionZ_;    
}

void PSI_Coordinate::CalcTrackToOpticalModuleVector( 
    double &x, double &y, double &z )
{
    if ( !trackToOpticalModuleVectorCalculated_ ) {
	trackToOpticalModuleVectorCalculated_ = true;
	trackToOpticalModuleVectorX_ = opticalModuleX_-trackX_;
	trackToOpticalModuleVectorY_ = opticalModuleY_-trackY_;
	trackToOpticalModuleVectorZ_ = opticalModuleZ_-trackZ_;	
    }
    x = trackToOpticalModuleVectorX_;
    y = trackToOpticalModuleVectorY_;
    z = trackToOpticalModuleVectorZ_;
}

//Get impact vector
void PSI_Coordinate::CalcRhoVector( 
    double &x,
    double &y,
    double &z ) 
{
    if ( !rhoVectorCalculated_ ) {
	double ttomx,ttomy,ttomz;
	CalcTrackToOpticalModuleVector(ttomx,ttomy,ttomz);
	double dirx,diry,dirz;
	CalcTrackDirection(dirx,diry,dirz);
	double onTrackDistance = CalcOnTrackDistance();
	rhoVectorCalculated_ = true;
	rhoVectorX_ = ttomx - onTrackDistance * dirx;
	rhoVectorY_ = ttomy - onTrackDistance * diry;
	rhoVectorZ_ = ttomz - onTrackDistance * dirz;       
    }
    x = rhoVectorX_;
    y = rhoVectorY_;
    z = rhoVectorZ_;
}

double PSI_Coordinate::CalcRho()
{
    if ( !rhoCalculated_ ) {
	double distance = CalcDistance();
	double onTrackDistance = CalcOnTrackDistance();
	rhoCalculated_ = true;
	rho_ = (distance-onTrackDistance)*(distance+onTrackDistance); 
	if (rho_>0) {
	    rho_=sqrt(rho_);
	} else { 
	    rho_=0.0;
	}
    }    
    return rho_;
}

double PSI_Coordinate::CalcDistanceTrackToEmission()
{    
    if ( !distanceTrackToEmissionCalculated_ ) {
	distanceTrackToEmissionCalculated_ = true ;
	//Muons have vertex different from emission point
	if ( trackType_ == 0 ) {
	    distanceTrackToEmission_=
		CalcOnTrackDistance()-CalcRho()/tan_cherenkov_;
	} else {
	    distanceTrackToEmission_=0.0;
	}
    }    
    return distanceTrackToEmission_;
}

double PSI_Coordinate::CalcDistanceEmissionToOpticalModule()
{
    if ( !distanceEmissionToOpticalModuleCalculated_ ) {
	distanceEmissionToOpticalModuleCalculated_=true;
	//Muons and shower have different distance
	if ( trackType_ == 0 ) {
	    distanceEmissionToOpticalModule_=CalcRho()/sin_cherenkov_;
	} else {
	    distanceEmissionToOpticalModule_=CalcDistance();
	}
    }   
    return distanceEmissionToOpticalModule_;
}

double PSI_Coordinate::CalcTGeo( )
{    
    if ( !tGeoCalculated_ ) {
	tGeoCalculated_=true;
	//Muons and shower have different tGeo
	if ( trackType_ == 0 ) {
	    tGeo_ =(CalcDistanceTrackToEmission()+
		    CalcDistanceEmissionToOpticalModule()*n_group_)
		/c_vac_m_ns_;
	} else {
	    tGeo_ = CalcDistance()*n_group_/c_vac_m_ns_;
	}
    }
    return tGeo_;
}

void PSI_Coordinate::MuonEnergyCorrection ( ) 
{
    if ( trackType_ != 0 ) {
	log_warn("Will not correct energy since track is not a muon");
	return;
    }

    if ( trackLength_ <= 0 ) {
	log_warn("Will not correct energy since muon length is <=0");
	return;
    }

    if ( trackEnergy_ <= 0 ) {
	log_warn("Will not correct energy since it is <=0");
	return;
    }

    double L = CalcDistanceTrackToEmission();
   
    if ( L < 0 ) {
	log_debug("Will not correct muon energy since point of"
		  " cherenkov emission occurs before track vertex");	
	return;
    }

    if ( L > trackLength_) {
	log_debug("Distance from vertex to point of cherenkov emission"
		  "exceeds track length, using track length instead");
	L = trackLength_;
    }
    
    // Using dE/dX=-a-b*E => E=(E0+a/b)*e^(-b*x)-a/b
    

    // b is assumed to be constant
    // Value taken from fit in MMC report (hep-ph/0407075) and multiplied
    // with 0.92 for ice density
    const double beff = 3.3e-4;

    // a is calculated from muon length (assuming E=0 when x=trackLength)
    const double aeff=(beff*trackEnergy_)/(exp(trackLength_*beff)-1);

    // fraction 
    const double f=aeff/beff;

    // Now use a and b and solution to dE/dX to calculate corrected 
    // track energy at nominal point of emission (x=L)
    double energy = (trackEnergy_+f)*exp(-beff*L)-f;
    
    if ( energy < 0 ) {
	log_debug("Energy is smaller than 0 GeV and thus invalid, "
		  "setting energy to 0 GeV");
	energy = 0;
    }
    
    log_debug("MuonEnergyCorrection track(Length,E0,E)=(%f,%f,%f)",
	      trackLength_,
	      trackEnergy_,
	      energy);
    
    trackEnergy_ = energy;
}

#ifndef PSI_DISABLE_ICE3

//Scales energy of hadronic shower to fraction of EM shower according to formula by M.Kowalski
void PSI_Coordinate::HadronEnergyCorrection(I3RandomServicePtr random)
{
  double had_frac_f = 1 + pow( (trackEnergy_/0.399), -0.130 )*( 0.467 - 1 );
  double logenergy = log10(trackEnergy_);
  double rms = 0.379 * pow(logenergy, -1.160);
  had_frac_f = random->Gaus(had_frac_f, rms);

  //To not get out of range
  if(had_frac_f < 0) had_frac_f = 0;
  if(had_frac_f > 1) had_frac_f = 1; 

  trackEnergy_ *= had_frac_f;

  //Photonics will scale to 0.802 of energy, so we prepare for that (temp)
  trackEnergy_ /= 0.802;
}

#endif

//Print 
void PSI_Coordinate::Print(
    ostream &o ) const
{
    o << "PSI_Coordinate:\n\tOpticalModule:" 
      << "\n\t\tx=" << opticalModuleX_
      << "\n\t\ty=" << opticalModuleY_
      << "\n\t\tz=" << opticalModuleZ_
      << "\n\t\to=" <<opticalModuleOrientation_
      << "\n\tTrack:"
      << "\n\t\tx=" << trackX_
      << "\n\t\ty=" << trackY_
      << "\n\t\tz=" << trackZ_ 
      << "\n\t\ttheta=" << trackTheta_
      << "\n\t\tphi=" << trackPhi_ 
      << "\n\t\tlength=" << trackLength_
      << "\n\t\tenergy=" << trackEnergy_
      << "\n\t\ttype=" <<  trackType_;
    if ( trackDirectionCalculated_ ) {
	o << "\n\ttrack direction (" 
	  << trackDirectionX_ << "," 
	  << trackDirectionY_ << "," 
	  << trackDirectionZ_ << ")";
    }    
    if ( trackToOpticalModuleVectorCalculated_ ) {
	o << "\n\t Track to OM ( " 
	  << trackToOpticalModuleVectorX_ << ","
	  << trackToOpticalModuleVectorY_ << ","
	  << trackToOpticalModuleVectorZ_ << ")";
    }
    if ( onTrackDistanceCalculated_ ) {
	o << "\n\t Distance Track to OM along track = " <<
	    onTrackDistance_ ;
    }    
    if ( distanceCalculated_ ) { 
	o << "Distance track to OM = " 
	  << distanceCalculated_;
    }
    if ( rhoCalculated_ ) { 
	o << "rho = " << rho_;
    }
    if ( rhoVectorCalculated_ ) {
	o << "rho vector = ("
	  << rhoVectorX_ << ","
	  << rhoVectorY_ << ","
	  << rhoVectorZ_ << ")";
    }
    if ( !distanceTrackToEmissionCalculated_ ) {
	o << "Distance track to emission" 
	  << distanceTrackToEmission_;
    }
    if ( !distanceEmissionToOpticalModuleCalculated_ ) {
	o <<  distanceEmissionToOpticalModule_;
    }
    if ( !tGeoCalculated_ ) {
	o << "tGeo =" << tGeo_;
    }
}

//Printing operator
ostream& operator<< (
    ostream &o, 
    const PSI_Coordinate &coordinate )
{
    coordinate.Print(o);
    return o;
}

//Printing operator for pointer
ostream& operator<< (
    ostream &o,
    const PSI_Coordinate *coordinate )
{
    coordinate->Print(o);
    return o;
}
