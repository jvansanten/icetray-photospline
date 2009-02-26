/**
 *@file
 *@brief Definition of basic Track-Optical Module coordinate class
 *
 *@author Thomas Burgess
 *
 * (c) the IceCube Collaboration
 *
 * $Revision: 1.3 $
 * $Author$
 * $Date$
 * $Id$
 */
 
#ifndef __PSI_Coordinate_h__
#define __PSI_Coordinate_h__

//Standard C/C++ includes
#include <cmath>
#include <ostream>
#include <string>
using std::string;
using std::ostream;

#ifndef PSI_DISABLE_ROOT
//Root includes
#include "TObject.h" //Needed for ClassDef
#endif

//Local includes
#include "PSInterface/PSI_Logging.h"

#ifndef PSI_DISABLE_ICE3
#include "phys-services/I3RandomService.h"
#endif

/**
 *@brief A track-OM coordinate class
 *
 * This is a base class for coordinates used in  
 * for example photonics or PTD PSInterface 
 * implementations. 
 *
 * Basic properties are calculated on construction and
 * can be accessed by GetFunctions. 
 *                                                         
 *@author Thomas Burgess                                   
 */
class PSI_Coordinate
#ifndef PSI_DISABLE_ROOT
: public TObject 
#endif
{

#ifndef PSI_DISABLE_ROOT
    ///ROOT class definition macro
    ClassDef(PSI_Coordinate, 2); ///A track-OM coordinate class
#endif

public:
    ///Name log for icetray    
    SET_LOGGER("PSI_Coordinate");

    ///Allow PSI_Coordinate_Photonics to access coordinate private members
    friend class PSI_Coordinate_Photonics;

    ///Allow PSI_Photonics to access coordinate private members
    friend class PSI_Photonics;
    friend class PSI_Photospline;

    ///Allow PSI_Coordinate_PTD to access coordinate private members
    friend class PSI_Coordinate_PTD;
        
    /**
     *@brief Default Constructor
     *
     *@param opticalModuleX            OM location X coordinate (m)
     *@param opticalModuleY            OM location Y coordinate (m)
     *@param opticalModuleZ            OM location Z coordinate (m)
     *@param opticalModuleOrientation  OM orientation (1.0 up, -1.0 down)
     *@param trackX                    Track vertex X coordinate (m)
     *@param trackY                    Track vertex Y coordinate (m)
     *@param trackZ                    Track vertex Z coordinate (m)
     *@param trackTheta                Track zenith angle (degrees)
     *@param trackPhi                  Track azimuthal angle (degrees)
     *@param trackLength               Track length (m, negative=infinite)
     *@param trackEnergy               Track energy (GeV)
     *@param trackType                 Track type (mu:0,e+e-:1,hadr:2)
     */
    PSI_Coordinate(
	    const double &opticalModuleX = 0,
	    const double &opticalModuleY = 0,
	    const double &opticalModuleZ = 0,
	    const double &opticalModuleOrientation = 0,
	    const double &trackX = 0,
	    const double &trackY = 0,
	    const double &trackZ = 0,
	    const double &trackTheta = 0,
	    const double &trackPhi = 0,
	    const double &trackLength = 0,
	    const double &trackEnergy = 0,
	    const int    &trackType = 0);

    /**
     *@brief Copy constructor
     *
     *@param coordinate coordinate to copy from
     */
    PSI_Coordinate( const PSI_Coordinate& coordinate );

    /**
     *@brief Virtual Destructor
     */
    virtual ~PSI_Coordinate();    

    /**
     *@brief Assignment operator
     *
     *@param coordinate coordinate to assign from
     *@return Assigned coordinate
     */
    const PSI_Coordinate& operator=( const PSI_Coordinate& coordinate );
    
    /**
     *@brief Printing operator for this and all derived clases
     *
     * Don't overload, to change behaviour in deriving classes
     * override Print which iscalled in this method.
     *
     *@param o   Stream to write to
     *@param coordinate PSInterface instance to print
     *
     *@result The written stream
     */
    friend std::ostream& operator<< (
	std::ostream &o, 
	const PSI_Coordinate &coordinate );

    /**
     *@brief Printing operator for this and all derived clases
     *
     *@param o   Stream to write to
     *@param coordinate PSInterface instance to print
     *
     *@result The written stream
     */
    friend std::ostream& operator<< (
	std::ostream &o, 
	const PSI_Coordinate *coordinate);

    /**
     *@brief Clears all values in the coordinate
     */
    virtual void Clear();

    /**
     *@param opticalModuleX            OM location X coordinate (m)
     *@param opticalModuleY            OM location Y coordinate (m)
     *@param opticalModuleZ            OM location Z coordinate (m)
     *@param opticalModuleOrientation  OM orientation (1.0 up, -1.0 down)
     *@param trackX                    Track vertex X coordinate (m)
     *@param trackY                    Track vertex Y coordinate (m)
     *@param trackZ                    Track vertex Z coordinate (m)
     *@param trackTheta                Track zenith angle (degrees)
     *@param trackPhi                  Track azimuthal angle (degrees)
     *@param trackLength               Track length (m, negative=infinite)
     *@param trackEnergy               Track energy (GeV)
     *@param trackType                 Track type (mu:0,e+e-:1,hadr:2)
     */
    virtual void Set(
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
	const int    &trackType);
      
    /**
     *@brief Calculate track direction vector
     *
     *@param x    direction x component result (m)
     *@param y    direction y component result (m)
     *@param z    direction z component result (m)
     */
    void CalcTrackDirection(    
	double &x, double &y, double &z);

    /**
     *@brief Calculate distance between OM location and track vertex
     *@return distance (m)
     */
    double CalcDistance();
    
    /**
     *@brief Calculate distance between track and OM along the track
     *
     * Calculates the distance from track vertex to optical module 
     * location projected along the track. In other words length along
     * track from vertex to point of closest approach. The result will
     * be 0 for track starting at the point of closest approach,
     * positive if the track starts before the OM and negative if
     * the track starts after the OM.
     *
     * This is the photonics coordinate l
     *
     *@return on-track distance (m)
     */
    double CalcOnTrackDistance();

    /**
     *@brief Calculate rho
     *
     * Rho is the impact parameter - ie the distance between
     * the optical module location and to the closest point on
     * track.
     *
     *@return rho (m)
     */
    double CalcRho();

    /**
     *@brief Calculate rho vector
     *
     *@param x   direction x component (m)
     *@param y   direction y component (m)
     *@param z   direction z component (m)
     */
    void CalcRhoVector(
	double &x, double &y, double &z);

    /**
     *@brief Calculate vector from track to optical module
     *
     *@param x    direction x component (m)
     *@param y    direction y component (m)
     *@param z    direction z component (m)
     */
    void CalcTrackToOpticalModuleVector( 
	double &x, double &y, double &z);

    /**
     *@brief Calculate distance from track vertex to cherenkov emission point
     *
     *@return distance (m)
     */
    double CalcDistanceTrackToEmission();

    /**
     *@brief Calculate distance point of cherenkov emission and OM
     *@return distance (m)
     */
    double CalcDistanceEmissionToOpticalModule();
    
    /**
     *@brief Calculate light time from track vertex to OM hit
     *@return time (ns)
     */
    double CalcTGeo();
    
    /**
     *@brief Get Optical module X coordinate
     */
    double GetOpticalModuleX() const { 
	return opticalModuleX_; }

    /**
     *@brief Get Optical module y coordinate
     */
    double GetOpticalModuleY() const { 
	return opticalModuleY_; }

    /**
     *@brief Get Optical module Z coordinate
     */
    double GetOpticalModuleZ() const { 
	return opticalModuleZ_; }

    /**
     *@brief Get Optical module orientation
     */
    double GetOpticalModuleOrientation() const { 
	return opticalModuleOrientation_; }
    
    /**
     *@brief Get Track vertex X coordinate
     */
    double GetTrackX() const { 
	return trackX_; }
    
    /**
     *@brief Get Track vertex Y coordinate
     */
    double GetTrackY() const { 
	return trackY_; }
    
    /**
     *@brief Get Track vertex Z coordinate
     */
    double GetTrackZ() const { 
	return trackZ_; }
    
    /**
     *@brief Get Track zenith angle 
     */
    double GetTrackTheta() const { 
	return trackTheta_; }
    
    /**
     *@brief Get Track azimuth angle
     */
    double GetTrackPhi() const { 
	return trackPhi_; }
    
    /**
     *@brief Get Track length
     */
    double GetTrackLength() const { 
	return trackLength_; }

    /**
     *@brief Get Track energy
     */
    double GetTrackEnergy() const { 
	return trackEnergy_; }
    
    /**
     *@brief Get Track type
     */
    int GetTrackType() const { 
	return trackType_; }

    /**
     *@brief Correct muon energy for energy losses
     */
    void MuonEnergyCorrection();

#ifndef PSI_DISABLE_ICE3
    /**
     *@brief Correct hadron energy
     */
    void HadronEnergyCorrection(I3RandomServicePtr);
#endif

protected:
    /**
     *@brief Print coordinate to stream
     *
     *@param o      stream to print to
     */
    virtual void Print( 
	ostream &o ) const;    

private:    

    /**
     *@brief Optical module X coordinate (m)
     */
    double opticalModuleX_;

    /**
     *@brief Optical module y coordinate (m)
     */
    double opticalModuleY_;

    /**
     *@brief Optical module Z coordinate (m)
     */
    double opticalModuleZ_;

    /**
     *@brief Optical module orientation (1.0 up, -1.0 down)
     */
    double opticalModuleOrientation_;
    
    /**
     *@brief Track vertex X coordinate (m)
     */
    double trackX_;
    
    /**
     *@brief Track vertex Y coordinate (m)
     */
    double trackY_;
    
    /**
     *@brief Track vertex Z coordinate (m)
     */
    double trackZ_;
    
    /**
     *@brief Track zenith angle (degrees)
     */
    double trackTheta_;
    
    /**
     *@brief Track azimuth angle (degrees)
     */
    double trackPhi_;
    
    /**
     *@brief Track length (m, negative=infinite)
     */
    double trackLength_;

    /**
     *@brief Track energy (GeV)
     */
    double trackEnergy_;
    
    /**
     *@brief Track type (mu:0,e+e-:1,hadr:2) 
     */
    int    trackType_;

    /**
     *@brief Theta in photon simulation system
     */
    double thetaRadians_;

    /**
     *@brief Phi in photon simulation system
     */
    double phiRadians_;

    /**
     *@brief Sin of theta in photon simulation system
     */
    double sinThetaRadians_;

    /**
     *@brief Cos of theta in photon simulation system
     */
    double cosThetaRadians_;

    /**
     *@brief Sin of phi in photon simulation system
     */
    double sinPhiRadians_;

    /**
     *@brief Cos of phi in photon simulation system
     */
    double cosPhiRadians_;
 
    /**
     *@brief True if track direction has been calculated
     */
    bool trackDirectionCalculated_;
    
    /**
     *@brief Track direction X-component
     */
    double trackDirectionX_;

    /**
     *@brief Track direction Y-component
     */
    double trackDirectionY_;

    /**
     *@brief Track direction Z-component
     */
    double trackDirectionZ_;

    /**
     *@brief True if track to OM vector is calculated
     */
    bool trackToOpticalModuleVectorCalculated_;

    /**
     *@brief Track to OM vector X-component
     */
    double trackToOpticalModuleVectorX_;

    /**
     *@brief Track to OM vector Y-component
     */
    double trackToOpticalModuleVectorY_;
    
    /**
     *@brief Track to OM vector Z-component
     */
    double trackToOpticalModuleVectorZ_;

    /**
     *@brief True if on track distance is calculated
     */
    bool onTrackDistanceCalculated_;
    
    /**
     *@brief On track distance
     */
    double onTrackDistance_;

    /**
     *@brief True if distance is calculated
     */
    bool distanceCalculated_;

    /**
     *@brief Distance track to OM
     */
    double distance_;

    /**
     *@brief True if rho is calculated
     */
    bool rhoCalculated_;

    /**
     *@brief Rho - impact parameter (m)
     */
    double rho_;

    /** 
     *@brief True if rho vector is calculated
     */
    bool rhoVectorCalculated_;

    /**
     *@brief Rho vector Z-component
     */
    double rhoVectorX_;

    /**
     *@brief Rho vector Z-component
     */
    double rhoVectorY_;

    /**
     *@brief Rho vector Z-component
     */
    double rhoVectorZ_;

    /** 
     *@brief True if distance from track vertex to emission is calculated
     */
    bool distanceTrackToEmissionCalculated_;

    /**
     *@brief Distance from track vertex to emission 
     */
    double distanceTrackToEmission_;

    /** 
     *@brief True if is calculated
     */
    bool distanceEmissionToOpticalModuleCalculated_;

    /**
     *@brief Distance from emission to OM
     */
    double distanceEmissionToOpticalModule_;

    /** 
     *@brief True if geometrical time is calculated
     */
    bool tGeoCalculated_;    
   
    /**
     *@brief Geometrical time
     */
    double tGeo_;

    /**
     *@brief speed of light in m/ns
     */
    double c_vac_m_ns_;

    /**
     *@brief phase index of refraction
     */
    double n_phase_;

    /**
     *@brief group index of refraction
     */
    double n_group_;

    /**
     *@brief tan of Cherenkov angle
     */
    double tan_cherenkov_;

    /**
     *@brief sin of Cherenkov angle
     */
    double sin_cherenkov_;

    /**
     *@brief cos of Cherenkov angle
     */
    double cos_cherenkov_;

    /**
     *@brief Chrenkov angle
     */
    double cherenkov_;
};

#endif
