/**
 *@file
 *@brief Definition of Photonics track-OM coordinate class
 *
 *@author Thomas Burgess
 *
 * (c) the IceCube Collaboration
 *
 * $Revision $
 * $Author $
 * $Date $
 * $Id $
 */

#ifndef __PSI_Coordinate_Photonics_h__
#define __PSI_Coordinate_Photonics_h__

//
// Inludes, namespaces, predeclarations
//

//Standard C/C++ includes
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
#include "PSInterface/PSI_Coordinate.h"
class PSI_Photonics;

/**
 *@brief A track-OM coordinate class for Photonics
 *
 * This is a coordinate that can calculate and hold Photonics specific 
 * information.It is used in PSI_Photonics 
 * 
 * A photonics coordinates is (zSrc,theta,l,phi,rho,dist) where z is
 * the depth, theta the track angle, l the length along the track from
 * vertex to emission, rho the impact parameter and dist the distance between
 * the track to the om. (zSrc,theta) are used to select which table to read.
 * (l,phi,rho) are used to aquire data from within a table. The coordinate 
 * information can be accessed by the corresponding GetCoordinate function 
 * for each coordinate.
 *
 * For the combination of (zSrc,theta) the is a table id, and for each 
 * of (l,phi,rho) there is a bin id.
 *
 * For finite muons there exists two sets of (zSrc,l,rho,dist) 
 * coordinates and bins, one releative to the vertex and one relative to the 
 * track end. Coordinates can be access by the GetStop function for each 
 * coordinate.
 *
 * All coordinates and bins are initialized to -1 before on construction.
 * To calculate the coordinates call CalculateCoordinate(). The bin 
 * information is filled when calling PSInterface::GetMeanAmplitude(). 
 *
 *@author Thomas Burgess                                   
 */
class PSI_Coordinate_Photonics : public PSI_Coordinate {

#ifndef PSI_DISABLE_ROOT
    ///ROOT class definition macro
    ClassDef(PSI_Coordinate_Photonics, 2); ///A track-OM coordinate class
#endif

public:
    ///Name log for icetray    
    SET_LOGGER("PSI_Coordinate_Photonics");

    ///Allow PSI_Photonics to access coordinate private members
    friend class PSI_Photonics;
    
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
    PSI_Coordinate_Photonics(	
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
    PSI_Coordinate_Photonics( 
	const PSI_Coordinate_Photonics& coordinate );

    /**
     *@brief Virtual Destructor
     */
    virtual ~PSI_Coordinate_Photonics();    

    /**
     *@brief Assignment operator
     *
     *@param coordinate coordinate to assign from
     *
     *@return Assigned coordinate
     */
    const PSI_Coordinate_Photonics& operator=( 
	const PSI_Coordinate_Photonics& coordinate );

    /**
     *@brief Clear photonics coordinate data
     */
    virtual void Clear();

    /**
     *@brief Calculate track-system x-axis
     *
     * This is perpendicular to the track direction
     * and defines the track phi axis
     *
     *@param x calculated pperp vector x-coordinate
     *@param y calculated pperp vector y-coordinate
     *@param z calculated pperp vector z-coordinate
     */
    void CalculatePPerp(double &x, double &y, double &z);

    /**
     *@brief Calculate Photonics coordinate
     */
    void CalculateCoordinate();

    /**
     *@brief Get Depth of the light source in amanda coordinate system
     *@return zSrc
     */
    float GetCoordinateZSrc() const { return photonics_zsrc_; }
    
    /**
     *@brief Get theta - zenith angle of the track (radians)
     *@return Theta
     */
    float GetCoordinateTheta() const { return photonics_theta_; }
    
    /**
     *@brief Get distance from source to emitter along the track
     *
     * (how much before (+) or behind (-) the emission
     * point the OM is located)
     *
     *@return L
     */
    float GetCoordinateL() const { return photonics_l_; }
    
    /**
     *@brief Get rho - perpendicular distance track to OM (Impact parameter)
     *@return Rho
     */
    float GetCoordinateRho() const { return photonics_rho_; }
 
    /**
     *@brief Get phi - angle between x-axis rotated into muon coordinate
     * system and rho or dist
     *@return  phi
     */
    float GetCoordinatePhi() const { return photonics_phi_; }
    
    /**
     *@brief Get distance track to OM
     *@return dist
     */
    float GetCoordinateDist() const { return photonics_dist_; }
 
    /**
     *@brief Get ZSrc for track end (Only meaningful for finite muons)
     *@see GetCoordinateZSrc()
     *@return ZSrc
     */
    float GetCoordinateStopZSrc() const { return photonics_stop_zsrc_; }
    
    /**
     *@brief Get L for track end (Only meaningful for finite muons)
     *@see GetCoordinateL()
     *@return L
     */
    float GetCoordinateStopL() const { return photonics_stop_l_; }
 
    /**
     *@brief Get dist for track end (Only meaningful for finite muons)
     *@see GetCoordinateDist()
     *@return dist
     */
    float GetCoordinateStopDist() const { return photonics_stop_dist_; }

    /**
     *@brief Get Z offset from origin
     *@return z origin
     */
    float GetZOrigin() const { return zOrigin_ ; }

    /**
     *@brief Get finite state
     *@return true if muon is finite, false otherwise
     */
    bool IsFinite() const { return finite_ ; }

    /**
     *@brief Set Z offset from origin
     *@param zOrigin z origin
     */
    void SetZOrigin(const float &zOrigin ) { zOrigin_=zOrigin; }

    /**
     *@brief Set finite state
     *@param finite true if muon is finite, false otherwise
     */
    void SetFinite( const bool &finite ) { finite_=finite; }

protected:
    /**
     *@brief Print coordinate to stream
     *
     *@param o      stream to print to
     */
    virtual void Print( 
	ostream &o ) const;    

private:    
    /**@brief photonics table set id*/
    int tableSetId_;

    /**@brief photonics table id*/
    int tableId_;

    /**@brief photonics table length bin id @see GetCoordinateL()*/
    int lBin_;

    /**@brief photonics table phi bin id @see GetCoordinatePhi()*/
    int phiBin_;

    /**@brief photonics table rho bin id @see GetCoordinateRho()*/
    int rhoBin_;

    /**
     *@brief photonics table id for track end
     * (Only meaningful for finite muons)
     */
    int stopTableId_;

    /**
     *@brief photonics table length bin id for track end @see GetCoordinateL()
     * (Only meaningful for finite muons)
     */
    int stopLBin_;

    /**
     *@brief photonics table rho bin id for track end @see GetCoordinateRho()
     * (Only meaningful for finite muons)
     */
    int stopRhoBin_;

    /**
     *@brief photonics table phi bin if for track end @see GetCoordinateZSrc()
     *( Only meaningful for finite muons)
     */
     int stopPhiBin_;

    /**
     *@brief zSrc
     *@see GetCoordinateZSrc()
     */
    float photonics_zsrc_;

    /**
     *@brief theta
     *@see GetCoordinateTheta()
     */
    float photonics_theta_;
    
    /**
     *@brief photonics l coordinate
     *@see GetCoordinateL()
     */
    float photonics_l_;

    /**
     *@brief photonics rho coordinate
     *@see GetCoordinateRho()
     */
    float photonics_rho_;

    /**
     *@brief photoncs phi coordinate
     *@see GetCoordinatePhi()
     */
    float photonics_phi_;

    /**
     *@brief dist 
     *@see GetCoordinateDist()
     */
    float photonics_dist_;

    /**
     *@brief photoncs z src coordinate for track end
     *@see GetCoordinateStopZSrc()
     */
    float photonics_stop_zsrc_;
    
    /**
     *@brief L for track end
     *@see GetCoordinateStopL()
     */
    float photonics_stop_l_;

    /**
     *@brief dist for end of track
     *@see GetCoordinateStopDist()
     */
    float photonics_stop_dist_;

    /**
     *@brief Offset in Z
     */
    float zOrigin_;    

    /**
     *@brief true if muon is finite 
     */
    bool finite_;

    /**
     *@brief true if coordinate is calculated
     */
    bool calculated_;
};

#endif
