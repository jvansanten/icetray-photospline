/**
 *@file
 *@brief Definition of PTD Track-OM coordinate class
 *
 *@author Thomas Burgess
 *
 * (c) the IceCube Collaboration
 *
 * $Revision: 1.2 $
 * $Author$
 * $Date$
 * $Id$
 */
 
#ifndef __PSI_Coordinate_PTD_h__
#define __PSI_Coordinate_PTD_h__

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
class PSI_Coordinate;
class PSI_PTD;
//PTD includes
extern "C" {
#include "ptdread.h"
}

/**
 *@brief A track-OM coordinate class for PTD
 *
 * This is a coordinate that can calculate and hold PTD specific information.
 * It is used in PSI_PTD. Before using the coordinate CaclulateCoordinate() 
 * should be called with the simParam that was initialized when loading tables.
 *                                                         
 *@author Thomas Burgess                                   
 */
class PSI_Coordinate_PTD : public PSI_Coordinate {

#ifndef PSI_DISABLE_ROOT
    ///ROOT class definition macro
    ClassDef(PSI_Coordinate_PTD, 2); ///A track-OM coordinate class
#endif

public:
    ///Name log for icetray    
    SET_LOGGER("PSI_Coordinate_PTD");

    ///Allow PSI_PTD to access coordinate private members
    friend class PSI_PTD;

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
     *@param trackLength               Track length (m, negative=infinite
     *@param trackEnergy               Track energy (GeV)
     *@param trackType                 Track type (mu:0,e+e-:1,hadr:2)
     *@param zOrigin                   Shift in z (m)
     */
    PSI_Coordinate_PTD(	
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
	const int    &trackType = 0,
	const double &zOrigin = 0);

    /**
     *@brief Copy constructor
     *
     *@param coordinate coordinate to copy from
     */
    PSI_Coordinate_PTD( 
	const PSI_Coordinate_PTD& coordinate );

    /**
     *@brief Virtual Destructor
     */
    virtual ~PSI_Coordinate_PTD();    

    /**
     *@brief Assignment operator
     *
     *@param coordinate coordinate to assign from
     *@return Assigned coordinate
     */
    const PSI_Coordinate_PTD& operator=( 
	const PSI_Coordinate_PTD& coordinate );

    /**
     *@brief Calculate PTD coordinate
     */
    bool CalculateCoordinate(
 	const simParams &simparams );

    /**
     *@brief Get PTD photon values (PTD internal coordinate)
     *@return handle to photonVals struct
     */
    const photonVals& GetPhotonVals() const { return photonVals_; } 
    
protected:
    /**
     *@brief Print coordinate to stream
     *
     *@param o      stream to print to
     */
    virtual void Print( 
	ostream &o ) const;    

    /**
     *@brief PTD coordinate data struct
     */
    photonVals photonVals_;

    /**
     *@brief Shift all z values by this amount
     *
     * default is 0.0
     */
     double zOrigin_;
};

#endif
