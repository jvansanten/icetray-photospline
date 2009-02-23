/**
 *@file
 *@brief Definition of simple photonics PSInterface
 *
 *@author Thomas Burgess
 *
 * (c) the IceCube Collaboration
 *
 * $Revision: 1.13 $
 * $Author$
 * $Date$
 * $Id$
 */
 
#ifndef __PSI_Photonics_h__
#define __PSI_Photonics_h__

//
// Inludes, namespaces, predeclarations
//

//Standard C/C++ 
#include <iostream>
using std::ostream;
#include <vector>
using std::vector; 
#include <string>
using std::string;

#ifndef PSI_DISABLE_ROOT
//Root includes
#include "TObject.h" //Needed for ClassDef
#endif

//Photonics C++ io 
#include "photonicsCPPio.h"

//Local includes
#include "PSInterface/PSI_Logging.h"
#include "PSInterface/PSInterface.h"
#include "PSInterface/PSI_Coordinate.h"
#include "PSInterface/PSI_Coordinate_Photonics.h"

/**
 *@brief Photonics PSInterface implementation
 *
 * Compatible with photonics 1.5 (starlight) or later
 *
 *@author Thomas Burgess
 */
class PSI_Photonics : public PSInterface {

#ifndef PSI_DISABLE_ROOT
    ///Root class definition macro
    ClassDef(PSI_Photonics,2); ///PhotonicsPSInterface implementation    
#endif

public:   

    ///Name log for icetray    
    SET_LOGGER("PSI_Photonics");

    /**
     *@brief Default constructor
     */
    PSI_Photonics();

    /**
     *@brief Destructor
     */
    virtual ~PSI_Photonics();

    /**
     *@brief Load photonics tables
     *
     *@param level      Photonics level (1 or 2)
     *@param fileName   Name of table list file
     *@param verbose    0 for not verbose, 1 for verbose
     * 
     *@return True if successful, false otherwise
     */
    bool LoadTables( 
 	const int& level,
	string fileName = "",
	const int& verbose = 0 );

    /**
     *@brief Clear photonic tables 
     *
     *@param level Photonics level (1 or 2)
     *
     *@return True if tables were cleared, false otherwise
     */
    bool ClearTables( 
	const int& level );
    
    /**
     *@brief Get total table memory usage for (bytes) from photonics
     *
     *@param level Photonics level (1 or 2)
     *
     *@return memory usage
     */
    unsigned long GetMemoryUsage( 
	const int& level);
    
    /**
     *@brief Define the angular region of interest for future requests.
     *
     * If tables are already loaded they are cleared and a new set 
     * is loaded if possible.If no tables are loaded, the region will
     * still be used at the next LoadTable request
     *
     *@param level Photonics level (1 or 2)
     *@param low   Lower angle 
     *@param high  Higher angle
     *
     *@return True on success, false otherwise
     */
    bool SetAngularSelection(
	const int& level,
	const float& low,
	const float& high);    

    /**
     *@brief Define the depth region of interest for future requests.
     *
     * If tables are already loaded they are cleared and a new set 
     * is loaded if possible.If no tables are loaded, the region will
     * still be used at the next LoadTable request
     *
     *@param level Photonics level (1 or 2)
     *@param low   Lower depth
     *@param high  Upper depth
     *
     *@return True on success, false otherwise
     */
    bool SetDepthSelection(
	const int& level,
	const float& low,
	const float& high);    

    /**
     *@brief Get photonics angular interval (may not be the same as requested)
     *
     * If tables are already loaded they are cleared and a new set 
     * is loaded if possible.If no tables are loaded, the region will
     * still be used at the next LoadTable request
     *
     *@param level Photonics level (1 or 2)
     *@param low   Resulting lower angle 
     *@param high  Resulting higher angle
     *
     *@return True on success, false otherwise
     */
    bool GetAngularInterval(
	const int& level,
	float& low,
	float& high);    

    /**
     *@brief Get photonics depth interval (may not be the same as requested)
     *
     * If tables are already loaded they are cleared and a new set 
     * is loaded if possible.If no tables are loaded, the region will
     * still be used at the next LoadTable request
     *
     *@param level Photonics level (1 or 2)
     *@param low   Resulting lower depth
     *@param high  Resulting upper depth
     *
     *@return True on success, false otherwise
     */
    bool GetDepthInterval(
	const int& level,
	float& low,
	float& high);    

    /**
     *@brief Determine memory usage of an angular selection request
     *
     * Before calling a minimal set (for example alow=0,high=0 must be loaded)
     *
     *@param level Photonics level (1 or 2)
     *@param low   Lower angle 
     *@param high  Higher angle
     *
     *@return memory usage on success, 0 otherwise
     */
    unsigned long GetAngularMemoryRequirement(
	const int& level,
	const float& low,
	const float& high) ;
    
    /**
     *@brief Set photonics verbosity
     *
     *@param level Photonics level (1 or 2)
     *@param verbosity verbosity
     *
     *@return True on success, false otherwise
     */
    bool SetPhotonicsVerbosity( 
	const int& level,
	const int &verbosity);

    /**
     *@brief Set group reference refractive indeces.
     *
     * Sets the group refractive index implicitly used for
     * residual time request (both for input and output) 
     *
     *@param level Photonics level (1 or 2)
     *@param ref_ng reference group refractive index
     *
     *@return True on success, false otherwise
     */    
    bool SetRefRefractive(
	const int& level,
	const double& ref_ng);

    /**
     *@brief Get level1 group and phase reference refractive indeces.
     *
     * Get the group and phase refractive indices implicitly used for
     * residual time request (both for input and output) 
     *
     *@param level Photonics level (1 or 2)
     *@param ref_ng reference group refractive index
     *@param ref_np reference phase refractive index
     *
     *@return True on success, false otherwise
     */    
    bool GetRefRefractive(const int& level,double &ref_ng,double &ref_np);
    
    /**
     *@brief Check if photonics tables are loaded
     *
     *@param level Photonics level (1 or 2)
     *
     *@return true if tables are loaded
     */
    bool IsTablesLoaded( 
	const int& level ) const;

    /**
     *@brief Set default driver file name for photonics
     *
     * This name will be used when LoadTables() is called without a file name.
     * After LoadTables() this name will be updated with the name of the file
     * that was loaded.
     *
     *@param level Photonics level (1 or 2)
     *@param fileName photonics driver file name
     *
     *@return True on success, false otherwise
     */
    bool SetDriverFileName( 
	const int& level,
	const string& fileName );

    /**
     *@brief Set default driver file name for photonics
     *
     *@param level Photonics level (1 or 2)
     *
     *@return Driver file name if on success, "" otherwise
     */
    string GetDriverFileName( 
	const int& level ) const;

    /**
     *@brief Get interpolation mode
     *
     *@return interpolation Interpolation mode
     */
    int GetInterpolationMode() const 
    { 
	return interpolationMode_; 
    }

    /**
     *@brief Set interpolation mode
     *
     *@param interpolation Interpolation mode (135 or 255 is recommended)
     */
    void SetInterpolationMode( 
	const int& interpolation );

protected:

    /**
     *@brief Print instance to stream
     *
     *@param o      stream to print to
     */
    virtual void Print(
        ostream &o ) const;

    /**
     *@brief Creates a PSI_Coordinate
     *
     *@see MakeCoordinate
     */
    virtual PSI_Coordinate* Coordinate( 
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
     *@brief Get a mean amplitude
     *     
     *@param amplitude    Amplitude to fill (-1 if failed)
     *@param coordinate   Coordinate for query (create it using MakeCoordinate)
     *
     *@return             true on success false otherwise
     *
     *@see GetMeanAmplitude
     */
    virtual bool MeanAmplitude(
	double &amplitude,
	PSI_Coordinate* coordinate);

    /**
     *@brief Gets a Time Delay
     *
     *@param timeDelay   Time delay to fill (-1 if failed)
     *@param random      Uniform random number [0,1] to sample timedelay 
     *                   distribution with
     *@param coordinate  Coordinate for query (in most cases coodinate 
     *                   must be used in GetMeanAmplitude first)
     *
     *@return             true on success false otherwise
     *
     *@see GetTimeDelay
     */
    virtual bool TimeDelay(
	double &timeDelay,
	const double &random,
	PSI_Coordinate* coordinate);

    /**
     *@brief Get a hit probability
     *
     *@param probability  Probability to fill (-1 if failed)
     *@param timeDelay    Time delay (>0) to sample probability distribution 
     *                    with
     *@param coordinate   Coordinate for query (in most cases coodinate 
     *                    must be used in GetMeanAmplitude first)
     *    
     *@return             true on success false otherwise
     *
     *@see GetProbability
     */
    virtual bool Probability(
	double &probability,
	const double &timeDelay,
	PSI_Coordinate* coordinate);

    /**
     *@brief Change the tables loaded state for photonics
     *
     * Warning invalid states will cause problems!
     *
     *@param level Photonics level (1 or 2)
     *@param state New state of tables loaded state
     *
     *@return true if tables are loaded
     */
    bool SetTablesLoaded( 
	const int& level,
	const bool &state );

    /**
     *@brief Photonics table reader library object.
     *
     * This class is provided by libphotonicsCPPio
     * we make it protected so that other psi instanses can in principle
     * share our table. 
     *
     * photonics_cppio supports level2 only. level2 can be used with any
     * type of light source; showers, muons (finite, infinite) etc...
     *
     */
    photonics_cppio photonics_cppio_obj_;

private:

    /**
     *@brief Copy constructor  
     *
     * Since there is no good way to copy the memory allocated
     * by photonics copying is prohibited and hence this method
     * is private and unimplemented.
     *
     *@param psinterface PSI_Photonics to copy from
     */
    PSI_Photonics( const PSI_Photonics& psinterface );

    /**
     *@brief Assignment operator 
     *
     * Since there is no good way to copy the memory allocated
     * by photonics copying is prohibited and hence this method
     * is private and unimplemented.
     *
     *@param psinterface  PSI_Photonics to assign from
     *
     *@return             assigned om
     */
    const PSI_Photonics& operator=( const PSI_Photonics& psinterface );

    /**
     *@brief Get Photonics level
     *
     *@param coordinate Coordinate used for query
     *@return photonics level or -1 on failure
     */
    int GetLevel( PSI_Coordinate_Photonics* coordinate) const;

    /**
     *@brief True if level2 tables are loaded (in current table set)
     */
    bool level1TablesLoaded_;

    /**
     *@brief True if level2 tables are loaded (in current table set)
     */
    bool level2TablesLoaded_;

    /**
     *@brief photonics level 1 driver file name (for current table set)
     */
    string level1DriverFileName_;

    /**
     *@brief photonics level 2 driver file name (for current table set)
     */
    string level2DriverFileName_;

    /**
     *@brief photonics interpolation mode (often 255 is a good value)
     */
    int interpolationMode_;
};

#endif
