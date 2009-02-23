/**
 *@file
 *@brief Definition PTD PSInterface implementation
 *
 * Thomas Burgess
 * (c) the IceCube Collaboration
 *
 * $Revision: 1.17 $
 * $Author$
 * $Date$
 * $Id$
 */
 
#ifndef __PSI_PTD_h__
#define __PSI_PTD_h__

//
// Inludes, namespaces, predeclarations
//

//Standard C/C++ 
#include <vector>
#include <string>
#include <iostream>
using std::ostream;
using std::vector; 
using std::string;

#ifndef PSI_DISABLE_ROOT
//Root includes
#include "TObject.h" //Needed for ClassDef
#endif

//Local includes
#include "PSInterface/PSI_Logging.h"
#include "PSInterface/PSI_Coordinate.h"
class PSInterface;
//PTD includes
extern "C" {
#include "ptdread.h"
}


/**
 *@brief PTD PSInterface implementation
 *
 * PTD supports multiple tables but lacks ability to unload tables.
 * The set or check which table is currently used for amplitude, 
 * timedelay and probability use GetTableId and SetTableId. Each time 
 * a new table is loaded the current table is automatically set to 
 * the new table. The properties of the current photon table can be 
 * accessed using GetSimParams which returns a PTD simParam struct.
 *
 * In order to use PTD like photonics several things has to be taken
 * care of manually (and they are not taken care of by PSI automatically). 
 * - Different particle types use different tables. One needs to manually 
 *   set table id to a table maching the particle.
 * - Muons are differential. Many queries along the track must be added 
 *   together to get a mean amplitude.
 * - Layered tables are not really layered. Layering is simulated by
 *   manually switching PTD table by optical module Z coordinate.
 *
 *@author Thomas Burgess
 */
class PSI_PTD : public PSInterface {

#ifndef PSI_DISABLE_ROOT
    ///Root class definition macro
    ClassDef(PSI_PTD,2); ///PTD PSInterface implementation
#endif

public:   
    ///Name log for icetray    
    SET_LOGGER("PSI_PTD");
    
    /**
     *@brief Default constructor
     */
    PSI_PTD();

    /**
     *@brief Destructor
     */
    virtual ~PSI_PTD();

    /**
     *@brief Load ptd table
     *
     * Will set the current table to the one loaded.
     * 
     *@param fileName   Name of table file to load
     *
     *@return true on success, false otherwise
     */
    bool LoadTable( 
	const string &fileName );

    /**
     *@brief set the ptd mode
     *
     *Valid choises are:
     * - 1    No interpolation
     * - 101  Interpolation
     * - 500  Differential in time, no interpolation
     * - 501  Differential in time, iterpolation
     *
     *@param mode Ptd mode
     *
     *@return true if valid mode selected otherwise false
     */
    bool SetPTDMode( const int &mode );

    /**
     *@brief Get current photon table
     *
     *@return current table id ( -1 if there are no tables)
     */
    int GetTableId() const { return currentTableId_; }

    
    /**
     *@brief Set current photon table
     *
     *@param id table id 
     *
     *@return true if success, false otherwise
     */
    bool SetTableId( const unsigned int& id );
    
    /**
     *@brief Print simparams struct information for current photon table
     */
    void PrintSimParams( ostream &o )const;
    
    /**
     *@brief Get simulation parameters of current photon table
     *
     * fatal if no tables are loaded
     *
     *@return Simulation params
     */
    const simParams& GetSimParams() const;

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
     *@brief Gets a mean amplitude
     *
     *@see GetMeanAmplitude
     */
    virtual bool MeanAmplitude(
	double &amplitude,
	PSI_Coordinate* coordinate);

    /**
     *@brief Gets a Time Delay
     *
     *@see GetTimeDelay
     */
    virtual bool TimeDelay(
	double &timeDelay,
	const double &random,
	PSI_Coordinate* coordinate);

    /**
     *@brief Gets a dummy Probability
     *
     *@see GetProbability
     *
     *@param probability  Probability to fill (-1 if failed)
     *@param timeDelay    Time delay (>0) to sample probability distribution
     *                    with
     *@param coordinate   Coordinate for query (in most cases coodinate
     *                    must be used in GetMeanAmplitude first)
     *
     *@return             true on success false otherwise
     */
    virtual bool Probability(
	double &probability,
	const double &timeDelay,
	PSI_Coordinate* coordinate);

private:
    /**
     *@brief Copy constructor  
     *
     * Since PTD uses globally allocated memory it doens't handle copying 
     * well and and hence this function is unimplemented and declared private
     *
     *@param psinterface PSI_PTD to copy from
     */
    PSI_PTD( const PSI_PTD& psinterface );

    /**
     *@brief Assignment operator 
     *
     * Since PTD uses globally allocated memory it doens't handle copying 
     * well and and hence this function is unimplemented and declared private
     *
     *@param psinterface  PSI_PTD to assign from
     *
     *@return             assigned om
     */
    const PSI_PTD& operator=( const PSI_PTD& psinterface );

    /**
     *@brief PTD simulation parameters 
     */
    vector<simParams> simParams_;

    /**
     *@brief PTD mode
     */
    int mode_;

    /**
     *@brief Current PTD photon table id
     */
    int currentTableId_;

    /**
     *@brief Vector with internal photon table id as returned by PTD
     */
    vector<int> PTDTableIds_;

    /**
     *@brief Vector with photon table filenames
     */
    vector<string> PTDTableFileNames_;
};

#endif
