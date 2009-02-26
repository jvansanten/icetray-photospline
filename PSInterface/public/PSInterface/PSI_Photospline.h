/**
 *@file
 *@brief Definition of simple photonics PSInterface
 *
 *@author Thomas Burgess
 *
 * (c) the IceCube Collaboration
 *
 * $Revision: 1.13 $
 * $Author: nwhitehorn $
 * $Date: 2009-02-23 10:54:11 -0600 (Mon, 23 Feb 2009) $
 * $Id: PSI_Photospline.h 52769 2009-02-23 16:54:11Z nwhitehorn $
 */
 
#ifndef __PSI_Photospline_h__
#define __PSI_Photospline_h__

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

//Local includes
#include "PSInterface/PSI_Logging.h"
#include "PSInterface/PSInterface.h"
#include "PSInterface/PSI_Coordinate.h"
#include "PSInterface/PSI_Coordinate_Photonics.h"

struct splinetable;

/**
 *@brief Photonics PSInterface implementation
 *
 * Compatible with photonics 1.5 (starlight) or later
 *
 *@author Thomas Burgess
 */
class PSI_Photospline : public PSInterface {

#ifndef PSI_DISABLE_ROOT
    ///Root class definition macro
    ClassDef(PSI_Photospline,2); ///PhotonicsPSInterface implementation    
#endif

public:   

    ///Name log for icetray    
    SET_LOGGER("PSI_Photospline");

    /**
     *@brief Default constructor
     */
    PSI_Photospline();

    /**
     *@brief Destructor
     */
    virtual ~PSI_Photospline();

    /**
     *@brief Load photonics tables
     *
     *@param fileName   Name of table FITS file
     * 
     *@return True if successful, false otherwise
     */
    bool LoadTables(std::string fileName = "");

    /**
     *@brief Clear photonic tables 
     *
     *@return True if tables were cleared, false otherwise
     */
    bool ClearTables();
    
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

private:

    /**
     *@brief Copy constructor  
     *
     * Since there is no good way to copy the memory allocated
     * by photonics copying is prohibited and hence this method
     * is private and unimplemented.
     *
     *@param psinterface PSI_Photospline to copy from
     */
    PSI_Photospline( const PSI_Photospline& psinterface );

    /**
     *@brief Assignment operator 
     *
     * Since there is no good way to copy the memory allocated
     * by photonics copying is prohibited and hence this method
     * is private and unimplemented.
     *
     *@param psinterface  PSI_Photospline to assign from
     *
     *@return             assigned om
     */
    const PSI_Photospline& operator=( const PSI_Photospline& psinterface );

    struct splinetable splinetable;
};

#endif
