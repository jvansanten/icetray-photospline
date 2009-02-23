/**
 *@file
 *@brief Definition of Minimal class implementing PSInterface
 *
 * Thomas Burgess
 * (c) the IceCube Collaboration
 *
 * $Revision: 1.5 $
 * $Author: burgess $
 * $Date: 2005-10-23 11:48:01 -0400 (Sun, 23 Oct 2005) $
 * $Id: PSI_FATInjector.h 11767 2005-10-23 15:48:01Z burgess $
 */
 
#ifndef __PSI_FATInjector_h__
#define __PSI_FATInjector_h__

//
// Inludes, namespaces, predeclarations
//

//C/C++ 
#include <ostream>
using std::ostream;

#ifndef PSI_DISABLE_ROOT
//Root 
#include "TObject.h" //Needed for ClassDef
#endif

//PSI
#include "PSInterface/PSInterface.h"
#include "PSInterface/PSI_Logging.h"
class PSI_Coordinate;


/**
 *@brief Minimal class implementing PSInterface
 *
 * This class is for use in testing and will not provide any useful physical 
 * results.
 *
 *@author Thomas Burgess
 */
class PSI_FATInjector : public PSInterface {
#ifndef PSI_DISABLE_ROOT
     ///Root class definition macro
     ClassDef(PSI_FATInjector,2); ///Minimal class implementing PSInterface
#endif
public:   
    ///Name log for icetray
    SET_LOGGER("PSI_FATInjector");
    
    /**
     *@brief Default constructor
     */
    PSI_FATInjector();

    /**
     *@brief Destructor
     */
    virtual ~PSI_FATInjector();

    /**
     *@brief Copy constructor 
     *
     *@param psinterface PSI_FATInjector to copy from
     */
    PSI_FATInjector ( const PSI_FATInjector& psinterface );

    /**
     *@brief Assignment operator 
     *
     *@param psinterface  PSI_FATInjector to assign from
     *
     *@return             assigned om
     */
    const PSI_FATInjector& operator=( const PSI_FATInjector& psinterface );
    
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
     *@brief Gets a DUMMY mean amplitude
     *
     * Will set amplitude to trackEnergy
     *
     *@see GetMeanAmplitude
     */
    virtual bool MeanAmplitude(
	double &amplitude,
	PSI_Coordinate* coordinate);

    /**
     *@brief Gets a DUMMY Time Delay
     *
     * Will set time delay to random
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
     * Will set probability to timeDelay
     *
     *@see GetProbability
     */
    virtual bool Probability(
	double &probability,
	const double &timeDelay,
	PSI_Coordinate* coordinate);
private:
};

#endif
