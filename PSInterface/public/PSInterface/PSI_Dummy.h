/**
 *@file
 *@brief Definition of Minimal class implementing PSInterface
 *
 * Thomas Burgess
 * (c) the IceCube Collaboration
 *
 * $Revision: 1.5 $
 * $Author$
 * $Date$
 * $Id$
 */
 
#ifndef __PSI_Dummy_h__
#define __PSI_Dummy_h__

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
class PSI_Dummy : public PSInterface {
#ifndef PSI_DISABLE_ROOT
     ///Root class definition macro
     ClassDef(PSI_Dummy,2); ///Minimal class implementing PSInterface
#endif
public:   
    ///Name log for icetray
    SET_LOGGER("PSI_Dummy");
    
    /**
     *@brief Default constructor
     */
    PSI_Dummy();

    /**
     *@brief Destructor
     */
    virtual ~PSI_Dummy();

    /**
     *@brief Copy constructor 
     *
     *@param psinterface PSI_Dummy to copy from
     */
    PSI_Dummy ( const PSI_Dummy& psinterface );

    /**
     *@brief Assignment operator 
     *
     *@param psinterface  PSI_Dummy to assign from
     *
     *@return             assigned om
     */
    const PSI_Dummy& operator=( const PSI_Dummy& psinterface );
    
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
