/**
 *@file
 *@brief Definition of PSInterface - Photon Simulation Interface
 *
 *@author Thomas Burgess
 *
 * (c) the IceCube Collaboration
 *
 * $Revision: 1.23 $
 * $Author$
 * $Date$
 * $Id$
 */

#ifndef __PSInterface_h__
#define __PSInterface_h__

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
#include "PSInterface/PSI_Logging.h"
class PSI_Coordinate;

/**
 *@brief PSInterface - Photon Simulation Interface
 *
 * PSInterface is the base class for all Photon Simulation Interface.
 *
 * It has 4 pure abstract methods which has to be implemented in subclasses.
 * Each method has a public visible method to access it:
 * - Coordinate() accessed by MakeCoordinate()
 * - MeanAmplitude() accessed by GetMeanAmplitude()
 * - TimeDelay() accessed by GetTimeDelay()
 * - Probability() accessed by GetProbability()
 *
 *@author Thomas Burgess
 */
class PSInterface {

#ifndef PSI_DISABLE_ROOT
    ///ROOT class definition macro
    ClassDef(PSInterface,2); ///PSInterface - Photon Simulation Interface
#endif

public:
    ///Name log for icetray    
    SET_LOGGER("PSInterface");

    /**
     *@brief Default constructor
     */
    PSInterface();

    /**
     *@brief Destructor
     */
    virtual ~PSInterface();

    /**
     *@brief Copy constructor 
     *
     *@param psinterface PSInterface to copy from
     */
    PSInterface ( const PSInterface& psinterface );

    /**
     *@brief Assignment operator 
     *
     *@param psinterface  PSInterface to assign from
     *
     *@return             assigned om
     */
    const PSInterface& operator=( const PSInterface& psinterface );
    
    /**
     *@brief Printing operator for this and all derived classes
     *
     * Don't overload, to change behaviour in deriving classes
     * override Print() which is called by this method.
     *
     *@param o    Stream to write to
     *@param psi  PSInterface instance to print
     *
     *@result     The written stream
     */
    friend ostream& operator<< (
        ostream &o,
        const PSInterface &psi );
    
    /**
     *@brief Printing operator for this and all derived classes
     *
     * This is the << for pointers to a PSInterface
     *
     * Don't overload, to change behaviour in deriving classes
     * override Print() which is called by this method.
     *
     *@param o    Stream to write to
     *@param psi  PSInterface instance to print
     *
     *@result     The written stream
     */
    friend ostream& operator<< (
        ostream &o,
        const PSInterface *psi );

    /**
     *@brief Create a coordinate
     *
     * Creates a new PSI coordinate (@see PSI_Coordinate)
     *
     * PSI does not own the coordinate and will not delete it after use.
     *
     *@param opticalModuleX            OM location X coordinate (m)
     *@param opticalModuleY            OM location Y coordinate (m)
     *@param opticalModuleZ            OM location Z coordinate (m)
     *@param opticalModuleOrientation  OM orientation (1.0 up, -1.0 down)
     *@param trackX                    Track vertex X coordinate (m)
     *@param trackY                    Track vertex Y coordinate (m)
     *@param trackZ                    Track vertex Z coordinate (m)
     *@param trackTheta                Track zenith angle (degrees)
     *@param trackPhi                  Track azimuth angle (degrees)
     *@param trackLength               Track length (m, negative=infinite)
     *@param trackEnergy               Track energy (GeV)
     *@param trackType                 Track type 
     * (depends on psi implementation, default is 
     * - 0 mu
     * - 1 e+e-
     * - 2 hadronic shower 
     *
     *@return                          Pointer to new PSI_Coordinate 
     *                                 if successful, 0 otherwise
     */
    PSI_Coordinate* MakeCoordinate( 
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
     */
    bool GetMeanAmplitude(
	double &amplitude,
	PSI_Coordinate* coordinate);

    /**
     *@brief Get hit time delay
     *
     *@param timeDelay   Time delay to fill (-1 if failed)
     *@param random      Uniform random number [0,1] to sample timedelay 
     *                   distribution with
     *@param coordinate  Coordinate for query (in most cases coodinate 
     *                   must be used in GetMeanAmplitude first)
     *    
     *@return            true on success false otherwise
     */
    bool GetTimeDelay(
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
     */
    bool GetProbability(
	double &probability,
	const double &timeDelay,
	PSI_Coordinate* coordinate);

    /**
     *@brief Return a Poisson Random Number
     *
     * From a uniform random number 'random'~[0,1] the probability
     * for P(X=k) is calculated for all k>=0 until the cumulant
     * is larger than 'random', the final k will be ~Po(lambda) and
     * is returned.
     *
     * The cumulant is calculated from probability P which in turn is
     * iteratively calculated as P(n+1)=P(n)*lambda/k where P(0)=e^(-lambda)
     *
     * If k exceeds maxPoisson_ the calculation will stop and
     * maxPoisson_ will be returned. (currently 10000)
     *
     * If random is out of the range [0,1] 0 will be returned.
     *
     * Code based on routine from Christian Walck
     *
     *@param lambda  Poisson parameter
     *@param random  Uniform random number
     *
     *@return        Poisson distributed random number, or -1 if failed
     */
    int PoissonRandomNumber( 
	const double &lambda,
	const double &random ) const;

  
    /** 
     *@brief Set angular selection for loading of photonics tables 
     *@param level       Photonics level affected 
     *@param low         Minimum zenith angle (degrees) 
     *@param high        Maximum zenith angle (degrees) 
     */ 
    virtual bool SetAngularSelection(int const& level, float const& low, float const& high);

    /** 
     *@brief Get lower limit of angular selection used for loading of photonics tables 
     *@return Minimum zenith angle (degrees) 
     */ 
    inline double GetAngularSelectionLow() { return angularSelectLow_; } 
    
    /** 
     *@brief Get upper limit of angular selection used for loading of photonics tables 
     *@return Maximum zenith angle (degrees) 
     */ 
    inline double GetAngularSelectionHigh() { return angularSelectHigh_; }

	
protected:
    /**
     *@brief Print instance to stream
     *
     *@param o      stream to print to
     */
    virtual void Print(
        ostream &o ) const;
    
    /**
     *@brief Implementation specifics for GetCoordinate
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
	const int    &trackType) = 0;

    /**
     *@brief Implementation specifics for GetMeanAmplitude()
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
	PSI_Coordinate* coordinate) = 0;

    /**
     *@brief  Implementation specifics for GetTimeDelay()
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
	PSI_Coordinate* coordinate) = 0;

    /**
     *@brief Implementation specifics for GetHitProbability()
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
	PSI_Coordinate* coordinate) = 0;


    
private:
    /**
     *@brief Maximal Poisson random number
     */
    static const unsigned maxPoisson_ = 100000;

protected:
   /**
    * Minimum zenith angle (degrees)
    */
    double angularSelectLow_;

   /**
    * Maximum zenith angle (degrees)
    */
    double angularSelectHigh_;
};

#endif
