/**
 *@file
 *@brief Definition of photonics Photorec PSInterface
 *
 *@author Thomas Burgess
 *
 * (c) the IceCube Collaboration
 *
 * $Revision: $
 * $Author:  $
 * $Date: $
 * $Id: $
 */
 
#ifndef __PSI_PhotonicsPhotorec_h__
#define __PSI_PhotonicsPhotorec_h__

//
// Inludes, namespaces, predeclarations
//

//Standard C/C++ 
#include <iostream>
using std::ostream;

#ifndef PSI_DISABLE_ROOT
//Root includes
#include "TObject.h" //Needed for ClassDef
#endif

//Local includes
#include "PSInterface/PSI_Logging.h"
#include "PSInterface/PSI_Coordinate_Photonics.h"
#include "PSInterface/PSI_Photonics.h"

/**
 *@brief Photonics Photorec PSInterface implementation
 *
 * The Photonics Photorec PSInterface can do anything a normal 
 * PSI_Photonics can do. Here is an example of typical photo rec usage.
<pre>
//Initialize
PSI_PhotonicsPhotorec rec;
rec.LoadTables(2,"level2_table.list");

//Create a photonics coordinate for a 
//downlooking optical module at (0,0,0) and a
//100 GeV muon track at (0,0,-10) with angles (45,180) for 10m
double OM_x=0, OM_y=0, OM_z=0, OM_o=-1;
double TR_x=0, TR_y=0, TR_z=-10, TR_th=45, TR_ph=180,
       TR_length=10, TR_energy=100;
int TR_type = 0;
PSI_Coordinate_Photonics* coord = 
    dynamic_cast<PSI_Coordinate_Photonics*>(
        rec.MakeCoordinate(
            OM_x,OM_y,OM_z,OM_o,
            TR_x,TR_y,TR_z,TR_th,TR_ph,
            TR_length,TR_energy,TR_type));
//Define timedelay for query
double timedelay = 25;
//Define amplitude and probability to fetch
double amplitude = -1;
double probability = -1;

if ( coord ! = 0 ) {
  //Make photorec query and print result
  rec.Get(c,timedelay,amplitude,probability);
}
cout << "Amplitude="<<amplitude<<" pdf="<<probability << endl;
</pre>
 *
 *@author Thomas Burgess
 */
class PSI_PhotonicsPhotorec : public PSI_Photonics {

#ifndef PSI_DISABLE_ROOT
    ///Root class definition macro
    ClassDef(PSI_PhotonicsPhotorec,2); ///Photonics Photo Rec PSInterface
#endif
    
public:   
    ///Name log for icetray    
    SET_LOGGER("PSI_PhotonicsPhotorec");

    /**
     *@brief Default constructor
     */
    PSI_PhotonicsPhotorec();

    /**
     *@brief Destructor
     */
    virtual ~PSI_PhotonicsPhotorec();

    /**
     *@brief Get Photorec
     *
     *@param coordinate   Coordinate used for request
     *@param delay        Input time delay for request
     *@param amplitude    Resulting amplitude, for minimum ionizing particle
     *@param probability  Resulting differential probability: 
     *                    dp/dt, integrates to 1
     *
     *@return true on success, false otherwise
     */
    bool Get( 
	PSI_Coordinate_Photonics* coordinate,	      
	const double &delay,
	double &amplitude,	
	double &probability);

protected:
    /**
     *@brief Print instance to stream
     *
     *@param o      stream to print to
     */
    virtual void Print(
        ostream &o ) const;

private:
    /**
     *@brief Copy constructor  
     *
     * Since there is no good way to copy the memory allocated
     * by photonics copying is prohibited and hence this method
     * is private and unimplemented.
     *
     *@param psinterface PSI_PhotonicsPhotorec to copy from
     */
    PSI_PhotonicsPhotorec( const PSI_PhotonicsPhotorec& psinterface );

    /**
     *@brief Assignment operator 
     *
     * Since there is no good way to copy the memory allocated
     * by photonics copying is prohibited and hence this method
     * is private and unimplemented.
     *
     *@param psinterface  PSI_PhotonicsPhotorec to assign from
     *
     *@return             assigned om
     */
    const PSI_PhotonicsPhotorec& operator=( 
	const PSI_PhotonicsPhotorec& psinterface );

    /**
     *@brief Get Mean Amplitude, unusable in photorec
     */
    virtual bool MeanAmplitude(
	double &amplitude,
	PSI_Coordinate* coordinate);

    /**
     *@brief  Get time delay, unusable in photorec
     */
    virtual bool TimeDelay(
	double &timeDelay,
	const double &random,
	PSI_Coordinate* coordinate);

    /**
     *@brief Get Hit Probability, unusable in photorec
     */
    virtual bool Probability(
	double &probability,
	const double &timeDelay,
	PSI_Coordinate* coordinate);

};

#endif
