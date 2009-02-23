/**
 *@file
 *@brief Implementation of PTD PSInterface 
 *
 * Thomas Burgess
 * (c) the IceCube Collaboration
 *
 * $Revision: 1.21 $
 * $Author$
 * $Date$
 * $Id$
 */

//Standard C/C++ includes
#include <iostream>
#include <vector>
using std::vector;
using std::ostream;

//Local includes
#include "PSInterface/PSInterface.h"
#include "PSInterface/PSI_Logging.h"
#include "PSInterface/PSI_PTD.h"
#include "PSInterface/PSI_Coordinate_PTD.h" 
#include "PSInterface/PSI_Coordinate.h" 
//PTD includes
extern "C" {
#include "ptdread.h"
}

#ifndef PSI_DISABLE_ROOT
//Root includes
#include "TObject.h" //Needed for ClassImp
///Root class implementation macro
ClassImp(PSI_PTD);
#endif

//Default constructor
PSI_PTD::PSI_PTD() :
    mode_(1),
    currentTableId_(-1)
{

}

//Destructor
PSI_PTD::~PSI_PTD()
{
}

//Load ptd table
bool PSI_PTD::LoadTable( 
    const string &fileName )
{
    photonVals *coord=0;
    const int mode = 0; //0 for load_tables in getphotons
    int PTDTableId = 0; //PTD internal table id
    simParams sim;

    sim.nFileId=0;
    sim.fTrackStepSize=0;
    sim.maxPhoton=0;
    sim.sourcePosition[0]=0;
    sim.sourcePosition[1]=0;
    sim.sourcePosition[2]=0;
    sim.maxAlfa=0;
    sim.maxScatter=0;
    sim.nCoordsys=0;
    sim.OMType=0;
    sim.scaModels=0;
    sim.BeamOn=0;
    sim.lambdaAbs=0; 
    sim.maxSimDist=0;
    sim.injectedPhotons=0;
    sim.dummy[0]=0;   
    sim.dummy[1]=0;   
    sim.dummy[2]=0;   
    sim.dummy[3]=0;   
    sim.ra_par1=0;  
    sim.th_par1=0;  
    sim.ti_par1=0;  
    sim.ph_par1=0;  
    sim.rh_par1=0;  
    sim.ra_bins=0;  
    sim.th_bins=0;  
    sim.ti_bins=0;  
    sim.ph_bins=0;  
    sim.rh_bins=0;  
    sim.ra_max=0;   
    sim.th_max=0;   
    sim.ti_max=0;   
    sim.ph_max=0;   
    sim.rh_max=0;   
    sim.ra_min=0;   
    sim.th_min=0; 
    sim.ti_min=0; 
    sim.ph_min=0; 
    sim.rh_min=0; 
    sim.ra_step=0;
    sim.th_step=0;
    sim.ti_step=0;
    sim.ph_step=0;
    sim.rh_step=0;
    sim.ra_over=0;
    sim.th_over=0;
    sim.ti_over=0;
    sim.ph_over=0;
    sim.rh_over=0;
    sim.po_bins=0;

    //Call getphotons
    int errorCode = 
	getphotons(fileName.c_str(), &PTDTableId, &mode, coord, &sim);
    if ( errorCode < 0 ) {
	log_error("Failed to load %s, ptd error %d", 
		  fileName.c_str(), errorCode);
	return false;
    }
    PTDTableFileNames_.push_back(fileName);
    simParams_.push_back(sim);
    PTDTableIds_.push_back(PTDTableId);

    log_debug("Loaded table %s", fileName.c_str());

    currentTableId_++;

    return true;
}

//set the ptd mode
bool PSI_PTD::SetPTDMode( 
    const int &mode )
{
    //Check mode
    if ( (mode==1)||(mode==101)||(mode==500)||(mode==501) ) {
	mode_=mode;
	log_debug("Set PTD mode %d", mode_);
	return true; 		
    } 

    log_warn("Invalid mode %d. Not doing anything. "
	     "(try 1, 101, 500, 501 or reading the docs...)",
	     mode);

    return false;
}
    
//Set current photon table
bool PSI_PTD::SetTableId( 
    const unsigned int& id )
{
    //Check id
    if ( id>=simParams_.size() ) {
	log_error("id %d to large when %zu tables",id,simParams_.size());
	return false;
    }
    
    currentTableId_=id;
    log_debug("Set current table id to %d", currentTableId_);

    return true;
}
    
//Get simulation parameters of current photon table
const simParams& PSI_PTD::GetSimParams() const
{
    //Check if there are tables
    if ( simParams_.size() == 0 ) {
	log_fatal("No PTD tables loaded!");
    }
    return simParams_[currentTableId_];
}


void PSI_PTD::PrintSimParams( ostream &o ) const
{
    o << "SimParams for table \"" << PTDTableFileNames_[currentTableId_] 
      << "\"\n"
	"\tnFileId "<< simParams_[currentTableId_].nFileId << "\n"
	"\tfTrackStepSize " << 
	simParams_[currentTableId_].fTrackStepSize << "\n"
	"\tmaxPhoton " << simParams_[currentTableId_].maxPhoton << "\n"
	"\tsourcePosition[0] " << 
	simParams_[currentTableId_].sourcePosition[0] << "\n"
	"\tsourcePosition[1] " << 
	simParams_[currentTableId_].sourcePosition[1] << "\n"
	"\tsourcePosition[2] " << 
	simParams_[currentTableId_].sourcePosition[2] << "\n"
	"\tmaxAlfa " << simParams_[currentTableId_].maxAlfa << "\n"
	"\tmaxScatter " << simParams_[currentTableId_].maxScatter << "\n"
	"\tnCoordsys " << simParams_[currentTableId_].nCoordsys << "\n"
	"\tOMType " << simParams_[currentTableId_].OMType << "\n"
	"\tscaModels " << simParams_[currentTableId_].scaModels << "\n"
	"\tBeamOn " << simParams_[currentTableId_].BeamOn << "\n"
	"\tlambdaAbs " << simParams_[currentTableId_].lambdaAbs << "\n"
	"\tmaxSimDist " << simParams_[currentTableId_].maxSimDist << "\n"
	"\tinjectedPhotons " << 
	simParams_[currentTableId_].injectedPhotons << "\n"
	"\tdummy[0] " << simParams_[currentTableId_].dummy[0] << "\n"
	"\tdummy[1] " << simParams_[currentTableId_].dummy[1] << "\n"
	"\tdummy[2] " << simParams_[currentTableId_].dummy[2] << "\n"
	"\tdummy[3] " << simParams_[currentTableId_].dummy[3] << "\n"
	"\tra_par1 " << simParams_[currentTableId_].ra_par1 << "\n"
	"\tth_par1 " << simParams_[currentTableId_].th_par1 << "\n"
	"\tti_par1 " << simParams_[currentTableId_].ti_par1 << "\n"
	"\tph_par1 " << simParams_[currentTableId_].ph_par1 << "\n"
	"\trh_par1 " << simParams_[currentTableId_].rh_par1 << "\n"
	"\tra_bins " << simParams_[currentTableId_].ra_bins << "\n"
	"\trh_bins " << simParams_[currentTableId_].th_bins << "\n"
	"\tti_bins " << simParams_[currentTableId_].ti_bins << "\n"
	"\tph_bins " << simParams_[currentTableId_].ph_bins << "\n"
	"\tth_bins " << simParams_[currentTableId_].rh_bins << "\n"
	"\tra_max " << simParams_[currentTableId_].ra_max << "\n"
	"\tth_max " << simParams_[currentTableId_].th_max << "\n"
	"\tti_max " << simParams_[currentTableId_].ti_max << "\n"
	"\tph_max " << simParams_[currentTableId_].ph_max << "\n"
	"\tth_max " << simParams_[currentTableId_].rh_max << "\n"
	"\tra_min " << simParams_[currentTableId_].ra_min << "\n"
	"\tth_min " << simParams_[currentTableId_].th_min << "\n"
	"\tti_min " << simParams_[currentTableId_].ti_min << "\n"
	"\tph_min " << simParams_[currentTableId_].ph_min << "\n"
	"\trh_min " << simParams_[currentTableId_].rh_min << "\n"
	"\tra_step " << simParams_[currentTableId_].ra_step << "\n"
	"\tth_step " << simParams_[currentTableId_].th_step << "\n"
	"\tti_step " << simParams_[currentTableId_].ti_step << "\n"
	"\tph_step " << simParams_[currentTableId_].ph_step << "\n"
	"\trh_step " << simParams_[currentTableId_].rh_step << "\n"
	"\tra_over " << simParams_[currentTableId_].ra_over << "\n"
	"\tth_over " << simParams_[currentTableId_].th_over << "\n"
	"\tti_over " << simParams_[currentTableId_].ti_over << "\n"
	"\tph_over " << simParams_[currentTableId_].ph_over << "\n"
	"\trh_over " << simParams_[currentTableId_].rh_over << "\n"
	"\tpo_bins " << simParams_[currentTableId_].po_bins << "\n"
	"\tpo_max " << simParams_[currentTableId_].po_max << "\n"
	"\tpo_min " << simParams_[currentTableId_].po_min << "\n"
	"\tpo_step " << simParams_[currentTableId_].po_step << "\n"
	"\tpo_over " << simParams_[currentTableId_].po_over << "\n";
}

//Print instance to stream
void PSI_PTD::Print(
    ostream &o ) const
{
    PSInterface::Print(o);
    o << "\nPSI_PTD\n\t"
	"NTable " << simParams_.size() << "\n\t"
	"current table id " << currentTableId_ << "\n\t"
	"mode " << mode_ <<
	"\nPhoton table files loaded:\n";
    for (unsigned int i=0;i<PTDTableFileNames_.size();i++)
	o << "\t" << PTDTableFileNames_[i] << "\n";
}

//Creates a PSI_Coordinate
PSI_Coordinate* PSI_PTD::Coordinate( 
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
    const int    &trackType)
{
    PSI_Coordinate_PTD* coord = new PSI_Coordinate_PTD(
	opticalModuleX,
	opticalModuleY,
	opticalModuleZ,
	opticalModuleOrientation,
	trackX,
	trackY,
	trackZ,
	trackTheta,
	trackPhi,
	trackLength,
	trackEnergy,
	trackType);
    
    return coord;
}
    
//Gets a  mean amplitude
bool PSI_PTD::MeanAmplitude(
    double &amplitude,
    PSI_Coordinate* coordinate)
{
    //Initialize amplitude
    amplitude = -1;
        
    //Check that at least one table is loaded
    if (simParams_.size()==0) {
	log_error("No tables loaded cannot get a mean amplitude");
	return false;
    }

    //Get and check coordinate
    PSI_Coordinate_PTD* coord = 
	dynamic_cast<PSI_Coordinate_PTD*>(coordinate);
    if ( coord == 0 ) {
	log_error("MeanAmplitude failed, "
		  "did not recieve a PTD coordinate");
	return false;
    }
    if ( ! coord->CalculateCoordinate( simParams_[currentTableId_] ) ) {
	log_error("Invalid coordinate - "
		  "probably outside PTD table boundaries");
	return false;
    }

    //Call PTD and get the mean amplitude
    const char* str="";
    int errorCode  = -1;
    if ( ( errorCode = getphotons(
	       str,&PTDTableIds_[currentTableId_],&mode_, 
	       &coord->photonVals_, &simParams_[currentTableId_] ) ) != 0 )  
    {
	log_error("MeanAmplitude failed, PTD error %d", errorCode);
	return false;
    }

    amplitude  = coord->photonVals_.amp;  
    if ( coord->GetTrackType()== 2  ) {
	amplitude*=0.802; //Ems/Hadron effective length ratio
    }	
    
    log_debug("Got MeanAmplitude=%f from PTD",amplitude);

    return true;
}

//Gets a  Time Delay
bool PSI_PTD::TimeDelay(
    double &timeDelay,
    const double &random,
    PSI_Coordinate* coordinate)
{
    //Initialize time delay
    timeDelay=-1;

    //Check that at least one table is loaded
    if (simParams_.size()==0) {
	log_error("No tables loaded cannot get a time delay");
	return false;
    }
    
    //Get and check coordinate
    PSI_Coordinate_PTD* coord = 
	dynamic_cast<PSI_Coordinate_PTD*>(coordinate);
    if ( coord == 0 ) {
	log_error("TimeDelay failed, "
		  "did not recieve a PTD coordinate");
	return false;
    }

    //Put probability in coordinate
    coord->photonVals_.random = random;
		
    //Read amplitude from PTD
    const char* str="";
    int errorCode = -1;
    if ( ( errorCode = getphotons(
	       str, &PTDTableIds_[currentTableId_],&mode_,
	       &coord->photonVals_, &simParams_[currentTableId_] ) ) != 0 )
    {
	log_error("TimeDelay failed, PTD error %d", errorCode);
	return false;
    } 

    timeDelay = coord->photonVals_.time;

    log_debug("Got time delay %f from PTD", timeDelay);

    return true;
}

//Gets a  Probability
bool PSI_PTD::Probability(
    double &probability,
    const double &timeDelay,
    PSI_Coordinate* coordinate)
{
    //Initialize probability
    probability = -1.0;

    //Check that at least one table is loaded
    if (simParams_.size()==0) {
	log_error("No tables loaded cannot get a probability");
	return false;
    }
    
    //Get and check coordinate
    PSI_Coordinate_PTD* coord = 
	dynamic_cast<PSI_Coordinate_PTD*>(coordinate);
    if ( coord == 0 ) {
	log_error("Probability failed, "
		  "did not recieve a PTD coordinate");
	return false;
    }

    //Put time delay in coordinate
    coord->photonVals_.time = timeDelay;
	    
    //Read amplitude from PTD
    const char* str="";
    int errorCode = -1;
    if ( ( errorCode =  getphotons(
	       str, &PTDTableIds_[currentTableId_],&mode_,
	       &coord->photonVals_, &simParams_[currentTableId_] ) ) != 0 )
    {
	log_error("Probability failed, PTD error %d", errorCode );
	return false;
    }
    
    probability = coord->photonVals_.random;	

    log_debug("Got probability %f from PTD", probability);

    return true;
}
