/**
 *@file
 *@brief Implementation of simple photonics PSInterface
 *
 *@author Thomas Burgess
 *
 * (c) the IceCube Collaboration
 *
 * $Revision: 1.36 $
 * $Author$
 * $Date$
 * $Id$
 */


//Standard C/C++ includes
#include <iostream>
using std::ostream;

//Photonics includes
extern "C" {
#include "photoamasim.h"
#include "level2.h"    
}
#include "photonicsCPPio.h"
#include "photoamasim.h"
#include "level2.h"    

//Local includes
#include "PSInterface/PSInterface.h"
#include "PSInterface/PSI_Logging.h"
#include "PSInterface/PSI_Photonics.h"
#include "PSInterface/PSI_Coordinate_Photonics.h" 
#include "PSInterface/PSI_Coordinate.h" 

//Check for photonics C++ method.
#ifndef PHOTONICS_CPPMODE_v001
#error "You seem not to have an up to date development photonics version: http://photonics.tsl.uu.se/icecube.php"
#endif

#ifndef PSI_DISABLE_ROOT
//Root includes
#include "TObject.h" //Needed for ClassImp
///Root class implementation macro
ClassImp(PSI_Photonics);
#endif


//Default constructor
PSI_Photonics::PSI_Photonics() :
    PSInterface(),
    photonics_cppio_obj_(),
    level1TablesLoaded_(false),
    level2TablesLoaded_(false),    
    level1DriverFileName_(""),
    level2DriverFileName_(""),
    interpolationMode_(135)
{
}

//Destructor
PSI_Photonics::~PSI_Photonics()
{
    if ( level1TablesLoaded_ ) {
	if (!ClearTables(1)) {
	    log_warn("Failed to clear level 1 tables");
	}
    }
    if ( level2TablesLoaded_ ) {
	if (!ClearTables(2)) {
	    log_warn("Failed to clear level 2 tables");
	}
    }
}

string PSI_Photonics::GetDriverFileName( 
    const int& level ) const
{
    if ( level==1 ) {
	return level1DriverFileName_;
    } else if ( level==2 ) {
	return level1DriverFileName_;
    }
    log_error("GetDriverFileName failed! "
	      "Invalid photonics level %d, select 1 or 2, returning \"\"",
	      level);
    return string("");
}


//Load photonics tables
bool PSI_Photonics::LoadTables( 
    const int& level ,
    string fileName,
    const int& verbose)
{

    //Check level
    if ( ( level != 1 ) && ( level != 2 ) ) {
	log_error("LoadTables failed! "
		  "Invalid photonics level %d, select 1 or 2",level);
	return 0;
    } 

    //Check if tables are allready loaded
    if ( IsTablesLoaded(level) ) {
	log_error("LoadTables failed! "
		  "Level %d tables allready loaded, clear tables first",
		  level);
	return false;
    }

    //Use default file name if neccessary
    if ( fileName.length()==0 ) {
	if ( GetDriverFileName(level).length()>0 ) {
	    fileName=level1DriverFileName_;
	    log_debug("Using %s as level %d driver file name", 
		      fileName.c_str(),level);
	} else {
	    log_error("LoadTables failed! "
		      "PSI_Photonics doesn't know the name "
		      "of the file you want to load");
	    return false;
	}
    }

    //Load tables
    if ( level==1 ) {
	if ( load_tables( fileName.c_str(), verbose ) == 0 ) {
	    //Loading failed
	    log_error("LoadTables failed! "
		      "Loading of level 1 driver file \"%s\" failed!\n", 
		      fileName.c_str());
	    return false;
	}
	//Set tables loaded flag and remember old file name
	level1TablesLoaded_ = true;
	level1DriverFileName_ = fileName;
    } else if ( level == 2 ) {
	if ( photonics_cppio_obj_.load_level2_tables( 
		 fileName.c_str(), verbose ) == 0 ) {
	    //Loading failed
	    log_error("LoadTables failed! "
		      "Loading of level 2 driver file \"%s\" failed!\n", 
		      fileName.c_str());
	    return false;
	}
	//Set tables loaded flag and remember old file name
	level2TablesLoaded_ = true;
	level2DriverFileName_ = fileName;
    }

    log_debug("Loaded photonics level %d driver file %s",
	      level, fileName.c_str());

    return true;
}

//Clear photonics tables 
bool PSI_Photonics::ClearTables(
    const int& level )
{
    //Check level
    if ( ( level != 1 ) && ( level != 2 ) ) {
	log_error("ClearTables failed! "
		  "Invalid photonics level %d, select 1 or 2",
		  level);
	return false;
    } 

    //Check if tables are loaded
    if (!IsTablesLoaded(level)) {
	log_error("ClearTables failed! "
		  "Trying to clear level %d tables, "
		  "but no tables loaded",level);
	return false;	
    }

    //Free tables
    if ( level == 1 ) {
	free_tablesets_level1(); 
	level1TablesLoaded_ = false;    
    } else if ( level == 2 ) {
	photonics_cppio_obj_.free_tablesets_level2(); 
	level2TablesLoaded_ = false;
    } 

    log_debug("Clearing level %d tables",level);
    return true;

}

//Get total table memory usage (bytes) from photonics
unsigned long PSI_Photonics::GetMemoryUsage(
    const int& level )
{
  unsigned long memuse;
  
    if ( level == 1 ) {
	return get_level1_memory_usage();
    } else if ( level == 2 ) {
	photonics_cppio_obj_.get_level2_memory_usage(&memuse);
    } 

    log_error("GetMemoryUsage failed! "
	      "Invalid photonics level %d, select 1 or 2",level);
    return 0;
}

bool PSI_Photonics::SetAngularSelection(
    const int& level,
    const float& low,
    const float& high) 
{
    if ( level == 1 ) {
	set_level1_angular_interval(low,high);
    } else if ( level == 2 ) {
	photonics_cppio_obj_.set_level2_angular_interval(low,high);
    } else {
	log_error("SetAngularSelection failed! "
		  "Invalid photonics level %d, select 1 or 2",
		  level);
	return false;
    }
    
    //Reload tables if neccesary
    if ( IsTablesLoaded(level) ) {
	ClearTables(level);
	if (!LoadTables(level)) {
	    log_error("SetAngularSelection failed! "
		      "Failed to reload tables");
	    return false;
	}
    }

    log_debug("Set level %d angular selection to %f<angle<%f",
	      level,low,high);
    // store angular selection so it can be queried by modules using service
    angularSelectLow_ = low;
    angularSelectHigh_ = high;

    return true;    
}


bool PSI_Photonics::SetDepthSelection(
    const int& level,
    const float& low,
    const float& high) 
{
    if ( level == 1 ) {
	set_level1_depth_interval(low,high);
    } else if ( level == 2 ) {
	photonics_cppio_obj_.set_level2_depth_interval(low,high);
    } else {
	log_error("SetDepthSelection failed! "
		  "Invalid photonics level %d, select 1 or 2",
		  level);
	return false;
    }
    
    //Reload tables if neccesary
    if ( IsTablesLoaded(level) ) {
	ClearTables(level);
	if (!LoadTables(level)) {
	    log_error("SetDepthSelection failed! "
		      "Failed to reload tables");
	    return false;
	}
    }

    log_debug("Set level %d depth selection to %f<depth<%f",
	      level,low,high);

    return true;    
}

bool PSI_Photonics::GetAngularInterval(
    const int& level,
    float& low,
    float& high) 
{
    if ( level == 1 ) {
	get_level1_angular_interval(&low,&high);
    } else if ( level == 2 ) {
	photonics_cppio_obj_.get_level2_angular_interval(&low,&high);
    } else {
	log_error("GetAngularInterval failed! "
		  "Invalid photonics level %d, select 1 or 2",
		  level);
	return false;
    }
    return true;    
}

bool PSI_Photonics::GetDepthInterval(
    const int& level,
    float& low,
    float& high) 
{
    if ( level == 1 ) {
	get_level1_depth_interval(&low,&high);
    } else if ( level == 2 ) {
	photonics_cppio_obj_.get_level2_depth_interval(&low,&high);
    } else {
	log_error("GetDepthInterval failed! "
		  "Invalid photonics level %d, select 1 or 2",
		  level);
	return false;
    }
    return true;    
}

//Determine memory usage of an angular selection request
unsigned long PSI_Photonics::GetAngularMemoryRequirement(
    const int &level,
    const float &low,
    const float &high) 
{
    //Check level
    if ( ( level != 1 ) && ( level != 2 ) ) {
	log_error("GetAngularMemoryRequirement failed! "
		  "Invalid photonics level %d, select 1 or 2",
		  level);
	return 0;
    } 

    //Check if tabels are loaded
    if ( !IsTablesLoaded(level) ) {
	log_error("GetAngularMemoryRequirement failed! "
		  "You must have at least a minimal set of loaded "
		  "in order to estimate memory usage");
	return 0;
    }

    float a=low;
    float b=high;

    //Get the memory usage
    unsigned long result = 0;
    if ( level == 1 ) { 
	result = get_level1_angular_memory_requirement(low,high);
    } else if ( level == 2 ) {
	photonics_cppio_obj_.get_level2_angular_memory_requirement(
	    a,b,&result);
    }

    return result;
}

//Set photonics verbosity level
bool PSI_Photonics::SetPhotonicsVerbosity( 
    const int &level,
    const int &verbosity)
{
    if ( level == 1 ) {
	set_level1_verbosity(verbosity);
    } else if ( level == 2 ) {
	photonics_cppio_obj_.set_level2_verbosity(verbosity);
    } else {
	log_error("SetPhotonicsVerbosity failed! "
		  "Invalid photonics level %d, select 1 or 2",level);
	return false;
    }

    log_debug("Set photonics verbosity for level %d to %d",
	      level,verbosity);
    
    return true;
}

//Get the photonics residual time convention for loadade tables
bool PSI_Photonics::GetRefRefractive(
    const int& level,
    double &ref_ng,
    double &ref_np) 
{
    if ( level == 1 ) {
	return get_level1_reference_refractive(&ref_ng,&ref_np);
    } else if ( level == 2 ){
	return photonics_cppio_obj_.get_level2_reference_refractive(
	    &ref_ng,&ref_np);
    } 
    log_error("GetRefRefractive failed! "
	      "Invalid photonics level %d, select 1 or 2",level);
    return false; 
}

//Set the photonics residual time convention (ngroup) for loadade tables
bool PSI_Photonics::SetRefRefractive(
    const int& level, 
    const double& ref_ng){
    if ( level == 1 ) {
	return set_level1_reference_refractive(ref_ng);
    } else if ( level == 2 ) {
	return photonics_cppio_obj_.set_level2_reference_refractive(ref_ng);
    }
    log_error("SetRefRefractive failed! "
	      "Invalid photonics level %d, select 1 or 2",level);
    return false; 
}


//Check if photonics tables are loaded
bool PSI_Photonics::IsTablesLoaded( 
    const int &level ) const
{
    if ( level == 1 ) {
	return level1TablesLoaded_;
    } else if ( level == 2 ) {
	return level2TablesLoaded_;
    } 
    log_error("IsTablesLoaded failed! "
	      "Invalid photonics level %d, select 1 or 2",level);
    return false;
}

//Set default driver file name for photonics 
bool PSI_Photonics::SetDriverFileName( 
    const int &level,
    const string& fileName )
{
    if ( level == 1 ) {
	level1DriverFileName_=fileName;
    } else if ( level == 2 ) {
	level2DriverFileName_=fileName;
    } else {
	log_error("SetDriverFileName failed! "
		  "Invalid photonics level %d, select 1 or 2",level);
	return false;
    }
    log_debug("Set level %d driver file name to %s",
	      level, fileName.c_str());
    return true;
}

//Set interpolation mode
void PSI_Photonics::SetInterpolationMode( 
    const int& interpolation )
{
    log_debug(
	"Photonics interpolation mode set to %d",
	interpolation);
    interpolationMode_=interpolation;
}

//Print instance to stream
void PSI_Photonics::Print(
    ostream &o ) const
{
    PSInterface::Print(o);
    o << "PSI_Photonics";
}
    
//Creates a PSI_Coordinate
PSI_Coordinate* PSI_Photonics::Coordinate( 
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
    PSI_Coordinate_Photonics* coordinate = new PSI_Coordinate_Photonics(
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
    return coordinate;
}

int PSI_Photonics::GetLevel( PSI_Coordinate_Photonics* coordinate) const
{
    if ( coordinate == 0 ) {
	log_error("GetLevel failed, null coordinate supplied");
	return -1;
    }

    if ( !( IsTablesLoaded(1) || IsTablesLoaded(2) ) ) {
	log_error("GetLevel failed, no tables loaded");
	return -1;
    }

    int level = -1;
    
    if ( ( coordinate->trackType_ == 1 ) || 
	 ( coordinate->trackType_ == 2 ) ||
	 ( coordinate->trackType_ == 9 ) ||
	 ( coordinate->trackType_ == 10 ) ||
	 ( coordinate->trackType_ == 11 ) ||
	 ( coordinate->trackType_ == 12 ) ){
	//Shower (point like emmission)
	if ( IsTablesLoaded(1) ) {
	    level = 1;
	} else if ( IsTablesLoaded(2) ) {
	  //log_info("Assuming level 2 tables handles showers");
	    level = 2;
	}
    } else if ( coordinate->trackType_ == 0 ) {
	//Muon (extended emmission)
	if  ( IsTablesLoaded(2) ) {
	    level = 2;
	}
    } else {
	log_error("GetLevel failed, unknown particle type %d",coordinate->trackType_);
    }
    
    return level;
}


//Gets a  mean amplitude
bool PSI_Photonics::MeanAmplitude(
    double &amplitude,
    PSI_Coordinate* coordinate)
{
    //Initialize amplitude
    amplitude = -1;

    //Get and check coordinate
    PSI_Coordinate_Photonics* coord = 
	dynamic_cast<PSI_Coordinate_Photonics*>(coordinate);
    if ( coord == 0 ) {
	log_trace("MeanAmplitude failed, "
		  "did not recieve a Photonics coordinate");
	return false;
    }
    
    //Check the photonics level
    int level = GetLevel(coord);
    if ( level<0 ) {
	log_trace("There was a problem determining photonics level");
	return false;
    }
    
    //Calculate the photonics coordinate
    coord->CalculateCoordinate();
    
    //Get the amplitude
    if ( level == 1 ) {
      int type(1);
      if(coord->trackType_ >= 9 && 
	 coord->trackType_ <= 12) type = coord->trackType_;
	// Use EM shower tables for now (type=1):
	amplitude = get_photon_density(
	    coord->photonics_zsrc_,
	    coord->photonics_theta_,
	    type,
	    coord->photonics_l_,
	    coord->photonics_rho_,
	    coord->photonics_phi_*180/M_PI,
	    interpolationMode_,
	    &coord->tableId_,
	    &coord->lBin_,
	    &coord->rhoBin_,
	    &coord->phiBin_);
    } else if ( level == 2 ) {
	amplitude = 
	    photonics_cppio_obj_.get_level2_general_photon_density(
		coord->photonics_theta_,
		coord->photonics_rho_,
		coord->photonics_phi_,
		coord->photonics_l_,
		coord->photonics_stop_l_,
		coord->photonics_zsrc_,
		coord->photonics_stop_zsrc_,
		interpolationMode_,
		&coord->tableSetId_,
		&coord->tableId_,
		&coord->stopTableId_,
		&coord->lBin_,
		&coord->stopLBin_,
		&coord->rhoBin_,
		&coord->stopRhoBin_,
		&coord->phiBin_,
		&coord->stopPhiBin_ );	
    } else {
	log_trace("MeanAmplitude failed, unknown photonics level %d",level);
	return false;
    }

    if ( !std::isnormal(amplitude) ) {
	log_debug("Photonics returned NAN or INF "
		 "amplitude, setting amplitude to 0");
	amplitude = 0;
	return false;
    } else if ( amplitude < 0 ) {
	log_debug(
	    "Photonics returned negative amplitude %f, "
	    "(probably because reading outside table limits) "
	    "setting amplitude to 0 ", amplitude);
	amplitude = 0;
	return false;
    } else if ( amplitude == 0 ) {
	log_debug("Photonics amplitude was 0, "
		  "there will be no timing information availible");
	return false;
    }
    
    if ( coord->trackType_ < 3  && 
	 coord->trackEnergy_ >= 0 ) {
	//Scale amplitude with energy
	double energy = coord->trackEnergy_;

	//If energy is too low, use 1 GeV so light is well behaved
	if ( energy < 1 ) { 
	    log_trace(
		"Energy %f<1 GeV, "
		"using 1 GeV (lowest possible) for light factor "
		"calculation", 
		energy);
	    energy = 1;
	}

	double lightFactor = 
	    photonics_cppio_obj_.light(coord->trackType_, energy );

	log_debug("Scaling amplitude %f amplitude with light factor %f",
		  amplitude, lightFactor);
	amplitude = amplitude*lightFactor;
    } else {
	log_trace(
	    "MeanAmplitude, will not scale amplitude with light factor."
	    "Reason: Either track energy is <0 or track type unknown."
	    "track type = %d, track energy = %f",
	    coord->trackType_, coord->trackEnergy_);
    }
    
    log_debug("Got MeanAmplitude=%f from Photonics",amplitude);
    
    return true;
}

//Gets a  Time Delay
bool PSI_Photonics::TimeDelay(
    double &timeDelay,
    const double &random,
    PSI_Coordinate* coordinate)
{
    //Initialize time delay
    timeDelay = -1;
    
    //Get and check coordinate
    PSI_Coordinate_Photonics* coord = 
	dynamic_cast< PSI_Coordinate_Photonics* >(coordinate);   
    if ( coord == 0 ) {
	log_trace("TimeDelay failed, "
		  "did not recieve a Photonics coordinate");
	return false;
    }

    //Check the photonics level
    int level = GetLevel(coord);
    if ( level<0 ) {
	log_trace("TimeDelay failed, "
		  "There was a problem determining photonics level");
	return false;
    }

    //Get the photon time delay
    if ( level == 1 ) {
	timeDelay = get_photon_delay(
	    random,
	    coord->tableId_,
	    coord->lBin_,
	    coord->rhoBin_,
	    coord->phiBin_);	
    } else if ( level == 2 ) {	
	timeDelay = photonics_cppio_obj_.get_level2_general_delay(
	    random,
	    coord->tableSetId_,
	    coord->tableId_,
	    coord->stopTableId_,
	    coord->lBin_,
	    coord->stopLBin_,
	    coord->rhoBin_,
	    coord->stopRhoBin_,
	    coord->phiBin_,
	    coord->stopPhiBin_ );
    } else {
	log_trace("TimeDelay failed, unknown photonics level %d, level",level);
	return false;
    }

    if ( timeDelay==-1 ) {
	log_debug(
	    "Photonics returned -1 timedelay, "
	    "perhaps something went wrong, "
	    "(PSI will not take any action on this condition)"); 
    }

    log_debug("Got time delay %f from Photonics", timeDelay);

    return true;
}

//Gets a  Probability
bool PSI_Photonics::Probability(
    double &probability,
    const double &timeDelay,
    PSI_Coordinate* coordinate)
{
    //Initialize probability
    probability = -1;

    //Get and check coordinate
    PSI_Coordinate_Photonics* coord = 
	dynamic_cast<PSI_Coordinate_Photonics*>(coordinate);
    if ( coord == 0 ) {
	log_trace("Probability failed, "
		  "did not recieve a Photonics coordinate");
	return false;
    }
    
    //Check the photonics level
    int level = GetLevel(coord);
    if ( level<0 ) {
	log_trace("Probability failed, there was a problem "
		  "determining photonics level");
	return false;
    }
    
    //Get the photon time delay
    if ( level == 1 ) {
	probability = get_delay_prob(
	    timeDelay,
	    coord->tableId_,
	    coord->lBin_,
	    coord->rhoBin_,
	    coord->phiBin_);	
    } else if (level == 2 ) {	
	probability = photonics_cppio_obj_.get_level2_general_delay_prob(
	    timeDelay,
	    coord->tableSetId_,
	    coord->tableId_,
	    coord->stopTableId_,
	    coord->lBin_,
	    coord->stopLBin_,
	    coord->rhoBin_,
	    coord->stopRhoBin_,
	    coord->phiBin_,
	    coord->stopPhiBin_ );
    } else {
	log_trace("Probability failed, unknown photonics level %d",level);
	return false;
    }

    
    //Check probability
    if ( probability < 0 ) {	
	log_debug(
	    "Photonics returned probability %f < 0, "
	    "setting probability to 0",probability);
	probability = 0;
    } else if ( probability > 1 ) {	
	log_debug(
	    "Photonics returned probability %f > 1, "
	    "setting probability to 1",probability);
	probability = 1;
    } else {
	log_debug("Got probability %f from Photonics", probability);
    }

    return true;
}

//Change the tables loaded state for level 2 photonics
bool PSI_Photonics::SetTablesLoaded( 
    const int &level,
    const bool &state )
{
    if (level==1) {
	level1TablesLoaded_ = state; 
    } else if (level==2){
	level2TablesLoaded_ = state; 
    } else {
	log_error(
	    "SetTablesLoaded failed! "
	    "Invalid photonics level %d, select 1 or 2",
	    level);
	return false;
    }
    string strstate="false";
    if (state) strstate="true";
    log_debug("Changed level 1 tables loaded status to %s", 
	      strstate.c_str());
    return true;
}
