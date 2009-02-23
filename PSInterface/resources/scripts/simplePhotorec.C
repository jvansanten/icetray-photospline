/**
 *@file
 *@brief Simple example for PSI_PhotonicsPhotorec
 *
 *@author Thomas Burgess
 * (c) the IceCube Collaboration
 */

/**
 * PSI_PhotonicsPhotorec simple test
 *
 * Run with ".x simplePhotorec.C" in root or  "root -q -l -b simplePhotorec.C"
 *  from commandline
 *
 *@author Thomas Burgess
 */
void simplePhotorec() {
    //Load the psi library
    if ( gSystem->Load("libPSInterface.so") < 0 ) {
        cerr << "Failed to load PSInterface, exiting" << endl;
	return;
    }    

    //Initialize
    PSI_PhotonicsPhotorec rec;
    if ( !rec.LoadTables(2,"level2_muon_photorec.list") ) {
	cerr << "Failed to load photonics tables, exiting" << endl;
	return;
    }
    
    //Create a photonics coordinate for a 
    //downlooking optical module at (0,0,0) and a
    //100 GeV muon track at (0,0,-10) with angles (45,180) for 10m
    double OM_x=0, OM_y=0, OM_z=0, OM_o=-1;
    double TR_x=0, TR_y=0, TR_z=-10, TR_th=45, TR_ph=180,
	TR_length=10, TR_energy=100;
    int TR_type = 0;

    PSI_Coordinate_Photonics* coordinate = 
	dynamic_cast<PSI_Coordinate_Photonics*>(
	    rec.MakeCoordinate(
		OM_x,OM_y,OM_z,OM_o,
		TR_x,TR_y,TR_z,TR_th,TR_ph,
		TR_length,TR_energy,TR_type));

    //Check coordinate
    if ( coordinate == 0 ) {
	cerr << "Failed to make coordinate" << endl;
	return;
    }
    
    //Define timedelay for query
    double timedelay = 25;
    //Define amplitude and probability to fetch
    double amplitude = -1;
    double probability = -1;
    
    //Make photorec query and print result
    if ( !rec.Get(coordinate,timedelay,amplitude,probability) ) {
	cerr << "Failed to get photorec result" << endl;
	return;
    }

    cout << "Amplitude="<<amplitude<<" pdf="<<probability << endl;
}
