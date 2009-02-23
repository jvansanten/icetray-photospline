/**
 * \file PSI_Coordinate_Photonics_Tests.cxx
 * \brief test for PSI_Coordinate_Photonics
 *  
 * \author Thomas Burgess 
 *
 * (c) the IceCube Collaboration
 */

#include <I3Test.h>
#include "PSInterface/PSI_Coordinate.h"
#include "PSInterface/PSI_Coordinate_Photonics.h"
#include<iostream>
using namespace std;
TEST_GROUP(PSI_Coordinate_Photonics);
//Return true if a=b within accuracy

inline bool testDouble( double a, double b )
{
    //Accuracy foe double should be ok upto ~10^15
    const double accuracy=1e7;
    return (round(a*accuracy)==round(b*accuracy));
}

TEST(Coordinate_Photonics)
{
    //Checks that base class Get functions return correct values
    //An error here means that the Coordinate constructor was 
    //called the wrong way in Coordinate_Photonics


    const double opticalModuleX = 0;
    const double opticalModuleY = 0;
    const double opticalModuleZ = 0;
    const double opticalModuleOrientation = -1;
    const double trackX = -10;
    const double trackY = 0;
    const double trackZ = 10;
    const double trackTheta = 0;
    const double trackPhi = 270;
    const double trackLength = 10;
    const double trackEnergy = 50;
    const int    trackType = 0; 

    PSI_Coordinate_Photonics coordinate_Photonics(
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
  
    ENSURE((opticalModuleX==coordinate_Photonics.GetOpticalModuleX()),
	   "There was a problem with Optical Module X");
    ENSURE((opticalModuleY==coordinate_Photonics.GetOpticalModuleY()),
	   "There was a problem with Optical Module Y");
    ENSURE((opticalModuleZ==coordinate_Photonics.GetOpticalModuleZ()),
	   "There was a problem with Optical Module Z");
    ENSURE((opticalModuleOrientation==
	    coordinate_Photonics.GetOpticalModuleOrientation()),
	   "There was a problem with Optical Module");
    ENSURE((trackX==coordinate_Photonics.GetTrackX()),
	   "There was a problem with Track X");
    ENSURE((trackY==coordinate_Photonics.GetTrackY()),
	   "There was a problem with Track Y");
    ENSURE((trackZ==coordinate_Photonics.GetTrackZ()),
	   "There was a problem with Track Z");
    ENSURE((trackTheta==coordinate_Photonics.GetTrackTheta()),
	   "There was a problem with Track Theta");
    ENSURE((trackPhi==coordinate_Photonics.GetTrackPhi()),
	   "There was a problem with Track Phi");
    ENSURE((trackLength==coordinate_Photonics.GetTrackLength()),
	   "There was a problem with Track Length");
    ENSURE((trackEnergy==coordinate_Photonics.GetTrackEnergy()),
	   "There was a problem with Track Energy");
    ENSURE((trackType==coordinate_Photonics.GetTrackType()),
	   "There was a problem with Track Type");
}

TEST(PPerp) 
{
    //OM and Track variables
    double OM_x, OM_y, OM_z, OM_o;
    double TR_x, TR_y, TR_z, TR_theta, TR_phi, TR_l, TR_e;
    int TR_type;

    //Downlooking OM at origin
    OM_x=0;OM_y=0; OM_z=0; OM_o=-1;

    // 1 GeV ems shower along x with y=0
    TR_x=10; TR_y=10; TR_z=10; TR_l=0; TR_e=1; TR_type=1;
    
    //Direction vector
    double x,y,z;

    //Perpendicular vector
    double px,py,pz;
    
    //PSI_Coordinate
    PSI_Coordinate_Photonics coord;

    //along -z
    TR_theta=0; TR_phi=0;
    coord=PSI_Coordinate_Photonics(
	OM_x, OM_y, OM_z, OM_o,
	TR_x, TR_y, TR_z, TR_theta, TR_phi, TR_l, TR_e, TR_type);    
    coord.CalcTrackDirection(x,y,z);   
    coord.CalculatePPerp(px,py,pz);   

    //Test that pperp is perpendicular to track direction
    ENSURE( testDouble(x*px,0),"PPerp X wrong, not perpendicular to track direction");
    ENSURE( testDouble(y*py,0),"PPerp Y wrong, not perpendicular to track direction");
    ENSURE( testDouble(z*pz,0),"PPerp Z wrong, not perpendicular to track direction");
}
