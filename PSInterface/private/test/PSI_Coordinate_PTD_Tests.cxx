/**
 * \file PSI_Coordinate_PTD_Tests.cxx
 * \test for PSI_Coordinate_Ptd
 *  
 *@author Thomas Burgess 
 *
 * (c) the IceCube Collaboration
 */

#include <I3Test.h>
#include "PSInterface/PSI_Coordinate.h"
#include "PSInterface/PSI_Coordinate_PTD.h"
#include<iostream>
using namespace std;
TEST_GROUP(PSI_Coordinate_PTD);

TEST(Coordinate_PTD)
{
    //Checks that base class Get functions return correct values
    //An error here means that the Coordinate constructor was 
    //called the wrong way in Coordinate_PTD

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

    PSI_Coordinate_PTD coordinate_PTD(
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
  
    ENSURE((opticalModuleX==coordinate_PTD.GetOpticalModuleX()),
	   "There was a problem with Optical Module X");
    ENSURE((opticalModuleY==coordinate_PTD.GetOpticalModuleY()),
	   "There was a problem with Optical Module Y");
    ENSURE((opticalModuleZ==coordinate_PTD.GetOpticalModuleZ()),
	   "There was a problem with Optical Module Z");
    ENSURE((opticalModuleOrientation==
	    coordinate_PTD.GetOpticalModuleOrientation()),
	   "There was a problem with Optical Module");
    ENSURE((trackX==coordinate_PTD.GetTrackX()),
	   "There was a problem with Track X");
    ENSURE((trackY==coordinate_PTD.GetTrackY()),
	   "There was a problem with Track Y");
    ENSURE((trackZ==coordinate_PTD.GetTrackZ()),
	   "There was a problem with Track Z");
    ENSURE((trackTheta==coordinate_PTD.GetTrackTheta()),
	   "There was a problem with Track Theta");
    ENSURE((trackPhi==coordinate_PTD.GetTrackPhi()),
	   "There was a problem with Track Phi");
    ENSURE((trackLength==coordinate_PTD.GetTrackLength()),
	   "There was a problem with Track Length");
    ENSURE((trackEnergy==coordinate_PTD.GetTrackEnergy()),
	   "There was a problem with Track Energy");
    ENSURE((trackType==coordinate_PTD.GetTrackType()),
	   "There was a problem with Track Type");
}
