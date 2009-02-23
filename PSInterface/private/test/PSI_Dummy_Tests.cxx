/**
 * \file PSI_Dummy_Tests.cxx
 * \brief test for PSI_Dummy
 *  
 *@author Thomas Burgess and Alex Olivas
 *
 * (c) the IceCube Collaboration
 */

#include <I3Test.h>

#include "PSInterface/PSI_Dummy.h"
#include "PSInterface/PSI_Coordinate.h"

TEST_GROUP(PSI_Dummy);

TEST(PSI_Dummy)
{
    //Creates, copies and deletes PSI_Dummy, should always work

    PSI_Dummy psi;        //Constructor
    PSI_Dummy psi2(psi);  //Copy constructor   
    psi=psi2;             //Assignment
    PSI_Dummy* psiptr = new PSI_Dummy();
    delete psiptr;        //Destructor
}


TEST(MakeCoordinate)
{
    // Makes a coordinate and tests if it is not = 0

    PSI_Dummy psi;

    //Initialize OM and Track variables
    double OM_x=0, OM_y=0, OM_z=0, OM_o=-1;
    double TR_x=10, TR_y=10, TR_z=10, TR_theta=90, TR_phi=0, TR_l=0, TR_e=10;
    int TR_type=1;
    
    //PSI_Coordinate
    PSI_Coordinate* coord = 0;
    
    coord=psi.MakeCoordinate(
	OM_x, OM_y, OM_z, OM_o,
	TR_x, TR_y, TR_z, TR_theta, TR_phi, TR_l, TR_e, TR_type);    

    ENSURE(coord!=0, "There was a problem making a coordinate.");

    delete coord;
}

TEST(GetMeanAmplitude)
{
    //Makes a psinterface, 
    // - makes a coordinate
    // - ensures coordinate!=0 
    // - gets amplitude 
    // - ensures that amplitude call was successful
    
    PSI_Dummy psi;

    //OM and Track variables
    double OM_x, OM_y, OM_z, OM_o;
    double TR_x, TR_y, TR_z, TR_theta, TR_phi, TR_l, TR_e;
    int TR_type;

    //Downlooking OM at origin
    OM_x=0;OM_y=0; OM_z=0; OM_o=-1;

    // 1 GeV ems shower at (10,10,10)
    TR_x=10; TR_y=10; TR_z=10; TR_l=0; TR_e=1; TR_type=1;
    
    //PSI_Coordinate
    PSI_Coordinate* coord = 0;

    //along -z
    TR_theta=0; TR_phi=0;
    coord=psi.MakeCoordinate(
	OM_x, OM_y, OM_z, OM_o,
	TR_x, TR_y, TR_z, TR_theta, TR_phi, TR_l, TR_e, TR_type);    

    ENSURE(coord!=0, "There was a problem making a coordinate."); 

    double amplitude = -1;
    bool getMeanAmplitude = psi.GetMeanAmplitude(amplitude,coord);
    ENSURE(getMeanAmplitude,"There was a problem getting a mean amplitude");
    
    delete coord;
}

TEST(GetTimeDelay)
{
    //Makes a psinterface, 
    // - makes a coordinate
    // - ensures coordinate!=0 
    // - gets amplitude 
    // - ensures that amplitude call was successful
    // - gets timedelay
    // - ensures that timedelay call was successful

    PSI_Dummy psi;

    //OM and Track variables
    double OM_x, OM_y, OM_z, OM_o;
    double TR_x, TR_y, TR_z, TR_theta, TR_phi, TR_l, TR_e;
    int TR_type;

    //Downlooking OM at origin
    OM_x=0;OM_y=0; OM_z=0; OM_o=-1;

    // 1 GeV ems shower at (10,10,10)
    TR_x=10; TR_y=10; TR_z=10; TR_l=0; TR_e=1; TR_type=1;
    
    //PSI_Coordinate
    PSI_Coordinate* coord = 0;

    //along -z
    TR_theta=0; TR_phi=0;
    coord=psi.MakeCoordinate(
	OM_x, OM_y, OM_z, OM_o,
	TR_x, TR_y, TR_z, TR_theta, TR_phi, TR_l, TR_e, TR_type);    

    ENSURE(coord!=0, "There was a problem making a coordinate."); 

    double amplitude = -1;
    bool getMeanAmplitude = psi.GetMeanAmplitude(amplitude,coord);
    ENSURE(getMeanAmplitude,"There was a problem getting a mean amplitude");
    
    double timeDelay=-1;
    double random = 0.5;
    bool getTimeDelay = psi.GetTimeDelay(timeDelay,random,coord);
    ENSURE(getTimeDelay,"There was a problem getting time delay");
    
    delete coord;
}

TEST(GetProbability)
{
        //Makes a psinterface, 
    // - makes a coordinate
    // - ensures coordinate!=0 
    // - gets amplitude 
    // - ensures that amplitude call was successful
    // - gets timedelay
    // - ensures that timedelay call was successful

    PSI_Dummy psi;

    //OM and Track variables
    double OM_x, OM_y, OM_z, OM_o;
    double TR_x, TR_y, TR_z, TR_theta, TR_phi, TR_l, TR_e;
    int TR_type;

    //Downlooking OM at origin
    OM_x=0;OM_y=0; OM_z=0; OM_o=-1;

    // 1 GeV ems shower at (10,10,10)
    TR_x=10; TR_y=10; TR_z=10; TR_l=0; TR_e=1; TR_type=1;
    
    //PSI_Coordinate
    PSI_Coordinate* coord = 0;

    //along -z
    TR_theta=0; TR_phi=0;
    coord=psi.MakeCoordinate(
	OM_x, OM_y, OM_z, OM_o,
	TR_x, TR_y, TR_z, TR_theta, TR_phi, TR_l, TR_e, TR_type);    

    ENSURE(coord!=0, "There was a problem making a coordinate."); 

    double amplitude = -1;
    bool getMeanAmplitude = psi.GetMeanAmplitude(amplitude,coord);
    ENSURE(getMeanAmplitude,"There was a problem getting a mean amplitude");
    
    double probability = -1;
    double timeDelay=10;
    bool getProbability = psi.GetProbability(probability,timeDelay,coord);
    ENSURE(getProbability,"There was a problem getting hit probability");

    delete coord;
}

