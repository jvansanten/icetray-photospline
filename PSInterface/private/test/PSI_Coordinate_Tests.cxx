/**
 * \file PSI_Coordinate_Tests.cxx
 * \test for PSI_Coordinate
 *  
 *@author Thomas Burgess 
 *
 * (c) the IceCube Collaboration
 */

#include <I3Test.h>
#include "PSInterface/PSI_Coordinate.h"
#include "PSInterface/PSI_Dummy.h"
#include "boost/filesystem/path.hpp"
#include "boost/filesystem/operations.hpp"
using boost::filesystem::exists;
using boost::filesystem::path;
using boost::filesystem::current_path;
using boost::filesystem::remove;

TEST_GROUP(PSI_Coordinate);

//Return true if a=b within accuracy
inline bool testDouble( double a, double b )
{
    //Accuracy foe double should be ok upto ~10^15
    const double accuracy=1e7;
    return (round(a*accuracy)==round(b*accuracy));
}

TEST(Coordinate)
{
    //Checks that Get funcions return what coordinate was contrutcted with
    //An error here means that the intialization of Coordinate or 
    //a Get function is wrong

    const double opticalModuleX = 0.0;
    const double opticalModuleY = 0.0;
    const double opticalModuleZ = 0.0;
    const double opticalModuleOrientation = -1.0;
    const double trackX = -10.0;
    const double trackY = 0.0;
    const double trackZ = 10.0;
    const double trackTheta = 0.0;
    const double trackPhi = 270.0;
    const double trackLength = 10.0;
    const double trackEnergy = 50.0;
    const int    trackType = 0; 

    PSI_Coordinate coordinate(
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
  
    ENSURE((opticalModuleX==coordinate.GetOpticalModuleX()),
	   "There was a problem with Optical Module X");
    ENSURE((opticalModuleY==coordinate.GetOpticalModuleY()),
	   "There was a problem with Optical Module Y");
    ENSURE((opticalModuleZ==coordinate.GetOpticalModuleZ()),
	   "There was a problem with Optical Module Z");
    ENSURE((opticalModuleOrientation==
	    coordinate.GetOpticalModuleOrientation()),
	   "There was a problem with Optical Module");
    ENSURE((trackX==coordinate.GetTrackX()),
	   "There was a problem with Track X");
    ENSURE((trackY==coordinate.GetTrackY()),
	   "There was a problem with Track Y");
    ENSURE((trackZ==coordinate.GetTrackZ()),
	   "There was a problem with Track Z");
    ENSURE((trackTheta==coordinate.GetTrackTheta()),
	   "There was a problem with Track Theta");
    ENSURE((trackPhi==coordinate.GetTrackPhi()),
	   "There was a problem with Track Phi");
    ENSURE((trackLength==coordinate.GetTrackLength()),
	   "There was a problem with Track Length");
    ENSURE((trackEnergy==coordinate.GetTrackEnergy()),
	   "There was a problem with Track Energy");
    ENSURE((trackType==coordinate.GetTrackType()),
	   "There was a problem with Track Type");
}

TEST(CalcTrackDirection) 
{  
    // The direction vector (x,y,z) for different track angles (theta,phi)
    // (using an OM at origin, a track vertex at (10,10,10))
    //
    // n |theta| phi | x | y | z |
    // --+-----+-----+---+---+---|
    // 1 |   0 |   0 | 0 | 0 |-1 | along -z (down)
    // 2 |  90 |   0 |-1 | 0 | 0 | along -x (horizontal)
    // 3 |  90 |  90 | 0 |-1 | 0 | along -y
    // 4 |  90 | 180 | 1 | 0 | 0 | along  x
    // 5 |  90 | 270 | 0 | 1 | 0 | along  y 
    // 6 | 180 |   0 | 0 | 0 | 1 | along  z (up)

    //OM and Track variables
    double OM_x, OM_y, OM_z, OM_o;
    double TR_x, TR_y, TR_z, TR_theta, TR_phi, TR_l, TR_e;
    int TR_type;

    //Downlooking OM at origin
    OM_x=0.0;OM_y=0.0; OM_z=0.0; OM_o=-1.0;

    // 1 GeV ems shower at (10,10,10)
    TR_x=10.0; TR_y=10.0; TR_z=10.0; TR_l=0.0; TR_e=1.0; TR_type=1;
    
    //Direction vector
    double x,y,z;
    
    //PSI_Coordinate
    PSI_Coordinate coord;

    //along -z
    TR_theta=0.0; TR_phi=0.0;
    coord=PSI_Coordinate(
	OM_x, OM_y, OM_z, OM_o,
	TR_x, TR_y, TR_z, TR_theta, TR_phi, TR_l, TR_e, TR_type);    
    coord.CalcTrackDirection(x,y,z);   
    ENSURE( testDouble(x,0.0),"X wrong for (theta,phi)=(0,0)");
    ENSURE( testDouble(y,0.0),"Y wrong for (theta,phi)=(0,0)");
    ENSURE( testDouble(z,-1.0),"Z wrong for (theta,phi)=(0,0)");

    //Horizontal along -x
    TR_theta=90.0; TR_phi=0.0;
    coord=PSI_Coordinate(
	OM_x, OM_y, OM_z, OM_o,
	TR_x, TR_y, TR_z, TR_theta, TR_phi, TR_l, TR_e, TR_type);    
    coord.CalcTrackDirection(x,y,z);   
    ENSURE( testDouble(x,-1.0),"X wrong for (theta,phi)=(90,0)");
    ENSURE( testDouble(y,0.0),"Y wrong for (theta,phi)=(90,0)");
    ENSURE( testDouble(z,0.0),"Z wrong for (theta,phi)=(90,0)");
    
    //Horizontal along -y
    TR_theta=90.0; TR_phi=90.0;
    coord=PSI_Coordinate(
	OM_x, OM_y, OM_z, OM_o,
	TR_x, TR_y, TR_z, TR_theta, TR_phi, TR_l, TR_e, TR_type);    
    coord.CalcTrackDirection(x,y,z);   
    ENSURE( testDouble(x,0.0),"X wrong for (theta,phi)=(90,90)");
    ENSURE( testDouble(y,-1.0),"Y wrong for (theta,phi)=(90,90)");
    ENSURE( testDouble(z,0.0),"Z wrong for (theta,phi)=(90,90)");

    //Horizontal along x
    TR_theta=90.0; TR_phi=180.0;
    coord=PSI_Coordinate(
	OM_x, OM_y, OM_z, OM_o,
	TR_x, TR_y, TR_z, TR_theta, TR_phi, TR_l, TR_e, TR_type);    
    coord.CalcTrackDirection(x,y,z);   
    ENSURE( testDouble(x,1.0),"X wrong for (theta,phi)=(90,180)");
    ENSURE( testDouble(y,0.0),"Y wrong for (theta,phi)=(90,180)");
    ENSURE( testDouble(z,0.0),"Z wrong for (theta,phi)=(90,180)");
 
    //Horizontal along y
    TR_theta=90.0; TR_phi=270.0;
    coord=PSI_Coordinate(
	OM_x, OM_y, OM_z, OM_o,
	TR_x, TR_y, TR_z, TR_theta, TR_phi, TR_l, TR_e, TR_type);    
    coord.CalcTrackDirection(x,y,z);   
    ENSURE( testDouble(x,0.0),"X wrong for (theta,phi)=(90,270)");
    ENSURE( testDouble(y,1.0),"Y wrong for (theta,phi)=(90,270)");
    ENSURE( testDouble(z,0.0),"Z wrong for (theta,phi)=(90,270)");

    //Vertical along z (up)
    TR_theta=180.0; TR_phi=0.0;
    coord=PSI_Coordinate(
	OM_x, OM_y, OM_z, OM_o,
	TR_x, TR_y, TR_z, TR_theta, TR_phi, TR_l, TR_e, TR_type);    
    coord.CalcTrackDirection(x,y,z);   
    ENSURE( testDouble(x,0.0),"X wrong for (theta,phi)=(180,0)");
    ENSURE( testDouble(y,0.0),"Y wrong for (theta,phi)=(180,0)");
    ENSURE( testDouble(z,1.0),"Z wrong for (theta,phi)=(180,0)");
}

TEST(CalcDistance) 
{
    PSI_Coordinate coord(
	0.0,0.0,0.0,-1.0,-10.0,0.0,10.0,90.0,0.0,0.0,10.0,0);

    //OM at (0,0,0), Track at (-10,0,10) hence
    //distance between track and OM is 
    //dist = sqrt((-10.0)^2+0^2+(10.0)^2) = 10.0*sqrt(2.0)

    ENSURE( testDouble(coord.CalcDistance(),10.0*sqrt(2.0)),
	    "Distance calculation wrong");
}

TEST(CalcRho) 
{
    PSI_Coordinate coord(
	0.0,0.0,0.0,-1.0,-10.0,0.0,10.0,90.0,0.0,0.0,10.0,0);

    //OM at (0,0,0), Track at (-10,0,10) heading(0,0) hence
    //shortest perpendicular distance between track and OM is 10 ,
    
    ENSURE( testDouble(coord.CalcRho(),10.0),
	    "Rho calculation wrong");
}

TEST(CalcRhoVector) 
{
    PSI_Coordinate coord(0.0,0.0,0.0,-1.0,-10.0,0.0,10.0,90.0,0.0,0.0,10.0,0);

    //OM at (0,0,0), Track at (-10,0,10) heading(0,0) hence
    //vector shortest perpendicular distance between track and OM is
    //0,0,-10
    double x,y,z;
    coord.CalcRhoVector(x,y,z);
    
    ENSURE( testDouble(x,0), "rho X wrong");
    ENSURE( testDouble(y,0), "rho Y wrong");
    ENSURE( testDouble(z,-10), "rho Z wrong");
}

TEST(CalcOnTrackDistance)
{  
    // OnTrackDistance l, 
    // for different track theta,x and z given a downlooking OM at (0,0,0)
    //                (if x the y,z=10)
    // n |theta| phi | xyz | l |
    // --+-----+-----+-----+---+---|
    //   |   0 |   0 |z -10|-10| along -z (vertical down)
    //   |   0 |   0 |z   0|  0| 
    //   |   0 |   0 |z  10| 10| 
    //   |  90 |   0 |x -10|-10| along -x (horizontal)
    //   |  90 |   0 |x   0|  0| 
    //   |  90 |   0 |x  10| 10| 
    //   |  90 |  90 |y -10|-10| along -y
    //   |  90 |  90 |y   0|  0|
    //   |  90 |  90 |y  10| 10|
    //   |  90 | 180 |x -10| 10| along +x
    //   |  90 | 180 |x   0|  0|             
    //   |  90 | 180 |x  10|-10|             
    //   |  90 | 270 |y -10| 10| along +y
    //   |  90 | 270 |y   0|  0|
    //   |  90 | 270 |y  10|-10|
    //   | 180 |   0 |z -10| 10| along +z (vertical up)
    //   | 180 |   0 |z   0|  0| 
    //   | 180 |   0 |z  10|-10| 


    //OM and Track variables
    double OM_x, OM_y, OM_z, OM_o;
    double TR_x, TR_y, TR_z, TR_theta, TR_phi, TR_l, TR_e;
    int TR_type;
    //Downlooking OM at origin
    OM_x=0.0; OM_y=0.0; OM_z=0.0; OM_o=-1.0;

    // 1 GeV ems shower along x with y=0
    TR_l=0.0; TR_e=1.0; TR_type=1;

    PSI_Coordinate coord;

    //Down
    TR_theta=0.0; TR_phi=0.0; TR_x=10.0; TR_y=10.0; TR_z=-10.0; 
    coord=PSI_Coordinate(
	OM_x, OM_y, OM_z, OM_o,
	TR_x, TR_y, TR_z, TR_theta, TR_phi, TR_l, TR_e, TR_type);    
    ENSURE( testDouble( coord.CalcOnTrackDistance(), -10.0), 
	    "OnTrackDistance wrong (theta,phi,x,y,z)=(0,0,10,10,-10)");

    TR_theta=0.0; TR_phi=0.0; TR_x=10.0; TR_y=10.0; TR_z=0.0; 
    coord=PSI_Coordinate(
	OM_x, OM_y, OM_z, OM_o,
	TR_x, TR_y, TR_z, TR_theta, TR_phi, TR_l, TR_e, TR_type);    
    ENSURE( testDouble( coord.CalcOnTrackDistance(), 0.0), 
	    "OnTrackDistance wrong (theta,phi,x,y,z)=(0,0,10,10,0)");

    TR_theta=0.0; TR_phi=0.0; TR_x=10.0; TR_y=10.0; TR_z=10.0; 
    coord=PSI_Coordinate(
	OM_x, OM_y, OM_z, OM_o,
	TR_x, TR_y, TR_z, TR_theta, TR_phi, TR_l, TR_e, TR_type);    
    ENSURE( testDouble( coord.CalcOnTrackDistance(), 10.0), 
	    "OnTrackDistance wrong (theta,phi,x,y,z)=(0,0,10,10,10)");

    //Horizontal along -x
    TR_theta=90.0; TR_phi=0.0; TR_x=-10.0; TR_y=10.0; TR_z=10.0; 
    coord=PSI_Coordinate(
	OM_x, OM_y, OM_z, OM_o,
	TR_x, TR_y, TR_z, TR_theta, TR_phi, TR_l, TR_e, TR_type);    
    ENSURE( testDouble( coord.CalcOnTrackDistance(), -10.0), 
	    "OnTrackDistance wrong (theta,phi,x,y,z)=(90,0,-10,10,10)");

    TR_theta=90.0; TR_phi=0.0; TR_x=0.0; TR_y=10.0; TR_z=10.0; 
    coord=PSI_Coordinate(
	OM_x, OM_y, OM_z, OM_o,
	TR_x, TR_y, TR_z, TR_theta, TR_phi, TR_l, TR_e, TR_type);    
    ENSURE( testDouble( coord.CalcOnTrackDistance(), 0.0), 
	    "OnTrackDistance wrong (theta,phi,x,y,z)=(90,0,0,10,10)");

    TR_theta=90.0; TR_phi=0.0; TR_x=10.0; TR_y=10.0; TR_z=10.0; 
    coord=PSI_Coordinate(
	OM_x, OM_y, OM_z, OM_o,
	TR_x, TR_y, TR_z, TR_theta, TR_phi, TR_l, TR_e, TR_type);    
    ENSURE( testDouble( coord.CalcOnTrackDistance(), 10.0), 
	    "OnTrackDistance wrong (theta,phi,x,y,z)=(90,0,10,10,10)");

    //Horizontal along -y
    TR_theta=90.0; TR_phi=90.0; TR_x=10.0; TR_y=-10.0; TR_z=10.0; 
    coord=PSI_Coordinate(
	OM_x, OM_y, OM_z, OM_o,
	TR_x, TR_y, TR_z, TR_theta, TR_phi, TR_l, TR_e, TR_type);    
    ENSURE( testDouble( coord.CalcOnTrackDistance(), -10.0), 
	    "OnTrackDistance wrong (theta,phi,x,y,z)=(90,90,10,-10,10)");

    TR_theta=90.0; TR_phi=90.0; TR_x=10.0; TR_y=0.0; TR_z=10.0; 
    coord=PSI_Coordinate(
	OM_x, OM_y, OM_z, OM_o,
	TR_x, TR_y, TR_z, TR_theta, TR_phi, TR_l, TR_e, TR_type);    
    ENSURE( testDouble( coord.CalcOnTrackDistance(), 0.0), 
	    "OnTrackDistance wrong (theta,phi,x,y,z)=(90,90,10,0,10)");

    TR_theta=90.0; TR_phi=90.0; TR_x=10.0; TR_y=10.0; TR_z=10.0; 
    coord=PSI_Coordinate(
	OM_x, OM_y, OM_z, OM_o,
	TR_x, TR_y, TR_z, TR_theta, TR_phi, TR_l, TR_e, TR_type);    
    ENSURE( testDouble( coord.CalcOnTrackDistance(), 10.0), 
	    "OnTrackDistance wrong (theta,phi,x,y,z)=(90,90,10,10,10)");


    //Horizontal along x
    TR_theta=90.0; TR_phi=180.0; TR_x=-10.0; TR_y=10.0; TR_z=10.0; 
    coord=PSI_Coordinate(
	OM_x, OM_y, OM_z, OM_o,
	TR_x, TR_y, TR_z, TR_theta, TR_phi, TR_l, TR_e, TR_type);    
    ENSURE( testDouble( coord.CalcOnTrackDistance(), 10.0), 
	    "OnTrackDistance wrong (theta,phi,x,y,z)=(90,180,-10,10,10)");

    TR_theta=90.0; TR_phi=180.0; TR_x=0.0; TR_y=10.0; TR_z=10.0; 
    coord=PSI_Coordinate(
	OM_x, OM_y, OM_z, OM_o,
	TR_x, TR_y, TR_z, TR_theta, TR_phi, TR_l, TR_e, TR_type);    
    ENSURE( testDouble( coord.CalcOnTrackDistance(), 0.0), 
	    "OnTrackDistance wrong (theta,phi,x,y,z)=(90,180,0,10,10)");

    TR_theta=90.0; TR_phi=180.0; TR_x=10.0; TR_y=10.0; TR_z=10.0; 
    coord=PSI_Coordinate(
	OM_x, OM_y, OM_z, OM_o,
	TR_x, TR_y, TR_z, TR_theta, TR_phi, TR_l, TR_e, TR_type);    
    ENSURE( testDouble( coord.CalcOnTrackDistance(), -10.0), 
	    "OnTrackDistance wrong (theta,phi,x,y,z)=(90,180,10,10,10)");

    //Horizontal along -y
    TR_theta=90.0; TR_phi=270.0; TR_x=10.0; TR_y=-10.0; TR_z=10.0; 
    coord=PSI_Coordinate(
	OM_x, OM_y, OM_z, OM_o,
	TR_x, TR_y, TR_z, TR_theta, TR_phi, TR_l, TR_e, TR_type);    
    ENSURE( testDouble( coord.CalcOnTrackDistance(), 10.0), 
	    "OnTrackDistance wrong (theta,phi,x,y,z)=(90,270,10,-10,10)");

    TR_theta=90.0; TR_phi=270.0; TR_x=10.0; TR_y=0.0; TR_z=10.0; 
    coord=PSI_Coordinate(
	OM_x, OM_y, OM_z, OM_o,
	TR_x, TR_y, TR_z, TR_theta, TR_phi, TR_l, TR_e, TR_type);    
    ENSURE( testDouble( coord.CalcOnTrackDistance(), 0.0), 
	    "OnTrackDistance wrong (theta,phi,x,y,z)=(90,270,10,0,10)");

    TR_theta=90.0; TR_phi=270.0; TR_x=10.0; TR_y=10.0; TR_z=10.0; 
    coord=PSI_Coordinate(
	OM_x, OM_y, OM_z, OM_o,
	TR_x, TR_y, TR_z, TR_theta, TR_phi, TR_l, TR_e, TR_type);    
    ENSURE( testDouble( coord.CalcOnTrackDistance(), -10.0), 
	    "OnTrackDistance wrong (theta,phi,x,y,z)=(90,270,10,10,10)");

    //Up
    TR_theta=180.0; TR_phi=0.0; TR_x=10.0; TR_y=10.0; TR_z=-10.0; 
    coord=PSI_Coordinate(
	OM_x, OM_y, OM_z, OM_o,
	TR_x, TR_y, TR_z, TR_theta, TR_phi, TR_l, TR_e, TR_type);    
    ENSURE( testDouble( coord.CalcOnTrackDistance(), 10.0), 
	    "OnTrackDistance wrong (theta,phi,x,y,z)=(180,0,10,10,-10)");

    TR_theta=180.0; TR_phi=0.0; TR_x=10.0; TR_y=10.0; TR_z=0.0; 
    coord=PSI_Coordinate(
	OM_x, OM_y, OM_z, OM_o,
	TR_x, TR_y, TR_z, TR_theta, TR_phi, TR_l, TR_e, TR_type);    
    ENSURE( testDouble( coord.CalcOnTrackDistance(), 0.0), 
	    "OnTrackDistance wrong (theta,phi,x,y,z)=(180,0,10,10,0)");

    TR_theta=180.0; TR_phi=0.0; TR_x=10.0; TR_y=10.0; TR_z=10.0; 
    coord=PSI_Coordinate(
	OM_x, OM_y, OM_z, OM_o,
	TR_x, TR_y, TR_z, TR_theta, TR_phi, TR_l, TR_e, TR_type);    
    ENSURE( testDouble( coord.CalcOnTrackDistance(), -10.0), 
	    "OnTrackDistance wrong (theta,phi,x,y,z)=(180,0,10,10,10)");
}


TEST(Coordinate_Set)
{

    shared_ptr<PSI_Dummy> psi( new PSI_Dummy);

    double omx=90,omy=80,omz=7;
    double vx=100,vy=90,vz=0;
    double ori = -1; // looking down
    double zenith = 90;
    double azimuth = -170;
    double length = -1; // infinite
    double energy = 1000; // GeV
    int ptype = 1; // cascade

    double foo = -500;
    double bar = +500;

    // initialize to some random values
    PSI_Coordinate* coord1 =
                    dynamic_cast<PSI_Coordinate*>(
                                 psi->MakeCoordinate( foo, foo, foo, 0, bar, bar, bar,
                                                           0, 0, 0, 0, 0 ) );
    // initialize to some useful values
    PSI_Coordinate* coord2 =
                    dynamic_cast<PSI_Coordinate*>(
                                 psi->MakeCoordinate( omx, omy, omz, ori, vx, vy, vz,
                                                           zenith, azimuth, length, energy, ptype ) );

    // coord2->SetFinite(false);

    // reset to some useful values
    coord1->Clear();
    coord1->Set( omx, omy, omz, ori, vx, vy, vz,
                 zenith, azimuth, length, energy, ptype );
    // coord1->SetFinite(false);
    // coord1->CalculateCoordinate();

    // these two should now be identical
    double omx1 = coord1->GetOpticalModuleX();
    double omx2 = coord2->GetOpticalModuleX();
    ENSURE( testDouble( omx1, omx2 ),
            "setting coordinate with Set does not give consistent omX");

    // these two should now be identical
    double trky1 = coord1->GetTrackY();
    double trky2 = coord2->GetTrackY();
    ENSURE( testDouble( trky1, trky2 ),
            "setting coordinate with Set does not give consistent trackY");

    // these two should now be identical
    double theta1 = coord1->GetTrackTheta();
    double theta2 = coord2->GetTrackTheta();
    ENSURE( testDouble( theta1, theta2 ),
            "setting coordinate with Set does not give consistent theta");

    // these two should now be identical
    double dist1 = coord1->CalcDistance();
    double dist2 = coord2->CalcDistance();
    ENSURE( testDouble( dist1, dist2 ),
            "setting coordinate with Set does not give consistent omX");

    // these two should now be identical
    double rho1 = coord1->CalcRho();
    double rho2 = coord2->CalcRho();
    ENSURE( testDouble( rho1, rho2 ),
            "setting coordinate with Set does not give consistent rho");

    // these two should now be identical
    double rl1 = coord1->CalcDistanceTrackToEmission();
    double rl2 = coord2->CalcDistanceTrackToEmission();
    ENSURE( testDouble( rl1, rl2 ),
            "setting coordinate with Set does not give consistent dist to Cherenkov emission point");

    // these two should now be identical
    double rcer1 = coord1->CalcDistanceEmissionToOpticalModule();
    double rcer2 = coord2->CalcDistanceEmissionToOpticalModule();
    ENSURE( testDouble( rcer1, rcer2 ),
            "setting coordinate with Set does not give consistent Cherenkov distance");

    // these two should now be identical
    double tgeo1 = coord1->CalcTGeo();
    double tgeo2 = coord2->CalcTGeo();
    ENSURE( testDouble( tgeo1, tgeo2 ),
            "setting coordinate with Set does not give consistent TGeo");

    delete coord1;
    delete coord2;
}

