/**
 * \file  test/PSI_PhotonicsPhotorec.cxx 
 * \brief Simple test making sure that multiple psi photonics
 *  objects can coexist.
 *
 *@author Johan Lundberg
 *
 * (c) the IceCube Collaboration
 */

#include <I3Test.h>
#include "PSInterface/PSI_Photonics.h"
#include "PSInterface/PSI_PhotonicsPhotorec.h"
#include<iostream>
#include <stdio.h>

using namespace std;
TEST_GROUP(PSI_Photonics_CPPio);
//Return true if a=b within accuracy

inline bool testDouble( double a, double b )
{
    //Accuracy foe double should be ok upto ~10^15
    const double accuracy=1e7;
    return (round(a*accuracy)==round(b*accuracy));
}

inline bool test( double a, double b )
{
    //Accuracy foe double should be ok upto ~10^15
    const double accuracy=1e7;
    return (round(a*accuracy)==round(b*accuracy));
}


TEST(PSI_Photonics_Multiple)                             
{

  double a=1.21,b=1.31,c,d;
  
  double e,f;
  
  PSI_Photonics psi_1;
  PSI_Photonics psi_2;
  
  psi_1.SetRefRefractive(2,a);
  psi_2.SetRefRefractive(2,b);

  psi_1.GetRefRefractive(2,c,e);
  psi_2.GetRefRefractive(2,d,f);

  ENSURE(((a==c)&&(b==d)),"Set and Get did not match");

}






















