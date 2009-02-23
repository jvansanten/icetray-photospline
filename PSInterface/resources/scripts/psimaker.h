/**
 *@file
 *@brief Initalizes a new PSInterface instance
 *
 *@author Thomas Burgess
 *
 * (c) the IceCube Collaboration
 *
 * $Revision:$
 * $Author:$
 * $Date:$
 * $Id:$
 */

#ifndef __psimaker_h__
#define __psimaker_h_

#ifndef __CINT__
#include <iostream>
using namespace std;
#include "PSInterface/PSInterface.h"
#endif

/**
 *@brief construct a new PSInterface statement
 *
 *@param interface What interface to construct, valid options are
 * - Dummy
 * - PTD
 * - Photonics 
 *@param photonicsLevel1TableList Name of photonics level 1 list to load
 * only used if interface==Photonics, if ="" no level 1 tables will be loaded
 *@param photonicsLevel2TableList Name of photonics level 2 list to load
 * (only used if interface==Photonics, if ="" no level 2 tables will be loaded
 *@param PTDTable                 Name of PTD table to load 
 * (only used if interface==PTD)
 *@param PTDMode PTD mode
 *
 *@author Thomas Burgess
 *
 *@return PSInterface on success, 0 otherwise
 */
PSInterface* psimaker(
    string interface, 
    string photonicsLevel1TableList, 
    string photonicsLevel2TableList,
    string PTDTable,
    int PTDMode );

#endif
