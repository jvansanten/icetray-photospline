/**
 *@file
 *@brief Direct hit drawing example for PSInterface
 *
 *@author Thomas Burgess
 *
 * (c) the IceCube Collaboration
 *
 * $Revision: $
 * $Author: $
 * $Date: $
 * $Id: $
 */

#ifndef __CINT__
#include <iostream>
using namespace std;
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TBenchmark.h"
#include "TH2F.h"
#include "PSInterface/PSInterface.h"
#include "PSInterface/PSI_Coordinate.h"
#include "psimaker.h"
#endif

/**
 *@brief PSInterface direct hits test
 *
 * Run with ".x directhits.C" in root or  "root -q -l -b directhits.C"
 *  from commandline
 * Remeber to set LD_LIBRARY_PATH so root can find
 * - libPSInterface.so
 * - liblevel2amasim.so
 * - libphotoamasim.so
 * - libphotonics.so
 * - libptd.so
 * (or copy the libraries(or softlink) to current directory)
 *
 * This macro can also be compiled if needed
 *<pre>
 * g++ directhits.C -I../../public/ \
 * `root-config --libs --cflags` -L../.. -lPSInterface \
 * -L${PTDLIB} -I${PTDINC} -lptd \ 
 * -I{PHOTONICSINC}  -L${PHOTONICSLIB} \ 
 * -llevel2amasim -lphotoamasim -lphotonics	\
 * -DENABLE_PSI_LOGGING -o directhits
 *</pre>
 *(PTDINC, PTDLIB, PHOTONICSINC and PHOTONICSLIB should be substituted for
 * library and include directories)
 *
 * This program quits as soon as it is done so one might miss any graphical
 * output when compiling the macro.
 * 
 *@author Thomas Burgess
 */
int directhits()
{
#ifdef __CINT__
    //Load the psi library and psimaker if using root interactive
    if ( gSystem->Load("libPSInterface.so") < 0 ) {
	return -1;
    }
    gROOT->ProcessLine(".L psimaker.C");
#endif   

    //Root style settings
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetCanvasBorderSize(0);
    gStyle->SetTitleYOffset(2);
    gStyle->SetPadRightMargin(0.14);
    gStyle->SetPadLeftMargin(0.16);
    gStyle->SetMarkerColor(1);
    gStyle->SetMarkerSize(0.2);
    gStyle->SetLineWidth(2);

    TCanvas* theCanvas=new TCanvas("the canvas","the canvas",650,600);
    theCanvas->cd();	


    gStyle->SetPalette(128);

    PSInterface* psi = 0;
    if ( ( psi = psimaker(
	       "Photonics", 
	       "level1_table.list",
	       "level2_table.list",
	       "tables5_bulk_mam/31ce5.9s21h2.pb.125.3302_011.nta",
	       500 ) ) == 0 ) {
	cerr << "Somehow no PSI was created, exiting" << endl;
	return 0;
    }

    double OM_x, OM_y, OM_z=0, OM_o=-1;	    
    double TR_x=0, TR_y=0, TR_z=-10, TR_th=90, TR_ph=0,
	TR_length=10, TR_energy=10;
    int TR_type = 1;

    const int nbins=149;          //Should be odd
    const double extent=50.0;    //width of plot

    TBenchmark bench;
    
    TH2F* directhits = new TH2F(
	"directhits","Direct Hits",
	nbins,-extent,extent,nbins,-extent,extent);
    bench.Start("MakeGraph");
    for (int i =0; i<nbins; i++) 
    {	
	//Downlooking optical module at (0,0,0)
	OM_x=((i-nbins/2+1)/((double)nbins-1))*extent*2;	
	for (int j =0; j<nbins; j++) {
	    OM_y=((j-nbins/2+1)/((double)nbins-1))*extent*2;
	    PSI_Coordinate* coord = 0;
	    
	    //Create a coordinate
	    if ( ( coord = psi->MakeCoordinate(
		OM_x,OM_y,OM_z,OM_o,
		TR_x,TR_y,TR_z,TR_th,TR_ph,
		TR_length,TR_energy,TR_type) ) == 0 ) 
	    {
		cerr << "Failed to make a coordinate" << endl;
		return -1;
	    }

	    double amplitude = 0;    //Hit amplitude
	    double time = 0;         //Time delay
	    double prob=0;           //Probabilty
	    double thresh=25;        //Max time
	    const double inc=0.1;    //step size

	    if (!psi->GetMeanAmplitude(amplitude,coord)) {
	   	cerr << "Failed to get mean ampitude" << endl; 
		return -1;
	    } 
 
	    while ( time < thresh ) {
		if (amplitude<=0) break; //Skip loop if no amplitude

		prob+=inc;	//Increase prob until >=1
		if ( prob>=1 ) break; 

		//Get time delay
 		if (!psi->GetTimeDelay(time,prob,coord))  {
		    cerr << "Failed to get probability" << endl;
		    return -1;
		}
		
		if (time>thresh) break; //Stop if threshold exceded
	    }
	    
	    amplitude*=prob; //Weight amplitude by probability
	    //Loop amplitude times and fill histogram (max 10 times)
	    for ( int hits=10; hits>0; --hits) 
		if (amplitude>hits) directhits->Fill(OM_x,OM_y,1); 
	    if (amplitude>0.5) 
		directhits->Fill(OM_x,OM_y,0.5); //Low probs
	    
	    delete coord;
	}
    }
    directhits->Draw("colz");

    bench.Stop("MakeGraph");
    bench.Show("MakeGraph");
    
    theCanvas->SaveAs("directhits.gif");
    return 0;
}

#ifndef __CINT__
int main() { return directhits(); }
#endif
