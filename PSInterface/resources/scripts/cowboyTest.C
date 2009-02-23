/**
 *@file
 *@brief Amplitude cowboy test example for PSInterface
 *
 *@author Thomas Burgess with Johan Lundberg
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
#include <string>
#include <vector>
using namespace std;
#include "TROOT.h"
#include "TRandom.h"
#include "TMath.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TBenchmark.h"
#include "TString.h"
#include "TH2F.h"

#include "PSInterface/PSInterface.h"
#include "PSInterface/PSI_Coordinate.h"
#include "psimaker.h"
#endif

/**
 *@brief Draw and save cowboy test histogram
 *
 *@param xTitle    Histogram X-axis title
 *@param yTitle    Histogram Y-axis title
 *@param histos    Vector with histograms to plot
 *@param min_x     Histogram X-axis min
 *@param max_x     Histogram X-axis max
 *@param min_y     Histogram Y-axis min
 *@param max_y     Histogram Y-axis max
 *@param logy     True if logscale Y-axis max
 *@param canvas    TCanvas to draw on
 *@param id        SubPad of canvas to draw on
 */
void drawHistograms (    
    string xTitle, string yTitle, TObjArray histos, 
    double min_x, double max_x, double min_y, double max_y, bool logy,
    TCanvas* canvas, int id)
{
    canvas->cd(id);
    if (logy) canvas->GetPad(id)->SetLogy();
    canvas->GetPad(id)->SetGrid();

    for ( int om = 0; om < histos.GetSize(); ++om ) 
    {	
      if ( (dynamic_cast<TH2F*>(histos[om])) == 0 ) break;

	int col = om;
	switch ( om ) {
	case 0: col = 5; break; //Yellow, -200
	case 1: col = 3; break; //Green,  -100
	case 2: col = 1; break; //Black,   0
	case 3: col = 2; break; //Red,     100
	case 4: col = 4; break; //Blue     200
	}
	TProfile* profile = ((TH2F*)histos[om])->ProfileX();
	profile->SetLineColor(col);
	profile->SetAxisRange(min_x,max_x,"X");
	profile->SetAxisRange(min_y,max_y,"Y");
	profile->SetXTitle(xTitle.c_str());
	profile->SetYTitle(yTitle.c_str());
	(om==0)?profile->Draw():profile->Draw("same");
	profile->Clear();
    }
}

/**
 *@brief Produce evenly distributed random angles on a sphere
 *@param theta result zenith angle (radians)
 *@param phi   result azimuth angle (radians)
 */
void SphereRandom( double& theta, double& phi){
    theta = acos(2.0*gRandom->Rndm()-1);  
    phi = 2.0*TMath::Pi()*gRandom->Rndm()-TMath::Pi();
}

/**
 *@brief Make tack randomly distributed within a sphere
 *
 *@param TR_x Resutling track x-coordinate 
 *@param TR_y Resutling track y-coordinate 
 *@param TR_z Resutling track z-coordinate 
 *@param TR_th Resutling track theta angle
 *@param TR_ph Resutling track phi angle
 *@param rCowboy Radius of sphere
 */
void MakeTrack( 
    double &TR_x, 
    double &TR_y, 
    double &TR_z, 
    double &TR_th, 
    double &TR_ph,
    double rCowboy ) 
{
    //Calculate random position on sphere   
    double theta, phi;
    SphereRandom(theta,phi);
    //Calculate random radius	
    const double r = rCowboy*TMath::Power(gRandom->Rndm(),1.0/3.0);
    //Calculate x,y,z from theta,phi,r
    TR_x = sin(theta)*cos(phi)*r;
    TR_y = sin(theta)*sin(phi)*r;
    TR_z = cos(theta)*r;
    //Create a random angle for the track
    SphereRandom(TR_th,TR_ph);
    TR_th*=TMath::RadToDeg();
    TR_ph*=TMath::RadToDeg();
}

/**
 *@brief Sets root style settings
 */
void RootStyle()
{
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
}

/**
 *@brief PSInterface direct hits test
 *
 * Run with ".x cowboyTest.C" in root or  "root -q -l -b cowboyTest.C"
 *  from commandline
 *
 * This script can be compiled using the sample scripts make file 
 * PSInterface/resources/scripts/Makefile
 *
 * Remeber to set LD_LIBRARY_PATH so root can find
 * - libPSInterface.so
 * - liblevel2amasim.so
 * - libphotoamasim.so
 * - libphotonics.so
 * - libptd.so
 * (or copy the libraries(or softlink) to current directory)
 *
 * This program quits as soon as it is done so one might miss any graphical
 * output when compiling the macro.
 *
 *@author Thomas Burgess with Johan Lundberg
 */
int cowboyTest()
{
#ifdef __CINT__
    //Load the psi library and psimaker if using root interactive
    if ( gSystem->Load("libPSInterface.so") < 0 ) {
        return -1;
    }
    gROOT->ProcessLine(".L psimaker.C");
#endif

    //Set root style
    RootStyle();

    TCanvas* canvas = new TCanvas("CowboyTest","Cowboy Test", 800,400);
    canvas->Divide(4,1);

    string interface = "PTD";

    PSInterface* psi = 0;
    if ( ( psi = psimaker(
               interface,
               "level1_table.list",
               "level2_table.list",
	       "tables5_bulk_mam/31ce5.9s21h2.pb.125.3302_011.nta",
               501 ) ) == 0 ) {
	cerr << "Somehow no PSI was created, exiting" << endl;
	return 0;
    }

    TObjArray rhoHistos;
    TObjArray distHistos;
    TObjArray timeHistos;
    TObjArray distTimeHistos;

    const double rCowboy=100.0;
    const double nTracks = 1000; //use 10000 or more
    const int nXbins = 50;         //Number of x bins
    const int nOMs = 5;            //Number of OMs
    const double OMZmax = 200.0 ;  //Maximum OM z coordinate
    const double zoffset =0.0;     //Shift everything in z by this amount    

    const double dmax = 300;        //Maximum distance (and rho)
    const double tmax = 5000;       //Maximum timedelay

    //Create array wity y bin ranges
    const int nYbins = 500;     //Number of bins
    const double amppicmin=1e-3;
    const double amppicmax=1e3;
    const Double_t ymin = amppicmin/1000.0;  //lower y
    const Double_t ymax = amppicmax*1.0;
    const Double_t binwidthy = 
	(TMath::Log10(ymax)-TMath::Log10(ymin))/nYbins;
    double ybins[nYbins+1];    //Array with bin ranges
    for (int i=0; i<=nYbins; ++i) {
	ybins[i] = ymin + 
	    TMath::Power(10,TMath::Log10(ymin)+i*binwidthy);
    }

    TBenchmark bench;
    bench.Start("processing");  

    //Loop over tracks
    for (int track=0; track<nTracks; ++track) {
        if (track%10000==0) cout << "Processing track " << track << endl;
	double TR_x, TR_y, TR_z, TR_th, TR_ph;
 	MakeTrack(TR_x, TR_y, TR_z, TR_th, TR_ph, rCowboy);
	
	//Loop over OMs
	for (int i=0; i< nOMs; ++i){ 
	    double OM_z = -OMZmax+zoffset;  
	    if (nOMs>1) OM_z += i*2*OMZmax/(nOMs-1);
	    
	    //Make a coordinate
	    PSI_Coordinate* coord = 0;	    
	    if ( ( coord = psi->MakeCoordinate(
		       0,0,OM_z,-1,
		       TR_x,TR_y,TR_z,TR_th,TR_ph,
		       0,1,1) ) == 0 )
            {
                cerr << "Failed to make a coordinate" << endl;
                return 0;
            }
	    
	    //Get an amplitude and timedelay
	    double amplitude = 0;
	    if (!psi->GetMeanAmplitude(amplitude,coord)) amplitude = 0;
	    double dt=0;
	    if ((amplitude==0)||
		(!psi->GetTimeDelay(dt,gRandom->Rndm(),coord))) dt = 0;
	    	    
	    if ( track == 0 ) {
		//First track, make histograms
		TString s="OM: "+i;
		rhoHistos.Add(new TH2F(s+" rho","MA vs Rho",
					     nXbins,0,dmax,nYbins,ybins));
		distHistos.Add(new TH2F(s+" Dist","MA vs Dist",
					      nXbins,0,dmax,nYbins,ybins));
		timeHistos.Add(new TH2F(s+" time","MA vs Time",
					      nXbins,0,tmax,nYbins,ybins));
  		distTimeHistos.Add(new TH2F(s+" dist time","Dist vs Time",
					    nXbins,0,dmax,nXbins,0,tmax));
	    }
	    //Fill histograms
	    ((TH2F*)rhoHistos[i])->Fill(coord->CalcRho(),amplitude);
	    ((TH2F*)distHistos[i])->Fill(coord->CalcDistance(),amplitude);
	    ((TH2F*)timeHistos[i])->Fill(dt,amplitude);
	    ((TH2F*)distTimeHistos[i])->Fill(coord->CalcDistance(),dt);

	    delete coord;
	}         
    }
    bench.Stop("processing");
    bench.Show("processing");    

    //Save pictures  if ( 

    drawHistograms ("Rho (m)","MeanAmplitude (NPE)",
 		    rhoHistos,0,dmax,amppicmin/100.0,amppicmax/10.0,
 		    true,canvas,1);
    drawHistograms ("Distance (m)","MeanAmplitude (NPE)",
 		    distHistos,0,dmax,amppicmin/100.0,amppicmax/10.0,
 		    true,canvas,2);
    drawHistograms ("TimeDelay (ns)","MeanAmplitude (NPE)",
		    timeHistos,0,tmax,amppicmin/100.0,amppicmax/10.0,
		    true,canvas,3);
    drawHistograms ("Distance (m)","TimeDelay (ns)",
		    distTimeHistos,0,dmax,0,tmax,
		    false,canvas,4);
    TString file(interface.c_str());
    file+=".gif";
    canvas->SaveAs(file);    
}

#ifndef __CINT__
int main() { return cowboyTest(); }
#endif
