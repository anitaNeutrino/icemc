/*
  Macro to give useful information about flares during an ANITA flight (or over any timescale): timing, which source, which type, etc.
  Can print unique, catalogued sources if specified, and plot them on a sky map.
*/
#include <iostream>
#include <libgen.h>
#include <map>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TH1F.h"
#include "TMath.h" 
#include "TObjString.h" 
#include "TTimeStamp.h"
#include "fava.h" 
#include "TH1.h"
#include "TGraph.h"
#include "TMarker.h"
#include "SkyMap.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TCanvas.h"
#include "TStyle.h"

void getFlareInfo()
{
  ///////////////////////////////SET THESE TO WHAT YOU DESIRE////////////////
  bool ignoreAssoc = true; // set to true to only look for known (associated) sources, IGNORE those that don't point to a catalogued astrophysical object
  bool uniqueSourcesOnly = true; // set to get a list of unique astrophyical objects associated with these blazar flares
  int anitaFlight = 4; // set to which flight you want. i.e. 3 for ANITA-3. Special option: 0 to use all data from all times (including those when anita wasn't flying) from FAVA
  bool skyMap = true; // print a sky map of all the unique catalogued objects flaring during the flight
  /////////////////////////////////////////////////////////////////////////////////////

  
  TFile *file = new TFile("fava.root"); 
  TTree *tree = (TTree*) file->Get("fava");
  
  FAVAEntry *fava = new FAVAEntry();
 
  tree->SetBranchAddress("fava",&fava);
  int a_tmin;
  int a_tmax;
  if(anitaFlight == 4)
    {
      a_tmin = 1480707643;
      a_tmax = 1483004563;
    }
  else if(anitaFlight == 3)
    {
      a_tmin = 1418938406; 
      a_tmax = 1420777814;
    }
  else if(anitaFlight == 0)
    {
      a_tmin = 1217864618;
      a_tmax = 1546271018;
    }
  else
    {
      cout << "You specified a flight which isn't supported! Using A4 flight times" << endl;
      a_tmin = 1480707643;
      a_tmax = 1483004563;
    }
 
  int currentWeek = 0;
  bool foundFlaresPreviously = false; // always should start as false
  int flareTotal = 0;
  int uniqueSources = 0;
  std::vector<TString> vFlareList;
  std::vector<double> vRA;
  std::vector<double> vDEC;
  std::vector<TString> vSourceClass;

  for (int i = 0; i < tree->GetEntries(); i++) 
    {
      tree->GetEntry(i);

      // If week times fall completely out of ANITA flight times, ignore
      if( (fava->unix_tmin < a_tmin && fava->unix_tmax < a_tmin) || (fava->unix_tmin > a_tmax && fava->unix_tmax > a_tmax))
	{
	  continue;
	}
      else
	{
	  if(foundFlaresPreviously == false)
	    {
	      cout << "Listing found flares:" << endl;
	    }
	  
	  if(fava->week > currentWeek)
	    {
	      std::cout << "___Week" << fava->week << "___" << std::endl;
	    }

	  currentWeek = fava->week;
	  
	  if(ignoreAssoc==true) //If true, only consider those which have an association
	    {
	      if(fava->association.GetString() != "None")
		{
		  std::cout << "Num:" << fava->num;
		  cout << " \t Source class:" << fava->source_class.GetString();
		  std::cout <<  "\t RA:" << fava->ra << "\t Dec:" << fava->dec;
		  cout << "\t Assoc:" << fava->association.GetString() << endl;
		  vFlareList.push_back(fava->association.GetString());
		  vSourceClass.push_back(fava->source_class.GetString());
		  vRA.push_back(fava->ra);
		  vDEC.push_back(fava->dec);
		  foundFlaresPreviously = true;
		  flareTotal++;
		}
	      
	    }

	  else // Consider all
	    {
	      std::cout << "Num:" << fava->num;
	      cout << " \t Source class:" << fava->source_class.GetString();
	      std::cout <<  "\t RA:" << fava->ra << "\t Dec:" << fava->dec;
	      cout << "\t Assoc:" << fava->association.GetString() << endl;
	      vFlareList.push_back(fava->association.GetString());
	      vSourceClass.push_back(fava->source_class.GetString());
	      vRA.push_back(fava->ra);
	      vDEC.push_back(fava->dec);
	      foundFlaresPreviously = true;
	      flareTotal++;
	    }	  
	  
	}

    }
  
  // Drawing vars
  TMarker *astroObject = new TMarker();
  SkyMap *skyMapOut = new SkyMap();
  TLegend *legend = new TLegend(0.8,0.8,0.99,0.99);
  TLegendEntry *le = new TLegendEntry();
  
  // Which source classes and how many of them do we see?
  // Include all object types for completeness
  bool bllPass = false; int bllCount = 0; // BL-Lac
  bool fsrqPass = false; int fsrqCount = 0; // flat spectrum radio quasar
  bool bcuPass = false; int bcuCount = 0; // blazar candidate (uncertain type)
  bool rdgPass = false; int rdgCount = 0; // radio galaxy
  bool hmbPass = false; int hmbCount = 0; // high mass binary
  bool psrPass = false; int psrCount = 0; //pulsar
  bool nlsy1Pass = false; int nlsy1Count = 0; // narrow line Seyfert 1
  bool pwnPass = false; int pwnCount = 0; //pulsar wind nebula
  bool agnPass = false; int agnCount = 0; //active galactic nuclei
  bool novPass = false; int novCount = 0; //nova
  
  // Print out unique sources
  if(uniqueSourcesOnly == true)
    {
      cout << "___Catalogued sources list___" << endl;
      for(unsigned int i = 0; i < vFlareList.size(); i++)
	{
	  unsigned int j;
	  for(j=0; j<i; j++)
	    {
	      if(vFlareList[i] == "None") // Don't print "none", and don't assign it as a catalogued source
		{
		  break;
		}
	      if(vFlareList[i] == vFlareList[j]) // Only unique
		{
		  break;
		}	  
	    }
	  if(i==j)
	    {
	      uniqueSources++;
	      cout << "   * " << vFlareList[i] << endl;

	      // We only want to draw each catalogued source on the map once
	      if(skyMap==true)
		{
		  astroObject->SetX(vRA[i]);
		  astroObject->SetY(vDEC[i]);
		  skyMapOut->addMarker(astroObject);
		  // Yes, this is inefficient, but TLegend doesn't work well otherwise
		  if(vSourceClass[i] == "bll")
		    {
		      astroObject->SetMarkerColor(1);
		      if(bllPass == false)
			{
			  le = legend->AddEntry(astroObject,"BLL","*");
			  le->SetTextColor(1);
			  bllPass = true;
			}
		      bllCount++;
		    }
		  else if(vSourceClass[i] == "fsrq")
		    {
		      astroObject->SetMarkerColor(2);
		      if(fsrqPass == false)
			{
			  le = legend->AddEntry(astroObject,"FSRQ","*");
			  le->SetTextColor(2);
			  fsrqPass = true;
			}
		      fsrqCount++;
		    }
		  
		  else if(vSourceClass[i] == "bcu")
		    {
		      astroObject->SetMarkerColor(3);
		      if(bcuPass == false)
			{
			  le = legend->AddEntry(astroObject,"BCU","*");
			  le->SetTextColor(3);
			  bcuPass = true;
			}
		      bcuCount++;
		    }
		  else if(vSourceClass[i] == "rdg")
		    {
		      astroObject->SetMarkerColor(4);
		      if(rdgPass == false)
			{
			  le = legend->AddEntry(astroObject,"RDG","*");
			  le->SetTextColor(4);
			  rdgPass = true;
			}
		      rdgCount++;
		    }
		  else if(vSourceClass[i] == "hmb")
		    {
		      astroObject->SetMarkerColor(15);
		      if(hmbPass == false)
			{
			  le = legend->AddEntry(astroObject,"HMB","*");
			  le->SetTextColor(15);
			  hmbPass = true;
			}
		      hmbCount++;
		    }
		  else if(vSourceClass[i] == "psr")
		    {
		      astroObject->SetMarkerColor(6);
		      if(psrPass == false)
			{
			  le = legend->AddEntry(astroObject,"PSR","*");
			  le->SetTextColor(6);
			  psrPass = true;
			}
		      psrCount++;
		    }
		  else if(vSourceClass[i] == "nlsy1")
		    {
		      astroObject->SetMarkerColor(7);
		      if(nlsy1Pass == false)
			{
			  le = legend->AddEntry(astroObject,"NLSY1","*");
			  le->SetTextColor(7);
			  nlsy1Pass = true;
			}
		      nlsy1Count++;
		    }
		  else if(vSourceClass[i] == "pwn")
		    {
		      astroObject->SetMarkerColor(8);
		      if(pwnPass == false)
			{
			  le = legend->AddEntry(astroObject,"PWN","*");
			  le->SetTextColor(8);
			  pwnPass = true;
			}
		      pwnCount++;
		    }
		  else if(vSourceClass[i] == "agn")
		    {
		      cout << "I AM HERE" << endl;
		      astroObject->SetMarkerColor(9);
		      if(agnPass == false)
			{
			  le = legend->AddEntry(astroObject,"AGN","*");
			  le->SetTextColor(9);
			  agnPass = true;
			}
		      agnCount++;
		    }
		  else if(vSourceClass[i] == "nov")
		    {
		      astroObject->SetMarkerColor(kOrange);
		      if(novPass == false)
			{
			  le = legend->AddEntry(astroObject,"NOV","*");
			  le->SetTextColor(kOrange);
			  novPass = true;
			}
		      novCount++;
		    }
		  else // not associated aren't listed here (or pulsar with no pulsations seen by LAT)
		    {
		      cout << "." << endl;
		    }
		  //Global object marker options: not source specific
		  astroObject->SetMarkerStyle(20);
		  astroObject->SetMarkerSize(1.75);
		}
	      
	    }      
	}
    }
  
  cout << "____________Summary_____________" << endl; 
  //cout << flareTotal << " flares occurred during the ANITA-4 flight:" << endl;
  if(uniqueSourcesOnly == true)
    {
      //cout << "   * From " <<  uniqueSources << " catalogued sources" << endl;
      cout << uniqueSources << " catalogued sources were flaring during the flight: " << endl;
      if(bllPass == true){cout << "   * " << bllCount << " BLLs" << endl;}
      if(fsrqPass == true){cout << "   * " << fsrqCount << " FSRQs" << endl;}
      if(bcuPass == true){cout << "   * " << bcuCount << " BCUs" << endl;}
      if(rdgPass == true){cout << "   * " << rdgCount << " RDGs" << endl;}
      if(hmbPass == true){cout << "   * " << hmbCount << " HMBs" << endl;}
      if(psrPass == true){cout << "   * " << psrCount << " PSRs" << endl;}
      if(nlsy1Pass == true){cout << "   * " << nlsy1Count << " NLSY1s" << endl;}
      if(pwnPass == true){cout << "   * " << pwnCount << " PWNs" << endl;}
      if(agnPass == true){cout << "   * " << agnCount << " AGNs" << endl;}
      if(novPass == true){cout << "   * " << novCount << " NOVs" << endl;}
    }

  if(skyMap==true)
    {
      TCanvas *c = new TCanvas();
      c->SetTitle("ANITA flaring sources map ");
      skyMapOut->Draw();
      legend->Draw();
    }
  return;
  
}
