// Rootify GRB tsvs for SWIFT catalog for all time

#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <string>
#include "TTree.h"
#include "TFile.h"

// This gets rid of unwanted chars: *, ? etc from data
void trim(std::string& s, const char t)
{
  s.erase(std::remove(s.begin(), s.end(), t), s.end());
}

void rootifyGRB()
{

  TFile file("./data/GRB.root","RECREATE");
  TTree grbTree("grbTree","Gamma Ray Burst SWIFT tree");

  // Scalar vars
  std::string name;
  std::string datePre;
  std::string letterIndex;
  std::string time;
  float RA; // BAT
  float dec; // BAT
  float T90; // BAT
  float fluence; // BAT
  float peakPhotonFluxOneSec; // BAT
  float photonIndex;
  std::string powerLawType;
  float gamma; // XRT

  // Immediate vars
  std::string tempIndex;
  std::string tempIndexFound;
  std::string na = "n/a";

  grbTree.Branch("name",&name);
  grbTree.Branch("datePre",&datePre);
  grbTree.Branch("letterIndex",&letterIndex);
  grbTree.Branch("time",&time);
  grbTree.Branch("RA",&RA,"RA/F");
  grbTree.Branch("dec",&dec,"dec/F");
  grbTree.Branch("T90",&T90,"T90/F");
  grbTree.Branch("fluence",&fluence,"fluence/F");
  grbTree.Branch("peakPhotonFluxOneSec",&peakPhotonFluxOneSec,"peakPhotonFluxOneSec/F");
  grbTree.Branch("photonIndex",&photonIndex,"photonIndex/F");
  grbTree.Branch("powerLawType",&powerLawType);
  grbTree.Branch("gamma",&gamma,"gamma/F");
  
  std::ifstream ifile("./data/swift_grb.tsv");
  std::string line;
  std::ifstream vifile;
  std::string vline;
  
  int lineNumber = 0;
  
  while (std::getline(ifile, line))
    {
      lineNumber++;
      if(lineNumber > 1)
	{
	  std::istringstream iss(line); // construct a string stream from line
	  // read the tokens from current line separated by comma
	  std::vector<std::string> tokens; // here we store the tokens
	  std::string token; // current token
	  // much better to use tab as a delimiter over , and ; (lists, field names, citations etc use these!)
	  while (std::getline(iss, token, '\t'))
	    {
	      //trim(token, '*');
	      tokens.push_back(token); // add the token to the vector
	    }

	  // map the tokens into our variables
	  name = tokens[0];
	  datePre = name.substr(0,6);
	  letterIndex = name.substr(6,7);
	  time = tokens[1];
	  // Some fields (except name time) have n/a for a full component, or only n/a for some vars...
	  // i.e. only a single part of SWIFT detected something, such as the XRT
	  // Assign -999 for non-string vars and "n/a" for strings which are n/a.
	  RA = (tokens[3] == na) ? -999 : atof(tokens[3].c_str());
	  dec = (tokens[4] == na) ? -999 : atof(tokens[4].c_str());
	  T90 = (tokens[6] == na) ? -999 : atof(tokens[6].c_str());
	  fluence = (tokens[7] == na) ? -999 : atof(tokens[7].c_str());
	  peakPhotonFluxOneSec = (tokens[9] == na) ? -999 : atof(tokens[9].c_str());
	  tempIndex = tokens[11]; // get power law string
	  tempIndexFound = tempIndex.substr(0, tempIndex.find(","));
	  gamma = (tokens[21] == na) ? -999 : atof(tokens[21].c_str());	  	  
	  
	  if(tempIndexFound == na)
	    {
	      photonIndex = -999;
	      powerLawType = "na";
	    }
	  else
	    {
	      photonIndex = atof(tempIndexFound.c_str());
	      powerLawType = tempIndex.substr(tempIndex.find(",") + 2); 
	    }
	  /*  
	  if(lineNumber < 6)
	    {
	      cout << "name = " << name << endl;
	      cout << "time = " << time << endl;
	      cout << "ra = " << RA << endl;
	      cout << "photonindex = " << photonIndex << endl;
	      cout << "power law = " << powerLawType << endl;
	      cout << "gamma " << gamma << endl;
	      cout << "---" << endl;
	    }
	  */	  		
	  grbTree.Fill();
	    
	}
    }
  
  grbTree.Write();
  file.Close();
  cout << "Processed SWIFT GRB catalog!" << endl;
  
  return;
}
