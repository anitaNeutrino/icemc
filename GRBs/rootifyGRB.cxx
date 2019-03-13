// Rootify GRB tsvs for SWIFT and IceCube catalogs for all time

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

bool RemoveExtraneousSpaces(char lhs, char rhs)
{
  return (lhs == rhs) && (lhs == ' ');
}

void rootifyGRBSwift()
{

  TFile file("./data/GRBSwift.root","RECREATE");
  TTree swiftTree("swiftTree","Gamma Ray Burst SWIFT tree");

  // Scalar vars
  std::string name;
  std::string date;
  std::string letterIndex;
  std::string time;
  int unixTime;
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

  swiftTree.Branch("name",&name);
  swiftTree.Branch("date",&date);
  swiftTree.Branch("letterIndex",&letterIndex);
  swiftTree.Branch("time",&time);
  swiftTree.Branch("unixTime",&unixTime, "unixTime/I");
  swiftTree.Branch("RA",&RA,"RA/F");
  swiftTree.Branch("dec",&dec,"dec/F");
  swiftTree.Branch("T90",&T90,"T90/F");
  swiftTree.Branch("fluence",&fluence,"fluence/F");
  swiftTree.Branch("peakPhotonFluxOneSec",&peakPhotonFluxOneSec,"peakPhotonFluxOneSec/F");
  swiftTree.Branch("photonIndex",&photonIndex,"photonIndex/F");
  swiftTree.Branch("powerLawType",&powerLawType);
  swiftTree.Branch("gamma",&gamma,"gamma/F");
  
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
	  // Some fields (except names) have n/a for a full component, or only n/a for some vars...
	  // i.e. only a single part of SWIFT detected something, such as the XRT
	  // Assign -999 for non-string vars and "n/a" for strings which are n/a.
	  date = name.substr(0,6);
	  letterIndex = name.substr(6,7); // skip milliseconds, only HH:MM:SS needed
	  time = (tokens[1] == na) ? "12:00:00" : tokens[1].substr(0,8); // if n/a, just give seconds in the middle of the day

	  struct tm tm;
	  time_t ts = 0;
	  memset(&tm, 0, sizeof(tm));
	  std::string humanTime = date + " " + time;
	  strptime(humanTime.c_str(), "%y%m%d %T", &tm);
	  ts = mktime(&tm);
	  unixTime = (int)ts;

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
	      cout << "date = " << date << endl;
	      cout << "time = " << time << endl;
	      cout << "unixTime = " << unixTime << endl;
	      cout << "ra = " << RA << endl;
	      cout << "photonindex = " << photonIndex << endl;
	      cout << "power law = " << powerLawType << endl;
	      cout << "gamma " << gamma << endl;
	      cout << "---" << endl;
	    }
	  */	  
	  swiftTree.Fill();
	    
	}
    }
  
  swiftTree.Write();
  file.Close();
  cout << "Processed SWIFT GRB catalog!" << endl;
  
  return;
}

void rootifyGRBIcecube()
{

  TFile file("./data/GRBIceCube.root","RECREATE");
  TTree iceCubeTree("iceCubeTree","Gamma Ray Burst IceCube tree");

  // Scalar vars
  std::string iceCubeName;
  std::string fermiName;
  std::string name;
  std::string date;
  int unixTriggerTime;
  float RA;
  float dec;
  float T90; // *duration* of T90 window
  int unixT90Time;
  float fluence;
  float redshift;
  float T100;

  // Intermediate vars
  std::string triggerTime;
  std::string T90Time;

  iceCubeTree.Branch("name",&name);
  iceCubeTree.Branch("date",&date);
  iceCubeTree.Branch("unixTriggerTime",&unixTriggerTime,"unixTriggerTime/I");
  iceCubeTree.Branch("RA",&RA,"RA/F");
  iceCubeTree.Branch("dec",&dec,"dec/F");
  iceCubeTree.Branch("T90",&T90,"T90/F");
  iceCubeTree.Branch("unixT90Time",&unixT90Time,"unixT90Time/I");
  iceCubeTree.Branch("fluence",&fluence,"fluence/F");
  iceCubeTree.Branch("redshift",&redshift,"redshift/F");
  iceCubeTree.Branch("T100",&T100,"T100/F");

  ////////////////////////////////////////////////

  std::ifstream ifile("./data/Summary_table.txt");
  std::string line;
  std::ifstream vifile;
  std::string vline;
  
  int lineNumber = 0;
  
  while (std::getline(ifile, line))
    {
      lineNumber++;
      if(lineNumber > 4 && lineNumber) // skip format, header, units
	{
	  std::istringstream iss(line); // construct a string stream from line
	  // read the tokens from current line separated by comma
	  std::vector<std::string> tokens; // here we store the tokens
	  std::string token; // current token

	  // note: used sed -i 's/ \{1,\}/,/g' file      on original file to replace any amount of space with a single ,
	  //sed -i 's/^[ \t]*//' file      on original file to remove leading whitespace
	  // sed -i 's/^,//' file        on original file to remove first comma if first character is comma
	  while (std::getline(iss, token, ','))
	    {
	      tokens.push_back(token); // add the token to the vector
	    }

	  // map the tokens into our variables
	  iceCubeName = tokens[0];
	  fermiName = tokens[1];

	  if(fermiName!="None" && iceCubeName=="None")
	    {
	      date = fermiName.substr(3,6);
	      name = fermiName;
	    }
	  else if(fermiName=="None" && iceCubeName!="None")
	    {
	      date = iceCubeName.substr(3,6);
	      name = iceCubeName;
	    }
	  else if(fermiName!="None" && iceCubeName!="None")
	    {
	      date = iceCubeName.substr(3,6); // use iceCubeName if both names exist
	      name = iceCubeName;
	    }
	  else
	    {
	      date = "noDate";
	      name= "noName";
	    }
	  
	  triggerTime = tokens[2];
	  

	  // trigger time
	  struct tm tm;
	  time_t ts = 0;
	  memset(&tm, 0, sizeof(tm));
	  std::string humanTime = date + " " + triggerTime;
	  strptime(humanTime.c_str(), "%y%m%d %T", &tm);
	  ts = mktime(&tm);
	  unixTriggerTime = (int)ts;
	  //
	  RA = atof(tokens[3].c_str());
	  dec = atof(tokens[4].c_str());
	  T90 = atof(tokens[6].c_str());
	  T90Time = tokens[8];
	  // same for T90 time
	  struct tm tm2;
	  time_t ts2 = 0;
	  memset(&tm2, 0, sizeof(tm2));
	  std::string humanTime2 = date + " " + T90Time;
	  strptime(humanTime2.c_str(), "%y%m%d %T", &tm2);
	  ts2 = mktime(&tm2);
	  unixT90Time = (int)ts2;
	  // if no T90 Start is NOT recorded (-999), also set the calculated unix time to -999...
	  if(T90Time == "-999"){unixT90Time = -999;}
	  //

	  fluence = atof(tokens[9].c_str());
	  redshift = atof(tokens[11].c_str());
	  T100 = atof(tokens[12].c_str());
	  
	  iceCubeTree.Fill();

	  // account for -999s 
	}
    }
  
  iceCubeTree.Write();
  file.Close();
  cout << "Processed IceCube GRB catalog!" << endl;
  
  return;
  
}

void rootifyGRB()
{
  
  rootifyGRBSwift();
  rootifyGRBIcecube();
  
  return;
  
}
