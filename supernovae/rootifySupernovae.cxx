// Rootify supernovae tsvs for both catalogs

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

void rootifySupernovaeOSC()
{

  TFile file("./data/supernovaeOSC.root","RECREATE");
  TTree oscTree("oscTree","Open Supernova Catalog tree");

  // Scalar vars
  std::string name;
  std::string alias;
  std::string discoverer;
  int discoverUnixTime;
  int maxUnixTime;
  float maxappmag;
  std::string host;
  float ra;
  float dec;
  float redshift;
  std::string claimedtype;

  oscTree.Branch("name",&name);
  oscTree.Branch("alias",&alias);
  oscTree.Branch("discoverer",&discoverer);
  oscTree.Branch("discoverUnixTime",&discoverUnixTime,"discoverUnixTime/I");
  oscTree.Branch("maxUnixTime",&maxUnixTime,"maxUnixTime/I");
  oscTree.Branch("maxappmag",&maxappmag,"maxappmag/F");
  oscTree.Branch("host",&host);
  oscTree.Branch("ra",&ra,"ra/F");
  oscTree.Branch("dec",&dec,"dec/F");
  oscTree.Branch("redshift",&redshift,"redshift/F");
  oscTree.Branch("claimedtype",&claimedtype);

  std::ifstream ifile("./data/supernovaeOSC.tsv");
  std::string line;
  std::ifstream vifile;
  std::string vline;
  
  int lineNumber = 0;
  
  while (std::getline(ifile, line))
    {
      lineNumber++;
      std::istringstream iss(line); // construct a string stream from line
      // read the tokens from current line separated by comma
      std::vector<std::string> tokens; // here we store the tokens
      std::string token; // current token
      // much better to use tab as a delimiter over , and ; (lists, field names, citations etc use these!)
      while (std::getline(iss, token, '\t'))
	{
	  trim(token, '*');
	  tokens.push_back(token); // add the token to the vector
	}

      // map the tokens into our variables
      name = tokens[0];
      alias = tokens[1];
      discoverer = tokens[2];
      discoverUnixTime= atoi(tokens[3].c_str());
      maxUnixTime = atoi(tokens[4].c_str());
      maxappmag = atof(tokens[5].c_str());
      host = tokens[6];
      ra = atof(tokens[7].c_str());
      dec = atof(tokens[8].c_str());
      redshift = atof(tokens[9].c_str());
      claimedtype = tokens[10];
	
      oscTree.Fill();

      if(false)
	{
	  std::cout << "" << std::endl;
	  std::cout << "___Example supernova___" << std::endl;
	  std::cout << "Name: " << name << std::endl;
	  std::cout << "alias: " << alias << std::endl;
	  std::cout << "discoverUnixTime: " << discoverUnixTime << std::endl;
	  std::cout << "maxUnixTime: " << maxUnixTime << std::endl;
	  std::cout << "ra: " << ra << std::endl;
	  std::cout << "dec: " << dec << std::endl;
	  std::cout << "type: " << claimedtype << std::endl;
	  std::cout << "redshift: " << redshift << std::endl;
	  std::cout << "_______________________" << std::endl;
	  std::cout << " " << std::endl;
	} 
    }

  oscTree.Write();
  file.Close();
  cout << "Processed OSC!" << endl;
  
  return;
}

void rootifySupernovaeTNS()
{

  TFile file("./data/supernovaeTNS.root","RECREATE");
  TTree tnsTree("tnsTree","Transient Name Server tree");

  // Scalar vars from tsv
  std::string id;
  std::string name;
  std::string raAlt;
  std::string decAlt;
  std::string fullObjType;
  float redshift;
  std::string hostName;
  float hostRedshift;
  std::string discoveryGroup;
  std::string classifyingGroup;
  std::string associatedGroup;
  std::string discoveryInternalName;
  std::string discoveryInstrument;
  std::string classifyingInstrument;
  std::string tnsAT;
  std::string publicSN;
  float discoveryMag;
  std::string discoveryMagFilter;
  std::string discoveryDate;

  // intermediate vars (not from tsv, not kept)
  int raHours; int raMins; float raSecs;
  int decDegs; int decMins; float decSecs;
  
  // New scalar vars
  float ra; // ra in 0 -> 360 deg format
  float dec; // dec in -90 -> 90 deg format
  int discoveryUnixTime; // converted from date
  std::string objSubtype; // I (SN I) or II (SN II) or more complex types (SLSN), etc
  std::string objSubtypeSuffix; // a (SNIa, SNIIb), b-pec (SNIb-pec) etc

  tnsTree.Branch("id",&id);
  tnsTree.Branch("name",&name);
  tnsTree.Branch("ra",&ra);
  tnsTree.Branch("dec",&dec);
  tnsTree.Branch("fullObjType",&fullObjType);
  tnsTree.Branch("objSubtype",&objSubtype);
  tnsTree.Branch("objSubtypeSuffix",&objSubtypeSuffix);
  tnsTree.Branch("redshift",&redshift);
  tnsTree.Branch("hostName",&hostName);
  tnsTree.Branch("hostRedshift",&hostRedshift);
  tnsTree.Branch("discoveryGroup",&discoveryGroup);
  tnsTree.Branch("classifyingGroup",&classifyingGroup);
  tnsTree.Branch("associatedGroup",&associatedGroup);
  tnsTree.Branch("discoveryInternalName",&discoveryInternalName);
  tnsTree.Branch("discoveryInstrument",&discoveryInstrument);
  tnsTree.Branch("classifyingInstrument",&classifyingInstrument);
  tnsTree.Branch("tnsAT",&tnsAT); // 
  tnsTree.Branch("publicSN",&publicSN);
  tnsTree.Branch("discoveryMag",&discoveryMag);
  tnsTree.Branch("discoveryMagFilter",&discoveryMagFilter);
  tnsTree.Branch("discoveryDate",&discoveryDate);
  tnsTree.Branch("discoveryUnixTime",&discoveryUnixTime);
  
  std::ifstream ifile("./data/supernovaeTNS.tsv");
  std::string line;
  std::ifstream vifile;
  std::string vline;
  
  int lineNumber = 0;
  
  while (std::getline(ifile, line))
    {
      lineNumber++;
      // skip header
      if(lineNumber > 1)
	{
	  std::istringstream iss(line); // construct a string stream from line
	  // read the tokens from current line separated by comma
	  std::vector<std::string> tokens; // here we store the tokens
	  std::string token; // current token
	  // much better to use % as a delimiter over , and ; (field names, citations etc use these!)
	  while (std::getline(iss, token, '\t'))
	    {
	      trim(token, '\"');
	      tokens.push_back(token); // add the token to the vector
	    }

	  // map the tokens into our variables
	  id = tokens[0];
	  name = tokens[1];
	  raAlt = tokens[2];
	  decAlt = tokens[3];
	  fullObjType = tokens[4];
	  //
	  redshift = atof(tokens[5].c_str());
	  hostName = tokens[6]; if(hostName == ""){hostName = "Unknown";}
	  hostRedshift =  atof(tokens[7].c_str()); if(hostRedshift == 0){hostRedshift = -999;} // if unknown (i.e. a blank field = 0, set to -999)
	  discoveryGroup = tokens[8]; if(discoveryGroup == ""){discoveryGroup = "Unknown";}
	  classifyingGroup = tokens[9]; if(classifyingGroup == ""){classifyingGroup = "Unknown";}
	  associatedGroup = tokens[10]; if(associatedGroup == ""){associatedGroup = "Unknown";}
	  discoveryInternalName = tokens[11]; if(discoveryInternalName == ""){discoveryInternalName = "Unknown";}
	  discoveryInstrument = tokens[12]; if(discoveryInstrument == ""){discoveryInstrument = "Unknown";}
	  classifyingInstrument = tokens[13]; if(classifyingInstrument == ""){classifyingInstrument = "Unknown";}
	  tnsAT = tokens[14];
	  publicSN = tokens[15];
	  //endProp = tokens[16]; // isn't used
	  discoveryMag = atof(tokens[17].c_str()); if(discoveryMag == 0){discoveryMag = -999;} 
	  discoveryMagFilter = tokens[18]; if(discoveryMagFilter == ""){discoveryMagFilter = "Unknown";}
	  discoveryDate = tokens[19]; if(discoveryDate == ""){discoveryDate = "Unknown";}

	  // intermediate vars
	  // RA
	  int pass = 0;
	  string jub;
	  string raAltSep[3] = {};
	  string supernovaString = "SN ";
	  //size_t prefix;
	  // Extract substrings
	  while(jub != raAlt)
	    {
	      jub = raAlt.substr(0,raAlt.find_first_of(":"));
	      raAlt = raAlt.substr(raAlt.find_first_of(":") + 1);
	      raAltSep[pass] = jub.c_str();
	      //cout << raAltSep[pass] << endl;
	      pass++;
	    }
	  ra = 15. *  (stof(raAltSep[0]) + stof(raAltSep[1])/60. + stof(raAltSep[2])/3600.);
	  // END RA

	  //cout << decAlt << endl;
	  // START DEC
	  int pass2 = 0;
	  string jub2;
	  string decAltSep[3] = {};
	  // Extract substrings
	  while(jub2 != decAlt)
	    {
	      jub2 = decAlt.substr(0,decAlt.find_first_of(":"));
	      decAlt = decAlt.substr(decAlt.find_first_of(":") + 1);
	      decAltSep[pass2] = jub2.c_str();
	      //cout << decAltSep[pass2] << endl;
	      pass2++;
	    }	  
	  // If negative attached to hour...
	  if (decAltSep[0].find('-') != std::string::npos)
	    {
	      //cout << "deg = " << stof(decAltSep[0]) << ", m = " << stof(decAltSep[1]) << ", s = " << stof(decAltSep[2]) << endl;
	      dec = -1.*(-stof(decAltSep[0]) + stof(decAltSep[1])/60. + stof(decAltSep[2])/3600.);
	    }
	  else
	    {
	      dec = stof(decAltSep[0]) + stof(decAltSep[1])/60. + stof(decAltSep[2])/3600.;
	    }				      
	  // END DEC

	  // START TIME
	  // YYYY-MM-DD HH:MM:SS -> unixTime
	  //cout << discoveryDate << endl;
	  struct tm tm;
	  time_t ts = 0;
	  memset(&tm, 0, sizeof(tm));
	  std::string humanTime = discoveryDate;
	  strptime(humanTime.c_str(), "%Y-%m-%d %T", &tm);
	  ts = mktime(&tm);
	  discoveryUnixTime = (int)ts;
	  //cout << discoveryUnixTime << endl;
	  // END TIME

	  // Supernovae have complex names sometimes...
	  // Also SNI is contained within the SNII string...
	  if(fullObjType.find("SN II") == 0)
	    {
	      objSubtype = "II";
	    }
	  else if(fullObjType.find("SN I") == 0)
	    {
	      objSubtype = "I";
	    }
	  else
	    {
	      objSubtype = fullObjType; // Non-standard SN, such as SLSN
	    }

	  // Keep suffix
	  // SN 1a-jub, suffix is a-jub
	  if(objSubtype!=fullObjType)
	    {
	      objSubtypeSuffix = fullObjType.substr(fullObjType.find(objSubtype)+1);
	    }
	  else
	    {
	      objSubtypeSuffix = "None";
	    }
	  /*
	  std::cout << "full = " << fullObjType << std::endl;
	  std::cout << "sub = " << objSubtype << std::endl;
	  std::cout << "suffix = " << objSubtypeSuffix << std::endl;
	  */
	  
	  tnsTree.Fill();
	}
    }
  tnsTree.Write();
  file.Close();
  cout << "Processed TNS!" << endl;
  
  return;
}


void rootifySupernovae()
{
  bool processOSC = true;
  bool processTNS = true;
  
  if(processOSC == true){cout << "Rootifying OSC..." << endl; rootifySupernovaeOSC(); cout << "" << endl;}
  if(processTNS == true){cout << "Rootifying TNS..." << endl; rootifySupernovaeTNS(); cout << "" << endl;}

  cout << "All done!" << endl;
}
