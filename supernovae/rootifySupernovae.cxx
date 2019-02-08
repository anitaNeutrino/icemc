// Rootify supernovae csvs

#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <string>

// This gets rid of the quote marks
void trim(std::string& s, const char t)
{
  s.erase(std::remove(s.begin(), s.end(), t), s.end());
}

void rootifySupernovae()
{

  TFile file("supernovae.root","RECREATE");
  TTree snTree("snTree","Open Supernova Catalog tree");

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

  snTree.Branch("name",&name);
  snTree.Branch("alias",&alias);
  snTree.Branch("discoverer",&discoverer);
  snTree.Branch("discoverUnixTime",&discoverUnixTime,"discoverUnixTime/I");
  snTree.Branch("maxUnixTime",&maxUnixTime,"maxUnixTime/I");
  snTree.Branch("maxappmag",&maxappmag,"maxappmag/F");
  snTree.Branch("host",&host);
  snTree.Branch("ra",&ra,"ra/F");
  snTree.Branch("dec",&dec,"dec/F");
  snTree.Branch("redshift",&redshift,"redshift/F");
  snTree.Branch("claimedtype",&claimedtype);

  std::ifstream ifile("supernovae.csv");
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
      // much better to use % as a delimiter over , and ; (field names, citations etc use these!)
      while (std::getline(iss, token, '%'))
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
	
      snTree.Fill();

      if(true)
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

  snTree.Write();
  file.Close();
  cout << "Done!" << endl;
  
  return;
}
// Find a way to handle those without info / bugs
///i.e. go to supernovae.root and draw / scan / show variables: do they make sense?
