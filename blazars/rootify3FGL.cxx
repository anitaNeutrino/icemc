// Rootify fava and light curve csvs

#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <string>
#include <vector>

void trim(std::string& s, const char t)
{
  //s.erase(0, s.find_first_not_of(t));
  //s.erase(s.find_last_not_of(t) + 1);
  s.erase(std::remove(s.begin(), s.end(), t), s.end());
}

void rootify3FGL()
{

  TFile f("3FGL.root","RECREATE");
  TTree catTree("catTree","3FGL catalog");

  // CATALOG VARS:

  // Scalar vars
  std::string sourceName;
  std::string association;
  double ra;
  double dec;
  std::string classType;
  std::string tevCat;
  std::string tevAssoc;

  catTree.Branch("sourceName",&sourceName);
  catTree.Branch("association",&association);
  catTree.Branch("ra",&ra,"ra/D");
  catTree.Branch("dec",&dec,"dec/D");
  catTree.Branch("classType",&classType);
  catTree.Branch("tevCat",&tevCat);
  catTree.Branch("tevAssoc",&tevAssoc);

  // Vector vars
  std::vector<int> met; int metTemp;
  std::vector<int> unixTime; int unixTimeTemp;
  std::vector<double> le_count; double le_countTemp;
  std::vector<double> le_bg; double le_bgTemp;
  std::vector<double> le_sig; double le_sigTemp;
  std::vector<double> he_count; double he_countTemp;
  std::vector<double> he_bg; double he_bgTemp;
  std::vector<double> he_sig; double he_sigTemp;

  bool doLightCurves = false;
  if(doLightCurves == true)
    {
      catTree.Branch("met",&met);
      catTree.Branch("unixTime",&unixTime);
      catTree.Branch("le_count",&le_count);
      catTree.Branch("le_bg",&le_bg);
      catTree.Branch("le_sig",&le_sig);
      catTree.Branch("he_count",&he_count);
      catTree.Branch("he_bg",&he_bg);
      catTree.Branch("he_sig",&he_sig);
    }
  // Load catalog data
  std::ifstream ifile("data_3fgl.csv");
  std::string line;
  std::ifstream vifile;
  std::string vline;

  int lineNumber = 0;

  while (std::getline(ifile, line))
    {
      lineNumber++;
      // skip header
      if(lineNumber > 1)   // if you only want, say, 10 sources, add && lineNumber <= 10 to the conditional
	{
	  std::istringstream iss(line); // construct a string stream from line
	  // read the tokens from current line separated by comma
	  std::vector<std::string> tokens; // here we store the tokens
	  std::string token; // current token
	  while (std::getline(iss, token, ','))
	    {
	      trim(token, '\"');
	      trim(token, '*');
	      tokens.push_back(token); // add the token to the vector
	    }

	  // map the tokens into our variables
	  sourceName = tokens[0];
	  association = tokens[1];
	  ra = atof(tokens[2].c_str());
	  dec = atof(tokens[3].c_str());
	  classType = tokens[4];
	  tevCat = tokens[5];
	  tevAssoc = tokens[6];

	  if(doLightCurves == true)
	    {	  
	      vifile = std::ifstream(TString::Format("./lightCurves/lightcurve_%s_%s.csv",tokens[2].c_str(),tokens[3].c_str()));
	      met.clear();
	      unixTime.clear();
	      le_count.clear();
	      le_bg.clear();
	      le_sig.clear();
	      he_count.clear();
	      he_bg.clear();
	      he_sig.clear();
	  
	      //vifile = std::ifstream("./lightCurves/lightcurve_21.3379_0.3612.csv");
	      while(std::getline(vifile,vline))
		{
		  // Grab the vectors from each corresponding light curve
		  std::istringstream viss(vline); // construct a string stream from line

		  // read the tokens from current line separated by comma
		  std::vector<std::string> vtokens; // here we store the tokens
		  std::string vtoken; // current token
		  while (std::getline(viss, vtoken, ','))
		    {
		      vtokens.push_back(vtoken); // add the token to the vector
		    }

		  metTemp = std::stoi(vtokens[0]); met.push_back(metTemp);
		  unixTimeTemp = std::stoi(vtokens[1]); unixTime.push_back(unixTimeTemp);
		  le_countTemp = std::stod(vtokens[2]); le_count.push_back(le_countTemp);
		  le_bgTemp = std::stod(vtokens[3]); le_bg.push_back(le_bgTemp);
		  le_sigTemp = std::stod(vtokens[4]); le_sig.push_back(le_sigTemp);
		  he_countTemp = std::stod(vtokens[5]); he_count.push_back(he_countTemp);
		  he_bgTemp = std::stod(vtokens[6]); he_bg.push_back(he_bgTemp);
		  he_sigTemp = std::stod(vtokens[7]); he_sig.push_back(he_sigTemp);
	  
		  //cout << "metTemp = " << metTemp << endl;
		  //cout << le_countTemp << endl;
	  
		}	  
	    }
	  catTree.Fill();
	  /*	  
		  std::cout << " " << std::endl;
		  std::cout << "Name: " << sourceName << std::endl;
		  std::cout << "dec: " << dec << std::endl;
		  std::cout << "ra: " << ra << std::endl;
		  std::cout << "type: " << classType << std::endl;
	  */
	  //std::cout << "association: " << sourceName << std::endl;
	  	  
	}
    }

  catTree.Write();
  f.Close();
  cout << "Done!" << endl;
  
  return;
}
