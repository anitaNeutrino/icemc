#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

void rootify3FGL()
{

  TFile f("3FGL.root","RECREATE");
  TTree catTree("catTree","3FGL catalog");

  // CATALOG VARS:
  //'Source_Name','DEJ2000','RAJ2000','GLON','GLAT','CLASS1','ASSOC1'
  
  std::string sourceName;
  double dec;
  double ra;
  double glon;
  double glat;
  std::string classType;
  std::string association;
  
  catTree.Branch("sourceName",&sourceName);
  catTree.Branch("dec",&dec,"dec/D");
  catTree.Branch("ra",&ra,"ra/D");
  catTree.Branch("glon",&glon,"glon/D");
  catTree.Branch("glat",&glat,"glat/D");
  catTree.Branch("classType",&classType);
  catTree.Branch("association",&association);
  
  // Load catalog data
  std::ifstream ifile("data_3fgl.csv");

  std::string line; 
  while (std::getline(ifile, line))
    {
      
      std::istringstream iss(line); // construct a string stream from line

      // read the tokens from current line separated by comma
      std::vector<std::string> tokens; // here we store the tokens
      std::string token; // current token
      while (std::getline(iss, token, ','))
        {
	  tokens.push_back(token); // add the token to the vector
        }

      // map the tokens into our variables
      sourceName = tokens[0];
      dec = std::stod(tokens[1]);
      ra = std::stod(tokens[2]);
      glon = std::stod(tokens[3]);
      glat = std::stod(tokens[4]);
      classType = tokens[5];
      association = tokens[6];

      catTree.Fill();
      /*
      std::cout << " " << std::endl;
      std::cout << "Name: " << sourceName << std::endl;
      std::cout << "dec: " << dec << std::endl;
      std::cout << "ra: " << ra << std::endl;
      std::cout << "glon: " << glon << std::endl;
      std::cout << "glat: " << glat << std::endl;
      std::cout << "type: " << classType << std::endl;
      std::cout << "association:" << association << std::endl;
      */
    }

  catTree.Write();
  std::cout << "Catalog written to tree." << std::endl;

  return;
}
