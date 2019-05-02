#include "bvv-util.h"

// cpp'ed version of https://stackoverflow.com/questions/9210528/split-string-with-delimiters-in-c
std::vector<std::string> str_split(char* a_str, const char a_delim)
{
  char delim[2];
  delim[0] = a_delim;
  delim[1] = 0;

  std::vector<std::string> result;

  char* token = strtok(a_str, delim);

  while (token) {
    result.push_back(token);
    token = strtok(0, delim);
  }

  return result;
}


void StrArrPrint(std::vector<std::string> str){
  for (int ind = 0; ind < int(str.size()); ind++) {
    std::cout << str[ind] << std::endl;
  }
}

