#include <vector>

#include "itensor/all.h"
#include "itensor/util/print_macro.h"

using namespace itensor;

/*
This file contains auxliary functions for automatically creating file names
*/

auto filename(auto VarNames, auto VarValsString) {
  string fname = VarNames[0] + "-" + VarValsString[0];
  for (unsigned int i = 1; i < VarNames.size(); i++)
    fname += "." + VarNames[i] + "-" + VarValsString[i];
  fname += ".dat";
  return fname;
}

// Exclude variables in positions from 0 to VarNames.size() in the filename
auto filenameExcludeParams(auto VarNames, auto VarValsString,
                           auto ExcludeList) {
  string fname;
  unsigned int j = 0;
  std::sort(ExcludeList.begin(), ExcludeList.end());
  if (ExcludeList[0] == 0) {
    while (ExcludeList[j] == j) j++;
    fname = VarNames[j] + "-" + VarValsString[j];
    for (unsigned int i = j + 1; i < VarNames.size(); i++) {
      if (i != ExcludeList[j])
        fname += "." + VarNames[i] + "-" + VarValsString[i];
      else
        j++;
    }
  } else {
    fname = VarNames[0] + "-" + VarValsString[0];
    for (unsigned int i = 1; i < VarNames.size(); i++) {
      if (i != ExcludeList[j])
        fname += "." + VarNames[i] + "-" + VarValsString[i];
      else
        j++;
    }
  }
  fname += ".dat";
  return fname;
}
