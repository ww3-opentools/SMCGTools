#ifndef CPP_PRINTVECTOR_H
#define CPP_PRINTVECTOR_H

#include <iostream>
using namespace std;
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>

void PrintVector(const std::vector<float>& Vctr, int NCel, int NTS, string FN, const string& Msgs = "");

#endif

