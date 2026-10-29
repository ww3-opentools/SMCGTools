#ifndef CPP_RESTRFILE_H
#define CPP_RESTRFILE_H

#include <iostream>
using namespace std;
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>

#include "WvSpectr.h"

void RestrRead( WvSpectr& Spc, int& NSt, string RstrFile );
void RestWrite( WvSpectr& Spc, int NTS, string Flnm, string& Msgs);

#endif

