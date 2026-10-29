
#pragma once

/* Program header for ctime stamp printing functions in C++.  JGLi20Aug2026  */
#include <iostream>
#include <iomanip>
using namespace std;
#include <fstream>
#include <string>
#include <ctime>

void printTime(const std::string& msg = "",
               std::ostream& ostrm = std::cout);

void printDateTime(const std::string& msg = "",
                   std::ostream& ostrm = std::cout);

void printElapsed(std::time_t start, std::ostream& ostrm = std::cout); 

