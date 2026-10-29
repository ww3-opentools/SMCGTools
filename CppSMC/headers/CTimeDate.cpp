/* Programs to print ctime stamps in model runs.   JGLi20Aug2026  */
#include <iostream>
#include <iomanip>
using namespace std;
#include <fstream>
#include <string>
#include <ctime>

void printTime(const std::string& msg = "",
               std::ostream& ostrm = std::cout) {
    std::time_t t = std::time(nullptr);
    if (!msg.empty()) ostrm << msg;
    ostrm << std::put_time(std::localtime(&t), "%H:%M:%S");
}

void printDateTime(const std::string& msg = "",
                   std::ostream& ostrm = std::cout) {
    std::time_t t = std::time(nullptr);
    if (!msg.empty()) ostrm << msg;
    ostrm << std::put_time(std::localtime(&t),
                        "%Y-%m-%d %H:%M:%S");
}

void printElapsed(std::time_t start, std::ostream& ostrm = std::cout) {
    std::time_t elapsed = std::time(nullptr) - start;
    ostrm <<" Total elapsed: "<< elapsed <<" s (";

    int h = elapsed / 3600;
    int m = (elapsed % 3600) / 60;
    int s = elapsed % 60;

    ostrm << h <<"h "<< m <<"m "<< s <<"s) \n";
}

