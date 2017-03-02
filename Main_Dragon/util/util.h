#ifndef __UTIL_H__
#define __UTIL_H__

#include <vector>
#include <string>

//using namespace std;

// lines_t is the type vector of strings
typedef vector< std::string > lines_t;
// function to retrieve the basename of a filename


string getfilebasename(string filename);

#endif
