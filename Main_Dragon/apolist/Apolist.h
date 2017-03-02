#ifndef __APOLIST__
#define __APOLIST__
//-----------------
//// Header
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>

//-----------------
//// Constant Definition
//
////-----------------
//// Type Definition
using namespace std;
  
/*
 ** Create a type definition for vector of string
 **/ 

typedef vector<string> Vecstring;

typedef struct{
  string status;  // status of parallelization
  int lineno;     // the start line number of the loop
  Vecstring info; // list of information/messages
} Lineinfo;

/**
 ** Structure of line of information
 **/

typedef vector<Lineinfo> Veclineinfo;

/**
 ** Type of vector of Lineinfo
 **/

typedef struct{
  Veclineinfo lineinfo;
} Proginfo;

/**
 ** Main datatype for information of a subprogram
 **/

typedef map<string, Proginfo> MapProginfo;

/**
 ** Function prototype
 **/ 

int loadFile(const char *filename, MapProginfo &mapProginfo);

#endif
