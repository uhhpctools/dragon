#include <stdlib.h>
#include <string.h>
#include "GraphVCG.h"

// --------------------------------------------
// default constructor
// --------------------------------------------
GraphVCG::GraphVCG() {
  strName=0;
}

// --------------------------------------------
// constructor with assigned graph name
// --------------------------------------------
GraphVCG::GraphVCG(const char *strGraph) 
{
  setName(strGraph);
}

// --------------------------------------------
// function to set name
// --------------------------------------------
void GraphVCG::setName(const char *strGraph)
{
  if (strName != NULL) {
    delete strName;
  }
  //try{
    strName = new char[strlen(strGraph)+1];
    strcpy(strName,strGraph);
    //}

}

// --------------------------------------------
// function to set name
// --------------------------------------------
inline char* GraphVCG::getName()
{
  return strName;
}


// --------------------------------------------
// destructor
// --------------------------------------------
GraphVCG::~GraphVCG()
{
  if (strName!=NULL) {
    delete strName;
  }
}
