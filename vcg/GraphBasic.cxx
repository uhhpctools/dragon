/*
 * Copyright (c) 2002
 * University of Houston.
 *
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.  University of Houston makes no
 * representations about the suitability of this software for any
 * purpose.  It is provided "as is" without express or implied warranty.
 */

#include <stdlib.h>
#include <string.h>
#include "GraphBasic.h"

// --------------------------------------------
// default constructor
// --------------------------------------------
GraphBasic::GraphBasic() {
  strName[0]='\0';
}

// --------------------------------------------
// constructor with assigned graph name
// --------------------------------------------
GraphBasic::GraphBasic(const char *strGraph) 
{
  setName(strGraph);
}

// --------------------------------------------
// function to set name
// --------------------------------------------
void GraphBasic::setName(const char *strGraph)
{
    strcpy(strName,strGraph);
}

// --------------------------------------------
// function to set name
// --------------------------------------------
char* GraphBasic::getName(void)
{
  return strName;
}

// --------------------------------------------
// function to set Label
// --------------------------------------------
void GraphBasic::setLabel(const char *strGraph)
{
    strcpy(strLabel,strGraph);
}

// --------------------------------------------
// function to set Label
// --------------------------------------------
char* GraphBasic::getLabel(void)
{
  return strLabel;
}


// --------------------------------------------
// destructor
// --------------------------------------------
GraphBasic::~GraphBasic()
{
}
