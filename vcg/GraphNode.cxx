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

#include <vector>
#include <stdlib.h>
#include "GraphNode.h"

// constructor with an edge
GraphNode::GraphNode(GraphEdge nNewEdge):GraphBasic()
{
  setEdge(nNewEdge);
}

// the same: allocate an immediate edge
GraphNode::GraphNode(int iPointTo):GraphBasic()
{
  setEdge(iPointTo);
}

// allocate the title, label
GraphNode::GraphNode(const char *sNodeName, const char *sNodeLabel, char shape='t', char color='0'):GraphBasic(sNodeName)
{
  strcpy(strLabel, sNodeLabel);
  cShape = shape;
  cColor = color;
}

// allocate the title plus the edge
GraphNode::GraphNode(const char *sNodeName, int iPointTo):GraphBasic(sNodeName)
{
  GraphEdge objEdge(iPointTo);
  setEdge(objEdge);  
}

// allocate the title, label, plus the edge
GraphNode::GraphNode(const char *sNodeName, const char *sNodeLabel, int iPointTo):GraphBasic(sNodeName)
{
  strcpy(strLabel, sNodeLabel);
  GraphEdge objEdge(iPointTo);
  setEdge(objEdge);  
}

// allocate a title and an edge
GraphNode::GraphNode(const char *sNodeName, GraphEdge objNewEdge):GraphBasic(sNodeName)
{
  setEdge(objNewEdge);
}

// destructor
GraphNode::~GraphNode()
{
  vEdge.clear();
}

// create a new edge
int GraphNode::setEdge(GraphEdge objNewEdge)
{
  vEdge.push_back(objNewEdge);
  return vEdge.size();
}

  
// create a new edge
int GraphNode::setEdge(char* title, char type='e')
{
  GraphEdge objEdge(title, type);
  return setEdge(objEdge);
}


// create a new edge
int GraphNode::setEdge(int iPointTo){
  GraphEdge objEdge(iPointTo);
  return setEdge(objEdge);
}

// create a new edge with its name
int GraphNode::setEdge(const char *sEdgeName, int iPointTo) {
  GraphEdge objEdge(sEdgeName, iPointTo);
  return setEdge(objEdge);
}

// get number of edges
int GraphNode::getNbEdges() {
  return vEdge.size();
}

// get an edge
GraphEdge GraphNode::getEdge(int iNum) {
  return vEdge[iNum];
}


