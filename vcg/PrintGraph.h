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

#ifndef __PRINTGRAPH_H__
#define __PRINTGRAPH_H__

#include <ostream.h>
#include <fstream.h>
#include <vector>
#include "GraphNode.h"

class PrintGraph {
 public:
  PrintGraph() {};
  PrintGraph(GraphNode objRoot, const char *sFilename="");
  PrintGraph(const char *sFilename);
  void addNode(GraphNode objNode);
  void addNode(const char *strNode);
  void addNode(const char *strNode, int iPointTo);
  void addNode(GraphNode *ptrNode);
  void setFile(const char *sFilename);
  int printNodes(const char* sFilename="", const char *sTitle="");
  ~PrintGraph();

 private:

  ofstream objOut;
  ofstream objTmp;
  vector<GraphNode> vNodes;
  unsigned char bHasBeenAllocated;

};

#endif
