
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

#ifndef __GRAPH_NODE_H__
#define __GRAPH_NODE_H__

#include <vector>
#include "GraphBasic.h"
#include "GraphEdge.h"

class GraphNode: public GraphBasic{
public:
  GraphNode(GraphEdge nNewEdge);
  GraphNode(int iPointTo);
  GraphNode(const char *sNodeName):GraphBasic(sNodeName){};
  GraphNode(const char *sNodeName, const char *sNodeLabel, char shape='t', char color='0');
  GraphNode(const char *sNodeName, int iPointTo);
  GraphNode(const char *sNodeName, GraphEdge objNewEdge);
  GraphNode(const char *sNodeName, const char *sNodeLabel, int iPointTo);
  ~GraphNode();

  int setEdge(GraphEdge nNewEdge);
  int setEdge(char* title, char type='e');
  int setEdge(const char *sEdgeName, int iPointTo);
  int setEdge(int iPointTo);
  int getNbEdges();
  GraphEdge getEdge(int iNum);
  void setShape(char shape){cShape = shape;}
  char getShape(){return cShape;}
  void setColor(char color){cColor = color;}
  char getColor(){return cColor;}


private:
  vector <GraphEdge> vEdge;
  char cShape;	// Shape of the node: 't': rectangle; 'r': rhomb; 'e': ellipse
  char cColor;	// Color of the node: 'b': blue; 'l': lightgreen; 'y': yellow; 'o': orange
};

typedef vector<GraphNode> VecGraphNode;

#endif
