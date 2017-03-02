#ifndef __GRAPH_EDGE_H__
#define __GRAPH_EDGE_H__

#include "GraphBasic.h"

class GraphEdge : public GraphBasic{
 public:
  GraphEdge(int iPoint);
  GraphEdge(const char* str, int iPoint);
  GraphEdge(const char* Succs, char type='e');
  void setPoint(int iTo);
  int getPoint(void);
  void setSuccs(const char* succs);
  char* getSuccs(void);
  void setType(char type);
  char getType(void);

 private:
  int iPointTo;
  char strSuccs[5];	// the title of successor
  char cType;	//	the type of edge : b: backedge, e: normal edge, t: bent edge
};

#endif
