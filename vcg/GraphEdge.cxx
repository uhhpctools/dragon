#include <string.h>

#include "GraphBasic.h"
#include "GraphEdge.h"

//------------------------
// constructor edge
//------------------------
GraphEdge::GraphEdge(int iPoint):GraphBasic() {
  setPoint(iPoint);
}

GraphEdge::GraphEdge(const char *Succs, char type='e'):GraphBasic() {
  setSuccs(Succs);
  cType = type;
}

//------------------------
// yet another constructor
//------------------------
GraphEdge::GraphEdge(const char *str, int iPoint):GraphBasic(str) {
  setPoint(iPoint);
}

//------------------------
// to set a point
//------------------------
inline void GraphEdge::setPoint(int iTo) {
  iPointTo = iTo;
}

//------------------------
// to get a point
//------------------------
int GraphEdge::getPoint(void) {
  return iPointTo;
}

//------------------------
// to set a title of successor
//------------------------
inline void GraphEdge::setSuccs(const char* succs) {
  strcpy(strSuccs,succs);
}

//------------------------
// to get a title of successor
//------------------------
char* GraphEdge::getSuccs(void) {
  return strSuccs;
}

//------------------------
// to set a type of edge
//------------------------
inline void GraphEdge::setType(char type) {
  cType = type;
}

//------------------------
// to get a type of edge
//------------------------
char GraphEdge::getType(void) {
  return cType;
}
