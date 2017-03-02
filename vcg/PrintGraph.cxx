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
#include <unistd.h>
#include <stdio.h>

#include <vector>
#include <iostream.h>
#include <fstream.h>
#include <string.h>

#include "GraphNode.h"
#include "PrintGraph.h"

// create an object with a node root and output file name
PrintGraph::PrintGraph(GraphNode objRoot, const char *sFilename="") {
  vNodes.push_back(objRoot);
  setFile(sFilename);
}

// destructor
PrintGraph::~PrintGraph() {
//  if (bHasBeenAllocated) {
//    std::cout = objTmp;
//  }
  vNodes.clear();
}

// set file name
void PrintGraph::setFile(const char *sFilename) {
  if (sFilename[0]!='\0') {
    objOut.open(sFilename, ofstream::out );
    // objTmp = std::cout;
//    std::cout = objOut;
    bHasBeenAllocated = (objOut != NULL);
  } 
}

// add a new node
void PrintGraph::addNode(GraphNode objNode) {
  vNodes.push_back(objNode);
}

//add a new node with its pointer
void PrintGraph::addNode(GraphNode *ptrNode) {
  vNodes.push_back(*ptrNode);
}

void PrintGraph::addNode(const char *strNode) {
  GraphNode objNode(strNode);
  addNode(strNode);
}

void PrintGraph::addNode(const char *strNode, int iPointTo) {
  GraphNode objNode(strNode, iPointTo);
  addNode(objNode);
}

int PrintGraph::printNodes(const char* sFilename="", const char *sTitle=""){
  unsigned char bOk=1;
  int i;
  setFile(sFilename);

  // header
  objOut << "graph: { " << endl;
  if (sTitle[0] != '\0')
    objOut << "title: \"" <<sTitle<<"\"" << endl;

  // list if nodes
  VecGraphNode::iterator iNodeIter;

  for(iNodeIter=vNodes.begin(); iNodeIter != vNodes.end(); iNodeIter++) {
    objOut << "node: { title:\"" << iNodeIter->getName() << "\"" << endl;
    objOut << "	    label: \"" << iNodeIter->getLabel() << "\"" ;
    switch(iNodeIter->getShape()){
    case 'r':
	objOut << " shape: rhomb ";
	break;
    case 'e':
	objOut << " shape: ellipse borderwidth: 5 ";
	break;
    }	

    switch(iNodeIter->getColor()){
    case 'y':
	objOut << " color: yellow";
	break;
    case 'b':
	objOut << " color: lightblue";
	break;
    case 'l':
	objOut << " color: lightgreen";
	break;
    case 'o':
	objOut << " color: orange";
	break;
    }
    objOut <<  "}" << endl;
  }

  // list of edges

  char *ptrName;
  for(iNodeIter=vNodes.begin(); iNodeIter != vNodes.end(); iNodeIter++) {

    for (int i=0;i<iNodeIter->getNbEdges();i++) {
      // edge 
      switch(iNodeIter->getEdge(i).getType()){
      case 'b':
      	objOut << "backedge: { ";
	break;
      case 't':
      	objOut << "bentnearedge: { ";
	break;
      default:
      	objOut << "edge: { ";
	break;
      }

      ptrName = iNodeIter->getEdge(i).getName();
      if (*ptrName != '\0')
	objOut <<"title: \"" << ptrName << "\" " ;
      objOut << "sourcename: \""<< iNodeIter->getName() << "\" targetname: \"" << iNodeIter->getEdge(i).getSuccs() << "\" }" << endl;
    }
  }

  // footer
  objOut << "}" << endl;
  objOut.close();
  
  //char argv[][20]={"xvcg","-psoutput",sFilename,"",'\0'};
//  execlp("xvcg", "-psoutput",sFilename);
//  if(execlp("xvcg", "xvcg", "vcg")<0)
//	printf("execute xvcg error!\n");
//  cout << "errorno: " << execlp("xvcg",sFilename) <<endl;
  
  return bHasBeenAllocated;
}


