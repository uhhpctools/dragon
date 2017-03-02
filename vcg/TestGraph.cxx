#include <stdlib.h>
#include <unistd.h>

#include "GraphNode.h"
#include "PrintGraph.h"

int main() 
{
  PrintGraph myPrint;
  GraphNode myNode("one",1);
  myPrint.addNode(myNode);
  cout << "test one" << endl;
  myPrint.addNode(new GraphNode("two",2));
  myPrint.addNode(new GraphNode("three",3));
  myPrint.addNode(new GraphNode("four",4));
  myPrint.addNode(new GraphNode("five",5));
  myPrint.addNode(new GraphNode("six"));
  myPrint.printNodes("vcg","toto");

  execlp("xvcg", "xvcg", "vcg");

}
