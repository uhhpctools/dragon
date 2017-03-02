#ifndef __GRAPH_Basic_H__
#define __GRAPH_Basic_H__

#include <stdlib.h>

#define MAX_TITLE_CHARS 80

class GraphBasic
{
 public:
  GraphBasic(const char *strGraph);
  GraphBasic();

  ~GraphBasic();

  void setName(const char *strGraph);
  char* getName(void);
  void setLabel(const char *strGraph);
  char* getLabel(void);
  
 private:
  char strName[5];
 
 protected:
  char strLabel[MAX_TITLE_CHARS];
 
};

#endif

