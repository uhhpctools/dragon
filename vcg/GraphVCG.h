#ifndef __GRAPH_VCG_H__
#define __GRAPH_VCG_H__

class GraphVCG
{
 public:
  GraphVCG(const char *strGraph);
  GraphVCG();

  ~GraphVCG();

  void setName(const char *strGraph);
  char* getName();
  
 private:
  char *strName;
 
 
};

#endif

