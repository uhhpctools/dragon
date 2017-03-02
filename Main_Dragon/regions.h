#include <string>
extern void create_InterRegionShell (Widget parent);
extern void CloseRegions(Widget, XtPointer,XtPointer);
extern void CreateRegions(void);

class pregions
{
 public:
    
  int pregion;
  string lower;
  string upper;
  string stride;

  pregions() {pregion =-1;
              lower="";
              upper="";
              stride="";};

};



class regions
{
 public:
  string procname;
  string variablename;
  string usage;
  string prstart;
  string prend; 
  string dimensions;
  int start;
  int end;
  regions() {procname =""; 
             variablename="";
             usage ="";
             prstart ="-1";
             prend = "-1";
             dimensions = "-1";
             int start =-1;
             int end =-1; 
  };

};
