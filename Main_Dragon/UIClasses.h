#ifndef UIClasses_h
#define UIClasses_h
#include <vector.h>
// The class for call graph 


//Yi Wen added on 07/01/2003. The global vlist array for
//call graph. 
//Each element in vlists is the actual procedure parameters set in some callsites 
extern vector<vector<int> > vlists;
extern vector<vector<int> > fake_edges;
extern bool compareVlist(int idx1,int idx2);

// C.H. Liao added on 3/4/2004, for feedback information of call graph
class CGEDGE
{
 public:
  unsigned int caller_id;
  unsigned int callee_id;
  unsigned int edge_id;
  unsigned int callsite_linenum;
  Widget Arc;
  float freq;
  CGEDGE();
  CGEDGE(unsigned int caller,
	 unsigned int callee,
	 unsigned int edge,
	 unsigned int line,
	 float frequency):
    caller_id(caller),
    callee_id(callee),
    edge_id(edge),
    callsite_linenum(line),
    freq(frequency)
    {}
  //  ~CGEDGE();

  CGEDGE& operator=(const CGEDGE * ed);
 void  Print(){
    cout<<caller_id<<"--"<<edge_id<<"-->"<<callee_id<< endl;
    cout<<"Callsite"<<callsite_linenum<<" Frequency:"<<freq<<endl;
  }

}; // end of class CGEDGE

class UIGraphNode
{
 public:
  UIGraphNode();
  ~UIGraphNode();
  void SetHasParent(bool value) { HasParent = value;};
  bool GetHasParent(){return HasParent;};
  void SetFilename(const char *);
  void SetName(char *);
  void SetLineNumber(unsigned int LineNum);
  void SetSuccs(vector<unsigned int>&);
  void SetSuccs();
  void SetSuccs(int *,int*, int *, int *, int *, int *,int*, int);
  void SetCallsites(int *,int);
  void SetFGSuccs(vector<long> Succs);
  int  GetNumFGSuccs() {return FGNumSuccs;};
  int  GetNumSuccs() {return NumSuccs;};

 
/* added by Yi Wen on 06/27/2003 for multiple edges*/
  int GetRealNumSuccs() const{return realTotalSuccs;};
  int GetFGSuccsID(int num) {return FGSuccs[num];};
  unsigned int GetSuccsID(int num){ 
	if(num<realTotalSuccs) return Succs[num];
	else{ cout<<"getSuccsID out of range\n"; return 0;} }
// added by Yi Wen ion 06/30/2003
  int GetEdgeID(int num) const{
	if(num<realTotalSuccs) return edge_idx[num];
	else{ cout<<"GetEdgeID out of range\n"; return -1;} }
  int GetCallsitesID(int num) const{
	if(num<realTotalSuccs) return callsites_idx[num];
	else{ cout<<"GetCallsitesID out of range\n"; return -1;} }
  int GetPreEdgeID(int num) const{
	if(num<realTotalSuccs) return pre_edge_idx[num];
	else{ cout<<"GetPreEdgeID out of range\n"; return -1;} }
  int GetVlistIndex(int num) const {
	if(num<realTotalSuccs) return vlist_index[num];
	else{ cout<<"GetVlistIndex out of range\n"; return -1;} }
  int GetFakeEdgeID(int num) const {
	if(num<realTotalSuccs) return fake_edges_idx[num];
	else{ cout<<"GetFakeEdgeID out of range\n"; return -1;} }
  int GetPreCallerID(int num) const {
	if(num<realTotalSuccs) return pre_caller_idx[num];
	else{ cout<<"GetPreCallerIndex out of range\n"; return -1;} }
  Widget GetCallgraphArcs(int num) const{
	if(num<realTotalSuccs) return cg_arc[num];
	else{ cout<<"GetCallgraphArcs out of range\n"; return NULL;} }
  void SetCallgraphArcs(int num,Widget arc) {cg_arc[num]=arc;}   

  void addParents(int p);
  //By Liao, for feedback information
  void addOutEdges(CGEDGE* ed){out_calledges.push_back(ed); };
  void addInEdges(CGEDGE* ed){in_calledges.push_back(ed);}
  void PrintFeedback();  //For debugging

  int GetParentsSize() const{return parents_idx.size();}
  int GetParents(int num) const{return parents_idx[num];}

  void SetIndex(int i){_index=i;}
  int GetIndex() const{return _index;}
//end of Yi added
  void SetUINode(Widget Wid) {UINode = Wid;};
  char * GetName(void) const{return Name;};
  Widget GetUINode(void){return UINode;};
  char * GetFilename() {return Filename;};
  unsigned int GetLineNumber(){return LineNum;};
  static unsigned int NumOfNodes;
  static unsigned int NumofMP;
  //By Liao
  //  static unsigned int CGfdback; 
  void Set_feedback(unsigned int wt,float freq, float cc)
{  
  weight=wt;
  caller_freq=freq;
  cycle=cc;
};
 


//Yi Wen added on 06/30/2003: to decide if the graph is constructed by new cg alg
  static bool NewCG;
//end of Yi Wen added
 

  bool GetVisited(void){return Visited;};
  void SetVisited(bool flag) {Visited = flag;};
  unsigned long GetFGNumOfNodes() {return FGNumOfNodes;};
  void SetFGNumOfNodes(unsigned long num){FGNumOfNodes = num;};  
  void AllocateFlowgraph(unsigned long num) {Flowgraph = new UIGraphNode[num];};
  UIGraphNode *GetFlowgraph() {return Flowgraph;};
  bool GetMP() {return Has_MP;};
  bool SetMP(bool flag) { Has_MP = flag;};
  int *FGSuccs;
  int FGNumSuccs;
   UIGraphNode *Flowgraph;
    Widget UINode;
 // Added by Liao for node feedback information
  vector <CGEDGE *> out_calledges;
  vector <CGEDGE *> in_calledges;
   unsigned int weight;
  float caller_freq;
  float cycle;
 private:
  int _index; //The index in GraphNode array
  bool Has_MP;
  bool Visited;
  bool HasParent;
  char *Filename;
  char *Name;
  unsigned int LineNum;
  int NumSuccs;
// added by Yi Wen on 06/27/2003 for multiple edges
  int realTotalSuccs;
  unsigned int *Succs;
// added by Yi Wen ion 06/30/2003
// add two variables for edge_index and previous edge index
//These two corresponding to Succs. e.g. Succs[0] has edge index _edge_idx[0]
//and its previous edge is _pre_edge_idx[0]  
  int *edge_idx,*callsites_idx,*pre_edge_idx, *vlist_index,*pre_caller_idx;
  int *fake_edges_idx;
  WidgetList cg_arc;
  vector<int> parents_idx;
 
 

  // end of feedback 
  int NumCallsites;
  unsigned int *Callsites; 
  
  unsigned long FGNumOfNodes; 
};




#endif
