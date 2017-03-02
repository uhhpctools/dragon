#include <X11/Xatom.h>
#include <X11/Intrinsic.h>
#include <X11/Shell.h>
#include <Xm/Xm.h> 
#include "UIClasses.h"

unsigned int UIGraphNode::NumOfNodes = 0;
unsigned int UIGraphNode::NumofMP = 0;
vector<vector<int> > vlists;
vector<vector<int> > fake_edges;
bool UIGraphNode::NewCG=True;

CGEDGE& 
CGEDGE::operator=(const CGEDGE *ed){
  caller_id=ed->caller_id;
  callee_id=ed->callee_id;
  edge_id=ed->edge_id;
  callsite_linenum=ed->callsite_linenum;
  freq=ed->freq;
 return *this;
}
void UIGraphNode::PrintFeedback()
{
  cout<<"PU:"<<Name<<"---------"<<endl;
  cout<<"Weight:"<<weight<<"freq"<<caller_freq<<"cycle"<<cycle<<endl;
  for(int i=0;i<out_calledges.size();i++)
    out_calledges[i]->Print();
}

UIGraphNode::UIGraphNode()
{
  Filename = NULL;
  Name = NULL;
  Succs = NULL;
  vlist_index=NULL;
  cg_arc=NULL;
  edge_idx=NULL;
  callsites_idx=NULL;
  pre_edge_idx=NULL; 
  pre_caller_idx=NULL; 
  LineNum = 0;
  FGNumSuccs = 0;
  NumSuccs = 0;
/* added by Yi Wen on 06/27/2003 for multiple edges*/
  realTotalSuccs=0;
  UINode = NULL;
  FGNumOfNodes = 0;
  Flowgraph = NULL;
  Visited = False;
  FGSuccs = NULL;
  Callsites = NULL;
  NumCallsites=0;
  HasParent = False;
}

UIGraphNode::~UIGraphNode()
{
  delete [] Filename; 
  delete [] Name; 
  delete [] Succs;
  delete [] vlist_index;
  delete [] cg_arc;
  delete [] edge_idx;
  delete [] callsites_idx;
  delete [] pre_edge_idx;
  delete [] pre_caller_idx;
  delete [] FGSuccs;
  delete [] Flowgraph;
  delete [] Callsites;
  for (int i=0;i<out_calledges.size();i++)
    if(!out_calledges[i]) delete out_calledges[i];
  for (int i=0;i<in_calledges.size();i++)
     if(!in_calledges[i])delete in_calledges[i];
}

void UIGraphNode::SetFilename(const char *Filename)
{
  if (UIGraphNode::Filename = new char [strlen (Filename)+1])
  {  strcpy(UIGraphNode::Filename, Filename); 

 
  } 
}

void UIGraphNode::SetName(char *Name)
{
  
  if( UIGraphNode::Name = new char [strlen (Name)+1])
    {
      strcpy(UIGraphNode::Name, Name);
      int len = strlen (UIGraphNode::Name);
  for (int i=0; i< len; i++)
    {
      if(UIGraphNode::Name[i]=='.' && UIGraphNode::Name[i+1]=='.')
	UIGraphNode::Name[i] = '\0';


    }  

    }
}

void UIGraphNode::SetLineNumber(unsigned int LineNum)
{
  UIGraphNode::LineNum = LineNum;

}


void UIGraphNode::SetSuccs(vector<unsigned int> &Succs)
{
  UIGraphNode::Succs = new unsigned int [Succs.size()];
  UIGraphNode::NumSuccs = Succs.size();
  for (int i = 0; i<Succs.size(); i++)
  {
    
     UIGraphNode::Succs[i] = Succs[i];   

  }
   
}
/* Yi Wen modified this function on 06/30/2003 for multiple edges
   The main change are: add a realTotalSuccs as total
   Made Succs has size total and fill in every elements.
   And also set edge_id and previous edge id
*/
void UIGraphNode::SetSuccs(int *Succs, int *edges, int* callsites, int *pre_edges, int* pre_caller, int* vlist_idx,int *fake_idx,int total)
{
  int *temp = new int [total];
/* added by Yi Wen on 06/27/2003 for multiple edges*/
  UIGraphNode::realTotalSuccs=total;
  // UIGraphNode::Succs = new unsigned int [total];
  //UIGraphNode::NumSuccs = total;

  for (int i = 0; i<total; i++)
  {
    temp[i] = Succs[i];
//	cout<<"Succs["<<i<<"]: "<<Succs[i]<<"\n";
    // UIGraphNode::Succs[i] = Succs[i];   

  }
   
  for (int i = 0; i<total-1; i++)
  {
    for (int j = i+1; j<total; j++)
      if (temp[i] == temp[j])
	   temp[j]=-1; 
  }
 
  int count = 0;

  for (int i = 0; i< total; i++)
    if (temp[i]!=-1) count++;

  UIGraphNode::NumSuccs = count;
	//cout<<"1\n";
  if (count > 0)
/* modified by Yi Wen on 06/27/2003 for multiple edges*/

    {   
	int index=0;
	int index2=total-1;
	UIGraphNode::Succs = new unsigned int [total];
         if(NewCG){  
		UIGraphNode::edge_idx = new int [total];
		UIGraphNode::callsites_idx = new int [total];
		UIGraphNode::pre_edge_idx = new int [total];
		UIGraphNode::pre_caller_idx = new int [total];
		UIGraphNode::vlist_index = new int [total];
		UIGraphNode::fake_edges_idx = new int [total];
	}
	//cout<<"2\n";
//initialize cg_index here:
       UIGraphNode::cg_arc = new Widget[total];
       for (int i=0; i<total; i++){
		if(temp[i]!=-1){
			UIGraphNode::Succs[index]=Succs[i];
         		if(NewCG){  
				UIGraphNode::edge_idx[index]=edges[i];
				UIGraphNode::callsites_idx[index]=callsites[i];
				UIGraphNode::pre_edge_idx[index]=pre_edges[i];
				UIGraphNode::pre_caller_idx[index]=pre_caller[i];
				UIGraphNode::vlist_index[index]=vlist_idx[i];
				UIGraphNode::fake_edges_idx[index]=fake_idx[i];
			}
			index++;
		}
		else{
			UIGraphNode::Succs[index2]=Succs[i];
         		if(NewCG){  
				UIGraphNode::edge_idx[index2]=edges[i];
				UIGraphNode::callsites_idx[index2]=callsites[i];
				UIGraphNode::pre_edge_idx[index2]=pre_edges[i];
				UIGraphNode::pre_caller_idx[index2]=pre_caller[i];
				UIGraphNode::vlist_index[index2]=vlist_idx[i];
				UIGraphNode::fake_edges_idx[index2]=fake_idx[i];
			}
			index2--;
		}
	  
	}	
    } 
//`  delete temp;

}


void UIGraphNode::addParents(int p)
{	
	bool find=False;
	for(int i=0;i<parents_idx.size();i++)
		find=(parents_idx[i]==p);
	//Debug by Liao
		if(!find)
		parents_idx.push_back(p);
}
/*return true is two vlists are the same, false otherwise
*/
bool compareVlist(int idx1,int idx2){
	if(vlists[idx1].size()!=vlists[idx2].size())
		return False;
	int size=vlists[idx1].size();
	cout<<"idx1: "<<idx1<<", idx2"<<idx2<<"\n";
	for(int i=0;i<size;i++){
		cout<<"1: "<<vlists[idx1][i]<<" ,2: "<<
		vlists[idx2][i]<<"\n";
		if(vlists[idx1][i]!=vlists[idx2][i])
			return False;
	}
	return True;
}

void UIGraphNode::SetCallsites(int *Succs, int total)
{
  UIGraphNode::Callsites = new unsigned int [total];
  UIGraphNode::NumCallsites = total;
  for (int i = 0; i<total; i++)
  {
    
     UIGraphNode::Callsites[i] = Succs[i];   

  }

  

   
}

void UIGraphNode::SetSuccs()
{
  bool *unique = new bool [NumCallsites];
  for (int i =0;i<NumCallsites;i++)
    unique[i] = True;  
     
  for (int i =0; i<NumCallsites;i++)
    {  
      // needs improvement algorithm
      unsigned int val = Callsites[i];
     
      for (int j=i+1; j<NumCallsites; j++)
      {
        if (val == Callsites[j])
          unique[j]=False;    
       
      } 


    }

  int total = 0;
  for(int i = 0; i<NumCallsites;i++)
    if(unique[i]==True) total++;

  if (total>0) { 
  UIGraphNode::Succs = new unsigned int [total];
  UIGraphNode::NumSuccs = total;
 

  for (int i = 0; i<NumCallsites; i++)
  {
    if(unique[i]==True)
     UIGraphNode::Succs[i] = Callsites[i];   

  }
  } 
  delete [] unique;
   
}



void UIGraphNode::SetFGSuccs(vector<long> Succs)
{
  UIGraphNode::FGSuccs = new int [Succs.size()];
  UIGraphNode::FGNumSuccs = Succs.size();
  for (int i = 0; i<Succs.size(); i++)
  {
    int num;
    num = Succs[i];
     UIGraphNode::FGSuccs[i] = Succs[i];   

  }
   
}












