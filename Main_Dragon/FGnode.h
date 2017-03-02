/************************************************
*						*
* 	a class of Control Flow Graph Node 	*
*	Lei Huang 07/31/02			*
*						*
************************************************/



#ifndef FGnode_h
#define FGnode_h


#include <X11/Xatom.h>
#include <X11/Intrinsic.h>
#include <X11/Shell.h>
#include <Xm/Xm.h> 
#include <vector.h>
#include <fstream.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>

// #include "opt_bb.h"
#include "feedback.h"
#define FALSE 0
#define TRUE  1

//----------------------
// Added by C.H. Liao, 3/2/04
// To support feedback information about control flow graph edges
// Using OPT_FB_EDGE as template here
//FOR DRAGON side
class FGEDGE {

 public:
  
  IDTYPE       identity; 
  IDTYPE       source;
  IDTYPE       destination;

  FB_EDGE_TYPE edge_type;
  FB_FREQ      freq;
  //--------------------------
  FGEDGE() {};
  FGEDGE(IDTYPE src, IDTYPE dst): 
    identity(0),
    source(src),
    destination(dst){
    edge_type=FB_EDGE_UNINIT;
    freq._type=FB_FREQ_TYPE_UNINIT;
    freq._value=-1.0;
  }
  FGEDGE(IDTYPE id, IDTYPE src, IDTYPE dst,
	       FB_EDGE_TYPE typ  = FB_EDGE_UNINIT,
	       FB_FREQ      freq = FB_FREQ_UNINIT ) :
    identity(id),
    source(src),
    destination(dst),
    edge_type(typ),
    freq(freq) {}
  // The size of edge is crucial for calculating the size of 
  //Control flow graph written into a binary file!

  unsigned GetSize(){ return 3*sizeof(IDTYPE)+sizeof(float)
+sizeof(FB_EDGE_TYPE)+sizeof(FB_FREQ_TYPE);
  }
  void Print(FILE *fp ) const 
    {
      char buffer[FGEDGE_TYPE_NAME_LENGTH];
      FB_EDGE_TYPE_sprintf( buffer, edge_type );

      fprintf( fp, "Edge[%3d]:  (%3d --> %3d) : freq = ",
	       identity, source, destination );
      freq.Print( fp );
      fprintf( fp, " : %s\n", buffer );
    } 

  // write the content of this object into a file
  void Print(ofstream &outfile){

 outfile.write((char *) &identity, sizeof(IDTYPE)); 	
 outfile.write((char *) &source, sizeof(IDTYPE));
 outfile.write((char *) &destination, sizeof(IDTYPE));

 outfile.write((char *) &edge_type, sizeof(FB_EDGE_TYPE));
 outfile.write((char *) &freq._type, sizeof(FB_FREQ_TYPE));
 outfile.write((char *) &freq._value, sizeof(float)); 

  }// END of Print

  FGEDGE (ifstream &fin) // Construct an object from file stream
    {
     
  FB_FREQ_TYPE f_type;
      float f_value;

  fin.read((char *)&identity, sizeof(IDTYPE));  // get ID
  fin.read((char *)&source, sizeof(IDTYPE));
  fin.read((char *)&destination, sizeof(IDTYPE));

  fin.read((char *) &edge_type, sizeof(FB_EDGE_TYPE));
  fin.read((char *) &f_type, sizeof(FB_FREQ_TYPE));
  fin.read((char *) &f_value, sizeof(float));

  freq=FB_FREQ(f_type, f_value);


    }// end of constructor from file stream

  void Print(){
    char buffer[FGEDGE_TYPE_NAME_LENGTH];
    FB_EDGE_TYPE_sprintf( buffer, edge_type );

    // cout<<"An edge of CFG:"<<endl;
    cout<<source<<"--"<<identity<<"-->"<<destination<<endl;
    cout<<"Edge type"<<buffer<<endl;
    freq.Sprintf(buffer);
    cout <<"Freqency:"<<buffer<<endl;
  }
};


class FGNODE
{
      private:

	 char Label[20];	
	 char Kind[20];
	 unsigned int ID;	   // basic block ID
//	 vector<unsigned int>preds;  // parents
//	 vector<unsigned int>succs;  // childen
	 char visited;
	 Boolean redundant;

      public:
	
	vector<unsigned int> line_nums;
	vector<int>preds;  // parents
	vector<int>succs;  // childen
	
	vector <FGEDGE *> in_edge_vec;  //Added by Liao
	vector <FGEDGE *> out_edge_vec;

	Widget UINode;
	FGNODE(){ visited = FALSE; redundant=FALSE; }
	FGNODE(unsigned int id, char* kd){
		ID = id;
		if(strcmp(kd,"GOTO")==0)
			strcpy(Label, "STMTS");
		else if(strcmp(kd, "DOSTEP")==0)
			strcpy(Label, "DOEND");
		else if(strcmp(kd, "WHILEEND")==0)
			strcpy(Label, "LOOP");
		else if(strcmp(kd, "LOGIF")==0)
			strcpy(Label, "BRANCH");
		else 
			strcpy(Label,kd);		
		strcpy(Kind,kd);
		visited = FALSE;
		redundant = FALSE;
	}
	void safe_clear();
	~FGNODE() { 	  
	  safe_clear();
	  
	};
	void Insert_line_num(int line_num){
		line_nums.push_back(line_num);
	}
	void Insert_succ(unsigned int id){
		Boolean found=FALSE;
		for(int i=0; i<succs.size(); i++)
			if(succs[i]==id){	// already exists
				found=TRUE;
				break;
			}
		if(!found)			
			succs.push_back(id);
	}
	void Insert_pred(unsigned int id){
		Boolean found=FALSE;
		for(int i=0; i<preds.size(); i++)
			if(preds[i]==id){	// already exists
				found=TRUE;
				break;
			}
		if(!found)			
			preds.push_back(id);
	}
	// Liao
	void Insert_succ(FGEDGE *o_vec){
	  out_edge_vec.push_back(o_vec);
	}
	void Insert_pred(FGEDGE *i_vec){
	  in_edge_vec.push_back(i_vec);
	}
	// functions to help change frequency after node deletion
	int changeOutEdge(IDTYPE old_node, IDTYPE new_node);
	FGEDGE * getOutEdge(IDTYPE src, IDTYPE dst);

	char* Get_Kind() {return Kind; }
	void  Set_Kind(char* kd) {strcpy(Kind,kd);}
	char* Get_Name() { return Label; }
	void  Set_Name(char* name) {strcpy(Label, name);}
	unsigned int Get_ID()   { return ID;    }
	char Is_Visited() { return visited; }
	void Set_Visited(char visit) { visited = visit; }
	int GetSizeofSuccs() { return succs.size(); }
	int GetSizeofPreds() { return preds.size(); }
	int GetSizeofLinenums() { return line_nums.size(); }
	int GetIDofSuccs(int num) { return succs[num]; }
	char Get_Redundant() { return redundant; }
	char Set_Redundant(Boolean rd) {redundant = rd;}
	void Print_pred();
	void Print_succ();
	void Print_line_num();
	void Print();
	void Print_pred(FILE *file);
	void Print_succ(FILE *file);
	void Print_line_num(FILE *file);
	void Print(FILE *file);
	void Print(ofstream &outfile);
	unsigned int GetSize();
	void SetUINode(Widget Wid) {UINode = Wid;}
	Widget GetUINode() {return UINode;}

	FGNODE& operator= (FGNODE &fnode);
};

/********************************************************
*							*
*	Flow Graph class contains Flow graph nodes 	*
*							*
********************************************************/

class FGRAPH
{

	private:
		char  pu_name[50];
		char  Finished;


	public:
		vector<FGNODE> FG_Nodes;
		FGRAPH(){Finished = FALSE;}
		FGRAPH(char* name){
//			pu_name = CXX_NEW(char[strlen(pu_name)+1],Malloc_Mem_Pool);
			strcpy(pu_name, name);
			Finished = FALSE;
		}

		~FGRAPH(){
//			CXX_DELETE(pu_name,Malloc_Mem_Pool);
		  for(int i=0;i<FG_Nodes.size();i++)
		    FG_Nodes[i].safe_clear();
		  
		}
		void Insert_FNode(FGNODE node);
		char If_Finished(){ return Finished; }
		void Set_Finished(char finish){ Finished = finish; }
		void Set_pu_name(char* name){
//			pu_name = CXX_NEW(char[strlen(pu_name)+1],Malloc_Mem_Pool);
			strcpy(pu_name, name);
		}
		void Print(){
			cout<< "\n ---- " << pu_name <<" ------------\n";
			for(int i=0; i<FG_Nodes.size(); i++){
				FG_Nodes[i].Print();
				
			}
		}
		void Print(FILE *file);
		void Print(ofstream &outfile);
		unsigned int GetSize();
		void Load_CFG(char* filename, char* pu);
		int GetSizeofFGnodes() { return FG_Nodes.size(); }
		void Modify_CFG();
		void Delete_Node(FGNODE &fnode);  // delete a specific BB 
		void PrintVCG();
};

//extern vector<FGRAPH, mempool_allocator<FGRAPH> > CFGRAPH;
//extern FGRAPH Curr_CFG;

extern FGRAPH cfg;
#endif
