/************************************************************************
*									*
*									*
* 		Control Flow graph class  07/31/2002  Lei Huang 	*
*									*
*	This source code contains some methods of Control Flow graph 	*
*	Nodes Class and Control Flow Graph Class. 			*
*									*
*	FGRAPH::Modify_CFG() is a method to modify the CFG that Open64  *
*	dumped so that the CFG could match source codes more precisely  *
*	and standardly.	It deletes some nodes without any line_no, and	*
*	modify loop CFG
*									*
************************************************************************/

#include <stdlib.h>
#include <unistd.h>

#include <vector.h>
#include "FGnode.h"
#include "vcg/GraphNode.h"
#include "vcg/PrintGraph.h"

const char *FB_EDGE_NAMES[]
    = { "------", "INCOMING", "OUTGOING", "ENTRY_OUTGOING",
	"BRANCH_TAKEN", "BRANCH_NOT_TAKEN",
	"LOOP_ZERO", "LOOP_POSITIVE", "LOOP_OUT", "LOOP_BACK",
	"LOOP_EXIT", "LOOP_ITERATE",
	"CIRCUIT_LEFT", "CIRCUIT_RIGHT", "CIRCUIT_NEITHER",
	"CALL_INCOMING", "CALL_OUTGOING", "CALL_INOUTSAME",
	"IO_OUTGOING", "IO_ESCAPE[1]", "IO_ESCAPE[2]", "IO_ESCAPE[3]",
	"SWITCH_DEFAULT" };

INT
 FB_EDGE_TYPE_sprintf( char *buffer, const FB_EDGE_TYPE fb_type )
{
 if ( fb_type < FB_EDGE_SWITCH_BASE ) {
    return sprintf( buffer, "%s", FB_EDGE_NAMES[fb_type] );
  } else {
    return sprintf( buffer, "SWITCH[%d]", fb_type - FB_EDGE_SWITCH_BASE );
  }
}

/*
FGNODE&
FGNODE::operator= (const BB_NODE *bb_node){

	ID = bb_node->Id();
	strcpy(Label, bb_node->Kind_name());
	return *this;
}
*/

// copy a Flow Graph Node
FGNODE&
FGNODE::operator= (FGNODE &fnode){

  FGEDGE *tempEd;
	ID = fnode.Get_ID();
	strcpy(Label, fnode.Get_Name());
	strcpy(Kind,  fnode.Get_Kind());

	line_nums.clear();
	for(int i=0; i<fnode.GetSizeofLinenums(); i++)
		Insert_line_num(fnode.line_nums[i]);

	preds.clear();
	in_edge_vec.clear();
	out_edge_vec.clear();
	for(int i=0; i<fnode.GetSizeofPreds(); i++)
	  {
		Insert_pred(fnode.preds[i]);

		tempEd=new FGEDGE();
		(*tempEd)=*(fnode.in_edge_vec[i]);
		Insert_pred(tempEd);
		  
	  }
	succs.clear();
	for(int i=0; i<fnode.GetSizeofSuccs(); i++)
	  {		
	    Insert_succ(fnode.succs[i]);

	    tempEd=new FGEDGE();
	    (*tempEd)=*(fnode.out_edge_vec[i]);
	    Insert_succ(tempEd);
	  }
	return *this;
}
void FGNODE::safe_clear(){
   for (int i=0;i<in_edge_vec.size();i++)
	   if(!in_edge_vec[i]) delete in_edge_vec[i];
	  for (int j=0;j<out_edge_vec.size();j++)
	    if(!out_edge_vec[j]) delete out_edge_vec[j];
}

void 
FGNODE::Print_pred(){
  for(int i=0; i<preds.size(); i++){
			cout<<preds[i]<<" ";
			in_edge_vec[i]->Print();
  }
}
	
void
FGNODE::Print_succ(){
  for(int i=0; i<succs.size(); i++){
	cout<<succs[i]<<" ";
	out_edge_vec[i]->Print();
  }
}
	
void
FGNODE::Print_line_num(){
		for(int i=0; i<line_nums.size(); i++)
			cout<<line_nums[i]<<" ";
}

void 
FGNODE::Print_pred(FILE *file){
		fprintf(file,"Preds: ");
		for(int i=0; i<preds.size(); i++)
		  {
			fprintf(file, "%i ", preds[i]);
			in_edge_vec[i]->Print(file);
		  }
		fprintf(file,"\n");
}
	
void
FGNODE::Print_succ(FILE *file){
		fprintf(file,"Succs: ");
		for(int i=0; i<succs.size(); i++){
			fprintf(file,"%i ", succs[i]);
			out_edge_vec[i]->Print(file);
		}
		fprintf(file,"\n");
}
	
void
FGNODE::Print_line_num(FILE *file){
		fprintf(file,"Line_no: ");
		for(int i=0; i<line_nums.size(); i++)
			fprintf(file, "%i ", line_nums[i]);
		fprintf(file,"\n");
}

	
void
FGNODE::Print(){
  cout <<"------------------------------"<<endl;
  cout<<"\nNode ID=="<<ID<<", Label=="<<Label<<endl;
		cout<<"\n--Predecessors--\n";
		Print_pred();
		cout<< "\n--Successors -->\n";
		Print_succ();
		cout<<"\n line no.:";
		Print_line_num();
		cout<<"End of one node\n";
}
	
void 
FGNODE::Print(FILE *file){
		fprintf(file,"\n");
		Print_pred(file);
		fprintf(file, "ID: %i Label: %s\n", ID, Label);
		Print_line_num(file);
		Print_succ(file);
		cout<<"\n";
}

// print a flow graph node to a file
void
FGNODE::Print(ofstream &outfile){
	int len;
	outfile.write((char *) &ID, sizeof(int));	// ID
	len = strlen(Label);
	outfile.write((char *) &len, sizeof(int));	// Label
	outfile.write(Label, len);
	len = line_nums.size();
	outfile.write((char *) &len, sizeof(int));  // number of line no.
	for(int i=0; i<len; i++)		// line no.
		outfile.write((char *) &line_nums[i], sizeof(int));
	len = preds.size();
	outfile.write((char *) &len, sizeof(int));  // number of predecessors
	for(int i=0; i<len; i++)		// predecessors
	  {
	    outfile.write((char *) &preds[i], sizeof(int));
	    in_edge_vec[i]->Print(outfile);
	  }
	len = succs.size();
	outfile.write((char *) &len, sizeof(int));  // number of successors
	for(int i=0; i<len; i++)	
	  {	// successors
		outfile.write((char *) &succs[i], sizeof(int));
		 out_edge_vec[i]->Print(outfile);
	  }
}
// Modified by Liao, to consider the size of edges for the node
unsigned
FGNODE::GetSize()
{
	unsigned int len=0, len_int, edge_size;
	len_int = sizeof(int);
	//FGEDGE ed;
        //edge_size=ed.GetSize();
	edge_size=3*sizeof(IDTYPE)+sizeof(float)+
	  sizeof(FB_EDGE_TYPE)+sizeof(FB_FREQ_TYPE);

	len += len_int;		// ID
	len += len_int;		// Label
	len += strlen(Label);   
	len += len_int;		// line_nums
	len += line_nums.size()*len_int;	
	len += len_int;		// predecessors
	len += preds.size()*(len_int+ edge_size);	
	len += len_int;		// successors
	len += succs.size()*(len_int+ edge_size);
	return len;
}
//
// change the destination node of out going edge from old_node to new_node
//Used for node deletion
// return 1 if success
// 0 otherwise
int
FGNODE::changeOutEdge(IDTYPE old_node, IDTYPE new_node)
{
  int retI=0;
  for (int i=0;i<out_edge_vec.size();i++)
    if (out_edge_vec[i]->destination == old_node) {
      out_edge_vec[i]->destination=new_node;
      retI=1;
      break; // should be only one edge from src --> dst in CFG
    }
  return retI;

}
// return the edge of src-->dst 
FGEDGE * 
FGNODE::getOutEdge(IDTYPE src, IDTYPE dst){
  for(int i=0;i<out_edge_vec.size();i++)
    if ((out_edge_vec[i]->source == src)&&
	(out_edge_vec[i]->destination == dst))
      return out_edge_vec[i];
}
//----------------------------------------------------------
// member functions for FGRAPH
//
// insert a flow graph node into control flow grap, sorted by the node's ID
void 
FGRAPH::Insert_FNode(FGNODE node){
	if(node.Get_ID()>GetSizeofFGnodes()){   
		for(int i=GetSizeofFGnodes(); i<node.Get_ID()-1; i++){
		// insert some temp fnodes to make sure the node is inserted in right position
			FGNODE ftmp;	
			FG_Nodes.push_back(ftmp);
		}
		FG_Nodes.push_back(node);
	}else{	// replace old node
		FG_Nodes[node.Get_ID()-1]=node;
	}
}


void 
FGRAPH::Print(FILE *file){
	fprintf(file, "\n PU_NAME: %s\n", pu_name);
	for(int i=0; i<FG_Nodes.size(); i++){
		FG_Nodes[i].Print(file);
	}
}


void
FGRAPH::Print(ofstream &outfile)
{
	int len;
	len = GetSize();	// get size of the control flow graph
	cout << "size === " << len <<"\n";
	outfile.write((char *) &len, sizeof(int));	
	len = strlen(pu_name);			// pu_name
	outfile.write((char *) &len, sizeof(int));
	outfile.write(pu_name,len);
	len = FG_Nodes.size();
	outfile.write((char *) &len, sizeof(int));      // Nodes number
	cout<<"\nsize of nodes "<<len<<"\n";
	for(int i=0; i<len; i++)		// store each nodes information
		FG_Nodes[i].Print(outfile);
}

unsigned int
FGRAPH::GetSize()   // calculate the length of a flow graph saved in a file
{
	unsigned int len=0, len_int;
	len_int = sizeof(int);
	len += len_int;
	len += strlen(pu_name);    // pu_name 
	len += len_int;		   // Nodes number
	for(int i=0; i<FG_Nodes.size(); i++){
		len += FG_Nodes[i].GetSize();
	}	
	return len;
}		
		
// load control flow graph for a .cfg file, which was saved when Open64 compiled the source code
void 
FGRAPH::Load_CFG(char* filename, char* pu)
{
  FGEDGE * tempEd;
  FG_Nodes.clear();//dangerous operation here!
  ifstream fin(filename,ios::in | ios::binary);
	
  if(fin.fail()){ 
    cout << "Error opening control flow graph file:" << filename << "\n";
    return;
  }

  //cout << "filename == " << filename << "\n";
  //cout << "pu from callgraph: " << pu <<endl;
	fin.seekg(ios::beg);
	while(!fin.eof()){	
		int size;
		fin.read((char *)&size,sizeof(int));
		if(fin.eof())
			break;
		int len;
		fin.read((char *)&len, sizeof(int));
		char* data = new char[len+1];
		fin.read(data, len);	// get pu's name
		data[len]='\0';
		//cout<< " pu name ===== "<<data <<"\n\n";
		char* data_ext = new char[len+3];
		strcpy(data_ext, data);
		strcat(data_ext, "..");
	     
		if(strcmp(data,pu)==0 || (strncmp(data_ext,pu,len+2)==0)){  
		  // find the pu's cfg

			Set_pu_name(pu);
			int id,num,size;
			fin.read((char *)&num, sizeof(int)); 
			// get CFG nodes number
			//cout << "\n size of nodes: " << num <<"\n";
			size=FG_Nodes.size();			

			for(int i=0; i<num; i++){
				fin.read((char *)&id, sizeof(int));  // get ID
				fin.read((char *)&len, sizeof(int));	
				fin.read(data,len);		// get Label;
				data[len]='\0';
				FGNODE fnode(id, data);		
				int size,line,pred,succ;
				fin.read((char *)&size, sizeof(int));  // get line_num size
				for(int j=0; j<size; j++){
					fin.read((char *)&line, sizeof(int)); // get line No.
					fnode.Insert_line_num(line);
				}

				//Step 1. read predecessors' information
				fin.read((char *)&size, sizeof(int)); 
 // get predecessors size
				for(int j=0; j<size; j++){
				  fin.read((char *)&pred, sizeof(int)); 
// get predecessors.
				  // printf("1. predecessors size: %d iter:%d\n",size,j);
				  tempEd=new FGEDGE(fin);
				  fnode.Insert_pred(pred);
				  fnode.Insert_pred(tempEd);
	   //Debug ,Liao
				  //printf("\npred:%d\n",pred);
				  //tempEd->Print();
				}
				//Step 2. read successors' information
				
				fin.read((char *)&size, sizeof(int));  // get successors size
				for(int j=0; j<size; j++){
					fin.read((char *)&succ, sizeof(int)); // get successors
					//printf("2.Successor:size%d iter:%d",size,j);
					tempEd=new FGEDGE(fin);	
					fnode.Insert_succ(succ);
					fnode.Insert_succ(tempEd);
					//printf("\nsucc:%d\n",succ);
					//tempEd->Print();
				}
//				if(id<=i+1)
					Insert_FNode(fnode);
			}
			delete data;
			delete data_ext;
			break;
		   }			
		   delete data;
		   delete data_ext;
		   fin.seekg(size-len-sizeof(int),ios::cur);
	}
	fin.close();
	
}

//Modified by Liao,3/6/2004
//add one more constraint for node deletion, since it is impossible to
// assign new frequency to edges after deletion if the deleted node has
// more than 1 successor.  For Example:
//    a --b ---c
//  e ___| |____d
// 
// if we can delete b
//    a --  ---c
//  e ___| |____d
// how to connect the remaining edges and assign frequency??

void
FGRAPH::Modify_CFG()  /* delete some BBs that don't have matched source code */
{
//	cout<<"begin"<<endl;

   for(int i=0; i<GetSizeofFGnodes(); i++){

     /* set one branch line no. in IF statements the same as the IF statement line no. if the IF statement has two branches but none of them has line_num. */

     if(strcmp(FG_Nodes[i].Get_Kind(), "LOGIF")==0 && FG_Nodes[i].GetSizeofSuccs()==2 && FG_Nodes[i].line_nums.size()>0){
	if(strcmp(FG_Nodes[FG_Nodes[i].succs[0]-1].Get_Kind(), "GOTO")==0 && FG_Nodes[FG_Nodes[i].succs[0]-1].line_nums.size()==0 
	&& strcmp(FG_Nodes[FG_Nodes[i].succs[1]-1].Get_Kind(), "GOTO")==0 && FG_Nodes[FG_Nodes[i].succs[1]-1].line_nums.size()==0){

		FG_Nodes[FG_Nodes[i].succs[0]-1].Insert_line_num(FG_Nodes[i].line_nums[FG_Nodes[i].line_nums.size()-1]);  /* set one branch line no. as the IF statement line no. */
//		cout << "set IF statement" <<endl;
		continue;
	}
     }


/* delete the nodes that has no line no.  */

     if(strcmp(FG_Nodes[i].Get_Kind(), "GOTO")==0 && 
FG_Nodes[i].GetSizeofSuccs() <= 1 
&&FG_Nodes[i].line_nums.size()==0){  // may be a redundant BB
//	cout << "delete a node" << i << endl;
	Delete_Node(FG_Nodes[i]);
	continue;
     } 

/* remove "DOEND" nodes (remove "DOSTART" nodes and replace "DOEND" with "DOSTART")  */

     if(strcmp(FG_Nodes[i].Get_Kind(), "DOSTART")==0  /* && FG_Nodes[i].GetSizeofPreds() == 1 */ 
	&& FG_Nodes[i].GetSizeofSuccs() == 1 && strcmp(FG_Nodes[FG_Nodes[i].succs[0]-1].Get_Kind(), "DOEND")==0){  

	FG_Nodes[FG_Nodes[i].succs[0]-1].Set_Name("DOSTART");
	for(int j=0; j<FG_Nodes[i].line_nums.size(); j++)
		FG_Nodes[FG_Nodes[i].succs[0]-1].Insert_line_num(FG_Nodes[i].line_nums[j]);
	Delete_Node(FG_Nodes[i]);
	continue;
     }



/* remove a statement node before DOSTART (with same line no.) */

     if(strcmp(FG_Nodes[i].Get_Kind(),"GOTO")==0 /*&& FG_Nodes[i].GetSizeofPreds() == 1 */&& FG_Nodes[i].GetSizeofSuccs() == 1 && FG_Nodes[i].line_nums.size()==1 
	&& ((strcmp(FG_Nodes[FG_Nodes[i].succs[0]-1].Get_Kind(), "DOSTART")==0) || (strcmp(FG_Nodes[FG_Nodes[i].succs[0]-1].Get_Kind(), "WHILEEND")==0)) 
	&& FG_Nodes[i].line_nums[0]==FG_Nodes[FG_Nodes[i].succs[0]-1].line_nums[0]){ 

	Delete_Node(FG_Nodes[i]);
	continue;
     }
		

/* if there is a statement before DOSTEP, which has only one line No. and the DOSTEP has no line No. then set the DOSTEP line no. as the statement line no. and remove the statement node */

     if(strcmp(FG_Nodes[i].Get_Kind(),"GOTO")==0 && FG_Nodes[i].GetSizeofPreds() == 1 && FG_Nodes[i].GetSizeofSuccs() == 1 && FG_Nodes[i].line_nums.size()==1 
	&& (strcmp(FG_Nodes[FG_Nodes[i].succs[0]-1].Get_Kind(), "DOSTEP")==0)  && FG_Nodes[FG_Nodes[i].succs[0]-1].line_nums.size()==0){ 

	FG_Nodes[FG_Nodes[i].succs[0]-1].Insert_line_num(FG_Nodes[i].line_nums[0]);
	Delete_Node(FG_Nodes[i]);
	continue;
     }

/* if there is a statement before DOSTEP, which has only one line No. and the DOSTEP has the same line No. with its successor DOSTART. then set the DOSTEP line no. as the statement line no. and remove the statement node 
	DOSTART <----------------	
	  |			 |
	  ---->STMTS---->DOSTEP---
*/

     if(strcmp(FG_Nodes[i].Get_Kind(),"GOTO")==0 && FG_Nodes[i].GetSizeofPreds() == 1 && FG_Nodes[i].GetSizeofSuccs() == 1 && FG_Nodes[i].line_nums.size()==1 
	&& (strcmp(FG_Nodes[FG_Nodes[i].succs[0]-1].Get_Kind(), "DOSTEP")==0)  && FG_Nodes[FG_Nodes[i].succs[0]-1].line_nums.size()==1)
	if ((strcmp(FG_Nodes[FG_Nodes[FG_Nodes[i].succs[0]-1].succs[0]-1].Get_Name(), "DOSTART")==0))  
	if (FG_Nodes[FG_Nodes[FG_Nodes[i].succs[0]-1].succs[0]-1].line_nums.size()>0 
	&& FG_Nodes[FG_Nodes[i].succs[0]-1].line_nums[0] == FG_Nodes[FG_Nodes[FG_Nodes[i].succs[0]-1].succs[0]-1].line_nums[FG_Nodes[FG_Nodes[FG_Nodes[i].succs[0]-1].succs[0]-1].line_nums.size()-1] ){ 

	FG_Nodes[FG_Nodes[i].succs[0]-1].line_nums[0]=FG_Nodes[i].line_nums[0];
	Delete_Node(FG_Nodes[i]);
	continue;
     }


/*   Give REGIONEXIT node line num if it doesn't have one  */
     if(strcmp(FG_Nodes[i].Get_Kind(), "REGIONEXIT")==0  &&  FG_Nodes[i].line_nums.size()==0){
	if(strcmp(FG_Nodes[FG_Nodes[i].preds[0]-1].Get_Name(),"DOSTART")==0){
		for(int j=0; j<FG_Nodes[FG_Nodes[i].preds[0]-1].GetSizeofPreds(); j++){
			int num=FG_Nodes[FG_Nodes[i].preds[0]-1].preds[j]-1;
			if(strcmp(FG_Nodes[num].Get_Kind(),"DOSTEP") == 0 && FG_Nodes[num].line_nums.size()>0){
				FG_Nodes[i].Insert_line_num(FG_Nodes[num].line_nums[FG_Nodes[num].line_nums.size()-1]);
				break;
			}
		}
	}
	continue;
     }



/*   change an EXIT node to STMTS/EXIT if the node has more than one line 	*/
     if(strcmp(FG_Nodes[i].Get_Kind(), "EXIT")==0  &&  FG_Nodes[i].line_nums.size()>1){
	FG_Nodes[i].Set_Name("STMTS/EXIT");
	continue;
     }


/**********************************************************************************/		
/* delete an EXIT node that has neither any line no. predecessors, nor successors */
/**********************************************************************************/

     if(strcmp(FG_Nodes[i].Get_Kind(), "EXIT")==0  &&  FG_Nodes[i].line_nums.size()==0){

	if(FG_Nodes[i].GetSizeofPreds() == 0 && FG_Nodes[i].GetSizeofSuccs() == 0){  // no predecessors and succssors
	
		Delete_Node(FG_Nodes[i]);
		continue;
	}

	int j;
	for(j=0; j< FG_Nodes[i].GetSizeofPreds(); j++){  // search the predecessors
		if(FG_Nodes[i].preds[j]<=GetSizeofFGnodes())	// if valid, then break
			break;
	}
	if(j == FG_Nodes[i].GetSizeofPreds()){	//  all predecessors are invalid		for		for
		for(j =0; j< FG_Nodes[i].GetSizeofSuccs(); j++){  // search the succssors
			if(FG_Nodes[i].succs[j]<=GetSizeofFGnodes())	// if valid, then break
				break;
		}
		if(j == FG_Nodes[i].GetSizeofSuccs()){		// all predecessors and succssors are invalid
			Delete_Node(FG_Nodes[i]);
			continue;
		}
			
	}

     } 

  }

				
}

// Modified by Liao, 3/6/2004
// Added one constraint for node deletion: only node with less than or
// equal to 1 successor.
// otherwise it is impossible to change edge frequency.
void
FGRAPH::Delete_Node(FGNODE &fnode)  // delete a specific Flow graph node *
{
  FGEDGE * tempEd;
  int predId,succId;
/*	int pred;
	pred=fnode.preds[0];
	for(int i=0; i<fnode.succs.size(); i++){
		FG_Nodes[pred-1].Insert_succ(fnode.succs[i]);
		FG_Nodes[fnode.succs[i]-1].Insert_pred(pred);
	}
	for(int i=1; i<fnode.preds.size(); i++){
		FG_Nodes[pred].Insert_succ(fnode.succs[i]);
	} */


//	if(fnode.GetSizeofPreds() == 1 && fnode.GetSizeofSuccs() == 1){
  fnode.Set_Redundant(TRUE);
//		cout << "redundant node: " << fnode.Get_ID() << "\n";

//step1. change successors of the node's predecessors to the node's successors 
  for(int i=0; i<fnode.GetSizeofPreds(); i++){
    if(fnode.preds[i]>GetSizeofFGnodes())
      continue; 
    predId= FG_Nodes[fnode.preds[i]-1].Get_ID();
    
    for(int j=0; j<FG_Nodes[fnode.preds[i]-1].GetSizeofSuccs(); j++){
      if(FG_Nodes[fnode.preds[i]-1].succs[j]==fnode.Get_ID()){
	if(fnode.GetSizeofSuccs()==0){ // in case fnode is the last node
	  FG_Nodes[fnode.preds[i]-1].succs[j]=-1;
	  break;
	}
	FG_Nodes[fnode.preds[i]-1].succs[j]=fnode.succs[0];  // get succs
	for(int k=1; k<fnode.GetSizeofSuccs(); k++)// should be only 1 succ.
	  {
	    succId=fnode.succs[k];
	  FG_Nodes[fnode.preds[i]-1].Insert_succ(fnode.succs[k]);
    	  FG_Nodes[fnode.preds[i]-1].changeOutEdge
	    (fnode.Get_ID(),fnode.succs[k]);  // change edge information ***
	  // ADD new edges to the successor
	  tempEd= new FGEDGE;
	  (*tempEd)= *( FG_Nodes[fnode.preds[i]-1].getOutEdge(predId,succId));
	  FG_Nodes[succId -1].Insert_pred(tempEd);
	  // Do we need to delete the edge from the deleted node?
	  // not likely, since the deleted node has flag of redundant already
	  }
	break;
      }
    }
  }
//Step 2. change predecessors of the node's successors to the node's predessors
		for(int i=0; i<fnode.GetSizeofSuccs(); i++){
			if(fnode.succs[i]>GetSizeofFGnodes())
				continue;
			for(int j=0; j<FG_Nodes[fnode.succs[i]-1].GetSizeofPreds(); j++){
				if(FG_Nodes[fnode.succs[i]-1].preds[j]==fnode.Get_ID()){
					FG_Nodes[fnode.succs[i]-1].preds[j]=fnode.preds[0];  // get preds
					for(int k=1; k<fnode.GetSizeofPreds(); k++)
						FG_Nodes[fnode.succs[i]-1].Insert_pred(fnode.preds[k]);				
					break;
				}
			}
		}
		
//	}
	
}

void
FGRAPH::PrintVCG(){
  PrintGraph myPrint;
/*  GraphNode myNode("one",1);
  myNode.setEdge(6);
  myPrint.addNode(myNode);
  cout << "test one" << endl;
  myPrint.addNode(new GraphNode("two",2));
  myPrint.addNode(new GraphNode("three",3));
  myPrint.addNode(new GraphNode("four",4));
  myPrint.addNode(new GraphNode("five",5));
  myPrint.addNode(new GraphNode("six"));
  myPrint.printNodes("vcg111","toto");
*/
	int len = FG_Nodes.size();

	for(int i=0; i<len; i++){		// store each nodes information

	     	if(FG_Nodes[i].Get_Redundant()){   // don't display redundant BBs
			continue;
		} 
		char title[10];
		sprintf(title,"%d",i+1);
		GraphNode myNode(title, FG_Nodes[i].Get_Name());
		// set the shape and color of the node
		if(strcmp(FG_Nodes[i].Get_Name(),"BRANCH")==0){
			myNode.setShape('r');	// rhomb shape of node
			myNode.setColor('l');
		}else if(strcmp(FG_Nodes[i].Get_Name(),"ENTRY")==0 || strcmp(FG_Nodes[i].Get_Name(),"EXIT")==0 || strcmp(FG_Nodes[i].Get_Name(),"STMTS/EXIT")==0 ){
			myNode.setShape('e');	// ellipse
			myNode.setColor('y');
		}else if(strcmp(FG_Nodes[i].Get_Name(),"DOSTART")==0 || strcmp(FG_Nodes[i].Get_Name(),"DOEND")==0 || strcmp(FG_Nodes[i].Get_Name(),"LOOP")==0){
			myNode.setColor('b');	// blue
                }else if(strcmp(FG_Nodes[i].Get_Name(),"REGIONSTART")==0 || strcmp(FG_Nodes[i].Get_Name(),"REGIONEXIT")==0){
                        myNode.setColor('o');	// orange	
		}

		for(int j=0; j<FG_Nodes[i].succs.size(); j++){
		// avoid the BB node exceeds the total size
	     		if(FG_Nodes[i].succs[j]>len || FG_Nodes[i].succs[j]<0) continue;     
			sprintf(title,"%d",FG_Nodes[i].succs[j]);
			if(strcmp(FG_Nodes[i].Get_Kind(),"LOGIF")==0)	// bent edge
				myNode.setEdge(title,'t');
			else if(strcmp(FG_Nodes[i].Get_Kind(),"DOSTEP")==0)  // backedge
				myNode.setEdge(title,'b');
			else
				myNode.setEdge(title);
		}

		myPrint.addNode(myNode);
	}
  myPrint.printNodes("cfg.vcg","Control Flow Graph");
  system("xvcg cfg.vcg&");
}

