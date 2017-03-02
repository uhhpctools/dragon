/*
** Generated by WorkShop Visual
*/
/*
**LIBS: -lXm -lXt -lX11
*/

#include <stdlib.h>
#include <X11/Xatom.h>
#include <X11/Intrinsic.h>
#include <X11/Shell.h>
#include <X11/Xatom.h>
#include <Xm/RepType.h>
#include <Xm/Protocols.h>
#include <X11/Shell.h>
#include <Xm/Xm.h>
#include <Xm/DialogS.h>
#include <Xm/Form.h>
#include <Xm/Label.h>
#include <Xm/PushB.h>
#include <Xm/Separator.h>
#include <XmGraph/Graph.h>
#include <XmGraph/Arc.h>
#include <Xm/TextF.h>
#include <Xm/Text.h>
#include "ddgui.h"
#include "dragon.h"
#include "graphwindow.h"
#include <string>
#include <vector>
#include <fstream.h>
#include <stdio.h>
#include <string.h>
Widget DDShell = (Widget) NULL;
Widget DDMainForm = (Widget) NULL;
Widget DDseparator1 = (Widget) NULL;
Widget DDlabel = (Widget) NULL;
Widget DDseparator2 = (Widget) NULL;
Widget DDGraphForm = (Widget) NULL;
Widget DDseparator3 = (Widget) NULL;
Widget DDLabelForm = (Widget) NULL;
Widget DDPicForm = (Widget) NULL;
Widget DDClose = (Widget) NULL;
Widget DDgraph = (Widget) NULL;
Widget DDMessageScrolledText = (Widget) NULL;
Widget DDMessageHorScrollBar = (Widget) NULL;
Widget DDMessageVertScrollBar = (Widget) NULL;
Widget DDMessageText = (Widget) NULL;

class DEP
{
   public: 

  DEP() {id=-1; lsc=-1; name=""; pointed=false; visited=false;}
   int id;
   int lsc;
   int line;
  int array;
   string name;
   vector<int> kids;
   bool pointed;
  bool visited;
  string exp; 
   Widget widget;
};

void DDPrintMessage(char Message[]);
extern void ShowDep(Widget, XtPointer,  XtPointer);
vector<DEP> depgraph; 
extern void CloseDDGraph(Widget CloseButtonWidget, XtPointer, XtPointer)
{

  //	XmGraphDestroyAllArcs(DDgraph);
    UnHighlightText();
    DDwindowopen = false;
   // ProcedureLoaded = rootnode;
   XtUnmanageChild(DDShell);
   
}
void ConstructDD(void);
void create_DDShell (Widget parent)
{
	Widget children[6];      /* Children to manage */
	Arg al[64];                    /* Arg List */
	register int ac = 0;           /* Arg Count */
	XtPointer tmp_value;             /* ditto */
	XmString xmstrings[16];    /* temporary storage for XmStrings */

	XtSetArg(al[ac], XmNallowShellResize, False); ac++;
	XtSetArg(al[ac], XmNtitle, "Data Dependence Information"); ac++;
         XtSetArg(al[ac], XmNheight, 400); ac++;
        XtSetArg(al[ac], XmNwidth,600); ac++;
	DDShell = XmCreateDialogShell ( parent, (char *) "DDShell", al, ac );

         Atom WmDeleteWindow = XmInternAtom(XtDisplay(DDShell),"WM_DELETE_WINDOW",False);
        XmAddWMProtocolCallback(DDShell, WmDeleteWindow, CloseDDGraph, NULL); 

	ac = 0;
	xmstrings[0] = XmStringCreateLtoR ( "Data Dependence Information", (XmStringCharSet)XmFONTLIST_DEFAULT_TAG );
	XtSetArg(al[ac], XmNdialogTitle, xmstrings[0]); ac++;
	XtSetArg(al[ac], XmNautoUnmanage, FALSE); ac++;
         XtSetArg(al[ac], XmNheight, 400); ac++;
        XtSetArg(al[ac], XmNwidth,600); ac++;
	DDMainForm = XmCreateForm ( DDShell, (char *) "DDMainForm", al, ac );
	ac = 0;
	XmStringFree ( xmstrings [ 0 ] );
	XtSetArg(al[ac], XmNseparatorType, XmSINGLE_DASHED_LINE); ac++;
	DDseparator1 = XmCreateSeparator ( DDMainForm, (char *) "DDseparator1", al, ac );
	ac = 0;
	xmstrings[0] = XmStringCreateLtoR ( "Array Data Dependence for: ", (XmStringCharSet)XmFONTLIST_DEFAULT_TAG );
	XtSetArg(al[ac], XmNlabelString, xmstrings[0]); ac++;
	DDlabel = XmCreateLabel ( DDMainForm, (char *) "DDlabel", al, ac );
	ac = 0;
	XmStringFree ( xmstrings [ 0 ] );
	DDseparator2 = XmCreateSeparator ( DDMainForm, (char *) "DDseparator2", al, ac );




        ac = 0;

        
        XtSetArg(al[ac], XmNheight, 100); ac++;
        XtSetArg(al[ac], XmNwidth,600); ac++;
	DDLabelForm = XmCreateForm ( DDMainForm, (char *) "DDLabelForm", al, ac );

        ac =0;
	//  XtSetArg(al[ac], XmNheight, 300); ac++;
        //XtSetArg(al[ac], XmNwidth,600); ac++;
	DDGraphForm = XmCreateForm ( DDMainForm, (char *) "DDGraphForm", al, ac );
        ac = 0; 
         XtSetArg(al[ac], XmNorientation,XmHORIZONTAL); ac++;
        XtSetArg(al[ac], XmNautoLayoutMode, XmNEVER); ac++;
        XtSetArg(al[ac], XmNtwinsVisible, TRUE); ac++;
        XtSetArg(al[ac], XmNeditable, True); ac++;
        DDgraph = XmCreateScrolledGraph(DDGraphForm,"DDgraph", al, ac);
  

          XtAddCallback (DDgraph,
    		     XmNnodeMovedCallback,
    		     (XtCallbackProc)MoveHandling,   NULL);
        
        XtAddCallback (DDgraph,
    		     XmNarcMovedCallback,
    		     (XtCallbackProc)MoveHandling,    NULL);
         
         /*  avoid creation of new arcs or nodes */
        XtAddCallback (DDgraph,
    		     XmNnewNodeCallback,
    		     (XtCallbackProc)AvoidNew,     NULL);
      
        XtAddCallback (DDgraph,
    		     XmNnewArcCallback,
    		     (XtCallbackProc)AvoidNew,      NULL);
   
       
        XtAddCallback (DDgraph,
    		     XmNselectNodeCallback,
    		     (XtCallbackProc)ShowDep,  NULL);
        ac = 0;

        XtSetArg(al[ac], XmNtopAttachment, XmATTACH_FORM);ac++;
        XtSetArg(al[ac], XmNbottomAttachment, XmATTACH_FORM);ac++;
        XtSetArg(al[ac], XmNleftAttachment, XmATTACH_FORM);ac++;
        XtSetArg(al[ac], XmNrightAttachment, XmATTACH_FORM);ac++;
        XtSetValues (XtParent(XtParent(DDgraph)), al,ac);
        
     


        ac =0;
 	DDseparator3 = XmCreateSeparator ( DDMainForm, (char *) "DDseparator3", al, ac );

	//  ac = 0;

        
        //XtSetArg(al[ac], XmNheight, 100); ac++;
        //XtSetArg(al[ac], XmNwidth,600); ac++;
	//	DDLabelForm = XmCreateForm ( DDMainForm, (char *) "DDLabelForm", al, ac );
        ac = 0;
	DDPicForm = XmCreateForm ( DDLabelForm, (char *) "DDPicForm", al, ac );
	xmstrings[0] = XmStringCreateLtoR ( "Close", (XmStringCharSet)XmFONTLIST_DEFAULT_TAG );
	XtSetArg(al[ac], XmNlabelString, xmstrings[0]); ac++;
	DDClose = XmCreatePushButton ( DDLabelForm, (char *) "DDClose", al, ac );

        XtAddCallback( DDClose, XmNactivateCallback, CloseDDGraph, (XtPointer) 0 );


        ac = 0;

        DDMessageText = XmCreateScrolledText (DDPicForm, (char *) "MessageText", al, ac );
	ac = 0;


	DDMessageScrolledText = XtParent ( DDMessageText );
	XtSetArg(al[ac], XmNhorizontalScrollBar, &DDMessageHorScrollBar); ac++;
	XtSetArg(al[ac], XmNverticalScrollBar, &DDMessageVertScrollBar); ac++;
	XtGetValues(DDMessageScrolledText, al, ac );
 

	ac = 0;
	XmStringFree ( xmstrings [ 0 ] );


	XtSetArg(al[ac], XmNtopAttachment, XmATTACH_FORM); ac++;
	XtSetArg(al[ac], XmNleftAttachment, XmATTACH_FORM); ac++;
	XtSetArg(al[ac], XmNrightAttachment, XmATTACH_FORM); ac++;
	XtSetValues ( DDseparator1,al, ac );
	ac = 0;

	XtSetArg(al[ac], XmNtopAttachment, XmATTACH_WIDGET); ac++;
	XtSetArg(al[ac], XmNtopWidget, DDseparator1); ac++;
	XtSetArg(al[ac], XmNleftAttachment, XmATTACH_FORM); ac++;
	XtSetArg(al[ac], XmNrightAttachment, XmATTACH_FORM); ac++;
	XtSetValues ( DDlabel,al, ac );
	ac = 0;

	XtSetArg(al[ac], XmNtopAttachment, XmATTACH_WIDGET); ac++;
	XtSetArg(al[ac], XmNtopWidget, DDlabel); ac++;
	XtSetArg(al[ac], XmNleftAttachment, XmATTACH_FORM); ac++;
	XtSetArg(al[ac], XmNrightAttachment, XmATTACH_FORM); ac++;
	XtSetValues ( DDseparator2,al, ac );
	ac = 0;

	XtSetArg(al[ac], XmNtopAttachment, XmATTACH_NONE); ac++;
	//	XtSetArg(al[ac], XmNtopWidget, DDseparator3); ac++;
	XtSetArg(al[ac], XmNbottomAttachment, XmATTACH_FORM); ac++;
	XtSetArg(al[ac], XmNleftAttachment, XmATTACH_FORM); ac++;
	XtSetArg(al[ac], XmNrightAttachment, XmATTACH_FORM); ac++;
	XtSetValues ( DDLabelForm,al, ac );
       
	ac = 0;

	//XtSetArg(al[ac], XmNtopAttachment, XmATTACH_WIDGET); ac++;
	//XtSetArg(al[ac], XmNtopWidget, DDGraphForm); ac++;
	//XtSetArg(al[ac], XmNbottomAttachment, XmATTACH_NONE); ac++;
        XtSetArg(al[ac], XmNtopAttachment, XmATTACH_NONE); ac++;
        XtSetArg(al[ac], XmNbottomAttachment, XmATTACH_WIDGET); ac++;
        XtSetArg(al[ac], XmNbottomWidget, DDLabelForm); ac++;
	XtSetArg(al[ac], XmNleftAttachment, XmATTACH_FORM); ac++;
	XtSetArg(al[ac], XmNrightAttachment, XmATTACH_FORM); ac++;
	XtSetValues ( DDseparator3,al, ac );
	ac = 0;

	


        XtSetArg(al[ac], XmNtopAttachment, XmATTACH_WIDGET); ac++;
	XtSetArg(al[ac], XmNtopWidget, DDseparator2); ac++;
	//	XtSetArg(al[ac], XmNbottomAttachment, XmATTACH_NONE); ac++;
        XtSetArg(al[ac], XmNbottomWidget, DDseparator3); ac++;
        XtSetArg(al[ac], XmNbottomAttachment, XmATTACH_WIDGET); ac++;
	XtSetArg(al[ac], XmNleftAttachment, XmATTACH_FORM); ac++;
	XtSetArg(al[ac], XmNrightAttachment, XmATTACH_FORM); ac++;
	XtSetValues ( DDGraphForm,al, ac );
       
	ac = 0;

	XtSetArg(al[ac], XmNtopAttachment, XmATTACH_NONE); ac++;
	XtSetArg(al[ac], XmNbottomAttachment, XmATTACH_FORM); ac++;
	XtSetArg(al[ac], XmNleftAttachment, XmATTACH_NONE); ac++;
	XtSetArg(al[ac], XmNrightAttachment, XmATTACH_FORM); ac++;
	XtSetValues ( DDClose,al, ac );
        
        ac=0;
	XtSetArg(al[ac], XmNtopAttachment, XmATTACH_FORM); ac++;
	XtSetArg(al[ac], XmNbottomAttachment, XmATTACH_FORM); ac++;
	XtSetArg(al[ac], XmNleftAttachment, XmATTACH_FORM); ac++;
        XtSetArg(al[ac], XmNrightWidget, DDClose); ac++;
        XtSetArg(al[ac], XmNrightAttachment, XmATTACH_WIDGET); ac++;
	XtSetValues ( DDPicForm,al, ac );

        ac=0;
	XtSetArg(al[ac], XmNtopAttachment, XmATTACH_FORM); ac++;
	XtSetArg(al[ac], XmNbottomAttachment, XmATTACH_FORM); ac++;
	XtSetArg(al[ac], XmNleftAttachment, XmATTACH_FORM); ac++;
        XtSetArg(al[ac], XmNrightAttachment, XmATTACH_FORM); ac++;  
        XtSetValues ( DDMessageScrolledText,al, ac );

	ac = 0;

	if ((children[ac] = DDPicForm) != (Widget) 0) { ac++; }
	if ((children[ac] = DDClose) != (Widget) 0) { ac++; }
	if (ac > 0) { XtManageChildren(children, ac); }
	ac = 0;
	if ((children[ac] = DDseparator1) != (Widget) 0) { ac++; }
	if ((children[ac] = DDlabel) != (Widget) 0) { ac++; }
	if ((children[ac] = DDseparator2) != (Widget) 0) { ac++; }
	if ((children[ac] = DDGraphForm) != (Widget) 0) { ac++; }
	if ((children[ac] = DDseparator3) != (Widget) 0) { ac++; }
	if ((children[ac] = DDLabelForm) != (Widget) 0) { ac++; }
       
	if (ac > 0) { XtManageChildren(children, ac); }
	ac = 0;
	DDwindowopen=true;
          if (DDgraph != (Widget) 0) { XtManageChild(DDgraph);} 
        XtManageChild(DDMainForm);
           ConstructDD();
          XmGraphLayout(DDgraph);
          XmUpdateDisplay(DDgraph);
       

      
	XtManageChild(DDShell);
        

}







void ConstructDD(void)
{
       
        Arg al[10];                    /* Arg List */
	register int ac = 0;
        bool debug = false;
        char *PU_Name, *FileName;
	int len=0, n;
        int totaltrue=0, totalanti=0, totaloutput=0, totalcall=0;
      
        bool graphfound= false;
        int graphsize = 0;       
        DEP dep;
        depgraph.clear(); 
          XmTextSetString(DDMessageText,"");
	len=strlen(UICallgraphRoot[ProcedureLoaded].GetName());
	PU_Name = new char[len+1];
	strcpy(PU_Name, UICallgraphRoot[ProcedureLoaded].GetName());
	len = strlen(UICallgraphRoot[ProcedureLoaded].GetFilename());
	FileName = new char[len+5];
 	strcpy(FileName,UICallgraphRoot[ProcedureLoaded].GetFilename());
	for(n=len-1; n>0; n--)
		if(FileName[n]!='.')
			FileName[n]='\0';
		else
			break;
	if(n>0)
		strcat(FileName,"dep");
	else{
		strcpy(FileName,UICallgraphRoot[ProcedureLoaded].GetFilename());
		strcat(FileName,".dep");
	}   

      
        fstream infile(FileName, ios::in | ios::bin);
       if (!infile) 
	{
            string message = "No dependence information file: "; 
            message = message + FileName;
            DDPrintMessage((char *)message.c_str());
          delete [] FileName;
        delete [] PU_Name;
         return; //add message.
         }
        char C;
      while((C=infile.peek()) != EOF) {
        infile.read((char *) &len, sizeof(int));
        
        char *retreived_name = new char[len];
        infile.read(retreived_name,len);
        if (debug) cout << "\nretreived function="<<retreived_name;
	// delete [] retreived_name;
        int total_nodes = 0;
        infile.read((char *) &total_nodes, sizeof(int));
        if (debug) cout << "\n total nodes="<<total_nodes;
       
         
         for (int i=0; i<total_nodes; i++)
	 {
            vector<int> kids;
	  int nodeid;
          infile.read((char *) &nodeid, sizeof(int));
	  if (debug) cout << "\nNode ID="<< nodeid;
     
          int line;
          infile.read((char *) &line, sizeof(int));
	  if (debug) cout << "\nline="<< line;

          int array; 
           infile.read((char *) &array, sizeof(int));
	  if (debug) cout << "\narray="<< array;

          int lsc;
          infile.read((char *) &lsc, sizeof(int));
	  if (debug) cout << "\nlsc="<< lsc;


          int len_name;
           infile.read((char *) &len_name, sizeof(int));
	  if (debug) cout << "\nlen_name="<< len_name;

          char *array_name = new char[len_name];
          infile.read(array_name, len_name);
          if (debug) cout << "\narray_name="<<array_name;
	  string name_string=array_name;
          delete [] array_name; 

          string exp_string="";
          int array_expressions=0;
          infile.read((char *) &array_expressions, sizeof(int));
	  if (debug) cout << "\n #array expressions="<< array_expressions; 
	 
          for (int exp=0; exp <array_expressions; exp++)
	    {  
               
	      int array_exp_len;
              infile.read((char *) &array_exp_len, sizeof (int));
              
              char *array_exp_str = new char [array_exp_len];  
               infile.read((char *) array_exp_str, array_exp_len);
               
	       exp_string = exp_string + array_exp_str;
               if (debug) cout << "\narray expressions="<< exp_string;
               delete [] array_exp_str;

            }

           int num_edges=0;
          infile.read((char *) &num_edges, sizeof(int));
	  if (debug) cout << "\n #num edges="<< num_edges;


          for (int edg=0; edg < num_edges; edg++)
	  {
               
             int edge=0;
             infile.read((char *) &edge, sizeof(int));
             kids.push_back(edge);
	     if (debug) cout <<"\nedge="<< edge;

	  }
	  //  cout << retreived_name<<":"<<PU_Name;   
	  if(strcmp(retreived_name, PU_Name)==0)
	    { 

	      if (debug)
                        cout << "\nAllocating element for: "<<PU_Name<<endl; 
             dep.id = nodeid;
             dep.name = name_string;
             dep.exp = exp_string;
             dep.kids = kids;
             dep.line = line;
             dep.array = array;
             dep.lsc = lsc; 
               
             depgraph.push_back(dep);
             graphfound=true;
             graphsize=total_nodes;
            } 

	 } // end of individual graph or nodes.

        delete [] retreived_name;
      }// end of file 

      if (!graphfound)
	{
	  string message = "There are no dependences in:  "; 
            message = message + PU_Name;
            DDPrintMessage((char *)message.c_str());
             delete [] FileName;
             delete [] PU_Name;
         return;

       
        }

      bool graph_empty = false;
      int nodes_in_graph =0;
       for(int i=0; i<graphsize; i++)
       {
           if(!depgraph[i].visited)             
              depgraph[i].visited = true; 

            int num_succs = depgraph[i].kids.size();

          if (num_succs > 0)
	    {
              depgraph[i].pointed = true;
              nodes_in_graph++;
              for (int j=0; j<num_succs; j++)
                { 
                  // find correct node with succs
                  int succ = depgraph[i].kids[j];
                  for (int k=0; k<depgraph.size(); k++) { if (depgraph[k].id == succ) { succ = k; break;}}
                   depgraph[succ].pointed = true;
                    nodes_in_graph++;
                   if(succ==i && depgraph[i].lsc==0 && depgraph[succ].lsc==0)
                     {
                       nodes_in_graph = nodes_in_graph-2;
                       depgraph[i].pointed = false;
                       depgraph[succ].pointed = false;
                     }
                  
                }
            }

       } 

       if (nodes_in_graph==0)
	 {
            string message = "There are no dependences in:  "; 
            message = message + PU_Name;
            DDPrintMessage((char *)message.c_str());
             delete [] FileName;
             delete [] PU_Name;
              return;

         } 

       for(int i=0; i < graphsize; i++) depgraph[i].visited = false; 

      for(int i=0; i<graphsize; i++) 
	{
          if(!depgraph[i].visited && depgraph[i].pointed) {
            depgraph[i].widget = CreateAnEdge(DDgraph,(char *)depgraph[i].name.c_str(),NULL,"test2");
            depgraph[i].visited = true;     
          }
          int num_succs = depgraph[i].kids.size();

          if (num_succs > 0)
	    {
              depgraph[i].pointed = true;
              for (int j=0; j<num_succs; j++)
                { 
                  // find correct node with succs
                  int succ = depgraph[i].kids[j];
                  for (int k=0; k<depgraph.size(); k++) { if (depgraph[k].id == succ) { succ = k; break;}}
                   depgraph[succ].pointed = true; 
                  if(!depgraph[succ].visited && depgraph[succ].pointed)
		    {
		      // depgraph[succ].widget = CreateAnEdgeDD(DDgraph, (char *)depgraph[succ].name.c_str(), depgraph[i].widget,NULL);
                      depgraph[succ].visited = true;
                     
                       if(depgraph[i].lsc==0 && depgraph[succ].lsc==1)
			 {
			   char label[]="anti";
                           char color[]="blue";
                             totalanti++;
                              depgraph[succ].widget = CreateAnEdgeDD(DDgraph, (char *)depgraph[succ].name.c_str(), depgraph[i].widget,label,color);
                         }
                          if(depgraph[i].lsc==1 && depgraph[succ].lsc==0)
			 {
                             char label[]="true";
                           char color[]="red";
                            totaltrue++;
                             depgraph[succ].widget = CreateAnEdgeDD(DDgraph, (char *)depgraph[succ].name.c_str(), depgraph[i].widget,label,color);
                         } 
                          if(depgraph[i].lsc==1 && depgraph[succ].lsc==1)
			 {
                            char label[]="out";
                           char color[]="green";
                           totaloutput++;
                            depgraph[succ].widget = CreateAnEdgeDD(DDgraph, (char *)depgraph[succ].name.c_str(), depgraph[i].widget,label,color);
                         }
                          if(depgraph[i].lsc==2 || depgraph[succ].lsc==2)
			{ 
                           char label[]="call";
                           char color[]="orange";
                          totalcall++;
                            depgraph[succ].widget = CreateAnEdgeDD(DDgraph, (char *)depgraph[succ].name.c_str(), depgraph[i].widget,label,color);
                        }

                    }
		  else if (depgraph[succ].pointed)
		     {
                        Arg al[5];
                        int ac = 0; 
                        Widget NewArc;  
                        XtSetArg(al[ac],XmNmapLabel,"True"); ac++;
                        XtSetArg(al[ac],XmNarcDirection, XmDIRECTED); ac++;
                        
                        //Widget NewArc = XmCreateAttachedArc(DDgraph,"NewArc",depgraph[i].widget, 
                        //                                    depgraph[succ].widget,
                        //                                    al,ac); 
                         

                          if(depgraph[i].lsc==0 && depgraph[succ].lsc==1)
			    {  
			      totalanti++;
                               NewArc = XmCreateAttachedArc(DDgraph,"anti",depgraph[i].widget, 
                                                            depgraph[succ].widget,
                                                            al,ac);
                                XtManageChild(NewArc); 
                              ColorWidget(NewArc,"blue",False);
                            }
                          if(depgraph[i].lsc==1 && depgraph[succ].lsc==0)
			    {
                              totaltrue++;
                               NewArc = XmCreateAttachedArc(DDgraph,"true",depgraph[i].widget, 
                                                            depgraph[succ].widget,
                                                            al,ac); 
                               XtManageChild(NewArc); 
                              ColorWidget(NewArc,"red",False);
                            }  
                          if(depgraph[i].lsc==1 && depgraph[succ].lsc==1)
                            {
			      totaloutput++;
                                NewArc = XmCreateAttachedArc(DDgraph,"out",depgraph[i].widget, 
                                                            depgraph[succ].widget,
                                                            al,ac); 
                                XtManageChild(NewArc); 
                              ColorWidget(NewArc,"green",False);
                            } 
                          if(depgraph[i].lsc==2 || depgraph[succ].lsc==2)
			    { 
			      totalcall++;
                                NewArc = XmCreateAttachedArc(DDgraph,"call",depgraph[i].widget, 
                                                            depgraph[succ].widget,
                                                            al,ac); 
                                XtManageChild(NewArc);  
                              ColorWidget(NewArc,"orange",False);
                            }   


                     }

                  
       
                } //end for succss
            }
        }

      // for(int i=0; i<graphsize; i++)
      //	   if (!depgraph[i].pointed) XtUnmanageChild(depgraph[i].widget); 
             

        if(debug) 
        cout << "filename from callgraph" << FileName; 
        cout.flush();


        ac = 0;
        char templabel[100];
        strcpy(templabel,"Array Dependence Information for : ");
        strcat(templabel,PU_Name);
        XmString tempstring = XmStringCreateLtoR (templabel, (XmStringCharSet)XmFONTLIST_DEFAULT_TAG );
	XtSetArg(al[ac], XmNlabelString, tempstring); ac++;
	XtSetValues(DDlabel, al, ac );
        XmStringFree ( tempstring);

	char ttrue[10], tanti[10], tout[10], tcall[10];
         sprintf(ttrue, "%d", totaltrue);
          sprintf(tanti, "%d", totalanti);
          sprintf(tout, "%d", totaloutput);
           sprintf(tcall, "%d", totalcall);
 
         string depmessage=""; 
      
         depmessage = "Data dependence summary for: ";
         depmessage += PU_Name;
         depmessage += "\n";
	 depmessage += ttrue;
	 depmessage += " true dependences.\n"; 
	 depmessage += tanti; 
	 depmessage += " anti dependences.\n"; 
	 depmessage += tout; 
         depmessage +=" output dependences.\n";
         depmessage += tcall; 
         depmessage +=" procedure dependences.";  
         DDPrintMessage((char *)depmessage.c_str());
          XmTextShowPosition(DDMessageText, 0);           



        delete [] FileName;
        delete [] PU_Name;





}

extern void ShowDep(Widget w, XtPointer client,  XtPointer Calldata)
{
     XmGraphCallbackStruct *calld = (XmGraphCallbackStruct *) Calldata;

  if (calld->reason == XmCR_SELECT_NODE && XtIsSubclass(calld->widget,xmPushButtonWidgetClass))
  {
    for (int i = 0; i<depgraph.size(); i++) 
    {  
      if(calld->widget == depgraph[i].widget)
	{
	  string message=""; 
           UnHighlightText();
	   HighlightLine(depgraph[i].line-1);
            
           switch (depgraph[i].lsc)
	     {
	       case 0:
	       message = "Reading from "; 
	       break;
	       case 1:
		 message = "Writing to ";
	       break;
	       case 2:
		 message = "Calling procedure ";
	       break;  
 
             } 

            switch (depgraph[i].array)
	     {
	       case 0:
	       message = message + "scalar "; 
	       break;
	       case 1:
		 message = message + "array ";
	       break;
	        
             } 

	    message = message + depgraph[i].name +" (line = ";
            char linenum[10];
            sprintf(linenum, "%d", depgraph[i].line-1); 
            message = message + linenum;
            message = message + ").";

	   DDPrintMessage((char *)message.c_str());




          
 
        } 
    }
  }
}


void UpdateDDGraph()
{

       Arg al[5];  
      int ac;
    
    
      for(int i=0; i<depgraph.size(); i++)
	{
           if(depgraph[i].pointed)
	       XtUnmanageChild(depgraph[i].widget);
        }  
		
        XtUnmanageChild(DDgraph);
	
        ac = 0;   
        XtSetArg(al[ac], XmNorientation,XmHORIZONTAL); ac++;
        XtSetArg(al[ac], XmNautoLayoutMode, XmALWAYS); ac++;
        XtSetArg(al[ac], XmNtwinsVisible, TRUE); ac++;
        XtSetArg(al[ac], XmNeditable, True); ac++;
        DDgraph = XmCreateScrolledGraph(DDGraphForm,"DDgraph", al, ac);

        /*  allow moving of nodes */
        XtAddCallback (DDgraph ,
    		     XmNnodeMovedCallback,
    		     (XtCallbackProc)MoveHandling,   NULL);
         /* avoid moving of arcs */
        XtAddCallback ( DDgraph ,
    		     XmNarcMovedCallback,
    		     (XtCallbackProc)MoveHandling,    NULL);
    
         /*  avoid creation of new arcs or nodes */
        XtAddCallback ( DDgraph        ,
    		     XmNnewNodeCallback,
    		     (XtCallbackProc)AvoidNew,     NULL);
      
        XtAddCallback (  DDgraph,
    		     XmNnewArcCallback,
    		     (XtCallbackProc)AvoidNew,      NULL);
   
       
        XtAddCallback (DDgraph,
    		     XmNselectNodeCallback,
    		     (XtCallbackProc)ShowDep,  NULL);
        /*
        XtAddCallback(Flowgraph,
    		    XmNselectArcCallback,
    		    (XtCallbackProc)CallSiteInfo,NULL);
        */




	ac = 0;

        XtSetArg(al[ac], XmNtopAttachment, XmATTACH_FORM);ac++;
        XtSetArg(al[ac], XmNbottomAttachment, XmATTACH_FORM);ac++;
        XtSetArg(al[ac], XmNleftAttachment, XmATTACH_FORM);ac++;
        XtSetArg(al[ac], XmNrightAttachment, XmATTACH_FORM);ac++;
        XtSetValues (XtParent(XtParent(DDgraph)), al,ac);
       

        //XmGraphUnmap(DDgraph);

        ConstructDD();

	// XmGraphMap(DDgraph);
       
         XtManageChild (DDgraph);
           XmGraphUnmap(DDgraph);
           XmGraphMap(DDgraph);


}
void DDPrintMessage(char Message[])
{

  char *TextToInsert = new char[strlen(Message)+2];
  strcpy(TextToInsert,Message);
  TextToInsert[strlen(Message)] = '\n';
  TextToInsert[strlen(Message)+1] = '\0';
  
  XmTextPosition CurPos;
  CurPos = XmTextGetInsertionPosition(DDMessageText);
  XmTextInsert(DDMessageText, CurPos, TextToInsert);

  CurPos = CurPos + strlen(TextToInsert);
  XmTextShowPosition(DDMessageText, CurPos);
  XmTextSetInsertionPosition(DDMessageText, CurPos);
  
  delete [] TextToInsert;


}
