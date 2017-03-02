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
#include <X11/Intrinsic.h>
#include <Xm/RepType.h>
#include <Xm/Protocols.h>

#include <Xm/Xm.h>
#include <Xm/DialogS.h>
#include <Xm/Form.h>
#include <Xm/Label.h>
#include <Xm/List.h>
#include <Xm/PushB.h>
#include <Xm/ScrollBar.h>
#include <Xm/Separator.h>

#include <string>
#include <stdio.h>
#include <mysql.h>
#include <vector>
#include <iostream.h>

#include "dragon.h"
#include "graphwindow.h"
#include "regions.h"


MYSQL_RES *result;
MYSQL_ROW row;
MYSQL *connection,mysql;
char selected_proc[100];

Widget InterRegionShell = (Widget) NULL;
Widget InterRegionForm = (Widget) NULL;
Widget list1 = (Widget)NULL;
Widget list2 = (Widget)NULL;
Widget list3 = (Widget)NULL;
Widget list4 = (Widget)NULL;
Widget list5 = (Widget)NULL;

void CreateRegions(void); 


void CloseRegions(Widget CloseButtonWidget, XtPointer, XtPointer)
{

  //	XmGraphDestroyAllArcs(Callgraph);
   UnHighlightText();
   Regionswindowopen = False;
   // ProcedureLoaded = rootnode;
   XtUnmanageChild(InterRegionShell);
   mysql_close(connection);
}

void ProcedureCallback(Widget CloseButtonWidget, XtPointer, XtPointer);
void ArrayNameCallback(Widget CloseButtonWidget, XtPointer, XtPointer);


void create_InterRegionShell (Widget parent)
{
	Widget children[14];      /* Children to manage */
	Arg al[64];                    /* Arg List */
	register int ac = 0;           /* Arg Count */
	XtPointer tmp_value;             /* ditto */
	Widget separator1 = (Widget)NULL;
	Widget label1 = (Widget)NULL;
	Widget separator2 = (Widget)NULL;
	Widget label2 = (Widget)NULL;
	Widget label3 = (Widget)NULL;
	Widget label4 = (Widget)NULL;
	Widget label5 = (Widget)NULL;
	Widget label6 = (Widget)NULL;
	Widget scrolledList1 = (Widget)NULL;
	Widget scrolledList2 = (Widget)NULL;
	Widget scrolledList3 = (Widget)NULL;
	Widget scrolledList4 = (Widget)NULL;
	Widget button2 = (Widget)NULL;

       
        XtSetArg(al[ac], XmNallowShellResize, TRUE); ac++;
	XtSetArg(al[ac], XmNtitle, "Interprocedural Array Regions Analysis"); ac++;
        XtSetArg(al[ac], XmNheight, 400); ac++;
        XtSetArg(al[ac], XmNwidth,700); ac++;
       
	InterRegionShell = XtCreatePopupShell ( (char *) "InterRegionShell", topLevelShellWidgetClass, parent, al, ac );

        Atom WmDeleteWindow = XmInternAtom(XtDisplay(InterRegionShell),"WM_DELETE_WINDOW",False);
        XmAddWMProtocolCallback(InterRegionShell, WmDeleteWindow, CloseRegions, NULL); 



	ac = 0;
	XtSetArg(al[ac], XmNautoUnmanage, FALSE); ac++;
	InterRegionForm = XmCreateForm ( InterRegionShell, (char *) "InterRegionForm", al, ac );
	ac = 0;
	XtSetArg(al[ac], XmNseparatorType, XmSINGLE_DASHED_LINE); ac++;
	separator1 = XmCreateSeparator ( InterRegionForm, (char *) "separator1", al, ac );
	ac = 0;
	label1 = XmCreateLabel ( InterRegionForm, (char *) "label1", al, ac );
	separator2 = XmCreateSeparator ( InterRegionForm, (char *) "separator2", al, ac );
	label2 = XmCreateLabel ( InterRegionForm, (char *) "label2", al, ac );
	label3 = XmCreateLabel ( InterRegionForm, (char *) "label3", al, ac );
	label4 = XmCreateLabel ( InterRegionForm, (char *) "label4", al, ac );
	label5 = XmCreateLabel ( InterRegionForm, (char *) "label5", al, ac );
	label6 = XmCreateLabel ( InterRegionForm, (char *) "label6", al, ac );
	list1 = XmCreateScrolledList ( InterRegionForm, (char *) "list1", al, ac );
	scrolledList1 = XtParent ( list1 );

	XtSetArg(al[ac], XmNlistSizePolicy, XmRESIZE_IF_POSSIBLE); ac++;
	list2 = XmCreateList ( InterRegionForm, (char *) "list2", al, ac );
	ac = 0;
	list3 = XmCreateScrolledList ( InterRegionForm, (char *) "list3", al, ac );
	scrolledList2 = XtParent ( list3 );

	list4 = XmCreateScrolledList ( InterRegionForm, (char *) "list4", al, ac );
	scrolledList3 = XtParent ( list4 );

	list5 = XmCreateScrolledList ( InterRegionForm, (char *) "list5", al, ac );
	scrolledList4 = XtParent ( list5 );

	button2 = XmCreatePushButton ( InterRegionForm, (char *) "button2", al, ac );
        XtAddCallback( button2, XmNactivateCallback, CloseRegions, (XtPointer) 0 ); 

	XtSetArg(al[ac], XmNtopAttachment, XmATTACH_FORM); ac++;
	XtSetArg(al[ac], XmNtopOffset, 3); ac++;
	XtSetArg(al[ac], XmNbottomAttachment, XmATTACH_NONE); ac++;
	XtSetArg(al[ac], XmNleftAttachment, XmATTACH_FORM); ac++;
	XtSetArg(al[ac], XmNrightAttachment, XmATTACH_FORM); ac++;
	XtSetValues ( separator1,al, ac );
	ac = 0;

	XtSetArg(al[ac], XmNtopAttachment, XmATTACH_WIDGET); ac++;
	XtSetArg(al[ac], XmNtopWidget, separator1); ac++;
	XtSetArg(al[ac], XmNleftAttachment, XmATTACH_FORM); ac++;
	XtSetArg(al[ac], XmNrightAttachment, XmATTACH_FORM); ac++;
	XtSetValues ( label1,al, ac );
	ac = 0;

	XtSetArg(al[ac], XmNtopAttachment, XmATTACH_WIDGET); ac++;
	XtSetArg(al[ac], XmNtopWidget, label1); ac++;
	XtSetArg(al[ac], XmNleftAttachment, XmATTACH_FORM); ac++;
	XtSetArg(al[ac], XmNrightAttachment, XmATTACH_FORM); ac++;
	XtSetValues ( separator2,al, ac );
	ac = 0;

	XtSetArg(al[ac], XmNtopAttachment, XmATTACH_WIDGET); ac++;
	XtSetArg(al[ac], XmNtopWidget, separator2); ac++;
	XtSetArg(al[ac], XmNleftAttachment, XmATTACH_FORM); ac++;
	XtSetValues ( label2,al, ac );
	ac = 0;

	XtSetArg(al[ac], XmNtopAttachment, XmATTACH_WIDGET); ac++;
	XtSetArg(al[ac], XmNtopWidget, separator2); ac++;
	XtSetArg(al[ac], XmNleftAttachment, XmATTACH_WIDGET); ac++;
	XtSetArg(al[ac], XmNleftWidget, label2); ac++;
	XtSetValues ( label3,al, ac );
	ac = 0;

	XtSetArg(al[ac], XmNtopAttachment, XmATTACH_WIDGET); ac++;
	XtSetArg(al[ac], XmNtopWidget, separator2); ac++;
	XtSetArg(al[ac], XmNleftAttachment, XmATTACH_WIDGET); ac++;
	XtSetArg(al[ac], XmNleftWidget, label3); ac++;
	XtSetValues ( label4,al, ac );
	ac = 0;

	XtSetArg(al[ac], XmNtopAttachment, XmATTACH_WIDGET); ac++;
	XtSetArg(al[ac], XmNtopWidget, separator2); ac++;
	XtSetArg(al[ac], XmNleftAttachment, XmATTACH_WIDGET); ac++;
	XtSetArg(al[ac], XmNleftWidget, label4); ac++;
	XtSetValues ( label5,al, ac );
	ac = 0;

	XtSetArg(al[ac], XmNtopAttachment, XmATTACH_WIDGET); ac++;
	XtSetArg(al[ac], XmNtopWidget, separator2); ac++;
	XtSetArg(al[ac], XmNleftAttachment, XmATTACH_WIDGET); ac++;
	XtSetArg(al[ac], XmNleftWidget, label5); ac++;
	XtSetArg(al[ac], XmNrightAttachment, XmATTACH_FORM); ac++;
	XtSetValues ( label6,al, ac );
	ac = 0;

	XtSetArg(al[ac], XmNtopAttachment, XmATTACH_FORM); ac++;
	XtSetArg(al[ac], XmNtopOffset, 44); ac++;
	XtSetArg(al[ac], XmNbottomAttachment, XmATTACH_WIDGET); ac++;
	XtSetArg(al[ac], XmNbottomWidget, button2); ac++;
	XtSetArg(al[ac], XmNleftAttachment, XmATTACH_FORM); ac++;
	XtSetArg(al[ac], XmNrightAttachment, XmATTACH_OPPOSITE_WIDGET); ac++;
	XtSetArg(al[ac], XmNrightWidget, label2); ac++;
	XtSetValues ( scrolledList1,al, ac );
	ac = 0;

	XtSetArg(al[ac], XmNtopAttachment, XmATTACH_WIDGET); ac++;
	XtSetArg(al[ac], XmNtopWidget, label3); ac++;
	XtSetArg(al[ac], XmNbottomAttachment, XmATTACH_WIDGET); ac++;
	XtSetArg(al[ac], XmNbottomWidget, button2); ac++;
	XtSetArg(al[ac], XmNleftAttachment, XmATTACH_WIDGET); ac++;
	XtSetArg(al[ac], XmNleftWidget, scrolledList1); ac++;
	XtSetArg(al[ac], XmNrightAttachment, XmATTACH_OPPOSITE_WIDGET); ac++;
	XtSetArg(al[ac], XmNrightWidget, label3); ac++;
	XtSetValues ( list2,al, ac );
	ac = 0;

	XtSetArg(al[ac], XmNtopAttachment, XmATTACH_WIDGET); ac++;
	XtSetArg(al[ac], XmNtopWidget, label4); ac++;
	XtSetArg(al[ac], XmNbottomAttachment, XmATTACH_WIDGET); ac++;
	XtSetArg(al[ac], XmNbottomWidget, button2); ac++;
	XtSetArg(al[ac], XmNleftAttachment, XmATTACH_WIDGET); ac++;
	XtSetArg(al[ac], XmNleftWidget, list2); ac++;
	XtSetArg(al[ac], XmNrightAttachment, XmATTACH_OPPOSITE_WIDGET); ac++;
	XtSetArg(al[ac], XmNrightWidget, label4); ac++;
	XtSetValues ( scrolledList2,al, ac );
	ac = 0;

	XtSetArg(al[ac], XmNtopAttachment, XmATTACH_WIDGET); ac++;
	XtSetArg(al[ac], XmNtopWidget, label5); ac++;
	XtSetArg(al[ac], XmNbottomAttachment, XmATTACH_WIDGET); ac++;
	XtSetArg(al[ac], XmNbottomWidget, button2); ac++;
	XtSetArg(al[ac], XmNleftAttachment, XmATTACH_WIDGET); ac++;
	XtSetArg(al[ac], XmNleftWidget, scrolledList2); ac++;
	XtSetArg(al[ac], XmNrightAttachment, XmATTACH_WIDGET); ac++;
	XtSetArg(al[ac], XmNrightWidget, label6); ac++;
	XtSetValues ( scrolledList3,al, ac );
	ac = 0;

	XtSetArg(al[ac], XmNtopAttachment, XmATTACH_WIDGET); ac++;
	XtSetArg(al[ac], XmNtopWidget, label6); ac++;
	XtSetArg(al[ac], XmNbottomAttachment, XmATTACH_WIDGET); ac++;
	XtSetArg(al[ac], XmNbottomWidget, button2); ac++;
	XtSetArg(al[ac], XmNleftAttachment, XmATTACH_WIDGET); ac++;
	XtSetArg(al[ac], XmNleftWidget, scrolledList3); ac++;
	XtSetArg(al[ac], XmNrightAttachment, XmATTACH_OPPOSITE_WIDGET); ac++;
	XtSetArg(al[ac], XmNrightWidget, label6); ac++;
	XtSetValues ( scrolledList4,al, ac );
	ac = 0;

	XtSetArg(al[ac], XmNtopAttachment, XmATTACH_NONE); ac++;
	XtSetArg(al[ac], XmNbottomAttachment, XmATTACH_FORM); ac++;
	XtSetArg(al[ac], XmNleftAttachment, XmATTACH_NONE); ac++;
	XtSetArg(al[ac], XmNrightAttachment, XmATTACH_FORM); ac++;
	XtSetValues ( button2,al, ac );

        XtAddCallback(list1, XmNbrowseSelectionCallback, ProcedureCallback, (XtPointer) 0 );
        XtAddCallback(list2, XmNbrowseSelectionCallback, ArrayNameCallback, (XtPointer) 0 ); 
         
	ac = 0;
	if (list1 != (Widget) 0) { XtManageChild(list1); }
	if (list3 != (Widget) 0) { XtManageChild(list3); }
	if (list4 != (Widget) 0) { XtManageChild(list4); }
	if (list5 != (Widget) 0) { XtManageChild(list5); }
	if ((children[ac] = separator1) != (Widget) 0) { ac++; }
	if ((children[ac] = label1) != (Widget) 0) { ac++; }
	if ((children[ac] = separator2) != (Widget) 0) { ac++; }
	if ((children[ac] = label2) != (Widget) 0) { ac++; }
	if ((children[ac] = label3) != (Widget) 0) { ac++; }
	if ((children[ac] = label4) != (Widget) 0) { ac++; }
	if ((children[ac] = label5) != (Widget) 0) { ac++; }
	if ((children[ac] = label6) != (Widget) 0) { ac++; }
	if ((children[ac] = list2) != (Widget) 0) { ac++; }
	if ((children[ac] = button2) != (Widget) 0) { ac++; }
	if (ac > 0) { XtManageChildren(children, ac); }
	ac = 0;
	if (InterRegionForm != (Widget) 0) { XtManageChild ( InterRegionForm); }
        if (InterRegionShell != (Widget) 0) { XtManageChild ( InterRegionShell); }

        CreateRegions();

        Regionswindowopen=true;
	
}

//-----------------------------
//DISPLAYS ELEMENTS IN THE GUI
//-----------------------------

void AddElement2List(char* addtext, Widget w)
{
    XmString xmstringt;
    xmstringt = XmStringCreateLtoR (addtext , (XmStringCharSet)XmFONTLIST_DEFAULT_TAG );
    XmListAddItem(w,xmstringt,0);
    XmStringFree(xmstringt);
  }
   
//---------------------------------------
//FUNCTION TO PRINT THE SQL QUERY RESULTS
//--------------------------------------- 
vector<string> Print_Query_Result(int num_fields, vector <string> proc_lst) 
{
     MYSQL_ROW row;
     int i, j=0;

     while ((row = mysql_fetch_row(result)) != NULL) 
     {
	 for(i = 0; i <num_fields; i++) 
	 {
	    proc_lst.push_back(row[i]);
	    j++;
	 }
	   
     }
      return proc_lst;
}


//-----------------------------------------------
//EVENT ON MOUSE CLICK ON A PROCEDURE/SCOPE NAME
//-----------------------------------------------
void ProcedureCallback(Widget w, XtPointer client_data, XtPointer call_data)
{
  bool debug = false;
  unsigned int num_fields;
  unsigned int i,j;
  vector <string> proc_lst2;
  char *SelectedListItem; //contains the string of the selected call;
  XmListCallbackStruct *acs = (XmListCallbackStruct *)call_data;
  XmStringGetLtoR(acs->item, XmFONTLIST_DEFAULT_TAG, &SelectedListItem);

  if (debug) cout <<SelectedListItem; 

  if(strcmp(SelectedListItem,"Global") == 0) 
      strcpy(selected_proc,"@");
  else 
      strcpy(selected_proc,SelectedListItem);

  string query_s="select distinct arr_var_name from array_reg_table WHERE proc_name=' "; 
  query_s+=selected_proc;
  query_s+="';";
  mysql_real_query(connection,(char*)query_s.c_str(), 555);
  result = mysql_store_result(connection);

  XmListDeleteAllItems(list2);
  XmListDeleteAllItems(list3);
  XmListDeleteAllItems(list4);
  XmListDeleteAllItems(list5);
	
  int num_rows = mysql_num_rows(result);
  if (num_rows == 0) 
          AddElement2List("No array vars!",list2);
  else 
  {
  	  num_fields = mysql_field_count(&mysql);	
	  proc_lst2 = Print_Query_Result(num_fields, proc_lst2);	 
	  for(int j=0; j<proc_lst2.size();j++) 
	  { 
	    AddElement2List((char *)proc_lst2[j].c_str(),list2);
	  }
  }
  mysql_free_result(result);

}

//--------------------------------------
//EVENT ON MOUSE CLICK ON AN ARRAY NAME
//--------------------------------------

void ArrayNameCallback(Widget w, XtPointer client_data, XtPointer call_data)
{
  bool debug = false;
  unsigned int num_fields;
  unsigned int i,j;
  vector <string> proc_lst3;
  char *query;
  char *SelectedListItem; //contains the string of the selected call;

  XmListCallbackStruct *acs = (XmListCallbackStruct *)call_data;
  XmStringGetLtoR(acs->item, XmFONTLIST_DEFAULT_TAG, &SelectedListItem);
  
   string query_s = "select DISTINCT num_dims,reg_usage from array_reg_table WHERE proc_name=' ";
  query_s+=selected_proc;
  query_s+="' AND arr_var_name='";
  query_s+=SelectedListItem;
  query_s+="';";
 
  mysql_real_query(connection, (char *)query_s.c_str(), 255);
  result = mysql_store_result(connection);
  num_fields = mysql_field_count(&mysql);
  proc_lst3 = Print_Query_Result(num_fields, proc_lst3);

  XmListDeleteAllItems(list3);
  XmListDeleteAllItems(list4);
  XmListDeleteAllItems(list5);

  for(int j=0; j<proc_lst3.size();j++) { 
     AddElement2List((char *)proc_lst3[j].c_str(),list3);
     AddElement2List((char *)proc_lst3[++j].c_str(),list4);
  }
  mysql_free_result(result);

    string query_s1 = "select SUM((preg_end_id - preg_start_id)+1) AS reg_access from array_reg_table WHERE proc_name=' ";
  query_s1+=selected_proc;
  query_s1+="' AND arr_var_name='";
  query_s1+=SelectedListItem;
  query_s1+="';";
 
  mysql_real_query(connection, (char *)query_s1.c_str(), 255);
  result = mysql_store_result(connection);
  row = mysql_fetch_row(result);  
  AddElement2List(row[0],list5);
  mysql_free_result(result);
}

//------------------------------------------
//C API TO CONNECT TO THE MYSQL DATABASE
//-------------------------------------------

void Mysql_C_API(char *datafile, char *current_procedure_name) {

        unsigned int num_fields;
        unsigned int i,j;
        unsigned long *lengths;

        mysql_init(&mysql);
        connection =mysql_real_connect(&mysql,"localhost","dragon","dragon","mysql",0,0,0);

        if(connection == NULL) {
            printf(mysql_error(&mysql));
        }

	mysql_real_query(connection,"CREATE DATABASE array_database;",255);
	mysql_real_query(connection,"USE array_database;",255);

	mysql_real_query(connection,"DROP TABLE IF EXISTS array_reg_table;",255);

	mysql_real_query(connection,"CREATE TABLE array_reg_table (file_id INT, reg_id int not null, PRIMARY KEY(file_id,reg_id), proc_name char(50), arr_var_name VARCHAR(50), reg_usage VARCHAR(50), preg_start_id INT, preg_end_id INT, num_dims INT );", 555);
	
        string temp = " LOAD DATA INFILE '";
        temp = temp + datafile;
        temp = temp + "' INTO TABLE array_reg_table FIELDS TERMINATED BY ',';";
     
	mysql_real_query(connection,(char *)temp.c_str(),255);
	
	XmListDeleteAllItems(list1);
	//AddElement2List("rhs_norm_",list1);
	AddElement2List(current_procedure_name,list1);
	AddElement2List("Global",list1);
	
	
}

//-------------------

void CreateRegions(void)
{
        char *PU_Name, *FileName;
	int len=0,n;
        bool debug = false;
   
        len=strlen(UICallgraphRoot[ProcedureLoaded].GetName());
	PU_Name = new char[len+1];
	strcpy(PU_Name, UICallgraphRoot[ProcedureLoaded].GetName());
	len = strlen(SelectedFile);
	FileName = new char[len+5];
 	strcpy(FileName,SelectedFile);
	for(n=len-1; n>0; n--)
		if(FileName[n]!='.')
			FileName[n]='\0';
		else
			break;
	if(n>0)
		strcat(FileName,"rgn");
	else{
		strcpy(FileName,SelectedFile);
		strcat(FileName,".rgn");
	}  

	if (debug)
	  {
	    cout << "\nProcedure Loaded: "<< FileName << " PU: "<<PU_Name; 
            
          }
        Mysql_C_API(FileName, PU_Name);
  
 // the filename that is going to read is the same as the project filename with extension .rgn
 // 1. Read File
 // 2. Add elements to columns.
 // 3. Fix the GUI;  
  
 
  delete [] FileName;
  delete [] PU_Name; 
}

