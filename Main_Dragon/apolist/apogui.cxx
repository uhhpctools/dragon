#include <stdlib.h>
#include <stdio.h>
#include <X11/Intrinsic.h>
#include <X11/Shell.h>
#include <Xm/Xm.h>
#include <Xm/DialogS.h>
#include <Xm/Form.h>
#include <Xm/Label.h>
#include <Xm/List.h>
#include <Xm/PushB.h>
#include <Xm/ScrollBar.h>
#include <Xm/Text.h>
#include <Xm/Protocols.h>
#include <Xm/PanedW.h>
#include <Xm/RowColumn.h>
#include "apogen.h"

static Widget wAPO = (Widget) NULL;
static Widget wList = (Widget) NULL;
static Widget wListLines = (Widget) NULL;
static Widget wText = (Widget) NULL;
static VecParStatus parstatus;
static int iPUcode = 0;

extern void showScrollPosition(unsigned long lLineNumber);

/*------------------------------------------
 * Close the window
 * */
void CloseAPOWindow(Widget CloseButtonWidget, XtPointer, XtPointer) {
   XtUnmanageChild(wAPO);
}

/* ------------------------------------
 * selection callback
 * ------------------------------------ */
void selectLineCB(Widget w, XtPointer client_data, XtPointer call_data)
{
   int iLine, i;
   XmListCallbackStruct *list_cbs = (XmListCallbackStruct*) call_data;
   char sTmp[1024];
   iLine = list_cbs->item_position  - 1;
   sprintf(sTmp, "Status: %s\n %s", parstatus[iPUcode].codeStatus[iLine].detailstatus.c_str(),
		   parstatus[iPUcode].codeStatus[iLine].info.c_str());
   XmTextSetString(wText, (char *) sTmp);
   showScrollPosition(parstatus[iPUcode].codeStatus[iLine].iLineNoStart);
}


/* ------------------------------------
 * selection callback
 * ------------------------------------ */
void selectCB(Widget w, XtPointer client_data, XtPointer call_data)
{
    int x, y;
    int mem_allocated;
    char sTmp[1024];
    XmString s;
    XmListCallbackStruct *list_cbs = (XmListCallbackStruct*) call_data;
    
    x = list_cbs->item_position  - 1; // we start from zero while list_cbs starts from 1
    iPUcode = x;
    XmListDeleteAllItems(wListLines);
    for(y=0;y<parstatus[x].codeStatus.size();y++) {
        sprintf(sTmp, "line %d-%d", parstatus[x].codeStatus[y].iLineNoStart,
		 	 parstatus[x].codeStatus[y].iLineNoEnd);
        s = XmStringCreate((char *) sTmp, XmSTRING_DEFAULT_CHARSET);
	XmListAddItem(wListLines, s, 0); 
    }
    XmTextSetString(wText,"");
}

/*-----------------------------------------------------------
 * Add array of items into the list
 * ------------------------------------------------------------------*/
void AddItems(VecParStatus vecparstatus)
{
    XmString s;
    int list_cnt;
    char sBuf[1024];

    for (list_cnt=0; list_cnt<vecparstatus.size(); list_cnt++)
    {
        s = XmStringCreate((char *)vecparstatus[list_cnt].PUname.c_str(), XmSTRING_DEFAULT_CHARSET);
        XmListAddItem(wList,s,0);
        XmStringFree(s);
    }
}

/*---------------------------------------------
 *  * User interface for OpenMP parallelization
 *   * This will list all loops (SNL) and their status
 *    * If the loops are not parallelizables, show the details
 *     * --------------------------------------------*/
Widget createAPOwindow(Widget parent) {
   register int ac = 0;           /* Arg Count */
   Arg al[64];                    /* Arg List */
   Widget wPanel;
   
   XtSetArg(al[ac], XmNallowShellResize, TRUE); ac++;
   XtSetArg(al[ac], XmNtitle, "Parallelization status"); ac++;
   XtSetArg(al[ac], XmNheight, 400); ac++;
   XtSetArg(al[ac], XmNwidth,700); ac++;

   wAPO = XtCreatePopupShell ( (char *) "APOshell", topLevelShellWidgetClass, parent, al, ac );

   Atom WmDeleteWindow = XmInternAtom(XtDisplay(wAPO),"WM_DELETE_WINDOW",False);
   XmAddWMProtocolCallback(wAPO, WmDeleteWindow, CloseAPOWindow, NULL);

    /* create the paned window container. */
   ac = 0;
   XtSetArg(al[0], XmNorientation, XmHORIZONTAL); ac++;
/*   XtSetArg(al[ac], XmNtopAttachment, XmATTACH_FORM);ac++;
   XtSetArg(al[ac], XmNbottomAttachment, XmATTACH_FORM);ac++;
   XtSetArg(al[ac], XmNleftAttachment, XmATTACH_FORM);ac++;
   XtSetArg(al[ac], XmNrightAttachment, XmATTACH_FORM);ac++;
   */
   XtSetArg(al[ac], XmNautoUnmanage, FALSE); ac++;
   XtSetArg(al[ac], XmNheight, 400); ac++;
   XtSetArg(al[ac], XmNwidth, 700); ac++;
   wPanel = XmCreateForm (wAPO, "panel",  al, ac); 
   //wPanel = XmCreateRowColumn (wAPO, "rowcolumn",  al, ac); 
   XtManageChild(wPanel);

    /* create a list widget */
   ac=0;
   XtSetArg(al[ac], XmNselectionPolicy, XmSINGLE_SELECT); ac++;
   XtSetArg(al[ac], XmNvisibleItemCount, 10); ac++;
   XtSetArg(al[ac], XmNtopAttachment, XmATTACH_FORM);ac++;
   XtSetArg(al[ac], XmNbottomAttachment, XmATTACH_FORM);ac++;
   XtSetArg(al[ac], XmNleftAttachment, XmATTACH_FORM);ac++;
   wList=XmCreateScrolledList(wPanel,"list",al,ac);
   XtManageChild(wList);
   XtAddCallback (wList, XmNsingleSelectionCallback, selectCB, NULL);

    /* create a list line widget */
   ac=0;
   XtSetArg(al[ac], XmNselectionPolicy, XmSINGLE_SELECT); ac++;
   XtSetArg(al[ac], XmNvisibleItemCount, 10); ac++;
   XtSetArg(al[ac], XmNtopAttachment, XmATTACH_FORM);ac++;
   XtSetArg(al[ac], XmNbottomAttachment, XmATTACH_FORM);ac++;
   XtSetArg(al[ac], XmNleftAttachment, XmATTACH_WIDGET);ac++;
   XtSetArg(al[ac], XmNleftWidget, wList);ac++;
   wListLines = XmCreateScrolledList(wPanel,"listline",al,ac);
   XtManageChild(wListLines);
   XtAddCallback (wListLines, XmNsingleSelectionCallback, selectLineCB, NULL);

    /* create a text widget */
   ac=0;
   XtSetArg(al[ac], XmNeditMode, XmMULTI_LINE_EDIT ); 	ac++;
   XtSetArg(al[ac], XmNeditable, False);		ac++;
   XtSetArg(al[ac], XmNtopAttachment, XmATTACH_FORM);ac++;
   XtSetArg(al[ac], XmNbottomAttachment, XmATTACH_FORM);ac++;
   XtSetArg(al[ac], XmNrightAttachment, XmATTACH_FORM);ac++;
   XtSetArg(al[ac], XmNleftAttachment, XmATTACH_WIDGET);ac++;
   XtSetArg(al[ac], XmNleftWidget, wListLines);ac++;
   
   wText=XmCreateScrolledText(wPanel,"textinfo",al,ac);
   XtManageChild(wText);

   return wAPO;
}

/* ---------------------------------
 * Show APO Windo
 * --------------*/
void ShowAPOwindow(Widget wParent, VecParStatus vecparstatus) {
  if(wAPO) {
  } else {
     createAPOwindow(wParent);
     XtManageChild(wAPO);
  }
  AddItems(vecparstatus);
  XtManageChild(wAPO);
  parstatus = vecparstatus;
}


