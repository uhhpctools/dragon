/*
 * $Id: Graph.h,v 1.1.1.1 2004/08/09 20:03:03 liao Exp $
 *
 * $Log: Graph.h,v $
 * Revision 1.1.1.1  2004/08/09 20:03:03  liao
 * Imported Dragon from Lurch
 *
 * Revision 1.1.1.1  2002/09/10 22:45:19  sgicomp
 * Create the first version of the Dragon
 *
 * Revision 1.3  1996/06/27 18:12:14  bob
 * Inserted typedef struct _XmGraphCallbackStruct.
 *
 * Revision 1.2  1996/06/14 21:08:38  markus
 * just a tiny change...
 *
 * Revision 1.1  1996/02/20 09:57:33  martin
 * Initial revision
 *
 */

/*
 *   File:        Graph.h
 *
 *   Project:     Motif -- OSF X User Environment Widget Set
 *
 *   Description: Public include file for the XmLabel class
 *
 *
 */ 

#ifndef _XmGraph_h
#define _XmGraph_h

/* Class record constant */

externalref WidgetClass xmGraphWidgetClass;
typedef struct _XmGraphClassRec		*XmGraphWidgetClass;
typedef struct _XmGraphRec		*XmGraphWidget;


#include <Xm/Xm.h>
#ifdef _XFWFINTERFACE
#include "Xfwf/Slider2.h"
#endif
#ifndef XmIsGraph
#define XmIsGraph (w)	XtIsSubclass(w, xmGraphWidgetClass)
#endif


/******************* CONVENIENCE FUNCTIONS **********************/

/* creates and returns an instance of the graph widget */

#ifdef __cplusplus                    /* do not leave open across includes */
extern "C" {                                  /* for C++ V2.0 */
#endif

#ifdef _XFWFINTERFACE
extern void XfwfResizeThumb_Interface(Widget wt, double parm1, double parm2);
extern void XfwfGetThumb_Interface(Widget, XfwfScrollInfo *);
extern void XfwfMoveThumb_Interface(Widget, double, double);
#endif

#ifndef _NO_PROTO

extern void	XmGraphUpdateSpringLayout (Widget w);
extern void	XmGraphSpringLayout (Widget w,Widget ww);
extern void	XmGraphSpringLayoutList (Widget w,Widget *wlist,int size);

extern void	XmGraphFixAllNodePos (Widget , Boolean );

extern void	XmGraphCenterAll(Widget);
extern void	XmGraphCenterAtPos(Widget,int,int);
extern void	XmGraphFixNodePos (Widget , Widget , Boolean );
extern void	XmGraphFixSelectedNodePos (Widget , Boolean );
extern WidgetList XmGraphGetConnectedNodes(Widget , Widget , int *);
extern void	XmGraphGetConsts (Widget ,double *len, double *springK,
				double *chargeK, double *loss);
extern int	XmGraphGetGridSize(Widget);
extern void	XmGraphSetSnapGrid (Widget , Boolean );
extern void	XmGraphSetGridSize(Widget,Cardinal);
extern void	XmGraphSetConsts (Widget,double len, double springK,
				double chargeK, double loss);
extern void     XmGraphZoom(Widget , float);
#else

extern void	XmGraphUpdateSpringLayout ();
extern void	XmGraphSpringLayout ();
extern void	XmGraphSpringLayoutList ();

extern void	XmGraphFixAllNodePos ();
extern void	XmGraphCenterAll();
extern void	XmGraphCenterAtPos();
extern void	XmGraphFixNodePos ();
extern void	XmGraphFixSelectedNodePos ();
extern WidgetList XmGraphGetConnectedNodes();
extern void	XmGraphGetConsts ();
extern int	XmGraphGetGridSize();
extern void	XmGraphSetSnapGrid ();
extern void	XmGraphSetGridSize();
extern void	XmGraphSetConsts ();
#endif

#ifndef _NO_PROTO
extern void XmGraphSetScrolled(Widget,int);

extern Widget	XmCreateGraph(Widget,String,ArgList,Cardinal);
extern Widget	XmCreateManagedGraph (Widget,String,ArgList,Cardinal);
extern Widget	XmCreateScrolledGraph(Widget,String,ArgList,Cardinal);

extern  void       XmGraphAddNode(Widget);
extern  void       XmGraphCenterAroundWidget(Widget, Widget);
extern  void       XmGraphDestroyAllArcs(Widget);
extern  void       XmGraphDestroyAllNodes(Widget);
extern  void       XmGraphDestroySelectedArcsOrNodes(Widget);
extern  void       XmGraphGetArcNodes(Widget, Widget, Widget *, Widget *);
extern  WidgetList XmGraphGetArcsBetweenNodes(Widget, Widget, Widget, int *);
extern  WidgetList XmGraphGetArcs(Widget, int *);
extern  void	   XmGraphGetNodeArcs(Widget, Widget,
				   WidgetList *, WidgetList *, int *, int *);
extern  WidgetList XmGraphGetNodes(Widget, int *);
extern  WidgetList XmGraphGetRoots(Widget, int *);
extern  WidgetList XmGraphGetSelectedArcs(Widget,int *);
extern  WidgetList XmGraphGetSelectedNodes(Widget,int *);
extern  Widget     XmGraphInputOverArc(Widget, int, int);
extern  void       XmGraphInsertRoots(Widget, WidgetList, int);
extern  Boolean    XmGraphIsPointInArc(Widget, int, int);
extern  Boolean    XmGraphIsSelectedArc(Widget, Widget);
extern  Boolean    XmGraphIsSelectedNode(Widget, Widget);
extern  void       XmGraphLayout(Widget);
extern  void	   XmGraphMoveAll(Widget, int, int);
extern  Boolean    XmGraphMoveArc(Widget, Widget, Widget, Widget);
extern  Boolean    XmGraphMoveNode(Widget,Widget, Position, Position);
extern  int        XmGraphNumArcs(Widget);
extern  int        XmGraphNumNodes(Widget);
extern  int        XmGraphNumNodeArcs(Widget, Widget, int * , int *);
extern  int        XmGraphNumRoots(Widget);
extern  int        XmGraphNumSelectedArcs(Widget);
extern  int        XmGraphNumSelectedNodes(Widget);
extern  void       XmGraphRelaySubgraph(Widget, Widget);
extern  void       XmGraphRemoveArcBetweenNodes(Widget, Widget, Widget);
extern  void       XmGraphRemoveRoots(Widget, WidgetList, int);
extern  void       XmGraphSelectNode(Widget,Widget);
extern  void       XmGraphSelectNodes(Widget, WidgetList, int);
extern  void       XmGraphSelectArc(Widget, Widget);
extern  void       XmGraphSelectArcs(Widget, WidgetList, int);
extern  void       XmGraphUnselectArc(Widget,Widget);
extern  void       XmGraphUnselectArcs(Widget, WidgetList, int);
extern  void       XmGraphUnselectNode(Widget,Widget);
extern  void       XmGraphUnselectNodes(Widget,WidgetList, int);
extern  void	   XmGraphUnmap(Widget);
extern  void	   XmGraphMap(Widget);

#else    /* _NO_PROTO */
extern void XmGraphSetScrolled();

extern Widget XmCreateGraph();
extern Widget XmCreateManagedGraph ();
extern Widget XmCreateScrolledGraph();
extern Widget XmCreateArc ();
extern Widget XmCreateAttachedArc ();

extern  void    XmGraphGetArcNodes();
extern  void    XmGraphGetNodeArcs();
extern  Widget  XmGraphInputOverArc();
extern  void    XmGraphInsertRoots();
extern  Boolean XmGraphIsSelectedArc();
extern  Boolean XmGraphMoveArc();
extern  Boolean XmGraphMoveNode();
extern  void    XmGraphSelectNode();
extern  void    XmGraphSelectNodes();
extern  void    XmGraphSelectArc();
extern  void    XmGraphSelectArcs();

extern  void	XmGraphCenterAroundWidget();
extern  void	XmGraphDestroyAllArcs();
extern  void	XmGraphDestroyAllNodes();
extern  void	XmGraphDestroySelectedArcsOrNodes();
extern  WidgetList	XmGraphGetArcsBetweenNodes();
extern  WidgetList	XmGraphGetArcs();
extern  WidgetList	XmGraphGetNodes();
extern  WidgetList	XmGraphGetRoots();
extern  WidgetList	XmGraphGetSelectedArcs();
extern  WidgetList	XmGraphGetSelectedNodes();
extern  Boolean XmGraphIsPointInArc();
extern  Boolean XmGraphIsSelectedNode();
extern  int	XmGraphNumArcs();
extern  int	XmGraphNumNodes();
extern  int	XmGraphNumNodeArcs();
extern  int	XmGraphNumRoots();
extern  int	XmGraphNumSelectedArcs();
extern  int	XmGraphNumSelectedNodes();
extern  void	XmGraphMoveAll();
extern  void	XmGraphLayout();
extern  void	XmGraphRelaySubgraph();
extern  void	XmGraphRemoveArcBetweenNodes();
extern  void	XmGraphRemoveRoots();
extern  void	XmGraphUnselectArc();
extern  void	XmGraphUnselectArcs();
extern  void	XmGraphUnselectNode();
extern  void	XmGraphUnselectNodes();
extern  void	XmGraphAddNode();
extern  void	XmGraphUnmap();
extern  void	XmGraphMap();
#endif /*_NO_PROTO */

#ifdef __cplusplus
}                                             /* for C++ V2.0 */
#endif

#define XmNsnapGridSize     "snapGridSize"
#define XmCSnapGridSize     "SnapGridSize"

#define XmNsnapGridOn     "snapGridOn"
#define XmCSnapGridOn     "SnapGridOn"

#define XmNtwinsVisible     "twinsVisible"
#define XmCTwinsVisible     "TwinsVisible"

#define XmNarcDrawMode      "arcDrawMode"
#define XmCArcDrawMode      "ArcDrawMode"


#define XmNautoLayoutMode   "autoLayoutMode"
#define XmCAutoLayoutMode   "AutoLayoutMode"

#define XmNreLayout	    "reLayout"
#define XmCReLayout	    "ReLayout"

#define XmNreorient         "reorient"
#define XmCReorient         "Reorient"

#define XmNchildSpacing     "childSpacing"
#define XmCChildSpacing     "ChildSpacing"

#define XmNsiblingSpacing   "siblingSpacing"
#define XmCSiblingSpacing   "SiblingSpacing"

#define XmNsaveGraph       "saveGraph"
#define XmCSaveGraph       "SaveGraph"

#define XmNsaveFileName     "saveFileName"
#define XmCSaveFileName     "SaveFileName"

	/* these have class XmCCallback  */

#define XmNnewArcCallback           "newArcCallback"

#define XmNnewNodeCallback	    "newNodeCallback"
#define XmNallowMultipleSelections  "allowMultipleSelections"
#define  XmCAllowMultipleSelections "AllowMultipleSelections"

#define XmNnodeMovedCallback	    "nodeMovedCallback"

#define XmNarcMovedCallback	    "arcMovedCallback"

#define XmNselectNodeCallback	    "selectNodeCallback"

#define XmNselectArcCallback	    "selectArcCallback"

#define XmNdeselectCallback	    "deselectCallback"

#define XmNselectSubgraphCallback   "selectSubgraphCallback"

#define XmNdeleteNodeCallback	    "deleteNodeCallback"

#define XmNdeleteArcCallback	    "deleteArcCallback"

#define XmNdefaultNodeClass "defaultNodeClass"

#define XmCDefaultNodeClass "DefaultNodeClass"

#define XmNdefaultLabel     "defaultLabel"

#define XmRArcDrawMode      "ArcDrawMode"

#define XmRAutoLayoutType     "AutoLayoutType"

#define XmNinteractiveArcDirection     "interactiveArcDirection"

#define XmCInteractiveArcDirection     "InteractiveArcDirection"

#define XmNmovableNodes   "movableNodes"
#define XmCMovableNodes    "MovableNodes"

#define XmNshowCrossingArcs  "showCrossingArcs"
#define XmCShowCrossingArcs  "ShowCrossingArcs"

#define XmNptrCursor  "ptrCursor"
#define XmNmotionCursor  "motionCursor"
#define XmNindicateCursor  "indicateCursor"
#define XmNindicateChildCursor  "indicateChildCursor"
#define XmNindicateParentCursor "indicateParentCursor"

#define XmNinitializeDataCallback  "initializeDataCallback"
#define XmNfreeDataCallback  "freeDataCallback"

#define XmNlayoutProc  "layoutProc"
#define XmCLayoutProc  "LayoutProc"


/****************************************************
 * Callback reasons
 *****************************************************/

#define XmCR_NEW_ARC			41
#define XmCR_NEW_NODE                   42
#define XmCR_NODE_MOVED			43
#define XmCR_ARC_MOVED			44
#define XmCR_SUBGRAPH_MOVED		45
#define XmCR_ARC_EDITED			46
#define XmCR_SELECT_NODE		47
#define XmCR_SELECT_ARC			48
#define XmCR_SELECT_SUBGRAPH		49
#define XmCR_DELETE_NODE		50
#define XmCR_DELETE_ARC			51

#ifndef XmCR_SELECT
#define XmCR_SELECT			52
#endif

#define XmCR_RELEASE			53
#define XmCR_NODE_DOUBLE_CLICK		54
#define XmCR_ARC_DOUBLE_CLICK		55
#define XmCR_DOUBLE_CLICK		56
#define XmCR_DESELECT_ARC		57
#define XmCR_DESELECT_NODE		58
#define XmCR_DESELECT		        59
#define XmCR_NODES_MOVED		60
#define XmCR_SELECT_NODES		61
#define XmCR_SELECT_ARCS		62
#define XmCR_SELECT_ARCS_AND_NODES		63

#define XmPOSITION_FIXED                 0
#define XmPOSITION_RELATIVE              1

typedef struct _XmGraphCallbackStruct {
    int            reason;
    XEvent        *event;
    Boolean        interactive;
    WidgetList     selected_widgets;
    int            num_selected_widgets;
    WidgetList     selected_arcs;
    int            num_selected_arcs;
    Widget         widget;
    Widget         old_to;    /* Used for move and edit arc callbacks */
    Widget         old_from;
    Widget         new_to;
    Widget         new_from;
    Boolean        doit;
} XmGraphCallbackStruct;

typedef enum  {XmNEVER, XmALWAYS, XmARCS_ONLY, XmNODES_ONLY, XmPARTIAL} autoLayoutType;

#endif /* _XmGraph_h */


/* DO NOT ADD ANYTHING AFTER THIS #endif */



