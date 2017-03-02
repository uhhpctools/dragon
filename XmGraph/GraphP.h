#ifndef _XmGraphP_h
#define _XmGraphP_h

/*
 * $Id: GraphP.h,v 1.1.1.1 2004/08/09 20:03:03 liao Exp $
 *
 * $Log: GraphP.h,v $
 * Revision 1.1.1.1  2004/08/09 20:03:03  liao
 * Imported Dragon from Lurch
 *
 * Revision 1.1.1.1  2002/09/10 22:45:19  sgicomp
 * Create the first version of the Dragon
 *
 * Revision 1.4  1996/07/01 15:34:34  markus
 * added an #include <Xm/Xm.h>
 * added an #include <Xm/XmP.h> and #include <Xm/ManagerP.h> for SUNOS 4.XX
 * changed include "Graph.h" to include "XmGraph/Graph.h"
 * added _XmGraphPart to the appropriate typedef struct...
 *
 * Revision 1.3  1996/06/27 18:55:28  bob
 * Checked typo.
 *
 * Revision 1.2  1996/06/27 18:15:05  bob
 * Changed typedef typedef struct _XmGraphPart.
 * Inserted #include <Xm/Xm.h> and #include <Xm/ManagerP.h>.
 * Changed #include "XmGraph/Graph.h" and #include "XmGrpah/ArcP.h".
 *
 * Revision 1.1  1996/02/20 09:57:37  martin
 * Initial revision
 *
 */

/*
 *  File:         GraphP.h
 *  Description:  Private header file for GraphWidget.
 *  Author:       Luis Miguel (luis@postgres.berkeley.edu)
 *  
 */

/* CYY modified 2/1/91 for  R3 and R4 compatible and ANSI_C fixes*/

#include <Xm/Xm.h>

#ifdef SUNOS4
#include <Xm/XmP.h>
#endif

//#ifndef SUNOS4
#include <Xm/ManagerP.h>
//#endif

#include "XmGraph/Graph.h"
#include "XmGraph/ArcP.h"


/* useful macros */
#ifndef MAX
#define MAX(a,b)        ((a) > (b) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a,b)        ((a) < (b) ? (a) : (b))
#endif

/* some constants */

#define  LARGE_NUM           987654321

typedef enum  {NORMAL, 
NODE_INDICATED, 
NODE_INDICATED_PENDING_CANCEL, 
NODES_INDICATED, 
NODES_INDICATED_PENDING_CANCEL, 
ARC_INDICATED, 
ARC_INDICATED_PENDING_CANCEL, 
ARCS_INDICATED, 
ARCS_INDICATED_PENDING_CANCEL, 
NODE_SELECTED_FOR_MOTION, 
ARC_SELECTED_FOR_MOTION,
MULTIPLE_NODES_SELECTED_FOR_MOTION, 
NODE_SELECTED,
NODES_SELECTED, 
ARC_SELECTED, 
ARCS_SELECTED,
POSSIBLE_DOUBLE_CLICK,
SUBGRAPH_INDICATED,
SUBGRAPH_INDICATED_PENDING_CANCEL,
MOVING_ARC, 
MOVING_ARCS, 
MOVING_NODE, 
MOVING_NODES, 
ADDING_NODE, 
ADDING_ARC,
ADDING_ARC_IN_PARENT,
REGION_INDICATED} ActionMode;

/* The node structure "surrounding" widgets */

typedef struct _node *NodePtr;
typedef NodePtr *NList;

typedef struct nodelist {
    NList      nodes;
    int        n_nodes;
    int        n_slots;
} NodeList;


typedef struct _node {
    /* all node and arc widgets */
    Widget      widget;
    ArcList     from_arcs;
    ArcList     to_arcs;
    
    /* the structure of the graph */
    NodeList    parents;
    NodeList    kids;
    int         visited;  /* the visited flag can contain a Boolean or count */
    
    /* tree processing needs these */
    int         level;
    NodePtr     tree_parent;
    NodeList    tree_kids;
    
    
} Node;

/* New fields for the GraphWidget widget class record */

typedef struct _XmGraphClassPart {
    XtPointer extension;
} XmGraphClassPart;

/* Full Class record declaration */

typedef struct _XmGraphClassRec {
    CoreClassPart	core_class;
    CompositeClassPart  composite_class;
    ConstraintClassPart	constraint_class;
    XmManagerClassPart  manager_class;
    XmGraphClassPart	graph_class;
} XmGraphClassRec;

externalref XmGraphClassRec xmGraphClassRec;

/* New Fields for the GraphWidget  widget record */


typedef struct _XmGraphPart {
    ActionMode     current_action;
    Time           last_button_time;
    int            double_click_interval;
    Cursor         indicate_child_cursor;
    Cursor         indicate_parent_cursor;
    Cursor         ptr_cursor;
    Cursor         motion_cursor;
    Cursor         indicate_cursor;
    Cursor         high_right_cursor;
    Cursor         high_left_cursor;
    Cursor         low_right_cursor;
    Cursor         low_left_cursor;
    Pixel          foreground;
    int            start_x, start_y, end_x, end_y;
    int            delta_x, delta_y;
    int            child_spacing;
    int            sibling_spacing;
    int            original_child_spacing;
    int            original_sibling_spacing;
    GC             gc;
    GC             xorgc;
    GC             cleargc;
    Boolean        edit_mode;
    Boolean        nodes_movable;
    autoLayoutType auto_layout_mode;
    Boolean        siblings_visible;
    Boolean        re_layout;
    Boolean	   re_orient;
    Boolean	   layed_out;
    Boolean        batch_drawing_mode;
    unsigned char  direction;
    unsigned char  arc_draw_mode;
    XmArcWidgetList     arcs;
    ArcList        selected_arcs;
    NodeList       selected_nodes;
    NodeList       user_roots;
    int            n_arcs;
    int            n_arc_slots;
    int            n_nodes;
    NodePtr        root;
    NodePtr        current_node;
    NodePtr        current_subgraph;
    Widget         current_arc;
    Widget         highlighted_arc;
    Widget         indicated_widget;
    WidgetClass    default_widget_class;
    Boolean allow_multiple_selections;
    void  (*layout_proc)();
    unsigned char  arc_style;
    
    XtCallbackList new_arc;  
    XtCallbackList new_node;  
    XtCallbackList node_moved;
    XtCallbackList arc_moved;
    XtCallbackList double_click;  
    XtCallbackList select_node;  
    XtCallbackList deselect;  
    XtCallbackList select_arc;  
    XtCallbackList select_subgraph;  
    XtCallbackList constraint_free;  
    XtCallbackList constraint_init;  

    Boolean         is_scrolled;
    Boolean         show_crossing_arcs;
    Dimension       clip_width;
    Dimension       clip_height;

    Boolean         deferedChangeManaged;	/* added to defer all calls
						 to ChangeManaged until the
						 graph is visible. */
    Boolean         deferedAdjustSize;

    Boolean         snap_on;
    int 	    snap_grid_size;
} XmGraphPart;


/**************************************************************************
 *
 * Full instance record declaration
 *
 **************************************************************************/

typedef struct _XmGraphRec {
    CorePart	    core;
    CompositePart   composite;
    ConstraintPart  constraint;
    XmManagerPart   manager;
    XmGraphPart     graph;
}  XmGraphRec;

/* Graph constraint record */

typedef struct _XmGraphConstraintsRec {
    XmManagerConstraintPart manager;
    int           old_mask;
    NodePtr       node;    
    int           delta_x;
    int           delta_y;
    int           old_x;
    int           old_y;
/* CYY added */
    Boolean 	  pos_fixed;
    double	  pos_x;
    double	  pos_y;
    double	  vel_x;
    double	  vel_y;
/* CYY added end */

    XtPointer       extension;
} XmGraphConstraintsRec, *XmGraphConstraintPtr;




#ifdef _R4_
#define NODEPTR(widget) ( (widget) ?	\
			 ( XtIsWidget(widget) ? \
			  ((XmGraphConstraintsRec *) ((Widget)(widget))->core.constraints)->node   :    \
			  ((XmGraphConstraintsRec *) ((Object)(widget))->object.constraints)->node)  :  \
			 NULL)
#define CONSTRAINTREC(widget) ( (widget) ?                                     \
			       ( XtIsWidget(widget)  ?		      \
				((XmGraphConstraintsRec *) ((Widget)(widget))->core.constraints)   :     \
				((XmGraphConstraintsRec *) ((Object)(widget))->object.constraints))  :  \
			       NULL)
#else

/* ObjectP and  CoreP are identical up to constraints */
#define NODEPTR(w) ( (w) ? ((XmGraphConstraintsRec *) ((w)->core.constraints))->node : NULL)

#define CONSTRAINTREC(w) ((w)? ( (XmGraphConstraintsRec *) w->core.constraints):NULL)


#endif

/* functions exported by Graph.c for internal use only */
#ifndef _NO_PROTO
extern void   _GetPoints	(XmGraphWidget, int,int,int,int,
						int,int,int,int,
				           	float*, float *,
						float *, float *); /* Arc.c */
extern void   _InitArcList	(XmGraphWidget);    /* Arc.c */
#else
extern void   _GetPoints	();	/* Arc.c */
extern void   _InitArcList	();     /* Arc.c */
#endif

#endif /* _XmGraphP_h  */
/* DON"T ADD ANYTHING AFTER THIS endif */


