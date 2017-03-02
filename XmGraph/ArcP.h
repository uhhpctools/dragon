#ifndef ARCP_H
#define ARCP_H

/*
 * $Id: ArcP.h,v 1.1.1.1 2004/08/09 20:03:03 liao Exp $
 *
 * $Log: ArcP.h,v $
 * Revision 1.1.1.1  2004/08/09 20:03:03  liao
 * Imported Dragon from Lurch
 *
 * Revision 1.1.1.1  2002/09/10 22:45:19  sgicomp
 * Create the first version of the Dragon
 *
 * Revision 1.3  1996/07/01 15:30:45  markus
 * inserted an #include for SUNOS 4.XX
 * changed #include "Arc.h" to #include "XmGraph/Arc.h"
 * added a "_Rubber" at the appropriate place of typedef struct...{}Rubber;
 *
 * Revision 1.2  1996/06/27 18:18:32  bob
 * Inserted typedef struct _Rubber.
 * Changed #include "XmGraph/Arc.h".
 * #ifndef ARCP_H...#define ARCP_H put at top of file.
 *
 * Revision 1.1  1996/02/20 09:57:21  martin
 * Initial revision
 *
 */

/* CYY modified 2/1/91 for  R3 and R4 compatible and ANSI_C fixes*/

#ifdef SUNOS4
#include <X11/IntrinsicP.h>
#endif

#include <Xm/Xm.h>

/* CYY added */
#include "XmGraph/Arc.h"

    
typedef struct arclist {
	XmArcWidgetList  arcs;
	int         n_arcs;
	int         n_slots;
    } ArcList;

typedef struct _XmArcClassPart
{
    XtWidgetProc         input_dispatch;
} XmArcClassPart;


typedef struct _XmArcClassRec {
    CoreClassPart	core_class;
    XmArcClassPart	arc_class;
} XmArcClassRec;

externalref XmArcClassRec xmArcClassRec;

typedef struct _XmArcPart {
    int		from_x,             /* cache the last x,y coords of the */
    		from_y,             /* head (to) and tail (from) of this */
    		to_x,               /* arc widget. Used to efficiently rubber */
    		to_y,               /* band arcs to moving nodes */
    		width,
    		dashes,
    		dash_offset;
    Widget	to;
    Widget	from;
    unsigned char direction,
    		  style,
    		  cap_style;
    Pixel      foreground;
    Pixel      highlightcolor;
    Boolean      highlight;
    Boolean      keephidden;
    Boolean      visible;
    XmFontList	font;  
    _XmString     label;
    Dimension    labelwidth, labelheight;
    Boolean      map_name;
    ArcList      *siblings;
    GC           gc;
    GC           highlight_gc;
    GC           current_gc;
    GC           clear_gc;
    int          delta;
    Region       region;
    XtCallbackList arm_callback;
    XtCallbackList activate_callback;
    XtCallbackList disarm_callback;
    XtCallbackList arc_edited;
    XtPointer        user_data;
    Boolean        armed;
    Boolean        up_to_date;
    int            rank;
} XmArcPart;


typedef struct _XmArcRec {
    CorePart	core;
    XmArcPart	arc;
} XmArcRec;

typedef struct _Rubber {
    Widget child;
    int isarc;
} Rubber;



#endif /* ARCP_H */
