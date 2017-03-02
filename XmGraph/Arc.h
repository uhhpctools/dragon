/*
 * $Id: Arc.h,v 1.1.1.1 2004/08/09 20:03:03 liao Exp $
 *
 * $Log: Arc.h,v $
 * Revision 1.1.1.1  2004/08/09 20:03:03  liao
 * Imported Dragon from Lurch
 *
 * Revision 1.1.1.1  2002/09/10 22:45:19  sgicomp
 * Create the first version of the Dragon
 *
 * Revision 1.2  1996/06/27 18:16:59  bob
 * #define XmRLineStyle            "LineStyle" #ifndef'ed for LINUX.
 *
 * Revision 1.1  1996/02/20 09:57:17  martin
 * Initial revision
 *
 */

/* CYY modified 2/1/91 for  R3 and R4 compatible and ANSI_C fixes*/
#ifndef _XmArc_h
#define _XmArc_h

externalref WidgetClass xmArcWidgetClass;
typedef struct _XmArcClassRec *XmArcWidgetClass;
typedef struct _XmArcRec      *XmArcWidget;
typedef XmArcWidget 	      *XmArcWidgetList;

#define  XmARC_BIT    (1<<29)

#define XmIsArc(w)  XtIsSubclass(w, xmArcWidgetClass)


#ifndef XmNhighlightColor
#define XmNhighlightColor	"highlightColor"
#endif

#ifndef XmCHighlightColor
#define XmCHighlightColor	"HighlightColor"
#endif

#define XmNlabel		"label"
#define XmNarcEditedCallback	"arcEditedCallback"
#define XmNaddArcCallback       "addArcCallback"
#define XmNto		        "to"
#define XmCTo	                "To"
#define XmNfrom		        "from"
#define XmCFrom	             	"From"
#define XmNmapLabel             "mapLabel"
#define XmCMapLabel             "MapLabel"
#define XmNdelta		"delta"
#define XmCDelta		"Delta"
#define XmCHighlight		"Highlight"
#define XmNkeephidden           "keephidden"
#define XmCKeephidden           "Keephidden"
#define XmNarcDirection		"arcDirection" /* XmCDirection class */
#define XmNstyle                "style"
#define XmCStyle                "Style"

/* weng */
/* #ifndef LINUX  */
#define XmRLineStyle            "LineStyle"
/* #endif */

#define XmNcapStyle             "capStyle"
#define XmCCapStyle             "CapStyle"
#define XmRCapStyle             "CapStyle"
#define XmNdashes               "dashes"
#define XmCDashes               "Dashes"
#define XmNdashOffset           "dashOffset"
#define XmCDashOffset           "DashOffset"
#define XmCDirection            "Direction"
#define XmRArcDirection         "ArcDirection"
#define XmNarcWidth             "arcWidth"
#define XmCArcWidth         	"ArcWidth"


#ifdef __cplusplus                    /* do not leave open across includes */
extern "C" {                                  /* for C++ V2.0 */
#endif

#ifndef _NO_PROTO

extern Region  _AddRegionToRegion (Region,Region);   /* Arc.c */
extern Boolean _ArcInRect	(XmArcWidget, XRectangle *); /* Graph.c */
extern void    _EraseArc 	(XmArcWidget arc);  /* Arc.c */
extern void    _SetupArcInternal (XmArcWidget arc);  /* Arc.c */
extern void    _XmUnhighlightArc (XmArcWidget);      /* Arc.c */
extern void    _XmHighlightArc 	(XmArcWidget);      /* Arc.c */
extern int     _sibling_rank 	(XmArcWidget);      /* Arc.c */
extern void    ComputeRegionsForArc (XmArcWidget);   /* Arc.c */
extern void    FreeArcRegions 	(XmArcWidget);      /* Arc.c */
extern void    XmArcKeepHidden  (XmArcWidget, Boolean); /*ME : new in Arc.c */
extern void    XmArcGetPos	(XmArcWidget, Position *, Position *, 
					      Position *, Position *);/*Arc.c*/
extern Widget  XmCreateArc	(Widget,String,ArgList,Cardinal); /* Graph.c */
extern Widget  XmCreateAttachedArc (Widget,String,Widget, Widget, 
					  ArgList,Cardinal); /* Graph.c */



#else /* _NO_PROTO */

extern Widget XmCreateArc();	 	/* Graph.c */
extern Widget XmCreateAttachedArc();	 /* Graph.c */
extern Region _AddRegionToRegion ();     /* Arc.c */
extern Boolean _ArcInRect	();      /* Graph.c */
extern void   _EraseArc 	();      /* Arc.c */
extern void   _XmUnhighlightArc ();      /* Arc.c */
extern void   _XmHighlightArc 	();      /* Arc.c */
extern int    _sibling_rank 	();      /* Arc.c */
extern void   ComputeRegionsForArc ();   /* Arc.c */
extern void   FreeArcRegions 	();      /* Arc.c */
extern void    XmArcKeepHidden  ();      /*ME : new in Arc.c */
#endif  /* _NO_PROTO */


#ifdef __cplusplus
}                                             /* for C++ V2.0 */
#endif


#define XmBIDIRECTED  0
#define XmDIRECTED    1
#define XmUNDIRECTED  2

#define XmLineSolid               0
#define XmLineOnOffDash           1
#define XmLineDoubleDash          2


#define XmCapNotLast              0
#define XmCapButt                 1
#define XmCapRound                2
#define XmCapProjecting           3

#endif
