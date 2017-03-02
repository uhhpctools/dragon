/* Generated by wbuild from "Frame.w"
** (generator version $Revision: 1.1.1.1 $ of $Date: 2004/08/09 20:03:03 $)
*/
#ifndef _XfwfFrame_H_
#define _XfwfFrame_H_
#include <Xfwf/Common.h>
typedef enum {
	    XfwfRaised, XfwfSunken, XfwfChiseled, XfwfLedged } FrameType;

typedef enum {XfwfAuto, XfwfColor, XfwfStipple} ShadowScheme;

typedef Pixmap  Bitmap;

void XfwfDrawFrame(
#if NeedFunctionPrototypes
Widget,int ,int ,int ,int ,FrameType ,int ,GC ,GC 
#endif
);
Boolean  cvtStringToFrameType(
#if NeedFunctionPrototypes
Display *,XrmValuePtr ,Cardinal *,XrmValuePtr ,XrmValuePtr ,XtPointer *
#endif
);
Boolean  cvtFrameTypeToString(
#if NeedFunctionPrototypes
Display *,XrmValuePtr ,Cardinal *,XrmValuePtr ,XrmValuePtr ,XtPointer *
#endif
);
Boolean  cvtStringToShadowScheme(
#if NeedFunctionPrototypes
Display *,XrmValuePtr ,Cardinal *,XrmValuePtr ,XrmValuePtr ,XtPointer *
#endif
);
Boolean  cvtShadowSchemeToString(
#if NeedFunctionPrototypes
Display *,XrmValuePtr ,Cardinal *,XrmValuePtr ,XrmValuePtr ,XtPointer *
#endif
);
#ifndef XtNcursor
#define XtNcursor "cursor"
#endif
#ifndef XtCCursor
#define XtCCursor "Cursor"
#endif
#ifndef XtRCursor
#define XtRCursor "Cursor"
#endif

#ifndef XtNframeType
#define XtNframeType "frameType"
#endif
#ifndef XtCFrameType
#define XtCFrameType "FrameType"
#endif
#ifndef XtRFrameType
#define XtRFrameType "FrameType"
#endif

#ifndef XtNframeWidth
#define XtNframeWidth "frameWidth"
#endif
#ifndef XtCFrameWidth
#define XtCFrameWidth "FrameWidth"
#endif
#ifndef XtRDimension
#define XtRDimension "Dimension"
#endif

#ifndef XtNouterOffset
#define XtNouterOffset "outerOffset"
#endif
#ifndef XtCOuterOffset
#define XtCOuterOffset "OuterOffset"
#endif
#ifndef XtRDimension
#define XtRDimension "Dimension"
#endif

#ifndef XtNinnerOffset
#define XtNinnerOffset "innerOffset"
#endif
#ifndef XtCInnerOffset
#define XtCInnerOffset "InnerOffset"
#endif
#ifndef XtRDimension
#define XtRDimension "Dimension"
#endif

#ifndef XtNshadowScheme
#define XtNshadowScheme "shadowScheme"
#endif
#ifndef XtCShadowScheme
#define XtCShadowScheme "ShadowScheme"
#endif
#ifndef XtRShadowScheme
#define XtRShadowScheme "ShadowScheme"
#endif

#ifndef XtNtopShadowColor
#define XtNtopShadowColor "topShadowColor"
#endif
#ifndef XtCTopShadowColor
#define XtCTopShadowColor "TopShadowColor"
#endif
#ifndef XtRPixel
#define XtRPixel "Pixel"
#endif

#ifndef XtNbottomShadowColor
#define XtNbottomShadowColor "bottomShadowColor"
#endif
#ifndef XtCBottomShadowColor
#define XtCBottomShadowColor "BottomShadowColor"
#endif
#ifndef XtRPixel
#define XtRPixel "Pixel"
#endif

#ifndef XtNtopShadowStipple
#define XtNtopShadowStipple "topShadowStipple"
#endif
#ifndef XtCTopShadowStipple
#define XtCTopShadowStipple "TopShadowStipple"
#endif
#ifndef XtRBitmap
#define XtRBitmap "Bitmap"
#endif

#ifndef XtNbottomShadowStipple
#define XtNbottomShadowStipple "bottomShadowStipple"
#endif
#ifndef XtCBottomShadowStipple
#define XtCBottomShadowStipple "BottomShadowStipple"
#endif
#ifndef XtRBitmap
#define XtRBitmap "Bitmap"
#endif

typedef struct _XfwfFrameClassRec *XfwfFrameWidgetClass;
typedef struct _XfwfFrameRec *XfwfFrameWidget;
externalref WidgetClass xfwfFrameWidgetClass;
#endif /*_XfwfFrame_H_*/
