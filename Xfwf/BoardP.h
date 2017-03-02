/* Generated by wbuild from "Board.w"
** (generator version $Revision: 1.1.1.1 $ of $Date: 2004/08/09 20:03:03 $)
*/
#ifndef _XfwfBoardP_H_
#define _XfwfBoardP_H_
#include <Xfwf/FrameP.h>
#include <Xfwf/Board.h>
#define MAGICNUM 12349 


typedef void (*set_abs_location_Proc)(
#if NeedFunctionPrototypes
Widget,unsigned  int ,int ,int ,int ,int 
#endif
);
#define XtInherit_set_abs_location ((set_abs_location_Proc) _XtInherit)

typedef struct {
/* methods */
set_abs_location_Proc set_abs_location;
/* class variables */
} XfwfBoardClassPart;

typedef struct _XfwfBoardClassRec {
CoreClassPart core_class;
CompositeClassPart composite_class;
XfwfCommonClassPart xfwfCommon_class;
XfwfFrameClassPart xfwfFrame_class;
XfwfBoardClassPart xfwfBoard_class;
} XfwfBoardClassRec;

typedef struct {
/* resources */
Position  abs_x;
float  rel_x;
Position  abs_y;
float  rel_y;
Position  abs_width;
float  rel_width;
Position  abs_height;
float  rel_height;
float  hunit;
float  vunit;
String  location;
/* private state */
} XfwfBoardPart;

typedef struct _XfwfBoardRec {
CorePart core;
CompositePart composite;
XfwfCommonPart xfwfCommon;
XfwfFramePart xfwfFrame;
XfwfBoardPart xfwfBoard;
} XfwfBoardRec;

externalref XfwfBoardClassRec xfwfBoardClassRec;

#endif /* _XfwfBoardP_H_ */
