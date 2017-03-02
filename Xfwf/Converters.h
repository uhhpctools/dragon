#ifndef _Converters_h
#define _Converters_h

#include "Xfwf/long.h"
#include "Xfwf/icon.h"
#include "Xfwf/choosecol.h"
#include "Xfwf/strarray.h"
#include "Xfwf/StrToPmap.h"
#include "Xfwf/Pen.h"

#define done_bert(type, value) \
    do {\
	if (to->addr != NULL) {\
	    if (to->size < sizeof(type)) {\
	        to->size = sizeof(type);\
	        return False;\
	    }\
	    *(type*)(to->addr) = (value);\
        } else {\
	    static type static_val;\
	    static_val = (value);\
	    to->addr = (XtPointer)&static_val;\
        }\
        to->size = sizeof(type);\
        return True;\
    } while (0)

#define done_bob(type, value) \
   {                                    \
       if (toVal->addr != NULL) {       \
          if (toVal->size < sizeof(type)) {\
             toVal->size = sizeof(type);\
             return False;              \
          }                             \
          *(type*)(toVal->addr) = (value);\
       }                                \
       else {                           \
          static type static_val;       \
          static_val = (value);         \
          toVal->addr = (XPointer)&static_val;\
       }                                \
       toVal->size = sizeof(type);      \
   }

#endif /* _Converters_h */
