#ifndef __GUI_UTIL__
#define __GUI_UTIL__

#include <Xm/Xm.h>
#include "util.h"
/**
 * vecToArrString: converting a vector of string into array of string
 * input:
 *   arr: vector of string
 *   maxlen: maximum length of string
 * output:
 *   arrStr: array of string
 * return:
 *   the number of elements or lines
 */ 
int vecToArrString(lines_t arr, int maxlen, String *&arrstr); 
#endif
