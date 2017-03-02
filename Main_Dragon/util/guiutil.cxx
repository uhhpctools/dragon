#include "guiutil.h"
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
int vecToArrString(lines_t arr, int maxlen, String *&arrstr) {
   int l = arr.size();
   arrstr = new String[l];
   for(int i=0;i<l;i++) {
      arrstr[i] = arr[i];

   }
   return l;
}
