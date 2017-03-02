#ifndef __READFILE_H__
#define __READFILE_H__
#include <string>
#include <vector>

typedef std::vector< std::string > lines_t;

// read_lines: read a file, and copy the content into a vector of strings
// input:
//   filename: the file name
// output:
//   lines: the content of the file (in vector)
//   maxlength: the maximum length of the lines
// return: 
//   the number of lines  
int read_lines(char *filename, lines_t &lines, int &maxlength);
std::string getfilebasename(string filename);
string getfilebasename(string filename, string &ext); 

/**
 * vecToArray: converting a vector of string into array of char*
 * input:
 *   arr: vector of string
 *   maxlen: maximum length of string
 * output:
 *   arrStr: array of chars
 * return:
 *   the number of elements or lines
 */ 
int vecToArray(lines_t arr, int maxlen, char ***&arrStr); 

#endif	
