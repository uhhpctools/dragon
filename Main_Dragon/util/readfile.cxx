/* Modified from: readfile.cpp
 */

#include <stdio.h>
#include <string>
#include <vector>
#include "readfile.h"
/**
 * Data definition
 */ 

	
/**
 * Function implementation
 */ 
// get the basic filename and its extension
string getfilebasename(string filename, string &ext) {
  string basename;
  int pos = filename.rfind(".");
  if(pos>0) {
      basename = filename.substr(0,pos);
      ext = filename.substr(pos+1);
  } else {
      basename = filename;
      ext = "";
  }
  return basename;
}

// get the basic filename
string getfilebasename(string filename) {
  string basename;
  int pos = filename.rfind(".");
  if(pos>0) {
      basename = filename.substr(0,pos);
  } else {
      basename = filename;
  }
  return basename;
}

// read_lines: read a file, and copy the content into a vector of strings
// input:
//   filename: the file name
// output:
//   lines: the content of the file (in vector)
//   maxlength: the maximum length of the lines
// return: 
//   the number of lines  
int read_lines(char *filename, lines_t &lines, int &maxlength)
{
	// verify the file first
	if(!filename){
		cerr <<"Error: Not valid filename "<<endl;
		return -1;
	}

	// open the file
	FILE *fp = fopen(filename, "r");
	if (!fp)
	{
		cerr <<"Error: Failed to open "<<filename<<endl;
		return -2;
	}

	const size_t buffer_size = 1024;
	char *buffer = new char[buffer_size];
	int linesizes = 0;
	std::string residue;

	maxlength = 0;
	// browse the content of the file until end of file or error
	// we use primitive c function because it is faster (theoritically) than
	//  the famous c++ ifstream
	while (!feof(fp) && !ferror(fp))
	{
		size_t bytes_read = fread(buffer, 1, buffer_size, fp);

		const char *buffer_end = buffer + bytes_read;
		
		const char *b = buffer;
		const char *e;
		while((e = (const char *)memchr(b, '\n', buffer_end - b)) != NULL)
		{
			std::string line(b, e);
			residue += line;
			lines.push_back(residue);
			// looking for the maximum length
			if(residue.size()>maxlength)
			   maxlength=residue.size();
			residue.clear();

			b = e + 1;
		}

		// Is there any residue?
		if (b < buffer_end)
		{
			// There's some residue
			residue += std::string(b, buffer_end);
		}
	}

	if (!residue.empty())
		lines.push_back(residue);

	delete [] buffer;

	fclose(fp);

	int n = 0;
	int linessize = lines.size();

	return linessize;
}

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
int vecToArray(lines_t arr, int maxlen, char ***&arrStr) {
   int l = arr.size();
   int i;
   arrStr = new char**[l];
   for(i=0;i<l;i++) {
     arrStr[i] =  new char*[3];
     arrStr[i][0] = new char[6];
     arrStr[i][1] = new char[128];
     arrStr[i][2] = new char[maxlen];
     sprintf(arrStr[i][0],"%d",i+1);
     strcpy(arrStr[i][1],"N/A");
     strcpy(arrStr[i][2], arr[i].c_str());
   }
   return l;
}

#ifdef __TEST__
int main(int argc, char *argv[]) 
{
   char ***arrStr;
   lines_t arr;
   int maxlen;
   int l = read_lines(argv[1],arr, maxlen);
   l = vecToArray(arr, maxlen, arrStr);
   cout <<argv[1]<<" has "<<l<<" lines"<<" with max len="<<maxlen<<endl;
//   arrStr = new char*[l];
   for(int i=0;i<l;i++) {
  //   arrStr[i] = new char[maxlen];
  //   strcpy(arrStr[i], arr[i].c_str());
     cout <<i<<": "<<arrStr[i]<<endl;
   }
   for(int i=0;i<l;i++)
      delete[] arrStr[i];
   delete[] arrStr;
}
#endif
