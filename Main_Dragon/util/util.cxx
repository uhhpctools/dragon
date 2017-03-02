#include <algorithm>
#include <vector>
#include <iterator>
#include <stdio.h>

template <class T> 
T * to_array (vector<T> & from)
{
  T * ret = new T [from.size()];
  copy (from.begin(), from.end(), ret);
  return ret;
}

#ifdef __TEST_UTIL__
typedef char* String;

int main() {
  vector <String> from;
  vector <int> de;
  int i;
  char buf[20];
  for(i=0;i<10;i++) {
     sprintf(buf,"%d ",i);
     from.push_back(buf);
     cout <<i<<"="<<buf<<endl;
     de.push_back(i);
  }
  String *ptr = to_array<String>(from);
  for(i=0;i<10;i++)
    cout <<from[i]<< " vs.  " <<ptr[i]<<endl;
  delete ptr;
}
#endif
