#ifndef __APOGEN__
#define __APOGEN__

#include <map>
#include <string>
#include <vector>

typedef struct {
  int iLineNoStart;
  int iLineNoEnd;
  char cStatus;
  string detailstatus;
  string info;
}CodeStatus;

typedef vector<CodeStatus> vecCodeStatus;

typedef struct {
  string filename;
  string PUname;
  vecCodeStatus codeStatus;
} ParStatus;

typedef vector<ParStatus> VecParStatus;
typedef map<int,char> MapLineCodeStatus;

void generateOpenMPCode(bool bApolist, vector<string> &vecSourceCode,
		                                 MapLineCodeStatus codestatus);

void apogen(); 
bool listParStatus(VecParStatus &vecparstatus, MapLineCodeStatus &maplinestatus);
#endif
