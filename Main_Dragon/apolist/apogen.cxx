#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include "../FGnode.h"
#include "../UIClasses.h"
#include "Apolist.h"
#include "../util/util.h"
#include "apogen.h"

extern UIGraphNode *UICallgraphRoot;     // main data
extern int ProcedureLoaded;              // to know which one is the current subroutine
extern void PrintMessage(char Message[]); // to print message in main window

/**************************
 * Function implementation
 */
void generateOpenMPCode(bool bApolist, vector<string> &vecSourceCode,
                 MapLineCodeStatus codestatus)
{
        // we generate only if we have information from Apolist flag
  if (bApolist) {
    MapLineCodeStatus::iterator iter;
    char *ptrFilename = UICallgraphRoot[0].GetFilename();
    int ilen = strlen(ptrFilename);
    // check the language
    // if *.c or *.cxx or *.C or *.CXX then it's a c/c++
    bool c   = (ptrFilename[ilen-1]=='c') || (ptrFilename[ilen-1]=='C')
               || (ptrFilename[ilen-3]=='c') || (ptrFilename[ilen-3]=='C');

    // try to find the parall status
    vector<string>::iterator itbegin=vecSourceCode.begin();
    int iPar = 0;
    for(iter=codestatus.begin(); iter!=codestatus.end(); iter++) {
      if(iter->second == 'A') {
        int iLine = (int)iter->first;
        if(c) { // c family
           vecSourceCode.insert(itbegin + iLine + iPar - 1, "#pragma omp parallel do");
        } else { // fortran family
           vecSourceCode.insert(itbegin + iLine + iPar - 1, "!$OMP PARALLEL DO");
        }
        iPar++;
      }
    }
  }
}

/**
 * List all parallel status for all subroutines
 */
bool listParStatus(VecParStatus &vecparstatus, MapLineCodeStatus &maplinestatus) {
//  map<string,int> mapPUname;
  FGRAPH cfg;
  int i=0;
  MapProginfo proginfo;
  for(i=0;i<UIGraphNode::NumOfNodes;i++) {
     ParStatus parStatus;
     string filename = UICallgraphRoot[i].GetFilename();
     string listfile;
     string PU_Name  = UICallgraphRoot[i].GetName();
     filename = getfilebasename(filename);
     listfile = filename + ".list"; // filename for -apolist
     filename = filename + ".cfg"; // filename for callgraph
     //     
     cfg.Load_CFG((char*)filename.c_str(), (char*)PU_Name.c_str());
     parStatus.filename = listfile;
     parStatus.PUname = PU_Name;
     if (loadFile(listfile.c_str(), proginfo) != 0)	// load apolist info file
	     return false;
	/**
	 * Now try to use   it.
	 **/ 
     int iloop = 0; // we need this special variable to be the indice of the looops
                  // since the cfg nodes may contain other nodes than loops (even enddo)
     char msg[2024];
     for(int l=0;l<cfg.GetSizeofFGnodes(); l++)
     {
      	   if(strcmp(cfg.FG_Nodes[l].Get_Name(),"DOSTART")==0 || 
   	        strcmp(cfg.FG_Nodes[l].Get_Name(),"LOOP")==0) {
	        /* here we find the loop node */
	 	strcpy(msg,cfg.FG_Nodes[l].Get_Name()); // get the name of the node
		CodeStatus objStatus;
		objStatus.iLineNoStart = cfg.FG_Nodes[l].line_nums[0];
		objStatus.iLineNoEnd   = cfg.FG_Nodes[l].line_nums[cfg.FG_Nodes[l].GetSizeofLinenums()-1];
	  	if(proginfo.find(parStatus.PUname)!=proginfo.end()){ // the subroutine exist
		  // make sure that the subprogram has loop information
		  // otherwise, a segmentation fault can be encountered
		  if(proginfo[parStatus.PUname].lineinfo.size()>iloop) {
		     objStatus.detailstatus = proginfo[parStatus.PUname].lineinfo[iloop].status.c_str();
		     // find the status code
		     if(objStatus.detailstatus[11]=='M')
			     objStatus.cStatus = 'M'; // manual
		     else if(objStatus.detailstatus[16]=='S')
			     objStatus.cStatus = 'S'; // synchronized
		     else if(objStatus.detailstatus[11]=='A')
			     objStatus.cStatus = 'A'; // automatic
		     else
			     objStatus.cStatus = 'N'; // not a parallel loop

		     objStatus.info = "";
	           // Print the message (if exist)
            	     for(int k=0;k<proginfo[parStatus.PUname].lineinfo[iloop].info.size();k++) {
			objStatus.info += proginfo[parStatus.PUname].lineinfo[iloop].info[k];
            	     } //for
		  } // if
	  	} // if
		parStatus.codeStatus.push_back(objStatus);
		maplinestatus[objStatus.iLineNoStart] = objStatus.cStatus;
	     	iloop++;
	     } // if
     } //for
     vecparstatus.push_back(parStatus);	 
  }
  return(i>0);
}

/**
 * Generate information based on apolist flag
 */ 
void apogen() {
  char *PU_Name, *FileName, *strlistname;
  int len, n;

	/*--------------------  get the file name first */
  len=strlen(UICallgraphRoot[ProcedureLoaded].GetName());
  PU_Name = new char[len+1];
  strcpy(PU_Name, UICallgraphRoot[ProcedureLoaded].GetName());
  len = strlen(UICallgraphRoot[ProcedureLoaded].GetFilename());
  FileName = new char[len+5];
  strlistname = new char[len+6];
  strcpy(FileName,UICallgraphRoot[ProcedureLoaded].GetFilename());
  for(n=len-1; n>0; n--)
	if(FileName[n]!='.')
		FileName[n]='\0';
	else
		break;
   if(n>0) {
	sprintf(strlistname,"%slist",FileName);
       	strcat(FileName,"cfg");
   } else{
	strcpy(FileName,UICallgraphRoot[ProcedureLoaded].GetFilename());
	sprintf(strlistname,"%s.list",FileName);
	strcat(FileName,".cfg");
   }

	/*-------------------- CFG creation */
   FGRAPH cfg;
   cfg.Load_CFG(FileName,PU_Name);	// load CFG for a file
   MapProginfo proginfo;
   loadFile(strlistname, proginfo);	// load apolist info file

	/**
	 * Now try to use   it.
	 **/ 
   int iloop = 0; // we need this special variable to be the indice of the looops
                  // since the cfg nodes may contain other nodes than loops (even enddo)
   char msg[1024];
   sprintf(msg,"Subprogram %s", PU_Name);
   PrintMessage(msg);
   for(int i=0;i<cfg.GetSizeofFGnodes(); i++)
   {
      if(strcmp(cfg.FG_Nodes[i].Get_Name(),"DOSTART")==0 || 
        strcmp(cfg.FG_Nodes[i].Get_Name(),"LOOP")==0) {
	  /* here we find the loop node */
	  strcpy(msg,cfg.FG_Nodes[i].Get_Name());
	  for (int j=0;j<cfg.FG_Nodes[i].GetSizeofLinenums();j++) {
	     sprintf(msg,"%s line %d",msg,cfg.FG_Nodes[i].line_nums[j]); // print the line numbers
	  }
	  string str=PU_Name;
	  if(proginfo.find(str)!=proginfo.end()){ // the subroutine exist
	    char *ptrMsg = new char[proginfo[str].lineinfo[iloop].status.size()+2];
	    strcpy(ptrMsg,proginfo[str].lineinfo[iloop].status.c_str());
	    sprintf(msg,"%s %s",msg,ptrMsg);
	    PrintMessage(msg); // print it into the main window
	    delete ptrMsg;
	    // Print the message (if exist)
            for(int k=0;k<proginfo[str].lineinfo[iloop].info.size();k++) {
              sprintf(msg,"\t%s",proginfo[str].lineinfo[iloop].info[k].c_str());
	      PrintMessage(msg);
            }
	  }
	  iloop++;
      }
    }
    delete FileName;
    delete strlistname;
}


#ifdef __TEST_APOGEN__
int main(int argc, char *argv[])
{
   if(argc>0) {
     VecParStatus parstatus;
     if(listParStatus(parstatus)) {
	for(int i=0;i<parstatus.size();i++) {
	   cout <<"Filename:"<<parstatus[i].filename<<endl;
	   cout <<"PU  name:"<<parstatus[i].PUname<<endl;
	   for(int j=0;j<parstatus[i].codeStatus.size();j++) {
	     cout <<"line  :"<<parstatus[i].codeStatus[j].iLineNoStart<<endl;
	     cout <<"status:"<<parstatus[i].codeStatus[j].detailstatus<<endl;
	     cout <<"info  :"<<parstatus[i].codeStatus[j].info<<endl;
	   }
	}
     }
   }
}
#endif
