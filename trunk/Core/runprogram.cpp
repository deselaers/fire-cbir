#include <iostream>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <sstream>
#include <errno.h>
#include "diag.hpp"
#include "runprogram.hpp"


using namespace std;

// this is much cooler than just calling "system" or "popen" :-)
// although calling system or popen would have been much easier...
// anyway, if I had thought about this, this file would be much shorter ;-)


void runInBackground(string commandline) {
  pid_t pid=vfork();
  if(pid==0) {
    DBG(20) << "Forking child..." << endl;
    run(commandline);
  } else if(pid>0) {
    // here is the parent process, nothing to do, just return
  } else {
    ERR << "Error forking child:" << strerror(errno) << endl;
  }
}

void run(string commandline) {
  DBG(20) << "Executing " << commandline << endl;
  istringstream iss(commandline);
  string command, tmp;
  iss >> command;
  char **opts=(char**)malloc(1*sizeof(char*));
  opts[0]=(char*)malloc(sizeof(char)*(command.size()+1));
  strcpy(opts[0],command.c_str());
  int len=1;
  while(!iss.eof()) {
    iss >> tmp;
    opts=(char**)realloc(opts,(len+1)*sizeof(char*));
    opts[len]=(char*)malloc(sizeof(char)*(tmp.size()+1));
    strcpy(opts[len],tmp.c_str());
    ++len;    
  }
  opts=(char**)realloc(opts,(len+1)*sizeof(char*));
  opts[len]=NULL;
  
  DBG(15) << "Running :";
  for(int i=0;i<len;++i) { BLINK(15) << opts[i] << " ";} BLINK(15) << endl;
  
  char *cmd=(char*)malloc(sizeof(char)*(command.size()+1));
  strcpy(cmd,command.c_str());
  int ERROR=execvp(cmd,opts);
  if(ERROR==-1) {
    ERR << "Execution failed: " << errno << ": " << strerror(errno) << endl;
  }
  
  for(int i=0;i<len+1;++i) {
    free(opts[i]);
  } 
  free(opts);
  free(cmd);
}


void waitForAllChildren() {
  DBG(10) << "Waiting for all children to finish... " ;
  pid_t pid;
  pid=1;
  while(pid>0) {
    pid=wait(NULL);
    DBG(20) << "Child " << pid << " terminated." << endl;
  }
  DBG(10) << "done" << endl;
}
