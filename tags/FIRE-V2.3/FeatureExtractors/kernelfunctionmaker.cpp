#include <string>
#include <sstream>
//#include "invariantfeaturehistogram.hpp"
#include "kernelfunction.hpp"
#include "monomial1kernelfunction.hpp"
#include "monomial2kernelfunction.hpp"
#include "monomial3kernelfunction.hpp"
#include "relationalkernelfunction.hpp"

using namespace std;

// helping function to retrieve the parameters of the kernelfunction
int getIntAfter(const string &pattern,const string& s) {
  uint pos=s.find(pattern);
  string nrstring="";
  for(uint i=pos+pattern.size();(s[i]>='0' && s[i]<='9')|| s[i]=='-';++i) {
    nrstring+=s[i];
  }
  istringstream istr(nrstring);
  
  int result=0;
  istr>>result;
  
  DBG(1000) << "after '" << pattern << "' there is a " << result << "='" << nrstring <<"' in '" << s <<"'." << endl;
  return result;
}


// creates a kernelfunction 
KernelFunction* getNewKernelFunction(string description) {
	
  KernelFunction* result=NULL;
	
  if(description.find("mon2:")==0) {
    
    int x1=getIntAfter("x1=",description);
    int y1=getIntAfter("y1=",description);
    int x2=getIntAfter("x2=",description);
    int y2=getIntAfter("y2=",description);
    DBG(25) << "Creating Monomial2KernelFunction" << endl;
    result=new Monomial2KernelFunction(x1,y1,x2,y2);
	  
  } else if(description.find("mon3:")==0) {

    int x1=getIntAfter("x1=",description);
    int y1=getIntAfter("y1=",description);
    int x2=getIntAfter("x2=",description);
    int y2=getIntAfter("y2=",description);
    int x3=getIntAfter("x3=",description);
    int y3=getIntAfter("y3=",description);
    DBG(25) << "Creating Monomial3KernelFunction" << endl;
    result=new Monomial3KernelFunction(x1,y1,x2,y2,x3,y3);

  } else if(description.find("mon1:")==0) {

    int x1=getIntAfter("x1=",description);
    int y1=getIntAfter("y1=",description);
    DBG(25) << "Creating Monomial1KernelFunction" << endl;
    result=new Monomial1KernelFunction(x1,y1);
    
  } else if(description.find("rel:")==0) {

    int x1=getIntAfter("x1=",description);
    int y1=getIntAfter("y1=",description);
    int x2=getIntAfter("x2=",description);
    int y2=getIntAfter("y2=",description);
    DBG(25) << "Creating RelationalKernelFunction" << endl;
    result=new RelationalKernelFunction(x1,y1,x2,y2);
    
  } else {
	  
    ERR << "No valid kernel function given" << endl;
	  
  }
  
  return result;
}
