#ifndef SCOPETIMER_H
#define SCOPETIMER_H

#ifdef WIN32
 #include <Windows.h>
 #undef min
 #undef max
 #define Li2Double(x) ((double)((x).HighPart) * 4.294967296E9 + (double)((x).LowPart))
#else
 #include <sys/time.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif
#include <map>



class RealTimerCL
{
private:
   double _t_begin, _t_end;
   double _time;

public:
   RealTimerCL(double time=0.);
   void Reset(double time= 0);
   void Start();
   void Stop();
   double GetTime();
   
   double timestamp();
};



//*******************************************************************
// S c o p e T i m e r C L
//   accumulates elapsed times between creation and 
//   destruction of objects with the same 'name'
//*******************************************************************
class ScopeTimerCL
{
private: 
   char* _name;
   RealTimerCL _timer;

public:
   ScopeTimerCL(char *name);
   ~ScopeTimerCL();
};
typedef ScopeTimerCL ScopeTimer;



#endif
