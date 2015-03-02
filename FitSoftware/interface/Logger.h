#ifndef SimpleFits_Logger_h
#define SimpleFits_Logger_h

#include <string.h>
#include <cstdio>
#include <iostream> 
#include <ostream>

class Logger {
 public:
  enum level{Fatal=0,Error=1,Warning=2,Verbose=3,Debug=4};
  static Logger* Instance(){if(instance==NULL) instance=new Logger(); return instance;}

  // set output stream 
  void Set_cout(){s=&std::cout;}
  void Set_cerr(){s=&std::cerr;}
  void Set_stream(std::ostream *stream){s=stream;}

  // Manipulate output levels
  void SetLevel(level _l){l=_l;}
  level Level(){return l;}
  std::ostream& Stream(){return (*s);} 

 private:
  Logger():l(Verbose){Set_cout();}
  virtual ~Logger(){};
  
  static Logger *instance;
  level l;
  std::ostream *s;
};

#define Logger(level) \
  if(Logger::Instance()->Level()<=level) \
    Logger::Instance()->Stream() << #level << " [File: " << __FILE__ << " Function: " <<  __func__ << " Line: " << __LINE__ << "] - "  

#endif

