#ifndef SimpleFits_Log_h
#define SimpleFits_Log_h

#include <string.h>
#include <cstdio>
#include <iostream> 
#include <ostream>

class Log {
 public:
  enum level{Fatal=0,Error=1,Warning=2,Verbose=3,Debug=4};
  static Log* Instance(){if(instance==NULL) instance=new Log(); return instance;}

  // set output stream 
  void Set_cout(){s=&std::cout;}
  void Set_cerr(){s=&std::cerr;}
  void Set_stream(std::ostream *stream){s=stream;}

  // Manipulate output levels
  void SetLevel(level _l){l=_l;}
  level Level(){return l;}
  std::ostream& Stream(){return (*s);} 

 private:
  Log():l(Verbose){Set_cout();}
  virtual ~Log(){};
  
  static Log *instance;
  level l;
  std::ostream *s;
};

#define Log(level) \
  if(Log::Instance()->Level()<=level) \
    Log::Instance()->Stream() << "#level: " << __FILE__ << " " <<  __func__ << " " << __LINE__ << " "  

#endif

