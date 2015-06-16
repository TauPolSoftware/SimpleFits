#ifndef SimpleFits_Logger_h
#define SimpleFits_Logger_h

#include <string.h>
#include <cstdio>
#include <iostream> 
#include <ostream>

class Logger {
 public:
   enum level{Fatal=0,Error=1,Warning=2,Info=3,Verbose=4,Debug=5};
  static Logger* Instance(){if(instance==NULL) instance=new Logger(); return instance;}

  // set output stream 
  void Set_cout(){s=&std::cout;}
  void Set_cerr(){s=&std::cerr;}
  void Set_stream(std::ostream *stream){s=stream;}

  // Manipulate output levels
  void SetLevel(level _l){l=_l;}
  level Level(){return l;}
  std::ostream& Stream(){return (*s);}

  static int levelColor(level l){
	  if (l == Fatal)	return 41; // red background
	  if (l == Error)	return 43; // yellow background
	  if (l == Warning)	return 31; // red font
	  if (l == Info)	return 34; // blue font
	  if (l == Verbose)	return 0; // nothing
	  if (l == Debug)	return 0; // nothing
	  return 30;
  }

 private:
  Logger():l(Verbose){Set_cout();}
  virtual ~Logger(){};
  
  static Logger *instance;
  level l;
  std::ostream *s;
};

#define Logger(level) \
  if(Logger::Instance()->Level()>=level) \
    Logger::Instance()->Stream() << "\033[1;" << Logger::levelColor(level) << "m" << #level << "\033[0m" << "[" << __FILE__ << " " <<  __func__ << "(..) l. " << __LINE__ << "] - "

#endif

