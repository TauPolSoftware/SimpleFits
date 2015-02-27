#ifndef ConstraintFitterInterator_h
#define ConstraintFitterInterator_h

#include "SimpleFits/FitSoftware/interface/Constraint.h"
#include <vector>
#include <string.h>
#include "TH1D.h"

class  ConstraintFitterInterator{
 public:
  ConstraintFitterInterator(){};
  virtual ~ConstraintFitterInterator(){};

  virtual bool FitwithConstraint();
  virtual void   AddConstraints(Constraint *_constraint);
  virtual void   SaveConstraintConvergence(std::string _file);
  inline unsigned int maxIterations(){return maxiterations_;} 
  virtual double Chi2()=0;
  virtual double NDF()=0;

 protected:
  virtual bool Fit()=0;
  
 private:
  unsigned int maxiterations_;
  std::vector<Constraint*> constraints_;
  std::vector<TH1D>       ConstraintConvergence_;
  
};
#endif


