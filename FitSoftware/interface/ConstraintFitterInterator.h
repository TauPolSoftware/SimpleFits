#ifndef SimpleFits_ConstraintFitterInterator_h
#define SimpleFits_ConstraintFitterInterator_h

#include "SimpleFits/FitSoftware/interface/Constraint.h"
#include "SimpleFits/FitSoftware/interface/FitterBase.h"
#include <vector>
#include <string.h>
#include "TH1D.h"

class  ConstraintFitInterator : public FitterBase {
 public:
  ConstraintFitterInterator()includeConstraint(true),monitor(false){};
  virtual ~ConstraintFitterInterator(){};

  virtual bool Fit();
  inline unsigned int maxIterations(){return maxiterations_;}
  bool Monitor(bool m){monitor=m;}
  virtual void SaveConstraintConvergence(std::string _file);

 protected:
  virtual bool FitwithConstraints()=0; // add internal fit function that the iteration can call 
  virtual void AddConstraint(Constraint* c);

 private:
  bool FitAndMonitor();
  void CorrectChi2();
  unsigned int maxiterations_;
  unsigned int iteration_;
  double uncorrChi2;
  bool monitor_;
  std::vector<Constraint*> constraints_;
  std::vector<TH1D>        ConstraintConvergence_;
  std::vector<TH1D>        ConstraintSignificance_;
  TH1D FitChi2; 
  TH1D FitUncorrectedChi2;
};
#endif


