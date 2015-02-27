#ifndef SimpleFits_HardConstraint_h
#define SimpleFits_HardConstraint_h

#include  "SimpleFits/FitSoftware/interface/Constraint.h"
#include <string.h>
#include "TMatrixT.h"
#include "TMatrixTSym.h"

class HardConstraint : public Constraint{
 public:
  HardConstraint(std::string _name, unsigned int _idx, TMatrixT<double> *_val, TMatrixTSym<double> *_cov, TMatrixT<double> *_dalpha, double _initalSigma, double _dx, unsigned int _maxIterations=5, double _scalefactor=0.75);
  virtual ~HardConstraint();

  virtual bool isSoft(){return false;}
  virtual bool isHard(){return true;}

  // important function for modifying the cov to iteratively increase the constriants
  virtual void UpdateCovariance();
  virtual bool isConverged();

};
#endif

