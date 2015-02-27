#ifndef SimpleFits_SoftConstraint_h
#define SimpleFits_SoftConstraint_h

#include  "SimpleFits/FitSoftware/interface/Constraint.h"
#include <string.h>
#include "TMatrixT.h"
#include "TMatrixTSym.h"

class SoftConstraint : public Constraint{
 public:
  SoftConstraint(std::string _name, unsigned int _idx, TMatrixT<double> *_val, TMatrixTSym<double> *_cov, TMatrixT<double> *_dalpha,double _value, double _initalSigma, double _sigma_min, unsigned int _maxIterations=5, double _scalefactor=0.75);
  virtual ~SoftConstraint();

  virtual void addCorrelation(int idx, double corr);

  // important function for modifying the cov to iteratively increase the constriants
  virtual void UpdateCovariance();
  virtual bool isConverged();

  // flags for defining behaviour of fit
  virtual bool isSoft(){return true;}
  virtual bool isHard(){return false;}

};
#endif

