#include  "SimpleFits/FitSoftware/interface/HardConstraint.h"

HardConstraint::HardConstraint(std::string _name, unsigned int _idx,TMatrixT<double> *_val,TMatrixTSym<double> *_cov,TMatrixT<double> *_dalpha,double _value, double _initalSigma, double _dx, unsigned int _maxIterations,double _scalefactor):
  Constraint(_name,_idx,_val,_cov,_dalpha,_value,_initalSigma,0,_dx,_maxIterations,_scalefactor)
{
}

HardConstraint::~HardConstraint(){}

bool HardConstraint::isConverged(){
  if(dx()>dalpha())return true; 
  return false;
}

void HardConstraint::UpdateCovariance(){
  if(maxIterations_>nIterations_ && !isConverged()){
    nIterations_++;
    ((*cov_)(idx_,idx_))*=scalefactor_*scalefactor_;
  }
}
