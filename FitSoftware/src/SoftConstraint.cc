#include  "SimpleFits/FitSoftware/interface/SoftConstraint.h"

SoftConstraint::SoftConstraint(std::string _name, unsigned int _idx,TMatrixT<double> *_val,TMatrixTSym<double> *_cov,TMatrixT<double> *_dalpha,double _value, double _initalSigma, double _sigma_min, unsigned int _maxIterations,double _scalefactor):
  Constraint(_name,_idx,_val,_cov,_dalpha,_value,_initalSigma,_sigma_min,0,_maxIterations,_scalefactor)
{

}

SoftConstraint::~SoftConstraint(){}

bool SoftConstraint::isConverged(){
  if(sigma_min()*sigma_min()==((*cov_)(idx_,idx_)))return true; 
  return false;
}

void SoftConstraint::UpdateCovariance(){
  if(maxIterations_>nIterations_ && !isConverged()){
    nIterations_++;
    if(sigma_min()*sigma_min()>((*cov_)(idx_,idx_))*scalefactor_*scalefactor_){
      // prevent sqrt(cov) from being below the requested minimum sigma
      scalefactor_=sigma_min()*sigma_min()/((*cov_)(idx_,idx_));
      ((*cov_)(idx_,idx_))=sigma_min()*sigma_min();
    }
    else{
      // else use regular scaling
      ((*cov_)(idx_,idx_))*=scalefactor_*scalefactor_;
    }
    for(unsigned int i=0;i<corridx_.size();i++) ((*cov_)(idx_,corridx_.at(i)))*=scalefactor_;
  }
}

void SoftConstraint::addCorrelation(int idx, double corr){
  if(idx<cov_->GetNrows()) corridx_.push_back(idx);
  cov_(idx_,idx)=sqrt(fabs(cov_(idx,idx)))*sqrt(fabs(cov_(idx_,idx_)))*corr;
  cov_(idx,idx_)=cov_(idx_,idx);
}
