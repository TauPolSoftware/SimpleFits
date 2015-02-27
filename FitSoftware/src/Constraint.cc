#include  "SimpleFits/FitSoftware/interface/Constraint.h"

Constraint::Constraint(std::string _name, unsigned int _idx,TMatrixT<double> *_val,TMatrixTSym<double> *_cov,TMatrixT<double> *_dalpha, double _initalSigma, double _sigma_min, double _dx, unsigned int _maxIterations,double _scalefactor):
  name_(_name),
  idx_(_idx),
  val_(_val),
  cov_(_cov),
  dalpha_(_dalpha),
  initalSigma_(_initalSigma),
  sigma_min_(_sigma_min),
  dx_(_dx),
  nIterations_(0),
  maxIterations_(_maxIterations),
  scalefactor_(_scalefactor)
{
  val_(idx_,0)=value;
  cov_(idx_,idx_)=initalSigma_;
}

Constraint::~Constraint(){}

bool Constraint::isConverged(){
  if(dx()>dalpha() || sigma_min()*sigma_min()==((*cov_)(idx_,idx_)))return true; 
  return false;
}

void Constraint::UpdateCovariance(){
  if(maxIterations_>nIterations_){
    nIterations_++;
    if(sigma_min()*sigma_min()>((*cov_)(idx_,idx_))*scalefactor_*scalefactor_){
      // prevent sqrt(cov) from being below the requested minimum sigma
      scalefactor_=sigma_min()*sigma_min()/((*cov_)(idx_,idx_));
      ((*cov_)(idx_,idx_))=sigma_min();
    }
    else{
      // else use regular scaling
      ((*cov_)(idx_,idx_))*=scalefactor_*scalefactor_;
    }
    for(unsigned int i=0;i<corridx_.size();i++) ((*cov_)(idx_,corridx_.at(i)));
  }
}
