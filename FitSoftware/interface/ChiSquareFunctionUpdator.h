#ifndef ChiSquareFunctionUpdator_h
#define ChiSquareFunctionUpdator_h

#include "Minuit2/FCNBase.h"
#include "TMatrixT.h"

template <class T>
class ChiSquareFunctionUpdator  : public ROOT::Minuit2::FCNBase {
 public:
  ChiSquareFunctionUpdator(T *Fitter){Fitter_=Fitter;}
  virtual ~ChiSquareFunctionUpdator(){};
  
  virtual double operator() (const std::vector<double> & x)const{
    TMatrixT<double> X(x.size(),1);
    for(unsigned int i=0; i<x.size();i++){X(i,0)=x.at(i);}
    return Fitter_->UpdateChisquare(X);
  }
  virtual double Up()const{return 1.0;}// Error definiton for Chi^2

 private:
  T *Fitter_;
  
};
#endif

