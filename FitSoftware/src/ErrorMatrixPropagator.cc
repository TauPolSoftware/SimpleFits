#include "SimpleFits/FitSoftware/interface/ErrorMatrixPropagator.h"
#include "math.h"
#include <iostream>

TMatrixTSym<double> ErrorMatrixPropagator::PropogateError(TMatrixT<double> (*f)(TMatrixT<double> &par),TMatrixT<double> inPar,TMatrixTSym<double> inCov, double epsilon, double errorEpsilonRatio){
  std::cout << "ErrorMatrixPropagator::PropogateError" << std::endl;
  TMatrixT<double> v=f(inPar);
  TMatrixT<double> Jacobian(inPar.GetNrows(),v.GetNrows());
  std::cout << "ErrorMatrixPropagator::PropogateError A" << std::endl;
  for(int i=0;i<inPar.GetNrows();i++){
    TMatrixT<double> ParPlusEpsilon=inPar;
    double error=sqrt(fabs(inCov(i,i)));
    double delta=epsilon;
    if(delta*errorEpsilonRatio<error) delta=error/errorEpsilonRatio;
    ParPlusEpsilon(i,0)+=delta;
    TMatrixT<double> vp=f(ParPlusEpsilon);
    for(int j=0;j<v.GetNrows();j++){Jacobian(i,j)=(vp(j,0)-v(j,0))/delta;}// Newtons approx.
  }
  std::cout << "ErrorMatrixPropagator::PropogateError almost ended" << std::endl;
  return inCov.SimilarityT(Jacobian);
}
