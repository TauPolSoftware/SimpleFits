#include "SimpleFits/FitSoftware/interface/ErrorMatrixPropagator.h"
#include "math.h"

TMatrixTSym<double> ErrorMatrixPropagator::PropogateError(TMatrixT<double> (*f)(TMatrixT<double> &par),TMatrixT<double> inPar,TMatrixTSym<double> inCov, double epsilon, double errorEpsilonRatio){
  TMatrixT<double> v=f(inPar);
  TMatrixT<double> Jacobian(inPar.GetNrows(),v.GetNrows());
  for(int i=0;i<inPar.GetNrows();i++){
    TMatrixT<double> ParPlusEpsilon=inPar;
    double error=sqrt(fabs(inCov(i,i)));
    double delta=epsilon;
    if(delta*errorEpsilonRatio<error) delta=error/errorEpsilonRatio;
    ParPlusEpsilon(i,0)+=delta;
    TMatrixT<double> vp=f(ParPlusEpsilon);
    for(int j=0;j<inPar.GetNrows();j++){Jacobian(i,j)=(vp(j,0)-v(j,0))/delta;}// Newtons approx.                                                                                                                                                
  }
  return inCov.SimilarityT(Jacobian);
}
