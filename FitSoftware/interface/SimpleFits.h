#ifndef SimpleFits_SimpleFits_h
#define SimpleFits_SimpleFits_h

#include "TMatrixT.h"
#include "TMatrixTSym.h"

class  SimpleFits{
 public:
  static SimpleFits* Get(){if(instance==NULL) instance=new SimpleFits(); return instance;}
  
  void Set_epsilonRatio(double _epsilonRatio){epsilonRatio_=_epsilonRatio;}
  double epsilonRatio(){return epsilonRatio_;}
  
  template <class T> 
    TMatrixTSym<double> PropogateError(T *t, TMatrixT<double> (T::*f)(TMatrixT<double> &par),TMatrixT<double> inPar,TMatrixTSym<double> inCov){
    TMatrixT<double> v=((t)->*(f))(inPar);
    TMatrixT<double> Jacobian(inPar.GetNrows(),v.GetNrows());
    for(int i=0;i<inPar.GetNrows();i++){
      TMatrixT<double> ParPlusEpsilon=inPar;
      double error=sqrt(fabs(inCov(i,i)));
      double delta=epsilonRatio_*error;
      ParPlusEpsilon(i,0)+=delta;
      TMatrixT<double> vp=((t)->*(f))(ParPlusEpsilon);
      for(int j=0;j<v.GetNrows();j++){Jacobian(i,j)=(vp(j,0)-v(j,0))/delta;}// Newtons approx.
    }
    TMatrixTSym<double> newCov=inCov.SimilarityT(Jacobian);
    return newCov;
  }
  
 private:
 SimpleFits():epsilonRatio_(0.001){};
  virtual ~SimpleFits(){};

  static SimpleFits *instance;
  double epsilonRatio_;
};
#endif


