#ifndef ErrorMatrixPropagator_h
#define ErrorMatrixPropagator_h

#include "TMatrixT.h"
#include "TMatrixTSym.h"

class  ErrorMatrixPropagator {
 public:
  ErrorMatrixPropagator(){};
  virtual ~ErrorMatrixPropagator(){};
  static TMatrixTSym<double> PropagateError(TMatrixT<double> (*function)(TMatrixT<double> &par),
                                            TMatrixT<double> inPar,TMatrixTSym<double> inCov,
                                            double epsilon=0.001, double errorEpsilonRatio=1000);
};
#endif


