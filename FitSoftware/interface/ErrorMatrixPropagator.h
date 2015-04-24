#ifndef ErrorMatrixPropagator_h
#define ErrorMatrixPropagator_h

#include "TMatrixT.h"
#include "TMatrixTSym.h"

namespace  ErrorMatrixPropagator {
  TMatrixTSym<double> PropagateError(void* ptr2Object, TMatrixT<double> (*ptr2Function)(void* ptr2Object, TMatrixT<double> &par),TMatrixT<double> inPar,TMatrixTSym<double> inCov, double epsilon=0.001, double errorEpsilonRatio=1000);
};
#endif


