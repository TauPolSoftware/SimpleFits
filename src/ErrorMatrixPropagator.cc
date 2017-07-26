#include "TauPolSoftware/SimpleFits/interface/ErrorMatrixPropagator.h"
#include "math.h"
#include <iostream>

TMatrixTSym<double> ErrorMatrixPropagator::PropagateError(TMatrixT<double> (*function)(TMatrixT<double> &par),
                                                          TMatrixT<double> inPar, TMatrixTSym<double> inCov,
                                                          double epsilon, double errorEpsilonRatio) {
  
  // https://de.wikipedia.org/wiki/Fehlerfortpflanzung#Generalisiertes_Fehlerfortpflanzungsgesetz
  // https://de.wikipedia.org/wiki/Jacobi-Matrix
  // https://root.cern.ch/doc/v608/classTMatrixTSym.html#ab26c079dab1c5e71c49934fbe4bd6c21
  TMatrixT<double> outPar = function(inPar);
  TMatrixT<double> jacobian(outPar.GetNrows(), inPar.GetNrows());
  for (int col = 0; col < jacobian.GetNcols(); col++) {
    TMatrixT<double> inParPlusEpsilon = inPar;
    double error = std::sqrt(std::fabs(inCov(col, col)));
    double delta = epsilon;
    if ((std::fabs(delta * errorEpsilonRatio) > error) && (error > 0.0)) {
      delta = std::fabs(error / errorEpsilonRatio);
    }
    inParPlusEpsilon(col, 0) += delta;
    TMatrixT<double> outParPlusEpsilon = function(inParPlusEpsilon);
    for(int row = 0; row < outPar.GetNrows(); row++) {
      jacobian(row, col) = (outParPlusEpsilon(row, 0) - outPar(row, 0)) / delta; // Newtons approx.
    }
  }
  TMatrixTSym<double> outCov = inCov.Similarity(jacobian);
  /// debug
  /*
  std::cout << "ErrorMatrixPropagator::PropogateError inCov = (" << inCov.GetNrows() << " x " << inCov.GetNcols() << ")" << std::endl;
  inCov.Print();
  std::cout << "ErrorMatrixPropagator::PropogateError jacobian = (" << jacobian.GetNrows() << " x " << jacobian.GetNcols() << ")" << std::endl;
  jacobian.Print();
  std::cout << "ErrorMatrixPropagator::PropogateError outCov = jacobian * inCov * jacobian^T = (" << outCov.GetNrows() << " x " << outCov.GetNcols() << ")" << std::endl;
  outCov.Print();
  */
  return outCov;
}
