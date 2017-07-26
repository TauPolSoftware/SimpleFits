#include "TauPolSoftware/SimpleFits/interface/ErrorMatrixPropagator.h"
#include "math.h"
#include <iostream>

TMatrixTSym<double> ErrorMatrixPropagator::PropagateError(TMatrixT<double> (*function)(TMatrixT<double> &par),
                                                          TMatrixT<double> inPar, TMatrixTSym<double> inCov,
                                                          double epsilon, double errorEpsilonRatio) {
  TMatrixT<double> outPar = function(inPar);
  TMatrixT<double> jacobian(inPar.GetNrows(), outPar.GetNrows());
  for (int row = 0; row < inPar.GetNrows(); row++) {
    TMatrixT<double> inParPlusEpsilon = inPar;
    double error = std::sqrt(std::fabs(inCov(row, row)));
    double delta = epsilon;
    if ((std::fabs(delta * errorEpsilonRatio) > error) && (error > 0.0)) {
      delta = std::fabs(error / errorEpsilonRatio);
    }
    inParPlusEpsilon(row, 0) += delta;
    TMatrixT<double> outParPlusEpsilon = function(inParPlusEpsilon);
    for(int col = 0; col < outPar.GetNrows(); col++) {
      jacobian(row, col) = (outParPlusEpsilon(col, 0) - outPar(col, 0)) / delta; // Newtons approx.
    }
  }
  TMatrixTSym<double> outCov = inCov.SimilarityT(jacobian);
  /// debug
  /*
  std::cout << "ErrorMatrixPropagator::PropogateError inCov:" << std::endl;
  inCov.Print();
  std::cout << "ErrorMatrixPropagator::PropogateError jacobian:" << std::endl;
  jacobian.Print();
  std::cout << "ErrorMatrixPropagator::PropogateError outCov:" << std::endl;
  outCov.Print();
  */
  return outCov;
}
