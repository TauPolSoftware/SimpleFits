#ifndef Particle_h
#define Particle_h

#include "TMatrixT.h"
#include "TMatrixTSym.h"

// Notes
// Store B in units of 1/GeV

class Particle {
public:
  Particle(TMatrixT<double> par_, TMatrixTSym<double> cov_, int pdgid_, double charge_, double b_);
  virtual ~Particle(){};

  virtual double Parameter(int i) const {
    if (0<=i && i<par.GetNrows()) {
      return par(i, 0);
    }
    return 0;
  }

  virtual double Covariance(int i, int j) {
    if (0<=i && i<cov.GetNrows() && 0<=j && j<cov.GetNrows()) {
      return cov(i, j);
    }
    return 0;
  }

  double BField() const { return b; }
  int PDGID() const { return pdgid; }
  double Charge() const { return charge; }
  double qB() const { return b * charge; }

  virtual int NParameters() = 0;

  TMatrixTSym<double> getCovMatrix() const {
	  return cov;
  }

  TMatrixT<double> getParMatrix() const {
	  return par;
  }

private:
  TMatrixT<double> par;
  TMatrixTSym<double> cov;
  int pdgid;
  double charge;
  double b;

};
#endif

