#include "TauPolSoftware/SimpleFits/interface/Particle.h"

Particle::Particle(TMatrixT<double> par_, TMatrixTSym<double> cov_, int pdgid_, double charge_, double b_):
  par(par_),
  cov(cov_),
  pdgid(pdgid_),
  charge(charge_),
  b(b_)
{

}
