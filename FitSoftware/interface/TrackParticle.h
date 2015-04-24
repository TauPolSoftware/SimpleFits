#ifndef TrackParticle_h
#define TrackParticle_h

#include "SimpleFits/FitSoftware/interface/Particle.h"
#include "TString.h"

class TrackParticle : public Particle {
 public:
  enum {kappa=0,lambda,phi,dxy,dz,NHelixPar};// 5 track helix Parameters

  TrackParticle();
  TrackParticle(const TrackParticle& another);
  TrackParticle(TMatrixT<double> par_, TMatrixTSym<double> cov_, int pdgid_, double mass_,double charge_, double b_);
  virtual ~TrackParticle(){};

  static TString Name(int i);
  int NParameters(){return NHelixPar;}

  double Mass() const{return mass;}
  
 private:
  double mass;
};
#endif


