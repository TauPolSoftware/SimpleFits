#ifndef TrackTools_h
#define TrackTools_h

#include "TMatrixT.h"
#include "TMatrixTSym.h"
#include "TVector3.h"
#include "TauPolSoftware/SimpleFits/interface/TrackParticle.h"
#include "TauPolSoftware/SimpleFits/interface/LorentzVectorParticle.h"

namespace  TrackTools {
  TVector3 PropogateToXPosition(TrackParticle &p,double &x);
  TVector3 PropogateToYPosition(TrackParticle &p,double &y);
  TVector3 PropogateToZPosition(TrackParticle &p,double &z);
  LorentzVectorParticle LorentzParticleAtPosition(TrackParticle &p,TVector3 &v);
};
#endif


