#ifndef TrackHelixVertexFitter_h
#define TrackHelixVertexFitter_h

// system include files
#include <TMatrixT.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TVectorD.h>

#include "SimpleFits/FitSoftware/interface/Particle.h"
#include "SimpleFits/FitSoftware/interface/TrackParticle.h"
#include "SimpleFits/FitSoftware/interface/LorentzVectorParticle.h"
#include "SimpleFits/FitSoftware/interface/ErrorMatrixPropagator.h"
#include "SimpleFits/FitSoftware/interface/TrackHelixVertexTools.h"

class TrackHelixVertexFitter : public TrackHelixVertexTools{
 public:
  TrackHelixVertexFitter(std::vector<TrackParticle> &particles_,TVector3 &vguess);
  virtual ~TrackHelixVertexFitter();

   virtual bool Fit()=0;
   virtual double UpdateChisquare(TMatrixT<double> inpar);
   virtual double ChiSquare(){return chi2;}
   virtual double NDF(){return ndf;}

 protected:
   bool isFit,isConfigure;
   TMatrixT<double> par;
   TMatrixTSym<double> parcov;
   virtual TString FreeParName(int Par);
   double chi2, ndf;

   std::vector<TrackParticle> particles;
   TMatrixT<double> val;
   TMatrixTSym<double> cov;
   TMatrixTSym<double> cov_inv;
   int nParticles;

};
#endif


