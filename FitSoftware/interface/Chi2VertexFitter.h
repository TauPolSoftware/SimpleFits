#ifndef Chi2VertexFitter_h
#define Chi2VertexFitter_h

#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/FCNBase.h"
#include "SimpleFits/FitSoftware/interface/TrackHelixVertexFitter.h"

class  Chi2VertexFitter : public TrackHelixVertexFitter {
 public:
  Chi2VertexFitter(std::vector<TrackParticle> &particles,TVector3 vguess,double nsigma_=4.0):TrackHelixVertexFitter(particles,vguess),nsigma(nsigma_){};
  virtual ~Chi2VertexFitter(){};

  virtual bool Fit();

 private:   
  ROOT::Minuit2::FunctionMinimum FindMinimum(ROOT::Minuit2::FCNBase &updator,ROOT::Minuit2::FunctionMinimum min);
  double nsigma;
};
#endif


