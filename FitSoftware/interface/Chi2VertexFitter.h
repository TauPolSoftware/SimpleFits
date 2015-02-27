#ifndef Chi2VertexFitter_h
#define Chi2VertexFitter_h

#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/FCNBase.h"
#include "SimpleFits/FitSoftware/interface/TrackHelixVertexTools.h"
#include "SimpleFits/FitSoftware/interface/FitterBase.h"

class  Chi2VertexFitter : public FitterBase, public TrackHelixVertexTools {
 public:
  Chi2VertexFitter(std::vector<TrackParticle> &particles,TVector3 vguess,double nsigma_=4.0);
  virtual ~Chi2VertexFitter(){};

  virtual bool Fit();
  virtual bool UpdateChisquare();

 private:   
  double nsigma;
};
#endif


