#ifndef Chi2VertexFitter_h
#define Chi2VertexFitter_h

#include <TMinuit.h>
#include "SimpleFits/FitSoftware/interface/TrackHelixVertexFitter.h"

class  Chi2VertexFitter : public TrackHelixVertexFitter {
 public:
  Chi2VertexFitter(std::vector<TrackParticle> &particles,TVector3 vguess):TrackHelixVertexFitter(particles,vguess){};
  virtual ~Chi2VertexFitter(){};

  virtual bool Fit();
   
};
#endif


