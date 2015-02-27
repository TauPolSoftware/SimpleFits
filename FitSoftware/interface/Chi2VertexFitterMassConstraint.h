#ifndef SimpleFits_Chi2VertexFitterMassConstraint_h
#define SimpleFits_Chi2VertexFitterMassConstraint_h

#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/FCNBase.h"
#include "SimpleFits/FitSoftware/interface/TrackHelixVertexTools.h"
#include "SimpleFits/FitSoftware/interface/ConstraintFitterInterator.h"

class  Chi2VertexFitterMassConstraint : public ConstraintFitterInterator, public TrackHelixVertexTools {
 public:
  enum constraintList{mass=0};
  
  Chi2VertexFitterMassConstraint(std::vector<TrackParticle> &particles,TVector3 vguess,double mass,double nsigma_=4.0);
  virtual ~Chi2VertexFitterMassConstraint(){};
  
  virtual bool UpdateChisquare();

 protected:
  virtual bool FitwithConstraints();

 private:   
  double nsigma;
};
#endif


