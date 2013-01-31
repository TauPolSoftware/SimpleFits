#ifndef TauA1NuConstrainedFitter_H
#define TauA1NuConstrainedFitter_H

#include "SimpleFits/FitSoftware/interface/LagrangeMultipliersFitter.h"
#include "TVector3.h"
#include "SimpleFits/FitSoftware/interface/MultiProngTauSolver.h"
#include "SimpleFits/FitSoftware/interface/LorentzVectorParticle.h"
#include "SimpleFits/FitSoftware/interface/ErrorMatrixPropagator.h"
#include <vector>

class TauA1NuConstrainedFitter : public LagrangeMultipliersFitter, public MultiProngTauSolver, private ErrorMatrixPropagator {
 public:
  TauA1NuConstrainedFitter(unsigned int ambiguity_,std::vector<LorentzVectorParticle> particles,TVector3 PVertex, TMatrixTSym<double> VertexCov,double mtau);
  virtual ~TauA1NuConstrainedFitter(){};

  enum Pars{tau_phi=0,tau_theta,a1_px,a1_py,a1_pz,a1_m,nu_px,nu_py,nu_pz,npar,norigpar=13};
  enum ExpandedPars{a1_vx=9,a1_vy,a1_vz,nexpandedpar};

  virtual double NConstraints(){return 3;}
  virtual double NDF(){return 0;}
  virtual int    NDaughters(){return 2;}

  std::vector<LorentzVectorParticle> GetReFitDaughters();
  LorentzVectorParticle GetMother();

 protected:
  virtual TVectorD Value(TVectorD &v);
    
 private:
  static TMatrixT<double> ComputeInitalPar(TMatrixT<double> &inpar);
  static TMatrixT<double> ComputeNuLorentzVectorPar(TMatrixT<double> &inpar);
  static TMatrixT<double> ComputeA1LorentzVectorPar(TMatrixT<double> &inpar);
  static TMatrixT<double> ComputeMotherLorentzVectorPar(TMatrixT<double> &inpar);
  void UpdateExpandedPar();
  void SetFitPar();

  double mtau_c;
  TMatrixT<double> exppar;
  TMatrixTSym<double> expcov;
  std::vector<LorentzVectorParticle> particles_;
};
#endif
