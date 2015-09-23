#ifndef DiTauConstrainedFitter_H
#define DiTauConstrainedFitter_H

#include "SimpleFits/FitSoftware/interface/LagrangeMultipliersFitter.h"
#include "TVector3.h"
#include "SimpleFits/FitSoftware/interface/MultiProngTauSolver.h"
#include "SimpleFits/FitSoftware/interface/LorentzVectorParticle.h"
#include "SimpleFits/FitSoftware/interface/ErrorMatrixPropagator.h"
#include "SimpleFits/FitSoftware/interface/TrackParticle.h"
#include "SimpleFits/FitSoftware/interface/PTObject.h"

#include <vector>

class DiTauConstrainedFitter : public LagrangeMultipliersFitter{
 public:
  DiTauConstrainedFitter(LorentzVectorParticle TauA1,TrackParticle MuTrack, double phiz, TVector3 PVertex, TMatrixTSym<double> VertexCov);
  DiTauConstrainedFitter(LorentzVectorParticle TauA1,TrackParticle MuTrack, double phiz, TVector3 PVertex, TMatrixTSym<double> VertexCov, double MassConstraint);
  DiTauConstrainedFitter(LorentzVectorParticle TauA1,TrackParticle MuTrack, PTObject ResPtEstimate, TVector3 PVertex, TMatrixTSym<double> VertexCov, double MassConstraint);
  virtual ~DiTauConstrainedFitter(){};

  void Configure(LorentzVectorParticle TauA1,TrackParticle MuTrack, TVector3 PVertex, TMatrixTSym<double> VertexCov);

  enum Pars{taua1_px=0,taua1_py,taua1_pz,taumu_px,taumu_py,taumu_pz,npar};
  enum ExpandedPars{z_vx=6,z_vy,z_vz,nexpandedpar};

  enum ParsTrunc{tau_px=0,tau_py,tau_pz,npartr};

  virtual double NConstraints(){ return 2;}
  virtual double NSoftConstraints(){if(!useFullRecoil_) return 3; return 2;}
  virtual double NDF(){return 1;}
  virtual int    NDaughters(){return 2;}
  void DebugFit();
  std::vector<LorentzVectorParticle> GetReFitDaughters();
  std::vector<LorentzVectorParticle> GetInitialDaughters(){return particles0_;};
  LorentzVectorParticle GetMother();
  LorentzVectorParticle GetInitMother(){return Init_Resonance_;};
  LorentzVectorParticle GetTauMuEstimate();
  double GetMassConstraint() const{return MassConstraint_;};
  void SetMassConstraint(double MassConstraint) const{MassConstraint_ = MassConstraint;};
  static double MassConstraint_;

  TMatrixD GetExppar() const{return exppar;}
  TMatrixDSym GetExpcov() const{return expcov;}

 protected:
  virtual TVectorD HardValue(TVectorD &va,TVectorD &vb);
  virtual TVectorD SoftValue(TVectorD &va,TVectorD &vb);

 private:

  static TMatrixT<double> ComputeInitalExpPar(TMatrixT<double> &inpar);
  static TMatrixT<double> ComputeExpParToPar(TMatrixT<double> &inpar);
  static TMatrixT<double> ComputeExpParToPara(TMatrixT<double> &inpar);
  static TMatrixT<double> ComputeExpParToParb(TMatrixT<double> &inpar);
  static TMatrixT<double> ComputeTauMuLorentzVectorPar(TMatrixT<double> &inpar);
  static TMatrixT<double> ComputeTauA1LorentzVectorPar(TMatrixT<double> &inpar);
  static TMatrixT<double> ComputeMotherLorentzVectorPar(TMatrixT<double> &inpar);

  void UpdateExpandedPar();
  void CovertParToObjects(TVectorD &va,TVectorD &vb,TLorentzVector &Taua1,TLorentzVector &Taumu, double &Zmass);

  LorentzVectorParticle  TauMuStartingPoint(TrackParticle MuTrack,LorentzVectorParticle TauA1, TVector3 PV,TMatrixTSym<double>  PVCov, TVector3 SV, TMatrixTSym<double>  SVCov);
  static TMatrixT<double> EstimateTauDirectionAdvanced(TMatrixT<double> &inpar);
  static TMatrixT<double> EstimateTauKinematic(TMatrixT<double> &inpar);
  TMatrixT<double> ConfigureParameters(TrackParticle MuTrack, std::pair<double, double> phiAngle);
  TMatrixT<double> ConfigureParameterErrors(TrackParticle MuTrack, std::pair<double, double> phiAngle);
  TMatrixT<double> ConfigureInitialAdvancedParameters(TrackParticle MuTrack, TVector3 PV, LorentzVectorParticle TauA1);
  TMatrixTSym<double> ConfigureInitialAdvancedParameterErrors(TrackParticle MuTrack, TMatrixTSym<double> PVCov,TMatrixTSym<double> TauA1Cov);
  TMatrixT<double> ConfigureKinematicParameters(TMatrixT<double>  TauMuDir, LorentzVectorParticle TauA1);
  TMatrixTSym<double> ConfigureKinematicParameterErrors(TMatrixTSym<double>  TauMuDirError, LorentzVectorParticle TauA1);
  TMatrixT<double> ComputeAngleCovarianceAnalytically(TrackParticle MuTrack, std::pair<double, double> phiAngle,  TVector3 PV, TVector3 SV, LorentzVectorParticle  TauA1);
  std::pair<double, double> EstimatePhiAngle( TVector3 dir, TVector3 dirE);

  LorentzVectorParticle TauMuStartingPointwithFullRecoil(TrackParticle MuTrack,LorentzVectorParticle TauA1, PTObject METminusNeutrino, TVector3 PV, TMatrixTSym<double>  PVCov, TVector3 SV, TMatrixTSym<double>  SVCov);
  TMatrixT<double> ConfigureKinematicParametersFullRecoil(TrackParticle MuTrack, TVector3 PV, LorentzVectorParticle TauA1, PTObject METminusNeutrino);
  TMatrixTSym<double> ConfigureKinematicParameterErrorsFullRecoil(TrackParticle MuTrack, TMatrixTSym<double> PVCov, LorentzVectorParticle TauA1, PTObject METminusNeutrino);
  static TMatrixT<double> EstimateTauKinematicFullRecoil(TMatrixT<double> &inpar);
  TMatrixDSym EstimateTauKinematicErrorFullRecoil(LorentzVectorParticle TauA1, TLorentzVector TauMu, PTObject ResPtEstimate);

  double CosThetaTauMu(TLorentzVector TauMu);

  TMatrixT<double> exppar;
  TMatrixTSym<double> expcov;
  TMatrixT<double> exppara;
  TMatrixTSym<double> expcova;
  TMatrixT<double> expparb;
  TMatrixTSym<double> expcovb;

  std::vector<LorentzVectorParticle> particles_, particles0_;

  LorentzVectorParticle Init_Resonance_;
  TrackParticle MuTrack_;
  TVector3 PV_;
  double ThetaForConstrTemporaryImplementation_;
  double phiz_;
  double RecoilX_, RecoilY_;
  PTObject ResPtEstimate_;
  bool debug;
  bool AnalyticalCovariance;
  int ConstraintMode;



};
#endif
