#ifndef DiTauConstrainedFitter_H
#define DiTauConstrainedFitter_H

#include "TauPolSoftware/SimpleFits/interface/LagrangeMultipliersFitter.h"
#include "TVector3.h"
#include "TauPolSoftware/SimpleFits/interface/MultiProngTauSolver.h"
#include "TauPolSoftware/SimpleFits/interface/LorentzVectorParticle.h"
#include "TauPolSoftware/SimpleFits/interface/ErrorMatrixPropagator.h"
#include "TauPolSoftware/SimpleFits/interface/TrackParticle.h"
#include "TauPolSoftware/SimpleFits/interface/PTObject.h"

#include <vector>

class ThreeProngOneProngFitter : public LagrangeMultipliersFitter{
 public:
  ThreeProngOneProngFitter(LorentzVectorParticle TauThreeProng, LorentzVectorParticle ThreeProng, LorentzVectorParticle OneProng, TrackParticle OneProngTrack, double phiz, TVector3 PVertex, TMatrixTSym<double> VertexCov);
  ThreeProngOneProngFitter(LorentzVectorParticle TauThreeProng, LorentzVectorParticle ThreeProng, LorentzVectorParticle OneProng, TrackParticle OneProngTrack, double phiz, TVector3 PVertex, TMatrixTSym<double> VertexCov, double MassConstraint);
  ThreeProngOneProngFitter(LorentzVectorParticle TauThreeProng, LorentzVectorParticle ThreeProng, LorentzVectorParticle OneProng, TrackParticle OneProngTrack, PTObject ResPtEstimate, TVector3 PVertex, TMatrixTSym<double> VertexCov, double MassConstraint);
  virtual ~ThreeProngOneProngFitter(){};

  void Configure(LorentzVectorParticle TauThreeProng, LorentzVectorParticle ThreeProng, LorentzVectorParticle OneProng, TrackParticle OneProngTrack, TVector3 PVertex, TMatrixTSym<double> VertexCov);

  enum Pars{tau3prong_px=0,tau3prong_py,tau3prong_pz,tau1prong_px,tau1prong_py,tau1prong_pz,npar};
  enum ExpandedPars{z_vx=npar,z_vy,z_vz,nexpandedpar};
  enum Lambdas{lambda_1=nexpandedpar,lambda_2,nexpandedparpluslambdas};

  enum ParsTrunc{tau_px=0,tau_py,tau_pz,npartr};

  // Number of Hard Constraints
  virtual double NConstraints(){ return 2;}

  // Number of Soft Constraints
  virtual double NSoftConstraints(){if(!useFullRecoil_) return 3; return 2;}

  // Number of degrees of freedom. NDF = NConstraints + NSoftConstraints - NPar(unmeasured Particles)
  virtual double NDF(){if(!useFullRecoil_) return 2; return 1;}

  // Number of daughter particles from the resonance/boson decay
  virtual int NDaughters(){return 2;}

  // Names of fitted parameters. Used to name parameters inside Minuit2
  virtual TString ParName(int par);

  void DebugFit();
  std::vector<LorentzVectorParticle> GetReFitDaughters();
  std::vector<LorentzVectorParticle> GetInitialDaughters(){return particles0_;};
  LorentzVectorParticle GetMother();
  LorentzVectorParticle GetInitMother(){return Init_Resonance_;};
  LorentzVectorParticle GetTauOneProngEstimate();

  double GetMassConstraint() const{return MassConstraint_;};
  void SetMassConstraint(double MassConstraint) const{MassConstraint_ = MassConstraint;};
  static double MassConstraint_;

  bool GetUseCollinearityTauOneProng() const{return useCollinearityTauOneProng_;};
  void SetUseCollinearityTauOneProng(const bool UseCollinearityTauOneProng) const{useCollinearityTauOneProng_ = UseCollinearityTauOneProng;};
  static bool useCollinearityTauOneProng_;

  TMatrixD GetExppar() const{return exppar;}
  TMatrixDSym GetExpcov() const{return expcov;}

  bool isConverged() override;
  bool Fit() override;

 protected:
  virtual TVectorD HardValue(TVectorD &va,TVectorD &vb);
  virtual TVectorD SoftValue(TVectorD &va,TVectorD &vb);

 private:

  static TMatrixT<double> ComputeInitalExpPar(TMatrixT<double> &inpar);
  static TMatrixT<double> ComputeExpParToPar(TMatrixT<double> &inpar);
  static TMatrixT<double> ComputeExpParToPara(TMatrixT<double> &inpar);
  static TMatrixT<double> ComputeExpParToParb(TMatrixT<double> &inpar);
  static TMatrixT<double> ComputeTauOneProngLorentzVectorPar(TMatrixT<double> &inpar);
  static TMatrixT<double> ComputeTauThreeProngLorentzVectorPar(TMatrixT<double> &inpar);
  static TMatrixT<double> ComputeMotherLorentzVectorPar(TMatrixT<double> &inpar);

  void UpdateExpandedPar();
  void CovertParToObjects(TVectorD &va,TVectorD &vb,TLorentzVector &TauThreeProng,TLorentzVector &Taumu);

  LorentzVectorParticle  TauOneProngStartingPoint(TrackParticle OneProngTrack,LorentzVectorParticle TauThreeProng, TVector3 PV,TMatrixTSym<double>  PVCov, TVector3 SV, TMatrixTSym<double>  SVCov);
  static TMatrixT<double> EstimateTauDirectionAdvanced(TMatrixT<double> &inpar);
  static TMatrixT<double> EstimateTauKinematic(TMatrixT<double> &inpar);
  TMatrixT<double> ConfigureParameters(TrackParticle OneProngTrack, std::pair<double, double> phiAngle);
  TMatrixT<double> ConfigureParameterErrors(TrackParticle OneProngTrack, std::pair<double, double> phiAngle);
  TMatrixT<double> ConfigureInitialAdvancedParameters(TrackParticle OneProngTrack, TVector3 PV, LorentzVectorParticle TauThreeProng);
  TMatrixTSym<double> ConfigureInitialAdvancedParameterErrors(TrackParticle OneProngTrack, TMatrixTSym<double> PVCov,TMatrixTSym<double> TauThreeProngCov);
  TMatrixT<double> ConfigureKinematicParameters(TMatrixT<double>  TauOneProngDir, LorentzVectorParticle TauThreeProng);
  TMatrixTSym<double> ConfigureKinematicParameterErrors(TMatrixTSym<double>  TauOneProngDirError, LorentzVectorParticle TauThreeProng);
  TMatrixT<double> ComputeAngleCovarianceAnalytically(TrackParticle OneProngTrack, std::pair<double, double> phiAngle,  TVector3 PV, TVector3 SV, LorentzVectorParticle  TauThreeProng);
  std::pair<double, double> EstimatePhiAngle( TVector3 dir, TVector3 dirE);

  LorentzVectorParticle TauOneProngStartingPointwithFullRecoil(TrackParticle OneProngTrack,LorentzVectorParticle TauThreeProng, PTObject METminusNeutrino, TVector3 PV, TMatrixTSym<double>  PVCov, TVector3 SV, TMatrixTSym<double>  SVCov);
  TMatrixT<double> ConfigureKinematicParametersFullRecoil(TrackParticle OneProngTrack, TVector3 PV, LorentzVectorParticle TauThreeProng, PTObject METminusNeutrino);
  TMatrixTSym<double> ConfigureKinematicParameterErrorsFullRecoil(TrackParticle OneProngTrack, TMatrixTSym<double> PVCov, LorentzVectorParticle TauThreeProng, PTObject METminusNeutrino);
  static TMatrixT<double> EstimateTauKinematicFullRecoil(TMatrixT<double> &inpar);
  TMatrixDSym EstimateTauKinematicErrorFullRecoil(LorentzVectorParticle TauThreeProng, TLorentzVector TauOneProng, PTObject ResPtEstimate);

  double CosThetaTauOneProng(TLorentzVector TauOneProng);

  // Expanded Parameters
  // exppar_ : 3 Momenta of tau->3prong (a1/3prong+neutral pions) and tau->1prong (rho/pi/mu), i.e. (px_1,py_1,pz_1,px_2,py_2,pz_2) and corresponding covariance expcov_
  TMatrixT<double> exppar_;
  TMatrixTSym<double> expcov_;

  std::vector<LorentzVectorParticle> particles_, particles0_;

  LorentzVectorParticle Init_Resonance_;
  TrackParticle OneProngTrack_;
  LorentzVectorParticle ThreeProng_;
  TVector3 PV_;
  double ThetaForConstrTemporaryImplementation_;
  double phiz_;
  double RecoilX_, RecoilY_;
  PTObject ResPtEstimate_;
  bool debug;
  bool AnalyticalCovariance_;
  int ConstraintMode_;

};
#endif
