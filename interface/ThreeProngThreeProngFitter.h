#ifndef ThreeProngThreeProngFitter_H
#define ThreeProngThreeProngFitter_H

#include "TauPolSoftware/SimpleFits/interface/LagrangeMultipliersFitter.h"
#include "TVector3.h"
#include "TauPolSoftware/SimpleFits/interface/MultiProngTauSolver.h"
#include "TauPolSoftware/SimpleFits/interface/LorentzVectorParticle.h"
#include "TauPolSoftware/SimpleFits/interface/ErrorMatrixPropagator.h"
#include "TauPolSoftware/SimpleFits/interface/TrackParticle.h"
#include "TauPolSoftware/SimpleFits/interface/PTObject.h"

#include <vector>

class ThreeProngThreeProngFitter : public LagrangeMultipliersFitter{
 public:
  ThreeProngThreeProngFitter(std::vector< LorentzVectorParticle > TauThreeProngs, std::vector< LorentzVectorParticle > ThreeProngs, TVector3 PVertex, TMatrixTSym<double> VertexCov);
  ThreeProngThreeProngFitter(std::vector< LorentzVectorParticle > TauThreeProngs, std::vector< LorentzVectorParticle > ThreeProngs, PTObject ResPtEstimate, TVector3 PVertex, TMatrixTSym<double> VertexCov);
  ThreeProngThreeProngFitter(std::vector< LorentzVectorParticle > TauThreeProngs, std::vector< LorentzVectorParticle > ThreeProngs, PTObject ResPtEstimate, TVector3 PVertex, TMatrixTSym<double> VertexCov, double MassConstraint);
  virtual ~ThreeProngThreeProngFitter(){};

  void Configure(std::vector< LorentzVectorParticle > TauThreeProngs, std::vector< LorentzVectorParticle > ThreeProngs, TMatrixTSym<double> VertexCov);

  // enum Pars{taua1_px=0,taua1_py,taua1_pz,taumu_px,taumu_py,taumu_pz,npar};
  enum Pars{tau1_px=0,tau1_py,tau1_pz,tau2_px,tau2_py,tau2_pz,npar};
  enum ExpandedPars{z_vx=npar,z_vy,z_vz,nexpandedpar};
  enum Lambdas{lambda_1=nexpandedpar,lambda_2,nexpandedparpluslambdas};

  enum ParsTrunc{tau_px=0,tau_py,tau_pz,npartr};

  virtual double NConstraints(){return 1;}
  virtual double NSoftConstraints(){if(!useFullRecoil_) return 0; return 2;}
  virtual double NDF(){return NConstraints() + NSoftConstraints();}
  virtual int    NDaughters(){return 2;}
  virtual TString ParName(int par);
  void DebugFit();
  std::vector<LorentzVectorParticle> GetReFitDaughters();
  std::vector<LorentzVectorParticle> GetInitialDaughters(){return particles0_;};
  LorentzVectorParticle GetMother();
  LorentzVectorParticle GetInitMother(){return Init_Resonance_;};
  // LorentzVectorParticle GetTauMuEstimate();
  double GetMassConstraint() const{return MassConstraint_;};
  void SetMassConstraint(double MassConstraint) const{MassConstraint_ = MassConstraint;};
  static double MassConstraint_;

  bool GetUseCollinearityTauMu() const{return useCollinearityTauMu_;};
  void SetUseCollinearityTauMu(const bool UseCollinearityTauMu) const{useCollinearityTauMu_ = UseCollinearityTauMu;};
  static bool useCollinearityTauMu_;

  TMatrixD GetExppar() const{return exppar_;}
  TMatrixDSym GetExpcov() const{return expcov_;}


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
  // static std::vector< TMatrixT<double> > ComputeTauThreeProngLorentzVectorPar(TMatrixT<double> &inpar);
  static TMatrixT<double> ComputeTau1LorentzVectorPar(TMatrixT<double> &inpar);
  static TMatrixT<double> ComputeTau2LorentzVectorPar(TMatrixT<double> &inpar);
  static TMatrixT<double> ComputeMotherLorentzVectorPar(TMatrixT<double> &inpar);

  void UpdateExpandedPar();
  void CovertParToObjects(TVectorD &va,TVectorD &vb,TLorentzVector &Tau1,TLorentzVector &Tau2);

  // LorentzVectorParticle  TauMuStartingPoint(TrackParticle MuTrack,std::vector< LorentzVectorParticle > TauThreeProngs, TVector3 PV,TMatrixTSym<double>  PVCov, TVector3 SV, TMatrixTSym<double>  SVCov);
  static TMatrixT<double> EstimateTauDirectionAdvanced(TMatrixT<double> &inpar);
  static TMatrixT<double> EstimateTauKinematic(TMatrixT<double> &inpar);
  TMatrixT<double> ConfigureParameters(TrackParticle MuTrack, std::pair<double, double> phiAngle);
  TMatrixT<double> ConfigureParameterErrors(TrackParticle MuTrack, std::pair<double, double> phiAngle);
  TMatrixT<double> ConfigureInitialAdvancedParameters(TrackParticle MuTrack, TVector3 PV, std::vector< LorentzVectorParticle > TauThreeProngs);
  TMatrixTSym<double> ConfigureInitialAdvancedParameterErrors(TrackParticle MuTrack, TMatrixTSym<double> PVCov,TMatrixTSym<double> TauA1Cov);
  TMatrixT<double> ConfigureKinematicParameters(TMatrixT<double>  TauMuDir, std::vector< LorentzVectorParticle > TauThreeProngs);
  TMatrixTSym<double> ConfigureKinematicParameterErrors(TMatrixTSym<double>  TauMuDirError, std::vector< LorentzVectorParticle > TauThreeProngs);
  TMatrixT<double> ComputeAngleCovarianceAnalytically(TrackParticle MuTrack, std::pair<double, double> phiAngle,  TVector3 PV, TVector3 SV, LorentzVectorParticle  TauA1);
  std::pair<double, double> EstimatePhiAngle( TVector3 dir, TVector3 dirE);

  // LorentzVectorParticle TauMuStartingPointwithFullRecoil(TrackParticle MuTrack,std::vector< LorentzVectorParticle > TauThreeProngs, PTObject METminusNeutrino, TVector3 PV, TMatrixTSym<double>  PVCov, TVector3 SV, TMatrixTSym<double>  SVCov);
  // TMatrixT<double> ConfigureKinematicParametersFullRecoil(TrackParticle MuTrack, TVector3 PV, std::vector< LorentzVectorParticle > TauThreeProngs, PTObject METminusNeutrino);
  // TMatrixTSym<double> ConfigureKinematicParameterErrorsFullRecoil(TrackParticle MuTrack, TMatrixTSym<double> PVCov, std::vector< LorentzVectorParticle > TauThreeProngs, PTObject METminusNeutrino);
  // static TMatrixT<double> EstimateTauKinematicFullRecoil(TMatrixT<double> &inpar);
  // TMatrixDSym EstimateTauKinematicErrorFullRecoil(std::vector< LorentzVectorParticle > TauThreeProngs, TLorentzVector TauMu, PTObject ResPtEstimate);

  // double CosThetaTauMu(TLorentzVector TauMu);

  TMatrixT<double> exppar_;
  TMatrixTSym<double> expcov_;
  // TMatrixT<double> exppara_;
  // TMatrixTSym<double> expcova_;
  // TMatrixT<double> expparb_;
  // TMatrixTSym<double> expcovb_;

  std::vector<LorentzVectorParticle> particles_, particles0_;

  LorentzVectorParticle Init_Resonance_;
  std::vector< LorentzVectorParticle > ThreeProngs_;
  TVector3 PV_;
  TMatrixTSym<double> PVCov_;
  double ThetaForConstrTemporaryImplementation_;
  double phiz_;
  double RecoilX_, RecoilY_;
  PTObject ResPtEstimate_;
  bool debug_;
  bool AnalyticalCovariance_;
  int ConstraintMode_;



};
#endif
