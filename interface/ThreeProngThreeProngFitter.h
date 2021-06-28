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

  void Configure();

  // enum Pars{taua1_px=0,taua1_py,taua1_pz,taumu_px,taumu_py,taumu_pz,npar};
  enum Pars{tau1_px=0,tau1_py,tau1_pz,tau2_px,tau2_py,tau2_pz,npar};
  enum ExpandedPars{z_vx=npar,z_vy,z_vz,nexpandedpar};
  enum Lambdas{lambda_1=nexpandedpar,lambda_2,nexpandedparpluslambdas};

  enum ParsTrunc{tau_px=0,tau_py,tau_pz,npartr};

  double NConstraints(){return 1;}
  double NSoftConstraints(){if(!useFullRecoil_) return 0; return 2;}
  double NDF(){return NConstraints() + NSoftConstraints();}
  int    NDaughters(){return 2;}
  TString ParName(int par);
  void DebugFit();
  std::vector<LorentzVectorParticle> GetReFitDaughters();
  std::vector<LorentzVectorParticle> GetInitialDaughters() const{return particles0_;};
  LorentzVectorParticle GetMother();
  LorentzVectorParticle GetInitMother(){return Init_Resonance_;};
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
  bool ApplyLagrangianConstraints() override;
  TVectorD ChiSquareUsingInitalPoint(TMatrixT<double> a, TMatrixT<double> b, TMatrixT<double> lambda, TMatrixTSym<double> V_f_inv) override;
  TVectorD HardValue(TVectorD &va,TVectorD &vb);
  TVectorD SoftValue(TVectorD &va,TVectorD &vb);

 private:
  static TMatrixT<double> ComputeInitalExpPar(TMatrixT<double> &inpar);
  static TMatrixT<double> ComputeExpParToPar(TMatrixT<double> &inpar);
  static TMatrixT<double> ComputeExpParToPara(TMatrixT<double> &inpar);
  static TMatrixT<double> ComputeExpParToParb(TMatrixT<double> &inpar);
  static TMatrixT<double> ComputeTau1LorentzVectorPar(TMatrixT<double> &inpar);
  static TMatrixT<double> ComputeTau2LorentzVectorPar(TMatrixT<double> &inpar);
  static TMatrixT<double> ComputeMotherLorentzVectorPar(TMatrixT<double> &inpar);

  void UpdateExpandedPar();
  void CovertParToObjects(TVectorD &va,TVectorD &vb,TLorentzVector &Tau1,TLorentzVector &Tau2);

  TMatrixT<double> exppar_;
  TMatrixTSym<double> expcov_;

  TMatrixT<double> z_;

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
