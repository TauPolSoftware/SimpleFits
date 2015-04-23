#ifndef DiTauConstrainedFitter_H
#define DiTauConstrainedFitter_H

#include "SimpleFits/FitSoftware/interface/LagrangeMultipliersFitter.h"
#include "TVector3.h"
#include "SimpleFits/FitSoftware/interface/MultiProngTauSolver.h"
#include "SimpleFits/FitSoftware/interface/LorentzVectorParticle.h"
#include "SimpleFits/FitSoftware/interface/ErrorMatrixPropagator.h"
#include "SimpleFits/FitSoftware/interface/TrackParticle.h"

#include <vector>

class DiTauConstrainedFitter : public LagrangeMultipliersFitter{
 public:
  DiTauConstrainedFitter(LorentzVectorParticle TauA1,TrackParticle MuTrack, double phiz, TVector3 PVertex, TMatrixTSym<double> VertexCov);
  virtual ~DiTauConstrainedFitter(){};
 
  enum Pars{taua1_px=0,taua1_py,taua1_pz,taumu_px,taumu_py,taumu_pz,npar};
  enum ExpandedPars{z_vx=6,z_vy,z_vz,nexpandedpar};

  enum ParsTrunc{tau_px=0,tau_py,tau_pz,npartr};


  virtual bool Fit();
  virtual double NConstraints(){return 2;}
  virtual double NSoftConstraints(){return 3;}
  virtual double NDF(){return 1;}
  virtual int    NDaughters(){return 2;}
  void DebugFit();
  std::vector<LorentzVectorParticle> GetReFitDaughters();
  LorentzVectorParticle GetMother();
  LorentzVectorParticle GetTauMuEstimate();
  static TVector3 TauMuPtBalanceEstimator(TMatrixT<double> Muon, TVector3 PV, TVector3 SV);
  static double Distance(TVector3 Location1, TVector3 Location2, TVector3 DirectionVector1, TVector3 DirectionVector2);
 protected:
  virtual TVectorD HardValue(TVectorD &va,TVectorD &vb);
  virtual TVectorD SoftValue(TVectorD &va,TVectorD &vb);
    
 private:

    LorentzVectorParticle EstimateTauMu(TVector3 PV,  TMatrixTSym<double>  PVCov, TVector3 SV, TMatrixTSym<double>  SVCov, TrackParticle MuTrack,LorentzVectorParticle TauA1, double ZMassR);


  static TMatrixT<double> ComputeInitalExpPar(TMatrixT<double> &inpar);
  static TMatrixT<double> ComputeExpParToPar(TMatrixT<double> &inpar);
  static TMatrixT<double> ComputeExpParToPara(TMatrixT<double> &inpar);
  static TMatrixT<double> ComputeExpParToParb(TMatrixT<double> &inpar);
  static TMatrixT<double> ComputeTauMuLorentzVectorPar(TMatrixT<double> &inpar);
  static TMatrixT<double> ComputeTauA1LorentzVectorPar(TMatrixT<double> &inpar);
  static TMatrixT<double> ComputeMotherLorentzVectorPar(TMatrixT<double> &inpar);

  std::vector<double>  ReturnDecayPoint(TrackParticle MuTrack, TLorentzVector MuonLV, TVector3 PV,TVector3 SV,TVector3 TauDir,TVector3 TauDirError);
  TVector3 IntersectionPoint(TVector3 PV,TVector3 SV,double xc, double yc, double r);
  LorentzVectorParticle IntersectionPointLinearApproximation(TVector3 PV,TVector3 SV, TVector3 MuonPoca, TLorentzVector MuonLV,TrackParticle MuTrack,TVector3 TauDirError,LorentzVectorParticle TauA1);

  LorentzVectorParticle ConvertTrackParticleToLorentzVectorParticle(TrackParticle MuTrack);

  void UpdateExpandedPar();
  void CovertParToObjects(TVectorD &va,TVectorD &vb,TLorentzVector &Taua1,TLorentzVector &Taumu, double &Zmass);
 
 

  LorentzVectorParticle  TauMuStartingPoint(TrackParticle MuTrack,LorentzVectorParticle TauA1, TVector3 PV,TMatrixTSym<double>  PVCov, TVector3 SV, TMatrixTSym<double>  SVCov );
  static TMatrixT<double> EstimateTauDirectionAdvanced(TMatrixT<double> &inpar);
  static TMatrixT<double> EstimateTauKinematic(TMatrixT<double> &inpar);
  TMatrixT<double> ConfigureParameters(TrackParticle MuTrack, std::pair<double, double> phiAngle);
  TMatrixT<double> ConfigureParameterErrors(TrackParticle MuTrack, std::pair<double, double> phiAngle);
  TMatrixT<double> ConfigureInitialAdvancedParameters(TrackParticle MuTrack, TVector3 PV, TVector3 SV);
  TMatrixTSym<double> ConfigureInitialAdvancedParameterErrors(TrackParticle MuTrack, TMatrixTSym<double>  PVCov,TMatrixTSym<double>  SVCov);
  TMatrixT<double> ConfigureKinematicParameters(TMatrixT<double>  TauMuDir, LorentzVectorParticle TauA1);
  TMatrixTSym<double> ConfigureKinematicParameterErrors(TMatrixTSym<double>  TauMuDirError, LorentzVectorParticle TauA1);
  TMatrixT<double> ComputeAngleCovarianceAnalytically(TrackParticle MuTrack, std::pair<double, double> phiAngle,  TVector3 PV, TVector3 SV, LorentzVectorParticle  TauA1);
  std::pair<double, double> EstimatePhiAngle( TVector3 dir, TVector3 dirE);


 

  TMatrixT<double> exppar;
  TMatrixTSym<double> expcov;
  TMatrixT<double> exppara;
  TMatrixTSym<double> expcova;
  TMatrixT<double> expparb;
  TMatrixTSym<double> expcovb;



  std::vector<LorentzVectorParticle> particles_;
  double ThetaForConstrTemporaryIMplementation_;
  double phiz_;
  bool debug;
  bool AnalyticalCovariance;
  int ConstraintMode;



};
#endif
