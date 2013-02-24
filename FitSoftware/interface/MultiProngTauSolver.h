#ifndef MultiProngTauSolver_h
#define MultiProngTauSolver_h

#include "TVector3.h" 
#include "TLorentzVector.h" 
#include "SimpleFits/FitSoftware/interface/LorentzVectorParticle.h"
#include "SimpleFits/FitSoftware/interface/PDGInfo.h"

class  MultiProngTauSolver {
 public:
  enum Ambiguity{zero,minus,plus};

  // constructor and Destructor
  MultiProngTauSolver(double mtau_):mtau(mtau_){};
  MultiProngTauSolver():mtau(PDGInfo::tau_mass()){};
  virtual ~MultiProngTauSolver(){};
      
  void quadratic(double &x_plus,double &x_minus,double a, double b, double c);
  void AnalyticESolver(TLorentzVector &nu_plus,TLorentzVector &nu_minus,TLorentzVector A1);
  void NumericalESolver(TLorentzVector &nu_plus,TLorentzVector &nu_minus,TLorentzVector A1);
  void SolvebyRotation(TVector3 TauDir,TLorentzVector A1, TLorentzVector &Tau_plus,TLorentzVector &Tau_minus,TLorentzVector &nu_plus,TLorentzVector &nu_minus,bool rotateback=true);
  bool SetTauDirectionatThetaGJMax(TLorentzVector a1, double &theta,double &phi,double scale=1.0);  
  static double ThetaGJMax(TLorentzVector a1, double Mtau);
  LorentzVectorParticle EstimateNu(LorentzVectorParticle &a1,TVector3 pv,int ambiguity,TLorentzVector &tau);
 private: 
  double mtau;
};

#endif
