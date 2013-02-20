#ifndef MultiProngTauSolver_h
#define MultiProngTauSolver_h

#include "TVector3.h" 
#include "TLorentzVector.h" 

class  MultiProngTauSolver {
 public:
  // constructor and Destructor
  MultiProngTauSolver(double mtau_):mtau(mtau_){};
  virtual ~MultiProngTauSolver(){};

  void quadratic(double &x_plus,double &x_minus,double a, double b, double c);
  void AnalyticESolver(TLorentzVector &nu_plus,TLorentzVector &nu_minus,TLorentzVector A1);
  void NumericalESolver(TLorentzVector &nu_plus,TLorentzVector &nu_minus,TLorentzVector A1);
  void SolvebyRotation(TVector3 TauDir,TLorentzVector A1, TLorentzVector &Tau_plus,TLorentzVector &Tau_minus,TLorentzVector &nu_plus,TLorentzVector &nu_minus,bool rotateback=true);
  bool SetTauDirectionatThetaGJMax(TLorentzVector a1, double &theta,double &phi);  

 private: 
  double mtau;
};

#endif
