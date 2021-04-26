#ifndef LagrangeMultipliersFitter_H
#define LagrangeMultipliersFitter_H

#include "TauPolSoftware/SimpleFits/interface/LorentzVectorParticle.h"
#include "TauPolSoftware/SimpleFits/interface/PDGInfo.h"
#include "TMatrixT.h"
#include "TVectorT.h"
#include "TMatrixTSym.h"
#include "TMatrixDEigen.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/FCNBase.h"
#include <vector>

class LagrangeMultipliersFitter{
 public:
  enum Position{pos_x=0,pos_y,pos_z,nposdim};
  enum Parameters{par_vx=0,par_vy,par_vz,par_px,par_py,par_pz,par_m,npardim};
  enum ConvergeProc{ConstraintMin=0,Chi2Min,Chi2AndConstaintMin};
  enum FittingProc{Standard=0,Minuit};

  LagrangeMultipliersFitter();
  virtual ~LagrangeMultipliersFitter(){};

  virtual void   SetWeight(double weight){weight_=weight;}
  virtual void   SetMaxDelta(double MaxDelta){MaxDelta_=MaxDelta;}
  virtual void   SetNIterMax(int Nitermax){nitermax_=Nitermax;}
  virtual void   SetEpsilon(double epsilon){epsilon_=epsilon;}

  virtual bool Fit();
  virtual bool isConverged();
  virtual bool isConfigured(){return isconfigured;}
  virtual double ChiSquare(){return chi2;}
  virtual TVectorD ChiSquareVector(){return chi2_vec;}
  virtual double CSum(){return delta;};
  virtual double NIter(){return niter;};
  virtual double NConstraints()=0;
  virtual double NSoftConstraints()=0;
  virtual double NDF()=0;
  virtual int    NDaughters()=0;
  virtual TString ParName(int par)=0;
  virtual int NPara(){return para_0.GetNrows();};
  virtual int NParb(){return parb_0.GetNrows();};

  virtual std::vector<LorentzVectorParticle> GetReFitDaughters()=0;
  virtual LorentzVectorParticle GetMother()=0;

  static TVectorT<double> convertToVector(TMatrixT<double> M);
  static TMatrixT<double> convertToMatrix(TVectorT<double> V);


   TMatrixT<double> MakeFullMatrix(TMatrixT<double> M11,TMatrixT<double> M12,TMatrixT<double> M21,TMatrixT<double> M22,TMatrixT<double> A,TMatrixT<double> B);
   TMatrixT<double> MakeFullVector(TMatrixT<double> V1,TMatrixT<double> V2,TMatrixT<double> V3);
   TMatrixT<double> solutionlambda(TMatrixT<double> M);
   TMatrixT<double> solutiona(TMatrixT<double> M);
   TMatrixT<double> solutionb(TMatrixT<double> M);
   void Print(TMatrixT<double> M);

  double UpdateChisquare(TVectorD a, TVectorD b);

 protected:
  virtual TVectorD HardValue(TVectorD &va,TVectorD &vb)=0;
  virtual TVectorD SoftValue(TVectorD &va,TVectorD &vb)=0;
  bool  ApplyLagrangianConstraints();

  TVectorD par_0; // parameter values for linearization point
  TVectorD par; // current parameter values
  TMatrixTSym<double> cov_0; //covariance matrix for linearization point (corresponding to par_0) 
  TMatrixTSym<double> cov; // current covariance matrix (corresponding to par) 

  TMatrixT<double> y_;

  //  a and b denotes taua1 and taumu parameters correspondingly
  TVectorD para_0; // parameter values for linearization point
  TVectorD para, paraprev; // current parameter values
  TMatrixTSym<double> cova_0; //covariance matrix for linearization point (corresponding to par_0) 
  TMatrixTSym<double> cova; // current covariance matrix (corresponding to par) 

  TVectorD parb_0; // parameter values for linearization point
  TVectorD parb, parbprev; // current parameter values
  TMatrixTSym<double> covb_0; //covariance matrix for linearization point (corresponding to par_0)
  TMatrixTSym<double> covb; // current covariance matrix (corresponding to par)

  TMatrixT<double> lambda_; // current lagrangian multipliers


  bool isconfigured;
  bool isFit;
  bool useFullRecoil_;
  int fittingMode_;

  // Configuration parameters
  double epsilon_,weight_,MaxDelta_,nitermax_,MaxParDelta_, MaxHCDelta_, MaxSCDelta_, MaxChi2Delta_, nCutStepmax_;

  // Fit variables
  double chi2,chi2prev,delta,niter,pardelta, pardeltaprev;
  TVectorD harddelta_vec, harddelta_vecprev, softdelta_vec, softdelta_vecprev, chi2_vec, chi2_vecprev;
  TMatrixD chi2s;

 private:
  TMatrixT<double> Derivative();

  TMatrixT<double> DerivativeHCa();
  TMatrixT<double> DerivativeHCb();
  TMatrixT<double> DerivativeSCa();
  TMatrixT<double> DerivativeSCb();
  TMatrixTSym<double> ComputeV_f(TMatrixTSym<double>  ca,TMatrixTSym<double>  cb,TVectorD pa,TVectorD pb);
  TMatrixTSym<double> ScaleMatrix(TMatrixTSym<double>  M, double scale);
  double ChiSquare(TMatrixT<double> delta_alpha,TMatrixT<double> lambda,TMatrixT<double> D,TMatrixT<double> d);
  TVectorD ChiSquareUsingInitalPoint(TMatrixT<double> y, TMatrixT<double> a,TMatrixT<double> b,TMatrixT<double> lambda,TMatrixTSym<double> V_f_inv);
  double ConstraintDelta(TVectorT<double> a,TVectorT<double>  b);
  TMatrixT<double> ComputeVariance();
  TMatrixT<double> ComputeVariancea();
  TMatrixT<double> ComputeVarianceb();

  // covariances and derivatives info
  TMatrixT<double> Fa;
  TMatrixT<double> Fb;
  TMatrixT<double> A;
  TMatrixT<double> B;

  TMatrixTSym<double> V_alpha0_inv;
  TMatrixTSym<double> V_f_inv_;
  TMatrixT<double> D;
  TMatrixTSym<double> V_D;
  double ScaleFactor;
  TMatrixT<double> V_corr_prev;


};
#endif
