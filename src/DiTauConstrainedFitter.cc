#include "TauPolSoftware/SimpleFits/interface/DiTauConstrainedFitter.h"
#include "TauPolSoftware/SimpleFits/interface/ChiSquareFunctionUpdator.h"
#include "TauPolSoftware/SimpleFits/interface/PDGInfo.h"
#include "TauPolSoftware/SimpleFits/interface/Logger.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnSimplex.h"
#include "Minuit2/CombinedMinimizer.h"
#include "Minuit2/MnMinimize.h"
#include "Minuit2/MnScan.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnHesse.h"
#include "Minuit2/MnContours.h"
#include "Minuit2/MnPlot.h"
#include "Minuit2/MinosError.h"
#include "Minuit2/ContoursError.h"
#include <iostream>

double DiTauConstrainedFitter::MassConstraint_ = 91.5;
bool DiTauConstrainedFitter::useCollinearityTauMu_ = false;

DiTauConstrainedFitter::DiTauConstrainedFitter(LorentzVectorParticle TauA1, LorentzVectorParticle A1, TrackParticle MuTrack, double phiz, TVector3 PVertex, TMatrixTSym<double> VertexCov){
  DiTauConstrainedFitter(TauA1, A1, MuTrack, phiz, PVertex, VertexCov, 91.5);
}

DiTauConstrainedFitter::DiTauConstrainedFitter(LorentzVectorParticle TauA1, LorentzVectorParticle A1, TrackParticle MuTrack, double phiz, TVector3 PVertex, TMatrixTSym<double> VertexCov, double MassConstraint):
  LagrangeMultipliersFitter()
{
  debug = false;
  AnalyticalCovariance =false;

  phiz_ = phiz;
  RecoilX_ = 0; //not used in this version
  RecoilY_ = 0; //not used in this version
  useFullRecoil_ = false;
  MassConstraint_ = MassConstraint;

  Configure(TauA1, A1, MuTrack, PVertex, VertexCov);
}

DiTauConstrainedFitter::DiTauConstrainedFitter(LorentzVectorParticle TauA1, LorentzVectorParticle A1, TrackParticle MuTrack, PTObject ResPtEstimate, TVector3 PVertex, TMatrixTSym<double> VertexCov, double MassConstraint):
  LagrangeMultipliersFitter()
{
  debug = false;
  AnalyticalCovariance =false;

  ResPtEstimate_ = ResPtEstimate;
  phiz_ = 0; //not used in this version
  RecoilX_ = ResPtEstimate.Par()(0,0);
  RecoilY_ = ResPtEstimate.Par()(1,0);
  useFullRecoil_ = true;
  MassConstraint_ = MassConstraint;

  Configure(TauA1, A1, MuTrack, PVertex, VertexCov);
}

void DiTauConstrainedFitter::Configure(LorentzVectorParticle TauA1, LorentzVectorParticle A1, TrackParticle MuTrack, TVector3 PVertex, TMatrixTSym<double> VertexCov){
  debug = false;
  AnalyticalCovariance =false;
  MuTrack_ = MuTrack;
  A1_ = A1;
  PV_ = PVertex;

  LorentzVectorParticle  TauMuGuess;
  LorentzVectorParticle  DiTau;
  if(!useFullRecoil_) TauMuGuess  = TauMuStartingPoint( MuTrack,TauA1,PVertex, VertexCov, TauA1.Vertex(),TauA1.VertexCov());
  else{
	TauMuGuess  = TauMuStartingPointwithFullRecoil(MuTrack,TauA1, ResPtEstimate_, PVertex, VertexCov, TauA1.Vertex(),TauA1.VertexCov());
  }

  Logger(Logger::Debug) << "TauA1 covariance: " << std::endl;
  if(Logger::Instance()->Level() == Logger::Debug){
	TauA1.getCovMatrix().Print();
  }
  Logger(Logger::Debug) << "TauMuGuess covariance: " << std::endl;
  if(Logger::Instance()->Level() == Logger::Debug){
	TauMuGuess.getCovMatrix().Print();
  }
  Logger(Logger::Debug) << "Passed: TauMuStartingPoint" << std::endl;

  ThetaForConstrTemporaryImplementation_=TauMuGuess.LV().Theta();
  particles_.push_back(TauA1);
  particles_.push_back(TauMuGuess);

  particles0_.push_back(TauA1);
  particles0_.push_back(TauMuGuess);

  Logger(Logger::Debug) << "(TauA1.LV() + TauMuGuess.LV()).M(): " << (TauA1.LV() + TauMuGuess.LV()).M() << std::endl;

  isconfigured=false;

  int size=particles_.size()*3;
  int sizeTrunc = 3;
  TMatrixT<double>    inpar(size,1);
  TMatrixTSym<double> incov(size);

  // Get primary vertex information
  if(VertexCov.GetNrows()!=LorentzVectorParticle::NVertex)return;

  //set input paramterts:  TauA1 - TauMu
  inpar(taua1_px,0)=TauA1.LV().Px();
  inpar(taua1_py,0)=TauA1.LV().Py();
  inpar(taua1_pz,0)=TauA1.LV().Pz();

  inpar(taumu_px,0)=TauMuGuess.LV().Px();
  inpar(taumu_py,0)=TauMuGuess.LV().Py();
  inpar(taumu_pz,0)=TauMuGuess.LV().Pz();

  int TauA1offset=0;
  int TauMuoffset=3;

  for(int i=0; i<3;i++){
    for(int j=0; j<3;j++){
      incov(i+TauA1offset,j+TauA1offset)=TauA1.Covariance(i+3,j+3);
      incov(i+TauMuoffset,j+TauMuoffset)=TauMuGuess.Covariance(i+3,j+3);
    }
  }

// store expanded par for computation of final par (assumes fit has neglegible impact on a1 correlations with vertex uncertainties)

  exppar.ResizeTo(size,1);
  expcov.ResizeTo(size,size);
  exppar=ComputeInitalExpPar(inpar);
  expcov=ErrorMatrixPropagator::PropagateError(&DiTauConstrainedFitter::ComputeInitalExpPar,inpar,incov);

  // store linearization point
  TMatrixT<double> PAR_0(size,1);
  par_0.ResizeTo(size);
  cov_0.ResizeTo(size,size);
  PAR_0=ComputeExpParToPar(exppar);
  for(int i=0; i<npar;i++)par_0(i)=PAR_0(i,0);
  cov_0=ErrorMatrixPropagator::PropagateError(&DiTauConstrainedFitter::ComputeExpParToPar,exppar,expcov);

  // set up inital point for fit (cov handled in Fit() function)
  par.ResizeTo(npar);
  par=par_0;


  TMatrixT<double> PARa_0(sizeTrunc,1);
  para_0.ResizeTo(sizeTrunc);
  cova_0.ResizeTo(sizeTrunc,sizeTrunc);
  PARa_0=ComputeExpParToPara(exppar);
  for(int i=0; i<npartr;i++)para_0(i)=PARa_0(i,0);

  cova_0=ErrorMatrixPropagator::PropagateError(&DiTauConstrainedFitter::ComputeExpParToPara,exppar,expcov);


  para.ResizeTo(npartr);
  para=para_0;
  TMatrixT<double> PARb_0(sizeTrunc,1);
  parb_0.ResizeTo(sizeTrunc);
  covb_0.ResizeTo(sizeTrunc,sizeTrunc);
  PARb_0=ComputeExpParToParb(exppar);
  for(int i=0; i<npartr;i++)parb_0(i)=PARb_0(i,0);
  y_.ResizeTo(npartr,1); y_ = convertToMatrix(para_0);

  covb_0=ErrorMatrixPropagator::PropagateError(&DiTauConstrainedFitter::ComputeExpParToParb,exppar,expcov);
  parb.ResizeTo(npartr);
  parb=parb_0;

  isconfigured=true;
  Init_Resonance_ = GetMother();


  Logger(Logger::Debug) << "exppar covariance: " << std::endl;
  if(Logger::Instance()->Level() == Logger::Debug){
	expcov.Print();
  }
  Logger(Logger::Debug) << "para_0 covariance: " << std::endl;
  if(Logger::Instance()->Level() == Logger::Debug){
	cova_0.Print();
  }
  Logger(Logger::Debug) << "parb covariance: " << std::endl;
  if(Logger::Instance()->Level() == Logger::Debug){
	covb_0.Print();
  }
}

TMatrixT<double> DiTauConstrainedFitter::ComputeInitalExpPar(TMatrixT<double> &inpar){
  TMatrixT<double> outpar(6,1);
  int offset=0;//LorentzVectorParticle::NVertex;// for TauA1
  outpar(taua1_px,0)=inpar(LorentzVectorParticle::vx+offset,0);
  outpar(taua1_py,0)=inpar(LorentzVectorParticle::vy+offset,0);
  outpar(taua1_pz,0)=inpar(LorentzVectorParticle::vz+offset,0);

  offset+=3;//LorentzVectorParticle::NLorentzandVertexPar; // for TauNu

  outpar(taumu_px,0)=inpar(LorentzVectorParticle::vx+offset,0);
  outpar(taumu_py,0)=inpar(LorentzVectorParticle::vy+offset,0);
  outpar(taumu_pz,0)=inpar(LorentzVectorParticle::vz+offset,0);

  return outpar;

}

bool DiTauConstrainedFitter::Fit(){
  if(cov.GetNrows()!=par_0.GetNrows()){
    cov.ResizeTo(par_0.GetNrows(),par_0.GetNrows());
    cov=cov_0;
  }
  if (fittingMode_ == FittingProc::Standard){
    // Standard Fitting method
    if(!isconfigured) return false;
    if(isFit)return isConverged();
    isFit=true;
    // chi2s.ResizeTo((int)nitermax_,4);
    chi2_vec.ResizeTo(3);
    harddelta_vec.ResizeTo(NConstraints());

    Logger(Logger::Debug) << "Start: ApplyLagrangianConstraints "<< std::endl;

    for(niter=0; niter < nitermax_; niter++){
      Logger(Logger::Debug) << "Start: ApplyLagrangianConstraints iteration no.: " << niter << std::endl;
      bool passed=ApplyLagrangianConstraints();

      bool cut_step(false);
      bool directing_step(false);
      // if(harddelta_vecprev.Norm1() < harddelta_vec.Norm1() && harddelta_vec.Norm1() > MaxHCDelta_) cut_step = true;//CutStep method
      if(chi2prev < chi2 && chi2prev > 0 && fabs(chi2prev-chi2) > MaxChi2Delta_) cut_step = true;
      if(chi2 < 0 && chi2prev > 0) directing_step = false;
      if(cut_step || directing_step){
        if(cut_step){
          Logger(Logger::Debug) << "Starting CutStep Method" << std::endl;
        }
        else if(directing_step){
          Logger(Logger::Debug) << "Starting DirectingStep Method" << std::endl;
        }
        TVectorD paraprev_(paraprev); //save copies because para/bprev are changed at the beginning of ApplyLagrangianConstraints
        TVectorD parbprev_(parbprev);
        TVectorD para_delta_(para_0 - paraprev_);
        TVectorD parb_delta_(parb_0 - parbprev_);
        TVectorD para_delta(para_delta_);
        TVectorD parb_delta(parb_delta_);
        TVectorD harddelta_vecprev_(harddelta_vecprev);
        double chi2CutStep(chi2prev);
        TVectorD chi2CutStep_vec(chi2_vecprev);
        // double chi2CutStepprev(chi2prev);
        if(Logger::Instance()->Level() == Logger::Debug){
          Logger(Logger::Debug) << "para_delta: " << std::endl;
          para_delta.Print();
          Logger(Logger::Debug) << "parb_delta: " << std::endl;
          parb_delta.Print();
        }
        for(unsigned i_CutStep = 0; i_CutStep <= nCutStepmax_; i_CutStep++){
          //double corrFactor = 1/pow(2., (double) i_CutStep + 1);
          double corrFactor = (double) i_CutStep/(double) nCutStepmax_;
          corrFactor = directing_step ? corrFactor : 1-corrFactor;
          // double corrFactor = directing_step ? i_CutStep/nCutStepmax_ : 1-pow(((double)i_CutStep/nCutStepmax_),2.);
          Logger(Logger::Debug) << "corrFactor: " << corrFactor << " after " << i_CutStep << " Steps" << std::endl;
          para_delta = para_delta_;
          parb_delta = parb_delta_;
          para_delta *= corrFactor;
          parb_delta *= corrFactor;
          paraprev = paraprev_;
          parbprev = parbprev_;
          para_0 = paraprev_ + para_delta;
          parb_0 = parbprev_ + parb_delta;
          // if(Logger::Instance()->Level() == Logger::Debug){
          // Logger(Logger::Debug) << "para_0: after i = " << i_CutStep << " Cutsteps: " << std::endl;
          // para_0.Print();
          // Logger(Logger::Debug) << "para_delta after i = " << i_CutStep << " Cutsteps: " << std::endl;
          // para_delta.Print();
          // Logger(Logger::Debug) << "parb_0: after i = " << i_CutStep << " Cutsteps: " << std::endl;
          // parb_0.Print();
          // Logger(Logger::Debug) << "parb_delta after i = " << i_CutStep << " Cutsteps: " << std::endl;
          // parb_delta.Print();
          // }
          passed = ApplyLagrangianConstraints();
          Logger(Logger::Debug) << "passed: " << passed << ", chi2_vec(0): " << chi2_vec(0) << ", chi2_vec(1): " << chi2_vec(1) << ", chi2_vec(2): " << chi2_vec(2) << std::endl;
          // if(cut_step && passed && harddelta_vecprev_.Norm1() > harddelta_vec.Norm1()){
          if(cut_step && passed && chi2 < chi2CutStep && chi2_vec(0) > 0 && chi2_vec(1) > 0 && chi2_vec(2) > -1*MaxHCDelta_){
            Logger(Logger::Debug) << "Stopped CutStep Method after "<< i_CutStep << " iterations" << std::endl;
            Logger(Logger::Debug) << "Chi2 before CutStep Method = "<< chi2CutStep << " and Chi2 after = " << chi2 << std::endl;
            break;
          }
          if(directing_step && passed && chi2 < chi2CutStep && chi2_vec(0) > 0 && chi2_vec(1) > 0 && chi2_vec(2) > -1*MaxHCDelta_){
            Logger(Logger::Debug) << "Stopped DirectingStep Method after " << i_CutStep << " iterations" << std::endl;
            Logger(Logger::Debug) << "Chi2 before DirectingStep Method = " << chi2CutStep << " and Chi2 after = " << chi2 << std::endl;
            break;
          }
        }
      }

      // chi2s(niter,0) = chi2_vec(0);
      // chi2s(niter,1) = harddelta_vec.Norm1();
      // chi2s(niter,2) = harddelta_vec(0);
      // chi2s(niter,3) = harddelta_vec(1);

      if(!passed){
        Logger(Logger::Error) << "Did not pass ApplyLagrangianConstraints(). Matrix inversion failed probably." << std::endl;
        return false;
      }

      // if(isConverged() && !cut_step && !directing_step){
      if(isConverged()){
        break;
      }
      if(niter==nitermax_-1){
        Logger(Logger::Debug) << "Reached Maximum number of iterations = " << niter << std::endl;
        if(harddelta_vec.Norm1() > MaxHCDelta_){
          Logger(Logger::Debug) << "Deviations from hard constraints = " << harddelta_vec.Norm1() << " greater than max. allowed = " << MaxHCDelta_ << std::endl;
          Logger(Logger::Debug) << "Deviations from hard constraints = " << harddelta_vec(0) << " and = " << harddelta_vec(1) << std::endl;
          }
        if(fabs(chi2prev-chi2) > MaxChi2Delta_){
          Logger(Logger::Debug) << "chi2prev-chi2 = " << chi2prev-chi2 << " greater than max. allowed = " << MaxChi2Delta_ << std::endl;
        }
        //chi2s.Print();
        return false;
        }
        if(harddelta_vecprev.Norm1()>=1e12) {
          Logger(Logger::Debug) << "Huge deviations from hard constraints --> overshoot during chi2 min. niter = "<< niter << " and delta = "<< harddelta_vecprev.Norm1() <<std::endl;
          return false;
        }
    }

    if(Logger::Instance()->Level() == Logger::Debug){
      Logger(Logger::Debug) << "para_0: " << std::endl;
      para_0.Print();
      Logger(Logger::Debug) << "parb_0: " << std::endl;
      parb_0.Print();
    }

     // ComputeVariancea();
     // ComputeVarianceb();

    return true;
  }
  else if (fittingMode_ == FittingProc::Minuit){
    // Fitting method Minuit2
    // TODO: actually implement properly. Only copypasted from Chi2VertexFitter atm
    if(isFit==true) return true;// do not refit
    if(!isconfigured) return false; // do not fit if configuration failed
    unsigned int nPar = para_0.GetNrows() + parb_0.GetNrows();
    DiTauConstrainedFitterChiSquareFunctionUpdator updator(this);
    ROOT::Minuit2::MnUserParameters MnPar;
    ROOT::Minuit2::MnUserCovariance MnCov(nPar);
    // Logger(Logger::Info) << "Starting Fit with FittingProc::Minuit" << std::endl;
    // Logger(Logger::Info) << "Setting up parameters a:" << std::endl;
    for(int i=0;i<para_0.GetNrows();i++){
      TString name=ParName(i);
      // if not limited (vhigh <= vlow)
      // MnPar.Add(name.Data(),para_0(i,0),sqrt(fabs(cova_0(i,i))),para_0(i,0)-nsigma*sqrt(fabs(cova_0(i,i))),para_0(i,0)+nsigma*sqrt(fabs(cova_0(i,i))));
      MnPar.Add(name.Data(),para_0(i));//,sqrt(fabs(cova_0(i,i))));
      for(int j=0;j<para_0.GetNrows();j++){
        MnCov(i,j) = cova_0(i,j);
      }
      // Logger(Logger::Info) << "\t" << name << std::endl;
    }
    // Logger(Logger::Info) << "Setting up parameters b" << std::endl;
    for(int i=0;i<parb_0.GetNrows();i++){
      int offset = para_0.GetNrows();
      TString name=ParName(i+offset);
      // if not limited (vhigh <= vlow)
      // MnPar.Add(name.Data(),parb_0(i,0),sqrt(fabs(covb_0(i,i))),parb_0(i,0)-nsigma*sqrt(fabs(covb_0(i,i))),parb_0(i,0)+nsigma*sqrt(fabs(covb_0(i,i))));
      MnPar.Add(name.Data(),parb_0(i)); //,sqrt(fabs(covb_0(i,i))));
      for(int j=0;j<parb_0.GetNrows();j++){
        MnCov(i+offset,j+offset) = covb_0(i,j);
      }
      // Logger(Logger::Info) << "\t" << name << std::endl;
    }
    // set limits for tau mu pz
    MnPar.SetLimits(5, -5.0*MuTrack_.P(), 5.0*MuTrack_.P());

    unsigned int max=10;
    // int numberofcalls=200+nPar*100+nPar*nPar*5;
    int numberofcalls=100000;
    double tolerance(1.0);
    double edmMin(0.001*updator.Up()*tolerance);

    // ROOT::Minuit2::MnPrint MnLogger;
    // MnLogger.SetLevel(0);
    // Logger(Logger::Info) << "Begin minimization" << std::endl;
    ROOT::Minuit2::MnMinimize minimize(updator, MnPar, MnCov);
    minimize.Fix(5);
    // ROOT::Minuit2::FunctionMinimum min= minimize(numberofcalls,tolerance);
    ROOT::Minuit2::FunctionMinimum min0 = minimize(numberofcalls,tolerance);
    minimize.Release(5);
    ROOT::Minuit2::FunctionMinimum min1 = minimize(numberofcalls,tolerance);
    minimize.RemoveLimits(5);
    ROOT::Minuit2::FunctionMinimum min = minimize(numberofcalls,tolerance);
    // for(unsigned int i=0;i<=max && min.Edm()>edmMin;i++){
    //   if(i==max) return false;
    //   min = minimize(i*numberofcalls,tolerance);
    // }
    // give return flag based on status
    if(min.IsAboveMaxEdm()){Logger(Logger::Error) << "Found Solution that is above EDM " << std::endl; return false;}
    if(!min.IsValid()){
      Logger(Logger::Error) << "Failed min.IsValid()" << std::endl;
      if(!min.HasValidParameters()){Logger(Logger::Error) << "Failed min.HasValidParameters()" << std::endl; }
      if(!min.HasValidCovariance()){Logger(Logger::Error) << "Failed min.HasValidCovariance()" << std::endl; }
      if(!min.HesseFailed()){Logger(Logger::Error) << "Failed min.HesseFailed()" << std::endl; }
      if(!min.HasReachedCallLimit()){Logger(Logger::Error) << "Failed min.HasReachedCallLimit()" << std::endl; }
      return false;
    }
    chi2=min.Fval();
    // Logger(Logger::Error) << "minimum with fixed tau mu pz: " << min0 << std::endl;
    // Logger(Logger::Error) << "minimum with limited tau mu pz: " << " 5.0*MuTrack_.P(): " << 5.0*MuTrack_.P() << min1 << std::endl;
    // Logger(Logger::Error) << "minimum with all free parameters: " << min << std::endl;
    // // Get output parameters
    // for(int i=0;i<par.GetNrows();i++){ par(i,0)=min.UserParameters().Value(i);}
    // // Get output covariance
    // for(int i=0;i<par.GetNrows();i++){
    //   for(int j=0;j<par.GetNrows();j++){parcov(i,j)=min.UserCovariance()(i,j);}
    // }
    // Logger(Logger::Info) << "Saving parameters a" << std::endl;
    for(int i=0;i<para_0.GetNrows();i++){
      para(i) = min.UserParameters().Value(i);
      // for(int j=0;j<para_0.GetNrows();j++){
      //   cova(i,j) = min.UserCovariance()(i,j);
      // }
    }
    // Logger(Logger::Info) << "Saving parameters b" << std::endl;
    for(int i=0;i<parb_0.GetNrows();i++){
      int offset = para_0.GetNrows();
      parb(i) = min.UserParameters().Value(i+offset);
      // for(int j=0;j<parb_0.GetNrows();j++){
      //   covb(i,j) = min.UserCovariance()(i+offset,j+offset);
      // }
    }
    isFit=true;
    return isFit;
  }
  return false;
}

TMatrixT<double> DiTauConstrainedFitter::ComputeExpParToPar(TMatrixT<double> &inpar){
  TMatrixT<double> outpar(npar,1);
  for(int i=0;i<npar;i++){outpar(i,0)=inpar(i,0);}
  return outpar;
}
TMatrixT<double> DiTauConstrainedFitter::ComputeExpParToPara(TMatrixT<double> &inpar){
  TMatrixT<double> outpar(npartr,1);
  for(int i=0;i<npartr;i++){outpar(i,0)=inpar(i,0);}
  return outpar;
}

TMatrixT<double> DiTauConstrainedFitter::ComputeExpParToParb(TMatrixT<double> &inpar){
  TMatrixT<double> outpar(npartr,1);
  int offset = 3;
  for(int i=0;i<npartr;i++){outpar(i,0)=inpar(i+offset,0);}
  return outpar;
}

TMatrixT<double> DiTauConstrainedFitter::ComputeTauMuLorentzVectorPar(TMatrixT<double> &inpar){
  //  start with index 3 to fill only momenta part
  TMatrixT<double> outpar(7,1);
  outpar(LorentzVectorParticle::vx,0)=0;
  outpar(LorentzVectorParticle::vy,0)=0;
  outpar(LorentzVectorParticle::vz,0)=0;

  outpar(LorentzVectorParticle::px,0)=inpar(3,0);
  outpar(LorentzVectorParticle::py,0)=inpar(4,0);
  outpar(LorentzVectorParticle::pz,0)=inpar(5,0);
  outpar(LorentzVectorParticle::m,0)=PDGInfo::tau_mass();

  return outpar;
}

TMatrixT<double> DiTauConstrainedFitter::ComputeTauA1LorentzVectorPar(TMatrixT<double> &inpar){
  //  start with index 3 to fill only momenta part
  TMatrixT<double> outpar(7,1);
  outpar(LorentzVectorParticle::vx,0)=0;
  outpar(LorentzVectorParticle::vy,0)=0;
  outpar(LorentzVectorParticle::vz,0)=0;
  outpar(LorentzVectorParticle::px,0)=inpar(0,0);
  outpar(LorentzVectorParticle::py,0)=inpar(1,0);
  outpar(LorentzVectorParticle::pz,0)=inpar(2,0);
  outpar(LorentzVectorParticle::m,0) =PDGInfo::tau_mass();

  return outpar;
}

TMatrixT<double> DiTauConstrainedFitter::ComputeMotherLorentzVectorPar(TMatrixT<double> &inpar){

  TMatrixT<double> outpar(7,1);
  TMatrixT<double> Taumupar=ComputeTauMuLorentzVectorPar(inpar);
  TMatrixT<double> Taua1par=ComputeTauA1LorentzVectorPar(inpar);
  outpar(LorentzVectorParticle::px,0)=Taumupar(LorentzVectorParticle::px,0)+Taua1par(LorentzVectorParticle::px,0);
  outpar(LorentzVectorParticle::py,0)=Taumupar(LorentzVectorParticle::py,0)+Taua1par(LorentzVectorParticle::py,0);
  outpar(LorentzVectorParticle::pz,0)=Taumupar(LorentzVectorParticle::pz,0)+Taua1par(LorentzVectorParticle::pz,0);

  double Etaumu2=pow(Taumupar(LorentzVectorParticle::px,0),2.0)+pow(Taumupar(LorentzVectorParticle::py,0),2.0)+pow(Taumupar(LorentzVectorParticle::pz,0),2.0)+pow(Taumupar(LorentzVectorParticle::m,0),2.0);
  double Etaua12=pow(Taua1par(LorentzVectorParticle::px,0),2.0)+pow(Taua1par(LorentzVectorParticle::py,0),2.0)+pow(Taua1par(LorentzVectorParticle::pz,0),2.0)+pow(Taua1par(LorentzVectorParticle::m,0),2.0);
  double P2=pow(outpar(LorentzVectorParticle::px,0),2.0)+pow(outpar(LorentzVectorParticle::py,0),2.0)+pow(outpar(LorentzVectorParticle::pz,0),2.0);
  outpar(LorentzVectorParticle::m,0)=sqrt(fabs(pow(sqrt(Etaumu2)+sqrt(Etaua12),2.0)-P2));

  return outpar;
}

void DiTauConstrainedFitter::UpdateExpandedPar(){
  for(int i=0; i<npartr;i++){
	  exppar(i,0)=para(i);
	  for(int j=0; j<npartr;j++){expcov(i,j)=cova_0(i,j);}
  }
  int offset = npartr;
  for(int i=0; i<npartr;i++){
	  exppar(i+offset,0)=parb(i);
	  for(int j=0; j<npartr;j++){expcov(i+offset,j+offset)=covb_0(i,j);}
  }

}

std::vector<LorentzVectorParticle> DiTauConstrainedFitter::GetReFitDaughters(){
  std::vector<LorentzVectorParticle> refitParticles;
  UpdateExpandedPar();

  double c(0),b(0);
  for(unsigned int i=0;i<particles_.size();i++){c+=particles_.at(i).Charge();b=particles_.at(i).BField();}
  TMatrixT<double> a1=ComputeTauA1LorentzVectorPar(exppar);
  TMatrixTSym<double> a1cov=ErrorMatrixPropagator::PropagateError(&DiTauConstrainedFitter::ComputeTauA1LorentzVectorPar,exppar,expcov);

   for(int i=0; i<LorentzVectorParticle::NVertex; i++){
     for(int j=0; j<LorentzVectorParticle::NVertex; j++){
       a1cov(i,j)=particles_.at(0).VertexCov()(i,j);
     }
   }

    refitParticles.push_back(LorentzVectorParticle(a1,a1cov,PDGInfo::tau_plus,c,b));

    TMatrixT<double> mu=ComputeTauMuLorentzVectorPar(exppar);
    mu(0,0)= particles_.at(1).Parameter(LorentzVectorParticle::vx);
    mu(1,0)= particles_.at(1).Parameter(LorentzVectorParticle::vy);
    mu(2,0)= particles_.at(1).Parameter(LorentzVectorParticle::vz);

    TMatrixTSym<double> mucov=ErrorMatrixPropagator::PropagateError(&DiTauConstrainedFitter::ComputeTauMuLorentzVectorPar,exppar,expcov);
    for(int i=0; i<LorentzVectorParticle::NVertex; i++){
      for(int j=0; j<LorentzVectorParticle::NVertex; j++){
	mucov(i,j)=particles_.at(1).VertexCov()(i,j);
      }
    }

    refitParticles.push_back(LorentzVectorParticle(mu,mucov,PDGInfo::tau_minus,0.0,b));

  return refitParticles;
}

LorentzVectorParticle
DiTauConstrainedFitter::GetMother(){
  UpdateExpandedPar();
  double c(0),b(0);
  for(unsigned int i=0;i<particles_.size();i++){c+=particles_.at(i).Charge();b=particles_.at(i).BField();}
  TMatrixT<double> m=ComputeMotherLorentzVectorPar(exppar);
  TMatrixTSym<double> mcov=ErrorMatrixPropagator::PropagateError(&DiTauConstrainedFitter::ComputeMotherLorentzVectorPar,exppar,expcov);

  return LorentzVectorParticle(m,mcov,PDGInfo::Z0,c,b);
}

TVectorD
DiTauConstrainedFitter::HardValue(TVectorD &va,TVectorD &vb,bool debug){
  TLorentzVector Taua1,Taumu;
  double ZMass;
  CovertParToObjects(va,vb,Taua1,Taumu,ZMass);

  TLorentzVector z=Taua1+Taumu;
  TVectorD d(NConstraints());

	//d(0) = pow(Taua1.E() + Taumu.E(), 2.) - (Taua1.Vect() + Taumu.Vect()).Mag2()- pow(MassConstraint_,2.);
	//d(1) = Taumu.Pz() - CosThetaTauMu(Taumu)*Taumu.P();
	d(0) = z.M() - MassConstraint_;
	d(1) = Taumu.Pz() - CosThetaTauMu(Taumu)*Taumu.P();
  return d;
}

TVectorD
DiTauConstrainedFitter::SoftValue(TVectorD &va,TVectorD &vb,bool debug){
  TLorentzVector Taua1,Taumu;
  double ZMass;
  CovertParToObjects(va,vb,Taua1,Taumu,ZMass);
  TVectorD d(NSoftConstraints());

  if(!useFullRecoil_){
	d(0) = Taua1.Px() + Taumu.Px();
	d(1) = Taua1.Py() + Taumu.Py();
	d(2) = ( (Taua1.Py() + Taumu.Py())/(Taua1.Px() + Taumu.Px())) -  tan(phiz_);
	Logger(Logger::Debug) << "SCVec: " << d(0) << ", " << d(1) << ", " << d(2) << std::endl;
	Logger(Logger::Debug) << "RecoilX_: " << RecoilX_ << ", RecoilY_: " << RecoilY_ << std::endl;
  }
  else{
	d(0) = Taua1.Px() + Taumu.Px() - RecoilX_;
	d(1) = Taua1.Py() + Taumu.Py() - RecoilY_;
	Logger(Logger::Debug) << "SCVec: " << d(0) << ", " << d(1) << std::endl;
	Logger(Logger::Debug) << "RecoilX_: " << RecoilX_ << ", RecoilY_: " << RecoilY_ << std::endl;
  }
  return d;
}

void DiTauConstrainedFitter::CovertParToObjects(TVectorD &va,TVectorD &vb,TLorentzVector &Taua1,TLorentzVector &Taumu, double &Zmass){
  // Taua1=particles_.at(0).LV();//TLorentzVector(v(taua1_px),v(taua1_py),v(taua1_pz),sqrt(1.777*1.777+v(taua1_px)*v(taua1_px)+v(taua1_py)*v(taua1_py)+v(taua1_pz)*v(taua1_pz)));
  // Taumu=TLorentzVector(v(taua1_px),v(taua1_py),v(taua1_pz),sqrt(1.777*1.777+v(taua1_px)*v(taua1_px)+v(taua1_py)*v(taua1_py)+v(taua1_pz)*v(taua1_pz)));
  Taua1=TLorentzVector(va(tau_px),va(tau_py),va(tau_pz),sqrt(PDGInfo::tau_mass()*PDGInfo::tau_mass()+va(tau_px)*va(tau_px)+va(tau_py)*va(tau_py)+va(tau_pz)*va(tau_pz)));
  Taumu=TLorentzVector(vb(tau_px),vb(tau_py),vb(tau_pz),sqrt(PDGInfo::tau_mass()*PDGInfo::tau_mass()+vb(tau_px)*vb(tau_px)+vb(tau_py)*vb(tau_py)+vb(tau_pz)*vb(tau_pz)));
  Zmass = 91.5;
}

LorentzVectorParticle
DiTauConstrainedFitter::TauMuStartingPoint(TrackParticle MuTrack,LorentzVectorParticle TauA1, TVector3 PV,TMatrixTSym<double>  PVCov, TVector3 SV, TMatrixTSym<double>  SVCov){

  TVector3 TauDir =   PV - SV;
  TVector3 TauDirError(sqrt(SVCov(0,0) + PVCov(0,0)),
		       sqrt(SVCov(1,1) + PVCov(1,1)),
		       sqrt(SVCov(2,2) + PVCov(2,2)));

  TMatrixT<double>    parameters;
  parameters.ResizeTo(5,1);
  TMatrixT<double>    parameterErrors;
  parameterErrors.ResizeTo(5,5);

  TMatrixT<double>    parametersAd;
  parametersAd.ResizeTo(9,1);
  TMatrixTSym<double>    parameterErrorsAd;
  parameterErrorsAd.ResizeTo(9,9);

  TMatrixT<double>    taumudirection;
  taumudirection.ResizeTo(2,1);
  TMatrixTSym<double>    taumudirectionError;
  taumudirectionError.ResizeTo(2,2);

  TMatrixT<double>    kinematicparameters;
  kinematicparameters.ResizeTo(5,1);
  TMatrixTSym<double>    kinematicparametererrors;
  kinematicparametererrors.ResizeTo(5,5);

  TMatrixT<double>    TauKin;
  TauKin.ResizeTo(3,1);

  TMatrixT<double>    TauKinErrorAnalytical;
  TauKinErrorAnalytical.ResizeTo(3,3);

  TMatrixT<double>    TauKinErrorNumerical;
  TauKinErrorNumerical.ResizeTo(3,3);

  parameters = ConfigureParameters(MuTrack, EstimatePhiAngle(TauDir,TauDirError));
  parameterErrors= ConfigureParameterErrors(MuTrack, EstimatePhiAngle(TauDir,TauDirError));

  parametersAd = ConfigureInitialAdvancedParameters(MuTrack, PV, TauA1);
  parameterErrorsAd= ConfigureInitialAdvancedParameterErrors(MuTrack, PVCov, TauA1.getCovMatrix());

  taumudirection = EstimateTauDirectionAdvanced(parametersAd);
  taumudirectionError = ErrorMatrixPropagator::PropagateError(&DiTauConstrainedFitter::EstimateTauDirectionAdvanced,ConfigureInitialAdvancedParameters(MuTrack, PV, TauA1),ConfigureInitialAdvancedParameterErrors(MuTrack, PVCov, TauA1.getCovMatrix()));

  kinematicparameters = ConfigureKinematicParameters(taumudirection,TauA1);
  kinematicparametererrors =  ConfigureKinematicParameterErrors(taumudirectionError,TauA1);

  TauKin=EstimateTauKinematic(kinematicparameters);
  TauKinErrorNumerical = ErrorMatrixPropagator::PropagateError(&DiTauConstrainedFitter::EstimateTauKinematic,ConfigureKinematicParameters(taumudirection,TauA1),ConfigureKinematicParameterErrors(taumudirectionError,TauA1));


  TauKinErrorAnalytical=ComputeAngleCovarianceAnalytically(MuTrack,EstimatePhiAngle(TauDir,TauDirError),PV,SV,TauA1);

  if(debug){
    std::cout<<"Tau Kinematic Paramteres 'TauKin' ====>"<<std::endl;
    TauKin.Print();

    std::cout<<"Numericaly  calculated tau kinematic parameter errors 'TauKinErrorNumerical'  ====> "<< std::endl;
    TauKinErrorNumerical.Print();

    std::cout<<"Inpute Kinematic Parameters 'kinematicparameters' ====>"<<std::endl;
    kinematicparameters.Print();

    std::cout<<"Input Kinematic Paramter Errors  'kinematicparametererrors' ====>"<< std::endl;
    kinematicparametererrors.Print();

    std::cout<<"TauMu direction paramters 'taumudirection'  ====>  "<< std::endl;
    taumudirection.Print();

    std::cout<<"TauMu direction paramter errors  'taumudirectionError'  ====>  "<< std::endl;
    taumudirectionError.Print();

    std::cout<<"Analitically  calculated tau kinematic parameter errors 'TauKinErrorNumerical'  ====>  "<< std::endl;
    TauKinErrorAnalytical.Print();

  }

  TMatrixT<double>    par(LorentzVectorParticle::NLorentzandVertexPar,1);
  TMatrixTSym<double> Cov(LorentzVectorParticle::NLorentzandVertexPar);

  par(LorentzVectorParticle::vx,0)=0; // fill zero vertex for now
  par(LorentzVectorParticle::vy,0)=0;
  par(LorentzVectorParticle::vz,0)=0;
  par(LorentzVectorParticle::px,0)=TauKin(0,0);
  par(LorentzVectorParticle::py,0)=TauKin(1,0);
  par(LorentzVectorParticle::pz,0)=TauKin(2,0);
  par(LorentzVectorParticle::m,0) =PDGInfo::tau_mass();

   for(int i=0; i<LorentzVectorParticle::NVertex; i++){
     for(int j=0; j<LorentzVectorParticle::NVertex; j++){
       if(AnalyticalCovariance){Cov(i+3,j+3)=TauKinErrorAnalytical(i,j);}
       else{Cov(i+3,j+3)=TauKinErrorNumerical(i,j);}
     }
   }

  return LorentzVectorParticle(par,Cov,PDGInfo::tau_minus,0,0);
}

LorentzVectorParticle
DiTauConstrainedFitter::TauMuStartingPointwithFullRecoil(TrackParticle MuTrack,LorentzVectorParticle TauA1, PTObject METminusNeutrino, TVector3 PV, TMatrixTSym<double>  PVCov, TVector3 SV, TMatrixTSym<double>  SVCov){

  TMatrixT<double>    kinematicparameters;
  kinematicparameters.ResizeTo(12,1);
  TMatrixTSym<double>    kinematicparametererrors;
  kinematicparametererrors.ResizeTo(12,12);

  TMatrixT<double>    TauKin;
  TauKin.ResizeTo(6,1);

  TMatrixT<double>    TauKinErrorNumerical;
  TauKinErrorNumerical.ResizeTo(6,6);
  TMatrixT<double>    TauKinErrorAnalytical;
  TauKinErrorAnalytical.ResizeTo(6,6);
  TMatrixTSym<double>    TauMuErrorAnalytical;
  TauMuErrorAnalytical.ResizeTo(3,3);

  kinematicparameters = ConfigureKinematicParametersFullRecoil(MuTrack, PV, TauA1, METminusNeutrino);
  kinematicparametererrors =  ConfigureKinematicParameterErrorsFullRecoil(MuTrack, PVCov, TauA1, METminusNeutrino);

  TauKin=EstimateTauKinematicFullRecoil(kinematicparameters);
  TauKinErrorNumerical = ErrorMatrixPropagator::PropagateError(&DiTauConstrainedFitter::EstimateTauKinematicFullRecoil,kinematicparameters,kinematicparametererrors);
  TLorentzVector TauMu(TauKin(3,0),TauKin(4,0),TauKin(5,0),PDGInfo::tau_mass());

  TauKinErrorAnalytical = TauKinErrorNumerical;
  TauMuErrorAnalytical = EstimateTauKinematicErrorFullRecoil(TauA1, TauMu, METminusNeutrino);

  if(Logger::Instance()->Level() == Logger::Debug){
	Logger(Logger::Debug) << "TauKinErrorAnalytical: " << std::endl;
	TauKinErrorAnalytical.Print();
	Logger(Logger::Debug) << "TauMuErrorAnalytical: " << std::endl;
	TauMuErrorAnalytical.Print();
  }

  TMatrixT<double>    par(LorentzVectorParticle::NLorentzandVertexPar,1);
  TMatrixTSym<double> Cov(LorentzVectorParticle::NLorentzandVertexPar);

  par(LorentzVectorParticle::vx,0)=PV_.X(); // fill tauh as "vertex" to save correlations
  par(LorentzVectorParticle::vy,0)=PV_.Y();
  par(LorentzVectorParticle::vz,0)=PV_.Z();
  par(LorentzVectorParticle::px,0)=TauKin(3,0);
  par(LorentzVectorParticle::py,0)=TauKin(4,0);
  par(LorentzVectorParticle::pz,0)=TauKin(5,0);
  par(LorentzVectorParticle::m,0) =PDGInfo::tau_mass();

  for(int i=0; i<LorentzVectorParticle::NVertex; i++){
	for(int j=0; j<LorentzVectorParticle::NVertex; j++){
	  Cov(i,j)=PVCov(i,j);
	  Cov(i+3,j+3)=TauMuErrorAnalytical(i,j);
	}
  }
  if(Logger::Instance()->Level() == Logger::Debug){
	Logger(Logger::Debug) << "Cov: " << std::endl;
	Cov.Print();
  }

  return LorentzVectorParticle(par,Cov,PDGInfo::tau_minus,0,0);
}

TMatrixT<double>
DiTauConstrainedFitter::EstimateTauDirectionAdvanced(TMatrixT<double> &inpar){
  TMatrixT<double>    outpar(2,1);

  double dxy   =inpar(0,0);
  double phi0  =inpar(1,0);
  double lam   =inpar(2,0);
  double dz    =inpar(3,0);

  TVector3 PV(inpar(4,0),inpar(5,0),inpar(6,0));
  TVector2 TauMuPt(-inpar(7,0), -inpar(8,0));

  double tanphi_tau = TauMuPt.Y()/TauMuPt.X();
  double a = tanphi_tau;
  double b = PV.Y() - a*PV.X();

  double tmu  = (dxy*(cos(phi0) + a*sin(phi0)) - b) / (a*cos(phi0) - sin(phi0)); //projection onto XY plane

  double xdoc = -dxy*sin(phi0) + tmu*cos(phi0);
  double ydoc = dxy*cos(phi0) + tmu*sin(phi0);
  double zdoc = dz + tmu*tan(lam);

  TVector3 PointGuess(xdoc, ydoc, zdoc);
  TVector3 TauMuDir = PointGuess - PV;

  outpar(0,0) = TauMuDir.Theta();
  outpar(1,0) = atan2(TauMuPt.Y(),TauMuPt.X());

  Logger(Logger::Debug) << "PV: " << inpar(4,0) << ", " << inpar(5,0) << ", " << inpar(6,0) << std::endl;
  Logger(Logger::Debug) << "dxy, phi0, lam, dz: " << dxy << ", " << phi0 << ", " << lam << ", " << dz << std::endl;
  Logger(Logger::Debug) << "PointGuess: " << xdoc << ", " << ydoc << ", " << zdoc << std::endl;
  Logger(Logger::Debug) << "DirGuess: " << TauMuDir.X() << ", " << TauMuDir.Y() << ", " << TauMuDir.Z() << std::endl;
  Logger(Logger::Debug) << "a, b, t: " << a << ", " << b << ", " << tmu <<  std::endl;
  Logger(Logger::Debug) << "TauH X, Y: " << inpar(7,0) << ", " << inpar(8,0) << std::endl;
  Logger(Logger::Debug) << "TauMuDir.Phi(): " << TauMuDir.Phi() << " TauH phi: " << atan2(inpar(8,0), inpar(7,0)) << std::endl;
  Logger(Logger::Debug) << "dPhi: " << TauMuDir.Phi() - atan2(inpar(8,0), inpar(7,0)) << "/pi= " << (TauMuDir.Phi() - atan2(inpar(8,0), inpar(7,0)))/TMath::Pi() << std::endl;

  return outpar;
}


TMatrixT<double>
DiTauConstrainedFitter::EstimateTauKinematic(TMatrixT<double> &inpar){
  TMatrixT<double>    outpar(3,1);
  TLorentzVector TauA1p4(inpar(2,0),
			 inpar(3,0),
			 inpar(4,0),
			 sqrt(pow(PDGInfo::tau_mass(),2) + pow(inpar(2,0),2) + pow(inpar(3,0),2) + pow(inpar(4,0),2)  ) );

  TVector3 TauMuDir(cos(inpar(1,0))*sin(inpar(0,0)), sin(inpar(1,0))*sin(inpar(0,0)), cos(inpar(0,0)));
  TVector3 P_Tauh(inpar(2,0), inpar(3,0), inpar(4,0));
  TLorentzVector P4_Tauh; P4_Tauh.SetXYZM(P_Tauh.X(),P_Tauh.Y(),P_Tauh.Z(),PDGInfo::tau_mass());
  double theta = P_Tauh.Angle(TauMuDir);
  double MassDiffsq = pow(MassConstraint_, 2.) - pow(PDGInfo::tau_mass(), 2.);
  double Denominator = pow(P4_Tauh.P()*sin(theta),2.) + pow(PDGInfo::tau_mass(),2.);
  double P_TauMu = (MassDiffsq*P4_Tauh.P()*cos(theta) + P4_Tauh.E()*sqrt( pow(MassDiffsq, 2.) - 4*pow(PDGInfo::tau_mass(),2.)*Denominator))/2/(Denominator);

  outpar(0,0) = P_TauMu*TauMuDir.X();
  outpar(1,0) = P_TauMu*TauMuDir.Y();
  outpar(2,0) = P_TauMu*TauMuDir.Z();

  //outpar(0,0) = P_Tauh.Pt()*cos(inpar(1,0));
  //outpar(1,0) = P_Tauh.Pt()*sin(inpar(1,0));
  //outpar(2,0) = P_Tauh.Pt()/tan(inpar(0,0));

  Logger(Logger::Debug) << "TauDir.Phi(): " << TauA1p4.Phi() << " TauMuDirNEW2.Phi(): " << inpar(1,0) << std::endl;
  Logger(Logger::Debug) << "TauA1 p3: " << TauA1p4.X() << ", " << TauA1p4.Y() << ", " << TauA1p4.Z() << std::endl;
  Logger(Logger::Debug) << "TauMu p3: " << outpar(0,0) << ", " << outpar(1,0) << ", " << outpar(2,0) << std::endl;


  return outpar;
}


TMatrixT<double>
DiTauConstrainedFitter::ConfigureParameters(TrackParticle MuTrack, std::pair<double, double> phiAngle){
   TMatrixT<double>    outpar(5,1);
   outpar(0,0) = MuTrack.Parameter(TrackParticle::dxy);
   outpar(1,0) = MuTrack.Parameter(TrackParticle::phi);
   outpar(2,0) = MuTrack.Parameter(TrackParticle::lambda);
   outpar(3,0) = MuTrack.Parameter(TrackParticle::dz);
   outpar(4,0) = phiAngle.first;

   return outpar;
 }

TMatrixT<double>
DiTauConstrainedFitter::ConfigureParameterErrors(TrackParticle MuTrack, std::pair<double, double> phiAngle){
    TMatrixT<double>  Cov;
    Cov.ResizeTo(5,5);

    Cov(0,0) = MuTrack.Covariance(TrackParticle::dxy,TrackParticle::dxy);
    Cov(0,1) = MuTrack.Covariance(TrackParticle::dxy,TrackParticle::phi);
    Cov(0,2) = MuTrack.Covariance(TrackParticle::dxy,TrackParticle::dz);
    Cov(0,3) = MuTrack.Covariance(TrackParticle::dxy,TrackParticle::lambda);
    Cov(0,4) = 0;

    Cov(1,0) = MuTrack.Covariance(TrackParticle::phi,TrackParticle::dxy);
    Cov(1,1) = MuTrack.Covariance(TrackParticle::phi,TrackParticle::phi);
    Cov(1,2) = MuTrack.Covariance(TrackParticle::phi,TrackParticle::dz);
    Cov(1,3) = MuTrack.Covariance(TrackParticle::phi,TrackParticle::lambda);
    Cov(1,4) = 0;

    Cov(2,0) = MuTrack.Covariance(TrackParticle::dz,TrackParticle::dxy);
    Cov(2,1) = MuTrack.Covariance(TrackParticle::dz,TrackParticle::phi);
    Cov(2,2) = MuTrack.Covariance(TrackParticle::dz,TrackParticle::dz);
    Cov(2,3) = MuTrack.Covariance(TrackParticle::dz,TrackParticle::lambda);
    Cov(2,4) = 0;

    Cov(3,0) = MuTrack.Covariance(TrackParticle::lambda,TrackParticle::dxy);
    Cov(3,1) = MuTrack.Covariance(TrackParticle::lambda,TrackParticle::phi);
    Cov(3,2) = MuTrack.Covariance(TrackParticle::lambda,TrackParticle::dz);
    Cov(3,3) = MuTrack.Covariance(TrackParticle::lambda,TrackParticle::lambda);
    Cov(3,4) = 0;

    Cov(4,0) = 0;
    Cov(4,1) = 0;
    Cov(4,2) = 0;
    Cov(4,3) = 0;
    Cov(4,4) = phiAngle.second;

    return Cov;
}


TMatrixT<double>
DiTauConstrainedFitter::ConfigureInitialAdvancedParameters(TrackParticle MuTrack,  TVector3 PV, LorentzVectorParticle TauA1){
  TMatrixT<double>    outpar;
  outpar.ResizeTo(9,1);
   outpar(0,0)  = MuTrack.Parameter(TrackParticle::dxy);
   outpar(1,0)  = MuTrack.Parameter(TrackParticle::phi);
   outpar(2,0)  = MuTrack.Parameter(TrackParticle::lambda);
   outpar(3,0)  = MuTrack.Parameter(TrackParticle::dz);
   outpar(4,0)  = PV.X();
   outpar(5,0)  = PV.Y();
   outpar(6,0)  = PV.Z();
   outpar(7,0) = TauA1.LV().X();
   outpar(8,0) = TauA1.LV().Y();

   return outpar;
 }

TMatrixTSym<double>
DiTauConstrainedFitter::ConfigureInitialAdvancedParameterErrors(TrackParticle MuTrack, TMatrixTSym<double> PVCov, TMatrixTSym<double> TauA1Cov){
  TMatrixTSym<double>  Cov;
    Cov.ResizeTo(9,9);

    Cov(0,0) = MuTrack.Covariance(TrackParticle::dxy,TrackParticle::dxy);
    Cov(0,1) = MuTrack.Covariance(TrackParticle::dxy,TrackParticle::phi);
    Cov(0,2) = MuTrack.Covariance(TrackParticle::dxy,TrackParticle::dz);
    Cov(0,3) = MuTrack.Covariance(TrackParticle::dxy,TrackParticle::lambda);

    Cov(1,0) = MuTrack.Covariance(TrackParticle::phi,TrackParticle::dxy);
    Cov(1,1) = MuTrack.Covariance(TrackParticle::phi,TrackParticle::phi);
    Cov(1,2) = MuTrack.Covariance(TrackParticle::phi,TrackParticle::dz);
    Cov(1,3) = MuTrack.Covariance(TrackParticle::phi,TrackParticle::lambda);

    Cov(2,0) = MuTrack.Covariance(TrackParticle::lambda,TrackParticle::dxy);
    Cov(2,1) = MuTrack.Covariance(TrackParticle::lambda,TrackParticle::phi);
    Cov(2,2) = MuTrack.Covariance(TrackParticle::lambda,TrackParticle::dz);
    Cov(2,3) = MuTrack.Covariance(TrackParticle::lambda,TrackParticle::lambda);

    Cov(3,0) = MuTrack.Covariance(TrackParticle::dz,TrackParticle::dxy);
    Cov(3,1) = MuTrack.Covariance(TrackParticle::dz,TrackParticle::phi);
    Cov(3,2) = MuTrack.Covariance(TrackParticle::dz,TrackParticle::dz);
    Cov(3,3) = MuTrack.Covariance(TrackParticle::dz,TrackParticle::lambda);

    for(int i=0; i<3; i++){
      for(int j=0; j<3; j++){
	Cov(i+4,j+4)=PVCov(i,j);
      }
    }
    for(int i=0; i<2; i++){
      for(int j=0; j<2; j++){
	Cov(i+7,j+7)=TauA1Cov(i+4,j+4);
      }
    }
    if(debug){
      std::cout<<"  Input paramters covariance matrix "<<std::endl;
      for(int str =0; str < Cov.GetNrows(); str++){
	for(int kol =0; kol < Cov.GetNcols(); kol++){
	  std::cout<<"  "<< Cov(str,kol)<<"  ";
	}
	std::cout<<std::endl;
      }
    }

    return Cov;
 }

TMatrixT<double>
DiTauConstrainedFitter::ConfigureKinematicParameters(TMatrixT<double>  TauMuDir, LorentzVectorParticle TauA1){
  TMatrixT<double>    outpar;
  outpar.ResizeTo(5,1);
  outpar(0,0)  = TauMuDir(0,0);
  outpar(1,0)  = TauMuDir(1,0);
  outpar(2,0)  = TauA1.LV().Px();
  outpar(3,0)  = TauA1.LV().Py();
  outpar(4,0)  = TauA1.LV().Pz();


  if(debug){
    std::cout<<"EstimateTauKinematic  TauA1  "<<TauA1.LV().Px()<<"  "<<TauA1.LV().Py()<<"  "<<TauA1.LV().Pz()<<std::endl;
    std::cout<<"ConfigureKinematicParameters:  outpar(0,0)  "<<outpar(0,0) <<std::endl;
    std::cout<<"ConfigureKinematicParameters:  outpar(1,0)  "<<outpar(1,0) <<std::endl;
    std::cout<<"ConfigureKinematicParameters:  outpar(2,0)  "<<outpar(2,0) <<std::endl;
    std::cout<<"ConfigureKinematicParameters:  outpar(3,0)  "<<outpar(3,0) <<std::endl;
    std::cout<<"ConfigureKinematicParameters:  outpar(4,0)  "<<outpar(4,0) <<std::endl;
  }
  return outpar;
 }

TMatrixTSym<double>
DiTauConstrainedFitter::ConfigureKinematicParameterErrors(TMatrixTSym<double>  TauMuDirError, LorentzVectorParticle TauA1){
  TMatrixTSym<double>    outpar;
  outpar.ResizeTo(5,5);

  for(int i=0; i<5; i++){
    for(int j=0; j<5; j++){

      if(i<2 && j < 2){
	outpar(i,j)=TauMuDirError(i,j);
      }else{
	outpar(i,j)=TauA1.Covariance(i+1,j+1);
      }
    }
  }

  return outpar;
}

std::pair<double, double>
DiTauConstrainedFitter::EstimatePhiAngle( TVector3 dir, TVector3 dirE){
    std::pair<double, double> outpar;
  double phi = atan2(dir.Y(),dir.X());
  double deltaphi = sqrt(dirE.Y()*dirE.Y()/dir.X()/dir.X()+
 			 dirE.X()*dirE.X()*dir.Y()*dir.Y()/dir.X()/dir.X()/dir.X()/dir.X() )*cos(phi)*cos(phi);
  outpar = std::make_pair(phi, deltaphi);
  return outpar;
}

//////////////////////////////////////////////////////////////////////
// Analytical calculation of TauMu Covariance; Has to be rewritten....
TMatrixT<double>
DiTauConstrainedFitter::ComputeAngleCovarianceAnalytically(TrackParticle MuTrack, std::pair<double, double> phiAngle,  TVector3 PV, TVector3 SV, LorentzVectorParticle  TauA1){
   TMatrixT<double>    outpar;
   outpar.ResizeTo(3,3);
   double dxy   =MuTrack.Parameter(TrackParticle::dxy);
   //double kappa =MuTrack.Parameter(TrackParticle::kappa);
   double phi0  =MuTrack.Parameter(TrackParticle::phi);
   double lam   =MuTrack.Parameter(TrackParticle::lambda);
   double dz    =MuTrack.Parameter(TrackParticle::dz);
   //double c     =MuTrack.Charge();
   TVector3 A1SV = -SV + PV;

   double phiAnot  = atan2(A1SV.Y(), A1SV.X());
   double xpoca2Anot = dxy*sin(phi0) - PV.X();
   double ypoca2Anot = -dxy*cos(phi0) - PV.Y();
   double aAnot = tan(phi0);
   double bAnot = ypoca2Anot - aAnot*xpoca2Anot;
   double r = sqrt( pow(bAnot/(tan(phiAnot) - aAnot )  ,2) + pow(bAnot*tan(phiAnot)/(tan(phiAnot) - aAnot) ,2));//(bAnot)/(tan(phiAnot) - tan(phi0))/cos(phiAnot);


   //double bz0 = fabs(dxy) - dz/lam ;
   double ZNeu = lam*(r +  dxy)  - dz;
   //double XNeu = r*cos(phiAnot);
   //double YNeu = r*sin(phiAnot);


   double br = tan(phi0)*cos(phiAnot) - sin(phiAnot);
   double Zc = ZNeu;

   double MinusSintheta = -sqrt(1  - Zc*Zc/(r*r + Zc*Zc) );

   double cosTheta = Zc/sqrt(r*r + Zc*Zc) ;
   double sinTheta = -MinusSintheta;
   if(fabs(phiAnot - TauA1.LV().Phi()) < 0.05) phiAnot = phiAnot - TMath::Pi();
   //-derivaticves
   double drdphi =2*dxy*sin(phi0)*(cos(phiAnot) + tan(phi0)*sin(phiAnot))/ (tan(phi0)*cos(phiAnot) - sin(phiAnot))/(tan(phi0)*cos(phiAnot) - sin(phiAnot));
   double drzdd = 2*sin(phi0)/br/Zc  - 2*dxy*sin(phi0)*lam/br/Zc/Zc;
   double drzdphi0 = (2*dxy*cos(phi0)/(tan(phi0)*cos(phiAnot) - sin(phiAnot))   - 2*dxy*sin(phi0)*cos(phiAnot)/cos(phi0)/cos(phi0)/(tan(phi0)*cos(phiAnot) - sin(phiAnot))/(tan(phi0)*cos(phiAnot) - sin(phiAnot)))/Zc;
   double drzdlam = -2*dxy*sin(phi0)*(r+dxy)/br/Zc/Zc;
   double drzdz0 = -2*dxy*sin(phi0)/br/Zc;
   double drzdphi = drdphi*(1 - r*lam/Zc)/Zc;//2*dxy*sin(phi0)*(cos(phiAnot) + tan(phi0)*sin(phiAnot))/ (tan(phi0)*cos(phiAnot) - sin(phiAnot))/(tan(phi0)*cos(phiAnot) - sin(phiAnot));
   double dcosThetadrz = -r/Zc/sqrt(pow(1 + r*r/Zc/Zc,3));
   //double dThetadPhi = -dcosThetadrz*drzdphi/sinTheta;
   //-derivaticves

   double cosTauTau2 = cos(TauA1.LV().Theta() + acos(cosTheta));
   //double sinTauTau2 = sin(TauA1.LV().Theta() + acos(cosTheta));

   TMatrixT<double>  HelixCov;
   HelixCov.ResizeTo(5,5);

   HelixCov(0,0) = MuTrack.Covariance(TrackParticle::dxy,TrackParticle::dxy);
   HelixCov(0,1) = MuTrack.Covariance(TrackParticle::dxy,TrackParticle::phi);
   HelixCov(0,2) = MuTrack.Covariance(TrackParticle::dxy,TrackParticle::dz);
   HelixCov(0,3) = MuTrack.Covariance(TrackParticle::dxy,TrackParticle::lambda);
   HelixCov(0,4) = 0;

   HelixCov(1,0) = MuTrack.Covariance(TrackParticle::phi,TrackParticle::dxy);
   HelixCov(1,1) = MuTrack.Covariance(TrackParticle::phi,TrackParticle::phi);
   HelixCov(1,2) = MuTrack.Covariance(TrackParticle::phi,TrackParticle::dz);
   HelixCov(1,3) = MuTrack.Covariance(TrackParticle::phi,TrackParticle::lambda);
   HelixCov(1,4) = 0;

   HelixCov(2,0) = MuTrack.Covariance(TrackParticle::dz,TrackParticle::dxy);
   HelixCov(2,1) = MuTrack.Covariance(TrackParticle::dz,TrackParticle::phi);
   HelixCov(2,2) = MuTrack.Covariance(TrackParticle::dz,TrackParticle::dz);
   HelixCov(2,3) = MuTrack.Covariance(TrackParticle::dz,TrackParticle::lambda);
   HelixCov(2,4) = 0;

   HelixCov(3,0) = MuTrack.Covariance(TrackParticle::lambda,TrackParticle::dxy);
   HelixCov(3,1) = MuTrack.Covariance(TrackParticle::lambda,TrackParticle::phi);
   HelixCov(3,2) = MuTrack.Covariance(TrackParticle::lambda,TrackParticle::dz);
   HelixCov(3,3) = MuTrack.Covariance(TrackParticle::lambda,TrackParticle::lambda);
   HelixCov(3,4) = 0;

   HelixCov(4,0) = 0;
   HelixCov(4,1) = 0;
   HelixCov(4,2) = 0;
   HelixCov(4,3) = 0;
   HelixCov(4,4) = phiAngle.second;//sqrt(TauDirError.Y()*TauDirError.Y()/A1SV.X()/A1SV.X()   + TauDirError.X()*TauDirError.X()*A1SV.Y()*A1SV.Y()/A1SV.X()/A1SV.X()/A1SV.X()/A1SV.X() ) *cos(phiAnot)*cos(phiAnot);

   TMatrixT<double>  DerivativesHelixToAngles;
   DerivativesHelixToAngles.ResizeTo(2,5);

   DerivativesHelixToAngles(0,0)  = 0;
   DerivativesHelixToAngles(0,1)  = 0;
   DerivativesHelixToAngles(0,2)  = 0;
   DerivativesHelixToAngles(0,3)  = 0;
   DerivativesHelixToAngles(0,4)  = 1;

   DerivativesHelixToAngles(1,0)  = dcosThetadrz*drzdd/MinusSintheta;
   DerivativesHelixToAngles(1,1)  = dcosThetadrz*drzdphi0/MinusSintheta;
   DerivativesHelixToAngles(1,2)  = dcosThetadrz*drzdlam/MinusSintheta;
   DerivativesHelixToAngles(1,3)  = dcosThetadrz*drzdz0/MinusSintheta;
   DerivativesHelixToAngles(1,4)  = dcosThetadrz*drzdphi/MinusSintheta;

   TMatrixT<double> DerivativesHelixToAnglesT=DerivativesHelixToAngles; DerivativesHelixToAnglesT.T();
   TMatrixT<double> CovAngleFrame=DerivativesHelixToAngles*HelixCov*DerivativesHelixToAnglesT;

   //double ZMassR = 91.5;
   double TauA1deltaP = sqrt(   (pow(TauA1.LV().Px(),2)*fabs(TauA1.Covariance(3,3))   +  pow(TauA1.LV().Py(),2)*fabs(TauA1.Covariance(4,4))  +  pow(TauA1.LV().Pz(),2)*fabs(TauA1.Covariance(5,5))    )/TauA1.LV().P()/TauA1.LV().P()  )  ;
   double TauMudeltaP_2 = TauA1deltaP*pow(MassConstraint_/2/TauA1.LV().P(),2)/(1-cosTauTau2);

   double TauMuP_2 =  MassConstraint_*MassConstraint_/2/(1-cosTauTau2)/TauA1.LV().P();


   TMatrixT<double> CovPPhiThetaFrame;
   CovPPhiThetaFrame.ResizeTo(3,3);
   //---------- dp dphi dtheta
   CovPPhiThetaFrame(0,0) = TauMudeltaP_2;
   CovPPhiThetaFrame(0,1) = 0;
   CovPPhiThetaFrame(0,2) = 0;

   CovPPhiThetaFrame(1,0) = 0;
   CovPPhiThetaFrame(1,1) = CovAngleFrame(0,0);
   CovPPhiThetaFrame(1,2) = CovAngleFrame(0,1);

   CovPPhiThetaFrame(2,0) = 0;
   CovPPhiThetaFrame(2,1) = CovAngleFrame(1,0);
   CovPPhiThetaFrame(2,2) = CovAngleFrame(1,1);


   TMatrixT<double> DirevativesPPhiThetaToPxPyPz;
   DirevativesPPhiThetaToPxPyPz.ResizeTo(3,3);

   DirevativesPPhiThetaToPxPyPz(0,0) =  cos(phiAnot)*sinTheta;
   DirevativesPPhiThetaToPxPyPz(0,1) = -TauMuP_2*sin(phiAnot)*sinTheta;
   DirevativesPPhiThetaToPxPyPz(0,2) =  TauMuP_2*cos(phiAnot)*cosTheta;

   DirevativesPPhiThetaToPxPyPz(1,0) = sin(phiAnot)*sinTheta;
   DirevativesPPhiThetaToPxPyPz(1,1) = TauMuP_2*cos(phiAnot)*sinTheta;
   DirevativesPPhiThetaToPxPyPz(1,2) = TauMuP_2*sin(phiAnot)*cosTheta;

   DirevativesPPhiThetaToPxPyPz(2,0) = cosTheta;
   DirevativesPPhiThetaToPxPyPz(2,1) = 0;
   DirevativesPPhiThetaToPxPyPz(2,2) =-TauMuP_2*sinTheta;


   TMatrixT<double> DirevativesPPhiThetaToPxPyPzT=DirevativesPPhiThetaToPxPyPz; DirevativesPPhiThetaToPxPyPzT.T();
   TMatrixT<double> CovPxPyPzFrame=DirevativesPPhiThetaToPxPyPz*CovPPhiThetaFrame*DirevativesPPhiThetaToPxPyPzT;


   outpar = CovPxPyPzFrame;
   return outpar;
 }


LorentzVectorParticle
DiTauConstrainedFitter::GetTauMuEstimate(){
  return particles_.at(1);
}


TMatrixT<double> DiTauConstrainedFitter::ConfigureKinematicParametersFullRecoil(TrackParticle MuTrack, TVector3 PV, LorentzVectorParticle TauA1, PTObject METMinusNeutrino){
  TMatrixT<double> outpar;
  outpar.ResizeTo(12,1);

  outpar(0,0)  = MuTrack.Parameter(TrackParticle::dxy);
  outpar(1,0)  = MuTrack.Parameter(TrackParticle::phi);
  outpar(2,0)  = MuTrack.Parameter(TrackParticle::lambda);
  outpar(3,0)  = MuTrack.Parameter(TrackParticle::dz);
  outpar(4,0)  = PV.X();
  outpar(5,0)  = PV.Y();
  outpar(6,0)  = PV.Z();
  outpar(7,0)  = TauA1.LV().X();
  outpar(8,0)  = TauA1.LV().Y();
  outpar(9,0)  = TauA1.LV().Z();
  outpar(10,0) = METMinusNeutrino.Par()(0,0);
  outpar(11,0) = METMinusNeutrino.Par()(1,0);

  return outpar;
}

TMatrixTSym<double> DiTauConstrainedFitter::ConfigureKinematicParameterErrorsFullRecoil(TrackParticle MuTrack, TMatrixTSym<double> PVCov, LorentzVectorParticle TauA1, PTObject METMinusNeutrino){
  TMatrixTSym<double>  Cov;
    Cov.ResizeTo(12,12);

    Cov(0,0) = MuTrack.Covariance(TrackParticle::dxy,TrackParticle::dxy);
    Cov(0,1) = MuTrack.Covariance(TrackParticle::dxy,TrackParticle::phi);
    Cov(0,2) = MuTrack.Covariance(TrackParticle::dxy,TrackParticle::dz);
    Cov(0,3) = MuTrack.Covariance(TrackParticle::dxy,TrackParticle::lambda);

    Cov(1,0) = MuTrack.Covariance(TrackParticle::phi,TrackParticle::dxy);
    Cov(1,1) = MuTrack.Covariance(TrackParticle::phi,TrackParticle::phi);
    Cov(1,2) = MuTrack.Covariance(TrackParticle::phi,TrackParticle::dz);
    Cov(1,3) = MuTrack.Covariance(TrackParticle::phi,TrackParticle::lambda);

    Cov(2,0) = MuTrack.Covariance(TrackParticle::lambda,TrackParticle::dxy);
    Cov(2,1) = MuTrack.Covariance(TrackParticle::lambda,TrackParticle::phi);
    Cov(2,2) = MuTrack.Covariance(TrackParticle::lambda,TrackParticle::dz);
    Cov(2,3) = MuTrack.Covariance(TrackParticle::lambda,TrackParticle::lambda);

    Cov(3,0) = MuTrack.Covariance(TrackParticle::dz,TrackParticle::dxy);
    Cov(3,1) = MuTrack.Covariance(TrackParticle::dz,TrackParticle::phi);
    Cov(3,2) = MuTrack.Covariance(TrackParticle::dz,TrackParticle::dz);
    Cov(3,3) = MuTrack.Covariance(TrackParticle::dz,TrackParticle::lambda);

    for(int i=0; i<3; i++){
      for(int j=0; j<3; j++){
        Cov(i+4,j+4)=PVCov(i,j);
        Cov(i+7,j+7)=TauA1.Covariance(i+4,j+4);
      }
    }
    for(int i=0; i<2; i++){
      for(int j=0; j<2; j++){
        Cov(i+10,j+10)=METMinusNeutrino.Cov()(i,j);
      }
    }
    return Cov;
}

TMatrixT<double> DiTauConstrainedFitter::EstimateTauKinematicFullRecoil(TMatrixT<double> &inpar){
  TMatrixT<double> outpar;
  outpar.ResizeTo(6,1);

  double dxy   =inpar(0,0);
  double phi0  =inpar(1,0);
  double lam   =inpar(2,0);
  double dz    =inpar(3,0);

  TVector3 PV(inpar(4,0),inpar(5,0),inpar(6,0));
  TVector3 P_Tauh(inpar(7,0), inpar(8,0), inpar(9,0));
  TLorentzVector P4_Tauh; P4_Tauh.SetXYZM(P_Tauh.X(),P_Tauh.Y(),P_Tauh.Z(),PDGInfo::tau_mass());
  TVector2 TauMuPt(inpar(10,0) - P_Tauh.X(), inpar(11,0) - P_Tauh.Y());
  TVector3 P_Taumu(inpar(10,0) - P_Tauh.X(), inpar(11,0) - P_Tauh.Y(), 0);

  double tanphi_tau = TauMuPt.Y()/TauMuPt.X();
  double a = tanphi_tau;
  double b = PV.Y() - a*PV.X();

  double tmu  = (dxy*(cos(phi0) + a*sin(phi0)) - b) / (a*cos(phi0) - sin(phi0)); //projection onto XY plane

  double xdoc = -dxy*sin(phi0) + tmu*cos(phi0);
  double ydoc = dxy*cos(phi0) + tmu*sin(phi0);
  double zdoc = dz + tmu*tan(lam);

  TVector3 PointGuess(xdoc, ydoc, zdoc);
  TVector3 TauMuDir = PointGuess - PV;

  double phitau = TauMuDir.Phi();
  double thetatau = TauMuDir.Theta();

  outpar(0,0) = P4_Tauh.X();
  outpar(1,0) = P4_Tauh.Y();
  outpar(2,0) = P4_Tauh.Z();

  if(useCollinearityTauMu_){
    outpar(3,0) = TauMuPt.X();
    outpar(4,0) = TauMuPt.Y();
    P_Taumu.SetZ(TauMuPt.Mod()*tan(lam));
    outpar(5,0) = P_Taumu.Z();
  }
  else{
    outpar(3,0) = TauMuPt.X();
    outpar(4,0) = TauMuPt.Y();
    P_Taumu.SetZ(TauMuPt.Mod()/tan(thetatau));
    outpar(5,0) = P_Taumu.Z();
  }

  return outpar;
}

TMatrixDSym DiTauConstrainedFitter::EstimateTauKinematicErrorFullRecoil(LorentzVectorParticle TauA1, TLorentzVector TauMu, PTObject ResPtEstimate){
  double TauMuPt = sqrt(pow(ResPtEstimate.X() - TauA1.LV().X(), 2.) + pow(ResPtEstimate.Y() - TauA1.LV().Y(), 2.));
  TMatrixDSym TauA1_ZPt_Cov(4);
  for(int i_row = 0; i_row<ResPtEstimate.Cov().GetNrows(); i_row++){
	for(int i_col = 0; i_col<ResPtEstimate.Cov().GetNcols(); i_col++){
	  TauA1_ZPt_Cov(i_row,i_col) = TauA1.Covariance(i_row + LorentzVectorParticle::px, i_col + LorentzVectorParticle::px);
	  TauA1_ZPt_Cov(i_row+2,i_col+2) = ResPtEstimate.Cov()(i_row, i_col);
	}
  }
  TMatrixD Jacobi(3,4);
  Jacobi(0,0) = -1.;
  Jacobi(0,2) = +1.;
  Jacobi(1,1) = -1.;
  Jacobi(1,3) = +1.;
  Jacobi(2,0) = -TauA1.LV().X()/TauMuPt*tan(TauMu.Theta());
  Jacobi(2,1) = +ResPtEstimate.X()/TauMuPt*tan(TauMu.Theta());
  Jacobi(2,2) = -TauA1.LV().Y()/TauMuPt*tan(TauMu.Theta());
  Jacobi(2,3) = +ResPtEstimate.Y()/TauMuPt*tan(TauMu.Theta());

  TMatrixDSym TauMuCov(TauA1_ZPt_Cov);
  TauMuCov.Similarity(Jacobi);

  if(Logger::Instance()->Level() == Logger::Debug){
	Logger(Logger::Debug) << "TauA1_ZPt_Cov: " << std::endl;
	TauA1_ZPt_Cov.Print();
	Logger(Logger::Debug) << "Jacobi: " << std::endl;
	Jacobi.Print();
	Logger(Logger::Debug) << "TauMuCov: " << std::endl;
	TauMuCov.Print();
  }

  return TauMuCov;
}

double DiTauConstrainedFitter::CosThetaTauMu(TLorentzVector TauMu){
  double dxy   = MuTrack_.Parameter(TrackParticle::dxy);
  double phi0  = MuTrack_.Parameter(TrackParticle::phi);
  double lam   = MuTrack_.Parameter(TrackParticle::lambda);
  double dz    = MuTrack_.Parameter(TrackParticle::dz);
  double tanphi_tau = TauMu.Y()/TauMu.X();
  double a = tanphi_tau;
  double b = PV_.Y() - a*PV_.X();

  double tmu  = (dxy*(cos(phi0) + a*sin(phi0)) - b) / (a*cos(phi0) - sin(phi0)); //projection onto XY plane

  double xdoc = -dxy*sin(phi0) + tmu*cos(phi0);
  double ydoc = dxy*cos(phi0) + tmu*sin(phi0);
  double zdoc = dz + tmu*tan(lam);

  TVector3 PointGuess(xdoc, ydoc, zdoc);
  TVector3 TauMuDir = PointGuess - PV_;

  if(useCollinearityTauMu_){
    return sin(lam);
  }
  else{
    return TauMuDir.CosTheta();
  }
}

bool DiTauConstrainedFitter::isConverged(){
  if(!useFullRecoil_){
//	if(pardelta<MaxParDelta_ /*&& harddelta_vec.Norm1() < MaxHCDelta_ && chi2prev-chi2 < MaxChi2Delta_  && chi2prev>chi2*/){
		if(/*pardelta<MaxParDelta_ &&*/ harddelta_vec.Norm1() < MaxHCDelta_ && fabs(chi2prev-chi2) < MaxChi2Delta_ /*&& softdelta_vec.Norm1() < MaxSCDelta_ && chi2prev>chi2*/
      && (pow(para(0), 2.0) + pow(para(1), 2.0)) > A1_.LV().Perp2()
      && (pow(parb(0), 2.0) + pow(parb(1), 2.0)) > pow(MuTrack_.Pt(), 2.0)
    ){
    //	Logger(Logger::Verbose) << "converged " << delta << " chi2 " <<  chi2 << " chi2prev " << chi2prev <<"  Maxdelta  " <<MaxDelta_ <<std::endl;
	  return true;
	}
  }
  else{
	//if(pardelta<MaxParDelta_ && harddelta_vec.Norm1() < MaxHCDelta_ && softdelta_vec.Norm1() < MaxSCDelta_){
	// if(/*pardelta<MaxParDelta_ &&*/ harddelta_vec.Norm1() < MaxHCDelta_ && fabs(chi2prev-chi2) < MaxChi2Delta_ /*&& chi2prev>chi2*/){
	if(/*pardelta<MaxParDelta_ && harddelta_vec.Norm1() < MaxHCDelta_ && */fabs(chi2prev-chi2) < MaxChi2Delta_ /*&& chi2prev>chi2*/
      && (pow(para(0), 2.0) + pow(para(1), 2.0)) > A1_.LV().Perp2()
      && (pow(parb(0), 2.0) + pow(parb(1), 2.0)) > pow(MuTrack_.Pt(), 2.0)
    ){
	  return true;
	}
	else{
	  Logger(Logger::Debug) << "Fit did not converge, because: " << std::endl;
	  if(pardelta > MaxParDelta_) Logger(Logger::Debug) << "pardelta  = " << pardelta <<" > MaxParDelta_ = " << MaxParDelta_ << std::endl;
	  if(chi2prev - chi2 > MaxChi2Delta_) Logger(Logger::Debug) << "chi2prev - chi2  = "<< chi2prev - chi2 << " > " << MaxChi2Delta_ << std::endl;
	  if(harddelta_vec.Sum() >= MaxHCDelta_) Logger(Logger::Debug) << "harddelta_vec.Norm1() = " << harddelta_vec.Norm1() <<" >= MaxHCDelta_ = " << MaxHCDelta_ << std::endl;
	  if(softdelta_vec.Sum() >= MaxSCDelta_) Logger(Logger::Debug) << "softdelta_vec.Norm1() = " << softdelta_vec.Norm1() <<" >= MaxSCDelta_ = " << MaxSCDelta_ << std::endl;
	}
  }
  // if(delta<MaxDelta_ /*&& chi2prev-chi2<1.0 && chi2prev>chi2*/){
  //   std::cout << "converged " << delta << " chi2 " <<  chi2 << " chi2prev " << chi2prev <<"  Maxdelta  " <<MaxDelta_ <<std::endl;

  //   return true;
  // }
   return false;
   // return true;
}

TString DiTauConstrainedFitter::ParName(int par){
  switch (par){
    case taua1_px: return "taua1_px";
    case taua1_py: return "taua1_py";
    case taua1_pz: return "taua1_pz";
    case taumu_px: return "taumu_px";
    case taumu_py: return "taumu_py";
    case taumu_pz: return "taumu_pz";
    case lambda_1: return "lambda_1";
    case lambda_2: return "lambda_2";
    // add more/remove lambdas if number of hard constraints exceeds number of lambdas or vice versa
    default:
      Logger(Logger::Warning) << "Paramater Index " << par << " out of bounds !!";
      return "";
  }
}