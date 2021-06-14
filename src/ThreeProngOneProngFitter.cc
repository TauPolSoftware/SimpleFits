#include "TauPolSoftware/SimpleFits/interface/ThreeProngOneProngFitter.h"
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

double ThreeProngOneProngFitter::MassConstraint_ = 125.2;
bool ThreeProngOneProngFitter::useCollinearityTauOneProng_ = false;

ThreeProngOneProngFitter::ThreeProngOneProngFitter(LorentzVectorParticle TauThreeProng, LorentzVectorParticle ThreeProng, LorentzVectorParticle OneProng, TrackParticle OneProngTrack, PTObject ResPtEstimate, TVector3 PVertex, TMatrixTSym<double> VertexCov):
  LagrangeMultipliersFitter()
{
  ThreeProngOneProngFitter(TauThreeProng, ThreeProng, OneProng, OneProngTrack, ResPtEstimate, PVertex, VertexCov, ThreeProngOneProngFitter::MassConstraint_);
}

ThreeProngOneProngFitter::ThreeProngOneProngFitter(LorentzVectorParticle TauThreeProng, LorentzVectorParticle ThreeProng, LorentzVectorParticle OneProng, TrackParticle OneProngTrack, PTObject ResPtEstimate, TVector3 PVertex, TMatrixTSym<double> VertexCov, double MassConstraint):
  LagrangeMultipliersFitter()
{
  debug = false;
  AnalyticalCovariance_ =false;

  ResPtEstimate_ = ResPtEstimate;
  phiz_ = 0; //not used in this version
  RecoilX_ = ResPtEstimate.Par()(0,0);
  RecoilY_ = ResPtEstimate.Par()(1,0);
  useFullRecoil_ = true;
  MassConstraint_ = MassConstraint;

  Configure(TauThreeProng, ThreeProng, OneProng, OneProngTrack, PVertex, VertexCov);
}

void ThreeProngOneProngFitter::Configure(LorentzVectorParticle TauThreeProng, LorentzVectorParticle ThreeProng, LorentzVectorParticle OneProng, TrackParticle OneProngTrack, TVector3 PVertex, TMatrixTSym<double> VertexCov){
  debug = false;
  AnalyticalCovariance_ =false;
  OneProngTrack_ = OneProngTrack;
  ThreeProng_ = ThreeProng;
  PV_ = PVertex;

  LorentzVectorParticle  TauOneProngGuess;
  LorentzVectorParticle  DiTau;
  if(!useFullRecoil_)
    TauOneProngGuess = TauOneProngStartingPoint( OneProngTrack,TauThreeProng,PVertex, VertexCov, TauThreeProng.Vertex(),TauThreeProng.VertexCov());
  else
    TauOneProngGuess = TauOneProngStartingPointwithFullRecoil(OneProngTrack,TauThreeProng, ResPtEstimate_, PVertex, VertexCov, TauThreeProng.Vertex(),TauThreeProng.VertexCov());

  Logger(Logger::Debug) << "TauThreeProng covariance: " << std::endl;
  if(Logger::Instance()->Level() == Logger::Debug){
    TauThreeProng.getCovMatrix().Print();
  }
  Logger(Logger::Debug) << "TauOneProngGuess covariance: " << std::endl;
  if(Logger::Instance()->Level() == Logger::Debug){
    TauOneProngGuess.getCovMatrix().Print();
  }
  Logger(Logger::Debug) << "Passed: TauOneProngStartingPoint" << std::endl;

  ThetaForConstrTemporaryImplementation_=TauOneProngGuess.LV().Theta();
  particles_.push_back(TauThreeProng);
  particles_.push_back(TauOneProngGuess);

  particles0_.push_back(TauThreeProng);
  particles0_.push_back(TauOneProngGuess);

  Logger(Logger::Debug) << "(TauThreeProng.LV() + TauOneProngGuess.LV()).M(): " << (TauThreeProng.LV() + TauOneProngGuess.LV()).M() << std::endl;

  isconfigured=false;

  int size=particles_.size()*3;
  int sizeTrunc = 3;
  TMatrixT<double>    inpar(size,1);
  TMatrixTSym<double> incov(size);

  // Get primary vertex information
  if(VertexCov.GetNrows()!=LorentzVectorParticle::NVertex)return;

  //set input paramterts:  TauThreeProng - TauOneProng
  inpar(tau3prong_px,0)=TauThreeProng.LV().Px();
  inpar(tau3prong_py,0)=TauThreeProng.LV().Py();
  inpar(tau3prong_pz,0)=TauThreeProng.LV().Pz();

  inpar(tau1prong_px,0)=TauOneProngGuess.LV().Px();
  inpar(tau1prong_py,0)=TauOneProngGuess.LV().Py();
  inpar(tau1prong_pz,0)=TauOneProngGuess.LV().Pz();

  int TauThreeProngoffset=0;
  int TauOneProngoffset=3;

  for(int i=0; i<3;i++){
    for(int j=0; j<3;j++){
      incov(i+TauThreeProngoffset,j+TauThreeProngoffset)=TauThreeProng.Covariance(i+3,j+3);
      incov(i+TauOneProngoffset,j+TauOneProngoffset)=TauOneProngGuess.Covariance(i+3,j+3);
    }
  }

  // store expanded par for computation of final par (assumes fit has neglegible impact on ThreeProng correlations with vertex uncertainties)
  exppar_.ResizeTo(size,1);
  expcov_.ResizeTo(size,size);
  exppar_=ComputeInitalExpPar(inpar);
  expcov_=ErrorMatrixPropagator::PropagateError(&ThreeProngOneProngFitter::ComputeInitalExpPar,inpar,incov);

  // store linearization point
  TMatrixT<double> PAR_0(size,1);
  par_0.ResizeTo(size);
  cov_0.ResizeTo(size,size);
  PAR_0=ComputeExpParToPar(exppar_);
  for(int i=0; i<npar;i++)par_0(i)=PAR_0(i,0);
  cov_0=ErrorMatrixPropagator::PropagateError(&ThreeProngOneProngFitter::ComputeExpParToPar,exppar_,expcov_);

  // set up inital point for fit (cov handled in Fit() function)
  par.ResizeTo(npar);
  par=par_0;


  TMatrixT<double> PARa_0(sizeTrunc,1);
  para_0.ResizeTo(sizeTrunc);
  cova_0.ResizeTo(sizeTrunc,sizeTrunc);
  PARa_0=ComputeExpParToPara(exppar_);
  for(int i=0; i<npartr;i++)para_0(i)=PARa_0(i,0);

  cova_0=ErrorMatrixPropagator::PropagateError(&ThreeProngOneProngFitter::ComputeExpParToPara,exppar_,expcov_);


  para.ResizeTo(npartr);
  para=para_0;
  TMatrixT<double> PARb_0(sizeTrunc,1);
  parb_0.ResizeTo(sizeTrunc);
  covb_0.ResizeTo(sizeTrunc,sizeTrunc);
  PARb_0=ComputeExpParToParb(exppar_);
  for(int i=0; i<npartr;i++)parb_0(i)=PARb_0(i,0);
  y_.ResizeTo(npartr,1); y_ = convertToMatrix(para_0);

  covb_0=ErrorMatrixPropagator::PropagateError(&ThreeProngOneProngFitter::ComputeExpParToParb,exppar_,expcov_);
  parb.ResizeTo(npartr);
  parb=parb_0;

  isconfigured=true;
  Init_Resonance_ = GetMother();


  Logger(Logger::Debug) << "exppar_ covariance: " << std::endl;
  if(Logger::Instance()->Level() == Logger::Debug){
	expcov_.Print();
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

TMatrixT<double> ThreeProngOneProngFitter::ComputeInitalExpPar(TMatrixT<double> &inpar){
  TMatrixT<double> outpar(6,1);
  int offset=0;//LorentzVectorParticle::NVertex;// for TauThreeProng
  outpar(tau3prong_px,0)=inpar(LorentzVectorParticle::vx+offset,0);
  outpar(tau3prong_py,0)=inpar(LorentzVectorParticle::vy+offset,0);
  outpar(tau3prong_pz,0)=inpar(LorentzVectorParticle::vz+offset,0);

  offset+=3;//LorentzVectorParticle::NLorentzandVertexPar; // for TauNu

  outpar(tau1prong_px,0)=inpar(LorentzVectorParticle::vx+offset,0);
  outpar(tau1prong_py,0)=inpar(LorentzVectorParticle::vy+offset,0);
  outpar(tau1prong_pz,0)=inpar(LorentzVectorParticle::vz+offset,0);

  return outpar;

}

bool ThreeProngOneProngFitter::Fit(){
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
    ThreeProngOneProngFitterChiSquareFunctionUpdator updator(this);
    ROOT::Minuit2::MnUserParameters MnPar;
    ROOT::Minuit2::MnUserCovariance MnCov(nPar);
    // Logger(Logger::Info) << "Starting Fit with FittingProc::Minuit" << std::endl;
    // Logger(Logger::Info) << "Setting up parameters a:" << std::endl;
    for(int i=0;i<para_0.GetNrows();i++){
      TString name=ParName(i);
      // if not limited (vhigh <= vlow)
      // MnPar.Add(name.Data(),para_0(i,0),sqrt(fabs(cova_0(i,i))),para_0(i,0)-nsigma*sqrt(fabs(cova_0(i,i))),para_0(i,0)+nsigma*sqrt(fabs(cova_0(i,i))));
      MnPar.Add(name.Data(),para_0(i),sqrt(fabs(cova_0(i,i))));
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
      MnPar.Add(name.Data(),parb_0(i),sqrt(fabs(covb_0(i,i))));
      for(int j=0;j<parb_0.GetNrows();j++){
        MnCov(i+offset,j+offset) = covb_0(i,j);
      }
      // Logger(Logger::Info) << "\t" << name << std::endl;
    }
    // set limits for tau mu pz
    MnPar.SetLimits(5, -5.0*OneProngTrack_.P(), 5.0*OneProngTrack_.P());

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
    // Logger(Logger::Error) << "minimum with limited tau mu pz: " << " 5.0*OneProngTrack_.P(): " << 5.0*OneProngTrack_.P() << min1 << std::endl;
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

TMatrixT<double> ThreeProngOneProngFitter::ComputeExpParToPar(TMatrixT<double> &inpar){
  TMatrixT<double> outpar(npar,1);
  for(int i=0;i<npar;i++){outpar(i,0)=inpar(i,0);}
  return outpar;
}
TMatrixT<double> ThreeProngOneProngFitter::ComputeExpParToPara(TMatrixT<double> &inpar){
  TMatrixT<double> outpar(npartr,1);
  for(int i=0;i<npartr;i++){outpar(i,0)=inpar(i,0);}
  return outpar;
}

TMatrixT<double> ThreeProngOneProngFitter::ComputeExpParToParb(TMatrixT<double> &inpar){
  TMatrixT<double> outpar(npartr,1);
  int offset = 3;
  for(int i=0;i<npartr;i++){outpar(i,0)=inpar(i+offset,0);}
  return outpar;
}

TMatrixT<double> ThreeProngOneProngFitter::ComputeTauOneProngLorentzVectorPar(TMatrixT<double> &inpar){
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

TMatrixT<double> ThreeProngOneProngFitter::ComputeTauThreeProngLorentzVectorPar(TMatrixT<double> &inpar){
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

TMatrixT<double> ThreeProngOneProngFitter::ComputeMotherLorentzVectorPar(TMatrixT<double> &inpar){

  TMatrixT<double> outpar(7,1);
  TMatrixT<double> TauOneProngpar=ComputeTauOneProngLorentzVectorPar(inpar);
  TMatrixT<double> TauThreeProngpar=ComputeTauThreeProngLorentzVectorPar(inpar);
  outpar(LorentzVectorParticle::px,0)=TauOneProngpar(LorentzVectorParticle::px,0)+TauThreeProngpar(LorentzVectorParticle::px,0);
  outpar(LorentzVectorParticle::py,0)=TauOneProngpar(LorentzVectorParticle::py,0)+TauThreeProngpar(LorentzVectorParticle::py,0);
  outpar(LorentzVectorParticle::pz,0)=TauOneProngpar(LorentzVectorParticle::pz,0)+TauThreeProngpar(LorentzVectorParticle::pz,0);

  double Etau1prong2=pow(TauOneProngpar(LorentzVectorParticle::px,0),2.0)+pow(TauOneProngpar(LorentzVectorParticle::py,0),2.0)+pow(TauOneProngpar(LorentzVectorParticle::pz,0),2.0)+pow(TauOneProngpar(LorentzVectorParticle::m,0),2.0);
  double EtauThreeProng2=pow(TauThreeProngpar(LorentzVectorParticle::px,0),2.0)+pow(TauThreeProngpar(LorentzVectorParticle::py,0),2.0)+pow(TauThreeProngpar(LorentzVectorParticle::pz,0),2.0)+pow(TauThreeProngpar(LorentzVectorParticle::m,0),2.0);
  double P2=pow(outpar(LorentzVectorParticle::px,0),2.0)+pow(outpar(LorentzVectorParticle::py,0),2.0)+pow(outpar(LorentzVectorParticle::pz,0),2.0);
  outpar(LorentzVectorParticle::m,0)=sqrt(fabs(pow(sqrt(Etau1prong2)+sqrt(EtauThreeProng2),2.0)-P2));

  return outpar;
}

void ThreeProngOneProngFitter::UpdateExpandedPar(){
  for(int i=0; i<npartr;i++){
	  exppar_(i,0)=para(i);
	  for(int j=0; j<npartr;j++){expcov_(i,j)=cova_0(i,j);}
  }
  int offset = npartr;
  for(int i=0; i<npartr;i++){
	  exppar_(i+offset,0)=parb(i);
	  for(int j=0; j<npartr;j++){expcov_(i+offset,j+offset)=covb_0(i,j);}
  }

}

std::vector<LorentzVectorParticle> ThreeProngOneProngFitter::GetReFitDaughters(){
  std::vector<LorentzVectorParticle> refitParticles;
  UpdateExpandedPar();

  TMatrixT<double> TauThreeProng=ComputeTauThreeProngLorentzVectorPar(exppar_);
  TMatrixTSym<double> TauThreeProngCov=ErrorMatrixPropagator::PropagateError(&ThreeProngOneProngFitter::ComputeTauThreeProngLorentzVectorPar,exppar_,expcov_);

  for(int i=0; i<LorentzVectorParticle::NVertex; i++){
    for(int j=0; j<LorentzVectorParticle::NVertex; j++){
      TauThreeProngCov(i,j)=particles_.at(0).VertexCov()(i,j);
    }
  }

  int charge = particles_.at(0).Charge();
  refitParticles.push_back(LorentzVectorParticle(TauThreeProng,TauThreeProngCov,charge*(PDGInfo::tau_plus),charge,particles_.at(0).BField()));

  TMatrixT<double> TauOneProng=ComputeTauOneProngLorentzVectorPar(exppar_);
  TauOneProng(0,0)= particles_.at(1).Parameter(LorentzVectorParticle::vx);
  TauOneProng(1,0)= particles_.at(1).Parameter(LorentzVectorParticle::vy);
  TauOneProng(2,0)= particles_.at(1).Parameter(LorentzVectorParticle::vz);

  TMatrixTSym<double> TauOneProngcov=ErrorMatrixPropagator::PropagateError(&ThreeProngOneProngFitter::ComputeTauOneProngLorentzVectorPar,exppar_,expcov_);
  for(int i=0; i<LorentzVectorParticle::NVertex; i++){
    for(int j=0; j<LorentzVectorParticle::NVertex; j++){
      TauOneProngcov(i,j)=particles_.at(1).VertexCov()(i,j);
    }
  }

  charge = particles_.at(1).Charge();
  refitParticles.push_back(LorentzVectorParticle(TauOneProng,TauOneProngcov,charge*(PDGInfo::tau_minus),charge,particles_.at(0).BField()));

  return refitParticles;
}

LorentzVectorParticle ThreeProngOneProngFitter::GetMother(){
  UpdateExpandedPar();
  double c(0),b(0);
  for(unsigned int i=0;i<particles_.size();i++){c+=particles_.at(i).Charge();b=particles_.at(i).BField();}
  TMatrixT<double> Mother=ComputeMotherLorentzVectorPar(exppar_);
  TMatrixTSym<double> MotherCov=ErrorMatrixPropagator::PropagateError(&ThreeProngOneProngFitter::ComputeMotherLorentzVectorPar,exppar_,expcov_);

  return LorentzVectorParticle(Mother,MotherCov,PDGInfo::Z0,c,b);
}

TVectorD ThreeProngOneProngFitter::HardValue(TVectorD &va,TVectorD &vb){
  TLorentzVector TauThreeProng,TauOneProng;
  CovertParToObjects(va,vb,TauThreeProng,TauOneProng);

  TLorentzVector z=TauThreeProng+TauOneProng;
  TVectorD d(NConstraints());

	//d(0) = pow(TauThreeProng.E() + TauOneProng.E(), 2.) - (TauThreeProng.Vect() + TauOneProng.Vect()).Mag2()- pow(MassConstraint_,2.);
	//d(1) = TauOneProng.Pz() - CosThetaTauOneProng(TauOneProng)*TauOneProng.P();
	d(0) = z.M() - MassConstraint_;
	d(1) = TauOneProng.Pz() - CosThetaTauOneProng(TauOneProng)*TauOneProng.P();
  return d;
}

TVectorD ThreeProngOneProngFitter::SoftValue(TVectorD &va,TVectorD &vb){
  TLorentzVector TauThreeProng,TauOneProng;
  CovertParToObjects(va,vb,TauThreeProng,TauOneProng);
  TVectorD d(NSoftConstraints());

  if(!useFullRecoil_){
    d(0) = TauThreeProng.Px() + TauOneProng.Px();
    d(1) = TauThreeProng.Py() + TauOneProng.Py();
    d(2) = ( (TauThreeProng.Py() + TauOneProng.Py())/(TauThreeProng.Px() + TauOneProng.Px())) -  tan(phiz_);
    Logger(Logger::Debug) << "SCVec: " << d(0) << ", " << d(1) << ", " << d(2) << std::endl;
    Logger(Logger::Debug) << "RecoilX_: " << RecoilX_ << ", RecoilY_: " << RecoilY_ << std::endl;
  }
  else{
    d(0) = TauThreeProng.Px() + TauOneProng.Px() - RecoilX_;
    d(1) = TauThreeProng.Py() + TauOneProng.Py() - RecoilY_;
    Logger(Logger::Debug) << "SCVec: " << d(0) << ", " << d(1) << std::endl;
    Logger(Logger::Debug) << "RecoilX_: " << RecoilX_ << ", RecoilY_: " << RecoilY_ << std::endl;
  }
  return d;
}

void ThreeProngOneProngFitter::CovertParToObjects(TVectorD &va,TVectorD &vb,TLorentzVector &TauThreeProng,TLorentzVector &TauOneProng){
  // TauThreeProng=particles_.at(0).LV();//TLorentzVector(v(tau3prong_px),v(tau3prong_py),v(tau3prong_pz),sqrt(1.777*1.777+v(tau3prong_px)*v(tau3prong_px)+v(tau3prong_py)*v(tau3prong_py)+v(tau3prong_pz)*v(tau3prong_pz)));
  // TauOneProng=TLorentzVector(v(tau3prong_px),v(tau3prong_py),v(tau3prong_pz),sqrt(1.777*1.777+v(tau3prong_px)*v(tau3prong_px)+v(tau3prong_py)*v(tau3prong_py)+v(tau3prong_pz)*v(tau3prong_pz)));
  TauThreeProng=TLorentzVector(va(tau_px),va(tau_py),va(tau_pz),sqrt(PDGInfo::tau_mass()*PDGInfo::tau_mass()+va(tau_px)*va(tau_px)+va(tau_py)*va(tau_py)+va(tau_pz)*va(tau_pz)));
  TauOneProng=TLorentzVector(vb(tau_px),vb(tau_py),vb(tau_pz),sqrt(PDGInfo::tau_mass()*PDGInfo::tau_mass()+vb(tau_px)*vb(tau_px)+vb(tau_py)*vb(tau_py)+vb(tau_pz)*vb(tau_pz)));
}

LorentzVectorParticle ThreeProngOneProngFitter::TauOneProngStartingPoint(TrackParticle OneProngTrack,LorentzVectorParticle TauThreeProng, TVector3 PV,TMatrixTSym<double>  PVCov, TVector3 SV, TMatrixTSym<double>  SVCov){

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

  TMatrixT<double>    tau1prongdirection;
  tau1prongdirection.ResizeTo(2,1);
  TMatrixTSym<double>    tau1prongdirectionError;
  tau1prongdirectionError.ResizeTo(2,2);

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

  parameters = ConfigureParameters(OneProngTrack, EstimatePhiAngle(TauDir,TauDirError));
  parameterErrors= ConfigureParameterErrors(OneProngTrack, EstimatePhiAngle(TauDir,TauDirError));

  parametersAd = ConfigureInitialAdvancedParameters(OneProngTrack, PV, TauThreeProng);
  parameterErrorsAd= ConfigureInitialAdvancedParameterErrors(OneProngTrack, PVCov, TauThreeProng.getCovMatrix());

  tau1prongdirection = EstimateTauDirectionAdvanced(parametersAd);
  tau1prongdirectionError = ErrorMatrixPropagator::PropagateError(&ThreeProngOneProngFitter::EstimateTauDirectionAdvanced,ConfigureInitialAdvancedParameters(OneProngTrack, PV, TauThreeProng),ConfigureInitialAdvancedParameterErrors(OneProngTrack, PVCov, TauThreeProng.getCovMatrix()));

  kinematicparameters = ConfigureKinematicParameters(tau1prongdirection,TauThreeProng);
  kinematicparametererrors =  ConfigureKinematicParameterErrors(tau1prongdirectionError,TauThreeProng);

  TauKin=EstimateTauKinematic(kinematicparameters);
  TauKinErrorNumerical = ErrorMatrixPropagator::PropagateError(&ThreeProngOneProngFitter::EstimateTauKinematic,ConfigureKinematicParameters(tau1prongdirection,TauThreeProng),ConfigureKinematicParameterErrors(tau1prongdirectionError,TauThreeProng));


  TauKinErrorAnalytical=ComputeAngleCovarianceAnalytically(OneProngTrack,EstimatePhiAngle(TauDir,TauDirError),PV,SV,TauThreeProng);

  if(debug){
    std::cout<<"Tau Kinematic Paramteres 'TauKin' ====>"<<std::endl;
    TauKin.Print();

    std::cout<<"Numericaly  calculated tau kinematic parameter errors 'TauKinErrorNumerical'  ====> "<< std::endl;
    TauKinErrorNumerical.Print();

    std::cout<<"Inpute Kinematic Parameters 'kinematicparameters' ====>"<<std::endl;
    kinematicparameters.Print();

    std::cout<<"Input Kinematic Paramter Errors  'kinematicparametererrors' ====>"<< std::endl;
    kinematicparametererrors.Print();

    std::cout<<"TauOneProng direction paramters 'tau1prongdirection'  ====>  "<< std::endl;
    tau1prongdirection.Print();

    std::cout<<"TauOneProng direction paramter errors  'tau1prongdirectionError'  ====>  "<< std::endl;
    tau1prongdirectionError.Print();

    std::cout<<"Analitically  calculated tau kinematic parameter errors 'TauKinErrorNumerical'  ====>  "<< std::endl;
    TauKinErrorAnalytical.Print();

  }

  TMatrixT<double>    par(LorentzVectorParticle::NLorentzandVertexPar,1);
  TMatrixTSym<double> Cov(LorentzVectorParticle::NLorentzandVertexPar);

  // par(LorentzVectorParticle::vx,0)=0; // fill zero vertex for now TODO: Why?
  // par(LorentzVectorParticle::vy,0)=0;
  // par(LorentzVectorParticle::vz,0)=0;
  par(LorentzVectorParticle::vx,0)=PV_.X();
  par(LorentzVectorParticle::vy,0)=PV_.Y();
  par(LorentzVectorParticle::vz,0)=PV_.Z();
  par(LorentzVectorParticle::px,0)=TauKin(0,0);
  par(LorentzVectorParticle::py,0)=TauKin(1,0);
  par(LorentzVectorParticle::pz,0)=TauKin(2,0);
  par(LorentzVectorParticle::m,0) =PDGInfo::tau_mass();

   for(int i=0; i<LorentzVectorParticle::NVertex; i++){
     for(int j=0; j<LorentzVectorParticle::NVertex; j++){
       if(AnalyticalCovariance_){Cov(i+3,j+3)=TauKinErrorAnalytical(i,j);}
       else{Cov(i+3,j+3)=TauKinErrorNumerical(i,j);}
     }
   }

  return LorentzVectorParticle(par,Cov,PDGInfo::tau_minus,0,0);
}

LorentzVectorParticle ThreeProngOneProngFitter::TauOneProngStartingPointwithFullRecoil(TrackParticle OneProngTrack,LorentzVectorParticle TauThreeProng, PTObject METminusNeutrino, TVector3 PV, TMatrixTSym<double>  PVCov, TVector3 SV, TMatrixTSym<double>  SVCov){

  TMatrixT<double> kinematicparameters;
  kinematicparameters.ResizeTo(12,1);
  TMatrixTSym<double> kinematicparametererrors;
  kinematicparametererrors.ResizeTo(12,12);

  TMatrixT<double> TauKin;
  TauKin.ResizeTo(6,1);

  TMatrixT<double> TauKinErrorNumerical;
  TauKinErrorNumerical.ResizeTo(6,6);
  TMatrixT<double> TauKinErrorAnalytical;
  TauKinErrorAnalytical.ResizeTo(6,6);
  TMatrixTSym<double> TauOneProngErrorAnalytical;
  TauOneProngErrorAnalytical.ResizeTo(3,3);

  kinematicparameters = ConfigureKinematicParametersFullRecoil(OneProngTrack, PV, TauThreeProng, METminusNeutrino);
  kinematicparametererrors =  ConfigureKinematicParameterErrorsFullRecoil(OneProngTrack, PVCov, TauThreeProng, METminusNeutrino);

  TauKin=EstimateTauKinematicFullRecoil(kinematicparameters);
  TauKinErrorNumerical = ErrorMatrixPropagator::PropagateError(&ThreeProngOneProngFitter::EstimateTauKinematicFullRecoil,kinematicparameters,kinematicparametererrors);
  TLorentzVector TauOneProng(TauKin(3,0),TauKin(4,0),TauKin(5,0),PDGInfo::tau_mass());

  TauKinErrorAnalytical = TauKinErrorNumerical;
  TauOneProngErrorAnalytical = EstimateTauKinematicErrorFullRecoil(TauThreeProng, TauOneProng, METminusNeutrino);

  if(Logger::Instance()->Level() == Logger::Debug){
	Logger(Logger::Debug) << "TauKinErrorAnalytical: " << std::endl;
	TauKinErrorAnalytical.Print();
	Logger(Logger::Debug) << "TauOneProngErrorAnalytical: " << std::endl;
	TauOneProngErrorAnalytical.Print();
  }

  TMatrixT<double> par(LorentzVectorParticle::NLorentzandVertexPar,1);
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
      Cov(i+3,j+3)=TauOneProngErrorAnalytical(i,j);
    }
  }
  if(Logger::Instance()->Level() == Logger::Debug){
	Logger(Logger::Debug) << "Cov: " << std::endl;
	Cov.Print();
  }

  return LorentzVectorParticle(par,Cov,PDGInfo::tau_minus,0,0);
}

TMatrixT<double> ThreeProngOneProngFitter::EstimateTauDirectionAdvanced(TMatrixT<double> &inpar){
  TMatrixT<double>    outpar(2,1);

  double dxy   =inpar(0,0);
  double phi0  =inpar(1,0);
  double lam   =inpar(2,0);
  double dz    =inpar(3,0);

  TVector3 PV(inpar(4,0),inpar(5,0),inpar(6,0));
  TVector2 TauOneProngPt(-inpar(7,0), -inpar(8,0));

  double tanphi_tau = TauOneProngPt.Y()/TauOneProngPt.X();
  double a = tanphi_tau;
  double b = PV.Y() - a*PV.X();

  double tmu  = (dxy*(cos(phi0) + a*sin(phi0)) - b) / (a*cos(phi0) - sin(phi0)); //projection onto XY plane

  double xdoc = -dxy*sin(phi0) + tmu*cos(phi0);
  double ydoc = dxy*cos(phi0) + tmu*sin(phi0);
  double zdoc = dz + tmu*tan(lam);

  TVector3 PointGuess(xdoc, ydoc, zdoc);
  TVector3 TauOneProngDir = PointGuess - PV;

  outpar(0,0) = TauOneProngDir.Theta();
  outpar(1,0) = atan2(TauOneProngPt.Y(),TauOneProngPt.X());

  Logger(Logger::Debug) << "PV: " << inpar(4,0) << ", " << inpar(5,0) << ", " << inpar(6,0) << std::endl;
  Logger(Logger::Debug) << "dxy, phi0, lam, dz: " << dxy << ", " << phi0 << ", " << lam << ", " << dz << std::endl;
  Logger(Logger::Debug) << "PointGuess: " << xdoc << ", " << ydoc << ", " << zdoc << std::endl;
  Logger(Logger::Debug) << "DirGuess: " << TauOneProngDir.X() << ", " << TauOneProngDir.Y() << ", " << TauOneProngDir.Z() << std::endl;
  Logger(Logger::Debug) << "a, b, t: " << a << ", " << b << ", " << tmu <<  std::endl;
  Logger(Logger::Debug) << "TauH X, Y: " << inpar(7,0) << ", " << inpar(8,0) << std::endl;
  Logger(Logger::Debug) << "TauOneProngDir.Phi(): " << TauOneProngDir.Phi() << " TauH phi: " << atan2(inpar(8,0), inpar(7,0)) << std::endl;
  Logger(Logger::Debug) << "dPhi: " << TauOneProngDir.Phi() - atan2(inpar(8,0), inpar(7,0)) << "/pi= " << (TauOneProngDir.Phi() - atan2(inpar(8,0), inpar(7,0)))/TMath::Pi() << std::endl;

  return outpar;
}

TMatrixT<double> ThreeProngOneProngFitter::EstimateTauKinematic(TMatrixT<double> &inpar){
  TMatrixT<double>    outpar(3,1);
  TLorentzVector TauThreeProngp4(inpar(2,0),
			 inpar(3,0),
			 inpar(4,0),
			 sqrt(pow(PDGInfo::tau_mass(),2) + pow(inpar(2,0),2) + pow(inpar(3,0),2) + pow(inpar(4,0),2)  ) );

  TVector3 TauOneProngDir(cos(inpar(1,0))*sin(inpar(0,0)), sin(inpar(1,0))*sin(inpar(0,0)), cos(inpar(0,0)));
  TVector3 P_Tauh(inpar(2,0), inpar(3,0), inpar(4,0));
  TLorentzVector P4_Tauh; P4_Tauh.SetXYZM(P_Tauh.X(),P_Tauh.Y(),P_Tauh.Z(),PDGInfo::tau_mass());
  double theta = P_Tauh.Angle(TauOneProngDir);
  double MassDiffsq = pow(MassConstraint_, 2.) - pow(PDGInfo::tau_mass(), 2.);
  double Denominator = pow(P4_Tauh.P()*sin(theta),2.) + pow(PDGInfo::tau_mass(),2.);
  double P_TauOneProng = (MassDiffsq*P4_Tauh.P()*cos(theta) + P4_Tauh.E()*sqrt( pow(MassDiffsq, 2.) - 4*pow(PDGInfo::tau_mass(),2.)*Denominator))/2/(Denominator);

  outpar(0,0) = P_TauOneProng*TauOneProngDir.X();
  outpar(1,0) = P_TauOneProng*TauOneProngDir.Y();
  outpar(2,0) = P_TauOneProng*TauOneProngDir.Z();

  //outpar(0,0) = P_Tauh.Pt()*cos(inpar(1,0));
  //outpar(1,0) = P_Tauh.Pt()*sin(inpar(1,0));
  //outpar(2,0) = P_Tauh.Pt()/tan(inpar(0,0));

  Logger(Logger::Debug) << "TauDir.Phi(): " << TauThreeProngp4.Phi() << " TauOneProngDirNEW2.Phi(): " << inpar(1,0) << std::endl;
  Logger(Logger::Debug) << "TauThreeProng p3: " << TauThreeProngp4.X() << ", " << TauThreeProngp4.Y() << ", " << TauThreeProngp4.Z() << std::endl;
  Logger(Logger::Debug) << "TauOneProng p3: " << outpar(0,0) << ", " << outpar(1,0) << ", " << outpar(2,0) << std::endl;


  return outpar;
}

TMatrixT<double> ThreeProngOneProngFitter::ConfigureParameters(TrackParticle OneProngTrack, std::pair<double, double> phiAngle){
   TMatrixT<double>    outpar(5,1);
   outpar(0,0) = OneProngTrack.Parameter(TrackParticle::dxy);
   outpar(1,0) = OneProngTrack.Parameter(TrackParticle::phi);
   outpar(2,0) = OneProngTrack.Parameter(TrackParticle::lambda);
   outpar(3,0) = OneProngTrack.Parameter(TrackParticle::dz);
   outpar(4,0) = phiAngle.first;

   return outpar;
 }

TMatrixT<double> ThreeProngOneProngFitter::ConfigureParameterErrors(TrackParticle OneProngTrack, std::pair<double, double> phiAngle){
    TMatrixT<double>  Cov;
    Cov.ResizeTo(5,5);

    Cov(0,0) = OneProngTrack.Covariance(TrackParticle::dxy,TrackParticle::dxy);
    Cov(0,1) = OneProngTrack.Covariance(TrackParticle::dxy,TrackParticle::phi);
    Cov(0,2) = OneProngTrack.Covariance(TrackParticle::dxy,TrackParticle::dz);
    Cov(0,3) = OneProngTrack.Covariance(TrackParticle::dxy,TrackParticle::lambda);
    Cov(0,4) = 0;

    Cov(1,0) = OneProngTrack.Covariance(TrackParticle::phi,TrackParticle::dxy);
    Cov(1,1) = OneProngTrack.Covariance(TrackParticle::phi,TrackParticle::phi);
    Cov(1,2) = OneProngTrack.Covariance(TrackParticle::phi,TrackParticle::dz);
    Cov(1,3) = OneProngTrack.Covariance(TrackParticle::phi,TrackParticle::lambda);
    Cov(1,4) = 0;

    Cov(2,0) = OneProngTrack.Covariance(TrackParticle::dz,TrackParticle::dxy);
    Cov(2,1) = OneProngTrack.Covariance(TrackParticle::dz,TrackParticle::phi);
    Cov(2,2) = OneProngTrack.Covariance(TrackParticle::dz,TrackParticle::dz);
    Cov(2,3) = OneProngTrack.Covariance(TrackParticle::dz,TrackParticle::lambda);
    Cov(2,4) = 0;

    Cov(3,0) = OneProngTrack.Covariance(TrackParticle::lambda,TrackParticle::dxy);
    Cov(3,1) = OneProngTrack.Covariance(TrackParticle::lambda,TrackParticle::phi);
    Cov(3,2) = OneProngTrack.Covariance(TrackParticle::lambda,TrackParticle::dz);
    Cov(3,3) = OneProngTrack.Covariance(TrackParticle::lambda,TrackParticle::lambda);
    Cov(3,4) = 0;

    Cov(4,0) = 0;
    Cov(4,1) = 0;
    Cov(4,2) = 0;
    Cov(4,3) = 0;
    Cov(4,4) = phiAngle.second;

    return Cov;
}

TMatrixT<double> ThreeProngOneProngFitter::ConfigureInitialAdvancedParameters(TrackParticle OneProngTrack,  TVector3 PV, LorentzVectorParticle TauThreeProng){
  TMatrixT<double>    outpar;
  outpar.ResizeTo(9,1);
   outpar(0,0)  = OneProngTrack.Parameter(TrackParticle::dxy);
   outpar(1,0)  = OneProngTrack.Parameter(TrackParticle::phi);
   outpar(2,0)  = OneProngTrack.Parameter(TrackParticle::lambda);
   outpar(3,0)  = OneProngTrack.Parameter(TrackParticle::dz);
   outpar(4,0)  = PV.X();
   outpar(5,0)  = PV.Y();
   outpar(6,0)  = PV.Z();
   outpar(7,0) = TauThreeProng.LV().X();
   outpar(8,0) = TauThreeProng.LV().Y();

   return outpar;
 }

TMatrixTSym<double>
ThreeProngOneProngFitter::ConfigureInitialAdvancedParameterErrors(TrackParticle OneProngTrack, TMatrixTSym<double> PVCov, TMatrixTSym<double> TauThreeProngCov){
  TMatrixTSym<double>  Cov;
    Cov.ResizeTo(9,9);

    Cov(0,0) = OneProngTrack.Covariance(TrackParticle::dxy,TrackParticle::dxy);
    Cov(0,1) = OneProngTrack.Covariance(TrackParticle::dxy,TrackParticle::phi);
    Cov(0,2) = OneProngTrack.Covariance(TrackParticle::dxy,TrackParticle::dz);
    Cov(0,3) = OneProngTrack.Covariance(TrackParticle::dxy,TrackParticle::lambda);

    Cov(1,0) = OneProngTrack.Covariance(TrackParticle::phi,TrackParticle::dxy);
    Cov(1,1) = OneProngTrack.Covariance(TrackParticle::phi,TrackParticle::phi);
    Cov(1,2) = OneProngTrack.Covariance(TrackParticle::phi,TrackParticle::dz);
    Cov(1,3) = OneProngTrack.Covariance(TrackParticle::phi,TrackParticle::lambda);

    Cov(2,0) = OneProngTrack.Covariance(TrackParticle::lambda,TrackParticle::dxy);
    Cov(2,1) = OneProngTrack.Covariance(TrackParticle::lambda,TrackParticle::phi);
    Cov(2,2) = OneProngTrack.Covariance(TrackParticle::lambda,TrackParticle::dz);
    Cov(2,3) = OneProngTrack.Covariance(TrackParticle::lambda,TrackParticle::lambda);

    Cov(3,0) = OneProngTrack.Covariance(TrackParticle::dz,TrackParticle::dxy);
    Cov(3,1) = OneProngTrack.Covariance(TrackParticle::dz,TrackParticle::phi);
    Cov(3,2) = OneProngTrack.Covariance(TrackParticle::dz,TrackParticle::dz);
    Cov(3,3) = OneProngTrack.Covariance(TrackParticle::dz,TrackParticle::lambda);

    for(int i=0; i<3; i++){
      for(int j=0; j<3; j++){
        Cov(i+4,j+4)=PVCov(i,j);
      }
    }
    for(int i=0; i<2; i++){
      for(int j=0; j<2; j++){
        Cov(i+7,j+7)=TauThreeProngCov(i+4,j+4);
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

TMatrixT<double> ThreeProngOneProngFitter::ConfigureKinematicParameters(TMatrixT<double>  TauOneProngDir, LorentzVectorParticle TauThreeProng){
  TMatrixT<double>    outpar;
  outpar.ResizeTo(5,1);
  outpar(0,0)  = TauOneProngDir(0,0);
  outpar(1,0)  = TauOneProngDir(1,0);
  outpar(2,0)  = TauThreeProng.LV().Px();
  outpar(3,0)  = TauThreeProng.LV().Py();
  outpar(4,0)  = TauThreeProng.LV().Pz();


  if(debug){
    std::cout<<"EstimateTauKinematic  TauThreeProng  "<<TauThreeProng.LV().Px()<<"  "<<TauThreeProng.LV().Py()<<"  "<<TauThreeProng.LV().Pz()<<std::endl;
    std::cout<<"ConfigureKinematicParameters:  outpar(0,0)  "<<outpar(0,0) <<std::endl;
    std::cout<<"ConfigureKinematicParameters:  outpar(1,0)  "<<outpar(1,0) <<std::endl;
    std::cout<<"ConfigureKinematicParameters:  outpar(2,0)  "<<outpar(2,0) <<std::endl;
    std::cout<<"ConfigureKinematicParameters:  outpar(3,0)  "<<outpar(3,0) <<std::endl;
    std::cout<<"ConfigureKinematicParameters:  outpar(4,0)  "<<outpar(4,0) <<std::endl;
  }
  return outpar;
 }

TMatrixTSym<double> ThreeProngOneProngFitter::ConfigureKinematicParameterErrors(TMatrixTSym<double>  TauOneProngDirError, LorentzVectorParticle TauThreeProng){
  TMatrixTSym<double> outpar;
  outpar.ResizeTo(5,5);

  for(int i=0; i<5; i++){
    for(int j=0; j<5; j++){
      if(i<2 && j < 2){
        outpar(i,j)=TauOneProngDirError(i,j);
      }else{
        outpar(i,j)=TauThreeProng.Covariance(i+1,j+1);
      }
    }
  }
  return outpar;
}

std::pair<double, double> ThreeProngOneProngFitter::EstimatePhiAngle( TVector3 dir, TVector3 dirE){
  std::pair<double, double> outpar;
  double phi = atan2(dir.Y(),dir.X());
  double deltaphi = sqrt(dirE.Y()*dirE.Y()/dir.X()/dir.X()+
                    dirE.X()*dirE.X()*dir.Y()*dir.Y()/dir.X()/dir.X()/dir.X()/dir.X() )*cos(phi)*cos(phi);
  outpar = std::make_pair(phi, deltaphi);
  return outpar;
}

//////////////////////////////////////////////////////////////////////
// Analytical calculation of TauOneProng Covariance; Has to be rewritten....
TMatrixT<double> ThreeProngOneProngFitter::ComputeAngleCovarianceAnalytically(TrackParticle OneProngTrack, std::pair<double, double> phiAngle,  TVector3 PV, TVector3 SV, LorentzVectorParticle  TauThreeProng){
  TMatrixT<double>    outpar;
  outpar.ResizeTo(3,3);
  double dxy   =OneProngTrack.Parameter(TrackParticle::dxy);
  //double kappa =OneProngTrack.Parameter(TrackParticle::kappa);
  double phi0  =OneProngTrack.Parameter(TrackParticle::phi);
  double lam   =OneProngTrack.Parameter(TrackParticle::lambda);
  double dz    =OneProngTrack.Parameter(TrackParticle::dz);
  //double c     =OneProngTrack.Charge();
  TVector3 ThreeProngSV = -SV + PV;

  double phiAnot  = atan2(ThreeProngSV.Y(), ThreeProngSV.X());
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
  if(fabs(phiAnot - TauThreeProng.LV().Phi()) < 0.05) phiAnot = phiAnot - TMath::Pi();
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

  double cosTauTau2 = cos(TauThreeProng.LV().Theta() + acos(cosTheta));
  //double sinTauTau2 = sin(TauThreeProng.LV().Theta() + acos(cosTheta));

  TMatrixT<double>  HelixCov;
  HelixCov.ResizeTo(5,5);

  HelixCov(0,0) = OneProngTrack.Covariance(TrackParticle::dxy,TrackParticle::dxy);
  HelixCov(0,1) = OneProngTrack.Covariance(TrackParticle::dxy,TrackParticle::phi);
  HelixCov(0,2) = OneProngTrack.Covariance(TrackParticle::dxy,TrackParticle::dz);
  HelixCov(0,3) = OneProngTrack.Covariance(TrackParticle::dxy,TrackParticle::lambda);
  HelixCov(0,4) = 0;

  HelixCov(1,0) = OneProngTrack.Covariance(TrackParticle::phi,TrackParticle::dxy);
  HelixCov(1,1) = OneProngTrack.Covariance(TrackParticle::phi,TrackParticle::phi);
  HelixCov(1,2) = OneProngTrack.Covariance(TrackParticle::phi,TrackParticle::dz);
  HelixCov(1,3) = OneProngTrack.Covariance(TrackParticle::phi,TrackParticle::lambda);
  HelixCov(1,4) = 0;

  HelixCov(2,0) = OneProngTrack.Covariance(TrackParticle::dz,TrackParticle::dxy);
  HelixCov(2,1) = OneProngTrack.Covariance(TrackParticle::dz,TrackParticle::phi);
  HelixCov(2,2) = OneProngTrack.Covariance(TrackParticle::dz,TrackParticle::dz);
  HelixCov(2,3) = OneProngTrack.Covariance(TrackParticle::dz,TrackParticle::lambda);
  HelixCov(2,4) = 0;

  HelixCov(3,0) = OneProngTrack.Covariance(TrackParticle::lambda,TrackParticle::dxy);
  HelixCov(3,1) = OneProngTrack.Covariance(TrackParticle::lambda,TrackParticle::phi);
  HelixCov(3,2) = OneProngTrack.Covariance(TrackParticle::lambda,TrackParticle::dz);
  HelixCov(3,3) = OneProngTrack.Covariance(TrackParticle::lambda,TrackParticle::lambda);
  HelixCov(3,4) = 0;

  HelixCov(4,0) = 0;
  HelixCov(4,1) = 0;
  HelixCov(4,2) = 0;
  HelixCov(4,3) = 0;
  HelixCov(4,4) = phiAngle.second;//sqrt(TauDirError.Y()*TauDirError.Y()/ThreeProngSV.X()/ThreeProngSV.X()   + TauDirError.X()*TauDirError.X()*ThreeProngSV.Y()*ThreeProngSV.Y()/ThreeProngSV.X()/ThreeProngSV.X()/ThreeProngSV.X()/ThreeProngSV.X() ) *cos(phiAnot)*cos(phiAnot);
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

  //double ZMassR = 125.2;
  double TauThreeProngdeltaP = sqrt(   (pow(TauThreeProng.LV().Px(),2)*fabs(TauThreeProng.Covariance(3,3))   +  pow(TauThreeProng.LV().Py(),2)*fabs(TauThreeProng.Covariance(4,4))  +  pow(TauThreeProng.LV().Pz(),2)*fabs(TauThreeProng.Covariance(5,5))    )/TauThreeProng.LV().P()/TauThreeProng.LV().P()  )  ;
  double TauOneProngdeltaP_2 = TauThreeProngdeltaP*pow(MassConstraint_/2/TauThreeProng.LV().P(),2)/(1-cosTauTau2);

  double TauOneProngP_2 =  MassConstraint_*MassConstraint_/2/(1-cosTauTau2)/TauThreeProng.LV().P();


  TMatrixT<double> CovPPhiThetaFrame;
  CovPPhiThetaFrame.ResizeTo(3,3);
  //---------- dp dphi dtheta
  CovPPhiThetaFrame(0,0) = TauOneProngdeltaP_2;
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
  DirevativesPPhiThetaToPxPyPz(0,1) = -TauOneProngP_2*sin(phiAnot)*sinTheta;
  DirevativesPPhiThetaToPxPyPz(0,2) =  TauOneProngP_2*cos(phiAnot)*cosTheta;

  DirevativesPPhiThetaToPxPyPz(1,0) = sin(phiAnot)*sinTheta;
  DirevativesPPhiThetaToPxPyPz(1,1) = TauOneProngP_2*cos(phiAnot)*sinTheta;
  DirevativesPPhiThetaToPxPyPz(1,2) = TauOneProngP_2*sin(phiAnot)*cosTheta;

  DirevativesPPhiThetaToPxPyPz(2,0) = cosTheta;
  DirevativesPPhiThetaToPxPyPz(2,1) = 0;
  DirevativesPPhiThetaToPxPyPz(2,2) =-TauOneProngP_2*sinTheta;


  TMatrixT<double> DirevativesPPhiThetaToPxPyPzT=DirevativesPPhiThetaToPxPyPz; DirevativesPPhiThetaToPxPyPzT.T();
  TMatrixT<double> CovPxPyPzFrame=DirevativesPPhiThetaToPxPyPz*CovPPhiThetaFrame*DirevativesPPhiThetaToPxPyPzT;


  outpar = CovPxPyPzFrame;
  return outpar;
 }

LorentzVectorParticle ThreeProngOneProngFitter::GetTauOneProngEstimate(){
  return particles_.at(1);
}

TMatrixT<double> ThreeProngOneProngFitter::ConfigureKinematicParametersFullRecoil(TrackParticle OneProngTrack, TVector3 PV, LorentzVectorParticle TauThreeProng, PTObject METMinusNeutrino){
  TMatrixT<double> outpar;
  outpar.ResizeTo(12,1);

  outpar(0,0)  = OneProngTrack.Parameter(TrackParticle::dxy);
  outpar(1,0)  = OneProngTrack.Parameter(TrackParticle::phi);
  outpar(2,0)  = OneProngTrack.Parameter(TrackParticle::lambda);
  outpar(3,0)  = OneProngTrack.Parameter(TrackParticle::dz);
  outpar(4,0)  = PV.X();
  outpar(5,0)  = PV.Y();
  outpar(6,0)  = PV.Z();
  outpar(7,0)  = TauThreeProng.LV().X();
  outpar(8,0)  = TauThreeProng.LV().Y();
  outpar(9,0)  = TauThreeProng.LV().Z();
  outpar(10,0) = METMinusNeutrino.Par()(0,0);
  outpar(11,0) = METMinusNeutrino.Par()(1,0);

  return outpar;
}

TMatrixTSym<double> ThreeProngOneProngFitter::ConfigureKinematicParameterErrorsFullRecoil(TrackParticle OneProngTrack, TMatrixTSym<double> PVCov, LorentzVectorParticle TauThreeProng, PTObject METMinusNeutrino){
  TMatrixTSym<double>  Cov;
    Cov.ResizeTo(12,12);

    Cov(0,0) = OneProngTrack.Covariance(TrackParticle::dxy,TrackParticle::dxy);
    Cov(0,1) = OneProngTrack.Covariance(TrackParticle::dxy,TrackParticle::phi);
    Cov(0,2) = OneProngTrack.Covariance(TrackParticle::dxy,TrackParticle::dz);
    Cov(0,3) = OneProngTrack.Covariance(TrackParticle::dxy,TrackParticle::lambda);

    Cov(1,0) = OneProngTrack.Covariance(TrackParticle::phi,TrackParticle::dxy);
    Cov(1,1) = OneProngTrack.Covariance(TrackParticle::phi,TrackParticle::phi);
    Cov(1,2) = OneProngTrack.Covariance(TrackParticle::phi,TrackParticle::dz);
    Cov(1,3) = OneProngTrack.Covariance(TrackParticle::phi,TrackParticle::lambda);

    Cov(2,0) = OneProngTrack.Covariance(TrackParticle::lambda,TrackParticle::dxy);
    Cov(2,1) = OneProngTrack.Covariance(TrackParticle::lambda,TrackParticle::phi);
    Cov(2,2) = OneProngTrack.Covariance(TrackParticle::lambda,TrackParticle::dz);
    Cov(2,3) = OneProngTrack.Covariance(TrackParticle::lambda,TrackParticle::lambda);

    Cov(3,0) = OneProngTrack.Covariance(TrackParticle::dz,TrackParticle::dxy);
    Cov(3,1) = OneProngTrack.Covariance(TrackParticle::dz,TrackParticle::phi);
    Cov(3,2) = OneProngTrack.Covariance(TrackParticle::dz,TrackParticle::dz);
    Cov(3,3) = OneProngTrack.Covariance(TrackParticle::dz,TrackParticle::lambda);

    for(int i=0; i<3; i++){
      for(int j=0; j<3; j++){
        Cov(i+4,j+4)=PVCov(i,j);
        Cov(i+7,j+7)=TauThreeProng.Covariance(i+4,j+4);
      }
    }
    for(int i=0; i<2; i++){
      for(int j=0; j<2; j++){
        Cov(i+10,j+10)=METMinusNeutrino.Cov()(i,j);
      }
    }
    return Cov;
}

TMatrixT<double> ThreeProngOneProngFitter::EstimateTauKinematicFullRecoil(TMatrixT<double> &inpar){
  TMatrixT<double> outpar;
  outpar.ResizeTo(6,1);

  double dxy   =inpar(0,0);
  double phi0  =inpar(1,0);
  double lam   =inpar(2,0);
  double dz    =inpar(3,0);

  TVector3 PV(inpar(4,0),inpar(5,0),inpar(6,0));
  TVector3 P_Tauh(inpar(7,0), inpar(8,0), inpar(9,0));
  TLorentzVector P4_Tauh; P4_Tauh.SetXYZM(P_Tauh.X(),P_Tauh.Y(),P_Tauh.Z(),PDGInfo::tau_mass());
  TVector2 TauOneProngPt(inpar(10,0) - P_Tauh.X(), inpar(11,0) - P_Tauh.Y());
  TVector3 P_TauOneProng(inpar(10,0) - P_Tauh.X(), inpar(11,0) - P_Tauh.Y(), 0);

  double tanphi_tau = TauOneProngPt.Y()/TauOneProngPt.X();
  double a = tanphi_tau;
  double b = PV.Y() - a*PV.X();

  double tmu  = (dxy*(cos(phi0) + a*sin(phi0)) - b) / (a*cos(phi0) - sin(phi0)); //projection onto XY plane

  double xdoc = -dxy*sin(phi0) + tmu*cos(phi0);
  double ydoc = dxy*cos(phi0) + tmu*sin(phi0);
  double zdoc = dz + tmu*tan(lam);

  TVector3 PointGuess(xdoc, ydoc, zdoc);
  TVector3 TauOneProngDir = PointGuess - PV;

  // double phitau = TauOneProngDir.Phi();
  double thetatau = TauOneProngDir.Theta();

  outpar(0,0) = P4_Tauh.X();
  outpar(1,0) = P4_Tauh.Y();
  outpar(2,0) = P4_Tauh.Z();

  if(useCollinearityTauOneProng_){
    outpar(3,0) = TauOneProngPt.X();
    outpar(4,0) = TauOneProngPt.Y();
    P_TauOneProng.SetZ(TauOneProngPt.Mod()*tan(lam));
    outpar(5,0) = P_TauOneProng.Z();
  }
  else{
    outpar(3,0) = TauOneProngPt.X();
    outpar(4,0) = TauOneProngPt.Y();
    P_TauOneProng.SetZ(TauOneProngPt.Mod()/tan(thetatau));
    outpar(5,0) = P_TauOneProng.Z();
  }

  return outpar;
}

TMatrixDSym ThreeProngOneProngFitter::EstimateTauKinematicErrorFullRecoil(LorentzVectorParticle TauThreeProng, TLorentzVector TauOneProng, PTObject ResPtEstimate){
  double TauOneProngPt = sqrt(pow(ResPtEstimate.X() - TauThreeProng.LV().X(), 2.) + pow(ResPtEstimate.Y() - TauThreeProng.LV().Y(), 2.));
  TMatrixDSym TauThreeProng_ZPt_Cov(4);
  for(int i_row = 0; i_row<ResPtEstimate.Cov().GetNrows(); i_row++){
    for(int i_col = 0; i_col<ResPtEstimate.Cov().GetNcols(); i_col++){
      TauThreeProng_ZPt_Cov(i_row,i_col) = TauThreeProng.Covariance(i_row + LorentzVectorParticle::px, i_col + LorentzVectorParticle::px);
      TauThreeProng_ZPt_Cov(i_row+2,i_col+2) = ResPtEstimate.Cov()(i_row, i_col);
    }
  }
  TMatrixD Jacobi(3,4);
  Jacobi(0,0) = -1.;
  Jacobi(0,2) = +1.;
  Jacobi(1,1) = -1.;
  Jacobi(1,3) = +1.;
  Jacobi(2,0) = -TauThreeProng.LV().X()/TauOneProngPt*tan(TauOneProng.Theta());
  Jacobi(2,1) = +ResPtEstimate.X()/TauOneProngPt*tan(TauOneProng.Theta());
  Jacobi(2,2) = -TauThreeProng.LV().Y()/TauOneProngPt*tan(TauOneProng.Theta());
  Jacobi(2,3) = +ResPtEstimate.Y()/TauOneProngPt*tan(TauOneProng.Theta());

  TMatrixDSym TauOneProngCov(TauThreeProng_ZPt_Cov);
  TauOneProngCov.Similarity(Jacobi);

  if(Logger::Instance()->Level() == Logger::Debug){
    Logger(Logger::Debug) << "TauThreeProng_ZPt_Cov: " << std::endl;
    TauThreeProng_ZPt_Cov.Print();
    Logger(Logger::Debug) << "Jacobi: " << std::endl;
    Jacobi.Print();
    Logger(Logger::Debug) << "TauOneProngCov: " << std::endl;
    TauOneProngCov.Print();
  }

  return TauOneProngCov;
}

double ThreeProngOneProngFitter::CosThetaTauOneProng(TLorentzVector TauOneProng){
  double dxy   = OneProngTrack_.Parameter(TrackParticle::dxy);
  double phi0  = OneProngTrack_.Parameter(TrackParticle::phi);
  double lam   = OneProngTrack_.Parameter(TrackParticle::lambda);
  double dz    = OneProngTrack_.Parameter(TrackParticle::dz);
  double tanphi_tau = TauOneProng.Y()/TauOneProng.X();
  double a = tanphi_tau;
  double b = PV_.Y() - a*PV_.X();

  double tmu  = (dxy*(cos(phi0) + a*sin(phi0)) - b) / (a*cos(phi0) - sin(phi0)); //projection onto XY plane

  double xdoc = -dxy*sin(phi0) + tmu*cos(phi0);
  double ydoc = dxy*cos(phi0) + tmu*sin(phi0);
  double zdoc = dz + tmu*tan(lam);

  TVector3 PointGuess(xdoc, ydoc, zdoc);
  TVector3 TauOneProngDir = PointGuess - PV_;

  if(useCollinearityTauOneProng_){
    return sin(lam);
  }
  else{
    return TauOneProngDir.CosTheta();
  }
}

bool ThreeProngOneProngFitter::isConverged(){
  if(!useFullRecoil_){
    // if(pardelta<MaxParDelta_ /*&& harddelta_vec.Norm1() < MaxHCDelta_ && chi2prev-chi2 < MaxChi2Delta_  && chi2prev>chi2*/){
    if(/*pardelta<MaxParDelta_ &&*/ harddelta_vec.Norm1() < MaxHCDelta_ && fabs(chi2prev-chi2) < MaxChi2Delta_ /*&& softdelta_vec.Norm1() < MaxSCDelta_ && chi2prev>chi2*/
      && (pow(para(0), 2.0) + pow(para(1), 2.0)) > ThreeProng_.LV().Perp2()
      && (pow(parb(0), 2.0) + pow(parb(1), 2.0)) > pow(OneProngTrack_.Pt(), 2.0)
    ){
      // Logger(Logger::Verbose) << "converged " << delta << " chi2 " <<  chi2 << " chi2prev " << chi2prev <<"  Maxdelta  " <<MaxDelta_ <<std::endl;
      return true;
    }
  }
  else{
    //if(pardelta<MaxParDelta_ && harddelta_vec.Norm1() < MaxHCDelta_ && softdelta_vec.Norm1() < MaxSCDelta_){
    // if(/*pardelta<MaxParDelta_ &&*/ harddelta_vec.Norm1() < MaxHCDelta_ && fabs(chi2prev-chi2) < MaxChi2Delta_ /*&& chi2prev>chi2*/){
    if(/*pardelta<MaxParDelta_ && harddelta_vec.Norm1() < MaxHCDelta_ && */fabs(chi2prev-chi2) < MaxChi2Delta_ /*&& chi2prev>chi2*/
      && (pow(para(0), 2.0) + pow(para(1), 2.0)) > ThreeProng_.LV().Perp2()
      && (pow(parb(0), 2.0) + pow(parb(1), 2.0)) > pow(OneProngTrack_.Pt(), 2.0)
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

TString ThreeProngOneProngFitter::ParName(int par){
  switch (par){
    case tau3prong_px: return "tau3prong_px";
    case tau3prong_py: return "tau3prong_py";
    case tau3prong_pz: return "tau3prong_pz";
    case tau1prong_px: return "tau1prong_px";
    case tau1prong_py: return "tau1prong_py";
    case tau1prong_pz: return "tau1prong_pz";
    case lambda_1: return "lambda_1";
    case lambda_2: return "lambda_2";
    // add more/remove lambdas if number of hard constraints exceeds number of lambdas or vice versa
    default:
      Logger(Logger::Warning) << "Paramater Index " << par << " out of bounds !!";
      return "";
  }
}