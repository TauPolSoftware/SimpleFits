#include "TauPolSoftware/SimpleFits/interface/ThreeProngThreeProngFitter.h"
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

double ThreeProngThreeProngFitter::MassConstraint_ = 125.2;
bool ThreeProngThreeProngFitter::useCollinearityTauMu_ = false;

ThreeProngThreeProngFitter::ThreeProngThreeProngFitter(std::vector< LorentzVectorParticle > TauThreeProngs, std::vector< LorentzVectorParticle > ThreeProngs, TVector3 PVertex, TMatrixTSym<double> VertexCov){
  ThreeProngThreeProngFitter(TauThreeProngs, ThreeProngs, PTObject(), PVertex, VertexCov, ThreeProngThreeProngFitter::MassConstraint_);
}

ThreeProngThreeProngFitter::ThreeProngThreeProngFitter(std::vector< LorentzVectorParticle > TauThreeProngs, std::vector< LorentzVectorParticle > ThreeProngs, PTObject ResPtEstimate, TVector3 PVertex, TMatrixTSym<double> VertexCov){
  ThreeProngThreeProngFitter(TauThreeProngs, ThreeProngs, ResPtEstimate, PVertex, VertexCov, ThreeProngThreeProngFitter::MassConstraint_);
}

ThreeProngThreeProngFitter::ThreeProngThreeProngFitter(std::vector< LorentzVectorParticle > TauThreeProngs, std::vector< LorentzVectorParticle > ThreeProngs, PTObject ResPtEstimate, TVector3 PVertex, TMatrixTSym<double> VertexCov, double MassConstraint):
  LagrangeMultipliersFitter(),
  particles_ (ThreeProngs),
  particles0_(ThreeProngs),
  PV_(PVertex),
  PVCov_(VertexCov),
  RecoilX_(ResPtEstimate.Par()(0,0)),
  RecoilY_(ResPtEstimate.Par()(1,0)),
  ResPtEstimate_(ResPtEstimate),
  debug_(false),
  AnalyticalCovariance_(false)
{
  MassConstraint_ = MassConstraint;
  if(ResPtEstimate.isValid()) useFullRecoil_ = true;
  Configure(TauThreeProngs, ThreeProngs, VertexCov);
}

void ThreeProngThreeProngFitter::Configure(std::vector< LorentzVectorParticle > TauThreeProngs, std::vector< LorentzVectorParticle > ThreeProngs, TMatrixTSym<double> VertexCov){

  Logger(Logger::Debug) << "Tau 1 covariance: " << std::endl;
  if(Logger::Instance()->Level() == Logger::Debug){
    particles0_.at(0).getCovMatrix().Print();
  }
  Logger(Logger::Debug) << "Tau 2 covariance: " << std::endl;
  if(Logger::Instance()->Level() == Logger::Debug){
    particles0_.at(1).getCovMatrix().Print();
  }

  Logger(Logger::Debug) << "(ThreeProngs.at(0).LV() + ThreeProngs.at(0).LV()).M(): " << (particles0_.at(0).LV() + particles0_.at(0).LV()).M() << std::endl;

  isconfigured=false;

  int size = particles_.size()*3;
  int sizeTrunc = 3;
  TMatrixT<double>    inpar(size,1);
  TMatrixTSym<double> incov(size);

  // Get primary vertex information
  if(VertexCov.GetNrows() != LorentzVectorParticle::NVertex)
    return;

  //set input paramterts:  TauA1 - TauMu
  inpar(tau1_px,0)=particles_.at(0).LV().Px();
  inpar(tau1_py,0)=particles_.at(0).LV().Py();
  inpar(tau1_pz,0)=particles_.at(0).LV().Pz();

  inpar(tau2_px,0)=particles_.at(1).LV().Px();
  inpar(tau2_py,0)=particles_.at(1).LV().Py();
  inpar(tau2_pz,0)=particles_.at(1).LV().Pz();

  int Tau1Offset=0;
  int Tau2Offset=3;

  for(int i=0; i<3;i++){
    for(int j=0; j<3;j++){
      incov(i+Tau1Offset,j+Tau1Offset)=particles_.at(0).Covariance(i+3,j+3);
      incov(i+Tau2Offset,j+Tau2Offset)=particles_.at(1).Covariance(i+3,j+3);
    }
  }

  // store expanded par for computation of final par (assumes fit has neglegible impact on a1 correlations with vertex uncertainties)
  exppar_.ResizeTo(size,1);
  expcov_.ResizeTo(size,size);
  exppar_=ComputeInitalExpPar(inpar);
  expcov_=ErrorMatrixPropagator::PropagateError(&ThreeProngThreeProngFitter::ComputeInitalExpPar,inpar,incov);

  // store linearization point
  TMatrixT<double> PAR_0(size,1);
  par_0.ResizeTo(size);
  cov_0.ResizeTo(size,size);
  PAR_0=ComputeExpParToPar(exppar_);
  for(int i=0; i<npar;i++)par_0(i)=PAR_0(i,0);
  cov_0=ErrorMatrixPropagator::PropagateError(&ThreeProngThreeProngFitter::ComputeExpParToPar,exppar_,expcov_);

  // set up inital point for fit (cov handled in Fit() function)
  par.ResizeTo(npar);
  par=par_0;

  TMatrixT<double> PARa_0(sizeTrunc,1);
  para_0.ResizeTo(sizeTrunc);
  cova_0.ResizeTo(sizeTrunc,sizeTrunc);
  PARa_0=ComputeExpParToPara(exppar_);
  for( int i = 0; i<npartr; i++ )
    para_0(i) = PARa_0(i,0);

  cova_0=ErrorMatrixPropagator::PropagateError(&ThreeProngThreeProngFitter::ComputeExpParToPara,exppar_,expcov_);

  para.ResizeTo(npartr);
  para=para_0;
  TMatrixT<double> PARb_0(sizeTrunc,1);
  parb_0.ResizeTo(sizeTrunc);
  covb_0.ResizeTo(sizeTrunc,sizeTrunc);
  PARb_0=ComputeExpParToParb(exppar_);
  for(int i=0; i<npartr;i++)parb_0(i)=PARb_0(i,0);
  y_.ResizeTo(npartr,1); y_ = convertToMatrix(para_0);

  covb_0=ErrorMatrixPropagator::PropagateError(&ThreeProngThreeProngFitter::ComputeExpParToParb,exppar_,expcov_);
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

TMatrixT<double> ThreeProngThreeProngFitter::ComputeInitalExpPar(TMatrixT<double> &inpar){
  TMatrixT<double> outpar(6,1);
  int offset=0;//LorentzVectorParticle::NVertex;// for TauA1
  outpar(tau1_px,0)=inpar(LorentzVectorParticle::vx+offset,0);
  outpar(tau1_py,0)=inpar(LorentzVectorParticle::vy+offset,0);
  outpar(tau1_pz,0)=inpar(LorentzVectorParticle::vz+offset,0);

  offset+=3;//LorentzVectorParticle::NLorentzandVertexPar; // for TauNu

  outpar(tau2_px,0)=inpar(LorentzVectorParticle::vx+offset,0);
  outpar(tau2_py,0)=inpar(LorentzVectorParticle::vy+offset,0);
  outpar(tau2_pz,0)=inpar(LorentzVectorParticle::vz+offset,0);

  return outpar;

}

bool ThreeProngThreeProngFitter::Fit(){
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
    ThreeProngThreeProngFitterChiSquareFunctionUpdator updator(this);
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
    // MnPar.SetLimits(5, -5.0*MuTrack_.P(), 5.0*MuTrack_.P());

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

TMatrixT<double> ThreeProngThreeProngFitter::ComputeExpParToPar(TMatrixT<double> &inpar){
  TMatrixT<double> outpar(npar,1);
  for(int i=0;i<npar;i++){outpar(i,0)=inpar(i,0);}
  return outpar;
}
TMatrixT<double> ThreeProngThreeProngFitter::ComputeExpParToPara(TMatrixT<double> &inpar){
  TMatrixT<double> outpar(npartr,1);
  for(int i=0;i<npartr;i++){outpar(i,0)=inpar(i,0);}  
  return outpar;
}

TMatrixT<double> ThreeProngThreeProngFitter::ComputeExpParToParb(TMatrixT<double> &inpar){
  TMatrixT<double> outpar(npartr,1);
  int offset = 3; 
  for(int i=0;i<npartr;i++){outpar(i,0)=inpar(i+offset,0);}
  return outpar;
}

TMatrixT<double> ThreeProngThreeProngFitter::ComputeTau2LorentzVectorPar(TMatrixT<double> &inpar){
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

TMatrixT<double> ThreeProngThreeProngFitter::ComputeTau1LorentzVectorPar(TMatrixT<double> &inpar){
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

TMatrixT<double> ThreeProngThreeProngFitter::ComputeMotherLorentzVectorPar(TMatrixT<double> &inpar){

  TMatrixT<double> outpar(7,1);
  TMatrixT<double> Taumupar=ComputeTau2LorentzVectorPar(inpar);
  TMatrixT<double> Taua1par=ComputeTau1LorentzVectorPar(inpar);
  outpar(LorentzVectorParticle::px,0)=Taumupar(LorentzVectorParticle::px,0)+Taua1par(LorentzVectorParticle::px,0);
  outpar(LorentzVectorParticle::py,0)=Taumupar(LorentzVectorParticle::py,0)+Taua1par(LorentzVectorParticle::py,0);
  outpar(LorentzVectorParticle::pz,0)=Taumupar(LorentzVectorParticle::pz,0)+Taua1par(LorentzVectorParticle::pz,0);

  double Etaumu2=pow(Taumupar(LorentzVectorParticle::px,0),2.0)+pow(Taumupar(LorentzVectorParticle::py,0),2.0)+pow(Taumupar(LorentzVectorParticle::pz,0),2.0)+pow(Taumupar(LorentzVectorParticle::m,0),2.0);
  double Etaua12=pow(Taua1par(LorentzVectorParticle::px,0),2.0)+pow(Taua1par(LorentzVectorParticle::py,0),2.0)+pow(Taua1par(LorentzVectorParticle::pz,0),2.0)+pow(Taua1par(LorentzVectorParticle::m,0),2.0);
  double P2=pow(outpar(LorentzVectorParticle::px,0),2.0)+pow(outpar(LorentzVectorParticle::py,0),2.0)+pow(outpar(LorentzVectorParticle::pz,0),2.0);
  outpar(LorentzVectorParticle::m,0)=sqrt(fabs(pow(sqrt(Etaumu2)+sqrt(Etaua12),2.0)-P2));

  return outpar;
}

void ThreeProngThreeProngFitter::UpdateExpandedPar(){
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

std::vector<LorentzVectorParticle> ThreeProngThreeProngFitter::GetReFitDaughters(){
  std::vector<LorentzVectorParticle> refitParticles;
  UpdateExpandedPar();

  double c(0),b(0);
  for(unsigned int i=0;i<particles_.size();i++){c+=particles_.at(i).Charge();b=particles_.at(i).BField();}
  TMatrixT<double> a1=ComputeTau1LorentzVectorPar(exppar_);
  TMatrixTSym<double> a1cov=ErrorMatrixPropagator::PropagateError(&ThreeProngThreeProngFitter::ComputeTau1LorentzVectorPar,exppar_,expcov_);

   for(int i=0; i<LorentzVectorParticle::NVertex; i++){
     for(int j=0; j<LorentzVectorParticle::NVertex; j++){
       a1cov(i,j)=particles_.at(0).VertexCov()(i,j);
     }
   }

    refitParticles.push_back(LorentzVectorParticle(a1,a1cov,PDGInfo::tau_plus,c,b));

    TMatrixT<double> mu=ComputeTau2LorentzVectorPar(exppar_);
    mu(0,0)= particles_.at(1).Parameter(LorentzVectorParticle::vx);
    mu(1,0)= particles_.at(1).Parameter(LorentzVectorParticle::vy);
    mu(2,0)= particles_.at(1).Parameter(LorentzVectorParticle::vz);

    TMatrixTSym<double> mucov=ErrorMatrixPropagator::PropagateError(&ThreeProngThreeProngFitter::ComputeTau2LorentzVectorPar,exppar_,expcov_);
    for(int i=0; i<LorentzVectorParticle::NVertex; i++){
      for(int j=0; j<LorentzVectorParticle::NVertex; j++){
        mucov(i,j)=particles_.at(1).VertexCov()(i,j);
      }
    }
    
    refitParticles.push_back(LorentzVectorParticle(mu,mucov,PDGInfo::tau_minus,0.0,b));
   
  return refitParticles; 
}

LorentzVectorParticle
ThreeProngThreeProngFitter::GetMother(){
  UpdateExpandedPar();
  double c(0),b(0);
  for(unsigned int i=0;i<particles_.size();i++){c+=particles_.at(i).Charge();b=particles_.at(i).BField();}
  TMatrixT<double> m=ComputeMotherLorentzVectorPar(exppar_);
  TMatrixTSym<double> mcov=ErrorMatrixPropagator::PropagateError(&ThreeProngThreeProngFitter::ComputeMotherLorentzVectorPar,exppar_,expcov_);

  return LorentzVectorParticle(m,mcov,PDGInfo::Z0,c,b);
}

TVectorD
ThreeProngThreeProngFitter::HardValue(TVectorD &va,TVectorD &vb){
  TLorentzVector Tau1,Tau2;
  CovertParToObjects(va,vb,Tau1,Tau2);

  TLorentzVector z=Tau1+Tau2;
  TVectorD d(NConstraints());

  //d(0) = pow(Tau1.E() + Tau2.E(), 2.) - (Tau1.Vect() + Tau2.Vect()).Mag2()- pow(MassConstraint_,2.);
  d(0) = z.M() - MassConstraint_;
  return d;
}

TVectorD
ThreeProngThreeProngFitter::SoftValue(TVectorD &va,TVectorD &vb){
  TLorentzVector Tau1,Tau2;
  CovertParToObjects(va,vb,Tau1,Tau2);
  TVectorD d(NSoftConstraints());

  d(0) = Tau1.Px() + Tau2.Px() - RecoilX_;
  d(1) = Tau1.Py() + Tau2.Py() - RecoilY_;
  Logger(Logger::Debug) << "SCVec: " << d(0) << ", " << d(1) << std::endl;
  Logger(Logger::Debug) << "RecoilX_: " << RecoilX_ << ", RecoilY_: " << RecoilY_ << std::endl;
  return d;
} 

void ThreeProngThreeProngFitter::CovertParToObjects(TVectorD &va,TVectorD &vb,TLorentzVector &Tau1,TLorentzVector &Tau2){
  // Tau1=particles_.at(0).LV();//TLorentzVector(v(tau1_px),v(tau1_py),v(tau1_pz),sqrt(1.777*1.777+v(tau1_px)*v(tau1_px)+v(tau1_py)*v(tau1_py)+v(tau1_pz)*v(tau1_pz)));
  // Tau2=TLorentzVector(v(tau1_px),v(tau1_py),v(tau1_pz),sqrt(1.777*1.777+v(tau1_px)*v(tau1_px)+v(tau1_py)*v(tau1_py)+v(tau1_pz)*v(tau1_pz)));
  Tau1=TLorentzVector(va(tau_px),va(tau_py),va(tau_pz),sqrt(PDGInfo::tau_mass()*PDGInfo::tau_mass()+va(tau_px)*va(tau_px)+va(tau_py)*va(tau_py)+va(tau_pz)*va(tau_pz)));
  Tau2=TLorentzVector(vb(tau_px),vb(tau_py),vb(tau_pz),sqrt(PDGInfo::tau_mass()*PDGInfo::tau_mass()+vb(tau_px)*vb(tau_px)+vb(tau_py)*vb(tau_py)+vb(tau_pz)*vb(tau_pz)));
}

TMatrixT<double> ThreeProngThreeProngFitter::EstimateTauKinematic(TMatrixT<double> &inpar){
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
ThreeProngThreeProngFitter::ConfigureParameters(TrackParticle MuTrack, std::pair<double, double> phiAngle){
   TMatrixT<double>    outpar(5,1);
   outpar(0,0) = MuTrack.Parameter(TrackParticle::dxy);
   outpar(1,0) = MuTrack.Parameter(TrackParticle::phi);
   outpar(2,0) = MuTrack.Parameter(TrackParticle::lambda);
   outpar(3,0) = MuTrack.Parameter(TrackParticle::dz);
   outpar(4,0) = phiAngle.first;

   return outpar;
 }

TMatrixT<double>
ThreeProngThreeProngFitter::ConfigureParameterErrors(TrackParticle MuTrack, std::pair<double, double> phiAngle){
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
ThreeProngThreeProngFitter::ConfigureInitialAdvancedParameters(TrackParticle MuTrack,  TVector3 PV, std::vector< LorentzVectorParticle > TauThreeProngs){
  TMatrixT<double>    outpar;
  outpar.ResizeTo(9,1);
   outpar(0,0)  = MuTrack.Parameter(TrackParticle::dxy);
   outpar(1,0)  = MuTrack.Parameter(TrackParticle::phi);
   outpar(2,0)  = MuTrack.Parameter(TrackParticle::lambda);
   outpar(3,0)  = MuTrack.Parameter(TrackParticle::dz);
   outpar(4,0)  = PV.X();
   outpar(5,0)  = PV.Y();
   outpar(6,0)  = PV.Z();
   outpar(7,0) = particles_.at(0).LV().X();
   outpar(8,0) = particles_.at(0).LV().Y();

   return outpar;
 }

TMatrixT<double>
ThreeProngThreeProngFitter::ConfigureKinematicParameters(TMatrixT<double>  TauMuDir, std::vector< LorentzVectorParticle > TauThreeProngs){
  TMatrixT<double>    outpar;
  outpar.ResizeTo(5,1);
  outpar(0,0)  = TauMuDir(0,0);
  outpar(1,0)  = TauMuDir(1,0);
  outpar(2,0)  = particles_.at(0).LV().Px();
  outpar(3,0)  = particles_.at(0).LV().Py();
  outpar(4,0)  = particles_.at(0).LV().Pz();


  if(debug_){
    std::cout<<"EstimateTauKinematic  TauA1  "<<particles_.at(0).LV().Px()<<"  "<<particles_.at(0).LV().Py()<<"  "<<particles_.at(0).LV().Pz()<<std::endl;
    std::cout<<"ConfigureKinematicParameters:  outpar(0,0)  "<<outpar(0,0) <<std::endl;
    std::cout<<"ConfigureKinematicParameters:  outpar(1,0)  "<<outpar(1,0) <<std::endl;
    std::cout<<"ConfigureKinematicParameters:  outpar(2,0)  "<<outpar(2,0) <<std::endl;
    std::cout<<"ConfigureKinematicParameters:  outpar(3,0)  "<<outpar(3,0) <<std::endl;
    std::cout<<"ConfigureKinematicParameters:  outpar(4,0)  "<<outpar(4,0) <<std::endl;
  }
  return outpar;
 }

TMatrixTSym<double>
ThreeProngThreeProngFitter::ConfigureKinematicParameterErrors(TMatrixTSym<double>  TauMuDirError, std::vector< LorentzVectorParticle > TauThreeProngs){
  TMatrixTSym<double>    outpar;
  outpar.ResizeTo(5,5);
  
  for(int i=0; i<5; i++){
    for(int j=0; j<5; j++){
      
      if(i<2 && j < 2){
	outpar(i,j)=TauMuDirError(i,j);
      }else{
	outpar(i,j)=particles_.at(0).Covariance(i+1,j+1);
      }
    }
  }

  return outpar;
}

std::pair<double, double>
ThreeProngThreeProngFitter::EstimatePhiAngle( TVector3 dir, TVector3 dirE){
    std::pair<double, double> outpar;
  double phi = atan2(dir.Y(),dir.X());
  double deltaphi = sqrt(dirE.Y()*dirE.Y()/dir.X()/dir.X()+
 			 dirE.X()*dirE.X()*dir.Y()*dir.Y()/dir.X()/dir.X()/dir.X()/dir.X() )*cos(phi)*cos(phi);
  outpar = std::make_pair(phi, deltaphi);
  return outpar;
}

bool ThreeProngThreeProngFitter::isConverged(){
  if(!useFullRecoil_){
//	if(pardelta<MaxParDelta_ /*&& harddelta_vec.Norm1() < MaxHCDelta_ && chi2prev-chi2 < MaxChi2Delta_  && chi2prev>chi2*/){
		if(/*pardelta<MaxParDelta_ &&*/ harddelta_vec.Norm1() < MaxHCDelta_ && fabs(chi2prev-chi2) < MaxChi2Delta_ /*&& softdelta_vec.Norm1() < MaxSCDelta_ && chi2prev>chi2*/
      && (pow(para(0), 2.0) + pow(para(1), 2.0)) > ThreeProngs_.at(0).LV().Perp2()
      && (pow(parb(0), 2.0) + pow(parb(1), 2.0)) > ThreeProngs_.at(1).LV().Perp2()
    ){
    //	Logger(Logger::Verbose) << "converged " << delta << " chi2 " <<  chi2 << " chi2prev " << chi2prev <<"  Maxdelta  " <<MaxDelta_ <<std::endl;
	  return true;
	}
  }
  else{
	//if(pardelta<MaxParDelta_ && harddelta_vec.Norm1() < MaxHCDelta_ && softdelta_vec.Norm1() < MaxSCDelta_){
	// if(/*pardelta<MaxParDelta_ &&*/ harddelta_vec.Norm1() < MaxHCDelta_ && fabs(chi2prev-chi2) < MaxChi2Delta_ /*&& chi2prev>chi2*/){
	if(/*pardelta<MaxParDelta_ && harddelta_vec.Norm1() < MaxHCDelta_ && */fabs(chi2prev-chi2) < MaxChi2Delta_ /*&& chi2prev>chi2*/
      && (pow(para(0), 2.0) + pow(para(1), 2.0)) > ThreeProngs_.at(0).LV().Perp2()
      && (pow(parb(0), 2.0) + pow(parb(1), 2.0)) > ThreeProngs_.at(1).LV().Perp2()
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

TString ThreeProngThreeProngFitter::ParName(int par){
  switch (par){
    case tau1_px: return "tau1_px";
    case tau1_py: return "tau1_py";
    case tau1_pz: return "tau1_pz";
    case tau2_px: return "tau2_px";
    case tau2_py: return "tau2_py";
    case tau2_pz: return "tau2_pz";
    case lambda_1: return "lambda_1";
    case lambda_2: return "lambda_2";
    // add more/remove lambdas if number of hard constraints exceeds number of lambdas or vice versa
    default:
      Logger(Logger::Warning) << "Paramater Index " << par << " out of bounds !!";
      return "";
  }
}