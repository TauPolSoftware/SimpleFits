#include "TauPolSoftware/SimpleFits/interface/ThreeProngThreeProngFitter.h"
#include "TauPolSoftware/SimpleFits/interface/ChiSquareFunctionUpdator.h"
#include "TauPolSoftware/SimpleFits/interface/PDGInfo.h"
#include "TauPolSoftware/SimpleFits/interface/Logger.h"
#include "TDecompBK.h"
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

ThreeProngThreeProngFitter::ThreeProngThreeProngFitter(std::vector< LorentzVectorParticle > TauThreeProngs, std::vector< LorentzVectorParticle > ThreeProngs, TVector3 PVertex, TMatrixTSym<double> VertexCov) :
  ThreeProngThreeProngFitter(TauThreeProngs, ThreeProngs, PTObject(), PVertex, VertexCov, ThreeProngThreeProngFitter::MassConstraint_)
{
}

ThreeProngThreeProngFitter::ThreeProngThreeProngFitter(std::vector< LorentzVectorParticle > TauThreeProngs, std::vector< LorentzVectorParticle > ThreeProngs, PTObject ResPtEstimate, TVector3 PVertex, TMatrixTSym<double> VertexCov) :
  ThreeProngThreeProngFitter(TauThreeProngs, ThreeProngs, ResPtEstimate, PVertex, VertexCov, ThreeProngThreeProngFitter::MassConstraint_)
{
}

ThreeProngThreeProngFitter::ThreeProngThreeProngFitter(std::vector< LorentzVectorParticle > TauThreeProngs, std::vector< LorentzVectorParticle > ThreeProngs, PTObject ResPtEstimate, TVector3 PVertex, TMatrixTSym<double> VertexCov, double MassConstraint):
  LagrangeMultipliersFitter(),
  particles_(TauThreeProngs),
  particles0_(TauThreeProngs),
  ThreeProngs_(ThreeProngs),
  PV_(PVertex),
  RecoilX_(ResPtEstimate.Par()(0,0)),
  RecoilY_(ResPtEstimate.Par()(1,0)),
  ResPtEstimate_(ResPtEstimate),
  debug_(false),
  AnalyticalCovariance_(false)
{
  PVCov_.ResizeTo(VertexCov); PVCov_ = VertexCov;
  MassConstraint_ = MassConstraint;
  if(ResPtEstimate.isValid()) useFullRecoil_ = true;
  Configure();
}

void ThreeProngThreeProngFitter::Configure(){

  Logger(Logger::Debug) << "Tau 1 covariance: " << std::endl;
  if(Logger::Instance()->Level() == Logger::Debug){
    particles0_.at(0).getParMatrix().Print();
    particles0_.at(0).getCovMatrix().Print();
  }
  Logger(Logger::Debug) << "Tau 2 covariance: " << std::endl;
  if(Logger::Instance()->Level() == Logger::Debug){
    particles0_.at(1).getParMatrix().Print();
    particles0_.at(1).getCovMatrix().Print();
  }

  Logger(Logger::Debug) << "(TauThreeProngs.at(0).LV() + TauThreeProngs.at(1).LV()).M(): " << (particles0_.at(0).LV() + particles0_.at(1).LV()).M() << std::endl;

  isconfigured=false;

  int size = particles_.size()*3;
  int sizeTrunc = 3;
  TMatrixT<double>    inpar(size,1);
  TMatrixTSym<double> incov(size);

  // Get primary vertex information
  if(PVCov_.GetNrows() != LorentzVectorParticle::NVertex)
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
  for(int i=0; i<npartr;i++)
    parb_0(i) = PARb_0(i,0);
  y_.ResizeTo(npartr,1); y_ = convertToMatrix(para_0);
  z_.ResizeTo(npartr,1); z_ = convertToMatrix(parb_0);

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
  TMatrixT<double> Tau2par=ComputeTau2LorentzVectorPar(inpar);
  TMatrixT<double> Tau1par=ComputeTau1LorentzVectorPar(inpar);

  outpar(LorentzVectorParticle::px,0)=Tau2par(LorentzVectorParticle::px,0)+Tau1par(LorentzVectorParticle::px,0);
  outpar(LorentzVectorParticle::py,0)=Tau2par(LorentzVectorParticle::py,0)+Tau1par(LorentzVectorParticle::py,0);
  outpar(LorentzVectorParticle::pz,0)=Tau2par(LorentzVectorParticle::pz,0)+Tau1par(LorentzVectorParticle::pz,0);

  double Etau2sq=pow(Tau2par(LorentzVectorParticle::px,0),2.0)+pow(Tau2par(LorentzVectorParticle::py,0),2.0)+pow(Tau2par(LorentzVectorParticle::pz,0),2.0)+pow(Tau2par(LorentzVectorParticle::m,0),2.0);
  double Etau1sq=pow(Tau1par(LorentzVectorParticle::px,0),2.0)+pow(Tau1par(LorentzVectorParticle::py,0),2.0)+pow(Tau1par(LorentzVectorParticle::pz,0),2.0)+pow(Tau1par(LorentzVectorParticle::m,0),2.0);
  double P2=pow(outpar(LorentzVectorParticle::px,0),2.0)+pow(outpar(LorentzVectorParticle::py,0),2.0)+pow(outpar(LorentzVectorParticle::pz,0),2.0);
  outpar(LorentzVectorParticle::m,0)=sqrt(fabs(pow(sqrt(Etau2sq)+sqrt(Etau1sq),2.0)-P2));

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

  TMatrixT<double> tau1=ComputeTau1LorentzVectorPar(exppar_);
  TMatrixTSym<double> tau1cov=ErrorMatrixPropagator::PropagateError(&ThreeProngThreeProngFitter::ComputeTau1LorentzVectorPar,exppar_,expcov_);

   for(int i=0; i<LorentzVectorParticle::NVertex; i++){
     for(int j=0; j<LorentzVectorParticle::NVertex; j++){
       tau1cov(i,j)=particles_.at(0).VertexCov()(i,j);
     }
   }

  refitParticles.push_back(LorentzVectorParticle(tau1,tau1cov,particles_.at(0).Charge()*(PDGInfo::tau_plus),particles_.at(0).Charge(),particles_.at(0).BField()));

  TMatrixT<double> tau2=ComputeTau2LorentzVectorPar(exppar_);
  tau2(0,0)= particles_.at(1).Parameter(LorentzVectorParticle::vx);
  tau2(1,0)= particles_.at(1).Parameter(LorentzVectorParticle::vy);
  tau2(2,0)= particles_.at(1).Parameter(LorentzVectorParticle::vz);

  TMatrixTSym<double> tau2cov=ErrorMatrixPropagator::PropagateError(&ThreeProngThreeProngFitter::ComputeTau2LorentzVectorPar,exppar_,expcov_);
  for(int i=0; i<LorentzVectorParticle::NVertex; i++){
    for(int j=0; j<LorentzVectorParticle::NVertex; j++){
      tau2cov(i,j)=particles_.at(1).VertexCov()(i,j);
    }
  }

  refitParticles.push_back(LorentzVectorParticle(tau2,tau2cov,particles_.at(1).Charge()*(PDGInfo::tau_plus),particles_.at(1).Charge(),particles_.at(1).BField()));

  return refitParticles;
}

LorentzVectorParticle ThreeProngThreeProngFitter::GetMother(){
  UpdateExpandedPar();
  double c(0),b(0);
  for(unsigned int i=0;i<particles_.size();i++){c+=particles_.at(i).Charge();b=particles_.at(i).BField();}
  TMatrixT<double> m=ComputeMotherLorentzVectorPar(exppar_);
  TMatrixTSym<double> mcov=ErrorMatrixPropagator::PropagateError(&ThreeProngThreeProngFitter::ComputeMotherLorentzVectorPar,exppar_,expcov_);

  return LorentzVectorParticle(m,mcov,PDGInfo::Z0,c,b);
}

TVectorD ThreeProngThreeProngFitter::HardValue(TVectorD &va,TVectorD &vb){
  TLorentzVector Tau1,Tau2;
  CovertParToObjects(va,vb,Tau1,Tau2);

  TLorentzVector z=Tau1+Tau2;
  TVectorD d(NConstraints());

  //d(0) = pow(Tau1.E() + Tau2.E(), 2.) - (Tau1.Vect() + Tau2.Vect()).Mag2()- pow(MassConstraint_,2.);
  d(0) = z.M() - MassConstraint_;
  return d;
}

TVectorD ThreeProngThreeProngFitter::SoftValue(TVectorD &va,TVectorD &vb){
  TVectorD d(NSoftConstraints());
  if(useFullRecoil_){
    TLorentzVector Tau1,Tau2;
    CovertParToObjects(va,vb,Tau1,Tau2);
    d(0) = Tau1.Px() + Tau2.Px() - RecoilX_;
    d(1) = Tau1.Py() + Tau2.Py() - RecoilY_;
    Logger(Logger::Debug) << "SCVec: " << d(0) << ", " << d(1) << std::endl;
    // Logger(Logger::Debug) << "RecoilX_: " << RecoilX_ << ", RecoilY_: " << RecoilY_ << std::endl;
  }
  return d;
}

void ThreeProngThreeProngFitter::CovertParToObjects(TVectorD &va,TVectorD &vb,TLorentzVector &Tau1,TLorentzVector &Tau2){
  // Tau1=particles_.at(0).LV();//TLorentzVector(v(tau1_px),v(tau1_py),v(tau1_pz),sqrt(1.777*1.777+v(tau1_px)*v(tau1_px)+v(tau1_py)*v(tau1_py)+v(tau1_pz)*v(tau1_pz)));
  // Tau2=TLorentzVector(v(tau1_px),v(tau1_py),v(tau1_pz),sqrt(1.777*1.777+v(tau1_px)*v(tau1_px)+v(tau1_py)*v(tau1_py)+v(tau1_pz)*v(tau1_pz)));
  Tau1=TLorentzVector(va(tau_px),va(tau_py),va(tau_pz),sqrt(PDGInfo::tau_mass()*PDGInfo::tau_mass()+va(tau_px)*va(tau_px)+va(tau_py)*va(tau_py)+va(tau_pz)*va(tau_pz)));
  Tau2=TLorentzVector(vb(tau_px),vb(tau_py),vb(tau_pz),sqrt(PDGInfo::tau_mass()*PDGInfo::tau_mass()+vb(tau_px)*vb(tau_px)+vb(tau_py)*vb(tau_py)+vb(tau_pz)*vb(tau_pz)));
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

bool ThreeProngThreeProngFitter::ApplyLagrangianConstraints(){
  //  if(V_a.GetNrows()!=para.GetNrows()) V_a.ResizeTo(para.GetNrows(),para.GetNrows());
  if(A.GetNrows()!=NConstraints() || A.GetNcols()!=para.GetNrows()) A.ResizeTo(NConstraints(),para.GetNrows());
  if(B.GetNrows()!=NConstraints() || B.GetNcols()!=parb.GetNrows()) B.ResizeTo(NConstraints(),parb.GetNrows());
  if(Fa.GetNrows()!=NSoftConstraints() || Fa.GetNcols()!=para.GetNrows()) Fa.ResizeTo(NSoftConstraints(),para.GetNrows());
  if(Fb.GetNrows()!=NSoftConstraints() || Fb.GetNcols()!=parb.GetNrows()) Fb.ResizeTo(NSoftConstraints(),parb.GetNrows());

  // Setup intial values
	harddelta_vecprev.ResizeTo(NConstraints());
	harddelta_vecprev = HardValue(para_0,parb_0);
	softdelta_vecprev.ResizeTo(NSoftConstraints());
	softdelta_vecprev = SoftValue(para_0,parb_0);
	paraprev.ResizeTo(para_0); paraprev=para_0;
	parbprev.ResizeTo(parb_0); parbprev=parb_0;

  // Setup initial values II

  TMatrixT<double> y(y_);
  TMatrixT<double> z(z_);
  TMatrixT<double> a0 = convertToMatrix(para_0);
  TMatrixT<double> b0 = convertToMatrix(parb_0);
  TMatrixT<double> c0 = convertToMatrix(HardValue(para_0,parb_0));
  TMatrixT<double> f0 = convertToMatrix(SoftValue(para_0,parb_0));
  Fa = DerivativeSCa();
  Fb = DerivativeSCb();
  A  = DerivativeHCa();
  B  = DerivativeHCb();
  TMatrixT<double> FaT = Fa; FaT.T();
  TMatrixT<double> FbT = Fb; FbT.T();
  TMatrixT<double> AT = A; AT.T();
  TMatrixT<double> BT = B; BT.T();

  TMatrixTSym<double> V_a = cova_0;
  TMatrixTSym<double> V_b = covb_0;
  //TMatrixTSym<double> V_f=V_a;//ComputeV_f(V_a,V_b,para_0, parb_0);
  TMatrixTSym<double> V_f; V_f.ResizeTo(NSoftConstraints(),NSoftConstraints());

  V_f = ComputeV_f(V_a, V_b, para_0, parb_0);

  // if(Logger::Instance()->Level() == Logger::Debug){
  //   Logger(Logger::Debug) << "Jacobi Matrices Fa and Fb " << std::endl;
  //   Fa.Print();
  //   Fb.Print();
  //   Logger(Logger::Debug) << "Jacobi Matrices A and B " << std::endl;
  //   A.Print();
  //   B.Print();
  // }

  V_f.SetTol(1.e-50);
  if(!useFullRecoil_) V_f.Similarity(Fa);

  //----  fill final matrix blocks
  TMatrixTSym<double> V_a_inv = V_a;
  if( fabs(V_a_inv.Determinant())  < 1e-25){
       std::cout << "Fit failed: unable to invert, matrix is singular " << " \n" << std::endl;
       return false;
  }
  V_a_inv.Invert();

  TMatrixTSym<double> V_b_inv = V_b;
  if( fabs(V_b_inv.Determinant())  < 1e-25){
       std::cout << "Fit failed: unable to invert, matrix is singular " << " \n" << std::endl;
       return false;
  }
  V_b_inv.Invert();

  TMatrixTSym<double> V_f_inv = V_f;

  double detVf = V_f_inv.Determinant();
  if( fabs(detVf)  < 1e-25){
       std::cout << "Fit failed: unable to invert, matrix is singular " << detVf << " \n" << std::endl;
       return false;
  }
  V_f_inv.Invert();
  V_f_inv_.ResizeTo(V_f_inv);
  V_f_inv_ = V_f_inv;

  TMatrixT<double> M11 = V_a_inv + FaT*V_f_inv*Fa;
  TMatrixT<double> M12 = FaT*V_f_inv*Fb;
  TMatrixT<double> M21 = FbT*V_f_inv*Fa;
  TMatrixT<double> M22 = V_b_inv + FbT*V_f_inv*Fb;

  TMatrixT<double> V1 = V_a_inv*y - FaT*V_f_inv*(f0 - Fa*a0 - Fb*b0);
  TMatrixT<double> V2 = V_b_inv*z - FbT*V_f_inv*(f0 - Fa*a0 - Fb*b0);
  TMatrixT<double> V3 = A*a0 + B*b0  - c0;
  TMatrixT<double> M = MakeFullMatrix(M11,M12,M21,M22,A,B);
  TMatrixT<double> V = MakeFullVector(V1,V2,V3);

  // TMatrixDEigen  MEig(V_f);
  // std::cout<<" EigenValues "<<std::endl;(MEig.GetEigenValues()).Print();
  // std::cout<<" EigenVectors "<<std::endl;(MEig.GetEigenVectors()).Print();

  double detM = M.Determinant();

  if(fabs(detM)>1e40 or fabs(detM)  < 1e-25){
       Logger(Logger::Error) << "Fit failed: unable to invert SYM  matrix LARGE Determinant or Singular Matrix" << detM << " \n" << std::endl;
       return false;
  }
  TMatrixT<double> M_inv = M; M_inv.Invert();
  // solve equations
  TMatrixT<double> res = M_inv*V;

  TMatrixT<double> par_a = solutiona(res);
  TMatrixT<double> par_b = solutionb(res);
  TMatrixT<double> lambda = solutionlambda(res);

  lambda_.ResizeTo(NConstraints(),1);
  lambda_ = lambda;

  para = convertToVector(par_a);
  parb = convertToVector(par_b);

  chi2prev = chi2;
  chi2_vecprev.ResizeTo(chi2_vec); chi2_vecprev = chi2_vec;

  TVectorD Currentchi2_vec = ChiSquareUsingInitalPoint(par_a,par_b,lambda,V_f_inv);
  double Curentchi2(Currentchi2_vec.Sum()), Currentdelta(ConstraintDelta(para,parb));

  TMatrixT<double> a_s = par_a;
  TMatrixT<double> b_s = par_b;

  a_s = par_a;
  b_s = par_b;
  TVectorD a_v = convertToVector(a_s);
  TVectorD b_v = convertToVector(b_s);
  softdelta_vec.ResizeTo(NSoftConstraints());
  harddelta_vec.ResizeTo(NConstraints());
  softdelta_vec = SoftValue(a_v,b_v);
  harddelta_vec = HardValue(a_v,b_v);

  // set chi2
  chi2=Curentchi2;
  chi2_vec.ResizeTo(Currentchi2_vec); chi2_vec = Currentchi2_vec;
  //set delta
  delta=Currentdelta;
  para = convertToVector(a_s);
  parb = convertToVector(b_s);
  pardelta=0;
  for(int l=0;l<par_a.GetNrows();l++){
    pardelta+=fabs(a0(l,0) - a_s(l,0));
    pardelta+=fabs(b0(l,0) - b_s(l,0));
  }

  para_0  = convertToVector(a_s);
  parb_0 =  convertToVector(b_s);

  return true;
}

TVectorD ThreeProngThreeProngFitter::ChiSquareUsingInitalPoint(TMatrixT<double> a,TMatrixT<double> b,TMatrixT<double> lambda,TMatrixTSym<double> V_f_inv){
  TMatrixTSym<double> V_alpha0 = cova_0;
  TMatrixTSym<double> V_alpha0_inv = cova_0;
  V_alpha0_inv.ResizeTo(cova_0.GetNrows(),cova_0.GetNrows());
  TDecompBK InverterA(V_alpha0);
  if(!InverterA.Decompose()){ // handle rare case where inversion is not possible (ie assume diagonal)
    std::cout << "LagrangeMultipliersFitter::ChiSquareUsingInitalPoint: Error non-invertable Matrix... Calculating under assumption that correlations can be neglected!!!" << std::endl;
    for(int j = 0; j<par.GetNrows(); j++){
      for(int i = 0; i<par.GetNrows(); i++){
        if(i==j)
          V_alpha0_inv(i,j) = 1.0/V_alpha0(i,j);
        else
          V_alpha0_inv(i,j) = 0.0;
      }
    }
  }
  else{
    V_alpha0_inv=InverterA.Invert();
  }

  TMatrixTSym<double> V_beta0 = covb_0;
  TMatrixTSym<double> V_beta0_inv = covb_0;
  V_beta0_inv.ResizeTo(covb_0.GetNrows(),covb_0.GetNrows());
  TDecompBK InverterB(V_beta0);
  if(!InverterB.Decompose()){ // handle rare case where inversion is not possible (ie assume diagonal)
    std::cout << "LagrangeMultipliersFitter::ChiSquareUsingInitalPoint: Error non-invertable Matrix... Calculating under assumption that correlations can be neglected!!!" << std::endl;
    for(int j=0; j<par.GetNrows(); j++){
      for(int i=0; i<par.GetNrows(); i++){
        if(i==j)
          V_beta0_inv(i,j) = 1.0/V_beta0(i,j);
        else
          V_beta0_inv(i,j) = 0.0;
      }
    }
  }
  else{
    V_beta0_inv=InverterB.Invert();
  }

  TMatrixT<double> lambdaT = lambda; lambdaT.T();
  TMatrixT<double> a0 = convertToMatrix(para_0);
  TMatrixT<double> b0 = convertToMatrix(parb_0);
  TMatrixT<double> da = y_-a;
  TMatrixT<double> db = z_-b;

  TMatrixT<double> daT = da;  daT.T();
  TMatrixT<double> dbT = db;  dbT.T();
  TMatrixT<double> chisquare_vara = daT*(V_alpha0_inv*da);
  TMatrixT<double> chisquare_varb = dbT*(V_beta0_inv*db);
  TVectorT<double> a_v = convertToVector(a);
  TVectorT<double> b_v = convertToVector(b);
  TMatrixT<double> f = convertToMatrix(SoftValue(a_v,b_v));
  TMatrixT<double> fT = f; fT.T();

  TVectorD chi2(3);
  chi2(0) = chisquare_vara(0,0) + chisquare_varb(0,0);
  chi2(1) = (fT*V_f_inv*f)(0,0);
  chi2(2) = 2*(lambdaT*convertToMatrix(HardValue(a_v,b_v)))(0,0);

  Logger(Logger::Debug) << "chi2 contributions: " << chi2(0) << " (orig) + " << chi2(1) << " (SC) + " << chi2(2) << " (HC) = " << chi2.Sum() << std::endl;
  Logger(Logger::Debug) << "hard constraints contributions: " << lambdaT(0,0)*HardValue(a_v,b_v)(0) << " (Mass) = " << chi2(2)/2 << std::endl;

  return chi2;
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
