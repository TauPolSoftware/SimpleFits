#include "TauPolSoftware/SimpleFits/interface/LagrangeMultipliersFitter.h"
#include "TauPolSoftware/SimpleFits/interface/Logger.h"
#include "TauPolSoftware/SimpleFits/interface/ChiSquareFunctionUpdator.h"
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

LagrangeMultipliersFitter::LagrangeMultipliersFitter():
  isconfigured(false),
  isFit(false),
  fittingMode_(FittingProc::Minuit),
  epsilon_(0.001),
  weight_(1.0),
  MaxDelta_(0.1),
  nitermax_(100),
  MaxParDelta_(10.),
  MaxHCDelta_(1.0),
  MaxSCDelta_(50.),
  MaxChi2Delta_(0.1),
  nCutStepmax_(100),
  chi2(1e10),
  niter(0),
  D(1,1),
  V_D(1,1)
{

}

bool LagrangeMultipliersFitter::Fit(){
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
    LagrangeMultipliersFitterChiSquareFunctionUpdator updator(this);
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
    MnPar.SetLimits(5, -500, 500);

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
    Logger(Logger::Error) << "minimum: " << min << std::endl;
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



bool LagrangeMultipliersFitter::ApplyLagrangianConstraints(){
  //  if(V_a.GetNrows()!=para.GetNrows()) V_a.ResizeTo(para.GetNrows(),para.GetNrows());
  if(A.GetNrows()!=NConstraints() || A.GetNcols()!=para.GetNrows()) A.ResizeTo(NConstraints(),para.GetNrows());
  if(B.GetNrows()!=NConstraints() || B.GetNcols()!=parb.GetNrows()) B.ResizeTo(NConstraints(),parb.GetNrows());
  if(Fa.GetNrows()!=NSoftConstraints() || Fa.GetNcols()!=para.GetNrows()) Fa.ResizeTo(NSoftConstraints(),para.GetNrows());
  if(Fb.GetNrows()!=NSoftConstraints() || Fb.GetNcols()!=parb.GetNrows()) Fb.ResizeTo(NSoftConstraints(),parb.GetNrows());



  // Setup intial values
	harddelta_vecprev.ResizeTo(NConstraints());
	harddelta_vecprev=HardValue(para_0,parb_0);
	softdelta_vecprev.ResizeTo(NSoftConstraints());
	softdelta_vecprev=SoftValue(para_0,parb_0);
	paraprev.ResizeTo(para_0); paraprev=para_0;
	parbprev.ResizeTo(parb_0); parbprev=parb_0;

  // TMatrixT<double> alpha_A=convertToMatrix(par);
  // TMatrixT<double> alpha_0=convertToMatrix(par_0);

  // TMatrixT<double> delta_alpha_A=alpha_A-alpha_0;
  // Setup initial values II

  TMatrixT<double> y(y_);
  TMatrixT<double> a0 = convertToMatrix(para_0);
  TMatrixT<double> b0 = convertToMatrix(parb_0);
  TMatrixT<double> c0 = convertToMatrix(HardValue(para_0,parb_0));
  TMatrixT<double> f0 = convertToMatrix(SoftValue(para_0,parb_0));
  Fa =DerivativeSCa();
  Fb =DerivativeSCb();
  A  =DerivativeHCa();
  B  =DerivativeHCb();
  TMatrixT<double> FaT=Fa; FaT.T();
  TMatrixT<double> FbT=Fb; FbT.T();
  TMatrixT<double> AT=A; AT.T();
  TMatrixT<double> BT=B; BT.T();


  TMatrixTSym<double> V_a=cova_0;
  TMatrixTSym<double> V_b=covb_0;
  //TMatrixTSym<double> V_f=V_a;//ComputeV_f(V_a,V_b,para_0, parb_0);
  TMatrixTSym<double> V_f; V_f.ResizeTo(NSoftConstraints(),NSoftConstraints());

  V_f=ComputeV_f(V_a,V_b,para_0, parb_0);
 

  if(Logger::Instance()->Level() == Logger::Debug){
	Logger(Logger::Debug) << "Jacobi Matrices Fa and Fb " << std::endl;
	Fa.Print();
	Fb.Print();
	Logger(Logger::Debug) << "Jacobi Matrices A and B " << std::endl;
	A.Print();
	B.Print();
  }


  V_f.SetTol(1.e-50);
  if(!useFullRecoil_) V_f.Similarity(Fa);


  //----  fill final matrix blocks

  TMatrixTSym<double> V_a_inv= V_a;
  double detVa = V_a_inv.Determinant();
  if( fabs(detVa)  < 1e-25){
       std::cout << "Fit failed: unable to invert, V_a matrix is singular. Determinant: " << detVa << "\n" << std::endl;
       V_a.Print();
       return false;
  } V_a_inv.Invert();




  TMatrixTSym<double> V_f_inv= V_f;

  double detVf = V_f_inv.Determinant();
  if( fabs(detVf)  < 1e-25){
       std::cout << "Fit failed: unable to invert, V_f matrix is singular Determinant: " << detVf << " \n" << std::endl;
       V_f.Print();
       return false;
  } V_f_inv.Invert();
  V_f_inv_.ResizeTo(V_f_inv);
  V_f_inv_ = V_f_inv;

  TMatrixT<double> M11 = V_a_inv + FaT*V_f_inv*Fa;
  TMatrixT<double> M12 = FaT*V_f_inv*Fb;
  TMatrixT<double> M21 = FbT*V_f_inv*Fa;
  TMatrixT<double> M22 = FbT*V_f_inv*Fb;

  TMatrixT<double> V1 = V_a_inv*y - FaT*V_f_inv*(f0 - Fa*a0 - Fb*b0);
  TMatrixT<double> V2 =  FbT*V_f_inv*(Fa*a0 + Fb*b0 - f0);
  TMatrixT<double> V3 = A*a0 + B*b0  - c0;
  TMatrixT<double> M = MakeFullMatrix(M11,M12,M21,M22,A,B);
  TMatrixT<double> V = MakeFullVector(V1,V2,V3);


  TMatrixDEigen  MEig(V_f);
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
  TMatrixT<double> lambda=solutionlambda(res);

  lambda_.ResizeTo(NConstraints(),1);
  lambda_=lambda;


  para  = convertToVector(par_a);
  parb  = convertToVector(par_b);


  // do while loop to see if the convergance criteria are satisfied
  //double s(1), stepscale(0.05);
  chi2prev=chi2;
  chi2_vecprev.ResizeTo(chi2_vec); chi2_vecprev = chi2_vec;

  TVectorD Currentchi2_vec = ChiSquareUsingInitalPoint(par_a,par_b,lambda,V_f_inv);
  double Curentchi2(Currentchi2_vec.Sum()), Currentdelta(ConstraintDelta(para,parb));

  TMatrixT<double> a_s=par_a;
  TMatrixT<double> b_s=par_b;
  // convergence in 2 step procedure to minimize chi2 within MaxDelta_ of the constriants
  // 1) Get within 5x MaxDelta_
  // 2) converge based on improving chi2 and constrianed delta
  //unsigned int Proc=ConstraintMin;
  //if(ConstraintDelta(para,parb)<5*MaxDelta_)Proc=Chi2AndConstaintMin;
  //int  NIter=(int)(1.0/stepscale);



  // for(int iter=0;iter<NIter;iter++){

  //   // compute safty cutoff for numberical constraint
  //   double diff=0;
  //   for(int l=0;l<par_a.GetNrows();l++){
  //     if(diff<y(l,0)-a_s(l,0) + b_s(l,0) - b0(l,0))diff=y(l,0)-a_s(l,0) +  b_s(l,0) - b0(l,0);
  //   }
  //   double delta_s=ConstraintDelta(convertToVector(a_s),convertToVector(b_s));
  //   std::cout<<"  iteration number ---> "<< iter <<std::endl;
  //   std::cout<<"diff  "<< diff <<std::endl;
  //   std::cout<<"delta_s  "<< delta_s << "  CurrentDelta  " <<Currentdelta<<std::endl;
  //   if(Proc==ConstraintMin){
  //     if(delta_s<Currentdelta || iter==NIter || diff<100*epsilon_){Curentchi2=ChiSquareUsingInitalPoint(y,a_s,b_s,lambda); Currentdelta=delta_s; ScaleFactor=s; break;}
  //   }

  //   else if(Proc==Chi2AndConstaintMin){
  //     double chi2_s=ChiSquareUsingInitalPoint(y,a_s,b_s,lambda);
  //     std::cout<<"chi2_s  "<< chi2_s << "  Curentchi2  " <<Curentchi2<<std::endl;

  //     if((delta_s<Currentdelta/*+MaxDelta_*/ && chi2_s<Curentchi2) || iter==NIter || diff<100*epsilon_){Curentchi2=chi2_s; Currentdelta=delta_s; ScaleFactor=s; break;}
  //   }
  //   s-=stepscale;
  //   //    alpha_s=alpha_A+s*(alpha-alpha_A);

  //   a_s = par_a + s*(y - par_a);
  //   b_s = par_b + s*(b0 - par_b);
  // }

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
    //pardelta=y(l,0)-a_s(l,0) +  b_s(l,0) - b0(l,0);
	pardelta+=fabs(a0(l,0) - a_s(l,0));
	pardelta+=fabs(b0(l,0) - b_s(l,0));
  }

  para_0  = convertToVector(a_s);
  parb_0 =  convertToVector(b_s);

  return true;
}
TMatrixD LagrangeMultipliersFitter::DerivativeHCa(){ // always evaluated at current par
  TMatrixD Derivatives(NConstraints(),para_0.GetNrows());
  TVectorD para_plus(para_0.GetNrows());
  TVectorD value(NConstraints());
  TVectorD value_plus(NConstraints());
  for(int j=0;j<para_0.GetNrows();j++){
    for(int i=0;i<para_0.GetNrows();i++){
      para_plus(i)=para_0(i);
       if(i==j) para_plus(i)=para_0(i)+epsilon_;
    }
    value=HardValue(para_0,parb_0);
    value_plus=HardValue(para_plus,parb_0);
    for(int i=0; i<NConstraints();i++){
      Derivatives(i,j)=(value_plus(i)-value(i))/epsilon_;
    }
  }
  return Derivatives;
}

TMatrixD LagrangeMultipliersFitter::DerivativeHCb(){ // always evaluated at current par
  TMatrixD Derivatives(NConstraints(),parb_0.GetNrows());
  TVectorD parb_plus(parb_0.GetNrows());
  TVectorD value(NConstraints());
  TVectorD value_plus(NConstraints());
  for(int j=0;j<parb_0.GetNrows();j++){
    for(int i=0;i<parb_0.GetNrows();i++){
      parb_plus(i)=parb_0(i);
         if(i==j) parb_plus(i)=parb_0(i)+epsilon_;
    }
    value=HardValue(para_0,parb_0);
    value_plus=HardValue(para_0,parb_plus);
    for(int i=0; i<NConstraints();i++){
       Derivatives(i,j)=(value_plus(i)-value(i))/epsilon_;
    }
  }
  return Derivatives;
}


TMatrixD LagrangeMultipliersFitter::DerivativeSCa(){ // always evaluated at current par
  TMatrixD Derivatives(NSoftConstraints(),para_0.GetNrows());
  TVectorD para_plus(para_0.GetNrows());

  TVectorD paraFD(para_0.GetNrows());

  if(!useFullRecoil_){
	paraFD(0) = para_0(0);
	paraFD(1) = para_0(1);
	paraFD(2) = atan(para_0(1)/para_0(0));
  }
  else paraFD = para_0;

  TVectorD paraFD_plus(para_0.GetNrows());

  TVectorD value(NSoftConstraints());
  TVectorD value_plus(NSoftConstraints());
  for(int j=0;j<paraFD.GetNrows();j++){
    for(int i=0;i<paraFD.GetNrows();i++){
      paraFD_plus(i)=paraFD(i);
       if(i==j) paraFD_plus(i)=paraFD(i)+epsilon_;
    }
    value=SoftValue(paraFD,parb_0);
    value_plus=SoftValue(paraFD_plus,parb_0);
    for(int i=0; i<NSoftConstraints();i++){
      Derivatives(i,j)=(value_plus(i)-value(i))/epsilon_;
    }
  }
  if(!useFullRecoil_) Derivatives(2,2) =1;
  return Derivatives;
}

TMatrixD LagrangeMultipliersFitter::DerivativeSCb(){ // always evaluated at current par
  TMatrixD Derivatives(NSoftConstraints(),parb_0.GetNrows());
  TVectorD parb_plus(parb_0.GetNrows());
  TVectorD value(NSoftConstraints());
  TVectorD value_plus(NSoftConstraints());
  for(int j=0;j<parb_0.GetNrows();j++){
    for(int i=0;i<parb_0.GetNrows();i++){
      parb_plus(i)=parb_0(i);
       if(i==j) parb_plus(i)=parb_0(i)+epsilon_;
    }
    value=SoftValue(para_0,parb_0);
    value_plus=SoftValue(para_0,parb_plus);
    for(int i=0; i<NSoftConstraints();i++){
      Derivatives(i,j)=(value_plus(i)-value(i))/epsilon_;
    }
  }
  if(!useFullRecoil_) Derivatives(2,2) =1;
   // Derivatives(0,2) =1/Derivatives(2,0);
   // Derivatives(1,2) =1/Derivatives(2,1);

  return Derivatives;
}


bool LagrangeMultipliersFitter::isConverged(){
  if(!useFullRecoil_){
//	if(pardelta<MaxParDelta_ /*&& harddelta_vec.Norm1() < MaxHCDelta_ && chi2prev-chi2 < MaxChi2Delta_  && chi2prev>chi2*/){
		if(/*pardelta<MaxParDelta_ &&*/ harddelta_vec.Norm1() < MaxHCDelta_ && fabs(chi2prev-chi2) < MaxChi2Delta_ /*&& softdelta_vec.Norm1() < MaxSCDelta_ && chi2prev>chi2*/){

    //	Logger(Logger::Verbose) << "converged " << delta << " chi2 " <<  chi2 << " chi2prev " << chi2prev <<"  Maxdelta  " <<MaxDelta_ <<std::endl;
	  return true;
	}
  }
  else{
	//if(pardelta<MaxParDelta_ && harddelta_vec.Norm1() < MaxHCDelta_ && softdelta_vec.Norm1() < MaxSCDelta_){
	// if(/*pardelta<MaxParDelta_ &&*/ harddelta_vec.Norm1() < MaxHCDelta_ && fabs(chi2prev-chi2) < MaxChi2Delta_ /*&& chi2prev>chi2*/){
	if(/*pardelta<MaxParDelta_ && harddelta_vec.Norm1() < MaxHCDelta_ && */fabs(chi2prev-chi2) < MaxChi2Delta_ /*&& chi2prev>chi2*/){
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



TVectorT<double> LagrangeMultipliersFitter::convertToVector(TMatrixT<double> M){
  TVectorT<double> V(M.GetNrows());
  for(int i=0; i<M.GetNrows();i++){
    V(i)=M(i,0);
    // std::cout<<"finPar   "<<V(i)<<std::endl;
  }
  return V;
}


TMatrixT<double> LagrangeMultipliersFitter::convertToMatrix(TVectorT<double> V){
  TMatrixT<double> M(V.GetNrows(),1);
  for(int i=0; i<M.GetNrows();i++){
    M(i,0)=V(i);
  }
  return M;
}


TMatrixT<double>
LagrangeMultipliersFitter::MakeFullVector(TMatrixT<double> V1,TMatrixT<double> V2,TMatrixT<double> V3){
 TMatrixT<double> M(V1.GetNrows()+V2.GetNrows()+V3.GetNrows(),1);

 int offset =0;
 for(int row =0; row < V1.GetNrows(); row++){  M(row,0) = V1(row,0); }

 offset =V1.GetNrows();
 for(int row =0; row < V2.GetNrows(); row++){  M(row+offset,0) = V2(row,0); }

 offset =V1.GetNrows()+V2.GetNrows();
 for(int row =0; row < V3.GetNrows(); row++){  M(row+offset,0) = V3(row,0); }

 return M;
}


TMatrixT<double>
LagrangeMultipliersFitter::MakeFullMatrix(TMatrixT<double> M11,TMatrixT<double> M12,TMatrixT<double> M21,TMatrixT<double> M22,TMatrixT<double> A,TMatrixT<double> B){

  TMatrixT<double> AT=A; AT.T();
  TMatrixT<double> BT=B; BT.T();
  TMatrixT<double> M(M11.GetNrows() +M21.GetNrows() +A.GetNrows(),M11.GetNcols() +M12.GetNcols() +AT.GetNcols());


  int offsetcols=0;int offsetrows=0;
  for(int row =0; row < M11.GetNrows(); row++){for(int col =0; col < M11.GetNcols(); col++){M(row + offsetrows,col +offsetcols)  = M11(row,col);}}
  offsetcols=0; offsetrows=M11.GetNrows();
  for(int row =0; row < M21.GetNrows(); row++){for(int col =0; col < M21.GetNcols(); col++){M(row + offsetrows,col +offsetcols)  = M21(row,col);}}
  offsetcols=0; offsetrows=M11.GetNrows() + M21.GetNrows();
  for(int row =0; row < A.GetNrows(); row++){for(int col =0; col < A.GetNcols(); col++){M(row + offsetrows,col +offsetcols)  = A(row,col);}}

  offsetcols=M11.GetNcols(); offsetrows=0;
  for(int row =0; row < M12.GetNrows(); row++){for(int col =0; col < M12.GetNcols(); col++){M(row + offsetrows,col +offsetcols)  = M12(row,col);}}
  offsetcols=M11.GetNcols(); offsetrows=M12.GetNrows();
  for(int row =0; row < M22.GetNrows(); row++){for(int col =0; col < M22.GetNcols(); col++){M(row + offsetrows,col +offsetcols)  = M22(row,col);}}
  offsetcols=M11.GetNcols(); offsetrows=M12.GetNrows() + M22.GetNrows();
  for(int row =0; row < B.GetNrows(); row++){for(int col =0; col < B.GetNcols(); col++){M(row + offsetrows,col +offsetcols)  = B(row,col);}}

  offsetcols=M11.GetNcols() + M12.GetNcols(); offsetrows=0;
  for(int row =0; row < AT.GetNrows(); row++){for(int col =0; col < AT.GetNcols(); col++){M(row + offsetrows,col +offsetcols)  = AT(row,col);}}
  offsetcols=M11.GetNcols() + M12.GetNcols(); offsetrows=AT.GetNrows();
  for(int row =0; row < BT.GetNrows(); row++){for(int col =0; col < BT.GetNcols(); col++){M(row + offsetrows,col +offsetcols)  = BT(row,col);}}
  offsetcols=M11.GetNcols() + M12.GetNcols(); offsetrows=AT.GetNrows() + BT.GetNrows();
  for(int row =0; row < NConstraints(); row++){for(int col =0; col < NConstraints(); col++){M(row + offsetrows,col +offsetcols)  = 0;}}

  return M;
}

TMatrixT<double>
LagrangeMultipliersFitter::solutiona(TMatrixT<double> M){
  TMatrixT<double> outpar(para_0.GetNrows(),1);
  for(int row =0; row < outpar.GetNrows(); row++){outpar(row,0)  = M(row,0);}
  return outpar;
}
TMatrixT<double>
LagrangeMultipliersFitter::solutionb(TMatrixT<double> M){
  int offsetrow = para_0.GetNrows(); //int offsetcol = para_0.GetNrows();
  TMatrixT<double> outpar(parb_0.GetNrows(),1);
  for(int row =0; row < outpar.GetNrows(); row++){outpar(row,0)  = M(row + offsetrow,0);}
  return outpar;
}

TMatrixT<double>
LagrangeMultipliersFitter::solutionlambda(TMatrixT<double> M){
  int offsetrow = para_0.GetNrows()+parb_0.GetNrows(); //int offsetcol = para_0.GetNrows() + parb_0.GetNrows();
  TMatrixT<double> outpar(NConstraints(),1);
  for(int row =0; row < outpar.GetNrows(); row++){outpar(row,0)  = M(row + offsetrow,0);}
  return outpar;
}
double LagrangeMultipliersFitter::ChiSquare(TMatrixT<double> delta_alpha,TMatrixT<double> lambda,TMatrixT<double> D, TMatrixT<double> d){
  TMatrixT<double> lambdaT=lambda; lambdaT.T();
  TMatrixT<double> chisquare=lambdaT*(D*delta_alpha+d);

  //  std::cout << "chi2  " <<chisquare(0,0) <<std::endl;
  double c2=chisquare(0,0);
  return c2;
}
TVectorD LagrangeMultipliersFitter::ChiSquareUsingInitalPoint(TMatrixT<double> a,TMatrixT<double> b,TMatrixT<double> lambda,TMatrixTSym<double> V_f_inv){
  // if(cova_0.GetNrows()!=V_alpha0_inv.GetNrows()){
  TMatrixTSym<double> V_alpha0(cova_0);
  TMatrixTSym<double> V_alpha0_inv(cova_0);
  TDecompBK Inverter(V_alpha0);
  if(!Inverter.Decompose()){ // handle rare case where inversion is not possible (ie assume diagonal)
    Logger(Logger::Error) << "non-invertable Matrix... Calculating under assumption that correlations can be neglected!!!" << std::endl;
    for(int j=0;j<cova_0.GetNrows();j++){
      for(int i=0;i<cova_0.GetNcols();i++){
        if(i==j) V_alpha0_inv(i,j)=1.0/V_alpha0(i,j);
        else V_alpha0_inv(i,j)=0.0;
      }
    }
  }
  else{
    V_alpha0_inv=Inverter.Invert();
  }
  // V_alpha0_inv=V_alpha0;V_alpha0_inv.Invert();

  // std::cout<<"V_alpha0  "<<std::endl;V_alpha0.Print();
  // std::cout<<"V_alpha0_inv  "<<std::endl;V_alpha0_inv.Print();
  TMatrixT<double> lambdaT=lambda; lambdaT.T();
  TMatrixT<double> a0=convertToMatrix(para_0);
  TMatrixT<double> b0=convertToMatrix(parb_0);
  TMatrixT<double> da=y_-a;

  TMatrixT<double> daT=da;  daT.T();
  TMatrixT<double> chisquare_var=daT*(V_alpha0_inv*da);
  TVectorT<double> a_v=convertToVector(a);
  TVectorT<double> b_v=convertToVector(b);
  TMatrixT<double> f = convertToMatrix(SoftValue(a_v,b_v));
  TMatrixT<double> fT = f; fT.T();
  TMatrixT<double> Fa =DerivativeSCa();
  TMatrixT<double> FaT=Fa; FaT.T();

  TMatrixT<double> chisquare_constraints(1,1);
  chisquare_constraints=lambdaT*convertToMatrix(HardValue(a_v,b_v)) + fT*V_f_inv*f;
  double c2=chisquare_var(0,0)+chisquare_constraints(0,0);

  TVectorD chi2(3);
  chi2(0) = chisquare_var(0,0);
  chi2(1) = (fT*V_f_inv*f)(0,0);
  chi2(2) = 2*(lambdaT*convertToMatrix(HardValue(a_v,b_v)))(0,0);

  Logger(Logger::Debug) << "chi2 comparison: " << c2 << ", " << chi2.Sum() << std::endl;
  Logger(Logger::Debug) << "chi2 contributions: " << chi2(0) << " (orig) + " << chi2(1) << " (SC) + " << chi2(2) << " (HC) = " << chi2.Sum() << std::endl;
  Logger(Logger::Debug) << "hard constraints contributions: " << lambdaT(0,0)*HardValue(a_v,b_v)(0) << " (Mass) + " << lambdaT(0,1)*HardValue(a_v,b_v)(1) << " (Theta) = " << chi2(2)/2 << std::endl;

  return chi2;

}

double LagrangeMultipliersFitter::UpdateChisquare(TVectorD a,TVectorD b){
  para_0 = a;
  parb_0 = b;
  TMatrixT<double> aM = convertToMatrix(a);
  TMatrixT<double> bM = convertToMatrix(b);
  ApplyLagrangianConstraints();
  // Logger(Logger::Info) << "y_.GetNrows(): " << y_.GetNrows() << " y_.GetNcols(): " << y_.GetNcols() << std::endl;
  // Logger(Logger::Info) << "aM.GetNrows(): " << aM.GetNrows() << " aM.GetNcols(): " << aM.GetNcols() << std::endl;
  // Logger(Logger::Info) << "bM.GetNrows(): " << bM.GetNrows() << " bM.GetNcols(): " << bM.GetNcols() << std::endl;
  // Logger(Logger::Info) << "lambda_.GetNrows(): " << lambda_.GetNrows() << " lambda_.GetNcols(): " << lambda_.GetNcols() << std::endl;
  // Logger(Logger::Info) << "V_f_inv_.GetNrows(): " << V_f_inv_.GetNrows() << " V_f_inv_.GetNcols(): " << V_f_inv_.GetNcols() << std::endl;
  TVectorD chi2vec=ChiSquareUsingInitalPoint(aM,bM,lambda_,V_f_inv_);
  return chi2vec.Sum();
}

// double LagrangeMultipliersFitter::ChiSquareUsingInitalPoint(TMatrixT<double> alpha,TMatrixT<double> lambda){
//   if(cov_0.GetNrows()!=V_alpha0_inv.GetNrows()){
//     TMatrixTSym<double> V_alpha0=cov_0;
//     V_alpha0_inv.ResizeTo(cov_0.GetNrows(),cov_0.GetNrows());
//     TDecompBK Inverter(V_alpha0);
//     if(!Inverter.Decompose()){ // handle rare case where inversion is not possible (ie assume diagonal)
//       std::cout << "LagrangeMultipliersFitter::ChiSquareUsingInitalPoint: Error non-invertable Matrix... Calculating under assumption that correlations can be neglected!!!" << std::endl;
//       for(int j=0;j<par.GetNrows();j++){
// 	for(int i=0;i<par.GetNrows();i++){
// 	  if(i==j) V_alpha0_inv(i,j)=1.0/V_alpha0(i,j);
// 	  else V_alpha0_inv(i,j)=0.0;
// 	}
//       }
//     }
//     else{
//       V_alpha0_inv=Inverter.Invert();
//     }
//   }

//   TMatrixT<double> lambdaT=lambda; lambdaT.T();
//   TMatrixT<double> alpha_0=convertToMatrix(par_0);
//   TMatrixT<double> dalpha=alpha-alpha_0;
//   TMatrixT<double> dalphaT=dalpha;  dalphaT.T();
//   TMatrixT<double> chisquare_var=dalphaT*(V_alpha0_inv*dalpha);
//   TVectorT<double> alpha_v=convertToVector(alpha);
//   //  std::cout<<"call value to compute chi2   "<<std::endl;
//   TMatrixT<double> chisquare_constraints=lambdaT*convertToMatrix(Value(alpha_v));
//   double c2=chisquare_var(0,0)+chisquare_constraints(0,0);
//   return c2;
// }

double LagrangeMultipliersFitter::ConstraintDelta(TVectorT<double> a,TVectorT<double>  b){
  //  std::cout<<"call value to compute ConstraintDelta   "<<std::endl;
  TVectorD dh_par=HardValue(a,b);
  TVectorD ds_par=SoftValue(a,b);
  //  TVectorD ds_par=converToMatrix(SoftValue(convertToVector(a),convertToVector(b)));
  double delta_d(0);
  double delta_dNew = dh_par.Norm1() + ds_par.Norm1();
  for(int i = 0; i<dh_par.GetNrows(); i++){
    delta_d+=fabs(dh_par(i)) + fabs(ds_par(i));
  }
  if(!useFullRecoil_) return delta_d;
  else return delta_dNew;
}

TMatrixT<double>  LagrangeMultipliersFitter::ComputeVarianceb(){

  // TMatrixTSym<double> V_0=covb_0;
  // TMatrixTSym<double> DTV_DD=V_0.SimilarityT(Fa);DTV_DD.Invert();
  // TMatrixT<double> FbT = Fb; FbT.T();
  // TMatrixT<double> CovCor=FbT*DTV_DD*Fb;CovCor.Invert();


  // for(int i=0; i<covb.GetNrows();i++){
  //   for(int j=0; j<=i;j++){
  //     covb(i,j)=CovCor(i,j);
  //   }
  // }

  //  return covb;

return covb_0;
}

TMatrixT<double>  LagrangeMultipliersFitter::ComputeVariancea(){



  // TMatrixTSym<double> V_0=cova_0;
  // TMatrixTSym<double> DTV_DD=V_0.SimilarityT(Fa);DTV_DD.Invert();
  // TMatrixT<double> FbT = Fb; FbT.T();
  // TMatrixT<double> Vz=FbT*DTV_DD*Fb;
  // TMatrixT<double> Vlambda = DTV_DD - DTV_DD*Vz*DTV_DD;
  // TMatrixT<double> Valpha = V_0 - V_0*(V_0.SimilarityT(Fa).T())*V_0;


  // for(int i=0; i<cova.GetNrows();i++){
  //   for(int j=0; j<=i;j++){
  //     cova(i,j)=Valpha(i,j);
  //   }
  // }

  return cova_0;
}
void  LagrangeMultipliersFitter::Print(TMatrixT<double> M){


   for(int str =0; str < M.GetNrows(); str++){
     for(int kol =0; kol < M.GetNcols(); kol++){
       std::cout<<"  "<< M(str,kol)<<"  ";

     }
     std::cout<<std::endl;
   }


}
TMatrixTSym<double> LagrangeMultipliersFitter::ComputeV_f(TMatrixTSym<double>  ca, TMatrixTSym<double>  cb, TVectorD pa,TVectorD pb){
  TMatrixTSym<double> Vf;
  if(!useFullRecoil_){
	Vf.ResizeTo(ca.GetNrows(), ca.GetNcols());
	Vf = ca;


	TMatrixTSym<double> Vfa;
	Vfa.ResizeTo(ca.GetNrows(), ca.GetNcols());
	Vfa = ca;

	TMatrixTSym<double> Vfb;
	Vfb.ResizeTo(cb.GetNrows(), cb.GetNcols());
	Vfb = cb;

	TMatrixDEigen  MEigca(ca);

	// double dgdx = -(pa(1) + pb(1))/(pow(pa(0) + pb(0),2) + pow(pa(1) + pb(1),2));
	// double dgdya = (pa(0) + pb(0))/(pow(pa(0) + pb(0),2) + pow(pa(1) + pb(1),2));
	// double dgdyb = dgdya;
	double dgdx=-(pa(1) + pb(1))/(pow(pa(0) + pb(0),2));
	double dgdya=pb(1)/(pa(0) + pb(0));
	double dgdyb=pa(1)/(pa(0) + pb(0));

	double deltaxxa = ca(0,0);
	double deltayya = ca(1,1);
	double deltaxya = ca(0,1);

	double deltaxxb = cb(0,0);
	double deltayyb = cb(1,1);
	double deltaxyb = cb(0,1);


	Vfa(2,0) = dgdx*deltaxxa + dgdya*deltaxya; Vfa(0,2) = Vfa(2,0);
	Vfa(2,1) = dgdx*deltaxya + dgdya*deltayya; Vfa(1,2) = Vfa(2,1);
	Vfa(2,2) = pow(dgdx,2)*deltaxxa + pow(dgdya,2)*deltayya + 2*dgdx*dgdya*deltaxya;

	Vfb(2,0) = dgdx*deltaxxb + dgdyb*deltaxyb; Vfb(0,2) = Vfb(2,0);
	Vfb(2,1) = dgdx*deltaxyb + dgdyb*deltayyb; Vfb(1,2) = Vfb(2,1);
	Vfb(2,2) = (pow(dgdx,2)*deltaxxb + pow(dgdyb,2)*deltayyb + 2*dgdx*dgdyb*deltaxyb);

	if(Logger::Instance()->Level() == Logger::Debug){
	  Logger(Logger::Debug) << "Vf matrix" << std::endl;
	  Vf.Print();
	}

	Vf = Vfa + Vfb;
  }
  else{
	Vf.ResizeTo(NSoftConstraints(),NSoftConstraints());

	TMatrixT<double> Jacobi(NSoftConstraints(),par.GetNrows());
	TMatrixTSym<double> fullcov = cov; //copy of the cov matrix as similarity overrides the matrix

	/*
	double Resonance_Pt = sqrt(pow(para_0(0) + parb_0(0), 2.) + pow(para_0(1) + parb_0(1), 2.));
	for(unsigned i=0; i<para_0.GetNrows();i++){
	  if(i==2) {
		Jacobi(0,i) = 0; //derivates of the z component are zero for the soft constraints
		Jacobi(1,i) = 0;
		Jacobi(0,i+3) = 0;
		Jacobi(1,i+3) = 0;
	  }
	  else{
		Jacobi(0,i) = (para_0(i) + parb_0(i))/Resonance_Pt;
		Jacobi(0,i+3) = (para_0(i) + parb_0(i))/Resonance_Pt;
	  }
	}
	Jacobi(1,0) = (para_0(1) + parb_0(1))/pow(para_0(0) + parb_0(0), 2.); Jacobi(1,3) = Jacobi(1,0);
	Jacobi(1,1) = parb_0(1)/(para_0(0) + parb_0(0));
	Jacobi(1,4) = para_0(1)/(para_0(0) + parb_0(0));
	*/
	Jacobi(0,0) = 1;
	Jacobi(0,3) = 1;
	Jacobi(1,1) = 1;
	Jacobi(1,4) = 1;
	Vf = fullcov.Similarity(Jacobi);
	if(Logger::Instance()->Level() == Logger::Debug){
	  Logger(Logger::Debug) << "Jacobi matrix" << std::endl;
	  Jacobi.Print();
	}
	if(Logger::Instance()->Level() == Logger::Debug){
	  Logger(Logger::Debug) << "Vf matrix" << std::endl;
	  Vf.Print();
	}
  }
  return Vf;

}

TMatrixTSym<double>  LagrangeMultipliersFitter::ScaleMatrix(TMatrixTSym<double>  M, double scale){
  TMatrixTSym<double> out(M.GetNrows());

  for(int i=0; i<out.GetNrows(); i++){
    for(int j=0; j<out.GetNcols(); j++){
      out(i,j)= M(i,j)*scale;
    }
  }


  return out;

}
