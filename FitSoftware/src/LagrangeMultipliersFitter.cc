// Code written by Ian M. Nugent
// minor modified by Vladimir Cherepanov on RWTH Aachen March 4 2014 (update 2)

#include "SimpleFits/FitSoftware/interface/LagrangeMultipliersFitter.h"
#include "TDecompBK.h"
#include <iostream>

LagrangeMultipliersFitter::LagrangeMultipliersFitter():
  isconfigured(false), 
  isFit(false),
  epsilon_(0.001),
  weight_(1.0),
  MaxDelta_(0.01),
  nitermax_(100),
  chi2(1e10),
  D(1,1), 
  V_D(1,1)
{

}

bool LagrangeMultipliersFitter::Fit(){
  //  std::cout<<"LagrangeMultipliersFitter  deb 1  "<<std::endl;
  if(cov.GetNrows()!=par_0.GetNrows()){
    //std::cout<<"LagrangeMultipliersFitter  deb 2  "<<std::endl;
    // set cov to cov_0 until value is computed
    cov.ResizeTo(par_0.GetNrows(),par_0.GetNrows());
    cov=cov_0;
 
  }
  if(!isconfigured) return false;
  if(isFit)return isConverged();
  isFit=true;
  niter=0;
  for(niter=0;niter<=nitermax_;niter++){
    //std::cout<<"before ApplyLagrangianConstraints:   chi2 delta  "<<chi2<<"   "<< delta<<"   probability  " <<  TMath::Prob(chi2,2)<<std::endl;
    bool passed=ApplyLagrangianConstraints();
    //std::cout<<"fit   chi2 delta  "<<chi2<<"   "<< delta<<"   probability  " <<  TMath::Prob(chi2,2)<<std::endl;
    if (!passed || (niter==nitermax_ && delta>=4.0*MaxDelta_)) {
      //std::cout << "Reached Maximum number of iterations..." << niter << " and delta "<< delta <<std::endl;

      return false;
      //return true;
    }
    if(isConverged()) break; 
  }
    //std::cout<<" Lagrangian constrained fitter =====> "<< chi2 <<std::endl;
  ComputeVariance();
//   std::cout<<"info before return true  fit   chi2 delta  "<<chi2<<"   "<< delta<<"   probability  " <<  TMath::Prob(chi2,2)<<std::endl;
//   std::cout<<"LagrangeMultipliersFitter::Fit  2 "<<  par(0) << "   " <<par(1) <<"   " << par(2) <<std::endl;
  return true;
}



bool LagrangeMultipliersFitter::ApplyLagrangianConstraints(){
  if(V_D.GetNrows()!=NConstraints()) V_D.ResizeTo(NConstraints(),NConstraints());
  if(D.GetNrows()!=NConstraints() || D.GetNcols()!=par.GetNrows()) D.ResizeTo(NConstraints(),par.GetNrows());

  //std::cout << "ApplyLagrangianConstraints debug: 1" <<std::endl;
  // Setup intial values
  TMatrixT<double> alpha_A=convertToMatrix(par);
  TMatrixT<double> alpha_0=convertToMatrix(par_0);
  TMatrixT<double> delta_alpha_A=alpha_A-alpha_0;

  //std::cout << "ApplyLagrangianConstraints debug: 2" <<std::endl;
  D=Derivative();

  //std::cout << "Print Derivatives matrix" <<std::endl;
  //D.Print();
  //std::cout << "ApplyLagrangianConstraints debug: 3" <<std::endl;

  // std::cout<<"call value to compute Value   "<<std::endl;
  TMatrixT<double> d=convertToMatrix(Value(par));
  TMatrixT<double> C=D*delta_alpha_A-d;
  TMatrixTSym<double> V_alpha0=cov_0;
  TMatrixTSym<double> V_D_inv=V_alpha0;

  //if(niter==0){
  //	std::cout << "VD matrix with initial cov" << std::endl;
  //	V_D_inv.Print();
  //}

  V_D_inv.Similarity(D);

  //if(niter==0){
  //	std::cout << "VD matrix after similarity" << std::endl;
  //	V_D_inv.Print();
  //}

  double det = V_D_inv.Determinant();

  //TMatrixTSym<double> X(V_D_inv);
  //std::cout << "PRINT X!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
  //std::cout << "X Det" << X.Determinant() <<  std::endl;
  //X.Print();
  //X.Invert();
  //X.Print();
  //if(niter==0){
  //	std::cout << "VD matrix det: " << det <<std::endl;
  //}

  TDecompBK Inverter(V_D_inv);
  if(fabs(det)>1e40){
       std::cout << "Fit failed: unable to invert SYM gain matrix LARGE Determinant" << det << " \n" << std::endl;
    return false;
  }
  if(!Inverter.Decompose()){
        std::cout << "Fit failed: unable to invert SYM gain matrix " << det << " \n" << std::endl;
        //std::cout << "Inverter.Decompose()" << Inverter.Decompose() << std::endl;
        //V_D_inv.Invert();
        //cov_0.Print();
        //d.Print();
        //D.Print();
        //V_D_inv.Print();
        //std::cout << "V_D_inv.Determinant()" << V_D_inv.Determinant() << std::endl;
    return false;
  }
  V_D=Inverter.Invert();
  
  //if(niter==0){
  //	std::cout << "VD matrix after inversion" << std::endl;
  //	V_D.Print();
  //}

  // solve equations
  TMatrixT<double> lambda=-1.0*V_D*C;
  TMatrixT<double> DT=D; DT.T();
  TMatrixT<double> alpha=alpha_0-V_alpha0*DT*lambda;


  TVectorD finPar=convertToVector(alpha);

  // do while loop to see if the convergence criteria are satisfied
  double s(1), stepscale(0.01);
  chi2prev=chi2;
  //std::cout<<"chi2 "<<chi2<< std::endl;
  double Curentchi2(ChiSquareUsingInitalPoint(alpha_A,lambda)), Currentdelta(ConstraintDelta(par));
  //std::cout<<"Curentchi2 "<<Curentchi2<< std::endl;
  TMatrixT<double> alpha_s=alpha;

  // convergence in 2 step procedure to minimize chi2 within MaxDelta_ of the constraints
  // 1) Get within 5x MaxDelta_
  // 2) converge based on improving chi2 and constrained delta
  unsigned int Proc=ConstraintMin;
  if(ConstraintDelta(par)<5*MaxDelta_)Proc=Chi2AndConstaintMin;
  int  NIter=(int)(1.0/stepscale);
  for(int iter=0;iter<NIter;iter++){
    // compute safety cutoff for numerical constraint
    double diff=0;
    for(int l=0;l<alpha_s.GetNrows();l++){
      if(diff<alpha_s(l,0)-alpha_A(l,0))diff=alpha_s(l,0)-alpha_A(l,0);
    }
    double delta_alpha_s=ConstraintDelta(convertToVector(alpha_s));
    if(Proc==ConstraintMin){
      if(delta_alpha_s<Currentdelta || iter==NIter || diff<100*epsilon_){Curentchi2=ChiSquareUsingInitalPoint(alpha_s,lambda); Currentdelta=delta_alpha_s; ScaleFactor=s; break;}
    }
    else if(Proc==Chi2AndConstaintMin){
      double chi2_s=ChiSquareUsingInitalPoint(alpha_s,lambda);
      //std::cout<<"chi2_s "<<chi2_s<< std::endl;
      if((delta_alpha_s<Currentdelta/*+MaxDelta_*/ && chi2_s<Curentchi2) || iter==NIter || diff<100*epsilon_){
    	Curentchi2=chi2_s;
    	Currentdelta=delta_alpha_s;
    	ScaleFactor=s;
    	//if(iter==NIter) std::cout << "Proc==Chi2AndConstraintMin broke, because: iter==NIter, iter: " << iter << std::endl;
    	//if(diff<1*epsilon_) std::cout << "Proc==Chi2AndConstraintMin broke, because: diff<1*epsilon_, diff: " << diff << std::endl;
    	break;}
    }
    s-=stepscale;
    alpha_s=alpha_A+s*(alpha-alpha_A);
  }


  //std::cout<<"Curentchi2 "<<Curentchi2<< std::endl;

  // set chi2
  chi2=Curentchi2;  
  //std::cout<<"chi2 "<<chi2<< std::endl;

  //set delta
  delta=Currentdelta;
  //   std::cout << "LagrangeMultipliersFitter Chi^2 " << chi2 << " delta " << Currentdelta << std::endl; 
  //correct finPar to new stepsize
  par=convertToVector(alpha_s);


  for(int i=0; i<alpha_s.GetNrows();i++){

    //   std::cout<<"correct finParfinPar   "<<alpha_s(i,0)<<std::endl;
  }
  return true;
}
TMatrixD LagrangeMultipliersFitter::Derivative(){ // alway evaluated at current par
  //std::cout<<" deb 1"<<std::endl;
  TMatrixD Derivatives(NConstraints(),par.GetNrows());
  TVectorD par_plus(par.GetNrows());
  TVectorD value(NConstraints());
  TVectorD value_plus(NConstraints());
  //std::cout<<" deb 2"<<std::endl;
  for(int j=0;j<par.GetNrows();j++){
    for(int i=0;i<par.GetNrows();i++){
      par_plus(i)=par(i);
            //std::cout<<" deb 3 "<<epsilon_ <<std::endl;
      if(i==j) par_plus(i)=par(i)+epsilon_;
    }
   //std::cout<<" deb 4"<<std::endl;
//   std::cout<<" paramters   "<<par(0) <<" "<<par(1) <<" "<<par(2) <<" "<<par(3) <<" "<<par(4) <<" "<<par(5) <<" "<<std::endl;
//   std::cout<<" paramters_plus   "<<par_plus(0) <<" "<<par_plus(1) <<" "<<par_plus(2) <<" "<<par_plus(3) <<" "<<par_plus(4) <<" "<<par_plus(5) <<" "<<std::endl;

    value=Value(par);
      //std::cout<<" deb 5"<<std::endl;
    value_plus=Value(par_plus);
    for(int i=0; i<NConstraints();i++){
        //std::cout<<" deb 6"<<std::endl;
      Derivatives(i,j)=(value_plus(i)-value(i))/epsilon_;
    }
  }
  return Derivatives;
}


bool LagrangeMultipliersFitter::isConverged(){
  if(delta<MaxDelta_ /*&& chi2prev-chi2<1.0 && chi2prev>chi2*/){
    //    std::cout << "converged " << delta << " chi2 " <<  chi2 << " chi2prev " << chi2prev <<"  Maxdelta  " <<MaxDelta_ <<std::endl; 
    return true;
  }
  return false;
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



double LagrangeMultipliersFitter::ChiSquare(TMatrixT<double> delta_alpha,TMatrixT<double> lambda,TMatrixT<double> D, TMatrixT<double> d){
  TMatrixT<double> lambdaT=lambda; lambdaT.T();
  TMatrixT<double> chisquare=lambdaT*(D*delta_alpha+d);

  std::cout << "chi2  " <<chisquare(0,0) <<std::endl;
  double c2=chisquare(0,0);
  return c2;
}

double LagrangeMultipliersFitter::ChiSquareUsingInitalPoint(TMatrixT<double> alpha,TMatrixT<double> lambda){
  if(cov_0.GetNrows()!=V_alpha0_inv.GetNrows()){
  	//std::cout << "cov_0 matrix" << std::endl;
  	//cov_0.Print();
    TMatrixTSym<double> V_alpha0=cov_0;
    V_alpha0_inv.ResizeTo(cov_0.GetNrows(),cov_0.GetNrows());
  	//std::cout << "V_alpha0_inv. matrix" << std::endl;
  	//V_alpha0_inv.Print();
    TDecompBK Inverter(V_alpha0);
    if(!Inverter.Decompose()){ // handle rare case where inversion is not possible (ie assume diagonal)
      std::cout << "LagrangeMultipliersFitter::ChiSquareUsingInitalPoint: Error non-invertable Matrix... Calculating under assumption that correlations can be neglected!!!" << std::endl;
      for(int j=0;j<par.GetNrows();j++){
	for(int i=0;i<par.GetNrows();i++){
	  if(i==j) V_alpha0_inv(i,j)=1.0/V_alpha0(i,j);
	  else V_alpha0_inv(i,j)=0.0;
	}
      }
    }
    else{
      V_alpha0_inv=Inverter.Invert();
    }
  }
	//std::cout << "V_alpha0_inv. matrix" << std::endl;
	//V_alpha0_inv.Print();

  TMatrixT<double> lambdaT=lambda; lambdaT.T();
  TMatrixT<double> alpha_0=convertToMatrix(par_0);
  TMatrixT<double> dalpha=alpha-alpha_0;
  TMatrixT<double> dalphaT=dalpha;  dalphaT.T();
  TMatrixT<double> chisquare_var=dalphaT*(V_alpha0_inv*dalpha);
  TVectorT<double> alpha_v=convertToVector(alpha);
  //  std::cout<<"call value to compute chi2   "<<std::endl;
  TMatrixT<double> chisquare_constraints=lambdaT*convertToMatrix(Value(alpha_v));
  double c2=chisquare_var(0,0)+chisquare_constraints(0,0);
  return c2;
}

double LagrangeMultipliersFitter::ConstraintDelta(TVectorT<double> par){
  //  std::cout<<"call value to compute ConstraintDelta   "<<std::endl;
  TVectorD d_par=Value(par);
  double delta_d(0);
  for(int i = 0; i<d_par.GetNrows(); i++){
    delta_d+=fabs(d_par(i));
  }
  return delta_d;
}


TMatrixT<double>  LagrangeMultipliersFitter::ComputeVariance(){
  TMatrixTSym<double> V_alpha0=cov_0;
  TMatrixTSym<double> DTV_DD=V_D.SimilarityT(D);
  TMatrixT<double> DTV_DDV=DTV_DD*V_alpha0;
  TMatrixT<double> VDTV_DDV=V_alpha0*DTV_DDV;
  TMatrixT<double> CovCor=VDTV_DDV;
  //CovCor*=ScaleFactor;
  if(V_corr_prev.GetNrows()!=V_alpha0.GetNrows()){
    V_corr_prev.ResizeTo(V_alpha0.GetNrows(),V_alpha0.GetNrows());
    V_corr_prev=CovCor;
  }
  else{
    V_corr_prev*=(1-ScaleFactor);
    CovCor+=V_corr_prev;
    V_corr_prev=CovCor;
  }
  
  TMatrixT<double> V_alpha = V_alpha0-CovCor;
  for(int i=0; i<cov.GetNrows();i++){
    for(int j=0; j<=i;j++){
      cov(i,j)=V_alpha(i,j);
    }
  }
  return cov;
}
