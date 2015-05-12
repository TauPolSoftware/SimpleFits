#include "SimpleFits/FitSoftware/interface/LagrangeMultipliersFitter.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"
#include "TDecompBK.h"
#include <iostream>

LagrangeMultipliersFitter::LagrangeMultipliersFitter():
  isconfigured(false),
  isFit(false),
  epsilon_(0.001),
  weight_(1.0),
  MaxDelta_(0.1),
  nitermax_(50),
  chi2(1e10),
  MaxParDelta_(0.1),
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
	  Logger(Logger::Debug) << "Debug 8.8" << std::endl;

    bool passed=ApplyLagrangianConstraints();
    
    //std::cout<<"global iteration ========== > "<< niter<<" fit   chi2 delta  "<<chi2<<"   "<< delta<<"   probability  " <<  TMath::Prob(chi2,1)<<std::endl; 
    if (!passed || (niter==nitermax_ && delta>=4.0*MaxDelta_)) {
      std::cout << "Reached Maximum number of iterations..." << niter << " and delta "<< delta <<std::endl;
       return false;
      //  return true;
    }
    if(isConverged()) break; 
  
  }

   ComputeVariancea();
   ComputeVarianceb();

  return true;
}



bool LagrangeMultipliersFitter::ApplyLagrangianConstraints(){
  //  if(V_a.GetNrows()!=para.GetNrows()) V_a.ResizeTo(para.GetNrows(),para.GetNrows());
  if(A.GetNrows()!=NConstraints() || A.GetNcols()!=para.GetNrows()) A.ResizeTo(NConstraints(),para.GetNrows());
  if(B.GetNrows()!=NConstraints() || B.GetNcols()!=parb.GetNrows()) B.ResizeTo(NConstraints(),parb.GetNrows());
  if(Fa.GetNrows()!=NSoftConstraints() || Fa.GetNcols()!=para.GetNrows()) Fa.ResizeTo(NSoftConstraints(),para.GetNrows());
  if(Fb.GetNrows()!=NSoftConstraints() || Fb.GetNcols()!=parb.GetNrows()) Fb.ResizeTo(NSoftConstraints(),parb.GetNrows());



  // Setup intial values

  // TMatrixT<double> alpha_A=convertToMatrix(par);
  // TMatrixT<double> alpha_0=convertToMatrix(par_0);

  // TMatrixT<double> delta_alpha_A=alpha_A-alpha_0;
  // Setup initial values II
 
  TMatrixT<double> y  = convertToMatrix(para_0);
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
  TMatrixTSym<double> V_f=V_a;//ComputeV_f(V_a,para_0, parb_0);
  V_f.SetTol(1.e-40);
  V_f.Similarity(Fa);

  //----  fill final matrix blocks 

  TMatrixTSym<double> V_a_inv= V_a; V_a_inv.Invert();



 

  TMatrixTSym<double> V_f_inv= V_f; V_f_inv.Invert();
  TMatrixT<double> M11 = V_a_inv + FaT*V_f_inv*Fa;
  TMatrixT<double> M12 = FaT*V_f_inv*Fb;
  TMatrixT<double> M21 = FbT*V_f_inv*Fa;
  TMatrixT<double> M22 = FbT*V_f_inv*Fb;

  TMatrixT<double> V1 = V_a_inv*y - FaT*V_f_inv*(f0 - Fa*a0 - Fb*b0);
  TMatrixT<double> V2 =  FbT*V_f_inv*(Fa*a0 + Fb*b0 - f0);
  TMatrixT<double> V3 = A*a0 + B*b0  - c0;
  TMatrixT<double> M = MakeFullMatrix(M11,M12,M21,M22,A,B);
  TMatrixT<double> V = MakeFullVector(V1,V2,V3);

  // V_f_inv.Print();



  double detM = M.Determinant();
  if(fabs(detM)>1e40){
       Logger(Logger::Error) << "Fit failed: unable to invert SYM gain matrix LARGE Determinant" << detM << " \n" << std::endl;
       return false;
  }
  TMatrixT<double> M_inv = M; M_inv.Invert();

 // solve equations
  TMatrixT<double> res = M_inv*V;
 

  TMatrixT<double> par_a = solutiona(res);
  TMatrixT<double> par_b = solutionb(res);
  TMatrixT<double> lambda=solutionlambda(res);


  para  = convertToVector(par_a);
  parb  = convertToVector(par_b);

  // std::cout<<" par_a  " <<std::endl;
  // par_a.Print();

  // std::cout<<" par_b  " <<std::endl;
  // par_b.Print();

  // do while loop to see if the convergance criteria are satisfied
  double s(1), stepscale(0.05);
  chi2prev=chi2;



  
  double Curentchi2(ChiSquareUsingInitalPoint(y,par_a,par_b,lambda)), Currentdelta(ConstraintDelta(para,parb));

  TMatrixT<double> a_s=par_a;
  TMatrixT<double> b_s=par_b;

  // convergence in 2 step procedure to minimize chi2 within MaxDelta_ of the constriants
  // 1) Get within 5x MaxDelta_
  // 2) converge based on improving chi2 and constrianed delta
  unsigned int Proc=ConstraintMin;
  if(ConstraintDelta(para,parb)<5*MaxDelta_)Proc=Chi2AndConstaintMin;
  int  NIter=(int)(1.0/stepscale);



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
  //set delta
  delta=Currentdelta;
  //    std::cout << "LagrangeMultipliersFitter Chi^2 " << chi2 << " delta " << Currentdelta << std::endl; 
  //correct finPar to new stepsize
  para = convertToVector(a_s);
  parb = convertToVector(b_s);
  for(int l=0;l<par_a.GetNrows();l++){
    pardelta=y(l,0)-a_s(l,0) +  b_s(l,0) - b0(l,0);
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

  paraFD(0) = para_0(0);
  paraFD(1) = para_0(1);
  paraFD(2) = atan(para_0(1)/para_0(0));

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
  Derivatives(2,2) =1;
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
   Derivatives(2,2) =1;
  return Derivatives;
}


bool LagrangeMultipliersFitter::isConverged(){
  if(pardelta<MaxParDelta_ /*&& chi2prev-chi2<1.0 && chi2prev>chi2*/){
	Logger(Logger::Verbose) << "converged " << delta << " chi2 " <<  chi2 << " chi2prev " << chi2prev <<"  Maxdelta  " <<MaxDelta_ <<std::endl;
    return true;
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
double LagrangeMultipliersFitter::ChiSquareUsingInitalPoint(TMatrixT<double> y, TMatrixT<double> a,TMatrixT<double> b,TMatrixT<double> lambda){
  // if(cova_0.GetNrows()!=V_alpha0_inv.GetNrows()){
  TMatrixTSym<double> V_alpha0=cova_0;
  V_alpha0_inv.ResizeTo(cova_0.GetNrows(),cova_0.GetNrows());
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
  
  TMatrixT<double> lambdaT=lambda; lambdaT.T();
  TMatrixT<double> a0=convertToMatrix(para_0);
  TMatrixT<double> b0=convertToMatrix(parb_0);
  TMatrixT<double> da=y-a;
  

  TMatrixT<double> daT=da;  daT.T();
  TMatrixT<double> chisquare_var=daT*(V_alpha0_inv*da);
  TVectorT<double> a_v=convertToVector(a);
  TVectorT<double> b_v=convertToVector(b);
  TMatrixT<double> f = convertToMatrix(SoftValue(a_v,b_v));
  TMatrixT<double> fT = f; fT.T();
  TMatrixT<double> Fa =DerivativeSCa();
  TMatrixT<double> FaT=Fa; Fa.T();

  //  
  TMatrixT<double> chisquare_constraints=lambdaT*convertToMatrix(HardValue(a_v,b_v));// + fT*(Fa*V_alpha0_inv*FaT)*f;
  double c2=chisquare_var(0,0)+chisquare_constraints(0,0);
  return c2;

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
  return delta_d;
}

TMatrixT<double>  LagrangeMultipliersFitter::ComputeVarianceb(){
 
  TMatrixTSym<double> V_0=cova_0;
  TMatrixTSym<double> DTV_DD=V_0.SimilarityT(Fa);DTV_DD.Invert();
  TMatrixT<double> FbT = Fb; FbT.T();
  TMatrixT<double> CovCor=FbT*DTV_DD*Fb;CovCor.Invert();


  for(int i=0; i<covb.GetNrows();i++){
    for(int j=0; j<=i;j++){
      covb(i,j)=CovCor(i,j);
    }
  }

  return covb;
}

TMatrixT<double>  LagrangeMultipliersFitter::ComputeVariancea(){


 
  TMatrixTSym<double> V_0=cova_0;
  TMatrixTSym<double> DTV_DD=V_0.SimilarityT(Fa);DTV_DD.Invert();
  TMatrixT<double> FbT = Fb; FbT.T();
  TMatrixT<double> Vz=FbT*DTV_DD*Fb;
  TMatrixT<double> Vlambda = DTV_DD - DTV_DD*Vz*DTV_DD;
  TMatrixT<double> Valpha = V_0 - V_0*(V_0.SimilarityT(Fa).T())*V_0;

  
  for(int i=0; i<cova.GetNrows();i++){
    for(int j=0; j<=i;j++){
      cova(i,j)=Valpha(i,j);
    }
  }

  return cova;
}
void  LagrangeMultipliersFitter::Print(TMatrixT<double> M){


   for(int str =0; str < M.GetNrows(); str++){
     for(int kol =0; kol < M.GetNcols(); kol++){
       std::cout<<"  "<< M(str,kol)<<"  ";

     }    
     std::cout<<std::endl;
   }


}
TMatrixTSym<double> LagrangeMultipliersFitter::ComputeV_f(TMatrixTSym<double>  cov, TVectorD para,TVectorD parb){
  TMatrixTSym<double> Vf;
  Vf.ResizeTo(cov.GetNrows(), cov.GetNcols());
  Vf = cov;

  double dgdx = -(para(1) + parb(1))/(pow(para(0) + parb(0),2) + pow(para(1) + parb(1),2));
  double dgdy = (para(0) + parb(0))/(pow(para(0) + parb(0),2) + pow(para(1) + parb(1),2));
  
  double deltaxx = cov(0,0);
  double deltayy = cov(1,1);
  double deltaxy = cov(0,1);

  Vf(2,0) = dgdx*deltaxx + dgdy*deltaxy; Vf(0,2) = Vf(2,0);
  Vf(2,1) = dgdx*deltaxy + dgdy*deltayy; Vf(1,2) = Vf(2,1);
  Vf(2,2) = pow(dgdx,2)*deltaxx + pow(dgdy,2)*deltayy + 2*dgdx*dgdy*deltaxy; 


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
