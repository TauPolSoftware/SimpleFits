#include "SimpleFits/FitSoftware/interface/TauA1NuConstrainedFitter.h"
#include "Validation/EventGenerator/interface/PdtPdgMini.h"
#include <iostream>

TauA1NuConstrainedFitter::TauA1NuConstrainedFitter(unsigned int ambiguity,LorentzVectorParticle A1,TVector3 PVertex, TMatrixTSym<double> VertexCov,double mtau):
  LagrangeMultipliersFitter(),
  MultiProngTauSolver(mtau),
  mtau_c(mtau),
  ambiguity_(ambiguity)
{
  TLorentzVector Tau(0,0,0,0);
  LorentzVectorParticle Nu=EstimateNu(A1,PVertex, ambiguity,Tau);
  particles_.push_back(A1);
  particles_.push_back(Nu);

  isconfigured=false;
  // setup 13 by 13 matrix
  int size=LorentzVectorParticle::NVertex+particles_.size()*LorentzVectorParticle::NLorentzandVertexPar;
  TMatrixT<double>    inpar(size,1);
  TMatrixTSym<double> incov(size);

  // Get primary vertex information
  if(VertexCov.GetNrows()!=LorentzVectorParticle::NVertex)return;
  inpar(LorentzVectorParticle::vx,0)=PVertex.X();
  inpar(LorentzVectorParticle::vy,0)=PVertex.Y();
  inpar(LorentzVectorParticle::vz,0)=PVertex.Z();  
  for(int i=0; i<LorentzVectorParticle::NVertex;i++){
    for(int j=0; j<LorentzVectorParticle::NVertex;j++)incov(i,j)=VertexCov(i,j);
  }
  int A1offset=LorentzVectorParticle::NVertex;
  int Nuoffset=LorentzVectorParticle::NLorentzandVertexPar+LorentzVectorParticle::NVertex;
  for(int i=0; i<LorentzVectorParticle::NLorentzandVertexPar;i++){
    inpar(i+A1offset,0)=A1.Parameter(i);
    inpar(i+Nuoffset,0)=Nu.Parameter(i);
    for(int j=0; j<LorentzVectorParticle::NLorentzandVertexPar;j++){
      incov(i+A1offset,j+A1offset)=A1.Covariance(i,j);
      incov(i+Nuoffset,j+Nuoffset)=Nu.Covariance(i,j);
    }
  }
  // store expanded par for computation of final par (assumes fit has neglegible impact on a1 correlations with vertex errors)
  exppar.ResizeTo(nexpandedpar,1);
  exppar=ComputeInitalExpPar(inpar);
  expcov.ResizeTo(nexpandedpar,nexpandedpar);
  expcov=ErrorMatrixPropagator::PropogateError(&TauA1NuConstrainedFitter::ComputeInitalExpPar,inpar,incov);
  // store linearization point
  TMatrixT<double> PAR_0(npar,1);
  par_0.ResizeTo(npar);
  cov_0.ResizeTo(npar,npar);
  PAR_0=ComputeExpParToPar(exppar);
  for(int i=0; i<npar;i++)par_0(i)=PAR_0(i,0);
  cov_0=ErrorMatrixPropagator::PropogateError(&TauA1NuConstrainedFitter::ComputeExpParToPar,exppar,expcov);

  for(int i=0; i<npar;i++){
    for(int j=0;j<npar;j++){cov_0(i,j)=expcov(i,j);}
  }
  // set up inital point for fit (cov handled in Fit() function)
  par.ResizeTo(npar);
  par=par_0;

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Check if Tau Direction is unphysical and if nessicary set the starting point to Theta_{GJ-Max} 
  /*
  TLorentzVector a1(par(a1_px),par(a1_py),par(a1_pz),sqrt(par(a1_m)*par(a1_m)+par(a1_px)*par(a1_px)+par(a1_py)*par(a1_py)+par(a1_pz)*par(a1_pz)));
  double phi(par(tau_phi)),theta(par(tau_theta));
  if(SetTauDirectionatThetaGJMax(a1,theta,phi,1/sqrt(2.0))){
    std::cout <<  "resetting phi and theta" << std::endl;
    TLorentzVector Tau_plus,Tau_minus,nu_plus,nu_minus;
    TVector3 TauDir; TauDir.SetMagThetaPhi(1.0,theta,phi);
    SolvebyRotation(TauDir,a1,Tau_plus,Tau_minus,nu_plus,nu_minus);
    par(tau_phi)=phi;
    par(tau_theta)=theta;
    if(ambiguity==plus){
      nu_plus.Print();
      par(nu_px)=nu_plus.Px();
      par(nu_py)=nu_plus.Py();
      par(nu_pz)=nu_plus.Pz();
    }
    if(ambiguity==minus){
      nu_minus.Print();
      par(nu_px)=nu_minus.Px();
      par(nu_py)=nu_minus.Py();
      par(nu_pz)=nu_minus.Pz();
    }
  }
  */
  isconfigured=true;  
}


TMatrixT<double> TauA1NuConstrainedFitter::ComputeInitalExpPar(TMatrixT<double> &inpar){
  TMatrixT<double> outpar(nexpandedpar,1);
  int offset=LorentzVectorParticle::NVertex;// for A1
  TVector3 pv(inpar(LorentzVectorParticle::vx,0),inpar(LorentzVectorParticle::vy,0),inpar(LorentzVectorParticle::vz,0));
  TVector3 sv(inpar(LorentzVectorParticle::vx+offset,0),inpar(LorentzVectorParticle::vy+offset,0),inpar(LorentzVectorParticle::vz+offset,0));
  TVector3 TauDir=sv-pv;
  outpar(tau_phi,0)=TauDir.Phi();
  outpar(tau_theta,0)=TauDir.Theta();
  outpar(a1_px,0)=inpar(LorentzVectorParticle::px+offset,0);
  outpar(a1_py,0)=inpar(LorentzVectorParticle::py+offset,0);
  outpar(a1_pz,0)=inpar(LorentzVectorParticle::pz+offset,0);
  outpar(a1_m,0)=inpar(LorentzVectorParticle::m+offset,0);
  outpar(a1_vx,0)=inpar(LorentzVectorParticle::vx+offset,0);
  outpar(a1_vy,0)=inpar(LorentzVectorParticle::vy+offset,0);
  outpar(a1_vz,0)=inpar(LorentzVectorParticle::vz+offset,0);
  offset+=LorentzVectorParticle::NLorentzandVertexPar; // for Nu
  outpar(nu_px,0)=inpar(LorentzVectorParticle::px+offset,0);
  outpar(nu_py,0)=inpar(LorentzVectorParticle::py+offset,0);
  outpar(nu_pz,0)=inpar(LorentzVectorParticle::pz+offset,0);
  return outpar;
}


TMatrixT<double> TauA1NuConstrainedFitter::ComputeExpParToPar(TMatrixT<double> &inpar){
  TMatrixT<double> outpar(npar,1);
  for(int i=0;i<npar;i++){outpar(i,0)=inpar(i,0);}
  return outpar;
}


TMatrixT<double> TauA1NuConstrainedFitter::ComputeNuLorentzVectorPar(TMatrixT<double> &inpar){
  TMatrixT<double> outpar(LorentzVectorParticle::NLorentzandVertexPar,1);
  outpar(LorentzVectorParticle::vx,0)=inpar(a1_vx,0);
  outpar(LorentzVectorParticle::vy,0)=inpar(a1_vy,0);
  outpar(LorentzVectorParticle::vz,0)=inpar(a1_vz,0);
  outpar(LorentzVectorParticle::px,0)=inpar(nu_px,0);
  outpar(LorentzVectorParticle::py,0)=inpar(nu_py,0);
  outpar(LorentzVectorParticle::pz,0)=inpar(nu_pz,0);
  outpar(LorentzVectorParticle::m,0)=0;
  return outpar;
}

TMatrixT<double> TauA1NuConstrainedFitter::ComputeA1LorentzVectorPar(TMatrixT<double> &inpar){
  TMatrixT<double> outpar(LorentzVectorParticle::NLorentzandVertexPar,1);
  outpar(LorentzVectorParticle::vx,0)=inpar(a1_vx,0);
  outpar(LorentzVectorParticle::vy,0)=inpar(a1_vy,0);
  outpar(LorentzVectorParticle::vz,0)=inpar(a1_vz,0);
  outpar(LorentzVectorParticle::px,0)=inpar(a1_px,0);
  outpar(LorentzVectorParticle::py,0)=inpar(a1_py,0);
  outpar(LorentzVectorParticle::pz,0)=inpar(a1_pz,0);
  outpar(LorentzVectorParticle::m,0)=inpar(a1_m,0);
  return outpar;
}

TMatrixT<double> TauA1NuConstrainedFitter::ComputeMotherLorentzVectorPar(TMatrixT<double> &inpar){
  TMatrixT<double> outpar(LorentzVectorParticle::NLorentzandVertexPar,1);
  TMatrixT<double> nupar=ComputeNuLorentzVectorPar(inpar);
  TMatrixT<double> a1par=ComputeA1LorentzVectorPar(inpar);
  for(int i=0;i<LorentzVectorParticle::NLorentzandVertexPar;i++){
    if(i==LorentzVectorParticle::m)continue;
    if(i<LorentzVectorParticle::NVertex){outpar(i,0)=a1par(i,0);}
    else{outpar(i,0)=nupar(i,0)+a1par(i,0);}
  }
  double Enu2=pow(nupar(LorentzVectorParticle::px,0),2.0)+pow(nupar(LorentzVectorParticle::py,0),2.0)+pow(nupar(LorentzVectorParticle::pz,0),2.0);
  double Ea12=pow(a1par(LorentzVectorParticle::px,0),2.0)+pow(a1par(LorentzVectorParticle::py,0),2.0)+pow(a1par(LorentzVectorParticle::pz,0),2.0)+pow(a1par(LorentzVectorParticle::m,0),2.0);
  double P2=pow(outpar(LorentzVectorParticle::px,0),2.0)+pow(outpar(LorentzVectorParticle::py,0),2.0)+pow(outpar(LorentzVectorParticle::pz,0),2.0);
  outpar(LorentzVectorParticle::m,0)=sqrt(fabs(pow(sqrt(Enu2)+sqrt(Ea12),2.0)-P2));
  return outpar;
}

void TauA1NuConstrainedFitter::UpdateExpandedPar(){
  // assumes changes to a1 correlation to vertex is small
  if(par.GetNrows()==npar && cov.GetNrows() && exppar.GetNrows()==npar && expcov.GetNrows()) return;
  for(int i=0; i<npar;i++){
    exppar(i,0)=par(i);
    for(int j=0; j<npar;j++){expcov(i,j)=cov(i,j);}
  }
}

std::vector<LorentzVectorParticle> TauA1NuConstrainedFitter::GetReFitDaughters(){
  std::vector<LorentzVectorParticle> refitParticles;
  UpdateExpandedPar();
  double c(0),b(0);
  for(unsigned int i=0;i<particles_.size();i++){c+=particles_.at(i).Charge();b=particles_.at(i).BField();}
  TMatrixT<double> a1=ComputeA1LorentzVectorPar(exppar);
  TMatrixTSym<double> a1cov=ErrorMatrixPropagator::PropogateError(&TauA1NuConstrainedFitter::ComputeA1LorentzVectorPar,exppar,expcov);
  refitParticles.push_back(LorentzVectorParticle(a1,a1cov,fabs(PdtPdgMini::a_1_plus)*c,c,b));
  TMatrixT<double> nu=ComputeNuLorentzVectorPar(exppar);
  TMatrixTSym<double> nucov=ErrorMatrixPropagator::PropogateError(&TauA1NuConstrainedFitter::ComputeNuLorentzVectorPar,exppar,expcov);
  refitParticles.push_back(LorentzVectorParticle(nu,nucov,PdtPdgMini::nu_tau,0.0,b));
  return refitParticles;
}

LorentzVectorParticle TauA1NuConstrainedFitter::GetMother(){
  UpdateExpandedPar();
  double c(0),b(0);
  for(unsigned int i=0;i<particles_.size();i++){c+=particles_.at(i).Charge();b=particles_.at(i).BField();}
  TMatrixT<double> m=ComputeMotherLorentzVectorPar(exppar);
  TMatrixTSym<double> mcov=ErrorMatrixPropagator::PropogateError(&TauA1NuConstrainedFitter::ComputeMotherLorentzVectorPar,exppar,expcov);
  return LorentzVectorParticle(m,mcov,-1.0*fabs(PdtPdgMini::tau_minus)*c,c,b);
}

TVectorD TauA1NuConstrainedFitter::Value(TVectorD &v){
  TLorentzVector a1,nu;
  double phi(0),theta(0);
  TVector3 TauDir;
  CovertParToObjects(v,a1,nu,phi,theta,TauDir);
  TLorentzVector a1_d=a1;
  TLorentzVector nu_d=nu;
  TLorentzVector Tau_plus,Tau_minus,nu_plus,nu_minus;
  SolvebyRotation(TauDir,a1_d,Tau_plus,Tau_minus,nu_plus,nu_minus,false);
  a1.RotateZ(-phi);
  a1.RotateY(-theta);
  nu.RotateZ(-phi);
  nu.RotateY(-theta);
  TLorentzVector nufixed(-a1.Px(),-a1.Py(),nu.Pz(),sqrt(a1.Pt()*a1.Pt()+nu.Pz()*nu.Pz()));
  TLorentzVector tau=a1+nufixed;
  TVectorD d(3);
  if(ambiguity_==minus){    d(0)=sqrt(pow(nu.Pz()-nu_minus.Pz(),4.0)+pow(tau.M2()-mtau_c*mtau_c,2.0));}
  else if(ambiguity_==plus){d(0)=sqrt(pow(nu.Pz()-nu_plus.Pz(),4.0)+pow(tau.M2()-mtau_c*mtau_c,2.0));}
  else {d(0) = tau.M2()-mtau_c*mtau_c;}
  d(1) = a1.Px()+nu.Px();
  d(2) = a1.Py()+nu.Py();
  return d;
}


void TauA1NuConstrainedFitter::CovertParToObjects(TVectorD &v,TLorentzVector &a1,TLorentzVector &nu,double &phi,double &theta,TVector3 &TauDir){
  a1=TLorentzVector(v(a1_px),v(a1_py),v(a1_pz),sqrt(v(a1_m)*v(a1_m)+v(a1_px)*v(a1_px)+v(a1_py)*v(a1_py)+v(a1_pz)*v(a1_pz)));
  nu=TLorentzVector(v(nu_px),v(nu_py),v(nu_pz),sqrt(v(nu_px)*v(nu_px)+v(nu_py)*v(nu_py)+v(nu_pz)*v(nu_pz)));
  phi=par(tau_phi);
  theta=par(tau_theta);
  TauDir.SetMagThetaPhi(1.0,theta,phi);
}


bool TauA1NuConstrainedFitter::Fit(){
  ////////////////////////////////////////////
  // Run Kicker to force +/- solution to avoid solution being stuck in the local minimum
  if(ambiguity_==minus || ambiguity_==plus){ 
    TLorentzVector a1,nu;
    double phi(0),theta(0);
    TVector3 TauDir;
    CovertParToObjects(par,a1,nu,phi,theta,TauDir);
    TLorentzVector Tau_plus,Tau_minus,nu_plus,nu_minus,nu_correct,nu_incorrect;
    if(ambiguity_==minus)SolvebyRotation(TauDir,a1,Tau_plus,Tau_minus,nu_incorrect,nu_correct,false);
    if(ambiguity_==plus)SolvebyRotation(TauDir,a1,Tau_plus,Tau_minus,nu_correct,nu_incorrect,false);
    nu.RotateZ(-phi);
    nu.RotateY(-theta);
    if(fabs(nu_incorrect.Pz()-nu.Pz())<fabs(nu_correct.Pz()-nu.Pz())){
      double pzkicked=nu_correct.Pz()-(nu_incorrect.Pz()-nu.Pz()); // minus sign is to make the kick a reflex about the ambiguity point 
      TLorentzVector nuKicked(nu.Px(),nu.Py(),pzkicked,sqrt(nu.Px()*nu.Px()+nu.Py()*nu.Py()+pzkicked*pzkicked));
      nuKicked.RotateY(-theta);
      nuKicked.RotateZ(-phi);
      par(nu_px)=nuKicked.Px();
      par(nu_py)=nuKicked.Py();
      par(nu_pz)=nuKicked.Pz();
    }
  }
  return LagrangeMultipliersFitter::Fit();
}

