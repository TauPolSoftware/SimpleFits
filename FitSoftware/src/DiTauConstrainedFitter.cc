#include "SimpleFits/FitSoftware/interface/DiTauConstrainedFitter.h"
#include "SimpleFits/FitSoftware/interface/PDGInfo.h"
#include <iostream>

double DiTauConstrainedFitter::MassConstraint_ = 91.5;

DiTauConstrainedFitter::DiTauConstrainedFitter(LorentzVectorParticle TauA1,TrackParticle MuTrack, TVector3 PVertex, TMatrixTSym<double> VertexCov):
  LagrangeMultipliersFitter()
{

  LorentzVectorParticle TauMuGuess  = EstimateTauMu( PVertex, VertexCov, TauA1.Vertex(),TauA1.VertexCov(), MuTrack,  TauA1);

  particles_.push_back(TauA1);
  particles_.push_back(TauMuGuess);
  particles0_.push_back(TauA1);
  particles0_.push_back(TauMuGuess);

  isconfigured=false;

  // setup 6 by 6 matrix
  int size=particles_.size()*3;
  TMatrixT<double>    inpar(size,1);
  TMatrixTSym<double> incov(size);
 
  // Get primary vertex information
  if(VertexCov.GetNrows()!=LorentzVectorParticle::NVertex)return;

  // set input parameters:  TauA1 - TauMu
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

 // store expanded par for computation of final par (assumes fit has negligible impact on a1 correlations with vertex uncertainties)

  exppar.ResizeTo(size,1);
  exppar=ComputeInitalExpPar(inpar);
  expcov.ResizeTo(size,size);
  expcov=ErrorMatrixPropagator::PropagateError(&DiTauConstrainedFitter::ComputeInitalExpPar,inpar,incov);

  // store linearization point
  TMatrixT<double> PAR_0(size,1);
  par_0.ResizeTo(size);
  cov_0.ResizeTo(size,size);
  PAR_0=ComputeExpParToPar(exppar);
  for(int i=0; i<npar;i++)par_0(i)=PAR_0(i,0);
  cov_0=ErrorMatrixPropagator::PropagateError(&DiTauConstrainedFitter::ComputeExpParToPar,exppar,expcov);

  for(int i=0; i<npar;i++){
	  for(int j=0;j<npar;j++){cov_0(i,j)=expcov(i,j);}
  }
  // set up initial point for fit (Covariance Matrix handled in Fit() function)
  par.ResizeTo(npar);
  par=par_0;

  isconfigured=true; 
}

TMatrixT<double> DiTauConstrainedFitter::ComputeInitalExpPar(TMatrixT<double> &inpar){
  // std::cout<<"deb ComputeInitalExpPar "<<std::endl;
  TMatrixT<double> outpar(6,1);
  int offset=0;//LorentzVectorParticle::NVertex;// for TauA1
  outpar(taua1_px,0)=inpar(LorentzVectorParticle::vx+offset,0);
  outpar(taua1_py,0)=inpar(LorentzVectorParticle::vy+offset,0);
  outpar(taua1_pz,0)=inpar(LorentzVectorParticle::vz+offset,0);

  offset+=3;//LorentzVectorParticle::NLorentzandVertexPar; // for TauMu

  outpar(taumu_px,0)=inpar(LorentzVectorParticle::vx+offset,0);
  outpar(taumu_py,0)=inpar(LorentzVectorParticle::vy+offset,0);
  outpar(taumu_pz,0)=inpar(LorentzVectorParticle::vz+offset,0);

  return outpar; 
}

TMatrixT<double> DiTauConstrainedFitter::ComputeExpParToPar(TMatrixT<double> &inpar){
  TMatrixT<double> outpar(npar,1);
  for(int i=0;i<npar;i++){outpar(i,0)=inpar(i,0);}
  return outpar;
}

TMatrixT<double> DiTauConstrainedFitter::ComputeTauMuLorentzVectorPar(TMatrixT<double> &inpar){
  //  start with index 3 to fill only momenta part
  TMatrixT<double> outpar(7,1);

   outpar(3,0)=inpar(3,0);
   outpar(4,0)=inpar(4,0);
   outpar(5,0)=inpar(5,0);
   outpar(6,0)=1.777;
  return outpar;
}

TMatrixT<double> DiTauConstrainedFitter::ComputeTauA1LorentzVectorPar(TMatrixT<double> &inpar){
  //  start with index 3 to fill only momenta part
  TMatrixT<double> outpar(7,1);
  outpar(3,0)=inpar(0,0);
  outpar(4,0)=inpar(1,0);
  outpar(5,0)=inpar(2,0);
  outpar(6,0) =1.777;
  return outpar;
}

TMatrixT<double> DiTauConstrainedFitter::ComputeMotherLorentzVectorPar(TMatrixT<double> &inpar){
  TMatrixT<double> outpar(7,1);
  TMatrixT<double> Taunupar=ComputeTauMuLorentzVectorPar(inpar);
  TMatrixT<double> Taua1par=ComputeTauA1LorentzVectorPar(inpar);

  outpar(3,0)=Taunupar(3,0)+Taua1par(3,0);
  outpar(4,0)=Taunupar(4,0)+Taua1par(4,0);
  outpar(5,0)=Taunupar(5,0)+Taua1par(5,0);

  double Etaumu2=pow(Taunupar(3,0),2.0)+pow(Taunupar(4,0),2.0)+pow(Taunupar(5,0),2.0)+pow(Taunupar(6,0),2.0);
  double Etaua12=pow(Taua1par(3,0),2.0)+pow(Taua1par(4,0),2.0)+pow(Taua1par(5,0),2.0)+pow(Taua1par(6,0),2.0);
  double P2=pow(outpar(3,0),2.0)+pow(outpar(4,0),2.0)+pow(outpar(5,0),2.0);
  outpar(6,0)=sqrt(fabs(pow(sqrt(Etaumu2)+sqrt(Etaua12),2.0)-P2));
  return outpar;
}

void DiTauConstrainedFitter::UpdateExpandedPar(){
  // assumes changes to a1 correlation to vertex is small
  //if(par.GetNrows()==npar && cov.GetNrows() && exppar.GetNrows()==npar && expcov.GetNrows()) return;
  for(int i=0; i<npar;i++){
	  exppar(i,0)=par(i);
	  for(int j=0; j<npar;j++){expcov(i,j)=cov(i,j);}
  }
}

std::vector<LorentzVectorParticle> DiTauConstrainedFitter::GetReFitDaughters(){
  std::vector<LorentzVectorParticle> refitParticles;
  UpdateExpandedPar();
  double c(0),b(0);
  for(unsigned int i=0;i<particles_.size();i++){c+=particles_.at(i).Charge();b=particles_.at(i).BField();}
  TMatrixT<double> TauA1=ComputeTauA1LorentzVectorPar(exppar);
  TMatrixTSym<double> TauA1Cov=ErrorMatrixPropagator::PropagateError(&DiTauConstrainedFitter::ComputeTauA1LorentzVectorPar,exppar,expcov);
  refitParticles.push_back(LorentzVectorParticle(TauA1,TauA1Cov,fabs(PDGInfo::tau_plus)*c,c,b));

  TMatrixT<double> TauMu=ComputeTauMuLorentzVectorPar(exppar);
  TauMu(0,0)= particles_.at(1).Parameter(LorentzVectorParticle::vx);
  TauMu(1,0)= particles_.at(1).Parameter(LorentzVectorParticle::vy);
  TauMu(2,0)= particles_.at(1).Parameter(LorentzVectorParticle::vz);
  TMatrixTSym<double> TauMuCov=ErrorMatrixPropagator::PropagateError(&DiTauConstrainedFitter::ComputeTauMuLorentzVectorPar,exppar,expcov);
  for(int i=0; i<LorentzVectorParticle::NVertex; i++){
	  for(int j=0; j<LorentzVectorParticle::NVertex; j++){
		  TauMuCov(i,j)=particles_.at(1).VertexCov()(i,j);
	  }
  }
  refitParticles.push_back(LorentzVectorParticle(TauMu,TauMuCov,PDGInfo::tau_minus,0.0,b));
  return refitParticles; 
}

LorentzVectorParticle DiTauConstrainedFitter::GetMother(){
  UpdateExpandedPar();
  double c(0),b(0);
  for(unsigned int i=0;i<particles_.size();i++){c+=particles_.at(i).Charge();b=particles_.at(i).BField();}
  TMatrixT<double> Mother=ComputeMotherLorentzVectorPar(exppar);
  TMatrixTSym<double> MotherCov=ErrorMatrixPropagator::PropagateError(&DiTauConstrainedFitter::ComputeMotherLorentzVectorPar,exppar,expcov);
  return LorentzVectorParticle(Mother,MotherCov,-1.0*fabs(PDGInfo::tau_minus)*c,c,b);
}

TVectorD DiTauConstrainedFitter::Value(TVectorD &v){
  TLorentzVector Taua1,Taumu;
  CovertParToObjects(v,Taua1,Taumu);

  TLorentzVector z=Taua1+Taumu;
  TVectorD d(3);

  d(0) = z.M() - MassConstraint_;
  d(1) = Taua1.Px() + Taumu.Px();
  d(2) = Taua1.Py() + Taumu.Py();

  return d;
}

void DiTauConstrainedFitter::CovertParToObjects(TVectorD &v,TLorentzVector &TauA1,TLorentzVector &TauMu){
  TauA1=TLorentzVector(v(taua1_px),v(taua1_py),v(taua1_pz),sqrt(1.777*1.777+v(taua1_px)*v(taua1_px)+v(taua1_py)*v(taua1_py)+v(taua1_pz)*v(taua1_pz)));
  TauMu=TLorentzVector(v(taumu_px),v(taumu_py),v(taumu_pz),sqrt(1.777*1.777+v(taumu_px)*v(taumu_px)+v(taumu_py)*v(taumu_py)+v(taumu_pz)*v(taumu_pz)));
}

 bool DiTauConstrainedFitter::Fit(){
//   std::cout<<" DiTauConstrainedFitter::Fit  "<<   std::endl;
   return LagrangeMultipliersFitter::Fit();
 }

LorentzVectorParticle
DiTauConstrainedFitter::EstimateTauMu(TVector3 PV,TMatrixTSym<double>  PVCov, TVector3 SV, TMatrixTSym<double>  SVCov,TrackParticle MuTrack, LorentzVectorParticle TauA1){

  TVector3 Point;
  // Define Tau direction
  TVector3 TauDir =  SV + PV;
  TVector3 TauDirError(sqrt(SVCov(0,0) + PVCov(0,0)),
			sqrt(SVCov(1,1) + PVCov(1,1)),
			sqrt(SVCov(2,2) + PVCov(2,2)));

  // Vector parameters
  double dxy   =MuTrack.Parameter(TrackParticle::dxy);
  double kappa =MuTrack.Parameter(TrackParticle::kappa);
  double phi0  =MuTrack.Parameter(TrackParticle::phi);
  double lam   =MuTrack.Parameter(TrackParticle::lambda);
  double dz    =MuTrack.Parameter(TrackParticle::dz);
  double c     =MuTrack.Charge();

  TVector3 ReferencePoint(c*dxy*sin(phi0) +PV.X() ,-c*dxy*cos(phi0)-PV.Y(),dz+PV.Z());
  TVector3 A1SV = -SV + PV;
  //------------------------------
  //all WRT PV

  //=====
  double phiAnot  = atan2(A1SV.Y(), A1SV.X());
  double xpoca2Anot = dxy*sin(phi0) - PV.X();
  double ypoca2Anot = -dxy*cos(phi0) - PV.Y();
  double aAnot = tan(phi0);
  double bAnot = ypoca2Anot - aAnot*xpoca2Anot;
  double r = sqrt( pow(bAnot/(tan(phiAnot) - aAnot )  ,2) + pow(bAnot*tan(phiAnot)/(tan(phiAnot) - aAnot) ,2));//(bAnot)/(tan(phiAnot) - tan(phi0))/cos(phiAnot);

  double bz0 = fabs(dxy) - dz/lam ;
  double ZNeu = lam*(r +  dxy)  - dz;
  double XNeu = r*cos(phiAnot);
  double YNeu = r*sin(phiAnot);
  //-------------------------------
  //------------------------------
  //all WRT 000

  double phiAnot1  = atan2(A1SV.Y(), A1SV.X());

  double xpoca2Anot1 = dxy*sin(phi0);
  double ypoca2Anot1 = -dxy*cos(phi0);
  double at = (SV-PV).Y()/(SV-PV).X(), bt =SV.Y() - at*SV.X() ;
  double am = tan(phi0), bm =ypoca2Anot1-am*xpoca2Anot1;

  double XNeu1 = (bt-bm)/(am-at);
  double YNeu1 = am*XNeu1 + bm;

  double r1 = sqrt(XNeu1*XNeu1  + YNeu1*YNeu1  );
  double phiWRT000 = atan2(YNeu1,XNeu1);

  double bz01 = fabs(dxy) - dz/lam ;
  double ZNeu1 = lam*(r1 +  dxy)  - dz;

  //---- Covariance of original parameters (dxy, phi0,dz,lambda)
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
  HelixCov(4,4) = sqrt(TauDirError.Y()*TauDirError.Y()/A1SV.X()/A1SV.X()   + TauDirError.X()*TauDirError.X()*A1SV.Y()*A1SV.Y()/A1SV.X()/A1SV.X()/A1SV.X()/A1SV.X() ) *cos(phiAnot)*cos(phiAnot);

  //-------------------
  double drdd = 2*sin(phi0)/(tan(phi0)*cos(phiAnot) - sin(phiAnot));
  double drdphi0 = 2*dxy*cos(phi0)/(tan(phi0)*cos(phiAnot) - sin(phiAnot))   - 2*dxy*sin(phi0)*cos(phiAnot)/cos(phi0)/cos(phi0)/(tan(phi0)*cos(phiAnot) - sin(phiAnot))/(tan(phi0)*cos(phiAnot) - sin(phiAnot));
  double drdphi =2*dxy*sin(phi0)*(cos(phiAnot) + tan(phi0)*sin(phiAnot))/ (tan(phi0)*cos(phiAnot) - sin(phiAnot))/(tan(phi0)*cos(phiAnot) - sin(phiAnot));
  TMatrixT<double>  DerivativesXYZToHelix;

  DerivativesXYZToHelix.ResizeTo(3,5);


  // Set derivatives matrix  (dxy,phi0,lambda,z,phi)
  //  std::cout<<" PARAMETERES:    drdd  "<< drdd <<" drdphi0  "<<drdphi0<<"  drdphi  "<<drdphi<< " lam  "<< lam <<std::endl;
  DerivativesXYZToHelix(0,0)  = drdd*cos(phiAnot);
  DerivativesXYZToHelix(0,1)  = drdphi0*cos(phiAnot);
  DerivativesXYZToHelix(0,2)  = 0;
  DerivativesXYZToHelix(0,3)  = 0;
  DerivativesXYZToHelix(0,4)  = drdphi*cos(phiAnot) - r*sin(phiAnot);

  DerivativesXYZToHelix(1,0)  = drdd*sin(phiAnot);
  DerivativesXYZToHelix(1,1)  = drdphi0*sin(phiAnot);
  DerivativesXYZToHelix(1,2)  = 0;
  DerivativesXYZToHelix(1,3)  = 0;
  DerivativesXYZToHelix(1,4)  = drdphi*sin(phiAnot) - r*cos(phiAnot);

  DerivativesXYZToHelix(2,0)  = lam;
  DerivativesXYZToHelix(2,1)  = lam*drdphi0;
  DerivativesXYZToHelix(2,2)  = r+dxy;
  DerivativesXYZToHelix(2,3)  = -1;
  DerivativesXYZToHelix(2,4)  = lam*drdphi;

  TVector3 NormalisedTauMuDirection;  TVector3 NormalisedTauMuDirectionWRT000;
  TVector3 NormalisedTauA1Direction;


  TVector3 PointWRTPV(XNeu,YNeu,ZNeu);
  TVector3 PointWRT000(XNeu1 - PV.X(),YNeu1 - PV.Y(),ZNeu1 - PV.Z());

  NormalisedTauMuDirection.SetX(PointWRTPV.X()/PointWRTPV.Mag());
  NormalisedTauMuDirection.SetY(PointWRTPV.Y()/PointWRTPV.Mag());
  NormalisedTauMuDirection.SetZ(PointWRTPV.Z()/PointWRTPV.Mag());

  NormalisedTauA1Direction.SetX(TauA1.LV().Px()/TauA1.LV().P());
  NormalisedTauA1Direction.SetY(TauA1.LV().Py()/TauA1.LV().P());
  NormalisedTauA1Direction.SetZ(TauA1.LV().Pz()/TauA1.LV().P());
  

  //  std::cout<<"chek new phi  "<<  atan2(YNeu,XNeu) <<" SV  "<< phiAnot<<" fabs(new Phi)   " <<fabs(phiAnot - atan2(YNeu,XNeu))<<"  " <<atan2(NormalisedTauMuDirection.Y(),  NormalisedTauMuDirection.X() )<<std::endl;

  double cosTauTau = NormalisedTauMuDirection.X()* NormalisedTauA1Direction.X() + NormalisedTauMuDirection.Y()* NormalisedTauA1Direction.Y() + NormalisedTauMuDirection.Z()* NormalisedTauA1Direction.Z();


  TMatrixT<double> DerivativesXYZToHelixT=DerivativesXYZToHelix; DerivativesXYZToHelixT.T();
  TMatrixT<double> CovXYZFrame=DerivativesXYZToHelix*HelixCov*DerivativesXYZToHelixT;
  
  //  double TauA1deltaP = sqrt(   (pow(TauA1.LV().Px()*TauA1.Covariance(3,3) ,2)  +  pow(TauA1.LV().Py()*TauA1.Covariance(4,4),2)  +  pow(TauA1.LV().Pz()*TauA1.Covariance(5,5),2)    )/TauA1.LV().P()/TauA1.LV().P()  )  ;
  double TauA1deltaP = sqrt(   (pow(TauA1.LV().Px(),2)*TauA1.Covariance(3,3)   +  pow(TauA1.LV().Py(),2)*TauA1.Covariance(4,4)  +  pow(TauA1.LV().Pz(),2)*TauA1.Covariance(5,5)    )/TauA1.LV().P()/TauA1.LV().P()  )  ;
  double TauMudeltaP = TauA1deltaP*pow(MassConstraint_/2/TauA1.LV().P(),2)/(1-cosTauTau);
  double TauMuP = MassConstraint_*MassConstraint_/2/(1-cosTauTau)/TauA1.LV().P();

  TLorentzVector TauMuEstimate;
  TauMuEstimate.SetXYZM(TauMuP*NormalisedTauMuDirection.X(),TauMuP*NormalisedTauMuDirection.Y(),TauMuP*NormalisedTauMuDirection.Z(),1.777);

//  std::cout<<"TauA1 covariance "<< sqrt(TauA1.Covariance(3,3))<<" "<<sqrt(TauA1.Covariance(4,4))<<" "<<sqrt(TauA1.Covariance(5,5))<<" TauA1 deltaPa  " <<TauA1deltaP <<" TauMu deltaPa  " <<TauMudeltaP<<std::endl;
//   std::cout<<" TauMuEstimate  "<< TauMuEstimate.Px()<<" "<<TauMuEstimate.Py()<<"  "<<TauMuEstimate.Pz()<<std::endl;
//   std::cout<<" TauA1P  "<< TauA1.LV().Px()<<" "<<TauA1.LV().Py()<<"  "<<TauA1.LV().Pz()<<std::endl;
//   std::cout<<" deltaPhi Mu - TauA1  "<< fabs(TauA1.LV().Phi() - TauMuEstimate.Phi())<<std::endl;
//   std::cout<<" deltaPhi TauMu - TauA1  "<< TauA1.LV().Phi() - phi0  <<std::endl;
//   std::cout<<" DiTauMass  "<< (TauA1.LV() + TauMuEstimate).M()  <<" generated value  " <<MassConstraint_ <<std::endl;



   //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>   Compute Covariance in terms of angles only a

  double br = tan(phi0)*cos(phiAnot) - sin(phiAnot);
  double Zc = ZNeu;

  double MinusSintheta = -sqrt(1  - Zc*Zc/(r*r + Zc*Zc) );

  double cosTheta = Zc/sqrt(r*r + Zc*Zc) ;
  double sinTheta = -MinusSintheta;
  if(fabs(phiAnot - TauA1.LV().Phi()) < 0.05) phiAnot = phiAnot - TMath::Pi();
  //-derivatives

  double drzdd = 2*sin(phi0)/br/Zc  - 2*dxy*sin(phi0)*lam/br/Zc/Zc;
  double drzdphi0 = (2*dxy*cos(phi0)/(tan(phi0)*cos(phiAnot) - sin(phiAnot))   - 2*dxy*sin(phi0)*cos(phiAnot)/cos(phi0)/cos(phi0)/(tan(phi0)*cos(phiAnot) - sin(phiAnot))/(tan(phi0)*cos(phiAnot) - sin(phiAnot)))/Zc;

  double drzdlam = -2*dxy*sin(phi0)*(r+dxy)/br/Zc/Zc;
  double drzdz0 = -2*dxy*sin(phi0)/br/Zc;

  double drzdphi = drdphi*(1 - r*lam/Zc)/Zc;//2*dxy*sin(phi0)*(cos(phiAnot) + tan(phi0)*sin(phiAnot))/ (tan(phi0)*cos(phiAnot) - sin(phiAnot))/(tan(phi0)*cos(phiAnot) - sin(phiAnot));
  double dcosThetadrz = -r/Zc/sqrt(pow(1 + r*r/Zc/Zc,3));

  double dThetadPhi = -dcosThetadrz*drzdphi/sinTheta;
  //-derivatives

  double cosTauTau2 = cos(TauA1.LV().Theta() + acos(cosTheta));
  double sinTauTau2 = sin(TauA1.LV().Theta() + acos(cosTheta));

  double TauMudeltaP_2 = TauA1deltaP*pow(MassConstraint_/2/TauA1.LV().P(),2)/(1-cosTauTau2);

  double TauMuP_2 =  MassConstraint_*MassConstraint_/2/(1-cosTauTau2)/TauA1.LV().P();

  double dPdTheta  = MassConstraint_*MassConstraint_/2/TauA1.LV().P()*sinTauTau2/(1-cosTauTau2)/(1-cosTauTau2);
  double dPdPhi  = dPdTheta*dThetadPhi;

  TLorentzVector TauMuEstimate2;
  TauMuEstimate2.SetXYZM(TauMuP_2*cos(phiAnot)*sinTheta,TauMuP_2*sin(phiAnot)*sinTheta,TauMuP_2*cosTheta,1.777);


  //double dPd

  TMatrixT<double>  DerivativesHelixToAngles;
  DerivativesHelixToAngles.ResizeTo(2,5);

  // Set derivatives matrix  (dxy,phi0,lambda,z,phi)
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

  TMatrixT<double> CovPPhiThetaFrame;
  CovPPhiThetaFrame.ResizeTo(3,3);
//---------- dp dphi dtheta
  CovPPhiThetaFrame(0,0) = TauMudeltaP;
  CovPPhiThetaFrame(0,1) = 0;
  CovPPhiThetaFrame(0,2) = 0;

  CovPPhiThetaFrame(1,0) = 0;
  CovPPhiThetaFrame(1,1) = CovAngleFrame(0,0);
  CovPPhiThetaFrame(1,2) = CovAngleFrame(0,1);

  CovPPhiThetaFrame(2,0) = 0;
  CovPPhiThetaFrame(2,1) = CovAngleFrame(1,0);
  CovPPhiThetaFrame(2,2) = CovAngleFrame(1,1);
//    include correlations

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

   //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   Compute Covariance in terms of angles only a

  TMatrixT<double>    par(LorentzVectorParticle::NLorentzandVertexPar,10);
  TMatrixTSym<double> Cov(LorentzVectorParticle::NLorentzandVertexPar);

  par(LorentzVectorParticle::vx,0)=XNeu;
  par(LorentzVectorParticle::vy,0)=YNeu;
  par(LorentzVectorParticle::vz,0)=ZNeu;
  par(LorentzVectorParticle::px,0)=TauMuEstimate2.Px();
  par(LorentzVectorParticle::py,0)=TauMuEstimate2.Py();
  par(LorentzVectorParticle::pz,0)=TauMuEstimate2.Pz();
  par(LorentzVectorParticle::m,0) =1.777;

  //---- FillVertexCov
  for(int i=0; i<LorentzVectorParticle::NVertex; i++){
    for(int j=0; j<LorentzVectorParticle::NVertex; j++){
      Cov(i,j)=CovXYZFrame(i,j);
    }
  }
  //---- FillVertexCov

  TMatrixT<double>    VectorEForCovariance;
  TMatrixT<double>    VectorRForCovariance;

  TMatrixT<double>    MatrixEForCovariance;
  TMatrixT<double>    FinCovariance;

  VectorEForCovariance.ResizeTo(3,1);
  VectorRForCovariance.ResizeTo(3,1);

  MatrixEForCovariance.ResizeTo(3,3);
  FinCovariance.ResizeTo(3,3);

  //VectorEForCovariance(0,0) = TauMuEstimate.P();
  //VectorEForCovariance(1,0) = TauMuEstimate.P();
  //VectorEForCovariance(2,0) = TauMuEstimate.P();

  //VectorRForCovariance(0,0) = XNeu;
  //VectorRForCovariance(1,0) = YNeu;
  //VectorRForCovariance(2,0) = ZNeu;

  for(int i=0; i<LorentzVectorParticle::NVertex; i++){
	  for(int j=0; j<LorentzVectorParticle::NVertex; j++){
		  MatrixEForCovariance(i,j)=0;
		  if(i==j)MatrixEForCovariance(i,j) = TauMudeltaP*(pow(XNeu,2)+pow(YNeu,2)+pow(ZNeu,2) );
	  }
  }
  FinCovariance = MatrixEForCovariance + TauMuEstimate.P()*TauMuEstimate.P()*CovXYZFrame;

  for(int i=0; i<LorentzVectorParticle::NVertex; i++){
	  for(int j=0; j<LorentzVectorParticle::NVertex; j++){
		  Cov(i+3,j+3)=CovPxPyPzFrame(i,j);
	  }
  }

    return LorentzVectorParticle(par,Cov,PDGInfo::tau_minus,0,0);
 }

LorentzVectorParticle DiTauConstrainedFitter::GetTauMuEstimate(){
return particles_.at(1);
}
