#include "SimpleFits/FitSoftware/interface/TrackHelixVertexFitter.h"
#include "TDecompBK.h"
#include <iostream>

TrackHelixVertexFitter::TrackHelixVertexFitter(std::vector<TrackParticle> &particles_,TVector3 vguess):
  isFit(false),
  nParticles(particles_.size()),
  nPar((NFreeTrackPar-NFreeVertexPar)*particles_.size()+NFreeVertexPar),
  nVal(TrackParticle::NHelixPar*particles_.size())
{
  std::cout << "TrackHelixVertexFitter::TrackHelixVertexFitter nParticles " << nParticles << " nPar " << nPar << " " << " nVal " << nVal << std::endl;
  particles=particles_;
  par.ResizeTo(nPar,1);
  parcov.ResizeTo(nPar,nPar);
  val.ResizeTo(nVal,1);
  cov.ResizeTo(nVal,nVal);
  for(unsigned int p=0; p<particles.size();p++){
    for(unsigned int j=0; j<TrackParticle::NHelixPar;j++){
      val(MeasuredValueIndex(j,p),0)=particles.at(p).Parameter(j);
      std::cout << "Input Par" << val(MeasuredValueIndex(j,p),0) << std::endl;
      for(unsigned int k=0; k<TrackParticle::NHelixPar;k++){
	cov(MeasuredValueIndex(j,p),MeasuredValueIndex(k,p))=particles.at(p).Covariance(j,k);
      }
    }
  }
  TDecompBK Inverter(cov);
  cov_inv.ResizeTo(nVal,nVal);
  cov_inv=Inverter.Invert();
  ndf=nVal-nPar;
  // Set Initial conditions within reason
  par(x0,0)  = vguess.X(); parcov(x0,x0)=pow(5.0,0);
  par(y0,0)  = vguess.Y(); parcov(y0,y0)=pow(5.0,0);
  par(z0,0)  = vguess.Z(); parcov(z0,z0)=pow(50.0,0);
  for(unsigned int p=0; p<particles.size();p++){
    par(FreeParIndex(kappa0,p),0)  = val(MeasuredValueIndex(TrackParticle::kappa,p),0);
    par(FreeParIndex(lambda0,p),0) = val(MeasuredValueIndex(TrackParticle::lambda,p),0);
    par(FreeParIndex(phi0,p),0)    = val(MeasuredValueIndex(TrackParticle::phi,p),0);
    //    
    parcov(FreeParIndex(kappa0,p),FreeParIndex(kappa0,p))   = cov(MeasuredValueIndex(TrackParticle::kappa,p),MeasuredValueIndex(TrackParticle::kappa,p));
    parcov(FreeParIndex(lambda0,p),FreeParIndex(lambda0,p)) = cov(MeasuredValueIndex(TrackParticle::lambda,p),MeasuredValueIndex(TrackParticle::lambda,p));
    parcov(FreeParIndex(phi0,p),FreeParIndex(phi0,p))       = cov(MeasuredValueIndex(TrackParticle::phi,p),MeasuredValueIndex(TrackParticle::phi,p));
  }
  std::cout << "TrackHelixVertexFitter::TrackHelixVertexFitter done" << std::endl;
}

TrackHelixVertexFitter::~TrackHelixVertexFitter(){}

double TrackHelixVertexFitter::UpdateChisquare(TMatrixT<double> inpar){
  std::cout << "TrackHelixVertexFitter::UpdateChisquare" << std::endl;
  TMatrixT<double> vprime=ComputePar(inpar);
  TMatrixT<double> dalpha=vprime-val;
  TMatrixT<double> dalphaT=dalpha;  dalphaT.T();
  TMatrixT<double> chisquare=dalphaT*(cov_inv*dalpha);
  std::cout << "TrackHelixVertexFitter::UpdateChisquare done" << std::endl;
  return chisquare(0,0);
}

std::vector<TrackParticle> TrackHelixVertexFitter::GetReFitTracks(){
  std::cout << "TrackHelixVertexFitter::GetReFitTracks" << std::endl;
  std::vector<TrackParticle> refitParticles;
  for(unsigned int p=0;p<particles.size();p++){
    TMatrixT<double> FreePar(NFreeTrackPar,1);
    TMatrixTSym<double> FreeParCov(NFreeTrackPar);
    for(int i=0;i<FreeParCov.GetNrows();i++){
      FreePar(i,0)=par(FreeParIndex(i,p),0);
      for(int j=0;j<FreeParCov.GetNrows();j++){
	FreeParCov(i,j)=parcov(FreeParIndex(i,p),FreeParIndex(j,p));
      }
    }
    TMatrixT<double>    TrackPar=ComputePar(FreePar);
    TMatrixTSym<double> TrackCov=ErrorMatrixPropagator::PropogateError(&TrackHelixVertexFitter::ComputePar,FreePar,FreeParCov);
    refitParticles.push_back(TrackParticle(TrackPar,TrackCov,particles.at(p).PDGID(),particles.at(p).Mass(),particles.at(p).Charge(),particles.at(p).BField()));
  }
  std::cout << "TrackHelixVertexFitter::GetReFitTracks done" << std::endl;
  return particles;
}



std::vector<LorentzVectorParticle> TrackHelixVertexFitter::GetReFitLorentzVectorParticles(){
  std::cout << "TrackHelixVertexFitter::GetReFitLorentzVectorParticles " <<std::endl;
  std::vector<LorentzVectorParticle> refitParticles;
  for(unsigned int p=0;p<particles.size();p++){
    TMatrixT<double>    FreePar(NFreeTrackPar+NExtraPar+MassOffSet,1);
    TMatrixTSym<double> FreeParCov(NFreeTrackPar+NExtraPar+MassOffSet);
    for(int i=0;i<NFreeTrackPar;i++){
      std::cout << "TrackHelixVertexFitter::GetReFitLorentzVectorParticles i" << i << " " << par.GetNrows() << std::endl;
      FreePar(i,0)=par(FreeParIndex(i,p),0);
      for(int j=0;j<NFreeTrackPar;j++){
        FreeParCov(i,j)=parcov(FreeParIndex(i,p),FreeParIndex(j,p));
      }
    }
    std::cout << "TrackHelixVertexFitter::GetReFitLorentzVectorParticles B" <<std::endl;
    FreePar(NFreeTrackPar+MassOffSet,0)=particles.at(p).Mass();
    FreePar(NFreeTrackPar+BField0,0)=particles.at(p).BField();
    std::cout << "TrackHelixVertexFitter::GetReFitLorentzVectorParticles A" <<std::endl;
    TMatrixT<double>    LVPar=ComputeLorentzVectorPar(FreePar);
    TMatrixTSym<double> LVCov=ErrorMatrixPropagator::PropogateError(&TrackHelixVertexFitter::ComputeLorentzVectorPar,FreePar,FreeParCov);
    refitParticles.push_back(LorentzVectorParticle(LVPar,LVCov,particles.at(p).PDGID(),particles.at(p).Charge(),particles.at(p).BField()));
  }
  std::cout << "TrackHelixVertexFitter::GetReFitLorentzVectorParticles done" <<std::endl;
  return refitParticles;
}

LorentzVectorParticle TrackHelixVertexFitter::GetMother(int pdgid){
  std::cout << "TrackHelixVertexFitter::GetMother" << std::endl;
  double c(0),b(0);
  TMatrixT<double>    FreePar(par.GetNrows()+NExtraPar+particles.size(),1);
  TMatrixTSym<double> FreeParCov(par.GetNrows()+NExtraPar+particles.size());
  std::cout << "TrackHelixVertexFitter::GetMother a" << std::endl; 
  for(int i=0;i<par.GetNrows();i++){
    FreePar(i,0)=par(i,0);
    for(int j=0;j<par.GetNrows();j++){FreeParCov(i,j)=parcov(i,j);}
  }
  std::cout << "TrackHelixVertexFitter::GetMother A" << std::endl;
  for(unsigned int p=0; p<particles.size();p++){
    b=particles.at(p).BField();
    c+=particles.at(p).Charge();
    FreePar(par.GetNrows()+MassOffSet+p,0)=particles.at(p).Mass();
  }
  std::cout << "TrackHelixVertexFitter::GetMother B" << std::endl;
  FreePar(par.GetNrows()+BField0,0)=b;
  std::cout << "TrackHelixVertexFitter::GetMother C" << std::endl;
  TMatrixT<double>    mpar=ComputeMotherLorentzVectorPar(FreePar);
  std::cout << "TrackHelixVertexFitter::GetMother D" << std::endl;
  TMatrixTSym<double> mcov=ErrorMatrixPropagator::PropogateError(&TrackHelixVertexFitter::ComputeMotherLorentzVectorPar,FreePar,FreeParCov);
  std::cout << "TrackHelixVertexFitter::GetMother done" << std::endl;
  return LorentzVectorParticle(mpar,mcov,pdgid,c,b);
}

TVector3 TrackHelixVertexFitter::GetVertex(){
  std::cout << "TrackHelixVertexFitter::GetVertex" << std::endl;
  return TVector3(par(FreeParIndex(x0,0),0),par(FreeParIndex(y0,0),0),par(FreeParIndex(z0,0),0));
}

TMatrixTSym<double> TrackHelixVertexFitter::GetVertexError(){
  std::cout << "TrackHelixVertexFitter::GetVertexError" << std::endl;
  TMatrixTSym<double> c(NFreeVertexPar);
  for(unsigned int i=0;i<NFreeVertexPar;i++){
    for(unsigned int j=0;j<NFreeVertexPar;j++){c(FreeParIndex(i,0),FreeParIndex(j,0));}
  }
  std::cout << "TrackHelixVertexFitter::GetVertexError done" << std::endl;
  return c;
}

void TrackHelixVertexFitter::Computedxydz(TMatrixT<double> &inpar,int p,double &c,double &lam,double &phi,double &x,double &y,double &z,double &s,double &dxy,double &dz){
  std::cout << "TrackHelixVertexFitter::Computedxydz p" << std::endl;
  std::cout << "Size " << inpar.GetNrows() << " FreeParIndex kappa0:" << FreeParIndex(kappa0,p) << " lambda0: " << FreeParIndex(lambda0,p) << " phi0: " << FreeParIndex(phi0,p) << " x0: " << FreeParIndex(x0,p) << " y0: " << FreeParIndex(y0,p) << " z0: " << FreeParIndex(z0,p) << std::endl;    
  c=inpar(FreeParIndex(kappa0,p),0);
  lam=inpar(FreeParIndex(lambda0,p),0);
  phi=inpar(FreeParIndex(phi0,p),0);
  x=inpar(FreeParIndex(x0,p),0);
  y=inpar(FreeParIndex(y0,p),0);
  z=inpar(FreeParIndex(z0,p),0);
  s=1/(2*c)*asin(2*c*(x*cos(phi)+y*sin(phi)));
  dxy=y*cos(phi)-x*sin(phi)-(1/c)*sin(c*s)*sin(c*s);
  dz=z-s*tan(lam);
  std::cout << "TrackHelixVertexFitter::Computedxydz done" << std::endl;
}

TMatrixT<double> TrackHelixVertexFitter::ComputePar(TMatrixT<double> &inpar){
  std::cout << "TrackHelixVertexFitter::ComputePar" << std::endl;
  int nparticles=(inpar.GetNrows()-NFreeVertexPar)/(NFreeTrackPar-NFreeVertexPar);
  TMatrixT<double> helices(nparticles*TrackParticle::NHelixPar,1);
  for(int p=0;p<nparticles;p++){
    TMatrixT<double> TrackPar=ComputeTrackPar(inpar,p);
    for(int i=0;i<TrackParticle::NHelixPar;i++){helices(MeasuredValueIndex(i,p),0)=TrackPar(i,0);}
  }
  std::cout << "TrackHelixVertexFitter::ComputePar done" << std::endl;
  return helices;
}

TMatrixT<double> TrackHelixVertexFitter::ComputeTrackPar(TMatrixT<double> &inpar, int p){
  std::cout << "TrackHelixVertexFitter::ComputeTrackPar" << std::endl; 
  TMatrixT<double> helix(TrackParticle::NHelixPar,1);
  // copy parameters that are 1 to 1
  double c,lam,phi,x,y,z,s,dxy,dz;
  TrackHelixVertexFitter::Computedxydz(inpar,p,c,lam,phi,x,y,z,s,dxy,dz);
  helix(TrackParticle::kappa,0)=c;
  helix(TrackParticle::lambda,0)=lam;
  helix(TrackParticle::phi,0)=phi;
  helix(TrackParticle::dxy,0)=dxy;
  helix(TrackParticle::dz,0)=dz;
  std::cout << "TrackHelixVertexFitter::ComputeTrackPar done" << std::endl;
  return helix;
}

TMatrixT<double> TrackHelixVertexFitter::ComputeLorentzVectorPar(TMatrixT<double> &inpar){
  std::cout << "TrackHelixVertexFitter::ComputeLorentzVectorPar" << std::endl;
  int np(0), parsize(0); ParSizeInfo(inpar,np,parsize,true);
  std::cout << "parsize " << parsize << " np " << np << std::endl;
  double B=inpar(parsize+BField0,0);
  double massHypothesis=inpar(parsize+MassOffSet,0);
  TMatrixT<double> LV(LorentzVectorParticle::NLorentzandVertexPar,1);
  double c,lam,phi,x,y,z,s,dxy,dz;
  int p=0;
  std::cout << "before Computedxydz" << std::endl;
  TrackHelixVertexFitter::Computedxydz(inpar,p,c,lam,phi,x,y,z,s,dxy,dz);
  std::cout << "after Computedxydz" << std::endl;
  LV(LorentzVectorParticle::px,0)=B*(1.0/c)*cos(2*s*c+phi);
  LV(LorentzVectorParticle::py,0)=B*(1.0/c)*sin(2*s*c+phi);
  LV(LorentzVectorParticle::pz,0)=B*(1.0/c)*tan(lam) ;
  LV(LorentzVectorParticle::m,0) =massHypothesis;
  LV(LorentzVectorParticle::vx,0)=x;
  LV(LorentzVectorParticle::vy,0)=y;
  LV(LorentzVectorParticle::vz,0)=z;
  std::cout << "TrackHelixVertexFitter::ComputeLorentzVectorPar done" << std::endl;
  return LV;
}

TMatrixT<double> TrackHelixVertexFitter::ComputeMotherLorentzVectorPar(TMatrixT<double> &inpar){
  std::cout << "TrackHelixVertexFitter::ComputeMotherLorentzVectorPar" << std::endl;
  TMatrixT<double> mother(LorentzVectorParticle::NLorentzandVertexPar,1);
  double E(0);
  int np(0), parsize(0); ParSizeInfo(inpar,np,parsize,true);
  std::cout << "parsize " << parsize << " np " << np << std::endl;
  for(int p=0;p<np;p++){
    std::cout << "p " << p << std::endl;
    TMatrixT<double> particlepar(NFreeTrackPar+NExtraPar+MassOffSet,1);
    for(int i=0;i<NFreeTrackPar;i++){std::cout << " " << inpar(FreeParIndex(i,p),0) << " " << std::endl;particlepar(i,0)=inpar(FreeParIndex(i,p),0);}
    std::cout << "size " << NFreeTrackPar+NExtraPar+MassOffSet << " " << particlepar.GetNrows() << std::endl; 
    particlepar(NFreeTrackPar+BField0,0)=inpar(parsize+BField0,0);
    particlepar(NFreeTrackPar+MassOffSet,0)=inpar(parsize+MassOffSet+p,0);
    std::cout << "Particle Filled" << std::endl;
    TMatrixT<double> daughter=TrackHelixVertexFitter::ComputeLorentzVectorPar(particlepar);
    std::cout << "have daugther 4 vector" << std::endl;
    mother(LorentzVectorParticle::px,0)+=daughter(LorentzVectorParticle::px,0);
    mother(LorentzVectorParticle::py,0)+=daughter(LorentzVectorParticle::py,0);
    mother(LorentzVectorParticle::pz,0)+=daughter(LorentzVectorParticle::pz,0);
    mother(LorentzVectorParticle::vx,0)=daughter(LorentzVectorParticle::vx,0);
    mother(LorentzVectorParticle::vy,0)=daughter(LorentzVectorParticle::vy,0);
    mother(LorentzVectorParticle::vz,0)=daughter(LorentzVectorParticle::vz,0);
    E+=sqrt((daughter(LorentzVectorParticle::px,0)*daughter(LorentzVectorParticle::px,0)+
	     daughter(LorentzVectorParticle::py,0)*daughter(LorentzVectorParticle::py,0)+
	     daughter(LorentzVectorParticle::pz,0)*daughter(LorentzVectorParticle::pz,0)+
	     daughter(LorentzVectorParticle::m,0)*daughter(LorentzVectorParticle::m,0)));
  }
  double P2=(mother(LorentzVectorParticle::px,0)*mother(LorentzVectorParticle::px,0)+
	     mother(LorentzVectorParticle::py,0)*mother(LorentzVectorParticle::py,0)+
	     mother(LorentzVectorParticle::pz,0)*mother(LorentzVectorParticle::pz,0));
  mother(LorentzVectorParticle::m,0)=(E*E-P2)/sqrt(fabs(E*E-P2));
  std::cout << "TrackHelixVertexFitter::ComputeMotherLorentzVectorPar done" << std::endl;
  return mother;
}

TString TrackHelixVertexFitter::FreeParName(int Par){
  int p(0);
  if(Par==x0)     return "x0";
  if(Par==y0)     return "y0";
  if(Par==z0)     return "z0";
  for(p=0;p<nParticles;p++){
    if((Par-NFreeVertexPar)<(p+1)*(NFreeTrackPar-NFreeVertexPar))break;
  }
  TString n;
  int index=Par-p*(NFreeTrackPar-NFreeVertexPar);
  std::cout << "particle " << p << " index " << index << std::endl; 
  if(index==kappa0)  n="kappa0";
  if(index==lambda0) n="lambda0";
  if(index==phi0)    n="phi0";
  n+="_particle";n+=p;
  return n;
}

void TrackHelixVertexFitter::ParSizeInfo(TMatrixT<double> &inpar, int &np, int &parsize, bool hasextras){
  if(hasextras)np=(inpar.GetNrows()-NFreeVertexPar-NExtraPar)/(NFreeTrackPar+MassOffSet-NFreeVertexPar);
  else np=(inpar.GetNrows()-NFreeVertexPar)/(NFreeTrackPar-NFreeVertexPar);
  parsize=np*(NFreeTrackPar-NFreeVertexPar)+NFreeVertexPar;
}
