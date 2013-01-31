#include "SimpleFits/FitSoftware/interface/TrackHelixVertexFitter.h"
#include "TDecompBK.h"


TrackHelixVertexFitter::TrackHelixVertexFitter(std::vector<TrackParticle> &particles_):
  isFit(false),
  nParticles(particles_.size()),
  nPar(3*particles_.size()+3),
  nVal(5*particles_.size())
{
  particles=particles_;
  par.ResizeTo(nPar,1);
  parcov.ResizeTo(nPar,1);
  val.ResizeTo(nVal,1);
  cov.ResizeTo(nVal,nPar);
  for(unsigned int i=0; i<particles.size();i++){
    for(unsigned int j=0; j<TrackParticle::NHelixPar;j++){
      val(MeasuredValueIndex(j,i),0)=particles.at(i).Parameter(j);
      for(unsigned int k=0; k<TrackParticle::NHelixPar;k++){
	cov(MeasuredValueIndex(j,i),MeasuredValueIndex(k,i))=particles.at(i).Covariance(j,k);
      }
    }
  }
  TDecompBK Inverter(cov);
  cov_inv.ResizeTo(nVal,nPar);
  cov_inv=Inverter.Invert();
  ndf=nVal-nPar;
}

TrackHelixVertexFitter::~TrackHelixVertexFitter(){}

double TrackHelixVertexFitter::UpdateChisquare(TMatrixT<double> inpar){
  TMatrixT<double> vprime=ComputePar(inpar);
  TMatrixT<double> dalpha=vprime-val;
  TMatrixT<double> dalphaT=dalpha;  dalphaT.T();
  TMatrixT<double> chisquare=dalphaT*(cov_inv*dalpha);
  return chisquare(0,0);
}

std::vector<TrackParticle> TrackHelixVertexFitter::GetReFitTracks(){
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
  return particles;
}



std::vector<LorentzVectorParticle> TrackHelixVertexFitter::GetReFitLorentzVectorParticles(){
  std::vector<LorentzVectorParticle> refitParticles;
  for(unsigned int p=0;p<particles.size();p++){
    TMatrixT<double>    FreePar(NFreeTrackPar+NExtraPar+particles.size(),1);
    TMatrixTSym<double> FreeParCov(NFreeTrackPar+NExtraPar+particles.size());
    for(int i=0;i<NFreeTrackPar;i++){
      FreePar(i,0)=par(FreeParIndex(i,p),0);
      for(int j=0;j<NFreeTrackPar;j++){
        FreeParCov(i,j)=parcov(FreeParIndex(i,p),FreeParIndex(j,p));
      }
    }
    FreePar(par.GetNrows()+BField0,0)=particles.at(p).BField();
    FreePar(par.GetNrows()+MassOffSet,0)=particles.at(p).Mass();
    TMatrixT<double>    LVPar=ComputeLorentzVectorPar(FreePar);
    TMatrixTSym<double> LVCov=ErrorMatrixPropagator::PropogateError(&TrackHelixVertexFitter::ComputeLorentzVectorPar,FreePar,FreeParCov);
    refitParticles.push_back(LorentzVectorParticle(LVPar,LVCov,particles.at(p).PDGID(),particles.at(p).Charge(),particles.at(p).BField()));
  }
  return refitParticles;
}

LorentzVectorParticle TrackHelixVertexFitter::GetMother(int pdgid){
  double c(0),b(0);
  TMatrixT<double>    FreePar(par.GetNrows()+NExtraPar+particles.size(),1);
  TMatrixTSym<double> FreeParCov(par.GetNrows()+NExtraPar+particles.size());
  for(int i=0;i<FreePar.GetNrows()-NExtraPar;i++){
    FreePar(i,0)=par(i,0);
    for(int j=0;j<FreePar.GetNrows()-NExtraPar;j++){FreeParCov(i,j)=parcov(i,j);}
  }
  for(unsigned int p=0; p<particles.size();p++){
    b=particles.at(p).BField();
    c+=particles.at(p).Charge();
    FreePar(par.GetNrows()+MassOffSet+p,0)=particles.at(p).Mass();
  }
  FreePar(par.GetNrows()+BField0,0)=b;
  TMatrixT<double>    mpar=ComputeMotherLorentzVectorPar(FreePar);
  TMatrixTSym<double> mcov=ErrorMatrixPropagator::PropogateError(&TrackHelixVertexFitter::ComputeMotherLorentzVectorPar,FreePar,FreeParCov);
  return LorentzVectorParticle(mpar,mcov,pdgid,c,b);
}

TVector3 TrackHelixVertexFitter::GetVertex(){
  return TVector3(par(FreeParIndex(x0,0),0),par(FreeParIndex(y0,0),0),par(FreeParIndex(z0,0),0));
}

TMatrixTSym<double> TrackHelixVertexFitter::GetVertexError(){
  TMatrixTSym<double> c(NFreeVertexPar);
  for(unsigned int i=0;i<NFreeVertexPar;i++){
    for(unsigned int j=0;j<NFreeVertexPar;j++){c(FreeParIndex(i,0),FreeParIndex(j,0));}
  }
  return c;
}

void TrackHelixVertexFitter::Computedxydz(TMatrixT<double> &inpar,int p,double &c,double &lam,double &phi,double &x,double &y,double &z,double &s,double &dxy,double &dz){

  c=inpar(FreeParIndex(kappa0,p),0);
  lam=inpar(FreeParIndex(lambda0,p),0);
  phi=inpar(FreeParIndex(phi0,p),0);
  x=inpar(FreeParIndex(x0,p),0);
  y=inpar(FreeParIndex(y0,p),0);
  z=inpar(FreeParIndex(z0,p),0);
  s=1/(2*c)*asin(2*c*(x*cos(phi)+y*sin(phi)));
  dxy=y*cos(phi)-x*sin(phi)-(1/c)*sin(c*s)*sin(c*s);
  dz=z-s*tan(lam);
}

TMatrixT<double> TrackHelixVertexFitter::ComputePar(TMatrixT<double> &inpar){
  int nparticles=(inpar.GetNrows()-NFreeVertexPar)/(NFreeTrackPar-NFreeVertexPar);
  TMatrixT<double> helices(nparticles*TrackParticle::NHelixPar,1);
  for(int p=0;p<nparticles;p++){
    TMatrixT<double> TrackPar=ComputeTrackPar(inpar,p);
    for(int i=0;i<TrackParticle::NHelixPar;i++){helices(MeasuredValueIndex(i,p),0)=TrackPar(i,0);}
  }
  return helices;
}

TMatrixT<double> TrackHelixVertexFitter::ComputeTrackPar(TMatrixT<double> &inpar, int p){
  TMatrixT<double> helix(TrackParticle::NHelixPar,1);
  // copy parameters that are 1 to 1
  double c,lam,phi,x,y,z,s,dxy,dz;
  TrackHelixVertexFitter::Computedxydz(inpar,p,c,lam,phi,x,y,z,s,dxy,dz);
  helix(TrackParticle::kappa,0)=c;
  helix(TrackParticle::lambda,0)=lam;
  helix(TrackParticle::phi,0)=phi;
  helix(TrackParticle::dxy,0)=dxy;
  helix(TrackParticle::dz,0)=dz;
  return helix;
}

TMatrixT<double> TrackHelixVertexFitter::ComputeLorentzVectorPar(TMatrixT<double> &inpar){
  int parsize=NFreeTrackPar+NFreeVertexPar;
  double B=inpar(parsize+BField0,0);
  double massHypothesis=inpar(parsize+MassOffSet,0);
  TMatrixT<double> LV(LorentzVectorParticle::NLorentzandVertexPar,1);
  double c,lam,phi,x,y,z,s,dxy,dz;
  int p=0;
  TrackHelixVertexFitter::Computedxydz(inpar,p,c,lam,phi,x,y,z,s,dxy,dz);
  LV(LorentzVectorParticle::px,0)=B*(1.0/c)*cos(2*s*c+phi);
  LV(LorentzVectorParticle::py,0)=B*(1.0/c)*sin(2*s*c+phi);
  LV(LorentzVectorParticle::pz,0)=B*(1.0/c)*tan(lam) ;
  LV(LorentzVectorParticle::m,0) =massHypothesis;
  LV(LorentzVectorParticle::vx,0)=x;
  LV(LorentzVectorParticle::vy,0)=y;
  LV(LorentzVectorParticle::vz,0)=z;
  return LV;
}

TMatrixT<double> TrackHelixVertexFitter::ComputeMotherLorentzVectorPar(TMatrixT<double> &inpar){
  TMatrixT<double> mother(LorentzVectorParticle::NLorentzandVertexPar,1);
  double E(0);
  int parsize=NFreeTrackPar+NFreeVertexPar;
  int np=(inpar.GetNrows()-NFreeVertexPar-NExtraPar)/(NFreeTrackPar+MassOffSet-NFreeVertexPar);
  for(int p=0;p<np;p++){
    TMatrixT<double> particlepar(NFreeTrackPar+NExtraPar+MassOffSet,1);
    for(int i=0;i<NFreeTrackPar;i++){particlepar(i,0)=inpar(FreeParIndex(i,p),0);}
    particlepar(particlepar.GetNrows()+BField0,0)=inpar(parsize+BField0,0);
    particlepar(particlepar.GetNrows()+MassOffSet,0)=inpar(parsize+MassOffSet+p,0);
    TMatrixT<double> daughter=TrackHelixVertexFitter::ComputeLorentzVectorPar(particlepar);
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
  return mother;
}
