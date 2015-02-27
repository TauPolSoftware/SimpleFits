#include "SimpleFits/FitSoftware/interface/TrackHelixVertexTools.h"
#include "SimpleFits/FitSoftware/interface/SimpleFits.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"
#include "TDecompBK.h"
#include <iostream>

TrackHelixVertexTools::TrackHelixVertexTools(std::vector<TrackParticle> &particles_,TVector3 &vguess):
  particles(particles_)
{
  int nConstHelicesPar=(NFreeTrackPar-NFreeVertexPar)*particles_.size()+NFreeVertexPar;
  constHelicesPar.ResizeTo(nConstHelicesPar,1);
  constHelicesCov.ResizeTo(nConstHelicesPar,nConstHelicesPar);
  int nHelicesPar=TrackParticle::NHelixPar*particles_.size(); 
  helicesPar.ResizeTo(nHelicesPar,1);
  helicesCov.ResizeTo(nHelicesPar,nHelicesPar);
  
  BField=0;
  for(unsigned int p=0; p<particles.size();p++){
    BField+=particles.at(p).BField();
  }
  if(particles.size()>0)BField/=particles.size();
  // Check that BField is uniform
  for(unsigned int p=0; p<particles.size();p++){
    if(fabs(BField-particles.at(p).BField())>0.00001) Log(Log::Error) << "B-Feild not uniform at 10^{-5}???? B = " << BField 
								      << "(Tesla) deltaB = " << fabs(BField-particles.at(p).BField()) 
								      << "(Tesla) Particle " << p << std::endl; 
  }
  
  for(unsigned int p=0; p<particles.size();p++){
    for(unsigned int j=0; j<TrackParticle::NHelixPar;j++){
      helicesPar(MeasuredHelixIndex(j,p),0)=particles.at(p).Parameter(j);
      for(unsigned int k=0; k<TrackParticle::NHelixPar;k++){
	helicesCov(MeasuredHelixIndex(j,p),MeasuredHelixIndex(k,p))=particles.at(p).Covariance(j,k);
      }
    }
  }

  // Set guess for constrained vertex of constrained helicesPar
  constHelicesPar(x0,0)  = vguess.X(); constHelicesCov(x0,x0)=pow(25.0,2.0);
  constHelicesPar(y0,0)  = vguess.Y(); constHelicesCov(y0,y0)=pow(25.0,2.0);
  constHelicesPar(z0,0)  = vguess.Z(); constHelicesCov(z0,z0)=pow(25.0,2.0);

  for(unsigned int p=0; p<particles.size();p++){
    constHelicesPar(ConstrainedHelixIndex(kappa0,p),0)  = helicesPar(MeasuredHelixIndex(TrackParticle::kappa,p),0);
    constHelicesPar(ConstrainedHelixIndex(lambda0,p),0) = helicesPar(MeasuredHelixIndex(TrackParticle::lambda,p),0);
    constHelicesPar(ConstrainedHelixIndex(phi0,p),0)    = helicesPar(MeasuredHelixIndex(TrackParticle::phi,p),0);
    //    
    constHelicesCov(ConstrainedHelixIndex(kappa0,p),ConstrainedHelixIndex(kappa0,p))   = helicesCov(MeasuredHelixIndex(TrackParticle::kappa,p),MeasuredHelixIndex(TrackParticle::kappa,p));
    constHelicesCov(ConstrainedHelixIndex(lambda0,p),ConstrainedHelixIndex(lambda0,p)) = helicesCov(MeasuredHelixIndex(TrackParticle::lambda,p),MeasuredHelixIndex(TrackParticle::lambda,p));
    constHelicesCov(ConstrainedHelixIndex(phi0,p),ConstrainedHelixIndex(phi0,p))       = helicesCov(MeasuredHelixIndex(TrackParticle::phi,p),MeasuredHelixIndex(TrackParticle::phi,p));
  }


  // Debugging Comments
  for(int i=0;i<helicesPar.GetNrows();i++) Log(Log::Debug) << "HelicesPar: " << helicesPar(i,0) << " " << TrackParticle::Name(i%TrackParticle::NHelixPar) << std::endl;
  for(int i=0;i<helicesCov.GetNrows();i++){
    for(int j=0;j<<helicesCov.GetNrows();j++)  Log(Log::Debug) << helicesCov(i,j) << " ";
    Log(Log::Debug) << std::endl;
  }
  for(int i=0;i<constHelicesPar.GetNrows();i++) Log(Log::Debug) << "constHelicesPar: " << constHelicesPar(i,0) << " " <<  ConstrainedHelicesName(i) << std::endl;
  for(int i=0;i<constHelicesCov.GetNrows();i++){
    for(int j=0;j<constHelicesCov.GetNrows();j++)  Log(Log::Debug) << constHelicesCov(i,j) << " ";   
    Log(Log::Debug) << std::endl;
  }
}

TrackHelixVertexTools::~TrackHelixVertexTools(){}

void TrackHelixVertexTools::ConstrainedHelices(TMatrixT<double> &par, TMatrixTSym<double> &cov){
  if(constHelicesPar.GetNrows()!=par.GetNrows())par.ResizeTo(constHelicesPar.GetNrows(),1);
  par=constHelicesPar;
  if(constHelicesCov.GetNrows()!=cov.GetNrows())cov.ResizeTo(constHelicesCov.GetNrows(),constHelicesCov.GetNrows());
  cov=constHelicesCov;
}

void TrackHelixVertexTools::SetConstrainedHelices(TMatrixT<double>& par,TMatrixTSym<double>& cov){
  if(constHelicesPar.GetNrows()==par.GetNrows() && constHelicesCov.GetNrows()==cov.GetNrows()){
    constHelicesPar=par;
    constHelicesCov=cov;
  }
  else{
    Log(Log::Fatal) << "SetConstrainedHelices invalid sizes constHelicesPar: " << constHelicesPar.GetNrows()
		    << " par: " << par.GetNrows() << " cov: " << cov.GetNrows() << std::endl;
    exit(1);
  }
}

void TrackHelixVertexTools::Helices(TMatrixT<double> &par, TMatrixTSym<double> &cov){
  if(helicesPar.GetNrows()!=par.GetNrows())par.ResizeTo(helicesPar.GetNrows(),1);
  par=helicesPar;
  if(helicesCov.GetNrows()!=cov.GetNrows())cov.ResizeTo(helicesCov.GetNrows(),helicesCov.GetNrows());
  cov=helicesCov;
}

std::vector<TrackParticle> TrackHelixVertexTools::ReFitTracks(){
  std::vector<TrackParticle> refitParticles;
  for(unsigned int p=0;p<particles.size();p++){
    TMatrixT<double> ConstrainedHelix(NFreeTrackPar,1);
    TMatrixTSym<double> ConstrainedHelix_Cov(NFreeTrackPar);
    for(int i=0;i<ConstrainedHelix_Cov.GetNrows();i++){
      ConstrainedHelix(i,0)=constHelicesPar(ConstrainedHelixIndex(i,p),0);
      for(int j=0;j<ConstrainedHelix_Cov.GetNrows();j++){
	ConstrainedHelix_Cov(i,j)=constHelicesCov(ConstrainedHelixIndex(i,p),ConstrainedHelixIndex(j,p));
      }
    }
    TMatrixT<double>    TrackPar=ConvertToHelicesNotation(ConstrainedHelix);
    TMatrixTSym<double> TrackCov=SimpleFits::Get()->PropogateError(this,&TrackHelixVertexTools::ConvertToHelicesNotation,ConstrainedHelix,ConstrainedHelix_Cov);
    refitParticles.push_back(TrackParticle(TrackPar,TrackCov,particles.at(p).PDGID(),particles.at(p).Mass(),particles.at(p).Charge(),particles.at(p).BField()));
  }
  return particles;
}



std::vector<LorentzVectorParticle> TrackHelixVertexTools::ReFitLorentzVectorParticles(){
  std::vector<LorentzVectorParticle> refitParticles;
  for(unsigned int p=0;p<particles.size();p++){
    TMatrixT<double>    ConstrainedHelix(NFreeTrackPar+1,1);
    TMatrixTSym<double> ConstrainedHelix_Cov(NFreeTrackPar+1);
    for(int i=0;i<NFreeTrackPar;i++){
      ConstrainedHelix(i,0)=constHelicesPar(ConstrainedHelixIndex(i,p),0);
      for(int j=0;j<NFreeTrackPar;j++){
        ConstrainedHelix_Cov(i,j)=constHelicesCov(ConstrainedHelixIndex(i,p),ConstrainedHelixIndex(j,p));
      }
    }
    ConstrainedHelix(NFreeTrackPar,0)=particles.at(p).Mass();
    TMatrixT<double>    LVPar=ConvertToLorentzVectorNotation(ConstrainedHelix);
    TMatrixTSym<double> LVCov=SimpleFits::Get()->PropogateError(this,&TrackHelixVertexTools::ConvertToLorentzVectorNotation,ConstrainedHelix,ConstrainedHelix_Cov);

    refitParticles.push_back(LorentzVectorParticle(LVPar,LVCov,particles.at(p).PDGID(),particles.at(p).Charge(),particles.at(p).BField()));
  }
  return refitParticles;
}

LorentzVectorParticle TrackHelixVertexTools::MotherLorentzVectorParticle(int pdgid){
  double c(0),b(0);
  TMatrixT<double>    ConstrainedHelix(constHelicesPar.GetNrows()+NExtraPar+particles.size(),1);
  TMatrixTSym<double> ConstrainedHelix_Cov(constHelicesPar.GetNrows()+NExtraPar+particles.size());
  for(int i=0;i<constHelicesPar.GetNrows();i++){
    ConstrainedHelix(i,0)=constHelicesPar(i,0);
    for(int j=0;j<constHelicesPar.GetNrows();j++){ConstrainedHelix_Cov(i,j)=constHelicesCov(i,j);}
  }
  for(unsigned int p=0; p<particles.size();p++){
    b=particles.at(p).BField();
    c+=particles.at(p).Charge();
  }
  TMatrixT<double>    mpar=ComputeMotherLorentzVectorPar(ConstrainedHelix);
  TMatrixTSym<double> mcov=SimpleFits::Get()->PropogateError(this,&TrackHelixVertexTools::ComputeMotherLorentzVectorPar,ConstrainedHelix,ConstrainedHelix_Cov);
  return LorentzVectorParticle(mpar,mcov,pdgid,c,b);
}

TVector3 TrackHelixVertexTools::Vertex(){
  return TVector3(constHelicesPar(ConstrainedHelixIndex(x0,0),0),constHelicesPar(ConstrainedHelixIndex(y0,0),0),constHelicesPar(ConstrainedHelixIndex(z0,0),0));
}

TMatrixTSym<double> TrackHelixVertexTools::VertexError(){
  TMatrixTSym<double> c(NFreeVertexPar);
  for(unsigned int i=0;i<NFreeVertexPar;i++){
    for(unsigned int j=0;j<NFreeVertexPar;j++){c(ConstrainedHelixIndex(i,0),ConstrainedHelixIndex(j,0))=constHelicesCov(ConstrainedHelixIndex(i,0),ConstrainedHelixIndex(j,0));}
  }
  return c;
}

TMatrixT<double> TrackHelixVertexTools::ConvertToHelicesNotation(TMatrixT<double> &inpar){
  int nparticles=(inpar.GetNrows()-NFreeVertexPar)/(NFreeTrackPar-NFreeVertexPar);
  TMatrixT<double> helicesPar(nparticles*TrackParticle::NHelixPar,1);
  for(int p=0;p<nparticles;p++){
    TMatrixT<double> helix=ConvertToHelixNotation(inpar,p);
    for(int i=0;i<TrackParticle::NHelixPar;i++){helicesPar(MeasuredHelixIndex(i,p),0)=helix(i,0);}
  }
  return helicesPar;
}

void TrackHelixVertexTools::ComputeHelixParameters(TMatrixT<double> &inpar,double &kappa,double &lam,double &phi,double &dxy, double &dz, double &s, unsigned int p){
  kappa=inpar(ConstrainedHelixIndex(kappa0,p),0);
  lam=inpar(ConstrainedHelixIndex(lambda0,p),0);
  phi=inpar(ConstrainedHelixIndex(phi0,p),0);
  double x=inpar(ConstrainedHelixIndex(x0,p),0);
  double y=inpar(ConstrainedHelixIndex(y0,p),0);
  double z=inpar(ConstrainedHelixIndex(z0,p),0);
  double v=(2.0*kappa*(x*cos(phi)+y*sin(phi)));
  double arcsinv=0;
  if(v>=1.0){arcsinv=TMath::Pi()/2;}
  else if(v<=-1.0){arcsinv=-TMath::Pi()/2;}
  else{arcsinv=asin(v);}
  s=1.0/(2.0*kappa)*arcsinv;//asin(2.0*kappa*(x*cos(phi)+y*sin(phi))); // check paper                
  dxy=y*cos(phi)-x*sin(phi)-(1/kappa)*sin(kappa*s)*sin(kappa*s);
  dz=z-s*tan(lam);

  Log(Log::Debug) << "TrackHelixVertexFitter::Computedxydz done" 
		  << "\n kappa0 " << kappa
		  << "\n lambda0 " << lam
		  << "\n phi0 " << phi
		  << "\n x0 " << x
		  << "\n y0 " << y
		  << "\n z0 " << z
		  << "\n Computation: "
		  << "\n arcsin " << asin(2*kappa*(x*cos(phi)+y*sin(phi))) 
		  << "\n cos(phi) " << cos(phi) 
		  << "\n sin(phi) " << sin(phi) 
		  << "\n F " << x*cos(phi)+y*sin(phi) 
		  << "\n s " << s 
		  << "\n v " << v 
		  << "\n dxy " << dxy 
		  << "\n dz "  << dz 
		  << std::endl;
}

TMatrixT<double> TrackHelixVertexTools::ConvertToHelixNotation(TMatrixT<double> &inpar, int p){
  TMatrixT<double> helix(TrackParticle::NHelixPar,1);
  double kappa, lam, phi, dxy, dz, s;
  ComputeHelixParameters(inpar,kappa,lam,phi,dxy,dz,s);
  helix(TrackParticle::kappa,0)=kappa;
  helix(TrackParticle::lambda,0)=lam;
  helix(TrackParticle::phi,0)=phi;
  helix(TrackParticle::dxy,0)=dxy;
  helix(TrackParticle::dz,0)=dz;
  return helix;
}

TMatrixT<double> TrackHelixVertexTools::ConvertToLorentzVectorNotation(TMatrixT<double> &inpar){
  TMatrixT<double> LV(LorentzVectorParticle::NLorentzandVertexPar,1);
  double kappa, lam, phi, dxy, dz, s;
  ComputeHelixParameters(inpar,kappa,lam,phi,dxy,dz,s);
  LV(LorentzVectorParticle::px,0)=BField*(1.0/fabs(kappa))*cos(2*s*kappa+phi);
  LV(LorentzVectorParticle::py,0)=BField*(1.0/fabs(kappa))*sin(2*s*kappa+phi);
  LV(LorentzVectorParticle::pz,0)=BField*(1.0/fabs(kappa))*tan(lam) ;
  LV(LorentzVectorParticle::m,0) =inpar(NFreeTrackPar,0);
  LV(LorentzVectorParticle::vx,0)=inpar(ConstrainedHelixIndex(x0),0);
  LV(LorentzVectorParticle::vy,0)=inpar(ConstrainedHelixIndex(y0),0);
  LV(LorentzVectorParticle::vz,0)=inpar(ConstrainedHelixIndex(z0),0);
  return LV;
}

TMatrixT<double> TrackHelixVertexTools::ComputeMotherLorentzVectorPar(TMatrixT<double> &inpar){
  TMatrixT<double> mother(LorentzVectorParticle::NLorentzandVertexPar,1);
  double E(0);
  for(unsigned int p=0;p<particles.size();p++){
    TMatrixT<double> particlepar(NFreeTrackPar+1,1);
    for(int i=0;i<NFreeTrackPar;i++){particlepar(i,0)=inpar(ConstrainedHelixIndex(i,p),0);}
    particlepar(NFreeTrackPar,0)=particles.at(p).Mass();
    TMatrixT<double> daughter=TrackHelixVertexTools::ConvertToLorentzVectorNotation(particlepar);
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

TString TrackHelixVertexTools::ConstrainedHelicesName(int Par){
  int p(0);
  if(Par==x0)     return "x0";
  if(Par==y0)     return "y0";
  if(Par==z0)     return "z0";
  for(p=0;p<(int)particles.size();p++){
    if((Par-NFreeVertexPar)<(p+1)*(NFreeTrackPar-NFreeVertexPar))break;
  }
  TString n;
  int index=Par-p*(NFreeTrackPar-NFreeVertexPar);
  if(index==kappa0)  n="kappa0";
  if(index==lambda0) n="lambda0";
  if(index==phi0)    n="phi0";
  n+="_particle";n+=p;
  return n;
}
