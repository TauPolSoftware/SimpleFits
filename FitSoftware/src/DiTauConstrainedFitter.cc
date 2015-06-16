#include "SimpleFits/FitSoftware/interface/DiTauConstrainedFitter.h"
#include "SimpleFits/FitSoftware/interface/PDGInfo.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"
#include <iostream>

double DiTauConstrainedFitter::MassConstraint_ = 91.5;

DiTauConstrainedFitter::DiTauConstrainedFitter(LorentzVectorParticle TauA1,TrackParticle MuTrack, double phiz, TVector3 PVertex, TMatrixTSym<double> VertexCov){
  DiTauConstrainedFitter(TauA1, MuTrack, phiz, PVertex, VertexCov, 91.5);
}

DiTauConstrainedFitter::DiTauConstrainedFitter(LorentzVectorParticle TauA1,TrackParticle MuTrack, double phiz, TVector3 PVertex, TMatrixTSym<double> VertexCov, double MassConstraint):
  LagrangeMultipliersFitter()
{
  debug = false;
  AnalyticalCovariance =false;

  phiz_ = phiz;
  ptz_ = 0;
  useFullRecoil_ = false;
  MassConstraint_ = MassConstraint;

  Configure(TauA1, MuTrack, PVertex, VertexCov);
}

DiTauConstrainedFitter::DiTauConstrainedFitter(LorentzVectorParticle TauA1,TrackParticle MuTrack, PTObject METminusNeutrino, TVector3 PVertex, TMatrixTSym<double> VertexCov, double MassConstraint):
  LagrangeMultipliersFitter()
{
  debug = false;
  AnalyticalCovariance =false;

  METminusNeutrino_ = METminusNeutrino;
  useFullRecoil_ = true;
  MassConstraint_ = MassConstraint;

  Configure(TauA1, MuTrack, PVertex, VertexCov);
}

void DiTauConstrainedFitter::Configure(LorentzVectorParticle TauA1,TrackParticle MuTrack, TVector3 PVertex, TMatrixTSym<double> VertexCov){
  debug = false;
  AnalyticalCovariance =false;

  LorentzVectorParticle  TauMuGuess;
  if(!useFullRecoil_) TauMuGuess  = TauMuStartingPoint( MuTrack,TauA1,PVertex, VertexCov, TauA1.Vertex(),TauA1.VertexCov());
  else{
	TauMuGuess  = TauMuStartingPointwithFullRecoil(MuTrack,TauA1, METminusNeutrino_, PVertex, VertexCov, TauA1.Vertex(),TauA1.VertexCov());
	phiz_ = (TauA1.LV() + TauMuGuess.LV()).Phi();
  }

  Logger(Logger::Debug) << "Passed: TauMuStartingPoint" << std::endl;

  RecoilX_ = TauA1.LV().X() + TauMuGuess.LV().X();
  RecoilY_ = TauA1.LV().Y() + TauMuGuess.LV().Y();
  ThetaForConstrTemporaryIMplementation_=TauMuGuess.LV().Theta();
  particles_.push_back(TauA1);
  particles_.push_back(TauMuGuess);

  particles0_.push_back(TauA1);
  particles0_.push_back(TauMuGuess);

  Logger(Logger::Debug) << "(TauA1.LV() + TauMuGuess.LV()).M(): " << (TauA1.LV() + TauMuGuess.LV()).M() << std::endl;


  isconfigured=false;

  int size=particles_.size()*3;  
  int sizeTrunc = 3; 
  TMatrixT<double>    inpar(size,1);
  TMatrixTSym<double> incov(size);
 
  // Get primary vertex information
  if(VertexCov.GetNrows()!=LorentzVectorParticle::NVertex)return;


  //set input paramterts:  TauA1 - TauMu
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

 // store expanded par for computation of final par (assumes fit has neglegible impact on a1 correlations with vertex uncertainties)

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

  // set up inital point for fit (cov handled in Fit() function)
  par.ResizeTo(npar);
  par=par_0;


  TMatrixT<double> PARa_0(sizeTrunc,1);
  para_0.ResizeTo(sizeTrunc);
  cova_0.ResizeTo(sizeTrunc,sizeTrunc);
  PARa_0=ComputeExpParToPara(exppar);
  for(int i=0; i<npartr;i++)para_0(i)=PARa_0(i,0);

  cova_0=ErrorMatrixPropagator::PropagateError(&DiTauConstrainedFitter::ComputeExpParToPara,exppar,expcov);


  para.ResizeTo(npartr);
  para=para_0;
  TMatrixT<double> PARb_0(sizeTrunc,1);
  parb_0.ResizeTo(sizeTrunc);
  covb_0.ResizeTo(sizeTrunc,sizeTrunc);
  PARb_0=ComputeExpParToParb(exppar);
  for(int i=0; i<npartr;i++)parb_0(i)=PARb_0(i,0);

  covb_0=ErrorMatrixPropagator::PropagateError(&DiTauConstrainedFitter::ComputeExpParToParb,exppar,expcov);
  parb.ResizeTo(npartr);
  parb=parb_0;

  isconfigured=true;
  Init_Resonance_ = GetMother();

}

TMatrixT<double> DiTauConstrainedFitter::ComputeInitalExpPar(TMatrixT<double> &inpar){
  TMatrixT<double> outpar(6,1);
  int offset=0;//LorentzVectorParticle::NVertex;// for TauA1
  outpar(taua1_px,0)=inpar(LorentzVectorParticle::vx+offset,0);
  outpar(taua1_py,0)=inpar(LorentzVectorParticle::vy+offset,0);
  outpar(taua1_pz,0)=inpar(LorentzVectorParticle::vz+offset,0);

  offset+=3;//LorentzVectorParticle::NLorentzandVertexPar; // for TauNu

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
TMatrixT<double> DiTauConstrainedFitter::ComputeExpParToPara(TMatrixT<double> &inpar){
  TMatrixT<double> outpar(npartr,1);
  for(int i=0;i<npartr;i++){outpar(i,0)=inpar(i,0);}  
  return outpar;
}

TMatrixT<double> DiTauConstrainedFitter::ComputeExpParToParb(TMatrixT<double> &inpar){
  TMatrixT<double> outpar(npartr,1);
  int offset = 3; 
  for(int i=0;i<npartr;i++){outpar(i,0)=inpar(i+offset,0);}
  return outpar;
}



TMatrixT<double> DiTauConstrainedFitter::ComputeTauMuLorentzVectorPar(TMatrixT<double> &inpar){
  //  start with index 3 to fill only momenta part
  TMatrixT<double> outpar(7,1);
  outpar(LorentzVectorParticle::vx,0)=0;
  outpar(LorentzVectorParticle::vy,0)=0;
  outpar(LorentzVectorParticle::vz,0)=0;
  
  outpar(LorentzVectorParticle::px,0)=inpar(3,0);
  outpar(LorentzVectorParticle::py,0)=inpar(4,0);
  outpar(LorentzVectorParticle::pz,0)=inpar(5,0);
  outpar(LorentzVectorParticle::m,0)=1.777;


  return outpar;
}

TMatrixT<double> DiTauConstrainedFitter::ComputeTauA1LorentzVectorPar(TMatrixT<double> &inpar){
  //  start with index 3 to fill only momenta part
  TMatrixT<double> outpar(7,1);
  outpar(LorentzVectorParticle::vx,0)=0;
  outpar(LorentzVectorParticle::vy,0)=0;
  outpar(LorentzVectorParticle::vz,0)=0;
  outpar(LorentzVectorParticle::px,0)=inpar(0,0);
  outpar(LorentzVectorParticle::py,0)=inpar(1,0);
  outpar(LorentzVectorParticle::pz,0)=inpar(2,0);
  outpar(LorentzVectorParticle::m,0) =1.777;  

  return outpar;
}

TMatrixT<double> DiTauConstrainedFitter::ComputeMotherLorentzVectorPar(TMatrixT<double> &inpar){

  TMatrixT<double> outpar(7,1);
  TMatrixT<double> Taumupar=ComputeTauMuLorentzVectorPar(inpar);
  TMatrixT<double> Taua1par=ComputeTauA1LorentzVectorPar(inpar);
  outpar(LorentzVectorParticle::px,0)=Taumupar(LorentzVectorParticle::px,0)+Taua1par(LorentzVectorParticle::px,0);
  outpar(LorentzVectorParticle::py,0)=Taumupar(LorentzVectorParticle::py,0)+Taua1par(LorentzVectorParticle::py,0);
  outpar(LorentzVectorParticle::pz,0)=Taumupar(LorentzVectorParticle::pz,0)+Taua1par(LorentzVectorParticle::pz,0);

  double Etaumu2=pow(Taumupar(LorentzVectorParticle::px,0),2.0)+pow(Taumupar(LorentzVectorParticle::py,0),2.0)+pow(Taumupar(LorentzVectorParticle::pz,0),2.0)+pow(Taumupar(LorentzVectorParticle::m,0),2.0);
  double Etaua12=pow(Taua1par(LorentzVectorParticle::px,0),2.0)+pow(Taua1par(LorentzVectorParticle::py,0),2.0)+pow(Taua1par(LorentzVectorParticle::pz,0),2.0)+pow(Taua1par(LorentzVectorParticle::m,0),2.0);
  double P2=pow(outpar(LorentzVectorParticle::px,0),2.0)+pow(outpar(LorentzVectorParticle::py,0),2.0)+pow(outpar(LorentzVectorParticle::pz,0),2.0);
  outpar(LorentzVectorParticle::m,0)=sqrt(fabs(pow(sqrt(Etaumu2)+sqrt(Etaua12),2.0)-P2));

  return outpar;
}

void DiTauConstrainedFitter::UpdateExpandedPar(){
  for(int i=0; i<npartr;i++){
	  exppar(i,0)=para(i);
	  for(int j=0; j<npartr;j++){expcov(i,j)=cova_0(i,j);}
  }
  int offset = npartr;
  for(int i=0; i<npartr;i++){
	  exppar(i+offset,0)=parb(i);
	  for(int j=0; j<npartr;j++){expcov(i+offset,j+offset)=covb_0(i,j);}
  }

}

std::vector<LorentzVectorParticle> DiTauConstrainedFitter::GetReFitDaughters(){
  std::vector<LorentzVectorParticle> refitParticles;
  UpdateExpandedPar();

  double c(0),b(0);
  for(unsigned int i=0;i<particles_.size();i++){c+=particles_.at(i).Charge();b=particles_.at(i).BField();}
  TMatrixT<double> a1=ComputeTauA1LorentzVectorPar(exppar);
  TMatrixTSym<double> a1cov=ErrorMatrixPropagator::PropagateError(&DiTauConstrainedFitter::ComputeTauA1LorentzVectorPar,exppar,expcov);
     
   for(int i=0; i<LorentzVectorParticle::NVertex; i++){
     for(int j=0; j<LorentzVectorParticle::NVertex; j++){
       a1cov(i,j)=particles_.at(0).VertexCov()(i,j);
     }
   }


    refitParticles.push_back(LorentzVectorParticle(a1,a1cov,PDGInfo::tau_plus,c,b));

    TMatrixT<double> mu=ComputeTauMuLorentzVectorPar(exppar);
    mu(0,0)= particles_.at(1).Parameter(LorentzVectorParticle::vx);
    mu(1,0)= particles_.at(1).Parameter(LorentzVectorParticle::vy);
    mu(2,0)= particles_.at(1).Parameter(LorentzVectorParticle::vz);
  
    TMatrixTSym<double> mucov=ErrorMatrixPropagator::PropagateError(&DiTauConstrainedFitter::ComputeTauMuLorentzVectorPar,exppar,expcov);
    for(int i=0; i<LorentzVectorParticle::NVertex; i++){
      for(int j=0; j<LorentzVectorParticle::NVertex; j++){
	mucov(i,j)=particles_.at(1).VertexCov()(i,j);
      }
    }
    
    refitParticles.push_back(LorentzVectorParticle(mu,mucov,PDGInfo::tau_minus,0.0,b));
   
  return refitParticles; 
}

LorentzVectorParticle 
DiTauConstrainedFitter::GetMother(){
  UpdateExpandedPar();
  double c(0),b(0);
  for(unsigned int i=0;i<particles_.size();i++){c+=particles_.at(i).Charge();b=particles_.at(i).BField();}
  TMatrixT<double> m=ComputeMotherLorentzVectorPar(exppar);
  TMatrixTSym<double> mcov=ErrorMatrixPropagator::PropagateError(&DiTauConstrainedFitter::ComputeMotherLorentzVectorPar,exppar,expcov);

  return LorentzVectorParticle(m,mcov,PDGInfo::Z0,c,b);
}

TVectorD 
DiTauConstrainedFitter::HardValue(TVectorD &va,TVectorD &vb){
  TLorentzVector Taua1,Taumu; 
  double ZMass;
  CovertParToObjects(va,vb,Taua1,Taumu,ZMass);


  TLorentzVector z=Taua1+Taumu; 
  TVectorD d(NConstraints());

  d(0) = z.M() - MassConstraint_;
  d(1) = Taumu.Pz()/Taumu.P() - cos(ThetaForConstrTemporaryIMplementation_);
  return d;
} 


TVectorD 
DiTauConstrainedFitter::SoftValue(TVectorD &va,TVectorD &vb){
  TLorentzVector Taua1,Taumu; 
  double ZMass;
  CovertParToObjects(va,vb,Taua1,Taumu,ZMass);
  TVectorD d(NSoftConstraints());


  if(!useFullRecoil_){
	d(0) = Taua1.Px() + Taumu.Px();
	d(1) = Taua1.Py() + Taumu.Py();
	d(2) = ( (Taua1.Py() + Taumu.Py())/(Taua1.Px() + Taumu.Px())) -  tan(phiz_);
  }
  else{
	d(0) = Taua1.Px() + Taumu.Px() - RecoilX_;
	d(1) = Taua1.Py() + Taumu.Py() - RecoilY_;
	d(2) = ( (Taua1.Py() + Taumu.Py())/(Taua1.Px() + Taumu.Px())) -  tan(phiz_);
  }
  return d;
} 

void DiTauConstrainedFitter::CovertParToObjects(TVectorD &va,TVectorD &vb,TLorentzVector &Taua1,TLorentzVector &Taumu, double &Zmass){
  // Taua1=particles_.at(0).LV();//TLorentzVector(v(taua1_px),v(taua1_py),v(taua1_pz),sqrt(1.777*1.777+v(taua1_px)*v(taua1_px)+v(taua1_py)*v(taua1_py)+v(taua1_pz)*v(taua1_pz)));
  // Taumu=TLorentzVector(v(taua1_px),v(taua1_py),v(taua1_pz),sqrt(1.777*1.777+v(taua1_px)*v(taua1_px)+v(taua1_py)*v(taua1_py)+v(taua1_pz)*v(taua1_pz)));
  Taua1=TLorentzVector(va(tau_px),va(tau_py),va(tau_pz),sqrt(1.777*1.777+va(tau_px)*va(tau_px)+va(tau_py)*va(tau_py)+va(tau_pz)*va(tau_pz)));
  Taumu=TLorentzVector(vb(tau_px),vb(tau_py),vb(tau_pz),sqrt(1.777*1.777+vb(tau_px)*vb(tau_px)+vb(tau_py)*vb(tau_py)+vb(tau_pz)*vb(tau_pz)));
  Zmass = 91.5;
}

   
bool DiTauConstrainedFitter::Fit(){
   return LagrangeMultipliersFitter::Fit();
}
  

LorentzVectorParticle  
DiTauConstrainedFitter::TauMuStartingPoint(TrackParticle MuTrack,LorentzVectorParticle TauA1, TVector3 PV,TMatrixTSym<double>  PVCov, TVector3 SV, TMatrixTSym<double>  SVCov){
  
  TVector3 TauDir =   PV - SV;
  TVector3 TauDirError(sqrt(SVCov(0,0) + PVCov(0,0)),
		       sqrt(SVCov(1,1) + PVCov(1,1)),
		       sqrt(SVCov(2,2) + PVCov(2,2)));

  TMatrixT<double>    parameters;
  parameters.ResizeTo(5,1);
  TMatrixT<double>    parameterErrors;
  parameterErrors.ResizeTo(5,5);

  TMatrixT<double>    parametersAd;
  parametersAd.ResizeTo(9,1);
  TMatrixTSym<double>    parameterErrorsAd;
  parameterErrorsAd.ResizeTo(9,9);

  TMatrixT<double>    taumudirection;
  taumudirection.ResizeTo(2,1);
  TMatrixTSym<double>    taumudirectionError;
  taumudirectionError.ResizeTo(2,2);
 
  TMatrixT<double>    kinematicparameters;
  kinematicparameters.ResizeTo(5,1);
  TMatrixTSym<double>    kinematicparametererrors;
  kinematicparametererrors.ResizeTo(5,5);
  
  TMatrixT<double>    TauKin;
  TauKin.ResizeTo(3,1);
  
  TMatrixT<double>    TauKinErrorAnalytical;
  TauKinErrorAnalytical.ResizeTo(3,3);
  
  TMatrixT<double>    TauKinErrorNumerical;
  TauKinErrorNumerical.ResizeTo(3,3);

  parameters = ConfigureParameters(MuTrack, EstimatePhiAngle(TauDir,TauDirError));
  parameterErrors= ConfigureParameterErrors(MuTrack, EstimatePhiAngle(TauDir,TauDirError));

  parametersAd = ConfigureInitialAdvancedParameters(MuTrack, PV, TauA1);
  parameterErrorsAd= ConfigureInitialAdvancedParameterErrors(MuTrack, PVCov, TauA1.getCovMatrix());

  taumudirection = EstimateTauDirectionAdvanced(parametersAd);
  taumudirectionError = ErrorMatrixPropagator::PropagateError(&DiTauConstrainedFitter::EstimateTauDirectionAdvanced,ConfigureInitialAdvancedParameters(MuTrack, PV, TauA1),ConfigureInitialAdvancedParameterErrors(MuTrack, PVCov, TauA1.getCovMatrix()));

  kinematicparameters = ConfigureKinematicParameters(taumudirection,TauA1);
  kinematicparametererrors =  ConfigureKinematicParameterErrors(taumudirectionError,TauA1);

  TauKin=EstimateTauKinematic(kinematicparameters);
  TauKinErrorNumerical = ErrorMatrixPropagator::PropagateError(&DiTauConstrainedFitter::EstimateTauKinematic,ConfigureKinematicParameters(taumudirection,TauA1),ConfigureKinematicParameterErrors(taumudirectionError,TauA1));


  TauKinErrorAnalytical=ComputeAngleCovarianceAnalytically(MuTrack,EstimatePhiAngle(TauDir,TauDirError),PV,SV,TauA1);

  /*
  TMatrixT<double> MuonTauA1(7,1);
  TMatrixTSym<double> MuonTauA1Cov(7);

  TMatrixT<double> TauKinNEW(3,1);
  TMatrixTSym<double> TauKinErrorNEW(3);

  MuonTauA1 = ConfigureMuTrackTauA1Parameters(MuTrack, TauA1);
  MuonTauA1Cov = ConfigureMuTrackTauA1Errors(MuTrack, TauA1);

  TauKinNEW = TauMuEstimator(MuonTauA1);
  TauKinErrorNEW = ErrorMatrixPropagator::PropagateError(&DiTauConstrainedFitter::TauMuEstimator, MuonTauA1, MuonTauA1Cov);


  Logger(Logger::Debug) << "TauKinNEW" << std::endl;
  TauKinNEW.Print();

  Logger(Logger::Debug) << "TauKinErrorNEW" << std::endl;
  TauKinErrorNEW.Print();
   */

  if(debug){
    std::cout<<"Tau Kinematic Paramteres 'TauKin' ====>"<<std::endl;
    TauKin.Print();

    std::cout<<"Numericaly  calculated tau kinematic parameter errors 'TauKinErrorNumerical'  ====> "<< std::endl;
    TauKinErrorNumerical.Print();

    std::cout<<"Inpute Kinematic Parameters 'kinematicparameters' ====>"<<std::endl;
    kinematicparameters.Print();

    std::cout<<"Input Kinematic Paramter Errors  'kinematicparametererrors' ====>"<< std::endl;
    kinematicparametererrors.Print();

    std::cout<<"TauMu direction paramters 'taumudirection'  ====>  "<< std::endl;
    taumudirection.Print();

    std::cout<<"TauMu direction paramter errors  'taumudirectionError'  ====>  "<< std::endl;
    taumudirectionError.Print();

    std::cout<<"Analitically  calculated tau kinematic parameter errors 'TauKinErrorNumerical'  ====>  "<< std::endl;
    TauKinErrorAnalytical.Print();
    
  }

  TMatrixT<double>    par(LorentzVectorParticle::NLorentzandVertexPar,1);
  TMatrixTSym<double> Cov(LorentzVectorParticle::NLorentzandVertexPar);

  par(LorentzVectorParticle::vx,0)=0; // fill zero vertex for now
  par(LorentzVectorParticle::vy,0)=0;
  par(LorentzVectorParticle::vz,0)=0;
  par(LorentzVectorParticle::px,0)=TauKin(0,0);
  par(LorentzVectorParticle::py,0)=TauKin(1,0);
  par(LorentzVectorParticle::pz,0)=TauKin(2,0);
  par(LorentzVectorParticle::m,0) =1.777;
  
   for(int i=0; i<LorentzVectorParticle::NVertex; i++){
     for(int j=0; j<LorentzVectorParticle::NVertex; j++){
       if(AnalyticalCovariance){Cov(i+3,j+3)=TauKinErrorAnalytical(i,j);}
       else{Cov(i+3,j+3)=TauKinErrorNumerical(i,j);}
     }
   }

  return LorentzVectorParticle(par,Cov,PDGInfo::tau_minus,0,0);
}

LorentzVectorParticle
DiTauConstrainedFitter::TauMuStartingPointwithFullRecoil(TrackParticle MuTrack,LorentzVectorParticle TauA1, PTObject METminusNeutrino, TVector3 PV, TMatrixTSym<double>  PVCov, TVector3 SV, TMatrixTSym<double>  SVCov){

  TVector3 TauDir =   PV - SV;
  TVector3 TauDirError(sqrt(SVCov(0,0) + PVCov(0,0)),
		       sqrt(SVCov(1,1) + PVCov(1,1)),
		       sqrt(SVCov(2,2) + PVCov(2,2)));

  TMatrixT<double>    kinematicparameters;
  kinematicparameters.ResizeTo(12,1);
  TMatrixTSym<double>    kinematicparametererrors;
  kinematicparametererrors.ResizeTo(12,12);

  TMatrixT<double>    taumuptparameters;
  taumuptparameters.ResizeTo(4,1);
  TMatrixTSym<double>    taumuptparametererrors;
  taumuptparametererrors.ResizeTo(4,4);

  TMatrixT<double>    TauMuPt;
  TauMuPt.ResizeTo(2,1);
  TMatrixTSym<double>    TauMuPterrors;
  TauMuPterrors.ResizeTo(2,2);

  TMatrixT<double>    TauKin;
  TauKin.ResizeTo(3,1);

  TMatrixT<double>    TauKinErrorAnalytical;
  TauKinErrorAnalytical.ResizeTo(3,3);

  TMatrixT<double>    TauKinErrorNumerical;
  TauKinErrorNumerical.ResizeTo(3,3);

  taumuptparameters = ConfigureTauMuPtParameters(MuTrack, METminusNeutrino);
  taumuptparametererrors = ConfigureTauMuPtParameterErrors(MuTrack, METminusNeutrino);

  TauMuPt = EstimateTauPt(taumuptparameters);
  TauMuPterrors = ErrorMatrixPropagator::PropagateError(&DiTauConstrainedFitter::EstimateTauPt,taumuptparameters,taumuptparametererrors);

  kinematicparameters = ConfigureKinematicParametersFullRecoil(MuTrack, PV, TauA1, TauMuPt);
  kinematicparametererrors =  ConfigureKinematicParameterErrorsFullRecoil(MuTrack, PVCov, TauA1, TauMuPterrors);

  TauKin=EstimateTauKinematicFullRecoil(kinematicparameters);
  TauKinErrorNumerical = ErrorMatrixPropagator::PropagateError(&DiTauConstrainedFitter::EstimateTauKinematicFullRecoil,kinematicparameters,kinematicparametererrors);


  TauKinErrorAnalytical=ComputeAngleCovarianceAnalytically(MuTrack,EstimatePhiAngle(TauDir,TauDirError),PV,SV,TauA1);

  TMatrixT<double>    par(LorentzVectorParticle::NLorentzandVertexPar,1);
  TMatrixTSym<double> Cov(LorentzVectorParticle::NLorentzandVertexPar);

  par(LorentzVectorParticle::vx,0)=0; // fill zero vertex for now
  par(LorentzVectorParticle::vy,0)=0;
  par(LorentzVectorParticle::vz,0)=0;
  par(LorentzVectorParticle::px,0)=TauKin(0,0);
  par(LorentzVectorParticle::py,0)=TauKin(1,0);
  par(LorentzVectorParticle::pz,0)=TauKin(2,0);
  par(LorentzVectorParticle::m,0) =1.777;

   for(int i=0; i<LorentzVectorParticle::NVertex; i++){
     for(int j=0; j<LorentzVectorParticle::NVertex; j++){
       if(AnalyticalCovariance){Cov(i+3,j+3)=TauKinErrorAnalytical(i,j);}
       else{Cov(i+3,j+3)=TauKinErrorNumerical(i,j);}
     }
   }

  return LorentzVectorParticle(par,Cov,PDGInfo::tau_minus,0,0);
}

TMatrixT<double> 
DiTauConstrainedFitter::EstimateTauDirectionAdvanced(TMatrixT<double> &inpar){
  TMatrixT<double>    outpar(2,1);

  double dxy   =inpar(0,0);
  double phi0  =inpar(1,0);
  double lam   =inpar(2,0);
  double dz    =inpar(3,0);

  TVector3 PV(inpar(4,0),inpar(5,0),inpar(6,0));
  TVector2 TauMuPt(-inpar(7,0), -inpar(8,0));
  
  double tanphi_tau = TauMuPt.Y()/TauMuPt.X();
  double a = tanphi_tau;
  double b = PV.Y() - a*PV.X();

  double tmu  = (dxy*(cos(phi0) + a*sin(phi0)) - b) / (a*cos(phi0) - sin(phi0)); //projection onto XY plane

  double xdoc = -dxy*sin(phi0) + tmu*cos(phi0);
  double ydoc = dxy*cos(phi0) + tmu*sin(phi0);
  double zdoc = dz + tmu*tan(lam);

  TVector3 PointGuess(xdoc, ydoc, zdoc);
  TVector3 TauMuDir = PointGuess - PV;

  outpar(0,0) = TauMuDir.Theta();
  outpar(1,0) = atan2(TauMuPt.Y(),TauMuPt.X());
 
  Logger(Logger::Debug) << "PV: " << inpar(4,0) << ", " << inpar(5,0) << ", " << inpar(6,0) << std::endl;
  Logger(Logger::Debug) << "dxy, phi0, lam, dz: " << dxy << ", " << phi0 << ", " << lam << ", " << dz << std::endl;
  Logger(Logger::Debug) << "PointGuess: " << xdoc << ", " << ydoc << ", " << zdoc << std::endl;
  Logger(Logger::Debug) << "DirGuess: " << TauMuDir.X() << ", " << TauMuDir.Y() << ", " << TauMuDir.Z() << std::endl;
  Logger(Logger::Debug) << "a, b, t: " << a << ", " << b << ", " << tmu <<  std::endl;
  Logger(Logger::Debug) << "TauH X, Y: " << inpar(7,0) << ", " << inpar(8,0) << std::endl;
  Logger(Logger::Debug) << "TauMuDir.Phi(): " << TauMuDir.Phi() << " TauH phi: " << atan2(inpar(8,0), inpar(7,0)) << std::endl;
  Logger(Logger::Debug) << "dPhi: " << TauMuDir.Phi() - atan2(inpar(8,0), inpar(7,0)) << "/pi= " << (TauMuDir.Phi() - atan2(inpar(8,0), inpar(7,0)))/TMath::Pi() << std::endl;

  return outpar; 
}


TMatrixT<double> 
DiTauConstrainedFitter::EstimateTauKinematic(TMatrixT<double> &inpar){
  TMatrixT<double>    outpar(3,1);
  TLorentzVector TauA1p4(inpar(2,0),
			 inpar(3,0),
			 inpar(4,0),
			 sqrt(pow(1.777,2) + pow(inpar(2,0),2) + pow(inpar(3,0),2) + pow(inpar(4,0),2)  ) );
  double ZMassR = 91.5;
  double cosTauTau = cos(TauA1p4.Theta() + inpar(0,0));
  double TauMuP =  MassConstraint_*MassConstraint_/2/(1-cosTauTau)/TauA1p4.P(); // neglecting tau mass w.r.t Z
  
  TVector3 TauMuDir(cos(inpar(1,0))*sin(inpar(0,0)), sin(inpar(1,0))*sin(inpar(0,0)), cos(inpar(0,0)));
  double TauMuPNew = MassConstraint_*MassConstraint_/2/(1-cos(TauA1p4.Angle(TauMuDir)))/TauA1p4.P();

  //Logger(Logger::Debug) << "TauMuP: " << TauMuP << std::endl;

  outpar(0,0) = TauMuP*cos(inpar(1,0))*sin(inpar(0,0));
  outpar(1,0) = TauMuP*sin(inpar(1,0))*sin(inpar(0,0));
  outpar(2,0) = TauMuP*cos(inpar(0,0));

  Logger(Logger::Debug) << "TauDir.Phi(): " << TauA1p4.Phi() << " TauMuDirNEW2.Phi(): " << inpar(1,0) << std::endl;
  Logger(Logger::Debug) << "TauA1 p3: " << TauA1p4.X() << ", " << TauA1p4.Y() << ", " << TauA1p4.Z() << std::endl;
  Logger(Logger::Debug) << "TauMu p3: " << outpar(0,0) << ", " << outpar(1,0) << ", " << outpar(2,0) << std::endl;


  return outpar; 
}


TMatrixT<double> 
DiTauConstrainedFitter::ConfigureParameters(TrackParticle MuTrack, std::pair<double, double> phiAngle){
   TMatrixT<double>    outpar(5,1);
   outpar(0,0) = MuTrack.Parameter(TrackParticle::dxy);
   outpar(1,0) = MuTrack.Parameter(TrackParticle::phi);
   outpar(2,0) = MuTrack.Parameter(TrackParticle::lambda);
   outpar(3,0) = MuTrack.Parameter(TrackParticle::dz);
   outpar(4,0) = phiAngle.first;

   return outpar;
 }

TMatrixT<double> 
DiTauConstrainedFitter::ConfigureParameterErrors(TrackParticle MuTrack, std::pair<double, double> phiAngle){
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
DiTauConstrainedFitter::ConfigureInitialAdvancedParameters(TrackParticle MuTrack,  TVector3 PV, LorentzVectorParticle TauA1){
  TMatrixT<double>    outpar;
  outpar.ResizeTo(9,1);
   outpar(0,0)  = MuTrack.Parameter(TrackParticle::dxy);
   outpar(1,0)  = MuTrack.Parameter(TrackParticle::phi);
   outpar(2,0)  = MuTrack.Parameter(TrackParticle::lambda);
   outpar(3,0)  = MuTrack.Parameter(TrackParticle::dz);
   outpar(4,0)  = PV.X();
   outpar(5,0)  = PV.Y();
   outpar(6,0)  = PV.Z();
   outpar(7,0) = TauA1.LV().X();
   outpar(8,0) = TauA1.LV().Y();

   return outpar;
 }

TMatrixTSym<double> 
DiTauConstrainedFitter::ConfigureInitialAdvancedParameterErrors(TrackParticle MuTrack, TMatrixTSym<double> PVCov, TMatrixTSym<double> TauA1Cov){
  TMatrixTSym<double>  Cov;
    Cov.ResizeTo(9,9);
 
    Cov(0,0) = MuTrack.Covariance(TrackParticle::dxy,TrackParticle::dxy);
    Cov(0,1) = MuTrack.Covariance(TrackParticle::dxy,TrackParticle::phi);
    Cov(0,2) = MuTrack.Covariance(TrackParticle::dxy,TrackParticle::dz);
    Cov(0,3) = MuTrack.Covariance(TrackParticle::dxy,TrackParticle::lambda);

    Cov(1,0) = MuTrack.Covariance(TrackParticle::phi,TrackParticle::dxy);
    Cov(1,1) = MuTrack.Covariance(TrackParticle::phi,TrackParticle::phi);
    Cov(1,2) = MuTrack.Covariance(TrackParticle::phi,TrackParticle::dz);
    Cov(1,3) = MuTrack.Covariance(TrackParticle::phi,TrackParticle::lambda);

    Cov(2,0) = MuTrack.Covariance(TrackParticle::dz,TrackParticle::dxy);
    Cov(2,1) = MuTrack.Covariance(TrackParticle::dz,TrackParticle::phi);
    Cov(2,2) = MuTrack.Covariance(TrackParticle::dz,TrackParticle::dz);
    Cov(2,3) = MuTrack.Covariance(TrackParticle::dz,TrackParticle::lambda);

    Cov(3,0) = MuTrack.Covariance(TrackParticle::lambda,TrackParticle::dxy);
    Cov(3,1) = MuTrack.Covariance(TrackParticle::lambda,TrackParticle::phi);
    Cov(3,2) = MuTrack.Covariance(TrackParticle::lambda,TrackParticle::dz);
    Cov(3,3) = MuTrack.Covariance(TrackParticle::lambda,TrackParticle::lambda);
    for(int i=0; i<3; i++){
      for(int j=0; j<3; j++){
	Cov(i+4,j+4)=PVCov(i,j);
      }
    }
    for(int i=0; i<2; i++){
      for(int j=0; j<2; j++){
	Cov(i+7,j+7)=TauA1Cov(i+4,j+4);
      }
    }
    if(debug){
      std::cout<<"  Input paramters covariance matrix "<<std::endl;
      for(int str =0; str < Cov.GetNrows(); str++){
	for(int kol =0; kol < Cov.GetNcols(); kol++){
	  std::cout<<"  "<< Cov(str,kol)<<"  ";
	}    
	std::cout<<std::endl;
      }
    }

    return Cov;
 }

TMatrixT<double> 
DiTauConstrainedFitter::ConfigureKinematicParameters(TMatrixT<double>  TauMuDir, LorentzVectorParticle TauA1){
  TMatrixT<double>    outpar;
  outpar.ResizeTo(5,1);
  outpar(0,0)  = TauMuDir(0,0);
  outpar(1,0)  = TauMuDir(1,0);
  outpar(2,0)  = TauA1.LV().Px();
  outpar(3,0)  = TauA1.LV().Py();
  outpar(4,0)  = TauA1.LV().Pz();


  if(debug){
    std::cout<<"EstimateTauKinematic  TauA1  "<<TauA1.LV().Px()<<"  "<<TauA1.LV().Py()<<"  "<<TauA1.LV().Pz()<<std::endl;
    std::cout<<"ConfigureKinematicParameters:  outpar(0,0)  "<<outpar(0,0) <<std::endl;
    std::cout<<"ConfigureKinematicParameters:  outpar(1,0)  "<<outpar(1,0) <<std::endl;
    std::cout<<"ConfigureKinematicParameters:  outpar(2,0)  "<<outpar(2,0) <<std::endl;
    std::cout<<"ConfigureKinematicParameters:  outpar(3,0)  "<<outpar(3,0) <<std::endl;
    std::cout<<"ConfigureKinematicParameters:  outpar(4,0)  "<<outpar(4,0) <<std::endl;
  }
  return outpar;
 }

TMatrixTSym<double> 
DiTauConstrainedFitter::ConfigureKinematicParameterErrors(TMatrixTSym<double>  TauMuDirError, LorentzVectorParticle TauA1){
  TMatrixTSym<double>    outpar;
  outpar.ResizeTo(5,5);
  
  for(int i=0; i<5; i++){
    for(int j=0; j<5; j++){
      
      if(i<2 && j < 2){
	outpar(i,j)=TauMuDirError(i,j);
      }else{
	outpar(i,j)=TauA1.Covariance(i+1,j+1);
      }
    }
  }
 
  return outpar;
}

std::pair<double, double> 
DiTauConstrainedFitter::EstimatePhiAngle( TVector3 dir, TVector3 dirE){
    std::pair<double, double> outpar;
  double phi = atan2(dir.Y(),dir.X());
  double deltaphi = sqrt(dirE.Y()*dirE.Y()/dir.X()/dir.X()+
 			 dirE.X()*dirE.X()*dir.Y()*dir.Y()/dir.X()/dir.X()/dir.X()/dir.X() )*cos(phi)*cos(phi);
  outpar = std::make_pair(phi, deltaphi);
  return outpar;
}

TMatrixT<double> DiTauConstrainedFitter::ConfigureMuTrackTauA1Parameters(TrackParticle Muon, LorentzVectorParticle TauA1){
  TMatrixT<double> outpar(7,1);
  outpar(0,0) = Muon.Parameter(TrackParticle::lambda);
  outpar(1,0) = Muon.Parameter(TrackParticle::phi);
  outpar(2,0) = Muon.Parameter(TrackParticle::dxy);
  outpar(3,0) = Muon.Parameter(TrackParticle::dz);
  outpar(4,0) = TauA1.Parameter(LorentzVectorParticle::px);
  outpar(5,0) = TauA1.Parameter(LorentzVectorParticle::py);
  outpar(6,0) = TauA1.Parameter(LorentzVectorParticle::pz);
  return outpar;
}
TMatrixTSym<double> DiTauConstrainedFitter::ConfigureMuTrackTauA1Errors(TrackParticle Muon, LorentzVectorParticle TauA1){
  TMatrixTSym<double> outCov(7);
  for(int i=0; i<outCov.GetNrows(); i++){
	for(int j=0; j<outCov.GetNcols(); j++){
	  if(i < 4 && j < 4){ 						 // only need 4 Muon track parameters (straight line approx., kappa is not needed)
		outCov(i,j) = Muon.Covariance(i+1, j+1); // +1 to exclude kappa
	  }
	  else{
		outCov(i,j) = TauA1.Covariance(i-1, j-1); // Parameters px,py,pz are 3,4,5 in enum
	  }
	}
  }
  return outCov;
}

//////////////////////////////////////////////////////////////////////
// Analytical calculation of TauMu Covariance; Has to be rewritten....
TMatrixT<double> 
DiTauConstrainedFitter::ComputeAngleCovarianceAnalytically(TrackParticle MuTrack, std::pair<double, double> phiAngle,  TVector3 PV, TVector3 SV, LorentzVectorParticle  TauA1){
   TMatrixT<double>    outpar;
   outpar.ResizeTo(3,3);
   double dxy   =MuTrack.Parameter(TrackParticle::dxy);
   double kappa =MuTrack.Parameter(TrackParticle::kappa);
   double phi0  =MuTrack.Parameter(TrackParticle::phi);
   double lam   =MuTrack.Parameter(TrackParticle::lambda);
   double dz    =MuTrack.Parameter(TrackParticle::dz);
   double c     =MuTrack.Charge();
   TVector3 A1SV = -SV + PV;
   
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


   double br = tan(phi0)*cos(phiAnot) - sin(phiAnot);
   double Zc = ZNeu;

   double MinusSintheta = -sqrt(1  - Zc*Zc/(r*r + Zc*Zc) );

   double cosTheta = Zc/sqrt(r*r + Zc*Zc) ;
   double sinTheta = -MinusSintheta;
   if(fabs(phiAnot - TauA1.LV().Phi()) < 0.05) phiAnot = phiAnot - TMath::Pi();
   //-derivaticves
   double drdphi =2*dxy*sin(phi0)*(cos(phiAnot) + tan(phi0)*sin(phiAnot))/ (tan(phi0)*cos(phiAnot) - sin(phiAnot))/(tan(phi0)*cos(phiAnot) - sin(phiAnot));
   double drzdd = 2*sin(phi0)/br/Zc  - 2*dxy*sin(phi0)*lam/br/Zc/Zc;
   double drzdphi0 = (2*dxy*cos(phi0)/(tan(phi0)*cos(phiAnot) - sin(phiAnot))   - 2*dxy*sin(phi0)*cos(phiAnot)/cos(phi0)/cos(phi0)/(tan(phi0)*cos(phiAnot) - sin(phiAnot))/(tan(phi0)*cos(phiAnot) - sin(phiAnot)))/Zc;
   double drzdlam = -2*dxy*sin(phi0)*(r+dxy)/br/Zc/Zc;
   double drzdz0 = -2*dxy*sin(phi0)/br/Zc;
   double drzdphi = drdphi*(1 - r*lam/Zc)/Zc;//2*dxy*sin(phi0)*(cos(phiAnot) + tan(phi0)*sin(phiAnot))/ (tan(phi0)*cos(phiAnot) - sin(phiAnot))/(tan(phi0)*cos(phiAnot) - sin(phiAnot));
   double dcosThetadrz = -r/Zc/sqrt(pow(1 + r*r/Zc/Zc,3));
   double dThetadPhi = -dcosThetadrz*drzdphi/sinTheta;
   //-derivaticves

   double cosTauTau2 = cos(TauA1.LV().Theta() + acos(cosTheta));
   double sinTauTau2 = sin(TauA1.LV().Theta() + acos(cosTheta));

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
   HelixCov(4,4) = phiAngle.second;//sqrt(TauDirError.Y()*TauDirError.Y()/A1SV.X()/A1SV.X()   + TauDirError.X()*TauDirError.X()*A1SV.Y()*A1SV.Y()/A1SV.X()/A1SV.X()/A1SV.X()/A1SV.X() ) *cos(phiAnot)*cos(phiAnot);

   TMatrixT<double>  DerivativesHelixToAngles;
   DerivativesHelixToAngles.ResizeTo(2,5);

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

   double ZMassR = 91.5;
   double TauA1deltaP = sqrt(   (pow(TauA1.LV().Px(),2)*fabs(TauA1.Covariance(3,3))   +  pow(TauA1.LV().Py(),2)*fabs(TauA1.Covariance(4,4))  +  pow(TauA1.LV().Pz(),2)*fabs(TauA1.Covariance(5,5))    )/TauA1.LV().P()/TauA1.LV().P()  )  ;
   double TauMudeltaP_2 = TauA1deltaP*pow(MassConstraint_/2/TauA1.LV().P(),2)/(1-cosTauTau2);

   double TauMuP_2 =  MassConstraint_*MassConstraint_/2/(1-cosTauTau2)/TauA1.LV().P();


   TMatrixT<double> CovPPhiThetaFrame;
   CovPPhiThetaFrame.ResizeTo(3,3);
   //---------- dp dphi dtheta
   CovPPhiThetaFrame(0,0) = TauMudeltaP_2;
   CovPPhiThetaFrame(0,1) = 0;
   CovPPhiThetaFrame(0,2) = 0;
   
   CovPPhiThetaFrame(1,0) = 0;
   CovPPhiThetaFrame(1,1) = CovAngleFrame(0,0);
   CovPPhiThetaFrame(1,2) = CovAngleFrame(0,1);

   CovPPhiThetaFrame(2,0) = 0;
   CovPPhiThetaFrame(2,1) = CovAngleFrame(1,0);
   CovPPhiThetaFrame(2,2) = CovAngleFrame(1,1);


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


   outpar = CovPxPyPzFrame;
   return outpar;
 }


LorentzVectorParticle 
DiTauConstrainedFitter::GetTauMuEstimate(){
  return particles_.at(1);
}

TVector3 DiTauConstrainedFitter::TauMuPtBalanceEstimator(TMatrixT<double> Muon, TVector3 PV, TVector3 SV){
	double dxy 		= Muon(0,0);	//= Muon.Parameter(TrackParticle::dxy);
	double phi0 	= Muon(1,0);	//= Muon.Parameter(TrackParticle::phi);
	double lambda 	= Muon(2,0);   	//= Muon.Parameter(TrackParticle::lambda);
	double dz 		= Muon(3,0);  //= Muon.Parameter(TrackParticle::dz);
	double cosphi0 	= cos(phi0);  double sinphi0 = sin(phi0);
	double coslambda= cos(lambda);double sinlambda = sin(lambda);
	double tanlambda= tan(lambda);

	TVector3 SVPV = SV - PV;
	TVector3 Muon0 = TVector3(-sinphi0*dxy,cosphi0*dxy,dz);
	TVector3 Taumu0 = TVector3(PV);
	TVector3 Delta0 = TVector3(Muon0 - Taumu0);

	double phitau = atan2(-SVPV.Y(),-SVPV.X());
	double cosphitau	= cos(phitau);	double sinphitau	= sin(phitau);
	//double cosphitau 	= -SVPV.X()/SVPV.Perp(); 		double sinphitau	= -SVPV.Y()/SVPV.Perp();

	double tmu = (Delta0.X()/cosphitau - Delta0.Y()/sinphitau)/(cosphi0/cosphitau - sinphi0/sinphitau);
	TVector3 MuonDir = TVector3(tmu*cosphi0, tmu*sinphi0, tmu*tanlambda);
	TVector3 Intersection = MuonDir + Muon0;
	TVector3 TauMuDirFromIntersection = Intersection - Taumu0;

	double distnew = Distance(Muon0, Taumu0, MuonDir, TauMuDirFromIntersection);

	TVector3 TauMuDir = TVector3();

	//(distplus < distminus) ? TauMuDir = TauMuDirPlus : TauMuDir = TauMuDirMinus;
	TauMuDir = TauMuDirFromIntersection;

	//TauMuDir.SetMag(1.);
	return TauMuDir;
}
double DiTauConstrainedFitter::Distance(TVector3 Location1, TVector3 Location2, TVector3 DirectionVector1, TVector3 DirectionVector2){
       TVector3 vecproduct = DirectionVector1.Cross(DirectionVector2);
       TVector3 relvec = TVector3(Location2 - Location1);
       double Mag = abs(vecproduct.Mag());
       double distance = (Mag == 0) ? -1. : abs(vecproduct.Dot(relvec))/Mag;
       return distance;
}
/*TMatrixT<double> DiTauConstrainedFitter::TauMuEstimator(TMatrixT<double> &inpar){
	double lambda = inpar(0,0);   	//= Muon.Parameter(TrackParticle::lambda);
	double phi0   = inpar(1,0);	    //= Muon.Parameter(TrackParticle::phi);
	double dxy    = inpar(2,0);	//= Muon.Parameter(TrackParticle::dxy);
	double dz     = inpar(3,0);  //= Muon.Parameter(TrackParticle::dz);
	double cosphi0 		= cos(phi0); 				double sinphi0 		= sin(phi0);

	TVector3 Tauh = TVector3(inpar(4,0),inpar(5,0),inpar(6,0));
	TVector3 MuonDir = TVector3(cosphi0*cos(lambda), sinphi0*cos(lambda),sin(lambda));
	double P_TauMu = pow(MassConstraint_, 2.)/2./Tauh.Mag()/(1 - cos(Tauh.Angle(MuonDir)));

	double cosPhiR = cos(static_ResRecoil_); double sinPhiR = sin(static_ResRecoil_);
	double Pt_TauMu = -(Tauh.X()/cosPhiR - Tauh.Y()/sinPhiR)/(cosphi0/cosPhiR - sinphi0/sinPhiR);
	double Pz_TauMu = sqrt(pow(P_TauMu, 2.) - pow(Pt_TauMu, 2.));

	TMatrixT<double>	outpar(3,1);
	//outpar(0,0) = Pt_TauMu*cosphi0;
	//outpar(1,0) = Pt_TauMu*sinphi0;
	//outpar(2,0) = Pz_TauMu;

	outpar(0,0) = P_TauMu*MuonDir.X();
	outpar(1,0) = P_TauMu*MuonDir.Y();
	outpar(2,0) = P_TauMu*MuonDir.Z();

	return outpar;
}
*/

TMatrixT<double> DiTauConstrainedFitter::ConfigureKinematicParametersFullRecoil(TrackParticle MuTrack, TVector3 PV, LorentzVectorParticle TauA1, TMatrixT<double> TauMuPt){
  TMatrixT<double> outpar;
  outpar.ResizeTo(12,1);

  outpar(0,0)  = MuTrack.Parameter(TrackParticle::dxy);
  outpar(1,0)  = MuTrack.Parameter(TrackParticle::phi);
  outpar(2,0)  = MuTrack.Parameter(TrackParticle::lambda);
  outpar(3,0)  = MuTrack.Parameter(TrackParticle::dz);
  outpar(4,0)  = PV.X();
  outpar(5,0)  = PV.Y();
  outpar(6,0)  = PV.Z();
  outpar(7,0)  = TauA1.LV().X();
  outpar(8,0)  = TauA1.LV().Y();
  outpar(9,0)  = TauA1.LV().Z();
  outpar(10,0) = TauMuPt(0,0);
  outpar(11,0) = TauMuPt(1,0);

  return outpar;
}

TMatrixTSym<double> DiTauConstrainedFitter::ConfigureKinematicParameterErrorsFullRecoil(TrackParticle MuTrack, TMatrixTSym<double> PVCov, LorentzVectorParticle TauA1, TMatrixTSym<double> TauMuPtCov){
  TMatrixTSym<double>  Cov;
    Cov.ResizeTo(12,12);

    Cov(0,0) = MuTrack.Covariance(TrackParticle::dxy,TrackParticle::dxy);
    Cov(0,1) = MuTrack.Covariance(TrackParticle::dxy,TrackParticle::phi);
    Cov(0,2) = MuTrack.Covariance(TrackParticle::dxy,TrackParticle::dz);
    Cov(0,3) = MuTrack.Covariance(TrackParticle::dxy,TrackParticle::lambda);

    Cov(1,0) = MuTrack.Covariance(TrackParticle::phi,TrackParticle::dxy);
    Cov(1,1) = MuTrack.Covariance(TrackParticle::phi,TrackParticle::phi);
    Cov(1,2) = MuTrack.Covariance(TrackParticle::phi,TrackParticle::dz);
    Cov(1,3) = MuTrack.Covariance(TrackParticle::phi,TrackParticle::lambda);

    Cov(2,0) = MuTrack.Covariance(TrackParticle::dz,TrackParticle::dxy);
    Cov(2,1) = MuTrack.Covariance(TrackParticle::dz,TrackParticle::phi);
    Cov(2,2) = MuTrack.Covariance(TrackParticle::dz,TrackParticle::dz);
    Cov(2,3) = MuTrack.Covariance(TrackParticle::dz,TrackParticle::lambda);

    Cov(3,0) = MuTrack.Covariance(TrackParticle::lambda,TrackParticle::dxy);
    Cov(3,1) = MuTrack.Covariance(TrackParticle::lambda,TrackParticle::phi);
    Cov(3,2) = MuTrack.Covariance(TrackParticle::lambda,TrackParticle::dz);
    Cov(3,3) = MuTrack.Covariance(TrackParticle::lambda,TrackParticle::lambda);

    for(int i=0; i<3; i++){
      for(int j=0; j<3; j++){
        Cov(i+4,j+4)=PVCov(i,j);
        Cov(i+7,j+7)=TauA1.Covariance(i+4,j+4);
      }
    }
    for(int i=0; i<2; i++){
      for(int j=0; j<2; j++){
        Cov(i+10,j+10)=TauMuPtCov(i,j);
      }
    }
    return Cov;
}

TMatrixT<double> DiTauConstrainedFitter::ConfigureTauMuPtParameters(TrackParticle Muon, PTObject METminusNeutrino){
  TMatrixT<double>  outpar;
  outpar.ResizeTo(4,1);

  outpar(0,0) = Muon.Parameter(TrackParticle::kappa)/Muon.BField();
  outpar(1,0) = Muon.Parameter(TrackParticle::phi);
  outpar(2,0) = METminusNeutrino.X();
  outpar(3,0) = METminusNeutrino.Y();

  return outpar;
}

TMatrixTSym<double> DiTauConstrainedFitter::ConfigureTauMuPtParameterErrors(TrackParticle Muon, PTObject METminusNeutrino){
  TMatrixTSym<double>  outpar;
  outpar.ResizeTo(4,4);

  outpar(0,0) = Muon.Covariance(TrackParticle::kappa, TrackParticle::kappa)/Muon.BField()/Muon.BField();
  outpar(1,0) = Muon.Covariance(TrackParticle::kappa, TrackParticle::phi)/Muon.BField();
  outpar(0,1) = Muon.Covariance(TrackParticle::phi, TrackParticle::kappa)/Muon.BField();
  outpar(1,1) = Muon.Covariance(TrackParticle::phi, TrackParticle::phi);

  for(unsigned i=0; i<2; i++){
	for(unsigned j=0; j<2; j++){
	  outpar(i+2,j+2) = METminusNeutrino.Cov()(i,j);
	}
  }
  return outpar;
}

TMatrixT<double> DiTauConstrainedFitter::EstimateTauPt(TMatrixT<double> &inpar){
  TMatrixT<double>  outpar;
  outpar.ResizeTo(2,1);

  double MetX = inpar(2,0);
  double MetY = inpar(3,0);

  double kappaoveralpha = inpar(0,0);
  double phi0 = inpar(1,0);

  outpar(0,0) = MetX + fabs(1/kappaoveralpha)*cos(phi0);
  outpar(1,0) = MetY + fabs(1/kappaoveralpha)*sin(phi0);

  return outpar;
}

TMatrixT<double> DiTauConstrainedFitter::EstimateTauKinematicFullRecoil(TMatrixT<double> &inpar){
  TMatrixT<double> outpar;
  outpar.ResizeTo(3,1);

  double dxy   =inpar(0,0);
  double phi0  =inpar(1,0);
  double lam   =inpar(2,0);
  double dz    =inpar(3,0);

  TVector3 PV(inpar(4,0),inpar(5,0),inpar(6,0));
  TVector3 P_Tauh(inpar(7,0), inpar(8,0), inpar(9,0));
  TVector2 TauMuPt(inpar(10,0),inpar(11,0));

  double tanphi_tau = TauMuPt.Y()/TauMuPt.X();
  double a = tanphi_tau;
  double b = PV.Y() - a*PV.X();

  double tmu  = (dxy*(cos(phi0) + a*sin(phi0)) - b) / (a*cos(phi0) - sin(phi0)); //projection onto XY plane

  double xdoc = -dxy*sin(phi0) + tmu*cos(phi0);
  double ydoc = dxy*cos(phi0) + tmu*sin(phi0);
  double zdoc = dz + tmu*tan(lam);

  TVector3 PointGuess(xdoc, ydoc, zdoc);
  TVector3 TauMuDir = PointGuess - PV;

  double P_TauMu = pow(MassConstraint_, 2.)/2./P_Tauh.Mag()/(1 - cos(P_Tauh.Angle(TauMuDir)));

  double phitau = atan2(TauMuPt.Y(),TauMuPt.X());
  double thetatau = TauMuDir.Theta();

  outpar(0,0) = P_TauMu*cos(phitau)*sin(thetatau);
  outpar(1,0) = P_TauMu*sin(phitau)*sin(thetatau);
  outpar(2,0) = P_TauMu*cos(thetatau);

  Logger(Logger::Debug) << "PV: " << inpar(4,0) << ", " << inpar(5,0) << ", " << inpar(6,0) << std::endl;
  Logger(Logger::Debug) << "dxy, phi0, lam, dz: " << dxy << ", " << phi0 << ", " << lam << ", " << dz << std::endl;
  Logger(Logger::Debug) << "PointGuess: " << xdoc << ", " << ydoc << ", " << zdoc << std::endl;
  Logger(Logger::Debug) << "DirGuess: " << TauMuDir.X() << ", " << TauMuDir.Y() << ", " << TauMuDir.Z() << std::endl;
  Logger(Logger::Debug) << "a, b, t: " << a << ", " << b << ", " << tmu <<  std::endl;
  Logger(Logger::Debug) << "TauH X, Y: " << inpar(7,0) << ", " << inpar(8,0) << std::endl;
  Logger(Logger::Debug) << "TauMuDir.Phi(): " << TauMuDir.Phi() << " TauH phi: " << atan2(inpar(8,0), inpar(7,0)) << std::endl;
  Logger(Logger::Debug) << "dPhi: " << TauMuDir.Phi() - atan2(inpar(8,0), inpar(7,0)) << "/pi= " << (TauMuDir.Phi() - atan2(inpar(8,0), inpar(7,0)))/TMath::Pi() << std::endl;

  return outpar;
}
