// Code written by Vladimir Cherepanov
// RWTH Aachen March 4 2014 (update 2)

#include "SimpleFits/FitSoftware/interface/DiTauConstrainedFitter.h"
#include "SimpleFits/FitSoftware/interface/PDGInfo.h"
#include <iostream>

double DiTauConstrainedFitter::MassConstraint_ = 91.5;
double DiTauConstrainedFitter::PtConstraint_ = 0;
/*double a(0);
TMatrixT<double> a1(5,1);
TMatrixTSym<double> a2(5,5);
TrackParticle DiTauConstrainedFitter::MuTrack_(a1, a2, (int)a, a, a, a);
*/

DiTauConstrainedFitter::DiTauConstrainedFitter(LorentzVectorParticle TauA1,TrackParticle MuTrack, TVector3 PVertex, TMatrixTSym<double> VertexCov):
  LagrangeMultipliersFitter()
{
  debug = false;
  AnalyticalCovariance = false;

  LorentzVectorParticle  TauMuGuess  = TauMuStartingPoint( MuTrack,TauA1,PVertex, VertexCov, TauA1.Vertex(),TauA1.VertexCov() );

  particles_.push_back(TauA1);
  particles_.push_back(TauMuGuess);

//   std::cout<<"DiTauConstrainedFitter  ZMass    "<< (TauA1.LV() + TauMuGuess.LV()).M()<<"  Pt  "<< (TauA1.LV() + TauMuGuess.LV()).Pt() << "  "<<TauA1.LV().Pt() - TauMuGuess.LV().Pt() <<std::endl;
//   std::cout<<"starting D value     "<< (TauA1.LV() + TauMuGuess.LV()).M() - 91.5 <<"  Pt  "<< TauA1.LV().Px()  +  TauMuGuess.LV().Px() << "  "<<TauA1.LV().Py()  +  TauMuGuess.LV().Py() <<std::endl;
//   std::cout<<"starting D value     "<< (TauA1.LV() + TauMuGuess.LV()).M() - 91.5 + TauA1.LV().Px()  +  TauMuGuess.LV().Px() + TauA1.LV().Py()  +  TauMuGuess.LV().Py() <<std::endl;

  isconfigured=false;

  int size=particles_.size()*3;   // 
  int sizeExpanded = 3; 
  TMatrixT<double>    inpar(size,1);
  TMatrixTSym<double> incov(size);
 
  // Get primary vertex information
  if(VertexCov.GetNrows()!=LorentzVectorParticle::NVertex)return;

  //set input parameters:  TauA1 - TauMu
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
  //std::cout<< "DiTauConstrainedFitter::DiTauConstrainedFitter: print incov" <<std::endl;
  //incov.Print();
 // store expanded par for computation of final par (assumes fit has neglegible impact on a1 correlations with vertex uncertainties)
  exppar.ResizeTo(size,1);
  exppar=ComputeInitalExpPar(inpar);
  expcov.ResizeTo(size,size);
  expcov=ErrorMatrixPropagator::PropagateError(&DiTauConstrainedFitter::ComputeInitalExpPar,inpar,incov);

  //std::cout<< "DiTauConstrainedFitter::DiTauConstrainedFitter: print expcov" <<std::endl;
  //expcov.Print();

  // store linearization point
  TMatrixT<double> PAR_0(3,1);
  par_0.ResizeTo(3);
  cov_0.ResizeTo(3,3);
  PAR_0=ComputeExpParToPar(exppar);
  for(int i=0; i<3;i++)par_0(i)=PAR_0(i,0);
  //std::cout<< "DiTauConstrainedFitter::DiTauConstrainedFitter: print cov_0" <<std::endl;
  //cov_0.Print();
  cov_0=ErrorMatrixPropagator::PropagateError(&DiTauConstrainedFitter::ComputeExpParToPar,exppar,expcov);
  //std::cout<< "DiTauConstrainedFitter::DiTauConstrainedFitter: print cov_0" <<std::endl;
  //cov_0.Print();
  // set up initial point for fit (cov handled in Fit() function)
  par.ResizeTo(3);
  par=par_0;

  isconfigured=true; 
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
  TMatrixT<double> outpar(3,1);
  //  for(int i=0;i<npar;i++){outpar(i,0)=inpar(i,0);}
  outpar(0,0) = inpar(3,0);
  outpar(1,0) = inpar(4,0);
  outpar(2,0) = inpar(5,0);
  return outpar;
}


TMatrixT<double> DiTauConstrainedFitter::ComputeTauMuLorentzVectorPar(TMatrixT<double> &inpar){
  //  start with index 3 to fill only momenta part
  TMatrixT<double> outpar(7,1);

   outpar(LorentzVectorParticle::px,0)=inpar(3,0);
   outpar(LorentzVectorParticle::py,0)=inpar(4,0);
   outpar(LorentzVectorParticle::pz,0)=inpar(5,0);
   outpar(LorentzVectorParticle::m,0)=1.777;

//    std::cout<<"--ComputeTauMuLorentzVectorPar  "<<std::endl; 
//    std::cout<<"--ComputeTauMuLorentzVectorPar  outpar(0,0)"<< outpar(3,0)<<std::endl; 
//    std::cout<<"--ComputeTauMuLorentzVectorPar  outpar(1,0)"<< outpar(4,0)<<std::endl; 
//    std::cout<<"--ComputeTauMuLorentzVectorPar  outpar(2,0)"<< outpar(5,0)<<std::endl; 
//    std::cout<<"--ComputeTauMuLorentzVectorPar  outpar(3,0)"<< outpar(6,0)<<std::endl; 

  return outpar;
}

TMatrixT<double> DiTauConstrainedFitter::ComputeTauA1LorentzVectorPar(TMatrixT<double> &inpar){
  //  start with index 3 to fill only momenta part
  TMatrixT<double> outpar(7,1);
  outpar(LorentzVectorParticle::px,0)=inpar(0,0);
  outpar(LorentzVectorParticle::py,0)=inpar(1,0);
  outpar(LorentzVectorParticle::pz,0)=inpar(2,0);
  outpar(LorentzVectorParticle::m,0) =1.777;
//    std::cout<<"ComputeTauA1LorentzVectorPar  "<<std::endl; 
//    std::cout<<"ComputeTauA1LorentzVectorPar  outpar(0,0)"<< outpar(3,0)<<std::endl; 
//    std::cout<<"ComputeTauA1LorentzVectorPar  outpar(1,0)"<< outpar(4,0)<<std::endl; 
//    std::cout<<"ComputeTauA1LorentzVectorPar  outpar(2,0)"<< outpar(5,0)<<std::endl; 
//    std::cout<<"ComputeTauA1LorentzVectorPar  outpar(3,0)"<< outpar(6,0)<<std::endl; 
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
  // assumes changes to a1 correlation to vertex is small
  //if(par.GetNrows()==3 && cov.GetNrows() && exppar.GetNrows()==npar && expcov.GetNrows()) return;
  for(int i=0; i<3;i++){
    exppar(i+3,0)=par(i);
    //    std::cout<<"  in loop  par "<<par(i)<<std::endl;
    for(int j=0; j<3;j++){expcov(i+3,j+3)=cov(i,j);}
  }
}

bool DiTauConstrainedFitter::Fit(){
   return LagrangeMultipliersFitter::Fit();
}

std::vector<LorentzVectorParticle> DiTauConstrainedFitter::GetReFitDaughters(){
  std::vector<LorentzVectorParticle> refitParticles;

  UpdateExpandedPar();

  double c(0),b(0);
  for(unsigned int i=0;i<particles_.size();i++){c+=particles_.at(i).Charge();b=particles_.at(i).BField();}
  TMatrixT<double> a1=ComputeTauA1LorentzVectorPar(exppar);

  TMatrixTSym<double> a1cov=ErrorMatrixPropagator::PropagateError(&DiTauConstrainedFitter::ComputeTauA1LorentzVectorPar,exppar,expcov);
  refitParticles.push_back(LorentzVectorParticle(a1,a1cov,PDGInfo::tau_plus,c,b));

//     std::cout<<"a1 cross2 ComputeTauA1LorentzVectorPar  "<<std::endl; 
//     std::cout<<"a1 cross2 ComputeTauA1LorentzVectorPar  outpar(0,0)"<< refitParticles.at(0).LV().Px()<<std::endl; 
//     std::cout<<"a1 cross2 ComputeTauA1LorentzVectorPar  outpar(1,0)"<< refitParticles.at(0).LV().Py()<<std::endl; 
//     std::cout<<"a1 cross2 ComputeTauA1LorentzVectorPar  outpar(2,0)"<< refitParticles.at(0).LV().Pz()<<std::endl; 
//     std::cout<<"a1 cross2 ComputeTauA1LorentzVectorPar  outpar(3,0)"<< refitParticles.at(0).LV().M()<<std::endl; 

  TMatrixT<double> nu=ComputeTauMuLorentzVectorPar(exppar);

  nu(0,0)= particles_.at(1).Parameter(LorentzVectorParticle::vx);
  nu(1,0)= particles_.at(1).Parameter(LorentzVectorParticle::vy);
  nu(2,0)= particles_.at(1).Parameter(LorentzVectorParticle::vz);
  TMatrixTSym<double> nucov=ErrorMatrixPropagator::PropagateError(&DiTauConstrainedFitter::ComputeTauMuLorentzVectorPar,exppar,expcov);
  for(int i=0; i<LorentzVectorParticle::NVertex; i++){
    for(int j=0; j<LorentzVectorParticle::NVertex; j++){
      nucov(i,j)=particles_.at(1).VertexCov()(i,j);
      
    }
  }
  
  refitParticles.push_back(LorentzVectorParticle(nu,nucov,PDGInfo::tau_minus,0.0,b));

  return refitParticles; 
}

LorentzVectorParticle DiTauConstrainedFitter::GetMother(){
  UpdateExpandedPar();
  double c(0),b(0);
  for(unsigned int i=0;i<particles_.size();i++){c+=particles_.at(i).Charge();b=particles_.at(i).BField();}
  TMatrixT<double> m=ComputeMotherLorentzVectorPar(exppar);
  TMatrixTSym<double> mcov=ErrorMatrixPropagator::PropagateError(&DiTauConstrainedFitter::ComputeMotherLorentzVectorPar,exppar,expcov);
  return LorentzVectorParticle(m,mcov,PDGInfo::Z0,c,b);
}

TVectorD DiTauConstrainedFitter::Value(TVectorD &v){
  TLorentzVector Taua1,Taumu; 
  double ZMass;
  ConvertParToObjects(v,Taua1,Taumu,ZMass);

  TLorentzVector z=Taua1+Taumu; 
  TVector3 Taumu_p0(par_0(0),par_0(1),par_0(2));

  TVectorD d(2);
  d(0) = z.M() - MassConstraint_;
  d(1) = Taua1.Pt() - Taumu.Pt();


  //TVectorD d(4);
  //d(0) = z.M() - MassConstraint_;
  //d(1) = Taua1.Px() + Taumu.Px();
  //d(2) = Taua1.Py() + Taumu.Py();
  //d(3) = Taumu_p0.CosTheta() - Taumu.CosTheta();

  return d;
}

void DiTauConstrainedFitter::ConvertParToObjects(TVectorD &v,TLorentzVector &Taua1,TLorentzVector &Taumu, double &Zmass){
  Taua1=particles_.at(0).LV();//TLorentzVector(v(taua1_px),v(taua1_py),v(taua1_pz),sqrt(1.777*1.777+v(taua1_px)*v(taua1_px)+v(taua1_py)*v(taua1_py)+v(taua1_pz)*v(taua1_pz)));
  Taumu=TLorentzVector(v(taua1_px),v(taua1_py),v(taua1_pz),sqrt(1.777*1.777+v(taua1_px)*v(taua1_px)+v(taua1_py)*v(taua1_py)+v(taua1_pz)*v(taua1_pz)));
  Zmass = 91.5;
}
  
LorentzVectorParticle DiTauConstrainedFitter::TauMuStartingPoint(TrackParticle MuTrack,LorentzVectorParticle TauA1, TVector3 PV,TMatrixTSym<double>  PVCov, TVector3 SV, TMatrixTSym<double>  SVCov ){
  
  //MuTrack_ = MuTrack;

  //std::cout << MuTrack_.Parameter(0) << std::endl;
  //std::cout << MuTrack.Parameter(0) << std::endl;

  //std::cout << "TauMuStartingPoint: PVCov" << std::endl;
  //PVCov.Print();
  //std::cout << "TauMuStartingPoint: SVCov" << std::endl;
  //SVCov.Print();

  TVector3 TauDir =   PV - SV;
  TVector3 TauDirError(sqrt(SVCov(0,0) + PVCov(0,0)),
		       sqrt(SVCov(1,1) + PVCov(1,1)),
		       sqrt(SVCov(2,2) + PVCov(2,2)));

  TMatrixT<double>    parameters;
  parameters.ResizeTo(5,1);
  TMatrixT<double>    parameterErrors;
  parameterErrors.ResizeTo(5,5);

  TMatrixT<double>    parametersAd;
  parametersAd.ResizeTo(10,1);
  TMatrixTSym<double>    parameterErrorsAd;
  parameterErrorsAd.ResizeTo(10,10);

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
  //std::cout << "TauMuStartingPoint: parameterErrors" << std::endl;
  //parameterErrors.Print();

  parametersAd = ConfigureInitialAdvancedParameters(MuTrack, PV,SV);
  parameterErrorsAd= ConfigureInitialAdvancedParameterErrors(MuTrack, PVCov,SVCov);

  //std::cout << "TauMuStartingPoint: parametersAd" << std::endl;
  //parametersAd.Print();
  //std::cout << "TauMuStartingPoint: parameterErrorsAd" << std::endl;
  //parameterErrorsAd.Print();

  taumudirection = EstimateTauDirectionAdvanced(parametersAd);
  taumudirectionError = ErrorMatrixPropagator::PropagateError(&DiTauConstrainedFitter::EstimateTauDirectionAdvanced,ConfigureInitialAdvancedParameters(MuTrack, PV,SV),ConfigureInitialAdvancedParameterErrors(MuTrack, PVCov,SVCov));

  //std::cout << "TauMuStartingPoint: taumudirection" << std::endl;
  //taumudirection.Print();
  //std::cout << "TauMuStartingPoint: taumudirectionError" << std::endl;
  //taumudirectionError.Print();

  kinematicparameters = ConfigureKinematicParameters(taumudirection,TauA1);
  kinematicparametererrors =  ConfigureKinematicParameterErrors(taumudirectionError,TauA1);

  //std::cout << "TauMuStartingPoint: kinematicparameters" << std::endl;
  //kinematicparameters.Print();
  //std::cout << "TauMuStartingPoint: kinematicparametererrors" << std::endl;
  //kinematicparametererrors.Print();

  TauKin=EstimateTauKinematic(kinematicparameters);
  TauKinErrorNumerical = ErrorMatrixPropagator::PropagateError(&DiTauConstrainedFitter::EstimateTauKinematic,ConfigureKinematicParameters(taumudirection,TauA1),ConfigureKinematicParameterErrors(taumudirectionError,TauA1));

  //std::cout << "TauMuStartingPoint: TauKin" << std::endl;
  //TauKin.Print();
  //std::cout << "TauMuStartingPoint: TauKinErrorNumerical" << std::endl;
  //TauKinErrorNumerical.Print();

  TauKinErrorAnalytical=ComputeAngleCovarianceAnalytically(MuTrack,EstimatePhiAngle(TauDir,TauDirError),PV,SV,TauA1);

  if(debug){
    std::cout<<"Tau Kinematic Paramteres 'TauKin' ====>"<<std::endl;
    for(int str =0; str < TauKin.GetNrows(); str++){
      for(int kol =0; kol < TauKin.GetNcols(); kol++){
	std::cout<<"  "<< TauKin(str,kol)<<"  ";
      }
      std::cout<<std::endl;
    }
    std::cout<<"Numericaly  calculated tau kinematic parameter errors 'TauKinErrorNumerical'  ====> "<< std::endl;
    for(int str =0; str < TauKinErrorNumerical.GetNrows(); str++){
      for(int kol =0; kol < TauKinErrorNumerical.GetNcols(); kol++){
	std::cout<<"  "<< TauKinErrorNumerical(str,kol)<<"  ";
      }
      std::cout<<std::endl;
    }
    std::cout<<"Inpute Kinematic Parameters 'kinematicparameters' ====>"<<std::endl;
    for(int str =0; str < kinematicparameters.GetNrows(); str++){
      for(int kol =0; kol < kinematicparameters.GetNcols(); kol++){
	std::cout<<"  "<< kinematicparameters(str,kol)<<"  ";
      }
      std::cout<<std::endl;
    }
    std::cout<<"Input Kinematic Paramter Errors  'kinematicparametererrors' ====>"<< std::endl;
    for(int str =0; str < kinematicparametererrors.GetNrows(); str++){
      for(int kol =0; kol < kinematicparametererrors.GetNcols(); kol++){
	std::cout<<"  "<< kinematicparametererrors(str,kol)<<"  ";
      }
      std::cout<<std::endl;
    }
    std::cout<<"TauMu direction paramters 'taumudirection'  ====>  "<< std::endl;
    for(int str =0; str < taumudirection.GetNrows(); str++){
      for(int kol =0; kol < taumudirection.GetNcols(); kol++){
	std::cout<<"  "<< taumudirection(str,kol)<<"  ";
      }
      std::cout<<std::endl;
    }
    std::cout<<"TauMu direction paramter errors  'taumudirectionError'  ====>  "<< std::endl;
    for(int str =0; str < taumudirectionError.GetNrows(); str++){
      for(int kol =0; kol < taumudirectionError.GetNcols(); kol++){
	std::cout<<"  "<< taumudirectionError(str,kol)<<"  ";
      }
      std::cout<<std::endl;
    }
    std::cout<<"Analitically  calculated tau kinematic parameter errors 'TauKinErrorNumerical'  ====>  "<< std::endl;
    for(int str =0; str < TauKinErrorAnalytical.GetNrows(); str++){
      for(int kol =0; kol < TauKinErrorAnalytical.GetNcols(); kol++){
	std::cout<<"  "<< TauKinErrorAnalytical(str,kol)<<"  ";
      }
      std::cout<<std::endl;
    }
  }

  TMatrixT<double>    par(LorentzVectorParticle::NLorentzandVertexPar,10);
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

TMatrixT<double> DiTauConstrainedFitter::EstimateTauDirectionAdvanced(TMatrixT<double> &inpar){
  TMatrixT<double>    outpar(2,1);
  
  double dxy   =inpar(0,0);
  double phi0  =inpar(1,0);
  double lam   =inpar(2,0);
  double dsz   =inpar(3,0);
  //  double phiAnot = inpar(4,0);

  TVector3 PV(inpar(4,0),inpar(5,0),inpar(6,0));
  TVector3 SV(inpar(7,0),inpar(8,0),inpar(9,0)); 
  
  TVector3 TauDir = PV - SV;
  double TauA1OppositePhi = atan2(TauDir.Y(), TauDir.X());

  double xdoc = dxy*sin(phi0) - PV.X(); //dxy*sin(phi0)
  double ydoc = -dxy*cos(phi0) - PV.Y(); //-dxy*cos(phi0)
  double zdoc = dsz; //this is right assuming ref. point is POCA

  double amuon = tan(phi0);
  double bmuon = ydoc  - amuon*xdoc;
  
  double r = sqrt( pow(bmuon/(tan(TauA1OppositePhi) - amuon )  ,2) + pow(bmuon*tan(TauA1OppositePhi)/(tan(TauA1OppositePhi) - amuon) ,2));//(bAnot)/(tan(phiAnot) - tan(phi0))/cos(phiAnot);
  
  double xPoint = r*cos(TauA1OppositePhi);
  double yPoint = r*sin(TauA1OppositePhi);
  
  //    double zPoint = zdoc + r*cos(phi0 +TauA1OppositePhi )*lam;
  double zPoint = (r - dxy)*lam - zdoc;
  
  double cosTheta2 = zPoint/sqrt(r*r + zPoint*zPoint) ;
  TVector3 PointGuess(xPoint, yPoint, zPoint);
  TVector3 TauMuDir = PointGuess - PV;
  
  /*
  //new approach:
  
  double x0 = dxy*sin(phi0);
  double y0 = -dxy*cos(phi0);
  double z0 = dsz;

  double alpha = -1./MuTrack_.BField();
  double radius = alpha/MuTrack_.Parameter(0);

  double x = x0 + radius*sin(phi0) - radius*sin(TauA1OppositePhi + phi0);
  double y = y0 - radius*cos(phi0) + radius*cos(TauA1OppositePhi + phi0);
  double z = z0 - radius*tan(lam)*TauA1OppositePhi;

  TVector3 PointGuessNew(x,y,z);
  TVector3 TauMuDirNew = PointGuessNew - PV;

  */
  outpar(0,0) = acos(cosTheta2);
  outpar(1,0) = TauDir.Phi();
   

  /*
  std::cout<<" TauA1OppositePhi: " << TauA1OppositePhi <<std::endl;
  std::cout<<" radius: " << radius << ", BField: "<< MuTrack_.BField() << ", kappa and charge: " << MuTrack_.Parameter(0) << " ," << MuTrack_.Charge() <<std::endl;
  std::cout<<" x = " << x0 << " + " << radius*sin(phi0) << " - " <<  radius*sin(TauA1OppositePhi + phi0) << std::endl;
  std::cout<<" x0, y0, z0  "<< x0 << ", " << y0 << ", "<< z0 <<std::endl;
  std::cout<<" xPoint, yPoint, zPoint  "<< xPoint << ", " << yPoint << ", "<< zPoint <<std::endl;
  std::cout<<" x, y, z  "<< x << ", " << y << ", "<< z <<std::endl;
  */
//    std::cout<<" TauDir.Phi()   "<< TauDir.Phi() <<"  PointGuess.Phi() "<< TauMuDir.Phi()<<std::endl;
//    std::cout<<" TauMuDir.Theta()   "<< TauMuDir.Theta() <<std::endl;
//    std::cout<<" cosTheta2   "<< cosTheta2 <<std::endl;
//    std::cout<<"  EstimateTauDirectionAdvanced  inpar "<<std::endl;
//    for(int str =0; str < inpar.GetNrows(); str++){
//      for(int kol =0; kol < inpar.GetNcols(); kol++){
//        std::cout<<"  "<< inpar(str,kol)<<"  ";
//      }    
//      std::cout<<std::endl;
//    }
//    std::cout<<"  EstimateTauDirectionAdvanced  outpar(0,0) "<<outpar(0,0)<<std::endl;
//    std::cout<<"  EstimateTauDirectionAdvanced  outpar(1,0) "<<outpar(1,0)<<std::endl;
   

  return outpar; 
}


TMatrixT<double> DiTauConstrainedFitter::EstimateTauKinematic(TMatrixT<double> &inpar){
  TMatrixT<double>    outpar(3,1);
  TLorentzVector TauA1p4(inpar(2,0),
			 inpar(3,0),
			 inpar(4,0),
			 sqrt(pow(1.777,2) + pow(inpar(2,0),2) + pow(inpar(3,0),2) + pow(inpar(4,0),2)  ) );

  outpar(0,0) = TauA1p4.Pt()*cos(inpar(1,0));
  outpar(1,0) = TauA1p4.Pt()*sin(inpar(1,0));
  outpar(2,0) = TauA1p4.Pt()/tan(inpar(0,0));

  /*
  double ZMassR = 91.5;
  double cosTauTau = cos(TauA1p4.Theta() + inpar(0,0));
  double TauMuP =  MassConstraint_*MassConstraint_/2/(1-cosTauTau)/TauA1p4.P(); // neglecting tau mass w.r.t Z
  
  outpar(0,0) = TauMuP*cos(inpar(1,0))*sin(inpar(0,0));
  outpar(1,0) = TauMuP*sin(inpar(1,0))*sin(inpar(0,0));
  outpar(2,0) = TauMuP*cos(inpar(0,0));
    */
  /*
  std::cout<<"  inpar(0,0)   "<<inpar(0,0) <<"  inpar(1,0)   " <<inpar(1,0) <<std::endl;
  std::cout<<"  TauA1p4.Phi()   "<<TauA1p4.Phi() <<"  TauA1p4.P()   " <<TauA1p4.P()<<std::endl;
  std::cout<<"  compare phis2  "<<fabs(TauA1p4.Phi() - inpar(1,0))<<std::endl;

  std::cout<<"  TauMuP   "<<TauMuP <<"  cosTauTau   " <<cosTauTau <<std::endl;

  std::cout<<"  EstimateTauKinematic   outpar  "<<std::endl;
  for(int str =0; str < outpar.GetNrows(); str++){
    for(int kol =0; kol < outpar.GetNcols(); kol++){
      std::cout<<"  "<< outpar(str,kol)<<"  ";
    }    
    std::cout<<std::endl;
  }
  */
  return outpar; 
}

TMatrixT<double> DiTauConstrainedFitter::ConfigureParameters(TrackParticle MuTrack, std::pair<double, double> phiAngle){
   TMatrixT<double>    outpar(5,1);
   outpar(0,0) = MuTrack.Parameter(TrackParticle::dxy);
   outpar(1,0) = MuTrack.Parameter(TrackParticle::phi);
   outpar(2,0) = MuTrack.Parameter(TrackParticle::lambda);
   outpar(3,0) = MuTrack.Parameter(TrackParticle::dz);
   outpar(4,0) = phiAngle.first;

   return outpar;
 }

TMatrixT<double> DiTauConstrainedFitter::ConfigureParameterErrors(TrackParticle MuTrack, std::pair<double, double> phiAngle){
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

TMatrixT<double> DiTauConstrainedFitter::ConfigureInitialAdvancedParameters(TrackParticle MuTrack,  TVector3 PV, TVector3 SV){
  TMatrixT<double>    outpar;
  outpar.ResizeTo(10,1);
   outpar(0,0)  = MuTrack.Parameter(TrackParticle::dxy);
   outpar(1,0)  = MuTrack.Parameter(TrackParticle::phi);
   outpar(2,0)  = MuTrack.Parameter(TrackParticle::lambda);
   outpar(3,0)  = MuTrack.Parameter(TrackParticle::dz);
   outpar(4,0)  = PV.X();
   outpar(5,0)  = PV.Y();
   outpar(6,0)  = PV.Z();
   outpar(7,0)  = SV.X();
   outpar(8,0)  = SV.Y();
   outpar(9,0)  = SV.Z();
   /*
   std::cout<<" ConfigureInitialAdvancedParameters   outpar "<<std::endl;
   for(int str =0; str < outpar.GetNrows(); str++){
     for(int kol =0; kol < outpar.GetNcols(); kol++){
       std::cout<<"  "<< outpar(str,kol)<<"  ";
     }    
     std::cout<<std::endl;
   }
   */
   return outpar;
}

TMatrixTSym<double> DiTauConstrainedFitter::ConfigureInitialAdvancedParameterErrors(TrackParticle MuTrack, TMatrixTSym<double>  PVCov,TMatrixTSym<double>  SVCov){
  TMatrixTSym<double>  Cov;
  Cov.ResizeTo(10,10);
 
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
  	Cov(i+7,j+7)=SVCov(i,j);
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

TMatrixT<double> DiTauConstrainedFitter::ConfigureKinematicParameters(TMatrixT<double>  TauMuDir, LorentzVectorParticle TauA1){
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

TMatrixTSym<double> DiTauConstrainedFitter::ConfigureKinematicParameterErrors(TMatrixTSym<double>  TauMuDirError, LorentzVectorParticle TauA1){
  TMatrixTSym<double>    outpar;
  outpar.ResizeTo(5,5);
  
  for(int i=0; i<5; i++){
    for(int j=0; j<5; j++){
      if(i<2 && j < 2){
    	outpar(i,j)=TauMuDirError(i,j);
      }
      else{
    	outpar(i,j)=TauA1.Covariance(i+1,j+1);
      }
    }
  }
  /*
  std::cout<<"  ConfigureKinematicParameterErrors "<<std::endl;
  for(int str =0; str < outpar.GetNrows(); str++){
    for(int kol =0; kol < outpar.GetNcols(); kol++){
      std::cout<<"  "<< outpar(str,kol)<<"  ";
    }    
    std::cout<<std::endl;
  }
  */
  return outpar;
}

std::pair<double, double> DiTauConstrainedFitter::EstimatePhiAngle( TVector3 dir, TVector3 dirE){
    std::pair<double, double> outpar;
  double phi = atan2(dir.Y(),dir.X());
  double deltaphi = sqrt(dirE.Y()*dirE.Y()/dir.X()/dir.X()+
 			 dirE.X()*dirE.X()*dir.Y()*dir.Y()/dir.X()/dir.X()/dir.X()/dir.X() )*cos(phi)*cos(phi);
  outpar = std::make_pair(phi, deltaphi);
  return outpar;
}

//////////////////////////////////////////////////////////////////////
// Analytical calculation of TauMu Covariance; Has to be rewritten....
TMatrixT<double> DiTauConstrainedFitter::ComputeAngleCovarianceAnalytically(TrackParticle MuTrack, std::pair<double, double> phiAngle,  TVector3 PV, TVector3 SV, LorentzVectorParticle  TauA1){
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


LorentzVectorParticle DiTauConstrainedFitter::GetTauMuEstimate(){
  return particles_.at(1);
}
