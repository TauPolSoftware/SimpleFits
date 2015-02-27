#include "SimpleFits/FitSoftware/interface/Chi2VertexFitter.h"
#include "SimpleFits/FitSoftware/interface/ChiSquareFunctionUpdator.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnStrategy.h"
#include "Minuit2/MnSimplex.h"
#include "Minuit2/CombinedMinimizer.h"
#include "Minuit2/MnMinimize.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnHesse.h"
#include "Minuit2/MnContours.h"
#include "Minuit2/MnPlot.h"
#include "Minuit2/MinosError.h"
#include "Minuit2/ContoursError.h"
#include <iostream>


Chi2VertexFitter::Chi2VertexFitter(std::vector<TrackParticle> &particles_,TVector3 &vguess, double mass, double initwidth, double x, double nsigma_):
  TrackHelixVertexTools(particles_,vguess),
  nsigma(nsigma_)
{
  ConstrainedHelices(par,parcov);
  Helices(trackval,trackcov);
  val.ReSizeTo(trackval.GetNrows()+1,1);
  cov.ReSizeTo(trackval.GetNrows()+1,trackval.GetNrows()+1);
  
  // Map track values to val
  for(int i=0;i<trackval.GetNrows();i++){
    val(i,0)=trackval(i,0);
    for(int j=0;j<trackval.GetNrows();j++)cov(i,j)=trackcov(i,j);
  }
  constraints_.push_back(new HardConstraint("Mass",val.GetNrows()-1,val,cov,dalpha,mass,initwidth,x));
  FitterBase::ReSet();
}

double Chi2VertexFitter::UpdateChisquare(TMatrixT<double> inpar){
  // Update track paramters
  TMatrixT<double> trackvprime=ConvertToHelicesNotation(inpar);
  TMatrixT<double> vprime=val; for(int i=0;i<trackvprime.GetNrows();i++)vprime(i,0)=trackvprime(i,0);

  // compute the mass
  if(includeConstraint) vprime(constraint_.at(mass)->Index(),0)=(ComputeMotherLorentzVectorPar(inpar)(LorentzVectorParticle::m,0));

  // compute chi2
  dalpha=vprime-val;
  TMatrixT<double> dalphaT=dalpha;  dalphaT.T();
  TMatrixT<double> chisquare=dalphaT*(cov_inv*dalpha);
  return chisquare(0,0);
}

bool Chi2VertexFitter::FitwithConstraints(){
  if(isFit==true) return true;// do not refit
  if(!isConfigure) return false; // do not fit if configuration failed
  ChiSquareFunctionUpdator<FitterBase> updator(this);
  ROOT::Minuit2::MnUserParameters MnPar;
  for(int i=0;i<par.GetNrows();i++){
    TString name=FreeParName(i);
    // if not limited (vhigh <= vlow)
    MnPar.Add(name.Data(),par(i,0),sqrt(fabs(parcov(i,i))),par(i,0)-nsigma*sqrt(fabs(parcov(i,i))),par(i,0)+nsigma*sqrt(fabs(parcov(i,i))));
  }
  
  unsigned int max=10;
  int numberofcalls=200+par.GetNrows()*100+par.GetNrows()*par.GetNrows()*5;
  double tolerance(0.01);
  double edmMin(0.001*updator.Up()*tolerance); 
  
  ROOT::Minuit2::MnHesse IntHesse;
  ROOT::Minuit2::MnUserParameterState intState=IntHesse(updator,MnPar);
  ROOT::Minuit2::MnStrategy Strategy(1);
  ROOT::Minuit2::MnMigrad minimize(updator,intState,Strategy);
  //ROOT::Minuit2::MnMinimize minimize(updator,MnPar);
  ROOT::Minuit2::FunctionMinimum min= minimize(numberofcalls,tolerance);
  for(unsigned int i=0;i<=max && min.Edm()>edmMin;i++){
    if(i==max) return false;
    min = minimize(i*numberofcalls,tolerance);
    Log(Log::Warning) << "Minimum not found: interating...." << i << std::endl;
  }
  ROOT::Minuit2::MnHesse Hesse;
  IntHesse(updator,min);
  
  // give return flag based on status
  if(min.IsAboveMaxEdm()){Log(Log::Error) << "Chi2VertexFitter::Fit() Found Vertex that is above EDM " << std::endl; return false;}
  if(!min.IsValid()){
    Log(Log::Error) << "Chi2VertexFitter::Fit(): Failed min.IsValid()" << std::endl; 
    if(!min.HasValidParameters()){Log(Log::Error) << "Chi2VertexFitter::Fit(): Failed min.HasValidParameters()" << std::endl; }
    if(!min.HasValidCovariance()){Log(Log::Error) << "Chi2VertexFitter::Fit(): Failed min.HasValidCovariance()" << std::endl; }
    if(!min.HesseFailed()){Log(Log::Error) << "Chi2VertexFitter::Fit(): Failed min.HesseFailed()" << std::endl; }
    if(!min.HasReachedCallLimit()){Log(Log::Error) << "Chi2VertexFitter::Fit(): Failed min.HasReachedCallLimit()" << std::endl; }
    return false;
  }
  chi2=min.Fval();
  // Get output parameters
  for(int i=0;i<par.GetNrows();i++){ par(i,0)=min.UserParameters().Value(i);}
  // Get output covariance
  for(int i=0;i<par.GetNrows();i++){
    for(int j=0;j<par.GetNrows();j++){parcov(i,j)=min.UserCovariance()(i,j);}
  } 

  TMatrixT<double> vprime=ConvertToHelicesNotation(par);
  TMatrixT<double> d_alpha=vprime-val;
  
  for(int i=0;i<par.GetNrows();i++) Log(Log::Verbose) << "par " << FreeParName(i) 
						      << " "  << i 
						      << " par " <<  vprime(i,0) 
						      << " val " << val(i,0) 
						      << " delta " <<  d_alpha(i,0) 
						      << " error " << sqrt(parcov(i,i)) 
						      << " sig " << dalpha(i,0)/sqrt(fabs(parcov(i,i))) 
						      << " original cov " << sqrt(cov(i,i)) << " " <<  d_alpha(i,0)/sqrt(cov(i,i))  <<std::endl; 
  Log(Log::Verbose) << "Chi2VertexFitter::Fit() " << chi2 << " " << ndf <<  std::endl;
  isFit=true;
  return isFit;
}
