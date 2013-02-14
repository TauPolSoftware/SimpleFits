#include "SimpleFits/FitSoftware/interface/Chi2VertexFitter.h"
#include "SimpleFits/FitSoftware/interface/ChiSquareFunctionUpdator.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnMigrad.h"
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

bool Chi2VertexFitter::Fit(){
  if(isFit==true) return true;// do not refit
  if(!isConfigure) return false; // do not fit if configuration failed
  std::cout << "Chi2VertexFitter::Fit ready to start" << std::endl;
  ChiSquareFunctionUpdator updator(this);
  ROOT::Minuit2::MnUserParameters MnPar;
  for(int i=0;i<par.GetNrows();i++){
    TString name=FreeParName(i);
    // if not limited (vhigh <= vlow)
    std::cout << "Adding Parameter " << i << " name " << name << " " << par(i,0) << std::endl; 
    MnPar.Add(name.Data(),par(i,0),sqrt(fabs(parcov(i,i))),par(i,0)-nsigma*sqrt(fabs(parcov(i,i))),par(i,0)+nsigma*sqrt(fabs(parcov(i,i))));
  }
  // create MIGRAD minimizer
  ROOT::Minuit2::MnMigrad minimize(updator,MnPar);
  //ROOT::Minuit2::MnMinimize minimize(updator,MnPar);
  //minimize.Fix(x0);
  //minimize.Fix(y0);
  //minimize.Fix(z0);
  ROOT::Minuit2::FunctionMinimum min = minimize(0,5e-2);

  // ROOT::Minuit2::MnHesse hesse(0);// low strategy
  //ROOT::Minuit2::MnUserParameterState hessemin=hesse(updator,min.UserParameters(),min.UserCovariance());

  /*for(unsigned int i=0;i<10 && !min.IsValid();i++){
    min=FindMinimum(updator,min);
    if(min.IsValid()) break;
    }*/
  //ROOT::Minuit2::MnMinos minos(updator,min);

  // give return flag based on status
  if(!min.IsValid()){
    std::cout << "Chi2VertexFitter::Fit(): Failed min.IsValid()" << std::endl; 
    if(!min.HasValidParameters()){std::cout << "Chi2VertexFitter::Fit(): Failed min.HasValidParameters()" << std::endl; }
    if(!min.HasValidCovariance()){std::cout << "Chi2VertexFitter::Fit(): Failed min.HasValidCovariance()" << std::endl; }
    if(!min.HesseFailed()){std::cout << "Chi2VertexFitter::Fit(): Failed min.HesseFailed()" << std::endl; }
    if(!min.HasReachedCallLimit()){std::cout << "Chi2VertexFitter::Fit(): Failed min.HasReachedCallLimit()" << std::endl; }
    //exit(1);
  }
  std::cout << "Function Value " << min.Fval() << std::endl;
  if(!min.HasValidCovariance()){std::cout << "Chi2VertexFitter::Fit(): Failed min.HasValidCovariance()" << std::endl; }
  if(!min.HesseFailed()){std::cout << "Chi2VertexFitter::Fit(): Failed min.HesseFailed()" << std::endl; }
  if(!min.HasReachedCallLimit()){std::cout << "Chi2VertexFitter::Fit(): Failed min.HasReachedCallLimit()" << std::endl; }

  // Get output parameters
  for(int i=0;i<par.GetNrows();i++){ par(i,0)=min.UserParameters().Value(i); std::cout << "par " << par(i,0) << std::endl;}
  // Get output covariance
  for(int i=0;i<par.GetNrows();i++){
    std::cout << "cov" << std::endl;
    for(int j=0;j<par.GetNrows();j++){parcov(i,j)=min.UserCovariance()(i,j); std::cout << parcov(i,j) << " ";}
    std::cout << std::endl;
  }



  isFit=true;
  return isFit;
}


ROOT::Minuit2::FunctionMinimum Chi2VertexFitter::FindMinimum(ROOT::Minuit2::FCNBase &updator,ROOT::Minuit2::FunctionMinimum min){
  //MnMinimize and CombinedMinimizer
  ROOT::Minuit2::MnSimplex simplex(updator,min.UserParameters());
  ROOT::Minuit2::FunctionMinimum minsimplex=simplex();
  ROOT::Minuit2::MnMigrad migrad(updator,minsimplex.UserParameters());
  ROOT::Minuit2::FunctionMinimum minmigrad = migrad(0,5e-2);
  return minmigrad;
}
/*
http://cmslxr.fnal.gov/lxr/source/RecoVertex/BeamSpotProducer/src/BSFitter.cc#018
658         FunctionMinimum fmin = migrad();
659         ff_minimum = fmin.Fval();
660 
661         reco::BeamSpot::CovarianceMatrix matrix;
662         for (int j = 0 ; j < 6 ; ++j) {
663                 for(int k = j ; k < 6 ; ++k) {
664                         matrix(j,k) = fmin.Error().Matrix()(j,k);
665                 }
666         }
667         
668         return reco::BeamSpot( reco::BeamSpot::Point(fmin.Parameters().Vec()(0),
669                                                      fmin.Parameters().Vec()(1),
670                                                      0.),
671                                0.,
672                                fmin.Parameters().Vec()(4),
673                                fmin.Parameters().Vec()(5),
674                                0.,
675                                matrix,
676                                fbeamtype );
*/
