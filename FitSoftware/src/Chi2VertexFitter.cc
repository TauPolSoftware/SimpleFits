#include "SimpleFits/FitSoftware/interface/Chi2VertexFitter.h"
#include "SimpleFits/FitSoftware/interface/ChiSquareFunctionUpdator.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnContours.h"
#include "Minuit2/MnPlot.h"
#include "Minuit2/MinosError.h"
#include "Minuit2/ContoursError.h"
#include <iostream>

bool Chi2VertexFitter::Fit(){
  if(isFit==true) return true;// do not refit
  ChiSquareFunctionUpdator updator(this);
  ROOT::Minuit2::MnUserParameters upar;
  for(int i=0;i<par.GetNrows();i++){
    TString name=FreeParName(i);
    // if not limited (vhigh <= vlow)
    std::cout << "Adding Parameter " << i << " name " << name << std::endl; 
    upar.Add(name.Data(),par(i,0),parcov(i,i),par(i,0)-10*parcov(i,i),par(i,0)+10*parcov(i,i));
  }
  // Setup function 

  // create MIGRAD minimizer
  ROOT::Minuit2::MnMigrad migrad(updator,upar);
  
  // Minimize
  ROOT::Minuit2::FunctionMinimum min = migrad();

  // give return flag based on status
  if(min.IsValid())             return 1;
  if(min.HasValidParameters())  return 2;
  if(min.HasValidCovariance())  return 3;
  if(min.HesseFailed())         return 4;
  if(min.HasReachedCallLimit()) return 5;
  for(int i=0;i<par.GetNrows();i++){
    par(i,0)=min.UserParameters().Value(i);
    for(int j=0;j<par.GetNrows();j++){std::cout << i << " " << j << std::endl; parcov(i,j)=min.UserCovariance()(i,j);}
  }
  isFit=true;
  return isFit;
}
