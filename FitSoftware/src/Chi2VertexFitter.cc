#include "SimpleFits/FitSoftware/interface/Chi2VertexFitter.h"
#include "SimpleFits/FitSoftware/interface/ChiSquareFunctionUpdator.h"
#include "TFitterMinuit.h"

bool Chi2VertexFitter::Fit(){
  if(isFit==true) return true;// do not refit
  TFitterMinuit minuit;
  ChiSquareFunctionUpdator updator(this);
  minuit.SetMinuitFCN(&updator);
  int j(0),p(0);
  for(int i=0;i<par.GetNrows();i++){
    TString name=TrackParticle::Name(j); name+=p;
    // if not limited (vhigh <= vlow)
    minuit.SetParameter(i,name,par(i,0),0.5*parcov(i,i),0,0);
    j++;
    if(j==TrackParticle::NHelixPar){j=0;p++;}
  }
  //minuit.SetPrintLevel(3);
  // create Minimizer (default is Migrad)
  minuit.CreateMinimizer();
  int iret = minuit.Minimize();
  if (iret != 0) { 
    return isFit; 
  }
  for(int i=0;i<par.GetNrows();i++){
    par(i,0)=minuit.GetParameter(i);
    for(int j=0;j<par.GetNrows();j++){parcov(i,j)=minuit.GetCovarianceMatrixElement(i,j);}
  }
  isFit=true;
  return isFit;
}


