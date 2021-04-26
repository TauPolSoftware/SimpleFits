#ifndef ChiSquareFunctionUpdator_h
#define ChiSquareFunctionUpdator_h

#include "Minuit2/FCNBase.h"
#include "TMatrixT.h"
#include "TauPolSoftware/SimpleFits/interface/TrackHelixVertexFitter.h"
#include "TauPolSoftware/SimpleFits/interface/LagrangeMultipliersFitter.h"
#include "TauPolSoftware/SimpleFits/interface/DiTauConstrainedFitter.h"
#include "TauPolSoftware/SimpleFits/interface/Logger.h"

class ChiSquareFunctionUpdator  : public ROOT::Minuit2::FCNBase {
 public:
  ChiSquareFunctionUpdator(TrackHelixVertexFitter *VF_){VF=VF_;}
  virtual ~ChiSquareFunctionUpdator(){};

  virtual double operator() (const std::vector<double> & x)const{
    TMatrixT<double> X(x.size(),1);
    for(unsigned int i=0; i<x.size();i++){X(i,0)=x.at(i);}
    return VF->UpdateChisquare(X);
  }
  virtual double Up()const{return 1.0;}// Error definiton for Chi^2

 private:
  TrackHelixVertexFitter *VF;

};

class LagrangeMultipliersFitterChiSquareFunctionUpdator  : public ROOT::Minuit2::FCNBase {
 public:
  LagrangeMultipliersFitterChiSquareFunctionUpdator(LagrangeMultipliersFitter *LMF_){LMF=LMF_;}
  virtual ~LagrangeMultipliersFitterChiSquareFunctionUpdator(){};

  virtual double operator() (const std::vector<double> & x)const{
    int offset = LMF->NPara();
    TVectorD a(LMF->NPara()),b(LMF->NParb());
    for(int i=0; i<LMF->NPara(); i++) a(i)=x.at(i);
    for(int i=0; i<LMF->NParb(); i++) b(i)=x.at(i+offset);
    return LMF->UpdateChisquare(a,b);
  }
  virtual double Up()const{return 1.0;}// Error definiton for Chi^2

 private:
  LagrangeMultipliersFitter *LMF;

};

class DiTauConstrainedFitterChiSquareFunctionUpdator  : public ROOT::Minuit2::FCNBase {
 public:
  DiTauConstrainedFitterChiSquareFunctionUpdator(DiTauConstrainedFitter *DTCF_){DTCF=DTCF_;}
  virtual ~DiTauConstrainedFitterChiSquareFunctionUpdator(){};

  virtual double operator() (const std::vector<double> & x)const{
    int offset = DTCF->NPara();
    TVectorD a(DTCF->NPara()),b(DTCF->NParb());
    for(int i=0; i<DTCF->NPara(); i++) a(i)=x.at(i);
    for(int i=0; i<DTCF->NParb(); i++) b(i)=x.at(i+offset);
    return DTCF->UpdateChisquare(a,b);
  }
  virtual double Up()const{return 1.0;}// Error definiton for Chi^2

 private:
  DiTauConstrainedFitter *DTCF;

};
#endif

