#include "SimpleFits/FitSoftware/interface/FitterBase.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"
#include <iostream>

FitterBase::~FitterBase(){

}


void FitterBase::ReSet(){
  isFit=false;
  TDecompBK Inverter(cov);
  double det = cov.Determinant();
  isConfigure=true;
  if(!Inverter.Decompose()){
    Log(Log::Fatal) << "TrackHelixVertexFitter::TrackHelixVertexFitter Fit failed: unable to invert SYM gain matrix " << det << " \n" << std::endl;
    isConfigure=false;
    return;
  }
  cov_inv.ResizeTo(val.GetNrows(),val.GetNrows());
  cov_inv=Inverter.Invert();
  dalpha.ReSizeTo(val.GetNrows());
  ndf=val.GetNrows()-par.GetNrows();
}
