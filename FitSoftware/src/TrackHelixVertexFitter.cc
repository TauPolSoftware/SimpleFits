#include "SimpleFits/FitSoftware/interface/TrackHelixVertexFitter.h"
#include "SimpleFits/FitSoftware/interface/SimpleFits.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"
#include "TDecompBK.h"
#include <iostream>

TrackHelixVertexFitter::TrackHelixVertexFitter(std::vector<TrackParticle> &particles_,TVector3 &vguess):
  TrackHelixVertexTools(particles_,vguess),
  isFit(false),
  isConfigure(false),
  nParticles(particles_.size())
{
  ConstrainedHelices(par,parcov);
  Helices(val,cov);

  TDecompBK Inverter(cov);
  double det = cov.Determinant();
  if(!Inverter.Decompose()){
    Log(Log::Fatal) << "TrackHelixVertexFitter::TrackHelixVertexFitter Fit failed: unable to invert SYM gain matrix " << det << " \n" << std::endl;
    exit(1);
  }
  cov_inv.ResizeTo(val.GetNrows(),val.GetNrows());
  cov_inv=Inverter.Invert();
  ndf=val.GetNrows()-par.GetNrows();
  isConfigure=true;
}

TrackHelixVertexFitter::~TrackHelixVertexFitter(){}

double TrackHelixVertexFitter::UpdateChisquare(TMatrixT<double> inpar){
  TMatrixT<double> vprime=ConvertToHelicesNotation(inpar);
  TMatrixT<double> dalpha=vprime-val;
  TMatrixT<double> dalphaT=dalpha;  dalphaT.T();
  TMatrixT<double> chisquare=dalphaT*(cov_inv*dalpha);
  return chisquare(0,0);
}

