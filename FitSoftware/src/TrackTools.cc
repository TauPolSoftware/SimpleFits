#include "SimpleFits/FitSoftware/interface/TrackTools.h"
#include "math.h"
#include <iostream>

TVector3 TrackTools::PropogateToXPosition(TrackParticle &p,double &x){
  double kappa=p.Parameter(TrackParticle::kappa);
  double lam=p.Parameter(TrackParticle::lambda);
  double phi=p.Parameter(TrackParticle::phi);
  double dxy=p.Parameter(TrackParticle::dxy);
  double dz=p.Parameter(TrackParticle::dz);
  double r=kappa/2.0;
  double s=(asin((x+(r+dxy)*sin(phi))/r)-phi)/(2*kappa);
  double y=-r*cos(2.0*s*kappa+phi)+(r+dxy)*cos(phi);
  double z=dz+s*tan(lam);
  return TVector3(x,y,z);
}

TVector3 TrackTools::PropogateToYPosition(TrackParticle &p,double &y){
  double kappa=p.Parameter(TrackParticle::kappa);
  double lam=p.Parameter(TrackParticle::lambda);
  double phi=p.Parameter(TrackParticle::phi);
  double dxy=p.Parameter(TrackParticle::dxy);
  double dz=p.Parameter(TrackParticle::dz);
  double r=kappa/2.0;
  double s=(acos(((r+dxy)*cos(phi)-y)/r)-phi)/(2*kappa);
  double x=r*sin(2.0*s*kappa+phi)-(r+dxy)*sin(phi);
  double z=dz+s*tan(lam);
  return TVector3(x,y,z);
}

TVector3 TrackTools::PropogateToZPosition(TrackParticle &p,double &z){
  double kappa=p.Parameter(TrackParticle::kappa);
  double lam=p.Parameter(TrackParticle::lambda);
  double phi=p.Parameter(TrackParticle::phi);
  double dxy=p.Parameter(TrackParticle::dxy);
  double dz=p.Parameter(TrackParticle::dz);
  double s=(z-dz)/tan(lam);
  double r=kappa/2.0;
  double x=r*sin(2.0*s*kappa+phi)-(r+dxy)*sin(phi);
  double y=-r*cos(2.0*s*kappa+phi)+(r+dxy)*cos(phi);
  return TVector3(x,y,z);
}
