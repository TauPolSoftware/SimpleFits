#include "TauPolSoftware/SimpleFits/interface/LorentzVectorParticle.h"


LorentzVectorParticle::LorentzVectorParticle():
  Particle(TMatrixT<double>(NLorentzandVertexPar,1),TMatrixTSym<double>(NLorentzandVertexPar),0,0,0)
{
}

LorentzVectorParticle::LorentzVectorParticle(const LorentzVectorParticle& other):
  Particle(other.getParMatrix(), other.getCovMatrix(), other.PDGID(), other.Charge(), other.BField())
{
}


LorentzVectorParticle::LorentzVectorParticle(TMatrixT<double> par_,TMatrixTSym<double> cov_,int pdgid_,double charge_,double b_):
  Particle(par_,cov_,pdgid_,charge_,b_)
{
}

/*LorentzVectorParticle& LorentzVectorParticle::operator=(const LorentzVectorParticle& other){
  return LorentzVectorParticle(other);
}
*/
TString LorentzVectorParticle::Name(int i){
  if(i==px)  return "px";
  if(i==py)  return "py";
  if(i==pz)  return "pz";
  if(i==m)   return "m";
  if(i==vx)  return "vx";
  if(i==vy)  return "vy";
  if(i==vz)  return "vz";
  return "invalid";
}

double LorentzVectorParticle::Parameter(int i){
  if(i==E)  return sqrt(pow(Particle::Parameter(m),2.0)+pow(Particle::Parameter(px),2.0)+pow(Particle::Parameter(py),2.0)+pow(Particle::Parameter(pz),2.0));
  if(i==p)  return sqrt(pow(Particle::Parameter(px),2.0)+pow(Particle::Parameter(py),2.0)+pow(Particle::Parameter(pz),2.0));
  if(i==pt) return sqrt(pow(Particle::Parameter(px),2.0)+pow(Particle::Parameter(py),2.0));
  return Particle::Parameter(i);
}

double LorentzVectorParticle::Covariance(int i,int j){
	if(i==E && j==E)
		return pow(Parameter(E), -2) * (
				pow(Parameter(px),2)*Particle::Covariance(px,px) + pow(Parameter(py),2)*Particle::Covariance(py,py) + pow(Parameter(pz),2)*Particle::Covariance(pz,pz) + Parameter(m)*Parameter(m)* Particle::Covariance(m,m) +
				2*Parameter(px)*Parameter(py)*Particle::Covariance(px,py) + 2*Parameter(px)*Parameter(pz)*Particle::Covariance(px,pz) + 2*Parameter(m)*Parameter(px)*Particle::Covariance(px,m) +
				2*Parameter(py)*Parameter(pz)*Particle::Covariance(py,pz) + 2*Parameter(m)*Parameter(py)*Particle::Covariance(py,m) +
				2*Parameter(m)*Parameter(pz)*Particle::Covariance(pz,m));
	if( i==E && (j==px || j==py || j==pz) )
		return -(Particle::Covariance(j,m)*Parameter(m) + Particle::Covariance(j,px)*Parameter(px) +
				 Particle::Covariance(j,py)*Parameter(py) + Particle::Covariance(j,pz)*Parameter(pz)) / Parameter(E);
	if (j==E && (i==px || i==py || i==pz))
		return Covariance(j,i);

	return Particle::Covariance(i,j);
}
