#ifndef LorentzVectorParticle_h
#define LorentzVectorParticle_h

#include "TauPolSoftware/SimpleFits/interface/Particle.h"
#include "TString.h"
#include "TMath.h"
#include "TLorentzVector.h"

class LorentzVectorParticle : public Particle {
 public:
  enum LorentzandVectorPar{vx=0,vy,vz,px,py,pz,m,NLorentzandVertexPar,E=-1,p=-2,pt=-3};// Lorentez+vertex parameters
  enum VertexInfo{NVertex=3};
  LorentzVectorParticle();
  LorentzVectorParticle(const LorentzVectorParticle& other);
  LorentzVectorParticle(TMatrixT<double> par_, TMatrixTSym<double> cov_, int pdgid_, double charge_, double b_);
  virtual ~LorentzVectorParticle(){};

  /*LorentzVectorParticle& operator=(const LorentzVectorParticle& other);*/

  static TString Name(int i);
  virtual int NParameters(){return NLorentzandVertexPar;}
  virtual double Parameter(int i);
  virtual double Covariance(int i,int j);
  virtual double Mass(){return Parameter(m);}
  TLorentzVector LV(){return TLorentzVector(Parameter(px),Parameter(py),Parameter(pz),Parameter(E));}
  TMatrixTSym<double> LVCov(){
	TMatrixTSym<double> lvcov(4);
    for(int i=px;i<=pz;i++){
      for(int j=px;j<pz;j++){lvcov(i,j)=Covariance(i,j);} // 3x3 matrix
      lvcov(4,i) = Covariance(E,i);
      lvcov(i,4) = lvcov(4,i);
    }
    lvcov(4,4) = Covariance(E,E);
    return lvcov;
  }
  TVector3 Vertex(){return TVector3(Parameter(vx),Parameter(vy),Parameter(vz));}
  TMatrixTSym<double> VertexCov(){
    TMatrixTSym<double> vcov(NVertex);
    for(int i=0;i<NVertex;i++){
      for(int j=0;j<NVertex;j++){vcov(i,j)=Covariance(i,j);}
    }
    return vcov; 
  }
};
#endif


