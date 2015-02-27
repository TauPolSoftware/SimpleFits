#ifndef SimpleFits_TrackHelixVertexTools_h
#define SimpleFits_TrackHelixVertexTools_h

// system include files
#include <TMatrixT.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TVectorD.h>

#include "SimpleFits/FitSoftware/interface/Particle.h"
#include "SimpleFits/FitSoftware/interface/TrackParticle.h"
#include "SimpleFits/FitSoftware/interface/LorentzVectorParticle.h"
#include "SimpleFits/FitSoftware/interface/ErrorMatrixPropagator.h"

////////////////////////////////////////////////////////////////////////////
//
// TrackHelixVertexTools
// Converts between the 5 parameter track format (Helix/Helices) 
// 3 parameter track format for tracks constrainted by a common vertex (ConstHelix/ConstHelices) 
// and the LorentzVector of the particles
//
///////////////////////////////////////////////////////////////////////////

class TrackHelixVertexTools{
 public:
  enum FreeVertexPar{x0=0,y0,z0,NFreeVertexPar};
  enum FreeTrackPar{kappa0=3,lambda0,phi0,NFreeTrackPar};
  enum ExtraPar{BField0=0,MassOffSet=1,NExtraPar=1};
  
  TrackHelixVertexTools(std::vector<TrackParticle> &particles_,TVector3 &vguess);
  virtual ~TrackHelixVertexTools();

   virtual LorentzVectorParticle MotherLorentzVectorParticle(int pdgid);
   virtual std::vector<TrackParticle> OriginalTracks(){return particles;}
   virtual std::vector<TrackParticle> ReFitTracks();
   virtual std::vector<LorentzVectorParticle> ReFitLorentzVectorParticles();

   virtual TVector3 Vertex();
   virtual TMatrixTSym<double> VertexError();

   // Helix information for Helices Constrainted to a common vertex
   void ConstrainedHelices(TMatrixT<double> &par, TMatrixTSym<double> &cov);
   void SetConstrainedHelices(TMatrixT<double>& par,TMatrixTSym<double>& cov);
   
   // Standard 5 parameter Helix Information
   void Helices(TMatrixT<double> &par, TMatrixTSym<double> &cov);
   
 protected:
   virtual TString ConstrainedHelicesName(int Par);
   virtual void ComputeHelixParameters(TMatrixT<double> &inpar,double &kappa,double &lam,double &phi,double &dxy, double &dz, double &s, unsigned int p=0);

   TMatrixT<double> ComputeMotherLorentzVectorPar(TMatrixT<double> &inpar);
   TMatrixT<double> ConvertToLorentzVectorNotation(TMatrixT<double> &inpar);
   TMatrixT<double> ConvertToHelicesNotation(TMatrixT<double> &inpar);
   TMatrixT<double> ConvertToHelixNotation(TMatrixT<double> &inpar,int p=0);

   inline int MeasuredHelixIndex(int TrackPar,int Particle=0){return TrackPar+Particle*TrackParticle::NHelixPar;}
   inline int ConstrainedHelixIndex(int Par,int Particle=0){
     if(Par==x0 || Par==y0 || Par==z0) return Par;
     return Par+Particle*(NFreeTrackPar-NFreeVertexPar);
   }

 private:
   std::vector<TrackParticle> particles; //Original TrackParticles  
   double BField;

   // Helix information for Helices Constrainted to a common vertex
   TMatrixT<double> constHelicesPar;
   TMatrixTSym<double> constHelicesCov;
   // Standard 5 parameter Helix Information
   TMatrixT<double> helicesPar;
   TMatrixTSym<double> helicesCov;

};
#endif


