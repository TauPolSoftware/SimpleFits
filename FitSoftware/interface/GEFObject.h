/*
 * GEFObject.h
 *
 *  Created on: Apr 16, 2015
 *      Author: zotz
 */

#ifndef GEFObject_h
#define GEFObject_h

#include "SimpleFits/FitSoftware/interface/LorentzVectorParticle.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"

class GEFObject{
  public:
	GEFObject();
	GEFObject(const GEFObject& other);
	GEFObject(std::vector< LorentzVectorParticle > InitDaughters,
			  LorentzVectorParticle InitResonance,
			  std::vector< LorentzVectorParticle > FitDaughters,
			  LorentzVectorParticle FitResonance,
			  bool fitconverged, double chi2, double csum, double Niterations, int Index);
	virtual ~GEFObject(){}
	
	double getChi2() const{ return chi2_;}
	double getCsum() const{ return csum_;} //returns the sum of all deviations from fit parameters to both hard and soft constraints
	bool Fitconverged() const{ return fitconverged_;}
	int getIndex() const{ return Index_;} //returns the ambiguity index resolved by the fit (0 = unphysical, 1 = minus, 2 = plus)
	LorentzVectorParticle getInitResonance() const{ return InitResonance_;}
	LorentzVectorParticle getInitTauH() const{ return InitTauH_;}
	LorentzVectorParticle getInitTauMu() const{ return InitTauMu_;}
	bool isValid() const{ return isvalid_;}
	double getNiterations() const{ return Niterations_;}
	LorentzVectorParticle getResonance() const{ return Resonance_;}
	LorentzVectorParticle getTauH() const{ return TauH_;}
	LorentzVectorParticle getTauMu() const{ return TauMu_;}
	
  protected:
	
  private:
	bool isvalid_;
	bool fitconverged_;
	double chi2_;
	double csum_;
	double Niterations_;
	int Index_;
	LorentzVectorParticle InitTauH_, InitTauMu_, InitResonance_;
	LorentzVectorParticle TauH_, TauMu_, Resonance_;
};



#endif /* GEFObject_h */
