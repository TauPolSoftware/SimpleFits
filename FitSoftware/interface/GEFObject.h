/*
 * GEFObject.h
 *
 *  Created on: Apr 16, 2015
 *      Author: zotz
 */

#ifndef GEFOBJECT_H_
#define GEFOBJECT_H_

class GEFObject{
  public:
	GEFObject();
	GEFObject(std::vector<LorentzVectorParticle> InitDaughters,
			  LorentzVectorParticle InitResonance,
			  std::vector<LorentzVectorParticle> FitDaughters,
			  LorentzVectorParticle FitResonance,
			  bool fitconverged, double chi2, double csum, double Niterations, int Index);
	virtual ~GEFObject(){}
	
	double getChi2() const{ return chi2_;}
	double getCsum() const{ return csum_;}
	bool Fitconverged() const{ return fitconverged_;}
	int getIndex() const{ return Index_;}
	const LorentzVectorParticle& getInitResonance() const{ return InitResonance_;}
	const LorentzVectorParticle& getInitTauH() const{ return InitTauH_;}
	const LorentzVectorParticle& getInitTauMu() const{ return InitTauMu_;}
	bool isIsvalid() const{ return isvalid_;}
	double getNiterations() const{ return Niterations_;}
	const LorentzVectorParticle& getResonance() const{ return Resonance_;}
	const LorentzVectorParticle& getTauH() const{ return TauH_;}
	const LorentzVectorParticle& getTauMu() const{ return TauMu_;}
	
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



#endif /* GEFOBJECT_H_ */
