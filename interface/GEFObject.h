/*
 * GEFObject.h
 *
 *  Created on: Apr 16, 2015
 *      Author: zotz#
 *
 *      Container object for all fit results
 */

#ifndef GEFObject_h
#define GEFObject_h

#include "TVectorD.h"
#include "TauPolSoftware/SimpleFits/interface/LorentzVectorParticle.h"
#include "TauPolSoftware/SimpleFits/interface/Logger.h"

class GEFObject{
  public:
	GEFObject();
	GEFObject(const GEFObject& other);
	GEFObject(std::vector< std::vector<LorentzVectorParticle> > InitDaughters,
			  std::vector<LorentzVectorParticle> InitResonance,
			  std::vector< std::vector<LorentzVectorParticle> > FitDaughters,
			  std::vector<LorentzVectorParticle> FitResonance,
			  bool fitconverged, std::vector<TVectorD> chi2, std::vector<double> csum, std::vector<double> Niterations, int Index);
	virtual ~GEFObject(){}
	
	double getChi2() const{return chi2_.at(Index_).Sum();} // returns full "chi2" of the solution picked by the fit
	TVectorD getChi2Vector() const{ return chi2_.at(Index_);} // returns a vector with each chi2 contribution from original chi2, soft and hard constraints (in this order) for the picked solution of the fit
	std::vector<TVectorD> getChi2Vectors() const{ return chi2_;} // returns a vector with each chi2 contribution from original chi2, soft and hard constraints (in this order) for each solution of the fit
	double getCsum() const{ return csum_.at(Index_);} //returns the sum of all deviations from fit parameters to both hard and soft constraints
	bool Fitconverged() const{ return fitconverged_;}
	int getIndex() const{ return Index_;} //returns the ambiguity index resolved by the fit (0 = unphysical, 1 = minus, 2 = plus)

	// return the input estimates of the fit particles for the resolved ambiguity
	LorentzVectorParticle getInitResonance() const{ return InitResonance_;}
	LorentzVectorParticle getInitTauH() const{ return InitTauH_;}
	LorentzVectorParticle getInitTauMu() const{ return InitTauMu_;}

	// return the input estimates of the fit particles for all possible ambiguities
	std::vector<LorentzVectorParticle> getInitResonances() const{ return InitResonances_;}
	std::vector<LorentzVectorParticle> getInitTauHs() const{ return InitTauHs_;}
	std::vector<LorentzVectorParticle> getInitTauMus() const{ return InitTauMus_;}

	// return the fit particles for the resolved ambiguity
	LorentzVectorParticle getResonance() const{ return Resonance_;}
	LorentzVectorParticle getTauH() const{ return TauH_;}
	LorentzVectorParticle getTauMu() const{ return TauMu_;}
	
	// return the fit particles for all possible ambiguities
	std::vector<LorentzVectorParticle> getResonances() const{ return InitResonances_;}
	std::vector<LorentzVectorParticle> getTauHs() const{ return InitTauHs_;}
	std::vector<LorentzVectorParticle> getTauMus() const{ return InitTauMus_;}

	bool isValid() const{ return isvalid_;}
	double getNiterations() const{ return Niterations_.at(Index_);}

	std::vector<double> getallCsums() const{ return csum_;}
	std::vector<double> getallNiterations() const{ return Niterations_;}

  protected:
	
  private:
	bool isvalid_;
	bool fitconverged_;
	std::vector<TVectorD> chi2_;
	std::vector<double> csum_, Niterations_;
	int Index_;
	LorentzVectorParticle InitTauH_, InitTauMu_, InitResonance_;
	LorentzVectorParticle TauH_, TauMu_, Resonance_;
	std::vector<LorentzVectorParticle> InitTauHs_, InitTauMus_, InitResonances_;
	std::vector<LorentzVectorParticle> TauHs_, TauMus_, Resonances_;
};



#endif /* GEFObject_h */
