/*
 * TPTRObject.h
 *
 *  Created on: Apr 15, 2015
 *      Author: zotz
 */

#ifndef TPTRObject_h
#define TPTRObject_h

#include <vector>
#include "SimpleFits/FitSoftware/interface/LorentzVectorParticle.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"

class TPTRObject{
  public:
	TPTRObject();
	TPTRObject(LorentzVectorParticle A1, std::vector<LorentzVectorParticle> Taus, std::vector<LorentzVectorParticle> Neutrinos, bool isambiguous, double RotSig, bool isvalid);
	virtual ~TPTRObject(){};

	std::vector<bool> CreateVectorFromAmbiguity();
	bool isAmbiguous() const;

	LorentzVectorParticle getNeutrinoZero() const;
	LorentzVectorParticle getNeutrinoMinus() const;
	LorentzVectorParticle getNeutrinoPlus() const;
	std::vector<LorentzVectorParticle> getNeutrinos() const;

	LorentzVectorParticle getTauZero() const;
	LorentzVectorParticle getTauMinus() const;
	LorentzVectorParticle getTauPlus() const;
	std::vector<LorentzVectorParticle> getTaus() const;

	LorentzVectorParticle getA1() const;

	double getRotSigma() const;
	double getRotationSignificance() const;

  private:
	bool isvalid_;
	bool isambiguous_; //true = 2 physical solutions, false = 1 unphysical solution projected to maximal GJ angle
	std::vector<LorentzVectorParticle> Taus_;
	std::vector<LorentzVectorParticle> Neutrinos_;
	LorentzVectorParticle A1_;
	double RotationSignificance_;
};


#endif /* TPTRObject_h */
