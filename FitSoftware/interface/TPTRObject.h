/*
 * TPTRObject.h
 *
 *  Created on: Apr 15, 2015
 *      Author: zotz
 */

#ifndef TPTROBJECT_H_
#define TPTROBJECT_H_

#include <vector>
#include "SimpleFits/FitSoftware/interface/LorentzVectorParticle.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"

class TPTRObject{
  public:
	TPTRObject();
	TPTRObject(LorentzVectorParticle A1, std::vector<LorentzVectorParticle> Taus, std::vector<LorentzVectorParticle> Neutrinos, bool isambiguous, bool isvalid);
	virtual ~TPTRObject(){};

	std::vector<bool> CreateVectorFromAmbiguity();

	const LorentzVectorParticle& getNeutrinoZero() const;
	const LorentzVectorParticle& getNeutrinoMinus() const;
	const LorentzVectorParticle& getNeutrinoPlus() const;
	const std::vector<LorentzVectorParticle>& getNeutrinos() const;

	const LorentzVectorParticle& getTauZero() const;
	const LorentzVectorParticle& getTauMinus() const;
	const LorentzVectorParticle& getTauPlus() const;
	const std::vector<LorentzVectorParticle>& getTaus() const;

	const LorentzVectorParticle& getA1() const;

  private:
	bool isvalid_;
	bool isambiguous_; //true = 2 physical solutions, false = 1 unphysical solution projected to maximal GJ angle
	std::vector<LorentzVectorParticle> Taus_;
	std::vector<LorentzVectorParticle> Neutrinos_;
	LorentzVectorParticle A1_;
};


#endif /* TPTROBJECT_H_ */
