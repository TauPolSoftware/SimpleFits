/*
 * GEFObject.cc
 *
 *  Created on: Apr 16, 2015
 *      Author: zotz
 */

#include "SimpleFits/FitSoftware/interface/GEFObject.h"

GEFObject::GEFObject(){
  isvalid_ = false;
}
GEFObject::GEFObject(std::vector<LorentzVectorParticle> InitDaughters,
					 LorentzVectorParticle InitResonance,
					 std::vector<LorentzVectorParticle> FitDaughters,
					 LorentzVectorParticle FitResonance,
					 bool fitconverged, double chi2, double csum, double Niterations, int Index){
  isvalid_ = true;
  chi2_ = chi2;
  csum_ = csum;
  Niterations_ = Niterations;
  fitconverged_ = fitconverged;
  Index_ = Index;
  InitTauH_ = InitDaughters.at(0);
  InitTauMu_ = InitDaughters.at(1);
  InitResonance_ = InitResonance;
  TauH_ = FitDaughters.at(0);
  TauMu_ = FitDaughters.at(1);
  Resonance_ = FitResonance;
};
