/*
 * GEFObject.cc
 *
 *  Created on: Apr 16, 2015
 *      Author: zotz
 */

#include "SimpleFits/FitSoftware/interface/GEFObject.h"

GEFObject::GEFObject(){
  isvalid_ = false;
  chi2_ = -1;
  csum_ = -1;
  Niterations_ = -1;
  fitconverged_ = false;
  Index_ = -1;
}

GEFObject::GEFObject(const GEFObject& other){
  isvalid_ = other.isValid();
  chi2_ = other.getChi2();
  csum_ = other.getCsum();
  Niterations_ = other.getNiterations();
  fitconverged_ = other.Fitconverged();
  Index_ = other.getIndex();
  InitTauH_ = other.getInitTauH();
  InitTauMu_ = other.getInitTauMu();
  InitResonance_ = other.getInitResonance();
  TauH_ = other.getTauH();
  TauMu_ = other.getTauMu();
  Resonance_ = other.getResonance();
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
