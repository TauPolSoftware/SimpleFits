/*
 * GEFObject.cc
 *
 *  Created on: Apr 16, 2015
 *      Author: zotz
 */

#include "TauPolSoftware/SimpleFits/interface/GEFObject.h"

GEFObject::GEFObject(){
  isvalid_ = false;
  fitconverged_ = false;
  Index_ = -1;
}

GEFObject::GEFObject(const GEFObject& other){
  isvalid_ = other.isValid();
  chi2_ = other.getChi2Vectors();
  csum_ = other.getallCsums();
  Niterations_ = other.getallNiterations();
  fitconverged_ = other.Fitconverged();
  Index_ = other.getIndex();
  InitTauH_ = other.getInitTauH();
  InitTauMu_ = other.getInitTauMu();
  InitResonance_ = other.getInitResonance();
  TauH_ = other.getTauH();
  TauMu_ = other.getTauMu();
  Resonance_ = other.getResonance();
  InitTauHs_ = other.getInitTauHs();
  InitTauMus_ = other.getInitTauMus();
  InitResonances_ = other.getInitResonances();
  TauHs_ = other.getTauHs();
  TauMus_ = other.getTauMus();
  Resonances_ = other.getResonances();
}

GEFObject::GEFObject(std::vector< std::vector<LorentzVectorParticle> > InitDaughters,
					 std::vector<LorentzVectorParticle> InitResonance,
					 std::vector< std::vector<LorentzVectorParticle> > FitDaughters,
					 std::vector<LorentzVectorParticle> FitResonance,
					 bool fitconverged, std::vector<TVectorD> chi2, std::vector<double> csum, std::vector<double> Niterations, int Index){
  isvalid_ = true;
  chi2_ = chi2;
  csum_ = csum;
  Niterations_ = Niterations;
  fitconverged_ = fitconverged;
  Index_ = Index;
  InitTauH_ = InitDaughters.at(Index).at(0);
  InitTauMu_ = InitDaughters.at(Index).at(1);
  InitResonance_ = InitResonance.at(Index);
  TauH_ = FitDaughters.at(Index).at(0);
  TauMu_ = FitDaughters.at(Index).at(1);
  Resonance_ = FitResonance.at(Index);
  for(unsigned i=0; i<3; i++){
	InitTauHs_.push_back(InitDaughters.at(i).at(0));
	InitTauMus_.push_back(InitDaughters.at(i).at(1));
	InitResonances_.push_back(InitResonance.at(i));
	TauHs_.push_back(FitDaughters.at(i).at(0));
	TauMus_.push_back(FitDaughters.at(i).at(1));
	Resonances_.push_back(FitResonance.at(i));
  }
};
