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
					 bool fitvalid, bool fitconverged, std::vector<TVectorD> chi2, std::vector<double> csum, std::vector<double> Niterations, int Index)
{
  for(unsigned i=0; i<InitDaughters.size(); i++){
    InitTauHs_.push_back(InitDaughters.at(i).at(0));
    InitTauMus_.push_back(InitDaughters.at(i).at(1));
    InitResonances_.push_back(InitResonance.at(i));
  }
  isvalid_ = fitvalid;
  chi2_ = chi2;
  csum_ = csum;
  Niterations_ = Niterations;
  fitconverged_ = fitconverged;
  // std::cout << "isvalid_: " << isvalid_ << "\n";
  // std::cout << "fitconverged_: " << fitconverged_ << "\n";
  // std::cout << "Index: " << Index << "\n";
  if(fitconverged){
    Index_ = Index;
    InitTauH_ = InitDaughters.at(Index_).at(0);
    InitTauMu_ = InitDaughters.at(Index_).at(1);
    InitResonance_ = InitResonance.at(Index_);
    TauH_ = FitDaughters.at(Index_).at(0);
    TauMu_ = FitDaughters.at(Index_).at(1);
    Taus_.push_back(TauH_);
    Taus_.push_back(TauMu_);
    Resonance_ = FitResonance.at(Index_);
    for(unsigned i=0; i<InitDaughters.size(); i++){
      TauHs_.push_back(FitDaughters.at(i).at(0));
      TauMus_.push_back(FitDaughters.at(i).at(1));
      Resonances_.push_back(FitResonance.at(i));
    }
  }
  else
  {
    Index_ = 0;
    InitTauH_ = LorentzVectorParticle();
    InitTauMu_ = LorentzVectorParticle();
    InitResonance_ = LorentzVectorParticle();
    TauH_ = LorentzVectorParticle();
    TauMu_ = LorentzVectorParticle();
    Taus_.push_back(TauH_);
    Taus_.push_back(TauMu_);
    Resonance_ = LorentzVectorParticle();
    for(unsigned i=0; i<InitDaughters.size(); i++){
      TauHs_.emplace_back(LorentzVectorParticle());
      TauMus_.emplace_back(LorentzVectorParticle());
      Resonances_.emplace_back(LorentzVectorParticle());
    }
  }
}