/*
 * TPTRObject.cc
 *
 *  Created on: Apr 16, 2015
 *      Author: zotz
 */

#include "SimpleFits/FitSoftware/interface/TPTRObject.h"

TPTRObject::TPTRObject(){
  isvalid_ = false;
}

TPTRObject::TPTRObject(LorentzVectorParticle A1, std::vector<LorentzVectorParticle> Taus, std::vector<LorentzVectorParticle> Neutrinos, bool isambiguous, double RotSig, bool isvalid){
  isvalid_ = isvalid;
  isambiguous_ = isambiguous;
  if(isvalid_){
	A1_ = A1;
	Taus_ = Taus;
	Neutrinos_ = Neutrinos;
	RotationSignificance_ = isambiguous_ ? 0 : RotSig;
  }
}

bool TPTRObject::isAmbiguous() const{
  if(isvalid_){
	return isambiguous_;
  }
  Logger(Logger::Error) << "TPTRObject is NOT valid!" << std::endl;
  return false;
}

LorentzVectorParticle TPTRObject::getNeutrinoMinus() const{
  if(isvalid_){
	if(isambiguous_){
	  return Neutrinos_.at(1);
	}
	else{
	  Logger(Logger::Error) << "TPTRObject has no ambiguity!" << std::endl;
	  return LorentzVectorParticle();
	}
  }
  Logger(Logger::Error) << "TPTRObject is NOT valid!" << std::endl;
  return LorentzVectorParticle();
}

LorentzVectorParticle TPTRObject::getNeutrinoPlus() const{
  if(isvalid_){
	if(isambiguous_){
	  return Neutrinos_.at(2);
	}
  }
  else{
	  Logger(Logger::Error) << "TPTRObject has no ambiguity!" << std::endl;
	  return LorentzVectorParticle();
  }
  Logger(Logger::Error) << "TPTRObject is NOT valid!" << std::endl;
  return LorentzVectorParticle();
}

std::vector<LorentzVectorParticle> TPTRObject::getNeutrinos() const{
  if(isvalid_){
	return Neutrinos_;
  }
  Logger(Logger::Error) << "TPTRObject is NOT valid!" << std::endl;
  std::vector<LorentzVectorParticle> tmp;
  return tmp;
}

LorentzVectorParticle TPTRObject::getNeutrinoZero() const{
  if(isvalid_){
	if(!isambiguous_){
	  return Neutrinos_.at(0);
	}
	else{
	  Logger(Logger::Error) << "TPTRObject has an ambiguity!" << std::endl;
	  return LorentzVectorParticle();
	}
  }
  Logger(Logger::Error) << "TPTRObject is NOT valid!" << std::endl;
  return LorentzVectorParticle();
}

LorentzVectorParticle TPTRObject::getTauMinus() const{
  if(isvalid_){
	if(isambiguous_){
	  return Taus_.at(1);
	}
	else{
	  Logger(Logger::Error) << "TPTRObject has no ambiguity!" << std::endl;
	  return LorentzVectorParticle();
	}
  }
  Logger(Logger::Error) << "TPTRObject is NOT valid!" << std::endl;
  return LorentzVectorParticle();
}

LorentzVectorParticle TPTRObject::getTauPlus() const{
  if(isvalid_){
	if(isambiguous_){
	  return Taus_.at(2);
	}
	else{
	  Logger(Logger::Error) << "TPTRObject has no ambiguity!" << std::endl;
	  return LorentzVectorParticle();
	}
  }
  Logger(Logger::Error) << "TPTRObject is NOT valid!" << std::endl;
  return LorentzVectorParticle();
}

std::vector<LorentzVectorParticle> TPTRObject::getTaus() const{
  if(isvalid_){
	return Taus_;
  }
  Logger(Logger::Error) << "TPTRObject is NOT valid!" << std::endl;
  std::vector<LorentzVectorParticle> tmp;
  return tmp;
}

LorentzVectorParticle TPTRObject::getTauZero() const{
  if(isvalid_){
	if(!isambiguous_){
	  return Taus_.at(0);
	}
	else{
	  Logger(Logger::Error) << "TPTRObject has an ambiguity!" << std::endl;
	  return LorentzVectorParticle();
	}
  }
  Logger(Logger::Error) << "TPTRObject is NOT valid!" << std::endl;
  return LorentzVectorParticle();
}

LorentzVectorParticle TPTRObject::getA1() const{
  if(isvalid_){
	return A1_;
  }
  Logger(Logger::Error) << "TPTRObject is NOT valid!" << std::endl;
  return LorentzVectorParticle();
}

double TPTRObject::getRotationSignificance() const{
  if(isvalid_){
	if(!isambiguous_) return RotationSignificance_;
	else{
	  Logger(Logger::Error) << "TPTRObject has an ambiguity and thus no rotation!" << std::endl;
	  return 0;
	}
  }
  else{
	  Logger(Logger::Error) << "TPTRObject is NOT valid!" << std::endl;
	  return 0;
  }
}

double TPTRObject::getRotSigma() const{
  return getRotationSignificance();
}

std::vector<bool> TPTRObject::CreateVectorFromAmbiguity(){
  std::vector<bool> vec;
  if(isambiguous_){
	vec.push_back(false);
	vec.push_back(true);
	vec.push_back(true);
  }
  else{
	vec.push_back(true);
	vec.push_back(false);
	vec.push_back(false);
  }
  return vec;
}
