/*
 * TPTRObject.cc
 *
 *  Created on: Apr 16, 2015
 *      Author: zotz
 */

#include "TPTRObject.h"

TPTRObject::TPTRObject(){
  isvalid_ = false;
}

TPTRObject::TPTRObject(LorentzVectorParticle A1, std::vector<LorentzVectorParticle> Taus, std::vector<LorentzVectorParticle> Neutrinos, bool isvalid, bool isambiguous){
  isvalid_ = isvalid;
  isambiguous_ = isambiguous;
  if(isvalid_){
	A1_ = A1;
	Taus_ = Taus;
	Neutrinos_ = Neutrinos;
  }
}

const LorentzVectorParticle& TPTRObject::getNeutrinoMinus() const{
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

const LorentzVectorParticle& TPTRObject::getNeutrinoPlus() const{
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

const std::vector<LorentzVectorParticle>& TPTRObject::getNeutrinos() const{
  if(isvalid_){
	return Neutrinos_;
  }
  Logger(Logger::Error) << "TPTRObject is NOT valid!" << std::endl;
  return LorentzVectorParticle();
}

const LorentzVectorParticle& TPTRObject::getNeutrinoZero() const{
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

const LorentzVectorParticle& TPTRObject::getTauMinus() const{
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

const LorentzVectorParticle& TPTRObject::getTauPlus() const{
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

const std::vector<LorentzVectorParticle>& TPTRObject::getTaus() const{
  if(isvalid_){
	return Taus_;
  }
  Logger(Logger::Error) << "TPTRObject is NOT valid!" << std::endl;
  return LorentzVectorParticle();
}

const LorentzVectorParticle& TPTRObject::getTauZero() const{
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

const LorentzVectorParticle& TPTRObject::getA1() const{
  if(isvalid_){
	return A1_;
  }
  Logger(Logger::Error) << "TPTRObject is NOT valid!" << std::endl;
  return LorentzVectorParticle();
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
