/*
 * GlobalEventFit.h
 *
 *  Created on: Apr 16, 2015
 *      Author: zotz
 */

#ifndef GlobalEventFit_h
#define GlobalEventFit_h

#include "SimpleFits/FitSoftware/interface/TPTRObject.h"
#include "SimpleFits/FitSoftware/interface/GEFObject.h"
#include "SimpleFits/FitSoftware/interface/LorentzVectorParticle.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"
#include "SimpleFits/FitSoftware/interface/TrackParticle.h"
#include "SimpleFits/FitSoftware/interface/PTObject.h"

class GlobalEventFit{
  public:
	GlobalEventFit(TrackParticle Muon, LorentzVectorParticle A1, double Phi_Res, TVector3 PV, TMatrixTSym<double> PVCov);
	GlobalEventFit(TrackParticle Muon, LorentzVectorParticle A1, PTObject METminusNeutrino, TVector3 PV, TMatrixTSym<double> PVCov);
	virtual ~GlobalEventFit();

	GEFObject Fit();
	bool AmbiguitySolverByChi2(std::vector<bool> A1Fit, std::vector<bool> EventFit, std::vector<double> Chi2s, int &IndexToReturn);
	std::vector<LorentzVectorParticle> FitDaughtersCorr(std::vector<LorentzVectorParticle> FitDaughters);
	const LorentzVectorParticle& getA1() const{ return A1_;}
	const GEFObject& getGEFObject() const{ return GEFObject_;}
	
	bool isConfigured() const{return isConfigured_;}
	bool isFit() const{return isFit_;}
	
	const TrackParticle getMuon() const{return Muon_;}
	const TVector3 getPV() const{return PV_;}
	const TMatrixTSym<double> getPVCov() const{return PVCov_;}
	const TVector3 getSV() const{return SV_;}
	const TMatrixTSym<double> getSVCov() const{return SVCov_;}
	const TPTRObject getTPTRObject() const{return TPTRObject_;}
	
	int getMaxIterations() const{return MaxIterations_;}
	void setMaxIterations(int maxIterations){
	  MaxIterations_ = maxIterations;
	  useDefaultMaxIterations_ = false;
	}
	double getMaxDelta() const{return MaxDelta_;}
	void setMaxDelta(double maxDelta){
	  MaxDelta_ = maxDelta;
	  useDefaultMaxDelta_ = false;
	}
	double getEpsilon() const{return Epsilon_;}
	void setEpsilon(double epsilon){
	  Epsilon_ = epsilon;
	  useDefaultEpsilon_ = false;
	}

	double getMassConstraint() const{return MassConstraint_;}
	void setMassConstraint(double MassConstraint){
	  MassConstraint_ = MassConstraint;
	  useDefaultMassConstraint_ = false;
	}

  protected:
	TPTRObject ThreeProngTauReco();
	bool IsAmbiguous(std::vector<bool> recostatus);
	PTObject SubtractNeutrinoFromMET(unsigned Ambiguity);

  private:
	bool isConfigured_;
	bool isFit_;
	TPTRObject TPTRObject_;
	GEFObject GEFObject_;
	TrackParticle Muon_;
	LorentzVectorParticle A1_;
	TVector3 PV_, SV_;
	TMatrixTSym<double> PVCov_, SVCov_;
	PTObject MET_;
	std::vector<PTObject> METminusNeutrino_;
	bool useMassConstraint_;
	double MassConstraint_;
	double Phi_Res_;
	int MaxIterations_;
	double MaxDelta_;
	double Epsilon_;
	bool useDefaultMaxIterations_, useDefaultMaxDelta_, useDefaultEpsilon_, useDefaultMassConstraint_;
	bool useFullRecoil_;
};


#endif /* GlobalEventFit_h */
