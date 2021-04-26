/*
 * GlobalEventFit.h
 *
 *  Created on: Apr 16, 2015
 *      Author: zotz
 */

#ifndef GlobalEventFit_h
#define GlobalEventFit_h

#include "TauPolSoftware/SimpleFits/interface/TPTRObject.h"
#include "TauPolSoftware/SimpleFits/interface/GEFObject.h"
#include "TauPolSoftware/SimpleFits/interface/LorentzVectorParticle.h"
#include "TauPolSoftware/SimpleFits/interface/Logger.h"
#include "TauPolSoftware/SimpleFits/interface/TrackParticle.h"
#include "TauPolSoftware/SimpleFits/interface/PTObject.h"

class GlobalEventFit{
  public:
	GlobalEventFit(TrackParticle Muon, LorentzVectorParticle A1, double Phi_Res, TVector3 PV, TMatrixTSym<double> PVCov);
	GlobalEventFit(TrackParticle Muon, LorentzVectorParticle A1, PTObject METminusNeutrino, TVector3 PV, TMatrixTSym<double> PVCov);
	virtual ~GlobalEventFit();

	GEFObject Fit();
	bool AmbiguitySolverByChi2(std::vector<bool> A1Fit, std::vector<bool> EventFit, std::vector<TVectorD> Chi2Vecs, int &IndexToReturn);
	bool AmbiguitySolverByChi2Minuit(std::vector<bool> A1Fit, std::vector<bool> EventFit, std::vector<double> Chi2s, int &IndexToReturn);
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
	std::vector<bool> getFitStatuses() const{return fitstatuses_;}

	double getMassConstraint() const{return MassConstraint_;}
	void setMassConstraint(double MassConstraint){
	  MassConstraint_ = MassConstraint;
	  useDefaultMassConstraint_ = false;
	}
	void SetCorrectPt(bool correct){correctPt_ = correct;}  //default is set to true. Can be set to false inside your analysis to prevent the corrections of the reconstructed pt of both taus to preserve the hard constraints imposed on the fit resonance
	void setUseCollinearityTauMu(bool useCollinearityTauMu){useCollinearityTauMu_ = useCollinearityTauMu;}

  protected:
	void Configure(TrackParticle Muon, LorentzVectorParticle A1, TVector3 PV, TMatrixTSym<double> PVCov);
	TPTRObject ThreeProngTauReconstruction();
	bool IsAmbiguous(std::vector<bool> recostatus);
	PTObject SubtractNeutrinoFromMET(unsigned Ambiguity);
	PTObject AddA1(PTObject MET);
	PTObject AddMuon(PTObject MET);

  private:
	bool isConfigured_;
	bool isFit_;
	bool isValid_;
	TPTRObject TPTRObject_;
	GEFObject GEFObject_;
	TrackParticle Muon_;
	LorentzVectorParticle A1_;
	TVector3 PV_, SV_;
	TMatrixTSym<double> PVCov_, SVCov_;
	PTObject MET_;
	std::vector<PTObject> METminusNeutrino_;
	std::vector<bool> fitstatuses_;
	TMatrixD FitPar_;
	TMatrixDSym FitCov_;
	bool useMassConstraint_;
	double MassConstraint_;
	double Phi_Res_;
	int MaxIterations_;
	double MaxDelta_;
	double Epsilon_;
	bool useDefaultMassConstraint_, useFullRecoil_, correctPt_, useCollinearityTauMu_;
};


#endif /* GlobalEventFit_h */
