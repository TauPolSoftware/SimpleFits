/*
 * GlobalEventFit.h
 *
 *  Created on: Apr 16, 2015
 *      Author: zotz
 */

#ifndef GLOBALEVENTFIT_H_
#define GLOBALEVENTFIT_H_

#include "TPTRObject.h"
#include "GEFObject.h"

class GlobalEventFit{
  public:
	GlobalEventFit(TrackParticle Muon, LorentzVectorParticle A1, TVector3 PV, TMatrixTSym<double> PVCov);
	virtual ~GlobalEventFit();

	void Configure();
	void ReConfigure(unsigned int MuonIndex,unsigned int TauIndex);
	GEFObject Fit();
	bool AmbiguitySolverByChi2(std::vector<bool> A1Fit, std::vector<bool> EventFit, std::vector<double> Chi2s, int &IndexToReturn);

	const LorentzVectorParticle& getA1() const{ return A1_;}
	const GEFObject& getGEFObject() const{ return GEFObject_;}
	
	bool isConfigured() const{return isConfigured_;}
	bool isFit() const{return isFit_;}
	
	const TrackParticle& getMuon() const{return Muon_;}
	const TVector3& getPv() const{return PV_;}
	const TMatrixTSym<double>& getPvCov() const{return PVCov_;}
	const TVector3& getSv() const{return SV_;}
	const TMatrixTSym<double>& getSvCov() const{return SVCov_;}
	const TPTRObject& getTPTRObject() const{return TPTRObject_;}

  protected:
	TPTRObject ThreeProngTauReco();
	bool IsAmbiguous(std::vector<bool> recostatus);

  private:
	bool isConfigured_;
	bool isFit_;
	TPTRObject TPTRObject_;
	GEFObject GEFObject_;
	TrackParticle Muon_;
	LorentzVectorParticle A1_;
	TVector3 PV_, SV_;
	TMatrixTSym<double> PVCov_, SVCov_;
};


#endif /* GLOBALEVENTFIT_H_ */
