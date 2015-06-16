/*
 * GlobalEventFit.cc
 *
 *  Created on: Apr 16, 2015
 *      Author: zotz
 *
 *
 *      This class offers an interface for the global event fit including the reconstruction of the hadronic tau and the fit of the di-tau system.
 *      After the class has been instantiated, the fit can be performed using the Fit() function, which returns a GEFObject that contains the fit results.
 */

#include "SimpleFits/FitSoftware/interface/GlobalEventFit.h"
#include "SimpleFits/FitSoftware/interface/TauA1NuConstrainedFitter.h"
#include "SimpleFits/FitSoftware/interface/DiTauConstrainedFitter.h"

GlobalEventFit::GlobalEventFit(TrackParticle Muon, LorentzVectorParticle A1, double Phi_Res, TVector3 PV, TMatrixTSym<double> PVCov){
	isConfigured_ = false;
	isFit_ = false;
	Muon_ = Muon;
	A1_= A1;
	PV_ = PV;
	PVCov_.ResizeTo(PVCov);
	PVCov_= PVCov;
	SV_ = A1.Vertex();
	SVCov_.ResizeTo(A1.VertexCov());
	SVCov_ = A1.VertexCov();
	Phi_Res_ = Phi_Res;

	TPTRObject_ = ThreeProngTauReco();

	useDefaultMaxIterations_ = true;
	useDefaultMaxDelta_ = true;
	useDefaultEpsilon_ = true;
	useDefaultMassConstraint_ = true;
	useFullRecoil_ = false;
}

GlobalEventFit::GlobalEventFit(TrackParticle Muon, LorentzVectorParticle A1, PTObject MET, TVector3 PV, TMatrixTSym<double> PVCov){
	isConfigured_ = false;
	isFit_ = false;
	Muon_ = Muon;
	A1_= A1;
	PV_ = PV;
	PVCov_.ResizeTo(PVCov);
	PVCov_= PVCov;
	SV_ = A1.Vertex();
	SVCov_.ResizeTo(A1.VertexCov());
	SVCov_ = A1.VertexCov();
	MET_ = MET;

	TPTRObject_ = ThreeProngTauReco();

	useDefaultMaxIterations_ = true;
	useDefaultMaxDelta_ = true;
	useDefaultEpsilon_ = true;
	useDefaultMassConstraint_ = true;
	useFullRecoil_ = true;

}

GlobalEventFit::~GlobalEventFit(){

}

// Is called in the constructor and determines whether the hadronic tau decay is ambiguous and calculates the possible four-momenta of the taus.
TPTRObject GlobalEventFit::ThreeProngTauReco(){
	std::vector<LorentzVectorParticle> Taus;
	std::vector<LorentzVectorParticle> Neutrinos;
	std::vector<bool> recostatus;

	if(A1_.getParMatrix().GetNrows() != LorentzVectorParticle::NLorentzandVertexPar){
	  Logger(Logger::Error) << "A1 is not a valid LorentzVectorParticle." << std::endl;
	  return TPTRObject();
	}

	for(unsigned int Ambiguity = 0; Ambiguity<3; Ambiguity++){
		TauA1NuConstrainedFitter TauA1NU(Ambiguity,A1_,PV_,PVCov_);
		recostatus.push_back(TauA1NU.Fit());
		if(recostatus.at(Ambiguity)){
			Taus.push_back(TauA1NU.GetMother());
			LorentzVectorParticle Nu = TauA1NU.GetReFitDaughters().at(1);
			Neutrinos.push_back(Nu);
		}
		else{
			Taus.push_back(LorentzVectorParticle());
			Neutrinos.push_back(LorentzVectorParticle());
		}
	}

	bool isambiguous(IsAmbiguous(recostatus));

	TPTRObject Results = TPTRObject(A1_, Taus, Neutrinos, isambiguous, true);
	isConfigured_ = true;
	return Results;
}

// Translates the vector of ambiguity into a single boolean
bool GlobalEventFit::IsAmbiguous(std::vector<bool> recostatus){
	if(recostatus.at(0) && !recostatus.at(1) && !recostatus.at(2)) return false;
	else if (!recostatus.at(0) && recostatus.at(1) && recostatus.at(2)) return true;
	else{
		Logger(Logger::Error) << "Three prong tau reconstruction failed." << std::endl;
		return false;
	}
}

// Performs the fit for every possible tau if ambiguous. Picks solution with lowest chi2.
GEFObject GlobalEventFit::Fit(){
	if(!isConfigured_) {
		Logger(Logger::Error) << "GlobalEventFit not configured." << std::endl;
		return GEFObject();
	}
	std::vector<bool> recostatus = TPTRObject_.CreateVectorFromAmbiguity();
	std::vector<LorentzVectorParticle> Taus = TPTRObject_.getTaus();
	std::vector< std::vector<LorentzVectorParticle> > InitDaughters, RefitDaughters;
	std::vector<LorentzVectorParticle> InitResonance, FitResonance;
	std::vector<PTObject> METMinusNeutrino;
	std::vector<double> Chi2s, Csums, Niterats;
	std::vector<bool> fitstatus;

	for(unsigned Ambiguity = 0; Ambiguity<recostatus.size(); Ambiguity ++){
		if(!recostatus.at(Ambiguity)){
			fitstatus.push_back(false);
			std::vector<LorentzVectorParticle> tmp;
			InitDaughters.push_back(tmp);
			RefitDaughters.push_back(tmp);
			InitResonance.push_back(LorentzVectorParticle());
			FitResonance.push_back(LorentzVectorParticle());
			METMinusNeutrino.push_back(PTObject());
			Chi2s.push_back(-1);
			Csums.push_back(-1);
			Niterats.push_back(-1);
			continue;
		}

		METMinusNeutrino.push_back(SubtractNeutrinoFromMET(Ambiguity));

		DiTauConstrainedFitter* ptr2DTCF = NULL;
		if(useFullRecoil_){
			if(useDefaultMassConstraint_){
				ptr2DTCF = new DiTauConstrainedFitter(Taus.at(Ambiguity), Muon_, METMinusNeutrino.at(Ambiguity), PV_, PVCov_, 91.5);
				Logger(Logger::Debug) << "Case 1: ptr2DTCF->GetMassConstraint(): " << ptr2DTCF->GetMassConstraint() << std::endl;
			}
			else{
			  ptr2DTCF = new DiTauConstrainedFitter(Taus.at(Ambiguity), Muon_, METMinusNeutrino.at(Ambiguity), PV_, PVCov_, MassConstraint_);
				Logger(Logger::Debug) << "Case 2: ptr2DTCF->GetMassConstraint(): " << ptr2DTCF->GetMassConstraint() << std::endl;
			}
		}
		else{
			if(useDefaultMassConstraint_){
				ptr2DTCF = new DiTauConstrainedFitter(Taus.at(Ambiguity), Muon_, Phi_Res_, PV_, PVCov_, 91.5);
				Logger(Logger::Debug) << "Case 3: ptr2DTCF->GetMassConstraint(): " << ptr2DTCF->GetMassConstraint() << std::endl;
			}
			else{
				ptr2DTCF = new DiTauConstrainedFitter(Taus.at(Ambiguity), Muon_, Phi_Res_, PV_, PVCov_, MassConstraint_);
				Logger(Logger::Debug) << "Case 4: ptr2DTCF->GetMassConstraint(): " << ptr2DTCF->GetMassConstraint() << std::endl;
			}
		}

		InitDaughters.push_back(ptr2DTCF->GetInitialDaughters());

		if(!useDefaultMaxDelta_) ptr2DTCF->SetMaxDelta(MaxDelta_);
		if(!useDefaultMaxIterations_) ptr2DTCF->SetNIterMax(MaxIterations_);
		if(!useDefaultEpsilon_) ptr2DTCF->SetEpsilon(Epsilon_);

		fitstatus.push_back(ptr2DTCF->Fit());
		if(fitstatus.at(Ambiguity) && ptr2DTCF->isConverged()){
			FitResonance.push_back(ptr2DTCF->GetMother());
			InitResonance.push_back(ptr2DTCF->GetInitMother()); //TODO: implementation and calculation of initial resonance inside DiTauConstrainedFitter
			RefitDaughters.push_back(ptr2DTCF->GetReFitDaughters());
			Chi2s.push_back(ptr2DTCF->ChiSquare());
			Niterats.push_back(ptr2DTCF->NIter());
			Csums.push_back(ptr2DTCF->CSum());
		}
		else{
			fitstatus.push_back(false);
			std::vector<LorentzVectorParticle> tmp;
			RefitDaughters.push_back(tmp);
			InitResonance.push_back(LorentzVectorParticle());
			FitResonance.push_back(LorentzVectorParticle());
			Chi2s.push_back(-1);
			Csums.push_back(-1);
			Niterats.push_back(-1);
		}
		delete ptr2DTCF;
	}
	int IndexToReturn(-1);
	if(AmbiguitySolverByChi2(recostatus, fitstatus, Chi2s, IndexToReturn)){
		GEFObject Results(InitDaughters.at(IndexToReturn),
			InitResonance.at(IndexToReturn),
			RefitDaughters.at(IndexToReturn),
			FitResonance.at(IndexToReturn),
			true, Chi2s.at(IndexToReturn), Csums.at(IndexToReturn), Niterats.at(IndexToReturn), IndexToReturn);
		isFit_ = true;
		return Results;
	}
	else{
		Logger(Logger::Verbose) << "Fit failed: Ambiguity was not solvable" << std::endl;
		return GEFObject();
	}
}

// Solves ambiguity by chi2
bool GlobalEventFit::AmbiguitySolverByChi2(std::vector<bool> A1Fit, std::vector<bool> EventFit, std::vector<double> Chi2s, int &IndexToReturn){

	if(EventFit.at(0) == true && EventFit.at(1) == false && EventFit.at(2) == false && Chi2s.at(0) > 0){IndexToReturn =0; return true;}
	if(EventFit.at(1) == true && EventFit.at(2) == false && Chi2s.at(1) > 0){ IndexToReturn = 1;return true;}
	if(EventFit.at(1) == false && EventFit.at(2) == true && Chi2s.at(2) > 0){ IndexToReturn = 2;return true;}

	if((A1Fit.at(1) == true && A1Fit.at(2) == true) && (EventFit.at(1) == true && EventFit.at(2) == true)){
		if(Chi2s.at(1) >= 0 && Chi2s.at(2) < 0){ IndexToReturn = 1;return true;}
		else if(Chi2s.at(1) < 0 && Chi2s.at(2) >= 0){ IndexToReturn = 2;return true;}
		else if(Chi2s.at(1) >= 0 && Chi2s.at(2) >= 0){
			if(Chi2s.at(1) < Chi2s.at(2)){ IndexToReturn = 1;return true;}
			if(Chi2s.at(1) > Chi2s.at(2)){ IndexToReturn = 2;return true;}
		}
	}
	return false;
}

PTObject GlobalEventFit::SubtractNeutrinoFromMET(unsigned Ambiguity){
  TMatrixT<double> METMinusNeutrinoPar; METMinusNeutrinoPar.ResizeTo(2,1);
  TMatrixTSym<double> METMinusNeutrinoCov; METMinusNeutrinoCov.ResizeTo(2,2);

  for(unsigned i=0; i<METMinusNeutrinoCov.GetNrows(); i++){
	METMinusNeutrinoPar(i,0) = MET_.Par()(i,0) - TPTRObject_.getNeutrinos().at(Ambiguity).getParMatrix()(i,0);
	Logger(Logger::Debug) << "METMinusNeutrinoPar(" << i << ",0) " << METMinusNeutrinoPar(i,0) << std::endl;
	for(unsigned j=0; j<METMinusNeutrinoCov.GetNcols(); j++){
	  METMinusNeutrinoCov(i,j) = MET_.Cov()(i,j) + TPTRObject_.getNeutrinos().at(Ambiguity).getCovMatrix()(i,j);
	  Logger(Logger::Debug) << "METMinusNeutrinoCov(" << i << "," << j << ") " << METMinusNeutrinoCov(i,j) << std::endl;
	}
  }

  PTObject METminusNeutrino(METMinusNeutrinoPar, METMinusNeutrinoCov);

  return METminusNeutrino;
}
