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

#include "TauPolSoftware/SimpleFits/interface/GlobalEventFit.h"
#include "TauPolSoftware/SimpleFits/interface/TauA1NuConstrainedFitter.h"
#include "TauPolSoftware/SimpleFits/interface/DiTauConstrainedFitter.h"
#include "TauPolSoftware/SimpleFits/interface/ThreeProngThreeProngFitter.h"

GlobalEventFit::GlobalEventFit(TrackParticle Muon, LorentzVectorParticle A1, double Phi_Res, TVector3 PV, TMatrixTSym<double> PVCov){
	Configure(Muon, A1, PV, PVCov);
	Phi_Res_ = Phi_Res;
	useFullRecoil_ = false;
}

GlobalEventFit::GlobalEventFit(TrackParticle Muon, LorentzVectorParticle A1, PTObject MET, TVector3 PV, TMatrixTSym<double> PVCov){
	Configure(Muon, A1, PV, PVCov);
	MET_ = MET;
	useFullRecoil_ = true;
}

GlobalEventFit::GlobalEventFit(std::vector<LorentzVectorParticle> A1s, PTObject MET, TVector3 PV, TMatrixTSym<double> PVCov){
	Configure(A1s, PV, PVCov);
	MET_ = MET;
	useFullRecoil_ = true;
}

GlobalEventFit::~GlobalEventFit(){

}

void GlobalEventFit::Configure(TrackParticle Muon, LorentzVectorParticle A1, TVector3 PV, TMatrixTSym<double> PVCov){
	isConfigured_ = false;
	isFit_ = false;
	isValid_ = false;
	Muon_ = Muon;
	A1_= A1;
	A1s_.push_back(A1);
	PV_ = PV;
	PVCov_.ResizeTo(PVCov);
	PVCov_= PVCov;
	SV_ = A1.Vertex();
	SVCov_.ResizeTo(A1.VertexCov());
	SVCov_ = A1.VertexCov();
	SVs_.push_back(A1.Vertex());
	SVCovs_.push_back(A1.VertexCov());

	// Logger::Instance()->SetLevel(Logger::level::Debug);

	ThreeProngTauReconstruction();
	if (isConfigured_) TPTRObject_ = TPTRObjects_.at(0);

	useDefaultMassConstraint_ = true;
	correctPt_ = false;
	useCollinearityTauMu_ = false;
}

void GlobalEventFit::Configure(std::vector<LorentzVectorParticle> A1s, TVector3 PV, TMatrixTSym<double> PVCov){
	isConfigured_ = false;
	isFit_ = false;
	isValid_ = false;
	A1_= A1s.at(0); // for backwards compatibility
	A1s_= A1s;
	PV_ = PV;
	PVCov_.ResizeTo(PVCov);
	PVCov_= PVCov;
	SV_ = A1s.at(0).Vertex();  // for backwards compatibility
	SVCov_.ResizeTo(A1s.at(0).VertexCov());
	SVCov_ = A1s.at(0).VertexCov();  // for backwards compatibility

	for(unsigned i=0; i<A1s_.size(); i++){
		SVs_.push_back(A1s.at(i).Vertex());
		TMatrixTSym<double> SVCov;
		SVCov.ResizeTo(A1s.at(i).VertexCov());
		SVCovs_.push_back(SVCov);
	}
	// Logger::Instance()->SetLevel(Logger::level::Debug);

	ThreeProngTauReconstruction();
	if (isConfigured_) TPTRObject_ = TPTRObjects_.at(0);

	useDefaultMassConstraint_ = true;
	correctPt_ = false;
	useCollinearityTauMu_ = false;
}


// Is called in the constructor and determines whether the hadronic tau decay is ambiguous and calculates the possible four-momenta of the taus.
void GlobalEventFit::ThreeProngTauReconstruction(){
	for(unsigned i=0; i<A1s_.size(); i++){
		std::vector<LorentzVectorParticle> Taus;
		std::vector<LorentzVectorParticle> Neutrinos;
		double RotationSignificance = 0.0;
		std::vector<bool> recostatus;

		if(A1s_.at(i).getParMatrix().GetNrows() != LorentzVectorParticle::NLorentzandVertexPar){
			Logger(Logger::Error) << "A1 is not a valid LorentzVectorParticle." << std::endl;
			TPTRObjects_.emplace_back(TPTRObject());
			return;
		}

		for(unsigned int Ambiguity = 0; Ambiguity<3; Ambiguity++){
			TauA1NuConstrainedFitter TauA1NU(Ambiguity,A1s_.at(i),PV_,PVCov_);
			recostatus.push_back(TauA1NU.Fit());
			if(recostatus.at(Ambiguity)){
				Logger(Logger::Debug) << "Ambiguity: " << Ambiguity << std::endl;
				Taus.push_back(TauA1NU.GetMother());
				LorentzVectorParticle Nu = TauA1NU.GetReFitDaughters().at(1);
				Neutrinos.push_back(Nu);
				if(Ambiguity == MultiProngTauSolver::zero){
				  RotationSignificance = TauA1NU.GetTauRotationSignificance();
					Logger(Logger::Debug) << "RotationSignificance: " << RotationSignificance << std::endl;
				}
				Logger(Logger::Debug) << "Tau par and covariance: " << std::endl;
				if(Logger::Instance()->Level() == Logger::Debug){
				  Taus.at(Ambiguity).getParMatrix().Print();
				  Taus.at(Ambiguity).getCovMatrix().Print();
				}
				Logger(Logger::Debug) << "Neutrino par and covariance: " << std::endl;
				if(Logger::Instance()->Level() == Logger::Debug){
				  Neutrinos.at(Ambiguity).getParMatrix().Print();
				  Neutrinos.at(Ambiguity).getCovMatrix().Print();
				}
			}
			else{
				Taus.push_back(LorentzVectorParticle());
				Neutrinos.push_back(LorentzVectorParticle());
			}
		}
		bool isambiguous(IsAmbiguous(recostatus));
		TPTRObjects_.emplace_back(TPTRObject(A1s_.at(i), Taus, Neutrinos, isambiguous, RotationSignificance, true));
	}
	isConfigured_ = true;
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

// Performs the fit for every possible tau if ambiguous. Picks solution with lowest chi2 > 0.
GEFObject GlobalEventFit::Fit(){
	if(!isConfigured_) {
		Logger(Logger::Error) << "GlobalEventFit is not configured." << std::endl;
		return GEFObject();
	}
	std::vector< std::vector<bool> > recostatus;
	std::vector< std::vector<LorentzVectorParticle> > Taus;

	std::vector< std::vector<LorentzVectorParticle> > InitDaughters, RefitDaughters;
	std::vector<LorentzVectorParticle> InitResonance, FitResonance;
	std::vector<PTObject> METMinusNeutrino;
	std::vector<TVectorD> Chi2Vecs;
	std::vector<double> Chi2s, Csums, Niterats;
	std::vector<bool> fitstatus;
	std::vector<bool> fitvalid;

	if(A1s_.size() == 1){
		recostatus.push_back(TPTRObject_.CreateVectorFromAmbiguity());
		Taus.push_back(TPTRObject_.getTaus());

		DiTauConstrainedFitter* ptr2Fitter = NULL;

		for(unsigned Ambiguity = 0; Ambiguity<recostatus.at(0).size(); Ambiguity ++){
			if(!recostatus.at(0).at(Ambiguity)){
				fitstatus.push_back(false);
				std::vector<LorentzVectorParticle> tmp;
				for(unsigned i=0; i<2; i++) tmp.push_back(LorentzVectorParticle());
				InitDaughters.push_back(tmp);
				RefitDaughters.push_back(tmp);
				InitResonance.push_back(LorentzVectorParticle());
				FitResonance.push_back(LorentzVectorParticle());
				METMinusNeutrino.push_back(PTObject());
				Chi2Vecs.push_back(TVectorD());
				Chi2s.push_back(-1);
				Csums.push_back(-1);
				Niterats.push_back(-1);
				continue;
			}

			if(useFullRecoil_){
				PTObject ZPtEst(MET_);
				AddA1(ZPtEst);
				AddMuon(ZPtEst);
				METMinusNeutrino.push_back(ZPtEst);
				//METMinusNeutrino.push_back(SubtractNeutrinoFromMET(Ambiguity));
				if(useDefaultMassConstraint_){
					ptr2Fitter = new DiTauConstrainedFitter(Taus.at(0).at(Ambiguity), A1_, Muon_, METMinusNeutrino.at(Ambiguity), PV_, PVCov_, 91.5);
					Logger(Logger::Debug) << "Case 1: With Recoil, Default MassConstraint: " << ptr2Fitter->GetMassConstraint() << std::endl;
				}
				else{
				  ptr2Fitter = new DiTauConstrainedFitter(Taus.at(0).at(Ambiguity), A1_, Muon_, METMinusNeutrino.at(Ambiguity), PV_, PVCov_, MassConstraint_);
					Logger(Logger::Debug) << "Case 2: With Recoil, User MassConstraint: " << ptr2Fitter->GetMassConstraint() << std::endl;
				}
			}
			else{
				METMinusNeutrino.push_back(PTObject());
				if(useDefaultMassConstraint_){
					ptr2Fitter = new DiTauConstrainedFitter(Taus.at(0).at(Ambiguity), A1_, Muon_, Phi_Res_, PV_, PVCov_, 91.5);
					Logger(Logger::Debug) << "Case 3: No Recoil, Default MassConstraint: " << ptr2Fitter->GetMassConstraint() << std::endl;
				}
				else{
					ptr2Fitter = new DiTauConstrainedFitter(Taus.at(0).at(Ambiguity), A1_, Muon_, Phi_Res_, PV_, PVCov_, MassConstraint_);
					Logger(Logger::Debug) << "Case 4: No Recoil, User MassConstraint: " << ptr2Fitter->GetMassConstraint() << std::endl;
				}
			}
			if(useCollinearityTauMu_){
				ptr2Fitter->SetUseCollinearityTauMu(true);
			}

			InitDaughters.push_back(ptr2Fitter->GetInitialDaughters());
			InitResonance.push_back(ptr2Fitter->GetInitMother()); //TODO: implementation and calculation of initial resonance inside DiTauConstrainedFitter

			fitvalid.push_back(ptr2Fitter->Fit());
			isValid_ = isValid_ || fitvalid.back();
			// fitstatus.push_back(fitvalid.back() && ptr2Fitter->isConverged());
			fitstatus.push_back(fitvalid.back());
			if(fitvalid.back()){
				FitResonance.push_back(ptr2Fitter->GetMother());
				RefitDaughters.push_back(ptr2Fitter->GetReFitDaughters());
				Chi2Vecs.push_back(ptr2Fitter->ChiSquareVector());
				Chi2s.push_back(ptr2Fitter->ChiSquare());
				Niterats.push_back(ptr2Fitter->NIter());
				Csums.push_back(ptr2Fitter->CSum());
				FitPar_.ResizeTo(ptr2Fitter->GetExppar()); FitPar_ = ptr2Fitter->GetExppar();
				FitCov_.ResizeTo(ptr2Fitter->GetExpcov()); FitCov_ = ptr2Fitter->GetExpcov();
			}
			else{
				std::vector<LorentzVectorParticle> tmp;
				for(unsigned i=0; i<2; i++) tmp.push_back(LorentzVectorParticle());
				RefitDaughters.push_back(tmp);
				FitResonance.push_back(LorentzVectorParticle());
				Chi2Vecs.push_back(TVectorD());
				Chi2s.push_back(-1);
				Csums.push_back(-1);
				Niterats.push_back(-1);
			}
			delete ptr2Fitter;
		}
	}
	else if(A1s_.size() == 2){
		for (size_t i = 0; i < A1s_.size(); i++) {
			recostatus.push_back(TPTRObjects_.at(i).CreateVectorFromAmbiguity());
			Taus.push_back(TPTRObjects_.at(i).getTaus());
		}
		// TODO finish the  recostatus and Taus stuff

		ThreeProngThreeProngFitter* ptr2Fitter = NULL;

		for(unsigned AmbiguityTau1 = 0; AmbiguityTau1 < recostatus.at(0).size(); AmbiguityTau1++){
			for(unsigned AmbiguityTau2 = 0; AmbiguityTau2 < recostatus.at(1).size(); AmbiguityTau2++){
				if(!recostatus.at(0).at(AmbiguityTau1) || !recostatus.at(1).at(AmbiguityTau2)){
					fitstatus.push_back(false);
					std::vector<LorentzVectorParticle> tmp;
					for(unsigned i=0; i<2; i++) tmp.push_back(LorentzVectorParticle());
					InitDaughters.push_back(tmp);
					RefitDaughters.push_back(tmp);
					InitResonance.push_back(LorentzVectorParticle());
					FitResonance.push_back(LorentzVectorParticle());
					METMinusNeutrino.push_back(PTObject());
					Chi2Vecs.push_back(TVectorD());
					Chi2s.push_back(-1);
					Csums.push_back(-1);
					Niterats.push_back(-1);
					continue;
				}

				std::vector<LorentzVectorParticle> TauThreeProngs {Taus.at(0).at(AmbiguityTau1), Taus.at(1).at(AmbiguityTau2)};
				if(useFullRecoil_){
					PTObject ResPtEst(MET_);
					AddA1s(ResPtEst);
					METMinusNeutrino.push_back(ResPtEst);
					//METMinusNeutrino.push_back(SubtractNeutrinoFromMET(AmbiguityTau1));
					if(useDefaultMassConstraint_){
						ptr2Fitter = new ThreeProngThreeProngFitter(TauThreeProngs, A1s_, ResPtEst, PV_, PVCov_);
						Logger(Logger::Debug) << "Case 1: With Recoil, Default MassConstraint: " << ptr2Fitter->GetMassConstraint() << std::endl;
					}
					else{
					  ptr2Fitter = new ThreeProngThreeProngFitter(TauThreeProngs, A1s_, ResPtEst, PV_, PVCov_);
						ptr2Fitter->SetMassConstraint(MassConstraint_);
						Logger(Logger::Debug) << "Case 2: With Recoil, User MassConstraint: " << ptr2Fitter->GetMassConstraint() << std::endl;
					}
				}
				else{
					METMinusNeutrino.push_back(PTObject());
					if(useDefaultMassConstraint_){
						ptr2Fitter = new ThreeProngThreeProngFitter(TauThreeProngs, A1s_, PV_, PVCov_);
						Logger(Logger::Debug) << "Case 3: No Recoil, Default MassConstraint: " << ptr2Fitter->GetMassConstraint() << std::endl;
					}
					else{
						ptr2Fitter = new ThreeProngThreeProngFitter(TauThreeProngs, A1s_, PV_, PVCov_);
						ptr2Fitter->SetMassConstraint(MassConstraint_);
						Logger(Logger::Debug) << "Case 4: No Recoil, User MassConstraint: " << ptr2Fitter->GetMassConstraint() << std::endl;
					}
				}
				if(useCollinearityTauMu_){
					ptr2Fitter->SetUseCollinearityTauMu(true);
				}

				ptr2Fitter->SetFittingMode(minimizer_);

				InitDaughters.push_back(ptr2Fitter->GetInitialDaughters());
				InitResonance.push_back(ptr2Fitter->GetInitMother());

				fitvalid.push_back(ptr2Fitter->Fit());
				isValid_ = isValid_ || fitvalid.back();
				// fitstatus.push_back(fitvalid.back() && ptr2Fitter->isConverged());
				fitstatus.push_back(fitvalid.back());
				if(fitvalid.back()){
					FitResonance.push_back(ptr2Fitter->GetMother());
					RefitDaughters.push_back(ptr2Fitter->GetReFitDaughters());
					Chi2Vecs.push_back(ptr2Fitter->ChiSquareVector());
					Chi2s.push_back(ptr2Fitter->ChiSquare());
					Niterats.push_back(ptr2Fitter->NIter());
					Csums.push_back(ptr2Fitter->CSum());
					FitPar_.ResizeTo(ptr2Fitter->GetExppar()); FitPar_ = ptr2Fitter->GetExppar();
					FitCov_.ResizeTo(ptr2Fitter->GetExpcov()); FitCov_ = ptr2Fitter->GetExpcov();
				}
				else{
					std::vector<LorentzVectorParticle> tmp;
					for(unsigned i=0; i<2; i++) tmp.emplace_back(LorentzVectorParticle());
					RefitDaughters.push_back(tmp);
					// InitResonance.push_back(LorentzVectorParticle());
					FitResonance.emplace_back(LorentzVectorParticle());
					Chi2Vecs.emplace_back(TVectorD());
					Chi2s.push_back(-1);
					Csums.push_back(-1);
					Niterats.push_back(-1);
				}
				delete ptr2Fitter;
			}
		}
	}

	fitstatuses_ = fitstatus;
	int IndexToReturn(-1);
	bool foundSolution(false);

	if(minimizer_ == LagrangeMultipliersFitter::FittingProc::Standard)
		foundSolution = AmbiguitySolverByChi2(recostatus, fitstatus, Chi2Vecs, IndexToReturn);
	else if(minimizer_ == LagrangeMultipliersFitter::FittingProc::Minuit)
		foundSolution = AmbiguitySolverByChi2Minuit(fitstatus, Chi2s, IndexToReturn);

	if(foundSolution){
		// std::vector<LorentzVectorParticle> CorrFitDaughters = FitDaughtersCorr(RefitDaughters.at(IndexToReturn.at(0)));
		// RefitDaughters.at(IndexToReturn.at(0)).at(0) = CorrFitDaughters.at(0);
		// RefitDaughters.at(IndexToReturn.at(0)).at(1) = CorrFitDaughters.at(1);
		isFit_ = true;
	}
	else{
		Logger(Logger::Verbose) << "AmbiguitySolver failed: Fit did not converge." << std::endl;
	}
	// Logger(Logger::Info) << "IndexToReturn: " << IndexToReturn << std::endl;
	return GEFObject(InitDaughters,
		InitResonance,
		RefitDaughters,
		FitResonance,
		isValid_, foundSolution, Chi2Vecs, Csums, Niterats, IndexToReturn);
}

// Solves ambiguity by chi2
bool GlobalEventFit::AmbiguitySolverByChi2(std::vector< std::vector<bool> > A1Fit, std::vector<bool> EventFit, std::vector<TVectorD> Chi2Vecs, int &IndexToReturn){
	if(A1s_.size() == 1){
		if(EventFit.at(0) == true && EventFit.at(1) == false && EventFit.at(2) == false && Chi2Vecs.at(0).GetNoElements() == 3){IndexToReturn =0; return true;}
		if(EventFit.at(1) == true && EventFit.at(2) == false && Chi2Vecs.at(1).GetNoElements() == 3){ IndexToReturn = 1;return true;}
		if(EventFit.at(1) == false && EventFit.at(2) == true && Chi2Vecs.at(2).GetNoElements() == 3){ IndexToReturn = 2;return true;}

		if((A1Fit.at(0).at(1) == true && A1Fit.at(0).at(2) == true) && (EventFit.at(1) == true && EventFit.at(2) == true)){

		  /*
		  if(Chi2s.at(1) >= 0 && Chi2s.at(2) < 0){ IndexToReturn = 1;return true;}
		  else if(Chi2s.at(1) < 0 && Chi2s.at(2) >= 0){ IndexToReturn = 2;return true;}
		  else if(Chi2s.at(1) >= 0 && Chi2s.at(2) >= 0){
		    if(Chi2s.at(1) < Chi2s.at(2)){ IndexToReturn = 1;return true;}
		    if(Chi2s.at(1) > Chi2s.at(2)){ IndexToReturn = 2;return true;}
		  }
		  */
		    if(Chi2Vecs.at(1).Sum() < Chi2Vecs.at(2).Sum()){ IndexToReturn = 1;return true;}
		    if(Chi2Vecs.at(1).Sum() > Chi2Vecs.at(2).Sum()){ IndexToReturn = 2;return true;}
		}
	}
	else if(A1s_.size() == 2){
		bool found_solution(false);
		double Chi2Min(9999);
		for(size_t i = 0; i < EventFit.size(); i++) {
			// Logger(Logger::Info) << "EventFit.at(" << i << "): " << EventFit.at(i) << std::endl;
			if(EventFit.at(i) && (Chi2Vecs.at(i).Sum() < Chi2Min)){
				Chi2Min = Chi2Vecs.at(i).Sum();
				IndexToReturn = i;
				found_solution = true;
			}
		}
		return found_solution;
	}
	// if((A1Fit.at(0).at(1) == true && A1Fit.at(0).at(2) == true) && (EventFit.at(1) == false && EventFit.at(2) == false)){
	//     if(Chi2Vecs.at(1).Sum() < Chi2Vecs.at(2).Sum()){ IndexToReturn = 1;return false;}
	//     if(Chi2Vecs.at(1).Sum() > Chi2Vecs.at(2).Sum()){ IndexToReturn = 2;return false;}
	// }
	return false;
}

// Solves ambiguity by chi2 for Minuit based minimization
bool GlobalEventFit::AmbiguitySolverByChi2Minuit(std::vector<bool> EventFit, std::vector<double> Chi2s, int &IndexToReturn){
	// Logger(Logger::Info) << "A1Fit.at(0): " << A1Fit.at(0) << std::endl;
	// Logger(Logger::Info) << "A1Fit.at(1): " << A1Fit.at(1) << std::endl;
	// Logger(Logger::Info) << "A1Fit.at(2): " << A1Fit.at(2) << std::endl;
	//
	// Logger(Logger::Info) << "EventFit.at(0): " << EventFit.at(0) << std::endl;
	// Logger(Logger::Info) << "EventFit.at(1): " << EventFit.at(1) << std::endl;
	// Logger(Logger::Info) << "EventFit.at(2): " << EventFit.at(2) << std::endl;
	if(A1s_.size() == 1){
		if(EventFit.at(0) == true && EventFit.at(1) == false && EventFit.at(2) == false){
			IndexToReturn = 0;
			return true;
		}
		else if(EventFit.at(1) == true && EventFit.at(2) == false){
			IndexToReturn = 1;
			return true;
		}
		else if(EventFit.at(1) == false && EventFit.at(2) == true){
			IndexToReturn = 2;
			return true;
		}
		else if(EventFit.at(1) == true && EventFit.at(2) == true){
			// Logger(Logger::Info) << "Chi2s.at(1): " << Chi2s.at(1) << std::endl;
			// Logger(Logger::Info) << "Chi2s.at(2): " << Chi2s.at(2) << std::endl;
			if(Chi2s.at(1) < Chi2s.at(2)){ IndexToReturn = 1;return true;}
			if(Chi2s.at(1) > Chi2s.at(2)){ IndexToReturn = 2;return true;}
		}
	}
	else if(A1s_.size() == 2){
		bool found_solution(false);
		double Chi2Min(9999);
		for(size_t i = 0; i < EventFit.size(); i++) {
			// Logger(Logger::Info) << "EventFit.at(" << i << "): " << EventFit.at(i) << std::endl;
			if(EventFit.at(i) && (Chi2s.at(i) < Chi2Min)){
				Chi2Min = Chi2s.at(i);
				IndexToReturn = i;
				found_solution = true;
			}
		}
		return found_solution;
	}
	return false;
}


std::vector<LorentzVectorParticle> GlobalEventFit::FitDaughtersCorr(std::vector<LorentzVectorParticle> FitDaughters){

  std::vector<LorentzVectorParticle> out;
  //  constants are taken from pol1 fit to profile of pt resolution vs post fit pt:

  //***************************************************************
  //--------------------------without recoil-----------------------
  //***************************************************************
  // **********Tau A1 ******************************
  // Minimizer is Linear
  // Chi2                      =      149.414
  // NDf                       =           48
  // p0                        =      35.9909   +/-   0.725711    
  // p1                        =    -0.833072   +/-   0.0172851  


  // **********Tau Mu ******************************
  // Minimizer is Linear
  // Chi2                      =      76.3996
  // NDf                       =           48
  // p0                        =      44.9105   +/-   0.728919    
  // p1                        =    -0.967167   +/-   0.0186608   


  //***************************************************************
  //--------------------------with recoil--------------------------
  //***************************************************************
  // **********Tau Mu ******************************
  // **********without ambiguity
  //****************************************
  //Minimizer is Linear
  //Chi2                      =      49.1458
  //NDf                       =           26
  //p0                        =      -40.692   +/-   0.388669
  //p1                        =     0.937738   +/-   0.0116998
  //
  // **********with ambiguity
  //****************************************
  //Minimizer is Linear
  //Chi2                      =      23.6312
  //NDf                       =           29
  //p0                        =     -41.4081   +/-   0.455316
  //p1                        =     0.947888   +/-   0.0165346
  //
  // **********Tau A1 ******************************
  // **********without ambiguity
  //****************************************
  //Minimizer is Linear
  //Chi2                      =      14.2012
  //NDf                       =           16
  //p0                        =     -21.9364   +/-   0.648501
  //p1                        =     0.571471   +/-   0.0139919
  //
  // **********with ambiguity
  //****************************************
  //Minimizer is Linear
  //Chi2                      =      78.1695
  //NDf                       =           33
  //p0                        =     -31.3707   +/-   0.705552
  //p1                        =     0.843426   +/-   0.0121233

  TLorentzVector TauA1 =FitDaughters.at(0).LV();
  TLorentzVector TauMu =FitDaughters.at(1).LV();

  TMatrixT<double> parTauA1 = FitDaughters.at(0).getParMatrix();
  TMatrixT<double> parTauMu = FitDaughters.at(1).getParMatrix();

  if(correctPt_){
	double p0tmu, p1tmu, p0ta1, p1ta1;
	if(!useFullRecoil_){
	  p0tmu = 44.9105;
	  p1tmu =-0.967167;

	  p0ta1 = 35.9909;
	  p1ta1 =-0.833072;
	}
	else{
	  if(TPTRObject_.isAmbiguous()){
		p0tmu = 41.4081;
		p1tmu = -0.947888;

		p0ta1 = 31.3707;
		p1ta1 = -0.843426;
	  }
	  else{
		p0tmu = 40.692;
		p1tmu = -0.937738;

		p0ta1 = 21.9364;
		p1ta1 = -0.571471;
	  }
	}
	TauA1.SetPerp(TauA1.Perp()*(1+ (p1ta1*TauA1.Perp() + p0ta1)/TauA1.Perp()));
	TauMu.SetPerp(TauMu.Perp()*(1+ (p1tmu*TauMu.Perp() + p0tmu)/TauMu.Perp()));
  }

  // -----  correct only LorenzVector  for now;   TODO: momentum covariance to be corrected later;
  parTauA1(LorentzVectorParticle::px,0) = TauA1.Px();
  parTauA1(LorentzVectorParticle::py,0) = TauA1.Py();
  parTauA1(LorentzVectorParticle::pz,0) = TauA1.Pz();
  parTauA1(LorentzVectorParticle::m,0)  = 1.777;

  parTauMu(LorentzVectorParticle::px,0) = TauMu.Px();
  parTauMu(LorentzVectorParticle::py,0) = TauMu.Py();
  parTauMu(LorentzVectorParticle::pz,0) = TauMu.Pz();
  parTauMu(LorentzVectorParticle::m,0)  = 1.777;

  out.push_back(LorentzVectorParticle(parTauA1,FitDaughters.at(0).getCovMatrix(),FitDaughters.at(0).PDGID(),FitDaughters.at(0).Charge(),FitDaughters.at(0).BField()));
  out.push_back(LorentzVectorParticle(parTauMu,FitDaughters.at(1).getCovMatrix(),FitDaughters.at(1).PDGID(),FitDaughters.at(1).Charge(),FitDaughters.at(1).BField()));

  return out;
}



PTObject GlobalEventFit::SubtractNeutrinoFromMET(unsigned Ambiguity){
  TMatrixT<double> METMinusNeutrinoPar; METMinusNeutrinoPar.ResizeTo(2,1);
  TMatrixTSym<double> METMinusNeutrinoCov; METMinusNeutrinoCov.ResizeTo(2,2);

  for(int i=0; i<METMinusNeutrinoCov.GetNrows(); i++){
	METMinusNeutrinoPar(i,0) = MET_.Par()(i,0) - TPTRObject_.getNeutrinos().at(Ambiguity).getParMatrix()(i,0);
	Logger(Logger::Debug) << "METMinusNeutrinoPar(" << i << ",0) " << METMinusNeutrinoPar(i,0) << std::endl;
	for(int j=0; j<METMinusNeutrinoCov.GetNcols(); j++){
	  METMinusNeutrinoCov(i,j) = MET_.Cov()(i,j) + TPTRObject_.getNeutrinos().at(Ambiguity).getCovMatrix()(i+3,j+3);
	  Logger(Logger::Debug) << "METMinusNeutrinoCov(" << i << "," << j << ") " << METMinusNeutrinoCov(i,j) << std::endl;
	}
  }

  Logger(Logger::Debug) << "MET covariance: " << std::endl;
  if(Logger::Instance()->Level() == Logger::Debug){
	MET_.Cov().Print();
  }
  Logger(Logger::Debug) << "Neutrino covariance: " << std::endl;
  if(Logger::Instance()->Level() == Logger::Debug){
	TPTRObject_.getNeutrinos().at(Ambiguity).getCovMatrix().Print();
  }
  Logger(Logger::Debug) << "METMinusNeutrinoCov covariance: " << std::endl;
  if(Logger::Instance()->Level() == Logger::Debug){
	METMinusNeutrinoCov.Print();
  }
  PTObject METminusNeutrino(METMinusNeutrinoPar, METMinusNeutrinoCov);

  return METminusNeutrino;
}

PTObject GlobalEventFit::AddA1(PTObject MET){
  TMatrixT<double> METPlusA1Par; METPlusA1Par.ResizeTo(2,1);
  TMatrixTSym<double> METPlusA1Cov; METPlusA1Cov.ResizeTo(2,2);

  for(int i=0; i<METPlusA1Cov.GetNrows(); i++){
	METPlusA1Par(i,0) = MET.Par()(i,0) + TPTRObject_.getA1().getParMatrix()(i,0);
	Logger(Logger::Debug) << "METMinusNeutrinoPar(" << i << ",0) " << METPlusA1Par(i,0) << std::endl;
	for(int j=0; j<METPlusA1Cov.GetNcols(); j++){
	  METPlusA1Cov(i,j) = MET.Cov()(i,j) + TPTRObject_.getA1().getCovMatrix()(i+3,j+3);
	  Logger(Logger::Debug) << "METMinusNeutrinoCov(" << i << "," << j << ") " << METPlusA1Cov(i,j) << std::endl;
	}
  }

  Logger(Logger::Debug) << "MET covariance: " << std::endl;
  if(Logger::Instance()->Level() == Logger::Debug){
	MET.Cov().Print();
  }
  Logger(Logger::Debug) << "A1 covariance: " << std::endl;
  if(Logger::Instance()->Level() == Logger::Debug){
	TPTRObject_.getA1().getCovMatrix().Print();
  }
  Logger(Logger::Debug) << "METPlusA1Cov covariance: " << std::endl;
  if(Logger::Instance()->Level() == Logger::Debug){
	METPlusA1Cov.Print();
  }
  PTObject METPlusA1(METPlusA1Par, METPlusA1Cov);

  return METPlusA1;
}

PTObject GlobalEventFit::AddA1s(PTObject MET){
  TMatrixT<double> METPlusA1Par; METPlusA1Par.ResizeTo(2,1);
  TMatrixTSym<double> METPlusA1Cov; METPlusA1Cov.ResizeTo(2,2);

  for (size_t i_a1 = 0; i_a1 < A1s_.size(); i_a1++) {
    for(int i=0; i<METPlusA1Cov.GetNrows(); i++){
      METPlusA1Par(i,0) = MET.Par()(i,0) + TPTRObjects_.at(i_a1).getA1().getParMatrix()(i,0);
      Logger(Logger::Debug) << "METMinusNeutrinoPar(" << i << ",0) " << METPlusA1Par(i,0) << std::endl;
      for(int j=0; j<METPlusA1Cov.GetNcols(); j++){
        METPlusA1Cov(i,j) = MET.Cov()(i,j) + TPTRObjects_.at(i_a1).getA1().getCovMatrix()(i+3,j+3);
        Logger(Logger::Debug) << "METMinusNeutrinoCov(" << i << "," << j << ") " << METPlusA1Cov(i,j) << std::endl;
      }
    }
  }

  Logger(Logger::Debug) << "MET covariance: " << std::endl;
  if(Logger::Instance()->Level() == Logger::Debug){
	MET.Cov().Print();
  }
  Logger(Logger::Debug) << "A1 covariance: " << std::endl;
  if(Logger::Instance()->Level() == Logger::Debug){
	TPTRObjects_.at(0).getA1().getCovMatrix().Print();
	TPTRObjects_.at(1).getA1().getCovMatrix().Print();
  }
  Logger(Logger::Debug) << "METPlusA1Cov covariance: " << std::endl;
  if(Logger::Instance()->Level() == Logger::Debug){
	METPlusA1Cov.Print();
  }
  PTObject METPlusA1(METPlusA1Par, METPlusA1Cov);

  return METPlusA1;
}

PTObject GlobalEventFit::AddMuon(PTObject MET){
  TMatrixT<double> METPlusMuonPar; METPlusMuonPar.ResizeTo(2,1);
  TMatrixTSym<double> METPlusMuonCov; METPlusMuonCov.ResizeTo(2,2);

  double kappa = Muon_.Parameter(TrackParticle::kappa);
  double alpha = Muon_.BField();
  double phi0 = Muon_.Parameter(TrackParticle::phi);
  METPlusMuonPar(0,0) = MET.Par()(0,0) + fabs(alpha/kappa)*cos(phi0);
  METPlusMuonPar(1,0) = MET.Par()(1,0) + fabs(alpha/kappa)*sin(phi0);

  // Logger(Logger::Info) << "fabs(alpha/kappa): " << fabs(alpha/kappa) << "\n";
  // Logger(Logger::Info) << "fabs(alpha/kappa)*cos(phi0): " << fabs(alpha/kappa)*cos(phi0) << "\n";
  // Logger(Logger::Info) << "fabs(alpha/kappa)*sin(phi0): " << fabs(alpha/kappa)*sin(phi0) << "\n";

  TMatrixDSym MuonKappaPhiCov(2);
  MuonKappaPhiCov(0,0) = Muon_.Covariance(TrackParticle::kappa,TrackParticle::kappa);
  MuonKappaPhiCov(0,1) = Muon_.Covariance(TrackParticle::kappa,TrackParticle::phi);
  MuonKappaPhiCov(1,0) = MuonKappaPhiCov(0,1);
  MuonKappaPhiCov(1,1) = Muon_.Covariance(TrackParticle::phi,TrackParticle::phi);

  TMatrixD Jacobi(2,2);
  Jacobi(0,0) = - alpha*kappa/fabs(pow(kappa, 3.))*cos(phi0);
  Jacobi(0,1) = - alpha/fabs(kappa)*sin(phi0);
  Jacobi(1,0) = - alpha*kappa/fabs(pow(kappa, 3.))*sin(phi0);
  Jacobi(1,1) = + alpha/fabs(kappa)*cos(phi0);

  TMatrixDSym MuonPtCov(MuonKappaPhiCov);
  MuonPtCov.Similarity(Jacobi);

  for(int i=0; i<METPlusMuonCov.GetNrows(); i++){
	for(int j=0; j<METPlusMuonCov.GetNcols(); j++){
	  METPlusMuonCov(i,j) = MET.Cov()(i,j) + MuonPtCov(i,j);
	}
  }

  Logger(Logger::Debug) << "MET covariance: " << std::endl;
  if(Logger::Instance()->Level() == Logger::Debug){
	MET.Cov().Print();
  }
  Logger(Logger::Debug) << "MuonPtCov covariance: " << std::endl;
  if(Logger::Instance()->Level() == Logger::Debug){
	MuonPtCov.Print();
  }
  Logger(Logger::Debug) << "METPlusA1Cov covariance: " << std::endl;
  if(Logger::Instance()->Level() == Logger::Debug){
	METPlusMuonCov.Print();
  }
  PTObject METPlusMuon(METPlusMuonPar, METPlusMuonCov);

  return METPlusMuon;
}
