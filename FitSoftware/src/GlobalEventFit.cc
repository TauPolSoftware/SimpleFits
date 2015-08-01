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
	Configure(Muon, A1, PV, PVCov);
	Phi_Res_ = Phi_Res;
	useFullRecoil_ = false;
}

GlobalEventFit::GlobalEventFit(TrackParticle Muon, LorentzVectorParticle A1, PTObject MET, TVector3 PV, TMatrixTSym<double> PVCov){
	Configure(Muon, A1, PV, PVCov);
	MET_ = MET;
	useFullRecoil_ = true;
}

GlobalEventFit::~GlobalEventFit(){

}

void GlobalEventFit::Configure(TrackParticle Muon, LorentzVectorParticle A1, TVector3 PV, TMatrixTSym<double> PVCov){
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

	TPTRObject_ = ThreeProngTauReconstruction();

	useDefaultMassConstraint_ = true;
	correctPt_ = true;
}

// Is called in the constructor and determines whether the hadronic tau decay is ambiguous and calculates the possible four-momenta of the taus.
TPTRObject GlobalEventFit::ThreeProngTauReconstruction(){
	std::vector<LorentzVectorParticle> Taus;
	std::vector<LorentzVectorParticle> Neutrinos;
	double RotationSignificance;
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
			if(Ambiguity == MultiProngTauSolver::zero){
			  RotationSignificance = TauA1NU.GetTauRotationSignificance();
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

	TPTRObject Results = TPTRObject(A1_, Taus, Neutrinos, isambiguous, RotationSignificance, true);
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

// Performs the fit for every possible tau if ambiguous. Picks solution with lowest chi2 > 0.
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
	std::vector<TVectorD> Chi2Vecs;
	std::vector<double> Chi2s, Csums, Niterats;
	std::vector<bool> fitstatus;
	DiTauConstrainedFitter* ptr2DTCF = NULL;

	for(unsigned Ambiguity = 0; Ambiguity<recostatus.size(); Ambiguity ++){
		if(!recostatus.at(Ambiguity)){
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
			METMinusNeutrino.push_back(SubtractNeutrinoFromMET(Ambiguity));
			if(useDefaultMassConstraint_){
				ptr2DTCF = new DiTauConstrainedFitter(Taus.at(Ambiguity), Muon_, METMinusNeutrino.at(Ambiguity), PV_, PVCov_, 91.5);
				Logger(Logger::Debug) << "Case 1: With Recoil, Default MassConstraint: " << ptr2DTCF->GetMassConstraint() << std::endl;
			}
			else{
			  ptr2DTCF = new DiTauConstrainedFitter(Taus.at(Ambiguity), Muon_, METMinusNeutrino.at(Ambiguity), PV_, PVCov_, MassConstraint_);
				Logger(Logger::Debug) << "Case 2: With Recoil, User MassConstraint: " << ptr2DTCF->GetMassConstraint() << std::endl;
			}
		}
		else{
			METMinusNeutrino.push_back(PTObject());
			if(useDefaultMassConstraint_){
				ptr2DTCF = new DiTauConstrainedFitter(Taus.at(Ambiguity), Muon_, Phi_Res_, PV_, PVCov_, 91.5);
				Logger(Logger::Debug) << "Case 3: No Recoil, Default MassConstraint: " << ptr2DTCF->GetMassConstraint() << std::endl;
			}
			else{
				ptr2DTCF = new DiTauConstrainedFitter(Taus.at(Ambiguity), Muon_, Phi_Res_, PV_, PVCov_, MassConstraint_);
				Logger(Logger::Debug) << "Case 4: No Recoil, User MassConstraint: " << ptr2DTCF->GetMassConstraint() << std::endl;
			}
		}

		InitDaughters.push_back(ptr2DTCF->GetInitialDaughters());

		fitstatus.push_back(ptr2DTCF->Fit() && ptr2DTCF->isConverged());
		if(fitstatus.at(Ambiguity)){
			FitResonance.push_back(ptr2DTCF->GetMother());
			InitResonance.push_back(ptr2DTCF->GetInitMother()); //TODO: implementation and calculation of initial resonance inside DiTauConstrainedFitter
			RefitDaughters.push_back(ptr2DTCF->GetReFitDaughters());
			Chi2Vecs.push_back(ptr2DTCF->ChiSquareVector());
			Chi2s.push_back(ptr2DTCF->ChiSquare());
			Niterats.push_back(ptr2DTCF->NIter());
			Csums.push_back(ptr2DTCF->CSum());
			FitPar_.ResizeTo(ptr2DTCF->GetExppar()); FitPar_ = ptr2DTCF->GetExppar();
			FitCov_.ResizeTo(ptr2DTCF->GetExpcov()); FitCov_ = ptr2DTCF->GetExpcov();
		}
		else{
			fitstatus.push_back(false);
			std::vector<LorentzVectorParticle> tmp;
			for(unsigned i=0; i<2; i++) tmp.push_back(LorentzVectorParticle());
			RefitDaughters.push_back(tmp);
			InitResonance.push_back(LorentzVectorParticle());
			FitResonance.push_back(LorentzVectorParticle());
			Chi2Vecs.push_back(TVectorD());
			Chi2s.push_back(-1);
			Csums.push_back(-1);
			Niterats.push_back(-1);
		}
		delete ptr2DTCF;
	}
	fitstatuses_ = fitstatus;
	int IndexToReturn(-1);
	if(AmbiguitySolverByChi2(recostatus, fitstatus, Chi2Vecs, IndexToReturn)){

	  std::vector<LorentzVectorParticle> CorrFitDaughters = FitDaughtersCorr(RefitDaughters.at(IndexToReturn));
	  RefitDaughters.at(IndexToReturn).at(0) = CorrFitDaughters.at(0);
	  RefitDaughters.at(IndexToReturn).at(1) = CorrFitDaughters.at(1);
	  GEFObject Results = GEFObject(InitDaughters,
					InitResonance,
					RefitDaughters,
					FitResonance,
					true, Chi2Vecs, Csums, Niterats, IndexToReturn);

		isFit_ = true;
		return Results;
	}
	else{
	  	Logger(Logger::Verbose) << "Fit failed: Ambiguity was not solvable" << std::endl;
		return GEFObject();
	}
}

// Solves ambiguity by chi2
bool GlobalEventFit::AmbiguitySolverByChi2(std::vector<bool> A1Fit, std::vector<bool> EventFit, std::vector<TVectorD> Chi2Vecs, int &IndexToReturn){

	if(EventFit.at(0) == true && EventFit.at(1) == false && EventFit.at(2) == false && Chi2Vecs.at(0).GetNoElements() == 3){IndexToReturn =0; return true;}
	if(EventFit.at(1) == true && EventFit.at(2) == false && Chi2Vecs.at(1).GetNoElements() == 3){ IndexToReturn = 1;return true;}
	if(EventFit.at(1) == false && EventFit.at(2) == true && Chi2Vecs.at(2).GetNoElements() == 3){ IndexToReturn = 2;return true;}

	if((A1Fit.at(1) == true && A1Fit.at(2) == true) && (EventFit.at(1) == true && EventFit.at(2) == true)){

	  /*
	  if(Chi2s.at(1) >= 0 && Chi2s.at(2) < 0){ IndexToReturn = 1;return true;}
	  else if(Chi2s.at(1) < 0 && Chi2s.at(2) >= 0){ IndexToReturn = 2;return true;}
	  else if(Chi2s.at(1) >= 0 && Chi2s.at(2) >= 0){
	    if(Chi2s.at(1) < Chi2s.at(2)){ IndexToReturn = 1;return true;}
	    if(Chi2s.at(1) > Chi2s.at(2)){ IndexToReturn = 2;return true;}
	  }
	  */
	    if(Chi2Vecs.at(1)(0) < Chi2Vecs.at(2)(0)){ IndexToReturn = 1;return true;}
	    if(Chi2Vecs.at(1)(0) > Chi2Vecs.at(2)(0)){ IndexToReturn = 2;return true;}
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

  for(unsigned i=0; i<METMinusNeutrinoCov.GetNrows(); i++){
	METMinusNeutrinoPar(i,0) = MET_.Par()(i,0) - TPTRObject_.getNeutrinos().at(Ambiguity).getParMatrix()(i,0);
	Logger(Logger::Debug) << "METMinusNeutrinoPar(" << i << ",0) " << METMinusNeutrinoPar(i,0) << std::endl;
	for(unsigned j=0; j<METMinusNeutrinoCov.GetNcols(); j++){
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
//GEFObject GlobalEventFit::Refit(){
//
//}
