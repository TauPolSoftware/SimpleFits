#include "SimpleFits/FitSoftware/interface/ConstraintFitterInterator.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"
#include "TFile.h"

ConstraintFitterInterator::~ConstraintFitterInterator(){
  for(unsigned int i=0;i<constraints_.size();i++){
    delete constraints_.at(i);
    constraints_.at(i)=NULL;
  }
  constraints_.clear();
}

bool ConstraintFitterInterator::Fit(){
  unsigned int nmax=1;
  for(unsigned int i=0;i<constraints_.size();i++){if(constraints_.at(i)->maxIterations()>nmax) nmax=constraints_.at(i)->maxIterations();}
  FitChi2=TH1D("FitChi2","FitChi2",nmax+1,-0.5,((float)nmax)+0.5);
  FitUncorrectedChi2=TH1D("FitUncorrectedChi2","FitUncorrectedChi2",nmax+1,-0.5,((float)nmax)+0.5);
  if(FitAndMonitor()){
    for(iteration=0;n<=nmax;n++){
      for(unsigned int i=0;i<constraints_.size();i++){
	if(!constraints_.at(i)->isConverged() && !constraints_.at(i)->atLimit()) constraints_.at(i)->UpdateCovariance();
      }
      if(!FitAndMonitor()) return false;
    }
  }
  else return false;
  return true;
}

void ConstraintFitterInterator::AddConstraint(Constraint* c){
  constraints_.push_back(c);
  ConstraintConvergence_.push_back(TH1D(c->name().data()+"_dx",c->name().data()+"_dx",,c->maxIterations()+1,-0.5,((float)c->maxIterations())+0.5));
  ConstraintSignificance_.push_back(TH1D(c->name().data()+"_sig",c->name().data()"_sig",c->maxIterations()+1,-0.5,((float)c->maxIterations())+0.5));
}

void ConstraintFitterInterator::SaveConstraintConvergence(std::string _file){
  if(monitor_){
    Log(Log::Verbose) << "Saving Constriants in " << _file << std::endl;
    TFile f(_file.data(),"NEW");
    for(unsigned int i=0;i<constraints_.size();i++){
      ConstraintConvergence_.at(i).Write(ConstraintConvergence_.at(i).Name());
      ConstraintSignificance_.at(i).Write(ConstraintSignificance_.at(i).Name());
    }
    FitChi2.Write(FitChi2.Name());
    FitUncorrectedChi2.Write(FitUncorrectedChi2.Name());
    f.Close();
  }
}

bool ConstraintFitterInterator::FitCorrectedAndMonitor(){
  FitwithConstraints();
  uncorrChi2=chi2;
  for(unsigned int i=0;i<constraints_.size();i++){if(constraints_.at(i)->isHard())chi2-=constraints_.at(i)->chi2();}
  if(monitor_){
    for(unsigned int i=0;i<constraints_.size();i++){
      ConstraintConvergence_.at(i).Fill(iteration_,constraints_.at(i)->dx());
      ConstraintSignificance_.at(i).Fill(iteration_,constraints_.at(i)->significance());
    }
    FitChi2.Fill(iteration_,chi2);
    FitUncorrectedChi2.Fill(iteration_,uncorrChi2);
  }
}

