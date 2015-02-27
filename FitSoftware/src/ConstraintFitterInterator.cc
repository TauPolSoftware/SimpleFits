#include "SimpleFits/FitSoftware/interface/ConstraintFitterInterator.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"
#include "TFile.h"

bool ConstraintFitterInterator::FitwithConstraint(){
  
  unsigned int nmax=1;
  for(unsigned int i=0;i<constraints_.size();i++){if(constraints_.at(i)->maxIterations()>nmax) nmax=constraints_.at(i)->maxIterations();}
  if(Fit()){
    for(unsigned int n=0;n<=nmax;n++){
      for(unsigned int i=0;i<constraints_.size();i++){
	if(!constraints_.at(i)->isConverged() && !constraints_.at(i)->atLimit()) constraints_.at(i)->UpdateCovariance();
      }
      if(!Fit()) return false;
    }
  }
  else return false;
  return true;
}


void ConstraintFitterInterator::AddConstraints(Constraint *_constraint){
  constraints_.push_back(_constraint);
  ConstraintConvergence_.push_back(TH1D(_constraint->name().data(),_constraint->name().data(),_constraint->maxIterations()+1,-0.5,((float)_constraint->maxIterations())+0.5));
}

void ConstraintFitterInterator::SaveConstraintConvergence(std::string _file){
  Log(Log::Verbose) << "Saving Constriants in " << _file << std::endl;
  TFile f(_file.data(),"NEW");
  for(unsigned int i=0;i<constraints_.size();i++)ConstraintConvergence_.at(i).Write(constraints_.at(i)->name().data());
  f.Close();
}
