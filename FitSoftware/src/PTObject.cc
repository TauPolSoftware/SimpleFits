/*
 * PTObject.cc
 *
 *  Created on: Jun 10, 2015
 *      Author: zotz
 */
#include "SimpleFits/FitSoftware/interface/PTObject.h"

  PTObject::PTObject(){
	TMatrixT<double> par; par.ResizeTo(2,1);
	par_.ResizeTo(par);
	par_=par;
	TMatrixTSym<double> cov; cov.ResizeTo(2,2);
	cov_.ResizeTo(cov);
	cov_=cov;
  }

  PTObject::PTObject(const PTObject& other){
	par_.ResizeTo(other.Par());
	par_ = other.Par();
	cov_.ResizeTo(other.Cov());
	cov_ = other.Cov();
  }

  PTObject::PTObject(TMatrixT<double> par, TMatrixTSym<double> cov){
	par_.ResizeTo(par);
	par_=par;
	cov_.ResizeTo(cov);
	cov_=cov;
  }
