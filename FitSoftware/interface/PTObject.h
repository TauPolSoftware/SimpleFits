/*
 * PTObject.h
 *
 *  Created on: Jun 10, 2015
 *      Author: zotz
 */

#ifndef PTObject_H_
#define PTObject_H_
#include "TMatrixT.h"
#include "TMatrixTSym.h"

class PTObject{

 public:
  PTObject();
  PTObject(const PTObject& other);
  PTObject(TMatrixT<double> par, TMatrixTSym<double> cov);
  virtual ~PTObject(){};

  double X(){return par_(0,0);}
  double Y(){return par_(1,0);}
  double Phi(){return X()== 0 && Y() == 0 ? 0 : atan2(Y(),X());}
  double Pt(){return sqrt(pow(X(),2.) + pow(Y(),2.));}
  TMatrixT<double> Par() const{return par_;}
  TMatrixTSym<double> Cov() const{return cov_;}

 private:
  TMatrixT<double> par_;
  TMatrixTSym<double> cov_;

};
#endif /* PTObject_H_ */
