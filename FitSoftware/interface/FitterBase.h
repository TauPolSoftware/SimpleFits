#ifndef SimpleFits_FitterBase_h
#define SimpleFits_FitterBase_h

// system include files
#include "TMatrixT.h"
#include "TMatrixTSym.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TVectorD.h"

class FitterBase{
 public:
  FitterBase(){};
  virtual ~FitterBase();

   virtual bool Fit()=0;
   virtual void ResetFit();
   virtual double UpdateChisquare(TMatrixT<double> inpar)=0;
   virtual double ChiSquare(){return chi2;}
   virtual double NDF(){return ndf;}

 protected:
   bool isFit,isConfigure;
   // Free parameters
   TMatrixT<double> par;
   TMatrixTSym<double> parcov;
   double chi2, ndf;

   // Measured values
   TMatrixT<double> val;
   TMatrixTSym<double> cov;
   TMatrixTSym<double> cov_inv;
   TMatrixTSym<double> dalpha;

};
#endif


