#ifndef SimpleFits_Constraint_h
#define SimpleFits_Constraint_h

#include <string.h>
#include "TMatrixT.h"
#include "TMatrixTSym.h"

class Constraint{
 public:
  Constraint(std::string _name, unsigned int _idx, TMatrixT<double> *_val, TMatrixTSym<double> *_cov, TMatrixT<double> *_dalpha, double _initalSigma, double _sigma_min, double _dx, unsigned int _maxIterations=5, double _scalefactor=0.75);
  virtual ~Constraint();

  inline std::string name(){return name_;}
  inline double constraint(){return ((*val_)(idx_,0));}
  inline double sigma(){return sqrt((*cov_)(idx_,idx_));}
  inline double dalpha(){return ((*dalpha_)(idx_,0));}
  inline double significance(){return ((*dalpha_)(idx_,0))/((*cov_)(idx_,idx_));}
  inline double dx(){return dx_;}
  inline double initalSigma(){return initalSigma_;}
  inline double sigma_min(){return sigma_min_;}
  inline bool   atLimit(){return (nIterations_>=maxIterations_);}
  inline unsigned int nIterations(){return nIterations_;}
  inline unsigned int maxIterations(){return maxIterations_;}
  virtual void addCorrelation(int idx){};// to be defined in the the inheriting Constraint class if correlations can be added.

  // important function for modifying the cov to iteratively increase the constriants
  virtual void UpdateCovariance()=0;
  virtual bool isConverged()=0;

  // important flags for determining when the chi2 is used
  virtual bool isSoft()=0; // always include the chi2
  virtual bool isHard()=0; // inlcude the chi2 when iterating, but final chi2 of fit does not include it 
                           // - since it is an exact constraint to a hypersurface  
  
 protected:
  std::string name_;
  unsigned int idx_;
  TMatrixT<double> *val_;
  TMatrixTSym<double> *cov_;
  TMatrixT<double> *dalpha_;
  double initalSigma_;
  double sigma_min_;
  double dx_;
  unsigned int nIterations_;
  unsigned int maxIterations_;
  std::vector<unsigned int> corridx_;
  double scalefactor_;

};
#endif

