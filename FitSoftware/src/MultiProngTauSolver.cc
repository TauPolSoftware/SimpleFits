#include "SimpleFits/FitSoftware/interface/MultiProngTauSolver.h"
#include <iostream>

void MultiProngTauSolver::quadratic(double &x_plus,double &x_minus,double a, double b, double c){
  double R=b*b-4*a*c;
  if(R<0){R=0;}
  x_minus=(-b-sqrt(R))/(2.0*a);
  x_plus=(-b+sqrt(R))/(2.0*a);
}

void MultiProngTauSolver::AnalyticESolver(TLorentzVector &nu1,TLorentzVector &nu2,TLorentzVector A1){
  double a=(A1.Pz()*A1.Pz())/(A1.E()*A1.E())-1.0;
  double K=(mtau*mtau-A1.M2()-2.0*A1.Pt()*A1.Pt())/(2.0*A1.E());
  double b=2.0*K*A1.Pz()/A1.E();
  double c=K*K-A1.Pt()*A1.Pt();
  double z1(0),z2(0);
  quadratic(z1,z2,a,b,c);
  nu1.SetPxPyPzE(-A1.Px(),-A1.Py(),z1,sqrt(z1*z1+A1.Pt()*A1.Pt()));
  nu2.SetPxPyPzE(-A1.Px(),-A1.Py(),z2,sqrt(z2*z2+A1.Pt()*A1.Pt()));
}

void MultiProngTauSolver::NumericalESolver(TLorentzVector &nu1,TLorentzVector &nu2,TLorentzVector A1){
  double rmin(-100), rmax(100), step(0.01), mtau2(pow(mtau,2.0)), z1(-999), z2(-999), zmin(-999), min(9999), prev(9999);
  double z=rmin;
  TLorentzVector nu,tau;
  for(int i=0;i<=(int)(rmax-rmin)/step;i++){
    nu.SetPxPyPzE(-A1.Px(),-A1.Py(),z,sqrt(z*z+A1.Pt()*A1.Pt()));
    tau=A1+nu;
    double m2=tau.M2();
    if(m2-mtau2<0 && prev-mtau2>=0) z1=z;
    if(m2-mtau2>0 && prev-mtau2<=0) z2=z;
    if(min>m2){ zmin=z; min=m2;}
    prev=m2;
    z+=step;
  }
  if(z1!=-999 && z2!=-999){
    nu1.SetPxPyPzE(-A1.Px(),-A1.Py(),z1,sqrt(z1*z1+A1.Pt()*A1.Pt()));
    nu2.SetPxPyPzE(-A1.Px(),-A1.Py(),z2,sqrt(z2*z2+A1.Pt()*A1.Pt()));
  }
  else{
    nu1.SetPxPyPzE(-A1.Px(),-A1.Py(),zmin,sqrt(zmin*zmin+A1.Pt()*A1.Pt()));
    nu2.SetPxPyPzE(-A1.Px(),-A1.Py(),zmin,sqrt(zmin*zmin+A1.Pt()*A1.Pt()));
  }
}

void MultiProngTauSolver::SolvebyRotation(TVector3 TauDir,TLorentzVector A1, TLorentzVector &Tau1,TLorentzVector &Tau2,
					  TLorentzVector &nu1,TLorentzVector &nu2, bool rotateback){
  TLorentzVector A1rot=A1;
  double phi(TauDir.Phi()),theta(TauDir.Theta());
  A1rot.RotateZ(-phi);
  A1rot.RotateY(-theta);
  /////////////////////////////////////////////////////
  //  NumericalESolver(nu1,nu2,A1rot); // for debugging AnalyticESolver (slow)
  /* std::cout << "NumericalESolver nu" << std::endl;
  nu1.Print();
  nu2.Print();
  AnalyticESolver(nu1,nu2,A1rot);
  std::cout << "AnalyticESolver nu" << std::endl;
  nu1.Print();
  nu2.Print();*/
  AnalyticESolver(nu1,nu2,A1rot);
  /////////////////////////////////////////////////////
  if(rotateback){
    nu1.RotateY(theta);
    nu1.RotateZ(phi);
    Tau1=A1+nu1;
    //
    nu2.RotateY(theta);
    nu2.RotateZ(phi);
    Tau2=A1+nu2;
  }
  else{
    Tau1=A1rot+nu1;
    Tau2=A1rot+nu2;
  }
  /*
  std::cout << "Solution A1" << std::endl;
  A1.Print();
  std::cout << "Solution nu" << std::endl;
  nu1.Print();
  nu2.Print();
  std::cout << "Solution tau" << std::endl;
  Tau1.Print();
  Tau2.Print();
  TauDir.Print();
  std::cout << " " << theta << " " << phi << std::endl;
  */
}

bool MultiProngTauSolver::SetTauDirectionatThetaGJMax(TLorentzVector a1, double &theta,double &phi){
  double thetaGJMax =asin(( mtau*mtau-a1.M2())/(2.0*mtau*fabs(a1.P())));
  double dtheta=(theta-a1.Theta());
  double dphi=(phi-a1.Phi());
  double dphitheta=sqrt(dtheta*dtheta+dphi*dphi);
  std::cout << " theta " << theta << " phi " << phi << "a1 theta " << a1.Theta() << "at phi " << a1.Phi()<< " thetaGJMax "  << thetaGJMax << " " << dphitheta << std::endl;
  if(thetaGJMax<dphitheta){
    theta=a1.Theta()+dtheta*thetaGJMax/dphitheta*0.98;
    phi=a1.Phi()+dphi*thetaGJMax/dphitheta*0.98;
    std::cout << "modified theta " << theta << " phi " << phi << " dr " << sqrt(pow(theta-a1.Theta(),2.0)+pow(phi-a1.Phi(),2.0)) <<  std::endl;
    return true;
  }
  return false;
}
