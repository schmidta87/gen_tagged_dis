#include "Cross_Sections.h"
#include <cmath>
#include "constants.h"

double F2(double xB, double QSq)
{
  const double C1p = 1.417;
  const double C2p =-0.108;
  const double C3p = 1.486;
  const double C4p =-5.979;
  const double C5p = 3.524;
  const double C6p =-0.011;
  const double C7p =-0.619;
  const double C8p = 1.385;
  const double C9p = 0.270;
  const double C10p=-2.179;
  const double C11p= 4.722;
  const double C12p=-4.363;

  const double onemx=1.-xB;
  const double F2thr = onemx*onemx*onemx*(onemx*(onemx*(onemx*(onemx*C5p + C4p) + C3p) + C2p) + C1p);
  const double A=1.22*exp(3.2*xB);
  const double lam1=C9p+C10p*xB+C11p*xB*xB + C12p*xB*xB*xB;
  const double lam2=(QSq>A)? 0. : C6p + C7p*xB * C8p*xB*xB;
  const double betap=1.;

  return betap*F2thr*(1. + lam1*log(QSq/A) + lam2*pow(log(QSq/A),2.));
}

double R(double xB, double QSq)
{
  double Theta=1. + 12.*QSq/(1.+QSq)*(0.015625 / (xB*xB+0.015625));
  return 0.0635/log(QSq/0.04) *Theta + 0.5747/QSq - 0.3534/(QSq*QSq+0.09);
}

double F1(double xB, double QSq)
{
  return F2(xB,QSq)/(2.*xB*(1.+R(xB,QSq))) * (1. + 4.*mP*mP*xB*xB/QSq);
}

double disCS_ft(double E1, double E3, double theta, double thisF1, double thisF2)
{
  const double nu = E1-E3;
  const double QSq=2.*E1*E3*(1.-cos(theta));
  const double xB=QSq/(2.*mP*nu);
  const double cosTo2 = cos(0.5*theta);
  const double sinTo2 = sin(0.5*theta);
  const double tanTo2 = sinTo2/cosTo2;
  
  return 4.*alpha*alpha*E3*E3*cosTo2*cosTo2/(QSq*QSq) * (thisF2/nu + 2./mP * thisF1*tanTo2*tanTo2);
}

double disCS(double s, double xB, double QSq, double thisF1, double thisF2)
{
  const double y = QSq/(xB*(s-mP*mP));

  return 4.*M_PI*alpha*alpha/(QSq*QSq) * ( (1.-y-mP*mP*QSq/(s-mP*mP)) * thisF2/xB + y*y*thisF1);
}

double alpha1(double a2){ //light cone fraction for the struck nucleon
 return 2 - a2;
}

double alpha2(double pz2, double E2 ){ //light cone fraction for the spectator
  return (pz2 + E2)/mH2;
}

double v1(double pPerp, double a1, double m1, double m2){
  return (0.5*mH2*a1*(mH2 - (pPerp*pPerp + 0.5*m2*m2)/(mH2*(2-a1))) - pPerp - m1*m1)/(m1*m1);
}


double csStruckNeutron(double QSq, double xB, double phi, double a2, double pPerp2, double E1){
  //Just pull the form of this function from Sara's paper? 
  // If so, what should I use as the density functions for the nucleons? 
  //Assuming that we use hbar = 1
  double W1 = F2(xB, QSq)*xB*mN*2/QSq;
  double W2 = F1(xB, QSq)/mN;
//  double theta = acos(QSq/(2*E1*(E1 - QSq/(2*mH2*xB)))-1);
  double theta = acos((2*E1*QSq/2*mH2*xB - 2*E1*E1 -QSq - (QSq/2*mH2*xB)*(QSq/2*mH2*xB)) /(2*E1*(E1 - QSq/(2*mH2*xB)))-1);
  double dCS = pow((alpha/(2*E1*sin(theta/2)*sin(theta/2))), 2) * (2*W1 *sin(theta/2)*sin(theta/2) + W2 * cos(theta/2)*cos(theta/2))*0.389;
  return dCS;

}

double csSW(double QSq, double xB, double phi, double a2, double pPerp2, double E1){
 //Looks like 1/GeV^2 is about 400 barns
  //What bothers me here is that Strikman and Weiss define this as a cross section for the interaction with the deuteron, not specifically one of the nucleons. And the F2 they use is for the deuteron, not just a nucleon
  double y = QSq/(xB*E1*mH2) ;//the 10 here is E of the electron beam
  double epsilon = (1- y- (xB*xB * y*y * mH2*mH2)/QSq)/(1- y +( y*y/2) -(xB*xB * y*y * mH2*mH2)/QSq);
  double dCS = (alpha*y*y)/(QSq*QSq*(1-epsilon))*(F2(xB, QSq))*(1/(2*pow(2*pi, 2)* a2 * mH2 )) *0.389;//Scaling to mb
  return dCS;

}
