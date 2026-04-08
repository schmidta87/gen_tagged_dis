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
