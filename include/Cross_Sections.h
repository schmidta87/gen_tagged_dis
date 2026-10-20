#ifndef __CROSS_SECTIONS_H__
#define __CROSS_SECTIONS_H__

double F1(double xB,double QSq);
double F2(double xB,double QSq);
double R(double xB,double QSq);

double disCS_ft(double E1, double E3, double theta, double thisF1, double thisF2);

double disCS(double s, double xB, double QSq, double thisF1, double thisF2);
double alpha1(double pz, double Es); //light cone fraction for the struck nucleon
double alpha2(double a1 ); //light cone fraction for the struck nucleon
double csStruckNeutron(double QSq, double xB, double phi, double a1, double pPerp, double E1);
//double v1(double Qsq, double xB, double a1, double M1);
double v1( double pPerp, double a1, double m1, double m2);
double csSW(double QSq, double xB, double phi, double a2, double pPerp2, double E1);

#endif

