#ifndef __CROSS_SECTIONS_H__
#define __CROSS_SECTIONS_H__

double F1(double xB,double QSq);
double F2(double xB,double QSq);
double R(double xB,double QSq);

double disCS_ft(double E1, double E3, double theta, double thisF1, double thisF2);

double disCS(double s, double xB, double QSq, double thisF1, double thisF2);


#endif

