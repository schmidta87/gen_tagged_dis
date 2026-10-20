#include "Cross_Sections.h"
#include <cmath>
#include "constants.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <time.h>


int plotV1(double pPerp, double a1, double m1, double m2){
    std::ofstream dataFile("v1Data.txt");
    if(!dataFile.is_open()) {
        std::cerr << "Could not open data file" << std::endl;
        return 1;
    }
    for(double i = 0; i <= 2; i+=0.01){
        double v = v1(pPerp, i, mP, mN);
        dataFile << i*1.0 << " " << v << '\n';   
    }
    dataFile.close();
    std::string filename = "v1Data.txt";

    std::string cmd =
    "gnuplot -persist -e \""
    "set xlabel 'a1'; "
    "set ylabel 'virtuality'; "
  //  "set logscale y;"
    "plot '" + filename + "' using 1:2 with lines title 'v1'\"";    std::system(cmd.c_str());
    
    return 0;
}

int plotCSsw(double Qsq, double xB, double phi, double a2, double pPerp2, double E1){
    std::ofstream dataFile("csData.txt");
    if(!dataFile.is_open()) {
        std::cerr << "Could not open data file" << std::endl;
        return 1;
    }
    for(double i = 0.1; i < 120; i++){
        double cs = csSW(i, xB, phi, a2, pPerp2, E1);
        dataFile << i*1.0 << " " << cs << '\n';   
    }
    dataFile.close();
    std::string filename = "csData.txt";

    std::string cmd =
    "gnuplot -persist -e \""
    "set xlabel 'Qsq'; "
    "set ylabel 'CSsw(mb)'; "
    "set logscale y;"
    "plot '" + filename + "' using 1:2 with lines title 'v1'\"";    std::system(cmd.c_str());
    
    return 0;
}
int plotCSgriff(double QSq, double xB, double phi, double a2, double pPerp2, double E1, double theta){
    std::ofstream dataFile("csDataGriff.txt");
    if(!dataFile.is_open()) {
        std::cerr << "Could not open data file" << std::endl;
        return 1;
    }
    for(double i = 0; i < 3.14; i+=0.1){
        double cs = csStruckNeutron(QSq, xB, phi, a2, pPerp2, E1);
        dataFile << i*1.0 << " " << cs << '\n';   
    }
    dataFile.close();
    std::string filename = "csDataGriff.txt";

    std::string cmd =
    "gnuplot -persist -e \""
    "set xlabel 'Qsq'; "
    "set ylabel 'CSgriff (mb)'; "
    "set logscale y;"
    "plot '" + filename + "' using 1:2 with lines title 'v1'\"";    
    std::system(cmd.c_str());
    
    return 0;
}

void roughMC(int N){
    /* generate some Qsq, xB, phi, and pPerp from some probability function P
        get the weights as the dCS/P, save this with the xB and phi kinematics to get the total cross section
    */
   //For a first pass I'm just going to use a normal distribution

   //There is a limit on the (QSq, xB) when using the rest frame momentum for all of this since QSq = 4EE'(1-cos(theta)) must work for some E' = E1-x, E' can be much larger in the collider frame than the rest frame
//   Q = 2*xB*mH2*(E1*E1) * (1 - cos( theta))/(E*(1-cos(theta) - xB*mH2 ))
//we can iterate through xB bins and generate random QSq from the max and min defined above
//How to define convergence of the MC? 
//FOr what 
   double xBmin = 0.1;
   double xBmax = 0.7;
   double QSqMin = 1;
   double QSqMax = 100;
   double pPerpMin = 0.001;
   double pPerpMax = 0.5;
   double xB, QSq, pPerp, phi, wCS;
   //Add in some a2 range
   //We can sample pPper from AV18 
   double a2 = 1.0; //Just arbitrary for now
   std::ofstream dataFile("mcData.dat");
   srand(time(NULL));
   for(int i =0; i < N; i++){
        phi = (rand() % 101)*0.01*2*3.1415;
        xB = (rand() % 101)*0.01*(xBmax - xBmin);
        QSq = (rand() % 101)*0.01*(QSqMax-QSqMin);
        pPerp = (rand() % 101)*0.01*2*(pPerpMax - pPerpMin);
     //   wCS = csStruckNeutron(QSq, xB, phi, a2, pPerp, 10 )* (QSqMax - QSqMin)*(xBmax - xBmin) * (pPerpMax - pPerpMin) * 2*3.1415;
        wCS = csSW(QSq, xB, phi, a2, pPerp, 10 )* (QSqMax - QSqMin)*(xBmax - xBmin) * (pPerpMax - pPerpMin) * 2*3.1415;
     
     dataFile << QSq << ' ' << xB << ' ' << wCS << '\n';  
   }
   dataFile.close();
}

int main()
{
    double QSq = 50.0;
    double xB  = 0.4;
    double Es  = 0.5;
    double a1 = 1.1;
    double a2 = 2-a1;
    double E1 = 10;
    double pPerp = 0.5;
    double phi = 0.3;
    double theta = 0.1;
    
    //plotV1(pPerp, a1, mP, mN);
    //plotCSsw(QSq, xB, phi, a2, pPerp, E1);
    plotCSsw(QSq, xB,phi, a2, pPerp, E1);
    roughMC(10000);
    return 0;
}
