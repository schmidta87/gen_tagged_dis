#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <TH2F.h>

void makePlot(){
       std::ifstream dataFile("mcData.dat");    
       double w,Q, x;
        TH2 * csH = new TH2F("csH", "MC results", 100, 0.0, 100, 100, 0.0, 0.7);
        while(dataFile >> Q >> x >> w){
            csH->Fill(Q,x,w);
        }
   
   dataFile.close();
   csH->Draw("COLZ");


}

void plotMC(){

    makePlot();
    
}