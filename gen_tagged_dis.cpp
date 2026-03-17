#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <unistd.h>

#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TLorentzVector.h"

#include "constants.h"

using namespace std;

int main(int argc, char ** argv)
{
  cout << "Hello world.\n";
  if (argc != 3)
    {
      cerr << "Wrong number of arguments.\n"
	   << "Instead use:\n\n"
	   << "   gen_tagged_dis /path/to/output/file [Nevents]\n\n"
	   << "Aborting...\n";
      return -1;
    }

  const int Nevents = atoi(argv[2]);
  TFile * outfile = new TFile(argv[1],"RECREATE");
  TTree * outtree = new TTree("genT","Generated Events Tree");

  // Let's make some branches
  double mom_e, theta_e, phi_e, mom_s, theta_s, phi_s, weight;
  int specID; // PDG code for spectator nucleon 
  outtree->Branch("mom_e",&mom_e,"mom_e/D");
  outtree->Branch("theta_e",&theta_e,"theta_e/D");
  outtree->Branch("phi_e",&phi_e,"phi_e/D");
  outtree->Branch("mom_s",&mom_s,"mom_s/D");
  outtree->Branch("theta_s",&theta_s,"theta_s/D");
  outtree->Branch("phi_s",&phi_s,"phi_s/D");
  outtree->Branch("weight",&weight,"weight/D");
  outtree->Branch("specID",&specID,"specID/I");

  // Create other helper functions
  TRandom3 myRand(0);

  // Figure out the boost kinematics
  const double pBeam_e=5.;
  const double pBeam_d=41.;
  const double EBeam_d=TMath::Sqrt(pBeam_d*pBeam_d + mH2*mH2);
  const double crab_angle=2.*M_PI/180.;  // 2 degrees converted to radians
  TLorentzVector v4beam_e(pBeam_e*TMath::Sin(crab_angle), 0., -pBeam_e*TMath::Cos(crab_angle),pBeam_e);
  TLorentzVector v4beam_d(0.,0.,pBeam_d,EBeam_d);

  TLorentzVector v4beam_e_ft=v4beam_e;
  v4beam_e_ft.Boost(-v4beam_d.BoostVector());
  cerr << "Kinematics:\n"
       << "    Deuteron: " << pBeam_d  << " GeV/c\n"
       << "    Electron: " << pBeam_e  << " GeV/c\n"
       << "This kinematic setting is equivalent to a fixed-target beam energy of " << v4beam_e_ft.P() << " GeV.\n";
   
  // Loop over the events
  for (int event=0; event <= Nevents ; event++)
    {
      weight=0;
      mom_e=0;
      theta_e=0;
      phi_e=0;
      mom_s=0;
      theta_s=0;
      phi_s=0;
      specID=pCode;

      // Generate random kinematics

      // Evaluate the cross section

      // Boost to the Collider frame
      
      outtree->Fill();
    }

  outtree->Write();
  outfile->Close();
  
  return 0;
}
