#include "wiser_pion.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Random/RandFlat.h"
#include <iostream>
#include "TFile.h"
#include "TTree.h"

using namespace std;

int main(int argc, char** argv){

  if( argc == 1 || (0 == strcmp("--help", argv[1]))) {
    cout << " usage: build/runWiser [options]" << endl
         << " --nEvents <sample number>" << endl
         << " --beamE <beam energy in GeV>" << endl
      //liquid H2 RadLen 890.4cm ... 20 cm tgt is 2.25% radiator      
         << " --radLen <radiation Length in procent/100>"<< endl
         << " --momMinMax <min and max accepted momentum in GeV>"<< endl
	 << " --thetaCentral <central theta value>"<<endl
      //theta range
         << " --thetaMinMax <theta min[rad]> <theta max[rad]>" << endl;
    return 1;
  }

  int nEvents;
  double beamE,radLen,thMin,thMax,momMin,momMax,thCentral(999);
  for(int i = 1; i < argc; i++) {
    if(0 == strcmp("--beamE", argv[i])) {
      beamE = atof(argv[i+1]);
    }
    if(0 == strcmp("--nEvents", argv[i])) {
      nEvents = atoi(argv[i+1]);
    }
    if(0 == strcmp("--radLen", argv[i])) {
      radLen = atof(argv[i+1]);
    }
    if(0 == strcmp("--thetaCentral", argv[i])) {
      thCentral = atof(argv[i+1])*CLHEP::deg;
    }
    if(0 == strcmp("--thetaMinMax", argv[i])) {
      thMin = atof(argv[i+1])*CLHEP::deg;
      thMax = atof(argv[i+2])*CLHEP::deg;
    }
    if(0 == strcmp("--momMinMax", argv[i])) {
      momMin = atof(argv[i+1]);
      momMax = atof(argv[i+2]);
    }
  }

  double th,pf,sigpip,sigpim,effXsect;
  TFile *fout=new TFile("o_wiser.root","RECREATE");
  TTree *t=new TTree("t","wiser tree for specific configuration");
  t->Branch("th",&th,"th/D");
  t->Branch("pf",&pf,"pf/D");
  t->Branch("sigpip",&sigpip,"sigpip/D");
  t->Branch("sigpim",&sigpim,"sigpim/D");
  t->Branch("effXsect",&effXsect,"effXsect/D");

  //const double alpha=1/137.036;//defined in wiser.h as 0.007299
  const double Z=1;
  const double A=1;
  const double intrad = 2.0*alpha*log(beamE/CLHEP::electron_mass_c2)/CLHEP::pi;
  const double V = 2*CLHEP::pi * (cos(thMin) - cos(thMax)) * beamE;

  for(int i=0;i<nEvents;i++){
    if(i%10000==1) cout<<"processed "<<i<<endl;

    if(thCentral == 999)
      th = CLHEP::RandFlat::shoot(thMin, thMax);
    else
      th=thCentral;
    
    pf = CLHEP::RandFlat::shoot(momMin  , momMax);
    //pf=beamE;
    
    sigpip = wiser_sigma(beamE, pf, th, radLen*4.0/3.0 , 0);
    sigpim = wiser_sigma(beamE, pf, th, radLen*4.0/3.0 , 1);
    // sigpip = wiser_sigma(beamE, pf, th, radLen*4.0/3.0 + intrad, 0)*CLHEP::nanobarn/CLHEP::GeV;
    // sigpim = wiser_sigma(beamE, pf, th, radLen*4.0/3.0 + intrad, 1)*CLHEP::nanobarn/CLHEP::GeV;

    effXsect = V * ( Z*sigpim + A/Z*sigpip );
    t->Fill();
  }

  fout->cd();
  t->Write();
  fout->Close();
    
  return 0;
}
