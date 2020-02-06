#include "../src/GlobalVariables.h"

const double MinMomentumPi = 0.1;// GeV/c
const double MinMomontumN = 0.14;// GeV/c
const double cdhL = 79.0;//cm
const double cdhLV = -14.3; //cm/ns
const double cdhR = 56.0; //cm
const double costheta_max = cdhL/2.0/sqrt(cdhLV*cdhLV/4.0+cdhR*cdhR);
const double constheta_min = -costheta_max;


void GenPiPinFake(){

	if (!gROOT->GetClass("TGenPhaseSpace")) gSystem->Load("libPhysics");

  static const double mom_beam   = 1.05;
  TLorentzVector* target = new TLorentzVector(0.0, 0.0, 0.0, dMass);
  TLorentzVector* beam = new TLorentzVector(0.0, 0.0, mom_beam, sqrt(mom_beam*mom_beam+kpMass*kpMass));
	TLorentzVector* W = new TLorentzVector(*target + *beam);
  
  const unsigned int EventNum=1e9; 
	TGenPhaseSpace event;
  TRandom3 *rand3 = new TRandom3();
  rand3->SetSeed(1);
  for(unsigned int ievt=0;ievt<EventNum;ievt++){
    double MissMass = rand3->Uniform(0,1.5);//MAX 1.5GeV
    //std::cout << "MissMass " << MissMass << std::endl;

    Double_t masses[4] = {piMass , piMass, nMass, MissMass};
    event.SetDecay(*W,4, masses);
    Double_t weight = event.Generate();
    TLorentzVector *LVec_pip = event.GetDecay(0);
    TLorentzVector *LVec_pim = event.GetDecay(1);
    TLorentzVector *LVec_n   = event.GetDecay(2);
    TLorentzVector *LVec_miss = event.GetDecay(3);
    bool flag_acc = checkAcceptance(LVec_pip,LVec_pim,LVec_n);
  }






















};

bool checkAcceptance(TLorentzVector* pip,TLorentzVector *pim,TLorentzVector *n){
  bool flag = true;

	// for pip
	double momentum_pip = pip->P();
	double costheta_pip = pip->CosTheta();	
	if(momentum_pip<=momentum_min) flag = false;
	if(costheta_pip<=costheta_min || costheta_pip<=costheta_k) flag = false;






}
