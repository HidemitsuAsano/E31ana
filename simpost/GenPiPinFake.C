#include "../src/GlobalVariables.h"

const double MinMomentumPi = 0.06;// GeV/c
const double MinMomentumN = 0.14;// GeV/c
const double cdhL = 79.0;//cm
const double cdhLV = -14.3; //cm/ns
const double cdhR = 56.0; //cm
const double costheta_max = cdhL/2.0/sqrt(cdhLV*cdhLV/4.0+cdhR*cdhR);
const double costheta_min = -costheta_max;

bool checkAcceptance(TLorentzVector* pip,TLorentzVector *pim,TLorentzVector *n);

void GenPiPinFake(){

	if (!gROOT->GetClass("TGenPhaseSpace")) gSystem->Load("libPhysics");

  static const double mom_beam   = 1.05;
  TLorentzVector* target = new TLorentzVector(0.0, 0.0, 0.0, dMass);
  TLorentzVector* beam = new TLorentzVector(0.0, 0.0, mom_beam, sqrt(mom_beam*mom_beam+kpMass*kpMass));
	TLorentzVector* W = new TLorentzVector(*target + *beam);
  
  const unsigned int EventNum=1e8; 
	TGenPhaseSpace event;
  TRandom3 *rand3 = new TRandom3();
  rand3->SetSeed(1);
  for(unsigned int ievt=0;ievt<EventNum;ievt++){
    if(ievt%1e7==0) std::cout << ievt << std::endl;
    double MissMass = rand3->Uniform(0,1.5);//MAX 1.5GeV
    //std::cout << "MissMass " << MissMass << std::endl;

    Double_t masses[4] = {piMass , piMass, nMass, MissMass};
    event.SetDecay(*W,4, masses);
    Double_t weight = event.Generate();
    TLorentzVector *LVec_pip = event.GetDecay(0);//lab
    TLorentzVector *LVec_pim = event.GetDecay(1);//lab
    TLorentzVector *LVec_n   = event.GetDecay(2);//lab
    TLorentzVector *LVec_miss = event.GetDecay(3);//lab
    bool flag_acc = checkAcceptance(LVec_pip,LVec_pim,LVec_n);
  }






















};

bool checkAcceptance(TLorentzVector* pip,TLorentzVector *pim,TLorentzVector *n){
	
  double momPiP = pip->P();
	double costhetaPiP = pip->CosTheta();	
	if(momPiP<=MinMomentumPi) return false;
	if(costhetaPiP<=costheta_min || costheta_max<=costhetaPiP) return false;

	double momPiM = pim->P();
	double costhetaPiM = pim->CosTheta();	
	if(momPiM<=MinMomentumPi) return false;
	if(costhetaPiM<=costheta_min || costheta_max<=costhetaPiM) return false;


  double momN = n->P();
  double costhetaN = n->CosTheta();
  if(momN<MinMomentumN) return false;
  if(costhetaN<=costheta_min || costheta_max<=costhetaN) return false;


  return true;
}
