#include "../src/GlobalVariables.h"
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TRandom3.h>
#include <TVector3.h>
#include <TGenPhaseSpace.h>
#include <TLorentzVector.h>

TLorentzVector mom_beam;   // 4-momentum(beam)
TLorentzVector mom_beam_Sp;   // 4-momentum(beam)
TLorentzVector mom_beam_Sm;   // 4-momentum(beam)
TLorentzVector mom_target; // 4-momentum(target)
TLorentzVector mom_pip;    // 4-momentum(pi+)
TLorentzVector mom_pim;    // 4-momentum(pi-)
TLorentzVector mom_n;      // 4-momentum(neutron)
TLorentzVector mom_n_beam;      // 4-momentum(neutron)
TLorentzVector mom_n_Sp;      // 4-momentum(neutron)
TLorentzVector mom_n_Sm;      // 4-momentum(neutron)
double NeutralBetaCDH; // veracity of neutral particle on CDH
double NeutralBetaCDH_beam; // veracity of neutral particle on CDH
double NeutralBetaCDH_vtx[2]; // veracity of neutral particle on CDH
double dE;   // energy deposit on CDH
int neutralseg;   // energy deposit on CDH
TVector3 vtx_reaction; //  vertex(reaction)   
TVector3 vtx_pip_beam; //  
TVector3 vtx_pim_beam; //   
TVector3 vtx_pip_cdc;
TVector3 vtx_pim_cdc;
TVector3 CA_pip;
TVector3 CA_pim;
TVector3 CDH_Pos;
TVector3 CDH_Pos_pip;
TVector3 CDH_Pos_pim;
int run_num;   // run number
int event_num; // event number
int block_num; // block number
int mc_nparticle;

//GEANT4 info after decay, interaction with detector, momentum is true
TLorentzVector mcmom_beam;  // generated 4-momentum(beam)
TLorentzVector mcmom_pip;   // generated 4-momentum(pi+)
TLorentzVector mcmom_pim;   // generated 4-momentum(pi-)
TLorentzVector mcmom_ncds;    // generated 4-momentum(neutron)
TLorentzVector mcmom_nmiss;    // generated 4-momentum(neutron)


//GEANT4 info generated particle before interaction
TLorentzVector react_nmiss;
TLorentzVector react_Sigma;//Sigma+ or Simga -
TLorentzVector react_pi;//pi- or pi +


TVector3 mc_vtx;
TLorentzVector kfMomBeamSpmode;   // 4-momentum(beam) after kinematical refit for pi- Sigma+
TLorentzVector kfMom_pip_Spmode;    // 4-momentum(pi+) after kinematical refit for pi- Sigma+
TLorentzVector kfMom_pim_Spmode;    // 4-momentum(pi-) after kinematical refit for pi- Sigma+
TLorentzVector kfMom_n_Spmode;      // 4-momentum(neutron) after kinematical refit for pi- Sigma+
double kf_chi2_Spmode;   // chi2 of kinematical refit
double kf_NDF_Spmode;    // NDF of kinematical refit
double kf_status_Spmode; // status of kinematical refit -> details can be found in this code
double kf_pvalue_Spmode; // p-value of kinematical refit
TLorentzVector kfMomBeamSmmode;   // 4-momentum(beam) after kinematical refit for pi+ Sigma-
TLorentzVector kfMom_pip_Smmode;    // 4-momentum(pi+) after kinematical refit for pi+ Sigma-
TLorentzVector kfMom_pim_Smmode;    // 4-momentum(pi-) after kinematical refit for pi+ Sigma-
TLorentzVector kfMom_n_Smmode;      // 4-momentum(neutron) after kinematical refit for pi+ Sigma-
double kf_chi2_Smmode;   // chi2 of kinematical refit
double kf_NDF_Smmode;    // NDF of kinematical refit
double kf_status_Smmode; // status of kinematical refit -> details can be found in this code
double kf_pvalue_Smmode; // p-value of kinematical refit
int kf_flag; // flag of correct pair reconstruction, etc
const double MinMomentumPi = 0.06;// GeV/c
const double MinMomentumN = 0.14;// GeV/c
const double MaxMomentumN = 1.0;
const double cdhL = 79.0;//cm
const double cdhLV = -14.3; //cm/ns
const double cdhR = 56.0; //cm
const double costheta_max = cdhL/2.0/sqrt(cdhLV*cdhLV/4.0+cdhR*cdhR);
const double costheta_min = -costheta_max;

bool checkAcceptance(TLorentzVector* pip,TLorentzVector *pim,TLorentzVector *n);

void GenPiPinFakeK0(int seed=1){

	if (!gROOT->GetClass("TGenPhaseSpace")) gSystem->Load("libPhysics");

  // std::string treefile_name("fakepippimn_pippimn.root");
  TFile *treefile = new TFile( Form("/gpfs/group/had/knucl/e15/asano/sim/K0nn_pippimn%d.root",seed), "recreate");
  std::cout << treefile->GetName() << std::endl;
  treefile->cd();
  TTree *npippimTree = new TTree( "EventTree", "EventTree");
  npippimTree->Branch( "mom_beam",   &mom_beam );//
  npippimTree->Branch( "mom_beam_Sp",  &mom_beam_Sp );//
  npippimTree->Branch( "mom_beam_Sm",  &mom_beam_Sm );
  npippimTree->Branch( "mom_target", &mom_target );
  npippimTree->Branch( "mom_pip", &mom_pip );
  npippimTree->Branch( "mom_pim", &mom_pim );
  npippimTree->Branch( "mom_n", &mom_n );
  npippimTree->Branch( "mom_n_beam", &mom_n_beam );
  npippimTree->Branch( "mom_n_Sp", &mom_n_Sp );
  npippimTree->Branch( "mom_n_Sm", &mom_n_Sm );
  npippimTree->Branch( "NeutralBetaCDH", &NeutralBetaCDH );
  npippimTree->Branch( "NeutralBetaCDH_beam", &NeutralBetaCDH_beam );
  npippimTree->Branch( "NeutralBetaCDH_vtx[2]", NeutralBetaCDH_vtx);
  npippimTree->Branch( "dE", &dE );
  npippimTree->Branch( "neutralseg", &neutralseg );
  npippimTree->Branch( "vtx_reaction", &vtx_reaction );
  npippimTree->Branch( "vtx_pip_beam", &vtx_pip_beam );
  npippimTree->Branch( "vtx_pim_beam", &vtx_pim_beam );
  npippimTree->Branch( "vtx_pip_cdc", &vtx_pip_cdc );
  npippimTree->Branch( "vtx_pim_cdc", &vtx_pim_cdc );
  npippimTree->Branch( "CA_pip",&CA_pip);
  npippimTree->Branch( "CA_pim",&CA_pim);
  npippimTree->Branch( "CDH_Pos",&CDH_Pos);
  npippimTree->Branch( "CDH_Pos_pim",&CDH_Pos_pim);
  npippimTree->Branch( "CDH_Pos_pip",&CDH_Pos_pip);
  //npippimTree->Branch( "run_num", &run_num );
  //npippimTree->Branch( "event_num", &event_num );
  //npippimTree->Branch( "block_num", &block_num );
  //npippimTree->Branch( "mc_nparticle",   &mc_nparticle );
  npippimTree->Branch( "mcmom_beam",   &mcmom_beam );
  npippimTree->Branch( "mcmom_pip", &mcmom_pip );
  npippimTree->Branch( "mcmom_pim", &mcmom_pim );
  npippimTree->Branch( "mcmom_ncds", &mcmom_ncds );
  npippimTree->Branch( "mcmom_nmiss", &mcmom_nmiss );
  npippimTree->Branch( "react_nmiss",&react_nmiss);
  npippimTree->Branch( "react_Sigma",&react_Sigma);
  npippimTree->Branch( "react_pi",&react_pi);
  npippimTree->Branch( "mc_vtx", &mc_vtx );
  npippimTree->Branch( "kfSpmode_mom_beam",   &kfMomBeamSpmode );
  npippimTree->Branch( "kfSpmode_mom_pip", &kfMom_pip_Spmode );
  npippimTree->Branch( "kfSpmode_mom_pim", &kfMom_pim_Spmode );
  npippimTree->Branch( "kfSpmode_mom_n", &kfMom_n_Spmode );
  npippimTree->Branch( "kfSpmode_chi2", &kf_chi2_Spmode );
  npippimTree->Branch( "kfSpmode_NDF", &kf_NDF_Spmode );
  npippimTree->Branch( "kfSpmode_status", &kf_status_Spmode );
  npippimTree->Branch( "kfSpmode_pvalue", &kf_pvalue_Spmode );
  npippimTree->Branch( "kfSmmode_mom_beam",   &kfMomBeamSmmode );
  npippimTree->Branch( "kfSmmode_mom_pip", &kfMom_pip_Smmode );
  npippimTree->Branch( "kfSmmode_mom_pim", &kfMom_pim_Smmode );
  npippimTree->Branch( "kfSmmode_mom_n", &kfMom_n_Smmode );
  npippimTree->Branch( "kfSmmode_chi2", &kf_chi2_Smmode );
  npippimTree->Branch( "kfSmmode_NDF", &kf_NDF_Smmode );
  npippimTree->Branch( "kfSmmode_status", &kf_status_Smmode );
  npippimTree->Branch( "kfSmmode_pvalue", &kf_pvalue_Smmode );
  npippimTree->Branch( "kf_flag", &kf_flag );


  const unsigned int EventNum=5e6; 
  TGenPhaseSpace event;
  TGenPhaseSpace event_2decay;
  TRandom3 *rand3 = new TRandom3(seed);
  //rand3->SetSeed(seed);
  //checking histograms.
  TH2D *h2pippimmom_gen = new TH2D("h2pippimmom_gen","pippimmom.",1000,0,2,1000,0,2);
  h2pippimmom_gen->SetXTitle("pip mom");
  h2pippimmom_gen->SetYTitle("pim mom");
  TH2D *h2pippimmom = new TH2D("h2pippimmom","pippimmom.",1000,0,2,1000,0,2);
  h2pippimmom->SetXTitle("pip mom");
  h2pippimmom->SetYTitle("pim mom");
  TH2D *h2nmissmom = new TH2D("h2nmissmom","nmissmom.",1000,0,2,1000,0,2);
  h2nmissmom->SetXTitle("nmom");
  h2nmissmom->SetYTitle("miss. mom.");
  
  TH2D* h2IMnpimIMnpip = new TH2D("h2IMnpimIMnpip","IMnpim vs IMnpip",1000,1,2,1000,1,2);
  TH2D* h2qIMnpipi = new TH2D("h2qIMnpipi","h2qIMnpipi",1000,1,2,1000,0,2);
  TH2D* h2qIMnpipi_w = new TH2D("h2qIMnpipi_w","h2qIMnpipi_w",1000,1,2,1000,0,2);
  TH1D* h1missmass = new TH1D("h1missmass","miss. mass",1500,0,1.5);
  TH1D* h1IMpippim = new TH1D("h1IMpippim","IMpippim",500,0,1);
  TH1D* h1Mompippim = new TH1D("h1Mompippim","Mompippim",500,0,1);
  TH1D* h1coslabpip = new TH1D("h1coslabpip","coslabpip",100,-1,1);
  TH1D* h1coslabpim = new TH1D("h1coslabpim","h1coslabpim",100,-1,1);

  static const double mombeam   = 1.00;
  TLorentzVector* target = new TLorentzVector(0.0, 0.0, 0.0, dMass);
  TLorentzVector* beam = new TLorentzVector(0.0, 0.0, mombeam, sqrt(mombeam*mombeam+kpMass*kpMass));
	TLorentzVector* W = new TLorentzVector(*target + *beam);
  
  unsigned int ievt=0;
  while(ievt<EventNum){
    if(ievt%500000==0) std::cout << ievt << std::endl;
    
    double k0Massres = k0Mass+gRandom->Gaus(0,0.0074845);//K0 mass reso.
    Double_t masses[3] = {k0Massres, nMass, nMass};
    event.SetDecay(*W,3, masses);
    Double_t weight = event.Generate();
    TLorentzVector *LVec_K0 = event.GetDecay(0);//lab
    TLorentzVector *LVec_n   = event.GetDecay(1);//lab
    TLorentzVector *LVec_miss = event.GetDecay(2);//lab
    Double_t masses_2decay[2] = {piMass,piMass};
    event_2decay.SetDecay(*LVec_K0,2, masses_2decay);
    Double_t weight_2body = event_2decay.Generate();
    TLorentzVector *LVec_pip = event_2decay.GetDecay(0);//lab
    TLorentzVector *LVec_pim = event_2decay.GetDecay(1);//lab
    h2pippimmom_gen->Fill((*LVec_pip).P(),(*LVec_pim).P());
    bool flag_acc = checkAcceptance(LVec_pip,LVec_pim,LVec_n);
    if(!flag_acc) continue;
    h1coslabpip->Fill((*LVec_pip).CosTheta());
    h1coslabpim->Fill((*LVec_pim).CosTheta());
    h1IMpippim->Fill((*LVec_pip+*LVec_pim).M());
    h1Mompippim->Fill((*LVec_pip+*LVec_pim).P());
    //fake beam mom. calculated from missing mom.
    TLorentzVector LVec_beam = *LVec_pip + *LVec_pim + *LVec_n - *LVec_miss - *target ;

    //fillTree
    mom_beam = LVec_beam;
    mom_target = *target;
    mom_pip = *LVec_pip;
    mom_pim = *LVec_pim;
    mom_n = *LVec_n;
    dE = 5.0;
    NeutralBetaCDH = 0.5;
    npippimTree->Fill();
    h2pippimmom->Fill((*LVec_pip).P(),(*LVec_pim).P());
    h2nmissmom->Fill((*LVec_n).P(),(*LVec_miss).P());
    h1missmass->Fill((*LVec_miss).M());
    TLorentzVector LVec_n_pip = *LVec_pip + *LVec_n;
    TLorentzVector LVec_n_pim = *LVec_pim + *LVec_n;
    h2IMnpimIMnpip->Fill(LVec_n_pip.M(),LVec_n_pim.M());
    TLorentzVector LVec_n_pip_pim = *LVec_pip + *LVec_pim + *LVec_n;
    TLorentzVector q = *beam - *LVec_miss;
    h2qIMnpipi->Fill(LVec_n_pip_pim.M(),q.P());
    h2qIMnpipi_w->Fill(LVec_n_pip_pim.M(),q.P(),weight);
    ievt++;
  }
  
  h2pippimmom->Write();
  h2nmissmom->Write();
  h1missmass->Write();
  npippimTree->Write();
  h2IMnpimIMnpip->Write();
  h2qIMnpipi->Write();
  h2qIMnpipi_w->Write();
  treefile->Write();
  treefile->Close();

  return;
};


bool checkAcceptance(TLorentzVector *pip,TLorentzVector *pim,TLorentzVector *n)
{

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
  if(momN>MaxMomentumN) return false;
  if(costhetaN<=costheta_min || costheta_max<=costhetaN) return false;


  return true;
}
