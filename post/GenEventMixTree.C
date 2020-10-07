#include <iostream>
#include <vector>
#include <string>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TH2.h>

#include "globalana.h"

TLorentzVector *LVec_beam2=NULL;   // 4-momentum(beam)
TLorentzVector *LVec_beam_Sp2=NULL;   // 4-momentum(beam),Sp mode assumption
TLorentzVector *LVec_beam_Sm2=NULL;   // 4-momentum(beam),Sm mode assumption
TLorentzVector *LVec_target2=NULL; // 4-momentum(target)
TLorentzVector *LVec_pip2=NULL;    // 4-momentum(pi+)
TLorentzVector *LVec_pim2=NULL;    // 4-momentum(pi-)
TLorentzVector *LVec_n2=NULL;      // 4-momentum(neutron)
TLorentzVector *LVec_n_beam2=NULL;      // 4-momentum(neutron)
TLorentzVector *LVec_n_Sp2=NULL;      // 4-momentum(neutron),Sp mode assumption
TLorentzVector *LVec_n_Sm2=NULL;      // 4-momentum(neutron),Sm mode assumption
TLorentzVector *mcmom_beam2=NULL;   // generated 4-momentum(beam)
TLorentzVector *mcmom_pip2=NULL;    // generated 4-momentum(pi+)
TLorentzVector *mcmom_pim2=NULL;    // generated 4-momentum(pi-)
TLorentzVector *mcmom_ncds2=NULL;      // generated 4-momentum(neutron)
TLorentzVector *mcmom_nmiss2=NULL;      // generated 4-momentum(neutron)
TLorentzVector *react_nmiss2=NULL;      // generated 4-momentum(neutron)
TLorentzVector *react_Sigma2=NULL;      // generated 4-momentum(neutron)
TLorentzVector *react_pi2=NULL;      // generated 4-momentum(neutron)
double NeutralBetaCDH2; // velocity of neutral particle on CDH
double NeutralBetaCDH_beam2; // velocity of neutral particle on CDH
double NeutralBetaCDH_vtx2[2]; // velocity of neutral particle on CDH,0: Spmode 1:Smmode
double dE2;   // energy deposit on CDH
int neutralseg2;
double mcncanvtxr2;
double mcncanvtxz2;
int mcncdsgen2;
int mcpattern2;
TVector3 *vtx_reaction2 = NULL; // vertex(reaction) 
TVector3 *vtx_pip_beam2 = NULL; //C.A.P of pip-beam beam side
TVector3 *vtx_pim_beam2 = NULL; //C.A.P of pim-beam beam side
TVector3 *vtx_pip_cdc2 = NULL;//C.A.P of pip-beam pip side
TVector3 *vtx_pim_cdc2 = NULL;//C.A.P of pim-beam pim side
TVector3 *CA_pip2 = NULL;//C.A.P of pip-pim pip side
TVector3 *CA_pim2 = NULL;//C.A.P of pip-pim pim side
TVector3 *CDH_Pos2 = NULL;
TVector3 *CDH_Pos_pim2 = NULL;
TVector3 *CDH_Pos_pip2 = NULL;


void GenEventMixTree(const char* filename = "evanaIMpisigma_npippim_v197.root")
{
  TFile *f = new TFile(filename);
  //TFile *f = new TFile("sim_piSpn_dE0_Al.root");
  TTree *tree = (TTree*)f->Get("EventTree");
  if(tree==0) {
    std::cout << "EventTree is not found " << std::endl;
    std::cout << "end " << std::endl;
    return ;
  }

  bool SimSpmode = (std::string(filename).find("Sp")!= std::string::npos);
  bool SimSmmode = (std::string(filename).find("Sm")!= std::string::npos);
  tree->SetBranchAddress( "mom_beam",   &LVec_beam );
  tree->SetBranchAddress( "mom_beam_Sp",   &LVec_beam_Sp );
  tree->SetBranchAddress( "mom_beam_Sm",   &LVec_beam_Sm );
  tree->SetBranchAddress( "mom_target", &LVec_target );
  tree->SetBranchAddress( "mom_pip", &LVec_pip );
  tree->SetBranchAddress( "mom_pim", &LVec_pim );
  tree->SetBranchAddress( "mom_n", &LVec_n );
  tree->SetBranchAddress( "mom_n_beam", &LVec_n_beam );//from v192
  tree->SetBranchAddress( "mom_n_Sp", &LVec_n_Sp );
  tree->SetBranchAddress( "mom_n_Sm", &LVec_n_Sm );
  tree->SetBranchAddress( "NeutralBetaCDH", &NeutralBetaCDH );
  tree->SetBranchAddress( "NeutralBetaCDH_beam", &NeutralBetaCDH_beam );//from v192
  tree->SetBranchAddress( "NeutralBetaCDH_vtx[2]", NeutralBetaCDH_vtx );
  tree->SetBranchAddress( "dE", &dE );
  tree->SetBranchAddress( "neutralseg", &neutralseg );  
  tree->SetBranchAddress( "vtx_reaction", &vtx_reaction );
  tree->SetBranchAddress( "vtx_pip_beam",&vtx_pip_beam);
  tree->SetBranchAddress( "vtx_pim_beam",&vtx_pim_beam);
  tree->SetBranchAddress( "vtx_pip_cdc",&vtx_pip_cdc);
  tree->SetBranchAddress( "vtx_pim_cdc",&vtx_pim_cdc);
  tree->SetBranchAddress( "CA_pip",&CA_pip);
  tree->SetBranchAddress( "CA_pim",&CA_pim);
  tree->SetBranchAddress( "CDH_Pos",&CDH_Pos);
  tree->SetBranchAddress( "CDH_Pos_pim",&CDH_Pos_pim);//from v193
  tree->SetBranchAddress( "CDH_Pos_pip",&CDH_Pos_pip);//from v193
  //tree->SetBranchAddress( "run_num", &run_num );
  //tree->SetBranchAddress( "event_num", &event_num );
  //tree->SetBranchAddress( "block_num", &block_num );
  if(SimSpmode || SimSmmode) {
    tree->SetBranchAddress( "mcncanvtxr", &mcncanvtxr);
    tree->SetBranchAddress( "mcncanvtxz", &mcncanvtxz);
    tree->SetBranchAddress( "mcncdsgen", &mcncdsgen);
    tree->SetBranchAddress( "mcpattern", &mcpattern);
    tree->SetBranchAddress( "mcmom_beam",  &mcmom_beam );
    tree->SetBranchAddress( "mcmom_pip", &mcmom_pip);
    tree->SetBranchAddress( "mcmom_pim", &mcmom_pim);
    tree->SetBranchAddress( "mcmom_ncds", &mcmom_ncds);
    tree->SetBranchAddress( "mcmom_nmiss", &mcmom_nmiss);
    tree->SetBranchAddress( "react_nmiss", &react_nmiss);
    tree->SetBranchAddress( "react_Sigma", &react_Sigma);
    tree->SetBranchAddress( "react_pi", &react_pi);
  }
   
  Int_t nevent = tree->GetEntries();
  std::cerr<<"# of events = "<<nevent<<std::endl;
  
  for ( Int_t i=0; i<nevent; i++ ) {
    tree->GetEvent(i);
    if(i%50000==0) std::cout << "Event# " << i << std::endl;
    *LVec_beam2 =  *LVec_beam;
    *LVec_beam_Sp2 = *LVec_beam_Sp;
    *LVec_beam_Sm2 = *LVec_beam_Sm;
    *LVec_target2 = *LVec_target;
    *LVec_pip2 = *LVec_pip;
    *LVec_pim2 = *LVec_pim;
    *LVec_n2 = *LVec_n;
    *LVec_n_beam2 = *LVec_n_beam;
    *LVec_n_Sp2 = *LVec_n_Sp;
    *LVec_n_Sm2 = *LVec_n_Sm;
    NeutralBetaCDH2 = NeutralBetaCDH;
    dE2 = dE;
    neutralseg2 = neutralseg;
    *CA_pip2 = *CA_pip;
    *CA_pim2 = *CA_pim;
    *CDH_Pos2 = *CDH_Pos;
    *CDH_Pos_pim2 = *CDH_Pos_pim;
    *CDH_Pos_pip2 = *CDH_Pos_pip;
    if(SimSpmode || SimSmmode){
      *mcmom_beam2 = *mcmom_beam;
      *mcmom_pip2 = *mcmom_pip;
      *mcmom_pim2 = *mcmom_pim;
      *mcmom_ncds2 = *mcmom_ncds;
      *mcmom_nmiss2 = *mcmom_nmiss;
      *react_nmiss2 = *react_nmiss;
      *react_Sigma2 = *react_Sigma;
      *react_pi2 = *react_pi;
    }

  }


  TTree *treeMIX = new TTree("EventTree","EventTreeMIX");
  treeMIX->Branch( "mom_beam",   &LVec_beam2 );//
  treeMIX->Branch( "mom_beam_Sp",  &LVec_beam_Sp2 );//
  treeMIX->Branch( "mom_beam_Sm",  &LVec_beam_Sm2 );
  treeMIX->Branch( "mom_target", &LVec_target2 );
  treeMIX->Branch( "mom_pip", &LVec_pip2 );
  treeMIX->Branch( "mom_pim", &LVec_pim2 );
  treeMIX->Branch( "mom_n", &LVec_n2 );
  treeMIX->Branch( "mom_n_beam", &LVec_n_beam2 );
  treeMIX->Branch( "mom_n_Sp", &LVec_n_Sp2 );
  treeMIX->Branch( "mom_n_Sm", &LVec_n_Sm2 );
  treeMIX->Branch( "NeutralBetaCDH", &NeutralBetaCDH2 );
  treeMIX->Branch( "NeutralBetaCDH_beam", &NeutralBetaCDH_beam2 );
  treeMIX->Branch( "NeutralBetaCDH_vtx[2]", NeutralBetaCDH_vtx2);
  treeMIX->Branch( "dE", &dE2 );
  treeMIX->Branch( "neutralseg", &neutralseg2 );
  treeMIX->Branch( "vtx_reaction", &vtx_reaction2 );
  treeMIX->Branch( "vtx_pip_beam", &vtx_pip_beam2 );
  treeMIX->Branch( "vtx_pim_beam", &vtx_pim_beam2 );
  treeMIX->Branch( "vtx_pip_cdc", &vtx_pip_cdc2 );
  treeMIX->Branch( "vtx_pim_cdc", &vtx_pim_cdc2 );
  treeMIX->Branch( "CA_pip",&CA_pip2);
  treeMIX->Branch( "CA_pim",&CA_pim2);
  treeMIX->Branch( "CDH_Pos",&CDH_Pos2);
  treeMIX->Branch( "CDH_Pos_pim",&CDH_Pos_pim2);
  treeMIX->Branch( "CDH_Pos_pip",&CDH_Pos_pip2);
  if(SimSpmode || SimSmmode) {
    treeMIX->Branch( "mcmom_beam",   &mcmom_beam2 );
    treeMIX->Branch( "mcmom_pip", &mcmom_pip2 );
    treeMIX->Branch( "mcmom_pim", &mcmom_pim2 );
    treeMIX->Branch( "mcmom_ncds", &mcmom_ncds2 );
    treeMIX->Branch( "mcmom_nmiss", &mcmom_nmiss2 );
    treeMIX->Branch( "react_nmiss",&react_nmiss2);
    treeMIX->Branch( "react_Sigma",&react_Sigma2);
    treeMIX->Branch( "react_pi",&react_pi2);
  }

 
  TString outname = std::string(filename);
  outname.Replace(std::string(filename).size()-5,5,"_MIX.root");


  TFile *fout = new TFile(outname.Data(),"RECREATE");
  fout->Print();
  fout->cd();
  treeMIX->Write();
  fout->Close();




}
