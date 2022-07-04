#include <iostream>
#include <vector>
#include <string>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TH2.h>

#include "globalana.h"
#include "anacuts.h"

TLorentzVector *LVec_beam2=NULL;   // 4-momentum(beam)
TLorentzVector *LVec_beam_Sp2=NULL;   // 4-momentum(beam),Sp mode assumption
TLorentzVector *LVec_beam_Sm2=NULL;   // 4-momentum(beam),Sm mode assumption
TLorentzVector *LVec_beam_K02=NULL;   // 4-momentum(beam),K0 mode assumption
TLorentzVector *LVec_target2=NULL; // 4-momentum(target)
TLorentzVector *LVec_pip2=NULL;    // 4-momentum(pi+)
TLorentzVector *LVec_pim2=NULL;    // 4-momentum(pi-)
TLorentzVector *LVec_n2=NULL;      // 4-momentum(neutron)
TLorentzVector *LVec_n_beam2=NULL;      // 4-momentum(neutron)
TLorentzVector *LVec_n_Sp2=NULL;      // 4-momentum(neutron),Sp mode assumption
TLorentzVector *LVec_n_Sm2=NULL;      // 4-momentum(neutron),Sm mode assumption
TLorentzVector *LVec_n_K02=NULL;      // 4-momentum(neutron),Sm mode assumption
TLorentzVector *mcmom_beam2=NULL;   // generated 4-momentum(beam)
TLorentzVector *mcmom_pip2=NULL;    // generated 4-momentum(pi+)
TLorentzVector *mcmom_pim2=NULL;    // generated 4-momentum(pi-)
TLorentzVector *mcmom_ncds2=NULL;      // generated 4-momentum(neutron)
TLorentzVector *mcmom_nmiss2=NULL;      // generated 4-momentum(neutron)
TLorentzVector *react_nmiss2=NULL;      // generated 4-momentum(neutron)
TLorentzVector *react_Sigma2=NULL;      // generated 4-momentum(neutron)
TLorentzVector *react_pi2=NULL;      // generated 4-momentum(neutron)
TLorentzVector *mc_vtx2=NULL;      // generated 4-momentum(neutron)
TLorentzVector *mc_disvtx2=NULL;      // generated 4-momentum(neutron)
double NeutralBetaCDH2; // velocity of neutral particle on CDH
double NeutralBetaCDH_beam2; // velocity of neutral particle on CDH
double NeutralBetaCDH_vtx2[3]; // velocity of neutral particle on CDH,0: Spmode 1:Smmode
double tofpim2;
double tofpip2;
double tofn2;
double dE2;   // energy deposit on CDH
int neutralseg2;
int nhitOutCDC2;
int ForwardCharge2;
double mcncanvtxr2;
double mcncanvtxz2;
int mcncdsgen2;
int mcpattern2;
TVector3 *vtx_reaction2 = NULL; // vertex(reaction) 
TVector3 *vtx_displaced2 = NULL; // vertex(reaction) 
TVector3 *vtx_pip_beam2 = NULL; //C.A.P of pip-beam beam side
TVector3 *vtx_pim_beam2 = NULL; //C.A.P of pim-beam beam side
TVector3 *vtx_pip_cdc2 = NULL;//C.A.P of pip-beam pip side
TVector3 *vtx_pim_cdc2 = NULL;//C.A.P of pim-beam pim side
TVector3 *CA_pip2 = NULL;//C.A.P of pip-pim pip side
TVector3 *CA_pim2 = NULL;//C.A.P of pip-pim pim side
TVector3 *CDH_Pos2 = NULL;
TVector3 *CDH_Pos_pim2 = NULL;
TVector3 *CDH_Pos_pip2 = NULL;


void GenEventMixTree(const char* filename = "evanaIMpisigma_npippim_v245.root")
{
  if(gROOT->GetVersionInt() < 60000){
    std::cout << "Use ROOT6 !!" << std::endl;
    return;
  }
  TFile *f = new TFile(filename,"READ");
  TTree *tree = (TTree*)f->Get("EventTree");
  if(tree==0) {
    std::cout << "EventTree is not found " << std::endl;
    std::cout << "end " << std::endl;
    return ;
  }

  bool SimSpmode = (std::string(filename).find("Sp")!= std::string::npos);
  bool SimSmmode = (std::string(filename).find("Sm")!= std::string::npos);
  bool SimK0nnmode = (std::string(filename).find("K0nn")!= std::string::npos);
  tree->SetBranchAddress( "mom_beam",   &LVec_beam );
  tree->SetBranchAddress( "mom_beam_Sp",   &LVec_beam_Sp );
  tree->SetBranchAddress( "mom_beam_Sm",   &LVec_beam_Sm );
  tree->SetBranchAddress( "mom_beam_K0",   &LVec_beam_K0 );
  tree->SetBranchAddress( "mom_target", &LVec_target );
  tree->SetBranchAddress( "mom_pip", &LVec_pip );
  tree->SetBranchAddress( "mom_pim", &LVec_pim );
  tree->SetBranchAddress( "mom_n", &LVec_n );
  tree->SetBranchAddress( "mom_n_beam", &LVec_n_beam );//from v192
  tree->SetBranchAddress( "mom_n_Sp", &LVec_n_Sp );
  tree->SetBranchAddress( "mom_n_Sm", &LVec_n_Sm );
  tree->SetBranchAddress( "mom_n_K0", &LVec_n_K0 );
  tree->SetBranchAddress( "NeutralBetaCDH", &NeutralBetaCDH );
  tree->SetBranchAddress( "NeutralBetaCDH_beam", &NeutralBetaCDH_beam );//from v192
  tree->SetBranchAddress( "NeutralBetaCDH_vtx[3]", NeutralBetaCDH_vtx );
  tree->SetBranchAddress( "tofpim",&tofpim);
  tree->SetBranchAddress( "tofpip",&tofpip);
  tree->SetBranchAddress( "tofn",&tofn);
  tree->SetBranchAddress( "dE", &dE );
  tree->SetBranchAddress( "neutralseg", &neutralseg );  
  tree->SetBranchAddress( "nhitOutCDC", &nhitOutCDC ); //charge veto by Outer 3 layer of 3cdc
  tree->SetBranchAddress( "ForwardCharge", &ForwardCharge);
  tree->SetBranchAddress( "vtx_reaction", &vtx_reaction );
  tree->SetBranchAddress( "vtx_displaced", &vtx_displaced );
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
  if(SimSpmode || SimSmmode || SimK0nnmode) {
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
    tree->SetBranchAddress( "mc_vtx", &mc_vtx);
    tree->SetBranchAddress( "mc_disvtx", &mc_disvtx);
  }
  
  TString outname = std::string(filename);
  //cut5 : remove missing n from event mixing
  outname.Replace(std::string(filename).size()-5,5,"_MIX_cut4.root");


  TFile *fout = new TFile(outname.Data(),"RECREATE");
  fout->Print();
  fout->cd();
   
  TTree *treeMIX = new TTree("EventTree","EventTreeMIX");
  treeMIX->Branch( "mom_beam",   &LVec_beam2 );//
  treeMIX->Branch( "mom_beam_Sp",  &LVec_beam_Sp2 );//
  treeMIX->Branch( "mom_beam_Sm",  &LVec_beam_Sm2 );
  treeMIX->Branch( "mom_beam_K0",  &LVec_beam_K02 );
  treeMIX->Branch( "mom_target", &LVec_target2 );
  treeMIX->Branch( "mom_pip", &LVec_pip2 );
  treeMIX->Branch( "mom_pim", &LVec_pim2 );
  treeMIX->Branch( "mom_n", &LVec_n2 );
  treeMIX->Branch( "mom_n_beam", &LVec_n_beam2 );
  treeMIX->Branch( "mom_n_Sp", &LVec_n_Sp2 );
  treeMIX->Branch( "mom_n_Sm", &LVec_n_Sm2 );
  treeMIX->Branch( "mom_n_K0", &LVec_n_K02 );
  treeMIX->Branch( "NeutralBetaCDH", &NeutralBetaCDH2 );
  treeMIX->Branch( "NeutralBetaCDH_beam", &NeutralBetaCDH_beam2 );
  treeMIX->Branch( "NeutralBetaCDH_vtx[3]", NeutralBetaCDH_vtx2);
  treeMIX->Branch( "tofpim",&tofpim2);
  treeMIX->Branch( "tofpip",&tofpip2);
  treeMIX->Branch( "tofn",&tofn2);
  treeMIX->Branch( "dE", &dE2 );
  treeMIX->Branch( "neutralseg", &neutralseg2 );
  treeMIX->Branch( "nhitOutCDC", &nhitOutCDC2 ); //charge veto by Outer 3 layer of 3cdc
  treeMIX->Branch( "ForwardCharge", &ForwardCharge2);
  treeMIX->Branch( "vtx_reaction", &vtx_reaction2 );
  treeMIX->Branch( "vtx_displaced", &vtx_displaced2 );
  treeMIX->Branch( "vtx_pip_beam", &vtx_pip_beam2 );
  treeMIX->Branch( "vtx_pim_beam", &vtx_pim_beam2 );
  treeMIX->Branch( "vtx_pip_cdc", &vtx_pip_cdc2 );
  treeMIX->Branch( "vtx_pim_cdc", &vtx_pim_cdc2 );
  treeMIX->Branch( "CA_pip",&CA_pip2);
  treeMIX->Branch( "CA_pim",&CA_pim2);
  treeMIX->Branch( "CDH_Pos",&CDH_Pos2);
  treeMIX->Branch( "CDH_Pos_pim",&CDH_Pos_pim2);
  treeMIX->Branch( "CDH_Pos_pip",&CDH_Pos_pip2);
  if(SimSpmode || SimSmmode || SimK0nnmode) {
    treeMIX->Branch( "mcncanvtxr",&mcncanvtxr2);
    treeMIX->Branch( "mcncanvtxz",&mcncanvtxz2);
    treeMIX->Branch( "mcncdsgen", &mcncdsgen2);
    treeMIX->Branch( "mcpattern", &mcpattern2);
    treeMIX->Branch( "mcmom_beam", &mcmom_beam2 );
    treeMIX->Branch( "mcmom_pip", &mcmom_pip2 );
    treeMIX->Branch( "mcmom_pim", &mcmom_pim2 );
    treeMIX->Branch( "mcmom_ncds", &mcmom_ncds2 );
    treeMIX->Branch( "mcmom_nmiss", &mcmom_nmiss2 );
    treeMIX->Branch( "react_nmiss",&react_nmiss2);
    treeMIX->Branch( "react_Sigma",&react_Sigma2);
    treeMIX->Branch( "react_pi",&react_pi2);
    treeMIX->Branch( "mc_vtx",&mc_vtx2);
    treeMIX->Branch( "mc_disvtx",&mc_disvtx2);
  }

  Int_t nevent = tree->GetEntries();
  std::cerr<<"# of events = "<<nevent<<std::endl;
  std::cout << "dE cut:" << anacuts::dE_MIN << std::endl;
  //mixed variables
  std::vector <TLorentzVector> vec_LVec_n;
  std::vector <TLorentzVector> vec_LVec_n_Sp;
  std::vector <TLorentzVector> vec_LVec_n_Sm;
  std::vector <TLorentzVector> vec_LVec_n_K0;
  std::vector <double> vec_dE;
  std::vector <double> vec_NeutralBetaCDH;
  std::vector <TVector3> vec_CDH_Pos;
  std::vector <int> vec_neutralseg;
  std::vector <double> vec_tofn;

  for ( Int_t i=0; i<nevent; i++ ) {
    tree->GetEvent(i);
    if( (*LVec_n).P()<anacuts::nmomcut) continue;

    bool NBetaOK=false;
    bool NdEOK=false;
  
    if(anacuts::beta_MIN<NeutralBetaCDH &&  NeutralBetaCDH<anacuts::beta_MAX  ) NBetaOK=true;
    if(anacuts::dE_MIN<dE) NdEOK=true;
  
    TVector3 diffpim = (*CDH_Pos)-(*CDH_Pos_pim);
    double diffPhinpim = (*CDH_Pos).Phi()-(*CDH_Pos_pim).Phi();
    if(diffPhinpim<-1.0*TMath::Pi()) diffPhinpim += 2.0*TMath::Pi();
    else if(diffPhinpim>1.0*TMath::Pi()) diffPhinpim -= 2.0*TMath::Pi();
    if(diffPhinpim<0){
      if( pow((diffPhinpim-anacuts::Isonpim_shift)/anacuts::Isonpim_phicut_left,2.0)+pow(diffpim.Z()/anacuts::Isonpim_zcut,2.0) <1 ) continue;
    }else{
      if( pow((diffPhinpim-anacuts::Isonpim_shift)/anacuts::Isonpim_phicut_right,2.0)+pow(diffpim.Z()/anacuts::Isonpim_zcut,2.0) <1 ) continue;
    }
    if( -anacuts::CDHwidthphi< diffPhinpim  && diffPhinpim < anacuts::CDHwidthphi ) continue;

    TVector3 diffpip = (*CDH_Pos)-(*CDH_Pos_pip);
    double diffPhinpip = (*CDH_Pos).Phi()-(*CDH_Pos_pip).Phi();
    if(diffPhinpip<-1.0*TMath::Pi()) diffPhinpip += 2.0*TMath::Pi();
    else if(diffPhinpip>1.0*TMath::Pi()) diffPhinpip -= 2.0*TMath::Pi();
    if(diffPhinpip<0){
      if( pow((diffPhinpip-anacuts::Isonpip_shift)/anacuts::Isonpip_phicut_left,2.0)+pow(diffpip.Z()/anacuts::Isonpip_zcut,2.0) <1 ) continue;
    }else{
      if( pow((diffPhinpip-anacuts::Isonpip_shift)/anacuts::Isonpip_phicut_right,2.0)+pow(diffpip.Z()/anacuts::Isonpip_zcut,2.0) <1 ) continue;
    }
    //for mixed events, avoid sharing same CDH segments
    if( -anacuts::CDHwidthphi< diffPhinpip  && diffPhinpip < anacuts::CDHwidthphi ) continue;
    
    //if((nhitOutCDC!=0) ) continue;

    if(NBetaOK && NdEOK){
      vec_LVec_n.push_back(*LVec_n);
      vec_LVec_n_Sp.push_back(*LVec_n_Sp);
      vec_LVec_n_Sm.push_back(*LVec_n_Sm);
      vec_LVec_n_K0.push_back(*LVec_n_K0);
      vec_NeutralBetaCDH.push_back(NeutralBetaCDH);
      vec_tofn.push_back(tofn);
      vec_dE.push_back(dE);
      vec_neutralseg.push_back(neutralseg);
      vec_CDH_Pos.push_back(*CDH_Pos);
    }
  }

  decltype(vec_LVec_n)::iterator last_n = vec_LVec_n.end();
  decltype(vec_LVec_n_Sp)::iterator last_n_Sp = vec_LVec_n_Sp.end();
  decltype(vec_LVec_n_Sm)::iterator last_n_Sm = vec_LVec_n_Sm.end();
  decltype(vec_LVec_n_K0)::iterator last_n_K0 = vec_LVec_n_K0.end();
  decltype(vec_NeutralBetaCDH)::iterator last_Beta =vec_NeutralBetaCDH.end();
  decltype(vec_tofn)::iterator last_tofn = vec_tofn.end();
  decltype(vec_dE)::iterator last_dE = vec_dE.end();
  decltype(vec_neutralseg)::iterator last_neutralseg = vec_neutralseg.end();
  decltype(vec_CDH_Pos)::iterator last_CDH_Pos = vec_CDH_Pos.end();
  
  const size_t nsize=vec_LVec_n.size();
  const int npipipairused = 20;
  for ( Int_t i=0; i<nevent*npipipairused; i++ ) {
    tree->GetEvent(i%(nevent-1));
    if(i%50000==0) std::cout << "Event# " << i << std::endl;
    //decrement iterators 
    last_n--;
    last_n_Sp--;
    last_n_Sm--;
    last_n_K0--;
    last_Beta--;
    last_tofn--;
    last_dE--;
    last_neutralseg--;
    last_CDH_Pos--;
    //std::cout << (*last_n).P() << std::endl;
    //std::cout << *last_Beta << std::endl;
    //std::cout << *last_tofn << std::endl;
    //std::cout << *last_dE << std::endl;
    //std::cout << *last_neutralseg << std::endl;
    //std::cout << (*last_CDH_Pos).Phi() << std::endl;

    *LVec_beam2 =  *LVec_beam;
    *LVec_beam_Sp2 = *LVec_beam_Sp;
    *LVec_beam_Sm2 = *LVec_beam_Sm;
    *LVec_beam_K02 = *LVec_beam_K0;
    *LVec_target2 = *LVec_target;
    *LVec_pip2 = *LVec_pip;
    *LVec_pim2 = *LVec_pim;
    *LVec_n2 = *last_n;//mixing param.
    *LVec_n_Sp2 = *last_n_Sp;//mixing param.
    *LVec_n_Sm2 = *last_n_Sm;//mixing param.
    *LVec_n_K02 = *last_n_K0;//mixing param.
    if((*LVec_n2).P() < 0.0001){
      std::cout << "Something is wrong " << std::endl;
      std::cout << (*LVec_n2).P() << std::endl;
    }
    *LVec_n_beam2 = *LVec_n_beam;
    NeutralBetaCDH2 = NeutralBetaCDH;
    //NeutralBetaCDH2 = *last_Beta;//mixing param.
    tofpim2 = tofpim;  
    tofpip2 = tofpip;  
    //tofn2   = tofn; //mixing param.   
    tofn2   = *last_tofn;    
    dE2 = dE;
    //dE2 = *last_dE; //mixing param.
    neutralseg2 = neutralseg;
    //neutralseg2 = *last_neutralseg; //mixing param.
    nhitOutCDC2 = nhitOutCDC;
    ForwardCharge2 = ForwardCharge;
    *vtx_reaction2 = *vtx_reaction;
    *vtx_displaced2 = *vtx_displaced;
    *vtx_pip_beam2 = *vtx_pip_beam;
    *vtx_pim_beam2 = *vtx_pim_beam;
    *vtx_pip_cdc2 = *vtx_pip_cdc;
    *vtx_pim_cdc2 = *vtx_pim_cdc;
    *CA_pip2 = *CA_pip;
    *CA_pim2 = *CA_pim;
    //*CDH_Pos2 = *CDH_Pos;
    *CDH_Pos2 = *last_CDH_Pos;//mixing param.
    *CDH_Pos_pim2 = *CDH_Pos_pim;
    *CDH_Pos_pip2 = *CDH_Pos_pip;
    if(SimSpmode || SimSmmode || SimK0nnmode){
      mcncanvtxr2 = mcncanvtxr;
      mcncanvtxz2 = mcncanvtxz;
      mcncdsgen2  = mcncdsgen;
      mcpattern2  = mcpattern;
      *mcmom_beam2 = *mcmom_beam;
      *mcmom_pip2 = *mcmom_pip;
      *mcmom_pim2 = *mcmom_pim;
      *mcmom_ncds2 = *mcmom_ncds;
      *mcmom_nmiss2 = *mcmom_nmiss;
      *react_nmiss2 = *react_nmiss;
      *react_Sigma2 = *react_Sigma;
      *react_pi2 = *react_pi;
      *vtx_reaction2 = *vtx_reaction;
      *vtx_displaced2 = *vtx_displaced;
    }
    
    //re-use mixing param.
    if(i!=0 &&  (i%(nsize-1)==0)){
     std::cout << "L." << __LINE__ <<  "  "  << i << std::endl;
     last_n+=(nsize-1);
     last_n_Sp+=(nsize-1);
     last_n_Sm+=(nsize-1);
     last_n_K0+=(nsize-1);
     last_Beta+=(nsize-1);
     last_tofn+=(nsize-1);
     last_dE+=(nsize-1);
     last_neutralseg+=(nsize-1);
     last_CDH_Pos+=(nsize-1);
    }
    treeMIX->Fill();
  }

  treeMIX->Write();
  fout->Close();

  std::cout << "end" << std::endl;

}
