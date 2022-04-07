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
TLorentzVector *LVec_target2=NULL; // 4-momentum(target)
TLorentzVector *LVec_pi2=NULL;    // 4-momentum(pi+)
int chargepi2;
TLorentzVector *LVec_n2=NULL;      // 4-momentum(neutron)
double NeutralBetaCDH2; // velocity of neutral particle on CDH
double tofpi2;
double tofn2;
double dE2;   // energy deposit on CDH
int neutralseg2;
int nhitOutCDC2;
int ForwardCharge2;
TVector3 *vtx_reaction2 = NULL; // vertex(reaction) 
TVector3 *vtx_pi_beam2 = NULL; //C.A.P of pip-beam beam side
TVector3 *vtx_pi_cdc2 = NULL;//C.A.P of pip-beam pip side
TVector3 *CDH_Pos2 = NULL;
TVector3 *CDH_Pos_pi2 = NULL;


void GenEventMixTree_h2(const char* filename = "evanaIMsigma_npi_h2_v14.root")
{
  if(gROOT->GetVersionInt() < 60000){
    std::cout << "Use ROOT6 " << std::endl;
    return;
  }
  TFile *f = new TFile(filename,"READ");
  //TFile *f = new TFile("sim_piSpn_dE0_Al.root");
  TTree *tree = (TTree*)f->Get("EventTree");
  if(tree==0) {
    std::cout << "EventTree is not found " << std::endl;
    std::cout << "end " << std::endl;
    return ;
  }

  tree->SetBranchAddress( "mom_beam",   &LVec_beam );
  tree->SetBranchAddress( "mom_target", &LVec_target );
  tree->SetBranchAddress( "mom_pi", &LVec_pi );
  tree->SetBranchAddress( "chargepi", &chargepi );
  tree->SetBranchAddress( "mom_n", &LVec_n );
  tree->SetBranchAddress( "NeutralBetaCDH", &NeutralBetaCDH );
  tree->SetBranchAddress( "tofpi",&tofpi);
  tree->SetBranchAddress( "tofn",&tofn);
  tree->SetBranchAddress( "dE", &dE );
  tree->SetBranchAddress( "neutralseg", &neutralseg );  
  tree->SetBranchAddress( "nhitOutCDC", &nhitOutCDC ); //charge veto by Outer 3 layer of 3cdc
  tree->SetBranchAddress( "ForwardCharge", &ForwardCharge);
  tree->SetBranchAddress( "vtx_reaction", &vtx_reaction );
  tree->SetBranchAddress( "vtx_pi_beam",&vtx_pi_beam);
  tree->SetBranchAddress( "vtx_pi_cdc",&vtx_pi_cdc);
  tree->SetBranchAddress( "CDH_Pos",&CDH_Pos);
  tree->SetBranchAddress( "CDH_Pos_pi",&CDH_Pos_pi);//
  //tree->SetBranchAddress( "run_num", &run_num );
  //tree->SetBranchAddress( "event_num", &event_num );
  //tree->SetBranchAddress( "block_num", &block_num );
  
  TString outname = std::string(filename);
  //cut5 : remove missing n from event mixing
  outname.Replace(std::string(filename).size()-5,5,"_MIX.root");


  TFile *fout = new TFile(outname.Data(),"RECREATE");
  fout->Print();
  fout->cd();
   
  TTree *treeMIX = new TTree("EventTree","EventTreeMIX");
  treeMIX->Branch( "mom_beam",   &LVec_beam2 );//
  treeMIX->Branch( "mom_target", &LVec_target2 );
  treeMIX->Branch( "mom_pi", &LVec_pi2 );
  treeMIX->Branch( "chargepi", &chargepi2 );
  treeMIX->Branch( "mom_n", &LVec_n2 );
  treeMIX->Branch( "NeutralBetaCDH", &NeutralBetaCDH2 );
  treeMIX->Branch( "tofpi",&tofpi2);
  treeMIX->Branch( "tofn",&tofn2);
  treeMIX->Branch( "dE", &dE2 );
  treeMIX->Branch( "neutralseg", &neutralseg2 );
  treeMIX->Branch( "nhitOutCDC", &nhitOutCDC2 ); //charge veto by Outer 3 layer of 3cdc
  treeMIX->Branch( "ForwardCharge", &ForwardCharge2);
  treeMIX->Branch( "vtx_reaction", &vtx_reaction2 );
  treeMIX->Branch( "vtx_pi_beam", &vtx_pi_beam2 );
  treeMIX->Branch( "vtx_pi_cdc", &vtx_pi_cdc2 );
  treeMIX->Branch( "CDH_Pos",&CDH_Pos2);
  treeMIX->Branch( "CDH_Pos_pi",&CDH_Pos_pi2);

  Int_t nevent = tree->GetEntries();
  std::cerr<<"# of events = "<<nevent<<std::endl;
  std::cout << "dE cut:" << anacuts::dE_MIN << std::endl;
  //mixed variables
  std::vector <TLorentzVector> vec_LVec_n;
  std::vector <double> vec_dE;
  std::vector <double> vec_NeutralBetaCDH;
  std::vector <TVector3> vec_CDH_Pos;
  std::vector <int> vec_neutralseg;
  std::vector <double> vec_tofn;

  for ( Int_t i=0; i<nevent; i++ ) {
    tree->GetEvent(i);
    //if( (*LVec_n).P()<anacuts::nmomcut) continue;

    bool NBetaOK=false;
    bool NdEOK=false;
  
    if(anacuts::beta_MIN<NeutralBetaCDH &&  NeutralBetaCDH<anacuts::beta_MAX  ) NBetaOK=true;
    if(anacuts::dE_MIN<dE) NdEOK=true;

    TVector3 diffpim = (*CDH_Pos)-(*CDH_Pos_pi);
    double diffPhinpim = (*CDH_Pos).Phi()-(*CDH_Pos_pi).Phi();
    if(diffPhinpim<-1.0*TMath::Pi()) diffPhinpim += 2.0*TMath::Pi();
    else if(diffPhinpim>1.0*TMath::Pi()) diffPhinpim -= 2.0*TMath::Pi();
    
    /*
    if(chargepi==0){
      if( pow((diffPhinpim-anacuts::Isonpim_shift)/anacuts::Isonpim_phicut,2.0)+pow(diffpim.Z()/anacuts::Isonpim_zcut,2.0) <1 ) continue;
    //for mixed events, avoid sharing same CDH segments
      if( -anacuts::CDHwidthphi< diffPhinpim  && diffPhinpim < anacuts::CDHwidthphi ) continue;
    }*/

    TVector3 diffpip = (*CDH_Pos)-(*CDH_Pos_pi);
    double diffPhinpip = (*CDH_Pos).Phi()-(*CDH_Pos_pi).Phi();
    if(diffPhinpip<-1.0*TMath::Pi()) diffPhinpip += 2.0*TMath::Pi();
    else if(diffPhinpip>1.0*TMath::Pi()) diffPhinpip -= 2.0*TMath::Pi();
    /*
    if(chargepi==1){
      if( pow((diffPhinpip-anacuts::Isonpip_shift)/anacuts::Isonpip_phicut,2.0)+pow(diffpip.Z()/anacuts::Isonpip_zcut,2.0) <1 ) continue;
      if( -anacuts::CDHwidthphi< diffPhinpip  && diffPhinpip < anacuts::CDHwidthphi ) continue;
    }*/

    if((chargepi==0)) {
      if( pow((diffPhinpim-anacuts::Isonpim_shift)/anacuts::Isonpim_phicut_left,2.0)+pow(diffpim.Z()/anacuts::Isonpim_zcut,2.0) <1 ){
        continue;
      }else{
        if( pow((diffPhinpim-anacuts::Isonpim_shift)/anacuts::Isonpim_phicut_right,2.0)+pow(diffpim.Z()/anacuts::Isonpim_zcut,2.0) <1 ){
          continue;
        }
      }
    }else if(chargepi==1) {
      //round cut
      if( pow((diffPhinpip-anacuts::Isonpip_shift)/anacuts::Isonpip_phicut_left,2.0)+pow(diffpip.Z()/anacuts::Isonpip_zcut,2.0) <1 ){
        continue;
      }else{
        if( pow((diffPhinpip-anacuts::Isonpip_shift)/anacuts::Isonpip_phicut_right,2.0)+pow(diffpip.Z()/anacuts::Isonpip_zcut,2.0) <1 ){
          continue;
        }
      }
    }
   

    if(NBetaOK && NdEOK){
      vec_LVec_n.push_back(*LVec_n);
      vec_NeutralBetaCDH.push_back(NeutralBetaCDH);
      vec_tofn.push_back(tofn);
      vec_dE.push_back(dE);
      vec_neutralseg.push_back(neutralseg);
      vec_CDH_Pos.push_back(*CDH_Pos);
    }
  }

  decltype(vec_LVec_n)::iterator last_n = vec_LVec_n.end();
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
    *LVec_target2 = *LVec_target;
    *LVec_pi2 = *LVec_pi;
    chargepi2 = chargepi;
    *LVec_n2 = *last_n;//mixing param.
    if((*LVec_n2).P() < 0.0001){
      std::cout << "Something is wrong " << std::endl;
      std::cout << (*LVec_n2).P() << std::endl;
    }
    NeutralBetaCDH2 = NeutralBetaCDH;
    //NeutralBetaCDH2 = *last_Beta;//mixing param.
    tofpi2 = tofpi;  
    tofn2   = *last_tofn;    
    dE2 = dE;
    //dE2 = *last_dE; //mixing param.
    neutralseg2 = neutralseg;
    //neutralseg2 = *last_neutralseg; //mixing param.
    nhitOutCDC2 = nhitOutCDC;
    ForwardCharge2 = ForwardCharge;
    *vtx_reaction2 = *vtx_reaction;
    *vtx_pi_beam2 = *vtx_pi_beam;
    *vtx_pi_cdc = *vtx_pi_cdc;
    //*CDH_Pos2 = *CDH_Pos;
    *CDH_Pos2 = *last_CDH_Pos;//mixing param.
    *CDH_Pos_pi2 = *CDH_Pos_pi;
    
    //re-use mixing param.
    if(i!=0 &&  (i%(nsize-1)==0)){
      std::cout << "L." << __LINE__ <<  "  "  << i << std::endl;
      last_n+=(nsize-1);
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
