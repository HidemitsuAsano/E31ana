#include <iostream>
#include <vector>
#include <string>
#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TMultiGraph.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TLegend.h>
#include <TPDF.h>
#include <TPaveText.h>
#include <TLatex.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TEfficiency.h>

#include "../src/GlobalVariables.h"
#include "anacuts.h"
#include "globalana.h"

const double pvalcut = 0.005;
const bool gridon=true;
const bool staton=false;
const bool UseKinFit = false;
const bool UseKinFitVal = true;
//const bool Sim1400Cut = false;

//0: diagonal cut
//1: 3 sigma cut
//2: 5 simga cut 
const unsigned int sigmacuttype=0;

//0:diagonal cut
//1:3 sigma cut
//2:5 sigma cut
const unsigned int sidebandtype=0;


void plot_IMpisigma(const char* filename="",const int qvalcutflag=0)
{

  //gROOT->SetStyle("Plain");
  if(staton)gStyle->SetOptStat(111111);
  else gStyle->SetOptStat(0);
  gStyle->SetOptFit(111111);
  gStyle->SetPadGridX(gridon);
  gStyle->SetPadGridY(gridon);
  gStyle->SetStatX(0.9);     
  gStyle->SetStatY(0.9);      
  gStyle->SetPalette(1);
  gStyle->SetStatBorderSize(1);
  gStyle->SetCanvasDefH(800); gStyle->SetCanvasDefW(1000);
  //gStyle->SetTitleFontSize(0.1);
  
  TH1::SetDefaultSumw2();

  std::cout << "infile " << filename <<std::endl;
  TString pdfname = std::string(filename);
  pdfname.Replace(std::string(filename).size()-4,5,"pdf");
  std::cout << "pdfname: " << pdfname << std::endl;
  std::cout << std::endl;
  std::cout << "Use Kin Fit Val ? " << std::endl; 

  if(UseKinFitVal) std::cout << "Yes" << std::endl;
  else             std::cout << "No"  << std::endl;
  
  std::cout << std::endl;
  std::cout << "Sigma selection type     " << sigmacuttype << std::endl;
  std::cout << "Side band selection type " << sidebandtype << std::endl;

  bool Spmode = (std::string(filename).find("Sp")!= std::string::npos);
  bool Smmode = (std::string(filename).find("Sm")!= std::string::npos);
  

  //= = = = pipipnn final-sample tree = = = =//
  
  TFile *f = new TFile(filename);
  //TFile *f = new TFile("sim_piSpn_dE0_Al.root");
  TTree *tree = (TTree*)f->Get("EventTree");
  if(tree==0){
    std::cout << "EventTree is not found " << std::endl;
    return ;
  }
  
  tree->SetBranchAddress( "mom_beam",   &LVec_beam );
  tree->SetBranchAddress( "mom_beam_Sp",   &LVec_beam_Sp );
  tree->SetBranchAddress( "mom_beam_Sm",   &LVec_beam_Sm );
  tree->SetBranchAddress( "mom_target", &LVec_target );
  tree->SetBranchAddress( "mom_pip", &LVec_pip );
  tree->SetBranchAddress( "mom_pim", &LVec_pim );
  tree->SetBranchAddress( "mom_n", &LVec_n );
  tree->SetBranchAddress( "mom_n_Sp", &LVec_n_Sp );
  tree->SetBranchAddress( "mom_n_Sm", &LVec_n_Sm );
  tree->SetBranchAddress( "NeutralBetaCDH", &NeutralBetaCDH );
  tree->SetBranchAddress( "NeutralBetaCDH_vtx[2]", NeutralBetaCDH_vtx );
  tree->SetBranchAddress( "dE", &dE );
  tree->SetBranchAddress( "vtx_reaction", &vtx_reaction );
  tree->SetBranchAddress( "vtx_pip_beam",&vtx_pip_beam);
  tree->SetBranchAddress( "vtx_pim_beam",&vtx_pim_beam);
  tree->SetBranchAddress( "vtx_pip_cdc",&vtx_pip_cdc);
  tree->SetBranchAddress( "vtx_pim_cdc",&vtx_pim_cdc);
  tree->SetBranchAddress( "CA_pip",&CA_pip);
  tree->SetBranchAddress( "CA_pim",&CA_pim);
  //tree->SetBranchAddress( "run_num", &run_num );
  //tree->SetBranchAddress( "event_num", &event_num );
  //tree->SetBranchAddress( "block_num", &block_num );
  if(Spmode || Smmode){
    tree->SetBranchAddress( "mcmom_beam",  &mcmom_beam );
    tree->SetBranchAddress( "mcmom_pip", &mcmom_pip);
    tree->SetBranchAddress( "mcmom_pim", &mcmom_pim);
    tree->SetBranchAddress( "mcmom_ncds", &mcmom_ncds);
    tree->SetBranchAddress( "mcmom_nmiss", &mcmom_nmiss);
  }
  if(UseKinFit){
    tree->SetBranchAddress( "kfSpmode_mom_beam",   &kfSpmode_mom_beam );
    tree->SetBranchAddress( "kfSpmode_mom_pip", &kfSpmode_mom_pip );
    tree->SetBranchAddress( "kfSpmode_mom_pim", &kfSpmode_mom_pim );
    tree->SetBranchAddress( "kfSpmode_mom_n", &kfSpmode_mom_n );
    tree->SetBranchAddress( "kfSpmode_chi2", &kfSpmode_chi2 );
    tree->SetBranchAddress( "kfSpmode_NDF", &kfSpmode_NDF );
    tree->SetBranchAddress( "kfSpmode_status", &kfSpmode_status );
    tree->SetBranchAddress( "kfSpmode_pvalue", &kfSpmode_pvalue );
    tree->SetBranchAddress( "kfSmmode_mom_beam",   &kfSmmode_mom_beam );
    tree->SetBranchAddress( "kfSmmode_mom_pip", &kfSmmode_mom_pip );
    tree->SetBranchAddress( "kfSmmode_mom_pim", &kfSmmode_mom_pim );
    tree->SetBranchAddress( "kfSmmode_mom_n", &kfSmmode_mom_n );
    tree->SetBranchAddress( "kfSmmode_chi2", &kfSmmode_chi2 );
    tree->SetBranchAddress( "kfSmmode_NDF", &kfSmmode_NDF );
    tree->SetBranchAddress( "kfSmmode_status", &kfSmmode_status );
    tree->SetBranchAddress( "kfSmmode_pvalue", &kfSmmode_pvalue );
    tree->SetBranchAddress( "kf_flag", &kf_flag );
  }
  
  // w/o kinematic fit 
  TH2F* dE_betainv_fid;//
  TH2F* dE_MMom_fid_beta_woK0;
  TH2F* dE_MMass_fid_beta_woK0;
  TH2F* MMom_MMass_woK0;
  TH2F* MMom_MMass_woK0_wSid;
  TH2F* IMnpim_IMnpip_dE_woK0;
  TH2F* IMnpim_IMnpip_dE_woK0_n;//
  TH2F* IMnpim_IMnpip_dE_woK0_n_cut;//Sp or Sm mode after selection
  TH2F* IMnpim_IMnpip_dE_woK0_n_Sp;//Spmode, avoiding crossing-point
  TH2F* IMnpim_IMnpip_dE_woK0_n_Sp_bg;//Spmode+background region, avoiding crossing-point
  TH2F* IMnpim_IMnpip_dE_woK0_n_Sm;//Smmode, avoiding crossing-point
  TH2F* IMnpim_IMnpip_dE_woK0_n_Sm_bg;//Smmode, avoiding crossing-point
  TH2F* IMnpim_IMnpip_dE_woK0_n_side;//Side band for Sp + Sm mode
  TH2F* IMnpim_IMnpip_dE_woK0_n_Sp_side[2];//low mass,high mass 
  TH2F* IMnpim_IMnpip_dE_woK0_n_Sm_side[2];//low mass,high mass
  TH2F* IMnpim_IMnpip_dE_woK0_n_Sm_side_cut[2];//low mass,high mass
  TH2F* MMnmiss_IMnpip_dE_woK0;
  TH2F* MMnmiss_IMnpim_dE_woK0;
  TH2F* MMnpip_MMnpim_woK0_wSid_n;
  TH2F* dE_IMnpim_woK0;
  TH2F* dE_IMnpim_woK0_n;
  TH2F* dE_IMnpip_woK0;
  TH2F* dE_IMnpip_woK0_n;
  TH2F* dE_IMnpipi_woK0_wSid_n;
  TH2F* Cosn_IMnpipi_woK0_wSid_n;
  TH2F* MMnmiss_IMnpipi_woK0_wSid_n;
  TH2F* q_IMnpipi_wSid_n;
  TH2F* q_IMnpipi_woK0_wSid_n;
  TH2F* q_IMnpip_gen;//fine bins, no cuts for separating S+/S-
  TH2F* q_IMnpim_gen;//fine bins, no cuts for separating S+/S-
  TH2F* q_IMnpipi_gen;//fine bins, no cuts for separating S+/S-
  TH2F* q_IMnpipi_wSid_n_acc;//fine bins, no cuts for separating S+/S-
  TH2F* q_IMnpipi_wSid_n_acc_reco;//fine bins, no cuts for separating S+/S- ,reconstructed value
  TH2F* q_IMnpipi_woK0_wSid_n_acc;//fine bins, no cuts for separationg S+/S-
  TH2F* q_IMnpipi_woK0_wSid_n_acc_reco;//fine bins, no cuts for separationg S+/S-, reconstructed value
  TH2F* q_IMnpipi_woK0_wSid_n_Sp;
  TH2F* q_IMnpipi_woK0_wSid_n_Sp_acc;//fine bins
  TH2F* q_IMnpipi_woK0_wSid_n_Sp_acc_reco;//fine bins
  TH2F* q_IMnpipi_woK0_wSid_n_Sp_side[3][2];//sideband type, low high side
  TH2F* q_IMnpipi_woK0_wSid_n_Sm;
  TH2F* q_IMnpipi_woK0_wSid_n_Sm_acc;//fine bins
  TH2F* q_IMnpipi_woK0_wSid_n_Sm_acc_reco;//fine bins
  TH2F* q_IMnpipi_woK0_wSid_n_Sm_side[3][2];//sideband type, low high side
  TH2F* IMnpip_IMnpipi_woK0_n;
  TH2F* IMnpim_IMnpipi_woK0_n;

  TH2F* q_IMnpipi_wSid_n_side;//side band method
  TH2F* q_IMnpipi_woK0_wSid_n_side;//side band method
  TH2F* nmom_IMnpipi_woK0_wSid_n;
  TH2F* q_pippim_n;

  const char smode[][4]={"Sp","Sm"};
  
  const int nbinIMnpipi = 100;//1-2 GeV/c^2
  const int nbinq = 25;//0-1.5 GeV/c
  const int nbinIMnpi = 400; //1-2 GeV/c^2
  const int nbinnmiss = 100; //0-1.5 GeV/c
  const int nbindE = 200;
  
  dE_betainv_fid = new TH2F(Form("dE_betainv_fid"),Form("dE_betainv_fid"),1000, 0, 50, nbindE, 0, 50);
  dE_betainv_fid->SetXTitle("1/#beta");
  dE_betainv_fid->SetYTitle("dE [MeVee]");
  
  dE_MMom_fid_beta_woK0 = new TH2F(Form("dE_MMom_fid_beta_woK0"),Form("dE_MMom_fid_beta_woK0"),100, 0, 1.5, nbindE, 0, 50);
  dE_MMom_fid_beta_woK0->SetXTitle("Missing Mom. [GeV/c]");
  dE_MMom_fid_beta_woK0->SetYTitle("dE [MeVee]");

  dE_MMass_fid_beta_woK0 = new TH2F(Form("dE_MMass_fid_beta_woK0"),Form("dE_MMass_fid_beta_woK0"), 140, 0.4, 1.8, nbindE, 0, 50);
  dE_MMass_fid_beta_woK0->SetXTitle("Missing mass [GeV/c^{2}]");
  dE_MMass_fid_beta_woK0->SetYTitle("dE [MeVee]");

  MMom_MMass_woK0 = new TH2F(Form("MMom_MMass_woK0"),Form("MMom_MMass_woK0"), 140, 0.4, 1.8, 100, 0, 1.5);
  MMom_MMass_woK0->SetXTitle("Missing Mass [GeV/c^{2}]");
  MMom_MMass_woK0->SetYTitle("Missing Mom. [GeV/c]");

  MMom_MMass_woK0_wSid = new TH2F(Form("MMom_MMass_woK0_wSid"),Form("MMom_MMass_woK0_wSid"), 140, 0.4, 1.8, 100, 0, 1.5);
  MMom_MMass_woK0_wSid->SetXTitle("Missing Mass [GeV/c^{2}]");
  MMom_MMass_woK0_wSid->SetYTitle("Missing Mom. [GeV/c]");
  
  IMnpim_IMnpip_dE_woK0 = new TH2F(Form("IMnpim_IMnpip_dE_woK0"), Form("IMnpim_IMnpip_dE_woK0"),nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_woK0->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_woK0->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
   
  MMnmiss_IMnpip_dE_woK0 = new TH2F("MMnmiss_IMnpip_dE_woK0", "MMnmiss_IMnpip_dE_woK0",nbinIMnpi,1.,2.0,nbinnmiss,0,1.5);
  MMnmiss_IMnpip_dE_woK0->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  MMnmiss_IMnpip_dE_woK0->SetYTitle("Miss Mass. [GeV/c^{2}]");

  MMnmiss_IMnpim_dE_woK0 = new TH2F("MMnmiss_IMnpim_dE_woK0", "MMnmiss_IMnpim_dE_woK0",nbinIMnpi,1.,2.0,nbinnmiss,0,1.5);
  MMnmiss_IMnpim_dE_woK0->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMnpim_dE_woK0->SetYTitle("Miss Mass. [GeV/c^{2}]");

  IMnpim_IMnpip_dE_woK0_n = new TH2F(Form("IMnpim_IMnpip_dE_woK0_n"),Form("IMnpim_IMnpip_dE_woK0_n"),nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_woK0_n->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_woK0_n->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  
  IMnpim_IMnpip_dE_woK0_n_side = new TH2F(Form("IMnpim_IMnpip_dE_woK0_n_side"),Form("IMnpim_IMnpip_dE_woK0_n_side"),nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_woK0_n_side->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_woK0_n_side->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  
  for(int i=0;i<2;i++){
    const char  lh[][6]={"low","high"};
    IMnpim_IMnpip_dE_woK0_n_Sp_side[i] = new TH2F(Form("IMnpim_IMnpip_dE_woK0_n_Sp_side_%s",lh[i]),Form("IMnpim_IMnpip_dE_woK0_n_Sp_side_%s",lh[i]), nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
    IMnpim_IMnpip_dE_woK0_n_Sp_side[i]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
    IMnpim_IMnpip_dE_woK0_n_Sp_side[i]->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
    
    IMnpim_IMnpip_dE_woK0_n_Sm_side[i] = new TH2F(Form("IMnpim_IMnpip_dE_woK0_n_Sm_side_%s",lh[i]),Form("IMnpim_IMnpip_dE_woK0_n_Sm_side_%s",lh[i]), nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
    IMnpim_IMnpip_dE_woK0_n_Sm_side[i]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
    IMnpim_IMnpip_dE_woK0_n_Sm_side[i]->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
    IMnpim_IMnpip_dE_woK0_n_Sm_side_cut[i] = new TH2F(Form("IMnpim_IMnpip_dE_woK0_n_Sm_side_%s_cut",lh[i]),Form("IMnpim_IMnpip_dE_woK0_n_Sm_side_%s_cut",lh[i]), nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
    IMnpim_IMnpip_dE_woK0_n_Sm_side_cut[i]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
    IMnpim_IMnpip_dE_woK0_n_Sm_side_cut[i]->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  }

  IMnpim_IMnpip_dE_woK0_n_cut = new TH2F(Form("IMnpim_IMnpip_dE_woK0_n_cut"),Form("IMnpim_IMnpip_dE_woK0_n_cut"), nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_woK0_n_cut->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_woK0_n_cut->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  
  IMnpim_IMnpip_dE_woK0_n_Sp = new TH2F(Form("IMnpim_IMnpip_dE_woK0_n_Sp"),Form("IMnpim_IMnpip_dE_woK0_n_Sp"), nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_woK0_n_Sp->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_woK0_n_Sp->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  
  IMnpim_IMnpip_dE_woK0_n_Sp_bg = new TH2F(Form("IMnpim_IMnpip_dE_woK0_n_Sp_bg"),Form("IMnpim_IMnpip_dE_woK0_n_Sp_bg"), nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_woK0_n_Sp_bg->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_woK0_n_Sp_bg->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  
  IMnpim_IMnpip_dE_woK0_n_Sm = new TH2F(Form("IMnpim_IMnpip_dE_woK0_n_Sm"),Form("IMnpim_IMnpip_dE_woK0_n_Sm"), nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_woK0_n_Sm->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_woK0_n_Sm->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  
  IMnpim_IMnpip_dE_woK0_n_Sm_bg = new TH2F(Form("IMnpim_IMnpip_dE_woK0_n_Sm_bg"),Form("IMnpim_IMnpip_dE_woK0_n_Sm_bg"), nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_woK0_n_Sm_bg->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_woK0_n_Sm_bg->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
    
  MMnpip_MMnpim_woK0_wSid_n = new TH2F(Form("MMnpip_MMnpim_woK0_wSid_n"),Form("MMnpip_MMnpim_woK0_wSid_n"),70, 1, 1.7, 70, 1, 1.7);
  MMnpip_MMnpim_woK0_wSid_n->SetXTitle("Miss. Mass(n#pi^{+}) [GeV/c^{2}]");
  MMnpip_MMnpim_woK0_wSid_n->SetYTitle("Miss. Mass(n#pi^{-}) [GeV/c^{2}]");
  
  dE_IMnpim_woK0 = new TH2F(Form("dE_IMnpim_woK0"),Form("dE_IMnpim_woK0"), nbinIMnpi, 1.0, 2.0, nbindE, 0, 50.);
  dE_IMnpim_woK0->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  dE_IMnpim_woK0->SetYTitle("dE [MeVee]");
  
  dE_IMnpim_woK0_n = new TH2F(Form("dE_IMnpim_woK0_n"),Form("dE_IMnpim_woK0_n"), nbinIMnpi, 1.0, 2.0, nbindE, 0, 50.);
  dE_IMnpim_woK0_n->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  dE_IMnpim_woK0_n->SetYTitle("dE [MeVee]");

  dE_IMnpip_woK0 = new TH2F(Form("dE_IMnpip_woK0"),Form("dE_IMnpip_woK0"), nbinIMnpi, 1.0, 2.0, nbindE, 0, 50.);
  dE_IMnpip_woK0->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  dE_IMnpip_woK0->SetYTitle("dE [MeVee]");
  
  dE_IMnpip_woK0_n = new TH2F(Form("dE_IMnpip_woK0_n"),Form("dE_IMnpip_woK0_n"), nbinIMnpi, 1.0, 2.0, nbindE, 0, 50.);
  dE_IMnpip_woK0_n->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  dE_IMnpip_woK0_n->SetYTitle("dE [MeVee]");

  dE_IMnpipi_woK0_wSid_n = new TH2F(Form("dE_IMnpipi_woK0_wSid_n"),Form("dE_IMnpipi_woK0_wSid_n"),nbinIMnpipi, 1, 2, nbindE, 0, 50);
  dE_IMnpipi_woK0_wSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  dE_IMnpipi_woK0_wSid_n->SetYTitle("dE [MeVee]");
    
  Cosn_IMnpipi_woK0_wSid_n = new TH2F(Form("Cosn_IMnpipi_woK0_wSid_n"),Form("dE_Cosn_IMnpipi_woK0_wSid_n"),100, 1, 2, 50, -1, 1);
  Cosn_IMnpipi_woK0_wSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  Cosn_IMnpipi_woK0_wSid_n->SetYTitle("cos#theta_{n} (CM)");
    
  MMnmiss_IMnpipi_woK0_wSid_n = new TH2F(Form("MMnmiss_IMnpipi_woK0_wSid_n"),Form("MMnmiss_IMnpipi_woK0_wSid_n"),nbinIMnpipi,1,2,100,0,1.5);
  MMnmiss_IMnpipi_woK0_wSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMnpipi_woK0_wSid_n->SetYTitle("Miss. Mass. [GeV/c^{2}]");
    
  q_IMnpipi_wSid_n = new TH2F(Form("q_IMnpipi_wSid_n"),Form("q_IMnpipi_wSid_n"),nbinIMnpipi,1,2,nbinq,0,1.5);
  q_IMnpipi_wSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wSid_n->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_IMnpipi_woK0_wSid_n = new TH2F(Form("q_IMnpipi_woK0_wSid_n"),Form("q_IMnpipi_woK0_wSid_n"), nbinIMnpipi,1,2, nbinq,0,1.5);
  q_IMnpipi_woK0_wSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_wSid_n->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_IMnpip_gen = new TH2F(Form("q_IMnpip_gen"),Form("q_IMnpip_gen"),500,1,2,300,0,1.5);
  q_IMnpip_gen->SetXTitle("true IM(n#pi^{+}) [GeV/c^{2}]");
  q_IMnpip_gen->SetYTitle("true Mom. Transfer [GeV/c]");
  
  q_IMnpim_gen = new TH2F(Form("q_IMnpim_gen"),Form("q_IMnpim_gen"),500,1,2,300,0,1.5);
  q_IMnpim_gen->SetXTitle("true IM(n#pi^{-}) [GeV/c^{2}]");
  q_IMnpim_gen->SetYTitle("true Mom. Transfer [GeV/c]");
  
  q_IMnpipi_gen = new TH2F(Form("q_IMnpipi_gen"),Form("q_IMnpipi_gen"),500,1,2,300,0,1.5);
  q_IMnpipi_gen->SetXTitle("true IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_gen->SetYTitle("true Mom. Transfer [GeV/c]");
  
  q_IMnpipi_wSid_n_acc = new TH2F(Form("q_IMnpipi_wSid_n_acc"),Form("q_IMnpipi_wSid_n_acc"),500,1,2,300,0,1.5);
  q_IMnpipi_wSid_n_acc->SetXTitle("true IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wSid_n_acc->SetYTitle("true Mom. Transfer [GeV/c]");
  
  q_IMnpipi_wSid_n_acc_reco = new TH2F(Form("q_IMnpipi_wSid_n_acc_reco"),Form("q_IMnpipi_wSid_n_acc_reco"),500,1,2,300,0,1.5);
  q_IMnpipi_wSid_n_acc_reco->SetXTitle("reco. IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wSid_n_acc_reco->SetYTitle("reco. Mom. Transfer [GeV/c]");
  
  q_IMnpipi_woK0_wSid_n_acc = new TH2F(Form("q_IMnpipi_woK0_wSid_n_acc"),Form("q_IMnpipi_woK0_wSid_n_acc"),500,1,2,300,0,1.5);
  q_IMnpipi_woK0_wSid_n_acc->SetXTitle("true IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_wSid_n_acc->SetYTitle("true Mom. Transfer [GeV/c]");
  
  q_IMnpipi_woK0_wSid_n_acc_reco = new TH2F(Form("q_IMnpipi_woK0_wSid_n_acc_reco"),Form("q_IMnpipi_woK0_wSid_n_acc_reco"),500,1,2,300,0,1.5);
  q_IMnpipi_woK0_wSid_n_acc_reco->SetXTitle("reco. IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_wSid_n_acc_reco->SetYTitle("reco. Mom. Transfer [GeV/c]");
  
  q_IMnpipi_woK0_wSid_n_Sp = new TH2F(Form("q_IMnpipi_woK0_wSid_n_Sp"),Form("q_IMnpipi_woK0_wSid_n_Sp"), nbinIMnpipi,1,2, nbinq,0,1.5);
  q_IMnpipi_woK0_wSid_n_Sp->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_wSid_n_Sp->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_IMnpipi_woK0_wSid_n_Sp_acc = new TH2F(Form("q_IMnpipi_woK0_wSid_n_Sp_acc"),Form("q_IMnpipi_woK0_wSid_n_Sp_acc"),500,1,2,300,0,1.5);
  q_IMnpipi_woK0_wSid_n_Sp_acc->SetXTitle("true IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_wSid_n_Sp_acc->SetYTitle("true Mom. Transfer [GeV/c]");
  
  q_IMnpipi_woK0_wSid_n_Sp_acc_reco = new TH2F(Form("q_IMnpipi_woK0_wSid_n_Sp_acc_reco"),Form("q_IMnpipi_woK0_wSid_n_Sp_acc_reco"),500,1,2,300,0,1.5);
  q_IMnpipi_woK0_wSid_n_Sp_acc_reco->SetXTitle("reco. IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_wSid_n_Sp_acc_reco->SetYTitle("reco. Mom. Transfer [GeV/c]");
  
  for(int itype=0;itype<3;itype++){
    for(int ilh=0;ilh<2;ilh++){
      const char  lh[][6]={"low","high"};
      q_IMnpipi_woK0_wSid_n_Sp_side[itype][ilh] = new TH2F(Form("q_IMnpipi_woK0_wSid_n_Sp_side_%d_%s",itype,lh[ilh]),Form("q_IMnpipi_woK0_wSid_n_Sp_side_%d_%s",itype,lh[ilh]), nbinIMnpipi,1,2, nbinq,0,1.5);
      q_IMnpipi_woK0_wSid_n_Sp_side[itype][ilh]->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
      q_IMnpipi_woK0_wSid_n_Sp_side[itype][ilh]->SetYTitle("Mom. Transfer [GeV/c]");
    }
  }

  q_IMnpipi_woK0_wSid_n_Sm = new TH2F(Form("q_IMnpipi_woK0_wSid_n_Sm"),Form("q_IMnpipi_woK0_wSid_n_Sm"), nbinIMnpipi,1,2, nbinq,0,1.5);
  q_IMnpipi_woK0_wSid_n_Sm->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_wSid_n_Sm->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_IMnpipi_woK0_wSid_n_Sm_acc = new TH2F(Form("q_IMnpipi_woK0_wSid_n_Sm_acc"),Form("q_IMnpipi_woK0_wSid_n_Sm_acc"),500,1,2,300,0,1.5);
  q_IMnpipi_woK0_wSid_n_Sm_acc->SetXTitle("true IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_wSid_n_Sm_acc->SetYTitle("true Mom. Transfer [GeV/c]");
  
  q_IMnpipi_woK0_wSid_n_Sm_acc_reco = new TH2F(Form("q_IMnpipi_woK0_wSid_n_Sm_acc_reco"),Form("q_IMnpipi_woK0_wSid_n_Sm_acc_reco"),500,1,2,300,0,1.5);
  q_IMnpipi_woK0_wSid_n_Sm_acc_reco->SetXTitle("reco. IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_wSid_n_Sm_acc_reco->SetYTitle("reco. Mom. Transfer [GeV/c]");
  
  for(int itype=0;itype<3;itype++){
    for(int ilh=0;ilh<2;ilh++){
      const char  lh[][6]={"low","high"};
      q_IMnpipi_woK0_wSid_n_Sm_side[itype][ilh] = new TH2F(Form("q_IMnpipi_woK0_wSid_n_Sm_side_%d_%s",itype,lh[ilh]),Form("q_IMnpipi_woK0_wSid_n_Sm_side_%d_%s",itype,lh[ilh]),nbinIMnpipi,1,2,nbinq,0,1.5);
      q_IMnpipi_woK0_wSid_n_Sm_side[itype][ilh]->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
      q_IMnpipi_woK0_wSid_n_Sm_side[itype][ilh]->SetYTitle("Mom. Transfer [GeV/c]");
    }
  }

  q_IMnpipi_wSid_n_side = new TH2F(Form("q_IMnpipi_wSid_n_side"),Form("q_IMnpipi_wSid_n_side"), nbinIMnpipi,1,2, nbinq,0,1.5);
  q_IMnpipi_wSid_n_side->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wSid_n_side->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_IMnpipi_woK0_wSid_n_side = new TH2F(Form("q_IMnpipi_woK0_wSid_n_side"),Form("q_IMnpipi_woK0_wSid_n_side"), nbinIMnpipi,1,2, nbinq,0,1.5);
  q_IMnpipi_woK0_wSid_n_side->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_wSid_n_side->SetYTitle("Mom. Transfer [GeV/c]");
  

  IMnpip_IMnpipi_woK0_n = new TH2F(Form("IMnpip_IMnpipi_woK0_n"),Form("IMnpip_IMnpipi_woK0_n"),nbinIMnpi,1.0,2.0, nbinIMnpipi,1.0,2.0);
  IMnpip_IMnpipi_woK0_n->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpip_IMnpipi_woK0_n->SetYTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");

  IMnpim_IMnpipi_woK0_n = new TH2F(Form("IMnpim_IMnpipi_woK0_n"),Form("IMnpim_IMnpipi_woK0_n"),nbinIMnpi,1.0,2.0, nbinIMnpipi,1.0,2.0);
  IMnpim_IMnpipi_woK0_n->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  IMnpim_IMnpipi_woK0_n->SetYTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");


  nmom_IMnpipi_woK0_wSid_n = new TH2F(Form("nmom_IMnpipi_woK0_wSid_n"),Form("nmom_IMnpipi_woK0_wSid_n"), nbinIMnpipi,1,2,100,0,1.0);
  nmom_IMnpipi_woK0_wSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  nmom_IMnpipi_woK0_wSid_n->SetYTitle("nmom  [GeV/c]");
   
  q_pippim_n = new TH2F("q_pippim_n","q_pippim_n",500,0,1, 100,0,1.5);
  q_pippim_n->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_pippim_n->SetYTitle("Mom. Transfer [GeV/c]");

  TH1F *nmom = new TH1F("nmom", "nmom", 50, 0, 1.0);
  nmom->SetXTitle("mom. [GeV/c]");
  
  TH2F *dE_nmom = new TH2F("dE_nmom", "dE_nmom", 50, 0, 1.0, 200 , 0, 50);
  dE_nmom->SetXTitle("mom. [GeV/c]");
  dE_nmom->SetYTitle("dE. [MeVee]");
  
  TH1F *mnmom = new TH1F("mnmom", "mnmom", 100, 0, 2.0);
  mnmom->SetXTitle("mom. [GeV/c]");
  
  TH1F *npipmom = new TH1F("npipmom", "npipmom", 150, 0, 3.0);
  npipmom->SetXTitle("mom. [GeV/c]");
  
  TH1F *npimmom = new TH1F("npimmom", "npimmom", 150, 0, 3.0);
  npimmom->SetXTitle("mom. [GeV/c]");
  
  //DCA analysis
  TH1F* DCA_pip_beam = new TH1F("DCA_pip_beam","DCA_pip_beam",3000,0,30);
  DCA_pip_beam->SetXTitle("DCA [cm]");
  TH1F* DCA_pim_beam = new TH1F("DCA_pim_beam","DCA_pim_beam",3000,0,30);
  DCA_pim_beam->SetXTitle("DCA [cm]");
  TH1F* DCA_pip_pim = new TH1F("DCA_pip_pim","DCA_pip_pim",3000,0,30);
  DCA_pip_pim->SetXTitle("DCA [cm]");
  

  //MC info. for resolution evaluation
  TH2F* q_IMnpipi_woK0_wSid_n_mc;
  q_IMnpipi_woK0_wSid_n_mc = new TH2F(Form("q_IMnpipi_woK0_wSid_n_mc"),Form("q_IMnpipi_woK0_wSid_n_mc"),nbinIMnpipi,1,2,nbinq,0,1.5);
  q_IMnpipi_woK0_wSid_n_mc->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_wSid_n_mc->SetYTitle("Mom. Transfer [GeV/c]");
  
  TH2F* q_IMnpipi_woK0_wSid_n_Sp_mc;
  q_IMnpipi_woK0_wSid_n_Sp_mc = new TH2F(Form("q_IMnpipi_woK0_wSid_n_Sp_mc"),Form("q_IMnpipi_woK0_wSid_n_Sp_mc"),nbinIMnpipi,1,2,nbinq,0,1.5);
  q_IMnpipi_woK0_wSid_n_Sp_mc->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_wSid_n_Sp_mc->SetYTitle("Mom. Transfer [GeV/c]");

  TH2F* q_IMnpipi_woK0_wSid_n_Sm_mc;
  q_IMnpipi_woK0_wSid_n_Sm_mc = new TH2F(Form("q_IMnpipi_woK0_wSid_n_Sm_mc"),Form("q_IMnpipi_woK0_wSid_n_Sm_mc"),nbinIMnpipi,1,2,nbinq,0,1.5);
  q_IMnpipi_woK0_wSid_n_Sm_mc->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_wSid_n_Sm_mc->SetYTitle("Mom. Transfer [GeV/c]");
  
  
  TH2F* diff_IMnpipi_woK0_wSid_n;
  diff_IMnpipi_woK0_wSid_n = new TH2F(Form("diff_IMnpipi_woK0_wSid_n"),Form("diff_IMnpipi_woK0_wSid_n"),nbinIMnpipi,1,2,100,-1,1);
  diff_IMnpipi_woK0_wSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  diff_IMnpipi_woK0_wSid_n->SetYTitle("reco. - gen.  [GeV/c^{2}]");

  TH2F* diff_IMnpipi_woK0_wSid_n_Sp;
  diff_IMnpipi_woK0_wSid_n_Sp = new TH2F(Form("diff_IMnpipi_woK0_wSid_n_Sp"),Form("diff_IMnpipi_woK0_wSid_n_Sp"),nbinIMnpipi,1,2,600,-0.3,0.3);
  diff_IMnpipi_woK0_wSid_n_Sp->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  diff_IMnpipi_woK0_wSid_n_Sp->SetYTitle("reco. - gen. [GeV/c^{2}]");
  
  TH2F* diff_q_woK0_wSid_n_Sp;
  diff_q_woK0_wSid_n_Sp = new TH2F(Form("diff_q_woK0_wSid_n_Sp"),Form("diff_q_woK0_wSid_n_Sp"),nbinq,0,1,600,-0.3,0.3);
  diff_q_woK0_wSid_n_Sp->SetXTitle("Mom. Transfer [GeV/c]");
  diff_q_woK0_wSid_n_Sp->SetYTitle("reco. - gen. [GeV/c^]");

  TH2F* diff_IMnpipi_woK0_wSid_n_Sm;
  diff_IMnpipi_woK0_wSid_n_Sm = new TH2F(Form("diff_IMnpipi_woK0_wSid_n_Sm"),Form("diff_IMnpipi_woK0_wSid_n_Sm"),nbinIMnpipi,1,2,600,-0.3,0.3);
  diff_IMnpipi_woK0_wSid_n_Sm->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  diff_IMnpipi_woK0_wSid_n_Sm->SetYTitle("reco. - gen.  [GeV/c^{2}]");

  TH2F* diff_q_woK0_wSid_n_Sm;
  diff_q_woK0_wSid_n_Sm = new TH2F(Form("diff_q_woK0_wSid_n_Sm"),Form("diff_q_woK0_wSid_n_Sm"),nbinq,0,1,600,-0.3,0.3);
  diff_q_woK0_wSid_n_Sm->SetXTitle("Mom. Transfer [GeV/c]");
  diff_q_woK0_wSid_n_Sm->SetYTitle("reco. - gen. [GeV/c^]");

  Int_t nevent = tree->GetEntries();
  std::cerr<<"# of events = "<<nevent<<std::endl;
  std::cout << "p-value cut:" << pvalcut << std::endl; 
  std::cout << "dE cut:" << anacuts::dE_MIN << std::endl; 
  TCanvas *cinfo = new TCanvas("cinfo","info");
  TPaveText *pt = new TPaveText(.05,.05,.95,.7);
  pt->AddText(Form("p-value cut: %f ",pvalcut));
  pt->AddText(Form("dE cut: %f " ,anacuts::dE_MIN));
  pt->AddText(Form("1/beta min.: %f ",1./anacuts::beta_MAX));
  pt->AddText(Form("1/beta max : %f ",1./anacuts::beta_MIN));
  pt->AddText(Form("K^{0} window : %0.3f - %0.3f",anacuts::pipi_MIN,anacuts::pipi_MAX )); 
  pt->AddText(Form("#Sigma^{+} window : %0.3f - %0.3f",anacuts::Sigmap_MIN,anacuts::Sigmap_MAX )); 
  pt->AddText(Form("#Sigma^{-} window : %0.3f - %0.3f",anacuts::Sigmam_MIN,anacuts::Sigmam_MAX )); 
  pt->AddText(Form("miss. n window : %0.3f - %0.3f",anacuts::neutron_MIN,anacuts::neutron_MAX )); 
  pt->Draw(); 
  //------------------------//
  //--- event roop start ---//
  //------------------------//
  for ( Int_t i=0; i<nevent; i++ ) {
    tree->GetEvent(i);
    if(i%50000==0) std::cout << "Event# " << i << std::endl; 
    TVector3 vtx_pip = *vtx_pip_cdc ;
    TVector3 vtx_pim = *vtx_pim_cdc ;
    // calc missing n //
    TLorentzVector LVec_n_miss = *LVec_target+*LVec_beam-*LVec_pip-*LVec_pim-*LVec_n;
    TLorentzVector LVec_n_miss_vtx[2];
    if(UseKinFit){
      if(!UseKinFitVal){
        LVec_n_miss_vtx[0] = *LVec_target+*LVec_beam_Sp-*LVec_pip-*LVec_pim-*LVec_n_Sp;
        LVec_n_miss_vtx[1] = *LVec_target+*LVec_beam_Sm-*LVec_pip-*LVec_pim-*LVec_n_Sm;
      }else{
        LVec_n_miss_vtx[0] = *LVec_target+*kfSpmode_mom_beam-*kfSpmode_mom_pip-*kfSpmode_mom_pim-*kfSpmode_mom_n;
        LVec_n_miss_vtx[1] = *LVec_target+*kfSmmode_mom_beam-*kfSmmode_mom_pip-*kfSmmode_mom_pim-*kfSmmode_mom_n;
      }
    }
    double nmiss_mass = LVec_n_miss.M();
    double nmiss_mass_vtx[2]={LVec_n_miss_vtx[0].M(),LVec_n_miss_vtx[1].M()};
    double nmiss_mom = LVec_n_miss.P();

    // calc cos(theta) of missing n //
    TVector3 boost = (*LVec_target+*LVec_beam).BoostVector();
    TVector3 boost_vtx[2] = {(*LVec_target+*LVec_beam_Sp).BoostVector(),(*LVec_target+*LVec_beam_Sm).BoostVector()} ;
    TLorentzVector LVec_n_miss_CM = LVec_n_miss;
    TLorentzVector LVec_beam_CM = *LVec_beam;
    LVec_n_miss_CM.Boost(-boost);
    LVec_beam_CM.Boost(-boost);
    double cos_n = LVec_n_miss_CM.Vect().Dot(LVec_beam_CM.Vect())/(LVec_n_miss_CM.Vect().Mag()*LVec_beam_CM.Vect().Mag());
    if(Spmode || Smmode){
      TVector3 boost_mc =  (*LVec_target+*mcmom_beam).BoostVector();
    }
    TLorentzVector qkn_mc;
    if(Spmode || Smmode){
      qkn_mc = *mcmom_beam-*mcmom_nmiss;
    }
    TLorentzVector qkn_vtx[2];
    if(UseKinFit){
      if(!UseKinFitVal){ 
        qkn_vtx[0] = *LVec_beam_Sp-LVec_n_miss_vtx[0];
        qkn_vtx[1] = *LVec_beam_Sm-LVec_n_miss_vtx[1];
      }else{
        qkn_vtx[0] = *kfSpmode_mom_beam-LVec_n_miss_vtx[0];
        qkn_vtx[1] = *kfSmmode_mom_beam-LVec_n_miss_vtx[1];
      }
    }
    // calc pi+pi- //
    TLorentzVector LVec_pip_pim = *LVec_pip+*LVec_pim;
    TLorentzVector LVec_pip_pim_mc;
    if(Spmode || Smmode){
      LVec_pip_pim_mc = *mcmom_pip+*mcmom_pim; 
    }
    // calc pi+n //
    TLorentzVector LVec_pip_n = *LVec_pip+*LVec_n;
    TLorentzVector LVec_pip_n_mc;
    if(Spmode || Smmode){
     LVec_pip_n_mc  = *mcmom_pip+*mcmom_ncds;
    }
    TLorentzVector LVec_pip_n_vtx[2];
    if(UseKinFit){
      if(!UseKinFitVal){
        LVec_pip_n_vtx[0] = *LVec_pip+*LVec_n_Sp;
        LVec_pip_n_vtx[1] = *LVec_pip+*LVec_n_Sm;
      }else{
        LVec_pip_n_vtx[0] = *kfSpmode_mom_pip+*kfSpmode_mom_n;
        LVec_pip_n_vtx[1] = *kfSmmode_mom_pip+*kfSmmode_mom_n;
      }
    }
    // calc pi-n //
    TLorentzVector LVec_pim_n = *LVec_pim+*LVec_n;
    TLorentzVector LVec_pim_n_mc;
    if(Spmode || Smmode){
      LVec_pim_n_mc = *mcmom_pim+*mcmom_ncds;
    }
    TLorentzVector LVec_pim_n_vtx[2];
    if(UseKinFit){
      if(!UseKinFitVal){
        LVec_pim_n_vtx[0]= *LVec_pim+*LVec_n_Sp;
        LVec_pim_n_vtx[1]= *LVec_pim+*LVec_n_Sm;
      }else{
        LVec_pim_n_vtx[0]= *kfSpmode_mom_pim+*kfSpmode_mom_n;
        LVec_pim_n_vtx[1]= *kfSmmode_mom_pim+*kfSmmode_mom_n;
      }
    }

    // calc missing Sp //
    TLorentzVector LVec_pip_n_miss = *LVec_target+*LVec_beam-*LVec_pim-*LVec_n;
    
    TLorentzVector LVec_pip_n_miss_mc;
    if(Spmode || Smmode){
      LVec_pip_n_miss_mc = *LVec_target+*mcmom_beam-*mcmom_pim-*mcmom_ncds;
    }
    TLorentzVector LVec_pip_n_miss_vtx[2];
    if(UseKinFit){
      if(!UseKinFitVal){
        LVec_pip_n_miss_vtx[0]= *LVec_target+*LVec_beam_Sp-*LVec_pim-*LVec_n_Sp;                          
        LVec_pip_n_miss_vtx[1]= *LVec_target+*LVec_beam_Sm-*LVec_pim-*LVec_n_Sm;
      }else{
        LVec_pip_n_miss_vtx[0]= *LVec_target+*kfSpmode_mom_beam-*kfSpmode_mom_pim-*kfSpmode_mom_n;                          
        LVec_pip_n_miss_vtx[1]= *LVec_target+*kfSmmode_mom_beam-*kfSmmode_mom_pim-*kfSmmode_mom_n;
      }
    }

    // calc missing Sm //
    TLorentzVector LVec_pim_n_miss = *LVec_target+*LVec_beam-*LVec_pip-*LVec_n;
    TLorentzVector LVec_pim_n_miss_mc; 
    if(Spmode || Smmode){
      LVec_pim_n_miss_mc = *LVec_target+*mcmom_beam-*mcmom_pip-*mcmom_ncds;
    }
    TLorentzVector LVec_pim_n_miss_vtx[2];
    if(UseKinFit){
      if(!UseKinFitVal){
        LVec_pim_n_miss_vtx[0] = *LVec_target+*LVec_beam_Sp-*LVec_pip-*LVec_n_Sp;
        LVec_pim_n_miss_vtx[1] = *LVec_target+*LVec_beam_Sm-*LVec_pip-*LVec_n_Sm;   
      }else{
        LVec_pim_n_miss_vtx[0] = *LVec_target+*kfSpmode_mom_beam-*kfSpmode_mom_pip-*kfSpmode_mom_n;
        LVec_pim_n_miss_vtx[1] = *LVec_target+*kfSmmode_mom_beam-*kfSmmode_mom_pip-*kfSmmode_mom_n;     
      }
    }
    
    // calc pi+pi-n //
    TLorentzVector LVec_pip_pim_n = *LVec_pip+*LVec_pim+*LVec_n;
    TLorentzVector LVec_pip_pim_n_mc; 
    if(Spmode || Smmode){
      LVec_pip_pim_n_mc = *mcmom_pip+*mcmom_pim+*mcmom_ncds;
    }
    TLorentzVector LVec_pip_pim_n_vtx[2];
    if(UseKinFit){
      if(!UseKinFitVal){
        LVec_pip_pim_n_vtx[0] = *LVec_pip+*LVec_pim+*LVec_n_Sp;
        LVec_pip_pim_n_vtx[1] = *LVec_pip+*LVec_pim+*LVec_n_Sm;
      }else{
        LVec_pip_pim_n_vtx[0] = *kfSpmode_mom_pip+*kfSpmode_mom_pim+*kfSpmode_mom_n;
        LVec_pip_pim_n_vtx[1] = *kfSmmode_mom_pip+*kfSmmode_mom_pim+*kfSmmode_mom_n;
      }
    }
    //if(Sim1400Cut && (Spmode || Smmode)){
    //  if(!(((0.85<cos_n) && (cos_n<=1)) || (LVec_pip_pim_n.M()<1.40))) continue;  
    // }
    TLorentzVector qkn = *LVec_beam-LVec_n_miss;
    TLorentzVector LVec_pip_pim_n_CM = LVec_pip_pim_n;
    LVec_pip_pim_n_CM.Boost(-boost);
    //double cos_X = LVec_pip_pim_n_CM.Vect().Dot(LVec_beam_CM.Vect())/(LVec_pip_pim_n_CM.Vect().Mag()*LVec_beam_CM.Vect().Mag());
    
    if( (qkn.P()>=anacuts::qvalcut) && (qvalcutflag==1) ) continue;
    if( (qkn.P()<anacuts::qvalcut) && (qvalcutflag==2) ) continue;
    //if(qkn.P()>0.70 ) continue;
    //if(LVec_pip_pim_n.M() < 1.45) continue;
    //if(LVec_pip_pim_n.M() > 1.55) continue;
    //if(dcapippim < 1) continue;
    //if(LVec_pip_pim_n.M()<1.45 ) continue;
    //double chi2 = kfSpmode_chi2<kfSmmode_chi2 ? kfSpmode_chi2:kfSmmode_chi2;
    double pvalue = -9999;
    if(UseKinFit) pvalue = kfSmmode_pvalue<kfSpmode_pvalue ? kfSpmode_pvalue:kfSmmode_pvalue;
    
    //Filling generated info.

    if(Spmode || Smmode){
      q_IMnpip_gen->Fill(LVec_pip_n_mc.M(),qkn_mc.P());
      q_IMnpim_gen->Fill(LVec_pim_n_mc.M(),qkn_mc.P());
      q_IMnpipi_gen->Fill(LVec_pip_pim_n_mc.M(),qkn_mc.P());
    }

    bool K0rejectFlag=false;
    bool MissNFlag=false;
    bool NBetaOK=false;
    bool NdEOK=false;
    bool SigmaPFlag=false;
    bool SigmaMFlag=false;
    bool SigmawidePFlag=false;
    bool SigmawideMFlag=false;
    bool SigmaPsideFlag[3]={false,false,false};
    bool SigmaPsideLowFlag=false;
    bool SigmaPsideHighFlag=false;
    bool SigmaMsideFlag[3]={false,false,false};
    bool SigmaMsideLowFlag=false;
    bool SigmaMsideHighFlag=false;
    
    //triangular cuts to minimize acceptance distortion
    bool SigmaCrossFlagTop=false;
    bool SigmaCrossFlagBottom=false;
    bool SigmaCrossFlagLeft=false;
    bool SigmaCrossFlagRight=false;
    
    //triangular cuts for side band of Sp mode low mass side
    bool SigmaCrossPsideLowFlagTop=false;
    bool SigmaCrossPsideLowFlagBottom=false;
    bool SigmaCrossPsideLowFlagLeft=false;
    bool SigmaCrossPsideLowFlagRight=false;
    
    //triangular cuts for side band of Sp mode high mass side
    bool SigmaCrossPsideHighFlagTop=false;
    bool SigmaCrossPsideHighFlagBottom=false;
    bool SigmaCrossPsideHighFlagLeft=false;
    bool SigmaCrossPsideHighFlagRight=false;

    //triangular cuts for side band of Sm mode low mass side
    bool SigmaCrossMsideLowFlagTop=false;
    bool SigmaCrossMsideLowFlagBottom=false;
    bool SigmaCrossMsideLowFlagLeft=false;
    bool SigmaCrossMsideLowFlagRight=false;
    
    //triangular cuts for side band of Sm mode high mass side
    bool SigmaCrossMsideHighFlagTop=false;
    bool SigmaCrossMsideHighFlagBottom=false;
    bool SigmaCrossMsideHighFlagLeft=false;
    bool SigmaCrossMsideHighFlagRight=false;
  
    double dca_pip_beam = (*vtx_pip_beam-*vtx_pip_cdc).Mag();
    double dca_pim_beam = (*vtx_pim_beam-*vtx_pim_cdc).Mag();
    double dca_pip_pim =(*CA_pip-*CA_pim).Mag();
    //-- neutron-ID, K0 and missing neutron selection --//
    if(anacuts::beta_MIN<NeutralBetaCDH &&  NeutralBetaCDH<anacuts::beta_MAX  ) NBetaOK=true;
    if(anacuts::dE_MIN<dE) NdEOK=true;
    double MassNPip= (*LVec_n+*LVec_pip).M();
    double MassNPim= (*LVec_n+*LVec_pim).M();

    //Sigma+ production in CDS
    //band cut for signal
    if( (anacuts::Sigmap_MIN<MassNPip && MassNPip<anacuts::Sigmap_MAX)) SigmaPFlag=true;
        
    //Sigma- production in CDS
    //band cut for signal
    if( (anacuts::Sigmam_MIN<MassNPim && MassNPim<anacuts::Sigmam_MAX)) SigmaMFlag=true;
    
    //Sigma+ production in CDS
    //
    if( (anacuts::Sigmap_MIN_wide<MassNPip && MassNPip<anacuts::Sigmap_MAX_wide)) SigmawidePFlag=true;
        
    //Sigma- production in CDS
    if( (anacuts::Sigmam_MIN_wide<MassNPim && MassNPim<anacuts::Sigmam_MAX_wide)) SigmawideMFlag=true;
      
    //
    //diagonal cuts
    //
    //Sigma cross region top
    if(  (MassNPim >=  (MassNPip - anacuts::Sigmap_center) + anacuts::Sigmam_center )  
      && (MassNPim > -(MassNPip - anacuts::Sigmap_center) + anacuts::Sigmam_center ) ) SigmaCrossFlagTop=true;
    
    //Sigma cross region bottom
    if(  (MassNPim <=  (MassNPip - anacuts::Sigmap_center) + anacuts::Sigmam_center )  
      && (MassNPim < -(MassNPip - anacuts::Sigmap_center) + anacuts::Sigmam_center ) ) SigmaCrossFlagBottom=true;

    //Sigma cross region right
    if(  (MassNPim <  (MassNPip - anacuts::Sigmap_center) + anacuts::Sigmam_center )  
      && (MassNPim >= -(MassNPip - anacuts::Sigmap_center) + anacuts::Sigmam_center ) ) SigmaCrossFlagRight=true;
   
    //Sigma cross region left
    if(  (MassNPim >  (MassNPip - anacuts::Sigmap_center) + anacuts::Sigmam_center )  
      && (MassNPim <= -(MassNPip - anacuts::Sigmap_center) + anacuts::Sigmam_center ) ) SigmaCrossFlagLeft=true;
    
    
    //0: diagonal cut
    //1: 3 sigma cut
    //2: 5 simga cut 
    bool SigmaPcutFlag[3] = {false,false,false};
    bool SigmaMcutFlag[3] = {false,false,false};
    
    //diagnal cut
    if(sigmacuttype==0){
      if(SigmaPFlag && !SigmaCrossFlagLeft && !SigmaCrossFlagRight){ 
        SigmaPcutFlag[0]=true;
      }
      if(SigmaMFlag && !SigmaCrossFlagTop  && !SigmaCrossFlagBottom){
        SigmaMcutFlag[0]=true;
      }
    }
    
    if(sigmacuttype==1){
      if(SigmaPFlag && !SigmaMFlag){ 
        SigmaPcutFlag[1]=true;
      }
      if(SigmaMFlag && !SigmaPFlag){
        SigmaMcutFlag[1]=true;
      }
    }

    if(sigmacuttype==2){
      if(SigmaPFlag && !SigmawideMFlag){ 
        SigmaPcutFlag[2]=true;
      }
      if(SigmaMFlag && !SigmawidePFlag){
        SigmaMcutFlag[2]=true;
      }
    }

    //Sigma cross side band Sp low mass top
    if(  (MassNPim >  (MassNPip - anacuts::Sigmap_sidelow_center) + anacuts::Sigmam_center)  
      && (MassNPim > -(MassNPip - anacuts::Sigmap_sidelow_center) + anacuts::Sigmam_center)) SigmaCrossPsideLowFlagTop=true;

    //Sigma cross side band Sp low mass bottom
    if(  (MassNPim <  (MassNPip - anacuts::Sigmap_sidelow_center) + anacuts::Sigmam_center)  
      && (MassNPim < -(MassNPip - anacuts::Sigmap_sidelow_center) + anacuts::Sigmam_center)) SigmaCrossPsideLowFlagBottom=true;
    
    //Sigma cross side band Sp low mass right 
    if(  (MassNPim <  (MassNPip - anacuts::Sigmap_sidelow_center) + anacuts::Sigmam_center)  
      && (MassNPim > -(MassNPip - anacuts::Sigmap_sidelow_center) + anacuts::Sigmam_center)) SigmaCrossPsideLowFlagRight=true;

    //Sigma cross side band Sp low mass left
    if(  (MassNPim >  (MassNPip - anacuts::Sigmap_sidelow_center) + anacuts::Sigmam_center)  
      && (MassNPim < -(MassNPip - anacuts::Sigmap_sidelow_center) + anacuts::Sigmam_center)) SigmaCrossPsideLowFlagLeft=true;

     
    //Sigma cross side band Sp high mass top
    if(  (MassNPim >  (MassNPip - anacuts::Sigmap_sidehigh_center) + anacuts::Sigmam_center)  
      && (MassNPim > -(MassNPip - anacuts::Sigmap_sidehigh_center) + anacuts::Sigmam_center)) SigmaCrossPsideHighFlagTop=true;

    //Sigma cross side band Sp high mass bottom
    if(  (MassNPim <  (MassNPip - anacuts::Sigmap_sidehigh_center) + anacuts::Sigmam_center)  
      && (MassNPim < -(MassNPip - anacuts::Sigmap_sidehigh_center) + anacuts::Sigmam_center)) SigmaCrossPsideHighFlagBottom=true;
    
    //Sigma cross side band Sp high mass right 
    if(  (MassNPim <  (MassNPip - anacuts::Sigmap_sidehigh_center) + anacuts::Sigmam_center)  
      && (MassNPim > -(MassNPip - anacuts::Sigmap_sidehigh_center) + anacuts::Sigmam_center)) SigmaCrossPsideHighFlagRight=true;

    //Sigma cross side band Sp high mass left
    if(  (MassNPim >  (MassNPip - anacuts::Sigmap_sidehigh_center) + anacuts::Sigmam_center)  
      && (MassNPim < -(MassNPip - anacuts::Sigmap_sidehigh_center) + anacuts::Sigmam_center)) SigmaCrossPsideHighFlagLeft=true;


    //Sigma cross side band Sm low mass top
    if(  (MassNPim >  (MassNPip - anacuts::Sigmap_center) + anacuts::Sigmam_sidelow_center)  
      && (MassNPim > -(MassNPip - anacuts::Sigmap_center) + anacuts::Sigmam_sidelow_center)) SigmaCrossMsideLowFlagTop=true;

    //Sigma cross side band Sp low mass bottom
    if(  (MassNPim <  (MassNPip - anacuts::Sigmap_center) + anacuts::Sigmam_sidelow_center)  
      && (MassNPim < -(MassNPip - anacuts::Sigmap_center) + anacuts::Sigmam_sidelow_center)) SigmaCrossMsideLowFlagBottom=true;
    
    //Sigma cross side band Sp low mass right 
    if(  (MassNPim <  (MassNPip - anacuts::Sigmap_center) + anacuts::Sigmam_sidelow_center)  
      && (MassNPim > -(MassNPip - anacuts::Sigmap_center) + anacuts::Sigmam_sidelow_center)) SigmaCrossMsideLowFlagRight=true;

    //Sigma cross side band Sp low mass left
    if(  (MassNPim >  (MassNPip - anacuts::Sigmap_center) + anacuts::Sigmam_sidelow_center)  
      && (MassNPim < -(MassNPip - anacuts::Sigmap_center) + anacuts::Sigmam_sidelow_center)) SigmaCrossMsideLowFlagLeft=true;

     
    //Sigma cross side band Sm high mass top
    if(  (MassNPim >  (MassNPip - anacuts::Sigmap_center) + anacuts::Sigmam_sidehigh_center)  
      && (MassNPim > -(MassNPip - anacuts::Sigmap_center) + anacuts::Sigmam_sidehigh_center)) SigmaCrossMsideHighFlagTop=true;

    //Sigma cross side band Sm high mass bottom
    if(  (MassNPim <  (MassNPip - anacuts::Sigmap_center) + anacuts::Sigmam_sidehigh_center)  
      && (MassNPim < -(MassNPip - anacuts::Sigmap_center) + anacuts::Sigmam_sidehigh_center)) SigmaCrossMsideHighFlagBottom=true;
    
    //Sigma cross side band Sm high mass right 
    if(  (MassNPim <  (MassNPip - anacuts::Sigmap_center) + anacuts::Sigmam_sidehigh_center)  
      && (MassNPim > -(MassNPip - anacuts::Sigmap_center) + anacuts::Sigmam_sidehigh_center)) SigmaCrossMsideHighFlagRight=true;

    //Sigma cross side band Sm high mass left
    if(  (MassNPim >  (MassNPip - anacuts::Sigmap_center) + anacuts::Sigmam_sidehigh_center)  
      && (MassNPim < -(MassNPip - anacuts::Sigmap_center) + anacuts::Sigmam_sidehigh_center)) SigmaCrossMsideHighFlagLeft=true;

    
    //Sigma+ production side band low mass side
    if( ((anacuts::Sigmap_sidelow_MIN<MassNPip) && 
         (MassNPip < anacuts::Sigmap_sidelow_MAX))) SigmaPsideLowFlag=true;
    
    //Sigma+ production side band high mass side
    if( ((anacuts::Sigmap_sidehigh_MIN<MassNPip) && 
         (MassNPip < anacuts::Sigmap_sidehigh_MAX))) SigmaPsideHighFlag=true;
     
    //Sigma+ production side band low or high mass side
    //
    //type 0: diagonal cut
    if( (SigmaPsideLowFlag  && !SigmaCrossPsideLowFlagLeft  && !SigmaCrossPsideLowFlagRight) 
        ||  (SigmaPsideHighFlag && !SigmaCrossPsideHighFlagLeft  && !SigmaCrossPsideHighFlagRight)
      ) SigmaPsideFlag[0]=true;
    
    //type 1: rectangle cut, avoiding signal region
    if( (SigmaPsideLowFlag || SigmaPsideHighFlag) && !SigmaMFlag) SigmaPsideFlag[1]=true;
    
    //type 2: rectangle cut, avoiding signal region + 2 sigma  
    if( (SigmaPsideLowFlag || SigmaPsideHighFlag) && !SigmawideMFlag) SigmaPsideFlag[2]=true;

    //Sigma- production side band low mass side
    if( ((anacuts::Sigmam_sidelow_MIN<MassNPim) && 
         (MassNPim <  anacuts::Sigmam_sidelow_MAX))) SigmaMsideLowFlag=true;
    
    //Sigma- production side band high mass side
    if( ((anacuts::Sigmam_sidehigh_MIN<MassNPim) && 
         (MassNPim <  anacuts::Sigmam_sidehigh_MAX))) SigmaMsideHighFlag=true;
    
    //Sigma- production side band low or high mass side
    if( (SigmaMsideLowFlag  && !SigmaCrossMsideLowFlagTop  && !SigmaCrossMsideLowFlagBottom) 
        ||  (SigmaMsideHighFlag && !SigmaCrossMsideHighFlagTop && !SigmaCrossMsideHighFlagBottom)
      ) SigmaMsideFlag[0]=true;

    //type 1: rectangle cut, avoiding signal region
    if( (SigmaMsideLowFlag || SigmaMsideHighFlag) &&  !SigmaPFlag) SigmaMsideFlag[1]=true;

    //type 2: rectangle cut, avoiding signal region + 2 sigma
    if( (SigmaMsideLowFlag || SigmaMsideHighFlag) &&  !SigmawidePFlag) SigmaMsideFlag[2]=true;


    if(anacuts::neutron_MIN<nmiss_mass && nmiss_mass<anacuts::neutron_MAX ) MissNFlag=true;
    
    //K0 rejection using original momentum
    if( (LVec_pip_pim.M()<anacuts::pipi_MIN || anacuts::pipi_MAX<LVec_pip_pim.M())) K0rejectFlag=true;
    
    bool isrecoPassed = (dE>0);
    //w/o kinfit
    if(K0rejectFlag){
      dE_betainv_fid->Fill(1./NeutralBetaCDH,dE);
    }
    if(K0rejectFlag && NBetaOK){
      dE_MMom_fid_beta_woK0->Fill(LVec_n_miss.P(),dE);
      dE_MMass_fid_beta_woK0->Fill(LVec_n_miss.M(),dE);
      dE_IMnpim_woK0->Fill(LVec_pim_n.M(),dE);
      dE_IMnpip_woK0->Fill(LVec_pip_n.M(),dE);
    }
    if(K0rejectFlag && NBetaOK && MissNFlag){
      dE_IMnpim_woK0_n->Fill(LVec_pim_n.M(),dE);
      dE_IMnpip_woK0_n->Fill(LVec_pip_n.M(),dE);
    }
    if(K0rejectFlag && NBetaOK && NdEOK){
      MMom_MMass_woK0->Fill(LVec_n_miss.M(),LVec_n_miss.P());
      IMnpim_IMnpip_dE_woK0->Fill(LVec_pip_n.M(),LVec_pim_n.M());
      nmom->Fill((*LVec_n).P());
      dE_nmom->Fill((*LVec_n).P(),dE);
      mnmom->Fill(nmiss_mom);
      npipmom->Fill(LVec_pip_n.P());
      npimmom->Fill(LVec_pim_n.P());
      MMnmiss_IMnpip_dE_woK0->Fill(LVec_pip_n.M(),nmiss_mass);
      MMnmiss_IMnpim_dE_woK0->Fill(LVec_pim_n.M(),nmiss_mass);
      if(SigmaPFlag || SigmaMFlag){
        MMom_MMass_woK0_wSid->Fill(LVec_n_miss.M(),LVec_n_miss.P());
      }
      MMnmiss_IMnpipi_woK0_wSid_n->Fill(LVec_pip_pim_n.M(), nmiss_mass);
    }
    if(K0rejectFlag && NBetaOK && NdEOK && MissNFlag){
      IMnpim_IMnpip_dE_woK0_n->Fill(LVec_pip_n.M(),LVec_pim_n.M());
      IMnpip_IMnpipi_woK0_n->Fill(LVec_pip_n.M(),LVec_pip_pim_n.M());
      IMnpim_IMnpipi_woK0_n->Fill(LVec_pim_n.M(),LVec_pip_pim_n.M());
      
      
      if(SigmaPsideFlag[sidebandtype]|| SigmaMsideFlag[sidebandtype]){
        IMnpim_IMnpip_dE_woK0_n_side->Fill(LVec_pip_n.M(),LVec_pim_n.M());
      }
      if(SigmaPsideLowFlag && SigmaPsideFlag[sidebandtype] ){
        IMnpim_IMnpip_dE_woK0_n_Sp_side[0]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
      }
      if(SigmaPsideHighFlag && SigmaPsideFlag[sidebandtype]){
        IMnpim_IMnpip_dE_woK0_n_Sp_side[1]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
      }
      if(SigmaMsideLowFlag && SigmaMsideFlag[sidebandtype]){
        IMnpim_IMnpip_dE_woK0_n_Sm_side[0]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
        IMnpim_IMnpip_dE_woK0_n_Sm_side_cut[0]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
      }
      if(SigmaMsideHighFlag && SigmaMsideFlag[sidebandtype]){
        IMnpim_IMnpip_dE_woK0_n_Sm_side[1]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
        IMnpim_IMnpip_dE_woK0_n_Sm_side_cut[1]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
      }
      
      if(SigmaPFlag || SigmaMFlag){
        MMnpip_MMnpim_woK0_wSid_n->Fill(LVec_pim_n_miss.M(),LVec_pip_n_miss.M());
        dE_IMnpipi_woK0_wSid_n->Fill(LVec_pip_pim_n.M(),dE);
        Cosn_IMnpipi_woK0_wSid_n->Fill(LVec_pip_pim_n.M(),cos_n);
        IMnpim_IMnpip_dE_woK0_n_cut->Fill(LVec_pip_n.M(),LVec_pim_n.M());
        //MMnmiss_IMnpipi_woK0_wSid_n->Fill(LVec_pip_pim_n.M(), nmiss_mass);
        q_IMnpipi_woK0_wSid_n_acc->Fill(LVec_pip_pim_n_mc.M(),qkn_mc.P());
        q_IMnpipi_woK0_wSid_n_acc_reco->Fill(LVec_pip_pim_n.M(),qkn.P());
        q_IMnpipi_woK0_wSid_n->Fill(LVec_pip_pim_n.M(),qkn.P());
        nmom_IMnpipi_woK0_wSid_n->Fill(LVec_pip_pim_n.M(),(*LVec_n).P());
        DCA_pip_beam->Fill( dca_pip_beam);
        DCA_pim_beam->Fill( dca_pim_beam );
        DCA_pip_pim->Fill(dca_pip_pim);
      }
       


      //0: diagonal cut
      //1: 3 sigma cut
      //2: 5 simga cut 
      if(SigmaPcutFlag[sigmacuttype]) {
        IMnpim_IMnpip_dE_woK0_n_Sp->Fill(LVec_pip_n.M(),LVec_pim_n.M());
        q_IMnpipi_woK0_wSid_n_Sp->Fill(LVec_pip_pim_n.M(),qkn.P());
        q_IMnpipi_woK0_wSid_n_Sp_acc->Fill(LVec_pip_pim_n_mc.M(),qkn_mc.P());
        q_IMnpipi_woK0_wSid_n_Sp_acc_reco->Fill(LVec_pip_pim_n.M(),qkn.P());
        if(Spmode){
          q_IMnpipi_woK0_wSid_n_Sp_mc->Fill(LVec_pip_pim_n_mc.M(),qkn_mc.P());
          diff_IMnpipi_woK0_wSid_n_Sp->Fill(LVec_pip_pim_n.M(),LVec_pip_pim_n.M()-LVec_pip_pim_n_mc.M());
          diff_q_woK0_wSid_n_Sp->Fill(qkn.P(),qkn.P()-qkn_mc.P());
        }
      }
      if(SigmaPcutFlag[sigmacuttype] || (!SigmaPFlag)){
        IMnpim_IMnpip_dE_woK0_n_Sp_bg->Fill(LVec_pip_n.M(),LVec_pim_n.M());
      }


      for(int itype=0;itype<3;itype++){
        if(SigmaPsideFlag[itype] && SigmaPsideLowFlag ){
          q_IMnpipi_woK0_wSid_n_Sp_side[itype][0]->Fill(LVec_pip_pim_n.M(),qkn.P());
        }
        if(SigmaPsideFlag[itype] && SigmaPsideHighFlag ){
          q_IMnpipi_woK0_wSid_n_Sp_side[itype][1]->Fill(LVec_pip_pim_n.M(),qkn.P());
        }
      }

      if(SigmaMcutFlag[sigmacuttype]){
        IMnpim_IMnpip_dE_woK0_n_Sm->Fill(LVec_pip_n.M(),LVec_pim_n.M());
        q_IMnpipi_woK0_wSid_n_Sm->Fill(LVec_pip_pim_n.M(),qkn.P());
        //q_IMnpipi_woK0_wSid_n_Sm_acc->Fill(LVec_pip_pim_n.M(),qkn.P());
        q_IMnpipi_woK0_wSid_n_Sm_acc->Fill(LVec_pip_pim_n_mc.M(),qkn_mc.P());
        q_IMnpipi_woK0_wSid_n_Sm_acc_reco->Fill(LVec_pip_pim_n.M(),qkn.P());
        if(Smmode){
          q_IMnpipi_woK0_wSid_n_Sm_mc->Fill(LVec_pip_pim_n_mc.M(),qkn_mc.P());
          diff_IMnpipi_woK0_wSid_n_Sm->Fill(LVec_pip_pim_n.M(),LVec_pip_pim_n.M()-LVec_pip_pim_n_mc.M());
          diff_q_woK0_wSid_n_Sm->Fill(qkn.P(),qkn.P()-qkn_mc.P());
        }
      }
      if(SigmaMcutFlag[sigmacuttype] || (!SigmaMFlag)){
        IMnpim_IMnpip_dE_woK0_n_Sm_bg->Fill(LVec_pip_n.M(),LVec_pim_n.M());
      }
      for(int itype=0;itype<3;itype++){
        if(SigmaMsideFlag[itype] && SigmaMsideLowFlag){
          q_IMnpipi_woK0_wSid_n_Sm_side[itype][0]->Fill(LVec_pip_pim_n.M(),qkn.P());
        }
        if(SigmaMsideFlag[itype] && SigmaMsideHighFlag){
          q_IMnpipi_woK0_wSid_n_Sm_side[itype][1]->Fill(LVec_pip_pim_n.M(),qkn.P());
        }
      }

    }
    if(K0rejectFlag && NBetaOK && NdEOK && MissNFlag){
      if(SigmaPsideFlag[sidebandtype] || SigmaMsideFlag[sidebandtype]){
        q_IMnpipi_woK0_wSid_n_side->Fill(LVec_pip_pim_n.M(),qkn.P());
      }
    }

    
    //including K0 
    if(NBetaOK && NdEOK && MissNFlag){
      q_pippim_n->Fill(LVec_pip_pim.M(),qkn.P());
      if(SigmaPFlag || SigmaMFlag){
        q_IMnpipi_wSid_n->Fill(LVec_pip_pim_n.M(),qkn.P());
        q_IMnpipi_wSid_n_acc->Fill(LVec_pip_pim_n_mc.M(),qkn_mc.P());
        q_IMnpipi_wSid_n_acc_reco->Fill(LVec_pip_pim_n.M(),qkn.P());
      }
      if(SigmaPsideFlag[sidebandtype] || SigmaMsideFlag[sidebandtype]){
        q_IMnpipi_wSid_n_side->Fill(LVec_pip_pim_n.M(),qkn.P());
      }
    }

	}//for ievt
   

  //----------
  //Drawing Part
  //---------
  

  TCanvas *cMMnmiss_IMnpip_dE_woK0 = new TCanvas("cMMnmiss_IMnpip_dE_woK0","MMnmiss_IMnpip_dE_woK0");
  MMnmiss_IMnpip_dE_woK0->GetXaxis()->SetRangeUser(1.05,1.5);
  MMnmiss_IMnpip_dE_woK0->GetYaxis()->SetRangeUser(0.6,1.5);
  MMnmiss_IMnpip_dE_woK0->Draw("colz");

  TCanvas *cMMnmiss_IMnpim_dE_woK0 = new TCanvas("cMMnmiss_IMnpim_dE_woK0","MMnmiss_IMnpim_dE_woK0");
  MMnmiss_IMnpim_dE_woK0->GetXaxis()->SetRangeUser(1.05,1.5);
  MMnmiss_IMnpim_dE_woK0->GetYaxis()->SetRangeUser(0.6,1.5);
  MMnmiss_IMnpim_dE_woK0->Draw("colz");
  
  //Sigma reconstruction 
  TCanvas *cIMnpim_IMnpip_dE_woK0_n = new TCanvas("cIMnpim_IMnpip_dE_woK0_n","IMnpim_IMnpip_dE_woK0_n");
  cIMnpim_IMnpip_dE_woK0_n->cd();
  IMnpim_IMnpip_dE_woK0_n->GetXaxis()->SetRangeUser(0,1.7);
  IMnpim_IMnpip_dE_woK0_n->GetYaxis()->SetRangeUser(0,1.7);
  IMnpim_IMnpip_dE_woK0_n->Draw("colz");
    
  //Sigma+ or Sigma- selection
  TCanvas *cIMnpim_IMnpip_dE_woK0_n_SpSm = new TCanvas("cIMnpim_IMnpip_dE_woK0_n_SpSm","IMnpim_IMnpip_dE_woK0_n_SpSM");
  IMnpim_IMnpip_dE_woK0_n_Sp->GetXaxis()->SetRangeUser(0,1.7);
  IMnpim_IMnpip_dE_woK0_n_Sp->GetYaxis()->SetRangeUser(0,1.7);
  IMnpim_IMnpip_dE_woK0_n_Sp->SetMaximum(IMnpim_IMnpip_dE_woK0_n->GetMaximum());
  IMnpim_IMnpip_dE_woK0_n_Sp->Draw("colz");
  IMnpim_IMnpip_dE_woK0_n_Sm->Draw("colsame");


  TCanvas *cIMnpim_IMnpip_dE_woK0_n_Sp = new TCanvas("cIMnpim_IMnpip_dE_woK0_n_Sp","IMnpim_IMnpip_dE_woK0_n_Sp");
  IMnpim_IMnpip_dE_woK0_n_Sp->SetMaximum(IMnpim_IMnpip_dE_woK0_n->GetMaximum());
  IMnpim_IMnpip_dE_woK0_n_Sp->Draw("colz");
  IMnpim_IMnpip_dE_woK0_n_Sp_side[0]->Draw("colsame");
  IMnpim_IMnpip_dE_woK0_n_Sp_side[1]->Draw("colsame");

  TCanvas *cIMnpim_IMnpip_dE_woK0_n_Sm = new TCanvas("cIMnpim_IMnpip_dE_woK0_n_Sm","IMnpim_IMnpip_dE_woK0_n_Sm");
  IMnpim_IMnpip_dE_woK0_n_Sm->SetMaximum(IMnpim_IMnpip_dE_woK0_n->GetMaximum());
  IMnpim_IMnpip_dE_woK0_n_Sm->GetXaxis()->SetRangeUser(0,1.7);
  IMnpim_IMnpip_dE_woK0_n_Sm->GetYaxis()->SetRangeUser(0,1.7);
  IMnpim_IMnpip_dE_woK0_n_Sm->Draw("colz");
  IMnpim_IMnpip_dE_woK0_n_Sm_side_cut[0]->Draw("colsame");
  IMnpim_IMnpip_dE_woK0_n_Sm_side_cut[1]->Draw("colsame");
  
  TCanvas *cIMnpim_IMnpip_dE_woK0_n_Sp_bg = new TCanvas("cIMnpim_IMnpip_dE_woK0_n_Sp_bg","IMnpim_IMnpip_dE_woK0_n_Sp_bg");
  IMnpim_IMnpip_dE_woK0_n_Sp_bg->SetMaximum(IMnpim_IMnpip_dE_woK0_n->GetMaximum());
  IMnpim_IMnpip_dE_woK0_n_Sp_bg->GetXaxis()->SetRangeUser(0,1.7);
  IMnpim_IMnpip_dE_woK0_n_Sp_bg->GetYaxis()->SetRangeUser(0,1.7);
  IMnpim_IMnpip_dE_woK0_n_Sp_bg->Draw("colz");

  TCanvas *cIMnpim_IMnpip_dE_woK0_n_Sm_bg = new TCanvas("cIMnpim_IMnpip_dE_woK0_n_Sm_bg","IMnpim_IMnpip_dE_woK0_n_Sm_bg");
  IMnpim_IMnpip_dE_woK0_n_Sm_bg->SetMaximum(IMnpim_IMnpip_dE_woK0_n->GetMaximum());
  IMnpim_IMnpip_dE_woK0_n_Sm_bg->GetXaxis()->SetRangeUser(0,1.7);
  IMnpim_IMnpip_dE_woK0_n_Sm_bg->GetYaxis()->SetRangeUser(0,1.7);
  IMnpim_IMnpip_dE_woK0_n_Sm_bg->Draw("colz");
  //
  TCanvas *cIMnpim_IMnpip_dE_woK0_n_px2 = new TCanvas("cIMnpim_IMnpip_dE_woK0_n_px2","IMnpim_IMnpip_dE_woK0_n_px2");
  cIMnpim_IMnpip_dE_woK0_n_px2->cd();
  //TH1D *IMnpim_IMnpip_dE_woK0_n_px = IMnpim_IMnpip_dE_woK0_n->ProjectionX();
  TH1D *IMnpim_IMnpip_dE_woK0_n_px = IMnpim_IMnpip_dE_woK0_n_Sp_bg->ProjectionX();
  IMnpim_IMnpip_dE_woK0_n_px->GetXaxis()->SetRangeUser(1.0,1.3);
  //IMnpim_IMnpip_dE_woK0_n_px->GetYaxis()->SetTitle("Counts/ 10 MeV/c^{2}");
  IMnpim_IMnpip_dE_woK0_n_px->GetYaxis()->CenterTitle();
  IMnpim_IMnpip_dE_woK0_n_px->Draw("EH");
  //
  //
  TH1D* IMnpim_IMnpip_dE_woK0_n_px_1 = (TH1D*) IMnpim_IMnpip_dE_woK0_n_px->Clone();
  int binlow  = IMnpim_IMnpip_dE_woK0_n_px_1->GetXaxis()->FindBin(anacuts::Sigmap_MIN);
  int binhigh = IMnpim_IMnpip_dE_woK0_n_px_1->GetXaxis()->FindBin(anacuts::Sigmap_MAX);
  double binlowcenter  = IMnpim_IMnpip_dE_woK0_n_px_1->GetXaxis()->GetBinCenter(binlow);
  double binhighcenter = IMnpim_IMnpip_dE_woK0_n_px_1->GetXaxis()->GetBinCenter(binhigh);
  IMnpim_IMnpip_dE_woK0_n_px_1->GetXaxis()->SetRangeUser(anacuts::Sigmap_MIN,anacuts::Sigmap_MAX);
  //IMnpim_IMnpip_dE_woK0_n_px_1->GetXaxis()->SetRangeUser(binlowcenter,binhighcenter);
  IMnpim_IMnpip_dE_woK0_n_px_1->SetFillColor(2);
  IMnpim_IMnpip_dE_woK0_n_px_1->SetFillStyle(3002);
  IMnpim_IMnpip_dE_woK0_n_px_1->Draw("HEsame");
  //TF1 *f1 = new TF1("f1","gaus(0)+pol1(3)",1.16,1.23);
  TF1 *f1 = new TF1("f1","gaus(0)+pol1(3)",anacuts::Sigmap_MIN,anacuts::Sigmap_MAX);
  f1->SetParameter(0,4.34555e+02);
  f1->SetParameter(1,1.18847e+00);
  f1->SetParameter(2, 5.10023e-03);
  f1->SetParameter(3,1.16667e-05);
  f1->SetParameter(4, 3.61475e+02);
  //IMnpim_IMnpip_dE_woK0_n_px->Fit("f1","","",1.16,1.23);
  IMnpim_IMnpip_dE_woK0_n_px->Fit("f1","","",anacuts::Sigmap_MIN,anacuts::Sigmap_MAX);
  TF1 *fbgSp = new TF1("fbgSp","pol1",binlowcenter-0.001,binhighcenter+0.001);
  fbgSp->SetLineColor(4);
  fbgSp->SetFillColor(4);
  fbgSp->SetFillStyle(3002);
  fbgSp->SetParameter(0,f1->GetParameter(3));
  fbgSp->SetParameter(1,f1->GetParameter(4));
  fbgSp->Draw("same");
  double bgsp = fbgSp->Integral(anacuts::Sigmap_MIN,anacuts::Sigmap_MAX);


  TH1D *IMnpim_IMnpip_dE_woK0_n_px_2 = (TH1D*)IMnpim_IMnpip_dE_woK0_n_px->Clone();
  IMnpim_IMnpip_dE_woK0_n_px_2->GetXaxis()->SetRangeUser(anacuts::Sigmap_sidelow_MIN,anacuts::Sigmap_sidelow_MAX);
  IMnpim_IMnpip_dE_woK0_n_px_2->SetFillColor(3);
  IMnpim_IMnpip_dE_woK0_n_px_2->Draw("HEsame");

  TH1D *IMnpim_IMnpip_dE_woK0_n_px_3 = (TH1D*)IMnpim_IMnpip_dE_woK0_n_px->Clone();
  IMnpim_IMnpip_dE_woK0_n_px_3->GetXaxis()->SetRangeUser(anacuts::Sigmap_sidehigh_MIN,anacuts::Sigmap_sidehigh_MAX);
  IMnpim_IMnpip_dE_woK0_n_px_3->SetFillColor(3);
  IMnpim_IMnpip_dE_woK0_n_px_3->Draw("HEsame");
  TH1D *IMnpim_IMnpip_dE_woK0_n_Sp_side_0_px = IMnpim_IMnpip_dE_woK0_n_Sp_side[0]->ProjectionX();
  IMnpim_IMnpip_dE_woK0_n_Sp_side_0_px->SetFillColor(4);
  IMnpim_IMnpip_dE_woK0_n_Sp_side_0_px->SetFillStyle(3009);
  IMnpim_IMnpip_dE_woK0_n_Sp_side_0_px->Draw("HEsame");
  TH1D *IMnpim_IMnpip_dE_woK0_n_Sp_side_1_px = IMnpim_IMnpip_dE_woK0_n_Sp_side[1]->ProjectionX();
  IMnpim_IMnpip_dE_woK0_n_Sp_side_1_px->SetFillColor(4);
  IMnpim_IMnpip_dE_woK0_n_Sp_side_1_px->SetFillStyle(3009);
  IMnpim_IMnpip_dE_woK0_n_Sp_side_1_px->Draw("HEsame");
  
  std::cout << "Sigma+ signal region:    " << IMnpim_IMnpip_dE_woK0_n_px_1->Integral() << std::endl;
  std::cout << "Sigma+ sideband low:     " << IMnpim_IMnpip_dE_woK0_n_px_2->Integral() << std::endl;
  std::cout << "Sigma+ sideband low cut: " << IMnpim_IMnpip_dE_woK0_n_Sp_side_0_px->Integral() << std::endl;
  std::cout << "Sigma+ sideband high:    " << IMnpim_IMnpip_dE_woK0_n_px_3->Integral() << std::endl;
  std::cout << "Sigma+ sideband high cut:" << IMnpim_IMnpip_dE_woK0_n_Sp_side_1_px->Integral() << std::endl;
  std::cout << "bg (Integral)           :" << bgsp       << std::endl;
  std::cout << "bg (Integral)/binw      :" << bgsp/IMnpim_IMnpip_dE_woK0_n_px_1->GetBinWidth(100) << std::endl;
  double trapezoidbgSp = (fbgSp->Eval(anacuts::Sigmap_MIN)+fbgSp->Eval(anacuts::Sigmap_MAX))*
                         (anacuts::Sigmap_MAX-anacuts::Sigmap_MIN)
                         /IMnpim_IMnpip_dE_woK0_n_px_1->GetBinWidth(100)/2.0;
  std::cout << "bg (trapezoid)          :" << trapezoidbgSp << std::endl;



  TCanvas *cIMnpim_IMnpip_dE_woK0_n_py2 = new TCanvas("cIMnpim_IMnpip_dE_woK0_n_py2","IMnpim_IMnpip_dE_woK0_n_py2");
  cIMnpim_IMnpip_dE_woK0_n_py2->cd();
 // TH1D *IMnpim_IMnpip_dE_woK0_n_py = IMnpim_IMnpip_dE_woK0_n->ProjectionY();
  TH1D *IMnpim_IMnpip_dE_woK0_n_py = IMnpim_IMnpip_dE_woK0_n_Sm_bg->ProjectionY();
  IMnpim_IMnpip_dE_woK0_n_py->GetXaxis()->SetRangeUser(1.0,1.3);
  //IMnpim_IMnpip_dE_woK0_n_py->GetYaxis()->SetTitle("Counts/ 10 MeV/c^{2}");
  IMnpim_IMnpip_dE_woK0_n_py->GetYaxis()->CenterTitle();
  IMnpim_IMnpip_dE_woK0_n_py->Draw("EH");
  
  TH1D *IMnpim_IMnpip_dE_woK0_n_py_1 = (TH1D*)IMnpim_IMnpip_dE_woK0_n_py->Clone();
  IMnpim_IMnpip_dE_woK0_n_py_1->GetXaxis()->SetRangeUser(anacuts::Sigmam_MIN,anacuts::Sigmam_MAX);
  IMnpim_IMnpip_dE_woK0_n_py_1->SetFillColor(2);
  IMnpim_IMnpip_dE_woK0_n_py_1->SetFillStyle(3002);
  IMnpim_IMnpip_dE_woK0_n_py_1->Draw("HEsame");
  //TF1 *fbgSm = new TF1("fbgSm","gaus(0)+pol2(3)",1.16,1.23);
  TF1 *fbgSm = new TF1("fbgSm","gaus(0)+pol2(3)",anacuts::Sigmam_MIN,anacuts::Sigmam_MAX);
  fbgSm->SetParameter(0,4.34555e+02);
  fbgSm->SetParameter(1,1.18847e+00);
  fbgSm->SetParameter(2, 5.10023e-03);
  fbgSm->SetParameter(3,1.16667e-05);
  fbgSm->SetParameter(4, 3.61475e+02);
  //IMnpim_IMnpip_dE_woK0_n_py->Fit("fbgSm","","",1.16,1.23);
  IMnpim_IMnpip_dE_woK0_n_py->Fit("fbgSm","","",anacuts::Sigmam_MIN,anacuts::Sigmam_MAX);
  TF1 *fpol2 = new TF1("fpol2","pol2",anacuts::Sigmam_MIN-0.001,anacuts::Sigmam_MAX+0.002);
  fpol2->SetLineColor(4);
  fpol2->SetFillColor(4);
  fpol2->SetFillStyle(3002);
  fpol2->SetParameter(0,fbgSm->GetParameter(3));
  fpol2->SetParameter(1,fbgSm->GetParameter(4));
  fpol2->SetParameter(2,fbgSm->GetParameter(5));
  fpol2->Draw("same");
  double bgsm = fpol2->Integral(anacuts::Sigmam_MIN,anacuts::Sigmam_MAX);
  
  TH1D *IMnpim_IMnpip_dE_woK0_n_py_2 = (TH1D*)IMnpim_IMnpip_dE_woK0_n_py->Clone();
  IMnpim_IMnpip_dE_woK0_n_py_2->GetXaxis()->SetRangeUser(anacuts::Sigmam_sidelow_MIN,anacuts::Sigmam_sidelow_MAX);
  IMnpim_IMnpip_dE_woK0_n_py_2->SetFillColor(3);
  IMnpim_IMnpip_dE_woK0_n_py_2->Draw("HEsame");
  TH1D *IMnpim_IMnpip_dE_woK0_n_py_3 = (TH1D*)IMnpim_IMnpip_dE_woK0_n_py->Clone();
  IMnpim_IMnpip_dE_woK0_n_py_3->GetXaxis()->SetRangeUser(anacuts::Sigmam_sidehigh_MIN,anacuts::Sigmam_sidehigh_MAX);
  IMnpim_IMnpip_dE_woK0_n_py_3->SetFillColor(3);
  IMnpim_IMnpip_dE_woK0_n_py_3->Draw("HEsame");
  
  TH1D *IMnpim_IMnpip_dE_woK0_n_Sm_side_0_py = IMnpim_IMnpip_dE_woK0_n_Sm_side[0]->ProjectionY();
  IMnpim_IMnpip_dE_woK0_n_Sm_side_0_py->SetFillColor(4);
  IMnpim_IMnpip_dE_woK0_n_Sm_side_0_py->SetFillStyle(3009);
  IMnpim_IMnpip_dE_woK0_n_Sm_side_0_py->Draw("HEsame");
  TH1D *IMnpim_IMnpip_dE_woK0_n_Sm_side_1_py = IMnpim_IMnpip_dE_woK0_n_Sm_side[1]->ProjectionY();
  IMnpim_IMnpip_dE_woK0_n_Sm_side_1_py->SetFillColor(4);
  IMnpim_IMnpip_dE_woK0_n_Sm_side_1_py->SetFillStyle(3009);
  IMnpim_IMnpip_dE_woK0_n_Sm_side_1_py->Draw("HEsame");

  std::cout << "Sigma- signal region:    " << IMnpim_IMnpip_dE_woK0_n_py_1->Integral() << std::endl;
  std::cout << "Sigma- sideband low:     " << IMnpim_IMnpip_dE_woK0_n_py_2->Integral() << std::endl;
  std::cout << "Sigma- sideband low cut: " << IMnpim_IMnpip_dE_woK0_n_Sm_side_0_py->Integral() << std::endl;
  std::cout << "Sigma- sideband high:    " << IMnpim_IMnpip_dE_woK0_n_py_3->Integral() << std::endl;
  std::cout << "Sigma- sideband high cut:" << IMnpim_IMnpip_dE_woK0_n_Sm_side_1_py->Integral() << std::endl;
  std::cout << "bg (Integral)           :" << bgsm << std::endl;
  std::cout << "bg (Integral)/binw      :" << bgsm/IMnpim_IMnpip_dE_woK0_n_py_1->GetBinWidth(100)  << std::endl;
  double trapezoidbgSm = (fbgSm->Eval(anacuts::Sigmam_MIN)+fbgSm->Eval(anacuts::Sigmam_MAX))*
                         (anacuts::Sigmam_MAX-anacuts::Sigmam_MIN)
                         /IMnpim_IMnpip_dE_woK0_n_py_1->GetBinWidth(100)/2.0;
  std::cout << "bg (trapezoid)          :" << trapezoidbgSm    << std::endl;
  
  /*
  TCanvas *cq_IMnpipi_woK0_wSid_n_px_side = new TCanvas("cq_IMnpipi_woK0_wSid_n_px_side","q_IMnpipi_woK0_wSid_n_px_side"); 
  TH1D *q_IMnpipi_wSid_n_px1  = q_IMnpipi_wSid_n->ProjectionX();
  //w/o K0
  //q_IMnpipi_woK0_wSid_n_px->GetYaxis()->SetTitle("Counts / 20 MeV/c^{2}");
  q_IMnpipi_woK0_wSid_n_px->GetYaxis()->CenterTitle();
  q_IMnpipi_woK0_wSid_n_px->Draw("EH");
  //including K0
  q_IMnpipi_wSid_n_px1->SetLineColor(5);
  //q_IMnpipi_wSid_n_px1->GetYaxis()->SetTitle("Counts / 20 MeV/c^{2}");
  q_IMnpipi_wSid_n_px1->GetYaxis()->CenterTitle();
  //q_IMnpipi_wSid_n_px1->Draw("HEsame");
  
  //sideband w/o K0
  TH1D *q_IMnpipi_woK0_wSid_n_side_px = q_IMnpipi_woK0_wSid_n_side->ProjectionX(); 
  q_IMnpipi_woK0_wSid_n_side_px->SetLineColor(4);
  //std::cout << "scale "  ;
  //std::cout << 0.8*(float)q_IMnpipi_woK0_wSid_n_px->Integral()/(float) q_IMnpipi_woK0_wSid_n_side_px->Integral() << std::endl;
  //q_IMnpipi_woK0_wSid_n_side_px->Scale(0.8*q_IMnpipi_woK0_wSid_n_px->Integral()/q_IMnpipi_woK0_wSid_n_side_px->Integral());
  //q_IMnpipi_woK0_wSid_n_side_px->Scale(1.23);
  q_IMnpipi_woK0_wSid_n_side_px->Draw("HEsame");

  //sideband including K0
  TH1D *q_IMnpipi_wSid_n_side_px = q_IMnpipi_wSid_n_side->ProjectionX(); 
  q_IMnpipi_wSid_n_side_px->SetLineColor(3);
  //q_IMnpipi_wSid_n_side_px->Scale(1.23);
  //q_IMnpipi_wSid_n_side_px->Draw("HEsame");
  */
  /*
  //TCanvas *csub = new TCanvas("csub","csub");
  //TH1D *q_IMnpipi_woK0_wSid_n_px_sub = (TH1D*) q_IMnpipi_woK0_wSid_n_px->Clone();
  //q_IMnpipi_woK0_wSid_n_px_sub->Add(q_IMnpipi_woK0_wSid_n_side_px,-1);
  //q_IMnpipi_woK0_wSid_n_px_sub->SetMinimum(0);
  //q_IMnpipi_woK0_wSid_n_px_sub->Draw("EH"); 
  //fkp->Draw("same");
  */ 
  
  TCanvas *cq_IMnpipi_woK0_wSid_n_Sp = new TCanvas("cq_IMnpipi_woK0_wSid_n_Sp","q_IMnpipi_woK0_wSid_n_Sp");
  cq_IMnpipi_woK0_wSid_n_Sp->cd();
  q_IMnpipi_woK0_wSid_n_Sp->Draw("colz");

  
  TCanvas *cq_IMnpipi_woK0_wSid_n_Sm = new TCanvas("cq_IMnpipi_woK0_wSid_n_Sm","q_IMnpipi_woK0_wSid_n_Sm");
  cq_IMnpipi_woK0_wSid_n_Sm->cd();
//  q_IMnpipi_woK0_wSid_n_Sm->SetMaximum(q_IMnpipi_woK0_wSid_n_Sp->GetMaximum());
  q_IMnpipi_woK0_wSid_n_Sm->Draw("colz");
  
  TCanvas *cq_IMnpipi_woK0_wSid_n_Sp_side_low = new TCanvas("cq_IMnpipi_woK0_wSid_n_Sp_side_low","q_IMnpipi_woK0_wSid_n_Sp_side_low");
  cq_IMnpipi_woK0_wSid_n_Sp_side_low->cd();
  q_IMnpipi_woK0_wSid_n_Sp_side[sidebandtype][0]->Draw("colz");
  
  TCanvas *cq_IMnpipi_woK0_wSid_n_Sp_side_high = new TCanvas("cq_IMnpipi_woK0_wSid_n_Sp_side_high","q_IMnpipi_woK0_wSid_n_Sp_side_high");
  cq_IMnpipi_woK0_wSid_n_Sp_side_high->cd();
  q_IMnpipi_woK0_wSid_n_Sp_side[sidebandtype][1]->Draw("colz");




  TCanvas *cq_IMnpipi_woK0_wSid_n_Sm_side_low = new TCanvas("cq_IMnpipi_woK0_wSid_n_Sm_side_low","q_IMnpipi_woK0_wSid_n_Sm_side_low");
  cq_IMnpipi_woK0_wSid_n_Sm_side_low->cd();
  q_IMnpipi_woK0_wSid_n_Sm_side[sidebandtype][0]->Draw("colz");

  TCanvas *cq_IMnpipi_woK0_wSid_n_Sm_side_high = new TCanvas("cq_IMnpipi_woK0_wSid_n_Sm_side_high","q_IMnpipi_woK0_wSid_n_Sm_side_high");
  cq_IMnpipi_woK0_wSid_n_Sm_side_high->cd();
  q_IMnpipi_woK0_wSid_n_Sm_side[sidebandtype][1]->Draw("colz");


  TCanvas *cq_IMnpipi_woK0_wSid_n_px_SpSm = new TCanvas("cq_IMnpipi_woK0_wSid_n_px_SpSm","q_IMnpipi_woK0_wSid_n_px_SpSm"); 
  TH1D *q_IMnpipi_woK0_wSid_n_px = q_IMnpipi_woK0_wSid_n->ProjectionX();
  q_IMnpipi_woK0_wSid_n_px->Draw("EH");
  TH1D* q_IMnpipi_woK0_wSid_n_Sp_px = (TH1D*) q_IMnpipi_woK0_wSid_n_Sp->ProjectionX();
  TH1D* q_IMnpipi_woK0_wSid_n_Sm_px = (TH1D*) q_IMnpipi_woK0_wSid_n_Sm->ProjectionX();
  q_IMnpipi_woK0_wSid_n_Sp_px->SetLineColor(2);
  q_IMnpipi_woK0_wSid_n_Sm_px->SetLineColor(3);
  //q_IMnpipi_wSid_n_px1->Draw("HEsame");
  q_IMnpipi_woK0_wSid_n_Sp_px->Draw("HEsame");
  q_IMnpipi_woK0_wSid_n_Sm_px->Draw("HEsame");
  TH1D* q_IMnpipi_woK0_wSid_n_wocross = (TH1D*)q_IMnpipi_woK0_wSid_n_Sp_px->Clone();
  q_IMnpipi_woK0_wSid_n_wocross->Add(q_IMnpipi_woK0_wSid_n_Sm_px,1);
  q_IMnpipi_woK0_wSid_n_wocross->SetLineColor(4);
  //q_IMnpipi_woK0_wSid_n_wocross->Draw("HEsame");
  

  TCanvas *cq_IMnpipi_woK0_wSid_n_px_Sp = new TCanvas("cq_IMnpipi_woK0_wSid_n_px_Sp","q_IMnpipi_woK0_wSid_n_px_Sp"); 
  //q_IMnpipi_woK0_wSid_n_Sp_px->SetMaximum(q_IMnpipi_woK0_wSid_n_Sm_px->GetMaximum());
  if(qvalcutflag==1 && !Spmode && !Smmode) q_IMnpipi_woK0_wSid_n_Sp_px->SetMaximum(160);
  if(qvalcutflag==2 && !Spmode && !Smmode) q_IMnpipi_woK0_wSid_n_Sp_px->SetMaximum(260);
  q_IMnpipi_woK0_wSid_n_Sp_px->Draw("HE");
  //TH1D* q_IMnpipi_woK0_wSid_n_Sp_side_px = (TH1D*) q_IMnpipi_woK0_wSid_n_Sp_side[sidebandtype]->ProjectionX();
  //q_IMnpipi_woK0_wSid_n_Sp_side_px->SetLineColor(5);
  //q_IMnpipi_woK0_wSid_n_Sp_side_px->Draw("HEsame");
  
  TH1D* q_IMnpipi_woK0_wSid_n_Sp_side_0_px = (TH1D*) q_IMnpipi_woK0_wSid_n_Sp_side[sidebandtype][0]->ProjectionX();
  q_IMnpipi_woK0_wSid_n_Sp_side_0_px->SetLineColor(kCyan);
  q_IMnpipi_woK0_wSid_n_Sp_side_0_px->SetMarkerStyle(20);
  q_IMnpipi_woK0_wSid_n_Sp_side_0_px->SetMarkerColor(kCyan);
  q_IMnpipi_woK0_wSid_n_Sp_side_0_px->Draw("Esame");
  
  TH1D* q_IMnpipi_woK0_wSid_n_Sp_side_1_px = (TH1D*) q_IMnpipi_woK0_wSid_n_Sp_side[sidebandtype][1]->ProjectionX();
  q_IMnpipi_woK0_wSid_n_Sp_side_1_px->SetLineColor(kCyan+2);
  q_IMnpipi_woK0_wSid_n_Sp_side_1_px->SetMarkerStyle(21);
  //q_IMnpipi_woK0_wSid_n_Sp_side_1_px->SetMarkerColor(kOrange+6);
  q_IMnpipi_woK0_wSid_n_Sp_side_1_px->Draw("Esame");
  
  TH1D* q_IMnpipi_woK0_wSid_nSp_side_px_sum = (TH1D*) q_IMnpipi_woK0_wSid_n_Sp_side_0_px->Clone();
  q_IMnpipi_woK0_wSid_nSp_side_px_sum->Add(q_IMnpipi_woK0_wSid_n_Sp_side_1_px);
  q_IMnpipi_woK0_wSid_nSp_side_px_sum->SetLineColor(kCyan+4);
  q_IMnpipi_woK0_wSid_nSp_side_px_sum->SetMarkerStyle(22);
  q_IMnpipi_woK0_wSid_nSp_side_px_sum->SetMarkerColor(kCyan+4);
  q_IMnpipi_woK0_wSid_nSp_side_px_sum->Draw("Esame");

  //TH1D* q_IMnpipi_woK0_wSid_n_Sp_side_px_1 = (TH1D*) q_IMnpipi_woK0_wSid_n_Sp_side[1]->ProjectionX();
  //q_IMnpipi_woK0_wSid_n_Sp_side_px_1->SetLineColor(kOrange+2);
  //q_IMnpipi_woK0_wSid_n_Sp_side_px_1->Draw("HEsame");
  //TH1D* q_IMnpipi_woK0_wSid_n_Sp_side_px_2 = (TH1D*) q_IMnpipi_woK0_wSid_n_Sp_side[2]->ProjectionX();
  //q_IMnpipi_woK0_wSid_n_Sp_side_px_2->SetLineColor(kOrange+4);
  //q_IMnpipi_woK0_wSid_n_Sp_side_px_2->Draw("HEsame");

  

  //sideband subtracted spectra for Sp mode
  TCanvas *cq_IMnpipi_woK0_wSid_n_Sp_sub = new TCanvas("cq_IMnpipi_woK0_wSid_n_Sp_sub","q_IMnpipi_woK0_wSid_n_Sp_sub");
  cq_IMnpipi_woK0_wSid_n_Sp_sub->cd();
  TH2F* q_IMnpipi_woK0_wSid_n_Sp_sub = (TH2F*)q_IMnpipi_woK0_wSid_n_Sp->Clone("Sp_sub");
  q_IMnpipi_woK0_wSid_n_Sp_sub->Add(q_IMnpipi_woK0_wSid_n_Sp_side[sidebandtype][0],-1);
  q_IMnpipi_woK0_wSid_n_Sp_sub->Add(q_IMnpipi_woK0_wSid_n_Sp_side[sidebandtype][1],-1);
  q_IMnpipi_woK0_wSid_n_Sp_sub->SetMinimum(0);
  q_IMnpipi_woK0_wSid_n_Sp_sub->Draw("colz");
  
  
  //overlay signal and sideband projected 1d histo to IM
  TCanvas *cq_IMnpipi_woK0_wSid_n_Sm_px = new TCanvas("cq_IMnpipi_woK0_wSid_n_Sm_px","q_IMnpipi_woK0_wSid_n_Sm_px"); 
  //q_IMnpipi_woK0_wSid_n_Sm_px->SetMaximum(q_IMnpipi_woK0_wSid_n_px->GetMaximum());
  if(qvalcutflag==1) q_IMnpipi_woK0_wSid_n_Sm_px->SetMaximum(160);
  if(qvalcutflag==2) q_IMnpipi_woK0_wSid_n_Sm_px->SetMaximum(260);
  q_IMnpipi_woK0_wSid_n_Sm_px->Draw("EH");
  TH1D* q_IMnpipi_woK0_wSid_n_Sm_side_0_px = (TH1D*) q_IMnpipi_woK0_wSid_n_Sm_side[sidebandtype][0]->ProjectionX();
  q_IMnpipi_woK0_wSid_n_Sm_side_0_px->SetLineColor(kPink+1);
  q_IMnpipi_woK0_wSid_n_Sm_side_0_px->SetMarkerStyle(20);
  q_IMnpipi_woK0_wSid_n_Sm_side_0_px->Draw("Esame");
  
  TH1D* q_IMnpipi_woK0_wSid_n_Sm_side_1_px = (TH1D*) q_IMnpipi_woK0_wSid_n_Sm_side[sidebandtype][1]->ProjectionX();
  q_IMnpipi_woK0_wSid_n_Sm_side_1_px->SetLineColor(kPink+3);
  q_IMnpipi_woK0_wSid_n_Sm_side_1_px->SetMarkerStyle(21);
  q_IMnpipi_woK0_wSid_n_Sm_side_1_px->SetMarkerColor(kPink+3);
  q_IMnpipi_woK0_wSid_n_Sm_side_1_px->Draw("Esame");
  
  TH1D* q_IMnpipi_woK0_wSid_n_Sm_side_px_sum = (TH1D*) q_IMnpipi_woK0_wSid_n_Sm_side_0_px->Clone();
  q_IMnpipi_woK0_wSid_n_Sm_side_px_sum->Add(q_IMnpipi_woK0_wSid_n_Sm_side_1_px);
  q_IMnpipi_woK0_wSid_n_Sm_side_px_sum->SetLineColor(kPink+10);
  q_IMnpipi_woK0_wSid_n_Sm_side_px_sum->SetMarkerStyle(22);
  q_IMnpipi_woK0_wSid_n_Sm_side_px_sum->SetMarkerColor(kPink+10);
  q_IMnpipi_woK0_wSid_n_Sm_side_px_sum->Draw("Esame");

  //TH1D* q_IMnpipi_woK0_wSid_n_Sm_side_px_1 = (TH1D*) q_IMnpipi_woK0_wSid_n_Sm_side[1]->ProjectionX();
  //q_IMnpipi_woK0_wSid_n_Sm_side_px_1->SetLineColor(kPink+3);
  //q_IMnpipi_woK0_wSid_n_Sm_side_px_1->Draw("HEsame");
  //TH1D* q_IMnpipi_woK0_wSid_n_Sm_side_px_2 = (TH1D*) q_IMnpipi_woK0_wSid_n_Sm_side[2]->ProjectionX();
  //q_IMnpipi_woK0_wSid_n_Sm_side_px_2->SetLineColor(kPink+5);
  //q_IMnpipi_woK0_wSid_n_Sm_side_px_2->Draw("HEsame");

  //sideband subtracted spectra for Sm mode
  TCanvas *cq_IMnpipi_woK0_wSid_n_Sm_sub = new TCanvas("cq_IMnpipi_woK0_wSid_n_Sm_sub","q_IMnpipi_woK0_wSid_n_Sm_sub");
  cq_IMnpipi_woK0_wSid_n_Sm_sub->cd();
  TH2F* q_IMnpipi_woK0_wSid_n_Sm_sub = (TH2F*)q_IMnpipi_woK0_wSid_n_Sm->Clone("Sm_sub");
  q_IMnpipi_woK0_wSid_n_Sm_sub->Add(q_IMnpipi_woK0_wSid_n_Sm_side[sidebandtype][0],-1);
  q_IMnpipi_woK0_wSid_n_Sm_sub->Add(q_IMnpipi_woK0_wSid_n_Sm_side[sidebandtype][1],-1);
  q_IMnpipi_woK0_wSid_n_Sm_sub->SetMinimum(0);
  q_IMnpipi_woK0_wSid_n_Sm_sub->Draw("colz");
  
  
  TCanvas *cq_IMnpipi_woK0_wSid_n_SpSm_sub_px = new TCanvas("cq_IMnpipi_woK0_wSid_n_SpSm_sub_px","q_IMnpipi_woK0_wSid_n_SpSm_sub_px");
  cq_IMnpipi_woK0_wSid_n_SpSm_sub_px->cd();
  TH1D* q_IMnpipi_woK0_wSid_n_Sp_sub_px = q_IMnpipi_woK0_wSid_n_Sp_sub->ProjectionX();
  q_IMnpipi_woK0_wSid_n_Sp_sub_px->SetLineColor(2);
  q_IMnpipi_woK0_wSid_n_Sp_sub_px->SetMinimum(0);
  q_IMnpipi_woK0_wSid_n_Sp_sub_px->Draw("HE");

  TH1D* q_IMnpipi_woK0_wSid_n_Sm_sub_px = q_IMnpipi_woK0_wSid_n_Sm_sub->ProjectionX();
  q_IMnpipi_woK0_wSid_n_Sm_sub_px->SetLineColor(3);
  q_IMnpipi_woK0_wSid_n_Sm_sub_px->Draw("HEsame");

  TH1D* SpSm_sub_sum_px = (TH1D*)q_IMnpipi_woK0_wSid_n_Sp_sub_px->Clone("SpSm_sub_px");
  SpSm_sub_sum_px->Add(q_IMnpipi_woK0_wSid_n_Sm_sub_px);
  SpSm_sub_sum_px->SetLineColor(4);
  SpSm_sub_sum_px->Draw("HEsame");


  //w/ K0
  TCanvas *cq_IMnpipi_wSid_n = new TCanvas("cq_IMnpipi_wSid_n","q_IMnpipi_wSid_n");
  cq_IMnpipi_wSid_n->cd();
  q_IMnpipi_wSid_n->Draw("colz");
  
  TCanvas *cq_IMnpipi_wSid_n_side = new TCanvas("cq_IMnpipi_wSid_n_side","q_IMnpipi_wSid_n_side");
  cq_IMnpipi_wSid_n_side->cd();
  q_IMnpipi_wSid_n_side->SetMaximum(q_IMnpipi_wSid_n->GetMaximum());
  q_IMnpipi_wSid_n_side->Draw("colz");


  TCanvas *cq_IMnpipi_woK0_wSid_n = new TCanvas("cq_IMnpipi_woK0_wSid_n","q_IMnpipi_woK0_wSid_n");
  cq_IMnpipi_woK0_wSid_n->cd();
  q_IMnpipi_woK0_wSid_n->Draw("colz");
  const double Kp_mass = pMass + kpMass;  
  TF1 *fkp = new TF1("f", "sqrt(((x*x-[0]*[0]-[1]*[1])/(2*[0]))*((x*x-[0]*[0]-[1]*[1])/(2*[0]))-[1]*[1])",Kp_mass-0.001,2);
  fkp->SetParameter(0,pMass);
  fkp->SetParameter(1,kpMass);
  //fkp->SetLineColor(4);
  fkp->SetLineWidth(4);
  fkp->SetLineStyle(4);
  fkp->SetLineColorAlpha(kPink, 0.35);
  fkp->Draw("same");
  //q_IMnpipi_woK0_wSid_n->Draw("colzsame");
  //fnu->cd();
  //TMultiGraph *mg = (TMultiGraph*)fnu->Get("mg"); 
  //mg->Draw("c");
  f->cd();
  //TCanvas *ctest = new TCanvas("ctest","ctezst");
  //TLatex *tex = new TLatex();
  //tex->SetTextSize(0.04);
  //tex->SetTextAngle(60);
  //tex->DrawLatexNDC( 0.90, 0.90, "M_{F}(q)" );
  //tex->SetTextAngle(0);

  //TCanvas *cq_IMnpipi_woK0_wSid_n_side = new TCanvas("cq_IMnpipi_woK0_wSid_n_side","q_IMnpipi_woK0_wSid_n_side");
  //cq_IMnpipi_woK0_wSid_n_side->cd();
  //cq_IMnpipi_woK0_wSid_n->cd(2);
  //q_IMnpipi_woK0_wSid_n_side->SetMaximum(q_IMnpipi_woK0_wSid_n->GetMaximum());
  //q_IMnpipi_woK0_wSid_n_side->Draw("colz");
  


  //TCanvas *cqsub2 = new TCanvas("cqsub2","cqsub2");
  //cqsub2->cd();
  //TH2F *q_IMnpipi_woK0_wSid_n_sub = (TH2F*) q_IMnpipi_woK0_wSid_n->Clone("sub");
  //q_IMnpipi_woK0_wSid_n_sub->Add(q_IMnpipi_woK0_wSid_n_side,-1);
  //q_IMnpipi_woK0_wSid_n_sub->SetMinimum(0);
  //q_IMnpipi_woK0_wSid_n_sub->Draw("colz");
  //fkp->Draw("same");
  //TFile *fnu = new TFile("NumericalRootFinder.root");
  //fnu->cd();
  //TMultiGraph *mg = (TMultiGraph*)fnu->Get("mg"); 
  //mg->Draw("c");
  //f->cd();
  //fnu->cd();
  //TMultiGraph *mg = (TMultiGraph*)fnu->Get("mg"); 
  //mg->Draw("c");
  //f->cd();

  TCanvas *cdE_betainv_fid = new TCanvas(Form("cdE_betainv_fid"), "dE_betainv_fid");
  cdE_betainv_fid->cd();
  dE_betainv_fid->GetXaxis()->SetRangeUser(0,10);
  dE_betainv_fid->Draw("colz");


  TCanvas *cdE_betainv_fid_px = new TCanvas("cdE_betainv_fid_px","dE_betainv_fid_px");
  int bin2mev = dE_betainv_fid->GetYaxis()->FindBin(2.0);
  int bin4mev = dE_betainv_fid->GetYaxis()->FindBin(4.0);
  int bin6mev = dE_betainv_fid->GetYaxis()->FindBin(6.0);

  TH1D* h1_nocut = dE_betainv_fid->ProjectionX("px");
  TH1D* h1_2mevcut = dE_betainv_fid->ProjectionX("px1",bin2mev,-1);
  TH1D* h1_4mevcut = dE_betainv_fid->ProjectionX("px2",bin4mev,-1);
  TH1D* h1_6mevcut = dE_betainv_fid->ProjectionX("px3",bin6mev,-1);

  h1_nocut->Draw("HE");
  h1_2mevcut->SetLineColor(2);
  h1_2mevcut->Draw("HEsame");
  h1_4mevcut->SetLineColor(3);
  h1_4mevcut->Draw("HEsame");
  h1_6mevcut->SetLineColor(4);
  h1_6mevcut->Draw("HEsame");

  //TCanvas *cIMnpim_IMnpip_dE_woK0 = new TCanvas(Form("cIMnpim_IMnpip_dE_woK0"),"IMnpim_IMnpip_dE_woK0");
  //cIMnpim_IMnpip_dE_woK0->cd();
  //IMnpim_IMnpip_dE_woK0->Draw("colz");

  /*
     TCanvas *cnmom_IMnpipi_woK0_wSid_n = new TCanvas(Form("cnmom_IMnpipi_woK0_wSid_n"),"nmom_IMnpipi_woK0_wSid_n");
     cnmom_IMnpipi_woK0_wSid_n->cd();
     nmom_IMnpipi_woK0_wSid_n->Draw("colz");

     TCanvas *cnmom_IMnpipi_woK0_wSid_n_py = new TCanvas(Form("cnmom_IMnpipi_woK0_wSid_n_py"),"nmom_IMnpipi_woK0_wSid_n_py");
     cnmom_IMnpipi_woK0_wSid_n_py->cd();
     TH1D* nmom_IMnpipi_woK0_wSid_n_py = nmom_IMnpipi_woK0_wSid_n->ProjectionY();
     nmom_IMnpipi_woK0_wSid_n_py->Draw();
     */

  //TCanvas *cMMnmiss_IMnpipi_woK0_wSid_n_px = new TCanvas("cMMnmiss_IMnpipi_woK0_wSid_n_px","MMnmiss_IMnpipi_woK0_wSid_n_px"); 
  //cMMnmiss_IMnpipi_woK0_wSid_n_px->cd();
  //MMnmiss_IMnpipi_woK0_wSid_n->ProjectionX()->Draw("");
  //TCanvas *cMMnmiss_IMnpipi_woK0_wSid_n_py = new TCanvas("cMMnmiss_IMnpipi_woK0_wSid_n_py","MMnmiss_IMnpipi_woK0_wSid_n_py"); 
  //cMMnmiss_IMnpipi_woK0_wSid_n_py->cd();
  //MMnmiss_IMnpipi_woK0_wSid_n->ProjectionY()->Draw("");

  TCanvas *cMMnmiss_IMnpipi_woK0_wSid_n = new TCanvas("cMMnmiss_IMnpipi_woK0_wSid_n","MMnmiss_IMnpipi_woK0_wSid_n");
  cMMnmiss_IMnpipi_woK0_wSid_n->cd();
  MMnmiss_IMnpipi_woK0_wSid_n->Draw("colz");

  //TCanvas *cq_pippim_n = new TCanvas("cq_pippim_n","q_pippim_n");
  //cq_pippim_n->cd();
  //q_pippim_n->Draw("colz");

  TCanvas *cq_pippim_n_px = new TCanvas("cq_pippim_n_px","q_pippim_n_px");
  cq_pippim_n_px->cd();
  TH1D* q_pippim_n_px = q_pippim_n->ProjectionX();
  q_pippim_n_px->GetXaxis()->SetRangeUser(0.4,0.6);
  q_pippim_n_px->GetYaxis()->SetTitle("Counts / 2 MeV/c^{2}");
  q_pippim_n_px->GetYaxis()->CenterTitle();
  q_pippim_n_px->Draw("EH");
  TH1D* q_pippim_n_px_cut = (TH1D*)q_pippim_n_px->Clone();
  q_pippim_n_px_cut->GetXaxis()->SetRangeUser(anacuts::pipi_MIN,anacuts::pipi_MAX);
  q_pippim_n_px_cut->SetFillColor(2);
  q_pippim_n_px_cut->Draw("HEsame");

  //TCanvas *cMMnpip_MMnpim_woK0_wSid_n = new TCanvas("cMMnpip_MMnpim_woK0_wSid_n","MMnpip_MMnpim_woK0_wSid_n");
  //cMMnpip_MMnpim_woK0_wSid_n->cd();
  //MMnpip_MMnpim_woK0_wSid_n->Draw("colz");
  
  //TCanvas *cMMom_MMass_woK0  = new TCanvas("cMMom_MMass_woK0","MMom_MMass_woK0");
  //cMMom_MMass_woK0->cd();
  //MMom_MMass_woK0->Draw("colz");
  
  //TCanvas *cMMom_MMass_woK0_px = new TCanvas("cMMom_MMass_woK0_px","MMom_MMass_woK0_px");
  //TH1D *MMom_MMass_woK0_px = MMom_MMass_woK0->ProjectionX();
  //MMom_MMass_woK0_px->Draw("");
  
  TCanvas *cMMom_MMass_woK0_wSid  = new TCanvas("cMMom_MMass_woK0_wSid","MMom_MMass_woK0_wSid");
  cMMom_MMass_woK0_wSid->cd();
  MMom_MMass_woK0_wSid->Draw("colz");
  

  TCanvas *cMMom_MMass_woK0_wSid_px  = new TCanvas("cMMom_MMass_woK0_wSid_px","MMom_MMass_woK0_wSid_px");
  cMMom_MMass_woK0_wSid_px->cd();
  TH1D *MMom_MMass_woK0_wSid_px = MMom_MMass_woK0_wSid->ProjectionX();
  MMom_MMass_woK0_wSid_px->GetYaxis()->SetTitle("Counts/ 10 MeV/c^{2}");
  MMom_MMass_woK0_wSid_px->GetYaxis()->CenterTitle();
  MMom_MMass_woK0_wSid_px->Draw("HE");
  TH1D *MMom_MMass_woK0_wSid_px_clone1= (TH1D*) MMom_MMass_woK0_wSid_px->Clone();
  MMom_MMass_woK0_wSid_px_clone1->GetXaxis()->SetRangeUser(anacuts::neutron_MIN,anacuts::neutron_MAX);
  MMom_MMass_woK0_wSid_px_clone1->SetFillColor(4);
  MMom_MMass_woK0_wSid_px_clone1->Draw("HEsame");
  


  //TCanvas *cMMom_MMass_woK0_px_sup = new TCanvas("cMMom_MMass_woK0_px_sup","MMom_MMass_woK0_px_sup");
  //MMom_MMass_woK0_px->Draw();
  TH1D *MMom_MMass_woK0_wSid_px_clone = (TH1D*)MMom_MMass_woK0_wSid_px->Clone(); 
  //MMom_MMass_woK0_wSid_px_clone->SetLineColor(4);
  //cMMom_MMass_woK0_px_sup->cd();
  MMom_MMass_woK0_wSid_px_clone->SetTitle("Missing Mass d(K^{-},#pi^{+}#pi^{-}n)\"X\"");
  //MMom_MMass_woK0_wSid_px_clone->Draw("");

  //TCanvas *cMMom_MMass_woK0_py = new TCanvas("cMMom_MMass_woK0_py","MMom_MMass_woK0_py");
  TH1D *MMom_MMass_woK0_py = MMom_MMass_woK0->ProjectionY();
  TH1D *MMom_MMass_woK0_wSid_py = MMom_MMass_woK0_wSid->ProjectionY();
  //MMom_MMass_woK0_py->Draw();
  TH1D *MMom_MMass_woK0_wSid_py_clone = (TH1D*)MMom_MMass_woK0_wSid_py->Clone();
  //MMom_MMass_woK0_wSid_py_clone->SetLineColor(4);
  //cMMom_MMass_woK0_py->cd();
  MMom_MMass_woK0_wSid_py_clone->SetTitle("Missing Mom. d(K^{-},#pi^{+}#pi^{-}n)\"X\"");
  //MMom_MMass_woK0_wSid_py_clone->Draw("");
  
  /*
  TCanvas *cdE_IMnpim_woK0_n = new TCanvas("cdE_IMnpim_woK0_n","dE_IMnpim_woK0_n");
  cdE_IMnpim_woK0_n->cd();
  dE_IMnpim_woK0_n->Draw("colz");

  TCanvas *cdE_IMnpip_woK0_n = new TCanvas("cdE_IMnpip_woK0_n","dE_IMnpip_woK0_n");
  cdE_IMnpip_woK0_n->cd();
  dE_IMnpip_woK0_n->Draw("colz");
  */
  
  //TCanvas *cIMnpim_IMnpip_dE_woK0_px = new TCanvas("cIMnpim_IMnpip_dE_woK0_px","IMnpim_IMnpip_dE_woK0_px");
  //cIMnpim_IMnpip_dE_woK0_px->cd();
  //TH1D *IMnpim_IMnpip_dE_woK0_px = IMnpim_IMnpip_dE_woK0->ProjectionX();
  //IMnpim_IMnpip_dE_woK0_px->Draw("EH");
  
  //TCanvas *cIMnpim_IMnpip_dE_woK0_py = new TCanvas("cIMnpim_IMnpip_dE_woK0_py","IMnpim_IMnpip_dE_woK0_py");
  //cIMnpim_IMnpip_dE_woK0_py->cd();
  //TH1D *IMnpim_IMnpip_dE_woK0_py = IMnpim_IMnpip_dE_woK0->ProjectionY();
  //IMnpim_IMnpip_dE_woK0_py->Draw("HE");

  /*
  TCanvas *c_DCA = new TCanvas("c_DCA","c_DCA");
  c_DCA->cd();
  DCA_pip_beam->Draw("");
  std::cout << DCA_pip_beam->Integral() << std::endl;
  DCA_pim_beam->SetLineColor(2);
  DCA_pim_beam->Draw("same");
  std::cout << DCA_pim_beam->Integral() << std::endl;
  DCA_pip_pim->SetLineColor(3);
  DCA_pip_pim->Draw("same");
  std::cout << DCA_pip_pim->Integral() << std::endl;
  
  
  TCanvas *c_DCA_Sp = new TCanvas("c_DCA_Sp","c_DCA_Sp");
  c_DCA_Sp->cd();
  
  TCanvas *c_DCA_Sm = new TCanvas("c_DCA_Sm","c_DCA_Sm");
  c_DCA_Sm->cd();
   */
  
  //acceptance calculation
  TFile *facc = new TFile("acc.root","READ");
  
  if(Spmode || Smmode){
    if(Spmode){
      std::cout << "This is Sigma+ mode sim." << std::endl;
    }else{
      std::cout << "This is Sigma- mode sim." << std::endl;
    }
    TString sacc = std::string(filename);
    sacc.Replace(sacc.Length()-5,10,"_acc.root");
    std::cout << "acc file Sp mode: " << sacc.Data() << std::endl;
    TFile *fsacc = new TFile(sacc.Data(),"RECREATE");

    TCanvas *cphase = new TCanvas("cphase","cphase");
    cphase->cd();
    q_IMnpipi_gen->Draw("colz");
    TCanvas *ceff = new TCanvas("ceff","ceff");
    TH2F *h2acc=NULL;
    if(Spmode) h2acc =  (TH2F*)q_IMnpipi_woK0_wSid_n_Sp_acc->Clone();
    else       h2acc =  (TH2F*)q_IMnpipi_woK0_wSid_n_Sm_acc->Clone();
    h2acc->Sumw2();
    q_IMnpipi_woK0_wSid_n_Sp->Sumw2();
    q_IMnpipi_woK0_wSid_n_Sm->Sumw2();
    q_IMnpipi_gen->Sumw2();
    q_IMnpipi_gen->Print("base");
    q_IMnpipi_woK0_wSid_n_Sp_acc->Print("base");
    std::cout << "calc. acc." << std::endl;
    TEfficiency *pEff ;

    //cleaning 
    /*
    for(int ibinx=0;ibinx<q_IMnpipi_gen->GetNbinsX();ibinx++){
      for(int ibiny=0;ibiny<q_IMnpipi_gen->GetNbinsY();ibiny++){
        int bingen =  q_IMnpipi_gen->GetBinContent(ibinx,ibiny);
        int binacc =  q_IMnpipi_woK0_wSid_n_Sp_acc_reco->GetBinContent(ibinx,ibiny);
        if(binacc>=bingen) q_IMnpipi_woK0_wSid_n_Sp_acc_reco->SetBinContent(ibinx,ibiny,0.0);
        binacc =  q_IMnpipi_woK0_wSid_n_Sm_acc->GetBinContent(ibinx,ibiny);
        if(binacc>=bingen) q_IMnpipi_woK0_wSid_n_Sm_acc->SetBinContent(ibinx,ibiny,0.0);
        binacc =  q_IMnpipi_wSid_n_acc->GetBinContent(ibinx,ibiny);
        if(binacc>=bingen) q_IMnpipi_wSid_n_acc->SetBinContent(ibinx,ibiny,0.0);
        binacc =  q_IMnpipi_woK0_wSid_n_acc->GetBinContent(ibinx,ibiny);
        if(binacc>=bingen) q_IMnpipi_woK0_wSid_n_acc->SetBinContent(ibinx,ibiny,0.0);
        
        if(bingen<100){
          q_IMnpipi_wSid_n_acc->SetBinContent(ibinx,ibiny,0.0); 
          q_IMnpipi_woK0_wSid_n_acc->SetBinContent(ibinx,ibiny,0.0); 
          q_IMnpipi_woK0_wSid_n_Sp_acc->SetBinContent(ibinx,ibiny,0.0); 
          q_IMnpipi_woK0_wSid_n_Sm_acc->SetBinContent(ibinx,ibiny,0.0); 
        }
      }
    }*/


    q_IMnpipi_gen->Write();
    q_IMnpipi_woK0_wSid_n_Sp_acc->Write();
    q_IMnpipi_woK0_wSid_n_Sp_acc_reco->Write();
    q_IMnpipi_woK0_wSid_n_Sm_acc->Write();
    q_IMnpipi_woK0_wSid_n_Sm_acc_reco->Write();
    q_IMnpipi_wSid_n_acc->Write();
    q_IMnpipi_wSid_n_acc_reco->Write();
    q_IMnpipi_woK0_wSid_n_acc->Write();
    q_IMnpipi_woK0_wSid_n_acc_reco->Write();
    
    fsacc->Close();
  }//Spmode or Smmode
  
  /*
  facc->cd();
  TH2F* acc_Sp_cal = (TH2F*)facc->Get("acc_Sp");
  if(acc_Sp_cal == NULL){
    std::cout << " acc_Sp is NULL " << std::endl;
  }
  TH2F* acc_Sm_cal = (TH2F*)facc->Get("acc_Sm");
  if(acc_Sm_cal == NULL){
    std::cout << " acc_Sm is NULL " << std::endl;
  }
  f->cd(); 
  std::cout << std::endl;
  std::cout << "calculation CS of Sp mode..." << std::endl;
  std::cout << std::endl;
  TCanvas *cq_IMnpipi_woK0_wSid_n_Sp_cs = new TCanvas("cq_IMnpipi_woK0_wSid_n_Sp_cs","q_IMnpipi_woK0_wSid_n_Sp_cs");
  cq_IMnpipi_woK0_wSid_n_Sp_cs->cd();
  TH2F *q_IMnpipi_woK0_wSid_n_Sp_cs = (TH2F*)q_IMnpipi_woK0_wSid_n_Sp->Clone("Sp_cs");
  q_IMnpipi_woK0_wSid_n_Sp_cs->Sumw2();
  acc_Sp_cal->Print("base");
  //acc_Sp_cal->RebinX(5);
  q_IMnpipi_woK0_wSid_n_Sp_cs->Print("base");
  q_IMnpipi_woK0_wSid_n_Sp_cs->Divide(acc_Sp_cal);
  q_IMnpipi_woK0_wSid_n_Sp_cs->SetMaximum(30000);
  q_IMnpipi_woK0_wSid_n_Sp_cs->Draw("colz");
  
  TCanvas *cq_IMnpipi_woK0_wSid_n_Sp_cs_px = new TCanvas("cq_IMnpipi_woK0_wSid_n_Sp_cs_px","q_IMnpipi_woK0_wSid_n_Sp_cs_px");
  TH1D *q_IMnpipi_woK0_wSid_n_Sp_cs_px = (TH1D*)q_IMnpipi_woK0_wSid_n_Sp_cs->ProjectionX();
  q_IMnpipi_woK0_wSid_n_Sp_cs_px->SetLineColor(2);
  q_IMnpipi_woK0_wSid_n_Sp_cs_px->Draw("HE");
  
  TCanvas *cq_IMnpipi_woK0_wSid_n_Sm_cs = new TCanvas("cq_IMnpipi_woK0_wSid_n_Sm_cs","q_IMnpipi_woK0_wSid_n_Sm_cs");
  cq_IMnpipi_woK0_wSid_n_Sm_cs->cd();
  TH2F *q_IMnpipi_woK0_wSid_n_Sm_cs = (TH2F*)q_IMnpipi_woK0_wSid_n_Sm->Clone("Sm_cs");
  q_IMnpipi_woK0_wSid_n_Sm_cs->Sumw2();
  //acc_Sm_cal->Print("base");
  //  q_IMnpipi_woK0_wSid_n_Sm_cs->Print("base");
  q_IMnpipi_woK0_wSid_n_Sm_cs->Divide(acc_Sm_cal);
  q_IMnpipi_woK0_wSid_n_Sm_cs->SetMaximum(30000);
  q_IMnpipi_woK0_wSid_n_Sm_cs->Draw("colz");

  //TCanvas *cq_IMnpipi_woK0_wSid_n_Sm_cs_px = new TCanvas("cq_IMnpipi_woK0_wSid_n_Sm_cs_px","q_IMnpipi_woK0_wSid_n_Sm_cs_px");
  TH1D *q_IMnpipi_woK0_wSid_n_Sm_cs_px = (TH1D*)q_IMnpipi_woK0_wSid_n_Sm_cs->ProjectionX();
  q_IMnpipi_woK0_wSid_n_Sm_cs_px->SetLineColor(3);
  q_IMnpipi_woK0_wSid_n_Sm_cs_px->Draw("HEsame");

  TH1D* cs_sum = (TH1D*)q_IMnpipi_woK0_wSid_n_Sp_cs_px->Clone();
  cs_sum->Add(q_IMnpipi_woK0_wSid_n_Sm_cs_px);
  cs_sum->SetLineColor(4);
  cs_sum->Draw("HEsame");

  TCanvas *cq_IMnpipi_woK0_wSid_n_Sp_sub_cs = new TCanvas("cq_IMnpipi_woK0_wSid_n_Sp_sub_cs","q_IMnpipi_woK0_wSid_n_Sp_sub_cs");
  cq_IMnpipi_woK0_wSid_n_Sp_sub_cs->cd();
  TH2F *q_IMnpipi_woK0_wSid_n_Sp_sub_cs = (TH2F*)q_IMnpipi_woK0_wSid_n_Sp_sub->Clone("Sp_cs_sub");
  q_IMnpipi_woK0_wSid_n_Sp_sub_cs->Sumw2();
  q_IMnpipi_woK0_wSid_n_Sp_sub_cs->Divide(acc_Sp_cal);
  q_IMnpipi_woK0_wSid_n_Sp_sub_cs->SetMaximum(30000);
  q_IMnpipi_woK0_wSid_n_Sp_sub_cs->Draw("colz");
  

  TCanvas *cq_IMnpipi_woK0_wSid_n_Sm_sub_cs = new TCanvas("cq_IMnpipi_woK0_wSid_n_Sm_sub_cs","q_IMnpipi_woK0_wSid_n_Sm_sub_cs");
  cq_IMnpipi_woK0_wSid_n_Sm_sub_cs->cd();
  TH2F *q_IMnpipi_woK0_wSid_n_Sm_sub_cs = (TH2F*)q_IMnpipi_woK0_wSid_n_Sm_sub->Clone("Sm_cs_sub");
  q_IMnpipi_woK0_wSid_n_Sm_sub_cs->Sumw2();
  q_IMnpipi_woK0_wSid_n_Sm_sub_cs->Divide(acc_Sm_cal);
  q_IMnpipi_woK0_wSid_n_Sm_sub_cs->SetMaximum(30000);
  q_IMnpipi_woK0_wSid_n_Sm_sub_cs->Draw("colz");

  TCanvas *cq_IMnpipi_woK0_wSid_n_SpSm_sub_cs = new TCanvas("cq_IMnpipi_woK0_wSid_n_SpSm_sub_cs","q_IMnpipi_woK0_wSid_n_SpSm_sub_cs");
  cq_IMnpipi_woK0_wSid_n_SpSm_sub_cs->cd();
  TH2F* q_IMnpipi_woK0_wSid_n_SpSm_sub_cs = (TH2F*)q_IMnpipi_woK0_wSid_n_Sp_sub_cs->Clone("SpSm_cs_sub");    
  q_IMnpipi_woK0_wSid_n_SpSm_sub_cs->Add(q_IMnpipi_woK0_wSid_n_Sm_sub_cs);
  q_IMnpipi_woK0_wSid_n_SpSm_sub_cs->SetMinimum(0);
  q_IMnpipi_woK0_wSid_n_SpSm_sub_cs->SetMaximum(30000);
  q_IMnpipi_woK0_wSid_n_SpSm_sub_cs->SetTitle("q_IMnpipi_woK0_wSid_n_SpSm_sub_cs"); 
  q_IMnpipi_woK0_wSid_n_SpSm_sub_cs->Draw("colz");
  
  gStyle->SetErrorX(0) ;
  TCanvas *cq_IMnpipi_woK0_wSid_n_Sp_sub_cs_px = new TCanvas("cq_IMnpipi_woK0_wSid_n_Sp_sub_cs_px","q_IMnpipi_woK0_wSid_n_Sp_sub_cs_px");
  TH1D *q_IMnpipi_woK0_wSid_n_Sp_sub_cs_px = (TH1D*)q_IMnpipi_woK0_wSid_n_Sp_sub_cs->ProjectionX();
  q_IMnpipi_woK0_wSid_n_Sp_sub_cs_px->SetLineColor(2);
  q_IMnpipi_woK0_wSid_n_Sp_sub_cs_px->SetMinimum(0);
  q_IMnpipi_woK0_wSid_n_Sp_sub_cs_px->GetYaxis()->SetTitle("C.S. (a.u.)");
  q_IMnpipi_woK0_wSid_n_Sp_sub_cs_px->GetYaxis()->CenterTitle();
  q_IMnpipi_woK0_wSid_n_Sp_sub_cs_px->SetMarkerStyle(20);
  q_IMnpipi_woK0_wSid_n_Sp_sub_cs_px->SetMarkerSize(1);
  q_IMnpipi_woK0_wSid_n_Sp_sub_cs_px->SetMarkerColor(2);
  q_IMnpipi_woK0_wSid_n_Sp_sub_cs_px->Draw("PE");

  //TCanvas *cq_IMnpipi_woK0_wSid_n_Sm_cs_px = new TCanvas("cq_IMnpipi_woK0_wSid_n_Sm_cs_px","q_IMnpipi_woK0_wSid_n_Sm_cs_px");
  TH1D *q_IMnpipi_woK0_wSid_n_Sm_sub_cs_px = (TH1D*)q_IMnpipi_woK0_wSid_n_Sm_sub_cs->ProjectionX();
  q_IMnpipi_woK0_wSid_n_Sm_sub_cs_px->SetLineColor(3);
  q_IMnpipi_woK0_wSid_n_Sm_sub_cs_px->Draw("PEsame");

  TH1D* cs_sub_sum = (TH1D*)q_IMnpipi_woK0_wSid_n_Sp_sub_cs_px->Clone();
  cs_sub_sum->Add(q_IMnpipi_woK0_wSid_n_Sm_sub_cs_px);
  cs_sub_sum->SetLineColor(4);
  cs_sub_sum->Draw("PEsame");
  */
  /*
  TCanvas *c1 = new TCanvas("c1","c1");
  c1->cd();
  q_IMnpipi_woK0_wSid_n_Sp_sub_cs_px->SetMarkerStyle(20);
  q_IMnpipi_woK0_wSid_n_Sp_sub_cs_px->SetMarkerSize(1);
  q_IMnpipi_woK0_wSid_n_Sp_sub_cs_px->SetMarkerColor(2);
  q_IMnpipi_woK0_wSid_n_Sp_sub_cs_px->SetMaximum(q_IMnpipi_woK0_wSid_n_Sp_sub_cs_px->GetMaximum()+50000.);
  q_IMnpipi_woK0_wSid_n_Sp_sub_cs_px->Draw("PE");


  TCanvas *c2 = new TCanvas("c2","c2");
  c2->cd();
  q_IMnpipi_woK0_wSid_n_Sm_sub_cs_px->SetMarkerStyle(20);
  q_IMnpipi_woK0_wSid_n_Sm_sub_cs_px->SetMarkerSize(1);
  q_IMnpipi_woK0_wSid_n_Sm_sub_cs_px->SetMinimum(0);
  q_IMnpipi_woK0_wSid_n_Sm_sub_cs_px->SetMaximum(q_IMnpipi_woK0_wSid_n_Sp_sub_cs_px->GetMaximum());
  q_IMnpipi_woK0_wSid_n_Sm_sub_cs_px->SetMarkerColor(3);
  q_IMnpipi_woK0_wSid_n_Sm_sub_cs_px->Draw("PE");
  
  TCanvas *c3 = new TCanvas("c3","c3");
  c3->cd();
  cs_sub_sum->SetMarkerStyle(20);
  cs_sub_sum->SetMarkerSize(1);
  cs_sub_sum->SetMarkerColor(4);
  cs_sub_sum->SetMinimum(0);
  cs_sub_sum->SetMaximum(q_IMnpipi_woK0_wSid_n_Sp_sub_cs_px->GetMaximum());
  cs_sub_sum->Draw("PE");
  */
   
  /*
  TCanvas *cq_IMnpipi_woK0_wSid_n_Sp_sub_cs_py = new TCanvas("cq_IMnpipi_woK0_wSid_n_Sp_sub_cs_py","q_IMnpipi_woK0_wSid_n_Sp_sub_cs_py");
  TH1D *q_IMnpipi_woK0_wSid_n_Sp_sub_cs_py = (TH1D*)q_IMnpipi_woK0_wSid_n_Sp_sub_cs->ProjectionY();
  q_IMnpipi_woK0_wSid_n_Sp_sub_cs_py->SetLineColor(2);
  q_IMnpipi_woK0_wSid_n_Sp_sub_cs_py->SetMinimum(0);
  q_IMnpipi_woK0_wSid_n_Sp_sub_cs_py->GetYaxis()->SetTitle("C.S. (a.u.)");
  q_IMnpipi_woK0_wSid_n_Sp_sub_cs_py->GetYaxis()->CenterTitle();
  q_IMnpipi_woK0_wSid_n_Sp_sub_cs_py->Draw("E");

  //TCanvas *cq_IMnpipi_woK0_wSid_n_Sm_cs_px = new TCanvas("cq_IMnpipi_woK0_wSid_n_Sm_cs_px","q_IMnpipi_woK0_wSid_n_Sm_cs_px");
  TH1D *q_IMnpipi_woK0_wSid_n_Sm_sub_cs_py = (TH1D*)q_IMnpipi_woK0_wSid_n_Sm_sub_cs->ProjectionY();
  q_IMnpipi_woK0_wSid_n_Sm_sub_cs_py->SetLineColor(3);
  q_IMnpipi_woK0_wSid_n_Sm_sub_cs_py->Draw("Esame");

  TH1D* cs_sub_sum_py = (TH1D*)q_IMnpipi_woK0_wSid_n_Sp_sub_cs_py->Clone();
  cs_sub_sum_py->Add(q_IMnpipi_woK0_wSid_n_Sm_sub_cs_py);
  cs_sub_sum_py->SetLineColor(4);
  cs_sub_sum_py->Draw("Esame");
  

  TCanvas *cq_IMnpipi_woK0_wSid_n_Sp_side_cs = new TCanvas("cq_IMnpipi_woK0_wSid_n_Sp_side_cs","q_IMnpipi_woK0_wSid_n_Sp_side_cs");
  cq_IMnpipi_woK0_wSid_n_Sp_side_cs->cd();
  TH2F *q_IMnpipi_woK0_wSid_n_Sp_side_cs = (TH2F*)q_IMnpipi_woK0_wSid_n_Sp_side[sidebandtype][0]->Clone("Sp_cs_side");
  q_IMnpipi_woK0_wSid_n_Sp_side_cs->Sumw2();
  q_IMnpipi_woK0_wSid_n_Sp_side_cs->Divide(acc_Sp_cal);
  q_IMnpipi_woK0_wSid_n_Sp_side_cs->SetMaximum(30000);
  q_IMnpipi_woK0_wSid_n_Sp_side_cs->Draw("colz");

  TCanvas *cq_IMnpipi_woK0_wSid_n_Sm_side_cs = new TCanvas("cq_IMnpipi_woK0_wSid_n_Sm_side_cs","q_IMnpipi_woK0_wSid_n_Sm_side_cs");
  cq_IMnpipi_woK0_wSid_n_Sm_side_cs->cd();
  TH2F *q_IMnpipi_woK0_wSid_n_Sm_side_cs = (TH2F*)q_IMnpipi_woK0_wSid_n_Sm_side[sidebandtype][0]->Clone("Sm_cs_side");
  q_IMnpipi_woK0_wSid_n_Sm_side_cs->Sumw2();
  q_IMnpipi_woK0_wSid_n_Sm_side_cs->Divide(acc_Sm_cal);
  q_IMnpipi_woK0_wSid_n_Sm_side_cs->SetMaximum(30000);
  q_IMnpipi_woK0_wSid_n_Sm_side_cs->Draw("colz");

  TCanvas *cq_IMnpipi_woK0_wSid_n_SpSm_side_cs = new TCanvas("cq_IMnpipi_woK0_wSid_n_SpSm_side_cs","q_IMnpipi_woK0_wSid_n_SpSm_side_cs");
  cq_IMnpipi_woK0_wSid_n_SpSm_side_cs->cd();
  TH2F* q_IMnpipi_woK0_wSid_n_SpSm_side_cs = (TH2F*)q_IMnpipi_woK0_wSid_n_Sp_side_cs->Clone("SpSm_cs_side");    
  q_IMnpipi_woK0_wSid_n_SpSm_side_cs->Scale(0.48);
  q_IMnpipi_woK0_wSid_n_SpSm_side_cs->Add(q_IMnpipi_woK0_wSid_n_Sm_side_cs);
  q_IMnpipi_woK0_wSid_n_SpSm_side_cs->Scale(1.48);
  q_IMnpipi_woK0_wSid_n_SpSm_side_cs->SetMinimum(0);
  q_IMnpipi_woK0_wSid_n_SpSm_side_cs->SetTitle("q_IMnpipi_woK0_wSid_n_SpSm_side_cs");
  q_IMnpipi_woK0_wSid_n_SpSm_side_cs->SetMaximum(60000);
  q_IMnpipi_woK0_wSid_n_SpSm_side_cs->Draw("colz");
  

  TCanvas *cIMnpip_IMnpipi_woK0_n = new TCanvas("cIMnpip_IMnpipi_woK0_n","IMnpip_IMnpipi_woK0_n");
  cIMnpip_IMnpipi_woK0_n->cd();
  IMnpip_IMnpipi_woK0_n->Draw("colz");
  
  TCanvas *cIMnpim_IMnpipi_woK0_n = new TCanvas("cIMnpim_IMnpipi_woK0_n","IMnpim_IMnpipi_woK0_n");
  cIMnpim_IMnpipi_woK0_n->cd();
  IMnpim_IMnpipi_woK0_n->Draw("colz");
  */

  //TCanvas *cq_IMnpipi_woK0_wSid_n_SpSm_side_cs_px = new TCanvas("cq_IMnpipi_woK0_wSid_n_SpSm_side_cs_px","q_IMnpipi_woK0_wSid_n_SpSm_side_cs_px");
  //cq_IMnpipi_woK0_wSid_n_SpSm_side_cs_px->cd();
  //TH1D* q_IMnpipi_woK0_wSid_n_SpSm_side_cs_px = q_IMnpipi_woK0_wSid_n_SpSm_side_cs->ProjectionX();
  //q_IMnpipi_woK0_wSid_n_SpSm_side_cs_px->Draw("E");
  
  //resolution evaluation
  if(Spmode){
    TCanvas *cdiff_IMnpipi_woK0_wSid_n_Sp = new TCanvas("cdiff_IMnpipi_woK0_wSid_n_Sp","diff_IMnpipi_woK0_wSid_n_Sp");
    cdiff_IMnpipi_woK0_wSid_n_Sp->cd();
    diff_IMnpipi_woK0_wSid_n_Sp->GetYaxis()->SetRangeUser(-0.1,0.1);
    diff_IMnpipi_woK0_wSid_n_Sp->Draw("colz");
    TProfile *pfxSp = (TProfile*)diff_IMnpipi_woK0_wSid_n_Sp->ProfileX("pfxSp",1,-1,"s");
    pfxSp->SetLineColor(2);
    pfxSp->SetMarkerStyle(33);
    pfxSp->Draw("same");
    TCanvas *cgaus = new TCanvas("cgaus","cgaus");
    double recomass[nbinIMnpipi];
    double cent[nbinIMnpipi];
    double cent_err[nbinIMnpipi];
    double sigma[nbinIMnpipi];
    double sigma_err[nbinIMnpipi];
    //cgaus->Divide(10,10);
    TH1D *px[nbinIMnpipi];
    for(int i=0;i<nbinIMnpipi;i++){
      px[i] = (TH1D*)diff_IMnpipi_woK0_wSid_n_Sp->ProjectionY(Form("px%d",i),i+1,i+2,"");
      recomass[i] = diff_IMnpipi_woK0_wSid_n_Sp->GetXaxis()->GetBinCenter(i+1);
      cgaus->cd(i+1);
      if(px[i]->GetEntries()>200){
        //px[i]->Draw("HE");
        px[i]->Fit("gaus","q"); 
        cent[i] = px[i]->GetFunction("gaus")->GetParameter(1);
        cent_err[i] = px[i]->GetFunction("gaus")->GetParError(1);
        sigma[i]= px[i]->GetFunction("gaus")->GetParameter(2);
        sigma_err[i]= px[i]->GetFunction("gaus")->GetParError(2);
      }else{
        cent[i]=-999.;
        cent_err[i]=-0.;
        sigma[i]=-999.;
        sigma_err[i]=0.;
      }
    }
    TCanvas *cfitmean = new TCanvas("cfitmean","fitmean");
    cfitmean->cd();
    TGraphErrors *grcent = new TGraphErrors(nbinIMnpipi,recomass,cent,0,cent_err); 
    grcent->SetTitle("gaussian mean");
    grcent->GetXaxis()->SetTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
    grcent->GetXaxis()->CenterTitle();
    grcent->GetYaxis()->SetTitle("mass center [GeV/c^{2}]");
    grcent->GetYaxis()->SetRangeUser(-0.01,0.03);
    grcent->GetYaxis()->CenterTitle();
    grcent->SetMarkerStyle(20);
    grcent->GetYaxis()->SetTitleOffset(1.5);
    grcent->Draw("APE");
    TCanvas *cfitsigma = new TCanvas("cfitsigma","fitsigma");
    cfitsigma->cd();
    TGraphErrors *grsigma = new TGraphErrors(nbinIMnpipi,recomass,sigma,0,sigma_err); 
    grsigma->SetTitle("mass resolution");
    grsigma->GetXaxis()->SetTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
    grsigma->GetXaxis()->CenterTitle();
    grsigma->GetYaxis()->SetTitle("mass resolution [GeV/c^{2}]");
    grsigma->GetYaxis()->CenterTitle();
    grsigma->GetYaxis()->SetTitleOffset(1.5);
    grsigma->GetYaxis()->SetRangeUser(0.0,0.03);
    grsigma->SetMarkerStyle(20);
    grsigma->Draw("APE");
    
    
    TCanvas *cdiff_q_woK0_wSid_n_Sp = new TCanvas("cdiff_q_woK0_wSid_n_Sp","diff_q_woK0_wSid_n_Sp");
    cdiff_q_woK0_wSid_n_Sp->cd();
    diff_q_woK0_wSid_n_Sp->GetYaxis()->SetRangeUser(-0.1,0.1);
    diff_q_woK0_wSid_n_Sp->Draw("colz");
    TProfile *pfxSp_q = (TProfile*)diff_q_woK0_wSid_n_Sp->ProfileX("pfxSp_q",1,-1,"s");
    pfxSp_q->SetLineColor(2);
    pfxSp_q->SetMarkerStyle(33);
    pfxSp_q->Draw("same");
    TCanvas *cgaus_q = new TCanvas("cgaus_q","cgaus_q");
    double recoq[nbinIMnpipi];
    double cent_q[nbinIMnpipi];
    double cent_q_err[nbinIMnpipi];
    double sigma_q[nbinIMnpipi];
    double sigma_q_err[nbinIMnpipi];
    //cgaus->Divide(10,10);
    TH1D *px_q[nbinIMnpipi];
    for(int i=0;i<nbinq;i++){
      px_q[i] = (TH1D*)diff_q_woK0_wSid_n_Sp->ProjectionY(Form("px_q%d",i),i+1,i+2,"");
      recoq[i] = diff_q_woK0_wSid_n_Sp->GetXaxis()->GetBinCenter(i+1);
      cgaus_q->cd(i+1);
      if(px_q[i]->GetEntries()>200){
        //px[i]->Draw("HE");
        px_q[i]->Fit("gaus","q"); 
        cent_q[i] = px_q[i]->GetFunction("gaus")->GetParameter(1);
        cent_q_err[i] = px_q[i]->GetFunction("gaus")->GetParError(1);
        sigma_q[i]= px_q[i]->GetFunction("gaus")->GetParameter(2);
        sigma_q_err[i]= px_q[i]->GetFunction("gaus")->GetParError(2);
      }else{
        cent_q[i]=-999.;
        cent_q_err[i]=-0.;
        sigma_q[i]=-999.;
        sigma_q_err[i]=0.;
      }
    }
    TCanvas *cfitmean_q = new TCanvas("cfitmean_q","fitmean_q");
    cfitmean_q->cd();
    TGraphErrors *grcent_q = new TGraphErrors(nbinq,recoq,cent_q,0,cent_q_err); 
    grcent_q->SetTitle("gaussian mean");
    grcent_q->GetXaxis()->SetTitle("Mom. Transfer [GeV/c]");
    grcent_q->GetXaxis()->CenterTitle();
    grcent_q->GetYaxis()->SetTitle("Mom. center [GeV/c]");
    grcent_q->GetYaxis()->SetRangeUser(-0.01,0.05);
    grcent_q->GetYaxis()->CenterTitle();
    grcent_q->SetMarkerStyle(20);
    grcent_q->GetYaxis()->SetTitleOffset(1.5);
    grcent_q->Draw("APE");
    TCanvas *cfitsigma_q = new TCanvas("cfitsigma_q","fitsigma_q");
    cfitsigma_q->cd();
    TGraphErrors *grsigma_q = new TGraphErrors(nbinq,recoq,sigma_q,0,sigma_q_err); 
    grsigma_q->SetTitle("Mom. resolution");
    grsigma_q->GetXaxis()->SetTitle("Mom. Transfer [GeV/c]");
    grsigma_q->GetXaxis()->CenterTitle();
    grsigma_q->GetYaxis()->SetTitle("Mom. resolution [GeV/c]");
    grsigma_q->GetYaxis()->CenterTitle();
    grsigma_q->GetYaxis()->SetTitleOffset(1.5);
    grsigma_q->GetYaxis()->SetRangeUser(0.0,0.05);
    grsigma_q->SetMarkerStyle(20);
    grsigma_q->Draw("APE");
  }
  
  if(Smmode){
    TCanvas *cdiff_IMnpipi_woK0_wSid_n_Sm = new TCanvas("cdiff_IMnpipi_woK0_wSid_n_Sm","diff_IMnpipi_woK0_wSid_n_Sm");
    cdiff_IMnpipi_woK0_wSid_n_Sm->cd();
    diff_IMnpipi_woK0_wSid_n_Sm->GetYaxis()->SetRangeUser(-0.1,0.1);
    diff_IMnpipi_woK0_wSid_n_Sm->Draw("colz");
    TProfile *pfxSm = diff_IMnpipi_woK0_wSid_n_Sm->ProfileX("pfxSm",1,-1,"s");
    pfxSm->SetLineColor(2);
    pfxSm->SetMarkerStyle(33);
    pfxSm->Draw("same");
    //TCanvas *cgaus = new TCanvas("cgaus","cgaus");
    double recomass[nbinIMnpipi];
    double cent[nbinIMnpipi];
    double cent_err[nbinIMnpipi];
    double sigma[nbinIMnpipi];
    double sigma_err[nbinIMnpipi];
    //cgaus->Divide(10,10);
    TH1D *px[nbinIMnpipi];
    for(int i=0;i<nbinIMnpipi;i++){
      px[i] = (TH1D*)diff_IMnpipi_woK0_wSid_n_Sm->ProjectionY(Form("px%d",i),i+1,i+2,"");
      recomass[i] = diff_IMnpipi_woK0_wSid_n_Sm->GetXaxis()->GetBinCenter(i+1);
      //cgaus->cd(i+1);
      if(px[i]->GetEntries()>200){
        //px[i]->Draw("HE");
        px[i]->Fit("gaus","q"); 
        cent[i] = px[i]->GetFunction("gaus")->GetParameter(1);
        cent_err[i] = px[i]->GetFunction("gaus")->GetParError(1);
        sigma[i]= px[i]->GetFunction("gaus")->GetParameter(2);
        sigma_err[i]= px[i]->GetFunction("gaus")->GetParError(2);
      }else{
        cent[i]=-999.;
        cent_err[i]=0;
        sigma[i]=-999.;
        sigma_err[i]=0;
      }
    }
    TCanvas *cfitmean = new TCanvas("cfitmean","fitmean");
    cfitmean->cd();
    TGraphErrors *grcent = new TGraphErrors(nbinIMnpipi,recomass,cent,0,cent_err); 
    grcent->SetTitle("gaussian mean");
    grcent->GetXaxis()->SetTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
    grcent->GetXaxis()->CenterTitle();
    grcent->GetYaxis()->SetTitle("mass center [GeV/c^{2}]");
    grcent->GetYaxis()->SetRangeUser(-0.02,0.03);
    grcent->GetYaxis()->CenterTitle();
    grcent->SetMarkerStyle(20);
    grcent->GetYaxis()->SetTitleOffset(1.5);
    grcent->Draw("APE");
    
    TCanvas *cfitsigma = new TCanvas("cfitsigma","fitsigma");
    cfitsigma->cd();
    TGraphErrors *grsigma = new TGraphErrors(nbinIMnpipi,recomass,sigma,0,sigma_err); 
    grsigma->SetTitle("mass resolution");
    grsigma->GetXaxis()->SetTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
    grsigma->GetXaxis()->CenterTitle();
    grsigma->GetYaxis()->SetTitle("mass resolution [GeV/c^{2}]");
    grsigma->GetYaxis()->CenterTitle();
    grsigma->GetYaxis()->SetTitleOffset(1.5);
    grsigma->GetYaxis()->SetRangeUser(0.0,0.03);
    grsigma->SetMarkerStyle(20);
    grsigma->Draw("APE");
    
    TCanvas *cdiff_q_woK0_wSid_n_Sm = new TCanvas("cdiff_q_woK0_wSid_n_Sm","diff_q_woK0_wSid_n_Sm");
    cdiff_q_woK0_wSid_n_Sm->cd();
    diff_q_woK0_wSid_n_Sm->GetYaxis()->SetRangeUser(-0.1,0.1);
    diff_q_woK0_wSid_n_Sm->Draw("colz");
    TProfile *pfxSm_q = (TProfile*)diff_q_woK0_wSid_n_Sm->ProfileX("pfxSm_q",1,-1,"s");
    pfxSm_q->SetLineColor(2);
    pfxSm_q->SetMarkerStyle(33);
    pfxSm_q->Draw("same");
    TCanvas *cgaus_q = new TCanvas("cgaus_q","cgaus_q");
    double recoq[nbinIMnpipi];
    double cent_q[nbinIMnpipi];
    double cent_q_err[nbinIMnpipi];
    double sigma_q[nbinIMnpipi];
    double sigma_q_err[nbinIMnpipi];
    //cgaus->Divide(10,10);
    TH1D *px_q[nbinIMnpipi];
    for(int i=0;i<nbinq;i++){
      px_q[i] = (TH1D*)diff_q_woK0_wSid_n_Sm->ProjectionY(Form("px_q%d",i),i+1,i+2,"");
      recoq[i] = diff_q_woK0_wSid_n_Sm->GetXaxis()->GetBinCenter(i+1);
      cgaus_q->cd(i+1);
      if(px_q[i]->GetEntries()>200){
        //px[i]->Draw("HE");
        px_q[i]->Fit("gaus","q"); 
        cent_q[i] = px_q[i]->GetFunction("gaus")->GetParameter(1);
        cent_q_err[i] = px_q[i]->GetFunction("gaus")->GetParError(1);
        sigma_q[i]= px_q[i]->GetFunction("gaus")->GetParameter(2);
        sigma_q_err[i]= px_q[i]->GetFunction("gaus")->GetParError(2);
      }else{
        cent_q[i]=-999.;
        cent_q_err[i]=-0.;
        sigma_q[i]=-999.;
        sigma_q_err[i]=0.;
      }
    }
    TCanvas *cfitmean_q = new TCanvas("cfitmean_q","fitmean_q");
    cfitmean_q->cd();
    TGraphErrors *grcent_q = new TGraphErrors(nbinq,recoq,cent_q,0,cent_q_err); 
    grcent_q->SetTitle("gaussian mean");
    grcent_q->GetXaxis()->SetTitle("Mom. Transfer [GeV/c]");
    grcent_q->GetXaxis()->CenterTitle();
    grcent_q->GetYaxis()->SetTitle("Mom. center [GeV/c]");
    grcent_q->GetYaxis()->SetRangeUser(-0.02,0.05);
    grcent_q->GetYaxis()->CenterTitle();
    grcent_q->SetMarkerStyle(20);
    grcent_q->GetYaxis()->SetTitleOffset(1.5);
    grcent_q->Draw("APE");
    TCanvas *cfitsigma_q = new TCanvas("cfitsigma_q","fitsigma_q");
    cfitsigma_q->cd();
    TGraphErrors *grsigma_q = new TGraphErrors(nbinq,recoq,sigma_q,0,sigma_q_err); 
    grsigma_q->SetTitle("Mom. resolution");
    grsigma_q->GetXaxis()->SetTitle("Mom. Transfer [GeV/c]");
    grsigma_q->GetXaxis()->CenterTitle();
    grsigma_q->GetYaxis()->SetTitle("Mom. resolution [GeV/c]");
    grsigma_q->GetYaxis()->CenterTitle();
    grsigma_q->GetYaxis()->SetTitleOffset(1.5);
    grsigma_q->GetYaxis()->SetRangeUser(0.0,0.05);
    grsigma_q->SetMarkerStyle(20);
    grsigma_q->Draw("APE");
  }

  
  //centering title of all histograms 
  f->cd();
  TIter nexthist(gDirectory->GetList());
  TH1F *h1 = nullptr;
  TH1D *h1d = nullptr;
  TH2F *h2 = nullptr;
  TObject *obj = nullptr;
  while( (obj = (TObject*)nexthist())!=nullptr  ){
    if(obj->InheritsFrom("TH1F")){
      h1 = (TH1F*) obj;
      h1->GetXaxis()->CenterTitle();
      //h1->GetXaxis()->SetTitleSize(0.05);
      //h1->GetXaxis()->SetTitleOffset(0.80);
      h1->GetYaxis()->SetTitleOffset(1.3);
    }
    if(obj->InheritsFrom("TH1D")){
      h1d = (TH1D*) obj;
      h1d->GetXaxis()->CenterTitle();
      //h1d->GetXaxis()->SetTitleSize(0.05);
      //h1d->GetXaxis()->SetTitleOffset(0.80);
      h1d->GetYaxis()->SetTitleOffset(1.3);
    }
    if(obj->InheritsFrom("TH2")){
      h2 = (TH2F*) obj;
      h2->GetXaxis()->CenterTitle();
      h2->GetYaxis()->CenterTitle();
      //h2->GetXaxis()->SetTitleSize(0.05);
      //h2->GetXaxis()->SetTitleOffset(0.80);
      //h2->GetYaxis()->SetTitleSize(0.05);
      h2->GetYaxis()->SetTitleOffset(1.1);
    }
  }

  TCanvas *c = nullptr;
  //TPDF *pdf = new TPDF(pdfname);
  TSeqCollection *SCol = gROOT->GetListOfCanvases();
  int size = SCol->GetSize();
  //TIter next(gROOT->GetListOfCanvases());
  TIter next(SCol);
  for(int i=0;i<size;i++){
  //while((c= (TCanvas*)next())){
    //pdf->NewPage();
    c= (TCanvas*)next();
    c->Draw();
    c->cd();
    TPaveText *pt;
    if(Spmode || Smmode){
      pt = new TPaveText(.80,0.90,0.98,0.99,"NDC");    
    }else{
      pt = new TPaveText(.80,0.90,0.98,0.99,"NDC");    
      //pt = new TPaveText(.74,.81,0.9,0.90,"NDC");
    }
    if(Spmode){
      pt->SetFillColor(kAzure-4);
      pt->AddText("MC #Sigma+#pi- mode");
    }
    else if(Smmode){
      pt->SetFillColor(kAzure-4);
      pt->AddText("MC #Sigma-#pi+ mode"); 
    }else{
      pt->AddText("Real Data");
      pt->SetFillColor(kCyan-9);
    }
    pt->SetBorderSize(1);
    pt->Draw();

    //gPad->SetLeftMargin(0.13);
    //gPad->SetBottomMargin(0.13);
    c->Modified();
    c->Update();
    if(i==0) c->Print(pdfname+"(",Form("Title:%s",c->GetTitle()));
    else if(i==size-1)c->Print(pdfname+")",Form("Title:%s",c->GetTitle())); 
    else c->Print(pdfname,Form("Title:%s",c->GetTitle())); 
  }
  //pdf->Close();
  std::cout << "closing pdf " << std::endl;
  //TString outname = std::string(filename);
  //outname.Replace(std::string(filename).size()-4,5,"out.root");
  //TFile *fout = new TFile(outname.Data(),"RECREATE");
  //fout->cd();
  //fout->Close();
  
}
