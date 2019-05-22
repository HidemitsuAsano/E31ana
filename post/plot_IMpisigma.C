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
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TLegend.h>
#include <TPDF.h>
#include <TPaveText.h>
#include <TLatex.h>
#include <TEfficiency.h>

#include "../src/GlobalVariables.h"
#include "anacuts.h"

const double pvalcut = 0.005;
const bool gridon=true;
const bool staton=true;
const bool UseKinFitVal = true;

void plot_IMpisigma(const char* filename="")
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
  gStyle->SetCanvasDefH(800); gStyle->SetCanvasDefW(900);
  //gStyle->SetTitleFontSize(0.1);
  

  std::cout << "infile " << filename <<std::endl;
  TString pdfname = std::string(filename);
  pdfname.Replace(std::string(filename).size()-4,5,"pdf");
  std::cout << "pdfname: " << pdfname << std::endl;
  std::cout << std::endl;
  std::cout << "Use Kin Fit Val ? " << std::endl;

  bool Spmode = (std::string(filename).find("Sp")!= std::string::npos);
  bool Smmode = (std::string(filename).find("Sm")!= std::string::npos);

  if(UseKinFitVal) std::cout << "Yes" << std::endl;
  else             std::cout << "No"  << std::endl;

  TH1::SetDefaultSumw2();
  //= = = = pipipnn final-sample tree = = = =//
  TLorentzVector *LVec_beam=nullptr;   // 4-momentum(beam)
  TLorentzVector *LVec_beam_Sp=nullptr;   // 4-momentum(beam),Sp mode assumption
  TLorentzVector *LVec_beam_Sm=nullptr;   // 4-momentum(beam),Sm mode assumption
  TLorentzVector *LVec_target=nullptr; // 4-momentum(target)
  TLorentzVector *LVec_pip=nullptr;    // 4-momentum(pi+)
  TLorentzVector *LVec_pim=nullptr;    // 4-momentum(pi-)
  TLorentzVector *LVec_n=nullptr;      // 4-momentum(neutron)
  TLorentzVector *LVec_n_Sp=nullptr;      // 4-momentum(neutron),Sp mode assumption
  TLorentzVector *LVec_n_Sm=nullptr;      // 4-momentum(neutron),Sm mode assumption
  TLorentzVector *mcmom_beam=nullptr;   // generated 4-momentum(beam)
  TLorentzVector *mcmom_pip=nullptr;    // generated 4-momentum(pi+)
  TLorentzVector *mcmom_pim=nullptr;    // generated 4-momentum(pi-)
  TLorentzVector *mcmom_ncds=nullptr;      // generated 4-momentum(neutron)
  TLorentzVector *mcmom_nmiss=nullptr;      // generated 4-momentum(neutron)
  double NeutralBetaCDH; // velocity of neutral particle on CDH
  double NeutralBetaCDH_vtx[2]; // velocity of neutral particle on CDH,0: Spmode 1:Smmode
  double dE;   // energy deposit on CDH
  TVector3 *vtx_reaction = nullptr; // vertex(reaction) 
  TVector3 *vtx_pip_beam = nullptr; //C.A.P of pip-beam beam side
  TVector3 *vtx_pim_beam = nullptr; //C.A.P of pim-beam beam side
  TVector3 *vtx_pip_cdc = nullptr;//C.A.P of pip-beam pip side
  TVector3 *vtx_pim_cdc = nullptr;//C.A.P of pim-beam pim side
  TVector3 *CA_pip = nullptr;//C.A.P of pip-pim pip side
  TVector3 *CA_pim = nullptr;//C.A.P of pip-pim pim side
  //int run_num;   // run number
  //int event_num; // event number
  //int block_num; // block number
  TLorentzVector *kfSpmode_mom_beam=nullptr;   // 4-momentum(beam) after kinematical refit for pi- Sigma+
  TLorentzVector *kfSpmode_mom_pip=nullptr;    // 4-momentum(pi+) after kinematical refit for pi- Sigma+
  TLorentzVector *kfSpmode_mom_pim=nullptr;    // 4-momentum(pi-) after kinematical refit for pi- Sigma+
  TLorentzVector *kfSpmode_mom_n=nullptr;      // 4-momentum(neutron) after kinematical refit for pi- Sigma+
  double kfSpmode_chi2;   // chi2 of kinematical refit
  double kfSpmode_NDF;    // NDF of kinematical refit
  double kfSpmode_status; // status of kinematical refit -> details can be found in this code
  double kfSpmode_pvalue; // p-value of kinematical refit
  TLorentzVector *kfSmmode_mom_beam=nullptr;   // 4-momentum(beam) after kinematical refit for pi+ Sigma-
  TLorentzVector *kfSmmode_mom_pip=nullptr;    // 4-momentum(pi+) after kinematical refit for pi+ Sigma-
  TLorentzVector *kfSmmode_mom_pim=nullptr;    // 4-momentum(pi-) after kinematical refit for pi+ Sigma-
  TLorentzVector *kfSmmode_mom_n=nullptr;      // 4-momentum(neutron) after kinematical refit for pi+ Sigma-
  double kfSmmode_chi2;   // chi2 of kinematical refit
  double kfSmmode_NDF;    // NDF of kinematical refit
  double kfSmmode_status; // status of kinematical refit -> details can be found in this code
  double kfSmmode_pvalue; // p-value of kinematical refit
  int kf_flag; // flag of correct pair reconstruction, etc

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
  TH2F* IMnpim_IMnpip_dE_woK0_n_Sm;//Smmode, avoiding crossing-point
  TH2F* IMnpim_IMnpip_dE_woK0_n_side;//Side band for Sp + Sm mode
  TH2F* IMnpim_IMnpip_dE_woK0_n_Sp_side[2];//low mass,high mass 
  TH2F* IMnpim_IMnpip_dE_woK0_n_Sp_side_cut[2];//low mass,high mass 
  TH2F* IMnpim_IMnpip_dE_woK0_n_Sp_side_sum;//low mass,high mass 
  TH2F* IMnpim_IMnpip_dE_woK0_n_Sm_side[2];//low mass,high mass
  TH2F* IMnpim_IMnpip_dE_woK0_n_Sm_side_cut[2];//low mass,high mass
  TH2F* IMnpim_IMnpip_dE_woK0_n_Sm_side_sum;//low mass,high mass 
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
  TH2F* q_IMnpipi_woK0_wSid_n_Sp;
  TH2F* q_IMnpipi_woK0_wSid_n_Sp_acc;//fine bins
  TH2F* q_IMnpipi_woK0_wSid_n_Sp_side;
  TH2F* q_IMnpipi_woK0_wSid_n_Sm;
  TH2F* q_IMnpipi_woK0_wSid_n_Sm_acc;//fine bins
  TH2F* q_IMnpipi_woK0_wSid_n_Sm_side;
  //TH2F* q_IMnpipi_woK0_wSid_n_wocross;

  TH2F* q_IMnpipi_wSid_n_side;//side band method
  TH2F* q_IMnpipi_woK0_wSid_n_side;//side band method
  TH2F* nmom_IMnpipi_woK0_wSid_n;
  TH2F* q_pippim_n;

  // w/ kinematic fit
  TH2F* dE_betainv_fid_kin[2];//
  TH2F* dE_MMom_fid_beta_woK0_kin[2];
  TH2F* dE_MMass_fid_beta_woK0_kin[2];
  TH2F* MMom_MMass_woK0_kin[2];
  TH2F* IMnpim_IMnpip_dE_woK0_kin[2];
  TH2F* MMnmiss_IMnpip_dE_woK0_kin[2];//not yet
  TH2F* MMnmiss_IMnpim_dE_woK0_kin[2];//not yet
  TH2F* MMnpip_MMnpim_woK0_kin[2];
  TH2F* dE_IMnpipi_woK0_kin[2];
  TH2F* Cosn_IMnpipi_woK0_kin[2];
  TH2F* MMnmiss_IMnpipi_woK0_kin[2];
  TH2F* q_IMnpipi_kin[2];
  TH2F* q_IMnpipi_woK0_kin[2];
  TH2F* nmom_IMnpipi_woK0_kin[2];
  TH2F* q_pippim_n_kin[2];
  const char smode[][4]={"Sp","Sm"};
  
  const int nbinIMnpipi = 100;//1-2 GeV/c^2
  const int nbinq = 25;//0-1.5 GeV/c
  const int nbinIMnpi = 200; //1-2 GeV/c^2
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
  
  IMnpim_IMnpip_dE_woK0_n_Sp_side_sum = new TH2F(Form("IMnpim_IMnpip_dE_woK0_n_Sp_side_sum"),Form("IMnpim_IMnpip_dE_woK0_n_Sp_side_sum"),nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_woK0_n_Sp_side_sum->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_woK0_n_Sp_side_sum->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  
  IMnpim_IMnpip_dE_woK0_n_Sm_side_sum = new TH2F(Form("IMnpim_IMnpip_dE_woK0_n_Sm_side_sum"),Form("IMnpim_IMnpip_dE_woK0_n_Sm_side_sum"), nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_woK0_n_Sm_side_sum->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_woK0_n_Sm_side_sum->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  
  const char  lh[][6]={"low","high"};
  for(int i=0;i<2;i++){
    IMnpim_IMnpip_dE_woK0_n_Sp_side[i] = new TH2F(Form("IMnpim_IMnpip_dE_woK0_n_Sp_side_%s",lh[i]),Form("IMnpim_IMnpip_dE_woK0_n_Sp_side_%s",lh[i]), nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
    IMnpim_IMnpip_dE_woK0_n_Sp_side[i]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
    IMnpim_IMnpip_dE_woK0_n_Sp_side[i]->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
    IMnpim_IMnpip_dE_woK0_n_Sp_side_cut[i] = new TH2F(Form("IMnpim_IMnpip_dE_woK0_n_Sp_side_%s_cut",lh[i]),Form("IMnpim_IMnpip_dE_woK0_n_Sp_side_%s_cut",lh[i]), nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
    IMnpim_IMnpip_dE_woK0_n_Sp_side_cut[i]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
    IMnpim_IMnpip_dE_woK0_n_Sp_side_cut[i]->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
    
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
  
  IMnpim_IMnpip_dE_woK0_n_Sm = new TH2F(Form("IMnpim_IMnpip_dE_woK0_n_Sm"),Form("IMnpim_IMnpip_dE_woK0_n_Sm"), nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_woK0_n_Sm->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_woK0_n_Sm->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
    
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
  MMnmiss_IMnpipi_woK0_wSid_n->SetYTitle("Miss Mass. [GeV/c^{2}]");
    
  q_IMnpipi_wSid_n = new TH2F(Form("q_IMnpipi_wSid_n"),Form("q_IMnpipi_wSid_n"),nbinIMnpipi,1,2,nbinq,0,1.5);
  q_IMnpipi_wSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wSid_n->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_IMnpipi_woK0_wSid_n = new TH2F(Form("q_IMnpipi_woK0_wSid_n"),Form("q_IMnpipi_woK0_wSid_n"), nbinIMnpipi,1,2, nbinq,0,1.5);
  q_IMnpipi_woK0_wSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_wSid_n->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_IMnpipi_woK0_wSid_n_Sp = new TH2F(Form("q_IMnpipi_woK0_wSid_n_Sp"),Form("q_IMnpipi_woK0_wSid_n_Sp"), nbinIMnpipi,1,2, nbinq,0,1.5);
  q_IMnpipi_woK0_wSid_n_Sp->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_wSid_n_Sp->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_IMnpipi_woK0_wSid_n_Sp_acc = new TH2F(Form("q_IMnpipi_woK0_wSid_n_Sp_acc"),Form("q_IMnpipi_woK0_wSid_n_Sp_acc"),500,1,2,300,0,1.5);
  q_IMnpipi_woK0_wSid_n_Sp_acc->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_wSid_n_Sp_acc->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_IMnpipi_woK0_wSid_n_Sp_side = new TH2F(Form("q_IMnpipi_woK0_wSid_n_Sp_side"),Form("q_IMnpipi_woK0_wSid_n_Sp_side"), nbinIMnpipi,1,2, nbinq,0,1.5);
  q_IMnpipi_woK0_wSid_n_Sp_side->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_wSid_n_Sp_side->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_IMnpipi_woK0_wSid_n_Sm = new TH2F(Form("q_IMnpipi_woK0_wSid_n_Sm"),Form("q_IMnpipi_woK0_wSid_n_Sm"), nbinIMnpipi,1,2, nbinq,0,1.5);
  q_IMnpipi_woK0_wSid_n_Sm->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_wSid_n_Sm->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_IMnpipi_woK0_wSid_n_Sm_acc = new TH2F(Form("q_IMnpipi_woK0_wSid_n_Sm_acc"),Form("q_IMnpipi_woK0_wSid_n_Sm_acc"),500,1,2,300,0,1.5);
  q_IMnpipi_woK0_wSid_n_Sm_acc->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_wSid_n_Sm_acc->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_IMnpipi_woK0_wSid_n_Sm_side = new TH2F(Form("q_IMnpipi_woK0_wSid_n_Sm_side"),Form("q_IMnpipi_woK0_wSid_n_Sm_side"),nbinIMnpipi,1,2,nbinq,0,1.5);
  q_IMnpipi_woK0_wSid_n_Sm_side->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_wSid_n_Sm_side->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_IMnpipi_wSid_n_side = new TH2F(Form("q_IMnpipi_wSid_n_side"),Form("q_IMnpipi_wSid_n_side"), nbinIMnpipi,1,2, nbinq,0,1.5);
  q_IMnpipi_wSid_n_side->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wSid_n_side->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_IMnpipi_woK0_wSid_n_side = new TH2F(Form("q_IMnpipi_woK0_wSid_n_side"),Form("q_IMnpipi_woK0_wSid_n_side"), nbinIMnpipi,1,2, nbinq,0,1.5);
  q_IMnpipi_woK0_wSid_n_side->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_wSid_n_side->SetYTitle("Mom. Transfer [GeV/c]");
  
  nmom_IMnpipi_woK0_wSid_n = new TH2F(Form("nmom_IMnpipi_woK0_wSid_n"),Form("nmom_IMnpipi_woK0_wSid_n"), nbinIMnpipi,1,2,100,0,1.0);
  nmom_IMnpipi_woK0_wSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  nmom_IMnpipi_woK0_wSid_n->SetYTitle("nmom  [GeV/c]");
   
  q_pippim_n = new TH2F("q_pippim_n","q_pippim_n",500,0,1, 100,0,1.5);
  q_pippim_n->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_pippim_n->SetYTitle("Mom. Transfer [GeV/c]");

  for(int imode=0;imode<2;imode++){
    dE_betainv_fid_kin[imode] = new TH2F(Form("dE_betainv_fid_kin_%s",smode[imode]),Form("dE_betainv_fid_kin_%s",smode[imode]),1000, 0, 50, nbindE, 0, 50);
    dE_betainv_fid_kin[imode]->SetXTitle("1/#beta");
    dE_betainv_fid_kin[imode]->SetYTitle("dE [MeVee]");

    
    dE_MMom_fid_beta_woK0_kin[imode] = new TH2F(Form("dE_MMom_fid_beta_woK0_kin_%s",smode[imode]),Form("dE_MMom_fid_beta_woK0_kin_%s",smode[imode]),100, 0, 1.5, nbindE, 0, 50);
    dE_MMom_fid_beta_woK0_kin[imode]->SetXTitle("Missing Mom. [GeV/c]");
    dE_MMom_fid_beta_woK0_kin[imode]->SetYTitle("dE [MeVee]");


    dE_MMass_fid_beta_woK0_kin[imode] = new TH2F(Form("dE_MMass_fid_beta_woK0_kin_%s",smode[imode]),Form("dE_MMass_fid_beta_woK0_kin_%s",smode[imode]), 140, 0.4, 1.8, nbindE, 0, 50);
    dE_MMass_fid_beta_woK0_kin[imode]->SetXTitle("Missing mass [GeV/c^{2}]");
    dE_MMass_fid_beta_woK0_kin[imode]->SetYTitle("dE [MeVee]");
  

    MMom_MMass_woK0_kin[imode] = new TH2F(Form("MMom_MMass_woK0_kin_%s",smode[imode]),Form("MMom_MMass_woK0_kin_%s",smode[imode]), 140, 0.4, 1.8, 100, 0, 1.5);
    MMom_MMass_woK0_kin[imode]->SetXTitle("Missing Mass [GeV/c^{2}]");
    MMom_MMass_woK0_kin[imode]->SetYTitle("Missing Mom. [GeV/c]");
    

    IMnpim_IMnpip_dE_woK0_kin[imode] = new TH2F(Form("IMnpim_IMnpip_dE_woK0_kin_%s",smode[imode]), Form("IMnpim_IMnpip_dE_woK0_kin_%s",smode[imode]), nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
    IMnpim_IMnpip_dE_woK0_kin[imode]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
    IMnpim_IMnpip_dE_woK0_kin[imode]->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  

    MMnpip_MMnpim_woK0_kin[imode] = new TH2F(Form("MMnpip_MMnpim_woK0_kin_%s",smode[imode]),Form("MMnpip_MMnpim_woK0_kin_%s",smode[imode]),70, 1, 1.7, 70, 1, 1.7);
    MMnpip_MMnpim_woK0_kin[imode]->SetXTitle("Miss. Mass(n#pi^{+}) [GeV/c^{2}]");
    MMnpip_MMnpim_woK0_kin[imode]->SetYTitle("Miss. Mass(n#pi^{-}) [GeV/c^{2}]");

 
    dE_IMnpipi_woK0_kin[imode] = new TH2F(Form("dE_IMnpipi_woK0_kin_%s",smode[imode]),Form("dE_IMnpipi_woK0_kin_%s",smode[imode]),100, 1, 2, nbindE, 0, 50);
    dE_IMnpipi_woK0_kin[imode]->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
    dE_IMnpipi_woK0_kin[imode]->SetYTitle("dE [MeVee]");

  
    Cosn_IMnpipi_woK0_kin[imode] = new TH2F(Form("Cosn_IMnpipi_woK0_kin_%s",smode[imode]),Form("dE_Cosn_IMnpipi_woK0_kin_%s",smode[imode]),100, 1, 2, 50, -1, 1);
    Cosn_IMnpipi_woK0_kin[imode]->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
    Cosn_IMnpipi_woK0_kin[imode]->SetYTitle("cos#theta_{n} (CM)");
  

    MMnmiss_IMnpipi_woK0_kin[imode] = new TH2F(Form("MMnmiss_IMnpipi_woK0_kin_%s",smode[imode]),Form("MMnmiss_IMnpipi_woK0_kin_%s",smode[imode]),100,1,2,100,0,1.5);
    MMnmiss_IMnpipi_woK0_kin[imode]->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
    MMnmiss_IMnpipi_woK0_kin[imode]->SetYTitle("Miss Mass. [GeV/c^{2}]");

    q_IMnpipi_kin[imode] = new TH2F(Form("q_IMnpipi_kin_%s",smode[imode]),Form("q_IMnpipi_kin_%s",smode[imode]),nbinIMnpipi,1,2, nbinq,0,1.5);
    q_IMnpipi_kin[imode]->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
    q_IMnpipi_kin[imode]->SetYTitle("Mom. Transfer [GeV/c]");
    
    q_IMnpipi_woK0_kin[imode] = new TH2F(Form("q_IMnpipi_woK0_kin_%s",smode[imode]),Form("q_IMnpipi_woK0_kin_%s",smode[imode]), nbinIMnpipi,1,2, nbinq,0,1.5);
    q_IMnpipi_woK0_kin[imode]->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
    q_IMnpipi_woK0_kin[imode]->SetYTitle("Mom. Transfer [GeV/c]");

    nmom_IMnpipi_woK0_kin[imode] = new TH2F(Form("nmom_IMnpipi_woK0_kin_%s",smode[imode]),Form("nmom_IMnpipi_woK0_kin_%s",smode[imode]), nbinIMnpipi,1,2,100,0,1.0);
    nmom_IMnpipi_woK0_kin[imode]->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
    nmom_IMnpipi_woK0_kin[imode]->SetYTitle("nmom  [GeV/c]");
  };

  TH2F *KFpvalue_vs = new TH2F("KFpvalue_vs", "KFpvalue_vs", 500, 0, 1, 500, 0, 1 );
  KFpvalue_vs->SetXTitle("P-value (#Sigma^{+} mode)");
  KFpvalue_vs->SetYTitle("P-value (#Sigma^{-} mode)");
  
  TH2F *KFchi2ndf_vs = new TH2F("KFchi2ndf_vs", "KFchi2ndf_vs", 300, 0, 300, 300, 0, 300 );
  KFchi2ndf_vs->SetXTitle("chi2/NDF (#Sigma^{+} mode kin. fit)");
  KFchi2ndf_vs->SetYTitle("chi2/NDF (#Sigma^{-} mode kin. fit)");

  TH1F *nmom = new TH1F("nmom", "nmom", 50, 0, 1.0);
  nmom->SetXTitle("mom. [GeV/c]");
  
  TH1F *nmom_kin = new TH1F("nmom_kin", "nmom_kin", 50, 0, 1.0);
  nmom_kin->SetXTitle("mom. [GeV/c]");
  
  TH2F *dE_nmom = new TH2F("dE_nmom", "dE_nmom", 50, 0, 1.0, 200 , 0, 50);
  dE_nmom->SetXTitle("mom. [GeV/c]");
  dE_nmom->SetYTitle("dE. [MeVee]");
  
  TH2F *dE_nmom_kin = new TH2F("dE_nmom_kin", "dE_nmom_kin", 50, 0, 1.0, 200 , 0, 50);
  dE_nmom_kin->SetXTitle("mom. [GeV/c]");
  dE_nmom_kin->SetYTitle("dE. [MeVee]");
  
  TH1F *mnmom = new TH1F("mnmom", "mnmom", 100, 0, 2.0);
  mnmom->SetXTitle("mom. [GeV/c]");
  
  TH1F *mnmom_kin = new TH1F("mnmom_kin", "mnmom_kin", 100, 0, 2.0);
  mnmom_kin->SetXTitle("mom. [GeV/c]");
  
  TH1F *npipmom = new TH1F("npipmom", "npipmom", 150, 0, 3.0);
  npipmom->SetXTitle("mom. [GeV/c]");
  
  TH1F *npipmom_kin = new TH1F("npipmom_kin", "npipmom_kin", 150, 0, 3.0);
  npipmom_kin->SetXTitle("mom. [GeV/c]");
  
  TH1F *npimmom = new TH1F("npimmom", "npimmom", 150, 0, 3.0);
  npimmom->SetXTitle("mom. [GeV/c]");
  
  TH1F *npimmom_kin = new TH1F("npimmom_kin", "npimmom_kin", 150, 0, 3.0);
  npimmom_kin->SetXTitle("mom. [GeV/c]");
  
  //DCA analysis
  TH1F* DCA_pip_beam = new TH1F("DCA_pip_beam","DCA_pip_beam",3000,0,30);
  DCA_pip_beam->SetXTitle("DCA [cm]");
  TH1F* DCA_pim_beam = new TH1F("DCA_pim_beam","DCA_pim_beam",3000,0,30);
  DCA_pim_beam->SetXTitle("DCA [cm]");
  TH1F* DCA_pip_pim = new TH1F("DCA_pip_pim","DCA_pip_pim",3000,0,30);
  DCA_pip_pim->SetXTitle("DCA [cm]");
  TH1F* DCA_pip_beam_kin[2];
  TH1F* DCA_pim_beam_kin[2];
  TH1F* DCA_pip_pim_kin[2];
  for(int imode=0;imode<2;imode++){
    DCA_pip_beam_kin[imode] = new TH1F(Form("DCA_pip_beam_kin_%s",smode[imode]),Form("DCA_pip_beam_kin_%s",smode[imode]),3000,0,30);
    DCA_pip_beam_kin[imode]->SetXTitle("DCA [cm]");
    DCA_pim_beam_kin[imode] = new TH1F(Form("DCA_pim_beam_kin_%s",smode[imode]),Form("DCA_pim_beam_kin_%s",smode[imode]),3000,0,30);
    DCA_pim_beam_kin[imode]->SetXTitle("DCA [cm]");
    DCA_pip_pim_kin[imode] = new TH1F(Form("DCA_pip_pim_kin_%s",smode[imode]),Form("DCA_pip_pim_kin_%s",smode[imode]),3000,0,30);
    DCA_pip_pim_kin[imode]->SetXTitle("DCA [cm]");
  }



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
    if(!UseKinFitVal){
      LVec_n_miss_vtx[0] = *LVec_target+*LVec_beam_Sp-*LVec_pip-*LVec_pim-*LVec_n_Sp;
      LVec_n_miss_vtx[1] = *LVec_target+*LVec_beam_Sm-*LVec_pip-*LVec_pim-*LVec_n_Sm;
    }else{
      LVec_n_miss_vtx[0] = *LVec_target+*kfSpmode_mom_beam-*kfSpmode_mom_pip-*kfSpmode_mom_pim-*kfSpmode_mom_n;
      LVec_n_miss_vtx[1] = *LVec_target+*kfSmmode_mom_beam-*kfSmmode_mom_pip-*kfSmmode_mom_pim-*kfSmmode_mom_n;
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
    TLorentzVector qkn = *LVec_beam-LVec_n_miss;
    TLorentzVector qkn_vtx[2];
    if(!UseKinFitVal){ 
      qkn_vtx[0] = *LVec_beam_Sp-LVec_n_miss_vtx[0];
      qkn_vtx[1] = *LVec_beam_Sm-LVec_n_miss_vtx[1];
    }else{
      qkn_vtx[0] = *kfSpmode_mom_beam-LVec_n_miss_vtx[0];
      qkn_vtx[1] = *kfSmmode_mom_beam-LVec_n_miss_vtx[1];
    }
    // calc pi+pi- //
    TLorentzVector LVec_pip_pim = *LVec_pip+*LVec_pim;
    TLorentzVector LVec_pip_pim_kinval[2];
    LVec_pip_pim_kinval[0]=*kfSpmode_mom_pip+*kfSpmode_mom_pim;
    LVec_pip_pim_kinval[1]=*kfSmmode_mom_pip+*kfSmmode_mom_pim;
     

    // calc pi+n //
    TLorentzVector LVec_pip_n = *LVec_pip+*LVec_n;
    TLorentzVector LVec_pip_n_vtx[2];
    if(!UseKinFitVal){
      LVec_pip_n_vtx[0] = *LVec_pip+*LVec_n_Sp;
      LVec_pip_n_vtx[1] = *LVec_pip+*LVec_n_Sm;
    }else{
      LVec_pip_n_vtx[0] = *kfSpmode_mom_pip+*kfSpmode_mom_n;
      LVec_pip_n_vtx[1] = *kfSmmode_mom_pip+*kfSmmode_mom_n;
    }
    // calc pi-n //
    TLorentzVector LVec_pim_n = *LVec_pim+*LVec_n;
    TLorentzVector LVec_pim_n_vtx[2];
    if(!UseKinFitVal){
      LVec_pim_n_vtx[0]= *LVec_pim+*LVec_n_Sp;
      LVec_pim_n_vtx[1]= *LVec_pim+*LVec_n_Sm;
    }else{
      LVec_pim_n_vtx[0]= *kfSpmode_mom_pim+*kfSpmode_mom_n;
      LVec_pim_n_vtx[1]= *kfSmmode_mom_pim+*kfSmmode_mom_n;
    }

    // calc missing Sp //
    TLorentzVector LVec_pip_n_miss = *LVec_target+*LVec_beam-*LVec_pim-*LVec_n;
    TLorentzVector LVec_pip_n_miss_vtx[2];
    if(!UseKinFitVal){
      LVec_pip_n_miss_vtx[0]= *LVec_target+*LVec_beam_Sp-*LVec_pim-*LVec_n_Sp;                          
      LVec_pip_n_miss_vtx[1]= *LVec_target+*LVec_beam_Sm-*LVec_pim-*LVec_n_Sm;
    }else{
      LVec_pip_n_miss_vtx[0]= *LVec_target+*kfSpmode_mom_beam-*kfSpmode_mom_pim-*kfSpmode_mom_n;                          
      LVec_pip_n_miss_vtx[1]= *LVec_target+*kfSmmode_mom_beam-*kfSmmode_mom_pim-*kfSmmode_mom_n;
    }

    // calc missing Sm //
    TLorentzVector LVec_pim_n_miss = *LVec_target+*LVec_beam-*LVec_pip-*LVec_n;
    TLorentzVector LVec_pim_n_miss_vtx[2];
    if(!UseKinFitVal){
      LVec_pim_n_miss_vtx[0] = *LVec_target+*LVec_beam_Sp-*LVec_pip-*LVec_n_Sp;
      LVec_pim_n_miss_vtx[1] = *LVec_target+*LVec_beam_Sm-*LVec_pip-*LVec_n_Sm;   
    }else{
      LVec_pim_n_miss_vtx[0] = *LVec_target+*kfSpmode_mom_beam-*kfSpmode_mom_pip-*kfSpmode_mom_n;
      LVec_pim_n_miss_vtx[1] = *LVec_target+*kfSmmode_mom_beam-*kfSmmode_mom_pip-*kfSmmode_mom_n;     
    }
    
    // calc pi+pi-n //
    TLorentzVector LVec_pip_pim_n = *LVec_pip+*LVec_pim+*LVec_n;
    TLorentzVector LVec_pip_pim_n_vtx[2];
    if(!UseKinFitVal){
      LVec_pip_pim_n_vtx[0] = *LVec_pip+*LVec_pim+*LVec_n_Sp;
      LVec_pip_pim_n_vtx[1] = *LVec_pip+*LVec_pim+*LVec_n_Sm;
    }else{
      LVec_pip_pim_n_vtx[0] = *kfSpmode_mom_pip+*kfSpmode_mom_pim+*kfSpmode_mom_n;
      LVec_pip_pim_n_vtx[1] = *kfSmmode_mom_pip+*kfSmmode_mom_pim+*kfSmmode_mom_n;
    }
    TLorentzVector LVec_pip_pim_n_CM = LVec_pip_pim_n;
    LVec_pip_pim_n_CM.Boost(-boost);
    //double cos_X = LVec_pip_pim_n_CM.Vect().Dot(LVec_beam_CM.Vect())/(LVec_pip_pim_n_CM.Vect().Mag()*LVec_beam_CM.Vect().Mag());

    //if(qkn.P()<anacuts::qvalcut || qkn.P()>0.70) continue;
    //if(qkn.P()>anacuts::qvalcut ) continue;
    //if(qkn.P()>0.70 ) continue;
    //if(LVec_pip_pim_n.M() < 1.45) continue;
    //if(LVec_pip_pim_n.M() > 1.55) continue;
    //if(dcapippim < 1) continue;
    //if(LVec_pip_pim_n.M()<1.45 ) continue;
    //double chi2 = kfSpmode_chi2<kfSmmode_chi2 ? kfSpmode_chi2:kfSmmode_chi2;
    double pvalue = kfSmmode_pvalue<kfSpmode_pvalue ? kfSpmode_pvalue:kfSmmode_pvalue;

    bool K0rejectFlag=false;
    bool MissNFlag=false;
    bool NBetaOK=false;
    bool NdEOK=false;
    bool SigmaPFlag=false;
    bool SigmaMFlag=false;
    bool SigmawidePFlag=false;
    bool SigmawideMFlag=false;
    bool SidebandFlag=false;
    bool SigmaPsideFlag=false;
    bool SigmaPsideLowFlag=false;
    bool SigmaPsideHighFlag=false;
    bool SigmaMsideFlag=false;
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
    //triangular cuts
    //
    //Sigma cross region top
    if(  (MassNPim >  (MassNPip - anacuts::Sigmap_center) + anacuts::Sigmam_center )  
      && (MassNPim > -(MassNPip - anacuts::Sigmap_center) + anacuts::Sigmam_center ) ) SigmaCrossFlagTop=true;
    
    //Sigma cross region bottom
    if(  (MassNPim <  (MassNPip - anacuts::Sigmap_center) + anacuts::Sigmam_center )  
      && (MassNPim < -(MassNPip - anacuts::Sigmap_center) + anacuts::Sigmam_center ) ) SigmaCrossFlagBottom=true;

    //Sigma cross region right
    if(  (MassNPim <  (MassNPip - anacuts::Sigmap_center) + anacuts::Sigmam_center )  
      && (MassNPim > -(MassNPip - anacuts::Sigmap_center) + anacuts::Sigmam_center ) ) SigmaCrossFlagRight=true;
   
    //Sigma cross region left
    if(  (MassNPim >  (MassNPip - anacuts::Sigmap_center) + anacuts::Sigmam_center )  
      && (MassNPim < -(MassNPip - anacuts::Sigmap_center) + anacuts::Sigmam_center ) ) SigmaCrossFlagLeft=true;
    

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
    if( (SigmaPsideLowFlag  && !SigmaCrossPsideLowFlagLeft  && !SigmaCrossPsideLowFlagRight) 
    ||  (SigmaPsideHighFlag && !SigmaCrossPsideHighFlagLeft  && !SigmaCrossPsideHighFlagRight)
    ) SigmaPsideFlag=true;

    //Sigma- production side band low mass side
    if( ((anacuts::Sigmam_sidelow_MIN<MassNPim) && 
         (MassNPim <  anacuts::Sigmam_sidelow_MAX))) SigmaMsideLowFlag=true;
    
    //Sigma- production side band high mass side
    if( ((anacuts::Sigmam_sidehigh_MIN<MassNPim) && 
         (MassNPim <  anacuts::Sigmam_sidehigh_MAX))) SigmaMsideHighFlag=true;
    
    //Sigma- production side band low or high mass side
    //if(SigmaMsideLowFlag || SigmaMsideHighFlag) SigmaMsideFlag=true;
    if( (SigmaMsideLowFlag  && !SigmaCrossMsideLowFlagTop  && !SigmaCrossMsideLowFlagBottom) 
    ||  (SigmaMsideHighFlag && !SigmaCrossMsideHighFlagTop && !SigmaCrossMsideHighFlagBottom)
    ) SigmaMsideFlag=true;

    if( (SigmaPsideFlag || SigmaMsideFlag)) SidebandFlag=true;

    if(anacuts::neutron_MIN<nmiss_mass && nmiss_mass<anacuts::neutron_MAX ) MissNFlag=true;
    
    //K0 rejection using original momentum
    if( (LVec_pip_pim.M()<anacuts::pipi_MIN || anacuts::pipi_MAX<LVec_pip_pim.M())) K0rejectFlag=true;
       
    if(K0rejectFlag && NBetaOK && NdEOK){
      KFpvalue_vs->Fill(kfSpmode_pvalue,kfSmmode_pvalue);
      KFchi2ndf_vs->Fill(kfSpmode_chi2/kfSpmode_NDF, kfSmmode_chi2/kfSmmode_NDF);
    }
    //w/o kinfit
    if(K0rejectFlag ){
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
      if(SidebandFlag){
        IMnpim_IMnpip_dE_woK0_n_side->Fill(LVec_pip_n.M(),LVec_pim_n.M());
      }
      if(SigmaPsideLowFlag && !SigmaCrossPsideLowFlagLeft && !SigmaCrossPsideLowFlagRight){
        IMnpim_IMnpip_dE_woK0_n_Sp_side[0]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
        //if(!SigmawideMFlag){
          IMnpim_IMnpip_dE_woK0_n_Sp_side_cut[0]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
        //}
      }
      if(SigmaPsideHighFlag && (SigmaCrossPsideHighFlagTop || SigmaCrossPsideHighFlagBottom)) {
        IMnpim_IMnpip_dE_woK0_n_Sp_side[1]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
        //if(!SigmawideMFlag){
          IMnpim_IMnpip_dE_woK0_n_Sp_side_cut[1]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
        //}
      }
      if(SigmaMsideLowFlag && (SigmaCrossMsideLowFlagLeft || SigmaCrossMsideLowFlagRight)){
        IMnpim_IMnpip_dE_woK0_n_Sm_side[0]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
        //if(!SigmawidePFlag){
          IMnpim_IMnpip_dE_woK0_n_Sm_side_cut[0]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
        //}
      }
      if(SigmaMsideHighFlag && (SigmaCrossMsideHighFlagLeft || SigmaCrossMsideHighFlagRight)){
        IMnpim_IMnpip_dE_woK0_n_Sm_side[1]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
        //if(!SigmawidePFlag){
          IMnpim_IMnpip_dE_woK0_n_Sm_side_cut[1]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
        //}
      }
      
      if(SigmaPFlag || SigmaMFlag){
        MMnpip_MMnpim_woK0_wSid_n->Fill(LVec_pim_n_miss.M(),LVec_pip_n_miss.M());
        dE_IMnpipi_woK0_wSid_n->Fill(LVec_pip_pim_n.M(),dE);
        Cosn_IMnpipi_woK0_wSid_n->Fill(LVec_pip_pim_n.M(),cos_n);
        IMnpim_IMnpip_dE_woK0_n_cut->Fill(LVec_pip_n.M(),LVec_pim_n.M());
        //MMnmiss_IMnpipi_woK0_wSid_n->Fill(LVec_pip_pim_n.M(), nmiss_mass);
        q_IMnpipi_woK0_wSid_n->Fill(LVec_pip_pim_n.M(),qkn.P());
        nmom_IMnpipi_woK0_wSid_n->Fill(LVec_pip_pim_n.M(),(*LVec_n).P());
        DCA_pip_beam->Fill( dca_pip_beam);
        DCA_pim_beam->Fill( dca_pim_beam );
        DCA_pip_pim->Fill(dca_pip_pim);
      }

      //if(SigmaPFlag && !SigmawideMFlag) {
      if(SigmaPFlag && !SigmaCrossFlagLeft && !SigmaCrossFlagRight) {
        IMnpim_IMnpip_dE_woK0_n_Sp->Fill(LVec_pip_n.M(),LVec_pim_n.M());
        q_IMnpipi_woK0_wSid_n_Sp->Fill(LVec_pip_pim_n.M(),qkn.P());
        q_IMnpipi_woK0_wSid_n_Sp_acc->Fill(LVec_pip_pim_n.M(),qkn.P());
      }
      if(SigmaPsideFlag && !SigmawideMFlag ){
     // if(SigmaPsideFlag && !SigmaMFlag) {
        IMnpim_IMnpip_dE_woK0_n_Sp_side_sum->Fill(LVec_pip_n.M(),LVec_pim_n.M());
        q_IMnpipi_woK0_wSid_n_Sp_side->Fill(LVec_pip_pim_n.M(),qkn.P());
      }
      //if(!SigmawidePFlag && SigmaMFlag){
      if( SigmaMFlag && !SigmaCrossFlagTop && !SigmaCrossFlagBottom){
        IMnpim_IMnpip_dE_woK0_n_Sm->Fill(LVec_pip_n.M(),LVec_pim_n.M());
        q_IMnpipi_woK0_wSid_n_Sm->Fill(LVec_pip_pim_n.M(),qkn.P());
        q_IMnpipi_woK0_wSid_n_Sm_acc->Fill(LVec_pip_pim_n.M(),qkn.P());
      }
      if(SigmaMsideFlag && !SigmawidePFlag ){
     // if(SigmaMsideFlag && !SigmaPFlag) {
        IMnpim_IMnpip_dE_woK0_n_Sm_side_sum->Fill(LVec_pip_n.M(),LVec_pim_n.M());
        q_IMnpipi_woK0_wSid_n_Sm_side->Fill(LVec_pip_pim_n.M(),qkn.P());
      }

    }
    if(K0rejectFlag && NBetaOK && NdEOK && MissNFlag){
      if(SidebandFlag){
        q_IMnpipi_woK0_wSid_n_side->Fill(LVec_pip_pim_n.M(),qkn.P());
      }
    }

    
    //including K0 
    if(NBetaOK && NdEOK && MissNFlag){
      q_pippim_n->Fill(LVec_pip_pim.M(),qkn.P());
      if(SigmaPFlag || SigmaMFlag){
        q_IMnpipi_wSid_n->Fill(LVec_pip_pim_n.M(),qkn.P());
      }
      if(SidebandFlag){
        q_IMnpipi_wSid_n_side->Fill(LVec_pip_pim_n.M(),qkn.P());
      }
    }

    // w/ kinfit
    if( (-1<kf_flag) && (pvalcut < pvalue) && K0rejectFlag){
      nmom_kin->Fill((*LVec_n).P());
      mnmom_kin->Fill(nmiss_mom);
      npipmom_kin->Fill(LVec_pip_n.P());
      npimmom_kin->Fill(LVec_pim_n.P());

      if( (kfSmmode_pvalue<kfSpmode_pvalue) && kfSpmode_status==0 ){
        dE_betainv_fid_kin[0]->Fill(1./NeutralBetaCDH_vtx[0],dE);
        if(NBetaOK){
          dE_MMom_fid_beta_woK0_kin[0]->Fill(LVec_n_miss_vtx[0].P(),dE);
          dE_MMass_fid_beta_woK0_kin[0]->Fill(LVec_n_miss_vtx[0].M(),dE);
        }
        if(NBetaOK && NdEOK){
          MMom_MMass_woK0_kin[0]->Fill(LVec_n_miss_vtx[0].M(),LVec_n_miss_vtx[0].P());
          IMnpim_IMnpip_dE_woK0_kin[0]->Fill(LVec_pip_n_vtx[0].M(),LVec_pim_n_vtx[0].M());
        }
        if(NBetaOK && NdEOK){
          MMnpip_MMnpim_woK0_kin[0]->Fill(LVec_pim_n_miss_vtx[0].M(),LVec_pip_n_miss_vtx[0].M());
          dE_IMnpipi_woK0_kin[0]->Fill(LVec_pip_pim_n_vtx[0].M(),dE);
          Cosn_IMnpipi_woK0_kin[0]->Fill(LVec_pip_pim_n_vtx[0].M(),cos_n);
          MMnmiss_IMnpipi_woK0_kin[0]->Fill(LVec_pip_pim_n_vtx[0].M(), nmiss_mass_vtx[0]);
          q_IMnpipi_woK0_kin[0]->Fill(LVec_pip_pim_n_vtx[0].M(),qkn_vtx[0].P());
          nmom_IMnpipi_woK0_kin[0]->Fill(LVec_pip_pim_n_vtx[0].M(),(*LVec_n_Sp).P());
          DCA_pip_beam_kin[0]->Fill( (*vtx_pip_beam-*vtx_pip_cdc).Mag());
          DCA_pim_beam_kin[0]->Fill( (*vtx_pim_beam-*vtx_pim_cdc).Mag());
          DCA_pip_pim_kin[0]->Fill( (*CA_pip-*CA_pim).Mag());
        }
      }else if( (kfSpmode_pvalue<kfSmmode_pvalue) && kfSmmode_status==0 ){//S- mode
        dE_betainv_fid_kin[1]->Fill(1./NeutralBetaCDH_vtx[1],dE);
        if(NBetaOK){
          dE_MMom_fid_beta_woK0_kin[1]->Fill(LVec_n_miss_vtx[1].P(),dE);
          dE_MMass_fid_beta_woK0_kin[1]->Fill(LVec_n_miss_vtx[1].M(),dE);
        }
        if(NBetaOK && NdEOK){
          MMom_MMass_woK0_kin[1]->Fill(LVec_n_miss_vtx[1].M(),LVec_n_miss_vtx[1].P());
          IMnpim_IMnpip_dE_woK0_kin[1]->Fill(LVec_pip_n_vtx[1].M(),LVec_pim_n_vtx[1].M());
          MMnpip_MMnpim_woK0_kin[1]->Fill(LVec_pim_n_miss_vtx[1].M(),LVec_pip_n_miss_vtx[1].M());
          dE_IMnpipi_woK0_kin[1]->Fill(LVec_pip_pim_n_vtx[1].M(),dE);
          Cosn_IMnpipi_woK0_kin[1]->Fill(LVec_pip_pim_n_vtx[1].M(),cos_n);
          MMnmiss_IMnpipi_woK0_kin[1]->Fill(LVec_pip_pim_n_vtx[1].M(), nmiss_mass_vtx[1]);
          q_IMnpipi_woK0_kin[1]->Fill(LVec_pip_pim_n_vtx[1].M(),qkn_vtx[1].P());
          nmom_IMnpipi_woK0_kin[1]->Fill(LVec_pip_pim_n_vtx[1].M(),(*LVec_n_Sm).P());
          DCA_pip_beam_kin[1]->Fill( (*vtx_pip_beam-*vtx_pip_cdc).Mag());
          DCA_pim_beam_kin[1]->Fill( (*vtx_pim_beam-*vtx_pim_cdc).Mag());
          DCA_pip_pim_kin[1]->Fill( (*CA_pip-*CA_pim).Mag());
        }
      }
    }//pval cut && Korejectflag
    
    if( (-1<kf_flag) && (pvalcut < pvalue) ){
      if( (kfSmmode_pvalue<kfSpmode_pvalue) && (kfSmmode_status==0) ){
        if(NBetaOK && NdEOK){ //&& MissNFlag)
          q_IMnpipi_kin[0]->Fill(LVec_pip_pim_n_vtx[0].M(),qkn_vtx[0].P());
        }
      }else if( (kfSpmode_pvalue<kfSmmode_pvalue) && (kfSmmode_status==0) ){//S- mode
        if(NBetaOK && NdEOK){ //&& MissNFlag)
          q_IMnpipi_kin[1]->Fill(LVec_pip_pim_n_vtx[1].M(),qkn_vtx[1].P());
        }
      }
    }   
	}//for ievt
   
  /*
  TCanvas *cKFpvalue_vs = new TCanvas(Form("cKFpvalue_vs"),"KFpvalue_vs");
  cKFpvalue_vs->cd();
  cKFpvalue_vs->SetLogy();
  TH1D *KFpvalue_vs_px = (TH1D*) KFpvalue_vs->ProjectionX();
  KFpvalue_vs_px->SetMinimum(1);
  KFpvalue_vs_px->SetXTitle("p-value");
  KFpvalue_vs_px->Draw();
  TH1D *KFpvalue_vs_py = (TH1D*) KFpvalue_vs->ProjectionY();
  KFpvalue_vs_py->SetLineColor(2);
  KFpvalue_vs_py->Draw("same");
  TLegend *legKFpvalue_vs = new TLegend(0.55,0.65,0.76,0.82);
  legKFpvalue_vs->AddEntry(KFpvalue_vs_px,"#Sigma^{+} mode");
  legKFpvalue_vs->AddEntry(KFpvalue_vs_py,"#Sigma^{-} mode");
  legKFpvalue_vs->SetFillColor(0);
  legKFpvalue_vs->Draw();

  int spbin = KFpvalue_vs_px->FindBin(pvalcut);
  int smbin = KFpvalue_vs_py->FindBin(pvalcut);
  std::cout << "Sp mode rejection ratio:" << KFpvalue_vs_px->Integral(0,spbin)/(KFpvalue_vs_px->Integral(0,201)) << std::endl;
  std::cout << "Sm mode rejection ratio:" << KFpvalue_vs_py->Integral(0,smbin)/(KFpvalue_vs_py->Integral(0,201)) << std::endl;

  //cumulative dist. of prob.
  TCanvas *cKFpvalue_vs_cum = new TCanvas(Form("cKFpvalue_vs_cum"),"KFpvalue_vs_cum");
  cKFpvalue_vs_cum->cd();
  TH1 *px_cum = KFpvalue_vs_px->GetCumulative();
  px_cum->Scale(1./KFpvalue_vs_px->GetEntries());
  px_cum->SetXTitle("p-value cut");
  px_cum->GetXaxis()->CenterTitle();
  px_cum->SetTitle("KF pvalue cumulative");
  px_cum->Draw();
  TH1 *py_cum = KFpvalue_vs_py->GetCumulative();
  py_cum->SetLineColor(2);
  py_cum->Scale(1./KFpvalue_vs_py->GetEntries());
  py_cum->Draw("same");
  TLegend *legKFpvalue_vs_cum = new TLegend(0.55,0.25,0.76,0.42);
  legKFpvalue_vs_cum->AddEntry(KFpvalue_vs_px,"#Sigma^{+} mode");
  legKFpvalue_vs_cum->AddEntry(KFpvalue_vs_py,"#Sigma^{-} mode");
  legKFpvalue_vs_cum->SetFillColor(0);
  legKFpvalue_vs_cum->Draw();

  TCanvas *cKFpvalue = new TCanvas(Form("cKFpvalue"),"KFpvalue");
  KFpvalue_vs->RebinX(5);
  KFpvalue_vs->RebinY(5);
  KFpvalue_vs->Draw("colz");
  */
    
  /*
  //chi2/ndf                                                                              
  TCanvas *cKFchi2ndf_vs = new TCanvas(Form("cKFchi2ndf_vs"),"KFchi2ndf_vs");             
  cKFchi2ndf_vs->cd();                                                                    
  cKFchi2ndf_vs->SetLogy();                                                               
  TH1D *KFchi2ndf_vs_px = (TH1D*) KFchi2ndf_vs->ProjectionX();                            
  //KFchi2ndf_vs_px->SetMinimum(1);                                                         
  KFchi2ndf_vs_px->GetXaxis()->SetTitle("chi2/ndf");                                      
  KFchi2ndf_vs_px->Draw();                                                                
  TH1D *KFchi2ndf_vs_py = (TH1D*) KFchi2ndf_vs->ProjectionY();                             
  KFchi2ndf_vs_py->SetLineColor(2);                                                       
  KFchi2ndf_vs_py->Draw("same");                                                          
  TLegend *legKFchi2ndf_vs = new TLegend(0.55,0.65,0.76,0.82);                            
  legKFchi2ndf_vs->AddEntry(KFchi2ndf_vs_px,"#Sigma^{+} mode");                                        
  legKFchi2ndf_vs->AddEntry(KFchi2ndf_vs_py,"#Sigma^{-} mode");                                        
  legKFchi2ndf_vs->SetFillColor(0);                                                       
  legKFchi2ndf_vs->Draw();                                                                

  //int spbin = KFchi2ndf_vs_px->FindBin(chi2cut);                                        
  //int smbin = KFchi2ndf_vs_py->FindBin(chi2cut);                                        
  //std::cout << "Sp mode rejection ratio:" << KFchi2ndf_vs_px->Integral(0,spbin)/(KFchi2n
  //std::cout << "Sm mode rejection ratio:" << KFchi2ndf_vs_py->Integral(0,smbin)/(KFchi2n

  //cumulative dist. of prob.                                                             
  TCanvas *cKFchi2ndf_vs_cum = new TCanvas(Form("cKFchi2ndf_vs_cum"),"KFchi2ndf_vs_cum"); 
  cKFchi2ndf_vs_cum->cd();                                                                
  TH1 *KFchi2ndf_vs_cum_px = KFchi2ndf_vs_px->GetCumulative();                            
  KFchi2ndf_vs_cum_px->Scale(1./KFchi2ndf_vs_px->GetEntries());                                        
  KFchi2ndf_vs_cum_px->SetXTitle("chi2ndf cut");                                          
  KFchi2ndf_vs_cum_px->GetXaxis()->CenterTitle();                                         
  KFchi2ndf_vs_cum_px->SetTitle("KF chi2ndf cumulative");                                 
  KFchi2ndf_vs_cum_px->Draw();                                                            
  TH1 *KFchi2ndf_vs_cum_py = KFchi2ndf_vs_py->GetCumulative();                            
  KFchi2ndf_vs_cum_py->SetLineColor(2);                                                   
  KFchi2ndf_vs_cum_py->Scale(1./KFchi2ndf_vs_py->GetEntries());                                        
  KFchi2ndf_vs_cum_py->Draw("same");                                                      
  TLegend *legKFchi2ndf_vs_cum = new TLegend(0.55,0.25,0.76,0.42);                        
  legKFchi2ndf_vs_cum->AddEntry(KFchi2ndf_vs_cum_px,"#Sigma^{+} mode");                                    
  legKFchi2ndf_vs_cum->AddEntry(KFchi2ndf_vs_cum_py,"#Sigma^{-} mode");                                    
  legKFchi2ndf_vs_cum->SetFillColor(0);                                                   
  legKFchi2ndf_vs_cum->Draw();                                                            

  TCanvas *cKFchi2ndf = new TCanvas(Form("cKFchi2ndf"),"KFchi2ndf");                      
  KFchi2ndf_vs->RebinX(5);                                                                
  KFchi2ndf_vs->RebinY(5);                                                                
  KFchi2ndf_vs->Draw("colz");                                                             

  TCanvas *cnmom = new TCanvas("cnmom","nmom");
  cnmom->cd();
  nmom->Draw("");
  gPad->SetLogy(0);

  TCanvas *cmnmom = new TCanvas("cmnmom","mnmom");
  cmnmom->cd();
  mnmom->Draw("");

  TCanvas *cnpipmom = new TCanvas("cnpipmom","npipmom");
  cnpipmom->cd();
  npipmom->Draw("");
  TCanvas *cnpimmom = new TCanvas("cnpimmom","npimmom");
  cnpimmom->cd();
  npimmom->Draw("");
  */
  TCanvas *cMMnmiss_IMnpip_dE_woK0 = new TCanvas("cMMnmiss_IMnpip_dE_woK0","MMnmiss_IMnpip_dE_woK0");
  MMnmiss_IMnpip_dE_woK0->GetXaxis()->SetRangeUser(1.05,1.5);
  MMnmiss_IMnpip_dE_woK0->GetYaxis()->SetRangeUser(0.6,1.5);
  MMnmiss_IMnpip_dE_woK0->Draw("colz");

  TCanvas *cMMnmiss_IMnpim_dE_woK0 = new TCanvas("cMMnmiss_IMnpim_dE_woK0","MMnmiss_IMnpim_dE_woK0");
  MMnmiss_IMnpim_dE_woK0->GetXaxis()->SetRangeUser(1.05,1.5);
  MMnmiss_IMnpim_dE_woK0->GetYaxis()->SetRangeUser(0.6,1.5);
  MMnmiss_IMnpim_dE_woK0->Draw("colz");
  
  /*
  TCanvas *cq_IMnpipi_wSid_n_px = new TCanvas("cq_IMnpipi_wSid_n_px","q_IMnpipi_wSid_n_px"); 
  cq_IMnpipi_wSid_n_px ->cd();
  TH1D *q_IMnpipi_wSid_n_px = q_IMnpipi_wSid_n->ProjectionX();
  TH1D *q_IMnpipi_kin_0_n_px = q_IMnpipi_kin[0]->ProjectionX();
  TH1D *q_IMnpipi_kin_1_n_px = q_IMnpipi_kin[1]->ProjectionX();
  q_IMnpipi_wSid_n_px->Draw("HE");
  q_IMnpipi_kin_0_n_px->SetLineColor(2);
  q_IMnpipi_kin_0_n_px->Draw("HEsame");
  q_IMnpipi_kin_1_n_px->SetLineColor(3);
  q_IMnpipi_kin_1_n_px->Draw("HEsame");
  TH1D *pxSum1 = (TH1D*) q_IMnpipi_kin_0_n_px->Clone();
  pxSum1->Add(q_IMnpipi_kin_1_n_px);
  pxSum1->SetLineColor(4);
  pxSum1->Draw("HEsame");
  */
  
  //TCanvas *cq_IMnpipi_woK0_wSid_n_px = new TCanvas("cq_IMnpipi_woK0_wSid_n_px","q_IMnpipi_woK0_wSid_n_px"); 
  //cq_IMnpipi_woK0_wSid_n_px ->cd();
  TH1D *q_IMnpipi_woK0_wSid_n_px = q_IMnpipi_woK0_wSid_n->ProjectionX();
  //TH1D *q_IMnpipi_woK0_kin_0_n_px = q_IMnpipi_woK0_kin[0]->ProjectionX();
  //TH1D *q_IMnpipi_woK0_kin_1_n_px = q_IMnpipi_woK0_kin[1]->ProjectionX();
  //q_IMnpipi_woK0_wSid_n_px->Draw("HE");
  //q_IMnpipi_woK0_kin_0_n_px->SetLineColor(2);
  //q_IMnpipi_woK0_kin_0_n_px->Draw("HEsame");
  //q_IMnpipi_woK0_kin_1_n_px->SetLineColor(3);
  //q_IMnpipi_woK0_kin_1_n_px->Draw("HEsame");
  //TH1D *pxSum = (TH1D*) q_IMnpipi_woK0_kin_0_n_px->Clone();
  //pxSum->Add(q_IMnpipi_woK0_kin_1_n_px);
  //pxSum->SetLineColor(4);
  //pxSum->Draw("HEsame");
  

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
  q_IMnpipi_wSid_n_px1->Draw("HEsame");
  
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
  TCanvas *csub = new TCanvas("csub","csub");
  TH1D *q_IMnpipi_woK0_wSid_n_px_sub = (TH1D*) q_IMnpipi_woK0_wSid_n_px->Clone();
  q_IMnpipi_woK0_wSid_n_px_sub->Add(q_IMnpipi_woK0_wSid_n_side_px,-1);
  q_IMnpipi_woK0_wSid_n_px_sub->SetMinimum(0);
  q_IMnpipi_woK0_wSid_n_px_sub->Draw("EH"); 
  //fkp->Draw("same");
  
  /*
  TCanvas *cq_IMnpipi_woK0_wSid_n_py = new TCanvas("cq_IMnpipi_woK0_wSid_n_py","q_IMnpipi_woK0_wSid_n_py"); 
  cq_IMnpipi_woK0_wSid_n_py->cd();
  TH1D *q_IMnpipi_wSid_n_py = q_IMnpipi_wSid_n->ProjectionY();
  TH1D *q_IMnpipi_woK0_wSid_n_py = q_IMnpipi_woK0_wSid_n->ProjectionY();
  TH1D *q_IMnpipi_woK0_kin_0_py = q_IMnpipi_woK0_kin[0]->ProjectionY();
  TH1D *q_IMnpipi_woK0_kin_1_py = q_IMnpipi_woK0_kin[1]->ProjectionY();
  q_IMnpipi_woK0_wSid_n_py->Draw("HE");
  q_IMnpipi_wSid_n_py->SetLineColor(5);
  q_IMnpipi_wSid_n_py->Draw("HEsame");
  q_IMnpipi_woK0_kin_1_py->SetLineColor(2);
  q_IMnpipi_woK0_kin_1_py->Draw("HEsame");
  q_IMnpipi_woK0_kin_0_py->SetLineColor(3);
  q_IMnpipi_woK0_kin_0_py->Draw("HEsame");
  TH1D *pySum = (TH1D*) q_IMnpipi_woK0_kin_0_py->Clone();
  pySum->Add(q_IMnpipi_woK0_kin_1_py);
  pySum->SetLineColor(4);
  pySum->Draw("HEsame");
  */


  TCanvas *cq_IMnpipi_woK0_wSid_n_px_SpSm = new TCanvas("cq_IMnpipi_woK0_wSid_n_px_SpSm","q_IMnpipi_woK0_wSid_n_px_SpSm"); 
  TH1D* q_IMnpipi_woK0_wSid_n_Sp_px = (TH1D*) q_IMnpipi_woK0_wSid_n_Sp->ProjectionX();
  TH1D* q_IMnpipi_woK0_wSid_n_Sm_px = (TH1D*) q_IMnpipi_woK0_wSid_n_Sm->ProjectionX();
  q_IMnpipi_woK0_wSid_n_Sp_px->SetLineColor(2);
  q_IMnpipi_woK0_wSid_n_Sm_px->SetLineColor(3);
  q_IMnpipi_woK0_wSid_n_px->Draw("EH");
  //q_IMnpipi_wSid_n_px1->Draw("HEsame");
  q_IMnpipi_woK0_wSid_n_Sp_px->Draw("HEsame");
  q_IMnpipi_woK0_wSid_n_Sm_px->Draw("HEsame");
  TH1D* q_IMnpipi_woK0_wSid_n_wocross = (TH1D*)q_IMnpipi_woK0_wSid_n_Sp_px->Clone();
  q_IMnpipi_woK0_wSid_n_wocross->Add(q_IMnpipi_woK0_wSid_n_Sm_px,1);
  q_IMnpipi_woK0_wSid_n_wocross->SetLineColor(4);
  q_IMnpipi_woK0_wSid_n_wocross->Draw("HEsame");

  TCanvas *cq_IMnpipi_woK0_wSid_n_px_Sp = new TCanvas("cq_IMnpipi_woK0_wSid_n_px_Sp","q_IMnpipi_woK0_wSid_n_px_Sp"); 
  //q_IMnpipi_woK0_wSid_n_Sp_px->SetMaximum(q_IMnpipi_woK0_wSid_n_px->GetMaximum());
  q_IMnpipi_woK0_wSid_n_Sp_px->Draw("HE");
  TH1D* q_IMnpipi_woK0_wSid_n_Sp_side_px = (TH1D*) q_IMnpipi_woK0_wSid_n_Sp_side->ProjectionX();
  q_IMnpipi_woK0_wSid_n_Sp_side_px->SetLineColor(5);
  q_IMnpipi_woK0_wSid_n_Sp_side_px->Draw("HEsame");
  
  TCanvas *csubSp = new TCanvas("csubSp","csubSp");
  csubSp->cd();
  TH2F* q_IMnpipi_woK0_wSid_n_Sp_sub = (TH2F*)q_IMnpipi_woK0_wSid_n_Sp->Clone("Sp_sub");
  q_IMnpipi_woK0_wSid_n_Sp_sub->Add(q_IMnpipi_woK0_wSid_n_Sp_side,-1);
  q_IMnpipi_woK0_wSid_n_Sp_sub->SetMinimum(0);
  q_IMnpipi_woK0_wSid_n_Sp_sub->Draw("colz");



  TCanvas *csubSm = new TCanvas("csubSm","csubSm");
  csubSm->cd();
  TH2F* q_IMnpipi_woK0_wSid_n_Sm_sub = (TH2F*)q_IMnpipi_woK0_wSid_n_Sm->Clone("Sm_sub");
  q_IMnpipi_woK0_wSid_n_Sm_sub->Add(q_IMnpipi_woK0_wSid_n_Sm_side,-1);
  q_IMnpipi_woK0_wSid_n_Sm_sub->SetMinimum(0);
  q_IMnpipi_woK0_wSid_n_Sm_sub->Draw("colz");
  
  
  TCanvas *csubSp_px = new TCanvas("csubSp_px","csubSp_px");
  csubSp_px->cd();
  TH1D* q_IMnpipi_woK0_wSid_n_Sp_sub_px = q_IMnpipi_woK0_wSid_n_Sp_sub->ProjectionX();
  q_IMnpipi_woK0_wSid_n_Sp_sub_px->SetLineColor(2);
  q_IMnpipi_woK0_wSid_n_Sp_sub_px->SetMinimum(0);
  q_IMnpipi_woK0_wSid_n_Sp_sub_px->Draw("E");

  //TCanvas *csubSm_px = new TCanvas("csubSm_px","csubSm_px");
  //csubSm_px->cd();
  TH1D* q_IMnpipi_woK0_wSid_n_Sm_sub_px = q_IMnpipi_woK0_wSid_n_Sm_sub->ProjectionX();
  q_IMnpipi_woK0_wSid_n_Sm_sub_px->SetLineColor(3);
  q_IMnpipi_woK0_wSid_n_Sm_sub_px->Draw("Esame");

  TH1D* SpSm_sub_px = (TH1D*)q_IMnpipi_woK0_wSid_n_Sp_sub_px->Clone("SpSm_sub_px");
  SpSm_sub_px->Add(q_IMnpipi_woK0_wSid_n_Sm_sub_px);
  SpSm_sub_px->SetLineColor(4);
  SpSm_sub_px->Draw("Esame");


  TCanvas *cq_IMnpipi_woK0_wSid_n_px_Sm = new TCanvas("cq_IMnpipi_woK0_wSid_n_px_Sm","q_IMnpipi_woK0_wSid_n_px_Sm"); 
  //q_IMnpipi_woK0_wSid_n_Sm_px->SetMaximum(q_IMnpipi_woK0_wSid_n_px->GetMaximum());
  q_IMnpipi_woK0_wSid_n_Sm_px->Draw("EH");
  TH1D* q_IMnpipi_woK0_wSid_n_Sm_side_px = (TH1D*) q_IMnpipi_woK0_wSid_n_Sm_side->ProjectionX();
  q_IMnpipi_woK0_wSid_n_Sm_side_px->SetLineColor(5);
  q_IMnpipi_woK0_wSid_n_Sm_side_px->Draw("HEsame");

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

   TCanvas *cq_IMnpipi_woK0_wSid_n_side = new TCanvas("cq_IMnpipi_woK0_wSid_n_side","q_IMnpipi_woK0_wSid_n_side");
  cq_IMnpipi_woK0_wSid_n_side->cd();
  //cq_IMnpipi_woK0_wSid_n->cd(2);
  q_IMnpipi_woK0_wSid_n_side->SetMaximum(q_IMnpipi_woK0_wSid_n->GetMaximum());
  q_IMnpipi_woK0_wSid_n_side->Draw("colz");
  
  TCanvas *cq_IMnpipi_woK0_wSid_n_Sp = new TCanvas("cq_IMnpipi_woK0_wSid_n_Sp","q_IMnpipi_woK0_wSid_n_Sp");
  cq_IMnpipi_woK0_wSid_n_Sp->cd();
  q_IMnpipi_woK0_wSid_n_Sp->Draw("colz");

  
  TCanvas *cq_IMnpipi_woK0_wSid_n_Sm = new TCanvas("cq_IMnpipi_woK0_wSid_n_Sm","q_IMnpipi_woK0_wSid_n_Sm");
  cq_IMnpipi_woK0_wSid_n_Sm->cd();
//  q_IMnpipi_woK0_wSid_n_Sm->SetMaximum(q_IMnpipi_woK0_wSid_n_Sp->GetMaximum());
  q_IMnpipi_woK0_wSid_n_Sm->Draw("colz");
  
  TCanvas *cq_IMnpipi_woK0_wSid_n_Sp_side = new TCanvas("cq_IMnpipi_woK0_wSid_n_Sp_side","q_IMnpipi_woK0_wSid_n_Sp_side");
  cq_IMnpipi_woK0_wSid_n_Sp_side->cd();
  q_IMnpipi_woK0_wSid_n_Sp_side->Draw("colz");


  
  TCanvas *cq_IMnpipi_woK0_wSid_n_Sm_side = new TCanvas("cq_IMnpipi_woK0_wSid_n_Sm_side","q_IMnpipi_woK0_wSid_n_Sm_side");
  cq_IMnpipi_woK0_wSid_n_Sm_side->cd();
  q_IMnpipi_woK0_wSid_n_Sm_side->Draw("colz");


  TCanvas *cqsub2 = new TCanvas("cqsub2","cqsub2");
  cqsub2->cd();
  TH2F *q_IMnpipi_woK0_wSid_n_sub = (TH2F*) q_IMnpipi_woK0_wSid_n->Clone("sub");
  q_IMnpipi_woK0_wSid_n_sub->Add(q_IMnpipi_woK0_wSid_n_side,-1);
  q_IMnpipi_woK0_wSid_n_sub->SetMinimum(0);
  q_IMnpipi_woK0_wSid_n_sub->Draw("colz");
  fkp->Draw("same");
  TFile *fnu = new TFile("NumericalRootFinder.root");
  fnu->cd();
  TMultiGraph *mg = (TMultiGraph*)fnu->Get("mg"); 
  mg->Draw("c");
  f->cd();

  //TCanvas *cq_IMnpipi_woK0_kin_sum = new TCanvas("cq_IMnpipi_woK0_kin_sum","q_IMnpipi_woK0_kin_sum");
  //cq_IMnpipi_woK0_kin_sum->cd();
  //TH2F *q_IMnpipi_woK0_kin_sum = (TH2F*) q_IMnpipi_woK0_kin[0]->Clone();  
  //q_IMnpipi_woK0_kin_sum->Add(q_IMnpipi_woK0_kin[1]);
  //q_IMnpipi_woK0_kin_sum->SetTitle("q_IMnpipi_woK0_kin_sum");
  //q_IMnpipi_woK0_kin_sum->GetZaxis()->SetRangeUser(0,80);
  //q_IMnpipi_woK0_kin_sum->Draw("colz");
  //TFile *fnu = new TFile("NumericalRootFinder_Spmode.root");
  //TFile *fnu = new TFile("NumericalRootFinder.root");
  fnu->cd();
  //TMultiGraph *mg = (TMultiGraph*)fnu->Get("mg"); 
  mg->Draw("c");
  f->cd();

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

  TCanvas *cMMnmiss_IMnpipi_woK0_wSid_n_px = new TCanvas("cMMnmiss_IMnpipi_woK0_wSid_n_px","MMnmiss_IMnpipi_woK0_wSid_n_px"); 
  cMMnmiss_IMnpipi_woK0_wSid_n_px->cd();
  MMnmiss_IMnpipi_woK0_wSid_n->ProjectionX()->Draw("");
  TCanvas *cMMnmiss_IMnpipi_woK0_wSid_n_py = new TCanvas("cMMnmiss_IMnpipi_woK0_wSid_n_py","MMnmiss_IMnpipi_woK0_wSid_n_py"); 
  cMMnmiss_IMnpipi_woK0_wSid_n_py->cd();
  MMnmiss_IMnpipi_woK0_wSid_n->ProjectionY()->Draw("");

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

  TCanvas *cdE_betainv_fid_kin[2];//           
  TCanvas *cdE_MMom_fid_beta_woK0_kin[2];      
  TCanvas *cdE_MMass_fid_beta_woK0_kin[2];     
  TCanvas *cMMom_MMass_woK0_kin[2];
  TCanvas *cMMom_MMass_woK0_kin_px[2];
  TCanvas *cMMom_MMass_woK0_kin_py[2];
  TH1D *MMom_MMass_woK0_kin_px[2];
  TH1D *MMom_MMass_woK0_kin_py[2];
  TCanvas *cIMnpim_IMnpip_dE_woK0_kin[2]; 
  TCanvas *cIMnpim_IMnpip_dE_woK0_kin_px[2];
  TH1D *IMnpim_IMnpip_dE_woK0_kin_px[2];
  TCanvas *cIMnpim_IMnpip_dE_woK0_kin_py[2]; 
  TH1D *IMnpim_IMnpip_dE_woK0_kin_py[2];

  TCanvas *cMMnpip_MMnpim_woK0_kin[2];  
  TCanvas *cdE_IMnpipi_woK0_kin[2];     
  TCanvas *cCosn_IMnpipi_n_kin[2];   
  TCanvas *cMMnmiss_IMnpipi_woK0_kin[2];
  TCanvas *cq_IMnpipi_woK0_kin[2];      
  TCanvas *cnmom_IMnpipi_woK0_kin[2];     
  
  
  for(int imode=0;imode<2;imode++){
    /*
    cdE_betainv_fid_kin[imode] = new TCanvas(Form("cdE_betainv_fid_kin_%s",smode[imode]), Form("dE_betainv_fid_kin_%s",smode[imode]));
    cdE_betainv_fid_kin[imode]->cd();
    dE_betainv_fid_kin[imode]->Draw("colz");


    cdE_MMom_fid_beta_woK0_kin[imode] = new TCanvas(Form("cdE_MMom_fid_beta_woK0_kin_%s",smode[imode]),Form("dE_MMom_fid_beta_woK0_kin_%s",smode[imode]));
    cdE_MMom_fid_beta_woK0_kin[imode]->cd();
    dE_MMom_fid_beta_woK0_kin[imode]->Draw("colz");

    cdE_MMass_fid_beta_woK0_kin[imode] = new TCanvas(Form("cdE_MMass_fid_beta_woK0_%s_kin",smode[imode]),Form("dE_MMass_fid_beta_woK0_%s_kin",smode[imode]));
    cdE_MMass_fid_beta_woK0_kin[imode]->cd(); 
    dE_MMass_fid_beta_woK0_kin[imode]->Draw("colz");

    cMMom_MMass_woK0_kin[imode] = new TCanvas(Form("cMMom_MMass_woK0_kin_%s",smode[imode]),Form("MMom_MMass_woK0_kin_%s",smode[imode]));
    cMMom_MMass_woK0_kin[imode]->cd();
    MMom_MMass_woK0_kin[imode]->Draw("colz");
    */
    //cMMom_MMass_woK0_kin_px[imode] = new TCanvas(Form("cMMom_MMass_woK0_kin_px_%s",smode[imode]),Form("MMom_MMass_woK0_kin_px_%s",smode[imode]));
    //cMMom_MMass_woK0_kin_px[imode]->cd();
    MMom_MMass_woK0_kin_px[imode] = MMom_MMass_woK0_kin[imode]->ProjectionX();
    //MMom_MMass_woK0_kin_px[imode]->Draw();
    
    //cMMom_MMass_woK0_kin_py[imode] = new TCanvas(Form("cMMom_MMass_woK0_kin_py_%s",smode[imode]),Form("MMom_MMass_woK0_kin_py_%s",smode[imode]));
    //cMMom_MMass_woK0_kin_py[imode]->cd();
    MMom_MMass_woK0_kin_py[imode] = MMom_MMass_woK0_kin[imode]->ProjectionY();
    //MMom_MMass_woK0_kin_py[imode]->Draw();
    
    
    //cIMnpim_IMnpip_dE_woK0_kin[imode] = new TCanvas(Form("cIMnpim_IMnpip_dE_woK0_kin_%s",smode[imode]),Form("IMnpim_IMnpip_dE_woK0_kin_%s",smode[imode]));
    //cIMnpim_IMnpip_dE_woK0_kin[imode]->cd();
    //IMnpim_IMnpip_dE_woK0_kin[imode]->Draw("colz");
    
    //cIMnpim_IMnpip_dE_woK0_kin_px[imode] = new TCanvas(Form("cIMnpim_IMnpip_dE_woK0_kin_px_%s",smode[imode]),Form("IMnpim_IMnpip_dE_woK0_kin_px_%s",smode[imode]));
    //cIMnpim_IMnpip_dE_woK0_kin_px[imode]->cd();
    IMnpim_IMnpip_dE_woK0_kin_px[imode] = IMnpim_IMnpip_dE_woK0_kin[imode]->ProjectionX();
    //IMnpim_IMnpip_dE_woK0_kin_px[imode]->Draw();
    
    //cIMnpim_IMnpip_dE_woK0_kin_py[imode] = new TCanvas(Form("cIMnpim_IMnpip_dE_woK0_kin_py_%s",smode[imode]),Form("IMnpim_IMnpip_dE_woK0_kin_py_%s",smode[imode]));
    //cIMnpim_IMnpip_dE_woK0_kin_py[imode]->cd();
    IMnpim_IMnpip_dE_woK0_kin_py[imode] = IMnpim_IMnpip_dE_woK0_kin[imode]->ProjectionY();
    //IMnpim_IMnpip_dE_woK0_kin_py[imode]->Draw();
    
    //cMMnpip_MMnpim_woK0_kin[imode] = new TCanvas(Form("cMMnpip_MMnpim_woK0_kin_%s",smode[imode]),Form("MMnpip_MMnpim_woK0_kin_%s",smode[imode]));
    //cMMnpip_MMnpim_woK0_kin[imode]->cd();
    //MMnpip_MMnpim_woK0_kin[imode]->Draw("colz");

    //cdE_IMnpipi_woK0_kin[imode] = new TCanvas(Form("cdE_IMnpipi_woK0_kin_%s",smode[imode]),Form("dE_IMnpipi_woK0_kin_%s",smode[imode]));
    //cdE_IMnpipi_woK0_kin[imode]->cd();
    //dE_IMnpipi_woK0_kin[imode]->Draw("colz");
    
    //cCosn_IMnpipi_n_kin[imode] = new TCanvas(Form("cCosn_IMnpipi_n_kin_%s",smode[imode]),Form("Cosn_IMnpipi_n_kin_%s",smode[imode]));
    //cCosn_IMnpipi_n_kin[imode]->cd();
    //Cosn_IMnpipi_woK0_kin[imode]->Draw("colz");

    //cMMnmiss_IMnpipi_woK0_kin[imode] = new TCanvas(Form("cMMnmiss_IMnpipi_woK0_kin_%s",smode[imode]),Form("MMnmiss_IMnpipi_woK0_kin_%s",smode[imode]));
    //cMMnmiss_IMnpipi_woK0_kin[imode]->cd();
    //MMnmiss_IMnpipi_woK0_kin[imode]->Draw("colz");
    
    //cq_IMnpipi_woK0_kin[imode] = new TCanvas(Form("cq_IMnpipi_woK0_kin_%s",smode[imode]),Form("q_IMnpipi_woK0_kin_%s",smode[imode]));
    //cq_IMnpipi_woK0_kin[imode]->cd();
    //q_IMnpipi_woK0_kin[imode]->Draw("colz");
    if(imode==0){
      TFile *fnup = new TFile("NumericalRootFinder_Spmode.root");
      fnup->cd();
      TMultiGraph *mgp = (TMultiGraph*)fnup->Get("mg"); 
      mgp->Draw("c");
      f->cd();
    }else{
      TFile *fnum = new TFile("NumericalRootFinder_Smmode.root");
      fnum->cd();
      TMultiGraph *mgm = (TMultiGraph*)fnum->Get("mg"); 
      mgm->Draw("c");
      f->cd();
    }
   
    if(imode==1)q_IMnpipi_woK0_kin[imode]->SetMaximum(q_IMnpipi_woK0_kin[0]->GetMaximum());
    
    //cnmom_IMnpipi_woK0_kin[imode] = new TCanvas(Form("cnmom_IMnpipi_woK0_kin_%s",smode[imode]),Form("nmom_IMnpipi_woK0_kin_%s",smode[imode]));
    //cnmom_IMnpipi_woK0_kin[imode]->cd();
    //nmom_IMnpipi_woK0_kin[imode]->Draw("colz");
  }
  

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
  TH1D *MMom_MMass_woK0_kin_px_clone[3];
  MMom_MMass_woK0_kin_px_clone[0] = (TH1D*)MMom_MMass_woK0_kin_px[0]->Clone();
  MMom_MMass_woK0_kin_px_clone[0]->SetLineColor(2);
  //MMom_MMass_woK0_kin_px_clone[0]->Draw("same");
  MMom_MMass_woK0_kin_px_clone[1] = (TH1D*)MMom_MMass_woK0_kin_px[1]->Clone();
  MMom_MMass_woK0_kin_px_clone[1]->SetLineColor(3);
  //MMom_MMass_woK0_kin_px_clone[1]->Draw("same");
  MMom_MMass_woK0_kin_px_clone[2] = (TH1D*)MMom_MMass_woK0_kin_px[0]->Clone();
  MMom_MMass_woK0_kin_px_clone[2]->Add(MMom_MMass_woK0_kin_px_clone[1]);
  MMom_MMass_woK0_kin_px_clone[2]->SetLineColor(4);
  //MMom_MMass_woK0_kin_px_clone[2]->Draw("same");
  std::cout << "miss n mean " << MMom_MMass_woK0_kin_px_clone[2]->GetMean() << std::endl;
  std::cout << "miss n rms " << MMom_MMass_woK0_kin_px_clone[2]->GetRMS() << std::endl;

  //TCanvas *cMMom_MMass_woK0_py = new TCanvas("cMMom_MMass_woK0_py","MMom_MMass_woK0_py");
  TH1D *MMom_MMass_woK0_py = MMom_MMass_woK0->ProjectionY();
  TH1D *MMom_MMass_woK0_wSid_py = MMom_MMass_woK0_wSid->ProjectionY();
  //MMom_MMass_woK0_py->Draw();
  TH1D *MMom_MMass_woK0_wSid_py_clone = (TH1D*)MMom_MMass_woK0_wSid_py->Clone();
  //MMom_MMass_woK0_wSid_py_clone->SetLineColor(4);
  //cMMom_MMass_woK0_py->cd();
  MMom_MMass_woK0_wSid_py_clone->SetTitle("Missing Mom. d(K^{-},#pi^{+}#pi^{-}n)\"X\"");
  //MMom_MMass_woK0_wSid_py_clone->Draw("");
  TH1D* MMom_MMass_woK0_kin_py_clone[3];
  MMom_MMass_woK0_kin_py_clone[0] = (TH1D*) MMom_MMass_woK0_kin_py[0]->Clone();
  MMom_MMass_woK0_kin_py_clone[0]->SetLineColor(2);
  //MMom_MMass_woK0_kin_py_clone[0]->Draw("same");
  MMom_MMass_woK0_kin_py_clone[1] = (TH1D*) MMom_MMass_woK0_kin_py[1]->Clone();
  MMom_MMass_woK0_kin_py_clone[1]->SetLineColor(3);
  //MMom_MMass_woK0_kin_py_clone[1]->Draw("same");
  MMom_MMass_woK0_kin_py_clone[2] = (TH1D*) MMom_MMass_woK0_kin_py[0]->Clone(); 
  MMom_MMass_woK0_kin_py_clone[2]->Add(MMom_MMass_woK0_kin_py_clone[1]);
  MMom_MMass_woK0_kin_py_clone[2]->SetLineColor(4);
  //MMom_MMass_woK0_kin_py_clone[2]->Draw("same");
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
  //TH1D * IMnpim_IMnpip_dE_woK0_n_px_sum =(TH1D*) IMnpim_IMnpip_dE_woK0_kin_px[0]->Clone();
  
  //TCanvas *cIMnpim_IMnpip_dE_woK0_py = new TCanvas("cIMnpim_IMnpip_dE_woK0_py","IMnpim_IMnpip_dE_woK0_py");
  //cIMnpim_IMnpip_dE_woK0_py->cd();
  //TH1D *IMnpim_IMnpip_dE_woK0_py = IMnpim_IMnpip_dE_woK0->ProjectionY();
  //IMnpim_IMnpip_dE_woK0_py->Draw("HE");

  TCanvas *cIMnpim_IMnpip_dE_woK0_n = new TCanvas("cIMnpim_IMnpip_dE_woK0_n","IMnpim_IMnpip_dE_woK0_n",800,800);
  cIMnpim_IMnpip_dE_woK0_n->cd();
  IMnpim_IMnpip_dE_woK0_n->GetXaxis()->SetRangeUser(0,1.7);
  IMnpim_IMnpip_dE_woK0_n->GetYaxis()->SetRangeUser(0,1.7);
  IMnpim_IMnpip_dE_woK0_n->Draw("colz");
  
  TCanvas *cIMnpim_IMnpip_dE_woK0_n_cut = new TCanvas("cIMnpim_IMnpip_dE_woK0_n_cut","IMnpim_IMnpip_dE_woK0_n_cut",800,800);
  cIMnpim_IMnpip_dE_woK0_n_cut->cd();
  IMnpim_IMnpip_dE_woK0_n_cut->SetMaximum(IMnpim_IMnpip_dE_woK0_n->GetMaximum());
  IMnpim_IMnpip_dE_woK0_n_cut->GetXaxis()->SetRangeUser(0,1.7);
  IMnpim_IMnpip_dE_woK0_n_cut->GetYaxis()->SetRangeUser(0,1.7);
  IMnpim_IMnpip_dE_woK0_n_cut->Draw("colz");
  
  //TCanvas *cIMnpim_IMnpip_dE_woK0_n_side = new TCanvas("cIMnpim_IMnpip_dE_woK0_n_side","IMnpim_IMnpip_dE_woK0_n_side",800,800);
 // cIMnpim_IMnpip_dE_woK0_n_side->cd();
  IMnpim_IMnpip_dE_woK0_n_side->SetMaximum(IMnpim_IMnpip_dE_woK0_n->GetMaximum());
  IMnpim_IMnpip_dE_woK0_n_side->GetXaxis()->SetRangeUser(0,1.7);
  IMnpim_IMnpip_dE_woK0_n_side->GetYaxis()->SetRangeUser(0,1.7);
  IMnpim_IMnpip_dE_woK0_n_side->Draw("colsame");


  TCanvas *cIMnpim_IMnpip_dE_woK0_n_SpSm = new TCanvas("cIMnpim_IMnpip_dE_woK0_n_SpSm","IMnpim_IMnpip_dE_woK0_n_SpSM",1000,1000);
  IMnpim_IMnpip_dE_woK0_n_Sp->GetXaxis()->SetRangeUser(0,1.7);
  IMnpim_IMnpip_dE_woK0_n_Sp->GetYaxis()->SetRangeUser(0,1.7);
  IMnpim_IMnpip_dE_woK0_n_Sp->Draw("colz");
  IMnpim_IMnpip_dE_woK0_n_Sm->Draw("colsame");


  TCanvas *cIMnpim_IMnpip_dE_woK0_n_Sp = new TCanvas("cIMnpim_IMnpip_dE_woK0_n_Sp","IMnpim_IMnpip_dE_woK0_n_Sp");
  IMnpim_IMnpip_dE_woK0_n_Sp->Draw("colz");
  IMnpim_IMnpip_dE_woK0_n_Sp_side_cut[0]->Draw("colsame");
  IMnpim_IMnpip_dE_woK0_n_Sp_side_cut[1]->Draw("colsame");
  
  TCanvas *cIMnpim_IMnpip_dE_woK0_n_Sm = new TCanvas("cIMnpim_IMnpip_dE_woK0_n_Sm","IMnpim_IMnpip_dE_woK0_n_Sm");
  IMnpim_IMnpip_dE_woK0_n_Sm->GetXaxis()->SetRangeUser(0,1.7);
  IMnpim_IMnpip_dE_woK0_n_Sm->GetYaxis()->SetRangeUser(0,1.7);
  IMnpim_IMnpip_dE_woK0_n_Sm->Draw("colz");
  IMnpim_IMnpip_dE_woK0_n_Sm_side_cut[0]->Draw("colsame");
  IMnpim_IMnpip_dE_woK0_n_Sm_side_cut[1]->Draw("colsame");
  
  //
  ///TCanvas *cIMnpim_IMnpip_dE_woK0_n_px = new TCanvas("cIMnpim_IMnpip_dE_woK0_n_px","IMnpim_IMnpip_dE_woK0_n_px");
  //cIMnpim_IMnpip_dE_woK0_n_px->cd();
  TH1D *IMnpim_IMnpip_dE_woK0_n_px = IMnpim_IMnpip_dE_woK0_n->ProjectionX();
  IMnpim_IMnpip_dE_woK0_n_px->GetYaxis()->SetTitle("Counts/ 10 MeV/c^{2}");
  IMnpim_IMnpip_dE_woK0_n_px->GetYaxis()->CenterTitle();
  // IMnpim_IMnpip_dE_woK0_n_px->Draw("HE");
  TH1D* IMnpim_IMnpip_dE_woK0_kin_px_clone[2];
  IMnpim_IMnpip_dE_woK0_kin_px_clone[0] = (TH1D*) IMnpim_IMnpip_dE_woK0_kin_px[0]->Clone();
  IMnpim_IMnpip_dE_woK0_kin_px_clone[0]->SetLineColor(2);
  //IMnpim_IMnpip_dE_woK0_kin_px_clone[0]->Draw("HEsame");
  std::cout << "IMpip " << IMnpim_IMnpip_dE_woK0_kin_px_clone[0]->GetMean() << std::endl;
  std::cout << "IMpip " << IMnpim_IMnpip_dE_woK0_kin_px_clone[0]->GetRMS() << std::endl;
  IMnpim_IMnpip_dE_woK0_kin_px_clone[1] = (TH1D*) IMnpim_IMnpip_dE_woK0_kin_px[1]->Clone();
  IMnpim_IMnpip_dE_woK0_kin_px_clone[1]->SetLineColor(3);
  //IMnpim_IMnpip_dE_woK0_kin_px_clone[1]->Draw("HEsame");
  
  //
  TCanvas *cIMnpim_IMnpip_dE_woK0_n_px2 = new TCanvas("cIMnpim_IMnpip_dE_woK0_n_px2","IMnpim_IMnpip_dE_woK0_n_px2");
  cIMnpim_IMnpip_dE_woK0_n_px2->cd();
  IMnpim_IMnpip_dE_woK0_n_px->GetXaxis()->SetRangeUser(1.0,1.3);
  IMnpim_IMnpip_dE_woK0_n_px->Draw("EH");
  
  TH1D *IMnpim_IMnpip_dE_woK0_n_px_1 = (TH1D*)IMnpim_IMnpip_dE_woK0_n_px->Clone();
  IMnpim_IMnpip_dE_woK0_n_px_1->GetXaxis()->SetRangeUser(anacuts::Sigmap_MIN,anacuts::Sigmap_MAX);
  IMnpim_IMnpip_dE_woK0_n_px_1->SetFillColor(2);
  IMnpim_IMnpip_dE_woK0_n_px_1->Draw("HEsame");
  TH1D *IMnpim_IMnpip_dE_woK0_n_px_2 = (TH1D*)IMnpim_IMnpip_dE_woK0_n_px->Clone();
  IMnpim_IMnpip_dE_woK0_n_px_2->GetXaxis()->SetRangeUser(anacuts::Sigmap_sidelow_MIN,anacuts::Sigmap_sidelow_MAX);
  IMnpim_IMnpip_dE_woK0_n_px_2->SetFillColor(3);
  IMnpim_IMnpip_dE_woK0_n_px_2->Draw("HEsame");
  TH1D *IMnpim_IMnpip_dE_woK0_n_px_3 = (TH1D*)IMnpim_IMnpip_dE_woK0_n_px->Clone();
  IMnpim_IMnpip_dE_woK0_n_px_3->GetXaxis()->SetRangeUser(anacuts::Sigmap_sidehigh_MIN,anacuts::Sigmap_sidehigh_MAX);
  IMnpim_IMnpip_dE_woK0_n_px_3->SetFillColor(3);
  IMnpim_IMnpip_dE_woK0_n_px_3->Draw("HEsame");
   
  TCanvas *cIMnpim_IMnpip_dE_woK0_n_Sp_py = new TCanvas("cIMnpim_IMnpip_dE_woK0_n_Sp_py","IMnpim_IMnpip_dE_woK0_n_Sp_py");
  cIMnpim_IMnpip_dE_woK0_n_Sp_py->cd();
  TH1D* IMnpim_IMnpip_dE_woK0_n_Sp_py = IMnpim_IMnpip_dE_woK0_n_Sp->ProjectionY();
  IMnpim_IMnpip_dE_woK0_n_Sp_py->Draw("EH");
  TH1D* IMnpim_IMnpip_dE_woK0_n_Sp_side_low_py = IMnpim_IMnpip_dE_woK0_n_Sp_side[0]->ProjectionY();
  IMnpim_IMnpip_dE_woK0_n_Sp_side_low_py->SetLineColor(2);
  std::cout << "Sp side low_low " << 
  IMnpim_IMnpip_dE_woK0_n_Sp_side_low_py->
  Integral(IMnpim_IMnpip_dE_woK0_n_Sp_side_low_py->
  FindBin(anacuts::Sigmam_MIN_wide),
  IMnpim_IMnpip_dE_woK0_n_Sp_side_low_py->FindBin(anacuts::Sigmam_MIN)-1) 
  <<std::endl;
  std::cout << "Sp side low_high " << 
  IMnpim_IMnpip_dE_woK0_n_Sp_side_low_py->
  Integral(IMnpim_IMnpip_dE_woK0_n_Sp_side_low_py->
  FindBin(anacuts::Sigmam_MAX),
  IMnpim_IMnpip_dE_woK0_n_Sp_side_low_py->FindBin(anacuts::Sigmam_MAX_wide)-1) 
  <<std::endl;
  //std::cout << "scaling factor Sp low_low " << (anacuts::Sigmam_MIN-anacuts::Sigmam_MIN_wide)/anacuts::Sigmam_MIN << std::endl;
  //std::cout << "scaling factor Sp low_high " << (anacuts::Sigmam_MAX_wide-anacuts::Sigmam_MAX)  << std::endl;
  
  IMnpim_IMnpip_dE_woK0_n_Sp_side_low_py->Draw("HEsame");
  TH1D* IMnpim_IMnpip_dE_woK0_n_Sp_side_high_py = IMnpim_IMnpip_dE_woK0_n_Sp_side[1]->ProjectionY();
  IMnpim_IMnpip_dE_woK0_n_Sp_side_high_py->SetLineColor(3);
  std::cout << "Sp side high_low " << 
  IMnpim_IMnpip_dE_woK0_n_Sp_side_high_py->
  Integral(IMnpim_IMnpip_dE_woK0_n_Sp_side_high_py->
  FindBin(anacuts::Sigmam_MIN_wide),
  IMnpim_IMnpip_dE_woK0_n_Sp_side_high_py->FindBin(anacuts::Sigmam_MIN)-1) 
  <<std::endl;
  std::cout << "Sp side high_high " << 
  IMnpim_IMnpip_dE_woK0_n_Sp_side_high_py->
  Integral(IMnpim_IMnpip_dE_woK0_n_Sp_side_high_py->
  FindBin(anacuts::Sigmam_MAX),
  IMnpim_IMnpip_dE_woK0_n_Sp_side_high_py->FindBin(anacuts::Sigmam_MAX_wide)-1) 
  <<std::endl;
  IMnpim_IMnpip_dE_woK0_n_Sp_side_high_py->Draw("HEsame");
  TH1D* IMnpim_IMnpip_dE_woK0_n_Sp_side_cut_low_py = IMnpim_IMnpip_dE_woK0_n_Sp_side_cut[0]->ProjectionY();
  IMnpim_IMnpip_dE_woK0_n_Sp_side_cut_low_py->SetLineColor(4);
  IMnpim_IMnpip_dE_woK0_n_Sp_side_cut_low_py->Draw("HEsame");
  TH1D* IMnpim_IMnpip_dE_woK0_n_Sp_side_cut_high_py = IMnpim_IMnpip_dE_woK0_n_Sp_side_cut[1]->ProjectionY();
  IMnpim_IMnpip_dE_woK0_n_Sp_side_cut_high_py->SetLineColor(5);
  IMnpim_IMnpip_dE_woK0_n_Sp_side_cut_high_py->Draw("HEsame");
  
  TCanvas *cIMnpim_IMnpip_dE_woK0_n_Sm_px = new TCanvas("cIMnpim_IMnpip_dE_woK0_n_Sm_px","IMnpim_IMnpip_dE_woK0_n_Sm_px");
  cIMnpim_IMnpip_dE_woK0_n_Sm_px->cd();
  TH1D* IMnpim_IMnpip_dE_woK0_n_Sm_px = IMnpim_IMnpip_dE_woK0_n_Sm->ProjectionX();
  IMnpim_IMnpip_dE_woK0_n_Sm_px->Draw("EH");
  TH1D* IMnpim_IMnpip_dE_woK0_n_Sm_side_low_px = IMnpim_IMnpip_dE_woK0_n_Sm_side[0]->ProjectionX();
  IMnpim_IMnpip_dE_woK0_n_Sm_side_low_px->SetLineColor(2);
  
  std::cout << "Sm side low_low " << 
  IMnpim_IMnpip_dE_woK0_n_Sm_side_low_px->
  Integral(IMnpim_IMnpip_dE_woK0_n_Sm_side_low_px->
  FindBin(anacuts::Sigmap_MIN_wide),
  IMnpim_IMnpip_dE_woK0_n_Sm_side_low_px->FindBin(anacuts::Sigmap_MIN)-1) 
  <<std::endl;
  std::cout << "Sm side low_high " << 
  IMnpim_IMnpip_dE_woK0_n_Sm_side_low_px->
  Integral(IMnpim_IMnpip_dE_woK0_n_Sm_side_low_px->
  FindBin(anacuts::Sigmap_MAX),
  IMnpim_IMnpip_dE_woK0_n_Sm_side_low_px->FindBin(anacuts::Sigmap_MAX_wide)-1) 
  <<std::endl;
  IMnpim_IMnpip_dE_woK0_n_Sm_side_low_px->Draw("HEsame");
  
  TH1D* IMnpim_IMnpip_dE_woK0_n_Sm_side_high_px = IMnpim_IMnpip_dE_woK0_n_Sm_side[1]->ProjectionX();
  IMnpim_IMnpip_dE_woK0_n_Sm_side_high_px->SetLineColor(3);
  
  std::cout << "Sm side high_low " << 
  IMnpim_IMnpip_dE_woK0_n_Sm_side_high_px->
  Integral(IMnpim_IMnpip_dE_woK0_n_Sm_side_high_px->
  FindBin(anacuts::Sigmap_MIN_wide),
  IMnpim_IMnpip_dE_woK0_n_Sm_side_high_px->FindBin(anacuts::Sigmap_MIN)-1) 
  <<std::endl;
  std::cout << "Sm side high_high " << 
  IMnpim_IMnpip_dE_woK0_n_Sm_side_high_px->
  Integral(IMnpim_IMnpip_dE_woK0_n_Sm_side_high_px->
  FindBin(anacuts::Sigmap_MAX),
  IMnpim_IMnpip_dE_woK0_n_Sm_side_high_px->FindBin(anacuts::Sigmap_MAX_wide)-1) 
  <<std::endl;
  
  IMnpim_IMnpip_dE_woK0_n_Sm_side_high_px->Draw("HEsame");
  TH1D* IMnpim_IMnpip_dE_woK0_n_Sm_side_cut_low_px = IMnpim_IMnpip_dE_woK0_n_Sm_side_cut[0]->ProjectionX();
  IMnpim_IMnpip_dE_woK0_n_Sm_side_cut_low_px->SetLineColor(4);
  IMnpim_IMnpip_dE_woK0_n_Sm_side_cut_low_px->Draw("HEsame");
  TH1D* IMnpim_IMnpip_dE_woK0_n_Sm_side_cut_high_px = IMnpim_IMnpip_dE_woK0_n_Sm_side_cut[1]->ProjectionX();
  IMnpim_IMnpip_dE_woK0_n_Sm_side_cut_high_px->SetLineColor(5);
  IMnpim_IMnpip_dE_woK0_n_Sm_side_cut_high_px->Draw("HEsame");


  //TCanvas *cIMnpim_IMnpip_dE_woK0_n_py = new TCanvas("cIMnpim_IMnpip_dE_woK0_n_py","IMnpim_IMnpip_dE_woK0_n_py");
  //cIMnpim_IMnpip_dE_woK0_n_py->cd();
  TH1D *IMnpim_IMnpip_dE_woK0_n_py = IMnpim_IMnpip_dE_woK0_n->ProjectionY();
  IMnpim_IMnpip_dE_woK0_n_py->GetYaxis()->SetTitle("Counts/ 10 MeV/c^{2}");
  IMnpim_IMnpip_dE_woK0_n_py->GetYaxis()->CenterTitle();
  IMnpim_IMnpip_dE_woK0_n_py->Draw("HE");
  TH1D *IMnpim_IMnpip_dE_woK0_kin_py_clone[2];
  IMnpim_IMnpip_dE_woK0_kin_py_clone[0] = (TH1D*)IMnpim_IMnpip_dE_woK0_kin_py[0]->Clone();
  IMnpim_IMnpip_dE_woK0_kin_py_clone[0]->SetLineColor(2);
  //IMnpim_IMnpip_dE_woK0_kin_py_clone[0]->Draw("HEsame");
  IMnpim_IMnpip_dE_woK0_kin_py_clone[1] = (TH1D*)IMnpim_IMnpip_dE_woK0_kin_py[1]->Clone();
  IMnpim_IMnpip_dE_woK0_kin_py_clone[1]->SetLineColor(3);
  //IMnpim_IMnpip_dE_woK0_kin_py_clone[1]->Draw("HEsame");
  std::cout << "IMpim " << IMnpim_IMnpip_dE_woK0_kin_py_clone[1]->GetMean() << std::endl;
  std::cout << "IMpim " << IMnpim_IMnpip_dE_woK0_kin_py_clone[1]->GetRMS() << std::endl;
   
  TCanvas *cIMnpim_IMnpip_dE_woK0_n_py2 = new TCanvas("cIMnpim_IMnpip_dE_woK0_n_py2","IMnpim_IMnpip_dE_woK0_n_py2");
  cIMnpim_IMnpip_dE_woK0_n_py2->cd();
  IMnpim_IMnpip_dE_woK0_n_py->GetXaxis()->SetRangeUser(1.0,1.3);
  IMnpim_IMnpip_dE_woK0_n_py->Draw("EH");
  
  
  TH1D *IMnpim_IMnpip_dE_woK0_n_py_1 = (TH1D*)IMnpim_IMnpip_dE_woK0_n_py->Clone();
  IMnpim_IMnpip_dE_woK0_n_py_1->GetXaxis()->SetRangeUser(anacuts::Sigmam_MIN,anacuts::Sigmam_MAX);
  IMnpim_IMnpip_dE_woK0_n_py_1->SetFillColor(2);
  IMnpim_IMnpip_dE_woK0_n_py_1->Draw("HEsame");
  
  TH1D *IMnpim_IMnpip_dE_woK0_n_py_2 = (TH1D*)IMnpim_IMnpip_dE_woK0_n_py->Clone();
  IMnpim_IMnpip_dE_woK0_n_py_2->GetXaxis()->SetRangeUser(anacuts::Sigmam_sidelow_MIN,anacuts::Sigmam_sidelow_MAX);
  IMnpim_IMnpip_dE_woK0_n_py_2->SetFillColor(3);
  IMnpim_IMnpip_dE_woK0_n_py_2->Draw("HEsame");
  TH1D *IMnpim_IMnpip_dE_woK0_n_py_3 = (TH1D*)IMnpim_IMnpip_dE_woK0_n_py->Clone();
  IMnpim_IMnpip_dE_woK0_n_py_3->GetXaxis()->SetRangeUser(anacuts::Sigmam_sidehigh_MIN,anacuts::Sigmam_sidehigh_MAX);
  IMnpim_IMnpip_dE_woK0_n_py_3->SetFillColor(3);
  IMnpim_IMnpip_dE_woK0_n_py_3->Draw("HEsame");
  
  /*
  TCanvas *cIMnpim_IMnpip_dE_woK0_n_Sp_side_sum = new TCanvas("cIMnpim_IMnpip_dE_woK0_n_Sp_side_sum","IMnpim_IMnpip_dE_woK0_n_Sp_side_sum");
  IMnpim_IMnpip_dE_woK0_n_Sp_side_sum->GetXaxis()->SetRangeUser(0,1.7);
  IMnpim_IMnpip_dE_woK0_n_Sp_side_sum->GetYaxis()->SetRangeUser(0,1.7);
  IMnpim_IMnpip_dE_woK0_n_Sp_side_sum->Draw("colz");

  TCanvas *cIMnpim_IMnpip_dE_woK0_n_Sm_side_sum = new TCanvas("cIMnpim_IMnpip_dE_woK0_n_Sm_side_sum","IMnpim_IMnpip_dE_woK0_n_Sm_side_sum");
  IMnpim_IMnpip_dE_woK0_n_Sm_side_sum->GetXaxis()->SetRangeUser(0,1.7);
  IMnpim_IMnpip_dE_woK0_n_Sm_side_sum->GetYaxis()->SetRangeUser(0,1.7);
  IMnpim_IMnpip_dE_woK0_n_Sm_side_sum->Draw("colz");
  */

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
  DCA_pip_beam_kin[0]->Draw("");
  DCA_pim_beam_kin[0]->SetLineColor(2);
  DCA_pim_beam_kin[0]->Draw("same");
  DCA_pip_pim_kin[0]->SetLineColor(3);
  DCA_pip_pim_kin[0]->Draw("same");
  std::cout << DCA_pip_beam_kin[0]->Integral() << std::endl;
  std::cout << DCA_pim_beam_kin[0]->Integral() << std::endl;
  std::cout << DCA_pip_pim_kin[0]->Integral() << std::endl;
  
  TCanvas *c_DCA_Sm = new TCanvas("c_DCA_Sm","c_DCA_Sm");
  c_DCA_Sm->cd();
  DCA_pip_beam_kin[1]->Draw("");
  DCA_pim_beam_kin[1]->SetLineColor(2);
  DCA_pim_beam_kin[1]->Draw("same");
  DCA_pip_pim_kin[1]->SetLineColor(3);
  DCA_pip_pim_kin[1]->Draw("same");
  std::cout << DCA_pip_beam_kin[1]->Integral() << std::endl;
  std::cout << DCA_pim_beam_kin[1]->Integral() << std::endl;
  std::cout << DCA_pip_pim_kin[1]->Integral() << std::endl;
   */
  
  //acceptance calculation
  
  TFile *facc = new TFile("acc.root","READ");
  
  if(Spmode || Smmode){
    if(Spmode){
      std::cout << "This is Sigma+ mode sim." << std::endl;
    }else{
      std::cout << "This is Sigma- mode sim." << std::endl;
    }
    //sphis = new TFile("../simpost/simIMpisigma_nSppim_DoraAir_v28.root","READ");
    //sphis = new TFile("../simpost/simIMpisigma_nSppim_DoraAir_v32.root","READ");
    TFile *genhis = new TFile("../simpost/simIMpisigma_nSppim_DoraAir_v28_v32.root","READ");
    std::cout << "file for generated info " ;
    std::cout << genhis->GetName() << std::endl;
    TString sacc = genhis->GetName();
    sacc.Replace(sacc.Length()-5,10,"_acc.root");
    std::cout << "acc file Sp mode: " << sacc.Data() << std::endl;
   
    //facc->SetName(sacc.Data());
    genhis->cd();
    TCanvas *cphase = new TCanvas("cphase","cphase");
    cphase->cd();
    TH2F *React_q_IMPiSigma = (TH2F*)genhis->Get("React_q_IMPiSigma");
    React_q_IMPiSigma->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
    React_q_IMPiSigma->SetYTitle("Mom. tranfer [GeV/c]");
    if(Spmode)React_q_IMPiSigma->SetTitle("Generated Events (#pi^{-}#Sigma^{+} mode) ");
    if(Smmode)React_q_IMPiSigma->SetTitle("Generated Events (#pi^{+}#Sigma^{-} mode) ");
    //React_q_IMPiSigma->RebinX(1);
    React_q_IMPiSigma->RebinY(12);
    React_q_IMPiSigma->Draw("colz");
   
    TCanvas *ceff = new TCanvas("ceff","ceff");
    //q_IMnpipi_woK0_wSid_n_Sp->RebinY(12);
    TH2F *h2acc=NULL;
    if(Spmode) h2acc =  (TH2F*)q_IMnpipi_woK0_wSid_n_Sp_acc->Clone();
    else       h2acc =  (TH2F*)q_IMnpipi_woK0_wSid_n_Sm_acc->Clone();
    h2acc->Sumw2();
    h2acc->RebinY(12);
    q_IMnpipi_woK0_wSid_n_Sp->Sumw2();
    q_IMnpipi_woK0_wSid_n_Sm->Sumw2();
    React_q_IMPiSigma->Sumw2();
    React_q_IMPiSigma->Print("base");
    q_IMnpipi_woK0_wSid_n_Sp_acc->RebinY(12);
    q_IMnpipi_woK0_wSid_n_Sm_acc->RebinY(12);
    q_IMnpipi_woK0_wSid_n_Sp_acc->Print("base");
    std::cout << "calc. acc." << std::endl;
    if(Spmode){
      h2acc->Divide(q_IMnpipi_woK0_wSid_n_Sp_acc,React_q_IMPiSigma,1.0,1.0,"b");
    }else{
      h2acc->Divide(q_IMnpipi_woK0_wSid_n_Sm_acc,React_q_IMPiSigma,1.0,1.0,"b");
    }
    h2acc->SetMaximum(0.02);
    h2acc->Draw("colz");
    TFile *fsacc = new TFile(sacc.Data(),"RECREATE");
    fsacc->cd();
    h2acc->SetName("acc_Sp");
    h2acc->SetTitle("acc_Sp");
    h2acc->Write();
    TH2F *acc_err = new TH2F("acc_err","acc_err",500,1,2,300,0,1.5);
    acc_err->RebinY(12);
    for(int ix=0;ix<h2acc->GetNbinsX();ix++){
      for(int iy=0;iy<h2acc->GetNbinsY();iy++){
        double err = h2acc->GetBinErrorUp(ix,iy);
        double cont = h2acc->GetBinContent(ix,iy);
        if(cont) acc_err->SetBinContent(ix,iy,err/cont);
      }
    }
    TCanvas *cacc_err = new TCanvas("cacc_err","acc_err");
    acc_err->Draw("colz");
    acc_err->Write();
    if(Spmode)React_q_IMPiSigma->SetName("React_q_IMPiSigma_Sp");
    else      React_q_IMPiSigma->SetName("React_q_IMPiSigma_Sm");
    React_q_IMPiSigma->Write();
    if(Spmode){
      q_IMnpipi_woK0_wSid_n_Sp_acc->Write();
    }else{
      q_IMnpipi_woK0_wSid_n_Sm_acc->Write();
    }
    fsacc->Close();
  }//Spmode or Smmode
  
  
  facc->cd();
  TH2F* acc_Sp_cal = (TH2F*)facc->Get("acc_Sp");
  if(acc_Sp_cal == NULL){
    std::cout << " acc_Sp is NULL " << std::endl;
  }
  /*
  for(int i=0;i<acc_Sp_cal->GetNbinsX();i++){
    for(int j=0;j<acc_Sp_cal->GetNbinsY();j++){
      double val = acc_Sp_cal->GetBinContent(i,j);
      if(val>0.02) {
        std::cout << "(x,y)" << i << " " << j  << std::endl;
      }
    }
  }*/
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
  TH2F *q_IMnpipi_woK0_wSid_n_Sp_side_cs = (TH2F*)q_IMnpipi_woK0_wSid_n_Sp_side->Clone("Sp_cs_side");
  q_IMnpipi_woK0_wSid_n_Sp_side_cs->Sumw2();
  q_IMnpipi_woK0_wSid_n_Sp_side_cs->Divide(acc_Sp_cal);
  q_IMnpipi_woK0_wSid_n_Sp_side_cs->SetMaximum(30000);
  q_IMnpipi_woK0_wSid_n_Sp_side_cs->Draw("colz");

  TCanvas *cq_IMnpipi_woK0_wSid_n_Sm_side_cs = new TCanvas("cq_IMnpipi_woK0_wSid_n_Sm_side_cs","q_IMnpipi_woK0_wSid_n_Sm_side_cs");
  cq_IMnpipi_woK0_wSid_n_Sm_side_cs->cd();
  TH2F *q_IMnpipi_woK0_wSid_n_Sm_side_cs = (TH2F*)q_IMnpipi_woK0_wSid_n_Sm_side->Clone("Sm_cs_side");
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

  //TCanvas *cq_IMnpipi_woK0_wSid_n_SpSm_side_cs_px = new TCanvas("cq_IMnpipi_woK0_wSid_n_SpSm_side_cs_px","q_IMnpipi_woK0_wSid_n_SpSm_side_cs_px");
  //cq_IMnpipi_woK0_wSid_n_SpSm_side_cs_px->cd();
  //TH1D* q_IMnpipi_woK0_wSid_n_SpSm_side_cs_px = q_IMnpipi_woK0_wSid_n_SpSm_side_cs->ProjectionX();
  //q_IMnpipi_woK0_wSid_n_SpSm_side_cs_px->Draw("E");


  
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
    }
    if(obj->InheritsFrom("TH1D")){
      h1d = (TH1D*) obj;
      h1d->GetXaxis()->CenterTitle();
      //h1d->GetXaxis()->SetTitleSize(0.05);
      //h1d->GetXaxis()->SetTitleOffset(0.80);
    }
    if(obj->InheritsFrom("TH2")){
      h2 = (TH2F*) obj;
      h2->GetXaxis()->CenterTitle();
      h2->GetYaxis()->CenterTitle();
      //h2->GetXaxis()->SetTitleSize(0.05);
      //h2->GetXaxis()->SetTitleOffset(0.80);
      //h2->GetYaxis()->SetTitleSize(0.05);
      //h2->GetYaxis()->SetTitleOffset(0.85);
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
      pt->AddText("MC #Sigma+ mode");
    }
    else if(Smmode){
      pt->SetFillColor(kAzure-4);
      pt->AddText("MC #Sigma- mode"); 
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
 // pdf->Close();
  std::cout << "closing pdf " << std::endl;
  //TString outname = std::string(filename);
  //outname.Replace(std::string(filename).size()-4,5,"out.root");
  //TFile *fout = new TFile(outname.Data(),"RECREATE");
  //fout->cd();
  //q_IMnpipi_woK0_wSid_n->Write();
  //fout->Close();
  
}
