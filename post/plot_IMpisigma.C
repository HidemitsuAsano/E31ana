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
#include <TPaveText.h>
#include <TLatex.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>

#include "../src/GlobalVariables.h"
#include "anacuts.h"
#include "globalana.h"

const double pvalcut = 0.005;
const bool gridon=true;
const bool staton=true;
const bool UseKinFit = false;
const bool UseKinFitVal = true;
const double lumi = 56.0e9; // rough value
const double Beamsurvival =  0.4; //rough value from Kawasaki ana
const double D2density = 4.83e23; //kawasaki ana
const double DAQeff = 0.77;//kawasaki ana
const double trigeff = 2.0; //pre-scale factor of CDH3 trigger
const int ngap=1;
const int nzone=2;//for wide range side band study
const unsigned int LOWside=0;
const unsigned int HIGHside=1;

//0: diagonal cut
//1: 3 sigma cut
//2: 5 simga cut 
const unsigned int sigmacuttype=0;

//0:diagonal cut
//1:3 sigma cut
//2:5 sigma cut
const unsigned int sidebandtype=0;

//color def. 
//Sp mode Signal :2 (red)
//Sm mode Signal :3 (green)
//Sp mode sideband Cyan+4
//Sp mode sideband low mass side :kCyan
//Sp mode sideband high mass side :kCyan+2
//Sm mode sideband : kPink+10
//Sm mode sideband low mass side kPink+1
//Sm mode sideband high mass side kPink+3
//neutron mass : 4 (blue)

void plot_IMpisigma(const char* filename="",const int qvalcutflag=0)
{
  gROOT->SetBatch(true);
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
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadLeftMargin(0.12);
  //gStyle->SetTitleFontSize(0.1);
  
  TH1::SetDefaultSumw2();

  std::cout << "infile " << filename <<std::endl;
  TString pdfname = std::string(filename);
  if(qvalcutflag==0) pdfname.Replace(std::string(filename).size()-4,5,"pdf");
  if(qvalcutflag==1) pdfname.Replace(std::string(filename).size()-5,8,"_1.pdf");
  if(qvalcutflag==2) pdfname.Replace(std::string(filename).size()-5,8,"_2.pdf");
  std::cout << "pdfname: " << pdfname << std::endl;
  std::cout << std::endl;

  std::cout << "Use Kin Fit ? " << std::endl; 
  if(UseKinFit) std::cout << "Yes" << std::endl;
  else             std::cout << "No"  << std::endl;
  
  std::cout << "Use Kin Fit Val ? " << std::endl; 
  if(UseKinFitVal) std::cout << "Yes" << std::endl;
  else             std::cout << "No"  << std::endl;
  
  std::cout << std::endl;
  std::cout << "Sigma selection type     " << sigmacuttype << std::endl;
  std::cout << "Side band selection type " << sidebandtype << std::endl;

  bool SimSpmode = (std::string(filename).find("Sp")!= std::string::npos);
  bool SimSmmode = (std::string(filename).find("Sm")!= std::string::npos);
  bool SimK0nn = (std::string(filename).find("K0nn")!= std::string::npos);
  bool SimK0n_ns = (std::string(filename).find("K0n_ns")!= std::string::npos);
  bool SimnpipiL = (std::string(filename).find("npipiL")!= std::string::npos);

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
  tree->SetBranchAddress( "CDH_Pos",&CDH_Pos);
  //tree->SetBranchAddress( "run_num", &run_num );
  //tree->SetBranchAddress( "event_num", &event_num );
  //tree->SetBranchAddress( "block_num", &block_num );
  if(SimSpmode || SimSmmode){
    tree->SetBranchAddress( "mcmom_beam",  &mcmom_beam );
    tree->SetBranchAddress( "mcmom_pip", &mcmom_pip);
    tree->SetBranchAddress( "mcmom_pim", &mcmom_pim);
    tree->SetBranchAddress( "mcmom_ncds", &mcmom_ncds);
    tree->SetBranchAddress( "mcmom_nmiss", &mcmom_nmiss);
    tree->SetBranchAddress( "react_nmiss", &react_nmiss);
    tree->SetBranchAddress( "react_Sigma", &react_Sigma);
    tree->SetBranchAddress( "react_pi", &react_pi);
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
  TH2F* dE_nmom_fid_beta;
  TH2F* dE_MMom_fid_beta;
  TH2F* dE_MMass_fid_beta;
  TH2F* dE_MMass_fid_beta_wSid;
  TH2F* dE_MMom_fid_beta_woK0;
  TH2F* dE_MMass_fid_beta_woK0;
  TH2F* dE_MMass_fid_beta_woK0_wSid;
  TH2F* MMom_MMass;
  TH2F* MMom_MMass_wSid;
  TH2F* MMom_MMass_woK0;
  TH2F* MMom_MMass_woK0_wSid;
  TH2F* IMnpim_IMnpip_dE;
  TH2F* IMnpim_IMnpip_dE_woK0;
  TH2F* IMnpim_IMnpip_dE_n;//
  TH2F* IMnpim_IMnpip_dE_woK0_n;//
  TH2F* IMnpip_CDHphi_dE_woK0_n;//
  TH2F* IMnpip_CDHz_dE_woK0_n;//
  TH2F* IMnpim_CDHphi_dE_woK0_n;//
  TH2F* IMnpim_CDHz_dE_woK0_n;//
  TH2F* IMnpip_DCApip_dE_woK0_n;//
  TH2F* IMnpip_DCApim_dE_woK0_n;//
  TH2F* IMnpim_DCApip_dE_woK0_n;//
  TH2F* IMnpim_DCApim_dE_woK0_n;//
  TH2F* IMnpim_IMnpip_dE_n_Sp;//Spmode, avoiding crossing-point
  TH2F* IMnpim_IMnpip_dE_woK0_n_Sp;//Spmode, avoiding crossing-point
  TH2F* IMnpim_IMnpip_dE_woK0_n_Sp_bg;//Spmode+background region, avoiding crossing-point
  TH2F* IMnpim_IMnpip_dE_n_Sm;//Smmode, avoiding crossing-point
  TH2F* IMnpim_IMnpip_dE_woK0_n_Sm;//Smmode, avoiding crossing-point
  TH2F* IMnpim_IMnpip_dE_woK0_n_Sm_bg;//Smmode, avoiding crossing-point
  TH2F* IMnpim_IMnpip_dE_n_side[ngap];//Side band for Sp + Sm mode
  TH2F* IMnpim_IMnpip_dE_woK0_n_side[ngap];//Side band for Sp + Sm mode
  TH2F* IMnpim_IMnpip_dE_n_Sp_side[2][ngap];//low mass,high mass 
  TH2F* IMnpim_IMnpip_dE_n_Sm_side[2][ngap];//low mass,high mass
  TH2F* IMnpim_IMnpip_dE_woK0_n_Sp_side[2][ngap];//low mass,high mass 
  TH2F* IMnpim_IMnpip_dE_woK0_n_Sm_side[2][ngap];//low mass,high mass
  TH2F* IMnpim_IMnpip_dE_n_Sp_sidewide[nzone];//scan from low mass side 
  TH2F* IMnpim_IMnpip_dE_n_Sm_sidewide[nzone];//scan from low mass side
  TH2F* IMnpim_IMnpip_dE_woK0_n_Sp_sidewide[nzone];//scan from low mass side 
  TH2F* IMnpim_IMnpip_dE_woK0_n_Sm_sidewide[nzone];//scan from low mass side
  TH2F* nmom_IMnpim_dE_n;
  TH2F* nmom_IMnpim_dE_woK0_n;
  TH2F* nmom_IMnpip_dE_n;
  TH2F* nmom_IMnpip_dE_woK0_n;
  TH2F* MMnmiss_IMpippim_dE;
  TH2F* MMnmiss_IMnpip_dE;
  TH2F* MMnmiss_IMnpim_dE;
  TH2F* MMnmiss_IMnpip_dE_woK0;
  TH2F* MMnmiss_IMnpim_dE_woK0;
  TH2F* MMnpip_MMnpim_n;
  TH2F* MMnpip_MMnpim_wSid_n;
  TH2F* MMnpip_MMnpim_woK0_n;
  TH2F* MMnpip_MMnpim_woK0_wSid_n;
  TH2F* dE_IMnpim;
  TH2F* dE_IMnpim_n;
  TH2F* dE_IMnpip;
  TH2F* dE_IMnpip_n;
  TH2F* dE_IMnpim_woK0;
  TH2F* dE_IMnpim_woK0_n;
  TH2F* dE_IMnpip_woK0;
  TH2F* dE_IMnpip_woK0_n;
  TH2F* dE_IMnpipi_wSid_n;
  TH2F* dE_IMnpipi_woK0_wSid_n;
  TH2F* Cosn_IMnpipi_wSid_n;
  TH2F* Cosn_IMnpipi_woK0_wSid_n;
  TH2F* MMnmiss_IMnpipi_wSid_n;
  TH2F* MMnmiss_IMnpipi_woK0_wSid;
  TH2F* MMnmiss_IMnpipi_woK0_wSid_Sp;
  TH2F* MMnmiss_IMnpipi_woK0_wSid_Sm;
  TH2F* q_IMnpipi_wSid_n;
  TH2F* q_IMnpipi_woK0_wSid_n;
  TH2F* q_IMpiSigma_gen;//fine bins
  TH2F* q_IMpiSigma_wSid_n_genacc;//fine bins,  reaction data
  TH2F* q_IMnpipi_wSid_n_acc;//fine bins, McData
  TH2F* q_IMnpipi_wSid_n_acc_reco;//fine bins, reconstructed value
  TH2F* q_IMpiSigma_woK0_wSid_n_genacc;//fine bins, reaction data 
  TH2F* q_IMnpipi_woK0_wSid_n_acc;//fine bins, no cuts for separationg S+/S-
  TH2F* q_IMnpipi_woK0_wSid_n_acc_reco;//fine bins, no cuts for separationg S+/S-, reconstructed value
  TH2F* q_IMnpipi_wSid_n_Sp;
  TH2F* q_IMnpipi_wK0_n;
  TH2F* q_IMnpipi_wK0_wSid_n_Sp;
  TH2F* q_IMnpipi_woK0_wSid_n_Sp;
  TH2F* q_IMpiSigma_wSid_n_Sp_genacc;//fine bins
  TH2F* q_IMpiSigma_woK0_wSid_n_Sp_genacc;//fine bins
  TH2F* q_IMnpipi_wSid_n_Sp_acc;//fine bins
  TH2F* q_IMnpipi_woK0_wSid_n_Sp_acc;//fine bins
  TH2F* q_IMnpipi_wSid_n_Sp_acc_reco;//fine bins
  TH2F* q_IMnpipi_woK0_wSid_n_Sp_acc_reco;//fine bins
  TH2F* q_IMnpipi_wSid_n_Sp_side[3][2][ngap];//sideband type, 
  TH2F* q_IMnpipi_woK0_wSid_n_Sp_side[3][2][ngap];//sideband type, 
                                                //low high side, 
                                                //distance to signal 
  TH2F* q_IMnpipi_wSid_n_Sp_sidewide[nzone];//
  TH2F* q_IMnpipi_woK0_wSid_n_Sp_sidewide[nzone];//
  TH2F* q_IMnpipi_wSid_n_Sm;
  TH2F* q_IMnpipi_wK0_wSid_n_Sm;
  TH2F* q_IMnpipi_woK0_wSid_n_Sm;
  TH2F* q_IMpiSigma_wSid_n_Sm_genacc;//fine bins, reaction data
  TH2F* q_IMpiSigma_woK0_wSid_n_Sm_genacc;//fine bins, reaction data
  TH2F* q_IMnpipi_wSid_n_Sm_acc;//fine bins, MCdata
  TH2F* q_IMnpipi_woK0_wSid_n_Sm_acc;//fine bins, MCdata
  TH2F* q_IMnpipi_wSid_n_Sm_acc_reco;//fine bins, reconstructed value
  TH2F* q_IMnpipi_woK0_wSid_n_Sm_acc_reco;//fine bins, reconstructed value
  TH2F* q_IMnpipi_wSid_n_Sm_side[3][2][ngap];//sideband type, 
  TH2F* q_IMnpipi_woK0_wSid_n_Sm_side[3][2][ngap];//sideband type, 
                                                //low high side
                                                //distance to signal
  TH2F* q_IMnpipi_wSid_n_Sm_sidewide[nzone];//
  TH2F* q_IMnpipi_woK0_wSid_n_Sm_sidewide[nzone];//
  TH2F* IMnpip_IMnpipi_n;
  TH2F* IMnpim_IMnpipi_n;
  TH2F* IMnpip_IMnpipi_wK0_n;
  TH2F* IMnpim_IMnpipi_wK0_n;
  TH2F* IMnpip_IMnpipi_woK0_n;
  TH2F* IMnpim_IMnpipi_woK0_n;

  TH2F* q_IMnpipi_wSid_n_side[ngap];//side band method
  TH2F* q_IMnpipi_woK0_wSid_n_side[ngap];//side band method
  TH2F* nmom_IMnpipi_wSid_n;
  TH2F* nmom_IMnpipi_woK0_wSid_n;
  TH2F* nmom_IMnpipi_woK0_wSid_n_Sp;
  TH2F* nmom_IMnpipi_woK0_wSid_n_Sm;
  //K0 study
  TH2F* nmom_IMnpipi_wK0_n;
  TH2F* nmom_cosnlab_K0_n;//ncds mom. vs cos ncdds lab. frame K0
  TH2F* nmom_IMpippim_n;
  TH2F* nmom_IMpippim_wSid_n;
  TH2F* mnmom_IMpippim_n;//missing neutron mom.
  TH2F* mnmom_IMpippim_wSid_n;//missing neutron mom.
  TH2F* q_IMpippim_n;
  TH2F* q_IMpippim_n_wSid;
  TH2F* IMpippim_IMnpipi_n;
  TH2F* IMpippim_IMnpipi_n_wSid;
  TH2F* IMpippim_IMnpip_n;
  TH2F* IMpippim_IMnpim_n;
  TH2F* Mompippim_IMnpipi_dE_wK0_n;
  TH2F* Mompippim_IMnpipi_dE_wK0_n_Sp;
  TH2F* Mompippim_IMnpipi_dE_wK0_n_Sm;

  const int nbinIMnpipi = 100;//1-2 GeV/c^2
  const int nbinq = 25;//0-1.5 GeV/c
  const int nbinIMnpi = 400; //1-2 GeV/c^2
  const int nbinnmiss = 100; //0-1.5 GeV/c
  const int nbindE = 200;
  const int nbinpippim = 500;
  
  dE_betainv_fid = new TH2F(Form("dE_betainv_fid"),Form("dE_betainv_fid"),1000, 0, 50, nbindE, 0, 50);
  dE_betainv_fid->SetXTitle("1/#beta");
  dE_betainv_fid->SetYTitle("dE [MeVee]");
  
  dE_nmom_fid_beta = new TH2F(Form("dE_nmom_fid_beta"),Form("dE_nmom_fid_beta"),100, 0, 1.5, nbindE, 0, 50);
  dE_nmom_fid_beta->SetXTitle("Mom. [GeV/c]");
  dE_nmom_fid_beta->SetYTitle("dE [MeVee]");
  
  dE_MMom_fid_beta = new TH2F(Form("dE_MMom_fid_beta"),Form("dE_MMom_fid_beta"),100, 0, 1.5, nbindE, 0, 50);
  dE_MMom_fid_beta->SetXTitle("Missing Mom. [GeV/c]");
  dE_MMom_fid_beta->SetYTitle("dE [MeVee]");
  
  dE_MMom_fid_beta_woK0 = new TH2F(Form("dE_MMom_fid_beta_woK0"),Form("dE_MMom_fid_beta_woK0"),100, 0, 1.5, nbindE, 0, 50);
  dE_MMom_fid_beta_woK0->SetXTitle("Missing Mom. [GeV/c]");
  dE_MMom_fid_beta_woK0->SetYTitle("dE [MeVee]");

  dE_MMass_fid_beta = new TH2F(Form("dE_MMass_fid_beta"),Form("dE_MMass_fid_beta"), 140, 0.4, 1.8, nbindE, 0, 50);
  dE_MMass_fid_beta->SetXTitle("Missing mass [GeV/c^{2}]");
  dE_MMass_fid_beta->SetYTitle("dE [MeVee]");
  
  dE_MMass_fid_beta_woK0 = new TH2F(Form("dE_MMass_fid_beta_woK0"),Form("dE_MMass_fid_beta_woK0"), 140, 0.4, 1.8, nbindE, 0, 50);
  dE_MMass_fid_beta_woK0->SetXTitle("Missing mass [GeV/c^{2}]");
  dE_MMass_fid_beta_woK0->SetYTitle("dE [MeVee]");
  
  dE_MMass_fid_beta_wSid = new TH2F(Form("dE_MMass_fid_beta_wSid"),Form("dE_MMass_fid_beta_wSid"), 140, 0.4, 1.8, nbindE, 0, 10);
  dE_MMass_fid_beta_wSid->SetXTitle("Missing mass [GeV/c^{2}]");
  dE_MMass_fid_beta_wSid->SetYTitle("dE [MeVee]");
  
  dE_MMass_fid_beta_woK0_wSid = new TH2F(Form("dE_MMass_fid_beta_woK0_wSid"),Form("dE_MMass_fid_beta_woK0_wSid"), 140, 0.4, 1.8, nbindE, 0, 10);
  dE_MMass_fid_beta_woK0_wSid->SetXTitle("Missing mass [GeV/c^{2}]");
  dE_MMass_fid_beta_woK0_wSid->SetYTitle("dE [MeVee]");

  MMom_MMass = new TH2F(Form("MMom_MMass"),Form("MMom_MMass"), 140, 0.4, 1.8, 100, 0, 1.5);
  MMom_MMass->SetXTitle("Missing Mass [GeV/c^{2}]");
  MMom_MMass->SetYTitle("Missing Mom. [GeV/c]");
  
  MMom_MMass_woK0 = new TH2F(Form("MMom_MMass_woK0"),Form("MMom_MMass_woK0"), 140, 0.4, 1.8, 100, 0, 1.5);
  MMom_MMass_woK0->SetXTitle("Missing Mass [GeV/c^{2}]");
  MMom_MMass_woK0->SetYTitle("Missing Mom. [GeV/c]");

  MMom_MMass_wSid = new TH2F(Form("MMom_MMass_wSid"),Form("MMom_MMass_wSid"), 140, 0.4, 1.8, 100, 0, 1.5);
  MMom_MMass_wSid->SetXTitle("Missing Mass [GeV/c^{2}]");
  MMom_MMass_wSid->SetYTitle("Missing Mom. [GeV/c]");
  
  MMom_MMass_woK0_wSid = new TH2F(Form("MMom_MMass_woK0_wSid"),Form("MMom_MMass_woK0_wSid"), 140, 0.4, 1.8, 100, 0, 1.5);
  MMom_MMass_woK0_wSid->SetXTitle("Missing Mass [GeV/c^{2}]");
  MMom_MMass_woK0_wSid->SetYTitle("Missing Mom. [GeV/c]");
  
  IMnpim_IMnpip_dE = new TH2F(Form("IMnpim_IMnpip_dE"), Form("IMnpim_IMnpip_dE"),nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  
  IMnpim_IMnpip_dE_woK0 = new TH2F(Form("IMnpim_IMnpip_dE_woK0"), Form("IMnpim_IMnpip_dE_woK0"),nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_woK0->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_woK0->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
   
  MMnmiss_IMpippim_dE = new TH2F("MMnmiss_IMpippim_dE", "MMnmiss_IMpippim_dE",nbinpippim,0.,1.0,nbinnmiss,0,1.5);
  MMnmiss_IMpippim_dE->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMpippim_dE->SetYTitle("Miss Mass. [GeV/c^{2}]");
  
  MMnmiss_IMnpip_dE = new TH2F("MMnmiss_IMnpip_dE", "MMnmiss_IMnpip_dE",nbinIMnpi,1.,2.0,nbinnmiss,0,1.5);
  MMnmiss_IMnpip_dE->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  MMnmiss_IMnpip_dE->SetYTitle("Miss Mass. [GeV/c^{2}]");
  
  MMnmiss_IMnpip_dE_woK0 = new TH2F("MMnmiss_IMnpip_dE_woK0", "MMnmiss_IMnpip_dE_woK0",nbinIMnpi,1.,2.0,nbinnmiss,0,1.5);
  MMnmiss_IMnpip_dE_woK0->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  MMnmiss_IMnpip_dE_woK0->SetYTitle("Miss Mass. [GeV/c^{2}]");

  MMnmiss_IMnpim_dE = new TH2F("MMnmiss_IMnpim_dE", "MMnmiss_IMnpim_dE",nbinIMnpi,1.,2.0,nbinnmiss,0,1.5);
  MMnmiss_IMnpim_dE->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMnpim_dE->SetYTitle("Miss Mass. [GeV/c^{2}]");
  
  MMnmiss_IMnpim_dE_woK0 = new TH2F("MMnmiss_IMnpim_dE_woK0", "MMnmiss_IMnpim_dE_woK0",nbinIMnpi,1.,2.0,nbinnmiss,0,1.5);
  MMnmiss_IMnpim_dE_woK0->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMnpim_dE_woK0->SetYTitle("Miss Mass. [GeV/c^{2}]");

  IMnpim_IMnpip_dE_n = new TH2F(Form("IMnpim_IMnpip_dE_n"),Form("IMnpim_IMnpip_dE_n"),nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_n->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_n->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  
  IMnpim_IMnpip_dE_woK0_n = new TH2F(Form("IMnpim_IMnpip_dE_woK0_n"),Form("IMnpim_IMnpip_dE_woK0_n"),nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_woK0_n->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_woK0_n->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  
  IMnpip_CDHphi_dE_woK0_n = new TH2F("IMnpip_CDHphi_dE_woK0_n","IMnpip_CDHphi_dE_woK0_n",100,-3.5,3.5,nbinIMnpi,1.,2.);
  IMnpip_CDHphi_dE_woK0_n->SetXTitle("CDHphi [radian]");
  IMnpip_CDHphi_dE_woK0_n->SetYTitle("IM(n#pi^{+} [GeV/c^{2}]");
  
  IMnpip_CDHz_dE_woK0_n = new TH2F("IMnpip_CDHz_dE_woK0_n","IMnpip_CDHz_dE_woK0_n",100,-50,50,nbinIMnpi,1.,2.);
  IMnpip_CDHz_dE_woK0_n->SetXTitle("CDHz [cm]");
  IMnpip_CDHz_dE_woK0_n->SetYTitle("IM(n#pi^{+} [GeV/c^{2}]");

  IMnpim_CDHphi_dE_woK0_n = new TH2F("IMnpim_CDHphi_dE_woK0_n","IMnpim_CDHphi_dE_woK0_n",100,-3.5,3.5,nbinIMnpi,1.,2.);
  IMnpim_CDHphi_dE_woK0_n->SetXTitle("CDHphi [radian]");
  IMnpim_CDHphi_dE_woK0_n->SetYTitle("IM(n#pi^{-} [GeV/c^{2}]");
  
  IMnpim_CDHz_dE_woK0_n = new TH2F("IMnpim_CDHz_dE_woK0_n","IMnpim_CDHz_dE_woK0_n",100,-50,50,nbinIMnpi,1.,2.);
  IMnpim_CDHz_dE_woK0_n->SetXTitle("CDHz [cm]");
  IMnpim_CDHz_dE_woK0_n->SetYTitle("IM(n#pi^{-} [GeV/c^{2}]");
  
  IMnpip_DCApip_dE_woK0_n = new TH2F("IMnpip_DCApip_dE_woK0_n","IMnpip_DCApip_dE_woK0_n",300,0,20,nbinIMnpi,1.,2.);
  IMnpip_DCApip_dE_woK0_n->SetXTitle("DCA #pi^{+} [cm]");
  IMnpip_DCApip_dE_woK0_n->SetYTitle("IM(n#pi^{+} [GeV/c^{2}]");
  
  IMnpip_DCApim_dE_woK0_n = new TH2F("IMnpip_DCApim_dE_woK0_n","IMnpip_DCApim_dE_woK0_n",300,0,20,nbinIMnpi,1.,2.);
  IMnpip_DCApim_dE_woK0_n->SetXTitle("DCA #pi^{-} [cm]");
  IMnpip_DCApim_dE_woK0_n->SetYTitle("IM(n#pi^{+} [GeV/c^{2}]");
 
  IMnpim_DCApip_dE_woK0_n = new TH2F("IMnpim_DCApip_dE_woK0_n","IMnpim_DCApip_dE_woK0_n",300,0,20,nbinIMnpi,1.,2.);
  IMnpim_DCApip_dE_woK0_n->SetXTitle("DCA #pi^{+} [cm]");
  IMnpim_DCApip_dE_woK0_n->SetYTitle("IM(n#pi^{-} [GeV/c^{2}]");
  
  IMnpim_DCApim_dE_woK0_n = new TH2F("IMnpim_DCApim_dE_woK0_n","IMnpim_DCApim_dE_woK0_n",300,0,20,nbinIMnpi,1.,2.);
  IMnpim_DCApim_dE_woK0_n->SetXTitle("DCA #pi^{-} [cm]");
  IMnpim_DCApim_dE_woK0_n->SetYTitle("IM(n#pi^{-} [GeV/c^{2}]");
  
  for(int igap=0;igap<ngap;igap++){
    IMnpim_IMnpip_dE_n_side[igap] = new TH2F(Form("IMnpim_IMnpip_dE_n_side_%d",igap),Form("IMnpim_IMnpip_dE_n_side_%d",igap),nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
    IMnpim_IMnpip_dE_n_side[igap]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
    IMnpim_IMnpip_dE_n_side[igap]->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
    
    IMnpim_IMnpip_dE_woK0_n_side[igap] = new TH2F(Form("IMnpim_IMnpip_dE_woK0_n_side_%d",igap),Form("IMnpim_IMnpip_dE_woK0_n_side_%d",igap),nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
    IMnpim_IMnpip_dE_woK0_n_side[igap]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
    IMnpim_IMnpip_dE_woK0_n_side[igap]->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  }

  for(int i=0;i<2;i++){
    for(int igap=0;igap<ngap;igap++){
      const char  lh[][6]={"low","high"};
      IMnpim_IMnpip_dE_n_Sp_side[i][igap] = new TH2F(Form("IMnpim_IMnpip_dE_n_Sp_side_%s_%d",lh[i],igap),Form("IMnpim_IMnpip_dE_n_Sp_side_%s_%d",lh[i],igap), nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
      IMnpim_IMnpip_dE_n_Sp_side[i][igap]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
      IMnpim_IMnpip_dE_n_Sp_side[i][igap]->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
    
      IMnpim_IMnpip_dE_n_Sm_side[i][igap] = new TH2F(Form("IMnpim_IMnpip_dE_n_Sm_side_%s_%d",lh[i],igap),Form("IMnpim_IMnpip_dE_n_Sm_side_%s_%d",lh[i],igap), nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
      IMnpim_IMnpip_dE_n_Sm_side[i][igap]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
      IMnpim_IMnpip_dE_n_Sm_side[i][igap]->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
 
      IMnpim_IMnpip_dE_woK0_n_Sp_side[i][igap] = new TH2F(Form("IMnpim_IMnpip_dE_woK0_n_Sp_side_%s_%d",lh[i],igap),Form("IMnpim_IMnpip_dE_woK0_n_Sp_side_%s_%d",lh[i],igap), nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
      IMnpim_IMnpip_dE_woK0_n_Sp_side[i][igap]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
      IMnpim_IMnpip_dE_woK0_n_Sp_side[i][igap]->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
    
      IMnpim_IMnpip_dE_woK0_n_Sm_side[i][igap] = new TH2F(Form("IMnpim_IMnpip_dE_woK0_n_Sm_side_%s_%d",lh[i],igap),Form("IMnpim_IMnpip_dE_woK0_n_Sm_side_%s_%d",lh[i],igap), nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
      IMnpim_IMnpip_dE_woK0_n_Sm_side[i][igap]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
      IMnpim_IMnpip_dE_woK0_n_Sm_side[i][igap]->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
    }
  }
  
  for(int izone=0;izone<nzone;izone++){
    IMnpim_IMnpip_dE_n_Sp_sidewide[izone] = new TH2F(Form("IMnpim_IMnpip_dE_n_Sp_side_wide_%d",izone),Form("IMnpim_IMnpip_dE_n_Sp_side_wide_%d",izone), nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
    IMnpim_IMnpip_dE_n_Sp_sidewide[izone]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
    IMnpim_IMnpip_dE_n_Sp_sidewide[izone]->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");

    IMnpim_IMnpip_dE_n_Sm_sidewide[izone] = new TH2F(Form("IMnpim_IMnpip_dE_n_Sm_side_wide_%d",izone),Form("IMnpim_IMnpip_dE_n_Sm_side_wide_%d",izone), nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
    IMnpim_IMnpip_dE_n_Sm_sidewide[izone]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
    IMnpim_IMnpip_dE_n_Sm_sidewide[izone]->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
    
    IMnpim_IMnpip_dE_woK0_n_Sp_sidewide[izone] = new TH2F(Form("IMnpim_IMnpip_dE_woK0_n_Sp_side_wide_%d",izone),Form("IMnpim_IMnpip_dE_woK0_n_Sp_side_wide_%d",izone), nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
    IMnpim_IMnpip_dE_woK0_n_Sp_sidewide[izone]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
    IMnpim_IMnpip_dE_woK0_n_Sp_sidewide[izone]->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");

    IMnpim_IMnpip_dE_woK0_n_Sm_sidewide[izone] = new TH2F(Form("IMnpim_IMnpip_dE_woK0_n_Sm_side_wide_%d",izone),Form("IMnpim_IMnpip_dE_woK0_n_Sm_side_wide_%d",izone), nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
    IMnpim_IMnpip_dE_woK0_n_Sm_sidewide[izone]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
    IMnpim_IMnpip_dE_woK0_n_Sm_sidewide[izone]->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  }

  IMnpim_IMnpip_dE_n_Sp = new TH2F(Form("IMnpim_IMnpip_dE_n_Sp"),Form("IMnpim_IMnpip_dE_n_Sp"), nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_n_Sp->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_n_Sp->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  
  IMnpim_IMnpip_dE_woK0_n_Sp = new TH2F(Form("IMnpim_IMnpip_dE_woK0_n_Sp"),Form("IMnpim_IMnpip_dE_woK0_n_Sp"), nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_woK0_n_Sp->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_woK0_n_Sp->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  
  IMnpim_IMnpip_dE_woK0_n_Sp_bg = new TH2F(Form("IMnpim_IMnpip_dE_woK0_n_Sp_bg"),Form("IMnpim_IMnpip_dE_woK0_n_Sp_bg"), nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_woK0_n_Sp_bg->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_woK0_n_Sp_bg->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  
  IMnpim_IMnpip_dE_n_Sm = new TH2F(Form("IMnpim_IMnpip_dE_n_Sm"),Form("IMnpim_IMnpip_dE_n_Sm"), nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_n_Sm->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_n_Sm->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  
  IMnpim_IMnpip_dE_woK0_n_Sm = new TH2F(Form("IMnpim_IMnpip_dE_woK0_n_Sm"),Form("IMnpim_IMnpip_dE_woK0_n_Sm"), nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_woK0_n_Sm->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_woK0_n_Sm->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  
  IMnpim_IMnpip_dE_woK0_n_Sm_bg = new TH2F(Form("IMnpim_IMnpip_dE_woK0_n_Sm_bg"),Form("IMnpim_IMnpip_dE_woK0_n_Sm_bg"), nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_woK0_n_Sm_bg->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_woK0_n_Sm_bg->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
    
  
  nmom_IMnpim_dE_n = new TH2F("nmom_IMnpim_dE_n","nmom_IMnpim_dE_n",nbinIMnpi,1.,2.,100,0.,1.0);
  nmom_IMnpim_dE_n->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  nmom_IMnpim_dE_n->SetYTitle("nmom. [GeV/c]");
  
  nmom_IMnpim_dE_woK0_n = new TH2F("nmom_IMnpim_dE_woK0_n","nmom_IMnpim_dE_woK0_n",nbinIMnpi,1.,2.,100,0.,1.0);
  nmom_IMnpim_dE_woK0_n->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  nmom_IMnpim_dE_woK0_n->SetYTitle("nmom. [GeV/c]");

  nmom_IMnpip_dE_n = new TH2F("nmom_IMnpip_dE_n","nmom_IMnpip_dE_n",nbinIMnpi,1.,2.,100,0.,1.0);
  nmom_IMnpip_dE_n->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  nmom_IMnpip_dE_n->SetYTitle("nmom. [GeV/c]");
  
  nmom_IMnpip_dE_woK0_n = new TH2F("nmom_IMnpip_dE_woK0_n","nmom_IMnpip_dE_woK0_n",nbinIMnpi,1.,2.,100,0.,1.0);
  nmom_IMnpip_dE_woK0_n->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  nmom_IMnpip_dE_woK0_n->SetYTitle("nmom. [GeV/c]");

  MMnpip_MMnpim_n = new TH2F(Form("MMnpip_MMnpim_n"),Form("MMnpip_MMnpim_n"),70, 1, 1.7, 70, 1, 1.7);
  MMnpip_MMnpim_n->SetXTitle("IM(n_{miss}#pi^{+}) [GeV/c^{2}]");
  MMnpip_MMnpim_n->SetYTitle("IM(n_{miss}#pi^{-}) [GeV/c^{2}]");
  
  MMnpip_MMnpim_woK0_n = new TH2F(Form("MMnpip_MMnpim_woK0_n"),Form("MMnpip_MMnpim_woK0_n"),70, 1, 1.7, 70, 1, 1.7);
  MMnpip_MMnpim_woK0_n->SetXTitle("IM(n_{miss}#pi^{+}) [GeV/c^{2}]");
  MMnpip_MMnpim_woK0_n->SetYTitle("IM(n_{miss}#pi^{-}) [GeV/c^{2}]");
  
  MMnpip_MMnpim_wSid_n = new TH2F(Form("MMnpip_MMnpim_wSid_n"),Form("MMnpip_MMnpim_wSid_n"),70, 1, 1.7, 70, 1, 1.7);
  MMnpip_MMnpim_wSid_n->SetXTitle("IM(n_{miss}#pi^{+}) [GeV/c^{2}]");
  MMnpip_MMnpim_wSid_n->SetYTitle("IM(n_{miss}#pi^{-}) [GeV/c^{2}]");
  
  MMnpip_MMnpim_woK0_wSid_n = new TH2F(Form("MMnpip_MMnpim_woK0_wSid_n"),Form("MMnpip_MMnpim_woK0_wSid_n"),70, 1, 1.7, 70, 1, 1.7);
  MMnpip_MMnpim_woK0_wSid_n->SetXTitle("IM(n_{miss}#pi^{+}) [GeV/c^{2}]");
  MMnpip_MMnpim_woK0_wSid_n->SetYTitle("IM(n_{miss}#pi^{-}) [GeV/c^{2}]");
  

  dE_IMnpim = new TH2F(Form("dE_IMnpim"),Form("dE_IMnpim"), nbinIMnpi, 1.0, 2.0, nbindE, 0, 50.);
  dE_IMnpim->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  dE_IMnpim->SetYTitle("dE [MeVee]");
  
  dE_IMnpim_woK0 = new TH2F(Form("dE_IMnpim_woK0"),Form("dE_IMnpim_woK0"), nbinIMnpi, 1.0, 2.0, nbindE, 0, 50.);
  dE_IMnpim_woK0->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  dE_IMnpim_woK0->SetYTitle("dE [MeVee]");
  
  dE_IMnpim_n = new TH2F(Form("dE_IMnpim_n"),Form("dE_IMnpim_n"), nbinIMnpi, 1.0, 2.0, nbindE, 0, 50.);
  dE_IMnpim_n->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  dE_IMnpim_n->SetYTitle("dE [MeVee]");
  
  dE_IMnpim_woK0_n = new TH2F(Form("dE_IMnpim_woK0_n"),Form("dE_IMnpim_woK0_n"), nbinIMnpi, 1.0, 2.0, nbindE, 0, 50.);
  dE_IMnpim_woK0_n->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  dE_IMnpim_woK0_n->SetYTitle("dE [MeVee]");

  dE_IMnpip = new TH2F(Form("dE_IMnpip"),Form("dE_IMnpip"), nbinIMnpi, 1.0, 2.0, nbindE, 0, 50.);
  dE_IMnpip->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  dE_IMnpip->SetYTitle("dE [MeVee]");
  
  dE_IMnpip_woK0 = new TH2F(Form("dE_IMnpip_woK0"),Form("dE_IMnpip_woK0"), nbinIMnpi, 1.0, 2.0, nbindE, 0, 50.);
  dE_IMnpip_woK0->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  dE_IMnpip_woK0->SetYTitle("dE [MeVee]");
  
  dE_IMnpip_n = new TH2F(Form("dE_IMnpip_n"),Form("dE_IMnpip_n"), nbinIMnpi, 1.0, 2.0, nbindE, 0, 50.);
  dE_IMnpip_n->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  dE_IMnpip_n->SetYTitle("dE [MeVee]");
  
  dE_IMnpip_woK0_n = new TH2F(Form("dE_IMnpip_woK0_n"),Form("dE_IMnpip_woK0_n"), nbinIMnpi, 1.0, 2.0, nbindE, 0, 50.);
  dE_IMnpip_woK0_n->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  dE_IMnpip_woK0_n->SetYTitle("dE [MeVee]");

  dE_IMnpipi_wSid_n = new TH2F(Form("dE_IMnpipi_wSid_n"),Form("dE_IMnpipi_wSid_n"),nbinIMnpipi, 1, 2, nbindE, 0, 50);
  dE_IMnpipi_wSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  dE_IMnpipi_wSid_n->SetYTitle("dE [MeVee]");
  
  dE_IMnpipi_woK0_wSid_n = new TH2F(Form("dE_IMnpipi_woK0_wSid_n"),Form("dE_IMnpipi_woK0_wSid_n"),nbinIMnpipi, 1, 2, nbindE, 0, 50);
  dE_IMnpipi_woK0_wSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  dE_IMnpipi_woK0_wSid_n->SetYTitle("dE [MeVee]");
    
  Cosn_IMnpipi_wSid_n = new TH2F(Form("Cosn_IMnpipi_wSid_n"),Form("dE_Cosn_IMnpipi_wSid_n"),100, 1, 2, 50, -1, 1);
  Cosn_IMnpipi_wSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  Cosn_IMnpipi_wSid_n->SetYTitle("cos#theta_{n} (CM)");
  
  Cosn_IMnpipi_woK0_wSid_n = new TH2F(Form("Cosn_IMnpipi_woK0_wSid_n"),Form("dE_Cosn_IMnpipi_woK0_wSid_n"),100, 1, 2, 50, -1, 1);
  Cosn_IMnpipi_woK0_wSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  Cosn_IMnpipi_woK0_wSid_n->SetYTitle("cos#theta_{n} (CM)");
    
  MMnmiss_IMnpipi_wSid_n = new TH2F(Form("MMnmiss_IMnpipi_wSid_n"),Form("MMnmiss_IMnpipi_wSid_n"),nbinIMnpipi,1,2,140,0.4,1.8);
  MMnmiss_IMnpipi_wSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMnpipi_wSid_n->SetYTitle("Miss. Mass. [GeV/c^{2}]");
  
  MMnmiss_IMnpipi_woK0_wSid = new TH2F(Form("MMnmiss_IMnpipi_woK0_wSid"),Form("MMnmiss_IMnpipi_woK0_wSid"),nbinIMnpipi,1,2,140,0.4,1.8);
  MMnmiss_IMnpipi_woK0_wSid->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMnpipi_woK0_wSid->SetYTitle("Miss. Mass. [GeV/c^{2}]");
  
  MMnmiss_IMnpipi_woK0_wSid_Sp = new TH2F(Form("MMnmiss_IMnpipi_woK0_wSid_Sp"),Form("MMnmiss_IMnpipi_woK0_wSid_Sp"),nbinIMnpipi,1,2,140,0.4,1.8);
  MMnmiss_IMnpipi_woK0_wSid_Sp->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMnpipi_woK0_wSid_Sp->SetYTitle("Miss. Mass. [GeV/c^{2}]");
  
  MMnmiss_IMnpipi_woK0_wSid_Sm = new TH2F(Form("MMnmiss_IMnpipi_woK0_wSid_Sm"),Form("MMnmiss_IMnpipi_woK0_wSid_Sm"),nbinIMnpipi,1,2,140,0.4,1.8);
  MMnmiss_IMnpipi_woK0_wSid_Sm->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMnpipi_woK0_wSid_Sm->SetYTitle("Miss. Mass. [GeV/c^{2}]");
    
  q_IMnpipi_wSid_n = new TH2F(Form("q_IMnpipi_wSid_n"),Form("q_IMnpipi_wSid_n"),nbinIMnpipi,1,2,nbinq,0,1.5);
  q_IMnpipi_wSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wSid_n->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_IMnpipi_woK0_wSid_n = new TH2F(Form("q_IMnpipi_woK0_wSid_n"),Form("q_IMnpipi_woK0_wSid_n"), nbinIMnpipi,1,2, nbinq,0,1.5);
  q_IMnpipi_woK0_wSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_wSid_n->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_IMpiSigma_gen = new TH2F(Form("q_IMpiSigma_gen"),Form("q_IMpiSigma_gen"),nbinpippim,1,2,300,0,1.5);
  if(SimSpmode){
    q_IMpiSigma_gen->SetXTitle("true IM(#Sigma^{+}#pi^{-}) [GeV/c^{2}]");
  }else if(SimSmmode){
    q_IMpiSigma_gen->SetXTitle("true IM(#Sigma^{-}#pi^{+}) [GeV/c^{2}]");
  }
  q_IMpiSigma_gen->SetYTitle("true Mom. Transfer [GeV/c]");
  
  q_IMpiSigma_wSid_n_genacc = new TH2F(Form("q_IMpiSigma_wSid_n_genacc"),Form("q_IMpiSigma_wSid_n_genacc"),nbinpippim,1,2,300,0,1.5);
  if(SimSpmode){
    q_IMpiSigma_wSid_n_genacc->SetXTitle("true IM(#Sigma^{+}#pi^{-}) [GeV/c^{2}]");
  }else if(SimSmmode){
    q_IMpiSigma_wSid_n_genacc->SetXTitle("true IM(#Sigma^{-}#pi^{+}) [GeV/c^{2}]");
  }
  q_IMpiSigma_wSid_n_genacc->SetYTitle("true Mom. Transfer [GeV/c]");

  q_IMnpipi_wSid_n_acc = new TH2F(Form("q_IMnpipi_wSid_n_acc"),Form("q_IMnpipi_wSid_n_acc"),nbinpippim,1,2,300,0,1.5);
  q_IMnpipi_wSid_n_acc->SetXTitle("true IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wSid_n_acc->SetYTitle("true Mom. Transfer [GeV/c]");
  
  q_IMnpipi_wSid_n_acc_reco = new TH2F(Form("q_IMnpipi_wSid_n_acc_reco"),Form("q_IMnpipi_wSid_n_acc_reco"),nbinpippim,1,2,300,0,1.5);
  q_IMnpipi_wSid_n_acc_reco->SetXTitle("reco. IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wSid_n_acc_reco->SetYTitle("reco. Mom. Transfer [GeV/c]");
  
  q_IMpiSigma_woK0_wSid_n_genacc = new TH2F(Form("q_IMpiSigma_woK0_wSid_n_genacc"),Form("q_IMpiSigma_woK0_wSid_n_genacc"),nbinpippim,1,2,300,0,1.5);
  if(SimSpmode){
    q_IMpiSigma_woK0_wSid_n_genacc->SetXTitle("true IM(#Sigma^{+}#pi^{-}) [GeV/c^{2}]");
  }else if(SimSmmode){
    q_IMpiSigma_woK0_wSid_n_genacc->SetXTitle("true IM(#Sigma^{-}#pi^{+}) [GeV/c^{2}]");
  }
  q_IMpiSigma_woK0_wSid_n_genacc->SetYTitle("true Mom. Transfer [GeV/c]");
  
  q_IMnpipi_woK0_wSid_n_acc = new TH2F(Form("q_IMnpipi_woK0_wSid_n_acc"),Form("q_IMnpipi_woK0_wSid_n_acc"),nbinpippim,1,2,300,0,1.5);
  q_IMnpipi_woK0_wSid_n_acc->SetXTitle("true IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_wSid_n_acc->SetYTitle("true Mom. Transfer [GeV/c]");
  
  q_IMnpipi_woK0_wSid_n_acc_reco = new TH2F(Form("q_IMnpipi_woK0_wSid_n_acc_reco"),Form("q_IMnpipi_woK0_wSid_n_acc_reco"),nbinpippim,1,2,300,0,1.5);
  q_IMnpipi_woK0_wSid_n_acc_reco->SetXTitle("reco. IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_wSid_n_acc_reco->SetYTitle("reco. Mom. Transfer [GeV/c]");
  
  q_IMnpipi_wSid_n_Sp = new TH2F(Form("q_IMnpipi_wSid_n_Sp"),Form("q_IMnpipi_wSid_n_Sp"), nbinIMnpipi,1,2, nbinq,0,1.5);
  q_IMnpipi_wSid_n_Sp->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wSid_n_Sp->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_IMnpipi_wK0_n = new TH2F(Form("q_IMnpipi_wK0_n"),Form("q_IMnpipi_wK0_n"), nbinIMnpipi,1,2, nbinq,0,1.5);
  q_IMnpipi_wK0_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wK0_n->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_IMnpipi_wK0_wSid_n_Sp = new TH2F(Form("q_IMnpipi_wK0_wSid_n_Sp"),Form("q_IMnpipi_wK0_wSid_n_Sp"), nbinIMnpipi,1,2, nbinq,0,1.5);
  q_IMnpipi_wK0_wSid_n_Sp->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wK0_wSid_n_Sp->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_IMnpipi_woK0_wSid_n_Sp = new TH2F(Form("q_IMnpipi_woK0_wSid_n_Sp"),Form("q_IMnpipi_woK0_wSid_n_Sp"), nbinIMnpipi,1,2, nbinq,0,1.5);
  q_IMnpipi_woK0_wSid_n_Sp->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_wSid_n_Sp->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_IMpiSigma_wSid_n_Sp_genacc = new TH2F(Form("q_IMpiSigma_wSid_n_Sp_genacc"),Form("q_IMpiSigma_wSid_n_Sp_genacc"),nbinpippim,1,2,300,0,1.5);
  if(SimSpmode){
    q_IMpiSigma_wSid_n_Sp_genacc->SetXTitle("true IM(#Sigma^{+}#pi^{-}) [GeV/c^{2}]");
  }else if(SimSmmode){
    q_IMpiSigma_wSid_n_Sp_genacc->SetXTitle("true IM(#Sigma^{-}#pi^{+}) [GeV/c^{2}]");
  }

  q_IMpiSigma_woK0_wSid_n_Sp_genacc = new TH2F(Form("q_IMpiSigma_woK0_wSid_n_Sp_genacc"),Form("q_IMpiSigma_woK0_wSid_n_Sp_genacc"),nbinpippim,1,2,300,0,1.5);
  if(SimSpmode){
    q_IMpiSigma_woK0_wSid_n_Sp_genacc->SetXTitle("true IM(#Sigma^{+}#pi^{-}) [GeV/c^{2}]");
  }else if(SimSmmode){
    q_IMpiSigma_woK0_wSid_n_Sp_genacc->SetXTitle("true IM(#Sigma^{-}#pi^{+}) [GeV/c^{2}]");
  }
  q_IMpiSigma_woK0_wSid_n_Sp_genacc->SetYTitle("true Mom. Transfer [GeV/c]");
  
  q_IMnpipi_wSid_n_Sp_acc = new TH2F(Form("q_IMnpipi_wSid_n_Sp_acc"),Form("q_IMnpipi_wSid_n_Sp_acc"),nbinpippim,1,2,300,0,1.5);
  q_IMnpipi_wSid_n_Sp_acc->SetXTitle("true IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wSid_n_Sp_acc->SetYTitle("true Mom. Transfer [GeV/c]");
  
  q_IMnpipi_woK0_wSid_n_Sp_acc = new TH2F(Form("q_IMnpipi_woK0_wSid_n_Sp_acc"),Form("q_IMnpipi_woK0_wSid_n_Sp_acc"),nbinpippim,1,2,300,0,1.5);
  q_IMnpipi_woK0_wSid_n_Sp_acc->SetXTitle("true IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_wSid_n_Sp_acc->SetYTitle("true Mom. Transfer [GeV/c]");
  
  q_IMnpipi_wSid_n_Sp_acc_reco = new TH2F(Form("q_IMnpipi_wSid_n_Sp_acc_reco"),Form("q_IMnpipi_wSid_n_Sp_acc_reco"),nbinpippim,1,2,300,0,1.5);
  q_IMnpipi_wSid_n_Sp_acc_reco->SetXTitle("reco. IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wSid_n_Sp_acc_reco->SetYTitle("reco. Mom. Transfer [GeV/c]");
  
  q_IMnpipi_woK0_wSid_n_Sp_acc_reco = new TH2F(Form("q_IMnpipi_woK0_wSid_n_Sp_acc_reco"),Form("q_IMnpipi_woK0_wSid_n_Sp_acc_reco"),nbinpippim,1,2,300,0,1.5);
  q_IMnpipi_woK0_wSid_n_Sp_acc_reco->SetXTitle("reco. IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_wSid_n_Sp_acc_reco->SetYTitle("reco. Mom. Transfer [GeV/c]");
  
  for(int itype=0;itype<3;itype++){
    for(int ilh=0;ilh<2;ilh++){
      for(int igap=0;igap<ngap;igap++){
        const char  lh[][6]={"low","high"};
        q_IMnpipi_wSid_n_Sp_side[itype][ilh][igap] = new TH2F(Form("q_IMnpipi_wSid_n_Sp_side_%d_%s_%d",itype,lh[ilh],igap),Form("q_IMnpipi_wSid_n_Sp_side_%d_%s_%d",itype,lh[ilh],igap), nbinIMnpipi,1,2, nbinq,0,1.5);
        q_IMnpipi_wSid_n_Sp_side[itype][ilh][igap]->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
        q_IMnpipi_wSid_n_Sp_side[itype][ilh][igap]->SetYTitle("Mom. Transfer [GeV/c]");
        
        q_IMnpipi_woK0_wSid_n_Sp_side[itype][ilh][igap] = new TH2F(Form("q_IMnpipi_woK0_wSid_n_Sp_side_%d_%s_%d",itype,lh[ilh],igap),Form("q_IMnpipi_woK0_wSid_n_Sp_side_%d_%s_%d",itype,lh[ilh],igap), nbinIMnpipi,1,2, nbinq,0,1.5);
        q_IMnpipi_woK0_wSid_n_Sp_side[itype][ilh][igap]->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
        q_IMnpipi_woK0_wSid_n_Sp_side[itype][ilh][igap]->SetYTitle("Mom. Transfer [GeV/c]");
      }
    }
  }

  for(int izone=0;izone<nzone;izone++){
    q_IMnpipi_wSid_n_Sp_sidewide[izone] = new TH2F(Form("q_IMnpipi_wSid_n_Sp_sidewide_%d",izone),Form("q_IMnpipi_wSid_n_Sp_sidewide_%d",izone), nbinIMnpipi,1,2, nbinq,0,1.5);
    q_IMnpipi_wSid_n_Sp_sidewide[izone]->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
    q_IMnpipi_wSid_n_Sp_sidewide[izone]->SetYTitle("Mom. Transfer [GeV/c]");
    q_IMnpipi_woK0_wSid_n_Sp_sidewide[izone] = new TH2F(Form("q_IMnpipi_woK0_wSid_n_Sp_sidewide_%d",izone),Form("q_IMnpipi_woK0_wSid_n_Sp_sidewide_%d",izone), nbinIMnpipi,1,2, nbinq,0,1.5);
    q_IMnpipi_woK0_wSid_n_Sp_sidewide[izone]->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
    q_IMnpipi_woK0_wSid_n_Sp_sidewide[izone]->SetYTitle("Mom. Transfer [GeV/c]");
    q_IMnpipi_wSid_n_Sm_sidewide[izone] = new TH2F(Form("q_IMnpipi_wSid_n_Sm_sidewide_%d",izone),Form("q_IMnpipi_wSid_n_Sm_sidewide_%d",izone), nbinIMnpipi,1,2, nbinq,0,1.5);
    q_IMnpipi_wSid_n_Sm_sidewide[izone]->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
    q_IMnpipi_wSid_n_Sm_sidewide[izone]->SetYTitle("Mom. Transfer [GeV/c]");
    q_IMnpipi_woK0_wSid_n_Sm_sidewide[izone] = new TH2F(Form("q_IMnpipi_woK0_wSid_n_Sm_sidewide_%d",izone),Form("q_IMnpipi_woK0_wSid_n_Sm_sidewide_%d",izone), nbinIMnpipi,1,2, nbinq,0,1.5);
    q_IMnpipi_woK0_wSid_n_Sm_sidewide[izone]->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
    q_IMnpipi_woK0_wSid_n_Sm_sidewide[izone]->SetYTitle("Mom. Transfer [GeV/c]");
  }

  q_IMnpipi_wSid_n_Sm = new TH2F(Form("q_IMnpipi_wSid_n_Sm"),Form("q_IMnpipi_wSid_n_Sm"), nbinIMnpipi,1,2, nbinq,0,1.5);
  q_IMnpipi_wSid_n_Sm->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wSid_n_Sm->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_IMnpipi_wK0_wSid_n_Sm = new TH2F(Form("q_IMnpipi_wK0_wSid_n_Sm"),Form("q_IMnpipi_wK0_wSid_n_Sm"), nbinIMnpipi,1,2, nbinq,0,1.5);
  q_IMnpipi_wK0_wSid_n_Sm->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wK0_wSid_n_Sm->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_IMnpipi_woK0_wSid_n_Sm = new TH2F(Form("q_IMnpipi_woK0_wSid_n_Sm"),Form("q_IMnpipi_woK0_wSid_n_Sm"), nbinIMnpipi,1,2, nbinq,0,1.5);
  q_IMnpipi_woK0_wSid_n_Sm->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_wSid_n_Sm->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_IMpiSigma_wSid_n_Sm_genacc = new TH2F(Form("q_IMpiSigma_wSid_n_Sm_genacc"),Form("q_IMpiSigma_wSid_n_Sm_genacc"),nbinpippim,1,2,300,0,1.5);
  if(SimSpmode){
    q_IMpiSigma_wSid_n_Sm_genacc->SetXTitle("true IM(#Sigma^{+}#pi^{-}) [GeV/c^{2}]");
  }else if(SimSmmode){
    q_IMpiSigma_wSid_n_Sm_genacc->SetXTitle("true IM(#Sigma^{-}#pi^{+}) [GeV/c^{2}]");
  }
  q_IMpiSigma_wSid_n_Sm_genacc->SetYTitle("true Mom. Transfer [GeV/c]");
  
  q_IMpiSigma_woK0_wSid_n_Sm_genacc = new TH2F(Form("q_IMpiSigma_woK0_wSid_n_Sm_genacc"),Form("q_IMpiSigma_woK0_wSid_n_Sm_genacc"),nbinpippim,1,2,300,0,1.5);
  if(SimSpmode){
    q_IMpiSigma_woK0_wSid_n_Sm_genacc->SetXTitle("true IM(#Sigma^{+}#pi^{-}) [GeV/c^{2}]");
  }else if(SimSmmode){
    q_IMpiSigma_woK0_wSid_n_Sm_genacc->SetXTitle("true IM(#Sigma^{-}#pi^{+}) [GeV/c^{2}]");
  }
  q_IMpiSigma_woK0_wSid_n_Sm_genacc->SetYTitle("true Mom. Transfer [GeV/c]");
  
  q_IMnpipi_wSid_n_Sm_acc = new TH2F(Form("q_IMnpipi_wSid_n_Sm_acc"),Form("q_IMnpipi_wSid_n_Sm_acc"),nbinpippim,1,2,300,0,1.5);
  q_IMnpipi_wSid_n_Sm_acc->SetXTitle("true IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wSid_n_Sm_acc->SetYTitle("true Mom. Transfer [GeV/c]");
  
  q_IMnpipi_woK0_wSid_n_Sm_acc = new TH2F(Form("q_IMnpipi_woK0_wSid_n_Sm_acc"),Form("q_IMnpipi_woK0_wSid_n_Sm_acc"),nbinpippim,1,2,300,0,1.5);
  q_IMnpipi_woK0_wSid_n_Sm_acc->SetXTitle("true IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_wSid_n_Sm_acc->SetYTitle("true Mom. Transfer [GeV/c]");
  
  q_IMnpipi_wSid_n_Sm_acc_reco = new TH2F(Form("q_IMnpipi_wSid_n_Sm_acc_reco"),Form("q_IMnpipi_wSid_n_Sm_acc_reco"),nbinpippim,1,2,300,0,1.5);
  q_IMnpipi_wSid_n_Sm_acc_reco->SetXTitle("reco. IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wSid_n_Sm_acc_reco->SetYTitle("reco. Mom. Transfer [GeV/c]");
  
  q_IMnpipi_woK0_wSid_n_Sm_acc_reco = new TH2F(Form("q_IMnpipi_woK0_wSid_n_Sm_acc_reco"),Form("q_IMnpipi_woK0_wSid_n_Sm_acc_reco"),nbinpippim,1,2,300,0,1.5);
  q_IMnpipi_woK0_wSid_n_Sm_acc_reco->SetXTitle("reco. IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_wSid_n_Sm_acc_reco->SetYTitle("reco. Mom. Transfer [GeV/c]");
  
  for(int itype=0;itype<3;itype++){
    for(int ilh=0;ilh<2;ilh++){
      for(int igap=0;igap<ngap;igap++){
        const char  lh[][6]={"low","high"};
        q_IMnpipi_wSid_n_Sm_side[itype][ilh][igap] = new TH2F(Form("q_IMnpipi_wSid_n_Sm_side_%d_%s_%d",itype,lh[ilh],igap),Form("q_IMnpipi_wSid_n_Sm_side_%d_%s_%d",itype,lh[ilh],igap),nbinIMnpipi,1,2,nbinq,0,1.5);
        q_IMnpipi_wSid_n_Sm_side[itype][ilh][igap]->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
        q_IMnpipi_wSid_n_Sm_side[itype][ilh][igap]->SetYTitle("Mom. Transfer [GeV/c]");
        q_IMnpipi_woK0_wSid_n_Sm_side[itype][ilh][igap] = new TH2F(Form("q_IMnpipi_woK0_wSid_n_Sm_side_%d_%s_%d",itype,lh[ilh],igap),Form("q_IMnpipi_woK0_wSid_n_Sm_side_%d_%s_%d",itype,lh[ilh],igap),nbinIMnpipi,1,2,nbinq,0,1.5);
        q_IMnpipi_woK0_wSid_n_Sm_side[itype][ilh][igap]->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
        q_IMnpipi_woK0_wSid_n_Sm_side[itype][ilh][igap]->SetYTitle("Mom. Transfer [GeV/c]");
      }
    }
  }
  
  for(int igap=0;igap<ngap;igap++){
    q_IMnpipi_wSid_n_side[igap] = new TH2F(Form("q_IMnpipi_wSid_n_side_%d",igap),Form("q_IMnpipi_wSid_n_side_%d",igap), nbinIMnpipi,1,2, nbinq,0,1.5);
    q_IMnpipi_wSid_n_side[igap]->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
    q_IMnpipi_wSid_n_side[igap]->SetYTitle("Mom. Transfer [GeV/c]");
  }

  for(int igap=0;igap<ngap;igap++){
    q_IMnpipi_woK0_wSid_n_side[igap] = new TH2F(Form("q_IMnpipi_woK0_wSid_n_side_%d",igap),Form("q_IMnpipi_woK0_wSid_n_side_%d",igap), nbinIMnpipi,1,2, nbinq,0,1.5);
    q_IMnpipi_woK0_wSid_n_side[igap]->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
    q_IMnpipi_woK0_wSid_n_side[igap]->SetYTitle("Mom. Transfer [GeV/c]");
  }
   

  IMnpip_IMnpipi_n = new TH2F(Form("IMnpip_IMnpipi_n"),Form("IMnpip_IMnpipi_n"),nbinIMnpipi,1.0,2.0, nbinIMnpi,1.0,2.0);
  IMnpip_IMnpipi_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  IMnpip_IMnpipi_n->SetYTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  
  IMnpip_IMnpipi_wK0_n = new TH2F(Form("IMnpip_IMnpipi_wK0_n"),Form("IMnpip_IMnpipi_wK0_n"),nbinIMnpipi,1.0,2.0, nbinIMnpi,1.0,2.0);
  IMnpip_IMnpipi_wK0_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  IMnpip_IMnpipi_wK0_n->SetYTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  
  IMnpip_IMnpipi_woK0_n = new TH2F(Form("IMnpip_IMnpipi_woK0_n"),Form("IMnpip_IMnpipi_woK0_n"),nbinIMnpipi,1.0,2.0, nbinIMnpi,1.0,2.0);
  IMnpip_IMnpipi_woK0_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  IMnpip_IMnpipi_woK0_n->SetYTitle("IM(n#pi^{+}) [GeV/c^{2}]");

  IMnpim_IMnpipi_n = new TH2F(Form("IMnpim_IMnpipi_n"),Form("IMnpim_IMnpipi_n"),nbinIMnpipi,1.0,2.0, nbinIMnpi,1.0,2.0);
  IMnpim_IMnpipi_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  IMnpim_IMnpipi_n->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  
  IMnpim_IMnpipi_wK0_n = new TH2F(Form("IMnpim_IMnpipi_wK0_n"),Form("IMnpim_IMnpipi_wK0_n"),nbinIMnpipi,1.0,2.0, nbinIMnpi,1.0,2.0);
  IMnpim_IMnpipi_wK0_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  IMnpim_IMnpipi_wK0_n->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  
  IMnpim_IMnpipi_woK0_n = new TH2F(Form("IMnpim_IMnpipi_woK0_n"),Form("IMnpim_IMnpipi_woK0_n"),nbinIMnpipi,1.0,2.0, nbinIMnpi,1.0,2.0);
  IMnpim_IMnpipi_woK0_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  IMnpim_IMnpipi_woK0_n->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  
  nmom_IMnpipi_wK0_n = new TH2F(Form("nmom_IMnpipi_wK0_n"),Form("nmom_IMnpipi_wK0_n"), nbinIMnpipi,1,2,100,0,1.0);
  nmom_IMnpipi_wK0_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  nmom_IMnpipi_wK0_n->SetYTitle("nmom  [GeV/c]");
  
  nmom_IMnpipi_wSid_n = new TH2F(Form("nmom_IMnpipi_wSid_n"),Form("nmom_IMnpipi_wSid_n"), nbinIMnpipi,1,2,100,0,1.0);
  nmom_IMnpipi_wSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  nmom_IMnpipi_wSid_n->SetYTitle("nmom  [GeV/c]");
  
  nmom_IMnpipi_woK0_wSid_n = new TH2F(Form("nmom_IMnpipi_woK0_wSid_n"),Form("nmom_IMnpipi_woK0_wSid_n"), nbinIMnpipi,1,2,100,0,1.0);
  nmom_IMnpipi_woK0_wSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  nmom_IMnpipi_woK0_wSid_n->SetYTitle("nmom  [GeV/c]");
  
  nmom_IMnpipi_woK0_wSid_n_Sp = new TH2F(Form("nmom_IMnpipi_woK0_wSid_n_Sp"),Form("nmom_IMnpipi_woK0_wSid_n_Sp"), nbinIMnpipi,1,2,100,0,1.0);
  nmom_IMnpipi_woK0_wSid_n_Sp->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  nmom_IMnpipi_woK0_wSid_n_Sp->SetYTitle("nmom  [GeV/c]");
   
  nmom_IMnpipi_woK0_wSid_n_Sm = new TH2F(Form("nmom_IMnpipi_woK0_wSid_n_Sm"),Form("nmom_IMnpipi_woK0_wSid_n_Sm"), nbinIMnpipi,1,2,100,0,1.0);
  nmom_IMnpipi_woK0_wSid_n_Sm->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  nmom_IMnpipi_woK0_wSid_n_Sm->SetYTitle("nmom  [GeV/c]");
   
  nmom_cosnlab_K0_n = new TH2F(Form("nmom_cosnlab_K0_n"),Form("nmom_cosnlab_K0_n"), 200,-1,1,100,0,1.0);
  nmom_cosnlab_K0_n->SetXTitle("cos_n (lab.)");
  nmom_cosnlab_K0_n->SetYTitle("nmom  [GeV/c]");
  
  nmom_IMpippim_n = new TH2F(Form("nmom_IMpippim_n"),Form("nmom_IMpippim_n"), nbinpippim,0,1,100,0,1.0);
  nmom_IMpippim_n->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  nmom_IMpippim_n->SetYTitle("nmom  [GeV/c]");
  
  nmom_IMpippim_wSid_n = new TH2F(Form("nmom_IMpippim_wSid_n"),Form("nmom_IMpippim_wSid_n"), nbinpippim,0,1,100,0,1.0);
  nmom_IMpippim_wSid_n->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  nmom_IMpippim_wSid_n->SetYTitle("nmom  [GeV/c]");
  
  mnmom_IMpippim_n = new TH2F(Form("mnmom_IMpippim_n"),Form("mnmom_IMpippim_n"), nbinpippim,0,1,200,0,2.0);
  mnmom_IMpippim_n->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  mnmom_IMpippim_n->SetYTitle("miss. nmom  [GeV/c]");
  
  mnmom_IMpippim_wSid_n = new TH2F(Form("mnmom_IMpippim_wSid_n"),Form("mnmom_IMpippim_wSid_n"), nbinpippim,0,1,200,0,2.0);
  mnmom_IMpippim_wSid_n->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  mnmom_IMpippim_wSid_n->SetYTitle("miss. nmom  [GeV/c]");
  
  q_IMpippim_n = new TH2F("q_IMpippim_n","q_IMpippim_n",nbinpippim,0,1, 100,0,1.5);
  q_IMpippim_n->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMpippim_n->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_IMpippim_n_wSid = new TH2F("q_IMpippim_n_wSid","q_IMpippim_n_wSid",nbinpippim,0,1, 100,0,1.5);
  q_IMpippim_n_wSid->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMpippim_n_wSid->SetYTitle("Mom. Transfer [GeV/c]");
  
  IMpippim_IMnpipi_n = new TH2F("IMpippim_IMnpipi_n","IMpippim_IMnpipi_n",100,1,2,nbinpippim,0,1);
  IMpippim_IMnpipi_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  IMpippim_IMnpipi_n->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  
  IMpippim_IMnpipi_n_wSid = new TH2F("IMpippim_IMnpipi_n_wSid","IMpippim_IMnpipi_n_wSid",100,1,2,nbinpippim,0,1);
  IMpippim_IMnpipi_n_wSid->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  IMpippim_IMnpipi_n_wSid->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");

  IMpippim_IMnpip_n = new TH2F("IMpippim_IMnpip_n","IMpippim_IMnpip_n",nbinIMnpi,1,2,nbinpippim,0,1);
  IMpippim_IMnpip_n->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMpippim_IMnpip_n->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  
  IMpippim_IMnpim_n = new TH2F("IMpippim_IMnpim_n","IMpippim_IMnpim_n",nbinIMnpi,1,2,nbinpippim,0,1);
  IMpippim_IMnpim_n->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  IMpippim_IMnpim_n->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  
  Mompippim_IMnpipi_dE_wK0_n = new TH2F("Mompippim_IMnpipi_dE_wK0_n","Mompippim_IMnpipi_dE_wK0_n",nbinIMnpipi,1,2,50,0,1);
  Mompippim_IMnpipi_dE_wK0_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  Mompippim_IMnpipi_dE_wK0_n->SetYTitle("Mom(#pi^{+}#pi^{-}) [GeV/c]");
  
  Mompippim_IMnpipi_dE_wK0_n_Sp = new TH2F("Mompippim_IMnpipi_dE_wK0_n_Sp","Mompippim_IMnpipi_dE_wK0_n_Sp",nbinIMnpipi,1,2,50,0,1);
  Mompippim_IMnpipi_dE_wK0_n_Sp->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  Mompippim_IMnpipi_dE_wK0_n_Sp->SetYTitle("Mom(#pi^{+}#pi^{-}) [GeV/c]");
  
  Mompippim_IMnpipi_dE_wK0_n_Sm = new TH2F("Mompippim_IMnpipi_dE_wK0_n_Sm","Mompippim_IMnpipi_dE_wK0_n_Sm",nbinIMnpipi,1,2,50,0,1);
  Mompippim_IMnpipi_dE_wK0_n_Sm->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  Mompippim_IMnpipi_dE_wK0_n_Sm->SetYTitle("Mom(#pi^{+}#pi^{-}) [GeV/c]");

  //
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
  TH2F* q_IMnpipi_wSid_n_mc;
  q_IMnpipi_wSid_n_mc = new TH2F(Form("q_IMnpipi_wSid_n_mc"),Form("q_IMnpipi_wSid_n_mc"),nbinIMnpipi,1,2,nbinq,0,1.5);
  q_IMnpipi_wSid_n_mc->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wSid_n_mc->SetYTitle("Mom. Transfer [GeV/c]");
  
  TH2F* q_IMnpipi_woK0_wSid_n_mc;
  q_IMnpipi_woK0_wSid_n_mc = new TH2F(Form("q_IMnpipi_woK0_wSid_n_mc"),Form("q_IMnpipi_woK0_wSid_n_mc"),nbinIMnpipi,1,2,nbinq,0,1.5);
  q_IMnpipi_woK0_wSid_n_mc->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_wSid_n_mc->SetYTitle("Mom. Transfer [GeV/c]");
  
  TH2F* q_IMnpipi_wSid_n_Sp_mc;
  q_IMnpipi_wSid_n_Sp_mc = new TH2F(Form("q_IMnpipi_wSid_n_Sp_mc"),Form("q_IMnpipi_wSid_n_Sp_mc"),nbinIMnpipi,1,2,nbinq,0,1.5);
  q_IMnpipi_wSid_n_Sp_mc->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wSid_n_Sp_mc->SetYTitle("Mom. Transfer [GeV/c]");
  
  TH2F* q_IMnpipi_woK0_wSid_n_Sp_mc;
  q_IMnpipi_woK0_wSid_n_Sp_mc = new TH2F(Form("q_IMnpipi_woK0_wSid_n_Sp_mc"),Form("q_IMnpipi_woK0_wSid_n_Sp_mc"),nbinIMnpipi,1,2,nbinq,0,1.5);
  q_IMnpipi_woK0_wSid_n_Sp_mc->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_wSid_n_Sp_mc->SetYTitle("Mom. Transfer [GeV/c]");

  TH2F* q_IMnpipi_wSid_n_Sm_mc;
  q_IMnpipi_wSid_n_Sm_mc = new TH2F(Form("q_IMnpipi_wSid_n_Sm_mc"),Form("q_IMnpipi_wSid_n_Sm_mc"),nbinIMnpipi,1,2,nbinq,0,1.5);
  q_IMnpipi_wSid_n_Sm_mc->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wSid_n_Sm_mc->SetYTitle("Mom. Transfer [GeV/c]");
  
  TH2F* q_IMnpipi_woK0_wSid_n_Sm_mc;
  q_IMnpipi_woK0_wSid_n_Sm_mc = new TH2F(Form("q_IMnpipi_woK0_wSid_n_Sm_mc"),Form("q_IMnpipi_woK0_wSid_n_Sm_mc"),nbinIMnpipi,1,2,nbinq,0,1.5);
  q_IMnpipi_woK0_wSid_n_Sm_mc->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_wSid_n_Sm_mc->SetYTitle("Mom. Transfer [GeV/c]");
  
  
  TH2F* diff_IMnpipi_wSid_n;
  diff_IMnpipi_wSid_n = new TH2F(Form("diff_IMnpipi_wSid_n"),Form("diff_IMnpipi_wSid_n"),nbinIMnpipi,1,2,100,-1,1);
  diff_IMnpipi_wSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  diff_IMnpipi_wSid_n->SetYTitle("reco. - gen.  [GeV/c^{2}]");
  
  TH2F* diff_IMnpipi_woK0_wSid_n;
  diff_IMnpipi_woK0_wSid_n = new TH2F(Form("diff_IMnpipi_woK0_wSid_n"),Form("diff_IMnpipi_woK0_wSid_n"),nbinIMnpipi,1,2,100,-1,1);
  diff_IMnpipi_woK0_wSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  diff_IMnpipi_woK0_wSid_n->SetYTitle("reco. - gen.  [GeV/c^{2}]");

  TH2F* diff_IMnpipi_wSid_n_Sp;
  diff_IMnpipi_wSid_n_Sp = new TH2F(Form("diff_IMnpipi_wSid_n_Sp"),Form("diff_IMnpipi_wSid_n_Sp"),nbinIMnpipi,1,2,600,-0.3,0.3);
  diff_IMnpipi_wSid_n_Sp->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  diff_IMnpipi_wSid_n_Sp->SetYTitle("reco. - gen. [GeV/c^{2}]");
  
  TH2F* diff_IMnpipi_woK0_wSid_n_Sp;
  diff_IMnpipi_woK0_wSid_n_Sp = new TH2F(Form("diff_IMnpipi_woK0_wSid_n_Sp"),Form("diff_IMnpipi_woK0_wSid_n_Sp"),nbinIMnpipi,1,2,600,-0.3,0.3);
  diff_IMnpipi_woK0_wSid_n_Sp->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  diff_IMnpipi_woK0_wSid_n_Sp->SetYTitle("reco. - gen. [GeV/c^{2}]");
  
  TH2F* diff_q_wSid_n_Sp;
  diff_q_wSid_n_Sp = new TH2F(Form("diff_q_wSid_n_Sp"),Form("diff_q_wSid_n_Sp"),nbinq,0,1,600,-0.3,0.3);
  diff_q_wSid_n_Sp->SetXTitle("Mom. Transfer [GeV/c]");
  diff_q_wSid_n_Sp->SetYTitle("reco. - gen. [GeV/c^]");

  TH2F* diff_q_woK0_wSid_n_Sp;
  diff_q_woK0_wSid_n_Sp = new TH2F(Form("diff_q_woK0_wSid_n_Sp"),Form("diff_q_woK0_wSid_n_Sp"),nbinq,0,1,600,-0.3,0.3);
  diff_q_woK0_wSid_n_Sp->SetXTitle("Mom. Transfer [GeV/c]");
  diff_q_woK0_wSid_n_Sp->SetYTitle("reco. - gen. [GeV/c^]");

  TH2F* diff_IMnpipi_wSid_n_Sm;
  diff_IMnpipi_wSid_n_Sm = new TH2F(Form("diff_IMnpipi_wSid_n_Sm"),Form("diff_IMnpipi_wSid_n_Sm"),nbinIMnpipi,1,2,600,-0.3,0.3);
  diff_IMnpipi_wSid_n_Sm->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  diff_IMnpipi_wSid_n_Sm->SetYTitle("reco. - gen.  [GeV/c^{2}]");
  
  TH2F* diff_IMnpipi_woK0_wSid_n_Sm;
  diff_IMnpipi_woK0_wSid_n_Sm = new TH2F(Form("diff_IMnpipi_woK0_wSid_n_Sm"),Form("diff_IMnpipi_woK0_wSid_n_Sm"),nbinIMnpipi,1,2,600,-0.3,0.3);
  diff_IMnpipi_woK0_wSid_n_Sm->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  diff_IMnpipi_woK0_wSid_n_Sm->SetYTitle("reco. - gen.  [GeV/c^{2}]");

  TH2F* diff_q_wSid_n_Sm;
  diff_q_wSid_n_Sm = new TH2F(Form("diff_q_wSid_n_Sm"),Form("diff_q_wSid_n_Sm"),nbinq,0,1,600,-0.3,0.3);
  diff_q_wSid_n_Sm->SetXTitle("Mom. Transfer [GeV/c]");
  diff_q_wSid_n_Sm->SetYTitle("reco. - gen. [GeV/c^]");
  
  TH2F* diff_q_woK0_wSid_n_Sm;
  diff_q_woK0_wSid_n_Sm = new TH2F(Form("diff_q_woK0_wSid_n_Sm"),Form("diff_q_woK0_wSid_n_Sm"),nbinq,0,1,600,-0.3,0.3);
  diff_q_woK0_wSid_n_Sm->SetXTitle("Mom. Transfer [GeV/c]");
  diff_q_woK0_wSid_n_Sm->SetYTitle("reco. - gen. [GeV/c^]");
  
  //reading TTree
  Int_t nevent = tree->GetEntries();
  std::cerr<<"# of events = "<<nevent<<std::endl;
  std::cout << "p-value cut:" << pvalcut << std::endl; 
  std::cout << "dE cut:" << anacuts::dE_MIN << std::endl; 
  TCanvas *cinfo = new TCanvas("cinfo","info");
  cinfo->cd();
  TPaveText *pt = new TPaveText(.05,.05,.95,.7);
  if(qvalcutflag == 1) {
    pt->AddText(Form("Mom. Transfer  < %0.2f ",anacuts::qvalcut));
  }else if(qvalcutflag == 2){
    pt->AddText(Form("Mom. Transfer  > %0.2f ",anacuts::qvalcut));
  }
  pt->AddText(Form("p-value cut: %f ",pvalcut));
  pt->AddText(Form("dE cut: %0.2f " ,anacuts::dE_MIN));
  pt->AddText(Form("1/beta min.: %f ",1./anacuts::beta_MAX));
  pt->AddText(Form("1/beta max : %f ",1./anacuts::beta_MIN));
  pt->AddText(Form("K^{0} window : %0.3f - %0.3f",anacuts::pipi_MIN,anacuts::pipi_MAX )); 
  pt->AddText(Form("#Sigma^{+} window : %0.3f - %0.3f (%0.1f sigma cut)",anacuts::Sigmap_MIN,anacuts::Sigmap_MAX,anacuts::Sigma_NSigmacut)); 
  pt->AddText(Form("#Sigma^{-} window : %0.3f - %0.3f (%0.1f sigma cut)",anacuts::Sigmam_MIN,anacuts::Sigmam_MAX,anacuts::Sigma_NSigmacut)); 
  pt->AddText(Form("miss. n window : %0.3f - %0.3f (%0.1f sigma cut)",anacuts::neutron_MIN,anacuts::neutron_MAX,anacuts::neutron_NSigmacut )); 
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
    TLorentzVector LVec_nmiss = *LVec_target+*LVec_beam-*LVec_pip-*LVec_pim-*LVec_n;
    TLorentzVector LVec_nmiss_vtx[2];
    if(UseKinFit){
      if(!UseKinFitVal){
        LVec_nmiss_vtx[0] = *LVec_target+*LVec_beam_Sp-*LVec_pip-*LVec_pim-*LVec_n_Sp;
        LVec_nmiss_vtx[1] = *LVec_target+*LVec_beam_Sm-*LVec_pip-*LVec_pim-*LVec_n_Sm;
      }else{
        LVec_nmiss_vtx[0] = *LVec_target+*kfSpmode_mom_beam-*kfSpmode_mom_pip-*kfSpmode_mom_pim-*kfSpmode_mom_n;
        LVec_nmiss_vtx[1] = *LVec_target+*kfSmmode_mom_beam-*kfSmmode_mom_pip-*kfSmmode_mom_pim-*kfSmmode_mom_n;
      }
    }
    double nmiss_mass = LVec_nmiss.M();
    double nmiss_mass_vtx[2]={LVec_nmiss_vtx[0].M(),LVec_nmiss_vtx[1].M()};
    double nmiss_mom = LVec_nmiss.P();

    // calc cos(theta) of missing n //
    TVector3 boost = (*LVec_target+*LVec_beam).BoostVector();
    TVector3 boost_vtx[2] = {(*LVec_target+*LVec_beam_Sp).BoostVector(),(*LVec_target+*LVec_beam_Sm).BoostVector()} ;
    TLorentzVector LVec_nmiss_CM = LVec_nmiss;
    TLorentzVector LVec_beam_CM = *LVec_beam;
    LVec_nmiss_CM.Boost(-boost);
    LVec_beam_CM.Boost(-boost);
    double cos_nmiss = LVec_nmiss_CM.Vect().Dot(LVec_beam_CM.Vect())/(LVec_nmiss_CM.Vect().Mag()*LVec_beam_CM.Vect().Mag());
    //cos(theta) of n_cds in lab
    double cos_nlab = LVec_n->Vect().Dot(LVec_beam->Vect())/(LVec_n->Vect().Mag()*LVec_beam->Vect().Mag());
    if(SimSpmode || SimSmmode){
      TVector3 boost_mc =  (*LVec_target+*mcmom_beam).BoostVector();
    }
    TLorentzVector qkn_mc;
    if(SimSpmode || SimSmmode){
      qkn_mc = *mcmom_beam-*mcmom_nmiss;
    }
    TLorentzVector qkn_vtx[2];
    if(UseKinFit){
      if(!UseKinFitVal){ 
        qkn_vtx[0] = *LVec_beam_Sp-LVec_nmiss_vtx[0];
        qkn_vtx[1] = *LVec_beam_Sm-LVec_nmiss_vtx[1];
      }else{
        qkn_vtx[0] = *kfSpmode_mom_beam-LVec_nmiss_vtx[0];
        qkn_vtx[1] = *kfSmmode_mom_beam-LVec_nmiss_vtx[1];
      }
    }
    // calc pi+pi- //
    TLorentzVector LVec_pip_pim = *LVec_pip+*LVec_pim;
    TLorentzVector LVec_pip_pim_mc;
    if(SimSpmode || SimSmmode){
      LVec_pip_pim_mc = *mcmom_pip+*mcmom_pim; 
    }
    // calc pi+n //
    TLorentzVector LVec_pip_n = *LVec_pip+*LVec_n;
    TLorentzVector LVec_pip_n_mc;
    if(SimSpmode || SimSmmode){
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
    if(SimSpmode || SimSmmode){
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
    if(SimSpmode || SimSmmode){
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
    if(SimSpmode || SimSmmode){
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
    if(SimSpmode || SimSmmode){
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
    TLorentzVector qkn = *LVec_beam-LVec_nmiss;
    TLorentzVector LVec_pip_pim_n_CM = LVec_pip_pim_n;
    LVec_pip_pim_n_CM.Boost(-boost);
    //double cos_X = LVec_pip_pim_n_CM.Vect().Dot(LVec_beam_CM.Vect())/(LVec_pip_pim_n_CM.Vect().Mag()*LVec_beam_CM.Vect().Mag());
    
    if( (qkn.P()>=anacuts::qvalcut) && (qvalcutflag==1) ) continue;
    if( (qkn.P()<anacuts::qvalcut) && (qvalcutflag==2) ) continue;
    if((*LVec_n).P()<0.10 ) continue;
    //if(LVec_pip_pim_n.M() < 1.45) continue;
    //if(LVec_pip_pim_n.M() > 1.55) continue;
    //if(dcapippim < 1) continue;
    //if(LVec_pip_pim_n.M()<1.45 ) continue;
    //double chi2 = kfSpmode_chi2<kfSmmode_chi2 ? kfSpmode_chi2:kfSmmode_chi2;
    double pvalue = -9999;
    if(UseKinFit) pvalue = kfSmmode_pvalue<kfSpmode_pvalue ? kfSpmode_pvalue:kfSmmode_pvalue;
    
    //Filling generated info.
    
    TLorentzVector LVec_Sigma_react;
    TLorentzVector LVec_pi_react;
    TLorentzVector LVec_piSigma_react;
    TLorentzVector qkn_react;
    if(SimSpmode || SimSmmode){
      LVec_Sigma_react = *react_Sigma;
      LVec_pi_react = *react_pi;
      LVec_piSigma_react = LVec_Sigma_react + LVec_pi_react;
      TLorentzVector LVec_beam_r = *LVec_beam;
      double px = (*LVec_beam).Px();
      double py = (*LVec_beam).Py();
      double pz = (*LVec_beam).Pz();
      double E = (*LVec_beam).E();
      TLorentzVector LVec_beam_unit;
      LVec_beam_unit.SetPx(px*1000.0);
      LVec_beam_unit.SetPy(py*1000.0);
      LVec_beam_unit.SetPz(pz*1000.0);
      LVec_beam_unit.SetE(E*1000.0);
      //std::cout << LVec_beam_r.P() << std::endl;
      //std::cout << LVec_beam_r.M() << std::endl;
      //std::cout << "beam p " << LVec_beam_unit.P() << std::endl;
      //std::cout << "beam m " << LVec_beam_unit.M() << std::endl;
      //std::cout << "nmiss P" <<  (*react_nmiss).P() << std::endl;
      qkn_react = LVec_beam_unit-*react_nmiss ;
      q_IMpiSigma_gen->Fill(LVec_piSigma_react.M()/1000.,qkn_react.P()/1000.);
    }

    bool K0rejectFlag=false;
    bool MissNFlag=false;
    bool NBetaOK=false;
    bool NdEOK=false;
    bool SigmaPFlag=false;
    bool SigmaMFlag=false;
    bool SigmawidePFlag=false;
    bool SigmawideMFlag=false;
    bool SigmaPsideFlag[3][ngap];//pattern, gap
    for(int ipat=0;ipat<3;ipat++){
      for(int igap=0;igap<ngap;igap++){
        SigmaPsideFlag[ipat][igap]=false;
      }
    }
    bool SigmaPsideLowFlag[ngap];
    bool SigmaPsideHighFlag[ngap];
    for(int igap=0;igap<ngap;igap++){
      SigmaPsideLowFlag[igap]=false;
      SigmaPsideHighFlag[igap]=false;
    }

    bool SigmaMsideFlag[3][ngap];//pattern,gap
    for(int ipat=0;ipat<3;ipat++){
      for(int igap=0;igap<ngap;igap++){
        SigmaMsideFlag[ipat][igap]=false;
      }
    }
    bool SigmaMsideLowFlag[ngap];
    bool SigmaMsideHighFlag[ngap];
    for(int igap=0;igap<ngap;igap++){
      SigmaMsideHighFlag[igap]=false;
      SigmaMsideLowFlag[igap]=false;
    }
    
    bool SigmaPsideWideFlag[nzone];
    bool SigmaMsideWideFlag[nzone];
    for(int izone=0;izone<nzone;izone++){
      SigmaPsideWideFlag[izone]=false;
      SigmaMsideWideFlag[izone]=false;
    }

    
    //triangular cuts to minimize acceptance distortion
    bool SigmaCrossFlagTop=false;
    bool SigmaCrossFlagBottom=false;
    bool SigmaCrossFlagLeft=false;
    bool SigmaCrossFlagRight=false;
    
    //triangular cuts for side band of Sp mode low mass side
    bool SigmaCrossPsideLowFlagTop[ngap];
    bool SigmaCrossPsideLowFlagBottom[ngap];
    bool SigmaCrossPsideLowFlagLeft[ngap];
    bool SigmaCrossPsideLowFlagRight[ngap];
    
    //triangular cuts for side band of Sp mode high mass side
    bool SigmaCrossPsideHighFlagTop[ngap];
    bool SigmaCrossPsideHighFlagBottom[ngap];
    bool SigmaCrossPsideHighFlagLeft[ngap];
    bool SigmaCrossPsideHighFlagRight[ngap];

    //triangular cuts for side band of Sm mode low mass side
    bool SigmaCrossMsideLowFlagTop[ngap];
    bool SigmaCrossMsideLowFlagBottom[ngap];
    bool SigmaCrossMsideLowFlagLeft[ngap];
    bool SigmaCrossMsideLowFlagRight[ngap];
    
    //triangular cuts for side band of Sm mode high mass side
    bool SigmaCrossMsideHighFlagTop[ngap];
    bool SigmaCrossMsideHighFlagBottom[ngap];
    bool SigmaCrossMsideHighFlagLeft[ngap];
    bool SigmaCrossMsideHighFlagRight[ngap];
    
    bool SigmaCrossPsideWideFlagTop[nzone];
    bool SigmaCrossPsideWideFlagBottom[nzone];
    bool SigmaCrossPsideWideFlagLeft[nzone];
    bool SigmaCrossPsideWideFlagRight[nzone];
    
    bool SigmaCrossMsideWideFlagTop[nzone];
    bool SigmaCrossMsideWideFlagBottom[nzone];
    bool SigmaCrossMsideWideFlagLeft[nzone];
    bool SigmaCrossMsideWideFlagRight[nzone];

    for(int igap=0;igap<ngap;igap++){
      SigmaCrossPsideLowFlagTop[igap] = false;   
      SigmaCrossPsideLowFlagBottom[igap] = false;
      SigmaCrossPsideLowFlagLeft[igap] = false;  
      SigmaCrossPsideLowFlagRight[igap]  = false; 

      SigmaCrossPsideHighFlagTop[igap] = false;     
      SigmaCrossPsideHighFlagBottom[igap] = false;  
      SigmaCrossPsideHighFlagLeft[igap] = false;    
      SigmaCrossPsideHighFlagRight[igap] = false;      
      
      SigmaCrossMsideLowFlagTop[igap] = false;    
      SigmaCrossMsideLowFlagBottom[igap] = false; 
      SigmaCrossMsideLowFlagLeft[igap] = false;   
      SigmaCrossMsideLowFlagRight[igap] = false;  
      
      SigmaCrossMsideHighFlagTop[igap] = false;     
      SigmaCrossMsideHighFlagBottom[igap] = false;  
      SigmaCrossMsideHighFlagLeft[igap] = false;    
      SigmaCrossMsideHighFlagRight[igap] = false;
    }

    for(int izone=0;izone<nzone;izone++){
      SigmaCrossPsideWideFlagTop[izone] = false;   
      SigmaCrossPsideWideFlagBottom[izone] = false;
      SigmaCrossPsideWideFlagLeft[izone] = false;  
      SigmaCrossPsideWideFlagRight[izone] = false; 
                                    
      SigmaCrossMsideWideFlagTop[izone] = false;   
      SigmaCrossMsideWideFlagBottom[izone] = false;
      SigmaCrossMsideWideFlagLeft[izone] =false;  
      SigmaCrossMsideWideFlagRight[izone] = false; 
    }

  
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

    const double IMnpibinWidth = 1./nbinIMnpi;//2.5 MeV 
    //Sigma cross side band Sp low mass top
    for(int igap=0;igap<ngap;igap++){
      if(  (MassNPim >  (MassNPip - (anacuts::Sigmap_sidelow_center-igap*IMnpibinWidth)) + anacuts::Sigmam_center)  
          && (MassNPim > -(MassNPip - (anacuts::Sigmap_sidelow_center-igap*IMnpibinWidth)) + anacuts::Sigmam_center)) SigmaCrossPsideLowFlagTop[igap]=true;

      //Sigma cross side band Sp low mass bottom
      if(  (MassNPim <  (MassNPip - (anacuts::Sigmap_sidelow_center-igap*IMnpibinWidth)) + anacuts::Sigmam_center)  
          && (MassNPim < -(MassNPip - (anacuts::Sigmap_sidelow_center-igap*IMnpibinWidth)) + anacuts::Sigmam_center)) SigmaCrossPsideLowFlagBottom[igap]=true;

      //Sigma cross side band Sp low mass right 
      if(  (MassNPim <  (MassNPip - (anacuts::Sigmap_sidelow_center-igap*IMnpibinWidth)) + anacuts::Sigmam_center)  
          && (MassNPim > -(MassNPip - (anacuts::Sigmap_sidelow_center-igap*IMnpibinWidth)) + anacuts::Sigmam_center)) SigmaCrossPsideLowFlagRight[igap]=true;

      //Sigma cross side band Sp low mass left
      if(  (MassNPim >  (MassNPip - (anacuts::Sigmap_sidelow_center-igap*IMnpibinWidth)) + anacuts::Sigmam_center)  
          && (MassNPim < -(MassNPip - (anacuts::Sigmap_sidelow_center-igap*IMnpibinWidth)) + anacuts::Sigmam_center)) SigmaCrossPsideLowFlagLeft[igap]=true;


      //Sigma cross side band Sp high mass top
      if(  (MassNPim >  (MassNPip - (anacuts::Sigmap_sidehigh_center+igap*IMnpibinWidth)) + anacuts::Sigmam_center)  
          && (MassNPim > -(MassNPip - (anacuts::Sigmap_sidehigh_center+igap*IMnpibinWidth)) + anacuts::Sigmam_center)) SigmaCrossPsideHighFlagTop[igap]=true;

      //Sigma cross side band Sp high mass bottom
      if(  (MassNPim <  (MassNPip - (anacuts::Sigmap_sidehigh_center+igap*IMnpibinWidth)) + anacuts::Sigmam_center)  
          && (MassNPim < -(MassNPip - (anacuts::Sigmap_sidehigh_center+igap*IMnpibinWidth)) + anacuts::Sigmam_center)) SigmaCrossPsideHighFlagBottom[igap]=true;

      //Sigma cross side band Sp high mass right 
      if(  (MassNPim <  (MassNPip - (anacuts::Sigmap_sidehigh_center+igap*IMnpibinWidth)) + anacuts::Sigmam_center)  
          && (MassNPim > -(MassNPip - (anacuts::Sigmap_sidehigh_center+igap*IMnpibinWidth)) + anacuts::Sigmam_center)) SigmaCrossPsideHighFlagRight[igap]=true;

      //Sigma cross side band Sp high mass left
      if(  (MassNPim >  (MassNPip - (anacuts::Sigmap_sidehigh_center+igap*IMnpibinWidth)) + anacuts::Sigmam_center)  
          && (MassNPim < -(MassNPip - (anacuts::Sigmap_sidehigh_center+igap*IMnpibinWidth)) + anacuts::Sigmam_center)) SigmaCrossPsideHighFlagLeft[igap]=true;


      //Sigma cross side band Sm low mass top
      if(  (MassNPim >  (MassNPip - anacuts::Sigmap_center) + (anacuts::Sigmam_sidelow_center-igap*IMnpibinWidth))  
          && (MassNPim > -(MassNPip - anacuts::Sigmap_center) + (anacuts::Sigmam_sidelow_center-igap*IMnpibinWidth))) SigmaCrossMsideLowFlagTop[igap]=true;

      //Sigma cross side band Sp low mass bottom
      if(  (MassNPim <  (MassNPip - anacuts::Sigmap_center) + (anacuts::Sigmam_sidelow_center-igap*IMnpibinWidth))  
          && (MassNPim < -(MassNPip - anacuts::Sigmap_center) + (anacuts::Sigmam_sidelow_center-igap*IMnpibinWidth))) SigmaCrossMsideLowFlagBottom[igap]=true;

      //Sigma cross side band Sp low mass right 
      if(  (MassNPim <  (MassNPip - anacuts::Sigmap_center) + (anacuts::Sigmam_sidelow_center-igap*IMnpibinWidth))  
          && (MassNPim > -(MassNPip - anacuts::Sigmap_center) + (anacuts::Sigmam_sidelow_center-igap*IMnpibinWidth))) SigmaCrossMsideLowFlagRight[igap]=true;

      //Sigma cross side band Sp low mass left
      if(  (MassNPim >  (MassNPip - anacuts::Sigmap_center) + (anacuts::Sigmam_sidelow_center-igap*IMnpibinWidth))  
          && (MassNPim < -(MassNPip - anacuts::Sigmap_center) + (anacuts::Sigmam_sidelow_center-igap*IMnpibinWidth))) SigmaCrossMsideLowFlagLeft[igap]=true;


      //Sigma cross side band Sm high mass top
      if(  (MassNPim >  (MassNPip - anacuts::Sigmap_center) + (anacuts::Sigmam_sidehigh_center+igap*IMnpibinWidth))  
          && (MassNPim > -(MassNPip - anacuts::Sigmap_center) + (anacuts::Sigmam_sidehigh_center+igap*IMnpibinWidth))) SigmaCrossMsideHighFlagTop[igap]=true;

      //Sigma cross side band Sm high mass bottom
      if(  (MassNPim <  (MassNPip - anacuts::Sigmap_center) + (anacuts::Sigmam_sidehigh_center+igap*IMnpibinWidth))  
          && (MassNPim < -(MassNPip - anacuts::Sigmap_center) + (anacuts::Sigmam_sidehigh_center+igap*IMnpibinWidth))) SigmaCrossMsideHighFlagBottom[igap]=true;

      //Sigma cross side band Sm high mass right 
      if(  (MassNPim <  (MassNPip - anacuts::Sigmap_center) + (anacuts::Sigmam_sidehigh_center+igap*IMnpibinWidth))  
          && (MassNPim > -(MassNPip - anacuts::Sigmap_center) + (anacuts::Sigmam_sidehigh_center+igap*IMnpibinWidth))) SigmaCrossMsideHighFlagRight[igap]=true;

      //Sigma cross side band Sm high mass left
      if(  (MassNPim >  (MassNPip - anacuts::Sigmap_center) + (anacuts::Sigmam_sidehigh_center+igap*IMnpibinWidth)) 
          && (MassNPim < -(MassNPip - anacuts::Sigmap_center) + (anacuts::Sigmam_sidehigh_center+igap*IMnpibinWidth))) SigmaCrossMsideHighFlagLeft[igap]=true;


      //Sigma+ production side band low mass side
      if( ((anacuts::Sigmap_sidelow_MIN-igap*IMnpibinWidth)<MassNPip) && 
            (MassNPip < (anacuts::Sigmap_sidelow_MAX-igap*IMnpibinWidth))) SigmaPsideLowFlag[igap]=true;

      //Sigma+ production side band high mass side
      if( ((anacuts::Sigmap_sidehigh_MIN+igap*IMnpibinWidth)<MassNPip) && 
            (MassNPip < (anacuts::Sigmap_sidehigh_MAX+igap*IMnpibinWidth))) SigmaPsideHighFlag[igap]=true;

      //Sigma+ production side band low or high mass side
      //
      //type 0: diagonal cut
      if( (SigmaPsideLowFlag[igap]  && !SigmaCrossPsideLowFlagLeft[igap]  && !SigmaCrossPsideLowFlagRight[igap]) 
          ||  (SigmaPsideHighFlag[igap] && !SigmaCrossPsideHighFlagLeft[igap]  && !SigmaCrossPsideHighFlagRight[igap])
        ) SigmaPsideFlag[0][igap]=true;

      //type 1: rectangle cut, avoiding signal region
      if( (SigmaPsideLowFlag[igap] || SigmaPsideHighFlag[igap]) && !SigmaMFlag) SigmaPsideFlag[1][igap]=true;

      //type 2: rectangle cut, avoiding signal region + 2 sigma  
      if( (SigmaPsideLowFlag[igap] || SigmaPsideHighFlag[igap]) && !SigmawideMFlag) SigmaPsideFlag[2][igap]=true;

      //Sigma- production side band low mass side
      if( ((anacuts::Sigmam_sidelow_MIN-igap*IMnpibinWidth)<MassNPim) && 
            (MassNPim < (anacuts::Sigmam_sidelow_MAX-igap*IMnpibinWidth))) SigmaMsideLowFlag[igap]=true;

      //Sigma- production side band high mass side
      if( ((anacuts::Sigmam_sidehigh_MIN+igap*IMnpibinWidth)<MassNPim) && 
           (MassNPim <  (anacuts::Sigmam_sidehigh_MAX+igap*IMnpibinWidth))) SigmaMsideHighFlag[igap]=true;
    
      //Sigma- production side band low or high mass side
      if( (SigmaMsideLowFlag[igap]  && !SigmaCrossMsideLowFlagTop[igap] && !SigmaCrossMsideLowFlagBottom[igap]) 
          ||  (SigmaMsideHighFlag[igap] && !SigmaCrossMsideHighFlagTop[igap] && !SigmaCrossMsideHighFlagBottom[igap])
        ) SigmaMsideFlag[0][igap]=true;

    
      //type 1: rectangle cut, avoiding signal region
      if( (SigmaMsideLowFlag[igap] || SigmaMsideHighFlag[igap]) &&  !SigmaPFlag) SigmaMsideFlag[1][igap]=true;

      //type 2: rectangle cut, avoiding signal region + 2 sigma
      if( (SigmaMsideLowFlag[igap] || SigmaMsideHighFlag[igap]) &&  !SigmawidePFlag) SigmaMsideFlag[2][igap]=true;

    }//for ngap

    
    for(int izone=0;izone<nzone;izone++){
      //Sigma + mode
      if( (MassNPip>anacuts::Sigmap_MIN+(int(izone-nzone/2.))*IMnpibinWidth)
        && (MassNPip<anacuts::Sigmap_MAX+(int(izone-nzone/2.))*IMnpibinWidth)) SigmaPsideWideFlag[izone] = true;
          

      if(  (MassNPim >  (MassNPip - (anacuts::Sigmap_center+(izone-nzone/2)*IMnpibinWidth)) + anacuts::Sigmam_center)  
          && (MassNPim > -(MassNPip - (anacuts::Sigmap_center+(izone-nzone/2)*IMnpibinWidth)) + anacuts::Sigmam_center)) SigmaCrossPsideWideFlagTop[izone]=true;

      if(  (MassNPim <  (MassNPip - (anacuts::Sigmap_center+(izone-nzone/2)*IMnpibinWidth)) + anacuts::Sigmam_center)  
          && (MassNPim < -(MassNPip - (anacuts::Sigmap_center+(izone-nzone/2)*IMnpibinWidth)) + anacuts::Sigmam_center)) SigmaCrossPsideWideFlagBottom[izone]=true;

      if(  (MassNPim <  (MassNPip - (anacuts::Sigmap_center+(izone-nzone/2)*IMnpibinWidth)) + anacuts::Sigmam_center)  
          && (MassNPim > -(MassNPip - (anacuts::Sigmap_center+(izone-nzone/2)*IMnpibinWidth)) + anacuts::Sigmam_center)) SigmaCrossPsideWideFlagRight[izone]=true;

      if(  (MassNPim >  (MassNPip - (anacuts::Sigmap_center+(izone-nzone/2)*IMnpibinWidth)) + anacuts::Sigmam_center)  
          && (MassNPim < -(MassNPip - (anacuts::Sigmap_center+(izone-nzone/2)*IMnpibinWidth)) + anacuts::Sigmam_center)) SigmaCrossPsideWideFlagLeft[izone]=true;
      
      //Sigma - mode
      if( (MassNPim>anacuts::Sigmam_MIN+(izone-nzone/2)*IMnpibinWidth)
        && (MassNPim<anacuts::Sigmam_MAX+(izone-nzone/2)*IMnpibinWidth)) SigmaMsideWideFlag[izone] = true;
      
      if(  (MassNPim >  (MassNPip - anacuts::Sigmap_center) + anacuts::Sigmam_center+(izone-nzone/2)*IMnpibinWidth)  
          && (MassNPim > -(MassNPip - anacuts::Sigmap_center) + anacuts::Sigmam_center+(izone-nzone/2)*IMnpibinWidth)) SigmaCrossMsideWideFlagTop[izone]=true;

      if(  (MassNPim <  (MassNPip - anacuts::Sigmap_center) + anacuts::Sigmam_center+(izone-nzone/2)*IMnpibinWidth)  
          && (MassNPim < -(MassNPip - anacuts::Sigmap_center) + anacuts::Sigmam_center+(izone-nzone/2)*IMnpibinWidth)) SigmaCrossMsideWideFlagBottom[izone]=true;

      if(  (MassNPim <  (MassNPip - anacuts::Sigmap_center) + anacuts::Sigmam_center+(izone-nzone/2)*IMnpibinWidth)  
          && (MassNPim > -(MassNPip - anacuts::Sigmap_center) + anacuts::Sigmam_center+(izone-nzone/2)*IMnpibinWidth)) SigmaCrossMsideWideFlagRight[izone]=true;

      if(  (MassNPim >  (MassNPip - anacuts::Sigmap_center) + anacuts::Sigmam_center+(izone-nzone/2)*IMnpibinWidth)  
          && (MassNPim < -(MassNPip - anacuts::Sigmap_center) + anacuts::Sigmam_center+(izone-nzone/2)*IMnpibinWidth)) SigmaCrossMsideWideFlagLeft[izone]=true;
    
    }//izone


    if(anacuts::neutron_MIN<nmiss_mass && nmiss_mass<anacuts::neutron_MAX ) MissNFlag=true;
    
    //K0 rejection using original momentum
    if( (LVec_pip_pim.M()<anacuts::pipi_MIN || anacuts::pipi_MAX<LVec_pip_pim.M())) K0rejectFlag=true;
    //---end of Flag definition-----------------------------------------------------

    //w/o kinfit
    //---including K0 --------------------------------------------------------------
    dE_betainv_fid->Fill(1./NeutralBetaCDH,dE);
    if(NBetaOK){
      dE_nmom_fid_beta->Fill((*LVec_n).P(),dE);
      dE_MMom_fid_beta->Fill(LVec_nmiss.P(),dE);
      dE_MMass_fid_beta->Fill(LVec_nmiss.M(),dE);
      if(SigmaPFlag || SigmaMFlag){
        dE_MMass_fid_beta_wSid->Fill(LVec_nmiss.M(),dE);
      }
      dE_IMnpim->Fill(LVec_pim_n.M(),dE);
      dE_IMnpip->Fill(LVec_pip_n.M(),dE);
    }
    if(NBetaOK && MissNFlag){
      dE_IMnpim_n->Fill(LVec_pim_n.M(),dE);
      dE_IMnpip_n->Fill(LVec_pip_n.M(),dE);
    }
    if(NBetaOK && NdEOK){
      MMom_MMass->Fill(LVec_nmiss.M(),LVec_nmiss.P());
      IMnpim_IMnpip_dE->Fill(LVec_pip_n.M(),LVec_pim_n.M());
      MMnmiss_IMnpip_dE->Fill(LVec_pip_n.M(),nmiss_mass);
      MMnmiss_IMnpim_dE->Fill(LVec_pim_n.M(),nmiss_mass);
      MMnmiss_IMpippim_dE->Fill(LVec_pip_pim.M(),nmiss_mass);
      if(SigmaPFlag || SigmaMFlag){
        MMom_MMass_wSid->Fill(LVec_nmiss.M(),LVec_nmiss.P());
      }
      MMnmiss_IMnpipi_wSid_n->Fill(LVec_pip_pim_n.M(), nmiss_mass);
    }

    if(NBetaOK && NdEOK && MissNFlag){
      IMnpim_IMnpip_dE_n->Fill(LVec_pip_n.M(),LVec_pim_n.M());
      IMnpip_IMnpipi_n->Fill(LVec_pip_pim_n.M(),LVec_pip_n.M());
      IMnpim_IMnpipi_n->Fill(LVec_pip_pim_n.M(),LVec_pim_n.M());

      for(int igap=0;igap<ngap;igap++){
        if(SigmaPsideFlag[sidebandtype][igap] || SigmaMsideFlag[sidebandtype][igap]){
          IMnpim_IMnpip_dE_n_side[igap]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
        }
        if(SigmaPsideLowFlag[igap] && SigmaPsideFlag[sidebandtype][igap]){
          IMnpim_IMnpip_dE_n_Sp_side[LOWside][igap]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
        }
        if(SigmaPsideHighFlag[igap] && SigmaPsideFlag[sidebandtype][igap]){
          IMnpim_IMnpip_dE_n_Sp_side[HIGHside][igap]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
        }
        if(SigmaMsideLowFlag[igap] && SigmaMsideFlag[sidebandtype][igap]){
          IMnpim_IMnpip_dE_n_Sm_side[LOWside][igap]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
        }
        if(SigmaMsideHighFlag[igap] && SigmaMsideFlag[sidebandtype][igap]){
          IMnpim_IMnpip_dE_n_Sm_side[HIGHside][igap]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
        }
      }

      for(int izone=0;izone<nzone;izone++){
        if(sidebandtype==0){
          if( !SigmaCrossPsideWideFlagLeft[izone] && !SigmaCrossPsideWideFlagRight[izone]  && SigmaPsideWideFlag[izone]){
            IMnpim_IMnpip_dE_n_Sp_sidewide[izone]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
            q_IMnpipi_wSid_n_Sp_sidewide[izone]->Fill(LVec_pip_pim_n.M(),qkn.P());
          }
          if( !SigmaCrossMsideWideFlagTop[izone] && !SigmaCrossMsideWideFlagBottom[izone]  && SigmaMsideWideFlag[izone]){
            IMnpim_IMnpip_dE_n_Sm_sidewide[izone]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
            q_IMnpipi_wSid_n_Sm_sidewide[izone]->Fill(LVec_pip_pim_n.M(),qkn.P());
          }
        }else if(sidebandtype==1){
          if( !SigmaMFlag  && SigmaPsideWideFlag[izone]){
            IMnpim_IMnpip_dE_n_Sp_sidewide[izone]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
            q_IMnpipi_wSid_n_Sp_sidewide[izone]->Fill(LVec_pip_pim_n.M(),qkn.P());
          }
          if( !SigmaPFlag  && SigmaMsideWideFlag[izone]){
            IMnpim_IMnpip_dE_n_Sm_sidewide[izone]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
            q_IMnpipi_wSid_n_Sm_sidewide[izone]->Fill(LVec_pip_pim_n.M(),qkn.P());
          }
        }else if(sidebandtype==2){
          if( !SigmawideMFlag  && SigmaPsideWideFlag[izone]){
            IMnpim_IMnpip_dE_n_Sp_sidewide[izone]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
            q_IMnpipi_wSid_n_Sp_sidewide[izone]->Fill(LVec_pip_pim_n.M(),qkn.P());
          }
          if( !SigmawidePFlag  && SigmaMsideWideFlag[izone]){
            IMnpim_IMnpip_dE_n_Sm_sidewide[izone]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
            q_IMnpipi_wSid_n_Sm_sidewide[izone]->Fill(LVec_pip_pim_n.M(),qkn.P());
          }
        }
      }//izone


      MMnpip_MMnpim_n->Fill(LVec_pim_n_miss.M(),LVec_pip_n_miss.M());
     
      //0: diagonal cut
      //1: 3 sigma cut
      //2: 5 simga cut 
      nmom->Fill((*LVec_n).P());
      mnmom->Fill(nmiss_mom);
      dE_nmom->Fill((*LVec_n).P(),dE);
      npipmom->Fill(LVec_pip_n.P());
      npimmom->Fill(LVec_pim_n.P());
      nmom_IMnpim_dE_n->Fill(LVec_pim_n.M(),(*LVec_n).P());
      nmom_IMnpip_dE_n->Fill(LVec_pip_n.M(),(*LVec_n).P());
      if(SigmaPcutFlag[sigmacuttype]) {
        IMnpim_IMnpip_dE_n_Sp->Fill(LVec_pip_n.M(),LVec_pim_n.M());
        q_IMnpipi_wSid_n_Sp->Fill(LVec_pip_pim_n.M(),qkn.P());
        q_IMpiSigma_wSid_n_Sp_genacc->Fill(LVec_piSigma_react.M()/1000.,qkn_react.P()/1000.);
        q_IMnpipi_wSid_n_Sp_acc->Fill(LVec_pip_pim_n_mc.M(),qkn_mc.P());
        q_IMnpipi_wSid_n_Sp_acc_reco->Fill(LVec_pip_pim_n.M(),qkn.P());
        if(SimSpmode){
          q_IMnpipi_wSid_n_Sp_mc->Fill(LVec_pip_pim_n_mc.M(),qkn_mc.P());
          diff_IMnpipi_wSid_n_Sp->Fill(LVec_pip_pim_n.M(),LVec_pip_pim_n.M()-LVec_pip_pim_n_mc.M());
          diff_q_wSid_n_Sp->Fill(qkn.P(),qkn.P()-qkn_mc.P());
        }
      }
      
      for(int itype=0;itype<3;itype++){
        for(int igap=0;igap<ngap;igap++){
          if(SigmaPsideFlag[itype][igap] && SigmaPsideLowFlag[igap] ){
            q_IMnpipi_wSid_n_Sp_side[itype][LOWside][igap]->Fill(LVec_pip_pim_n.M(),qkn.P());
          }
          if(SigmaPsideFlag[itype][igap] && SigmaPsideHighFlag[igap] ){
            q_IMnpipi_wSid_n_Sp_side[itype][HIGHside][igap]->Fill(LVec_pip_pim_n.M(),qkn.P());
          }
        }
      }
      
      if(SigmaMcutFlag[sigmacuttype]){
        IMnpim_IMnpip_dE_n_Sm->Fill(LVec_pip_n.M(),LVec_pim_n.M());
        q_IMnpipi_wSid_n_Sm->Fill(LVec_pip_pim_n.M(),qkn.P());
        q_IMpiSigma_wSid_n_Sm_genacc->Fill(LVec_piSigma_react.M()/1000.,qkn_react.P()/1000.);
        q_IMnpipi_wSid_n_Sm_acc->Fill(LVec_pip_pim_n_mc.M(),qkn_mc.P());
        q_IMnpipi_wSid_n_Sm_acc_reco->Fill(LVec_pip_pim_n.M(),qkn.P());
        if(SimSmmode){
          q_IMnpipi_wSid_n_Sm_mc->Fill(LVec_pip_pim_n_mc.M(),qkn_mc.P());
          diff_IMnpipi_wSid_n_Sm->Fill(LVec_pip_pim_n.M(),LVec_pip_pim_n.M()-LVec_pip_pim_n_mc.M());
          diff_q_wSid_n_Sm->Fill(qkn.P(),qkn.P()-qkn_mc.P());
        }
      }
      for(int igap=0;igap<ngap;igap++){
        for(int itype=0;itype<3;itype++){
          if(SigmaMsideFlag[itype][igap] && SigmaMsideLowFlag[igap]){
            q_IMnpipi_wSid_n_Sm_side[itype][LOWside][igap]->Fill(LVec_pip_pim_n.M(),qkn.P());
          }
          if(SigmaMsideFlag[itype][igap] && SigmaMsideHighFlag[igap]){
            q_IMnpipi_wSid_n_Sm_side[itype][HIGHside][igap]->Fill(LVec_pip_pim_n.M(),qkn.P());
          }
        }
        if(NBetaOK && NdEOK && MissNFlag){
          if(SigmaPsideFlag[sidebandtype][igap] || SigmaMsideFlag[sidebandtype][igap]){
            q_IMnpipi_wSid_n_side[igap]->Fill(LVec_pip_pim_n.M(),qkn.P());
          }
        }
      }
      
      q_IMpippim_n->Fill(LVec_pip_pim.M(),qkn.P());
      IMpippim_IMnpipi_n->Fill(LVec_pip_pim_n.M(), LVec_pip_pim.M());
      IMpippim_IMnpip_n->Fill(LVec_pip_n.M(),LVec_pip_pim.M());
      IMpippim_IMnpim_n->Fill(LVec_pim_n.M(),LVec_pip_pim.M());
      if(!K0rejectFlag){
        nmom_cosnlab_K0_n->Fill(cos_nlab,(*LVec_n).P());
        nmom_IMnpipi_wK0_n->Fill(LVec_pip_pim_n.M(),(*LVec_n).P());
      }
      nmom_IMpippim_n->Fill(LVec_pip_pim.M(),(*LVec_n).P());
      mnmom_IMpippim_n->Fill(LVec_pip_pim.M(),(LVec_nmiss).P());
      if(SigmaPFlag || SigmaMFlag){
        MMnpip_MMnpim_wSid_n->Fill(LVec_pim_n_miss.M(),LVec_pip_n_miss.M());
        q_IMpippim_n_wSid->Fill(LVec_pip_pim.M(),qkn.P());
        Cosn_IMnpipi_wSid_n->Fill(LVec_pip_pim_n.M(),cos_nmiss);
        dE_IMnpipi_wSid_n->Fill(LVec_pip_pim_n.M(),dE);
        IMpippim_IMnpipi_n_wSid->Fill(LVec_pip_pim_n.M(), LVec_pip_pim.M());
        nmom_IMpippim_wSid_n->Fill(LVec_pip_pim.M(),(*LVec_n).P());
        mnmom_IMpippim_wSid_n->Fill(LVec_pip_pim.M(),(LVec_nmiss).P());
        q_IMnpipi_wSid_n->Fill(LVec_pip_pim_n.M(),qkn.P());
        q_IMpiSigma_wSid_n_genacc->Fill(LVec_piSigma_react.M()/1000.,qkn_react.P()/1000.);
        q_IMnpipi_wSid_n_acc->Fill(LVec_pip_pim_n_mc.M(),qkn_mc.P());
        q_IMnpipi_wSid_n_acc_reco->Fill(LVec_pip_pim_n.M(),qkn.P());
        nmom_IMnpipi_wSid_n->Fill(LVec_pip_pim_n.M(),(*LVec_n).P());
      }
      for(int igap=0;igap<ngap;igap++){
        if(SigmaPsideFlag[sidebandtype][igap] || SigmaMsideFlag[sidebandtype][igap]){
          q_IMnpipi_wSid_n_side[igap]->Fill(LVec_pip_pim_n.M(),qkn.P());
        }
      }
    }//NBetaOK && NdEOK && MissNFlag
    //---including K0 end----------------------------------------------------------------------
    
    //selection K0
    if(!K0rejectFlag && NBetaOK && NdEOK && MissNFlag){
      Mompippim_IMnpipi_dE_wK0_n->Fill(LVec_pip_pim_n.M(),LVec_pip_pim.P());
      q_IMnpipi_wK0_n->Fill(LVec_pip_pim_n.M(),qkn.P());
      if(SigmaPcutFlag[sigmacuttype]){
        Mompippim_IMnpipi_dE_wK0_n_Sp->Fill(LVec_pip_pim_n.M(),LVec_pip_pim.P());
        q_IMnpipi_wK0_wSid_n_Sp->Fill(LVec_pip_pim_n.M(),qkn.P());
      }
      IMnpip_IMnpipi_wK0_n->Fill(LVec_pip_pim_n.M(),LVec_pip_n.M());
      IMnpim_IMnpipi_wK0_n->Fill(LVec_pip_pim_n.M(),LVec_pim_n.M());
    
      if(SigmaMcutFlag[sigmacuttype]){
        Mompippim_IMnpipi_dE_wK0_n_Sm->Fill(LVec_pip_pim_n.M(),LVec_pip_pim.P());
        q_IMnpipi_wK0_wSid_n_Sm->Fill(LVec_pip_pim_n.M(),qkn.P());
      }
    }
    //---rejecting K0--------------------------------------------------------------------------
    //K0 rejection
    if(K0rejectFlag && NBetaOK){
      dE_MMom_fid_beta_woK0->Fill(LVec_nmiss.P(),dE);
      dE_MMass_fid_beta_woK0->Fill(LVec_nmiss.M(),dE);
      if(SigmaPFlag || SigmaMFlag){
        dE_MMass_fid_beta_woK0_wSid->Fill(LVec_nmiss.M(),dE);
      }
      dE_IMnpim_woK0->Fill(LVec_pim_n.M(),dE);
      dE_IMnpip_woK0->Fill(LVec_pip_n.M(),dE);
    }
    if(K0rejectFlag && NBetaOK && MissNFlag){
      dE_IMnpim_woK0_n->Fill(LVec_pim_n.M(),dE);
      dE_IMnpip_woK0_n->Fill(LVec_pip_n.M(),dE);
    }
    if(K0rejectFlag && NBetaOK && NdEOK){
      MMom_MMass_woK0->Fill(LVec_nmiss.M(),LVec_nmiss.P());
      IMnpim_IMnpip_dE_woK0->Fill(LVec_pip_n.M(),LVec_pim_n.M());
      MMnmiss_IMnpip_dE_woK0->Fill(LVec_pip_n.M(),nmiss_mass);
      MMnmiss_IMnpim_dE_woK0->Fill(LVec_pim_n.M(),nmiss_mass);
      if(SigmaPFlag || SigmaMFlag){
        MMom_MMass_woK0_wSid->Fill(LVec_nmiss.M(),LVec_nmiss.P());
        MMnmiss_IMnpipi_woK0_wSid->Fill(LVec_pip_pim_n.M(), nmiss_mass);
      }
      if(SigmaPcutFlag[sigmacuttype]){
        MMnmiss_IMnpipi_woK0_wSid_Sp->Fill(LVec_pip_pim_n.M(), nmiss_mass);
      }
      if(SigmaMcutFlag[sigmacuttype]){
        MMnmiss_IMnpipi_woK0_wSid_Sm->Fill(LVec_pip_pim_n.M(), nmiss_mass);
      }
    }
    if(K0rejectFlag && NBetaOK && NdEOK && MissNFlag){
      IMnpim_IMnpip_dE_woK0_n->Fill(LVec_pip_n.M(),LVec_pim_n.M());
      IMnpip_CDHphi_dE_woK0_n->Fill((*CDH_Pos).Phi(),LVec_pip_n.M());
      IMnpip_CDHz_dE_woK0_n->Fill((*CDH_Pos).z(),LVec_pip_n.M());
      IMnpip_DCApip_dE_woK0_n->Fill(dca_pip_beam,LVec_pip_n.M());
      IMnpip_DCApim_dE_woK0_n->Fill(dca_pim_beam,LVec_pip_n.M());
      IMnpim_DCApip_dE_woK0_n->Fill(dca_pip_beam,LVec_pim_n.M());
      IMnpim_DCApim_dE_woK0_n->Fill(dca_pim_beam,LVec_pim_n.M());
      IMnpim_CDHphi_dE_woK0_n->Fill((*CDH_Pos).Phi(),LVec_pim_n.M());
      IMnpim_CDHz_dE_woK0_n->Fill((*CDH_Pos).z(),LVec_pim_n.M());
      IMnpip_IMnpipi_woK0_n->Fill(LVec_pip_pim_n.M(),LVec_pip_n.M());
      IMnpim_IMnpipi_woK0_n->Fill(LVec_pip_pim_n.M(),LVec_pim_n.M());
      
      for(int igap=0;igap<ngap;igap++){
        if(SigmaPsideFlag[sidebandtype][igap] || SigmaMsideFlag[sidebandtype][igap]){
          IMnpim_IMnpip_dE_woK0_n_side[igap]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
        }
        if(SigmaPsideLowFlag[igap] && SigmaPsideFlag[sidebandtype][igap]){
          IMnpim_IMnpip_dE_woK0_n_Sp_side[LOWside][igap]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
        }
        if(SigmaPsideHighFlag[igap] && SigmaPsideFlag[sidebandtype][igap]){
          IMnpim_IMnpip_dE_woK0_n_Sp_side[HIGHside][igap]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
        }
        if(SigmaMsideLowFlag[igap] && SigmaMsideFlag[sidebandtype][igap]){
          IMnpim_IMnpip_dE_woK0_n_Sm_side[LOWside][igap]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
        }
        if(SigmaMsideHighFlag[igap] && SigmaMsideFlag[sidebandtype][igap]){
          IMnpim_IMnpip_dE_woK0_n_Sm_side[HIGHside][igap]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
        }
      }
      
      for(int izone=0;izone<nzone;izone++){
        if(sidebandtype==0){
          if( !SigmaCrossPsideWideFlagLeft[izone] && !SigmaCrossPsideWideFlagRight[izone]  && SigmaPsideWideFlag[izone]){
            IMnpim_IMnpip_dE_woK0_n_Sp_sidewide[izone]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
            q_IMnpipi_woK0_wSid_n_Sp_sidewide[izone]->Fill(LVec_pip_pim_n.M(),qkn.P());
          }
          if( !SigmaCrossMsideWideFlagTop[izone] && !SigmaCrossMsideWideFlagBottom[izone]  && SigmaMsideWideFlag[izone]){
            IMnpim_IMnpip_dE_woK0_n_Sm_sidewide[izone]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
            q_IMnpipi_woK0_wSid_n_Sm_sidewide[izone]->Fill(LVec_pip_pim_n.M(),qkn.P());
          }
        }else if(sidebandtype==1){
          if( !SigmaMFlag  && SigmaPsideWideFlag[izone]){
            IMnpim_IMnpip_dE_woK0_n_Sp_sidewide[izone]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
            q_IMnpipi_woK0_wSid_n_Sp_sidewide[izone]->Fill(LVec_pip_pim_n.M(),qkn.P());
          }
          if( !SigmaPFlag  && SigmaMsideWideFlag[izone]){
            IMnpim_IMnpip_dE_woK0_n_Sm_sidewide[izone]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
            q_IMnpipi_woK0_wSid_n_Sm_sidewide[izone]->Fill(LVec_pip_pim_n.M(),qkn.P());
          }
        }else if(sidebandtype==2){
          if( !SigmawideMFlag  && SigmaPsideWideFlag[izone]){
            IMnpim_IMnpip_dE_woK0_n_Sp_sidewide[izone]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
            q_IMnpipi_woK0_wSid_n_Sp_sidewide[izone]->Fill(LVec_pip_pim_n.M(),qkn.P());
          }
          if( !SigmawidePFlag  && SigmaMsideWideFlag[izone]){
            IMnpim_IMnpip_dE_woK0_n_Sm_sidewide[izone]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
            q_IMnpipi_woK0_wSid_n_Sm_sidewide[izone]->Fill(LVec_pip_pim_n.M(),qkn.P());
          }
        }
      }
      
      MMnpip_MMnpim_woK0_n->Fill(LVec_pim_n_miss.M(),LVec_pip_n_miss.M());
      nmom_IMnpim_dE_woK0_n->Fill(LVec_pim_n.M(),(*LVec_n).P());
      nmom_IMnpip_dE_woK0_n->Fill(LVec_pip_n.M(),(*LVec_n).P());
      if(SigmaPFlag || SigmaMFlag){
        MMnpip_MMnpim_woK0_wSid_n->Fill(LVec_pim_n_miss.M(),LVec_pip_n_miss.M());
        dE_IMnpipi_woK0_wSid_n->Fill(LVec_pip_pim_n.M(),dE);
        Cosn_IMnpipi_woK0_wSid_n->Fill(LVec_pip_pim_n.M(),cos_nmiss);
        q_IMpiSigma_woK0_wSid_n_genacc->Fill(LVec_piSigma_react.M()/1000.,qkn_react.P()/1000.);
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
        nmom_IMnpipi_woK0_wSid_n_Sp->Fill(LVec_pip_pim_n.M(),(*LVec_n).P());
        q_IMpiSigma_woK0_wSid_n_Sp_genacc->Fill(LVec_piSigma_react.M()/1000.,qkn_react.P()/1000.);
        q_IMnpipi_woK0_wSid_n_Sp_acc->Fill(LVec_pip_pim_n_mc.M(),qkn_mc.P());
        q_IMnpipi_woK0_wSid_n_Sp_acc_reco->Fill(LVec_pip_pim_n.M(),qkn.P());
        if(SimSpmode){
          q_IMnpipi_woK0_wSid_n_Sp_mc->Fill(LVec_pip_pim_n_mc.M(),qkn_mc.P());
          diff_IMnpipi_woK0_wSid_n_Sp->Fill(LVec_pip_pim_n.M(),LVec_pip_pim_n.M()-LVec_pip_pim_n_mc.M());
          diff_q_woK0_wSid_n_Sp->Fill(qkn.P(),qkn.P()-qkn_mc.P());
        }
      }
      if(SigmaPcutFlag[sigmacuttype] || (!SigmaPFlag)){
        IMnpim_IMnpip_dE_woK0_n_Sp_bg->Fill(LVec_pip_n.M(),LVec_pim_n.M());
      }


      for(int itype=0;itype<3;itype++){
        for(int igap=0;igap<ngap;igap++){
          if(SigmaPsideFlag[itype][igap] && SigmaPsideLowFlag[igap] ){
            q_IMnpipi_woK0_wSid_n_Sp_side[itype][LOWside][igap]->Fill(LVec_pip_pim_n.M(),qkn.P());
          }
          if(SigmaPsideFlag[itype][igap] && SigmaPsideHighFlag[igap] ){
            q_IMnpipi_woK0_wSid_n_Sp_side[itype][HIGHside][igap]->Fill(LVec_pip_pim_n.M(),qkn.P());
          }
        }
      }

      if(SigmaMcutFlag[sigmacuttype]){
        IMnpim_IMnpip_dE_woK0_n_Sm->Fill(LVec_pip_n.M(),LVec_pim_n.M());
        q_IMnpipi_woK0_wSid_n_Sm->Fill(LVec_pip_pim_n.M(),qkn.P());
        nmom_IMnpipi_woK0_wSid_n_Sm->Fill(LVec_pip_pim_n.M(),(*LVec_n).P());
        q_IMpiSigma_woK0_wSid_n_Sm_genacc->Fill(LVec_piSigma_react.M()/1000.,qkn_react.P()/1000.);
        q_IMnpipi_woK0_wSid_n_Sm_acc->Fill(LVec_pip_pim_n_mc.M(),qkn_mc.P());
        q_IMnpipi_woK0_wSid_n_Sm_acc_reco->Fill(LVec_pip_pim_n.M(),qkn.P());
        if(SimSmmode){
          q_IMnpipi_woK0_wSid_n_Sm_mc->Fill(LVec_pip_pim_n_mc.M(),qkn_mc.P());
          diff_IMnpipi_woK0_wSid_n_Sm->Fill(LVec_pip_pim_n.M(),LVec_pip_pim_n.M()-LVec_pip_pim_n_mc.M());
          diff_q_woK0_wSid_n_Sm->Fill(qkn.P(),qkn.P()-qkn_mc.P());
        }
      }
      if(SigmaMcutFlag[sigmacuttype] || (!SigmaMFlag)){
        IMnpim_IMnpip_dE_woK0_n_Sm_bg->Fill(LVec_pip_n.M(),LVec_pim_n.M());
      }
      for(int igap=0;igap<ngap;igap++){
        for(int itype=0;itype<3;itype++){
          if(SigmaMsideFlag[itype][igap] && SigmaMsideLowFlag[igap]){
            q_IMnpipi_woK0_wSid_n_Sm_side[itype][LOWside][igap]->Fill(LVec_pip_pim_n.M(),qkn.P());
          }
          if(SigmaMsideFlag[itype][igap] && SigmaMsideHighFlag[igap]){
            q_IMnpipi_woK0_wSid_n_Sm_side[itype][HIGHside][igap]->Fill(LVec_pip_pim_n.M(),qkn.P());
          }
        }
        if(SigmaPsideFlag[sidebandtype][igap] || SigmaMsideFlag[sidebandtype][igap]){
          q_IMnpipi_woK0_wSid_n_side[igap]->Fill(LVec_pip_pim_n.M(),qkn.P());
        }
      }//for igap
    }//if K0rejectFlag && NBetaOK && NdEOK && MissNFlag0
	  //---removing K0 END---------------------------------------------- 
  
  }//for ievt
  //--- Filling Histogram END -------------------------------------------------- 


  //----------------------------------------------------------------------------
  //---Drawing Part
  //----------------------------------------------------------------------------

  TCanvas *cMMnmiss_IMnpip_dE = new TCanvas("cMMnmiss_IMnpip_dE","MMnmiss_IMnpip_dE");
  cMMnmiss_IMnpip_dE->cd();
  MMnmiss_IMnpip_dE->GetXaxis()->SetRangeUser(1.05,1.5);
  MMnmiss_IMnpip_dE->GetYaxis()->SetRangeUser(0.6,1.5);
  MMnmiss_IMnpip_dE->Draw("colz");
  
  TCanvas *cMMnmiss_IMnpip_dE_woK0 = new TCanvas("cMMnmiss_IMnpip_dE_woK0","MMnmiss_IMnpip_dE_woK0");
  cMMnmiss_IMnpip_dE_woK0->cd(); 
  MMnmiss_IMnpip_dE_woK0->GetXaxis()->SetRangeUser(1.05,1.5);
  MMnmiss_IMnpip_dE_woK0->GetYaxis()->SetRangeUser(0.6,1.5);
  MMnmiss_IMnpip_dE_woK0->Draw("colz");

  TCanvas *cMMnmiss_IMnpim_dE = new TCanvas("cMMnmiss_IMnpim_dE","MMnmiss_IMnpim_dE");
  cMMnmiss_IMnpim_dE->cd();
  MMnmiss_IMnpim_dE->GetXaxis()->SetRangeUser(1.05,1.5);
  MMnmiss_IMnpim_dE->GetYaxis()->SetRangeUser(0.6,1.5);
  MMnmiss_IMnpim_dE->Draw("colz");

  TCanvas *cMMnmiss_IMpippim_dE = new TCanvas("cMMnmiss_IMpippim_dE","MMnmiss_IMpippim_dE");
  cMMnmiss_IMpippim_dE->cd();
  MMnmiss_IMpippim_dE->Draw("colz");
  
  TCanvas *cMMnmiss_IMpippim_dE_px = new TCanvas("cMMnmiss_IMpippim_dE_px","MMnmiss_IMpippim_dE_px");
  cMMnmiss_IMpippim_dE_px->cd();
  TH1D* MMnmiss_IMpippim_dE_px = (TH1D*)MMnmiss_IMpippim_dE->ProjectionX();
  MMnmiss_IMpippim_dE_px->Draw("HE");
  TH1D* MMnmiss_IMpippim_dE_px_cut = (TH1D*)MMnmiss_IMpippim_dE_px->Clone("MMnmiss_IMpippim_dE_px_cut");
  MMnmiss_IMpippim_dE_px_cut->GetXaxis()->SetRangeUser(anacuts::pipi_MIN,anacuts::pipi_MAX);
  MMnmiss_IMpippim_dE_px_cut->SetFillColor(kViolet);
  MMnmiss_IMpippim_dE_px_cut->Draw("HEsame");
  
  TCanvas *cMMnmiss_IMpippim_dE_py = new TCanvas("cMMnmiss_IMpippim_dE_py","MMnmiss_IMpippim_dE_py");
  cMMnmiss_IMpippim_dE_py->cd();
  TH1D* MMnmiss_IMpippim_dE_py = (TH1D*)MMnmiss_IMpippim_dE->ProjectionY("MMnmiss_IMpippim_dE_py");
  //MMnmiss_IMpippim_dE_py->Draw("HE");
  
  TH1D* MMnmiss_IMpippim_dE_py_cut = (TH1D*)MMnmiss_IMpippim_dE->ProjectionY("MMnmiss_IMpippim_dE_py_cut",
    MMnmiss_IMpippim_dE->GetXaxis()->FindBin(anacuts::pipi_MIN),
    MMnmiss_IMpippim_dE->GetXaxis()->FindBin(anacuts::pipi_MAX)
  );
  MMnmiss_IMpippim_dE_py_cut->SetFillColor(kViolet);
  MMnmiss_IMpippim_dE_py_cut->Draw("HE");
  
  TCanvas *cMMnmiss_IMnpim_dE_woK0 = new TCanvas("cMMnmiss_IMnpim_dE_woK0","MMnmiss_IMnpim_dE_woK0");
  cMMnmiss_IMnpim_dE_woK0->cd();
  MMnmiss_IMnpim_dE_woK0->GetXaxis()->SetRangeUser(1.05,1.5);
  MMnmiss_IMnpim_dE_woK0->GetYaxis()->SetRangeUser(0.6,1.5);
  MMnmiss_IMnpim_dE_woK0->Draw("colz");
  
  //Sigma reconstruction 
  TCanvas *cIMnpim_IMnpip_dE_n = new TCanvas("cIMnpim_IMnpip_dE_n","IMnpim_IMnpip_dE_n");
  cIMnpim_IMnpip_dE_n->cd();
  IMnpim_IMnpip_dE_n->GetXaxis()->SetRangeUser(0,1.7);
  IMnpim_IMnpip_dE_n->GetYaxis()->SetRangeUser(0,1.7);
  IMnpim_IMnpip_dE_n->Draw("colz");
  
  TCanvas *cIMnpim_IMnpip_dE_woK0_n = new TCanvas("cIMnpim_IMnpip_dE_woK0_n","IMnpim_IMnpip_dE_woK0_n");
  cIMnpim_IMnpip_dE_woK0_n->cd();
  IMnpim_IMnpip_dE_woK0_n->GetXaxis()->SetRangeUser(0,1.7);
  IMnpim_IMnpip_dE_woK0_n->GetYaxis()->SetRangeUser(0,1.7);
  IMnpim_IMnpip_dE_woK0_n->Draw("colz");

  TCanvas *cIMnpim_IMnpip_dE_woK0_n_px = new TCanvas("cIMnpim_IMnpip_dE_woK0_n_px","IMnpim_IMnpip_dE_woK0_n_px");
  cIMnpim_IMnpip_dE_woK0_n_px->cd();
  TH1D* IMnpim_IMnpip_dE_woK0_n_px = (TH1D*)IMnpim_IMnpip_dE_woK0_n->ProjectionX("IMnpim_IMnpip_dE_woK0_n_px");
  IMnpim_IMnpip_dE_woK0_n_px->Draw("HE");
  TH1D* IMnpim_IMnpip_dE_woK0_n_px_cut = (TH1D*)IMnpim_IMnpip_dE_woK0_n_px->Clone("IMnpim_IMnpip_dE_woK0_n_px_cut");
  IMnpim_IMnpip_dE_woK0_n_px_cut->SetFillColor(2);
  IMnpim_IMnpip_dE_woK0_n_px_cut->GetXaxis()->SetRangeUser(anacuts::Sigmap_MIN,anacuts::Sigmap_MAX);
  IMnpim_IMnpip_dE_woK0_n_px_cut->Draw("HEsame");

  TCanvas *cIMnpim_IMnpip_dE_woK0_n_py = new TCanvas("cIMnpim_IMnpip_dE_woK0_n_py","IMnpim_IMnpip_dE_woK0_n_py");
  cIMnpim_IMnpip_dE_woK0_n_py->cd();
  TH1D* IMnpim_IMnpip_dE_woK0_n_py = (TH1D*)IMnpim_IMnpip_dE_woK0_n->ProjectionY("IMnpim_IMnpip_dE_woK0_n_py");
  IMnpim_IMnpip_dE_woK0_n_py->Draw("HE");
  TH1D* IMnpim_IMnpip_dE_woK0_n_py_cut = (TH1D*)IMnpim_IMnpip_dE_woK0_n_py->Clone("IMnpim_IMnpip_dE_woK0_n_py_cut");
  IMnpim_IMnpip_dE_woK0_n_py_cut->SetFillColor(3);
  IMnpim_IMnpip_dE_woK0_n_py_cut->GetXaxis()->SetRangeUser(anacuts::Sigmam_MIN,anacuts::Sigmam_MAX);
  IMnpim_IMnpip_dE_woK0_n_py_cut->Draw("HEsame");
   
  TCanvas *cIMnpip_CDHphi_dE_woK0_n = new TCanvas("cIMnpip_CDHphi_dE_woK0_n","IMnpip_CDHphi_dE_woK0_n");
  cIMnpip_CDHphi_dE_woK0_n->cd();
  IMnpip_CDHphi_dE_woK0_n->Draw("colz");
  
  TCanvas *cIMnpip_CDHz_dE_woK0_n = new TCanvas("cIMnpip_CDHz_dE_woK0_n","IMnpip_CDHz_dE_woK0_n");
  cIMnpip_CDHz_dE_woK0_n->cd();
  IMnpip_CDHz_dE_woK0_n->RebinX(5);
  IMnpip_CDHz_dE_woK0_n->Draw("colz");
  
  TCanvas *cIMnpim_CDHphi_dE_woK0_n = new TCanvas("cIMnpim_CDHphi_dE_woK0_n","IMnpim_CDHphi_dE_woK0_n");
  cIMnpim_CDHphi_dE_woK0_n->cd();
  IMnpim_CDHphi_dE_woK0_n->Draw("colz");
  
  TCanvas *cIMnpim_CDHphi_dE_woK0_n_px = new TCanvas("cIMnpim_CDHphi_dE_woK0_n_px","IMnpim_CDHphi_dE_woK0_n_px");
  cIMnpim_CDHphi_dE_woK0_n_px->cd();
  int bin105 = IMnpim_CDHphi_dE_woK0_n->GetYaxis()->FindBin(1.05);
  int bin115 = IMnpim_CDHphi_dE_woK0_n->GetYaxis()->FindBin(1.15);
  TH1D *IMnpim_CDHphi_dE_woK0_n_px = (TH1D*) IMnpim_CDHphi_dE_woK0_n->ProjectionX("IMnpim_CDHphi_dE_woK0_n_px",bin105,bin115);
  IMnpim_CDHphi_dE_woK0_n_px->Draw("HE");

  TCanvas *cIMnpim_CDHz_dE_woK0_n = new TCanvas("cIMnpim_CDHz_dE_woK0_n","IMnpim_CDHz_dE_woK0_n");
  cIMnpim_CDHz_dE_woK0_n->cd();
  IMnpim_CDHz_dE_woK0_n->RebinX(5);
  IMnpim_CDHz_dE_woK0_n->Draw("colz");
  
  TCanvas *cIMnpim_CDHz_dE_woK0_n_px = new TCanvas("cIMnpim_CDHz_dE_woK0_n_px","IMnpim_CDHz_dE_woK0_n_px");
  cIMnpim_CDHz_dE_woK0_n_px->cd();
  TH1D *IMnpim_CDHz_dE_woK0_n_px = (TH1D*)IMnpim_CDHz_dE_woK0_n->ProjectionX("IMnpim_CDHz_dE_woK0_n_px",bin105,bin115);
  IMnpim_CDHz_dE_woK0_n_px->Draw("HE");

  
  TCanvas *cIMnpip_DCApip_dE_woK0_n = new TCanvas("cIMnpip_DCApip_dE_woK0_n","IMnpip_DCApip_dE_woK0_n");
  cIMnpip_DCApip_dE_woK0_n->cd();
  cIMnpip_DCApip_dE_woK0_n->SetLogx();
  IMnpip_DCApip_dE_woK0_n->Draw("colz");

  TCanvas *cIMnpip_DCApim_dE_woK0_n = new TCanvas("cIMnpip_DCApim_dE_woK0_n","IMnpip_DCApim_dE_woK0_n");
  cIMnpip_DCApim_dE_woK0_n->cd();
  cIMnpip_DCApim_dE_woK0_n->SetLogx();
  IMnpip_DCApim_dE_woK0_n->Draw("colz");

  TCanvas *cIMnpim_DCApip_dE_woK0_n = new TCanvas("cIMnpim_DCApip_dE_woK0_n","IMnpim_DCApip_dE_woK0_n");
  cIMnpim_DCApip_dE_woK0_n->cd();
  cIMnpim_DCApip_dE_woK0_n->SetLogx();
  IMnpim_DCApip_dE_woK0_n->Draw("colz");

  TCanvas *cIMnpim_DCApim_dE_woK0_n = new TCanvas("cIMnpim_DCApim_dE_woK0_n","IMnpim_DCApim_dE_woK0_n");
  cIMnpim_DCApim_dE_woK0_n->cd();
  cIMnpim_DCApim_dE_woK0_n->SetLogx();
  IMnpim_DCApim_dE_woK0_n->Draw("colz");

  //Sigma+ or Sigma- selection
  TCanvas *cIMnpim_IMnpip_dE_n_SpSm = new TCanvas("cIMnpim_IMnpip_dE_n_SpSm","IMnpim_IMnpip_dE_n_SpSm");
  cIMnpim_IMnpip_dE_n_SpSm->cd();
  IMnpim_IMnpip_dE_n_Sp->GetXaxis()->SetRangeUser(0,1.7);
  IMnpim_IMnpip_dE_n_Sp->GetYaxis()->SetRangeUser(0,1.7);
  IMnpim_IMnpip_dE_n_Sp->SetMaximum(IMnpim_IMnpip_dE_woK0_n->GetMaximum());
  IMnpim_IMnpip_dE_n_Sp->Draw("colz");
  IMnpim_IMnpip_dE_n_Sm->Draw("colsame");
  
  TCanvas *cIMnpim_IMnpip_dE_woK0_n_SpSm = new TCanvas("cIMnpim_IMnpip_dE_woK0_n_SpSm","IMnpim_IMnpip_dE_woK0_n_SpSm");
  cIMnpim_IMnpip_dE_woK0_n_SpSm->cd();
  IMnpim_IMnpip_dE_woK0_n_Sp->GetXaxis()->SetRangeUser(0,1.7);
  IMnpim_IMnpip_dE_woK0_n_Sp->GetYaxis()->SetRangeUser(0,1.7);
  IMnpim_IMnpip_dE_woK0_n_Sp->SetMaximum(IMnpim_IMnpip_dE_woK0_n->GetMaximum());
  IMnpim_IMnpip_dE_woK0_n_Sp->Draw("colz");
  IMnpim_IMnpip_dE_woK0_n_Sm->Draw("colsame");
  


  
  TCanvas *cIMnpim_IMnpip_dE_woK0_n_Sp[ngap];
  for(int igap=0;igap<ngap;igap++){
    cIMnpim_IMnpip_dE_woK0_n_Sp[igap] = new TCanvas(Form("cIMnpim_IMnpip_dE_woK0_n_Sp_%d",igap),Form("IMnpim_IMnpip_dE_woK0_n_Sp_%d",igap));   
    cIMnpim_IMnpip_dE_woK0_n_Sp[igap]->cd();
    IMnpim_IMnpip_dE_woK0_n_Sp->SetMaximum(IMnpim_IMnpip_dE_woK0_n->GetMaximum());
    IMnpim_IMnpip_dE_woK0_n_Sp->Draw("colz");
    IMnpim_IMnpip_dE_woK0_n_Sp_side[LOWside][igap]->Draw("colsame");
    IMnpim_IMnpip_dE_woK0_n_Sp_side[HIGHside][igap]->Draw("colsame");
    if(SimSpmode || SimSmmode) break;
  }

  TCanvas *cIMnpim_IMnpip_dE_woK0_n_Sm[ngap];
  for(int igap=0;igap<ngap;igap++){
    cIMnpim_IMnpip_dE_woK0_n_Sm[igap] = new TCanvas(Form("cIMnpim_IMnpip_dE_woK0_n_Sm_%d",igap),Form("IMnpim_IMnpip_dE_woK0_n_Sm_%d",igap));
    cIMnpim_IMnpip_dE_woK0_n_Sm[igap]->cd();
    IMnpim_IMnpip_dE_woK0_n_Sm->SetMaximum(IMnpim_IMnpip_dE_woK0_n->GetMaximum());
    IMnpim_IMnpip_dE_woK0_n_Sm->GetXaxis()->SetRangeUser(0,1.7);
    IMnpim_IMnpip_dE_woK0_n_Sm->GetYaxis()->SetRangeUser(0,1.7);
    IMnpim_IMnpip_dE_woK0_n_Sm->Draw("colz");
    IMnpim_IMnpip_dE_woK0_n_Sm_side[LOWside][igap]->Draw("colsame");
    IMnpim_IMnpip_dE_woK0_n_Sm_side[HIGHside][igap]->Draw("colsame");
    if(SimSpmode || SimSmmode) break;
  }

  TCanvas *cIMnpim_IMnpip_dE_n_Sp_wide[nzone];
  TH1D* q_IMnpipi_wSid_n_Sp_sidewide_px[nzone];
  for(int izone=0;izone<nzone;izone++){
    cIMnpim_IMnpip_dE_n_Sp_wide[izone] = new TCanvas(Form("cIMnpim_IMnpip_dE_n_Sp_wide_%d",izone),Form("IMnpim_IMnpip_dE_n_Sp_wide_%d",izone));   
    cIMnpim_IMnpip_dE_n_Sp_wide[izone]->Divide(2,2);
    cIMnpim_IMnpip_dE_n_Sp_wide[izone]->cd(1);
    IMnpim_IMnpip_dE_n_Sp->SetMaximum(IMnpim_IMnpip_dE_woK0_n->GetMaximum());
    IMnpim_IMnpip_dE_n_Sp->GetXaxis()->SetRangeUser(0,1.7);
    IMnpim_IMnpip_dE_n_Sp->GetYaxis()->SetRangeUser(0,1.7);
    IMnpim_IMnpip_dE_n_Sp->Draw("colz");
    IMnpim_IMnpip_dE_n_Sp_sidewide[izone]->Draw("colsame");
    cIMnpim_IMnpip_dE_n_Sp_wide[izone]->cd(2);
    TH1D* q_IMnpipi_wSid_n_Sp_px = (TH1D*) q_IMnpipi_wSid_n_Sp->ProjectionX();
    q_IMnpipi_wSid_n_Sp_px->Draw("HE");
    q_IMnpipi_wSid_n_Sp_sidewide_px[izone] = (TH1D*)q_IMnpipi_wSid_n_Sp_sidewide[izone]->ProjectionX(Form("Sp_sidewide_px_%d",izone));
    q_IMnpipi_wSid_n_Sp_sidewide_px[izone]->SetLineColor(kCyan);
    q_IMnpipi_wSid_n_Sp_sidewide_px[izone]->SetMarkerColor(kCyan);
    q_IMnpipi_wSid_n_Sp_sidewide_px[izone]->SetMarkerStyle(22);
    q_IMnpipi_wSid_n_Sp_sidewide_px[izone]->Draw("HEsame");
    cIMnpim_IMnpip_dE_n_Sp_wide[izone]->cd(3);
    q_IMnpipi_wSid_n_Sp->Draw("colz");
    cIMnpim_IMnpip_dE_n_Sp_wide[izone]->cd(4);
    q_IMnpipi_wSid_n_Sp_sidewide[izone]->SetMaximum(q_IMnpipi_wSid_n_Sp->GetMaximum());
    q_IMnpipi_wSid_n_Sp_sidewide[izone]->Draw("colz"); 
    
    if(SimSpmode || SimSmmode) break;
  }

  TCanvas *cIMnpim_IMnpip_dE_woK0_n_Sp_wide[nzone];
  TH1D* q_IMnpipi_woK0_wSid_n_Sp_sidewide_px[nzone];
  for(int izone=0;izone<nzone;izone++){
    cIMnpim_IMnpip_dE_woK0_n_Sp_wide[izone] = new TCanvas(Form("cIMnpim_IMnpip_dE_woK0_n_Sp_wide_%d",izone),Form("IMnpim_IMnpip_dE_woK0_n_Sp_wide_%d",izone));   
    cIMnpim_IMnpip_dE_woK0_n_Sp_wide[izone]->Divide(2,2);
    cIMnpim_IMnpip_dE_woK0_n_Sp_wide[izone]->cd(1);
    IMnpim_IMnpip_dE_woK0_n_Sp->SetMaximum(IMnpim_IMnpip_dE_woK0_n->GetMaximum());
    IMnpim_IMnpip_dE_woK0_n_Sp->GetXaxis()->SetRangeUser(0,1.7);
    IMnpim_IMnpip_dE_woK0_n_Sp->GetYaxis()->SetRangeUser(0,1.7);
    IMnpim_IMnpip_dE_woK0_n_Sp->Draw("colz");
    IMnpim_IMnpip_dE_woK0_n_Sp_sidewide[izone]->Draw("colsame");
    cIMnpim_IMnpip_dE_woK0_n_Sp_wide[izone]->cd(2);
    TH1D* q_IMnpipi_woK0_wSid_n_Sp_px = (TH1D*) q_IMnpipi_woK0_wSid_n_Sp->ProjectionX();
    q_IMnpipi_woK0_wSid_n_Sp_px->Draw("HE");
    q_IMnpipi_woK0_wSid_n_Sp_sidewide_px[izone] = (TH1D*)q_IMnpipi_woK0_wSid_n_Sp_sidewide[izone]->ProjectionX(Form("Sp_sidewide_woK0_px_%d",izone));
    q_IMnpipi_woK0_wSid_n_Sp_sidewide_px[izone]->SetLineColor(kCyan);
    q_IMnpipi_woK0_wSid_n_Sp_sidewide_px[izone]->SetMarkerColor(kCyan);
    q_IMnpipi_woK0_wSid_n_Sp_sidewide_px[izone]->SetMarkerStyle(22);
    q_IMnpipi_woK0_wSid_n_Sp_sidewide_px[izone]->Draw("HEsame");
    cIMnpim_IMnpip_dE_woK0_n_Sp_wide[izone]->cd(3);
    q_IMnpipi_woK0_wSid_n_Sp->Draw("colz");
    cIMnpim_IMnpip_dE_woK0_n_Sp_wide[izone]->cd(4);
    q_IMnpipi_woK0_wSid_n_Sp_sidewide[izone]->SetMaximum(q_IMnpipi_woK0_wSid_n_Sp->GetMaximum());
    q_IMnpipi_woK0_wSid_n_Sp_sidewide[izone]->Draw("colz"); 
    
    if(SimSpmode || SimSmmode) break;
  }
  
  TCanvas *cIMnpim_IMnpip_dE_n_Sm_wide[nzone];
  TH1D* q_IMnpipi_wSid_n_Sm_sidewide_px[nzone];
  for(int izone=0;izone<nzone;izone++){
    cIMnpim_IMnpip_dE_n_Sm_wide[izone] = new TCanvas(Form("cIMnpim_IMnpip_dE_n_Sm_wide_%d",izone),Form("IMnpim_IMnpip_dE_n_Sm_wide_%d",izone));
    cIMnpim_IMnpip_dE_n_Sm_wide[izone]->Divide(2,2);
    cIMnpim_IMnpip_dE_n_Sm_wide[izone]->cd(1);
    IMnpim_IMnpip_dE_n_Sm->SetMaximum(IMnpim_IMnpip_dE_n->GetMaximum());
    IMnpim_IMnpip_dE_n_Sm->GetXaxis()->SetRangeUser(0,1.7);
    IMnpim_IMnpip_dE_n_Sm->GetYaxis()->SetRangeUser(0,1.7);
    IMnpim_IMnpip_dE_n_Sm->Draw("colz");
    IMnpim_IMnpip_dE_n_Sm_sidewide[izone]->Draw("colsame");
    cIMnpim_IMnpip_dE_n_Sm_wide[izone]->cd(2);
    TH1D* q_IMnpipi_wSid_n_Sm_px = (TH1D*) q_IMnpipi_wSid_n_Sm->ProjectionX();
    q_IMnpipi_wSid_n_Sm_px->Draw("HE");
    q_IMnpipi_wSid_n_Sm_sidewide_px[izone] = (TH1D*)q_IMnpipi_wSid_n_Sm_sidewide[izone]->ProjectionX(Form("Sm_sidewide_px_%d",izone));
    q_IMnpipi_wSid_n_Sm_sidewide_px[izone]->SetLineColor(kPink);
    q_IMnpipi_wSid_n_Sm_sidewide_px[izone]->SetMarkerColor(kPink);
    q_IMnpipi_wSid_n_Sm_sidewide_px[izone]->SetMarkerStyle(22);
    q_IMnpipi_wSid_n_Sm_sidewide_px[izone]->Draw("HEsame");
    cIMnpim_IMnpip_dE_n_Sm_wide[izone]->cd(3);
    q_IMnpipi_wSid_n_Sm->Draw("colz");
    cIMnpim_IMnpip_dE_n_Sm_wide[izone]->cd(4);
    q_IMnpipi_wSid_n_Sm_sidewide[izone]->SetMaximum(q_IMnpipi_woK0_wSid_n_Sm->GetMaximum());
    q_IMnpipi_wSid_n_Sm_sidewide[izone]->Draw("colz"); 

    if(SimSpmode || SimSmmode) break;
  }

  TCanvas *cIMnpim_IMnpip_dE_woK0_n_Sm_wide[nzone];
  TH1D* q_IMnpipi_woK0_wSid_n_Sm_sidewide_px[nzone];
  for(int izone=0;izone<nzone;izone++){
    cIMnpim_IMnpip_dE_woK0_n_Sm_wide[izone] = new TCanvas(Form("cIMnpim_IMnpip_dE_woK0_n_Sm_wide_%d",izone),Form("IMnpim_IMnpip_dE_woK0_n_Sm_wide_%d",izone));
    cIMnpim_IMnpip_dE_woK0_n_Sm_wide[izone]->Divide(2,2);
    cIMnpim_IMnpip_dE_woK0_n_Sm_wide[izone]->cd(1);
    IMnpim_IMnpip_dE_woK0_n_Sm->SetMaximum(IMnpim_IMnpip_dE_woK0_n->GetMaximum());
    IMnpim_IMnpip_dE_woK0_n_Sm->GetXaxis()->SetRangeUser(0,1.7);
    IMnpim_IMnpip_dE_woK0_n_Sm->GetYaxis()->SetRangeUser(0,1.7);
    IMnpim_IMnpip_dE_woK0_n_Sm->Draw("colz");
    IMnpim_IMnpip_dE_woK0_n_Sm_sidewide[izone]->Draw("colsame");
    cIMnpim_IMnpip_dE_woK0_n_Sm_wide[izone]->cd(2);
    TH1D* q_IMnpipi_woK0_wSid_n_Sm_px = (TH1D*) q_IMnpipi_woK0_wSid_n_Sm->ProjectionX();
    q_IMnpipi_woK0_wSid_n_Sm_px->Draw("HE");
    q_IMnpipi_woK0_wSid_n_Sm_sidewide_px[izone] = (TH1D*)q_IMnpipi_woK0_wSid_n_Sm_sidewide[izone]->ProjectionX(Form("Sm_sidewide_woK0_px_%d",izone));
    q_IMnpipi_woK0_wSid_n_Sm_sidewide_px[izone]->SetLineColor(kPink);
    q_IMnpipi_woK0_wSid_n_Sm_sidewide_px[izone]->SetMarkerColor(kPink);
    q_IMnpipi_woK0_wSid_n_Sm_sidewide_px[izone]->SetMarkerStyle(22);
    q_IMnpipi_woK0_wSid_n_Sm_sidewide_px[izone]->Draw("HEsame");
    cIMnpim_IMnpip_dE_woK0_n_Sm_wide[izone]->cd(3);
    q_IMnpipi_woK0_wSid_n_Sm->Draw("colz");
    cIMnpim_IMnpip_dE_woK0_n_Sm_wide[izone]->cd(4);
    q_IMnpipi_woK0_wSid_n_Sm_sidewide[izone]->SetMaximum(q_IMnpipi_woK0_wSid_n_Sm->GetMaximum());
    q_IMnpipi_woK0_wSid_n_Sm_sidewide[izone]->Draw("colz"); 

    if(SimSpmode || SimSmmode) break;
  }



  TCanvas *cIMnpim_IMnpip_dE_woK0_n_Sp_bg = new TCanvas("cIMnpim_IMnpip_dE_woK0_n_Sp_bg","IMnpim_IMnpip_dE_woK0_n_Sp_bg");
  cIMnpim_IMnpip_dE_woK0_n_Sp_bg->cd();
  IMnpim_IMnpip_dE_woK0_n_Sp_bg->SetMaximum(IMnpim_IMnpip_dE_woK0_n->GetMaximum());
  IMnpim_IMnpip_dE_woK0_n_Sp_bg->GetXaxis()->SetRangeUser(0,1.7);
  IMnpim_IMnpip_dE_woK0_n_Sp_bg->GetYaxis()->SetRangeUser(0,1.7);
  IMnpim_IMnpip_dE_woK0_n_Sp_bg->Draw("colz");

  TCanvas *cIMnpim_IMnpip_dE_woK0_n_Sm_bg = new TCanvas("cIMnpim_IMnpip_dE_woK0_n_Sm_bg","IMnpim_IMnpip_dE_woK0_n_Sm_bg");
  cIMnpim_IMnpip_dE_woK0_n_Sm_bg->cd();
  IMnpim_IMnpip_dE_woK0_n_Sm_bg->SetMaximum(IMnpim_IMnpip_dE_woK0_n->GetMaximum());
  IMnpim_IMnpip_dE_woK0_n_Sm_bg->GetXaxis()->SetRangeUser(0,1.7);
  IMnpim_IMnpip_dE_woK0_n_Sm_bg->GetYaxis()->SetRangeUser(0,1.7);
  IMnpim_IMnpip_dE_woK0_n_Sm_bg->Draw("colz");
  //
  TCanvas *cIMnpim_IMnpip_dE_woK0_n_px2 = new TCanvas("cIMnpim_IMnpip_dE_woK0_n_px2","IMnpim_IMnpip_dE_woK0_n_px2");
  cIMnpim_IMnpip_dE_woK0_n_px2->cd();
  //TH1D *IMnpim_IMnpip_dE_woK0_n_px = IMnpim_IMnpip_dE_woK0_n->ProjectionX();
  TH1D *IMnpim_IMnpip_dE_woK0_n_Sp_bg_px = IMnpim_IMnpip_dE_woK0_n_Sp_bg->ProjectionX();
  IMnpim_IMnpip_dE_woK0_n_Sp_bg_px->GetXaxis()->SetRangeUser(1.0,1.3);
  //IMnpim_IMnpip_dE_woK0_n_px->GetYaxis()->SetTitle("Counts/ 10 MeV/c^{2}");
  IMnpim_IMnpip_dE_woK0_n_Sp_bg_px->GetYaxis()->CenterTitle();
  IMnpim_IMnpip_dE_woK0_n_Sp_bg_px->Draw("EH");
  //
  //
  TH1D* IMnpim_IMnpip_dE_woK0_n_px_1 = (TH1D*) IMnpim_IMnpip_dE_woK0_n_Sp_bg_px->Clone();
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
  IMnpim_IMnpip_dE_woK0_n_Sp_bg_px->Fit("f1","","",anacuts::Sigmap_MIN,anacuts::Sigmap_MAX);
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
  //TH1D *IMnpim_IMnpip_dE_woK0_n_Sp_side_0_px = IMnpim_IMnpip_dE_woK0_n_Sp_side[0]->ProjectionX();
  //IMnpim_IMnpip_dE_woK0_n_Sp_side_0_px->SetFillColor(4);
  //IMnpim_IMnpip_dE_woK0_n_Sp_side_0_px->SetFillStyle(3009);
  //IMnpim_IMnpip_dE_woK0_n_Sp_side_0_px->Draw("HEsame");
  //TH1D *IMnpim_IMnpip_dE_woK0_n_Sp_side_1_px = IMnpim_IMnpip_dE_woK0_n_Sp_side[1]->ProjectionX();
  //IMnpim_IMnpip_dE_woK0_n_Sp_side_1_px->SetFillColor(4);
  //IMnpim_IMnpip_dE_woK0_n_Sp_side_1_px->SetFillStyle(3009);
  //IMnpim_IMnpip_dE_woK0_n_Sp_side_1_px->Draw("HEsame");
  
  std::cout << "Sigma+ signal region:    " << IMnpim_IMnpip_dE_woK0_n_px_1->Integral() << std::endl;
  std::cout << "Sigma+ sideband low:     " << IMnpim_IMnpip_dE_woK0_n_px_2->Integral() << std::endl;
  ///std::cout << "Sigma+ sideband low cut: " << IMnpim_IMnpip_dE_woK0_n_Sp_side_0_px->Integral() << std::endl;
  std::cout << "Sigma+ sideband high:    " << IMnpim_IMnpip_dE_woK0_n_px_3->Integral() << std::endl;
  //std::cout << "Sigma+ sideband high cut:" << IMnpim_IMnpip_dE_woK0_n_Sp_side_1_px->Integral() << std::endl;
  std::cout << "bg (Integral)           :" << bgsp       << std::endl;
  std::cout << "bg (Integral)/binw      :" << bgsp/IMnpim_IMnpip_dE_woK0_n_px_1->GetBinWidth(100) << std::endl;
  double trapezoidbgSp = (fbgSp->Eval(anacuts::Sigmap_MIN)+fbgSp->Eval(anacuts::Sigmap_MAX))*
                         (anacuts::Sigmap_MAX-anacuts::Sigmap_MIN)
                         /IMnpim_IMnpip_dE_woK0_n_px_1->GetBinWidth(100)/2.0;
  std::cout << "bg (trapezoid)          :" << trapezoidbgSp << std::endl;



  TCanvas *cIMnpim_IMnpip_dE_woK0_n_py2 = new TCanvas("cIMnpim_IMnpip_dE_woK0_n_py2","IMnpim_IMnpip_dE_woK0_n_py2");
  cIMnpim_IMnpip_dE_woK0_n_py2->cd();
 // TH1D *IMnpim_IMnpip_dE_woK0_n_py = IMnpim_IMnpip_dE_woK0_n->ProjectionY();
  TH1D *IMnpim_IMnpip_dE_woK0_n_Sm_bg_py = IMnpim_IMnpip_dE_woK0_n_Sm_bg->ProjectionY();
  IMnpim_IMnpip_dE_woK0_n_Sm_bg_py->GetXaxis()->SetRangeUser(1.0,1.3);
  //IMnpim_IMnpip_dE_woK0_n_py->GetYaxis()->SetTitle("Counts/ 10 MeV/c^{2}");
  IMnpim_IMnpip_dE_woK0_n_Sm_bg_py->GetYaxis()->CenterTitle();
  IMnpim_IMnpip_dE_woK0_n_Sm_bg_py->Draw("EH");
  
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
  IMnpim_IMnpip_dE_woK0_n_Sm_bg_py->Fit("fbgSm","","",anacuts::Sigmam_MIN,anacuts::Sigmam_MAX);
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
  
  //TH1D *IMnpim_IMnpip_dE_woK0_n_Sm_side_0_py = IMnpim_IMnpip_dE_woK0_n_Sm_side[0]->ProjectionY();
  //IMnpim_IMnpip_dE_woK0_n_Sm_side_0_py->SetFillColor(4);
  //IMnpim_IMnpip_dE_woK0_n_Sm_side_0_py->SetFillStyle(3009);
  //IMnpim_IMnpip_dE_woK0_n_Sm_side_0_py->Draw("HEsame");
  //TH1D *IMnpim_IMnpip_dE_woK0_n_Sm_side_1_py = IMnpim_IMnpip_dE_woK0_n_Sm_side[1]->ProjectionY();
  //IMnpim_IMnpip_dE_woK0_n_Sm_side_1_py->SetFillColor(4);
  //IMnpim_IMnpip_dE_woK0_n_Sm_side_1_py->SetFillStyle(3009);
  //IMnpim_IMnpip_dE_woK0_n_Sm_side_1_py->Draw("HEsame");

  std::cout << "Sigma- signal region:    " << IMnpim_IMnpip_dE_woK0_n_py_1->Integral() << std::endl;
  std::cout << "Sigma- sideband low:     " << IMnpim_IMnpip_dE_woK0_n_py_2->Integral() << std::endl;
  //std::cout << "Sigma- sideband low cut: " << IMnpim_IMnpip_dE_woK0_n_Sm_side_0_py->Integral() << std::endl;
  std::cout << "Sigma- sideband high:    " << IMnpim_IMnpip_dE_woK0_n_py_3->Integral() << std::endl;
  //std::cout << "Sigma- sideband high cut:" << IMnpim_IMnpip_dE_woK0_n_Sm_side_1_py->Integral() << std::endl;
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
  
  const double Kp_mass = pMass + kpMass;  
  TF1 *fkp = new TF1("f", "sqrt(((x*x-[0]*[0]-[1]*[1])/(2*[0]))*((x*x-[0]*[0]-[1]*[1])/(2*[0]))-[1]*[1])",Kp_mass-0.0001,2);
  fkp->SetParameter(0,nMass);
  fkp->SetParameter(1,kpMass);
  //fkp->SetLineColor(4);
  fkp->SetLineWidth(5);
  fkp->SetLineStyle(4);
  fkp->SetLineColorAlpha(kPink, 0.35);
  TCanvas *cq_IMnpipi_wSid_n_Sp = new TCanvas("cq_IMnpipi_wSid_n_Sp","q_IMnpipi_wSid_n_Sp");
  cq_IMnpipi_wSid_n_Sp->cd();
  q_IMnpipi_wSid_n_Sp->Draw("colz");
  fkp->Draw("same");
  TCanvas *cq_IMnpipi_woK0_wSid_n_Sp = new TCanvas("cq_IMnpipi_woK0_wSid_n_Sp","q_IMnpipi_woK0_wSid_n_Sp");
  cq_IMnpipi_woK0_wSid_n_Sp->cd();
  q_IMnpipi_woK0_wSid_n_Sp->Draw("colz");
  fkp->Draw("same");
  TFile *fnuSp = new TFile("NumericalRootFinder_Spmode.root");
  fnuSp->cd();
  TMultiGraph *mg = (TMultiGraph*)fnuSp->Get("mg"); 
  mg->Draw("c");
  f->cd();
  

  TCanvas *cq_IMnpipi_wK0_n = new TCanvas("cq_IMnpipi_wK0_n","q_IMnpipi_wK0_n");
  cq_IMnpipi_wK0_n->cd();
  q_IMnpipi_wK0_n->GetYaxis()->SetRangeUser(0,1.0);
  q_IMnpipi_wK0_n->Draw("colz");
  fkp->Draw("same");

  TCanvas *cq_IMnpipi_wK0_n_px = new TCanvas("cq_IMnpipi_wK0_n_px","q_IMnpipi_wK0_n_px");
  cq_IMnpipi_wK0_n_px->cd();
  TH1D* q_IMnpipi_wK0_n_px = (TH1D*)q_IMnpipi_wK0_n->ProjectionX("q_IMnpipi_wK0_n_px");
  q_IMnpipi_wK0_n_px->SetFillColor(kViolet);
  //q_IMnpipi_wK0_n_px->SetMaximum(2200);
  q_IMnpipi_wK0_n_px->Draw("HE");
  
  TCanvas *cq_IMnpipi_wK0_n_py = new TCanvas("cq_IMnpipi_wK0_n_py","q_IMnpipi_wK0_n_py");
  cq_IMnpipi_wK0_n_py->cd();
  TH1D* q_IMnpipi_wK0_n_py = (TH1D*)q_IMnpipi_wK0_n->ProjectionY("q_IMnpipi_wK0_n_py");
  q_IMnpipi_wK0_n_py->SetFillColor(kViolet);
  //q_IMnpipi_wK0_n_py->SetMaximum(4800);
  q_IMnpipi_wK0_n_py->Draw("HE");



  TCanvas *cq_IMnpipi_wK0_wSid_n_Sp = new TCanvas("cq_IMnpipi_wK0_wSid_n_Sp","q_IMnpipi_wK0_wSid_n_Sp");
  cq_IMnpipi_wK0_wSid_n_Sp->cd();
  q_IMnpipi_wK0_wSid_n_Sp->Draw("colz");
  fkp->Draw("same");
  
  TCanvas *cq_IMnpipi_wSid_n_Sm = new TCanvas("cq_IMnpipi_wSid_n_Sm","q_IMnpipi_wSid_n_Sm");
  cq_IMnpipi_wSid_n_Sm->cd();
//  q_IMnpipi_woK0_wSid_n_Sm->SetMaximum(q_IMnpipi_woK0_wSid_n_Sp->GetMaximum());
  q_IMnpipi_wSid_n_Sm->Draw("colz");
  fkp->Draw("same");



  TCanvas *cq_IMnpipi_woK0_wSid_n_Sm = new TCanvas("cq_IMnpipi_woK0_wSid_n_Sm","q_IMnpipi_woK0_wSid_n_Sm");
  cq_IMnpipi_woK0_wSid_n_Sm->cd();
//  q_IMnpipi_woK0_wSid_n_Sm->SetMaximum(q_IMnpipi_woK0_wSid_n_Sp->GetMaximum());
  q_IMnpipi_woK0_wSid_n_Sm->Draw("colz");
  fkp->Draw("same");
  TFile *fnuSm = new TFile("NumericalRootFinder_Smmode.root");
  fnuSm->cd();
  TMultiGraph *mg2 = (TMultiGraph*)fnuSm->Get("mg"); 
  mg2->Draw("c");
  f->cd();

  TCanvas *cq_IMnpipi_wK0_wSid_n_Sm = new TCanvas("cq_IMnpipi_wK0_wSid_n_Sm","q_IMnpipi_wK0_wSid_n_Sm");
  cq_IMnpipi_wK0_wSid_n_Sm->cd();
//  q_IMnpipi_woK0_wSid_n_Sm->SetMaximum(q_IMnpipi_wK0_wSid_n_Sp->GetMaximum());
  q_IMnpipi_wK0_wSid_n_Sm->Draw("colz");
  fkp->Draw("same");
  
  TCanvas *cq_IMnpipi_wSid_n_Sp_side_low = new TCanvas("cq_IMnpipi_wSid_n_Sp_side_low","q_IMnpipi_wSid_n_Sp_side_low");
  cq_IMnpipi_wSid_n_Sp_side_low->cd();
  q_IMnpipi_wSid_n_Sp_side[sidebandtype][LOWside][0]->Draw("colz");
  TCanvas *cq_IMnpipi_woK0_wSid_n_Sp_side_low = new TCanvas("cq_IMnpipi_woK0_wSid_n_Sp_side_low","q_IMnpipi_woK0_wSid_n_Sp_side_low");
  cq_IMnpipi_woK0_wSid_n_Sp_side_low->cd();
  q_IMnpipi_woK0_wSid_n_Sp_side[sidebandtype][LOWside][0]->Draw("colz");
  
  TCanvas *cq_IMnpipi_wSid_n_Sp_side_high = new TCanvas("cq_IMnpipi_wSid_n_Sp_side_high","q_IMnpipi_wSid_n_Sp_side_high");
  cq_IMnpipi_wSid_n_Sp_side_high->cd();
  q_IMnpipi_wSid_n_Sp_side[sidebandtype][HIGHside][0]->Draw("colz");
  TCanvas *cq_IMnpipi_woK0_wSid_n_Sp_side_high = new TCanvas("cq_IMnpipi_woK0_wSid_n_Sp_side_high","q_IMnpipi_woK0_wSid_n_Sp_side_high");
  cq_IMnpipi_woK0_wSid_n_Sp_side_high->cd();
  q_IMnpipi_woK0_wSid_n_Sp_side[sidebandtype][HIGHside][0]->Draw("colz");



  TCanvas *cq_IMnpipi_wSid_n_Sm_side_low = new TCanvas("cq_IMnpipi_wSid_n_Sm_side_low","q_IMnpipi_wSid_n_Sm_side_low");
  cq_IMnpipi_wSid_n_Sm_side_low->cd();
  q_IMnpipi_wSid_n_Sm_side[sidebandtype][LOWside][0]->Draw("colz");
  TCanvas *cq_IMnpipi_woK0_wSid_n_Sm_side_low = new TCanvas("cq_IMnpipi_woK0_wSid_n_Sm_side_low","q_IMnpipi_woK0_wSid_n_Sm_side_low");
  cq_IMnpipi_woK0_wSid_n_Sm_side_low->cd();
  q_IMnpipi_woK0_wSid_n_Sm_side[sidebandtype][LOWside][0]->Draw("colz");

  TCanvas *cq_IMnpipi_wSid_n_Sm_side_high = new TCanvas("cq_IMnpipi_wSid_n_Sm_side_high","q_IMnpipi_wSid_n_Sm_side_high");
  cq_IMnpipi_wSid_n_Sm_side_high->cd();
  q_IMnpipi_wSid_n_Sm_side[sidebandtype][HIGHside][0]->Draw("colz");
  TCanvas *cq_IMnpipi_woK0_wSid_n_Sm_side_high = new TCanvas("cq_IMnpipi_woK0_wSid_n_Sm_side_high","q_IMnpipi_woK0_wSid_n_Sm_side_high");
  cq_IMnpipi_woK0_wSid_n_Sm_side_high->cd();
  q_IMnpipi_woK0_wSid_n_Sm_side[sidebandtype][HIGHside][0]->Draw("colz");

  TCanvas *cq_IMnpipi_wSid_n_px_SpSm = new TCanvas("cq_IMnpipi_wSid_n_px_SpSm","q_IMnpipi_wSid_n_px_SpSm"); 
  cq_IMnpipi_wSid_n_px_SpSm->cd();
  TH1D *q_IMnpipi_wSid_n_px = q_IMnpipi_wSid_n->ProjectionX();
  q_IMnpipi_wSid_n_px->Draw("EH");
  TH1D* q_IMnpipi_wSid_n_Sp_px = (TH1D*) q_IMnpipi_wSid_n_Sp->ProjectionX();
  TH1D* q_IMnpipi_wSid_n_Sm_px = (TH1D*) q_IMnpipi_wSid_n_Sm->ProjectionX();
  q_IMnpipi_wSid_n_Sp_px->SetLineColor(2);
  q_IMnpipi_wSid_n_Sm_px->SetLineColor(3);
  //q_IMnpipi_wSid_n_px1->Draw("HEsame");
  q_IMnpipi_wSid_n_Sp_px->Draw("HEsame");
  q_IMnpipi_wSid_n_Sm_px->Draw("HEsame");
  //TH1D* q_IMnpipi_wSid_n_sum = (TH1D*)q_IMnpipi_wSid_n_Sp_px->Clone();
  //q_IMnpipi_wSid_n_sum->Add(q_IMnpipi_wSid_n_Sm_px,1);
  //q_IMnpipi_wSid_n_sum->Draw("HEsame");

  TCanvas *cq_IMnpipi_woK0_wSid_n_px_SpSm = new TCanvas("cq_IMnpipi_woK0_wSid_n_px_SpSm","q_IMnpipi_woK0_wSid_n_px_SpSm"); 
  cq_IMnpipi_woK0_wSid_n_px_SpSm->cd();
  TH1D *q_IMnpipi_woK0_wSid_n_px = q_IMnpipi_woK0_wSid_n->ProjectionX();
  q_IMnpipi_woK0_wSid_n_px->Draw("EH");
  TH1D* q_IMnpipi_woK0_wSid_n_Sp_px = (TH1D*) q_IMnpipi_woK0_wSid_n_Sp->ProjectionX();
  TH1D* q_IMnpipi_woK0_wSid_n_Sm_px = (TH1D*) q_IMnpipi_woK0_wSid_n_Sm->ProjectionX();
  q_IMnpipi_woK0_wSid_n_Sp_px->SetLineColor(2);
  q_IMnpipi_woK0_wSid_n_Sm_px->SetLineColor(3);
  //q_IMnpipi_wSid_n_px1->Draw("HEsame");
  q_IMnpipi_woK0_wSid_n_Sp_px->Draw("HEsame");
  q_IMnpipi_woK0_wSid_n_Sm_px->Draw("HEsame");
  //TH1D* q_IMnpipi_woK0_wSid_n_sum = (TH1D*)q_IMnpipi_woK0_wSid_n_Sp_px->Clone();
  // q_IMnpipi_woK0_wSid_n_sum->Add(q_IMnpipi_woK0_wSid_n_Sm_px,1);
  //q_IMnpipi_woK0_wSid_n_sum->SetLineColor(4);
  //q_IMnpipi_woK0_wSid_n_wocross->Draw("HEsame");
  
  TCanvas *cq_IMnpipi_wwoK0_wSid_n_px_SpSm = new TCanvas("cq_IMnpipi_wwoK0_wSid_n_px_SpSm","q_IMnpipi_wwoK0_wSid_n_px_SpSm"); 
  cq_IMnpipi_wwoK0_wSid_n_px_SpSm->cd();
  cq_IMnpipi_wwoK0_wSid_n_px_SpSm->Divide(2,1);
  cq_IMnpipi_wwoK0_wSid_n_px_SpSm->cd(1);
  TH1D *q_IMnpipi_wSid_n_Sp_px_clone = (TH1D*)q_IMnpipi_wSid_n_Sp_px->Clone();
  TH1D *q_IMnpipi_woK0_wSid_n_Sp_px_clone = (TH1D*)q_IMnpipi_woK0_wSid_n_Sp_px->Clone();
  q_IMnpipi_wSid_n_Sp_px_clone->SetMarkerColor(1);
  q_IMnpipi_wSid_n_Sp_px_clone->SetLineColor(1);
  q_IMnpipi_wSid_n_Sp_px_clone->SetMarkerStyle(20);
  q_IMnpipi_wSid_n_Sp_px_clone->Draw("E");
  q_IMnpipi_woK0_wSid_n_Sp_px_clone->SetMarkerStyle(21);
  q_IMnpipi_woK0_wSid_n_Sp_px_clone->SetMarkerColor(2);
  q_IMnpipi_woK0_wSid_n_Sp_px_clone->Draw("Esame");
  
  cq_IMnpipi_wwoK0_wSid_n_px_SpSm->cd(2);
  TH1D *q_IMnpipi_wSid_n_Sm_px_clone = (TH1D*)q_IMnpipi_wSid_n_Sm_px->Clone();
  TH1D *q_IMnpipi_woK0_wSid_n_Sm_px_clone = (TH1D*)q_IMnpipi_woK0_wSid_n_Sm_px->Clone();
  q_IMnpipi_wSid_n_Sm_px_clone->SetMarkerColor(1);
  q_IMnpipi_wSid_n_Sm_px_clone->SetLineColor(1);
  q_IMnpipi_wSid_n_Sm_px_clone->SetMarkerStyle(20);
  q_IMnpipi_wSid_n_Sm_px_clone->Draw("E");
  q_IMnpipi_woK0_wSid_n_Sm_px_clone->SetMarkerStyle(21);
  q_IMnpipi_woK0_wSid_n_Sm_px_clone->SetMarkerColor(3);
  q_IMnpipi_woK0_wSid_n_Sm_px_clone->Draw("HEsame");

  
  TCanvas *cq_IMnpipi_wSid_n_px_Sp[ngap];
  TH1D* q_IMnpipi_wSid_n_Sp_side_low_px[ngap];
  TH1D* q_IMnpipi_wSid_n_Sp_side_high_px[ngap];
  TH1D* q_IMnpipi_wSid_nSp_side_px_sum[ngap];
  for(int igap=0;igap<ngap;igap++){
    cq_IMnpipi_wSid_n_px_Sp[igap]  = new TCanvas(Form("cq_IMnpipi_wSid_n_px_Sp_%d",igap),Form("q_IMnpipi_wSid_n_px_Sp_%d",igap)); 
    cq_IMnpipi_wSid_n_px_Sp[igap]->Divide(2,2); 
    cq_IMnpipi_wSid_n_px_Sp[igap]->cd(1); 
    IMnpim_IMnpip_dE_n_Sp->SetMaximum(IMnpim_IMnpip_dE_n->GetMaximum());
    IMnpim_IMnpip_dE_n_Sp->Draw("colz");
    IMnpim_IMnpip_dE_n_Sp_side[LOWside][igap]->Draw("colsame");
    IMnpim_IMnpip_dE_n_Sp_side[HIGHside][igap]->Draw("colsame");

    cq_IMnpipi_wSid_n_px_Sp[igap]->cd(2); 
    //q_IMnpipi_woK0_wSid_n_Sp_px->SetMaximum(q_IMnpipi_woK0_wSid_n_Sm_px->GetMaximum());
    if(qvalcutflag==1 && !SimSpmode && !SimSmmode) q_IMnpipi_wSid_n_Sp_px->SetMaximum(160);
    if(qvalcutflag==2 && !SimSpmode && !SimSmmode) q_IMnpipi_wSid_n_Sp_px->SetMaximum(260);
    q_IMnpipi_wSid_n_Sp_px->Draw("HE");
  
    q_IMnpipi_wSid_n_Sp_side_low_px[igap] = (TH1D*) q_IMnpipi_wSid_n_Sp_side[sidebandtype][LOWside][igap]->ProjectionX();
    q_IMnpipi_wSid_n_Sp_side_low_px[igap]->SetLineColor(kCyan);
    q_IMnpipi_wSid_n_Sp_side_low_px[igap]->SetMarkerStyle(20);
    q_IMnpipi_wSid_n_Sp_side_low_px[igap]->SetMarkerColor(kCyan);
    // q_IMnpipi_woK0_wSid_n_Sp_side_low_px[igap]->Draw("Esame");
  
    //q_IMnpipi_woK0_wSid_n_Sp_side_low_px[igap]->Draw("E");
    q_IMnpipi_wSid_n_Sp_side_high_px[igap] = (TH1D*) q_IMnpipi_wSid_n_Sp_side[sidebandtype][HIGHside][igap]->ProjectionX();
    q_IMnpipi_wSid_n_Sp_side_high_px[igap]->SetLineColor(kCyan+2);
    q_IMnpipi_wSid_n_Sp_side_high_px[igap]->SetMarkerStyle(21);
    //q_IMnpipi_woK0_wSid_n_Sp_side_high_px[igap]->Draw("Esame");
  
    q_IMnpipi_wSid_nSp_side_px_sum[igap] = (TH1D*) q_IMnpipi_wSid_n_Sp_side_low_px[igap]->Clone();
    q_IMnpipi_wSid_nSp_side_px_sum[igap]->Add(q_IMnpipi_wSid_n_Sp_side_high_px[igap]);
    q_IMnpipi_wSid_nSp_side_px_sum[igap]->SetLineColor(kCyan+4);
    q_IMnpipi_wSid_nSp_side_px_sum[igap]->SetMarkerStyle(22);
    q_IMnpipi_wSid_nSp_side_px_sum[igap]->SetMarkerColor(kCyan+4);
    q_IMnpipi_wSid_nSp_side_px_sum[igap]->Draw("Esame");
    cq_IMnpipi_wSid_n_px_Sp[igap]->cd(3); 
    q_IMnpipi_wSid_nSp_side_px_sum[igap]->Draw("E");
    q_IMnpipi_wSid_n_Sp_side_low_px[igap]->Draw("Esame");
    q_IMnpipi_wSid_n_Sp_side_high_px[igap]->Draw("Esame");
    if(SimSpmode || SimSmmode) break;
  }

  TCanvas *cq_IMnpipi_woK0_wSid_n_px_Sp[ngap];
  TH1D* q_IMnpipi_woK0_wSid_n_Sp_side_low_px[ngap];
  TH1D* q_IMnpipi_woK0_wSid_n_Sp_side_high_px[ngap];
  TH1D* q_IMnpipi_woK0_wSid_nSp_side_px_sum[ngap];
  for(int igap=0;igap<ngap;igap++){
    cq_IMnpipi_woK0_wSid_n_px_Sp[igap]  = new TCanvas(Form("cq_IMnpipi_woK0_wSid_n_px_Sp_%d",igap),Form("q_IMnpipi_woK0_wSid_n_px_Sp_%d",igap)); 
    cq_IMnpipi_woK0_wSid_n_px_Sp[igap]->Divide(2,2); 
    cq_IMnpipi_woK0_wSid_n_px_Sp[igap]->cd(1); 
    IMnpim_IMnpip_dE_woK0_n_Sp->SetMaximum(IMnpim_IMnpip_dE_woK0_n->GetMaximum());
    IMnpim_IMnpip_dE_woK0_n_Sp->Draw("colz");
    IMnpim_IMnpip_dE_woK0_n_Sp_side[LOWside][igap]->Draw("colsame");
    IMnpim_IMnpip_dE_woK0_n_Sp_side[HIGHside][igap]->Draw("colsame");

    cq_IMnpipi_woK0_wSid_n_px_Sp[igap]->cd(2); 
    //q_IMnpipi_woK0_wSid_n_Sp_px->SetMaximum(q_IMnpipi_woK0_wSid_n_Sm_px->GetMaximum());
    if(qvalcutflag==1 && !SimSpmode && !SimSmmode) q_IMnpipi_woK0_wSid_n_Sp_px->SetMaximum(160);
    if(qvalcutflag==2 && !SimSpmode && !SimSmmode) q_IMnpipi_woK0_wSid_n_Sp_px->SetMaximum(260);
    q_IMnpipi_woK0_wSid_n_Sp_px->Draw("HE");
  
    q_IMnpipi_woK0_wSid_n_Sp_side_low_px[igap] = (TH1D*) q_IMnpipi_woK0_wSid_n_Sp_side[sidebandtype][LOWside][igap]->ProjectionX();
    q_IMnpipi_woK0_wSid_n_Sp_side_low_px[igap]->SetLineColor(kCyan);
    q_IMnpipi_woK0_wSid_n_Sp_side_low_px[igap]->SetMarkerStyle(20);
    q_IMnpipi_woK0_wSid_n_Sp_side_low_px[igap]->SetMarkerColor(kCyan);
    // q_IMnpipi_woK0_wSid_n_Sp_side_low_px[igap]->Draw("Esame");
  
    //q_IMnpipi_woK0_wSid_n_Sp_side_low_px[igap]->Draw("E");
    q_IMnpipi_woK0_wSid_n_Sp_side_high_px[igap] = (TH1D*) q_IMnpipi_woK0_wSid_n_Sp_side[sidebandtype][HIGHside][igap]->ProjectionX();
    q_IMnpipi_woK0_wSid_n_Sp_side_high_px[igap]->SetLineColor(kCyan+2);
    q_IMnpipi_woK0_wSid_n_Sp_side_high_px[igap]->SetMarkerStyle(21);
    //q_IMnpipi_woK0_wSid_n_Sp_side_high_px[igap]->Draw("Esame");
  
    q_IMnpipi_woK0_wSid_nSp_side_px_sum[igap] = (TH1D*) q_IMnpipi_woK0_wSid_n_Sp_side_low_px[igap]->Clone();
    q_IMnpipi_woK0_wSid_nSp_side_px_sum[igap]->Add(q_IMnpipi_woK0_wSid_n_Sp_side_high_px[igap]);
    q_IMnpipi_woK0_wSid_nSp_side_px_sum[igap]->SetLineColor(kCyan+4);
    q_IMnpipi_woK0_wSid_nSp_side_px_sum[igap]->SetMarkerStyle(22);
    q_IMnpipi_woK0_wSid_nSp_side_px_sum[igap]->SetMarkerColor(kCyan+4);
    q_IMnpipi_woK0_wSid_nSp_side_px_sum[igap]->Draw("Esame");
    cq_IMnpipi_woK0_wSid_n_px_Sp[igap]->cd(3); 
    q_IMnpipi_woK0_wSid_nSp_side_px_sum[igap]->Draw("E");
    q_IMnpipi_woK0_wSid_n_Sp_side_low_px[igap]->Draw("Esame");
    q_IMnpipi_woK0_wSid_n_Sp_side_high_px[igap]->Draw("Esame");
    if(SimSpmode || SimSmmode) break;
  }
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
  q_IMnpipi_woK0_wSid_n_Sp_sub->Add(q_IMnpipi_woK0_wSid_n_Sp_side[sidebandtype][LOWside][0],-1);
  q_IMnpipi_woK0_wSid_n_Sp_sub->Add(q_IMnpipi_woK0_wSid_n_Sp_side[sidebandtype][HIGHside][0],-1);
  q_IMnpipi_woK0_wSid_n_Sp_sub->SetMinimum(0);
  q_IMnpipi_woK0_wSid_n_Sp_sub->Draw("colz");
  
  TCanvas *cq_IMnpipi_wSid_n_Sm_px[ngap];
  TH1D* q_IMnpipi_wSid_n_Sm_side_low_px[ngap];
  TH1D* q_IMnpipi_wSid_n_Sm_side_high_px[ngap];
  TH1D* q_IMnpipi_wSid_n_Sm_side_px_sum[ngap];
  for(int igap=0;igap<ngap;igap++){
    cq_IMnpipi_wSid_n_Sm_px[igap] = new TCanvas(Form("cq_IMnpipi_wSid_n_Sm_px_%d",igap),Form("q_IMnpipi_wSid_n_Sm_px_%d",igap)); 
    cq_IMnpipi_wSid_n_Sm_px[igap]->Divide(2,2);  
    cq_IMnpipi_wSid_n_Sm_px[igap]->cd(1);  
    IMnpim_IMnpip_dE_n_Sm->SetMaximum(IMnpim_IMnpip_dE_n->GetMaximum());
    IMnpim_IMnpip_dE_n_Sm->GetXaxis()->SetRangeUser(0,1.7);
    IMnpim_IMnpip_dE_n_Sm->GetYaxis()->SetRangeUser(0,1.7);
    IMnpim_IMnpip_dE_n_Sm->Draw("colz");
    IMnpim_IMnpip_dE_n_Sm_side[LOWside][igap]->Draw("colsame");
    IMnpim_IMnpip_dE_n_Sm_side[HIGHside][igap]->Draw("colsame");
 
    cq_IMnpipi_wSid_n_Sm_px[igap]->cd(2);  
   //q_IMnpipi_woK0_wSid_n_Sm_px->SetMaximum(q_IMnpipi_woK0_wSid_n_px->GetMaximum());
    if(qvalcutflag==1) q_IMnpipi_wSid_n_Sm_px->SetMaximum(160);
    if(qvalcutflag==2) q_IMnpipi_wSid_n_Sm_px->SetMaximum(260);
    q_IMnpipi_wSid_n_Sm_px->Draw("EH");
    q_IMnpipi_wSid_n_Sm_side_low_px[igap] = (TH1D*) q_IMnpipi_wSid_n_Sm_side[sidebandtype][LOWside][igap]->ProjectionX();
    q_IMnpipi_wSid_n_Sm_side_low_px[igap]->SetLineColor(kPink+1);
    q_IMnpipi_wSid_n_Sm_side_low_px[igap]->SetMarkerStyle(20);
    //q_IMnpipi_woK0_wSid_n_Sm_side_low_px[igap]->Draw("Esame");
  
    q_IMnpipi_wSid_n_Sm_side_high_px[igap] = (TH1D*) q_IMnpipi_wSid_n_Sm_side[sidebandtype][HIGHside][igap]->ProjectionX();
    q_IMnpipi_wSid_n_Sm_side_high_px[igap]->SetLineColor(kPink+3);
    q_IMnpipi_wSid_n_Sm_side_high_px[igap]->SetMarkerStyle(21);
    q_IMnpipi_wSid_n_Sm_side_high_px[igap]->SetMarkerColor(kPink+3);
    //q_IMnpipi_woK0_wSid_n_Sm_side_high_px[igap]->Draw("Esame");
  
    q_IMnpipi_wSid_n_Sm_side_px_sum[igap] = (TH1D*) q_IMnpipi_wSid_n_Sm_side_low_px[igap]->Clone();
    q_IMnpipi_wSid_n_Sm_side_px_sum[igap]->Add(q_IMnpipi_wSid_n_Sm_side_high_px[igap]);
    q_IMnpipi_wSid_n_Sm_side_px_sum[igap]->SetLineColor(kPink+10);
    q_IMnpipi_wSid_n_Sm_side_px_sum[igap]->SetMarkerStyle(22);
    q_IMnpipi_wSid_n_Sm_side_px_sum[igap]->SetMarkerColor(kPink+10);
    q_IMnpipi_wSid_n_Sm_side_px_sum[igap]->Draw("Esame");
    cq_IMnpipi_wSid_n_Sm_px[igap]->cd(3);  
    q_IMnpipi_wSid_n_Sm_side_px_sum[igap]->Draw("E");
    q_IMnpipi_wSid_n_Sm_side_low_px[igap]->Draw("Esame");
    q_IMnpipi_wSid_n_Sm_side_high_px[igap]->Draw("Esame");
    if(SimSpmode || SimSmmode) break;
  }
  
  //overlay signal and sideband projected 1d histo to IM
  TCanvas *cq_IMnpipi_woK0_wSid_n_Sm_px[ngap];
  TH1D* q_IMnpipi_woK0_wSid_n_Sm_side_low_px[ngap];
  TH1D* q_IMnpipi_woK0_wSid_n_Sm_side_high_px[ngap];
  TH1D* q_IMnpipi_woK0_wSid_n_Sm_side_px_sum[ngap];
  for(int igap=0;igap<ngap;igap++){
    cq_IMnpipi_woK0_wSid_n_Sm_px[igap] = new TCanvas(Form("cq_IMnpipi_woK0_wSid_n_Sm_px_%d",igap),Form("q_IMnpipi_woK0_wSid_n_Sm_px_%d",igap)); 
    cq_IMnpipi_woK0_wSid_n_Sm_px[igap]->Divide(2,2);  
    cq_IMnpipi_woK0_wSid_n_Sm_px[igap]->cd(1);  
    IMnpim_IMnpip_dE_woK0_n_Sm->SetMaximum(IMnpim_IMnpip_dE_woK0_n->GetMaximum());
    IMnpim_IMnpip_dE_woK0_n_Sm->GetXaxis()->SetRangeUser(0,1.7);
    IMnpim_IMnpip_dE_woK0_n_Sm->GetYaxis()->SetRangeUser(0,1.7);
    IMnpim_IMnpip_dE_woK0_n_Sm->Draw("colz");
    IMnpim_IMnpip_dE_woK0_n_Sm_side[LOWside][igap]->Draw("colsame");
    IMnpim_IMnpip_dE_woK0_n_Sm_side[HIGHside][igap]->Draw("colsame");
 
    cq_IMnpipi_woK0_wSid_n_Sm_px[igap]->cd(2);  
   //q_IMnpipi_woK0_wSid_n_Sm_px->SetMaximum(q_IMnpipi_woK0_wSid_n_px->GetMaximum());
    if(qvalcutflag==1) q_IMnpipi_woK0_wSid_n_Sm_px->SetMaximum(160);
    if(qvalcutflag==2) q_IMnpipi_woK0_wSid_n_Sm_px->SetMaximum(260);
    q_IMnpipi_woK0_wSid_n_Sm_px->Draw("EH");
    q_IMnpipi_woK0_wSid_n_Sm_side_low_px[igap] = (TH1D*) q_IMnpipi_woK0_wSid_n_Sm_side[sidebandtype][LOWside][igap]->ProjectionX();
    q_IMnpipi_woK0_wSid_n_Sm_side_low_px[igap]->SetLineColor(kPink+1);
    q_IMnpipi_woK0_wSid_n_Sm_side_low_px[igap]->SetMarkerStyle(20);
    //q_IMnpipi_woK0_wSid_n_Sm_side_low_px[igap]->Draw("Esame");
  
    q_IMnpipi_woK0_wSid_n_Sm_side_high_px[igap] = (TH1D*) q_IMnpipi_woK0_wSid_n_Sm_side[sidebandtype][HIGHside][igap]->ProjectionX();
    q_IMnpipi_woK0_wSid_n_Sm_side_high_px[igap]->SetLineColor(kPink+3);
    q_IMnpipi_woK0_wSid_n_Sm_side_high_px[igap]->SetMarkerStyle(21);
    q_IMnpipi_woK0_wSid_n_Sm_side_high_px[igap]->SetMarkerColor(kPink+3);
    //q_IMnpipi_woK0_wSid_n_Sm_side_high_px[igap]->Draw("Esame");
  
    q_IMnpipi_woK0_wSid_n_Sm_side_px_sum[igap] = (TH1D*) q_IMnpipi_woK0_wSid_n_Sm_side_low_px[igap]->Clone();
    q_IMnpipi_woK0_wSid_n_Sm_side_px_sum[igap]->Add(q_IMnpipi_woK0_wSid_n_Sm_side_high_px[igap]);
    q_IMnpipi_woK0_wSid_n_Sm_side_px_sum[igap]->SetLineColor(kPink+10);
    q_IMnpipi_woK0_wSid_n_Sm_side_px_sum[igap]->SetMarkerStyle(22);
    q_IMnpipi_woK0_wSid_n_Sm_side_px_sum[igap]->SetMarkerColor(kPink+10);
    q_IMnpipi_woK0_wSid_n_Sm_side_px_sum[igap]->Draw("Esame");
    cq_IMnpipi_woK0_wSid_n_Sm_px[igap]->cd(3);  
    q_IMnpipi_woK0_wSid_n_Sm_side_px_sum[igap]->Draw("E");
    q_IMnpipi_woK0_wSid_n_Sm_side_low_px[igap]->Draw("Esame");
    q_IMnpipi_woK0_wSid_n_Sm_side_high_px[igap]->Draw("Esame");
    if(SimSpmode || SimSmmode) break;
  }
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
  q_IMnpipi_woK0_wSid_n_Sm_sub->Add(q_IMnpipi_woK0_wSid_n_Sm_side[sidebandtype][LOWside][0],-1);
  q_IMnpipi_woK0_wSid_n_Sm_sub->Add(q_IMnpipi_woK0_wSid_n_Sm_side[sidebandtype][HIGHside][0],-1);
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
  
  TCanvas *cq_IMnpipi_wSid_n_side[ngap];
  for(int igap=0;igap<ngap;igap++){
    cq_IMnpipi_wSid_n_side[igap] = new TCanvas(Form("cq_IMnpipi_wSid_n_side_%d",igap),Form("q_IMnpipi_wSid_n_side_%d",igap));
    cq_IMnpipi_wSid_n_side[igap]->cd();
    q_IMnpipi_wSid_n_side[igap]->SetMaximum(q_IMnpipi_wSid_n->GetMaximum());
    q_IMnpipi_wSid_n_side[igap]->Draw("colz");
  }

  TCanvas *cq_IMnpipi_woK0_wSid_n = new TCanvas("cq_IMnpipi_woK0_wSid_n","q_IMnpipi_woK0_wSid_n");
  cq_IMnpipi_woK0_wSid_n->cd();
  q_IMnpipi_woK0_wSid_n->Draw("colz");
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
  cdE_betainv_fid_px->cd();
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
  TLegend *ldE_betainv_fid_px = new TLegend(0.6,0.7,0.9,0.9);
  ldE_betainv_fid_px->AddEntry(h1_nocut,"< 1 MeVee","f");
  ldE_betainv_fid_px->AddEntry(h1_2mevcut,"< 2 MeVee","f");
  ldE_betainv_fid_px->AddEntry(h1_4mevcut,"< 4 MeVee","f");
  ldE_betainv_fid_px->AddEntry(h1_6mevcut,"< 6 MeVee","f");
  ldE_betainv_fid_px->Draw();

  TCanvas *cdE_nmom_fid_beta = new TCanvas("cdE_nmom_fid_beta","dE_nmom_fid_beta");
  cdE_nmom_fid_beta->cd();
  dE_nmom_fid_beta->Draw("colz");
  cdE_nmom_fid_beta->SetLogz();

  TCanvas *cdE_MMom_fid_beta = new TCanvas("cdE_MMom_fid_beta","dE_MMom_fid_beta");
  cdE_MMom_fid_beta->cd();
  dE_MMom_fid_beta->Draw("colz");

  TCanvas *cdE_MMom_fid_beta_woK0 = new TCanvas("cdE_MMom_fid_beta_woK0","dE_MMom_fid_beta_woK0");
  cdE_MMom_fid_beta_woK0->cd();
  dE_MMom_fid_beta_woK0->Draw("colz");
  
  TCanvas *cdE_MMass_fid_beta_woK0 = new TCanvas("cdE_MMass_fid_beta_woK0","dE_MMass_fid_beta_woK0");
  cdE_MMass_fid_beta_woK0->cd();
  dE_MMass_fid_beta_woK0->Draw("colz");

  TCanvas *cdE_MMass_fid_beta_woK0_wSid = new TCanvas("cdE_MMass_fid_beta_woK0_wSid","dE_MMass_fid_beta_woK0_wSid");
  cdE_MMass_fid_beta_woK0_wSid->cd();
  dE_MMass_fid_beta_woK0_wSid->Draw("colz");
   


  TCanvas *cnmom_IMnpip_dE_n = new TCanvas("cnmom_IMnpip_dE_n","nmom_IMnpip_dE_n");
  cnmom_IMnpip_dE_n->cd();
  nmom_IMnpip_dE_n->Draw("colz");
  
  TCanvas *cnmom_IMnpim_dE_n = new TCanvas("cnmom_IMnpim_dE_n","nmom_IMnpim_dE_n");
  cnmom_IMnpim_dE_n->cd();
  nmom_IMnpim_dE_n->Draw("colz");
  
  TCanvas *cnmom_IMnpip_dE_woK0_n = new TCanvas("cnmom_IMnpip_dE_woK0_n","nmom_IMnpip_dE_woK0_n");
  cnmom_IMnpip_dE_woK0_n->cd();
  nmom_IMnpip_dE_woK0_n->Draw("colz");
  
  TCanvas *cnmom_IMnpim_dE_woK0_n = new TCanvas("cnmom_IMnpim_dE_woK0_n","nmom_IMnpim_dE_woK0_n");
  cnmom_IMnpim_dE_woK0_n->cd();
  nmom_IMnpim_dE_woK0_n->Draw("colz");

  //TCanvas *cIMnpim_IMnpip_dE_woK0 = new TCanvas(Form("cIMnpim_IMnpip_dE_woK0"),"IMnpim_IMnpip_dE_woK0");
  //cIMnpim_IMnpip_dE_woK0->cd();
  //IMnpim_IMnpip_dE_woK0->Draw("colz");

  TCanvas *cnmom_IMnpipi_wK0_n = new TCanvas(Form("cnmom_IMnpipi_wK0_n"),"nmom_IMnpipi_wK0_n");
  cnmom_IMnpipi_wK0_n->cd();
  nmom_IMnpipi_wK0_n->Draw("colz");

  TCanvas *cnmom_IMnpipi_wK0_n_px = new TCanvas(Form("cnmom_IMnpipi_wK0_n_px"),"nmom_IMnpipi_wK0_n_px");
  cnmom_IMnpipi_wK0_n_px->cd();
  TH1D* nmom_IMnpipi_wK0_n_px = nmom_IMnpipi_wK0_n->ProjectionX();
  nmom_IMnpipi_wK0_n_px->SetFillColor(kViolet);
  nmom_IMnpipi_wK0_n_px->Draw("HE");
  
  TCanvas *cnmom_IMnpipi_wK0_n_py = new TCanvas(Form("cnmom_IMnpipi_wK0_n_py"),"nmom_IMnpipi_wK0_n_py");
  cnmom_IMnpipi_wK0_n_py->cd();
  TH1D* nmom_IMnpipi_wK0_n_py = nmom_IMnpipi_wK0_n->ProjectionY();
  nmom_IMnpipi_wK0_n_py->SetFillColor(kViolet);
  //nmom_IMnpipi_wK0_n_py->SetMaximum(3300);
  nmom_IMnpipi_wK0_n_py->Draw("HE");
  
  TCanvas *cnmom_IMnpipi_woK0_wSid_n = new TCanvas(Form("cnmom_IMnpipi_woK0_wSid_n"),"nmom_IMnpipi_woK0_wSid_n");
  cnmom_IMnpipi_woK0_wSid_n->cd();
  nmom_IMnpipi_woK0_wSid_n->Draw("colz");

  TCanvas *cnmom_IMnpipi_woK0_wSid_n_px = new TCanvas(Form("cnmom_IMnpipi_woK0_wSid_n_px"),"nmom_IMnpipi_woK0_wSid_n_px");
  cnmom_IMnpipi_woK0_wSid_n_px->cd();
  TH1D* nmom_IMnpipi_woK0_wSid_n_px = nmom_IMnpipi_woK0_wSid_n->ProjectionX();
  nmom_IMnpipi_woK0_wSid_n_px->Draw();
  
  TCanvas *cnmom_IMnpipi_woK0_wSid_n_py = new TCanvas(Form("cnmom_IMnpipi_woK0_wSid_n_py"),"nmom_IMnpipi_woK0_wSid_n_py");
  cnmom_IMnpipi_woK0_wSid_n_py->cd();
  TH1D* nmom_IMnpipi_woK0_wSid_n_py = nmom_IMnpipi_woK0_wSid_n->ProjectionY();
  nmom_IMnpipi_woK0_wSid_n_py->Draw();


  TCanvas *cnmom_IMnpipi_woK0_wSid_n_Sp = new TCanvas(Form("cnmom_IMnpipi_woK0_wSid_n_Sp"),"nmom_IMnpipi_woK0_wSid_n_Sp");
  cnmom_IMnpipi_woK0_wSid_n_Sp->cd();
  nmom_IMnpipi_woK0_wSid_n_Sp->Draw("colz");

  TCanvas *cnmom_IMnpipi_woK0_wSid_n_Sp_px = new TCanvas(Form("cnmom_IMnpipi_woK0_wSid_n_Sp_px"),"nmom_IMnpipi_woK0_wSid_n_Sp_px");
  cnmom_IMnpipi_woK0_wSid_n_Sp_px->cd();
  TH1D* nmom_IMnpipi_woK0_wSid_n_Sp_px = nmom_IMnpipi_woK0_wSid_n_Sp->ProjectionX();
  nmom_IMnpipi_woK0_wSid_n_Sp_px->SetLineColor(2);
  nmom_IMnpipi_woK0_wSid_n_Sp_px->Draw("HE");
  
  TCanvas *cnmom_IMnpipi_woK0_wSid_n_Sp_py = new TCanvas(Form("cnmom_IMnpipi_woK0_wSid_n_Sp_py"),"nmom_IMnpipi_woK0_wSid_n_Sp_py");
  cnmom_IMnpipi_woK0_wSid_n_Sp_py->cd();
  TH1D* nmom_IMnpipi_woK0_wSid_n_Sp_py = nmom_IMnpipi_woK0_wSid_n_Sp->ProjectionY();
  nmom_IMnpipi_woK0_wSid_n_Sp_py->SetLineColor(2);
  nmom_IMnpipi_woK0_wSid_n_Sp_py->Draw("HE");

  TCanvas *cnmom_IMnpipi_woK0_wSid_n_Sm = new TCanvas(Form("cnmom_IMnpipi_woK0_wSid_n_Sm"),"nmom_IMnpipi_woK0_wSid_n_Sm");
  cnmom_IMnpipi_woK0_wSid_n_Sm->cd();
  nmom_IMnpipi_woK0_wSid_n_Sm->Draw("colz");

  TCanvas *cnmom_IMnpipi_woK0_wSid_n_Sm_px = new TCanvas(Form("cnmom_IMnpipi_woK0_wSid_n_Sm_px"),"nmom_IMnpipi_woK0_wSid_n_Sm_px");
  cnmom_IMnpipi_woK0_wSid_n_Sm_px->cd();
  TH1D* nmom_IMnpipi_woK0_wSid_n_Sm_px = nmom_IMnpipi_woK0_wSid_n_Sm->ProjectionX();
  nmom_IMnpipi_woK0_wSid_n_Sm_px->SetLineColor(3);
  nmom_IMnpipi_woK0_wSid_n_Sm_px->Draw("HE");
  
  TCanvas *cnmom_IMnpipi_woK0_wSid_n_Sm_py = new TCanvas(Form("cnmom_IMnpipi_woK0_wSid_n_Sm_py"),"nmom_IMnpipi_woK0_wSid_n_Sm_py");
  cnmom_IMnpipi_woK0_wSid_n_Sm_py->cd();
  TH1D* nmom_IMnpipi_woK0_wSid_n_Sm_py = nmom_IMnpipi_woK0_wSid_n_Sm->ProjectionY();
  nmom_IMnpipi_woK0_wSid_n_Sm_py->SetLineColor(3);
  nmom_IMnpipi_woK0_wSid_n_Sm_py->Draw("HE");


  TCanvas *cnmom_IMnpipi_woK0_wSid_n_SpSm_py = new TCanvas("cnmom_IMnpipi_woK0_wSid_n_SpSm_py","cnmom_IMnpipi_woK0_wSid_n_SpSm_py");
  cnmom_IMnpipi_woK0_wSid_n_SpSm_py->cd();
  nmom_IMnpipi_woK0_wSid_n_Sm_py->Draw("HE");
  nmom_IMnpipi_woK0_wSid_n_Sp_py->Draw("HEsame");
  
  TCanvas *cnmom_cosnlab_K0_n = new TCanvas(Form("cnmom_cosnlab_K0_n"),"nmom_cosnlab_K0_n");
  cnmom_cosnlab_K0_n->cd();
  nmom_cosnlab_K0_n->Draw("colz");
  
  TCanvas *cnmom_IMpippim_n = new TCanvas(Form("cnmom_IMpippim_n"),"nmom_IMpippim_n");
  cnmom_IMpippim_n->cd();
  nmom_IMpippim_n->Draw("colz");

  TCanvas *cnmom_IMpippim_wSid_n = new TCanvas(Form("cnmom_IMpippim_wSid_n"),"nmom_IMpippim_wSid_n");
  cnmom_IMpippim_wSid_n->cd();
  nmom_IMpippim_wSid_n->Draw("colz");

  TCanvas *cnmom_IMpippim_wSid_n_py = new TCanvas(Form("cnmom_IMpippim_wSid_n_py"),"nmom_IMpippim_wSid_n_py");
  cnmom_IMpippim_wSid_n_py->cd();
  TH1D* nmom_IMpippim_wSid_n_py = (TH1D*)nmom_IMpippim_wSid_n->ProjectionY();
  // mom_IMpippim_wSid_n_py->SetFillColor(kViolet);
  nmom_IMpippim_wSid_n_py->Draw("HE");
  TH1D* nmom_IMpippim_wSid_n_py_cut = (TH1D*)nmom_IMpippim_wSid_n->ProjectionY("nmom_IMpippim_wSid_n_py_cut",
                                                                               nmom_IMpippim_wSid_n->GetXaxis()->FindBin(anacuts::pipi_MIN),
                                                                               nmom_IMpippim_wSid_n->GetXaxis()->FindBin(anacuts::pipi_MAX));
  nmom_IMpippim_wSid_n_py_cut->SetFillColor(kViolet);
  nmom_IMpippim_wSid_n_py_cut->Draw("HEsame");

  TCanvas *cmnmom_IMpippim_n = new TCanvas(Form("cmnmom_IMpippim_n"),"mnmom_IMpippim_n");
  cmnmom_IMpippim_n->cd();
  mnmom_IMpippim_n->Draw("colz");
  
  TCanvas *cmnmom_IMpippim_n_px = new TCanvas(Form("cmnmom_IMpippim_n_px"),"mnmom_IMpippim_n_px");
  cmnmom_IMpippim_n_px->cd();
  TH1D* mnmom_IMpippim_n_px = (TH1D*)mnmom_IMpippim_n->ProjectionX();
  mnmom_IMpippim_n_px->Draw("HE");
  TH1D* mnmom_IMpippim_n_px_cut= (TH1D*)mnmom_IMpippim_n_px->Clone();
  mnmom_IMpippim_n_px_cut->GetXaxis()->SetRangeUser(anacuts::pipi_MIN,anacuts::pipi_MAX);
  mnmom_IMpippim_n_px_cut->SetFillColor(kViolet);
  mnmom_IMpippim_n_px_cut->Draw("HEsame");

  TCanvas *cmnmom_IMpippim_n_py = new TCanvas(Form("cmnmom_IMpippim_n_py"),"mnmom_IMpippim_n_py");
  cmnmom_IMpippim_n_py->cd();
  TH1D* mnmom_IMpippim_n_py = (TH1D*)mnmom_IMpippim_n->ProjectionY();
  mnmom_IMpippim_n_py->Draw("HE");
  
  //TCanvas *cmnmom_IMpippim_n_py_cut = new TCanvas(Form("cmnmom_IMpippim_n_py_cut"),"mnmom_IMpippim_n_py_cut");
  //cmnmom_IMpippim_n_py_cut->cd();
  TH1D* mnmom_IMpippim_n_py_cut = (TH1D*)mnmom_IMpippim_n->ProjectionY("mnmom_IMpippim_n_py_cut",
                                                                        mnmom_IMpippim_n->GetXaxis()->FindBin(anacuts::pipi_MIN),
                                                                        mnmom_IMpippim_n->GetXaxis()->FindBin(anacuts::pipi_MAX));
  //mnmom_IMpippim_n_py_cut->SetLineColor(2);
  mnmom_IMpippim_n_py_cut->SetFillColor(kViolet);
  mnmom_IMpippim_n_py_cut->Draw("HEsame");
  TLegend *lmnmom_IMpippim_n_py  = new TLegend(0.15,0.7,0.45,0.9);
  lmnmom_IMpippim_n_py->AddEntry(mnmom_IMpippim_n_py,Form("missing mom. "),"f");
  lmnmom_IMpippim_n_py->AddEntry(mnmom_IMpippim_n_py_cut,Form("missing mom. (K0 region selected)"),"f");
  lmnmom_IMpippim_n_py->Draw();

  TCanvas *cmnmom_IMpippim_wSid_n = new TCanvas(Form("cmnmom_IMpippim_wSid_n"),"mnmom_IMpippim_wSid_n");
  cmnmom_IMpippim_wSid_n->cd();
  mnmom_IMpippim_wSid_n->Draw("colz");
  
  //TCanvas *cMMnmiss_IMnpipi_woK0_wSid_px = new TCanvas("cMMnmiss_IMnpipi_woK0_wSid_px","MMnmiss_IMnpipi_woK0_wSid_px"); 
  //cMMnmiss_IMnpipi_woK0_wSid_px->cd();
  //MMnmiss_IMnpipi_woK0_wSid->ProjectionX()->Draw("");
  //TCanvas *cMMnmiss_IMnpipi_woK0_wSid_py = new TCanvas("cMMnmiss_IMnpipi_woK0_wSid_py","MMnmiss_IMnpipi_woK0_wSid_py"); 
  //cMMnmiss_IMnpipi_woK0_wSid_py->cd();
  //MMnmiss_IMnpipi_woK0_wSid->ProjectionY()->Draw("");
  
  TCanvas *cMMnmiss_IMnpipi_woK0_wSid_SpSm_py = new TCanvas("cMMnmiss_IMnpipi_woK0_wSid_SpSm_py","MMnmiss_IMnpipi_woK0_wSid_SpSm_py"); 
  cMMnmiss_IMnpipi_woK0_wSid_SpSm_py->cd();
  TH1D* MMnmiss_IMnpipi_woK0_wSid_py = (TH1D*)MMnmiss_IMnpipi_woK0_wSid->ProjectionY();
  MMnmiss_IMnpipi_woK0_wSid_py->Draw("HE");
  TH1D* MMnmiss_IMnpipi_woK0_wSid_Sp_py = (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sp->ProjectionY();
  MMnmiss_IMnpipi_woK0_wSid_Sp_py->SetLineColor(2);
  MMnmiss_IMnpipi_woK0_wSid_Sp_py->Draw("HEsame");
  TH1D* MMnmiss_IMnpipi_woK0_wSid_Sm_py = (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sm->ProjectionY();
  MMnmiss_IMnpipi_woK0_wSid_Sm_py->SetLineColor(3);
  MMnmiss_IMnpipi_woK0_wSid_Sm_py->Draw("HEsame");


  TCanvas *cMMnmiss_IMnpipi_wSid_n = new TCanvas("cMMnmiss_IMnpipi_wSid_n","MMnmiss_IMnpipi_wSid_n");
  cMMnmiss_IMnpipi_wSid_n->cd();
  MMnmiss_IMnpipi_wSid_n->Draw("colz");
  
  TCanvas *cMMnmiss_IMnpipi_woK0_wSid = new TCanvas("cMMnmiss_IMnpipi_woK0_wSid","MMnmiss_IMnpipi_woK0_wSid");
  cMMnmiss_IMnpipi_woK0_wSid->cd();
  MMnmiss_IMnpipi_woK0_wSid->Draw("colz");

  TCanvas *cq_IMpippim_n = new TCanvas("cq_IMpippim_n","q_IMpippim_n");
  cq_IMpippim_n->cd();
  q_IMpippim_n->GetYaxis()->SetRangeUser(0,1.0);
  q_IMpippim_n->Draw("colz");

  TCanvas *cq_IMpippim_n_px = new TCanvas("cq_IMpippim_n_px","q_IMpippim_n_px");
  cq_IMpippim_n_px->cd();
  TH1D* q_IMpippim_n_px = q_IMpippim_n->ProjectionX();
  q_IMpippim_n_px->GetXaxis()->SetRangeUser(0.4,0.6);
  q_IMpippim_n_px->GetYaxis()->SetTitle("Counts / 2 MeV/c^{2}");
  q_IMpippim_n_px->GetYaxis()->CenterTitle();
  q_IMpippim_n_px->Draw("EH");
  TH1D* q_IMpippim_n_px_cut = (TH1D*)q_IMpippim_n_px->Clone();
  q_IMpippim_n_px_cut->GetXaxis()->SetRangeUser(anacuts::pipi_MIN,anacuts::pipi_MAX);
  q_IMpippim_n_px_cut->SetFillColor(kViolet);
  q_IMpippim_n_px_cut->Draw("HEsame");
  
  TCanvas *cq_IMpippim_n_wSid = new TCanvas("cq_IMpippim_n_wSid","q_IMpippim_n_wSid");
  cq_IMpippim_n_wSid->cd();
  q_IMpippim_n_wSid->Draw("colz");
  
  TCanvas *cMompippim_IMnpipi_dE_wK0_n = new TCanvas("cMompippim_IMnpipi_dE_wK0_n","Mompippim_IMnpipi_dE_wK0_n");
  cMompippim_IMnpipi_dE_wK0_n->cd();
  Mompippim_IMnpipi_dE_wK0_n->Draw("colz");
  fkp->Draw("same");
  TCanvas *cMompippim_IMnpipi_dE_wK0_n_py = new TCanvas("cMompippim_IMnpipi_dE_wK0_n_py","Mompippim_IMnpipi_dE_wK0_n_py");
  cMompippim_IMnpipi_dE_wK0_n_py->cd();
  TH1D* Mompippim_IMnpipi_dE_wK0_n_py = (TH1D*)Mompippim_IMnpipi_dE_wK0_n->ProjectionY();
  Mompippim_IMnpipi_dE_wK0_n_py->SetFillColor(kViolet);
  //Mompippim_IMnpipi_dE_wK0_n_py->SetMaximum(1800);
  Mompippim_IMnpipi_dE_wK0_n_py->Draw("HE");

  
  TCanvas *cMompippim_IMnpipi_dE_wK0_n_Sp = new TCanvas("cMompippim_IMnpipi_dE_wK0_n_Sp","Mompippim_IMnpipi_dE_wK0_n_Sp");
  cMompippim_IMnpipi_dE_wK0_n_Sp->cd();
  Mompippim_IMnpipi_dE_wK0_n_Sp->Draw("colz");
  
  TCanvas *cMompippim_IMnpipi_dE_wK0_n_Sp_py = new TCanvas("cMompippim_IMnpipi_dE_wK0_n_Sp_py","Mompippim_IMnpipi_dE_wK0_n_Sp_py");
  cMompippim_IMnpipi_dE_wK0_n_Sp_py->cd();
  TH1D* Mompippim_IMnpipi_dE_wK0_n_Sp_py = (TH1D*)Mompippim_IMnpipi_dE_wK0_n_Sp->ProjectionY();
  Mompippim_IMnpipi_dE_wK0_n_Sp_py->SetFillColor(kViolet);
  Mompippim_IMnpipi_dE_wK0_n_Sp_py->Draw("HE");

  TCanvas *cMompippim_IMnpipi_dE_wK0_n_Sm = new TCanvas("cMompippim_IMnpipi_dE_wK0_n_Sm","Mompippim_IMnpipi_dE_wK0_n_Sm");
  cMompippim_IMnpipi_dE_wK0_n_Sm->cd();
  Mompippim_IMnpipi_dE_wK0_n_Sm->Draw("colz");
  
  TCanvas *cMompippim_IMnpipi_dE_wK0_n_Sm_py = new TCanvas("cMompippim_IMnpipi_dE_wK0_n_Sm_py","Mompippim_IMnpipi_dE_wK0_n_Sm_py");
  cMompippim_IMnpipi_dE_wK0_n_Sm_py->cd();
  TH1D* Mompippim_IMnpipi_dE_wK0_n_Sm_py = (TH1D*)Mompippim_IMnpipi_dE_wK0_n_Sm->ProjectionY();
  Mompippim_IMnpipi_dE_wK0_n_Sm_py->SetFillColor(kViolet);
  Mompippim_IMnpipi_dE_wK0_n_Sm_py->Draw("HE");

  /*
  TCanvas *cq_IMpippim_wSid_n_Sp = new TCanvas("cq_IMpippim_wSid_n_Sp","q_IMpippim_wSid_n_Sp");
  cq_IMpippim_wSid_n_Sp->cd();
  q_IMpippim_wSid_n_Sp->Draw("colz");

  TCanvas *cq_IMpippim_wSid_n_Sm = new TCanvas("cq_IMpippim_wSid_n_Sm","q_IMpippim_wSid_n_Sm");
  cq_IMpippim_wSid_n_Sm->cd();
  q_IMpippim_wSid_n_Sm->Draw("colz");
  */

  //TCanvas *cq_IMnpipi_wSid_n_Sm_px = new TCanvas("cq_IMnpipi_wSid_n_Sm_px","q_IMnpipi_wSid_n_Sm_px");
  //cq_IMnpipi_wSid_n_Sm_px->cd();
  //TCanvas *cq_IMnpipi_wSid_n_SpSm_px = new TCanvas("cq_IMnpipi_wSid_n_SpSm_px","q_IMpippim_wSid_n_SpSm_px");
  //cq_IMnpipi_wSid_n_SpSm_px->cd();
  //TH1D* q_IMnpipi_wSid_n_px = (TH1D*)q_IMnpipi_wSid_n->ProjectionX("q_IMnpipi_wSid_n_px");
  //q_IMnpipi_wSid_n_px->Draw("HE");
  //TH1D* q_IMnpipi_wSid_n_Sp_px = (TH1D*)q_IMnpipi_wSid_n_Sp->ProjectionX("q_IMnpipi_wSid_n_Sp_px");
  //q_IMnpipi_wSid_n_Sp_px->SetLineColor(2);
  //q_IMnpipi_wSid_n_Sp_px->Draw("Esame");
  //TH1D* q_IMnpipi_wSid_n_Sm_px = (TH1D*)q_IMnpipi_wSid_n_Sm->ProjectionX("q_IMnpipi_wSid_n_Sm_px");
  //q_IMnpipi_wSid_n_Sm_px->SetLineColor(3);
  //q_IMnpipi_wSid_n_Sm_px->Draw("HEsame");

  
  TCanvas *cIMpippim_IMnpipi_n = new TCanvas("cIMpippim_IMnpipi_n","IMpippim_IMnpipi_n");
  cIMpippim_IMnpipi_n->cd();
  IMpippim_IMnpipi_n->Draw("colz");

  TCanvas *cIMpippim_IMnpipi_n_wSid = new TCanvas("cIMpippim_IMnpipi_n_wSid","IMpippim_IMnpipi_n_wSid");
  cIMpippim_IMnpipi_n_wSid->cd();
  IMpippim_IMnpipi_n_wSid->Draw("colz");

  TCanvas *cIMpippim_IMnpip_n = new TCanvas("cIMpippim_IMnpip_n","IMpippim_IMnpip_n");
  cIMpippim_IMnpip_n->cd();
  IMpippim_IMnpip_n->Draw("colz");
  
  TCanvas *cIMpippim_IMnpim_n = new TCanvas("cIMpippim_IMnpim_n","IMpippim_IMnpim_n");
  cIMpippim_IMnpim_n->cd();
  IMpippim_IMnpim_n->Draw("colz");
  
  TCanvas *cMMnpip_MMnpim_woK0_n = new TCanvas("cMMnpip_MMnpim_woK0_n","MMnpip_MMnpim_woK0_n");
  cMMnpip_MMnpim_woK0_n->cd();
  MMnpip_MMnpim_woK0_n->Draw("colz");
  
  TCanvas *cMMnpip_MMnpim_wwoK0_n_px = new TCanvas("cMMnpip_MMnpim_wwoK0_n_px","cMMnpip_MMnpim_wwoK0_n_px");
  cMMnpip_MMnpim_wwoK0_n_px->cd();
  TH1D* MMnpip_MMnpim_n_px = (TH1D*)MMnpip_MMnpim_n->ProjectionX();
  MMnpip_MMnpim_n_px->Draw("HE");
  TH1D* MMnpip_MMnpim_woK0_n_px = (TH1D*)MMnpip_MMnpim_woK0_n->ProjectionX();
  MMnpip_MMnpim_woK0_n_px->SetLineColor(kViolet);
  MMnpip_MMnpim_woK0_n_px->Draw("HEsame");

  TCanvas *cMMnpip_MMnpim_wwoK0_n_py = new TCanvas("cMMnpip_MMnpim_wwoK0_n_py","cMMnpip_MMnpim_wwoK0_n_py");
  cMMnpip_MMnpim_wwoK0_n_py->cd();
  TH1D* MMnpip_MMnpim_n_py = (TH1D*)MMnpip_MMnpim_n->ProjectionY();
  MMnpip_MMnpim_n_py->Draw("HE");
  TH1D* MMnpip_MMnpim_woK0_n_py = (TH1D*)MMnpip_MMnpim_woK0_n->ProjectionY();
  MMnpip_MMnpim_woK0_n_py->SetLineColor(kViolet);
  MMnpip_MMnpim_woK0_n_py->Draw("HEsame");


  TCanvas *cMMnpip_MMnpim_woK0_wSid_n = new TCanvas("cMMnpip_MMnpim_woK0_wSid_n","MMnpip_MMnpim_woK0_wSid_n");
  cMMnpip_MMnpim_woK0_wSid_n->cd();
  MMnpip_MMnpim_woK0_wSid_n->Draw("colz");
  
  TCanvas *cMMnpip_MMnpim_woK0_wSid_n_px = new TCanvas("cMMnpip_MMnpim_woK0_wSid_n_px","cMMnpip_MMnpim_woK0_wSid_n_px");
  cMMnpip_MMnpim_woK0_wSid_n_px->cd();
  TH1D* MMnpip_MMnpim_woK0_wSid_n_px = (TH1D*)MMnpip_MMnpim_woK0_wSid_n->ProjectionX();
  MMnpip_MMnpim_woK0_wSid_n_px->Draw("HE");

  TCanvas *cMMnpip_MMnpim_woK0_wSid_n_py = new TCanvas("cMMnpip_MMnpim_woK0_wSid_n_py","cMMnpip_MMnpim_woK0_wSid_n_py");
  cMMnpip_MMnpim_woK0_wSid_n_py->cd();
  TH1D* MMnpip_MMnpim_woK0_wSid_n_py = (TH1D*)MMnpip_MMnpim_woK0_wSid_n->ProjectionY();
  MMnpip_MMnpim_woK0_wSid_n_py->Draw("HE");


  //TCanvas *cMMom_MMass_woK0  = new TCanvas("cMMom_MMass_woK0","MMom_MMass_woK0");
  //cMMom_MMass_woK0->cd();
  //MMom_MMass_woK0->Draw("colz");
  
  //TCanvas *cMMom_MMass_woK0_px = new TCanvas("cMMom_MMass_woK0_px","MMom_MMass_woK0_px");
  //TH1D *MMom_MMass_woK0_px = MMom_MMass_woK0->ProjectionX();
  //MMom_MMass_woK0_px->Draw("");
  
  TCanvas *cMMom_MMass_wSid  = new TCanvas("cMMom_MMass_wSid","MMom_MMass_wSid");
  cMMom_MMass_wSid->cd();
  MMom_MMass_wSid->Draw("colz");
  
  TCanvas *cMMom_MMass_woK0_wSid  = new TCanvas("cMMom_MMass_woK0_wSid","MMom_MMass_woK0_wSid");
  cMMom_MMass_woK0_wSid->cd();
  MMom_MMass_woK0_wSid->Draw("colz");
  

  TCanvas *cMMom_MMass_woK0_wSid_px  = new TCanvas("cMMom_MMass_woK0_wSid_px","MMom_MMass_woK0_wSid_px");
  cMMom_MMass_woK0_wSid_px->cd();
  //TH1D *MMom_MMass_wSid_px = MMom_MMass_wSid->ProjectionX();
  TH1D *MMom_MMass_woK0_wSid_px = MMom_MMass_woK0_wSid->ProjectionX();
  MMom_MMass_woK0_wSid_px->GetYaxis()->SetTitle("Counts/ 10 MeV/c^{2}");
  MMom_MMass_woK0_wSid_px->GetYaxis()->CenterTitle();
  //MMom_MMass_wSid_px->Draw("HE");
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

  TCanvas *cMMom_MMass_woK0_py = new TCanvas("cMMom_MMass_woK0_py","MMom_MMass_woK0_py");
  TH1D *MMom_MMass_woK0_py = MMom_MMass_woK0->ProjectionY();
  TH1D *MMom_MMass_woK0_wSid_py = MMom_MMass_woK0_wSid->ProjectionY();
  MMom_MMass_woK0_py->Draw();
  TH1D *MMom_MMass_woK0_wSid_py_clone = (TH1D*)MMom_MMass_woK0_wSid_py->Clone();
  MMom_MMass_woK0_wSid_py_clone->SetLineColor(4);
  cMMom_MMass_woK0_py->cd();
  MMom_MMass_woK0_wSid_py_clone->SetTitle("Missing Mom. d(K^{-},#pi^{+}#pi^{-}n)\"X\"");
  MMom_MMass_woK0_wSid_py_clone->Draw("");
  
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
  TFile *facc_Sp = new TFile("../simpost/acc_Sp.root","READ");
  TFile *facc_Sm = new TFile("../simpost/acc_Sm.root","READ");
  
  if(SimSpmode || SimSmmode){
    if(SimSpmode){
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
    q_IMpiSigma_gen->Draw("colz");
    TCanvas *ceff = new TCanvas("ceff","ceff");
    TH2F *h2acc=NULL;
    if(SimSpmode) h2acc =  (TH2F*)q_IMnpipi_woK0_wSid_n_Sp_acc->Clone();
    else       h2acc =  (TH2F*)q_IMnpipi_woK0_wSid_n_Sm_acc->Clone();
    h2acc->Sumw2();
    q_IMnpipi_woK0_wSid_n_Sp->Sumw2();
    q_IMnpipi_woK0_wSid_n_Sm->Sumw2();
    q_IMpiSigma_gen->Sumw2();
    q_IMpiSigma_gen->Print("base");
    q_IMnpipi_woK0_wSid_n_Sp_acc->Print("base");
    std::cout << "calc. acc." << std::endl;
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


    q_IMpiSigma_gen->Write();
    q_IMpiSigma_wSid_n_genacc->Write();
    q_IMnpipi_wSid_n_acc->Write();
    q_IMnpipi_wSid_n_acc_reco->Write();
    q_IMpiSigma_woK0_wSid_n_genacc->Write();
    q_IMnpipi_woK0_wSid_n_acc->Write();
    q_IMnpipi_woK0_wSid_n_acc_reco->Write();
    q_IMpiSigma_woK0_wSid_n_Sp_genacc->Write();
    q_IMnpipi_woK0_wSid_n_Sp_acc->Write();
    q_IMnpipi_wSid_n_Sp_acc->Write();
    q_IMnpipi_wSid_n_Sp_acc_reco->Write();
    q_IMnpipi_woK0_wSid_n_Sp_acc_reco->Write();
    q_IMpiSigma_woK0_wSid_n_Sm_genacc->Write();
    q_IMnpipi_wSid_n_Sm_acc->Write();
    q_IMnpipi_wSid_n_Sm_acc_reco->Write();
    q_IMnpipi_woK0_wSid_n_Sm_acc->Write();
    q_IMnpipi_woK0_wSid_n_Sm_acc_reco->Write();
    
    fsacc->Close();
  }//SimSpmode or SimSmmode
  
  
  //TH2F* acc_Sp = (TH2F*)facc_Sp->Get("eff_q_IMpiSigma_woK0_wSid_n_SpSm");

  //including K0, true (q,IM)
  //TH2F* acc_Sp = (TH2F*)facc_Sp->Get("eff_q_IMpiSigma_wSid_n");
  TH2F* acc_Sp = (TH2F*)facc_Sp->Get("eff_q_IMpiSigma_wSid_n_reco");
  //w/o K0, true (q,IM)
  TH2F* acc_Sp_woK0 = (TH2F*)facc_Sp->Get("eff_q_IMpiSigma_woK0_wSid_n");
  //w/o K0, Sp selection
  TH2F* acc_SpSel_woK0 = (TH2F*)facc_Sp->Get("eff_q_IMpiSigma_woK0_wSid_n_SpSm_reco");
  //TH2F* acc_SpSel_woK0 = (TH2F*)facc_Sp->Get("eff_q_IMpiSigma_woK0_wSid_n_SpSm");
  TH2F* acc_Sp_err = (TH2F*)facc_Sp->Get("heff_err");
  if(acc_Sp == NULL){
    std::cout << " acc_Sp is NULL " << std::endl;
    return;
  }
  if(acc_Sp_woK0 == NULL){
    std::cout << " acc_Sp_woK0 is NULL " << std::endl;
    return;
  }
  if(acc_SpSel_woK0 == NULL){
    std::cout << " acc_SpSel_woK0 is NULL " << std::endl;
    return;
  }
  if(acc_Sp_err == NULL){
    std::cout << " acc_Sp_err is NULL " << std::endl;
    return;
  }

  //TH2F* acc_Sm = (TH2F*)facc_Sm->Get("eff_q_IMpiSigma_woK0_wSid_n_SpSm");
  TH2F* acc_Sm = (TH2F*)facc_Sm->Get("eff_q_IMpiSigma_wSid_n");
  TH2F* acc_Sm_woK0 = (TH2F*)facc_Sm->Get("eff_q_IMpiSigma_woK0_wSid_n");
  TH2F* acc_SmSel_woK0 = (TH2F*)facc_Sm->Get("eff_q_IMpiSigma_woK0_wSid_n_SpSm_reco");
  TH2F* acc_Sm_err = (TH2F*)facc_Sm->Get("heff_err");
  if(acc_Sm == NULL){
    std::cout << " acc_Sm is NULL " << std::endl;
    return;
  }
  if(acc_Sm_woK0 == NULL){
    std::cout << " acc_Sm_woK0 is NULL " << std::endl;
    return;
  }
  if(acc_SmSel_woK0 == NULL){
    std::cout << " acc_SmSel_woK0 is NULL " << std::endl;
    return;
  }
  if(acc_Sm_err == NULL){
    std::cout << " acc_Sm is NULL " << std::endl;
    return;
  }
  f->cd(); 
  std::cout << std::endl;
  std::cout << "calculation CS of Sp mode..." << std::endl;
  std::cout << std::endl;
  const double efflumi = lumi*Beamsurvival*D2density*DAQeff*trigeff*10e-30;//micro barn(^-1)
  std::cout << "eff. lumi " << efflumi << std::endl;
  

  //Sp mode cross section
  TCanvas *cq_IMnpipi_wSid_n_cs = new TCanvas("cq_IMnpipi_wSid_n_cs","q_IMnpipi_wSid_n_cs");
  cq_IMnpipi_wSid_n_cs->cd();
  TH2F *q_IMnpipi_wSid_n_cs = (TH2F*)q_IMnpipi_wSid_n->Clone("Sp_cs");
  q_IMnpipi_wSid_n_cs->Sumw2();
  acc_Sp->Print("base");
  //cleaning up
  for(int ibinx=0;ibinx<acc_Sp_err->GetNbinsX();ibinx++){
    for(int ibiny=0;ibiny<acc_Sp_err->GetNbinsY();ibiny++){
      double precision =  acc_Sp_err->GetBinContent(ibinx,ibiny);
      if(precision>0.20){ 
        q_IMnpipi_wSid_n_cs->SetBinContent(ibinx,ibiny,0.0);
        q_IMnpipi_wSid_n_cs->SetBinError(ibinx,ibiny,0.0);
      }
    }
  }
  q_IMnpipi_wSid_n_cs->Print("base");
  q_IMnpipi_wSid_n_cs->Divide(acc_Sp);
  if(!(SimSpmode || SimSmmode))q_IMnpipi_wSid_n_cs->Scale(1./efflumi);
  q_IMnpipi_wSid_n_cs->Draw("colz");
  fkp->Draw("same");
  
  
  TCanvas *cq_IMnpipi_woK0_wSid_n_Sp_cs = new TCanvas("cq_IMnpipi_woK0_wSid_n_Sp_cs","q_IMnpipi_woK0_wSid_n_Sp_cs");
  cq_IMnpipi_woK0_wSid_n_Sp_cs->cd();
  TH2F *q_IMnpipi_woK0_wSid_n_Sp_cs = (TH2F*)q_IMnpipi_woK0_wSid_n_Sp->Clone("Sp_cs");
  q_IMnpipi_woK0_wSid_n_Sp_cs->Sumw2();
  acc_SpSel_woK0->Print("base");
  //cleaning up
  for(int ibinx=0;ibinx<acc_Sp_err->GetNbinsX();ibinx++){
    for(int ibiny=0;ibiny<acc_Sp_err->GetNbinsY();ibiny++){
      double precision =  acc_Sp_err->GetBinContent(ibinx,ibiny);
      if(precision>0.20){ 
        q_IMnpipi_woK0_wSid_n_Sp_cs->SetBinContent(ibinx,ibiny,0.0);
        q_IMnpipi_woK0_wSid_n_Sp_cs->SetBinError(ibinx,ibiny,0.0);
      }
    }
  }
  q_IMnpipi_woK0_wSid_n_Sp_cs->Print("base");
  q_IMnpipi_woK0_wSid_n_Sp_cs->Divide(acc_SpSel_woK0);
  //q_IMnpipi_woK0_wSid_n_Sp_cs->SetMaximum(30000);
  if(!(SimSpmode || SimSmmode))q_IMnpipi_woK0_wSid_n_Sp_cs->Scale(1./efflumi);
  q_IMnpipi_woK0_wSid_n_Sp_cs->Draw("colz");
  fkp->Draw("same");



  TCanvas *cq_IMnpipi_woK0_wSid_n_Sp_cs_px = new TCanvas("cq_IMnpipi_woK0_wSid_n_Sp_cs_px","q_IMnpipi_woK0_wSid_n_Sp_cs_px");
  cq_IMnpipi_woK0_wSid_n_Sp_cs_px->cd();
  TH1D *q_IMnpipi_woK0_wSid_n_Sp_cs_px = (TH1D*)q_IMnpipi_woK0_wSid_n_Sp_cs->ProjectionX();
  q_IMnpipi_woK0_wSid_n_Sp_cs_px->SetLineColor(2);
 // q_IMnpipi_woK0_wSid_n_Sp_cs_px->SetYTitle("a.u");
  q_IMnpipi_woK0_wSid_n_Sp_cs_px->SetYTitle("#mu b");
  q_IMnpipi_woK0_wSid_n_Sp_cs_px->Draw("HE");
  
  TCanvas *cq_IMnpipi_woK0_wSid_n_Sp_side_cs[ngap];
  TH2F *q_IMnpipi_woK0_wSid_n_Sp_side_cs_low[ngap];
  TH2F *q_IMnpipi_woK0_wSid_n_Sp_side_cs_high[ngap];
  TH2F *q_IMnpipi_woK0_wSid_n_Sp_side_cs[ngap];
  for(int igap=0;igap<ngap;igap++){
    cq_IMnpipi_woK0_wSid_n_Sp_side_cs[igap] = new TCanvas(Form("cq_IMnpipi_woK0_wSid_n_Sp_side_cs_%d",igap),Form("q_IMnpipi_woK0_wSid_n_Sp_side_cs_%d",igap));
    cq_IMnpipi_woK0_wSid_n_Sp_side_cs[igap]->cd();
    q_IMnpipi_woK0_wSid_n_Sp_side_cs_low[igap] = (TH2F*)q_IMnpipi_woK0_wSid_n_Sp_side[sidebandtype][LOWside][igap]->Clone(Form("Sp_cs_side_low_%d",igap));
    q_IMnpipi_woK0_wSid_n_Sp_side_cs_high[igap] = (TH2F*)q_IMnpipi_woK0_wSid_n_Sp_side[sidebandtype][HIGHside][igap]->Clone(Form("Sp_cs_side_high_%d",igap));
    q_IMnpipi_woK0_wSid_n_Sp_side_cs_low[igap]->Sumw2();
    q_IMnpipi_woK0_wSid_n_Sp_side_cs_high[igap]->Sumw2();
    q_IMnpipi_woK0_wSid_n_Sp_side_cs_low[igap]->Divide(acc_SpSel_woK0);
    q_IMnpipi_woK0_wSid_n_Sp_side_cs_high[igap]->Divide(acc_SpSel_woK0);
    if(!(SimSpmode || SimSmmode)){
      q_IMnpipi_woK0_wSid_n_Sp_side_cs_low[igap]->Scale(1./efflumi);
      q_IMnpipi_woK0_wSid_n_Sp_side_cs_high[igap]->Scale(1./efflumi);
    }
    q_IMnpipi_woK0_wSid_n_Sp_side_cs[igap] = (TH2F*)q_IMnpipi_woK0_wSid_n_Sp_side[sidebandtype][LOWside][igap]->Clone(Form("Sp_cs_side_sum_%d",igap));
    q_IMnpipi_woK0_wSid_n_Sp_side_cs[igap]->Add(q_IMnpipi_woK0_wSid_n_Sp_side[sidebandtype][HIGHside][igap]);
    q_IMnpipi_woK0_wSid_n_Sp_side_cs[igap]->Sumw2();
    q_IMnpipi_woK0_wSid_n_Sp_side_cs[igap]->Divide(acc_SpSel_woK0);
    if(!(SimSpmode || SimSmmode)){
      q_IMnpipi_woK0_wSid_n_Sp_side_cs[igap]->Scale(1./efflumi);
    }
    q_IMnpipi_woK0_wSid_n_Sp_side_cs[igap]->SetMaximum(q_IMnpipi_woK0_wSid_n_Sp_cs->GetMaximum());
    q_IMnpipi_woK0_wSid_n_Sp_side_cs[igap]->Draw("colz");
  }

  TCanvas *cq_IMnpipi_woK0_wSid_n_Sp_side_cs_px[ngap];
  TH1D* q_IMnpipi_woK0_wSid_n_Sp_side_cs_low_px[ngap];
  TH1D* q_IMnpipi_woK0_wSid_n_Sp_side_cs_high_px[ngap];
  TH1D* q_IMnpipi_woK0_wSid_n_Sp_side_cs_px[ngap];
  for(int igap=0;igap<ngap;igap++){
    cq_IMnpipi_woK0_wSid_n_Sp_side_cs_px[igap] = new TCanvas(Form("cq_IMnpipi_woK0_wSid_n_Sp_side_cs_px_%d",igap),Form("q_IMnpipi_woK0_wSid_n_Sp_side_cs_px_%d",igap));
    cq_IMnpipi_woK0_wSid_n_Sp_side_cs_px[igap]->cd();
    q_IMnpipi_woK0_wSid_n_Sp_cs_px->Draw("HE");
    q_IMnpipi_woK0_wSid_n_Sp_side_cs_low_px[igap] = (TH1D*)q_IMnpipi_woK0_wSid_n_Sp_side_cs_low[igap]->ProjectionX();
    q_IMnpipi_woK0_wSid_n_Sp_side_cs_low_px[igap]->SetLineColor(kCyan);  
    q_IMnpipi_woK0_wSid_n_Sp_side_cs_low_px[igap]->SetMarkerStyle(20);   
    q_IMnpipi_woK0_wSid_n_Sp_side_cs_low_px[igap]->SetMarkerColor(kCyan);
    q_IMnpipi_woK0_wSid_n_Sp_side_cs_low_px[igap]->Draw("Esame");
  
    q_IMnpipi_woK0_wSid_n_Sp_side_cs_high_px[igap] = (TH1D*)q_IMnpipi_woK0_wSid_n_Sp_side_cs_high[igap]->ProjectionX();
    q_IMnpipi_woK0_wSid_n_Sp_side_cs_high_px[igap]->SetLineColor(kCyan+2);  
    q_IMnpipi_woK0_wSid_n_Sp_side_cs_high_px[igap]->SetMarkerStyle(21);   
    q_IMnpipi_woK0_wSid_n_Sp_side_cs_high_px[igap]->Draw("Esame");
  
    q_IMnpipi_woK0_wSid_n_Sp_side_cs_px[igap] = (TH1D*)q_IMnpipi_woK0_wSid_n_Sp_side_cs[igap]->ProjectionX();
    q_IMnpipi_woK0_wSid_n_Sp_side_cs_px[igap]->SetLineColor(kCyan+4);      
    q_IMnpipi_woK0_wSid_n_Sp_side_cs_px[igap]->SetMarkerStyle(22);         
    q_IMnpipi_woK0_wSid_n_Sp_side_cs_px[igap]->SetMarkerColor(kCyan+4);    
    q_IMnpipi_woK0_wSid_n_Sp_side_cs_px[igap]->Draw("Esame");              
  }

  
  //Sm mode Cross section
  TCanvas *cq_IMnpipi_woK0_wSid_n_Sm_cs = new TCanvas("cq_IMnpipi_woK0_wSid_n_Sm_cs","q_IMnpipi_woK0_wSid_n_Sm_cs");
  cq_IMnpipi_woK0_wSid_n_Sm_cs->cd();
  TH2F *q_IMnpipi_woK0_wSid_n_Sm_cs = (TH2F*)q_IMnpipi_woK0_wSid_n_Sm->Clone("Sm_cs");
  q_IMnpipi_woK0_wSid_n_Sm_cs->Sumw2();
  acc_Sm->Print("base");
  //cleaning up
  for(int ibinx=0;ibinx<acc_Sm_err->GetNbinsX();ibinx++){
    for(int ibiny=0;ibiny<acc_Sm_err->GetNbinsY();ibiny++){
      double precision =  acc_Sm_err->GetBinContent(ibinx,ibiny);
      //double cont = q_IMnpipi_woK0_wSid_n_Sm_cs->GetBinContent(ibinx,ibiny);
      if(precision>0.20 ){ 
        q_IMnpipi_woK0_wSid_n_Sm_cs->SetBinContent(ibinx,ibiny,0.0);
        q_IMnpipi_woK0_wSid_n_Sm_cs->SetBinError(ibinx,ibiny,0.0);
      }
    }
  }
  q_IMnpipi_woK0_wSid_n_Sm_cs->Divide(acc_SmSel_woK0);
  if(!(SimSpmode || SimSmmode))q_IMnpipi_woK0_wSid_n_Sm_cs->Scale(1.0/efflumi);
  q_IMnpipi_woK0_wSid_n_Sm_cs->Draw("colz");
  fkp->Draw("same");
  


  TCanvas *cq_IMnpipi_woK0_wSid_n_Sm_cs_px = new TCanvas("cq_IMnpipi_woK0_wSid_n_Sm_cs_px","q_IMnpipi_woK0_wSid_n_Sm_cs_px");
  cq_IMnpipi_woK0_wSid_n_Sm_cs_px->cd();
  TH1D *q_IMnpipi_woK0_wSid_n_Sm_cs_px = (TH1D*)q_IMnpipi_woK0_wSid_n_Sm_cs->ProjectionX();
  q_IMnpipi_woK0_wSid_n_Sm_cs_px->SetLineColor(3);
  //q_IMnpipi_woK0_wSid_n_Sm_cs_px->SetYTitle("a.u");
  q_IMnpipi_woK0_wSid_n_Sm_cs_px->SetYTitle("#mu b");
  q_IMnpipi_woK0_wSid_n_Sm_cs_px->Draw("HE");


  TCanvas *cq_IMnpipi_woK0_wSid_n_SpSm_cs_px = new TCanvas("cq_IMnpipi_woK0_wSid_n_SpSm_cs_px","q_IMnpipi_woK0_wSid_n_SpSm_cs_px");
  cq_IMnpipi_woK0_wSid_n_SpSm_cs_px->cd();
  q_IMnpipi_woK0_wSid_n_Sp_cs_px->Draw("HE");
  q_IMnpipi_woK0_wSid_n_Sm_cs_px->Draw("HEsame");
  
  TCanvas *cq_IMnpipi_woK0_wSid_n_Sm_side_cs[ngap];
  TH2F *q_IMnpipi_woK0_wSid_n_Sm_side_cs_low[ngap];
  TH2F *q_IMnpipi_woK0_wSid_n_Sm_side_cs_high[ngap];
  TH2F *q_IMnpipi_woK0_wSid_n_Sm_side_cs[ngap];
  for(int igap=0;igap<ngap;igap++){
    cq_IMnpipi_woK0_wSid_n_Sm_side_cs[igap] = new TCanvas(Form("cq_IMnpipi_woK0_wSid_n_Sm_side_cs_%d",igap),Form("q_IMnpipi_woK0_wSid_n_Sm_side_cs_%d",igap));
    cq_IMnpipi_woK0_wSid_n_Sm_side_cs[igap]->cd();
    q_IMnpipi_woK0_wSid_n_Sm_side_cs_low[igap] = (TH2F*)q_IMnpipi_woK0_wSid_n_Sm_side[sidebandtype][LOWside][igap]->Clone(Form("Sm_cs_side_low_%d",igap));
    q_IMnpipi_woK0_wSid_n_Sm_side_cs_high[igap] = (TH2F*)q_IMnpipi_woK0_wSid_n_Sm_side[sidebandtype][HIGHside][igap]->Clone(Form("Sm_cs_side_high_%d",igap));
    //q_IMnpipi_woK0_wSid_n_Sm_side_cs->Add(q_IMnpipi_woK0_wSid_n_Sm_side[sidebandtype][HIGHside][0]);
    q_IMnpipi_woK0_wSid_n_Sm_side_cs_low[igap]->Sumw2();
    q_IMnpipi_woK0_wSid_n_Sm_side_cs_high[igap]->Sumw2();
    q_IMnpipi_woK0_wSid_n_Sm_side_cs_low[igap]->Divide(acc_SmSel_woK0);
    q_IMnpipi_woK0_wSid_n_Sm_side_cs_high[igap]->Divide(acc_SmSel_woK0);
    if(!(SimSpmode || SimSmmode)){
      q_IMnpipi_woK0_wSid_n_Sm_side_cs_low[igap]->Scale(1.0/efflumi);
      q_IMnpipi_woK0_wSid_n_Sm_side_cs_high[igap]->Scale(1.0/efflumi);
    }
    q_IMnpipi_woK0_wSid_n_Sm_side_cs[igap] = (TH2F*)q_IMnpipi_woK0_wSid_n_Sm_side[sidebandtype][LOWside][igap]->Clone(Form("Sm_cs_side_sum_%d",igap));
    q_IMnpipi_woK0_wSid_n_Sm_side_cs[igap]->Add(q_IMnpipi_woK0_wSid_n_Sm_side[sidebandtype][HIGHside][igap]);
    q_IMnpipi_woK0_wSid_n_Sm_side_cs[igap]->Sumw2();
    q_IMnpipi_woK0_wSid_n_Sm_side_cs[igap]->Divide(acc_SmSel_woK0);
    if(!(SimSpmode || SimSmmode)){
      q_IMnpipi_woK0_wSid_n_Sm_side_cs[igap]->Scale(1./efflumi);
    }
    q_IMnpipi_woK0_wSid_n_Sm_side_cs[igap]->SetMaximum(q_IMnpipi_woK0_wSid_n_Sm_cs->GetMaximum());
    q_IMnpipi_woK0_wSid_n_Sm_side_cs[igap]->Draw("colz");
  }

  TCanvas *cq_IMnpipi_woK0_wSid_n_Sm_side_cs_px[ngap];
  TH1D* q_IMnpipi_woK0_wSid_n_Sm_side_cs_low_px[ngap];
  TH1D* q_IMnpipi_woK0_wSid_n_Sm_side_cs_high_px[ngap];
  TH1D* q_IMnpipi_woK0_wSid_n_Sm_side_cs_px[ngap];
  for(int igap=0;igap<ngap;igap++){
    cq_IMnpipi_woK0_wSid_n_Sm_side_cs_px[igap] = new TCanvas(Form("cq_IMnpipi_woK0_wSid_n_Sm_side_cs_px_%d",igap),Form("q_IMnpipi_woK0_wSid_n_Sm_side_cs_px_%d",igap));
    cq_IMnpipi_woK0_wSid_n_Sm_side_cs_px[igap]->cd();
    q_IMnpipi_woK0_wSid_n_Sm_cs_px->Draw("HE");
    q_IMnpipi_woK0_wSid_n_Sm_side_cs_low_px[igap] = (TH1D*)q_IMnpipi_woK0_wSid_n_Sm_side_cs_low[igap]->ProjectionX();
    q_IMnpipi_woK0_wSid_n_Sm_side_cs_low_px[igap]->SetLineColor(kPink+1);
    q_IMnpipi_woK0_wSid_n_Sm_side_cs_low_px[igap]->SetMarkerStyle(20);
    q_IMnpipi_woK0_wSid_n_Sm_side_cs_low_px[igap]->Draw("Esame");
  
    q_IMnpipi_woK0_wSid_n_Sm_side_cs_high_px[igap] = (TH1D*)q_IMnpipi_woK0_wSid_n_Sm_side_cs_high[igap]->ProjectionX();
    q_IMnpipi_woK0_wSid_n_Sm_side_cs_high_px[igap]->SetLineColor(kPink+3);
    q_IMnpipi_woK0_wSid_n_Sm_side_cs_high_px[igap]->SetMarkerStyle(21);
    q_IMnpipi_woK0_wSid_n_Sm_side_cs_high_px[igap]->SetMarkerColor(kPink+3);
    q_IMnpipi_woK0_wSid_n_Sm_side_cs_high_px[igap]->Draw("Esame");
  
    q_IMnpipi_woK0_wSid_n_Sm_side_cs_px[igap] = (TH1D*)q_IMnpipi_woK0_wSid_n_Sm_side_cs[igap]->ProjectionX();
    q_IMnpipi_woK0_wSid_n_Sm_side_cs_px[igap]->SetLineColor(kPink+10);     
    q_IMnpipi_woK0_wSid_n_Sm_side_cs_px[igap]->SetMarkerStyle(22);         
    q_IMnpipi_woK0_wSid_n_Sm_side_cs_px[igap]->SetMarkerColor(kPink+10);   
    q_IMnpipi_woK0_wSid_n_Sm_side_cs_px[igap]->Draw("Esame");              
  }

  TCanvas *cIMnpim_IMnpip_dE_woK0_n_Sp_wide_cs[nzone];
  TH2D* q_IMnpipi_woK0_wSid_n_Sp_sidewide_cs[nzone];
  TH1D* q_IMnpipi_woK0_wSid_n_Sp_sidewide_cs_px[nzone];
  for(int izone=0;izone<nzone;izone++){
    cIMnpim_IMnpip_dE_woK0_n_Sp_wide_cs[izone] = new TCanvas(Form("cIMnpim_IMnpip_dE_woK0_n_Sp_wide_cs_%d",izone),Form("IMnpim_IMnpip_dE_woK0_n_Sp_wide_cs_%d",izone));   
    cIMnpim_IMnpip_dE_woK0_n_Sp_wide_cs[izone]->Divide(2,2);
    cIMnpim_IMnpip_dE_woK0_n_Sp_wide_cs[izone]->cd(1);
    IMnpim_IMnpip_dE_woK0_n_Sp->SetMaximum(IMnpim_IMnpip_dE_woK0_n->GetMaximum());
    IMnpim_IMnpip_dE_woK0_n_Sp->GetXaxis()->SetRangeUser(0,1.7);
    IMnpim_IMnpip_dE_woK0_n_Sp->GetYaxis()->SetRangeUser(0,1.7);
    IMnpim_IMnpip_dE_woK0_n_Sp->Draw("colz");
    IMnpim_IMnpip_dE_woK0_n_Sp_sidewide[izone]->Draw("colsame");
    cIMnpim_IMnpip_dE_woK0_n_Sp_wide_cs[izone]->cd(2);
    q_IMnpipi_woK0_wSid_n_Sp_cs_px->Draw("HE");
    q_IMnpipi_woK0_wSid_n_Sp_sidewide_cs[izone] = (TH2D*) q_IMnpipi_woK0_wSid_n_Sp_sidewide[izone]->Clone(Form("Sp_sidewide_cs_%d",izone));
    q_IMnpipi_woK0_wSid_n_Sp_sidewide_cs[izone]->Sumw2();
    q_IMnpipi_woK0_wSid_n_Sp_sidewide_cs[izone]->Divide(acc_SpSel_woK0);
    if(!(SimSpmode || SimSmmode)){
      q_IMnpipi_woK0_wSid_n_Sp_sidewide_cs[izone]->Scale(1./efflumi);
    }
    q_IMnpipi_woK0_wSid_n_Sp_sidewide_cs_px[izone] = (TH1D*)q_IMnpipi_woK0_wSid_n_Sp_sidewide_cs[izone]->ProjectionX();
    q_IMnpipi_woK0_wSid_n_Sp_sidewide_cs_px[izone]->SetLineColor(kCyan);
    q_IMnpipi_woK0_wSid_n_Sp_sidewide_cs_px[izone]->SetMarkerStyle(22);
    q_IMnpipi_woK0_wSid_n_Sp_sidewide_cs_px[izone]->SetMarkerColor(kCyan);
    q_IMnpipi_woK0_wSid_n_Sp_sidewide_cs_px[izone]->Draw("Esame");
    cIMnpim_IMnpip_dE_woK0_n_Sp_wide_cs[izone]->cd(3);
    q_IMnpipi_woK0_wSid_n_Sp_cs->Draw("colz");
    cIMnpim_IMnpip_dE_woK0_n_Sp_wide_cs[izone]->cd(4);
    q_IMnpipi_woK0_wSid_n_Sp_sidewide_cs[izone]->SetMaximum(q_IMnpipi_woK0_wSid_n_Sp_cs->GetMaximum());
    q_IMnpipi_woK0_wSid_n_Sp_sidewide_cs[izone]->Draw("colz");
    if(SimSpmode || SimSmmode) break;
  }

  TCanvas *cIMnpim_IMnpip_dE_woK0_n_Sm_wide_cs[nzone];
  TH2D* q_IMnpipi_woK0_wSid_n_Sm_sidewide_cs[nzone];
  TH1D* q_IMnpipi_woK0_wSid_n_Sm_sidewide_cs_px[nzone];
  for(int izone=0;izone<nzone;izone++){
    cIMnpim_IMnpip_dE_woK0_n_Sm_wide_cs[izone] = new TCanvas(Form("cIMnpim_IMnpip_dE_woK0_n_Sm_wide_cs_%d",izone),Form("IMnpim_IMnpip_dE_woK0_n_Sm_wide_cs_%d",izone));   
    cIMnpim_IMnpip_dE_woK0_n_Sm_wide_cs[izone]->Divide(2,2);
    cIMnpim_IMnpip_dE_woK0_n_Sm_wide_cs[izone]->cd(1);
    IMnpim_IMnpip_dE_woK0_n_Sm->SetMaximum(IMnpim_IMnpip_dE_woK0_n->GetMaximum());
    IMnpim_IMnpip_dE_woK0_n_Sm->GetXaxis()->SetRangeUser(0,1.7);
    IMnpim_IMnpip_dE_woK0_n_Sm->GetYaxis()->SetRangeUser(0,1.7);
    IMnpim_IMnpip_dE_woK0_n_Sm->Draw("colz");
    IMnpim_IMnpip_dE_woK0_n_Sm_sidewide[izone]->Draw("colsame");
    cIMnpim_IMnpip_dE_woK0_n_Sm_wide_cs[izone]->cd(2);
    q_IMnpipi_woK0_wSid_n_Sm_cs_px->Draw("HE");
    q_IMnpipi_woK0_wSid_n_Sm_sidewide_cs[izone] = (TH2D*) q_IMnpipi_woK0_wSid_n_Sm_sidewide[izone]->Clone(Form("Sm_sidewide_cs_%d",izone));
    q_IMnpipi_woK0_wSid_n_Sm_sidewide_cs[izone]->Sumw2();
    q_IMnpipi_woK0_wSid_n_Sm_sidewide_cs[izone]->Divide(acc_SmSel_woK0);
    if(!(SimSpmode || SimSmmode)){
      q_IMnpipi_woK0_wSid_n_Sm_sidewide_cs[izone]->Scale(1./efflumi);
    }
    q_IMnpipi_woK0_wSid_n_Sm_sidewide_cs_px[izone] = (TH1D*)q_IMnpipi_woK0_wSid_n_Sm_sidewide_cs[izone]->ProjectionX();
    q_IMnpipi_woK0_wSid_n_Sm_sidewide_cs_px[izone]->SetLineColor(kPink);
    q_IMnpipi_woK0_wSid_n_Sm_sidewide_cs_px[izone]->SetMarkerStyle(22);
    q_IMnpipi_woK0_wSid_n_Sm_sidewide_cs_px[izone]->SetMarkerColor(kPink);
    q_IMnpipi_woK0_wSid_n_Sm_sidewide_cs_px[izone]->Draw("Esame");

    cIMnpim_IMnpip_dE_woK0_n_Sm_wide_cs[izone]->cd(3);
    q_IMnpipi_woK0_wSid_n_Sm_cs->Draw("colz");

    cIMnpim_IMnpip_dE_woK0_n_Sm_wide_cs[izone]->cd(4);
    q_IMnpipi_woK0_wSid_n_Sm_sidewide_cs[izone]->SetMaximum(q_IMnpipi_woK0_wSid_n_Sm_cs->GetMaximum());
    q_IMnpipi_woK0_wSid_n_Sm_sidewide_cs[izone]->Draw("colz");
    
    if(SimSpmode || SimSmmode) break;
  }


  TCanvas *ccs_sum = new TCanvas("ccs_sum","ccs_sum");
  ccs_sum->cd();
  TH1D* cs_sum = (TH1D*)q_IMnpipi_woK0_wSid_n_Sp_cs_px->Clone();
  cs_sum->Add(q_IMnpipi_woK0_wSid_n_Sm_cs_px);
  cs_sum->SetLineColor(4);
  cs_sum->SetYTitle("#mu b");
  cs_sum->GetYaxis()->CenterTitle();
  cs_sum->Draw("HE");
  q_IMnpipi_woK0_wSid_n_Sm_cs_px->Draw("HEsame");
  q_IMnpipi_woK0_wSid_n_Sp_cs_px->Draw("HEsame");
  
  /*
  TCanvas *cq_IMnpipi_woK0_wSid_n_Sp_sub_cs = new TCanvas("cq_IMnpipi_woK0_wSid_n_Sp_sub_cs","q_IMnpipi_woK0_wSid_n_Sp_sub_cs");
  cq_IMnpipi_woK0_wSid_n_Sp_sub_cs->cd();
  TH2F *q_IMnpipi_woK0_wSid_n_Sp_sub_cs = (TH2F*)q_IMnpipi_woK0_wSid_n_Sp_sub->Clone("Sp_cs_sub");
  q_IMnpipi_woK0_wSid_n_Sp_sub_cs->Sumw2();
  q_IMnpipi_woK0_wSid_n_Sp_sub_cs->Divide(acc_Sp);
  q_IMnpipi_woK0_wSid_n_Sp_sub_cs->SetMaximum(30000);
  q_IMnpipi_woK0_wSid_n_Sp_sub_cs->Draw("colz");
  
  TCanvas *cq_IMnpipi_woK0_wSid_n_Sm_sub_cs = new TCanvas("cq_IMnpipi_woK0_wSid_n_Sm_sub_cs","q_IMnpipi_woK0_wSid_n_Sm_sub_cs");
  cq_IMnpipi_woK0_wSid_n_Sm_sub_cs->cd();
  TH2F *q_IMnpipi_woK0_wSid_n_Sm_sub_cs = (TH2F*)q_IMnpipi_woK0_wSid_n_Sm_sub->Clone("Sm_cs_sub");
  q_IMnpipi_woK0_wSid_n_Sm_sub_cs->Sumw2();
  q_IMnpipi_woK0_wSid_n_Sm_sub_cs->Divide(acc_Sm);
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
  */  
  
  TCanvas *cIMnpip_IMnpipi_wK0_n = new TCanvas("cIMnpip_IMnpipi_wK0_n","IMnpip_IMnpipi_wK0_n");
  cIMnpip_IMnpipi_wK0_n->cd();
  IMnpip_IMnpipi_wK0_n->GetYaxis()->SetRangeUser(1,1.8);
  IMnpip_IMnpipi_wK0_n->Draw("colz");
  
  TCanvas *cIMnpim_IMnpipi_wK0_n = new TCanvas("cIMnpim_IMnpipi_wK0_n","IMnpim_IMnpipi_wK0_n");
  cIMnpim_IMnpipi_wK0_n->cd();
  IMnpim_IMnpipi_wK0_n->GetYaxis()->SetRangeUser(1,1.8);
  IMnpim_IMnpipi_wK0_n->Draw("colz");

  TCanvas *cIMnpip_IMnpipi_woK0_n = new TCanvas("cIMnpip_IMnpipi_woK0_n","IMnpip_IMnpipi_woK0_n");
  cIMnpip_IMnpipi_woK0_n->cd();
  IMnpip_IMnpipi_woK0_n->GetYaxis()->SetRangeUser(1,1.8);
  IMnpip_IMnpipi_woK0_n->Draw("colz");
  
  TCanvas *cIMnpim_IMnpipi_woK0_n = new TCanvas("cIMnpim_IMnpipi_woK0_n","IMnpim_IMnpipi_woK0_n");
  cIMnpim_IMnpipi_woK0_n->cd();
  IMnpim_IMnpipi_woK0_n->GetYaxis()->SetRangeUser(1,1.8);
  IMnpim_IMnpipi_woK0_n->Draw("colz");

  TCanvas *cIMnpip_IMnpipi_n = new TCanvas("cIMnpip_IMnpipi_n","IMnpip_IMnpipi_n");
  cIMnpip_IMnpipi_n->cd();
  IMnpip_IMnpipi_n->GetYaxis()->SetRangeUser(1,1.8);
  IMnpip_IMnpipi_n->Draw("colz");
  
  TCanvas *cIMnpim_IMnpipi_n = new TCanvas("cIMnpim_IMnpipi_n","IMnpim_IMnpipi_n");
  cIMnpim_IMnpipi_n->cd();
  IMnpim_IMnpipi_n->GetYaxis()->SetRangeUser(1,1.8);
  IMnpim_IMnpipi_n->Draw("colz");

  //TCanvas *cq_IMnpipi_woK0_wSid_n_SpSm_side_cs_px = new TCanvas("cq_IMnpipi_woK0_wSid_n_SpSm_side_cs_px","q_IMnpipi_woK0_wSid_n_SpSm_side_cs_px");
  //cq_IMnpipi_woK0_wSid_n_SpSm_side_cs_px->cd();
  //TH1D* q_IMnpipi_woK0_wSid_n_SpSm_side_cs_px = q_IMnpipi_woK0_wSid_n_SpSm_side_cs->ProjectionX();
  //q_IMnpipi_woK0_wSid_n_SpSm_side_cs_px->Draw("E");
  
  //resolution evaluation
  if(SimSpmode){
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
  }//if Spmode
  
  if(SimSmmode){
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
  }//if SimSmmode

  
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
      h1->GetYaxis()->SetTitleOffset(1.4);
    }
    if(obj->InheritsFrom("TH1D")){
      h1d = (TH1D*) obj;
      h1d->GetXaxis()->CenterTitle();
      //h1d->GetXaxis()->SetTitleSize(0.05);
      //h1d->GetXaxis()->SetTitleOffset(0.80);
      h1d->GetYaxis()->SetTitleOffset(1.4);
    }
    if(obj->InheritsFrom("TH2")){
      h2 = (TH2F*) obj;
      h2->GetXaxis()->CenterTitle();
      h2->GetYaxis()->CenterTitle();
      //h2->GetXaxis()->SetTitleSize(0.05);
      //h2->GetXaxis()->SetTitleOffset(0.80);
      //h2->GetYaxis()->SetTitleSize(0.05);
      h2->GetYaxis()->SetTitleOffset(1.3);
    }
  }
  TCanvas *c = nullptr;
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
    if(SimSpmode || SimSmmode){
      pt = new TPaveText(.80,0.90,0.98,0.99,"NDC");    
    }else{
      pt = new TPaveText(.76,0.90,0.90,0.98,"NDC");    
      //pt = new TPaveText(.74,.81,0.9,0.90,"NDC");
    }
    if(SimSpmode){
      pt->SetFillColor(kAzure-4);
      pt->AddText("MC #Sigma+#pi- mode");
    }
    else if(SimSmmode){
      pt->SetFillColor(kAzure-4);
      pt->AddText("MC #Sigma-#pi+ mode"); 
    }
    else if(SimK0nn){
      pt->SetFillColor(kAzure-4);
      pt->AddText("MC K0nn"); 
    }
    else if(SimK0n_ns){
      pt->SetFillColor(kAzure-4);
      pt->AddText("MC K0n_ns"); 
    }
    else if(SimnpipiL){
      pt->SetFillColor(kAzure-4);
      pt->AddText("MC n#pi+#pi-#Lambda");
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
  std::cout << "closing pdf " << std::endl;
  
  TIter nexthist2(gDirectory->GetList());
  TString outname = std::string(filename);
  outname.Replace(std::string(filename).size()-5,5,"_out.root");
  TFile *fout = new TFile(outname.Data(),"RECREATE");
  while( (obj = (TObject*)nexthist2())!=nullptr  ){
    obj->Write();
  }
  fout->cd();
  fout->Close();
}

