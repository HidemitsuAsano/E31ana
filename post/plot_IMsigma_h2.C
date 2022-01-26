//asano memo 
//macro for analysis and histogram of p(K-,Sigma+/-)"pi-/+" 

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
const bool UseKinFitVal = false;
const unsigned int LOWside=0;
const unsigned int HIGHside=1;

//0: diagonal cut for S+/S- separation
//1: 3 sigma cut
//2: 5 simga cut
const unsigned int sigmacuttype=0;

//0:diagonal cut
//1:3 sigma cut
//2:5 sigma cut
const unsigned int sidebandtype=0;

//v198 def, 
//0 no isolation
//1 round cut
//2 revert round cut
const unsigned int IsolationFlag=1;

const bool CDCChargeVetoFlag=true;

const bool ForwardVetoFlag=false;

//check GEANT4 info. to reject fake neutron events
//maybe, also forward Sigma events should be rejected ?
const bool SimRejectFake = true;

const bool RejectStoppedSigma = true;

//color def.

void plot_IMsigma_h2(const char* filename="", const int qvalcutflag=0)
{
  gROOT->Reset();
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
  gStyle->SetCanvasDefH(800);
  gStyle->SetCanvasDefW(1000);
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadLeftMargin(0.12);

  TH1::SetDefaultSumw2();

  std::cout << "infile " << filename <<std::endl;
  TString pdfname = std::string(filename);
  if(qvalcutflag==0) pdfname.Replace(std::string(filename).size()-5,6,".pdf");
  std::cout << "pdfname: " << pdfname << std::endl;
  std::cout << std::endl;

  std::cout << "IsolationFlag ? " << IsolationFlag << std::endl;
  std::cout << "CDCChargeVetoFlag ? " << CDCChargeVetoFlag << std::endl;
  std::cout << "ForwardVetoFlag ? " << ForwardVetoFlag << std::endl;

  std::cout << std::endl;
  std::cout << "SimRejectFake ? "  << std::endl;
  if(SimRejectFake) std::cout << "Yes" << std::endl;
  else              std::cout << "No"  << std::endl;

  bool SimSpmode = (std::string(filename).find("Sp")!= std::string::npos);
  bool SimSmmode = (std::string(filename).find("Sm")!= std::string::npos);
  bool SimK0nmode = (std::string(filename).find("K0n")!= std::string::npos);//K0n
  bool RealDatamode = (std::string(filename).find("evanaIMsigma_npi_h2")!=std::string::npos);
  bool MIXmode = (std::string(filename).find("MIX")!=std::string::npos);
  
  if(RealDatamode) std::cout << "Real Data analysis" << std::endl;
  if(MIXmode) std::cout << "mixed events " << std::endl;

  TFile *f = new TFile(filename);
  //TFile *f = new TFile("sim_piSpn_dE0_Al.root");
  TTree *tree = (TTree*)f->Get("EventTree");
  if(tree==0) {
    std::cout << "EventTree is not found " << std::endl;
    std::cout << "end " << std::endl;
    return ;
  }
  //f->Print();
  tree->SetBranchAddress( "mom_beam",   &LVec_beam );
  tree->SetBranchAddress( "mom_target", &LVec_target );
  tree->SetBranchAddress( "mom_pi", &LVec_pi );
  tree->SetBranchAddress( "chargepi",&chargepi);
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
  //tree->SetBranchAddress( "CA_pi",&CA_pi);
  tree->SetBranchAddress( "CDH_Pos",&CDH_Pos);
  tree->SetBranchAddress( "CDH_Pos_pi",&CDH_Pos_pi);//
  if(SimSpmode || SimSmmode || SimK0nmode) {
    tree->SetBranchAddress( "mcncanvtxr", &mcncanvtxr);
    tree->SetBranchAddress( "mcncanvtxz", &mcncanvtxz);
    tree->SetBranchAddress( "mcncdsgen", &mcncdsgen);
    tree->SetBranchAddress( "mcpattern", &mcpattern);
    tree->SetBranchAddress( "mcmom_beam",  &mcmom_beam );
    tree->SetBranchAddress( "mcmom_pi", &mcmom_pi);//TODO :add chargepi val. to sim reco code
    tree->SetBranchAddress( "mcchargepi", &mcchargepi);//TODO :add chargepi val. to sim reco code
    tree->SetBranchAddress( "mcmom_ncds", &mcmom_ncds);
    tree->SetBranchAddress( "mcmom_pimiss", &mcmom_pimiss);
    tree->SetBranchAddress( "react_pimiss", &react_pimiss);
    tree->SetBranchAddress( "react_Sigma", &react_Sigma);
    tree->SetBranchAddress( "react_pi", &react_pi);
    tree->SetBranchAddress( "mc_vtx",&mc_vtx);
    tree->SetBranchAddress( "mc_disvtx",&mc_disvtx);
  }
  //std::cout << __LINE__ << std::endl;
  
  //std::cout << __LINE__ << std::endl;
  f->cd();
  TH1I* NHitCDCOut;
  TH1I* IsForwardCharge;
  TH2F* CDHphi_betainv_fid;
  TH2F* CDHz_betainv_fid;
  TH2F* CDHz_nmom_fid;
  TH2F* dE_betainv_fid;//
  TH2F* dE_nmom_fid_beta;
  TH2F* dE_MMom_fid_beta;
  TH2F* dE_MMass_fid_beta;
  TH2F* dE_MMass_fid_beta_wSid;
  TH2F* MMom_MMass;
  TH2F* MMom_MMass_wSid;
  //TH2F* MMom_MMass_wSid_n;
  TH1F* MMnpim_react; 

  //
  TH2F* nmom_IMnpim_mc;
  TH2F* nmom_IMnpip_mc;
  TH2F* vtxr_generation_ncan_wSid_n_mc;//
  TH2F* vtxr_diffmom_npip_ncan_wSid_n_mc;
  TH2F* vtxr_diffmom_npim_ncan_wSid_n_mc;
  TH2F* vtxr_diffMass_npip_ncan_wSid_n_mc;
  TH2F* vtxr_diffMass_npim_ncan_wSid_n_mc;
  TH2F* vtxr_vtxz_mc;
  TH2F* vtxr_vtxz_mc_pat2;
  TH2F* vtxr_vtxz_mc_pat7;
  TH2F* vtxr_vtxz_ncan_wSid_n_mc;
  TH2F* vtxr_vtxz_ncan_wSid_n_mc_pat2;
  TH2F* vtxr_vtxz_ncan_wSid_n_mc_pat7;
  TH2F* generation_diffmom_npip_ncan_wSid_n_mc;
  TH2F* generation_diffmom_npim_ncan_wSid_n_mc;
  TH2F* generation_diffMass_npip_ncan_wSid_n_mc;
  TH2F* generation_diffMass_npim_ncan_wSid_n_mc;
  
  //GEANT reaction data - mcData matching 
  TH1F* diff_IMnpim_reactmc;//store mcData - reaction data
  TH2F* IMnpim_MMnpim_mc_wSid_n_fake;//store mcData
  TH2F* IMnpim_MMnpim_mc_wSid_n_fake_pat2;//store mcData
  TH2F* IMnpim_MMnpim_mc_wSid_n_fake_pat7;//store mcData
  TH2F* IMnpim_MMnpim_mc_vtx_pat2;//store mcData 
  TH2F* IMnpim_MMnpim_mc_vtx_pat7;//store mcData
  TH2F* diff2D_IMnpim_Momnpim_wSid_n_reactmc;//store mcData - reaction data
  TH2F* diff2D_IMnpim_Momnpim_wSid_n_reactmc_fake1;//store mcData - reaction data
  TH1F* diff_IMnpip_reactmc;//store mcData - reaction data
  TH2F* IMnpip_MMnpip_mc;//store mcData 
  TH2F* IMnpip_MMnpip_mc_wSid_n_fake;//store mcData 
  TH2F* IMnpip_MMnpip_mc_wSid_n_fake_pat2;//store mcData 
  TH2F* IMnpip_MMnpip_mc_wSid_n_fake_pat7;//store mcData 
  TH2F* IMnpip_MMnpip_mc_vtx;//store mcData 
  TH2F* IMnpip_MMnpip_mc_vtx_pat2;//store mcData 
  TH2F* IMnpip_MMnpip_mc_vtx_pat7;//store mcData 
  TH2F* diff2D_IMnpip_Momnpip_wSid_n_reactmc;//store mcData - reaction data
  TH2F* diff2D_IMnpip_Momnpip_wSid_n_reactmc_fake1;//store mcData - reaction data
  TH1F* diff_nmiss_reactmc;//store mcData - reaction data mom.
  TH1F* diff_nmiss_wSid_n_reactmc;//store mcData - reaction data mom.
  TH1F* diff_cosnmiss_reactmc;//store mcData - reaction data mom.
  TH1F* diff_cosnmiss_wSid_n_reactmc;//store mcData - reaction data mom.
  TH2F* diff2D_MMnpim_IMnpim_reactmc;
  TH2F* diff2D_MMnpimissIMnpim_wSid_n_reactmc;
  TH2F* diff2D_MMnpip_IMnpip_reactmc;
  TH2F* diff2D_MMnpimissIMnpip_wSid_n_reactmc;
  
  //GEANT mcData - reco. data matching 
  //fake1: fake event flag based on mcData - reaction data matching
  TH2F* diff2D_MMnpimissIMnpim_recomc_wSid_n;//store reconstructed data - mcData
  TH2F* diff2D_MMnpimissIMnpip_recomc_wSid_n;//store reconstructed data - mcData
  TH2F* diff2D_MMnpimissIMnpim_recomc_wSid_n_fake1;//store reconstructed data - mcData
  TH2F* diff2D_MMnpimissIMnpip_recomc_wSid_n_fake1;//store reconstructed data - mcData
  TH2F* diff2D_nmom_IMnpim_recomc_wSid_n;//store reconstructed data - mcData
  TH1F* diffMMom_recomc_wSid_n;//store reconstructed data - mcData
  TH2F* diff2D_nmom_IMnpip_recomc_wSid_n;//store reconstructed data - mcData
  TH2F* diff2D_nmom_IMnpim_recomc_wSid_n_fake1;//store reconstructed data - mcData
  TH2F* diff2D_nmom_IMnpip_recomc_wSid_n_fake1;//store reconstructed data - mcData
  

  const unsigned int nwbin = 3;
  const int nbinIMnpi = 500; //1-2 GeV/c^2
  const int nbinnmiss = 150; //0-1.5 GeV/c
  const double nmisslow = 0.0;
  const double nmisshigh = 1.5;
  const int nbindE = 200;
  const int nbinnmom = 200;
  const int nbinmomnpi = 150;
  
  //bin width to be adjusted to resolution
  TH2F* IMnpim_IMnpim_mc_dE_n;//
  TH2F* IMnpim_IMnpim_mc_dE_n_vtx;//
  TH2F* IMnpim_IMnpim_mc_dE_n_vtx_pat2;//Sigma decay neutron
  TH2F* IMnpim_IMnpim_mc_dE_n_vtx_pat7;//Initial neutron
  TH2F* IMnpip_IMnpip_mc_dE_n;//
  TH2F* IMnpip_IMnpip_mc_dE_n_vtx;//
  TH2F* IMnpip_IMnpip_mc_dE_n_vtx_pat2;//
  TH2F* IMnpip_IMnpip_mc_dE_n_vtx_pat7;//
  //TH2F* pipmom_IMnpip;//no nuetron ID by CDS
  //TH2F* pipmom_IMnpip_dE;//+nuetron ID by CDS
  //TH2F* pipmom_IMnpip_dE_n;//missing nuetron ID
  //TH2F* pipmom_MMnpimissdE;//missing nuetron ID
  //TH2F* pimmom_IMnpim;//no nuetron ID
  //TH2F* pimmom_IMnpim_dE;//+nuetron ID by CDS
  //TH2F* pimmom_IMnpim_dE_n;//missing nuetron ID
  //TH2F* pimmom_MMnpimissdE;//missing nuetron ID
  TH2F* pipmom_Momnpip_wSid_n;
  TH2F* pipmom_Momnpip_wSid_n_Sp;
  TH2F* pimmom_Momnpim_wSid_n;
  TH2F* pimmom_Momnpim_wSid_n_Sm;
  TH2F* IMnpip_DCApipibeam;//DCApipibeam distance of the center point of pi+ and pi- to beam in X-Y plane
  TH2F* IMnpip_DCApipibeam_n;
  TH2F* IMnpim_DCApipibeam;//DCApipibeam distance of the center point of pi+ and pi- to beam in X-Y plane
  TH2F* IMnpim_DCApipibeam_n;
  TH2F* nmom_IMnpip_dE_n;
  TH2F* nmom_IMnpip_dE_n_Sp;
  TH2F* nmom_IMnpim_dE_n;
  TH2F* nmom_IMnpim_dE_n_Sm;
  TH2F* nmom_Momnpip;
  TH2F* nmom_Momnpip_n_Sp;
  TH2F* nmom_Momnpim;
  TH2F* nmom_Momnpim_n_Sm;
  
  //w dE cuts
  TH2F* MMnpi_IMnpip;
  TH2F* MM2npi_IMnpip;
  //TH2F* MMnpi_IMnpip_fake;
  TH2F* MMnpi_IMnpim;
  TH2F* MM2npi_IMnpim;
  //TH2F* MMnpi_IMnpim_fake;//for m
  TH2F* Cospicm_IMnpip;//cos pi cm 
  TH2F* Cospicm_IMnpim;//cos pi cm
  TH2F* MMn_IMnpip;//for checking forward K0bar
  TH2F* MMn_IMnpim;//for checking forward K0bar
  TH2F* MMpi_IMnpip;//checking forward sigma
  TH2F* MMpi_IMnpim;//checking forward sigma
  TH2F* MMn_IMnpip_pi;//for checking forward K0bar
  TH2F* MMn_IMnpim_pi;//for checking forward K0bar
  TH2F* MMpi_IMnpip_pi;//checking forward sigma
  TH2F* MMpi_IMnpim_pi;//checking forward sigma
  TH2F* Cospicm_IMnpip_pi;//cos pi cm 
  TH2F* Cospicm_IMnpim_pi;//cos pi cm
  
  TH2F* dE_CDHphi;
  TH2F* dE_CDHz;
  TH2F* dE_IMnpim;
  TH2F* dE_IMnpip;
  TH2F* dE_IMnpim_pi;
  TH2F* dE_IMnpip_pi;
  //low high side,
  //distance to signal
  TH2F* nmom_CDHphi;
  TH2F* nmom_cosn_wSid_n;//cosn = CDH neutron angle to beam axis.
  TH2F* nmom_cosnmiss_wSid_n;//

  TH2F* nmom_MMnpimisswSid;
  TH2F* nmom_MMnpimisswSid_fake;
  TH2F* nmom_MMnpimisswSid_n;
  //TH2F* diffnmom_diffdca_n;
  //TH2F* diffnmom_diffdcar_n;
  TH2F* diffnmom_diffdcaz_n;
  TH2F* vtxr_vtxz_n;
  TH2F* disvtxr_disvtxz_n;
  TH2F* mcvtxr_mcvtxz_n_mc;
  TH2F* mcdisvtxr_mcdisvtxz_n_mc;
  TH2F* diffnmom_diffdisvtx_CApip_n;
  TH2F* diffnmom_diffdisvtx_CApip_r_n;
  TH2F* diffnmom_diffdisvtx_CApip_z_n;
  TH2F* diffnmom_diffdisvtx_CApim_n;
  TH2F* diffnmom_diffdisvtx_CApim_r_n;
  TH2F* diffnmom_diffdisvtx_CApim_z_n;
  TH2F* diffnmom_diffdisvtx_cdcbeam_pip_n;
  TH2F* diffnmom_diffdisvtx_cdcbeam_pip_r_n;
  TH2F* diffnmom_diffdisvtx_cdcbeam_pip_z_n;
  TH2F* diffnmom_diffdisvtx_cdcbeam_pim_n;
  TH2F* diffnmom_diffdisvtx_cdcbeam_pim_r_n;
  TH2F* diffnmom_diffdisvtx_cdcbeam_pim_z_n;

  TH2F* diff2d_CDC_CDH_pim;//neutral hit - pim track projected position
  TH2F* diff2d_CDC_CDH_pip;//neutral hit - pip track projected position
  TH2F* diff2d_CDC_CDH_pim_phi_tof;//neutral hit - pim track projected position (phi) vs tof(ncan - pim)
  TH2F* diff2d_CDC_CDH_pip_phi_tof;//neutral hit - pip track projected position (phi) vs tof(ncan - pip)
  TH2F* diff2d_CDC_CDH_pim_z_tof;//neutral hit - pim track projected positio (z) vs tof (ncan - -pim)
  TH2F* diff2d_CDC_CDH_pip_z_tof;//neutral hit - pip track projected position (z) vs tof (ncan -pip)
  TH2F* dE_diffphi_CDC_CDH_pim;
  TH2F* dE_diffphi_CDC_CDH_pip;
  //TH2F* MMnpimissdE;
  TH2F* MMnpimissdiffphi_CDC_CDH_pim;
  TH2F* MMnpimissdiffphi_CDC_CDH_pip;
  TH2F* pimmom_diffphi_CDC_CDH_pim;
  TH2F* pipmom_diffphi_CDC_CDH_pip;
  TH2F* pimmom_diffz_CDC_CDH_pim;
  TH2F* pipmom_diffz_CDC_CDH_pip;
  TH2F* nmom_diffphi_CDC_CDH_pim;
  TH2F* nmom_diffphi_CDC_CDH_pip;
  TH2F* nmom_diffz_CDC_CDH_pim;
  TH2F* nmom_diffz_CDC_CDH_pip;

  std::cout << __LINE__ << std::endl;

  
  NHitCDCOut = new TH1I("NHitCDCOut","NHitCDCOut",10,0,10);
  NHitCDCOut->SetXTitle("Nhits of CDC Outer 2 layers");
  
  IsForwardCharge = new TH1I("IsForwardCharge","IsForwardCharge",2,0,2);

  CDHphi_betainv_fid = new TH2F("CDHphi_betainv_fid","CDHphi_betainv_fid",100,0,10,100,-TMath::Pi(),TMath::Pi());
  CDHphi_betainv_fid->SetXTitle("1/#beta");
  CDHphi_betainv_fid->SetYTitle("CDH phi");

  CDHz_betainv_fid = new TH2F("CDHz_betainv_fid","CDHz_betainv_fid",100,0,10,100,-50,50);
  CDHz_betainv_fid->SetXTitle("1/#beta");
  CDHz_betainv_fid->SetYTitle("CDH z [cm]");

  CDHz_nmom_fid = new TH2F("CDHz_nmom_fid","CDHz_nmom_fid",nbinnmom,0,1,100,-50,50);
  CDHz_nmom_fid->SetXTitle("nmom. [GeV/c]");
  CDHz_nmom_fid->SetYTitle("CDH z [cm]");

  dE_betainv_fid = new TH2F("dE_betainv_fid","dE_betainv_fid",1000, 0, 50, nbindE, 0, 50);
  dE_betainv_fid->SetXTitle("1/#beta");
  dE_betainv_fid->SetYTitle("dE [MeVee]");

  dE_nmom_fid_beta = new TH2F("dE_nmom_fid_beta","dE_nmom_fid_beta",1000, 0, 1.5, 500, 0, 50);
  dE_nmom_fid_beta->SetXTitle("Mom. [GeV/c]");
  dE_nmom_fid_beta->SetYTitle("dE [MeVee]");

  //dE_nmom_fid_beta_wK0 = new TH2F("dE_nmom_fid_beta_wK0","dE_nmom_fid_beta_wK0",100, 0, 1.5, 500, 0, 50);
  //dE_nmom_fid_beta_wK0->SetXTitle("Mom. [GeV/c]");
  //dE_nmom_fid_beta_wK0->SetYTitle("dE [MeVee]");

  dE_MMom_fid_beta = new TH2F(Form("dE_MMom_fid_beta"),Form("dE_MMom_fid_beta"),100, 0, 1.5, nbindE, 0, 50);
  dE_MMom_fid_beta->SetXTitle("Missing Mom. [GeV/c]");
  dE_MMom_fid_beta->SetYTitle("dE [MeVee]");

  dE_MMass_fid_beta = new TH2F(Form("dE_MMass_fid_beta"),Form("dE_MMass_fid_beta"), nbinnmiss, nmisslow, nmisshigh, nbindE, 0, 50);
  dE_MMass_fid_beta->SetXTitle("Missing mass [GeV/c^{2}]");
  dE_MMass_fid_beta->SetYTitle("dE [MeVee]");

  dE_MMass_fid_beta_wSid = new TH2F(Form("dE_MMass_fid_beta_wSid"),Form("dE_MMass_fid_beta_wSid"), nbinnmiss, nmisslow, nmisshigh, nbindE, 0, 10);
  dE_MMass_fid_beta_wSid->SetXTitle("Missing mass [GeV/c^{2}]");
  dE_MMass_fid_beta_wSid->SetYTitle("dE [MeVee]");

  MMom_MMass = new TH2F("MMom_MMass","MMom_MMass", nbinnmiss, nmisslow, nmisshigh, 200, 0, 2.0);
  MMom_MMass->SetXTitle("Missing Mass [GeV/c^{2}]");
  MMom_MMass->SetYTitle("Missing Mom. [GeV/c]");

  MMom_MMass_wSid = new TH2F(Form("MMom_MMass_wSid"),Form("MMom_MMass_wSid"), nbinnmiss, nmisslow, nmisshigh, 200, 0, 2.0);
  MMom_MMass_wSid->SetXTitle("Missing Mass [GeV/c^{2}]");
  MMom_MMass_wSid->SetYTitle("Missing Mom. [GeV/c]");
  
  //MMom_MMass_wSid_n = new TH2F(Form("MMom_MMass_wSid_n"),Form("MMom_MMass_wSid_n"), nbinnmiss, nmisslow, nmisshigh, 200, 0, 2.0);
  //MMom_MMass_wSid_n->SetXTitle("Missing Mass [GeV/c^{2}]");
  //MMom_MMass_wSid_n->SetYTitle("Missing Mom. [GeV/c]");

  //MMnpim_react = new TH1F("MMnpimissreact","MMnpimissreact",100,0,1.5);
  //MMnpim_react->SetXTitle("react. Missing Mass [GeV/c^{2}");
  
  diff_IMnpim_reactmc = new TH1F("diff_IMnpim_reactmc","diff_IMnpim_reactmc",100,-1.0,1.0);
  diff_IMnpim_reactmc->SetXTitle("diff. IMnpim (MCData - Reac.) [GeV/c^{2}]");
  
  
  IMnpim_MMnpim_mc_wSid_n_fake = new TH2F("IMnpim_MMnpim_mc_wSid_n_fake","IMnpim_MMnpim_mc_wSid_n_fake",140,1.,1.7,140,1.0,1.7);
  IMnpim_MMnpim_mc_wSid_n_fake->SetYTitle("true  IM(n_{CDS}#pi^{-}) [GeV/c^{2}]");
  IMnpim_MMnpim_mc_wSid_n_fake->SetXTitle("true  IM(n_{miss}#pi^{-}) [GeV/c^{2}]");
  
  IMnpim_MMnpim_mc_wSid_n_fake_pat2 = new TH2F("IMnpim_MMnpim_mc_wSid_n_fake_pat2","IMnpim_MMnpim_mc_wSid_n_fake_pat2",140,1.,1.7,140,1.0,1.7);
  IMnpim_MMnpim_mc_wSid_n_fake_pat2->SetYTitle("true  IM(n_{CDS}#pi^{-}) [GeV/c^{2}]");
  IMnpim_MMnpim_mc_wSid_n_fake_pat2->SetXTitle("true  IM(n_{miss}#pi^{-}) [GeV/c^{2}]");
  
  IMnpim_MMnpim_mc_wSid_n_fake_pat7 = new TH2F("IMnpim_MMnpim_mc_wSid_n_fake_pat7","IMnpim_MMnpim_mc_wSid_n_fake_pat7",140,1.,1.7,140,1.0,1.7);
  IMnpim_MMnpim_mc_wSid_n_fake_pat7->SetYTitle("true  IM(n_{CDS}#pi^{-}) [GeV/c^{2}]");
  IMnpim_MMnpim_mc_wSid_n_fake_pat7->SetXTitle("true  IM(n_{miss}#pi^{-}) [GeV/c^{2}]");
  
  
  IMnpim_MMnpim_mc_vtx_pat2 = new TH2F("IMnpim_MMnpim_mc_vtx_pat2","IMnpim_MMnpim_mc_vtx_pat2",140,1.,1.7,140,1.0,1.7);
  IMnpim_MMnpim_mc_vtx_pat2->SetYTitle("true  IM(n_{CDS}#pi^{-})  [GeV/c^{2}]");
  IMnpim_MMnpim_mc_vtx_pat2->SetXTitle("true  IM(n_{miss}#pi^{-}) [GeV/c^{2}]");

  IMnpim_MMnpim_mc_vtx_pat7 = new TH2F("IMnpim_MMnpim_mc_vtx_pat7","IMnpim_MMnpim_mc_vtx_pat7",140,1.,1.7,140,1.0,1.7);
  IMnpim_MMnpim_mc_vtx_pat7->SetYTitle("true  IM(n_{CDS}#pi^{-})  [GeV/c^{2}]");
  IMnpim_MMnpim_mc_vtx_pat7->SetXTitle("true  IM(n_{miss}#pi^{-}) [GeV/c^{2}]");


  diff2D_IMnpim_Momnpim_wSid_n_reactmc = new TH2F("diff2D_IMnpim_Momnpim_wSid_n_reactmc","diff2D_IMnpim_Momnpim_wSid_n_reactmc",100,-1.0,1.0,100,-1.0,1.0);
  diff2D_IMnpim_Momnpim_wSid_n_reactmc->SetXTitle("diff. Mom.(n#pi^{-}) (MCData - Reac.) [GeV/c]");
  diff2D_IMnpim_Momnpim_wSid_n_reactmc->SetYTitle("diff. IM(n#pi^{-}) (MCData - Reac.) [GeV/c^{2}]");
  
  diff2D_IMnpim_Momnpim_wSid_n_reactmc_fake1 = new TH2F("diff2D_IMnpim_Momnpim_wSid_n_reactmc_fake1","diff2D_IMnpim_Momnpim_wSid_n_reactmc_fake1",100,-1.0,1.0,100,-1.0,1.0);
  diff2D_IMnpim_Momnpim_wSid_n_reactmc_fake1->SetXTitle("diff. Mom.(n#pi^{-}) (MCData - Reac.) [GeV/c]");
  diff2D_IMnpim_Momnpim_wSid_n_reactmc_fake1->SetYTitle("diff. IM(n#pi^{-}) (MCData - Reac.) [GeV/c^{2}]");
  
  diff_IMnpip_reactmc = new TH1F("diff_IMnpip_reactmc","diff_IMnpip_reactmc",1000,-1.0,1.0);
  diff_IMnpip_reactmc->SetXTitle("diff. IMnpip (MCData - Reac.) [GeV/c^{2}]");
  
  IMnpip_MMnpip_mc = new TH2F("IMnpip_MMnpip_mc","IMnpip_MMnpip_mc",140,1,1.7,140,1.0,1.7);
  IMnpip_MMnpip_mc->SetYTitle("true IM(n_{CDS}#pi^{+}) [GeV/c^{2}]");
  IMnpip_MMnpip_mc->SetXTitle("true IM(n_{miss}#pi^{+}) [GeV/c^{2}]");
  
  IMnpip_MMnpip_mc_wSid_n_fake = new TH2F("IMnpip_MMnpip_mc_wSid_n_fake","IMnpip_MMnpip_mc_wSid_n_fake",140,1,1.7,140,1.0,1.7);
  IMnpip_MMnpip_mc_wSid_n_fake->SetYTitle("true IM(n_{CDS}#pi^{+}) [GeV/c^{2}]");
  IMnpip_MMnpip_mc_wSid_n_fake->SetXTitle("true IM(n_{miss}#pi^{+}) [GeV/c^{2}]");
  
  IMnpip_MMnpip_mc_wSid_n_fake_pat2 = new TH2F("IMnpip_MMnpip_mc_wSid_n_fake_pat2","IMnpip_MMnpip_mc_wSid_n_fake_pat2",140,1,1.7,140,1.0,1.7);
  IMnpip_MMnpip_mc_wSid_n_fake_pat2->SetYTitle("true IM(n_{CDS}#pi^{+}) [GeV/c^{2}]");
  IMnpip_MMnpip_mc_wSid_n_fake_pat2->SetXTitle("true IM(n_{miss}#pi^{+}) [GeV/c^{2}]");
  
  IMnpip_MMnpip_mc_wSid_n_fake_pat7 = new TH2F("IMnpip_MMnpip_mc_wSid_n_fake_pat7","IMnpip_MMnpip_mc_wSid_n_fake_pat7",140,1,1.7,140,1.0,1.7);
  IMnpip_MMnpip_mc_wSid_n_fake_pat7->SetYTitle("true IM(n_{CDS}#pi^{+}) [GeV/c^{2}]");
  IMnpip_MMnpip_mc_wSid_n_fake_pat7->SetXTitle("true IM(n_{miss}#pi^{+}) [GeV/c^{2}]");
  
  
  IMnpip_MMnpip_mc_vtx = new TH2F("IMnpip_MMnpip_mc_vtx","IMnpip_MMnpip_mc_vtx",140,1,1.7,140,1.,1.7);
  IMnpip_MMnpip_mc_vtx->SetYTitle("true IM(n_{CDS}#pi^{+}) [GeV/c^{2}]");
  IMnpip_MMnpip_mc_vtx->SetXTitle("true IM(n_{miss}#pi^{+} ) [GeV/c^{2}]");
  
  IMnpip_MMnpip_mc_vtx_pat2 = new TH2F("IMnpip_MMnpip_mc_vtx_pat2","IMnpip_MMnpip_mc_vtx_pat2",140,1,1.7,140,1.,1.7);
  IMnpip_MMnpip_mc_vtx_pat2->SetYTitle("true IM(n_{CDS}#pi^{+}) [GeV/c^{2}]");
  IMnpip_MMnpip_mc_vtx_pat2->SetXTitle("true IM(n_{miss}#pi^{+} ) [GeV/c^{2}]");
  
  IMnpip_MMnpip_mc_vtx_pat7 = new TH2F("IMnpip_MMnpip_mc_vtx_pat7","IMnpip_MMnpip_mc_vtx_pat7",140,1,1.7,140,1.,1.7);
  IMnpip_MMnpip_mc_vtx_pat7->SetYTitle("true IM(n_{CDS}#pi^{+}) [GeV/c^{2}]");
  IMnpip_MMnpip_mc_vtx_pat7->SetXTitle("true IM(n_{miss}#pi^{+} ) [GeV/c^{2}]");
  
  
  
  diff2D_IMnpip_Momnpip_wSid_n_reactmc = new TH2F("diff2D_IMnpip_Momnpip_wSid_n_reactmc","diff2D_IMnpip_Momnpip_wSid_n_reactmc",100,-1.0,1.0,100,-1.0,1.0);
  diff2D_IMnpip_Momnpip_wSid_n_reactmc->SetXTitle("diff. Mom.(n#pi^{+}) (MCData - Reac.) [GeV/c]");
  diff2D_IMnpip_Momnpip_wSid_n_reactmc->SetYTitle("diff. IM(n#pi^{+}) (MCData - Reac.) [GeV/c^{2}]");

  diff2D_IMnpip_Momnpip_wSid_n_reactmc_fake1 = new TH2F("diff2D_IMnpip_Momnpip_wSid_n_reactmc_fake1","diff2D_IMnpip_Momnpip_wSid_n_reactmc_fake1",100,-1.0,1.0,100,-1.0,1.0);
  diff2D_IMnpip_Momnpip_wSid_n_reactmc_fake1->SetXTitle("diff. Mom.(n#pi^{+}) (MCData - Reac.) [GeV/c]");
  diff2D_IMnpip_Momnpip_wSid_n_reactmc_fake1->SetYTitle("diff. IM(n#pi^{+}) (MCData - Reac.) [GeV/c^{2}]");
  
  diff_nmiss_reactmc = new TH1F("diff_nmiss_reactmc","diff_nmiss_reactmc",1500,-1.5,1.5);
  diff_nmiss_reactmc->SetXTitle("diff. nmiss (MCData - Reac.) [GeV/c]");
  
  diff_nmiss_wSid_n_reactmc = new TH1F("diff_nmiss_wSid_n_reactmc","diff_nmiss_wSid_n_reactmc",1500,-1.5,1.5);
  diff_nmiss_wSid_n_reactmc->SetXTitle("diff. nmiss (MCData - Reac.) [GeV/c]");
  
  diff_cosnmiss_reactmc = new TH1F("diff_cosnmiss_reactmc","diff_cosnmiss_reactmc",2000,-1.0,1.0);
  diff_cosnmiss_reactmc->SetXTitle("diff. cos. nmiss (MCData - Reac.) [radian]");
  
  diff_cosnmiss_wSid_n_reactmc = new TH1F("diff_cosnmiss_wSid_n_reactmc","diff_cosnmiss_wSid_n_reactmc",200,-1.0,1.0);
  diff_cosnmiss_wSid_n_reactmc->SetXTitle("diff. cos. nmiss (MCData - Reac.) [radian]");
  
  diff2D_MMnpim_IMnpim_reactmc = new TH2F("diff2D_MMnpim_IMnpim_reactmc","diff2D_MMnpim_IMnpim_reactmc",200,-1.0,1.0,150,-1.5,1.5);
  diff2D_MMnpim_IMnpim_reactmc->SetXTitle("diff. IMnpim (MCData - Reac.) [GeV/c^{2}]");
  diff2D_MMnpim_IMnpim_reactmc->SetYTitle("diff. Missing Mass (MCData - Reac.) [GeV/c^{2}]");
   
  diff2D_MMnpimissIMnpim_wSid_n_reactmc = new TH2F("diff2D_MMnpimissIMnpim_wSid_n_reactmc","diff2D_MMnpimissIMnpim_wSid_n_reactmc",200,-1.0,1.0,150,-1.5,1.5);
  diff2D_MMnpimissIMnpim_wSid_n_reactmc->SetXTitle("diff. IMnpim (MCData - Reac.) [GeV/c^{2}]");
  diff2D_MMnpimissIMnpim_wSid_n_reactmc->SetYTitle("diff. Missing Mass (MCData - Reac.) [GeV/c^{2}]");

  diff2D_MMnpip_IMnpip_reactmc = new TH2F("diff2D_MMnpip_IMnpip_reactmc","diff2D_MMnpip_IMnpip_reactmc",200,-1.0,1.0,150,-1.5,1.5);
  diff2D_MMnpip_IMnpip_reactmc->SetXTitle("diff. IMnpip (MCData - Reac.) [GeV/c^{2}]");
  diff2D_MMnpip_IMnpip_reactmc->SetYTitle("diff. Missing Mass (MCData - Reac.) [GeV/c^{2}]");
   
  diff2D_MMnpimissIMnpip_wSid_n_reactmc = new TH2F("diff2D_MMnpimissIMnpip_wSid_n_reactmc","diff2D_MMnpimissIMnpip_wSid_n_reactmc",200,-1.0,1.0,150,-1.5,1.5);
  diff2D_MMnpimissIMnpip_wSid_n_reactmc->SetXTitle("diff. IMnpip (MCData - Reac.) [GeV/c^{2}]");
  diff2D_MMnpimissIMnpip_wSid_n_reactmc->SetYTitle("diff. Missing Mass (MCData - Reac.) [GeV/c^{2}]");

  diff2D_MMnpimissIMnpim_recomc_wSid_n = new TH2F("diff2D_MMnpimissIMnpim_recomc_wSid_n","diff2D_MMnpimissIMnpim_recomc_wSid_n",200,-1.0,1.0,150,-1.5,1.5);
  diff2D_MMnpimissIMnpim_recomc_wSid_n->SetXTitle("diff. IMnpim (reco. - MCData) [GeV/c^{2}]");
  diff2D_MMnpimissIMnpim_recomc_wSid_n->SetYTitle("diff. Missing Mass (reco. - MCData) [GeV/c^{2}]");
  
  diff2D_MMnpimissIMnpim_recomc_wSid_n_fake1 = new TH2F("diff2D_MMnpimissIMnpim_recomc_wSid_n_fake1","diff2D_MMnpimissIMnpim_recomc_wSid_n_fake1",200,-1.0,1.0,150,-1.5,1.5);
  diff2D_MMnpimissIMnpim_recomc_wSid_n_fake1->SetXTitle("diff. IMnpim (reco. - MCData) [GeV/c^{2}]");
  diff2D_MMnpimissIMnpim_recomc_wSid_n_fake1->SetYTitle("diff. Missing Mass (reco. - MCData) [GeV/c^{2}]");

  diff2D_MMnpimissIMnpip_recomc_wSid_n = new TH2F("diff2D_MMnpimissIMnpip_recomc_wSid_n","diff2D_MMnpimissIMnpip_recomc_wSid_n",200,-1.0,1.0,150,-1.5,1.5);
  diff2D_MMnpimissIMnpip_recomc_wSid_n->SetXTitle("diff. IMnpip (reco. - MCData) [GeV/c^{2}]");
  diff2D_MMnpimissIMnpip_recomc_wSid_n->SetYTitle("diff. Missing Mass (reco. - MCData) [GeV/c^{2}]");
  
  diff2D_MMnpimissIMnpip_recomc_wSid_n_fake1 = new TH2F("diff2D_MMnpimissIMnpip_recomc_wSid_n_fake1","diff2D_MMnpimissIMnpip_recomc_wSid_n_fake1",200,-1.0,1.0,150,-1.5,1.5);
  diff2D_MMnpimissIMnpip_recomc_wSid_n_fake1->SetXTitle("diff. IMnpip (reco. - MCData) [GeV/c^{2}]");
  diff2D_MMnpimissIMnpip_recomc_wSid_n_fake1->SetYTitle("diff. Missing Mass (reco. - MCData) [GeV/c^{2}]");

  diff2D_nmom_IMnpim_recomc_wSid_n = new TH2F("diff2D_nmom_IMnpim_recomc_wSid_n","diff2D_nmom_IMnpim_recomc_wSid_n",2000,-0.2,0.2,400,-0.4,0.4);
  diff2D_nmom_IMnpim_recomc_wSid_n->SetXTitle("diff. IMnpim (reco. - MCData) [GeV/c^{2}]");
  diff2D_nmom_IMnpim_recomc_wSid_n->SetYTitle("diff. n_{CDS} mom. (reco. - MCData) [GeV/c]");
  
  diff2D_nmom_IMnpim_recomc_wSid_n_fake1 = new TH2F("diff2D_nmom_IMnpim_recomc_wSid_n_fake1","diff2D_nmom_IMnpim_recomc_wSid_n_fake1",200,-1.0,1.0,150,-1.5,1.5);
  diff2D_nmom_IMnpim_recomc_wSid_n_fake1->SetXTitle("diff. IMnpim (reco. - MCData) [GeV/c^{2}]");
  diff2D_nmom_IMnpim_recomc_wSid_n_fake1->SetYTitle("diff. n_{CDS} mom. (reco. - MCData) [GeV/c]");

  if(SimSpmode || SimSmmode || SimK0nmode){
    diff2D_nmom_IMnpip_recomc_wSid_n = new TH2F("diff2D_nmom_IMnpip_recomc_wSid_n","diff2D_nmom_IMnpip_recomc_wSid_n",2000,-0.2,0.2,400,-0.4,0.4);
    diff2D_nmom_IMnpip_recomc_wSid_n->SetXTitle("diff. IMnpip (reco. - MCData) [GeV/c^{2}]");
    diff2D_nmom_IMnpip_recomc_wSid_n->SetYTitle("diff. n_{CDS} mom. (reco. - MCData) [GeV/c]");

    diff2D_nmom_IMnpip_recomc_wSid_n_fake1 = new TH2F("diff2D_nmom_IMnpip_recomc_wSid_n_fake1","diff2D_nmom_IMnpip_recomc_wSid_n_fake1",200,-1.0,1.0,150,-1.5,1.5);
    diff2D_nmom_IMnpip_recomc_wSid_n_fake1->SetXTitle("diff. IMnpip (reco. - MCData) [GeV/c^{2}]");
    diff2D_nmom_IMnpip_recomc_wSid_n_fake1->SetYTitle("diff. n_{CDS} mom. (reco. - MCData) [GeV/c]");

    diffMMom_recomc_wSid_n = new TH1F("diffMMom_recomc_wSid_n","diffMMom_recomc_wSid_n",1000,-1.0,1.0);
  
    nmom_IMnpim_mc = new TH2F("nmom_IMnpim_mc", "nmom_IMnpim_mc",nbinIMnpi, 1, 2.0, nbinnmom,0,1.0);
    nmom_IMnpim_mc->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
    nmom_IMnpim_mc->SetYTitle("n_{CDS} mom. [GeV/c]");

    nmom_IMnpip_mc = new TH2F("nmom_IMnpip_mc", "nmom_IMnpip_mc",nbinIMnpi, 1, 2.0, nbinnmom,0,1.0);
    nmom_IMnpip_mc->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
    nmom_IMnpip_mc->SetYTitle("n_{CDS} mom. [GeV/c]");

    vtxr_generation_ncan_wSid_n_mc = new TH2F("vtxr_generation_ncan_wSid_n_mc","vtxr_generation_ncan_wSid_n_mc",10,0,10, 480,0,120.0);
    vtxr_generation_ncan_wSid_n_mc->SetXTitle("n_{CDS} generation");
    vtxr_generation_ncan_wSid_n_mc->SetYTitle("n_{CDS} origin in R [cm]");

    vtxr_diffmom_npip_ncan_wSid_n_mc = new TH2F("vtxr_diffmom_npip_ncan_wSid_n_mc","vtxr_diffmom_npip_ncan_wSid_n_mc",200,-0.1,0.1,480,0,120.0);
    vtxr_diffmom_npip_ncan_wSid_n_mc->SetXTitle("diff. mom. (n#pi^{+}) MCdata - React.");
    vtxr_diffmom_npip_ncan_wSid_n_mc->SetYTitle("n_{CDS} origin in R [cm] ");

    vtxr_diffMass_npip_ncan_wSid_n_mc = new TH2F("vtxr_diffMass_npip_ncan_wSid_n_mc","vtxr_diffMass_npip_ncan_wSid_n_mc",200,-0.1,0.1,480,0,120.0);
    vtxr_diffMass_npip_ncan_wSid_n_mc->SetXTitle("diff. Mass. (n#pi^{+}) MCdata - React.");
    vtxr_diffMass_npip_ncan_wSid_n_mc->SetYTitle("n_{CDS} origin in R [cm] ");

    vtxr_diffmom_npim_ncan_wSid_n_mc = new TH2F("vtxr_diffmom_npim_ncan_wSid_n_mc","vtxr_diffmom_npim_ncan_wSid_n_mc",200,-0.1,0.1,480,0,120.0);
    vtxr_diffmom_npim_ncan_wSid_n_mc->SetXTitle("diff. mom. (n#pi^{-}) MCdata - React.");
    vtxr_diffmom_npim_ncan_wSid_n_mc->SetYTitle("n_{CDS} origin in R [cm] ");

    vtxr_diffMass_npim_ncan_wSid_n_mc = new TH2F("vtxr_diffMass_npim_ncan_wSid_n_mc","vtxr_diffMass_npim_ncan_wSid_n_mc",200,-0.1,0.1,480,0,120.0);
    vtxr_diffMass_npim_ncan_wSid_n_mc->SetXTitle("diff. Mass. (n#pi^{-}) MCdata - React.");
    vtxr_diffMass_npim_ncan_wSid_n_mc->SetYTitle("n_{CDS} origin in R [cm] ");

    vtxr_vtxz_mc = new TH2F("vtxr_vtxz_mc","vtxr_vtxz_mc",75,-75,75,480,0,120.0);
    vtxr_vtxz_mc->SetXTitle("n_{CDS} origin in Z [cm] ");
    vtxr_vtxz_mc->SetYTitle("n_{CDS} origin in R [cm] ");

    vtxr_vtxz_mc_pat2 = new TH2F("vtxr_vtxz_mc_pat2","vtxr_vtxz_mc_pat2",75,-75,75,480,0,120.0);
    vtxr_vtxz_mc_pat2->SetXTitle("n_{CDS} origin in Z [cm] ");
    vtxr_vtxz_mc_pat2->SetYTitle("n_{CDS} origin in R [cm] ");

    vtxr_vtxz_mc_pat7 = new TH2F("vtxr_vtxz_mc_pat7","vtxr_vtxz_mc_pat7",75,-75,75,480,0,120.0);
    vtxr_vtxz_mc_pat7->SetXTitle("n_{CDS} origin in Z [cm] ");
    vtxr_vtxz_mc_pat7->SetYTitle("n_{CDS} origin in R [cm] ");

    vtxr_vtxz_ncan_wSid_n_mc = new TH2F("vtxr_vtxz_ncan_wSid_n_mc","vtxr_vtxz_ncan_wSid_n_mc",75,-75,75,480,0,120.0);
    vtxr_vtxz_ncan_wSid_n_mc->SetXTitle("n_{CDS} origin in Z [cm] ");
    vtxr_vtxz_ncan_wSid_n_mc->SetYTitle("n_{CDS} origin in R [cm] ");

    vtxr_vtxz_ncan_wSid_n_mc_pat2 = new TH2F("vtxr_vtxz_ncan_wSid_n_mc_pat2","vtxr_vtxz_ncan_wSid_n_mc_pat2",75,-75,75,480,0,120.0);
    vtxr_vtxz_ncan_wSid_n_mc_pat2->SetXTitle("n_{CDS} origin in Z [cm] ");
    vtxr_vtxz_ncan_wSid_n_mc_pat2->SetYTitle("n_{CDS} origin in R [cm] ");

    vtxr_vtxz_ncan_wSid_n_mc_pat7 = new TH2F("vtxr_vtxz_ncan_wSid_n_mc_pat7","vtxr_vtxz_ncan_wSid_n_mc_pat7",75,-75,75,480,0,120.0);
    vtxr_vtxz_ncan_wSid_n_mc_pat7->SetXTitle("n_{CDS} origin in Z [cm] ");
    vtxr_vtxz_ncan_wSid_n_mc_pat7->SetYTitle("n_{CDS} origin in R [cm] ");

    generation_diffmom_npip_ncan_wSid_n_mc = new TH2F("generation_diffmom_npip_ncan_wSid_n_mc","generation_diffmom_npip_ncan_wSid_n_mc",200,-0.1,0.1,10,0,10);
    generation_diffmom_npip_ncan_wSid_n_mc->SetXTitle("diff. mom. (n#pi^{+}) MCdata - React.");
    generation_diffmom_npip_ncan_wSid_n_mc->SetYTitle("n_{CDS} generation");

    generation_diffMass_npip_ncan_wSid_n_mc = new TH2F("generation_diffMass_npip_ncan_wSid_n_mc","generation_diffMass_npip_ncan_wSid_n_mc",200,-0.1,0.1,10,0,10);
    generation_diffMass_npip_ncan_wSid_n_mc->SetXTitle("diff. Mass. (n#pi^{+}) MCdata - React.");
    generation_diffMass_npip_ncan_wSid_n_mc->SetYTitle("n_{CDS} generation");

    generation_diffmom_npim_ncan_wSid_n_mc = new TH2F("generation_diffmom_npim_ncan_wSid_n_mc","generation_diffmom_npim_ncan_wSid_n_mc",200,-0.1,0.1,10,0,10);
    generation_diffmom_npim_ncan_wSid_n_mc->SetXTitle("diff. mom. (n#pi^{-}) MCdata - React.");
    generation_diffmom_npim_ncan_wSid_n_mc->SetYTitle("n_{CDS} generation");

    generation_diffMass_npim_ncan_wSid_n_mc = new TH2F("generation_diffMass_npim_ncan_wSid_n_mc","generation_diffMass_npim_ncan_wSid_n_mc",200,-0.1,0.1,10,0,10);
    generation_diffMass_npim_ncan_wSid_n_mc->SetXTitle("diff. Mass. (n#pi^{-}) MCdata - React.");
    generation_diffMass_npim_ncan_wSid_n_mc->SetYTitle("n_{CDS} generation");
  }

  IMnpim_IMnpim_mc_dE_n = new TH2F("IMnpim_IMnpim_mc_dE_n","IMnpim_IMnpim_mc_dE_n",nbinIMnpi, 1, 2.0, nbinIMnpi, -1.0, 1.0);
  IMnpim_IMnpim_mc_dE_n->SetXTitle("reco. IM(n#pi^{-}) [GeV/c^{2}]");
  IMnpim_IMnpim_mc_dE_n->SetYTitle("reco. - true IM(n#pi^{-}) [GeV/c^{2}]");
  
  IMnpim_IMnpim_mc_dE_n_vtx = new TH2F("IMnpim_IMnpim_mc_dE_n_vtx","IMnpim_IMnpim_mc_dE_n_vtx",nbinIMnpi, 1, 2.0, nbinIMnpi, -1.0, 1.0);
  IMnpim_IMnpim_mc_dE_n_vtx->SetXTitle("reco. IM(n#pi^{-}) [GeV/c^{2}]");
  IMnpim_IMnpim_mc_dE_n_vtx->SetYTitle("reco. true IM(n#pi^{-}) [GeV/c^{2}]");
  
  IMnpim_IMnpim_mc_dE_n_vtx_pat2 = new TH2F("IMnpim_IMnpim_mc_dE_n_vtx_pat2","IMnpim_IMnpim_mc_dE_n_vtx_pat2",nbinIMnpi, 1, 2.0, nbinIMnpi, -1.0, 1.0);
  IMnpim_IMnpim_mc_dE_n_vtx_pat2->SetXTitle("reco. IM(n#pi^{-}) [GeV/c^{2}]");
  IMnpim_IMnpim_mc_dE_n_vtx_pat2->SetYTitle("reco. - true IM(n#pi^{-}) [GeV/c^{2}]");
  
  IMnpim_IMnpim_mc_dE_n_vtx_pat7 = new TH2F("IMnpim_IMnpim_mc_dE_n_vtx_pat7","IMnpim_IMnpim_mc_dE_n_vtx_pat7",nbinIMnpi, 1, 2.0, nbinIMnpi, -1.0, 1.0);
  IMnpim_IMnpim_mc_dE_n_vtx_pat7->SetXTitle("reco. IM(n#pi^{-}) [GeV/c^{2}]");
  IMnpim_IMnpim_mc_dE_n_vtx_pat7->SetYTitle("reco. - true IM(n#pi^{-}) [GeV/c^{2}]");
  
  IMnpip_IMnpip_mc_dE_n = new TH2F("IMnpip_IMnpip_mc_dE_n","IMnpip_IMnpip_mc_dE_n",nbinIMnpi, 1, 2.0, nbinIMnpi, -1.0, 1.0);
  IMnpip_IMnpip_mc_dE_n->SetXTitle("reco. IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpip_IMnpip_mc_dE_n->SetYTitle("reco. - true IM(n#pi^{+}) [GeV/c^{2}]");
  
  IMnpip_IMnpip_mc_dE_n_vtx = new TH2F("IMnpip_IMnpip_mc_dE_n_vtx","IMnpip_IMnpip_mc_dE_n_vtx",nbinIMnpi, 1, 2.0, nbinIMnpi, -1.0, 1.0);
  IMnpip_IMnpip_mc_dE_n_vtx->SetXTitle("reco. IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpip_IMnpip_mc_dE_n_vtx->SetYTitle("reco. - true IM(n#pi^{+}) [GeV/c^{2}]");
  
  IMnpip_IMnpip_mc_dE_n_vtx_pat2 = new TH2F("IMnpip_IMnpip_mc_dE_n_vtx_pat2","IMnpip_IMnpip_mc_dE_n_vtx_pat2",nbinIMnpi, 1, 2.0, nbinIMnpi, -1.0, 1.0);
  IMnpip_IMnpip_mc_dE_n_vtx_pat2->SetXTitle("reco. IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpip_IMnpip_mc_dE_n_vtx_pat2->SetYTitle("reco. - true IM(n#pi^{+}) [GeV/c^{2}]");
  
  IMnpip_IMnpip_mc_dE_n_vtx_pat7 = new TH2F("IMnpip_IMnpip_mc_dE_n_vtx_pat7","IMnpip_IMnpip_mc_dE_n_vtx_pat7",nbinIMnpi, 1, 2.0, nbinIMnpi, -1.0, 1.0);
  IMnpip_IMnpip_mc_dE_n_vtx_pat7->SetXTitle("reco. IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpip_IMnpip_mc_dE_n_vtx_pat7->SetYTitle("reco. - true IM(n#pi^{+}) [GeV/c^{2}]");
  
  //pimmom_IMnpim = new TH2F("pimmom_IMnpim","pimmom_IMnpim",nbinIMnpi,1.,2.,200,0,1);
  //pimmom_IMnpim->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}] ");
  //pimmom_IMnpim->SetYTitle("Mom(#pi^{-}) [GeV/c] ");

  //pimmom_IMnpim_dE = new TH2F("pimmom_IMnpim_dE","pimmom_IMnpim_dE",nbinIMnpi,1.,2.,200,0,1);
  //pimmom_IMnpim_dE->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}] ");
  //pimmom_IMnpim_dE->SetYTitle("Mom(#pi^{-}) [GeV/c] ");

  //pimmom_IMnpim_dE_n = new TH2F("pimmom_IMnpim_dE_n","pimmom_IMnpim_dE_n",nbinIMnpi,1.,2.,200,0,1);
  //pimmom_IMnpim_dE_n->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}] ");
  //pimmom_IMnpim_dE_n->SetYTitle("Mom(#pi^{-}) [GeV/c] ");

  //pimmom_MMnpimissdE = new TH2F("pimmom_MMnpimissdE","pimmom_MMnpimissdE",nbinnmiss, nmisslow, nmisshigh,200,0,1);
  //pimmom_MMnpimissdE->SetXTitle("Miss Mass. [GeV/c^{2}]");
  //pimmom_MMnpimissdE->SetYTitle("Mom(#pi^{-}) [GeV/c] ");

  pipmom_Momnpip_wSid_n = new TH2F("pipmom_Momnpip_wSid_n","pipmom_Momnpip_wSid_n",200,0.,1.0,200,0.,1.0);
  pipmom_Momnpip_wSid_n->SetXTitle("Mom(n#pi^{+}) [GeV/c]");
  pipmom_Momnpip_wSid_n->SetYTitle("Mom(#pi^{+}) [GeV/c]");
  
  pipmom_Momnpip_wSid_n_Sp = new TH2F("pipmom_Momnpip_wSid_n_Sp","pipmom_Momnpip_wSid_n_Sp",200,0.,1.0,200,0.,1.0);
  pipmom_Momnpip_wSid_n_Sp->SetXTitle("Mom(n#pi^{+}) [GeV/c]");
  pipmom_Momnpip_wSid_n_Sp->SetYTitle("Mom(#pi^{+}) [GeV/c]");
  
  pimmom_Momnpim_wSid_n = new TH2F("pimmom_Momnpim_wSid_n","pimmom_Momnpim_wSid_n",200,0,1,200,0.,1.0);
  pimmom_Momnpim_wSid_n->SetXTitle("Mom(n#pi^{-}) [GeV/c]");
  pimmom_Momnpim_wSid_n->SetYTitle("Mom(#pi^{-}) [GeV/c]");
  
  pimmom_Momnpim_wSid_n_Sm = new TH2F("pimmom_Momnpim_wSid_n_Sm","pimmom_Momnpim_wSid_n_Sm",200,0.,1.0,200,0.,1.0);
  pimmom_Momnpim_wSid_n_Sm->SetXTitle("Mom(n#pi^{-}) [GeV/c]");
  pimmom_Momnpim_wSid_n_Sm->SetYTitle("Mom(#pi^{-}) [GeV/c]");

  IMnpip_DCApipibeam = new TH2F("IMnpip_DCApipibeam","IMnpip_DCApipibeam",300,0,20,nbinIMnpi,1.,2.);
  IMnpip_DCApipibeam->SetXTitle("DCA #pi#pi-beam [cm]");
  IMnpip_DCApipibeam->SetYTitle("IM(n#pi^{+}) [GeV/c^{2}] ");

  IMnpip_DCApipibeam_n = new TH2F("IMnpip_DCApipibeam_n","IMnpip_DCApipibeam_n",300,0,20,nbinIMnpi,1.,2.);
  IMnpip_DCApipibeam_n->SetXTitle("DCA #pi#pi-beam [cm]");
  IMnpip_DCApipibeam_n->SetYTitle("IM(n#pi^{+}) [GeV/c^{2}] ");

  IMnpim_DCApipibeam = new TH2F("IMnpim_DCApipibeam","IMnpim_DCApipibeam",300,0,20,nbinIMnpi,1.,2.);
  IMnpim_DCApipibeam->SetXTitle("DCA #pi#pi-beam [cm]");
  IMnpim_DCApipibeam->SetYTitle("IM(n#pi^{+}) [GeV/c^{2}] ");

  IMnpim_DCApipibeam_n = new TH2F("IMnpim_DCApipibeam_n","IMnpim_DCApipibeam_n",300,0,20,nbinIMnpi,1.,2.);
  IMnpim_DCApipibeam_n->SetXTitle("DCA #pi#pi-beam [cm]");
  IMnpim_DCApipibeam_n->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}] ");

  std::cout << __LINE__ << std::endl;
  
  nmom_IMnpim_dE_n = new TH2F("nmom_IMnpim_dE_n","nmom_IMnpim_dE_n",nbinIMnpi,1.,2.,nbinnmom,0.,1.0);
  nmom_IMnpim_dE_n->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  nmom_IMnpim_dE_n->SetYTitle("nmom. [GeV/c]");

  nmom_IMnpim_dE_n_Sm = new TH2F("nmom_IMnpim_dE_n_Sm","nmom_IMnpim_dE_n_Sm",nbinIMnpi,1.,2.,nbinnmom,0.,1.0);
  nmom_IMnpim_dE_n_Sm->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  nmom_IMnpim_dE_n_Sm->SetYTitle("nmom. [GeV/c]");
  
  nmom_Momnpim = new TH2F("nmom_Momnpim","nmom_Momnpim",100,0.,1.,nbinnmom,0.,1.0);
  nmom_Momnpim->SetXTitle("Mom.(n#pi^{-}) [GeV/c]");
  nmom_Momnpim->SetYTitle("nmom. [GeV/c]");
  
  nmom_Momnpim_n_Sm = new TH2F("nmom_Momnpim_n_Sm","nmom_Momnpim_n_Sm",100,0.,1.,nbinnmom,0.,1.0);
  nmom_Momnpim_n_Sm->SetXTitle("Mom.(n#pi^{-}) [GeV/c]");
  nmom_Momnpim_n_Sm->SetYTitle("nmom. [GeV/c]");

  nmom_IMnpip_dE_n = new TH2F("nmom_IMnpip_dE_n","nmom_IMnpip_dE_n",nbinIMnpi,1.,2.,nbinnmom,0.,1.0);
  nmom_IMnpip_dE_n->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  nmom_IMnpip_dE_n->SetYTitle("nmom. [GeV/c]");

  nmom_IMnpip_dE_n_Sp = new TH2F("nmom_IMnpip_dE_n_Sp","nmom_IMnpip_dE_n_Sp",nbinIMnpi,1.,2.,nbinnmom,0.,1.0);
  nmom_IMnpip_dE_n_Sp->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  nmom_IMnpip_dE_n_Sp->SetYTitle("nmom. [GeV/c]");

  nmom_Momnpip = new TH2F("nmom_Momnpip","nmom_Momnpip",100,0.,1.,nbinnmom,0.,1.0);
  nmom_Momnpip->SetXTitle("Mom.(n#pi^{+}) [GeV/c]");
  nmom_Momnpip->SetYTitle("nmom. [GeV/c]");

  nmom_Momnpip_n_Sp = new TH2F("nmom_Momnpip_n_Sp","nmom_Momnpip_n_Sp",100,0.,1.,nbinnmom,0.,1.0);
  nmom_Momnpip_n_Sp->SetXTitle("Mom.(n#pi^{+}) [GeV/c]");
  nmom_Momnpip_n_Sp->SetYTitle("nmom. [GeV/c]");
  
  MMnpi_IMnpip = new TH2F("MMnpi_IMnpip","MMnpi_IMnpip",700,1.,1.7,200,-1.,1.);
  MMnpi_IMnpip->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  MMnpi_IMnpip->SetYTitle("Miss. Mass (n#pi^{+}) [GeV/c^{2}]");
  
  MM2npi_IMnpip = new TH2F("MM2npi_IMnpip","MM2npi_IMnpip",700,1.,1.7,200,-1,1.);
  MM2npi_IMnpip->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  MM2npi_IMnpip->SetYTitle("Miss. Mass^{2} (n#pi^{+}) [(GeV/c^{2})^{2}]");
  
  MMnpi_IMnpim = new TH2F("MMnpi_IMnpim","MMnpi_IMnpim",700,1.,1.7,200,-1.,1.);
  MMnpi_IMnpim->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  MMnpi_IMnpim->SetYTitle("Miss. Mass (n#pi^{-}) [GeV/c^{2}]");  
  
  MM2npi_IMnpim = new TH2F("MM2npi_IMnpim","MM2npi_IMnpim",700,1.,1.7,200,-1.,1.);
  MM2npi_IMnpim->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  MM2npi_IMnpim->SetYTitle("Miss. Mass^{2} (n#pi^{-}) [(GeV/c^{2})^{2}]");  

  Cospicm_IMnpip = new TH2F("Cospicm_IMnpip","Cospicm_IMnpip",700,1,1.7 ,50,0,1);
  Cospicm_IMnpip->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  Cospicm_IMnpip->SetYTitle("Cos.#theta_{CM} miss-#pi" );
  
  Cospicm_IMnpim = new TH2F("Cospicm_IMnpim","Cospicm_IMnpim",700,1,1.7 ,50,0,1);
  Cospicm_IMnpim->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  Cospicm_IMnpim->SetYTitle("Cos.#theta_{CM} miss-#pi" );
  
  Cospicm_IMnpip_pi = new TH2F("Cospicm_IMnpip_pi","Cospicm_IMnpip_pi",700,1,1.7 ,50,0,1);
  Cospicm_IMnpip_pi->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  Cospicm_IMnpip_pi->SetYTitle("Cos.#theta_{CM} miss-#pi" );
  
  Cospicm_IMnpim_pi = new TH2F("Cospicm_IMnpim_pi","Cospicm_IMnpim_pi",700,1,1.7 ,50,0,1);
  Cospicm_IMnpim_pi->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  Cospicm_IMnpim_pi->SetYTitle("Cos.#theta_{CM} miss-#pi" );
  
  MMn_IMnpip = new TH2F("MMn_IMnpip","MMn_IMnpip",700,1.,1.7,100,0.,1.);
  MMn_IMnpip->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  MMn_IMnpip->SetYTitle("Miss. Mass (n) [GeV/c^{2}]");  
  
  MMn_IMnpim = new TH2F("MMn_IMnpim","MMn_IMnpim",700,1.,1.7,100,0.,1.);
  MMn_IMnpim->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  MMn_IMnpim->SetYTitle("Miss. Mass (n) [GeV/c^{2}]");  

  MMpi_IMnpip = new TH2F("MMpi_IMnpip","MMpi_IMnpip",700,1.,1.7,100,1,2.);
  MMpi_IMnpip->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  MMpi_IMnpip->SetYTitle("Miss. Mass (#pi^{+}) [GeV/c^{2}]");  
  
  MMpi_IMnpim = new TH2F("MMpi_IMnpim","MMpi_IMnpim",700,1.,1.7,100,1.,2.);
  MMpi_IMnpim->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  MMpi_IMnpim->SetYTitle("Miss. Mass (#pi^{-}) [GeV/c^{2}]");  

  MMn_IMnpip_pi = new TH2F("MMn_IMnpip_pi","MMn_IMnpip_pi",700,1.,1.7,100,0.,1.);
  MMn_IMnpip_pi->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  MMn_IMnpip_pi->SetYTitle("Miss. Mass (n) [GeV/c^{2}]");  
  
  MMn_IMnpim_pi = new TH2F("MMn_IMnpim_pi","MMn_IMnpim_pi",700,1.,1.7,100,0.,1.);
  MMn_IMnpim_pi->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  MMn_IMnpim_pi->SetYTitle("Miss. Mass (n) [GeV/c^{2}]");  

  MMpi_IMnpip_pi = new TH2F("MMpi_IMnpip_pi","MMpi_IMnpip_pi",700,1.,1.7,100,1,2.);
  MMpi_IMnpip_pi->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  MMpi_IMnpip_pi->SetYTitle("Miss. Mass (#pi^{+}) [GeV/c^{2}]");  
  
  MMpi_IMnpim_pi = new TH2F("MMpi_IMnpim_pi","MMpi_IMnpim_pi",700,1.,1.7,100,1.,2.);
  MMpi_IMnpim_pi->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  MMpi_IMnpim_pi->SetYTitle("Miss. Mass (#pi^{-}) [GeV/c^{2}]");  
  
  dE_CDHphi = new TH2F(Form("dE_CDHphi"),Form("dE_CDHphi"),100,-3.14,3.14, nbindE,0,50);
  dE_CDHphi->SetXTitle("CDH phi");
  dE_CDHphi->SetYTitle("dE [MeVee]");

  dE_CDHz = new TH2F(Form("dE_CDHz"),"dE_CDHz",100,-50,50,nbindE,0,50);
  dE_CDHz->SetXTitle("CDH z [cm]");
  dE_CDHz->SetYTitle("dE [MeVee]");

  dE_IMnpim = new TH2F(Form("dE_IMnpim"),Form("dE_IMnpim"), nbinIMnpi, 1.0, 2.0, nbindE, 0, 50.);
  dE_IMnpim->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  dE_IMnpim->SetYTitle("dE [MeVee]");

  dE_IMnpim_pi = new TH2F(Form("dE_IMnpim_pi"),Form("dE_IMnpim_pi"), nbinIMnpi, 1.0, 2.0, nbindE, 0, 50.);
  dE_IMnpim_pi->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  dE_IMnpim_pi->SetYTitle("dE [MeVee]");

  dE_IMnpip = new TH2F(Form("dE_IMnpip"),Form("dE_IMnpip"), nbinIMnpi, 1.0, 2.0, nbindE, 0, 50.);
  dE_IMnpip->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  dE_IMnpip->SetYTitle("dE [MeVee]");

  dE_IMnpip_pi = new TH2F(Form("dE_IMnpip_pi"),Form("dE_IMnpip_pi"), nbinIMnpi, 1.0, 2.0, nbindE, 0, 50.);
  dE_IMnpip_pi->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  dE_IMnpip_pi->SetYTitle("dE [MeVee]");

  nmom_CDHphi = new TH2F("nmom_CDHphi","nmom_CDHphi",100,-3.14,3.14,nbinnmom,0,1.0);
  nmom_CDHphi->SetXTitle("CDH phi");
  nmom_CDHphi->SetYTitle("nmom [GeV/c]");

  nmom_cosn_wSid_n = new TH2F("nmom_cosn_wSid_n","nmom_cosn_wSid_n",100,-1.0,1.0,nbinnmom,0,1.0);
  nmom_cosn_wSid_n->SetXTitle("nCDS cos#theta_{LAB}");
  nmom_cosn_wSid_n->SetYTitle("nCDS mom [GeV/c]");

  std::cout << __LINE__ << std::endl;

  nmom_cosnmiss_wSid_n = new TH2F("nmom_cosnmiss_wSid_n","nmom_cosnmiss_wSid_n",100,-1.0,1.0,nbinnmom,0,1.0);

  nmom_MMnpimisswSid = new TH2F("nmom_MMnpimisswSid","nmom_MMnpimisswSid", nbinnmiss, nmisslow, nmisshigh,nbinnmom,0,1.0);
  nmom_MMnpimisswSid->SetXTitle("Miss. Mass [GeV/c^{2}]");
  nmom_MMnpimisswSid->SetYTitle("nmom  [GeV/c]");
  
  nmom_MMnpimisswSid_fake = new TH2F("nmom_MMnpimisswSid_fake","nmom_MMnpimisswSid_fake", nbinnmiss, nmisslow, nmisshigh,nbinnmom,0,1.0);
  nmom_MMnpimisswSid_fake->SetXTitle("Miss. Mass [GeV/c^{2}]");
  nmom_MMnpimisswSid_fake->SetYTitle("nmom  [GeV/c]");
  
  nmom_MMnpimisswSid_n = new TH2F("nmom_MMnpimisswSid_n","nmom_MMnpimisswSid_n", nbinnmiss, nmisslow, nmisshigh,nbinnmom,0,1.0);
  nmom_MMnpimisswSid_n->SetXTitle("Miss. Mass [GeV/c^{2}]");
  nmom_MMnpimisswSid_n->SetYTitle("nmom  [GeV/c]");
  
  //diffnmom_diffdca_n = new TH2F("diffnmom_diffdca_n","diffnmom_diffdca_n",1000,-2,2,1000,-0.2,0.2);
  //diffnmom_diffdca_n->SetXTitle("reco - true n_{CDS} vertex. [cm]");
  //diffnmom_diffdca_n->SetYTitle("reco - true n_{CDS} mom. [GeV/c]");

  //diffnmom_diffdcar_n = new TH2F("diffnmom_diffdcar_n","diffnmom_diffdcar_n",1000,-2,2,1000,-0.2,0.2);
  //diffnmom_diffdcar_n->SetXTitle("reco - true n_{CDS} vertex in r. [cm]");
  //diffnmom_diffdcar_n->SetYTitle("reco - true n_{CDS} mom. [GeV/c]");
  
  diffnmom_diffdcaz_n = new TH2F("diffnmom_diffdcaz_n","diffnmom_diffdcaz_n",1000,-4,4,1000,-0.2,0.2);
  diffnmom_diffdcaz_n->SetXTitle("reco - true n_{CDS} vertex in z. [cm]");
  diffnmom_diffdcaz_n->SetYTitle("reco - true n_{CDS} mom. [GeV/c]");
  
  vtxr_vtxz_n = new TH2F("vtxr_vtxz_n","vtxr_vtxz_n",1000,-12.5,12.5,300,0,30);
  vtxr_vtxz_n->SetYTitle("reco. vtx X [cm]  ");
  vtxr_vtxz_n->SetXTitle("reco. vtx Z [cm]  ");

  disvtxr_disvtxz_n = new TH2F("disvtxr_disvtxz_n","disvtxr_disvtxz_n",1000,-12.5,12.5,300,0,30);
  disvtxr_disvtxz_n->SetYTitle("reco. vtx X [cm]  ");
  disvtxr_disvtxz_n->SetXTitle("reco. vtx Z [cm]  ");
  
  mcvtxr_mcvtxz_n_mc = new TH2F("mcvtxr_mcvtxz_n_mc","mcvtxr_mcvtxz_n_mc",1000,-12.5,12.5,300,0,30);
  mcvtxr_mcvtxz_n_mc->SetYTitle("mc vtx X [cm]  ");
  mcvtxr_mcvtxz_n_mc->SetXTitle("mc vtx Z [cm]  ");
  
  mcdisvtxr_mcdisvtxz_n_mc = new TH2F("mcdisvtxr_mcdisvtxz_n_mc","mcdisvtxr_mcdisvtxz_n_mc",1000,-12.5,12.5,300,0,30);
  mcdisvtxr_mcdisvtxz_n_mc->SetYTitle("mc vtx X [cm]  ");
  mcdisvtxr_mcdisvtxz_n_mc->SetXTitle("mc vtx Z [cm]  ");

  diffnmom_diffdisvtx_CApip_n = new TH2F("diffnmom_diffdisvtx_CApip_n","diffnmom_diffdisvtx_CApip_n",100,0,4,100,-0.2,0.2);
  diffnmom_diffdisvtx_CApip_n->SetXTitle("diff vtx [cm] ");
  diffnmom_diffdisvtx_CApip_n->SetYTitle("reco. - true nmom [GeV/c]");

  diffnmom_diffdisvtx_CApip_r_n = new TH2F("diffnmom_diffdisvtx_CApip_r_n","diffnmom_diffdisvtx_CApip_r_n",100,0,4,100,-0.2,0.2);
  diffnmom_diffdisvtx_CApip_r_n->SetXTitle("diff vtx_r [cm] ");
  diffnmom_diffdisvtx_CApip_r_n->SetYTitle("reco. - true nmom [GeV/c]");

  diffnmom_diffdisvtx_CApip_z_n = new TH2F("diffnmom_diffdisvtx_CApip_z_n","diffnmom_diffdisvtx_CApip_z_n",100,-4,4,100,-0.2,0.2);
  diffnmom_diffdisvtx_CApip_z_n->SetXTitle("diff vtx_z [cm] ");
  diffnmom_diffdisvtx_CApip_z_n->SetYTitle("reco. - true nmom [GeV/c]");

  diffnmom_diffdisvtx_CApim_n = new TH2F("diffnmom_diffdisvtx_CApim_n","diffnmom_diffdisvtx_CApim_n",100,0,4,100,-0.2,0.2);
  diffnmom_diffdisvtx_CApim_n->SetXTitle("diff vtx [cm] ");
  diffnmom_diffdisvtx_CApim_n->SetYTitle("reco. - true nmom [GeV/c]");

  diffnmom_diffdisvtx_CApim_r_n = new TH2F("diffnmom_diffdisvtx_CApim_r_n","diffnmom_diffdisvtx_CApim_r_n",100,0,4,100,-0.2,0.2);
  diffnmom_diffdisvtx_CApim_r_n->SetXTitle("diff vtx_r [cm] ");
  diffnmom_diffdisvtx_CApim_r_n->SetYTitle("reco. - true nmom [GeV/c]");

  diffnmom_diffdisvtx_CApim_z_n = new TH2F("diffnmom_diffdisvtx_CApim_z_n","diffnmom_diffdisvtx_CApim_z_n",100,-4,4,100,-0.2,0.2);
  diffnmom_diffdisvtx_CApim_z_n->SetXTitle("diff vtx_z [cm] ");
  diffnmom_diffdisvtx_CApim_z_n->SetYTitle("reco. - true nmom [GeV/c]");
    
  diffnmom_diffdisvtx_cdcbeam_pip_n = new TH2F("diffnmom_diffdisvtx_cdcbeam_pip_n","diffnmom_diffdisvtx_cdcbeam_pip_n",100,0,4,100,-0.2,0.2);
  diffnmom_diffdisvtx_cdcbeam_pip_n->SetXTitle("diff vtx [cm] ");
  diffnmom_diffdisvtx_cdcbeam_pip_n->SetYTitle("reco. - true nmom [GeV/c]");

  diffnmom_diffdisvtx_cdcbeam_pip_r_n = new TH2F("diffnmom_diffdisvtx_cdcbeam_pip_r_n","diffnmom_diffdisvtx_cdcbeam_pip_r_n",100,0,4,100,-0.2,0.2);
  diffnmom_diffdisvtx_cdcbeam_pip_r_n->SetXTitle("diff vtx_r [cm] ");
  diffnmom_diffdisvtx_cdcbeam_pip_r_n->SetYTitle("reco. - true nmom [GeV/c]");

  diffnmom_diffdisvtx_cdcbeam_pip_z_n = new TH2F("diffnmom_diffdisvtx_cdcbeam_pip_z_n","diffnmom_diffdisvtx_cdcbeam_pip_z_n",100,-4,4,100,-0.2,0.2);
  diffnmom_diffdisvtx_cdcbeam_pip_z_n->SetXTitle("diff vtx_z [cm] ");
  diffnmom_diffdisvtx_cdcbeam_pip_z_n->SetYTitle("reco. - true nmom [GeV/c]");
  
  diffnmom_diffdisvtx_cdcbeam_pim_n = new TH2F("diffnmom_diffdisvtx_cdcbeam_pim_n","diffnmom_diffdisvtx_cdcbeam_pim_n",100,0,4,100,-0.2,0.2);
  diffnmom_diffdisvtx_cdcbeam_pim_n->SetXTitle("diff vtx [cm] ");
  diffnmom_diffdisvtx_cdcbeam_pim_n->SetYTitle("reco. - true nmom [GeV/c]");

  diffnmom_diffdisvtx_cdcbeam_pim_r_n = new TH2F("diffnmom_diffdisvtx_cdcbeam_pim_r_n","diffnmom_diffdisvtx_cdcbeam_pim_r_n",100,0,4,100,-0.2,0.2);
  diffnmom_diffdisvtx_cdcbeam_pim_r_n->SetXTitle("diff vtx_r [cm] ");
  diffnmom_diffdisvtx_cdcbeam_pim_r_n->SetYTitle("reco. - true nmom [GeV/c]");

  diffnmom_diffdisvtx_cdcbeam_pim_z_n = new TH2F("diffnmom_diffdisvtx_cdcbeam_pim_z_n","diffnmom_diffdisvtx_cdcbeam_pim_z_n",100,-4,4,100,-0.2,0.2);
  diffnmom_diffdisvtx_cdcbeam_pim_z_n->SetXTitle("diff vtx_z [cm] ");
  diffnmom_diffdisvtx_cdcbeam_pim_z_n->SetYTitle("reco. - true nmom [GeV/c]");

  diff2d_CDC_CDH_pim = new TH2F("diff2d_CDC_CDH_pim","diff2d_CDC_CDH_pim",100,-1.*TMath::Pi(),TMath::Pi(),100,-100,100);
  diff2d_CDC_CDH_pim->SetXTitle("CDH hit - #pi^{-} track (phi) [radian]");
  diff2d_CDC_CDH_pim->SetYTitle("CDH hit - #pi^{-} track (z) [cm]");
  
  diff2d_CDC_CDH_pim_phi_tof = new TH2F("diff2d_CDC_CDH_pim_phi_tof","diff2d_CDC_CDH_pim_phi_tof",100,-1.*TMath::Pi(),TMath::Pi(),120,-10,50);
  diff2d_CDC_CDH_pim_phi_tof->SetXTitle("CDH hit - #pi^{-} track (phi) [radian]");
  diff2d_CDC_CDH_pim_phi_tof->SetYTitle("CDH hit - #pi^{-} TOF [nsec]");
  
  diff2d_CDC_CDH_pim_z_tof = new TH2F("diff2d_CDC_CDH_pim_z_tof","diff2d_CDC_CDH_pim_z_tof",100,-100,100,120,-10,50);
  diff2d_CDC_CDH_pim_z_tof->SetXTitle("CDH hit - #pi^{-} track (z) [cm]");
  diff2d_CDC_CDH_pim_z_tof->SetYTitle("CDH hit - #pi^{-} TOF [nsec]");

  diff2d_CDC_CDH_pip = new TH2F("diff2d_CDC_CDH_pip","diff2d_CDC_CDH_pip",100,-1.*TMath::Pi(),TMath::Pi(),100,-100,100);
  diff2d_CDC_CDH_pip->SetXTitle("CDH hit - #pi^{+} track (phi) [radian]");
  diff2d_CDC_CDH_pip->SetYTitle("CDH hit - #pi^{+} track (z) [cm]");
  
  diff2d_CDC_CDH_pip_phi_tof = new TH2F("diff2d_CDC_CDH_pip_phi_tof","diff2d_CDC_CDH_pip_phi_tof",100,-1.*TMath::Pi(),TMath::Pi(),120,-10,50);
  diff2d_CDC_CDH_pip_phi_tof->SetXTitle("CDH hit - #pi^{+} track (phi) [radian]");
  diff2d_CDC_CDH_pip_phi_tof->SetYTitle("CDH hit - #pi^{+} TOF [nsec]");
  
  diff2d_CDC_CDH_pip_z_tof = new TH2F("diff2d_CDC_CDH_pip_z_tof","diff2d_CDC_CDH_pip_z_tof",100,-100,100,120,-10,50);
  diff2d_CDC_CDH_pip_z_tof->SetXTitle("CDH hit - #pi^{-} track (z) [cm]");
  diff2d_CDC_CDH_pip_z_tof->SetYTitle("CDH hit - #pi^{-} TOF [nsec]");
  
  dE_diffphi_CDC_CDH_pim = new TH2F("dE_diffphi_CDC_CDH_pim","dE_diffphi_CDC_CDH_pim",100,-1.*TMath::Pi(),TMath::Pi(),nbindE,0,50);
  dE_diffphi_CDC_CDH_pim->SetXTitle("CDH hit - #pi^{-} track (phi) [radian]");
  dE_diffphi_CDC_CDH_pim->SetYTitle("CDH dE [MeVee]");

  dE_diffphi_CDC_CDH_pip = new TH2F("dE_diffphi_CDC_CDH_pip","dE_diffphi_CDC_CDH_pip",100,-1.*TMath::Pi(),TMath::Pi(),nbindE,0,50);
  dE_diffphi_CDC_CDH_pip->SetXTitle("CDH hit - #pi^{+} track (phi) [radian]");
  dE_diffphi_CDC_CDH_pip->SetYTitle("CDH dE [MeVee]");

  //MMnpimissdE = new TH2F("MMnpimissdE","MMnpimissdE",nbindE,0,50,100,0.4,1.9);
  //MMnpimissdE->SetXTitle("CDH dE [MeVee]");
  //MMnpimissdE->SetYTitle("Miss. Mass [GeV/c^{2}]");

  MMnpimissdiffphi_CDC_CDH_pim
    = new TH2F("MMnpimissdiffphi_CDC_CDH_pim","MMnpimissdiffphi_CDC_CDH_pim",100,-1.*TMath::Pi(),TMath::Pi(),100,0.4,1.9);
  MMnpimissdiffphi_CDC_CDH_pim->SetXTitle("CDH hit - #pi^{-} track (phi) [radian]");
  MMnpimissdiffphi_CDC_CDH_pim->SetYTitle("Miss. Mass [GeV/c^{2}]");

  MMnpimissdiffphi_CDC_CDH_pip
    = new TH2F("MMnpimissdiffphi_CDC_CDH_pip","MMnpimissdiffphi_CDC_CDH_pip",100,-1.*TMath::Pi(),TMath::Pi(),100,0.4,1.9);
  MMnpimissdiffphi_CDC_CDH_pip->SetXTitle("CDH hit - #pi^{+} track (phi) [radian]");
  MMnpimissdiffphi_CDC_CDH_pip->SetYTitle("Miss. Mass [GeV/c^{2}]");

  pimmom_diffphi_CDC_CDH_pim = new TH2F("pimmom_diffphi_CDC_CDH_pim","pimmom_diffphi_CDC_CDH_pim",100,-1.*TMath::Pi(),TMath::Pi(), 200,0,1);
  pimmom_diffphi_CDC_CDH_pim->SetXTitle("CDH hit - #pi^{-} track (phi) [radian]");
  pimmom_diffphi_CDC_CDH_pim->SetYTitle("#pi^{-} mom. [GeV/c]");
  
  pipmom_diffphi_CDC_CDH_pip = new TH2F("pipmom_diffphi_CDC_CDH_pip","pipmom_diffphi_CDC_CDH_pip",100,-1.*TMath::Pi(),TMath::Pi(), 200,0,1);
  pipmom_diffphi_CDC_CDH_pip->SetXTitle("CDH hit - #pi^{+} track (phi) [radian]");
  pipmom_diffphi_CDC_CDH_pip->SetYTitle("#pi^{+} mom. [GeV/c]");
  
  pimmom_diffz_CDC_CDH_pim = new TH2F("pimmom_diffz_CDC_CDH_pim","pimmom_diffz_CDC_CDH_pim",100,-100,100, 200,0,1);
  pimmom_diffz_CDC_CDH_pim->SetXTitle("CDH hit - #pi^{-} track (z) [cm]");
  pimmom_diffz_CDC_CDH_pim->SetYTitle("#pi^{-} mom. [GeV/c]");

  pipmom_diffz_CDC_CDH_pip = new TH2F("pipmom_diffz_CDC_CDH_pip","pipmom_diffz_CDC_CDH_pip",100,-100,100, 200,0,1);
  pipmom_diffz_CDC_CDH_pip->SetXTitle("CDH hit - #pi^{+} track (z) [cm]");
  pipmom_diffz_CDC_CDH_pip->SetYTitle("#pi^{+} mom. [GeV/c]");

  nmom_diffphi_CDC_CDH_pim = new TH2F("nmom_diffphi_CDC_CDH_pim","nmom_diffphi_CDC_CDH_pim",100,-1.*TMath::Pi(),TMath::Pi(),100,0,1);
  nmom_diffphi_CDC_CDH_pim->SetXTitle("CDH hit - #pi^{-} track (phi) [radian]");
  nmom_diffphi_CDC_CDH_pim->SetYTitle("nmom.  [GeV/c]");

  nmom_diffphi_CDC_CDH_pip = new TH2F("nmom_diffphi_CDC_CDH_pip","nmom_diffphi_CDC_CDH_pip",100,-1.*TMath::Pi(),TMath::Pi(),100,0,1);
  nmom_diffphi_CDC_CDH_pip->SetXTitle("CDH hit - #pi^{+} track (phi) [radian]");
  nmom_diffphi_CDC_CDH_pip->SetYTitle("nmom.  [GeV/c]");

  nmom_diffz_CDC_CDH_pim = new TH2F("nmom_diffz_CDC_CDH_pim","nmom_diffz_CDC_CDH_pim",100,-100,100,100,0,1);
  nmom_diffz_CDC_CDH_pim->SetXTitle("CDH hit - #pi^{-} track (z) [cm]");
  nmom_diffz_CDC_CDH_pim->SetYTitle("nmom.  [GeV/c]");

  nmom_diffz_CDC_CDH_pip = new TH2F("nmom_diffz_CDC_CDH_pip","nmom_diffz_CDC_CDH_pip",100,-100,100,100,0,1);
  nmom_diffz_CDC_CDH_pip->SetXTitle("CDH hit - #pi^{+} track (z) [cm]");
  nmom_diffz_CDC_CDH_pip->SetYTitle("nmom.  [GeV/c]");

  //
  TH1F *nmom = new TH1F("nmom", "nmom", 50, 0, 1.0);
  nmom->SetXTitle("mom. [GeV/c]");

  TH2F *dE_nmom = new TH2F("dE_nmom", "dE_nmom", 50, 0, 1.0, 200, 0, 50);
  dE_nmom->SetXTitle("mom. [GeV/c]");
  dE_nmom->SetYTitle("dE. [MeVee]");

  TH1F *mnmom = new TH1F("mnmom", "mnmom", 100, 0, 2.0);
  mnmom->SetXTitle("mom. [GeV/c]");

  TH1F *npipmom = new TH1F("npipmom", "npipmom", 150, 0, 3.0);
  npipmom->SetXTitle("mom. [GeV/c]");

  TH1F *npimmom = new TH1F("npimmom", "npimmom", 150, 0, 3.0);
  npimmom->SetXTitle("mom. [GeV/c]");

  //DCA analysis
  TH1F* DCA_pip_beam = new TH1F("DCA_pip_beam","DCA_pip_beam",300,0,30);
  DCA_pip_beam->SetXTitle("DCA [cm]");
  TH1F* DCA_pim_beam = new TH1F("DCA_pim_beam","DCA_pim_beam",300,0,30);
  DCA_pim_beam->SetXTitle("DCA [cm]");
  TH1F* DCA_pip_pim = new TH1F("DCA_pip_pim","DCA_pip_pim",300,0,30);
  DCA_pip_pim->SetXTitle("DCA [cm]");

  std::cout << __LINE__ << std::endl;

  //MC info. for resolution evaluation

  TH2F* Vtx_ZX = new TH2F("Vtx_ZX","Vtx_ZX",1000,-25,25,500,-12.5,12.5);
  TH2F* Vtx_XY = new TH2F("Vtx_XY","Vtx_XY",500,-12.5,12.5,500,-12.5,12.5);

  std::cout << __LINE__ << std::endl;


  //reading TTree
  Int_t nevent = tree->GetEntries();
  std::cerr<<"# of events = "<<nevent<<std::endl;
  std::cout << "p-value cut:" << pvalcut << std::endl;
  std::cout << "dE cut:" << anacuts::dE_MIN << std::endl;
  TCanvas *cinfo = new TCanvas("cinfo","info");
  cinfo->cd();
  TPaveText *pt = new TPaveText(.05,.05,.95,.7);
  
  pt->AddText(Form("p-value cut: %f ",pvalcut));
  pt->AddText(Form("dE cut: %0.2f ",anacuts::dE_MIN));
  pt->AddText(Form("1/beta min.: %f ",1./anacuts::beta_MAX));
  pt->AddText(Form("1/beta max : %f ",1./anacuts::beta_MIN));
  pt->AddText(Form("K^{0} window : %0.3f - %0.3f",anacuts::pipi_MIN,anacuts::pipi_MAX ));
  pt->AddText(Form("#Sigma^{+} window : %0.3f - %0.3f (%0.1f sigma cut)",anacuts::Sigmap_MIN,anacuts::Sigmap_MAX,anacuts::Sigma_NSigmacut));
  pt->AddText(Form("#Sigma^{-} window : %0.3f - %0.3f (%0.1f sigma cut)",anacuts::Sigmam_MIN,anacuts::Sigmam_MAX,anacuts::Sigma_NSigmacut));
  pt->AddText(Form("miss. n window : %0.3f - %0.3f (%0.1f sigma cut)",anacuts::neutron_MIN,anacuts::neutron_MAX,anacuts::neutron_NSigmacut ));
  pt->AddText(Form("CDS_neutron mom. cut : %0.3f > GeV/c ",anacuts::nmomcut ));
  pt->Draw();

  //------------------------//
  //--- event roop start ---//
  //------------------------//
  for ( Int_t i=0; i<nevent; i++ ) {
    tree->GetEvent(i);
    if(i%50000==0) std::cout << "Event# " << i << std::endl;
    //TVector3 vtx_pi = *vtx_pi_cdc ;
    // calc missing npi+pi- mass //
    TLorentzVector LVec_pim;
    TLorentzVector LVec_pip;
    if(chargepi==1) LVec_pip = *LVec_pi;
    if(chargepi==0) LVec_pim = *LVec_pi;

    TLorentzVector LVec_nmiss = *LVec_target+*LVec_beam-*LVec_n;
    TLorentzVector LVec_pimiss = *LVec_target+*LVec_beam-*LVec_pi;
    TLorentzVector LVec_npimiss = *LVec_target+*LVec_beam-*LVec_pi-*LVec_n;
    double npimiss_mass = LVec_npimiss.M();
    double npimiss_mass2 = LVec_npimiss.M2();
    //if(npimiss_mass<0) continue;
    double npimiss_mom = LVec_npimiss.P();
    TLorentzVector LVec_pip_n;
    TLorentzVector LVec_pim_n;
    if(chargepi==1){
      LVec_pip_n = *LVec_n+*LVec_pi;
    }else if(chargepi==0){
      LVec_pim_n = *LVec_n+*LVec_pi;
    }
    //std::cout << __LINE__ << std::endl;
    // calc cos(theta) of missing npi+pi- //
    TVector3 boost = (*LVec_target+*LVec_beam).BoostVector();
    double cos_pimisslab = LVec_npimiss.Vect().Dot((*LVec_beam).Vect())/(LVec_npimiss.Vect().Mag()*(*LVec_beam).Vect().Mag());
    double pimissthetalab = acos(cos_pimisslab);
    TLorentzVector LVec_npimiss_CM = LVec_npimiss;
    TLorentzVector LVec_beam_CM = *LVec_beam;
    LVec_npimiss_CM.Boost(-boost);
    LVec_beam_CM.Boost(-boost);
    double cos_pimissCM = LVec_npimiss_CM.Vect().Dot(LVec_beam_CM.Vect())/(LVec_npimiss_CM.Vect().Mag()*LVec_beam_CM.Vect().Mag());
    //cos(theta) of n_cds in lab
    double cos_ncdslab =  LVec_n->Vect().Dot(LVec_beam->Vect())/(LVec_n->Vect().Mag()*LVec_beam->Vect().Mag());
    TLorentzVector LVec_n_CM = *LVec_n;
    LVec_n_CM.Boost(-boost);
    double cos_ncdsCM = LVec_n_CM.Vect().Dot(LVec_beam_CM.Vect())/(LVec_n_CM.Vect().Mag()*LVec_beam_CM.Vect().Mag());

    double cos_pi = (*LVec_pi).Vect().Dot((*LVec_beam).Vect())/((*LVec_pi).Vect().Mag()*(*LVec_beam).Vect().Mag());

    TLorentzVector LVec_pi_CM = *LVec_pi;
    LVec_pi_CM.Boost(-boost);
    double cos_piCM = LVec_pi_CM.Vect().Dot(LVec_beam_CM.Vect())/(LVec_pi_CM.Vect().Mag()*LVec_beam_CM.Vect().Mag());

    //std::cout << __LINE__ << std::endl;
    if(SimSpmode || SimSmmode || SimK0nmode) {
      TVector3 boost_mc =  (*LVec_target+*mcmom_beam).BoostVector();
    }
    TLorentzVector mcmom_pimiss_calc;
    if(SimSpmode || SimSmmode || SimK0nmode) {
      mcmom_pimiss_calc = *LVec_target+*mcmom_beam-*mcmom_pi-*mcmom_ncds;
    }

    // calc pi+n //
    TLorentzVector LVec_pi_n = *LVec_pi+*LVec_n;
    double phi_npi = (*LVec_n).Phi()-(*LVec_pi).Phi();
    double phi_pi = (*LVec_pi).Phi();
    double phi_n = (*LVec_n).Phi();

    if(phi_npi<-1.0*TMath::Pi()) phi_npi += 2.0*TMath::Pi();
    else if(phi_npi>1.0*TMath::Pi()) phi_npi -= 2.0*TMath::Pi();


    TLorentzVector LVec_pip_n_mc;
    TLorentzVector LVec_pim_n_mc;
    if(SimSpmode || SimSmmode || SimK0nmode) {
      if(mcchargepi==1){
        LVec_pip_n_mc  = *mcmom_pi+*mcmom_ncds;
      }else if(mcchargepi==0){
        LVec_pim_n_mc  = *mcmom_pi+*mcmom_ncds;
      }
    }

    TLorentzVector LVec_npimiss_mc;
    if(SimSpmode || SimSmmode || SimK0nmode) {
      LVec_npimiss_mc = *LVec_target+*mcmom_beam-*mcmom_pi-*mcmom_ncds;
    }

    //update the momentum 
    if( (SimSpmode || SimSmmode) && 
        SimRejectFake && 
      (mcpattern!=2)) continue;
    if(SimK0nmode && SimRejectFake && (mcpattern!=7) ) continue;

    //std::cout << __LINE__ << std::endl;
    TLorentzVector LVec_Sigma_react;
    TLorentzVector LVec_pi_react;//missing pi 
    TLorentzVector LVec_piSigma_react;
    bool IsFakeN1 = false;
    bool IsFakeN2 = false;
    bool IsFakebyVTX = false;
    if(SimSpmode || SimSmmode || SimK0nmode) {
      if( (mcncanvtxr>58) || (fabs(mcncanvtxz)>40) ){
      //if( (mcncanvtxr>15) || (fabs(mcncanvtxz)>15) ){
        IsFakebyVTX = true;
      }
      LVec_Sigma_react = *react_Sigma;
      LVec_pi_react = *react_pi;
      LVec_piSigma_react = LVec_Sigma_react + LVec_pi_react;
      TLorentzVector LVec_beam_r = *LVec_beam;
      double px = (*LVec_beam).Px();
      double py = (*LVec_beam).Py();
      double pz = (*LVec_beam).Pz();
      double E = (*LVec_beam).E();
      TLorentzVector LVec_beam_unit;//adjust unit [GeV/c^2, GeV/c]
      LVec_beam_unit.SetPx(px*1000.0);
      LVec_beam_unit.SetPy(py*1000.0);
      LVec_beam_unit.SetPz(pz*1000.0);
      LVec_beam_unit.SetE(E*1000.0);
      //std::cout << LVec_beam_r.P() << std::endl;
      //std::cout << LVec_beam_r.M() << std::endl;
      //std::cout << "beam p " << LVec_beam_unit.P() << std::endl;
      //std::cout << "beam m " << LVec_beam_unit.M() << std::endl;
      //std::cout << "nmiss P" <<  (*react_pimiss).P() << std::endl;
      //std::cout << __LINE__ << std::endl;
      vtxr_vtxz_mc->Fill(mcncanvtxz,mcncanvtxr);
      if(mcpattern==2){
        vtxr_vtxz_mc_pat2->Fill(mcncanvtxz,mcncanvtxr);
      }
      else if(mcpattern==7){
        vtxr_vtxz_mc_pat7->Fill(mcncanvtxz,mcncanvtxr);
      }
      
      if(!IsFakebyVTX  &&  (*mcmom_pimiss).P()>0.01 && (*mcmom_pimiss).P()<100.0 && (*mcmom_ncds).P()>0.01 && (*mcmom_pi).P()>0.01) {
        //true neutron from Sigma decay case
        if(mcpattern==2){
          IMnpim_MMnpim_mc_vtx_pat2->Fill(LVec_npimiss_mc.M(),LVec_pim_n_mc.M());
          IMnpip_MMnpip_mc_vtx_pat2->Fill(LVec_npimiss_mc.M(),LVec_pip_n_mc.M());
        }
        //true initial neutron case
        if(mcpattern==7){
          IMnpim_MMnpim_mc_vtx_pat7->Fill(LVec_npimiss_mc.M(),LVec_pim_n_mc.M());
          IMnpip_MMnpip_mc_vtx_pat7->Fill(LVec_npimiss_mc.M(),LVec_pip_n_mc.M());
        }
      }
      
      diff_nmiss_reactmc->Fill((*mcmom_pimiss).P()-(*react_pimiss).P()/1000.);
      diff2D_MMnpim_IMnpim_reactmc->Fill(LVec_pim_n_mc.M()-LVec_Sigma_react.M()/1000.,(*mcmom_pimiss).M()-(*react_pimiss).M()/1000.);
      diff2D_MMnpip_IMnpip_reactmc->Fill(LVec_pip_n_mc.M()-LVec_Sigma_react.M()/1000.,(*mcmom_pimiss).M()-(*react_pimiss).M()/1000.);
      nmom_IMnpip_mc->Fill((*mcmom_ncds).P(), LVec_pip_n_mc.M());
      nmom_IMnpim_mc->Fill((*mcmom_ncds).P(), LVec_pim_n_mc.M());
      //if(SimSpmode){
      //  if(fabs(LVec_pip_n_mc.M()-LVec_Sigma_react.M()/1000.) > 0.002) IsFakeN1 = true;
      //  if(fabs((*mcmom_pimiss).M()-(*react_pimiss).M()/1000.) > 0.002) IsFakeN1 = true;
      //}
      //if(SimSmmode){
      //  if(fabs(LVec_pim_n_mc.M()-LVec_Sigma_react.M()/1000.) > 0.002) IsFakeN1 = true;
      //  if(fabs((*mcmom_pimiss).M()-(*react_pimiss).M()/1000.) > 0.002) IsFakeN1 = true;
      //}
      //added angle check
      //if(!IsFakeN1){
      //  diff_cosnmiss_reactmc->Fill((*mcmom_pimiss).CosTheta()-(*react_pimiss).CosTheta());
      //  if(fabs((*mcmom_pimiss).CosTheta()-(*react_pimiss).CosTheta())>0.002) IsFakeN1 = true;
      //}
      //std::cout << __LINE__ << std::endl;
      if(SimSpmode){
        if(fabs(LVec_pip_n_mc.M()-LVec_Sigma_react.M()/1000.)>0.02) IsFakeN2=true;
        //std::cout << __LINE__ << std::endl;
        //if( (LVec_pip_n.M() -  LVec_pip_n_mc.M())< -0.012 || 0.010< (LVec_pip_n.M() -  LVec_pip_n_mc.M())) IsFakeN2=true;
        //if( fabs(LVec_pip_nmiss_mc.M()-1.18937)<0.01) IsFakeN2=true;
        //if( (diffIMnpip_recomc<-0.012) || (0.010<diffIMnpip_recomc)) IsFakeN2=true;
        //if(diffnpip_recomc.P()>0.10) IsFakeN2 = true;
      }
      if(SimSmmode){
        if(fabs(LVec_pim_n_mc.M()-LVec_Sigma_react.M()/1000.)>0.02) IsFakeN2=true;
        //if( (LVec_pim_n.M() -  LVec_pim_n_mc.M())< -0.012 || 0.010< (LVec_pim_n.M() -  LVec_pim_n_mc.M())) IsFakeN2=true;
        //if( fabs(LVec_pim_nmiss_mc.M()-1.197449)<0.01) IsFakeN2=true;
        //if(diffnpim_recomc.P()>0.10) IsFakeN2 = true;
      }

    }

    bool MissPiFlag=false;
    bool MissNwideFlag=false;
    bool NBetaOK=false;
    bool NdEOK=false;
    bool SigmaPFlag=false;
    bool SigmaMFlag=false;
    bool SigmawidePFlag=false;
    bool SigmawideMFlag=false;
    bool SigmaPMissNViciFlag=false;//select vicinity of the signal
    bool SigmaMMissNViciFlag=false;//select vicinity of the signal
    bool SigmaPMissNViciextFlag=false;//select vicinity of the signal
    bool SigmaMMissNViciextFlag=false;//select vicinity of the signal

    

    if( (LVec_pip_n.P()<anacuts::SigmaPMomCut) ){
      Vtx_ZX->Fill((*vtx_pi_cdc).Z(),(*vtx_pi_cdc).X());
      Vtx_XY->Fill((*vtx_pi_cdc).X(),(*vtx_pi_cdc).Y());
    }
    double dca_pip_beam = (*vtx_pi_beam-*vtx_pi_cdc).Mag();
    double dca_pim_beam = (*vtx_pi_beam-*vtx_pi_cdc).Mag();
    //double dca_pip_pim =(*CA_pip-*CA_pim).Mag();


    //-- neutron-ID, K0 and missing neutron selection --//
    if(anacuts::beta_MIN<NeutralBetaCDH &&  NeutralBetaCDH<anacuts::beta_MAX  ) NBetaOK=true;
    if(anacuts::dE_MIN<dE) NdEOK=true;
    double MassNPip= (*LVec_n+*LVec_pi).M();
    double MassNPim= (*LVec_n+*LVec_pi).M();

    TVector3 diffpim;
    double diffPhinpim = 0;
    double difftofnpim = 0;
    if(chargepi==0){
      diffpim = (*CDH_Pos)-(*CDH_Pos_pi);
      diffPhinpim = (*CDH_Pos).Phi()-(*CDH_Pos_pi).Phi();
      difftofnpim = tofn - tofpi;
    }
    if(diffPhinpim<-1.0*TMath::Pi()) diffPhinpim += 2.0*TMath::Pi();
    else if(diffPhinpim>1.0*TMath::Pi()) diffPhinpim -= 2.0*TMath::Pi();
    
    TVector3 diffpip; 
    double diffPhinpip = 0;
    double difftofnpip = 0;
    if(chargepi==1){
      diffpip = (*CDH_Pos)-(*CDH_Pos_pi);
      diffPhinpip = (*CDH_Pos).Phi()-(*CDH_Pos_pi).Phi();
      difftofnpip = tofn - tofpi;
    }
    if(diffPhinpip<-1.0*TMath::Pi()) diffPhinpip += 2.0*TMath::Pi();
    else if(diffPhinpip>1.0*TMath::Pi()) diffPhinpip -= 2.0*TMath::Pi();
    
    if(IsolationFlag==1 && (chargepi==0)) {
      //round cut
      if( pow((diffPhinpim-anacuts::Isonpim_shift)/anacuts::Isonpim_phicut,2.0)+pow(diffpim.Z()/anacuts::Isonpim_zcut,2.0) <1 ) continue;
      //for mixed events, avoid sharing same CDH segments
      if( -anacuts::CDHwidthphi< diffPhinpim  && diffPhinpim < anacuts::CDHwidthphi ) continue;
    }else if(IsolationFlag==2){ 
      //round cut wide
      if( pow((diffPhinpim-anacuts::Isonpim_shift)/anacuts::Isonpim_phicutwide,2.0)+pow(diffpim.Z()/anacuts::Isonpim_zcutwide,2.0) <1 ) continue;
      //for mixed events, avoid sharing same CDH segments
      if( -anacuts::CDHwidthphi< diffPhinpim  && diffPhinpim < anacuts::CDHwidthphi ) continue;
    } else if(IsolationFlag==3) {
      if( pow((diffPhinpim-anacuts::Isonpim_shift)/anacuts::Isonpim_phicut,2.0)+pow(diffpim.Z()/anacuts::Isonpim_zcut,2.0) >=1 ) continue;
      //for mixed events, avoid sharing same CDH segments
      if( -anacuts::CDHwidthphi< diffPhinpim  && diffPhinpim < anacuts::CDHwidthphi ) continue;
    }

    if(IsolationFlag==1 && (chargepi==1) ) {
      //round cut
      if( pow((diffPhinpip-anacuts::Isonpip_shift)/anacuts::Isonpip_phicut,2.0)+pow(diffpip.Z()/anacuts::Isonpip_zcut,2.0) <1 ) continue;
      //for mixed events, avoid sharing same CDH segments
      if( -anacuts::CDHwidthphi< diffPhinpip  && diffPhinpip < anacuts::CDHwidthphi ) continue;
    } else if(IsolationFlag==3) {
      if( pow((diffPhinpip-anacuts::Isonpip_shift)/anacuts::Isonpip_phicut,2.0)+pow(diffpip.Z()/anacuts::Isonpip_zcut,2.0)>=1 ) continue;
      //for mixed events, avoid sharing same CDH segments
      if( -anacuts::CDHwidthphi< diffPhinpip  && diffPhinpip < anacuts::CDHwidthphi ) continue;
    }

    if(CDCChargeVetoFlag && (nhitOutCDC!=0) ) continue;
    
    if(ForwardVetoFlag && ForwardCharge) continue;
    
    //std::cout << __LINE__ << std::endl;
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
    

    //vicinity of Sigma+ & NMiss events
    if(pow(((MassNPip - anacuts::Sigmap_center)/5.0/anacuts::Sigmap_sigma),2.0) +
       pow(((npimiss_mass - anacuts::neutron_center)/5.0/anacuts::neutron_sigma),2.0) < 1.0){
       if(npimiss_mass < 1.05){
         SigmaPMissNViciFlag=true;
       }
    }
    
    //vicinity of Sigma+ &Nmiss events +extended area
    if( (npimiss_mass >=anacuts::neutron_center) && (pow(((MassNPip - anacuts::Sigmap_center)/5.0/anacuts::Sigmap_sigma),2.0) +
       pow(((npimiss_mass - anacuts::neutron_center)/5.0/anacuts::neutron_sigma),2.0) < 1.0)){
       if(npimiss_mass < 1.05){
         SigmaPMissNViciextFlag=true;
       }
    }
    
    if( (npimiss_mass < anacuts::neutron_center) && (fabs( MassNPip - anacuts::Sigmap_center) < 5.0*anacuts::Sigmap_sigma)){
      SigmaPMissNViciextFlag=true;
    }


    //vicinity of Sigma- & NMiss events
    if( (npimiss_mass >= anacuts::neutron_center) && 
        (pow(((MassNPim - anacuts::Sigmam_center)/5.0/anacuts::Sigmam_sigma),2.0) +
        pow(((npimiss_mass - anacuts::neutron_center)/3.0/anacuts::neutron_sigma),2.0) < 1.0)){
        SigmaMMissNViciFlag=true;
        SigmaMMissNViciextFlag=true;
    }
    if( (npimiss_mass < anacuts::neutron_center) && 
        (pow(((MassNPim - anacuts::Sigmam_center)/5.0/anacuts::Sigmam_sigma),2.0) +
        pow(((npimiss_mass - anacuts::neutron_center)/5.0/anacuts::neutron_sigma),2.0) < 1.0)){
        SigmaMMissNViciFlag=true;
    }

    
    if( (npimiss_mass < anacuts::neutron_center) && (fabs( MassNPim - anacuts::Sigmam_center) < 5.0*anacuts::Sigmam_sigma)){
      SigmaMMissNViciextFlag=true;
    }


    if(anacuts::Miss2Pi_MIN<npimiss_mass2 && npimiss_mass2<anacuts::Miss2Pi_MAX ) MissPiFlag=true;
    if(anacuts::neutron_MIN_wide<npimiss_mass && npimiss_mass<anacuts::neutron_MAX_wide ) MissNwideFlag=true;

    //momentum update

    //if( (*LVec_n).P()<anacuts::nmomcut) continue;
    if(RejectStoppedSigma){
      if(LVec_pip_n.P()<anacuts::SigmaPMomCut && chargepi==1) continue;
      if(LVec_pim_n.P()<anacuts::SigmaMMomCut && chargepi==0) continue;
    }

    //std::cout << __LINE__ << std::endl;
    double weight = 1.0;
    if(MIXmode){
      weight = 4.24608060240400029e-02;
      if(SimSpmode){
        weight *=0.72;
        weight *=6.45779095649856028e-01;
      }
      if(SimSmmode){
        weight *=6.56913599999999986e-01;
        weight *=6.51604999999999879e-01;
      }
      if(SimK0nmode){
        weight *=5.69470101983999943e-01;
        weight *=7.37894527999999994e-01;
      }
    }
    static bool isState = false;
    if(!isState){
      if(SimSpmode) std::cout << "Sim Sp mode " << std::endl;
      if(SimSmmode) std::cout << "Sim Sm mode " << std::endl;
      if(SimK0nmode) std::cout << "Sim K0n mode " << std::endl;
      std::cout  << " weighting factor " << weight << std::endl;
      isState = true;
    }

    //---end of Flag definition-----------------------------------------------------
    //w/o kinfit
    //---including K0 --------------------------------------------------------------

    //std::cout << __LINE__ << std::endl;
    NHitCDCOut->Fill(nhitOutCDC,weight);
    IsForwardCharge->Fill(ForwardCharge,weight);
    CDHphi_betainv_fid->Fill(1./NeutralBetaCDH,(*CDH_Pos).Phi());
    CDHz_betainv_fid->Fill(1./NeutralBetaCDH,(*CDH_Pos).z());
    dE_betainv_fid->Fill(1./NeutralBetaCDH,dE);
    //pipmom_IMnpip->Fill(LVec_pip_n.M(),(LVec_pip).P());
    //pimmom_IMnpim->Fill(LVec_pim_n.M(),(LVec_pim).P());
    //std::cout << __LINE__ << std::endl;
    if(NBetaOK) {
      dE_nmom_fid_beta->Fill((*LVec_n).P(),dE);
      dE_MMom_fid_beta->Fill(LVec_npimiss.P(),dE);
      dE_MMass_fid_beta->Fill(LVec_npimiss.M(),dE);
      if(SigmaPFlag || SigmaMFlag) {
        dE_MMass_fid_beta_wSid->Fill(LVec_npimiss.M(),dE);
      }
      dE_CDHphi->Fill((*CDH_Pos).Phi(),dE);
      dE_CDHz->Fill((*CDH_Pos).z(),dE);
      dE_IMnpim->Fill(LVec_pim_n.M(),dE);
      dE_IMnpip->Fill(LVec_pip_n.M(),dE);
    }
    if(NBetaOK && MissPiFlag) {
      dE_IMnpim_pi->Fill(LVec_pim_n.M(),dE,weight);
      dE_IMnpip_pi->Fill(LVec_pip_n.M(),dE,weight);
    }
    if(NBetaOK && NdEOK) {

      CDHz_nmom_fid->Fill((*LVec_n).P(),(*CDH_Pos).z());
      MMom_MMass->Fill(LVec_npimiss.M(),LVec_npimiss.P(),weight);
      nmom_Momnpip->Fill(LVec_pip_n.P(),(*LVec_n).P());
      nmom_Momnpim->Fill(LVec_pim_n.P(),(*LVec_n).P());
      if(SigmaPFlag || SigmaMFlag) {
        MMom_MMass_wSid->Fill(LVec_npimiss.M(),LVec_npimiss.P(),weight);
        nmom_MMnpimisswSid->Fill(npimiss_mass,(*LVec_n).P(),weight);
        
        if(SimSpmode || SimSmmode || SimK0nmode){
          if(!( (mcpattern==2)  ||  (mcpattern==7))) {// || IsFakebyVTX )
            nmom_MMnpimisswSid_fake->Fill(npimiss_mass,(*LVec_n).P(),weight);
          }
        }
      }
      //pipmom_IMnpip_dE->Fill(LVec_pip_n.M(),(LVec_pip).P());
      //pimmom_IMnpim_dE->Fill(LVec_pim_n.M(),(LVec_pim).P());
      //pipmom_MMnpimissdE->Fill(npimiss_mass,(LVec_pip).P());
      //pimmom_MMnpimissdE->Fill(npimiss_mass,(LVec_pim).P());
      nmom_CDHphi->Fill((*CDH_Pos).Phi(),(*LVec_n).P());
      if(chargepi==0){//pi-
        diff2d_CDC_CDH_pim->Fill(diffPhinpim,diffpim.z());
        diff2d_CDC_CDH_pim_phi_tof->Fill(diffPhinpim,difftofnpim);
        diff2d_CDC_CDH_pim_z_tof->Fill(diffpim.z(),difftofnpim);
        dE_diffphi_CDC_CDH_pim->Fill(diffPhinpim,dE);
        MMnpimissdiffphi_CDC_CDH_pim->Fill(diffPhinpim, npimiss_mass);
        pimmom_diffphi_CDC_CDH_pim->Fill(diffPhinpim,(LVec_pim).P());
        pimmom_diffz_CDC_CDH_pim->Fill(diffpim.z(),(LVec_pim).P());
        nmom_diffphi_CDC_CDH_pim->Fill(diffPhinpim,(*LVec_n).P());
        nmom_diffz_CDC_CDH_pim->Fill(diffpim.z(),(*LVec_n).P());
        MMnpi_IMnpim->Fill(LVec_pim_n.M(),npimiss_mass);
        MM2npi_IMnpim->Fill(LVec_pim_n.M(),npimiss_mass2);
        Cospicm_IMnpim->Fill(LVec_pim_n.M(),cos_pimissCM);
        MMn_IMnpim->Fill(LVec_pim_n.M(),LVec_nmiss.M());
        MMpi_IMnpim->Fill(LVec_pim_n.M(),LVec_pimiss.M());
      }
      if(chargepi==1){//pi+
        diff2d_CDC_CDH_pip->Fill(diffPhinpip,diffpip.z());
        diff2d_CDC_CDH_pip_phi_tof->Fill(diffPhinpip,difftofnpip);
        diff2d_CDC_CDH_pip_z_tof->Fill(diffpip.z(),difftofnpip);
        dE_diffphi_CDC_CDH_pip->Fill(diffPhinpip,dE);
      //MMnpimissdE->Fill(dE, npimiss_mass);
        MMnpimissdiffphi_CDC_CDH_pip->Fill(diffPhinpip, npimiss_mass);
        pipmom_diffphi_CDC_CDH_pip->Fill(diffPhinpip,(LVec_pip).P());
        pipmom_diffz_CDC_CDH_pip->Fill(diffpip.z(),(LVec_pip).P());
        nmom_diffphi_CDC_CDH_pip->Fill(diffPhinpip,(*LVec_n).P());
        nmom_diffz_CDC_CDH_pip->Fill(diffpip.z(),(*LVec_n).P());
        MMnpi_IMnpip->Fill(LVec_pip_n.M(),npimiss_mass);
        MM2npi_IMnpip->Fill(LVec_pip_n.M(),npimiss_mass2);
        Cospicm_IMnpip->Fill(LVec_pip_n.M(),cos_pimissCM);
        MMn_IMnpip->Fill(LVec_pip_n.M(),LVec_nmiss.M());
        MMpi_IMnpip->Fill(LVec_pip_n.M(),LVec_pimiss.M());
      }
    } //if(NBetaOK && NdEOK) 

    //std::cout << __LINE__ << std::endl;


    //std::cout << __LINE__ << std::endl;
    //std::cout << __LINE__ << std::endl;
    if(NBetaOK && NdEOK && MissPiFlag) {
      nmom->Fill((*LVec_n).P());
      mnmom->Fill(npimiss_mom);
      dE_nmom->Fill((*LVec_n).P(),dE);
      if(chargepi==1)npipmom->Fill(LVec_pip_n.P());
      if(chargepi==0)npimmom->Fill(LVec_pim_n.P());
       
      if(chargepi==0){//pi-
        nmom_IMnpim_dE_n->Fill(LVec_pim_n.M(),(*LVec_n).P(),weight);
        MMn_IMnpim_pi->Fill(LVec_pim_n.M(),LVec_nmiss.M());
        MMpi_IMnpim_pi->Fill(LVec_pim_n.M(),LVec_pimiss.M());
        Cospicm_IMnpim_pi->Fill(LVec_pim_n.M(),cos_pimissCM);
      }else if(chargepi==1){
        MMn_IMnpip_pi->Fill(LVec_pip_n.M(),LVec_nmiss.M());
        MMpi_IMnpip_pi->Fill(LVec_pip_n.M(),LVec_pimiss.M());
        nmom_IMnpip_dE_n->Fill(LVec_pip_n.M(),(*LVec_n).P(),weight);
        Cospicm_IMnpip_pi->Fill(LVec_pip_n.M(),cos_pimissCM);
      }

      if(SimSpmode || SimSmmode){
        IMnpim_IMnpim_mc_dE_n->Fill(LVec_pim_n.M(),LVec_pim_n.M()-LVec_pim_n_mc.M(),weight);
        IMnpip_IMnpip_mc_dE_n->Fill(LVec_pip_n.M(),LVec_pip_n.M()-LVec_pip_n_mc.M(),weight);
        if((mcpattern==2) || (mcpattern==7)){
          IMnpim_IMnpim_mc_dE_n_vtx->Fill(LVec_pim_n.M(),LVec_pim_n.M()-LVec_pim_n_mc.M(),weight);
          IMnpip_IMnpip_mc_dE_n_vtx->Fill(LVec_pip_n.M(),LVec_pip_n.M()-LVec_pip_n_mc.M(),weight);
          if(mcpattern==2){
            IMnpim_IMnpim_mc_dE_n_vtx_pat2->Fill(LVec_pim_n.M(),LVec_pim_n.M()-LVec_pim_n_mc.M(),weight);
            IMnpip_IMnpip_mc_dE_n_vtx_pat2->Fill(LVec_pip_n.M(),LVec_pip_n.M()-LVec_pip_n_mc.M(),weight);
          }else if(mcpattern==7){
            IMnpim_IMnpim_mc_dE_n_vtx_pat7->Fill(LVec_pim_n.M(),LVec_pim_n.M()-LVec_pim_n_mc.M(),weight);
            IMnpip_IMnpip_mc_dE_n_vtx_pat7->Fill(LVec_pip_n.M(),LVec_pip_n.M()-LVec_pip_n_mc.M(),weight);
          }
        }
      }

      if(SigmaPFlag) {
        nmom_IMnpip_dE_n_Sp->Fill(LVec_pip_n.M(),(*LVec_n).P(),weight);
        nmom_Momnpip_n_Sp->Fill(LVec_pip_n.P(),(*LVec_n).P(),weight);
      }
      //std::cout << __LINE__ << std::endl;

      if(SigmaMFlag) {
        nmom_IMnpim_dE_n_Sm->Fill(LVec_pim_n.M(),(*LVec_n).P(),weight);
        nmom_Momnpim_n_Sm->Fill(LVec_pim_n.P(),(*LVec_n).P());
      }
      
      if(SigmaPFlag || SigmaMFlag) {
        //MMom_MMass_wSid_n->Fill(LVec_npipimiss.M(),LVec_npipimiss.P(),weight);
        nmom_cosn_wSid_n->Fill(cos_ncdslab,(*LVec_n).P(),weight);
        nmom_cosnmiss_wSid_n->Fill(cos_pimissCM,(*LVec_n).P(),weight);
        pipmom_Momnpip_wSid_n->Fill(LVec_pip_n.P(),(LVec_pip).P(),weight);
        pimmom_Momnpim_wSid_n->Fill(LVec_pim_n.P(),(LVec_pim).P(),weight);

        nmom_MMnpimisswSid_n->Fill(npimiss_mass,(*LVec_n).P(),weight);
        if(SigmaPFlag){
          pipmom_Momnpip_wSid_n_Sp->Fill(LVec_pip_n.P(),(LVec_pip).P() ,weight);
        }

        if(SigmaMFlag){
          pimmom_Momnpim_wSid_n_Sm->Fill(LVec_pim_n.P(),(LVec_pim).P(),weight);
        }

        //reaction data - mcData matching
        bool IsMissMassNOK = false;
        bool IsMcNMassOK = false;
        if(SimSpmode || SimSmmode || SimK0nmode){
          double diffIMnpim_reactmc = LVec_pim_n.M()- LVec_pim_n_mc.M();
          double diffIMnpip_reactmc = LVec_pip_n.M()- LVec_pip_n_mc.M();
          double diffMMnpimissrecomc = npimiss_mass - (*mcmom_pimiss).M();
          
          //MCdata missing mass check
          if(fabs((*mcmom_pimiss).M()-nMass)<0.01) IsMissMassNOK = true;
          //MCdata ncds check
          if(fabs((*mcmom_ncds).M()-nMass)<0.01) IsMcNMassOK = true;
          
          
          if(IsMissMassNOK && IsMcNMassOK){
            diff2D_IMnpim_Momnpim_wSid_n_reactmc->Fill(LVec_pim_n_mc.P()-LVec_Sigma_react.P()/1000.,LVec_pim_n_mc.M()-LVec_Sigma_react.M()/1000.);
            diff2D_IMnpip_Momnpip_wSid_n_reactmc->Fill(LVec_pip_n_mc.P()-LVec_Sigma_react.P()/1000.,LVec_pip_n_mc.M()-LVec_Sigma_react.M()/1000.);
            diff_nmiss_wSid_n_reactmc->Fill((*mcmom_pimiss).P()-(*react_pimiss).P()/1000.);
            diff2D_MMnpimissIMnpim_wSid_n_reactmc->Fill(LVec_pim_n_mc.M()-LVec_Sigma_react.M()/1000.,(*mcmom_pimiss).M()-(*react_pimiss).M()/1000.);
            diff2D_MMnpimissIMnpip_wSid_n_reactmc->Fill(LVec_pip_n_mc.M()-LVec_Sigma_react.M()/1000.,(*mcmom_pimiss).M()-(*react_pimiss).M()/1000.);
            //Sigma mass check and missing neurtom mom. check. 
            //note: Reaction data does not have the info. of the neutron from Sigma
            if(SimSpmode){
              if(fabs(LVec_pip_n_mc.M()-LVec_Sigma_react.M()/1000.) > 0.002) IsFakeN1 = true;
              if(fabs((*mcmom_pimiss).P()-(*react_pimiss).P()/1000.) > 0.002) IsFakeN1 = true;
            }
            if(SimSmmode){
              if(fabs(LVec_pim_n_mc.M()-LVec_Sigma_react.M()/1000.) > 0.002) IsFakeN1 = true;
              if(fabs((*mcmom_pimiss).P()-(*react_pimiss).P()/1000.) > 0.002) IsFakeN1 = true;
            }
            //added angle check->why angle ? ->mom. check. 
            //if(!IsFakeN1){
              diff2D_IMnpim_Momnpim_wSid_n_reactmc_fake1->Fill(LVec_pim_n_mc.P()-LVec_Sigma_react.P()/1000.,LVec_pim_n_mc.M()-LVec_Sigma_react.M()/1000.);
              diff2D_IMnpip_Momnpip_wSid_n_reactmc_fake1->Fill(LVec_pip_n_mc.P()-LVec_Sigma_react.P()/1000.,LVec_pip_n_mc.M()-LVec_Sigma_react.M()/1000.);
              diff_cosnmiss_wSid_n_reactmc->Fill((*mcmom_pimiss).CosTheta()-(*react_pimiss).CosTheta());
              if(fabs((*mcmom_pimiss).CosTheta()-(*react_pimiss).CosTheta())>0.002) IsFakeN1 = true;
            //}
          }
          
          //std::cout << __LINE__ << std::endl;

          double diffIMnpim_recomc = LVec_pim_n.M()- LVec_pim_n_mc.M();
          double diffIMnpip_recomc = LVec_pip_n.M()- LVec_pip_n_mc.M();
          TLorentzVector diffnpip_recomc = LVec_pip_n - LVec_pip_n_mc; 
          TLorentzVector diffnpim_recomc = LVec_pim_n - LVec_pim_n_mc; 
          TLorentzVector diffMMom_recomc = LVec_npimiss - *mcmom_pimiss;
          double diffnmom_recomc = (*LVec_n).P() - (*mcmom_ncds).P();
          double diffnpip_mcreact = LVec_pip_n_mc.P() - LVec_Sigma_react.P()/1000.0;
          double diffnpim_mcreact = LVec_pim_n_mc.P() - LVec_Sigma_react.P()/1000.0;
          double diffMassnpip_mcreact = LVec_pip_n_mc.M() - LVec_Sigma_react.M()/1000.0;
          double diffMassnpim_mcreact = LVec_pim_n_mc.M() - LVec_Sigma_react.M()/1000.0;
          //double dcareco = (*vtx_displaced-*vtx_reaction).Mag();
          //double dcarreco = (*vtx_displaced-*vtx_reaction).Perp();
          //double dcazreco = (*vtx_displaced-*vtx_reaction).z();
          double dcatrue = (*mc_disvtx-*mc_vtx).Mag();
          double dcartrue = (*mc_disvtx-*mc_vtx).Perp();
          double dcaztrue = (*mc_disvtx-*mc_vtx).z();
          /*
          double diffdisvtx_CApip = (*CA_pip-*mc_disvtx).Mag();
          double diffdisvtx_CApip_r = (*CA_pip-*mc_disvtx).Perp();
          double diffdisvtx_CApip_z = (*CA_pip-*mc_disvtx).z();
          double diffdisvtx_CApim = (*CA_pim-*mc_disvtx).Mag();
          double diffdisvtx_CApim_r = (*CA_pim-*mc_disvtx).Perp();
          double diffdisvtx_CApim_z = (*CA_pim-*mc_disvtx).z();
          double diffdisvtx_cdcbeam_pip = (*vtx_pip_cdc-*mc_disvtx).Mag();
          double diffdisvtx_cdcbeam_pip_r = (*vtx_pip_cdc-*mc_disvtx).Perp();
          double diffdisvtx_cdcbeam_pip_z = (*vtx_pip_cdc-*mc_disvtx).z();
          double diffdisvtx_cdcbeam_pim = (*vtx_pim_cdc-*mc_disvtx).Mag();
          double diffdisvtx_cdcbeam_pim_r = (*vtx_pim_cdc-*mc_disvtx).Perp();
          double diffdisvtx_cdcbeam_pim_z = (*vtx_pim_cdc-*mc_disvtx).z();
          */
          vtxr_vtxz_n->Fill((*vtx_reaction).z(),(*vtx_reaction).Perp());
          //disvtxr_disvtxz_n->Fill((*vtx_displaced).z(),(*vtx_displaced).Perp());
          //disvtxr_disvtxz_n->Fill((*vtx_pip_cdc).z(),(*vtx_pip_cdc).Perp());
          mcvtxr_mcvtxz_n_mc->Fill((*mc_vtx).z(),(*mc_vtx).Perp() );
          mcdisvtxr_mcdisvtxz_n_mc->Fill((*mc_disvtx).z(),(*mc_disvtx).Perp() );
          diff2D_MMnpimissIMnpim_recomc_wSid_n->Fill(diffIMnpim_recomc,diffMMnpimissrecomc);
          diff2D_MMnpimissIMnpip_recomc_wSid_n->Fill(diffIMnpip_recomc,diffMMnpimissrecomc);
          diff2D_nmom_IMnpim_recomc_wSid_n->Fill(diffIMnpim_recomc,diffnmom_recomc);
          diff2D_nmom_IMnpip_recomc_wSid_n->Fill(diffIMnpip_recomc,diffnmom_recomc);
          diffMMom_recomc_wSid_n->Fill(diffMMom_recomc.P());
          vtxr_generation_ncan_wSid_n_mc->Fill(mcncdsgen,mcncanvtxr);
          //diffnmom_diffdca_n->Fill(dcareco-dcatrue,diffnmom_recomc);
          //diffnmom_diffdcar_n->Fill(dcarreco-dcartrue,diffnmom_recomc);
          //diffnmom_diffdcaz_n->Fill(dcazreco-dcaztrue,diffnmom_recomc);
          /*
          diffnmom_diffdisvtx_CApip_n->Fill(diffdisvtx_CApip,diffnmom_recomc);
          diffnmom_diffdisvtx_CApip_r_n->Fill(diffdisvtx_CApip_r,diffnmom_recomc);
          diffnmom_diffdisvtx_CApip_z_n->Fill(diffdisvtx_CApip_z,diffnmom_recomc);
          diffnmom_diffdisvtx_CApim_n->Fill(diffdisvtx_CApim,diffnmom_recomc);
          diffnmom_diffdisvtx_CApim_r_n->Fill(diffdisvtx_CApim_r,diffnmom_recomc);
          diffnmom_diffdisvtx_CApim_z_n->Fill(diffdisvtx_CApim_z,diffnmom_recomc);

          diffnmom_diffdisvtx_cdcbeam_pip_n->Fill(diffdisvtx_cdcbeam_pip,diffnmom_recomc);
          diffnmom_diffdisvtx_cdcbeam_pip_r_n->Fill(diffdisvtx_cdcbeam_pip_r,diffnmom_recomc);
          diffnmom_diffdisvtx_cdcbeam_pip_z_n->Fill(diffdisvtx_cdcbeam_pip_z,diffnmom_recomc);
          diffnmom_diffdisvtx_cdcbeam_pim_n->Fill(diffdisvtx_cdcbeam_pim,diffnmom_recomc);
          diffnmom_diffdisvtx_cdcbeam_pim_r_n->Fill(diffdisvtx_cdcbeam_pim_r,diffnmom_recomc);
          diffnmom_diffdisvtx_cdcbeam_pim_z_n->Fill(diffdisvtx_cdcbeam_pim_z,diffnmom_recomc);
          */
          
          if(!IsFakebyVTX){
            vtxr_diffmom_npip_ncan_wSid_n_mc->Fill(diffnpip_mcreact,mcncanvtxr);
            vtxr_diffmom_npim_ncan_wSid_n_mc->Fill(diffnpim_mcreact,mcncanvtxr);
            vtxr_diffMass_npip_ncan_wSid_n_mc->Fill(diffMassnpip_mcreact,mcncanvtxr);
            vtxr_diffMass_npim_ncan_wSid_n_mc->Fill(diffMassnpim_mcreact,mcncanvtxr);
            generation_diffmom_npip_ncan_wSid_n_mc->Fill(diffnpip_mcreact,mcncdsgen);
            generation_diffmom_npim_ncan_wSid_n_mc->Fill(diffnpim_mcreact,mcncdsgen);
            generation_diffMass_npip_ncan_wSid_n_mc->Fill(diffMassnpip_mcreact,mcncdsgen);
            generation_diffMass_npim_ncan_wSid_n_mc->Fill(diffMassnpim_mcreact,mcncdsgen);
          }

          vtxr_vtxz_ncan_wSid_n_mc->Fill(mcncanvtxz,mcncanvtxr);
          if(mcpattern==2)vtxr_vtxz_ncan_wSid_n_mc_pat2->Fill(mcncanvtxz,mcncanvtxr);
          else if(mcpattern==7)vtxr_vtxz_ncan_wSid_n_mc_pat7->Fill(mcncanvtxz,mcncanvtxr);
          
          if(!IsFakebyVTX){
            diff2D_MMnpimissIMnpim_recomc_wSid_n_fake1->Fill(diffIMnpim_recomc,diffMMnpimissrecomc);
            diff2D_MMnpimissIMnpip_recomc_wSid_n_fake1->Fill(diffIMnpip_recomc,diffMMnpimissrecomc);
            if(IsMcNMassOK){
              diff2D_nmom_IMnpim_recomc_wSid_n_fake1->Fill(diffIMnpim_recomc,diffnmom_recomc);
              diff2D_nmom_IMnpip_recomc_wSid_n_fake1->Fill(diffIMnpip_recomc,diffnmom_recomc);
            }
          }
          if(!( (mcpattern==2)  ||  (mcpattern==7))) {// || IsFakebyVTX )
            IMnpim_MMnpim_mc_wSid_n_fake->Fill(LVec_npimiss_mc.M(),LVec_pim_n_mc.M());
            IMnpip_MMnpip_mc_wSid_n_fake->Fill(LVec_npimiss_mc.M(),LVec_pip_n_mc.M());
          }
        }//SimSpmode or SimSmmode
      }//wSid_n
    }
    //std::cout << __LINE__ << std::endl;
  }//for ievt
  //--- Filling Histogram END --------------------------------------------------


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

  //centering title of all histograms
  f->cd();
  TIter nexthist(gDirectory->GetList());
  TH1F *h1 = nullptr;
  TH1D *h1d = nullptr;
  TH2F *h2 = nullptr;
  TObject *obj = nullptr;
  while( (obj = (TObject*)nexthist())!=nullptr  ) {
    if(obj->InheritsFrom("TH1F")) {
      h1 = (TH1F*) obj;
      h1->GetXaxis()->CenterTitle();
      //h1->GetXaxis()->SetTitleSize(0.05);
      //h1->GetXaxis()->SetTitleOffset(0.80);
      h1->GetYaxis()->SetTitleOffset(1.4);
    }
    if(obj->InheritsFrom("TH1D")) {
      h1d = (TH1D*) obj;
      h1d->GetXaxis()->CenterTitle();
      h1d->GetYaxis()->CenterTitle();
      //h1d->GetXaxis()->SetTitleSize(0.05);
      //h1d->GetXaxis()->SetTitleOffset(0.80);
      h1d->GetYaxis()->SetTitleOffset(1.4);
    }
    if(obj->InheritsFrom("TH2")) {
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
  for(int i=0; i<size; i++) {
    //while((c= (TCanvas*)next()))
    //pdf->NewPage();
    c= (TCanvas*)next();
    c->Draw();
    c->cd();
    TPaveText *pt;
    if(SimSpmode || SimSmmode) {
      pt = new TPaveText(.80,0.90,0.98,0.99,"NDC");
    } else {
      pt = new TPaveText(.76,0.90,0.90,0.98,"NDC");
      //pt = new TPaveText(.74,.81,0.9,0.90,"NDC");
    }
    if(SimSpmode) {
      pt->SetFillColor(kAzure-4);
      pt->AddText("MC #Sigma+#pi- mode");
    }
    else if(SimSmmode) {
      pt->SetFillColor(kAzure-4);
      pt->AddText("MC #Sigma-#pi+ mode");
    }
    else if(SimK0nmode) {
      pt->SetFillColor(kAzure-4);
      pt->AddText("MC K0nn");
    }else {
      pt->AddText("Real Data");
      pt->SetFillColor(kCyan-9);
    }
    pt->SetBorderSize(1);
    pt->Draw();

    c->Modified();
    c->Update();
    if(i==0) c->Print(pdfname+"(",Form("Title:%s",c->GetTitle()));
    else if(i==size-1)c->Print(pdfname+")",Form("Title:%s",c->GetTitle()));
    else c->Print(pdfname,Form("Title:%s",c->GetTitle()));
  }
  std::cout << "closing pdf " << std::endl;

  TIter nexthist2(gDirectory->GetList());
  TString outname = std::string(filename);
  if(IsolationFlag==0) outname.Replace(std::string(filename).size()-5,5,"_out_noiso.root");
  else if(IsolationFlag==1) outname.Replace(std::string(filename).size()-5,5,"_out_iso.root");
  else if(IsolationFlag==2) outname.Replace(std::string(filename).size()-5,5,"_out_isowide.root");
  else if(IsolationFlag==3) outname.Replace(std::string(filename).size()-5,5,"_out_isorev.root");
  
  if(SimRejectFake && (SimSpmode || SimSmmode || SimK0nmode)){
    outname.Replace(std::string(outname).size()-5,5,"_rej.root");
  }
  
  if(RejectStoppedSigma && (RealDatamode || SimSpmode || SimSmmode || SimK0nmode)){
    outname.Replace(std::string(outname).size()-5,5,"_nostop.root");
  }
    
  TFile *fout = new TFile(outname.Data(),"RECREATE");
  fout->Print();
  fout->cd();
  while( (obj = (TObject*)nexthist2())!=nullptr) {
    obj->Write();
  }
  fout->cd();
  fout->Close();

}

