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
#include <TPaveText.h>
#include <TLatex.h>
#include <TEfficiency.h>
#include <TGeoManager.h>

#include "../src/GlobalVariables.h"
#include "anacuts.h"

const double pvalcut = 0.005;
//const double pvalcut = 1.0e-5;
const bool gridon=true;
const bool staton=false;
const bool UseKinFitVal = false;
const double ForwardAngle=0.996;
const double MallCut= -0.2;
const int ncut=6;

int GetID(const TVector3 &pos,TGeoManager *geom);
//mode 0: Sigma+ ,1: Sigma-
void plot_IMLambdaPim(const char* filename="", const int qvalcutflag=0)
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
//  gStyle->SetMarkerStyle(20);
//  gStyle->SetMarkerSize(1.2);
  gStyle->SetCanvasDefH(800);
  gStyle->SetCanvasDefW(900);
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadLeftMargin(0.12);
  //gStyle->SetTitleFontSize(0.1);


  std::cout << "infile " << filename <<std::endl;
  TString pdfname = std::string(filename);
  if(qvalcutflag==0) pdfname.Replace(std::string(filename).size()-4,5,"pdf");
  if(qvalcutflag==1) pdfname.Replace(std::string(filename).size()-5,8,"_1.pdf");
  if(qvalcutflag==2) pdfname.Replace(std::string(filename).size()-5,8,"_2.pdf");
  std::cout << "pdfname: " << pdfname << std::endl;
  std::cout << std::endl;
  std::cout << "Use Kin Fit Val ? " << std::endl;

  bool SimMode = (std::string(filename).find("sim")!= std::string::npos);

  if(UseKinFitVal) std::cout << "Yes" << std::endl;
  else             std::cout << "No"  << std::endl;

  TH1::SetDefaultSumw2();
  //= = = = pppi-pi- final-sample tree = = = =//
  TLorentzVector *LVec_beam = NULL;   // 4-momentum(beam)
  TLorentzVector *LVec_target = NULL; // 4-momentum(target)
  TLorentzVector *LVec_pim1 = NULL;    // 4-momentum(pi+)
  TLorentzVector *LVec_pim2 = NULL;    // 4-momentum(pi-)
  TLorentzVector *LVec_p = NULL;      // 4-momentum(proton)
  TLorentzVector *LVec_p2 = NULL;      // 4-momentum(proton)
  TLorentzVector *mcmom_beam = NULL;   // true 4-momentum(beam)
  TLorentzVector *mcmom_pim1 = NULL;    // true 4-momentum(pi-)
  TLorentzVector *mcmom_pim2 = NULL;    // true 4-momentum(pi-),always decay pion from Lambda
  TLorentzVector *mcmom_p = NULL;      // true 4-momentum(proton),always decay proton from Lambda
  TLorentzVector *mcmom_pmiss = NULL;      // true 4-momentum(proton)
  TLorentzVector *react_pmiss = NULL;     //generated initial proton
  TLorentzVector *react_Lambda = NULL;   // generated Lambda
  TLorentzVector *react_pim = NULL; // generated pim
  TVector3 *vtx_reaction = NULL; // vertex(reaction)
  TVector3 *vtx_displaced = NULL; // vertex(reaction)
  TVector3 *vtx_pim1_beam = NULL; //C.A.P of pip-beam beam side
  TVector3 *vtx_pim2_beam = NULL; //C.A.P of pim-beam beam side
  TVector3 *vtx_p_beam = NULL; //C.A.P of pim-beam beam side
  TVector3 *vtx_pim1_cdc = NULL;//C.A.P of pip-beam pip side
  TVector3 *vtx_pim2_cdc = NULL;//C.A.P of pim-beam pim side
  TVector3 *vtx_p_cdc = NULL;//C.A.P of pim-beam pim side
  TVector3 *CA_pim1_pim1p = NULL;//C.A.P of pip-pim pip side
  TVector3 *CA_p_pim1p = NULL;//C.A.P of pip-pim pim side
  TVector3 *CA_pim2_pim2p = NULL;//C.A.P of pip-pim pip side
  TVector3 *CA_p_pim2p = NULL;//C.A.P of pip-pim pim side
  TVector3 *CA_pim1_pim1pim2 = NULL;//C.A.P of pip-pim pim side
  TVector3 *CA_pim2_pim1pim2 = NULL;//C.A.P of pip-pim pim side
  TVector3 *vtx_Lcan_p_pim1 = NULL;
  TVector3 *vtx_Lcan_p_pim2 = NULL;
  TVector3 *vtx_Lcan_p2_pim1 = NULL;
  TVector3 *vtx_Lcan_p2_pim2 = NULL;
  int ForwardCharge=0;
  //int run_num;   // run number
  //int event_num; // event number
  //int block_num; // block number
  TLorentzVector *kf_mom_beam = NULL;   // 4-momentum(beam) after kinematical refit for pi- Sigma+
  TLorentzVector *kf_mom_pim1 = NULL;    // 4-momentum(pi+) after kinematical refit for pi- Sigma+
  TLorentzVector *kf_mom_pim2 = NULL;    // 4-momentum(pi-) after kinematical refit for pi- Sigma+
  TLorentzVector *kf_mom_proton = NULL;      // 4-momentum(neutron) after kinematical refit for pi- Sigma+
  double kf_chi2;   // chi2 of kinematical refit
  double kf_NDF;    // NDF of kinematical refit
  double kf_status; // status of kinematical refit -> details can be found in this code
  double kf_pvalue; // p-value of kinematical refit
  int kf_flag; // flag of correct pair reconstruction, etc

  //= = = = pipipnn final-sample tree = = = =//

  TFile *f = new TFile(filename);
  //TFile *f = new TFile("sim_piSpn_dE0_Al.root");
  TTree *tree = (TTree*)f->Get("EventTree");
  if(tree==0) {
    std::cout << "EventTree is not found " << std::endl;
    return ;
  }

  tree->SetBranchAddress( "mom_beam",   &LVec_beam );
  tree->SetBranchAddress( "mom_target", &LVec_target );
  tree->SetBranchAddress( "mom_pim1", &LVec_pim1 );
  tree->SetBranchAddress( "mom_pim2", &LVec_pim2 );
  tree->SetBranchAddress( "mom_p", &LVec_p );
  tree->SetBranchAddress( "mom_p2", &LVec_p2 );
  tree->SetBranchAddress( "vtx_reaction", &vtx_reaction );
  tree->SetBranchAddress( "vtx_displaced", &vtx_displaced );
  tree->SetBranchAddress( "vtx_pim1_beam",&vtx_pim1_beam);
  tree->SetBranchAddress( "vtx_pim2_beam",&vtx_pim2_beam);
  tree->SetBranchAddress( "vtx_p_beam",&vtx_p_beam);
  tree->SetBranchAddress( "vtx_pim1_cdc",&vtx_pim1_cdc);
  tree->SetBranchAddress( "vtx_pim2_cdc",&vtx_pim2_cdc);
  tree->SetBranchAddress( "vtx_p_beam", &vtx_p_beam);
  tree->SetBranchAddress( "vtx_p_cdc", &vtx_p_cdc);
  tree->SetBranchAddress( "CA_pim1_pim1p", &CA_pim1_pim1p);
  tree->SetBranchAddress( "CA_p_pim1p", &CA_p_pim1p);
  tree->SetBranchAddress( "CA_pim2_pim2p", &CA_pim2_pim2p);
  tree->SetBranchAddress( "CA_p_pim2p", &CA_p_pim2p);
  tree->SetBranchAddress( "CA_pim1_pim1pim2", &CA_pim1_pim1pim2);
  tree->SetBranchAddress( "CA_pim2_pim1pim2", &CA_pim2_pim1pim2);
  tree->SetBranchAddress( "vtx_Lcan_p_pim1",&vtx_Lcan_p_pim1);
  tree->SetBranchAddress( "vtx_Lcan_p_pim2",&vtx_Lcan_p_pim2);
  tree->SetBranchAddress( "vtx_Lcan_p2_pim1",&vtx_Lcan_p2_pim1);
  tree->SetBranchAddress( "vtx_Lcan_p2_pim2",&vtx_Lcan_p2_pim2);
  tree->SetBranchAddress( "ForwardCharge",&ForwardCharge);
  //tree->SetBranchAddress( "run_num", &run_num );
  //tree->SetBranchAddress( "event_num", &event_num );
  //tree->SetBranchAddress( "block_num", &block_num );
  if(SimMode) {
    tree->SetBranchAddress( "mcmom_beam",  &mcmom_beam );
    tree->SetBranchAddress( "mcmom_pim1", &mcmom_pim1);
    tree->SetBranchAddress( "mcmom_pim2", &mcmom_pim2);//always decay pi- from Lambda
    tree->SetBranchAddress( "mcmom_p", &mcmom_p);//always decay pi- from Lambda
    tree->SetBranchAddress( "mcmom_pmiss", &mcmom_pmiss);
    tree->SetBranchAddress( "react_pmiss",&react_pmiss);
    tree->SetBranchAddress( "react_Lambda",&react_Lambda);
    tree->SetBranchAddress( "react_pim",&react_pim);
  }

  //tree->SetBranchAddress( "kfSpmode_mom_beam",   &kfSpmode_mom_beam );
  //tree->SetBranchAddress( "kfSpmode_mom_pip", &kfSpmode_mom_pip );
  //tree->SetBranchAddress( "kfSpmode_mom_pim", &kfSpmode_mom_pim );
  //tree->SetBranchAddress( "kfSpmode_mom_n", &kfSpmode_mom_n );
  //tree->SetBranchAddress( "kfSpmode_chi2", &kfSpmode_chi2 );
  //tree->SetBranchAddress( "kfSpmode_NDF", &kfSpmode_NDF );
  //tree->SetBranchAddress( "kfSpmode_status", &kfSpmode_status );
  //tree->SetBranchAddress( "kfSpmode_pvalue", &kfSpmode_pvalue );
  //tree->SetBranchAddress( "kfSmmode_mom_beam",   &kfSmmode_mom_beam );
  //tree->SetBranchAddress( "kfSmmode_mom_pip", &kfSmmode_mom_pip );
  //tree->SetBranchAddress( "kfSmmode_mom_pim", &kfSmmode_mom_pim );
  //tree->SetBranchAddress( "kfSmmode_mom_n", &kfSmmode_mom_n );
  //tree->SetBranchAddress( "kfSmmode_chi2", &kfSmmode_chi2 );
  //tree->SetBranchAddress( "kfSmmode_NDF", &kfSmmode_NDF );
  //tree->SetBranchAddress( "kfSmmode_status", &kfSmmode_status );
  //tree->SetBranchAddress( "kfSmmode_pvalue", &kfSmmode_pvalue );
  //tree->SetBranchAddress( "kf_flag", &kf_flag );


  // w/o kinematic fit
  f->cd();
  TH2F* MMom_MMass;
  TH2F* MMom_MMass_2;
  TH2F* q_MMass;
  TH2F* q_MMass_forward;
  TH2F* MMom_MMass_p;
  TH2F* MMom_MMass_p2;
  TH2F* MMom_MMass_p_wL;
  TH2F* MMom_MMass_p_wL_nop2;
  TH2F* MMom_MMass_p_wL_wp2;
  TH2F* q_PMom;
  TH2F* q_MMom;
  TH2F* q_MMomCosTheta;
  TH2F* q_MMomCosTheta_nop2;
  TH2F* q_MMomCosTheta_wp2;
  TH2F* q_pmisscos;
  TH2F* q_PMom_p_wL;
  TH2F* q_LMom_p_wL;
  TH2F* q_MMom_p_wL;
  TH2F* q_MMomCosTheta_p_wL;
  TH2F* q_MMomCosTheta_p_wL_nop2;
  TH2F* q_MMomCosTheta_p_wL_wp2;
  TH2F* q2_MMom2CosTheta_p2_wL_wp2;
  TH2F* q_pmiss_p_wL;
  TH2F* q_P2Mom;
  TH2F* q_P2Mom_p2_wL;
  TH2F* q_L2Mom_p2_wL;
  TH2F* q_PMom_p_wL_sum;
  TH2F* MMom_PMom;
  TH2F* MMom_PMom_2;
  TH2F* IMppim1_IMppim2;
  TH2F* IMppim1_IMp2pim1;
  TH2F* IMppim2_IMp2pim2;
  TH2F* IMp2pim1_IMp2pim2;
  TH2F* IMppim1_IMppim2_p;
  TH2F* IMp2pim1_IMp2pim2_p2;
  TH2F* IMppim1_IMppim2_p_wL;
  TH2F* IMp2pim1_IMp2pim2_p2_wL;
  TH2F* IMMissppim1_IMMissppim2_p_wL;
  TH2F* MMass_IMppim1;
  TH2F* MMass_IMppim2;
  TH2F* MMass_IMppim_wL;
  TH2F* MMass_IMppim_p_wL;
  TH2F* MMass_IMp2pim_p2_wL;
  TH2F* MMass_IMppipi_wL;
  TH2F* MMass_IMp2pipi_wL;
  TH2F* MMass_IMppipi_wL_sum;
  TH2F* MMass_IMppipi_wL_sum_forward;
  TH1F* MpimfromL_p_sum;
  TH2F* Mppim1_Mppim2_p_sum;
  TH2F* Mppim1_Mppim2_p_wL_sum;
  TH1F* MppimnotL_sum;
  TH2F* MML_MmomL_p_wL;
  TH2F* q_IMppim_mc;
  TH2F* q_IMppipi_mc;
  TH2F* q_IMppipi_p_wL_mc;
  TH2F* q_IMppipi_p_wL_wp2_mc;
  TH2F* q_IMppipi_p_wL_nop2_mc;
  TH2F* q_IMppipi_p;
  TH2F* q_IMp2pipi_p2;
  TH2F* q_IMppipi_p_wL;
  TH2F* q_IMp2pipi_p2_wL;
  TH2F* CosTheta_IMppipi_p_wL;
  TH2F* CosTheta_IMppipi_p_wL_nocombi;
  TH2F* CosThetaCM_IMppipi_p_wL_nocombi;
  TH2F* CosTheta_IMppipi_p_wL_nocombi_mc;
  TH2F* CosTheta_IMp2pipi_p2_wL;
  TH2F* q_IMppipi_p_wL_sum[ncut];//0:default, 1 half low, 2 half high, 3 sigma0 region, 4 wide range
  TH2F* q_IMppipi_p_wL_sum_nop2[ncut];//0:default, 1 half low, 2 half high, 3 sigma0 region, 4 wide range
  TH2F* q_IMpppipi_p_wL_sum_nop2_addptest;//0:default, 1 half low, 2 half high, 3 sigma0 region, 4 wide range
  TH2F* q_IMppipi_p_wL_sum_wp2[ncut];//0:default, 1 half low, 2 half high, 3 sigma0 region, 4 wide range
  TH2F* q_IMppipi_p_wL_sum_nocombi[ncut];
  TH2F* q_IMppipi_p_wL_sum_nocombi_nop2[ncut];
  TH2F* q_IMppipi_p_wL_sum_nocombi_wp2[ncut];
  TH2F* q_IMppipi_p_wL_sum_nocombi_mc[ncut];
  TH2F* q_IMppipi_p_wL_sum_nocombi_nop2_mc[ncut];
  TH2F* q_IMppipi_p_wL_sum_nocombi_wp2_mc[ncut];
  TH2F* q_IMppipi_pboth_wL;
  TH2F* q_IMppipi_p_wL_wp2;
  TH2F* CosTheta_IMppipi_p_wL_sum;
  TH2F* CosThetaCM_IMppipi_p_wL_sum;
  TH2F* CosTheta_IMppipi_p_wL_sum_mc;//mc true val
  TH2F* q_IMppipi_p_wL_sum_forward;
  TH2F* q_IMppipi_p_wL_sum_forward_plus;
  TH2F* q_IMppipi_p_wL_sum_forward_minus;
  TH2F* CosTheta_IMppipi_p_wL_sum_forward;
  TH2F* q_IMppipi_p_wL_sum_fp;
  TH2F* CosTheta_IMppipi_p_wL_sum_fp;
  TH1F* IMpL_p_wL_wp2;
  TH2F* IMmisspL_IMppipi_p_wL;
  TH2F* IMmisspL_q_p_wL;
 
  TH1F* DCA_pim1_beam = new TH1F("DCA_pim1_beam","DCA_pim1_beam",300,0,30);
  DCA_pim1_beam->SetXTitle("DCA #pi^{-}1 [cm]");
  TH1F* DCA_pim2_beam = new TH1F("DCA_pim2_beam","DCA_pim2_beam",300,0,30);
  DCA_pim2_beam->SetXTitle("DCA #pi^{-}2 [cm]");
  TH1F* DCA_pim1_pim2 = new TH1F("DCA_pim1_pim2","DCA_pim1_pim2",300,0,30);
  DCA_pim1_pim2->SetXTitle("DCA #pi^{-}1#pi^{-}2");
  TH1F* BeamMom = new TH1F("BeamMom","Beam Mom.",100,0.9,1.1);
  BeamMom->SetXTitle("Beam Mom. [GeV/c]");
  TH1F* BeamMomCosTheta = new TH1F("BeamMomCosTheta","Beam Mom. CosTheta",10000,0,1);
  BeamMomCosTheta->SetXTitle("Beam Mom. CosTheta");

  TH2D* Vtx_ZX_Lcan = new TH2D("Vtx_ZX_Lcan","Vtx_ZX_Lcan",500,-25,25,250,-12.5,12.5);
  Vtx_ZX_Lcan->SetXTitle("vertex Z [cm]");
  Vtx_ZX_Lcan->SetYTitle("vertex X [cm]");
  TH2D* Vtx_ZY_Lcan = new TH2D("Vtx_ZY_Lcan","Vtx_ZY_Lcan",500,-25,25,250,-12.5,12.5);
  Vtx_ZY_Lcan->SetXTitle("vertex Z [cm]");
  Vtx_ZY_Lcan->SetYTitle("vertex Y [cm]");
  TH2D* Vtx_XY_Lcan = new TH2D("Vtx_XY_Lcan","Vtx_XY_Lcan",250,-12.5,12.5,250,-12.5,12.5);
  Vtx_XY_Lcan->SetXTitle("vertex X [cm]");
  Vtx_XY_Lcan->SetYTitle("vertex Y [cm]");
  TH2D* Vtx_ZX_Lcan_fid = new TH2D("Vtx_ZX_Lcan_fid","Vtx_ZX_Lcan_fid",500,-25,25,250,-12.5,12.5);
  Vtx_ZX_Lcan_fid->SetXTitle("vertex Z [cm]");
  Vtx_ZX_Lcan_fid->SetYTitle("vertex X [cm]");
  TH2D* Vtx_ZY_Lcan_fid = new TH2D("Vtx_ZY_Lcan_fid","Vtx_ZY_Lcan_fid",500,-25,25,250,-12.5,12.5);
  Vtx_ZY_Lcan_fid->SetXTitle("vertex Z [cm]");
  Vtx_ZY_Lcan_fid->SetYTitle("vertex Y [cm]");
  TH2D* Vtx_XY_Lcan_fid = new TH2D("Vtx_XY_Lcan_fid","Vtx_XY_Lcan_fid",250,-12.5,12.5,250,-12.5,12.5);
  Vtx_XY_Lcan_fid->SetXTitle("vertex X [cm]");
  Vtx_XY_Lcan_fid->SetYTitle("vertex Y [cm]");


  // w/ kinematic fit
  //TH2F* MMom_MMass_fid_kin;
  //TH2F* IMppim1_IMppim2_kin;
  //TH2F* MMom_PMom_kin[];
  //TH2F* q_IMppipi_kin[2];

  const int nbinIMppipi = 60;//1.2-2.1 GeV/c^2 15 MeV.
  //const int nbinIMppipi = 160;//1.2-2.1 GeV/c^2
  const int nbinIMppipicos = 60;//1.2-2.1 GeV/c^2
  const double IMppipilow = 1.2;//1.2-2.1 GeV/c^2
  const double IMppipihigh = 2.1;//1.2-2.1 GeV/c^2
  const int nbinq = 30;//0-1.5 GeV/c 
  //const int nbinq = 100;//0-1.5 GeV/c
  const int nbinIMppi = 2000; //1-2 GeV/c^2
  const int nbinpmiss = 150; //0-1.5 GeV/c
  const double pmisslow = 0.0;
  const double pmisshigh = 1.5;

  MMom_MMass = new TH2F("MMom_MMass","MMom_MMass", nbinpmiss, pmisslow, pmisshigh, 200, 0, 2.0);
  MMom_MMass->SetXTitle("Missing Mass [GeV/c^{2}]");
  MMom_MMass->SetYTitle("Missing Mom. [GeV/c]");

  MMom_MMass_2 = new TH2F("MMom_MMass_2","MMom_MMass_2", nbinpmiss, pmisslow, pmisshigh, 200, 0, 2.0);
  MMom_MMass_2->SetXTitle("Missing Mass [GeV/c^{2}]");
  MMom_MMass_2->SetYTitle("Missing Mom. [GeV/c]");

  q_MMass = new TH2F("q_MMass","q_MMass", nbinpmiss, pmisslow, pmisshigh, 200, 0, 2.0);
  q_MMass->SetXTitle("Missing Mass [GeV/c^{2}]");
  q_MMass->SetYTitle("Mom. Transfer. [GeV/c]");

  q_MMass_forward = new TH2F("q_MMass_forward","q_MMass_forward", nbinpmiss, pmisslow, pmisshigh, 200, 0, 2.0);
  q_MMass_forward->SetXTitle("Missing Mass [GeV/c^{2}]");
  q_MMass_forward->SetYTitle("Mom. Transfer. [GeV/c]");

  MMom_MMass_p = new TH2F("MMom_MMass_p","MMom_MMass_p", nbinpmiss, pmisslow, pmisshigh, 200, 0, 2.0);
  MMom_MMass_p->SetXTitle("Missing Mass [GeV/c^{2}]");
  MMom_MMass_p->SetYTitle("Missing Mom. [GeV/c]");

  MMom_MMass_p2 = new TH2F("MMom_MMass_p2","MMom_MMass_p2", nbinpmiss, pmisslow, pmisshigh, 200, 0, 2.0);
  MMom_MMass_p2->SetXTitle("Missing Mass [GeV/c^{2}]");
  MMom_MMass_p2->SetYTitle("Missing Mom. [GeV/c]");

  MMom_MMass_p_wL = new TH2F("MMom_MMass_p_wL","MMom_MMass_p_wL", nbinpmiss, pmisslow, pmisshigh, 200, 0, 2.0);
  MMom_MMass_p_wL->SetXTitle("Missing Mass [GeV/c^{2}]");
  MMom_MMass_p_wL->SetYTitle("Missing Mom. [GeV/c]");

  MMom_MMass_p_wL_nop2 = new TH2F("MMom_MMass_p_wL_nop2","MMom_MMass_p_wL_nop2", nbinpmiss, pmisslow, pmisshigh, 200, 0, 2.0);
  MMom_MMass_p_wL_nop2->SetXTitle("Missing Mass [GeV/c^{2}]");
  MMom_MMass_p_wL_nop2->SetYTitle("Missing Mom. [GeV/c]");
  
  MMom_MMass_p_wL_wp2 = new TH2F("MMom_MMass_p_wL_wp2","MMom_MMass_p_wL_wp2", nbinpmiss, pmisslow, pmisshigh, 200, 0, 2.0);
  MMom_MMass_p_wL_wp2->SetXTitle("Missing Mass [GeV/c^{2}]");
  MMom_MMass_p_wL_wp2->SetYTitle("Missing Mom. [GeV/c]");
  
  q_PMom = new TH2F("q_PMom","q_PMom",200,0,2,200,0,2.0);
  q_PMom->SetXTitle("P Mom. [GeV/c]");
  q_PMom->SetYTitle("Mom. Transfer. [GeV/c]");

  q_MMom = new TH2F("q_MMom","q_MMom",200,0,2,200,0,2.0);
  q_MMom->SetXTitle("Miss Mom. [GeV/c]");
  q_MMom->SetYTitle("Mom. Transfer. [GeV/c]");

  q_MMomCosTheta = new TH2F("q_MMomCosTheta","q_MMomCosTheta",200,0,1,200,0,2.0);
  q_MMomCosTheta->SetXTitle("Miss Mom. CosTheta ");
  q_MMomCosTheta->SetYTitle("Mom. Transfer. [GeV/c]");

  q_MMomCosTheta_nop2 = new TH2F("q_MMomCosTheta_nop2","q_MMomCosTheta_nop2",200,0,1,200,0,2.0);
  q_MMomCosTheta_nop2->SetXTitle("Miss Mom. CosTheta ");
  q_MMomCosTheta_nop2->SetYTitle("Mom. Transfer. [GeV/c]");

  q_MMomCosTheta_wp2 = new TH2F("q_MMomCosTheta_wp2","q_MMomCosTheta_wp2",200,0,1,200,0,2.0);
  q_MMomCosTheta_wp2->SetXTitle("Miss Mom. CosTheta ");
  q_MMomCosTheta_wp2->SetYTitle("Mom. Transfer. [GeV/c]");

  q_PMom_p_wL = new TH2F("q_PMom_p_wL","q_PMom_p_wL",200,0,2,200,0,2.0);
  q_PMom_p_wL->SetXTitle("P Mom. [GeV/c]");
  q_PMom_p_wL->SetYTitle("Mom. Transfer. [GeV/c]");

  q_LMom_p_wL = new TH2F("q_LMom_p_wL","q_LMom_p_wL",200,0,2,200,0,2.0);
  q_LMom_p_wL->SetXTitle("#Lambda Mom. [GeV/c]");
  q_LMom_p_wL->SetYTitle("Mom. Transfer. [GeV/c]");

  q_MMom_p_wL = new TH2F("q_MMom_p_wL","q_MMom_p_wL",200,0,2,200,0,2.0);
  q_MMom_p_wL->SetXTitle("Miss Mom. [GeV/c]");
  q_MMom_p_wL->SetYTitle("Mom. Transfer. [GeV/c]");

  q_MMomCosTheta_p_wL = new TH2F("q_MMomCosTheta_p_wL","q_MMomCosTheta_p_wL",200,0,1,200,0,2.0);
  q_MMomCosTheta_p_wL->SetXTitle("Miss Mom. [GeV/c]");
  q_MMomCosTheta_p_wL->SetYTitle("Mom. Transfer. CosTheta ");

  q_MMomCosTheta_p_wL_nop2 = new TH2F("q_MMomCosTheta_p_wL_nop2","q_MMomCosTheta_p_wL_nop2",200,0,1,200,0,2.0);
  q_MMomCosTheta_p_wL_nop2->SetXTitle("Miss Mom. [GeV/c]");
  q_MMomCosTheta_p_wL_nop2->SetYTitle("Mom. Transfer. CosTheta ");

  q_MMomCosTheta_p_wL_wp2 = new TH2F("q_MMomCosTheta_p_wL_wp2","q_MMomCosTheta_p_wL_wp2",200,0,1,200,0,2.0);
  q_MMomCosTheta_p_wL_wp2->SetXTitle("Miss Mom. [GeV/c]");
  q_MMomCosTheta_p_wL_wp2->SetYTitle("Mom. Transfer. CosTheta ");

  q2_MMom2CosTheta_p2_wL_wp2 = new TH2F("q2_MMom2CosTheta_p2_wL_wp2","q2_MMom2CosTheta_p2_wL_wp2",200,0,1,200,0,2.0);
  q2_MMom2CosTheta_p2_wL_wp2->SetXTitle("Miss Mom. [GeV/c]");
  q2_MMom2CosTheta_p2_wL_wp2->SetYTitle("Mom. Transfer. CosTheta ");


  q_P2Mom = new TH2F("q_P2Mom","q_P2Mom",200,0,2,200,0,2.0);
  q_P2Mom->SetXTitle("P_{2} Mom. [GeV/c]");
  q_P2Mom->SetYTitle("Mom. Transfer. [GeV/c]");

  q_P2Mom_p2_wL = new TH2F("q_P2Mom_p2_wL","q_P2Mom_p2_wL",200,0,2,200,0,2.0);
  q_P2Mom_p2_wL->SetXTitle("P_{2} Mom. [GeV/c]");
  q_P2Mom_p2_wL->SetYTitle("Mom. Transfer. [GeV/c]");

  q_PMom_p_wL_sum = new TH2F("q_PMom_p_wL_sum","q_PMom_p_wL_sum",200,0,2,200,0,2.0);
  q_PMom_p_wL_sum->SetXTitle("P_{2} Mom. [GeV/c]");
  q_PMom_p_wL_sum->SetYTitle("Mom. Transfer. [GeV/c]");

  MMom_PMom = new TH2F("MMom_PMom","MMom_PMom",200,0,2,200,0,2.0);
  MMom_PMom->SetXTitle("P Mom. [GeV/c]");
  MMom_PMom->SetYTitle("Missing Mom. [GeV/c]");

  MMom_PMom_2 = new TH2F("MMom_PMom_2","MMom_PMom_2",200,0,2,200,0,2.0);
  MMom_PMom_2->SetXTitle("P Mom. [GeV/c]");
  MMom_PMom_2->SetYTitle("Missing Mom. [GeV/c]");

  IMppim1_IMppim2 = new TH2F("IMppim1_IMppim2","IMppim1_IMppim2",nbinIMppi,1.,2.0,nbinIMppi,1.,2.0);
  IMppim1_IMppim2->SetXTitle("IM(p#pi^{-}2) [GeV/c^{2}]");
  IMppim1_IMppim2->SetYTitle("IM(p#pi^{-}1) [GeV/c^{2}]");

  IMppim1_IMp2pim1 = new TH2F("IMppim1_IMp2pim1","IMppim1_IMp2pim1",nbinIMppi,1.,2.0,nbinIMppi,1.,2.0);
  IMppim1_IMp2pim1->SetXTitle("IM(p_{2}#pi^{-}1) [GeV/c^{2}]");
  IMppim1_IMp2pim1->SetYTitle("IM(p#pi^{-}1) [GeV/c^{2}]");

  IMppim2_IMp2pim2 = new TH2F("IMppim2_IMp2pim2","IMppim2_IMp2pim2",nbinIMppi,1.,2.0,nbinIMppi,1.,2.0);
  IMppim2_IMp2pim2->SetXTitle("IM(p_{2}#pi^{-}2) [GeV/c^{2}]");
  IMppim2_IMp2pim2->SetYTitle("IM(p#pi^{-}2) [GeV/c^{2}]");

  IMp2pim1_IMp2pim2 = new TH2F("IMp2pim1_IMp2pim2","IMp2pim1_IMp2pim2",nbinIMppi,1.,2.0,nbinIMppi,1.,2.0);
  IMp2pim1_IMp2pim2->SetXTitle("IM(p_{2}#pi^{-}2) [GeV/c^{2}]");
  IMp2pim1_IMp2pim2->SetYTitle("IM(p_{2}#pi^{-}1) [GeV/c^{2}]");

  IMppim1_IMppim2_p = new TH2F("IMppim1_IMppim2_p","IMppim1_IMppim2_p",nbinIMppi,1.,2.0,nbinIMppi,1.,2.0);
  IMppim1_IMppim2_p->SetXTitle("IM(p#pi^{-}2) [GeV/c^{2}]");
  IMppim1_IMppim2_p->SetYTitle("IM(p#pi^{-}1) [GeV/c^{2}]");

  IMp2pim1_IMp2pim2_p2 = new TH2F("IMp2pim1_IMp2pim2_p2","IMp2pim1_IMp2pim2_p2",nbinIMppi,1.,2.0,nbinIMppi,1.,2.0);
  IMp2pim1_IMp2pim2_p2->SetXTitle("IM(p#pi^{-}2) [GeV/c^{2}]");
  IMp2pim1_IMp2pim2_p2->SetYTitle("IM(p#pi^{-}1) [GeV/c^{2}]");

  IMppim1_IMppim2_p_wL = new TH2F("IMppim1_IMppim2_p_wL","IMppim1_IMppim2_p_wL",nbinIMppi,1.,2.0,nbinIMppi,1.,2.0);
  IMppim1_IMppim2_p_wL->SetXTitle("IM(p#pi^{-}2) [GeV/c^{2}]");
  IMppim1_IMppim2_p_wL->SetYTitle("IM(p#pi^{-}1) [GeV/c^{2}]");

  IMp2pim1_IMp2pim2_p2_wL = new TH2F("IMp2pim1_IMp2pim2_p2_wL","IMp2pim1_IMp2pim2_p2_wL",nbinIMppi,1.,2.0,nbinIMppi,1.,2.0);
  IMp2pim1_IMp2pim2_p2_wL->SetXTitle("IM(p_{2}#pi^{-}2) [GeV/c^{2}]");
  IMp2pim1_IMp2pim2_p2_wL->SetYTitle("IM(p_{2}#pi^{-}1) [GeV/c^{2}]");

  IMMissppim1_IMMissppim2_p_wL = new TH2F("IMMissppim1_IMMissppim2_p_wL","IMMissppim1_IMMissppim2_p_wL",nbinIMppi,1.,2.0,nbinIMppi,1.,2.0);
  IMMissppim1_IMMissppim2_p_wL->SetXTitle("IM(missp#pi^{-}2) [GeV/c^{2}]");
  IMMissppim1_IMMissppim2_p_wL->SetYTitle("IM(missp#pi^{-}1) [GeV/c^{2}]");



  MMass_IMppim1 = new TH2F("MMass_IMppim1","MMass_IMppim1",nbinIMppi,1.,2.0,nbinpmiss, pmisslow, pmisshigh);
  MMass_IMppim1->SetXTitle("IM(p#pi^{-}) [GeV/c^{2}]");
  MMass_IMppim1->SetYTitle("Missing Mass [GeV/c^{2}]");

  MMass_IMppim2 = new TH2F("MMass_IMppim2","MMass_IMppim2",nbinIMppi,1.,2.0,nbinpmiss, pmisslow, pmisshigh);
  MMass_IMppim2->SetXTitle("IM(p#pi^{-}) [GeV/c^{2}]");
  MMass_IMppim2->SetYTitle("Missing Mass [GeV/c^{2}]");

  MMass_IMppim_wL = new TH2F("MMass_IMppim_wL","MMass_IMppim_wL",nbinIMppi,1.,2.0,nbinpmiss, pmisslow, pmisshigh);
  MMass_IMppim_wL->SetXTitle("IM(p#pi^{-}) [GeV/c^{2}]");
  MMass_IMppim_wL->SetYTitle("Missing Mass [GeV/c^{2}]");

  MMass_IMppim_p_wL = new TH2F("MMass_IMppim_p_wL","MMass_IMppim_p_wL",nbinIMppi,1.,2.0,nbinpmiss, pmisslow, pmisshigh);
  MMass_IMppim_p_wL->SetXTitle("IM(p#pi^{-}) [GeV/c^{2}]");
  MMass_IMppim_p_wL->SetYTitle("Missing Mass [GeV/c^{2}]");

  MMass_IMp2pim_p2_wL = new TH2F("MMass_IMp2pim_p2_wL","MMass_IMp2pim_p2_wL",nbinIMppi,1.,2.0,nbinpmiss, pmisslow, pmisshigh);
  MMass_IMp2pim_p2_wL->SetXTitle("IM(p#pi^{-}) [GeV/c^{2}]");
  MMass_IMp2pim_p2_wL->SetYTitle("Missing Mass [GeV/c^{2}]");

  MMass_IMppipi_wL = new TH2F("MMass_IMppipi_wL","MMass_IMppipi_wL",nbinIMppipi,IMppipilow,IMppipihigh,nbinpmiss, pmisslow, pmisshigh);
  MMass_IMppipi_wL->SetXTitle("IM(p#pi^{-}#pi^{+}) [GeV/c^{2}]");
  MMass_IMppipi_wL->SetYTitle("Missing Mass [GeV/c^{2}]");

  MMass_IMp2pipi_wL = new TH2F("MMass_IMp2pipi_wL","MMass_IMp2pipi_wL",nbinIMppipi,IMppipilow,IMppipihigh,nbinpmiss, pmisslow, pmisshigh);
  MMass_IMp2pipi_wL->SetXTitle("IM(p#pi^{-}#pi^{+}) [GeV/c^{2}]");
  MMass_IMp2pipi_wL->SetYTitle("Missing Mass [GeV/c^{2}]");

  MMass_IMppipi_wL_sum = new TH2F("MMass_IMppipi_wL_sum","MMass_IMppipi_wL_sum",nbinIMppipi,IMppipilow,IMppipihigh,nbinpmiss, pmisslow, pmisshigh);
  MMass_IMppipi_wL_sum->SetXTitle("IM(p#pi^{-}#pi^{+}) [GeV/c^{2}]");
  MMass_IMppipi_wL_sum->SetYTitle("Missing Mass [GeV/c^{2}]");

  MMass_IMppipi_wL_sum_forward = new TH2F("MMass_IMppipi_wL_sum_forward","MMass_IMppipi_wL_sum_forward",nbinIMppipi,IMppipilow,IMppipihigh,nbinpmiss, pmisslow, pmisshigh);
  MMass_IMppipi_wL_sum_forward->SetXTitle("IM(p#pi^{-}#pi^{+}) [GeV/c^{2}]");
  MMass_IMppipi_wL_sum_forward->SetYTitle("Missing Mass [GeV/c^{2}]");

  MpimfromL_p_sum = new TH1F("MpimfromL_p_sum","MpimfromL_p_sum",nbinIMppi/5.,1.,2.0);
  MpimfromL_p_sum->SetXTitle("Miss. Mass.  d(K^{-},#pi^{-}_{#Lambda}) [GeV/c^{2}]");
  
  Mppim1_Mppim2_p_sum = new TH2F("Mppim1_Mppim2_p_sum","Mppim1_Mppim2_p_sum",nbinIMppi/5.,1.,2.0,nbinIMppi/5.,1.,2.0);
  Mppim1_Mppim2_p_sum->SetXTitle("Miss. Mass.  d(K^{-},p#pi^{-}_{2}) [GeV/c^{2}]");
  Mppim1_Mppim2_p_sum->SetYTitle("Miss. Mass.  d(K^{-},p#pi^{-}_{1}) [GeV/c^{2}]");
  
  Mppim1_Mppim2_p_wL_sum = new TH2F("Mppim1_Mppim2_p_wL_sum","Mppim1_Mppim2_p_wL_sum",nbinIMppi/5.,1.,2.0,nbinIMppi/5.,1.,2.0);
  Mppim1_Mppim2_p_wL_sum->SetXTitle("Miss. Mass.  d(K^{-},p#pi^{-}_{2}) [GeV/c^{2}]");
  Mppim1_Mppim2_p_wL_sum->SetYTitle("Miss. Mass.  d(K^{-},p#pi^{-}_{1}) [GeV/c^{2}]");
  
  MppimnotL_sum = new TH1F("MppimnotL_sum","MppimnotL_sum",nbinIMppi/5.,1.,2.0);
  MppimnotL_sum->SetXTitle("Miss. Mass.  d(K^{-},p_{CDS}#pi^{-}_{not#Lambda}) [GeV/c^{2}]");
  
  MML_MmomL_p_wL = new TH2F("MML_MmomL_p_wL","MML_MmomL_p_wL",nbinIMppi/5.,1.,2.0,200,0.,2.0);
  MML_MmomL_p_wL->SetXTitle("Miss. Mass.  d(K^{-},L) [GeV/c^{2}]");
  MML_MmomL_p_wL->SetYTitle("Miss. Mom.  d(K^{-},L) [GeV/c^{2}]");

  q_IMppim_mc = new TH2F("q_IMppim_mc","q_IMppim_mc",nbinIMppi,1.,2., nbinq,0,1.5);
  q_IMppim_mc->SetXTitle("true IM(p#pi^{-}) [GeV/c^{2}]");
  q_IMppim_mc->SetYTitle("true Mom. Transfer [GeV/c]");

  q_IMppipi_mc = new TH2F("q_IMppipi_mc","q_IMppipi_mc",nbinIMppipi,IMppipilow,IMppipihigh, nbinq,0,1.5);
  q_IMppipi_mc->SetXTitle("true IM(p#pi^{-}#pi^{-}) [GeV/c^{2}]");
  q_IMppipi_mc->SetYTitle("true Mom. Transfer [GeV/c]");
  
  q_IMppipi_p_wL_mc = new TH2F("q_IMppipi_p_wL_mc","q_IMppipi_p_wL_mc",nbinIMppipi,IMppipilow,IMppipihigh, nbinq,0,1.5);
  q_IMppipi_p_wL_mc->SetXTitle("true IM(#Lambda#pi^{-}) [GeV/c^{2}]");
  q_IMppipi_p_wL_mc->SetYTitle("true Mom. Transfer [GeV/c]");
  
  q_IMppipi_p_wL_wp2_mc = new TH2F("q_IMppipi_p_wL_wp2_mc","q_IMppipi_p_wL_wp2_mc",nbinIMppipi,IMppipilow,IMppipihigh, nbinq,0,1.5);
  q_IMppipi_p_wL_wp2_mc->SetXTitle("true IM(#Lambda#pi^{-}) [GeV/c^{2}]");
  q_IMppipi_p_wL_wp2_mc->SetYTitle("true Mom. Transfer [GeV/c]");
  
  q_IMppipi_p_wL_nop2_mc = new TH2F("q_IMppipi_p_wL_nop2_mc","q_IMppipi_p_wL_nop2_mc",nbinIMppipi,IMppipilow,IMppipihigh, nbinq,0,1.5);
  q_IMppipi_p_wL_nop2_mc->SetXTitle("true IM(#Lambda#pi^{-}) [GeV/c^{2}]");
  q_IMppipi_p_wL_nop2_mc->SetYTitle("true Mom. Transfer [GeV/c]");

  q_IMppipi_p = new TH2F("q_IMppipi_p","q_IMppipi_p",nbinIMppipi,IMppipilow,IMppipihigh, nbinq,0,1.5);
  q_IMppipi_p->SetXTitle("IM(p#pi^{-}#pi^{-}) [GeV/c^{2}]");
  q_IMppipi_p->SetYTitle("Mom. Transfer [GeV/c]");


  q_IMp2pipi_p2 = new TH2F("q_IMp2pipi_p2","q_IMp2pipi_p2",nbinIMppipi,IMppipilow,IMppipihigh, nbinq,0,1.5);
  q_IMp2pipi_p2->SetXTitle("IM(p#pi^{-}#pi^{-}) [GeV/c^{2}]");
  q_IMp2pipi_p2->SetYTitle("Mom. Transfer [GeV/c]");

  q_IMppipi_p_wL = new TH2F("q_IMppipi_p_wL","q_IMppipi_p_wL",nbinIMppipi,IMppipilow,IMppipihigh, nbinq,0,1.5);
  q_IMppipi_p_wL->SetXTitle("IM(p#pi^{-}#pi^{-}) [GeV/c^{2}]");
  q_IMppipi_p_wL->SetYTitle("Mom. Transfer [GeV/c]");

  CosTheta_IMppipi_p_wL = new TH2F("CosTheta_IMppipi_p_wL","CosTheta_IMppipi_p_wL",nbinIMppipicos,IMppipilow,IMppipihigh, 2000,-1,1);
  CosTheta_IMppipi_p_wL->SetXTitle("IM(p#pi^{-}#pi^{-}) [GeV/c^{2}]");
  CosTheta_IMppipi_p_wL->SetYTitle("miss. p CosTheta");

  CosTheta_IMppipi_p_wL_nocombi = new TH2F("CosTheta_IMppipi_p_wL_nocombi","CosTheta_IMppipi_p_wL_nocombi",nbinIMppipicos,IMppipilow,IMppipihigh, 2000,-1,1);
  CosTheta_IMppipi_p_wL_nocombi->SetXTitle("IM(p#pi^{-}#pi^{-}) [GeV/c^{2}]");
  CosTheta_IMppipi_p_wL_nocombi->SetYTitle("miss. p CosTheta");
  
  CosThetaCM_IMppipi_p_wL_nocombi = new TH2F("CosThetaCM_IMppipi_p_wL_nocombi","CosThetaCM_IMppipi_p_wL_nocombi",nbinIMppipicos,IMppipilow,IMppipihigh, 2000,-1,1);
  CosThetaCM_IMppipi_p_wL_nocombi->SetXTitle("IM(p#pi^{-}#pi^{-}) [GeV/c^{2}]");
  CosThetaCM_IMppipi_p_wL_nocombi->SetYTitle("miss. p CosTheta");
  

  CosTheta_IMppipi_p_wL_nocombi_mc = new TH2F("CosTheta_IMppipi_p_wL_nocombi_mc","CosTheta_IMppipi_p_wL_nocombi_mc",nbinIMppipicos,IMppipilow,IMppipihigh, 2000,-1,1);
  CosTheta_IMppipi_p_wL_nocombi_mc->SetXTitle("IM(p#pi^{-}#pi^{-}) [GeV/c^{2}]");
  CosTheta_IMppipi_p_wL_nocombi_mc->SetYTitle("miss. p CosTheta");

  q_IMp2pipi_p2_wL = new TH2F("q_IMp2pipi_p2_wL","q_IMp2pipi_p2_wL",nbinIMppipi,IMppipilow,IMppipihigh, nbinq,0,1.5);
  q_IMp2pipi_p2_wL->SetXTitle("IM(p#pi^{-}#pi^{-}) [GeV/c^{2}]");
  q_IMp2pipi_p2_wL->SetYTitle("Mom. Transfer [GeV/c]");

  CosTheta_IMp2pipi_p2_wL = new TH2F("CosTheta_IMp2pipi_p2_wL","CosTheta_IMp2pipi_p2_wL",nbinIMppipicos,IMppipilow,IMppipihigh, 2000,-1,1);
  CosTheta_IMp2pipi_p2_wL->SetXTitle("IM(p#pi^{-}#pi^{-}) [GeV/c^{2}]");
  CosTheta_IMp2pipi_p2_wL->SetYTitle("miss. p CosTheta");
  
  for(int icut=0;icut<ncut;icut++){
    q_IMppipi_p_wL_sum[icut] = new TH2F(Form("q_IMppipi_p_wL_sum%d",icut),Form("q_IMppipi_p_wL_sum%d",icut),nbinIMppipi,IMppipilow,IMppipihigh, nbinq,0,1.5);
    q_IMppipi_p_wL_sum[icut]->SetXTitle("IM(#Lambda#pi^{-}) [GeV/c^{2}]");
    q_IMppipi_p_wL_sum[icut]->SetYTitle("Mom. Transfer [GeV/c]");
    
    q_IMppipi_p_wL_sum_nop2[icut] = new TH2F(Form("q_IMppipi_p_wL_sum_nop2%d",icut),Form("q_IMppipi_p_wL_sum_nop2%d",icut),nbinIMppipi,IMppipilow,IMppipihigh, nbinq,0,1.5);
    q_IMppipi_p_wL_sum_nop2[icut]->SetXTitle("IM(#Lambda#pi^{-}) [GeV/c^{2}]");
    q_IMppipi_p_wL_sum_nop2[icut]->SetYTitle("Mom. Transfer [GeV/c]");
    
    q_IMppipi_p_wL_sum_wp2[icut] = new TH2F(Form("q_IMppipi_p_wL_sum_wp2%d",icut),Form("q_IMppipi_p_wL_sum_wp2%d",icut),nbinIMppipi,IMppipilow,IMppipihigh, nbinq,0,1.5);
    q_IMppipi_p_wL_sum_wp2[icut]->SetXTitle("IM(#Lambda#pi^{-}) [GeV/c^{2}]");
    q_IMppipi_p_wL_sum_wp2[icut]->SetYTitle("Mom. Transfer [GeV/c]");

    q_IMppipi_p_wL_sum_nocombi[icut] = new TH2F(Form("q_IMppipi_p_wL_sum_nocombi%d",icut),Form("q_IMppipi_p_wL_sum_nocombi%d",icut),nbinIMppipi,IMppipilow,IMppipihigh, nbinq,0,1.5);
    q_IMppipi_p_wL_sum_nocombi[icut]->SetXTitle("IM(#Lambda#pi^{-}) [GeV/c^{2}]");
    q_IMppipi_p_wL_sum_nocombi[icut]->SetYTitle("Mom. Transfer [GeV/c]");
    
    q_IMppipi_p_wL_sum_nocombi_mc[icut] = new TH2F(Form("q_IMppipi_p_wL_sum_nocombi_mc%d",icut),Form("q_IMppipi_p_wL_sum_nocombi_mc%d",icut),nbinIMppipi,IMppipilow,IMppipihigh, nbinq,0,1.5);
    q_IMppipi_p_wL_sum_nocombi_mc[icut]->SetXTitle("true IM(#Lambda#pi^{-}) [GeV/c^{2}]");
    q_IMppipi_p_wL_sum_nocombi_mc[icut]->SetYTitle("true Mom. Transfer [GeV/c]");

    q_IMppipi_p_wL_sum_nocombi_nop2[icut] 
    = new TH2F(Form("q_IMppipi_p_wL_sum_nocombi_nop2%d",icut),Form("q_IMppipi_p_wL_sum_nocombi_nop2%d",icut),nbinIMppipi,IMppipilow,IMppipihigh, nbinq,0,1.5);
    q_IMppipi_p_wL_sum_nocombi_nop2[icut]->SetXTitle("IM(#Lambda#pi^{-}) [GeV/c^{2}]");
    q_IMppipi_p_wL_sum_nocombi_nop2[icut]->SetYTitle("Mom. Transfer [GeV/c]");
    
    q_IMppipi_p_wL_sum_nocombi_nop2_mc[icut] 
    = new TH2F(Form("q_IMppipi_p_wL_sum_nocombi_nop2_mc%d",icut),Form("q_IMppipi_p_wL_sum_nocombi_nop2_mc%d",icut),nbinIMppipi,IMppipilow,IMppipihigh, nbinq,0,1.5);
    q_IMppipi_p_wL_sum_nocombi_nop2_mc[icut]->SetXTitle("true IM(#Lambda#pi^{-}) [GeV/c^{2}]");
    q_IMppipi_p_wL_sum_nocombi_nop2_mc[icut]->SetYTitle("true Mom. Transfer [GeV/c]");

    q_IMppipi_p_wL_sum_nocombi_wp2[icut]
    = new TH2F(Form("q_IMppipi_p_wL_sum_nocombi_wp2%d",icut),Form("q_IMppipi_p_wL_sum_nocombi_wp2%d",icut),nbinIMppipi,IMppipilow,IMppipihigh, nbinq,0,1.5);
    q_IMppipi_p_wL_sum_nocombi_wp2[icut]->SetXTitle("IM(#Lambda#pi^{-}) [GeV/c^{2}]");
    q_IMppipi_p_wL_sum_nocombi_wp2[icut]->SetYTitle("Mom. Transfer [GeV/c]");
    
    q_IMppipi_p_wL_sum_nocombi_wp2_mc[icut]
    = new TH2F(Form("q_IMppipi_p_wL_sum_nocombi_wp2_mc%d",icut),Form("q_IMppipi_p_wL_sum_nocombi_wp2_mc%d",icut),nbinIMppipi,IMppipilow,IMppipihigh, nbinq,0,1.5);
    q_IMppipi_p_wL_sum_nocombi_wp2_mc[icut]->SetXTitle("true IM(#Lambda#pi^{-}) [GeV/c^{2}]");
    q_IMppipi_p_wL_sum_nocombi_wp2_mc[icut]->SetYTitle("true Mom. Transfer [GeV/c]");
  }//icut

  q_IMpppipi_p_wL_sum_nop2_addptest = new TH2F(Form("q_IMpppipi_p_wL_sum_nop2_addptest"),Form("q_IMpppipi_p_wL_sum_nop2_addptest"),120,IMppipilow+1.0,IMppipihigh+1.0, nbinq,0,1.5);
  q_IMpppipi_p_wL_sum_nop2_addptest->SetXTitle("IM(#Lambda#pi^{-}) [GeV/c^{2}]");
  q_IMpppipi_p_wL_sum_nop2_addptest->SetYTitle("Mom. Transfer [GeV/c]");

  q_IMppipi_pboth_wL = new TH2F("q_IMppipi_pboth_wL","q_IMppipi_pboth_wL",nbinIMppipi,IMppipilow,IMppipihigh, nbinq,0,1.5);
  q_IMppipi_pboth_wL->SetXTitle("IM(#Lambda#pi^{-}) [GeV/c^{2}]");
  q_IMppipi_pboth_wL->SetYTitle("Mom. Transfer [GeV/c]");

  q_IMppipi_p_wL_wp2 = new TH2F("q_IMppipi_p_wL_wp2","q_IMppipi_p_wL_wp2",nbinIMppipi,IMppipilow,IMppipihigh, nbinq,0,1.5);
  q_IMppipi_p_wL_wp2->SetXTitle("IM(#Lambda#pi^{-}) [GeV/c^{2}]");
  q_IMppipi_p_wL_wp2->SetYTitle("Mom. Transfer [GeV/c]");

  CosTheta_IMppipi_p_wL_sum = new TH2F("CosTheta_IMppipi_p_wL_sum","CosTheta_IMppipi_p_wL_sum",nbinIMppipicos,IMppipilow,IMppipihigh, 2000,-1.,1.);
  CosTheta_IMppipi_p_wL_sum->SetXTitle("IM(#Lambda#pi^{-}) [GeV/c^{2}]");
  CosTheta_IMppipi_p_wL_sum->SetYTitle("miss. p CosTheta");
  
  CosThetaCM_IMppipi_p_wL_sum = new TH2F("CosThetaCM_IMppipi_p_wL_sum","CosThetaCM_IMppipi_p_wL_sum",nbinIMppipicos,IMppipilow,IMppipihigh, 2000,-1.,1.);
  CosThetaCM_IMppipi_p_wL_sum->SetXTitle("IM(#Lambda#pi^{-}) [GeV/c^{2}]");
  CosThetaCM_IMppipi_p_wL_sum->SetYTitle("miss. p CosTheta_{CM}");
  
  CosTheta_IMppipi_p_wL_sum_mc = new TH2F("CosTheta_IMppipi_p_wL_sum_mc","CosTheta_IMppipi_p_wL_sum_mc",nbinIMppipicos,IMppipilow,IMppipihigh, 2000,-1.,1.);
  CosTheta_IMppipi_p_wL_sum_mc->SetXTitle("true IM(#Lambda#pi^{-}) [GeV/c^{2}]");
  CosTheta_IMppipi_p_wL_sum_mc->SetYTitle("true miss. p CosTheta");

  q_IMppipi_p_wL_sum_forward = new TH2F("q_IMppipi_p_wL_sum_forward","q_IMppipi_p_wL_sum_forward",nbinIMppipi,IMppipilow,IMppipihigh, nbinq,0,1.5);
  q_IMppipi_p_wL_sum_forward->SetXTitle("IM(#Lambda#pi^{-}) [GeV/c^{2}]");
  q_IMppipi_p_wL_sum_forward->SetYTitle("Mom. Transfer [GeV/c]");

  q_IMppipi_p_wL_sum_forward_plus = new TH2F("q_IMppipi_p_wL_sum_forward_plus","q_IMppipi_p_wL_sum_forward_plus",nbinIMppipi,IMppipilow,IMppipihigh, nbinq,0,1.5);
  q_IMppipi_p_wL_sum_forward_plus->SetXTitle("IM(#Lambda#pi^{-}) [GeV/c^{2}]");
  q_IMppipi_p_wL_sum_forward_plus->SetYTitle("Mom. Transfer [GeV/c]");

  q_IMppipi_p_wL_sum_forward_minus = new TH2F("q_IMppipi_p_wL_sum_forward_minus","q_IMppipi_p_wL_sum_forward_minus",nbinIMppipi,IMppipilow,IMppipihigh, nbinq,0,1.5);
  q_IMppipi_p_wL_sum_forward_minus->SetXTitle("IM(#Lambda#pi^{-}) [GeV/c^{2}]");
  q_IMppipi_p_wL_sum_forward_minus->SetYTitle("Mom. Transfer [GeV/c]");

  CosTheta_IMppipi_p_wL_sum_forward = new TH2F("CosTheta_IMppipi_p_wL_sum_forward","CosTheta_IMppipi_p_wL_sum_forward",nbinIMppipicos,IMppipilow,IMppipihigh, 2000,-1,1);
  CosTheta_IMppipi_p_wL_sum_forward->SetXTitle("IM(#Lambda#pi^{-}) [GeV/c^{2}]");
  CosTheta_IMppipi_p_wL_sum_forward->SetYTitle("miss. p CosTheta");

  //just require at least 1 hist on proton counter, no proton-ID by TOF
  q_IMppipi_p_wL_sum_fp = new TH2F("q_IMppipi_p_wL_sum_fp","q_IMppipi_p_wL_sum_fp",nbinIMppipi,IMppipilow,IMppipihigh, nbinq,0,1.5);
  q_IMppipi_p_wL_sum_fp->SetXTitle("IM(#Lambda#pi^{-}) [GeV/c^{2}]");
  q_IMppipi_p_wL_sum_fp->SetYTitle("Mom. Transfer [GeV/c]");

  IMpL_p_wL_wp2 = new TH1F("IMpL_p_wL_wp2","IMpL_p_wL_wp2",300,0,3);

  IMmisspL_IMppipi_p_wL = new TH2F("IMmisspL_IMppipi_p_wL","IMmisspL_IMppipi_p_wL",nbinIMppipi,IMppipilow,IMppipihigh,300,0,3);
  IMmisspL_q_p_wL = new TH2F("IMmisspL_q_p_wL","IMmisspL_q_p_wL",nbinq,0,1.5,300,0,3);
  
  TH2F* q_MMpppipi_p_wL;//k-d -> Lambda pim proton missing mass
  q_MMpppipi_p_wL = new TH2F("q_MMpppipi_p_wL","q_MMpppipi_p_wL",200,-1,1,nbinq,0.,1.5);
  q_MMpppipi_p_wL->SetXTitle("Miss. Mass (K^{-}d#rightarrow p#Lambda#pi^{-}) [GeV/c^{2}]");
  q_MMpppipi_p_wL->SetYTitle("Mom. Transfer [GeV/c]");
  TH2F* q_MMpppipi_p_wL_mc;//k-d -> Lambda pim proton missing mass
  q_MMpppipi_p_wL_mc = new TH2F("q_MMpppipi_p_wL_mc","q_MMpppipi_p_wL_mc",200,-1,1,nbinq,0.,1.5);
  q_MMpppipi_p_wL_mc->SetXTitle("Miss. Mass (K^{-}d#rightarrow p#Lambda#pi^{-}) [GeV/c^{2}]");
  q_MMpppipi_p_wL_mc->SetYTitle("Mom. Transfer [GeV/c]");

  TH1D* pim1cos = new TH1D("pim1cos","pim1cos",10000,-1,1);
  TH1D* pim2cos = new TH1D("pim2cos","pim2cos",10000,-1,1);
  TH1D* pcos = new TH1D("pcos","pcos",10000,-1,1);
  TH1D* ptheta = new TH1D("ptheta","ptheta",1000,-1.0*TMath::Pi(),TMath::Pi());
  TH1D* p2theta = new TH1D("p2theta","p2theta",1000,-1.0*TMath::Pi(),TMath::Pi());
  TH1D* pphi = new TH1D("pphi","pphi",1000,-1.0*TMath::Pi(),TMath::Pi());
  TH1D* p2phi = new TH1D("p2phi","p2phi",1000,-1.0*TMath::Pi(),TMath::Pi());
  TH1D* p2cos = new TH1D("p2cos","p2cos",10000,-1,1);
  TH1D* ppmisscostheta = new TH1D("ppmisscostheta","ppmisscostheta",10000,-1,1);
  TH1D* pp2phi = new TH1D("pp2phi","pp2phi",100,-1.0*TMath::Pi(),TMath::Pi());
  TH1D* pp2theta = new TH1D("pp2theta","pp2theta",100,-0.5*TMath::Pi(),0.5*TMath::Pi());
  TH2D* diffpp2theta_mom = new TH2D("diffpp2theta_mom","diffpp2theta_mom",1000,-1,1,100,-0.5*TMath::Pi(),0.5*TMath::Pi());
  TH2D* diffpmissp2theta_mom = new TH2D("diffpmissp2theta_mom","diffpmissp2theta_mom",1000,-1,1,100,-0.5*TMath::Pi(),0.5*TMath::Pi());
  diffpmissp2theta_mom->SetXTitle("diff mom. pmiss - p2 [GeV/c^{2}]");
  diffpmissp2theta_mom->SetYTitle("diff theta pmiss - p2 ");
  TH2D* pmisstheta_diffp2pmissmom = new TH2D("pmisstheta_diffp2pmissmom","pmisstheta_diffp2pmissmom",1000,-1,1,100,-1,1.0);
  pmisstheta_diffp2pmissmom->SetXTitle("diff mom. pmiss - p2 [GeV/c^{2}]");
  pmisstheta_diffp2pmissmom->SetYTitle("cos theta pmiss ");

 // TH1D* diffpmissmom_mc = new TH1D("diffpmissmom_mc","diffpmissmom_mc",1000,-1,1);
 // TH1D* diffpmissmom_nop2_mc = new TH1D("diffpmissmom_nop2_mc","diffpmissmom_nop2_mc",1000,-1,1);
 // TH1D* diffpmissmom_wp2_mc = new TH1D("diffpmissmom_wp2_mc","diffpmissmom_wp2_mc",1000,-1,1);

  TH1D* pmisscos = new TH1D("pmisscos","pmisscos",10000,-1,1);
  TH1D* pmisscos_fp = new TH1D("pmisscos_fp","pmisscos_fp",10000,-1,1);
  TH1D* pmisscos_wL_p = new TH1D("pmisscos_wL_p","pmisscos_wL_p",10000,-1,1);
  TH1D* pmisscos_wL_p_fp = new TH1D("pmisscos_wL_p_fp","pmisscos_wL_p_fp",10000,-1,1);
  TH1D* p2misscos = new TH1D("p2misscos","p2misscos",10000,-1,1);
  TH2F* react_q_IMLambdaPim_1 = new TH2F("react_q_IMLambdaPim_1","react_q_IMLambdaPim_1",nbinIMppipi,IMppipilow,IMppipihigh, nbinq,0,1.5);
  TH2F* react_q_IMLambdaPim_2 = new TH2F("react_q_IMLambdaPim_2","react_q_IMLambdaPim_2",nbinIMppipi,IMppipilow,IMppipihigh, nbinq,0,1.5);
  TH2F* q_diffpmom_mc = new TH2F("q_diffpmom_mc","q_diffpmom_mc",1100,-1,10,nbinq,0,1.5);
  TH1I* hp2flag = new TH1I("hp2flag","hp2flag",2,0,2);
  TH1D* MMass_wL_or = new TH1D("MMass_wL_or","MMass_wL_or",nbinpmiss, pmisslow, pmisshigh);
  MMass_wL_or->SetXTitle("Missing Mass [GeV/c^{2}]");
  TH1D* IMppim_p_or = new TH1D("IMppim_p_or","IMppim_p_or",10000,1,2);
  IMppim_p_or->SetXTitle("IM p#pi^{-} [GeV/c^{2}]");


  TH1D* MMass_wL_1 = new TH1D("MMass_wL_1","MMass_wL_1",nbinpmiss, pmisslow, pmisshigh);
  MMass_wL_1->SetXTitle("Missing Mass [GeV/c^{2}]");
  TH1D* MMass_wL_2 = new TH1D("MMass_wL_2","MMass_wL_2",nbinpmiss, pmisslow, pmisshigh);
  MMass_wL_2->SetXTitle("Missing Mass [GeV/c^{2}]");
  TH1D* MMass_woL = new TH1D("MMass_woL","MMass_woL",nbinpmiss, pmisslow, pmisshigh);
  MMass_woL->SetXTitle("Missing Mass [GeV/c^{2}]");

  TH2D* diffpcos_pcos = new TH2D("diffpcos_pcos","diffpcos_pcos",1000,0,1,2000,-0.1,0.1);
  diffpcos_pcos->SetXTitle("Miss. P cos#theta");
  diffpcos_pcos->SetYTitle("reco. - true  Miss. P cos#theta");

  TH2D* diffpcos_pcos_wp2 = new TH2D("diffpcos_pcos_wp2","diffpcos_pcos_wp2",1000,0,1,2000,-0.1,0.1);
  diffpcos_pcos_wp2->SetXTitle("Miss. P cos#theta");
  diffpcos_pcos_wp2->SetYTitle("reco. - true  Miss. P cos#theta");
  
  TH1D* diffLmom_mc = new TH1D("diffLmom_mc","diffLmom_mc",1000,-1,1);

  TH2D* diffIMppipi_IMppipi = new TH2D("diffIMppipi_IMppipi","diffIMppipi_IMppipi",nbinIMppipi,IMppipilow,IMppipihigh,100,-0.1,0.1);
  TH2D* diffIMppipi_IMppipi_f = new TH2D("diffIMppipi_IMppipi_f","diffIMppipi_IMppipi_f",nbinIMppipi,IMppipilow,IMppipihigh,100,-0.1,0.1);
  TH2D* diffq_q = new TH2D("diffq_q","diffq_q",nbinq,0,1.5,1000,-0.1,0.1);
  diffq_q->SetXTitle("Mom. Transfer [GeV/c]");
  diffq_q->SetYTitle("reco. - true Mom. Transfer [GeV/c]");
  TH2D* diffq_IMppipi = new TH2D("diffq_IMppipi","diffq_IMppipi", nbinIMppipi,IMppipilow,IMppipihigh,500,-0.5,0.5);
  diffq_IMppipi->SetXTitle("IM(#Lambda#pi^{-}) [GeV/c^{2}]");
  diffq_IMppipi->SetYTitle("reco. - true Mom. Transfer [GeV/c]");
  TH2D* diffpmissMom_diffLpimMom = new TH2D("diffpmissMom_diffLpimMom","diffpmissMom_diffLpimMom",1000,-1.0,1.0,1000,-1.0,1.0);
  diffpmissMom_diffLpimMom->SetXTitle("reco. - true #Lambda#pi^{-} Mom. [GeV/c]");
  diffpmissMom_diffLpimMom->SetYTitle("reco. - true Miss. Mom. [GeV/c]");
  TH2D* diffpmissMom_diffpcdsMom_nop2 = new TH2D("diffpmissMom_diffpcdsMom_nop2","diffpmissMom_diffpcdsMom_nop2",1000,0.,2.0,1000,0.,2.0);
  diffpmissMom_diffpcdsMom_nop2->SetXTitle("reco. pcds - true p miss Mom. [GeV/c]");
  diffpmissMom_diffpcdsMom_nop2->SetYTitle("reco. - true Miss. Mom. [GeV/c]");
  TH2D* diffpmissMom_diffpcdsMom_wp2 = new TH2D("diffpmissMom_diffpcdsMom_wp2","diffpmissMom_diffpcdsMom_wp2",1000,0,2.0,1000,0.,2.0);
  diffpmissMom_diffpcdsMom_wp2->SetXTitle("reco. pcds - true p miss Mom. [GeV/c]");
  diffpmissMom_diffpcdsMom_wp2->SetYTitle("reco. - true Miss. Mom. [GeV/c]");

  TH2F* pcos_Lambdacos = new TH2F("pcos_Lambdacos","pcos_Lambdacos",1000,-1,1,1000,-1,1);
  pcos_Lambdacos->SetXTitle("cos. #Lambda");
  pcos_Lambdacos->SetYTitle("cos. missp");
   
  TH2F* ppimom_IMppi_notL = new TH2F("ppimom_IMppi_notL","ppimom_IMppi_notL",1000,0,2,1000,0,2);
  //TH2F* CosTheta_diffpmiss = new TH2F("CosTheta_diffpmiss","CosTheta_diffpmiss",  );
 // 
  
  Int_t nevent = tree->GetEntries();
  std::cerr<<"# of events = "<<nevent<<std::endl;

  TCanvas *cinfo = new TCanvas("cinfo","info");
  TPaveText *pt = new TPaveText(.05,.05,.95,.7);
  pt->AddText(Form("#Lambda window : %0.3f - %0.3f",anacuts::Sigmap_MIN,anacuts::Sigmap_MAX ));
  pt->AddText(Form("miss. p window : %0.3f - %0.3f",anacuts::neutron_MIN,anacuts::neutron_MAX ));
  pt->Draw();
  //------------------------//
  //--- event roop start ---//
  //------------------------//
  TGeoManager *k18br_geom = new TGeoManager("test","test");
  Double_t hadron_hall_x = 20.0*m/2.0;
  Double_t hadron_hall_y = 5.0*m/2.0;
  Double_t hadron_hall_z = 50.0*m/2.0;
  Double_t hadron_hall_org[3] = {0., 0., 0.};
  TGeoMaterial *matVacuum    = new TGeoMaterial("Vacuum",    1.008,  1,  1.e-25);
  TGeoMedium *Vacuum;
  Vacuum       = new TGeoMedium("Vacuum",       0,  matVacuum);
  TGeoVolume *hadron_hall;
  hadron_hall = k18br_geom->MakeBox("hadron_hall",
                                    Vacuum,
                                    hadron_hall_x,
                                    hadron_hall_y,
                                    hadron_hall_z);
  k18br_geom->SetTopVolume(hadron_hall);


  TGeoVolume *targetCell = k18br_geom->MakeTube("targetCell",
                           Vacuum,
                           0.0,3.43*cm,12.524*cm/2.0);
  TGeoRotation *targetCell_rot = new TGeoRotation();
  targetCell_rot->RotateX(0.2*degree);
  targetCell_rot->RotateY(0.8*degree);
  TGeoCombiTrans *targetCell_trans = new TGeoCombiTrans(
    -0.4*cm,
    -0.5*cm,
    -5.5*cm,
    targetCell_rot);
  hadron_hall->AddNode(targetCell,CID_TarCell,targetCell_trans);

  Double_t tgtsize_r = 3.0*cm;
  Double_t tgtsize_z = 12.5*cm/2.0;
  TGeoVolume *tgt = k18br_geom->MakeTube("Fiducial",
                                         Vacuum,
                                         0.0,//inner r
                                         tgtsize_r,
                                         tgtsize_z);
  Double_t tgt_pos_x =  0.0*m;
  Double_t tgt_pos_y =  0.0*m;
  Double_t tgt_pos_z =  0.0*m;
  Double_t tgt_angle =  0.0*degree;
  TGeoRotation *tgt_rot = new TGeoRotation();
  tgt_rot->RotateY(tgt_angle);
  TGeoCombiTrans *tgt_trans = new TGeoCombiTrans(tgt_pos_x,
      tgt_pos_y,
      tgt_pos_z,
      tgt_rot);
  targetCell->AddNode(tgt,CID_Fiducial, tgt_trans);
  k18br_geom->CloseGeometry();
  TGeoNode *node = k18br_geom->FindNode(0,0,0);
  std::cout << node->GetNumber() << std::endl;
  for ( Int_t i=0; i<nevent; i++ ) {
    tree->GetEvent(i);
    if(i%10000==0) std::cout << "Event# " << i << std::endl;
    TVector3 vtx_pim1 = *vtx_pim1_cdc ;
    TVector3 vtx_pim2 = *vtx_pim2_cdc ;
    // calc missing p //
    TLorentzVector LVec_p_miss = *LVec_target+*LVec_beam-*LVec_pim1-*LVec_pim2-*LVec_p;
    double pmiss_mass = LVec_p_miss.M();
    double pmiss_mom = LVec_p_miss.P();
    bool p2flag=false;
    if((*LVec_p2).P()>0.1) p2flag = true;
    hp2flag->Fill(p2flag);
    TLorentzVector LVec_p2_miss = *LVec_target+*LVec_beam-*LVec_pim1-*LVec_pim2-*LVec_p2;
    double p2miss_mass = LVec_p2_miss.M();
    double p2miss_mom = LVec_p2_miss.P();
     
    TLorentzVector LVec_allmiss = *LVec_target+*LVec_beam-*LVec_pim1-*LVec_pim2-*LVec_p-*LVec_p2;
    // calc cos(theta) of missing p //
    TVector3 boost = ((*LVec_target)+(*LVec_beam)).BoostVector();
    TLorentzVector qkn = *LVec_beam-LVec_p_miss;
    TLorentzVector qkn2 = *LVec_beam-LVec_p2_miss;
    TLorentzVector LVec_p_miss_CM = LVec_p_miss;
    TLorentzVector LVec_p2_miss_CM = LVec_p2_miss;
    TLorentzVector LVec_beam_CM = *LVec_beam;
    LVec_p_miss_CM.Boost(-boost);
    LVec_p2_miss_CM.Boost(-boost);
    LVec_beam_CM.Boost(-boost);
    double cos_p_CM = LVec_p_miss_CM.Vect().Dot(LVec_beam_CM.Vect())/(LVec_p_miss_CM.Vect().Mag()*LVec_beam_CM.Vect().Mag());
    double cos_p2_CM = LVec_p2_miss_CM.Vect().Dot(LVec_beam_CM.Vect())/(LVec_p2_miss_CM.Vect().Mag()*LVec_beam_CM.Vect().Mag());
    TLorentzVector qkn_CM = LVec_beam_CM-LVec_p_miss_CM;
    // calc pi-pi- //
    TLorentzVector LVec_pim1_pim2 = *LVec_pim1+*LVec_pim2;


    // calc pi+p //
    TLorentzVector LVec_pim1_p = *LVec_pim1+*LVec_p;
    TLorentzVector LVec_pim2_p = *LVec_pim2+*LVec_p;
    TLorentzVector LVec_pim1_p2 = *LVec_pim1+*LVec_p2;
    TLorentzVector LVec_pim2_p2 = *LVec_pim2+*LVec_p2;

    // calc missing Lambda //
    TLorentzVector LVec_pim1_p_miss = *LVec_target+*LVec_beam-*LVec_pim1-*LVec_p;
    TLorentzVector LVec_pim2_p_miss = *LVec_target+*LVec_beam-*LVec_pim2-*LVec_p;
    TLorentzVector LVec_pim1_p2_miss = *LVec_target+*LVec_beam-*LVec_pim1-*LVec_p2;
    TLorentzVector LVec_pim2_p2_miss = *LVec_target+*LVec_beam-*LVec_pim2-*LVec_p2;
    TLorentzVector LVec_Miss_p_pim1 = *LVec_target+*LVec_beam-*LVec_pim1-*LVec_p;
    TLorentzVector LVec_Miss_p_pim2 = *LVec_target+*LVec_beam-*LVec_pim2-*LVec_p;
    TLorentzVector LVec_Miss_p2_pim1 = *LVec_target+*LVec_beam-*LVec_pim1-*LVec_p2;
    TLorentzVector LVec_Miss_p2_pim2 = *LVec_target+*LVec_beam-*LVec_pim2-*LVec_p2;
    TLorentzVector LVec_Miss_pim1 = *LVec_target+*LVec_beam-*LVec_pim1;
    TLorentzVector LVec_Miss_pim2 = *LVec_target+*LVec_beam-*LVec_pim2;
    // calc pi+pi-n //
    TLorentzVector LVec_pim1_pim2_p = *LVec_pim1+*LVec_pim2+*LVec_p;
    TLorentzVector LVec_pim1_pim2_p2 = *LVec_pim1+*LVec_pim2+*LVec_p2;
    TLorentzVector LVec_pim1_pim2_p_CM = LVec_pim1_pim2_p;
    LVec_pim1_pim2_p_CM.Boost(-boost);
    double cos_X = LVec_pim1_pim2_p_CM.Vect().Dot(LVec_beam_CM.Vect())/(LVec_pim1_pim2_p_CM.Vect().Mag()*LVec_beam_CM.Vect().Mag());

    TLorentzVector LVec_Lambda_mc ;
    TLorentzVector LVec_Lambda_reactcor;
    TLorentzVector LVec_LambdaPim_mc ;
    TLorentzVector LVec_LambdaPim_reactcor ;
    TLorentzVector qkn_mc ; 
    TLorentzVector LVec_pmiss_reactcor;
    
    if(SimMode){
      LVec_Lambda_mc = *mcmom_p+*mcmom_pim2;
      LVec_Lambda_reactcor = (*react_Lambda)*0.001;
      LVec_LambdaPim_mc = *mcmom_p+*mcmom_pim1+*mcmom_pim2;
      double m = (*react_pmiss).M()/1000.;
      TVector3 vec = (*react_pmiss).Vect()*0.001;
      LVec_pmiss_reactcor.SetVectM(vec,m);
      qkn_mc = LVec_pmiss_reactcor -*mcmom_beam;     
      TLorentzVector LVec_LambdaPim_react ;
      LVec_LambdaPim_react = *react_Lambda+*react_pim;
      LVec_LambdaPim_react = *react_Lambda+*react_pim;
      double mm = LVec_LambdaPim_react.M()/1000.;
      TVector3 vvec = LVec_LambdaPim_react.Vect()*0.001;
      LVec_LambdaPim_reactcor.SetVectM(vvec,mm);
    }

    if( (qkn.P()>=anacuts::qvalcut) && (qvalcutflag==1) ) continue;
    if( (qkn.P()<anacuts::qvalcut) && (qvalcutflag==2) ) continue;
    //double chi2 = kfSpmode_chi2<kfSmmode_chi2 ? kfSpmode_chi2:kfSmmode_chi2;
    //double pvalue = kfSmmode_pvalue<kfSpmode_pvalue ? kfSpmode_pvalue:kfSmmode_pvalue;

    pim1cos->Fill((*LVec_pim1).CosTheta());
    pim2cos->Fill((*LVec_pim2).CosTheta());
    pcos->Fill((*LVec_p).CosTheta());
    ptheta->Fill((*LVec_p).Theta());
    pphi->Fill((*LVec_p).Phi());
    pmisscos->Fill(LVec_p_miss.CosTheta());
    if(ForwardCharge)pmisscos_fp->Fill(LVec_p_miss.CosTheta());
    ppmisscostheta->Fill(cos(LVec_p_miss.Theta()-(*LVec_p).Theta()));
    if(p2flag) {
      p2cos->Fill((*LVec_p2).CosTheta());
      p2theta->Fill((*LVec_p2).Theta());
      p2phi->Fill((*LVec_p2).Phi());
      p2misscos->Fill(LVec_p2_miss.CosTheta());
      pp2phi->Fill((*LVec_p2).Phi()-(*LVec_p).Phi());
      pp2theta->Fill((*LVec_p2).Theta()-(*LVec_p).Theta());
      diffpp2theta_mom->Fill((*LVec_p2).P()-(*LVec_p).P(),(*LVec_p2).Theta()-(*LVec_p).Theta());
      diffpmissp2theta_mom->Fill((LVec_p_miss).P()-(*LVec_p2).P(),(LVec_p_miss).Theta()-(*LVec_p2).Theta());
      pmisstheta_diffp2pmissmom->Fill((LVec_p_miss).P()-(*LVec_p2).P(),(LVec_p_miss).CosTheta());
    }
    BeamMom->Fill((*LVec_beam).P());
    BeamMomCosTheta->Fill((*LVec_beam).CosTheta());
    
    //0:default, 1 half low, 2 half high, 3 sigma0 region, 4 wide range ,5 +/- 1 sigma
    bool MissPFlag[ncut]={false,false,false,false,false,false};
    bool MissP2Flag[ncut]={false,false,false,false,false,false};
    bool LambdaFlag=false;
    bool LambdaFlag_1=false;
    bool LambdaFlag_2=false;
    bool Lambda2Flag=false;
    bool Lambda2Flag_1=false;
    bool Lambda2Flag_2=false;
    bool MissLFlag=false;
    bool MissL2Flag=false;
    //-- neutron-ID, K0 and missing neutron selection --//

    double dca_pim1_beam = ((*vtx_pim1_beam)-(*vtx_pim2_cdc)).Mag();
    double dca_pim2_beam = ((*vtx_pim2_beam)-(*vtx_pim2_cdc)).Mag();
    double dca_pim1_pim2 =((*CA_pim1_pim1pim2)-(*CA_pim2_pim1pim2)).Mag();

    //Lambda production in CDS
    // if( (anacuts::Lambda_MIN<LVec_pim1_p.M() && LVec_pim1_p.M()<anacuts::Lambda_MAX)){
    if( (anacuts::Lambda_MIN_narrow<LVec_pim1_p.M() && LVec_pim1_p.M()<anacuts::Lambda_MAX_narrow)) {
      LambdaFlag=true;
      LambdaFlag_1=true;
    }
//    if( (anacuts::Lambda_MIN<LVec_pim2_p.M() && LVec_pim2_p.M()<anacuts::Lambda_MAX)){
    if( (anacuts::Lambda_MIN_narrow<LVec_pim2_p.M() && LVec_pim2_p.M()<anacuts::Lambda_MAX_narrow)) {
      LambdaFlag=true;
      LambdaFlag_2=true;
    }
//    if(p2flag && (anacuts::Lambda_MIN<LVec_pim1_p2.M() && LVec_pim1_p2.M()<anacuts::Lambda_MAX)){
    if(p2flag && (anacuts::Lambda_MIN_narrow<LVec_pim1_p2.M() && LVec_pim1_p2.M()<anacuts::Lambda_MAX_narrow)) {
      Lambda2Flag=true;
      Lambda2Flag_1=true;
    }
//    if(p2flag && (anacuts::Lambda_MIN<LVec_pim2_p2.M() && LVec_pim2_p2.M()<anacuts::Lambda_MAX)){
    if(p2flag && (anacuts::Lambda_MIN_narrow<LVec_pim2_p2.M() && LVec_pim2_p2.M()<anacuts::Lambda_MAX_narrow)) {
      Lambda2Flag=true;
      Lambda2Flag_2=true;
    }

    if(anacuts::Proton_MIN_narrow<pmiss_mass && pmiss_mass<anacuts::Proton_MAX_narrow ) MissPFlag[0]=true;
    if(p2flag && (anacuts::Proton_MIN_narrow<p2miss_mass) && p2miss_mass<anacuts::Proton_MAX_narrow ) MissP2Flag[0]=true;
    if(anacuts::Proton_MIN<pmiss_mass && pmiss_mass<anacuts::Proton_MAX ) MissPFlag[4]=true;
    if(p2flag && (anacuts::Proton_MIN<p2miss_mass) && p2miss_mass<anacuts::Proton_MAX ) MissP2Flag[4]=true;
    //study for S0pi-proton reaction
    //half low
    if(anacuts::Proton_MIN_narrow<pmiss_mass && pmiss_mass<pMass ) MissPFlag[1]=true;
    if(p2flag && (anacuts::Proton_MIN_narrow<p2miss_mass) && p2miss_mass<pMass ) MissP2Flag[1]=true;
    //half high
    if(pMass <= pmiss_mass && pmiss_mass<anacuts::Proton_MAX_narrow ) MissPFlag[2]=true;
    if(p2flag && (pMass<=p2miss_mass) && p2miss_mass<anacuts::Proton_MAX_narrow ) MissP2Flag[2]=true;
    //Sigma0 region
    if(anacuts::Proton_MAX_narrow < pmiss_mass && pmiss_mass < 1.05) MissPFlag[3]=true;
    if(anacuts::Proton_MAX_narrow < p2miss_mass && p2miss_mass < 1.05) MissP2Flag[3]=true;
    //+/-1 sigma cut 
    if(anacuts::Proton_MIN_tiny<pmiss_mass && pmiss_mass<anacuts::Proton_MAX_tiny ) MissPFlag[5]=true;
    if(p2flag && anacuts::Proton_MIN_tiny<p2miss_mass && p2miss_mass<anacuts::Proton_MAX_tiny ) MissP2Flag[5]=true;


    if(anacuts::Lambda_MIN-0.02<LVec_Miss_p_pim1.M() && LVec_Miss_p_pim1.M()<anacuts::Lambda_MAX+0.02) MissLFlag=true;
    if(anacuts::Lambda_MIN-0.02<LVec_Miss_p_pim2.M() && LVec_Miss_p_pim2.M()<anacuts::Lambda_MAX+0.02) MissLFlag=true;
    if(anacuts::Lambda_MIN-0.02<LVec_Miss_p2_pim1.M() && LVec_Miss_p2_pim1.M()<anacuts::Lambda_MAX+0.02) MissL2Flag=true;
    if(anacuts::Lambda_MIN-0.02<LVec_Miss_p2_pim2.M() && LVec_Miss_p2_pim2.M()<anacuts::Lambda_MAX+0.02) MissL2Flag=true;

    if(p2flag) {
      IMp2pim1_IMp2pim2->Fill(LVec_pim2_p2.M(),LVec_pim1_p2.M());
    }
    bool LvtxOK=false;
    if(GetID((*vtx_pim1_beam),k18br_geom)!=CID_Fiducial) continue;
    if(GetID((*vtx_pim2_beam),k18br_geom)!=CID_Fiducial) continue;
    if(GetID((*vtx_p_beam),k18br_geom)!=CID_Fiducial) continue;
    //if(MissLFlag) continue;
    //if(MissL2Flag) continue;

    if(LambdaFlag_1) {
      Vtx_ZX_Lcan->Fill((*vtx_Lcan_p_pim1).Z(),(*vtx_Lcan_p_pim1).X());
      Vtx_ZY_Lcan->Fill((*vtx_Lcan_p_pim1).Z(),(*vtx_Lcan_p_pim1).Y());
      Vtx_XY_Lcan->Fill((*vtx_Lcan_p_pim1).X(),(*vtx_Lcan_p_pim1).Y());
      if(GetID((*vtx_Lcan_p_pim1),k18br_geom)==CID_Fiducial) {
        Vtx_ZX_Lcan_fid->Fill((*vtx_Lcan_p_pim1).Z(),(*vtx_Lcan_p_pim1).X());
        Vtx_ZY_Lcan_fid->Fill((*vtx_Lcan_p_pim1).Z(),(*vtx_Lcan_p_pim1).Y());
        Vtx_XY_Lcan_fid->Fill((*vtx_Lcan_p_pim1).X(),(*vtx_Lcan_p_pim1).Y());
        LvtxOK=true;
      }
    } else if(LambdaFlag_2) {
      Vtx_ZX_Lcan->Fill((*vtx_Lcan_p_pim2).Z(),(*vtx_Lcan_p_pim2).X());
      Vtx_ZY_Lcan->Fill((*vtx_Lcan_p_pim2).Z(),(*vtx_Lcan_p_pim2).Y());
      Vtx_XY_Lcan->Fill((*vtx_Lcan_p_pim2).X(),(*vtx_Lcan_p_pim2).Y());
      if(GetID((*vtx_Lcan_p_pim2),k18br_geom)==CID_Fiducial) {
        Vtx_ZX_Lcan_fid->Fill((*vtx_Lcan_p_pim2).Z(),(*vtx_Lcan_p_pim2).X());
        Vtx_ZY_Lcan_fid->Fill((*vtx_Lcan_p_pim2).Z(),(*vtx_Lcan_p_pim2).Y());
        Vtx_XY_Lcan_fid->Fill((*vtx_Lcan_p_pim2).X(),(*vtx_Lcan_p_pim2).Y());
        LvtxOK=true;
      }
    } else if(Lambda2Flag_1) {
      Vtx_ZX_Lcan->Fill((*vtx_Lcan_p2_pim1).Z(),(*vtx_Lcan_p2_pim1).X());
      Vtx_ZY_Lcan->Fill((*vtx_Lcan_p2_pim1).Z(),(*vtx_Lcan_p2_pim1).Y());
      Vtx_XY_Lcan->Fill((*vtx_Lcan_p2_pim1).X(),(*vtx_Lcan_p2_pim1).Y());
      if(GetID((*vtx_Lcan_p2_pim1),k18br_geom)==CID_Fiducial) {
        Vtx_ZX_Lcan_fid->Fill((*vtx_Lcan_p2_pim1).Z(),(*vtx_Lcan_p2_pim1).X());
        Vtx_ZY_Lcan_fid->Fill((*vtx_Lcan_p2_pim1).Z(),(*vtx_Lcan_p2_pim1).Y());
        Vtx_XY_Lcan_fid->Fill((*vtx_Lcan_p2_pim1).X(),(*vtx_Lcan_p2_pim1).Y());
        LvtxOK=true;
      }
    } else if(Lambda2Flag_2) {
      Vtx_ZX_Lcan->Fill((*vtx_Lcan_p2_pim2).Z(),(*vtx_Lcan_p2_pim2).X());
      Vtx_ZY_Lcan->Fill((*vtx_Lcan_p2_pim2).Z(),(*vtx_Lcan_p2_pim2).Y());
      Vtx_XY_Lcan->Fill((*vtx_Lcan_p2_pim2).X(),(*vtx_Lcan_p2_pim2).Y());
      if(GetID((*vtx_Lcan_p2_pim2),k18br_geom)==CID_Fiducial) {
        Vtx_ZX_Lcan_fid->Fill((*vtx_Lcan_p2_pim2).Z(),(*vtx_Lcan_p2_pim2).X());
        Vtx_ZY_Lcan_fid->Fill((*vtx_Lcan_p2_pim2).Z(),(*vtx_Lcan_p2_pim2).Y());
        Vtx_XY_Lcan_fid->Fill((*vtx_Lcan_p2_pim2).X(),(*vtx_Lcan_p2_pim2).Y());
        LvtxOK=true;
      }
    }

    if(SimMode) {
      if(p2flag)q_IMppipi_p_wL_wp2->Fill(LVec_pim1_pim2_p.M(),qkn.P());
      q_IMppim_mc->Fill((*react_Lambda).M(),qkn_mc.P());
      q_IMppipi_mc->Fill(LVec_LambdaPim_reactcor.M(),qkn_mc.P());
      
      /*
      if(LambdaFlag) {
        //TLorentzVector TL_p_diff = LVec_p_miss -*mcmom_pmiss;
        TLorentzVector TL_p_diff = LVec_p_miss -LVec_pmiss_reactcor;
        //diffpmissmom_mc->Fill(TL_p_diff.P());;
        if(p2flag) {
          //TLorentzVector TL_p2_diff = *LVec_p2 - *mcmom_pmiss;
          TLorentzVector TL_p2_diff = *LVec_p2 - LVec_pmiss_reactcor;
          //diffpmissmom_wp2_mc->Fill(TL_p2_diff.P());;
        }else{
          //diffpmissmom_nop2_mc->Fill(TL_p_diff.P());;
        }
      }*/
    }

    MMom_MMass->Fill(pmiss_mass,pmiss_mom);
    q_PMom->Fill((*LVec_p).P(),qkn.P());
    q_MMom->Fill(pmiss_mom,qkn.P());
    q_MMomCosTheta->Fill(LVec_p_miss.CosTheta(),qkn.P());
    if(!p2flag)q_MMomCosTheta_nop2->Fill(LVec_p_miss.CosTheta(),qkn.P());
    else if(p2flag)q_MMomCosTheta_wp2->Fill(LVec_p_miss.CosTheta(),qkn.P());
    MMom_PMom->Fill((*LVec_p).P(),pmiss_mom);
    if(p2flag) {
      q_P2Mom->Fill((*LVec_p2).P(),qkn2.P());
      MMom_MMass_2->Fill(p2miss_mass,p2miss_mom);
      MMom_PMom_2->Fill((*LVec_p2).P(),p2miss_mom);
    }
    IMppim1_IMppim2->Fill(LVec_pim2_p.M(),LVec_pim1_p.M());
    if(p2flag) {
      IMppim1_IMp2pim1->Fill(LVec_pim1_p2.M(),LVec_pim1_p.M());
      IMppim2_IMp2pim2->Fill(LVec_pim2_p2.M(),LVec_pim2_p.M());
    }
    MMass_IMppim1->Fill(LVec_pim1_p.M(),pmiss_mass);
    MMass_IMppim2->Fill(LVec_pim2_p.M(),pmiss_mass);

    if(!LvtxOK) continue;
    if(MissPFlag[0]) {
      MMom_MMass_p->Fill(pmiss_mass,pmiss_mom);
      IMppim1_IMppim2_p->Fill(LVec_pim2_p.M(),LVec_pim1_p.M());
      q_IMppipi_p->Fill(LVec_pim1_pim2_p.M(),qkn.P());
      Mppim1_Mppim2_p_sum->Fill(LVec_Miss_p_pim2.M(),LVec_Miss_p_pim1.M() );
      if(SimMode) {
        if(!p2flag)react_q_IMLambdaPim_1->Fill(LVec_LambdaPim_reactcor.M(),qkn_mc.P());

        TLorentzVector TL_p_diff = *LVec_p - LVec_pmiss_reactcor;
        q_diffpmom_mc->Fill(TL_p_diff.P(),qkn_mc.P());
      }
    }

    if(MissP2Flag[0]) {
      MMom_MMass_p2->Fill(p2miss_mass,p2miss_mom);
      IMp2pim1_IMp2pim2_p2->Fill(LVec_pim2_p.M(),LVec_pim1_p.M());
      q_IMp2pipi_p2->Fill(LVec_pim1_pim2_p2.M(),qkn2.P());
      Mppim1_Mppim2_p_sum->Fill(LVec_Miss_p_pim2.M(),LVec_Miss_p_pim1.M() );
      if(SimMode) {
        react_q_IMLambdaPim_2->Fill(LVec_LambdaPim_reactcor.M(),qkn_mc.P());
      }
    }

    if(LambdaFlag) {
      if((anacuts::Lambda_MIN<LVec_pim1_p.M() && LVec_pim1_p.M()<anacuts::Lambda_MAX)) {
        MMass_IMppim_wL->Fill(LVec_pim1_p.M(),pmiss_mass);
      } else if((anacuts::Lambda_MIN<LVec_pim2_p.M() && LVec_pim2_p.M()<anacuts::Lambda_MAX)) {
        MMass_IMppim_wL->Fill(LVec_pim2_p.M(),pmiss_mass);
      }
    }

    if(LambdaFlag) {
      MMass_wL_or->Fill(LVec_p_miss.M());
      MMass_wL_1->Fill(LVec_p_miss.M());
      MMass_IMppipi_wL->Fill(LVec_pim1_pim2_p.M(),pmiss_mass);
      MMass_IMppipi_wL_sum->Fill(LVec_pim1_pim2_p.M(),pmiss_mass);
      q_MMass->Fill(pmiss_mass,qkn.P());
      if(LVec_p_miss.CosTheta()>ForwardAngle) {
        q_MMass_forward->Fill(pmiss_mass,qkn.P());
        MMass_IMppipi_wL_sum_forward->Fill(LVec_pim1_pim2_p.M(),pmiss_mass);
      }
      if(p2flag){
        if(MissPFlag[0]) {
          if(LambdaFlag_1){
            TLorentzVector LVec_ppi_notL = (*LVec_pim2)+(*LVec_p2);
            ppimom_IMppi_notL->Fill(LVec_ppi_notL.M(),LVec_ppi_notL.P());
          }else if(LambdaFlag_2){
            TLorentzVector LVec_ppi_notL = (*LVec_pim1)+(*LVec_p2);
            ppimom_IMppi_notL->Fill(LVec_ppi_notL.M(),LVec_ppi_notL.P());
          }
        }
      }
    } else if(Lambda2Flag) {
      MMass_wL_or->Fill(LVec_p2_miss.M());
      MMass_wL_2->Fill(LVec_p2_miss.M());
      MMass_IMp2pipi_wL->Fill(LVec_pim1_pim2_p2.M(),p2miss_mass);
      MMass_IMppipi_wL_sum->Fill(LVec_pim1_pim2_p2.M(),p2miss_mass);
      q_MMass->Fill(p2miss_mass,qkn2.P());
      if(LVec_p_miss.CosTheta()>ForwardAngle) {
        q_MMass_forward->Fill(p2miss_mass,qkn2.P());
        MMass_IMppipi_wL_sum_forward->Fill(LVec_pim1_pim2_p2.M(),p2miss_mass);
      }
      if(p2flag){
        if(MissP2Flag[0]) {
          if(Lambda2Flag_1){
            TLorentzVector LVec_ppi_notL = (*LVec_pim2)+(*LVec_p);
            ppimom_IMppi_notL->Fill(LVec_ppi_notL.M(),LVec_ppi_notL.P());
          }else if(Lambda2Flag_2){
            TLorentzVector LVec_ppi_notL = (*LVec_pim1)+(*LVec_p);
            ppimom_IMppi_notL->Fill(LVec_ppi_notL.M(),LVec_ppi_notL.P());
          }
        }
      }
    } else {
      MMass_woL->Fill(LVec_p_miss.M());
    }


    //for acceptance map without combinatorial background
    if(SimMode){
      double diffLpimMom = LVec_pim1_pim2_p.P()-LVec_LambdaPim_reactcor.P();
      double diffLpim2Mom = LVec_pim1_pim2_p2.P()-LVec_LambdaPim_reactcor.P();
      TLorentzVector LVec_diffpmiss = LVec_p_miss-LVec_pmiss_reactcor;
      TLorentzVector LVec_diffp2miss = LVec_p2_miss-LVec_pmiss_reactcor;
      TLorentzVector LVec_diffpcdsMom = *LVec_p-LVec_pmiss_reactcor;
      TLorentzVector LVec_diffp2cdsMom = *LVec_p2-LVec_pmiss_reactcor;
      TLorentzVector LVec_diffpim1p = LVec_pim1_p-LVec_Lambda_reactcor;
      TLorentzVector LVec_diffpim2p = LVec_pim2_p-LVec_Lambda_reactcor;
      TLorentzVector LVec_diffpim1p2 = LVec_pim1_p2-LVec_Lambda_reactcor;
      TLorentzVector LVec_diffpim2p2 = LVec_pim2_p2-LVec_Lambda_reactcor;
      double diffpmissMom = LVec_p_miss.P()-LVec_pmiss_reactcor.P();
      double diffp2missMom = LVec_p2_miss.P()-LVec_pmiss_reactcor.P();
      double diffpcdsMom = (*LVec_p).P()-LVec_pmiss_reactcor.P();
      double diffp2cdsMom = (*LVec_p2).P()-LVec_pmiss_reactcor.P();
      if(MissPFlag[0] && LambdaFlag) {
        if(LVec_diffpcdsMom.P()>LVec_diffpmiss.P()){
          if(p2flag)diffpmissMom_diffpcdsMom_wp2->Fill(LVec_diffpcdsMom.P(),LVec_diffpmiss.P());
          else diffpmissMom_diffpcdsMom_nop2->Fill(LVec_diffpcdsMom.P(),LVec_diffpmiss.P());
          diffpmissMom_diffLpimMom->Fill(diffLpimMom,diffpmissMom);
          q_IMppipi_p_wL_sum_nocombi[0]->Fill(LVec_pim1_pim2_p.M(),qkn.P());
          q_IMppipi_p_wL_sum_nocombi_mc[0]->Fill(LVec_LambdaPim_reactcor.M(),qkn_mc.P());
          if(p2flag&& (LVec_allmiss.M()>MallCut) ){ // && LVec_pmiss_reactcor.CosTheta()>0.60 ){
            q_IMppipi_p_wL_sum_nocombi_wp2[0]->Fill(LVec_pim1_pim2_p.M(),qkn.P());
            q_IMppipi_p_wL_sum_nocombi_wp2_mc[0]->Fill(LVec_LambdaPim_reactcor.M(),qkn_mc.P());
          }else if(!p2flag){
            q_IMppipi_p_wL_sum_nocombi_nop2[0]->Fill(LVec_pim1_pim2_p.M(),qkn.P());
            q_IMppipi_p_wL_sum_nocombi_nop2_mc[0]->Fill(LVec_LambdaPim_reactcor.M(),qkn_mc.P());
          }
          if(LVec_diffpim1p.P()<LVec_diffpim2p.P() ){
            diffLmom_mc->Fill(LVec_pim1_p.P()-LVec_Lambda_reactcor.P());
          }else{
            diffLmom_mc->Fill(LVec_pim2_p.P()-LVec_Lambda_reactcor.P());
          }
        }
      }else if(MissP2Flag[0] && Lambda2Flag) {
        if( (LVec_diffp2cdsMom.P()>LVec_diffp2miss.P()) &&  (LVec_allmiss.M()>-0.2) ){ //&& (LVec_pmiss_reactcor.CosTheta()>0.60) ){
          diffpmissMom_diffpcdsMom_wp2->Fill(LVec_diffp2cdsMom.P(),LVec_diffp2miss.P());
          diffpmissMom_diffLpimMom->Fill(diffLpim2Mom,diffp2missMom);
          q_IMppipi_p_wL_sum_nocombi[0]->Fill(LVec_pim1_pim2_p2.M(),qkn2.P());
          q_IMppipi_p_wL_sum_nocombi_mc[0]->Fill(LVec_LambdaPim_reactcor.M(),qkn_mc.P());
          if( LVec_allmiss.M()>MallCut   ){ 
            q_IMppipi_p_wL_sum_nocombi_wp2[0]->Fill(LVec_pim1_pim2_p2.M(),qkn2.P());
            q_IMppipi_p_wL_sum_nocombi_wp2_mc[0]->Fill(LVec_LambdaPim_reactcor.M(),qkn_mc.P());
          }
        }
        if(LVec_diffpim1p2.P()<LVec_diffpim2p2.P() ){
          diffLmom_mc->Fill(LVec_pim1_p2.P()-LVec_Lambda_reactcor.P());
        }else{
          diffLmom_mc->Fill(LVec_pim2_p2.P()-LVec_Lambda_reactcor.P());
        }
      }
    }

    if(MissPFlag[0] && LambdaFlag) {
      q_IMppipi_p_wL->Fill(LVec_pim1_pim2_p.M(),qkn.P());
      CosTheta_IMppipi_p_wL->Fill(LVec_pim1_pim2_p.M(),LVec_p_miss.CosTheta());
      if(SimMode){
        TLorentzVector LVec_diffpcdsMom = *LVec_p-LVec_pmiss_reactcor;
        TLorentzVector LVec_diffpmiss = LVec_p_miss-LVec_pmiss_reactcor;
        if(LVec_diffpcdsMom.P()>LVec_diffpmiss.P()){
          CosTheta_IMppipi_p_wL_nocombi->Fill(LVec_pim1_pim2_p.M(),LVec_p_miss.CosTheta());
          CosThetaCM_IMppipi_p_wL_nocombi->Fill(LVec_pim1_pim2_p.M(),cos_p_CM);
          CosTheta_IMppipi_p_wL_nocombi_mc->Fill(LVec_LambdaPim_reactcor.M(),LVec_pmiss_reactcor.CosTheta());
        }
      }
      q_IMppipi_p_wL_sum[0]->Fill(LVec_pim1_pim2_p.M(),qkn.P());
      
      if(p2flag){
        q_MMpppipi_p_wL->Fill(LVec_allmiss.M(),qkn.P());
        MMom_MMass_p_wL_wp2->Fill(pmiss_mass,pmiss_mom);
        if(SimMode){
          TLorentzVector LVec_allmiss_mc = *LVec_target+*LVec_beam-LVec_LambdaPim_reactcor-LVec_pmiss_reactcor;
          q_MMpppipi_p_wL_mc->Fill(LVec_allmiss_mc.M(),qkn_mc.P());
        }
      }

      if(p2flag &&  (LVec_allmiss.M()>MallCut)   )q_IMppipi_p_wL_sum_wp2[0]->Fill(LVec_pim1_pim2_p.M(),qkn.P());
      else if(!p2flag){
        q_IMppipi_p_wL_sum_nop2[0]->Fill(LVec_pim1_pim2_p.M(),qkn.P());
        TLorentzVector testp;
        TVector3 ptest(0,0,0);
        testp.SetVectM(ptest,pMass);
        TLorentzVector LVec_pim1_pim2_p_testp = LVec_pim1_pim2_p+testp;
        q_IMpppipi_p_wL_sum_nop2_addptest->Fill(LVec_pim1_pim2_p_testp.M(),qkn.P());
      }
      CosTheta_IMppipi_p_wL_sum->Fill(LVec_pim1_pim2_p.M(),LVec_p_miss.CosTheta());
      CosThetaCM_IMppipi_p_wL_sum->Fill(LVec_pim1_pim2_p.M(),cos_p_CM);
      CosTheta_IMppipi_p_wL_sum_mc->Fill(LVec_LambdaPim_reactcor.M(),LVec_pmiss_reactcor.CosTheta());
      if(LVec_p_miss.CosTheta()>ForwardAngle) {
        q_IMppipi_p_wL_sum_forward->Fill(LVec_pim1_pim2_p.M(),qkn.P());
        CosTheta_IMppipi_p_wL_sum_forward->Fill(LVec_pim1_pim2_p.M(),LVec_p_miss.CosTheta());
      }
      if(LVec_p_miss.CosTheta()>ForwardAngle+0.001) {
        q_IMppipi_p_wL_sum_forward_plus->Fill(LVec_pim1_pim2_p.M(),qkn.P());
      }
      if(LVec_p_miss.CosTheta()>ForwardAngle-0.001) {
        q_IMppipi_p_wL_sum_forward_minus->Fill(LVec_pim1_pim2_p.M(),qkn.P());
      }
      if(ForwardCharge)q_IMppipi_p_wL_sum_fp->Fill(LVec_pim1_pim2_p.M(),qkn.P());
      pmisscos_wL_p->Fill(LVec_p_miss.CosTheta());
      if(ForwardCharge)pmisscos_wL_p_fp->Fill(LVec_p_miss.CosTheta());
      q_PMom_p_wL->Fill((*LVec_p).P(),qkn.P());
      if(LambdaFlag_1) {
        q_LMom_p_wL->Fill((LVec_pim1_p.P()),qkn.P());
        IMppim_p_or->Fill(LVec_pim1_p.M());
      } else if(LambdaFlag_2) {
        q_LMom_p_wL->Fill((LVec_pim2_p.P()),qkn.P());
        IMppim_p_or->Fill(LVec_pim2_p.M());
      }
      q_MMom_p_wL->Fill(pmiss_mom,qkn.P());
      q_MMomCosTheta_p_wL->Fill(LVec_p_miss.CosTheta(),qkn.P());
      if(p2flag)q_MMomCosTheta_p_wL_wp2->Fill(LVec_p_miss.CosTheta(),qkn.P());
      else if(!p2flag)q_MMomCosTheta_p_wL_nop2->Fill(LVec_p_miss.CosTheta(),qkn.P());
      q_PMom_p_wL_sum->Fill((*LVec_p).P(),qkn.P());
      if((anacuts::Lambda_MIN<LVec_pim1_p.M() && LVec_pim1_p.M()<anacuts::Lambda_MAX)) {
        MMass_IMppim_p_wL->Fill(LVec_pim1_p.M(),pmiss_mass);
      } else if((anacuts::Lambda_MIN<LVec_pim2_p.M() && LVec_pim2_p.M()<anacuts::Lambda_MAX)) {
        MMass_IMppim_p_wL->Fill(LVec_pim2_p.M(),pmiss_mass);
      }
      MMom_MMass_p_wL->Fill(pmiss_mass,pmiss_mom);
      IMppim1_IMppim2_p_wL->Fill(LVec_pim2_p.M(),LVec_pim1_p.M());
      DCA_pim1_beam->Fill( dca_pim1_beam );
      DCA_pim2_beam->Fill( dca_pim2_beam );
      DCA_pim1_pim2->Fill( dca_pim1_pim2 );

      if(LambdaFlag_1) {
        pcos_Lambdacos->Fill(LVec_pim1_p.CosTheta(),LVec_p_miss.CosTheta());
        TLorentzVector LVec_misspL = LVec_p_miss+LVec_pim1_p;
        IMmisspL_IMppipi_p_wL->Fill(LVec_pim1_pim2_p.M(),LVec_misspL.M());
        IMmisspL_q_p_wL->Fill(qkn.P(),LVec_misspL.M());
        TLorentzVector LVec_L_miss = *LVec_target+*LVec_beam-*LVec_pim1-*LVec_p;
        if(!p2flag && LVec_pim1_pim2_p.M()>1.8)MML_MmomL_p_wL->Fill(LVec_L_miss.M(),LVec_L_miss.P() );
        MppimnotL_sum->Fill(LVec_Miss_p_pim2.M() );
      } else if(LambdaFlag_2) {
        pcos_Lambdacos->Fill(LVec_pim2_p.CosTheta(),LVec_p_miss.CosTheta());
        TLorentzVector LVec_misspL = LVec_p_miss+LVec_pim2_p;
        IMmisspL_IMppipi_p_wL->Fill(qkn.P(),LVec_misspL.M());
        TLorentzVector LVec_L_miss = *LVec_target+*LVec_beam-*LVec_pim2-*LVec_p;
        if(!p2flag && LVec_pim1_pim2_p.M()>1.8)MML_MmomL_p_wL->Fill(LVec_L_miss.M(),LVec_L_miss.P() );
        MppimnotL_sum->Fill(LVec_Miss_p_pim1.M() );
      }

      TLorentzVector LVec_Missp_pim1 = LVec_p_miss+*LVec_pim1;
      TLorentzVector LVec_Missp_pim2 = LVec_p_miss+*LVec_pim2;
      IMMissppim1_IMMissppim2_p_wL->Fill( LVec_Missp_pim2.M(), LVec_Missp_pim1.M());
      Mppim1_Mppim2_p_wL_sum->Fill(LVec_Miss_p_pim2.M(),LVec_Miss_p_pim1.M() );
      if(p2flag) {
        if(LambdaFlag_1) {
          TLorentzVector LVec_p2L = *LVec_p2+LVec_pim1_p;
          IMpL_p_wL_wp2->Fill(LVec_p2L.M());
        } else if(LambdaFlag_2) {
          TLorentzVector LVec_p2L = *LVec_p2+LVec_pim2_p;
          IMpL_p_wL_wp2->Fill(LVec_p2L.M());
        }
      }
      if(SimMode) {
        double diffpcos = LVec_p_miss.CosTheta() - (*react_pmiss).CosTheta();
        if(LVec_pim1_pim2_p.M()<1.6) {
          diffpcos_pcos->Fill(LVec_p_miss.CosTheta(),diffpcos);
        }
        double diffmass = LVec_pim1_pim2_p.M()-LVec_LambdaPim_reactcor.M();
        diffIMppipi_IMppipi->Fill(LVec_pim1_pim2_p.M(),diffmass);
        if(LVec_p_miss.CosTheta()>ForwardAngle) {
          diffIMppipi_IMppipi_f->Fill(LVec_pim1_pim2_p.M(),diffmass);
        }
        double diffq = qkn.P()-qkn_mc.P();
        diffq_q->Fill(qkn.P(),diffq);
        diffq_IMppipi->Fill(LVec_pim1_pim2_p.M(), diffq);
        q_IMppipi_p_wL_mc->Fill(LVec_LambdaPim_reactcor.M(),qkn_mc.P());
        if(p2flag) q_IMppipi_p_wL_wp2_mc->Fill(LVec_LambdaPim_reactcor.M(),qkn_mc.P());
        else q_IMppipi_p_wL_nop2_mc->Fill(LVec_LambdaPim_reactcor.M(),qkn_mc.P());
      }
    }else if(MissP2Flag[0] && Lambda2Flag ) {
      q_IMp2pipi_p2_wL->Fill(LVec_pim1_pim2_p2.M(),qkn2.P());
      CosTheta_IMp2pipi_p2_wL->Fill(LVec_pim1_pim2_p2.M(),LVec_p2_miss.CosTheta());
      q_IMppipi_p_wL_sum[0]->Fill(LVec_pim1_pim2_p2.M(),qkn2.P());
      CosTheta_IMppipi_p_wL_sum->Fill(LVec_pim1_pim2_p2.M(),LVec_p2_miss.CosTheta());
      CosThetaCM_IMppipi_p_wL_sum->Fill(LVec_pim1_pim2_p2.M(),cos_p2_CM);
      if(SimMode){
        TLorentzVector LVec_diffp2cdsMom = *LVec_p2-LVec_pmiss_reactcor;
        TLorentzVector LVec_diffp2miss = LVec_p2_miss-LVec_pmiss_reactcor;
        if(LVec_diffp2cdsMom.P()>LVec_diffp2miss.P()){
          CosThetaCM_IMppipi_p_wL_nocombi->Fill(LVec_pim1_pim2_p2.M(),cos_p2_CM);
        }
      }
      if(LVec_p2_miss.CosTheta()>ForwardAngle) {
        q_IMppipi_p_wL_sum_forward->Fill(LVec_pim1_pim2_p2.M(),qkn2.P());
        CosTheta_IMppipi_p_wL_sum_forward->Fill(LVec_pim1_pim2_p2.M(),LVec_p2_miss.CosTheta());
      }
      if(LVec_p_miss.CosTheta()>ForwardAngle+0.001) {
        q_IMppipi_p_wL_sum_forward_plus->Fill(LVec_pim1_pim2_p2.M(),qkn2.P());
      }
      if(LVec_p_miss.CosTheta()>ForwardAngle-0.001) {
        q_IMppipi_p_wL_sum_forward_minus->Fill(LVec_pim1_pim2_p2.M(),qkn2.P());
      }
      q_MMpppipi_p_wL->Fill(LVec_allmiss.M(),qkn2.P());
      MMom_MMass_p_wL_wp2->Fill(p2miss_mass,p2miss_mom);
      if(SimMode){
        TLorentzVector LVec_allmiss_mc = *LVec_target+*LVec_beam-LVec_LambdaPim_reactcor-LVec_pmiss_reactcor;
        q_MMpppipi_p_wL_mc->Fill(LVec_allmiss_mc.M(),qkn_mc.P());
      }

      if(ForwardCharge)q_IMppipi_p_wL_sum_fp->Fill(LVec_pim1_pim2_p2.M(),qkn2.P());
      q_P2Mom_p2_wL->Fill((*LVec_p2).P(),qkn2.P());
      q_PMom_p_wL_sum->Fill((*LVec_p2).P(),qkn2.P());
      if(LVec_allmiss.M()>MallCut){
        q2_MMom2CosTheta_p2_wL_wp2->Fill(LVec_p2_miss.CosTheta(),qkn2.P());
        q_IMppipi_p_wL_sum_wp2[0]->Fill(LVec_pim1_pim2_p.M(),qkn.P());
        if((anacuts::Lambda_MIN<LVec_pim1_p2.M() && LVec_pim1_p2.M()<anacuts::Lambda_MAX)) {
          MMass_IMp2pim_p2_wL->Fill(LVec_pim1_p2.M(),p2miss_mass);
        } else if((anacuts::Lambda_MIN<LVec_pim2_p2.M() && LVec_pim2_p2.M()<anacuts::Lambda_MAX)) {
          MMass_IMp2pim_p2_wL->Fill(LVec_pim2_p2.M(),p2miss_mass);
        }
      }
      //MMom_MMass_p_wL->Fill(pmiss_mass,pmiss_mom);
      IMp2pim1_IMp2pim2_p2_wL->Fill(LVec_pim2_p2.M(),LVec_pim1_p2.M());
      TLorentzVector LVec_Missp_pim1 = LVec_p_miss+*LVec_pim1;
      TLorentzVector LVec_Missp_pim2 = LVec_p_miss+*LVec_pim2;
      IMMissppim1_IMMissppim2_p_wL->Fill( LVec_Missp_pim2.M(), LVec_Missp_pim1.M());
      Mppim1_Mppim2_p_wL_sum->Fill(LVec_Miss_p2_pim2.M(),LVec_Miss_p2_pim1.M() );
     
      
      if(Lambda2Flag_1) {
        pcos_Lambdacos->Fill(LVec_pim1_p2.CosTheta(),LVec_p2_miss.CosTheta());
        TLorentzVector LVec_misspL2 = LVec_p_miss+LVec_pim1_p2;
        IMmisspL_IMppipi_p_wL->Fill(LVec_pim1_pim2_p2.M(), LVec_misspL2.M());
        IMmisspL_q_p_wL->Fill(qkn2.P(), LVec_misspL2.M());
        TLorentzVector LVec_pL2 = *LVec_p+LVec_pim1_p2;
        IMpL_p_wL_wp2->Fill(LVec_pL2.M());
        IMppim_p_or->Fill(LVec_pim1_p2.M());
        MppimnotL_sum->Fill(LVec_Miss_p2_pim2.M() );
      } else if(Lambda2Flag_2) {
        pcos_Lambdacos->Fill(LVec_pim2_p2.CosTheta(),LVec_p2_miss.CosTheta());
        TLorentzVector LVec_misspL2 = LVec_p_miss+LVec_pim2_p2;
        IMmisspL_IMppipi_p_wL->Fill(LVec_pim1_pim2_p2.M(), LVec_misspL2.M());
        TLorentzVector LVec_pL2 = *LVec_p+LVec_pim2_p2;
        IMpL_p_wL_wp2->Fill(LVec_pL2.M());
        IMmisspL_q_p_wL->Fill(qkn2.P(), LVec_misspL2.M());
        IMppim_p_or->Fill(LVec_pim2_p2.M());
        MppimnotL_sum->Fill(LVec_Miss_p2_pim1.M() );
      }

      if(SimMode) {
        double diffpcos = LVec_p2_miss.CosTheta()- (*react_pmiss).CosTheta();

        if(LVec_pim1_pim2_p.M()<1.6) {
          diffpcos_pcos->Fill((*react_pmiss).CosTheta(),diffpcos);
        }
        double diffmass = LVec_pim1_pim2_p.M()-LVec_LambdaPim_reactcor.M();
        diffIMppipi_IMppipi->Fill(LVec_pim1_pim2_p2.M(),diffmass);
        if(LVec_p2_miss.CosTheta()>ForwardAngle) {
          diffIMppipi_IMppipi_f->Fill(LVec_pim1_pim2_p2.M(),diffmass);
        }
        double diffq = qkn2.P()-qkn_mc.P();
        diffq_q->Fill(qkn2.P(),diffq);
        diffq_IMppipi->Fill(LVec_pim1_pim2_p.M(), diffq);
        q_IMppipi_p_wL_mc->Fill(LVec_LambdaPim_reactcor.M(),qkn_mc.P());
        q_IMppipi_p_wL_wp2_mc->Fill(LVec_LambdaPim_reactcor.M(),qkn_mc.P());
      }
    }

    if(MissPFlag[0] && MissP2Flag[0]) {
      if(LambdaFlag_1 || LambdaFlag_2) {
        q_IMppipi_pboth_wL->Fill(LVec_pim1_pim2_p.M(),qkn.P());
      } else if(Lambda2Flag_1 || Lambda2Flag_2)  {
        q_IMppipi_pboth_wL->Fill(LVec_pim1_pim2_p2.M(),qkn2.P());
      }
    }

    for(int icut=1;icut<ncut;icut++){
      if(MissPFlag[icut] && LambdaFlag) {
        q_IMppipi_p_wL_sum[icut]->Fill(LVec_pim1_pim2_p.M(),qkn.P());
        if(p2flag && LVec_allmiss.M()>MallCut )q_IMppipi_p_wL_sum_wp2[icut]->Fill(LVec_pim1_pim2_p.M(),qkn.P());
        else if(!p2flag) q_IMppipi_p_wL_sum_nop2[icut]->Fill(LVec_pim1_pim2_p.M(),qkn.P());
        if(SimMode){
          TLorentzVector LVec_diffpmiss = LVec_p_miss-LVec_pmiss_reactcor;
          TLorentzVector LVec_diffp2miss = LVec_p2_miss-LVec_pmiss_reactcor;
          TLorentzVector LVec_diffpcdsMom = *LVec_p-LVec_pmiss_reactcor;
          TLorentzVector LVec_diffp2cdsMom = *LVec_p2-LVec_pmiss_reactcor;
          if(LVec_diffpcdsMom.P()>LVec_diffpmiss.P()){
            q_IMppipi_p_wL_sum_nocombi[icut]->Fill(LVec_pim1_pim2_p.M(),qkn.P());
            if(p2flag && ((LVec_allmiss.M()>MallCut))){
              q_IMppipi_p_wL_sum_nocombi_wp2[icut]->Fill(LVec_pim1_pim2_p.M(),qkn.P());
              q_IMppipi_p_wL_sum_nocombi_wp2_mc[icut]->Fill(LVec_LambdaPim_reactcor.M(),qkn_mc.P());
            }else if(!p2flag){
              q_IMppipi_p_wL_sum_nocombi_nop2[icut]->Fill(LVec_pim1_pim2_p.M(),qkn.P());
              q_IMppipi_p_wL_sum_nocombi_nop2_mc[icut]->Fill(LVec_LambdaPim_reactcor.M(),qkn_mc.P());
            }
          }
        }
      }else if(MissP2Flag[icut] && Lambda2Flag){
        if(LVec_allmiss.M()>MallCut){
          q_IMppipi_p_wL_sum_wp2[icut]->Fill(LVec_pim1_pim2_p2.M(),qkn2.P());
        }
        
        if(SimMode){
          TLorentzVector LVec_diffpmiss = LVec_p_miss-LVec_pmiss_reactcor;
          TLorentzVector LVec_diffp2miss = LVec_p2_miss-LVec_pmiss_reactcor;
          TLorentzVector LVec_diffpcdsMom = *LVec_p-LVec_pmiss_reactcor;
          TLorentzVector LVec_diffp2cdsMom = *LVec_p2-LVec_pmiss_reactcor;
          if(LVec_diffp2cdsMom.P()>LVec_diffp2miss.P() && ((LVec_allmiss.M()>MallCut))){
            q_IMppipi_p_wL_sum_nocombi[icut]->Fill(LVec_pim1_pim2_p2.M(),qkn2.P());
            q_IMppipi_p_wL_sum_nocombi_wp2[icut]->Fill(LVec_pim1_pim2_p2.M(),qkn2.P());
            q_IMppipi_p_wL_sum_nocombi_wp2_mc[icut]->Fill(LVec_LambdaPim_reactcor.M(),qkn_mc.P());
          }
        }//SimMode
      }
    }


  }//for ievt

  /*
  const double Kp_mass = pMass + kpMass;
  TF1 *fkp = new TF1("f", "sqrt(((x*x-[0]*[0]-[1]*[1])/(2*[0]))*((x*x-[0]*[0]-[1]*[1])/(2*[0]))-[1]*[1])",Kp_mass-0.001,2);
  fkp->SetParameter(0,pMass);
  fkp->SetParameter(1,kpMass);
  //fkp->SetLineColor(4);
  fkp->SetLineWidth(4);
  fkp->SetLineStyle(4);
  fkp->SetLineColorAlpha(kPink, 0.35);
  fkp->Draw("same");
  */
  //TCanvas *ctest = new TCanvas("ctest","ctest");
  //q_IMppipi_p_wL->Draw("colz");

  f->cd();
  TIter nexthist(gDirectory->GetList());
  //TIter nexthist(gROOT->GetList());
  TH1F *h1 = NULL;
  TH1D *h1d = NULL;
  TH2F *h2 = NULL;
  TObject *obj = NULL;
  TString outname = std::string(filename);
  outname.Replace(std::string(filename).size()-5,5,"_out.root");
  if(qvalcutflag==1) outname.Replace(std::string(outname).size()-5,5,"_qlo.root");
  if(qvalcutflag==2) outname.Replace(std::string(outname).size()-5,5,"_qhi.root");
  if(qvalcutflag==3) outname.Replace(std::string(outname).size()-5,5,"_theta15.root");
  TFile *fout = new TFile(outname.Data(),"RECREATE");
  while( (obj = (TObject*)nexthist())!=NULL  ) {
    if(obj->InheritsFrom("TH1F")) {
      h1 = (TH1F*) obj;
      h1->GetXaxis()->CenterTitle();
      //h1->GetXaxis()->SetTitleSize(0.05);
      //h1->GetXaxis()->SetTitleOffset(0.80);
    }
    if(obj->InheritsFrom("TH1D")) {
      h1d = (TH1D*) obj;
      h1d->GetXaxis()->CenterTitle();
      //h1d->GetXaxis()->SetTitleSize(0.05);
      //h1d->GetXaxis()->SetTitleOffset(0.80);
    }
    if(obj->InheritsFrom("TH2")) {
      h2 = (TH2F*) obj;
      h2->GetXaxis()->CenterTitle();
      h2->GetYaxis()->CenterTitle();
      //h2->GetXaxis()->SetTitleSize(0.05);
      //h2->GetXaxis()->SetTitleOffset(0.80);
      //h2->GetYaxis()->SetTitleSize(0.05);
      //h2->GetYaxis()->SetTitleOffset(0.85);
    }
    obj->Write();
  }
  fout->Close();

  TCanvas *c = NULL;
  //TPDF *pdf = new TPDF(pdfname);
  TSeqCollection *SCol = gROOT->GetListOfCanvases();
  int size = SCol->GetSize();
  //TIter next(gROOT->GetListOfCanvases());
  TIter next(SCol);
  for(int i=0; i<size; i++) {
    //while((c= (TCanvas*)next())){
    //pdf->NewPage();
    c= (TCanvas*)next();
    c->Draw();
    c->cd();
    TPaveText *pt;
    pt = new TPaveText(.80,0.90,0.98,0.99,"NDC");
    pt->AddText("Real Data");
    pt->SetFillColor(kCyan-9);
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

  //f->cd();
  /*
  std::cout << "closing pdf " << std::endl;
  TString outname = std::string(filename);
  outname.Replace(std::string(filename).size()-5,5,"_out.root");
  if(qvalcutflag==1) outname.Replace(std::string(outname).size()-5,5,"_qlo.root");
  if(qvalcutflag==2) outname.Replace(std::string(outname).size()-5,5,"_qhi.root");
  if(qvalcutflag==3) outname.Replace(std::string(outname).size()-5,5,"_theta15.root");
  TFile *fout = new TFile(outname.Data(),"RECREATE");
  fout->Print();
  fout->cd();
  TIter nexthist2(gDirectory->GetList());
  while( (obj = (TObject*)nexthist2())!=NULL) {
    obj->Print();
    obj->Write();
  }
  fout->cd();
  fout->Close();
  */
}

int GetID(const TVector3 &pos,TGeoManager *k18br_geom)
{
  TGeoNode *node=k18br_geom->FindNode(pos.X(),pos.Y(),pos.Z());
  if(!node) return 0;
  int id=node->GetNumber();
  return id;
}


