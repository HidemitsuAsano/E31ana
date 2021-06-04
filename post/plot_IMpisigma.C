//asano memo 

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
//#include "weightfunc.h"
//#include "weightfuncGSp.h"
//#include "weightfuncGSm.h"

const double pvalcut = 0.005;
const bool gridon=true;
const bool staton=true;
const bool UseKinFit = false;
const bool UseKinFitVal = false;
const double lumi = 56.0e9; // rough value
const double Beamsurvival =  0.4; //rough value from Kawasaki ana
const double D2density = 4.83e23; //kawasaki ana
const double DAQeff = 0.77;//kawasaki ana
const double trigeff = 2.0; //pre-scale factor of CDH3 trigger
const int ngap=1;
const int nzone=2;//for wide range side band study
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

//BG flag -> handle BG region for _woSid_won histograms
//0 : all BG  = excludes missing neutron and Sigma+/- (cross cut)
//1 : BG near sigal region = addional excluded region for BG
const int  BGFlag_woSid_won=0;

const bool IsMCweighting = false;

//check GEANT4 info. to reject fake neutron events
//maybe, also forward Sigma events should be rejected ?
const bool SimRejectFake = true;

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

void plot_IMpisigma(const char* filename="", const int qvalcutflag=0)
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
  if(qvalcutflag==1) pdfname.Replace(std::string(filename).size()-5,8,"_qlo.pdf");
  if(qvalcutflag==2) pdfname.Replace(std::string(filename).size()-5,8,"_qhi.pdf");
  if(qvalcutflag==3) pdfname.Replace(std::string(filename).size()-5,8,"_theta15.pdf");
  std::cout << "pdfname: " << pdfname << std::endl;
  std::cout << std::endl;

  std::cout << "Use Kin Fit ? " << std::endl;
  if(UseKinFit) std::cout << "Yes" << std::endl;
  else             std::cout << "No"  << std::endl;

  std::cout << "Use Kin Fit Val of Mass ? " << std::endl;
  if(UseKinFitVal) std::cout << "Yes" << std::endl;
  else             std::cout << "No"  << std::endl;

  std::cout << std::endl;
  std::cout << "Sigma selection type     " << sigmacuttype << std::endl;
  std::cout << "Side band selection type " << sidebandtype << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;

  std::cout << "IsolationFlag ? " << IsolationFlag << std::endl;
  std::cout << "CDCChargeVetoFlag ? " << CDCChargeVetoFlag << std::endl;
  std::cout << "ForwardVetoFlag ? " << ForwardVetoFlag << std::endl;


  std::cout << "MC weighting ? " << std::endl;
  if(IsMCweighting) std::cout << "Yes" << std::endl;
  else              std::cout << "No"  << std::endl;
  
  std::cout << std::endl;
  std::cout << "SimRejectFake ? "  << std::endl;
  if(SimRejectFake) std::cout << "Yes" << std::endl;
  else              std::cout << "No"  << std::endl;

  bool SimSpmode = (std::string(filename).find("Sp")!= std::string::npos);
  bool SimSmmode = (std::string(filename).find("Sm")!= std::string::npos);
  bool SimK0nnmode = (std::string(filename).find("K0nn")!= std::string::npos);//K0nn 
  bool SimK0n_nsmode = (std::string(filename).find("K0n_ns")!= std::string::npos);
  bool SimK0_nntsmode = (std::string(filename).find("K0_nnts")!= std::string::npos);
  bool SimnpipiLmode = (std::string(filename).find("npipiL")!= std::string::npos);
  bool SimnS0pippimmode = (std::string(filename).find("nS0pippim")!= std::string::npos);
  bool SimSppi0mode = (std::string(filename).find("Sppimpi0")!=std::string::npos);
  bool SimSmpi0mode = (std::string(filename).find("Smpippi0")!=std::string::npos);
  bool SimSp_nsmode = (std::string(filename).find("Sppim_ns")!=std::string::npos);
  bool SimSm_nsmode = (std::string(filename).find("Smpip_ns")!=std::string::npos);
  bool SimFakemode  = (std::string(filename).find("fakepippim_")!=std::string::npos);
  bool SimFakeK0mode  = (std::string(filename).find("fakepippimK0_")!=std::string::npos);
  bool SimFakemode_gSp = (std::string(filename).find("fakemcgSp")!=std::string::npos);
  bool SimFakemode_gSm = (std::string(filename).find("fakemcgSm")!=std::string::npos);
  bool RealDatamode = (std::string(filename).find("evanaIMpisigma")!=std::string::npos);
  bool MIXmode = (std::string(filename).find("MIX")!=std::string::npos);
  //= = = = pipipnn final-sample tree = = = =//
  if(SimFakemode) std::cout << "fake pi+pi-nX mode "  << std::endl; 
  if(SimFakeK0mode) std::cout << "fake K0barnX mode "  << std::endl; 
  if(SimFakemode_gSp) std::cout << "fake for GEANT nSp mode" << std::endl;
  if(SimFakemode_gSm) std::cout << "fake for GEANT nSm mode" << std::endl;
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
  tree->SetBranchAddress( "tofpim",&tofpim);
  tree->SetBranchAddress( "tofpip",&tofpip);
  tree->SetBranchAddress( "tofn",&tofn);
  tree->SetBranchAddress( "dE", &dE );
  tree->SetBranchAddress( "neutralseg", &neutralseg );  
  tree->SetBranchAddress( "nhitOutCDC", &nhitOutCDC ); //charge veto by Outer 3 layer of 3cdc
  tree->SetBranchAddress( "ForwardCharge", &ForwardCharge);
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
  }
  if(UseKinFit) {
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
  //std::cout << __LINE__ << std::endl;
  //weight function of BG evaluation for MC
  
  //std::cout << __LINE__ << std::endl;
  f->cd();
  // w/o kinematic fit
  TH1I* NHitCDCOut;
  TH1I* IsForwardCharge;
  TH2F* CDHphi_betainv_fid;
  TH2F* CDHz_betainv_fid;
  TH2F* CDHz_nmom_fid;
  TH2F* dE_betainv_fid;//
  TH2F* dE_nmom_fid_beta;
  //TH2F* dE_nmom_fid_beta_wK0;
  TH2F* dE_MMom_fid_beta;
  TH2F* dE_MMass_fid_beta;
  TH2F* dE_MMass_fid_beta_wSid;
  //TH2F* dE_MMom_fid_beta_woK0;
  //TH2F* dE_MMass_fid_beta_woK0;
  //TH2F* dE_MMass_fid_beta_woK0_wSid;
  TH2F* MMom_MMass;
  TH2F* MMom_MMass_wSid;
  //TH2F* MMom_MMass_wSid_n;
  TH2F* MMom_MMass_woK0;
  TH2F* MMom_MMass_woK0_wSid;
  //TH2F* MMom_MMass_woK0_woSidn;
  TH2F* MMom_MMass_woK0_woSid_won;
  TH2F* MMom_MMass_wK0_woSidn_won;
  TH2F* MMom_MMass_wK0_woSid_won;
  TH1F* MMnmiss_react; 

  //
  TH2F* IMnpim_IMnpip_mc;//store mcData node
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
  TH2F* IMnpim_MMnpim_mc;//store mcData
  TH2F* IMnpim_MMnpim_mc_wSid_n_fake;//store mcData
  TH2F* IMnpim_MMnpim_mc_wSid_n_fake_pat2;//store mcData
  TH2F* IMnpim_MMnpim_mc_wSid_n_fake_pat7;//store mcData
  TH2F* IMnpim_MMnpim_mc_vtx;//store mcData 
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
  TH2F* diff2D_MMnmiss_IMnpim_reactmc;
  TH2F* diff2D_MMnmiss_IMnpim_wSid_n_reactmc;
  TH2F* diff2D_MMnmiss_IMnpip_reactmc;
  TH2F* diff2D_MMnmiss_IMnpip_wSid_n_reactmc;
  
  //GEANT mcData - reco. data matching 
  //fake1: fake event flag based on mcData - reaction data matching
  //TH2F* diff2D_MMnmiss_IMnpim_recomc_woK0_wSid_n;//store reconstructed data - mcData
  //TH2F* diff2D_MMnmiss_IMnpip_recomc_woK0_wSid_n;//store reconstructed data - mcData
  TH2F* diff2D_MMnmiss_IMnpim_recomc_woK0_wSid_n_fake1;//store reconstructed data - mcData
  TH2F* diff2D_MMnmiss_IMnpip_recomc_woK0_wSid_n_fake1;//store reconstructed data - mcData
  TH2F* diff2D_MMnmiss_IMnpim_recomc_wK0_wSid_n;//store reconstructed data - mcData
  TH2F* diff2D_MMnmiss_IMnpip_recomc_wK0_wSid_n;//store reconstructed data - mcData
  TH2F* diff2D_MMnmiss_IMnpim_recomc_wK0_wSid_n_fake1;//store reconstructed data - mcData
  TH2F* diff2D_MMnmiss_IMnpip_recomc_wK0_wSid_n_fake1;//store reconstructed data - mcData
  TH2F* diff2D_MMnmiss_IMnpim_recomc_wSid_n;//store reconstructed data - mcData
  TH2F* diff2D_MMnmiss_IMnpip_recomc_wSid_n;//store reconstructed data - mcData
  TH2F* diff2D_MMnmiss_IMnpim_recomc_wSid_n_fake1;//store reconstructed data - mcData
  TH2F* diff2D_MMnmiss_IMnpip_recomc_wSid_n_fake1;//store reconstructed data - mcData
  //TH2F* diff2D_nmom_IMnpim_recomc_woK0_wSid_n;//store reconstructed data - mcData
  //TH2F* diff2D_nmom_IMnpip_recomc_woK0_wSid_n;//store reconstructed data - mcData
  TH2F* diff2D_nmom_IMnpim_recomc_woK0_wSid_n_fake1;//store reconstructed data - mcData
  TH2F* diff2D_nmom_IMnpip_recomc_woK0_wSid_n_fake1;//store reconstructed data - mcData
  TH2F* diff2D_nmom_IMnpim_recomc_wK0_wSid_n;//store reconstructed data - mcData
  TH2F* diff2D_nmom_IMnpip_recomc_wK0_wSid_n;//store reconstructed data - mcData
  TH2F* diff2D_nmom_IMnpim_recomc_wK0_wSid_n_fake1;//store reconstructed data - mcData
  TH2F* diff2D_nmom_IMnpip_recomc_wK0_wSid_n_fake1;//store reconstructed data - mcData
  TH2F* diff2D_nmom_IMnpim_recomc_wSid_n;//store reconstructed data - mcData
  TH2F* diffMomnpim_Momnpip_recomc_wSid_n;//store reconstructed data - mcData
  TH1F* diffMMom_recomc_wSid_n;//store reconstructed data - mcData
  TH2F* diff2D_nmom_IMnpip_recomc_wSid_n;//store reconstructed data - mcData
  TH2F* diff2D_nmom_IMnpim_recomc_wSid_n_fake1;//store reconstructed data - mcData
  TH2F* diff2D_nmom_IMnpip_recomc_wSid_n_fake1;//store reconstructed data - mcData
  

  const unsigned int nwbin = 3;
  TH2F* IMnpim_IMnpip_dE;
  TH2F* IMnpim_IMnpip_dE_woK0;
  TH2F* IMnpim_IMnpip_dE_woK0_woSp_vici;//select vicinity of the signal
  TH2F* IMnpim_IMnpip_dE_woK0_woSm_vici;//select viciniti of the signal
  TH2F* IMnpim_IMnpip_dE_woK0_woSp_viciext;//select vicinity of the signal
  TH2F* IMnpim_IMnpip_dE_woK0_woSm_viciext;//select viciniti of the signal
  TH2F* IMnpim_IMnpip_dE_wK0;
  TH2F* IMnpim_IMnpip_dE_wK0_woSid_n;
  TH2F* IMnpim_IMnpip_dE_wK0_woSid_n_45rot;
  TH2F* IMnpim_IMnpip_dE_wK0_woSid_n_45rot2;
  TH2F* IMnpim_IMnpip_dE_wK0_woSid_n_45rot3;
  TH2F* IMnpim_IMnpip_dE_wK0_n_45rot3;
  TH2F* IMnpim_IMnpip_dE_wK0_woSid_n_wbin[nwbin];
  TH2F* IMnpim_IMnpip_dE_woK0_woSid_n;
  TH2F* IMnpim_IMnpip_dE_n;//
  const unsigned int nbintemplate = 100;
  //TH2F* IMnpim_IMnpip_dE_n_bin[nbintemplate];//
  //TH2F* IMnpim_IMnpip_dE_n_reg[100];// for template fit
  TH2F* IMnpim_IMnpim_mc_dE_n;//
  TH2F* IMnpim_IMnpim_mc_dE_n_vtx;//
  TH2F* IMnpim_IMnpim_mc_dE_n_vtx_pat2;//Sigma decay neutron
  TH2F* IMnpim_IMnpim_mc_dE_n_vtx_pat7;//Initial neutron
  TH2F* IMnpip_IMnpip_mc_dE_n;//
  TH2F* IMnpip_IMnpip_mc_dE_n_vtx;//
  TH2F* IMnpip_IMnpip_mc_dE_n_vtx_pat2;//
  TH2F* IMnpip_IMnpip_mc_dE_n_vtx_pat7;//
  TH2F* IMnpim_IMnpip_dE_n_fake;//for GEANT4 sim.
  TH2F* IMnpim_IMnpip_dE_wSid_n;
  //TH2F* IMnpim_IMnpip_dE_wSid_n_bin[nbintemplate];
  TH2F* IMnpim_IMnpip_dE_wK0orwSid_n;
  TH2F* IMnpim_IMnpip_dE_wK0orwSid_n_bin[nbintemplate];
  TH2F* IMnpim_IMnpip_dE_wK0orwSid_n_wbin[nwbin];
  TH2F* IMnpim_IMnpip_dE_wK0orwSid_n_woSp;
  TH2F* IMnpim_IMnpip_dE_wK0orwSid_n_woSm;
  TH2F* IMnpim_IMnpip_dE_wSid_n_Sp;
  TH2F* IMnpim_IMnpip_dE_wSid_n_Sp_bin[nbintemplate];
  TH2F* IMnpim_IMnpip_dE_wSid_n_Sp_wbin[nwbin];
  TH2F* IMnpim_IMnpip_dE_wSid_n_Sm;
  TH2F* IMnpim_IMnpip_dE_wSid_n_Sm_bin[nbintemplate];
  TH2F* IMnpim_IMnpip_dE_wSid_n_Sm_wbin[nwbin];
  TH2F* IMnpim_IMnpip_dE_wSid_n_fake;//for GEANT4 sim.
  TH2F* IMnpim_IMnpip_dE_wSid_n_fake_pat2;//for GEANT4 sim.
  TH2F* IMnpim_IMnpip_dE_wSid_n_fake_pat7;//for GEANT4 sim.
  TH2F* IMnpim_IMnpip_dE_woK0_n;//
  TH2F* IMnpim_IMnpip_dE_woK0_wSid_n;//
  TH2F* IMnpim_IMnpip_dE_woK0_wSid_n_woSm;//
  TH2F* IMnpim_IMnpip_dE_woK0_wSid_n_woSm_wbin[nwbin];//
  TH2F* IMnpim_IMnpip_dE_woK0_wSid_n_woSp;//
  TH2F* IMnpim_IMnpip_dE_woK0_wSid_n_woSp_wbin[nwbin];//
  //TH2F* IMnpim_IMnpip_dE_woK0_wSid_n_bin[nbintemplate];//
  TH2F* IMnpim_IMnpip_dE_wK0_wSid_n;//
  TH2F* IMnpim_IMnpip_dE_wK0_wSid_n_Sp;//
  TH2F* IMnpim_IMnpip_dE_wK0_wSid_n_Sp_wbin[nwbin];//resolution based binning
  TH2F* IMnpim_IMnpip_dE_wK0_wSid_n_Sm;//
  TH2F* IMnpim_IMnpip_dE_wK0_wSid_n_Sm_wbin[nwbin];//resolution based binning
  TH2F* IMnpim_IMnpip_dE_wK0_wSid_n_SpSm;//
  TH2F* IMnpim_IMnpip_dE_woSid;//
  //TH2F* IMnpim_IMnpip_dE_woSid_won;//
  TH2F* IMnpim_IMnpip_dE_woK0_woSid;//
  TH2F* IMnpim_IMnpip_dE_woK0_woSid_won;//
  TH2F* IMnpim_IMnpip_dE_wK0_n;//
  TH2F* IMnpim_IMnpip_dE_wK0_woSid_won;//
  TH2F* IMnpim_IMnpip_dE_wK0_woSidn_won;//
  //TH2F* IMnpip_CDHphi_dE_woK0_n;//
  //TH2F* IMnpip_CDHz_dE_woK0_n;//
  //TH2F* IMnpim_CDHphi_dE_woK0_n;//
  //TH2F* IMnpim_CDHz_dE_woK0_n;//
  //TH2F* pipmom_IMpippim;//no nuetron ID
  //TH2F* pipmom_IMpippim_dE;//+nuetron ID by CDS
  //TH2F* pipmom_IMpippim_dE_n;//missing nuetron ID
  //TH2F* pimmom_IMpippim;//no nuetron ID
  //TH2F* pimmom_IMpippim_dE;//+nuetron ID by CDS
  //TH2F* pimmom_IMpippim_dE_n;//missing nuetron ID
  //TH2F* pipmom_IMnpip;//no nuetron ID by CDS
  //TH2F* pipmom_IMnpip_dE;//+nuetron ID by CDS
  //TH2F* pipmom_IMnpip_dE_n;//missing nuetron ID
  //TH2F* pipmom_MMnmiss_dE;//missing nuetron ID
  //TH2F* pipmom_MMnmiss_dE_woK0_woSidn;//missing nuetron ID
  //TH2F* pipmom_MMnmiss_dE_woK0_woSid_won;//missing nuetron ID
  //TH2F* pipmom_MMnmiss_dE_wK0_woSid_won;//missing nuetron ID
  //TH2F* pipmom_MMnmiss_dE_wK0_woSidn_won;//missing nuetron ID
  //TH2F* pipmom_pimmom_dE_woK0_woSidn;
  //TH2F* pipmom_pimmom_dE_woK0_woSid_won;
  //TH2F* pipmom_pimmom_dE_wK0_woSid_won;
  //TH2F* pimmom_IMnpim;//no nuetron ID
  //TH2F* pimmom_IMnpim_dE;//+nuetron ID by CDS
  //TH2F* pimmom_IMnpim_dE_n;//missing nuetron ID
  //TH2F* pimmom_MMnmiss_dE;//missing nuetron ID
  //TH2F* pimmom_MMnmiss_dE_woK0_woSidn;//missing nuetron ID
  //TH2F* pimmom_MMnmiss_dE_woK0_woSid_won;//missing nuetron ID
  //TH2F* pimmom_MMnmiss_dE_wK0_woSid_won;//missing nuetron ID
  //TH2F* pimmom_MMnmiss_dE_wK0_woSidn_won;//missing nuetron ID
  TH2F* pipmom_Momnpip_wSid_n;
  TH2F* pimmom_Momnpim_wSid_n;
  TH2F* IMpippim_DCApipibeam;//DCApipibeam distance of the center point of pi+ and pi- to beam in X-Y plane
  TH2F* IMpippim_DCApipibeam_n;
  TH2F* IMpippim_DCApipibeam_wK0_n;
  TH2F* IMpippim_DCApipibeam_wK0_woSid_n;
  TH2F* IMnpip_DCApipibeam;//DCApipibeam distance of the center point of pi+ and pi- to beam in X-Y plane
  TH2F* IMnpip_DCApipibeam_n;
  TH2F* IMnpip_DCApipibeam_woK0_n;
  TH2F* IMnpip_DCApipibeam_woK0_n_Sp;
  TH2F* IMnpim_DCApipibeam;//DCApipibeam distance of the center point of pi+ and pi- to beam in X-Y plane
  TH2F* IMnpim_DCApipibeam_n;
  TH2F* IMnpim_DCApipibeam_woK0_n;
  TH2F* IMnpim_DCApipibeam_woK0_n_Sm;
  TH2F* MMnmiss_DCApipibeam_wK0;
  TH2F* MMnmiss_DCApipibeam_woK0_wSid;
  TH2F* IMnpip_DCApip_dE_woK0_n;//
  TH2F* IMnpip_DCApim_dE_woK0_n;//
  TH2F* IMnpim_DCApip_dE_woK0_n;//
  TH2F* IMnpim_DCApim_dE_woK0_n;//
  TH2F* IMnpim_IMnpip_dE_n_Sp;//Spmode, avoiding crossing-point
  TH2F* IMnpim_IMnpip_dE_woK0_n_Sp;//Spmode, avoiding crossing-point
  //TH2F* IMnpim_IMnpip_dE_woK0_n_Sp_bg;//Spmode+background region, avoiding crossing-point
  TH2F* IMnpim_IMnpip_dE_n_Sm;//Smmode, avoiding crossing-point
  TH2F* IMnpim_IMnpip_dE_woK0_n_Sm;//Smmode, avoiding crossing-point
  //TH2F* IMnpim_IMnpip_dE_woK0_n_Sm_bg;//Smmode, avoiding crossing-point
  //TH2F* IMnpim_IMnpip_dE_woK0_woSidn;
  //TH2F* IMnpim_IMnpip_dE_n_side[ngap];//Side band for Sp + Sm mode
  //TH2F* IMnpim_IMnpip_dE_woK0_n_side[ngap];//Side band for Sp + Sm mode
  //TH2F* IMnpim_IMnpip_dE_n_Sp_side[2][ngap];//low mass,high mass
  //TH2F* IMnpim_IMnpip_dE_n_Sm_side[2][ngap];//low mass,high mass
  //TH2F* IMnpim_IMnpip_dE_woK0_n_Sp_side[2][ngap];//low mass,high mass
  //TH2F* IMnpim_IMnpip_dE_woK0_n_Sm_side[2][ngap];//low mass,high mass
  //TH2F* IMnpim_IMnpip_dE_n_Sp_sidewide[nzone];//scan from low mass side
  //TH2F* IMnpim_IMnpip_dE_n_Sm_sidewide[nzone];//scan from low mass side
  //TH2F* IMnpim_IMnpip_dE_woK0_n_Sp_sidewide[nzone];//scan from low mass side
  //TH2F* IMnpim_IMnpip_dE_woK0_n_Sm_sidewide[nzone];//scan from low mass side
  TH2F* nmom_IMnpip_dE_n;
  TH2F* nmom_IMnpip_dE_n_Sp;
  TH2F* nmom_IMnpip_dE_woK0_n;
  TH2F* nmom_IMnpip_dE_woK0_woSid_won;
  TH2F* nmom_IMnpip_dE_wK0_woSid_won;
  TH2F* nmom_IMnpim_dE_n;
  TH2F* nmom_IMnpim_dE_n_Sm;
  TH2F* nmom_IMnpim_dE_woK0_n;
  TH2F* nmom_IMnpim_dE_woK0_woSid_won;
  TH2F* nmom_IMnpim_dE_wK0_woSid_won;
  TH2F* nmom_Momnpip_woK0_n_Sp;
  TH2F* nmom_Momnpim_woK0_n_Sm;
  TH2F* MMnmiss_IMpippim_dE;
  TH2F* MMnmiss_IMpippim_dE_pat2;
  TH2F* MMnmiss_IMpippim_dE_pat7;
  TH2F* MMnmiss_IMpippim_dE_viciSp;
  TH2F* MMnmiss_IMpippim_dE_viciSm;
  TH2F* MMnmiss_IMpippim_dE_viciK0;
  TH2F* MMnmiss_IMpippim_dE_wSid;
  TH2F* MMnmiss_IMpippim_dE_wSid_fake;//for GEANT4 sim.
  TH2F* MMnmiss_IMpippim_dE_wSid_n;
  TH2F* MMnmiss_IMpippim_dE_wSid_n_fake;//for GEANT4 sim.
  TH2F* MMnmiss_IMpippim_dE_wSid_n_pat2;//for GEANT4 sim.
  TH2F* MMnmiss_IMpippim_dE_wSid_n_pat7;//for GEANT4 sim.
  TH2F* MMnmiss_IMpippim_dE_woK0;
  TH2F* MMnmiss_IMpippim_dE_woK0_wSid;
  TH2F* MMnmiss_IMpippim_dE_woK0_wSid_n;
  TH2F* MMnmiss_IMpippim_dE_wK0_wSid;
  TH2F* MMnmiss_IMpippim_dE_wK0_wSid_n;
  TH2F* Momnpim_Momnpip_dE_wSid_n;
  TH2F* Momnpim_Momnpip_dE_woK0_wSid_n;
  TH2F* Momnpim_Momnpip_dE_wK0_wSid_n;
  TH2F* Momnpim_Mompippim_dE_wSid_n;
  TH2F* Momnpim_Mompippim_dE_woK0_wSid_n;
  TH2F* Momnpim_Mompippim_dE_wK0_wSid_n;
  TH2F* Momnpip_Mompippim_dE_wSid_n;
  TH2F* Momnpip_Mompippim_dE_woK0_wSid_n;
  TH2F* Momnpip_Mompippim_dE_wK0_wSid_n;

  TH2F* MMnmiss_IMpippim_dE_woK0_woSid;
  TH2F* MMnmiss_IMpippim_dE_woSid;
  //TH2F* MMnmiss_IMpippim_dE_woSid_won;
  TH2F* MMnmiss_IMpippim_dE_woK0_woSid_won;
  //TH2F* MMnmiss_IMpippim_dE_woK0_woSidn;// (miss n & Sigma+/-) is rejected
  TH2F* MMnmiss_IMpippim_dE_wK0_woSid_won;
  TH2F* MMnmiss_IMpippim_dE_wK0_woSid_n_wbin[nwbin];
  TH2F* MMnmiss_IMpippim_dE_wK0_woSidn_won;
  //TH2F* MMnmiss_Mompippim_dE_woK0_woSidn;// (miss n & Sigma+/-) is rejected
  TH2F* MMnmiss_Mompippim_dE_woK0_woSid_won;// (miss n & Sigma+/-) is rejected
  TH2F* MMnmiss_Mompippim_dE_wK0_woSid_won;// (miss n & Sigma+/-) is rejected
  TH2F* MMnmiss_Mompippim_dE_wK0_woSidn_won;// (miss n & Sigma+/-) is rejected
  TH2F* MMnmiss_IMnpip_dE;
  TH2F* MMnmiss_IMnpip_dE_fake;
  TH2F* MMnmiss_IMnpip_dE_woSid;
  //TH2F* MMnmiss_IMnpip_dE_woSid_won;
  TH2F* MMnmiss_IMnpim_dE;
  TH2F* MMnmiss_IMnpim_dE_fake;//for mc
  TH2F* MMnmiss_IMnpim_dE_woSid;
  //TH2F* MMnmiss_IMnpim_dE_woSid_won;
  TH2F* MMnmiss_IMnpip_dE_woK0;
  TH2F* MMnmiss_IMnpim_dE_woK0;
  TH2F* MMnmiss_IMnpip_dE_woK0_woSm;
  TH2F* MMnmiss_IMnpip_dE_woK0_wSid_n_woSm_wbin[nwbin];
  TH2F* MMnmiss_IMnpip_dE_woK0_woSm_vici;//select vicinity of the signal 
  TH2F* MMnmiss_IMnpip_dE_woK0_woSm_viciext;//select vicinity of the signal, +extented
  TH2F* MMnmiss_IMnpip_dE_wK0_woSm;
  TH2F* MMnmiss_IMnpip_dE_woK0_woSm_cross;
  TH2F* MMnmiss_IMnpim_dE_woK0_woSp;
  TH2F* MMnmiss_IMnpim_dE_woK0_wSid_n_woSp_wbin[nwbin];
  TH2F* MMnmiss_IMnpim_dE_woK0_woSp_vici;
  TH2F* MMnmiss_IMnpim_dE_woK0_woSp_viciext;
  TH2F* MMnmiss_IMnpim_dE_wK0_woSp;
  TH2F* MMnmiss_IMnpim_dE_woK0_woSp_cross;
  TH2F* MMnmiss_IMnpip_dE_woK0_woSm_n;
  TH2F* MMnmiss_IMnpip_dE_wK0_woSm_n;
  TH2F* MMnmiss_IMnpim_dE_woK0_woSp_n;
  TH2F* MMnmiss_IMnpim_dE_wK0_woSp_n;
  TH2F* MMnmiss_IMnpip_dE_woK0_woSid;// (miss n & Sigma+/-) is rejected
  //TH2F* MMnmiss_IMnpip_dE_woK0_woSidn;// (miss n & Sigma+/-) is rejected
  TH2F* MMnmiss_IMnpip_dE_woK0_woSid_won;// (miss n & Sigma+/-) is rejected
  TH2F* MMnmiss_Momnpip_dE_woK0_woSid_won;// (miss n & Sigma+/-) is rejected
  //TH2F* MMnmiss_IMnpip_dE_woK0_woSidn_cross;// (miss n & Sigma+/-) is rejected + (miss n select) + (Sigma select);
  //TH2F* MMnmiss_IMnpip_dE_woK0_woSidn_cross_Sp;// (miss n & Sigma+/-) is rejected + (miss n select) + (SigmaP select);
  TH2F* MMnmiss_IMnpip_dE_wK0_woSid_won;// (miss n & Sigma+/-) is rejected
  TH2F* MMnmiss_Momnpip_dE_wK0_woSid_won;// (miss n & Sigma+/-) is rejected
  TH2F* MMnmiss_IMnpip_dE_wK0_woSidn_won;// (miss n & Sigma+/-) is rejected
  TH2F* MMnmiss_IMnpim_dE_woK0_woSid;// (miss n & Sigma+/-) is rejected
  //TH2F* MMnmiss_IMnpim_dE_woK0_woSidn;// (miss n & Sigma+/-) is rejected
  TH2F* MMnmiss_IMnpim_dE_woK0_woSid_won;// (miss n & Sigma+/-) is rejected
  TH2F* MMnmiss_Momnpim_dE_woK0_woSid_won;// (miss n & Sigma+/-) is rejected
  //TH2F* MMnmiss_IMnpim_dE_woK0_woSidn_cross;// (miss n & Sigma+/-) is rejected
  //TH2F* MMnmiss_IMnpim_dE_woK0_woSidn_cross_Sm;// (miss n & Sigma+/-) is rejected
  TH2F* MMnmiss_IMnpim_dE_wK0_woSid_won;// (miss n & Sigma+/-) is rejected
  TH2F* MMnmiss_Momnpim_dE_wK0_woSid_won;// (miss n & Sigma+/-) is rejected
  TH2F* MMnmiss_IMnpim_dE_wK0_woSidn_won;// (miss n & Sigma+/-) is rejected
  //TH2F* Momnpim_Momnpip_dE_woSid_won;
  TH2F* Momnpim_Momnpip_dE_woK0_woSid_won;
  TH2F* Momnpim_Momnpip_dE_wK0_woSid_won;
  //TH2F* Momnpim_Mompippim_dE_woSid_won;
  TH2F* Momnpim_Mompippim_dE_woK0_woSid_won;
  TH2F* Momnpim_Mompippim_dE_wK0_woSid_won;
  TH2F* Momnpip_Mompippim_dE_woK0_woSid_won;
  TH2F* Momnpip_Mompippim_dE_wK0_woSid_won;
  TH2F* MMnpim_MMnpip;
  TH2F* MMnpim_MMnpip_n;
  TH2F* MMnpim_MMnpip_woSid_n;
  TH2F* MMnpim_MMnpip_woK0_woSid_n;
  TH2F* MMnpim_MMnpip_mc;
  TH2F* MMnpim_MMnpip_wSid_n;
  TH2F* MMnpim_MMnpip_wSid_n_mc;
  TH2F* MMnpim_MMnpip_woK0_n;
  TH2F* MMnpim_MMnpip_woK0_wSid_n;
  TH2F* dE_CDHphi;
  TH2F* dE_CDHz;
  TH2F* dE_IMnpim;
  TH2F* dE_IMnpim_n;
  TH2F* dE_IMnpip;
  TH2F* dE_IMnpip_n;
  //TH2F* dE_IMnpim_woK0;
  //TH2F* dE_IMnpim_woK0_n;
  //TH2F* dE_IMnpip_woK0;
  //TH2F* dE_IMnpip_woK0_n;
  //TH2F* dE_IMnpipi_wSid_n;
  //TH2F* dE_IMnpipi_woK0_wSid_n;
  TH2F* Cosn_IMnpipi_wSid_n;
  TH2F* Cosn_IMnpipi_woK0_wSid_n;
  TH2F* MMnmiss_IMnpipi_wSid;//MM= missing Mass
  TH2F* MMnmiss_IMnpipi_wK0_wSid;
  TH2F* MMnmiss_IMnpipi_woK0_wSid;
  TH2F* MMnmiss_IMnpipi_woK0_wSid_Sp;
  TH2F* MMnmiss_IMnpipi_woK0_wSid_Sm;
  TH2F* MMnmiss_IMnpipi_wK0_woSid_won;
  TH2F* MMnmiss_Momnpipi_wK0_woSid_won;
  TH2F* MMnmiss_IMnpipi_wK0_woSidn_won;
  //TH2F* MMnmiss_IMnpipi_woK0_woSidn;
  TH2F* MMnmiss_IMnpipi_woK0_woSid_won;
  TH2F* MMnmiss_Momnpipi_woK0_woSid_won;
  TH2F* q_IMnpipi_wSid_n;
  const int nthetacut=36; //180/5=36
  TH2F* q_IMnpipi_wSid_n_thetacut[nthetacut];
  TH2F* q_IMnpipi_wSid_n_fake;
  TH2F* q_IMnpipi_wSid_n_fake_pat2;
  TH2F* q_IMnpipi_wSid_n_fake_pat7;
  TH2F* q_IMnpipi_wSid_n_wocross;
  TH2F* q_IMnpipi_woK0_wSid_n;
  TH2F* q_IMnpipi_woK0_wSid_n_woSp;
  TH2F* q_IMnpipi_woK0_wSid_n_woSm;
  TH2F* q_IMnpipi_wK0_wSid_n;
  TH2F* q_IMnpipi_wK0orwSid_n;
  TH2F* q_IMnpipi_wK0_woSid_n;
  TH2F* q_IMpiSigma_gen;//fine bins
  TH2F* q_IMpiSigma_wSid_n_genacc;//fine bins,  reaction data
  TH2F* q_IMnpipi_wSid_n_acc;//fine bins, McData
  TH2F* q_IMnpipi_wSid_n_acc_reco;//fine bins, reconstructed value
  TH2F* q_IMpiSigma_woK0_wSid_n_genacc;//fine bins, reaction data
  TH2F* q_IMnpipi_woK0_wSid_n_acc;//fine bins, no cuts for separationg S+/S-
  TH2F* q_IMnpipi_woK0_wSid_n_acc_reco;//fine bins, no cuts for separationg S+/S-, reconstructed value
  TH2F* q_IMnpipi_wSid_n_Sp;
  TH2F* q_IMnpipi_wSid_n_Sp_Stop;
  TH2F* q_IMnpipi_wSid_n_Sp_NoStop;
  TH2F* q_IMnpipi_wSid_n_SpSm;//Sigma+ like & Sigma- like mode
  TH2F* q_IMnpipi_wSid_n_Sp_woSm;
  TH2F* q_IMnpipi_wK0_n;
  TH2F* q_IMnpipi_woK0_woSid;
  TH2F* q_IMnpipi_woK0_woSid_n;
  TH2F* q_IMnpipi_woK0_woSid_won;
  //TH2F* q_nmom_woSid_won;
  TH2F* q_nmom_woK0_woSid_won;
  TH2F* q_nmom_wK0_woSid_won;
  TH2F* q_nmom_wK0_woSid_n_wbin[nwbin];
  TH2F* q_nmom_wSid_n;
  TH2F* q_nmom_woK0_wSid_n;
  TH2F* q_nmom_woK0_wSid_n_woSm_wbin[nwbin];
  TH2F* q_nmom_woK0_wSid_n_woSp_wbin[nwbin];
  TH2F* q_nmom_wK0_wSid_n;
  TH2F* q_IMnpipi_wK0_wSid_n_SpSm;
  TH2F* q_IMnpipi_wK0_wSid_n_Sp;
  TH2F* q_IMnpipi_woK0_wSid_n_Sp;
  TH2F* q_IMpiSigma_wSid_n_Sp_genacc;//fine bins
  TH2F* q_IMpiSigma_woK0_wSid_n_Sp_genacc;//fine bins
  TH2F* q_IMnpipi_wSid_n_Sp_acc;//fine bins
  TH2F* q_IMnpipi_woK0_wSid_n_Sp_acc;//fine bins
  TH2F* q_IMnpipi_wSid_n_Sp_acc_reco;//fine bins
  TH2F* q_IMnpipi_woK0_wSid_n_Sp_acc_reco;//fine bins
  //TH2F* q_IMnpipi_wSid_n_Sp_side[3][2][ngap];//sideband type,
  //TH2F* q_IMnpipi_woK0_wSid_n_Sp_side[3][2][ngap];//sideband type,
  //low high side,
  //distance to signal
  //TH2F* q_IMnpipi_wSid_n_Sp_sidewide[nzone];//
  //TH2F* q_IMnpipi_woK0_wSid_n_Sp_sidewide[nzone];//
  TH2F* q_IMnpipi_wSid_n_Sm;
  TH2F* q_IMnpipi_wSid_n_Sm_Stop;
  TH2F* q_IMnpipi_wSid_n_Sm_NoStop;
  TH2F* q_IMnpipi_wSid_n_Sm_woSp;
  TH2F* q_IMnpipi_wK0_wSid_n_Sm;
  TH2F* q_IMnpipi_woK0_wSid_n_Sm;
  TH2F* q_IMpiSigma_wSid_n_Sm_genacc;//fine bins, reaction data
  TH2F* q_IMpiSigma_woK0_wSid_n_Sm_genacc;//fine bins, reaction data
  TH2F* q_IMnpipi_wSid_n_Sm_acc;//fine bins, MCdata
  TH2F* q_IMnpipi_woK0_wSid_n_Sm_acc;//fine bins, MCdata
  TH2F* q_IMnpipi_wSid_n_Sm_acc_reco;//fine bins, reconstructed value
  TH2F* q_IMnpipi_woK0_wSid_n_Sm_acc_reco;//fine bins, reconstructed value
  //TH2F* q_IMnpipi_wSid_n_Sm_side[3][2][ngap];//sideband type,
  //TH2F* q_IMnpipi_woK0_wSid_n_Sm_side[3][2][ngap];//sideband type,
  //low high side
  //distance to signal
  //TH2F* q_IMnpipi_wSid_n_Sm_sidewide[nzone];//
  //TH2F* q_IMnpipi_woK0_wSid_n_Sm_sidewide[nzone];//
  TH2F* q_IMnpipi_wK0_woSid_won;
  TH2F* q_IMnpipi_wK0_woSidn_won;
  //TH2F* q_IMnpipi_woK0_woSidn;
  //TH2F* q_IMnpipi_woK0_woSidn_cross;
  //TH2F* q_IMnpipi_woK0_woSidn_cross_Sp;
  //TH2F* q_IMnpipi_woK0_woSidn_cross_Sm;
  TH2F* IMnpip_IMnpipi_n;
  TH2F* IMnpim_IMnpipi_n;
  TH2F* IMnpip_IMnpipi_wK0_n;
  TH2F* IMnpim_IMnpipi_wK0_n;
  TH2F* IMnpip_IMnpipi_woK0_n;
  TH2F* IMnpim_IMnpipi_woK0_n;

  //TH2F* q_IMnpipi_wSid_n_side[ngap];//side band method
  //TH2F* q_IMnpipi_woK0_wSid_n_side[ngap];//side band method
  TH2F* pipmom_IMnpipi_wSid_n;
  TH2F* pimmom_IMnpipi_wSid_n;
  TH2F* pipmom_IMnpipi_woK0_wSid_n;
  TH2F* pimmom_IMnpipi_woK0_wSid_n;
  TH2F* nmom_CDHphi;
  TH2F* nmom_cosn_wK0;//cosn = CDH neutron angle to beam axis.
  //TH2F* nmom_cosn_wK0_n;//
  TH2F* nmom_cosn_wSid_n;//cosn = CDH neutron angle to beam axis.
  TH2F* nmom_cosn_woK0_wSid_n;//cosn = CDH neutron angle to beam axis.
  //TH2F* nmom_cosn_woK0_woSidn;//cosn = CDH neutron angle to beam axis.
  TH2F* nmom_cosn_woK0_woSid_won;//cosn = CDH neutron angle to beam axis.
  TH2F* nmom_cosn_wK0_woSid_won;//cosn = CDH neutron angle to beam axis.
  TH2F* nmom_cosn_wK0_n_forward;//w/ forward angle selection
  TH2F* nmom_cosK0_wK0;
  TH2F* nmom_cosK0_wK0_n;
  TH2F* nmom_cosK0n_wK0;// cosK0n = CDH neutron angle to K0
  TH2F* nmom_cosK0n_wK0_n;//
  TH2F* nmom_cosnmiss_wK0_n;//
  TH2F* nmom_cosnmiss_wSid_n;//
  TH2F* nmom_cosnmiss_woK0_wSid_n;//
  TH2F* nmom_cosnnmiss_wK0_n;//

  //TH2F* nmom_pipmom_woK0_woSidn;
  TH2F* nmom_pipmom_woK0_woSid_won;
  TH2F* nmom_pipmom_wK0_woSid_won;
  //TH2F* nmom_pimmom_woK0_woSidn;
  TH2F* nmom_pimmom_woK0_woSid_won;
  TH2F* nmom_pimmom_wK0_woSid_won;
  //TH2F* nmom_cospim_woK0_woSidn;//cosn = CDH neutron angle to beam axis.
  TH2F* nmom_cospim_woK0_woSid_won;//cosn = CDH neutron angle to beam axis.
  TH2F* nmom_coslabpim_wK0_woSid_won;//cosn = CDH neutron angle to beam axis.
  TH2F* nmom_cospim_wK0_woSid_won;//cosn = CDH neutron angle to beam axis.
  //TH2F* nmom_cospip_woK0_woSidn;//cosn = CDH neutron angle to beam axis.
  TH2F* nmom_cospip_woK0_woSid_won;//cosn = CDH neutron angle to beam axis.
  TH2F* nmom_cospippim_woK0_woSid_won;//cosn = CDH neutron angle to beam axis.
  TH2F* nmom_cospip_wK0_woSid_won;//cosn = CDH neutron angle to beam axis.
  TH2F* nmom_cospippim_wK0_woSid_won;//cosn = CDH neutron angle to beam axis.
  //TH2F* nmom_phinpim_woK0_woSidn;//cosn = CDH neutron angle to beam axis.
  TH2F* nmom_phinpim_woK0_woSid_won;//cosn = CDH neutron angle to beam axis.
  TH2F* nmom_phinpim_wK0_woSid_won;//cosn = CDH neutron angle to beam axis.
  //TH2F* nmom_phinpip_woK0_woSidn;//cosn = CDH neutron angle to beam axis.
  TH2F* nmom_phinpip_woK0_woSid_won;//cosn = CDH neutron angle to beam axis.
  TH2F* nmom_phinpip_wK0_woSid_won;//cosn = CDH neutron angle to beam axis.
  //TH2F* nmom_phipim_woK0_woSidn;//cosn = CDH neutron angle to beam axis.
  TH2F* nmom_phipim_woK0_woSid_won;//cosn = CDH neutron angle to beam axis.
  TH2F* nmom_phipim_wK0_woSid_won;//cosn = CDH neutron angle to beam axis.
  //TH2F* nmom_phipip_woK0_woSidn;//cosn = CDH neutron angle to beam axis.
  TH2F* nmom_phipip_woK0_woSid_won;//cosn = CDH neutron angle to beam axis.
  TH2F* nmom_phipip_wK0_woSid_won;//cosn = CDH neutron angle to beam axis.
  //TH2F* nmom_phin_woK0_woSidn;//cosn = CDH neutron angle to beam axis.
  TH2F* nmom_phin_woK0_woSid_won;//cosn = CDH neutron angle to beam axis.
  TH2F* nmom_phin_wK0_woSid_won;//cosn = CDH neutron angle to beam axis.
  TH2F* K0mom_cosK0_wK0;
  TH2F* K0mom_cosK0_wK0_n;
  TH2F* nmissmom_cosnmiss_wK0;// cosnmiss = missing neutron angle to beam axis.
  TH2F* nmissmom_cosnmiss_wK0_n;//
  //TH2F* nmissmom_cosK0nmiss_wK0;// cosnmiss = missing neutron angle to K0.
  TH2F* nmissmom_cosK0nmiss_wK0_n;//
  TH2F* nmom_K0mom;//K0 selection,
  TH2F* nmom_K0mom_n;//K0 selection, missing neutron cut
  TH2F* nmom_K0mom_woSid_n;//K0 selection, missing neutron cut
  TH2F* nmom_nmissmom_wK0;//K0 selection
  TH2F* nmom_nmissmom_wK0_n;//K0 selection, missing neutron cut
  TH2F* nmom_nmissmom_wK0_wSid;//Sigma selection, missing neutron cut
  TH2F* nmom_nmissmom_woK0_wSid;//Sigma selection, missing neutron cut
  TH2F* nmom_nmissmom_wK0_wSid_n;//Sigma selection, missing neutron cut
  TH2F* nmom_nmissmom_woK0_wSid_n;//Sigma selection, missing neutron cut
  TH2F* nmissmom_K0mom;
  TH2F* nmissmom_K0mom_n;


  TH2F* nmom_IMnpipi_wSid_n;
  //TH2F* nmom_IMnpipi_woK0_wSid_n;
  TH2F* nmom_IMnpipi_woK0_woSid_won;
  TH2F* nmom_IMnpipi_woK0_wSid_n_Sp;
  TH2F* nmom_IMnpipi_woK0_wSid_n_Sm;
  TH2F* nmom_MMnmiss_wSid;
  TH2F* nmom_MMnmiss_wSid_fake;
  TH2F* nmom_MMnmiss_wSid_n;
  TH2F* nmom_MMnmiss_woK0_wSid;
  TH2F* nmom_MMnmiss_woK0_wSid_Sp;//
  TH2F* nmom_MMnmiss_woK0_wSid_Sm;//
  TH2F* nmom_MMnmiss_woK0_woSid;
  //TH2F* nmom_MMnmiss_woK0_woSidn;
  TH2F* nmom_MMnmiss_woK0_woSid_won;
  //K0 study
  TH2F* nmom_IMnpipi_wK0_n;
  TH2F* nmom_IMnpipi_wK0_wSid_n;
  TH2F* nmom_IMnpipi_wK0_woSid_won;
  TH2F* nmom_IMnpipi_wK0_woSid_n;
  TH2F* nmom_MMnmiss_wK0;
  TH2F* nmom_MMnmiss_wK0_woSid;
  TH2F* nmom_MMnmiss_wK0_woSid_won;
  TH2F* nmom_MMnmiss_wK0_woSidn_won;
  TH2F* nmom_cosnlab_K0_n;//ncds mom. vs cos ncdds lab. frame K0
  TH2F* nmom_IMpippim;
  TH2F* nmom_MK0bar2;
  TH2F* nmom_IMpippim_n;
  TH2F* nmom_IMpippim_wSid_n;
  TH2F* nmom_IMpippim_woK0_woSid_won;
  TH2F* nmom_IMpippim_wK0_woSid_won;
  TH2F* mnmom_IMpippim_n;//missing neutron mom.
  TH2F* mnmom_IMpippim_wSid_n;//missing neutron mom.
  TH2F* q_IMpippim_n;
  TH2F* q_IMpippim_wSid_n;
  TH2F* q_IMpippim_woK0_woSid_won;
  TH2F* q_IMpippim_wK0_woSid_won;
  TH2F* q_IMnpip_n;
  TH2F* q_IMnpip_n_wSid;
  TH2F* q_IMnpip_woK0_woSid_won;
  TH2F* q_IMnpip_wK0_woSid_won;
  TH2F* q_IMnpim_n;
  TH2F* q_IMnpim_n_wSid;
  TH2F* q_IMnpim_woK0_woSid_won;
  TH2F* q_IMnpim_wK0_woSid_won;
  //TH2F* q_MMnmiss;
  //TH2F* q_MMnmiss_wSid;
  //TH2F* q_MMnmiss_n_wSid;
  //TH2F* q_MMnmiss_woK0_woSid_won;
  //TH2F* q_MMnmiss_wK0_woSid_won;
  //TH2F* q_nmom_n;
  //TH2F* q_nmom_n_wSid;
  TH2F* IMpippim_IMnpipi_n;
  TH2F* IMpippim_IMnpipi_n_wSid;
  TH2F* IMpippim_IMnpipi_n_wSid_woSp;
  TH2F* IMpippim_IMnpipi_n_wSid_woSm;
  TH2F* IMpippim_IMnpip_vici;
  TH2F* IMpippim_IMnpip_vici_woSm;
  TH2F* IMpippim_IMnpip_n;
  //TH2F* IMpippim_IMnpip_n_bin[nbintemplate];
  TH2F* IMpippim_IMnpip_n_woSm;
  //TH2F* IMpippim_IMnpip_n_woSm_bin[nbintemplate];
  TH2F* IMpippim_IMnpip_n_woSmdia;
  //TH2F* IMpippim_IMnpip_n_woSmdia_bin[nbintemplate];
  TH2F* IMpippim_IMnpip_wSid_n_Sp;
  TH2F* IMpippim_IMnpip_wSid_n_Sp_bin[nbintemplate];
  TH2F* IMpippim_IMnpip_wSid_n_Sp_wbin[nwbin];
  TH2F* IMpippim_IMnpip_wSid_n_woSm;
  //TH2F* IMpippim_IMnpip_wSid_n_woSm_bin[nbintemplate];
  TH2F* IMpippim_IMnpip_wK0_n;
  TH2F* IMpippim_IMnpip_wK0_n_bin[nbintemplate];
  TH2F* IMpippim_IMnpip_wK0_n_wbin[nwbin];
  TH2F* IMpippim_IMnpip_wK0_n_woSmdia;
  //TH2F* IMpippim_IMnpip_wK0_n_woSmdia_bin[nbintemplate];
  TH2F* IMpippim_IMnpip_wK0orwSid_n;
  TH2F* IMpippim_IMnpip_wK0orwSid_n_bin[nbintemplate];
  TH2F* IMpippim_IMnpip_wK0orwSid_n_wbin[nwbin];
  TH2F* IMpippim_IMnpip_wK0orwSid_n_woSmdia;
  //TH2F* IMpippim_IMnpip_wK0orwSid_n_woSmdia_bin[nbintemplate];
  TH2F* IMpippim_IMnpip_woK0_wSid_n;
  TH2F* IMpippim_IMnpip_woK0_wSid_n_woSp_wbin[nwbin];//for Sm template
  TH2F* IMpippim_IMnpip_woK0_wSid_n_woSm_wbin[nwbin];//for Sp template
  TH2F* IMpippim_IMnpim_vici;
  TH2F* IMpippim_IMnpim_vici_woSp;
  TH2F* IMpippim_IMnpim_n;
  //TH2F* IMpippim_IMnpim_n_bin[nbintemplate];
  TH2F* IMpippim_IMnpim_n_woSp;
  //TH2F* IMpippim_IMnpim_n_woSp_bin[nbintemplate];
  TH2F* IMpippim_IMnpim_n_woSpdia;
  //TH2F* IMpippim_IMnpim_n_woSpdia_bin[nbintemplate];
  TH2F* IMpippim_IMnpim_wSid_n_Sm;
  TH2F* IMpippim_IMnpim_wSid_n_Sm_bin[nbintemplate];
  TH2F* IMpippim_IMnpim_wSid_n_Sm_wbin[nwbin];
  TH2F* IMpippim_IMnpim_wSid_n_woSp;
  //TH2F* IMpippim_IMnpim_wSid_n_woSp_bin[nbintemplate];
  TH2F* IMpippim_IMnpim_wK0_n;
  TH2F* IMpippim_IMnpim_wK0_n_bin[nbintemplate];
  TH2F* IMpippim_IMnpim_wK0_n_wbin[nwbin];
  TH2F* IMpippim_IMnpim_wK0_n_woSpdia;
  //TH2F* IMpippim_IMnpim_wK0_n_woSpdia_bin[nbintemplate];
  TH2F* IMpippim_IMnpim_wK0orwSid_n;
  TH2F* IMpippim_IMnpim_wK0orwSid_n_bin[nbintemplate];
  TH2F* IMpippim_IMnpim_wK0orwSid_n_wbin[nwbin];
  TH2F* IMpippim_IMnpim_wK0orwSid_n_woSpdia;
  //TH2F* IMpippim_IMnpim_wK0orwSid_n_woSpdia_bin[nbintemplate];
  TH2F* IMpippim_IMnpim_woK0_wSid_n;
  TH2F* IMpippim_IMnpim_woK0_wSid_n_woSp_wbin[nwbin];
  TH2F* IMpippim_IMnpim_woK0_wSid_n_woSm_wbin[nwbin];
  TH2F* IMpippim_IMnpip_n_cross;
  TH2F* IMpippim_IMnpim_n_cross;
  TH2F* IMpippim_IMnpip_woK0_woSid_won;
  TH2F* IMpippim_IMnpip_wK0_woSid_n_wbin[nwbin];
  TH2F* IMpippim_IMnpim_wK0_woSid_n_wbin[nwbin];
  TH2F* IMpippim_IMnpip_wK0_woSid_won;
  TH2F* IMpippim_IMnpim_woK0_woSid_won;
  TH2F* IMpippim_IMnpim_wK0_woSid_won;

  TH2F* Mompippim_IMnpipi_dE_wK0_n;
  TH2F* Mompippim_nmom_dE_wK0_n;
  TH2F* Mompippim_IMnpipi_dE_wK0_n_Sp;
  TH2F* Mompippim_IMnpipi_dE_wK0_n_Sm;
  TH2F* diff2d_Phipippim_Phinpim;
  //TH2F* diff2d_Zpippim_Znpim;
  //TH2F* diff2d_Phipippim_Phinpim_woSid_won;
  TH2F* diff2d_Zpippim_Znpim_woSid_won;
  TH2F* diff2d_CDC_CDH_pim;//neutral hit - pim track projected position
  TH2F* diff2d_CDC_CDH_pip;//neutral hit - pip track projected position
  TH2F* diff2d_CDC_CDH_pim_phi_tof;//neutral hit - pim track projected position (phi) vs tof(ncan - pim)
  TH2F* diff2d_CDC_CDH_pip_phi_tof;//neutral hit - pip track projected position (phi) vs tof(ncan - pip)
  TH2F* diff2d_CDC_CDH_pim_z_tof;//neutral hit - pim track projected positio (z) vs tof (ncan - -pim)
  TH2F* diff2d_CDC_CDH_pip_z_tof;//neutral hit - pip track projected position (z) vs tof (ncan -pip)
  //TH2F* diff2d_CDC_CDH_pim_woSid_won;//neutral hit - pim track projected position
  //TH2F* diff2d_CDC_CDH_pip_woSid_won;//neutral hit - pip track projected position
  TH2F* diff2d_CDC_CDH_pim_phi_tof_woSid_won;//neutral hit - pim track projected position
  TH2F* diff2d_CDC_CDH_pip_phi_tof_woSid_won;//neutral hit - pip track projected position
  TH2F* diff2d_CDC_CDH_pim_z_tof_woSid_won;//neutral hit - pim track projected position
  TH2F* diff2d_CDC_CDH_pip_z_tof_woSid_won;//neutral hit - pip track projected position
  TH2F* diff2d_CDC_CDH_pim_wK0;//neutral hit - pim track projected position
  TH2F* diff2d_CDC_CDH_pip_wK0;//neutral hit - pip track projected position
  TH2F* diff2d_CDC_CDH_pim_woK0_wSid;//neutral hit - pim track projected position
  TH2F* diff2d_CDC_CDH_pip_woK0_wSid;//neutral hit - pip track projected position
  TH2F* diff2d_CDC_CDH_pim_woK0_wSid_n;//neutral hit - pim track projected position
  TH2F* diff2d_CDC_CDH_pip_woK0_wSid_n;//neutral hit - pip track projected position
  TH2F* dE_diffphi_CDC_CDH_pim;
  TH2F* dE_diffphi_CDC_CDH_pip;
  //TH2F* MMnmiss_dE;
  TH2F* MMnmiss_diffphi_CDC_CDH_pim;
  TH2F* MMnmiss_diffphi_CDC_CDH_pip;
  TH2F* MMnmiss_diffphi_CDC_CDH_pim_woK0_wSid;
  TH2F* MMnmiss_diffphi_CDC_CDH_pip_woK0_wSid;
  TH2F* MMnmiss_diffz_CDC_CDH_pim_woK0_wSid;
  TH2F* MMnmiss_diffz_CDC_CDH_pip_woK0_wSid;
  TH2F* pimmom_diffphi_CDC_CDH_pim;
  TH2F* pimmom_diffphi_CDC_CDH_pim_woSid_won;
  TH2F* pipmom_diffphi_CDC_CDH_pip;
  TH2F* pipmom_diffphi_CDC_CDH_pip_woSid_won;
  TH2F* pimmom_diffz_CDC_CDH_pim;
  TH2F* pipmom_diffz_CDC_CDH_pip;
  TH2F* nmom_diffphi_CDC_CDH_pim;
  TH2F* nmom_diffphi_CDC_CDH_pip;
  TH2F* nmom_diffz_CDC_CDH_pim;
  TH2F* nmom_diffz_CDC_CDH_pip;
  TH2F* pimmom_diffphi_CDC_CDH_pim_woK0_wSid;
  TH2F* pipmom_diffphi_CDC_CDH_pip_woK0_wSid;
  TH2F* pimmom_diffz_CDC_CDH_pim_woK0_wSid;
  TH2F* pipmom_diffz_CDC_CDH_pip_woK0_wSid;
  TH2F* nmom_diffphi_CDC_CDH_pim_woK0_wSid;
  TH2F* nmom_diffphi_CDC_CDH_pip_woK0_wSid;
  TH2F* nmom_diffz_CDC_CDH_pim_woK0_wSid;
  TH2F* nmom_diffz_CDC_CDH_pip_woK0_wSid;

  std::cout << __LINE__ << std::endl;

  const int nbinIMnpipi = 80;//1.2-2 GeV/c^2
  const float IMnpipilow = 1.2;
  const float IMnpipihi = 2.0;
  const int nbinq = 100;// 25;//0-1.5 GeV/c
  const int nbinIMnpi = 500; //1-2 GeV/c^2
  const int nbinnmiss = 150; //0-1.5 GeV/c
  const double nmisslow = 0.0;
  const double nmisshigh = 1.5;
  const int nbindE = 200;
  const int nbinpippim = 900;//0-0.9 GeV/c^2
  const int nbinnmom = 100;
  const int nbinmomnpi = 150;
  const int nbinmompipi = 150;
  
  const float wbinlow[nwbin] = {1.0,1.40,1.52};
  const float wbinhigh[nwbin] = {1.40,1.52,2.00};
  //bin width to be adjusted to resolution
  const double IMnpip_wbinlow = anacuts::Sigmap_center-17.0*2.0*anacuts::Sigmap_sigma;
  const double IMnpip_wbinhigh = anacuts::Sigmap_center+47.0*2.0*anacuts::Sigmap_sigma;
  const int nbinIMnpip_wbin = 32;
  const double IMnpim_wbinlow = anacuts::Sigmam_center-17.0*2.0*anacuts::Sigmam_sigma;
  const double IMnpim_wbinhigh = anacuts::Sigmam_center+41.0*2.0*anacuts::Sigmam_sigma;
  const int nbinIMnpim_wbin = 29;
  const double IMpippim_wbinlow = anacuts::K0_center-19.0*2.0*anacuts::K0_sigma;
  const double IMpippim_wbinhigh = anacuts::K0_center+27.0*2.0*anacuts::K0_sigma;
  const int nbinIMpippim_wbin = 23;
  
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

  dE_betainv_fid = new TH2F(Form("dE_betainv_fid"),Form("dE_betainv_fid"),100, 0, 50, nbindE, 0, 50);
  dE_betainv_fid->SetXTitle("1/#beta");
  dE_betainv_fid->SetYTitle("dE [MeVee]");

  dE_nmom_fid_beta = new TH2F("dE_nmom_fid_beta","dE_nmom_fid_beta",100, 0, 1.5, 500, 0, 50);
  dE_nmom_fid_beta->SetXTitle("Mom. [GeV/c]");
  dE_nmom_fid_beta->SetYTitle("dE [MeVee]");

  //dE_nmom_fid_beta_wK0 = new TH2F("dE_nmom_fid_beta_wK0","dE_nmom_fid_beta_wK0",100, 0, 1.5, 500, 0, 50);
  //dE_nmom_fid_beta_wK0->SetXTitle("Mom. [GeV/c]");
  //dE_nmom_fid_beta_wK0->SetYTitle("dE [MeVee]");

  dE_MMom_fid_beta = new TH2F(Form("dE_MMom_fid_beta"),Form("dE_MMom_fid_beta"),100, 0, 1.5, nbindE, 0, 50);
  dE_MMom_fid_beta->SetXTitle("Missing Mom. [GeV/c]");
  dE_MMom_fid_beta->SetYTitle("dE [MeVee]");

  //dE_MMom_fid_beta_woK0 = new TH2F(Form("dE_MMom_fid_beta_woK0"),Form("dE_MMom_fid_beta_woK0"),100, 0, 1.5, nbindE, 0, 50);
  //dE_MMom_fid_beta_woK0->SetXTitle("Missing Mom. [GeV/c]");
  //dE_MMom_fid_beta_woK0->SetYTitle("dE [MeVee]");

  dE_MMass_fid_beta = new TH2F(Form("dE_MMass_fid_beta"),Form("dE_MMass_fid_beta"), nbinnmiss, nmisslow, nmisshigh, nbindE, 0, 50);
  dE_MMass_fid_beta->SetXTitle("Missing mass [GeV/c^{2}]");
  dE_MMass_fid_beta->SetYTitle("dE [MeVee]");

  //dE_MMass_fid_beta_woK0 = new TH2F(Form("dE_MMass_fid_beta_woK0"),Form("dE_MMass_fid_beta_woK0"), nbinnmiss, nmisslow, nmisshigh, nbindE, 0, 50);
  //dE_MMass_fid_beta_woK0->SetXTitle("Missing mass [GeV/c^{2}]");
  //dE_MMass_fid_beta_woK0->SetYTitle("dE [MeVee]");

  dE_MMass_fid_beta_wSid = new TH2F(Form("dE_MMass_fid_beta_wSid"),Form("dE_MMass_fid_beta_wSid"), nbinnmiss, nmisslow, nmisshigh, nbindE, 0, 10);
  dE_MMass_fid_beta_wSid->SetXTitle("Missing mass [GeV/c^{2}]");
  dE_MMass_fid_beta_wSid->SetYTitle("dE [MeVee]");

  //dE_MMass_fid_beta_woK0_wSid = new TH2F(Form("dE_MMass_fid_beta_woK0_wSid"),Form("dE_MMass_fid_beta_woK0_wSid"), nbinnmiss, nmisslow, nmisshigh, nbindE, 0, 10);
  //dE_MMass_fid_beta_woK0_wSid->SetXTitle("Missing mass [GeV/c^{2}]");
  //dE_MMass_fid_beta_woK0_wSid->SetYTitle("dE [MeVee]");

  MMom_MMass = new TH2F(Form("MMom_MMass"),Form("MMom_MMass"), nbinnmiss, nmisslow, nmisshigh, 200, 0, 2.0);
  MMom_MMass->SetXTitle("Missing Mass [GeV/c^{2}]");
  MMom_MMass->SetYTitle("Missing Mom. [GeV/c]");

  MMom_MMass_woK0 = new TH2F(Form("MMom_MMass_woK0"),Form("MMom_MMass_woK0"), nbinnmiss, nmisslow, nmisshigh, 200, 0, 2.0);
  MMom_MMass_woK0->SetXTitle("Missing Mass [GeV/c^{2}]");
  MMom_MMass_woK0->SetYTitle("Missing Mom. [GeV/c]");

  MMom_MMass_wSid = new TH2F(Form("MMom_MMass_wSid"),Form("MMom_MMass_wSid"), nbinnmiss, nmisslow, nmisshigh, 200, 0, 2.0);
  MMom_MMass_wSid->SetXTitle("Missing Mass [GeV/c^{2}]");
  MMom_MMass_wSid->SetYTitle("Missing Mom. [GeV/c]");
  
  //MMom_MMass_wSid_n = new TH2F(Form("MMom_MMass_wSid_n"),Form("MMom_MMass_wSid_n"), nbinnmiss, nmisslow, nmisshigh, 200, 0, 2.0);
  //MMom_MMass_wSid_n->SetXTitle("Missing Mass [GeV/c^{2}]");
  //MMom_MMass_wSid_n->SetYTitle("Missing Mom. [GeV/c]");

  MMom_MMass_woK0_wSid = new TH2F("MMom_MMass_woK0_wSid","MMom_MMass_woK0_wSid", nbinnmiss, nmisslow, nmisshigh, 200, 0, 2.0);
  MMom_MMass_woK0_wSid->SetXTitle("Missing Mass [GeV/c^{2}]");
  MMom_MMass_woK0_wSid->SetYTitle("Missing Mom. [GeV/c]");

  //MMom_MMass_woK0_woSidn = new TH2F("MMom_MMass_woK0_woSidn","MMom_MMass_woK0_woSidn", nbinnmiss, nmisslow, nmisshigh, 200, 0, 2.0);
  //MMom_MMass_woK0_woSidn->SetXTitle("Missing Mass [GeV/c^{2}]");
  //MMom_MMass_woK0_woSidn->SetYTitle("Missing Mom. [GeV/c]");

  MMom_MMass_woK0_woSid_won = new TH2F("MMom_MMass_woK0_woSid_won","MMom_MMass_woK0_woSid_won", nbinnmiss, nmisslow, nmisshigh, 200, 0, 2.0);
  MMom_MMass_woK0_woSid_won->SetXTitle("Missing Mass [GeV/c^{2}]");
  MMom_MMass_woK0_woSid_won->SetYTitle("Missing Mom. [GeV/c]");

  MMom_MMass_wK0_woSid_won = new TH2F("MMom_MMass_wK0_woSid_won","MMom_MMass_wK0_woSid_won", nbinnmiss, nmisslow, nmisshigh, 200, 0, 2.0);
  MMom_MMass_wK0_woSid_won->SetXTitle("Missing Mass [GeV/c^{2}]");
  MMom_MMass_wK0_woSid_won->SetYTitle("Missing Mom. [GeV/c]");
  
  MMnmiss_react = new TH1F("MMnmiss_react","MMnmiss_react",100,0,1.5);
  MMnmiss_react->SetXTitle("react. Missing Mass [GeV/c^{2}");


  MMom_MMass_wK0_woSidn_won = new TH2F("MMom_MMass_wK0_woSidn_won","MMom_MMass_wK0_woSidn_won", nbinnmiss, nmisslow, nmisshigh, 200, 0, 2.0);
  MMom_MMass_wK0_woSidn_won->SetXTitle("Missing Mass [GeV/c^{2}]");
  MMom_MMass_wK0_woSidn_won->SetYTitle("Missing Mom. [GeV/c]");
  
  diff_IMnpim_reactmc = new TH1F("diff_IMnpim_reactmc","diff_IMnpim_reactmc",100,-1.0,1.0);
  diff_IMnpim_reactmc->SetXTitle("diff. IMnpim (MCData - Reac.) [GeV/c^{2}]");
  
  IMnpim_MMnpim_mc = new TH2F("IMnpim_MMnpim_mc","IMnpim_MMnpim_mc",140,1.,1.7,140,1.0,1.7);
  IMnpim_MMnpim_mc->SetYTitle("true  IM(n_{CDS}#pi^{-}) [GeV/c^{2}]");
  IMnpim_MMnpim_mc->SetXTitle("true  IM(n_{miss}#pi^{-}) [GeV/c^{2}]");
  
  IMnpim_MMnpim_mc_wSid_n_fake = new TH2F("IMnpim_MMnpim_mc_wSid_n_fake","IMnpim_MMnpim_mc_wSid_n_fake",140,1.,1.7,140,1.0,1.7);
  IMnpim_MMnpim_mc_wSid_n_fake->SetYTitle("true  IM(n_{CDS}#pi^{-}) [GeV/c^{2}]");
  IMnpim_MMnpim_mc_wSid_n_fake->SetXTitle("true  IM(n_{miss}#pi^{-}) [GeV/c^{2}]");
  
  IMnpim_MMnpim_mc_wSid_n_fake_pat2 = new TH2F("IMnpim_MMnpim_mc_wSid_n_fake_pat2","IMnpim_MMnpim_mc_wSid_n_fake_pat2",140,1.,1.7,140,1.0,1.7);
  IMnpim_MMnpim_mc_wSid_n_fake_pat2->SetYTitle("true  IM(n_{CDS}#pi^{-}) [GeV/c^{2}]");
  IMnpim_MMnpim_mc_wSid_n_fake_pat2->SetXTitle("true  IM(n_{miss}#pi^{-}) [GeV/c^{2}]");
  
  IMnpim_MMnpim_mc_wSid_n_fake_pat7 = new TH2F("IMnpim_MMnpim_mc_wSid_n_fake_pat7","IMnpim_MMnpim_mc_wSid_n_fake_pat7",140,1.,1.7,140,1.0,1.7);
  IMnpim_MMnpim_mc_wSid_n_fake_pat7->SetYTitle("true  IM(n_{CDS}#pi^{-}) [GeV/c^{2}]");
  IMnpim_MMnpim_mc_wSid_n_fake_pat7->SetXTitle("true  IM(n_{miss}#pi^{-}) [GeV/c^{2}]");
  
  
  IMnpim_MMnpim_mc_vtx = new TH2F("IMnpim_MMnpim_mc_vtx","IMnpim_MMnpim_mc_vtx",140,1.,1.7,140,1.0,1.7);
  IMnpim_MMnpim_mc_vtx->SetYTitle("true  IM(n_{CDS}#pi^{-})  [GeV/c^{2}]");
  IMnpim_MMnpim_mc_vtx->SetXTitle("true  IM(n_{miss}#pi^{-}) [GeV/c^{2}]");
  
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
  
  diff2D_MMnmiss_IMnpim_reactmc = new TH2F("diff2D_MMnmiss_IMnpim_reactmc","diff2D_MMnmiss_IMnpim_reactmc",200,-1.0,1.0,150,-1.5,1.5);
  diff2D_MMnmiss_IMnpim_reactmc->SetXTitle("diff. IMnpim (MCData - Reac.) [GeV/c^{2}]");
  diff2D_MMnmiss_IMnpim_reactmc->SetYTitle("diff. Missing Mass (MCData - Reac.) [GeV/c^{2}]");
   
  diff2D_MMnmiss_IMnpim_wSid_n_reactmc = new TH2F("diff2D_MMnmiss_IMnpim_wSid_n_reactmc","diff2D_MMnmiss_IMnpim_wSid_n_reactmc",200,-1.0,1.0,150,-1.5,1.5);
  diff2D_MMnmiss_IMnpim_wSid_n_reactmc->SetXTitle("diff. IMnpim (MCData - Reac.) [GeV/c^{2}]");
  diff2D_MMnmiss_IMnpim_wSid_n_reactmc->SetYTitle("diff. Missing Mass (MCData - Reac.) [GeV/c^{2}]");

  diff2D_MMnmiss_IMnpip_reactmc = new TH2F("diff2D_MMnmiss_IMnpip_reactmc","diff2D_MMnmiss_IMnpip_reactmc",200,-1.0,1.0,150,-1.5,1.5);
  diff2D_MMnmiss_IMnpip_reactmc->SetXTitle("diff. IMnpip (MCData - Reac.) [GeV/c^{2}]");
  diff2D_MMnmiss_IMnpip_reactmc->SetYTitle("diff. Missing Mass (MCData - Reac.) [GeV/c^{2}]");
   
  diff2D_MMnmiss_IMnpip_wSid_n_reactmc = new TH2F("diff2D_MMnmiss_IMnpip_wSid_n_reactmc","diff2D_MMnmiss_IMnpip_wSid_n_reactmc",200,-1.0,1.0,150,-1.5,1.5);
  diff2D_MMnmiss_IMnpip_wSid_n_reactmc->SetXTitle("diff. IMnpip (MCData - Reac.) [GeV/c^{2}]");
  diff2D_MMnmiss_IMnpip_wSid_n_reactmc->SetYTitle("diff. Missing Mass (MCData - Reac.) [GeV/c^{2}]");

  //diff2D_MMnmiss_IMnpim_recomc_woK0_wSid_n = new TH2F("diff2D_MMnmiss_IMnpim_recomc_woK0_wSid_n","diff2D_MMnmiss_IMnpim_recomc_woK0_wSid_n",200,-1.0,1.0,150,-1.5,1.5);
  //diff2D_MMnmiss_IMnpim_recomc_woK0_wSid_n->SetXTitle("diff. IMnpim (reco. - MCData) [GeV/c^{2}]");
  //diff2D_MMnmiss_IMnpim_recomc_woK0_wSid_n->SetYTitle("diff. Missing Mass (reco. - MCData) [GeV/c^{2}]");

  diff2D_MMnmiss_IMnpim_recomc_woK0_wSid_n_fake1 = new TH2F("diff2D_MMnmiss_IMnpim_recomc_woK0_wSid_n_fake1","diff2D_MMnmiss_IMnpim_recomc_woK0_wSid_n_fake1",200,-1.0,1.0,150,-1.5,1.5);
  diff2D_MMnmiss_IMnpim_recomc_woK0_wSid_n_fake1->SetXTitle("diff. IMnpim (reco. - MCData) [GeV/c^{2}]");
  diff2D_MMnmiss_IMnpim_recomc_woK0_wSid_n_fake1->SetYTitle("diff. Missing Mass (reco. - MCData) [GeV/c^{2}]");
  
  //diff2D_MMnmiss_IMnpip_recomc_woK0_wSid_n = new TH2F("diff2D_MMnmiss_IMnpip_recomc_woK0_wSid_n","diff2D_MMnmiss_IMnpip_recomc_woK0_wSid_n",200,-1.0,1.0,150,-1.5,1.5);
  //diff2D_MMnmiss_IMnpip_recomc_woK0_wSid_n->SetXTitle("diff. IMnpip (reco. - MCData) [GeV/c^{2}]");
  //diff2D_MMnmiss_IMnpip_recomc_woK0_wSid_n->SetYTitle("diff. Missing Mass (reco. - MCData) [GeV/c^{2}]");
  
  diff2D_MMnmiss_IMnpip_recomc_woK0_wSid_n_fake1 = new TH2F("diff2D_MMnmiss_IMnpip_recomc_woK0_wSid_n_fake1","diff2D_MMnmiss_IMnpip_recomc_woK0_wSid_n_fake1",200,-1.0,1.0,150,-1.5,1.5);
  diff2D_MMnmiss_IMnpip_recomc_woK0_wSid_n_fake1->SetXTitle("diff. IMnpip (reco. - MCData) [GeV/c^{2}]");
  diff2D_MMnmiss_IMnpip_recomc_woK0_wSid_n_fake1->SetYTitle("diff. Missing Mass (reco. - MCData) [GeV/c^{2}]");
  
  diff2D_MMnmiss_IMnpim_recomc_wK0_wSid_n = new TH2F("diff2D_MMnmiss_IMnpim_recomc_wK0_wSid_n","diff2D_MMnmiss_IMnpim_recomc_wK0_wSid_n",200,-1.0,1.0,150,-1.5,1.5);
  diff2D_MMnmiss_IMnpim_recomc_wK0_wSid_n->SetXTitle("diff. IMnpim (reco. - MCData) [GeV/c^{2}]");
  diff2D_MMnmiss_IMnpim_recomc_wK0_wSid_n->SetYTitle("diff. Missing Mass (reco. - MCData) [GeV/c^{2}]");
  
  diff2D_MMnmiss_IMnpim_recomc_wK0_wSid_n_fake1 = new TH2F("diff2D_MMnmiss_IMnpim_recomc_wK0_wSid_n_fake1","diff2D_MMnmiss_IMnpim_recomc_wK0_wSid_n_fake1",200,-1.0,1.0,150,-1.5,1.5);
  diff2D_MMnmiss_IMnpim_recomc_wK0_wSid_n_fake1->SetXTitle("diff. IMnpim (reco. - MCData) [GeV/c^{2}]");
  diff2D_MMnmiss_IMnpim_recomc_wK0_wSid_n_fake1->SetYTitle("diff. Missing Mass (reco. - MCData) [GeV/c^{2}]");

  diff2D_MMnmiss_IMnpip_recomc_wK0_wSid_n = new TH2F("diff2D_MMnmiss_IMnpip_recomc_wK0_wSid_n","diff2D_MMnmiss_IMnpip_recomc_wK0_wSid_n",200,-1.0,1.0,150,-1.5,1.5);
  diff2D_MMnmiss_IMnpip_recomc_wK0_wSid_n->SetXTitle("diff. IMnpip (reco. - MCData) [GeV/c^{2}]");
  diff2D_MMnmiss_IMnpip_recomc_wK0_wSid_n->SetYTitle("diff. Missing Mass (reco. - MCData) [GeV/c^{2}]");
  
  diff2D_MMnmiss_IMnpip_recomc_wK0_wSid_n_fake1 = new TH2F("diff2D_MMnmiss_IMnpip_recomc_wK0_wSid_n_fake1","diff2D_MMnmiss_IMnpip_recomc_wK0_wSid_n_fake1",200,-1.0,1.0,150,-1.5,1.5);
  diff2D_MMnmiss_IMnpip_recomc_wK0_wSid_n_fake1->SetXTitle("diff. IMnpip (reco. - MCData) [GeV/c^{2}]");
  diff2D_MMnmiss_IMnpip_recomc_wK0_wSid_n_fake1->SetYTitle("diff. Missing Mass (reco. - MCData) [GeV/c^{2}]");

  diff2D_MMnmiss_IMnpim_recomc_wSid_n = new TH2F("diff2D_MMnmiss_IMnpim_recomc_wSid_n","diff2D_MMnmiss_IMnpim_recomc_wSid_n",200,-1.0,1.0,150,-1.5,1.5);
  diff2D_MMnmiss_IMnpim_recomc_wSid_n->SetXTitle("diff. IMnpim (reco. - MCData) [GeV/c^{2}]");
  diff2D_MMnmiss_IMnpim_recomc_wSid_n->SetYTitle("diff. Missing Mass (reco. - MCData) [GeV/c^{2}]");
  
  diff2D_MMnmiss_IMnpim_recomc_wSid_n_fake1 = new TH2F("diff2D_MMnmiss_IMnpim_recomc_wSid_n_fake1","diff2D_MMnmiss_IMnpim_recomc_wSid_n_fake1",200,-1.0,1.0,150,-1.5,1.5);
  diff2D_MMnmiss_IMnpim_recomc_wSid_n_fake1->SetXTitle("diff. IMnpim (reco. - MCData) [GeV/c^{2}]");
  diff2D_MMnmiss_IMnpim_recomc_wSid_n_fake1->SetYTitle("diff. Missing Mass (reco. - MCData) [GeV/c^{2}]");

  diff2D_MMnmiss_IMnpip_recomc_wSid_n = new TH2F("diff2D_MMnmiss_IMnpip_recomc_wSid_n","diff2D_MMnmiss_IMnpip_recomc_wSid_n",200,-1.0,1.0,150,-1.5,1.5);
  diff2D_MMnmiss_IMnpip_recomc_wSid_n->SetXTitle("diff. IMnpip (reco. - MCData) [GeV/c^{2}]");
  diff2D_MMnmiss_IMnpip_recomc_wSid_n->SetYTitle("diff. Missing Mass (reco. - MCData) [GeV/c^{2}]");
  
  diff2D_MMnmiss_IMnpip_recomc_wSid_n_fake1 = new TH2F("diff2D_MMnmiss_IMnpip_recomc_wSid_n_fake1","diff2D_MMnmiss_IMnpip_recomc_wSid_n_fake1",200,-1.0,1.0,150,-1.5,1.5);
  diff2D_MMnmiss_IMnpip_recomc_wSid_n_fake1->SetXTitle("diff. IMnpip (reco. - MCData) [GeV/c^{2}]");
  diff2D_MMnmiss_IMnpip_recomc_wSid_n_fake1->SetYTitle("diff. Missing Mass (reco. - MCData) [GeV/c^{2}]");

  //diff2D_nmom_IMnpim_recomc_woK0_wSid_n = new TH2F("diff2D_nmom_IMnpim_recomc_woK0_wSid_n","diff2D_nmom_IMnpim_recomc_woK0_wSid_n",200,-1.0,1.0,150,-1.5,1.5);
  //diff2D_nmom_IMnpim_recomc_woK0_wSid_n->SetXTitle("diff. IMnpim (reco. - MCData) [GeV/c^{2}]");
  //diff2D_nmom_IMnpim_recomc_woK0_wSid_n->SetYTitle("diff. n_{CDS} mom. (reco. - MCData) [GeV/c]");
  
  diff2D_nmom_IMnpim_recomc_woK0_wSid_n_fake1 = new TH2F("diff2D_nmom_IMnpim_recomc_woK0_wSid_n_fake1","diff2D_nmom_IMnpim_recomc_woK0_wSid_n_fake1",200,-1.0,1.0,150,-1.5,1.5);
  diff2D_nmom_IMnpim_recomc_woK0_wSid_n_fake1->SetXTitle("diff. IMnpim (reco. - MCData) [GeV/c^{2}]");
  diff2D_nmom_IMnpim_recomc_woK0_wSid_n_fake1->SetYTitle("diff. n_{CDS} mom. (reco. - MCData) [GeV/c]");

  //diff2D_nmom_IMnpip_recomc_woK0_wSid_n = new TH2F("diff2D_nmom_IMnpip_recomc_woK0_wSid_n","diff2D_nmom_IMnpip_recomc_woK0_wSid_n",200,-1.0,1.0,150,-1.5,1.5);
  //diff2D_nmom_IMnpip_recomc_woK0_wSid_n->SetXTitle("diff. IMnpip (reco. - MCData) [GeV/c^{2}]");
  //diff2D_nmom_IMnpip_recomc_woK0_wSid_n->SetYTitle("diff. n_{CDS} mom. (reco. - MCData) [GeV/c]");
  
  diff2D_nmom_IMnpip_recomc_woK0_wSid_n_fake1 = new TH2F("diff2D_nmom_IMnpip_recomc_woK0_wSid_n_fake1","diff2D_nmom_IMnpip_recomc_woK0_wSid_n_fake1",200,-1.0,1.0,150,-1.5,1.5);
  diff2D_nmom_IMnpip_recomc_woK0_wSid_n_fake1->SetXTitle("diff. IMnpip (reco. - MCData) [GeV/c^{2}]");
  diff2D_nmom_IMnpip_recomc_woK0_wSid_n_fake1->SetYTitle("diff. n_{CDS} mom. (reco. - MCData) [GeV/c]");

  diff2D_nmom_IMnpim_recomc_wK0_wSid_n = new TH2F("diff2D_nmom_IMnpim_recomc_wK0_wSid_n","diff2D_nmom_IMnpim_recomc_wK0_wSid_n",200,-1.0,1.0,150,-1.5,1.5);
  diff2D_nmom_IMnpim_recomc_wK0_wSid_n->SetXTitle("diff. IMnpim (reco. - MCData) [GeV/c^{2}]");
  diff2D_nmom_IMnpim_recomc_wK0_wSid_n->SetYTitle("diff. n_{CDS} mom. (reco. - MCData) [GeV/c]");
  
  diff2D_nmom_IMnpim_recomc_wK0_wSid_n_fake1 = new TH2F("diff2D_nmom_IMnpim_recomc_wK0_wSid_n_fake1","diff2D_nmom_IMnpim_recomc_wK0_wSid_n_fake1",200,-1.0,1.0,150,-1.5,1.5);
  diff2D_nmom_IMnpim_recomc_wK0_wSid_n_fake1->SetXTitle("diff. IMnpim (reco. - MCData) [GeV/c^{2}]");
  diff2D_nmom_IMnpim_recomc_wK0_wSid_n_fake1->SetYTitle("diff. n_{CDS} mom. (reco. - MCData) [GeV/c]");

  diff2D_nmom_IMnpip_recomc_wK0_wSid_n = new TH2F("diff2D_nmom_IMnpip_recomc_wK0_wSid_n","diff2D_nmom_IMnpip_recomc_wK0_wSid_n",200,-1.0,1.0,150,-1.5,1.5);
  diff2D_nmom_IMnpip_recomc_wK0_wSid_n->SetXTitle("diff. IMnpip (reco. - MCData) [GeV/c^{2}]");
  diff2D_nmom_IMnpip_recomc_wK0_wSid_n->SetYTitle("diff. n_{CDS} mom. (reco. - MCData) [GeV/c]");
   
  diff2D_nmom_IMnpip_recomc_wK0_wSid_n_fake1 = new TH2F("diff2D_nmom_IMnpip_recomc_wK0_wSid_n_fake1","diff2D_nmom_IMnpip_recomc_wK0_wSid_n_fake1",200,-1.0,1.0,150,-1.5,1.5);
  diff2D_nmom_IMnpip_recomc_wK0_wSid_n_fake1->SetXTitle("diff. IMnpip (reco. - MCData) [GeV/c^{2}]");
  diff2D_nmom_IMnpip_recomc_wK0_wSid_n_fake1->SetYTitle("diff. n_{CDS} mom. (reco. - MCData) [GeV/c]");
  
  diff2D_nmom_IMnpim_recomc_wSid_n = new TH2F("diff2D_nmom_IMnpim_recomc_wSid_n","diff2D_nmom_IMnpim_recomc_wSid_n",200,-1.0,1.0,150,-1.5,1.5);
  diff2D_nmom_IMnpim_recomc_wSid_n->SetXTitle("diff. IMnpim (reco. - MCData) [GeV/c^{2}]");
  diff2D_nmom_IMnpim_recomc_wSid_n->SetYTitle("diff. n_{CDS} mom. (reco. - MCData) [GeV/c]");
  

  diff2D_nmom_IMnpim_recomc_wSid_n_fake1 = new TH2F("diff2D_nmom_IMnpim_recomc_wSid_n_fake1","diff2D_nmom_IMnpim_recomc_wSid_n_fake1",200,-1.0,1.0,150,-1.5,1.5);
  diff2D_nmom_IMnpim_recomc_wSid_n_fake1->SetXTitle("diff. IMnpim (reco. - MCData) [GeV/c^{2}]");
  diff2D_nmom_IMnpim_recomc_wSid_n_fake1->SetYTitle("diff. n_{CDS} mom. (reco. - MCData) [GeV/c]");

  diff2D_nmom_IMnpip_recomc_wSid_n = new TH2F("diff2D_nmom_IMnpip_recomc_wSid_n","diff2D_nmom_IMnpip_recomc_wSid_n",200,-1.0,1.0,150,-1.5,1.5);
  diff2D_nmom_IMnpip_recomc_wSid_n->SetXTitle("diff. IMnpip (reco. - MCData) [GeV/c^{2}]");
  diff2D_nmom_IMnpip_recomc_wSid_n->SetYTitle("diff. n_{CDS} mom. (reco. - MCData) [GeV/c]");

  diffMomnpim_Momnpip_recomc_wSid_n = new TH2F("diffMomnpim_Momnpip_recomc_wSid_n","diffMomnpim_Momnpip_recomc_wSid_n",100,-1.0,1.0,100,-1.0,1.0);
  diffMomnpim_Momnpip_recomc_wSid_n->SetXTitle("diff. Mom(n#pi^{+}) [GeV/c]");
  diffMomnpim_Momnpip_recomc_wSid_n->SetYTitle("diff. Mom(n#pi^{-}) [GeV/c]");

  diff2D_nmom_IMnpip_recomc_wSid_n_fake1 = new TH2F("diff2D_nmom_IMnpip_recomc_wSid_n_fake1","diff2D_nmom_IMnpip_recomc_wSid_n_fake1",200,-1.0,1.0,150,-1.5,1.5);
  diff2D_nmom_IMnpip_recomc_wSid_n_fake1->SetXTitle("diff. IMnpip (reco. - MCData) [GeV/c^{2}]");
  diff2D_nmom_IMnpip_recomc_wSid_n_fake1->SetYTitle("diff. n_{CDS} mom. (reco. - MCData) [GeV/c]");

  diffMMom_recomc_wSid_n = new TH1F("diffMMom_recomc_wSid_n","diffMMom_recomc_wSid_n",1000,-1.0,1.0);
  
  IMnpim_IMnpip_mc = new TH2F("IMnpim_IMnpip_mc", "IMnpim_IMnpip_mc",nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_mc->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_mc->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");

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

  IMnpim_IMnpip_dE = new TH2F(Form("IMnpim_IMnpip_dE"), Form("IMnpim_IMnpip_dE"),nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  
  const int nbinIMnpim_bin = 120;
  const double nbinIMnpim_low = anacuts::Sigmam_center - 30.0*anacuts::Sigmam_sigma;
  //const double nbinIMnpim_low = anacuts::Sigmam_center - 30.0*anacuts::Sigmam_sigma;
  const double nbinIMnpim_high = anacuts::Sigmam_center + 90.0*anacuts::Sigmam_sigma;
  //const double nbinIMnpim_high = anacuts::Sigmam_center + 90.0*anacuts::Sigmam_sigma;
  const int nbinIMnpip_bin = 120;
  const double nbinIMnpip_low = anacuts::Sigmap_center - 26.0*anacuts::Sigmap_sigma;
  const double nbinIMnpip_high = anacuts::Sigmap_center + 94.0*anacuts::Sigmap_sigma;
  
  /*
  const int nbinIMnpim_bin2 = 40;
  const double nbinIMnpim_low2 = anacuts::Sigmam_center - 10.0*3.0*anacuts::Sigmam_sigma;
  //const double nbinIMnpim_low = anacuts::Sigmam_center - 30.0*anacuts::Sigmam_sigma;
  const double nbinIMnpim_high2 = anacuts::Sigmam_center + 10.0*3.0*anacuts::Sigmam_sigma;
  //const double nbinIMnpim_high = anacuts::Sigmam_center + 90.0*anacuts::Sigmam_sigma;
  const int nbinIMnpip_bin2= 40;
  const double nbinIMnpip_low2 = anacuts::Sigmap_center - 27.0*anacuts::Sigmap_sigma;
  const double nbinIMnpip_high2 = anacuts::Sigmap_center + 93.0*anacuts::Sigmap_sigma;
  */

  IMnpim_IMnpip_dE_wSid_n = new TH2F("IMnpim_IMnpip_dE_wSid_n", "IMnpim_IMnpip_dE_wSid_n",
      nbinIMnpip_bin, nbinIMnpip_low, nbinIMnpip_high, nbinIMnpim_bin, nbinIMnpim_low, nbinIMnpim_high);
      //,nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_wSid_n->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_wSid_n->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  
  /*
  for(unsigned int ibin=0;ibin<nbintemplate;ibin++){
    float binlow=1.0+(float)ibin*1./nbintemplate;
    float binhigh=1.0+(float)(ibin+1.0)*1./nbintemplate;
    IMnpim_IMnpip_dE_wSid_n_bin[ibin] = new TH2F(Form("IMnpim_IMnpip_wSid_n_bin%d",ibin),Form("IMnpim_IMnpip_wSid_n %0.2f-%0.2f",binlow,binhigh)
        ,nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
    IMnpim_IMnpip_dE_wSid_n_bin[ibin]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
    IMnpim_IMnpip_dE_wSid_n_bin[ibin]->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  }*/
  
  IMnpim_IMnpip_dE_wK0orwSid_n = new TH2F("IMnpim_IMnpip_dE_wK0orwSid_n", "IMnpim_IMnpip_dE_wK0orwSid_n",
                                 nbinIMnpip_bin, nbinIMnpip_low, nbinIMnpip_high, nbinIMnpim_bin, nbinIMnpim_low, nbinIMnpim_high);
                              //nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_wK0orwSid_n->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_wK0orwSid_n->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  
  for(unsigned int ibin=0;ibin<nbintemplate;ibin++){
    float binlow=1.0+(float)ibin*1./nbintemplate;
    float binhigh=1.0+(float)(ibin+1.0)*1./nbintemplate;
    IMnpim_IMnpip_dE_wK0orwSid_n_bin[ibin] = new TH2F(Form("IMnpim_IMnpip_dE_wK0orwSid_n_bin%d",ibin),Form("IMnpim_IMnpip_dE_wK0orwSid_n %0.2f-%0.2f",binlow,binhigh)
        ,nbinIMnpi/4, 1, 2.0, nbinIMnpi/4, 1, 2.0);
    IMnpim_IMnpip_dE_wK0orwSid_n_bin[ibin]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
    IMnpim_IMnpip_dE_wK0orwSid_n_bin[ibin]->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  }
  
  
  for(unsigned int iwbin=0;iwbin<nwbin;iwbin++){
    IMnpim_IMnpip_dE_wK0orwSid_n_wbin[iwbin] = new TH2F(Form("IMnpim_IMnpip_wK0orwSid_n_wbin%d",iwbin),
                                                        Form("IMnpim_IMnpip_wK0orwSid_n %0.2f-%0.2f",wbinlow[iwbin],wbinhigh[iwbin]),
                                                        nbinIMnpip_wbin, IMnpip_wbinlow, IMnpip_wbinhigh, nbinIMnpim_wbin, IMnpim_wbinlow, IMnpim_wbinhigh);
    IMnpim_IMnpip_dE_wK0orwSid_n_wbin[iwbin]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
    IMnpim_IMnpip_dE_wK0orwSid_n_wbin[iwbin]->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  }
  
  IMnpim_IMnpip_dE_wK0orwSid_n_woSp = new TH2F("IMnpim_IMnpip_dE_wK0orwSid_n_woSp", "IMnpim_IMnpip_dE_wK0orwSid_n_woSp",
                                 nbinIMnpip_bin, nbinIMnpip_low, nbinIMnpip_high, nbinIMnpim_bin, nbinIMnpim_low, nbinIMnpim_high);
                               //nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_wK0orwSid_n_woSp->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_wK0orwSid_n_woSp->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  
  IMnpim_IMnpip_dE_wK0orwSid_n_woSm = new TH2F("IMnpim_IMnpip_dE_wK0orwSid_n_woSm", "IMnpim_IMnpip_dE_wK0orwSid_n_woSm",
                                 nbinIMnpip_bin, nbinIMnpip_low, nbinIMnpip_high, nbinIMnpim_bin, nbinIMnpim_low, nbinIMnpim_high);
                              //nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_wK0orwSid_n_woSm->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_wK0orwSid_n_woSm->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  
  IMnpim_IMnpip_dE_wSid_n_Sp = new TH2F("IMnpim_IMnpip_dE_wSid_n_Sp", "IMnpim_IMnpip_dE_wSid_n_Sp",nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_wSid_n_Sp->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_wSid_n_Sp->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  
  for(unsigned int ibin=0;ibin<nbintemplate;ibin++){
    float binlow=1.0+(float)ibin*1./nbintemplate;
    float binhigh=1.0+(float)(ibin+1.0)*1./nbintemplate;
    IMnpim_IMnpip_dE_wSid_n_Sp_bin[ibin] = new TH2F(Form("IMnpim_IMnpip_wSid_n_Sp_bin%d",ibin),Form("IMnpim_IMnpip_wSid_n_Sp %0.2f-%0.2f",binlow,binhigh)
        ,nbinIMnpi/4, 1, 2.0, nbinIMnpi/4, 1, 2.0);
    IMnpim_IMnpip_dE_wSid_n_Sp_bin[ibin]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
    IMnpim_IMnpip_dE_wSid_n_Sp_bin[ibin]->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  }
  

  for(unsigned int ibin=0;ibin<nwbin;ibin++){
    IMnpim_IMnpip_dE_wSid_n_Sp_wbin[ibin] = new TH2F(Form("IMnpim_IMnpip_wSid_n_Sp_wbin%d",ibin),Form("IMnpim_IMnpip_wSid_n_Sp %0.2f-%0.2f",wbinlow[ibin],wbinhigh[ibin])
        ,nbinIMnpip_wbin, IMnpip_wbinlow , IMnpip_wbinhigh, nbinIMnpim_wbin, IMnpim_wbinlow, IMnpim_wbinhigh);
    IMnpim_IMnpip_dE_wSid_n_Sp_wbin[ibin]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
    IMnpim_IMnpip_dE_wSid_n_Sp_wbin[ibin]->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  }
  
  IMnpim_IMnpip_dE_wSid_n_Sm = new TH2F("IMnpim_IMnpip_dE_wSid_n_Sm", "IMnpim_IMnpip_dE_wSid_n_Sm",nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_wSid_n_Sm->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_wSid_n_Sm->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  
  for(unsigned int ibin=0;ibin<nbintemplate;ibin++){
    float binlow=1.0+(float)ibin*1./nbintemplate;
    float binhigh=1.0+(float)(ibin+1.0)*1./nbintemplate;
    IMnpim_IMnpip_dE_wSid_n_Sm_bin[ibin] = new TH2F(Form("IMnpim_IMnpip_wSid_n_Sm_bin%d",ibin),Form("IMnpim_IMnpip_wSid_n_Sm %0.2f-%0.2f",binlow,binhigh)
        ,nbinIMnpi/4, 1, 2.0, nbinIMnpi/4, 1, 2.0);
    IMnpim_IMnpip_dE_wSid_n_Sm_bin[ibin]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
    IMnpim_IMnpip_dE_wSid_n_Sm_bin[ibin]->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  }
  
  for(unsigned int ibin=0;ibin<nwbin;ibin++){
    IMnpim_IMnpip_dE_wSid_n_Sm_wbin[ibin] = new TH2F(Form("IMnpim_IMnpip_wSid_n_Sm_wbin%d",ibin),Form("IMnpim_IMnpip_wSid_n_Sm %0.2f-%0.2f",wbinlow[ibin],wbinhigh[ibin])
        ,nbinIMnpip_wbin, IMnpip_wbinlow , IMnpip_wbinhigh, nbinIMnpim_wbin, IMnpim_wbinlow, IMnpim_wbinhigh);
    IMnpim_IMnpip_dE_wSid_n_Sm_wbin[ibin]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
    IMnpim_IMnpip_dE_wSid_n_Sm_wbin[ibin]->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  }
  
  IMnpim_IMnpip_dE_wSid_n_fake = new TH2F("IMnpim_IMnpip_dE_wSid_n_fake", "IMnpim_IMnpip_dE_wSid_n_fake",nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_wSid_n_fake->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_wSid_n_fake->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  
  IMnpim_IMnpip_dE_wSid_n_fake_pat2 = new TH2F("IMnpim_IMnpip_dE_wSid_n_fake_pat2", "IMnpim_IMnpip_dE_wSid_n_fake_pat2",nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_wSid_n_fake_pat2->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_wSid_n_fake_pat2->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");

  IMnpim_IMnpip_dE_wSid_n_fake_pat7 = new TH2F("IMnpim_IMnpip_dE_wSid_n_fake_pat7", "IMnpim_IMnpip_dE_wSid_n_fake_pat7",nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_wSid_n_fake_pat7->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_wSid_n_fake_pat7->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");


  IMnpim_IMnpip_dE_woK0 = new TH2F("IMnpim_IMnpip_dE_woK0", "IMnpim_IMnpip_dE_woK0",nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_woK0->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_woK0->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  
  IMnpim_IMnpip_dE_woK0_woSp_vici = new TH2F("IMnpim_IMnpip_dE_woK0_woSp_vici","IMnpim_IMnpip_dE_woK0_woSp_vici",nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_woK0_woSp_vici->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_woK0_woSp_vici->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");

  IMnpim_IMnpip_dE_woK0_woSm_vici = new TH2F("IMnpim_IMnpip_dE_woK0_woSm_vici","IMnpim_IMnpip_dE_woK0_woSm_vici",nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_woK0_woSm_vici->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_woK0_woSm_vici->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  
  IMnpim_IMnpip_dE_woK0_woSp_viciext = new TH2F("IMnpim_IMnpip_dE_woK0_woSp_viciext","IMnpim_IMnpip_dE_woK0_woSp_viciext",nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_woK0_woSp_viciext->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_woK0_woSp_viciext->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");

  IMnpim_IMnpip_dE_woK0_woSm_viciext = new TH2F("IMnpim_IMnpip_dE_woK0_woSm_viciext","IMnpim_IMnpip_dE_woK0_woSm_viciext",nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_woK0_woSm_viciext->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_woK0_woSm_viciext->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");

  IMnpim_IMnpip_dE_wK0 = new TH2F("IMnpim_IMnpip_dE_wK0", "IMnpim_IMnpip_dE_wK0",nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_wK0->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_wK0->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
   

  IMnpim_IMnpip_dE_wK0_woSid_n = new TH2F("IMnpim_IMnpip_dE_wK0_woSid_n","IMnpim_IMnpip_dE_wK0_woSid_n",
      nbinIMnpip_bin, nbinIMnpip_low, nbinIMnpip_high, nbinIMnpim_bin, nbinIMnpim_low, nbinIMnpim_high);
     //nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_wK0_woSid_n->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_wK0_woSid_n->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  
  IMnpim_IMnpip_dE_wK0_woSid_n_45rot = new TH2F("IMnpim_IMnpip_dE_wK0_woSid_n_45rot","IMnpim_IMnpip_dE_wK0_woSid_n_45rot",60, 1.6, 2.2, 100, -0.5,0.5);
  IMnpim_IMnpip_dE_wK0_woSid_n_45rot->SetXTitle("(IM(n#pi^{+})+IM(n#pi^{-}))/#sqrt{2} [GeV/c^{2}]");
  IMnpim_IMnpip_dE_wK0_woSid_n_45rot->SetYTitle("(IM(n#pi^{+})-IM(n#pi^{-}))/#sqrt{2} [GeV/c^{2}]");
  
  IMnpim_IMnpip_dE_wK0_woSid_n_45rot2 = new TH2F("IMnpim_IMnpip_dE_wK0_woSid_n_45rot2","IMnpim_IMnpip_dE_wK0_woSid_n_45rot2",100, -0.5,0.5,60, 1.6, 2.2);
  IMnpim_IMnpip_dE_wK0_woSid_n_45rot2->SetXTitle("(IM(n#pi^{+})-IM(n#pi^{-}))/#sqrt{2} [GeV/c^{2}]");
  IMnpim_IMnpip_dE_wK0_woSid_n_45rot2->SetYTitle("(IM(n#pi^{+})+IM(n#pi^{-}))/#sqrt{2} [GeV/c^{2}]");
  
  IMnpim_IMnpip_dE_wK0_woSid_n_45rot3 = new TH2F("IMnpim_IMnpip_dE_wK0_woSid_n_45rot3","IMnpim_IMnpip_dE_wK0_woSid_n_45rot3",100, -0.5,0.5,120, 0.9, 1.5);
  IMnpim_IMnpip_dE_wK0_woSid_n_45rot3->SetXTitle("(IM(n#pi^{+})-IM(n#pi^{-}))/#sqrt{2} [GeV/c^{2}]");
  //IMnpim_IMnpip_dE_wK0_woSid_n_45rot3->SetYTitle("(IM(n#pi^{+})+IM(n#pi^{-}))/#sqrt{2} [GeV/c^{2}]");
  
  IMnpim_IMnpip_dE_wK0_n_45rot3 = new TH2F("IMnpim_IMnpip_dE_wK0_n_45rot3","IMnpim_IMnpip_dE_wK0_n_45rot3",100, -0.5,0.5,120, 0.9, 1.5);
  IMnpim_IMnpip_dE_wK0_n_45rot3->SetXTitle("(IM(n#pi^{+})-IM(n#pi^{-}))/#sqrt{2} [GeV/c^{2}]");
  //IMnpim_IMnpip_dE_wK0_n_45rot3->SetYTitle("(IM(n#pi^{+})+IM(n#pi^{-}))/#sqrt{2} [GeV/c^{2}]");
  
  for(unsigned int iwbin=0;iwbin<nwbin;iwbin++){
    IMnpim_IMnpip_dE_wK0_woSid_n_wbin[iwbin] = new TH2F(Form("IMnpim_IMnpip_dE_wK0_woSid_n_wbin%d",iwbin),
                                                        Form("IMnpim_IMnpip_dE_wK0_woSid_n %0.2f-%0.2f",wbinlow[iwbin],wbinhigh[iwbin]),
                                                        nbinIMnpim_wbin,IMnpim_wbinlow,IMnpim_wbinhigh,nbinIMnpip_wbin,IMnpip_wbinlow,IMnpip_wbinhigh);
    IMnpim_IMnpip_dE_wK0_woSid_n_wbin[iwbin]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
    IMnpim_IMnpip_dE_wK0_woSid_n_wbin[iwbin]->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  }

  IMnpim_IMnpip_dE_woK0_woSid_n = new TH2F("IMnpim_IMnpip_dE_woK0_woSid_n","IMnpim_IMnpip_dE_woK0_woSid_n",nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_woK0_woSid_n->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_woK0_woSid_n->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  
  MMnmiss_IMpippim_dE = new TH2F("MMnmiss_IMpippim_dE", "MMnmiss_IMpippim_dE",nbinpippim,0.,0.9,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMpippim_dE->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMpippim_dE->SetYTitle("Miss Mass. [GeV/c^{2}]");
  
  MMnmiss_IMpippim_dE_pat2 = new TH2F("MMnmiss_IMpippim_dE_pat2", "MMnmiss_IMpippim_dE_pat2",nbinpippim,0.,0.9,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMpippim_dE_pat2->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMpippim_dE_pat2->SetYTitle("Miss Mass. [GeV/c^{2}]");
  
  MMnmiss_IMpippim_dE_pat7 = new TH2F("MMnmiss_IMpippim_dE_pat7", "MMnmiss_IMpippim_dE_pat7",nbinpippim,0.,0.9,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMpippim_dE_pat7->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMpippim_dE_pat7->SetYTitle("Miss Mass. [GeV/c^{2}]");
  
  MMnmiss_IMpippim_dE_viciSp = new TH2F("MMnmiss_IMpippim_dE_viciSp", "MMnmiss_IMpippim_dE_viciSp",nbinpippim,0.,0.9,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMpippim_dE_viciSp->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMpippim_dE_viciSp->SetYTitle("Miss Mass. [GeV/c^{2}]");
 
  MMnmiss_IMpippim_dE_viciSm = new TH2F("MMnmiss_IMpippim_dE_viciSm", "MMnmiss_IMpippim_dE_viciSm",nbinpippim,0.,0.9,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMpippim_dE_viciSm->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMpippim_dE_viciSm->SetYTitle("Miss Mass. [GeV/c^{2}]");
  
  MMnmiss_IMpippim_dE_viciK0 = new TH2F("MMnmiss_IMpippim_dE_viciK0", "MMnmiss_IMpippim_dE_viciK0",nbinpippim,0.,0.9,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMpippim_dE_viciK0->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMpippim_dE_viciK0->SetYTitle("Miss Mass. [GeV/c^{2}]");

  MMnmiss_IMpippim_dE_wSid = new TH2F("MMnmiss_IMpippim_dE_wSid", "MMnmiss_IMpippim_dE_wSid",nbinpippim,0.,0.9,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMpippim_dE_wSid->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMpippim_dE_wSid->SetYTitle("Miss Mass. [GeV/c^{2}]");
  
  MMnmiss_IMpippim_dE_wSid_fake = new TH2F("MMnmiss_IMpippim_dE_wSid_fake", "MMnmiss_IMpippim_dE_wSid_fake",nbinpippim,0.,0.9,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMpippim_dE_wSid_fake->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMpippim_dE_wSid_fake->SetYTitle("Miss Mass. [GeV/c^{2}]");

  MMnmiss_IMpippim_dE_wSid_n = new TH2F("MMnmiss_IMpippim_dE_wSid_n", "MMnmiss_IMpippim_dE_wSid_n",nbinpippim,0.,0.9,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMpippim_dE_wSid_n->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMpippim_dE_wSid_n->SetYTitle("Miss Mass. [GeV/c^{2}]");
  
  MMnmiss_IMpippim_dE_wSid_n_fake = new TH2F("MMnmiss_IMpippim_dE_wSid_n_fake", "MMnmiss_IMpippim_dE_wSid_n_fake",nbinpippim,0.,0.9,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMpippim_dE_wSid_n_fake->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMpippim_dE_wSid_n_fake->SetYTitle("Miss Mass. [GeV/c^{2}]");
  
  MMnmiss_IMpippim_dE_wSid_n_pat2 = new TH2F("MMnmiss_IMpippim_dE_wSid_n_pat2", "MMnmiss_IMpippim_dE_wSid_n_pat2",nbinpippim,0.,0.9,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMpippim_dE_wSid_n_pat2->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMpippim_dE_wSid_n_pat2->SetYTitle("Miss Mass. [GeV/c^{2}]");
  
  MMnmiss_IMpippim_dE_wSid_n_pat7 = new TH2F("MMnmiss_IMpippim_dE_wSid_n_pat7", "MMnmiss_IMpippim_dE_wSid_n_pat7",nbinpippim,0.,0.9,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMpippim_dE_wSid_n_pat7->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMpippim_dE_wSid_n_pat7->SetYTitle("Miss Mass. [GeV/c^{2}]");



  MMnmiss_IMpippim_dE_woK0 = new TH2F("MMnmiss_IMpippim_dE_woK0", "MMnmiss_IMpippim_dE_woK0",nbinpippim,0.,0.9,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMpippim_dE_woK0->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMpippim_dE_woK0->SetYTitle("Miss Mass. [GeV/c^{2}]");

  MMnmiss_IMpippim_dE_woK0_wSid = new TH2F("MMnmiss_IMpippim_dE_woK0_wSid", "MMnmiss_IMpippim_dE_woK0_wSid",nbinpippim,0.,0.9,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMpippim_dE_woK0_wSid->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMpippim_dE_woK0_wSid->SetYTitle("Miss Mass. [GeV/c^{2}]");
  
  MMnmiss_IMpippim_dE_woK0_wSid_n = new TH2F("MMnmiss_IMpippim_dE_woK0_wSid_n", "MMnmiss_IMpippim_dE_woK0_wSid_n",nbinpippim,0.,0.9,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMpippim_dE_woK0_wSid_n->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMpippim_dE_woK0_wSid_n->SetYTitle("Miss Mass. [GeV/c^{2}]");
  
  MMnmiss_IMpippim_dE_wK0_wSid = new TH2F("MMnmiss_IMpippim_dE_wK0_wSid", "MMnmiss_IMpippim_dE_wK0_wSid",nbinpippim,0.,0.9,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMpippim_dE_wK0_wSid->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMpippim_dE_wK0_wSid->SetYTitle("Miss Mass. [GeV/c^{2}]");
  
  MMnmiss_IMpippim_dE_wK0_wSid_n = new TH2F("MMnmiss_IMpippim_dE_wK0_wSid_n", "MMnmiss_IMpippim_dE_wK0_wSid_n",nbinpippim,0.,0.9,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMpippim_dE_wK0_wSid_n->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMpippim_dE_wK0_wSid_n->SetYTitle("Miss Mass. [GeV/c^{2}]");
   
  Momnpim_Momnpip_dE_wSid_n = new TH2F("Momnpim_Momnpip_dE_wSid_n","Momnpim_Momnpip_dE_wSid_n",nbinmomnpi,0,1.5,nbinmomnpi,0,1.5);
  Momnpim_Momnpip_dE_wSid_n->SetXTitle("Mom(n#pi^{+}) [GeV/c]");
  Momnpim_Momnpip_dE_wSid_n->SetYTitle("Mom(n#pi^{-}) [GeV/c]");

  Momnpim_Momnpip_dE_woK0_wSid_n = new TH2F("Momnpim_Momnpip_dE_woK0_wSid_n","Momnpim_Momnpip_dE_woK0_wSid_n",nbinmomnpi,0,1.5,nbinmomnpi,0,1.5);
  Momnpim_Momnpip_dE_woK0_wSid_n->SetXTitle("Mom(n#pi^{+}) [GeV/c]");
  Momnpim_Momnpip_dE_woK0_wSid_n->SetYTitle("Mom(n#pi^{-}) [GeV/c]");
  
  Momnpim_Momnpip_dE_wK0_wSid_n = new TH2F("Momnpim_Momnpip_dE_wK0_wSid_n","Momnpim_Momnpip_dE_wK0_wSid_n",nbinmomnpi,0,1.5,nbinmomnpi,0,1.5);
  Momnpim_Momnpip_dE_wK0_wSid_n->SetXTitle("Mom(n#pi^{+}) [GeV/c]");
  Momnpim_Momnpip_dE_wK0_wSid_n->SetYTitle("Mom(n#pi^{-}) [GeV/c]");
  
  Momnpim_Mompippim_dE_wSid_n = new TH2F("Momnpim_Mompippim_dE_wSid_n","Momnpim_Mompippim_dE_wSid_n",nbinmompipi,0,1.5,nbinmomnpi,0,1.5);
  Momnpim_Mompippim_dE_wSid_n->SetXTitle("Mom(#pi^{+}#pi^{-}) [GeV/c]");
  Momnpim_Mompippim_dE_wSid_n->SetYTitle("Mom(n#pi^{-}) [GeV/c]");

  Momnpim_Mompippim_dE_woK0_wSid_n = new TH2F("Momnpim_Mompippim_dE_woK0_wSid_n","Momnpim_Mompippim_dE_woK0_wSid_n",nbinmompipi,0,1.5,nbinmomnpi,0,1.5);
  Momnpim_Mompippim_dE_woK0_wSid_n->SetXTitle("Mom(#pi^{+}#pi^{-}) [GeV/c]");
  Momnpim_Mompippim_dE_woK0_wSid_n->SetYTitle("Mom(n#pi^{-}) [GeV/c]");
  
  Momnpim_Mompippim_dE_wK0_wSid_n = new TH2F("Momnpim_Mompippim_dE_wK0_wSid_n","Momnpim_Mompippim_dE_wK0_wSid_n",nbinmompipi,0,1.5,nbinmomnpi,0,1.5);
  Momnpim_Mompippim_dE_wK0_wSid_n->SetXTitle("Mom(#pi^{+}#pi^{-}) [GeV/c]");
  Momnpim_Mompippim_dE_wK0_wSid_n->SetYTitle("Mom(n#pi^{-}) [GeV/c]");
  
  Momnpip_Mompippim_dE_wSid_n = new TH2F("Momnpip_Mompippim_dE_wSid_n","Momnpip_Mompippim_dE_wSid_n",nbinmompipi,0,1.5,nbinmomnpi,0,1.5);
  Momnpip_Mompippim_dE_wSid_n->SetXTitle("Mom(#pi^{+}#pi^{-}) [GeV/c]");
  Momnpip_Mompippim_dE_wSid_n->SetYTitle("Mom(n#pi^{+}) [GeV/c]");

  Momnpip_Mompippim_dE_woK0_wSid_n = new TH2F("Momnpip_Mompippim_dE_woK0_wSid_n","Momnpip_Mompippim_dE_woK0_wSid_n",nbinmompipi,0,1.5,nbinmomnpi,0,1.5);
  Momnpip_Mompippim_dE_woK0_wSid_n->SetXTitle("Mom(#pi^{+}#pi^{-}) [GeV/c]");
  Momnpip_Mompippim_dE_woK0_wSid_n->SetYTitle("Mom(n#pi^{+}) [GeV/c]");
  
  Momnpip_Mompippim_dE_wK0_wSid_n = new TH2F("Momnpip_Mompippim_dE_wK0_wSid_n","Momnpip_Mompippim_dE_wK0_wSid_n",nbinmompipi,0,1.5,nbinmomnpi,0,1.5);
  Momnpip_Mompippim_dE_wK0_wSid_n->SetXTitle("Mom(#pi^{+}#pi^{-}) [GeV/c]");
  Momnpip_Mompippim_dE_wK0_wSid_n->SetYTitle("Mom(n#pi^{+}) [GeV/c]");

  MMnmiss_IMpippim_dE_woK0_woSid = new TH2F("MMnmiss_IMpippim_dE_woK0_woSid", "MMnmiss_IMpippim_dE_woK0_woSid",nbinpippim,0.,0.9,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMpippim_dE_woK0_woSid->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMpippim_dE_woK0_woSid->SetYTitle("Miss Mass. [GeV/c^{2}]");

  MMnmiss_IMpippim_dE_woSid = new TH2F("MMnmiss_IMpippim_dE_woSid", "MMnmiss_IMpippim_dE_woSid",nbinpippim,0.,0.9,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMpippim_dE_woSid->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMpippim_dE_woSid->SetYTitle("Miss Mass. [GeV/c^{2}]");
  
  //MMnmiss_IMpippim_dE_woSid_won = new TH2F("MMnmiss_IMpippim_dE_woSid_won", "MMnmiss_IMpippim_dE_woSid_won",nbinpippim,0.,0.9,nbinnmiss, nmisslow, nmisshigh);
  //MMnmiss_IMpippim_dE_woSid_won->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  //MMnmiss_IMpippim_dE_woSid_won->SetYTitle("Miss Mass. [GeV/c^{2}]");
  
  MMnmiss_IMpippim_dE_woK0_woSid_won = new TH2F("MMnmiss_IMpippim_dE_woK0_woSid_won", "MMnmiss_IMpippim_dE_woK0_woSid_won",nbinpippim,0.,0.9,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMpippim_dE_woK0_woSid_won->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMpippim_dE_woK0_woSid_won->SetYTitle("Miss Mass. [GeV/c^{2}]");

  MMnmiss_IMpippim_dE_wK0_woSid_won = new TH2F("MMnmiss_IMpippim_dE_wK0_woSid_won", "MMnmiss_IMpippim_dE_wK0_woSid_won",nbinpippim,0.,0.9,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMpippim_dE_wK0_woSid_won->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMpippim_dE_wK0_woSid_won->SetYTitle("Miss Mass. [GeV/c^{2}]");

  for(unsigned int iwbin=0;iwbin<nwbin;iwbin++){
    MMnmiss_IMpippim_dE_wK0_woSid_n_wbin[iwbin] = new TH2F(Form("MMnmiss_IMpippim_wK0_woSid_n_wbin%d",iwbin)
                                                          ,Form("MMnmiss_IMpippim_wK0_woSid_n %0.2f-%0.2f",wbinlow[iwbin],wbinhigh[iwbin])
                                                          ,nbinpippim,0.,0.9,nbinnmiss, nmisslow, nmisshigh);
    MMnmiss_IMpippim_dE_wK0_woSid_n_wbin[iwbin]->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
    MMnmiss_IMpippim_dE_wK0_woSid_n_wbin[iwbin]->SetYTitle("Miss Mass. [GeV/c^{2}]");
  }
  
  MMnmiss_IMpippim_dE_wK0_woSidn_won = new TH2F("MMnmiss_IMpippim_dE_wK0_woSidn_won", "MMnmiss_IMpippim_dE_wK0_woSidn_won",nbinpippim,0.,0.9,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMpippim_dE_wK0_woSidn_won->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMpippim_dE_wK0_woSidn_won->SetYTitle("Miss Mass. [GeV/c^{2}]");

  //MMnmiss_Mompippim_dE_woK0_woSidn = new TH2F("MMnmiss_Mompippim_dE_woK0_woSidn", "MMnmiss_Mompippim_dE_woK0_woSidn",nbinmompipi,0.,1.5,nbinnmiss, nmisslow, nmisshigh);
  //MMnmiss_Mompippim_dE_woK0_woSidn->SetXTitle("Mom.(#pi^{+}#pi^{-}) [GeV/c]");
  //MMnmiss_Mompippim_dE_woK0_woSidn->SetYTitle("Miss Mass. [GeV/c^{2}]");

  MMnmiss_Mompippim_dE_woK0_woSid_won = new TH2F("MMnmiss_Mompippim_dE_woK0_woSid_won", "MMnmiss_Mompippim_dE_woK0_woSid_won",nbinmompipi,0.,1.5,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_Mompippim_dE_woK0_woSid_won->SetXTitle("Mom.(#pi^{+}#pi^{-}) [GeV/c]");
  MMnmiss_Mompippim_dE_woK0_woSid_won->SetYTitle("Miss Mass. [GeV/c^{2}]");

  MMnmiss_Mompippim_dE_wK0_woSid_won = new TH2F("MMnmiss_Mompippim_dE_wK0_woSid_won", "MMnmiss_Mompippim_dE_wK0_woSid_won",nbinmompipi,0.,1.5,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_Mompippim_dE_wK0_woSid_won->SetXTitle("Mom.(#pi^{+}#pi^{-}) [GeV/c]");
  MMnmiss_Mompippim_dE_wK0_woSid_won->SetYTitle("Miss Mass. [GeV/c^{2}]");

  MMnmiss_Mompippim_dE_wK0_woSidn_won = new TH2F("MMnmiss_Mompippim_dE_wK0_woSidn_won", "MMnmiss_Mompippim_dE_wK0_woSidn_won",nbinmompipi,0.,1.5,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_Mompippim_dE_wK0_woSidn_won->SetXTitle("Mom.(#pi^{+}#pi^{-}) [GeV/c]");
  MMnmiss_Mompippim_dE_wK0_woSidn_won->SetYTitle("Miss Mass. [GeV/c^{2}]");

  MMnmiss_IMnpip_dE = new TH2F("MMnmiss_IMnpip_dE", "MMnmiss_IMnpip_dE",nbinIMnpi,1.,2.0,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMnpip_dE->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  MMnmiss_IMnpip_dE->SetYTitle("Miss Mass. [GeV/c^{2}]");
  
  MMnmiss_IMnpip_dE_fake = new TH2F("MMnmiss_IMnpip_dE_fake", "MMnmiss_IMnpip_dE_fake",nbinIMnpi,1.,2.0,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMnpip_dE_fake->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  MMnmiss_IMnpip_dE_fake->SetYTitle("Miss Mass. [GeV/c^{2}]");
  
  MMnmiss_IMnpip_dE_woSid = new TH2F("MMnmiss_IMnpip_dE_woSid", "MMnmiss_IMnpip_dE_woSid",nbinIMnpi,1.,2.0,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMnpip_dE_woSid->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  MMnmiss_IMnpip_dE_woSid->SetYTitle("Miss Mass. [GeV/c^{2}]");
  
  //MMnmiss_IMnpip_dE_woSid_won = new TH2F("MMnmiss_IMnpip_dE_woSid_won", "MMnmiss_IMnpip_dE_woSid_won",nbinIMnpi,1.,2.0,nbinnmiss, nmisslow, nmisshigh);
  //MMnmiss_IMnpip_dE_woSid_won->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  //MMnmiss_IMnpip_dE_woSid_won->SetYTitle("Miss Mass. [GeV/c^{2}]");

  MMnmiss_IMnpip_dE_woK0 = new TH2F("MMnmiss_IMnpip_dE_woK0", "MMnmiss_IMnpip_dE_woK0",nbinIMnpi,1.,2.0,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMnpip_dE_woK0->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  MMnmiss_IMnpip_dE_woK0->SetYTitle("Miss Mass. [GeV/c^{2}]");

  MMnmiss_IMnpip_dE_woK0_woSm = new TH2F("MMnmiss_IMnpip_dE_woK0_woSm", "MMnmiss_IMnpip_dE_woK0_woSm",nbinIMnpi,1.,2.0,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMnpip_dE_woK0_woSm->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  MMnmiss_IMnpip_dE_woK0_woSm->SetYTitle("Miss Mass. [GeV/c^{2}]");

  for(unsigned int iwbin=0;iwbin<nwbin;iwbin++){
    MMnmiss_IMnpip_dE_woK0_wSid_n_woSm_wbin[iwbin] = new TH2F(Form("MMnmiss_IMnpip_woK0_wSid_n_woSm_wbin%d",iwbin)
                                                      ,Form("MMnmiss_IMnpip_woK0_wSid_n_woSm %0.2f-%0.2f",wbinlow[iwbin],wbinhigh[iwbin])
                                                      ,nbinIMnpi,1.,2.0,nbinnmiss, nmisslow, nmisshigh);
    MMnmiss_IMnpip_dE_woK0_wSid_n_woSm_wbin[iwbin]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
    MMnmiss_IMnpip_dE_woK0_wSid_n_woSm_wbin[iwbin]->SetYTitle("Miss Mass. [GeV/c^{2}]");
  }


  MMnmiss_IMnpip_dE_woK0_woSm_vici = new TH2F("MMnmiss_IMnpip_dE_woK0_woSm_vici", "MMnmiss_IMnpip_dE_woK0_woSm_vici",nbinIMnpi,1.,2.0,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMnpip_dE_woK0_woSm_vici->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  MMnmiss_IMnpip_dE_woK0_woSm_vici->SetYTitle("Miss Mass. [GeV/c^{2}]");

  MMnmiss_IMnpip_dE_woK0_woSm_viciext = new TH2F("MMnmiss_IMnpip_dE_woK0_woSm_viciext", "MMnmiss_IMnpip_dE_woK0_woSm_viciext",nbinIMnpi,1.,2.0,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMnpip_dE_woK0_woSm_viciext->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  MMnmiss_IMnpip_dE_woK0_woSm_viciext->SetYTitle("Miss Mass. [GeV/c^{2}]");


  MMnmiss_IMnpip_dE_wK0_woSm = new TH2F("MMnmiss_IMnpip_dE_wK0_woSm", "MMnmiss_IMnpip_dE_wK0_woSm",nbinIMnpi,1.,2.0,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMnpip_dE_wK0_woSm->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  MMnmiss_IMnpip_dE_wK0_woSm->SetYTitle("Miss Mass. [GeV/c^{2}]");

  MMnmiss_IMnpip_dE_woK0_woSm_cross = new TH2F("MMnmiss_IMnpip_dE_woK0_woSm_cross", "MMnmiss_IMnpip_dE_woK0_woSm_cross",nbinIMnpi,1.,2.0,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMnpip_dE_woK0_woSm_cross->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  MMnmiss_IMnpip_dE_woK0_woSm_cross->SetYTitle("Miss Mass. [GeV/c^{2}]");

  MMnmiss_IMnpip_dE_woK0_woSm_n = new TH2F("MMnmiss_IMnpip_dE_woK0_woSm_n", "MMnmiss_IMnpip_dE_woK0_woSm_n",nbinIMnpi,1.,2.0,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMnpip_dE_woK0_woSm_n->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  MMnmiss_IMnpip_dE_woK0_woSm_n->SetYTitle("Miss Mass. [GeV/c^{2}]");
  
  MMnmiss_IMnpip_dE_wK0_woSm_n = new TH2F("MMnmiss_IMnpip_dE_wK0_woSm_n", "MMnmiss_IMnpip_dE_wK0_woSm_n",nbinIMnpi,1.,2.0,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMnpip_dE_wK0_woSm_n->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  MMnmiss_IMnpip_dE_wK0_woSm_n->SetYTitle("Miss Mass. [GeV/c^{2}]");

  MMnmiss_IMnpip_dE_woK0_woSid = new TH2F("MMnmiss_IMnpip_dE_woK0_woSid", "MMnmiss_IMnpip_dE_woK0_woSid",nbinIMnpi,1.,2.0,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMnpip_dE_woK0_woSid->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  MMnmiss_IMnpip_dE_woK0_woSid->SetYTitle("Miss Mass. [GeV/c^{2}]");


  //MMnmiss_IMnpip_dE_woK0_woSidn = new TH2F("MMnmiss_IMnpip_dE_woK0_woSidn", "MMnmiss_IMnpip_dE_woK0_woSidn",nbinIMnpi,1.,2.0,nbinnmiss, nmisslow, nmisshigh);
  //MMnmiss_IMnpip_dE_woK0_woSidn->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  //MMnmiss_IMnpip_dE_woK0_woSidn->SetYTitle("Miss Mass. [GeV/c^{2}]");

  MMnmiss_IMnpip_dE_woK0_woSid_won = new TH2F("MMnmiss_IMnpip_dE_woK0_woSid_won", "MMnmiss_IMnpip_dE_woK0_woSid_won",nbinIMnpi,1.,2.0,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMnpip_dE_woK0_woSid_won->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  MMnmiss_IMnpip_dE_woK0_woSid_won->SetYTitle("Miss Mass. [GeV/c^{2}]");
  
  MMnmiss_Momnpip_dE_woK0_woSid_won = new TH2F("MMnmiss_Momnpip_dE_woK0_woSid_won", "MMnmiss_Momnpip_dE_woK0_woSid_won",300,0.,3.0,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_Momnpip_dE_woK0_woSid_won->SetXTitle("Mom(n#pi^{+}) [GeV/c]");
  MMnmiss_Momnpip_dE_woK0_woSid_won->SetYTitle("Miss Mass. [GeV/c^{2}]");

  //MMnmiss_IMnpip_dE_woK0_woSidn_cross = new TH2F("MMnmiss_IMnpip_dE_woK0_woSidn_cross", "MMnmiss_IMnpip_dE_woK0_woSidn_cross",nbinIMnpi,1.,2.0,nbinnmiss, nmisslow, nmisshigh);
  //MMnmiss_IMnpip_dE_woK0_woSidn_cross->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  //MMnmiss_IMnpip_dE_woK0_woSidn_cross->SetYTitle("Miss Mass. [GeV/c^{2}]");

  //MMnmiss_IMnpip_dE_woK0_woSidn_cross_Sp = new TH2F("MMnmiss_IMnpip_dE_woK0_woSidn_cross_Sp", "MMnmiss_IMnpip_dE_woK0_woSidn_cross_Sp",nbinIMnpi,1.,2.0,nbinnmiss, nmisslow, nmisshigh);
  //MMnmiss_IMnpip_dE_woK0_woSidn_cross_Sp->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  //MMnmiss_IMnpip_dE_woK0_woSidn_cross_Sp->SetYTitle("Miss Mass. [GeV/c^{2}]");

  MMnmiss_IMnpip_dE_wK0_woSid_won = new TH2F("MMnmiss_IMnpip_dE_wK0_woSid_won", "MMnmiss_IMnpip_dE_wK0_woSid_won",nbinIMnpi,1.,2.0,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMnpip_dE_wK0_woSid_won->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  MMnmiss_IMnpip_dE_wK0_woSid_won->SetYTitle("Miss Mass. [GeV/c^{2}]");
  
  MMnmiss_Momnpip_dE_wK0_woSid_won = new TH2F("MMnmiss_Momnpip_dE_wK0_woSid_won", "MMnmiss_Momnpip_dE_wK0_woSid_won",300,0.,3.0,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_Momnpip_dE_wK0_woSid_won->SetXTitle("Mom(n#pi^{+}) [GeV/c]");
  MMnmiss_Momnpip_dE_wK0_woSid_won->SetYTitle("Miss Mass. [GeV/c^{2}]");

  MMnmiss_IMnpip_dE_wK0_woSidn_won = new TH2F("MMnmiss_IMnpip_dE_wK0_woSidn_won", "MMnmiss_IMnpip_dE_wK0_woSidn_won",nbinIMnpi,1.,2.0,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMnpip_dE_wK0_woSidn_won->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  MMnmiss_IMnpip_dE_wK0_woSidn_won->SetYTitle("Miss Mass. [GeV/c^{2}]");

  MMnmiss_IMnpim_dE = new TH2F("MMnmiss_IMnpim_dE", "MMnmiss_IMnpim_dE",nbinIMnpi,1.,2.0,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMnpim_dE->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMnpim_dE->SetYTitle("Miss Mass. [GeV/c^{2}]");
  
  MMnmiss_IMnpim_dE_fake = new TH2F("MMnmiss_IMnpim_dE_fake", "MMnmiss_IMnpim_dE_fake",nbinIMnpi,1.,2.0,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMnpim_dE_fake->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMnpim_dE_fake->SetYTitle("Miss Mass. [GeV/c^{2}]");

  MMnmiss_IMnpim_dE_woSid = new TH2F("MMnmiss_IMnpim_dE_woSid", "MMnmiss_IMnpim_dE_woSid",nbinIMnpi,1.,2.0,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMnpim_dE_woSid->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMnpim_dE_woSid->SetYTitle("Miss Mass. [GeV/c^{2}]");
  
  //MMnmiss_IMnpim_dE_woSid_won = new TH2F("MMnmiss_IMnpim_dE_woSid_won", "MMnmiss_IMnpim_dE_woSid_won",nbinIMnpi,1.,2.0,nbinnmiss, nmisslow, nmisshigh);
  //MMnmiss_IMnpim_dE_woSid_won->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  //MMnmiss_IMnpim_dE_woSid_won->SetYTitle("Miss Mass. [GeV/c^{2}]");
  
  MMnmiss_IMnpim_dE_woK0 = new TH2F("MMnmiss_IMnpim_dE_woK0", "MMnmiss_IMnpim_dE_woK0",nbinIMnpi,1.,2.0,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMnpim_dE_woK0->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMnpim_dE_woK0->SetYTitle("Miss Mass. [GeV/c^{2}]");

  MMnmiss_IMnpim_dE_woK0_woSp = new TH2F("MMnmiss_IMnpim_dE_woK0_woSp", "MMnmiss_IMnpim_dE_woK0_woSp",nbinIMnpi,1.,2.0,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMnpim_dE_woK0_woSp->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMnpim_dE_woK0_woSp->SetYTitle("Miss Mass. [GeV/c^{2}]");
  
  for(unsigned int iwbin=0;iwbin<nwbin;iwbin++){
    MMnmiss_IMnpim_dE_woK0_wSid_n_woSp_wbin[iwbin] = new TH2F(Form("MMnmiss_IMnpim_woK0_wSid_n_woSp_wbin%d",iwbin)
                                                      ,Form("MMnmiss_IMnpim_woK0_wSid_n_woSp %0.2f-%0.2f",wbinlow[iwbin],wbinhigh[iwbin])
                                                      ,nbinIMnpi,1.,2.0,nbinnmiss, nmisslow, nmisshigh);
    MMnmiss_IMnpim_dE_woK0_wSid_n_woSp_wbin[iwbin]->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
    MMnmiss_IMnpim_dE_woK0_wSid_n_woSp_wbin[iwbin]->SetYTitle("Miss Mass. [GeV/c^{2}]");
  }
  
  MMnmiss_IMnpim_dE_woK0_woSp_vici = new TH2F("MMnmiss_IMnpim_dE_woK0_woSp_vici", "MMnmiss_IMnpim_dE_woK0_woSp_vici",nbinIMnpi,1.,2.0,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMnpim_dE_woK0_woSp_vici->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMnpim_dE_woK0_woSp_vici->SetYTitle("Miss Mass. [GeV/c^{2}]");
  
  MMnmiss_IMnpim_dE_woK0_woSp_viciext = new TH2F("MMnmiss_IMnpim_dE_woK0_woSp_viciext", "MMnmiss_IMnpim_dE_woK0_woSp_viciext",nbinIMnpi,1.,2.0,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMnpim_dE_woK0_woSp_viciext->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMnpim_dE_woK0_woSp_viciext->SetYTitle("Miss Mass. [GeV/c^{2}]");

  MMnmiss_IMnpim_dE_wK0_woSp = new TH2F("MMnmiss_IMnpim_dE_wK0_woSp", "MMnmiss_IMnpim_dE_wK0_woSp",nbinIMnpi,1.,2.0,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMnpim_dE_wK0_woSp->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMnpim_dE_wK0_woSp->SetYTitle("Miss Mass. [GeV/c^{2}]");

  MMnmiss_IMnpim_dE_woK0_woSp_cross = new TH2F("MMnmiss_IMnpim_dE_woK0_woSp_cross", "MMnmiss_IMnpim_dE_woK0_woSp_cross",nbinIMnpi,1.,2.0,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMnpim_dE_woK0_woSp_cross->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMnpim_dE_woK0_woSp_cross->SetYTitle("Miss Mass. [GeV/c^{2}]");

  MMnmiss_IMnpim_dE_woK0_woSp_n = new TH2F("MMnmiss_IMnpim_dE_woK0_woSp_n", "MMnmiss_IMnpim_dE_woK0_woSp_n",nbinIMnpi,1.,2.0,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMnpim_dE_woK0_woSp_n->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMnpim_dE_woK0_woSp_n->SetYTitle("Miss Mass. [GeV/c^{2}]");
  
  MMnmiss_IMnpim_dE_wK0_woSp_n = new TH2F("MMnmiss_IMnpim_dE_wK0_woSp_n", "MMnmiss_IMnpim_dE_wK0_woSp_n",nbinIMnpi,1.,2.0,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMnpim_dE_wK0_woSp_n->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMnpim_dE_wK0_woSp_n->SetYTitle("Miss Mass. [GeV/c^{2}]");

  MMnmiss_IMnpim_dE_wK0_woSid_won = new TH2F("MMnmiss_IMnpim_dE_wK0_woSid_won", "MMnmiss_IMnpim_dE_wK0_woSid_won",nbinIMnpi,1.,2.0,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMnpim_dE_wK0_woSid_won->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMnpim_dE_wK0_woSid_won->SetYTitle("Miss Mass. [GeV/c^{2}]");
  
  MMnmiss_Momnpim_dE_wK0_woSid_won = new TH2F("MMnmiss_Momnpim_dE_wK0_woSid_won", "MMnmiss_Momnpim_dE_wK0_woSid_won",300,0.,3.0,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_Momnpim_dE_wK0_woSid_won->SetXTitle("Mom(n#pi^{-}) [GeV/c]");
  MMnmiss_Momnpim_dE_wK0_woSid_won->SetYTitle("Miss Mass. [GeV/c^{2}]");

  MMnmiss_IMnpim_dE_wK0_woSidn_won = new TH2F("MMnmiss_IMnpim_dE_wK0_woSidn_won", "MMnmiss_IMnpim_dE_wK0_woSidn_won",nbinIMnpi,1.,2.0,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMnpim_dE_wK0_woSidn_won->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMnpim_dE_wK0_woSidn_won->SetYTitle("Miss Mass. [GeV/c^{2}]");
 
  //Momnpim_Momnpip_dE_woSid_won = new TH2F("Momnpim_Momnpip_dE_woSid_won","Momnpim_Momnpip_dE_woSid_won",nbinmomnpi,0,1.5,nbinmomnpi,0,1.5);
  //Momnpim_Momnpip_dE_woSid_won->SetXTitle("Mom(n#pi^{+}) [GeV/c]");
  //Momnpim_Momnpip_dE_woSid_won->SetYTitle("Mom(n#pi^{-}) [GeV/c]");

  Momnpim_Momnpip_dE_woK0_woSid_won = new TH2F("Momnpim_Momnpip_dE_woK0_woSid_won","Momnpim_Momnpip_dE_woK0_woSid_won",nbinmomnpi,0,1.5,nbinmomnpi,0,1.5);
  Momnpim_Momnpip_dE_woK0_woSid_won->SetXTitle("Mom(n#pi^{+}) [GeV/c]");
  Momnpim_Momnpip_dE_woK0_woSid_won->SetYTitle("Mom(n#pi^{-}) [GeV/c]");

  Momnpim_Momnpip_dE_wK0_woSid_won = new TH2F("Momnpim_Momnpip_dE_wK0_woSid_won","Momnpim_Momnpip_dE_wK0_woSid_won",nbinmomnpi,0,1.5,nbinmomnpi,0,1.5);
  Momnpim_Momnpip_dE_wK0_woSid_won->SetXTitle("Mom(n#pi^{+}) [GeV/c]");
  Momnpim_Momnpip_dE_wK0_woSid_won->SetYTitle("Mom(n#pi^{-}) [GeV/c]");

  //Momnpim_Mompippim_dE_woSid_won = new TH2F("Momnpim_Mompippim_dE_woSid_won","Momnpim_Mompippim_dE_woSid_won",nbinmompipi,0,1.5,nbinmomnpi,0,1.5);
  //Momnpim_Mompippim_dE_woSid_won->SetXTitle("Mom(#pi^{+}#pi^{-}) [GeV/c]");
  //Momnpim_Mompippim_dE_woSid_won->SetYTitle("Mom(n#pi^{-}) [GeV/c]");
  
  Momnpim_Mompippim_dE_woK0_woSid_won = new TH2F("Momnpim_Mompippim_dE_woK0_woSid_won","Momnpim_Mompippim_dE_woK0_woSid_won",nbinmompipi,0,1.5,nbinmomnpi,0,1.5);
  Momnpim_Mompippim_dE_woK0_woSid_won->SetXTitle("Mom(#pi^{+}#pi^{-}) [GeV/c]");
  Momnpim_Mompippim_dE_woK0_woSid_won->SetYTitle("Mom(n#pi^{-}) [GeV/c]");

  Momnpim_Mompippim_dE_wK0_woSid_won = new TH2F("Momnpim_Mompippim_dE_wK0_woSid_won","Momnpim_Mompippim_dE_wK0_woSid_won",nbinmompipi,0,1.5,nbinmomnpi,0,1.5);
  Momnpim_Mompippim_dE_wK0_woSid_won->SetXTitle("Mom(#pi^{+}#pi^{-}) [GeV/c]");
  Momnpim_Mompippim_dE_wK0_woSid_won->SetYTitle("Mom(n#pi^{-}) [GeV/c]");

  Momnpip_Mompippim_dE_woK0_woSid_won = new TH2F("Momnpip_Mompippim_dE_woK0_woSid_won","Momnpip_Mompippim_dE_woK0_woSid_won",nbinmompipi,0,1.5,nbinmomnpi,0,1.5);
  Momnpip_Mompippim_dE_woK0_woSid_won->SetXTitle("Mom(#pi^{+}#pi^{-}) [GeV/c]");
  Momnpip_Mompippim_dE_woK0_woSid_won->SetYTitle("Mom(n#pi^{+}) [GeV/c]");

  Momnpip_Mompippim_dE_wK0_woSid_won = new TH2F("Momnpip_Mompippim_dE_wK0_woSid_won","Momnpip_Mompippim_dE_wK0_woSid_won",nbinmompipi,0,1.5,nbinmomnpi,0,1.5);
  Momnpip_Mompippim_dE_wK0_woSid_won->SetXTitle("Mom(#pi^{+}#pi^{-}) [GeV/c]");
  Momnpip_Mompippim_dE_wK0_woSid_won->SetYTitle("Mom(n#pi^{+}) [GeV/c]");

  MMnmiss_IMnpim_dE_woK0_woSid = new TH2F("MMnmiss_IMnpim_dE_woK0_woSid", "MMnmiss_IMnpim_dE_woK0_woSid",nbinIMnpi,1.,2.0,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMnpim_dE_woK0_woSid->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMnpim_dE_woK0_woSid->SetYTitle("Miss Mass. [GeV/c^{2}]");

  //MMnmiss_IMnpim_dE_woK0_woSidn = new TH2F("MMnmiss_IMnpim_dE_woK0_woSidn", "MMnmiss_IMnpim_dE_woK0_woSidn",nbinIMnpi,1.,2.0,nbinnmiss, nmisslow, nmisshigh);
  //MMnmiss_IMnpim_dE_woK0_woSidn->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  //MMnmiss_IMnpim_dE_woK0_woSidn->SetYTitle("Miss Mass. [GeV/c^{2}]");

  MMnmiss_IMnpim_dE_woK0_woSid_won = new TH2F("MMnmiss_IMnpim_dE_woK0_woSid_won", "MMnmiss_IMnpim_dE_woK0_woSid_won",nbinIMnpi,1.,2.0,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMnpim_dE_woK0_woSid_won->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMnpim_dE_woK0_woSid_won->SetYTitle("Miss Mass. [GeV/c^{2}]");
  
  MMnmiss_Momnpim_dE_woK0_woSid_won = new TH2F("MMnmiss_Momnpim_dE_woK0_woSid_won", "MMnmiss_Momnpim_dE_woK0_woSid_won",300,0.,3.0,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_Momnpim_dE_woK0_woSid_won->SetXTitle("Mom(n#pi^{-}) [GeV/c]");
  MMnmiss_Momnpim_dE_woK0_woSid_won->SetYTitle("Miss Mass. [GeV/c^{2}]");

  //MMnmiss_IMnpim_dE_woK0_woSidn_cross = new TH2F("MMnmiss_IMnpim_dE_woK0_woSidn_cross", "MMnmiss_IMnpim_dE_woK0_woSidn_cross",nbinIMnpi,1.,2.0,nbinnmiss, nmisslow, nmisshigh);
  //MMnmiss_IMnpim_dE_woK0_woSidn_cross->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  //MMnmiss_IMnpim_dE_woK0_woSidn_cross->SetYTitle("Miss Mass. [GeV/c^{2}]");

  //MMnmiss_IMnpim_dE_woK0_woSidn_cross_Sm = new TH2F("MMnmiss_IMnpim_dE_woK0_woSidn_cross_Sm", "MMnmiss_IMnpim_dE_woK0_woSidn_cross_Sm",nbinIMnpi,1.,2.0,nbinnmiss, nmisslow, nmisshigh);
  //MMnmiss_IMnpim_dE_woK0_woSidn_cross_Sm->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  //MMnmiss_IMnpim_dE_woK0_woSidn_cross_Sm->SetYTitle("Miss Mass. [GeV/c^{2}]");

  //MMnmiss_IMpippim_dE_woK0_woSidn = new TH2F("MMnmiss_IMpippim_dE_woK0_woSidn", "MMnmiss_IMpippim_dE_woK0_woSidn",nbinpippim,0.,0.9,nbinnmiss, nmisslow, nmisshigh);
  //MMnmiss_IMpippim_dE_woK0_woSidn->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  //MMnmiss_IMpippim_dE_woK0_woSidn->SetYTitle("Miss Mass. [GeV/c^{2}]");

  IMnpim_IMnpip_dE_n = new TH2F("IMnpim_IMnpip_dE_n","IMnpim_IMnpip_dE_n",nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_n->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_n->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  
  /*
  for(unsigned int ibin=0;ibin<nbintemplate;ibin++){
    float binlow=1.0+(float)ibin*1./nbintemplate;
    float binhigh=1.0+(float)(ibin+1.0)*1./nbintemplate;
    IMnpim_IMnpip_dE_n_bin[ibin] = new TH2F(Form("IMnpim_IMnpip_dE_n_bin%d",ibin),Form("IMnpim_IMnpip_dE_n %0.2f-%0.2f",binlow,binhigh)
          ,nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
    IMnpim_IMnpip_dE_n_bin[ibin]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
    IMnpim_IMnpip_dE_n_bin[ibin]->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  }
  */

  
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
  
  IMnpim_IMnpip_dE_n_fake = new TH2F("IMnpim_IMnpip_dE_n_fake","IMnpim_IMnpip_dE_n_fake",nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_n_fake->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_n_fake->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");

  IMnpim_IMnpip_dE_woK0_n = new TH2F(Form("IMnpim_IMnpip_dE_woK0_n"),Form("IMnpim_IMnpip_dE_woK0_n"),nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_woK0_n->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_woK0_n->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  
  IMnpim_IMnpip_dE_woK0_wSid_n = new TH2F(Form("IMnpim_IMnpip_dE_woK0_wSid_n"),Form("IMnpim_IMnpip_dE_woK0_wSid_n"),nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_woK0_wSid_n->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_woK0_wSid_n->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  
  IMnpim_IMnpip_dE_woK0_wSid_n_woSm = new TH2F("IMnpim_IMnpip_dE_woK0_wSid_n_woSm","IMnpim_IMnpip_dE_woK0_wSid_n_woSm",
                                     //,nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
                                    nbinIMnpip_bin, nbinIMnpip_low, nbinIMnpip_high, nbinIMnpim_bin, nbinIMnpim_low, nbinIMnpim_high);
  IMnpim_IMnpip_dE_woK0_wSid_n_woSm->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_woK0_wSid_n_woSm->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
   
  for(unsigned int iwbin=0;iwbin<nwbin;iwbin++){
    IMnpim_IMnpip_dE_woK0_wSid_n_woSm_wbin[iwbin] = new TH2F(
        Form("IMnpim_IMnpip_woK0_wSid_n_woSm_wbin%d",iwbin),
        Form("IMnpim_IMnpip_woK0_wSid_n_woSm %0.2f-%0.2f",wbinlow[iwbin],wbinhigh[iwbin]),
          nbinIMnpip_wbin, IMnpip_wbinlow, IMnpip_wbinhigh, nbinIMnpim_wbin, IMnpim_wbinlow, IMnpim_wbinhigh);
    IMnpim_IMnpip_dE_woK0_wSid_n_woSm_wbin[iwbin]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
    IMnpim_IMnpip_dE_woK0_wSid_n_woSm_wbin[iwbin]->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  }

  IMnpim_IMnpip_dE_woK0_wSid_n_woSp = new TH2F("IMnpim_IMnpip_dE_woK0_wSid_n_woSp","IMnpim_IMnpip_dE_woK0_wSid_n_woSp",
                                    nbinIMnpip_bin, nbinIMnpip_low, nbinIMnpip_high, nbinIMnpim_bin, nbinIMnpim_low, nbinIMnpim_high);
                                   //nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_woK0_wSid_n_woSp->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_woK0_wSid_n_woSp->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  
  for(unsigned int iwbin=0;iwbin<nwbin;iwbin++){
    IMnpim_IMnpip_dE_woK0_wSid_n_woSp_wbin[iwbin] = new TH2F(
        Form("IMnpim_IMnpip_woK0_wSid_n_woSp_wbin%d",iwbin),
        Form("IMnpim_IMnpip_woK0_wSid_n_woSp %0.2f-%0.2f",wbinlow[iwbin],wbinhigh[iwbin]),
        nbinIMnpip_wbin, IMnpip_wbinlow , IMnpip_wbinhigh, nbinIMnpim_wbin, IMnpim_wbinlow, IMnpim_wbinhigh);
    IMnpim_IMnpip_dE_woK0_wSid_n_woSp_wbin[iwbin]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
    IMnpim_IMnpip_dE_woK0_wSid_n_woSp_wbin[iwbin]->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  }
  /*
  for(unsigned int ibin=0;ibin<nbintemplate;ibin++){
    float binlow=1.0+(float)ibin*1./nbintemplate;
    float binhigh=1.0+(float)(ibin+1.0)*1./nbintemplate;
    IMnpim_IMnpip_dE_woK0_wSid_n_bin[ibin] = new TH2F(Form("IMnpim_IMnpip_woK0_wSid_n_bin%d",ibin),Form("IMnpim_IMnpip_woK0_wSid_n %0.2f-%0.2f",binlow,binhigh)
        ,nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
    IMnpim_IMnpip_dE_woK0_wSid_n_bin[ibin]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
    IMnpim_IMnpip_dE_woK0_wSid_n_bin[ibin]->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  }*/


  IMnpim_IMnpip_dE_woSid = new TH2F("IMnpim_IMnpip_dE_woSid","IMnpim_IMnpip_dE_woSid",nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_woSid->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_woSid->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  
  //IMnpim_IMnpip_dE_woSid_won = new TH2F("IMnpim_IMnpip_dE_woSid_won","IMnpim_IMnpip_dE_woSid_won",nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  //IMnpim_IMnpip_dE_woSid_won->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  //IMnpim_IMnpip_dE_woSid_won->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");

  IMnpim_IMnpip_dE_woK0_woSid = new TH2F("IMnpim_IMnpip_dE_woK0_woSid","IMnpim_IMnpip_dE_woK0_woSid",nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_woK0_woSid->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_woK0_woSid->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");

  IMnpim_IMnpip_dE_woK0_woSid_won = new TH2F("IMnpim_IMnpip_dE_woK0_woSid_won","IMnpim_IMnpip_dE_woK0_woSid_won",nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_woK0_woSid_won->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_woK0_woSid_won->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");

  IMnpim_IMnpip_dE_wK0_n = new TH2F("IMnpim_IMnpip_dE_wK0_n","IMnpim_IMnpip_dE_wK0_n",
                                   nbinIMnpip_bin, nbinIMnpip_low, nbinIMnpip_high, nbinIMnpim_bin, nbinIMnpim_low, nbinIMnpim_high);
  // nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_wK0_n->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_wK0_n->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");

  IMnpim_IMnpip_dE_wK0_wSid_n = new TH2F("IMnpim_IMnpip_dE_wK0_wSid_n","IMnpim_IMnpip_dE_wK0_wSid_n",nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_wK0_wSid_n->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_wK0_wSid_n->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  
  IMnpim_IMnpip_dE_wK0_wSid_n_Sp = new TH2F("IMnpim_IMnpip_dE_wK0_wSid_n_Sp","IMnpim_IMnpip_dE_wK0_wSid_n_Sp",nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_wK0_wSid_n_Sp->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_wK0_wSid_n_Sp->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  
  for(unsigned int iwbin=0;iwbin<nwbin;iwbin++){
    IMnpim_IMnpip_dE_wK0_wSid_n_Sp_wbin[iwbin] = new TH2F(Form("IMnpim_IMnpip_dE_wK0_wSid_n_Sp_wbin%d",iwbin),
                                                     Form("IMnpim_IMnpip_dE_wK0_wSid_n_Sp %0.2f-%0.2f",wbinlow[iwbin],wbinhigh[iwbin]),
                                                        nbinIMnpim_wbin,IMnpim_wbinlow,IMnpim_wbinhigh,nbinIMnpip_wbin,IMnpip_wbinlow,IMnpip_wbinhigh);
    IMnpim_IMnpip_dE_wK0_wSid_n_Sp_wbin[iwbin]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
    IMnpim_IMnpip_dE_wK0_wSid_n_Sp_wbin[iwbin]->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  }


  IMnpim_IMnpip_dE_wK0_wSid_n_Sm = new TH2F("IMnpim_IMnpip_dE_wK0_wSid_n_Sm","IMnpim_IMnpip_dE_wK0_wSid_n_Sm",nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_wK0_wSid_n_Sm->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_wK0_wSid_n_Sm->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");

  for(unsigned int iwbin=0;iwbin<nwbin;iwbin++){
    IMnpim_IMnpip_dE_wK0_wSid_n_Sm_wbin[iwbin] = new TH2F(Form("IMnpim_IMnpip_dE_wK0_wSid_n_Sm_wbin%d",iwbin),
                                                     Form("IMnpim_IMnpip_dE_wK0_wSid_n_Sm %0.2f-%0.2f",wbinlow[iwbin],wbinhigh[iwbin]),
                                                        nbinIMnpim_wbin,IMnpim_wbinlow,IMnpim_wbinhigh,nbinIMnpip_wbin,IMnpip_wbinlow,IMnpip_wbinhigh);
    IMnpim_IMnpip_dE_wK0_wSid_n_Sm_wbin[iwbin]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
    IMnpim_IMnpip_dE_wK0_wSid_n_Sm_wbin[iwbin]->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  }
  

  IMnpim_IMnpip_dE_wK0_wSid_n_SpSm = new TH2F(Form("IMnpim_IMnpip_dE_wK0_wSid_n_SpSm"),Form("IMnpim_IMnpip_dE_wK0_wSid_n_SpSm"),nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_wK0_wSid_n_SpSm->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_wK0_wSid_n_SpSm->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  
  IMnpim_IMnpip_dE_wK0_woSid_won = new TH2F(Form("IMnpim_IMnpip_dE_wK0_woSid_won"),Form("IMnpim_IMnpip_dE_wK0_woSid_won"),nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_wK0_woSid_won->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_wK0_woSid_won->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");

  IMnpim_IMnpip_dE_wK0_woSidn_won = new TH2F(Form("IMnpim_IMnpip_dE_wK0_woSidn_won"),Form("IMnpim_IMnpip_dE_wK0_woSidn_won"),nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_wK0_woSidn_won->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_wK0_woSidn_won->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");

  //IMnpip_CDHphi_dE_woK0_n = new TH2F("IMnpip_CDHphi_dE_woK0_n","IMnpip_CDHphi_dE_woK0_n",100,-3.5,3.5,nbinIMnpi,1.,2.);
  //IMnpip_CDHphi_dE_woK0_n->SetXTitle("CDHphi [radian]");
  //IMnpip_CDHphi_dE_woK0_n->SetYTitle("IM(n#pi^{+}) [GeV/c^{2}]");

  //IMnpip_CDHz_dE_woK0_n = new TH2F("IMnpip_CDHz_dE_woK0_n","IMnpip_CDHz_dE_woK0_n",100,-50,50,nbinIMnpi,1.,2.);
  //IMnpip_CDHz_dE_woK0_n->SetXTitle("CDHz [cm]");
  //IMnpip_CDHz_dE_woK0_n->SetYTitle("IM(n#pi^{+}) [GeV/c^{2}]");

  //IMnpim_CDHphi_dE_woK0_n = new TH2F("IMnpim_CDHphi_dE_woK0_n","IMnpim_CDHphi_dE_woK0_n",100,-3.5,3.5,nbinIMnpi,1.,2.);
  //IMnpim_CDHphi_dE_woK0_n->SetXTitle("CDHphi [radian]");
  //IMnpim_CDHphi_dE_woK0_n->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");

  //IMnpim_CDHz_dE_woK0_n = new TH2F("IMnpim_CDHz_dE_woK0_n","IMnpim_CDHz_dE_woK0_n",100,-50,50,nbinIMnpi,1.,2.);
  //IMnpim_CDHz_dE_woK0_n->SetXTitle("CDHz [cm]");
  //IMnpim_CDHz_dE_woK0_n->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");

  //pipmom_IMpippim= new TH2F("pipmom_IMpippim","pipmom_IMpippim",nbinpippim,0.,0.9,200,0,1);
  //pipmom_IMpippim->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}] ");
  //pipmom_IMpippim->SetYTitle("Mom(#pi^{+}) [GeV/c] ");

  //pipmom_IMpippim_dE= new TH2F("pipmom_IMpippim_dE","pipmom_IMpippim_dE",nbinpippim,0.,0.9,200,0,1);
  //pipmom_IMpippim_dE->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}] ");
  //pipmom_IMpippim_dE->SetYTitle("Mom(#pi^{+}) [GeV/c] ");

  //pipmom_IMpippim_dE_n = new TH2F("pipmom_IMpippim_dE_n","pipmom_IMpippim_dE_n",nbinpippim,0.,0.9,200,0,1);
  //pipmom_IMpippim_dE_n->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}] ");
  //pipmom_IMpippim_dE_n->SetYTitle("Mom(#pi^{+}) [GeV/c] ");

  //pimmom_IMpippim= new TH2F("pimmom_IMpippim","pimmom_IMpippim",nbinpippim,0.,0.9,200,0,1);
  //pimmom_IMpippim->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}] ");
  //pimmom_IMpippim->SetYTitle("Mom(#pi^{-}) [GeV/c] ");

  //pimmom_IMpippim_dE= new TH2F("pimmom_IMpippim_dE","pimmom_IMpippim_dE",nbinpippim,0.,0.9,200,0,1);
  //pimmom_IMpippim_dE->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}] ");
  //pimmom_IMpippim_dE->SetYTitle("Mom(#pi^{-}) [GeV/c] ");

  //pimmom_IMpippim_dE_n = new TH2F("pimmom_IMpippim_dE_n","pimmom_IMpippim_dE_n",nbinpippim,0.,0.9,200,0,1);
  //pimmom_IMpippim_dE_n->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}] ");
  //pimmom_IMpippim_dE_n->SetYTitle("Mom(#pi^{-}) [GeV/c] ");


  //pipmom_IMnpip= new TH2F("pipmom_IMnpip","pipmom_IMnpip",nbinIMnpi,1.,2.,200,0,1);
  //pipmom_IMnpip->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}] ");
  //pipmom_IMnpip->SetYTitle("Mom(#pi^{+}) [GeV/c] ");

  //pipmom_IMnpip_dE= new TH2F("pipmom_IMnpip_dE","pipmom_IMnpip_dE",nbinIMnpi,1.,2.,200,0,1);
  //pipmom_IMnpip_dE->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}] ");
  //pipmom_IMnpip_dE->SetYTitle("Mom(#pi^{+}) [GeV/c] ");

  //pipmom_IMnpip_dE_n = new TH2F("pipmom_IMnpip_dE_n","pipmom_IMnpip_dE_n",nbinIMnpi,1.,2.,200,0,1);
  //pipmom_IMnpip_dE_n->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}] ");
  //pipmom_IMnpip_dE_n->SetYTitle("Mom(#pi^{+}) [GeV/c] ");

  //pipmom_MMnmiss_dE = new TH2F("pipmom_MMnmiss_dE","pipmom_MMnmiss_dE",nbinnmiss, nmisslow, nmisshigh,200,0,1);
  //pipmom_MMnmiss_dE->SetXTitle("Miss Mass. [GeV/c^{2}]");
  //pipmom_MMnmiss_dE->SetYTitle("Mom(#pi^{+}) [GeV/c] ");

  //pipmom_MMnmiss_dE_woK0_woSidn = new TH2F("pipmom_MMnmiss_dE_woK0_woSidn","pipmom_MMnmiss_dE_woK0_woSidn",nbinnmiss, nmisslow, nmisshigh,200,0,1);
  //pipmom_MMnmiss_dE_woK0_woSidn->SetXTitle("Miss Mass. [GeV/c^{2}]");
  //pipmom_MMnmiss_dE_woK0_woSidn->SetYTitle("Mom(#pi^{+}) [GeV/c] ");

  //pipmom_MMnmiss_dE_woK0_woSid_won = new TH2F("pipmom_MMnmiss_dE_woK0_woSid_won","pipmom_MMnmiss_dE_woK0_woSid_won",nbinnmiss, nmisslow, nmisshigh,200,0,1);
  //pipmom_MMnmiss_dE_woK0_woSid_won->SetXTitle("Miss Mass. [GeV/c^{2}]");
  //pipmom_MMnmiss_dE_woK0_woSid_won->SetYTitle("Mom(#pi^{+}) [GeV/c] ");

  //pipmom_MMnmiss_dE_wK0_woSid_won = new TH2F("pipmom_MMnmiss_dE_wK0_woSid_won","pipmom_MMnmiss_dE_wK0_woSid_won",nbinnmiss, nmisslow, nmisshigh,200,0,1);
  //pipmom_MMnmiss_dE_wK0_woSid_won->SetXTitle("Miss Mass. [GeV/c^{2}]");
  //pipmom_MMnmiss_dE_wK0_woSid_won->SetYTitle("Mom(#pi^{+}) [GeV/c] ");

  //pipmom_MMnmiss_dE_wK0_woSidn_won = new TH2F("pipmom_MMnmiss_dE_wK0_woSidn_won","pipmom_MMnmiss_dE_wK0_woSidn_won",nbinnmiss, nmisslow, nmisshigh,200,0,1);
  //pipmom_MMnmiss_dE_wK0_woSidn_won->SetXTitle("Miss Mass. [GeV/c^{2}]");
  //pipmom_MMnmiss_dE_wK0_woSidn_won->SetYTitle("Mom(#pi^{+}) [GeV/c] ");

  //pipmom_pimmom_dE_woK0_woSidn = new TH2F("pipmom_pimmom_dE_woK0_woSidn","pipmom_pimmom_dE_woK0_woSidn",100,0,1.0,100,0,1.0);
  //pipmom_pimmom_dE_woK0_woSidn->SetXTitle("Mom(#pi^{-}) [GeV/c] ");
  //pipmom_pimmom_dE_woK0_woSidn->SetYTitle("Mom(#pi^{+}) [GeV/c] ");

  //pipmom_pimmom_dE_woK0_woSid_won = new TH2F("pipmom_pimmom_dE_woK0_woSid_won","pipmom_pimmom_dE_woK0_woSid_won",100,0,1.0,100,0,1.0);
  //pipmom_pimmom_dE_woK0_woSid_won->SetXTitle("Mom(#pi^{-}) [GeV/c] ");
  //pipmom_pimmom_dE_woK0_woSid_won->SetYTitle("Mom(#pi^{+}) [GeV/c] ");

  //pipmom_pimmom_dE_wK0_woSid_won = new TH2F("pipmom_pimmom_dE_wK0_woSid_won","pipmom_pimmom_dE_wK0_woSid_won",100,0,1.0,100,0,1.0);
  //pipmom_pimmom_dE_wK0_woSid_won->SetXTitle("Mom(#pi^{-}) [GeV/c] ");
  //pipmom_pimmom_dE_wK0_woSid_won->SetYTitle("Mom(#pi^{+}) [GeV/c] ");



  //pimmom_IMnpim = new TH2F("pimmom_IMnpim","pimmom_IMnpim",nbinIMnpi,1.,2.,200,0,1);
  //pimmom_IMnpim->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}] ");
  //pimmom_IMnpim->SetYTitle("Mom(#pi^{-}) [GeV/c] ");

  //pimmom_IMnpim_dE = new TH2F("pimmom_IMnpim_dE","pimmom_IMnpim_dE",nbinIMnpi,1.,2.,200,0,1);
  //pimmom_IMnpim_dE->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}] ");
  //pimmom_IMnpim_dE->SetYTitle("Mom(#pi^{-}) [GeV/c] ");

  //pimmom_IMnpim_dE_n = new TH2F("pimmom_IMnpim_dE_n","pimmom_IMnpim_dE_n",nbinIMnpi,1.,2.,200,0,1);
  //pimmom_IMnpim_dE_n->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}] ");
  //pimmom_IMnpim_dE_n->SetYTitle("Mom(#pi^{-}) [GeV/c] ");

  //pimmom_MMnmiss_dE = new TH2F("pimmom_MMnmiss_dE","pimmom_MMnmiss_dE",nbinnmiss, nmisslow, nmisshigh,200,0,1);
  //pimmom_MMnmiss_dE->SetXTitle("Miss Mass. [GeV/c^{2}]");
  //pimmom_MMnmiss_dE->SetYTitle("Mom(#pi^{-}) [GeV/c] ");

  //pimmom_MMnmiss_dE_woK0_woSidn = new TH2F("pimmom_MMnmiss_dE_woK0_woSidn","pimmom_MMnmiss_dE_woK0_woSidn",nbinnmiss, nmisslow, nmisshigh,200,0,1);
  //pimmom_MMnmiss_dE_woK0_woSidn->SetXTitle("Miss Mass. [GeV/c^{2}]");
  //pimmom_MMnmiss_dE_woK0_woSidn->SetYTitle("Mom(#pi^{-}) [GeV/c] ");

  //pimmom_MMnmiss_dE_woK0_woSid_won = new TH2F("pimmom_MMnmiss_dE_woK0_woSid_won","pimmom_MMnmiss_dE_woK0_woSid_won",nbinnmiss, nmisslow, nmisshigh,200,0,1);
  //pimmom_MMnmiss_dE_woK0_woSid_won->SetXTitle("Miss Mass. [GeV/c^{2}]");
  //pimmom_MMnmiss_dE_woK0_woSid_won->SetYTitle("Mom(#pi^{-}) [GeV/c] ");

  //pimmom_MMnmiss_dE_wK0_woSid_won = new TH2F("pimmom_MMnmiss_dE_wK0_woSid_won","pimmom_MMnmiss_dE_wK0_woSid_won",nbinnmiss, nmisslow, nmisshigh,200,0,1);
  //pimmom_MMnmiss_dE_wK0_woSid_won->SetXTitle("Miss Mass. [GeV/c^{2}]");
  //pimmom_MMnmiss_dE_wK0_woSid_won->SetYTitle("Mom(#pi^{-}) [GeV/c] ");

  //pimmom_MMnmiss_dE_wK0_woSidn_won = new TH2F("pimmom_MMnmiss_dE_wK0_woSidn_won","pimmom_MMnmiss_dE_wK0_woSidn_won",nbinnmiss, nmisslow, nmisshigh,200,0,1);
  //pimmom_MMnmiss_dE_wK0_woSidn_won->SetXTitle("Miss Mass. [GeV/c^{2}]");
  //pimmom_MMnmiss_dE_wK0_woSidn_won->SetYTitle("Mom(#pi^{-}) [GeV/c] ");
  
  pipmom_Momnpip_wSid_n = new TH2F("pipmom_Momnpip_wSid_n","pipmom_Momnpip_wSid_n",200,0,1,200,0.,1.0);
  pipmom_Momnpip_wSid_n->SetXTitle("Mom(n#pi^{+}) [GeV/c]");
  pipmom_Momnpip_wSid_n->SetYTitle("Mom(#pi^{+}) [GeV/c]");
  
  pimmom_Momnpim_wSid_n = new TH2F("pimmom_Momnpim_wSid_n","pimmom_Momnpim_wSid_n",200,0,1,200,0.,1.0);
  pimmom_Momnpim_wSid_n->SetXTitle("Mom(n#pi^{-}) [GeV/c]");
  pimmom_Momnpim_wSid_n->SetYTitle("Mom(#pi^{-}) [GeV/c]");

  IMpippim_DCApipibeam = new TH2F("IMpippim_DCApipibeam","IMpippim_DCApipibeam",300,0,20,nbinpippim,0.,0.9);
  IMpippim_DCApipibeam->SetXTitle("DCA #pi#pi-beam [cm]");
  IMpippim_DCApipibeam->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}] ");

  IMpippim_DCApipibeam_n = new TH2F("IMpippim_DCApipibeam_n","IMpippim_DCApipibeam_n",300,0,20,nbinpippim,0.,0.9);
  IMpippim_DCApipibeam_n->SetXTitle("DCA #pi#pi-beam [cm]");
  IMpippim_DCApipibeam_n->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}] ");

  IMpippim_DCApipibeam_wK0_n = new TH2F("IMpippim_DCApipibeam_wK0_n","IMpippim_DCApipibeam_wK0_n",300,0,20,nbinpippim,0.,0.9);
  IMpippim_DCApipibeam_wK0_n->SetXTitle("DCA #pi#pi-beam [cm]");
  IMpippim_DCApipibeam_wK0_n->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}] ");

  IMpippim_DCApipibeam_wK0_woSid_n = new TH2F("IMpippim_DCApipibeam_wK0_woSid_n","IMpippim_DCApipibeam_wK0_woSid_n",300,0,20,nbinpippim,0.,0.9);
  IMpippim_DCApipibeam_wK0_woSid_n->SetXTitle("DCA #pi#pi-beam [cm]");
  IMpippim_DCApipibeam_wK0_woSid_n->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}] ");

  IMnpip_DCApipibeam = new TH2F("IMnpip_DCApipibeam","IMnpip_DCApipibeam",300,0,20,nbinIMnpi,1.,2.);
  IMnpip_DCApipibeam->SetXTitle("DCA #pi#pi-beam [cm]");
  IMnpip_DCApipibeam->SetYTitle("IM(n#pi^{+}) [GeV/c^{2}] ");

  IMnpip_DCApipibeam_n = new TH2F("IMnpip_DCApipibeam_n","IMnpip_DCApipibeam_n",300,0,20,nbinIMnpi,1.,2.);
  IMnpip_DCApipibeam_n->SetXTitle("DCA #pi#pi-beam [cm]");
  IMnpip_DCApipibeam_n->SetYTitle("IM(n#pi^{+}) [GeV/c^{2}] ");

  IMnpip_DCApipibeam_woK0_n = new TH2F("IMnpip_DCApipibeam_woK0_n","IMnpip_DCApipibeam_woK0_n",300,0,20,nbinIMnpi,1.,2.);
  IMnpip_DCApipibeam_woK0_n->SetXTitle("DCA #pi#pi-beam [cm]");
  IMnpip_DCApipibeam_woK0_n->SetYTitle("IM(n#pi^{+}) [GeV/c^{2}] ");

  IMnpip_DCApipibeam_woK0_n_Sp = new TH2F("IMnpip_DCApipibeam_woK0_n_Sp","IMnpip_DCApipibeam_woK0_n_Sp",300,0,20,nbinIMnpi,1.,2.);
  IMnpip_DCApipibeam_woK0_n_Sp->SetXTitle("DCA #pi#pi-beam [cm]");
  IMnpip_DCApipibeam_woK0_n_Sp->SetYTitle("IM(n#pi^{+}) [GeV/c^{2}] ");

  IMnpim_DCApipibeam = new TH2F("IMnpim_DCApipibeam","IMnpim_DCApipibeam",300,0,20,nbinIMnpi,1.,2.);
  IMnpim_DCApipibeam->SetXTitle("DCA #pi#pi-beam [cm]");
  IMnpim_DCApipibeam->SetYTitle("IM(n#pi^{+}) [GeV/c^{2}] ");

  IMnpim_DCApipibeam_n = new TH2F("IMnpim_DCApipibeam_n","IMnpim_DCApipibeam_n",300,0,20,nbinIMnpi,1.,2.);
  IMnpim_DCApipibeam_n->SetXTitle("DCA #pi#pi-beam [cm]");
  IMnpim_DCApipibeam_n->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}] ");

  IMnpim_DCApipibeam_woK0_n = new TH2F("IMnpim_DCApipibeam_woK0_n","IMnpim_DCApipibeam_woK0_n",300,0,20,nbinIMnpi,1.,2.);
  IMnpim_DCApipibeam_woK0_n->SetXTitle("DCA #pi#pi-beam [cm]");
  IMnpim_DCApipibeam_woK0_n->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}] ");

  IMnpim_DCApipibeam_woK0_n_Sm = new TH2F("IMnpim_DCApipibeam_woK0_n_Sm","IMnpim_DCApipibeam_woK0_n_Sm",300,0,20,nbinIMnpi,1.,2.);
  IMnpim_DCApipibeam_woK0_n_Sm->SetXTitle("DCA #pi#pi-beam [cm]");
  IMnpim_DCApipibeam_woK0_n_Sm->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}] ");

  MMnmiss_DCApipibeam_wK0 = new TH2F("MMnmiss_DCApipibeam_wK0","MMnmiss_DCApipibeam_wK0",300,0,20,100,0.4,1.9);
  MMnmiss_DCApipibeam_wK0->SetXTitle("DCA #pi#pi-beam [cm]");
  MMnmiss_DCApipibeam_wK0->SetYTitle("Miss. Mass [GeV/c^{2}]");

  MMnmiss_DCApipibeam_woK0_wSid = new TH2F("MMnmiss_DCApipibeam_woK0_wSid","MMnmiss_DCApipibeam_woK0_wSid",300,0,20,100,0.4,1.9);
  MMnmiss_DCApipibeam_woK0_wSid->SetXTitle("DCA #pi#pi-beam [cm]");
  MMnmiss_DCApipibeam_woK0_wSid->SetYTitle("Miss. Mass [GeV/c^{2}]");

  IMnpip_DCApip_dE_woK0_n = new TH2F("IMnpip_DCApip_dE_woK0_n","IMnpip_DCApip_dE_woK0_n",300,0,20,nbinIMnpi,1.,2.);
  IMnpip_DCApip_dE_woK0_n->SetXTitle("DCA #pi^{+} [cm]");
  IMnpip_DCApip_dE_woK0_n->SetYTitle("IM(n#pi^{+}) [GeV/c^{2}]");

  IMnpip_DCApim_dE_woK0_n = new TH2F("IMnpip_DCApim_dE_woK0_n","IMnpip_DCApim_dE_woK0_n",300,0,20,nbinIMnpi,1.,2.);
  IMnpip_DCApim_dE_woK0_n->SetXTitle("DCA #pi^{-} [cm]");
  IMnpip_DCApim_dE_woK0_n->SetYTitle("IM(n#pi^{+}) [GeV/c^{2}]");

  IMnpim_DCApip_dE_woK0_n = new TH2F("IMnpim_DCApip_dE_woK0_n","IMnpim_DCApip_dE_woK0_n",300,0,20,nbinIMnpi,1.,2.);
  IMnpim_DCApip_dE_woK0_n->SetXTitle("DCA #pi^{+} [cm]");
  IMnpim_DCApip_dE_woK0_n->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");

  IMnpim_DCApim_dE_woK0_n = new TH2F("IMnpim_DCApim_dE_woK0_n","IMnpim_DCApim_dE_woK0_n",300,0,20,nbinIMnpi,1.,2.);
  IMnpim_DCApim_dE_woK0_n->SetXTitle("DCA #pi^{-} [cm]");
  IMnpim_DCApim_dE_woK0_n->SetYTitle("IM(n#pi^{-} [GeV/c^{2}]");
  
  /*
  for(int igap=0; igap<ngap; igap++) {
    IMnpim_IMnpip_dE_n_side[igap] = new TH2F(Form("IMnpim_IMnpip_dE_n_side_%d",igap),Form("IMnpim_IMnpip_dE_n_side_%d",igap),nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
    IMnpim_IMnpip_dE_n_side[igap]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
    IMnpim_IMnpip_dE_n_side[igap]->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");

    IMnpim_IMnpip_dE_woK0_n_side[igap] = new TH2F(Form("IMnpim_IMnpip_dE_woK0_n_side_%d",igap),Form("IMnpim_IMnpip_dE_woK0_n_side_%d",igap),nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
    IMnpim_IMnpip_dE_woK0_n_side[igap]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
    IMnpim_IMnpip_dE_woK0_n_side[igap]->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  }
  */
  std::cout << __LINE__ << std::endl;
  for(int i=0; i<2; i++) {
    for(int igap=0; igap<ngap; igap++) {
      const char  lh[][6]= {"low","high"};
      //IMnpim_IMnpip_dE_n_Sp_side[i][igap] = new TH2F(Form("IMnpim_IMnpip_dE_n_Sp_side_%s_%d",lh[i],igap),Form("IMnpim_IMnpip_dE_n_Sp_side_%s_%d",lh[i],igap), nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
      //IMnpim_IMnpip_dE_n_Sp_side[i][igap]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
      //IMnpim_IMnpip_dE_n_Sp_side[i][igap]->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");

      //IMnpim_IMnpip_dE_n_Sm_side[i][igap] = new TH2F(Form("IMnpim_IMnpip_dE_n_Sm_side_%s_%d",lh[i],igap),Form("IMnpim_IMnpip_dE_n_Sm_side_%s_%d",lh[i],igap), nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
      //IMnpim_IMnpip_dE_n_Sm_side[i][igap]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
      //IMnpim_IMnpip_dE_n_Sm_side[i][igap]->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");

      //IMnpim_IMnpip_dE_woK0_n_Sp_side[i][igap] = new TH2F(Form("IMnpim_IMnpip_dE_woK0_n_Sp_side_%s_%d",lh[i],igap),Form("IMnpim_IMnpip_dE_woK0_n_Sp_side_%s_%d",lh[i],igap), nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
      //IMnpim_IMnpip_dE_woK0_n_Sp_side[i][igap]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
      //IMnpim_IMnpip_dE_woK0_n_Sp_side[i][igap]->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");

      //IMnpim_IMnpip_dE_woK0_n_Sm_side[i][igap] = new TH2F(Form("IMnpim_IMnpip_dE_woK0_n_Sm_side_%s_%d",lh[i],igap),Form("IMnpim_IMnpip_dE_woK0_n_Sm_side_%s_%d",lh[i],igap), nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
      //IMnpim_IMnpip_dE_woK0_n_Sm_side[i][igap]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
      //IMnpim_IMnpip_dE_woK0_n_Sm_side[i][igap]->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
    }
  }
  
  /*
  for(int izone=0; izone<nzone; izone++) {
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
  }*/

  IMnpim_IMnpip_dE_n_Sp = new TH2F(Form("IMnpim_IMnpip_dE_n_Sp"),Form("IMnpim_IMnpip_dE_n_Sp"), nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_n_Sp->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_n_Sp->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");

  IMnpim_IMnpip_dE_woK0_n_Sp = new TH2F(Form("IMnpim_IMnpip_dE_woK0_n_Sp"),Form("IMnpim_IMnpip_dE_woK0_n_Sp"), nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_woK0_n_Sp->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_woK0_n_Sp->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");

  //IMnpim_IMnpip_dE_woK0_n_Sp_bg = new TH2F(Form("IMnpim_IMnpip_dE_woK0_n_Sp_bg"),Form("IMnpim_IMnpip_dE_woK0_n_Sp_bg"), nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  //IMnpim_IMnpip_dE_woK0_n_Sp_bg->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  //IMnpim_IMnpip_dE_woK0_n_Sp_bg->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");

  IMnpim_IMnpip_dE_n_Sm = new TH2F(Form("IMnpim_IMnpip_dE_n_Sm"),Form("IMnpim_IMnpip_dE_n_Sm"), nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_n_Sm->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_n_Sm->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");

  IMnpim_IMnpip_dE_woK0_n_Sm = new TH2F(Form("IMnpim_IMnpip_dE_woK0_n_Sm"),Form("IMnpim_IMnpip_dE_woK0_n_Sm"), nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  IMnpim_IMnpip_dE_woK0_n_Sm->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_woK0_n_Sm->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");

  //IMnpim_IMnpip_dE_woK0_n_Sm_bg = new TH2F(Form("IMnpim_IMnpip_dE_woK0_n_Sm_bg"),Form("IMnpim_IMnpip_dE_woK0_n_Sm_bg"), nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  //IMnpim_IMnpip_dE_woK0_n_Sm_bg->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  //IMnpim_IMnpip_dE_woK0_n_Sm_bg->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");

  //IMnpim_IMnpip_dE_woK0_woSidn = new TH2F(Form("IMnpim_IMnpip_dE_woK0_woSidn"),Form("IMnpim_IMnpip_dE_woK0_woSidn"), nbinIMnpi, 1, 2.0, nbinIMnpi, 1, 2.0);
  //IMnpim_IMnpip_dE_woK0_woSidn->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  //IMnpim_IMnpip_dE_woK0_woSidn->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");

  nmom_IMnpim_dE_n = new TH2F("nmom_IMnpim_dE_n","nmom_IMnpim_dE_n",nbinIMnpi,1.,2.,nbinnmom,0.,1.0);
  nmom_IMnpim_dE_n->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  nmom_IMnpim_dE_n->SetYTitle("nmom. [GeV/c]");

  nmom_IMnpim_dE_n_Sm = new TH2F("nmom_IMnpim_dE_n_Sm","nmom_IMnpim_dE_n_Sm",nbinIMnpi,1.,2.,nbinnmom,0.,1.0);
  nmom_IMnpim_dE_n_Sm->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  nmom_IMnpim_dE_n_Sm->SetYTitle("nmom. [GeV/c]");

  nmom_Momnpim_woK0_n_Sm = new TH2F("nmom_Momnpim_woK0_n_Sm","nmom_Momnpim_woK0_n_Sm",100,0.,1.,nbinnmom,0.,1.0);
  nmom_Momnpim_woK0_n_Sm->SetXTitle("Mom.(n#pi^{-}) [GeV/c]");
  nmom_Momnpim_woK0_n_Sm->SetYTitle("nmom. [GeV/c]");

  nmom_IMnpim_dE_woK0_n = new TH2F("nmom_IMnpim_dE_woK0_n","nmom_IMnpim_dE_woK0_n",nbinIMnpi,1.,2.,nbinnmom,0.,1.0);
  nmom_IMnpim_dE_woK0_n->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  nmom_IMnpim_dE_woK0_n->SetYTitle("nmom. [GeV/c]");

  nmom_IMnpip_dE_n = new TH2F("nmom_IMnpip_dE_n","nmom_IMnpip_dE_n",nbinIMnpi,1.,2.,nbinnmom,0.,1.0);
  nmom_IMnpip_dE_n->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  nmom_IMnpip_dE_n->SetYTitle("nmom. [GeV/c]");

  nmom_IMnpip_dE_n_Sp = new TH2F("nmom_IMnpip_dE_n_Sp","nmom_IMnpip_dE_n_Sp",nbinIMnpi,1.,2.,nbinnmom,0.,1.0);
  nmom_IMnpip_dE_n_Sp->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  nmom_IMnpip_dE_n_Sp->SetYTitle("nmom. [GeV/c]");

  nmom_IMnpip_dE_woK0_n = new TH2F("nmom_IMnpip_dE_woK0_n","nmom_IMnpip_dE_woK0_n",nbinIMnpi,1.,2.,nbinnmom,0.,1.0);
  nmom_IMnpip_dE_woK0_n->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  nmom_IMnpip_dE_woK0_n->SetYTitle("nmom. [GeV/c]");
  
  nmom_IMnpip_dE_woK0_woSid_won = new TH2F("nmom_IMnpip_dE_woK0_woSid_won","nmom_IMnpip_dE_woK0_woSid_won",nbinIMnpi,1.,2.,nbinnmom,0.,1.0);
  nmom_IMnpip_dE_woK0_woSid_won->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  nmom_IMnpip_dE_woK0_woSid_won->SetYTitle("nmom. [GeV/c]");
  
  nmom_IMnpip_dE_wK0_woSid_won = new TH2F("nmom_IMnpip_dE_wK0_woSid_won","nmom_IMnpip_dE_wK0_woSid_won",nbinIMnpi,1.,2.,nbinnmom,0.,1.0);
  nmom_IMnpip_dE_wK0_woSid_won->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  nmom_IMnpip_dE_wK0_woSid_won->SetYTitle("nmom. [GeV/c]");

  nmom_Momnpip_woK0_n_Sp = new TH2F("nmom_Momnpip_woK0_n_Sp","nmom_Momnpip_woK0_n_Sp",100,0.,1.,nbinnmom,0.,1.0);
  nmom_Momnpip_woK0_n_Sp->SetXTitle("Mom.(n#pi^{+}) [GeV/c]");
  nmom_Momnpip_woK0_n_Sp->SetYTitle("nmom. [GeV/c]");

  MMnpim_MMnpip = new TH2F("MMnpim_MMnpip","MMnpim_MMnpip",280, 1, 1.7, 280, 1, 1.7);
  MMnpim_MMnpip->SetXTitle("Miss. Mass (K^{-}d #rightarrow n#pi^{+}) [GeV/c^{2}]");
  MMnpim_MMnpip->SetYTitle("Miss. Mass (K^{-}d #rightarrow n#pi^{-}) [GeV/c^{2}]");
  
  MMnpim_MMnpip_n = new TH2F("MMnpim_MMnpip_n","MMnpim_MMnpip_n",280, 1, 1.7, 280, 1, 1.7);
  MMnpim_MMnpip_n->SetXTitle("Miss. Mass (K^{-}d #rightarrow n#pi^{+}) [GeV/c^{2}]");
  MMnpim_MMnpip_n->SetYTitle("Miss. Mass (K^{-}d #rightarrow n#pi^{-}) [GeV/c^{2}]");

  MMnpim_MMnpip_woSid_n = new TH2F("MMnpim_MMnpip_woSid_n","MMnpim_MMnpip_woSid_n",280, 1, 1.7, 280, 1, 1.7);
  MMnpim_MMnpip_woSid_n->SetXTitle("Miss. Mass (K^{-}d #rightarrow n#pi^{+}) [GeV/c^{2}]");
  MMnpim_MMnpip_woSid_n->SetYTitle("Miss. Mass (K^{-}d #rightarrow n#pi^{-}) [GeV/c^{2}]");

  MMnpim_MMnpip_woK0_woSid_n = new TH2F("MMnpim_MMnpip_woK0_woSid_n","MMnpim_MMnpip_woK0_woSid_n",280, 1, 1.7, 280, 1, 1.7);
  MMnpim_MMnpip_woK0_woSid_n->SetXTitle("Miss. Mass (K^{-}d #rightarrow n#pi^{+}) [GeV/c^{2}]");
  MMnpim_MMnpip_woK0_woSid_n->SetYTitle("Miss. Mass (K^{-}d #rightarrow n#pi^{-}) [GeV/c^{2}]");

  MMnpim_MMnpip_mc = new TH2F("MMnpim_MMnpip_mc","MMnpim_MMnpip_mc",280, 1, 1.7, 280, 1, 1.7);
  MMnpim_MMnpip_mc->SetXTitle("true IM(n_{miss}#pi^{+}) [GeV/c^{2}]");
  MMnpim_MMnpip_mc->SetYTitle("true IM(n_{miss}#pi^{-}) [GeV/c^{2}]");

  MMnpim_MMnpip_woK0_n = new TH2F("MMnpim_MMnpip_woK0_n","MMnpim_MMnpip_woK0_n",140, 1, 1.7, 140, 1, 1.7);
  MMnpim_MMnpip_woK0_n->SetXTitle("Miss. Mass (K^{-}d #rightarrow n#pi^{+}) [GeV/c^{2}]");
  MMnpim_MMnpip_woK0_n->SetYTitle("Miss. Mass (K^{-}d #rightarrow n#pi^{-}) [GeV/c^{2}]");

  MMnpim_MMnpip_wSid_n = new TH2F("MMnpim_MMnpip_wSid_n","MMnpim_MMnpip_wSid_n",140, 1, 1.7, 140, 1, 1.7);
  MMnpim_MMnpip_wSid_n->SetXTitle("Miss. Mass (K^{-}d #rightarrow n#pi^{+}) [GeV/c^{2}]");
  MMnpim_MMnpip_wSid_n->SetYTitle("Miss. Mass (K^{-}d #rightarrow n#pi^{-}) [GeV/c^{2}]");
  
  MMnpim_MMnpip_wSid_n_mc = new TH2F("MMnpim_MMnpip_wSid_n_mc","MMnpim_MMnpip_wSid_n_mc",280, 1, 1.7, 280, 1, 1.7);
  MMnpim_MMnpip_wSid_n_mc->SetXTitle("true IM(n_{miss}#pi^{+}) [GeV/c^{2}]");
  MMnpim_MMnpip_wSid_n_mc->SetYTitle("true IM(n_{miss}#pi^{-}) [GeV/c^{2}]");

  MMnpim_MMnpip_woK0_wSid_n = new TH2F("MMnpim_MMnpip_woK0_wSid_n","MMnpim_MMnpip_woK0_wSid_n",140, 1, 1.7, 140, 1, 1.7);
  MMnpim_MMnpip_woK0_wSid_n->SetXTitle("Miss. Mass (K^{-}d #rightarrow n#pi^{+}) [GeV/c^{2}]");
  MMnpim_MMnpip_woK0_wSid_n->SetYTitle("Miss. Mass (K^{-}d #rightarrow n#pi^{-}) [GeV/c^{2}]");

  dE_CDHphi = new TH2F(Form("dE_CDHphi"),Form("dE_CDHphi"),100,-3.14,3.14, nbindE,0,50);
  dE_CDHphi->SetXTitle("CDH phi");
  dE_CDHphi->SetYTitle("dE [MeVee]");

  dE_CDHz = new TH2F(Form("dE_CDHz"),"dE_CDHz",100,-50,50,nbindE,0,50);
  dE_CDHz->SetXTitle("CDH z [cm]");
  dE_CDHz->SetYTitle("dE [MeVee]");

  dE_IMnpim = new TH2F(Form("dE_IMnpim"),Form("dE_IMnpim"), nbinIMnpi, 1.0, 2.0, nbindE, 0, 50.);
  dE_IMnpim->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  dE_IMnpim->SetYTitle("dE [MeVee]");

  //dE_IMnpim_woK0 = new TH2F(Form("dE_IMnpim_woK0"),Form("dE_IMnpim_woK0"), nbinIMnpi, 1.0, 2.0, nbindE, 0, 50.);
  //dE_IMnpim_woK0->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  //dE_IMnpim_woK0->SetYTitle("dE [MeVee]");

  dE_IMnpim_n = new TH2F(Form("dE_IMnpim_n"),Form("dE_IMnpim_n"), nbinIMnpi, 1.0, 2.0, nbindE, 0, 50.);
  dE_IMnpim_n->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  dE_IMnpim_n->SetYTitle("dE [MeVee]");

  //dE_IMnpim_woK0_n = new TH2F(Form("dE_IMnpim_woK0_n"),Form("dE_IMnpim_woK0_n"), nbinIMnpi, 1.0, 2.0, nbindE, 0, 50.);
  //dE_IMnpim_woK0_n->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  //dE_IMnpim_woK0_n->SetYTitle("dE [MeVee]");
  
  nmom_IMnpim_dE_woK0_woSid_won = new TH2F("nmom_IMnpim_dE_woK0_woSid_won","nmom_IMnpim_dE_woK0_woSid_won",nbinIMnpi,1.,2.,nbinnmom,0.,1.0);
  nmom_IMnpim_dE_woK0_woSid_won->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  nmom_IMnpim_dE_woK0_woSid_won->SetYTitle("nmom. [GeV/c]");
  
  nmom_IMnpim_dE_wK0_woSid_won = new TH2F("nmom_IMnpim_dE_wK0_woSid_won","nmom_IMnpim_dE_wK0_woSid_won",nbinIMnpi,1.,2.,nbinnmom,0.,1.0);
  nmom_IMnpim_dE_wK0_woSid_won->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  nmom_IMnpim_dE_wK0_woSid_won->SetYTitle("nmom. [GeV/c]");

  dE_IMnpip = new TH2F(Form("dE_IMnpip"),Form("dE_IMnpip"), nbinIMnpi, 1.0, 2.0, nbindE, 0, 50.);
  dE_IMnpip->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  dE_IMnpip->SetYTitle("dE [MeVee]");

  //dE_IMnpip_woK0 = new TH2F(Form("dE_IMnpip_woK0"),Form("dE_IMnpip_woK0"), nbinIMnpi, 1.0, 2.0, nbindE, 0, 50.);
  //dE_IMnpip_woK0->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  //dE_IMnpip_woK0->SetYTitle("dE [MeVee]");

  dE_IMnpip_n = new TH2F(Form("dE_IMnpip_n"),Form("dE_IMnpip_n"), nbinIMnpi, 1.0, 2.0, nbindE, 0, 50.);
  dE_IMnpip_n->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  dE_IMnpip_n->SetYTitle("dE [MeVee]");

  //dE_IMnpip_woK0_n = new TH2F(Form("dE_IMnpip_woK0_n"),Form("dE_IMnpip_woK0_n"), nbinIMnpi, 1.0, 2.0, nbindE, 0, 50.);
  //dE_IMnpip_woK0_n->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  //dE_IMnpip_woK0_n->SetYTitle("dE [MeVee]");

  //dE_IMnpipi_wSid_n = new TH2F(Form("dE_IMnpipi_wSid_n"),Form("dE_IMnpipi_wSid_n"),nbinIMnpipi,IMnpipilow,IMnpipihi, nbindE, 0, 50);
  //dE_IMnpipi_wSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  //dE_IMnpipi_wSid_n->SetYTitle("dE [MeVee]");

  //dE_IMnpipi_woK0_wSid_n = new TH2F(Form("dE_IMnpipi_woK0_wSid_n"),Form("dE_IMnpipi_woK0_wSid_n"),nbinIMnpipi,IMnpipilow,IMnpipihi, nbindE, 0, 50);
  //dE_IMnpipi_woK0_wSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  //dE_IMnpipi_woK0_wSid_n->SetYTitle("dE [MeVee]");

  Cosn_IMnpipi_wSid_n = new TH2F(Form("Cosn_IMnpipi_wSid_n"),Form("dE_Cosn_IMnpipi_wSid_n"),100, 1, 2, 50, -1, 1);
  Cosn_IMnpipi_wSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  Cosn_IMnpipi_wSid_n->SetYTitle("cos#theta_{n} (CM)");

  Cosn_IMnpipi_woK0_wSid_n = new TH2F(Form("Cosn_IMnpipi_woK0_wSid_n"),Form("dE_Cosn_IMnpipi_woK0_wSid_n"),100, 1, 2, 50, -1, 1);
  Cosn_IMnpipi_woK0_wSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  Cosn_IMnpipi_woK0_wSid_n->SetYTitle("cos#theta_{n} (CM)");

  MMnmiss_IMnpipi_wSid = new TH2F("MMnmiss_IMnpipi_wSid","MMnmiss_IMnpipi_wSid",nbinIMnpipi,IMnpipilow,IMnpipihi,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMnpipi_wSid->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMnpipi_wSid->SetYTitle("Miss. Mass. [GeV/c^{2}]");

  MMnmiss_IMnpipi_wK0_wSid = new TH2F("MMnmiss_IMnpipi_wK0_wSid","MMnmiss_IMnpipi_wK0_wSid",nbinIMnpipi,IMnpipilow,IMnpipihi,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMnpipi_wK0_wSid->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMnpipi_wK0_wSid->SetYTitle("Miss. Mass. [GeV/c^{2}]");
  
  MMnmiss_IMnpipi_woK0_wSid = new TH2F("MMnmiss_IMnpipi_woK0_wSid","MMnmiss_IMnpipi_woK0_wSid",nbinIMnpipi,IMnpipilow,IMnpipihi,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMnpipi_woK0_wSid->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMnpipi_woK0_wSid->SetYTitle("Miss. Mass. [GeV/c^{2}]");

  MMnmiss_IMnpipi_woK0_wSid_Sp = new TH2F("MMnmiss_IMnpipi_woK0_wSid_Sp","MMnmiss_IMnpipi_woK0_wSid_Sp",nbinIMnpipi,IMnpipilow,IMnpipihi,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMnpipi_woK0_wSid_Sp->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMnpipi_woK0_wSid_Sp->SetYTitle("Miss. Mass. [GeV/c^{2}]");

  MMnmiss_IMnpipi_woK0_wSid_Sm = new TH2F("MMnmiss_IMnpipi_woK0_wSid_Sm","MMnmiss_IMnpipi_woK0_wSid_Sm",nbinIMnpipi,IMnpipilow,IMnpipihi,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMnpipi_woK0_wSid_Sm->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMnpipi_woK0_wSid_Sm->SetYTitle("Miss. Mass. [GeV/c^{2}]");

  MMnmiss_IMnpipi_wK0_woSid_won = new TH2F("MMnmiss_IMnpipi_wK0_woSid_won","MMnmiss_IMnpipi_wK0_woSid_won",nbinIMnpipi,IMnpipilow,IMnpipihi,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMnpipi_wK0_woSid_won->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMnpipi_wK0_woSid_won->SetYTitle("Miss. Mass. [GeV/c^{2}]");

  MMnmiss_Momnpipi_wK0_woSid_won = new TH2F("MMnmiss_Momnpipi_wK0_woSid_won","MMnmiss_Momnpipi_wK0_woSid_won",200,0,2,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_Momnpipi_wK0_woSid_won->SetXTitle("Mom.(n#pi^{+}#pi^{-}) [GeV/c]");
  MMnmiss_Momnpipi_wK0_woSid_won->SetYTitle("Miss. Mass. [GeV/c^{2}]");

  MMnmiss_IMnpipi_wK0_woSidn_won = new TH2F("MMnmiss_IMnpipi_wK0_woSidn_won","MMnmiss_IMnpipi_wK0_woSidn_won",nbinIMnpipi,IMnpipilow,IMnpipihi,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMnpipi_wK0_woSidn_won->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMnpipi_wK0_woSidn_won->SetYTitle("Miss. Mass. [GeV/c^{2}]");

  //MMnmiss_IMnpipi_woK0_woSidn = new TH2F("MMnmiss_IMnpipi_woK0_woSidn","MMnmiss_IMnpipi_woK0_woSidn",nbinIMnpipi,IMnpipilow,IMnpipihi,nbinnmiss, nmisslow, nmisshigh);
  //MMnmiss_IMnpipi_woK0_woSidn->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  //MMnmiss_IMnpipi_woK0_woSidn->SetYTitle("Miss. Mass. [GeV/c^{2}]");

  MMnmiss_IMnpipi_woK0_woSid_won = new TH2F("MMnmiss_IMnpipi_woK0_woSid_won","MMnmiss_IMnpipi_woK0_woSid_won",nbinIMnpipi,IMnpipilow,IMnpipihi,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_IMnpipi_woK0_woSid_won->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMnpipi_woK0_woSid_won->SetYTitle("Miss. Mass. [GeV/c^{2}]");

  MMnmiss_Momnpipi_woK0_woSid_won = new TH2F("MMnmiss_Momnpipi_woK0_woSid_won","MMnmiss_Momnpipi_woK0_woSid_won",200,0,2,nbinnmiss, nmisslow, nmisshigh);
  MMnmiss_Momnpipi_woK0_woSid_won->SetXTitle("Mom.(n#pi^{+}#pi^{-}) [GeV/c]");
  MMnmiss_Momnpipi_woK0_woSid_won->SetYTitle("Miss. Mass. [GeV/c^{2}]");

  q_IMnpipi_wSid_n = new TH2F("q_IMnpipi_wSid_n","q_IMnpipi_wSid_n",nbinIMnpipi,IMnpipilow,IMnpipihi,nbinq,0,1.5);
  q_IMnpipi_wSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wSid_n->SetYTitle("Mom. Transfer [GeV/c]");

  
  for(int ithcut=0;ithcut<nthetacut;ithcut++){
    int cuttheta = 5*(ithcut+1);
    q_IMnpipi_wSid_n_thetacut[ithcut] = new TH2F(Form("q_IMnpipi_wSid_n_thetacut%d",ithcut),Form("q_IMnpipi_wSid_n #theta_{lab} below %d^{#circ}",cuttheta),
                                                nbinIMnpipi,IMnpipilow,IMnpipihi,nbinq,0,1.5); 
    q_IMnpipi_wSid_n_thetacut[ithcut]->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");      
    q_IMnpipi_wSid_n_thetacut[ithcut]->SetYTitle("Mom. Transfer [GeV/c]");                
  }

  q_IMnpipi_wSid_n_fake = new TH2F("q_IMnpipi_wSid_n_fake","q_IMnpipi_wSid_n_fake",nbinIMnpipi,IMnpipilow,IMnpipihi,nbinq,0,1.5);
  q_IMnpipi_wSid_n_fake->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wSid_n_fake->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_IMnpipi_wSid_n_fake_pat2 = new TH2F("q_IMnpipi_wSid_n_fake_pat2","q_IMnpipi_wSid_n_fake_pat2",nbinIMnpipi,IMnpipilow,IMnpipihi,nbinq,0,1.5);
  q_IMnpipi_wSid_n_fake_pat2->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wSid_n_fake_pat2->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_IMnpipi_wSid_n_fake_pat7 = new TH2F("q_IMnpipi_wSid_n_fake_pat7","q_IMnpipi_wSid_n_fake_pat7",nbinIMnpipi,IMnpipilow,IMnpipihi,nbinq,0,1.5);
  q_IMnpipi_wSid_n_fake_pat7->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wSid_n_fake_pat7->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_IMnpipi_wSid_n_wocross = new TH2F("q_IMnpipi_wSid_n_wocross","q_IMnpipi_wSid_n_wocross",nbinIMnpipi,IMnpipilow,IMnpipihi,nbinq,0,1.5);
  q_IMnpipi_wSid_n_wocross->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wSid_n_wocross->SetYTitle("Mom. Transfer [GeV/c]");

  q_IMnpipi_woK0_wSid_n = new TH2F("q_IMnpipi_woK0_wSid_n","q_IMnpipi_woK0_wSid_n", nbinIMnpipi,IMnpipilow,IMnpipihi, nbinq,0,1.5);
  q_IMnpipi_woK0_wSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_wSid_n->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_IMnpipi_woK0_wSid_n_woSp = new TH2F("q_IMnpipi_woK0_wSid_n_woSp","q_IMnpipi_woK0_wSid_n_woSp", nbinIMnpipi,IMnpipilow,IMnpipihi, nbinq,0,1.5);
  q_IMnpipi_woK0_wSid_n_woSp->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_wSid_n_woSp->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_IMnpipi_woK0_wSid_n_woSm = new TH2F("q_IMnpipi_woK0_wSid_n_woSm","q_IMnpipi_woK0_wSid_n_woSm", nbinIMnpipi,IMnpipilow,IMnpipihi, nbinq,0,1.5);
  q_IMnpipi_woK0_wSid_n_woSm->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_wSid_n_woSm->SetYTitle("Mom. Transfer [GeV/c]");

  q_IMnpipi_wK0_wSid_n = new TH2F("q_IMnpipi_wK0_wSid_n","q_IMnpipi_wK0_wSid_n", nbinIMnpipi,IMnpipilow,IMnpipihi, nbinq,0,1.5);
  q_IMnpipi_wK0_wSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wK0_wSid_n->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_IMnpipi_wK0orwSid_n = new TH2F("q_IMnpipi_wK0orwSid_n","q_IMnpipi_wK0orwSid_n", nbinIMnpipi,IMnpipilow,IMnpipihi, nbinq,0,1.5);
  q_IMnpipi_wK0orwSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wK0orwSid_n->SetYTitle("Mom. Transfer [GeV/c]");

  q_IMnpipi_wK0_woSid_n = new TH2F("q_IMnpipi_wK0_woSid_n","q_IMnpipi_wK0_woSid_n", nbinIMnpipi,IMnpipilow,IMnpipihi, nbinq,0,1.5);
  q_IMnpipi_wK0_woSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wK0_woSid_n->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_IMpiSigma_gen = new TH2F(Form("q_IMpiSigma_gen"),Form("q_IMpiSigma_gen"),nbinIMnpipi,IMnpipilow,IMnpipihi,nbinq,0,1.5);
  if(SimSpmode) {
    q_IMpiSigma_gen->SetXTitle("true IM(#Sigma^{+}#pi^{-}) [GeV/c^{2}]");
  } else if(SimSmmode) {
    q_IMpiSigma_gen->SetXTitle("true IM(#Sigma^{-}#pi^{+}) [GeV/c^{2}]");
  }else if(SimK0nnmode){
    q_IMpiSigma_gen->SetXTitle("true IM(K^{0}n) [GeV/c^{2}]");
  }
  q_IMpiSigma_gen->SetYTitle("true Mom. Transfer [GeV/c]");

  q_IMpiSigma_wSid_n_genacc = new TH2F(Form("q_IMpiSigma_wSid_n_genacc"),Form("q_IMpiSigma_wSid_n_genacc"),nbinIMnpipi,IMnpipilow,IMnpipihi,nbinq,0,1.5);
  if(SimSpmode) {
    q_IMpiSigma_wSid_n_genacc->SetXTitle("true IM(#Sigma^{+}#pi^{-}) [GeV/c^{2}]");
  } else if(SimSmmode) {
    q_IMpiSigma_wSid_n_genacc->SetXTitle("true IM(#Sigma^{-}#pi^{+}) [GeV/c^{2}]");
  }else if(SimK0nnmode){
    q_IMpiSigma_wSid_n_genacc->SetXTitle("true IM(K^{0}n) [GeV/c^{2}]");
  }
  q_IMpiSigma_wSid_n_genacc->SetYTitle("true Mom. Transfer [GeV/c]");

  q_IMnpipi_wSid_n_acc = new TH2F(Form("q_IMnpipi_wSid_n_acc"),Form("q_IMnpipi_wSid_n_acc"),nbinIMnpipi,IMnpipilow,IMnpipihi,nbinq,0,1.5);
  q_IMnpipi_wSid_n_acc->SetXTitle("true IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wSid_n_acc->SetYTitle("true Mom. Transfer [GeV/c]");

  q_IMnpipi_wSid_n_acc_reco = new TH2F(Form("q_IMnpipi_wSid_n_acc_reco"),Form("q_IMnpipi_wSid_n_acc_reco"),nbinIMnpipi,IMnpipilow,IMnpipihi,nbinq,0,1.5);
  q_IMnpipi_wSid_n_acc_reco->SetXTitle("reco. IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wSid_n_acc_reco->SetYTitle("reco. Mom. Transfer [GeV/c]");

  q_IMpiSigma_woK0_wSid_n_genacc = new TH2F(Form("q_IMpiSigma_woK0_wSid_n_genacc"),Form("q_IMpiSigma_woK0_wSid_n_genacc"),nbinIMnpipi,IMnpipilow,IMnpipihi,nbinq,0,1.5);
  if(SimSpmode) {
    q_IMpiSigma_woK0_wSid_n_genacc->SetXTitle("true IM(#Sigma^{+}#pi^{-}) [GeV/c^{2}]");
  } else if(SimSmmode) {
    q_IMpiSigma_woK0_wSid_n_genacc->SetXTitle("true IM(#Sigma^{-}#pi^{+}) [GeV/c^{2}]");
  } else if(SimK0nnmode){
    q_IMpiSigma_woK0_wSid_n_genacc->SetXTitle("true IM(K^{0}n) [GeV/c^{2}]");
  }
  q_IMpiSigma_woK0_wSid_n_genacc->SetYTitle("true Mom. Transfer [GeV/c]");

  q_IMnpipi_woK0_wSid_n_acc = new TH2F(Form("q_IMnpipi_woK0_wSid_n_acc"),Form("q_IMnpipi_woK0_wSid_n_acc"),nbinIMnpipi,IMnpipilow,IMnpipihi,nbinq,0,1.5);
  q_IMnpipi_woK0_wSid_n_acc->SetXTitle("true IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_wSid_n_acc->SetYTitle("true Mom. Transfer [GeV/c]");

  q_IMnpipi_woK0_wSid_n_acc_reco = new TH2F(Form("q_IMnpipi_woK0_wSid_n_acc_reco"),Form("q_IMnpipi_woK0_wSid_n_acc_reco"),nbinIMnpipi,IMnpipilow,IMnpipihi,nbinq,0,1.5);
  q_IMnpipi_woK0_wSid_n_acc_reco->SetXTitle("reco. IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_wSid_n_acc_reco->SetYTitle("reco. Mom. Transfer [GeV/c]");

  q_IMnpipi_wSid_n_Sp = new TH2F("q_IMnpipi_wSid_n_Sp","q_IMnpipi_wSid_n_Sp", nbinIMnpipi,IMnpipilow,IMnpipihi, nbinq,0,1.5);
  q_IMnpipi_wSid_n_Sp->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wSid_n_Sp->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_IMnpipi_wSid_n_Sp_Stop = new TH2F("q_IMnpipi_wSid_n_Sp_Stop","q_IMnpipi_wSid_n_Sp_Stop Mom(#Sigma^{+}) < 30 MeV", nbinIMnpipi,IMnpipilow,IMnpipihi, nbinq,0,1.5);
  q_IMnpipi_wSid_n_Sp_Stop->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wSid_n_Sp_Stop->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_IMnpipi_wSid_n_Sp_NoStop = new TH2F("q_IMnpipi_wSid_n_Sp_NoStop","q_IMnpipi_wSid_n_Sp_NoStop Mom(#Sigma^{+}) >= 30 MeV", nbinIMnpipi,IMnpipilow,IMnpipihi, nbinq,0,1.5);
  q_IMnpipi_wSid_n_Sp_NoStop->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wSid_n_Sp_NoStop->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_IMnpipi_wSid_n_SpSm = new TH2F("q_IMnpipi_wSid_n_SpSm","q_IMnpipi_wSid_n_SpSm", nbinIMnpipi,IMnpipilow,IMnpipihi, nbinq,0,1.5);
  q_IMnpipi_wSid_n_SpSm->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wSid_n_SpSm->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_IMnpipi_wSid_n_Sp_woSm = new TH2F("q_IMnpipi_wSid_n_Sp_woSm","q_IMnpipi_wSid_n_Sp_woSm", nbinIMnpipi,IMnpipilow,IMnpipihi, nbinq,0,1.5);
  q_IMnpipi_wSid_n_Sp_woSm->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wSid_n_Sp_woSm->SetYTitle("Mom. Transfer [GeV/c]");

  q_IMnpipi_wK0_n = new TH2F("q_IMnpipi_wK0_n","q_IMnpipi_wK0_n", nbinIMnpipi,IMnpipilow,IMnpipihi, nbinq,0,1.5);
  q_IMnpipi_wK0_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wK0_n->SetYTitle("Mom. Transfer [GeV/c]");

  q_IMnpipi_woK0_woSid = new TH2F("q_IMnpipi_woK0_woSid","q_IMnpipi_woK0_woSid", nbinIMnpipi,IMnpipilow,IMnpipihi, nbinq,0,1.5);
  q_IMnpipi_woK0_woSid->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_woSid->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_IMnpipi_woK0_woSid_n = new TH2F("q_IMnpipi_woK0_woSid_n","q_IMnpipi_woK0_woSid_n", nbinIMnpipi,IMnpipilow,IMnpipihi, nbinq,0,1.5);
  q_IMnpipi_woK0_woSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_woSid_n->SetYTitle("Mom. Transfer [GeV/c]");

  q_IMnpipi_woK0_woSid_won = new TH2F(Form("q_IMnpipi_woK0_woSid_won"),Form("q_IMnpipi_woK0_woSid_won"), nbinIMnpipi,IMnpipilow,IMnpipihi, nbinq,0,1.5);
  q_IMnpipi_woK0_woSid_won->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_woSid_won->SetYTitle("Mom. Transfer [GeV/c]");
  
  //q_nmom_woSid_won = new TH2F("q_nmom_woSid_won","q_nmom_woSid_won", nbinnmom,0,1, nbinq,0,1.5);
  //q_nmom_woSid_won->SetXTitle("nmom [GeV/c]");
  //q_nmom_woSid_won->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_nmom_woK0_woSid_won = new TH2F("q_nmom_woK0_woSid_won","q_nmom_woK0_woSid_won", nbinnmom,0,1, nbinq,0,1.5);
  q_nmom_woK0_woSid_won->SetXTitle("nmom [GeV/c]");
  q_nmom_woK0_woSid_won->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_nmom_wSid_n = new TH2F("q_nmom_wSid_n","q_nmom_wSid_n", nbinnmom,0,1, nbinq,0,1.5);
  q_nmom_wSid_n->SetXTitle("nmom [GeV/c]");
  q_nmom_wSid_n->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_nmom_woK0_wSid_n = new TH2F("q_nmom_woK0_wSid_n","q_nmom_woK0_wSid_n", nbinnmom,0,1, nbinq,0,1.5);
  q_nmom_woK0_wSid_n->SetXTitle("nmom [GeV/c]");
  q_nmom_woK0_wSid_n->SetYTitle("Mom. Transfer [GeV/c]");
  
  for(unsigned int iwbin=0;iwbin<nwbin;iwbin++){
    q_nmom_woK0_wSid_n_woSm_wbin[iwbin] = new TH2F(
        Form("q_nmom_woK0_wSid_n_woSm_wbin%d",iwbin),
        Form("q_nmom_woK0_wSid_n_woSm %0.2f-%0.2f",wbinlow[iwbin],wbinhigh[iwbin]),
        nbinnmom,0,1, nbinq,0,1.5);
    q_nmom_woK0_wSid_n_woSm_wbin[iwbin]->SetXTitle("nmom [GeV/c]");
    q_nmom_woK0_wSid_n_woSm_wbin[iwbin]->SetYTitle("Mom. Transfer [GeV/c]");
    
    q_nmom_woK0_wSid_n_woSp_wbin[iwbin] = new TH2F(
        Form("q_nmom_woK0_wSid_n_woSp_wbin%d",iwbin),
        Form("q_nmom_woK0_wSid_n_woSp %0.2f-%0.2f",wbinlow[iwbin],wbinhigh[iwbin]),
        nbinnmom,0,1, nbinq,0,1.5);
    q_nmom_woK0_wSid_n_woSp_wbin[iwbin]->SetXTitle("nmom [GeV/c]");
    q_nmom_woK0_wSid_n_woSp_wbin[iwbin]->SetYTitle("Mom. Transfer [GeV/c]");
  }


  q_nmom_wK0_wSid_n = new TH2F("q_nmom_wK0_wSid_n","q_nmom_wK0_wSid_n", nbinnmom,0,1, nbinq,0,1.5);
  q_nmom_wK0_wSid_n->SetXTitle("nmom [GeV/c]");
  q_nmom_wK0_wSid_n->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_IMnpipi_wK0_wSid_n_SpSm = new TH2F(Form("q_IMnpipi_wK0_wSid_n_SpSm"),Form("q_IMnpipi_wK0_wSid_n_SpSm"), nbinIMnpipi,IMnpipilow,IMnpipihi, nbinq,0,1.5);
  q_IMnpipi_wK0_wSid_n_SpSm->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wK0_wSid_n_SpSm->SetYTitle("Mom. Transfer [GeV/c]");

  q_IMnpipi_wK0_wSid_n_Sp = new TH2F(Form("q_IMnpipi_wK0_wSid_n_Sp"),Form("q_IMnpipi_wK0_wSid_n_Sp"), nbinIMnpipi,IMnpipilow,IMnpipihi, nbinq,0,1.5);
  q_IMnpipi_wK0_wSid_n_Sp->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wK0_wSid_n_Sp->SetYTitle("Mom. Transfer [GeV/c]");

  q_IMnpipi_woK0_wSid_n_Sp = new TH2F(Form("q_IMnpipi_woK0_wSid_n_Sp"),Form("q_IMnpipi_woK0_wSid_n_Sp"), nbinIMnpipi,IMnpipilow,IMnpipihi, nbinq,0,1.5);
  q_IMnpipi_woK0_wSid_n_Sp->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_wSid_n_Sp->SetYTitle("Mom. Transfer [GeV/c]");

  q_IMpiSigma_wSid_n_Sp_genacc = new TH2F(Form("q_IMpiSigma_wSid_n_Sp_genacc"),Form("q_IMpiSigma_wSid_n_Sp_genacc"),nbinIMnpipi,IMnpipilow,IMnpipihi,nbinq,0,1.5);
  if(SimSpmode) {
    q_IMpiSigma_wSid_n_Sp_genacc->SetXTitle("true IM(#Sigma^{+}#pi^{-}) [GeV/c^{2}]");
  } else if(SimSmmode) {
    q_IMpiSigma_wSid_n_Sp_genacc->SetXTitle("true IM(#Sigma^{-}#pi^{+}) [GeV/c^{2}]");
  }

  q_IMpiSigma_woK0_wSid_n_Sp_genacc = new TH2F(Form("q_IMpiSigma_woK0_wSid_n_Sp_genacc"),Form("q_IMpiSigma_woK0_wSid_n_Sp_genacc"),nbinIMnpipi,IMnpipilow,IMnpipihi,nbinq,0,1.5);
  if(SimSpmode) {
    q_IMpiSigma_woK0_wSid_n_Sp_genacc->SetXTitle("true IM(#Sigma^{+}#pi^{-}) [GeV/c^{2}]");
  } else if(SimSmmode) {
    q_IMpiSigma_woK0_wSid_n_Sp_genacc->SetXTitle("true IM(#Sigma^{-}#pi^{+}) [GeV/c^{2}]");
  }
  q_IMpiSigma_woK0_wSid_n_Sp_genacc->SetYTitle("true Mom. Transfer [GeV/c]");

  q_IMnpipi_wSid_n_Sp_acc = new TH2F(Form("q_IMnpipi_wSid_n_Sp_acc"),Form("q_IMnpipi_wSid_n_Sp_acc"),nbinIMnpipi,IMnpipilow,IMnpipihi,nbinq,0,1.5);
  q_IMnpipi_wSid_n_Sp_acc->SetXTitle("true IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wSid_n_Sp_acc->SetYTitle("true Mom. Transfer [GeV/c]");

  q_IMnpipi_woK0_wSid_n_Sp_acc = new TH2F(Form("q_IMnpipi_woK0_wSid_n_Sp_acc"),Form("q_IMnpipi_woK0_wSid_n_Sp_acc"),nbinIMnpipi,IMnpipilow,IMnpipihi,nbinq,0,1.5);
  q_IMnpipi_woK0_wSid_n_Sp_acc->SetXTitle("true IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_wSid_n_Sp_acc->SetYTitle("true Mom. Transfer [GeV/c]");

  q_IMnpipi_wSid_n_Sp_acc_reco = new TH2F(Form("q_IMnpipi_wSid_n_Sp_acc_reco"),Form("q_IMnpipi_wSid_n_Sp_acc_reco"),nbinIMnpipi,IMnpipilow,IMnpipihi,nbinq,0,1.5);
  q_IMnpipi_wSid_n_Sp_acc_reco->SetXTitle("reco. IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wSid_n_Sp_acc_reco->SetYTitle("reco. Mom. Transfer [GeV/c]");

  q_IMnpipi_woK0_wSid_n_Sp_acc_reco = new TH2F(Form("q_IMnpipi_woK0_wSid_n_Sp_acc_reco"),Form("q_IMnpipi_woK0_wSid_n_Sp_acc_reco"),nbinIMnpipi,IMnpipilow,IMnpipihi,nbinq,0,1.5);
  q_IMnpipi_woK0_wSid_n_Sp_acc_reco->SetXTitle("reco. IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_wSid_n_Sp_acc_reco->SetYTitle("reco. Mom. Transfer [GeV/c]");
  
  /*
  for(int itype=0; itype<3; itype++) {
    for(int ilh=0; ilh<2; ilh++) {
      for(int igap=0; igap<ngap; igap++) {
        const char  lh[][6]= {"low","high"};
        q_IMnpipi_wSid_n_Sp_side[itype][ilh][igap] = new TH2F(Form("q_IMnpipi_wSid_n_Sp_side_%d_%s_%d",itype,lh[ilh],igap),Form("q_IMnpipi_wSid_n_Sp_side_%d_%s_%d",itype,lh[ilh],igap), nbinIMnpipi,IMnpipilow,IMnpipihi, nbinq,0,1.5);
        q_IMnpipi_wSid_n_Sp_side[itype][ilh][igap]->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
        q_IMnpipi_wSid_n_Sp_side[itype][ilh][igap]->SetYTitle("Mom. Transfer [GeV/c]");

        q_IMnpipi_woK0_wSid_n_Sp_side[itype][ilh][igap] = new TH2F(Form("q_IMnpipi_woK0_wSid_n_Sp_side_%d_%s_%d",itype,lh[ilh],igap),Form("q_IMnpipi_woK0_wSid_n_Sp_side_%d_%s_%d",itype,lh[ilh],igap), nbinIMnpipi,IMnpipilow,IMnpipihi, nbinq,0,1.5);
        q_IMnpipi_woK0_wSid_n_Sp_side[itype][ilh][igap]->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
        q_IMnpipi_woK0_wSid_n_Sp_side[itype][ilh][igap]->SetYTitle("Mom. Transfer [GeV/c]");
      }
    }
  }*/
  
  /*
  for(int izone=0; izone<nzone; izone++) {
    q_IMnpipi_wSid_n_Sp_sidewide[izone] = new TH2F(Form("q_IMnpipi_wSid_n_Sp_sidewide_%d",izone),Form("q_IMnpipi_wSid_n_Sp_sidewide_%d",izone), nbinIMnpipi,IMnpipilow,IMnpipihi, nbinq,0,1.5);
    q_IMnpipi_wSid_n_Sp_sidewide[izone]->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
    q_IMnpipi_wSid_n_Sp_sidewide[izone]->SetYTitle("Mom. Transfer [GeV/c]");
    q_IMnpipi_woK0_wSid_n_Sp_sidewide[izone] = new TH2F(Form("q_IMnpipi_woK0_wSid_n_Sp_sidewide_%d",izone),Form("q_IMnpipi_woK0_wSid_n_Sp_sidewide_%d",izone), nbinIMnpipi,IMnpipilow,IMnpipihi, nbinq,0,1.5);
    q_IMnpipi_woK0_wSid_n_Sp_sidewide[izone]->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
    q_IMnpipi_woK0_wSid_n_Sp_sidewide[izone]->SetYTitle("Mom. Transfer [GeV/c]");
    q_IMnpipi_wSid_n_Sm_sidewide[izone] = new TH2F(Form("q_IMnpipi_wSid_n_Sm_sidewide_%d",izone),Form("q_IMnpipi_wSid_n_Sm_sidewide_%d",izone), nbinIMnpipi,IMnpipilow,IMnpipihi, nbinq,0,1.5);
    q_IMnpipi_wSid_n_Sm_sidewide[izone]->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
    q_IMnpipi_wSid_n_Sm_sidewide[izone]->SetYTitle("Mom. Transfer [GeV/c]");
    q_IMnpipi_woK0_wSid_n_Sm_sidewide[izone] = new TH2F(Form("q_IMnpipi_woK0_wSid_n_Sm_sidewide_%d",izone),Form("q_IMnpipi_woK0_wSid_n_Sm_sidewide_%d",izone), nbinIMnpipi,IMnpipilow,IMnpipihi, nbinq,0,1.5);
    q_IMnpipi_woK0_wSid_n_Sm_sidewide[izone]->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
    q_IMnpipi_woK0_wSid_n_Sm_sidewide[izone]->SetYTitle("Mom. Transfer [GeV/c]");
  }*/

  q_IMnpipi_wSid_n_Sm = new TH2F("q_IMnpipi_wSid_n_Sm","q_IMnpipi_wSid_n_Sm", nbinIMnpipi,IMnpipilow,IMnpipihi, nbinq,0,1.5);
  q_IMnpipi_wSid_n_Sm->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wSid_n_Sm->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_IMnpipi_wSid_n_Sm_Stop = new TH2F("q_IMnpipi_wSid_n_Sm_Stop","q_IMnpipi_wSid_n_Sm_Stop Mom(#Sigma^{-}) < 30 MeV", nbinIMnpipi,IMnpipilow,IMnpipihi, nbinq,0,1.5);
  q_IMnpipi_wSid_n_Sm_Stop->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wSid_n_Sm_Stop->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_IMnpipi_wSid_n_Sm_NoStop = new TH2F("q_IMnpipi_wSid_n_Sm_NoStop","q_IMnpipi_wSid_n_Sm_NoStop Mom(#Sigma^{-}) => 30 MeV", nbinIMnpipi,IMnpipilow,IMnpipihi, nbinq,0,1.5);
  q_IMnpipi_wSid_n_Sm_NoStop->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wSid_n_Sm_NoStop->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_IMnpipi_wSid_n_Sm_woSp = new TH2F("q_IMnpipi_wSid_n_Sm_woSp","q_IMnpipi_wSid_n_Sm_woSp", nbinIMnpipi,IMnpipilow,IMnpipihi, nbinq,0,1.5);
  q_IMnpipi_wSid_n_Sm_woSp->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wSid_n_Sm_woSp->SetYTitle("Mom. Transfer [GeV/c]");

  q_IMnpipi_wK0_wSid_n_Sm = new TH2F(Form("q_IMnpipi_wK0_wSid_n_Sm"),Form("q_IMnpipi_wK0_wSid_n_Sm"), nbinIMnpipi,IMnpipilow,IMnpipihi, nbinq,0,1.5);
  q_IMnpipi_wK0_wSid_n_Sm->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wK0_wSid_n_Sm->SetYTitle("Mom. Transfer [GeV/c]");

  q_IMnpipi_woK0_wSid_n_Sm = new TH2F(Form("q_IMnpipi_woK0_wSid_n_Sm"),Form("q_IMnpipi_woK0_wSid_n_Sm"), nbinIMnpipi,IMnpipilow,IMnpipihi, nbinq,0,1.5);
  q_IMnpipi_woK0_wSid_n_Sm->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_wSid_n_Sm->SetYTitle("Mom. Transfer [GeV/c]");

  q_IMnpipi_wK0_woSid_won = new TH2F("q_IMnpipi_wK0_woSid_won","q_IMnpipi_wK0_woSid_won", nbinIMnpipi,IMnpipilow,IMnpipihi, nbinq,0,1.5);
  q_IMnpipi_wK0_woSid_won->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wK0_woSid_won->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_nmom_wK0_woSid_won = new TH2F("q_nmom_wK0_woSid_won","q_nmom_wK0_woSid_won", nbinnmom,0,1, nbinq,0,1.5);
  q_nmom_wK0_woSid_won->SetXTitle("nmom [GeV/c]");
  q_nmom_wK0_woSid_won->SetYTitle("Mom. Transfer [GeV/c]");
  
  for(unsigned int iwbin=0;iwbin<nwbin;iwbin++){
    q_nmom_wK0_woSid_n_wbin[iwbin] = new TH2F(
         Form("q_nmom_wK0_woSid_n_wbin%d",iwbin),
         Form("q_nmom_wK0_woSid_n %0.2f-%0.2f",wbinlow[iwbin],wbinhigh[iwbin]),
         nbinnmom,0,1, nbinq,0,1.5);
    q_nmom_wK0_woSid_n_wbin[iwbin]->SetXTitle("nmom [GeV/c]");
    q_nmom_wK0_woSid_n_wbin[iwbin]->SetYTitle("Mom. Transfer [GeV/c]");
  }
  
  q_IMnpipi_wK0_woSidn_won = new TH2F("q_IMnpipi_wK0_woSidn_won","q_IMnpipi_wK0_woSidn_won", nbinIMnpipi,IMnpipilow,IMnpipihi, nbinq,0,1.5);
  q_IMnpipi_wK0_woSidn_won->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wK0_woSidn_won->SetYTitle("Mom. Transfer [GeV/c]");

  //q_IMnpipi_woK0_woSidn = new TH2F("q_IMnpipi_woK0_woSidn","q_IMnpipi_woK0_woSidn", nbinIMnpipi,IMnpipilow,IMnpipihi, nbinq,0,1.5);
  //q_IMnpipi_woK0_woSidn->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  //q_IMnpipi_woK0_woSidn->SetYTitle("Mom. Transfer [GeV/c]");

  //q_IMnpipi_woK0_woSidn_cross = new TH2F("q_IMnpipi_woK0_woSidn_cross","q_IMnpipi_woK0_woSidn_cross", nbinIMnpipi,IMnpipilow,IMnpipihi, nbinq,0,1.5);
  //q_IMnpipi_woK0_woSidn_cross->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  //q_IMnpipi_woK0_woSidn_cross->SetYTitle("Mom. Transfer [GeV/c]");

  //q_IMnpipi_woK0_woSidn_cross_Sp = new TH2F("q_IMnpipi_woK0_woSidn_cross_Sp","q_IMnpipi_woK0_woSidn_cross_Sp", nbinIMnpipi,IMnpipilow,IMnpipihi, nbinq,0,1.5);
  //q_IMnpipi_woK0_woSidn_cross_Sp->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  //q_IMnpipi_woK0_woSidn_cross_Sp->SetYTitle("Mom. Transfer [GeV/c]");

  //q_IMnpipi_woK0_woSidn_cross_Sm = new TH2F("q_IMnpipi_woK0_woSidn_cross_Sm","q_IMnpipi_woK0_woSidn_cross_Sm", nbinIMnpipi,IMnpipilow,IMnpipihi, nbinq,0,1.5);
  //q_IMnpipi_woK0_woSidn_cross_Sm->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  //q_IMnpipi_woK0_woSidn_cross_Sm->SetYTitle("Mom. Transfer [GeV/c]");

  q_IMpiSigma_wSid_n_Sm_genacc = new TH2F(Form("q_IMpiSigma_wSid_n_Sm_genacc"),Form("q_IMpiSigma_wSid_n_Sm_genacc"),nbinIMnpipi,IMnpipilow,IMnpipihi,nbinq,0,1.5);
  if(SimSpmode) {
    q_IMpiSigma_wSid_n_Sm_genacc->SetXTitle("true IM(#Sigma^{+}#pi^{-}) [GeV/c^{2}]");
  } else if(SimSmmode) {
    q_IMpiSigma_wSid_n_Sm_genacc->SetXTitle("true IM(#Sigma^{-}#pi^{+}) [GeV/c^{2}]");
  }
  q_IMpiSigma_wSid_n_Sm_genacc->SetYTitle("true Mom. Transfer [GeV/c]");

  q_IMpiSigma_woK0_wSid_n_Sm_genacc = new TH2F(Form("q_IMpiSigma_woK0_wSid_n_Sm_genacc"),Form("q_IMpiSigma_woK0_wSid_n_Sm_genacc"),nbinIMnpipi,IMnpipilow,IMnpipihi,nbinq,0,1.5);
  if(SimSpmode) {
    q_IMpiSigma_woK0_wSid_n_Sm_genacc->SetXTitle("true IM(#Sigma^{+}#pi^{-}) [GeV/c^{2}]");
  } else if(SimSmmode) {
    q_IMpiSigma_woK0_wSid_n_Sm_genacc->SetXTitle("true IM(#Sigma^{-}#pi^{+}) [GeV/c^{2}]");
  }
  q_IMpiSigma_woK0_wSid_n_Sm_genacc->SetYTitle("true Mom. Transfer [GeV/c]");

  q_IMnpipi_wSid_n_Sm_acc = new TH2F(Form("q_IMnpipi_wSid_n_Sm_acc"),Form("q_IMnpipi_wSid_n_Sm_acc"),nbinIMnpipi,IMnpipilow,IMnpipihi,nbinq,0,1.5);
  q_IMnpipi_wSid_n_Sm_acc->SetXTitle("true IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wSid_n_Sm_acc->SetYTitle("true Mom. Transfer [GeV/c]");

  q_IMnpipi_woK0_wSid_n_Sm_acc = new TH2F(Form("q_IMnpipi_woK0_wSid_n_Sm_acc"),Form("q_IMnpipi_woK0_wSid_n_Sm_acc"),nbinIMnpipi,IMnpipilow,IMnpipihi,nbinq,0,1.5);
  q_IMnpipi_woK0_wSid_n_Sm_acc->SetXTitle("true IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_wSid_n_Sm_acc->SetYTitle("true Mom. Transfer [GeV/c]");

  q_IMnpipi_wSid_n_Sm_acc_reco = new TH2F(Form("q_IMnpipi_wSid_n_Sm_acc_reco"),Form("q_IMnpipi_wSid_n_Sm_acc_reco"),nbinIMnpipi,IMnpipilow,IMnpipihi,nbinq,0,1.5);
  q_IMnpipi_wSid_n_Sm_acc_reco->SetXTitle("reco. IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wSid_n_Sm_acc_reco->SetYTitle("reco. Mom. Transfer [GeV/c]");

  q_IMnpipi_woK0_wSid_n_Sm_acc_reco = new TH2F(Form("q_IMnpipi_woK0_wSid_n_Sm_acc_reco"),Form("q_IMnpipi_woK0_wSid_n_Sm_acc_reco"),nbinIMnpipi,IMnpipilow,IMnpipihi,nbinq,0,1.5);
  q_IMnpipi_woK0_wSid_n_Sm_acc_reco->SetXTitle("reco. IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_wSid_n_Sm_acc_reco->SetYTitle("reco. Mom. Transfer [GeV/c]");
  
  /*
  for(int itype=0; itype<3; itype++) {
    for(int ilh=0; ilh<2; ilh++) {
      for(int igap=0; igap<ngap; igap++) {
        const char  lh[][6]= {"low","high"};
        q_IMnpipi_wSid_n_Sm_side[itype][ilh][igap] = new TH2F(Form("q_IMnpipi_wSid_n_Sm_side_%d_%s_%d",itype,lh[ilh],igap),Form("q_IMnpipi_wSid_n_Sm_side_%d_%s_%d",itype,lh[ilh],igap),nbinIMnpipi,IMnpipilow,IMnpipihi,nbinq,0,1.5);
        q_IMnpipi_wSid_n_Sm_side[itype][ilh][igap]->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
        q_IMnpipi_wSid_n_Sm_side[itype][ilh][igap]->SetYTitle("Mom. Transfer [GeV/c]");
        q_IMnpipi_woK0_wSid_n_Sm_side[itype][ilh][igap] = new TH2F(Form("q_IMnpipi_woK0_wSid_n_Sm_side_%d_%s_%d",itype,lh[ilh],igap),Form("q_IMnpipi_woK0_wSid_n_Sm_side_%d_%s_%d",itype,lh[ilh],igap),nbinIMnpipi,IMnpipilow,IMnpipihi,nbinq,0,1.5);
        q_IMnpipi_woK0_wSid_n_Sm_side[itype][ilh][igap]->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
        q_IMnpipi_woK0_wSid_n_Sm_side[itype][ilh][igap]->SetYTitle("Mom. Transfer [GeV/c]");
      }
    }
  }*/
  
  /*
  for(int igap=0; igap<ngap; igap++) {
    q_IMnpipi_wSid_n_side[igap] = new TH2F(Form("q_IMnpipi_wSid_n_side_%d",igap),Form("q_IMnpipi_wSid_n_side_%d",igap), nbinIMnpipi,IMnpipilow,IMnpipihi, nbinq,0,1.5);
    q_IMnpipi_wSid_n_side[igap]->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
    q_IMnpipi_wSid_n_side[igap]->SetYTitle("Mom. Transfer [GeV/c]");
  }

  for(int igap=0; igap<ngap; igap++) {
    q_IMnpipi_woK0_wSid_n_side[igap] = new TH2F(Form("q_IMnpipi_woK0_wSid_n_side_%d",igap),Form("q_IMnpipi_woK0_wSid_n_side_%d",igap), nbinIMnpipi,IMnpipilow,IMnpipihi, nbinq,0,1.5);
    q_IMnpipi_woK0_wSid_n_side[igap]->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
    q_IMnpipi_woK0_wSid_n_side[igap]->SetYTitle("Mom. Transfer [GeV/c]");
  }*/


  IMnpip_IMnpipi_n = new TH2F(Form("IMnpip_IMnpipi_n"),Form("IMnpip_IMnpipi_n"),nbinIMnpipi,IMnpipilow,IMnpipihi, nbinIMnpi,1.0,2.0);
  IMnpip_IMnpipi_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  IMnpip_IMnpipi_n->SetYTitle("IM(n#pi^{+}) [GeV/c^{2}]");

  IMnpip_IMnpipi_wK0_n = new TH2F(Form("IMnpip_IMnpipi_wK0_n"),Form("IMnpip_IMnpipi_wK0_n"),nbinIMnpipi,IMnpipilow,IMnpipihi, nbinIMnpi,1.0,2.0);
  IMnpip_IMnpipi_wK0_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  IMnpip_IMnpipi_wK0_n->SetYTitle("IM(n#pi^{+}) [GeV/c^{2}]");

  IMnpip_IMnpipi_woK0_n = new TH2F(Form("IMnpip_IMnpipi_woK0_n"),Form("IMnpip_IMnpipi_woK0_n"),nbinIMnpipi,IMnpipilow,IMnpipihi, nbinIMnpi,1.0,2.0);
  IMnpip_IMnpipi_woK0_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  IMnpip_IMnpipi_woK0_n->SetYTitle("IM(n#pi^{+}) [GeV/c^{2}]");

  IMnpim_IMnpipi_n = new TH2F(Form("IMnpim_IMnpipi_n"),Form("IMnpim_IMnpipi_n"),nbinIMnpipi,IMnpipilow,IMnpipihi, nbinIMnpi,1.0,2.0);
  IMnpim_IMnpipi_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  IMnpim_IMnpipi_n->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");

  IMnpim_IMnpipi_wK0_n = new TH2F(Form("IMnpim_IMnpipi_wK0_n"),Form("IMnpim_IMnpipi_wK0_n"),nbinIMnpipi,IMnpipilow,IMnpipihi, nbinIMnpi,1.0,2.0);
  IMnpim_IMnpipi_wK0_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  IMnpim_IMnpipi_wK0_n->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");

  IMnpim_IMnpipi_woK0_n = new TH2F(Form("IMnpim_IMnpipi_woK0_n"),Form("IMnpim_IMnpipi_woK0_n"),nbinIMnpipi,IMnpipilow,IMnpipihi, nbinIMnpi,1.0,2.0);
  IMnpim_IMnpipi_woK0_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  IMnpim_IMnpipi_woK0_n->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");

  pipmom_IMnpipi_wSid_n = new TH2F("pipmom_IMnpipi_wSid_n","pipmom_IMnpipi_wSid_n",nbinIMnpipi,IMnpipilow,IMnpipihi,200,0,1.0);
  pipmom_IMnpipi_wSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  pipmom_IMnpipi_wSid_n->SetYTitle("Mom(#pi^{+}) [GeV/c] ");

  pimmom_IMnpipi_wSid_n = new TH2F("pimmom_IMnpipi_wSid_n","pimmom_IMnpipi_wSid_n",nbinIMnpipi,IMnpipilow,IMnpipihi,200,0,1.0);
  pimmom_IMnpipi_wSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  pimmom_IMnpipi_wSid_n->SetYTitle("Mom(#pi^{-}) [GeV/c] ");

  pipmom_IMnpipi_woK0_wSid_n = new TH2F("pipmom_IMnpipi_woK0_wSid_n","pipmom_IMnpipi_woK0_wSid_n",nbinIMnpipi,IMnpipilow,IMnpipihi,200,0,1.0);
  pipmom_IMnpipi_woK0_wSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  pipmom_IMnpipi_woK0_wSid_n->SetYTitle("Mom(#pi^{+}) [GeV/c] ");

  pimmom_IMnpipi_woK0_wSid_n = new TH2F("pimmom_IMnpipi_woK0_wSid_n","pimmom_IMnpipi_woK0_wSid_n",nbinIMnpipi,IMnpipilow,IMnpipihi,200,0,1.0);
  pimmom_IMnpipi_woK0_wSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  pimmom_IMnpipi_woK0_wSid_n->SetYTitle("Mom(#pi^{-}) [GeV/c] ");

  nmom_IMnpipi_wK0_n = new TH2F("nmom_IMnpipi_wK0_n","nmom_IMnpipi_wK0_n", nbinIMnpipi,IMnpipilow,IMnpipihi,nbinnmom,0,1.0);
  nmom_IMnpipi_wK0_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  nmom_IMnpipi_wK0_n->SetYTitle("nmom  [GeV/c]");

  nmom_IMnpipi_wK0_wSid_n = new TH2F(Form("nmom_IMnpipi_wK0_wSid_n"),Form("nmom_IMnpipi_wK0_wSid_n"), nbinIMnpipi,IMnpipilow,IMnpipihi,nbinnmom,0,1.0);
  nmom_IMnpipi_wK0_wSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  nmom_IMnpipi_wK0_wSid_n->SetYTitle("nmom  [GeV/c]");
  
  nmom_IMnpipi_wK0_woSid_won = new TH2F("nmom_IMnpipi_wK0_woSid_won","nmom_IMnpipi_wK0_woSid_won", nbinIMnpipi,IMnpipilow,IMnpipihi,nbinnmom,0,1.0);
  nmom_IMnpipi_wK0_woSid_won->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  nmom_IMnpipi_wK0_woSid_won->SetYTitle("nmom  [GeV/c]");

  nmom_IMnpipi_wK0_woSid_n = new TH2F("nmom_IMnpipi_wK0_woSid_n","nmom_IMnpipi_wK0_woSid_n", nbinIMnpipi,IMnpipilow,IMnpipihi,nbinnmom,0,1.0);
  nmom_IMnpipi_wK0_woSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  nmom_IMnpipi_wK0_woSid_n->SetYTitle("nmom  [GeV/c]");

  nmom_MMnmiss_wK0 = new TH2F(Form("nmom_MMnmiss_wK0"),Form("nmom_MMnmiss_wK0"), 100,0.4,1.9,nbinnmom,0,1.0);
  nmom_MMnmiss_wK0->SetXTitle("Miss. Mass [GeV/c^{2}]");
  nmom_MMnmiss_wK0->SetYTitle("nmom  [GeV/c]");

  nmom_MMnmiss_wK0_woSid = new TH2F("nmom_MMnmiss_wK0_woSid","nmom_MMnmiss_wK0_woSid", 100,0.4,1.9,nbinnmom,0,1.0);
  nmom_MMnmiss_wK0_woSid->SetXTitle("Miss. Mass [GeV/c^{2}]");
  nmom_MMnmiss_wK0_woSid->SetYTitle("nmom  [GeV/c]");

  nmom_MMnmiss_wK0_woSid_won = new TH2F("nmom_MMnmiss_wK0_woSid_won","nmom_MMnmiss_wK0_woSid_won", 100,0.4,1.9,nbinnmom,0,1.0);
  nmom_MMnmiss_wK0_woSid_won->SetXTitle("Miss. Mass [GeV/c^{2}]");
  nmom_MMnmiss_wK0_woSid_won->SetYTitle("nmom  [GeV/c]");

  nmom_MMnmiss_wK0_woSidn_won = new TH2F("nmom_MMnmiss_wK0_woSidn_won","nmom_MMnmiss_wK0_woSidn_won", 100,0.4,1.9,nbinnmom,0,1.0);
  nmom_MMnmiss_wK0_woSidn_won->SetXTitle("Miss. Mass [GeV/c^{2}]");
  nmom_MMnmiss_wK0_woSidn_won->SetYTitle("nmom  [GeV/c]");

  nmom_CDHphi = new TH2F("nmom_CDHphi","nmom_CDHphi",100,-3.14,3.14,nbinnmom,0,1.0);
  nmom_CDHphi->SetXTitle("CDH phi");
  nmom_CDHphi->SetYTitle("nmom [GeV/c]");

  nmom_cosn_wK0 = new TH2F("nmom_cosn_wK0","nmom_cosn_wK0",100,-1.0,1.0,nbinnmom,0,1.0);
  nmom_cosn_wK0->SetXTitle("nCDS mom [GeV/c]");
  nmom_cosn_wK0->SetYTitle("nCDS cos#theta_{LAB}");


  //nmom_cosn_wK0_n = new TH2F("nmom_cosn_wK0_n","nmom_cosn_wK0_n",100,-1.0,1.0,nbinnmom,0,1.0);

  nmom_cosn_wSid_n = new TH2F("nmom_cosn_wSid_n","nmom_cosn_wSid_n",100,-1.0,1.0,nbinnmom,0,1.0);
  nmom_cosn_wSid_n->SetXTitle("nCDS cos#theta_{LAB}");
  nmom_cosn_wSid_n->SetYTitle("nCDS mom [GeV/c]");

  nmom_cosn_woK0_wSid_n = new TH2F("nmom_cosn_woK0_wSid_n","nmom_cosn_woK0_wSid_n",100,-1.0,1.0,nbinnmom,0,1.0);
  nmom_cosn_woK0_wSid_n->SetXTitle("nCDS cos#theta_{LAB}");
  nmom_cosn_woK0_wSid_n->SetYTitle("nCDS mom [GeV/c]");

  //nmom_cosn_woK0_woSidn = new TH2F("nmom_cosn_woK0_woSidn","nmom_cosn_woK0_woSidn",100,-1.0,1.0,nbinnmom,0,1.0);
  //nmom_cosn_woK0_woSidn->SetXTitle("nCDS cos#theta_{LAB}");
  //nmom_cosn_woK0_woSidn->SetYTitle("nCDS mom [GeV/c]");

  nmom_cosn_woK0_woSid_won = new TH2F("nmom_cosn_woK0_woSid_won","nmom_cosn_woK0_woSid_won",100,-1.0,1.0,nbinnmom,0,1.0);
  nmom_cosn_woK0_woSid_won->SetXTitle("nCDS cos#theta_{LAB}");
  nmom_cosn_woK0_woSid_won->SetYTitle("nCDS mom [GeV/c]");

  nmom_cosn_wK0_woSid_won = new TH2F("nmom_cosn_wK0_woSid_won","nmom_cosn_wK0_woSid_won",100,-1.0,1.0,nbinnmom,0,1.0);
  nmom_cosn_wK0_woSid_won->SetXTitle("nCDS cos#theta_{LAB}");
  nmom_cosn_wK0_woSid_won->SetYTitle("nCDS mom [GeV/c]");

  //nmom_cospip_woK0_woSidn = new TH2F("nmom_cospip_woK0_woSidn","nmom_cospip_woK0_woSidn",100,-1.0,1.0,nbinnmom,0,1.0);
  //nmom_cospip_woK0_woSidn->SetXTitle("#pi^{+} cos#theta_{LAB}");
  //nmom_cospip_woK0_woSidn->SetYTitle("nCDS mom [GeV/c]");

  nmom_cospip_woK0_woSid_won = new TH2F("nmom_cospip_woK0_woSid_won","nmom_cospip_woK0_woSid_won",100,-1.0,1.0,nbinnmom,0,1.0);
  nmom_cospip_woK0_woSid_won->SetXTitle("#pi^{+} cos#theta_{LAB}");
  nmom_cospip_woK0_woSid_won->SetYTitle("nCDS mom [GeV/c]");

  nmom_cospippim_woK0_woSid_won = new TH2F("nmom_cospippim_woK0_woSid_won","nmom_cospippim_woK0_woSid_won",100,-1.0,1.0,nbinnmom,0,1.0);
  nmom_cospippim_woK0_woSid_won->SetXTitle("#pi^{+}-#pi^{-} cos#theta_{LAB}");
  nmom_cospippim_woK0_woSid_won->SetYTitle("nCDS mom [GeV/c]");

  nmom_cospip_wK0_woSid_won = new TH2F("nmom_cospip_wK0_woSid_won","nmom_cospip_wK0_woSid_won",100,-1.0,1.0,nbinnmom,0,1.0);
  nmom_cospip_wK0_woSid_won->SetXTitle("#pi^{+} cos#theta_{LAB}");
  nmom_cospip_wK0_woSid_won->SetYTitle("nCDS mom [GeV/c]");

  nmom_cospippim_wK0_woSid_won = new TH2F("nmom_cospippim_wK0_woSid_won","nmom_cospippim_wK0_woSid_won",100,-1.0,1.0,nbinnmom,0,1.0);
  nmom_cospippim_wK0_woSid_won->SetXTitle("#pi^{+}-#pi^{-} cos#theta_{LAB}");
  nmom_cospippim_wK0_woSid_won->SetYTitle("nCDS mom [GeV/c]");

  //nmom_cospim_woK0_woSidn = new TH2F("nmom_cospim_woK0_woSidn","nmom_cospim_woK0_woSidn",100,-1.0,1.0,nbinnmom,0,1.0);
  //nmom_cospim_woK0_woSidn->SetXTitle("#pi^{-} cos#theta_{LAB}");
  //nmom_cospim_woK0_woSidn->SetYTitle("nCDS mom [GeV/c]");

  nmom_cospim_woK0_woSid_won = new TH2F("nmom_cospim_woK0_woSid_won","nmom_cospim_woK0_woSid_won",100,-1.0,1.0,nbinnmom,0,1.0);
  nmom_cospim_woK0_woSid_won->SetXTitle("#pi^{-} cos#theta_{LAB}");
  nmom_cospim_woK0_woSid_won->SetYTitle("nCDS mom [GeV/c]");

  nmom_coslabpim_wK0_woSid_won = new TH2F("nmom_coslabpim_wK0_woSid_won","nmom_coslabpim_wK0_woSid_won",100,-1.0,1.0,nbinnmom,0,1.0);
  nmom_coslabpim_wK0_woSid_won->SetXTitle("#pi^{-} cos#theta_{LAB}");
  nmom_coslabpim_wK0_woSid_won->SetYTitle("nCDS mom [GeV/c]");

  nmom_cospim_wK0_woSid_won = new TH2F("nmom_cospim_wK0_woSid_won","nmom_cospim_wK0_woSid_won",100,-1.0,1.0,nbinnmom,0,1.0);
  nmom_cospim_wK0_woSid_won->SetXTitle("#pi^{-} cos#theta_{LAB}");
  nmom_cospim_wK0_woSid_won->SetYTitle("nCDS mom [GeV/c]");

  //nmom_phinpip_woK0_woSidn = new TH2F("nmom_phinpip_woK0_woSidn","nmom_phinpip_woK0_woSidn",100,-1.0*TMath::Pi(),TMath::Pi(),nbinnmom,0.,1.0);
  //nmom_phinpip_woK0_woSidn->SetXTitle("#Delta#phi (nCDS-#pi^{+}) [radian]");
  //nmom_phinpip_woK0_woSidn->SetYTitle("nCDS mom [GeV/c]");

  nmom_phinpip_woK0_woSid_won = new TH2F("nmom_phinpip_woK0_woSid_won","nmom_phinpip_woK0_woSid_won",100,-1.0*TMath::Pi(),TMath::Pi(),nbinnmom,0.,1.0);
  nmom_phinpip_woK0_woSid_won->SetXTitle("#Delta#phi (nCDS-#pi^{+}) [radian]");
  nmom_phinpip_woK0_woSid_won->SetYTitle("nCDS mom [GeV/c]");

  nmom_phinpip_wK0_woSid_won = new TH2F("nmom_phinpip_wK0_woSid_won","nmom_phinpip_wK0_woSid_won",100,-1.0*TMath::Pi(),TMath::Pi(),nbinnmom,0.,1.0);
  nmom_phinpip_wK0_woSid_won->SetXTitle("#Delta#phi (nCDS-#pi^{+}) [radian]");
  nmom_phinpip_wK0_woSid_won->SetYTitle("nCDS mom [GeV/c]");

  //nmom_phinpim_woK0_woSidn = new TH2F("nmom_phinpim_woK0_woSidn","nmom_phinpim_woK0_woSidn",100,-1.0*TMath::Pi(),TMath::Pi(),nbinnmom,0.,1.0);
  //nmom_phinpim_woK0_woSidn->SetXTitle("#Delta#phi (nCDS-#pi^{-}) [radian]");
  //nmom_phinpim_woK0_woSidn->SetYTitle("nCDS mom [GeV/c]");

  nmom_phinpim_woK0_woSid_won = new TH2F("nmom_phinpim_woK0_woSid_won","nmom_phinpim_woK0_woSid_won",100,-1.0*TMath::Pi(),TMath::Pi(),nbinnmom,0.,1.0);
  nmom_phinpim_woK0_woSid_won->SetXTitle("#Delta#phi (nCDS-#pi^{-}) [radian]");
  nmom_phinpim_woK0_woSid_won->SetYTitle("nCDS mom [GeV/c]");

  nmom_phinpim_wK0_woSid_won = new TH2F("nmom_phinpim_wK0_woSid_won","nmom_phinpim_wK0_woSid_won",100,-1.0*TMath::Pi(),TMath::Pi(),nbinnmom,0.,1.0);
  nmom_phinpim_wK0_woSid_won->SetXTitle("#Delta#phi (nCDS-#pi^{-}) [radian]");
  nmom_phinpim_wK0_woSid_won->SetYTitle("nCDS mom [GeV/c]");

  //nmom_phipip_woK0_woSidn = new TH2F("nmom_phipip_woK0_woSidn","nmom_phipip_woK0_woSidn",100,-1.0*TMath::Pi(),TMath::Pi(),nbinnmom,0.,1.0);
  //nmom_phipip_woK0_woSidn->SetXTitle("#phi #pi^{+} [radian]");
  //nmom_phipip_woK0_woSidn->SetYTitle("nCDS mom [GeV/c]");

  nmom_phipip_woK0_woSid_won = new TH2F("nmom_phipip_woK0_woSid_won","nmom_phipip_woK0_woSid_won",100,-1.0*TMath::Pi(),TMath::Pi(),nbinnmom,0.,1.0);
  nmom_phipip_woK0_woSid_won->SetXTitle("phi #pi^{+} [radian]");
  nmom_phipip_woK0_woSid_won->SetYTitle("nCDS mom [GeV/c]");

  nmom_phipip_wK0_woSid_won = new TH2F("nmom_phipip_wK0_woSid_won","nmom_phipip_wK0_woSid_won",100,-1.0*TMath::Pi(),TMath::Pi(),nbinnmom,0.,1.0);
  nmom_phipip_wK0_woSid_won->SetXTitle("#phi #pi^{+} [radian]");
  nmom_phipip_wK0_woSid_won->SetYTitle("nCDS mom [GeV/c]");

  //nmom_phipim_woK0_woSidn = new TH2F("nmom_phipim_woK0_woSidn","nmom_phipim_woK0_woSidn",100,-1.0*TMath::Pi(),TMath::Pi(),nbinnmom,0.,1.0);
  //nmom_phipim_woK0_woSidn->SetXTitle("#phi #pi^{-} [radian]");
  //nmom_phipim_woK0_woSidn->SetYTitle("nCDS mom [GeV/c]");

  nmom_phipim_woK0_woSid_won = new TH2F("nmom_phipim_woK0_woSid_won","nmom_phipim_woK0_woSid_won",100,-1.0*TMath::Pi(),TMath::Pi(),nbinnmom,0.,1.0);
  nmom_phipim_woK0_woSid_won->SetXTitle("#phi #pi^{-} [radian]");
  nmom_phipim_woK0_woSid_won->SetYTitle("nCDS mom [GeV/c]");

  nmom_phipim_wK0_woSid_won = new TH2F("nmom_phipim_wK0_woSid_won","nmom_phipim_wK0_woSid_won",100,-1.0*TMath::Pi(),TMath::Pi(),nbinnmom,0.,1.0);
  nmom_phipim_wK0_woSid_won->SetXTitle("#phi #pi^{-} [radian]");
  nmom_phipim_wK0_woSid_won->SetYTitle("nCDS mom [GeV/c]");

  //nmom_phin_woK0_woSidn = new TH2F("nmom_phin_woK0_woSidn","nmom_phin_woK0_woSidn",100,-1.0*TMath::Pi(),TMath::Pi(),nbinnmom,0.,1.0);
  //nmom_phin_woK0_woSidn->SetXTitle("#phi nCDS [radian]");
  //nmom_phin_woK0_woSidn->SetYTitle("nCDS mom [GeV/c]");

  nmom_phin_woK0_woSid_won = new TH2F("nmom_phin_woK0_woSid_won","nmom_phin_woK0_woSid_won",100,-1.0*TMath::Pi(),TMath::Pi(),nbinnmom,0.,1.0);
  nmom_phin_woK0_woSid_won->SetXTitle("#phi nCDS [radian]");
  nmom_phin_woK0_woSid_won->SetYTitle("nCDS mom [GeV/c]");

  nmom_phin_wK0_woSid_won = new TH2F("nmom_phin_wK0_woSid_won","nmom_phin_wK0_woSid_won",100,-1.0*TMath::Pi(),TMath::Pi(),nbinnmom,0.,1.0);
  nmom_phin_wK0_woSid_won->SetXTitle("#phi nCDS [radian]");
  nmom_phin_wK0_woSid_won->SetYTitle("nCDS mom [GeV/c]");


  //nmom_pipmom_woK0_woSidn = new TH2F("nmom_pipmom_woK0_woSidn","nmom_pipmom_woK0_woSidn",200,0,1.0,nbinnmom,0,1.0);
  //nmom_pipmom_woK0_woSidn->SetXTitle("#pi^{+} mom [GeV/c]");
  //nmom_pipmom_woK0_woSidn->SetYTitle("nCDS mom [GeV/c]");

  nmom_pipmom_woK0_woSid_won = new TH2F("nmom_pipmom_woK0_woSid_won","nmom_pipmom_woK0_woSid_won",200,0,1.0,nbinnmom,0,1.0);
  nmom_pipmom_woK0_woSid_won->SetXTitle("#pi^{+} mom [GeV/c]");
  nmom_pipmom_woK0_woSid_won->SetYTitle("nCDS mom [GeV/c]");

  nmom_pipmom_wK0_woSid_won = new TH2F("nmom_pipmom_wK0_woSid_won","nmom_pipmom_wK0_woSid_won",200,0,1.0,nbinnmom,0,1.0);
  nmom_pipmom_wK0_woSid_won->SetXTitle("#pi^{+} mom [GeV/c]");
  nmom_pipmom_wK0_woSid_won->SetYTitle("nCDS mom [GeV/c]");

  
  //nmom_pimmom_woK0_woSidn = new TH2F("nmom_pimmom_woK0_woSidn","nmom_pimmom_woK0_woSidn",200,0,1.0,nbinnmom,0,1.0);
  //nmom_pimmom_woK0_woSidn->SetXTitle("#pi^{-} mom [GeV/c]");
  //nmom_pimmom_woK0_woSidn->SetYTitle("nCDS mom [GeV/c]");

  nmom_pimmom_woK0_woSid_won = new TH2F("nmom_pimmom_woK0_woSid_won","nmom_pimmom_woK0_woSid_won",200,0,1.0,nbinnmom,0,1.0);
  nmom_pimmom_woK0_woSid_won->SetXTitle("#pi^{-} mom [GeV/c]");
  nmom_pimmom_woK0_woSid_won->SetYTitle("nCDS mom [GeV/c]");

  nmom_pimmom_wK0_woSid_won = new TH2F("nmom_pimmom_wK0_woSid_won","nmom_pimmom_wK0_woSid_won",200,0,1.0,nbinnmom,0,1.0);
  nmom_pimmom_wK0_woSid_won->SetXTitle("#pi^{-} mom [GeV/c]");
  nmom_pimmom_wK0_woSid_won->SetYTitle("nCDS mom [GeV/c]");
  
  std::cout << __LINE__ << std::endl;

  nmom_cosn_wK0_n_forward = new TH2F("nmom_cosn_wK0_n_forward","nmom_cosn_wK0_n_forward",100,-1.0,1.0,nbinnmom,0,1.0);

  nmom_cosK0_wK0 = new TH2F("nmom_cosK0_wK0","nmom_cosK0_wK0",100,-1.0,1.0,nbinnmom,0,1.0);

  nmom_cosK0_wK0_n = new TH2F("nmom_cosK0_wK0_n","nmom_cosK0_wK0_n",100,-1.0,1.0,nbinnmom,0,1.0);

  nmom_cosK0n_wK0 = new TH2F("nmom_cosK0n_wK0","nmom_cosK0n_wK0",100,-1.0,1.0,nbinnmom,0,1.0);

  nmom_cosK0n_wK0_n = new TH2F("nmom_cosK0n_wK0_n","nmom_cosK0n_wK0_n",100,-1.0,1.0,nbinnmom,0,1.0);

  nmom_cosnmiss_wK0_n = new TH2F("nmom_cosnmiss_wK0_n","nmom_cosnmiss_wK0_n",100,-1.0,1.0,nbinnmom,0,1.0);

  nmom_cosnmiss_wSid_n = new TH2F("nmom_cosnmiss_wSid_n","nmom_cosnmiss_wSid_n",100,-1.0,1.0,nbinnmom,0,1.0);

  nmom_cosnmiss_woK0_wSid_n = new TH2F("nmom_cosnmiss_woK0_wSid_n","nmom_cosnmiss_woK0_wSid_n",100,-1.0,1.0,nbinnmom,0,1.0);

  nmom_cosnnmiss_wK0_n = new TH2F("nmom_cosnnmiss_wK0_n","nmom_cosnnmiss_wK0_n",100,-1.0,1.0,nbinnmom,0,1.0);

  K0mom_cosK0_wK0 = new TH2F("K0mom_cosK0_wK0","K0mom_cosK0_wK0",100,-1.0,1.0,100,0,1.0);

  K0mom_cosK0_wK0_n = new TH2F("K0mom_cosK0_wK0_n","K0mom_cosK0_wK0_n",100,-1.0,1.0,100,0,1.0);

  nmissmom_cosnmiss_wK0 = new TH2F("nmissmom_cosnmiss_wK0","nmissmom_cosnmiss_wK0",100,-1.0,1.0,100,0,1.5);

  nmissmom_cosnmiss_wK0_n = new TH2F("nmissmom_cosnmiss_wK0_n","nmissmom_cosnmiss_wK0_n",100,-1.0,1.0,100,0,1.5);

  nmissmom_cosK0nmiss_wK0_n = new TH2F("nmissmom_cosK0nmiss_wK0_n","nmissmom_cosK0nmiss_wK0_n",100,-1.0,1.0,100,0,1.5);

  nmom_K0mom = new TH2F("nmom_K0mom","nmom_K0mom",100,0,1.0,nbinnmom,0,1.0);
  nmom_K0mom->SetXTitle("K^{0} mom [GeV/c]");  
  nmom_K0mom->SetYTitle("n_{CDS} mom [GeV/c]");

  nmom_K0mom_n = new TH2F("nmom_K0mom_n","nmom_K0mom_n",100,0,1.0,nbinnmom,0,1.0);
  nmom_K0mom_n->SetXTitle("K^{0} mom [GeV/c]");  
  nmom_K0mom_n->SetYTitle("n_{CDS} mom [GeV/c]");
  
  nmom_K0mom_woSid_n = new TH2F("nmom_K0mom_woSid_n","nmom_K0mom_woSid_n",100,0,1.0,nbinnmom,0,1.0);
  nmom_K0mom_woSid_n->SetXTitle("K^{0} mom [GeV/c]");  
  nmom_K0mom_woSid_n->SetYTitle("n_{CDS} mom [GeV/c]");

  nmom_nmissmom_wK0 = new TH2F("nmom_nmissmom_wK0","nmom_nmissmom_wK0",100,0,1.5,nbinnmom,0,1.0);//

  nmom_nmissmom_wK0_wSid = new TH2F("nmom_nmissmom_wK0_wSid","nmom_nmissmom_wK0_wSid",100,0,1.5,nbinnmom,0,1.0);//

  nmom_nmissmom_woK0_wSid = new TH2F("nmom_nmissmom_woK0_wSid","nmom_nmissmom_woK0_wSid",100,0,1.5,nbinnmom,0,1.0);//

  nmom_nmissmom_wK0_n = new TH2F("nmom_nmissmom_wK0_n","nmom_nmissmom_wK0_n",100,0,1.5,nbinnmom,0,1.0);//

  nmom_nmissmom_wK0_wSid_n = new TH2F("nmom_nmissmom_wK0_wSid_n","nmom_nmissmom_wK0_wSid_n",100,0,1.5,nbinnmom,0,1.0);//

  nmom_nmissmom_woK0_wSid_n = new TH2F("nmom_nmissmom_woK0_wSid_n","nmom_nmissmom_woK0_wSid_n",100,0,1.5,nbinnmom,0,1.0);

  nmissmom_K0mom = new TH2F("nmissmom_K0mom","nmissmom_K0mom",100,0,1.0,100,0,1.5);

  nmissmom_K0mom_n = new TH2F("nmissmom_K0mom_n","nmissmom_K0mom_n",100,0,1.0,100,0,1.5);

  nmom_IMnpipi_wSid_n = new TH2F(Form("nmom_IMnpipi_wSid_n"),Form("nmom_IMnpipi_wSid_n"), nbinIMnpipi,IMnpipilow,IMnpipihi,nbinnmom,0,1.0);
  nmom_IMnpipi_wSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  nmom_IMnpipi_wSid_n->SetYTitle("nmom  [GeV/c]");

  //nmom_IMnpipi_woK0_wSid_n = new TH2F(Form("nmom_IMnpipi_woK0_wSid_n"),Form("nmom_IMnpipi_woK0_wSid_n"), nbinIMnpipi,IMnpipilow,IMnpipihi,nbinnmom,0,1.0);
  //nmom_IMnpipi_woK0_wSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  //nmom_IMnpipi_woK0_wSid_n->SetYTitle("nmom  [GeV/c]");
  
  nmom_IMnpipi_woK0_woSid_won = new TH2F(Form("nmom_IMnpipi_woK0_woSid_won"),Form("nmom_IMnpipi_woK0_woSid_won"), nbinIMnpipi,IMnpipilow,IMnpipihi,nbinnmom,0,1.0);
  nmom_IMnpipi_woK0_woSid_won->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  nmom_IMnpipi_woK0_woSid_won->SetYTitle("nmom  [GeV/c]");

  nmom_IMnpipi_woK0_wSid_n_Sp = new TH2F(Form("nmom_IMnpipi_woK0_wSid_n_Sp"),Form("nmom_IMnpipi_woK0_wSid_n_Sp"), nbinIMnpipi,IMnpipilow,IMnpipihi,nbinnmom,0,1.0);
  nmom_IMnpipi_woK0_wSid_n_Sp->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  nmom_IMnpipi_woK0_wSid_n_Sp->SetYTitle("nmom  [GeV/c]");

  nmom_IMnpipi_woK0_wSid_n_Sm = new TH2F(Form("nmom_IMnpipi_woK0_wSid_n_Sm"),Form("nmom_IMnpipi_woK0_wSid_n_Sm"), nbinIMnpipi,IMnpipilow,IMnpipihi,nbinnmom,0,1.0);
  nmom_IMnpipi_woK0_wSid_n_Sm->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  nmom_IMnpipi_woK0_wSid_n_Sm->SetYTitle("nmom  [GeV/c]");

  nmom_MMnmiss_wSid = new TH2F("nmom_MMnmiss_wSid","nmom_MMnmiss_wSid", nbinnmiss, nmisslow, nmisshigh,nbinnmom,0,1.0);
  nmom_MMnmiss_wSid->SetXTitle("Miss. Mass [GeV/c^{2}]");
  nmom_MMnmiss_wSid->SetYTitle("nmom  [GeV/c]");
  
  nmom_MMnmiss_wSid_fake = new TH2F("nmom_MMnmiss_wSid_fake","nmom_MMnmiss_wSid_fake", nbinnmiss, nmisslow, nmisshigh,nbinnmom,0,1.0);
  nmom_MMnmiss_wSid_fake->SetXTitle("Miss. Mass [GeV/c^{2}]");
  nmom_MMnmiss_wSid_fake->SetYTitle("nmom  [GeV/c]");
  
  nmom_MMnmiss_wSid_n = new TH2F("nmom_MMnmiss_wSid_n","nmom_MMnmiss_wSid_n", nbinnmiss, nmisslow, nmisshigh,nbinnmom,0,1.0);
  nmom_MMnmiss_wSid_n->SetXTitle("Miss. Mass [GeV/c^{2}]");
  nmom_MMnmiss_wSid_n->SetYTitle("nmom  [GeV/c]");
  
  nmom_MMnmiss_woK0_wSid = new TH2F(Form("nmom_MMnmiss_woK0_wSid"),Form("nmom_MMnmiss_woK0_wSid"), nbinnmiss, nmisslow, nmisshigh,nbinnmom,0,1.0);
  nmom_MMnmiss_woK0_wSid->SetXTitle("Miss. Mass [GeV/c^{2}]");
  nmom_MMnmiss_woK0_wSid->SetYTitle("nmom  [GeV/c]");

  nmom_MMnmiss_woK0_woSid = new TH2F("nmom_MMnmiss_woK0_woSid","nmom_MMnmiss_woK0_woSid", nbinnmiss, nmisslow, nmisshigh,nbinnmom,0,1.0);
  nmom_MMnmiss_woK0_woSid->SetXTitle("Miss. Mass [GeV/c^{2}]");
  nmom_MMnmiss_woK0_woSid->SetYTitle("nmom  [GeV/c]");

  //nmom_MMnmiss_woK0_woSidn = new TH2F("nmom_MMnmiss_woK0_woSidn","nmom_MMnmiss_woK0_woSidn", nbinnmiss, nmisslow, nmisshigh,nbinnmom,0,1.0);
  //nmom_MMnmiss_woK0_woSidn->SetXTitle("Miss. Mass [GeV/c^{2}]");
  //nmom_MMnmiss_woK0_woSidn->SetYTitle("nmom  [GeV/c]");

  nmom_MMnmiss_woK0_woSid_won = new TH2F("nmom_MMnmiss_woK0_woSid_won","nmom_MMnmiss_woK0_woSid_won", nbinnmiss, nmisslow, nmisshigh,nbinnmom,0,1.0);
  nmom_MMnmiss_woK0_woSid_won->SetXTitle("Miss. Mass [GeV/c^{2}]");
  nmom_MMnmiss_woK0_woSid_won->SetYTitle("nmom  [GeV/c]");

  nmom_cosnlab_K0_n = new TH2F(Form("nmom_cosnlab_K0_n"),Form("nmom_cosnlab_K0_n"), 200,-1,1,nbinnmom,0,1.0);
  nmom_cosnlab_K0_n->SetXTitle("cos_n (lab.)");
  nmom_cosnlab_K0_n->SetYTitle("nmom  [GeV/c]");

  nmom_IMpippim = new TH2F(Form("nmom_IMpippim"),Form("nmom_IMpippim"), nbinpippim,0.,0.9,nbinnmom,0,1.0);
  nmom_IMpippim->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  nmom_IMpippim->SetYTitle("nmom  [GeV/c]");

  nmom_MK0bar2 = new TH2F(Form("nmom_MK0bar2"),Form("nmom_MK0bar2"), 150,-0.5,1,nbinnmom,0,1.0);
  nmom_MK0bar2->SetXTitle("(M(#bar{K^{0}})^{2} [(GeV/c^{2})^{2}]");
  nmom_MK0bar2->SetYTitle("nmom  [GeV/c]");

  nmom_IMpippim_n = new TH2F(Form("nmom_IMpippim_n"),Form("nmom_IMpippim_n"), nbinpippim,0.,0.9,nbinnmom,0,1.0);
  nmom_IMpippim_n->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  nmom_IMpippim_n->SetYTitle("nmom  [GeV/c]");

  nmom_IMpippim_wSid_n = new TH2F("nmom_IMpippim_wSid_n","nmom_IMpippim_wSid_n", nbinpippim,0.,0.9,nbinnmom,0,1.0);
  nmom_IMpippim_wSid_n->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  nmom_IMpippim_wSid_n->SetYTitle("nmom  [GeV/c]");
  
  nmom_IMpippim_woK0_woSid_won = new TH2F("nmom_IMpippim_woK0_woSid_won","nmom_IMpippim_woK0_woSid_won", nbinpippim,0.,0.9,nbinnmom,0,1.0);
  nmom_IMpippim_woK0_woSid_won->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  nmom_IMpippim_woK0_woSid_won->SetYTitle("nmom  [GeV/c]");
  
  nmom_IMpippim_wK0_woSid_won = new TH2F("nmom_IMpippim_wK0_woSid_won","nmom_IMpippim_wK0_woSid_won", nbinpippim,0.,0.9,nbinnmom,0,1.0);
  nmom_IMpippim_wK0_woSid_won->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  nmom_IMpippim_wK0_woSid_won->SetYTitle("nmom  [GeV/c]");

  mnmom_IMpippim_n = new TH2F(Form("mnmom_IMpippim_n"),Form("mnmom_IMpippim_n"), nbinpippim,0.,0.9,100,0,2.0);
  mnmom_IMpippim_n->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  mnmom_IMpippim_n->SetYTitle("miss. nmom  [GeV/c]");

  mnmom_IMpippim_wSid_n = new TH2F(Form("mnmom_IMpippim_wSid_n"),Form("mnmom_IMpippim_wSid_n"), nbinpippim,0.,0.9,100,0,2.0);
  mnmom_IMpippim_wSid_n->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  mnmom_IMpippim_wSid_n->SetYTitle("miss. nmom  [GeV/c]");

  q_IMpippim_n = new TH2F("q_IMpippim_n","q_IMpippim_n",nbinpippim,0.,0.9, 100,0,1.5);
  q_IMpippim_n->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMpippim_n->SetYTitle("Mom. Transfer [GeV/c]");

  q_IMpippim_wSid_n = new TH2F("q_IMpippim_wSid_n","q_IMpippim_wSid_n",nbinpippim,0.,0.9, 100,0,1.5);
  q_IMpippim_wSid_n->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMpippim_wSid_n->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_IMpippim_woK0_woSid_won = new TH2F("q_IMpippim_woK0_woSid_won","q_IMpippim_woK0_woSid_won",nbinpippim,0.,0.9, 100,0,1.5);
  q_IMpippim_woK0_woSid_won->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMpippim_woK0_woSid_won->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_IMpippim_wK0_woSid_won = new TH2F("q_IMpippim_wK0_woSid_won","q_IMpippim_wK0_woSid_won",nbinpippim,0.,0.9, 100,0,1.5);
  q_IMpippim_wK0_woSid_won->SetXTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMpippim_wK0_woSid_won->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_IMnpip_n = new TH2F("q_IMnpip_n","q_IMnpip_n",nbinIMnpi, 1.0, 2.0, 100,0,1.5);
  q_IMnpip_n->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  q_IMnpip_n->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_IMnpip_n_wSid = new TH2F("q_IMnpip_n_wSid","q_IMnpip_n_wSid",nbinIMnpi, 1.0, 2.0, 100,0,1.5);
  q_IMnpip_n_wSid->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  q_IMnpip_n_wSid->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_IMnpip_woK0_woSid_won = new TH2F("q_IMnpip_woK0_woSid_won","q_IMnpip_woK0_woSid_won",nbinIMnpi, 1.0, 2.0, 100,0,1.5);
  q_IMnpip_woK0_woSid_won->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  q_IMnpip_woK0_woSid_won->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_IMnpip_wK0_woSid_won = new TH2F("q_IMnpip_wK0_woSid_won","q_IMnpip_wK0_woSid_won",nbinIMnpi, 1.0, 2.0, 100,0,1.5);
  q_IMnpip_wK0_woSid_won->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  q_IMnpip_wK0_woSid_won->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_IMnpim_n = new TH2F("q_IMnpim_n","q_IMnpim_n",nbinIMnpi, 1.0, 2.0, 100,0,1.5);
  q_IMnpim_n->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  q_IMnpim_n->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_IMnpim_n_wSid = new TH2F("q_IMnpim_n_wSid","q_IMnpim_n_wSid",nbinIMnpi, 1.0, 2.0, 100,0,1.5);
  q_IMnpim_n_wSid->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  q_IMnpim_n_wSid->SetYTitle("Mom. Transfer [GeV/c]");
   
 
  q_IMnpim_woK0_woSid_won = new TH2F("q_IMnpim_woK0_woSid_won","q_IMnpim_woK0_woSid_won",nbinIMnpi, 1.0, 2.0, 100,0,1.5);
  q_IMnpim_woK0_woSid_won->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  q_IMnpim_woK0_woSid_won->SetYTitle("Mom. Transfer [GeV/c]");
  
  q_IMnpim_wK0_woSid_won = new TH2F("q_IMnpim_wK0_woSid_won","q_IMnpim_wK0_woSid_won",nbinIMnpi, 1.0, 2.0, 100,0,1.5);
  q_IMnpim_wK0_woSid_won->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  q_IMnpim_wK0_woSid_won->SetYTitle("Mom. Transfer [GeV/c]");
  
  //q_MMnmiss = new TH2F("q_MMnmiss","q_MMnmiss",nbinnmiss, nmisslow, nmisshigh,100,0.,1.5);
  //q_MMnmiss->SetXTitle("Miss. Mass. [GeV/c^{2}]");
  //q_MMnmiss->SetYTitle("Mom. Transfer [GeV/c]");

  //q_MMnmiss_wSid = new TH2F("q_MMnmiss_wSid","q_MMnmiss_wSid",nbinnmiss, nmisslow, nmisshigh,100,0.,1.5);
  //q_MMnmiss_wSid->SetXTitle("Miss. Mass. [GeV/c^{2}]");
  //q_MMnmiss_wSid->SetYTitle("Mom. Transfer [GeV/c]");
  
  //q_MMnmiss_n_wSid = new TH2F("q_MMnmiss_n_wSid","q_MMnmiss_n_wSid",nbinnmiss, nmisslow, nmisshigh,100,0.,1.5);
  //q_MMnmiss_n_wSid->SetXTitle("Miss. Mass. [GeV/c^{2}]");
  //q_MMnmiss_n_wSid->SetYTitle("Mom. Transfer [GeV/c]");
  
  //q_MMnmiss_woK0_woSid_won = new TH2F("q_MMnmiss_woK0_woSid_won","q_MMnmiss_woK0_woSid_won",nbinnmiss, nmisslow, nmisshigh,100,0.,1.5);
  //q_MMnmiss_woK0_woSid_won->SetXTitle("Miss. Mass. [GeV/c^{2}]");
  //q_MMnmiss_woK0_woSid_won->SetYTitle("Mom. Transfer [GeV/c]");
  
  //q_MMnmiss_wK0_woSid_won = new TH2F("q_MMnmiss_wK0_woSid_won","q_MMnmiss_wK0_woSid_won",nbinnmiss, nmisslow, nmisshigh,100,0.,1.5);
  //q_MMnmiss_wK0_woSid_won->SetXTitle("Miss. Mass. [GeV/c^{2}]");
  //q_MMnmiss_wK0_woSid_won->SetYTitle("Mom. Transfer [GeV/c]");
  
  //q_nmom_n = new TH2F("q_nmom_n","q_nmom_n",nbinnmom,0,1,100,0.,1.5);
  //q_nmom_n->SetXTitle("nmom  [GeV/c]");
  //q_nmom_n->SetYTitle("Mom. Transfer [GeV/c]");

  //q_nmom_n_wSid = new TH2F("q_nmom_n_wSid","q_nmom_n_wSid",nbinnmom,0,1,100,0.,1.5);
  //q_nmom_n_wSid->SetXTitle("nmom  [GeV/c]");
  //q_nmom_n_wSid->SetYTitle("Mom. Transfer [GeV/c]");
  
  IMpippim_IMnpipi_n = new TH2F("IMpippim_IMnpipi_n","IMpippim_IMnpipi_n",100,1,2,nbinpippim,0.,0.9);
  IMpippim_IMnpipi_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  IMpippim_IMnpipi_n->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");

  IMpippim_IMnpipi_n_wSid = new TH2F("IMpippim_IMnpipi_n_wSid","IMpippim_IMnpipi_n_wSid",100,1,2,nbinpippim,0.,0.9);
  IMpippim_IMnpipi_n_wSid->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  IMpippim_IMnpipi_n_wSid->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  
  IMpippim_IMnpipi_n_wSid_woSp = new TH2F("IMpippim_IMnpipi_n_wSid_woSp","IMpippim_IMnpipi_n_wSid_woSp",100,1,2,nbinpippim,0.,0.9);
  IMpippim_IMnpipi_n_wSid_woSp->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  IMpippim_IMnpipi_n_wSid_woSp->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");

  IMpippim_IMnpipi_n_wSid_woSm = new TH2F("IMpippim_IMnpipi_n_wSid_woSm","IMpippim_IMnpipi_n_wSid_woSm",100,1,2,nbinpippim,0.,0.9);
  IMpippim_IMnpipi_n_wSid_woSm->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  IMpippim_IMnpipi_n_wSid_woSm->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");


  IMpippim_IMnpip_vici = new TH2F("IMpippim_IMnpip_vici","IMpippim_IMnpip_vici",nbinIMnpi,1,2,nbinpippim,0.,0.9);
  IMpippim_IMnpip_vici->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMpippim_IMnpip_vici->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  
  IMpippim_IMnpip_vici_woSm = new TH2F("IMpippim_IMnpip_vici_woSm","IMpippim_IMnpip_vici_woSm",nbinIMnpi,1,2,nbinpippim,0.,0.9);
  IMpippim_IMnpip_vici_woSm->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMpippim_IMnpip_vici_woSm->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");

  IMpippim_IMnpip_n = new TH2F("IMpippim_IMnpip_n","IMpippim_IMnpip_n",nbinIMnpi,1,2,nbinpippim,0.,0.9);
  IMpippim_IMnpip_n->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMpippim_IMnpip_n->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  
  /*
  for(unsigned int ibin=0;ibin<nbintemplate;ibin++){
    float binlow=1.0+(float)ibin*1./nbintemplate;
    float binhigh=1.0+(float)(ibin+1.0)*1./nbintemplate;
    IMpippim_IMnpip_n_bin[ibin] = new TH2F(Form("IMpippim_IMnpip_n_bin%d",ibin),Form("IMpippim_IMnpip_n %0.2f-%0.2f",binlow,binhigh)
        ,nbinIMnpi,1,2,nbinpippim,0.,0.9);
    IMpippim_IMnpip_n_bin[ibin]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
    IMpippim_IMnpip_n_bin[ibin]->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  }*/
  
  IMpippim_IMnpip_n_woSm = new TH2F("IMpippim_IMnpip_n_woSm","IMpippim_IMnpip_n_woSm",nbinIMnpi,1,2,nbinpippim,0.,0.9);
  IMpippim_IMnpip_n_woSm->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMpippim_IMnpip_n_woSm->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  
  /*
  for(unsigned int ibin=0;ibin<nbintemplate;ibin++){
    float binlow=1.0+(float)ibin*1./nbintemplate;
    float binhigh=1.0+(float)(ibin+1.0)*1./nbintemplate;
    IMpippim_IMnpip_n_woSm_bin[ibin] = new TH2F(Form("IMpippim_IMnpip_n_woSm_bin%d",ibin),Form("IMpippim_IMnpip_n_woSm %0.2f-%0.2f",binlow,binhigh)
        ,nbinIMnpi,1,2,nbinpippim,0.,0.9);
    IMpippim_IMnpip_n_woSm_bin[ibin]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
    IMpippim_IMnpip_n_woSm_bin[ibin]->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  }*/
  
  //select SigmaP by diagonal cut
  IMpippim_IMnpip_n_woSmdia = new TH2F("IMpippim_IMnpip_n_woSmdia","IMpippim_IMnpip_n_woSmdia",nbinIMnpi,1,2,nbinpippim,0.,0.9);
  IMpippim_IMnpip_n_woSmdia->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMpippim_IMnpip_n_woSmdia->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  
  /*
  for(unsigned int ibin=0;ibin<nbintemplate;ibin++){
    float binlow=1.0+(float)ibin*1./nbintemplate;
    float binhigh=1.0+(float)(ibin+1.0)*1./nbintemplate;
    IMpippim_IMnpip_n_woSmdia_bin[ibin] = new TH2F(Form("IMpippim_IMnpip_n_woSmdia_bin%d",ibin),Form("IMpippim_IMnpip_n_woSmdia %0.2f-%0.2f",binlow,binhigh)
        ,nbinIMnpi,1,2,nbinpippim,0.,0.9);
    IMpippim_IMnpip_n_woSmdia_bin[ibin]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
    IMpippim_IMnpip_n_woSmdia_bin[ibin]->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  }*/

  IMpippim_IMnpip_wSid_n_Sp = new TH2F("IMpippim_IMnpip_wSid_n_Sp","IMpippim_IMnpip_wSid_n_Sp",nbinIMnpi,1,2,nbinpippim,0.,0.9);
  IMpippim_IMnpip_wSid_n_Sp->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMpippim_IMnpip_wSid_n_Sp->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
 
  for(unsigned int ibin=0;ibin<nbintemplate;ibin++){
    float binlow=1.0+(float)ibin*1./nbintemplate;
    float binhigh=1.0+(float)(ibin+1.0)*1./nbintemplate;
    IMpippim_IMnpip_wSid_n_Sp_bin[ibin] = new TH2F(Form("IMpippim_IMnpip_wSid_n_Sp_bin%d",ibin),
                                                   Form("IMpippim_IMnpip_wSid_n_Sp %0.2f-%0.2f",binlow,binhigh),
                                                   nbinIMnpip_wbin,IMnpip_wbinlow,IMnpip_wbinhigh,nbinIMpippim_wbin,IMpippim_wbinlow,IMpippim_wbinhigh);
    IMpippim_IMnpip_wSid_n_Sp_bin[ibin]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
    IMpippim_IMnpip_wSid_n_Sp_bin[ibin]->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  }
  
  for(unsigned int ibin=0;ibin<nwbin;ibin++){
    IMpippim_IMnpip_wSid_n_Sp_wbin[ibin] = new TH2F(Form("IMpippim_IMnpip_wSid_n_Sp_wbin%d",ibin),Form("IMpippim_IMnpip_wSid_n_Sp %0.2f-%0.2f",wbinlow[ibin],wbinhigh[ibin])
        ,nbinIMnpip_wbin,IMnpip_wbinlow,IMnpip_wbinhigh,nbinIMpippim_wbin,IMpippim_wbinlow,IMpippim_wbinhigh);
    IMpippim_IMnpip_wSid_n_Sp_wbin[ibin]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
    IMpippim_IMnpip_wSid_n_Sp_wbin[ibin]->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  }
 
  IMpippim_IMnpip_wSid_n_woSm = new TH2F("IMpippim_IMnpip_wSid_n_woSm","IMpippim_IMnpip_wSid_n_woSm",nbinIMnpi,1,2,nbinpippim,0.,0.9);
  IMpippim_IMnpip_wSid_n_woSm->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMpippim_IMnpip_wSid_n_woSm->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  
  /*
  for(unsigned int ibin=0;ibin<nbintemplate;ibin++){
    float binlow=1.0+(float)ibin*1./nbintemplate;
    float binhigh=1.0+(float)(ibin+1.0)*1./nbintemplate;
    IMpippim_IMnpip_wSid_n_woSm_bin[ibin] = new TH2F(Form("IMpippim_IMnpip_wSid_n_woSm_bin%d",ibin),Form("IMpippim_IMnpip_wSid_n_woSm %0.2f-%0.2f",binlow,binhigh)
        ,nbinIMnpi,1,2,nbinpippim,0.,0.9);
    IMpippim_IMnpip_wSid_n_woSm_bin[ibin]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
    IMpippim_IMnpip_wSid_n_woSm_bin[ibin]->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  }*/
  
  IMpippim_IMnpip_wK0_n = new TH2F("IMpippim_IMnpip_wK0_n","IMpippim_IMnpip_wK0_n",nbinIMnpi,1,2,nbinpippim,0.,0.9);
  IMpippim_IMnpip_wK0_n->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMpippim_IMnpip_wK0_n->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
 
  for(unsigned int ibin=0;ibin<nbintemplate;ibin++){
    float binlow=1.0+(float)ibin*1./nbintemplate;
    float binhigh=1.0+(float)(ibin+1.0)*1./nbintemplate;
    IMpippim_IMnpip_wK0_n_bin[ibin] = new TH2F(Form("IMpippim_IMnpip_wK0_n_bin%d",ibin),
        Form("IMpippim_IMnpip_wK0_n %0.2f-%0.2f",binlow,binhigh),
        nbinIMnpip_wbin,IMnpip_wbinlow,IMnpip_wbinhigh,nbinIMpippim_wbin,IMpippim_wbinlow,IMpippim_wbinhigh);
    IMpippim_IMnpip_wK0_n_bin[ibin]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
    IMpippim_IMnpip_wK0_n_bin[ibin]->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  }
  
  for(unsigned int ibin=0;ibin<nwbin;ibin++){
    IMpippim_IMnpip_wK0_n_wbin[ibin] = new TH2F(Form("IMpippim_IMnpip_wK0_n_wbin%d",ibin),Form("IMpippim_IMnpip_wK0_n %0.2f-%0.2f",wbinlow[ibin],wbinhigh[ibin])
        ,nbinIMnpip_wbin,IMnpip_wbinlow,IMnpip_wbinhigh,nbinIMpippim_wbin,IMpippim_wbinlow,IMpippim_wbinhigh);
    IMpippim_IMnpip_wK0_n_wbin[ibin]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
    IMpippim_IMnpip_wK0_n_wbin[ibin]->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  }
  
  IMpippim_IMnpip_wK0_n_woSmdia = new TH2F("IMpippim_IMnpip_wK0_n_woSmdia","IMpippim_IMnpip_wK0_n_woSmdia",nbinIMnpi,1,2,nbinpippim,0.,0.9);
  IMpippim_IMnpip_wK0_n_woSmdia->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMpippim_IMnpip_wK0_n_woSmdia->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  
  /*
  for(unsigned int ibin=0;ibin<nbintemplate;ibin++){
    float binlow=1.0+(float)ibin*1./nbintemplate;
    float binhigh=1.0+(float)(ibin+1.0)*1./nbintemplate;
    IMpippim_IMnpip_wK0_n_woSmdia_bin[ibin] = new TH2F(Form("IMpippim_IMnpip_wK0_n_woSmdia_bin%d",ibin),Form("IMpippim_IMnpip_wK0_n_woSmdia %0.2f-%0.2f",binlow,binhigh)
        ,nbinIMnpi,1,2,nbinpippim,0.,0.9);
    IMpippim_IMnpip_wK0_n_woSmdia_bin[ibin]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
    IMpippim_IMnpip_wK0_n_woSmdia_bin[ibin]->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  }*/
  
  IMpippim_IMnpip_wK0orwSid_n = new TH2F("IMpippim_IMnpip_wK0orwSid_n","IMpippim_IMnpip_wK0orwSid_n",nbinIMnpi,1,2,nbinpippim,0.,0.9);
  IMpippim_IMnpip_wK0orwSid_n->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMpippim_IMnpip_wK0orwSid_n->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
 
  for(unsigned int ibin=0;ibin<nbintemplate;ibin++){
    float binlow=1.0+(float)ibin*1./nbintemplate;
    float binhigh=1.0+(float)(ibin+1.0)*1./nbintemplate;
    IMpippim_IMnpip_wK0orwSid_n_bin[ibin] = new TH2F(Form("IMpippim_IMnpip_wK0orwSid_n_bin%d",ibin),Form("IMpippim_IMnpip_wK0orwSid_n %0.2f-%0.2f",binlow,binhigh)
        ,nbinIMnpi/4,1,2,nbinpippim/5,0.,0.9);
    IMpippim_IMnpip_wK0orwSid_n_bin[ibin]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
    IMpippim_IMnpip_wK0orwSid_n_bin[ibin]->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  }
  
  for(unsigned int iwbin=0;iwbin<nwbin;iwbin++){
    IMpippim_IMnpip_wK0orwSid_n_wbin[iwbin] = new TH2F(Form("IMpippim_IMnpip_wK0orwSid_n_wbin%d",iwbin),
                                                       Form("IMpippim_IMnpip_wK0orwSid_n %0.2f-%0.2f",wbinlow[iwbin],wbinhigh[iwbin]),
                                                       nbinIMnpip_wbin,IMnpip_wbinlow,IMnpip_wbinhigh,nbinIMpippim_wbin,IMpippim_wbinlow,IMpippim_wbinhigh);
    IMpippim_IMnpip_wK0orwSid_n_wbin[iwbin]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
    IMpippim_IMnpip_wK0orwSid_n_wbin[iwbin]->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  }
 
  IMpippim_IMnpip_wK0orwSid_n_woSmdia = new TH2F("IMpippim_IMnpip_wK0orwSid_n_woSmdia","IMpippim_IMnpip_wK0orwSid_n_woSmdia",nbinIMnpi,1,2,nbinpippim,0.,0.9);
  IMpippim_IMnpip_wK0orwSid_n_woSmdia->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMpippim_IMnpip_wK0orwSid_n_woSmdia->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  
  /*
  for(unsigned int ibin=0;ibin<nbintemplate;ibin++){
    float binlow=1.0+(float)ibin*1./nbintemplate;
    float binhigh=1.0+(float)(ibin+1.0)*1./nbintemplate;
    IMpippim_IMnpip_wK0orwSid_n_woSmdia_bin[ibin] = new TH2F(Form("IMpippim_IMnpip_wK0orwSid_n_woSmdia_bin%d",ibin),Form("IMpippim_IMnpip_wK0orwSid_n_woSmdia %0.2f-%0.2f",binlow,binhigh)
        ,nbinIMnpi,1,2,nbinpippim,0.,0.9);
    IMpippim_IMnpip_wK0orwSid_n_woSmdia_bin[ibin]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
    IMpippim_IMnpip_wK0orwSid_n_woSmdia_bin[ibin]->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  }*/

  IMpippim_IMnpip_woK0_wSid_n= new TH2F("IMpippim_IMnpip_woK0_wSid_n","IMpippim_IMnpip_woK0_wSid_n",nbinIMnpi,1,2,nbinpippim,0.,0.9);
  IMpippim_IMnpip_woK0_wSid_n->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMpippim_IMnpip_woK0_wSid_n->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  
  for(unsigned int iwbin=0;iwbin<nwbin;iwbin++){
    IMpippim_IMnpip_woK0_wSid_n_woSp_wbin[iwbin]= new TH2F(Form("IMpippim_IMnpip_woK0_wSid_n_woSp_wbin%d",iwbin)
                                                          ,Form("IMpippim_IMnpip_woK0_wSid_n_woSp %0.2f-%0.2f",wbinlow[iwbin],wbinhigh[iwbin])
                                                          ,nbinIMnpip_wbin,IMnpip_wbinlow,IMnpip_wbinhigh,nbinIMpippim_wbin,IMpippim_wbinlow,IMpippim_wbinhigh);
    IMpippim_IMnpip_woK0_wSid_n_woSp_wbin[iwbin]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
    IMpippim_IMnpip_woK0_wSid_n_woSp_wbin[iwbin]->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
    
    IMpippim_IMnpip_woK0_wSid_n_woSm_wbin[iwbin]= new TH2F(Form("IMpippim_IMnpip_woK0_wSid_n_woSm_wbin%d",iwbin)
                                                          ,Form("IMpippim_IMnpip_woK0_wSid_n_woSm %0.2f-%0.2f",wbinlow[iwbin],wbinhigh[iwbin])
                                                          ,nbinIMnpip_wbin,IMnpip_wbinlow,IMnpip_wbinhigh,nbinIMpippim_wbin,IMpippim_wbinlow,IMpippim_wbinhigh);
    IMpippim_IMnpip_woK0_wSid_n_woSm_wbin[iwbin]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
    IMpippim_IMnpip_woK0_wSid_n_woSm_wbin[iwbin]->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  }
  
  
  IMpippim_IMnpip_n_cross = new TH2F("IMpippim_IMnpip_n_cross","IMpippim_IMnpip_n_cross",nbinIMnpi,1,2,nbinpippim,0.,0.9);
  IMpippim_IMnpip_n_cross->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMpippim_IMnpip_n_cross->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  
  IMpippim_IMnpim_vici = new TH2F("IMpippim_IMnpim_vici","IMpippim_IMnpim_vici",nbinIMnpi,1,2,nbinpippim,0.,0.9);
  IMpippim_IMnpim_vici->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  IMpippim_IMnpim_vici->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  
  IMpippim_IMnpim_vici_woSp = new TH2F("IMpippim_IMnpim_vici_woSp","IMpippim_IMnpim_vici_woSp",nbinIMnpi,1,2,nbinpippim,0.,0.9);
  IMpippim_IMnpim_vici_woSp->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  IMpippim_IMnpim_vici_woSp->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");

  IMpippim_IMnpim_n = new TH2F("IMpippim_IMnpim_n","IMpippim_IMnpim_n",nbinIMnpi,1,2,nbinpippim,0.,0.9);
  IMpippim_IMnpim_n->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  IMpippim_IMnpim_n->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  
  /*
  for(unsigned int ibin=0;ibin<nbintemplate;ibin++){
    float binlow=1.0+(float)ibin*1./nbintemplate;
    float binhigh=1.0+(float)(ibin+1.0)*1./nbintemplate;
    IMpippim_IMnpim_n_bin[ibin] = new TH2F(Form("IMpippim_IMnpim_n_bin%d",ibin),Form("IMpippim_IMnpim_n %0.2f-%0.2f",binlow,binhigh)
        ,nbinIMnpi,1,2,nbinpippim,0.,0.9);
    IMpippim_IMnpim_n_bin[ibin]->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
    IMpippim_IMnpim_n_bin[ibin]->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  }*/
  
  IMpippim_IMnpim_n_woSp = new TH2F("IMpippim_IMnpim_n_woSp","IMpippim_IMnpim_n_woSp",nbinIMnpi,1,2,nbinpippim,0.,0.9);
  IMpippim_IMnpim_n_woSp->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  IMpippim_IMnpim_n_woSp->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  
  /*
  for(unsigned int ibin=0;ibin<nbintemplate;ibin++){
    float binlow=1.0+(float)ibin*1./nbintemplate;
    float binhigh=1.0+(float)(ibin+1.0)*1./nbintemplate;
    IMpippim_IMnpim_n_woSp_bin[ibin] = new TH2F(Form("IMpippim_IMnpim_n_woSp_bin%d",ibin),Form("IMpippim_IMnpim_n_woSp %0.2f-%0.2f",binlow,binhigh)
        ,nbinIMnpi,1,2,nbinpippim,0.,0.9);
    IMpippim_IMnpim_n_woSp_bin[ibin]->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
    IMpippim_IMnpim_n_woSp_bin[ibin]->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  }*/
  
  IMpippim_IMnpim_n_woSpdia = new TH2F("IMpippim_IMnpim_n_woSpdia","IMpippim_IMnpim_n_woSpdia",nbinIMnpi,1,2,nbinpippim,0.,0.9);
  IMpippim_IMnpim_n_woSpdia->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  IMpippim_IMnpim_n_woSpdia->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  
  /*
  for(unsigned int ibin=0;ibin<nbintemplate;ibin++){
    float binlow=1.0+(float)ibin*1./nbintemplate;
    float binhigh=1.0+(float)(ibin+1.0)*1./nbintemplate;
    IMpippim_IMnpim_n_woSpdia_bin[ibin] = new TH2F(Form("IMpippim_IMnpim_n_woSpdia_bin%d",ibin),Form("IMpippim_IMnpim_n_woSpdia %0.2f-%0.2f",binlow,binhigh)
        ,nbinIMnpi,1,2,nbinpippim,0.,0.9);
    IMpippim_IMnpim_n_woSpdia_bin[ibin]->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
    IMpippim_IMnpim_n_woSpdia_bin[ibin]->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  }*/
  
  IMpippim_IMnpim_wSid_n_Sm = new TH2F("IMpippim_IMnpim_wSid_n_Sm","IMpippim_IMnpim_wSid_n_Sm",nbinIMnpi,1,2,nbinpippim,0.,0.9);
  IMpippim_IMnpim_wSid_n_Sm->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  IMpippim_IMnpim_wSid_n_Sm->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  
  for(unsigned int ibin=0;ibin<nbintemplate;ibin++){
    float binlow=1.0+(float)ibin*1./nbintemplate;
    float binhigh=1.0+(float)(ibin+1.0)*1./nbintemplate;
    IMpippim_IMnpim_wSid_n_Sm_bin[ibin] = new TH2F(Form("IMpippim_IMnpim_wSid_n_Sm_bin%d",ibin),Form("IMpippim_IMnpim_wSid_n_Sm %0.2f-%0.2f",binlow,binhigh)
        ,nbinIMnpim_wbin,IMnpim_wbinlow,IMnpim_wbinhigh,nbinIMpippim_wbin,IMpippim_wbinlow,IMpippim_wbinhigh);
    IMpippim_IMnpim_wSid_n_Sm_bin[ibin]->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
    IMpippim_IMnpim_wSid_n_Sm_bin[ibin]->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  }

  for(unsigned int ibin=0;ibin<nwbin;ibin++){
    IMpippim_IMnpim_wSid_n_Sm_wbin[ibin] = new TH2F(Form("IMpippim_IMnpim_wSid_n_Sm_wbin%d",ibin),Form("IMpippim_IMnpim_wSid_n_Sm %0.2f-%0.2f",wbinlow[ibin],wbinhigh[ibin])
        ,nbinIMnpim_wbin,IMnpim_wbinlow,IMnpim_wbinhigh,nbinIMpippim_wbin,IMpippim_wbinlow,IMpippim_wbinhigh);
    IMpippim_IMnpim_wSid_n_Sm_wbin[ibin]->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
    IMpippim_IMnpim_wSid_n_Sm_wbin[ibin]->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  }

  IMpippim_IMnpim_wSid_n_woSp = new TH2F("IMpippim_IMnpim_wSid_n_woSp","IMpippim_IMnpim_wSid_n_woSp",nbinIMnpi,1,2,nbinpippim,0.,0.9);
  IMpippim_IMnpim_wSid_n_woSp->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  IMpippim_IMnpim_wSid_n_woSp->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  
  /*
  for(unsigned int ibin=0;ibin<nbintemplate;ibin++){
    float binlow=1.0+(float)ibin*1./nbintemplate;
    float binhigh=1.0+(float)(ibin+1.0)*1./nbintemplate;
    IMpippim_IMnpim_wSid_n_woSp_bin[ibin] = new TH2F(Form("IMpippim_IMnpim_wSid_n_woSp_bin%d",ibin),Form("IMpippim_IMnpim_wSid_n_woSp %0.2f-%0.2f",binlow,binhigh)
        ,nbinIMnpi,1,2,nbinpippim,0.,0.9);
    IMpippim_IMnpim_wSid_n_woSp_bin[ibin]->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
    IMpippim_IMnpim_wSid_n_woSp_bin[ibin]->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  }*/

  IMpippim_IMnpim_wK0_n = new TH2F("IMpippim_IMnpim_wK0_n","IMpippim_IMnpim_wK0_n",nbinIMnpi,1,2,nbinpippim,0.,0.9);
  IMpippim_IMnpim_wK0_n->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  IMpippim_IMnpim_wK0_n->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  
  for(unsigned int ibin=0;ibin<nbintemplate;ibin++){
    float binlow=1.0+(float)ibin*1./nbintemplate;
    float binhigh=1.0+(float)(ibin+1.0)*1./nbintemplate;
    IMpippim_IMnpim_wK0_n_bin[ibin] = new TH2F(Form("IMpippim_IMnpim_wK0_n_bin%d",ibin),Form("IMpippim_IMnpim_wK0_n %0.2f-%0.2f",binlow,binhigh)
        ,nbinIMnpi/4,1,2,nbinpippim/5,0.,0.9);
    IMpippim_IMnpim_wK0_n_bin[ibin]->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
    IMpippim_IMnpim_wK0_n_bin[ibin]->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  }
  
  for(unsigned int ibin=0;ibin<nwbin;ibin++){
    IMpippim_IMnpim_wK0_n_wbin[ibin] = new TH2F(Form("IMpippim_IMnpim_wK0_n_wbin%d",ibin),Form("IMpippim_IMnpim_wK0_n %0.2f-%0.2f",wbinlow[ibin],wbinhigh[ibin])
        ,nbinIMnpim_wbin,IMnpim_wbinlow,IMnpim_wbinhigh,nbinIMpippim_wbin,IMpippim_wbinlow,IMpippim_wbinhigh);
    IMpippim_IMnpim_wK0_n_wbin[ibin]->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
    IMpippim_IMnpim_wK0_n_wbin[ibin]->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  }
  
  IMpippim_IMnpim_wK0_n_woSpdia = new TH2F("IMpippim_IMnpim_wK0_n_woSpdia","IMpippim_IMnpim_wK0_n_woSpdia",nbinIMnpi,1,2,nbinpippim,0.,0.9);
  IMpippim_IMnpim_wK0_n_woSpdia->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  IMpippim_IMnpim_wK0_n_woSpdia->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  
  /*
  for(unsigned int ibin=0;ibin<nbintemplate;ibin++){
    float binlow=1.0+(float)ibin*1./nbintemplate;
    float binhigh=1.0+(float)(ibin+1.0)*1./nbintemplate;
    IMpippim_IMnpim_wK0_n_woSpdia_bin[ibin] = new TH2F(Form("IMpippim_IMnpim_wK0_n_woSpdia_bin%d",ibin),Form("IMpippim_IMnpim_wK0_n_woSpdia %0.2f-%0.2f",binlow,binhigh)
        ,nbinIMnpi,1,2,nbinpippim,0.,0.9);
    IMpippim_IMnpim_wK0_n_woSpdia_bin[ibin]->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
    IMpippim_IMnpim_wK0_n_woSpdia_bin[ibin]->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  }*/

  IMpippim_IMnpim_wK0orwSid_n = new TH2F("IMpippim_IMnpim_wK0orwSid_n","IMpippim_IMnpim_wK0orwSid_n",nbinIMnpi,1,2,nbinpippim,0.,0.9);
  IMpippim_IMnpim_wK0orwSid_n->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  IMpippim_IMnpim_wK0orwSid_n->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  
  for(unsigned int ibin=0;ibin<nbintemplate;ibin++){
    float binlow=1.0+(float)ibin*1./nbintemplate;
    float binhigh=1.0+(float)(ibin+1.0)*1./nbintemplate;
    IMpippim_IMnpim_wK0orwSid_n_bin[ibin] = new TH2F(Form("IMpippim_IMnpim_wK0orwSid_n_bin%d",ibin),Form("IMpippim_IMnpim_wK0orwSid_n %0.2f-%0.2f",binlow,binhigh)
        ,nbinIMnpi/4,1,2,nbinpippim/5,0.,0.9);
    IMpippim_IMnpim_wK0orwSid_n_bin[ibin]->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
    IMpippim_IMnpim_wK0orwSid_n_bin[ibin]->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  }

  for(unsigned int iwbin=0;iwbin<nwbin;iwbin++){
    IMpippim_IMnpim_wK0orwSid_n_wbin[iwbin] = new TH2F(Form("IMpippim_IMnpim_wK0orwSid_n_wbin%d",iwbin),
                                                       Form("IMpippim_IMnpim_wK0orwSid_n %0.2f-%0.2f",wbinlow[iwbin],wbinhigh[iwbin]),
                                                       nbinIMnpim_wbin,IMnpim_wbinlow,IMnpim_wbinhigh,nbinIMpippim_wbin,IMpippim_wbinlow,IMpippim_wbinhigh);
    IMpippim_IMnpim_wK0orwSid_n_wbin[iwbin]->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
    IMpippim_IMnpim_wK0orwSid_n_wbin[iwbin]->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  }

  IMpippim_IMnpim_wK0orwSid_n_woSpdia = new TH2F("IMpippim_IMnpim_wK0orwSid_n_woSpdia","IMpippim_IMnpim_wK0orwSid_n_woSpdia",nbinIMnpi,1,2,nbinpippim,0.,0.9);
  IMpippim_IMnpim_wK0orwSid_n_woSpdia->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  IMpippim_IMnpim_wK0orwSid_n_woSpdia->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  
  /*
  for(unsigned int ibin=0;ibin<nbintemplate;ibin++){
    float binlow=1.0+(float)ibin*1./nbintemplate;
    float binhigh=1.0+(float)(ibin+1.0)*1./nbintemplate;
    IMpippim_IMnpim_wK0orwSid_n_woSpdia_bin[ibin] = new TH2F(Form("IMpippim_IMnpim_wK0orwSid_n_woSpdia_bin%d",ibin),Form("IMpippim_IMnpim_wK0orwSid_n_woSpdia %0.2f-%0.2f",binlow,binhigh)
        ,nbinIMnpi,1,2,nbinpippim,0.,0.9);
    IMpippim_IMnpim_wK0orwSid_n_woSpdia_bin[ibin]->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
    IMpippim_IMnpim_wK0orwSid_n_woSpdia_bin[ibin]->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  }*/


  IMpippim_IMnpim_woK0_wSid_n = new TH2F("IMpippim_IMnpim_woK0_wSid_n","IMpippim_IMnpim_woK0_wSid_n",nbinIMnpi,1,2,nbinpippim,0.,0.9);
  IMpippim_IMnpim_woK0_wSid_n->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  IMpippim_IMnpim_woK0_wSid_n->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");

  for(unsigned int iwbin=0;iwbin<nwbin;iwbin++){
    IMpippim_IMnpim_woK0_wSid_n_woSp_wbin[iwbin] = new TH2F(
      Form("IMpippim_IMnpim_woK0_wSid_n_woSp_wbin%d",iwbin),
      Form("IMpippim_IMnpim_woK0_wSid_n_woSp %0.2f-%0.2f",wbinlow[iwbin],wbinhigh[iwbin]),
      nbinIMnpim_wbin,IMnpim_wbinlow,IMnpim_wbinhigh,nbinIMpippim_wbin,IMpippim_wbinlow,IMpippim_wbinhigh);
    IMpippim_IMnpim_woK0_wSid_n_woSp_wbin[iwbin]->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
    IMpippim_IMnpim_woK0_wSid_n_woSp_wbin[iwbin]->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
    
    IMpippim_IMnpim_woK0_wSid_n_woSm_wbin[iwbin] = new TH2F(
      Form("IMpippim_IMnpim_woK0_wSid_n_woSm_wbin%d",iwbin),
      Form("IMpippim_IMnpim_woK0_wSid_n_woSm %0.2f-%0.2f",wbinlow[iwbin],wbinhigh[iwbin]),
      nbinIMnpim_wbin,IMnpim_wbinlow,IMnpim_wbinhigh,nbinIMpippim_wbin,IMpippim_wbinlow,IMpippim_wbinhigh);
    IMpippim_IMnpim_woK0_wSid_n_woSm_wbin[iwbin]->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
    IMpippim_IMnpim_woK0_wSid_n_woSm_wbin[iwbin]->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  }

  IMpippim_IMnpim_n_cross = new TH2F("IMpippim_IMnpim_n_cross","IMpippim_IMnpim_n_cross",nbinIMnpi,1,2,nbinpippim,0.,0.9);
  IMpippim_IMnpim_n_cross->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  IMpippim_IMnpim_n_cross->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  
  IMpippim_IMnpip_woK0_woSid_won = new TH2F("IMpippim_IMnpip_woK0_woSid_won","IMpippim_IMnpip_woK0_woSid_won",nbinIMnpi,1,2,nbinpippim,0.,0.9);
  IMpippim_IMnpip_woK0_woSid_won->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMpippim_IMnpip_woK0_woSid_won->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");

  IMpippim_IMnpim_woK0_woSid_won = new TH2F("IMpippim_IMnpim_woK0_woSid_won","IMpippim_IMnpim_woK0_woSid_won",nbinIMnpi,1,2,nbinpippim,0.,0.9);
  IMpippim_IMnpim_woK0_woSid_won->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  IMpippim_IMnpim_woK0_woSid_won->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  
  for(unsigned int iwbin=0;iwbin<nwbin;iwbin++){
    IMpippim_IMnpip_wK0_woSid_n_wbin[iwbin] = new TH2F(Form("IMpippim_IMnpip_wK0_woSid_n_wbin%d",iwbin),
                                                       Form("IMpippim_IMnpip_wK0_woSid_n %0.2f-%0.2f",wbinlow[iwbin],wbinhigh[iwbin]),
                                                       nbinIMnpip_wbin,IMnpip_wbinlow,IMnpip_wbinhigh,nbinIMpippim_wbin,IMpippim_wbinlow,IMpippim_wbinhigh);
    IMpippim_IMnpip_wK0_woSid_n_wbin[iwbin]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
    IMpippim_IMnpip_wK0_woSid_n_wbin[iwbin]->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
    
    IMpippim_IMnpim_wK0_woSid_n_wbin[iwbin] = new TH2F(Form("IMpippim_IMnpim_wK0_woSid_n_wbin%d",iwbin),
                                                       Form("IMpippim_IMnpim_wK0_woSid_n %0.2f-%0.2f",wbinlow[iwbin],wbinhigh[iwbin]),
                                                       nbinIMnpim_wbin,IMnpim_wbinlow,IMnpim_wbinhigh,nbinIMpippim_wbin,IMpippim_wbinlow,IMpippim_wbinhigh);
    IMpippim_IMnpim_wK0_woSid_n_wbin[iwbin]->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
    IMpippim_IMnpim_wK0_woSid_n_wbin[iwbin]->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  }


  IMpippim_IMnpip_wK0_woSid_won = new TH2F("IMpippim_IMnpip_wK0_woSid_won","IMpippim_IMnpip_wK0_woSid_won",nbinIMnpi,1,2,nbinpippim,0.,0.9);
  IMpippim_IMnpip_wK0_woSid_won->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMpippim_IMnpip_wK0_woSid_won->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");

  IMpippim_IMnpim_wK0_woSid_won = new TH2F("IMpippim_IMnpim_wK0_woSid_won","IMpippim_IMnpim_wK0_woSid_won",nbinIMnpi,1,2,nbinpippim,0.,0.9);
  IMpippim_IMnpim_wK0_woSid_won->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  IMpippim_IMnpim_wK0_woSid_won->SetYTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");

  Mompippim_IMnpipi_dE_wK0_n = new TH2F("Mompippim_IMnpipi_dE_wK0_n","Mompippim_IMnpipi_dE_wK0_n",nbinIMnpipi,IMnpipilow,IMnpipihi,50,0,1);
  Mompippim_IMnpipi_dE_wK0_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  Mompippim_IMnpipi_dE_wK0_n->SetYTitle("Mom(#pi^{+}#pi^{-}) [GeV/c]");

  Mompippim_nmom_dE_wK0_n = new TH2F("Mompippim_nmom_dE_wK0_n","Mompippim_nmom_dE_wK0_n",nbinnmom,0.,1.0,50,0,1);
  Mompippim_nmom_dE_wK0_n->SetXTitle("nmom [GeV/c]");
  Mompippim_nmom_dE_wK0_n->SetYTitle("Mom(#pi^{+}#pi^{-}) [GeV/c]");

  Mompippim_IMnpipi_dE_wK0_n_Sp = new TH2F("Mompippim_IMnpipi_dE_wK0_n_Sp","Mompippim_IMnpipi_dE_wK0_n_Sp",nbinIMnpipi,IMnpipilow,IMnpipihi,50,0,1);
  Mompippim_IMnpipi_dE_wK0_n_Sp->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  Mompippim_IMnpipi_dE_wK0_n_Sp->SetYTitle("Mom(#pi^{+}#pi^{-}) [GeV/c]");

  Mompippim_IMnpipi_dE_wK0_n_Sm = new TH2F("Mompippim_IMnpipi_dE_wK0_n_Sm","Mompippim_IMnpipi_dE_wK0_n_Sm",nbinIMnpipi,IMnpipilow,IMnpipihi,50,0,1);
  Mompippim_IMnpipi_dE_wK0_n_Sm->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  Mompippim_IMnpipi_dE_wK0_n_Sm->SetYTitle("Mom(#pi^{+}#pi^{-}) [GeV/c]");

  diff2d_Phipippim_Phinpim = new TH2F("diff2d_Phipippim_Phinpim","diff2d_Phipippim_Phinpim",36,-1.*TMath::Pi(),TMath::Pi(),36,-1.*TMath::Pi(),TMath::Pi());
  diff2d_Phipippim_Phinpim->SetXTitle("n_{can.} - #pi^{-} hit (phi) [radian]");
  diff2d_Phipippim_Phinpim->SetYTitle("#pi^{+} - #pi^{-} hit (phi) [radian]");

  //diff2d_Zpippim_Znpim = new TH2F("diff2d_Zpippim_Znpim","diff2d_Zpippim_Znpim",100,-100,100,100,-100,100);
  //diff2d_Zpippim_Znpim->SetXTitle("n_{can.} - #pi^{-} hit (z) [cm]");
  //diff2d_Zpippim_Znpim->SetYTitle("#pi^{+} - #pi^{-} hit (z) [cm]");
  
  //diff2d_Phipippim_Phinpim_woSid_won = new TH2F("diff2d_Phipippim_Phinpim_woSid_won","diff2d_Phipippim_Phinpim_woSid_won",36,-1.*TMath::Pi(),TMath::Pi(),36,-1.*TMath::Pi(),TMath::Pi());
  //diff2d_Phipippim_Phinpim_woSid_won->SetXTitle("n_{fake} - #pi^{-} hit (phi) [radian]");
  //diff2d_Phipippim_Phinpim_woSid_won->SetYTitle("#pi^{+} - #pi^{-} hit (phi) [radian]");

  diff2d_Zpippim_Znpim_woSid_won = new TH2F("diff2d_Zpippim_Znpim_woSid_won","diff2d_Zpippim_Znpim_woSid_won",100,-100,100,100,-100,100);
  diff2d_Zpippim_Znpim_woSid_won->SetXTitle("n_{fake} - #pi^{-} hit (z) [cm]");
  diff2d_Zpippim_Znpim_woSid_won->SetYTitle("#pi^{+} - #pi^{-} hit (z) [cm]");

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
  
  //diff2d_CDC_CDH_pim_woSid_won = new TH2F("diff2d_CDC_CDH_pim_woSid_won","diff2d_CDC_CDH_pim_woSid_won",100,-1.*TMath::Pi(),TMath::Pi(),100,-100,100);
  //diff2d_CDC_CDH_pim_woSid_won->SetXTitle("CDH hit - #pi^{-} track (phi) [radian]");
  //diff2d_CDC_CDH_pim_woSid_won->SetYTitle("CDH hit - #pi^{-} track (z) [cm]");
  
  diff2d_CDC_CDH_pim_phi_tof_woSid_won = new TH2F("diff2d_CDC_CDH_pim_phi_tof_woSid_won","diff2d_CDC_CDH_pim_phi_tof_woSid_won",100,-1.*TMath::Pi(),TMath::Pi(),120,-10,50);
  diff2d_CDC_CDH_pim_phi_tof_woSid_won->SetXTitle("CDH hit - #pi^{-} track (phi) [radian]");
  diff2d_CDC_CDH_pim_phi_tof_woSid_won->SetYTitle("CDH hit - #pi^{-} TOF [nsec]");
  
  diff2d_CDC_CDH_pim_z_tof_woSid_won = new TH2F("diff2d_CDC_CDH_pim_z_tof_woSid_won","diff2d_CDC_CDH_pim_z_tof_woSid_won",100,-100,100,120,-10,50);
  diff2d_CDC_CDH_pim_z_tof_woSid_won->SetXTitle("CDH hit - #pi^{-} track (z) [cm]");
  diff2d_CDC_CDH_pim_z_tof_woSid_won->SetYTitle("CDH hit - #pi^{-} TOF [nsec]");

  //diff2d_CDC_CDH_pip_woSid_won = new TH2F("diff2d_CDC_CDH_pip_woSid_won","diff2d_CDC_CDH_pip_woSid_won",100,-1.*TMath::Pi(),TMath::Pi(),100,-100,100);
  //diff2d_CDC_CDH_pip_woSid_won->SetXTitle("CDH hit - #pi^{+} track (phi) [radian]");
  //diff2d_CDC_CDH_pip_woSid_won->SetYTitle("CDH hit - #pi^{+} track (z) [cm]");
  
  diff2d_CDC_CDH_pip_phi_tof_woSid_won = new TH2F("diff2d_CDC_CDH_pip_phi_tof_woSid_won","diff2d_CDC_CDH_pip_phi_tof_woSid_won",100,-1.*TMath::Pi(),TMath::Pi(),120,-10,50);
  diff2d_CDC_CDH_pip_phi_tof_woSid_won->SetXTitle("CDH hit - #pi^{+} track (phi) [radian]");
  diff2d_CDC_CDH_pip_phi_tof_woSid_won->SetYTitle("CDH hit - #pi^{+} TOF [nsec]");
  
  diff2d_CDC_CDH_pip_z_tof_woSid_won = new TH2F("diff2d_CDC_CDH_pip_z_tof_woSid_won","diff2d_CDC_CDH_pip_z_tof_woSid_won",100,-100,100,120,-10,50);
  diff2d_CDC_CDH_pip_z_tof_woSid_won->SetXTitle("CDH hit - #pi^{-} track (z) [cm]");
  diff2d_CDC_CDH_pip_z_tof_woSid_won->SetYTitle("CDH hit - #pi^{-} TOF [nsec]");

  diff2d_CDC_CDH_pim_wK0 = new TH2F("diff2d_CDC_CDH_pim_wK0","diff2d_CDC_CDH_pim_wK0",100,-1.*TMath::Pi(),TMath::Pi(),100,-100,100);
  diff2d_CDC_CDH_pim_wK0->SetXTitle("CDH hit - #pi^{-} track (phi) [radian]");
  diff2d_CDC_CDH_pim_wK0->SetYTitle("CDH hit - #pi^{-} track (z) [cm]");

  diff2d_CDC_CDH_pip_wK0 = new TH2F("diff2d_CDC_CDH_pip_wK0","diff2d_CDC_CDH_pip_wK0",100,-1.*TMath::Pi(),TMath::Pi(),100,-100,100);
  diff2d_CDC_CDH_pip_wK0->SetXTitle("CDH hit - #pi^{+} track (phi) [radian]");
  diff2d_CDC_CDH_pip_wK0->SetYTitle("CDH hit - #pi^{+} track (z) [cm]");

  diff2d_CDC_CDH_pim_woK0_wSid = new TH2F("diff2d_CDC_CDH_pim_woK0_wSid","diff2d_CDC_CDH_pim_woK0_wSid",100,-1.*TMath::Pi(),TMath::Pi(),100,-100,100);
  diff2d_CDC_CDH_pim_woK0_wSid->SetXTitle("CDH hit - #pi^{-} track (phi) [radian]");
  diff2d_CDC_CDH_pim_woK0_wSid->SetYTitle("CDH hit - #pi^{-} track (z) [cm]");

  diff2d_CDC_CDH_pip_woK0_wSid = new TH2F("diff2d_CDC_CDH_pip_woK0_wSid","diff2d_CDC_CDH_pip_woK0_wSid",100,-1.*TMath::Pi(),TMath::Pi(),100,-100,100);
  diff2d_CDC_CDH_pip_woK0_wSid->SetXTitle("CDH hit - #pi^{+} track (phi) [radian]");
  diff2d_CDC_CDH_pip_woK0_wSid->SetYTitle("CDH hit - #pi^{+} track (z) [cm]");

  diff2d_CDC_CDH_pim_woK0_wSid_n = new TH2F("diff2d_CDC_CDH_pim_woK0_wSid_n","diff2d_CDC_CDH_pim_woK0_wSid_n",100,-1.*TMath::Pi(),TMath::Pi(),100,-100,100);
  diff2d_CDC_CDH_pim_woK0_wSid_n->SetXTitle("CDH hit - #pi^{-} track (phi) [radian]");
  diff2d_CDC_CDH_pim_woK0_wSid_n->SetYTitle("CDH hit - #pi^{-} track (z) [cm]");

  diff2d_CDC_CDH_pip_woK0_wSid_n = new TH2F("diff2d_CDC_CDH_pip_woK0_wSid_n","diff2d_CDC_CDH_pip_woK0_wSid_n",100,-1.*TMath::Pi(),TMath::Pi(),100,-100,100);
  diff2d_CDC_CDH_pip_woK0_wSid_n->SetXTitle("CDH hit - #pi^{+} track (phi) [radian]");
  diff2d_CDC_CDH_pip_woK0_wSid_n->SetYTitle("CDH hit - #pi^{+} track (z) [cm]");

  dE_diffphi_CDC_CDH_pim = new TH2F("dE_diffphi_CDC_CDH_pim","dE_diffphi_CDC_CDH_pim",100,-1.*TMath::Pi(),TMath::Pi(),nbindE,0,50);
  dE_diffphi_CDC_CDH_pim->SetXTitle("CDH hit - #pi^{-} track (phi) [radian]");
  dE_diffphi_CDC_CDH_pim->SetYTitle("CDH dE [MeVee]");

  dE_diffphi_CDC_CDH_pip = new TH2F("dE_diffphi_CDC_CDH_pip","dE_diffphi_CDC_CDH_pip",100,-1.*TMath::Pi(),TMath::Pi(),nbindE,0,50);
  dE_diffphi_CDC_CDH_pip->SetXTitle("CDH hit - #pi^{+} track (phi) [radian]");
  dE_diffphi_CDC_CDH_pip->SetYTitle("CDH dE [MeVee]");

  //MMnmiss_dE = new TH2F("MMnmiss_dE","MMnmiss_dE",nbindE,0,50,100,0.4,1.9);
  //MMnmiss_dE->SetXTitle("CDH dE [MeVee]");
  //MMnmiss_dE->SetYTitle("Miss. Mass [GeV/c^{2}]");

  MMnmiss_diffphi_CDC_CDH_pim
    = new TH2F("MMnmiss_diffphi_CDC_CDH_pim","MMnmiss_diffphi_CDC_CDH_pim",100,-1.*TMath::Pi(),TMath::Pi(),100,0.4,1.9);
  MMnmiss_diffphi_CDC_CDH_pim->SetXTitle("CDH hit - #pi^{-} track (phi) [radian]");
  MMnmiss_diffphi_CDC_CDH_pim->SetYTitle("Miss. Mass [GeV/c^{2}]");

  MMnmiss_diffphi_CDC_CDH_pip
    = new TH2F("MMnmiss_diffphi_CDC_CDH_pip","MMnmiss_diffphi_CDC_CDH_pip",100,-1.*TMath::Pi(),TMath::Pi(),100,0.4,1.9);
  MMnmiss_diffphi_CDC_CDH_pip->SetXTitle("CDH hit - #pi^{+} track (phi) [radian]");
  MMnmiss_diffphi_CDC_CDH_pip->SetYTitle("Miss. Mass [GeV/c^{2}]");

  MMnmiss_diffphi_CDC_CDH_pim_woK0_wSid
    = new TH2F("MMnmiss_diffphi_CDC_CDH_pim_woK0_wSid","MMnmiss_diffphi_CDC_CDH_pim_woK0_wSid",100,-1.*TMath::Pi(),TMath::Pi(),100,0.4,1.9);
  MMnmiss_diffphi_CDC_CDH_pim_woK0_wSid->SetXTitle("CDH hit - #pi^{-} track (phi) [radian]");
  MMnmiss_diffphi_CDC_CDH_pim_woK0_wSid->SetYTitle("Miss. Mass [GeV/c^{2}]");

  MMnmiss_diffphi_CDC_CDH_pip_woK0_wSid
    = new TH2F("MMnmiss_diffphi_CDC_CDH_pip_woK0_wSid","MMnmiss_diffphi_CDC_CDH_pip_woK0_wSid",100,-1.*TMath::Pi(),TMath::Pi(),100,0.4,1.9);
  MMnmiss_diffphi_CDC_CDH_pip_woK0_wSid->SetXTitle("CDH hit - #pi^{+} track (phi) [radian]");
  MMnmiss_diffphi_CDC_CDH_pip_woK0_wSid->SetYTitle("Miss. Mass [GeV/c^{2}]");

  MMnmiss_diffz_CDC_CDH_pim_woK0_wSid
    = new TH2F("MMnmiss_diffz_CDC_CDH_pim_woK0_wSid","MMnmiss_diffz_CDC_CDH_pim_woK0_wSid",100,-100,100,100,0.4,1.9);
  MMnmiss_diffz_CDC_CDH_pim_woK0_wSid->SetXTitle("CDH hit - #pi^{-} track (z) [cm]");
  MMnmiss_diffz_CDC_CDH_pim_woK0_wSid->SetYTitle("Miss. Mass [GeV/c^{2}]");

  MMnmiss_diffz_CDC_CDH_pip_woK0_wSid
    = new TH2F("MMnmiss_diffz_CDC_CDH_pip_woK0_wSid","MMnmiss_diffz_CDC_CDH_pip_woK0_wSid",100,-100,100,100,0.4,1.9);
  MMnmiss_diffz_CDC_CDH_pip_woK0_wSid->SetXTitle("CDH hit - #pi^{+} track (z) [cm]");
  MMnmiss_diffz_CDC_CDH_pip_woK0_wSid->SetYTitle("Miss. Mass [GeV/c^{2}]");

  pimmom_diffphi_CDC_CDH_pim = new TH2F("pimmom_diffphi_CDC_CDH_pim","pimmom_diffphi_CDC_CDH_pim",100,-1.*TMath::Pi(),TMath::Pi(), 200,0,1);
  pimmom_diffphi_CDC_CDH_pim->SetXTitle("CDH hit - #pi^{-} track (phi) [radian]");
  pimmom_diffphi_CDC_CDH_pim->SetYTitle("#pi^{-} mom. [GeV/c]");
  
  pimmom_diffphi_CDC_CDH_pim_woSid_won = new TH2F("pimmom_diffphi_CDC_CDH_pim_woSid_won","pimmom_diffphi_CDC_CDH_pim_woSid_won",100,-1.*TMath::Pi(),TMath::Pi(), 200,0,1);
  pimmom_diffphi_CDC_CDH_pim_woSid_won->SetXTitle("CDH hit - #pi^{-} track (phi) [radian]");
  pimmom_diffphi_CDC_CDH_pim_woSid_won->SetYTitle("#pi^{-} mom. [GeV/c]");

  pipmom_diffphi_CDC_CDH_pip = new TH2F("pipmom_diffphi_CDC_CDH_pip","pipmom_diffphi_CDC_CDH_pip",100,-1.*TMath::Pi(),TMath::Pi(), 200,0,1);
  pipmom_diffphi_CDC_CDH_pip->SetXTitle("CDH hit - #pi^{+} track (phi) [radian]");
  pipmom_diffphi_CDC_CDH_pip->SetYTitle("#pi^{+} mom. [GeV/c]");
  
  pipmom_diffphi_CDC_CDH_pip_woSid_won = new TH2F("pipmom_diffphi_CDC_CDH_pip_woSid_won","pipmom_diffphi_CDC_CDH_pip_woSid_won",100,-1.*TMath::Pi(),TMath::Pi(), 200,0,1);
  pipmom_diffphi_CDC_CDH_pip_woSid_won->SetXTitle("CDH hit - #pi^{+} track (phi) [radian]");
  pipmom_diffphi_CDC_CDH_pip_woSid_won->SetYTitle("#pi^{+} mom. [GeV/c]");

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

  pimmom_diffphi_CDC_CDH_pim_woK0_wSid = new TH2F("pimmom_diffphi_CDC_CDH_pim_woK0_wSid","pimmom_diffphi_CDC_CDH_pim_woK0_wSid",100,-1.*TMath::Pi(),TMath::Pi(), 200,0,1);
  pimmom_diffphi_CDC_CDH_pim_woK0_wSid->SetXTitle("CDH hit - #pi^{-} track (phi) [radian]");
  pimmom_diffphi_CDC_CDH_pim_woK0_wSid->SetYTitle("#pi^{-} mom. [GeV/c]");

  pipmom_diffphi_CDC_CDH_pip_woK0_wSid = new TH2F("pipmom_diffphi_CDC_CDH_pip_woK0_wSid","pipmom_diffphi_CDC_CDH_pip_woK0_wSid",100,-1.*TMath::Pi(),TMath::Pi(), 200,0,1);
  pipmom_diffphi_CDC_CDH_pip_woK0_wSid->SetXTitle("CDH hit - #pi^{+} track (phi) [radian]");
  pipmom_diffphi_CDC_CDH_pip_woK0_wSid->SetYTitle("#pi^{+} mom. [GeV/c]");


  pimmom_diffz_CDC_CDH_pim_woK0_wSid = new TH2F("pimmom_diffz_CDC_CDH_pim_woK0_wSid","pimmom_diffz_CDC_CDH_pim_woK0_wSid",100,-100,100, 200,0,1);
  pimmom_diffz_CDC_CDH_pim_woK0_wSid->SetXTitle("CDH hit - #pi^{-} track (z) [cm]");
  pimmom_diffz_CDC_CDH_pim_woK0_wSid->SetYTitle("#pi^{-} mom. [GeV/c]");

  pipmom_diffz_CDC_CDH_pip_woK0_wSid = new TH2F("pipmom_diffz_CDC_CDH_pip_woK0_wSid","pipmom_diffz_CDC_CDH_pip_woK0_wSid",100,-100,100, 200,0,1);
  pipmom_diffz_CDC_CDH_pip_woK0_wSid->SetXTitle("CDH hit - #pi^{+} track (z) [cm]");
  pipmom_diffz_CDC_CDH_pip_woK0_wSid->SetYTitle("#pi^{+} mom. [GeV/c]");


  nmom_diffphi_CDC_CDH_pim_woK0_wSid = new TH2F("nmom_diffphi_CDC_CDH_pim_woK0_wSid","nmom_diffphi_CDC_CDH_pim_woK0_wSid",100,-1.*TMath::Pi(),TMath::Pi(),100,0,1);
  nmom_diffphi_CDC_CDH_pim_woK0_wSid->SetXTitle("CDH hit - #pi^{-} track (phi) [radian]");
  nmom_diffphi_CDC_CDH_pim_woK0_wSid->SetYTitle("nmom.  [GeV/c]");

  nmom_diffphi_CDC_CDH_pip_woK0_wSid = new TH2F("nmom_diffphi_CDC_CDH_pip_woK0_wSid","nmom_diffphi_CDC_CDH_pip_woK0_wSid",100,-1.*TMath::Pi(),TMath::Pi(),100,0,1);
  nmom_diffphi_CDC_CDH_pip_woK0_wSid->SetXTitle("CDH hit - #pi^{+} track (phi) [radian]");
  nmom_diffphi_CDC_CDH_pip_woK0_wSid->SetYTitle("nmom.  [GeV/c]");


  nmom_diffz_CDC_CDH_pim_woK0_wSid = new TH2F("nmom_diffz_CDC_CDH_pim_woK0_wSid","nmom_diffz_CDC_CDH_pim_woK0_wSid",100,-100,100,100,0,1);
  nmom_diffz_CDC_CDH_pim_woK0_wSid->SetXTitle("CDH hit - #pi^{-} track (z) [cm]");
  nmom_diffz_CDC_CDH_pim_woK0_wSid->SetYTitle("nmom.  [GeV/c]");

  nmom_diffz_CDC_CDH_pip_woK0_wSid = new TH2F("nmom_diffz_CDC_CDH_pip_woK0_wSid","nmom_diffz_CDC_CDH_pip_woK0_wSid",100,-100,100,100,0,1);
  nmom_diffz_CDC_CDH_pip_woK0_wSid->SetXTitle("CDH hit - #pi^{+} track (z) [cm]");
  nmom_diffz_CDC_CDH_pip_woK0_wSid->SetYTitle("nmom.  [GeV/c]");
  
 
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
  TH2F* q_IMnpipi_wSid_n_mc;
  q_IMnpipi_wSid_n_mc = new TH2F(Form("q_IMnpipi_wSid_n_mc"),Form("q_IMnpipi_wSid_n_mc"),nbinIMnpipi,IMnpipilow,IMnpipihi,nbinq,0,1.5);
  q_IMnpipi_wSid_n_mc->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wSid_n_mc->SetYTitle("Mom. Transfer [GeV/c]");

  TH2F* q_IMnpipi_woK0_wSid_n_mc;
  q_IMnpipi_woK0_wSid_n_mc = new TH2F(Form("q_IMnpipi_woK0_wSid_n_mc"),Form("q_IMnpipi_woK0_wSid_n_mc"),nbinIMnpipi,IMnpipilow,IMnpipihi,nbinq,0,1.5);
  q_IMnpipi_woK0_wSid_n_mc->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_wSid_n_mc->SetYTitle("Mom. Transfer [GeV/c]");

  TH2F* q_IMnpipi_wSid_n_Sp_mc;
  q_IMnpipi_wSid_n_Sp_mc = new TH2F(Form("q_IMnpipi_wSid_n_Sp_mc"),Form("q_IMnpipi_wSid_n_Sp_mc"),nbinIMnpipi,IMnpipilow,IMnpipihi,nbinq,0,1.5);
  q_IMnpipi_wSid_n_Sp_mc->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wSid_n_Sp_mc->SetYTitle("Mom. Transfer [GeV/c]");

  TH2F* q_IMnpipi_woK0_wSid_n_Sp_mc;
  q_IMnpipi_woK0_wSid_n_Sp_mc = new TH2F(Form("q_IMnpipi_woK0_wSid_n_Sp_mc"),Form("q_IMnpipi_woK0_wSid_n_Sp_mc"),nbinIMnpipi,IMnpipilow,IMnpipihi,nbinq,0,1.5);
  q_IMnpipi_woK0_wSid_n_Sp_mc->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_wSid_n_Sp_mc->SetYTitle("Mom. Transfer [GeV/c]");

  TH2F* q_IMnpipi_wSid_n_Sm_mc;
  q_IMnpipi_wSid_n_Sm_mc = new TH2F("q_IMnpipi_wSid_n_Sm_mc","q_IMnpipi_wSid_n_Sm_mc",nbinIMnpipi,IMnpipilow,IMnpipihi,nbinq,0,1.5);
  q_IMnpipi_wSid_n_Sm_mc->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wSid_n_Sm_mc->SetYTitle("Mom. Transfer [GeV/c]");

  TH2F* q_IMnpipi_woK0_wSid_n_Sm_mc = NULL;
  q_IMnpipi_woK0_wSid_n_Sm_mc = new TH2F("q_IMnpipi_woK0_wSid_n_Sm_mc","q_IMnpipi_woK0_wSid_n_Sm_mc",nbinIMnpipi,IMnpipilow,IMnpipihi,nbinq,0,1.5);
  q_IMnpipi_woK0_wSid_n_Sm_mc->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_wSid_n_Sm_mc->SetYTitle("Mom. Transfer [GeV/c]");


  TH2F* diff_IMnpipi_wSid_n = NULL;
  diff_IMnpipi_wSid_n = new TH2F("diff_IMnpipi_wSid_n","diff_IMnpipi_wSid_n",nbinIMnpipi,IMnpipilow,IMnpipihi,100,-1,1);
  diff_IMnpipi_wSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  diff_IMnpipi_wSid_n->SetYTitle("reco. - gen.  [GeV/c^{2}]");

  TH2F* diff_IMnpipi_woK0_wSid_n = NULL;
  diff_IMnpipi_woK0_wSid_n = new TH2F("diff_IMnpipi_woK0_wSid_n","diff_IMnpipi_woK0_wSid_n",nbinIMnpipi,IMnpipilow,IMnpipihi,100,-1,1);
  diff_IMnpipi_woK0_wSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  diff_IMnpipi_woK0_wSid_n->SetYTitle("reco. - gen.  [GeV/c^{2}]");

  TH2F* diff_IMnpipi_wSid_n_Sp = NULL;
  diff_IMnpipi_wSid_n_Sp = new TH2F("diff_IMnpipi_wSid_n_Sp","diff_IMnpipi_wSid_n_Sp",nbinIMnpipi,IMnpipilow,IMnpipihi,600,-0.3,0.3);
  diff_IMnpipi_wSid_n_Sp->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  diff_IMnpipi_wSid_n_Sp->SetYTitle("reco. - gen. [GeV/c^{2}]");

  TH2F* diff_IMnpipi_woK0_wSid_n_Sp = NULL;
  diff_IMnpipi_woK0_wSid_n_Sp = new TH2F("diff_IMnpipi_woK0_wSid_n_Sp","diff_IMnpipi_woK0_wSid_n_Sp",nbinIMnpipi,IMnpipilow,IMnpipihi,600,-0.3,0.3);
  diff_IMnpipi_woK0_wSid_n_Sp->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  diff_IMnpipi_woK0_wSid_n_Sp->SetYTitle("reco. - gen. [GeV/c^{2}]");

  TH2F* diff_q_wSid_n_Sp = NULL;
  diff_q_wSid_n_Sp = new TH2F("diff_q_wSid_n_Sp","diff_q_wSid_n_Sp",nbinq,0,1,600,-0.3,0.3);
  diff_q_wSid_n_Sp->SetXTitle("Mom. Transfer [GeV/c]");
  diff_q_wSid_n_Sp->SetYTitle("reco. - gen. [GeV/c^]");

  TH2F* diff_q_woK0_wSid_n_Sp = NULL;
  diff_q_woK0_wSid_n_Sp = new TH2F("diff_q_woK0_wSid_n_Sp","diff_q_woK0_wSid_n_Sp",nbinq,0,1,600,-0.3,0.3);
  diff_q_woK0_wSid_n_Sp->SetXTitle("Mom. Transfer [GeV/c]");
  diff_q_woK0_wSid_n_Sp->SetYTitle("reco. - gen. [GeV/c^]");

  TH2F* diff_IMnpipi_wSid_n_Sm = NULL;
  diff_IMnpipi_wSid_n_Sm = new TH2F("diff_IMnpipi_wSid_n_Sm","diff_IMnpipi_wSid_n_Sm",nbinIMnpipi,IMnpipilow,IMnpipihi,600,-0.3,0.3);
  diff_IMnpipi_wSid_n_Sm->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  diff_IMnpipi_wSid_n_Sm->SetYTitle("reco. - gen.  [GeV/c^{2}]");

  TH2F* diff_IMnpipi_woK0_wSid_n_Sm = NULL;
  diff_IMnpipi_woK0_wSid_n_Sm = new TH2F("diff_IMnpipi_woK0_wSid_n_Sm","diff_IMnpipi_woK0_wSid_n_Sm",nbinIMnpipi,IMnpipilow,IMnpipihi,600,-0.3,0.3);
  diff_IMnpipi_woK0_wSid_n_Sm->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  diff_IMnpipi_woK0_wSid_n_Sm->SetYTitle("reco. - gen.  [GeV/c^{2}]");

  TH2F* diff_q_wSid_n_Sm = NULL;
  diff_q_wSid_n_Sm = new TH2F(Form("diff_q_wSid_n_Sm"),Form("diff_q_wSid_n_Sm"),nbinq,0,1,600,-0.3,0.3);
  diff_q_wSid_n_Sm->SetXTitle("Mom. Transfer [GeV/c]");
  diff_q_wSid_n_Sm->SetYTitle("reco. - gen. [GeV/c^]");

  TH2F* diff_q_woK0_wSid_n_Sm = NULL;
  diff_q_woK0_wSid_n_Sm = new TH2F(Form("diff_q_woK0_wSid_n_Sm"),Form("diff_q_woK0_wSid_n_Sm"),nbinq,0,1,600,-0.3,0.3);
  diff_q_woK0_wSid_n_Sm->SetXTitle("Mom. Transfer [GeV/c]");
  diff_q_woK0_wSid_n_Sm->SetYTitle("reco. - gen. [GeV/c^]");


  std::cout << __LINE__ << std::endl;


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
  } else if(qvalcutflag == 2) {
    pt->AddText(Form("Mom. Transfer  > %0.2f ",anacuts::qvalcut));
  }else if(qvalcutflag ==3){
    pt->AddText(Form("nmiss thetacut 15 "));
  }
  
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
    TVector3 vtx_pip = *vtx_pip_cdc ;
    TVector3 vtx_pim = *vtx_pim_cdc ;
    // calc missing npi+pi- mass //
    TLorentzVector LVec_npipimiss = *LVec_target+*LVec_beam-*LVec_pip-*LVec_pim-*LVec_n;
    TLorentzVector LVec_npipimiss_vtx[2];
    if(UseKinFit) {
      if(!UseKinFitVal) {
        LVec_npipimiss_vtx[0] = *LVec_target+*LVec_beam_Sp-*LVec_pip-*LVec_pim-*LVec_n_Sp;
        LVec_npipimiss_vtx[1] = *LVec_target+*LVec_beam_Sm-*LVec_pip-*LVec_pim-*LVec_n_Sm;
      } else {
        LVec_npipimiss_vtx[0] = *LVec_target+*kfSpmode_mom_beam-*kfSpmode_mom_pip-*kfSpmode_mom_pim-*kfSpmode_mom_n;
        LVec_npipimiss_vtx[1] = *LVec_target+*kfSmmode_mom_beam-*kfSmmode_mom_pip-*kfSmmode_mom_pim-*kfSmmode_mom_n;
      }
    }
    double nmiss_mass = LVec_npipimiss.M();
    //if(nmiss_mass<0) continue;
    double nmiss_mass_vtx[2]= {LVec_npipimiss_vtx[0].M(),LVec_npipimiss_vtx[1].M()};
    double nmiss_mom = LVec_npipimiss.P();

    // calc cos(theta) of missing npi+pi- //
    TVector3 boost = (*LVec_target+*LVec_beam).BoostVector();
    TVector3 boost_vtx[2] = {(*LVec_target+*LVec_beam_Sp).BoostVector(),(*LVec_target+*LVec_beam_Sm).BoostVector()} ;
    double cos_nmisslab = LVec_npipimiss.Vect().Dot((*LVec_beam).Vect())/(LVec_npipimiss.Vect().Mag()*(*LVec_beam).Vect().Mag());
    double nmissthetalab = acos(cos_nmisslab);
    if(qvalcutflag==3 && (nmissthetalab>15.0/180.0*TMath::Pi())) continue;
    TLorentzVector LVec_npipimiss_CM = LVec_npipimiss;
    TLorentzVector LVec_beam_CM = *LVec_beam;
    LVec_npipimiss_CM.Boost(-boost);
    LVec_beam_CM.Boost(-boost);
    double cos_nmissCM = LVec_npipimiss_CM.Vect().Dot(LVec_beam_CM.Vect())/(LVec_npipimiss_CM.Vect().Mag()*LVec_beam_CM.Vect().Mag());
    //cos(theta) of n_cds in lab
    double cos_ncdslab = (SimFakemode || SimFakeK0mode) ?  (*LVec_n).CosTheta() : LVec_n->Vect().Dot(LVec_beam->Vect())/(LVec_n->Vect().Mag()*LVec_beam->Vect().Mag());
    TLorentzVector LVec_n_CM = *LVec_n;
    LVec_n_CM.Boost(-boost);
    double cos_ncdsCM = LVec_n_CM.Vect().Dot(LVec_beam_CM.Vect())/(LVec_n_CM.Vect().Mag()*LVec_beam_CM.Vect().Mag());

    double cos_pim = (SimFakemode || SimFakeK0mode) ? (*LVec_pim).CosTheta() : (*LVec_pim).Vect().Dot((*LVec_beam).Vect())/((*LVec_pim).Vect().Mag()*(*LVec_beam).Vect().Mag());
    double cos_pip = (SimFakemode || SimFakeK0mode) ? (*LVec_pip).CosTheta() : (*LVec_pip).Vect().Dot((*LVec_beam).Vect())/((*LVec_pip).Vect().Mag()*(*LVec_beam).Vect().Mag());
    double cos_pippim = (*LVec_pip).Vect().Dot((*LVec_pim).Vect())/((*LVec_pip).Vect().Mag()*(*LVec_pim).Vect().Mag());


    TLorentzVector LVec_pim_CM = *LVec_pim;
    LVec_pim_CM.Boost(-boost);
    double cos_pimCM = LVec_pim_CM.Vect().Dot(LVec_beam_CM.Vect())/(LVec_pim_CM.Vect().Mag()*LVec_beam_CM.Vect().Mag());

    TLorentzVector LVec_pip_CM = *LVec_pip;
    LVec_pip_CM.Boost(-boost);
    double cos_pipCM = LVec_pip_CM.Vect().Dot(LVec_beam_CM.Vect())/(LVec_pip_CM.Vect().Mag()*LVec_beam_CM.Vect().Mag());

    if(SimSpmode || SimSmmode || SimK0nnmode) {
      TVector3 boost_mc =  (*LVec_target+*mcmom_beam).BoostVector();
    }
    TLorentzVector qkn_mc;
    TLorentzVector mcmom_nmiss_calc;
    if(SimSpmode || SimSmmode || SimK0nnmode) {
      qkn_mc = *mcmom_beam-*mcmom_nmiss;
      mcmom_nmiss_calc = *LVec_target+*mcmom_beam-*mcmom_pip-*mcmom_pim-*mcmom_ncds;
    }
    TLorentzVector qkn_vtx[2];
    if(UseKinFit) {
      if(!UseKinFitVal) {
        qkn_vtx[0] = *LVec_beam_Sp-LVec_npipimiss_vtx[0];
        qkn_vtx[1] = *LVec_beam_Sm-LVec_npipimiss_vtx[1];
      } else {
        qkn_vtx[0] = *kfSpmode_mom_beam-LVec_npipimiss_vtx[0];
        qkn_vtx[1] = *kfSmmode_mom_beam-LVec_npipimiss_vtx[1];
      }
    }
    // calc pi+pi- //
    TLorentzVector LVec_pip_pim = *LVec_pip+*LVec_pim;
    TLorentzVector LVec_pip_pim_mc;
    if(SimSpmode || SimSmmode || SimK0nnmode) {
      LVec_pip_pim_mc = *mcmom_pip+*mcmom_pim;
    }
    TLorentzVector LVec_pip_pim_CM = LVec_pip_pim;
    LVec_pip_pim_CM.Boost(-boost);
    double cos_K0CM = LVec_pip_pim_CM.Vect().Dot(LVec_beam_CM.Vect())/(LVec_pip_pim_CM.Vect().Mag()*LVec_beam_CM.Vect().Mag());

    //K0 - nCDS cos CM frame
    double cos_K0_ncds_CM = LVec_pip_pim_CM.Vect().Dot(LVec_n_CM.Vect())/(LVec_pip_pim_CM.Vect().Mag()*LVec_n_CM.Vect().Mag());
    //missing particle - nCDS cos CM frame
    double cos_nmiss_ncds_CM = LVec_npipimiss_CM.Vect().Dot(LVec_n_CM.Vect())/(LVec_npipimiss_CM.Vect().Mag()*LVec_n_CM.Vect().Mag());

    //missing particle - K0 cos CM frame
    double cos_nmiss_pippim_CM = LVec_npipimiss_CM.Vect().Dot(LVec_pip_pim_CM.Vect())/(LVec_npipimiss_CM.Vect().Mag()*LVec_pip_pim_CM.Vect().Mag());

    // calc pi+n //
    TLorentzVector LVec_pip_n = *LVec_pip+*LVec_n;
    double phi_npip = (*LVec_n).Phi()-(*LVec_pip).Phi();
    double phi_pip = (*LVec_pip).Phi();
    double phi_n = (*LVec_n).Phi();

    if(phi_npip<-1.0*TMath::Pi()) phi_npip += 2.0*TMath::Pi();
    else if(phi_npip>1.0*TMath::Pi()) phi_npip -= 2.0*TMath::Pi();


    TLorentzVector LVec_pip_n_mc;
    if(SimSpmode || SimSmmode || SimK0nnmode) {
      LVec_pip_n_mc  = *mcmom_pip+*mcmom_ncds;
    }
    TLorentzVector LVec_pip_n_vtx[2];
    if(UseKinFit) {
      if(!UseKinFitVal) {
        LVec_pip_n_vtx[0] = *LVec_pip+*LVec_n_Sp;
        LVec_pip_n_vtx[1] = *LVec_pip+*LVec_n_Sm;
      } else {
        LVec_pip_n_vtx[0] = *kfSpmode_mom_pip+*kfSpmode_mom_n;
        LVec_pip_n_vtx[1] = *kfSmmode_mom_pip+*kfSmmode_mom_n;
      }
    }
    // calc pi-n //
    TLorentzVector LVec_pim_n = *LVec_pim+*LVec_n;
    double phi_npim = (*LVec_n).Phi()-(*LVec_pim).Phi();
    double phi_pim = (*LVec_pim).Phi();

    if(phi_npim<-1.0*TMath::Pi()) phi_npim += 2.0*TMath::Pi();
    else if(phi_npim>1.0*TMath::Pi()) phi_npim -= 2.0*TMath::Pi();
    TLorentzVector LVec_pim_n_mc;
    if(SimSpmode || SimSmmode || SimK0nnmode) {
      LVec_pim_n_mc = *mcmom_pim+*mcmom_ncds;
    }
    TLorentzVector LVec_pim_n_vtx[2];
    if(UseKinFit) {
      if(!UseKinFitVal) {
        LVec_pim_n_vtx[0]= *LVec_pim+*LVec_n_Sp;
        LVec_pim_n_vtx[1]= *LVec_pim+*LVec_n_Sm;
      } else {
        LVec_pim_n_vtx[0]= *kfSpmode_mom_pim+*kfSpmode_mom_n;
        LVec_pim_n_vtx[1]= *kfSmmode_mom_pim+*kfSmmode_mom_n;
      }
    }
    //std::cout << __LINE__ << std::endl;

    // calc missing pip+neutron //
    TLorentzVector LVec_pipmiss_nmiss = *LVec_target+*LVec_beam-*LVec_pim-*LVec_n;

    TLorentzVector LVec_pipmiss_nmiss_mc;
    if(SimSpmode || SimSmmode || SimK0nnmode) {
      LVec_pipmiss_nmiss_mc = *LVec_target+*mcmom_beam-*mcmom_pim-*mcmom_ncds;
    }
    TLorentzVector LVec_pipmiss_nmiss_vtx[2];
    if(UseKinFit) {
      if(!UseKinFitVal) {
        LVec_pipmiss_nmiss_vtx[0]= *LVec_target+*LVec_beam_Sp-*LVec_pim-*LVec_n_Sp;
        LVec_pipmiss_nmiss_vtx[1]= *LVec_target+*LVec_beam_Sm-*LVec_pim-*LVec_n_Sm;
      } else {
        LVec_pipmiss_nmiss_vtx[0]= *LVec_target+*kfSpmode_mom_beam-*kfSpmode_mom_pim-*kfSpmode_mom_n;
        LVec_pipmiss_nmiss_vtx[1]= *LVec_target+*kfSmmode_mom_beam-*kfSmmode_mom_pim-*kfSmmode_mom_n;
      }
    }

    // calc missing pim+neutron //
    TLorentzVector LVec_pimmiss_nmiss = *LVec_target+*LVec_beam-*LVec_pip-*LVec_n;
    TLorentzVector LVec_pimmiss_nmiss_mc;
    if(SimSpmode || SimSmmode || SimK0nnmode) {
      LVec_pimmiss_nmiss_mc = *LVec_target+*mcmom_beam-*mcmom_pip-*mcmom_ncds;
    }
    TLorentzVector LVec_pimmiss_nmiss_vtx[2];
    if(UseKinFit) {
      if(!UseKinFitVal) {
        LVec_pimmiss_nmiss_vtx[0] = *LVec_target+*LVec_beam_Sp-*LVec_pip-*LVec_n_Sp;
        LVec_pimmiss_nmiss_vtx[1] = *LVec_target+*LVec_beam_Sm-*LVec_pip-*LVec_n_Sm;
      } else {
        LVec_pimmiss_nmiss_vtx[0] = *LVec_target+*kfSpmode_mom_beam-*kfSpmode_mom_pip-*kfSpmode_mom_n;
        LVec_pimmiss_nmiss_vtx[1] = *LVec_target+*kfSmmode_mom_beam-*kfSmmode_mom_pip-*kfSmmode_mom_n;
      }
    }


    //calc pip_cds + missing neutron
    TLorentzVector LVec_pip_nmiss = *LVec_pip+LVec_npipimiss;
    TLorentzVector LVec_pip_nmiss_mc;
    if(SimSpmode || SimSmmode || SimK0nnmode){
      LVec_pip_nmiss_mc = *mcmom_pip+*mcmom_nmiss;
    }
    //calc pim_cds + missing neutron
    TLorentzVector LVec_pim_nmiss = *LVec_pim+LVec_npipimiss;
    TLorentzVector LVec_pim_nmiss_mc;
    if(SimSpmode || SimSmmode || SimK0nnmode){
      LVec_pim_nmiss_mc = *mcmom_pim+*mcmom_nmiss;
    }

    // calc pi+pi-n //
    TLorentzVector LVec_pip_pim_n = *LVec_pip+*LVec_pim+*LVec_n;
    TVector3 P_specn(0.,0.,0.);
    TLorentzVector LVec_specn;
    LVec_specn.SetVectM(P_specn,nMass);
    TLorentzVector LVec_K0bar = *LVec_pip+*LVec_pim+*LVec_n-LVec_specn;
    TLorentzVector LVec_pip_pim_n_mc;
    if(SimSpmode || SimSmmode || SimK0nnmode) {
      LVec_pip_pim_n_mc = *mcmom_pip+*mcmom_pim+*mcmom_ncds;
    }
    TLorentzVector LVec_pip_pim_n_vtx[2];
    if(UseKinFit) {
      if(!UseKinFitVal) {
        LVec_pip_pim_n_vtx[0] = *LVec_pip+*LVec_pim+*LVec_n_Sp;
        LVec_pip_pim_n_vtx[1] = *LVec_pip+*LVec_pim+*LVec_n_Sm;
      } else {
        LVec_pip_pim_n_vtx[0] = *kfSpmode_mom_pip+*kfSpmode_mom_pim+*kfSpmode_mom_n;
        LVec_pip_pim_n_vtx[1] = *kfSmmode_mom_pip+*kfSmmode_mom_pim+*kfSmmode_mom_n;
      }
    }
    TLorentzVector qkn = *LVec_beam-LVec_npipimiss;
    TLorentzVector LVec_pip_pim_n_CM = LVec_pip_pim_n;
    LVec_pip_pim_n_CM.Boost(-boost);
    //double cos_X = LVec_pip_pim_n_CM.Vect().Dot(LVec_beam_CM.Vect())/(LVec_pip_pim_n_CM.Vect().Mag()*LVec_beam_CM.Vect().Mag());


    if( (qkn.P()>=anacuts::qvalcut) && (qvalcutflag==1) ) continue;
    if( (qkn.P()<anacuts::qvalcut) && (qvalcutflag==2) ) continue;
    if( (*LVec_n).P()<anacuts::nmomcut) continue;
    //double chi2 = kfSpmode_chi2<kfSmmode_chi2 ? kfSpmode_chi2:kfSmmode_chi2;
    double pvalue = -9999;
    if(UseKinFit) pvalue = kfSmmode_pvalue<kfSpmode_pvalue ? kfSpmode_pvalue:kfSmmode_pvalue;
    if( (SimSpmode || SimSmmode) && 
        SimRejectFake && 
      (mcpattern!=2)) continue;
    if(SimK0nnmode && SimRejectFake && (mcpattern!=7) ) continue;
    //Filling generated info.

    //std::cout << __LINE__ << std::endl;
    TLorentzVector LVec_Sigma_react;
    TLorentzVector LVec_pi_react;
    TLorentzVector LVec_piSigma_react;
    TLorentzVector qkn_react;
    bool IsFakeN1 = false;
    bool IsFakeN2 = false;
    bool IsFakebyVTX = false;
    if(SimSpmode || SimSmmode || SimK0nnmode) {
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
      //std::cout << "nmiss P" <<  (*react_nmiss).P() << std::endl;
      qkn_react = LVec_beam_unit-*react_nmiss ;
      q_IMpiSigma_gen->Fill(LVec_piSigma_react.M()/1000.,qkn_react.P()/1000.);
      diff_IMnpim_reactmc->Fill(LVec_pim_n_mc.M()-LVec_Sigma_react.M()/1000.);
      IMnpim_MMnpim_mc->Fill(LVec_pim_nmiss_mc.M(),LVec_pim_n_mc.M());
      diff_IMnpip_reactmc->Fill(LVec_pip_n_mc.M()-LVec_Sigma_react.M()/1000.);
      IMnpip_MMnpip_mc->Fill(LVec_pip_nmiss_mc.M(),LVec_pip_n_mc.M());
      
      vtxr_vtxz_mc->Fill(mcncanvtxz,mcncanvtxr);
      if(mcpattern==2){
        vtxr_vtxz_mc_pat2->Fill(mcncanvtxz,mcncanvtxr);
      }
      else if(mcpattern==7){
        vtxr_vtxz_mc_pat7->Fill(mcncanvtxz,mcncanvtxr);
      }
      
      if(!IsFakebyVTX  &&  (*mcmom_nmiss).P()>0.01 && (*mcmom_nmiss).P()<100.0 && (*mcmom_ncds).P()>0.01 && (*mcmom_pip).P()>0.01 && (*mcmom_pim).P()>0.01) {
        IMnpim_MMnpim_mc_vtx->Fill(LVec_pim_nmiss_mc.M(),LVec_pim_n_mc.M());
        IMnpip_MMnpip_mc_vtx->Fill(LVec_pip_nmiss_mc.M(),LVec_pip_n_mc.M());
        //true neutron from Sigma decay case
        if(mcpattern==2){
          IMnpim_MMnpim_mc_vtx_pat2->Fill(LVec_pim_nmiss_mc.M(),LVec_pim_n_mc.M());
          IMnpip_MMnpip_mc_vtx_pat2->Fill(LVec_pip_nmiss_mc.M(),LVec_pip_n_mc.M());
        }
        //true initial neutron case
        if(mcpattern==7){
          IMnpim_MMnpim_mc_vtx_pat7->Fill(LVec_pim_nmiss_mc.M(),LVec_pim_n_mc.M());
          IMnpip_MMnpip_mc_vtx_pat7->Fill(LVec_pip_nmiss_mc.M(),LVec_pip_n_mc.M());
        }
      }
     
      
      diff_nmiss_reactmc->Fill((*mcmom_nmiss).P()-(*react_nmiss).P()/1000.);
      diff2D_MMnmiss_IMnpim_reactmc->Fill(LVec_pim_n_mc.M()-LVec_Sigma_react.M()/1000.,(*mcmom_nmiss).M()-(*react_nmiss).M()/1000.);
      diff2D_MMnmiss_IMnpip_reactmc->Fill(LVec_pip_n_mc.M()-LVec_Sigma_react.M()/1000.,(*mcmom_nmiss).M()-(*react_nmiss).M()/1000.);
      IMnpim_IMnpip_mc->Fill(LVec_pip_n_mc.M(), LVec_pim_n_mc.M());
      nmom_IMnpip_mc->Fill((*mcmom_ncds).P(), LVec_pip_n_mc.M());
      nmom_IMnpim_mc->Fill((*mcmom_ncds).P(), LVec_pim_n_mc.M());
      //if(SimSpmode){
      //  if(fabs(LVec_pip_n_mc.M()-LVec_Sigma_react.M()/1000.) > 0.002) IsFakeN1 = true;
      //  if(fabs((*mcmom_nmiss).M()-(*react_nmiss).M()/1000.) > 0.002) IsFakeN1 = true;
      //}
      //if(SimSmmode){
      //  if(fabs(LVec_pim_n_mc.M()-LVec_Sigma_react.M()/1000.) > 0.002) IsFakeN1 = true;
      //  if(fabs((*mcmom_nmiss).M()-(*react_nmiss).M()/1000.) > 0.002) IsFakeN1 = true;
      //}
      //added angle check
      //if(!IsFakeN1){
      //  diff_cosnmiss_reactmc->Fill((*mcmom_nmiss).CosTheta()-(*react_nmiss).CosTheta());
      //  if(fabs((*mcmom_nmiss).CosTheta()-(*react_nmiss).CosTheta())>0.002) IsFakeN1 = true;
      //}
      if(SimSpmode){
        if(fabs(LVec_pip_n_mc.M()-LVec_Sigma_react.M()/1000.)>0.02) IsFakeN2=true;
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

      if(!IsFakebyVTX)MMnpim_MMnpip_mc->Fill(LVec_pip_nmiss_mc.M(),LVec_pim_nmiss_mc.M());
    }

    bool K0Flag=false;
    bool K0rejectFlag=false;
    bool MissNFlag=false;
    bool MissNwideFlag=false;
    bool NBetaOK=false;
    bool NdEOK=false;
    bool SigmaPFlag=false;
    bool SigmaMFlag=false;
    bool SigmawidePFlag=false;
    bool SigmawideMFlag=false;
    bool SigmaPMissNViciFlag=false;//select vicinity of the signal
    bool SigmaMMissNViciFlag=false;//select vicinity of the signal
    bool K0MissNViciFlag=false;//select vicinity of the signal
    bool SigmaPMissNViciextFlag=false;//select vicinity of the signal
    bool SigmaMMissNViciextFlag=false;//select vicinity of the signal
    bool SigmaPsideFlag[3][ngap];//pattern, gap
    for(int ipat=0; ipat<3; ipat++) {
      for(int igap=0; igap<ngap; igap++) {
        SigmaPsideFlag[ipat][igap]=false;
      }
    }
    bool SigmaPsideLowFlag[ngap];
    bool SigmaPsideHighFlag[ngap];
    for(int igap=0; igap<ngap; igap++) {
      SigmaPsideLowFlag[igap]=false;
      SigmaPsideHighFlag[igap]=false;
    }

    bool SigmaMsideFlag[3][ngap];//pattern,gap
    for(int ipat=0; ipat<3; ipat++) {
      for(int igap=0; igap<ngap; igap++) {
        SigmaMsideFlag[ipat][igap]=false;
      }
    }
    bool SigmaMsideLowFlag[ngap];
    bool SigmaMsideHighFlag[ngap];
    for(int igap=0; igap<ngap; igap++) {
      SigmaMsideHighFlag[igap]=false;
      SigmaMsideLowFlag[igap]=false;
    }

    bool SigmaPsideWideFlag[nzone];
    bool SigmaMsideWideFlag[nzone];
    for(int izone=0; izone<nzone; izone++) {
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

    for(int igap=0; igap<ngap; igap++) {
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

    for(int izone=0; izone<nzone; izone++) {
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

    TVector3 pippim_CA = 0.5*(*CA_pip+*CA_pim);
    TVector3 beam_pippim_CA = 0.5*(*vtx_pip_beam+*vtx_pim_beam);
    double dca_pipibeam = (pippim_CA-beam_pippim_CA).Perp();

    //-- neutron-ID, K0 and missing neutron selection --//
    if(anacuts::beta_MIN<NeutralBetaCDH &&  NeutralBetaCDH<anacuts::beta_MAX  ) NBetaOK=true;
    if(anacuts::dE_MIN<dE) NdEOK=true;
    double MassNPip= (*LVec_n+*LVec_pip).M();
    double MassNPim= (*LVec_n+*LVec_pim).M();

    TVector3 diffpim = (*CDH_Pos)-(*CDH_Pos_pim);
    double diffPhinpim = (*CDH_Pos).Phi()-(*CDH_Pos_pim).Phi();
    const double difftofnpim = tofn - tofpim;
    if(diffPhinpim<-1.0*TMath::Pi()) diffPhinpim += 2.0*TMath::Pi();
    else if(diffPhinpim>1.0*TMath::Pi()) diffPhinpim -= 2.0*TMath::Pi();
    if(IsolationFlag==1) {
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

    TVector3 diffpip = (*CDH_Pos)-(*CDH_Pos_pip);
    double diffPhinpip = (*CDH_Pos).Phi()-(*CDH_Pos_pip).Phi();
    const double difftofnpip = tofn - tofpip;
    if(diffPhinpip<-1.0*TMath::Pi()) diffPhinpip += 2.0*TMath::Pi();
    else if(diffPhinpip>1.0*TMath::Pi()) diffPhinpip -= 2.0*TMath::Pi();
    if(IsolationFlag==1 || IsolationFlag==2) {
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
    
    TVector3 diffpippim = (*CDH_Pos_pip)-(*CDH_Pos_pim);
    double diffPhipippim = (*CDH_Pos_pip).Phi()-(*CDH_Pos_pim).Phi();
    if(diffPhipippim<-1.0*TMath::Pi()) diffPhipippim +=2.0*TMath::Pi();
    else if(diffPhipippim>1.0*TMath::Pi()) diffPhipippim -= 2.0*TMath::Pi();
    diff2d_Phipippim_Phinpim->Fill(diffPhinpim,diffPhipippim);
    //diff2d_Zpippim_Znpim->Fill(diffpim.Z(),diffpippim.Z());

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
       pow(((nmiss_mass - anacuts::neutron_center)/5.0/anacuts::neutron_sigma),2.0) < 1.0){
       if(nmiss_mass < 1.05){
         SigmaPMissNViciFlag=true;
       }
    }
    
    //vicinity of Sigma+ &Nmiss events +extended area
    if( (nmiss_mass >=anacuts::neutron_center) && (pow(((MassNPip - anacuts::Sigmap_center)/5.0/anacuts::Sigmap_sigma),2.0) +
       pow(((nmiss_mass - anacuts::neutron_center)/5.0/anacuts::neutron_sigma),2.0) < 1.0)){
       if(nmiss_mass < 1.05){
         SigmaPMissNViciextFlag=true;
       }
    }
    
    if( (nmiss_mass < anacuts::neutron_center) && (fabs( MassNPip - anacuts::Sigmap_center) < 5.0*anacuts::Sigmap_sigma)){
      SigmaPMissNViciextFlag=true;
    }


    //vicinity of Sigma- & NMiss events
    if( (nmiss_mass >= anacuts::neutron_center) && 
        (pow(((MassNPim - anacuts::Sigmam_center)/5.0/anacuts::Sigmam_sigma),2.0) +
        pow(((nmiss_mass - anacuts::neutron_center)/3.0/anacuts::neutron_sigma),2.0) < 1.0)){
        SigmaMMissNViciFlag=true;
        SigmaMMissNViciextFlag=true;
    }
    if( (nmiss_mass < anacuts::neutron_center) && 
        (pow(((MassNPim - anacuts::Sigmam_center)/5.0/anacuts::Sigmam_sigma),2.0) +
        pow(((nmiss_mass - anacuts::neutron_center)/5.0/anacuts::neutron_sigma),2.0) < 1.0)){
        SigmaMMissNViciFlag=true;
    }

    
    if( (nmiss_mass < anacuts::neutron_center) && (fabs( MassNPim - anacuts::Sigmam_center) < 5.0*anacuts::Sigmam_sigma)){
      SigmaMMissNViciextFlag=true;
    }

    
    //vicinity of K0 & NMiss events
    if( pow(((nmiss_mass - anacuts::neutron_center)/5.0/anacuts::neutron_sigma),2.0) + 
        pow(((LVec_pip_pim.M() - anacuts::K0_center)/5.0/anacuts::K0_sigma),2.0) <1.0){
      K0MissNViciFlag=true;
    }



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
    if(sigmacuttype==0) {
      if(SigmaPFlag && !SigmaCrossFlagLeft && !SigmaCrossFlagRight) {
        SigmaPcutFlag[0]=true;
      }
      if(SigmaMFlag && !SigmaCrossFlagTop  && !SigmaCrossFlagBottom) {
        SigmaMcutFlag[0]=true;
      }
    }

    if(sigmacuttype==1) {
      if(SigmaPFlag && !SigmaMFlag) {
        SigmaPcutFlag[1]=true;
      }
      if(SigmaMFlag && !SigmaPFlag) {
        SigmaMcutFlag[1]=true;
      }
    }

    if(sigmacuttype==2) {
      if(SigmaPFlag && !SigmawideMFlag) {
        SigmaPcutFlag[2]=true;
      }
      if(SigmaMFlag && !SigmawidePFlag) {
        SigmaMcutFlag[2]=true;
      }
    }

    const double IMnpibinWidth = 1./nbinIMnpi;//2.5 MeV
    //Sigma cross side band Sp low mass top
    for(int igap=0; igap<ngap; igap++) {
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


    for(int izone=0; izone<nzone; izone++) {
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
    if(anacuts::neutron_MIN_wide<nmiss_mass && nmiss_mass<anacuts::neutron_MAX_wide ) MissNwideFlag=true;

    //K0 rejection using original momentum
    if( (LVec_pip_pim.M()<anacuts::pipi_MIN || anacuts::pipi_MAX<LVec_pip_pim.M())) K0rejectFlag=true;
    //if( (LVec_pip_pim.M()<anacuts::pipi_MIN_narrow || anacuts::pipi_MAX_narrow<LVec_pip_pim.M())) K0rejectFlag=true;
    if( (anacuts::pipi_MIN_narrow < LVec_pip_pim.M())  && (LVec_pip_pim.M() < anacuts::pipi_MAX_narrow)) K0Flag=true;

    bool IsBGregion = false;
    if(BGFlag_woSid_won==0){
      if(!SigmawidePFlag && !SigmawideMFlag && !MissNwideFlag){
        IsBGregion = true;
      }
    }
    //exclude pi+pi-lambda+n_true, pi+pi-pi0+lambda+n_true and so on
    //pick up n_fake events as much as possible
    else if(BGFlag_woSid_won==1){
      if(!SigmawidePFlag && !SigmawideMFlag && !MissNwideFlag ){
        IsBGregion = true;
      }
      if( anacuts::neutron_MAX_wide <=nmiss_mass) IsBGregion = false;
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
      if(SimK0nnmode){
        weight *=5.69470101983999943e-01;
        weight *=7.37894527999999994e-01;
      }
    }
    static bool isState = false;
    if(!isState){
      if(SimSpmode) std::cout << "Sim Sp mode " << std::endl;
      if(SimSmmode) std::cout << "Sim Sm mode " << std::endl;
      if(SimK0nnmode) std::cout << "Sim K0nn mode " << std::endl;
      std::cout  << " weighting factor " << weight << std::endl;
      isState = true;
    }

    int binnum = (LVec_pip_pim_n.M()-1.0)*nbintemplate ;
    int wbinnum=0;
    if(1.40<=LVec_pip_pim_n.M() && LVec_pip_pim_n.M()<1.52) wbinnum=1;
    else if(1.52<= LVec_pip_pim_n.M()) wbinnum=2;
    else wbinnum=0;

    if(IsMCweighting) {
      if(!SimFakemode_gSp && !SimFakemode_gSm){//mc for real data
        if(SimFakemode) { //w/o K0
          if(!isState){
            std::cout << "MC for real data" << std::endl;
            isState = true;
          }
          //weight *= fweight_nmom_v353->Eval((*LVec_n).P()); 
          //weight *= fweight_IMpippim_v364->Eval(LVec_pip_pim.M());
          //weight *= fweight_q_v366->Eval(qkn.P()); 
          //weight *= fweight_IMnpim_v367->Eval(LVec_pim_n.M());
          //weight *= fweight_MMnmiss_v368->Eval(nmiss_mass);
          //weight *= fweight_IMnpip_v369->Eval(LVec_pip_n.M());
        }else if(SimFakeK0mode) { //wK0
          //weight *= fweight_q_wK0_v377->Eval((qkn.P()));
          //weight *= fweight_MMnmiss_wK0_v378->Eval(nmiss_mass);
          //weight *= fweight_nmom_wK0_v379->Eval((*LVec_n).P());
          //weight *= fweight_IMnpip_wK0_v380->Eval(LVec_pip_n.M());
          //weight *= fweight_IMnpim_wK0_v381->Eval(LVec_pim_n.M());
        }
      }else if(SimFakemode_gSp){
        if(!isState){
          std::cout << "MC for geant Sp" << std::endl;
          isState = true;
        }
        if(SimFakemode) { //w/o K0
          //weight *= fweight_nmom_vSp23->Eval((*LVec_n).P()); 
          //weight *= fweight_IMpippim_vSp22->Eval(LVec_pip_pim.M());
          //weight *= fweight_q_vSp20->Eval(qkn.P()); 
          //weight *= fweight_MMnmiss_vSp21->Eval(nmiss_mass);
          //weight *= fweight_IMnpim_vSp25->Eval(LVec_pim_n.M());
          //weight *= fweight_IMnpip_vSp24->Eval(LVec_pip_n.M());
          //weight *= fweight_Mompippim_vSp26->Eval(LVec_pip_pim.P());
          //weight *= fweight_Momnpip_vSp27->Eval(LVec_pip_n.P());
          //weight *= fweight_Momnpim_vSp28->Eval(LVec_pim_n.P());
        }else if(SimFakeK0mode) { //wK0
          //weight = 0;
          //weight *= fweight_q_wK0_v377->Eval((qkn.P()));
          //weight *= fweight_MMnmiss_wK0_v378->Eval(nmiss_mass);
          //weight *= fweight_nmom_wK0_v379->Eval((*LVec_n).P());
          //weight *= fweight_IMnpip_wK0_v380->Eval(LVec_pip_n.M());
          //weight *= fweight_IMnpim_wK0_v381->Eval(LVec_pim_n.M());
        }
      }else if(SimFakemode_gSm){
        if(!isState){
          std::cout << "MC for geant Sm" << std::endl;
          isState = true;
        }
        if(SimFakemode) { //w/o K0
          //weight *= fweight_IMnpip_vSm11->Interpolate(LVec_pip_n.M());
          //weight *= fweight_IMnpim_vSm12->Interpolate(LVec_pim_n.M());
          //weight *= fweight_IMnpip_vSm16->Interpolate(LVec_pip_n.M());
          //weight *= fweight_IMnpim_vSm17->Interpolate(LVec_pim_n.M());
          //weight *= fweight_IMnpip_vSm20->Interpolate(LVec_pip_n.M());
          //weight *= fweight_IMnpim_vSm21->Interpolate(LVec_pim_n.M());
          
          //weight *= fweight_q_vSm23->Eval(qkn.P()); 
          //weight *= fweight_MMnmiss_vSm24->Eval(nmiss_mass);
          //weight *= fweight_nmom_vSm25->Eval((*LVec_n).P()); 
          //weight *= fweight_IMpippim_vSm26->Eval(LVec_pip_pim.M());
          //weight *= fweight_IMnpip_vSm27->Eval(LVec_pip_n.M());
          //weight *= fweight_IMnpim_vSm28->Eval(LVec_pim_n.M());
        }else if(SimFakeK0mode) { //wK0
          //weight = 0;
          //weight *= fweight_q_wK0_v377->Eval((qkn.P()));
          //weight *= fweight_MMnmiss_wK0_v378->Eval(nmiss_mass);
          //weight *= fweight_nmom_wK0_v379->Eval((*LVec_n).P());
          //weight *= fweight_IMnpip_wK0_v380->Eval(LVec_pip_n.M());
          //weight *= fweight_IMnpim_wK0_v381->Eval(LVec_pim_n.M());
        }
      }
    }//MCweighting
    //---end of Flag definition-----------------------------------------------------
    //w/o kinfit
    //---including K0 --------------------------------------------------------------

    //std::cout << __LINE__ << std::endl;
    NHitCDCOut->Fill(nhitOutCDC,weight);
    IsForwardCharge->Fill(ForwardCharge,weight);
    CDHphi_betainv_fid->Fill(1./NeutralBetaCDH,(*CDH_Pos).Phi());
    CDHz_betainv_fid->Fill(1./NeutralBetaCDH,(*CDH_Pos).z());
    dE_betainv_fid->Fill(1./NeutralBetaCDH,dE);
    //pipmom_IMpippim->Fill(LVec_pip_pim.M(),(*LVec_pip).P());
    //pimmom_IMpippim->Fill(LVec_pip_pim.M(),(*LVec_pim).P());
    //pipmom_IMnpip->Fill(LVec_pip_n.M(),(*LVec_pip).P());
    //pimmom_IMnpim->Fill(LVec_pim_n.M(),(*LVec_pim).P());
    if(NBetaOK) {
      dE_nmom_fid_beta->Fill((*LVec_n).P(),dE);
      //if(K0Flag) dE_nmom_fid_beta_wK0->Fill((*LVec_n).P(),dE);
      dE_MMom_fid_beta->Fill(LVec_npipimiss.P(),dE);
      dE_MMass_fid_beta->Fill(LVec_npipimiss.M(),dE);
      if(SigmaPFlag || SigmaMFlag) {
        dE_MMass_fid_beta_wSid->Fill(LVec_npipimiss.M(),dE);
      }
      dE_CDHphi->Fill((*CDH_Pos).Phi(),dE);
      dE_CDHz->Fill((*CDH_Pos).z(),dE);
      dE_IMnpim->Fill(LVec_pim_n.M(),dE);
      dE_IMnpip->Fill(LVec_pip_n.M(),dE);
    }
    if(NBetaOK && MissNFlag) {
      dE_IMnpim_n->Fill(LVec_pim_n.M(),dE,weight);
      dE_IMnpip_n->Fill(LVec_pip_n.M(),dE,weight);
    }
    if(NBetaOK && NdEOK) {

      CDHz_nmom_fid->Fill((*LVec_n).P(),(*CDH_Pos).z());
      MMom_MMass->Fill(LVec_npipimiss.M(),LVec_npipimiss.P(),weight);
      IMnpim_IMnpip_dE->Fill(LVec_pip_n.M(),LVec_pim_n.M(),weight);
      MMnmiss_IMnpip_dE->Fill(LVec_pip_n.M(),nmiss_mass,weight);
      MMnmiss_IMnpim_dE->Fill(LVec_pim_n.M(),nmiss_mass,weight);
      nmom_IMpippim->Fill(LVec_pip_pim.M(),(*LVec_n).P(),weight);
      nmom_MK0bar2->Fill(LVec_K0bar.M2(),(*LVec_n).P(),weight);
      MMnmiss_IMpippim_dE->Fill(LVec_pip_pim.M(),nmiss_mass,weight);
      if(mcpattern==2){
        MMnmiss_IMpippim_dE_pat2->Fill(LVec_pip_pim.M(),nmiss_mass,weight);
      }
      if(mcpattern==7){
        MMnmiss_IMpippim_dE_pat7->Fill(LVec_pip_pim.M(),nmiss_mass,weight);
      }
      MMnpim_MMnpip->Fill(LVec_pimmiss_nmiss.M(),LVec_pipmiss_nmiss.M());
      if(SimSpmode || SimSmmode || SimK0nnmode){
        if(IsFakebyVTX || IsFakeN2){
          MMnmiss_IMnpip_dE_fake->Fill(LVec_pip_n.M(),nmiss_mass,weight);
          MMnmiss_IMnpim_dE_fake->Fill(LVec_pim_n.M(),nmiss_mass,weight);
        }
      }
      if(SigmaPFlag || SigmaMFlag) {
        MMom_MMass_wSid->Fill(LVec_npipimiss.M(),LVec_npipimiss.P(),weight);
        MMnmiss_IMpippim_dE_wSid->Fill(LVec_pip_pim.M(),nmiss_mass,weight);
        MMnmiss_IMnpipi_wSid->Fill(LVec_pip_pim_n.M(), nmiss_mass,weight);
        //q_MMnmiss_wSid->Fill(nmiss_mass,qkn.P(),weight);
        nmom_MMnmiss_wSid->Fill(nmiss_mass,(*LVec_n).P(),weight);
        
        if(SimSpmode || SimSmmode || SimK0nnmode){
          //if(IsFakebyVTX || IsFakeN2){
          if(!( (mcpattern==2)  ||  (mcpattern==7))) {// || IsFakebyVTX ){
            MMnmiss_IMpippim_dE_wSid_fake->Fill(LVec_pip_pim.M(),nmiss_mass);
            nmom_MMnmiss_wSid_fake->Fill(nmiss_mass,(*LVec_n).P(),weight);
          }
        }
      }
      //pipmom_IMpippim_dE->Fill(LVec_pip_pim.M(),(*LVec_pip).P());
      //pimmom_IMpippim_dE->Fill(LVec_pip_pim.M(),(*LVec_pim).P());
      //pipmom_IMnpip_dE->Fill(LVec_pip_n.M(),(*LVec_pip).P());
      //pimmom_IMnpim_dE->Fill(LVec_pim_n.M(),(*LVec_pim).P());
      //pipmom_MMnmiss_dE->Fill(nmiss_mass,(*LVec_pip).P());
      //pimmom_MMnmiss_dE->Fill(nmiss_mass,(*LVec_pim).P());
      nmom_CDHphi->Fill((*CDH_Pos).Phi(),(*LVec_n).P());
      IMpippim_DCApipibeam->Fill(dca_pipibeam,LVec_pip_pim.M());
      IMnpip_DCApipibeam->Fill(dca_pipibeam,LVec_pip_n.M());
      IMnpim_DCApipibeam->Fill(dca_pipibeam,LVec_pim_n.M());
      diff2d_CDC_CDH_pim->Fill(diffPhinpim,diffpim.z());
      diff2d_CDC_CDH_pim_phi_tof->Fill(diffPhinpim,difftofnpim);
      diff2d_CDC_CDH_pim_z_tof->Fill(diffpim.z(),difftofnpim);
      diff2d_CDC_CDH_pip->Fill(diffPhinpip,diffpip.z());
      diff2d_CDC_CDH_pip_phi_tof->Fill(diffPhinpip,difftofnpip);
      diff2d_CDC_CDH_pip_z_tof->Fill(diffpip.z(),difftofnpip);
      dE_diffphi_CDC_CDH_pim->Fill(diffPhinpim,dE);
      dE_diffphi_CDC_CDH_pip->Fill(diffPhinpip,dE);
      //MMnmiss_dE->Fill(dE, nmiss_mass);
      MMnmiss_diffphi_CDC_CDH_pim->Fill(diffPhinpim, nmiss_mass);
      MMnmiss_diffphi_CDC_CDH_pip->Fill(diffPhinpip, nmiss_mass);
      pimmom_diffphi_CDC_CDH_pim->Fill(diffPhinpim,(*LVec_pim).P());
      pipmom_diffphi_CDC_CDH_pip->Fill(diffPhinpip,(*LVec_pip).P());
      pimmom_diffz_CDC_CDH_pim->Fill(diffpim.z(),(*LVec_pim).P());
      pipmom_diffz_CDC_CDH_pip->Fill(diffpip.z(),(*LVec_pip).P());
      nmom_diffphi_CDC_CDH_pim->Fill(diffPhinpim,(*LVec_n).P());
      nmom_diffphi_CDC_CDH_pip->Fill(diffPhinpip,(*LVec_n).P());
      nmom_diffz_CDC_CDH_pim->Fill(diffpim.z(),(*LVec_n).P());
      nmom_diffz_CDC_CDH_pip->Fill(diffpip.z(),(*LVec_n).P());
      //q_MMnmiss->Fill(nmiss_mass,qkn.P(),weight);
      
      if(SigmaMMissNViciFlag){
        IMpippim_IMnpim_vici->Fill(LVec_pim_n.M(),LVec_pip_pim.M(),weight);
        MMnmiss_IMpippim_dE_viciSm->Fill(LVec_pip_pim.M(),nmiss_mass,weight);
        if(!SigmawidePFlag){
          IMpippim_IMnpim_vici_woSp->Fill(LVec_pim_n.M(),LVec_pip_pim.M(),weight);
        }
      }
      if(SigmaPMissNViciFlag){
        IMpippim_IMnpip_vici->Fill(LVec_pip_n.M(),LVec_pip_pim.M(),weight);
        MMnmiss_IMpippim_dE_viciSp->Fill(LVec_pip_pim.M(),nmiss_mass,weight);
        if(!SigmawideMFlag){
          IMpippim_IMnpip_vici_woSm->Fill(LVec_pip_n.M(),LVec_pip_pim.M(),weight);
        }
      }

      if(K0MissNViciFlag){
        MMnmiss_IMpippim_dE_viciK0->Fill(LVec_pip_pim.M(),nmiss_mass,weight);
      }

      if(!SigmawidePFlag && !SigmawideMFlag ) { 
        IMnpim_IMnpip_dE_woSid->Fill(LVec_pip_n.M(),LVec_pim_n.M(),weight);
        MMnmiss_IMnpip_dE_woSid->Fill(LVec_pip_n.M(),nmiss_mass,weight);
        MMnmiss_IMnpim_dE_woSid->Fill(LVec_pim_n.M(),nmiss_mass,weight);
        MMnmiss_IMpippim_dE_woSid->Fill(LVec_pip_pim.M(),nmiss_mass,weight);
      }
      
      if(IsBGregion){
        //diff2d_CDC_CDH_pim_woSid_won->Fill(diffPhinpim,diffpim.z());
        diff2d_CDC_CDH_pim_phi_tof_woSid_won->Fill(diffPhinpim,difftofnpim);
        diff2d_CDC_CDH_pim_z_tof_woSid_won->Fill(diffpim.z(),difftofnpim);
        //diff2d_CDC_CDH_pip_woSid_won->Fill(diffPhinpip,diffpip.z());
        diff2d_CDC_CDH_pip_phi_tof_woSid_won->Fill(diffPhinpip,difftofnpip);
        diff2d_CDC_CDH_pip_z_tof_woSid_won->Fill(diffpip.z(),difftofnpip);
        //diff2d_Phipippim_Phinpim_woSid_won->Fill(diffPhinpim,diffPhipippim);
        diff2d_Zpippim_Znpim_woSid_won->Fill(diffpim.Z(),diffpippim.Z());

        //IMnpim_IMnpip_dE_woSid_won->Fill(LVec_pip_n.M(),LVec_pim_n.M(),weight);
        //MMnmiss_IMpippim_dE_woSid_won->Fill(LVec_pip_pim.M(),nmiss_mass,weight);
        //q_nmom_woSid_won->Fill((*LVec_n).P(),qkn.P(),weight);
        //MMnmiss_IMnpip_dE_woSid_won->Fill(LVec_pip_n.M(),nmiss_mass,weight);
        //MMnmiss_IMnpim_dE_woSid_won->Fill(LVec_pim_n.M(),nmiss_mass,weight);
        //Momnpim_Momnpip_dE_woSid_won->Fill(LVec_pip_n.P(),LVec_pim_n.P(),weight);
        //Momnpim_Mompippim_dE_woSid_won->Fill(LVec_pip_pim.P(),LVec_pim_n.P(),weight);
        pimmom_diffphi_CDC_CDH_pim_woSid_won->Fill(diffPhinpim,(*LVec_pim).P());
        pipmom_diffphi_CDC_CDH_pip_woSid_won->Fill(diffPhinpip,(*LVec_pip).P());
      }


      if(K0Flag) {
        IMnpim_IMnpip_dE_wK0->Fill(LVec_pip_n.M(),LVec_pim_n.M(),weight);
        nmom_cosn_wK0->Fill(cos_ncdslab,(*LVec_n).P(),weight);
        nmom_cosK0_wK0->Fill(cos_K0CM,(*LVec_n).P(),weight);
        nmom_cosK0n_wK0->Fill(cos_K0_ncds_CM,(*LVec_n).P(),weight);
        K0mom_cosK0_wK0->Fill(cos_K0CM,LVec_pip_pim.P(),weight);
        nmissmom_cosnmiss_wK0->Fill(cos_nmissCM,LVec_npipimiss.P(),weight);
        nmom_K0mom->Fill(LVec_pip_pim.P(),(*LVec_n).P(),weight);
        nmom_nmissmom_wK0->Fill(LVec_npipimiss.P(),(*LVec_n).P(),weight);
        nmissmom_K0mom->Fill(LVec_pip_pim.P(),LVec_npipimiss.P(),weight);
        MMnmiss_DCApipibeam_wK0->Fill(dca_pipibeam,nmiss_mass,weight);
        nmom_MMnmiss_wK0->Fill(nmiss_mass,(*LVec_n).P(),weight);
        diff2d_CDC_CDH_pim_wK0->Fill(diffPhinpim,diffpim.z(),weight);
        diff2d_CDC_CDH_pip_wK0->Fill(diffPhinpip,diffpip.z(),weight);
        if( !(SigmawidePFlag && MissNwideFlag) && !(SigmawideMFlag && MissNwideFlag)) {
          if(!MissNwideFlag) {
            MMnmiss_IMnpip_dE_wK0_woSidn_won->Fill(LVec_pip_n.M(),nmiss_mass,weight);
            MMnmiss_IMnpim_dE_wK0_woSidn_won->Fill(LVec_pim_n.M(),nmiss_mass,weight);
            MMnmiss_IMpippim_dE_wK0_woSidn_won->Fill(LVec_pip_pim.M(),nmiss_mass,weight);
            IMnpim_IMnpip_dE_wK0_woSidn_won->Fill(LVec_pip_n.M(),LVec_pim_n.M(),weight);
            MMnmiss_IMnpipi_wK0_woSidn_won->Fill(LVec_pip_pim_n.M(),nmiss_mass,weight);
            q_IMnpipi_wK0_woSidn_won->Fill(LVec_pip_pim_n.M(),qkn.P(),weight);
            nmom_MMnmiss_wK0_woSidn_won->Fill(nmiss_mass,(*LVec_n).P(),weight);
            MMnmiss_Mompippim_dE_wK0_woSidn_won->Fill(LVec_pip_pim.P(),nmiss_mass,weight);
            MMom_MMass_wK0_woSidn_won->Fill(LVec_npipimiss.M(),LVec_npipimiss.P(),weight);
            //pipmom_MMnmiss_dE_wK0_woSidn_won->Fill(nmiss_mass,(*LVec_pip).P(),weight);
            //pimmom_MMnmiss_dE_wK0_woSidn_won->Fill(nmiss_mass,(*LVec_pim).P(),weight);
          }
        }
        //17th,Mar. 2020.
        //applied narrow window to remove \pi^{\pm}\Sigma^{\mp} final states from K0nn events
        //if(!SigmawidePFlag && !SigmawideMFlag && !K0rejectFlag_narrow)  
        if(!SigmawidePFlag && !SigmawideMFlag ) { 
          nmom_MMnmiss_wK0_woSid->Fill(nmiss_mass,(*LVec_n).P(),weight);
        }
        if(IsBGregion) {
          nmom_IMnpipi_wK0_woSid_won->Fill(LVec_pip_pim_n.M(),(*LVec_n).P(),weight);
          MMnmiss_IMnpip_dE_wK0_woSid_won->Fill(LVec_pip_n.M(),nmiss_mass,weight);
          MMnmiss_Momnpip_dE_wK0_woSid_won->Fill(LVec_pip_n.P(),nmiss_mass,weight);
          MMnmiss_IMnpim_dE_wK0_woSid_won->Fill(LVec_pim_n.M(),nmiss_mass,weight);
          MMnmiss_Momnpim_dE_wK0_woSid_won->Fill(LVec_pim_n.P(),nmiss_mass,weight);
          Momnpim_Momnpip_dE_wK0_woSid_won->Fill(LVec_pip_n.P(),LVec_pim_n.P(),weight);
          Momnpim_Mompippim_dE_wK0_woSid_won->Fill(LVec_pip_pim.P(),LVec_pim_n.P(),weight);
          Momnpip_Mompippim_dE_wK0_woSid_won->Fill(LVec_pip_pim.P(),LVec_pip_n.P(),weight);
          MMnmiss_IMpippim_dE_wK0_woSid_won->Fill(LVec_pip_pim.M(),nmiss_mass,weight);
          IMnpim_IMnpip_dE_wK0_woSid_won->Fill( LVec_pip_n.M(), LVec_pim_n.M(),weight);
          IMpippim_IMnpip_wK0_woSid_won->Fill(LVec_pip_n.M(),LVec_pip_pim.M(),weight);
          IMpippim_IMnpim_wK0_woSid_won->Fill(LVec_pim_n.M(),LVec_pip_pim.M(),weight);
          MMnmiss_IMnpipi_wK0_woSid_won->Fill(LVec_pip_pim_n.M(),nmiss_mass,weight);
          MMnmiss_Momnpipi_wK0_woSid_won->Fill(LVec_pip_pim_n.P(),nmiss_mass,weight);
          q_IMpippim_wK0_woSid_won->Fill(LVec_pip_pim.M(),qkn.P(),weight);
          q_IMnpipi_wK0_woSid_won->Fill(LVec_pip_pim_n.M(),qkn.P(),weight);
          q_nmom_wK0_woSid_won->Fill((*LVec_n).P(),qkn.P(),weight);
          q_IMnpip_wK0_woSid_won->Fill(LVec_pip_n.M(),qkn.P(),weight);
          q_IMnpim_wK0_woSid_won->Fill(LVec_pim_n.M(),qkn.P(),weight);
          //q_MMnmiss_wK0_woSid_won->Fill(nmiss_mass,qkn.P(),weight);
          nmom_MMnmiss_wK0_woSid_won->Fill(nmiss_mass,(*LVec_n).P(),weight);
          MMnmiss_Mompippim_dE_wK0_woSid_won->Fill(LVec_pip_pim.P(),nmiss_mass,weight);
          MMom_MMass_wK0_woSid_won->Fill(LVec_npipimiss.M(),LVec_npipimiss.P(),weight);
          //pipmom_MMnmiss_dE_wK0_woSid_won->Fill(nmiss_mass,(*LVec_pip).P(),weight);
          //pipmom_pimmom_dE_wK0_woSid_won->Fill((*LVec_pim).P(),(*LVec_pip).P(),weight);
          //pimmom_MMnmiss_dE_wK0_woSid_won->Fill(nmiss_mass,(*LVec_pim).P(),weight);
          nmom_cosn_wK0_woSid_won->Fill(cos_ncdslab,(*LVec_n).P(),weight);
          nmom_cospip_wK0_woSid_won->Fill(cos_pip,(*LVec_n).P(),weight);
          nmom_cospippim_wK0_woSid_won->Fill(cos_pippim,(*LVec_n).P(),weight);
          nmom_coslabpim_wK0_woSid_won->Fill(cos_pim,(*LVec_n).P(),weight);
          nmom_cospim_wK0_woSid_won->Fill(cos_pim,(*LVec_n).P(),weight);
          nmom_phinpip_wK0_woSid_won->Fill(phi_npip,(*LVec_n).P(),weight);
          nmom_phinpim_wK0_woSid_won->Fill(phi_npim,(*LVec_n).P(),weight);
          nmom_phipip_wK0_woSid_won->Fill(phi_pip,(*LVec_n).P(),weight);
          nmom_phipim_wK0_woSid_won->Fill(phi_pim,(*LVec_n).P(),weight);
          nmom_phin_wK0_woSid_won->Fill(phi_n,(*LVec_n).P(),weight);
          nmom_pipmom_wK0_woSid_won->Fill((*LVec_pip).P(),(*LVec_n).P(),weight);
          nmom_pimmom_wK0_woSid_won->Fill((*LVec_pim).P(),(*LVec_n).P(),weight);
          nmom_IMnpip_dE_wK0_woSid_won->Fill(LVec_pip_n.M(),(*LVec_n).P(),weight);
          nmom_IMnpim_dE_wK0_woSid_won->Fill(LVec_pim_n.M(),(*LVec_n).P(),weight);
          nmom_IMpippim_wK0_woSid_won->Fill(LVec_pip_pim.M(),(*LVec_n).P(),weight);
        }//BG region
        if(SigmaPFlag || SigmaMFlag) {
          nmom_nmissmom_wK0_wSid->Fill(LVec_npipimiss.P(),(*LVec_n).P(),weight);
          MMnmiss_IMnpipi_wK0_wSid->Fill(LVec_pip_pim_n.M(), nmiss_mass,weight);
          MMnmiss_IMpippim_dE_wK0_wSid->Fill(LVec_pip_pim.M(),nmiss_mass,weight);
        }
        if(!SigmawideMFlag) {
          MMnmiss_IMnpip_dE_wK0_woSm->Fill(LVec_pip_n.M(),nmiss_mass,weight);
        }
        if(!SigmawidePFlag) {
          MMnmiss_IMnpim_dE_wK0_woSp->Fill(LVec_pim_n.M(),nmiss_mass,weight);
        }
      }//wK0
    }



    //std::cout << __LINE__ << std::endl;
    if(NBetaOK && NdEOK && MissNFlag) {
      IMnpim_IMnpip_dE_n->Fill(LVec_pip_n.M(),LVec_pim_n.M(),weight);
      IMnpip_IMnpipi_n->Fill(LVec_pip_pim_n.M(),LVec_pip_n.M(),weight);
      IMnpim_IMnpipi_n->Fill(LVec_pip_pim_n.M(),LVec_pim_n.M(),weight);
      
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
        }else{
          IMnpim_IMnpip_dE_n_fake->Fill(LVec_pip_n.M(),LVec_pim_n.M());
        }
      }

      //pipmom_IMpippim_dE_n->Fill(LVec_pip_pim.M(),(*LVec_pip).P(),weight);
      //pimmom_IMpippim_dE_n->Fill(LVec_pip_pim.M(),(*LVec_pim).P(),weight);
      //pipmom_IMnpip_dE_n->Fill(LVec_pip_n.M(),(*LVec_pip).P(),weight);
      //pimmom_IMnpim_dE_n->Fill(LVec_pim_n.M(),(*LVec_pim).P(),weight);

      for(int igap=0; igap<ngap; igap++) {
        if(SigmaPsideFlag[sidebandtype][igap] || SigmaMsideFlag[sidebandtype][igap]) {
          //IMnpim_IMnpip_dE_n_side[igap]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
        }
        if(SigmaPsideLowFlag[igap] && SigmaPsideFlag[sidebandtype][igap]) {
          //IMnpim_IMnpip_dE_n_Sp_side[LOWside][igap]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
        }
        if(SigmaPsideHighFlag[igap] && SigmaPsideFlag[sidebandtype][igap]) {
          //IMnpim_IMnpip_dE_n_Sp_side[HIGHside][igap]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
        }
        if(SigmaMsideLowFlag[igap] && SigmaMsideFlag[sidebandtype][igap]) {
         // IMnpim_IMnpip_dE_n_Sm_side[LOWside][igap]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
        }
        if(SigmaMsideHighFlag[igap] && SigmaMsideFlag[sidebandtype][igap]) {
         // IMnpim_IMnpip_dE_n_Sm_side[HIGHside][igap]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
        }
      }
      
      /*
      for(int izone=0; izone<nzone; izone++) {
        if(sidebandtype==0) {
          if( !SigmaCrossPsideWideFlagLeft[izone] && !SigmaCrossPsideWideFlagRight[izone]  && SigmaPsideWideFlag[izone]) {
            //IMnpim_IMnpip_dE_n_Sp_sidewide[izone]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
            //q_IMnpipi_wSid_n_Sp_sidewide[izone]->Fill(LVec_pip_pim_n.M(),qkn.P());
          }
          if( !SigmaCrossMsideWideFlagTop[izone] && !SigmaCrossMsideWideFlagBottom[izone]  && SigmaMsideWideFlag[izone]) {
            //IMnpim_IMnpip_dE_n_Sm_sidewide[izone]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
            //q_IMnpipi_wSid_n_Sm_sidewide[izone]->Fill(LVec_pip_pim_n.M(),qkn.P());
          }
        } else if(sidebandtype==1) {
          if( !SigmaMFlag  && SigmaPsideWideFlag[izone]) {
            IMnpim_IMnpip_dE_n_Sp_sidewide[izone]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
            q_IMnpipi_wSid_n_Sp_sidewide[izone]->Fill(LVec_pip_pim_n.M(),qkn.P());
          }
          if( !SigmaPFlag  && SigmaMsideWideFlag[izone]) {
            IMnpim_IMnpip_dE_n_Sm_sidewide[izone]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
            q_IMnpipi_wSid_n_Sm_sidewide[izone]->Fill(LVec_pip_pim_n.M(),qkn.P());
          }
        } else if(sidebandtype==2) {
          if( !SigmawideMFlag  && SigmaPsideWideFlag[izone]) {
            IMnpim_IMnpip_dE_n_Sp_sidewide[izone]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
            q_IMnpipi_wSid_n_Sp_sidewide[izone]->Fill(LVec_pip_pim_n.M(),qkn.P());
          }
          if( !SigmawidePFlag  && SigmaMsideWideFlag[izone]) {
            IMnpim_IMnpip_dE_n_Sm_sidewide[izone]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
            q_IMnpipi_wSid_n_Sm_sidewide[izone]->Fill(LVec_pip_pim_n.M(),qkn.P());
          }
        }
      }//izone
      */

      MMnpim_MMnpip_n->Fill(LVec_pimmiss_nmiss.M(),LVec_pipmiss_nmiss.M(),weight);
      //MMnpim_MMnpip_n->Fill(LVec_pip_nmiss.M(),LVec_pim_nmiss.M());
      if(!SigmaPFlag  && !SigmaMFlag){
        //MMnpim_MMnpip_woSid_n->Fill(LVec_pip_nmiss.M(),LVec_pim_nmiss.M());
        MMnpim_MMnpip_woSid_n->Fill(LVec_pimmiss_nmiss.M(),LVec_pipmiss_nmiss.M(),weight);
        if(K0rejectFlag){
          MMnpim_MMnpip_woK0_woSid_n->Fill(LVec_pimmiss_nmiss.M(),LVec_pipmiss_nmiss.M(),weight);
          IMnpim_IMnpip_dE_woK0_woSid_n->Fill(LVec_pip_n.M(),LVec_pim_n.M(),weight);
          q_IMnpipi_woK0_woSid_n->Fill(LVec_pip_pim_n.M(),qkn.P(),weight);
        }
      }


      //0: diagonal cut
      //1: 3 sigma cut
      //2: 5 simga cut
      nmom->Fill((*LVec_n).P());
      mnmom->Fill(nmiss_mom);
      dE_nmom->Fill((*LVec_n).P(),dE);
      npipmom->Fill(LVec_pip_n.P());
      npimmom->Fill(LVec_pim_n.P());
      nmom_IMnpim_dE_n->Fill(LVec_pim_n.M(),(*LVec_n).P(),weight);
      nmom_IMnpip_dE_n->Fill(LVec_pip_n.M(),(*LVec_n).P(),weight);
      if(SigmaPFlag) {
        nmom_IMnpip_dE_n_Sp->Fill(LVec_pip_n.M(),(*LVec_n).P(),weight);
        q_IMnpipi_wSid_n_Sp->Fill(LVec_pip_pim_n.M(),qkn.P(),weight);
        if(LVec_pip_n.P()<0.03){
          q_IMnpipi_wSid_n_Sp_Stop->Fill(LVec_pip_pim_n.M(),qkn.P(),weight);
        }else{
          q_IMnpipi_wSid_n_Sp_NoStop->Fill(LVec_pip_pim_n.M(),qkn.P(),weight);
        }
        q_IMpiSigma_wSid_n_Sp_genacc->Fill(LVec_piSigma_react.M()/1000.,qkn_react.P()/1000.);
        q_IMnpipi_wSid_n_Sp_acc->Fill(LVec_pip_pim_n_mc.M(),qkn_mc.P());
        q_IMnpipi_wSid_n_Sp_acc_reco->Fill(LVec_pip_pim_n.M(),qkn.P());
        IMnpim_IMnpip_dE_n_Sp->Fill(LVec_pip_n.M(),LVec_pim_n.M(),weight);
        if(!SigmawideMFlag){
          q_IMnpipi_wSid_n_Sp_woSm->Fill(LVec_pip_pim_n.M(),qkn.P(),weight);
        }
        if(SimSpmode) {
          q_IMnpipi_wSid_n_Sp_mc->Fill(LVec_pip_pim_n_mc.M(),qkn_mc.P());
          diff_IMnpipi_wSid_n_Sp->Fill(LVec_pip_pim_n.M(),LVec_pip_pim_n.M()-LVec_pip_pim_n_mc.M(),weight);
          diff_q_wSid_n_Sp->Fill(qkn.P(),qkn.P()-qkn_mc.P(),weight);
        }
      }

      for(int itype=0; itype<3; itype++) {
        for(int igap=0; igap<ngap; igap++) {
          if(SigmaPsideFlag[itype][igap] && SigmaPsideLowFlag[igap] ) {
            //q_IMnpipi_wSid_n_Sp_side[itype][LOWside][igap]->Fill(LVec_pip_pim_n.M(),qkn.P());
          }
          if(SigmaPsideFlag[itype][igap] && SigmaPsideHighFlag[igap] ) {
            //q_IMnpipi_wSid_n_Sp_side[itype][HIGHside][igap]->Fill(LVec_pip_pim_n.M(),qkn.P());
          }
        }
      }

      if(SigmaMFlag) {
        nmom_IMnpim_dE_n_Sm->Fill(LVec_pim_n.M(),(*LVec_n).P(),weight);
        q_IMnpipi_wSid_n_Sm->Fill(LVec_pip_pim_n.M(),qkn.P(),weight);
        if(LVec_pim_n.P()<0.03){
          q_IMnpipi_wSid_n_Sm_Stop->Fill(LVec_pip_pim_n.M(),qkn.P(),weight);
        }else{
          q_IMnpipi_wSid_n_Sm_NoStop->Fill(LVec_pip_pim_n.M(),qkn.P(),weight);
        }
        q_IMpiSigma_wSid_n_Sm_genacc->Fill(LVec_piSigma_react.M()/1000.,qkn_react.P()/1000.);
        q_IMnpipi_wSid_n_Sm_acc->Fill(LVec_pip_pim_n_mc.M(),qkn_mc.P());
        q_IMnpipi_wSid_n_Sm_acc_reco->Fill(LVec_pip_pim_n.M(),qkn.P());
        IMnpim_IMnpip_dE_n_Sm->Fill(LVec_pip_n.M(),LVec_pim_n.M(),weight);
        if(!SigmawidePFlag){
          q_IMnpipi_wSid_n_Sm_woSp->Fill(LVec_pip_pim_n.M(),qkn.P(),weight);
        }
        
        if(SimSmmode) {
          q_IMnpipi_wSid_n_Sm_mc->Fill(LVec_pip_pim_n_mc.M(),qkn_mc.P());
          diff_IMnpipi_wSid_n_Sm->Fill(LVec_pip_pim_n.M(),LVec_pip_pim_n.M()-LVec_pip_pim_n_mc.M(),weight);
          diff_q_wSid_n_Sm->Fill(qkn.P(),qkn.P()-qkn_mc.P(),weight);
        }
      }
      
      if(SigmaPFlag && SigmaMFlag){
        q_IMnpipi_wSid_n_SpSm->Fill(LVec_pip_pim_n.M(),qkn.P(),weight);
      }


      for(int igap=0; igap<ngap; igap++) {
        for(int itype=0; itype<3; itype++) {
          if(SigmaMsideFlag[itype][igap] && SigmaMsideLowFlag[igap]) {
            //q_IMnpipi_wSid_n_Sm_side[itype][LOWside][igap]->Fill(LVec_pip_pim_n.M(),qkn.P());
          }
          if(SigmaMsideFlag[itype][igap] && SigmaMsideHighFlag[igap]) {
            //q_IMnpipi_wSid_n_Sm_side[itype][HIGHside][igap]->Fill(LVec_pip_pim_n.M(),qkn.P());
          }
        }
        if(NBetaOK && NdEOK && MissNFlag) {
          if(SigmaPsideFlag[sidebandtype][igap] || SigmaMsideFlag[sidebandtype][igap]) {
            //q_IMnpipi_wSid_n_side[igap]->Fill(LVec_pip_pim_n.M(),qkn.P());
          }
        }
      }

      q_IMpippim_n->Fill(LVec_pip_pim.M(),qkn.P(),weight);
      IMpippim_IMnpipi_n->Fill(LVec_pip_pim_n.M(), LVec_pip_pim.M(),weight);
      IMpippim_IMnpip_n->Fill(LVec_pip_n.M(),LVec_pip_pim.M(),weight);
      //IMpippim_IMnpip_n_bin[binnum]->Fill(LVec_pip_n.M(),LVec_pip_pim.M(),weight);
      IMpippim_IMnpim_n->Fill(LVec_pim_n.M(),LVec_pip_pim.M(),weight);
      //IMpippim_IMnpim_n_bin[binnum]->Fill(LVec_pim_n.M(),LVec_pip_pim.M(),weight);
      if(!SigmawideMFlag){
        IMpippim_IMnpip_n_woSm->Fill(LVec_pip_n.M(),LVec_pip_pim.M(),weight);
        //IMpippim_IMnpip_n_woSm_bin[binnum]->Fill(LVec_pip_n.M(),LVec_pip_pim.M(),weight);
      }
      if(!SigmawidePFlag){
        IMpippim_IMnpim_n_woSp->Fill(LVec_pim_n.M(),LVec_pip_pim.M(),weight);
        //IMpippim_IMnpim_n_woSp_bin[binnum]->Fill(LVec_pim_n.M(),LVec_pip_pim.M(),weight);
      }
      
      if(!SigmaMcutFlag[sigmacuttype]) {
        IMpippim_IMnpip_n_woSmdia->Fill(LVec_pip_n.M(),LVec_pip_pim.M(),weight);
        //IMpippim_IMnpip_n_woSmdia_bin[binnum]->Fill(LVec_pip_n.M(),LVec_pip_pim.M(),weight);
      }
      if(!SigmaPcutFlag[sigmacuttype])  {
        IMpippim_IMnpim_n_woSpdia->Fill(LVec_pim_n.M(),LVec_pip_pim.M(),weight);
        //IMpippim_IMnpim_n_woSpdia_bin[binnum]->Fill(LVec_pim_n.M(),LVec_pip_pim.M(),weight);
      }

      
      q_IMnpip_n->Fill(LVec_pip_n.M(),qkn.P(),weight);
      q_IMnpim_n->Fill(LVec_pim_n.M(),qkn.P(),weight);
      //q_nmom_n->Fill((*LVec_n).P(),qkn.P(),weight);
      if(K0Flag || SigmaPFlag) {
        IMpippim_IMnpip_n_cross->Fill(LVec_pip_n.M(),LVec_pip_pim.M(),weight);
      }
      if(K0Flag || SigmaMFlag) {
        IMpippim_IMnpim_n_cross->Fill(LVec_pim_n.M(),LVec_pip_pim.M(),weight);
      }

      nmom_IMpippim_n->Fill(LVec_pip_pim.M(),(*LVec_n).P(),weight);
      mnmom_IMpippim_n->Fill(LVec_pip_pim.M(),(LVec_npipimiss).P(),weight);
      IMpippim_DCApipibeam_n->Fill(dca_pipibeam,LVec_pip_pim.M(),weight);
      IMnpip_DCApipibeam_n->Fill(dca_pipibeam,LVec_pip_n.M(),weight);
      IMnpim_DCApipibeam_n->Fill(dca_pipibeam,LVec_pim_n.M(),weight);
      if(SigmaPFlag || SigmaMFlag) {
        //MMom_MMass_wSid_n->Fill(LVec_npipimiss.M(),LVec_npipimiss.P(),weight);
        //MMnpim_MMnpip_wSid_n->Fill(LVec_pip_nmiss.M(),LVec_pim_nmiss.M(),weight);
        MMnpim_MMnpip_wSid_n->Fill(LVec_pimmiss_nmiss.M(),LVec_pipmiss_nmiss.M());
        MMnmiss_IMpippim_dE_wSid_n->Fill(LVec_pip_pim.M(),nmiss_mass ,weight);
        Momnpim_Momnpip_dE_wSid_n->Fill(LVec_pip_n.P(),LVec_pim_n.P(),weight);
        Momnpim_Mompippim_dE_wSid_n->Fill(LVec_pip_pim.P(),LVec_pim_n.P(),weight);
        Momnpip_Mompippim_dE_wSid_n->Fill(LVec_pip_pim.P(),LVec_pip_n.P(),weight);
        q_IMpippim_wSid_n->Fill(LVec_pip_pim.M(),qkn.P(),weight);
        IMnpim_IMnpip_dE_wSid_n->Fill(LVec_pip_n.M(),LVec_pim_n.M(),weight);
        //IMnpim_IMnpip_dE_wSid_n_bin[binnum]->Fill(LVec_pip_n.M(),LVec_pim_n.M(),weight);
        Cosn_IMnpipi_wSid_n->Fill(LVec_pip_pim_n.M(),cos_nmissCM,weight);
        //dE_IMnpipi_wSid_n->Fill(LVec_pip_pim_n.M(),dE,weight);
        IMpippim_IMnpipi_n_wSid->Fill(LVec_pip_pim_n.M(), LVec_pip_pim.M(),weight);
        nmom_IMpippim_wSid_n->Fill(LVec_pip_pim.M(),(*LVec_n).P(),weight);
        nmom_cosn_wSid_n->Fill(cos_ncdslab,(*LVec_n).P(),weight);
        nmom_cosnmiss_wSid_n->Fill(cos_nmissCM,(*LVec_n).P(),weight);
        mnmom_IMpippim_wSid_n->Fill(LVec_pip_pim.M(),(LVec_npipimiss).P(),weight);
        q_IMnpipi_wSid_n->Fill(LVec_pip_pim_n.M(),qkn.P(),weight);
        pipmom_Momnpip_wSid_n->Fill(LVec_pip_n.P(),(*LVec_pip).P() ,weight);
        pimmom_Momnpim_wSid_n->Fill(LVec_pip_n.P(),(*LVec_pim).P(),weight);

        for(int icut=0;icut<nthetacut;icut++){
          if(nmissthetalab< TMath::Pi()/nthetacut*(icut+1)) {
            q_IMnpipi_wSid_n_thetacut[icut]->Fill(LVec_pip_pim_n.M(),qkn.P(),weight);
          }           
        }
        q_nmom_wSid_n->Fill((*LVec_n).P(), qkn.P(),weight);
        q_IMpiSigma_wSid_n_genacc->Fill(LVec_piSigma_react.M()/1000.,qkn_react.P()/1000.,weight);
        q_IMnpipi_wSid_n_acc->Fill(LVec_pip_pim_n_mc.M(),qkn_mc.P(),weight);
        q_IMnpipi_wSid_n_acc_reco->Fill(LVec_pip_pim_n.M(),qkn.P(),weight);
        q_IMnpip_n_wSid->Fill(LVec_pip_n.M(),qkn.P(),weight);
        q_IMnpim_n_wSid->Fill(LVec_pim_n.M(),qkn.P(),weight);
        //q_MMnmiss_n_wSid->Fill(nmiss_mass,qkn.P(),weight);
        //q_nmom_n_wSid->Fill((*LVec_n).P(),qkn.P(),weight);
        nmom_IMnpipi_wSid_n->Fill(LVec_pip_pim_n.M(),(*LVec_n).P(),weight);
        nmom_MMnmiss_wSid_n->Fill(nmiss_mass,(*LVec_n).P(),weight);
        pipmom_IMnpipi_wSid_n->Fill(LVec_pip_pim_n.M(),(*LVec_pip).P(),weight);
        pimmom_IMnpipi_wSid_n->Fill(LVec_pip_pim_n.M(),(*LVec_pim).P(),weight);
        if(SigmaPFlag){
          IMpippim_IMnpip_wSid_n_Sp->Fill(LVec_pip_n.M(),LVec_pip_pim.M(),weight);
          IMpippim_IMnpip_wSid_n_Sp_bin[binnum]->Fill(LVec_pip_n.M(),LVec_pip_pim.M(),weight);
          IMpippim_IMnpip_wSid_n_Sp_wbin[wbinnum]->Fill(LVec_pip_n.M(),LVec_pip_pim.M(),weight);
          IMnpim_IMnpip_dE_wSid_n_Sp->Fill(LVec_pip_n.M(),LVec_pim_n.M(),weight);
          IMnpim_IMnpip_dE_wSid_n_Sp_bin[binnum]->Fill(LVec_pip_n.M(),LVec_pim_n.M(),weight);
          IMnpim_IMnpip_dE_wSid_n_Sp_wbin[wbinnum]->Fill(LVec_pip_n.M(),LVec_pim_n.M(),weight);
        }

        if(SigmaMFlag){
          IMpippim_IMnpim_wSid_n_Sm->Fill(LVec_pim_n.M(),LVec_pip_pim.M(),weight);
          IMpippim_IMnpim_wSid_n_Sm_bin[binnum]->Fill(LVec_pim_n.M(),LVec_pip_pim.M(),weight);
          IMpippim_IMnpim_wSid_n_Sm_wbin[wbinnum]->Fill(LVec_pim_n.M(),LVec_pip_pim.M(),weight);
          IMnpim_IMnpip_dE_wSid_n_Sm->Fill(LVec_pip_n.M(),LVec_pim_n.M(),weight);
          IMnpim_IMnpip_dE_wSid_n_Sm_bin[binnum]->Fill(LVec_pip_n.M(),LVec_pim_n.M(),weight);
          IMnpim_IMnpip_dE_wSid_n_Sm_wbin[wbinnum]->Fill(LVec_pip_n.M(),LVec_pim_n.M(),weight);
        }

        if(!SigmawideMFlag){
          IMpippim_IMnpip_wSid_n_woSm->Fill(LVec_pip_n.M(),LVec_pip_pim.M(),weight);
          //IMpippim_IMnpip_wSid_n_woSm_bin[binnum]->Fill(LVec_pip_n.M(),LVec_pip_pim.M(),weight);
          IMpippim_IMnpipi_n_wSid_woSm->Fill(LVec_pip_pim_n.M(), LVec_pip_pim.M(),weight);
        }
        if(!SigmawidePFlag){
          IMpippim_IMnpim_wSid_n_woSp->Fill(LVec_pim_n.M(),LVec_pip_pim.M(),weight);
          //IMpippim_IMnpim_wSid_n_woSp_bin[binnum]->Fill(LVec_pim_n.M(),LVec_pip_pim.M(),weight);
          IMpippim_IMnpipi_n_wSid_woSp->Fill(LVec_pip_pim_n.M(), LVec_pip_pim.M(),weight);
        }
        
        if(!(SigmawidePFlag && SigmawideMFlag)){
          q_IMnpipi_wSid_n_wocross->Fill(LVec_pip_pim_n.M(),qkn.P(),weight);
        }

        //reaction data - mcData matching
        bool IsMissMassNOK = false;
        bool IsMcNMassOK = false;
        if(SimSpmode || SimSmmode || SimK0nnmode){
          if(!IsFakebyVTX)MMnpim_MMnpip_wSid_n_mc->Fill(LVec_pip_nmiss_mc.M(),LVec_pim_nmiss_mc.M(),weight);
          double diffIMnpim_reactmc = LVec_pim_n.M()- LVec_pim_n_mc.M();
          double diffIMnpip_reactmc = LVec_pip_n.M()- LVec_pip_n_mc.M();
          double diffMMnmiss_recomc = nmiss_mass - (*mcmom_nmiss).M();
          
          //MCdata missing mass check
          if(fabs((*mcmom_nmiss).M()-nMass)<0.01) IsMissMassNOK = true;
          //MCdata ncds check
          if(fabs((*mcmom_ncds).M()-nMass)<0.01) IsMcNMassOK = true;
          
          
          if(IsMissMassNOK && IsMcNMassOK){
            diff2D_IMnpim_Momnpim_wSid_n_reactmc->Fill(LVec_pim_n_mc.P()-LVec_Sigma_react.P()/1000.,LVec_pim_n_mc.M()-LVec_Sigma_react.M()/1000.);
            diff2D_IMnpip_Momnpip_wSid_n_reactmc->Fill(LVec_pip_n_mc.P()-LVec_Sigma_react.P()/1000.,LVec_pip_n_mc.M()-LVec_Sigma_react.M()/1000.);
            diff_nmiss_wSid_n_reactmc->Fill((*mcmom_nmiss).P()-(*react_nmiss).P()/1000.);
            diff2D_MMnmiss_IMnpim_wSid_n_reactmc->Fill(LVec_pim_n_mc.M()-LVec_Sigma_react.M()/1000.,(*mcmom_nmiss).M()-(*react_nmiss).M()/1000.);
            diff2D_MMnmiss_IMnpip_wSid_n_reactmc->Fill(LVec_pip_n_mc.M()-LVec_Sigma_react.M()/1000.,(*mcmom_nmiss).M()-(*react_nmiss).M()/1000.);
            //Sigma mass check and missing neurtom mom. check. 
            //note: Reaction data does not have the info. of the neutron from Sigma
            if(SimSpmode){
              if(fabs(LVec_pip_n_mc.M()-LVec_Sigma_react.M()/1000.) > 0.002) IsFakeN1 = true;
              if(fabs((*mcmom_nmiss).P()-(*react_nmiss).P()/1000.) > 0.002) IsFakeN1 = true;
            }
            if(SimSmmode){
              if(fabs(LVec_pim_n_mc.M()-LVec_Sigma_react.M()/1000.) > 0.002) IsFakeN1 = true;
              if(fabs((*mcmom_nmiss).P()-(*react_nmiss).P()/1000.) > 0.002) IsFakeN1 = true;
            }
            //added angle check->why angle ? ->mom. check. 
            //if(!IsFakeN1){
              diff2D_IMnpim_Momnpim_wSid_n_reactmc_fake1->Fill(LVec_pim_n_mc.P()-LVec_Sigma_react.P()/1000.,LVec_pim_n_mc.M()-LVec_Sigma_react.M()/1000.);
              diff2D_IMnpip_Momnpip_wSid_n_reactmc_fake1->Fill(LVec_pip_n_mc.P()-LVec_Sigma_react.P()/1000.,LVec_pip_n_mc.M()-LVec_Sigma_react.M()/1000.);
              diff_cosnmiss_wSid_n_reactmc->Fill((*mcmom_nmiss).CosTheta()-(*react_nmiss).CosTheta());
              if(fabs((*mcmom_nmiss).CosTheta()-(*react_nmiss).CosTheta())>0.002) IsFakeN1 = true;
            //}
          }
          
          //std::cout << __LINE__ << std::endl;

          double diffIMnpim_recomc = LVec_pim_n.M()- LVec_pim_n_mc.M();
          double diffIMnpip_recomc = LVec_pip_n.M()- LVec_pip_n_mc.M();
          TLorentzVector diffnpip_recomc = LVec_pip_n - LVec_pip_n_mc; 
          TLorentzVector diffnpim_recomc = LVec_pim_n - LVec_pim_n_mc; 
          TLorentzVector diffMMom_recomc = LVec_npipimiss - *mcmom_nmiss;
          double diffnmom_recomc = (*LVec_n).P() - (*mcmom_ncds).P();
          double diffnpip_mcreact = LVec_pip_n_mc.P() - LVec_Sigma_react.P()/1000.0;
          double diffnpim_mcreact = LVec_pim_n_mc.P() - LVec_Sigma_react.P()/1000.0;
          double diffMassnpip_mcreact = LVec_pip_n_mc.M() - LVec_Sigma_react.M()/1000.0;
          double diffMassnpim_mcreact = LVec_pim_n_mc.M() - LVec_Sigma_react.M()/1000.0;
          diff2D_MMnmiss_IMnpim_recomc_wSid_n->Fill(diffIMnpim_recomc,diffMMnmiss_recomc);
          diff2D_MMnmiss_IMnpip_recomc_wSid_n->Fill(diffIMnpip_recomc,diffMMnmiss_recomc);
          diff2D_nmom_IMnpim_recomc_wSid_n->Fill(diffIMnpim_recomc,diffnmom_recomc);
          diff2D_nmom_IMnpip_recomc_wSid_n->Fill(diffIMnpip_recomc,diffnmom_recomc);
          diffMomnpim_Momnpip_recomc_wSid_n->Fill(diffnpip_recomc.P(),diffnpim_recomc.P());
          diffMMom_recomc_wSid_n->Fill(diffMMom_recomc.P());
          vtxr_generation_ncan_wSid_n_mc->Fill(mcncdsgen,mcncanvtxr);
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
          //if(!IsFakeN1){
          if(!IsFakebyVTX){
            diff2D_MMnmiss_IMnpim_recomc_wSid_n_fake1->Fill(diffIMnpim_recomc,diffMMnmiss_recomc);
            diff2D_MMnmiss_IMnpip_recomc_wSid_n_fake1->Fill(diffIMnpip_recomc,diffMMnmiss_recomc);
            if(IsMcNMassOK){
              diff2D_nmom_IMnpim_recomc_wSid_n_fake1->Fill(diffIMnpim_recomc,diffnmom_recomc);
              diff2D_nmom_IMnpip_recomc_wSid_n_fake1->Fill(diffIMnpip_recomc,diffnmom_recomc);
            }
          }
          //if(SimSpmode){
          //  if( (diffIMnpip_recomc<-0.012) || (0.010<diffIMnpip_recomc)) IsFakeN2=true;
            //if(diffnpip_recomc.P()>0.10) IsFakeN2 = true;
          //}
          //if(SimSmmode){
          //  if( (diffIMnpim_recomc<-0.012) || (0.012<diffIMnpim_recomc)) IsFakeN2=true;
            //if(diffnpim_recomc.P()>0.10) IsFakeN2 = true;
          //}
          //if(IsFakeN1 || !IsMissMassNOK || !IsMcNMassOK || IsFakeN2){
          //if(!IsMissMassNOK || !IsMcNMassOK || IsFakeN2){
          //if( (mcncdsgen!=3) || (mcncanvtxr>20)  ){
          //if( (mcncanvtxr>58)  && (mcncdsgen!=3) ){
          //if(IsFakebyVTX  || IsFakeN2){
          //if( !((mcpattern==2))) { //|| (mcpattern==7))) { // || IsFakebyVTX ){
          //if( IsFakebyVTX ){
            //if(  (SimSpmode &&  ( ( LVec_pip_nmiss_mc.M() < 1.189) || (1.190 < LVec_pip_nmiss_mc.M())) )
            //   ||(SimSmmode &&  ( ( LVec_pim_nmiss_mc.M() < 1.197) || (1.198 < LVec_pim_nmiss_mc.M())) 
          if( mcpattern==2 ) {
            MMnmiss_IMpippim_dE_wSid_n_pat2->Fill(LVec_pip_pim.M(),nmiss_mass);
          }
          if( mcpattern==7 ) {
            MMnmiss_IMpippim_dE_wSid_n_pat7->Fill(LVec_pip_pim.M(),nmiss_mass);
          }
          if(!( (mcpattern==2)  ||  (mcpattern==7))) {// || IsFakebyVTX ){
            q_IMnpipi_wSid_n_fake->Fill(LVec_pip_pim_n.M(),qkn.P());
            MMnmiss_IMpippim_dE_wSid_n_fake->Fill(LVec_pip_pim.M(),nmiss_mass);
            IMnpim_IMnpip_dE_wSid_n_fake->Fill(LVec_pip_n.M(),LVec_pim_n.M());

            IMnpim_MMnpim_mc_wSid_n_fake->Fill(LVec_pim_nmiss_mc.M(),LVec_pim_n_mc.M());
            IMnpip_MMnpip_mc_wSid_n_fake->Fill(LVec_pip_nmiss_mc.M(),LVec_pip_n_mc.M());
          }
        }//SimSpmode or SimSmmode
      }//wSid_n

      for(int igap=0; igap<ngap; igap++) {
        if(SigmaPsideFlag[sidebandtype][igap] || SigmaMsideFlag[sidebandtype][igap]) {
          //q_IMnpipi_wSid_n_side[igap]->Fill(LVec_pip_pim_n.M(),qkn.P());
        }
      }

      if( K0Flag || (SigmaPFlag || SigmaMFlag)){
          IMpippim_IMnpip_wK0orwSid_n->Fill(LVec_pip_n.M(),LVec_pip_pim.M(),weight);
          IMpippim_IMnpip_wK0orwSid_n_bin[binnum]->Fill(LVec_pip_n.M(),LVec_pip_pim.M(),weight);
          IMpippim_IMnpip_wK0orwSid_n_wbin[wbinnum]->Fill(LVec_pip_n.M(),LVec_pip_pim.M(),weight);
          IMpippim_IMnpim_wK0orwSid_n->Fill(LVec_pim_n.M(),LVec_pip_pim.M(),weight);
          IMpippim_IMnpim_wK0orwSid_n_bin[binnum]->Fill(LVec_pim_n.M(),LVec_pip_pim.M(),weight);
          IMpippim_IMnpim_wK0orwSid_n_wbin[wbinnum]->Fill(LVec_pim_n.M(),LVec_pip_pim.M(),weight);
          IMnpim_IMnpip_dE_wK0orwSid_n->Fill(LVec_pip_n.M(),LVec_pim_n.M(),weight);
          IMnpim_IMnpip_dE_wK0orwSid_n_bin[binnum]->Fill(LVec_pip_n.M(),LVec_pim_n.M(),weight);
          IMnpim_IMnpip_dE_wK0orwSid_n_wbin[wbinnum]->Fill(LVec_pip_n.M(),LVec_pim_n.M(),weight);
          q_IMnpipi_wK0orwSid_n->Fill(LVec_pip_pim_n.M(),qkn.P(),weight);
          if(!SigmawidePFlag){
            IMnpim_IMnpip_dE_wK0orwSid_n_woSp->Fill(LVec_pip_n.M(),LVec_pim_n.M(),weight);
          }
          if(!SigmawideMFlag){
            IMnpim_IMnpip_dE_wK0orwSid_n_woSm->Fill(LVec_pip_n.M(),LVec_pim_n.M(),weight);
          }

          if(!SigmaMcutFlag[sigmacuttype]) {
            IMpippim_IMnpip_wK0orwSid_n_woSmdia->Fill(LVec_pip_n.M(),LVec_pip_pim.M(),weight);
            //IMpippim_IMnpip_wK0orwSid_n_woSmdia_bin[binnum]->Fill(LVec_pip_n.M(),LVec_pip_pim.M(),weight);
          }
          if(!SigmaPcutFlag[sigmacuttype]) {
            IMpippim_IMnpim_wK0orwSid_n_woSpdia->Fill(LVec_pim_n.M(),LVec_pip_pim.M(),weight);
            //IMpippim_IMnpim_wK0orwSid_n_woSpdia_bin[binnum]->Fill(LVec_pim_n.M(),LVec_pip_pim.M(),weight);
          }
      }
        }//NBetaOK && NdEOK && MissNFlag
    //---including K0 end----------------------------------------------------------------------

    //std::cout << __LINE__ << std::endl;
    //selection K0
    if(K0Flag && NBetaOK && NdEOK && MissNFlag) {
      //nmom_cosn_wK0_n->Fill(cos_ncdslab,(*LVec_n).P(),weight);
      if(cos_nmissCM>0.95 ) {
        nmom_cosn_wK0_n_forward->Fill(cos_ncdslab,(*LVec_n).P(),weight);
      }
      IMnpim_IMnpip_dE_wK0_n->Fill(LVec_pip_n.M(),LVec_pim_n.M(),weight);
      nmom_cosK0_wK0_n->Fill(cos_K0CM,(*LVec_n).P(),weight);
      nmom_cosK0n_wK0_n->Fill(cos_K0_ncds_CM,(*LVec_n).P(),weight);
      nmom_cosnmiss_wK0_n->Fill(cos_nmissCM,(*LVec_n).P(),weight);
      nmom_cosnnmiss_wK0_n->Fill(cos_nmiss_ncds_CM,(*LVec_n).P(),weight);
      K0mom_cosK0_wK0_n->Fill(cos_K0CM,LVec_pip_pim.P(),weight);
      nmissmom_cosnmiss_wK0_n->Fill(cos_nmissCM,LVec_npipimiss.P(),weight);
      nmissmom_cosK0nmiss_wK0_n->Fill(cos_nmiss_pippim_CM,LVec_npipimiss.P(),weight);
      nmom_K0mom_n->Fill(LVec_pip_pim.P(),(*LVec_n).P(),weight);
      nmom_nmissmom_wK0_n->Fill(LVec_npipimiss.P(),(*LVec_n).P(),weight);
      nmissmom_K0mom_n->Fill(LVec_pip_pim.P(),LVec_npipimiss.P(),weight);
      nmom_cosnlab_K0_n->Fill(cos_ncdslab,(*LVec_n).P(),weight);
      nmom_IMnpipi_wK0_n->Fill(LVec_pip_pim_n.M(),(*LVec_n).P(),weight);
      IMpippim_DCApipibeam_wK0_n->Fill(dca_pipibeam,LVec_pip_pim.M(),weight);

      IMpippim_IMnpip_wK0_n->Fill(LVec_pip_n.M(),LVec_pip_pim.M(),weight);
      IMpippim_IMnpip_wK0_n_bin[binnum]->Fill(LVec_pip_n.M(),LVec_pip_pim.M(),weight);
      IMpippim_IMnpip_wK0_n_wbin[wbinnum]->Fill(LVec_pip_n.M(),LVec_pip_pim.M(),weight);
      IMpippim_IMnpim_wK0_n->Fill(LVec_pim_n.M(),LVec_pip_pim.M(),weight);
      IMpippim_IMnpim_wK0_n_wbin[wbinnum]->Fill(LVec_pim_n.M(),LVec_pip_pim.M(),weight);

      if(!SigmaMcutFlag[sigmacuttype]) {
        IMpippim_IMnpip_wK0_n_woSmdia->Fill(LVec_pip_n.M(),LVec_pip_pim.M(),weight);
        //IMpippim_IMnpip_wK0_n_woSmdia_bin[binnum]->Fill(LVec_pip_n.M(),LVec_pip_pim.M(),weight);
      }
      if(!SigmaPcutFlag[sigmacuttype]) {
        IMpippim_IMnpim_wK0_n_woSpdia->Fill(LVec_pim_n.M(),LVec_pip_pim.M(),weight);
        //IMpippim_IMnpim_wK0_n_woSpdia_bin[binnum]->Fill(LVec_pim_n.M(),LVec_pip_pim.M(),weight);
      }



      if(SigmaPFlag || SigmaMFlag) {
        q_nmom_wK0_wSid_n->Fill((*LVec_n).P(), qkn.P(),weight);
        IMnpim_IMnpip_dE_wK0_wSid_n->Fill(LVec_pip_n.M(),LVec_pim_n.M(),weight);
        nmom_IMnpipi_wK0_wSid_n->Fill(LVec_pip_pim_n.M(),(*LVec_n).P(),weight);
        nmom_nmissmom_wK0_wSid_n->Fill(LVec_npipimiss.P(),(*LVec_n).P(),weight);
        MMnmiss_IMpippim_dE_wK0_wSid_n->Fill(LVec_pip_pim.M(),nmiss_mass,weight);
        Momnpim_Momnpip_dE_wK0_wSid_n->Fill(LVec_pip_n.P(),LVec_pim_n.P(),weight);
        Momnpim_Mompippim_dE_wK0_wSid_n->Fill(LVec_pip_pim.P(),LVec_pim_n.P(),weight);
        Momnpip_Mompippim_dE_wK0_wSid_n->Fill(LVec_pip_pim.P(),LVec_pip_n.P(),weight);
          
        if(SigmaPFlag){
          IMnpim_IMnpip_dE_wK0_wSid_n_Sp->Fill(LVec_pip_n.M(),LVec_pim_n.M(),weight);
          IMnpim_IMnpip_dE_wK0_wSid_n_Sp_wbin[wbinnum]->Fill(LVec_pip_n.M(),LVec_pim_n.M(),weight);
        }
        
        if(SigmaMFlag){
          IMnpim_IMnpip_dE_wK0_wSid_n_Sm->Fill(LVec_pip_n.M(),LVec_pim_n.M(),weight);
          IMnpim_IMnpip_dE_wK0_wSid_n_Sm_wbin[wbinnum]->Fill(LVec_pip_n.M(),LVec_pim_n.M(),weight);
        }


        if(SimSpmode || SimSmmode || SimK0nnmode){
          double diffIMnpim_recomc = LVec_pim_n.M()- LVec_pim_n_mc.M();
          double diffIMnpip_recomc = LVec_pip_n.M()- LVec_pip_n_mc.M();
          double diffMMnmiss_recomc = nmiss_mass - (*mcmom_nmiss).M();
          double diffnmom_recomc = (*LVec_n).P() - (*mcmom_ncds).P();
          diff2D_MMnmiss_IMnpim_recomc_wK0_wSid_n->Fill(diffIMnpim_recomc,diffMMnmiss_recomc);
          diff2D_MMnmiss_IMnpip_recomc_wK0_wSid_n->Fill(diffIMnpip_recomc,diffMMnmiss_recomc);
          diff2D_nmom_IMnpim_recomc_wK0_wSid_n->Fill(diffIMnpim_recomc,diffnmom_recomc);
          diff2D_nmom_IMnpip_recomc_wK0_wSid_n->Fill(diffIMnpip_recomc,diffnmom_recomc);
          if(!IsFakeN1){
            diff2D_MMnmiss_IMnpim_recomc_wK0_wSid_n_fake1->Fill(diffIMnpim_recomc,diffMMnmiss_recomc);
            diff2D_MMnmiss_IMnpip_recomc_wK0_wSid_n_fake1->Fill(diffIMnpip_recomc,diffMMnmiss_recomc);
            diff2D_nmom_IMnpim_recomc_wK0_wSid_n_fake1->Fill(diffIMnpim_recomc,diffnmom_recomc);
            diff2D_nmom_IMnpip_recomc_wK0_wSid_n_fake1->Fill(diffIMnpip_recomc,diffnmom_recomc);
          }
        }
      }

      double xx = 1./sqrt(2.0)*(LVec_pip_n.M()-LVec_pim_n.M());
      double yy = 1./sqrt(2.0)*(LVec_pip_n.M()+LVec_pim_n.M());
      //double yy2 = yy-(cosh(1.96*xx)-1.0);
      double yy2 = yy-(sqrt(6.76*xx*xx+2.725)-1.0);
      IMnpim_IMnpip_dE_wK0_n_45rot3->Fill(xx,yy2,weight);
      
      if(!SigmawidePFlag && !SigmawideMFlag) {
        IMpippim_DCApipibeam_wK0_woSid_n->Fill(dca_pipibeam,LVec_pip_pim.M(),weight);
        nmom_IMnpipi_wK0_woSid_n->Fill(LVec_pip_pim_n.M(),(*LVec_n).P(),weight);
        q_IMnpipi_wK0_woSid_n->Fill(LVec_pip_pim_n.M(),qkn.P(),weight);
        IMnpim_IMnpip_dE_wK0_woSid_n->Fill(LVec_pip_n.M(),LVec_pim_n.M(),weight);
        IMnpim_IMnpip_dE_wK0_woSid_n_45rot->Fill(1./sqrt(2.0)*(LVec_pip_n.M()+LVec_pim_n.M()),
                                                 1./sqrt(2.0)*(LVec_pip_n.M()-LVec_pim_n.M()),
                                                 weight);
        IMnpim_IMnpip_dE_wK0_woSid_n_45rot2->Fill(xx,
                                                  yy,
                                                 weight);
        IMnpim_IMnpip_dE_wK0_woSid_n_45rot3->Fill(xx,yy2,weight);

        IMnpim_IMnpip_dE_wK0_woSid_n_wbin[wbinnum]->Fill(LVec_pip_n.M(),LVec_pim_n.M(),weight);
        nmom_K0mom_woSid_n->Fill(LVec_pip_pim.P(),(*LVec_n).P(),weight);
        q_nmom_wK0_woSid_n_wbin[wbinnum]->Fill((*LVec_n).P(),qkn.P(),weight);
        IMpippim_IMnpip_wK0_woSid_n_wbin[wbinnum]->Fill(LVec_pip_n.M(),LVec_pip_pim.M(),weight);
        IMpippim_IMnpim_wK0_woSid_n_wbin[wbinnum]->Fill(LVec_pim_n.M(),LVec_pip_pim.M(),weight);
        MMnmiss_IMpippim_dE_wK0_woSid_n_wbin[wbinnum]->Fill(LVec_pip_pim.M(),nmiss_mass,weight);
      }

      Mompippim_IMnpipi_dE_wK0_n->Fill(LVec_pip_pim_n.M(),LVec_pip_pim.P(),weight);
      Mompippim_nmom_dE_wK0_n->Fill((*LVec_n).P(),LVec_pip_pim.P(),weight);
      q_IMnpipi_wK0_n->Fill(LVec_pip_pim_n.M(),qkn.P(),weight);
      
      if(SigmaPFlag){
        q_IMnpipi_wK0_wSid_n_Sp->Fill(LVec_pip_pim_n.M(),qkn.P(),weight);
      }

      if(SigmaMFlag){
        q_IMnpipi_wK0_wSid_n_Sm->Fill(LVec_pip_pim_n.M(),qkn.P(),weight);
      }

      if(SigmaPFlag || SigmaMFlag) {
        q_IMnpipi_wK0_wSid_n->Fill(LVec_pip_pim_n.M(),qkn.P(),weight);
      }
      
      if(SigmaPFlag && SigmaMFlag) {
        q_IMnpipi_wK0_wSid_n_SpSm->Fill(LVec_pip_pim_n.M(),qkn.P(),weight);
        IMnpim_IMnpip_dE_wK0_wSid_n_SpSm->Fill(LVec_pip_n.M(),LVec_pim_n.M(),weight);
      }

      if(SigmaPcutFlag[sigmacuttype]) {
        Mompippim_IMnpipi_dE_wK0_n_Sp->Fill(LVec_pip_pim_n.M(),LVec_pip_pim.P(),weight);
      }
      
      IMnpip_IMnpipi_wK0_n->Fill(LVec_pip_pim_n.M(),LVec_pip_n.M(),weight);
      IMnpim_IMnpipi_wK0_n->Fill(LVec_pip_pim_n.M(),LVec_pim_n.M(),weight);

      if(SigmaMcutFlag[sigmacuttype]) {
        Mompippim_IMnpipi_dE_wK0_n_Sm->Fill(LVec_pip_pim_n.M(),LVec_pip_pim.P(),weight);
      }
      if(!SigmawideMFlag) {
        MMnmiss_IMnpip_dE_wK0_woSm_n->Fill(LVec_pip_n.M(),nmiss_mass,weight);
      }
      if(!SigmawidePFlag) {
        MMnmiss_IMnpim_dE_wK0_woSp_n->Fill(LVec_pim_n.M(),nmiss_mass,weight);
      }
    }

    //std::cout << __LINE__ << std::endl;
    //---rejecting K0--------------------------------------------------------------------------
    //K0 rejection
    if(K0rejectFlag && NBetaOK) {
      //dE_MMom_fid_beta_woK0->Fill(LVec_npipimiss.P(),dE,weight);
      //dE_MMass_fid_beta_woK0->Fill(LVec_npipimiss.M(),dE,weight);
      if(SigmaPFlag || SigmaMFlag) {
        //dE_MMass_fid_beta_woK0_wSid->Fill(LVec_npipimiss.M(),dE,weight);
      }
      //dE_IMnpim_woK0->Fill(LVec_pim_n.M(),dE,weight);
      //dE_IMnpip_woK0->Fill(LVec_pip_n.M(),dE,weight);
    }
    if(K0rejectFlag && NBetaOK && MissNFlag) {
      //dE_IMnpim_woK0_n->Fill(LVec_pim_n.M(),dE,weight);
      //dE_IMnpip_woK0_n->Fill(LVec_pip_n.M(),dE,weight);
    }
    if(K0rejectFlag && NBetaOK && NdEOK) {
      MMom_MMass_woK0->Fill(LVec_npipimiss.M(),LVec_npipimiss.P(),weight);
      IMnpim_IMnpip_dE_woK0->Fill(LVec_pip_n.M(),LVec_pim_n.M(),weight);
      MMnmiss_IMnpip_dE_woK0->Fill(LVec_pip_n.M(),nmiss_mass,weight);
      MMnmiss_IMnpim_dE_woK0->Fill(LVec_pim_n.M(),nmiss_mass,weight);
      if(!SigmawidePFlag) {
        MMnmiss_IMnpim_dE_woK0_woSp->Fill(LVec_pim_n.M(),nmiss_mass,weight);
        if(MissNFlag && SigmaMFlag) {
          MMnmiss_IMnpim_dE_woK0_woSp_cross->Fill(LVec_pim_n.M(),nmiss_mass,weight);
        }
        if(SigmaMMissNViciFlag){
          MMnmiss_IMnpim_dE_woK0_woSp_vici->Fill(LVec_pim_n.M(),nmiss_mass,weight);
          IMnpim_IMnpip_dE_woK0_woSp_vici->Fill(LVec_pip_n.M(),LVec_pim_n.M(),weight);
        }
        if(SigmaMMissNViciextFlag){
          MMnmiss_IMnpim_dE_woK0_woSp_viciext->Fill(LVec_pim_n.M(),nmiss_mass,weight);
          IMnpim_IMnpip_dE_woK0_woSp_viciext->Fill(LVec_pip_n.M(),LVec_pim_n.M(),weight);
        }
      }
      if(!SigmawideMFlag) {
        MMnmiss_IMnpip_dE_woK0_woSm->Fill(LVec_pip_n.M(),nmiss_mass,weight);
        if(MissNFlag && SigmaPFlag) {
          MMnmiss_IMnpip_dE_woK0_woSm_cross->Fill(LVec_pip_n.M(),nmiss_mass,weight);
        }
        if(SigmaPMissNViciFlag){
          MMnmiss_IMnpip_dE_woK0_woSm_vici->Fill(LVec_pip_n.M(),nmiss_mass,weight);
          IMnpim_IMnpip_dE_woK0_woSm_vici->Fill(LVec_pip_n.M(),LVec_pim_n.M(),weight);
        }
        if(SigmaPMissNViciextFlag){
          MMnmiss_IMnpip_dE_woK0_woSm_viciext->Fill(LVec_pip_n.M(),nmiss_mass,weight);
          IMnpim_IMnpip_dE_woK0_woSm_viciext->Fill(LVec_pip_n.M(),LVec_pim_n.M(),weight);
        }
      }
      //std::cout << __LINE__ << std::endl;
      /*
      if( !(SigmawidePFlag && MissNwideFlag) && !(SigmawideMFlag && MissNwideFlag)) {
        MMnmiss_IMnpip_dE_woK0_woSidn->Fill(LVec_pip_n.M(),nmiss_mass,weight);
        MMnmiss_IMnpim_dE_woK0_woSidn->Fill(LVec_pim_n.M(),nmiss_mass,weight);
        MMnmiss_IMpippim_dE_woK0_woSidn->Fill(LVec_pip_pim.M(),nmiss_mass,weight);
        IMnpim_IMnpip_dE_woK0_woSidn->Fill(LVec_pip_n.M(),LVec_pim_n.M(),weight);
        MMnmiss_IMnpipi_woK0_woSidn->Fill(LVec_pip_pim_n.M(),nmiss_mass,weight);
        q_IMnpipi_woK0_woSidn->Fill(LVec_pip_pim_n.M(),qkn.P(),weight);
        nmom_MMnmiss_woK0_woSidn->Fill(nmiss_mass,(*LVec_n).P(),weight);
        MMnmiss_Mompippim_dE_woK0_woSidn->Fill(LVec_pip_pim.P(),nmiss_mass,weight);
        MMom_MMass_woK0_woSidn->Fill(LVec_npipimiss.M(),LVec_npipimiss.P(),weight);
        pipmom_MMnmiss_dE_woK0_woSidn->Fill(nmiss_mass,(*LVec_pip).P(),weight);
        pipmom_pimmom_dE_woK0_woSidn->Fill((*LVec_pim).P(), (*LVec_pip).P(),weight);
        pimmom_MMnmiss_dE_woK0_woSidn->Fill(nmiss_mass,(*LVec_pim).P(),weight);
        nmom_cosn_woK0_woSidn->Fill(cos_ncdslab,(*LVec_n).P(),weight);
        nmom_cospip_woK0_woSidn->Fill(cos_pip,(*LVec_n).P(),weight);
        nmom_cospim_woK0_woSidn->Fill(cos_pim,(*LVec_n).P(),weight);
        nmom_phinpip_woK0_woSidn->Fill(phi_npip,(*LVec_n).P(),weight);
        nmom_phinpim_woK0_woSidn->Fill(phi_npim,(*LVec_n).P(),weight);
        nmom_phipip_woK0_woSidn->Fill(phi_pip,(*LVec_n).P(),weight);
        nmom_phipim_woK0_woSidn->Fill(phi_pim,(*LVec_n).P(),weight);
        nmom_phin_woK0_woSidn->Fill(phi_n,(*LVec_n).P(),weight);
        nmom_pipmom_woK0_woSidn->Fill((*LVec_pip).P(),(*LVec_n).P(),weight);
        nmom_pimmom_woK0_woSidn->Fill((*LVec_pim).P(),(*LVec_n).P(),weight);
        if(MissNFlag && (SigmaPFlag || SigmaMFlag)) {
          MMnmiss_IMnpip_dE_woK0_woSidn_cross->Fill(LVec_pip_n.M(),nmiss_mass,weight);
          MMnmiss_IMnpim_dE_woK0_woSidn_cross->Fill(LVec_pim_n.M(),nmiss_mass,weight);
          q_IMnpipi_woK0_woSidn_cross->Fill(LVec_pip_pim_n.M(),qkn.P(),weight);
          if(SigmaPcutFlag[sigmacuttype]) {
            MMnmiss_IMnpip_dE_woK0_woSidn_cross_Sp->Fill(LVec_pip_n.M(),nmiss_mass,weight);
            q_IMnpipi_woK0_woSidn_cross_Sp->Fill(LVec_pip_pim_n.M(),qkn.P(),weight);
          }
          if(SigmaMcutFlag[sigmacuttype]) {
            MMnmiss_IMnpim_dE_woK0_woSidn_cross_Sm->Fill(LVec_pim_n.M(),nmiss_mass,weight);
            q_IMnpipi_woK0_woSidn_cross_Sm->Fill(LVec_pip_pim_n.M(),qkn.P(),weight);
          }
        }
      }*/

      MMnmiss_IMpippim_dE_woK0->Fill(LVec_pip_pim.M(),nmiss_mass,weight);
      if(SigmaPFlag || SigmaMFlag) {
        MMom_MMass_woK0_wSid->Fill(LVec_npipimiss.M(),LVec_npipimiss.P(),weight);
        MMnmiss_IMnpipi_woK0_wSid->Fill(LVec_pip_pim_n.M(), nmiss_mass,weight);
        MMnmiss_IMpippim_dE_woK0_wSid->Fill(LVec_pip_pim.M(),nmiss_mass,weight);
        MMnmiss_DCApipibeam_woK0_wSid->Fill(dca_pipibeam,nmiss_mass,weight);
        nmom_MMnmiss_woK0_wSid->Fill(nmiss_mass,(*LVec_n).P(),weight);
        nmom_nmissmom_woK0_wSid->Fill(LVec_npipimiss.P(),(*LVec_n).P(),weight);
        diff2d_CDC_CDH_pim_woK0_wSid->Fill(diffPhinpim,diffpim.z(),weight);
        diff2d_CDC_CDH_pip_woK0_wSid->Fill(diffPhinpip,diffpip.z(),weight);
        MMnmiss_diffphi_CDC_CDH_pim_woK0_wSid->Fill(diffPhinpim, nmiss_mass,weight);
        MMnmiss_diffphi_CDC_CDH_pip_woK0_wSid->Fill(diffPhinpip, nmiss_mass,weight);
        MMnmiss_diffz_CDC_CDH_pim_woK0_wSid->Fill(diffpim.z(), nmiss_mass,weight);
        MMnmiss_diffz_CDC_CDH_pip_woK0_wSid->Fill(diffpip.z(), nmiss_mass,weight);
        pimmom_diffphi_CDC_CDH_pim_woK0_wSid->Fill(diffPhinpim,(*LVec_pim).P(),weight);
        pipmom_diffphi_CDC_CDH_pip_woK0_wSid->Fill(diffPhinpip,(*LVec_pip).P(),weight);
        nmom_diffphi_CDC_CDH_pim_woK0_wSid->Fill(diffPhinpim,(*LVec_n).P(),weight);
        nmom_diffphi_CDC_CDH_pip_woK0_wSid->Fill(diffPhinpip,(*LVec_n).P(),weight);
        nmom_diffz_CDC_CDH_pim_woK0_wSid->Fill(diffpim.z(),(*LVec_n).P(),weight);
        nmom_diffz_CDC_CDH_pip_woK0_wSid->Fill(diffpip.z(),(*LVec_n).P(),weight);
      }
      //std::cout << __LINE__ << std::endl;
      if(!SigmawidePFlag && !SigmawideMFlag) {
        nmom_MMnmiss_woK0_woSid->Fill(nmiss_mass,(*LVec_n).P(),weight);
        MMnmiss_IMnpip_dE_woK0_woSid->Fill(LVec_pip_n.M(),nmiss_mass,weight);
        MMnmiss_IMnpim_dE_woK0_woSid->Fill(LVec_pim_n.M(),nmiss_mass,weight);
        MMnmiss_IMpippim_dE_woK0_woSid->Fill(LVec_pip_pim.M(),nmiss_mass,weight);
        IMnpim_IMnpip_dE_woK0_woSid->Fill(LVec_pip_n.M(),LVec_pim_n.M(),weight);
        q_IMnpipi_woK0_woSid->Fill(LVec_pip_pim_n.M(),qkn.P(),weight);
      }
      if(IsBGregion) {
        nmom_IMnpipi_woK0_woSid_won->Fill(LVec_pip_pim_n.M(),(*LVec_n).P(),weight);
        MMnmiss_IMnpip_dE_woK0_woSid_won->Fill(LVec_pip_n.M(),nmiss_mass,weight);
        MMnmiss_Momnpip_dE_woK0_woSid_won->Fill(LVec_pip_n.P(),nmiss_mass,weight);
        MMnmiss_IMnpim_dE_woK0_woSid_won->Fill(LVec_pim_n.M(),nmiss_mass,weight);
        MMnmiss_Momnpim_dE_woK0_woSid_won->Fill(LVec_pim_n.P(),nmiss_mass,weight);
        Momnpim_Momnpip_dE_woK0_woSid_won->Fill(LVec_pip_n.P(),LVec_pim_n.P(),weight);
        Momnpim_Mompippim_dE_woK0_woSid_won->Fill(LVec_pip_pim.P(),LVec_pim_n.P(),weight);
        Momnpip_Mompippim_dE_woK0_woSid_won->Fill(LVec_pip_pim.P(),LVec_pip_n.P(),weight);
        MMnmiss_IMpippim_dE_woK0_woSid_won->Fill(LVec_pip_pim.M(),nmiss_mass,weight);
        IMnpim_IMnpip_dE_woK0_woSid_won->Fill(LVec_pip_n.M(),LVec_pim_n.M(),weight);
        IMpippim_IMnpip_woK0_woSid_won->Fill(LVec_pip_n.M(),LVec_pip_pim.M(),weight);
        IMpippim_IMnpim_woK0_woSid_won->Fill(LVec_pim_n.M(),LVec_pip_pim.M(),weight);
        MMnmiss_IMnpipi_woK0_woSid_won->Fill(LVec_pip_pim_n.M(),nmiss_mass,weight);
        MMnmiss_Momnpipi_woK0_woSid_won->Fill(LVec_pip_pim_n.P(),nmiss_mass,weight);
        q_IMnpipi_woK0_woSid_won->Fill(LVec_pip_pim_n.M(),qkn.P(),weight);
        q_IMpippim_woK0_woSid_won->Fill(LVec_pip_pim.M(),qkn.P(),weight);
        q_nmom_woK0_woSid_won->Fill((*LVec_n).P(),qkn.P(),weight);
        q_IMnpip_woK0_woSid_won->Fill(LVec_pip_n.M(),qkn.P(),weight);
        q_IMnpim_woK0_woSid_won->Fill(LVec_pim_n.M(),qkn.P(),weight);
        //q_MMnmiss_woK0_woSid_won->Fill(nmiss_mass,qkn.P(),weight);
        nmom_MMnmiss_woK0_woSid_won->Fill(nmiss_mass,(*LVec_n).P(),weight);
        MMnmiss_Mompippim_dE_woK0_woSid_won->Fill(LVec_pip_pim.P(),nmiss_mass,weight);
        MMom_MMass_woK0_woSid_won->Fill(LVec_npipimiss.M(),LVec_npipimiss.P(),weight);
        //pipmom_MMnmiss_dE_woK0_woSid_won->Fill(nmiss_mass,(*LVec_pip).P(),weight);
        //pipmom_pimmom_dE_woK0_woSid_won->Fill((*LVec_pim).P(), (*LVec_pip).P(),weight);
        //pimmom_MMnmiss_dE_woK0_woSid_won->Fill(nmiss_mass,(*LVec_pim).P(),weight);
        nmom_cosn_woK0_woSid_won->Fill(cos_ncdslab,(*LVec_n).P(),weight);
        nmom_cospip_woK0_woSid_won->Fill(cos_pip,(*LVec_n).P(),weight);
        nmom_cospim_woK0_woSid_won->Fill(cos_pim,(*LVec_n).P(),weight);
        nmom_cospippim_woK0_woSid_won->Fill(cos_pippim,(*LVec_n).P(),weight);
        nmom_phinpip_woK0_woSid_won->Fill(phi_npip,(*LVec_n).P(),weight);
        nmom_phinpim_woK0_woSid_won->Fill(phi_npim,(*LVec_n).P(),weight);
        nmom_phipip_woK0_woSid_won->Fill(phi_pip,(*LVec_n).P(),weight);
        nmom_phipim_woK0_woSid_won->Fill(phi_pim,(*LVec_n).P(),weight);
        nmom_phin_woK0_woSid_won->Fill(phi_n,(*LVec_n).P(),weight);
        nmom_pipmom_woK0_woSid_won->Fill((*LVec_pip).P(),(*LVec_n).P(),weight);
        nmom_pimmom_woK0_woSid_won->Fill((*LVec_pim).P(),(*LVec_n).P(),weight);
        nmom_IMnpip_dE_woK0_woSid_won->Fill(LVec_pip_n.M(),(*LVec_n).P(),weight);
        nmom_IMnpim_dE_woK0_woSid_won->Fill(LVec_pim_n.M(),(*LVec_n).P(),weight);
        nmom_IMpippim_woK0_woSid_won->Fill(LVec_pip_pim.M(),(*LVec_n).P(),weight);
      }
      if(SigmaPcutFlag[sigmacuttype]) {
        MMnmiss_IMnpipi_woK0_wSid_Sp->Fill(LVec_pip_pim_n.M(), nmiss_mass,weight);
      }
      if(SigmaMcutFlag[sigmacuttype]) {
        MMnmiss_IMnpipi_woK0_wSid_Sm->Fill(LVec_pip_pim_n.M(), nmiss_mass,weight);
      }
    }//K0reject, NbetaOK ,NdEOk

    //std::cout << __LINE__ << std::endl;
    if(K0rejectFlag && NBetaOK && NdEOK && MissNFlag) {
      IMnpim_IMnpip_dE_woK0_n->Fill(LVec_pip_n.M(),LVec_pim_n.M(),weight);
      //IMnpip_CDHphi_dE_woK0_n->Fill((*CDH_Pos).Phi(),LVec_pip_n.M(),weight);
      //IMnpip_CDHz_dE_woK0_n->Fill((*CDH_Pos).z(),LVec_pip_n.M(),weight);
      IMnpip_DCApip_dE_woK0_n->Fill(dca_pip_beam,LVec_pip_n.M(),weight);
      IMnpip_DCApim_dE_woK0_n->Fill(dca_pim_beam,LVec_pip_n.M(),weight);
      IMnpim_DCApip_dE_woK0_n->Fill(dca_pip_beam,LVec_pim_n.M(),weight);
      IMnpim_DCApim_dE_woK0_n->Fill(dca_pim_beam,LVec_pim_n.M(),weight);
      //IMnpim_CDHphi_dE_woK0_n->Fill((*CDH_Pos).Phi(),LVec_pim_n.M(),weight);
      //IMnpim_CDHz_dE_woK0_n->Fill((*CDH_Pos).z(),LVec_pim_n.M(),weight);
      IMnpip_IMnpipi_woK0_n->Fill(LVec_pip_pim_n.M(),LVec_pip_n.M(),weight);
      IMnpim_IMnpipi_woK0_n->Fill(LVec_pip_pim_n.M(),LVec_pim_n.M(),weight);
      IMnpip_DCApipibeam_woK0_n->Fill(dca_pipibeam,LVec_pip_n.M(),weight);
      IMnpim_DCApipibeam_woK0_n->Fill(dca_pipibeam,LVec_pim_n.M(),weight);
      for(int igap=0; igap<ngap; igap++) {
        if(SigmaPsideFlag[sidebandtype][igap] || SigmaMsideFlag[sidebandtype][igap]) {
          //IMnpim_IMnpip_dE_woK0_n_side[igap]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
        }
        if(SigmaPsideLowFlag[igap] && SigmaPsideFlag[sidebandtype][igap]) {
          //IMnpim_IMnpip_dE_woK0_n_Sp_side[LOWside][igap]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
        }
        if(SigmaPsideHighFlag[igap] && SigmaPsideFlag[sidebandtype][igap]) {
          //IMnpim_IMnpip_dE_woK0_n_Sp_side[HIGHside][igap]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
        }
        if(SigmaMsideLowFlag[igap] && SigmaMsideFlag[sidebandtype][igap]) {
          //IMnpim_IMnpip_dE_woK0_n_Sm_side[LOWside][igap]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
        }
        if(SigmaMsideHighFlag[igap] && SigmaMsideFlag[sidebandtype][igap]) {
          //IMnpim_IMnpip_dE_woK0_n_Sm_side[HIGHside][igap]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
        }
      }
      
      /*
      for(int izone=0; izone<nzone; izone++) {
        if(sidebandtype==0) {
          if( !SigmaCrossPsideWideFlagLeft[izone] && !SigmaCrossPsideWideFlagRight[izone]  && SigmaPsideWideFlag[izone]) {
            //IMnpim_IMnpip_dE_woK0_n_Sp_sidewide[izone]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
            //q_IMnpipi_woK0_wSid_n_Sp_sidewide[izone]->Fill(LVec_pip_pim_n.M(),qkn.P());
          }
          if( !SigmaCrossMsideWideFlagTop[izone] && !SigmaCrossMsideWideFlagBottom[izone]  && SigmaMsideWideFlag[izone]) {
            IMnpim_IMnpip_dE_woK0_n_Sm_sidewide[izone]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
            q_IMnpipi_woK0_wSid_n_Sm_sidewide[izone]->Fill(LVec_pip_pim_n.M(),qkn.P());
          }
        } else if(sidebandtype==1) {
          if( !SigmaMFlag  && SigmaPsideWideFlag[izone]) {
            IMnpim_IMnpip_dE_woK0_n_Sp_sidewide[izone]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
            q_IMnpipi_woK0_wSid_n_Sp_sidewide[izone]->Fill(LVec_pip_pim_n.M(),qkn.P());
          }
          if( !SigmaPFlag  && SigmaMsideWideFlag[izone]) {
            IMnpim_IMnpip_dE_woK0_n_Sm_sidewide[izone]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
            q_IMnpipi_woK0_wSid_n_Sm_sidewide[izone]->Fill(LVec_pip_pim_n.M(),qkn.P());
          }
        } else if(sidebandtype==2) {
          if( !SigmawideMFlag  && SigmaPsideWideFlag[izone]) {
            IMnpim_IMnpip_dE_woK0_n_Sp_sidewide[izone]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
            q_IMnpipi_woK0_wSid_n_Sp_sidewide[izone]->Fill(LVec_pip_pim_n.M(),qkn.P());
          }
          if( !SigmawidePFlag  && SigmaMsideWideFlag[izone]) {
            IMnpim_IMnpip_dE_woK0_n_Sm_sidewide[izone]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
            q_IMnpipi_woK0_wSid_n_Sm_sidewide[izone]->Fill(LVec_pip_pim_n.M(),qkn.P());
          }
        }
      }
      */

      //std::cout << __LINE__ << std::endl;
      //MMnpim_MMnpip_woK0_n->Fill(LVec_pip_nmiss.M(),LVec_pim_nmiss.M(),weight);
      MMnpim_MMnpip_woK0_n->Fill(LVec_pimmiss_nmiss.M(),LVec_pipmiss_nmiss.M());
      nmom_IMnpim_dE_woK0_n->Fill(LVec_pim_n.M(),(*LVec_n).P(),weight);
      nmom_IMnpip_dE_woK0_n->Fill(LVec_pip_n.M(),(*LVec_n).P(),weight);
      if(!SigmawideMFlag) {
        MMnmiss_IMnpip_dE_woK0_woSm_n->Fill(LVec_pip_n.M(),nmiss_mass,weight);
      }
      if(!SigmawidePFlag) {
        MMnmiss_IMnpim_dE_woK0_woSp_n->Fill(LVec_pim_n.M(),nmiss_mass,weight);
      }
      if(SigmaPFlag || SigmaMFlag) {
        IMnpim_IMnpip_dE_woK0_wSid_n->Fill(LVec_pip_n.M(),LVec_pim_n.M(),weight);
        //IMnpim_IMnpip_dE_woK0_wSid_n_bin[binnum]->Fill(LVec_pip_n.M(),LVec_pim_n.M(),weight);
        //MMnpim_MMnpip_woK0_wSid_n->Fill(LVec_pip_nmiss.M(),LVec_pim_nmiss.M(),weight);
        MMnpim_MMnpip_woK0_wSid_n->Fill(LVec_pimmiss_nmiss.M(),LVec_pipmiss_nmiss.M());
        MMnmiss_IMpippim_dE_woK0_wSid_n->Fill(LVec_pip_pim.M(),nmiss_mass,weight);
        Momnpim_Momnpip_dE_woK0_wSid_n->Fill(LVec_pip_n.P(),LVec_pim_n.P(),weight);
        Momnpim_Mompippim_dE_woK0_wSid_n->Fill(LVec_pip_pim.P(),LVec_pim_n.P(),weight);
        Momnpip_Mompippim_dE_woK0_wSid_n->Fill(LVec_pip_pim.P(),LVec_pip_n.P(),weight);
        //dE_IMnpipi_woK0_wSid_n->Fill(LVec_pip_pim_n.M(),dE,weight);
        Cosn_IMnpipi_woK0_wSid_n->Fill(LVec_pip_pim_n.M(),cos_nmissCM,weight);
        q_IMpiSigma_woK0_wSid_n_genacc->Fill(LVec_piSigma_react.M()/1000.,qkn_react.P()/1000.);
        q_IMnpipi_woK0_wSid_n_acc->Fill(LVec_pip_pim_n_mc.M(),qkn_mc.P(),weight);
        q_IMnpipi_woK0_wSid_n_acc_reco->Fill(LVec_pip_pim_n.M(),qkn.P(),weight);
        q_IMnpipi_woK0_wSid_n->Fill(LVec_pip_pim_n.M(),qkn.P(),weight);
        q_nmom_woK0_wSid_n->Fill((*LVec_n).P(), qkn.P(),weight);
        //nmom_IMnpipi_woK0_wSid_n->Fill(LVec_pip_pim_n.M(),(*LVec_n).P(),weight);
        IMpippim_IMnpip_woK0_wSid_n->Fill(LVec_pip_n.M(),LVec_pip_pim.M(),weight);
        if(!SigmawidePFlag){
          MMnmiss_IMnpim_dE_woK0_wSid_n_woSp_wbin[wbinnum]->Fill(LVec_pim_n.M(),nmiss_mass,weight);
          IMpippim_IMnpip_woK0_wSid_n_woSp_wbin[wbinnum]->Fill(LVec_pip_n.M(),LVec_pip_pim.M(),weight);
          q_nmom_woK0_wSid_n_woSp_wbin[wbinnum]->Fill((*LVec_n).P(), qkn.P(),weight);
          IMpippim_IMnpim_woK0_wSid_n_woSp_wbin[wbinnum]->Fill(LVec_pim_n.M(),LVec_pip_pim.M(),weight);
          q_IMnpipi_woK0_wSid_n_woSp->Fill(LVec_pip_pim_n.M(),qkn.P(),weight);
          IMnpim_IMnpip_dE_woK0_wSid_n_woSp->Fill(LVec_pip_n.M(),LVec_pim_n.M(),weight);
          IMnpim_IMnpip_dE_woK0_wSid_n_woSp_wbin[wbinnum]->Fill(LVec_pip_n.M(),LVec_pim_n.M(),weight);
        }
        if(!SigmawideMFlag){
          MMnmiss_IMnpip_dE_woK0_wSid_n_woSm_wbin[wbinnum]->Fill(LVec_pip_n.M(),nmiss_mass,weight);
          IMpippim_IMnpip_woK0_wSid_n_woSm_wbin[wbinnum]->Fill(LVec_pip_n.M(),LVec_pip_pim.M(),weight);
          q_nmom_woK0_wSid_n_woSm_wbin[wbinnum]->Fill((*LVec_n).P(), qkn.P(),weight);
          IMpippim_IMnpim_woK0_wSid_n_woSm_wbin[wbinnum]->Fill(LVec_pim_n.M(),LVec_pip_pim.M(),weight);
          q_IMnpipi_woK0_wSid_n_woSm->Fill(LVec_pip_pim_n.M(),qkn.P(),weight);
          IMnpim_IMnpip_dE_woK0_wSid_n_woSm->Fill(LVec_pip_n.M(),LVec_pim_n.M(),weight);
          IMnpim_IMnpip_dE_woK0_wSid_n_woSm_wbin[wbinnum]->Fill(LVec_pip_n.M(),LVec_pim_n.M(),weight);
        }

        IMpippim_IMnpim_woK0_wSid_n->Fill(LVec_pim_n.M(),LVec_pip_pim.M(),weight);
        nmom_cosn_woK0_wSid_n->Fill(cos_ncdslab,(*LVec_n).P(),weight);
        nmom_cosnmiss_woK0_wSid_n->Fill(cos_nmissCM,(*LVec_n).P(),weight);
        DCA_pip_beam->Fill( dca_pip_beam,weight);
        DCA_pim_beam->Fill( dca_pim_beam,weight );
        DCA_pip_pim->Fill(dca_pip_pim,weight);
        pipmom_IMnpipi_woK0_wSid_n->Fill(LVec_pip_pim_n.M(),(*LVec_pip).P(),weight);
        pimmom_IMnpipi_woK0_wSid_n->Fill(LVec_pip_pim_n.M(),(*LVec_pim).P(),weight);
        nmom_nmissmom_woK0_wSid_n->Fill(LVec_npipimiss.P(),(*LVec_n).P(),weight);
        diff2d_CDC_CDH_pim_woK0_wSid_n->Fill(diffPhinpim,diffpim.z(),weight);
        diff2d_CDC_CDH_pip_woK0_wSid_n->Fill(diffPhinpip,diffpip.z(),weight);
        
        //std::cout << __LINE__ << std::endl;
        if(SimSpmode || SimSmmode || SimK0nnmode){
          //std::cout << __LINE__ << std::endl;
          double diffIMnpim_recomc = LVec_pim_n.M()- LVec_pim_n_mc.M();
          double diffIMnpip_recomc = LVec_pip_n.M()- LVec_pip_n_mc.M();
          double diffMMnmiss_recomc = nmiss_mass - (*mcmom_nmiss).M();
          double diffnmom_recomc = (*LVec_n).P() - (*mcmom_ncds).P();
          //diff2D_MMnmiss_IMnpim_recomc_woK0_wSid_n->Fill(diffIMnpim_recomc,diffMMnmiss_recomc);
          //diff2D_MMnmiss_IMnpip_recomc_woK0_wSid_n->Fill(diffIMnpip_recomc,diffMMnmiss_recomc);
          //diff2D_nmom_IMnpim_recomc_woK0_wSid_n->Fill(diffIMnpim_recomc,diffnmom_recomc);
          //diff2D_nmom_IMnpip_recomc_woK0_wSid_n->Fill(diffIMnpip_recomc,diffnmom_recomc);
          if(!IsFakeN1){
            diff2D_MMnmiss_IMnpim_recomc_woK0_wSid_n_fake1->Fill(diffIMnpim_recomc,diffMMnmiss_recomc);
            diff2D_MMnmiss_IMnpip_recomc_woK0_wSid_n_fake1->Fill(diffIMnpip_recomc,diffMMnmiss_recomc);
            diff2D_nmom_IMnpim_recomc_woK0_wSid_n_fake1->Fill(diffIMnpim_recomc,diffnmom_recomc);
            diff2D_nmom_IMnpip_recomc_woK0_wSid_n_fake1->Fill(diffIMnpip_recomc,diffnmom_recomc);
          }
        }
      }


      //0: diagonal cut
      //1: 3 sigma cut
      //2: 5 simga cut
      if(SigmaPcutFlag[sigmacuttype]) {
        IMnpim_IMnpip_dE_woK0_n_Sp->Fill(LVec_pip_n.M(),LVec_pim_n.M(),weight);
        q_IMnpipi_woK0_wSid_n_Sp->Fill(LVec_pip_pim_n.M(),qkn.P(),weight);
        nmom_IMnpipi_woK0_wSid_n_Sp->Fill(LVec_pip_pim_n.M(),(*LVec_n).P(),weight);
        q_IMpiSigma_woK0_wSid_n_Sp_genacc->Fill(LVec_piSigma_react.M()/1000.,qkn_react.P()/1000.);
        q_IMnpipi_woK0_wSid_n_Sp_acc->Fill(LVec_pip_pim_n_mc.M(),qkn_mc.P());
        q_IMnpipi_woK0_wSid_n_Sp_acc_reco->Fill(LVec_pip_pim_n.M(),qkn.P());
        IMnpip_DCApipibeam_woK0_n_Sp->Fill(dca_pipibeam,LVec_pip_n.M(),weight);
        nmom_Momnpip_woK0_n_Sp->Fill(LVec_pip_n.P(),(*LVec_n).P(),weight);
        if(SimSpmode) {
          q_IMnpipi_woK0_wSid_n_Sp_mc->Fill(LVec_pip_pim_n_mc.M(),qkn_mc.P());
          diff_IMnpipi_woK0_wSid_n_Sp->Fill(LVec_pip_pim_n.M(),LVec_pip_pim_n.M()-LVec_pip_pim_n_mc.M(),weight);
          diff_q_woK0_wSid_n_Sp->Fill(qkn.P(),qkn.P()-qkn_mc.P(),weight);
        }
      }
      //if(SigmaPcutFlag[sigmacuttype] || (!SigmaPFlag)) {
      //  IMnpim_IMnpip_dE_woK0_n_Sp_bg->Fill(LVec_pip_n.M(),LVec_pim_n.M(),weight);
      //}

      /*
      for(int itype=0; itype<3; itype++) {
        for(int igap=0; igap<ngap; igap++) {
          if(SigmaPsideFlag[itype][igap] && SigmaPsideLowFlag[igap] ) {
            q_IMnpipi_woK0_wSid_n_Sp_side[itype][LOWside][igap]->Fill(LVec_pip_pim_n.M(),qkn.P());
          }
          if(SigmaPsideFlag[itype][igap] && SigmaPsideHighFlag[igap] ) {
            q_IMnpipi_woK0_wSid_n_Sp_side[itype][HIGHside][igap]->Fill(LVec_pip_pim_n.M(),qkn.P());
          }
        }
      }*/

      if(SigmaMcutFlag[sigmacuttype]) {
        IMnpim_IMnpip_dE_woK0_n_Sm->Fill(LVec_pip_n.M(),LVec_pim_n.M());
        q_IMnpipi_woK0_wSid_n_Sm->Fill(LVec_pip_pim_n.M(),qkn.P());
        nmom_IMnpipi_woK0_wSid_n_Sm->Fill(LVec_pip_pim_n.M(),(*LVec_n).P());
        q_IMpiSigma_woK0_wSid_n_Sm_genacc->Fill(LVec_piSigma_react.M()/1000.,qkn_react.P()/1000.);
        q_IMnpipi_woK0_wSid_n_Sm_acc->Fill(LVec_pip_pim_n_mc.M(),qkn_mc.P());
        q_IMnpipi_woK0_wSid_n_Sm_acc_reco->Fill(LVec_pip_pim_n.M(),qkn.P());
        IMnpim_DCApipibeam_woK0_n_Sm->Fill(dca_pipibeam,LVec_pim_n.M());
        nmom_Momnpim_woK0_n_Sm->Fill(LVec_pim_n.P(),(*LVec_n).P());
        if(SimSmmode) {
          q_IMnpipi_woK0_wSid_n_Sm_mc->Fill(LVec_pip_pim_n_mc.M(),qkn_mc.P());
          diff_IMnpipi_woK0_wSid_n_Sm->Fill(LVec_pip_pim_n.M(),LVec_pip_pim_n.M()-LVec_pip_pim_n_mc.M());
          diff_q_woK0_wSid_n_Sm->Fill(qkn.P(),qkn.P()-qkn_mc.P());
        }
      }
      //if(SigmaMcutFlag[sigmacuttype] || (!SigmaMFlag)) {
      //  IMnpim_IMnpip_dE_woK0_n_Sm_bg->Fill(LVec_pip_n.M(),LVec_pim_n.M(),weight);
      //}

      /*
      for(int igap=0; igap<ngap; igap++) {
        for(int itype=0; itype<3; itype++) {
          if(SigmaMsideFlag[itype][igap] && SigmaMsideLowFlag[igap]) {
            q_IMnpipi_woK0_wSid_n_Sm_side[itype][LOWside][igap]->Fill(LVec_pip_pim_n.M(),qkn.P());
          }
          if(SigmaMsideFlag[itype][igap] && SigmaMsideHighFlag[igap]) {
            q_IMnpipi_woK0_wSid_n_Sm_side[itype][HIGHside][igap]->Fill(LVec_pip_pim_n.M(),qkn.P());
          }
        }
        if(SigmaPsideFlag[sidebandtype][igap] || SigmaMsideFlag[sidebandtype][igap]) {
          q_IMnpipi_woK0_wSid_n_side[igap]->Fill(LVec_pip_pim_n.M(),qkn.P());
        }
      }//for igap
      */
    }//if K0rejectFlag && NBetaOK && NdEOK && MissNFlag0
    //---removing K0 END----------------------------------------------
    //if(i> 1.00e+06) break;
  }//for ievt
  //--- Filling Histogram END --------------------------------------------------

  //std::cout << __LINE__ << std::endl;

  //----------------------------------------------------------------------------
  //---Drawing Part
  //----------------------------------------------------------------------------
  TCanvas *cMMnmiss_IMnpip_dE = new TCanvas("cMMnmiss_IMnpip_dE","MMnmiss_IMnpip_dE");
  cMMnmiss_IMnpip_dE->cd();
  //MMnmiss_IMnpip_dE->GetXaxis()->SetRangeUser(1.05,1.5);
  //MMnmiss_IMnpip_dE->GetYaxis()->SetRangeUser(0.6,1.5);
  MMnmiss_IMnpip_dE->Draw("colz");

  TCanvas *cMMnmiss_IMnpip_dE_woK0 = new TCanvas("cMMnmiss_IMnpip_dE_woK0","MMnmiss_IMnpip_dE_woK0");
  cMMnmiss_IMnpip_dE_woK0->cd();
  //MMnmiss_IMnpip_dE_woK0->GetXaxis()->SetRangeUser(1.05,1.5);
  //MMnmiss_IMnpip_dE_woK0->GetYaxis()->SetRangeUser(0.6,1.5);
  MMnmiss_IMnpip_dE_woK0->Draw("colz");

  TCanvas *cMMnmiss_IMnpim_dE = new TCanvas("cMMnmiss_IMnpim_dE","MMnmiss_IMnpim_dE");
  cMMnmiss_IMnpim_dE->cd();
  //MMnmiss_IMnpim_dE->GetXaxis()->SetRangeUser(1.05,1.5);
  //MMnmiss_IMnpim_dE->GetYaxis()->SetRangeUser(0.6,1.5);
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
  //MMnmiss_IMnpim_dE_woK0->GetXaxis()->SetRangeUser(1.05,1.5);
  //MMnmiss_IMnpim_dE_woK0->GetYaxis()->SetRangeUser(0.6,1.5);
  MMnmiss_IMnpim_dE_woK0->Draw("colz");

  //Sigma reconstruction
  TCanvas *cIMnpim_IMnpip_dE_n = new TCanvas("cIMnpim_IMnpip_dE_n","IMnpim_IMnpip_dE_n");
  cIMnpim_IMnpip_dE_n->cd();
  //IMnpim_IMnpip_dE_n->GetXaxis()->SetRangeUser(0,1.7);
  //IMnpim_IMnpip_dE_n->GetYaxis()->SetRangeUser(0,1.7);
  IMnpim_IMnpip_dE_n->Draw("colz");

  TCanvas *cIMnpim_IMnpip_dE_woK0_n = new TCanvas("cIMnpim_IMnpip_dE_woK0_n","IMnpim_IMnpip_dE_woK0_n");
  cIMnpim_IMnpip_dE_woK0_n->cd();
  //IMnpim_IMnpip_dE_woK0_n->GetXaxis()->SetRangeUser(0,1.7);
  //IMnpim_IMnpip_dE_woK0_n->GetYaxis()->SetRangeUser(0,1.7);
  IMnpim_IMnpip_dE_woK0_n->Draw("colz");

  TCanvas *cIMnpim_IMnpip_dE_woK0_n_px = new TCanvas("cIMnpim_IMnpip_dE_woK0_n_px","IMnpim_IMnpip_dE_woK0_n_px");
  cIMnpim_IMnpip_dE_woK0_n_px->cd();
  TH1D* IMnpim_IMnpip_dE_woK0_n_px = (TH1D*)IMnpim_IMnpip_dE_woK0_n->ProjectionX("IMnpim_IMnpip_dE_woK0_n_px");
  IMnpim_IMnpip_dE_woK0_n_px->Draw("HE");
  TH1D* IMnpim_IMnpip_dE_woK0_n_px_cut = (TH1D*)IMnpim_IMnpip_dE_woK0_n_px->Clone("IMnpim_IMnpip_dE_woK0_n_px_cut");
  IMnpim_IMnpip_dE_woK0_n_px_cut->SetFillColor(2);
  IMnpim_IMnpip_dE_woK0_n_px_cut->GetXaxis()->SetRangeUser(anacuts::Sigmap_MIN,anacuts::Sigmap_MAX);
  IMnpim_IMnpip_dE_woK0_n_px_cut->Draw("HEsame");
// TH1D* IMnpim_IMnpip_dE_woK0_n_px_Sm = (TH1D*)IMnpim_IMnpip_dE_woK0_n->ProjectionX("IMnpim_IMnpip_dE_woK0_n_px_Sm",
//     IMnpim_IMnpip_dE_woK0_n->GetYaxis()->FindBin(anacuts::Sigmam_MIN),
//     IMnpim_IMnpip_dE_woK0_n->GetYaxis()->FindBin(anacuts::Sigmam_MAX));
// IMnpim_IMnpip_dE_woK0_n_px_Sm->SetFillColor(3);
// IMnpim_IMnpip_dE_woK0_n_px_Sm->Draw("HEsame");


  TCanvas *cIMnpim_IMnpip_dE_woK0_n_py = new TCanvas("cIMnpim_IMnpip_dE_woK0_n_py","IMnpim_IMnpip_dE_woK0_n_py");
  cIMnpim_IMnpip_dE_woK0_n_py->cd();
  TH1D* IMnpim_IMnpip_dE_woK0_n_py = (TH1D*)IMnpim_IMnpip_dE_woK0_n->ProjectionY("IMnpim_IMnpip_dE_woK0_n_py");
  IMnpim_IMnpip_dE_woK0_n_py->Draw("HE");
  TH1D* IMnpim_IMnpip_dE_woK0_n_py_cut = (TH1D*)IMnpim_IMnpip_dE_woK0_n_py->Clone("IMnpim_IMnpip_dE_woK0_n_py_cut");
  IMnpim_IMnpip_dE_woK0_n_py_cut->SetFillColor(3);
  IMnpim_IMnpip_dE_woK0_n_py_cut->GetXaxis()->SetRangeUser(anacuts::Sigmam_MIN,anacuts::Sigmam_MAX);
  IMnpim_IMnpip_dE_woK0_n_py_cut->Draw("HEsame");
//  TH1D* IMnpim_IMnpip_dE_woK0_n_py_Sp = (TH1D*)IMnpim_IMnpip_dE_woK0_n->ProjectionY("IMnpim_IMnpip_dE_woK0_n_py_Sp",
//      IMnpim_IMnpip_dE_woK0_n->GetXaxis()->FindBin(anacuts::Sigmap_MIN),
//      IMnpim_IMnpip_dE_woK0_n->GetXaxis()->FindBin(anacuts::Sigmap_MAX));
//  IMnpim_IMnpip_dE_woK0_n_py_Sp->SetFillColor(2);
//  IMnpim_IMnpip_dE_woK0_n_py_Sp->Draw("HEsame");

  //TCanvas *cCDHphi_betainv_fid = new TCanvas("cCDHphi_betainv_fid","CDHphi_betainv_fid");
  //cCDHphi_betainv_fid->cd();
  //CDHphi_betainv_fid->Draw("colz");

  //TCanvas *cCDHz_betainv_fid = new TCanvas("cCDHz_betainv_fid","CDHz_betainv_fid");
  //cCDHz_betainv_fid->cd();
  //CDHz_betainv_fid->Draw("colz");

  //TCanvas *cnmom_CDHphi = new TCanvas("cnmom_CDHphi","nmom_CDHphi");
  //cnmom_CDHphi->cd();
  //nmom_CDHphi->Draw("colz");

  //TCanvas *cIMnpip_CDHphi_dE_woK0_n = new TCanvas("cIMnpip_CDHphi_dE_woK0_n","IMnpip_CDHphi_dE_woK0_n");
  //cIMnpip_CDHphi_dE_woK0_n->cd();
  //IMnpip_CDHphi_dE_woK0_n->Draw("colz");

  //TCanvas *cIMnpip_CDHz_dE_woK0_n = new TCanvas("cIMnpip_CDHz_dE_woK0_n","IMnpip_CDHz_dE_woK0_n");
  //cIMnpip_CDHz_dE_woK0_n->cd();
  //IMnpip_CDHz_dE_woK0_n->RebinX(5);
  //IMnpip_CDHz_dE_woK0_n->Draw("colz");

  //TCanvas *cpipmom_IMpippim_dE_n = new TCanvas("cpipmom_IMpippim_dE_n","pipmom_IMpippim_dE_n");
  //cpipmom_IMpippim_dE_n->cd();
  //pipmom_IMpippim_dE_n->Draw("colz");


  //TCanvas *cpipmom_IMnpip_dE_n = new TCanvas("cpipmom_IMnpip_dE_n","pipmom_IMnpip_dE_n");
  //cpipmom_IMnpip_dE_n->cd();
  //pipmom_IMnpip_dE_n->Draw("colz");
  
  /*
  TCanvas *cpipmom_IMnpip_dE_n_py = new TCanvas("cpipmom_IMnpip_dE_n_py","pipmom_IMnpip_dE_n_py");
  cpipmom_IMnpip_dE_n_py->cd();
  TH1D *pipmom_IMnpip_py = (TH1D*)pipmom_IMnpip->ProjectionY();
  TH1D *pipmom_IMnpip_dE_py = (TH1D*)pipmom_IMnpip_dE->ProjectionY();
  TH1D *pipmom_IMnpip_dE_n_py = (TH1D*)pipmom_IMnpip_dE_n->ProjectionY();
  pipmom_IMnpip_py->Draw("HE");
  pipmom_IMnpip_dE_py->SetLineColor(2);
  pipmom_IMnpip_dE_py->Draw("HEsame");
  pipmom_IMnpip_dE_n_py->SetLineColor(3);
  pipmom_IMnpip_dE_n_py->Draw("HEsame");
  */

  //TCanvas *cpipmom_MMnmiss_dE = new TCanvas("cpipmom_MMnmiss_dE","pipmom_MMnmiss_dE");
  //cpipmom_MMnmiss_dE->cd();
  //pipmom_MMnmiss_dE->Draw("colz");

  //TCanvas *cpimmom_IMpippim_dE_n = new TCanvas("cpimmom_IMpippim_dE_n","pimmom_IMpippim_dE_n");
  //cpimmom_IMpippim_dE_n->cd();
  //pimmom_IMpippim_dE_n->Draw("colz");

  //TCanvas *cpimmom_IMnpim_dE_n = new TCanvas("cpimmom_IMnpim_dE_n","pimmom_IMnpim_dE_n");
  //cpimmom_IMnpim_dE_n->cd();
  //pimmom_IMnpim_dE_n->Draw("colz");
  
  /*
  TCanvas *cpimmom_IMnpim_dE_n_py = new TCanvas("cpimmom_IMnpim_dE_n_py","pimmom_IMnpim_dE_n_py");
  cpimmom_IMnpim_dE_n_py->cd();
  TH1D *pimmom_IMnpim_py = (TH1D*)pimmom_IMnpim->ProjectionY();
  TH1D *pimmom_IMnpim_dE_py = (TH1D*)pimmom_IMnpim_dE->ProjectionY();
  TH1D *pimmom_IMnpim_dE_n_py = (TH1D*)pimmom_IMnpim_dE_n->ProjectionY();
  pimmom_IMnpim_py->Draw("HE");
  pimmom_IMnpim_dE_py->SetLineColor(2);
  pimmom_IMnpim_dE_py->Draw("HEsame");
  pimmom_IMnpim_dE_n_py->SetLineColor(3);
  pimmom_IMnpim_dE_n_py->Draw("HEsame");
  */

  //TCanvas *cIMnpim_CDHphi_dE_woK0_n = new TCanvas("cIMnpim_CDHphi_dE_woK0_n","IMnpim_CDHphi_dE_woK0_n");
  //cIMnpim_CDHphi_dE_woK0_n->cd();
  //IMnpim_CDHphi_dE_woK0_n->Draw("colz");

  //TCanvas *cIMnpim_CDHphi_dE_woK0_n_px = new TCanvas("cIMnpim_CDHphi_dE_woK0_n_px","IMnpim_CDHphi_dE_woK0_n_px");
  //cIMnpim_CDHphi_dE_woK0_n_px->cd();
  //int bin105 = IMnpim_CDHphi_dE_woK0_n->GetYaxis()->FindBin(1.05);
  //int bin115 = IMnpim_CDHphi_dE_woK0_n->GetYaxis()->FindBin(1.15);
  //TH1D *IMnpim_CDHphi_dE_woK0_n_px = (TH1D*) IMnpim_CDHphi_dE_woK0_n->ProjectionX("IMnpim_CDHphi_dE_woK0_n_px",bin105,bin115);
  //IMnpim_CDHphi_dE_woK0_n_px->Draw("HE");

  //TCanvas *cIMnpim_CDHz_dE_woK0_n = new TCanvas("cIMnpim_CDHz_dE_woK0_n","IMnpim_CDHz_dE_woK0_n");
  //cIMnpim_CDHz_dE_woK0_n->cd();
  //IMnpim_CDHz_dE_woK0_n->RebinX(5);
  //IMnpim_CDHz_dE_woK0_n->Draw("colz");

  //TCanvas *cIMnpim_CDHz_dE_woK0_n_px = new TCanvas("cIMnpim_CDHz_dE_woK0_n_px","IMnpim_CDHz_dE_woK0_n_px");
  //cIMnpim_CDHz_dE_woK0_n_px->cd();
  //TH1D *IMnpim_CDHz_dE_woK0_n_px = (TH1D*)IMnpim_CDHz_dE_woK0_n->ProjectionX("IMnpim_CDHz_dE_woK0_n_px",bin105,bin115);
  //IMnpim_CDHz_dE_woK0_n_px->Draw("HE");


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


   /*

  TCanvas *cIMnpim_IMnpip_dE_woK0_n_Sp[ngap];
  for(int igap=0; igap<ngap; igap++) {
    cIMnpim_IMnpip_dE_woK0_n_Sp[igap] = new TCanvas(Form("cIMnpim_IMnpip_dE_woK0_n_Sp_%d",igap),Form("IMnpim_IMnpip_dE_woK0_n_Sp_%d",igap));
    cIMnpim_IMnpip_dE_woK0_n_Sp[igap]->cd();
    IMnpim_IMnpip_dE_woK0_n_Sp->SetMaximum(IMnpim_IMnpip_dE_woK0_n->GetMaximum());
    IMnpim_IMnpip_dE_woK0_n_Sp->Draw("colz");
    //IMnpim_IMnpip_dE_woK0_n_Sp_side[LOWside][igap]->Draw("colsame");
    //IMnpim_IMnpip_dE_woK0_n_Sp_side[HIGHside][igap]->Draw("colsame");
    if(SimSpmode || SimSmmode) break;
  }

  TCanvas *cIMnpim_IMnpip_dE_woK0_n_Sm[ngap];
  for(int igap=0; igap<ngap; igap++) {
    cIMnpim_IMnpip_dE_woK0_n_Sm[igap] = new TCanvas(Form("cIMnpim_IMnpip_dE_woK0_n_Sm_%d",igap),Form("IMnpim_IMnpip_dE_woK0_n_Sm_%d",igap));
    cIMnpim_IMnpip_dE_woK0_n_Sm[igap]->cd();
    IMnpim_IMnpip_dE_woK0_n_Sm->SetMaximum(IMnpim_IMnpip_dE_woK0_n->GetMaximum());
    IMnpim_IMnpip_dE_woK0_n_Sm->GetXaxis()->SetRangeUser(0,1.7);
    IMnpim_IMnpip_dE_woK0_n_Sm->GetYaxis()->SetRangeUser(0,1.7);
    IMnpim_IMnpip_dE_woK0_n_Sm->Draw("colz");
    //IMnpim_IMnpip_dE_woK0_n_Sm_side[LOWside][igap]->Draw("colsame");
    //IMnpim_IMnpip_dE_woK0_n_Sm_side[HIGHside][igap]->Draw("colsame");
    if(SimSpmode || SimSmmode) break;
  }
  */

  /*
  TCanvas *cIMnpim_IMnpip_dE_n_Sp_wide[nzone];
  //TH1D* q_IMnpipi_wSid_n_Sp_sidewide_px[nzone];
  for(int izone=0; izone<nzone; izone++) {
    cIMnpim_IMnpip_dE_n_Sp_wide[izone] = new TCanvas(Form("cIMnpim_IMnpip_dE_n_Sp_wide_%d",izone),Form("IMnpim_IMnpip_dE_n_Sp_wide_%d",izone));
    cIMnpim_IMnpip_dE_n_Sp_wide[izone]->Divide(2,2);
    cIMnpim_IMnpip_dE_n_Sp_wide[izone]->cd(1);
    IMnpim_IMnpip_dE_n_Sp->SetMaximum(IMnpim_IMnpip_dE_woK0_n->GetMaximum());
    IMnpim_IMnpip_dE_n_Sp->GetXaxis()->SetRangeUser(0,1.7);
    IMnpim_IMnpip_dE_n_Sp->GetYaxis()->SetRangeUser(0,1.7);
    IMnpim_IMnpip_dE_n_Sp->Draw("colz");
    //IMnpim_IMnpip_dE_n_Sp_sidewide[izone]->Draw("colsame");
    cIMnpim_IMnpip_dE_n_Sp_wide[izone]->cd(2);
    TH1D* q_IMnpipi_wSid_n_Sp_px = (TH1D*) q_IMnpipi_wSid_n_Sp->ProjectionX();
    q_IMnpipi_wSid_n_Sp_px->Draw("HE");
    //q_IMnpipi_wSid_n_Sp_sidewide_px[izone] = (TH1D*)q_IMnpipi_wSid_n_Sp_sidewide[izone]->ProjectionX(Form("Sp_sidewide_px_%d",izone));
    //q_IMnpipi_wSid_n_Sp_sidewide_px[izone]->SetLineColor(kCyan);
    //q_IMnpipi_wSid_n_Sp_sidewide_px[izone]->SetMarkerColor(kCyan);
    //q_IMnpipi_wSid_n_Sp_sidewide_px[izone]->SetMarkerStyle(22);
    //q_IMnpipi_wSid_n_Sp_sidewide_px[izone]->Draw("HEsame");
    cIMnpim_IMnpip_dE_n_Sp_wide[izone]->cd(3);
    q_IMnpipi_wSid_n_Sp->Draw("colz");
    cIMnpim_IMnpip_dE_n_Sp_wide[izone]->cd(4);
    //q_IMnpipi_wSid_n_Sp_sidewide[izone]->SetMaximum(q_IMnpipi_wSid_n_Sp->GetMaximum());
    //q_IMnpipi_wSid_n_Sp_sidewide[izone]->Draw("colz");

    if(SimSpmode || SimSmmode) break;
  }*/
  
  /*
  TCanvas *cIMnpim_IMnpip_dE_woK0_n_Sp_wide[nzone];
  //TH1D* q_IMnpipi_woK0_wSid_n_Sp_sidewide_px[nzone];
  for(int izone=0; izone<nzone; izone++) {
    cIMnpim_IMnpip_dE_woK0_n_Sp_wide[izone] = new TCanvas(Form("cIMnpim_IMnpip_dE_woK0_n_Sp_wide_%d",izone),Form("IMnpim_IMnpip_dE_woK0_n_Sp_wide_%d",izone));
    cIMnpim_IMnpip_dE_woK0_n_Sp_wide[izone]->Divide(2,2);
    cIMnpim_IMnpip_dE_woK0_n_Sp_wide[izone]->cd(1);
    IMnpim_IMnpip_dE_woK0_n_Sp->SetMaximum(IMnpim_IMnpip_dE_woK0_n->GetMaximum());
    IMnpim_IMnpip_dE_woK0_n_Sp->GetXaxis()->SetRangeUser(0,1.7);
    IMnpim_IMnpip_dE_woK0_n_Sp->GetYaxis()->SetRangeUser(0,1.7);
    IMnpim_IMnpip_dE_woK0_n_Sp->Draw("colz");
    //IMnpim_IMnpip_dE_woK0_n_Sp_sidewide[izone]->Draw("colsame");
    cIMnpim_IMnpip_dE_woK0_n_Sp_wide[izone]->cd(2);
    TH1D* q_IMnpipi_woK0_wSid_n_Sp_px = (TH1D*) q_IMnpipi_woK0_wSid_n_Sp->ProjectionX();
    q_IMnpipi_woK0_wSid_n_Sp_px->Draw("HE");
    //q_IMnpipi_woK0_wSid_n_Sp_sidewide_px[izone] = (TH1D*)q_IMnpipi_woK0_wSid_n_Sp_sidewide[izone]->ProjectionX(Form("Sp_sidewide_woK0_px_%d",izone));
    //q_IMnpipi_woK0_wSid_n_Sp_sidewide_px[izone]->SetLineColor(kCyan);
    //q_IMnpipi_woK0_wSid_n_Sp_sidewide_px[izone]->SetMarkerColor(kCyan);
    //q_IMnpipi_woK0_wSid_n_Sp_sidewide_px[izone]->SetMarkerStyle(22);
    //q_IMnpipi_woK0_wSid_n_Sp_sidewide_px[izone]->Draw("HEsame");
    cIMnpim_IMnpip_dE_woK0_n_Sp_wide[izone]->cd(3);
    q_IMnpipi_woK0_wSid_n_Sp->Draw("colz");
    cIMnpim_IMnpip_dE_woK0_n_Sp_wide[izone]->cd(4);
    //q_IMnpipi_woK0_wSid_n_Sp_sidewide[izone]->SetMaximum(q_IMnpipi_woK0_wSid_n_Sp->GetMaximum());
    //q_IMnpipi_woK0_wSid_n_Sp_sidewide[izone]->Draw("colz");

    if(SimSpmode || SimSmmode) break;
  }

  TCanvas *cIMnpim_IMnpip_dE_n_Sm_wide[nzone];
  //TH1D* q_IMnpipi_wSid_n_Sm_sidewide_px[nzone];
  for(int izone=0; izone<nzone; izone++) {
    cIMnpim_IMnpip_dE_n_Sm_wide[izone] = new TCanvas(Form("cIMnpim_IMnpip_dE_n_Sm_wide_%d",izone),Form("IMnpim_IMnpip_dE_n_Sm_wide_%d",izone));
    cIMnpim_IMnpip_dE_n_Sm_wide[izone]->Divide(2,2);
    cIMnpim_IMnpip_dE_n_Sm_wide[izone]->cd(1);
    IMnpim_IMnpip_dE_n_Sm->SetMaximum(IMnpim_IMnpip_dE_n->GetMaximum());
    IMnpim_IMnpip_dE_n_Sm->GetXaxis()->SetRangeUser(0,1.7);
    IMnpim_IMnpip_dE_n_Sm->GetYaxis()->SetRangeUser(0,1.7);
    IMnpim_IMnpip_dE_n_Sm->Draw("colz");
    //IMnpim_IMnpip_dE_n_Sm_sidewide[izone]->Draw("colsame");
    cIMnpim_IMnpip_dE_n_Sm_wide[izone]->cd(2);
    TH1D* q_IMnpipi_wSid_n_Sm_px = (TH1D*) q_IMnpipi_wSid_n_Sm->ProjectionX();
    q_IMnpipi_wSid_n_Sm_px->Draw("HE");
    //q_IMnpipi_wSid_n_Sm_sidewide_px[izone] = (TH1D*)q_IMnpipi_wSid_n_Sm_sidewide[izone]->ProjectionX(Form("Sm_sidewide_px_%d",izone));
    //q_IMnpipi_wSid_n_Sm_sidewide_px[izone]->SetLineColor(kPink);
    //q_IMnpipi_wSid_n_Sm_sidewide_px[izone]->SetMarkerColor(kPink);
    //q_IMnpipi_wSid_n_Sm_sidewide_px[izone]->SetMarkerStyle(22);
    //q_IMnpipi_wSid_n_Sm_sidewide_px[izone]->Draw("HEsame");
    cIMnpim_IMnpip_dE_n_Sm_wide[izone]->cd(3);
    q_IMnpipi_wSid_n_Sm->Draw("colz");
    cIMnpim_IMnpip_dE_n_Sm_wide[izone]->cd(4);
    //q_IMnpipi_wSid_n_Sm_sidewide[izone]->SetMaximum(q_IMnpipi_woK0_wSid_n_Sm->GetMaximum());
    //q_IMnpipi_wSid_n_Sm_sidewide[izone]->Draw("colz");

    if(SimSpmode || SimSmmode) break;
  }

  TCanvas *cIMnpim_IMnpip_dE_woK0_n_Sm_wide[nzone];
  //TH1D* q_IMnpipi_woK0_wSid_n_Sm_sidewide_px[nzone];
  for(int izone=0; izone<nzone; izone++) {
    cIMnpim_IMnpip_dE_woK0_n_Sm_wide[izone] = new TCanvas(Form("cIMnpim_IMnpip_dE_woK0_n_Sm_wide_%d",izone),Form("IMnpim_IMnpip_dE_woK0_n_Sm_wide_%d",izone));
    cIMnpim_IMnpip_dE_woK0_n_Sm_wide[izone]->Divide(2,2);
    cIMnpim_IMnpip_dE_woK0_n_Sm_wide[izone]->cd(1);
    IMnpim_IMnpip_dE_woK0_n_Sm->SetMaximum(IMnpim_IMnpip_dE_woK0_n->GetMaximum());
    IMnpim_IMnpip_dE_woK0_n_Sm->GetXaxis()->SetRangeUser(0,1.7);
    IMnpim_IMnpip_dE_woK0_n_Sm->GetYaxis()->SetRangeUser(0,1.7);
    IMnpim_IMnpip_dE_woK0_n_Sm->Draw("colz");
    //IMnpim_IMnpip_dE_woK0_n_Sm_sidewide[izone]->Draw("colsame");
    cIMnpim_IMnpip_dE_woK0_n_Sm_wide[izone]->cd(2);
    TH1D* q_IMnpipi_woK0_wSid_n_Sm_px = (TH1D*) q_IMnpipi_woK0_wSid_n_Sm->ProjectionX();
    q_IMnpipi_woK0_wSid_n_Sm_px->Draw("HE");
    //q_IMnpipi_woK0_wSid_n_Sm_sidewide_px[izone] = (TH1D*)q_IMnpipi_woK0_wSid_n_Sm_sidewide[izone]->ProjectionX(Form("Sm_sidewide_woK0_px_%d",izone));
    //q_IMnpipi_woK0_wSid_n_Sm_sidewide_px[izone]->SetLineColor(kPink);
    //q_IMnpipi_woK0_wSid_n_Sm_sidewide_px[izone]->SetMarkerColor(kPink);
    //q_IMnpipi_woK0_wSid_n_Sm_sidewide_px[izone]->SetMarkerStyle(22);
    //q_IMnpipi_woK0_wSid_n_Sm_sidewide_px[izone]->Draw("HEsame");
    cIMnpim_IMnpip_dE_woK0_n_Sm_wide[izone]->cd(3);
    q_IMnpipi_woK0_wSid_n_Sm->Draw("colz");
    cIMnpim_IMnpip_dE_woK0_n_Sm_wide[izone]->cd(4);
    //q_IMnpipi_woK0_wSid_n_Sm_sidewide[izone]->SetMaximum(q_IMnpipi_woK0_wSid_n_Sm->GetMaximum());
    //q_IMnpipi_woK0_wSid_n_Sm_sidewide[izone]->Draw("colz");

    if(SimSpmode || SimSmmode) break;
  }
  */

  /*
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

  //std::cout << "Sigma+ signal region:    " << IMnpim_IMnpip_dE_woK0_n_px_1->Integral() << std::endl;
  //std::cout << "Sigma+ sideband low:     " << IMnpim_IMnpip_dE_woK0_n_px_2->Integral() << std::endl;
  ///std::cout << "Sigma+ sideband low cut: " << IMnpim_IMnpip_dE_woK0_n_Sp_side_0_px->Integral() << std::endl;
  //std::cout << "Sigma+ sideband high:    " << IMnpim_IMnpip_dE_woK0_n_px_3->Integral() << std::endl;
  //std::cout << "Sigma+ sideband high cut:" << IMnpim_IMnpip_dE_woK0_n_Sp_side_1_px->Integral() << std::endl;
  //std::cout << "bg (Integral)           :" << bgsp       << std::endl;
  //std::cout << "bg (Integral)/binw      :" << bgsp/IMnpim_IMnpip_dE_woK0_n_px_1->GetBinWidth(100) << std::endl;
  double trapezoidbgSp = (fbgSp->Eval(anacuts::Sigmap_MIN)+fbgSp->Eval(anacuts::Sigmap_MAX))*
                         (anacuts::Sigmap_MAX-anacuts::Sigmap_MIN)
                         /IMnpim_IMnpip_dE_woK0_n_px_1->GetBinWidth(100)/2.0;
  //std::cout << "bg (trapezoid)          :" << trapezoidbgSp << std::endl;



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

  //std::cout << "Sigma- signal region:    " << IMnpim_IMnpip_dE_woK0_n_py_1->Integral() << std::endl;
  //std::cout << "Sigma- sideband low:     " << IMnpim_IMnpip_dE_woK0_n_py_2->Integral() << std::endl;
  //std::cout << "Sigma- sideband low cut: " << IMnpim_IMnpip_dE_woK0_n_Sm_side_0_py->Integral() << std::endl;
  //std::cout << "Sigma- sideband high:    " << IMnpim_IMnpip_dE_woK0_n_py_3->Integral() << std::endl;
  //std::cout << "Sigma- sideband high cut:" << IMnpim_IMnpip_dE_woK0_n_Sm_side_1_py->Integral() << std::endl;
  //std::cout << "bg (Integral)           :" << bgsm << std::endl;
  //std::cout << "bg (Integral)/binw      :" << bgsm/IMnpim_IMnpip_dE_woK0_n_py_1->GetBinWidth(100)  << std::endl;
  double trapezoidbgSm = (fbgSm->Eval(anacuts::Sigmam_MIN)+fbgSm->Eval(anacuts::Sigmam_MAX))*
                         (anacuts::Sigmam_MAX-anacuts::Sigmam_MIN)
                         /IMnpim_IMnpip_dE_woK0_n_py_1->GetBinWidth(100)/2.0;
  //std::cout << "bg (trapezoid)          :" << trapezoidbgSm    << std::endl;
  */

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

  /*
  TCanvas *cq_IMnpipi_wK0_n = new TCanvas("cq_IMnpipi_wK0_n","q_IMnpipi_wK0_n");
  cq_IMnpipi_wK0_n->cd();
  //q_IMnpipi_wK0_n->GetYaxis()->SetRangeUser(0,1.0);
  q_IMnpipi_wK0_n->Draw("colz");
  fkp->Draw("same");
  */

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

  //TCanvas *cq_IMnpipi_wSid_n_Sp_side_low = new TCanvas("cq_IMnpipi_wSid_n_Sp_side_low","q_IMnpipi_wSid_n_Sp_side_low");
  //cq_IMnpipi_wSid_n_Sp_side_low->cd();
  //q_IMnpipi_wSid_n_Sp_side[sidebandtype][LOWside][0]->Draw("colz");
  //TCanvas *cq_IMnpipi_woK0_wSid_n_Sp_side_low = new TCanvas("cq_IMnpipi_woK0_wSid_n_Sp_side_low","q_IMnpipi_woK0_wSid_n_Sp_side_low");
  //cq_IMnpipi_woK0_wSid_n_Sp_side_low->cd();
  //q_IMnpipi_woK0_wSid_n_Sp_side[sidebandtype][LOWside][0]->Draw("colz");

  //TCanvas *cq_IMnpipi_wSid_n_Sp_side_high = new TCanvas("cq_IMnpipi_wSid_n_Sp_side_high","q_IMnpipi_wSid_n_Sp_side_high");
  //cq_IMnpipi_wSid_n_Sp_side_high->cd();
  //q_IMnpipi_wSid_n_Sp_side[sidebandtype][HIGHside][0]->Draw("colz");
  //TCanvas *cq_IMnpipi_woK0_wSid_n_Sp_side_high = new TCanvas("cq_IMnpipi_woK0_wSid_n_Sp_side_high","q_IMnpipi_woK0_wSid_n_Sp_side_high");
  //cq_IMnpipi_woK0_wSid_n_Sp_side_high->cd();
  //q_IMnpipi_woK0_wSid_n_Sp_side[sidebandtype][HIGHside][0]->Draw("colz");



  //TCanvas *cq_IMnpipi_wSid_n_Sm_side_low = new TCanvas("cq_IMnpipi_wSid_n_Sm_side_low","q_IMnpipi_wSid_n_Sm_side_low");
  //cq_IMnpipi_wSid_n_Sm_side_low->cd();
  //q_IMnpipi_wSid_n_Sm_side[sidebandtype][LOWside][0]->Draw("colz");
  //TCanvas *cq_IMnpipi_woK0_wSid_n_Sm_side_low = new TCanvas("cq_IMnpipi_woK0_wSid_n_Sm_side_low","q_IMnpipi_woK0_wSid_n_Sm_side_low");
  //cq_IMnpipi_woK0_wSid_n_Sm_side_low->cd();
  //q_IMnpipi_woK0_wSid_n_Sm_side[sidebandtype][LOWside][0]->Draw("colz");

  //TCanvas *cq_IMnpipi_wSid_n_Sm_side_high = new TCanvas("cq_IMnpipi_wSid_n_Sm_side_high","q_IMnpipi_wSid_n_Sm_side_high");
  //cq_IMnpipi_wSid_n_Sm_side_high->cd();
  //q_IMnpipi_wSid_n_Sm_side[sidebandtype][HIGHside][0]->Draw("colz");
  //TCanvas *cq_IMnpipi_woK0_wSid_n_Sm_side_high = new TCanvas("cq_IMnpipi_woK0_wSid_n_Sm_side_high","q_IMnpipi_woK0_wSid_n_Sm_side_high");
  //cq_IMnpipi_woK0_wSid_n_Sm_side_high->cd();
  //q_IMnpipi_woK0_wSid_n_Sm_side[sidebandtype][HIGHside][0]->Draw("colz");

  TCanvas *cq_IMnpipi_wSid_n_px_SpSm = new TCanvas("cq_IMnpipi_wSid_n_px_SpSm","q_IMnpipi_wSid_n_px_SpSm");
  cq_IMnpipi_wSid_n_px_SpSm->cd();
  TH1D *q_IMnpipi_wSid_n_px = q_IMnpipi_wSid_n->ProjectionX();
  q_IMnpipi_wSid_n_px->Draw("EH");
  TH1D* q_IMnpipi_wSid_n_Sp_px = (TH1D*) q_IMnpipi_wSid_n_Sp->ProjectionX();
  TH1D* q_IMnpipi_wSid_n_Sm_px = (TH1D*) q_IMnpipi_wSid_n_Sm->ProjectionX();
  q_IMnpipi_wSid_n_Sp_px->SetLineColor(2);
  q_IMnpipi_wSid_n_Sm_px->SetLineColor(3);
  q_IMnpipi_wSid_n_Sp_px->Draw("HEsame");
  q_IMnpipi_wSid_n_Sm_px->Draw("HEsame");

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

  
  /*
  TCanvas *cq_IMnpipi_wSid_n_px_Sp[ngap];
  TH1D* q_IMnpipi_wSid_n_Sp_side_low_px[ngap];
  TH1D* q_IMnpipi_wSid_n_Sp_side_high_px[ngap];
  TH1D* q_IMnpipi_wSid_nSp_side_px_sum[ngap];
  for(int igap=0; igap<ngap; igap++) {
    cq_IMnpipi_wSid_n_px_Sp[igap]  = new TCanvas(Form("cq_IMnpipi_wSid_n_px_Sp_%d",igap),Form("q_IMnpipi_wSid_n_px_Sp_%d",igap));
    cq_IMnpipi_wSid_n_px_Sp[igap]->Divide(2,2);
    cq_IMnpipi_wSid_n_px_Sp[igap]->cd(1);
    IMnpim_IMnpip_dE_n_Sp->SetMaximum(IMnpim_IMnpip_dE_n->GetMaximum());
    IMnpim_IMnpip_dE_n_Sp->Draw("colz");
    //IMnpim_IMnpip_dE_n_Sp_side[LOWside][igap]->Draw("colsame");
    //IMnpim_IMnpip_dE_n_Sp_side[HIGHside][igap]->Draw("colsame");

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
  }*/
  
  /*
  TCanvas *cq_IMnpipi_woK0_wSid_n_px_Sp[ngap];
  TH1D* q_IMnpipi_woK0_wSid_n_Sp_side_low_px[ngap];
  TH1D* q_IMnpipi_woK0_wSid_n_Sp_side_high_px[ngap];
  TH1D* q_IMnpipi_woK0_wSid_nSp_side_px_sum[ngap];
  for(int igap=0; igap<ngap; igap++) {
    cq_IMnpipi_woK0_wSid_n_px_Sp[igap]  = new TCanvas(Form("cq_IMnpipi_woK0_wSid_n_px_Sp_%d",igap),Form("q_IMnpipi_woK0_wSid_n_px_Sp_%d",igap));
    cq_IMnpipi_woK0_wSid_n_px_Sp[igap]->Divide(2,2);
    cq_IMnpipi_woK0_wSid_n_px_Sp[igap]->cd(1);
    IMnpim_IMnpip_dE_woK0_n_Sp->SetMaximum(IMnpim_IMnpip_dE_woK0_n->GetMaximum());
    IMnpim_IMnpip_dE_woK0_n_Sp->Draw("colz");
    //IMnpim_IMnpip_dE_woK0_n_Sp_side[LOWside][igap]->Draw("colsame");
    //IMnpim_IMnpip_dE_woK0_n_Sp_side[HIGHside][igap]->Draw("colsame");

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
  */

  /*
  //sideband subtracted spectra for Sp mode
  TCanvas *cq_IMnpipi_woK0_wSid_n_Sp_sub = new TCanvas("cq_IMnpipi_woK0_wSid_n_Sp_sub","q_IMnpipi_woK0_wSid_n_Sp_sub");
  cq_IMnpipi_woK0_wSid_n_Sp_sub->cd();
  TH2F* q_IMnpipi_woK0_wSid_n_Sp_sub = (TH2F*)q_IMnpipi_woK0_wSid_n_Sp->Clone("Sp_sub");
  q_IMnpipi_woK0_wSid_n_Sp_sub->Add(q_IMnpipi_woK0_wSid_n_Sp_side[sidebandtype][LOWside][0],-1);
  q_IMnpipi_woK0_wSid_n_Sp_sub->Add(q_IMnpipi_woK0_wSid_n_Sp_side[sidebandtype][HIGHside][0],-1);
  q_IMnpipi_woK0_wSid_n_Sp_sub->SetMinimum(0);
  q_IMnpipi_woK0_wSid_n_Sp_sub->Draw("colz");
  */ 

  /*
  TCanvas *cq_IMnpipi_wSid_n_Sm_px[ngap];
  TH1D* q_IMnpipi_wSid_n_Sm_side_low_px[ngap];
  TH1D* q_IMnpipi_wSid_n_Sm_side_high_px[ngap];
  TH1D* q_IMnpipi_wSid_n_Sm_side_px_sum[ngap];
  for(int igap=0; igap<ngap; igap++) {
    cq_IMnpipi_wSid_n_Sm_px[igap] = new TCanvas(Form("cq_IMnpipi_wSid_n_Sm_px_%d",igap),Form("q_IMnpipi_wSid_n_Sm_px_%d",igap));
    cq_IMnpipi_wSid_n_Sm_px[igap]->Divide(2,2);
    cq_IMnpipi_wSid_n_Sm_px[igap]->cd(1);
    IMnpim_IMnpip_dE_n_Sm->SetMaximum(IMnpim_IMnpip_dE_n->GetMaximum());
    IMnpim_IMnpip_dE_n_Sm->GetXaxis()->SetRangeUser(0,1.7);
    IMnpim_IMnpip_dE_n_Sm->GetYaxis()->SetRangeUser(0,1.7);
    IMnpim_IMnpip_dE_n_Sm->Draw("colz");
    //IMnpim_IMnpip_dE_n_Sm_side[LOWside][igap]->Draw("colsame");
    //IMnpim_IMnpip_dE_n_Sm_side[HIGHside][igap]->Draw("colsame");

    cq_IMnpipi_wSid_n_Sm_px[igap]->cd(2);
    //q_IMnpipi_woK0_wSid_n_Sm_px->SetMaximum(q_IMnpipi_woK0_wSid_n_px->GetMaximum());
    if(qvalcutflag==1) q_IMnpipi_wSid_n_Sm_px->SetMaximum(160);
    if(qvalcutflag==2) q_IMnpipi_wSid_n_Sm_px->SetMaximum(260);
    q_IMnpipi_wSid_n_Sm_px->Draw("EH");
    //q_IMnpipi_wSid_n_Sm_side_low_px[igap] = (TH1D*) q_IMnpipi_wSid_n_Sm_side[sidebandtype][LOWside][igap]->ProjectionX();
    //q_IMnpipi_wSid_n_Sm_side_low_px[igap]->SetLineColor(kPink+1);
    //q_IMnpipi_wSid_n_Sm_side_low_px[igap]->SetMarkerStyle(20);
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
  }*/
  
  /*
  //overlay signal and sideband projected 1d histo to IM
  TCanvas *cq_IMnpipi_woK0_wSid_n_Sm_px[ngap];
  TH1D* q_IMnpipi_woK0_wSid_n_Sm_side_low_px[ngap];
  TH1D* q_IMnpipi_woK0_wSid_n_Sm_side_high_px[ngap];
  TH1D* q_IMnpipi_woK0_wSid_n_Sm_side_px_sum[ngap];
  for(int igap=0; igap<ngap; igap++) {
    cq_IMnpipi_woK0_wSid_n_Sm_px[igap] = new TCanvas(Form("cq_IMnpipi_woK0_wSid_n_Sm_px_%d",igap),Form("q_IMnpipi_woK0_wSid_n_Sm_px_%d",igap));
    cq_IMnpipi_woK0_wSid_n_Sm_px[igap]->Divide(2,2);
    cq_IMnpipi_woK0_wSid_n_Sm_px[igap]->cd(1);
    IMnpim_IMnpip_dE_woK0_n_Sm->SetMaximum(IMnpim_IMnpip_dE_woK0_n->GetMaximum());
    IMnpim_IMnpip_dE_woK0_n_Sm->GetXaxis()->SetRangeUser(0,1.7);
    IMnpim_IMnpip_dE_woK0_n_Sm->GetYaxis()->SetRangeUser(0,1.7);
    IMnpim_IMnpip_dE_woK0_n_Sm->Draw("colz");
    //IMnpim_IMnpip_dE_woK0_n_Sm_side[LOWside][igap]->Draw("colsame");
    //IMnpim_IMnpip_dE_woK0_n_Sm_side[HIGHside][igap]->Draw("colsame");

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
  */

  /*
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
  */

  //w/ K0
  TCanvas *cq_IMnpipi_wSid_n = new TCanvas("cq_IMnpipi_wSid_n","q_IMnpipi_wSid_n");
  cq_IMnpipi_wSid_n->cd();
  q_IMnpipi_wSid_n->Draw("colz");
  
  /*
  TCanvas *cq_IMnpipi_wSid_n_side[ngap];
  for(int igap=0; igap<ngap; igap++) {
    cq_IMnpipi_wSid_n_side[igap] = new TCanvas(Form("cq_IMnpipi_wSid_n_side_%d",igap),Form("q_IMnpipi_wSid_n_side_%d",igap));
    cq_IMnpipi_wSid_n_side[igap]->cd();
    q_IMnpipi_wSid_n_side[igap]->SetMaximum(q_IMnpipi_wSid_n->GetMaximum());
    q_IMnpipi_wSid_n_side[igap]->Draw("colz");
  }
  */

  TCanvas *cq_IMnpipi_woK0_wSid_n = new TCanvas("cq_IMnpipi_woK0_wSid_n","q_IMnpipi_woK0_wSid_n");
  cq_IMnpipi_woK0_wSid_n->cd();
  q_IMnpipi_woK0_wSid_n->Draw("colz");
  fkp->Draw("same");
  //q_IMnpipi_woK0_wSid_n->Draw("colzsame");
  //fnu->cd();
  //TMultiGraph *mg = (TMultiGraph*)fnu->Get("mg");
  //mg->Draw("c");
  f->cd();

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

  //TCanvas *cdE_nmom_fid_beta_wK0 = new TCanvas("cdE_nmom_fid_beta_wK0","dE_nmom_fid_beta_wK0");
  //cdE_nmom_fid_beta_wK0->cd();
  //dE_nmom_fid_beta_wK0->Draw("colz");
  //cdE_nmom_fid_beta_wK0->SetLogz();

  TCanvas *cdE_MMom_fid_beta = new TCanvas("cdE_MMom_fid_beta","dE_MMom_fid_beta");
  cdE_MMom_fid_beta->cd();
  dE_MMom_fid_beta->Draw("colz");

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
  TH1D* nmom_IMnpipi_wK0_n_px = nmom_IMnpipi_wK0_n->ProjectionX("nmom_IMnpipi_wK0_n_px");
  nmom_IMnpipi_wK0_n_px->SetFillColor(kViolet);
  nmom_IMnpipi_wK0_n_px->Draw("HE");

  TCanvas *cnmom_IMnpipi_wK0_n_py = new TCanvas(Form("cnmom_IMnpipi_wK0_n_py"),"nmom_IMnpipi_wK0_n_py");
  cnmom_IMnpipi_wK0_n_py->cd();
  TH1D* nmom_IMnpipi_wK0_n_py = nmom_IMnpipi_wK0_n->ProjectionY("nmom_IMnpipi_wK0_n_py");
  nmom_IMnpipi_wK0_n_py->SetFillColor(kViolet);
  //nmom_IMnpipi_wK0_n_py->SetMaximum(3300);
  nmom_IMnpipi_wK0_n_py->Draw("HE");

  //TCanvas *cnmom_IMnpipi_woK0_wSid_n = new TCanvas(Form("cnmom_IMnpipi_woK0_wSid_n"),"nmom_IMnpipi_woK0_wSid_n");
  //cnmom_IMnpipi_woK0_wSid_n->cd();
  //nmom_IMnpipi_woK0_wSid_n->Draw("colz");

  //TCanvas *cnmom_IMnpipi_woK0_wSid_n_px = new TCanvas(Form("cnmom_IMnpipi_woK0_wSid_n_px"),"nmom_IMnpipi_woK0_wSid_n_px");
  //cnmom_IMnpipi_woK0_wSid_n_px->cd();
  //TH1D* nmom_IMnpipi_woK0_wSid_n_px = nmom_IMnpipi_woK0_wSid_n->ProjectionX("nmom_IMnpipi_woK0_wSid_n_px");
  //nmom_IMnpipi_woK0_wSid_n_px->Draw();

  //TCanvas *cnmom_IMnpipi_woK0_wSid_n_py = new TCanvas(Form("cnmom_IMnpipi_woK0_wSid_n_py"),"nmom_IMnpipi_woK0_wSid_n_py");
  //cnmom_IMnpipi_woK0_wSid_n_py->cd();
  //TH1D* nmom_IMnpipi_woK0_wSid_n_py = nmom_IMnpipi_woK0_wSid_n->ProjectionY("nmom_IMnpipi_woK0_wSid_n_py");
  //nmom_IMnpipi_woK0_wSid_n_py->Draw();


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

  TCanvas *cpipmom_IMnpipi_wSid_n = new TCanvas("cpipmom_IMnpipi_wSid_n","pipmom_IMnpipi_wSid_n");
  cpipmom_IMnpipi_wSid_n->cd();
  pipmom_IMnpipi_wSid_n->Draw("colz");

  TCanvas *cpimmom_IMnpipi_wSid_n = new TCanvas("cpimmom_IMnpipi_wSid_n","pimmom_IMnpipi_wSid_n");
  cpimmom_IMnpipi_wSid_n->cd();
  pimmom_IMnpipi_wSid_n->Draw("colz");

  TCanvas *cpipmom_IMnpipi_woK0_wSid_n = new TCanvas("cpipmom_IMnpipi_woK0_wSid_n","pipmom_IMnpipi_woK0_wSid_n");
  cpipmom_IMnpipi_woK0_wSid_n->cd();
  pipmom_IMnpipi_woK0_wSid_n->Draw("colz");

  TCanvas *cpimmom_IMnpipi_woK0_wSid_n = new TCanvas("cpimmom_IMnpipi_woK0_wSid_n","pimmom_IMnpipi_woK0_wSid_n");
  cpimmom_IMnpipi_woK0_wSid_n->cd();
  pimmom_IMnpipi_woK0_wSid_n->Draw("colz");

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
  TH1D* mnmom_IMpippim_n_py = (TH1D*)mnmom_IMpippim_n->ProjectionY("pytest");
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


  TCanvas *cMMnmiss_IMnpipi_wSid = new TCanvas("cMMnmiss_IMnpipi_wSid","MMnmiss_IMnpipi_wSid");
  cMMnmiss_IMnpipi_wSid->cd();
  MMnmiss_IMnpipi_wSid->Draw("colz");

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

  TCanvas *cq_IMpippim_wSid_n = new TCanvas("cq_IMpippim_wSid_n","q_IMpippim_wSid_n");
  cq_IMpippim_wSid_n->cd();
  q_IMpippim_wSid_n->Draw("colz");

  TCanvas *cMompippim_IMnpipi_dE_wK0_n = new TCanvas("cMompippim_IMnpipi_dE_wK0_n","Mompippim_IMnpipi_dE_wK0_n");
  cMompippim_IMnpipi_dE_wK0_n->cd();
  Mompippim_IMnpipi_dE_wK0_n->Draw("colz");
  fkp->Draw("same");

  TCanvas *cMompippim_nmom_dE_wK0_n = new TCanvas("cMompippim_nmom_dE_wK0_n","Mompippim_nmom_dE_wK0_n");
  cMompippim_nmom_dE_wK0_n->cd();
  Mompippim_nmom_dE_wK0_n->Draw("colz");

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

  TCanvas *cMMnpim_MMnpip_woK0_n = new TCanvas("cMMnpim_MMnpip_woK0_n","MMnpim_MMnpip_woK0_n");
  cMMnpim_MMnpip_woK0_n->cd();
  MMnpim_MMnpip_woK0_n->Draw("colz");

  TCanvas *cMMnpip_MMnpim_wwoK0_n_px = new TCanvas("cMMnpip_MMnpim_wwoK0_n_px","cMMnpip_MMnpim_wwoK0_n_px");
  cMMnpip_MMnpim_wwoK0_n_px->cd();
  TH1D* MMnpim_MMnpip_n_px = (TH1D*)MMnpim_MMnpip_n->ProjectionX();
  MMnpim_MMnpip_n_px->Draw("HE");
  TH1D* MMnpim_MMnpip_woK0_n_px = (TH1D*)MMnpim_MMnpip_woK0_n->ProjectionX();
  MMnpim_MMnpip_woK0_n_px->SetLineColor(kViolet);
  MMnpim_MMnpip_woK0_n_px->Draw("HEsame");

  TCanvas *cMMnpip_MMnpim_wwoK0_n_py = new TCanvas("cMMnpip_MMnpim_wwoK0_n_py","cMMnpip_MMnpim_wwoK0_n_py");
  cMMnpip_MMnpim_wwoK0_n_py->cd();
  TH1D* MMnpim_MMnpip_n_py = (TH1D*)MMnpim_MMnpip_n->ProjectionY();
  MMnpim_MMnpip_n_py->Draw("HE");
  TH1D* MMnpim_MMnpip_woK0_n_py = (TH1D*)MMnpim_MMnpip_woK0_n->ProjectionY();
  MMnpim_MMnpip_woK0_n_py->SetLineColor(kViolet);
  MMnpim_MMnpip_woK0_n_py->Draw("HEsame");


  TCanvas *cMMnpim_MMnpip_woK0_wSid_n = new TCanvas("cMMnpim_MMnpip_woK0_wSid_n","MMnpim_MMnpip_woK0_wSid_n");
  cMMnpim_MMnpip_woK0_wSid_n->cd();
  MMnpim_MMnpip_woK0_wSid_n->Draw("colz");

  TCanvas *cMMnpim_MMnpip_woK0_wSid_n_px = new TCanvas("cMMnpim_MMnpip_woK0_wSid_n_px","cMMnpim_MMnpip_woK0_wSid_n_px");
  cMMnpim_MMnpip_woK0_wSid_n_px->cd();
  TH1D* MMnpim_MMnpip_woK0_wSid_n_px = (TH1D*)MMnpim_MMnpip_woK0_wSid_n->ProjectionX();
  MMnpim_MMnpip_woK0_wSid_n_px->Draw("HE");

  TCanvas *cMMnpim_MMnpip_woK0_wSid_n_py = new TCanvas("cMMnpim_MMnpip_woK0_wSid_n_py","cMMnpim_MMnpip_woK0_wSid_n_py");
  cMMnpim_MMnpip_woK0_wSid_n_py->cd();
  TH1D* MMnpim_MMnpip_woK0_wSid_n_py = (TH1D*)MMnpim_MMnpip_woK0_wSid_n->ProjectionY();
  MMnpim_MMnpip_woK0_wSid_n_py->Draw("HE");


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

  /*
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
  */


  //TCanvas *cMMom_MMass_woK0_px_sup = new TCanvas("cMMom_MMass_woK0_px_sup","MMom_MMass_woK0_px_sup");
  //MMom_MMass_woK0_px->Draw();
  //TH1D *MMom_MMass_woK0_wSid_px_clone = (TH1D*)MMom_MMass_woK0_wSid_px->Clone();
  //MMom_MMass_woK0_wSid_px_clone->SetLineColor(4);
  //cMMom_MMass_woK0_px_sup->cd();
  //MMom_MMass_woK0_wSid_px_clone->SetTitle("Missing Mass d(K^{-},#pi^{+}#pi^{-}n)\"X\"");
  //MMom_MMass_woK0_wSid_px_clone->Draw("");
  /*
  TCanvas *cMMom_MMass_woK0_py = new TCanvas("cMMom_MMass_woK0_py","MMom_MMass_woK0_py");
  TH1D *MMom_MMass_woK0_py = MMom_MMass_woK0->ProjectionY();
  TH1D *MMom_MMass_woK0_wSid_py = MMom_MMass_woK0_wSid->ProjectionY();
  MMom_MMass_woK0_py->Draw();
  TH1D *MMom_MMass_woK0_wSid_py_clone = (TH1D*)MMom_MMass_woK0_wSid_py->Clone();
  MMom_MMass_woK0_wSid_py_clone->SetLineColor(4);
  cMMom_MMass_woK0_py->cd();
  MMom_MMass_woK0_wSid_py_clone->SetTitle("Missing Mom. d(K^{-},#pi^{+}#pi^{-}n)\"X\"");
  MMom_MMass_woK0_wSid_py_clone->Draw("");
  */
  //acceptance calculation
  TFile *facc_Sp = new TFile("../simpost/acc_Sp.root","READ");
  TFile *facc_Sm = new TFile("../simpost/acc_Sm.root","READ");

  if(SimSpmode || SimSmmode) {
    if(SimSpmode) {
      std::cout << "This is Sigma+ mode sim." << std::endl;
    } else {
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
  TH2F* acc_Sp_woK0 = (TH2F*)facc_Sp->Get("eff_q_IMpiSigma_woK0_wSid_n_reco");
  //w/o K0, Sp selection
  TH2F* acc_SpSel_woK0 = (TH2F*)facc_Sp->Get("eff_q_IMpiSigma_woK0_wSid_n_SpSm_reco");
  //TH2F* acc_SpSel_woK0 = (TH2F*)facc_Sp->Get("eff_q_IMpiSigma_woK0_wSid_n_SpSm");
  TH2F* acc_Sp_err = (TH2F*)facc_Sp->Get("heff_err");
  if(acc_Sp == NULL) {
    std::cout << " acc_Sp is NULL " << std::endl;
    return;
  }
  if(acc_Sp_woK0 == NULL) {
    std::cout << " acc_Sp_woK0 is NULL " << std::endl;
    //return;
  }
  if(acc_SpSel_woK0 == NULL) {
    std::cout << " acc_SpSel_woK0 is NULL " << std::endl;
    return;
  }
  if(acc_Sp_err == NULL) {
    std::cout << " acc_Sp_err is NULL " << std::endl;
    return;
  }

  //TH2F* acc_Sm = (TH2F*)facc_Sm->Get("eff_q_IMpiSigma_woK0_wSid_n_SpSm");
  TH2F* acc_Sm = (TH2F*)facc_Sm->Get("eff_q_IMpiSigma_wSid_n_reco");
  TH2F* acc_Sm_woK0 = (TH2F*)facc_Sm->Get("eff_q_IMpiSigma_woK0_wSid_n_reco");
  TH2F* acc_SmSel_woK0 = (TH2F*)facc_Sm->Get("eff_q_IMpiSigma_woK0_wSid_n_SpSm_reco");
  TH2F* acc_Sm_err = (TH2F*)facc_Sm->Get("heff_err");
  if(acc_Sm == NULL) {
    std::cout << " acc_Sm is NULL " << std::endl;
    return;
  }
  if(acc_Sm_woK0 == NULL) {
    std::cout << " acc_Sm_woK0 is NULL " << std::endl;
    return;
  }
  if(acc_SmSel_woK0 == NULL) {
    std::cout << " acc_SmSel_woK0 is NULL " << std::endl;
    return;
  }
  if(acc_Sm_err == NULL) {
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
  for(int ibinx=0; ibinx<acc_Sp_err->GetNbinsX(); ibinx++) {
    for(int ibiny=0; ibiny<acc_Sp_err->GetNbinsY(); ibiny++) {
      double precision =  acc_Sp_err->GetBinContent(ibinx,ibiny);
      if(precision>0.20) {
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
  for(int ibinx=0; ibinx<acc_Sp_err->GetNbinsX(); ibinx++) {
    for(int ibiny=0; ibiny<acc_Sp_err->GetNbinsY(); ibiny++) {
      double precision =  acc_Sp_err->GetBinContent(ibinx,ibiny);
      if(precision>0.20) {
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
  
  /*
  TCanvas *cq_IMnpipi_woK0_wSid_n_Sp_side_cs[ngap];
  TH2F *q_IMnpipi_woK0_wSid_n_Sp_side_cs_low[ngap];
  TH2F *q_IMnpipi_woK0_wSid_n_Sp_side_cs_high[ngap];
  TH2F *q_IMnpipi_woK0_wSid_n_Sp_side_cs[ngap];
  for(int igap=0; igap<ngap; igap++) {
    cq_IMnpipi_woK0_wSid_n_Sp_side_cs[igap] = new TCanvas(Form("cq_IMnpipi_woK0_wSid_n_Sp_side_cs_%d",igap),Form("q_IMnpipi_woK0_wSid_n_Sp_side_cs_%d",igap));
    cq_IMnpipi_woK0_wSid_n_Sp_side_cs[igap]->cd();
    q_IMnpipi_woK0_wSid_n_Sp_side_cs_low[igap] = (TH2F*)q_IMnpipi_woK0_wSid_n_Sp_side[sidebandtype][LOWside][igap]->Clone(Form("Sp_cs_side_low_%d",igap));
    q_IMnpipi_woK0_wSid_n_Sp_side_cs_high[igap] = (TH2F*)q_IMnpipi_woK0_wSid_n_Sp_side[sidebandtype][HIGHside][igap]->Clone(Form("Sp_cs_side_high_%d",igap));
    q_IMnpipi_woK0_wSid_n_Sp_side_cs_low[igap]->Sumw2();
    q_IMnpipi_woK0_wSid_n_Sp_side_cs_high[igap]->Sumw2();
    q_IMnpipi_woK0_wSid_n_Sp_side_cs_low[igap]->Divide(acc_SpSel_woK0);
    q_IMnpipi_woK0_wSid_n_Sp_side_cs_high[igap]->Divide(acc_SpSel_woK0);
    if(!(SimSpmode || SimSmmode)) {
      q_IMnpipi_woK0_wSid_n_Sp_side_cs_low[igap]->Scale(1./efflumi);
      q_IMnpipi_woK0_wSid_n_Sp_side_cs_high[igap]->Scale(1./efflumi);
    }
    q_IMnpipi_woK0_wSid_n_Sp_side_cs[igap] = (TH2F*)q_IMnpipi_woK0_wSid_n_Sp_side[sidebandtype][LOWside][igap]->Clone(Form("Sp_cs_side_sum_%d",igap));
    q_IMnpipi_woK0_wSid_n_Sp_side_cs[igap]->Add(q_IMnpipi_woK0_wSid_n_Sp_side[sidebandtype][HIGHside][igap]);
    q_IMnpipi_woK0_wSid_n_Sp_side_cs[igap]->Sumw2();
    q_IMnpipi_woK0_wSid_n_Sp_side_cs[igap]->Divide(acc_SpSel_woK0);
    if(!(SimSpmode || SimSmmode)) {
      q_IMnpipi_woK0_wSid_n_Sp_side_cs[igap]->Scale(1./efflumi);
    }
    q_IMnpipi_woK0_wSid_n_Sp_side_cs[igap]->SetMaximum(q_IMnpipi_woK0_wSid_n_Sp_cs->GetMaximum());
    q_IMnpipi_woK0_wSid_n_Sp_side_cs[igap]->Draw("colz");
  }

  TCanvas *cq_IMnpipi_woK0_wSid_n_Sp_side_cs_px[ngap];
  TH1D* q_IMnpipi_woK0_wSid_n_Sp_side_cs_low_px[ngap];
  TH1D* q_IMnpipi_woK0_wSid_n_Sp_side_cs_high_px[ngap];
  TH1D* q_IMnpipi_woK0_wSid_n_Sp_side_cs_px[ngap];
  for(int igap=0; igap<ngap; igap++) {
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
  */

  //Sm mode Cross section
  TCanvas *cq_IMnpipi_woK0_wSid_n_Sm_cs = new TCanvas("cq_IMnpipi_woK0_wSid_n_Sm_cs","q_IMnpipi_woK0_wSid_n_Sm_cs");
  cq_IMnpipi_woK0_wSid_n_Sm_cs->cd();
  TH2F *q_IMnpipi_woK0_wSid_n_Sm_cs = (TH2F*)q_IMnpipi_woK0_wSid_n_Sm->Clone("Sm_cs");
  q_IMnpipi_woK0_wSid_n_Sm_cs->Sumw2();
  acc_Sm->Print("base");
  //cleaning up
  for(int ibinx=0; ibinx<acc_Sm_err->GetNbinsX(); ibinx++) {
    for(int ibiny=0; ibiny<acc_Sm_err->GetNbinsY(); ibiny++) {
      double precision =  acc_Sm_err->GetBinContent(ibinx,ibiny);
      //double cont = q_IMnpipi_woK0_wSid_n_Sm_cs->GetBinContent(ibinx,ibiny);
      if(precision>0.20 ) {
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
  
  /*
  TCanvas *cq_IMnpipi_woK0_wSid_n_Sm_side_cs[ngap];
  TH2F *q_IMnpipi_woK0_wSid_n_Sm_side_cs_low[ngap];
  TH2F *q_IMnpipi_woK0_wSid_n_Sm_side_cs_high[ngap];
  TH2F *q_IMnpipi_woK0_wSid_n_Sm_side_cs[ngap];
  for(int igap=0; igap<ngap; igap++) {
    cq_IMnpipi_woK0_wSid_n_Sm_side_cs[igap] = new TCanvas(Form("cq_IMnpipi_woK0_wSid_n_Sm_side_cs_%d",igap),Form("q_IMnpipi_woK0_wSid_n_Sm_side_cs_%d",igap));
    cq_IMnpipi_woK0_wSid_n_Sm_side_cs[igap]->cd();
    q_IMnpipi_woK0_wSid_n_Sm_side_cs_low[igap] = (TH2F*)q_IMnpipi_woK0_wSid_n_Sm_side[sidebandtype][LOWside][igap]->Clone(Form("Sm_cs_side_low_%d",igap));
    q_IMnpipi_woK0_wSid_n_Sm_side_cs_high[igap] = (TH2F*)q_IMnpipi_woK0_wSid_n_Sm_side[sidebandtype][HIGHside][igap]->Clone(Form("Sm_cs_side_high_%d",igap));
    //q_IMnpipi_woK0_wSid_n_Sm_side_cs->Add(q_IMnpipi_woK0_wSid_n_Sm_side[sidebandtype][HIGHside][0]);
    q_IMnpipi_woK0_wSid_n_Sm_side_cs_low[igap]->Sumw2();
    q_IMnpipi_woK0_wSid_n_Sm_side_cs_high[igap]->Sumw2();
    q_IMnpipi_woK0_wSid_n_Sm_side_cs_low[igap]->Divide(acc_SmSel_woK0);
    q_IMnpipi_woK0_wSid_n_Sm_side_cs_high[igap]->Divide(acc_SmSel_woK0);
    if(!(SimSpmode || SimSmmode)) {
      q_IMnpipi_woK0_wSid_n_Sm_side_cs_low[igap]->Scale(1.0/efflumi);
      q_IMnpipi_woK0_wSid_n_Sm_side_cs_high[igap]->Scale(1.0/efflumi);
    }
    q_IMnpipi_woK0_wSid_n_Sm_side_cs[igap] = (TH2F*)q_IMnpipi_woK0_wSid_n_Sm_side[sidebandtype][LOWside][igap]->Clone(Form("Sm_cs_side_sum_%d",igap));
    q_IMnpipi_woK0_wSid_n_Sm_side_cs[igap]->Add(q_IMnpipi_woK0_wSid_n_Sm_side[sidebandtype][HIGHside][igap]);
    q_IMnpipi_woK0_wSid_n_Sm_side_cs[igap]->Sumw2();
    q_IMnpipi_woK0_wSid_n_Sm_side_cs[igap]->Divide(acc_SmSel_woK0);
    if(!(SimSpmode || SimSmmode)) {
      q_IMnpipi_woK0_wSid_n_Sm_side_cs[igap]->Scale(1./efflumi);
    }
    q_IMnpipi_woK0_wSid_n_Sm_side_cs[igap]->SetMaximum(q_IMnpipi_woK0_wSid_n_Sm_cs->GetMaximum());
    q_IMnpipi_woK0_wSid_n_Sm_side_cs[igap]->Draw("colz");
  }

  TCanvas *cq_IMnpipi_woK0_wSid_n_Sm_side_cs_px[ngap];
  TH1D* q_IMnpipi_woK0_wSid_n_Sm_side_cs_low_px[ngap];
  TH1D* q_IMnpipi_woK0_wSid_n_Sm_side_cs_high_px[ngap];
  TH1D* q_IMnpipi_woK0_wSid_n_Sm_side_cs_px[ngap];
  for(int igap=0; igap<ngap; igap++) {
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
  }*/
  
  /*
  TCanvas *cIMnpim_IMnpip_dE_woK0_n_Sp_wide_cs[nzone];
  //TH2D* q_IMnpipi_woK0_wSid_n_Sp_sidewide_cs[nzone];
  //TH1D* q_IMnpipi_woK0_wSid_n_Sp_sidewide_cs_px[nzone];
  for(int izone=0; izone<nzone; izone++) {
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
    if(!(SimSpmode || SimSmmode)) {
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
  */

  /*
  TCanvas *cIMnpim_IMnpip_dE_woK0_n_Sm_wide_cs[nzone];
  TH2D* q_IMnpipi_woK0_wSid_n_Sm_sidewide_cs[nzone];
  TH1D* q_IMnpipi_woK0_wSid_n_Sm_sidewide_cs_px[nzone];
  for(int izone=0; izone<nzone; izone++) {
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
    if(!(SimSpmode || SimSmmode)) {
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
  }*/


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
    else if(SimK0nnmode) {
      pt->SetFillColor(kAzure-4);
      pt->AddText("MC K0nn");
    }
    else if(SimK0n_nsmode) {
      pt->SetFillColor(kAzure-4);
      pt->AddText("MC K0n_ns");
    }
    else if(SimK0_nntsmode) {
      pt->SetFillColor(kAzure-4);
      pt->AddText("MC K0_nnts");
    }
    else if(SimnpipiLmode) {
      pt->SetFillColor(kAzure-4);
      pt->AddText("MC n#pi+#pi-#Lambda");
    }
    else if(SimnS0pippimmode) {
      pt->SetFillColor(kAzure-4);
      pt->AddText("MC n#Sigma0#pi#pi-");
    }
    else if(SimSppi0mode) {
      pt->SetFillColor(kAzure-4);
      pt->AddText("MC n#Sigma+#pi+#pi-#pi0");
    }
    else if(SimSmpi0mode) {
      pt->SetFillColor(kAzure-4);
      pt->AddText("MC n#Sigma-#pi+#pi-#pi0");
    } else if(SimSp_nsmode) {
      pt->SetFillColor(kAzure-4);
      pt->AddText("MC #Sigma+#pi+#pi- 1NA");
    } else if(SimSm_nsmode) {
      pt->SetFillColor(kAzure-4);
      pt->AddText("MC #Sigma-#pi+#pi- 1NA");
    } else if(SimFakemode) {
      pt->SetFillColor(kAzure-4);
      pt->AddText("MC Fake");
    } else if(SimFakeK0mode) {
      pt->SetFillColor(kAzure-4);
      pt->AddText("MC Fake");
    }
    else {
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
  if(IsolationFlag==0) outname.Replace(std::string(filename).size()-5,5,"_out_noiso.root");
  else if(IsolationFlag==1) outname.Replace(std::string(filename).size()-5,5,"_out_iso.root");
  else if(IsolationFlag==2) outname.Replace(std::string(filename).size()-5,5,"_out_isowide.root");
  else if(IsolationFlag==3) outname.Replace(std::string(filename).size()-5,5,"_out_isorev.root");
  //outname.Replace(std::string(filename).size()-5,5,"_outncutK015.root");
  
  if(qvalcutflag==1) outname.Replace(std::string(outname).size()-5,5,"_qlo.root");
  if(qvalcutflag==2) outname.Replace(std::string(outname).size()-5,5,"_qhi.root");
  if(qvalcutflag==3) outname.Replace(std::string(outname).size()-5,5,"_theta15.root");
  
  if(SimRejectFake && (SimSpmode || SimSmmode || SimK0nnmode)){
    outname.Replace(std::string(outname).size()-5,5,"_rej.root");
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

