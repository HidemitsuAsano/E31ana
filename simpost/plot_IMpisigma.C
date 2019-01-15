#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

#include <TApplication.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <TSystem.h>
#include <TBranch.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TRint.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TGraphErrors.h> 
#include <TDatabasePDG.h>
#include <TRandom.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TColor.h>
#include <TProfile.h>
#include <TFractionFitter.h>
#include <TPaletteAxis.h>

#include "../src/IMPiSigmaAnaPar.h"


const double pvalcut = 0.005;
const bool gridon=true;
const bool staton=true;

//mode 0: Sigma+ ,1: Sigma- 
void plot_IMpisigma(const char* filename="",const int mode=0)
{
  gROOT->SetStyle("Plain");
  if(staton)gStyle->SetOptStat(111111);
  else gStyle->SetOptStat(0);
  gStyle->SetOptFit(111111);
  gStyle->SetPadGridX(gridon);
  gStyle->SetPadGridY(gridon);
  
  std::string outfilename = string(filename);
  outfilename.insert(outfilename.size()-5,"_post");
  std::cout << outfilename << std::endl;
  //--- color style ---//
  
  //= = = = pipipnn final-sample tree = = = =//
  TLorentzVector *LVec_beam=nullptr;   // 4-momentum(beam)
  TLorentzVector *LVec_target=nullptr; // 4-momentum(target)
  TLorentzVector *LVec_pip=nullptr;    // 4-momentum(pi+)
  TLorentzVector *LVec_pim=nullptr;    // 4-momentum(pi-)
  TLorentzVector *LVec_n=nullptr;      // 4-momentum(neutron)
  double NeutralBetaCDH; // veracity of neutral particle on CDH
  double dE;   // energy deposit on CDH
  TVector3 *vtx_reaction = nullptr; // vertex(reaction)
  int run_num;   // run number
  int event_num; // event number
  int block_num; // block number
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
  
  tree->SetBranchAddress( "mom_target", &LVec_target );
  
  tree->SetBranchAddress( "mom_pip", &LVec_pip );
  tree->SetBranchAddress( "mom_pim", &LVec_pim );
  tree->SetBranchAddress( "mom_n", &LVec_n );
  tree->SetBranchAddress( "beta", &NeutralBetaCDH );
  tree->SetBranchAddress( "dE", &dE );
  tree->SetBranchAddress( "vtx_reaction", &vtx_reaction );
  tree->SetBranchAddress( "run_num", &run_num );
  tree->SetBranchAddress( "event_num", &event_num );
  //tree->SetBranchAddress( "block_num", &block_num );
  tree->SetBranchAddress( "kfMomBeamSpmode",   &kfSpmode_mom_beam );
  tree->SetBranchAddress( "kfMom_pip_Spmode", &kfSpmode_mom_pip );
  tree->SetBranchAddress( "kfMom_pim_Spmode", &kfSpmode_mom_pim );
  tree->SetBranchAddress( "kfMom_n_Spmode", &kfSpmode_mom_n );
  tree->SetBranchAddress( "kf_chi2_Spmode", &kfSpmode_chi2 );
  tree->SetBranchAddress( "kf_NDF_Spmode", &kfSpmode_NDF );
  tree->SetBranchAddress( "kf_status_Spmode", &kfSpmode_status );
  tree->SetBranchAddress( "kf_pvalue_Spmode", &kfSpmode_pvalue );
  tree->SetBranchAddress( "kfMomBeamSmmode",   &kfSmmode_mom_beam );
  tree->SetBranchAddress( "kfMom_pip_Smmode", &kfSmmode_mom_pip );
  tree->SetBranchAddress( "kfMom_pim_Smmode", &kfSmmode_mom_pim );
  tree->SetBranchAddress( "kfMom_n_Smmode", &kfSmmode_mom_n );
  tree->SetBranchAddress( "kf_chi2_Smmode", &kfSmmode_chi2 );
  tree->SetBranchAddress( "kf_NDF_Smmode", &kfSmmode_NDF );
  tree->SetBranchAddress( "kf_status_Smmode", &kfSmmode_status );
  tree->SetBranchAddress( "kf_pvalue_Smmode", &kfSmmode_pvalue );
  tree->SetBranchAddress( "kf_flag", &kf_flag );
  

  // w/o kinematic fit 
  TH2F* dE_betainv_fid;//
  TH2F* dE_MMom_fid_beta_woK0;
  TH2F* dE_MMass_fid_beta_woK0;
  TH2F* MMom_MMass_fid_beta_dE_woK0;
  TH2F* IMnpim_IMnpip_dE_woK0;
  TH2F* IMnpim_IMnpip_dE_woK0_n;
  TH2F* MMnpip_MMnpim_woK0_wSid_n;
  TH2F* dE_IMnpipi_woK0_wSid_n;
  TH2F* Cosn_IMnpipi_woK0_wSid_n;
  TH2F* MMnmiss_IMnpipi_woK0_wSid_n;
  TH2F* q_IMnpipi_woK0_wSid_n;
  // w/ kinematic fit
  TH2F* dE_betainv_fid_kin[2];//
  TH2F* dE_MMom_fid_beta_woK0_kin[2];
  TH2F* dE_MMass_fid_beta_woK0_kin[2];
  TH2F* MMom_MMass_fid_beta_dE_woK0_kin[2];
  TH2F* IMnpim_IMnpip_dE_woK0_kin[2];
  TH2F* IMnpim_IMnpip_dE_woK0_n_kin[2];
  TH2F* MMnpip_MMnpim_woK0_wSid_n_kin[2];
  TH2F* dE_IMnpipi_woK0_wSid_n_kin[2];
  TH2F* Cosn_IMnpipi_woK0_wSid_n_kin[2];
  TH2F* MMnmiss_IMnpipi_woK0_wSid_n_kin[2];
  TH2F* q_IMnpipi_woK0_wSid_n_kin[2];
  const char smode[][4]={"Sp","Sm"};


  dE_betainv_fid = new TH2F(Form("dE_betainv_fid"),Form("dE_betainv_fid"),1000, 0, 50, 200, 0, 50);
  dE_betainv_fid->SetXTitle("1/#beta");
  dE_betainv_fid->SetYTitle("dE [MeVee]");
  dE_betainv_fid->GetXaxis()->CenterTitle();
  dE_betainv_fid->GetYaxis()->CenterTitle();
  
  dE_MMom_fid_beta_woK0 = new TH2F(Form("dE_MMom_fid_beta_woK0"),Form("dE_MMom_fid_beta_woK0"),100, 0, 1.5, 200, 0, 50);
  dE_MMom_fid_beta_woK0->SetXTitle("Missing Mom. [GeV/c]");
  dE_MMom_fid_beta_woK0->SetYTitle("dE [MeVee]");
  dE_MMom_fid_beta_woK0->GetXaxis()->CenterTitle();
  dE_MMom_fid_beta_woK0->GetYaxis()->CenterTitle();

  dE_MMass_fid_beta_woK0 = new TH2F(Form("dE_MMass_fid_beta_woK0"),Form("dE_MMass_fid_beta_woK0"), 140, 0.4, 1.8, 200, 0, 50);
  dE_MMass_fid_beta_woK0->SetXTitle("Missing mass [GeV/c^{2}]");
  dE_MMass_fid_beta_woK0->SetYTitle("dE [MeVee]");
  dE_MMass_fid_beta_woK0->GetXaxis()->CenterTitle();
  dE_MMass_fid_beta_woK0->GetYaxis()->CenterTitle();

  MMom_MMass_fid_beta_dE_woK0 = new TH2F(Form("MMom_MMass_fid_beta_dE_woK0"),Form("MMom_MMass_fid_beta_dE_woK0"), 140, 0.4, 1.8, 100, 0, 1.5);
  MMom_MMass_fid_beta_dE_woK0->SetXTitle("Missing Mass [GeV/c^{2}]");
  MMom_MMass_fid_beta_dE_woK0->SetYTitle("Missing Mom. [GeV/c]");
  MMom_MMass_fid_beta_dE_woK0->GetXaxis()->CenterTitle();
  MMom_MMass_fid_beta_dE_woK0->GetYaxis()->CenterTitle();
  
  IMnpim_IMnpip_dE_woK0 = new TH2F(Form("IMnpim_IMnpip_dE_woK0"), Form("IMnpim_IMnpip_dE_woK0"),200, 1, 2.0, 200, 1, 2.0);
  IMnpim_IMnpip_dE_woK0->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_woK0->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_woK0->GetXaxis()->CenterTitle();
  IMnpim_IMnpip_dE_woK0->GetXaxis()->CenterTitle();
    
  IMnpim_IMnpip_dE_woK0_n = new TH2F(Form("IMnpim_IMnpip_dE_woK0_n"),Form("IMnpim_IMnpip_dE_woK0_n"),200, 1, 2.0, 200, 1, 2.0);
  IMnpim_IMnpip_dE_woK0_n->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_woK0_n->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_woK0_n->GetXaxis()->CenterTitle();
  IMnpim_IMnpip_dE_woK0_n->GetYaxis()->CenterTitle();
    
  MMnpip_MMnpim_woK0_wSid_n = new TH2F(Form("MMnpip_MMnpim_woK0_wSid_n"),Form("MMnpip_MMnpim_woK0_wSid_n"),70, 1, 1.7, 70, 1, 1.7);
  MMnpip_MMnpim_woK0_wSid_n->SetXTitle("Miss. Mass(n#pi^{+}) [GeV/c^{2}]");
  MMnpip_MMnpim_woK0_wSid_n->SetYTitle("Miss. Mass(n#pi^{-}) [GeV/c^{2}]");
  MMnpip_MMnpim_woK0_wSid_n->GetXaxis()->CenterTitle();
  MMnpip_MMnpim_woK0_wSid_n->GetYaxis()->CenterTitle();
  
  dE_IMnpipi_woK0_wSid_n = new TH2F(Form("dE_IMnpipi_woK0_wSid_n"),Form("dE_IMnpipi_woK0_wSid_n"),100, 1, 2, 200, 0, 50);
  dE_IMnpipi_woK0_wSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  dE_IMnpipi_woK0_wSid_n->SetYTitle("dE [MeVee]");
  dE_IMnpipi_woK0_wSid_n->GetXaxis()->CenterTitle();
  dE_IMnpipi_woK0_wSid_n->GetYaxis()->CenterTitle();
    
  Cosn_IMnpipi_woK0_wSid_n = new TH2F(Form("Cosn_IMnpipi_woK0_wSid_n"),Form("dE_Cosn_IMnpipi_woK0_wSid_n"),100, 1, 2, 50, -1, 1);
  Cosn_IMnpipi_woK0_wSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  Cosn_IMnpipi_woK0_wSid_n->SetYTitle("cos#theta_{n} (CM)");
  Cosn_IMnpipi_woK0_wSid_n->GetXaxis()->CenterTitle();
  Cosn_IMnpipi_woK0_wSid_n->GetYaxis()->CenterTitle();
    
  MMnmiss_IMnpipi_woK0_wSid_n = new TH2F(Form("MMnmiss_IMnpipi_woK0_wSid_n"),Form("MMnmiss_IMnpipi_woK0_wSid_n"),100,1,2,100,0,1.5);
  MMnmiss_IMnpipi_woK0_wSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMnpipi_woK0_wSid_n->SetYTitle("Miss Mass. [GeV/c^{2}]");
  MMnmiss_IMnpipi_woK0_wSid_n->GetXaxis()->CenterTitle();
  MMnmiss_IMnpipi_woK0_wSid_n->GetYaxis()->CenterTitle();
    
  q_IMnpipi_woK0_wSid_n = new TH2F(Form("q_IMnpipi_woK0_wSid_n"),Form("q_IMnpipi_woK0_wSid_n"),100,1,2,300,0,1.5);
  q_IMnpipi_woK0_wSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_wSid_n->SetYTitle("Mom. Transfer [GeV/c]");
  q_IMnpipi_woK0_wSid_n->GetXaxis()->CenterTitle();
  q_IMnpipi_woK0_wSid_n->GetYaxis()->CenterTitle();
    

  for(int imode=0;imode<2;imode++){
    dE_betainv_fid_kin[imode] = new TH2F(Form("dE_betainv_fid_kin_%s",smode[imode]),Form("dE_betainv_fid_kin_%s",smode[imode]),200, 0, 10, 200, 0, 50);
    dE_betainv_fid_kin[imode]->SetXTitle("1/#beta");
    dE_betainv_fid_kin[imode]->SetYTitle("dE [MeVee]");
    dE_betainv_fid_kin[imode]->GetXaxis()->CenterTitle();
    dE_betainv_fid_kin[imode]->GetYaxis()->CenterTitle();

    
    dE_MMom_fid_beta_woK0_kin[imode] = new TH2F(Form("dE_MMom_fid_beta_woK0_kin_%s",smode[imode]),Form("dE_MMom_fid_beta_woK0_kin_%s",smode[imode]),100, 0, 1.5, 200, 0, 50);
    dE_MMom_fid_beta_woK0_kin[imode]->SetXTitle("Missing Mom. [GeV/c]");
    dE_MMom_fid_beta_woK0_kin[imode]->SetYTitle("dE [MeVee]");
    dE_MMom_fid_beta_woK0_kin[imode]->GetXaxis()->CenterTitle();
    dE_MMom_fid_beta_woK0_kin[imode]->GetYaxis()->CenterTitle();


    dE_MMass_fid_beta_woK0_kin[imode] = new TH2F(Form("dE_MMass_fid_beta_woK0_kin_%s",smode[imode]),Form("dE_MMass_fid_beta_woK0_kin_%s",smode[imode]), 140, 0.4, 1.8, 200, 0, 50);
    dE_MMass_fid_beta_woK0_kin[imode]->SetXTitle("Missing mass [GeV/c^{2}]");
    dE_MMass_fid_beta_woK0_kin[imode]->SetYTitle("dE [MeVee]");
    dE_MMass_fid_beta_woK0_kin[imode]->GetXaxis()->CenterTitle();
    dE_MMass_fid_beta_woK0_kin[imode]->GetYaxis()->CenterTitle();
  

    MMom_MMass_fid_beta_dE_woK0_kin[imode] = new TH2F(Form("MMom_MMass_fid_beta_dE_woK0_kin_%s",smode[imode]),Form("MMom_MMass_fid_beta_dE_woK0_kin_%s",smode[imode]), 140, 0.4, 1.8, 100, 0, 1.5);
    MMom_MMass_fid_beta_dE_woK0_kin[imode]->SetXTitle("Missing Mass [GeV/c^{2}]");
    MMom_MMass_fid_beta_dE_woK0_kin[imode]->SetYTitle("Missing Mom. [GeV/c]");
    MMom_MMass_fid_beta_dE_woK0_kin[imode]->GetXaxis()->CenterTitle();
    MMom_MMass_fid_beta_dE_woK0_kin[imode]->GetYaxis()->CenterTitle();
    

    IMnpim_IMnpip_dE_woK0_kin[imode] = new TH2F(Form("IMnpim_IMnpip_dE_woK0_kin_%s",smode[imode]), Form("IMnpim_IMnpip_dE_woK0_kin_%s",smode[imode]),200, 1, 2.0, 200, 1, 2.0);
    IMnpim_IMnpip_dE_woK0_kin[imode]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
    IMnpim_IMnpip_dE_woK0_kin[imode]->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
    IMnpim_IMnpip_dE_woK0_kin[imode]->GetXaxis()->CenterTitle();
    IMnpim_IMnpip_dE_woK0_kin[imode]->GetXaxis()->CenterTitle();
  

    IMnpim_IMnpip_dE_woK0_n_kin[imode] = new TH2F(Form("IMnpim_IMnpip_dE_woK0_n_kin_%s",smode[imode]),Form("IMnpim_IMnpip_dE_woK0_n_kin_%s",smode[imode]),200, 1, 2.0, 200, 1, 2.0);
    IMnpim_IMnpip_dE_woK0_n_kin[imode]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
    IMnpim_IMnpip_dE_woK0_n_kin[imode]->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
    IMnpim_IMnpip_dE_woK0_n_kin[imode]->GetXaxis()->CenterTitle();
    IMnpim_IMnpip_dE_woK0_n_kin[imode]->GetYaxis()->CenterTitle();

  
    MMnpip_MMnpim_woK0_wSid_n_kin[imode] = new TH2F(Form("MMnpip_MMnpim_woK0_wSid_n_kin_%s",smode[imode]),Form("MMnpip_MMnpim_woK0_wSid_n_kin_%s",smode[imode]),70, 1, 1.7, 70, 1, 1.7);
    MMnpip_MMnpim_woK0_wSid_n_kin[imode]->SetXTitle("Miss. Mass(n#pi^{+}) [GeV/c^{2}]");
    MMnpip_MMnpim_woK0_wSid_n_kin[imode]->SetYTitle("Miss. Mass(n#pi^{-}) [GeV/c^{2}]");
    MMnpip_MMnpim_woK0_wSid_n_kin[imode]->GetXaxis()->CenterTitle();
    MMnpip_MMnpim_woK0_wSid_n_kin[imode]->GetYaxis()->CenterTitle();

 
    dE_IMnpipi_woK0_wSid_n_kin[imode] = new TH2F(Form("dE_IMnpipi_woK0_wSid_n_kin_%s",smode[imode]),Form("dE_IMnpipi_woK0_wSid_n_kin_%s",smode[imode]),100, 1, 2, 200, 0, 50);
    dE_IMnpipi_woK0_wSid_n_kin[imode]->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
    dE_IMnpipi_woK0_wSid_n_kin[imode]->SetYTitle("dE [MeVee]");
    dE_IMnpipi_woK0_wSid_n_kin[imode]->GetXaxis()->CenterTitle();
    dE_IMnpipi_woK0_wSid_n_kin[imode]->GetYaxis()->CenterTitle();

  
    Cosn_IMnpipi_woK0_wSid_n_kin[imode] = new TH2F(Form("Cosn_IMnpipi_woK0_wSid_n_kin_%s",smode[imode]),Form("dE_Cosn_IMnpipi_woK0_wSid_n_kin_%s",smode[imode]),100, 1, 2, 50, -1, 1);
    Cosn_IMnpipi_woK0_wSid_n_kin[imode]->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
    Cosn_IMnpipi_woK0_wSid_n_kin[imode]->SetYTitle("cos#theta_{n} (CM)");
    Cosn_IMnpipi_woK0_wSid_n_kin[imode]->GetXaxis()->CenterTitle();
    Cosn_IMnpipi_woK0_wSid_n_kin[imode]->GetYaxis()->CenterTitle();
  

    MMnmiss_IMnpipi_woK0_wSid_n_kin[imode] = new TH2F(Form("MMnmiss_IMnpipi_woK0_wSid_n_kin_%s",smode[imode]),Form("MMnmiss_IMnpipi_woK0_wSid_n_kin_%s",smode[imode]),100,1,2,100,0,1.5);
    MMnmiss_IMnpipi_woK0_wSid_n_kin[imode]->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
    MMnmiss_IMnpipi_woK0_wSid_n_kin[imode]->SetYTitle("Miss Mass. [GeV/c^{2}]");
    MMnmiss_IMnpipi_woK0_wSid_n_kin[imode]->GetXaxis()->CenterTitle();
    MMnmiss_IMnpipi_woK0_wSid_n_kin[imode]->GetYaxis()->CenterTitle();

    
    q_IMnpipi_woK0_wSid_n_kin[imode] = new TH2F(Form("q_IMnpipi_woK0_wSid_n_kin_%s",smode[imode]),Form("q_IMnpipi_woK0_wSid_n_kin_%s",smode[imode]),100,1,2,300,0,1.5);
    q_IMnpipi_woK0_wSid_n_kin[imode]->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
    q_IMnpipi_woK0_wSid_n_kin[imode]->SetYTitle("Mom. Transfer [GeV/c]");
    q_IMnpipi_woK0_wSid_n_kin[imode]->GetXaxis()->CenterTitle();
    q_IMnpipi_woK0_wSid_n_kin[imode]->GetYaxis()->CenterTitle();
  };

  TH2F *KFpvalue_vs = new TH2F("KFpvalue_vs", "KFpvalue_vs", 500, 0, 1, 500, 0, 1 );
  KFpvalue_vs->SetXTitle("P-value (#Sigma^{+} mode)");
  KFpvalue_vs->SetYTitle("P-value (#Sigma^{-} mode)");
  KFpvalue_vs->GetXaxis()->CenterTitle();
  KFpvalue_vs->GetYaxis()->CenterTitle();

  TH1F *nmom = new TH1F("nmom", "nmom", 50, 0, 1.0);
  nmom->SetXTitle("mom. [GeV/c]");
  nmom->GetXaxis()->CenterTitle();
  TH1F *nmom_kin = new TH1F("nmom_kin", "nmom_kin", 50, 0, 1.0);
  nmom_kin->SetXTitle("mom. [GeV/c]");
  nmom_kin->GetXaxis()->CenterTitle();
  TH1F *mnmom = new TH1F("mnmom", "mnmom", 100, 0, 2.0);
  mnmom->SetXTitle("mom. [GeV/c]");
  mnmom->GetXaxis()->CenterTitle();
  TH1F *mnmom_kin = new TH1F("mnmom_kin", "mnmom_kin", 100, 0, 2.0);
  mnmom_kin->SetXTitle("mom. [GeV/c]");
  mnmom_kin->GetXaxis()->CenterTitle();
  TH1F *npipmom = new TH1F("npipmom", "npipmom", 150, 0, 3.0);
  npipmom->SetXTitle("mom. [GeV/c]");
  npipmom->GetXaxis()->CenterTitle();
  TH1F *npipmom_kin = new TH1F("npipmom_kin", "npipmom_kin", 150, 0, 3.0);
  npipmom_kin->SetXTitle("mom. [GeV/c]");
  npipmom_kin->GetXaxis()->CenterTitle();
  TH1F *npimmom = new TH1F("npimmom", "npimmom", 150, 0, 3.0);
  npimmom->SetXTitle("mom. [GeV/c]");
  npimmom->GetXaxis()->CenterTitle();
  TH1F *npimmom_kin = new TH1F("npimmom_kin", "npimmom_kin", 150, 0, 3.0);
  npimmom_kin->SetXTitle("mom. [GeV/c]");
  npimmom_kin->GetXaxis()->CenterTitle();

  Int_t nevent = tree->GetEntries();
  std::cerr<<"# of events = "<<nevent<<std::endl;
  //------------------------//
  //--- event roop start ---//
  //------------------------//
  for ( Int_t i=0; i<nevent; i++ ) {
    tree->GetEvent(i);
    //if(i%1000==0) std::cout << "Event# " << i << std::endl; 
    
    // calc missing n //
    TLorentzVector LVec_n_miss = *LVec_target+*LVec_beam-*LVec_pip-*LVec_pim-*LVec_n;
    double nmiss_mass = LVec_n_miss.M();
    double nmiss_mom = LVec_n_miss.P();

    // calc cos(theta) of missing n //
    TVector3 boost = (*LVec_target+*LVec_beam).BoostVector();
    TLorentzVector LVec_n_miss_CM = LVec_n_miss;
    TLorentzVector LVec_beam_CM = *LVec_beam;
    LVec_n_miss_CM.Boost(-boost);
    LVec_beam_CM.Boost(-boost);
    double cos_n = LVec_n_miss_CM.Vect().Dot(LVec_beam_CM.Vect())/(LVec_n_miss_CM.Vect().Mag()*LVec_beam_CM.Vect().Mag());
    TLorentzVector qkn = *LVec_beam-LVec_n_miss;
   
    // calc pi+pi- //
    TLorentzVector LVec_pip_pim = *LVec_pip+*LVec_pim;

    // calc pi+n //
    TLorentzVector LVec_pip_n = *LVec_pip+*LVec_n;
    // calc pi-n //
    TLorentzVector LVec_pim_n = *LVec_pim+*LVec_n;

    // calc missing Sp //
    TLorentzVector LVec_pip_n_miss = *LVec_target+*LVec_beam-*LVec_pim-*LVec_n;
    npipmom->Fill(LVec_pip_n.P());
    
    // calc missing Sm //
    TLorentzVector LVec_pim_n_miss = *LVec_target+*LVec_beam-*LVec_pip-*LVec_n;
    npimmom->Fill(LVec_pim_n.P());

    // calc pi+pi-n //
    TLorentzVector LVec_pip_pim_n = *LVec_pip+*LVec_pim+*LVec_n;
    TLorentzVector LVec_pip_pim_n_CM = LVec_pip_pim_n;
    LVec_pip_pim_n_CM.Boost(-boost);
    double cos_X = LVec_pip_pim_n_CM.Vect().Dot(LVec_beam_CM.Vect())/(LVec_pip_pim_n_CM.Vect().Mag()*LVec_beam_CM.Vect().Mag());


    double chi2 = kfSpmode_chi2<kfSmmode_chi2 ? kfSpmode_chi2:kfSmmode_chi2;
    double pvalue = kfSmmode_pvalue<kfSpmode_pvalue ? kfSpmode_pvalue:kfSmmode_pvalue;
    

    bool K0rejectFlag=false;
    bool MissNFlag=false;
    bool NBetaOK=false;
    bool NdEOK=false;
    bool SigmaPFlag=false;
    bool SigmaMFlag=false;
    //-- neutron-ID, K0 and missing neutron selection --//

    if(NeutralBetaCDH<anacuts::beta_MAX) NBetaOK=true;
    if(anacuts::dE_MIN<dE) NdEOK=true;
    
    //Sigma+ production in CDS
    if( (anacuts::Sigmap_MIN<(*LVec_n+*LVec_pip).M() && (*LVec_n+*LVec_pip).M()<anacuts::Sigmap_MAX)) SigmaPFlag=true;
        
    //Sigma- production in CDS
    if( (anacuts::Sigmam_MIN<(*LVec_n+*LVec_pim).M() && (*LVec_n+*LVec_pim).M()<anacuts::Sigmam_MAX)) SigmaMFlag=true;
     
    if(anacuts::neutron_MIN<nmiss_mass && nmiss_mass<anacuts::neutron_MAX ) MissNFlag=true;
    
    
    
    //K0 rejection using original momentum
    if( (LVec_pip_pim.M()<anacuts::pipi_MIN || anacuts::pipi_MAX<LVec_pip_pim.M())) K0rejectFlag=true;
     
    if(K0rejectFlag && NBetaOK && NdEOK) KFpvalue_vs->Fill(kfSpmode_pvalue,kfSmmode_pvalue);
    //w/o kinfit
    if(K0rejectFlag ){
      dE_betainv_fid->Fill(1./NeutralBetaCDH,dE);
      if(NBetaOK){
        dE_MMom_fid_beta_woK0->Fill(LVec_n_miss.P(),dE);
        dE_MMass_fid_beta_woK0->Fill(LVec_n_miss.M(),dE);
      }
      if(NBetaOK && NdEOK){
        MMom_MMass_fid_beta_dE_woK0->Fill(LVec_n_miss.M(),LVec_n_miss.P());
        IMnpim_IMnpip_dE_woK0->Fill(LVec_pip_n.M(),LVec_pim_n.M());
        nmom->Fill((*LVec_n).P());
        mnmom->Fill(nmiss_mom);
      }
      if(NBetaOK && NdEOK && MissNFlag){
        IMnpim_IMnpip_dE_woK0_n->Fill(LVec_pip_n.M(),LVec_pim_n.M());
        if(SigmaPFlag || SigmaMFlag){
          MMnpip_MMnpim_woK0_wSid_n->Fill(LVec_pim_n_miss.M(),LVec_pip_n_miss.M());
          dE_IMnpipi_woK0_wSid_n->Fill(LVec_pip_pim_n.M(),dE);
          Cosn_IMnpipi_woK0_wSid_n->Fill(LVec_pip_pim_n.M(),cos_n);
          MMnmiss_IMnpipi_woK0_wSid_n->Fill(LVec_pip_pim_n.M(), nmiss_mass);
          q_IMnpipi_woK0_wSid_n->Fill(LVec_pip_pim_n.M(),qkn.P());
        }
      }
    }
    
    // w/ kinfit
    if( -1<kf_flag && 0.01<pvalue && K0rejectFlag ){
      nmom_kin->Fill((*LVec_n).P());
      mnmom_kin->Fill(nmiss_mom);
      npipmom_kin->Fill(LVec_pip_n.P());
      npimmom_kin->Fill(LVec_pim_n.P());

      if( kfSmmode_pvalue<kfSpmode_pvalue ){
        dE_betainv_fid_kin[0]->Fill(1./NeutralBetaCDH,dE);
        if(NBetaOK){
          dE_MMom_fid_beta_woK0_kin[0]->Fill(LVec_n_miss.P(),dE);
          dE_MMass_fid_beta_woK0_kin[0]->Fill(LVec_n_miss.M(),dE);
        }
        if(NBetaOK && NdEOK){
          MMom_MMass_fid_beta_dE_woK0_kin[0]->Fill(LVec_n_miss.M(),LVec_n_miss.P());
          IMnpim_IMnpip_dE_woK0_kin[0]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
        }
        if(NBetaOK && NdEOK && MissNFlag){
          IMnpim_IMnpip_dE_woK0_n_kin[0]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
          MMnpip_MMnpim_woK0_wSid_n_kin[0]->Fill(LVec_pim_n_miss.M(),LVec_pip_n_miss.M());
          dE_IMnpipi_woK0_wSid_n_kin[0]->Fill(LVec_pip_pim_n.M(),dE);
          Cosn_IMnpipi_woK0_wSid_n_kin[0]->Fill(LVec_pip_pim_n.M(),cos_n);
          MMnmiss_IMnpipi_woK0_wSid_n_kin[0]->Fill(LVec_pip_pim_n.M(), nmiss_mass);
          q_IMnpipi_woK0_wSid_n_kin[0]->Fill(LVec_pip_pim_n.M(),qkn.P());
        }
      
      
      }else{//S- mode
        dE_betainv_fid_kin[1]->Fill(1./NeutralBetaCDH,dE);
        if(NBetaOK){
          dE_MMom_fid_beta_woK0_kin[1]->Fill(LVec_n_miss.P(),dE);
          dE_MMass_fid_beta_woK0_kin[1]->Fill(LVec_n_miss.M(),dE);
        }
        if(NBetaOK && NdEOK){
          MMom_MMass_fid_beta_dE_woK0_kin[1]->Fill(LVec_n_miss.M(),LVec_n_miss.P());
          IMnpim_IMnpip_dE_woK0_kin[1]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
        }
        if(NBetaOK && NdEOK && MissNFlag){
          IMnpim_IMnpip_dE_woK0_n_kin[1]->Fill(LVec_pip_n.M(),LVec_pim_n.M());
          MMnpip_MMnpim_woK0_wSid_n_kin[1]->Fill(LVec_pim_n_miss.M(),LVec_pip_n_miss.M());
          dE_IMnpipi_woK0_wSid_n_kin[1]->Fill(LVec_pip_pim_n.M(),dE);
          Cosn_IMnpipi_woK0_wSid_n_kin[1]->Fill(LVec_pip_pim_n.M(),cos_n);
          MMnmiss_IMnpipi_woK0_wSid_n_kin[1]->Fill(LVec_pip_pim_n.M(), nmiss_mass);
          q_IMnpipi_woK0_wSid_n_kin[1]->Fill(LVec_pip_pim_n.M(),qkn.P());
        }
      }
    }
    
	}//for ievt
   
  TCanvas *cKFpvalue_vs = new TCanvas(Form("cKFpvalue_vs"),"KFpvalue_vs");
  cKFpvalue_vs->cd();
  TH1D *px = (TH1D*) KFpvalue_vs->ProjectionX();
  px->SetMinimum(1);
  px->SetXTitle("p-value");
  px->Draw();
  TH1D *py = (TH1D*) KFpvalue_vs->ProjectionY();
  py->SetLineColor(2);
  py->Draw("same");
  TLegend *legKFpvalue_vs = new TLegend(0.55,0.65,0.76,0.82);
  legKFpvalue_vs->AddEntry(px,"#Sigma^{+} mode");
  legKFpvalue_vs->AddEntry(py,"#Sigma^{-} mode");
  legKFpvalue_vs->SetFillColor(0);
  legKFpvalue_vs->Draw();
  gPad->SetLogy();
  int spbin = px->FindBin(pvalcut);
  //std::cout << spbin << std::endl;
  int smbin = py->FindBin(pvalcut);
  std::cout << "Sp mode rejection ratio:" << px->Integral(0,spbin)/(px->Integral(0,201)) << std::endl;
  std::cout << "Sm mode rejection ratio:" << py->Integral(0,smbin)/(py->Integral(0,201)) << std::endl;

  //cumulative dist. of prob.
  TCanvas *cKFpvalue_vs_cum = new TCanvas(Form("cKFpvalue_vs_cum"),"KFpvalue_vs_cum");
  cKFpvalue_vs_cum->cd();
  TH1 *px_cum = px->GetCumulative();
  px_cum->Scale(1./px->GetEntries());
  px_cum->SetXTitle("p-value cut");
  px_cum->GetXaxis()->CenterTitle();
  px_cum->SetTitle("KF pvalue cumulative");
  px_cum->Draw();
  TH1 *py_cum = py->GetCumulative();
  py_cum->SetLineColor(2);
  py_cum->Scale(1./py->GetEntries());
  py_cum->Draw("same");
  TLegend *legKFpvalue_vs_cum = new TLegend(0.55,0.25,0.76,0.42);
  legKFpvalue_vs_cum->AddEntry(px,"#Sigma^{+} mode");
  legKFpvalue_vs_cum->AddEntry(py,"#Sigma^{-} mode");
  legKFpvalue_vs_cum->SetFillColor(0);
  legKFpvalue_vs_cum->Draw();
  
  
  TCanvas *cnmom = new TCanvas("cnmom","cnmom");
  cnmom->cd();
  nmom->Draw("");
  gPad->SetLogy(0);

  TCanvas *cmnmom = new TCanvas("cmnmom","cmnmom");
  cmnmom->cd();
  mnmom->Draw("");
  
  TCanvas *cnpipmom = new TCanvas("cnpipmom","cnpipmom");
  cnpipmom->cd();
  npipmom->Draw("");
  
  TCanvas *cnpimmom = new TCanvas("cnpimmom","cnpimmom");
  cnpimmom->cd();
  npimmom->Draw("");
  
  TCanvas *c_IMnpipi_woK0_wSid_n = new TCanvas("c_IMnpipi_woK0_wSid_n","c_IMnpipi_woK0_wSid_n"); 
  c_IMnpipi_woK0_wSid_n->cd();
  TH1D *pxSp = q_IMnpipi_woK0_wSid_n_kin[0]->ProjectionX();
  TH1D *pxSm = q_IMnpipi_woK0_wSid_n_kin[1]->ProjectionX();
  pxSm->SetLineColor(2);
  pxSm->Rebin(2);
  pxSm->Draw("HE");
  pxSp->Rebin(2);
  pxSp->Draw("HEsame");

  TCanvas *cdE_betainv_fid = new TCanvas("cdE_betainv_fid","cdE_betainv_fid");
  cdE_betainv_fid->cd();
  dE_betainv_fid->Draw("colz");

  TCanvas *cdE_betainv_fid_kin[2];//           
  TCanvas *cdE_MMom_fid_beta_woK0[2];      
  TCanvas *cdE_MMass_fid_beta_woK0[2];     
  TCanvas *cMMom_MMass_fid_beta_dE_woK0[2];
  TCanvas *cIMnpim_IMnpip_dE_woK0[2];      
  TCanvas *cIMnpim_IMnpip_dE_woK0_n[2];    
  TCanvas *cMMnpip_MMnpim_woK0_wSid_n[2];  
  TCanvas *cdE_IMnpipi_woK0_wSid_n[2];     
  TCanvas *cCosn_IMnpipi_woK0_wSid_n[2];   
  TCanvas *cMMnmiss_IMnpipi_woK0_wSid_n[2];
  TCanvas *cq_IMnpipi_woK0_wSid_n[2];      

  for(int imode=0;imode<2;imode++){
    cdE_betainv_fid_kin[imode] = new TCanvas(Form("cdE_betainv_fid_kin_%s",smode[imode]), "");
    cdE_betainv_fid_kin[imode]->cd();
    dE_betainv_fid_kin[imode]->Draw("colz");

    cdE_MMom_fid_beta_woK0[imode] = new TCanvas(Form("cdE_MMom_fid_beta_woK0_%s",smode[imode]),"");
    cdE_MMom_fid_beta_woK0[imode]->cd();
    dE_MMom_fid_beta_woK0_kin[imode]->Draw("colz");

    cdE_MMass_fid_beta_woK0[imode] = new TCanvas(Form("cdE_MMass_fid_beta_woK0_%s",smode[imode]),"");
    cdE_MMass_fid_beta_woK0[imode]->cd(); 
    dE_MMass_fid_beta_woK0_kin[imode]->Draw("colz");

    cMMom_MMass_fid_beta_dE_woK0[imode] = new TCanvas(Form("cMMom_MMass_fid_beta_dE_woK0_%s",smode[imode]),"");
    cMMom_MMass_fid_beta_dE_woK0[imode]->cd();
    MMom_MMass_fid_beta_dE_woK0_kin[imode]->Draw("colz");
    
    //cIMnpim_IMnpip_dE_woK0[imode] = new TCanvas(Form("cIMnpim_IMnpip_dE_woK0_%s",smode[imode]),"");
    //cIMnpim_IMnpip_dE_woK0[imode]->cd();
    //IMnpim_IMnpip_dE_woK0[imode]->Draw("colz");
    
    cIMnpim_IMnpip_dE_woK0_n[imode] = new TCanvas(Form("cIMnpim_IMnpip_dE_woK0_n_%s",smode[imode]),"");
    cIMnpim_IMnpip_dE_woK0_n[imode]->cd();
    IMnpim_IMnpip_dE_woK0_n_kin[imode]->Draw("colz");

    cMMnpip_MMnpim_woK0_wSid_n[imode] = new TCanvas(Form("cMMnpip_MMnpim_woK0_wSid_n_%s",smode[imode]),"");
    cMMnpip_MMnpim_woK0_wSid_n[imode]->cd();
    MMnpip_MMnpim_woK0_wSid_n_kin[imode]->Draw("colz");

    cdE_IMnpipi_woK0_wSid_n[imode] = new TCanvas(Form("cdE_IMnpipi_woK0_wSid_n_%s",smode[imode]),"");
    cdE_IMnpipi_woK0_wSid_n[imode]->cd();
    dE_IMnpipi_woK0_wSid_n_kin[imode]->Draw("colz");
    
    cCosn_IMnpipi_woK0_wSid_n[imode] = new TCanvas(Form("cCosn_IMnpipi_woK0_wSid_n_%s",smode[imode]),"");
    cCosn_IMnpipi_woK0_wSid_n[imode]->cd();
    Cosn_IMnpipi_woK0_wSid_n_kin[imode]->Draw("colz");

    cMMnmiss_IMnpipi_woK0_wSid_n[imode] = new TCanvas(Form("cMMnmiss_IMnpipi_woK0_wSid_n_%s",smode[imode]),"");
    cMMnmiss_IMnpipi_woK0_wSid_n[imode]->cd();
    MMnmiss_IMnpipi_woK0_wSid_n_kin[imode]->Draw("colz");

    cq_IMnpipi_woK0_wSid_n[imode] = new TCanvas(Form("cq_IMnpipi_woK0_wSid_n_%s",smode[imode]),"");
    cq_IMnpipi_woK0_wSid_n[imode]->cd();
    q_IMnpipi_woK0_wSid_n_kin[imode]->Draw("colz");
  }



  TFile *fout = new TFile(outfilename.c_str(),"RECREATE");
  fout->cd();
  for(int imode=0;imode<2;imode++){
    dE_betainv_fid_kin[imode]->Write();
    dE_MMom_fid_beta_woK0_kin[imode]->Write();
    dE_MMass_fid_beta_woK0_kin[imode]->Write();
    MMom_MMass_fid_beta_dE_woK0_kin[imode]->Write();
    IMnpim_IMnpip_dE_woK0_kin[imode]->Write();
    IMnpim_IMnpip_dE_woK0_n_kin[imode]->Write();
    MMnpip_MMnpim_woK0_wSid_n_kin[imode]->Write();
    dE_IMnpipi_woK0_wSid_n_kin[imode]->Write();
    Cosn_IMnpipi_woK0_wSid_n_kin[imode]->Write();
    MMnmiss_IMnpipi_woK0_wSid_n_kin[imode]->Write();
    q_IMnpipi_woK0_wSid_n_kin[imode]->Write();
    nmom->Write();
    mnmom->Write();
    npipmom->Write();
    npimmom->Write();
    nmom_kin->Write();
    mnmom_kin->Write();
    npipmom_kin->Write();
    npimmom_kin->Write();
  };

  
  fout->Close();
  
}
