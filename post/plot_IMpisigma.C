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
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TLegend.h>
#include <TPDF.h>
#include <TPaveText.h>

#include "../src/IMPiSigmaAnaPar.h"

const double pvalcut = 0.005;
const double dEcut = 2.0;
//const double pvalcut = 1.0e-5;
const bool gridon=true;
const bool staton=false;

//mode 0: Sigma+ ,1: Sigma- 
void plot_IMpisigma(const char* filename="",const int mode=0)
{
  std::cout << "p-value cut:" << pvalcut << std::endl; 
  std::cout << "dE cut:" << dEcut << std::endl; 
  TCanvas *ctext = new TCanvas("ctext","ctext");
  TPaveText *pt = new TPaveText(.05,.05,.95,.7);
  pt->AddText(Form("p-value cut: %f ",pvalcut));
  pt->AddText(Form("dE cut: %f " ,dEcut));
  pt->AddText(Form("1/beta min.: %f ",1./anacuts::beta_MAX));
  pt->Draw(); 

  //gROOT->SetStyle("Plain");
  if(staton)gStyle->SetOptStat(111111);
  else gStyle->SetOptStat(0);
  gStyle->SetOptFit(111111);
  gStyle->SetPadGridX(gridon);
  gStyle->SetPadGridY(gridon);
  std::cout << "infile " << filename <<std::endl;
  std::string outfilename = std::string(filename);
  outfilename.insert(outfilename.size()-5,"_post");
  std::cout << "outfilename: " << outfilename << std::endl;
  TString pdfname = outfilename;
  pdfname.Replace(outfilename.size()-4,5,"pdf");
  std::cout << "pdfname: " << pdfname << std::endl;

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
  TVector3 *vtx_pip = nullptr; // vertex (pip)
  TVector3 *vtx_pim = nullptr; // vertex (pim)
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
  tree->SetBranchAddress( "mom_target", &LVec_target );
  tree->SetBranchAddress( "mom_pip", &LVec_pip );
  tree->SetBranchAddress( "mom_pim", &LVec_pim );
  tree->SetBranchAddress( "mom_n", &LVec_n );
  tree->SetBranchAddress( "NeutralBetaCDH", &NeutralBetaCDH );
  tree->SetBranchAddress( "dE", &dE );
  tree->SetBranchAddress( "vtx_reaction", &vtx_reaction );
  tree->SetBranchAddress( "vtx_pip",&vtx_pip);
  tree->SetBranchAddress( "vtx_pim",&vtx_pim);
  //tree->SetBranchAddress( "run_num", &run_num );
  //tree->SetBranchAddress( "event_num", &event_num );
  //tree->SetBranchAddress( "block_num", &block_num );
  //tree->SetBranchAddress( "kfSpmode_mom_beam",   &kfSpmode_mom_beam );
  //tree->SetBranchAddress( "kfSpmode_mom_pip", &kfSpmode_mom_pip );
  //tree->SetBranchAddress( "kfSpmode_mom_pim", &kfSpmode_mom_pim );
  //tree->SetBranchAddress( "kfSpmode_mom_n", &kfSpmode_mom_n );
  tree->SetBranchAddress( "kfSpmode_chi2", &kfSpmode_chi2 );
  tree->SetBranchAddress( "kfSpmode_NDF", &kfSpmode_NDF );
  tree->SetBranchAddress( "kfSpmode_status", &kfSpmode_status );
  tree->SetBranchAddress( "kfSpmode_pvalue", &kfSpmode_pvalue );
  //tree->SetBranchAddress( "kfSmmode_mom_beam",   &kfSmmode_mom_beam );
  //tree->SetBranchAddress( "kfSmmode_mom_pip", &kfSmmode_mom_pip );
  //tree->SetBranchAddress( "kfSmmode_mom_pim", &kfSmmode_mom_pim );
  //tree->SetBranchAddress( "kfSmmode_mom_n", &kfSmmode_mom_n );
  tree->SetBranchAddress( "kfSmmode_chi2", &kfSmmode_chi2 );
  tree->SetBranchAddress( "kfSmmode_NDF", &kfSmmode_NDF );
  tree->SetBranchAddress( "kfSmmode_status", &kfSmmode_status );
  tree->SetBranchAddress( "kfSmmode_pvalue", &kfSmmode_pvalue );
  tree->SetBranchAddress( "kf_flag", &kf_flag );
  
  
  // w/o kinematic fit 
  TH2F* dE_betainv_fid;//
  TH2F* dE_MMom_fid_beta_woK0;
  TH2F* dE_MMass_fid_beta_woK0;
  TH2F* MMom_MMass_fid_beta_dE_woK0;
  TH2F* MMom_MMass_fid_beta_dE_woK0_wSid;
  TH2F* IMnpim_IMnpip_dE_woK0;
  TH2F* IMnpim_IMnpip_dE_woK0_n;
  TH2F* MMnpip_MMnpim_woK0_wSid_n;
  TH2F* dE_IMnpim_woK0;
  TH2F* dE_IMnpim_woK0_n;
  TH2F* dE_IMnpip_woK0;
  TH2F* dE_IMnpip_woK0_n;
  TH2F* dE_IMnpipi_woK0_wSid_n;
  TH2F* Cosn_IMnpipi_woK0_wSid_n;
  TH2F* MMnmiss_IMnpipi_woK0_wSid_n;
  TH2F* q_IMnpipi_woK0_wSid_n;
  TH2F* nmom_IMnpipi_woK0_wSid_n;
  // w/ kinematic fit
  TH2F* dE_betainv_fid_kin[2];//
  TH2F* dE_MMom_fid_beta_woK0_kin[2];
  TH2F* dE_MMass_fid_beta_woK0_kin[2];
  TH2F* MMom_MMass_fid_beta_dE_woK0_kin[2];
  TH2F* IMnpim_IMnpip_dE_woK0_kin[2];
  TH2F* MMnpip_MMnpim_woK0_kin[2];
  TH2F* dE_IMnpipi_woK0_kin[2];
  TH2F* Cosn_IMnpipi_woK0_kin[2];
  TH2F* MMnmiss_IMnpipi_woK0_kin[2];
  TH2F* q_IMnpipi_woK0_kin[2];
  TH2F* nmom_IMnpipi_woK0_kin[2];
  const char smode[][4]={"Sp","Sm"};
  
  
  dE_betainv_fid = new TH2F(Form("dE_betainv_fid"),Form("dE_betainv_fid"),1000, 0, 50, 200, 0, 50);
  dE_betainv_fid->SetXTitle("1/#beta");
  dE_betainv_fid->SetYTitle("dE [MeVee]");
  
  dE_MMom_fid_beta_woK0 = new TH2F(Form("dE_MMom_fid_beta_woK0"),Form("dE_MMom_fid_beta_woK0"),100, 0, 1.5, 200, 0, 50);
  dE_MMom_fid_beta_woK0->SetXTitle("Missing Mom. [GeV/c]");
  dE_MMom_fid_beta_woK0->SetYTitle("dE [MeVee]");

  dE_MMass_fid_beta_woK0 = new TH2F(Form("dE_MMass_fid_beta_woK0"),Form("dE_MMass_fid_beta_woK0"), 140, 0.4, 1.8, 200, 0, 50);
  dE_MMass_fid_beta_woK0->SetXTitle("Missing mass [GeV/c^{2}]");
  dE_MMass_fid_beta_woK0->SetYTitle("dE [MeVee]");

  MMom_MMass_fid_beta_dE_woK0 = new TH2F(Form("MMom_MMass_fid_beta_dE_woK0"),Form("MMom_MMass_fid_beta_dE_woK0"), 140, 0.4, 1.8, 100, 0, 1.5);
  MMom_MMass_fid_beta_dE_woK0->SetXTitle("Missing Mass [GeV/c^{2}]");
  MMom_MMass_fid_beta_dE_woK0->SetYTitle("Missing Mom. [GeV/c]");


  MMom_MMass_fid_beta_dE_woK0_wSid = new TH2F(Form("MMom_MMass_fid_beta_dE_woK0_wSid"),Form("MMom_MMass_fid_beta_dE_woK0_wSid"), 140, 0.4, 1.8, 100, 0, 1.5);
  MMom_MMass_fid_beta_dE_woK0_wSid->SetXTitle("Missing Mass [GeV/c^{2}]");
  MMom_MMass_fid_beta_dE_woK0_wSid->SetYTitle("Missing Mom. [GeV/c]");

  
  IMnpim_IMnpip_dE_woK0 = new TH2F(Form("IMnpim_IMnpip_dE_woK0"), Form("IMnpim_IMnpip_dE_woK0"),200, 1, 2.0, 200, 1, 2.0);
  IMnpim_IMnpip_dE_woK0->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_woK0->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
    
  IMnpim_IMnpip_dE_woK0_n = new TH2F(Form("IMnpim_IMnpip_dE_woK0_n"),Form("IMnpim_IMnpip_dE_woK0_n"),200, 1, 2.0, 200, 1, 2.0);
  IMnpim_IMnpip_dE_woK0_n->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  IMnpim_IMnpip_dE_woK0_n->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
    
  MMnpip_MMnpim_woK0_wSid_n = new TH2F(Form("MMnpip_MMnpim_woK0_wSid_n"),Form("MMnpip_MMnpim_woK0_wSid_n"),70, 1, 1.7, 70, 1, 1.7);
  MMnpip_MMnpim_woK0_wSid_n->SetXTitle("Miss. Mass(n#pi^{+}) [GeV/c^{2}]");
  MMnpip_MMnpim_woK0_wSid_n->SetYTitle("Miss. Mass(n#pi^{-}) [GeV/c^{2}]");
  
  dE_IMnpim_woK0 = new TH2F(Form("dE_IMnpim_woK0"),Form("dE_IMnpim_woK0"), 200, 1.0, 2.0, 200, 0, 50.);
  dE_IMnpim_woK0->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  dE_IMnpim_woK0->SetYTitle("dE [MeVee]");
  
  dE_IMnpim_woK0_n = new TH2F(Form("dE_IMnpim_woK0_n"),Form("dE_IMnpim_woK0_n"), 200, 1.0, 2.0, 200, 0, 50.);
  dE_IMnpim_woK0_n->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  dE_IMnpim_woK0_n->SetYTitle("dE [MeVee]");

  dE_IMnpip_woK0 = new TH2F(Form("dE_IMnpip_woK0"),Form("dE_IMnpip_woK0"), 200, 1.0, 2.0, 200, 0, 50.);
  dE_IMnpip_woK0->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  dE_IMnpip_woK0->SetYTitle("dE [MeVee]");
  
  dE_IMnpip_woK0_n = new TH2F(Form("dE_IMnpip_woK0_n"),Form("dE_IMnpip_woK0_n"), 200, 1.0, 2.0, 200, 0, 50.);
  dE_IMnpip_woK0_n->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  dE_IMnpip_woK0_n->SetYTitle("dE [MeVee]");

  dE_IMnpipi_woK0_wSid_n = new TH2F(Form("dE_IMnpipi_woK0_wSid_n"),Form("dE_IMnpipi_woK0_wSid_n"),100, 1, 2, 200, 0, 50);
  dE_IMnpipi_woK0_wSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  dE_IMnpipi_woK0_wSid_n->SetYTitle("dE [MeVee]");
    
  Cosn_IMnpipi_woK0_wSid_n = new TH2F(Form("Cosn_IMnpipi_woK0_wSid_n"),Form("dE_Cosn_IMnpipi_woK0_wSid_n"),100, 1, 2, 50, -1, 1);
  Cosn_IMnpipi_woK0_wSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  Cosn_IMnpipi_woK0_wSid_n->SetYTitle("cos#theta_{n} (CM)");
    
  MMnmiss_IMnpipi_woK0_wSid_n = new TH2F(Form("MMnmiss_IMnpipi_woK0_wSid_n"),Form("MMnmiss_IMnpipi_woK0_wSid_n"),100,1,2,100,0,1.5);
  MMnmiss_IMnpipi_woK0_wSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMnpipi_woK0_wSid_n->SetYTitle("Miss Mass. [GeV/c^{2}]");
    
  q_IMnpipi_woK0_wSid_n = new TH2F(Form("q_IMnpipi_woK0_wSid_n"),Form("q_IMnpipi_woK0_wSid_n"),100,1,2,300,0,1.5);
  q_IMnpipi_woK0_wSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_woK0_wSid_n->SetYTitle("Mom. Transfer [GeV/c]");
  
  nmom_IMnpipi_woK0_wSid_n = new TH2F(Form("nmom_IMnpipi_woK0_wSid_n"),Form("nmom_IMnpipi_woK0_wSid_n"),100,1,2,100,0,1.0);
  nmom_IMnpipi_woK0_wSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  nmom_IMnpipi_woK0_wSid_n->SetYTitle("nmom  [GeV/c]");
    

  for(int imode=0;imode<2;imode++){
    dE_betainv_fid_kin[imode] = new TH2F(Form("dE_betainv_fid_kin_%s",smode[imode]),Form("dE_betainv_fid_kin_%s",smode[imode]),1000, 0, 50, 200, 0, 50);
    dE_betainv_fid_kin[imode]->SetXTitle("1/#beta");
    dE_betainv_fid_kin[imode]->SetYTitle("dE [MeVee]");

    
    dE_MMom_fid_beta_woK0_kin[imode] = new TH2F(Form("dE_MMom_fid_beta_woK0_kin_%s",smode[imode]),Form("dE_MMom_fid_beta_woK0_kin_%s",smode[imode]),100, 0, 1.5, 200, 0, 50);
    dE_MMom_fid_beta_woK0_kin[imode]->SetXTitle("Missing Mom. [GeV/c]");
    dE_MMom_fid_beta_woK0_kin[imode]->SetYTitle("dE [MeVee]");


    dE_MMass_fid_beta_woK0_kin[imode] = new TH2F(Form("dE_MMass_fid_beta_woK0_kin_%s",smode[imode]),Form("dE_MMass_fid_beta_woK0_kin_%s",smode[imode]), 140, 0.4, 1.8, 200, 0, 50);
    dE_MMass_fid_beta_woK0_kin[imode]->SetXTitle("Missing mass [GeV/c^{2}]");
    dE_MMass_fid_beta_woK0_kin[imode]->SetYTitle("dE [MeVee]");
  

    MMom_MMass_fid_beta_dE_woK0_kin[imode] = new TH2F(Form("MMom_MMass_fid_beta_dE_woK0_kin_%s",smode[imode]),Form("MMom_MMass_fid_beta_dE_woK0_kin_%s",smode[imode]), 140, 0.4, 1.8, 100, 0, 1.5);
    MMom_MMass_fid_beta_dE_woK0_kin[imode]->SetXTitle("Missing Mass [GeV/c^{2}]");
    MMom_MMass_fid_beta_dE_woK0_kin[imode]->SetYTitle("Missing Mom. [GeV/c]");
    

    IMnpim_IMnpip_dE_woK0_kin[imode] = new TH2F(Form("IMnpim_IMnpip_dE_woK0_kin_%s",smode[imode]), Form("IMnpim_IMnpip_dE_woK0_kin_%s",smode[imode]),200, 1, 2.0, 200, 1, 2.0);
    IMnpim_IMnpip_dE_woK0_kin[imode]->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
    IMnpim_IMnpip_dE_woK0_kin[imode]->SetYTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  

    MMnpip_MMnpim_woK0_kin[imode] = new TH2F(Form("MMnpip_MMnpim_woK0_kin_%s",smode[imode]),Form("MMnpip_MMnpim_woK0_kin_%s",smode[imode]),70, 1, 1.7, 70, 1, 1.7);
    MMnpip_MMnpim_woK0_kin[imode]->SetXTitle("Miss. Mass(n#pi^{+}) [GeV/c^{2}]");
    MMnpip_MMnpim_woK0_kin[imode]->SetYTitle("Miss. Mass(n#pi^{-}) [GeV/c^{2}]");

 
    dE_IMnpipi_woK0_kin[imode] = new TH2F(Form("dE_IMnpipi_woK0_kin_%s",smode[imode]),Form("dE_IMnpipi_woK0_kin_%s",smode[imode]),100, 1, 2, 200, 0, 50);
    dE_IMnpipi_woK0_kin[imode]->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
    dE_IMnpipi_woK0_kin[imode]->SetYTitle("dE [MeVee]");

  
    Cosn_IMnpipi_woK0_kin[imode] = new TH2F(Form("Cosn_IMnpipi_woK0_kin_%s",smode[imode]),Form("dE_Cosn_IMnpipi_woK0_kin_%s",smode[imode]),100, 1, 2, 50, -1, 1);
    Cosn_IMnpipi_woK0_kin[imode]->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
    Cosn_IMnpipi_woK0_kin[imode]->SetYTitle("cos#theta_{n} (CM)");
  

    MMnmiss_IMnpipi_woK0_kin[imode] = new TH2F(Form("MMnmiss_IMnpipi_woK0_kin_%s",smode[imode]),Form("MMnmiss_IMnpipi_woK0_kin_%s",smode[imode]),100,1,2,100,0,1.5);
    MMnmiss_IMnpipi_woK0_kin[imode]->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
    MMnmiss_IMnpipi_woK0_kin[imode]->SetYTitle("Miss Mass. [GeV/c^{2}]");

    
    q_IMnpipi_woK0_kin[imode] = new TH2F(Form("q_IMnpipi_woK0_kin_%s",smode[imode]),Form("q_IMnpipi_woK0_kin_%s",smode[imode]),100,1,2,300,0,1.5);
    q_IMnpipi_woK0_kin[imode]->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
    q_IMnpipi_woK0_kin[imode]->SetYTitle("Mom. Transfer [GeV/c]");


    nmom_IMnpipi_woK0_kin[imode] = new TH2F(Form("nmom_IMnpipi_woK0_kin_%s",smode[imode]),Form("nmom_IMnpipi_woK0_kin_%s",smode[imode]),100,1,2,100,0,1.0);
    nmom_IMnpipi_woK0_kin[imode]->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
    nmom_IMnpipi_woK0_kin[imode]->SetYTitle("nmom  [GeV/c]");
  };

  TH2F *KFpvalue_vs = new TH2F("KFpvalue_vs", "KFpvalue_vs", 500, 0, 1, 500, 0, 1 );
  KFpvalue_vs->SetXTitle("P-value (#Sigma^{+} mode)");
  KFpvalue_vs->SetYTitle("P-value (#Sigma^{-} mode)");
  
  TH2F *KFchi2ndf_vs = new TH2F("KFchi2ndf_vs", "KFchi2ndf_vs", 100, 0, 100, 100, 0, 100 );
  KFchi2ndf_vs->SetXTitle("chi2/NDF (#Sigma^{+} mode)");
  KFchi2ndf_vs->SetYTitle("chi2/NDF (#Sigma^{-} mode)");

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
  
  //centering title of all histograms 
  TIter nexthist(gDirectory->GetList());
  TH1F *h1 = nullptr;
  TH2F *h2 = nullptr;
  TObject *obj = nullptr;
  while( (obj = (TObject*)nexthist())!=nullptr  ){
    if(obj->InheritsFrom("TH1")){
      h1 = (TH1F*) obj;
      h1->GetXaxis()->CenterTitle();
    }
    if(obj->InheritsFrom("TH2")){
      h2 = (TH2F*) obj;
      h2->GetXaxis()->CenterTitle();
      h2->GetYaxis()->CenterTitle();
    }
  }

  Int_t nevent = tree->GetEntries();
  std::cerr<<"# of events = "<<nevent<<std::endl;
  //------------------------//
  //--- event roop start ---//
  //------------------------//
  for ( Int_t i=0; i<nevent; i++ ) {
    tree->GetEvent(i);
    if(i%10000==0) std::cout << "Event# " << i << std::endl; 
    
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
    
    // calc missing Sm //
    TLorentzVector LVec_pim_n_miss = *LVec_target+*LVec_beam-*LVec_pip-*LVec_n;

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
    //if(anacuts::dE_MIN<dE) NdEOK=true;
    if(dEcut<dE) NdEOK=true;
   
    //Sigma+ production in CDS
    if( (anacuts::Sigmap_MIN<(*LVec_n+*LVec_pip).M() && (*LVec_n+*LVec_pip).M()<anacuts::Sigmap_MAX)) SigmaPFlag=true;
        
    //Sigma- production in CDS
    if( (anacuts::Sigmam_MIN<(*LVec_n+*LVec_pim).M() && (*LVec_n+*LVec_pim).M()<anacuts::Sigmam_MAX)) SigmaMFlag=true;
    
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
      if(NBetaOK){
        dE_MMom_fid_beta_woK0->Fill(LVec_n_miss.P(),dE);
        dE_MMass_fid_beta_woK0->Fill(LVec_n_miss.M(),dE);
        dE_IMnpim_woK0->Fill(LVec_pim_n.M(),dE);
        dE_IMnpip_woK0->Fill(LVec_pip_n.M(),dE);
        if(MissNFlag){
          dE_IMnpim_woK0_n->Fill(LVec_pim_n.M(),dE);
          dE_IMnpip_woK0_n->Fill(LVec_pip_n.M(),dE);
        }
      }
      if(NBetaOK && NdEOK){
        MMom_MMass_fid_beta_dE_woK0->Fill(LVec_n_miss.M(),LVec_n_miss.P());
        IMnpim_IMnpip_dE_woK0->Fill(LVec_pip_n.M(),LVec_pim_n.M());
        nmom->Fill((*LVec_n).P());
        dE_nmom->Fill((*LVec_n).P(),dE);
        mnmom->Fill(nmiss_mom);
        npipmom->Fill(LVec_pip_n.P());
        npimmom->Fill(LVec_pim_n.P());
        if(SigmaPFlag || SigmaMFlag){
          MMom_MMass_fid_beta_dE_woK0_wSid->Fill(LVec_n_miss.M(),LVec_n_miss.P());
        }
      }
      if(NBetaOK && NdEOK && MissNFlag){
        IMnpim_IMnpip_dE_woK0_n->Fill(LVec_pip_n.M(),LVec_pim_n.M());
        if(SigmaPFlag || SigmaMFlag){
          MMnpip_MMnpim_woK0_wSid_n->Fill(LVec_pim_n_miss.M(),LVec_pip_n_miss.M());
          dE_IMnpipi_woK0_wSid_n->Fill(LVec_pip_pim_n.M(),dE);
          Cosn_IMnpipi_woK0_wSid_n->Fill(LVec_pip_pim_n.M(),cos_n);
          MMnmiss_IMnpipi_woK0_wSid_n->Fill(LVec_pip_pim_n.M(), nmiss_mass);
          q_IMnpipi_woK0_wSid_n->Fill(LVec_pip_pim_n.M(),qkn.P());
          nmom_IMnpipi_woK0_wSid_n->Fill(LVec_pip_pim_n.M(),(*LVec_n).P());
        }
      }
    }

    // w/ kinfit
    if( -1<kf_flag && pvalcut <pvalue && K0rejectFlag ){
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
        if(NBetaOK && NdEOK){// && MissNFlag){
          MMnpip_MMnpim_woK0_kin[0]->Fill(LVec_pim_n_miss.M(),LVec_pip_n_miss.M());
          dE_IMnpipi_woK0_kin[0]->Fill(LVec_pip_pim_n.M(),dE);
          Cosn_IMnpipi_woK0_kin[0]->Fill(LVec_pip_pim_n.M(),cos_n);
          MMnmiss_IMnpipi_woK0_kin[0]->Fill(LVec_pip_pim_n.M(), nmiss_mass);
          q_IMnpipi_woK0_kin[0]->Fill(LVec_pip_pim_n.M(),qkn.P());
          nmom_IMnpipi_woK0_kin[0]->Fill(LVec_pip_pim_n.M(),(*LVec_n).P());
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
        if(NBetaOK && NdEOK){ //&& MissNFlag){
          MMnpip_MMnpim_woK0_kin[1]->Fill(LVec_pim_n_miss.M(),LVec_pip_n_miss.M());
          dE_IMnpipi_woK0_kin[1]->Fill(LVec_pip_pim_n.M(),dE);
          Cosn_IMnpipi_woK0_kin[1]->Fill(LVec_pip_pim_n.M(),cos_n);
          MMnmiss_IMnpipi_woK0_kin[1]->Fill(LVec_pip_pim_n.M(), nmiss_mass);
          q_IMnpipi_woK0_kin[1]->Fill(LVec_pip_pim_n.M(),qkn.P());
          nmom_IMnpipi_woK0_kin[1]->Fill(LVec_pip_pim_n.M(),(*LVec_n).P());
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

  TCanvas *cKFpvalue = new TCanvas(Form("cKFpvalue"),"KFpvalue");
  KFpvalue_vs->RebinX(5);
  KFpvalue_vs->RebinY(5);
  KFpvalue_vs->Draw("colz");


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
  
  TCanvas *c_IMnpipi_woK0 = new TCanvas("c_IMnpipi_woK0","c_IMnpipi_woK0"); 
  c_IMnpipi_woK0->cd();
  TH1D *pxSbe = q_IMnpipi_woK0_wSid_n->ProjectionX();
  TH1D *pxSp = q_IMnpipi_woK0_kin[0]->ProjectionX();
  TH1D *pxSm = q_IMnpipi_woK0_kin[1]->ProjectionX();
  pxSbe->Rebin(2);
  pxSbe->Draw("HE");
  pxSm->SetLineColor(2);
  pxSm->Rebin(2);
  pxSm->Draw("HEsame");
  pxSp->Rebin(2);
  pxSp->SetLineColor(3);
  pxSp->Draw("HEsame");
  TH1D *pxSum = (TH1D*) pxSp->Clone();
  pxSum->Add(pxSm);
  pxSum->SetLineColor(4);
  pxSum->Draw("HEsame");


  TCanvas *cdE_betainv_fid = new TCanvas(Form("cdE_betainv_fid"), "");
  cdE_betainv_fid->cd();
  dE_betainv_fid->Draw("colz");

  
  TCanvas  *cIMnpim_IMnpip_dE_woK0 = new TCanvas(Form("cIMnpim_IMnpip_dE_woK0"),"");
  cIMnpim_IMnpip_dE_woK0->cd();
  IMnpim_IMnpip_dE_woK0->Draw("colz");
  
  /*
  TCanvas  *cIMnpim_IMnpip_dE_woK0_px = new TCanvas(Form("cIMnpim_IMnpip_dE_woK0_px"),"");
  cIMnpim_IMnpip_dE_woK0_px->cd();
  TH1D* IMnpim_IMnpip_dE_woK0_px = IMnpim_IMnpip_dE_woK0->ProjectionX();
  IMnpim_IMnpip_dE_woK0_px->Draw("");
  
  TCanvas  *cIMnpim_IMnpip_dE_woK0_py = new TCanvas(Form("cIMnpim_IMnpip_dE_woK0_py"),"");
  cIMnpim_IMnpip_dE_woK0_py->cd();
  TH1D* IMnpim_IMnpip_dE_woK0_py = IMnpim_IMnpip_dE_woK0->ProjectionY();
  IMnpim_IMnpip_dE_woK0_py->Draw("");
  
  */
  TCanvas *cnmom_IMnpipi_woK0_wSid_n = new TCanvas(Form("cnmom_IMnpipi_woK0_wSid_n"),"");
  cnmom_IMnpipi_woK0_wSid_n->cd();
  nmom_IMnpipi_woK0_wSid_n->Draw("colz");
  
  TCanvas *cnmom_IMnpipi_woK0_wSid_n_py = new TCanvas(Form("cnmom_IMnpipi_woK0_wSid_n_py"),"");
  cnmom_IMnpipi_woK0_wSid_n_py->cd();
  TH1D* nmom_IMnpipi_woK0_wSid_n_py = nmom_IMnpipi_woK0_wSid_n->ProjectionY();
  nmom_IMnpipi_woK0_wSid_n_py->Draw();

  TCanvas *cMMnmiss_IMnpipi_woK0_wSid_n = new TCanvas(Form("cMMnmiss_IMnpipi_woK0_wSid_n"),"");
  cMMnmiss_IMnpipi_woK0_wSid_n->cd();
  MMnmiss_IMnpipi_woK0_wSid_n->Draw("colz");

  TCanvas *cMMnmiss_IMnpipi_woK0_wSid_n_px = new TCanvas(Form("cMMnmiss_IMnpipi_woK0_wSid_n_px"),""); 
  cMMnmiss_IMnpipi_woK0_wSid_n_px->cd();
  MMnmiss_IMnpipi_woK0_wSid_n->ProjectionX()->Draw("");
  TCanvas *cMMnmiss_IMnpipi_woK0_wSid_n_py = new TCanvas(Form("cMMnmiss_IMnpipi_woK0_wSid_n_py"),""); 
  cMMnmiss_IMnpipi_woK0_wSid_n_py->cd();
  MMnmiss_IMnpipi_woK0_wSid_n->ProjectionY()->Draw("");


  TCanvas *cdE_betainv_fid_kin[2];//           
  TCanvas *cdE_MMom_fid_beta_woK0_kin[2];      
  TCanvas *cdE_MMass_fid_beta_woK0_kin[2];     
  TCanvas *cMMom_MMass_fid_beta_dE_woK0_kin[2];
  TCanvas *cMMom_MMass_fid_beta_dE_woK0_kin_px[2];
  TH1D *MMom_MMass_fid_beta_dE_woK0_kin_px[2];
  TCanvas *cMMom_MMass_fid_beta_dE_woK0_kin_py[2];
  TH1D *MMom_MMass_fid_beta_dE_woK0_kin_py[2];
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
    cdE_betainv_fid_kin[imode] = new TCanvas(Form("cdE_betainv_fid_kin_%s",smode[imode]), "");
    cdE_betainv_fid_kin[imode]->cd();
    dE_betainv_fid_kin[imode]->Draw("colz");


    cdE_MMom_fid_beta_woK0_kin[imode] = new TCanvas(Form("cdE_MMom_fid_beta_woK0_kin_%s",smode[imode]),"");
    cdE_MMom_fid_beta_woK0_kin[imode]->cd();
    dE_MMom_fid_beta_woK0_kin[imode]->Draw("colz");

    cdE_MMass_fid_beta_woK0_kin[imode] = new TCanvas(Form("cdE_MMass_fid_beta_woK0_%s_kin",smode[imode]),"");
    cdE_MMass_fid_beta_woK0_kin[imode]->cd(); 
    dE_MMass_fid_beta_woK0_kin[imode]->Draw("colz");

    cMMom_MMass_fid_beta_dE_woK0_kin[imode] = new TCanvas(Form("cMMom_MMass_fid_beta_dE_woK0_kin_%s",smode[imode]),"");
    cMMom_MMass_fid_beta_dE_woK0_kin[imode]->cd();
    MMom_MMass_fid_beta_dE_woK0_kin[imode]->Draw("colz");

    cMMom_MMass_fid_beta_dE_woK0_kin_px[imode] = new TCanvas(Form("cMMom_MMass_fid_beta_dE_woK0_kin_px_%s",smode[imode]),"");
    cMMom_MMass_fid_beta_dE_woK0_kin_px[imode]->cd();
    MMom_MMass_fid_beta_dE_woK0_kin_px[imode] = MMom_MMass_fid_beta_dE_woK0_kin[imode]->ProjectionX();
    MMom_MMass_fid_beta_dE_woK0_kin_px[imode]->Draw();

    cMMom_MMass_fid_beta_dE_woK0_kin_py[imode] = new TCanvas(Form("cMMom_MMass_fid_beta_dE_woK0_kin_py_%s",smode[imode]),"");
    cMMom_MMass_fid_beta_dE_woK0_kin_py[imode]->cd();
    MMom_MMass_fid_beta_dE_woK0_kin_py[imode] = MMom_MMass_fid_beta_dE_woK0_kin[imode]->ProjectionY();
    MMom_MMass_fid_beta_dE_woK0_kin_py[imode]->Draw();
    

    cIMnpim_IMnpip_dE_woK0_kin[imode] = new TCanvas(Form("cIMnpim_IMnpip_dE_woK0_kin_%s",smode[imode]),"");
    cIMnpim_IMnpip_dE_woK0_kin[imode]->cd();
    IMnpim_IMnpip_dE_woK0_kin[imode]->Draw("colz");
    
    cIMnpim_IMnpip_dE_woK0_kin_px[imode] = new TCanvas(Form("cIMnpim_IMnpip_dE_woK0_kin_px_%s",smode[imode]),"");
    cIMnpim_IMnpip_dE_woK0_kin_px[imode]->cd();
    IMnpim_IMnpip_dE_woK0_kin_px[imode] = IMnpim_IMnpip_dE_woK0_kin[imode]->ProjectionX();
    IMnpim_IMnpip_dE_woK0_kin_px[imode]->Draw();
    
    cIMnpim_IMnpip_dE_woK0_kin_py[imode] = new TCanvas(Form("cIMnpim_IMnpip_dE_woK0_kin_py_%s",smode[imode]),"");
    cIMnpim_IMnpip_dE_woK0_kin_py[imode]->cd();
    IMnpim_IMnpip_dE_woK0_kin_py[imode] = IMnpim_IMnpip_dE_woK0_kin[imode]->ProjectionY();
    IMnpim_IMnpip_dE_woK0_kin_py[imode]->Draw();
    
    cMMnpip_MMnpim_woK0_kin[imode] = new TCanvas(Form("cMMnpip_MMnpim_woK0_kin_%s",smode[imode]),"");
    cMMnpip_MMnpim_woK0_kin[imode]->cd();
    MMnpip_MMnpim_woK0_kin[imode]->Draw("colz");

    cdE_IMnpipi_woK0_kin[imode] = new TCanvas(Form("cdE_IMnpipi_woK0_kin_%s",smode[imode]),"");
    cdE_IMnpipi_woK0_kin[imode]->cd();
    dE_IMnpipi_woK0_kin[imode]->Draw("colz");
    
    cCosn_IMnpipi_n_kin[imode] = new TCanvas(Form("cCosn_IMnpipi_n_kin_%s",smode[imode]),"");
    cCosn_IMnpipi_n_kin[imode]->cd();
    Cosn_IMnpipi_woK0_kin[imode]->Draw("colz");

    cMMnmiss_IMnpipi_woK0_kin[imode] = new TCanvas(Form("cMMnmiss_IMnpipi_woK0_kin_%s",smode[imode]),"");
    cMMnmiss_IMnpipi_woK0_kin[imode]->cd();
    MMnmiss_IMnpipi_woK0_kin[imode]->Draw("colz");

    cq_IMnpipi_woK0_kin[imode] = new TCanvas(Form("cq_IMnpipi_woK0_kin_%s",smode[imode]),"");
    cq_IMnpipi_woK0_kin[imode]->cd();
    q_IMnpipi_woK0_kin[imode]->RebinX(2);
    q_IMnpipi_woK0_kin[imode]->RebinY(6);
    q_IMnpipi_woK0_kin[imode]->Draw("colz");
    
    cnmom_IMnpipi_woK0_kin[imode] = new TCanvas(Form("cnmom_IMnpipi_woK0_kin_%s",smode[imode]),"");
    cnmom_IMnpipi_woK0_kin[imode]->cd();
    nmom_IMnpipi_woK0_kin[imode]->RebinX(2);
    //nmom_IMnpipi_woK0_kin[imode]->RebinY(6);
    nmom_IMnpipi_woK0_kin[imode]->Draw("colz");
  }

  //TCanvas *cMMnpip_MMnpim_woK0_wSid_n = new TCanvas(Form("cMMnpip_MMnpim_woK0_wSid_n"),"");
  //cMMnpip_MMnpim_woK0_wSid_n->cd();
  //MMnpip_MMnpim_woK0_wSid_n->Draw("colz");
  
  TCanvas *cMMom_MMass_fid_beta_dE_woK0  = new TCanvas(Form("cMMom_MMass_fid_beta_dE_woK0"),"");
  cMMom_MMass_fid_beta_dE_woK0->cd();
  MMom_MMass_fid_beta_dE_woK0->Draw("colz");
  
  TCanvas *cMMom_MMass_fid_beta_dE_woK0_px = new TCanvas(Form("cMMom_MMass_fid_beta_dE_woK0_px"),"");
  TH1D *MMom_MMass_fid_beta_dE_woK0_px = MMom_MMass_fid_beta_dE_woK0->ProjectionX();
  MMom_MMass_fid_beta_dE_woK0_px->Draw("");
  
  TCanvas *cMMom_MMass_fid_beta_dE_woK0_wSid  = new TCanvas(Form("cMMom_MMass_fid_beta_dE_woK0_wSid"),"");
  cMMom_MMass_fid_beta_dE_woK0_wSid->cd();
  MMom_MMass_fid_beta_dE_woK0_wSid->Draw("colz");
  

  TCanvas *cMMom_MMass_fid_beta_dE_woK0_wSid_px  = new TCanvas(Form("cMMom_MMass_fid_beta_dE_woK0_wSid_px"),"");
  cMMom_MMass_fid_beta_dE_woK0_wSid_px->cd();
  TH1D *MMom_MMass_fid_beta_dE_woK0_wSid_px = MMom_MMass_fid_beta_dE_woK0_wSid->ProjectionX();
  MMom_MMass_fid_beta_dE_woK0_wSid_px->Draw("");


  TCanvas *cMMom_MMass_fid_beta_dE_woK0_wSid_py  = new TCanvas(Form("cMMom_MMass_fid_beta_dE_woK0_wSid_py"),"");
  cMMom_MMass_fid_beta_dE_woK0_wSid_py->cd();
  TH1D *MMom_MMass_fid_beta_dE_woK0_wSid_py = MMom_MMass_fid_beta_dE_woK0_wSid->ProjectionY();
  MMom_MMass_fid_beta_dE_woK0_wSid_py->Draw("");
   

  TCanvas *cMMom_MMass_fid_beta_dE_woK0_px_sup = new TCanvas("cMMom_MMass_fid_beta_dE_woK0_px_sup","");
  cMMom_MMass_fid_beta_dE_woK0_px_sup->cd();
  MMom_MMass_fid_beta_dE_woK0_px->Draw();
  TH1D *MMom_MMass_fid_beta_dE_woK0_wSid_px_clone = (TH1D*)MMom_MMass_fid_beta_dE_woK0_wSid_px->Clone(); 
  MMom_MMass_fid_beta_dE_woK0_wSid_px_clone->SetLineColor(4);
  cMMom_MMass_fid_beta_dE_woK0_px_sup->cd();
  MMom_MMass_fid_beta_dE_woK0_wSid_px_clone->Draw("same");
  TH1D *MMom_MMass_fid_beta_dE_woK0_kin_px_clone[2];
  MMom_MMass_fid_beta_dE_woK0_kin_px_clone[0] = (TH1D*)MMom_MMass_fid_beta_dE_woK0_kin_px[0]->Clone();
  MMom_MMass_fid_beta_dE_woK0_kin_px_clone[0]->SetLineColor(2);
  MMom_MMass_fid_beta_dE_woK0_kin_px_clone[0]->Draw("same");
  MMom_MMass_fid_beta_dE_woK0_kin_px_clone[1] = (TH1D*)MMom_MMass_fid_beta_dE_woK0_kin_px[1]->Clone();
  MMom_MMass_fid_beta_dE_woK0_kin_px_clone[1]->SetLineColor(3);
  MMom_MMass_fid_beta_dE_woK0_kin_px_clone[1]->Draw("same");

  TCanvas *cMMom_MMass_fid_beta_dE_woK0_py = new TCanvas("cMMom_MMass_fid_beta_dE_woK0_py","");
  cMMom_MMass_fid_beta_dE_woK0_py->cd();
  TH1D *MMom_MMass_fid_beta_dE_woK0_py = MMom_MMass_fid_beta_dE_woK0->ProjectionY();
  MMom_MMass_fid_beta_dE_woK0_py->Draw();
  //test
  TH1D *MMom_MMass_fid_beta_dE_woK0_wSid_py_clone = (TH1D*)MMom_MMass_fid_beta_dE_woK0_wSid_py->Clone();
  MMom_MMass_fid_beta_dE_woK0_wSid_py_clone->SetLineColor(4);
  cMMom_MMass_fid_beta_dE_woK0_py->cd();
  MMom_MMass_fid_beta_dE_woK0_wSid_py_clone->Draw("same");
  TH1D* MMom_MMass_fid_beta_dE_woK0_kin_py_clone[2];
  MMom_MMass_fid_beta_dE_woK0_kin_py_clone[0] = (TH1D*) MMom_MMass_fid_beta_dE_woK0_kin_py[0]->Clone();
  MMom_MMass_fid_beta_dE_woK0_kin_py_clone[0]->SetLineColor(2);
  MMom_MMass_fid_beta_dE_woK0_kin_py_clone[0]->Draw("same");
  MMom_MMass_fid_beta_dE_woK0_kin_py_clone[1] = (TH1D*) MMom_MMass_fid_beta_dE_woK0_kin_py[1]->Clone();
  MMom_MMass_fid_beta_dE_woK0_kin_py_clone[1]->SetLineColor(3);
  MMom_MMass_fid_beta_dE_woK0_kin_py_clone[1]->Draw("same");


  
  TCanvas *cIMnpim_IMnpip_dE_woK0_n = new TCanvas(Form("cIMnpim_IMnpip_dE_woK0_n"),"");
  cIMnpim_IMnpip_dE_woK0_n->cd();
  IMnpim_IMnpip_dE_woK0_n->Draw("colz");
  
  TCanvas *cIMnpim_IMnpip_dE_woK0_px = new TCanvas(Form("cIMnpim_IMnpip_dE_woK0_px"),"");
  cIMnpim_IMnpip_dE_woK0_px->cd();
  TH1D *IMnpim_IMnpip_dE_woK0_px = IMnpim_IMnpip_dE_woK0->ProjectionX();
  IMnpim_IMnpip_dE_woK0_px->Draw();
  //TH1D * IMnpim_IMnpip_dE_woK0_n_px_sum =(TH1D*) IMnpim_IMnpip_dE_woK0_kin_px[0]->Clone();
  
  TCanvas *cIMnpim_IMnpip_dE_woK0_py = new TCanvas(Form("cIMnpim_IMnpip_dE_woK0_py"),"");
  cIMnpim_IMnpip_dE_woK0_py->cd();
  TH1D *IMnpim_IMnpip_dE_woK0_py = IMnpim_IMnpip_dE_woK0->ProjectionY();
  IMnpim_IMnpip_dE_woK0_py->Draw();

  TCanvas *cIMnpim_IMnpip_dE_woK0_n_px = new TCanvas(Form("cIMnpim_IMnpip_dE_woK0_n_px"),"");
  cIMnpim_IMnpip_dE_woK0_n_px->cd();
  TH1D *IMnpim_IMnpip_dE_woK0_n_px = IMnpim_IMnpip_dE_woK0_n->ProjectionX();
  IMnpim_IMnpip_dE_woK0_n_px->Draw();
  TH1D* IMnpim_IMnpip_dE_woK0_kin_px_clone[2];
  IMnpim_IMnpip_dE_woK0_kin_px_clone[0] = (TH1D*) IMnpim_IMnpip_dE_woK0_kin_px[0]->Clone();
  IMnpim_IMnpip_dE_woK0_kin_px_clone[0]->SetLineColor(2);
  IMnpim_IMnpip_dE_woK0_kin_px_clone[0]->Draw("same");
  IMnpim_IMnpip_dE_woK0_kin_px_clone[1] = (TH1D*) IMnpim_IMnpip_dE_woK0_kin_px[1]->Clone();
  IMnpim_IMnpip_dE_woK0_kin_px_clone[1]->SetLineColor(3);
  IMnpim_IMnpip_dE_woK0_kin_px_clone[1]->Draw("same");
  //TH1D * IMnpim_IMnpip_dE_woK0_n_px_sum =(TH1D*) IMnpim_IMnpip_dE_woK0_kin_px[0]->Clone();
  
  TCanvas *cIMnpim_IMnpip_dE_woK0_n_py = new TCanvas(Form("cIMnpim_IMnpip_dE_woK0_n_py"),"");
  cIMnpim_IMnpip_dE_woK0_n_py->cd();
  TH1D *IMnpim_IMnpip_dE_woK0_n_py = IMnpim_IMnpip_dE_woK0_n->ProjectionY();
  IMnpim_IMnpip_dE_woK0_n_py->Draw();
  TH1D *IMnpim_IMnpip_dE_woK0_kin_py_clone[2];
  IMnpim_IMnpip_dE_woK0_kin_py_clone[0] = (TH1D*)IMnpim_IMnpip_dE_woK0_kin_py[0]->Clone();
  IMnpim_IMnpip_dE_woK0_kin_py_clone[0]->SetLineColor(2);
  IMnpim_IMnpip_dE_woK0_kin_py_clone[0]->Draw("same");
  IMnpim_IMnpip_dE_woK0_kin_py_clone[1] = (TH1D*)IMnpim_IMnpip_dE_woK0_kin_py[1]->Clone();
  IMnpim_IMnpip_dE_woK0_kin_py_clone[1]->SetLineColor(3);
  IMnpim_IMnpip_dE_woK0_kin_py_clone[1]->Draw("same");


  TFile *fout = new TFile(outfilename.c_str(),"RECREATE");
  fout->cd();

  /*
  for(int imode=0;imode<2;imode++){
    dE_betainv_fid_kin[imode]->Write();
    dE_MMom_fid_beta_woK0_kin[imode]->Write();
    dE_MMass_fid_beta_woK0_kin[imode]->Write();
    MMom_MMass_fid_beta_dE_woK0_kin[imode]->Write();
    IMnpim_IMnpip_dE_woK0_kin[imode]->Write();
    MMnpip_MMnpim_woK0_kin[imode]->Write();
    dE_IMnpipi_woK0_kin[imode]->Write();
    Cosn_IMnpipi_woK0_kin[imode]->Write();
    MMnmiss_IMnpipi_woK0_kin[imode]->Write();
    q_IMnpipi_woK0_kin[imode]->Write();
    nmom_IMnpipi_woK0_kin[imode]->Write();
  };

  q_IMnpipi_woK0_wSid_n->Write();
  nmom->Write();
  mnmom->Write();
  npipmom->Write();
  npimmom->Write();*/
  
  TCanvas *c = nullptr;
  TPDF *pdf = new TPDF(pdfname);
  TIter next(gROOT->GetListOfCanvases());
  while((c= (TCanvas*)next())){
    pdf->NewPage();
    c->Draw();
    c->cd();
    TPaveText *pt = new TPaveText(.74,.81,0.9,0.90,"NDC");
    pt->AddText("Real Data");
    pt->SetFillColor(kCyan-9);
    pt->SetBorderSize(1);
    pt->Draw();
    c->Modified();
    c->Update();
  }
  pdf->Close();
  std::cout << "closing pdf " << std::endl;
  
  fout->Close();
  
}
