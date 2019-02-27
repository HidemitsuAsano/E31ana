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
#include <TPaveText.h>

#include "../src/IMPiSigmaAnaPar.h"

const double pvalcut = 0.005;
const double dEcut = 2.0;
//const double pvalcut = 1.0e-5;
const bool gridon=true;
const bool staton=true;

//mode 0: Sigma+ ,1: Sigma- 
void plot_IMpisigma(const char* filename="")
{

  //gROOT->SetStyle("Plain");
  if(staton)gStyle->SetOptStat("emruo");
  else gStyle->SetOptStat(0);
  gStyle->SetOptFit(111111);
  gStyle->SetPadGridX(gridon);
  gStyle->SetPadGridY(gridon);
  //gStyle->SetStatX(0.98);
  gStyle->SetStatX(0.9);
  gStyle->SetStatY(0.9);
  gStyle->SetStatBorderSize(1);

  std::cout << "infile " << filename <<std::endl;
  //std::string outfilename = std::string(filename);
  //outfilename.insert(outfilename.size()-5,"_post");
  //std::cout << "outfilename: " << outfilename << std::endl;
  TString pdfname = std::string(filename);
  pdfname.Replace(std::string(filename).size()-4,5,"pdf");
  std::cout << "pdfname: " << pdfname << std::endl;

  bool Spmode = (std::string(filename).find("Sp")!= std::string::npos);
  if(Spmode){
    std::cout << "This is Sigma+ mode sim." << std::endl;
  }else{
    std::cout << "This is Sigma- mode sim." << std::endl;
  }
  
  //= = = = pipipnn final-sample tree = = = =//
  TLorentzVector *LVec_beam=nullptr;   // 4-momentum(beam)
  TLorentzVector *LVec_target=nullptr; // 4-momentum(target)
  TLorentzVector *LVec_pip=nullptr;    // 4-momentum(pi+)
  TLorentzVector *LVec_pim=nullptr;    // 4-momentum(pi-)
  TLorentzVector *LVec_n=nullptr;      // 4-momentum(neutron)
  TLorentzVector *mcmom_beam=nullptr;   // generated 4-momentum(beam)
  TLorentzVector *mcmom_pip=nullptr;    // generated 4-momentum(pi+)
  TLorentzVector *mcmom_pim=nullptr;    // generated 4-momentum(pi-)
  TLorentzVector *mcmom_ncds=nullptr;      // generated 4-momentum(neutron)
  TLorentzVector *mcmom_nmiss=nullptr;      // generated 4-momentum(neutron)
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
  tree->SetBranchAddress( "NeutralBetaCDH", &NeutralBetaCDH );//<- from v32.
  //tree->SetBranchAddress( "beta", &NeutralBetaCDH );
  tree->SetBranchAddress( "dE", &dE );
  tree->SetBranchAddress( "vtx_reaction", &vtx_reaction );
  tree->SetBranchAddress( "vtx_pip",&vtx_pip);
  tree->SetBranchAddress( "vtx_pim",&vtx_pim);
  //tree->SetBranchAddress( "run_num", &run_num );
  //tree->SetBranchAddress( "event_num", &event_num );
  //tree->SetBranchAddress( "block_num", &block_num );
  tree->SetBranchAddress( "mcmom_beam",  &mcmom_beam );
  tree->SetBranchAddress( "mcmom_pip", &mcmom_pip);
  tree->SetBranchAddress( "mcmom_pim", &mcmom_pim);
  tree->SetBranchAddress( "mcmom_ncds", &mcmom_ncds);
  tree->SetBranchAddress( "mcmom_nmiss", &mcmom_nmiss);
  tree->SetBranchAddress( "kfSpmode_mom_beam",   &kfSpmode_mom_beam );
  tree->SetBranchAddress( "kfSpmode_mom_pip", &kfSpmode_mom_pip );
  tree->SetBranchAddress( "kfSpmode_mom_pim", &kfSpmode_mom_pim );
  tree->SetBranchAddress( "kfSpmode_mom_n", &kfSpmode_mom_n );
  //tree->SetBranchAddress( "kf_chi2_Spmode", &kfSpmode_chi2 );     
  //tree->SetBranchAddress( "kf_NDF_Spmode", &kfSpmode_NDF );       
  //tree->SetBranchAddress( "kf_status_Spmode", &kfSpmode_status ); 
  //tree->SetBranchAddress( "kf_pvalue_Spmode", &kfSpmode_pvalue ); 
  tree->SetBranchAddress( "kfSpmode_chi2", &kfSpmode_chi2 );
  tree->SetBranchAddress( "kfSpmode_NDF", &kfSpmode_NDF );
  tree->SetBranchAddress( "kfSpmode_status", &kfSpmode_status );
  tree->SetBranchAddress( "kfSpmode_pvalue", &kfSpmode_pvalue );
  tree->SetBranchAddress( "kfSmmode_mom_beam",   &kfSmmode_mom_beam );
  tree->SetBranchAddress( "kfSmmode_mom_pip", &kfSmmode_mom_pip );
  tree->SetBranchAddress( "kfSmmode_mom_pim", &kfSmmode_mom_pim );
  tree->SetBranchAddress( "kfSmmode_mom_n", &kfSmmode_mom_n );
  //tree->SetBranchAddress( "kf_chi2_Smmode", &kfSmmode_chi2 );       
  //tree->SetBranchAddress( "kf_NDF_Smmode", &kfSmmode_NDF );         
  //tree->SetBranchAddress( "kf_status_Smmode", &kfSmmode_status );   
  //tree->SetBranchAddress( "kf_pvalue_Smmode", &kfSmmode_pvalue );   
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
  TH2F* nmom_IMnpipi_woK0_wSid_n;\
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
  TH2F* q_IMnpipi_kin[2];
  TH2F* q_IMnpipi_woK0_kin[2];
  TH2F* nmom_IMnpipi_woK0_kin[2];
  const char smode[][4]={"Sp","Sm"};
   
  /*
  //correlation plots again
  const int BIN = 4000;
  const double cov_MAX=0.4;
  TH1F* cov[kin::npart][4][4];
  TH2F* cov_mom[kin::npart][4][4];
  for( int ip=0; ip<kin::npart; ip++ ){
    for( int ii=0; ii<4; ii++ ){
      for( int jj=0; jj<4; jj++ ){
        TH1F* cov[ip][ii][jj] = new TH1F(Form("cov_%d_%d_%d", ip, ii, jj), BIN, -cov_MAX, cov_MAX);
        TH2F* cov_mom[ip][ii][jj] = new TH2F(Form("cov_mom_%d_%d_%d", ip, ii, jj), 500., -cov_MAX, cov_MAX,200,0,2.0);
      }
    }
  }//for ip
  */ 

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

  MMnmiss_IMnpip_dE_woK0 = new TH2F("MMnmiss_IMnpip_dE_woK0", "MMnmiss_IMnpip_dE_woK0",200,1.,2.0,100,0,1.5);
  MMnmiss_IMnpip_dE_woK0->SetXTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  MMnmiss_IMnpip_dE_woK0->SetYTitle("Miss Mass. [GeV/c^{2}]");

  MMnmiss_IMnpim_dE_woK0 = new TH2F("MMnmiss_IMnpim_dE_woK0", "MMnmiss_IMnpim_dE_woK0",200,1.,2.0,100,0,1.5);
  MMnmiss_IMnpim_dE_woK0->SetXTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  MMnmiss_IMnpim_dE_woK0->SetYTitle("Miss Mass. [GeV/c^{2}]");
  
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
    
  q_IMnpipi_wSid_n = new TH2F(Form("q_IMnpipi_wSid_n"),Form("q_IMnpipi_wSid_n"),100,1,2,300,0,1.5);
  q_IMnpipi_wSid_n->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
  q_IMnpipi_wSid_n->SetYTitle("Mom. Transfer [GeV/c]");
  
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

    q_IMnpipi_kin[imode] = new TH2F(Form("q_IMnpipi_kin_%s",smode[imode]),Form("q_IMnpipi_kin_%s",smode[imode]),100,1,2,300,0,1.5);
    q_IMnpipi_kin[imode]->SetXTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
    q_IMnpipi_kin[imode]->SetYTitle("Mom. Transfer [GeV/c]");
    
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
  
  //centering title of all histograms 
  TIter nexthist(gDirectory->GetList());
  TH1F *h1 = nullptr;
  TH2F *h2 = nullptr;
  TObject *obj = nullptr;
  while( (obj = (TObject*)nexthist())!=nullptr ){
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
  std::cout<<"# of events = "<<nevent<<std::endl;
  
  std::cout << "p-value cut:" << pvalcut << std::endl; 
  std::cout << "dE cut:" << dEcut << std::endl; 
  TCanvas *cinfo = new TCanvas("cinfo","info");
  TPaveText *pt = new TPaveText(.05,.05,.95,.7);
  pt->AddText(Form("p-value cut: %0.3f ",pvalcut));
  pt->AddText(Form("dE cut: %0.3f " ,dEcut));
  pt->AddText(Form("1/beta min.: %0.7f ",1./anacuts::beta_MAX));
  pt->AddText(Form("K0 window : %0.3f - %0.3f",anacuts::pipi_MIN,anacuts::pipi_MAX )); 
  pt->AddText(Form("Sp window : %0.3f - %0.3f",anacuts::Sigmap_MIN,anacuts::Sigmap_MAX )); 
  pt->AddText(Form("Sm window : %0.3f - %0.3f",anacuts::Sigmam_MIN,anacuts::Sigmam_MAX )); 
  pt->AddText(Form("# of events %d",nevent));
  pt->Draw(); 
  //------------------------//
  //--- event roop start ---//
  //------------------------//
  for ( Int_t i=0; i<nevent; i++ ) {
    tree->GetEvent(i);
    if(i%50000==0) std::cout << "Event# " << i << std::endl; 
    
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
      MMom_MMass_fid_beta_dE_woK0->Fill(LVec_n_miss.M(),LVec_n_miss.P());
      IMnpim_IMnpip_dE_woK0->Fill(LVec_pip_n.M(),LVec_pim_n.M());
      nmom->Fill((*LVec_n).P());
      dE_nmom->Fill((*LVec_n).P(),dE);
      mnmom->Fill(nmiss_mom);
      npipmom->Fill(LVec_pip_n.P());
      npimmom->Fill(LVec_pim_n.P());
      MMnmiss_IMnpip_dE_woK0->Fill(LVec_pip_n.M(),nmiss_mass);
      MMnmiss_IMnpim_dE_woK0->Fill(LVec_pim_n.M(),nmiss_mass);
      if(SigmaPFlag || SigmaMFlag){
        MMom_MMass_fid_beta_dE_woK0_wSid->Fill(LVec_n_miss.M(),LVec_n_miss.P());
      }
    }
    if(K0rejectFlag && NBetaOK && NdEOK && MissNFlag){
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
    
    //including K0 
    if(NBetaOK && NdEOK && MissNFlag){
      if(SigmaPFlag || SigmaMFlag){
        q_IMnpipi_wSid_n->Fill(LVec_pip_pim_n.M(),qkn.P());
      }
    }

    // w/ kinfit
    if( (-1<kf_flag) && (pvalcut < pvalue) && K0rejectFlag ){
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
        if(NBetaOK && NdEOK){// && MissNFlag)
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
        if(NBetaOK && NdEOK){ //&& MissNFlag)
          MMnpip_MMnpim_woK0_kin[1]->Fill(LVec_pim_n_miss.M(),LVec_pip_n_miss.M());
          dE_IMnpipi_woK0_kin[1]->Fill(LVec_pip_pim_n.M(),dE);
          Cosn_IMnpipi_woK0_kin[1]->Fill(LVec_pip_pim_n.M(),cos_n);
          MMnmiss_IMnpipi_woK0_kin[1]->Fill(LVec_pip_pim_n.M(), nmiss_mass);
          q_IMnpipi_woK0_kin[1]->Fill(LVec_pip_pim_n.M(),qkn.P());
          nmom_IMnpipi_woK0_kin[1]->Fill(LVec_pip_pim_n.M(),(*LVec_n).P());
        }
      }
    }//pvalcut && K0rejectflag
    
    if( (-1<kf_flag) && (pvalcut < pvalue) ){
      if( kfSmmode_pvalue<kfSpmode_pvalue ){
        if(NBetaOK && NdEOK){ //&& MissNFlag)
          q_IMnpipi_kin[0]->Fill(LVec_pip_pim_n.M(),qkn.P());
        }
      }else{
        if(NBetaOK && NdEOK){ //&& MissNFlag)
          q_IMnpipi_kin[1]->Fill(LVec_pip_pim_n.M(),qkn.P());
        }
      }
    }   
	}//for ievt
   

  TCanvas *cKFpvalue_vs = new TCanvas(Form("cKFpvalue_vs"),"KFpvalue_vs");
  cKFpvalue_vs->cd();
  cKFpvalue_vs->SetLogy();
  TH1D *KFpvalue_vs_px = (TH1D*) KFpvalue_vs->ProjectionX();
  KFpvalue_vs_px->SetMinimum(1);
  KFpvalue_vs_px->SetXTitle("p-value");
  TH1D *KFpvalue_vs_py = (TH1D*) KFpvalue_vs->ProjectionY();
  KFpvalue_vs_py->SetLineColor(2);
  if(Spmode){
    KFpvalue_vs_px->Draw();
    KFpvalue_vs_py->Draw("same");
  }else{
    KFpvalue_vs_py->Draw("");
    KFpvalue_vs_px->Draw("same");
  }
  TLegend *legKFpvalue_vs = new TLegend(0.55,0.65,0.76,0.82);
  legKFpvalue_vs->AddEntry(KFpvalue_vs_px,"#Sigma^{+} mode");
  legKFpvalue_vs->AddEntry(KFpvalue_vs_py,"#Sigma^{-} mode");
  legKFpvalue_vs->SetFillColor(0);
  legKFpvalue_vs->Draw();

  int spbin = KFpvalue_vs_px->FindBin(pvalcut);
  //std::cout << spbin << std::endl;
  int smbin = KFpvalue_vs_py->FindBin(pvalcut);
  std::cout << "Sp mode rejection ratio:" << KFpvalue_vs_px->Integral(0,spbin)/(KFpvalue_vs_px->Integral(0,201)) << std::endl;
  std::cout << "Sm mode rejection ratio:" << KFpvalue_vs_py->Integral(0,smbin)/(KFpvalue_vs_py->Integral(0,201)) << std::endl;

  //cumulative dist. of prob.
  TCanvas *cKFpvalue_vs_cum = new TCanvas(Form("cKFpvalue_vs_cum"),"KFpvalue_vs_cum");
  cKFpvalue_vs_cum->cd();
  TH1 *KFpvalue_vs_px_cum = KFpvalue_vs_px->GetCumulative();
  KFpvalue_vs_px_cum->Scale(1./KFpvalue_vs_px->GetEntries());
  KFpvalue_vs_px_cum->SetXTitle("p-value cut");
  KFpvalue_vs_px_cum->GetXaxis()->CenterTitle();
  KFpvalue_vs_px_cum->SetTitle("KF pvalue cumulative");
  TH1 *KFpvalue_vs_py_cum = KFpvalue_vs_py->GetCumulative();
  KFpvalue_vs_py_cum->SetLineColor(2);
  KFpvalue_vs_py_cum->Scale(1./KFpvalue_vs_py->GetEntries());
  if(Spmode){
    KFpvalue_vs_px_cum->Draw();
    KFpvalue_vs_py_cum->Draw("same");
  }else{
    KFpvalue_vs_py_cum->Draw("");
    KFpvalue_vs_px_cum->Draw("same");
  }

  TLegend *legKFpvalue_vs_cum = new TLegend(0.55,0.25,0.76,0.42);
  legKFpvalue_vs_cum->AddEntry(KFpvalue_vs_px_cum,"#Sigma^{+} mode");
  legKFpvalue_vs_cum->AddEntry(KFpvalue_vs_py_cum,"#Sigma^{-} mode");
  legKFpvalue_vs_cum->SetFillColor(0);
  legKFpvalue_vs_cum->Draw();

  TCanvas *cKFpvalue = new TCanvas(Form("cKFpvalue"),"KFpvalue");
  KFpvalue_vs->RebinX(5);
  KFpvalue_vs->RebinY(5);
  KFpvalue_vs->Draw("colz");
  
  //chi2/ndf
  TCanvas *cKFchi2ndf_vs = new TCanvas(Form("cKFchi2ndf_vs"),"KFchi2ndf_vs");
  cKFchi2ndf_vs->cd();
  //cKFchi2ndf_vs->SetLogy();
  TH1D *KFchi2ndf_vs_px = (TH1D*) KFchi2ndf_vs->ProjectionX();
  KFchi2ndf_vs_px->SetMinimum(1);
  KFchi2ndf_vs_px->GetXaxis()->SetTitle("chi2/ndf");
  TH1D *KFchi2ndf_vs_py = (TH1D*) KFchi2ndf_vs->ProjectionY();
  KFchi2ndf_vs_py->SetLineColor(2);
  KFchi2ndf_vs_py->GetXaxis()->SetTitle("chi2/ndf");
  if(Spmode){
    KFchi2ndf_vs_px->Draw();
    KFchi2ndf_vs_py->Draw("same");
  }else{
    KFchi2ndf_vs_py->Draw();
    KFchi2ndf_vs_px->Draw("same");
  }

  TLegend *legKFchi2ndf_vs = new TLegend(0.55,0.65,0.76,0.82);
  legKFchi2ndf_vs->AddEntry(KFchi2ndf_vs_px,"#Sigma^{+} mode");
  legKFchi2ndf_vs->AddEntry(KFchi2ndf_vs_py,"#Sigma^{-} mode");
  legKFchi2ndf_vs->SetFillColor(0);
  legKFchi2ndf_vs->Draw();

  //int spbin = KFchi2ndf_vs_px->FindBin(chi2cut);
  //int smbin = KFchi2ndf_vs_py->FindBin(chi2cut);
  //std::cout << "Sp mode rejection ratio:" << KFchi2ndf_vs_px->Integral(0,spbin)/(KFchi2ndf_vs_px->Integral(0,201)) << std::endl;
  //std::cout << "Sm mode rejection ratio:" << KFchi2ndf_vs_py->Integral(0,smbin)/(KFchi2ndf_vs_py->Integral(0,201)) << std::endl;

  //cumulative dist. of prob.
  TCanvas *cKFchi2ndf_vs_cum = new TCanvas(Form("cKFchi2ndf_vs_cum"),"KFchi2ndf_vs_cum");
  cKFchi2ndf_vs_cum->cd();
  TH1 *KFchi2ndf_vs_cum_px = KFchi2ndf_vs_px->GetCumulative();
  KFchi2ndf_vs_cum_px->Scale(1./KFchi2ndf_vs_px->GetEntries());
  KFchi2ndf_vs_cum_px->SetXTitle("chi2ndf cut");
  KFchi2ndf_vs_cum_px->GetXaxis()->CenterTitle();
  KFchi2ndf_vs_cum_px->SetTitle("KF chi2ndf cumulative");
  TH1 *KFchi2ndf_vs_cum_py = KFchi2ndf_vs_py->GetCumulative();
  KFchi2ndf_vs_cum_py->SetLineColor(2);
  KFchi2ndf_vs_cum_py->Scale(1./KFchi2ndf_vs_py->GetEntries());
  KFchi2ndf_vs_cum_py->SetXTitle("chi2ndf cut");
  KFchi2ndf_vs_cum_py->SetTitle("KF chi2ndf cumulative");
  if(Spmode){
    KFchi2ndf_vs_cum_px->Draw();
    KFchi2ndf_vs_cum_py->Draw("same");
  }else{
    KFchi2ndf_vs_cum_py->Draw();
    KFchi2ndf_vs_cum_px->Draw("same");
  }



  TLegend *legKFchi2ndf_vs_cum = new TLegend(0.55,0.25,0.76,0.42);
  legKFchi2ndf_vs_cum->AddEntry(KFchi2ndf_vs_px,"#Sigma^{+} mode");
  legKFchi2ndf_vs_cum->AddEntry(KFchi2ndf_vs_py,"#Sigma^{-} mode");
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

  TCanvas *cMMnmiss_IMnpip_dE_woK0 = new TCanvas("cMMnmiss_IMnpip_dE_woK0","MMnmiss_IMnpip_dE_woK0");
  MMnmiss_IMnpip_dE_woK0->RebinX(2);
  MMnmiss_IMnpip_dE_woK0->RebinY(2);
  MMnmiss_IMnpip_dE_woK0->GetXaxis()->SetRangeUser(1.05,1.5);
  MMnmiss_IMnpip_dE_woK0->GetYaxis()->SetRangeUser(0.6,1.5);
  MMnmiss_IMnpip_dE_woK0->Draw("colz");

  TCanvas *cMMnmiss_IMnpim_dE_woK0 = new TCanvas("cMMnmiss_IMnpim_dE_woK0","MMnmiss_IMnpim_dE_woK0");
  MMnmiss_IMnpim_dE_woK0->RebinX(2);
  MMnmiss_IMnpim_dE_woK0->RebinY(2);
  MMnmiss_IMnpim_dE_woK0->GetXaxis()->SetRangeUser(1.05,1.5);
  MMnmiss_IMnpim_dE_woK0->GetYaxis()->SetRangeUser(0.6,1.5);
  MMnmiss_IMnpim_dE_woK0->Draw("colz");


  TCanvas *cq_IMnpipi_woK0_wSid_n_px = new TCanvas("cq_IMnpipi_woK0_wSid_n_px","q_IMnpipi_woK0_wSid_n_px"); 
  cq_IMnpipi_woK0_wSid_n_px->cd();
  TH1D *q_IMnpipi_wSid_n_px = q_IMnpipi_wSid_n->ProjectionX();
  TH1D *q_IMnpipi_woK0_wSid_n_px = q_IMnpipi_woK0_wSid_n->ProjectionX();
  TH1D *q_IMnpipi_woK0_kin_0_px = q_IMnpipi_woK0_kin[0]->ProjectionX();
  TH1D *q_IMnpipi_woK0_kin_1_px = q_IMnpipi_woK0_kin[1]->ProjectionX();
  //q_IMnpipi_woK0_wSid_n_px->Rebin(2);
  q_IMnpipi_woK0_wSid_n_px->Draw("HE");
  //q_IMnpipi_wSid_n_px->Rebin(2);
  q_IMnpipi_wSid_n_px->SetLineColor(5);
  q_IMnpipi_wSid_n_px->Draw("HEsame");
  q_IMnpipi_woK0_kin_0_px->SetLineColor(2);
  //q_IMnpipi_woK0_kin_1_px->Rebin(2);
  q_IMnpipi_woK0_kin_0_px->Draw("HEsame");
  //q_IMnpipi_woK0_kin_0_px->Rebin(2);
  q_IMnpipi_woK0_kin_1_px->SetLineColor(3);
  q_IMnpipi_woK0_kin_1_px->Draw("HEsame");
  TH1D *pxSum = (TH1D*) q_IMnpipi_woK0_kin_0_px->Clone();
  pxSum->Add(q_IMnpipi_woK0_kin_1_px);
  pxSum->SetLineColor(4);
  pxSum->Draw("HEsame");
  
  TCanvas *cq_IMnpipi_woK0_wSid_n_py = new TCanvas("cq_IMnpipi_woK0_wSid_n_py","q_IMnpipi_woK0_wSid_n_py"); 
  cq_IMnpipi_woK0_wSid_n_py->cd();
  TH1D *q_IMnpipi_wSid_n_py = q_IMnpipi_wSid_n->ProjectionY();
  TH1D *q_IMnpipi_woK0_wSid_n_py = q_IMnpipi_woK0_wSid_n->ProjectionY();
  TH1D *q_IMnpipi_woK0_kin_0_py = q_IMnpipi_woK0_kin[0]->ProjectionY();
  TH1D *q_IMnpipi_woK0_kin_1_py = q_IMnpipi_woK0_kin[1]->ProjectionY();
  //q_IMnpipi_woK0_wSid_n_py->Rebin(2);
  q_IMnpipi_woK0_wSid_n_py->Draw("HE");
  //q_IMnpipi_wSid_n_py->Rebin(2);
  q_IMnpipi_wSid_n_py->SetLineColor(5);
  q_IMnpipi_wSid_n_py->Draw("HEsame");
  q_IMnpipi_woK0_kin_1_py->SetLineColor(2);
  //q_IMnpipi_woK0_kin_1_py->Rebin(2);
  q_IMnpipi_woK0_kin_1_py->Draw("HEsame");
  //q_IMnpipi_woK0_kin_0_py->Rebin(2);
  q_IMnpipi_woK0_kin_0_py->SetLineColor(3);
  q_IMnpipi_woK0_kin_0_py->Draw("HEsame");
  TH1D *pySum = (TH1D*) q_IMnpipi_woK0_kin_0_py->Clone();
  pySum->Add(q_IMnpipi_woK0_kin_1_py);
  pySum->SetLineColor(4);
  pySum->Draw("HEsame");

  TCanvas *cq_IMnpipi_wSid_n = new TCanvas("cq_IMnpipi_wSid_n","q_IMnpipi_wSid_n");
  cq_IMnpipi_wSid_n->cd();
  q_IMnpipi_wSid_n->RebinX(2);
  q_IMnpipi_wSid_n->RebinY(6);
  q_IMnpipi_wSid_n->Draw("colz");

  TCanvas *cq_IMnpipi_woK0_wSid_n = new TCanvas("cq_IMnpipi_woK0_wSid_n","q_IMnpipi_woK0_wSid_n");
  cq_IMnpipi_woK0_wSid_n->cd();
  q_IMnpipi_woK0_wSid_n->RebinX(2);
  q_IMnpipi_woK0_wSid_n->RebinY(6);
  q_IMnpipi_woK0_wSid_n->Draw("colz");


  TCanvas *cdE_betainv_fid = new TCanvas(Form("cdE_betainv_fid"), "dE_betainv_fid");
  cdE_betainv_fid->cd();
  dE_betainv_fid->Draw("colz");

  
  TCanvas *cnmom_IMnpipi_woK0_wSid_n = new TCanvas(Form("cnmom_IMnpipi_woK0_wSid_n"),"nmom_IMnpipi_woK0_wSid_n");
  cnmom_IMnpipi_woK0_wSid_n->cd();
  nmom_IMnpipi_woK0_wSid_n->Draw("colz");
  
  TCanvas *cnmom_IMnpipi_woK0_wSid_n_py = new TCanvas(Form("cnmom_IMnpipi_woK0_wSid_n_py"),"nmom_IMnpipi_woK0_wSid_n_py");
  cnmom_IMnpipi_woK0_wSid_n_py->cd();
  TH1D* nmom_IMnpipi_woK0_wSid_n_py = nmom_IMnpipi_woK0_wSid_n->ProjectionY();
  nmom_IMnpipi_woK0_wSid_n_py->Draw();

  TCanvas *cMMnmiss_IMnpipi_woK0_wSid_n = new TCanvas(Form("cMMnmiss_IMnpipi_woK0_wSid_n"),"MMnmiss_IMnpipi_woK0_wSid_n");
  cMMnmiss_IMnpipi_woK0_wSid_n->cd();
  MMnmiss_IMnpipi_woK0_wSid_n->Draw("colz");

  TCanvas *cMMnmiss_IMnpipi_woK0_wSid_n_px = new TCanvas(Form("cMMnmiss_IMnpipi_woK0_wSid_n_px"),"MMnmiss_IMnpipi_woK0_wSid_n_px"); 
  cMMnmiss_IMnpipi_woK0_wSid_n_px->cd();
  MMnmiss_IMnpipi_woK0_wSid_n->ProjectionX()->Draw("");
  TCanvas *cMMnmiss_IMnpipi_woK0_wSid_n_py = new TCanvas(Form("cMMnmiss_IMnpipi_woK0_wSid_n_py"),"MMnmiss_IMnpipi_woK0_wSid_n_py"); 
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
    cdE_betainv_fid_kin[imode] = new TCanvas(Form("cdE_betainv_fid_kin_%s",smode[imode]),Form("dE_betainv_fid_kin_%s",smode[imode]));
    cdE_betainv_fid_kin[imode]->cd();
    dE_betainv_fid_kin[imode]->Draw("colz");


    cdE_MMom_fid_beta_woK0_kin[imode] = new TCanvas(Form("cdE_MMom_fid_beta_woK0_kin_%s",smode[imode]),Form("dE_MMom_fid_beta_woK0_kin_%s",smode[imode]));
    cdE_MMom_fid_beta_woK0_kin[imode]->cd();
    dE_MMom_fid_beta_woK0_kin[imode]->Draw("colz");

    cdE_MMass_fid_beta_woK0_kin[imode] = new TCanvas(Form("cdE_MMass_fid_beta_woK0_%s_kin",smode[imode]),Form("dE_MMass_fid_beta_woK0_%s_kin",smode[imode]));
    cdE_MMass_fid_beta_woK0_kin[imode]->cd(); 
    dE_MMass_fid_beta_woK0_kin[imode]->Draw("colz");

    cMMom_MMass_fid_beta_dE_woK0_kin[imode] = new TCanvas(Form("cMMom_MMass_fid_beta_dE_woK0_kin_%s",smode[imode]),Form("MMom_MMass_fid_beta_dE_woK0_kin_%s",smode[imode]));
    cMMom_MMass_fid_beta_dE_woK0_kin[imode]->cd();
    MMom_MMass_fid_beta_dE_woK0_kin[imode]->Draw("colz");

    cMMom_MMass_fid_beta_dE_woK0_kin_px[imode] = new TCanvas(Form("cMMom_MMass_fid_beta_dE_woK0_kin_px_%s",smode[imode]),Form("MMom_MMass_fid_beta_dE_woK0_kin_px_%s",smode[imode]));
    cMMom_MMass_fid_beta_dE_woK0_kin_px[imode]->cd();
    MMom_MMass_fid_beta_dE_woK0_kin_px[imode] = MMom_MMass_fid_beta_dE_woK0_kin[imode]->ProjectionX();
    MMom_MMass_fid_beta_dE_woK0_kin_px[imode]->Draw();

    cMMom_MMass_fid_beta_dE_woK0_kin_py[imode] = new TCanvas(Form("cMMom_MMass_fid_beta_dE_woK0_kin_py_%s",smode[imode]),Form("MMom_MMass_fid_beta_dE_woK0_kin_py_%s",smode[imode]));
    cMMom_MMass_fid_beta_dE_woK0_kin_py[imode]->cd();
    MMom_MMass_fid_beta_dE_woK0_kin_py[imode] = MMom_MMass_fid_beta_dE_woK0_kin[imode]->ProjectionY();
    MMom_MMass_fid_beta_dE_woK0_kin_py[imode]->Draw();
    

    cIMnpim_IMnpip_dE_woK0_kin[imode] = new TCanvas(Form("cIMnpim_IMnpip_dE_woK0_kin_%s",smode[imode]),Form("IMnpim_IMnpip_dE_woK0_kin_%s",smode[imode]));
    cIMnpim_IMnpip_dE_woK0_kin[imode]->cd();
    IMnpim_IMnpip_dE_woK0_kin[imode]->Draw("colz");
    
    cIMnpim_IMnpip_dE_woK0_kin_px[imode] = new TCanvas(Form("cIMnpim_IMnpip_dE_woK0_kin_px_%s",smode[imode]),Form("IMnpim_IMnpip_dE_woK0_kin_px_%s",smode[imode]));
    cIMnpim_IMnpip_dE_woK0_kin_px[imode]->cd();
    IMnpim_IMnpip_dE_woK0_kin_px[imode] = IMnpim_IMnpip_dE_woK0_kin[imode]->ProjectionX();
    IMnpim_IMnpip_dE_woK0_kin_px[imode]->Draw();
    
    cIMnpim_IMnpip_dE_woK0_kin_py[imode] = new TCanvas(Form("cIMnpim_IMnpip_dE_woK0_kin_py_%s",smode[imode]),Form("IMnpim_IMnpip_dE_woK0_kin_py_%s",smode[imode]));
    cIMnpim_IMnpip_dE_woK0_kin_py[imode]->cd();
    IMnpim_IMnpip_dE_woK0_kin_py[imode] = IMnpim_IMnpip_dE_woK0_kin[imode]->ProjectionY();
    IMnpim_IMnpip_dE_woK0_kin_py[imode]->Draw();
    
    cMMnpip_MMnpim_woK0_kin[imode] = new TCanvas(Form("cMMnpip_MMnpim_woK0_kin_%s",smode[imode]),Form("MMnpip_MMnpim_woK0_kin_%s",smode[imode]));
    cMMnpip_MMnpim_woK0_kin[imode]->cd();
    MMnpip_MMnpim_woK0_kin[imode]->Draw("colz");

    cdE_IMnpipi_woK0_kin[imode] = new TCanvas(Form("cdE_IMnpipi_woK0_kin_%s",smode[imode]),Form("dE_IMnpipi_woK0_kin_%s",smode[imode]));
    cdE_IMnpipi_woK0_kin[imode]->cd();
    dE_IMnpipi_woK0_kin[imode]->Draw("colz");
    
    cCosn_IMnpipi_n_kin[imode] = new TCanvas(Form("cCosn_IMnpipi_n_kin_%s",smode[imode]),Form("Cosn_IMnpipi_n_kin_%s",smode[imode]));
    cCosn_IMnpipi_n_kin[imode]->cd();
    Cosn_IMnpipi_woK0_kin[imode]->Draw("colz");

    cMMnmiss_IMnpipi_woK0_kin[imode] = new TCanvas(Form("cMMnmiss_IMnpipi_woK0_kin_%s",smode[imode]),Form("MMnmiss_IMnpipi_woK0_kin_%s",smode[imode]));
    cMMnmiss_IMnpipi_woK0_kin[imode]->cd();
    MMnmiss_IMnpipi_woK0_kin[imode]->Draw("colz");

    cq_IMnpipi_woK0_kin[imode] = new TCanvas(Form("cq_IMnpipi_woK0_kin_%s",smode[imode]),Form("q_IMnpipi_woK0_kin_%s",smode[imode]));
    cq_IMnpipi_woK0_kin[imode]->cd();
    q_IMnpipi_woK0_kin[imode]->RebinX(2);
    q_IMnpipi_woK0_kin[imode]->RebinY(6);
    q_IMnpipi_woK0_kin[imode]->Draw("colz");
    
    cnmom_IMnpipi_woK0_kin[imode] = new TCanvas(Form("cnmom_IMnpipi_woK0_kin_%s",smode[imode]),Form("nmom_IMnpipi_woK0_kin_%s",smode[imode]));
    cnmom_IMnpipi_woK0_kin[imode]->cd();
    nmom_IMnpipi_woK0_kin[imode]->RebinX(2);
    //nmom_IMnpipi_woK0_kin[imode]->RebinY(6);
    nmom_IMnpipi_woK0_kin[imode]->Draw("colz");
  }

  TCanvas *cMMnpip_MMnpim_woK0_wSid_n = new TCanvas(Form("cMMnpip_MMnpim_woK0_wSid_n"),Form("MMnpip_MMnpim_woK0_wSid_n"));
  cMMnpip_MMnpim_woK0_wSid_n->cd();
  MMnpip_MMnpim_woK0_wSid_n->Draw("colz");
  
  TCanvas *cMMom_MMass_fid_beta_dE_woK0  = new TCanvas(Form("cMMom_MMass_fid_beta_dE_woK0"),"MMom_MMass_fid_beta_dE_woK0");
  cMMom_MMass_fid_beta_dE_woK0->cd();
  MMom_MMass_fid_beta_dE_woK0->Draw("colz");
  
  TCanvas *cMMom_MMass_fid_beta_dE_woK0_px = new TCanvas(Form("cMMom_MMass_fid_beta_dE_woK0_px"),"MMom_MMass_fid_beta_dE_woK0_px");
  TH1D *MMom_MMass_fid_beta_dE_woK0_px = MMom_MMass_fid_beta_dE_woK0->ProjectionX();
  MMom_MMass_fid_beta_dE_woK0_px->Draw("");
  
  TCanvas *cMMom_MMass_fid_beta_dE_woK0_wSid  = new TCanvas(Form("cMMom_MMass_fid_beta_dE_woK0_wSid"),"MMom_MMass_fid_beta_dE_woK0_wSid");
  cMMom_MMass_fid_beta_dE_woK0_wSid->cd();
  MMom_MMass_fid_beta_dE_woK0_wSid->Draw("colz");
  

  TCanvas *cMMom_MMass_fid_beta_dE_woK0_wSid_px  = new TCanvas(Form("cMMom_MMass_fid_beta_dE_woK0_wSid_px"),"MMom_MMass_fid_beta_dE_woK0_wSid_px");
  cMMom_MMass_fid_beta_dE_woK0_wSid_px->cd();
  TH1D *MMom_MMass_fid_beta_dE_woK0_wSid_px = MMom_MMass_fid_beta_dE_woK0_wSid->ProjectionX();
  MMom_MMass_fid_beta_dE_woK0_wSid_px->Draw("");


  TCanvas *cMMom_MMass_fid_beta_dE_woK0_wSid_py  = new TCanvas(Form("cMMom_MMass_fid_beta_dE_woK0_wSid_py"),"MMom_MMass_fid_beta_dE_woK0_wSid_py");
  cMMom_MMass_fid_beta_dE_woK0_wSid_py->cd();
  TH1D *MMom_MMass_fid_beta_dE_woK0_wSid_py = MMom_MMass_fid_beta_dE_woK0_wSid->ProjectionY();
  MMom_MMass_fid_beta_dE_woK0_wSid_py->Draw("");
   

  TCanvas *cMMom_MMass_fid_beta_dE_woK0_px_sup = new TCanvas("cMMom_MMass_fid_beta_dE_woK0_px_sup","MMom_MMass_fid_beta_dE_woK0_px_sup");
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

  TCanvas *cMMom_MMass_fid_beta_dE_woK0_py = new TCanvas("cMMom_MMass_fid_beta_dE_woK0_py","MMom_MMass_fid_beta_dE_woK0_py");
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


  TCanvas  *cIMnpim_IMnpip_dE_woK0 = new TCanvas(Form("cIMnpim_IMnpip_dE_woK0"),"IMnpim_IMnpip_dE_woK0");
  cIMnpim_IMnpip_dE_woK0->cd();
  IMnpim_IMnpip_dE_woK0->Draw("colz");
  
  
  TCanvas *cIMnpim_IMnpip_dE_woK0_n = new TCanvas(Form("cIMnpim_IMnpip_dE_woK0_n"),"IMnpim_IMnpip_dE_woK0_n");
  cIMnpim_IMnpip_dE_woK0_n->cd();
  IMnpim_IMnpip_dE_woK0_n->Draw("colz");
  
  TCanvas *cIMnpim_IMnpip_dE_woK0_px = new TCanvas(Form("cIMnpim_IMnpip_dE_woK0_px"),"IMnpim_IMnpip_dE_woK0_px");
  cIMnpim_IMnpip_dE_woK0_px->cd();
  TH1D *IMnpim_IMnpip_dE_woK0_px = IMnpim_IMnpip_dE_woK0->ProjectionX();
  IMnpim_IMnpip_dE_woK0_px->Draw();
  //TH1D * IMnpim_IMnpip_dE_woK0_n_px_sum =(TH1D*) IMnpim_IMnpip_dE_woK0_kin_px[0]->Clone();
  
  TCanvas *cIMnpim_IMnpip_dE_woK0_py = new TCanvas(Form("cIMnpim_IMnpip_dE_woK0_py"),"IMnpim_IMnpip_dE_woK0_py");
  cIMnpim_IMnpip_dE_woK0_py->cd();
  TH1D *IMnpim_IMnpip_dE_woK0_py = IMnpim_IMnpip_dE_woK0->ProjectionY();
  IMnpim_IMnpip_dE_woK0_py->Draw();

  TCanvas *cIMnpim_IMnpip_dE_woK0_n_px = new TCanvas(Form("cIMnpim_IMnpip_dE_woK0_n_px"),"IMnpim_IMnpip_dE_woK0_n_px");
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
  
  TCanvas *cIMnpim_IMnpip_dE_woK0_n_py = new TCanvas(Form("cIMnpim_IMnpip_dE_woK0_n_py"),"IMnpim_IMnpip_dE_woK0_n_py");
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


  //TFile *fout = new TFile(outfilename.c_str(),"RECREATE");
  //fout->cd();

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
  TSeqCollection *SCol = gROOT->GetListOfCanvases();
  int size = SCol->GetSize();
  TIter next(SCol);
  for(int i=0;i<size;i++){
    c= (TCanvas*)next();
    c->Draw();
    c->cd();
    TPaveText *pt = new TPaveText(.80,0.90,0.98,0.99,"NDC");
    if(Spmode) pt->AddText("MC #Sigma+ mode");
    else pt->AddText("MC #Sigma- mode");
    pt->SetFillColor(kAzure-4);
    pt->SetBorderSize(1);
    pt->Draw();
    c->Modified();
    c->Update();
    if(i==0) c->Print(pdfname+"(",Form("Title:%s",c->GetTitle()));
    else if(i==size-1)c->Print(pdfname+")",Form("Title:%s",c->GetTitle())); 
    else c->Print(pdfname,Form("Title:%s",c->GetTitle())); 
  }
  
  //fout->Close();
  
}
