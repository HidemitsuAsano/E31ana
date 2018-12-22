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


#define CM 1

void plot_IMpisigma(const char* filename="")
{
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(111111);
  
  std::string outfilename = string(filename);
  outfilename.insert(outfilename.size()-5,"_post");
  std::cout << outfilename << std::endl;
  //--- color style ---//
  
  //= = = = pipipnn final-sample tree = = = =//
  TLorentzVector *mom_beam = 0;   // 4-momentum(beam)
  TLorentzVector *mom_target = 0; // 4-momentum(target)
  TLorentzVector *mom_pip = 0;    // 4-momentum(pi+)
  TLorentzVector *mom_pim = 0;    // 4-momentum(pi-)
  TLorentzVector *mom_n = 0;      // 4-momentum(neutron)
  double beta; // veracity of neutral particle on CDH
  double dE;   // energy deposit on CDH
  TVector3 *vtx_reaction = 0; // vertex(reaction)
  int run_num;   // run number
  int event_num; // event number
  int block_num; // block number
  TLorentzVector *kfSpmode_mom_beam = 0;   // 4-momentum(beam) after kinematical refit for pi- Sigma+
  TLorentzVector *kfSpmode_mom_pip = 0;    // 4-momentum(pi+) after kinematical refit for pi- Sigma+
  TLorentzVector *kfSpmode_mom_pim = 0;    // 4-momentum(pi-) after kinematical refit for pi- Sigma+
  TLorentzVector *kfSpmode_mom_n = 0;      // 4-momentum(neutron) after kinematical refit for pi- Sigma+
  double kfSpmode_chi2;   // chi2 of kinematical refit
  double kfSpmode_NDF;    // NDF of kinematical refit
  double kfSpmode_status; // status of kinematical refit -> details can be found in this code
  double kfSpmode_pvalue; // p-value of kinematical refit
  TLorentzVector *kfSmmode_mom_beam = 0;   // 4-momentum(beam) after kinematical refit for pi+ Sigma-
  TLorentzVector *kfSmmode_mom_pip = 0;    // 4-momentum(pi+) after kinematical refit for pi+ Sigma-
  TLorentzVector *kfSmmode_mom_pim = 0;    // 4-momentum(pi-) after kinematical refit for pi+ Sigma-
  TLorentzVector *kfSmmode_mom_n = 0;      // 4-momentum(neutron) after kinematical refit for pi+ Sigma-
  double kfSmmode_chi2;   // chi2 of kinematical refit
  double kfSmmode_NDF;    // NDF of kinematical refit
  double kfSmmode_status; // status of kinematical refit -> details can be found in this code
  double kfSmmode_pvalue; // p-value of kinematical refit
  int kf_flag; // flag of correct pair reconstruction, etc
  //= = = = pipipnn final-sample tree = = = =//
  
  TFile *f = new TFile(filename);
  //TFile *f = new TFile("sim_piSpn_dE0_Al.root");
  TTree *tree = (TTree*)f->Get("EventTree");

  tree->SetBranchAddress( "mom_beam",   &mom_beam );
  tree->SetBranchAddress( "mom_target", &mom_target );
  tree->SetBranchAddress( "mom_pip", &mom_pip );
  tree->SetBranchAddress( "mom_pim", &mom_pim );
  tree->SetBranchAddress( "mom_n", &mom_n );
  tree->SetBranchAddress( "beta", &beta );
  tree->SetBranchAddress( "dE", &dE );
  tree->SetBranchAddress( "vtx_reaction", &vtx_reaction );
  tree->SetBranchAddress( "run_num", &run_num );
  tree->SetBranchAddress( "event_num", &event_num );
  tree->SetBranchAddress( "block_num", &block_num );
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

  TLegend *leg;
  TLine *line;
  double ymax;
  TH1F *htmp1;
  char com[32];

  const double beta_MAX = 0.728786; // p = 1.0 GeV/c for neutron & 1/beta = 1.372
  //const double dE_MIN = 5.0; // 8.0MeVee * 3cm / 5cm;
  const double dE_MIN = 1.8;
 
  const double COSN_CUT = 0.2;
  
  const double pipi_MIN = 0.485;
  const double pipi_MAX = 0.510;
  const double ppi_MIN = 1.1075;
  const double ppi_MAX = 1.1225;

  const double neutron_MIN = 0.85;
  const double neutron_MAX = 1.03;

  const double Sigmap_MIN = 1.18;
  const double Sigmap_MAX = 1.20;
  const double Sigmam_MIN = 1.19;
  const double Sigmam_MAX = 1.21;
  //const double Sigma_width = 0.02;
  const double Sigma_width = 0.01;

  const double CHI2 = 10; // maximum value
  const double PVALUE[9] = {1E-60, 0.0001, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.5}; // minimum value
  
  TH2F* h2_dE_betainv_fid_beta_dE[2];
  TH2F* h2_dE_MMom_fid_beta_woK0[2];
  TH2F* h2_dE_MMass_fid_beta_woK0[2];
  TH2F* h2_MMom_MMass_fid_beta_dE_woK0[2];
  TH2F* h2_IMnpim_IMnpip_dE_woK0[2];
  TH2F* h2_IMnpim_IMnpip_dE_woK0_n[2];
  TH2F* h2_MMnpip_MMnpim_woK0_wSid_n[2];
  TH2F* h2_dE_IMnpipi_woK0_wSid_n[2];
  TH2F* h2_Cosn_IMnpipi_woK0_wSid_n[2];
  TH2F* h2_MMnmiss_IMnpipi_woK0_wSid_n[2];
  TH2F* h2_q_IMnpipi_woK0_wSid_n[2];

  for(int imode=0;imode<2;imode++){
    h2_dE_betainv_fid_beta_dE[imode] = new TH2F(Form("h2_dE_betainv_fid_beta_dE%d",imode),Form("dE_betainv_fid_beta_dE%d",imode),200, 0, 10, 200, 0, 50);
    h2_dE_betainv_fid_beta_dE[imode]->SetXTitle("1/#beta");
    h2_dE_betainv_fid_beta_dE[imode]->SetYTitle("dE [MeVee]");
    h2_dE_betainv_fid_beta_dE[imode]->GetXaxis()->CenterTitle();
    h2_dE_betainv_fid_beta_dE[imode]->GetYaxis()->CenterTitle();

  
    h2_dE_MMom_fid_beta_woK0[imode] = new TH2F(Form("h2_dE_MMom_fid_beta_woK0%d",imode),Form("dE_MMom_fid_beta_woK0%d",imode),100, 0, 1.5, 200, 0, 50);
    h2_dE_MMom_fid_beta_woK0[imode]->SetXTitle("Missing Mom. [GeV/c]");
    h2_dE_MMom_fid_beta_woK0[imode]->SetYTitle("dE [MeVee]");
    h2_dE_MMom_fid_beta_woK0[imode]->GetXaxis()->CenterTitle();
    h2_dE_MMom_fid_beta_woK0[imode]->GetYaxis()->CenterTitle();

    h2_dE_MMass_fid_beta_woK0[imode] = new TH2F(Form("h2_dE_MMass_fid_beta_woK0%d",imode),Form("dE_MMass_fid_beta_woK0%d",imode), 140, 0.4, 1.8, 200, 0, 50);
    h2_dE_MMass_fid_beta_woK0[imode]->SetXTitle("Missing mass [GeV/c^{2}]");
    h2_dE_MMass_fid_beta_woK0[imode]->SetYTitle("dE [MeVee]");
    h2_dE_MMass_fid_beta_woK0[imode]->GetXaxis()->CenterTitle();
    h2_dE_MMass_fid_beta_woK0[imode]->GetYaxis()->CenterTitle();
  

    h2_MMom_MMass_fid_beta_dE_woK0[imode] = new TH2F(Form("h2_MMom_MMass_fid_beta_dE_woK0%d",imode),Form("MMom_MMass_fid_beta_dE_woK0%d",imode), 140, 0.4, 1.8, 100, 0, 1.5);
    h2_MMom_MMass_fid_beta_dE_woK0[imode]->SetXTitle("Missing Mass [GeV/c^{2}]");
    h2_MMom_MMass_fid_beta_dE_woK0[imode]->SetYTitle("Missing Mom. [GeV/c]");
    h2_MMom_MMass_fid_beta_dE_woK0[imode]->GetXaxis()->CenterTitle();
    h2_MMom_MMass_fid_beta_dE_woK0[imode]->GetYaxis()->CenterTitle();

    h2_IMnpim_IMnpip_dE_woK0[imode] = new TH2F(Form("h2_IMnpim_IMnpip_dE_woK0%d",imode), Form("IMnpim_IMnpip_dE_woK0%d",imode),140, 1, 1.7, 140, 1, 1.7);
    h2_IMnpim_IMnpip_dE_woK0[imode]->SetXTitle("IM(n#pi^{+} [GeV/c^{2}]");
    h2_IMnpim_IMnpip_dE_woK0[imode]->SetYTitle("IM(n#pi^{-} [GeV/c^{2}]");
    h2_IMnpim_IMnpip_dE_woK0[imode]->GetXaxis()->CenterTitle();
    h2_IMnpim_IMnpip_dE_woK0[imode]->GetXaxis()->CenterTitle();
  

    h2_IMnpim_IMnpip_dE_woK0_n[imode] = new TH2F(Form("h2_IMnpim_IMnpip_dE_woK0_n%d",imode),Form("IMnpim_IMnpip_dE_woK0_n%d",imode),140, 1, 1.7, 140, 1, 1.7);
    h2_IMnpim_IMnpip_dE_woK0_n[imode]->SetXTitle("IM(n#pi^{+} [GeV/c^{2}]");
    h2_IMnpim_IMnpip_dE_woK0_n[imode]->SetYTitle("IM(n#pi^{-} [GeV/c^{2}]");
    h2_IMnpim_IMnpip_dE_woK0_n[imode]->GetXaxis()->CenterTitle();
    h2_IMnpim_IMnpip_dE_woK0_n[imode]->GetYaxis()->CenterTitle();

  
    h2_MMnpip_MMnpim_woK0_wSid_n[imode] = new TH2F(Form("h2_MMnpip_MMnpim_woK0_wSid_n%d",imode),Form("MMnpip_MMnpim_woK0_wSid_n%d",imode),70, 1, 1.7, 70, 1, 1.7);
    h2_MMnpip_MMnpim_woK0_wSid_n[imode]->SetXTitle("Miss. Mass{n#pi^{+}} [GeV/c^2]");
    h2_MMnpip_MMnpim_woK0_wSid_n[imode]->SetYTitle("Miss. Mass{n#pi^{-}} [GeV/c^2]");
    h2_MMnpip_MMnpim_woK0_wSid_n[imode]->GetXaxis()->CenterTitle();
    h2_MMnpip_MMnpim_woK0_wSid_n[imode]->GetYaxis()->CenterTitle();

 
    h2_dE_IMnpipi_woK0_wSid_n[imode] = new TH2F(Form("h2_dE_IMnpipi_woK0_wSid_n%d",imode),Form("dE_IMnpipi_woK0_wSid_n%d",imode),100, 1, 2, 200, 0, 50);
    h2_dE_IMnpipi_woK0_wSid_n[imode]->SetXTitle("IM(n#pi^{+}#pi^{-} [GeV/c^{2}]");
    h2_dE_IMnpipi_woK0_wSid_n[imode]->SetYTitle("dE [MeVee]");
    h2_dE_IMnpipi_woK0_wSid_n[imode]->GetXaxis()->CenterTitle();
    h2_dE_IMnpipi_woK0_wSid_n[imode]->GetYaxis()->CenterTitle();

  
    h2_Cosn_IMnpipi_woK0_wSid_n[imode] = new TH2F(Form("h2_Cosn_IMnpipi_woK0_wSid_n%d",imode),Form("dE_Cosn_IMnpipi_woK0_wSid_n%d",imode),100, 1, 2, 50, -1, 1);
    h2_Cosn_IMnpipi_woK0_wSid_n[imode]->SetXTitle("IM(n#pi^{+}#pi^{-} [GeV/c^{2}]");
    h2_Cosn_IMnpipi_woK0_wSid_n[imode]->SetYTitle("cos#theta_{n}");
    h2_Cosn_IMnpipi_woK0_wSid_n[imode]->GetXaxis()->CenterTitle();
    h2_Cosn_IMnpipi_woK0_wSid_n[imode]->GetYaxis()->CenterTitle();
  
    h2_MMnmiss_IMnpipi_woK0_wSid_n[imode] = new TH2F(Form("h2_MMnmiss_IMnpipi_woK0_wSid_n%d",imode),Form("MMnmiss_IMnpipi_woK0_wSid_n%d",imode),100,1,2,100,0,1.5);
    h2_MMnmiss_IMnpipi_woK0_wSid_n[imode]->SetXTitle("IM(n#pi^{+}#pi^{-} [GeV/c^{2}]");
    h2_MMnmiss_IMnpipi_woK0_wSid_n[imode]->SetYTitle("Miss Mom. [GeV/c]");
    h2_MMnmiss_IMnpipi_woK0_wSid_n[imode]->GetXaxis()->CenterTitle();
    h2_MMnmiss_IMnpipi_woK0_wSid_n[imode]->GetYaxis()->CenterTitle();

    h2_q_IMnpipi_woK0_wSid_n[imode] = new TH2F(Form("h2_q_IMnpipi_woK0_wSid_n%d",imode),Form("q_IMnpipi_woK0_wSid_n%d",imode),100,1,2,300,0,1.5);
    h2_q_IMnpipi_woK0_wSid_n[imode]->SetXTitle("IM(n#pi^{+}#pi^{-} [GeV/c^{2}]");
    h2_q_IMnpipi_woK0_wSid_n[imode]->SetYTitle("Mom. Transfer [GeV/c]");
    h2_q_IMnpipi_woK0_wSid_n[imode]->GetXaxis()->CenterTitle();
    h2_q_IMnpipi_woK0_wSid_n[imode]->GetYaxis()->CenterTitle();
  };

  TH2F *KFpvalue_vs = new TH2F("KFpvalue_vs", "KFpvalue_vs", 100, 0, 1, 100, 0, 1 );
  

  TH1F *pmom = new TH1F("pmom", "pmom", 50, 0, 1.0);
  TH1F *nmom = new TH1F("nmom", "nmom", 50, 0, 1.0);
  TH1F *mnmom = new TH1F("mnmom", "mnmom", 50, 0, 1.0);
  TH1F *npipmom = new TH1F("npipmom", "npipmom", 50, 0, 1.0);
  TH1F *npimmom = new TH1F("npimmom", "npimmom", 50, 0, 1.0);

  Int_t nevent = tree->GetEntries();
  std::cerr<<"# of events = "<<nevent<<std::endl;
  //------------------------//
  //--- event roop start ---//
  //------------------------//
  for ( Int_t i=0; i<nevent; i++ ) {
    tree->GetEvent(i);
    //std::cerr<<i<<" "<<run_num<<" "<<event_num<<" "<<block_num<<std::endl;

    //** calc missing n **//
    TLorentzVector mom_n_miss = *mom_target+*mom_beam-*mom_pip-*mom_pim-*mom_n;
    double miss_mass = mom_n_miss.M();
    
    // calc cos(theta) of missing n //
    TVector3 boost = (*mom_target+*mom_beam).BoostVector();
    TLorentzVector mom_n_miss_CM = mom_n_miss;
    TLorentzVector mom_beam_CM = *mom_beam;
    mom_n_miss_CM.Boost(-boost);
    mom_beam_CM.Boost(-boost);
    double cos_n = mom_n_miss_CM.Vect().Dot(mom_beam_CM.Vect())/(mom_n_miss_CM.Vect().Mag()*mom_beam_CM.Vect().Mag());
    TLorentzVector qkn = *mom_beam-mom_n_miss;
   
    //** calc pi+pi- **//
    TLorentzVector mom_pip_pim = *mom_pip+*mom_pim;

    //** calc pi+n **//
    TLorentzVector mom_pip_n = *mom_pip+*mom_n;

    //** calc pi-n **//
    TLorentzVector mom_pim_n = *mom_pim+*mom_n;

    //** calc missing Sp **//
    TLorentzVector mom_pip_n_miss = *mom_target+*mom_beam-*mom_pim-*mom_n;

    //** calc missing Sm **//
    TLorentzVector mom_pim_n_miss = *mom_target+*mom_beam-*mom_pip-*mom_n;


    //** calc pi+pi-n **//
    TLorentzVector mom_pip_pim_n = *mom_pip+*mom_pim+*mom_n;
    TLorentzVector mom_pip_pim_n_CM = mom_pip_pim_n;
    mom_pip_pim_n_CM.Boost(-boost);
    double cos_X = mom_pip_pim_n_CM.Vect().Dot(mom_beam_CM.Vect())/(mom_pip_pim_n_CM.Vect().Mag()*mom_beam_CM.Vect().Mag());



    double chi2 = kfSpmode_chi2<kfSmmode_chi2 ? kfSpmode_chi2:kfSmmode_chi2;
    double pvalue = kfSmmode_pvalue<kfSpmode_pvalue ? kfSpmode_pvalue:kfSmmode_pvalue;

    bool K0rejectFlag=false;
    bool MissNFlag=false;
    bool SigmaPFlag=false;
    bool SigmaMFlag=false;
    bool NBetaOK=false;
    bool NdEOK=false;
    //-- neutron-ID, K0 and missing neutron selection --//
    //-- original

    if( beta<beta_MAX && dE_MIN<dE &&
    (mom_pip_pim.M()<pipi_MIN || pipi_MAX<mom_pip_pim.M())
    ){
   
      if( -1<kf_flag && 0.05<pvalue ){
        if( kfSmmode_pvalue<kfSpmode_pvalue ){
          h2_IMnpim_IMnpip_Sp->Fill(mom_pip_n.M(),mom_pim_n.M() );
          h2_MMom_MMass_fid_beta_dE_Sp->Fill(mom_n_miss.M(), mom_n_miss.P());
        }else{
          h2_IMnpim_IMnpip_Sm->Fill(mom_pip_n.M(),mom_pim_n.M() );
          h2_MMom_MMass_fid_beta_dE_Sm->Fill(mom_n_miss.M(), mom_n_miss.P());
        }
      }
    }
	}
  
  TFile *fout = new TFile(outfilename.c_str(),"RECREATE");
  fout->cd();
  h2_IMnpim_IMnpip_Sp->Write();
  h2_IMnpim_IMnpip_Sm->Write();
  h2_MMom_MMass_fid_beta_dE_Sp->Write();
  h2_MMom_MMass_fid_beta_dE_Sm->Write();

  fout->Close();
  
}
