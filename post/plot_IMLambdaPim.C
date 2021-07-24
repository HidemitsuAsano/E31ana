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

#include "../src/GlobalVariables.h"
#include "anacuts.h"

const double pvalcut = 0.005;
//const double pvalcut = 1.0e-5;
const bool gridon=true;
const bool staton=false;
const bool UseKinFitVal = false;

//mode 0: Sigma+ ,1: Sigma- 
void plot_IMLambdaPim(const char* filename="", const int qvalcutflag=0)
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
  std::cout << "Use Kin Fit Val ? " << std::endl;

  bool SimMode = (std::string(filename).find("sim")!= std::string::npos);

  if(UseKinFitVal) std::cout << "Yes" << std::endl;
  else             std::cout << "No"  << std::endl;

  TH1::SetDefaultSumw2();
  //--- color style ---//
  //= = = = pipipnn final-sample tree = = = =//
  TLorentzVector *LVec_beam=nullptr;   // 4-momentum(beam)
  TLorentzVector *LVec_target=nullptr; // 4-momentum(target)
  TLorentzVector *LVec_pim1=nullptr;    // 4-momentum(pi+)
  TLorentzVector *LVec_pim2=nullptr;    // 4-momentum(pi-)
  TLorentzVector *LVec_p=nullptr;      // 4-momentum(proton)
  TLorentzVector *mcmom_beam=nullptr;   // generated 4-momentum(beam)
  TLorentzVector *mcmom_pim1=nullptr;    // generated 4-momentum(pi+)
  TLorentzVector *mcmom_pim2=nullptr;    // generated 4-momentum(pi-)
  TLorentzVector *mcmom_p=nullptr;      // generated 4-momentum(neutron)
  TLorentzVector *mcmom_pmiss=nullptr;      // generated 4-momentum(neutron)
  TLorentzVector *react_pmiss=nullptr;
  TLorentzVector *react_Lambda=nullptr;
  TLorentzVector *react_pim=nullptr;
  TVector3 *vtx_reaction = nullptr; // vertex(reaction) 
  TVector3 *vtx_displaced = nullptr; // vertex(reaction) 
  TVector3 *vtx_pim1_beam = nullptr; //C.A.P of pip-beam beam side
  TVector3 *vtx_pim2_beam = nullptr; //C.A.P of pim-beam beam side
  TVector3 *vtx_p_beam = nullptr; //C.A.P of pim-beam beam side
  TVector3 *vtx_pim1_cdc = nullptr;//C.A.P of pip-beam pip side
  TVector3 *vtx_pim2_cdc = nullptr;//C.A.P of pim-beam pim side
  TVector3 *vtx_p_cdc = nullptr;//C.A.P of pim-beam pim side
  TVector3 *CA_pim1_pim1p = nullptr;//C.A.P of pip-pim pip side
  TVector3 *CA_p_pim1p = nullptr;//C.A.P of pip-pim pim side
  TVector3 *CA_pim2_pim2p = nullptr;//C.A.P of pip-pim pip side
  TVector3 *CA_p_pim2p = nullptr;//C.A.P of pip-pim pim side
  TVector3 *CA_pim1_pim1pim2 = nullptr;//C.A.P of pip-pim pim side
  TVector3 *CA_pim2_pim1pim2 = nullptr;//C.A.P of pip-pim pim side
  //int run_num;   // run number
  //int event_num; // event number
  //int block_num; // block number
  TLorentzVector *kf_mom_beam=nullptr;   // 4-momentum(beam) after kinematical refit for pi- Sigma+
  TLorentzVector *kf_mom_pim1=nullptr;    // 4-momentum(pi+) after kinematical refit for pi- Sigma+
  TLorentzVector *kf_mom_pim2=nullptr;    // 4-momentum(pi-) after kinematical refit for pi- Sigma+
  TLorentzVector *kf_mom_proton=nullptr;      // 4-momentum(neutron) after kinematical refit for pi- Sigma+
  double kf_chi2;   // chi2 of kinematical refit
  double kf_NDF;    // NDF of kinematical refit
  double kf_status; // status of kinematical refit -> details can be found in this code
  double kf_pvalue; // p-value of kinematical refit
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
  tree->SetBranchAddress( "mom_pim1", &LVec_pim1 );
  tree->SetBranchAddress( "mom_pim2", &LVec_pim2 );
  tree->SetBranchAddress( "mom_p", &LVec_p );
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

  //tree->SetBranchAddress( "run_num", &run_num );
  //tree->SetBranchAddress( "event_num", &event_num );
  //tree->SetBranchAddress( "block_num", &block_num );
  //if(Spmode ){
  //  tree->SetBranchAddress( "mcmom_beam",  &mcmom_beam );
  //  tree->SetBranchAddress( "mcmom_pip", &mcmom_pip);
  //  tree->SetBranchAddress( "mcmom_pim", &mcmom_pim);
  //  tree->SetBranchAddress( "mcmom_ncds", &mcmom_ncds);
  //  tree->SetBranchAddress( "mcmom_nmiss", &mcmom_nmiss);
  //}

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
  TH2F* MMom_MMass_p;
  TH2F* MMom_MMass_p_wL;
  TH2F* MMom_PMom_p;
  TH2F* IMppim1_IMppim2;
  TH2F* IMppim1_IMppim2_p;
  TH2F* IMppim1_IMppim2_p_wL;
  TH2F* MMass_IMppim1;
  TH2F* MMass_IMppim2;
  TH2F* MMass_IMppim_wL;
  TH2F* MMass_IMppim_p_wL;
  TH2F* q_IMppipi_p;
  TH2F* q_IMppipi_p_wL;
  TH1F* DCA_pim1_beam = new TH1F("DCA_pim1_beam","DCA_pim1_beam",300,0,30);
  DCA_pim1_beam->SetXTitle("DCA #pi^{-}1 [cm]");
  TH1F* DCA_pim2_beam = new TH1F("DCA_pim2_beam","DCA_pim2_beam",300,0,30);
  DCA_pim2_beam->SetXTitle("DCA #pi^{-}2 [cm]");
  TH1F* DCA_pim1_pim2 = new TH1F("DCA_pim1_pim2","DCA_pim1_pim2",300,0,30);
  DCA_pim1_pim2->SetXTitle("DCA #pi^{-}1#pi^{-}2");

  // w/ kinematic fit
  //TH2F* MMom_MMass_fid_kin;
  //TH2F* IMppim1_IMppim2_kin;
  //TH2F* MMom_PMom_kin[];
  //TH2F* q_IMppipi_kin[2];
  
  const int nbinIMppipi = 160;//1.2-2 GeV/c^2
  const double IMppipilow = 1.2;//1.2-2 GeV/c^2
  const double IMppipihigh = 2.0;//1.2-2 GeV/c^2
  const int nbinq = 100;//0-1.5 GeV/c
  const int nbinIMppi = 500; //1-2 GeV/c^2
  const int nbinpmiss = 150; //0-1.5 GeV/c
  const double pmisslow = 0.0;
  const double pmisshigh = 1.5;
  
  MMom_MMass = new TH2F("MMom_MMass","MMom_MMass", nbinpmiss, pmisslow, pmisshigh, 200, 0, 2.0);
  MMom_MMass->SetXTitle("Missing Mass [GeV/c^{2}]");
  MMom_MMass->SetYTitle("Missing Mom. [GeV/c]");

  MMom_MMass_p = new TH2F("MMom_MMass_p","MMom_MMass_p", nbinpmiss, pmisslow, pmisshigh, 200, 0, 2.0);
  MMom_MMass_p->SetXTitle("Missing Mass [GeV/c^{2}]");
  MMom_MMass_p->SetYTitle("Missing Mom. [GeV/c]");
   
  MMom_MMass_p_wL = new TH2F("MMom_MMass_p_wL","MMom_MMass_p_wL", nbinpmiss, pmisslow, pmisshigh, 200, 0, 2.0);
  MMom_MMass_p_wL->SetXTitle("Missing Mass [GeV/c^{2}]");
  MMom_MMass_p_wL->SetYTitle("Missing Mom. [GeV/c]");
   
  IMppim1_IMppim2 = new TH2F("IMppim1_IMppim2","IMppim1_IMppim2",nbinIMppi,1.,2.0,nbinIMppi,1.,2.0);
  IMppim1_IMppim2->SetXTitle("IM(p#pi^{-}2) [GeV/c^{2}]");
  IMppim1_IMppim2->SetYTitle("IM(p#pi^{-}1) [GeV/c^{2}]");

  IMppim1_IMppim2_p = new TH2F("IMppim1_IMppim2_p","IMppim1_IMppim2_p",nbinIMppi,1.,2.0,nbinIMppi,1.,2.0);
  IMppim1_IMppim2_p->SetXTitle("IM(p#pi^{-}2) [GeV/c^{2}]");
  IMppim1_IMppim2_p->SetYTitle("IM(p#pi^{-}1) [GeV/c^{2}]");

  IMppim1_IMppim2_p_wL = new TH2F("IMppim1_IMppim2_p_wL","IMppim1_IMppim2_p_wL",nbinIMppi,1.,2.0,nbinIMppi,1.,2.0);
  IMppim1_IMppim2_p_wL->SetXTitle("IM(p#pi^{-}2) [GeV/c^{2}]");
  IMppim1_IMppim2_p_wL->SetYTitle("IM(p#pi^{-}1) [GeV/c^{2}]");

  MMass_IMppim1 = new TH2F("MMass_IMppim1","MMass_IMppim1",nbinIMppi,1.,2.0,nbinpmiss, pmisslow, pmisshigh);
  MMass_IMppim1->SetXTitle("IM(p#pi^{-} [GeV/c^{2}]");
  MMass_IMppim1->SetYTitle("Missing Mass [GeV/c^{2}]");
  
  MMass_IMppim2 = new TH2F("MMass_IMppim2","MMass_IMppim2",nbinIMppi,1.,2.0,nbinpmiss, pmisslow, pmisshigh);
  MMass_IMppim2->SetXTitle("IM(p#pi^{-} [GeV/c^{2}]");
  MMass_IMppim2->SetYTitle("Missing Mass [GeV/c^{2}]");
  
  MMass_IMppim_wL = new TH2F("MMass_IMppim_wL","MMass_IMppim_wL",nbinIMppi,1.,2.0,nbinpmiss, pmisslow, pmisshigh);
  MMass_IMppim_wL->SetXTitle("IM(p#pi^{-} [GeV/c^{2}]");
  MMass_IMppim_wL->SetYTitle("Missing Mass [GeV/c^{2}]");

  MMass_IMppim_p_wL = new TH2F("MMass_IMppim_p_wL","MMass_IMppim_p_wL",nbinIMppi,1.,2.0,nbinpmiss, pmisslow, pmisshigh);
  MMass_IMppim_p_wL->SetXTitle("IM(p#pi^{-} [GeV/c^{2}]");
  MMass_IMppim_p_wL->SetYTitle("Missing Mass [GeV/c^{2}]");
  
  q_IMppipi_p = new TH2F("q_IMppipi_p","q_IMppipi_p",nbinIMppipi,IMppipilow,IMppipihigh, nbinq,0,1.5);
  q_IMppipi_p->SetXTitle("IM(p#pi^{-}#pi^{-}) [GeV/c^{2}]");
  q_IMppipi_p->SetYTitle("Mom. Transfer [GeV/c]");

  q_IMppipi_p_wL = new TH2F("q_IMppipi_p_wL","q_IMppipi_p_wL",nbinIMppipi,IMppipilow,IMppipihigh, nbinq,0,1.5);
  q_IMppipi_p_wL->SetXTitle("IM(p#pi^{-}#pi^{-}) [GeV/c^{2}]");
  q_IMppipi_p_wL->SetYTitle("Mom. Transfer [GeV/c]");

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
  for ( Int_t i=0; i<nevent; i++ ) {
    tree->GetEvent(i);
    if(i%50000==0) std::cout << "Event# " << i << std::endl; 
    TVector3 vtx_pim1 = *vtx_pim1_cdc ;
    TVector3 vtx_pim2 = *vtx_pim2_cdc ;
    // calc missing p //
    TLorentzVector LVec_p_miss = *LVec_target+*LVec_beam-*LVec_pim1-*LVec_pim2-*LVec_p;
    double pmiss_mass = LVec_p_miss.M();
    double pmiss_mom = LVec_p_miss.P();

    // calc cos(theta) of missing p //
    TVector3 boost = (*LVec_target+*LVec_beam).BoostVector();
    TLorentzVector LVec_p_miss_CM = LVec_p_miss;
    TLorentzVector LVec_beam_CM = *LVec_beam;
    LVec_p_miss_CM.Boost(-boost);
    LVec_beam_CM.Boost(-boost);
    double cos_p = LVec_p_miss_CM.Vect().Dot(LVec_beam_CM.Vect())/(LVec_p_miss_CM.Vect().Mag()*LVec_beam_CM.Vect().Mag());
    TLorentzVector qkn = *LVec_beam-LVec_p_miss;
    
    // calc pi-pi- //
    TLorentzVector LVec_pim1_pim2 = *LVec_pim1+*LVec_pim2;
     

    // calc pi+p //
    TLorentzVector LVec_pim1_p = *LVec_pim1+*LVec_p;
    TLorentzVector LVec_pim2_p = *LVec_pim2+*LVec_p;

    // calc missing Lambda //
    TLorentzVector LVec_pim1_p_miss = *LVec_target+*LVec_beam-*LVec_pim1-*LVec_p;
    TLorentzVector LVec_pim2_p_miss = *LVec_target+*LVec_beam-*LVec_pim2-*LVec_p;
    
    // calc pi+pi-n //
    TLorentzVector LVec_pim1_pim2_p = *LVec_pim1+*LVec_pim2+*LVec_p;
    TLorentzVector LVec_pim1_pim2_p_CM = LVec_pim1_pim2_p;
    LVec_pim1_pim2_p_CM.Boost(-boost);
    double cos_X = LVec_pim1_pim2_p_CM.Vect().Dot(LVec_beam_CM.Vect())/(LVec_pim1_pim2_p_CM.Vect().Mag()*LVec_beam_CM.Vect().Mag());

    if( (qkn.P()>=anacuts::qvalcut) && (qvalcutflag==1) ) continue;
    if( (qkn.P()<anacuts::qvalcut) && (qvalcutflag==2) ) continue;
    
    //double chi2 = kfSpmode_chi2<kfSmmode_chi2 ? kfSpmode_chi2:kfSmmode_chi2;
    //double pvalue = kfSmmode_pvalue<kfSpmode_pvalue ? kfSpmode_pvalue:kfSmmode_pvalue;

    bool MissPFlag=false;
    bool LambdaFlag=false;
    //-- neutron-ID, K0 and missing neutron selection --//

    double dca_pim1_beam = (*vtx_pim1_beam-*vtx_pim2_cdc).Mag();
    double dca_pim2_beam = (*vtx_pim2_beam-*vtx_pim2_cdc).Mag();
    double dca_pim1_pim2 =(*CA_pim1_pim1pim2-*CA_pim2_pim1pim2).Mag();
   
    //Lambda production in CDS
    if( (anacuts::Lambda_MIN<LVec_pim1_p.M() && LVec_pim1_p.M()<anacuts::Lambda_MAX)) LambdaFlag=true;
    if( (anacuts::Lambda_MIN<LVec_pim2_p.M() && LVec_pim2_p.M()<anacuts::Lambda_MAX)) LambdaFlag=true;
        
    if(anacuts::Proton_MIN<pmiss_mass && pmiss_mass<anacuts::Proton_MAX ) MissPFlag=true;
    
    MMom_MMass->Fill(pmiss_mass,pmiss_mom);
    IMppim1_IMppim2->Fill(LVec_pim2_p.M(),LVec_pim1_p.M());
    MMass_IMppim1->Fill(LVec_pim1_p.M(),pmiss_mass);
    MMass_IMppim2->Fill(LVec_pim2_p.M(),pmiss_mass);
    if(MissPFlag){
      MMom_MMass_p->Fill(pmiss_mass,pmiss_mom);
      IMppim1_IMppim2_p->Fill(LVec_pim2_p.M(),LVec_pim1_p.M());
      q_IMppipi_p->Fill(LVec_pim1_pim2_p.M(),qkn.P());
    }
    if(LambdaFlag){
      if((anacuts::Lambda_MIN<LVec_pim1_p.M() && LVec_pim1_p.M()<anacuts::Lambda_MAX)){
        MMass_IMppim_wL->Fill(LVec_pim1_p.M(),pmiss_mass);
      }else if((anacuts::Lambda_MIN<LVec_pim2_p.M() && LVec_pim2_p.M()<anacuts::Lambda_MAX)){
        MMass_IMppim_wL->Fill(LVec_pim2_p.M(),pmiss_mass);
      }
    }
    if(MissPFlag && LambdaFlag){
      if((anacuts::Lambda_MIN<LVec_pim1_p.M() && LVec_pim1_p.M()<anacuts::Lambda_MAX)){
        MMass_IMppim_p_wL->Fill(LVec_pim1_p.M(),pmiss_mass);
      }else if((anacuts::Lambda_MIN<LVec_pim2_p.M() && LVec_pim2_p.M()<anacuts::Lambda_MAX)){
        MMass_IMppim_p_wL->Fill(LVec_pim2_p.M(),pmiss_mass);
      }
      MMom_MMass_p_wL->Fill(pmiss_mass,pmiss_mom);
      IMppim1_IMppim2_p_wL->Fill(LVec_pim2_p.M(),LVec_pim1_p.M());
      q_IMppipi_p_wL->Fill(LVec_pim1_pim2_p.M(),qkn.P());
      DCA_pim1_beam->Fill( dca_pim1_beam );
      DCA_pim2_beam->Fill( dca_pim2_beam );
      DCA_pim1_pim2->Fill( dca_pim1_pim2 );
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
  TCanvas *ctest = new TCanvas("ctest","ctest");
  q_IMppipi_p_wL->Draw("colz");
  
  f->cd();
  TIter nexthist(gDirectory->GetList());
  TH1F *h1 = nullptr;
  TH1D *h1d = nullptr;
  TH2F *h2 = nullptr;
  TObject *obj = nullptr;
  TString outname = std::string(filename);
  outname.Replace(std::string(filename).size()-5,5,"_out.root");
  if(qvalcutflag==1) outname.Replace(std::string(outname).size()-5,5,"_qlo.root");
  if(qvalcutflag==2) outname.Replace(std::string(outname).size()-5,5,"_qhi.root");
  if(qvalcutflag==3) outname.Replace(std::string(outname).size()-5,5,"_theta15.root");
  TFile *fout = new TFile(outname.Data(),"RECREATE");
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
    obj->Write();
  }
  fout->Close();
  
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
  while( (obj = (TObject*)nexthist2())!=nullptr) {
    obj->Print();
    obj->Write();
  }
  fout->cd();
  fout->Close();
  */
}
