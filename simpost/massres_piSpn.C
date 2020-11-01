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

#define DEBUG 0
#define CM 1

int main()
{
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(111111);

  //--- color style ---//
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
  //--- color style ---//

  TDatabasePDG *pdg = new TDatabasePDG();
  pdg->ReadPDGTable("./pdg_table.txt");
  double n_mass   = pdg->GetParticle("neutron")->Mass();
  double L_mass   = pdg->GetParticle("Lambda0")->Mass();
  double Sp_mass  = pdg->GetParticle("Sigma+")->Mass();
  double Sm_mass  = pdg->GetParticle("Sigma-")->Mass();
  double Kpp_mass   = 2*pdg->GetParticle("proton")->Mass()+pdg->GetParticle("K-")->Mass();
  double piSp_mass1 = pdg->GetParticle("pi+")->Mass()+pdg->GetParticle("Sigma-")->Mass()+pdg->GetParticle("proton")->Mass();
  double piSp_mass2 = pdg->GetParticle("pi-")->Mass()+pdg->GetParticle("Sigma+")->Mass()+pdg->GetParticle("proton")->Mass();
  double piLn_mass  = pdg->GetParticle("pi+")->Mass()+pdg->GetParticle("Lambda0")->Mass()+pdg->GetParticle("neutron")->Mass();
  double Kp_mass    = pdg->GetParticle("proton")->Mass()+pdg->GetParticle("K-")->Mass();
  double piS_mass1 = pdg->GetParticle("pi+")->Mass()+pdg->GetParticle("Sigma-")->Mass();
  double piS_mass2 = pdg->GetParticle("pi-")->Mass()+pdg->GetParticle("Sigma+")->Mass();
  double piL_mass  = pdg->GetParticle("pi+")->Mass()+pdg->GetParticle("Lambda0")->Mass();

  
  //= = = = pipipnn final-sample tree = = = =//
  TLorentzVector *mom_beam = 0;   // 4-momentum(beam)
  TLorentzVector *mom_target = 0; // 4-momentum(target)
  TLorentzVector *mom_pip = 0;    // 4-momentum(pi+)
  TLorentzVector *mom_pim = 0;    // 4-momentum(pi-)
  TLorentzVector *mom_p = 0;      // 4-momentum(proton)
  TLorentzVector *mom_n = 0;      // 4-momentum(neutron)
  double beta; // veracity of neutral particle on CDH
  double dE;   // energy deposit on CDH
  TVector3 *vtx_reaction = 0; // vertex(reaction)
  int run_num;   // run number
  int event_num; // event number
  int block_num; // block number
  TLorentzVector *kf1mom_beam = 0;   // 4-momentum(beam) after kinematical refit for pi- Sigma+
  TLorentzVector *kf1mom_pip = 0;    // 4-momentum(pi+) after kinematical refit for pi- Sigma+
  TLorentzVector *kf1mom_pim = 0;    // 4-momentum(pi-) after kinematical refit for pi- Sigma+
  TLorentzVector *kf1mom_p = 0;      // 4-momentum(proton) after kinematical refit for pi- Sigma+
  TLorentzVector *kf1mom_n = 0;      // 4-momentum(neutron) after kinematical refit for pi- Sigma+
  TLorentzVector *mcmom_beam = 0;   // generated 4-momentum(beam)
  TLorentzVector *mcmom_pip = 0;    // generated 4-momentum(pi+)
  TLorentzVector *mcmom_pim = 0;    // generated 4-momentum(pi-)
  TLorentzVector *mcmom_p = 0;      // generated 4-momentum(proton)
  TLorentzVector *mcmom_n = 0;      // generated 4-momentum(neutron)
  double kf1_chi2;   // chi2 of kinematical refit
  double kf1_NDF;    // NDF of kinematical refit
  double kf1_status; // status of kinematical refit -> details can be found in this code
  double kf1_pvalue; // p-value of kinematical refit
  TLorentzVector *kf2mom_beam = 0;   // 4-momentum(beam) after kinematical refit for pi+ Sigma-
  TLorentzVector *kf2mom_pip = 0;    // 4-momentum(pi+) after kinematical refit for pi+ Sigma-
  TLorentzVector *kf2mom_pim = 0;    // 4-momentum(pi-) after kinematical refit for pi+ Sigma-
  TLorentzVector *kf2mom_p = 0;      // 4-momentum(proton) after kinematical refit for pi+ Sigma-
  TLorentzVector *kf2mom_n = 0;      // 4-momentum(neutron) after kinematical refit for pi+ Sigma-
  double kf2_chi2;   // chi2 of kinematical refit
  double kf2_NDF;    // NDF of kinematical refit
  double kf2_status; // status of kinematical refit -> details can be found in this code
  double kf2_pvalue; // p-value of kinematical refit
  int kf_flag; // flag of correct pair reconstruction, etc
  //= = = = pipipnn final-sample tree = = = =//
  
  //TFile *f = new TFile("run65_piLnn.root");
  //TFile *f = new TFile("sim_piSpn.root");
  //TFile *f = new TFile("sim_piSpn_dE0_Al.root");
  //TFile *f = new TFile("L1405pn.root");
  TFile *f = new TFile("test.root");
  TTree *tree = (TTree*)f->Get("EventTree");

  tree->SetBranchAddress( "mom_beam",   &mom_beam );
  tree->SetBranchAddress( "mom_target", &mom_target );
  tree->SetBranchAddress( "mom_pip", &mom_pip );
  tree->SetBranchAddress( "mom_pim", &mom_pim );
  tree->SetBranchAddress( "mom_p", &mom_p );
  tree->SetBranchAddress( "mom_n", &mom_n );
  tree->SetBranchAddress( "beta", &beta );
  tree->SetBranchAddress( "dE", &dE );
  tree->SetBranchAddress( "vtx_reaction", &vtx_reaction );
  tree->SetBranchAddress( "run_num", &run_num );
  tree->SetBranchAddress( "event_num", &event_num );
  tree->SetBranchAddress( "block_num", &block_num );
  tree->SetBranchAddress( "mcmom_beam",   &mcmom_beam );
  tree->SetBranchAddress( "mcmom_pip", &mcmom_pip );
  tree->SetBranchAddress( "mcmom_pim", &mcmom_pim );
  tree->SetBranchAddress( "mcmom_p", &mcmom_p );
  tree->SetBranchAddress( "mcmom_n", &mcmom_n );
  tree->SetBranchAddress( "kf1mom_beam",   &kf1mom_beam );
  tree->SetBranchAddress( "kf1mom_pip", &kf1mom_pip );
  tree->SetBranchAddress( "kf1mom_pim", &kf1mom_pim );
  tree->SetBranchAddress( "kf1mom_p", &kf1mom_p );
  tree->SetBranchAddress( "kf1mom_n", &kf1mom_n );
  tree->SetBranchAddress( "kf1_chi2", &kf1_chi2 );
  tree->SetBranchAddress( "kf1_NDF", &kf1_NDF );
  tree->SetBranchAddress( "kf1_status", &kf1_status );
  tree->SetBranchAddress( "kf1_pvalue", &kf1_pvalue );
  tree->SetBranchAddress( "kf2mom_beam",   &kf2mom_beam );
  tree->SetBranchAddress( "kf2mom_pip", &kf2mom_pip );
  tree->SetBranchAddress( "kf2mom_pim", &kf2mom_pim );
  tree->SetBranchAddress( "kf2mom_p", &kf2mom_p );
  tree->SetBranchAddress( "kf2mom_n", &kf2mom_n );
  tree->SetBranchAddress( "kf2_chi2", &kf2_chi2 );
  tree->SetBranchAddress( "kf2_NDF", &kf2_NDF );
  tree->SetBranchAddress( "kf2_status", &kf2_status );
  tree->SetBranchAddress( "kf2_pvalue", &kf2_pvalue );
  tree->SetBranchAddress( "kf_flag", &kf_flag );

  TLegend *leg;
  TLine *line;
  double ymax;
  TH1F *htmp1;
  char com[32];

  const double beta_MAX = 0.728786; // p = 1.0 GeV/c for neutron & 1/beta = 1.372
  //const double dE_MIN = 5.0; // 8.0MeVee * 3cm / 5cm;
  const double dE_MIN = 2.0;

  const double pipi_MIN = 0.485;
  const double pipi_MAX = 0.510;
  const double ppi_MIN = 1.1075;
  const double ppi_MAX = 1.1225;

  const double neutron_MIN = 0.85;
  const double neutron_MAX = 1.03;

  const double CHI2 = 10; // maximum value
  const double PVALUE[9] = {1E-60, 0.0001, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.5}; // minimum value

  const int NBINS_orig = 20;
  int NBINS = NBINS_orig;

  TH1F *IMnpipi[NBINS];
  TH1F *IMnppipi[NBINS];

  TH2F *Diff_mom_n[4];

  const int    MISSN_BIN = 70;
  const double MISSN_MIN = 0.4;
  const double MISSN_MAX = 1.8;

  const int    SIGMA_BIN = 70;
  const double SIGMA_MIN = 1.0;
  const double SIGMA_MAX = 1.7;

#if 0
  const int    PISIGMA_BIN = 50;
  const double PISIGMA_MIN = 1.0;
  const double PISIGMA_MAX = 2.0;
#else
  const int    PISIGMA_BIN = 66;
  const double PISIGMA_MIN = 1.0;
  const double PISIGMA_MAX = 1.99;
#endif

  const int    PISIGMAP_BIN = 50;
  const double PISIGMAP_MIN = 2.0;
  const double PISIGMAP_MAX = 3.0;

  const int    COSN_BIN = 40;
  const double COSN_MIN = -1.0;
  const double COSN_MAX = 1.0;

  for( int x=0; x<NBINS; x++ ){
    sprintf( com, "IMnpipi_%d", x );
    IMnpipi[x] = new TH1F( com, com, 160, -80, 80 ); // MeV
    sprintf( com, "IMnppipi_%d", x );
    IMnppipi[x] = new TH1F( com, com, 160, -80, 80 ); // MeV
  }

  int TGT = 4;

  for( int i=0; i<4; i++ ){
    sprintf(com, "Diff_mom_n_%d (cut_%d)", i, TGT);
    Diff_mom_n[i] = new TH2F(com,com,100,-0.2,0.2,100,0,0.8);
    sprintf(com, "Gen_mom [GeV/c]");
    Diff_mom_n[i]->SetYTitle(com);
    sprintf(com, "(Rec_mom-Gen_mom)_%d [GeV/c] (cut_%d)", i, TGT);
    Diff_mom_n[i]->SetXTitle(com);
  }

    
  double sta[2][NBINS];
  double end[2][NBINS];
  double mass[2][NBINS];
  double mass_e[2][NBINS];
  double mass_cen[2][NBINS];
  double mass_cen_e[2][NBINS];
  double mass_res[2][NBINS];
  double mass_res_e[2][NBINS];
  for( int x=0; x<NBINS; x++ ){
    sta[0][x] = PISIGMA_MIN+(PISIGMA_MAX-PISIGMA_MIN)/NBINS*x;
    end[0][x] = PISIGMA_MIN+(PISIGMA_MAX-PISIGMA_MIN)/NBINS*(x+1);
    mass[0][x] = PISIGMA_MIN+(PISIGMA_MAX-PISIGMA_MIN)/NBINS*(x+0.5);
    mass_e[0][x] = (PISIGMA_MAX-PISIGMA_MIN)/NBINS/2;
    sta[1][x] = PISIGMAP_MIN+(PISIGMAP_MAX-PISIGMAP_MIN)/NBINS*x;
    end[1][x] = PISIGMAP_MIN+(PISIGMAP_MAX-PISIGMAP_MIN)/NBINS*(x+1);
    mass[1][x] = PISIGMAP_MIN+(PISIGMAP_MAX-PISIGMAP_MIN)/NBINS*(x+0.5);
    mass_e[1][x] = (PISIGMAP_MAX-PISIGMAP_MIN)/NBINS/2;
    for( int y=0; y<2; y++ ){
      mass_cen[y][x] = 0;
      mass_cen_e[y][x] = 0;
      mass_res[y][x] = 0;
      mass_res_e[y][x] = 0;
    }
  }

  Int_t nevent = tree->GetEntries();
  std::cerr<<"# of events = "<<nevent<<std::endl;
  //------------------------//
  //--- event roop start ---//
  //------------------------//
  for ( Int_t i=0; i<nevent; i++ ) {
    tree->GetEvent(i);
    //std::cerr<<i<<" "<<run_num<<" "<<event_num<<" "<<block_num<<std::endl;

    //** calc pi-p **//
    TLorentzVector mom_pim_p = *mom_pim+*mom_p;

    //** calc pi+pi- **//
    TLorentzVector mom_pip_pim = *mom_pip+*mom_pim;

    //** calc pi+n **//
    TLorentzVector mom_pip_n = *mom_pip+*mom_n;

    //** calc pi-n **//
    TLorentzVector mom_pim_n = *mom_pim+*mom_n;

    //** calc pi+pi-n **//
    TLorentzVector mom_pip_pim_n = *mom_pip+*mom_pim+*mom_n;
    TLorentzVector mcmom_pip_pim_n = *mcmom_pip+*mcmom_pim+*mcmom_n;
    double diff_pip_pim_n = (mom_pip_pim_n.M()-mcmom_pip_pim_n.M())*1000; // MeV

    //** calc pi+pi-np **//
    TLorentzVector mom_pip_pim_n_p = *mom_pip+*mom_pim+*mom_n+*mom_p;
    TLorentzVector mcmom_pip_pim_n_p = *mcmom_pip+*mcmom_pim+*mcmom_n+*mcmom_p;
    double diff_pip_pim_n_p = (mom_pip_pim_n_p.M()-mcmom_pip_pim_n_p.M())*1000; // MeV

    //** calc detected-n **//
    TLorentzVector diff_n = *mom_n-*mcmom_n;
    
    //-- neutron-ID, K0 & Lambda subtraction --//
    if( beta<beta_MAX && dE_MIN<dE &&
	(mom_pip_pim.M()<pipi_MIN || pipi_MAX<mom_pip_pim.M()) &&
	(mom_pim_p.M()<ppi_MIN || ppi_MAX<mom_pim_p.M()) ){
      if( -1<kf_flag ){
	if( (kf1_status==0 && kf2_pvalue<kf1_pvalue && PVALUE[TGT]<kf1_pvalue) ||
	    (kf2_status==0 && kf1_pvalue<kf2_pvalue && PVALUE[TGT]<kf2_pvalue) ){
	  for( int x=0; x<NBINS; x++ ){
	    if( sta[0][x]<mcmom_pip_pim_n.M() && mcmom_pip_pim_n.M()<end[0][x] ){
	      IMnpipi[x]->Fill(diff_pip_pim_n);
	    }
	    if( sta[1][x]<mcmom_pip_pim_n_p.M() && mcmom_pip_pim_n_p.M()<end[1][x] ){
	      IMnppipi[x]->Fill(diff_pip_pim_n_p);
	    }
	  }
	  Diff_mom_n[0]->Fill(diff_n.X(),mcmom_n->P());
	  Diff_mom_n[1]->Fill(diff_n.Y(),mcmom_n->P());
	  Diff_mom_n[2]->Fill(diff_n.Z(),mcmom_n->P());
	  Diff_mom_n[3]->Fill(mom_n->E()-mcmom_n->E(),mcmom_n->P());
	  //std::cerr<<mom_n->E()-mcmom_n->E()<<std::endl;
	}
      }
    } // if( beta<beta_MAX && dE_MIN<dE &&

  } // for ( Int_t i=0; i<nevent; i++ ) {


  //** ++ ** ++ ** ++ ** ++ ** ++ ** ++ ** ++ ** ++ **//

  TH1F *his1, *his2, *his3, *his4, *his5, *his6;
  double width;
#if 0
  TCanvas *c0;
  c0 = new TCanvas("c0", "", 600, 600);
  c0->Divide(4,5);
  for( int x=0; x<NBINS; x++ ){
    c0->cd(x+1);
    IMnpipi[x]->Draw();
    IMnpipi[x]->SetXTitle("measured-generated mass-resolution [MeV/c^{2}]");
    if( 200<IMnpipi[x]->GetEntries() ){
      IMnpipi[x]->Fit("gaus","q");
      IMnpipi[x]->GetFunction("gaus")->SetLineColor(4);
      mass_cen[0][x] = IMnpipi[x]->GetFunction("gaus")->GetParameter(1);
      mass_cen_e[0][x] = IMnpipi[x]->GetFunction("gaus")->GetParError(1);
      mass_res[0][x] = IMnpipi[x]->GetFunction("gaus")->GetParameter(2);
      mass_res_e[0][x] = IMnpipi[x]->GetFunction("gaus")->GetParError(2);
    }
  }
  c0->Print("tmp.pdf(");

#if DEBUG
  for( int x=0; x<NBINS; x++ ){
    std::cerr<<sta[0][x]<<" - "<<end[0][x]<<" | " 
	     <<mass_cen[0][x]<<" +/- "<<mass_cen_e[0][x]<<" , "
	     <<mass_res[0][x]<<" +/- "<<mass_res_e[0][x]<<std::endl;
  }
#endif

  for( int x=0; x<NBINS; x++ ){
    if( mass_cen[0][x]==0 ){
      for( int y=x; y<NBINS; y++ ){
	sta[0][y] = sta[0][y+1];
	end[0][y] = end[0][y+1];
	mass[0][y] = mass[0][y+1];
	mass_cen[0][y] = mass_cen[0][y+1];
	mass_cen_e[0][y] = mass_cen_e[0][y+1];
	mass_res[0][y] = mass_res[0][y+1];
	mass_res_e[0][y] = mass_res_e[0][y+1];
      }
      x--;
      NBINS--;
    }
  }

#if DEBUG
  std::cerr<<" %%%%%%%%%%%%% "<<std::endl;
  for( int x=0; x<NBINS; x++ ){
    std::cerr<<sta[0][x]<<" - "<<end[0][x]<<" | " 
	     <<mass_cen[0][x]<<" +/- "<<mass_cen_e[0][x]<<" , "
	     <<mass_res[0][x]<<" +/- "<<mass_res_e[0][x]<<std::endl;
  }  
#endif

  TGraphErrors *gr1[2];
  TGraphErrors *gr2[2];
  gr1[0] = new TGraphErrors( NBINS, mass[0], mass_cen[0], mass_e[0], mass_cen_e[0] );
  gr2[0] = new TGraphErrors( NBINS, mass[0], mass_res[0], mass_e[0], mass_res_e[0] );


  TCanvas *c1;
  c1 = new TCanvas("c1", "", 600, 600);
  gr1[0]->GetXaxis()->SetTitle("IM(n#pi^{+}#pi^{-}) [MeV/c^{2}]");
  gr1[0]->GetYaxis()->SetTitle("measured-generated mass-center [MeV/c^{2}]");
  gr1[0]->SetMarkerStyle(21);
  gr1[0]->SetMarkerColor(1);
  gr1[0]->SetMarkerSize(1);
  //gr1[0]->GetXaxis()->SetNdivisions(6);
  //gr1[0]->GetXaxis()->SetLabelOffset(99);
  //gr1[0]->GetYaxis()->SetTitle("[%]");
  //gr1[0]->GetYaxis()->SetTitleOffset(1.3);
  //gr1[0]->GetYaxis()->SetRangeUser(-2.0, 4.0);
  //gr1[0]->Draw("P");
  gr1[0]->GetXaxis()->SetLimits(PISIGMA_MIN, PISIGMA_MAX);
  gr1[0]->SetMinimum(-20.0);
  gr1[0]->SetMaximum(20.0);
  gr1[0]->SetTitle("mass center");
  gr1[0]->Draw("AP");
  line = new TLine(PISIGMA_MIN, 0, PISIGMA_MAX, 0);
  line->SetLineColor(2);
  line->SetLineStyle(2);
  line->Draw();
  line = new TLine(Kp_mass, -20, Kp_mass, 20);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->Draw();
  line = new TLine(piS_mass1, -20, piS_mass1, 20);
  line->SetLineColor(4);
  line->SetLineStyle(2);
  line->Draw();
  line = new TLine(piS_mass2, -20, piS_mass2, 20);
  line->SetLineColor(4);
  line->SetLineStyle(2);
  line->Draw();
  //c1->SetGrid();
  c1->Print("tmp.pdf");

  TCanvas *c2;
  c2 = new TCanvas("c2", "", 600, 600);
  gr2[0]->GetXaxis()->SetTitle("IM(n#pi^{+}#pi^{-}) [MeV/c^{2}]");
  gr2[0]->GetYaxis()->SetTitle("measured-generated mass-resolution [MeV/c^{2}]");
  gr2[0]->SetMarkerStyle(21);
  gr2[0]->SetMarkerColor(1);
  gr2[0]->SetMarkerSize(1);
  //gr2[0]->GetXaxis()->SetNdivisions(6);
  //gr2[0]->GetXaxis()->SetLabelOffset(99);
  //gr2[0]->GetYaxis()->SetTitle("[%]");
  //gr2[0]->GetYaxis()->SetTitleOffset(1.3);
  //gr2[0]->GetYaxis()->SetRangeUser(-2.0, 4.0);
  gr2[0]->GetXaxis()->SetLimits(PISIGMA_MIN, PISIGMA_MAX);
  gr2[0]->SetMinimum(0.0);
  gr2[0]->SetMaximum(40.0);
  gr2[0]->SetTitle("mass resolution");
  gr2[0]->Draw("AP");
  line = new TLine(Kp_mass, 0, Kp_mass, 40);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->Draw();
  line = new TLine(piS_mass1, 0, piS_mass1, 40);
  line->SetLineColor(4);
  line->SetLineStyle(2);
  line->Draw();
  line = new TLine(piS_mass2, 0, piS_mass2, 40);
  line->SetLineColor(4);
  line->SetLineStyle(2);
  line->Draw();
  //c2->SetGrid();
  gr2[0]->Fit("pol3","","",piS_mass1,1.9);
  gr2[0]->GetFunction("pol3")->SetLineColor(4);
  c2->Print("tmp.pdf)");




  //=============//

  NBINS = NBINS_orig;

  TCanvas *d0;
  d0 = new TCanvas("d0", "", 600, 600);
  d0->Divide(4,5);
  for( int x=0; x<NBINS; x++ ){
    d0->cd(x+1);
    IMnppipi[x]->Draw();
    IMnppipi[x]->SetXTitle("measured-generated mass-resolution [MeV/c^{2}]");
    if( 200<IMnppipi[x]->GetEntries() ){
      IMnppipi[x]->Fit("gaus","q");
      IMnppipi[x]->GetFunction("gaus")->SetLineColor(4);
      mass_cen[1][x] = IMnppipi[x]->GetFunction("gaus")->GetParameter(1);
      mass_cen_e[1][x] = IMnppipi[x]->GetFunction("gaus")->GetParError(1);
      mass_res[1][x] = IMnppipi[x]->GetFunction("gaus")->GetParameter(2);
      mass_res_e[1][x] = IMnppipi[x]->GetFunction("gaus")->GetParError(2);
    }
  }
  d0->Print("tmp2.pdf(");

#if DEBUG
  for( int x=0; x<NBINS; x++ ){
    std::cerr<<sta[1][x]<<" - "<<end[1][x]<<" | " 
	     <<mass_cen[1][x]<<" +/- "<<mass_cen_e[1][x]<<" , "
	     <<mass_res[1][x]<<" +/- "<<mass_res_e[1][x]<<std::endl;
  }
#endif

  for( int x=0; x<NBINS; x++ ){
    if( mass_cen[1][x]==0 ){
      for( int y=x; y<NBINS; y++ ){
	sta[1][y] = sta[1][y+1];
	end[1][y] = end[1][y+1];
	mass[1][y] = mass[1][y+1];
	mass_cen[1][y] = mass_cen[1][y+1];
	mass_cen_e[1][y] = mass_cen_e[1][y+1];
	mass_res[1][y] = mass_res[1][y+1];
	mass_res_e[1][y] = mass_res_e[1][y+1];
      }
      x--;
      NBINS--;
    }
  }

#if DEBUG
  std::cerr<<" %%%%%%%%%%%%% "<<std::endl;
  for( int x=0; x<NBINS; x++ ){
    std::cerr<<sta[1][x]<<" - "<<end[1][x]<<" | " 
	     <<mass_cen[1][x]<<" +/- "<<mass_cen_e[1][x]<<" , "
	     <<mass_res[1][x]<<" +/- "<<mass_res_e[1][x]<<std::endl;
  }  
#endif

  gr1[1] = new TGraphErrors( NBINS, mass[1], mass_cen[1], mass_e[1], mass_cen_e[1] );
  gr2[1] = new TGraphErrors( NBINS, mass[1], mass_res[1], mass_e[1], mass_res_e[1] );


  TCanvas *d1;
  d1 = new TCanvas("d1", "", 600, 600);
  gr1[1]->GetXaxis()->SetTitle("IM(np#pi^{+}#pi^{-}) [MeV/c^{2}]");
  gr1[1]->GetYaxis()->SetTitle("measured-generated mass-center [MeV/c^{2}]");
  gr1[1]->SetMarkerStyle(21);
  gr1[1]->SetMarkerColor(1);
  gr1[1]->SetMarkerSize(1);
  //gr1[1]->GetXaxis()->SetNdivisions(6);
  //gr1[1]->GetXaxis()->SetLabelOffset(99);
  //gr1[1]->GetYaxis()->SetTitle("[%]");
  //gr1[1]->GetYaxis()->SetTitleOffset(1.3);
  //gr1[1]->GetYaxis()->SetRangeUser(-2.0, 4.0);
  //gr1[1]->Draw("P");
  gr1[1]->GetXaxis()->SetLimits(PISIGMAP_MIN, PISIGMAP_MAX);
  gr1[1]->SetMinimum(-20.0);
  gr1[1]->SetMaximum(20.0);
  gr1[1]->SetTitle("mass center");
  gr1[1]->Draw("AP");
  line = new TLine(PISIGMAP_MIN, 0, PISIGMAP_MAX, 0);
  line->SetLineColor(2);
  line->SetLineStyle(2);
  line->Draw();
  line = new TLine(Kpp_mass, -20, Kpp_mass, 20);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->Draw();
  line = new TLine(piSp_mass1, -20, piSp_mass1, 20);
  line->SetLineColor(4);
  line->SetLineStyle(2);
  line->Draw();
  line = new TLine(piSp_mass2, -20, piSp_mass2, 20);
  line->SetLineColor(4);
  line->SetLineStyle(2);
  line->Draw();
  //d1->SetGrid();
  d1->Print("tmp2.pdf");

  TCanvas *d2;
  d2 = new TCanvas("d2", "", 600, 600);
  gr2[1]->GetXaxis()->SetTitle("IM(np#pi^{+}#pi^{-}) [MeV/c^{2}]");
  gr2[1]->GetYaxis()->SetTitle("measured-generated mass-resolution [MeV/c^{2}]");
  gr2[1]->SetMarkerStyle(21);
  gr2[1]->SetMarkerColor(1);
  gr2[1]->SetMarkerSize(1);
  //gr2[1]->GetXaxis()->SetNdivisions(6);
  //gr2[1]->GetXaxis()->SetLabelOffset(99);
  //gr2[1]->GetYaxis()->SetTitle("[%]");
  //gr2[1]->GetYaxis()->SetTitleOffset(1.3);
  //gr2[1]->GetYaxis()->SetRangeUser(-2.0, 4.0);
  gr2[1]->GetXaxis()->SetLimits(PISIGMAP_MIN, PISIGMAP_MAX);
  gr2[1]->SetMinimum(0.0);
  gr2[1]->SetMaximum(40.0);
  gr2[1]->SetTitle("mass resolution");
  gr2[1]->Draw("AP");
  line = new TLine(Kpp_mass, 0, Kpp_mass, 40);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->Draw();
  line = new TLine(piSp_mass1, 0, piSp_mass1, 40);
  line->SetLineColor(4);
  line->SetLineStyle(2);
  line->Draw();
  line = new TLine(piSp_mass2, 0, piSp_mass2, 40);
  line->SetLineColor(4);
  line->SetLineStyle(2);
  line->Draw();
  //d2->SetGrid();
  d2->Print("tmp2.pdf)");
#endif

  TCanvas *d3;
  d3 = new TCanvas("d3", "", 600, 600);
  d3->Divide(2,2);
  d3->cd(1); Diff_mom_n[0]->Draw("colz");
  d3->cd(2); Diff_mom_n[1]->Draw("colz");
  d3->cd(3); Diff_mom_n[2]->Draw("colz");
  d3->cd(4); Diff_mom_n[3]->Draw("colz");
  d3->Print("tmp3.pdf");
}
