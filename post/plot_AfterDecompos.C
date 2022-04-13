//input : output files of the macro Fit2DK0.C 
//      : acceptance map
//
//output : decoposed spectra q vs IM(npi+pi-) and acceptance corrected spectra

const bool showBG = true;
const double UncertCut = 0.25;
#include "anacuts.h"

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
#include <TParameter.h>

#include "../src/GlobalVariables.h"

bool RemoveNegative = true;
bool RebinMode = true;
bool Sidefar=false;
bool FitNoWeight=true;

const int Version = 241;
const int versionSigma = 156;//SIM version
const int versionK0 = 30;//SIM version

void plot_AfterDecompos(const int dEcut=2,const int sysud=0)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetMarkerStyle(20); 
  gStyle->SetMarkerSize(1.2); 
  gROOT->ForceStyle();
  gROOT->SetBatch();

  TFile *fr[4] = {NULL};
  //Because the statistics is limited, we divide data into q<0.35 and q>0.35 and decompose K0 & Sigma+ & Simga- 
  fr[0] = TFile::Open(Form("evanaIMpisigma_npippim_v%d_out_dE%d_iso_nostop_sys%d_sub.root",Version,dEcut,sysud),"READ");
  fr[1] = TFile::Open(Form("evanaIMpisigma_npippim_v%d_out_dE%d_iso_qlo_nostop_sys%d_sub.root",Version,dEcut,sysud),"READ");
  fr[2] = TFile::Open(Form("evanaIMpisigma_npippim_v%d_out_dE%d_iso_qhi_nostop_sys%d_sub.root",Version,dEcut,sysud),"READ");
  fr[3] = TFile::Open(Form("evanaIMpisigma_npippim_v%d_out_dE%d_iso_theta15_nostop_sys%d_sub.root",Version,dEcut,sysud),"READ");
  
  std::cout << "FILE OPEN " << std::endl;
 

  TFile *fdeco[2]={NULL,NULL};
  fdeco[0] = TFile::Open(Form("fout_qlo_v%d_dE%d_sys%d.root",Version,2,0),"READ");
  fdeco[1] = TFile::Open(Form("fout_qhi_v%d_dE%d_sys%d.root",Version,2,0),"READ");
  fdeco[0]->Print();
  fdeco[1]->Print();

  TFile *flumi = TFile::Open("InteLumi.root");
  TParameter<double>*IntegLumi = (TParameter<double>*)flumi->Get("IntegLumi");
  TParameter<double>*Err = (TParameter<double>*)flumi->Get("Err");
  double lumi = IntegLumi->GetVal();
  double lumierr = Err->GetVal();
  const double trigScale = 0.5;
  std::cout << "Lumi:  " << lumi << std::endl;
  std::cout << "Err:   " << lumierr << std::endl;
  
  //gROOT->SetBatch(1);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TH1::SetDefaultSumw2();
  //gStyle->SetErrorX(0.);  

  const unsigned int nbinIMnpipi = 80;
  const int nqcut=4;
  const int qstart=0;

  //for the overlap of S+ & S- & K0 counting 
  TH2F* q_IMnpipi_wK0_wSid_n_SpSm[nqcut];
  TH1D* IMnpipi_wK0_wSid_n_SpSm[nqcut];
  TCanvas *cq_IMnpipi_wK0_wSid_n_SpSm[nqcut];
  TCanvas *cIMnpipi_Sp[nqcut];
  TCanvas *cq_IMnpipi_wSid_n_Sp[nqcut];
  TCanvas *cq_IMnpipi_wK0_wSid_n_Sp[nqcut];
  TCanvas *cIMnpipi_Sm[nqcut];
  TCanvas *cq_IMnpipi_wSid_n_Sm[nqcut];
  TCanvas *cq_IMnpipi_wK0_wSid_n_Sm[nqcut];
  TCanvas *cq_IMnpipi_wSid_n_SpSm[nqcut];
  TCanvas *cIMnpipi_K0[nqcut];
  TCanvas *cq_IMnpipi_K0_n[nqcut];
  double OverlapCount[nbinIMnpipi][nqcut]; 
  const unsigned int nwbin = 3;
  TH2D* IMnpim_IMnpip_dE_wK0_wSid_n_Sp_bin[nbinIMnpipi][nqcut];
  TH2D* IMnpim_IMnpip_dE_wK0_wSid_n_Sm_bin[nbinIMnpipi][nqcut];
  TH2D* IMnpim_IMnpip_dE_wSid_n_SpSm[nqcut];
  TH2D* IMnpim_IMnpip_dE_wSid_n_SpSm_bin[nbinIMnpipi][nqcut];
  TH2D* IMnpim_IMnpip_dE_wK0_wSid_n_Sp_wbin[nwbin][nqcut];
  TH2D* IMnpim_IMnpip_dE_wK0_wSid_n_Sm_wbin[nwbin][nqcut];
  TH2D* q_IMnpipi_wSid_n_Sp[nqcut];
  TH1D* IMnpipi_wSid_n_Sp[nqcut];
  TH2D* q_IMnpipi_wK0_wSid_n_Sp[nqcut];
  TH1D* IMnpipi_wK0_wSid_n_Sp[nqcut];
  TH1D* IMnpipi_wK0_wSid_n_Sp2[nqcut];
  TH2D* q_IMnpipi_wSid_n_Sm[nqcut];
  TH1D* IMnpipi_wSid_n_Sm[nqcut];
  TH2D* q_IMnpipi_wK0_wSid_n_Sm[nqcut];
  TH1D* IMnpipi_wK0_wSid_n_Sm[nqcut];
  TH1D* IMnpipi_wK0_wSid_n_Sm2[nqcut];
  TH2D* q_IMnpipi_wSid_n_SpSm[nqcut];
  TH1D* IMnpipi_wSid_n_SpSm[nqcut];
  TH2D* q_IMnpipi_wK0_n[nqcut];
  TH1D* IMnpipi_wK0_n[nqcut];

  const char cqcut[][10]= {"all","qlo","qhi","theta"};
  std::cout << __LINE__ << std::endl;
  for(int iq=qstart;iq<nqcut;iq++){
    //get hists from root files
    q_IMnpipi_wK0_wSid_n_SpSm[iq] = (TH2F*)fr[iq]->Get("q_IMnpipi_wK0_wSid_n_SpSm");
    q_IMnpipi_wK0_wSid_n_SpSm[iq]->RebinX(3);
    q_IMnpipi_wSid_n_Sp[iq] = (TH2D*)fr[iq]->Get("q_IMnpipi_wSid_n_Sp");
    q_IMnpipi_wSid_n_Sp[iq]->RebinX(3);
    q_IMnpipi_wK0_wSid_n_Sp[iq] = (TH2D*)fr[iq]->Get("q_IMnpipi_wK0_wSid_n_Sp");
    q_IMnpipi_wK0_wSid_n_Sp[iq]->RebinX(3);
    q_IMnpipi_wSid_n_Sm[iq] = (TH2D*)fr[iq]->Get("q_IMnpipi_wSid_n_Sm");
    q_IMnpipi_wSid_n_Sm[iq]->RebinX(3);
    q_IMnpipi_wK0_wSid_n_Sm[iq] = (TH2D*)fr[iq]->Get("q_IMnpipi_wK0_wSid_n_Sm");
    q_IMnpipi_wK0_wSid_n_Sm[iq]->RebinX(3);
    q_IMnpipi_wSid_n_SpSm[iq] = (TH2D*)fr[iq]->Get("q_IMnpipi_wSid_n_SpSm");
    q_IMnpipi_wSid_n_SpSm[iq]->RebinX(3);
    q_IMnpipi_wK0_n[iq] = (TH2D*)fr[iq]->Get("q_IMnpipi_wK0_n");
    q_IMnpipi_wK0_n[iq]->RebinX(3);

    //Draw plots for comfirmation
    cq_IMnpipi_wSid_n_Sp[iq] = new TCanvas(Form("cq_IMnpipi_wSid_n_Sp_%s",cqcut[iq]),Form("cq_IMnpipi_wSid_n_Sp_%s",cqcut[iq]));
    q_IMnpipi_wSid_n_Sp[iq]->Draw("colz");
    cq_IMnpipi_wK0_wSid_n_Sp[iq] = new TCanvas(Form("cq_IMnpipi_wK0_wSid_n_Sp_%s",cqcut[iq]),Form("cq_IMnpipi_wK0_wSid_n_Sp_%s",cqcut[iq]));
    q_IMnpipi_wK0_wSid_n_Sp[iq]->Draw("colz");
    
    cIMnpipi_Sp[iq] = new TCanvas(Form("IMnpipi_Sp_%s",cqcut[iq]),Form("IMnpipi_Sp_%s",cqcut[iq]),1000,800);
    IMnpipi_wSid_n_Sp[iq] = (TH1D*)q_IMnpipi_wSid_n_Sp[iq]->ProjectionX(Form("IMnpipi_wSid_n_Sp_%d",iq));
    IMnpipi_wSid_n_Sp[iq]->SetTitle(Form("IM(n#pi^{+}#pi^{-}) #Sigma^{+}#pi^{-} like mode %s",cqcut[iq]));
    IMnpipi_wSid_n_Sp[iq]->SetYTitle("Counts/15 MeV");
    IMnpipi_wSid_n_Sp[iq]->GetYaxis()->CenterTitle();
    IMnpipi_wSid_n_Sp[iq]->Draw("E");
    IMnpipi_wK0_wSid_n_Sp[iq] = (TH1D*)q_IMnpipi_wK0_wSid_n_Sp[iq]->ProjectionX(Form("IMnpipi_wK0_wSid_n_Sp_%d",iq));
    IMnpipi_wK0_wSid_n_Sp[iq]->SetLineColor(2);
    IMnpipi_wK0_wSid_n_Sp[iq]->SetMarkerColor(2);
    IMnpipi_wK0_wSid_n_Sp[iq]->Draw("Esame");
    IMnpipi_wSid_n_SpSm[iq] = (TH1D*)q_IMnpipi_wSid_n_SpSm[iq]->ProjectionX(Form("IMnpipi_wSid_n_SpSm_%d",iq));
    IMnpipi_wSid_n_SpSm[iq]->SetLineColor(3);
    IMnpipi_wSid_n_SpSm[iq]->SetMarkerColor(3);
    IMnpipi_wSid_n_SpSm[iq]->Draw("Esame");
    
    TLegend *lSpoverlap = new TLegend(0.6,0.7,0.9,0.9);
    lSpoverlap->AddEntry(IMnpipi_wSid_n_Sp[iq],"#Sigma^{+} like events (total)","l");
    lSpoverlap->AddEntry(IMnpipi_wK0_wSid_n_Sp[iq],"w/ #bar{K}^{0} overlap","l");
    lSpoverlap->AddEntry(IMnpipi_wSid_n_SpSm[iq],"w/ #Sigma^{-} overlap","l");
    lSpoverlap->Draw();

    cq_IMnpipi_wSid_n_Sm[iq] = new TCanvas(Form("cq_IMnpipi_wSid_n_Sm_%s",cqcut[iq]),Form("cq_IMnpipi_wSid_n_Sm_%s",cqcut[iq]));
    q_IMnpipi_wSid_n_Sm[iq]->Draw("colz");
    cq_IMnpipi_wK0_wSid_n_Sm[iq] = new TCanvas(Form("cq_IMnpipi_wK0_wSid_n_Sm_%s",cqcut[iq]),Form("cq_IMnpipi_wK0_wSid_n_Sm_%s",cqcut[iq]));
    q_IMnpipi_wK0_wSid_n_Sm[iq]->Draw("colz");
    cIMnpipi_Sm[iq] = new TCanvas(Form("IMnpipi_Sm_%s",cqcut[iq]),Form("IMnpipi_Sm_%s",cqcut[iq]),1000,800);
    IMnpipi_wSid_n_Sm[iq] = (TH1D*)q_IMnpipi_wSid_n_Sm[iq]->ProjectionX(Form("IMnpipi_wSid_n_Sm_%d",iq));
    IMnpipi_wSid_n_Sm[iq]->SetTitle(Form("IM(n#pi^{+}#pi^{-}) #Sigma^{-}#pi^{+} like mode %s",cqcut[iq]));
    IMnpipi_wSid_n_Sm[iq]->SetYTitle("Counts/15 MeV");
    IMnpipi_wSid_n_Sm[iq]->GetYaxis()->CenterTitle();
    IMnpipi_wSid_n_Sm[iq]->Draw("E");
    IMnpipi_wK0_wSid_n_Sm[iq] = (TH1D*)q_IMnpipi_wK0_wSid_n_Sm[iq]->ProjectionX(Form("IMnpipi_wK0_wSid_n_Sm_%d",iq));
    IMnpipi_wK0_wSid_n_Sm[iq]->SetLineColor(2);
    IMnpipi_wK0_wSid_n_Sm[iq]->SetMarkerColor(2);
    IMnpipi_wK0_wSid_n_Sm[iq]->Draw("Esame");
    IMnpipi_wSid_n_SpSm[iq]->Draw("Esame");
    TLegend *lSmoverlap = new TLegend(0.6,0.7,0.9,0.9);
    lSmoverlap->AddEntry(IMnpipi_wSid_n_Sm[iq],"#Sigma^{-} like events (total)","l");
    lSmoverlap->AddEntry(IMnpipi_wK0_wSid_n_Sm[iq],"w/ #bar{K}^{0} overlap","l");
    lSmoverlap->AddEntry(IMnpipi_wSid_n_SpSm[iq],"w/ #Sigma^{+} overlap","l");
    lSmoverlap->Draw();
    
    cq_IMnpipi_wSid_n_SpSm[iq] = new TCanvas(Form("cq_IMnpipi_wSid_n_SpSm_%s",cqcut[iq]),Form("cq_IMnpipi_wSid_n_SpSm_%s",cqcut[iq]));
    q_IMnpipi_wSid_n_SpSm[iq]->Draw("colz");


    cIMnpipi_K0[iq] = new TCanvas(Form("IMnpipi_K0_%s",cqcut[iq]),Form("IMnpipi_K0_%s",cqcut[iq]),1000,800);
    IMnpipi_wK0_n[iq] = (TH1D*)q_IMnpipi_wK0_n[iq]->ProjectionX(Form("IMnpipi_wK0_n_%d",iq));
    IMnpipi_wK0_n[iq]->SetTitle(Form("IM(n#pi^{+}#pi^{-}) #bar{K}^{0}n like mode %s",cqcut[iq]));
    IMnpipi_wK0_n[iq]->SetYTitle("Counts/15 MeV");
    IMnpipi_wK0_n[iq]->GetYaxis()->CenterTitle();
    IMnpipi_wK0_n[iq]->Draw("E");
    IMnpipi_wK0_wSid_n_Sp2[iq] = (TH1D*)IMnpipi_wK0_wSid_n_Sp[iq]->Clone(Form("IMnpipi_wK0_wSid_n_Sp2_%d",iq));
    IMnpipi_wK0_wSid_n_Sm2[iq] = (TH1D*)IMnpipi_wK0_wSid_n_Sm[iq]->Clone(Form("IMnpipi_wK0_wSid_n_Sm2_%d",iq));
    IMnpipi_wK0_wSid_n_Sp2[iq]->SetLineColor(2);
    IMnpipi_wK0_wSid_n_Sp2[iq]->SetMarkerColor(2);
    IMnpipi_wK0_wSid_n_Sm2[iq]->SetLineColor(3);
    IMnpipi_wK0_wSid_n_Sm2[iq]->SetMarkerColor(3);
    IMnpipi_wK0_wSid_n_Sp2[iq]->Draw("Esame");
    IMnpipi_wK0_wSid_n_Sm2[iq]->Draw("Esame");
    
    TLegend *lK0overlap = new TLegend(0.6,0.7,0.9,0.9);
    lK0overlap->AddEntry(IMnpipi_wK0_n[iq],"#bar{K}^{0} like events (total)","l");
    lK0overlap->AddEntry(IMnpipi_wK0_wSid_n_Sp2[iq],"w/ #Sigma^{+} overlap","l");
    lK0overlap->AddEntry(IMnpipi_wK0_wSid_n_Sm2[iq],"w/ #Sigma^{-} overlap","l");
    lK0overlap->Draw();


    cq_IMnpipi_wK0_wSid_n_SpSm[iq] = new TCanvas(Form("q_IMnpipi_wK0_wSid_n_SpSm_%s",cqcut[iq]),Form("q_IMnpipi_wK0_wSid_n_SpSm_%s",cqcut[iq]));
    IMnpipi_wK0_wSid_n_SpSm[iq] = (TH1D*)q_IMnpipi_wK0_wSid_n_SpSm[iq]->ProjectionX(Form("IMnpipi_wK0_wSid_n_SpSm_%d",iq));
    IMnpipi_wK0_wSid_n_SpSm[iq]->SetTitle(Form("IMnpipi K0&Sp&Sm %s",cqcut[iq]));
    IMnpipi_wK0_wSid_n_SpSm[iq]->SetYTitle("Counts/15 MeV");
    IMnpipi_wK0_wSid_n_SpSm[iq]->GetYaxis()->CenterTitle();
    IMnpipi_wK0_wSid_n_SpSm[iq]->Draw("E");
     

    IMnpim_IMnpip_dE_wSid_n_SpSm[iq] = (TH2D*)fr[iq]->Get("IMnpim_IMnpip_dE_wSid_n_SpSm");
    //Sigma+ & Sigma- & K0 overlap
    for(int ibin=0;ibin<nbinIMnpipi;ibin++){
      IMnpim_IMnpip_dE_wK0_wSid_n_Sp_bin[ibin][iq] = (TH2D*)fr[iq]->Get(Form("IMnpim_IMnpip_dE_wK0_wSid_n_Sp_bin%d",ibin));
      IMnpim_IMnpip_dE_wK0_wSid_n_Sm_bin[ibin][iq] = (TH2D*)fr[iq]->Get(Form("IMnpim_IMnpip_dE_wK0_wSid_n_Sm_bin%d",ibin));
      IMnpim_IMnpip_dE_wSid_n_SpSm_bin[ibin][iq] = (TH2D*)fr[iq]->Get(Form("IMnpim_IMnpip_dE_wSid_n_SpSm_bin%d",ibin));
      OverlapCount[ibin][iq] = IMnpipi_wK0_wSid_n_SpSm[iq]->GetBinContent(ibin+1);
      if(OverlapCount[ibin][iq]<0.0) OverlapCount[ibin][iq]=0.0;
    }
    for(int iwbin=0;iwbin<nwbin;iwbin++){
      IMnpim_IMnpip_dE_wK0_wSid_n_Sp_wbin[iwbin][iq] = (TH2D*)fr[iq]->Get(Form("IMnpim_IMnpip_dE_wK0_wSid_n_Sp_wbin%d",iwbin));
      IMnpim_IMnpip_dE_wK0_wSid_n_Sm_wbin[iwbin][iq] = (TH2D*)fr[iq]->Get(Form("IMnpim_IMnpip_dE_wK0_wSid_n_Sm_wbin%d",iwbin));
    }
  }//iq

  //You have only qlo(=0) and qhi(1) decomposition results so far.
  TH2D* IMnpim_IMnpip_K0inter[2];
  TH2D* IMnpim_IMnpip_K0inter_sysup[2];
  TH2D* IMnpim_IMnpip_K0inter_sysdown[2];
  TGraphErrors *gr_SpONnpim_fin_pol1[2];
  TGraphErrors *gr_SmONnpip_fin_pol1[2];
  TGraphErrors *gr_Spratio_SporK0_sys[3][2];
  TGraphErrors *gr_Smratio_SmorK0_sys[3][2];

  //only qlo(=0) and qhi(1) decomposition results so far
  for(int iq=0;iq<2;iq++){
    IMnpim_IMnpip_K0inter[iq] = (TH2D*)fdeco[iq]->Get("h2K0inter_3fine");
    IMnpim_IMnpip_K0inter_sysup[iq] = (TH2D*)fdeco[iq]->Get("h2K0inter_3fine_sysup");
    IMnpim_IMnpip_K0inter_sysdown[iq] = (TH2D*)fdeco[iq]->Get("h2K0inter_3fine_sysdown");
    gr_SpONnpim_fin_pol1[iq] = (TGraphErrors*)fdeco[iq]->Get("gr_SpONnpim_fin_pol1_final");
    gr_SmONnpip_fin_pol1[iq] = (TGraphErrors*)fdeco[iq]->Get("gr_SmONnpip_fin_pol1_final");
    for(int isys=0;isys<3;isys++){
      std::cout << "iq: " << iq << " sys "  << isys-1 << std::endl;
      gr_Spratio_SporK0_sys[isys][iq] = (TGraphErrors*)fdeco[iq]->Get(Form("gr_Spratio_SporK0_sys%d",isys-1));
      gr_Spratio_SporK0_sys[isys][iq]->Print();
      gr_Smratio_SmorK0_sys[isys][iq] = (TGraphErrors*)fdeco[iq]->Get(Form("gr_Smratio_SmorK0_sys%d",isys-1));
      gr_Smratio_SmorK0_sys[isys][iq]->Print();
    }
  }
  
  std::cout << __LINE__ << std::endl;
  //1st index 0 : all
  //          1 : qlo
  //          2 : qhigh
  //          3 : theta_n small
  //0 sysdown
  //1 normal
  //2 sysup
  TH2D* q_IMnpipi_K0orSp_ToSp[4][3];
  TH2D* q_IMnpipi_K0orSp_ToK0[4][3];
  TH1D* IMnpipi_K0orSp_ToSp[4][3];
  TH1D* IMnpipi_K0orSp_ToK0[4][3];
  
  //nwbin:0 no overlap
  //nwbin:1 K0 overlap
  //nwbin:2 K0 overlap and Sigma- overlap
  const double wbinlow[nwbin] = {1.0,1.40,1.52};
  const double wbinhigh[nwbin] = {1.40,1.52,2.00};
  
  //q-low,q-hi,ntheta_small
  for(int iq=0;iq<4;iq++){
    for(int isys=0;isys<3;isys++){
      q_IMnpipi_K0orSp_ToSp[iq][isys] = (TH2D*)q_IMnpipi_wK0_wSid_n_Sp[iq]->Clone(Form("q_IMnpipi_K0orSp_ToSp%d_sys%d",iq,isys-1));
      q_IMnpipi_K0orSp_ToK0[iq][isys] = (TH2D*)q_IMnpipi_wK0_wSid_n_Sp[iq]->Clone(Form("q_IMnpipi_K0orSp_ToK0%d_sys%d",iq,isys-1));
      q_IMnpipi_K0orSp_ToSp[iq][isys]->Reset();
      q_IMnpipi_K0orSp_ToK0[iq][isys]->Reset();
      //std::cout << "bin low "  << spbinlow[iwbin-1][iq] << std::endl;
      //std::cout << "bin high " <<  spbinhi[iwbin-1][iq] << std::endl;

      double NevtWide=0.0;
      const int q350bin = q_IMnpipi_wK0_wSid_n_Sp[iq]->GetYaxis()->FindBin(0.35);
      for(int iwbin=1;iwbin<nwbin;iwbin++){
        int wbinl = q_IMnpipi_K0orSp_ToSp[iq][isys]->GetXaxis()->FindBin(wbinlow[iwbin]);
        int wbinh = q_IMnpipi_K0orSp_ToSp[iq][isys]->GetXaxis()->FindBin(wbinhigh[iwbin]);
        for(int ix=wbinl;ix<wbinh;ix++){
          for(int iqbin=0;iqbin<q_IMnpipi_wK0_wSid_n_Sp[iq]->GetNbinsY();iqbin++){
            int qlowhigh = 0;
            if(q350bin <= iqbin ) qlowhigh=1;
            double nevt = q_IMnpipi_wK0_wSid_n_Sp[iq]->GetBinContent(ix,iqbin);
            double err = q_IMnpipi_wK0_wSid_n_Sp[iq]->GetBinError(ix,iqbin);
            double ToSp =0.0;
            double ToK0 =0.0;
            double ToSperr =0.0;
            double ToK0err =0.0;
            double ratioToSp = 0.0;
            if(iwbin==1) ratioToSp = gr_Spratio_SporK0_sys[qlowhigh][isys]->Eval(1)  ;
            else if(iwbin==2) ratioToSp = gr_Spratio_SporK0_sys[qlowhigh][isys]->Eval(2)  ;
            double ratioToK0 = 1.0-ratioToSp;
            //if(isys==1 && iq==1 && qlowhigh==1){
            //}

            ToSp = nevt*ratioToSp;
            ToK0 = nevt*ratioToK0;
            ToSperr = err*ratioToSp;
            ToK0err = err*ratioToK0;
            q_IMnpipi_K0orSp_ToSp[iq][isys]->SetBinContent(ix,iqbin,ToSp);
            q_IMnpipi_K0orSp_ToSp[iq][isys]->SetBinError(ix,iqbin,ToSperr);
            q_IMnpipi_K0orSp_ToK0[iq][isys]->SetBinContent(ix,iqbin,ToK0);
            q_IMnpipi_K0orSp_ToK0[iq][isys]->SetBinError(ix,iqbin,ToK0err);
          }//iqbin
        }//ix
      }//iwbin
      IMnpipi_K0orSp_ToSp[iq][isys] = (TH1D*)q_IMnpipi_K0orSp_ToSp[iq][isys]->ProjectionX(Form("IMnpipi_K0orSp_ToSp%d_sys%d",iq,isys-1));
      IMnpipi_K0orSp_ToK0[iq][isys] = (TH1D*)q_IMnpipi_K0orSp_ToK0[iq][isys]->ProjectionX(Form("IMnpipi_K0orSp_ToK0%d_sys%d",iq,isys-1));
    }//isys
  }//iq

  //Index 0 : qlow
  //Index 1 : qhigh
  //Index 2 : theta_n small
  TH2D* q_IMnpipi_K0orSm_ToSm[4][3];
  TH2D* q_IMnpipi_K0orSm_ToK0[4][3];
  TH1D* IMnpipi_K0orSm_ToSm[4][3];
  TH1D* IMnpipi_K0orSm_ToK0[4][3];

  for(int iq=0;iq<4;iq++){
    for(int isys=0;isys<3;isys++){
      q_IMnpipi_K0orSm_ToSm[iq][isys] = (TH2D*)q_IMnpipi_wK0_wSid_n_Sm[iq]->Clone(Form("q_IMnpipi_K0orSm_ToSm%d_sys%d",iq,isys-1));
      q_IMnpipi_K0orSm_ToK0[iq][isys] = (TH2D*)q_IMnpipi_wK0_wSid_n_Sm[iq]->Clone(Form("q_IMnpipi_K0orSm_ToK0%d_sys%d",iq,isys-1));
      q_IMnpipi_K0orSm_ToSm[iq][isys]->Reset();
      q_IMnpipi_K0orSm_ToK0[iq][isys]->Reset();
      for(int iwbin=1;iwbin<nwbin;iwbin++){
        int wbinl = q_IMnpipi_K0orSm_ToSm[iq][isys]->GetXaxis()->FindBin(wbinlow[iwbin]);
        int wbinh = q_IMnpipi_K0orSm_ToSm[iq][isys]->GetXaxis()->FindBin(wbinhigh[iwbin]);
        std::cout << wbinl << "  " << wbinh << std::endl;
        double NevtWide=0.0;
        const int q350bin = q_IMnpipi_wK0_wSid_n_Sm[iq]->GetYaxis()->FindBin(0.35);
        for(int ix=wbinl;ix<wbinh;ix++){
          for(int iqbin=0;iqbin<q_IMnpipi_wK0_wSid_n_Sm[iq]->GetNbinsY();iqbin++){
            double nevt = q_IMnpipi_wK0_wSid_n_Sm[iq]->GetBinContent(ix,iqbin);
            double nerr = q_IMnpipi_wK0_wSid_n_Sm[iq]->GetBinError(ix,iqbin);
            int qlowhigh = 0;
            if(q350bin <= iqbin ) qlowhigh=1;
            double ToSm = 0.0;
            double ToK0 = 0.0;
            double ToSmerr = 0.0;
            double ToK0err = 0.0;
            double ratioToSm = 0; 
            if(iwbin==1) ratioToSm = gr_Smratio_SmorK0_sys[qlowhigh][isys]->Eval(1);
            else if(iwbin==2) ratioToSm = gr_Smratio_SmorK0_sys[qlowhigh][isys]->Eval(2);
            //if(isys==1 && iq==1 && qlowhigh==1){
            //  std::cout << "iwbin: "  << iwbin <<  " ratio "  << ratioToSm << std::endl;
            //}
            double ratioToK0 = 1.0-ratioToSm;
            ToSm = nevt*ratioToSm;
            ToK0 = nevt*ratioToK0;
            ToSmerr = nerr*ratioToSm;
            ToK0err = nerr*ratioToK0;
              
            q_IMnpipi_K0orSm_ToSm[iq][isys]->SetBinContent(ix,iqbin,ToSm);
            q_IMnpipi_K0orSm_ToSm[iq][isys]->SetBinError(ix,iqbin,ToSmerr);
            q_IMnpipi_K0orSm_ToK0[iq][isys]->SetBinContent(ix,iqbin,ToK0);
            q_IMnpipi_K0orSm_ToK0[iq][isys]->SetBinError(ix,iqbin,ToK0err);
          }//iqbin
        }//ix
      }//iwbin
      IMnpipi_K0orSm_ToSm[iq][isys] = (TH1D*)q_IMnpipi_K0orSm_ToSm[iq][isys]->ProjectionX(Form("IMnpipi_K0orSm_ToSm%d_sys%d",iq,isys-1));
      IMnpipi_K0orSm_ToK0[iq][isys] = (TH1D*)q_IMnpipi_K0orSm_ToK0[iq][isys]->ProjectionX(Form("IMnpipi_K0orSm_ToK0%d_sys%d",iq,isys-1));
    }//isys
  }//iq


  
  TCanvas *cIMnpipi_wK0_Sp_afterDeco[4];
  for(int iq=0;iq<4;iq++){
    cIMnpipi_wK0_Sp_afterDeco[iq] = new TCanvas(Form("cIMnpipi_wK0_Sp_Deco_%s",cqcut[iq]),Form("cIMnpipi_wK0_Sp_Deco_%s",cqcut[iq]),1000,800);
    IMnpipi_wK0_wSid_n_Sp[iq]->SetTitle(Form("IM(n#pi^{+}#pi^{-}) #bar{K}^{0}n or #Sigma^{+}#pi^{-} %s",cqcut[iq]));
    IMnpipi_wK0_wSid_n_Sp[iq]->Draw("E");
    IMnpipi_K0orSp_ToSp[iq][1]->SetLineColor(3);
    IMnpipi_K0orSp_ToSp[iq][1]->SetMarkerColor(3);
    IMnpipi_K0orSp_ToSp[iq][1]->Draw("Esame");
    IMnpipi_K0orSp_ToK0[iq][1]->SetLineColor(4);
    IMnpipi_K0orSp_ToK0[iq][1]->SetMarkerColor(4);
    IMnpipi_K0orSp_ToK0[iq][1]->Draw("Esame");

    TLegend *l = new TLegend(0.6,0.7,0.9,0.9);
    l->AddEntry(IMnpipi_wK0_wSid_n_Sp[iq],"overlap","l");
    l->AddEntry(IMnpipi_K0orSp_ToSp[iq][1],"ToSp","l");
    l->AddEntry(IMnpipi_K0orSp_ToK0[iq][1],"ToK0","l");
    l->Draw();
  }

  TCanvas *cIMnpipi_wK0_Sm_afterDeco[4];
  for(int iq=0;iq<4;iq++){
    cIMnpipi_wK0_Sm_afterDeco[iq] = new TCanvas(Form("cIMnpipi_wK0_Sm_Deco_%s",cqcut[iq]),Form("cIMnpipi_wK0_Sm_Deco_%s",cqcut[iq]),1000,800);
    IMnpipi_wK0_wSid_n_Sm[iq]->SetTitle(Form("IM(n#pi^{+}#pi^{-}) #bar{K}^{0}n or #Sigma^{-}#pi^{+} %s",cqcut[iq]));
    IMnpipi_wK0_wSid_n_Sm[iq]->Draw("E");
    IMnpipi_K0orSm_ToSm[iq][1]->SetLineColor(3);
    IMnpipi_K0orSm_ToSm[iq][1]->SetMarkerColor(3);
    IMnpipi_K0orSm_ToSm[iq][1]->Draw("Esame");
    IMnpipi_K0orSm_ToK0[iq][1]->SetLineColor(4);
    IMnpipi_K0orSm_ToK0[iq][1]->SetMarkerColor(4);
    IMnpipi_K0orSm_ToK0[iq][1]->Draw("Esame");
    TLegend *l = new TLegend(0.6,0.7,0.9,0.9);
    l->AddEntry(IMnpipi_wK0_wSid_n_Sm[iq],"overlap","l");
    l->AddEntry(IMnpipi_K0orSm_ToSm[iq][1],"ToSm","l");
    l->AddEntry(IMnpipi_K0orSm_ToK0[iq][1],"ToK0","l");
    l->Draw();
  }
  

  TH2D* q_IMnpipi_SporSm_ToSp[4][3];//iqcut,isys
  TH2D* q_IMnpipi_SporSm_ToSm[4][3];//iqcut,isys
  TH1D* IMnpipi_SporSm_ToSp[4][3];//iqcut,isys
  TH1D* IMnpipi_SporSm_ToSm[4][3];//iqcut,isys
  double nSp_SporSm[2][3];//q-region (low or hi), no wbin, because SporSm events are localized 
  double nSm_SporSm[2][3];//q-region (low or hi), no wbin
  for(int iqlowhigh=0;iqlowhigh<2;iqlowhigh++){
    double nSporSm = IMnpim_IMnpip_dE_wSid_n_SpSm[iqlowhigh+1]->Integral();
    double nSp = 0.0;
    double nSm = 0.0;
    if(nSporSm>0.0){
      for(int ix=0;ix<IMnpim_IMnpip_dE_wSid_n_SpSm[iqlowhigh+1]->GetNbinsX();ix++){
        for(int iy=0;iy<IMnpim_IMnpip_dE_wSid_n_SpSm[iqlowhigh+1]->GetNbinsY();iy++){
          double cont = IMnpim_IMnpip_dE_wSid_n_SpSm[iqlowhigh+1]->GetBinContent(ix,iy);
          if(cont>0.0){
            double bincent_npip = IMnpim_IMnpip_dE_wSid_n_SpSm[iqlowhigh+1]->GetXaxis()->GetBinCenter(ix);
            double bincent_npim = IMnpim_IMnpip_dE_wSid_n_SpSm[iqlowhigh+1]->GetYaxis()->GetBinCenter(iy);
            double nSp_bin = gr_SpONnpim_fin_pol1[iqlowhigh]->Eval(bincent_npim);
            double nSm_bin = gr_SmONnpip_fin_pol1[iqlowhigh]->Eval(bincent_npip);
            nSp += nSp_bin;//center
            nSm += nSm_bin;//center
          }//cont>0
        }//iy
      }//ix
    }//if nSporSm
    double nSp_bin_err = gr_SpONnpim_fin_pol1[iqlowhigh]->GetErrorY(29);
    double nSm_bin_err = gr_SmONnpip_fin_pol1[iqlowhigh]->GetErrorY(29);
    nSp_SporSm[iqlowhigh][0] = nSp-nSp_bin_err;//down Sp num
    nSm_SporSm[iqlowhigh][0] = nSm+nSm_bin_err;//up Sm num
    nSp_SporSm[iqlowhigh][1] = nSp;
    nSm_SporSm[iqlowhigh][1] = nSm;
    nSp_SporSm[iqlowhigh][2] = nSp+nSp_bin_err;//up Sp num
    nSm_SporSm[iqlowhigh][2] = nSm-nSm_bin_err;//down Sm num
  
  }

  for(int isys=0;isys<3;isys++){
    for(int iq=0;iq<4;iq++){
      q_IMnpipi_SporSm_ToSp[iq][isys] = (TH2D*)q_IMnpipi_wSid_n_SpSm[iq]->Clone(Form("q_IMnpipi_SporSm_ToSp%d_sys%d",iq,isys-1));
      q_IMnpipi_SporSm_ToSm[iq][isys] = (TH2D*)q_IMnpipi_wSid_n_SpSm[iq]->Clone(Form("q_IMnpipi_SporSm_ToSm%d_sys%d",iq,isys-1));
      for(int ibin=0;ibin<q_IMnpipi_SporSm_ToSp[iq][isys]->GetNbinsX();ibin++){
        for(int iqbin=0;iqbin<q_IMnpipi_SporSm_ToSp[iq][isys]->GetNbinsY();iqbin++){
          double evt = q_IMnpipi_SporSm_ToSp[iq][isys]->GetBinContent(ibin,iqbin);
          double err = q_IMnpipi_SporSm_ToSp[iq][isys]->GetBinError(ibin,iqbin);
          const int q350bin = q_IMnpipi_SporSm_ToSp[iq][isys]->GetYaxis()->FindBin(0.35);
          int qlowhigh = 0;
          if(q350bin <= iqbin ) qlowhigh=1;
          if((nSp_SporSm[qlowhigh][isys]+nSm_SporSm[qlowhigh][isys])>0.0){
            double evtToSp = evt*nSp_SporSm[qlowhigh][isys]/(nSp_SporSm[qlowhigh][isys]+nSm_SporSm[qlowhigh][isys]);
            double evtToSm = evt*nSm_SporSm[qlowhigh][isys]/(nSp_SporSm[qlowhigh][isys]+nSm_SporSm[qlowhigh][isys]);
            double evtToSperr = err*nSp_SporSm[qlowhigh][isys]/(nSp_SporSm[qlowhigh][isys]+nSm_SporSm[qlowhigh][isys]);
            double evtToSmerr = err*nSm_SporSm[qlowhigh][isys]/(nSp_SporSm[qlowhigh][isys]+nSm_SporSm[qlowhigh][isys]);
            q_IMnpipi_SporSm_ToSp[iq][isys]->SetBinContent(ibin,iqbin,evtToSp);
            q_IMnpipi_SporSm_ToSp[iq][isys]->SetBinError(ibin,iqbin,evtToSperr);
            q_IMnpipi_SporSm_ToSm[iq][isys]->SetBinContent(ibin,iqbin,evtToSm);
            q_IMnpipi_SporSm_ToSm[iq][isys]->SetBinError(ibin,iqbin,evtToSmerr);
            if((nSp_SporSm[qlowhigh][isys]+nSm_SporSm[qlowhigh][isys])<0.00001){
              std::cout << "nSp + nSm " << nSp_SporSm[qlowhigh][isys]+nSm_SporSm[qlowhigh][isys] << std::endl;
              std::cout << "evt " << evt << std::endl;
              std::cout << "err " << err << std::endl;
            }
          }//nSp+nSm>0
        }//iqbin
      }//ibin
    IMnpipi_SporSm_ToSp[iq][isys] = (TH1D*)q_IMnpipi_SporSm_ToSp[iq][isys]->ProjectionX(Form("IMnpipi_SporSm_ToSp%d_sys%d",iq,isys-1));
    IMnpipi_SporSm_ToSm[iq][isys] = (TH1D*)q_IMnpipi_SporSm_ToSm[iq][isys]->ProjectionX(Form("IMnpipi_SporSm_ToSm%d_sys%d",iq,isys-1));
    }//iq
  }//isys
  
  TCanvas *cIMnpipi_SpSm_decosys[4];
  for(int iq=0;iq<4;iq++){
    cIMnpipi_SpSm_decosys[iq] = new TCanvas(Form("cIMnpipi_SpSm_decosys_%d",iq),Form("cIMnpipi_SpSm_decosys_%d",iq),1000,800);
    IMnpipi_SporSm_ToSp[iq][1]->Draw();
    IMnpipi_SporSm_ToSp[iq][0]->SetLineColor(2);
    IMnpipi_SporSm_ToSp[iq][0]->Draw("same");
    IMnpipi_SporSm_ToSp[iq][2]->SetLineColor(3);
    IMnpipi_SporSm_ToSp[iq][2]->Draw("same");
  }


  std::cout << __LINE__ << std::endl;
  TCanvas *cIMnpipi_SpSm_afterDeco[4];
  for(int iq=0;iq<4;iq++){
    cIMnpipi_SpSm_afterDeco[iq] = new TCanvas(Form("cIMnpipi_SpSm_Deco_%s",cqcut[iq+1]),Form("cIMnpipi_SpSm_Deco_%s",cqcut[iq]),1000,800);
    IMnpipi_wSid_n_SpSm[iq]->SetTitle(Form("IM(n#pi^{+}#pi^{-}) #Sigma^{+}#pi^{-} or #Sigma^{-}#pi^{+} %s",cqcut[iq]));
    IMnpipi_wSid_n_SpSm[iq]->SetYTitle("Counts/ 15 MeV");
    IMnpipi_wSid_n_SpSm[iq]->GetYaxis()->CenterTitle();
    IMnpipi_wSid_n_SpSm[iq]->Draw("E");
    IMnpipi_SporSm_ToSp[iq][1]->SetLineColor(4);
    IMnpipi_SporSm_ToSp[iq][1]->SetMarkerColor(4);
    IMnpipi_SporSm_ToSp[iq][1]->Draw("Esame");
    IMnpipi_SporSm_ToSm[iq][1]->SetLineColor(5);
    IMnpipi_SporSm_ToSm[iq][1]->SetMarkerColor(5);
    IMnpipi_SporSm_ToSm[iq][1]->Draw("Esame");
    TLegend *l = new TLegend(0.6,0.7,0.9,0.9);
    l->AddEntry(IMnpipi_wSid_n_SpSm[iq+1],"overlap","l");
    l->AddEntry(IMnpipi_SporSm_ToSp[iq][1],"ToSp","l");
    l->AddEntry(IMnpipi_SporSm_ToSm[iq][1],"ToSm","l");
    l->Draw();
  }
  std::cout << __LINE__ << std::endl;

  //treatment of K0 & Sigma+ & Sigma- overlap
  TH2D* q_IMnpipi_K0SpSm_ToK0[4][3];
  TH2D* q_IMnpipi_K0SpSm_ToSp[4][3];
  TH2D* q_IMnpipi_K0SpSm_ToSm[4][3];
  TH1D* IMnpipi_K0SpSm_ToK0[4][3];
  TH1D* IMnpipi_K0SpSm_ToSp[4][3];
  TH1D* IMnpipi_K0SpSm_ToSm[4][3];
  double OverlapToSp[2][3];// qlow-high, isys
  double OverlapToSm[2][3];// qlow-high, isys
  double OverlapToK0[2][3];// qlow-high, isys
  
  for(int iqlowhi=0;iqlowhi<2;iqlowhi++){
    for(int isys=0;isys<3;isys++){
      int wbinl = IMnpipi_SporSm_ToSp[iqlowhi+1][isys]->GetXaxis()->FindBin(wbinlow[1]);//1.40
      int wbinh = IMnpipi_K0orSp_ToSp[iqlowhi+1][isys]->GetXaxis()->FindBin(wbinhigh[1]);//1.52
      OverlapToSp[iqlowhi][isys] = IMnpipi_SporSm_ToSp[iqlowhi+1][isys]->Integral(wbinl,wbinh);
      OverlapToSp[iqlowhi][isys] += IMnpipi_K0orSp_ToSp[iqlowhi+1][isys]->Integral(wbinl,wbinh);

      OverlapToSm[iqlowhi][isys] = IMnpipi_SporSm_ToSm[iqlowhi+1][isys]->Integral(wbinl,wbinh);
      OverlapToSm[iqlowhi][isys] += IMnpipi_K0orSm_ToSm[iqlowhi+1][isys]->Integral(wbinl,wbinh);

      OverlapToK0[iqlowhi][isys] = IMnpipi_K0orSm_ToK0[iqlowhi+1][isys]->Integral(wbinl,wbinh);
      OverlapToK0[iqlowhi][isys] += IMnpipi_K0orSp_ToK0[iqlowhi+1][isys]->Integral(wbinl,wbinh);

      std::cout << " iq " << iqlowhi << " isys" << isys << std::endl;
      std::cout <<  "Sp " << OverlapToSp[iqlowhi][isys] << std::endl;
      std::cout <<  "Sm " << OverlapToSm[iqlowhi][isys] << std::endl;
      std::cout <<  "K0 " << OverlapToK0[iqlowhi][isys] << std::endl;
    }
  }
  
  for(int iq=0;iq<4;iq++){
    for(int isys=0;isys<3;isys++){
      q_IMnpipi_K0SpSm_ToK0[iq][isys] = (TH2D*)q_IMnpipi_wK0_wSid_n_SpSm[iq]->Clone(Form("q_IMnpipi_K0SpSm_ToK0%d",iq));
      q_IMnpipi_K0SpSm_ToSp[iq][isys] = (TH2D*)q_IMnpipi_wK0_wSid_n_SpSm[iq]->Clone(Form("q_IMnpipi_K0SpSm_ToSp%d",iq));
      q_IMnpipi_K0SpSm_ToSm[iq][isys] = (TH2D*)q_IMnpipi_wK0_wSid_n_SpSm[iq]->Clone(Form("q_IMnpipi_K0SpSm_ToSm%d",iq));
      q_IMnpipi_K0SpSm_ToK0[iq][isys]->Reset();
      q_IMnpipi_K0SpSm_ToSp[iq][isys]->Reset();
      q_IMnpipi_K0SpSm_ToSm[iq][isys]->Reset();
      int wbinl = q_IMnpipi_K0SpSm_ToK0[iq][isys]->GetXaxis()->FindBin(wbinlow[1]);//1.40
      int wbinh = q_IMnpipi_K0SpSm_ToK0[iq][isys]->GetXaxis()->FindBin(wbinhigh[1]);//1.52
      for(int ix=wbinl;ix<wbinh;ix++){
        for(int iy=0;iy<q_IMnpipi_wK0_wSid_n_SpSm[iq]->GetNbinsY();iy++){
          double cont = q_IMnpipi_wK0_wSid_n_SpSm[iq]->GetBinContent(ix,iy);
          double err = q_IMnpipi_wK0_wSid_n_SpSm[iq]->GetBinError(ix,iy);
          int qlowhi = 0;
          const int q350bin = q_IMnpipi_wK0_wSid_n_SpSm[iq]->GetYaxis()->FindBin(0.35);
          if(q350bin <= iy ) qlowhi=1;
          double ToK0 = cont*OverlapToK0[qlowhi][isys]/(OverlapToSp[qlowhi][isys]+OverlapToSm[qlowhi][isys]+OverlapToK0[qlowhi][isys]);
          double ToK0err = err*OverlapToK0[qlowhi][isys]/(OverlapToSp[qlowhi][isys]+OverlapToSm[qlowhi][isys]+OverlapToK0[qlowhi][isys]);
          q_IMnpipi_K0SpSm_ToK0[iq][isys]->SetBinContent(ix,iy,ToK0);
          q_IMnpipi_K0SpSm_ToK0[iq][isys]->SetBinError(ix,iy,ToK0err);
          
          double ToSp = cont*OverlapToSp[qlowhi][isys]/(OverlapToSp[qlowhi][isys]+OverlapToSm[qlowhi][isys]+OverlapToK0[qlowhi][isys]);
          double ToSperr = err*OverlapToSp[qlowhi][isys]/(OverlapToSp[qlowhi][isys]+OverlapToSm[qlowhi][isys]+OverlapToK0[qlowhi][isys]);
          q_IMnpipi_K0SpSm_ToSp[iq][isys]->SetBinContent(ix,iy,ToSp);
          q_IMnpipi_K0SpSm_ToSp[iq][isys]->SetBinError(ix,iy,ToSperr);
          
          double ToSm = cont*OverlapToSm[qlowhi][isys]/(OverlapToSp[qlowhi][isys]+OverlapToSm[qlowhi][isys]+OverlapToK0[qlowhi][isys]);
          double ToSmerr = err*OverlapToSm[qlowhi][isys]/(OverlapToSp[qlowhi][isys]+OverlapToSm[qlowhi][isys]+OverlapToK0[qlowhi][isys]);
          q_IMnpipi_K0SpSm_ToSm[iq][isys]->SetBinContent(ix,iy,ToSm);
          q_IMnpipi_K0SpSm_ToSm[iq][isys]->SetBinError(ix,iy,ToSmerr);
        }//iy
      }//ix
      IMnpipi_K0SpSm_ToK0[iq][isys] = (TH1D*)q_IMnpipi_K0SpSm_ToK0[iq][isys]->ProjectionX(Form("IMnpipi_K0SpSm_ToK0_%d_sys%d",iq,isys));
      IMnpipi_K0SpSm_ToSp[iq][isys] = (TH1D*)q_IMnpipi_K0SpSm_ToSp[iq][isys]->ProjectionX(Form("IMnpipi_K0SpSm_ToSp_%d_sys%d",iq,isys));
      IMnpipi_K0SpSm_ToSm[iq][isys] = (TH1D*)q_IMnpipi_K0SpSm_ToSm[iq][isys]->ProjectionX(Form("IMnpipi_K0SpSm_ToSm_%d_sys%d",iq,isys));
    }//isys
  }//iq

  TCanvas *cK0SpSm[4];
  for(int iq=0;iq<4;iq++){
    cK0SpSm[iq] = new TCanvas(Form("cK0SpSm%d",iq),Form("cK0SpSm%d",iq));
    IMnpipi_wK0_wSid_n_SpSm[iq]->Draw("E");
    IMnpipi_K0SpSm_ToK0[iq][1]->SetLineColor(4);
    IMnpipi_K0SpSm_ToK0[iq][1]->SetMarkerColor(4);
    IMnpipi_K0SpSm_ToK0[iq][1]->Draw("Esame");
    IMnpipi_K0SpSm_ToSp[iq][1]->SetLineColor(2);
    IMnpipi_K0SpSm_ToSp[iq][1]->SetMarkerColor(2);
    IMnpipi_K0SpSm_ToSp[iq][1]->Draw("Esame");
    IMnpipi_K0SpSm_ToSm[iq][1]->SetLineColor(3);
    IMnpipi_K0SpSm_ToSm[iq][1]->SetMarkerColor(3);
    IMnpipi_K0SpSm_ToSm[iq][1]->Draw("Esame");
    TLegend *lo = new TLegend(0.6,0.7,0.9,0.9);
    lo->AddEntry(IMnpipi_K0SpSm_ToK0[iq][1],"To K0 mode","l");
    lo->AddEntry(IMnpipi_K0SpSm_ToSp[iq][1],"To Sp mode","l");
    lo->AddEntry(IMnpipi_K0SpSm_ToSm[iq][1],"To Sm mode","l");
    lo->Draw();
  }


  TCanvas *cq_IMnpipi_Sp_afterDeco[4];
  TH2D* q_IMnpipi_Sp_sum[4][3];
  for(int iq=0;iq<4;iq++){
    for(int isys=0;isys<3;isys++){
      q_IMnpipi_Sp_sum[iq][isys] = (TH2D*)q_IMnpipi_wSid_n_Sp[iq]->Clone(Form("q_IMnpipi_Sp_sum%d_sys%d",iq,isys-1));
      q_IMnpipi_Sp_sum[iq][isys]->Add(q_IMnpipi_K0orSp_ToK0[iq][isys],-1.0);
      q_IMnpipi_Sp_sum[iq][isys]->Add(q_IMnpipi_SporSm_ToSm[iq][isys],-1.0);
      q_IMnpipi_Sp_sum[iq][isys]->Add(q_IMnpipi_K0SpSm_ToSp[iq][isys],1.0);//added over-subtraction component
    }//isys
    cq_IMnpipi_Sp_afterDeco[iq] = new TCanvas(Form("cq_IMnpipi_Sp_afterDeco%d",iq),Form("cq_IMnpipi_Sp_afterDeco%d",iq));
    cq_IMnpipi_Sp_afterDeco[iq]->cd();
    //q_IMnpipi_Sp_sum[iq][1]->SetMinimum(0);
    q_IMnpipi_Sp_sum[iq][1]->Draw("colz");
  }

  TCanvas *cq_IMnpipi_Sm_afterDeco[4];
  TH2D* q_IMnpipi_Sm_sum[4][3];
  for(int iq=0;iq<4;iq++){
    for(int isys=0;isys<3;isys++){
      q_IMnpipi_Sm_sum[iq][isys] = (TH2D*)q_IMnpipi_wSid_n_Sm[iq]->Clone(Form("q_IMnpipi_Sm_sum%d_sys%d",iq,isys-1));
      q_IMnpipi_Sm_sum[iq][isys]->Add(q_IMnpipi_K0orSm_ToK0[iq][isys],-1.0);
      q_IMnpipi_Sm_sum[iq][isys]->Add(q_IMnpipi_SporSm_ToSp[iq][isys],-1.0);
      q_IMnpipi_Sm_sum[iq][isys]->Add(q_IMnpipi_K0SpSm_ToSm[iq][isys],1.0);//added over-subtraction component
    }//isys
    cq_IMnpipi_Sm_afterDeco[iq] = new TCanvas(Form("cq_IMnpipi_Sm_afterDeco%d",iq),Form("cq_IMnpipi_Sm_afterDeco%d",iq));
    cq_IMnpipi_Sm_afterDeco[iq]->cd();
    //q_IMnpipi_Sm_sum[iq][1]->SetMinimum(0);
    q_IMnpipi_Sm_sum[iq][1]->Draw("colz");
  }

  std::cout << __LINE__ << std::endl;

  TCanvas *cq_IMnpipi_K0_afterDeco[4];
  TH2D* q_IMnpipi_K0_sum[4][3];
  for(int iq=0;iq<4;iq++){
    for(int isys=0;isys<3;isys++){
      q_IMnpipi_K0_sum[iq][isys] = (TH2D*)q_IMnpipi_wK0_n[iq]->Clone(Form("q_IMnpipi_K0_sum%d_sys%d",iq,isys-1));
      q_IMnpipi_K0_sum[iq][isys]->Add(q_IMnpipi_K0orSp_ToSp[iq][isys],-1.0);
      q_IMnpipi_K0_sum[iq][isys]->Add(q_IMnpipi_K0orSm_ToSm[iq][isys],-1.0);
      q_IMnpipi_K0_sum[iq][isys]->Add(q_IMnpipi_K0SpSm_ToK0[iq][isys],1.0);//added over-subtraction component
    }
    cq_IMnpipi_K0_afterDeco[iq] = new TCanvas(Form("cq_IMnpipi_K0_afterDeco%d",iq),Form("cq_IMnpipi_K0_afterDeco%d",iq));
    cq_IMnpipi_K0_afterDeco[iq]->cd();
    //q_IMnpipi_K0_sum[iq][1]->SetMinimum(0);
    q_IMnpipi_K0_sum[iq][1]->Draw("colz");
  }
  std::cout << __LINE__ << std::endl;



  TCanvas *cIMnpipi_Summary_Sp[4];
  TH1D* IMnpipi_Sp_noK0_noSm[4][3];//iq,isys
  for(int iq=0;iq<4;iq++){
    cIMnpipi_Summary_Sp[iq] = new TCanvas(Form("cIMnpipi_Summary_Sp_%s",cqcut[iq]),Form("cIMnpipi_Summary_Sp_%s",cqcut[iq]),1000,800);
    //all sum
    IMnpipi_wSid_n_Sp[iq]->SetYTitle("Counts/ 15 MeV");
    IMnpipi_wSid_n_Sp[iq]->GetYaxis()->CenterTitle();
    IMnpipi_wSid_n_Sp[iq]->Draw("E");
    for(int isys=0;isys<3;isys++){
      IMnpipi_Sp_noK0_noSm[iq][isys] = (TH1D*)IMnpipi_wSid_n_Sp[iq]->Clone(Form("IMnpipi_Sp_noK0_noSm%d_sys%d",iq,isys));
      IMnpipi_Sp_noK0_noSm[iq][isys]->Add(IMnpipi_K0orSp_ToK0[iq][isys],-1);
      IMnpipi_Sp_noK0_noSm[iq][isys]->Add(IMnpipi_SporSm_ToSm[iq][isys],-1);
      IMnpipi_Sp_noK0_noSm[iq][isys]->Add(IMnpipi_K0SpSm_ToSp[iq][isys],1);//add over-subtraction component
    }
    IMnpipi_Sp_noK0_noSm[iq][1]->SetLineColor(2);
    IMnpipi_Sp_noK0_noSm[iq][1]->SetMarkerColor(2);
    IMnpipi_Sp_noK0_noSm[iq][1]->Draw("Esame");
    TLegend *lSp = new TLegend(0.6,0.7,0.9,0.9);
    lSp->AddEntry(IMnpipi_wSid_n_Sp[iq],"Sigma+ like mode","l");
    lSp->AddEntry(IMnpipi_Sp_noK0_noSm[iq][1],"Sigma+ after deco. ","l");
    lSp->Draw();
  }
  std::cout << __LINE__ << std::endl;
  //Systematic graph
  TGraphAsymmErrors *gDecoErrorSp[4];
  for(int iq=0;iq<4;iq++){
    gDecoErrorSp[iq] = new TGraphAsymmErrors(IMnpipi_Sp_noK0_noSm[iq][1]);
    for(int ip=0;ip<(gDecoErrorSp[iq]->GetN());ip++){
      double valdown = IMnpipi_Sp_noK0_noSm[iq][0]->GetBinContent(ip+1);
      double valdef = IMnpipi_Sp_noK0_noSm[iq][1]->GetBinContent(ip+1);
      double valup = IMnpipi_Sp_noK0_noSm[iq][2]->GetBinContent(ip+1);

      double yeh = (valdown-valdef);
      double yel = (valup-valdef);
      //std::cout << ip << "  " << yeh  << " "  <<  yel << std::endl;
      if(yeh>yel){
        gDecoErrorSp[iq]->SetPointEYhigh(ip,fabs(yeh));
        gDecoErrorSp[iq]->SetPointEYlow(ip,fabs(yel));
      }else{
        gDecoErrorSp[iq]->SetPointEYhigh(ip,fabs(yel));
        gDecoErrorSp[iq]->SetPointEYlow(ip,fabs(yeh));
      }
    }
    cIMnpipi_Summary_Sp[iq]->cd();
    gDecoErrorSp[iq]->SetFillStyle(0);
    gDecoErrorSp[iq]->SetLineColor(3);
    gDecoErrorSp[iq]->Draw("5");
  }
  


  //summary plot with the systematic error of 
  TCanvas *cIMnpipi_Summary_Sp_re[4];
  for(int iq=0;iq<4;iq++){
    cIMnpipi_Summary_Sp_re[iq] = new TCanvas(Form("cIMnpipi_Summary_Sp_re%s",cqcut[iq]),Form("cIMnpipi_Summary_Sp_re%s",cqcut[iq]),1000,800);
    //all sum
    IMnpipi_Sp_noK0_noSm[iq][1]->SetYTitle("Counts/ 15 MeV");
    IMnpipi_Sp_noK0_noSm[iq][1]->GetYaxis()->CenterTitle();
    IMnpipi_Sp_noK0_noSm[iq][1]->Draw("E");
    gDecoErrorSp[iq]->Draw("5");
  }

  TCanvas *cIMnpipi_Summary_Sm[4];
  TH1D* IMnpipi_Sm_noK0_noSp[4][3];
  for(int iq=0;iq<4;iq++){
    cIMnpipi_Summary_Sm[iq] = new TCanvas(Form("cIMnpipi_Summary_Sm_%s",cqcut[iq]),Form("cIMnpipi_Summary_Sm_%s",cqcut[iq]),1000,800);
    //all sum
    IMnpipi_wSid_n_Sm[iq]->SetYTitle("Counts/ 15 MeV");
    IMnpipi_wSid_n_Sm[iq]->GetYaxis()->CenterTitle();
    IMnpipi_wSid_n_Sm[iq]->Draw("E");
    for(int isys=0;isys<3;isys++){
      IMnpipi_Sm_noK0_noSp[iq][isys] = (TH1D*)IMnpipi_wSid_n_Sm[iq]->Clone(Form("IMnpipi_Sm_noK0_noSp%d_sys%d",iq,isys));
      IMnpipi_Sm_noK0_noSp[iq][isys]->Add(IMnpipi_K0orSm_ToK0[iq][isys],-1);
      IMnpipi_Sm_noK0_noSp[iq][isys]->Add(IMnpipi_SporSm_ToSp[iq][isys],-1);
      IMnpipi_Sm_noK0_noSp[iq][isys]->Add(IMnpipi_K0SpSm_ToSm[iq][isys],1);//add over-subtraction component
    }
    IMnpipi_Sm_noK0_noSp[iq][1]->SetLineColor(2);
    IMnpipi_Sm_noK0_noSp[iq][1]->SetMarkerColor(2);
    IMnpipi_Sm_noK0_noSp[iq][1]->Draw("Esame");
    TLegend *lSm = new TLegend(0.6,0.7,0.9,0.9);
    lSm->AddEntry(IMnpipi_wSid_n_Sm[iq],"Sigma- like mode","l");
    lSm->AddEntry(IMnpipi_Sm_noK0_noSp[iq][1],"Sigma- after deco. ","l");
    lSm->Draw();
  }
  
  //Systematic graph
  TGraphAsymmErrors *gDecoErrorSm[4];
  for(int iq=0;iq<4;iq++){
    gDecoErrorSm[iq] = new TGraphAsymmErrors(IMnpipi_Sm_noK0_noSp[iq][1]);
    for(int ip=0;ip<(gDecoErrorSm[iq]->GetN());ip++){
      double valdown = IMnpipi_Sm_noK0_noSp[iq][0]->GetBinContent(ip+1);
      double valdef = IMnpipi_Sm_noK0_noSp[iq][1]->GetBinContent(ip+1);
      double valup = IMnpipi_Sm_noK0_noSp[iq][2]->GetBinContent(ip+1);

      double yeh = (valdown-valdef);
      double yel = (valup-valdef);
      //std::cout << ip << "  " << yeh  << " "  <<  yel << std::endl;
      if(yeh>yel){
        gDecoErrorSm[iq]->SetPointEYhigh(ip,fabs(yeh));
        gDecoErrorSm[iq]->SetPointEYlow(ip,fabs(yel));
      }else{
        gDecoErrorSm[iq]->SetPointEYhigh(ip,fabs(yel));
        gDecoErrorSm[iq]->SetPointEYlow(ip,fabs(yeh));
      }
    }
    cIMnpipi_Summary_Sm[iq]->cd();
    gDecoErrorSm[iq]->SetFillStyle(0);
    gDecoErrorSm[iq]->SetLineColor(3);
    gDecoErrorSm[iq]->Draw("5");
  }

  TCanvas *cIMnpipi_Summary_Sm_re[4];
  for(int iq=0;iq<4;iq++){
    cIMnpipi_Summary_Sm_re[iq] = new TCanvas(Form("cIMnpipi_Summary_Sm_re%s",cqcut[iq]),Form("cIMnpipi_Summary_Sm_re%s",cqcut[iq]),1000,800);
    //all sum
    IMnpipi_Sm_noK0_noSp[iq][1]->SetYTitle("Counts/ 15 MeV");
    IMnpipi_Sm_noK0_noSp[iq][1]->GetYaxis()->CenterTitle();;
    IMnpipi_Sm_noK0_noSp[iq][1]->Draw("E");
    gDecoErrorSm[iq]->Draw("5");
  }


  TCanvas *cIMnpipi_Summary_K0[4];
  TH1D* IMnpipi_K0_noSp_noSm[4][3];
  for(int iq=0;iq<4;iq++){
    cIMnpipi_Summary_K0[iq] = new TCanvas(Form("cIMnpipi_Summary_K0_%s",cqcut[iq]),Form("cIMnpipi_Summary_K0_%s",cqcut[iq]),1000,800);
    //all sum
    IMnpipi_wK0_n[iq]->SetYTitle("Counts/ 15 MeV");
    IMnpipi_wK0_n[iq]->GetYaxis()->CenterTitle();
    IMnpipi_wK0_n[iq]->Draw("E");
    for(int isys=0;isys<3;isys++){
      IMnpipi_K0_noSp_noSm[iq][isys] = (TH1D*)IMnpipi_wK0_n[iq]->Clone(Form("IMnpipi_K0_noSp_noSm%d_%d",iq,isys));
      IMnpipi_K0_noSp_noSm[iq][isys]->Add(IMnpipi_K0orSp_ToSp[iq][isys],-1);
      IMnpipi_K0_noSp_noSm[iq][isys]->Add(IMnpipi_K0orSm_ToSm[iq][isys],-1);
      IMnpipi_K0_noSp_noSm[iq][isys]->Add(IMnpipi_K0SpSm_ToK0[iq][isys],1);//add over-subtraction component
    }
    IMnpipi_K0_noSp_noSm[iq][1]->SetLineColor(2);
    IMnpipi_K0_noSp_noSm[iq][1]->SetMarkerColor(2);
    IMnpipi_K0_noSp_noSm[iq][1]->Draw("Esame");
    TLegend *lK0 = new TLegend(0.6,0.7,0.9,0.9);
    lK0->AddEntry(IMnpipi_wK0_n[iq],"K0nn like mode","l");
    lK0->AddEntry(IMnpipi_K0_noSp_noSm[iq][1],"K0nn after deco. ","l");
    lK0->Draw();
  }
  TGraphAsymmErrors *gDecoErrorK0[4];
  for(int iq=0;iq<4;iq++){
    gDecoErrorK0[iq] = new TGraphAsymmErrors(IMnpipi_K0_noSp_noSm[iq][1]);
    for(int ip=0;ip<(gDecoErrorK0[iq]->GetN());ip++){
      double valdown = IMnpipi_K0_noSp_noSm[iq][0]->GetBinContent(ip+1);
      double valdef = IMnpipi_K0_noSp_noSm[iq][1]->GetBinContent(ip+1);
      double valup = IMnpipi_K0_noSp_noSm[iq][2]->GetBinContent(ip+1);

      double yeh = (valdown-valdef);
      double yel = (valup-valdef);
      //std::cout << ip << "  " << yeh  << " "  <<  yel << std::endl;
      if(yeh>yel){
        gDecoErrorK0[iq]->SetPointEYhigh(ip,fabs(yeh));
        gDecoErrorK0[iq]->SetPointEYlow(ip,fabs(yel));
      }else{
        gDecoErrorK0[iq]->SetPointEYhigh(ip,fabs(yel));
        gDecoErrorK0[iq]->SetPointEYlow(ip,fabs(yeh));
      }
    }
    cIMnpipi_Summary_K0[iq]->cd();
    gDecoErrorK0[iq]->SetFillStyle(0);
    gDecoErrorK0[iq]->SetLineColor(3);
    gDecoErrorK0[iq]->Draw("5");
  }

  TCanvas *cIMnpipi_Summary_K0_re[4];
  for(int iq=0;iq<4;iq++){
    cIMnpipi_Summary_K0_re[iq] = new TCanvas(Form("cIMnpipi_Summary_K0_re%s",cqcut[iq]),Form("cIMnpipi_Summary_K0_re%s",cqcut[iq]),1000,800);
    //all sum
    IMnpipi_K0_noSp_noSm[iq][1]->SetYTitle("Counts/ 15 MeV");
    IMnpipi_K0_noSp_noSm[iq][1]->GetYaxis()->CenterTitle();;
    IMnpipi_K0_noSp_noSm[iq][1]->Draw("E");
    gDecoErrorK0[iq]->Draw("5");
  }

  //
  //Acceptance correction
  //

  std::cout << __LINE__ << std::endl;
  //TFile *facc = TFile::Open("../simpost/accmap.root");
  TFile *facc = TFile::Open(Form("../simpost/accmapv%d_%d_dE%d.root",versionSigma,versionK0,dEcut));
  TH2D *q_IMnpipi_Sp_accp= (TH2D*)facc->Get("q_IMnpipi_Sp_accp_0");
  TH2D *q_IMnpipi_Sm_accp= (TH2D*)facc->Get("q_IMnpipi_Sm_accp_0");
  TH2D *q_IMnpipi_K0_accp= (TH2D*)facc->Get("q_IMnpipi_K0_accp_0");
  TH2D *q_IMnpipi_Sp_accperr= (TH2D*)facc->Get("q_IMnpipi_Sp_accerr_0");
  TH2D *q_IMnpipi_Sm_accperr= (TH2D*)facc->Get("q_IMnpipi_Sm_accerr_0");
  TH2D *q_IMnpipi_K0_accperr= (TH2D*)facc->Get("q_IMnpipi_K0_accerr_0");
  
  TCanvas *cSpacc = new TCanvas("cSpacc","cSpacc",1600,800);
  cSpacc->Divide(2,1);
  cSpacc->cd(1);
  q_IMnpipi_Sp_accp->Draw("colz");
  cSpacc->cd(2);
  q_IMnpipi_Sp_accperr->Draw("colz");

  TCanvas *cSmacc = new TCanvas("cSmacc","cSmacc",1600,800);
  cSmacc->Divide(2,1);
  cSmacc->cd(1);
  q_IMnpipi_Sm_accp->Draw("colz");
  cSmacc->cd(2);
  q_IMnpipi_Sm_accperr->Draw("colz");

  TCanvas *cK0acc = new TCanvas("cK0acc","cK0acc",1600,800);
  cK0acc->Divide(2,1);
  cK0acc->cd(1);
  q_IMnpipi_K0_accp->Draw("colz");
  cK0acc->cd(2);
  q_IMnpipi_K0_accperr->Draw("colz");


  TH2D* q_IMnpipi_Sp_cs[4][3];//iq,isys
  TH2D* q_IMnpipi_Sm_cs[4][3];//iq,isys
  TH2D* q_IMnpipi_K0_cs[4][3];//iq,isys
  TH2D* q_IMnpipi_Sp_cserr[4][3];//iq,isys
  TH2D* q_IMnpipi_Sm_cserr[4][3];//iq,isys
  TH2D* q_IMnpipi_K0_cserr[4][3];//iq,isys
  const int qcut350 =q_IMnpipi_Sp_sum[0][1]->GetYaxis()->FindBin(0.35);
  const int qcut600 =q_IMnpipi_Sp_sum[0][1]->GetYaxis()->FindBin(0.60);
  for(int iq=0;iq<4;iq++){
    for(int isys=0;isys<3;isys++){
      q_IMnpipi_Sp_cs[iq][isys] = (TH2D*)q_IMnpipi_Sp_sum[iq][isys]->Clone(Form("q_IMnpipi_Sp_cs%d_sys%d",iq,isys-1));
      q_IMnpipi_Sm_cs[iq][isys] = (TH2D*)q_IMnpipi_Sm_sum[iq][isys]->Clone(Form("q_IMnpipi_Sm_cs%d_sys%d",iq,isys-1));
      q_IMnpipi_K0_cs[iq][isys] = (TH2D*)q_IMnpipi_K0_sum[iq][isys]->Clone(Form("q_IMnpipi_K0_cs%d_sys%d",iq,isys-1));
      q_IMnpipi_Sp_cserr[iq][isys] = (TH2D*)q_IMnpipi_Sp_sum[iq][isys]->Clone(Form("q_IMnpipi_Sp_cserr%d_sys%d",iq,isys-1));
      q_IMnpipi_Sm_cserr[iq][isys] = (TH2D*)q_IMnpipi_Sm_sum[iq][isys]->Clone(Form("q_IMnpipi_Sm_cserr%d_sys%d",iq,isys-1));
      q_IMnpipi_K0_cserr[iq][isys] = (TH2D*)q_IMnpipi_K0_sum[iq][isys]->Clone(Form("q_IMnpipi_K0_cserr%d_sys%d",iq,isys-1));
      q_IMnpipi_Sp_cs[iq][isys]->SetTitle(Form("q_IMnpipi_Sp_cs_%s_sys%d",cqcut[iq],isys-1));
      q_IMnpipi_Sm_cs[iq][isys]->SetTitle(Form("q_IMnpipi_Sm_cs_%s_sys%d",cqcut[iq],isys-1));
      q_IMnpipi_K0_cs[iq][isys]->SetTitle(Form("q_IMnpipi_K0_cs_%s_sys%d",cqcut[iq],isys-1));
      q_IMnpipi_Sp_cserr[iq][isys]->SetTitle(Form("q_IMnpipi_Sp_cserr_%s_sys%d",cqcut[iq],isys-1));
      q_IMnpipi_Sm_cserr[iq][isys]->SetTitle(Form("q_IMnpipi_Sm_cserr_%s_sys%d",cqcut[iq],isys-1));
      q_IMnpipi_K0_cserr[iq][isys]->SetTitle(Form("q_IMnpipi_K0_cserr_%s_sts%d",cqcut[iq],isys-1));

      for(int ix=0;ix<=q_IMnpipi_Sp_cs[iq][isys]->GetNbinsX();ix++){
        for(int iy=0;iy<=q_IMnpipi_Sp_cs[iq][isys]->GetNbinsY();iy++){
          double contSp =q_IMnpipi_Sp_cs[iq][isys]->GetBinContent(ix,iy);
          double contSperr =q_IMnpipi_Sp_cs[iq][isys]->GetBinError(ix,iy);
          double accSp =q_IMnpipi_Sp_accp->GetBinContent(ix,iy);
          double accerrSp =q_IMnpipi_Sp_accperr->GetBinContent(ix,iy);
          double csSp = 0.0; 
          double csSperr = 0.0; 
          double binwidthM = q_IMnpipi_Sp_cs[iq][isys]->ProjectionX()->GetBinWidth(1)*1000.0;
          double binwidthq = q_IMnpipi_Sp_cs[iq][isys]->ProjectionY()->GetBinWidth(1)*1000.0;

          if(accSp>0.0){
            csSp = contSp/accSp/binwidthM/binwidthq/trigScale/lumi ;
            csSperr = contSperr/accSp/binwidthM/binwidthq/trigScale/lumi;
          }
          double contSm =q_IMnpipi_Sm_cs[iq][isys]->GetBinContent(ix,iy);
          double contSmerr =q_IMnpipi_Sm_cs[iq][isys]->GetBinError(ix,iy);
          double accSm =q_IMnpipi_Sm_accp->GetBinContent(ix,iy);
          double accerrSm =q_IMnpipi_Sm_accperr->GetBinContent(ix,iy);
          double csSm = 0.0; 
          double csSmerr = 0.0; 
          if(accSm>0.0){
            csSm = contSm/accSm/binwidthM/binwidthq/trigScale/lumi;
            csSmerr = contSmerr/accSm/binwidthM/binwidthq/trigScale/lumi;
          }

          double contK0 =q_IMnpipi_K0_cs[iq][isys]->GetBinContent(ix,iy);
          double contK0err =q_IMnpipi_K0_cs[iq][isys]->GetBinError(ix,iy);
          double accK0 =q_IMnpipi_K0_accp->GetBinContent(ix,iy);
          double accerrK0 =q_IMnpipi_K0_accperr->GetBinContent(ix,iy);
          double csK0 = 0.0; 
          double csK0err = 0.0; 
          if(accK0>0.0){
            csK0 = contK0/accK0/binwidthM/binwidthq/trigScale/lumi;
            csK0err = contK0err/accK0/binwidthM/binwidthq/trigScale/lumi;
          }

          if(accerrSp<UncertCut){
            q_IMnpipi_Sp_cs[iq][isys]->SetBinContent(ix,iy,csSp);
            q_IMnpipi_Sp_cs[iq][isys]->SetBinError(ix,iy,csSperr);
            q_IMnpipi_Sp_cserr[iq][isys]->SetBinContent(ix,iy,csSperr/csSp);
          }else{
            q_IMnpipi_Sp_cs[iq][isys]->SetBinContent(ix,iy,0.);
            q_IMnpipi_Sp_cs[iq][isys]->SetBinError(ix,iy,0.);
            q_IMnpipi_Sp_cserr[iq][isys]->SetBinError(ix,iy,0.);
          }
        
          if(iy>qcut600){
            q_IMnpipi_Sp_cs[iq][isys]->SetBinContent(ix,iy,0.);
            q_IMnpipi_Sp_cs[iq][isys]->SetBinError(ix,iy,0.);
            q_IMnpipi_Sp_cserr[iq][isys]->SetBinContent(ix,iy,0.);
            q_IMnpipi_Sp_cserr[iq][isys]->SetBinError(ix,iy,0.);
          }
          if(accerrSm<UncertCut){
            q_IMnpipi_Sm_cs[iq][isys]->SetBinContent(ix,iy,csSm);
            q_IMnpipi_Sm_cs[iq][isys]->SetBinError(ix,iy,csSmerr);
            q_IMnpipi_Sm_cserr[iq][isys]->SetBinError(ix,iy,csSmerr/csSm);
          }else{
            q_IMnpipi_Sm_cs[iq][isys]->SetBinContent(ix,iy,0.);
            q_IMnpipi_Sm_cs[iq][isys]->SetBinError(ix,iy,0.);
          }
        
          if(iy>qcut600){
            q_IMnpipi_Sm_cs[iq][isys]->SetBinContent(ix,iy,0.);
            q_IMnpipi_Sm_cs[iq][isys]->SetBinError(ix,iy,0.);
            q_IMnpipi_Sm_cserr[iq][isys]->SetBinContent(ix,iy,0.);
            q_IMnpipi_Sm_cserr[iq][isys]->SetBinError(ix,iy,0.);
          }
        
          if(accerrK0<UncertCut){
            q_IMnpipi_K0_cs[iq][isys]->SetBinContent(ix,iy,csK0);
            q_IMnpipi_K0_cs[iq][isys]->SetBinError(ix,iy,csK0err);
            q_IMnpipi_K0_cserr[iq][isys]->SetBinError(ix,iy,csK0err/csK0);
          }else{
            q_IMnpipi_K0_cs[iq][isys]->SetBinContent(ix,iy,0.);
            q_IMnpipi_K0_cs[iq][isys]->SetBinError(ix,iy,0.);
          }
        
          if(iy>qcut600){
            q_IMnpipi_K0_cs[iq][isys]->SetBinContent(ix,iy,0.);
            q_IMnpipi_K0_cs[iq][isys]->SetBinError(ix,iy,0.);
            q_IMnpipi_K0_cserr[iq][isys]->SetBinContent(ix,iy,0.);
            q_IMnpipi_K0_cserr[iq][isys]->SetBinError(ix,iy,0.);
          }
        }//iy
      }//ix
    }//isys
  }//iq



  TCanvas *ccsSp[4];
  TCanvas *ccsSppro[4];
  TCanvas *ccsSm[4];
  TCanvas *ccsSmpro[4];
  TCanvas *ccsK0[4];
  TCanvas *ccsK0pro[4];
  
  //display each 1D CS
  TH1D* IMnpipi_Sp_cs_single[4][3];//iq,isys
  TH1D* IMnpipi_Sm_cs_single[4][3];//iq,isys
  TH1D* IMnpipi_K0_cs_single[4][3];//iq,isys
  for(int iq=0;iq<4;iq++){
    ccsSp[iq] = new TCanvas(Form("ccsSp%d",iq),Form("ccsSp%d",iq),1600,800);
    ccsSp[iq]->Divide(2,1);
    ccsSp[iq]->cd(1);
    q_IMnpipi_Sp_cs[iq][1]->Draw("colz");
    ccsSp[iq]->cd(2);
    q_IMnpipi_Sp_cserr[iq][1]->Draw("colz");
    
    double binwidthq = q_IMnpipi_Sp_cs[iq][1]->ProjectionY()->GetBinWidth(1)*1000.0;
    for(int isys=0;isys<3;isys++){
      IMnpipi_Sp_cs_single[iq][isys]= (TH1D*)q_IMnpipi_Sp_cs[iq][isys]->ProjectionX(Form("IMnpipi_Sp_cs_single%d_sys%d",iq,isys-1),1,qcut600);
      IMnpipi_Sp_cs_single[iq][isys]->SetYTitle("d#sigma/dM [#mu b (MeV/c^{2})]");
      IMnpipi_Sp_cs_single[iq][isys]->GetYaxis()->CenterTitle();
      IMnpipi_Sp_cs_single[iq][isys]->Scale(binwidthq);
    }
    ccsSppro[iq] = new TCanvas(Form("ccsSppro%d",iq),Form("ccsSppro%d",iq),800,800);

    //IMnpipi_Sp_cs_single[iq][1]->SetMinimum(0);
    IMnpipi_Sp_cs_single[iq][1]->Draw("E");

    ccsSm[iq] = new TCanvas(Form("ccsSm%d",iq),Form("ccsSm%d",iq),1600,800);
    ccsSm[iq]->Divide(2,1);
    ccsSm[iq]->cd(1);
    q_IMnpipi_Sm_cs[iq][1]->Draw("colz");
    ccsSm[iq]->cd(2);
    q_IMnpipi_Sm_cserr[iq][1]->Draw("colz");
    
    for(int isys=0;isys<3;isys++){
      IMnpipi_Sm_cs_single[iq][isys]= (TH1D*)q_IMnpipi_Sm_cs[iq][isys]->ProjectionX(Form("IMnpipi_Sm_cs_single%d_sys%d",iq,isys-1),1,qcut600);
      IMnpipi_Sm_cs_single[iq][isys]->SetYTitle("d#sigma/dM [#mu b (MeV/c^{2})]");
      IMnpipi_Sm_cs_single[iq][isys]->GetYaxis()->CenterTitle();
      IMnpipi_Sm_cs_single[iq][isys]->Scale(binwidthq);
    }

    ccsSmpro[iq] = new TCanvas(Form("ccsSmpro%d",iq),Form("ccsSmpro%d",iq),800,800);
    IMnpipi_Sm_cs_single[iq][1]->SetMinimum(0);
    IMnpipi_Sm_cs_single[iq][1]->Draw("E");

    ccsK0[iq] = new TCanvas(Form("ccsK0%d",iq),Form("ccsK0%d",iq),1600,800);
    ccsK0[iq]->Divide(2,1);
    ccsK0[iq]->cd(1);
    q_IMnpipi_K0_cs[iq][1]->Draw("colz");
    ccsK0[iq]->cd(2);
    q_IMnpipi_K0_cserr[iq][1]->Draw("colz");

    for(int isys=0;isys<3;isys++){
      IMnpipi_K0_cs_single[iq][isys]= (TH1D*)q_IMnpipi_K0_cs[iq][isys]->ProjectionX(Form("IMnpipi_K0_cs_single%d_sys%d",iq,isys-1),1,qcut600);
      IMnpipi_K0_cs_single[iq][isys]->SetYTitle("d#sigma/dM [#mu b (MeV/c^{2})]");
      IMnpipi_K0_cs_single[iq][isys]->GetYaxis()->CenterTitle();
      IMnpipi_K0_cs_single[iq][isys]->Scale(binwidthq);
    }
    ccsK0pro[iq] = new TCanvas(Form("ccsK0pro%d",iq),Form("ccsK0pro%d",iq),800,800);
    IMnpipi_K0_cs_single[iq][1]->SetMinimum(0);
    IMnpipi_K0_cs_single[iq][1]->Draw("E");
  }
  
  //systematic error of decomposition
  TGraphAsymmErrors *gDecoErrorSp_CS[4];
  TGraphAsymmErrors *gDecoErrorSm_CS[4];
  TGraphAsymmErrors *gDecoErrorK0_CS[4];
  for(int iq=0;iq<4;iq++){
    gDecoErrorSp_CS[iq] = new TGraphAsymmErrors(IMnpipi_Sp_cs_single[iq][1]);
    gDecoErrorSm_CS[iq] = new TGraphAsymmErrors(IMnpipi_Sm_cs_single[iq][1]);
    gDecoErrorK0_CS[iq] = new TGraphAsymmErrors(IMnpipi_K0_cs_single[iq][1]);
    for(int ip=0;ip<( gDecoErrorSp_CS[iq]->GetN());ip++){
      double valdown = IMnpipi_Sp_cs_single[iq][0]->GetBinContent(ip+1);
      double valdef = IMnpipi_Sp_cs_single[iq][1]->GetBinContent(ip+1);
      double valup = IMnpipi_Sp_cs_single[iq][2]->GetBinContent(ip+1);

      double yeh = valdown-valdef;
      double yel = valup-valdef;
     
      if(yeh>yel){ 
        gDecoErrorSp_CS[iq]->SetPointEYhigh(ip,fabs(yeh));
        gDecoErrorSp_CS[iq]->SetPointEYlow(ip,fabs(yel));
      }else{
        gDecoErrorSp_CS[iq]->SetPointEYhigh(ip,fabs(yel));
        gDecoErrorSp_CS[iq]->SetPointEYlow(ip,fabs(yeh));
      }
    }
    for(int ip=0;ip<( gDecoErrorSm_CS[iq]->GetN());ip++){
      double valdown = IMnpipi_Sm_cs_single[iq][0]->GetBinContent(ip+1);
      double valdef = IMnpipi_Sm_cs_single[iq][1]->GetBinContent(ip+1);
      double valup = IMnpipi_Sm_cs_single[iq][2]->GetBinContent(ip+1);

      double yeh = valdown-valdef;
      double yel = valup-valdef;
     
      if(yeh>yel){ 
        gDecoErrorSm_CS[iq]->SetPointEYhigh(ip,fabs(yeh));
        gDecoErrorSm_CS[iq]->SetPointEYlow(ip,fabs(yel));
      }else{
        gDecoErrorSm_CS[iq]->SetPointEYhigh(ip,fabs(yel));
        gDecoErrorSm_CS[iq]->SetPointEYlow(ip,fabs(yeh));
      }
    }
    for(int ip=0;ip<( gDecoErrorK0_CS[iq]->GetN());ip++){
      double valdown = IMnpipi_K0_cs_single[iq][0]->GetBinContent(ip+1);
      double valdef = IMnpipi_K0_cs_single[iq][1]->GetBinContent(ip+1);
      double valup = IMnpipi_K0_cs_single[iq][2]->GetBinContent(ip+1);

      double yeh = valdown-valdef;
      double yel = valup-valdef;
     
      if(yeh>yel){ 
        gDecoErrorK0_CS[iq]->SetPointEYhigh(ip,fabs(yeh));
        gDecoErrorK0_CS[iq]->SetPointEYlow(ip,fabs(yel));
      }else{
        gDecoErrorK0_CS[iq]->SetPointEYhigh(ip,fabs(yel));
        gDecoErrorK0_CS[iq]->SetPointEYlow(ip,fabs(yeh));
      }
    }
    ccsSppro[iq]->cd();
    gDecoErrorSp_CS[iq]->SetFillStyle(0);
    gDecoErrorSp_CS[iq]->SetLineColor(3);
    gDecoErrorSp_CS[iq]->Draw("5");
    ccsSmpro[iq]->cd();
    gDecoErrorSm_CS[iq]->SetFillStyle(0);
    gDecoErrorSm_CS[iq]->SetLineColor(3);
    gDecoErrorSm_CS[iq]->Draw("5");
    ccsK0pro[iq]->cd();
    gDecoErrorK0_CS[iq]->SetFillStyle(0);
    gDecoErrorK0_CS[iq]->SetLineColor(3);
    gDecoErrorK0_CS[iq]->Draw("5");
  }



  //Sigma+ Sigma- charge sum
  TCanvas *csum[4];
  TH2D* q_IMnpipi_SpSmSum[4][3];//iq,isys
  TH1D* IMnpipi_SpSmSum[4][3];//iq,isys
  for(int iq=0;iq<4;iq++){
    for(int isys=0;isys<3;isys++){
      q_IMnpipi_SpSmSum[iq][isys] = (TH2D*)q_IMnpipi_Sp_cs[iq][isys]->Clone(Form("q_IMnpipi_SpSmSum%d_sys%d",iq,isys-1));
      q_IMnpipi_SpSmSum[iq][isys]->Add(q_IMnpipi_Sm_cs[iq][isys],1.0);
      q_IMnpipi_SpSmSum[iq][isys]->SetTitle(Form("q_IMnpipi_SpSmSum_%s_sys%d",cqcut[iq],isys-1));
      IMnpipi_SpSmSum[iq][isys] = (TH1D*)q_IMnpipi_SpSmSum[iq][isys]->ProjectionX(Form("IMnpipi_SpSmSum%d_sys%d",iq,isys-1),1,qcut600);
    }
    csum[iq]  = new TCanvas(Form("csum%d",iq),Form("csum%d",iq),1600,800);
    csum[iq]->Divide(2,1);
    csum[iq]->cd(1);
    q_IMnpipi_SpSmSum[iq][1]->SetMinimum(0);
    q_IMnpipi_SpSmSum[iq][1]->Draw("colz");
    csum[iq]->cd(2);
    IMnpipi_SpSmSum[iq][1]->SetTitle("#Sigma^{+} #Sigma^{-} charge Sum");
    IMnpipi_SpSmSum[iq][1]->SetMinimum(-0.005);
    IMnpipi_SpSmSum[iq][1]->Draw("E");
  }

  TCanvas *csum2D = new TCanvas("csum2D","csum2D");
  csum2D->cd();
  q_IMnpipi_SpSmSum[1][1]->Draw("col");
  q_IMnpipi_SpSmSum[2][1]->Draw("colsame");



  TCanvas *csub = new TCanvas("csub","csub",1600,800);
  csub->Divide(2,1);
  TH2D* q_IMnpipi_SpSmSub[4][3];
  for(int iq=0;iq<4;iq++){
    for(int isys=0;isys<3;isys++){
      q_IMnpipi_SpSmSub[iq][isys] = (TH2D*)q_IMnpipi_Sp_cs[iq][isys]->Clone(Form("q_IMnpipi_SpSmSub%d_sys%d",iq,isys-1));
      q_IMnpipi_SpSmSub[iq][isys]->Add(q_IMnpipi_Sm_cs[iq][isys],-1.0);
    }
  }
  csub->cd(1);
  q_IMnpipi_SpSmSub[1][1]->SetTitle("q_IMnpipi_SpSmSub");
  q_IMnpipi_SpSmSub[1][1]->Draw("colz");
  csub->cd(2);
  q_IMnpipi_SpSmSub[1][1]->ProjectionX("sum_px",1,qcut600)->Draw("E");
   
  TCanvas *ccomp[4];//iq,isys
  TH1D* IMnpipi_Sp_cs[4][3];
  TH1D* IMnpipi_Sm_cs[4][3];
  TH1D* IMnpipi_K0_cs[4][3];
  
  for(int iq=0;iq<4;iq++){
    for(int isys=0;isys<3;isys++){
      IMnpipi_Sp_cs[iq][isys] = (TH1D*)IMnpipi_Sp_cs_single[iq][isys]->Clone(Form("IMnpipi_Sp_cs%d_sys%d",iq,isys-1));
      IMnpipi_Sm_cs[iq][isys] = (TH1D*)IMnpipi_Sm_cs_single[iq][isys]->Clone(Form("IMnpipi_Sm_cs%d_sys%d",iq,isys-1));
      IMnpipi_K0_cs[iq][isys] = (TH1D*)IMnpipi_K0_cs_single[iq][isys]->Clone(Form("IMnpipi_K0_cs%d_sys%d",iq,isys-1));
      IMnpipi_Sp_cs[iq][isys]->SetLineColor(3);
      IMnpipi_Sp_cs[iq][isys]->SetMarkerColor(3);
      IMnpipi_Sm_cs[iq][isys]->SetLineColor(4);
      IMnpipi_Sm_cs[iq][isys]->SetMarkerColor(4);
      IMnpipi_K0_cs[iq][isys]->SetLineColor(2);
      IMnpipi_K0_cs[iq][isys]->SetMarkerColor(2);
    }
    ccomp[iq] = new TCanvas(Form("ccomp%d",iq),Form("ccomp%d",iq),1000,800);
    //IMnpipi_K0_cs[iq]->RebinX(2);
    //IMnpipi_Sp_cs[iq]->RebinX(2);
    //IMnpipi_Sm_cs[iq]->RebinX(2);
    //IMnpipi_K0_cs[iq]->Draw("E");
    //IMnpipi_Sp_cs[iq]->Draw("Esame");
    IMnpipi_Sp_cs[iq][1]->SetMarkerStyle(20);
    IMnpipi_Sp_cs[iq][1]->SetMinimum(-0.1);
    IMnpipi_Sp_cs[iq][1]->Draw("E");
    IMnpipi_Sp_cs[iq][1]->SetYTitle("d#sigma/dM [#mu b (MeV/c^{2})]");
    IMnpipi_Sp_cs[iq][1]->GetYaxis()->CenterTitle();
    IMnpipi_Sm_cs[iq][1]->Draw("Esame");
    //IMnpipi_K0_cs[iq]->Draw("Esame");
  
    TLegend *lcs = new TLegend(0.6,0.7,0.9,0.9);
    lcs->AddEntry(IMnpipi_Sp_cs[iq][1],"#Sigma^{+} mode","l");
    lcs->AddEntry(IMnpipi_Sm_cs[iq][1],"#Sigma^{-} mode","l");
    //lcs->AddEntry(IMnpipi_K0_cs[iq],"#bar{K}^{0} mode","l");
    lcs->Draw();
  }

  TIter nexthist(gDirectory->GetList());
  TH1F *h1 = NULL;
  TH1D *h1d = NULL;
  TH2F *h2 = NULL;
  TObject *obj = NULL;
  while( (obj = (TObject*)nexthist())!=NULL  ) {
    if(obj->InheritsFrom("TH1F")) {
      h1 = (TH1F*) obj;
      h1->GetXaxis()->CenterTitle();
      h1->GetYaxis()->CenterTitle();
      //h1->GetXaxis()->SetTitleSize(0.05);
      //h1->GetXaxis()->SetTitleOffset(0.80);
      //h1->GetYaxis()->SetTitleOffset(1.4);
    }
    if(obj->InheritsFrom("TH1D")) {
      h1d = (TH1D*) obj;
      h1d->GetXaxis()->CenterTitle();
      h1d->GetYaxis()->CenterTitle();
      //h1d->GetXaxis()->SetTitleSize(0.05);
      //h1d->GetXaxis()->SetTitleOffset(0.80);
      //h1d->GetYaxis()->SetTitleOffset(1.4);
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

  TCanvas *c = NULL;
  TSeqCollection *SCol = gROOT->GetListOfCanvases();
  int size = SCol->GetSize();
  TIter next(SCol);
  TString pdfname = Form("AfterDecompos_v%d_dE%d_sysud%d.pdf",Version,dEcut,sysud);
  for(int i=0;i<size;i++){
    //pdf->NewPage();
    c= (TCanvas*)next();
    c->Modified();
    c->Update();
    c->Draw();
    c->cd();
    //inside the canvas
    //TPaveText *pt = new TPaveText(.74,.81,0.9,0.90,"NDC");
    c->Modified();
    c->Update();
    //std::cout << c->GetName() << std::endl;
    //make 1 pdf file
    if(i==0) c->Print(pdfname+"(",Form("pdf Title:%s",c->GetTitle()));
    else if(i==size-1)c->Print(pdfname+")",Form("pdf Title:%s",c->GetTitle())); 
    else c->Print(pdfname,Form("pdf Title:%s",c->GetTitle())); 
    //make separated pdf files
    //c->Print(Form("pdf/%s.pdf",c->GetTitle()));
  }
  
  std::cout << __LINE__ << std::endl;

  TFile* fout = new TFile(Form("cs_pisigma_v%d_dE%d_sys%d.root",Version,dEcut,sysud),"RECREATE");
  fout->Print();
  fout->cd();
  for(int iq=0;iq<4;iq++){
    for(int isys=0;isys<3;isys++){
      q_IMnpipi_Sp_cs[iq][isys]->Write();
      q_IMnpipi_Sm_cs[iq][isys]->Write();
      q_IMnpipi_K0_cs[iq][isys]->Write();
      q_IMnpipi_SpSmSum[iq][isys]->Write();
      IMnpipi_Sp_cs[iq][isys]->Write();
      IMnpipi_Sm_cs[iq][isys]->Write();
      IMnpipi_K0_cs[iq][isys]->Write();
      IMnpipi_SpSmSum[iq][isys]->Write();
    }
    gDecoErrorSp_CS[iq]->Write();
    gDecoErrorSm_CS[iq]->Write();
    gDecoErrorK0_CS[iq]->Write();
  }
  fout->Close();

  return;

}
