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

const int Version = 239;

void plot_AfterDecompos(const int dEcut=2,const int sysud=0)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetMarkerStyle(20); 
  gStyle->SetMarkerSize(1.2); 
  gROOT->ForceStyle();

  TFile *fr[4] = {NULL};
  //Because the statistics is limited, we just divide data into q<0.35 and q>0.35.
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
  gStyle->SetErrorX(0.);  

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
    IMnpipi_wSid_n_Sp[iq]->Draw("HE");
    IMnpipi_wK0_wSid_n_Sp[iq] = (TH1D*)q_IMnpipi_wK0_wSid_n_Sp[iq]->ProjectionX(Form("IMnpipi_wK0_wSid_n_Sp_%d",iq));
    IMnpipi_wK0_wSid_n_Sp[iq]->SetLineColor(2);
    IMnpipi_wK0_wSid_n_Sp[iq]->Draw("HEsame");
    IMnpipi_wSid_n_SpSm[iq] = (TH1D*)q_IMnpipi_wSid_n_SpSm[iq]->ProjectionX(Form("IMnpipi_wSid_n_SpSm_%d",iq));
    IMnpipi_wSid_n_SpSm[iq]->SetLineColor(3);
    IMnpipi_wSid_n_SpSm[iq]->Draw("HEsame");
    
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
    IMnpipi_wSid_n_Sm[iq]->Draw("HE");
    IMnpipi_wK0_wSid_n_Sm[iq] = (TH1D*)q_IMnpipi_wK0_wSid_n_Sm[iq]->ProjectionX(Form("IMnpipi_wK0_wSid_n_Sm_%d",iq));
    IMnpipi_wK0_wSid_n_Sm[iq]->SetLineColor(2);
    IMnpipi_wK0_wSid_n_Sm[iq]->Draw("HEsame");
    IMnpipi_wSid_n_SpSm[iq]->Draw("HEsame");
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
    IMnpipi_wK0_n[iq]->Draw("HE");
    IMnpipi_wK0_wSid_n_Sp2[iq] = (TH1D*)IMnpipi_wK0_wSid_n_Sp[iq]->Clone(Form("IMnpipi_wK0_wSid_n_Sp2_%d",iq));
    IMnpipi_wK0_wSid_n_Sm2[iq] = (TH1D*)IMnpipi_wK0_wSid_n_Sm[iq]->Clone(Form("IMnpipi_wK0_wSid_n_Sm2_%d",iq));
    IMnpipi_wK0_wSid_n_Sp2[iq]->SetLineColor(2);
    IMnpipi_wK0_wSid_n_Sm2[iq]->SetLineColor(3);
    IMnpipi_wK0_wSid_n_Sp2[iq]->Draw("HEsame");
    IMnpipi_wK0_wSid_n_Sm2[iq]->Draw("HEsame");
    
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
    IMnpipi_wK0_wSid_n_SpSm[iq]->Draw("HISTE");
     

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
    IMnpipi_wK0_wSid_n_Sp[iq]->Draw("HE");
    IMnpipi_K0orSp_ToSp[iq][1]->SetLineColor(3);
    IMnpipi_K0orSp_ToSp[iq][1]->Draw("HEsame");
    IMnpipi_K0orSp_ToK0[iq][1]->SetLineColor(4);
    IMnpipi_K0orSp_ToK0[iq][1]->Draw("HEsame");

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
    IMnpipi_wK0_wSid_n_Sm[iq]->Draw("HE");
    IMnpipi_K0orSm_ToSm[iq][1]->SetLineColor(3);
    IMnpipi_K0orSm_ToSm[iq][1]->Draw("HEsame");
    IMnpipi_K0orSm_ToK0[iq][1]->SetLineColor(4);
    IMnpipi_K0orSm_ToK0[iq][1]->Draw("HEsame");
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
    nSp_SporSm[iqlowhigh][0] = nSp-nSp_bin_err;
    nSm_SporSm[iqlowhigh][0] = nSm-nSm_bin_err;
    nSp_SporSm[iqlowhigh][1] = nSp;
    nSm_SporSm[iqlowhigh][1] = nSm;
    nSp_SporSm[iqlowhigh][2] = nSp+nSp_bin_err;
    nSm_SporSm[iqlowhigh][2] = nSm+nSm_bin_err;
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
  
  std::cout << __LINE__ << std::endl;
  TCanvas *cIMnpipi_SpSm_afterDeco[4];
  for(int iq=0;iq<4;iq++){
    cIMnpipi_SpSm_afterDeco[iq] = new TCanvas(Form("cIMnpipi_SpSm_Deco_%s",cqcut[iq+1]),Form("cIMnpipi_SpSm_Deco_%s",cqcut[iq]),1000,800);
    IMnpipi_wSid_n_SpSm[iq]->SetTitle(Form("IM(n#pi^{+}#pi^{-}) #Sigma^{+}#pi^{-} or #Sigma^{-}#pi^{+} %s",cqcut[iq]));
    IMnpipi_wSid_n_SpSm[iq]->SetYTitle("Counts/ 15 MeV");
    IMnpipi_wSid_n_SpSm[iq]->GetYaxis()->CenterTitle();
    IMnpipi_wSid_n_SpSm[iq]->Draw("HE");
    IMnpipi_SporSm_ToSp[iq][1]->SetLineColor(4);
    IMnpipi_SporSm_ToSp[iq][1]->Draw("HEsame");
    IMnpipi_SporSm_ToSm[iq][1]->SetLineColor(5);
    IMnpipi_SporSm_ToSm[iq][1]->Draw("HEsame");
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
    IMnpipi_wK0_wSid_n_SpSm[iq]->Draw("HE");
    IMnpipi_K0SpSm_ToK0[iq][1]->SetLineColor(4);
    IMnpipi_K0SpSm_ToK0[iq][1]->Draw("HEsame");
    IMnpipi_K0SpSm_ToSp[iq][1]->SetLineColor(2);
    IMnpipi_K0SpSm_ToSp[iq][1]->Draw("HEsame");
    IMnpipi_K0SpSm_ToSm[iq][1]->SetLineColor(3);
    IMnpipi_K0SpSm_ToSm[iq][1]->Draw("HEsame");
    TLegend *lo = new TLegend(0.6,0.7,0.9,0.9);
    lo->AddEntry(IMnpipi_K0SpSm_ToK0[iq][1],"To K0 mode","l");
    lo->AddEntry(IMnpipi_K0SpSm_ToSp[iq][1],"To Sp mode","l");
    lo->AddEntry(IMnpipi_K0SpSm_ToSm[iq][1],"To Sm mode","l");
    lo->Draw();
  }


  TCanvas *cq_IMnpipi_Sp_afterDeco[4];
  TH2D* q_IMnpipi_Sp_sum[4];
  for(int iq=0;iq<4;iq++){
    q_IMnpipi_Sp_sum[iq] = (TH2D*)q_IMnpipi_wSid_n_Sp[iq]->Clone(Form("q_IMnpipi_Sp_sum%d",iq));
    q_IMnpipi_Sp_sum[iq]->Add(q_IMnpipi_K0orSp_ToK0[iq][1],-1.0);
    q_IMnpipi_Sp_sum[iq]->Add(q_IMnpipi_SporSm_ToSm[iq][1],-1.0);
    q_IMnpipi_Sp_sum[iq]->Add(q_IMnpipi_K0SpSm_ToSp[iq][1],1.0);//added over-subtraction component
    cq_IMnpipi_Sp_afterDeco[iq] = new TCanvas(Form("cq_IMnpipi_Sp_afterDeco%d",iq),Form("cq_IMnpipi_Sp_afterDeco%d",iq));
    cq_IMnpipi_Sp_afterDeco[iq]->cd();
    q_IMnpipi_Sp_sum[iq]->SetMinimum(0);
    q_IMnpipi_Sp_sum[iq]->Draw("colz");
  }

  TCanvas *cq_IMnpipi_Sm_afterDeco[4];
  TH2D* q_IMnpipi_Sm_sum[4];
  for(int iq=0;iq<4;iq++){
    q_IMnpipi_Sm_sum[iq] = (TH2D*)q_IMnpipi_wSid_n_Sm[iq]->Clone(Form("q_IMnpipi_Sm_sum%d",iq));
    q_IMnpipi_Sm_sum[iq]->Add(q_IMnpipi_K0orSm_ToK0[iq][1],-1.0);
    q_IMnpipi_Sm_sum[iq]->Add(q_IMnpipi_SporSm_ToSp[iq][1],-1.0);
    q_IMnpipi_Sm_sum[iq]->Add(q_IMnpipi_K0SpSm_ToSm[iq][1],1.0);//added over-subtraction component
    cq_IMnpipi_Sm_afterDeco[iq] = new TCanvas(Form("cq_IMnpipi_Sm_afterDeco%d",iq),Form("cq_IMnpipi_Sm_afterDeco%d",iq));
    cq_IMnpipi_Sm_afterDeco[iq]->cd();
    q_IMnpipi_Sm_sum[iq]->SetMinimum(0);
    q_IMnpipi_Sm_sum[iq]->Draw("colz");
  }

  std::cout << __LINE__ << std::endl;

  TCanvas *cq_IMnpipi_K0_afterDeco[4];
  TH2D* q_IMnpipi_K0_sum[4];
  for(int iq=0;iq<4;iq++){
    q_IMnpipi_K0_sum[iq] = (TH2D*)q_IMnpipi_wK0_n[iq]->Clone(Form("q_IMnpipi_K0_sum%d",iq));
    q_IMnpipi_K0_sum[iq]->Add(q_IMnpipi_K0orSp_ToSp[iq][1],-1.0);
    q_IMnpipi_K0_sum[iq]->Add(q_IMnpipi_K0orSm_ToSm[iq][1],-1.0);
    q_IMnpipi_K0_sum[iq]->Add(q_IMnpipi_K0SpSm_ToK0[iq][1],1.0);//added over-subtraction component
    cq_IMnpipi_K0_afterDeco[iq] = new TCanvas(Form("cq_IMnpipi_K0_afterDeco%d",iq),Form("cq_IMnpipi_K0_afterDeco%d",iq));
    cq_IMnpipi_K0_afterDeco[iq]->cd();
    q_IMnpipi_K0_sum[iq]->SetMinimum(0);
    q_IMnpipi_K0_sum[iq]->Draw("colz");
  }
  std::cout << __LINE__ << std::endl;



  TCanvas *cIMnpipi_Summary_Sp[4];
  TH1D* IMnpipi_Sp_noK0_noSm[4];
  for(int iq=0;iq<4;iq++){
    cIMnpipi_Summary_Sp[iq] = new TCanvas(Form("cIMnpipi_Summary_Sp_%s",cqcut[iq]),Form("cIMnpipi_Summary_Sp_%s",cqcut[iq]),1000,800);
    //all sum
    IMnpipi_wSid_n_Sp[iq]->SetYTitle("Counts/ 15 MeV");
    IMnpipi_wSid_n_Sp[iq]->GetYaxis()->CenterTitle();
    IMnpipi_wSid_n_Sp[iq]->Draw("E");
    IMnpipi_Sp_noK0_noSm[iq] = (TH1D*)IMnpipi_wSid_n_Sp[iq]->Clone(Form("IMnpipi_Sp_noK0_noSm%d",iq));
    IMnpipi_Sp_noK0_noSm[iq]->Add(IMnpipi_K0orSp_ToK0[iq][1],-1);
    IMnpipi_Sp_noK0_noSm[iq]->Add(IMnpipi_SporSm_ToSm[iq][1],-1);
    IMnpipi_Sp_noK0_noSm[iq]->Add(IMnpipi_K0SpSm_ToSp[iq][1],1);//add over-subtraction component
    IMnpipi_Sp_noK0_noSm[iq]->SetLineColor(2);
    IMnpipi_Sp_noK0_noSm[iq]->Draw("Esame");
    TLegend *lSp = new TLegend(0.6,0.7,0.9,0.9);
    lSp->AddEntry(IMnpipi_wSid_n_Sp[iq],"Sigma+ like mode","l");
    lSp->AddEntry(IMnpipi_Sp_noK0_noSm[iq],"Sigma+ after deco. ","l");
    lSp->Draw();
  }
  
  //summary plot with the systematic error of 
  TCanvas *cIMnpipi_Summary_Sp_re[4];
  for(int iq=0;iq<4;iq++){
    cIMnpipi_Summary_Sp_re[iq] = new TCanvas(Form("cIMnpipi_Summary_Sp_re%s",cqcut[iq]),Form("cIMnpipi_Summary_Sp_re%s",cqcut[iq]),1000,800);
    //all sum
    IMnpipi_Sp_noK0_noSm[iq]->SetYTitle("Counts/ 15 MeV");
    IMnpipi_Sp_noK0_noSm[iq]->GetYaxis()->CenterTitle();
    IMnpipi_Sp_noK0_noSm[iq]->Draw("E");
  }

  TCanvas *cIMnpipi_Summary_Sm[4];
  TH1D* IMnpipi_Sm_noK0_noSp[4];
  for(int iq=0;iq<4;iq++){
    cIMnpipi_Summary_Sm[iq] = new TCanvas(Form("cIMnpipi_Summary_Sm_%s",cqcut[iq]),Form("cIMnpipi_Summary_Sm_%s",cqcut[iq]),1000,800);
    //all sum
    IMnpipi_wSid_n_Sm[iq]->SetYTitle("Counts/ 15 MeV");
    IMnpipi_wSid_n_Sm[iq]->GetYaxis()->CenterTitle();
    IMnpipi_wSid_n_Sm[iq]->Draw("E");
    IMnpipi_Sm_noK0_noSp[iq] = (TH1D*)IMnpipi_wSid_n_Sm[iq]->Clone(Form("IMnpipi_Sm_noK0_noSp%d",iq));
    IMnpipi_Sm_noK0_noSp[iq]->Add(IMnpipi_K0orSm_ToK0[iq][1],-1);
    IMnpipi_Sm_noK0_noSp[iq]->Add(IMnpipi_SporSm_ToSp[iq][1],-1);
    IMnpipi_Sm_noK0_noSp[iq]->Add(IMnpipi_K0SpSm_ToSm[iq][1],1);//add over-subtraction component
    IMnpipi_Sm_noK0_noSp[iq]->SetLineColor(2);
    IMnpipi_Sm_noK0_noSp[iq]->Draw("Esame");
    TLegend *lSm = new TLegend(0.6,0.7,0.9,0.9);
    lSm->AddEntry(IMnpipi_wSid_n_Sm[iq],"Sigma- like mode","l");
    lSm->AddEntry(IMnpipi_Sm_noK0_noSp[iq],"Sigma- after deco. ","l");
    lSm->Draw();
  }

  TCanvas *cIMnpipi_Summary_Sm_re[4];
  for(int iq=0;iq<4;iq++){
    cIMnpipi_Summary_Sm_re[iq] = new TCanvas(Form("cIMnpipi_Summary_Sm_re%s",cqcut[iq]),Form("cIMnpipi_Summary_Sm_re%s",cqcut[iq]),1000,800);
    //all sum
    IMnpipi_Sm_noK0_noSp[iq]->SetYTitle("Counts/ 15 MeV");
    IMnpipi_Sm_noK0_noSp[iq]->GetYaxis()->CenterTitle();;
    IMnpipi_Sm_noK0_noSp[iq]->Draw("E");
  }


  TCanvas *cIMnpipi_Summary_K0[4];
  TH1D* IMnpipi_K0_noSp_noSm[4];
  for(int iq=0;iq<4;iq++){
    cIMnpipi_Summary_K0[iq] = new TCanvas(Form("cIMnpipi_Summary_K0_%s",cqcut[iq]),Form("cIMnpipi_Summary_K0_%s",cqcut[iq]),1000,800);
    //all sum
    IMnpipi_wK0_n[iq]->SetYTitle("Counts/ 15 MeV");
    IMnpipi_wK0_n[iq]->GetYaxis()->CenterTitle();
    IMnpipi_wK0_n[iq]->Draw("E");
    IMnpipi_K0_noSp_noSm[iq] = (TH1D*)IMnpipi_wK0_n[iq]->Clone(Form("IMnpipi_K0_noSp_noSm%d",iq));
    IMnpipi_K0_noSp_noSm[iq]->Add(IMnpipi_K0orSp_ToSp[iq][1],-1);
    IMnpipi_K0_noSp_noSm[iq]->Add(IMnpipi_K0orSm_ToSm[iq][1],-1);
    IMnpipi_K0_noSp_noSm[iq]->Add(IMnpipi_K0SpSm_ToK0[iq][1],1);//add over-subtraction component
    IMnpipi_K0_noSp_noSm[iq]->SetLineColor(2);
    IMnpipi_K0_noSp_noSm[iq]->Draw("Esame");
    TLegend *lK0 = new TLegend(0.6,0.7,0.9,0.9);
    lK0->AddEntry(IMnpipi_wK0_n[iq],"K0nn like mode","l");
    lK0->AddEntry(IMnpipi_K0_noSp_noSm[iq],"K0nn after deco. ","l");
    lK0->Draw();
  }

  TCanvas *cIMnpipi_Summary_K0_re[4];
  for(int iq=0;iq<4;iq++){
    cIMnpipi_Summary_K0_re[iq] = new TCanvas(Form("cIMnpipi_Summary_K0_re%s",cqcut[iq]),Form("cIMnpipi_Summary_K0_re%s",cqcut[iq]),1000,800);
    //all sum
    IMnpipi_K0_noSp_noSm[iq]->SetYTitle("Counts/ 15 MeV");
    IMnpipi_K0_noSp_noSm[iq]->GetYaxis()->CenterTitle();;
    IMnpipi_K0_noSp_noSm[iq]->Draw("E");
  }


  //
  //Acceptance correction
  //

  std::cout << __LINE__ << std::endl;
  //TFile *facc = TFile::Open("../simpost/accmap.root");
  TFile *facc = TFile::Open("../simpost/accmapv152_27.root");
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


  TH2D* q_IMnpipi_Sp_cs[3];
  TH2D* q_IMnpipi_Sm_cs[3];
  TH2D* q_IMnpipi_K0_cs[3];
  TH2D* q_IMnpipi_Sp_cserr[3];
  TH2D* q_IMnpipi_Sm_cserr[3];
  TH2D* q_IMnpipi_K0_cserr[3];
  const int qcut650 =q_IMnpipi_Sp_sum[0]->GetYaxis()->FindBin(0.65);
  for(int iq=0;iq<3;iq++){
    q_IMnpipi_Sp_cs[iq] = (TH2D*)q_IMnpipi_Sp_sum[iq]->Clone(Form("q_IMnpipi_Sp_cs%d",iq));
    q_IMnpipi_Sm_cs[iq] = (TH2D*)q_IMnpipi_Sm_sum[iq]->Clone(Form("q_IMnpipi_Sm_cs%d",iq));
    q_IMnpipi_K0_cs[iq] = (TH2D*)q_IMnpipi_K0_sum[iq]->Clone(Form("q_IMnpipi_K0_cs%d",iq));
    q_IMnpipi_Sp_cs[iq]->SetTitle(Form("q_IMnpipi_Sp_cs_%d",iq));
    q_IMnpipi_Sm_cs[iq]->SetTitle(Form("q_IMnpipi_Sm_cs_%d",iq));
    q_IMnpipi_K0_cs[iq]->SetTitle(Form("q_IMnpipi_K0_cs_%d",iq));
    
    for(int ix=0;ix<q_IMnpipi_Sp_cs[iq]->GetNbinsX();ix++){
      for(int iy=0;iy<q_IMnpipi_Sp_cs[iq]->GetNbinsY();iy++){
        double contSp =q_IMnpipi_Sp_cs[iq]->GetBinContent(ix,iy);
        double contSperr =q_IMnpipi_Sp_cs[iq]->GetBinError(ix,iy);
        double accSp =q_IMnpipi_Sp_accp->GetBinContent(ix,iy);
        double accerrSp =q_IMnpipi_Sp_accperr->GetBinContent(ix,iy);
        double csSp = 0.0; 
        double csSperr = 0.0; 
        double binwidth = q_IMnpipi_Sp_cs[iq]->ProjectionX()->GetBinWidth(1)*1000.0;
        if(accSp>0.0){
          csSp = contSp/accSp/binwidth/trigScale/lumi ;
          csSperr = contSperr/accSp/binwidth/trigScale/lumi;
        }
        double contSm =q_IMnpipi_Sm_cs[iq]->GetBinContent(ix,iy);
        double contSmerr =q_IMnpipi_Sm_cs[iq]->GetBinError(ix,iy);
        double accSm =q_IMnpipi_Sm_accp->GetBinContent(ix,iy);
        double accerrSm =q_IMnpipi_Sm_accperr->GetBinContent(ix,iy);
        double csSm = 0.0; 
        double csSmerr = 0.0; 
        if(accSm>0.0){
          csSm = contSm/accSm/binwidth/trigScale/lumi;
          csSmerr = contSmerr/accSm/binwidth/trigScale/lumi;
        }

        double contK0 =q_IMnpipi_K0_cs[iq]->GetBinContent(ix,iy);
        double contK0err =q_IMnpipi_K0_cs[iq]->GetBinError(ix,iy);
        double accK0 =q_IMnpipi_K0_accp->GetBinContent(ix,iy);
        double accerrK0 =q_IMnpipi_K0_accperr->GetBinContent(ix,iy);
        double csK0 = 0.0; 
        double csK0err = 0.0; 
        if(accK0>0.0){
          csK0 = contK0/accK0/binwidth/trigScale/lumi;
          csK0err = contK0err/accK0/binwidth/trigScale/lumi;
        }

        if(accerrSp<UncertCut){
          q_IMnpipi_Sp_cs[iq]->SetBinContent(ix,iy,csSp);
          q_IMnpipi_Sp_cs[iq]->SetBinError(ix,iy,csSperr);
        }else{
          q_IMnpipi_Sp_cs[iq]->SetBinContent(ix,iy,0.);
          q_IMnpipi_Sp_cs[iq]->SetBinError(ix,iy,0.);
        }
        
        //if(iy>qcut650){
        //  q_IMnpipi_Sp_cs[iq]->SetBinContent(ix,iy,0.);
        //  q_IMnpipi_Sp_cs[iq]->SetBinError(ix,iy,0.);
        //}
        if(accerrSm<UncertCut){
          q_IMnpipi_Sm_cs[iq]->SetBinContent(ix,iy,csSm);
          q_IMnpipi_Sm_cs[iq]->SetBinError(ix,iy,csSmerr);
        }else{
          q_IMnpipi_Sm_cs[iq]->SetBinContent(ix,iy,0.);
          q_IMnpipi_Sm_cs[iq]->SetBinError(ix,iy,0.);
        }
        
        //if(iy>qcut650){
        //  q_IMnpipi_Sm_cs[iq]->SetBinContent(ix,iy,0.);
        //  q_IMnpipi_Sm_cs[iq]->SetBinError(ix,iy,0.);
        //}
        
        if(accerrK0<UncertCut){
          q_IMnpipi_K0_cs[iq]->SetBinContent(ix,iy,csK0);
          q_IMnpipi_K0_cs[iq]->SetBinError(ix,iy,csK0err);
        }else{
          q_IMnpipi_K0_cs[iq]->SetBinContent(ix,iy,0.);
          q_IMnpipi_K0_cs[iq]->SetBinError(ix,iy,0.);
        }
        
        //if(iy>qcut650){
        //  q_IMnpipi_K0_cs[iq]->SetBinContent(ix,iy,0.);
        //  q_IMnpipi_K0_cs[iq]->SetBinError(ix,iy,0.);
        //}
      }
    }
  }

  TCanvas *ccsSp[3];
  TCanvas *ccsSm[3];
  TCanvas *ccsK0[3];
  TH1D* IMnpipi_Sp_cs_single[3];
  TH1D* IMnpipi_Sm_cs_single[3];
  TH1D* IMnpipi_K0_cs_single[3];
  for(int iq=0;iq<3;iq++){
    ccsSp[iq] = new TCanvas(Form("ccsSp%d",iq),Form("ccsSp%d",iq),1600,800);
    ccsSp[iq]->Divide(2,1);
    ccsSp[iq]->cd(1);
    q_IMnpipi_Sp_cs[iq]->Draw("colz");
    ccsSp[iq]->cd(2);
    IMnpipi_Sp_cs_single[iq]= (TH1D*)q_IMnpipi_Sp_cs[iq]->ProjectionX(Form("IMnpipi_Sp_cs_single%d",iq));
    IMnpipi_Sp_cs_single[iq]->SetYTitle("d#sigma/dM [#mu b (MeV/c^{2})]");
    IMnpipi_Sp_cs_single[iq]->GetYaxis()->CenterTitle();
    IMnpipi_Sp_cs_single[iq]->Draw("HE");
    ccsSm[iq] = new TCanvas(Form("ccsSm%d",iq),Form("ccsSm%d",iq),1600,800);
    ccsSm[iq]->Divide(2,1);
    ccsSm[iq]->cd(1);
    q_IMnpipi_Sm_cs[iq]->Draw("colz");
    ccsSm[iq]->cd(2);
    IMnpipi_Sm_cs_single[iq]= (TH1D*)q_IMnpipi_Sm_cs[iq]->ProjectionX(Form("IMnpipi_Sm_cs_single%d",iq));
    IMnpipi_Sm_cs_single[iq]->SetYTitle("d#sigma/dM [#mu b (MeV/c^{2})]");
    IMnpipi_Sm_cs_single[iq]->GetYaxis()->CenterTitle();
    IMnpipi_Sm_cs_single[iq]->Draw("HE");

    ccsK0[iq] = new TCanvas(Form("ccsK0%d",iq),Form("ccsK0%d",iq),1600,800);
    ccsK0[iq]->Divide(2,1);
    ccsK0[iq]->cd(1);
    q_IMnpipi_K0_cs[iq]->Draw("colz");
    ccsK0[iq]->cd(2);
    IMnpipi_K0_cs_single[iq]= (TH1D*)q_IMnpipi_K0_cs[iq]->ProjectionX(Form("IMnpipi_K0_cs_single%d",iq));
    IMnpipi_K0_cs_single[iq]->SetYTitle("d#sigma/dM [#mu b (MeV/c^{2})]");
    IMnpipi_K0_cs_single[iq]->GetYaxis()->CenterTitle();
    IMnpipi_K0_cs_single[iq]->Draw("HE");
  }
  
  
  TCanvas *csum[3];
  TH2D* q_IMnpipi_SpSmSum[3];
  for(int iq=0;iq<3;iq++){
    csum[iq]  = new TCanvas(Form("csum%d",iq),Form("csum%d",iq),1600,800);
    csum[iq]->Divide(2,1);
    q_IMnpipi_SpSmSum[iq] = (TH2D*)q_IMnpipi_Sp_cs[iq]->Clone(Form("q_IMnpipi_SpSmSum%d",iq));
    q_IMnpipi_SpSmSum[iq]->Add(q_IMnpipi_Sm_cs[iq],1.0);
    csum[iq]->cd(1);
    q_IMnpipi_SpSmSum[iq]->SetTitle(Form("q_IMnpipi_SpSmSum_%d",iq));
    q_IMnpipi_SpSmSum[iq]->SetMinimum(0);
    q_IMnpipi_SpSmSum[iq]->Draw("colz");
    csum[iq]->cd(2);
    q_IMnpipi_SpSmSum[iq]->SetTitle("#Sigma^{+} #Sigma^{-} charge Sum");
    q_IMnpipi_SpSmSum[iq]->ProjectionX()->Draw("HE");
  }

  TCanvas *csum2D = new TCanvas("csum2D","csum2D");
  csum2D->cd();
  q_IMnpipi_SpSmSum[0]->Draw("col");
  q_IMnpipi_SpSmSum[1]->Draw("colsame");





  TCanvas *csub = new TCanvas("csub","csub",1600,800);
  csub->Divide(2,1);
  TH2D* q_IMnpipi_SpSmSub[3];
  for(int iq=0;iq<3;iq++){
    q_IMnpipi_SpSmSub[iq] = (TH2D*)q_IMnpipi_Sp_cs[iq]->Clone(Form("q_IMnpipi_SpSmSub%d",iq));
    q_IMnpipi_SpSmSub[iq]->Add(q_IMnpipi_Sm_cs[iq],-1.0);
  }
  csub->cd(1);
  q_IMnpipi_SpSmSub[1]->SetTitle("q_IMnpipi_SpSmSub");
  q_IMnpipi_SpSmSub[1]->Draw("colz");
  csub->cd(2);
  q_IMnpipi_SpSmSub[1]->ProjectionX()->Draw("HE");
   
  TCanvas *ccomp[3];
  TH1D* IMnpipi_Sp_cs[3];
  TH1D* IMnpipi_Sm_cs[3];
  TH1D* IMnpipi_K0_cs[3];
  
  for(int iq=0;iq<3;iq++){
    ccomp[iq] = new TCanvas(Form("ccomp%d",iq),Form("ccomp%d",iq),1000,800);
    IMnpipi_Sp_cs[iq] = (TH1D*)q_IMnpipi_Sp_cs[iq]->ProjectionX(Form("IMnpipi_Sp_cs%d",iq));
    IMnpipi_Sm_cs[iq] = (TH1D*)q_IMnpipi_Sm_cs[iq]->ProjectionX(Form("IMnpipi_Sm_cs%d",iq));
    IMnpipi_K0_cs[iq] = (TH1D*)q_IMnpipi_K0_cs[iq]->ProjectionX(Form("IMnpipi_K0_cs%d",iq));
    IMnpipi_Sp_cs[iq]->SetLineColor(3);
    IMnpipi_Sm_cs[iq]->SetLineColor(4);
    IMnpipi_K0_cs[iq]->SetLineColor(2);
    //IMnpipi_K0_cs[iq]->RebinX(2);
    //IMnpipi_Sp_cs[iq]->RebinX(2);
    //IMnpipi_Sm_cs[iq]->RebinX(2);
    //IMnpipi_K0_cs[iq]->Draw("HE");
    //IMnpipi_Sp_cs[iq]->Draw("HEsame");
    IMnpipi_Sp_cs[iq]->SetMarkerStyle(20);
    IMnpipi_Sp_cs[iq]->Draw("E");
    IMnpipi_Sp_cs[iq]->SetYTitle("d#sigma/dM [#mu b (MeV/c^{2})]");
    IMnpipi_Sp_cs[iq]->GetYaxis()->CenterTitle();
    IMnpipi_Sm_cs[iq]->Draw("Esame");
    //IMnpipi_K0_cs[iq]->Draw("HEsame");
  
    TLegend *lcs = new TLegend(0.6,0.7,0.9,0.9);
    lcs->AddEntry(IMnpipi_Sp_cs[iq],"#Sigma^{+} mode","l");
    lcs->AddEntry(IMnpipi_Sm_cs[iq],"#Sigma^{-} mode","l");
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
  TString pdfname = "AfterDecompos.pdf";
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
  
  TFile* fout = new TFile("cs_pisigma.root","RECREATE");
  fout->Print();
  fout->cd();
  for(int iq=0;iq<3;iq++){
    q_IMnpipi_Sp_cs[iq]->Write();
    q_IMnpipi_Sm_cs[iq]->Write();
    q_IMnpipi_K0_cs[iq]->Write();
    IMnpipi_Sp_cs[iq]->Write();
    IMnpipi_Sm_cs[iq]->Write();
    IMnpipi_K0_cs[iq]->Write();
  }
  fout->Close();

  return;

}
