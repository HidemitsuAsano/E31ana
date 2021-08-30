const bool showBG = true;
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

#include "../src/GlobalVariables.h"

bool RemoveNegative = true;
bool RebinMode = true;
bool Sidefar=false;
bool FitNoWeight=true;


void plot_AfterDecompos()
{
  TFile *fr[4]={NULL};
  
  fr[0] = TFile::Open("evanaIMpisigma_npippim_v229_out_iso_nostop_sub.root");
  fr[1] = TFile::Open("evanaIMpisigma_npippim_v229_out_iso_qlo_nostop_sub.root");
  fr[2] = TFile::Open("evanaIMpisigma_npippim_v229_out_iso_qhi_nostop_sub.root");
  fr[3] = TFile::Open("evanaIMpisigma_npippim_v229_out_iso_theta15_nostop_sub.root");
  fr[0]->Print();
  fr[1]->Print();
  fr[2]->Print();
  fr[3]->Print();
  TFile *fdeco[2]={NULL};
  fdeco[0] = TFile::Open("deco_qlo.root","READ");
  fdeco[1] = TFile::Open("deco_qhi.root","READ");
  fdeco[0]->Print();
  fdeco[1]->Print();

  TFile *flumi = TFile::Open("InteLumi.root");
  TParameter<double>*IntegLumi = flumi->Get("IntegLumi");
  TParameter<double>*Err = flumi->Get("Err");
  double lumi = IntegLumi->GetVal();
  double lumierr = Err->GetVal();
  const double trigScale = 2.0;
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
  const int qstart=1;

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
    //get
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

    //Draw
    cq_IMnpipi_wSid_n_Sp[iq] = new TCanvas(Form("cq_IMnpipi_wSid_n_Sp_%s",cqcut[iq]),Form("cq_IMnpipi_wSid_n_Sp_%s",cqcut[iq]));
    q_IMnpipi_wSid_n_Sp[iq]->Draw("colz");
    cq_IMnpipi_wK0_wSid_n_Sp[iq] = new TCanvas(Form("cq_IMnpipi_wK0_wSid_n_Sp_%s",cqcut[iq]),Form("cq_IMnpipi_wK0_wSid_n_Sp_%s",cqcut[iq]));
    q_IMnpipi_wK0_wSid_n_Sp[iq]->Draw("colz");
    
    cIMnpipi_Sp[iq] = new TCanvas(Form("IMnpipi_Sp_%s",cqcut[iq]),Form("IMnpipi_Sp_%s",cqcut[iq]),800,800);
    IMnpipi_wSid_n_Sp[iq] = (TH1D*)q_IMnpipi_wSid_n_Sp[iq]->ProjectionX(Form("IMnpipi_wSid_n_Sp_%d",iq));
    IMnpipi_wSid_n_Sp[iq]->Draw("HE");
    IMnpipi_wK0_wSid_n_Sp[iq] = (TH1D*)q_IMnpipi_wK0_wSid_n_Sp[iq]->ProjectionX(Form("IMnpipi_wK0_wSid_n_Sp_%d",iq));
    IMnpipi_wK0_wSid_n_Sp[iq]->SetLineColor(2);
    IMnpipi_wK0_wSid_n_Sp[iq]->Draw("HEsame");
    IMnpipi_wSid_n_SpSm[iq] = (TH1D*)q_IMnpipi_wSid_n_SpSm[iq]->ProjectionX(Form("IMnpipi_wSid_n_SpSm_%d",iq));
    IMnpipi_wSid_n_SpSm[iq]->SetLineColor(3);
    IMnpipi_wSid_n_SpSm[iq]->Draw("HEsame");
    
    TLegend *lSpoverlap = new TLegend(0.6,0.7,0.9,0.9);
    lSpoverlap->AddEntry(IMnpipi_wSid_n_Sp[iq],"#Sigma+ like mode","l");
    lSpoverlap->AddEntry(IMnpipi_wK0_wSid_n_Sp[iq],"K0 overlap","l");
    lSpoverlap->AddEntry(IMnpipi_wSid_n_SpSm[iq],"#Sigma- overlap","l");
    lSpoverlap->Draw();

    cq_IMnpipi_wSid_n_Sm[iq] = new TCanvas(Form("cq_IMnpipi_wSid_n_Sm_%s",cqcut[iq]),Form("cq_IMnpipi_wSid_n_Sm_%s",cqcut[iq]));
    q_IMnpipi_wSid_n_Sm[iq]->Draw("colz");
    cq_IMnpipi_wK0_wSid_n_Sm[iq] = new TCanvas(Form("cq_IMnpipi_wK0_wSid_n_Sm_%s",cqcut[iq]),Form("cq_IMnpipi_wK0_wSid_n_Sm_%s",cqcut[iq]));
    q_IMnpipi_wK0_wSid_n_Sm[iq]->Draw("colz");
    cIMnpipi_Sm[iq] = new TCanvas(Form("IMnpipi_Sm_%s",cqcut[iq]),Form("IMnpipi_Sm_%s",cqcut[iq]),800,800);
    IMnpipi_wSid_n_Sm[iq] = (TH1D*)q_IMnpipi_wSid_n_Sm[iq]->ProjectionX(Form("IMnpipi_wSid_n_Sm_%d",iq));
    IMnpipi_wSid_n_Sm[iq]->Draw("HE");
    IMnpipi_wK0_wSid_n_Sm[iq] = (TH1D*)q_IMnpipi_wK0_wSid_n_Sm[iq]->ProjectionX(Form("IMnpipi_wK0_wSid_n_Sm_%d",iq));
    IMnpipi_wK0_wSid_n_Sm[iq]->SetLineColor(2);
    IMnpipi_wK0_wSid_n_Sm[iq]->Draw("HEsame");
    IMnpipi_wSid_n_SpSm[iq]->Draw("HEsame");
    TLegend *lSmoverlap = new TLegend(0.6,0.7,0.9,0.9);
    lSmoverlap->AddEntry(IMnpipi_wSid_n_Sm[iq],"#Sigma- like mode","l");
    lSmoverlap->AddEntry(IMnpipi_wK0_wSid_n_Sm[iq],"K0 overlap","l");
    lSmoverlap->AddEntry(IMnpipi_wSid_n_SpSm[iq],"#Sigma+ overlap","l");
    lSmoverlap->Draw();
    
    cq_IMnpipi_wSid_n_SpSm[iq] = new TCanvas(Form("cq_IMnpipi_wSid_n_SpSm_%s",cqcut[iq]),Form("cq_IMnpipi_wSid_n_SpSm_%s",cqcut[iq]));
    q_IMnpipi_wSid_n_SpSm[iq]->Draw("colz");


    cIMnpipi_K0[iq] = new TCanvas(Form("IMnpipi_K0_%s",cqcut[iq]),Form("IMnpipi_K0_%s",cqcut[iq]),800,800);
    IMnpipi_wK0_n[iq] = (TH1D*)q_IMnpipi_wK0_n[iq]->ProjectionX(Form("IMnpipi_wK0_n_%d",iq));
    IMnpipi_wK0_n[iq]->Draw("HE");
    IMnpipi_wK0_wSid_n_Sp2[iq] = (TH1D*)IMnpipi_wK0_wSid_n_Sp[iq]->Clone(Form("IMnpipi_wK0_wSid_n_Sp2_%d",iq));
    IMnpipi_wK0_wSid_n_Sm2[iq] = (TH1D*)IMnpipi_wK0_wSid_n_Sm[iq]->Clone(Form("IMnpipi_wK0_wSid_n_Sm2_%d",iq));
    IMnpipi_wK0_wSid_n_Sp2[iq]->SetLineColor(2);
    IMnpipi_wK0_wSid_n_Sm2[iq]->SetLineColor(3);
    IMnpipi_wK0_wSid_n_Sp2[iq]->Draw("HEsame");
    IMnpipi_wK0_wSid_n_Sm2[iq]->Draw("HEsame");
    
    TLegend *lK0overlap = new TLegend(0.6,0.7,0.9,0.9);
    lK0overlap->AddEntry(IMnpipi_wK0_n[iq],"K0 like mode","l");
    lK0overlap->AddEntry(IMnpipi_wK0_wSid_n_Sp2[iq],"#Sigma+ overlap","l");
    lK0overlap->AddEntry(IMnpipi_wK0_wSid_n_Sm2[iq],"#Sigma- overlap","l");
    lK0overlap->Draw();


    cq_IMnpipi_wK0_wSid_n_SpSm[iq] = new TCanvas(Form("q_IMnpipi_wK0_wSid_n_SpSm_%s",cqcut[iq]),Form("q_IMnpipi_wK0_wSid_n_SpSm_%s",cqcut[iq]));
    IMnpipi_wK0_wSid_n_SpSm[iq] = (TH1D*)q_IMnpipi_wK0_wSid_n_SpSm[iq]->ProjectionX(Form("IMnpipi_wK0_wSid_n_SpSm_%d",iq));
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
  }

  //only qlo(=0) and qhi(1) decomposition results;
  TH2D* IMnpim_IMnpip_K0inter[2];
  TGraphErrors *gr_SpONnpim_fin_pol1[2];
  TGraphErrors *gr_SmONnpip_fin_pol1[2];
  TGraphErrors *gr_SpONnpim_fin_3rd[2];
  TGraphErrors *gr_SmONnpip_fin_3rd[2];

  //only qlo(=0) and qhi(1) decomposition results;
  for(int iq=0;iq<2;iq++){
    IMnpim_IMnpip_K0inter[iq] = (TH2D*)fdeco[iq]->Get("h2K0inter_3fine");
    gr_SpONnpim_fin_pol1[iq] = (TGraphErrors*)fdeco[iq]->Get("gr_SpONnpim_fin_pol1");
    gr_SmONnpip_fin_pol1[iq] = (TGraphErrors*)fdeco[iq]->Get("gr_SmONnpip_fin_pol1");
    gr_SpONnpim_fin_3rd[iq] = (TGraphErrors*)fdeco[iq]->Get("gr_SpONnpim_fin_3rd");
    gr_SmONnpip_fin_3rd[iq] = (TGraphErrors*)fdeco[iq]->Get("gr_SmONnpip_fin_3rd");
  }
  
  //have data only in Sigma+/- region
  TCanvas *cK0inter = new TCanvas("cK0inter","cK0inter",1600,800);
  cK0inter->Divide(2,1);
  cK0inter->cd(1);
  IMnpim_IMnpip_K0inter[0]->Draw("colz");
  cK0inter->cd(2);
  IMnpim_IMnpip_K0inter[1]->Draw("colz");
  

  TH2D* q_IMnpipi_K0orSp_ToSp[2];
  TH2D* q_IMnpipi_K0orSp_ToK0[2];
  TH1D* IMnpipi_K0orSp_ToSp[2];
  TH1D* IMnpipi_K0orSp_ToK0[2];
   
  //check overlap btw wide-binning IMnpip or IMnpim
  TCanvas *cIMnpim_Spmode_wbinoverlap[2];
  TCanvas *cIMnpip_Smmode_wbinoverlap[2];
  TH1D* IMnpim_Sp_wbin[nwbin][2];
  TH1D* IMnpip_Sm_wbin[nwbin][2];
  
  for(int iq=0;iq<2;iq++){
    cIMnpim_Spmode_wbinoverlap[iq] = new TCanvas(Form("cIMnpim_Spmode_wbinoverlap%d",iq),Form("cIMnpim_Spmode_wbinoverlap%d",iq));
    for(int iwbin=0;iwbin<nwbin;iwbin++){
      IMnpim_Sp_wbin[iwbin][iq] 
        = (TH1D*)IMnpim_IMnpip_dE_wK0_wSid_n_Sp_wbin[iwbin][iq+1]->ProjectionY(Form("IMnpim_Sp_wbin%d%d",iwbin,iq));
      IMnpim_Sp_wbin[iwbin][iq]->SetLineColor(iwbin);
    }
    IMnpim_Sp_wbin[1][iq]->SetTitle(Form("IMnpim K0orSp mode %d",iq));
    IMnpim_Sp_wbin[1][iq]->Draw("HE");
    IMnpim_Sp_wbin[2][iq]->Draw("HEsame");
    
    cIMnpip_Smmode_wbinoverlap[iq] = new TCanvas(Form("cIMnpip_Smmode_wbinoverlap%d",iq),Form("cIMnpip_Smmode_wbinoverlap%d",iq));
    for(int iwbin=0;iwbin<nwbin;iwbin++){
      IMnpip_Sm_wbin[iwbin][iq] 
        = (TH1D*)IMnpim_IMnpip_dE_wK0_wSid_n_Sm_wbin[iwbin][iq+1]->ProjectionX(Form("IMnpip_Sm_wbin%d%d",iwbin,iq));
      IMnpip_Sm_wbin[iwbin][iq]->SetLineColor(iwbin);
    }
    IMnpip_Sm_wbin[1][iq]->SetTitle(Form("IMnpip K0orSm mode %d",iq));
    IMnpip_Sm_wbin[1][iq]->Draw("HE");
    IMnpip_Sm_wbin[2][iq]->Draw("HEsame");
  }

  
  const float wbinlow[nwbin] = {1.0,1.40,1.52};
  const float wbinhigh[nwbin] = {1.40,1.52,2.00};
  const int spbinlow[2][2]={{21,39},
                            {21,38}};
  const int spbinhi[2][2]={{40,87},
                           {42,114}};
  for(int iq=0;iq<2;iq++){
    q_IMnpipi_K0orSp_ToSp[iq] = (TH2D*)q_IMnpipi_wK0_wSid_n_Sp[iq+1]->Clone(Form("q_IMnpipi_K0orSp_ToSp%d",iq));
    q_IMnpipi_K0orSp_ToK0[iq] = (TH2D*)q_IMnpipi_wK0_wSid_n_Sp[iq+1]->Clone(Form("q_IMnpipi_K0orSp_ToK0%d",iq));
    q_IMnpipi_K0orSp_ToSp[iq]->Reset();
    q_IMnpipi_K0orSp_ToK0[iq]->Reset();
    for(int iwbin=1;iwbin<3;iwbin++){//wbin=0 : no overlap region,so start from 1
      double nK0orSp = IMnpim_IMnpip_dE_wK0_wSid_n_Sp_wbin[iwbin][iq+1]->Integral();
      double nSp = 0.0;
      double nK0 = 0.0;
      if(nK0orSp>0){
        for(int ix=0;ix<IMnpim_IMnpip_dE_wK0_wSid_n_Sp_wbin[iwbin][iq+1]->GetNbinsX();ix++){
          for(int iy=spbinlow[iwbin-1][iq];iy<=spbinhi[iwbin-1][iq];iy++){
            double cont = IMnpim_IMnpip_dE_wK0_wSid_n_Sp_wbin[iwbin][iq+1]->GetBinContent(ix,iy);
            double nK0_bin = IMnpim_IMnpip_K0inter[iq]->GetBinContent(ix,iy);     
            if(cont<0.) nK0_bin = 0.0;
            double nSp_bin = cont-nK0_bin;
            if(nSp_bin<0.) nSp_bin=0.0;
            nSp += nSp_bin;
            nK0 += nK0_bin;
            if(cont>0.0 && iq==0){
              //std::cout << "cont " << cont << std::endl;
              //std::cout << "nK0_bin " << nK0_bin << std::endl;
              //std::cout << "nSp_bin " << nSp_bin << std::endl;
            }
          }
        }
      }//if nK0orSp
      //std::cout << "bin low "  << spbinlow[iwbin-1][iq] << std::endl;
      //std::cout << "bin high " <<  spbinhi[iwbin-1][iq] << std::endl;
      //std::cout << "nSp " << nSp << std::endl;
      //std::cout << "nK0 " << nK0 << std::endl;
      int wbinl = q_IMnpipi_K0orSp_ToSp[iq]->GetXaxis()->FindBin(wbinlow[iwbin]);
      int wbinh = q_IMnpipi_K0orSp_ToSp[iq]->GetXaxis()->FindBin(wbinhigh[iwbin]);
     
      double NevtWide=0.0;
      for(int ibin=wbinl;ibin<=wbinh;ibin++){
        for(int iqbin=0;iqbin<q_IMnpipi_wK0_wSid_n_Sp[iq+1]->GetNbinsY();iqbin++){
          double nevt = q_IMnpipi_wK0_wSid_n_Sp[iq+1]->GetBinContent(ibin,iqbin);
          //if(nevt>0.0){
            double ToSp = nevt*nSp/(nSp+nK0);
            double ToK0 = nevt*nK0/(nSp+nK0);
            if(iwbin==1 && iq==0){
            //if(isnan(ToSp)){
            //  std::cout << "iq: "   << iq << std::endl;
            //  std::cout << "ibin:"  << ibin << std::endl;
            //  std::cout << "iqbin:"  << iqbin << std::endl;
            //  std::cout << "nevt: " << nevt << std::endl;
            //  std::cout << "nK0orSp " << nK0orSp << std::endl;
            //  std::cout << "ToK0: " << nK0 << std::endl;
            }
              q_IMnpipi_K0orSp_ToSp[iq]->SetBinContent(ibin,iqbin,ToSp);
              q_IMnpipi_K0orSp_ToSp[iq]->SetBinError(ibin,iqbin,0);
              q_IMnpipi_K0orSp_ToK0[iq]->SetBinContent(ibin,iqbin,ToK0);
              q_IMnpipi_K0orSp_ToK0[iq]->SetBinError(ibin,iqbin,0);
           // }
          //}
          //if(nevt<=0.0 || nK0orSp <=0.0){ 
          //  q_IMnpipi_K0orSp_ToSp[iq]->SetBinContent(ibin,iqbin,0);
          //  q_IMnpipi_K0orSp_ToSp[iq]->SetBinError(ibin,iqbin,0);
          //  q_IMnpipi_K0orSp_ToK0[iq]->SetBinContent(ibin,iqbin,0);
          //  q_IMnpipi_K0orSp_ToK0[iq]->SetBinError(ibin,iqbin,0);
          //}
          //}
        }
      }//ibin
    }//iwbin
    IMnpipi_K0orSp_ToSp[iq] = (TH1D*)q_IMnpipi_K0orSp_ToSp[iq]->ProjectionX(Form("IMnpipi_K0orSp_ToSp%d",iq));
    IMnpipi_K0orSp_ToK0[iq] = (TH1D*)q_IMnpipi_K0orSp_ToK0[iq]->ProjectionX(Form("IMnpipi_K0orSp_ToK0%d",iq));
  }//iq
  


  TH2D* q_IMnpipi_K0orSm_ToSm[2];
  TH2D* q_IMnpipi_K0orSm_ToK0[2];
  TH1D* IMnpipi_K0orSm_ToSm[2];
  TH1D* IMnpipi_K0orSm_ToK0[2];
  const int smbinlow[2][2]={{16,34},
                            {17,34}};
  const int smbinhi[2][2]={{39,98},
                           {40,116}};
  for(int iq=0;iq<2;iq++){
    q_IMnpipi_K0orSm_ToSm[iq] = (TH2D*)q_IMnpipi_wK0_wSid_n_Sm[iq+1]->Clone(Form("q_IMnpipi_K0orSm_ToSm%d",iq));
    q_IMnpipi_K0orSm_ToK0[iq] = (TH2D*)q_IMnpipi_wK0_wSid_n_Sm[iq+1]->Clone(Form("q_IMnpipi_K0orSm_ToK0%d",iq));
    q_IMnpipi_K0orSm_ToSm[iq]->Reset();
    q_IMnpipi_K0orSm_ToK0[iq]->Reset();
    for(int iwbin=1;iwbin<3;iwbin++){//wbin=0 : no overlap region, so start from 1
      double nK0orSm = IMnpim_IMnpip_dE_wK0_wSid_n_Sm_wbin[iwbin][iq+1]->Integral();
      double nSm = 0.0;
      double nK0 = 0.0;
      if(nK0orSm>0){
        for(int ix=smbinlow[iwbin-1][iq];ix<=smbinhi[iwbin-1][iq];ix++){
          for(int iy=0;iy<IMnpim_IMnpip_dE_wK0_wSid_n_Sm_wbin[iwbin][iq+1]->GetNbinsY();iy++){
            double cont = IMnpim_IMnpip_dE_wK0_wSid_n_Sm_wbin[iwbin][iq+1]->GetBinContent(ix,iy);
            double nK0_bin = IMnpim_IMnpip_K0inter[iq]->GetBinContent(ix,iy);     
            if(cont<0.00001) nK0_bin = 0.0;
            double nSm_bin = cont-nK0_bin;
            if(nSm_bin<0.) nSm_bin=0.0;
            nSm += nSm_bin;
            nK0 += nK0_bin;
          }
        }
      }//if nK0orSm
      int wbinl = q_IMnpipi_K0orSm_ToSm[iq]->GetXaxis()->FindBin(wbinlow[iwbin]);
      int wbinh = q_IMnpipi_K0orSm_ToSm[iq]->GetXaxis()->FindBin(wbinhigh[iwbin]);
      double NevtWide=0.0;
      for(int ibin=wbinl;ibin<=wbinh;ibin++){
        for(int iqbin=0;iqbin<q_IMnpipi_wK0_wSid_n_Sm[iq+1]->GetNbinsY();iqbin++){
          double nevt = q_IMnpipi_wK0_wSid_n_Sm[iq+1]->GetBinContent(ibin,iqbin);
          //if(nevt>0.0){
            double ToSm = nevt*nSm/(nSm+nK0);
            double ToK0 = nevt*nK0/(nSm+nK0);
            //if(isnan(ToSm)){
           
            //}else{
              q_IMnpipi_K0orSm_ToSm[iq]->SetBinContent(ibin,iqbin,ToSm);
              q_IMnpipi_K0orSm_ToSm[iq]->SetBinError(ibin,iqbin,0);
              q_IMnpipi_K0orSm_ToK0[iq]->SetBinContent(ibin,iqbin,ToK0);
              q_IMnpipi_K0orSm_ToK0[iq]->SetBinError(ibin,iqbin,0);
           // }
          //}
          //if(nevt<=0.0 || nK0orSm <=0.0){ 
          //  q_IMnpipi_K0orSm_ToSm[iq]->SetBinContent(ibin,iqbin,0);
          //  q_IMnpipi_K0orSm_ToSm[iq]->SetBinError(ibin,iqbin,0);
          //  q_IMnpipi_K0orSm_ToK0[iq]->SetBinContent(ibin,iqbin,0);
          //  q_IMnpipi_K0orSm_ToK0[iq]->SetBinError(ibin,iqbin,0);
          //}
        }
      }//for iqbin
    }
    IMnpipi_K0orSm_ToSm[iq] = (TH1D*)q_IMnpipi_K0orSm_ToSm[iq]->ProjectionX(Form("IMnpipi_K0orSm_ToSm%d",iq));
    IMnpipi_K0orSm_ToK0[iq] = (TH1D*)q_IMnpipi_K0orSm_ToK0[iq]->ProjectionX(Form("IMnpipi_K0orSm_ToK0%d",iq));
  }


  
  TCanvas *cIMnpipi_wK0_Sp_afterDeco[2];
  for(int iq=0;iq<2;iq++){
    cIMnpipi_wK0_Sp_afterDeco[iq] = new TCanvas(Form("cIMnpipi_wK0_Sp_Deco_%s",cqcut[iq+1]),Form("cIMnpipi_wK0_Sp_Deco_%s",cqcut[iq+1]),800,800);
    IMnpipi_wK0_wSid_n_Sp[iq+1]->Draw("HE");
    IMnpipi_K0orSp_ToSp[iq]->SetLineColor(3);
    IMnpipi_K0orSp_ToSp[iq]->Draw("HEsame");
    IMnpipi_K0orSp_ToK0[iq]->SetLineColor(4);
    IMnpipi_K0orSp_ToK0[iq]->Draw("HEsame");
  }

  TCanvas *cIMnpipi_wK0_Sm_afterDeco[2];
  for(int iq=0;iq<2;iq++){
    cIMnpipi_wK0_Sm_afterDeco[iq] = new TCanvas(Form("cIMnpipi_wK0_Sm_Deco_%s",cqcut[iq+1]),Form("cIMnpipi_wK0_Sm_Deco_%s",cqcut[iq+1]),800,800);
    IMnpipi_wK0_wSid_n_Sm[iq+1]->Draw("HE");
    IMnpipi_K0orSm_ToSm[iq]->SetLineColor(3);
    IMnpipi_K0orSm_ToSm[iq]->Draw("HEsame");
    IMnpipi_K0orSm_ToK0[iq]->SetLineColor(4);
    IMnpipi_K0orSm_ToK0[iq]->Draw("HEsame");
  }
  

  TH2D* q_IMnpipi_SporSm_ToSp[2];
  TH2D* q_IMnpipi_SporSm_ToSm[2];
  TH1D* IMnpipi_SporSm_ToSp[2];
  TH1D* IMnpipi_SporSm_ToSm[2];
  for(int iq=0;iq<2;iq++){
    q_IMnpipi_SporSm_ToSp[iq] = (TH2D*)q_IMnpipi_wSid_n_SpSm[iq+1]->Clone(Form("q_IMnpipi_SporSm_ToSp%d",iq));
    q_IMnpipi_SporSm_ToSm[iq] = (TH2D*)q_IMnpipi_wSid_n_SpSm[iq+1]->Clone(Form("q_IMnpipi_SporSm_ToSm%d",iq));
    double nSporSm = IMnpim_IMnpip_dE_wSid_n_SpSm[iq+1]->Integral();
    double nSp = 0.;
    double nSm = 0.;
    if(nSporSm>0.0){
      for(int ix=0;ix<IMnpim_IMnpip_dE_wSid_n_SpSm[iq+1]->GetNbinsX();ix++){
        for(int iy=0;iy<IMnpim_IMnpip_dE_wSid_n_SpSm[iq+1]->GetNbinsY();iy++){
          double cont = IMnpim_IMnpip_dE_wSid_n_SpSm[iq+1]->GetBinContent(ix,iy);
          if(cont>0.0){
            double bincent_npip = IMnpim_IMnpip_dE_wSid_n_SpSm[iq+1]->GetXaxis()->GetBinCenter(ix);
            double bincent_npim = IMnpim_IMnpip_dE_wSid_n_SpSm[iq+1]->GetYaxis()->GetBinCenter(iy);
            double nSp_bin = gr_SpONnpim_fin_pol1[iq]->Eval(bincent_npim);
            double nSm_bin = gr_SmONnpip_fin_pol1[iq]->Eval(bincent_npip);
            //std::cout << "ibin:" << ibin << std::endl;
            //std::cout << "cont:" << cont << std::endl;
            //std::cout << "nSp " << nSp << std::endl;
            //std::cout << "nSp_bin" << nSp_bin << std::endl;
            //std::cout << "nSm_bin" << nSm_bin << std::endl;
            nSp += nSp_bin;
            nSm += nSm_bin;
            //std::cout<< "nSm_net_bin" << nSm_net_bin << std::endl;
            //std::cout<< "nSporSm" << nSporSm << std::endl;
          }//cont>0
        }//iy
      }//ix
    }//if nSporSm
    
    for(int ibin=0;ibin<q_IMnpipi_SporSm_ToSp[iq]->GetNbinsX();ibin++){
      for(int iqbin=0;iqbin<q_IMnpipi_SporSm_ToSp[iq]->GetNbinsY();iqbin++){
        double evt = q_IMnpipi_SporSm_ToSp[iq]->GetBinContent(ibin,iqbin);
        if((nSp+nSm)>0.0){
          double evtToSp = evt*nSp/(nSp+nSm);
          double evtToSm = evt*nSm/(nSp+nSm);
          q_IMnpipi_SporSm_ToSp[iq]->SetBinContent(ibin,iqbin,evtToSp);
          q_IMnpipi_SporSm_ToSp[iq]->SetBinError(ibin,iqbin,0);
          q_IMnpipi_SporSm_ToSm[iq]->SetBinContent(ibin,iqbin,evtToSm);
          q_IMnpipi_SporSm_ToSm[iq]->SetBinError(ibin,iqbin,0);
          if(nSp+nSm<0.00001){
            std::cout << "nSp + nSm " << nSp+nSm << std::endl;
            std::cout << "evt " << evt << std::endl;
          }
        }
      }
    }
    IMnpipi_SporSm_ToSp[iq] = (TH1D*)q_IMnpipi_SporSm_ToSp[iq]->ProjectionX(Form("IMnpipi_SporSm_ToSp%d",iq));
    IMnpipi_SporSm_ToSm[iq] = (TH1D*)q_IMnpipi_SporSm_ToSm[iq]->ProjectionX(Form("IMnpipi_SporSm_ToSm%d",iq));
  }
  
  TCanvas *cq_IMnpipi_Sp_afterDeco[2];
  TH2D* q_IMnpipi_Sp_sum[2];
  for(int iq=0;iq<2;iq++){
    q_IMnpipi_Sp_sum[iq] = (TH2D*)q_IMnpipi_wSid_n_Sp[iq+1]->Clone(Form("q_IMnpipi_Sp_sum%d",iq));
    q_IMnpipi_Sp_sum[iq]->Add(q_IMnpipi_K0orSp_ToK0[iq],-1.0);
    q_IMnpipi_Sp_sum[iq]->Add(q_IMnpipi_SporSm_ToSm[iq],-1.0);
    cq_IMnpipi_Sp_afterDeco[iq] = new TCanvas(Form("cq_IMnpipi_Sp_afterDeco%d",iq),Form("cq_IMnpipi_Sp_afterDeco%d",iq));
    cq_IMnpipi_Sp_afterDeco[iq]->cd();
    q_IMnpipi_Sp_sum[iq]->SetMinimum(0);
    q_IMnpipi_Sp_sum[iq]->Draw("colz");
  }

  TCanvas *cq_IMnpipi_Sm_afterDeco[2];
  TH2D* q_IMnpipi_Sm_sum[2];
  for(int iq=0;iq<2;iq++){
    q_IMnpipi_Sm_sum[iq] = (TH2D*)q_IMnpipi_wSid_n_Sm[iq+1]->Clone(Form("q_IMnpipi_Sm_sum%d",iq));
    q_IMnpipi_Sm_sum[iq]->Add(q_IMnpipi_K0orSm_ToK0[iq],-1.0);
    q_IMnpipi_Sm_sum[iq]->Add(q_IMnpipi_SporSm_ToSp[iq],-1.0);
    cq_IMnpipi_Sm_afterDeco[iq] = new TCanvas(Form("cq_IMnpipi_Sm_afterDeco%d",iq),Form("cq_IMnpipi_Sm_afterDeco%d",iq));
    cq_IMnpipi_Sm_afterDeco[iq]->cd();
    q_IMnpipi_Sm_sum[iq]->SetMinimum(0);
    q_IMnpipi_Sm_sum[iq]->Draw("colz");
  }


  TCanvas *cq_IMnpipi_K0_afterDeco[2];
  TH2D* q_IMnpipi_K0_sum[2];
  for(int iq=0;iq<2;iq++){
    q_IMnpipi_K0_sum[iq] = (TH2D*)q_IMnpipi_wK0_n[iq+1]->Clone(Form("q_IMnpipi_K0_sum%d",iq));
    q_IMnpipi_K0_sum[iq]->Add(q_IMnpipi_K0orSp_ToK0[iq],-1.0);
    q_IMnpipi_K0_sum[iq]->Add(q_IMnpipi_K0orSm_ToK0[iq],-1.0);
    cq_IMnpipi_K0_afterDeco[iq] = new TCanvas(Form("cq_IMnpipi_K0_afterDeco%d",iq),Form("cq_IMnpipi_K0_afterDeco%d",iq));
    cq_IMnpipi_K0_afterDeco[iq]->cd();
    q_IMnpipi_K0_sum[iq]->SetMinimum(0);
    q_IMnpipi_K0_sum[iq]->Draw("colz");
  }


  TCanvas *cIMnpipi_SpSm_afterDeco[2];
  for(int iq=0;iq<2;iq++){
    cIMnpipi_SpSm_afterDeco[iq] = new TCanvas(Form("cIMnpipi_SpSm_Deco_%s",cqcut[iq+1]),Form("cIMnpipi_SpSm_Deco_%s",cqcut[iq+1]),800,800);
    IMnpipi_wSid_n_SpSm[iq+1]->Draw("HE");
    IMnpipi_SporSm_ToSp[iq]->SetLineColor(4);
    IMnpipi_SporSm_ToSp[iq]->Draw("HEsame");
    IMnpipi_SporSm_ToSm[iq]->SetLineColor(5);
    IMnpipi_SporSm_ToSm[iq]->Draw("HEsame");
  }

  TCanvas *cIMnpipi_Summary_Sp[2];
  TH1D* IMnpipi_Sp_noK0_noSm[2];
  for(int iq=0;iq<2;iq++){
    cIMnpipi_Summary_Sp[iq] = new TCanvas(Form("cIMnpipi_Summary_Sp_%s",cqcut[iq+1]),Form("cIMnpipi_Summary_Sp_%s",cqcut[iq+1]),800,800);
    //all sum
    IMnpipi_wSid_n_Sp[iq+1]->Draw("HE");
    IMnpipi_Sp_noK0_noSm[iq] = (TH1D*)IMnpipi_wSid_n_Sp[iq+1]->Clone(Form("IMnpipi_Sp_noK0_noSm%d",iq));
    IMnpipi_Sp_noK0_noSm[iq]->Add(IMnpipi_K0orSp_ToK0[iq],-1);
    IMnpipi_Sp_noK0_noSm[iq]->Add(IMnpipi_SporSm_ToSm[iq],-1);
    IMnpipi_Sp_noK0_noSm[iq]->SetLineColor(2);
    IMnpipi_Sp_noK0_noSm[iq]->Draw("HEsame");
  }

  TCanvas *cIMnpipi_Summary_Sp_re[2];
  for(int iq=0;iq<2;iq++){
    cIMnpipi_Summary_Sp_re[iq] = new TCanvas(Form("cIMnpipi_Summary_Sp_re%s",cqcut[iq+1]),Form("cIMnpipi_Summary_Sp_re%s",cqcut[iq+1]),800,800);
    //all sum
    IMnpipi_Sp_noK0_noSm[iq]->Draw("HE");
  }

  TCanvas *cIMnpipi_Summary_Sm[2];
  TH1D* IMnpipi_Sm_noK0_noSp[2];
  for(int iq=0;iq<2;iq++){
    cIMnpipi_Summary_Sm[iq] = new TCanvas(Form("cIMnpipi_Summary_Sm_%s",cqcut[iq+1]),Form("cIMnpipi_Summary_Sm_%s",cqcut[iq+1]),800,800);
    //all sum
    IMnpipi_wSid_n_Sm[iq+1]->Draw("HE");
    IMnpipi_Sm_noK0_noSp[iq] = (TH1D*)IMnpipi_wSid_n_Sm[iq+1]->Clone(Form("IMnpipi_Sm_noK0_noSp%d",iq));
    IMnpipi_Sm_noK0_noSp[iq]->Add(IMnpipi_K0orSm_ToK0[iq],-1);
    IMnpipi_Sm_noK0_noSp[iq]->Add(IMnpipi_SporSm_ToSp[iq],-1);
    IMnpipi_Sm_noK0_noSp[iq]->SetLineColor(2);
    IMnpipi_Sm_noK0_noSp[iq]->Draw("HEsame");
  }

  TCanvas *cIMnpipi_Summary_Sm_re[2];
  for(int iq=0;iq<2;iq++){
    cIMnpipi_Summary_Sm_re[iq] = new TCanvas(Form("cIMnpipi_Summary_Sm_re%s",cqcut[iq+1]),Form("cIMnpipi_Summary_Sm_re%s",cqcut[iq+1]),800,800);
    //all sum
    IMnpipi_Sm_noK0_noSp[iq]->Draw("HE");
  }
   
  TFile *facc = TFile::Open("../simpost/accmap.root");
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


  TH2D* q_IMnpipi_Sp_cs[2];
  TH2D* q_IMnpipi_Sm_cs[2];
  TH2D* q_IMnpipi_K0_cs[2];
  TH2D* q_IMnpipi_Sp_cserr[2];
  TH2D* q_IMnpipi_Sm_cserr[2];
  TH2D* q_IMnpipi_K0_cserr[2];
  const int qcut650 =q_IMnpipi_Sp_sum[0]->GetYaxis()->FindBin(0.65);
  for(int iq=0;iq<2;iq++){
    q_IMnpipi_Sp_cs[iq] = (TH2D*)q_IMnpipi_Sp_sum[iq]->Clone(Form("q_IMnpipi_Sp_cs%d",iq));
    q_IMnpipi_Sm_cs[iq] = (TH2D*)q_IMnpipi_Sm_sum[iq]->Clone(Form("q_IMnpipi_Sm_cs%d",iq));
    q_IMnpipi_K0_cs[iq] = (TH2D*)q_IMnpipi_K0_sum[iq]->Clone(Form("q_IMnpipi_K0_cs%d",iq));
    for(int ix=0;ix<q_IMnpipi_Sp_cs[iq]->GetNbinsX();ix++){
      for(int iy=0;iy<q_IMnpipi_Sp_cs[iq]->GetNbinsY();iy++){
        double contSp =q_IMnpipi_Sp_cs[iq]->GetBinContent(ix,iy);
        double contSperr =q_IMnpipi_Sp_cs[iq]->GetBinError(ix,iy);
        double accSp =q_IMnpipi_Sp_accp->GetBinContent(ix,iy);
        double accerrSp =q_IMnpipi_Sp_accperr->GetBinContent(ix,iy);
        double csSp = 0.0; 
        double csSperr = 0.0; 
        if(accSp>0.0){
          csSp = contSp/accSp ;
          csSperr = contSperr/accSp;
        }
        double contSm =q_IMnpipi_Sm_cs[iq]->GetBinContent(ix,iy);
        double contSmerr =q_IMnpipi_Sm_cs[iq]->GetBinError(ix,iy);
        double accSm =q_IMnpipi_Sm_accp->GetBinContent(ix,iy);
        double accerrSm =q_IMnpipi_Sm_accperr->GetBinContent(ix,iy);
        double csSm = 0.0; 
        double csSmerr = 0.0; 
        if(accSm>0.0){
          csSm = contSm/accSm;
          csSmerr = contSmerr/accSm;
        }

        double contK0 =q_IMnpipi_K0_cs[iq]->GetBinContent(ix,iy);
        double contK0err =q_IMnpipi_K0_cs[iq]->GetBinError(ix,iy);
        double accK0 =q_IMnpipi_K0_accp->GetBinContent(ix,iy);
        double accerrK0 =q_IMnpipi_K0_accperr->GetBinContent(ix,iy);
        double csK0 = 0.0; 
        double csK0err = 0.0; 
        if(accK0>0.0){
          csK0 = contK0/accK0;
          csK0err = contK0err/accK0;
        }

        if(accerrSp<0.25){
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
        if(accerrSm<0.25){
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
        
        if(accerrK0<0.25){
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

  TCanvas *ccsSp[2];
  TCanvas *ccsSm[2];
  TCanvas *ccsK0[2];
  for(int iq=0;iq<2;iq++){
    ccsSp[iq] = new TCanvas(Form("ccsSp%d",iq),Form("ccsSp%d",iq),1600,800);
    ccsSp[iq]->Divide(2,1);
    ccsSp[iq]->cd(1);
    q_IMnpipi_Sp_cs[iq]->Draw("colz");
    ccsSp[iq]->cd(2);
    q_IMnpipi_Sp_cs[iq]->ProjectionX()->Draw("HE");
    ccsSm[iq] = new TCanvas(Form("ccsSm%d",iq),Form("ccsSm%d",iq),1600,800);
    ccsSm[iq]->Divide(2,1);
    ccsSm[iq]->cd(1);
    q_IMnpipi_Sm_cs[iq]->Draw("colz");
    ccsSm[iq]->cd(2);
    q_IMnpipi_Sm_cs[iq]->ProjectionX()->Draw("HE");
    ccsK0[iq] = new TCanvas(Form("ccsK0%d",iq),Form("ccsK0%d",iq),1600,800);
    ccsK0[iq]->Divide(2,1);
    ccsK0[iq]->cd(1);
    q_IMnpipi_K0_cs[iq]->Draw("colz");
    ccsK0[iq]->cd(2);
    q_IMnpipi_K0_cs[iq]->ProjectionX()->Draw("HE");
  }
  
  
  TCanvas *csum[2];
  TH2D* q_IMnpipi_SpSmSum[2];
  for(int iq=0;iq<2;iq++){
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
  TH2D* q_IMnpipi_SpSmSub[2];
  for(int iq=0;iq<2;iq++){
    q_IMnpipi_SpSmSub[iq] = (TH2D*)q_IMnpipi_Sp_cs[iq]->Clone(Form("q_IMnpipi_SpSmSub%d",iq));
    q_IMnpipi_SpSmSub[iq]->Add(q_IMnpipi_Sm_cs[iq],-1.0);
  }
  csub->cd(1);
  q_IMnpipi_SpSmSub[1]->SetTitle("q_IMnpipi_SpSmSub");
  q_IMnpipi_SpSmSub[1]->Draw("colz");
  csub->cd(2);
  q_IMnpipi_SpSmSub[1]->ProjectionX()->Draw("HE");
   
  TCanvas *ccomp[2];
  TH1D* IMnpipi_Sp_cs[2];
  TH1D* IMnpipi_Sm_cs[2];
  TH1D* IMnpipi_K0_cs[2];
  
  for(int iq=0;iq<2;iq++){
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
    IMnpipi_Sp_cs[iq]->Draw("HE");
    IMnpipi_Sm_cs[iq]->Draw("HEsame");
  }


  TCanvas *c = NULL;
  TSeqCollection *SCol = gROOT->GetListOfCanvases();
  int size = SCol->GetSize();
  TIter next(SCol);
  TString pdfname = "AfterDecompos.pdf";
  for(int i=0;i<size;i++){
    //pdf->NewPage();
    c= (TCanvas*)next();
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
  

  return;

}
