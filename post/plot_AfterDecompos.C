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
  
  fr[0] = TFile::Open("evanaIMpisigma_npippim_v202_out_iso_sub.root");
  fr[1] = TFile::Open("evanaIMpisigma_npippim_v202_out_iso_qlo_sub.root");
  fr[2] = TFile::Open("evanaIMpisigma_npippim_v202_out_iso_qhi_sub.root");
  fr[3] = TFile::Open("evanaIMpisigma_npippim_v202_out_iso_theta15_sub.root");
  fr[0]->Print();
  fr[1]->Print();
  fr[2]->Print();
  fr[3]->Print();
  TFile *fdeco[2]={NULL};
  fdeco[0] = TFile::Open("deco_qlo.root","READ");
  fdeco[1] = TFile::Open("deco_qhi.root","READ");
  fdeco[0]->Print();
  fdeco[1]->Print();

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
    q_IMnpipi_wSid_n_Sp[iq] = (TH2D*)fr[iq]->Get("q_IMnpipi_wSid_n_Sp");
    q_IMnpipi_wK0_wSid_n_Sp[iq] = (TH2D*)fr[iq]->Get("q_IMnpipi_wK0_wSid_n_Sp");
    q_IMnpipi_wSid_n_Sm[iq] = (TH2D*)fr[iq]->Get("q_IMnpipi_wSid_n_Sm");
    q_IMnpipi_wK0_wSid_n_Sm[iq] = (TH2D*)fr[iq]->Get("q_IMnpipi_wK0_wSid_n_Sm");
    q_IMnpipi_wSid_n_SpSm[iq] = (TH2D*)fr[iq]->Get("q_IMnpipi_wSid_n_SpSm");
    q_IMnpipi_wK0_n[iq] = (TH2D*)fr[iq]->Get("q_IMnpipi_wK0_n");
    
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

    //Sigma+ & Sigma- & K0 overlap
    for(unsigned int ibin=0;ibin<nbinIMnpipi;ibin++){
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

  return;
  
  const int spbinlow[2][2]={{7,11},
                            {7,11}};
  const int spbinhi[2][2]={{12,23},
                           {12,30}};

  for(int iq=0;iq<2;iq++){
    q_IMnpipi_K0orSp_ToSp[iq] = (TH2D*)q_IMnpipi_wK0_wSid_n_Sp[iq+1]->Clone(Form("q_IMnpipi_K0orSp_ToSp%d",iq));
    q_IMnpipi_K0orSp_ToK0[iq] = (TH2D*)q_IMnpipi_wK0_wSid_n_Sp[iq+1]->Clone(Form("q_IMnpipi_K0orSp_ToK0%d",iq));
    //IMnpipi_K0orSp_ToSp[iq] = (TH1D*)IMnpipi_wK0_wSid_n_Sp[iq+1]->Clone(Form("IMnpipi_K0orSp_ToSp%d",iq));
    //IMnpipi_K0orSp_ToK0[iq] = (TH1D*)IMnpipi_wK0_wSid_n_Sp[iq+1]->Clone(Form("IMnpipi_K0orSp_ToK0%d",iq));
    //for(int ibin=0;ibin<nbinIMnpipi;ibin++){
    for(int ibin=0;ibin<2;ibin++){
      //double nK0orSp = IMnpim_IMnpip_dE_wK0_wSid_n_Sp_bin[ibin][iq+1]->Integral();
      double nK0orSp = IMnpim_IMnpip_dE_wK0_wSid_n_Sp_wbin[ibin][iq+1]->Integral();
      double nSp = 0.0;
      double nK0 = 0.0;
      if(nK0orSp>0){
        for(int ix=0;ix<IMnpim_IMnpip_dE_wK0_wSid_n_Sp_wbin[ibin][iq+1]->GetNbinsX();ix++){
          for(int iy=spbinlow[ibin][iq];iy<IMnpim_IMnpip_dE_wK0_wSid_n_Sp_wbin[ibin][iq+1]->GetNbinsY();iy++){
            double cont = IMnpim_IMnpip_dE_wK0_wSid_n_Sp_wbin[ibin][iq+1]->GetBinContent(ix,iy);
            double nK0_bin = IMnpim_IMnpip_K0inter[iq]->GetBinContent(ix,iy);     
            double nSp_bin = cont-nK0_bin;
            if(nSp_bin<0.) nSp_bin=0.0;
            nSp += nSp_bin;
            nK0 += nK0_bin;
          }
        }
      }//if nK0orSp
      for(int iqbin=0;iqbin<q_IMnpipi_K0orSp_ToSp[iq]->GetNbinsY();iqbin++){
        double nevt = q_IMnpipi_K0orSp_ToSp[iq]->GetBinContent(ibin+1,iqbin);
        if(nevt>0.0){
          double ToSp = nevt*nSp/(nSp+nK0);
          double ToK0 = nevt*nK0/(nSp+nK0);
          if(isnan(ToSp)){
            std::cout << "iq: "   << iq << std::endl;
            std::cout << "ibin:"  << ibin << std::endl;
            std::cout << "iqbin:"  << iqbin << std::endl;
            std::cout << "nevt: " << nevt << std::endl;
            std::cout << "nK0orSp " << nK0orSp << std::endl;
            std::cout << "ToK0: " << nK0 << std::endl;
          }else{
            q_IMnpipi_K0orSp_ToSp[iq]->SetBinContent(ibin+1,iqbin,ToSp);
            q_IMnpipi_K0orSp_ToSp[iq]->SetBinError(ibin+1,iqbin,0);
            q_IMnpipi_K0orSp_ToK0[iq]->SetBinContent(ibin+1,iqbin,ToK0);
            q_IMnpipi_K0orSp_ToK0[iq]->SetBinError(ibin+1,iqbin,0);
          }
        }//if nevt>0.0
        if(nevt<=0.0 || nK0orSp <=0.0){ 
          q_IMnpipi_K0orSp_ToSp[iq]->SetBinContent(ibin+1,iqbin,0);
          q_IMnpipi_K0orSp_ToSp[iq]->SetBinError(ibin+1,iqbin,0);
          q_IMnpipi_K0orSp_ToK0[iq]->SetBinContent(ibin+1,iqbin,0);
          q_IMnpipi_K0orSp_ToK0[iq]->SetBinError(ibin+1,iqbin,0);
        }
      }
      //IMnpipi_K0orSp_ToSp[iq]->SetBinContent(ibin+1,nK0orSp);
      //IMnpipi_K0orSp_ToSp[iq]->SetBinError(ibin+1,0);
      //IMnpipi_K0orSp_ToK0[iq]->SetBinContent(ibin+1,nSp - nK0orSp);
      //IMnpipi_K0orSp_ToK0[iq]->SetBinError(ibin+1,0);
      IMnpipi_K0orSp_ToSp[iq] = (TH1D*)q_IMnpipi_K0orSp_ToSp[iq]->ProjectionX(Form("IMnpipi_K0orSp_ToSp%d",iq));
      IMnpipi_K0orSp_ToK0[iq] = (TH1D*)q_IMnpipi_K0orSp_ToK0[iq]->ProjectionX(Form("IMnpipi_K0orSp_ToK0%d",iq));
    }//for ibin
  }


  TH2D* q_IMnpipi_K0orSm_ToSm[2];
  TH2D* q_IMnpipi_K0orSm_ToK0[2];
  TH1D* IMnpipi_K0orSm_ToSm[2];
  TH1D* IMnpipi_K0orSm_ToK0[2];
  const int smbinlow[2][2]={{7,11},
                            {7,11}};
  const int smbinhi[2][2]={{12,28},
                           {12,30}};
  for(int iq=0;iq<2;iq++){
    q_IMnpipi_K0orSm_ToSm[iq] = (TH2D*)q_IMnpipi_wK0_wSid_n_Sm[iq+1]->Clone(Form("q_IMnpipi_K0orSm_ToSm%d",iq));
    q_IMnpipi_K0orSm_ToK0[iq] = (TH2D*)q_IMnpipi_wK0_wSid_n_Sm[iq+1]->Clone(Form("q_IMnpipi_K0orSm_ToK0%d",iq));
    //IMnpipi_K0orSm_ToSm[iq] = (TH1D*)IMnpipi_wK0_wSid_n_Sm[iq+1]->Clone(Form("IMnpipi_K0orSm_ToSm%d",iq));
    //IMnpipi_K0orSm_ToK0[iq] = (TH1D*)IMnpipi_wK0_wSid_n_Sm[iq+1]->Clone(Form("IMnpipi_K0orSm_ToK0%d",iq));
    for(int ibin=0;ibin<nbinIMnpipi;ibin++){
      double nK0orSm = IMnpim_IMnpip_dE_wK0_wSid_n_Sm_bin[ibin][iq+1]->Integral();
      double nSm = 0.0;
      double nK0 = 0.0;
      if(nK0orSm>0){
        for(int ix=0;ix<IMnpim_IMnpip_dE_wK0_wSid_n_Sm_bin[ibin][iq+1]->GetNbinsX();ix++){
          for(int iy=0;iy<IMnpim_IMnpip_dE_wK0_wSid_n_Sm_bin[ibin][iq+1]->GetNbinsY();iy++){
            double cont = IMnpim_IMnpip_dE_wK0_wSid_n_Sm_bin[ibin][iq+1]->GetBinContent(ix,iy);
            double nK0_bin = IMnpim_IMnpip_K0inter[iq]->GetBinContent(ix,iy);     
            double nSm_bin = cont-nK0_bin;
            if(nSm_bin<0.) nSm_bin=0.0;
            nSm += nSm_bin;
            nK0 += nK0_bin;
          }
        }
      }//if nK0orSm
      for(int iqbin=0;iqbin<q_IMnpipi_K0orSm_ToSm[iq]->GetNbinsY();iqbin++){
        double nevt = q_IMnpipi_K0orSm_ToSm[iq]->GetBinContent(ibin+1,iqbin);
        //if(nevt>0.0){
          double ToSm = nevt*nSm/(nSm+nK0);
          double ToK0 = nevt*nK0/(nSm+nK0);
          if(isnan(ToSm)){
           
          }else{
            q_IMnpipi_K0orSm_ToSm[iq]->SetBinContent(ibin+1,iqbin,ToSm);
            q_IMnpipi_K0orSm_ToSm[iq]->SetBinError(ibin+1,iqbin,0);
            q_IMnpipi_K0orSm_ToK0[iq]->SetBinContent(ibin+1,iqbin,ToK0);
            q_IMnpipi_K0orSm_ToK0[iq]->SetBinError(ibin+1,iqbin,0);
          }
        //}
        //if(nevt<=0.0 || nK0orSm <=0.0){ 
        if( nK0orSm <=0.0){ 
          q_IMnpipi_K0orSm_ToSm[iq]->SetBinContent(ibin+1,iqbin,0);
          q_IMnpipi_K0orSm_ToSm[iq]->SetBinError(ibin+1,iqbin,0);
          q_IMnpipi_K0orSm_ToK0[iq]->SetBinContent(ibin+1,iqbin,0);
          q_IMnpipi_K0orSm_ToK0[iq]->SetBinError(ibin+1,iqbin,0);
        }
      }//for iqbin
      IMnpipi_K0orSm_ToSm[iq] = (TH1D*)q_IMnpipi_K0orSm_ToSm[iq]->ProjectionX(Form("IMnpipi_K0orSm_ToSm%d",iq));
      IMnpipi_K0orSm_ToK0[iq] = (TH1D*)q_IMnpipi_K0orSm_ToK0[iq]->ProjectionX(Form("IMnpipi_K0orSm_ToK0%d",iq));
    }
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

  TH1D* IMnpipi_SporSm_ToSp[2];
  TH1D* IMnpipi_SporSm_ToSm[2];
  for(int iq=0;iq<2;iq++){
    IMnpipi_SporSm_ToSp[iq] = (TH1D*)IMnpipi_wSid_n_SpSm[iq+1]->Clone(Form("IMnpipi_SporSm_ToSp%d",iq));
    IMnpipi_SporSm_ToSm[iq] = (TH1D*)IMnpipi_wSid_n_SpSm[iq+1]->Clone(Form("IMnpipi_SporSm_ToSm%d",iq));
    for(int ibin=0;ibin<nbinIMnpipi;ibin++){
      double nSporSm = IMnpim_IMnpip_dE_wSid_n_SpSm_bin[ibin][iq+1]->Integral();
      double nSp = nSporSm;
      if(nSporSm>0.0){
        for(int ix=0;ix<IMnpim_IMnpip_dE_wSid_n_SpSm_bin[ibin][iq+1]->GetNbinsX();ix++){
          for(int iy=0;iy<IMnpim_IMnpip_dE_wSid_n_SpSm_bin[ibin][iq+1]->GetNbinsY();iy++){
            double cont = IMnpim_IMnpip_dE_wSid_n_SpSm_bin[ibin][iq+1]->GetBinContent(ix,iy);
            if(cont>0.0){
              double bincent_npip = IMnpim_IMnpip_dE_wSid_n_SpSm_bin[ibin][iq+1]->GetXaxis()->GetBinCenter(ix);
              double bincent_npim = IMnpim_IMnpip_dE_wSid_n_SpSm_bin[ibin][iq+1]->GetYaxis()->GetBinCenter(iy);
              double nSp_bin = gr_SpONnpim_fin_pol1[iq]->Eval(bincent_npim);
              double nSm_bin = gr_SmONnpip_fin_pol1[iq]->Eval(bincent_npip);
              //std::cout << "ibin:" << ibin << std::endl;
              //std::cout << "cont:" << cont << std::endl;
              //std::cout << "nSp " << nSp << std::endl;
              //std::cout << "nSp_bin" << nSp_bin << std::endl;
              //std::cout << "nSm_bin" << nSm_bin << std::endl;
              double nSm_net_bin = cont*nSm_bin/(nSp_bin+nSm_bin);
              nSporSm -=nSm_net_bin;
              //std::cout<< "nSm_net_bin" << nSm_net_bin << std::endl;
              //std::cout<< "nSporSm" << nSporSm << std::endl;
            }
          }
        }
      }//if nSporSm
      IMnpipi_SporSm_ToSp[iq]->SetBinContent(ibin+1,nSporSm);
      IMnpipi_SporSm_ToSp[iq]->SetBinError(ibin+1,0);
      IMnpipi_SporSm_ToSm[iq]->SetBinContent(ibin+1,nSp-nSporSm);
      IMnpipi_SporSm_ToSm[iq]->SetBinError(ibin+1,0);
    }//for ibin
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
