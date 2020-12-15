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

void HistToRorateGraph(TH1D* h1, TGraphErrors &gr)
{
  for(int ibin=1;ibin<=h1->GetNbinsX();ibin++){
    double cont = h1->GetBinContent(ibin);
    double err  = h1->GetBinError(ibin);
    double bincenter = h1->GetBinCenter(ibin);
    gr.SetPoint(ibin,cont,bincenter);
    gr.SetPointError(ibin,err,0);
    if(cont < -1000) std::cout << ibin << " " << cont << std::endl;
  }
  gr.GetYaxis()->SetRangeUser(h1->GetXaxis()->GetXmin(),h1->GetXaxis()->GetXmax());
  gr.GetYaxis()->SetLabelSize(0);
  gr.SetMarkerStyle(20);
  return;
}

void K0SigmaTemp()
{
  TFile *fr[3]={NULL};
  TFile *fmix[3]={NULL};
  
  
  fr[0] = TFile::Open("evanaIMpisigma_npippim_v202_out_iso.root");
  fmix[1] = TFile::Open("evanaIMpisigma_npippim_v202_MIX_cut4_out_iso.root");
  fr[1] = TFile::Open("evanaIMpisigma_npippim_v202_out_iso_qlo.root");
  fmix[1] = TFile::Open("evanaIMpisigma_npippim_v202_MIX_cut4_out_iso_qlo.root");
  fr[2] = TFile::Open("evanaIMpisigma_npippim_v202_out_iso_qhi.root");
  fmix[2] = TFile::Open("evanaIMpisigma_npippim_v202_MIX_cut4_out_iso_qhi.root");
  fr[0]->Print() ;
  fmix[0]->Print();
  fr[1]->Print() ;
  fmix[1]->Print();
  fr[2]->Print() ;
  fmix[2]->Print();
  
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0.);  

  const unsigned int nbintemplate = 50;
  TH2F* IMpippim_IMnpip_n_data[3];
  TH2F* IMpippim_IMnpim_n_data[3];
  TH2F* IMpippim_IMnpip_n_bin_data[nbintemplate][3];
  TH2F* IMpippim_IMnpim_n_bin_data[nbintemplate][3];
  TH2F* IMpippim_IMnpip_n_mix[3];
  TH2F* IMpippim_IMnpim_n_mix[3];
  TH2F* IMpippim_IMnpip_n_bin_mix[nbintemplate][3];
  TH2F* IMpippim_IMnpim_n_bin_mix[nbintemplate][3];
  TH2F* IMpippim_IMnpip_n_sub[3];
  TH2F* IMpippim_IMnpim_n_sub[3];
  TH2F* IMpippim_IMnpip_n_bin_sub[nbintemplate][3];
  TH2F* IMpippim_IMnpim_n_bin_sub[nbintemplate][3];
  TH2F* IMpippim_IMnpip_n_woSm_data[3];
  TH2F* IMpippim_IMnpim_n_woSp_data[3];
  TH2F* IMpippim_IMnpip_n_woSm_bin_data[nbintemplate][3];
  TH2F* IMpippim_IMnpim_n_woSp_bin_data[nbintemplate][3];
  TH2F* IMpippim_IMnpip_n_woSm_mix[3];
  TH2F* IMpippim_IMnpim_n_woSp_mix[3];
  TH2F* IMpippim_IMnpip_n_woSm_bin_mix[nbintemplate][3];
  TH2F* IMpippim_IMnpim_n_woSp_bin_mix[nbintemplate][3];
  TH2F* IMpippim_IMnpip_n_woSm_sub[3];
  TH2F* IMpippim_IMnpim_n_woSp_sub[3];
  TH2F* IMpippim_IMnpip_n_woSm_bin_sub[nbintemplate][3];
  TH2F* IMpippim_IMnpim_n_woSp_bin_sub[nbintemplate][3];
  TH2F* IMpippim_IMnpip_n_woSmdia_data[3];
  TH2F* IMpippim_IMnpim_n_woSpdia_data[3];
  TH2F* IMpippim_IMnpip_n_woSmdia_bin_data[nbintemplate][3];
  TH2F* IMpippim_IMnpim_n_woSpdia_bin_data[nbintemplate][3];
  TH2F* IMpippim_IMnpip_n_woSmdia_mix[3];
  TH2F* IMpippim_IMnpim_n_woSpdia_mix[3];
  TH2F* IMpippim_IMnpip_n_woSmdia_bin_mix[nbintemplate][3];
  TH2F* IMpippim_IMnpim_n_woSpdia_bin_mix[nbintemplate][3];
  TH2F* IMpippim_IMnpip_n_woSmdia_sub[3];
  TH2F* IMpippim_IMnpim_n_woSpdia_sub[3];
  TH2F* IMpippim_IMnpip_n_woSmdia_bin_sub[nbintemplate][3];
  TH2F* IMpippim_IMnpim_n_woSpdia_bin_sub[nbintemplate][3];
  
  const char cqcut[][6]= {"qlo","qhi"};
  for(int iqcut=0;iqcut<3;iqcut++){
    IMpippim_IMnpip_n_data[iqcut] = (TH2F*)fr[iqcut]->Get("IMpippim_IMnpip_n");
    IMpippim_IMnpim_n_data[iqcut] = (TH2F*)fr[iqcut]->Get("IMpippim_IMnpim_n");
    IMpippim_IMnpip_n_mix[iqcut] = (TH2F*)fmix[iqcut]->Get("IMpippim_IMnpip_n");
    IMpippim_IMnpim_n_mix[iqcut] = (TH2F*)fmix[iqcut]->Get("IMpippim_IMnpim_n");
    IMpippim_IMnpip_n_sub[iqcut] = (TH2F*)IMpippim_IMnpip_n_data[iqcut]->Clone(Form("IMpippim_IMnpip_n_sub_%s",cqcut[iqcut]));
    IMpippim_IMnpip_n_sub[iqcut]->Add(IMpippim_IMnpip_n_mix[iqcut],-1.0); 
    IMpippim_IMnpip_n_sub[iqcut]->SetTitle(Form("IMpippim_IMnpip_n_sub_%s",cqcut[iqcut]));
    IMpippim_IMnpim_n_sub[iqcut] = (TH2F*)IMpippim_IMnpim_n_data[iqcut]->Clone(Form("IMpippim_IMnpim_n_sub_%s",cqcut[iqcut]));
    IMpippim_IMnpim_n_sub[iqcut]->Add(IMpippim_IMnpim_n_mix[iqcut],-1.0); 
    IMpippim_IMnpim_n_sub[iqcut]->SetTitle(Form("IMpippim_IMnpim_n_sub_%s",cqcut[iqcut]));
    
    IMpippim_IMnpip_n_woSm_data[iqcut] = (TH2F*)fr[iqcut]->Get("IMpippim_IMnpip_n_woSm");
    IMpippim_IMnpim_n_woSp_data[iqcut] = (TH2F*)fr[iqcut]->Get("IMpippim_IMnpim_n_woSp");
    IMpippim_IMnpip_n_woSm_mix[iqcut] = (TH2F*)fmix[iqcut]->Get("IMpippim_IMnpip_n_woSm");
    IMpippim_IMnpim_n_woSp_mix[iqcut] = (TH2F*)fmix[iqcut]->Get("IMpippim_IMnpim_n_woSp");
    IMpippim_IMnpip_n_woSm_sub[iqcut] = (TH2F*)IMpippim_IMnpip_n_woSm_data[iqcut]->Clone(Form("IMpippim_IMnpip_n_woSm_sub_%s",cqcut[iqcut]));
    IMpippim_IMnpip_n_woSm_sub[iqcut]->Add(IMpippim_IMnpip_n_woSm_mix[iqcut],-1.0); 
    IMpippim_IMnpip_n_woSm_sub[iqcut]->SetTitle(Form("IMpippim_IMnpip_n_woSm_sub_%s",cqcut[iqcut]));
    IMpippim_IMnpim_n_woSp_sub[iqcut] = (TH2F*)IMpippim_IMnpim_n_woSp_data[iqcut]->Clone(Form("IMpippim_IMnpim_n_woSp_sub_%s",cqcut[iqcut]));
    IMpippim_IMnpim_n_woSp_sub[iqcut]->Add(IMpippim_IMnpim_n_woSp_mix[iqcut],-1.0); 
    IMpippim_IMnpim_n_woSp_sub[iqcut]->SetTitle(Form("IMpippim_IMnpim_n_woSp_sub_%s",cqcut[iqcut]));
    
    IMpippim_IMnpip_n_woSmdia_data[iqcut] = (TH2F*)fr[iqcut]->Get("IMpippim_IMnpip_n_woSmdia");
    IMpippim_IMnpim_n_woSpdia_data[iqcut] = (TH2F*)fr[iqcut]->Get("IMpippim_IMnpim_n_woSpdia");
    IMpippim_IMnpip_n_woSmdia_mix[iqcut] = (TH2F*)fmix[iqcut]->Get("IMpippim_IMnpip_n_woSmdia");
    IMpippim_IMnpim_n_woSpdia_mix[iqcut] = (TH2F*)fmix[iqcut]->Get("IMpippim_IMnpim_n_woSpdia");
    IMpippim_IMnpip_n_woSmdia_sub[iqcut] = (TH2F*)IMpippim_IMnpip_n_woSmdia_data[iqcut]->Clone(Form("IMpippim_IMnpip_n_woSmdia_sub_%s",cqcut[iqcut]));
    IMpippim_IMnpip_n_woSmdia_sub[iqcut]->Add(IMpippim_IMnpip_n_woSmdia_mix[iqcut],-1.0); 
    IMpippim_IMnpip_n_woSmdia_sub[iqcut]->SetTitle(Form("IMpippim_IMnpip_n_woSmdia_sub_%s",cqcut[iqcut]));
    IMpippim_IMnpim_n_woSpdia_sub[iqcut] = (TH2F*)IMpippim_IMnpim_n_woSpdia_data[iqcut]->Clone(Form("IMpippim_IMnpim_n_woSpdia_sub_%s",cqcut[iqcut]));
    IMpippim_IMnpim_n_woSpdia_sub[iqcut]->Add(IMpippim_IMnpim_n_woSpdia_mix[iqcut],-1.0); 
    IMpippim_IMnpim_n_woSpdia_sub[iqcut]->SetTitle(Form("IMpippim_IMnpim_n_woSpdia_sub_%s",cqcut[iqcut]));
  }
  

  for(unsigned int ibin=0;ibin<nbintemplate;ibin++){
    IMpippim_IMnpip_n_bin_data[ibin] = (TH2F*)fr->Get(Form("IMpippim_IMnpip_n_bin%d",ibin));
    IMpippim_IMnpim_n_bin_data[ibin] = (TH2F*)fr->Get(Form("IMpippim_IMnpim_n_bin%d",ibin));
    IMpippim_IMnpip_n_bin_mix[ibin] = (TH2F*)fmix->Get(Form("IMpippim_IMnpip_n_bin%d",ibin));
    IMpippim_IMnpim_n_bin_mix[ibin] = (TH2F*)fmix->Get(Form("IMpippim_IMnpim_n_bin%d",ibin));
    IMpippim_IMnpip_n_bin_sub[ibin] = (TH2F*)IMpippim_IMnpip_n_bin_data[ibin]->Clone(Form("IMpippim_IMnpip_n_sub_bin%d",ibin));
    IMpippim_IMnpip_n_bin_sub[ibin]->Add(IMpippim_IMnpip_n_bin_mix[ibin],-1.0);
    IMpippim_IMnpim_n_bin_sub[ibin] = (TH2F*)IMpippim_IMnpim_n_bin_data[ibin]->Clone(Form("IMpippim_IMnpim_n_sub_bin%d",ibin));
    IMpippim_IMnpim_n_bin_sub[ibin]->Add(IMpippim_IMnpim_n_bin_mix[ibin],-1.0);
    
    IMpippim_IMnpip_n_woSm_bin_data[ibin] = (TH2F*)fr->Get(Form("IMpippim_IMnpip_n_woSm_bin%d",ibin));
    IMpippim_IMnpim_n_woSp_bin_data[ibin] = (TH2F*)fr->Get(Form("IMpippim_IMnpim_n_woSp_bin%d",ibin));
    IMpippim_IMnpip_n_woSm_bin_mix[ibin] = (TH2F*)fmix->Get(Form("IMpippim_IMnpip_n_woSm_bin%d",ibin));
    IMpippim_IMnpim_n_woSp_bin_mix[ibin] = (TH2F*)fmix->Get(Form("IMpippim_IMnpim_n_woSp_bin%d",ibin));
    IMpippim_IMnpip_n_woSm_bin_sub[ibin] = (TH2F*)IMpippim_IMnpip_n_woSm_bin_data[ibin]->Clone(Form("IMpippim_IMnpip_n_woSm_sub_bin%d",ibin));
    IMpippim_IMnpip_n_woSm_bin_sub[ibin]->Add(IMpippim_IMnpip_n_woSm_bin_mix[ibin],-1.0);
    IMpippim_IMnpim_n_woSp_bin_sub[ibin] = (TH2F*)IMpippim_IMnpim_n_woSp_bin_data[ibin]->Clone(Form("IMpippim_IMnpim_n_woSp_sub_bin%d",ibin));
    IMpippim_IMnpim_n_woSp_bin_sub[ibin]->Add(IMpippim_IMnpim_n_woSp_bin_mix[ibin],-1.0);

    IMpippim_IMnpip_n_woSmdia_bin_data[ibin] = (TH2F*)fr->Get(Form("IMpippim_IMnpip_n_woSmdia_bin%d",ibin));
    IMpippim_IMnpim_n_woSpdia_bin_data[ibin] = (TH2F*)fr->Get(Form("IMpippim_IMnpim_n_woSpdia_bin%d",ibin));
    IMpippim_IMnpip_n_woSmdia_bin_mix[ibin] = (TH2F*)fmix->Get(Form("IMpippim_IMnpip_n_woSmdia_bin%d",ibin));
    IMpippim_IMnpim_n_woSpdia_bin_mix[ibin] = (TH2F*)fmix->Get(Form("IMpippim_IMnpim_n_woSpdia_bin%d",ibin));
    IMpippim_IMnpip_n_woSmdia_bin_sub[ibin] = (TH2F*)IMpippim_IMnpip_n_woSmdia_bin_data[ibin]->Clone(Form("IMpippim_IMnpip_n_woSmdia_sub_bin%d",ibin));
    IMpippim_IMnpip_n_woSmdia_bin_sub[ibin]->Add(IMpippim_IMnpip_n_woSmdia_bin_mix[ibin],-1.0);
    IMpippim_IMnpim_n_woSpdia_bin_sub[ibin] = (TH2F*)IMpippim_IMnpim_n_woSpdia_bin_data[ibin]->Clone(Form("IMpippim_IMnpim_n_woSpdia_sub_bin%d",ibin));
    IMpippim_IMnpim_n_woSpdia_bin_sub[ibin]->Add(IMpippim_IMnpim_n_woSpdia_bin_mix[ibin],-1.0);
  }
  
  TCanvas *cIMpippim_IMnpip_n_sub = new TCanvas("cIMpippim_IMnpip_n_sub","cIMpippim_IMnpip_n_sub",800,800);
  cIMpippim_IMnpip_n_sub->Divide(2,2,0,0);
  cIMpippim_IMnpip_n_sub->cd(3);
  IMpippim_IMnpip_n_sub->RebinY(2);
  IMpippim_IMnpip_n_sub->Draw("colz");
  cIMpippim_IMnpip_n_sub->cd(1);
  TH1D* IMnpip_n_sub = (TH1D*)IMpippim_IMnpip_n_sub->ProjectionX("IMnpip_n_sub");
  IMnpip_n_sub->Draw("HE");
  cIMpippim_IMnpip_n_sub->cd(4);
  TH1D* IMpippim_n_sub = (TH1D*)IMpippim_IMnpip_n_sub->ProjectionY("IMpippim_n_sub");
  TGraphErrors *gr_IMpippim_n_sub = new TGraphErrors();
  HistToRorateGraph(IMpippim_n_sub,*gr_IMpippim_n_sub);
  gr_IMpippim_n_sub->Draw("AP");
  
  TCanvas *cIMpippim_IMnpip_n_woSmdia_sub = new TCanvas("cIMpippim_IMnpip_n_woSmdia_sub","cIMpippim_IMnpip_n_woSmdia_sub",800,800);
  cIMpippim_IMnpip_n_woSmdia_sub->Divide(2,2,0,0);
  cIMpippim_IMnpip_n_woSmdia_sub->cd(3);
  IMpippim_IMnpip_n_woSmdia_sub->RebinY(2);
  IMpippim_IMnpip_n_woSmdia_sub->Draw("colz");
  cIMpippim_IMnpip_n_woSmdia_sub->cd(1);
  TH1D* IMnpip_n_woSmdia_sub = (TH1D*)IMpippim_IMnpip_n_woSmdia_sub->ProjectionX("IMnpip_n_woSmdia_sub");
  IMnpip_n_woSmdia_sub->Draw("HE");
  cIMpippim_IMnpip_n_woSmdia_sub->cd(4);
  TH1D* IMpippim_n_woSmdia_sub = (TH1D*)IMpippim_IMnpip_n_woSmdia_sub->ProjectionY("IMpippim_n_woSmdia_sub");
  TGraphErrors *gr_IMpippim_n_woSmdia_sub = new TGraphErrors();
  HistToRorateGraph(IMpippim_n_woSmdia_sub,*gr_IMpippim_n_woSmdia_sub);
  gr_IMpippim_n_woSmdia_sub->Draw("AP");
  
  TCanvas *cIMpippim_IMnpim_n_sub = new TCanvas("cIMpippim_IMnpim_n_sub","cIMpippim_IMnpim_n_sub",800,800);
  cIMpippim_IMnpim_n_sub->Divide(2,2,0,0);
  cIMpippim_IMnpim_n_sub->cd(3);
  IMpippim_IMnpim_n_sub->RebinY(2);
  IMpippim_IMnpim_n_sub->Draw("colz");
  cIMpippim_IMnpim_n_sub->cd(1);
  TH1D* IMnpim_n_sub = (TH1D*)IMpippim_IMnpim_n_sub->ProjectionX("IMnpim_n_sub");
  IMnpim_n_sub->Draw("HE");
  cIMpippim_IMnpim_n_sub->cd(4);
  gr_IMpippim_n_sub->Draw("AP");
  
  TCanvas *cIMpippim_IMnpim_n_woSpdia_sub = new TCanvas("cIMpippim_IMnpim_n_woSpdia_sub","cIMpippim_IMnpim_n_woSpdia_sub",800,800);
  cIMpippim_IMnpim_n_woSpdia_sub->Divide(2,2,0,0);
  cIMpippim_IMnpim_n_woSpdia_sub->cd(3);
  IMpippim_IMnpim_n_woSpdia_sub->RebinY(2);
  IMpippim_IMnpim_n_woSpdia_sub->Draw("colz");
  cIMpippim_IMnpim_n_woSpdia_sub->cd(1);
  TH1D* IMnpim_n_woSpdia_sub = (TH1D*)IMpippim_IMnpim_n_woSpdia_sub->ProjectionX("IMnpim_n_woSpdia_sub");
  IMnpim_n_woSpdia_sub->Draw("HE");
  cIMpippim_IMnpim_n_woSpdia_sub->cd(4);
  TH1D* IMpippim_n_woSpdia_sub = (TH1D*)IMpippim_IMnpim_n_woSpdia_sub->ProjectionY("IMpippim_n_woSpdia_sub");
  TGraphErrors *gr_IMpippim_n_woSpdia_sub = new TGraphErrors();
  HistToRorateGraph(IMpippim_n_woSpdia_sub,*gr_IMpippim_n_woSpdia_sub);
  gr_IMpippim_n_woSpdia_sub->Draw("AP");
  
  if(!gROOT->IsBatch()){
    TCanvas *cIMpippim_IMnpip_n_sub_bin22 = new TCanvas("cIMpippim_IMnpip_n_sub_bin22","cIMpippim_IMnpip_n_sub_bin22",800,800);
    cIMpippim_IMnpip_n_sub_bin22->Divide(2,2,0,0);
    cIMpippim_IMnpip_n_sub_bin22->cd(3);
    IMpippim_IMnpip_n_bin_sub[22]->RebinX(4);
    IMpippim_IMnpip_n_bin_sub[22]->RebinY(4);
    IMpippim_IMnpip_n_bin_sub[22]->Draw("colz");
    cIMpippim_IMnpip_n_sub_bin22->cd(1);
    TH1D* IMnpip_n_sub_bin22 = (TH1D*)IMpippim_IMnpip_n_bin_sub[22]->ProjectionX("IMnpip_n_sub_bin22");
    IMnpip_n_sub_bin22->Draw("E");
    cIMpippim_IMnpip_n_sub_bin22->cd(4);
    TH1D* IMpippim_n_sub_bin22 = (TH1D*)IMpippim_IMnpip_n_bin_sub[22]->ProjectionY("IMpippim_n_sub_bin22");
    TGraphErrors *gr_IMpippim_n_sub_bin22 = new TGraphErrors();
    HistToRorateGraph(IMpippim_n_sub_bin22,*gr_IMpippim_n_sub_bin22);
    gr_IMpippim_n_sub_bin22->Draw("AP");
    
    TCanvas *cIMpippim_IMnpip_n_sub_qhi_bin22 = new TCanvas("cIMpippim_IMnpip_n_sub_qhi_bin22","cIMpippim_IMnpip_n_sub_qhi_bin22",800,800);
    cIMpippim_IMnpip_n_sub_qhi_bin22->Divide(2,2,0,0);
    cIMpippim_IMnpip_n_sub_qhi_bin22->cd(3);
    IMpippim_IMnpip_n_bin_qhi_sub[22]->RebinX(4);
    IMpippim_IMnpip_n_bin_qhi_sub[22]->RebinY(4);
    IMpippim_IMnpip_n_bin_qhi_sub[22]->Draw("colz");
    cIMpippim_IMnpip_n_sub_qhi_bin22->cd(1);
    TH1D* IMnpip_n_sub_qhi_bin22 = (TH1D*)IMpippim_IMnpip_n_bin_qhi_sub[22]->ProjectionX("IMnpip_n_sub_qhi_bin22");
    IMnpip_n_sub_qhi_bin22->Draw("E");
    cIMpippim_IMnpip_n_sub_qhi_bin22->cd(4);
    TH1D* IMpippim_n_sub_qhi_bin22 = (TH1D*)IMpippim_IMnpip_n_bin_qhi_sub[22]->ProjectionY("IMpippim_n_sub_qhi_bin22");
    TGraphErrors *gr_IMpippim_n_sub_qhi_bin22 = new TGraphErrors();
    HistToRorateGraph(IMpippim_n_sub_qhi_bin22,*gr_IMpippim_n_sub_qhi_bin22);
    gr_IMpippim_n_sub_qhi_bin22->Draw("AP");
    
    TCanvas *cIMpippim_IMnpip_n_woSmdia_sub_bin22 = new TCanvas("cIMpippim_IMnpip_n_woSmdia_sub_bin22","cIMpippim_IMnpip_n_woSmdia_sub_bin22",800,800);
    cIMpippim_IMnpip_n_woSmdia_sub_bin22->Divide(2,2,0,0);
    cIMpippim_IMnpip_n_woSmdia_sub_bin22->cd(3);
    IMpippim_IMnpip_n_woSmdia_bin_sub[22]->RebinX(4);
    IMpippim_IMnpip_n_woSmdia_bin_sub[22]->RebinY(4);
    IMpippim_IMnpip_n_woSmdia_bin_sub[22]->Draw("colz");
    cIMpippim_IMnpip_n_woSmdia_sub_bin22->cd(1);
    TH1D* IMnpip_n_woSmdia_sub_bin22 = (TH1D*)IMpippim_IMnpip_n_woSmdia_bin_sub[22]->ProjectionX("IMnpip_n_woSmdia_sub_bin22");
    IMnpip_n_woSmdia_sub_bin22->Draw("E");
    cIMpippim_IMnpip_n_woSmdia_sub_bin22->cd(4);
    TH1D* IMpippim_n_woSmdia_sub_bin22 = (TH1D*)IMpippim_IMnpip_n_woSmdia_bin_sub[22]->ProjectionY("IMpippim_n_woSmdia_sub_bin22");
    TGraphErrors *gr_IMpippim_n_woSmdia_sub_bin22 = new TGraphErrors();
    HistToRorateGraph(IMpippim_n_woSmdia_sub_bin22,*gr_IMpippim_n_woSmdia_sub_bin22);
    gr_IMpippim_n_woSmdia_sub_bin22->Draw("AP");
    
    TCanvas *cIMpippim_IMnpip_n_woSmdia_sub_qhi_bin22 = new TCanvas("cIMpippim_IMnpip_n_woSmdia_sub_qhi_bin22","cIMpippim_IMnpip_n_woSmdia_sub_qhi_bin22",800,800);
    cIMpippim_IMnpip_n_woSmdia_sub_qhi_bin22->Divide(2,2,0,0);
    cIMpippim_IMnpip_n_woSmdia_sub_qhi_bin22->cd(3);
    IMpippim_IMnpip_n_woSmdia_bin_qhi_sub[22]->RebinX(4);
    IMpippim_IMnpip_n_woSmdia_bin_qhi_sub[22]->RebinY(4);
    IMpippim_IMnpip_n_woSmdia_bin_qhi_sub[22]->Draw("colz");
    cIMpippim_IMnpip_n_woSmdia_sub_qhi_bin22->cd(1);
    TH1D* IMnpip_n_woSmdia_sub_qhi_bin22 = (TH1D*)IMpippim_IMnpip_n_woSmdia_bin_qhi_sub[22]->ProjectionX("IMnpip_n_woSmdia_sub_qhi_bin22");
    IMnpip_n_woSmdia_sub_qhi_bin22->Draw("E");
    cIMpippim_IMnpip_n_woSmdia_sub_qhi_bin22->cd(4);
    TH1D* IMpippim_n_woSmdia_sub_qhi_bin22 = (TH1D*)IMpippim_IMnpip_n_woSmdia_bin_qhi_sub[22]->ProjectionY("IMpippim_n_woSmdia_sub_qhi_bin22");
    TGraphErrors *gr_IMpippim_n_woSmdia_sub_qhi_bin22 = new TGraphErrors();
    HistToRorateGraph(IMpippim_n_woSmdia_sub_qhi_bin22,*gr_IMpippim_n_woSmdia_sub_qhi_bin22);
    gr_IMpippim_n_woSmdia_sub_qhi_bin22->Draw("AP");
    
    TCanvas *cIMpippim_IMnpip_n_woSmdia_sub_qlo_bin22 = new TCanvas("cIMpippim_IMnpip_n_woSmdia_sub_qlo_bin22","cIMpippim_IMnpip_n_woSmdia_sub_qlo_bin22",800,800);
    cIMpippim_IMnpip_n_woSmdia_sub_qlo_bin22->Divide(2,2,0,0);
    cIMpippim_IMnpip_n_woSmdia_sub_qlo_bin22->cd(3);
    IMpippim_IMnpip_n_woSmdia_bin_qlo_sub[22]->RebinX(4);
    IMpippim_IMnpip_n_woSmdia_bin_qlo_sub[22]->RebinY(4);
    IMpippim_IMnpip_n_woSmdia_bin_qlo_sub[22]->Draw("colz");
    cIMpippim_IMnpip_n_woSmdia_sub_qlo_bin22->cd(1);
    TH1D* IMnpip_n_woSmdia_sub_qlo_bin22 = (TH1D*)IMpippim_IMnpip_n_woSmdia_bin_qlo_sub[22]->ProjectionX("IMnpip_n_woSmdia_sub_qlo_bin22");
    IMnpip_n_woSmdia_sub_qlo_bin22->Draw("E");
    cIMpippim_IMnpip_n_woSmdia_sub_qlo_bin22->cd(4);
    TH1D* IMpippim_n_woSmdia_sub_qlo_bin22 = (TH1D*)IMpippim_IMnpip_n_woSmdia_bin_qlo_sub[22]->ProjectionY("IMpippim_n_woSmdia_sub_qlo_bin22");
    TGraphErrors *gr_IMpippim_n_woSmdia_sub_qlo_bin22 = new TGraphErrors();
    HistToRorateGraph(IMpippim_n_woSmdia_sub_qlo_bin22,*gr_IMpippim_n_woSmdia_sub_qlo_bin22);
    gr_IMpippim_n_woSmdia_sub_qlo_bin22->Draw("AP");
  
  
    TCanvas *cIMpippim_IMnpim_n_woSpdia_sub_bin22 = new TCanvas("cIMpippim_IMnpim_n_woSpdia_sub_bin22","cIMpippim_IMnpim_n_woSpdia_sub_bin22",800,800);
    cIMpippim_IMnpim_n_woSpdia_sub_bin22->Divide(2,2,0,0);
    cIMpippim_IMnpim_n_woSpdia_sub_bin22->cd(3);
    IMpippim_IMnpim_n_woSpdia_bin_sub[22]->RebinX(4);
    IMpippim_IMnpim_n_woSpdia_bin_sub[22]->RebinY(4);
    IMpippim_IMnpim_n_woSpdia_bin_sub[22]->Draw("colz");
    cIMpippim_IMnpim_n_woSpdia_sub_bin22->cd(1);
    TH1D* IMnpim_n_woSpdia_sub_bin22 = (TH1D*)IMpippim_IMnpim_n_woSpdia_bin_sub[22]->ProjectionX("IMnpim_n_woSpdia_sub_bin22");
    IMnpim_n_woSpdia_sub_bin22->Draw("E");
    cIMpippim_IMnpim_n_woSpdia_sub_bin22->cd(4);
    TH1D* IMpippim_n_woSpdia_sub_bin22 = (TH1D*)IMpippim_IMnpim_n_woSpdia_bin_sub[22]->ProjectionY("IMpippim_n_woSpdia_sub_bin22");
    TGraphErrors *gr_IMpippim_n_woSpdia_sub_bin22 = new TGraphErrors();
    HistToRorateGraph(IMpippim_n_woSpdia_sub_bin22,*gr_IMpippim_n_woSpdia_sub_bin22);
    gr_IMpippim_n_woSpdia_sub_bin22->Draw("AP");
    
    TCanvas *cIMpippim_IMnpim_n_woSpdia_sub_qhi_bin22 = new TCanvas("cIMpippim_IMnpim_n_woSpdia_sub_qhi_bin22","cIMpippim_IMnpim_n_woSpdia_sub_qhi_bin22",800,800);
    cIMpippim_IMnpim_n_woSpdia_sub_qhi_bin22->Divide(2,2,0,0);
    cIMpippim_IMnpim_n_woSpdia_sub_qhi_bin22->cd(3);
    IMpippim_IMnpim_n_woSpdia_bin_qhi_sub[22]->RebinX(4);
    IMpippim_IMnpim_n_woSpdia_bin_qhi_sub[22]->RebinY(4);
    IMpippim_IMnpim_n_woSpdia_bin_qhi_sub[22]->Draw("colz");
    cIMpippim_IMnpim_n_woSpdia_sub_qhi_bin22->cd(1);
    TH1D* IMnpim_n_woSpdia_sub_qhi_bin22 = (TH1D*)IMpippim_IMnpim_n_woSpdia_bin_qhi_sub[22]->ProjectionX("IMnpim_n_woSpdia_sub_qhi_bin22");
    IMnpim_n_woSpdia_sub_qhi_bin22->Draw("E");
    cIMpippim_IMnpim_n_woSpdia_sub_qhi_bin22->cd(4);
    TH1D* IMpippim_n_woSpdia_sub_qhi_bin22 = (TH1D*)IMpippim_IMnpim_n_woSpdia_bin_qhi_sub[22]->ProjectionY("IMpippim_n_woSpdia_sub_qhi_bin22");
    TGraphErrors *gr_IMpippim_n_woSpdia_sub_qhi_bin22 = new TGraphErrors();
    HistToRorateGraph(IMpippim_n_woSpdia_sub_qhi_bin22,*gr_IMpippim_n_woSpdia_sub_qhi_bin22);
    gr_IMpippim_n_woSpdia_sub_qhi_bin22->Draw("AP");
    
    TCanvas *cIMpippim_IMnpim_n_woSpdia_sub_qlo_bin22 = new TCanvas("cIMpippim_IMnpim_n_woSpdia_sub_qlo_bin22","cIMpippim_IMnpim_n_woSpdia_sub_qlo_bin22",800,800);
    cIMpippim_IMnpim_n_woSpdia_sub_qlo_bin22->Divide(2,2,0,0);
    cIMpippim_IMnpim_n_woSpdia_sub_qlo_bin22->cd(3);
    IMpippim_IMnpim_n_woSpdia_bin_qlo_sub[22]->RebinX(4);
    IMpippim_IMnpim_n_woSpdia_bin_qlo_sub[22]->RebinY(4);
    IMpippim_IMnpim_n_woSpdia_bin_qlo_sub[22]->Draw("colz");
    cIMpippim_IMnpim_n_woSpdia_sub_qlo_bin22->cd(1);
    TH1D* IMnpim_n_woSpdia_sub_qlo_bin22 = (TH1D*)IMpippim_IMnpim_n_woSpdia_bin_qlo_sub[22]->ProjectionX("IMnpim_n_woSpdia_sub_qlo_bin22");
    IMnpim_n_woSpdia_sub_qlo_bin22->Draw("E");
    cIMpippim_IMnpim_n_woSpdia_sub_qlo_bin22->cd(4);
    TH1D* IMpippim_n_woSpdia_sub_qlo_bin22 = (TH1D*)IMpippim_IMnpim_n_woSpdia_bin_qlo_sub[22]->ProjectionY("IMpippim_n_woSpdia_sub_qlo_bin22");
    TGraphErrors *gr_IMpippim_n_woSpdia_sub_qlo_bin22 = new TGraphErrors();
    HistToRorateGraph(IMpippim_n_woSpdia_sub_qlo_bin22,*gr_IMpippim_n_woSpdia_sub_qlo_bin22);
    gr_IMpippim_n_woSpdia_sub_qlo_bin22->Draw("AP");
  
  
  }




  TCanvas *c = NULL;
  TSeqCollection *SCol = gROOT->GetListOfCanvases();
  int size = SCol->GetSize();
  TIter next(SCol);
  TString pdfname = "K0SigmaTemp.pdf";
  for(int i=0;i<size;i++){
    //pdf->NewPage();
    c= (TCanvas*)next();
    c->Draw();
    c->cd();
    //inside the canvas
    //TPaveText *pt = new TPaveText(.74,.81,0.9,0.90,"NDC");
    c->Modified();
    c->Update();
    std::cout << c->GetName() << std::endl;
    //make 1 pdf file
    if(i==0) c->Print(pdfname+"(",Form("pdf Title:%s",c->GetTitle()));
    else if(i==size-1)c->Print(pdfname+")",Form("pdf Title:%s",c->GetTitle())); 
    else c->Print(pdfname,Form("pdf Title:%s",c->GetTitle())); 
    
    //make separated pdf files
    //c->Print(Form("pdf/%s.pdf",c->GetTitle()));
  }
  

  TIter nexthist2(gDirectory->GetList());
  TFile *fout = new TFile("K0SigmaTemp.root","RECREATE");
  fout->Print();
  fout->cd();
  TObject *obj = nullptr;
  while( (obj = (TObject*)nexthist2())!=nullptr) {
    obj->Write();
  }
  fout->cd();
  fout->Close();


  return;
  

}
