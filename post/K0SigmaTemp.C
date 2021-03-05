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
bool RebinMode = false;
bool Sidefar=false;

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
  TFile *fr[4]={NULL};
  TFile *fmix[4]={NULL};
  
  fr[0] = TFile::Open("evanaIMpisigma_npippim_v202_out_iso.root");
  fmix[0] = TFile::Open("evanaIMpisigma_npippim_v202_MIX_cut4_out_iso.root");
  fr[1] = TFile::Open("evanaIMpisigma_npippim_v202_out_iso_qlo.root");
  fmix[1] = TFile::Open("evanaIMpisigma_npippim_v202_MIX_cut4_out_iso_qlo.root");
  fr[2] = TFile::Open("evanaIMpisigma_npippim_v202_out_iso_qhi.root");
  fmix[2] = TFile::Open("evanaIMpisigma_npippim_v202_MIX_cut4_out_iso_qhi.root");
  fr[3] = TFile::Open("evanaIMpisigma_npippim_v202_out_iso_theta15.root");
  fmix[3] = TFile::Open("evanaIMpisigma_npippim_v202_MIX_cut4_out_iso_theta15.root");
  fr[0]->Print();
  fmix[0]->Print();
  fr[1]->Print();
  fmix[1]->Print();
  fr[2]->Print();
  fmix[2]->Print();
  fr[3]->Print();
  fmix[3]->Print();
  
  gROOT->SetBatch(1);
  gStyle->SetPalette(1);
  //gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  //TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0.);  

  const unsigned int nbintemplate = 100;
  const int nqcut=4;
  const int qstart=0;
  TH2F* IMpippim_IMnpip_n_data[nqcut];
  TH2F* IMpippim_IMnpim_n_data[nqcut];
  TH2F* IMpippim_IMnpip_n_bin_data[nbintemplate][nqcut];
  TH2F* IMpippim_IMnpim_n_bin_data[nbintemplate][nqcut];
  TH2F* IMpippim_IMnpip_n_mix[nqcut];
  TH2F* IMpippim_IMnpim_n_mix[nqcut];
  TH2F* IMpippim_IMnpip_n_bin_mix[nbintemplate][nqcut];
  TH2F* IMpippim_IMnpim_n_bin_mix[nbintemplate][nqcut];
  TH2F* IMpippim_IMnpip_n_sub[nqcut];
  TH2F* IMpippim_IMnpim_n_sub[nqcut];
  TH2F* IMpippim_IMnpip_n_bin_sub[nbintemplate][nqcut];
  TH2F* IMpippim_IMnpim_n_bin_sub[nbintemplate][nqcut];
  TH2F* IMpippim_IMnpip_wSid_n_Sp_data[nqcut];
  TH2F* IMpippim_IMnpip_wSid_n_Sp_bin_data[nbintemplate][nqcut];
  TH2F* IMpippim_IMnpip_wSid_n_Sp_mix[nqcut];
  TH2F* IMpippim_IMnpip_wSid_n_Sp_bin_mix[nbintemplate][nqcut];
  TH2F* IMpippim_IMnpip_wSid_n_Sp_sub[nqcut];
  TH2F* IMpippim_IMnpip_wSid_n_Sp_bin_sub[nbintemplate][nqcut];
  TH2F* IMpippim_IMnpim_wSid_n_Sm_data[nqcut];
  TH2F* IMpippim_IMnpim_wSid_n_Sm_bin_data[nbintemplate][nqcut];
  TH2F* IMpippim_IMnpim_wSid_n_Sm_mix[nqcut];
  TH2F* IMpippim_IMnpim_wSid_n_Sm_bin_mix[nbintemplate][nqcut];
  TH2F* IMpippim_IMnpim_wSid_n_Sm_sub[nqcut];
  TH2F* IMpippim_IMnpim_wSid_n_Sm_bin_sub[nbintemplate][nqcut];
  TH2F* IMpippim_IMnpip_wK0_n_data[nqcut];
  TH2F* IMpippim_IMnpip_wK0_n_bin_data[nbintemplate][nqcut];
  TH2F* IMpippim_IMnpip_wK0_n_mix[nqcut];
  TH2F* IMpippim_IMnpip_wK0_n_bin_mix[nbintemplate][nqcut];
  TH2F* IMpippim_IMnpip_wK0_n_sub[nqcut];
  TH2F* IMpippim_IMnpip_wK0_n_bin_sub[nbintemplate][nqcut];
  TH2F* IMpippim_IMnpim_wK0_n_data[nqcut];
  TH2F* IMpippim_IMnpim_wK0_n_bin_data[nbintemplate][nqcut];
  TH2F* IMpippim_IMnpim_wK0_n_mix[nqcut];
  TH2F* IMpippim_IMnpim_wK0_n_bin_mix[nbintemplate][nqcut];
  TH2F* IMpippim_IMnpim_wK0_n_sub[nqcut];
  TH2F* IMpippim_IMnpim_wK0_n_bin_sub[nbintemplate][nqcut];
  TH2F* IMpippim_IMnpip_wK0orwSid_n_data[nqcut];
  TH2F* IMpippim_IMnpim_wK0orwSid_n_data[nqcut];
  TH2F* IMpippim_IMnpip_wK0orwSid_n_bin_data[nbintemplate][nqcut];
  TH2F* IMpippim_IMnpim_wK0orwSid_n_bin_data[nbintemplate][nqcut];
  TH2F* IMpippim_IMnpip_wK0orwSid_n_mix[nqcut];
  TH2F* IMpippim_IMnpim_wK0orwSid_n_mix[nqcut];
  TH2F* IMpippim_IMnpip_wK0orwSid_n_bin_mix[nbintemplate][nqcut];
  TH2F* IMpippim_IMnpim_wK0orwSid_n_bin_mix[nbintemplate][nqcut];
  TH2F* IMpippim_IMnpip_wK0orwSid_n_sub[nqcut];
  TH2F* IMpippim_IMnpim_wK0orwSid_n_sub[nqcut];
  TH2F* IMpippim_IMnpip_wK0orwSid_n_bin_sub[nbintemplate][nqcut];
  TH2F* IMpippim_IMnpim_wK0orwSid_n_bin_sub[nbintemplate][nqcut];
  TH2F* IMnpim_IMnpip_n_data[nqcut];
  TH2F* IMnpim_IMnpip_n_bin_data[nbintemplate][nqcut];
  TH2F* IMnpim_IMnpip_n_mix[nqcut];
  TH2F* IMnpim_IMnpip_n_bin_mix[nbintemplate][nqcut];
  TH2F* IMnpim_IMnpip_n_sub[nqcut];
  TH2F* IMnpim_IMnpip_n_bin_sub[nbintemplate][nqcut];
  TH2F* IMnpim_IMnpip_wSid_n_data[nqcut];
  TH2F* IMnpim_IMnpip_wSid_n_bin_data[nbintemplate][nqcut];
  TH2F* IMnpim_IMnpip_wSid_n_mix[nqcut];
  TH2F* IMnpim_IMnpip_wSid_n_bin_mix[nbintemplate][nqcut];
  TH2F* IMnpim_IMnpip_wSid_n_sub[nqcut];
  TH2F* IMnpim_IMnpip_wSid_n_bin_sub[nbintemplate][nqcut];
  TH2F* IMnpim_IMnpip_wK0orwSid_n_data[nqcut];
  TH2F* IMnpim_IMnpip_wK0orwSid_n_bin_data[nbintemplate][nqcut];
  TH2F* IMnpim_IMnpip_wK0orwSid_n_mix[nqcut];
  TH2F* IMnpim_IMnpip_wK0orwSid_n_bin_mix[nbintemplate][nqcut];
  TH2F* IMnpim_IMnpip_wK0orwSid_n_sub[nqcut];
  TH2F* IMnpim_IMnpip_wK0orwSid_n_bin_sub[nbintemplate][nqcut];
  TH2F* IMnpim_IMnpip_wSid_n_Sp_data[nqcut];
  TH2F* IMnpim_IMnpip_wSid_n_Sp_bin_data[nbintemplate][nqcut];
  TH2F* IMnpim_IMnpip_wSid_n_Sp_mix[nqcut];
  TH2F* IMnpim_IMnpip_wSid_n_Sp_bin_mix[nbintemplate][nqcut];
  TH2F* IMnpim_IMnpip_wSid_n_Sp_sub[nqcut];
  TH2F* IMnpim_IMnpip_wSid_n_Sp_bin_sub[nbintemplate][nqcut];
  TH2F* IMnpim_IMnpip_wSid_n_Sm_data[nqcut];
  TH2F* IMnpim_IMnpip_wSid_n_Sm_bin_data[nbintemplate][nqcut];
  TH2F* IMnpim_IMnpip_wSid_n_Sm_mix[nqcut];
  TH2F* IMnpim_IMnpip_wSid_n_Sm_bin_mix[nbintemplate][nqcut];
  TH2F* IMnpim_IMnpip_wSid_n_Sm_sub[nqcut];
  TH2F* IMnpim_IMnpip_wSid_n_Sm_bin_sub[nbintemplate][nqcut];
  
  //for the overlap of S+ & S- & K0 counting 
  TH2F* q_IMnpipi_wK0_wSid_n_SpSm_data[nqcut];
  TH2F* q_IMnpipi_wK0_wSid_n_SpSm_mix[nqcut];
  TH2F* q_IMnpipi_wK0_wSid_n_SpSm_sub[nqcut];
  TH1D* IMnpipi_wK0_wSid_n_SpSm_sub[nqcut];
  TCanvas *cq_IMnpipi_wK0_wSid_n_SpSm_sub[nqcut];
  double OverlapCount[nbintemplate][nqcut]; 

  TCanvas *cIMpippim_IMnpip_n_all[nqcut];
  TH1D* IMnpip_wK0_n_sub[nqcut];
  TH1D* IMpippim_wSid_n_Sp_sub[nqcut];
  TCanvas *cIMpippim_IMnpim_n_all[nqcut];
  TH1D* IMnpim_wK0_n_sub[nqcut];
  TH1D* IMpippim_wSid_n_Sm_sub[nqcut];
  TCanvas *cIMnpim_IMnpip_n_all[nqcut];
  TH1D* IMnpip_wSid_n_Sm_sub[nqcut];
  TH1D* IMnpim_wSid_n_Sp_sub[nqcut];
  const char cqcut[][10]= {"all","qlo","qhi","theta"};
  for(int iqcut=qstart;iqcut<nqcut;iqcut++){
    IMpippim_IMnpip_n_data[iqcut] = (TH2F*)fr[iqcut]->Get("IMpippim_IMnpip_n");
    IMpippim_IMnpip_n_mix[iqcut] = (TH2F*)fmix[iqcut]->Get("IMpippim_IMnpip_n");
    IMpippim_IMnpip_n_sub[iqcut] = (TH2F*)IMpippim_IMnpip_n_data[iqcut]->Clone(Form("IMpippim_IMnpip_n_sub_%s",cqcut[iqcut]));
    IMpippim_IMnpip_n_sub[iqcut]->Add(IMpippim_IMnpip_n_mix[iqcut],-1.0); 
    IMpippim_IMnpip_n_sub[iqcut]->SetTitle(Form("IMpippim_IMnpip_n_sub_%s",cqcut[iqcut]));
    if(RemoveNegative)IMpippim_IMnpip_n_sub[iqcut]->SetMinimum(0);
    
    IMpippim_IMnpim_n_data[iqcut] = (TH2F*)fr[iqcut]->Get("IMpippim_IMnpim_n");
    IMpippim_IMnpim_n_mix[iqcut] = (TH2F*)fmix[iqcut]->Get("IMpippim_IMnpim_n");
    IMpippim_IMnpim_n_sub[iqcut] = (TH2F*)IMpippim_IMnpim_n_data[iqcut]->Clone(Form("IMpippim_IMnpim_n_sub_%s",cqcut[iqcut]));
    IMpippim_IMnpim_n_sub[iqcut]->Add(IMpippim_IMnpim_n_mix[iqcut],-1.0); 
    IMpippim_IMnpim_n_sub[iqcut]->SetTitle(Form("IMpippim_IMnpim_n_sub_%s",cqcut[iqcut]));
    if(RemoveNegative)IMpippim_IMnpim_n_sub[iqcut]->SetMinimum(0);
    
    IMpippim_IMnpip_wSid_n_Sp_data[iqcut] = (TH2F*)fr[iqcut]->Get("IMpippim_IMnpip_wSid_n_Sp");
    IMpippim_IMnpip_wSid_n_Sp_mix[iqcut] = (TH2F*)fmix[iqcut]->Get("IMpippim_IMnpip_wSid_n_Sp");
    IMpippim_IMnpip_wSid_n_Sp_sub[iqcut] = (TH2F*)IMpippim_IMnpip_wSid_n_Sp_data[iqcut]->Clone(Form("IMpippim_IMnpip_wSid_n_Sp_sub_%s",cqcut[iqcut]));
    IMpippim_IMnpip_wSid_n_Sp_sub[iqcut]->Add(IMpippim_IMnpip_wSid_n_Sp_mix[iqcut],-1.0);
    IMpippim_IMnpip_wSid_n_Sp_sub[iqcut]->SetTitle(Form("IMpippim_IMnpip_wSid_n_Sp_sub_%s",cqcut[iqcut]));
    if(RemoveNegative)IMpippim_IMnpip_wSid_n_Sp_sub[iqcut]->SetMinimum(0);

    IMpippim_IMnpim_wSid_n_Sm_data[iqcut] = (TH2F*)fr[iqcut]->Get("IMpippim_IMnpim_wSid_n_Sm");
    IMpippim_IMnpim_wSid_n_Sm_mix[iqcut] = (TH2F*)fmix[iqcut]->Get("IMpippim_IMnpim_wSid_n_Sm");
    IMpippim_IMnpim_wSid_n_Sm_sub[iqcut] = (TH2F*)IMpippim_IMnpim_wSid_n_Sm_data[iqcut]->Clone(Form("IMpippim_IMnpim_wSid_n_Sm_sub_%s",cqcut[iqcut]));
    IMpippim_IMnpim_wSid_n_Sm_sub[iqcut]->Add(IMpippim_IMnpim_wSid_n_Sm_mix[iqcut],-1.0);
    IMpippim_IMnpim_wSid_n_Sm_sub[iqcut]->SetTitle(Form("IMpippim_IMnpim_wSid_n_Sm_sub_%s",cqcut[iqcut]));
    if(RemoveNegative)IMpippim_IMnpim_wSid_n_Sm_sub[iqcut]->SetMinimum(0);
   
    IMpippim_IMnpip_wK0_n_data[iqcut] = (TH2F*)fr[iqcut]->Get("IMpippim_IMnpip_wK0_n");
    IMpippim_IMnpip_wK0_n_mix[iqcut] = (TH2F*)fmix[iqcut]->Get("IMpippim_IMnpip_wK0_n");
    IMpippim_IMnpip_wK0_n_sub[iqcut] = (TH2F*)IMpippim_IMnpip_wK0_n_data[iqcut]->Clone(Form("IMpippim_IMnpip_wK0_n_sub_%s",cqcut[iqcut]));
    IMpippim_IMnpip_wK0_n_sub[iqcut]->Add(IMpippim_IMnpip_wK0_n_mix[iqcut],-1.0);
    IMpippim_IMnpip_wK0_n_sub[iqcut]->SetTitle(Form("IMpippim_IMnpip_wK0_n_sub_%s",cqcut[iqcut]));
    if(RemoveNegative)IMpippim_IMnpip_wK0_n_sub[iqcut]->SetMinimum(0);

    IMpippim_IMnpim_wK0_n_data[iqcut] = (TH2F*)fr[iqcut]->Get("IMpippim_IMnpim_wK0_n");
    IMpippim_IMnpim_wK0_n_mix[iqcut] = (TH2F*)fmix[iqcut]->Get("IMpippim_IMnpim_wK0_n");
    IMpippim_IMnpim_wK0_n_sub[iqcut] = (TH2F*)IMpippim_IMnpim_wK0_n_data[iqcut]->Clone(Form("IMpippim_IMnpim_wK0_n_sub_%s",cqcut[iqcut]));
    IMpippim_IMnpim_wK0_n_sub[iqcut]->Add(IMpippim_IMnpim_wK0_n_mix[iqcut],-1.0);
    IMpippim_IMnpim_wK0_n_sub[iqcut]->SetTitle(Form("IMpippim_IMnpim_wK0_n_sub_%s",cqcut[iqcut]));
    if(RemoveNegative)IMpippim_IMnpim_wK0_n_sub[iqcut]->SetMinimum(0);

    IMpippim_IMnpip_wK0orwSid_n_data[iqcut] = (TH2F*)fr[iqcut]->Get("IMpippim_IMnpip_wK0orwSid_n");
    IMpippim_IMnpip_wK0orwSid_n_mix[iqcut] = (TH2F*)fmix[iqcut]->Get("IMpippim_IMnpip_wK0orwSid_n");
    IMpippim_IMnpip_wK0orwSid_n_sub[iqcut] = (TH2F*)IMpippim_IMnpip_wK0orwSid_n_data[iqcut]->Clone(Form("IMpippim_IMnpip_wK0orwSid_n_sub_%s",cqcut[iqcut]));
    IMpippim_IMnpip_wK0orwSid_n_sub[iqcut]->Add(IMpippim_IMnpip_wK0orwSid_n_mix[iqcut],-1.0); 
    IMpippim_IMnpip_wK0orwSid_n_sub[iqcut]->SetTitle(Form("IMpippim_IMnpip_wK0orwSid_n_sub_%s",cqcut[iqcut]));
    if(RemoveNegative)IMpippim_IMnpip_wK0orwSid_n_sub[iqcut]->SetMinimum(0);
    cIMpippim_IMnpip_n_all[iqcut] = new TCanvas(Form("cIMpippim_IMnpip_n_all_%s",cqcut[iqcut]),Form("cIMpippim_IMnpip_n_all_%s",cqcut[iqcut]));
    cIMpippim_IMnpip_n_all[iqcut]->Divide(2,2);
    cIMpippim_IMnpip_n_all[iqcut]->cd(3);
    IMpippim_IMnpip_wK0orwSid_n_sub[iqcut]->RebinX(4);
    IMpippim_IMnpip_wK0orwSid_n_sub[iqcut]->RebinY(5);
    IMpippim_IMnpip_wK0orwSid_n_sub[iqcut]->GetXaxis()->SetRangeUser(1,1.7);
    IMpippim_IMnpip_wK0orwSid_n_sub[iqcut]->GetYaxis()->SetRangeUser(0.2,1.0);
    IMpippim_IMnpip_wK0orwSid_n_sub[iqcut]->Draw("colz");
    cIMpippim_IMnpip_n_all[iqcut]->cd(1);
    IMnpip_wK0_n_sub[iqcut] = (TH1D*)IMpippim_IMnpip_wK0_n_sub[iqcut]->ProjectionX(Form("IMnpip_wK0_n_sub_%d",iqcut));
    IMnpip_wK0_n_sub[iqcut]->SetTitle(Form("IMnpip_wK0_n_sub_%d",iqcut));
    IMnpip_wK0_n_sub[iqcut]->GetXaxis()->SetRangeUser(1,1.7);
    IMnpip_wK0_n_sub[iqcut]->Draw("HIST");
    cIMpippim_IMnpip_n_all[iqcut]->cd(4);
    IMpippim_wSid_n_Sp_sub[iqcut]
    = (TH1D*)IMpippim_IMnpip_wSid_n_Sp_sub[iqcut]->ProjectionY(Form("IMpippim_wSid_n_Sp_sub_%d",iqcut));
    IMpippim_wSid_n_Sp_sub[iqcut]->SetTitle(Form("IMpippim_wSid_n_Sp_sub_%d",iqcut));
    IMpippim_wSid_n_Sp_sub[iqcut]->GetXaxis()->SetRangeUser(0.2,1.0);
    IMpippim_wSid_n_Sp_sub[iqcut]->Draw("HIST");

    IMpippim_IMnpim_wK0orwSid_n_data[iqcut] = (TH2F*)fr[iqcut]->Get("IMpippim_IMnpim_wK0orwSid_n");
    IMpippim_IMnpim_wK0orwSid_n_mix[iqcut] = (TH2F*)fmix[iqcut]->Get("IMpippim_IMnpim_wK0orwSid_n");
    IMpippim_IMnpim_wK0orwSid_n_sub[iqcut] = (TH2F*)IMpippim_IMnpim_wK0orwSid_n_data[iqcut]->Clone(Form("IMpippim_IMnpim_wK0orwSid_n_sub_%s",cqcut[iqcut]));
    IMpippim_IMnpim_wK0orwSid_n_sub[iqcut]->Add(IMpippim_IMnpim_wK0orwSid_n_mix[iqcut],-1.0); 
    IMpippim_IMnpim_wK0orwSid_n_sub[iqcut]->SetTitle(Form("IMpippim_IMnpim_wK0orwSid_n_sub_%s",cqcut[iqcut]));
    if(RemoveNegative)IMpippim_IMnpim_wK0orwSid_n_sub[iqcut]->SetMinimum(0);

    IMpippim_IMnpim_wK0orwSid_n_data[iqcut] = (TH2F*)fr[iqcut]->Get("IMpippim_IMnpim_wK0orwSid_n");
    IMpippim_IMnpim_wK0orwSid_n_mix[iqcut] = (TH2F*)fmix[iqcut]->Get("IMpippim_IMnpim_wK0orwSid_n");
    IMpippim_IMnpim_wK0orwSid_n_sub[iqcut] = (TH2F*)IMpippim_IMnpim_wK0orwSid_n_data[iqcut]->Clone(Form("IMpippim_IMnpim_wK0orwSid_n_sub_%s",cqcut[iqcut]));
    IMpippim_IMnpim_wK0orwSid_n_sub[iqcut]->Add(IMpippim_IMnpim_wK0orwSid_n_mix[iqcut],-1.0); 
    IMpippim_IMnpim_wK0orwSid_n_sub[iqcut]->SetTitle(Form("IMpippim_IMnpim_wK0orwSid_n_sub_%s",cqcut[iqcut]));
    if(RemoveNegative)IMpippim_IMnpim_wK0orwSid_n_sub[iqcut]->SetMinimum(0);
    cIMpippim_IMnpim_n_all[iqcut] = new TCanvas(Form("cIMpippim_IMnpim_n_all_%s",cqcut[iqcut]),Form("cIMpippim_IMnpim_n_all_%s",cqcut[iqcut]));
    cIMpippim_IMnpim_n_all[iqcut]->Divide(2,2);
    cIMpippim_IMnpim_n_all[iqcut]->cd(3);
    IMpippim_IMnpim_wK0orwSid_n_sub[iqcut]->RebinX(4);
    IMpippim_IMnpim_wK0orwSid_n_sub[iqcut]->RebinY(5);
    IMpippim_IMnpim_wK0orwSid_n_sub[iqcut]->GetXaxis()->SetRangeUser(1,1.7);
    IMpippim_IMnpim_wK0orwSid_n_sub[iqcut]->GetYaxis()->SetRangeUser(0.2,1.0);
    IMpippim_IMnpim_wK0orwSid_n_sub[iqcut]->Draw("colz");
    cIMpippim_IMnpim_n_all[iqcut]->cd(1);
    IMnpim_wK0_n_sub[iqcut] = (TH1D*)IMpippim_IMnpim_wK0_n_sub[iqcut]->ProjectionX(Form("IMnpim_wK0_n_sub_%d",iqcut));
    IMnpim_wK0_n_sub[iqcut]->SetTitle(Form("IMnpim_wK0_n_sub_%d",iqcut));
    IMnpim_wK0_n_sub[iqcut]->GetXaxis()->SetRangeUser(1,1.7);
    IMnpim_wK0_n_sub[iqcut]->Draw("HIST");
    cIMpippim_IMnpim_n_all[iqcut]->cd(4);
    IMpippim_wSid_n_Sm_sub[iqcut]
    = (TH1D*)IMpippim_IMnpim_wSid_n_Sm_sub[iqcut]->ProjectionY(Form("IMpippim_wSid_n_Sm_sub_%d",iqcut));
    IMpippim_wSid_n_Sm_sub[iqcut]->SetTitle(Form("IMpippim_wSid_n_Sm_sub_%d",iqcut));
    IMpippim_wSid_n_Sm_sub[iqcut]->GetXaxis()->SetRangeUser(0.2,1.0);
    IMpippim_wSid_n_Sm_sub[iqcut]->Draw("HIST");

    IMnpim_IMnpip_n_data[iqcut] = (TH2F*)fr[iqcut]->Get("IMnpim_IMnpip_dE_n");
    IMnpim_IMnpip_n_mix[iqcut] = (TH2F*)fmix[iqcut]->Get("IMnpim_IMnpip_dE_n");
    IMnpim_IMnpip_n_sub[iqcut] = (TH2F*)IMnpim_IMnpip_n_data[iqcut]->Clone(Form("IMnpim_IMnpip_n_%s",cqcut[iqcut]));
    IMnpim_IMnpip_n_sub[iqcut]->Add(IMnpim_IMnpip_n_mix[iqcut],-1);
    IMnpim_IMnpip_n_sub[iqcut]->SetTitle(Form("IMnpim_IMnpip_n_%s",cqcut[iqcut]));
    if(RemoveNegative)IMnpim_IMnpip_n_sub[iqcut]->SetMinimum(0);
    
    IMnpim_IMnpip_wSid_n_data[iqcut] = (TH2F*)fr[iqcut]->Get("IMnpim_IMnpip_dE_wSid_n");
    IMnpim_IMnpip_wSid_n_mix[iqcut] = (TH2F*)fmix[iqcut]->Get("IMnpim_IMnpip_dE_wSid_n");
    IMnpim_IMnpip_wSid_n_sub[iqcut] = (TH2F*)IMnpim_IMnpip_wSid_n_data[iqcut]->Clone(Form("IMnpim_IMnpip_wSid_n_%s",cqcut[iqcut]));
    IMnpim_IMnpip_wSid_n_sub[iqcut]->Add(IMnpim_IMnpip_wSid_n_mix[iqcut],-1);
    IMnpim_IMnpip_wSid_n_sub[iqcut]->SetTitle(Form("IMnpim_IMnpip_wSid_n_%s",cqcut[iqcut]));
    if(RemoveNegative)IMnpim_IMnpip_wSid_n_sub[iqcut]->SetMinimum(0);
    
    IMnpim_IMnpip_wK0orwSid_n_data[iqcut] = (TH2F*)fr[iqcut]->Get("IMnpim_IMnpip_dE_wK0orwSid_n");
    IMnpim_IMnpip_wK0orwSid_n_mix[iqcut] = (TH2F*)fmix[iqcut]->Get("IMnpim_IMnpip_dE_wK0orwSid_n");
    IMnpim_IMnpip_wK0orwSid_n_sub[iqcut] = (TH2F*)IMnpim_IMnpip_wK0orwSid_n_data[iqcut]->Clone(Form("IMnpim_IMnpip_wK0orwSid_n_%s",cqcut[iqcut]));
    IMnpim_IMnpip_wK0orwSid_n_sub[iqcut]->Add(IMnpim_IMnpip_wK0orwSid_n_mix[iqcut],-1);
    IMnpim_IMnpip_wK0orwSid_n_sub[iqcut]->SetTitle(Form("IMnpim_IMnpip_wK0orwSid_n_%s",cqcut[iqcut]));
    if(RemoveNegative)IMnpim_IMnpip_wK0orwSid_n_sub[iqcut]->SetMinimum(0);

    IMnpim_IMnpip_wSid_n_Sm_data[iqcut] = (TH2F*)fr[iqcut]->Get("IMnpim_IMnpip_dE_wSid_n_Sm");
    IMnpim_IMnpip_wSid_n_Sm_mix[iqcut] = (TH2F*)fmix[iqcut]->Get("IMnpim_IMnpip_dE_wSid_n_Sm");
    IMnpim_IMnpip_wSid_n_Sm_sub[iqcut] = (TH2F*)IMnpim_IMnpip_wSid_n_Sm_data[iqcut]->Clone(Form("IMnpim_IMnpip_wSid_n_Sm_%s",cqcut[iqcut]));
    IMnpim_IMnpip_wSid_n_Sm_sub[iqcut]->Add(IMnpim_IMnpip_wSid_n_Sm_mix[iqcut],-1);
    IMnpim_IMnpip_wSid_n_Sm_sub[iqcut]->SetTitle(Form("IMnpim_IMnpip_wSid_n_Sm_%s",cqcut[iqcut]));
    if(RemoveNegative)IMnpim_IMnpip_wSid_n_Sm_sub[iqcut]->SetMinimum(0);
    
    IMnpim_IMnpip_wSid_n_Sp_data[iqcut] = (TH2F*)fr[iqcut]->Get("IMnpim_IMnpip_dE_wSid_n_Sp");
    IMnpim_IMnpip_wSid_n_Sp_mix[iqcut] = (TH2F*)fmix[iqcut]->Get("IMnpim_IMnpip_dE_wSid_n_Sp");
    IMnpim_IMnpip_wSid_n_Sp_sub[iqcut] = (TH2F*)IMnpim_IMnpip_wSid_n_Sp_data[iqcut]->Clone(Form("IMnpim_IMnpip_wSid_n_Sp_%s",cqcut[iqcut]));
    IMnpim_IMnpip_wSid_n_Sp_sub[iqcut]->Add(IMnpim_IMnpip_wSid_n_Sp_mix[iqcut],-1);
    IMnpim_IMnpip_wSid_n_Sp_sub[iqcut]->SetTitle(Form("IMnpim_IMnpip_wSid_n_Sp_%s",cqcut[iqcut]));
    if(RemoveNegative)IMnpim_IMnpip_wSid_n_Sp_sub[iqcut]->SetMinimum(0);

    IMnpim_IMnpip_wSid_n_Sm_data[iqcut] = (TH2F*)fr[iqcut]->Get("IMnpim_IMnpip_dE_wSid_n_Sm");
    IMnpim_IMnpip_wSid_n_Sm_mix[iqcut] = (TH2F*)fmix[iqcut]->Get("IMnpim_IMnpip_dE_wSid_n_Sm");
    IMnpim_IMnpip_wSid_n_Sm_sub[iqcut] = (TH2F*)IMnpim_IMnpip_wSid_n_Sm_data[iqcut]->Clone(Form("IMnpim_IMnpip_wSid_n_Sm_%s",cqcut[iqcut]));
    IMnpim_IMnpip_wSid_n_Sm_sub[iqcut]->Add(IMnpim_IMnpip_wSid_n_Sm_mix[iqcut],-1);
    IMnpim_IMnpip_wSid_n_Sm_sub[iqcut]->SetTitle(Form("IMnpim_IMnpip_wSid_n_Sm_%s",cqcut[iqcut]));
    if(RemoveNegative)IMnpim_IMnpip_wSid_n_Sm_sub[iqcut]->SetMinimum(0);
    
    cIMnpim_IMnpip_n_all[iqcut] = new TCanvas(Form("cIMnpim_IMnpip_n_all_%s",cqcut[iqcut]),Form("cIMnpim_IMnpip_n_all_%s",cqcut[iqcut]));
    cIMnpim_IMnpip_n_all[iqcut]->Divide(2,2);
    cIMnpim_IMnpip_n_all[iqcut]->cd(3);
    IMnpim_IMnpip_wK0orwSid_n_sub[iqcut]->RebinX(4);
    IMnpim_IMnpip_wK0orwSid_n_sub[iqcut]->RebinY(4);
    IMnpim_IMnpip_wK0orwSid_n_sub[iqcut]->GetXaxis()->SetRangeUser(1,1.7);
    IMnpim_IMnpip_wK0orwSid_n_sub[iqcut]->GetYaxis()->SetRangeUser(1,1.7);
    IMnpim_IMnpip_wK0orwSid_n_sub[iqcut]->Draw("colz");
    cIMnpim_IMnpip_n_all[iqcut]->cd(1);
    IMnpip_wSid_n_Sm_sub[iqcut] = (TH1D*)IMnpim_IMnpip_wSid_n_Sm_sub[iqcut]->ProjectionX(Form("IMnpip_wSid_n_Sm_sub_%d",iqcut));
    IMnpip_wSid_n_Sm_sub[iqcut]->SetTitle(Form("IMnpip_wSid_n_Sm_sub_%d",iqcut));
    IMnpip_wSid_n_Sm_sub[iqcut]->GetXaxis()->SetRangeUser(1,1.7);
    IMnpip_wSid_n_Sm_sub[iqcut]->Draw("HIST");
    cIMnpim_IMnpip_n_all[iqcut]->cd(4);
    IMnpim_wSid_n_Sp_sub[iqcut]
    = (TH1D*)IMnpim_IMnpip_wSid_n_Sp_sub[iqcut]->ProjectionY(Form("IMnpim_wSid_n_Sp_sub_%d",iqcut));
    IMnpim_wSid_n_Sp_sub[iqcut]->SetTitle(Form("IMnpim_wSid_n_Sp_sub_%d",iqcut));
    IMnpim_wSid_n_Sp_sub[iqcut]->GetXaxis()->SetRangeUser(1,1.7);
    IMnpim_wSid_n_Sp_sub[iqcut]->Draw("HIST");
    

    q_IMnpipi_wK0_wSid_n_SpSm_data[iqcut] = (TH2F*)fr[iqcut]->Get("q_IMnpipi_wK0_wSid_n_SpSm");
    q_IMnpipi_wK0_wSid_n_SpSm_mix[iqcut] = (TH2F*)fmix[iqcut]->Get("q_IMnpipi_wK0_wSid_n_SpSm");
    q_IMnpipi_wK0_wSid_n_SpSm_sub[iqcut] = (TH2F*)q_IMnpipi_wK0_wSid_n_SpSm_data[iqcut]->Clone(Form("q_IMnpipi_wK0_wSid_n_SpSm_sub_%s",cqcut[iqcut]));
    q_IMnpipi_wK0_wSid_n_SpSm_sub[iqcut]->Add(q_IMnpipi_wK0_wSid_n_SpSm_mix[iqcut],-1.0);
    q_IMnpipi_wK0_wSid_n_SpSm_sub[iqcut]->SetTitle(Form("q_IMnpipi_wK0_wSid_n_SpSm_%s",cqcut[iqcut]));
    cq_IMnpipi_wK0_wSid_n_SpSm_sub[iqcut] = new TCanvas(Form("q_IMnpipi_wK0_wSid_n_SpSm_%s",cqcut[iqcut]),Form("q_IMnpipi_wK0_wSid_n_SpSm_%s",cqcut[iqcut]));
    IMnpipi_wK0_wSid_n_SpSm_sub[iqcut] = (TH1D*)q_IMnpipi_wK0_wSid_n_SpSm_sub[iqcut]->ProjectionX(Form("IMnpipi_wK0_wSid_n_SpSm_%d",iqcut));
    IMnpipi_wK0_wSid_n_SpSm_sub[iqcut]->Draw("HISTE");
    
    for(unsigned int ibin=0;ibin<nbintemplate;ibin++){
      OverlapCount[ibin][iqcut] = IMnpipi_wK0_wSid_n_SpSm_sub[iqcut]->GetBinContent(ibin+1);
      if(OverlapCount[ibin][iqcut]<0.0) OverlapCount[ibin][iqcut]=0.0;
    }

    for(unsigned int ibin=0;ibin<nbintemplate;ibin++){
      std::cout << iqcut << "  " << ibin << std::endl;
      IMpippim_IMnpip_n_bin_data[ibin][iqcut] = (TH2F*)fr[iqcut]->Get(Form("IMpippim_IMnpip_n_bin%d",ibin));
      IMpippim_IMnpip_n_bin_mix[ibin][iqcut] = (TH2F*)fmix[iqcut]->Get(Form("IMpippim_IMnpip_n_bin%d",ibin));
      IMpippim_IMnpip_n_bin_sub[ibin][iqcut] = (TH2F*)IMpippim_IMnpip_n_bin_data[ibin][iqcut]->Clone(Form("IMpippim_IMnpip_n_sub_bin%d_%s",ibin,cqcut[iqcut]));
      IMpippim_IMnpip_n_bin_sub[ibin][iqcut]->Add(IMpippim_IMnpip_n_bin_mix[ibin][iqcut],-1.0);
      if(RemoveNegative)IMpippim_IMnpip_n_bin_sub[ibin][iqcut]->SetMinimum(0);

      IMpippim_IMnpim_n_bin_data[ibin][iqcut] = (TH2F*)fr[iqcut]->Get(Form("IMpippim_IMnpim_n_bin%d",ibin));
      IMpippim_IMnpim_n_bin_mix[ibin][iqcut] = (TH2F*)fmix[iqcut]->Get(Form("IMpippim_IMnpim_n_bin%d",ibin));
      IMpippim_IMnpim_n_bin_sub[ibin][iqcut] = (TH2F*)IMpippim_IMnpim_n_bin_data[ibin][iqcut]->Clone(Form("IMpippim_IMnpim_n_sub_bin%d_%s",ibin,cqcut[iqcut]));
      IMpippim_IMnpim_n_bin_sub[ibin][iqcut]->Add(IMpippim_IMnpim_n_bin_mix[ibin][iqcut],-1.0);
      if(RemoveNegative)IMpippim_IMnpim_n_bin_sub[ibin][iqcut]->SetMinimum(0);

      IMpippim_IMnpip_wSid_n_Sp_bin_data[ibin][iqcut] = (TH2F*)fr[iqcut]->Get(Form("IMpippim_IMnpip_wSid_n_Sp_bin%d",ibin));
      IMpippim_IMnpip_wSid_n_Sp_bin_mix[ibin][iqcut] = (TH2F*)fmix[iqcut]->Get(Form("IMpippim_IMnpip_wSid_n_Sp_bin%d",ibin));
      IMpippim_IMnpip_wSid_n_Sp_bin_sub[ibin][iqcut]
      = (TH2F*)IMpippim_IMnpip_wSid_n_Sp_bin_data[ibin][iqcut]->Clone(Form("IMpippim_IMnpip_wSid_n_Sp_bin%d_sub%s",ibin,cqcut[iqcut]));
      IMpippim_IMnpip_wSid_n_Sp_bin_sub[ibin][iqcut]->Add(IMpippim_IMnpip_wSid_n_Sp_bin_mix[ibin][iqcut],-1.0);
      IMpippim_IMnpip_wSid_n_Sp_bin_sub[ibin][iqcut]->SetTitle(Form("IMpippim_IMnpip_wSid_n_Sp_bin%d_sub%s",ibin,cqcut[iqcut]));
      IMpippim_IMnpip_wSid_n_Sp_bin_sub[ibin][iqcut]->GetXaxis()->SetRangeUser(1,1.7);
      IMpippim_IMnpip_wSid_n_Sp_bin_sub[ibin][iqcut]->GetYaxis()->SetRangeUser(0.2,1.0);
      if(RemoveNegative)IMpippim_IMnpip_wSid_n_Sp_bin_sub[ibin][iqcut]->SetMinimum(0);

      IMpippim_IMnpim_wSid_n_Sm_bin_data[ibin][iqcut] = (TH2F*)fr[iqcut]->Get(Form("IMpippim_IMnpim_wSid_n_Sm_bin%d",ibin));
      IMpippim_IMnpim_wSid_n_Sm_bin_mix[ibin][iqcut] = (TH2F*)fmix[iqcut]->Get(Form("IMpippim_IMnpim_wSid_n_Sm_bin%d",ibin));
      IMpippim_IMnpim_wSid_n_Sm_bin_sub[ibin][iqcut]
      = (TH2F*)IMpippim_IMnpim_wSid_n_Sm_bin_data[ibin][iqcut]->Clone(Form("IMpippim_IMnpim_wSid_n_Sm_bin%d_sub%s",ibin,cqcut[iqcut]));
      IMpippim_IMnpim_wSid_n_Sm_bin_sub[ibin][iqcut]->Add(IMpippim_IMnpim_wSid_n_Sm_bin_mix[ibin][iqcut],-1.0);
      IMpippim_IMnpim_wSid_n_Sm_bin_sub[ibin][iqcut]->SetTitle(Form("IMpippim_IMnpim_wSid_n_Sm_bin%d_sub%s",ibin,cqcut[iqcut]));
      IMpippim_IMnpim_wSid_n_Sm_bin_sub[ibin][iqcut]->GetXaxis()->SetRangeUser(1,1.7);
      IMpippim_IMnpim_wSid_n_Sm_bin_sub[ibin][iqcut]->GetYaxis()->SetRangeUser(0.2,1.0);
      if(RemoveNegative)IMpippim_IMnpim_wSid_n_Sm_bin_sub[ibin][iqcut]->SetMinimum(0);

      IMpippim_IMnpip_wK0_n_bin_data[ibin][iqcut] = (TH2F*)fr[iqcut]->Get(Form("IMpippim_IMnpip_wK0_n_bin%d",ibin));
      IMpippim_IMnpip_wK0_n_bin_mix[ibin][iqcut] = (TH2F*)fmix[iqcut]->Get(Form("IMpippim_IMnpip_wK0_n_bin%d",ibin));
      IMpippim_IMnpip_wK0_n_bin_sub[ibin][iqcut] = (TH2F*)IMpippim_IMnpip_wK0_n_bin_data[ibin][iqcut]->Clone(Form("IMpippim_IMnpip_wK0_n_bin%d_sub%s",ibin,cqcut[iqcut]));
      IMpippim_IMnpip_wK0_n_bin_sub[ibin][iqcut]->Add(IMpippim_IMnpip_wK0_n_bin_mix[ibin][iqcut],-1.0);
      IMpippim_IMnpip_wK0_n_bin_sub[ibin][iqcut]->SetTitle(Form("IMpippim_IMnpip_wK0_n_bin%d_sub%s",ibin,cqcut[iqcut]));
      IMpippim_IMnpip_wK0_n_bin_sub[ibin][iqcut]->GetXaxis()->SetRangeUser(1,1.7);
      IMpippim_IMnpip_wK0_n_bin_sub[ibin][iqcut]->GetYaxis()->SetRangeUser(0.2,1.0);
      if(RemoveNegative)IMpippim_IMnpip_wK0_n_bin_sub[ibin][iqcut]->SetMinimum(0);

      IMpippim_IMnpim_wK0_n_bin_data[ibin][iqcut] = (TH2F*)fr[iqcut]->Get(Form("IMpippim_IMnpim_wK0_n_bin%d",ibin));
      IMpippim_IMnpim_wK0_n_bin_mix[ibin][iqcut] = (TH2F*)fmix[iqcut]->Get(Form("IMpippim_IMnpim_wK0_n_bin%d",ibin));
      IMpippim_IMnpim_wK0_n_bin_sub[ibin][iqcut] = (TH2F*)IMpippim_IMnpim_wK0_n_bin_data[ibin][iqcut]->Clone(Form("IMpippim_IMnpim_wK0_n_bin%d_sub%s",ibin,cqcut[iqcut]));
      IMpippim_IMnpim_wK0_n_bin_sub[ibin][iqcut]->Add(IMpippim_IMnpim_wK0_n_bin_mix[ibin][iqcut],-1.0);
      IMpippim_IMnpim_wK0_n_bin_sub[ibin][iqcut]->SetTitle(Form("IMpippim_IMnpim_wK0_n_bin%d_sub%s",ibin,cqcut[iqcut]));
      IMpippim_IMnpim_wK0_n_bin_sub[ibin][iqcut]->GetXaxis()->SetRangeUser(1,1.7);
      IMpippim_IMnpim_wK0_n_bin_sub[ibin][iqcut]->GetYaxis()->SetRangeUser(0.2,1.0);
      if(RemoveNegative)IMpippim_IMnpim_wK0_n_bin_sub[ibin][iqcut]->SetMinimum(0);

      IMpippim_IMnpip_wK0orwSid_n_bin_data[ibin][iqcut] = (TH2F*)fr[iqcut]->Get(Form("IMpippim_IMnpip_wK0orwSid_n_bin%d",ibin));
      IMpippim_IMnpip_wK0orwSid_n_bin_mix[ibin][iqcut] = (TH2F*)fmix[iqcut]->Get(Form("IMpippim_IMnpip_wK0orwSid_n_bin%d",ibin));
      IMpippim_IMnpip_wK0orwSid_n_bin_sub[ibin][iqcut]
      = (TH2F*)IMpippim_IMnpip_wK0orwSid_n_bin_data[ibin][iqcut]->Clone(Form("IMpippim_IMnpip_wK0orwSid_n_sub_bin%d_%s",ibin,cqcut[iqcut]));
      IMpippim_IMnpip_wK0orwSid_n_bin_sub[ibin][iqcut]->Add(IMpippim_IMnpip_wK0orwSid_n_bin_mix[ibin][iqcut],-1.0);
      if(RemoveNegative)IMpippim_IMnpip_wK0orwSid_n_bin_sub[ibin][iqcut]->SetMinimum(0);
      
      IMpippim_IMnpim_wK0orwSid_n_bin_data[ibin][iqcut] = (TH2F*)fr[iqcut]->Get(Form("IMpippim_IMnpim_wK0orwSid_n_bin%d",ibin));
      IMpippim_IMnpim_wK0orwSid_n_bin_mix[ibin][iqcut] = (TH2F*)fmix[iqcut]->Get(Form("IMpippim_IMnpim_wK0orwSid_n_bin%d",ibin));
      IMpippim_IMnpim_wK0orwSid_n_bin_sub[ibin][iqcut]
      = (TH2F*)IMpippim_IMnpim_wK0orwSid_n_bin_data[ibin][iqcut]->Clone(Form("IMpippim_IMnpim_wK0orwSid_n_sub_bin%d_%s",ibin,cqcut[iqcut]));
      IMpippim_IMnpim_wK0orwSid_n_bin_sub[ibin][iqcut]->Add(IMpippim_IMnpim_wK0orwSid_n_bin_mix[ibin][iqcut],-1.0);
      if(RemoveNegative)IMpippim_IMnpim_wK0orwSid_n_bin_sub[ibin][iqcut]->SetMinimum(0);

      IMnpim_IMnpip_wSid_n_bin_data[ibin][iqcut] = (TH2F*)fr[iqcut]->Get(Form("IMnpim_IMnpip_wSid_n_bin%d",ibin));
      IMnpim_IMnpip_wSid_n_bin_mix[ibin][iqcut] = (TH2F*)fmix[iqcut]->Get(Form("IMnpim_IMnpip_wSid_n_bin%d",ibin));
      IMnpim_IMnpip_wSid_n_bin_sub[ibin][iqcut] 
      = (TH2F*)IMnpim_IMnpip_wSid_n_bin_data[ibin][iqcut]->Clone(Form("IMnpim_IMnpip_wSid_n_bin%d_%s",ibin,cqcut[iqcut]));
      IMnpim_IMnpip_wSid_n_bin_sub[ibin][iqcut]->Add(IMnpim_IMnpip_wSid_n_bin_mix[ibin][iqcut],-1);
      float binlow=1.0+(float)ibin*1./nbintemplate;
      float binhigh=1.0+(float)(ibin+1.0)*1./nbintemplate;
      IMnpim_IMnpip_wSid_n_bin_sub[ibin][iqcut]->SetTitle(Form("IMnpim_IMnpip_wSid_n_bin%d_%s %0.2f-%0.2f",ibin,cqcut[iqcut],binlow,binhigh));
      if(RemoveNegative)IMnpim_IMnpip_wSid_n_bin_sub[ibin][iqcut]->SetMinimum(0);
      
      IMnpim_IMnpip_wK0orwSid_n_bin_data[ibin][iqcut] = (TH2F*)fr[iqcut]->Get(Form("IMnpim_IMnpip_dE_wK0orwSid_n_bin%d",ibin));
      IMnpim_IMnpip_wK0orwSid_n_bin_mix[ibin][iqcut] = (TH2F*)fmix[iqcut]->Get(Form("IMnpim_IMnpip_dE_wK0orwSid_n_bin%d",ibin));
      IMnpim_IMnpip_wK0orwSid_n_bin_sub[ibin][iqcut] 
      = (TH2F*)IMnpim_IMnpip_wK0orwSid_n_bin_data[ibin][iqcut]->Clone(Form("IMnpim_IMnpip_wK0orwSid_n_bin%d_%s",ibin,cqcut[iqcut]));
      IMnpim_IMnpip_wK0orwSid_n_bin_sub[ibin][iqcut]->Add(IMnpim_IMnpip_wK0orwSid_n_bin_mix[ibin][iqcut],-1);
      IMnpim_IMnpip_wK0orwSid_n_bin_sub[ibin][iqcut]->SetTitle(Form("IMnpim_IMnpip_wK0orwSid_n_bin%d_%s %0.2f-%0.2f",ibin,cqcut[iqcut],binlow,binhigh));
      if(RemoveNegative)IMnpim_IMnpip_wK0orwSid_n_bin_sub[ibin][iqcut]->SetMinimum(0);
      
      IMnpim_IMnpip_wSid_n_Sp_bin_data[ibin][iqcut] = (TH2F*)fr[iqcut]->Get(Form("IMnpim_IMnpip_wSid_n_Sp_bin%d",ibin));
      IMnpim_IMnpip_wSid_n_Sp_bin_mix[ibin][iqcut] = (TH2F*)fmix[iqcut]->Get(Form("IMnpim_IMnpip_wSid_n_Sp_bin%d",ibin));
      IMnpim_IMnpip_wSid_n_Sp_bin_sub[ibin][iqcut] 
      = (TH2F*)IMnpim_IMnpip_wSid_n_Sp_bin_data[ibin][iqcut]->Clone(Form("IMnpim_IMnpip_wSid_n_Sp_bin%d_%s",ibin,cqcut[iqcut]));
      IMnpim_IMnpip_wSid_n_Sp_bin_sub[ibin][iqcut]->Add(IMnpim_IMnpip_wSid_n_Sp_bin_mix[ibin][iqcut],-1);
      IMnpim_IMnpip_wSid_n_Sp_bin_sub[ibin][iqcut]->SetTitle(Form("IMnpim_IMnpip_wSid_n_Sp_bin%d_%s",ibin,cqcut[iqcut]));
      IMnpim_IMnpip_wSid_n_Sp_bin_sub[ibin][iqcut]->GetXaxis()->SetRangeUser(1,1.7);
      IMnpim_IMnpip_wSid_n_Sp_bin_sub[ibin][iqcut]->GetYaxis()->SetRangeUser(1,1.7);
      if(RemoveNegative)IMnpim_IMnpip_wSid_n_Sp_bin_sub[ibin][iqcut]->SetMinimum(0);

      IMnpim_IMnpip_wSid_n_Sm_bin_data[ibin][iqcut] = (TH2F*)fr[iqcut]->Get(Form("IMnpim_IMnpip_wSid_n_Sm_bin%d",ibin));
      IMnpim_IMnpip_wSid_n_Sm_bin_mix[ibin][iqcut] = (TH2F*)fmix[iqcut]->Get(Form("IMnpim_IMnpip_wSid_n_Sm_bin%d",ibin));
      IMnpim_IMnpip_wSid_n_Sm_bin_sub[ibin][iqcut]
      = (TH2F*)IMnpim_IMnpip_wSid_n_Sm_bin_data[ibin][iqcut]->Clone(Form("IMnpim_IMnpip_wSid_n_Sm_bin%d_%s",ibin,cqcut[iqcut]));
      IMnpim_IMnpip_wSid_n_Sm_bin_sub[ibin][iqcut]->Add(IMnpim_IMnpip_wSid_n_Sm_bin_mix[ibin][iqcut],-1);
      IMnpim_IMnpip_wSid_n_Sm_bin_sub[ibin][iqcut]->SetTitle(Form("IMnpim_IMnpip_wSid_n_Sm_bin%d_%s",ibin,cqcut[iqcut]));
      IMnpim_IMnpip_wSid_n_Sm_bin_sub[ibin][iqcut]->GetXaxis()->SetRangeUser(1,1.7);
      IMnpim_IMnpip_wSid_n_Sm_bin_sub[ibin][iqcut]->GetYaxis()->SetRangeUser(1,1.7);
      if(RemoveNegative)IMnpim_IMnpip_wSid_n_Sm_bin_sub[ibin][iqcut]->SetMinimum(0);
    }
  }
  
  
  int ngroup = 3;
  TH1D* IMnpipi_overlapdeco_K0[ngroup][nqcut];
  TH1D* IMnpipi_overlapdeco_Sp[ngroup][nqcut];
  TH1D* IMnpipi_overlapdeco_Sm[ngroup][nqcut];

  TCanvas *cIMpippim_IMnpip_n_sub_bin[nbintemplate][nqcut];
  TH1D* IMnpip_n_sub_bin[nbintemplate][nqcut];
  TH1D* IMpippim_n_sub_bin_1[nbintemplate][nqcut];
  TGraphErrors *gr_IMpippim_n_sub_bin_1[nbintemplate][nqcut];
  TCanvas *cIMpippim_IMnpip_n_sub_bin_cut[nbintemplate][nqcut];
  TH1D* IMnpip_wK0_n_sub_bin[nbintemplate][nqcut];
  TH1D* IMnpip_wK0_n_sub_bin_select[nbintemplate][nqcut];
  TH1D* IMnpip_wK0_n_sub_bin_sidelo[nbintemplate][nqcut];
  TH1D* IMnpip_wK0_n_sub_bin_sidehi[nbintemplate][nqcut];
  TH1D* IMpippim_wSid_n_Sp_sub_bin[nbintemplate][nqcut];
  TH1D* IMpippim_wSid_n_Sp_sub_bin_select[nbintemplate][nqcut];
  TH1D* IMpippim_wSid_n_Sp_sub_bin_sidelo[nbintemplate][nqcut];
  TH1D* IMpippim_wSid_n_Sp_sub_bin_sidelo2[nbintemplate][nqcut];//2 sigma to 6 sigma away when there is no higher side
  TH1D* IMpippim_wSid_n_Sp_sub_bin_sidehi[nbintemplate][nqcut];
  //IM(pi+pi-) vs IM(npi-) correlations
  TCanvas *cIMpippim_IMnpim_n_sub_bin[nbintemplate][nqcut];
  TH1D* IMnpim_n_sub_bin[nbintemplate][nqcut];
  TH1D* IMpippim_n_sub_bin_2[nbintemplate][nqcut];
  TGraphErrors *gr_IMpippim_n_sub_bin_2[nbintemplate][nqcut];
  TCanvas *cIMpippim_IMnpim_n_sub_bin_cut[nbintemplate][nqcut];
  TH1D* IMnpim_wK0_n_sub_bin[nbintemplate][nqcut];
  TH1D* IMnpim_wK0_n_sub_bin_select[nbintemplate][nqcut];
  TH1D* IMnpim_wK0_n_sub_bin_sidelo[nbintemplate][nqcut];
  TH1D* IMnpim_wK0_n_sub_bin_sidehi[nbintemplate][nqcut];
  TH1D* IMpippim_wSid_n_Sm_sub_bin[nbintemplate][nqcut];
  TH1D* IMpippim_wSid_n_Sm_sub_bin_select[nbintemplate][nqcut];
  TH1D* IMpippim_wSid_n_Sm_sub_bin_sidelo[nbintemplate][nqcut];
  TH1D* IMpippim_wSid_n_Sm_sub_bin_sidelo2[nbintemplate][nqcut];
  TH1D* IMpippim_wSid_n_Sm_sub_bin_sidehi[nbintemplate][nqcut];
  
  TCanvas *cIMnpim_IMnpip_n_sub_bin[nbintemplate][nqcut];
  TH1D* IMnpip_n_sub_bin_2[nbintemplate][nqcut];
  TH1D* IMnpim_n_sub_bin_2[nbintemplate][nqcut];
  TGraphErrors *gr_IMnpim_n_sub_bin_2[nbintemplate][nqcut];
  TCanvas *cIMnpim_IMnpip_n_sub_bin_cut[nbintemplate][nqcut];
  TH1D* IMnpip_wSid_n_Sm_sub_bin[nbintemplate][nqcut];
  TH1D* IMnpip_wSid_n_Sm_sub_bin_select[nbintemplate][nqcut];
  TH1D* IMnpip_wSid_n_Sm_sub_bin_sidelo[nbintemplate][nqcut];//2 sigma
  TH1D* IMnpip_wSid_n_Sm_sub_bin_sidelo2[nbintemplate][nqcut];//4 sigma
  TH1D* IMnpip_wSid_n_Sm_sub_bin_sidehi[nbintemplate][nqcut];//2 sigma
  TH1D* IMnpip_wSid_n_Sm_sub_bin_sidehi2[nbintemplate][nqcut];//4 simga
  TH1D* IMnpim_wSid_n_Sp_sub_bin[nbintemplate][nqcut];
  TH1D* IMnpim_wSid_n_Sp_sub_bin_select[nbintemplate][nqcut];
  TH1D* IMnpim_wSid_n_Sp_sub_bin_sidelo[nbintemplate][nqcut];//2 sigma
  TH1D* IMnpim_wSid_n_Sp_sub_bin_sidelo2[nbintemplate][nqcut];//4 sigma
  TH1D* IMnpim_wSid_n_Sp_sub_bin_sidehi[nbintemplate][nqcut];//2 sigma
  TH1D* IMnpim_wSid_n_Sp_sub_bin_sidehi2[nbintemplate][nqcut];//4 sigma
  
  for(int ig=0;ig<ngroup;ig++){
    for(int iq=qstart;iq<nqcut;iq++){
      if(ig!=2)IMnpipi_overlapdeco_K0[ig][iq] = new TH1D(Form("IMnpipi_overlapdeco_K0_g%d_%d",ig,iq),Form("IMnpipi_overlapdeco_K0_g%d_%d",ig,iq),100,1,2);
      if(ig!=1)IMnpipi_overlapdeco_Sp[ig][iq] = new TH1D(Form("IMnpipi_overlapdeco_Sp_g%d_%d",ig,iq),Form("IMnpipi_overlapdeco_Sp_g%d_%d",ig,iq),100,1,2);
      if(ig!=0)IMnpipi_overlapdeco_Sm[ig][iq] = new TH1D(Form("IMnpipi_overlapdeco_Sm_g%d_%d",ig,iq),Form("IMnpipi_overlapdeco_Sm_g%d_%d",ig,iq),100,1,2);
    }
  }

  for(unsigned int ibin=39;ibin<nbintemplate;ibin++){
    for(int iq=qstart;iq<nqcut;iq++){
      /*
      //first canvas to display IM(pi+pi-) vs IM(npi+) and their simple projection just to see the signal and other background 
      //no complicated cut at this moment
      cIMpippim_IMnpip_n_sub_bin[ibin][iq] = new TCanvas(Form("cIMpippim_IMnpip_n_sub_bin%d_%d",ibin,iq),Form("cIMpippim_IMnpip_n_sub_bin%d_%d",ibin,iq),800,800);
      cIMpippim_IMnpip_n_sub_bin[ibin][iq]->Divide(2,2,0,0);
      cIMpippim_IMnpip_n_sub_bin[ibin][iq]->cd(3);
      IMpippim_IMnpip_n_bin_sub[ibin][iq]->RebinX(4);
      IMpippim_IMnpip_n_bin_sub[ibin][iq]->RebinY(5);
      IMpippim_IMnpip_n_bin_sub[ibin][iq]->GetXaxis()->SetRangeUser(1,1.7);
      IMpippim_IMnpip_n_bin_sub[ibin][iq]->GetYaxis()->SetRangeUser(0.2,1.0);
      IMpippim_IMnpip_n_bin_sub[ibin][iq]->Draw("colz");
      cIMpippim_IMnpip_n_sub_bin[ibin][iq]->cd(1);
      IMnpip_n_sub_bin[ibin][iq] = (TH1D*)IMpippim_IMnpip_n_bin_sub[ibin][iq]->ProjectionX(Form("IMnpip_n_sub_bin%d_%d",ibin,iq));
      IMnpip_n_sub_bin[ibin][iq]->SetTitle(Form("IMnpip_n_sub_bin%d_%d",ibin,iq));
      IMnpip_n_sub_bin[ibin][iq]->Draw("H");
      cIMpippim_IMnpip_n_sub_bin[ibin][iq]->cd(4);
      IMpippim_n_sub_bin_1[ibin][iq] = (TH1D*)IMpippim_IMnpip_n_bin_sub[ibin][iq]->ProjectionY(Form("IMpippim_n_sub_bin_1_%d_%d",ibin,iq));
      gr_IMpippim_n_sub_bin_1[ibin][iq] = new TGraphErrors();
      HistToRorateGraph(IMpippim_n_sub_bin_1[ibin][iq],*gr_IMpippim_n_sub_bin_1[ibin][iq]);
      gr_IMpippim_n_sub_bin_1[ibin][iq]->Draw("AP");
      */
      //2nd canvas to display the calculation of decomposion of K0, Sigma+,Sigma- 
      //2D plot in canvas(3) has cut to select K0 or Sigma+
      //projection plot in canvas(1) select Sigma+ events candidate + (Sigma+ & K0 ) overlap events
      //projection plot in canvas(4) select K0     events candidate + (Sigma+ & K0 ) overlap events 
      cIMpippim_IMnpip_n_sub_bin_cut[ibin][iq] 
      = new TCanvas(Form("cIMpippim_IMnpip_n_sub_bin_cut_%d_%d",ibin,iq),Form("cIMpippim_IMnpip_n_sub_bin_cut%d_%d",ibin,iq),800,800);
      cIMpippim_IMnpip_n_sub_bin_cut[ibin][iq]->Divide(2,2);
      cIMpippim_IMnpip_n_sub_bin_cut[ibin][iq]->cd(3);
      
      IMpippim_IMnpip_wK0orwSid_n_bin_sub[ibin][iq]->RebinX(4);
      IMpippim_IMnpip_wK0orwSid_n_bin_sub[ibin][iq]->RebinY(5);
      IMpippim_IMnpip_wK0orwSid_n_bin_sub[ibin][iq]->GetXaxis()->SetRangeUser(1,1.7);
      IMpippim_IMnpip_wK0orwSid_n_bin_sub[ibin][iq]->GetYaxis()->SetRangeUser(0.2,1.0);
      IMpippim_IMnpip_wK0orwSid_n_bin_sub[ibin][iq]->Draw("colz");
      cIMpippim_IMnpip_n_sub_bin_cut[ibin][iq]->cd(1);
      IMnpip_wK0_n_sub_bin[ibin][iq] = (TH1D*)IMpippim_IMnpip_wK0_n_bin_sub[ibin][iq]->ProjectionX(Form("IMnpip_wK0_n_sub_bin%d_%d",ibin,iq));
      IMnpip_wK0_n_sub_bin[ibin][iq]->SetTitle(Form("IMnpip_wK0_n_sub_bin%d_%d",ibin,iq));
      if(RebinMode)IMnpip_wK0_n_sub_bin[ibin][iq]->RebinX(5);
      IMnpip_wK0_n_sub_bin[ibin][iq]->Draw("HIST");
      IMnpip_wK0_n_sub_bin_select[ibin][iq] = (TH1D*)IMnpip_wK0_n_sub_bin[ibin][iq]->Clone(Form("IMnpip_wK0_n_sub_bin%d_%d_select",ibin,iq));
      IMnpip_wK0_n_sub_bin_select[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmap_MIN,anacuts::Sigmap_MAX);
      IMnpip_wK0_n_sub_bin_select[ibin][iq]->SetLineColor(2);
      //IMnpip_wK0_n_sub_bin_select[ibin][iq]->SetFillColor(2);
      IMnpip_wK0_n_sub_bin_select[ibin][iq]->Draw("HISTsame");
      IMnpip_wK0_n_sub_bin_sidelo[ibin][iq] = (TH1D*)IMnpip_wK0_n_sub_bin[ibin][iq]->Clone(Form("IMnpip_wK0_n_sub_bin%d_%d_sidelo",ibin,iq));
      IMnpip_wK0_n_sub_bin_sidehi[ibin][iq] = (TH1D*)IMnpip_wK0_n_sub_bin[ibin][iq]->Clone(Form("IMnpip_wK0_n_sub_bin%d_%d_sidehi",ibin,iq));
      IMnpip_wK0_n_sub_bin_sidelo[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmap_MIN-4.0*anacuts::Sigmap_sigma, anacuts::Sigmap_MIN-2.0*anacuts::Sigmap_sigma);
      IMnpip_wK0_n_sub_bin_sidehi[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmap_MAX+2.0*anacuts::Sigmap_sigma, anacuts::Sigmap_MAX+4.0*anacuts::Sigmap_sigma);
      IMnpip_wK0_n_sub_bin_sidelo[ibin][iq]->SetLineColor(4);
      //IMnpip_wK0_n_sub_bin_sidelo[ibin][iq]->SetFillColor(4);
      IMnpip_wK0_n_sub_bin_sidehi[ibin][iq]->SetLineColor(4);
      //IMnpip_wK0_n_sub_bin_sidehi[ibin][iq]->SetFillColor(4);
      //double err_IMnpip_wK0=0.0;
      if(ibin==44 || ibin==45){
        IMnpip_wK0_n_sub_bin_sidehi[ibin][iq]->Draw("HISTsame");
      }else{
        IMnpip_wK0_n_sub_bin_sidelo[ibin][iq]->Draw("HISTsame");
        IMnpip_wK0_n_sub_bin_sidehi[ibin][iq]->Draw("HISTsame");
      }
      cIMpippim_IMnpip_n_sub_bin_cut[ibin][iq]->cd(4);
      IMpippim_wSid_n_Sp_sub_bin[ibin][iq]
      = (TH1D*)IMpippim_IMnpip_wSid_n_Sp_bin_sub[ibin][iq]->ProjectionY(Form("IMpippim_wSid_n_Sp_sub_bin%d_%d",ibin,iq));
      IMpippim_wSid_n_Sp_sub_bin[ibin][iq]->SetTitle(Form("IMpippim_wSid_n_Sp_sub_bin%d_%d",ibin,iq));
      if(RebinMode)IMpippim_wSid_n_Sp_sub_bin[ibin][iq]->RebinX(5);
      IMpippim_wSid_n_Sp_sub_bin[ibin][iq]->Draw("HIST");
      IMpippim_wSid_n_Sp_sub_bin_select[ibin][iq] = (TH1D*)IMpippim_wSid_n_Sp_sub_bin[ibin][iq]->Clone(Form("IMpippim_wSid_n_Sp_sub_bin_%d_%d_select",ibin,iq));
      IMpippim_wSid_n_Sp_sub_bin_select[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::pipi_MIN_narrow,anacuts::pipi_MAX_narrow);
      IMpippim_wSid_n_Sp_sub_bin_select[ibin][iq]->SetLineColor(2);
      //IMpippim_wSid_n_Sp_sub_bin_select[ibin][iq]->SetFillColor(2);
      IMpippim_wSid_n_Sp_sub_bin_select[ibin][iq]->Draw("HISTsame");
      IMpippim_wSid_n_Sp_sub_bin_sidelo[ibin][iq] = (TH1D*)IMpippim_wSid_n_Sp_sub_bin[ibin][iq]->Clone(Form("IMpippim_wSid_n_Sp_sub_bin_%d_%d_sidelo",ibin,iq));
      IMpippim_wSid_n_Sp_sub_bin_sidelo2[ibin][iq] = (TH1D*)IMpippim_wSid_n_Sp_sub_bin[ibin][iq]->Clone(Form("IMpippim_wSid_n_Sp_sub_bin_%d_%d_sidelo2",ibin,iq));
      IMpippim_wSid_n_Sp_sub_bin_sidehi[ibin][iq] = (TH1D*)IMpippim_wSid_n_Sp_sub_bin[ibin][iq]->Clone(Form("IMpippim_wSid_n_Sp_sub_bin_%d_%d_sidehi",ibin,iq));
      IMpippim_wSid_n_Sp_sub_bin_sidelo[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::pipi_MIN_narrow-4.0*anacuts::K0_sigma,anacuts::pipi_MIN_narrow-2.0*anacuts::K0_sigma);
      //IMpippim_wSid_n_Sp_sub_bin_sidelo2[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::pipi_MIN_narrow-6.0*anacuts::K0_sigma,anacuts::pipi_MIN_narrow-2.0*anacuts::K0_sigma);
      if(ibin==45){
        IMpippim_wSid_n_Sp_sub_bin_sidelo2[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::pipi_MIN_narrow-14.0*anacuts::K0_sigma,anacuts::pipi_MIN_narrow-10.0*anacuts::K0_sigma);
      }else{
        IMpippim_wSid_n_Sp_sub_bin_sidelo2[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::pipi_MIN_narrow-6.0*anacuts::K0_sigma,anacuts::pipi_MIN_narrow-2.0*anacuts::K0_sigma);
      }
      IMpippim_wSid_n_Sp_sub_bin_sidehi[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::pipi_MAX_narrow+2.0*anacuts::K0_sigma,anacuts::pipi_MAX_narrow+4.0*anacuts::K0_sigma);
      IMpippim_wSid_n_Sp_sub_bin_sidelo[ibin][iq]->SetLineColor(4);
      //IMpippim_wSid_n_Sp_sub_bin_sidelo[ibin][iq]->SetFillColor(4);
      IMpippim_wSid_n_Sp_sub_bin_sidelo2[ibin][iq]->SetLineColor(4);
      //IMpippim_wSid_n_Sp_sub_bin_sidelo2[ibin][iq]->SetFillColor(4);
      IMpippim_wSid_n_Sp_sub_bin_sidehi[ibin][iq]->SetLineColor(4);
      //IMpippim_wSid_n_Sp_sub_bin_sidehi[ibin][iq]->SetFillColor(4);

      if(ibin>51){
        IMpippim_wSid_n_Sp_sub_bin_sidelo[ibin][iq]->Draw("HISTsame");
        IMpippim_wSid_n_Sp_sub_bin_sidehi[ibin][iq]->Draw("HISTsame");
      }else{
        IMpippim_wSid_n_Sp_sub_bin_sidelo2[ibin][iq]->Draw("HISTsame");
      }
      
      //gr_IMpippim_wSid_n_Sp_sub_bin[ibin][iq] = new TGraphErrors();
      //HistToRorateGraph(IMpippim_wSid_n_Sp_sub_bin[ibin][iq],*gr_IMpippim_wSid_n_Sp_sub_bin[ibin][iq]);
      //gr_IMpippim_wSid_n_Sp_sub_bin[ibin][iq]->Draw("AP");
      //TBox *box = new TBox(0,anacuts::pipi_MIN_narrow,2000,anacuts::pipi_MAX_narrow);
      //box->SetFillColor(4);
      //box->SetFillStyle(3002);
      //box->Draw();

      cIMpippim_IMnpip_n_sub_bin_cut[ibin][iq]->cd(2);
      TPaveText *pt = new TPaveText(.05,.05,.95,.7);
      double inteSp = IMpippim_wSid_n_Sp_sub_bin_select[ibin][iq]->Integral();  //->Integral(binpipi_MIN,binpipi_MAX);
      double inteSpsidelo = IMpippim_wSid_n_Sp_sub_bin_sidelo[ibin][iq]->Integral();  //>Integral(binpipi_MIN_4sigma,binpipi_MIN-1);
      double inteSpsidelo2 = IMpippim_wSid_n_Sp_sub_bin_sidelo2[ibin][iq]->Integral();  //>Integral(binpipi_MIN_4sigma,binpipi_MIN-1);
      double inteSpsidehi = IMpippim_wSid_n_Sp_sub_bin_sidehi[ibin][iq]->Integral();  //>Integral(binpipi_MIN_4sigma,binpipi_MIN-1);
      double Spestimated=0.0;
      double Spoverlap = 0.0; 
      if(ibin>51) Spestimated = inteSpsidelo+inteSpsidehi;
      else        Spestimated = inteSpsidelo2;
      Spoverlap = inteSp-Spestimated;
      if(Spestimated<0.0) Spestimated=0.0;
      if(Spoverlap<0.0) {
        Spoverlap = 0.0;
        Spestimated=0.0;
      }
      double inteK0 = IMnpip_wK0_n_sub_bin_select[ibin][iq]->Integral();  //->Integral(binnpip_MIN,binnpip_MAX);
      double inteK0sidelo = IMnpip_wK0_n_sub_bin_sidelo[ibin][iq]->Integral();  //(binnpip_MIN_2sigma,binnpip_MIN-1);
      double inteK0sidehi = IMnpip_wK0_n_sub_bin_sidehi[ibin][iq]->Integral();//(binnpip_MAX+1,binnpip_MAX_2sigma);
      double K0estimated = 0.0;
      if(ibin==44 || ibin==45) K0estimated = inteK0sidehi*2.0;
      else K0estimated = inteK0sidelo + inteK0sidehi;
      double K0overlap = inteK0 - K0estimated;
      if(K0estimated<0.0) K0estimated = 0.0;
      if(K0overlap<0.0){
        K0overlap = 0.0;
        K0estimated = 0.0;
      }
      pt->AddText(Form("IM(n#pi^{-}#pi^{+})  %0.2f-%0.2f %s",1.0+ibin*1.0/nbintemplate,1.0+(ibin+1.0)/nbintemplate,cqcut[iq])); 
      pt->AddText(Form("Sigma+ count     %0.2f ",inteSp));
      if(ibin>51){
        pt->AddText(Form("Sigma+ side low  %0.2f ",inteSpsidelo));
        pt->AddText(Form("Sigma+ side high %0.2f ",inteSpsidehi));
      }else{
        pt->AddText(Form("Sigma+ side low (4 sigma)  %0.2f ",inteSpsidelo2));
      }
      pt->AddText(Form("Sigma+ estimated %0.2f ", Spestimated));
      pt->AddText(Form("K0 count %0.2f ",inteK0));
      if(ibin==44 || ibin==45){
        pt->AddText(Form("K0 side high*2  %0.2f ",inteK0sidehi*2));
      }else{
        pt->AddText(Form("K0 side low   %0.2f ",inteK0sidelo));
        pt->AddText(Form("K0 side high  %0.2f ",inteK0sidehi));
      }
      pt->AddText(Form("K0 estimated %0.2f ", K0estimated));
      pt->Draw();
      IMnpipi_overlapdeco_K0[0][iq]->SetBinContent(ibin,K0estimated);
      IMnpipi_overlapdeco_K0[0][iq]->SetBinError(ibin,sqrt(K0estimated));
      IMnpipi_overlapdeco_Sp[0][iq]->SetBinContent(ibin,Spestimated);
      IMnpipi_overlapdeco_Sp[0][iq]->SetBinError(ibin,sqrt(Spestimated));
      
      /*
      //first canvas to display IM(pi+pi-) vs IM(npi-) and their simple projection just to see the signal and other background 
      //no complicated cut at this moment
      cIMpippim_IMnpim_n_sub_bin[ibin][iq] = new TCanvas(Form("cIMpippim_IMnpim_n_sub_bin%d_%d",ibin,iq),Form("cIMpippim_IMnpim_n_sub_bin%d_%d",ibin,iq),800,800);
      cIMpippim_IMnpim_n_sub_bin[ibin][iq]->Divide(2,2,0,0);
      cIMpippim_IMnpim_n_sub_bin[ibin][iq]->cd(3);
      IMpippim_IMnpim_n_bin_sub[ibin][iq]->RebinX(4);
      IMpippim_IMnpim_n_bin_sub[ibin][iq]->RebinY(5);
      IMpippim_IMnpim_n_bin_sub[ibin][iq]->GetXaxis()->SetRangeUser(1,1.7);
      IMpippim_IMnpim_n_bin_sub[ibin][iq]->GetYaxis()->SetRangeUser(0.2,1.0);
      IMpippim_IMnpim_n_bin_sub[ibin][iq]->Draw("colz");
      cIMpippim_IMnpim_n_sub_bin[ibin][iq]->cd(1);
      IMnpim_n_sub_bin[ibin][iq] = (TH1D*)IMpippim_IMnpim_n_bin_sub[ibin][iq]->ProjectionX(Form("IMnpim_n_sub_bin%d_%d",ibin,iq));
      IMnpim_n_sub_bin[ibin][iq]->SetTitle(Form("IMnpim_n_sub_bin%d_%d",ibin,iq));
      IMnpim_n_sub_bin[ibin][iq]->Draw("HE");
      cIMpippim_IMnpim_n_sub_bin[ibin][iq]->cd(4);
      IMpippim_n_sub_bin_2[ibin][iq] = (TH1D*)IMpippim_IMnpim_n_bin_sub[ibin][iq]->ProjectionY(Form("IMpippim_n_sub_bin_2_%d_%d",ibin,iq));
      gr_IMpippim_n_sub_bin_2[ibin][iq] = new TGraphErrors();
      HistToRorateGraph(IMpippim_n_sub_bin_2[ibin][iq],*gr_IMpippim_n_sub_bin_2[ibin][iq]);
      gr_IMpippim_n_sub_bin_2[ibin][iq]->Draw("AP");
      */

      //2nd canvas to display the calculation of decomposion of K0, Sigma+,Sigma- 
      //2D plot in canvas(3) has cut to select K0 or Sigma-
      //projection plot in canvas(1) select Sigma- events candidate + (Sigma- & K0 ) overlap events
      //projection plot in canvas(4) select K0     events candidate + (Sigma- & K0 ) overlap events 
      cIMpippim_IMnpim_n_sub_bin_cut[ibin][iq]
      = new TCanvas(Form("cIMpippim_IMnpim_n_sub_bin_cut_%d_%d",ibin,iq),Form("cIMpippim_IMnpim_n_sub_bin_cut%d_%d",ibin,iq),800,800);
      cIMpippim_IMnpim_n_sub_bin_cut[ibin][iq]->Divide(2,2);
      cIMpippim_IMnpim_n_sub_bin_cut[ibin][iq]->cd(3);

      IMpippim_IMnpim_wK0orwSid_n_bin_sub[ibin][iq]->RebinX(4);
      IMpippim_IMnpim_wK0orwSid_n_bin_sub[ibin][iq]->RebinY(5);
      IMpippim_IMnpim_wK0orwSid_n_bin_sub[ibin][iq]->GetXaxis()->SetRangeUser(1,1.7);
      IMpippim_IMnpim_wK0orwSid_n_bin_sub[ibin][iq]->GetYaxis()->SetRangeUser(0.2,1.0);
      IMpippim_IMnpim_wK0orwSid_n_bin_sub[ibin][iq]->Draw("colz");
      cIMpippim_IMnpim_n_sub_bin_cut[ibin][iq]->cd(1);
      IMnpim_wK0_n_sub_bin[ibin][iq] = (TH1D*)IMpippim_IMnpim_wK0_n_bin_sub[ibin][iq]->ProjectionX(Form("IMnpim_wK0_n_sub_bin%d_%d",ibin,iq));
      IMnpim_wK0_n_sub_bin[ibin][iq]->SetTitle(Form("IMnpim_wK0_n_sub_bin%d_%d",ibin,iq));
      if(RebinMode) IMnpim_wK0_n_sub_bin[ibin][iq]->RebinX(5);
      IMnpim_wK0_n_sub_bin[ibin][iq]->Draw("HIST");
      IMnpim_wK0_n_sub_bin_select[ibin][iq] = (TH1D*)IMnpim_wK0_n_sub_bin[ibin][iq]->Clone(Form("IMnpim_wK0_n_sub_bin%d_%d_select",ibin,iq));
      IMnpim_wK0_n_sub_bin_select[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmam_MIN,anacuts::Sigmam_MAX);
      IMnpim_wK0_n_sub_bin_select[ibin][iq]->SetLineColor(2);
      //IMnpim_wK0_n_sub_bin_select[ibin][iq]->SetFillColor(2);
      IMnpim_wK0_n_sub_bin_select[ibin][iq]->Draw("HISTsame");
      IMnpim_wK0_n_sub_bin_sidelo[ibin][iq] = (TH1D*)IMnpim_wK0_n_sub_bin[ibin][iq]->Clone(Form("IMnpim_wK0_n_sub_bin%d_%d_selectlo",ibin,iq));
      IMnpim_wK0_n_sub_bin_sidehi[ibin][iq] = (TH1D*)IMnpim_wK0_n_sub_bin[ibin][iq]->Clone(Form("IMnpim_wK0_n_sub_bin%d_%d_selecthi",ibin,iq));
      IMnpim_wK0_n_sub_bin_sidelo[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmam_MIN-4.0*anacuts::Sigmam_sigma,anacuts::Sigmam_MIN-2.0*anacuts::Sigmam_sigma);
      IMnpim_wK0_n_sub_bin_sidehi[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmam_MAX+2.0*anacuts::Sigmam_sigma,anacuts::Sigmam_MAX+4.0*anacuts::Sigmam_sigma);
      IMnpim_wK0_n_sub_bin_sidelo[ibin][iq]->SetLineColor(4);
      //IMnpim_wK0_n_sub_bin_sidelo[ibin][iq]->SetFillColor(4);
      IMnpim_wK0_n_sub_bin_sidehi[ibin][iq]->SetLineColor(4);
      //IMnpim_wK0_n_sub_bin_sidehi[ibin][iq]->SetFillColor(4);
      if(ibin==44 || ibin==45){
        IMnpim_wK0_n_sub_bin_sidehi[ibin][iq]->Draw("HISTsame");
      }else{
        IMnpim_wK0_n_sub_bin_sidelo[ibin][iq]->Draw("HISTsame");
        IMnpim_wK0_n_sub_bin_sidehi[ibin][iq]->Draw("HISTsame");
      }
      
      cIMpippim_IMnpim_n_sub_bin_cut[ibin][iq]->cd(4);
      IMpippim_wSid_n_Sm_sub_bin[ibin][iq]
      = (TH1D*)IMpippim_IMnpim_wSid_n_Sm_bin_sub[ibin][iq]->ProjectionY(Form("IMpippim_wSid_n_Sm_sub_bin%d_%d",ibin,iq));
      IMpippim_wSid_n_Sm_sub_bin[ibin][iq]->SetTitle(Form("IMpippim_wSid_n_Sm_sub_bin%d_%d",ibin,iq));
      if(RebinMode)IMpippim_wSid_n_Sm_sub_bin[ibin][iq]->RebinX(5);
      IMpippim_wSid_n_Sm_sub_bin[ibin][iq]->Draw("HIST");
      IMpippim_wSid_n_Sm_sub_bin_select[ibin][iq] = (TH1D*)IMpippim_wSid_n_Sm_sub_bin[ibin][iq]->Clone(Form("IMpippim_wSid_n_Sm_sub_bin_%d_%d_select",ibin,iq));
      IMpippim_wSid_n_Sm_sub_bin_select[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::pipi_MIN_narrow,anacuts::pipi_MAX_narrow);
      IMpippim_wSid_n_Sm_sub_bin_select[ibin][iq]->SetLineColor(2);
      //IMpippim_wSid_n_Sm_sub_bin_select[ibin][iq]->SetFillColor(2);
      IMpippim_wSid_n_Sm_sub_bin_select[ibin][iq]->Draw("HISTsame");
      IMpippim_wSid_n_Sm_sub_bin_sidelo[ibin][iq] = (TH1D*)IMpippim_wSid_n_Sm_sub_bin[ibin][iq]->Clone(Form("IMpippim_wSid_n_Sm_sub_bin_%d_%d_sidelo",ibin,iq));
      IMpippim_wSid_n_Sm_sub_bin_sidelo2[ibin][iq] = (TH1D*)IMpippim_wSid_n_Sm_sub_bin[ibin][iq]->Clone(Form("IMpippim_wSid_n_Sm_sub_bin_%d_%d_sidelo2",ibin,iq));
      IMpippim_wSid_n_Sm_sub_bin_sidehi[ibin][iq] = (TH1D*)IMpippim_wSid_n_Sm_sub_bin[ibin][iq]->Clone(Form("IMpippim_wSid_n_Sm_sub_bin_%d_%d_sidehi",ibin,iq));
      IMpippim_wSid_n_Sm_sub_bin_sidelo[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::pipi_MIN_narrow-4.0*anacuts::K0_sigma,anacuts::pipi_MIN_narrow-2.0*anacuts::K0_sigma);
      if(ibin==45){
        IMpippim_wSid_n_Sm_sub_bin_sidelo2[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::pipi_MIN_narrow-14.0*anacuts::K0_sigma,anacuts::pipi_MIN_narrow-10.0*anacuts::K0_sigma);
      }else{
        IMpippim_wSid_n_Sm_sub_bin_sidelo2[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::pipi_MIN_narrow-6.0*anacuts::K0_sigma,anacuts::pipi_MIN_narrow-2.0*anacuts::K0_sigma);
      }
        
      IMpippim_wSid_n_Sm_sub_bin_sidehi[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::pipi_MAX_narrow+2.0*anacuts::K0_sigma,anacuts::pipi_MAX_narrow+4.0*anacuts::K0_sigma);
      IMpippim_wSid_n_Sm_sub_bin_sidelo[ibin][iq]->SetLineColor(4);
      //IMpippim_wSid_n_Sm_sub_bin_sidelo[ibin][iq]->SetFillColor(4);
      IMpippim_wSid_n_Sm_sub_bin_sidelo2[ibin][iq]->SetLineColor(4);
      //IMpippim_wSid_n_Sm_sub_bin_sidelo2[ibin][iq]->SetFillColor(4);
      IMpippim_wSid_n_Sm_sub_bin_sidehi[ibin][iq]->SetLineColor(4);
      //IMpippim_wSid_n_Sm_sub_bin_sidehi[ibin][iq]->SetFillColor(4);
      if(ibin>51){
        IMpippim_wSid_n_Sm_sub_bin_sidelo[ibin][iq]->Draw("HISTsame");
        IMpippim_wSid_n_Sm_sub_bin_sidehi[ibin][iq]->Draw("HISTsame");
      }else{
        IMpippim_wSid_n_Sm_sub_bin_sidelo2[ibin][iq]->Draw("HISTsame");
      }

      //gr_IMpippim_wSid_n_Sm_sub_bin[ibin][iq] = new TGraphErrors();
      //HistToRorateGraph(IMpippim_wSid_n_Sm_sub_bin[ibin][iq],*gr_IMpippim_wSid_n_Sm_sub_bin[ibin][iq]);
      //gr_IMpippim_wSid_n_Sm_sub_bin[ibin][iq]->Draw("AP");
      
      cIMpippim_IMnpim_n_sub_bin_cut[ibin][iq]->cd(2);
      TPaveText *pt2 = new TPaveText(.05,.05,.95,.7);
      double inteSm_g2 = IMpippim_wSid_n_Sm_sub_bin_select[ibin][iq]->Integral();
      double inteSmsidelo_g2 = IMpippim_wSid_n_Sm_sub_bin_sidelo[ibin][iq]->Integral();
      double inteSmsidelo2_g2 = IMpippim_wSid_n_Sm_sub_bin_sidelo2[ibin][iq]->Integral();
      double inteSmsidehi_g2 = IMpippim_wSid_n_Sm_sub_bin_sidehi[ibin][iq]->Integral();
      double Smoverlap_g2 = 0.0;
      double Smestimated_g2 = 0.0;
      if(ibin>51) Smestimated_g2 = inteSmsidelo_g2+inteSmsidehi_g2;
      else        Smestimated_g2 = inteSmsidelo2_g2;
      Smoverlap_g2 = inteSm_g2-Smestimated_g2;
      if(Smestimated_g2<0.0) Smestimated_g2 = 0.0;
      if(Smoverlap_g2<0.0){
        Smoverlap_g2 = 0.0; 
        Smestimated_g2 = 0.0;
      }
      double inteK0_g2 = IMnpim_wK0_n_sub_bin_select[ibin][iq]->Integral();
      double inteK0sidelo_g2 = IMnpim_wK0_n_sub_bin_sidelo[ibin][iq]->Integral();
      double inteK0sidehi_g2 = IMnpim_wK0_n_sub_bin_sidehi[ibin][iq]->Integral();
      double K0estimated_g2 = 0.0;
      if(ibin==44 || ibin==45) K0estimated_g2 = inteK0sidehi_g2*2.0;
      else K0estimated_g2 = inteK0sidelo_g2 + inteK0sidehi_g2;
      double K0overlap_g2 = inteK0_g2 - K0estimated_g2;
      if(K0estimated_g2 <0.0) K0estimated_g2 = 0.0;
      if(K0overlap_g2<0.0) {
         K0overlap_g2 = 0.0;
         K0estimated_g2 = 0.0;
      }
      pt2->AddText(Form("IM(n#pi^{-}#pi^{+})  %0.2f-%0.2f %s",1.0+ibin*1.0/nbintemplate,1.0+(ibin+1.0)/nbintemplate,cqcut[iq])); 
      pt2->AddText(Form("Sigma- count %0.2f ",inteSm_g2));// 
      if(ibin>51){
        pt2->AddText(Form("Sigma- side low  %0.2f ",inteSmsidelo_g2));
        pt2->AddText(Form("Sigma- side high %0.2f ",inteSmsidehi_g2));
      }else{
        pt2->AddText(Form("Sigma- side low (4 sigma)  %0.2f ",inteSmsidelo2_g2));
      }
      pt2->AddText(Form("Sigma- estimated  %0.2f ",Smestimated_g2)); //Sm in (K0 ^ Sm)     
      pt2->AddText(Form("K0 count %0.2f ",inteK0_g2));
      if(ibin==44 || ibin==45){
        pt2->AddText(Form("K0 side high*2  %0.2f ",inteK0sidehi_g2*2.0));
      }else{
        pt2->AddText(Form("K0 side low   %0.2f ",inteK0sidelo_g2));
        pt2->AddText(Form("K0 side high  %0.2f ",inteK0sidehi_g2));
      }
      pt2->AddText(Form("K0 estimated  %0.2f ", K0estimated_g2));
      pt2->Draw();
      IMnpipi_overlapdeco_K0[1][iq]->SetBinContent(ibin,K0estimated_g2);
      IMnpipi_overlapdeco_K0[1][iq]->SetBinError(ibin,sqrt(K0estimated_g2));
      IMnpipi_overlapdeco_Sm[1][iq]->SetBinContent(ibin,Smestimated_g2);
      IMnpipi_overlapdeco_Sm[1][iq]->SetBinError(ibin,sqrt(Smestimated_g2));
      
      /*
      //first canvas to display IM(npi-) vs IM(npi+) and their simple projection just to see the signal and other background 
      //no complicated cut at this moment
      cIMnpim_IMnpip_n_sub_bin[ibin][iq] = new TCanvas(Form("cIMnpim_IMnpip_n_sub_bin%d_%d",ibin,iq),Form("cIMnpim_IMnpip_n_sub_bin%d_%d",ibin,iq),800,800);
      cIMnpim_IMnpip_n_sub_bin[ibin][iq]->Divide(2,2,0,0);
      cIMnpim_IMnpip_n_sub_bin[ibin][iq]->cd(3);
      IMnpim_IMnpip_n_bin_sub[ibin][iq]->RebinX(4);
      IMnpim_IMnpip_n_bin_sub[ibin][iq]->RebinY(4);
      IMnpim_IMnpip_n_bin_sub[ibin][iq]->Draw("colz");
      cIMnpim_IMnpip_n_sub_bin[ibin][iq]->cd(1);
      IMnpip_n_sub_bin_2[ibin][iq] = (TH1D*)IMnpim_IMnpip_n_bin_sub[ibin][iq]->ProjectionX(Form("IMnpip_n_sub_bin_2_%d_%d",ibin,iq));
      IMnpip_n_sub_bin_2[ibin][iq]->SetTitle(Form("IMnpip_n_sub_bin_2_%d_%d",ibin,iq));
      IMnpip_n_sub_bin_2[ibin][iq]->Draw("HE");

      cIMnpim_IMnpip_n_sub_bin[ibin][iq]->cd(4);
      IMnpim_n_sub_bin_2[ibin][iq] = (TH1D*)IMnpim_IMnpip_n_bin_sub[ibin][iq]->ProjectionY(Form("IMnpim_n_sub_bin_2_%d_%d",ibin,iq));
      gr_IMnpim_n_sub_bin_2[ibin][iq] = new TGraphErrors();
      HistToRorateGraph(IMnpim_n_sub_bin_2[ibin][iq],*gr_IMnpim_n_sub_bin_2[ibin][iq]);
      gr_IMnpim_n_sub_bin_2[ibin][iq]->Draw("AP");
      */

      //2nd canvas to display the calculation of decomposion of K0, Sigma+,Sigma- 
      //2D plot in canvas(3) has cut to select Sigma- or Sigma+
      //projection plot in canvas(1) select Sigma+ events candidate + (Sigma- & Sigma+ ) overlap events
      //projection plot in canvas(4) select Sigma- events candidate + (Sigma- & Sigma- ) overlap events 
      cIMnpim_IMnpip_n_sub_bin_cut[ibin][iq]
      = new TCanvas(Form("cIMnpim_IMnpip_n_sub_bin_cut_%d_%d",ibin,iq),Form("cIMnpim_IMnpip_n_sub_bin_cut%d_%d",ibin,iq),800,800);
      cIMnpim_IMnpip_n_sub_bin_cut[ibin][iq]->Divide(2,2);
      cIMnpim_IMnpip_n_sub_bin_cut[ibin][iq]->cd(3);

      IMnpim_IMnpip_wK0orwSid_n_bin_sub[ibin][iq]->RebinX(4);
      IMnpim_IMnpip_wK0orwSid_n_bin_sub[ibin][iq]->RebinY(4);
      IMnpim_IMnpip_wK0orwSid_n_bin_sub[ibin][iq]->GetXaxis()->SetRangeUser(1,1.7);
      IMnpim_IMnpip_wK0orwSid_n_bin_sub[ibin][iq]->GetYaxis()->SetRangeUser(1,1.7);
      IMnpim_IMnpip_wK0orwSid_n_bin_sub[ibin][iq]->Draw("colz");
      
      cIMnpim_IMnpip_n_sub_bin_cut[ibin][iq]->cd(1);
      IMnpip_wSid_n_Sm_sub_bin[ibin][iq] = (TH1D*)IMnpim_IMnpip_wSid_n_Sm_bin_sub[ibin][iq]->ProjectionX(Form("IMnpip_wSid_n_Sm_sub_bin%d_%d",ibin,iq));
      IMnpip_wSid_n_Sm_sub_bin[ibin][iq]->SetTitle(Form("IMnpip_wSid_n_Sm_sub_bin%d_%d",ibin,iq));
      if(RebinMode){
        IMnpip_wSid_n_Sm_sub_bin[ibin][iq]->RebinX(5);
      }
      IMnpip_wSid_n_Sm_sub_bin[ibin][iq]->Draw("HIST");

      IMnpip_wSid_n_Sm_sub_bin_select[ibin][iq] = (TH1D*)IMnpip_wSid_n_Sm_sub_bin[ibin][iq]->Clone(Form("IMnpip_wSid_n_Sm_sub_bin%d_%d_select",ibin,iq));
      IMnpip_wSid_n_Sm_sub_bin_select[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmap_MIN,anacuts::Sigmap_MAX);
      IMnpip_wSid_n_Sm_sub_bin_select[ibin][iq]->SetLineColor(2);
      //IMnpip_wSid_n_Sm_sub_bin_select[ibin][iq]->SetFillColor(2);
      IMnpip_wSid_n_Sm_sub_bin_select[ibin][iq]->Draw("HISTsame");
      IMnpip_wSid_n_Sm_sub_bin_sidelo[ibin][iq] = (TH1D*)IMnpip_wSid_n_Sm_sub_bin[ibin][iq]->Clone(Form("IMnpip_wSid_n_Sm_sub_bin%d_%d_sidelo",ibin,iq));
      IMnpip_wSid_n_Sm_sub_bin_sidelo2[ibin][iq] = (TH1D*)IMnpip_wSid_n_Sm_sub_bin[ibin][iq]->Clone(Form("IMnpip_wSid_n_Sm_sub_bin%d_%d_sidelo2",ibin,iq));
      IMnpip_wSid_n_Sm_sub_bin_sidehi[ibin][iq] = (TH1D*)IMnpip_wSid_n_Sm_sub_bin[ibin][iq]->Clone(Form("IMnpip_wSid_n_Sm_sub_bin%d_%d_sidehi",ibin,iq));
      IMnpip_wSid_n_Sm_sub_bin_sidehi2[ibin][iq] = (TH1D*)IMnpip_wSid_n_Sm_sub_bin[ibin][iq]->Clone(Form("IMnpip_wSid_n_Sm_sub_bin%d_%d_sidehi2",ibin,iq));
      IMnpip_wSid_n_Sm_sub_bin_sidelo[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmap_MIN-4.0*anacuts::Sigmap_sigma,anacuts::Sigmap_MIN-2.0*anacuts::Sigmap_sigma);
      IMnpip_wSid_n_Sm_sub_bin_sidelo2[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmap_MIN-6.0*anacuts::Sigmap_sigma,anacuts::Sigmap_MIN-2.0*anacuts::Sigmap_sigma);
      if(ibin==42){
        IMnpip_wSid_n_Sm_sub_bin_sidehi[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmap_MAX+1.0*anacuts::Sigmap_sigma,anacuts::Sigmap_MAX+3.0*anacuts::Sigmap_sigma);
      }else{
        IMnpip_wSid_n_Sm_sub_bin_sidehi[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmap_MAX+2.0*anacuts::Sigmap_sigma,anacuts::Sigmap_MAX+4.0*anacuts::Sigmap_sigma);
      }
      IMnpip_wSid_n_Sm_sub_bin_sidehi2[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmap_MAX+2.0*anacuts::Sigmap_sigma,anacuts::Sigmap_MAX+6.0*anacuts::Sigmap_sigma);
      IMnpip_wSid_n_Sm_sub_bin_sidelo[ibin][iq]->SetLineColor(4);
      //IMnpip_wSid_n_Sm_sub_bin_sidelo[ibin][iq]->SetFillColor(4);
      IMnpip_wSid_n_Sm_sub_bin_sidelo2[ibin][iq]->SetLineColor(4);
      //IMnpip_wSid_n_Sm_sub_bin_sidelo2[ibin][iq]->SetFillColor(4);
      IMnpip_wSid_n_Sm_sub_bin_sidehi[ibin][iq]->SetLineColor(4);
      //IMnpip_wSid_n_Sm_sub_bin_sidehi[ibin][iq]->SetFillColor(4);
      IMnpip_wSid_n_Sm_sub_bin_sidehi2[ibin][iq]->SetLineColor(4);
      //IMnpip_wSid_n_Sm_sub_bin_sidehi2[ibin][iq]->SetFillColor(4);
      if(ibin<42){
        IMnpip_wSid_n_Sm_sub_bin_sidelo2[ibin][iq]->Draw("HISTsame");
      }else if(ibin==42){
        IMnpip_wSid_n_Sm_sub_bin_sidehi[ibin][iq]->Draw("HISTsame");
      }else if(42 < ibin && ibin<44){
        IMnpip_wSid_n_Sm_sub_bin_sidelo[ibin][iq]->Draw("HISTsame");
        IMnpip_wSid_n_Sm_sub_bin_sidehi[ibin][iq]->Draw("HISTsame");
      }else if(44 <= ibin){
        IMnpip_wSid_n_Sm_sub_bin_sidehi2[ibin][iq]->Draw("HISTsame");
      }
      cIMnpim_IMnpip_n_sub_bin_cut[ibin][iq]->cd(4);
      IMnpim_wSid_n_Sp_sub_bin[ibin][iq]
      = (TH1D*)IMnpim_IMnpip_wSid_n_Sp_bin_sub[ibin][iq]->ProjectionY(Form("IMnpim_wSid_n_Sp_sub_bin%d_%d",ibin,iq));
      IMnpim_wSid_n_Sp_sub_bin[ibin][iq]->SetTitle(Form("IMnpim_wSid_n_Sp_sub_bin%d_%d",ibin,iq));
      if(RebinMode)IMnpim_wSid_n_Sp_sub_bin[ibin][iq]->RebinX(5);
      IMnpim_wSid_n_Sp_sub_bin[ibin][iq]->Draw("HIST");
      IMnpim_wSid_n_Sp_sub_bin_select[ibin][iq] = (TH1D*)IMnpim_wSid_n_Sp_sub_bin[ibin][iq]->Clone(Form("IMnpim_wSid_n_Sp_sub_bin%d_%d_select",ibin,iq));
      IMnpim_wSid_n_Sp_sub_bin_select[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmam_MIN,anacuts::Sigmam_MAX);
      IMnpim_wSid_n_Sp_sub_bin_select[ibin][iq]->SetLineColor(2);
      //IMnpim_wSid_n_Sp_sub_bin_select[ibin][iq]->SetFillColor(2);
      IMnpim_wSid_n_Sp_sub_bin_select[ibin][iq]->Draw("HISTsame");
      IMnpim_wSid_n_Sp_sub_bin_sidelo[ibin][iq] = (TH1D*)IMnpim_wSid_n_Sp_sub_bin[ibin][iq]->Clone(Form("IMnpim_wSid_n_Sp_sub_bin%d_%d_sidelo",ibin,iq));
      IMnpim_wSid_n_Sp_sub_bin_sidelo2[ibin][iq] = (TH1D*)IMnpim_wSid_n_Sp_sub_bin[ibin][iq]->Clone(Form("IMnpim_wSid_n_Sp_sub_bin%d_%d_sidelo2",ibin,iq));
      IMnpim_wSid_n_Sp_sub_bin_sidehi[ibin][iq] = (TH1D*)IMnpim_wSid_n_Sp_sub_bin[ibin][iq]->Clone(Form("IMnpim_wSid_n_Sp_sub_bin%d_%d_sidehi",ibin,iq));
      IMnpim_wSid_n_Sp_sub_bin_sidehi2[ibin][iq] = (TH1D*)IMnpim_wSid_n_Sp_sub_bin[ibin][iq]->Clone(Form("IMnpim_wSid_n_Sp_sub_bin%d_%d_sidehi2",ibin,iq));
      IMnpim_wSid_n_Sp_sub_bin_sidelo[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmam_MIN-4.0*anacuts::Sigmam_sigma,anacuts::Sigmam_MIN-2.0*anacuts::Sigmam_sigma);
      IMnpim_wSid_n_Sp_sub_bin_sidelo2[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmam_MIN-6.0*anacuts::Sigmam_sigma,anacuts::Sigmam_MIN-2.0*anacuts::Sigmam_sigma);
      if(ibin==42){
        IMnpim_wSid_n_Sp_sub_bin_sidehi[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmam_MAX+1.0*anacuts::Sigmam_sigma,anacuts::Sigmam_MAX+3.0*anacuts::Sigmam_sigma);
      }else{
        IMnpim_wSid_n_Sp_sub_bin_sidehi[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmam_MAX+2.0*anacuts::Sigmam_sigma,anacuts::Sigmam_MAX+4.0*anacuts::Sigmam_sigma);
      }  
      IMnpim_wSid_n_Sp_sub_bin_sidehi2[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmam_MAX+2.0*anacuts::Sigmam_sigma,anacuts::Sigmam_MAX+6.0*anacuts::Sigmam_sigma);
      IMnpim_wSid_n_Sp_sub_bin_sidelo[ibin][iq]->SetLineColor(4);
      //IMnpim_wSid_n_Sp_sub_bin_sidelo[ibin][iq]->SetFillColor(4);
      IMnpim_wSid_n_Sp_sub_bin_sidelo2[ibin][iq]->SetLineColor(4);
      //IMnpim_wSid_n_Sp_sub_bin_sidelo2[ibin][iq]->SetFillColor(4);
      IMnpim_wSid_n_Sp_sub_bin_sidehi[ibin][iq]->SetLineColor(4);
      //IMnpim_wSid_n_Sp_sub_bin_sidehi[ibin][iq]->SetFillColor(4);
      IMnpim_wSid_n_Sp_sub_bin_sidehi2[ibin][iq]->SetLineColor(4);
      //IMnpim_wSid_n_Sp_sub_bin_sidehi2[ibin][iq]->SetFillColor(4);
      if(ibin <42){
        IMnpim_wSid_n_Sp_sub_bin_sidelo2[ibin][iq]->Draw("HISTsame");
      }else if( ibin == 42){
        IMnpim_wSid_n_Sp_sub_bin_sidehi[ibin][iq]->Draw("HISTsame");
      }else if(42 < ibin && ibin<44){
        IMnpim_wSid_n_Sp_sub_bin_sidelo[ibin][iq]->Draw("HISTsame");
        IMnpim_wSid_n_Sp_sub_bin_sidehi[ibin][iq]->Draw("HISTsame");
      }else if(44 <= ibin){
        IMnpim_wSid_n_Sp_sub_bin_sidehi2[ibin][iq]->Draw("HISTsame");
      }

      cIMnpim_IMnpip_n_sub_bin_cut[ibin][iq]->cd(2);
      TPaveText *pt3 = new TPaveText(.05,.05,.95,.7);
      double inteSp_g3 = IMnpim_wSid_n_Sp_sub_bin_select[ibin][iq]->Integral();
      double inteSpsidelo_g3 = IMnpim_wSid_n_Sp_sub_bin_sidelo[ibin][iq]->Integral();
      double inteSpsidelo2_g3 = IMnpim_wSid_n_Sp_sub_bin_sidelo2[ibin][iq]->Integral();
      double inteSpsidehi_g3 = IMnpim_wSid_n_Sp_sub_bin_sidehi[ibin][iq]->Integral();
      double inteSpsidehi2_g3 = IMnpim_wSid_n_Sp_sub_bin_sidehi2[ibin][iq]->Integral();
      double Spestimated_g3 = 0.0;
      double Spoverlap_g3 = 0.0;
      
      if(ibin<42) Spestimated_g3 = inteSpsidelo2_g3;
      else if(ibin==42) Spestimated_g3 = inteSpsidehi_g3*2.0;
      else if(42<ibin && ibin<44) Spestimated_g3 = inteSpsidelo_g3+inteSpsidehi_g3; 
      else Spestimated_g3 = inteSpsidehi2_g3;     
      if(Spestimated_g3<0.0) Spestimated_g3 = 0.0;
      Spoverlap_g3 = inteSp_g3 - Spestimated_g3;
      if(Spoverlap_g3<0.0){
        Spoverlap_g3 = 0.0;
        Spestimated_g3 = 0.0;
      }
      double inteSm_g3 = IMnpip_wSid_n_Sm_sub_bin_select[ibin][iq]->Integral();
      double inteSmsidelo_g3 = IMnpip_wSid_n_Sm_sub_bin_sidelo[ibin][iq]->Integral();
      double inteSmsidelo2_g3 = IMnpip_wSid_n_Sm_sub_bin_sidelo2[ibin][iq]->Integral();
      double inteSmsidehi_g3 = IMnpip_wSid_n_Sm_sub_bin_sidehi[ibin][iq]->Integral();
      double inteSmsidehi2_g3 = IMnpip_wSid_n_Sm_sub_bin_sidehi2[ibin][iq]->Integral();
      double Smestimated_g3 = 0.0;
      double Smoverlap_g3 = 0.0;
      
      if(ibin<42) Smestimated_g3 = inteSmsidelo2_g3;
      else if(ibin==42) Smestimated_g3 = inteSmsidehi_g3*2.0;
      else if(42< ibin &&  ibin<44) Smestimated_g3 = inteSmsidehi_g3 + inteSmsidelo_g3;
      else Smestimated_g3 = inteSmsidehi2_g3;
      if(Smestimated_g3<0.0) Smestimated_g3 = 0.0;
      Smoverlap_g3 = inteSm_g3 - Smestimated_g3;
      if(Smoverlap_g3<0.0) {
        Smoverlap_g3 = 0.0;
        Smestimated_g3 = 0.0;
      }
      pt3->AddText(Form("IM(n#pi^{-}#pi^{+})  %0.2f-%0.2f %s",1.0+ibin*1.0/nbintemplate,1.0+(ibin+1.0)/nbintemplate,cqcut[iq])); 
      pt3->AddText(Form("Sigma+ count %0.2f ",inteSp_g3));
      if(ibin<42){
        pt3->AddText(Form("Sigma+ side low (4 sigma)   %0.2f ",inteSpsidelo2_g3));
      }else if(ibin==42){
        pt3->AddText(Form("Sigma+ side high*2   %0.2f ",inteSpsidehi_g3*2.0));
      }else if(42 < ibin && ibin<44){
        pt3->AddText(Form("Sigma+ side low   %0.2f ",inteSpsidelo_g3));
        pt3->AddText(Form("Sigma+ side high  %0.2f ",inteSpsidehi_g3));
      }else{
        pt3->AddText(Form("Sigma+ side high (4 sigma)   %0.2f ",inteSpsidehi2_g3));
      } 
      pt3->AddText(Form("Sigma+ estimated %0.2f ", Spestimated_g3));
      pt3->AddText(Form("Sigma- count %0.2f ",inteSm_g3));
      if(ibin <42){
        pt3->AddText(Form("Sigma- side low (4 sigma)   %0.2f ",inteSmsidelo2_g3));
      }else if(ibin==42){
        pt3->AddText(Form("Sigma- side high*2   %0.2f ",inteSmsidehi_g3*2.0));
      }else if(42 < ibin && ibin<44){
        pt3->AddText(Form("Sigma- side low   %0.2f ",inteSmsidelo_g3));
        pt3->AddText(Form("Sigma- side high   %0.2f ",inteSmsidehi_g3));
      }else{
        pt3->AddText(Form("Sigma- side high (4 sigma)  %0.2f ",inteSmsidehi2_g3));
      }
      pt3->AddText(Form("Sigma- estimated %0.2f ",Smestimated_g3)); 
      pt3->Draw();

      IMnpipi_overlapdeco_Sm[2][iq]->SetBinContent(ibin,Smestimated_g3);
      IMnpipi_overlapdeco_Sm[2][iq]->SetBinError(ibin,sqrt(Smestimated_g3));
      IMnpipi_overlapdeco_Sp[2][iq]->SetBinContent(ibin,Spestimated_g3);
      IMnpipi_overlapdeco_Sp[2][iq]->SetBinError(ibin,sqrt(Spestimated_g3));
    }
  }
   
  TCanvas *csum_Sp[ngroup][nqcut];
  TCanvas *csum_Sm[ngroup][nqcut];
  TCanvas *csum_K0[ngroup][nqcut];
  for(int ig=0;ig<ngroup;ig++){
    for(int iq=qstart;iq<nqcut;iq++){
      if(ig!=2){
        csum_K0[ig][iq] = new TCanvas(Form("csum_K0_g%d_q%d",ig,iq),Form("csum_K0_g%d_q%d",ig,iq));
        IMnpipi_overlapdeco_K0[ig][iq]->Draw("H");
      }
      if(ig!=1){  
        csum_Sp[ig][iq] = new TCanvas(Form("csum_Sp_g%d_q%d",ig,iq),Form("csum_Sp_g%d_q%d",ig,iq));
        IMnpipi_overlapdeco_Sp[ig][iq]->Draw("H");
      }
      if(ig!=0){
        csum_Sm[ig][iq] = new TCanvas(Form("csum_Sm_g%d_q%d",ig,iq),Form("csum_Sm_g%d_q%d",ig,iq));
        IMnpipi_overlapdeco_Sm[ig][iq]->Draw("H");
      }
    }
  }

  TCanvas *cratio_SpSm[nqcut];//Sp/Sm  ->group 2
  TCanvas *cratio_SpK0[nqcut];//Sp/K0  ->group 0
  TCanvas *cratio_SmK0[nqcut];//Sm/K0  ->group 1
  TH1D* IMnpipi_SpSm_ratio[nqcut];
  TH1D* IMnpipi_SpK0_ratio[nqcut];
  TH1D* IMnpipi_SmK0_ratio[nqcut];

  for(int iq=qstart;iq<nqcut;iq++){
    IMnpipi_SpSm_ratio[iq] = new TH1D(Form("IMnpipi_SpSm_ratio_%d",iq),Form("IMnpipi_SpSm_ratio_%s",cqcut[iq]),100,1,2);
    IMnpipi_SpK0_ratio[iq] = new TH1D(Form("IMnpipi_SpK0_ratio_%d",iq),Form("IMnpipi_SpK0_ratio_%s",cqcut[iq]),100,1,2);
    IMnpipi_SmK0_ratio[iq] = new TH1D(Form("IMnpipi_SmK0_ratio_%d",iq),Form("IMnpipi_SmK0_ratio_%s",cqcut[iq]),100,1,2);
  
    IMnpipi_SpSm_ratio[iq] = (TH1D*)IMnpipi_overlapdeco_Sp[2][iq]->Clone(Form("IMnpipi_SpSm_ratio_%d",iq));
    IMnpipi_SpSm_ratio[iq]->Divide(IMnpipi_overlapdeco_Sm[2][iq]);
    cratio_SpSm[iq] = new TCanvas(Form("cratrio_SpSm_%s",cqcut[iq]));
    IMnpipi_SpSm_ratio[iq]->SetTitle(Form("IMnpipi_Sp/Sm_ratio_%s",cqcut[iq]));
    IMnpipi_SpSm_ratio[iq]->Draw("HIST");
    
    IMnpipi_SpK0_ratio[iq] = (TH1D*)IMnpipi_overlapdeco_Sp[0][iq]->Clone(Form("IMnpipi_SpK0_ratio_%d",iq));
    IMnpipi_SpK0_ratio[iq]->Divide(IMnpipi_overlapdeco_K0[0][iq]);
    cratio_SpK0[iq] = new TCanvas(Form("cratrio_SpK0_%s",cqcut[iq]));
    IMnpipi_SpK0_ratio[iq]->SetTitle(Form("IMnpipi_Sp/K0_ratio_%s",cqcut[iq]));
    IMnpipi_SpK0_ratio[iq]->Draw("HIST");

    IMnpipi_SmK0_ratio[iq] = (TH1D*)IMnpipi_overlapdeco_Sm[1][iq]->Clone(Form("IMnpipi_SmK0_ratio_%d",iq));
    IMnpipi_SmK0_ratio[iq]->Divide(IMnpipi_overlapdeco_K0[1][iq]);
    cratio_SmK0[iq] = new TCanvas(Form("cratrio_SmK0_%s",cqcut[iq]));
    IMnpipi_SmK0_ratio[iq]->SetTitle(Form("IMnpipi_Sm/K0_ratio_%s",cqcut[iq]));
    IMnpipi_SmK0_ratio[iq]->Draw("HIST");
  }

  
  //merging some specific bins
  TH2D* IMpippim_IMnpip_wK0orwSid_n_merge[10][nqcut];
  TH2D* IMpippim_IMnpim_wK0orwSid_n_merge[10][nqcut];
  TH2D* IMnpim_IMnpip_wSid_n_merge[10][nqcut];
  TH1D* IMnpip_wK0_merge[10][nqcut];
  TH1D* IMnpip_wK0_merge_select[10][nqcut];
  TH1D* IMnpip_wK0_merge_lo[10][nqcut];
  TH1D* IMnpip_wK0_merge_hi[10][nqcut];
  TH1D* IMpippim_Sp_merge[10][nqcut];
  TH1D* IMpippim_Sp_merge_select[10][nqcut];
  TH1D* IMpippim_Sp_merge_lo[10][nqcut];
  TH1D* IMpippim_Sp_merge_hi[10][nqcut];
  TH1D* IMnpim_wK0_merge[10][nqcut];
  TH1D* IMnpim_wK0_merge_select[10][nqcut];
  TH1D* IMnpim_wK0_merge_lo[10][nqcut];
  TH1D* IMnpim_wK0_merge_hi[10][nqcut];
  TH1D* IMpippim_Sm_merge[10][nqcut];
  TH1D* IMpippim_Sm_merge_select[10][nqcut];
  TH1D* IMpippim_Sm_merge_lo[10][nqcut];
  TH1D* IMpippim_Sm_merge_hi[10][nqcut];
  TH1D* IMnpip_Sm_merge[10][nqcut];
  TH1D* IMnpip_Sm_merge_select[10][nqcut];
  TH1D* IMnpip_Sm_merge_lo[10][nqcut];
  TH1D* IMnpip_Sm_merge_hi[10][nqcut];
  TH1D* IMnpim_Sp_merge[10][nqcut];
  TH1D* IMnpim_Sp_merge_select[10][nqcut];
  TH1D* IMnpim_Sp_merge_lo[10][nqcut];
  TH1D* IMnpim_Sp_merge_hi[10][nqcut];
  
  TCanvas *cIMpippim_IMnpip_merge[10][nqcut];
  TCanvas *cIMpippim_IMnpim_merge[10][nqcut];
  TCanvas *cIMnpim_IMnpip_merge[10][nqcut];
  
  //merge
  //index 0 
  //1.40-1.52
  //
  //index 1 
  //1.52-2.00
  const int initbin[10]={40,52};
  const int startbin[10]={41,53};
  const int endbin[10]={52,100};
  const float startval[10]={1.40,1.52};
  const float endval[10]={1.52,2.00};
  for(int imerge=0;imerge<2;imerge++){
    for(int iqcut=0;iqcut<nqcut;iqcut++){
      //index 0 group
      //use M = 1.40-1.41 GeV bin for cloneing 
      IMpippim_IMnpip_wK0orwSid_n_merge[imerge][iqcut]  
        = (TH2D*)IMpippim_IMnpip_wK0orwSid_n_bin_sub[initbin[imerge]][iqcut]->Clone(Form("IMpippim_IMnpip_wK0orwSid_merge%d_%d",imerge,iqcut));
      IMpippim_IMnpip_wK0orwSid_n_merge[imerge][iqcut]->SetTitle(Form("IMpippim_IMnpip_wK0orwSid %0.2f-%0.2f %s",startval[imerge],endval[imerge],cqcut[iqcut]));
      IMpippim_IMnpim_wK0orwSid_n_merge[imerge][iqcut] 
        = (TH2D*)IMpippim_IMnpim_wK0orwSid_n_bin_sub[initbin[imerge]][iqcut]->Clone(Form("IMpippim_IMnpim_wK0orwSid_merge%d_%d",imerge,iqcut));
      IMpippim_IMnpim_wK0orwSid_n_merge[imerge][iqcut]->SetTitle(Form("IMpippim_IMnpim_wK0orwSid %0.2f-%0.2f %s",startval[imerge],endval[imerge],cqcut[iqcut]));
      IMnpim_IMnpip_wSid_n_merge[imerge][iqcut]
        = (TH2D*)IMnpim_IMnpip_wK0orwSid_n_bin_sub[initbin[imerge]][iqcut]->Clone(Form("IMnpim_IMnpip_wSid_n_merge%d_%d",imerge,iqcut));
      IMnpim_IMnpip_wSid_n_merge[imerge][iqcut]->SetTitle(Form("IMnpim_IMnpip_wK0orwSid  %0.2f-%0.2f %s",startval[imerge],endval[imerge],cqcut[iqcut]));
      IMnpip_wK0_merge[imerge][iqcut] = (TH1D*)IMnpip_wK0_n_sub_bin[initbin[imerge]][iqcut]->Clone(Form("IMnpip_wK0_merge%d_%d",imerge,iqcut));
      IMpippim_Sp_merge[imerge][iqcut] = (TH1D*)IMpippim_wSid_n_Sp_sub_bin[initbin[imerge]][iqcut]->Clone(Form("IMpippim_Sp_merge%d_%d",imerge,iqcut));
      IMnpim_wK0_merge[imerge][iqcut] = (TH1D*)IMnpim_wK0_n_sub_bin[initbin[imerge]][iqcut]->Clone(Form("IMnpim_wK0_merge%d_%d",imerge,iqcut));
      IMpippim_Sm_merge[imerge][iqcut] = (TH1D*)IMpippim_wSid_n_Sm_sub_bin[initbin[imerge]][iqcut]->Clone(Form("IMpippim_Sm_merge%d_%d",imerge,iqcut));
      IMnpip_Sm_merge[imerge][iqcut] = (TH1D*)IMnpip_wSid_n_Sm_sub_bin[initbin[imerge]][iqcut]->Clone(Form("IMnpip_Sm_merge%d_%d",imerge,iqcut));
      IMnpim_Sp_merge[imerge][iqcut] = (TH1D*)IMnpim_wSid_n_Sp_sub_bin[initbin[imerge]][iqcut]->Clone(Form("IMnpim_Sp_merge%d_%d",imerge,iqcut));

      for(int ibin=startbin[imerge];ibin<endbin[imerge];ibin++){
        IMpippim_IMnpip_wK0orwSid_n_merge[imerge][iqcut]->Add(IMpippim_IMnpip_wK0orwSid_n_bin_sub[ibin][iqcut]);
        IMpippim_IMnpim_wK0orwSid_n_merge[imerge][iqcut]->Add(IMpippim_IMnpim_wK0orwSid_n_bin_sub[ibin][iqcut]);
        IMnpim_IMnpip_wSid_n_merge[imerge][iqcut]->Add(IMnpim_IMnpip_wK0orwSid_n_bin_sub[ibin][iqcut]);
        IMnpip_wK0_merge[imerge][iqcut]->Add(IMnpip_wK0_n_sub_bin[ibin][iqcut]);
        IMpippim_Sp_merge[imerge][iqcut]->Add(IMpippim_wSid_n_Sp_sub_bin[ibin][iqcut]);
        IMnpim_wK0_merge[imerge][iqcut]->Add(IMnpim_wK0_n_sub_bin[ibin][iqcut]);
        IMpippim_Sm_merge[imerge][iqcut]->Add(IMpippim_wSid_n_Sm_sub_bin[ibin][iqcut]);
        IMnpip_Sm_merge[imerge][iqcut]->Add(IMnpip_wSid_n_Sm_sub_bin[ibin][iqcut]);
        IMnpim_Sp_merge[imerge][iqcut]->Add(IMnpim_wSid_n_Sp_sub_bin[ibin][iqcut]);
      }
      
      //Draw in Canvas 
      cIMpippim_IMnpip_merge[imerge][iqcut] = new TCanvas(Form("cIMpippim_IMnpip_merge%d_%d",imerge,iqcut),Form("cIMpippim_IMnpip_merge%d_%d",imerge,iqcut));
      cIMpippim_IMnpip_merge[imerge][iqcut]->Divide(2,2);
      cIMpippim_IMnpip_merge[imerge][iqcut]->cd(3);
      if(RemoveNegative)IMpippim_IMnpip_wK0orwSid_n_merge[imerge][iqcut]->SetMinimum(0);
      IMpippim_IMnpip_wK0orwSid_n_merge[imerge][iqcut]->Draw("colz");
      cIMpippim_IMnpip_merge[imerge][iqcut]->cd(1);
      IMnpip_wK0_merge[imerge][iqcut]->Draw("HIST");
      IMnpip_wK0_merge_select[imerge][iqcut] = (TH1D*)IMnpip_wK0_merge[imerge][iqcut]->Clone(Form("IMnpip_wK0_merge_select_%d_%d",imerge,iqcut));
      IMnpip_wK0_merge_lo[imerge][iqcut] = (TH1D*)IMnpip_wK0_merge[imerge][iqcut]->Clone(Form("IMnpip_wK0_merge_lo_%d_%d",imerge,iqcut));
      IMnpip_wK0_merge_hi[imerge][iqcut] = (TH1D*)IMnpip_wK0_merge[imerge][iqcut]->Clone(Form("IMnpip_wK0_merge_hi_%d_%d",imerge,iqcut));
      IMnpip_wK0_merge_select[imerge][iqcut]->SetLineColor(2);
      IMnpip_wK0_merge_lo[imerge][iqcut]->SetLineColor(4);
      IMnpip_wK0_merge_hi[imerge][iqcut]->SetLineColor(4);
      IMnpip_wK0_merge_select[imerge][iqcut]->GetXaxis()->SetRangeUser(anacuts::Sigmap_MIN,anacuts::Sigmap_MAX);
      if(Sidefar){
        IMnpip_wK0_merge_lo[imerge][iqcut]->GetXaxis()->SetRangeUser(anacuts::Sigmap_MIN-4.0*anacuts::Sigmap_sigma, anacuts::Sigmap_MIN-2.0*anacuts::Sigmap_sigma);
        IMnpip_wK0_merge_hi[imerge][iqcut]->GetXaxis()->SetRangeUser(anacuts::Sigmap_MAX+2.0*anacuts::Sigmap_sigma, anacuts::Sigmap_MAX+4.0*anacuts::Sigmap_sigma);
      }else{
        IMnpip_wK0_merge_lo[imerge][iqcut]->GetXaxis()->SetRangeUser(anacuts::Sigmap_MIN-2.0*anacuts::Sigmap_sigma, anacuts::Sigmap_MIN);
        IMnpip_wK0_merge_hi[imerge][iqcut]->GetXaxis()->SetRangeUser(anacuts::Sigmap_MAX, anacuts::Sigmap_MAX+2.0*anacuts::Sigmap_sigma);
      }
      IMnpip_wK0_merge_select[imerge][iqcut]->Draw("HISTsame");
      IMnpip_wK0_merge_lo[imerge][iqcut]->Draw("HISTsame");
      IMnpip_wK0_merge_hi[imerge][iqcut]->Draw("HISTsame");
      cIMpippim_IMnpip_merge[imerge][iqcut]->cd(4);
      IMpippim_Sp_merge[imerge][iqcut]->Draw("HIST");
      IMpippim_Sp_merge_select[imerge][iqcut] = (TH1D*)IMpippim_Sp_merge[imerge][iqcut]->Clone(Form("IMpippim_Sp_merge_select_%d_%d",imerge,iqcut));
      IMpippim_Sp_merge_lo[imerge][iqcut] = (TH1D*)IMpippim_Sp_merge[imerge][iqcut]->Clone(Form("IMpippim_Sp_merge_lo_%d_%d",imerge,iqcut));
      IMpippim_Sp_merge_hi[imerge][iqcut] = (TH1D*)IMpippim_Sp_merge[imerge][iqcut]->Clone(Form("IMpippim_Sp_merge_hi_%d_%d",imerge,iqcut));
      IMpippim_Sp_merge_select[imerge][iqcut]->SetLineColor(2);
      IMpippim_Sp_merge_lo[imerge][iqcut]->SetLineColor(4);
      IMpippim_Sp_merge_hi[imerge][iqcut]->SetLineColor(4);
      IMpippim_Sp_merge_select[imerge][iqcut]->GetXaxis()->SetRangeUser(anacuts::pipi_MIN_narrow,anacuts::pipi_MAX_narrow);
      if(Sidefar){
        IMpippim_Sp_merge_lo[imerge][iqcut]->GetXaxis()->SetRangeUser(anacuts::pipi_MIN_narrow-4.0*anacuts::K0_sigma,anacuts::pipi_MIN_narrow-2.0*anacuts::K0_sigma);
        IMpippim_Sp_merge_hi[imerge][iqcut]->GetXaxis()->SetRangeUser(anacuts::pipi_MAX_narrow+2.0*anacuts::K0_sigma,anacuts::pipi_MAX_narrow+4.0*anacuts::K0_sigma);
      }else{
        IMpippim_Sp_merge_lo[imerge][iqcut]->GetXaxis()->SetRangeUser(anacuts::pipi_MIN_narrow-2.0*anacuts::K0_sigma,anacuts::pipi_MIN_narrow-0.0*anacuts::K0_sigma);
        IMpippim_Sp_merge_hi[imerge][iqcut]->GetXaxis()->SetRangeUser(anacuts::pipi_MAX_narrow+0.0*anacuts::K0_sigma,anacuts::pipi_MAX_narrow+2.0*anacuts::K0_sigma);
      }
      IMpippim_Sp_merge_select[imerge][iqcut]->Draw("HISTsame");
      IMpippim_Sp_merge_lo[imerge][iqcut]->Draw("HISTsame");
      IMpippim_Sp_merge_hi[imerge][iqcut]->Draw("HISTsame");

      cIMpippim_IMnpim_merge[imerge][iqcut] = new TCanvas(Form("cIMpippim_IMnpim_merge%d_%d",imerge,iqcut),Form("cIMpippim_IMnpim_merge%d_%d",imerge,iqcut));
      cIMpippim_IMnpim_merge[imerge][iqcut]->Divide(2,2);
      cIMpippim_IMnpim_merge[imerge][iqcut]->cd(3);
      if(RemoveNegative)IMpippim_IMnpim_wK0orwSid_n_merge[imerge][iqcut]->SetMinimum(0);
      IMpippim_IMnpim_wK0orwSid_n_merge[imerge][iqcut]->Draw("colz");
      cIMpippim_IMnpim_merge[imerge][iqcut]->cd(1);
      IMnpim_wK0_merge[imerge][iqcut]->Draw("HIST");
      IMnpim_wK0_merge_select[imerge][iqcut] = (TH1D*)IMnpim_wK0_merge[imerge][iqcut]->Clone(Form("IMnpim_wK0_merge_select_%d_%d",imerge,iqcut));
      IMnpim_wK0_merge_lo[imerge][iqcut] = (TH1D*)IMnpim_wK0_merge[imerge][iqcut]->Clone(Form("IMnpim_wK0_merge_lo_%d_%d",imerge,iqcut));
      IMnpim_wK0_merge_hi[imerge][iqcut] = (TH1D*)IMnpim_wK0_merge[imerge][iqcut]->Clone(Form("IMnpim_wK0_merge_hi_%d_%d",imerge,iqcut));
      IMnpim_wK0_merge_select[imerge][iqcut]->SetLineColor(2);
      IMnpim_wK0_merge_lo[imerge][iqcut]->SetLineColor(4);
      IMnpim_wK0_merge_hi[imerge][iqcut]->SetLineColor(4);
      IMnpim_wK0_merge_select[imerge][iqcut]->GetXaxis()->SetRangeUser(anacuts::Sigmam_MIN,anacuts::Sigmam_MAX);
      if(Sidefar){
        IMnpim_wK0_merge_lo[imerge][iqcut]->GetXaxis()->SetRangeUser(anacuts::Sigmam_MIN-4.0*anacuts::Sigmam_sigma, anacuts::Sigmam_MIN-2.0*anacuts::Sigmam_sigma);
        IMnpim_wK0_merge_hi[imerge][iqcut]->GetXaxis()->SetRangeUser(anacuts::Sigmam_MAX+2.0*anacuts::Sigmam_sigma, anacuts::Sigmam_MAX+4.0*anacuts::Sigmam_sigma);
      }else{
        IMnpim_wK0_merge_lo[imerge][iqcut]->GetXaxis()->SetRangeUser(anacuts::Sigmam_MIN-2.0*anacuts::Sigmam_sigma, anacuts::Sigmam_MIN-0.0*anacuts::Sigmam_sigma);
        IMnpim_wK0_merge_hi[imerge][iqcut]->GetXaxis()->SetRangeUser(anacuts::Sigmam_MAX+0.0*anacuts::Sigmam_sigma, anacuts::Sigmam_MAX+2.0*anacuts::Sigmam_sigma);
      }
      IMnpim_wK0_merge_select[imerge][iqcut]->Draw("HISTsame");
      IMnpim_wK0_merge_lo[imerge][iqcut]->Draw("HISTsame");
      IMnpim_wK0_merge_hi[imerge][iqcut]->Draw("HISTsame");
      cIMpippim_IMnpim_merge[imerge][iqcut]->cd(4);
      IMpippim_Sm_merge[imerge][iqcut]->Draw("HIST");
      IMpippim_Sm_merge_select[imerge][iqcut] = (TH1D*)IMpippim_Sm_merge[imerge][iqcut]->Clone(Form("IMpippim_Sm_merge_select_%d_%d",imerge,iqcut));
      IMpippim_Sm_merge_lo[imerge][iqcut] = (TH1D*)IMpippim_Sm_merge[imerge][iqcut]->Clone(Form("IMpippim_Sm_merge_lo_%d_%d",imerge,iqcut));
      IMpippim_Sm_merge_hi[imerge][iqcut] = (TH1D*)IMpippim_Sm_merge[imerge][iqcut]->Clone(Form("IMpippim_Sm_merge_hi_%d_%d",imerge,iqcut));
      IMpippim_Sm_merge_select[imerge][iqcut]->SetLineColor(2);
      IMpippim_Sm_merge_lo[imerge][iqcut]->SetLineColor(4);
      IMpippim_Sm_merge_hi[imerge][iqcut]->SetLineColor(4);
      IMpippim_Sm_merge_select[imerge][iqcut]->GetXaxis()->SetRangeUser(anacuts::pipi_MIN_narrow,anacuts::pipi_MAX_narrow);
      if(Sidefar){
        IMpippim_Sm_merge_lo[imerge][iqcut]->GetXaxis()->SetRangeUser(anacuts::pipi_MIN_narrow-4.0*anacuts::K0_sigma,anacuts::pipi_MIN_narrow-2.0*anacuts::K0_sigma);
        IMpippim_Sm_merge_hi[imerge][iqcut]->GetXaxis()->SetRangeUser(anacuts::pipi_MAX_narrow+2.0*anacuts::K0_sigma,anacuts::pipi_MAX_narrow+4.0*anacuts::K0_sigma);
      }else{
        IMpippim_Sm_merge_lo[imerge][iqcut]->GetXaxis()->SetRangeUser(anacuts::pipi_MIN_narrow-2.0*anacuts::K0_sigma,anacuts::pipi_MIN_narrow-0.0*anacuts::K0_sigma);
        IMpippim_Sm_merge_hi[imerge][iqcut]->GetXaxis()->SetRangeUser(anacuts::pipi_MAX_narrow+0.0*anacuts::K0_sigma,anacuts::pipi_MAX_narrow+2.0*anacuts::K0_sigma);
      }
        
      IMpippim_Sm_merge_select[imerge][iqcut]->Draw("HISTsame");
      IMpippim_Sm_merge_lo[imerge][iqcut]->Draw("HISTsame");
      IMpippim_Sm_merge_hi[imerge][iqcut]->Draw("HISTsame");
      

      cIMnpim_IMnpip_merge[imerge][iqcut] = new TCanvas(Form("cIMnpim_IMnpip_merge%d_%d",imerge,iqcut),Form("cIMnpim_IMnpip_merge%d_%d",imerge,iqcut));
      cIMnpim_IMnpip_merge[imerge][iqcut]->Divide(2,2);
      cIMnpim_IMnpip_merge[imerge][iqcut]->cd(3);
      if(RemoveNegative)IMnpim_IMnpip_wSid_n_merge[imerge][iqcut]->SetMinimum(0);
      IMnpim_IMnpip_wSid_n_merge[imerge][iqcut]->Draw("colz");
      cIMnpim_IMnpip_merge[imerge][iqcut]->cd(1);
      IMnpip_Sm_merge[imerge][iqcut]->Draw("HIST");
      IMnpip_Sm_merge_select[imerge][iqcut] = (TH1D*)IMnpip_Sm_merge[imerge][iqcut]->Clone(Form("IMnpip_Sm_merge_select_%d_%d",imerge,iqcut));
      IMnpip_Sm_merge_lo[imerge][iqcut] = (TH1D*)IMnpip_Sm_merge[imerge][iqcut]->Clone(Form("IMnpip_Sm_merge_lo_%d_%d",imerge,iqcut));
      IMnpip_Sm_merge_hi[imerge][iqcut] = (TH1D*)IMnpip_Sm_merge[imerge][iqcut]->Clone(Form("IMnpip_Sm_merge_hi_%d_%d",imerge,iqcut));
      IMnpip_Sm_merge_select[imerge][iqcut]->SetLineColor(2);
      IMnpip_Sm_merge_lo[imerge][iqcut]->SetLineColor(4);
      IMnpip_Sm_merge_hi[imerge][iqcut]->SetLineColor(4);
      IMnpip_Sm_merge_select[imerge][iqcut]->GetXaxis()->SetRangeUser(anacuts::Sigmap_MIN,anacuts::Sigmap_MAX);
      if(Sidefar){
        IMnpip_Sm_merge_lo[imerge][iqcut]->GetXaxis()->SetRangeUser(anacuts::Sigmap_MIN-4.0*anacuts::Sigmap_sigma, anacuts::Sigmap_MIN-2.0*anacuts::Sigmap_sigma);
        IMnpip_Sm_merge_hi[imerge][iqcut]->GetXaxis()->SetRangeUser(anacuts::Sigmap_MAX+2.0*anacuts::Sigmap_sigma, anacuts::Sigmap_MAX+4.0*anacuts::Sigmap_sigma);
      }
      IMnpip_Sm_merge_select[imerge][iqcut]->Draw("HISTsame");
      IMnpip_Sm_merge_lo[imerge][iqcut]->Draw("HISTsame");
      IMnpip_Sm_merge_hi[imerge][iqcut]->Draw("HISTsame");
      cIMnpim_IMnpip_merge[imerge][iqcut]->cd(4);
      IMnpim_Sp_merge[imerge][iqcut]->Draw("HIST");
      IMnpim_Sp_merge_select[imerge][iqcut] = (TH1D*)IMnpim_Sp_merge[imerge][iqcut]->Clone(Form("IMnpim_Sp_merge_select_%d_%d",imerge,iqcut));
      IMnpim_Sp_merge_lo[imerge][iqcut] = (TH1D*)IMnpim_Sp_merge[imerge][iqcut]->Clone(Form("IMnpim_Sp_merge_lo_%d_%d",imerge,iqcut));
      IMnpim_Sp_merge_hi[imerge][iqcut] = (TH1D*)IMnpim_Sp_merge[imerge][iqcut]->Clone(Form("IMnpim_Sp_merge_hi_%d_%d",imerge,iqcut));
      IMnpim_Sp_merge_select[imerge][iqcut]->SetLineColor(2);
      IMnpim_Sp_merge_lo[imerge][iqcut]->SetLineColor(4);
      IMnpim_Sp_merge_hi[imerge][iqcut]->SetLineColor(4);
      IMnpim_Sp_merge_select[imerge][iqcut]->GetXaxis()->SetRangeUser(anacuts::Sigmam_MIN,anacuts::Sigmam_MAX);
      IMnpim_Sp_merge_lo[imerge][iqcut]->GetXaxis()->SetRangeUser(anacuts::Sigmam_MIN-4.0*anacuts::Sigmam_sigma, anacuts::Sigmam_MIN-2.0*anacuts::Sigmam_sigma);
      IMnpim_Sp_merge_hi[imerge][iqcut]->GetXaxis()->SetRangeUser(anacuts::Sigmam_MAX+2.0*anacuts::Sigmam_sigma, anacuts::Sigmam_MAX+4.0*anacuts::Sigmam_sigma);
      IMnpim_Sp_merge_select[imerge][iqcut]->Draw("HISTsame");
      IMnpim_Sp_merge_lo[imerge][iqcut]->Draw("HISTsame");
      IMnpim_Sp_merge_hi[imerge][iqcut]->Draw("HISTsame");
    }
  }



  TCanvas *c = NULL;
  TSeqCollection *SCol = gROOT->GetListOfCanvases();
  int size = SCol->GetSize();
  TIter next(SCol);
  TString pdfname = "K0SigmaTemp.pdf";
  if(RebinMode) pdfname = "K0SigmaTemp_Rebin.pdf";
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
  
  if(!RebinMode){
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
  }

  return;
  

}
