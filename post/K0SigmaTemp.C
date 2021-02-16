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
  fmix[0] = TFile::Open("evanaIMpisigma_npippim_v202_MIX_cut4_out_iso.root");
  fr[1] = TFile::Open("evanaIMpisigma_npippim_v202_out_iso_qlo.root");
  fmix[1] = TFile::Open("evanaIMpisigma_npippim_v202_MIX_cut4_out_iso_qlo.root");
  fr[2] = TFile::Open("evanaIMpisigma_npippim_v202_out_iso_qhi.root");
  fmix[2] = TFile::Open("evanaIMpisigma_npippim_v202_MIX_cut4_out_iso_qhi.root");
  fr[0]->Print();
  fmix[0]->Print();
  fr[1]->Print();
  fmix[1]->Print();
  fr[2]->Print();
  fmix[2]->Print();
  
  gROOT->SetBatch(1);
  gStyle->SetPalette(1);
  //gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0.);  

  const unsigned int nbintemplate = 100;
  const int nqcut=3;
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
  

  const char cqcut[][6]= {"all","qlo","qhi"};
  for(int iqcut=0;iqcut<nqcut;iqcut++){
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
    
    IMpippim_IMnpim_wK0orwSid_n_data[iqcut] = (TH2F*)fr[iqcut]->Get("IMpippim_IMnpim_wK0orwSid_n");
    IMpippim_IMnpim_wK0orwSid_n_mix[iqcut] = (TH2F*)fmix[iqcut]->Get("IMpippim_IMnpim_wK0orwSid_n");
    IMpippim_IMnpim_wK0orwSid_n_sub[iqcut] = (TH2F*)IMpippim_IMnpim_wK0orwSid_n_data[iqcut]->Clone(Form("IMpippim_IMnpim_wK0orwSid_n_sub_%s",cqcut[iqcut]));
    IMpippim_IMnpim_wK0orwSid_n_sub[iqcut]->Add(IMpippim_IMnpim_wK0orwSid_n_mix[iqcut],-1.0); 
    IMpippim_IMnpim_wK0orwSid_n_sub[iqcut]->SetTitle(Form("IMpippim_IMnpim_wK0orwSid_n_sub_%s",cqcut[iqcut]));
    if(RemoveNegative)IMpippim_IMnpim_wK0orwSid_n_sub[iqcut]->SetMinimum(0);

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

    std::cout << iqcut  << std::endl;
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
    for(int iq=0;iq<nqcut;iq++){
      IMnpipi_overlapdeco_K0[ig][iq] = new TH1D(Form("IMnpipi_overlapdeco_K0_g%d_%d",ig,iq),Form("IMnpipi_overlapdeco_K0_g%d_%d",ig,iq),100,1,2);
      IMnpipi_overlapdeco_Sp[ig][iq] = new TH1D(Form("IMnpipi_overlapdeco_Sp_g%d_%d",ig,iq),Form("IMnpipi_overlapdeco_Sp_g%d_%d",ig,iq),100,1,2);
      IMnpipi_overlapdeco_Sm[ig][iq] = new TH1D(Form("IMnpipi_overlapdeco_Sm_g%d_%d",ig,iq),Form("IMnpipi_overlapdeco_Sm_g%d_%d",ig,iq),100,1,2);
    }
  }

  for(unsigned int ibin=39;ibin<nbintemplate;ibin++){
    for(int iq=0;iq<nqcut;iq++){
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
      IMnpip_wK0_n_sub_bin[ibin][iq]->Draw("H");
      IMnpip_wK0_n_sub_bin_select[ibin][iq] = (TH1D*)IMnpip_wK0_n_sub_bin[ibin][iq]->Clone(Form("IMnpip_wK0_n_sub_bin%d_%d_select",ibin,iq));
      IMnpip_wK0_n_sub_bin_select[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmap_MIN,anacuts::Sigmap_MAX);
      IMnpip_wK0_n_sub_bin_select[ibin][iq]->SetLineColor(2);
      IMnpip_wK0_n_sub_bin_select[ibin][iq]->SetFillColor(2);
      IMnpip_wK0_n_sub_bin_select[ibin][iq]->Draw("Hsame");
      IMnpip_wK0_n_sub_bin_sidelo[ibin][iq] = (TH1D*)IMnpip_wK0_n_sub_bin[ibin][iq]->Clone(Form("IMnpip_wK0_n_sub_bin%d_%d_sidelo",ibin,iq));
      IMnpip_wK0_n_sub_bin_sidehi[ibin][iq] = (TH1D*)IMnpip_wK0_n_sub_bin[ibin][iq]->Clone(Form("IMnpip_wK0_n_sub_bin%d_%d_sidehi",ibin,iq));
      IMnpip_wK0_n_sub_bin_sidelo[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmap_MIN-4.0*anacuts::Sigmap_sigma, anacuts::Sigmap_MIN-2.0*anacuts::Sigmap_sigma);
      IMnpip_wK0_n_sub_bin_sidehi[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmap_MAX+2.0*anacuts::Sigmap_sigma, anacuts::Sigmap_MAX+4.0*anacuts::Sigmap_sigma);
      IMnpip_wK0_n_sub_bin_sidelo[ibin][iq]->SetLineColor(4);
      IMnpip_wK0_n_sub_bin_sidelo[ibin][iq]->SetFillColor(4);
      IMnpip_wK0_n_sub_bin_sidehi[ibin][iq]->SetLineColor(4);
      IMnpip_wK0_n_sub_bin_sidehi[ibin][iq]->SetFillColor(4);
      IMnpip_wK0_n_sub_bin_sidelo[ibin][iq]->Draw("Hsame");
      IMnpip_wK0_n_sub_bin_sidehi[ibin][iq]->Draw("Hsame");

      cIMpippim_IMnpip_n_sub_bin_cut[ibin][iq]->cd(4);
      IMpippim_wSid_n_Sp_sub_bin[ibin][iq]
      = (TH1D*)IMpippim_IMnpip_wSid_n_Sp_bin_sub[ibin][iq]->ProjectionY(Form("IMpippim_wSid_n_Sp_sub_bin%d_%d",ibin,iq));
      IMpippim_wSid_n_Sp_sub_bin[ibin][iq]->Draw("H");
      IMpippim_wSid_n_Sp_sub_bin_select[ibin][iq] = (TH1D*)IMpippim_wSid_n_Sp_sub_bin[ibin][iq]->Clone(Form("IMpippim_wSid_n_Sp_sub_bin_%d_%d_select",ibin,iq));
      IMpippim_wSid_n_Sp_sub_bin_select[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::pipi_MIN_narrow,anacuts::pipi_MAX_narrow);
      IMpippim_wSid_n_Sp_sub_bin_select[ibin][iq]->SetLineColor(2);
      IMpippim_wSid_n_Sp_sub_bin_select[ibin][iq]->SetFillColor(2);
      IMpippim_wSid_n_Sp_sub_bin_select[ibin][iq]->Draw("Hsame");
      IMpippim_wSid_n_Sp_sub_bin_sidelo[ibin][iq] = (TH1D*)IMpippim_wSid_n_Sp_sub_bin[ibin][iq]->Clone(Form("IMpippim_wSid_n_Sp_sub_bin_%d_%d_sidelo",ibin,iq));
      IMpippim_wSid_n_Sp_sub_bin_sidelo2[ibin][iq] = (TH1D*)IMpippim_wSid_n_Sp_sub_bin[ibin][iq]->Clone(Form("IMpippim_wSid_n_Sp_sub_bin_%d_%d_sidelo2",ibin,iq));
      IMpippim_wSid_n_Sp_sub_bin_sidehi[ibin][iq] = (TH1D*)IMpippim_wSid_n_Sp_sub_bin[ibin][iq]->Clone(Form("IMpippim_wSid_n_Sp_sub_bin_%d_%d_sidehi",ibin,iq));
      IMpippim_wSid_n_Sp_sub_bin_sidelo[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::pipi_MIN_narrow-4.0*anacuts::K0_sigma,anacuts::pipi_MIN_narrow-2.0*anacuts::K0_sigma);
      IMpippim_wSid_n_Sp_sub_bin_sidelo2[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::pipi_MIN_narrow-6.0*anacuts::K0_sigma,anacuts::pipi_MIN_narrow-2.0*anacuts::K0_sigma);
      IMpippim_wSid_n_Sp_sub_bin_sidehi[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::pipi_MAX_narrow+2.0*anacuts::K0_sigma,anacuts::pipi_MAX_narrow+4.0*anacuts::K0_sigma);
      IMpippim_wSid_n_Sp_sub_bin_sidelo[ibin][iq]->SetLineColor(4);
      IMpippim_wSid_n_Sp_sub_bin_sidelo[ibin][iq]->SetFillColor(4);
      IMpippim_wSid_n_Sp_sub_bin_sidelo2[ibin][iq]->SetLineColor(4);
      IMpippim_wSid_n_Sp_sub_bin_sidelo2[ibin][iq]->SetFillColor(4);
      IMpippim_wSid_n_Sp_sub_bin_sidehi[ibin][iq]->SetLineColor(4);
      IMpippim_wSid_n_Sp_sub_bin_sidehi[ibin][iq]->SetFillColor(4);
      if(ibin>47){
        IMpippim_wSid_n_Sp_sub_bin_sidelo[ibin][iq]->Draw("Hsame");
        IMpippim_wSid_n_Sp_sub_bin_sidehi[ibin][iq]->Draw("Hsame");
      }else{
        IMpippim_wSid_n_Sp_sub_bin_sidelo2[ibin][iq]->Draw("Hsame");
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
      double inteK0 = IMpippim_wSid_n_Sp_sub_bin_select[ibin][iq]->Integral();  //->Integral(binpipi_MIN,binpipi_MAX);
      double inteK0sidelo = IMpippim_wSid_n_Sp_sub_bin_sidelo[ibin][iq]->Integral();  //>Integral(binpipi_MIN_4sigma,binpipi_MIN-1);
      double inteK0sidelo2 = IMpippim_wSid_n_Sp_sub_bin_sidelo2[ibin][iq]->Integral();  //>Integral(binpipi_MIN_4sigma,binpipi_MIN-1);
      double inteK0sidehi = IMpippim_wSid_n_Sp_sub_bin_sidehi[ibin][iq]->Integral();  //>Integral(binpipi_MIN_4sigma,binpipi_MIN-1);
      double K0net = 0.0; 
      if(ibin>47) K0net = inteK0-inteK0sidelo-inteK0sidehi;
      else        K0net = inteK0-inteK0sidelo2;
      if(K0net<0.0) K0net = 0.0;
      double inteSp = IMnpip_wK0_n_sub_bin_select[ibin][iq]->Integral();  //->Integral(binnpip_MIN,binnpip_MAX);
      double inteSpsidelo = IMnpip_wK0_n_sub_bin_sidelo[ibin][iq]->Integral();  //(binnpip_MIN_2sigma,binnpip_MIN-1);
      double inteSpsidehi = IMnpip_wK0_n_sub_bin_sidehi[ibin][iq]->Integral();//(binnpip_MAX+1,binnpip_MAX_2sigma);
      double Spnet = inteSp - inteSpsidelo - inteSpsidehi;
      if(Spnet<0.0) Spnet = 0.0;
      
      pt->AddText(Form("IM(n#pi^{-}#pi^{+})  %0.2f-%0.2f",1.0+ibin*1.0/nbintemplate,1.0+(ibin+1.0)/nbintemplate)); 
      pt->AddText(Form("K0 count     %0.2f ",inteK0));
      if(ibin>47){
        pt->AddText(Form("K0 side low  %0.2f ",inteK0sidelo));
        pt->AddText(Form("K0 side high %0.2f ",inteK0sidehi));
      }else{
        pt->AddText(Form("K0 side low  %0.2f ",inteK0sidelo2));
      }
      pt->AddText(Form("K0 net (model)  %0.2f    ", K0net));
      pt->AddText(Form("Sigma+ count %0.2f ",inteSp));
      pt->AddText(Form("Sigma+ side low   %0.2f ",inteSpsidelo));
      pt->AddText(Form("Sigma+ side high  %0.2f ",inteSpsidehi));
      pt->AddText(Form("Sigma+ net (model)  %0.2f ", inteSp-inteSpsidelo-inteSpsidehi));
      pt->Draw();
      IMnpipi_overlapdeco_K0[0][iq]->Fill(1.0+ibin*1.0/nbintemplate,K0net);
      IMnpipi_overlapdeco_Sp[0][iq]->Fill(1.0+ibin*1.0/nbintemplate,Spnet);
      
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
      IMnpim_wK0_n_sub_bin[ibin][iq] ->SetTitle(Form("IMnpim_wK0_n_sub_bin%d_%d",ibin,iq));
      IMnpim_wK0_n_sub_bin[ibin][iq]->Draw("H");
      IMnpim_wK0_n_sub_bin_select[ibin][iq] = (TH1D*)IMnpim_wK0_n_sub_bin[ibin][iq]->Clone(Form("IMnpim_wK0_n_sub_bin%d_%d_select",ibin,iq));
      IMnpim_wK0_n_sub_bin_select[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmam_MIN,anacuts::Sigmam_MAX);
      IMnpim_wK0_n_sub_bin_select[ibin][iq]->SetLineColor(2);
      IMnpim_wK0_n_sub_bin_select[ibin][iq]->SetFillColor(2);
      IMnpim_wK0_n_sub_bin_select[ibin][iq]->Draw("Hsame");
      IMnpim_wK0_n_sub_bin_sidelo[ibin][iq] = (TH1D*)IMnpim_wK0_n_sub_bin[ibin][iq]->Clone(Form("IMnpim_wK0_n_sub_bin%d_%d_selectlo",ibin,iq));
      IMnpim_wK0_n_sub_bin_sidehi[ibin][iq] = (TH1D*)IMnpim_wK0_n_sub_bin[ibin][iq]->Clone(Form("IMnpim_wK0_n_sub_bin%d_%d_selecthi",ibin,iq));
      IMnpim_wK0_n_sub_bin_sidelo[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmam_MIN-4.0*anacuts::Sigmam_sigma,anacuts::Sigmam_MIN-2.0*anacuts::Sigmam_sigma);
      IMnpim_wK0_n_sub_bin_sidehi[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmam_MAX+2.0*anacuts::Sigmam_sigma,anacuts::Sigmam_MAX+4.0*anacuts::Sigmam_sigma);
      IMnpim_wK0_n_sub_bin_sidelo[ibin][iq]->SetLineColor(4);
      IMnpim_wK0_n_sub_bin_sidelo[ibin][iq]->SetFillColor(4);
      IMnpim_wK0_n_sub_bin_sidehi[ibin][iq]->SetLineColor(4);
      IMnpim_wK0_n_sub_bin_sidehi[ibin][iq]->SetFillColor(4);
      IMnpim_wK0_n_sub_bin_sidelo[ibin][iq]->Draw("Hsame");
      IMnpim_wK0_n_sub_bin_sidehi[ibin][iq]->Draw("Hsame");
      
      
      cIMpippim_IMnpim_n_sub_bin_cut[ibin][iq]->cd(4);
      IMpippim_wSid_n_Sm_sub_bin[ibin][iq]
      = (TH1D*)IMpippim_IMnpim_wSid_n_Sm_bin_sub[ibin][iq]->ProjectionY(Form("IMpippim_wSid_n_Sm_sub_bin%d_%d",ibin,iq));
      
      IMpippim_wSid_n_Sm_sub_bin[ibin][iq]->Draw("H");
      IMpippim_wSid_n_Sm_sub_bin_select[ibin][iq] = (TH1D*)IMpippim_wSid_n_Sm_sub_bin[ibin][iq]->Clone(Form("IMpippim_wSid_n_Sm_sub_bin_%d_%d_select",ibin,iq));
      IMpippim_wSid_n_Sm_sub_bin_select[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::pipi_MIN_narrow,anacuts::pipi_MAX_narrow);
      IMpippim_wSid_n_Sm_sub_bin_select[ibin][iq]->SetLineColor(2);
      IMpippim_wSid_n_Sm_sub_bin_select[ibin][iq]->SetFillColor(2);
      IMpippim_wSid_n_Sm_sub_bin_select[ibin][iq]->Draw("Hsame");
      IMpippim_wSid_n_Sm_sub_bin_sidelo[ibin][iq] = (TH1D*)IMpippim_wSid_n_Sm_sub_bin[ibin][iq]->Clone(Form("IMpippim_wSid_n_Sm_sub_bin_%d_%d_sidelo",ibin,iq));
      IMpippim_wSid_n_Sm_sub_bin_sidelo2[ibin][iq] = (TH1D*)IMpippim_wSid_n_Sm_sub_bin[ibin][iq]->Clone(Form("IMpippim_wSid_n_Sm_sub_bin_%d_%d_sidelo2",ibin,iq));
      IMpippim_wSid_n_Sm_sub_bin_sidehi[ibin][iq] = (TH1D*)IMpippim_wSid_n_Sm_sub_bin[ibin][iq]->Clone(Form("IMpippim_wSid_n_Sm_sub_bin_%d_%d_sidehi",ibin,iq));
      IMpippim_wSid_n_Sm_sub_bin_sidelo[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::pipi_MIN_narrow-4.0*anacuts::K0_sigma,anacuts::pipi_MIN_narrow-2.0*anacuts::K0_sigma);
      IMpippim_wSid_n_Sm_sub_bin_sidelo2[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::pipi_MIN_narrow-6.0*anacuts::K0_sigma,anacuts::pipi_MIN_narrow-2.0*anacuts::K0_sigma);
      IMpippim_wSid_n_Sm_sub_bin_sidehi[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::pipi_MAX_narrow+2.0*anacuts::K0_sigma,anacuts::pipi_MAX_narrow+4.0*anacuts::K0_sigma);
      IMpippim_wSid_n_Sm_sub_bin_sidelo[ibin][iq]->SetLineColor(4);
      IMpippim_wSid_n_Sm_sub_bin_sidelo[ibin][iq]->SetFillColor(4);
      IMpippim_wSid_n_Sm_sub_bin_sidelo2[ibin][iq]->SetLineColor(3);
      IMpippim_wSid_n_Sm_sub_bin_sidelo2[ibin][iq]->SetFillColor(3);
      IMpippim_wSid_n_Sm_sub_bin_sidehi[ibin][iq]->SetLineColor(4);
      IMpippim_wSid_n_Sm_sub_bin_sidehi[ibin][iq]->SetFillColor(4);
      if(ibin>47){
        IMpippim_wSid_n_Sm_sub_bin_sidelo[ibin][iq]->Draw("Hsame");
        IMpippim_wSid_n_Sm_sub_bin_sidehi[ibin][iq]->Draw("Hsame");
      }else{
        IMpippim_wSid_n_Sm_sub_bin_sidelo2[ibin][iq]->Draw("Hsame");
      }

      //gr_IMpippim_wSid_n_Sm_sub_bin[ibin][iq] = new TGraphErrors();
      //HistToRorateGraph(IMpippim_wSid_n_Sm_sub_bin[ibin][iq],*gr_IMpippim_wSid_n_Sm_sub_bin[ibin][iq]);
      //gr_IMpippim_wSid_n_Sm_sub_bin[ibin][iq]->Draw("AP");
      
      cIMpippim_IMnpim_n_sub_bin_cut[ibin][iq]->cd(2);
      TPaveText *pt2 = new TPaveText(.05,.05,.95,.7);
      double inteK0_g2 = IMpippim_wSid_n_Sm_sub_bin_select[ibin][iq]->Integral();
      double inteK0sidelo_g2 = IMpippim_wSid_n_Sm_sub_bin_sidelo[ibin][iq]->Integral();
      double inteK0sidelo2_g2 = IMpippim_wSid_n_Sm_sub_bin_sidelo2[ibin][iq]->Integral();
      double inteK0sidehi_g2 = IMpippim_wSid_n_Sm_sub_bin_sidehi[ibin][iq]->Integral();
      double K0net_g2 = 0.0;
      if(ibin>47) K0net_g2 = inteK0_g2-inteK0sidelo_g2-inteK0sidehi_g2;
      else        K0net_g2 = inteK0_g2-inteK0sidelo2_g2;
      if(K0net_g2<0.0) K0net_g2 = 0.0; 
      double inteSm_g2 = IMnpim_wK0_n_sub_bin_select[ibin][iq]->Integral();
      double inteSmsidelo_g2 = IMnpim_wK0_n_sub_bin_sidelo[ibin][iq]->Integral();
      double inteSmsidehi_g2 = IMnpim_wK0_n_sub_bin_sidehi[ibin][iq]->Integral();
      double Smnet_g2 = inteSm_g2- inteSmsidelo_g2 - inteSmsidehi_g2;
      if(Smnet_g2<0.0) Smnet_g2 = 0.0;

      pt2->AddText(Form("IM(n#pi^{-}#pi^{+})  %0.2f-%0.2f",1.0+ibin*1.0/nbintemplate,1.0+(ibin+1.0)/nbintemplate)); 
      pt2->AddText(Form("K0 count %0.2f ",inteK0_g2));
      if(ibin>47){
        pt2->AddText(Form("K0 side low  %0.2f ",inteK0sidelo_g2));
        pt2->AddText(Form("K0 side high %0.2f ",inteK0sidehi_g2));
      }else{
        pt2->AddText(Form("K0 side low  %0.2f ",inteK0sidelo2_g2));
      }
      pt2->AddText(Form("K0 net  (model)  %0.2f ",K0net_g2));
      pt2->AddText(Form("Sigma- count %0.2f ",inteSm_g2));
      pt2->AddText(Form("Sigma- side low   %0.2f ",inteSmsidelo_g2));
      pt2->AddText(Form("Sigma- side high  %0.2f ",inteSmsidehi_g2));
      pt2->AddText(Form("Sigma- net (model)  %0.2f ", Smnet_g2));
      pt2->Draw();
      IMnpipi_overlapdeco_K0[1][iq]->Fill(1.0+ibin*1.0/nbintemplate,K0net_g2);
      IMnpipi_overlapdeco_Sm[1][iq]->Fill(1.0+ibin*1.0/nbintemplate,Smnet_g2);
      
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
      IMnpip_wSid_n_Sm_sub_bin[ibin][iq]->Draw("H");
      IMnpip_wSid_n_Sm_sub_bin_select[ibin][iq] = (TH1D*)IMnpip_wSid_n_Sm_sub_bin[ibin][iq]->Clone(Form("IMnpip_wSid_n_Sm_sub_bin%d_%d_select",ibin,iq));
      IMnpip_wSid_n_Sm_sub_bin_select[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmap_MIN,anacuts::Sigmap_MAX);
      IMnpip_wSid_n_Sm_sub_bin_select[ibin][iq]->SetLineColor(2);
      IMnpip_wSid_n_Sm_sub_bin_select[ibin][iq]->SetFillColor(2);
      IMnpip_wSid_n_Sm_sub_bin_select[ibin][iq]->Draw("Hsame");
      IMnpip_wSid_n_Sm_sub_bin_sidelo[ibin][iq] = (TH1D*)IMnpip_wSid_n_Sm_sub_bin[ibin][iq]->Clone(Form("IMnpip_wSid_n_Sm_sub_bin%d_%d_sidelo",ibin,iq));
      IMnpip_wSid_n_Sm_sub_bin_sidelo2[ibin][iq] = (TH1D*)IMnpip_wSid_n_Sm_sub_bin[ibin][iq]->Clone(Form("IMnpip_wSid_n_Sm_sub_bin%d_%d_sidelo2",ibin,iq));
      IMnpip_wSid_n_Sm_sub_bin_sidehi[ibin][iq] = (TH1D*)IMnpip_wSid_n_Sm_sub_bin[ibin][iq]->Clone(Form("IMnpip_wSid_n_Sm_sub_bin%d_%d_sidehi",ibin,iq));
      IMnpip_wSid_n_Sm_sub_bin_sidehi2[ibin][iq] = (TH1D*)IMnpip_wSid_n_Sm_sub_bin[ibin][iq]->Clone(Form("IMnpip_wSid_n_Sm_sub_bin%d_%d_sidehi2",ibin,iq));
      IMnpip_wSid_n_Sm_sub_bin_sidelo[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmap_MIN-4.0*anacuts::Sigmap_sigma,anacuts::Sigmap_MIN-2.0*anacuts::Sigmap_sigma);
      IMnpip_wSid_n_Sm_sub_bin_sidelo2[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmap_MIN-6.0*anacuts::Sigmap_sigma,anacuts::Sigmap_MIN-2.0*anacuts::Sigmap_sigma);
      IMnpip_wSid_n_Sm_sub_bin_sidehi[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmap_MAX+2.0*anacuts::Sigmap_sigma,anacuts::Sigmap_MAX+4.0*anacuts::Sigmap_sigma);
      IMnpip_wSid_n_Sm_sub_bin_sidehi2[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmap_MAX+2.0*anacuts::Sigmap_sigma,anacuts::Sigmap_MAX+6.0*anacuts::Sigmap_sigma);
      IMnpip_wSid_n_Sm_sub_bin_sidelo[ibin][iq]->SetLineColor(4);
      IMnpip_wSid_n_Sm_sub_bin_sidelo[ibin][iq]->SetFillColor(4);
      IMnpip_wSid_n_Sm_sub_bin_sidelo2[ibin][iq]->SetLineColor(4);
      IMnpip_wSid_n_Sm_sub_bin_sidelo2[ibin][iq]->SetFillColor(4);
      IMnpip_wSid_n_Sm_sub_bin_sidehi[ibin][iq]->SetLineColor(4);
      IMnpip_wSid_n_Sm_sub_bin_sidehi[ibin][iq]->SetFillColor(4);
      IMnpip_wSid_n_Sm_sub_bin_sidehi2[ibin][iq]->SetLineColor(4);
      IMnpip_wSid_n_Sm_sub_bin_sidehi2[ibin][iq]->SetFillColor(4);
      if(ibin>42){
        IMnpip_wSid_n_Sm_sub_bin_sidelo[ibin][iq]->Draw("Hsame");
        IMnpip_wSid_n_Sm_sub_bin_sidehi[ibin][iq]->Draw("Hsame");
      }else{
        IMnpip_wSid_n_Sm_sub_bin_sidelo2[ibin][iq]->Draw("Hsame");
      }

      cIMnpim_IMnpip_n_sub_bin_cut[ibin][iq]->cd(4);
      IMnpim_wSid_n_Sp_sub_bin[ibin][iq]
      = (TH1D*)IMnpim_IMnpip_wSid_n_Sp_bin_sub[ibin][iq]->ProjectionY(Form("IMnpim_wSid_n_Sp_sub_bin%d_%d",ibin,iq));
      IMnpim_wSid_n_Sp_sub_bin[ibin][iq]->Draw("H");
      IMnpim_wSid_n_Sp_sub_bin_select[ibin][iq] = (TH1D*)IMnpim_wSid_n_Sp_sub_bin[ibin][iq]->Clone(Form("IMnpim_wSid_n_Sp_sub_bin%d_%d_select",ibin,iq));
      IMnpim_wSid_n_Sp_sub_bin_select[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmam_MIN,anacuts::Sigmam_MAX);
      IMnpim_wSid_n_Sp_sub_bin_select[ibin][iq]->SetLineColor(2);
      IMnpim_wSid_n_Sp_sub_bin_select[ibin][iq]->SetFillColor(2);
      IMnpim_wSid_n_Sp_sub_bin_select[ibin][iq]->Draw("Hsame");
      IMnpim_wSid_n_Sp_sub_bin_sidelo[ibin][iq] = (TH1D*)IMnpim_wSid_n_Sp_sub_bin[ibin][iq]->Clone(Form("IMnpim_wSid_n_Sp_sub_bin%d_%d_sidelo",ibin,iq));
      IMnpim_wSid_n_Sp_sub_bin_sidelo2[ibin][iq] = (TH1D*)IMnpim_wSid_n_Sp_sub_bin[ibin][iq]->Clone(Form("IMnpim_wSid_n_Sp_sub_bin%d_%d_sidelo2",ibin,iq));
      IMnpim_wSid_n_Sp_sub_bin_sidehi[ibin][iq] = (TH1D*)IMnpim_wSid_n_Sp_sub_bin[ibin][iq]->Clone(Form("IMnpim_wSid_n_Sp_sub_bin%d_%d_sidehi",ibin,iq));
      IMnpim_wSid_n_Sp_sub_bin_sidehi2[ibin][iq] = (TH1D*)IMnpim_wSid_n_Sp_sub_bin[ibin][iq]->Clone(Form("IMnpim_wSid_n_Sp_sub_bin%d_%d_sidehi2",ibin,iq));
      IMnpim_wSid_n_Sp_sub_bin_sidelo[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmam_MIN-4.0*anacuts::Sigmam_sigma,anacuts::Sigmam_MIN-2.0*anacuts::Sigmam_sigma);
      IMnpim_wSid_n_Sp_sub_bin_sidelo2[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmam_MIN-6.0*anacuts::Sigmam_sigma,anacuts::Sigmam_MIN-2.0*anacuts::Sigmam_sigma);
      IMnpim_wSid_n_Sp_sub_bin_sidehi[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmam_MAX+2.0*anacuts::Sigmam_sigma,anacuts::Sigmam_MAX+4.0*anacuts::Sigmam_sigma);
      IMnpim_wSid_n_Sp_sub_bin_sidehi2[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmam_MAX+2.0*anacuts::Sigmam_sigma,anacuts::Sigmam_MAX+6.0*anacuts::Sigmam_sigma);
      IMnpim_wSid_n_Sp_sub_bin_sidelo[ibin][iq]->SetLineColor(4);
      IMnpim_wSid_n_Sp_sub_bin_sidelo[ibin][iq]->SetFillColor(4);
      IMnpim_wSid_n_Sp_sub_bin_sidelo2[ibin][iq]->SetLineColor(4);
      IMnpim_wSid_n_Sp_sub_bin_sidelo2[ibin][iq]->SetFillColor(4);
      IMnpim_wSid_n_Sp_sub_bin_sidehi[ibin][iq]->SetLineColor(4);
      IMnpim_wSid_n_Sp_sub_bin_sidehi[ibin][iq]->SetFillColor(4);
      IMnpim_wSid_n_Sp_sub_bin_sidehi2[ibin][iq]->SetLineColor(4);
      IMnpim_wSid_n_Sp_sub_bin_sidehi2[ibin][iq]->SetFillColor(4);
      if(ibin>42){
        IMnpim_wSid_n_Sp_sub_bin_sidelo[ibin][iq]->Draw("Hsame");
        IMnpim_wSid_n_Sp_sub_bin_sidehi[ibin][iq]->Draw("Hsame");
      }else{
        IMnpim_wSid_n_Sp_sub_bin_sidelo2[ibin][iq]->Draw("Hsame");
      }

      cIMnpim_IMnpip_n_sub_bin_cut[ibin][iq]->cd(2);
      TPaveText *pt3 = new TPaveText(.05,.05,.95,.7);
      double inteSm_g3 = IMnpim_wSid_n_Sp_sub_bin_select[ibin][iq]->Integral();
      double inteSmsidelo_g3 = IMnpim_wSid_n_Sp_sub_bin_sidelo[ibin][iq]->Integral();
      double inteSmsidelo2_g3 = IMnpim_wSid_n_Sp_sub_bin_sidelo2[ibin][iq]->Integral();
      double inteSmsidehi_g3 = IMnpim_wSid_n_Sp_sub_bin_sidehi[ibin][iq]->Integral();
      double inteSmsidehi2_g3 = IMnpim_wSid_n_Sp_sub_bin_sidehi2[ibin][iq]->Integral();
      double Smnet_g3 = 0;
      if(ibin>42)Smnet_g3 = inteSm_g3-inteSmsidelo_g3-inteSmsidehi_g3;
      else Smnet_g3 = inteSm_g3 - inteSmsidelo2_g3; 
      if(Smnet_g3<0.0) Smnet_g3 = 0.0;

      double inteSp_g3 = IMnpip_wSid_n_Sm_sub_bin_select[ibin][iq]->Integral();
      double inteSpsidelo_g3 = IMnpip_wSid_n_Sm_sub_bin_sidelo[ibin][iq]->Integral();
      double inteSpsidelo2_g3 = IMnpip_wSid_n_Sm_sub_bin_sidelo2[ibin][iq]->Integral();
      double inteSpsidehi_g3 = IMnpip_wSid_n_Sm_sub_bin_sidehi[ibin][iq]->Integral();
      double inteSpsidehi2_g3 = IMnpip_wSid_n_Sm_sub_bin_sidehi2[ibin][iq]->Integral();
      double Spnet_g3 = 0;
      if(ibin>42) Spnet_g3 = inteSp_g3- inteSpsidelo_g3 - inteSpsidehi_g3;
      else Spnet_g3 = inteSp_g3 - inteSpsidelo2_g3;
      if(Spnet_g3<0.0) Spnet_g3 = 0.0;

      pt3->AddText(Form("IM(n#pi^{-}#pi^{+})  %0.2f-%0.2f",1.0+ibin*1.0/nbintemplate,1.0+(ibin+1.0)/nbintemplate)); 
      pt3->AddText(Form("Sigma- count %0.2f ",inteSm_g3));
      if(ibin>42){
        pt3->AddText(Form("Sigma- side low   %0.2f ",inteSmsidelo_g3));
        pt3->AddText(Form("Sigma- side high  %0.2f ",inteSmsidehi_g3));
      }else{
        pt3->AddText(Form("Sigma- side low2   %0.2f ",inteSmsidelo2_g3));
      } 
      pt3->AddText(Form("Sigma- net (model) %0.2f ", Smnet_g3));
      pt3->AddText(Form("Sigma+ count %0.2f ",inteSp_g3));
      if(ibin>42){
        pt3->AddText(Form("Sigma+ side low   %0.2f ",inteSpsidelo_g3));
        pt3->AddText(Form("Sigma+ side high   %0.2f ",inteSpsidehi_g3));
      }else{
        pt3->AddText(Form("Sigma+ side low2   %0.2f ",inteSpsidelo2_g3));
      }
      pt3->AddText(Form("Sigma+ net (model) %0.2f ",Spnet_g3)); 
      pt3->Draw();
      IMnpipi_overlapdeco_Sm[2][iq]->Fill(1.0+ibin*1.0/nbintemplate,Smnet_g3);
      IMnpipi_overlapdeco_Sp[2][iq]->Fill(1.0+ibin*1.0/nbintemplate,Spnet_g3);
    }
  }
   
  TCanvas *csum_Sp[ngroup][nqcut];
  TCanvas *csum_Sm[ngroup][nqcut];
  TCanvas *csum_K0[ngroup][nqcut];
  for(int ig=0;ig<ngroup;ig++){
    for(int iq=0;iq<nqcut;iq++){
    csum_K0[ig][iq] = new TCanvas(Form("csum_K0_g%d_q%d",ig,iq),Form("csum_K0_g%d_q%d",ig,iq));
    IMnpipi_overlapdeco_K0[ig][iq]->Draw("H");
    csum_Sp[ig][iq] = new TCanvas(Form("csum_Sp_g%d_q%d",ig,iq),Form("csum_Sp_g%d_q%d",ig,iq));
    IMnpipi_overlapdeco_Sp[ig][iq]->Draw("H");
    csum_Sm[ig][iq] = new TCanvas(Form("csum_Sm_g%d_q%d",ig,iq),Form("csum_Sm_g%d_q%d",ig,iq));
    IMnpipi_overlapdeco_Sm[ig][iq]->Draw("H");
    }
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
