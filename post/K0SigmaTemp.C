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
  TH2F* IMpippim_IMnpip_wSid_n_Sp_data[3];
  TH2F* IMpippim_IMnpip_wSid_n_Sp_bin_data[nbintemplate][3];
  TH2F* IMpippim_IMnpip_wSid_n_Sp_mix[3];
  TH2F* IMpippim_IMnpip_wSid_n_Sp_bin_mix[nbintemplate][3];
  TH2F* IMpippim_IMnpip_wSid_n_Sp_sub[3];
  TH2F* IMpippim_IMnpip_wSid_n_Sp_bin_sub[nbintemplate][3];
  TH2F* IMpippim_IMnpim_wSid_n_Sm_data[3];
  TH2F* IMpippim_IMnpim_wSid_n_Sm_bin_data[nbintemplate][3];
  TH2F* IMpippim_IMnpim_wSid_n_Sm_mix[3];
  TH2F* IMpippim_IMnpim_wSid_n_Sm_bin_mix[nbintemplate][3];
  TH2F* IMpippim_IMnpim_wSid_n_Sm_sub[3];
  TH2F* IMpippim_IMnpim_wSid_n_Sm_bin_sub[nbintemplate][3];
  TH2F* IMpippim_IMnpip_wK0_n_data[3];
  TH2F* IMpippim_IMnpip_wK0_n_bin_data[nbintemplate][3];
  TH2F* IMpippim_IMnpip_wK0_n_mix[3];
  TH2F* IMpippim_IMnpip_wK0_n_bin_mix[nbintemplate][3];
  TH2F* IMpippim_IMnpip_wK0_n_sub[3];
  TH2F* IMpippim_IMnpip_wK0_n_bin_sub[nbintemplate][3];
  TH2F* IMpippim_IMnpim_wK0_n_data[3];
  TH2F* IMpippim_IMnpim_wK0_n_bin_data[nbintemplate][3];
  TH2F* IMpippim_IMnpim_wK0_n_mix[3];
  TH2F* IMpippim_IMnpim_wK0_n_bin_mix[nbintemplate][3];
  TH2F* IMpippim_IMnpim_wK0_n_sub[3];
  TH2F* IMpippim_IMnpim_wK0_n_bin_sub[nbintemplate][3];
  TH2F* IMpippim_IMnpip_wK0orwSid_n_data[3];
  TH2F* IMpippim_IMnpim_wK0orwSid_n_data[3];
  TH2F* IMpippim_IMnpip_wK0orwSid_n_bin_data[nbintemplate][3];
  TH2F* IMpippim_IMnpim_wK0orwSid_n_bin_data[nbintemplate][3];
  TH2F* IMpippim_IMnpip_wK0orwSid_n_mix[3];
  TH2F* IMpippim_IMnpim_wK0orwSid_n_mix[3];
  TH2F* IMpippim_IMnpip_wK0orwSid_n_bin_mix[nbintemplate][3];
  TH2F* IMpippim_IMnpim_wK0orwSid_n_bin_mix[nbintemplate][3];
  TH2F* IMpippim_IMnpip_wK0orwSid_n_sub[3];
  TH2F* IMpippim_IMnpim_wK0orwSid_n_sub[3];
  TH2F* IMpippim_IMnpip_wK0orwSid_n_bin_sub[nbintemplate][3];
  TH2F* IMpippim_IMnpim_wK0orwSid_n_bin_sub[nbintemplate][3];
  TH2F* IMnpim_IMnpip_n_data[3];
  TH2F* IMnpim_IMnpip_n_bin_data[nbintemplate][3];
  TH2F* IMnpim_IMnpip_n_mix[3];
  TH2F* IMnpim_IMnpip_n_bin_mix[nbintemplate][3];
  TH2F* IMnpim_IMnpip_n_sub[3];
  TH2F* IMnpim_IMnpip_n_bin_sub[nbintemplate][3];
  TH2F* IMnpim_IMnpip_wSid_n_data[3];
  TH2F* IMnpim_IMnpip_wSid_n_bin_data[nbintemplate][3];
  TH2F* IMnpim_IMnpip_wSid_n_mix[3];
  TH2F* IMnpim_IMnpip_wSid_n_bin_mix[nbintemplate][3];
  TH2F* IMnpim_IMnpip_wSid_n_sub[3];
  TH2F* IMnpim_IMnpip_wSid_n_bin_sub[nbintemplate][3];
  TH2F* IMnpim_IMnpip_wSid_n_Sp_data[3];
  TH2F* IMnpim_IMnpip_wSid_n_Sp_bin_data[nbintemplate][3];
  TH2F* IMnpim_IMnpip_wSid_n_Sp_mix[3];
  TH2F* IMnpim_IMnpip_wSid_n_Sp_bin_mix[nbintemplate][3];
  TH2F* IMnpim_IMnpip_wSid_n_Sp_sub[3];
  TH2F* IMnpim_IMnpip_wSid_n_Sp_bin_sub[nbintemplate][3];
  TH2F* IMnpim_IMnpip_wSid_n_Sm_data[3];
  TH2F* IMnpim_IMnpip_wSid_n_Sm_bin_data[nbintemplate][3];
  TH2F* IMnpim_IMnpip_wSid_n_Sm_mix[3];
  TH2F* IMnpim_IMnpip_wSid_n_Sm_bin_mix[nbintemplate][3];
  TH2F* IMnpim_IMnpip_wSid_n_Sm_sub[3];
  TH2F* IMnpim_IMnpip_wSid_n_Sm_bin_sub[nbintemplate][3];
  

  const char cqcut[][6]= {"all","qlo","qhi"};
  const int nqcut=1;
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
      if(RemoveNegative)IMpippim_IMnpip_wSid_n_Sp_bin_sub[ibin][iqcut]->SetMinimum(0);

      IMpippim_IMnpim_wSid_n_Sm_bin_data[ibin][iqcut] = (TH2F*)fr[iqcut]->Get(Form("IMpippim_IMnpim_wSid_n_Sm_bin%d",ibin));
      IMpippim_IMnpim_wSid_n_Sm_bin_mix[ibin][iqcut] = (TH2F*)fmix[iqcut]->Get(Form("IMpippim_IMnpim_wSid_n_Sm_bin%d",ibin));
      IMpippim_IMnpim_wSid_n_Sm_bin_sub[ibin][iqcut]
      = (TH2F*)IMpippim_IMnpim_wSid_n_Sm_bin_data[ibin][iqcut]->Clone(Form("IMpippim_IMnpim_wSid_n_Sm_bin%d_sub%s",ibin,cqcut[iqcut]));
      IMpippim_IMnpim_wSid_n_Sm_bin_sub[ibin][iqcut]->Add(IMpippim_IMnpim_wSid_n_Sm_bin_mix[ibin][iqcut],-1.0);
      IMpippim_IMnpim_wSid_n_Sm_bin_sub[ibin][iqcut]->SetTitle(Form("IMpippim_IMnpim_wSid_n_Sm_bin%d_sub%s",ibin,cqcut[iqcut]));
      if(RemoveNegative)IMpippim_IMnpim_wSid_n_Sm_bin_sub[ibin][iqcut]->SetMinimum(0);

      IMpippim_IMnpip_wK0_n_bin_data[ibin][iqcut] = (TH2F*)fr[iqcut]->Get(Form("IMpippim_IMnpip_wK0_n_bin%d",ibin));
      IMpippim_IMnpip_wK0_n_bin_mix[ibin][iqcut] = (TH2F*)fmix[iqcut]->Get(Form("IMpippim_IMnpip_wK0_n_bin%d",ibin));
      IMpippim_IMnpip_wK0_n_bin_sub[ibin][iqcut] = (TH2F*)IMpippim_IMnpip_wK0_n_bin_data[ibin][iqcut]->Clone(Form("IMpippim_IMnpip_wK0_n_bin%d_sub%s",ibin,cqcut[iqcut]));
      IMpippim_IMnpip_wK0_n_bin_sub[ibin][iqcut]->Add(IMpippim_IMnpip_wK0_n_bin_mix[ibin][iqcut],-1.0);
      IMpippim_IMnpip_wK0_n_bin_sub[ibin][iqcut]->SetTitle(Form("IMpippim_IMnpip_wK0_n_bin%d_sub%s",ibin,cqcut[iqcut]));
      if(RemoveNegative)IMpippim_IMnpip_wK0_n_bin_sub[ibin][iqcut]->SetMinimum(0);

      IMpippim_IMnpim_wK0_n_bin_data[ibin][iqcut] = (TH2F*)fr[iqcut]->Get(Form("IMpippim_IMnpim_wK0_n_bin%d",ibin));
      IMpippim_IMnpim_wK0_n_bin_mix[ibin][iqcut] = (TH2F*)fmix[iqcut]->Get(Form("IMpippim_IMnpim_wK0_n_bin%d",ibin));
      IMpippim_IMnpim_wK0_n_bin_sub[ibin][iqcut] = (TH2F*)IMpippim_IMnpim_wK0_n_bin_data[ibin][iqcut]->Clone(Form("IMpippim_IMnpim_wK0_n_bin%d_sub%s",ibin,cqcut[iqcut]));
      IMpippim_IMnpim_wK0_n_bin_sub[ibin][iqcut]->Add(IMpippim_IMnpim_wK0_n_bin_mix[ibin][iqcut],-1.0);
      IMpippim_IMnpim_wK0_n_bin_sub[ibin][iqcut]->SetTitle(Form("IMpippim_IMnpim_wK0_n_bin%d_sub%s",ibin,cqcut[iqcut]));
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
      
      IMnpim_IMnpip_wSid_n_Sp_bin_data[ibin][iqcut] = (TH2F*)fr[iqcut]->Get(Form("IMnpim_IMnpip_wSid_n_Sp_bin%d",ibin));
      IMnpim_IMnpip_wSid_n_Sp_bin_mix[ibin][iqcut] = (TH2F*)fmix[iqcut]->Get(Form("IMnpim_IMnpip_wSid_n_Sp_bin%d",ibin));
      IMnpim_IMnpip_wSid_n_Sp_bin_sub[ibin][iqcut] 
      = (TH2F*)IMnpim_IMnpip_wSid_n_Sp_bin_data[ibin][iqcut]->Clone(Form("IMnpim_IMnpip_wSid_n_Sp_bin%d_%s",ibin,cqcut[iqcut]));
      IMnpim_IMnpip_wSid_n_Sp_bin_sub[ibin][iqcut]->Add(IMnpim_IMnpip_wSid_n_Sp_bin_mix[ibin][iqcut],-1);
      IMnpim_IMnpip_wSid_n_Sp_bin_sub[ibin][iqcut]->SetTitle(Form("IMnpim_IMnpip_wSid_n_Sp_bin%d_%s",ibin,cqcut[iqcut]));
      if(RemoveNegative)IMnpim_IMnpip_wSid_n_Sp_bin_sub[ibin][iqcut]->SetMinimum(0);

      IMnpim_IMnpip_wSid_n_Sm_bin_data[ibin][iqcut] = (TH2F*)fr[iqcut]->Get(Form("IMnpim_IMnpip_wSid_n_Sm_bin%d",ibin));
      IMnpim_IMnpip_wSid_n_Sm_bin_mix[ibin][iqcut] = (TH2F*)fmix[iqcut]->Get(Form("IMnpim_IMnpip_wSid_n_Sm_bin%d",ibin));
      IMnpim_IMnpip_wSid_n_Sm_bin_sub[ibin][iqcut]
      = (TH2F*)IMnpim_IMnpip_wSid_n_Sm_bin_data[ibin][iqcut]->Clone(Form("IMnpim_IMnpip_wSid_n_Sm_bin%d_%s",ibin,cqcut[iqcut]));
      IMnpim_IMnpip_wSid_n_Sm_bin_sub[ibin][iqcut]->Add(IMnpim_IMnpip_wSid_n_Sm_bin_mix[ibin][iqcut],-1);
      IMnpim_IMnpip_wSid_n_Sm_bin_sub[ibin][iqcut]->SetTitle(Form("IMnpim_IMnpip_wSid_n_Sm_bin%d_%s",ibin,cqcut[iqcut]));
      if(RemoveNegative)IMnpim_IMnpip_wSid_n_Sm_bin_sub[ibin][iqcut]->SetMinimum(0);
    }
  }

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
  //TGraphErrors *gr_IMpippim_wSid_n_Sp_sub_bin[nbintemplate][nqcut];
  //TGraphErrors *gr_IMpippim_wSid_n_Sp_sub_bin[nbintemplate][nqcut];
  
  for(unsigned int ibin=40;ibin<nbintemplate;ibin++){
    for(int iq=0;iq<nqcut;iq++){
      //first canvas to display IM(pi+pi-) vs IM(npi+) and their simple projection just to see the signal and other background 
      //no complicated cut at this moment
      cIMpippim_IMnpip_n_sub_bin[ibin][iq] = new TCanvas(Form("cIMpippim_IMnpip_n_sub_bin%d_%d",ibin,iq),Form("cIMpippim_IMnpip_n_sub_bin%d_%d",ibin,iq),800,800);
      cIMpippim_IMnpip_n_sub_bin[ibin][iq]->Divide(2,2,0,0);
      cIMpippim_IMnpip_n_sub_bin[ibin][iq]->cd(3);
      IMpippim_IMnpip_n_bin_sub[ibin][iq]->RebinX(4);
      IMpippim_IMnpip_n_bin_sub[ibin][iq]->RebinY(5);
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
      
      //2nd canvas to display the calculation of decomposion of K0, Sigma+,Sigma- 
      //2D plot in canvas(3) has cut to select K0 or Sigma+
      //projection plot in canvas(1) select Sigma+ events candidate + (Sigma+ & K0 ) overlap events
      //projection plot in canvas(4) select K0     events candidate + (Sigma+ & K0 ) overlap events 
      cIMpippim_IMnpip_n_sub_bin_cut[ibin][iq] 
      = new TCanvas(Form("cIMpippim_IMnpip_n_sub_bin_cut_%d_%d",ibin,iq),Form("cIMpippim_IMnpip_n_sub_bin_cut%d_%d",ibin,iq),800,800);
      cIMpippim_IMnpip_n_sub_bin_cut[ibin][iq]->Divide(2,2,0,0);
      cIMpippim_IMnpip_n_sub_bin_cut[ibin][iq]->cd(3);
      
      IMpippim_IMnpip_wK0orwSid_n_bin_sub[ibin][iq]->RebinX(4);
      IMpippim_IMnpip_wK0orwSid_n_bin_sub[ibin][iq]->RebinY(5);
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
      IMpippim_wSid_n_Sp_sub_bin_sidelo[ibin][iq]->SetLineColor(3);
      IMpippim_wSid_n_Sp_sub_bin_sidelo[ibin][iq]->SetFillColor(3);
      IMpippim_wSid_n_Sp_sub_bin_sidelo2[ibin][iq]->SetLineColor(4);
      IMpippim_wSid_n_Sp_sub_bin_sidelo2[ibin][iq]->SetFillColor(4);
      IMpippim_wSid_n_Sp_sub_bin_sidehi[ibin][iq]->SetLineColor(3);
      IMpippim_wSid_n_Sp_sub_bin_sidehi[ibin][iq]->SetFillColor(3);
      IMpippim_wSid_n_Sp_sub_bin_sidelo[ibin][iq]->Draw("Hsame");
      IMpippim_wSid_n_Sp_sub_bin_sidehi[ibin][iq]->Draw("Hsame");
      
      //gr_IMpippim_wSid_n_Sp_sub_bin[ibin][iq] = new TGraphErrors();
      //HistToRorateGraph(IMpippim_wSid_n_Sp_sub_bin[ibin][iq],*gr_IMpippim_wSid_n_Sp_sub_bin[ibin][iq]);
      //gr_IMpippim_wSid_n_Sp_sub_bin[ibin][iq]->Draw("AP");
      //TBox *box = new TBox(0,anacuts::pipi_MIN_narrow,2000,anacuts::pipi_MAX_narrow);
      //box->SetFillColor(4);
      //box->SetFillStyle(3002);
      //box->Draw();

      cIMpippim_IMnpip_n_sub_bin_cut[ibin][iq]->cd(2);
      TPaveText *pt = new TPaveText(.05,.05,.95,.7);
      int binpipi_MIN = IMpippim_wSid_n_Sp_sub_bin[ibin][iq]->GetXaxis()->FindBin(anacuts::pipi_MIN_narrow);
      int binpipi_MAX = IMpippim_wSid_n_Sp_sub_bin[ibin][iq]->GetXaxis()->FindBin(anacuts::pipi_MAX_narrow);
      int binpipi_MIN_4sigma = IMpippim_wSid_n_Sp_sub_bin[ibin][iq]->GetXaxis()->FindBin(anacuts::pipi_MIN_narrow-4.0*anacuts::K0_sigma);
      double inteK0 = IMpippim_wSid_n_Sp_sub_bin_select[ibin][iq]->Integral();  //->Integral(binpipi_MIN,binpipi_MAX);
      double inteK0sidelo = IMpippim_wSid_n_Sp_sub_bin_sidelo[ibin][iq]->Integral();  //>Integral(binpipi_MIN_4sigma,binpipi_MIN-1);
      double inteK0sidehi = IMpippim_wSid_n_Sp_sub_bin_sidehi[ibin][iq]->Integral();  //>Integral(binpipi_MIN_4sigma,binpipi_MIN-1);
      double K0net = inteK0-inteK0sidelo-inteK0sidehi;
      int binnpip_MIN = IMnpip_wK0_n_sub_bin[ibin][iq]->GetXaxis()->FindBin(anacuts::Sigmap_MIN);
      int binnpip_MIN_2sigma = IMnpip_wK0_n_sub_bin[ibin][iq]->GetXaxis()->FindBin(anacuts::Sigmap_MIN-2.0*anacuts::Sigmap_sigma);
      int binnpip_MAX = IMnpip_wK0_n_sub_bin[ibin][iq]->GetXaxis()->FindBin(anacuts::Sigmap_MAX);
      int binnpip_MAX_2sigma = IMnpip_wK0_n_sub_bin[ibin][iq]->GetXaxis()->FindBin(anacuts::Sigmap_MAX+2.0*anacuts::Sigmap_sigma);
      double inteSp = IMnpip_wK0_n_sub_bin_select[ibin][iq]->Integral();  //->Integral(binnpip_MIN,binnpip_MAX);
      double inteSpsidelo = IMnpip_wK0_n_sub_bin_sidelo[ibin][iq]->Integral();  //(binnpip_MIN_2sigma,binnpip_MIN-1);
      double inteSpsidehi = IMnpip_wK0_n_sub_bin_sidehi[ibin][iq]->Integral();//(binnpip_MAX+1,binnpip_MAX_2sigma);
      double Spnet = inteSp - inteSpsidelo - inteSpsidehi;
      
      pt->AddText(Form("IM(n#pi^{-#}pi^{+}  %0.2f-%0.2f",1.0+ibin*1.0/nbintemplate,1.0+(ibin+1.0)/nbintemplate)); 
      pt->AddText(Form("K0 count %0.2f ",inteK0));
      pt->AddText(Form("K0 side low  %0.2f ",inteK0sidelo));
      pt->AddText(Form("K0 side high %0.2f ",inteK0sidehi));
      pt->AddText(Form("K0 net   %0.2f    ", inteK0-inteK0sidelo-inteK0sidehi));
      pt->AddText(Form("Sigma+ count %0.2f ",inteSp));
      pt->AddText(Form("Sigma+ side low   %0.2f ",inteSpsidelo));
      pt->AddText(Form("Sigma+ side high  %0.2f ",inteSpsidehi));
      pt->AddText(Form("Sigma+ net  %0.2f ", inteSp-inteSpsidelo-inteSpsidehi));
      pt->AddText(Form("K0 ratio in cross  %0.2f ",K0net/(K0net+Spnet))); 
      pt->Draw();
    }
  }

  TCanvas *cIMpippim_IMnpim_n_sub_bin[nbintemplate][nqcut];
  TH1D* IMnpim_n_sub_bin[nbintemplate][nqcut];
  TH1D* IMpippim_n_sub_bin_2[nbintemplate][nqcut];
  TGraphErrors *gr_IMpippim_n_sub_bin_2[nbintemplate][nqcut];
  TCanvas *cIMpippim_IMnpim_n_sub_bin_cut[nbintemplate][nqcut];
  TH1D* IMnpim_wK0_n_sub_bin[nbintemplate][nqcut];
  TH1D* IMnpim_wK0_n_sub_bin_select[nbintemplate][nqcut];
  TH1D* IMpippim_wSid_n_Sm_sub_bin[nbintemplate][nqcut];
  TGraphErrors *gr_IMpippim_wSid_n_Sm_sub_bin[nbintemplate][nqcut];
  
  for(unsigned int ibin=40;ibin<nbintemplate;ibin++){
    for(int iq=0;iq<nqcut;iq++){
      //first canvas to display IM(pi+pi-) vs IM(npi-) and their simple projection just to see the signal and other background 
      //no complicated cut at this moment
      cIMpippim_IMnpim_n_sub_bin[ibin][iq] = new TCanvas(Form("cIMpippim_IMnpim_n_sub_bin%d_%d",ibin,iq),Form("cIMpippim_IMnpim_n_sub_bin%d_%d",ibin,iq),800,800);
      cIMpippim_IMnpim_n_sub_bin[ibin][iq]->Divide(2,2,0,0);
      cIMpippim_IMnpim_n_sub_bin[ibin][iq]->cd(3);
      IMpippim_IMnpim_n_bin_sub[ibin][iq]->RebinX(4);
      IMpippim_IMnpim_n_bin_sub[ibin][iq]->RebinY(5);
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

      //2nd canvas to display the calculation of decomposion of K0, Sigma+,Sigma- 
      //2D plot in canvas(3) has cut to select K0 or Sigma-
      //projection plot in canvas(1) select Sigma- events candidate + (Sigma- & K0 ) overlap events
      //projection plot in canvas(4) select K0     events candidate + (Sigma- & K0 ) overlap events 
      cIMpippim_IMnpim_n_sub_bin_cut[ibin][iq]
      = new TCanvas(Form("cIMpippim_IMnpim_n_sub_bin_cut_%d_%d",ibin,iq),Form("cIMpippim_IMnpim_n_sub_bin_cut%d_%d",ibin,iq),800,800);
      cIMpippim_IMnpim_n_sub_bin_cut[ibin][iq]->Divide(2,2,0,0);
      cIMpippim_IMnpim_n_sub_bin_cut[ibin][iq]->cd(3);

      IMpippim_IMnpim_wK0orwSid_n_bin_sub[ibin][iq]->RebinX(4);
      IMpippim_IMnpim_wK0orwSid_n_bin_sub[ibin][iq]->RebinY(5);
      IMpippim_IMnpim_wK0orwSid_n_bin_sub[ibin][iq]->Draw("colz");
      cIMpippim_IMnpim_n_sub_bin_cut[ibin][iq]->cd(1);
      IMnpim_wK0_n_sub_bin[ibin][iq] = (TH1D*)IMpippim_IMnpim_wK0_n_bin_sub[ibin][iq]->ProjectionX(Form("IMnpim_wK0_n_sub_bin%d_%d",ibin,iq));
      IMnpim_wK0_n_sub_bin[ibin][iq] ->SetTitle(Form("IMnpim_wK0_n_sub_bin%d_%d",ibin,iq));
      IMnpim_wK0_n_sub_bin[ibin][iq]->Draw("HE");
      IMnpim_wK0_n_sub_bin_select[ibin][iq] = (TH1D*)IMnpim_wK0_n_sub_bin[ibin][iq]->Clone(Form("IMnpim_wK0_n_sub_bin%d_%d_select",ibin,iq));
      IMnpim_wK0_n_sub_bin_select[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmam_MIN,anacuts::Sigmam_MAX);
      IMnpim_wK0_n_sub_bin_select[ibin][iq]->SetLineColor(2);
      IMnpim_wK0_n_sub_bin_select[ibin][iq]->SetFillColor(2);
      IMnpim_wK0_n_sub_bin_select[ibin][iq]->Draw("HEsame");
      
      
      cIMpippim_IMnpim_n_sub_bin_cut[ibin][iq]->cd(4);
      IMpippim_wSid_n_Sm_sub_bin[ibin][iq]
      = (TH1D*)IMpippim_IMnpim_wSid_n_Sm_bin_sub[ibin][iq]->ProjectionY(Form("IMpippim_wSid_n_Sm_sub_bin%d_%d",ibin,iq));
      
      gr_IMpippim_wSid_n_Sm_sub_bin[ibin][iq] = new TGraphErrors();
      HistToRorateGraph(IMpippim_wSid_n_Sm_sub_bin[ibin][iq],*gr_IMpippim_wSid_n_Sm_sub_bin[ibin][iq]);
      gr_IMpippim_wSid_n_Sm_sub_bin[ibin][iq]->Draw("AP");
      
      cIMpippim_IMnpim_n_sub_bin_cut[ibin][iq]->cd(2);
      TPaveText *pt = new TPaveText(.05,.05,.95,.7);
      int binpipi_MIN = IMpippim_wSid_n_Sm_sub_bin[ibin][iq]->GetXaxis()->FindBin(anacuts::pipi_MIN_narrow);
      int binpipi_MAX = IMpippim_wSid_n_Sm_sub_bin[ibin][iq]->GetXaxis()->FindBin(anacuts::pipi_MAX_narrow);
      int binpipi_MIN_4sigma = IMpippim_wSid_n_Sm_sub_bin[ibin][iq]->GetXaxis()->FindBin(anacuts::pipi_MIN_narrow-4.0*anacuts::K0_sigma);
      double inteK0 = IMpippim_wSid_n_Sm_sub_bin[ibin][iq]->Integral(binpipi_MIN,binpipi_MAX);
      double inteK0sidelo = IMpippim_wSid_n_Sm_sub_bin[ibin][iq]->Integral(binpipi_MIN_4sigma,binpipi_MIN-1);
      double K0net = inteK0-inteK0sidelo;
      int binnpim_MIN = IMnpim_wK0_n_sub_bin[ibin][iq]->GetXaxis()->FindBin(anacuts::Sigmam_MIN);
      int binnpim_MIN_2sigma = IMnpim_wK0_n_sub_bin[ibin][iq]->GetXaxis()->FindBin(anacuts::Sigmam_MIN-2.0*anacuts::Sigmam_sigma);
      int binnpim_MAX = IMnpim_wK0_n_sub_bin[ibin][iq]->GetXaxis()->FindBin(anacuts::Sigmam_MAX);
      int binnpim_MAX_2sigma = IMnpim_wK0_n_sub_bin[ibin][iq]->GetXaxis()->FindBin(anacuts::Sigmam_MAX+2.0*anacuts::Sigmam_sigma);
      double inteSm = IMnpim_wK0_n_sub_bin[ibin][iq]->Integral(binnpim_MIN,binnpim_MAX);
      double inteSmsidelo = IMnpim_wK0_n_sub_bin[ibin][iq]->Integral(binnpim_MIN_2sigma,binnpim_MIN-1);
      double inteSmsidehi = IMnpim_wK0_n_sub_bin[ibin][iq]->Integral(binnpim_MAX+1,binnpim_MAX_2sigma);
      double Smnet = inteSm- inteSmsidelo - inteSmsidehi;

      pt->AddText(Form("K0 count %0.2f ",inteK0));
      pt->AddText(Form("K0 net   %0.2f ",K0net));
      pt->AddText(Form("Sigma+ count %0.2f ",inteSm));
      pt->AddText(Form("Sigma+ net   %0.2f ",Smnet));
      pt->AddText(Form("K0 ratio in cross %0.2f ",K0net/(K0net+Smnet))); 
      pt->Draw();
    }
  }

  std::cout << "Sigma- vs Sigma+ decompositon " << std::endl;

  TCanvas *cIMnpim_IMnpip_n_sub_bin[nbintemplate][nqcut];
  TH1D* IMnpip_n_sub_bin_2[nbintemplate][nqcut];
  TH1D* IMnpim_n_sub_bin_2[nbintemplate][nqcut];
  TGraphErrors *gr_IMnpim_n_sub_bin_2[nbintemplate][nqcut];
  
  TCanvas *cIMnpim_IMnpip_n_sub_bin_cut[nbintemplate][nqcut];
  TH1D* IMnpip_wSid_n_Sm_sub_bin[nbintemplate][nqcut];
  TH1D* IMnpip_wSid_n_Sm_sub_bin_select[nbintemplate][nqcut];
  TH1D* IMnpim_wSid_n_Sp_sub_bin[nbintemplate][nqcut];
  TGraphErrors *gr_IMnpim_wSid_n_Sp_sub_bin[nbintemplate][nqcut];
  
  for(unsigned int ibin=0;ibin<nbintemplate;ibin++){
    for(int iq=0;iq<nqcut;iq++){
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
      cIMnpim_IMnpip_n_sub_bin_cut[ibin][iq]->Divide(2,2,0,0);
      cIMnpim_IMnpip_n_sub_bin_cut[ibin][iq]->cd(3);

      IMnpim_IMnpip_wSid_n_bin_sub[ibin][iq]->RebinX(4);
      IMnpim_IMnpip_wSid_n_bin_sub[ibin][iq]->RebinY(4);
      IMnpim_IMnpip_wSid_n_bin_sub[ibin][iq]->Draw("colz");
      
      cIMnpim_IMnpip_n_sub_bin_cut[ibin][iq]->cd(1);
      IMnpip_wSid_n_Sm_sub_bin[ibin][iq] = (TH1D*)IMnpim_IMnpip_wSid_n_Sm_bin_sub[ibin][iq]->ProjectionX(Form("IMnpip_wSid_n_Sm_sub_bin%d_%d",ibin,iq));
      IMnpip_wSid_n_Sm_sub_bin[ibin][iq]->SetTitle(Form("IMnpip_wSid_n_Sm_sub_bin%d_%d",ibin,iq));
      IMnpip_wSid_n_Sm_sub_bin[ibin][iq]->Draw("HE");
      IMnpip_wSid_n_Sm_sub_bin_select[ibin][iq] = (TH1D*)IMnpip_wSid_n_Sm_sub_bin[ibin][iq]->Clone(Form("IMnpip_wSid_n_Sm_sub_bin%d_%d_select",ibin,iq));
      IMnpip_wSid_n_Sm_sub_bin_select[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmap_MIN,anacuts::Sigmap_MAX);
      IMnpip_wSid_n_Sm_sub_bin_select[ibin][iq]->SetLineColor(2);
      IMnpip_wSid_n_Sm_sub_bin_select[ibin][iq]->SetFillColor(2);
      IMnpip_wSid_n_Sm_sub_bin_select[ibin][iq]->Draw("HEsame");
      
      cIMnpim_IMnpip_n_sub_bin_cut[ibin][iq]->cd(4);
      IMnpim_wSid_n_Sp_sub_bin[ibin][iq]
      = (TH1D*)IMnpim_IMnpip_wSid_n_Sp_bin_sub[ibin][iq]->ProjectionY(Form("IMnpim_wSid_n_Sp_sub_bin%d_%d",ibin,iq));
      
      gr_IMnpim_wSid_n_Sp_sub_bin[ibin][iq] = new TGraphErrors();
      HistToRorateGraph(IMnpim_wSid_n_Sp_sub_bin[ibin][iq],*gr_IMnpim_wSid_n_Sp_sub_bin[ibin][iq]);
      gr_IMnpim_wSid_n_Sp_sub_bin[ibin][iq]->Draw("AP");
      
      cIMnpim_IMnpip_n_sub_bin_cut[ibin][iq]->cd(2);
      TPaveText *pt = new TPaveText(.05,.05,.95,.7);
      int binnpim_MIN = IMnpim_wSid_n_Sp_sub_bin[ibin][iq]->GetXaxis()->FindBin(anacuts::Sigmam_MIN);
      int binnpim_MAX = IMnpim_wSid_n_Sp_sub_bin[ibin][iq]->GetXaxis()->FindBin(anacuts::Sigmam_MAX);
      int binnpim_MIN_4sigma = IMnpim_wSid_n_Sp_sub_bin[ibin][iq]->GetXaxis()->FindBin(anacuts::Sigmam_MIN-4.0*anacuts::Sigmam_sigma);
      double inteSm = IMnpim_wSid_n_Sp_sub_bin[ibin][iq]->Integral(binnpim_MIN,binnpim_MAX);
      double inteSmsidelo = IMnpim_wSid_n_Sp_sub_bin[ibin][iq]->Integral(binnpim_MIN_4sigma,binnpim_MIN-1);
      double Smnet = inteSm-inteSmsidelo;
      int binnpip_MIN = IMnpip_wSid_n_Sm_sub_bin[ibin][iq]->GetXaxis()->FindBin(anacuts::Sigmap_MIN);
      int binnpip_MIN_2sigma = IMnpip_wSid_n_Sm_sub_bin[ibin][iq]->GetXaxis()->FindBin(anacuts::Sigmap_MIN-2.0*anacuts::Sigmap_sigma);
      int binnpip_MAX = IMnpip_wSid_n_Sm_sub_bin[ibin][iq]->GetXaxis()->FindBin(anacuts::Sigmap_MAX);
      int binnpip_MAX_2sigma = IMnpip_wSid_n_Sm_sub_bin[ibin][iq]->GetXaxis()->FindBin(anacuts::Sigmap_MAX+2.0*anacuts::Sigmap_sigma);
      double inteSp = IMnpip_wSid_n_Sm_sub_bin[ibin][iq]->Integral(binnpip_MIN,binnpip_MAX);
      double inteSpsidelo = IMnpip_wSid_n_Sm_sub_bin[ibin][iq]->Integral(binnpip_MIN_2sigma,binnpip_MIN-1);
      double inteSpsidehi = IMnpip_wSid_n_Sm_sub_bin[ibin][iq]->Integral(binnpip_MAX+1,binnpip_MAX_2sigma);
      double Spnet = inteSp- inteSpsidelo - inteSpsidehi;

      pt->AddText(Form("Sigma- count %0.2f ",inteSm));
      pt->AddText(Form("Sigma- net   %0.2f ",Smnet));
      pt->AddText(Form("Sigma+ count %0.2f ",inteSp));
      pt->AddText(Form("Sigma+ net   %0.2f ",Spnet));
      pt->AddText(Form("K0 ratio in cross %0.2f ",Smnet/(Smnet+Spnet))); 
      pt->Draw();
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
