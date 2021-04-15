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
bool Sidefar=false;
bool FitNoWeight=true;

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

void K0SigmaTempSpline()
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
  
  gROOT->SetBatch(1);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0.);  

  const unsigned int nbintemplate = 100;
  const int nqcut=4;
  const int qstart=0;
  TH2F* IMpippim_IMnpip_wSid_n_Sp_sub[nqcut];
  TH2F* IMpippim_IMnpip_wSid_n_Sp_bin_sub[nbintemplate][nqcut];
  TH2F* IMpippim_IMnpim_wSid_n_Sm_sub[nqcut];
  TH2F* IMpippim_IMnpim_wSid_n_Sm_bin_sub[nbintemplate][nqcut];
  TH2F* IMpippim_IMnpip_wK0_n_sub[nqcut];
  TH2F* IMpippim_IMnpip_wK0_n_bin_sub[nbintemplate][nqcut];
  TH2F* IMpippim_IMnpim_wK0_n_sub[nqcut];
  TH2F* IMpippim_IMnpim_wK0_n_bin_sub[nbintemplate][nqcut];
  TH2F* IMpippim_IMnpip_wK0orwSid_n_sub[nqcut];
  TH2F* IMpippim_IMnpim_wK0orwSid_n_sub[nqcut];
  TH2F* IMpippim_IMnpip_wK0orwSid_n_bin_sub[nbintemplate][nqcut];
  TH2F* IMpippim_IMnpim_wK0orwSid_n_bin_sub[nbintemplate][nqcut];
  TH2F* IMnpim_IMnpip_wK0orwSid_n_sub[nqcut];
  TH2F* IMnpim_IMnpip_wK0orwSid_n_bin_sub[nbintemplate][nqcut];
  TH2F* IMnpim_IMnpip_wSid_n_Sp_sub[nqcut];
  TH2F* IMnpim_IMnpip_wSid_n_Sp_bin_sub[nbintemplate][nqcut];
  TH2F* IMnpim_IMnpip_wSid_n_Sm_sub[nqcut];
  TH2F* IMnpim_IMnpip_wSid_n_Sm_bin_sub[nbintemplate][nqcut];
  
  //for the overlap of S+ & S- & K0 counting 
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
  for(int iq=qstart;iq<nqcut;iq++){
    IMpippim_IMnpip_wSid_n_Sp_sub[iq] = (TH2F*)fr[iq]->Get("IMpippim_IMnpip_wSid_n_Sp");
    if(RemoveNegative)IMpippim_IMnpip_wSid_n_Sp_sub[iq]->SetMinimum(0);

    IMpippim_IMnpim_wSid_n_Sm_sub[iq] = (TH2F*)fr[iq]->Get("IMpippim_IMnpim_wSid_n_Sm");
    if(RemoveNegative)IMpippim_IMnpim_wSid_n_Sm_sub[iq]->SetMinimum(0);
   
    IMpippim_IMnpip_wK0_n_sub[iq] = (TH2F*)fr[iq]->Get("IMpippim_IMnpip_wK0_n");
    if(RemoveNegative)IMpippim_IMnpip_wK0_n_sub[iq]->SetMinimum(0);

    IMpippim_IMnpim_wK0_n_sub[iq] = (TH2F*)fr[iq]->Get("IMpippim_IMnpim_wK0_n");
    if(RemoveNegative)IMpippim_IMnpim_wK0_n_sub[iq]->SetMinimum(0);

    IMpippim_IMnpip_wK0orwSid_n_sub[iq] = (TH2F*)fr[iq]->Get("IMpippim_IMnpip_wK0orwSid_n");
    if(RemoveNegative)IMpippim_IMnpip_wK0orwSid_n_sub[iq]->SetMinimum(0);
    
    cIMpippim_IMnpip_n_all[iq] = new TCanvas(Form("cIMpippim_IMnpip_n_all_%s",cqcut[iq]),Form("cIMpippim_IMnpip_n_all_%s",cqcut[iq]));
    cIMpippim_IMnpip_n_all[iq]->Divide(2,2);
    cIMpippim_IMnpip_n_all[iq]->cd(3);
    //IMpippim_IMnpip_wK0orwSid_n_sub[iq]->RebinX(4);
    //IMpippim_IMnpip_wK0orwSid_n_sub[iq]->RebinY(5);
    IMpippim_IMnpip_wK0orwSid_n_sub[iq]->GetXaxis()->SetRangeUser(1,1.8);
    IMpippim_IMnpip_wK0orwSid_n_sub[iq]->GetYaxis()->SetRangeUser(0.2,0.9);
    IMpippim_IMnpip_wK0orwSid_n_sub[iq]->Draw("colz");
    cIMpippim_IMnpip_n_all[iq]->cd(1);
    IMnpip_wK0_n_sub[iq] = (TH1D*)IMpippim_IMnpip_wK0_n_sub[iq]->ProjectionX(Form("IMnpip_wK0_n_sub_%d",iq));
    IMnpip_wK0_n_sub[iq]->SetTitle(Form("IMnpip_wK0_n_sub_%d",iq));
    IMnpip_wK0_n_sub[iq]->GetXaxis()->SetRangeUser(1,1.8);
    IMnpip_wK0_n_sub[iq]->Draw("HE");
    cIMpippim_IMnpip_n_all[iq]->cd(4);
    IMpippim_wSid_n_Sp_sub[iq]
    = (TH1D*)IMpippim_IMnpip_wSid_n_Sp_sub[iq]->ProjectionY(Form("IMpippim_wSid_n_Sp_sub_%d",iq));
    IMpippim_wSid_n_Sp_sub[iq]->SetTitle(Form("IMpippim_wSid_n_Sp_sub_%d",iq));
    IMpippim_wSid_n_Sp_sub[iq]->GetXaxis()->SetRangeUser(0.2,0.9);
    IMpippim_wSid_n_Sp_sub[iq]->Draw("HE");

    IMpippim_IMnpim_wK0orwSid_n_sub[iq] = (TH2F*)fr[iq]->Get("IMpippim_IMnpim_wK0orwSid_n");
    if(RemoveNegative)IMpippim_IMnpim_wK0orwSid_n_sub[iq]->SetMinimum(0);

    IMpippim_IMnpim_wK0orwSid_n_sub[iq] = (TH2F*)fr[iq]->Get("IMpippim_IMnpim_wK0orwSid_n");
    if(RemoveNegative)IMpippim_IMnpim_wK0orwSid_n_sub[iq]->SetMinimum(0);
    cIMpippim_IMnpim_n_all[iq] = new TCanvas(Form("cIMpippim_IMnpim_n_all_%s",cqcut[iq]),Form("cIMpippim_IMnpim_n_all_%s",cqcut[iq]));
    cIMpippim_IMnpim_n_all[iq]->Divide(2,2);
    cIMpippim_IMnpim_n_all[iq]->cd(3);
    //IMpippim_IMnpim_wK0orwSid_n_sub[iq]->RebinX(4);
    //IMpippim_IMnpim_wK0orwSid_n_sub[iq]->RebinY(5);
    IMpippim_IMnpim_wK0orwSid_n_sub[iq]->GetXaxis()->SetRangeUser(1,1.8);
    IMpippim_IMnpim_wK0orwSid_n_sub[iq]->GetYaxis()->SetRangeUser(0.2,0.9);
    IMpippim_IMnpim_wK0orwSid_n_sub[iq]->Draw("colz");
    cIMpippim_IMnpim_n_all[iq]->cd(1);
    IMnpim_wK0_n_sub[iq] = (TH1D*)IMpippim_IMnpim_wK0_n_sub[iq]->ProjectionX(Form("IMnpim_wK0_n_sub_%d",iq));
    IMnpim_wK0_n_sub[iq]->SetTitle(Form("IMnpim_wK0_n_sub_%d",iq));
    IMnpim_wK0_n_sub[iq]->GetXaxis()->SetRangeUser(1,1.8);
    IMnpim_wK0_n_sub[iq]->Draw("HE");
    cIMpippim_IMnpim_n_all[iq]->cd(4);
    IMpippim_wSid_n_Sm_sub[iq]
    = (TH1D*)IMpippim_IMnpim_wSid_n_Sm_sub[iq]->ProjectionY(Form("IMpippim_wSid_n_Sm_sub_%d",iq));
    IMpippim_wSid_n_Sm_sub[iq]->SetTitle(Form("IMpippim_wSid_n_Sm_sub_%d",iq));
    IMpippim_wSid_n_Sm_sub[iq]->GetXaxis()->SetRangeUser(0.2,0.9);
    IMpippim_wSid_n_Sm_sub[iq]->Draw("HE");

    IMnpim_IMnpip_wK0orwSid_n_sub[iq] = (TH2F*)fr[iq]->Get("IMnpim_IMnpip_dE_wK0orwSid_n");
    if(RemoveNegative)IMnpim_IMnpip_wK0orwSid_n_sub[iq]->SetMinimum(0);

    IMnpim_IMnpip_wSid_n_Sm_sub[iq] = (TH2F*)fr[iq]->Get("IMnpim_IMnpip_dE_wSid_n_Sm");
    if(RemoveNegative)IMnpim_IMnpip_wSid_n_Sm_sub[iq]->SetMinimum(0);
    
    IMnpim_IMnpip_wSid_n_Sp_sub[iq] = (TH2F*)fr[iq]->Get("IMnpim_IMnpip_dE_wSid_n_Sp");
    if(RemoveNegative)IMnpim_IMnpip_wSid_n_Sp_sub[iq]->SetMinimum(0);

    IMnpim_IMnpip_wSid_n_Sm_sub[iq] = (TH2F*)fr[iq]->Get("IMnpim_IMnpip_dE_wSid_n_Sm");
    if(RemoveNegative)IMnpim_IMnpip_wSid_n_Sm_sub[iq]->SetMinimum(0);
    
    cIMnpim_IMnpip_n_all[iq] = new TCanvas(Form("cIMnpim_IMnpip_n_all_%s",cqcut[iq]),Form("cIMnpim_IMnpip_n_all_%s",cqcut[iq]));
    cIMnpim_IMnpip_n_all[iq]->Divide(2,2);
    cIMnpim_IMnpip_n_all[iq]->cd(3);
    //IMnpim_IMnpip_wK0orwSid_n_sub[iq]->RebinX(4);
    //IMnpim_IMnpip_wK0orwSid_n_sub[iq]->RebinY(4);
    IMnpim_IMnpip_wK0orwSid_n_sub[iq]->GetXaxis()->SetRangeUser(1,1.8);
    IMnpim_IMnpip_wK0orwSid_n_sub[iq]->GetYaxis()->SetRangeUser(1,1.8);
    IMnpim_IMnpip_wK0orwSid_n_sub[iq]->Draw("colz");
    cIMnpim_IMnpip_n_all[iq]->cd(1);
    IMnpip_wSid_n_Sm_sub[iq] = (TH1D*)IMnpim_IMnpip_wSid_n_Sm_sub[iq]->ProjectionX(Form("IMnpip_wSid_n_Sm_sub_%d",iq));
    IMnpip_wSid_n_Sm_sub[iq]->SetTitle(Form("IMnpip_wSid_n_Sm_sub_%d",iq));
    IMnpip_wSid_n_Sm_sub[iq]->GetXaxis()->SetRangeUser(1,1.8);
    IMnpip_wSid_n_Sm_sub[iq]->Draw("HIST");
    cIMnpim_IMnpip_n_all[iq]->cd(4);
    IMnpim_wSid_n_Sp_sub[iq]
    = (TH1D*)IMnpim_IMnpip_wSid_n_Sp_sub[iq]->ProjectionY(Form("IMnpim_wSid_n_Sp_sub_%d",iq));
    IMnpim_wSid_n_Sp_sub[iq]->SetTitle(Form("IMnpim_wSid_n_Sp_sub_%d",iq));
    IMnpim_wSid_n_Sp_sub[iq]->GetXaxis()->SetRangeUser(1,1.8);
    IMnpim_wSid_n_Sp_sub[iq]->Draw("HIST");
    

    q_IMnpipi_wK0_wSid_n_SpSm_sub[iq] = (TH2F*)fr[iq]->Get("q_IMnpipi_wK0_wSid_n_SpSm");
    cq_IMnpipi_wK0_wSid_n_SpSm_sub[iq] = new TCanvas(Form("q_IMnpipi_wK0_wSid_n_SpSm_%s",cqcut[iq]),Form("q_IMnpipi_wK0_wSid_n_SpSm_%s",cqcut[iq]));
    IMnpipi_wK0_wSid_n_SpSm_sub[iq] = (TH1D*)q_IMnpipi_wK0_wSid_n_SpSm_sub[iq]->ProjectionX(Form("IMnpipi_wK0_wSid_n_SpSm_%d",iq));
    IMnpipi_wK0_wSid_n_SpSm_sub[iq]->Draw("HISTE");
    
    for(unsigned int ibin=0;ibin<nbintemplate;ibin++){
      OverlapCount[ibin][iq] = IMnpipi_wK0_wSid_n_SpSm_sub[iq]->GetBinContent(ibin+1);
      if(OverlapCount[ibin][iq]<0.0) OverlapCount[ibin][iq]=0.0;
    }

    for(unsigned int ibin=0;ibin<nbintemplate;ibin++){
      std::cout << iq << "  " << ibin << std::endl;

      IMpippim_IMnpip_wSid_n_Sp_bin_sub[ibin][iq] = (TH2F*)fr[iq]->Get(Form("IMpippim_IMnpip_wSid_n_Sp_bin%d",ibin));
      if(RemoveNegative)IMpippim_IMnpip_wSid_n_Sp_bin_sub[ibin][iq]->SetMinimum(0);

      IMpippim_IMnpim_wSid_n_Sm_bin_sub[ibin][iq] = (TH2F*)fr[iq]->Get(Form("IMpippim_IMnpim_wSid_n_Sm_bin%d",ibin));
      if(RemoveNegative)IMpippim_IMnpim_wSid_n_Sm_bin_sub[ibin][iq]->SetMinimum(0);

      IMpippim_IMnpip_wK0_n_bin_sub[ibin][iq] = (TH2F*)fr[iq]->Get(Form("IMpippim_IMnpip_wK0_n_bin%d",ibin));
      if(RemoveNegative)IMpippim_IMnpip_wK0_n_bin_sub[ibin][iq]->SetMinimum(0);

      IMpippim_IMnpim_wK0_n_bin_sub[ibin][iq] = (TH2F*)fr[iq]->Get(Form("IMpippim_IMnpim_wK0_n_bin%d",ibin));
      if(RemoveNegative)IMpippim_IMnpim_wK0_n_bin_sub[ibin][iq]->SetMinimum(0);

      IMpippim_IMnpip_wK0orwSid_n_bin_sub[ibin][iq] = (TH2F*)fr[iq]->Get(Form("IMpippim_IMnpip_wK0orwSid_n_bin%d",ibin));
      if(RemoveNegative)IMpippim_IMnpip_wK0orwSid_n_bin_sub[ibin][iq]->SetMinimum(0);
      
      IMpippim_IMnpim_wK0orwSid_n_bin_sub[ibin][iq] = (TH2F*)fr[iq]->Get(Form("IMpippim_IMnpim_wK0orwSid_n_bin%d",ibin));
      if(RemoveNegative)IMpippim_IMnpim_wK0orwSid_n_bin_sub[ibin][iq]->SetMinimum(0);

      float binlow=1.0+(float)ibin*1./nbintemplate;
      float binhigh=1.0+(float)(ibin+1.0)*1./nbintemplate;
      
      IMnpim_IMnpip_wK0orwSid_n_bin_sub[ibin][iq] = (TH2F*)fr[iq]->Get(Form("IMnpim_IMnpip_dE_wK0orwSid_n_bin%d",ibin));
      if(RemoveNegative)IMnpim_IMnpip_wK0orwSid_n_bin_sub[ibin][iq]->SetMinimum(0);
      
      IMnpim_IMnpip_wSid_n_Sp_bin_sub[ibin][iq] = (TH2F*)fr[iq]->Get(Form("IMnpim_IMnpip_wSid_n_Sp_bin%d",ibin));
      if(RemoveNegative)IMnpim_IMnpip_wSid_n_Sp_bin_sub[ibin][iq]->SetMinimum(0);

      IMnpim_IMnpip_wSid_n_Sm_bin_sub[ibin][iq] = (TH2F*)fr[iq]->Get(Form("IMnpim_IMnpip_wSid_n_Sm_bin%d",ibin));
      if(RemoveNegative)IMnpim_IMnpip_wSid_n_Sm_bin_sub[ibin][iq]->SetMinimum(0);
    }
  }
  
  
  std::cout << __LINE__ << std::endl;
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
  TH1D* IMnpip_wK0_n_sub_bin_lo[nbintemplate][nqcut];
  TH1D* IMnpip_wK0_n_sub_bin_hi[nbintemplate][nqcut];
  TH1D* IMpippim_wSid_n_Sp_sub_bin[nbintemplate][nqcut];
  TH1D* IMpippim_wSid_n_Sp_sub_bin_select[nbintemplate][nqcut];
  TH1D* IMpippim_wSid_n_Sp_sub_bin_lo[nbintemplate][nqcut];
  TH1D* IMpippim_wSid_n_Sp_sub_bin_lo2[nbintemplate][nqcut];//2 sigma to 6 sigma away when there is no higher side
  TH1D* IMpippim_wSid_n_Sp_sub_bin_hi[nbintemplate][nqcut];
  //IM(pi+pi-) vs IM(npi-) correlations
  TCanvas *cIMpippim_IMnpim_n_sub_bin[nbintemplate][nqcut];
  TH1D* IMnpim_n_sub_bin[nbintemplate][nqcut];
  TH1D* IMpippim_n_sub_bin_2[nbintemplate][nqcut];
  TGraphErrors *gr_IMpippim_n_sub_bin_2[nbintemplate][nqcut];
  TCanvas *cIMpippim_IMnpim_n_sub_bin_cut[nbintemplate][nqcut];
  TH1D* IMnpim_wK0_n_sub_bin[nbintemplate][nqcut];
  TH1D* IMnpim_wK0_n_sub_bin_select[nbintemplate][nqcut];
  TH1D* IMnpim_wK0_n_sub_bin_lo[nbintemplate][nqcut];
  TH1D* IMnpim_wK0_n_sub_bin_hi[nbintemplate][nqcut];
  TH1D* IMpippim_wSid_n_Sm_sub_bin[nbintemplate][nqcut];
  TH1D* IMpippim_wSid_n_Sm_sub_bin_select[nbintemplate][nqcut];
  TH1D* IMpippim_wSid_n_Sm_sub_bin_lo[nbintemplate][nqcut];
  TH1D* IMpippim_wSid_n_Sm_sub_bin_lo2[nbintemplate][nqcut];
  TH1D* IMpippim_wSid_n_Sm_sub_bin_hi[nbintemplate][nqcut];
  
  TCanvas *cIMnpim_IMnpip_n_sub_bin[nbintemplate][nqcut];
  TH1D* IMnpip_n_sub_bin_2[nbintemplate][nqcut];
  TH1D* IMnpim_n_sub_bin_2[nbintemplate][nqcut];
  TGraphErrors *gr_IMnpim_n_sub_bin_2[nbintemplate][nqcut];
  TCanvas *cIMnpim_IMnpip_n_sub_bin_cut[nbintemplate][nqcut];
  TH1D* IMnpip_wSid_n_Sm_sub_bin[nbintemplate][nqcut];
  TH1D* IMnpip_wSid_n_Sm_sub_bin_select[nbintemplate][nqcut];
  TH1D* IMnpip_wSid_n_Sm_sub_bin_lo[nbintemplate][nqcut];//2 sigma
  TH1D* IMnpip_wSid_n_Sm_sub_bin_lo2[nbintemplate][nqcut];//4 sigma
  TH1D* IMnpip_wSid_n_Sm_sub_bin_hi[nbintemplate][nqcut];//2 sigma
  TH1D* IMnpip_wSid_n_Sm_sub_bin_hi2[nbintemplate][nqcut];//4 simga
  TH1D* IMnpim_wSid_n_Sp_sub_bin[nbintemplate][nqcut];
  TH1D* IMnpim_wSid_n_Sp_sub_bin_select[nbintemplate][nqcut];
  TH1D* IMnpim_wSid_n_Sp_sub_bin_lo[nbintemplate][nqcut];//2 sigma
  TH1D* IMnpim_wSid_n_Sp_sub_bin_lo2[nbintemplate][nqcut];//4 sigma
  TH1D* IMnpim_wSid_n_Sp_sub_bin_hi[nbintemplate][nqcut];//2 sigma
  TH1D* IMnpim_wSid_n_Sp_sub_bin_hi2[nbintemplate][nqcut];//4 sigma
  
  for(int ig=0;ig<ngroup;ig++){
    for(int iq=qstart;iq<nqcut;iq++){
      if(ig!=2)IMnpipi_overlapdeco_K0[ig][iq] = new TH1D(Form("IMnpipi_overlapdeco_K0_g%d_%d",ig,iq),Form("IMnpipi_overlapdeco_K0_g%d_%d",ig,iq),100,1,2);
      if(ig!=1)IMnpipi_overlapdeco_Sp[ig][iq] = new TH1D(Form("IMnpipi_overlapdeco_Sp_g%d_%d",ig,iq),Form("IMnpipi_overlapdeco_Sp_g%d_%d",ig,iq),100,1,2);
      if(ig!=0)IMnpipi_overlapdeco_Sm[ig][iq] = new TH1D(Form("IMnpipi_overlapdeco_Sm_g%d_%d",ig,iq),Form("IMnpipi_overlapdeco_Sm_g%d_%d",ig,iq),100,1,2);
    }
  }

  for(unsigned int ibin=39;ibin<nbintemplate;ibin++){
    for(int iq=qstart;iq<nqcut;iq++){
      //std::cout << ibin << " " << iq << std::endl;
      //first canvas to display IM(pi+pi-) vs IM(npi+) and their simple projection just to see the signal and other background 
      //no complicated cut at this moment
      //2nd canvas to display the calculation of decomposion of K0, Sigma+,Sigma- 
      //2D plot in canvas(3) has cut to select K0 or Sigma+
      //projection plot in canvas(1) select Sigma+ events candidate + (Sigma+ & K0 ) overlap events
      //projection plot in canvas(4) select K0     events candidate + (Sigma+ & K0 ) overlap events 
      cIMpippim_IMnpip_n_sub_bin_cut[ibin][iq] 
      = new TCanvas(Form("cIMpippim_IMnpip_n_sub_bin_cut_%d_%d",ibin,iq),Form("cIMpippim_IMnpip_n_sub_bin_cut%d_%d",ibin,iq),800,800);
      cIMpippim_IMnpip_n_sub_bin_cut[ibin][iq]->Divide(2,2);
      cIMpippim_IMnpip_n_sub_bin_cut[ibin][iq]->cd(3);
      
      //IMpippim_IMnpip_wK0orwSid_n_bin_sub[ibin][iq]->RebinX(4);
      //IMpippim_IMnpip_wK0orwSid_n_bin_sub[ibin][iq]->RebinY(5);
      IMpippim_IMnpip_wK0orwSid_n_bin_sub[ibin][iq]->GetXaxis()->SetRangeUser(1,1.8);
      IMpippim_IMnpip_wK0orwSid_n_bin_sub[ibin][iq]->GetYaxis()->SetRangeUser(0.2,0.9);
      IMpippim_IMnpip_wK0orwSid_n_bin_sub[ibin][iq]->Draw("colz");
      cIMpippim_IMnpip_n_sub_bin_cut[ibin][iq]->cd(1);
      IMnpip_wK0_n_sub_bin[ibin][iq] = (TH1D*)IMpippim_IMnpip_wK0_n_bin_sub[ibin][iq]->ProjectionX(Form("IMnpip_wK0_n_sub_bin%d_%d",ibin,iq));
      IMnpip_wK0_n_sub_bin[ibin][iq]->SetTitle(Form("IMnpip_wK0_n_sub_bin%d_%d",ibin,iq));
      IMnpip_wK0_n_sub_bin[ibin][iq]->GetXaxis()->SetRangeUser(1,1.8);
      IMnpip_wK0_n_sub_bin[ibin][iq]->Draw("HIST");
      IMnpip_wK0_n_sub_bin_select[ibin][iq] = (TH1D*)IMnpip_wK0_n_sub_bin[ibin][iq]->Clone(Form("IMnpip_wK0_n_sub_bin%d_%d_select",ibin,iq));
      IMnpip_wK0_n_sub_bin_select[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmap_MIN,anacuts::Sigmap_MAX);
      IMnpip_wK0_n_sub_bin_select[ibin][iq]->SetLineColor(2);
      //IMnpip_wK0_n_sub_bin_select[ibin][iq]->SetFillColor(2);
      IMnpip_wK0_n_sub_bin_select[ibin][iq]->Draw("HISTsame");
      //std::cout << __LINE__ << std::endl;
      //IMnpip_wK0_n_sub_bin_lo[ibin][iq] = (TH1D*)IMnpip_wK0_n_sub_bin[ibin][iq]->Clone(Form("IMnpip_wK0_n_sub_bin%d_%d_lo",ibin,iq));
      //IMnpip_wK0_n_sub_bin_hi[ibin][iq] = (TH1D*)IMnpip_wK0_n_sub_bin[ibin][iq]->Clone(Form("IMnpip_wK0_n_sub_bin%d_%d_hi",ibin,iq));
      //IMnpip_wK0_n_sub_bin_lo[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmap_MIN-4.0*anacuts::Sigmap_sigma, anacuts::Sigmap_MIN-2.0*anacuts::Sigmap_sigma);
      //IMnpip_wK0_n_sub_bin_hi[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmap_MAX+2.0*anacuts::Sigmap_sigma, anacuts::Sigmap_MAX+4.0*anacuts::Sigmap_sigma);
      //IMnpip_wK0_n_sub_bin_lo[ibin][iq]->SetLineColor(4);
      //IMnpip_wK0_n_sub_bin_lo[ibin][iq]->SetFillColor(4);
      //IMnpip_wK0_n_sub_bin_hi[ibin][iq]->SetLineColor(4);
      //IMnpip_wK0_n_sub_bin_hi[ibin][iq]->SetFillColor(4);
      //double err_IMnpip_wK0=0.0;
      //if(ibin==44 || ibin==45){
      //  IMnpip_wK0_n_sub_bin_hi[ibin][iq]->Draw("HISTsame");
      //}else{
      //  IMnpip_wK0_n_sub_bin_lo[ibin][iq]->Draw("HISTsame");
      //  IMnpip_wK0_n_sub_bin_hi[ibin][iq]->Draw("HISTsame");
      //}
      cIMpippim_IMnpip_n_sub_bin_cut[ibin][iq]->cd(4);
      IMpippim_wSid_n_Sp_sub_bin[ibin][iq]
      = (TH1D*)IMpippim_IMnpip_wSid_n_Sp_bin_sub[ibin][iq]->ProjectionY(Form("IMpippim_wSid_n_Sp_sub_bin%d_%d",ibin,iq));
      IMpippim_wSid_n_Sp_sub_bin[ibin][iq]->SetTitle(Form("IMpippim_wSid_n_Sp_sub_bin%d_%d",ibin,iq));
      IMpippim_wSid_n_Sp_sub_bin[ibin][iq]->GetXaxis()->SetRangeUser(0.2,0.9);
      IMpippim_wSid_n_Sp_sub_bin[ibin][iq]->Draw("HIST");
      IMpippim_wSid_n_Sp_sub_bin_select[ibin][iq] = (TH1D*)IMpippim_wSid_n_Sp_sub_bin[ibin][iq]->Clone(Form("IMpippim_wSid_n_Sp_sub_bin_%d_%d_select",ibin,iq));
      IMpippim_wSid_n_Sp_sub_bin_select[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::pipi_MIN_narrow,anacuts::pipi_MAX_narrow);
      IMpippim_wSid_n_Sp_sub_bin_select[ibin][iq]->SetLineColor(2);
      //IMpippim_wSid_n_Sp_sub_bin_select[ibin][iq]->SetFillColor(2);
      IMpippim_wSid_n_Sp_sub_bin_select[ibin][iq]->Draw("HISTsame");
      //IMpippim_wSid_n_Sp_sub_bin_lo[ibin][iq] = (TH1D*)IMpippim_wSid_n_Sp_sub_bin[ibin][iq]->Clone(Form("IMpippim_wSid_n_Sp_sub_bin_%d_%d_lo",ibin,iq));
      //IMpippim_wSid_n_Sp_sub_bin_lo2[ibin][iq] = (TH1D*)IMpippim_wSid_n_Sp_sub_bin[ibin][iq]->Clone(Form("IMpippim_wSid_n_Sp_sub_bin_%d_%d_lo2",ibin,iq));
      //IMpippim_wSid_n_Sp_sub_bin_hi[ibin][iq] = (TH1D*)IMpippim_wSid_n_Sp_sub_bin[ibin][iq]->Clone(Form("IMpippim_wSid_n_Sp_sub_bin_%d_%d_hi",ibin,iq));
      //IMpippim_wSid_n_Sp_sub_bin_lo[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::pipi_MIN_narrow-4.0*anacuts::K0_sigma,anacuts::pipi_MIN_narrow-2.0*anacuts::K0_sigma);
      //IMpippim_wSid_n_Sp_sub_bin_lo2[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::pipi_MIN_narrow-6.0*anacuts::K0_sigma,anacuts::pipi_MIN_narrow-2.0*anacuts::K0_sigma);
      //if(ibin==45){
      //  IMpippim_wSid_n_Sp_sub_bin_lo2[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::pipi_MIN_narrow-14.0*anacuts::K0_sigma,anacuts::pipi_MIN_narrow-10.0*anacuts::K0_sigma);
      //}else{
      //  IMpippim_wSid_n_Sp_sub_bin_lo2[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::pipi_MIN_narrow-6.0*anacuts::K0_sigma,anacuts::pipi_MIN_narrow-2.0*anacuts::K0_sigma);
      //}
      //IMpippim_wSid_n_Sp_sub_bin_hi[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::pipi_MAX_narrow+2.0*anacuts::K0_sigma,anacuts::pipi_MAX_narrow+4.0*anacuts::K0_sigma);
      //IMpippim_wSid_n_Sp_sub_bin_lo[ibin][iq]->SetLineColor(4);
      //IMpippim_wSid_n_Sp_sub_bin_lo[ibin][iq]->SetFillColor(4);
      //IMpippim_wSid_n_Sp_sub_bin_lo2[ibin][iq]->SetLineColor(4);
      //IMpippim_wSid_n_Sp_sub_bin_lo2[ibin][iq]->SetFillColor(4);
      //IMpippim_wSid_n_Sp_sub_bin_hi[ibin][iq]->SetLineColor(4);
      //IMpippim_wSid_n_Sp_sub_bin_hi[ibin][iq]->SetFillColor(4);

      //if(ibin>51){
      //  IMpippim_wSid_n_Sp_sub_bin_lo[ibin][iq]->Draw("HISTsame");
      //  IMpippim_wSid_n_Sp_sub_bin_hi[ibin][iq]->Draw("HISTsame");
      //}else{
      //  IMpippim_wSid_n_Sp_sub_bin_lo2[ibin][iq]->Draw("HISTsame");
      //}
      
      //gr_IMpippim_wSid_n_Sp_sub_bin[ibin][iq] = new TGraphErrors();
      //HistToRorateGraph(IMpippim_wSid_n_Sp_sub_bin[ibin][iq],*gr_IMpippim_wSid_n_Sp_sub_bin[ibin][iq]);
      //gr_IMpippim_wSid_n_Sp_sub_bin[ibin][iq]->Draw("AP");
      //TBox *box = new TBox(0,anacuts::pipi_MIN_narrow,2000,anacuts::pipi_MAX_narrow);
      //box->SetFillColor(4);
      //box->SetFillStyle(3002);
      //box->Draw();
       
      /*
      cIMpippim_IMnpip_n_sub_bin_cut[ibin][iq]->cd(2);
      TPaveText *pt = new TPaveText(.05,.05,.95,.7);
      double inteSp = IMpippim_wSid_n_Sp_sub_bin_select[ibin][iq]->Integral();  //->Integral(binpipi_MIN,binpipi_MAX);
      double inteSplo = IMpippim_wSid_n_Sp_sub_bin_lo[ibin][iq]->Integral();  //>Integral(binpipi_MIN_4sigma,binpipi_MIN-1);
      double inteSplo2 = IMpippim_wSid_n_Sp_sub_bin_lo2[ibin][iq]->Integral();  //>Integral(binpipi_MIN_4sigma,binpipi_MIN-1);
      double inteSphi = IMpippim_wSid_n_Sp_sub_bin_hi[ibin][iq]->Integral();  //>Integral(binpipi_MIN_4sigma,binpipi_MIN-1);
      double Spestimated=0.0;
      double Spoverlap = 0.0; 
      if(ibin>51) Spestimated = inteSplo+inteSphi;
      else        Spestimated = inteSplo2;
      Spoverlap = inteSp-Spestimated;
      if(Spestimated<0.0) Spestimated=0.0;
      if(Spoverlap<0.0) {
        Spoverlap = 0.0;
        Spestimated=0.0;
      }
      double inteK0 = IMnpip_wK0_n_sub_bin_select[ibin][iq]->Integral();  //->Integral(binnpip_MIN,binnpip_MAX);
      double inteK0lo = IMnpip_wK0_n_sub_bin_lo[ibin][iq]->Integral();  //(binnpip_MIN_2sigma,binnpip_MIN-1);
      double inteK0hi = IMnpip_wK0_n_sub_bin_hi[ibin][iq]->Integral();//(binnpip_MAX+1,binnpip_MAX_2sigma);
      double K0estimated = 0.0;
      if(ibin==44 || ibin==45) K0estimated = inteK0hi*2.0;
      else K0estimated = inteK0lo + inteK0hi;
      double K0overlap = inteK0 - K0estimated;
      if(K0estimated<0.0) K0estimated = 0.0;
      if(K0overlap<0.0){
        K0overlap = 0.0;
        K0estimated = 0.0;
      }
      pt->AddText(Form("IM(n#pi^{-}#pi^{+})  %0.2f-%0.2f %s",1.0+ibin*1.0/nbintemplate,1.0+(ibin+1.0)/nbintemplate,cqcut[iq])); 
      pt->AddText(Form("Sigma+ count     %0.2f ",inteSp));
      if(ibin>51){
        pt->AddText(Form("Sigma+ side low  %0.2f ",inteSplo));
        pt->AddText(Form("Sigma+ side high %0.2f ",inteSphi));
      }else{
        pt->AddText(Form("Sigma+ side low (4 sigma)  %0.2f ",inteSplo2));
      }
      pt->AddText(Form("Sigma+ estimated %0.2f ", Spestimated));
      pt->AddText(Form("K0 count %0.2f ",inteK0));
      if(ibin==44 || ibin==45){
        pt->AddText(Form("K0 side high*2  %0.2f ",inteK0hi*2));
      }else{
        pt->AddText(Form("K0 side low   %0.2f ",inteK0lo));
        pt->AddText(Form("K0 side high  %0.2f ",inteK0hi));
      }
      pt->AddText(Form("K0 estimated %0.2f ", K0estimated));
      pt->Draw();
      IMnpipi_overlapdeco_K0[0][iq]->SetBinContent(ibin,K0estimated);
      IMnpipi_overlapdeco_K0[0][iq]->SetBinError(ibin,sqrt(K0estimated));
      IMnpipi_overlapdeco_Sp[0][iq]->SetBinContent(ibin,Spestimated);
      IMnpipi_overlapdeco_Sp[0][iq]->SetBinError(ibin,sqrt(Spestimated));
      */
      /*
      //first canvas to display IM(pi+pi-) vs IM(npi-) and their simple projection just to see the signal and other background 
      //no complicated cut at this moment
      cIMpippim_IMnpim_n_sub_bin[ibin][iq] = new TCanvas(Form("cIMpippim_IMnpim_n_sub_bin%d_%d",ibin,iq),Form("cIMpippim_IMnpim_n_sub_bin%d_%d",ibin,iq),800,800);
      cIMpippim_IMnpim_n_sub_bin[ibin][iq]->Divide(2,2,0,0);
      cIMpippim_IMnpim_n_sub_bin[ibin][iq]->cd(3);
      IMpippim_IMnpim_n_bin_sub[ibin][iq]->RebinX(4);
      IMpippim_IMnpim_n_bin_sub[ibin][iq]->RebinY(5);
      IMpippim_IMnpim_n_bin_sub[ibin][iq]->GetXaxis()->SetRangeUser(1,1.8);
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
      //std::cout << __LINE__ << std::endl;
      cIMpippim_IMnpim_n_sub_bin_cut[ibin][iq]
      = new TCanvas(Form("cIMpippim_IMnpim_n_sub_bin_cut_%d_%d",ibin,iq),Form("cIMpippim_IMnpim_n_sub_bin_cut%d_%d",ibin,iq),800,800);
      cIMpippim_IMnpim_n_sub_bin_cut[ibin][iq]->Divide(2,2);
      cIMpippim_IMnpim_n_sub_bin_cut[ibin][iq]->cd(3);

      //IMpippim_IMnpim_wK0orwSid_n_bin_sub[ibin][iq]->RebinX(4);
      //IMpippim_IMnpim_wK0orwSid_n_bin_sub[ibin][iq]->RebinY(5);
      IMpippim_IMnpim_wK0orwSid_n_bin_sub[ibin][iq]->GetXaxis()->SetRangeUser(1,1.8);
      IMpippim_IMnpim_wK0orwSid_n_bin_sub[ibin][iq]->GetYaxis()->SetRangeUser(0.2,0.9);
      IMpippim_IMnpim_wK0orwSid_n_bin_sub[ibin][iq]->Draw("colz");
      cIMpippim_IMnpim_n_sub_bin_cut[ibin][iq]->cd(1);
      IMnpim_wK0_n_sub_bin[ibin][iq] = (TH1D*)IMpippim_IMnpim_wK0_n_bin_sub[ibin][iq]->ProjectionX(Form("IMnpim_wK0_n_sub_bin%d_%d",ibin,iq));
      IMnpim_wK0_n_sub_bin[ibin][iq]->SetTitle(Form("IMnpim_wK0_n_sub_bin%d_%d",ibin,iq));
      IMnpim_wK0_n_sub_bin[ibin][iq]->GetXaxis()->SetRangeUser(1,1.8);
      IMnpim_wK0_n_sub_bin[ibin][iq]->Draw("HIST");
      IMnpim_wK0_n_sub_bin_select[ibin][iq] = (TH1D*)IMnpim_wK0_n_sub_bin[ibin][iq]->Clone(Form("IMnpim_wK0_n_sub_bin%d_%d_select",ibin,iq));
      IMnpim_wK0_n_sub_bin_select[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmam_MIN,anacuts::Sigmam_MAX);
      IMnpim_wK0_n_sub_bin_select[ibin][iq]->SetLineColor(2);
      //IMnpim_wK0_n_sub_bin_select[ibin][iq]->SetFillColor(2);
      IMnpim_wK0_n_sub_bin_select[ibin][iq]->Draw("HISTsame");
      //IMnpim_wK0_n_sub_bin_lo[ibin][iq] = (TH1D*)IMnpim_wK0_n_sub_bin[ibin][iq]->Clone(Form("IMnpim_wK0_n_sub_bin%d_%d_selectlo",ibin,iq));
      //IMnpim_wK0_n_sub_bin_hi[ibin][iq] = (TH1D*)IMnpim_wK0_n_sub_bin[ibin][iq]->Clone(Form("IMnpim_wK0_n_sub_bin%d_%d_selecthi",ibin,iq));
      //IMnpim_wK0_n_sub_bin_lo[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmam_MIN-4.0*anacuts::Sigmam_sigma,anacuts::Sigmam_MIN-2.0*anacuts::Sigmam_sigma);
      //IMnpim_wK0_n_sub_bin_hi[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmam_MAX+2.0*anacuts::Sigmam_sigma,anacuts::Sigmam_MAX+4.0*anacuts::Sigmam_sigma);
      //IMnpim_wK0_n_sub_bin_lo[ibin][iq]->SetLineColor(4);
      //IMnpim_wK0_n_sub_bin_lo[ibin][iq]->SetFillColor(4);
      //IMnpim_wK0_n_sub_bin_hi[ibin][iq]->SetLineColor(4);
      //IMnpim_wK0_n_sub_bin_hi[ibin][iq]->SetFillColor(4);
      //if(ibin==44 || ibin==45){
      //  IMnpim_wK0_n_sub_bin_hi[ibin][iq]->Draw("HISTsame");
      //}else{
      //  IMnpim_wK0_n_sub_bin_lo[ibin][iq]->Draw("HISTsame");
      //  IMnpim_wK0_n_sub_bin_hi[ibin][iq]->Draw("HISTsame");
      //}
      
      //std::cout << __LINE__ << std::endl;
      cIMpippim_IMnpim_n_sub_bin_cut[ibin][iq]->cd(4);
      IMpippim_wSid_n_Sm_sub_bin[ibin][iq]
      = (TH1D*)IMpippim_IMnpim_wSid_n_Sm_bin_sub[ibin][iq]->ProjectionY(Form("IMpippim_wSid_n_Sm_sub_bin%d_%d",ibin,iq));
      IMpippim_wSid_n_Sm_sub_bin[ibin][iq]->SetTitle(Form("IMpippim_wSid_n_Sm_sub_bin%d_%d",ibin,iq));
      IMpippim_wSid_n_Sm_sub_bin[ibin][iq]->Draw("HIST");
      IMpippim_wSid_n_Sm_sub_bin_select[ibin][iq] = (TH1D*)IMpippim_wSid_n_Sm_sub_bin[ibin][iq]->Clone(Form("IMpippim_wSid_n_Sm_sub_bin_%d_%d_select",ibin,iq));
      IMpippim_wSid_n_Sm_sub_bin_select[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::pipi_MIN_narrow,anacuts::pipi_MAX_narrow);
      IMpippim_wSid_n_Sm_sub_bin_select[ibin][iq]->SetLineColor(2);
      //IMpippim_wSid_n_Sm_sub_bin_select[ibin][iq]->SetFillColor(2);
      IMpippim_wSid_n_Sm_sub_bin_select[ibin][iq]->Draw("HISTsame");
      //IMpippim_wSid_n_Sm_sub_bin_lo[ibin][iq] = (TH1D*)IMpippim_wSid_n_Sm_sub_bin[ibin][iq]->Clone(Form("IMpippim_wSid_n_Sm_sub_bin_%d_%d_lo",ibin,iq));
      //IMpippim_wSid_n_Sm_sub_bin_lo2[ibin][iq] = (TH1D*)IMpippim_wSid_n_Sm_sub_bin[ibin][iq]->Clone(Form("IMpippim_wSid_n_Sm_sub_bin_%d_%d_lo2",ibin,iq));
      //IMpippim_wSid_n_Sm_sub_bin_hi[ibin][iq] = (TH1D*)IMpippim_wSid_n_Sm_sub_bin[ibin][iq]->Clone(Form("IMpippim_wSid_n_Sm_sub_bin_%d_%d_hi",ibin,iq));
      //IMpippim_wSid_n_Sm_sub_bin_lo[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::pipi_MIN_narrow-4.0*anacuts::K0_sigma,anacuts::pipi_MIN_narrow-2.0*anacuts::K0_sigma);
      //if(ibin==45){
      //  IMpippim_wSid_n_Sm_sub_bin_lo2[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::pipi_MIN_narrow-14.0*anacuts::K0_sigma,anacuts::pipi_MIN_narrow-10.0*anacuts::K0_sigma);
      //}else{
      //  IMpippim_wSid_n_Sm_sub_bin_lo2[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::pipi_MIN_narrow-6.0*anacuts::K0_sigma,anacuts::pipi_MIN_narrow-2.0*anacuts::K0_sigma);
      //}
        
      //IMpippim_wSid_n_Sm_sub_bin_hi[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::pipi_MAX_narrow+2.0*anacuts::K0_sigma,anacuts::pipi_MAX_narrow+4.0*anacuts::K0_sigma);
      //IMpippim_wSid_n_Sm_sub_bin_lo[ibin][iq]->SetLineColor(4);
      //IMpippim_wSid_n_Sm_sub_bin_lo[ibin][iq]->SetFillColor(4);
      //IMpippim_wSid_n_Sm_sub_bin_lo2[ibin][iq]->SetLineColor(4);
      //IMpippim_wSid_n_Sm_sub_bin_lo2[ibin][iq]->SetFillColor(4);
      //IMpippim_wSid_n_Sm_sub_bin_hi[ibin][iq]->SetLineColor(4);
      //IMpippim_wSid_n_Sm_sub_bin_hi[ibin][iq]->SetFillColor(4);
      //if(ibin>51){
      //  IMpippim_wSid_n_Sm_sub_bin_lo[ibin][iq]->Draw("HISTsame");
      //  IMpippim_wSid_n_Sm_sub_bin_hi[ibin][iq]->Draw("HISTsame");
      //}else{
      //  IMpippim_wSid_n_Sm_sub_bin_lo2[ibin][iq]->Draw("HISTsame");
      //}

      //gr_IMpippim_wSid_n_Sm_sub_bin[ibin][iq] = new TGraphErrors();
      //HistToRorateGraph(IMpippim_wSid_n_Sm_sub_bin[ibin][iq],*gr_IMpippim_wSid_n_Sm_sub_bin[ibin][iq]);
      //gr_IMpippim_wSid_n_Sm_sub_bin[ibin][iq]->Draw("AP");
      /*
      cIMpippim_IMnpim_n_sub_bin_cut[ibin][iq]->cd(2);
      TPaveText *pt2 = new TPaveText(.05,.05,.95,.7);
      double inteSm_g2 = IMpippim_wSid_n_Sm_sub_bin_select[ibin][iq]->Integral();
      double inteSmlo_g2 = IMpippim_wSid_n_Sm_sub_bin_lo[ibin][iq]->Integral();
      double inteSmlo2_g2 = IMpippim_wSid_n_Sm_sub_bin_lo2[ibin][iq]->Integral();
      double inteSmhi_g2 = IMpippim_wSid_n_Sm_sub_bin_hi[ibin][iq]->Integral();
      double Smoverlap_g2 = 0.0;
      double Smestimated_g2 = 0.0;
      if(ibin>51) Smestimated_g2 = inteSmlo_g2+inteSmhi_g2;
      else        Smestimated_g2 = inteSmlo2_g2;
      Smoverlap_g2 = inteSm_g2-Smestimated_g2;
      if(Smestimated_g2<0.0) Smestimated_g2 = 0.0;
      if(Smoverlap_g2<0.0){
        Smoverlap_g2 = 0.0; 
        Smestimated_g2 = 0.0;
      }
      double inteK0_g2 = IMnpim_wK0_n_sub_bin_select[ibin][iq]->Integral();
      double inteK0lo_g2 = IMnpim_wK0_n_sub_bin_lo[ibin][iq]->Integral();
      double inteK0hi_g2 = IMnpim_wK0_n_sub_bin_hi[ibin][iq]->Integral();
      double K0estimated_g2 = 0.0;
      if(ibin==44 || ibin==45) K0estimated_g2 = inteK0hi_g2*2.0;
      else K0estimated_g2 = inteK0lo_g2 + inteK0hi_g2;
      double K0overlap_g2 = inteK0_g2 - K0estimated_g2;
      if(K0estimated_g2 <0.0) K0estimated_g2 = 0.0;
      if(K0overlap_g2<0.0) {
         K0overlap_g2 = 0.0;
         K0estimated_g2 = 0.0;
      }
      pt2->AddText(Form("IM(n#pi^{-}#pi^{+})  %0.2f-%0.2f %s",1.0+ibin*1.0/nbintemplate,1.0+(ibin+1.0)/nbintemplate,cqcut[iq])); 
      pt2->AddText(Form("Sigma- count %0.2f ",inteSm_g2));// 
      if(ibin>51){
        pt2->AddText(Form("Sigma- side low  %0.2f ",inteSmlo_g2));
        pt2->AddText(Form("Sigma- side high %0.2f ",inteSmhi_g2));
      }else{
        pt2->AddText(Form("Sigma- side low (4 sigma)  %0.2f ",inteSmlo2_g2));
      }
      pt2->AddText(Form("Sigma- estimated  %0.2f ",Smestimated_g2)); //Sm in (K0 ^ Sm)     
      pt2->AddText(Form("K0 count %0.2f ",inteK0_g2));
      if(ibin==44 || ibin==45){
        pt2->AddText(Form("K0 side high*2  %0.2f ",inteK0hi_g2*2.0));
      }else{
        pt2->AddText(Form("K0 side low   %0.2f ",inteK0lo_g2));
        pt2->AddText(Form("K0 side high  %0.2f ",inteK0hi_g2));
      }
      pt2->AddText(Form("K0 estimated  %0.2f ", K0estimated_g2));
      pt2->Draw();
      IMnpipi_overlapdeco_K0[1][iq]->SetBinContent(ibin,K0estimated_g2);
      IMnpipi_overlapdeco_K0[1][iq]->SetBinError(ibin,sqrt(K0estimated_g2));
      IMnpipi_overlapdeco_Sm[1][iq]->SetBinContent(ibin,Smestimated_g2);
      IMnpipi_overlapdeco_Sm[1][iq]->SetBinError(ibin,sqrt(Smestimated_g2));
      */     

      //std::cout << __LINE__ << std::endl;
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
      IMnpim_IMnpip_wK0orwSid_n_bin_sub[ibin][iq]->GetXaxis()->SetRangeUser(1,1.8);
      IMnpim_IMnpip_wK0orwSid_n_bin_sub[ibin][iq]->GetYaxis()->SetRangeUser(1,1.8);
      IMnpim_IMnpip_wK0orwSid_n_bin_sub[ibin][iq]->Draw("colz");
      
      cIMnpim_IMnpip_n_sub_bin_cut[ibin][iq]->cd(1);
      IMnpip_wSid_n_Sm_sub_bin[ibin][iq] = (TH1D*)IMnpim_IMnpip_wSid_n_Sm_bin_sub[ibin][iq]->ProjectionX(Form("IMnpip_wSid_n_Sm_sub_bin%d_%d",ibin,iq));
      IMnpip_wSid_n_Sm_sub_bin[ibin][iq]->SetTitle(Form("IMnpip_wSid_n_Sm_sub_bin%d_%d",ibin,iq));
      IMnpip_wSid_n_Sm_sub_bin[ibin][iq]->GetXaxis()->SetRangeUser(1,1.8);
      IMnpip_wSid_n_Sm_sub_bin[ibin][iq]->Draw("HIST");

      IMnpip_wSid_n_Sm_sub_bin_select[ibin][iq] = (TH1D*)IMnpip_wSid_n_Sm_sub_bin[ibin][iq]->Clone(Form("IMnpip_wSid_n_Sm_sub_bin%d_%d_select",ibin,iq));
      IMnpip_wSid_n_Sm_sub_bin_select[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmap_MIN,anacuts::Sigmap_MAX);
      IMnpip_wSid_n_Sm_sub_bin_select[ibin][iq]->SetLineColor(2);
      //IMnpip_wSid_n_Sm_sub_bin_select[ibin][iq]->SetFillColor(2);
      IMnpip_wSid_n_Sm_sub_bin_select[ibin][iq]->Draw("HISTsame");
      //std::cout << __LINE__ << std::endl;
      cIMnpim_IMnpip_n_sub_bin_cut[ibin][iq]->cd(4);
      IMnpim_wSid_n_Sp_sub_bin[ibin][iq]
      = (TH1D*)IMnpim_IMnpip_wSid_n_Sp_bin_sub[ibin][iq]->ProjectionY(Form("IMnpim_wSid_n_Sp_sub_bin%d_%d",ibin,iq));
      IMnpim_wSid_n_Sp_sub_bin[ibin][iq]->SetTitle(Form("IMnpim_wSid_n_Sp_sub_bin%d_%d",ibin,iq));
      IMnpim_wSid_n_Sp_sub_bin[ibin][iq]->GetXaxis()->SetRangeUser(1,1.8);
      IMnpim_wSid_n_Sp_sub_bin[ibin][iq]->Draw("HIST");
      IMnpim_wSid_n_Sp_sub_bin_select[ibin][iq] = (TH1D*)IMnpim_wSid_n_Sp_sub_bin[ibin][iq]->Clone(Form("IMnpim_wSid_n_Sp_sub_bin%d_%d_select",ibin,iq));
      IMnpim_wSid_n_Sp_sub_bin_select[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmam_MIN,anacuts::Sigmam_MAX);
      IMnpim_wSid_n_Sp_sub_bin_select[ibin][iq]->SetLineColor(2);
      //IMnpim_wSid_n_Sp_sub_bin_select[ibin][iq]->SetFillColor(2);
      IMnpim_wSid_n_Sp_sub_bin_select[ibin][iq]->Draw("HISTsame");
      //IMnpim_wSid_n_Sp_sub_bin_lo[ibin][iq] = (TH1D*)IMnpim_wSid_n_Sp_sub_bin[ibin][iq]->Clone(Form("IMnpim_wSid_n_Sp_sub_bin%d_%d_lo",ibin,iq));
      //IMnpim_wSid_n_Sp_sub_bin_lo2[ibin][iq] = (TH1D*)IMnpim_wSid_n_Sp_sub_bin[ibin][iq]->Clone(Form("IMnpim_wSid_n_Sp_sub_bin%d_%d_lo2",ibin,iq));
      //IMnpim_wSid_n_Sp_sub_bin_hi[ibin][iq] = (TH1D*)IMnpim_wSid_n_Sp_sub_bin[ibin][iq]->Clone(Form("IMnpim_wSid_n_Sp_sub_bin%d_%d_hi",ibin,iq));
      //IMnpim_wSid_n_Sp_sub_bin_hi2[ibin][iq] = (TH1D*)IMnpim_wSid_n_Sp_sub_bin[ibin][iq]->Clone(Form("IMnpim_wSid_n_Sp_sub_bin%d_%d_hi2",ibin,iq));
      //IMnpim_wSid_n_Sp_sub_bin_lo[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmam_MIN-4.0*anacuts::Sigmam_sigma,anacuts::Sigmam_MIN-2.0*anacuts::Sigmam_sigma);
      //IMnpim_wSid_n_Sp_sub_bin_lo2[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmam_MIN-6.0*anacuts::Sigmam_sigma,anacuts::Sigmam_MIN-2.0*anacuts::Sigmam_sigma);
      //if(ibin==42){
      //  IMnpim_wSid_n_Sp_sub_bin_hi[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmam_MAX+1.0*anacuts::Sigmam_sigma,anacuts::Sigmam_MAX+3.0*anacuts::Sigmam_sigma);
      //}else{
      //  IMnpim_wSid_n_Sp_sub_bin_hi[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmam_MAX+2.0*anacuts::Sigmam_sigma,anacuts::Sigmam_MAX+4.0*anacuts::Sigmam_sigma);
      //}  
      //IMnpim_wSid_n_Sp_sub_bin_hi2[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmam_MAX+2.0*anacuts::Sigmam_sigma,anacuts::Sigmam_MAX+6.0*anacuts::Sigmam_sigma);
      //IMnpim_wSid_n_Sp_sub_bin_lo[ibin][iq]->SetLineColor(4);
      //IMnpim_wSid_n_Sp_sub_bin_lo[ibin][iq]->SetFillColor(4);
      //IMnpim_wSid_n_Sp_sub_bin_lo2[ibin][iq]->SetLineColor(4);
      //IMnpim_wSid_n_Sp_sub_bin_lo2[ibin][iq]->SetFillColor(4);
      //IMnpim_wSid_n_Sp_sub_bin_hi[ibin][iq]->SetLineColor(4);
      //IMnpim_wSid_n_Sp_sub_bin_hi[ibin][iq]->SetFillColor(4);
      //IMnpim_wSid_n_Sp_sub_bin_hi2[ibin][iq]->SetLineColor(4);
      //IMnpim_wSid_n_Sp_sub_bin_hi2[ibin][iq]->SetFillColor(4);
      //if(ibin <42){
      //  IMnpim_wSid_n_Sp_sub_bin_lo2[ibin][iq]->Draw("HISTsame");
      //}else if( ibin == 42){
      //  IMnpim_wSid_n_Sp_sub_bin_hi[ibin][iq]->Draw("HISTsame");
      //}else if(42 < ibin && ibin<44){
      //  IMnpim_wSid_n_Sp_sub_bin_lo[ibin][iq]->Draw("HISTsame");
      //  IMnpim_wSid_n_Sp_sub_bin_hi[ibin][iq]->Draw("HISTsame");
      //}else if(44 <= ibin){
      //  IMnpim_wSid_n_Sp_sub_bin_hi2[ibin][iq]->Draw("HISTsame");
      //}
      
      /*
      cIMnpim_IMnpip_n_sub_bin_cut[ibin][iq]->cd(2);
      TPaveText *pt3 = new TPaveText(.05,.05,.95,.7);
      double inteSp_g3 = IMnpim_wSid_n_Sp_sub_bin_select[ibin][iq]->Integral();
      double inteSplo_g3 = IMnpim_wSid_n_Sp_sub_bin_lo[ibin][iq]->Integral();
      double inteSplo2_g3 = IMnpim_wSid_n_Sp_sub_bin_lo2[ibin][iq]->Integral();
      double inteSphi_g3 = IMnpim_wSid_n_Sp_sub_bin_hi[ibin][iq]->Integral();
      double inteSphi2_g3 = IMnpim_wSid_n_Sp_sub_bin_hi2[ibin][iq]->Integral();
      double Spestimated_g3 = 0.0;
      double Spoverlap_g3 = 0.0;
      
      if(ibin<42) Spestimated_g3 = inteSplo2_g3;
      else if(ibin==42) Spestimated_g3 = inteSphi_g3*2.0;
      else if(42<ibin && ibin<44) Spestimated_g3 = inteSplo_g3+inteSphi_g3; 
      else Spestimated_g3 = inteSphi2_g3;     
      if(Spestimated_g3<0.0) Spestimated_g3 = 0.0;
      Spoverlap_g3 = inteSp_g3 - Spestimated_g3;
      if(Spoverlap_g3<0.0){
        Spoverlap_g3 = 0.0;
        Spestimated_g3 = 0.0;
      }
      double inteSm_g3 = IMnpip_wSid_n_Sm_sub_bin_select[ibin][iq]->Integral();
      double inteSmlo_g3 = IMnpip_wSid_n_Sm_sub_bin_lo[ibin][iq]->Integral();
      double inteSmlo2_g3 = IMnpip_wSid_n_Sm_sub_bin_lo2[ibin][iq]->Integral();
      double inteSmhi_g3 = IMnpip_wSid_n_Sm_sub_bin_hi[ibin][iq]->Integral();
      double inteSmhi2_g3 = IMnpip_wSid_n_Sm_sub_bin_hi2[ibin][iq]->Integral();
      double Smestimated_g3 = 0.0;
      double Smoverlap_g3 = 0.0;
      
      if(ibin<42) Smestimated_g3 = inteSmlo2_g3;
      else if(ibin==42) Smestimated_g3 = inteSmhi_g3*2.0;
      else if(42< ibin &&  ibin<44) Smestimated_g3 = inteSmhi_g3 + inteSmlo_g3;
      else Smestimated_g3 = inteSmhi2_g3;
      if(Smestimated_g3<0.0) Smestimated_g3 = 0.0;
      Smoverlap_g3 = inteSm_g3 - Smestimated_g3;
      if(Smoverlap_g3<0.0) {
        Smoverlap_g3 = 0.0;
        Smestimated_g3 = 0.0;
      }
      pt3->AddText(Form("IM(n#pi^{-}#pi^{+})  %0.2f-%0.2f %s",1.0+ibin*1.0/nbintemplate,1.0+(ibin+1.0)/nbintemplate,cqcut[iq])); 
      pt3->AddText(Form("Sigma+ count %0.2f ",inteSp_g3));
      if(ibin<42){
        pt3->AddText(Form("Sigma+ side low (4 sigma)   %0.2f ",inteSplo2_g3));
      }else if(ibin==42){
        pt3->AddText(Form("Sigma+ side high*2   %0.2f ",inteSphi_g3*2.0));
      }else if(42 < ibin && ibin<44){
        pt3->AddText(Form("Sigma+ side low   %0.2f ",inteSplo_g3));
        pt3->AddText(Form("Sigma+ side high  %0.2f ",inteSphi_g3));
      }else{
        pt3->AddText(Form("Sigma+ side high (4 sigma)   %0.2f ",inteSphi2_g3));
      } 
      pt3->AddText(Form("Sigma+ estimated %0.2f ", Spestimated_g3));
      pt3->AddText(Form("Sigma- count %0.2f ",inteSm_g3));
      if(ibin <42){
        pt3->AddText(Form("Sigma- side low (4 sigma)   %0.2f ",inteSmlo2_g3));
      }else if(ibin==42){
        pt3->AddText(Form("Sigma- side high*2   %0.2f ",inteSmhi_g3*2.0));
      }else if(42 < ibin && ibin<44){
        pt3->AddText(Form("Sigma- side low   %0.2f ",inteSmlo_g3));
        pt3->AddText(Form("Sigma- side high   %0.2f ",inteSmhi_g3));
      }else{
        pt3->AddText(Form("Sigma- side high (4 sigma)  %0.2f ",inteSmhi2_g3));
      }
      pt3->AddText(Form("Sigma- estimated %0.2f ",Smestimated_g3)); 
      pt3->Draw();

      IMnpipi_overlapdeco_Sm[2][iq]->SetBinContent(ibin,Smestimated_g3);
      IMnpipi_overlapdeco_Sm[2][iq]->SetBinError(ibin,sqrt(Smestimated_g3));
      IMnpipi_overlapdeco_Sp[2][iq]->SetBinContent(ibin,Spestimated_g3);
      IMnpipi_overlapdeco_Sp[2][iq]->SetBinError(ibin,sqrt(Spestimated_g3));
      */
    }
  }
   
  std::cout << __LINE__ << std::endl;
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

  std::cout << __LINE__ << std::endl;
  TCanvas *cratio_SpSm[nqcut];//Sp/Sm  ->group 2
  TCanvas *cratio_SpK0[nqcut];//Sp/K0  ->group 0
  TCanvas *cratio_SmK0[nqcut];//Sm/K0  ->group 1
  TH1D* IMnpipi_SpoverSm_ratio[nqcut];
  TH1D* IMnpipi_SpoverK0_ratio[nqcut];
  TH1D* IMnpipi_SmoverK0_ratio[nqcut];

  for(int iq=qstart;iq<nqcut;iq++){
    IMnpipi_SpoverSm_ratio[iq] = new TH1D(Form("IMnpipi_SpoverSm_ratio_%d",iq),Form("IMnpipi_SpoverSm_ratio_%s",cqcut[iq]),100,1,2);
    IMnpipi_SpoverK0_ratio[iq] = new TH1D(Form("IMnpipi_SpoverK0_ratio_%d",iq),Form("IMnpipi_SpoverK0_ratio_%s",cqcut[iq]),100,1,2);
    IMnpipi_SmoverK0_ratio[iq] = new TH1D(Form("IMnpipi_SmoverK0_ratio_%d",iq),Form("IMnpipi_SmoverK0_ratio_%s",cqcut[iq]),100,1,2);
  
    IMnpipi_SpoverSm_ratio[iq] = (TH1D*)IMnpipi_overlapdeco_Sp[2][iq]->Clone(Form("IMnpipi_SpoverSm_ratio_%d",iq));
    IMnpipi_SpoverSm_ratio[iq]->Divide(IMnpipi_overlapdeco_Sm[2][iq]);
    cratio_SpSm[iq] = new TCanvas(Form("cratrio_SpSm_%s",cqcut[iq]));
    IMnpipi_SpoverSm_ratio[iq]->SetTitle(Form("IMnpipi_Sp/Sm_ratio_%s",cqcut[iq]));
    IMnpipi_SpoverSm_ratio[iq]->Draw("HIST");
    
    IMnpipi_SpoverK0_ratio[iq] = (TH1D*)IMnpipi_overlapdeco_Sp[0][iq]->Clone(Form("IMnpipi_SpoverK0_ratio_%d",iq));
    IMnpipi_SpoverK0_ratio[iq]->Divide(IMnpipi_overlapdeco_K0[0][iq]);
    cratio_SpK0[iq] = new TCanvas(Form("cratrio_SpK0_%s",cqcut[iq]));
    IMnpipi_SpoverK0_ratio[iq]->SetTitle(Form("IMnpipi_Sp/K0_ratio_%s",cqcut[iq]));
    IMnpipi_SpoverK0_ratio[iq]->Draw("HIST");

    IMnpipi_SmoverK0_ratio[iq] = (TH1D*)IMnpipi_overlapdeco_Sm[1][iq]->Clone(Form("IMnpipi_SmoverK0_ratio_%d",iq));
    IMnpipi_SmoverK0_ratio[iq]->Divide(IMnpipi_overlapdeco_K0[1][iq]);
    cratio_SmK0[iq] = new TCanvas(Form("cratrio_SmK0_%s",cqcut[iq]));
    IMnpipi_SmoverK0_ratio[iq]->SetTitle(Form("IMnpipi_Sm/K0_ratio_%s",cqcut[iq]));
    IMnpipi_SmoverK0_ratio[iq]->Draw("HIST");
  }

  
  //merging some specific bins
  TH2D* IMpippim_IMnpip_wK0orwSid_n_merge[10][nqcut];
  TH2D* IMpippim_IMnpim_wK0orwSid_n_merge[10][nqcut];
  TH2D* IMnpim_IMnpip_wSid_n_merge[10][nqcut];
  
  //tuned bins 2d histograms
  TH2F* IMnpim_IMnpip_wSid_n_Sp_wbin[3][nqcut];
  TH2F* IMnpim_IMnpip_wSid_n_Sm_wbin[3][nqcut];
  TH2F* IMpippim_IMnpip_wSid_n_Sp_wbin[3][nqcut];
  TH2F* IMpippim_IMnpip_wK0_n_wbin[3][nqcut];
  TH2F* IMpippim_IMnpim_wSid_n_Sm_wbin[3][nqcut];
  TH2F* IMpippim_IMnpim_wK0_n_wbin[3][nqcut];
  for(int iq=0;iq<nqcut;iq++){
    for(int imerge=0;imerge<3;imerge++){
      IMnpim_IMnpip_wSid_n_Sp_wbin[imerge][iq]= (TH2F*)fr[iq]->Get(Form("IMnpim_IMnpip_wSid_n_Sp_wbin%d",imerge));

      IMnpim_IMnpip_wSid_n_Sm_wbin[imerge][iq]= (TH2F*)fr[iq]->Get(Form("IMnpim_IMnpip_wSid_n_Sm_wbin%d",imerge));

      IMpippim_IMnpip_wSid_n_Sp_wbin[imerge][iq]= (TH2F*)fr[iq]->Get(Form("IMpippim_IMnpip_wSid_n_Sp_wbin%d",imerge));

      IMpippim_IMnpip_wK0_n_wbin[imerge][iq]= (TH2F*)fr[iq]->Get(Form("IMpippim_IMnpip_wK0_n_wbin%d",imerge));

      IMpippim_IMnpim_wSid_n_Sm_wbin[imerge][iq]= (TH2F*)fr[iq]->Get(Form("IMpippim_IMnpim_wSid_n_Sm_wbin%d",imerge));

      IMpippim_IMnpim_wK0_n_wbin[imerge][iq]= (TH2F*)fr[iq]->Get(Form("IMpippim_IMnpim_wK0_n_wbin%d",imerge));
    }
  }

  //projection hists
  TH1D* IMnpip_wK0_merge[10][nqcut];
  TH1D* IMnpip_wK0_merge_select[10][nqcut];
  TH1D* IMnpip_wK0_merge_lohi[10][nqcut];
  TH1D* IMpippim_Sp_merge[10][nqcut];
  TH1D* IMpippim_Sp_merge_select[10][nqcut];
  TH1D* IMpippim_Sp_merge_lohi[10][nqcut];
  TH1D* IMnpim_wK0_merge[10][nqcut];
  TH1D* IMnpim_wK0_merge_select[10][nqcut];
  TH1D* IMnpim_wK0_merge_lo[10][nqcut];
  TH1D* IMnpim_wK0_merge_hi[10][nqcut];
  TH1D* IMnpim_wK0_merge_lohi[10][nqcut];
  TH1D* IMpippim_Sm_merge[10][nqcut];
  TH1D* IMpippim_Sm_merge_select[10][nqcut];
  TH1D* IMpippim_Sm_merge_lo[10][nqcut];
  TH1D* IMpippim_Sm_merge_hi[10][nqcut];
  TH1D* IMpippim_Sm_merge_lohi[10][nqcut];
  TH1D* IMnpip_Sm_merge[10][nqcut];
  TH1D* IMnpip_Sm_merge_select[10][nqcut];
  TH1D* IMnpip_Sm_merge_lo[10][nqcut];
  TH1D* IMnpip_Sm_merge_hi[10][nqcut];
  TH1D* IMnpip_Sm_merge_lohi[10][nqcut];
  TH1D* IMnpim_Sp_merge[10][nqcut];
  TH1D* IMnpim_Sp_merge_select[10][nqcut];
  TH1D* IMnpim_Sp_merge_lo[10][nqcut];
  TH1D* IMnpim_Sp_merge_hi[10][nqcut];
  TH1D* IMnpim_Sp_merge_lohi[10][nqcut];
  
  TCanvas *cIMpippim_IMnpip_merge[10][nqcut];
  TCanvas *cIMpippim_IMnpim_merge[10][nqcut];
  TCanvas *cIMnpim_IMnpip_merge[10][nqcut];
  
  TGraphErrors *gSpwK0[10][nqcut];
  TGraphErrors *gK0wSp[10][nqcut];
  TGraphErrors *gSmwK0[10][nqcut];
  TGraphErrors *gK0wSm[10][nqcut];
  TGraphErrors *gSpwSm[10][nqcut];
  TGraphErrors *gSmwSp[10][nqcut];
  
  TSpline3 *s3SpwK0[10][nqcut];
  TSpline3 *s3K0wSp[10][nqcut];
  TSpline3 *s3SmwK0[10][nqcut];
  TSpline3 *s3K0wSm[10][nqcut];
  TSpline3 *s3SpwSm[10][nqcut];
  TSpline3 *s3SmwSp[10][nqcut];

  TSpline5 *s5SpwK0[10][nqcut];
  TSpline5 *s5K0wSp[10][nqcut];
  TSpline5 *s5SmwK0[10][nqcut];
  TSpline5 *s5K0wSm[10][nqcut];
  TSpline5 *s5SpwSm[10][nqcut];
  TSpline5 *s5SmwSp[10][nqcut];

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
  TH1D* IMnpipi_Sp_ratio_wK0_merge[nqcut];
  TH1D* IMnpipi_Sp_ratio_wSm_merge[nqcut];
  TH1D* IMnpipi_Sm_ratio_wK0_merge[nqcut];
  TH1D* IMnpipi_Sm_ratio_wSp_merge[nqcut];
  TH1D* IMnpipi_K0_ratio_wSp_merge[nqcut];
  TH1D* IMnpipi_K0_ratio_wSm_merge[nqcut];
  for(int iq=0;iq<nqcut;iq++){
    IMnpipi_Sp_ratio_wK0_merge[iq] = new TH1D(Form("IMnpipi_Sp_ratio_wK0_merge_%d",iq),Form("IMnpipi_Sp_ratio_wK0_merge_%s",cqcut[iq]),100,1,2);
    IMnpipi_Sp_ratio_wSm_merge[iq] = new TH1D(Form("IMnpipi_Sp_ratio_wSm_merge_%d",iq),Form("IMnpipi_Sp_ratio_wSm_merge_%s",cqcut[iq]),100,1,2);
    IMnpipi_Sm_ratio_wK0_merge[iq] = new TH1D(Form("IMnpipi_Sm_ratio_wK0_merge_%d",iq),Form("IMnpipi_Sm_ratio_wK0_merge_%s",cqcut[iq]),100,1,2);
    IMnpipi_Sm_ratio_wSp_merge[iq] = new TH1D(Form("IMnpipi_Sm_ratio_wSp_merge_%d",iq),Form("IMnpipi_Sm_ratio_wSp_merge_%s",cqcut[iq]),100,1,2);
    IMnpipi_K0_ratio_wSp_merge[iq] = new TH1D(Form("IMnpipi_K0_ratio_wSp_merge_%d",iq),Form("IMnpipi_K0_ratio_wSp_merge_%s",cqcut[iq]),100,1,2);
    IMnpipi_K0_ratio_wSm_merge[iq] = new TH1D(Form("IMnpipi_K0_ratio_wSm_merge_%d",iq),Form("IMnpipi_K0_ratio_wSm_merge_%s",cqcut[iq]),100,1,2);
  }

  std::cout << __LINE__ << std::endl;
  
  for(int imerge=0;imerge<2;imerge++){
    for(int iq=0;iq<nqcut;iq++){
      IMpippim_IMnpip_wK0orwSid_n_merge[imerge][iq]  
        = (TH2D*)IMpippim_IMnpip_wK0orwSid_n_bin_sub[initbin[imerge]][iq]->Clone(Form("IMpippim_IMnpip_wK0orwSid_merge%d_%d",imerge,iq));
      IMpippim_IMnpip_wK0orwSid_n_merge[imerge][iq]->SetTitle(Form("IMpippim_IMnpip_wK0orwSid %0.2f-%0.2f %s",startval[imerge],endval[imerge],cqcut[iq]));
      IMpippim_IMnpim_wK0orwSid_n_merge[imerge][iq] 
        = (TH2D*)IMpippim_IMnpim_wK0orwSid_n_bin_sub[initbin[imerge]][iq]->Clone(Form("IMpippim_IMnpim_wK0orwSid_merge%d_%d",imerge,iq));
      IMpippim_IMnpim_wK0orwSid_n_merge[imerge][iq]->SetTitle(Form("IMpippim_IMnpim_wK0orwSid %0.2f-%0.2f %s",startval[imerge],endval[imerge],cqcut[iq]));
      IMnpim_IMnpip_wSid_n_merge[imerge][iq]
        = (TH2D*)IMnpim_IMnpip_wK0orwSid_n_bin_sub[initbin[imerge]][iq]->Clone(Form("IMnpim_IMnpip_wSid_n_merge%d_%d",imerge,iq));
      
      IMnpip_wK0_merge[imerge][iq] = (TH1D*)IMpippim_IMnpip_wK0_n_wbin[imerge+1][iq]->ProjectionX(Form("IMnpip_wK0_merge%d_%d",imerge,iq));
      IMpippim_Sp_merge[imerge][iq] = (TH1D*)IMpippim_IMnpip_wSid_n_Sp_wbin[imerge+1][iq]->ProjectionY(Form("IMpippim_Sp_merge%d_%d",imerge,iq));
      IMnpim_wK0_merge[imerge][iq] = (TH1D*)IMpippim_IMnpim_wK0_n_wbin[imerge+1][iq]->ProjectionX(Form("IMnpim_wK0_merge%d_%d",imerge,iq));
      IMpippim_Sm_merge[imerge][iq] = (TH1D*)IMpippim_IMnpim_wSid_n_Sm_wbin[imerge+1][iq]->ProjectionY(Form("IMpippim_Sm_merge%d_%d",imerge,iq));
      IMnpip_Sm_merge[imerge][iq] = (TH1D*)IMnpim_IMnpip_wSid_n_Sm_wbin[imerge+1][iq]->ProjectionX(Form("IMnpip_Sm_merge%d_%d",imerge,iq));
      IMnpim_Sp_merge[imerge][iq] = (TH1D*)IMnpim_IMnpip_wSid_n_Sp_wbin[imerge+1][iq]->ProjectionY(Form("IMnpim_Sp_merge%d_%d",imerge,iq));

      std::cout << __LINE__ << std::endl;
      for(int ibin=startbin[imerge];ibin<endbin[imerge];ibin++){
        IMpippim_IMnpip_wK0orwSid_n_merge[imerge][iq]->Add(IMpippim_IMnpip_wK0orwSid_n_bin_sub[ibin][iq]);
        IMpippim_IMnpim_wK0orwSid_n_merge[imerge][iq]->Add(IMpippim_IMnpim_wK0orwSid_n_bin_sub[ibin][iq]);
        IMnpim_IMnpip_wSid_n_merge[imerge][iq]->Add(IMnpim_IMnpip_wK0orwSid_n_bin_sub[ibin][iq]);
      }
      std::cout << __LINE__ << std::endl;
      //Draw in Canvas 
      cIMpippim_IMnpip_merge[imerge][iq] = new TCanvas(Form("cIMpippim_IMnpip_merge%d_%d",imerge,iq),Form("cIMpippim_IMnpip_merge%d_%d",imerge,iq));
      cIMpippim_IMnpip_merge[imerge][iq]->Divide(2,2);
      cIMpippim_IMnpip_merge[imerge][iq]->cd(3);
      if(RemoveNegative)IMpippim_IMnpip_wK0orwSid_n_merge[imerge][iq]->SetMinimum(0);
      IMpippim_IMnpip_wK0orwSid_n_merge[imerge][iq]->Draw("colz");
      cIMpippim_IMnpip_merge[imerge][iq]->cd(1);
      IMnpip_wK0_merge[imerge][iq]->Draw("HE");
      IMnpip_wK0_merge_select[imerge][iq] = (TH1D*)IMnpip_wK0_merge[imerge][iq]->Clone(Form("IMnpip_wK0_merge_select_%d_%d",imerge,iq));
      IMnpip_wK0_merge_select[imerge][iq]->SetLineColor(2);
      IMnpip_wK0_merge_select[imerge][iq]->GetXaxis()->SetRangeUser(
          IMnpip_wK0_merge_select[imerge][iq]->GetBinLowEdge(9),
          IMnpip_wK0_merge_select[imerge][iq]->GetBinLowEdge(10)
          );
      IMnpip_wK0_merge_select[imerge][iq]->Draw("HEsame");
      IMnpip_wK0_merge_lohi[imerge][iq] = (TH1D*)IMnpip_wK0_merge[imerge][iq]->Clone(Form("IMnpip_wK0_merge_lohi_%d_%d",imerge,iq));
      IMnpip_wK0_merge_lohi[imerge][iq]->SetBinContent(9,0);
      IMnpip_wK0_merge_lohi[imerge][iq]->SetBinError(9,0);
      IMnpip_wK0_merge_lohi[imerge][iq]->SetLineColor(4);
      IMnpip_wK0_merge_lohi[imerge][iq]->Draw("HEsame");
      gSpwK0[imerge][iq] = new TGraphErrors(IMnpip_wK0_merge_lohi[imerge][iq]);
      gSpwK0[imerge][iq]->RemovePoint(8);
      gSpwK0[imerge][iq]->Draw("P");
      s3SpwK0[imerge][iq] = new TSpline3(Form("s3SpwK0%d%d",imerge,iq),gSpwK0[imerge][iq]);
      s3SpwK0[imerge][iq]->SetLineColor(5);
      s3SpwK0[imerge][iq]->Draw("same");
      s5SpwK0[imerge][iq] = new TSpline5(Form("s5SpwK0%d%d",imerge,iq),gSpwK0[imerge][iq]);
      s5SpwK0[imerge][iq]->SetLineColor(6);
      s5SpwK0[imerge][iq]->Draw("same");
      double count_SpwK0 = IMnpip_wK0_merge_select[imerge][iq]->GetBinContent(9);
      double err_SpwK0 = IMnpip_wK0_merge_select[imerge][iq]->GetBinError(9);
      double bincen_SpwK0 = IMnpip_wK0_merge_select[imerge][iq]->GetXaxis()->GetBinCenter(9);
      double bg_SpwK0 =  s3SpwK0[imerge][iq]->Eval(bincen_SpwK0);
      double bg_SpwK0err =  0.0;
      cIMpippim_IMnpip_merge[imerge][iq]->cd(4);
      IMpippim_Sp_merge[imerge][iq]->Draw("HE");
      IMpippim_Sp_merge_select[imerge][iq] = (TH1D*)IMpippim_Sp_merge[imerge][iq]->Clone(Form("IMpippim_Sp_merge_select_%d_%d",imerge,iq));
      IMpippim_Sp_merge_select[imerge][iq]->SetLineColor(2);
      //IMpippim_Sp_merge_lo[imerge][iq]->SetLineColor(4);
      //IMpippim_Sp_merge_hi[imerge][iq]->SetLineColor(4);
      IMpippim_Sp_merge_select[imerge][iq]->GetXaxis()->SetRangeUser(
          IMpippim_Sp_merge_select[imerge][iq]->GetBinLowEdge(10),
          IMpippim_Sp_merge_select[imerge][iq]->GetBinLowEdge(11));
      IMpippim_Sp_merge_select[imerge][iq]->Draw("HEsame");
      IMpippim_Sp_merge_lohi[imerge][iq] = (TH1D*)IMpippim_Sp_merge[imerge][iq]->Clone(Form("IMpippim_Sp_merge_lohi_%d_%d",imerge,iq));
      IMpippim_Sp_merge_lohi[imerge][iq]->SetBinContent(10,0);
      IMpippim_Sp_merge_lohi[imerge][iq]->SetBinError(10,0);
      IMpippim_Sp_merge_lohi[imerge][iq]->SetLineColor(4);
      IMpippim_Sp_merge_lohi[imerge][iq]->Draw("HEsame");
      gK0wSp[imerge][iq] = new TGraphErrors(IMpippim_Sp_merge_lohi[imerge][iq]);
      gK0wSp[imerge][iq]->RemovePoint(9);
      gK0wSp[imerge][iq]->Draw("P");
      s3K0wSp[imerge][iq] = new TSpline3(Form("s3K0wSp%d%d",imerge,iq),gK0wSp[imerge][iq]);
      s3K0wSp[imerge][iq]->SetLineColor(5);
      s3K0wSp[imerge][iq]->Draw("same");
      s5K0wSp[imerge][iq] = new TSpline5(Form("s5K0wSp%d%d",imerge,iq),gK0wSp[imerge][iq]);
      s5K0wSp[imerge][iq]->SetLineColor(6);
      s5K0wSp[imerge][iq]->Draw("same");

      double count_K0wSp = IMpippim_Sp_merge_select[imerge][iq]->GetBinContent(10);
      double err_K0wSp = IMpippim_Sp_merge_select[imerge][iq]->GetBinError(10);
      double bincen_K0wSp = IMpippim_Sp_merge_select[imerge][iq]->GetXaxis()->GetBinCenter(10);
      double bg_K0wSp =  s3K0wSp[imerge][iq]->Eval(bincen_K0wSp);
      double bg_K0wSperr = 0.0;//constfitK0wSp[imerge][iq]->GetParError(0);

      cIMpippim_IMnpip_merge[imerge][iq]->cd(2);
      TPaveText *pt = new TPaveText(.05,.05,.95,.7);
      //pt->AddText(Form("Sigma+ on K0 count %0.2f",count_SpwK0)); 
      pt->AddText(Form("Sigma+ & K0 overlap %0.1f +/- %0.1f",count_SpwK0,err_K0wSp)); 
      double errSpwK0 = sqrt(err_SpwK0*err_SpwK0+bg_SpwK0err*bg_SpwK0err);
      pt->AddText(Form("Sigma+ estimated %0.1f +/- %0.1f",count_SpwK0 - bg_SpwK0,errSpwK0 )); 
      //pt->AddText(Form("K0 on Sigma+ count %0.2f",count_K0wSp)); 
      double errK0wSp = sqrt(err_K0wSp*err_K0wSp+bg_K0wSperr*bg_K0wSperr);
      pt->AddText(Form("K0 estimated %0.1f +/- %0.1f",count_SpwK0 - bg_K0wSp,errK0wSp )); 
      pt->AddText(Form("Sigma+/K0 ratio %0.1f +/- %0.1f", (count_SpwK0-bg_SpwK0)/(count_K0wSp-bg_K0wSp),
          sqrt(errSpwK0*errSpwK0/pow(count_K0wSp-bg_K0wSp,2)+pow(count_SpwK0-bg_SpwK0 ,2)/pow(count_K0wSp-bg_K0wSp,4)*pow(errK0wSp,2))));
      pt->Draw();
      double Sp_ratio_wK0 = (count_SpwK0-bg_SpwK0)/((count_SpwK0-bg_SpwK0)+(count_K0wSp-bg_K0wSp));
      double K0_ratio_wSp = (count_K0wSp-bg_K0wSp)/((count_SpwK0-bg_SpwK0)+(count_K0wSp-bg_K0wSp));
      for(int ibin=initbin[imerge];ibin<endbin[imerge];ibin++){
        IMnpipi_Sp_ratio_wK0_merge[iq]->SetBinContent(ibin,Sp_ratio_wK0);
        IMnpipi_Sp_ratio_wK0_merge[iq]->SetBinError(ibin,0);
        IMnpipi_K0_ratio_wSp_merge[iq]->SetBinContent(ibin,K0_ratio_wSp);
        IMnpipi_K0_ratio_wSp_merge[iq]->SetBinError(ibin,0);
      }

      std::cout << __LINE__ << std::endl;

      //Sigma- and K0 overlap 
      cIMpippim_IMnpim_merge[imerge][iq] = new TCanvas(Form("cIMpippim_IMnpim_merge%d_%d",imerge,iq),Form("cIMpippim_IMnpim_merge%d_%d",imerge,iq));
      cIMpippim_IMnpim_merge[imerge][iq]->Divide(2,2);
      cIMpippim_IMnpim_merge[imerge][iq]->cd(3);
      if(RemoveNegative)IMpippim_IMnpim_wK0orwSid_n_merge[imerge][iq]->SetMinimum(0);
      IMpippim_IMnpim_wK0orwSid_n_merge[imerge][iq]->Draw("colz");
      cIMpippim_IMnpim_merge[imerge][iq]->cd(1);
      IMnpim_wK0_merge[imerge][iq]->Draw("HE");
      IMnpim_wK0_merge_select[imerge][iq] = (TH1D*)IMnpim_wK0_merge[imerge][iq]->Clone(Form("IMnpim_wK0_merge_select_%d_%d",imerge,iq));
      IMnpim_wK0_merge_select[imerge][iq]->SetLineColor(2);
      IMnpim_wK0_merge_select[imerge][iq]->GetXaxis()->SetRangeUser(
          IMnpim_wK0_merge_select[imerge][iq]->GetBinLowEdge(9),
          IMnpim_wK0_merge_select[imerge][iq]->GetBinLowEdge(10)); 

      IMnpim_wK0_merge_select[imerge][iq]->Draw("HEsame");
      IMnpim_wK0_merge_lohi[imerge][iq] = (TH1D*)IMnpim_wK0_merge[imerge][iq]->Clone(Form("IMnpim_wK0_merge_lohi_%d_%d",imerge,iq));
      IMnpim_wK0_merge_lohi[imerge][iq]->SetBinContent(9,0);
      IMnpim_wK0_merge_lohi[imerge][iq]->SetBinError(9,0);
      IMnpim_wK0_merge_lohi[imerge][iq]->SetLineColor(4);
      IMnpim_wK0_merge_lohi[imerge][iq]->Draw("HEsame");
      gSmwK0[imerge][iq] = new TGraphErrors(IMnpim_wK0_merge_lohi[imerge][iq]);
      gSmwK0[imerge][iq]->RemovePoint(8);
      gSmwK0[imerge][iq]->Draw("P");
      s3SmwK0[imerge][iq] = new TSpline3(Form("s3SmwK0%d%d",imerge,iq),gSmwK0[imerge][iq]);
      s3SmwK0[imerge][iq]->SetLineColor(5);
      s3SmwK0[imerge][iq]->Draw("same");
      s5SmwK0[imerge][iq] = new TSpline5(Form("s5SmwK0%d%d",imerge,iq),gSmwK0[imerge][iq]);
      s5SmwK0[imerge][iq]->SetLineColor(6);
      s5SmwK0[imerge][iq]->Draw("same");

      double count_SmwK0 = IMnpim_wK0_merge_select[imerge][iq]->GetBinContent(9);
      double err_SmwK0 = IMnpim_wK0_merge_select[imerge][iq]->GetBinError(9);
      double bincen_SmwK0 = IMnpim_wK0_merge_select[imerge][iq]->GetXaxis()->GetBinCenter(9);
      double bg_SmwK0 =  s3SmwK0[imerge][iq]->Eval(bincen_SmwK0);
      double bg_SmwK0err =  0.0;//constfitSmwK0[imerge][iq]->GetParError(0);
      
      
      cIMpippim_IMnpim_merge[imerge][iq]->cd(4);
      IMpippim_Sm_merge[imerge][iq]->Draw("HE");
      IMpippim_Sm_merge_select[imerge][iq] = (TH1D*)IMpippim_Sm_merge[imerge][iq]->Clone(Form("IMpippim_Sm_merge_select_%d_%d",imerge,iq));
      IMpippim_Sm_merge_select[imerge][iq]->SetLineColor(2);
      IMpippim_Sm_merge_select[imerge][iq]->GetXaxis()->SetRangeUser(
          IMpippim_Sm_merge_select[imerge][iq]->GetBinLowEdge(10),
          IMpippim_Sm_merge_select[imerge][iq]->GetBinLowEdge(11));
      IMpippim_Sm_merge_select[imerge][iq]->Draw("HEsame");
      IMpippim_Sm_merge_lohi[imerge][iq] = (TH1D*)IMpippim_Sm_merge[imerge][iq]->Clone(Form("IMpippim_Sm_merge_lohi_%d_%d",imerge,iq));
      IMpippim_Sm_merge_lohi[imerge][iq]->SetLineColor(4);
      IMpippim_Sm_merge_lohi[imerge][iq]->SetBinContent(10,0);
      IMpippim_Sm_merge_lohi[imerge][iq]->SetBinError(10,0);
      IMpippim_Sm_merge_lohi[imerge][iq]->Draw("HEsame");
      gK0wSm[imerge][iq] = new TGraphErrors(IMpippim_Sm_merge_lohi[imerge][iq]);
      gK0wSm[imerge][iq]->RemovePoint(9);
      gK0wSm[imerge][iq]->Draw("P");
      s3K0wSm[imerge][iq] = new TSpline3(Form("s3K0wSm%d%d",imerge,iq),gK0wSm[imerge][iq]);
      s3K0wSm[imerge][iq]->SetLineColor(5);
      s3K0wSm[imerge][iq]->Draw("same");
      s5K0wSm[imerge][iq] = new TSpline5(Form("s5K0wSm%d%d",imerge,iq),gK0wSm[imerge][iq]);
      s5K0wSm[imerge][iq]->SetLineColor(6);
      s5K0wSm[imerge][iq]->Draw("same");
      double count_K0wSm = IMpippim_Sm_merge_select[imerge][iq]->GetBinContent(10);
      double err_K0wSm = IMpippim_Sm_merge_select[imerge][iq]->GetBinError(10);
      double bincen_K0wSm = IMpippim_Sm_merge_select[imerge][iq]->GetXaxis()->GetBinCenter(10);
      double bg_K0wSm =  s3K0wSm[imerge][iq]->Eval(bincen_K0wSm);
      double bg_K0wSmerr = 0.0;//  constfitK0wSm[imerge][iq]->GetParError(0);
      
      cIMpippim_IMnpim_merge[imerge][iq]->cd(2);
      TPaveText *pt1 = new TPaveText(.05,.05,.95,.7);
      pt1->AddText(Form("Sigma- & K0 overlap %0.1f +/- %0.1f",count_SmwK0,err_SmwK0)); 
      double errSmwK0 = sqrt(err_SmwK0*err_SmwK0 + bg_SmwK0*bg_SmwK0);
      //pt1->AddText(Form("Sigma- es %0.1f",bg_SmwK0)); 
      pt1->AddText(Form("Sigma- estimated %0.1f +/- %0.1f",count_SmwK0-bg_SmwK0,errSmwK0)); 
      //pt1->AddText(Form("K0 on Sigma- count %0.1f",count_K0wSm)); 
      double errK0wSm = sqrt(err_K0wSm*err_K0wSm + bg_K0wSm*bg_K0wSm);
      pt1->AddText(Form("K0 estimated %0.1f",count_K0wSm-bg_K0wSm,errK0wSm)); 
      //pt1->AddText(Form("Sigma-/K0 ratio %0.2f", (count_SmwK0-bg_SmwK0)/(count_K0wSm-bg_K0wSm)));
      pt1->AddText(Form("Sigma-/K0 ratio %0.1f +/- %0.1f", (count_SmwK0-bg_SmwK0)/(count_K0wSm-bg_K0wSm),
          sqrt(errSmwK0*errSmwK0/pow(count_K0wSm-bg_K0wSm,2)+pow(count_SmwK0-bg_SmwK0 ,2)/pow(count_K0wSm-bg_K0wSm,4)*pow(errK0wSm,2))));
      pt1->Draw();
      
      double Sm_ratio_wK0 = (count_SmwK0-bg_SmwK0)/((count_SmwK0-bg_SmwK0)+(count_K0wSm-bg_K0wSm));
      double K0_ratio_wSm = (count_K0wSm-bg_K0wSm)/((count_SmwK0-bg_SmwK0)+(count_K0wSm-bg_K0wSm));
      for(int ibin=initbin[imerge];ibin<endbin[imerge];ibin++){
        IMnpipi_Sm_ratio_wK0_merge[iq]->SetBinContent(ibin,Sm_ratio_wK0);
        IMnpipi_Sm_ratio_wK0_merge[iq]->SetBinError(ibin,0);
        IMnpipi_K0_ratio_wSm_merge[iq]->SetBinContent(ibin,K0_ratio_wSm);
        IMnpipi_K0_ratio_wSm_merge[iq]->SetBinError(ibin,0);
      }

      cIMnpim_IMnpip_merge[imerge][iq] = new TCanvas(Form("cIMnpim_IMnpip_merge%d_%d",imerge,iq),Form("cIMnpim_IMnpip_merge%d_%d",imerge,iq));
      cIMnpim_IMnpip_merge[imerge][iq]->Divide(2,2);
      cIMnpim_IMnpip_merge[imerge][iq]->cd(3);
      if(RemoveNegative)IMnpim_IMnpip_wSid_n_merge[imerge][iq]->SetMinimum(0);
      IMnpim_IMnpip_wSid_n_merge[imerge][iq]->GetXaxis()->SetRangeUser(1,1.4);
      IMnpim_IMnpip_wSid_n_merge[imerge][iq]->GetYaxis()->SetRangeUser(1,1.4);
      IMnpim_IMnpip_wSid_n_merge[imerge][iq]->Draw("colz");
      cIMnpim_IMnpip_merge[imerge][iq]->cd(1);
      IMnpip_Sm_merge[imerge][iq]->GetXaxis()->SetRangeUser(1,1.4);
      IMnpip_Sm_merge[imerge][iq]->SetMinimum(0);
      IMnpip_Sm_merge[imerge][iq]->Draw("HE");
      IMnpip_Sm_merge_select[imerge][iq] = (TH1D*)IMnpip_Sm_merge[imerge][iq]->Clone(Form("IMnpip_Sm_merge_select_%d_%d",imerge,iq));
      IMnpip_Sm_merge_select[imerge][iq]->SetLineColor(2);
      IMnpip_Sm_merge_select[imerge][iq]->GetXaxis()->SetRangeUser(
          IMnpip_Sm_merge_select[imerge][iq]->GetBinLowEdge(9),
          IMnpip_Sm_merge_select[imerge][iq]->GetBinLowEdge(10));
      
      IMnpip_Sm_merge_select[imerge][iq]->Draw("HEsame");
      IMnpip_Sm_merge_lohi[imerge][iq] = (TH1D*)IMnpip_Sm_merge[imerge][iq]->Clone(Form("IMnpip_Sm_merge_lohi_%d_%d",imerge,iq));
      IMnpip_Sm_merge_lohi[imerge][iq]->SetLineColor(4);
      IMnpip_Sm_merge_lohi[imerge][iq]->SetBinContent(9,0);
      IMnpip_Sm_merge_lohi[imerge][iq]->SetBinError(9,0);
      IMnpip_Sm_merge_lohi[imerge][iq]->Draw("HEsame");
      gSpwSm[imerge][iq] = new TGraphErrors(IMnpip_Sm_merge_lohi[imerge][iq]);
      gSpwSm[imerge][iq]->RemovePoint(8);
      gSpwSm[imerge][iq]->Draw("P");
      s3SpwSm[imerge][iq] = new TSpline3(Form("s3SpwSm%d%d",imerge,iq),gSpwSm[imerge][iq]);
      s3SpwSm[imerge][iq]->SetLineColor(5);
      s3SpwSm[imerge][iq]->Draw("same");
      s5SpwSm[imerge][iq] = new TSpline5(Form("s5SpwSm%d%d",imerge,iq),gSpwSm[imerge][iq]);
      s5SpwSm[imerge][iq]->SetLineColor(6);
      s5SpwSm[imerge][iq]->Draw("same");
      double count_SpwSm = IMnpip_Sm_merge_select[imerge][iq]->GetBinContent(9);
      double err_SpwSm = IMnpip_Sm_merge_select[imerge][iq]->GetBinError(9);
      double bincen_SpwSm = IMnpip_Sm_merge_select[imerge][iq]->GetXaxis()->GetBinCenter(9);
      double bg_SpwSm =  s3SpwSm[imerge][iq]->Eval(bincen_SpwSm);
      double bg_SpwSmerr = 0.0;// constfitSpwSm[imerge][iq]->GetParError(0);

      cIMnpim_IMnpip_merge[imerge][iq]->cd(4);
      IMnpim_Sp_merge[imerge][iq]->GetXaxis()->SetRangeUser(1,1.4);
      IMnpim_Sp_merge[imerge][iq]->SetMinimum(0);
      IMnpim_Sp_merge[imerge][iq]->Draw("HE");
      IMnpim_Sp_merge_select[imerge][iq] = (TH1D*)IMnpim_Sp_merge[imerge][iq]->Clone(Form("IMnpim_Sp_merge_select_%d_%d",imerge,iq));
      IMnpim_Sp_merge_select[imerge][iq]->SetLineColor(2);
      IMnpim_Sp_merge_select[imerge][iq]->GetXaxis()->SetRangeUser(
          IMnpim_Sp_merge_select[imerge][iq]->GetBinLowEdge(9),
          IMnpim_Sp_merge_select[imerge][iq]->GetBinLowEdge(10));
      
      IMnpim_Sp_merge_select[imerge][iq]->Draw("HEsame");
      IMnpim_Sp_merge_lohi[imerge][iq] = (TH1D*)IMnpim_Sp_merge[imerge][iq]->Clone(Form("IMnpim_Sp_merge_lohi_%d_%d",imerge,iq));
      IMnpim_Sp_merge_lohi[imerge][iq]->SetLineColor(4);
      IMnpim_Sp_merge_lohi[imerge][iq]->SetBinContent(9,0);
      IMnpim_Sp_merge_lohi[imerge][iq]->SetBinError(9,0);
      IMnpim_Sp_merge_lohi[imerge][iq]->Draw("HEsame");
      gSmwSp[imerge][iq] = new TGraphErrors(IMnpim_Sp_merge_lohi[imerge][iq]);
      gSmwSp[imerge][iq]->RemovePoint(8);
      gSmwSp[imerge][iq]->Draw("P");
      s3SmwSp[imerge][iq] = new TSpline3(Form("s3SmwSp%d%d",imerge,iq),gSmwSp[imerge][iq]);
      s3SmwSp[imerge][iq]->SetLineColor(5);
      s3SmwSp[imerge][iq]->Draw("same");
      s5SmwSp[imerge][iq] = new TSpline5(Form("s5SmwSp%d%d",imerge,iq),gSmwSp[imerge][iq]);
      s5SmwSp[imerge][iq]->SetLineColor(6);
      s5SmwSp[imerge][iq]->Draw("same");
      double count_SmwSp = IMnpim_Sp_merge_select[imerge][iq]->GetBinContent(9);
      double err_SmwSp = IMnpim_Sp_merge_select[imerge][iq]->GetBinError(9);
      double bincen_SmwSp = IMnpim_Sp_merge_select[imerge][iq]->GetXaxis()->GetBinCenter(9);
      double bg_SmwSp =  s3SmwSp[imerge][iq]->Eval(bincen_SmwSp);
      double bg_SmwSperr = 0.0;// constfitSmwSp[imerge][iq]->GetParError(0);
      cIMnpim_IMnpip_merge[imerge][iq]->cd(2);
      TPaveText *pt2 = new TPaveText(.05,.05,.95,.7);
      pt2->AddText(Form("Sigma+ & Sigma- overlap %0.1f +/- %0.1f",count_SpwSm,err_SpwSm)); 
      double errSpwSm = sqrt(err_SmwSp*err_SmwSp+bg_SpwSmerr*bg_SpwSmerr);
      pt2->AddText(Form("Sigma+ estimated %0.1f +/- %0.1f",count_SpwSm-bg_SpwSm,errSpwSm)); 
      //pt2->AddText(Form("Simga- on Sigma+ count %0.2f",count_SmwSp)); 
      double errSmwSp = sqrt(err_SpwSm*err_SpwSm+bg_SmwSperr*bg_SmwSperr);
      pt2->AddText(Form("Sigma- estimated %0.1f +/- %0.1f",count_SmwSp-bg_SmwSp,errSmwSp)); 
      //pt2->AddText(Form("Sigma+/Sigma- ratio %0.1f +/- %0.1f", (count_SpwSm-bg_SpwSm)/(count_SmwSp-bg_SmwSp)));
      pt2->AddText(Form("Sigma+/Sigma- ratio  %0.1f +/- %0.1f", (count_SpwSm-bg_SpwSm)/(count_SmwSp-bg_SmwSp),
          sqrt(errSpwSm*errSpwSm/pow(count_SmwSp-bg_SmwSp,2)+pow(count_SpwSm-bg_SpwSm,2)/pow(count_SmwSp-bg_SmwSp,4)*pow(errSmwSp,2))));
      pt2->Draw();

      double Sp_ratio_wSm = (count_SpwSm-bg_SpwSm)/((count_SpwSm-bg_SpwSm)+(count_SmwSp-bg_SmwSp));
      double Sm_ratio_wSp = (count_SmwSp-bg_SmwSp)/((count_SmwSp-bg_SmwSp)+(count_SpwSm-bg_SpwSm));
      if(imerge==1){
        Sp_ratio_wSm = 0.0;
        Sm_ratio_wSp = 0.0;
      }
      for(int ibin=initbin[imerge];ibin<endbin[imerge];ibin++){
        IMnpipi_Sp_ratio_wSm_merge[iq]->SetBinContent(ibin,Sp_ratio_wSm);
        IMnpipi_Sp_ratio_wSm_merge[iq]->SetBinError(ibin,0);
        IMnpipi_Sm_ratio_wSp_merge[iq]->SetBinContent(ibin,Sm_ratio_wSp);
        IMnpipi_Sm_ratio_wSp_merge[iq]->SetBinError(ibin,0);
      }
    }
  }
  
  TCanvas *cIMnpipi_Sp_ratio_wK0_merge[nqcut];
  TCanvas *cIMnpipi_K0_ratio_wSp_merge[nqcut];
  TCanvas *cIMnpipi_Sm_ratio_wK0_merge[nqcut];
  TCanvas *cIMnpipi_K0_ratio_wSm_merge[nqcut];
  TCanvas *cIMnpipi_Sp_ratio_wSm_merge[nqcut];
  TCanvas *cIMnpipi_Sm_ratio_wSp_merge[nqcut];
  
  for(int iq=0;iq<nqcut;iq++){
    cIMnpipi_Sp_ratio_wK0_merge[iq] = new TCanvas(Form("cIMnpipi_Sp_ratio_wK0_merge_%d",iq),Form("cIMnpipi_Sp_ratio_wK0_merge_%d",iq));
    IMnpipi_Sp_ratio_wK0_merge[iq]->Draw("HIST");

    cIMnpipi_K0_ratio_wSp_merge[iq] = new TCanvas(Form("cIMnpipi_K0_ratio_wSp_merge_%d",iq),Form("cIMnpipi_K0_ratio_wSp_merge_%d",iq));
    IMnpipi_K0_ratio_wSp_merge[iq]->Draw("HIST");
    
    cIMnpipi_Sm_ratio_wK0_merge[iq] = new TCanvas(Form("cIMnpipi_Sm_ratio_wK0_merge_%d",iq),Form("cIMnpipi_Sm_ratio_wK0_merge_%d",iq));
    IMnpipi_Sm_ratio_wK0_merge[iq]->Draw("HIST");

    cIMnpipi_K0_ratio_wSm_merge[iq] = new TCanvas(Form("cIMnpipi_K0_ratio_wSm_merge_%d",iq),Form("cIMnpipi_K0_ratio_wSm_merge_%d",iq));
    IMnpipi_K0_ratio_wSm_merge[iq]->Draw("HIST");

    cIMnpipi_Sp_ratio_wSm_merge[iq] = new TCanvas(Form("cIMnpipi_Sp_ratio_wSm_merge_%d",iq),Form("cIMnpipi_Sp_ratio_wSm_merge_%d",iq));
    IMnpipi_Sp_ratio_wSm_merge[iq]->Draw("HIST");

    cIMnpipi_Sm_ratio_wSp_merge[iq] = new TCanvas(Form("cIMnpipi_Sm_ratio_wSp_merge_%d",iq),Form("cIMnpipi_Sm_ratio_wSp_merge_%d",iq));
    IMnpipi_Sm_ratio_wSp_merge[iq]->Draw("HIST");
  }


  TCanvas *c = NULL;
  TSeqCollection *SCol = gROOT->GetListOfCanvases();
  int size = SCol->GetSize();
  TIter next(SCol);
  TString pdfname = "K0SigmaDeco.pdf";
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
  
  TIter nexthist2(gDirectory->GetList());
  TFile *fout;
  fout = new TFile("K0SigmaDeco.root","RECREATE");
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
