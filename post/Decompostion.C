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


void Decomposition()
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

  //correlation plots including overlap region
  TH2F* IMpippim_IMnpip_wSid_n_Sp[nqcut];
  TH2F* IMpippim_IMnpip_wSid_n_Sp_bin[nbintemplate][nqcut];
  TH2F* IMpippim_IMnpim_wSid_n_Sm[nqcut];
  TH2F* IMpippim_IMnpim_wSid_n_Sm_bin[nbintemplate][nqcut];
  TH2F* IMpippim_IMnpip_wK0_n[nqcut];
  TH2F* IMpippim_IMnpip_wK0_n_bin[nbintemplate][nqcut];
  TH2F* IMpippim_IMnpim_wK0_n[nqcut];
  TH2F* IMpippim_IMnpim_wK0_n_bin[nbintemplate][nqcut];
  TH2F* IMpippim_IMnpip_wK0orwSid_n[nqcut];
  TH2F* IMpippim_IMnpim_wK0orwSid_n[nqcut];
  TH2F* IMpippim_IMnpip_wK0orwSid_n_bin[nbintemplate][nqcut];
  TH2F* IMpippim_IMnpim_wK0orwSid_n_bin[nbintemplate][nqcut];
  TH2F* IMnpim_IMnpip_wK0orwSid_n_bin[nbintemplate][nqcut];
  TH2F* IMnpim_IMnpip_wSid_n_Sp_bin[nbintemplate][nqcut];
  TH2F* IMnpim_IMnpip_wSid_n_Sm_bin[nbintemplate][nqcut];
  
  


  //for the overlap of S+ & S- & K0 counting 
  TH2F* q_IMnpipi_wK0_wSid_n_SpSm[nqcut];
  TH1D* IMnpipi_wK0_wSid_n_SpSm[nqcut];
  TCanvas *cq_IMnpipi_wK0_wSid_n_SpSm[nqcut];
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
    IMpippim_IMnpip_wK0orwSid_n_sub[iq]->RebinX(4);
    IMpippim_IMnpip_wK0orwSid_n_sub[iq]->RebinY(5);
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
    IMpippim_IMnpim_wK0orwSid_n_sub[iq]->RebinX(4);
    IMpippim_IMnpim_wK0orwSid_n_sub[iq]->RebinY(5);
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
    IMnpim_IMnpip_wK0orwSid_n_sub[iq]->RebinX(4);
    IMnpim_IMnpip_wK0orwSid_n_sub[iq]->RebinY(4);
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
      
      IMpippim_IMnpip_wK0orwSid_n_bin_sub[ibin][iq]->RebinX(4);
      IMpippim_IMnpip_wK0orwSid_n_bin_sub[ibin][iq]->RebinY(5);
      IMpippim_IMnpip_wK0orwSid_n_bin_sub[ibin][iq]->GetXaxis()->SetRangeUser(1,1.8);
      IMpippim_IMnpip_wK0orwSid_n_bin_sub[ibin][iq]->GetYaxis()->SetRangeUser(0.2,0.9);
      IMpippim_IMnpip_wK0orwSid_n_bin_sub[ibin][iq]->Draw("colz");
      cIMpippim_IMnpip_n_sub_bin_cut[ibin][iq]->cd(1);
      IMnpip_wK0_n_sub_bin[ibin][iq] = (TH1D*)IMpippim_IMnpip_wK0_n_bin_sub[ibin][iq]->ProjectionX(Form("IMnpip_wK0_n_sub_bin%d_%d",ibin,iq));
      IMnpip_wK0_n_sub_bin[ibin][iq]->SetTitle(Form("IMnpip_wK0_n_sub_bin%d_%d",ibin,iq));
      if(RebinMode)IMnpip_wK0_n_sub_bin[ibin][iq]->RebinX(10);
      IMnpip_wK0_n_sub_bin[ibin][iq]->GetXaxis()->SetRangeUser(1,1.8);
      IMnpip_wK0_n_sub_bin[ibin][iq]->Draw("HIST");
      IMnpip_wK0_n_sub_bin_select[ibin][iq] = (TH1D*)IMnpip_wK0_n_sub_bin[ibin][iq]->Clone(Form("IMnpip_wK0_n_sub_bin%d_%d_select",ibin,iq));
      IMnpip_wK0_n_sub_bin_select[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmap_MIN,anacuts::Sigmap_MAX);
      IMnpip_wK0_n_sub_bin_select[ibin][iq]->SetLineColor(2);
      IMnpip_wK0_n_sub_bin_select[ibin][iq]->Draw("HISTsame");
      cIMpippim_IMnpip_n_sub_bin_cut[ibin][iq]->cd(4);
      IMpippim_wSid_n_Sp_sub_bin[ibin][iq]
      = (TH1D*)IMpippim_IMnpip_wSid_n_Sp_bin_sub[ibin][iq]->ProjectionY(Form("IMpippim_wSid_n_Sp_sub_bin%d_%d",ibin,iq));
      IMpippim_wSid_n_Sp_sub_bin[ibin][iq]->SetTitle(Form("IMpippim_wSid_n_Sp_sub_bin%d_%d",ibin,iq));
      if(RebinMode)IMpippim_wSid_n_Sp_sub_bin[ibin][iq]->RebinX(30);
      IMpippim_wSid_n_Sp_sub_bin[ibin][iq]->GetXaxis()->SetRangeUser(0.2,0.9);
      IMpippim_wSid_n_Sp_sub_bin[ibin][iq]->Draw("HIST");
      IMpippim_wSid_n_Sp_sub_bin_select[ibin][iq] = (TH1D*)IMpippim_wSid_n_Sp_sub_bin[ibin][iq]->Clone(Form("IMpippim_wSid_n_Sp_sub_bin_%d_%d_select",ibin,iq));
      IMpippim_wSid_n_Sp_sub_bin_select[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::pipi_MIN_narrow,anacuts::pipi_MAX_narrow);
      IMpippim_wSid_n_Sp_sub_bin_select[ibin][iq]->SetLineColor(2);
      IMpippim_wSid_n_Sp_sub_bin_select[ibin][iq]->Draw("HISTsame");
       
      //2nd canvas to display the calculation of decomposion of K0, Sigma+,Sigma- 
      //2D plot in canvas(3) has cut to select K0 or Sigma-
      //projection plot in canvas(1) select Sigma- events candidate + (Sigma- & K0 ) overlap events
      //projection plot in canvas(4) select K0     events candidate + (Sigma- & K0 ) overlap events 
      //std::cout << __LINE__ << std::endl;
      cIMpippim_IMnpim_n_sub_bin_cut[ibin][iq]
      = new TCanvas(Form("cIMpippim_IMnpim_n_sub_bin_cut_%d_%d",ibin,iq),Form("cIMpippim_IMnpim_n_sub_bin_cut%d_%d",ibin,iq),800,800);
      cIMpippim_IMnpim_n_sub_bin_cut[ibin][iq]->Divide(2,2);
      cIMpippim_IMnpim_n_sub_bin_cut[ibin][iq]->cd(3);

      IMpippim_IMnpim_wK0orwSid_n_bin_sub[ibin][iq]->RebinX(4);
      IMpippim_IMnpim_wK0orwSid_n_bin_sub[ibin][iq]->RebinY(5);
      IMpippim_IMnpim_wK0orwSid_n_bin_sub[ibin][iq]->GetXaxis()->SetRangeUser(1,1.8);
      IMpippim_IMnpim_wK0orwSid_n_bin_sub[ibin][iq]->GetYaxis()->SetRangeUser(0.2,0.9);
      IMpippim_IMnpim_wK0orwSid_n_bin_sub[ibin][iq]->Draw("colz");
      cIMpippim_IMnpim_n_sub_bin_cut[ibin][iq]->cd(1);
      IMnpim_wK0_n_sub_bin[ibin][iq] = (TH1D*)IMpippim_IMnpim_wK0_n_bin_sub[ibin][iq]->ProjectionX(Form("IMnpim_wK0_n_sub_bin%d_%d",ibin,iq));
      IMnpim_wK0_n_sub_bin[ibin][iq]->SetTitle(Form("IMnpim_wK0_n_sub_bin%d_%d",ibin,iq));
      if(RebinMode) IMnpim_wK0_n_sub_bin[ibin][iq]->RebinX(10);
      IMnpim_wK0_n_sub_bin[ibin][iq]->GetXaxis()->SetRangeUser(1,1.8);
      IMnpim_wK0_n_sub_bin[ibin][iq]->Draw("HIST");
      IMnpim_wK0_n_sub_bin_select[ibin][iq] = (TH1D*)IMnpim_wK0_n_sub_bin[ibin][iq]->Clone(Form("IMnpim_wK0_n_sub_bin%d_%d_select",ibin,iq));
      IMnpim_wK0_n_sub_bin_select[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmam_MIN,anacuts::Sigmam_MAX);
      IMnpim_wK0_n_sub_bin_select[ibin][iq]->SetLineColor(2);
      IMnpim_wK0_n_sub_bin_select[ibin][iq]->Draw("HISTsame");
      cIMpippim_IMnpim_n_sub_bin_cut[ibin][iq]->cd(4);
      IMpippim_wSid_n_Sm_sub_bin[ibin][iq]
      = (TH1D*)IMpippim_IMnpim_wSid_n_Sm_bin_sub[ibin][iq]->ProjectionY(Form("IMpippim_wSid_n_Sm_sub_bin%d_%d",ibin,iq));
      IMpippim_wSid_n_Sm_sub_bin[ibin][iq]->SetTitle(Form("IMpippim_wSid_n_Sm_sub_bin%d_%d",ibin,iq));
      if(RebinMode)IMpippim_wSid_n_Sm_sub_bin[ibin][iq]->RebinX(30);
      IMpippim_wSid_n_Sm_sub_bin[ibin][iq]->Draw("HIST");
      IMpippim_wSid_n_Sm_sub_bin_select[ibin][iq] = (TH1D*)IMpippim_wSid_n_Sm_sub_bin[ibin][iq]->Clone(Form("IMpippim_wSid_n_Sm_sub_bin_%d_%d_select",ibin,iq));
      IMpippim_wSid_n_Sm_sub_bin_select[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::pipi_MIN_narrow,anacuts::pipi_MAX_narrow);
      IMpippim_wSid_n_Sm_sub_bin_select[ibin][iq]->SetLineColor(2);
      IMpippim_wSid_n_Sm_sub_bin_select[ibin][iq]->Draw("HISTsame");

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
      if(RebinMode){
        IMnpip_wSid_n_Sm_sub_bin[ibin][iq]->RebinX(10);
      }
      IMnpip_wSid_n_Sm_sub_bin[ibin][iq]->GetXaxis()->SetRangeUser(1,1.8);
      IMnpip_wSid_n_Sm_sub_bin[ibin][iq]->Draw("HIST");

      IMnpip_wSid_n_Sm_sub_bin_select[ibin][iq] = (TH1D*)IMnpip_wSid_n_Sm_sub_bin[ibin][iq]->Clone(Form("IMnpip_wSid_n_Sm_sub_bin%d_%d_select",ibin,iq));
      IMnpip_wSid_n_Sm_sub_bin_select[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmap_MIN,anacuts::Sigmap_MAX);
      IMnpip_wSid_n_Sm_sub_bin_select[ibin][iq]->SetLineColor(2);
      IMnpip_wSid_n_Sm_sub_bin_select[ibin][iq]->Draw("HISTsame");
      //std::cout << __LINE__ << std::endl;
      cIMnpim_IMnpip_n_sub_bin_cut[ibin][iq]->cd(4);
      IMnpim_wSid_n_Sp_sub_bin[ibin][iq]
      = (TH1D*)IMnpim_IMnpip_wSid_n_Sp_bin_sub[ibin][iq]->ProjectionY(Form("IMnpim_wSid_n_Sp_sub_bin%d_%d",ibin,iq));
      IMnpim_wSid_n_Sp_sub_bin[ibin][iq]->SetTitle(Form("IMnpim_wSid_n_Sp_sub_bin%d_%d",ibin,iq));
      if(RebinMode)IMnpim_wSid_n_Sp_sub_bin[ibin][iq]->RebinX(10);
      IMnpim_wSid_n_Sp_sub_bin[ibin][iq]->GetXaxis()->SetRangeUser(1,1.8);
      IMnpim_wSid_n_Sp_sub_bin[ibin][iq]->Draw("HIST");
      IMnpim_wSid_n_Sp_sub_bin_select[ibin][iq] = (TH1D*)IMnpim_wSid_n_Sp_sub_bin[ibin][iq]->Clone(Form("IMnpim_wSid_n_Sp_sub_bin%d_%d_select",ibin,iq));
      IMnpim_wSid_n_Sp_sub_bin_select[ibin][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmam_MIN,anacuts::Sigmam_MAX);
      IMnpim_wSid_n_Sp_sub_bin_select[ibin][iq]->SetLineColor(2);
      IMnpim_wSid_n_Sp_sub_bin_select[ibin][iq]->Draw("HISTsame");
      
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
  
  const int nwbin=3;
  //tuned bins 2d histograms
  TH2F* IMnpim_IMnpip_wSid_n_Sp_wbin[nwbin][nqcut];
  TH2F* IMnpim_IMnpip_wSid_n_Sm_wbin[nwbin][nqcut];
  TH2F* IMpippim_IMnpip_wSid_n_Sp_wbin[nwbin][nqcut];
  TH2F* IMpippim_IMnpip_wK0_n_wbin[nwbin][nqcut];
  TH2F* IMpippim_IMnpim_wSid_n_Sm_wbin[nwbin][nqcut];
  TH2F* IMpippim_IMnpim_wK0_n_wbin[nwbin][nqcut];
  for(int iq=0;iq<nqcut;iq++){
    for(int imerge=0;imerge<nwbin;imerge++){
      IMnpim_IMnpip_wSid_n_Sp_wbin[imerge][iq]= (TH2F*)fr[iq]->Get(Form("IMnpim_IMnpip_wSid_n_Sp_wbin%d",imerge));

      IMnpim_IMnpip_wSid_n_Sm_wbin[imerge][iq]= (TH2F*)fr[iq]->Get(Form("IMnpim_IMnpip_wSid_n_Sm_wbin%d",imerge));

      IMpippim_IMnpip_wSid_n_Sp_wbin[imerge][iq]= (TH2F*)fr[iq]->Get(Form("IMpippim_IMnpip_wSid_n_Sp_wbin%d",imerge));

      IMpippim_IMnpip_wK0_n_wbin[imerge][iq]= (TH2F*)fr[iq]->Get(Form("IMpippim_IMnpip_wK0_n_wbin%d",imerge));

      IMpippim_IMnpim_wSid_n_Sm_wbin[imerge][iq]= (TH2F*)fr[iq]->Get(Form("IMpippim_IMnpim_wSid_n_Sm_wbin%d",imerge));

      IMpippim_IMnpim_wK0_n_wbin[imerge][iq]= (TH2F*)fr[iq]->Get(Form("IMpippim_IMnpim_wK0_n_wbin%d",imerge));
    }
  }

  //correlation plots excluding overlap region
  TH2F* IMnpim_IMnpip_wK0_woSid_n_wbin[nwbin][nqcut];
  TH2F* IMnpim_IMnpip_woK0_wSid_n_woSm_wbin[nwbin][nqcut];//
  TH2F* IMnpim_IMnpip_woK0_wSid_n_woSp_wbin[nwbin][nqcut];//
  TH2F* MMnmiss_IMpippim_wK0_woSid_n_wbin[nwbin][nqcut];
  TH2F* MMnmiss_IMnpip_woK0_wSid_n_woSm_wbin[nwbin][nqcut];
  TH2F* MMnmiss_IMnpim_woK0_wSid_n_woSp_wbin[nwbin][nqcut];
  TH2F* q_nmom_wK0_woSid_n_wbin[nwbin][nqcut];
  TH2F* q_nmom_woK0_wSid_n_woSm_wbin[nwbin][nqcut];
  TH2F* q_nmom_woK0_wSid_n_woSp_wbin[nwbin][nqcut];
  TH2F* IMpippim_IMnpip_woK0_wSid_n_woSp_wbin[nwbin][nqcut];//for Sm template
  TH2F* IMpippim_IMnpip_woK0_wSid_n_woSm_wbin[nwbin][nqcut];//for Sp template
  TH2F* IMpippim_IMnpim_woK0_wSid_n_woSp_wbin[nwbin][nqcut];
  TH2F* IMpippim_IMnpim_woK0_wSid_n_woSm_wbin[nwbin][nqcut];
  TH2F* IMpippim_IMnpip_wK0_woSid_n_wbin[nwbin][nqcut];
  TH2F* IMpippim_IMnpim_wK0_woSid_n_wbin[nwbin][nqcut];
  TH2F* q_IMnpipi_wK0_woSid_n[nqcut];
  TH2F* q_IMnpipi_woK0_wSid_n_woSp[nqcut];
  TH2F* q_IMnpipi_woK0_wSid_n_woSm[nqcut];
  for(int iq=0;iq<nqcut;iq++){
    for(int imerge=0;imerge;nwbin;imerge++){
      IMnpim_IMnpip_wK0_woSid_n_wbin


    }
  }



  //projection hists
  TH1D* IMnpip_wK0_merge[10][nqcut];
  TH1D* IMnpip_wK0_merge_select[10][nqcut];
  TH1D* IMnpip_wK0_merge_lo[10][nqcut];
  TH1D* IMnpip_wK0_merge_hi[10][nqcut];
  TH1D* IMnpip_wK0_merge_lohi[10][nqcut];
  TH1D* IMpippim_Sp_merge[10][nqcut];
  TH1D* IMpippim_Sp_merge_select[10][nqcut];
  TH1D* IMpippim_Sp_merge_lo[10][nqcut];
  TH1D* IMpippim_Sp_merge_hi[10][nqcut];
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
  
  TF1 *constfitSpwK0[10][nqcut];
  TF1 *constfitK0wSp[10][nqcut];
  TF1 *constfitSmwK0[10][nqcut];
  TF1 *constfitK0wSm[10][nqcut];
  TF1 *constfitSpwSm[10][nqcut];
  TF1 *constfitSmwSp[10][nqcut];
  const double frSpwK0s[2][nqcut]={{1.1,1.1,1.1,1.1},
                                   {1.1,1.1,1.1,1.1}};
  const double frSpwK0e[2][nqcut]={{1.32,1.32,1.30,1.32},
                                   {1.32,1.30,1.32,1.32}};
  const double frK0wSps[2][nqcut]={{0.28,0.28,0.30,0.28},
                                   {0.30,0.30,0.30,0.30}};
  const double frK0wSpe[2][nqcut]={{0.57,0.56,0.57,0.57},
                                   {0.70,0.70,0.70,0.57}};
  const double frSmwK0s[2][nqcut]={{1.1,1.12,1.1,1.1},
                                   {1.17,1.17,1.17,1.17}};
  const double frSmwK0e[2][nqcut]={{1.32,1.3,1.27,1.27},
                                   {1.27,1.27,1.27,1.24}};
  const double frK0wSms[2][nqcut]={{0.28,0.30,0.28,0.28},
                                   {0.30,0.30,0.30,0.30}};
  const double frK0wSme[2][nqcut]={{0.57,0.57,0.57,0.57},
                                   {0.70,0.70,0.70,0.60}};
  const double frSpwSms[2][nqcut]={{1.14,1.14,1.12,1.14},
                                   {1.14,1.14,1.12,1.14}};
  const double frSpwSme[2][nqcut]={{1.28,1.28,1.30,1.28},
                                   {1.28,1.28,1.30,1.30}};
  const double frSmwSps[2][nqcut]={{1.14,1.14,1.12,1.14},
                                   {1.14,1.14,1.12,1.14}};
  const double frSmwSpe[2][nqcut]={{1.30,1.30,1.30,1.30},
                                   {1.30,1.30,1.30,1.30}};
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
      if(!RebinMode){
        IMnpip_wK0_merge_select[imerge][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmap_MIN,anacuts::Sigmap_MAX);
      }else{
        IMnpip_wK0_merge_select[imerge][iq]->GetXaxis()->SetRangeUser(
            IMnpip_wK0_merge_select[imerge][iq]->GetBinLowEdge(9),
            IMnpip_wK0_merge_select[imerge][iq]->GetBinLowEdge(10)
            );
      }
      IMnpip_wK0_merge_lohi[imerge][iq] = (TH1D*)IMnpip_wK0_merge[imerge][iq]->Clone(Form("IMnpip_wK0_merge_lohi_%d_%d",imerge,iq));
      IMnpip_wK0_merge_lohi[imerge][iq]->SetBinContent(9,0);
      IMnpip_wK0_merge_lohi[imerge][iq]->SetBinError(9,0);
      IMnpip_wK0_merge_lohi[imerge][iq]->SetLineColor(4);
      IMnpip_wK0_merge_lohi[imerge][iq]->Draw("HEsame");
      constfitSpwK0[imerge][iq] = new TF1(Form("fSpwK0%d%d",imerge,iq),"pol0");
      constfitSpwK0[imerge][iq]->SetLineColor(3);
      IMnpip_wK0_merge_select[imerge][iq]->Draw("HEsame");
      //IMnpip_wK0_merge_lo[imerge][iq]->Draw("HEsame");
      //IMnpip_wK0_merge_hi[imerge][iq]->Draw("HEsame");
      IMnpip_wK0_merge_lohi[imerge][iq]->Draw("HEsame");
      if(!FitNoWeight)IMnpip_wK0_merge_lohi[imerge][iq]->Fit(constfitSpwK0[imerge][iq],"","",frSpwK0s[imerge][iq],frSpwK0e[imerge][iq]);
      else IMnpip_wK0_merge_lohi[imerge][iq]->Fit(constfitSpwK0[imerge][iq],"w","",frSpwK0s[imerge][iq],frSpwK0e[imerge][iq]);
      double count_SpwK0 = IMnpip_wK0_merge_select[imerge][iq]->GetBinContent(9);
      double err_SpwK0 = IMnpip_wK0_merge_select[imerge][iq]->GetBinError(9);
      double bg_SpwK0 =  constfitSpwK0[imerge][iq]->GetParameter(0);
      double bg_SpwK0err =  constfitSpwK0[imerge][iq]->GetParError(0);
      cIMpippim_IMnpip_merge[imerge][iq]->cd(4);
      IMpippim_Sp_merge[imerge][iq]->Draw("HE");
      IMpippim_Sp_merge_select[imerge][iq] = (TH1D*)IMpippim_Sp_merge[imerge][iq]->Clone(Form("IMpippim_Sp_merge_select_%d_%d",imerge,iq));
      //IMpippim_Sp_merge_lo[imerge][iq] = (TH1D*)IMpippim_Sp_merge[imerge][iq]->Clone(Form("IMpippim_Sp_merge_lo_%d_%d",imerge,iq));
      //IMpippim_Sp_merge_hi[imerge][iq] = (TH1D*)IMpippim_Sp_merge[imerge][iq]->Clone(Form("IMpippim_Sp_merge_hi_%d_%d",imerge,iq));
      IMpippim_Sp_merge_select[imerge][iq]->SetLineColor(2);
      //IMpippim_Sp_merge_lo[imerge][iq]->SetLineColor(4);
      //IMpippim_Sp_merge_hi[imerge][iq]->SetLineColor(4);
      if(!RebinMode){
        IMpippim_Sp_merge_select[imerge][iq]->GetXaxis()->SetRangeUser(anacuts::pipi_MIN_narrow,anacuts::pipi_MAX_narrow);
      }else{
        IMpippim_Sp_merge_select[imerge][iq]->GetXaxis()->SetRangeUser(
            IMpippim_Sp_merge_select[imerge][iq]->GetBinLowEdge(10),
            IMpippim_Sp_merge_select[imerge][iq]->GetBinLowEdge(11));
      }
      IMpippim_Sp_merge_select[imerge][iq]->Draw("HEsame");
      IMpippim_Sp_merge_lohi[imerge][iq] = (TH1D*)IMpippim_Sp_merge[imerge][iq]->Clone(Form("IMpippim_Sp_merge_lohi_%d_%d",imerge,iq));
      IMpippim_Sp_merge_lohi[imerge][iq]->SetBinContent(10,0);
      IMpippim_Sp_merge_lohi[imerge][iq]->SetBinError(10,0);
      IMpippim_Sp_merge_lohi[imerge][iq]->SetLineColor(4);
      IMpippim_Sp_merge_lohi[imerge][iq]->Draw("HEsame");
      constfitK0wSp[imerge][iq] = new TF1(Form("fK0wSp%d%d",imerge,iq),"pol0");
      constfitK0wSp[imerge][iq]->SetLineColor(3);
      if(!FitNoWeight)IMpippim_Sp_merge_lohi[imerge][iq]->Fit(constfitK0wSp[imerge][iq],"","",frK0wSps[imerge][iq],frK0wSpe[imerge][iq]);
      else            IMpippim_Sp_merge_lohi[imerge][iq]->Fit(constfitK0wSp[imerge][iq],"w","",frK0wSps[imerge][iq],frK0wSpe[imerge][iq]);
      double count_K0wSp = IMpippim_Sp_merge_select[imerge][iq]->GetBinContent(10);
      double err_K0wSp = IMpippim_Sp_merge_select[imerge][iq]->GetBinError(10);
      double bg_K0wSp =  constfitK0wSp[imerge][iq]->GetParameter(0);
      double bg_K0wSperr = constfitK0wSp[imerge][iq]->GetParError(0);

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
      if(!RebinMode){
        IMnpim_wK0_merge_select[imerge][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmam_MIN,anacuts::Sigmam_MAX);
      }else{
        IMnpim_wK0_merge_select[imerge][iq]->GetXaxis()->SetRangeUser(
            IMnpim_wK0_merge_select[imerge][iq]->GetBinLowEdge(9),
            IMnpim_wK0_merge_select[imerge][iq]->GetBinLowEdge(10)); 
      }

      IMnpim_wK0_merge_select[imerge][iq]->Draw("HEsame");
      IMnpim_wK0_merge_lohi[imerge][iq] = (TH1D*)IMnpim_wK0_merge[imerge][iq]->Clone(Form("IMnpim_wK0_merge_lohi_%d_%d",imerge,iq));
      IMnpim_wK0_merge_lohi[imerge][iq]->SetBinContent(9,0);
      IMnpim_wK0_merge_lohi[imerge][iq]->SetBinError(9,0);
      IMnpim_wK0_merge_lohi[imerge][iq]->SetLineColor(4);
      IMnpim_wK0_merge_lohi[imerge][iq]->Draw("HEsame");
      constfitSmwK0[imerge][iq] = new TF1(Form("fSmwK0%d%d",imerge,iq),"pol0");
      constfitSmwK0[imerge][iq]->SetLineColor(3);
      if(!FitNoWeight)IMnpim_wK0_merge_lohi[imerge][iq]->Fit(constfitSmwK0[imerge][iq],"","",frSmwK0s[imerge][iq],frSmwK0e[imerge][iq]);
      else            IMnpim_wK0_merge_lohi[imerge][iq]->Fit(constfitSmwK0[imerge][iq],"w","",frSmwK0s[imerge][iq],frSmwK0e[imerge][iq]);
      double count_SmwK0 = IMnpim_wK0_merge_select[imerge][iq]->GetBinContent(9);
      double err_SmwK0 = IMnpim_wK0_merge_select[imerge][iq]->GetBinError(9);
      double bg_SmwK0 =  constfitSmwK0[imerge][iq]->GetParameter(0);
      double bg_SmwK0err =  constfitSmwK0[imerge][iq]->GetParError(0);
      
      
      cIMpippim_IMnpim_merge[imerge][iq]->cd(4);
      IMpippim_Sm_merge[imerge][iq]->Draw("HE");
      IMpippim_Sm_merge_select[imerge][iq] = (TH1D*)IMpippim_Sm_merge[imerge][iq]->Clone(Form("IMpippim_Sm_merge_select_%d_%d",imerge,iq));
      IMpippim_Sm_merge_select[imerge][iq]->SetLineColor(2);
      if(!RebinMode){
        IMpippim_Sm_merge_select[imerge][iq]->GetXaxis()->SetRangeUser(anacuts::pipi_MIN_narrow,anacuts::pipi_MAX_narrow);
      }else{
        IMpippim_Sm_merge_select[imerge][iq]->GetXaxis()->SetRangeUser(
            IMpippim_Sm_merge_select[imerge][iq]->GetBinLowEdge(10),
            IMpippim_Sm_merge_select[imerge][iq]->GetBinLowEdge(11));
      }
        
      IMpippim_Sm_merge_select[imerge][iq]->Draw("HEsame");
      IMpippim_Sm_merge_lohi[imerge][iq] = (TH1D*)IMpippim_Sm_merge[imerge][iq]->Clone(Form("IMpippim_Sm_merge_lohi_%d_%d",imerge,iq));
      IMpippim_Sm_merge_lohi[imerge][iq]->SetLineColor(4);
      IMpippim_Sm_merge_lohi[imerge][iq]->SetBinContent(10,0);
      IMpippim_Sm_merge_lohi[imerge][iq]->SetBinError(10,0);
      IMpippim_Sm_merge_lohi[imerge][iq]->Draw("HEsame");
      constfitK0wSm[imerge][iq] = new TF1(Form("fK0wSm%d%d",imerge,iq),"pol0");
      constfitK0wSm[imerge][iq]->SetLineColor(3);
      if(!FitNoWeight) IMpippim_Sm_merge_lohi[imerge][iq]->Fit(constfitK0wSm[imerge][iq],"","",frK0wSms[imerge][iq],frK0wSme[imerge][iq]);
      else             IMpippim_Sm_merge_lohi[imerge][iq]->Fit(constfitK0wSm[imerge][iq],"w","",frK0wSms[imerge][iq],frK0wSme[imerge][iq]);
      double count_K0wSm = IMpippim_Sm_merge_select[imerge][iq]->GetBinContent(10);
      double err_K0wSm = IMpippim_Sm_merge_select[imerge][iq]->GetBinError(10);
      double bg_K0wSm =  constfitK0wSm[imerge][iq]->GetParameter(0);
      double bg_K0wSmerr =  constfitK0wSm[imerge][iq]->GetParError(0);
      
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
      if(!RebinMode){
        IMnpip_Sm_merge_select[imerge][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmap_MIN,anacuts::Sigmap_MAX);
      }else{
        IMnpip_Sm_merge_select[imerge][iq]->GetXaxis()->SetRangeUser(
            IMnpip_Sm_merge_select[imerge][iq]->GetBinLowEdge(9),
            IMnpip_Sm_merge_select[imerge][iq]->GetBinLowEdge(10));
      }
      
      IMnpip_Sm_merge_select[imerge][iq]->Draw("HEsame");
      IMnpip_Sm_merge_lohi[imerge][iq] = (TH1D*)IMnpip_Sm_merge[imerge][iq]->Clone(Form("IMnpip_Sm_merge_lohi_%d_%d",imerge,iq));
      IMnpip_Sm_merge_lohi[imerge][iq]->SetLineColor(4);
      IMnpip_Sm_merge_lohi[imerge][iq]->SetBinContent(9,0);
      IMnpip_Sm_merge_lohi[imerge][iq]->SetBinError(9,0);
      IMnpip_Sm_merge_lohi[imerge][iq]->Draw("HEsame");
      constfitSpwSm[imerge][iq] = new TF1(Form("fSpwSm%d%d",imerge,iq),"pol0");
      constfitSpwSm[imerge][iq]->SetLineColor(3);
      if(!FitNoWeight)IMnpip_Sm_merge_lohi[imerge][iq]->Fit(constfitSpwSm[imerge][iq],"","",frSpwSms[imerge][iq],frSpwSme[imerge][iq]);
      else            IMnpip_Sm_merge_lohi[imerge][iq]->Fit(constfitSpwSm[imerge][iq],"w","",frSpwSms[imerge][iq],frSpwSme[imerge][iq]);
      double count_SpwSm = IMnpip_Sm_merge_select[imerge][iq]->GetBinContent(9);
      double err_SpwSm = IMnpip_Sm_merge_select[imerge][iq]->GetBinError(9);
      double bg_SpwSm =  constfitSpwSm[imerge][iq]->GetParameter(0);
      double bg_SpwSmerr =  constfitSpwSm[imerge][iq]->GetParError(0);

      cIMnpim_IMnpip_merge[imerge][iq]->cd(4);
      IMnpim_Sp_merge[imerge][iq]->GetXaxis()->SetRangeUser(1,1.4);
      IMnpim_Sp_merge[imerge][iq]->SetMinimum(0);
      IMnpim_Sp_merge[imerge][iq]->Draw("HE");
      IMnpim_Sp_merge_select[imerge][iq] = (TH1D*)IMnpim_Sp_merge[imerge][iq]->Clone(Form("IMnpim_Sp_merge_select_%d_%d",imerge,iq));
      IMnpim_Sp_merge_select[imerge][iq]->SetLineColor(2);
      if(!RebinMode){
        IMnpim_Sp_merge_select[imerge][iq]->GetXaxis()->SetRangeUser(anacuts::Sigmam_MIN,anacuts::Sigmam_MAX);
      }else{
        IMnpim_Sp_merge_select[imerge][iq]->GetXaxis()->SetRangeUser(
            IMnpim_Sp_merge_select[imerge][iq]->GetBinLowEdge(9),
            IMnpim_Sp_merge_select[imerge][iq]->GetBinLowEdge(10));
      }
      
      IMnpim_Sp_merge_select[imerge][iq]->Draw("HEsame");
      IMnpim_Sp_merge_lohi[imerge][iq] = (TH1D*)IMnpim_Sp_merge[imerge][iq]->Clone(Form("IMnpim_Sp_merge_lohi_%d_%d",imerge,iq));
      IMnpim_Sp_merge_lohi[imerge][iq]->SetLineColor(4);
      IMnpim_Sp_merge_lohi[imerge][iq]->SetBinContent(9,0);
      IMnpim_Sp_merge_lohi[imerge][iq]->SetBinError(9,0);
      IMnpim_Sp_merge_lohi[imerge][iq]->Draw("HEsame");
      constfitSmwSp[imerge][iq] = new TF1(Form("fSmwSp%d%d",imerge,iq),"pol0");
      constfitSmwSp[imerge][iq]->SetLineColor(3);
      if(!FitNoWeight)IMnpim_Sp_merge_lohi[imerge][iq]->Fit(constfitSmwSp[imerge][iq],"","",frSmwSps[imerge][iq],frSmwSpe[imerge][iq]);
      else            IMnpim_Sp_merge_lohi[imerge][iq]->Fit(constfitSmwSp[imerge][iq],"w","",frSmwSps[imerge][iq],frSmwSpe[imerge][iq]);
      double count_SmwSp = IMnpim_Sp_merge_select[imerge][iq]->GetBinContent(9);
      double err_SmwSp = IMnpim_Sp_merge_select[imerge][iq]->GetBinError(9);
      double bg_SmwSp =  constfitSmwSp[imerge][iq]->GetParameter(0);
      double bg_SmwSperr =  constfitSmwSp[imerge][iq]->GetParError(0);
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
  TString pdfname = "K0SigmaTemp.pdf";
  if(RebinMode) pdfname = "K0SigmaTemp_Rebin.pdf";
  for(int i=0;i<size;i++){
    //pdf->NewPage();
    c= (TCanvas*)next();
    c->Draw();
    c->cd();
    //inside the canvas
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
  if(!RebinMode) fout = new TFile("Decomposition.root","RECREATE");
  else           fout = new TFile("Decomposition_Rebin.root","RECREATE");
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
