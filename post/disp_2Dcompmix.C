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

void disp_2Dcompmix()
{
  TFile *fr = TFile::Open("evanaIMpisigma_npippim_v202_out_iso.root");
  TFile *fmix = TFile::Open("evanaIMpisigma_npippim_v202_MIX_cut4_out_iso.root");
  fr->Print() ;
  fmix->Print();
   
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0.);  


  //just plot each 2d row histograms before scaling
  TH2D* MMnmiss_IMnpip_woK0_woSm_data = (TH2D*)fr->Get("MMnmiss_IMnpip_dE_woK0_woSm");
  TCanvas *cMMnmiss_IMnpip_woK0_woSm_data = new TCanvas("cMMnmiss_IMnpip_woK0_woSm_data","cMMnmiss_IMnpip_woK0_woSm_data",800,800);
  MMnmiss_IMnpip_woK0_woSm_data->Draw("colz");

  TH2D* MMnmiss_IMnpim_woK0_woSp_data = (TH2D*)fr->Get("MMnmiss_IMnpim_dE_woK0_woSp");
  TCanvas *cMMnmiss_IMnpim_woK0_woSp_data = new TCanvas("cMMnmiss_IMnpim_woK0_woSp_data","cMMnmiss_IMnpim_woK0_woSp_data",800,800);
  MMnmiss_IMnpim_woK0_woSp_data->Draw("colz");

  TH2D* MMnmiss_IMnpip_woK0_woSm_mix = (TH2D*)fmix->Get("MMnmiss_IMnpip_dE_woK0_woSm");
  TCanvas *cMMnmiss_IMnpip_woK0_woSm_mix = new TCanvas("cMMnmiss_IMnpip_dE_woK0_woSm_mix","cMMnmiss_IMnpip_woK0_woSm_mix",800,800);
  MMnmiss_IMnpip_woK0_woSm_mix->Draw("colz");
  
  TH2D* MMnmiss_IMnpim_woK0_woSp_mix = (TH2D*)fmix->Get("MMnmiss_IMnpim_dE_woK0_woSp");
  TCanvas *cMMnmiss_IMnpim_woK0_woSp_mix = new TCanvas("cMMnmiss_IMnpim_woK0_woSp_mix","cMMnmiss_IMnpim_woK0_woSp_mix",800,800);
  MMnmiss_IMnpim_woK0_woSp_mix->Draw("colz");
  
  TH2D* IMnpim_IMnpip_woK0_woSp_vici_data = (TH2D*)fr->Get("IMnpim_IMnpip_dE_woK0_woSp_vici");
  TCanvas *cIMnpim_IMnpip_woK0_woSp_vici_data = new TCanvas("cIMnpim_IMnpip_woK0_woSp_vici_data","cIMnpim_IMnpip_woK0_woSp_vici_data",800,800);
  IMnpim_IMnpip_woK0_woSp_vici_data->Draw("colz");

  TH2D* IMnpim_IMnpip_woK0_woSm_vici_data = (TH2D*)fr->Get("IMnpim_IMnpip_dE_woK0_woSm_vici");
  TCanvas *cIMnpim_IMnpip_woK0_woSm_vici_data = new TCanvas("cIMnpim_IMnpip_woK0_woSm_vici_data","cIMnpim_IMnpip_woK0_woSm_vici_data",800,800);
  IMnpim_IMnpip_woK0_woSm_vici_data->Draw("colz");
  //

  TH2D* MMnmiss_IMnpip_woK0_woSm_vici_data = (TH2D*)fr->Get("MMnmiss_IMnpip_dE_woK0_woSm_vici");
  TH2D* MMnmiss_IMnpim_woK0_woSp_vici_data = (TH2D*)fr->Get("MMnmiss_IMnpim_dE_woK0_woSp_vici");
  TH2D* MMnmiss_IMnpip_woK0_woSm_viciext_data = (TH2D*)fr->Get("MMnmiss_IMnpip_dE_woK0_woSm_viciext");
  TH2D* MMnmiss_IMnpim_woK0_woSp_viciext_data = (TH2D*)fr->Get("MMnmiss_IMnpim_dE_woK0_woSp_viciext");
  TH2D* MMnmiss_IMnpip_woK0_woSm_vici_mix = (TH2D*)fmix->Get("MMnmiss_IMnpip_dE_woK0_woSm_vici");
  TH2D* MMnmiss_IMnpim_woK0_woSp_vici_mix = (TH2D*)fmix->Get("MMnmiss_IMnpim_dE_woK0_woSp_vici");
  TH2D* MMnmiss_IMnpip_woK0_woSm_viciext_mix = (TH2D*)fmix->Get("MMnmiss_IMnpip_dE_woK0_woSm_viciext");
  TH2D* MMnmiss_IMnpim_woK0_woSp_viciext_mix = (TH2D*)fmix->Get("MMnmiss_IMnpim_dE_woK0_woSp_viciext");
  
  //
  //compare data and mixed events
  //
  TCanvas *cMMnmiss_IMnpip_woK0_woSm_comp = new TCanvas("cMMnmiss_IMnpip_woK0_woSm_comp","cMMnmiss_IMnpip_woK0_woSm_data_comp",800,800);
  cMMnmiss_IMnpip_woK0_woSm_comp->Divide(2,2,0,0);
  cMMnmiss_IMnpip_woK0_woSm_comp->cd(3);
  TH2D* MMnmiss_IMnpip_woK0_woSm_data_zoom = (TH2D*)MMnmiss_IMnpip_woK0_woSm_data->Clone("MMnmiss_IMnpip_woK0_woSm_data_zoom");
  TH2D* MMnmiss_IMnpip_woK0_woSm_mix_zoom = (TH2D*)MMnmiss_IMnpip_woK0_woSm_mix->Clone("MMnmiss_IMnpip_woK0_woSm_mix_zoom");
  MMnmiss_IMnpip_woK0_woSm_data_zoom->GetXaxis()->SetRangeUser(1.1,1.3);
  MMnmiss_IMnpip_woK0_woSm_data_zoom->GetYaxis()->SetRangeUser(0.7,1.1);
  MMnmiss_IMnpip_woK0_woSm_vici_data->GetXaxis()->SetRangeUser(1.1,1.3);
  MMnmiss_IMnpip_woK0_woSm_vici_data->GetYaxis()->SetRangeUser(0.7,1.1);
  MMnmiss_IMnpip_woK0_woSm_vici_mix->GetXaxis()->SetRangeUser(1.1,1.3);
  MMnmiss_IMnpip_woK0_woSm_vici_mix->GetYaxis()->SetRangeUser(0.7,1.1);
  MMnmiss_IMnpip_woK0_woSm_data_zoom->SetTitle("");
  //MMnmiss_IMnpip_woK0_woSm_data_zoom->GetXaxis()->SetLabelOffset(0.005);
  //MMnmiss_IMnpip_woK0_woSm_data_zoom->GetYaxis()->SetLabelOffset(0.005);
  //MMnmiss_IMnpip_woK0_woSm_data_zoom->Draw("colz");
  MMnmiss_IMnpip_woK0_woSm_vici_data->Draw("colz");

  cMMnmiss_IMnpip_woK0_woSm_comp->cd(1);
  TH1D* IMnpip_woK0_woSm_vici_data = (TH1D*)MMnmiss_IMnpip_woK0_woSm_vici_data->ProjectionX("IMnpip_woK0_woSm_vici_data");
  TH1D* IMnpip_woK0_woSm_vici_mix  = (TH1D*)MMnmiss_IMnpip_woK0_woSm_vici_mix->ProjectionX("IMnpip_woK0_woSm_vici_mix");
  IMnpip_woK0_woSm_vici_data->Draw("HE");
  IMnpip_woK0_woSm_vici_mix->SetLineColor(2);
  IMnpip_woK0_woSm_vici_mix->Draw("HEsame");
  
  cMMnmiss_IMnpip_woK0_woSm_comp->cd(4);
  TH1D* MMnmiss_woK0_woSm_vici_data = (TH1D*)MMnmiss_IMnpip_woK0_woSm_vici_data->ProjectionY("MMnmiss_woK0_woSm_vici_data");
  TH1D* MMnmiss_woK0_woSm_vici_mix  = (TH1D*)MMnmiss_IMnpip_woK0_woSm_vici_mix->ProjectionY("MMnmiss_woK0_woSm_vici_mix");
  TGraphErrors *gr_MMnmiss_woK0_woSm_vici_data = new TGraphErrors();
  TGraphErrors *gr_MMnmiss_woK0_woSm_vici_mix = new TGraphErrors();
  HistToRorateGraph(MMnmiss_woK0_woSm_vici_data,*gr_MMnmiss_woK0_woSm_vici_data);
  HistToRorateGraph(MMnmiss_woK0_woSm_vici_mix,*gr_MMnmiss_woK0_woSm_vici_mix);
  gr_MMnmiss_woK0_woSm_vici_data->Draw("AP");
  gr_MMnmiss_woK0_woSm_vici_mix->SetMarkerColor(2);
  gr_MMnmiss_woK0_woSm_vici_mix->Draw("P");
  
  //
  //compare data and mixed events (ext hist.)
  //
  TCanvas *cMMnmiss_IMnpip_woK0_woSm_comp2 = new TCanvas("cMMnmiss_IMnpip_woK0_woSm_comp2","cMMnmiss_IMnpip_woK0_woSm_data_comp2",800,800);
  cMMnmiss_IMnpip_woK0_woSm_comp2->Divide(2,2,0,0);
  cMMnmiss_IMnpip_woK0_woSm_comp2->cd(3);
  TH2D* MMnmiss_IMnpip_woK0_woSm_data_zoom2 = (TH2D*)MMnmiss_IMnpip_woK0_woSm_data->Clone("MMnmiss_IMnpip_woK0_woSm_data_zoom2");
  TH2D* MMnmiss_IMnpip_woK0_woSm_mix_zoom2 = (TH2D*)MMnmiss_IMnpip_woK0_woSm_mix->Clone("MMnmiss_IMnpip_woK0_woSm_mix_zoom2");
  MMnmiss_IMnpip_woK0_woSm_data_zoom2->GetXaxis()->SetRangeUser(1.1,1.3);
  MMnmiss_IMnpip_woK0_woSm_data_zoom2->GetYaxis()->SetRangeUser(0.0,1.1);
  MMnmiss_IMnpip_woK0_woSm_viciext_data->GetXaxis()->SetRangeUser(1.1,1.3);
  MMnmiss_IMnpip_woK0_woSm_viciext_data->GetYaxis()->SetRangeUser(0.0,1.1);
  MMnmiss_IMnpip_woK0_woSm_viciext_mix->GetXaxis()->SetRangeUser(1.1,1.3);
  MMnmiss_IMnpip_woK0_woSm_viciext_mix->GetYaxis()->SetRangeUser(0.0,1.1);
  MMnmiss_IMnpip_woK0_woSm_data_zoom2->SetTitle("");
  //MMnmiss_IMnpip_woK0_woSm_data_zoom2->GetXaxis()->SetLabelOffset(0.005);
  //MMnmiss_IMnpip_woK0_woSm_data_zoom2->GetYaxis()->SetLabelOffset(0.005);
  //MMnmiss_IMnpip_woK0_woSm_data_zoom2->Draw("colz");
  MMnmiss_IMnpip_woK0_woSm_viciext_data->Draw("colz");

  cMMnmiss_IMnpip_woK0_woSm_comp2->cd(1);
  TH1D* IMnpip_woK0_woSm_viciext_data = (TH1D*)MMnmiss_IMnpip_woK0_woSm_viciext_data->ProjectionX("IMnpip_woK0_woSm_viciext_data");
  TH1D* IMnpip_woK0_woSm_viciext_mix  = (TH1D*)MMnmiss_IMnpip_woK0_woSm_viciext_mix->ProjectionX("IMnpip_woK0_woSm_viciext_mix");
  IMnpip_woK0_woSm_viciext_data->Draw("HE");
  IMnpip_woK0_woSm_viciext_mix->SetLineColor(2);
  IMnpip_woK0_woSm_viciext_mix->Draw("HEsame");
  
  cMMnmiss_IMnpip_woK0_woSm_comp2->cd(4);
  TH1D* MMnmiss_woK0_woSm_viciext_data = (TH1D*)MMnmiss_IMnpip_woK0_woSm_viciext_data->ProjectionY("MMnmiss_woK0_woSm_viciext_data");
  TH1D* MMnmiss_woK0_woSm_viciext_mix  = (TH1D*)MMnmiss_IMnpip_woK0_woSm_viciext_mix->ProjectionY("MMnmiss_woK0_woSm_viciext_mix");

  TGraphErrors *gr_MMnmiss_woK0_woSm_viciext_data = new TGraphErrors();
  TGraphErrors *gr_MMnmiss_woK0_woSm_viciext_mix = new TGraphErrors();
  HistToRorateGraph(MMnmiss_woK0_woSm_viciext_data,*gr_MMnmiss_woK0_woSm_viciext_data);
  HistToRorateGraph(MMnmiss_woK0_woSm_viciext_mix,*gr_MMnmiss_woK0_woSm_viciext_mix);
  gr_MMnmiss_woK0_woSm_viciext_data->Draw("AP");
  gr_MMnmiss_woK0_woSm_viciext_mix->SetMarkerColor(2);
  gr_MMnmiss_woK0_woSm_viciext_mix->Draw("P");

  TCanvas *cMMnmiss_IMnpip_woK0_woSm_data_sub2 = new TCanvas("cMMnmiss_IMnpip_woK0_woSm_data_sub2","cMMnmiss_IMnpip_woK0_woSm_data_sub2",800,800);
  TH2D* MMnmiss_IMnpip_woK0_woSm_data_sub = (TH2D*)MMnmiss_IMnpip_woK0_woSm_data->Clone("MMnmiss_IMnpip_woK0_woSm_data_sub");
  //gStyle->SetPalette(53);
  MMnmiss_IMnpip_woK0_woSm_data_sub->Add(MMnmiss_IMnpip_woK0_woSm_mix,-1.0);
  TH2D* MMnmiss_IMnpip_woK0_woSm_data_subRebin = (TH2D*)MMnmiss_IMnpip_woK0_woSm_data_sub->Clone("MMnmiss_IMnpip_woK0_woSm_data_subRebin");
  MMnmiss_IMnpip_woK0_woSm_data_subRebin->RebinX(2);
  MMnmiss_IMnpip_woK0_woSm_data_subRebin->RebinY(2);
  MMnmiss_IMnpip_woK0_woSm_data_subRebin->Draw("colz");

  //and subtract
  TCanvas *cMMnmiss_IMnpip_woK0_woSm_data_sub = new TCanvas("cMMnmiss_IMnpip_woK0_woSm_data_sub","cMMnmiss_IMnpip_woK0_woSm_data_sub",800,800);
  cMMnmiss_IMnpip_woK0_woSm_data_sub->Divide(2,2,0,0);
  cMMnmiss_IMnpip_woK0_woSm_data_sub->cd(3);
  //MMnmiss_IMnpip_woK0_woSm_data_sub->RebinX(2);
  //MMnmiss_IMnpip_woK0_woSm_data_sub->RebinY(2);
  MMnmiss_IMnpip_woK0_woSm_data_sub->SetTitle("");
  MMnmiss_IMnpip_woK0_woSm_data_sub->Draw("colz");
  cMMnmiss_IMnpip_woK0_woSm_data_sub->cd(1);
  TH1D* IMnpip_woK0_woSm_data_sub = (TH1D*)MMnmiss_IMnpip_woK0_woSm_data_sub->ProjectionX("IMnpip_woK0_woSm_data_sub");
  IMnpip_woK0_woSm_data_sub->Draw("HE");
  cMMnmiss_IMnpip_woK0_woSm_data_sub->cd(4);
  TH1D* MMnmiss_woK0_woSm_data_sub = (TH1D*)MMnmiss_IMnpip_woK0_woSm_data_sub->ProjectionY("MMnmiss_woK0_woSm_data_sub");
  TGraphErrors *gr_MMnmiss_woK0_woSm_data_sub = new TGraphErrors();
  HistToRorateGraph(MMnmiss_woK0_woSm_data_sub,*gr_MMnmiss_woK0_woSm_data_sub);
  gr_MMnmiss_woK0_woSm_data_sub->Draw("AP");
  //gr_MMnmiss_woK0_woSm_data_sub->Print("all");
   
 


  //draw vicinity subtracted
  TCanvas *cMMnmiss_IMnpip_woK0_woSm_vici_data_sub = new TCanvas("cMMnmiss_IMnpip_woK0_woSm_vici_data_sub","cMMnmiss_IMnpip_woK0_woSm_vici_data_sub",800,800);
  cMMnmiss_IMnpip_woK0_woSm_vici_data_sub->Divide(2,2,0,0);
  cMMnmiss_IMnpip_woK0_woSm_vici_data_sub->cd(3);
  TH2D* MMnmiss_IMnpip_woK0_woSm_vici_data_sub = (TH2D*)MMnmiss_IMnpip_woK0_woSm_vici_data->Clone("MMnmiss_IMnpip_woK0_woSm_vici_data_sub");
  MMnmiss_IMnpip_woK0_woSm_vici_data_sub->Add(MMnmiss_IMnpip_woK0_woSm_vici_mix,-1.0);
  MMnmiss_IMnpip_woK0_woSm_vici_data_sub->SetTitle("");
  //MMnmiss_IMnpip_woK0_woSm_vici_data_sub->GetXaxis()->SetLabelOffset(0.005);
  //MMnmiss_IMnpip_woK0_woSm_vici_data_sub->GetYaxis()->SetLabelOffset(0.005);
  MMnmiss_IMnpip_woK0_woSm_vici_data_sub->Draw("colz");
  
  cMMnmiss_IMnpip_woK0_woSm_vici_data_sub->cd(1);
  TH1D* IMnpip_woK0_woSm_vici_data_sub = (TH1D*)MMnmiss_IMnpip_woK0_woSm_vici_data_sub->ProjectionX("IMnpip_woK0_woSm_vici_data_sub");
  IMnpip_woK0_woSm_vici_data_sub->Draw("HE");
  
  
  cMMnmiss_IMnpip_woK0_woSm_vici_data_sub->cd(4);
  TH1D* MMnmiss_woK0_woSm_vici_data_sub = (TH1D*)MMnmiss_IMnpip_woK0_woSm_vici_data_sub->ProjectionY("MMnmiss_woK0_woSm_vici_data_sub");
  TGraphErrors *gr_MMnmiss_woK0_woSm_vici_data_sub = new TGraphErrors();
  HistToRorateGraph(MMnmiss_woK0_woSm_vici_data_sub,*gr_MMnmiss_woK0_woSm_vici_data_sub);
  gr_MMnmiss_woK0_woSm_vici_data_sub->Draw("AP");
  
  //
  //Sigma- mode
  //
  //real data

  //
  //compare data and mixed events
  //
  TCanvas *cMMnmiss_IMnpim_woK0_woSp_comp = new TCanvas("cMMnmiss_IMnpim_woK0_woSp_comp","cMMnmiss_IMnpim_woK0_woSp_data_comp",800,800);
  cMMnmiss_IMnpim_woK0_woSp_comp->Divide(2,2,0,0);
  cMMnmiss_IMnpim_woK0_woSp_comp->cd(3);
  TH2D* MMnmiss_IMnpim_woK0_woSp_data_zoom = (TH2D*)MMnmiss_IMnpim_woK0_woSp_data->Clone("MMnmiss_IMnpim_woK0_woSp_data_zoom");
  TH2D* MMnmiss_IMnpim_woK0_woSp_mix_zoom = (TH2D*)MMnmiss_IMnpim_woK0_woSp_mix->Clone("MMnmiss_IMnpim_woK0_woSp_mix_zoom");
  MMnmiss_IMnpim_woK0_woSp_data_zoom->GetXaxis()->SetRangeUser(1.1,1.3);
  MMnmiss_IMnpim_woK0_woSp_data_zoom->GetYaxis()->SetRangeUser(0.7,1.1);
  MMnmiss_IMnpim_woK0_woSp_vici_data->GetXaxis()->SetRangeUser(1.1,1.3);
  MMnmiss_IMnpim_woK0_woSp_vici_data->GetYaxis()->SetRangeUser(0.7,1.1);
  MMnmiss_IMnpim_woK0_woSp_vici_mix->GetXaxis()->SetRangeUser(1.1,1.3);
  MMnmiss_IMnpim_woK0_woSp_vici_mix->GetYaxis()->SetRangeUser(0.7,1.1);
  MMnmiss_IMnpim_woK0_woSp_data_zoom->SetTitle("");
  //MMnmiss_IMnpim_woK0_woSp_data_zoom->GetXaxis()->SetLabelOffset(0.005);
  //MMnmiss_IMnpim_woK0_woSp_data_zoom->GetYaxis()->SetLabelOffset(0.005);
  //MMnmiss_IMnpim_woK0_woSp_data_zoom->Draw("colz");
  MMnmiss_IMnpim_woK0_woSp_vici_data->Draw("colz");

  cMMnmiss_IMnpim_woK0_woSp_comp->cd(1);
  TH1D* IMnpim_woK0_woSp_vici_data = (TH1D*)MMnmiss_IMnpim_woK0_woSp_vici_data->ProjectionX("IMnpim_woK0_woSp_vici_data");
  TH1D* IMnpim_woK0_woSp_vici_mix  = (TH1D*)MMnmiss_IMnpim_woK0_woSp_vici_mix->ProjectionX("IMnpim_woK0_woSp_vici_mix");
  IMnpim_woK0_woSp_vici_data->Draw("HE");
  IMnpim_woK0_woSp_vici_mix->SetLineColor(2);
  IMnpim_woK0_woSp_vici_mix->Draw("HEsame");
  
  cMMnmiss_IMnpim_woK0_woSp_comp->cd(4);
  TH1D* MMnmiss_woK0_woSp_vici_data = (TH1D*)MMnmiss_IMnpim_woK0_woSp_vici_data->ProjectionY("MMnmiss_woK0_woSp_vici_data");
  TH1D* MMnmiss_woK0_woSp_vici_mix  = (TH1D*)MMnmiss_IMnpim_woK0_woSp_vici_mix->ProjectionY("MMnmiss_woK0_woSp_vici_mix");
  TGraphErrors *gr_MMnmiss_woK0_woSp_vici_data = new TGraphErrors();
  TGraphErrors *gr_MMnmiss_woK0_woSp_vici_mix = new TGraphErrors();
  HistToRorateGraph(MMnmiss_woK0_woSp_vici_data,*gr_MMnmiss_woK0_woSp_vici_data);
  HistToRorateGraph(MMnmiss_woK0_woSp_vici_mix,*gr_MMnmiss_woK0_woSp_vici_mix);
  gr_MMnmiss_woK0_woSp_vici_data->Draw("AP");
  gr_MMnmiss_woK0_woSp_vici_mix->SetMarkerColor(2);
  gr_MMnmiss_woK0_woSp_vici_mix->Draw("P");
  


  //
  //compare data and mixed events (ext hist.)
  //
  TCanvas *cMMnmiss_IMnpim_woK0_woSp_comp2 = new TCanvas("cMMnmiss_IMnpim_woK0_woSp_comp2","cMMnmiss_IMnpim_woK0_woSp_data_comp2",800,800);
  cMMnmiss_IMnpim_woK0_woSp_comp2->Divide(2,2,0,0);
  cMMnmiss_IMnpim_woK0_woSp_comp2->cd(3);
  TH2D* MMnmiss_IMnpim_woK0_woSp_data_zoom2 = (TH2D*)MMnmiss_IMnpim_woK0_woSp_data->Clone("MMnmiss_IMnpim_woK0_woSp_data_zoom2");
  TH2D* MMnmiss_IMnpim_woK0_woSp_mix_zoom2 = (TH2D*)MMnmiss_IMnpim_woK0_woSp_mix->Clone("MMnmiss_IMnpim_woK0_woSp_mix_zoom2");
  MMnmiss_IMnpim_woK0_woSp_data_zoom2->GetXaxis()->SetRangeUser(1.1,1.3);
  MMnmiss_IMnpim_woK0_woSp_data_zoom2->GetYaxis()->SetRangeUser(0.0,1.1);
  MMnmiss_IMnpim_woK0_woSp_viciext_data->GetXaxis()->SetRangeUser(1.1,1.3);
  MMnmiss_IMnpim_woK0_woSp_viciext_data->GetYaxis()->SetRangeUser(0.0,1.1);
  MMnmiss_IMnpim_woK0_woSp_viciext_mix->GetXaxis()->SetRangeUser(1.1,1.3);
  MMnmiss_IMnpim_woK0_woSp_viciext_mix->GetYaxis()->SetRangeUser(0.0,1.1);
  MMnmiss_IMnpim_woK0_woSp_data_zoom2->SetTitle("");
  //MMnmiss_IMnpim_woK0_woSp_data_zoom2->GetXaxis()->SetLabelOffset(0.005);
  //MMnmiss_IMnpim_woK0_woSp_data_zoom2->GetYaxis()->SetLabelOffset(0.005);
  //MMnmiss_IMnpim_woK0_woSp_data_zoom2->Draw("colz");
  MMnmiss_IMnpim_woK0_woSp_viciext_data->Draw("colz");

  cMMnmiss_IMnpim_woK0_woSp_comp2->cd(1);
  TH1D* IMnpim_woK0_woSp_viciext_data = (TH1D*)MMnmiss_IMnpim_woK0_woSp_viciext_data->ProjectionX("IMnpim_woK0_woSp_viciext_data");
  TH1D* IMnpim_woK0_woSp_viciext_mix  = (TH1D*)MMnmiss_IMnpim_woK0_woSp_viciext_mix->ProjectionX("IMnpim_woK0_woSp_viciext_mix");
  IMnpim_woK0_woSp_viciext_data->Draw("HE");
  IMnpim_woK0_woSp_viciext_mix->SetLineColor(2);
  IMnpim_woK0_woSp_viciext_mix->Draw("HEsame");
  
  cMMnmiss_IMnpim_woK0_woSp_comp2->cd(4);
  TH1D* MMnmiss_woK0_woSp_viciext_data = (TH1D*)MMnmiss_IMnpim_woK0_woSp_viciext_data->ProjectionY("MMnmiss_woK0_woSp_viciext_data");
  TH1D* MMnmiss_woK0_woSp_viciext_mix  = (TH1D*)MMnmiss_IMnpim_woK0_woSp_viciext_mix->ProjectionY("MMnmiss_woK0_woSp_viciext_mix");

  TGraphErrors *gr_MMnmiss_woK0_woSp_viciext_data = new TGraphErrors();
  TGraphErrors *gr_MMnmiss_woK0_woSp_viciext_mix = new TGraphErrors();
  HistToRorateGraph(MMnmiss_woK0_woSp_viciext_data,*gr_MMnmiss_woK0_woSp_viciext_data);
  HistToRorateGraph(MMnmiss_woK0_woSp_viciext_mix,*gr_MMnmiss_woK0_woSp_viciext_mix);
  gr_MMnmiss_woK0_woSp_viciext_data->Draw("AP");
  gr_MMnmiss_woK0_woSp_viciext_mix->SetMarkerColor(2);
  gr_MMnmiss_woK0_woSp_viciext_mix->Draw("P");

  TCanvas *cMMnmiss_IMnpim_woK0_woSp_data_sub2 = new TCanvas("cMMnmiss_IMnpim_woK0_woSp_data_sub2","cMMnmiss_IMnpim_woK0_woSp_data_sub2",800,800);
  TH2D* MMnmiss_IMnpim_woK0_woSp_data_sub = (TH2D*)MMnmiss_IMnpim_woK0_woSp_data->Clone("MMnmiss_IMnpim_woK0_woSp_data_sub");
  //gStyle->SetPalette(53);
  MMnmiss_IMnpim_woK0_woSp_data_sub->Add(MMnmiss_IMnpim_woK0_woSp_mix,-1.0);
  TH2D* MMnmiss_IMnpim_woK0_woSp_data_subRebin = (TH2D*)MMnmiss_IMnpim_woK0_woSp_data_sub->Clone("MMnmiss_IMnpim_woK0_woSp_data_subRebin");
  MMnmiss_IMnpim_woK0_woSp_data_subRebin->RebinX(2);
  MMnmiss_IMnpim_woK0_woSp_data_subRebin->RebinY(2);
  MMnmiss_IMnpim_woK0_woSp_data_subRebin->Draw("colz");

  //and subtract
  TCanvas *cMMnmiss_IMnpim_woK0_woSp_data_sub = new TCanvas("cMMnmiss_IMnpim_woK0_woSp_data_sub","cMMnmiss_IMnpim_woK0_woSp_data_sub",800,800);
  cMMnmiss_IMnpim_woK0_woSp_data_sub->Divide(2,2,0,0);
  cMMnmiss_IMnpim_woK0_woSp_data_sub->cd(3);
  //MMnmiss_IMnpim_woK0_woSp_data_sub->RebinX(2);
  //MMnmiss_IMnpim_woK0_woSp_data_sub->RebinY(2);
  MMnmiss_IMnpim_woK0_woSp_data_sub->SetTitle("");
  MMnmiss_IMnpim_woK0_woSp_data_sub->Draw("colz");
  cMMnmiss_IMnpim_woK0_woSp_data_sub->cd(1);
  TH1D* IMnpim_woK0_woSp_data_sub = (TH1D*)MMnmiss_IMnpim_woK0_woSp_data_sub->ProjectionX("IMnpim_woK0_woSp_data_sub");
  IMnpim_woK0_woSp_data_sub->Draw("HE");
  cMMnmiss_IMnpim_woK0_woSp_data_sub->cd(4);
  TH1D* MMnmiss_woK0_woSp_data_sub = (TH1D*)MMnmiss_IMnpim_woK0_woSp_data_sub->ProjectionY("MMnmiss_woK0_woSp_data_sub");
  TGraphErrors *gr_MMnmiss_woK0_woSp_data_sub = new TGraphErrors();
  HistToRorateGraph(MMnmiss_woK0_woSp_data_sub,*gr_MMnmiss_woK0_woSp_data_sub);
  gr_MMnmiss_woK0_woSp_data_sub->Draw("AP");
  //gr_MMnmiss_woK0_woSp_data_sub->Print("all");
   
 


  //draw vicinity subtracted
  TCanvas *cMMnmiss_IMnpim_woK0_woSp_vici_data_sub = new TCanvas("cMMnmiss_IMnpim_woK0_woSp_vici_data_sub","cMMnmiss_IMnpim_woK0_woSp_vici_data_sub",800,800);
  cMMnmiss_IMnpim_woK0_woSp_vici_data_sub->Divide(2,2,0,0);
  cMMnmiss_IMnpim_woK0_woSp_vici_data_sub->cd(3);
  TH2D* MMnmiss_IMnpim_woK0_woSp_vici_data_sub = (TH2D*)MMnmiss_IMnpim_woK0_woSp_vici_data->Clone("MMnmiss_IMnpim_woK0_woSp_vici_data_sub");
  MMnmiss_IMnpim_woK0_woSp_vici_data_sub->Add(MMnmiss_IMnpim_woK0_woSp_vici_mix,-1.0);
  MMnmiss_IMnpim_woK0_woSp_vici_data_sub->SetTitle("");
  //MMnmiss_IMnpim_woK0_woSp_vici_data_sub->GetXaxis()->SetLabelOffset(0.005);
  //MMnmiss_IMnpim_woK0_woSp_vici_data_sub->GetYaxis()->SetLabelOffset(0.005);
  MMnmiss_IMnpim_woK0_woSp_vici_data_sub->Draw("colz");
  
  cMMnmiss_IMnpim_woK0_woSp_vici_data_sub->cd(1);
  TH1D* IMnpim_woK0_woSp_vici_data_sub = (TH1D*)MMnmiss_IMnpim_woK0_woSp_vici_data_sub->ProjectionX("IMnpim_woK0_woSp_vici_data_sub");
  IMnpim_woK0_woSp_vici_data_sub->Draw("HE");
  
  
  cMMnmiss_IMnpim_woK0_woSp_vici_data_sub->cd(4);
  TH1D* MMnmiss_woK0_woSp_vici_data_sub = (TH1D*)MMnmiss_IMnpim_woK0_woSp_vici_data_sub->ProjectionY("MMnmiss_woK0_woSp_vici_data_sub");
  TGraphErrors *gr_MMnmiss_woK0_woSp_vici_data_sub = new TGraphErrors();
  HistToRorateGraph(MMnmiss_woK0_woSp_vici_data_sub,*gr_MMnmiss_woK0_woSp_vici_data_sub);
  gr_MMnmiss_woK0_woSp_vici_data_sub->Draw("AP");
  
  TH2D* IMnpim_IMnpip_n_data = (TH2D*)fr->Get("IMnpim_IMnpip_dE_n");
  TH2D* IMnpim_IMnpip_n_mix = (TH2D*)fmix->Get("IMnpim_IMnpip_dE_n");
  
  TCanvas *cIMnpim_IMnpip_n_sub = new TCanvas("cIMnpim_IMnpip_n_sub","cIMnpim_IMnpip_n_sub");
  TH2D* IMnpim_IMnpip_n_sub = (TH2D*)IMnpim_IMnpip_n_data->Clone("IMnpim_IMnpip_n_sub");
  IMnpim_IMnpip_n_sub->RebinX(2);
  IMnpim_IMnpip_n_sub->RebinY(2);
  IMnpim_IMnpip_n_mix->RebinX(2);
  IMnpim_IMnpip_n_mix->RebinY(2);
  IMnpim_IMnpip_n_sub->Add(IMnpim_IMnpip_n_mix,-1.);
  IMnpim_IMnpip_n_sub->SetTitle("IMnpim_IMnpip_n_sub");
  IMnpim_IMnpip_n_sub->Draw("colz");
  
  TH2D* IMnpim_IMnpip_wK0_n_data = (TH2D*)fr->Get("IMnpim_IMnpip_dE_wK0_n");
  TH2D* IMnpim_IMnpip_wK0_n_mix  = (TH2D*)fmix->Get("IMnpim_IMnpip_dE_wK0_n");
  TCanvas *cIMnpim_IMnpip_wK0_n_sub = new TCanvas("cIMnpim_IMnpip_wK0_n_sub","cIMnpim_IMnpip_wK0_n_sub");
  TH2D* IMnpim_IMnpip_wK0_n_sub = (TH2D*)IMnpim_IMnpip_wK0_n_data->Clone("IMnpim_IMnpip_wK0_n_sub");
  IMnpim_IMnpip_wK0_n_sub->RebinX(2);
  IMnpim_IMnpip_wK0_n_sub->RebinY(2);
  IMnpim_IMnpip_wK0_n_mix->RebinX(2);
  IMnpim_IMnpip_wK0_n_mix->RebinY(2);
  IMnpim_IMnpip_wK0_n_sub->Add(IMnpim_IMnpip_wK0_n_mix,-1.);
  IMnpim_IMnpip_wK0_n_sub->SetTitle("IMnpim_IMnpip_wK0_n_sub");
  IMnpim_IMnpip_wK0_n_sub->Draw("colz");
 
  TH2D* IMnpim_IMnpip_wK0_wSid_n_data = (TH2D*)fr->Get("IMnpim_IMnpip_dE_wK0_wSid_n");
  TH2D* IMnpim_IMnpip_wK0_wSid_n_mix  = (TH2D*)fmix->Get("IMnpim_IMnpip_dE_wK0_wSid_n");
  TCanvas *cIMnpim_IMnpip_wK0_wSid_n_sub = new TCanvas("cIMnpim_IMnpip_wK0_wSid_n_sub","cIMnpim_IMnpip_wK0_wSid_n_sub");
  TH2D* IMnpim_IMnpip_wK0_wSid_n_sub = (TH2D*)IMnpim_IMnpip_wK0_wSid_n_data->Clone("IMnpim_IMnpip_wK0_wSid_n_sub");
  IMnpim_IMnpip_wK0_wSid_n_sub->RebinX(2);
  IMnpim_IMnpip_wK0_wSid_n_sub->RebinY(2);
  IMnpim_IMnpip_wK0_wSid_n_mix->RebinX(2);
  IMnpim_IMnpip_wK0_wSid_n_mix->RebinY(2);
  IMnpim_IMnpip_wK0_wSid_n_sub->Add(IMnpim_IMnpip_wK0_wSid_n_mix,-1.);
  IMnpim_IMnpip_wK0_wSid_n_sub->SetTitle("IMnpim_IMnpip_wK0_wSid_n_sub");
  IMnpim_IMnpip_wK0_wSid_n_sub->Draw("colz");



  TH2D* q_IMnpipi_woK0_wSid_n_woSm_data = (TH2D*)fr->Get("q_IMnpipi_woK0_wSid_n_woSm");
  TH2D* q_IMnpipi_woK0_wSid_n_woSm_mix  = (TH2D*)fmix->Get("q_IMnpipi_woK0_wSid_n_woSm");
  

  TCanvas *cq_IMnpipi_woK0_wSid_n_woSm_comp= new TCanvas("cq_IMnpipi_woK0_wSid_n_woSm_comp","cq_IMnpipi_woK0_wSid_n_woSm_comp");
  cq_IMnpipi_woK0_wSid_n_woSm_comp->cd();
  TH1D* IMnpipi_woK0_wSid_n_woSm_data = (TH1D*)q_IMnpipi_woK0_wSid_n_woSm_data->ProjectionX("IMnpipi_woK0_wSid_n_woSm_data");
  IMnpipi_woK0_wSid_n_woSm_data->Draw("HE");
  TH1D* IMnpipi_woK0_wSid_n_woSm_mix = (TH1D*)q_IMnpipi_woK0_wSid_n_woSm_mix->ProjectionX("IMnpipi_woK0_wSid_n_woSm_mix");
  IMnpipi_woK0_wSid_n_woSm_mix->SetLineColor(2);
  IMnpipi_woK0_wSid_n_woSm_mix->Draw("HEsame");

  TCanvas *cq_IMnpipi_woK0_wSid_n_woSm_comp_0 = new TCanvas("cq_IMnpipi_woK0_wSid_n_woSm_comp_0","cq_IMnpipi_woK0_wSid_n_woSm_comp_0");
  cq_IMnpipi_woK0_wSid_n_woSm_comp_0->cd();
  const int bin350 = q_IMnpipi_woK0_wSid_n_woSm_data->GetYaxis()->FindBin(0.35);
  TH1D* IMnpipi_woK0_wSid_n_woSm_data_0 = (TH1D*)q_IMnpipi_woK0_wSid_n_woSm_data->ProjectionX("IMnpipi_woK0_wSid_n_woSm_data_0",1,bin350-1);
  IMnpipi_woK0_wSid_n_woSm_data_0->SetTitle("IMnpipi_woK0_wSid_n_woSm_data_0");
  IMnpipi_woK0_wSid_n_woSm_data_0->Draw("HE");
  TH1D* IMnpipi_woK0_wSid_n_woSm_mix_0 = (TH1D*)q_IMnpipi_woK0_wSid_n_woSm_mix->ProjectionX("IMnpipi_woK0_wSid_n_woSm_mix_0",1,bin350-1);
  IMnpipi_woK0_wSid_n_woSm_mix_0->SetLineColor(2);
  IMnpipi_woK0_wSid_n_woSm_mix_0->Draw("HEsame");
  
  TCanvas *cq_IMnpipi_woK0_wSid_n_woSm_sub_0 = new TCanvas("cq_IMnpipi_woK0_wSid_n_woSm_sub_0","cq_IMnpipi_woK0_wSid_n_woSm_sub_0");
  TH1D* IMnpipi_woK0_wSid_n_woSm_0_sub = (TH1D*)IMnpipi_woK0_wSid_n_woSm_data_0->Clone("IMnpipi_woK0_wSid_n_woSm_0_sub");
  IMnpipi_woK0_wSid_n_woSm_0_sub->Add(IMnpipi_woK0_wSid_n_woSm_mix_0,-1);
  IMnpipi_woK0_wSid_n_woSm_0_sub->Draw("HE");

  TCanvas *cq_IMnpipi_woK0_wSid_n_woSm_comp_350 = new TCanvas("cq_IMnpipi_woK0_wSid_n_woSm_comp_350","cq_IMnpipi_woK0_wSid_n_woSm_comp_350");
  cq_IMnpipi_woK0_wSid_n_woSm_comp_350->cd();
  TH1D* IMnpipi_woK0_wSid_n_woSm_data_350 = (TH1D*)q_IMnpipi_woK0_wSid_n_woSm_data->ProjectionX("IMnpipi_woK0_wSid_n_woSm_data_350",bin350,500);
  IMnpipi_woK0_wSid_n_woSm_data_350->SetTitle("IMnpipi_woK0_wSid_n_woSm_data_350");
  IMnpipi_woK0_wSid_n_woSm_data_350->Draw("HE");
  TH1D* IMnpipi_woK0_wSid_n_woSm_mix_350 = (TH1D*)q_IMnpipi_woK0_wSid_n_woSm_mix->ProjectionX("IMnpipi_woK0_wSid_n_woSm_mix_350",bin350,500);
  IMnpipi_woK0_wSid_n_woSm_mix_350->SetLineColor(2);
  IMnpipi_woK0_wSid_n_woSm_mix_350->Draw("HEsame");
  
  TCanvas *cq_IMnpipi_woK0_wSid_n_woSm_sub_350 = new TCanvas("cq_IMnpipi_woK0_wSid_n_woSm_sub_350","cq_IMnpipi_woK0_wSid_n_woSm_sub_350");
  TH1D* IMnpipi_woK0_wSid_n_woSm_350_sub = (TH1D*)IMnpipi_woK0_wSid_n_woSm_data_350->Clone("IMnpipi_woK0_wSid_n_woSm_350_sub");
  IMnpipi_woK0_wSid_n_woSm_350_sub->Add(IMnpipi_woK0_wSid_n_woSm_mix_350,-1);
  IMnpipi_woK0_wSid_n_woSm_350_sub->Draw("HE");

  TH2D* q_IMnpipi_woK0_wSid_n_woSp_data = (TH2D*)fr->Get("q_IMnpipi_woK0_wSid_n_woSp");
  TH2D* q_IMnpipi_woK0_wSid_n_woSp_mix  = (TH2D*)fmix->Get("q_IMnpipi_woK0_wSid_n_woSp");
  
  TCanvas *cq_IMnpipi_woK0_wSid_n_woSp_comp_0 = new TCanvas("cq_IMnpipi_woK0_wSid_n_woSp_comp_0","cq_IMnpipi_woK0_wSid_n_woSp_comp_0");
  cq_IMnpipi_woK0_wSid_n_woSp_comp_0->cd();
  TH1D* IMnpipi_woK0_wSid_n_woSp_data_0 = (TH1D*)q_IMnpipi_woK0_wSid_n_woSp_data->ProjectionX("IMnpipi_woK0_wSid_n_woSp_data_0",1,bin350-1);
  IMnpipi_woK0_wSid_n_woSp_data_0->SetTitle("IMnpipi_woK0_wSid_n_woSp_data_0");
  IMnpipi_woK0_wSid_n_woSp_data_0->Draw("HE");
  TH1D* IMnpipi_woK0_wSid_n_woSp_mix_0 = (TH1D*)q_IMnpipi_woK0_wSid_n_woSp_mix->ProjectionX("IMnpipi_woK0_wSid_n_woSp_mix_0",1,bin350-1);
  IMnpipi_woK0_wSid_n_woSp_mix_0->SetLineColor(2);
  IMnpipi_woK0_wSid_n_woSp_mix_0->Draw("HEsame");

  TCanvas *cq_IMnpipi_woK0_wSid_n_woSp_sub_0 = new TCanvas("cq_IMnpipi_woK0_wSid_n_woSp_sub_0","cq_IMnpipi_woK0_wSid_n_woSp_sub_0");
  TH1D* IMnpipi_woK0_wSid_n_woSp_0_sub = (TH1D*)IMnpipi_woK0_wSid_n_woSp_data_0->Clone("IMnpipi_woK0_wSid_n_woSp_0_sub");
  IMnpipi_woK0_wSid_n_woSp_0_sub->Add(IMnpipi_woK0_wSid_n_woSp_mix_0,-1);
  IMnpipi_woK0_wSid_n_woSp_0_sub->Draw("HE");

  TCanvas *cq_IMnpipi_woK0_wSid_n_woSp_comp_350 = new TCanvas("cq_IMnpipi_woK0_wSid_n_woSp_comp_350","cq_IMnpipi_woK0_wSid_n_woSp_comp_350");
  cq_IMnpipi_woK0_wSid_n_woSp_comp_350->cd();
  TH1D* IMnpipi_woK0_wSid_n_woSp_data_350 = (TH1D*)q_IMnpipi_woK0_wSid_n_woSp_data->ProjectionX("IMnpipi_woK0_wSid_n_woSp_data_350",bin350,500);
  IMnpipi_woK0_wSid_n_woSp_data_350->SetTitle("IMnpipi_woK0_wSid_n_woSp_data_350");
  IMnpipi_woK0_wSid_n_woSp_data_350->Draw("HE");
  TH1D* IMnpipi_woK0_wSid_n_woSp_mix_350 = (TH1D*)q_IMnpipi_woK0_wSid_n_woSp_mix->ProjectionX("IMnpipi_woK0_wSid_n_woSp_mix_350",bin350,500);
  IMnpipi_woK0_wSid_n_woSp_mix_350->SetLineColor(2);
  IMnpipi_woK0_wSid_n_woSp_mix_350->Draw("HEsame");

  TCanvas *cq_IMnpipi_woK0_wSid_n_woSp_sub_350 = new TCanvas("cq_IMnpipi_woK0_wSid_n_woSp_sub_350","cq_IMnpipi_woK0_wSid_n_woSp_sub_350");
  TH1D* IMnpipi_woK0_wSid_n_woSp_350_sub = (TH1D*)IMnpipi_woK0_wSid_n_woSp_data_350->Clone("IMnpipi_woK0_wSid_n_woSp_350_sub");
  IMnpipi_woK0_wSid_n_woSp_350_sub->Add(IMnpipi_woK0_wSid_n_woSp_mix_350,-1);
  IMnpipi_woK0_wSid_n_woSp_350_sub->Draw("HE");



  
  TH2D* q_IMnpipi_wSid_n_Sp_data = (TH2D*)fr->Get("q_IMnpipi_wSid_n_Sp");
  TH2D* q_IMnpipi_wSid_n_Sp_mix  = (TH2D*)fmix->Get("q_IMnpipi_wSid_n_Sp");

  TCanvas *cq_IMnpipi_wSid_n_Sp_comp= new TCanvas("cq_IMnpipi_wSid_n_Sp_comp","cq_IMnpipi_wSid_n_Sp_comp");
  cq_IMnpipi_wSid_n_Sp_comp->cd();
  TH1D* IMnpipi_wSid_n_Sp_data = (TH1D*)q_IMnpipi_wSid_n_Sp_data->ProjectionX("IMnpipi_wSid_n_Sp_data");
  IMnpipi_wSid_n_Sp_data->Draw("HE");
  TH1D* IMnpipi_wSid_n_Sp_mix = (TH1D*)q_IMnpipi_wSid_n_Sp_mix->ProjectionX("IMnpipi_wSid_n_Sp_mix");
  IMnpipi_wSid_n_Sp_mix->SetLineColor(2);
  IMnpipi_wSid_n_Sp_mix->Draw("HEsame");

  TCanvas *cq_IMnpipi_wSid_n_Sp_comp_0 = new TCanvas("cq_IMnpipi_wSid_n_Sp_comp_0","cq_IMnpipi_wSid_n_Sp_comp_0");
  cq_IMnpipi_wSid_n_Sp_comp_0->cd();
  //const int bin350 = q_IMnpipi_wSid_n_Sp_data->GetYaxis()->FindBin(0.35);
  TH1D* IMnpipi_wSid_n_Sp_data_0 = (TH1D*)q_IMnpipi_wSid_n_Sp_data->ProjectionX("IMnpipi_wSid_n_Sp_data_0",1,bin350-1);
  IMnpipi_wSid_n_Sp_data_0->SetTitle("IMnpipi_wSid_n_Sp_data_0");
  IMnpipi_wSid_n_Sp_data_0->Draw("HE");
  TH1D* IMnpipi_wSid_n_Sp_mix_0 = (TH1D*)q_IMnpipi_wSid_n_Sp_mix->ProjectionX("IMnpipi_wSid_n_Sp_mix_0",1,bin350-1);
  IMnpipi_wSid_n_Sp_mix_0->SetLineColor(2);
  IMnpipi_wSid_n_Sp_mix_0->Draw("HEsame");


  TCanvas *cq_IMnpipi_wSid_n_Sp_comp_350 = new TCanvas("cq_IMnpipi_wSid_n_Sp_comp_350","cq_IMnpipi_wSid_n_Sp_comp_350");
  cq_IMnpipi_wSid_n_Sp_comp_350->cd();
  TH1D* IMnpipi_wSid_n_Sp_data_350 = (TH1D*)q_IMnpipi_wSid_n_Sp_data->ProjectionX("IMnpipi_wSid_n_Sp_data_350",bin350,500);
  IMnpipi_wSid_n_Sp_data_350->SetTitle("IMnpipi_wSid_n_Sp_data_350");
  IMnpipi_wSid_n_Sp_data_350->Draw("HE");
  TH1D* IMnpipi_wSid_n_Sp_mix_350 = (TH1D*)q_IMnpipi_wSid_n_Sp_mix->ProjectionX("IMnpipi_wSid_n_Sp_mix_350",bin350,500);
  IMnpipi_wSid_n_Sp_mix_350->SetLineColor(2);
  IMnpipi_wSid_n_Sp_mix_350->Draw("HEsame");

  TCanvas *cq_IMnpipi_wSid_n_Sp_sub = new TCanvas("cq_IMnpipi_wSid_n_Sp_sub","cq_IMnpipi_wSid_n_Sp_sub");
  TH2D *q_IMnpipi_wSid_n_Sp_sub = (TH2D*)q_IMnpipi_wSid_n_Sp_data->Clone("q_IMnpipi_wSid_n_Sp_sub");
  q_IMnpipi_wSid_n_Sp_sub->RebinX(2);
  q_IMnpipi_wSid_n_Sp_mix->RebinX(2);
  q_IMnpipi_wSid_n_Sp_sub->Add(q_IMnpipi_wSid_n_Sp_mix,-1);
  q_IMnpipi_wSid_n_Sp_sub->Draw("colz");

  
  TCanvas *cIMnpipi_wSid_n_Sp_comp_0 = new TCanvas("cIMnpipi_wSid_n_Sp_sub_0","cIMnpipi_wSid_n_Sp_sub_0");
  cIMnpipi_wSid_n_Sp_comp_0->cd();
  TH1D* IMnpipi_wSid_n_Sp_data_0_sub = (TH1D*) IMnpipi_wSid_n_Sp_data_0->Clone("IMnpipi_wSid_n_Sp_data_0_sub");
  IMnpipi_wSid_n_Sp_data_0_sub->Add(IMnpipi_wSid_n_Sp_mix_0,-1);
  IMnpipi_wSid_n_Sp_data_0_sub->Draw("HE");

  TCanvas *cIMnpipi_wSid_n_Sp_comp_350 = new TCanvas("cIMnpipi_wSid_n_Sp_sub_350","cIMnpipi_wSid_n_Sp_sub_350");
  cIMnpipi_wSid_n_Sp_comp_350->cd();
  TH1D* IMnpipi_wSid_n_Sp_data_350_sub = (TH1D*) IMnpipi_wSid_n_Sp_data_350->Clone("IMnpipi_wSid_n_Sp_data_350_sub");
  IMnpipi_wSid_n_Sp_data_350_sub->Add(IMnpipi_wSid_n_Sp_mix_350,-1);
  IMnpipi_wSid_n_Sp_data_350_sub->Draw("HE");


  TH2D* q_IMnpipi_wSid_n_Sm_data = (TH2D*)fr->Get("q_IMnpipi_wSid_n_Sm");
  TH2D* q_IMnpipi_wSid_n_Sm_mix  = (TH2D*)fmix->Get("q_IMnpipi_wSid_n_Sm");
  
  TCanvas *cq_IMnpipi_wSid_n_Sm_comp_0 = new TCanvas("cq_IMnpipi_wSid_n_Sm_comp_0","cq_IMnpipi_wSid_n_Sm_comp_0");
  cq_IMnpipi_wSid_n_Sm_comp_0->cd();
  TH1D* IMnpipi_wSid_n_Sm_data_0 = (TH1D*)q_IMnpipi_wSid_n_Sm_data->ProjectionX("IMnpipi_wSid_n_Sm_data_0",1,bin350-1);
  IMnpipi_wSid_n_Sm_data_0->SetTitle("IMnpipi_wSid_n_Sm_data_0");
  IMnpipi_wSid_n_Sm_data_0->Draw("HE");
  TH1D* IMnpipi_wSid_n_Sm_mix_0 = (TH1D*)q_IMnpipi_wSid_n_Sm_mix->ProjectionX("IMnpipi_wSid_n_Sm_mix_0",1,bin350-1);
  IMnpipi_wSid_n_Sm_mix_0->SetLineColor(2);
  IMnpipi_wSid_n_Sm_mix_0->Draw("HEsame");


  TCanvas *cq_IMnpipi_wSid_n_Sm_comp_350 = new TCanvas("cq_IMnpipi_wSid_n_Sm_comp_350","cq_IMnpipi_wSid_n_Sm_comp_350");
  cq_IMnpipi_wSid_n_Sm_comp_350->cd();
  TH1D* IMnpipi_wSid_n_Sm_data_350 = (TH1D*)q_IMnpipi_wSid_n_Sm_data->ProjectionX("IMnpipi_wSid_n_Sm_data_350",bin350,500);
  IMnpipi_wSid_n_Sm_data_350->SetTitle("IMnpipi_wSid_n_Sm_data_350");
  IMnpipi_wSid_n_Sm_data_350->Draw("HE");
  TH1D* IMnpipi_wSid_n_Sm_mix_350 = (TH1D*)q_IMnpipi_wSid_n_Sm_mix->ProjectionX("IMnpipi_wSid_n_Sm_mix_350",bin350,500);
  IMnpipi_wSid_n_Sm_mix_350->SetLineColor(2);
  IMnpipi_wSid_n_Sm_mix_350->Draw("HEsame");

  TCanvas *cq_IMnpipi_wSid_n_Sm_sub = new TCanvas("cq_IMnpipi_wSid_n_Sm_sub","cq_IMnpipi_wSid_n_Sm_sub");
  TH2D *q_IMnpipi_wSid_n_Sm_sub = (TH2D*)q_IMnpipi_wSid_n_Sm_data->Clone("q_IMnpipi_wSid_n_Sm_sub");
  q_IMnpipi_wSid_n_Sm_sub->RebinX(2);
  q_IMnpipi_wSid_n_Sm_mix->RebinX(2);
  q_IMnpipi_wSid_n_Sm_sub->Add(q_IMnpipi_wSid_n_Sm_mix,-1);
  q_IMnpipi_wSid_n_Sm_sub->Draw("colz");

  
  TCanvas *cIMnpipi_wSid_n_Sm_comp_0 = new TCanvas("cIMnpipi_wSid_n_Sm_sub_0","cIMnpipi_wSid_n_Sm_sub_0");
  cIMnpipi_wSid_n_Sm_comp_0->cd();
  TH1D* IMnpipi_wSid_n_Sm_data_0_sub = (TH1D*) IMnpipi_wSid_n_Sm_data_0->Clone("IMnpipi_wSid_n_Sm_data_0_sub");
  IMnpipi_wSid_n_Sm_data_0_sub->Add(IMnpipi_wSid_n_Sm_mix_0,-1);
  IMnpipi_wSid_n_Sm_data_0_sub->Draw("HE");

  TCanvas *cIMnpipi_wSid_n_Sm_comp_350 = new TCanvas("cIMnpipi_wSid_n_Sm_sub_350","cIMnpipi_wSid_n_Sm_sub_350");
  cIMnpipi_wSid_n_Sm_comp_350->cd();
  TH1D* IMnpipi_wSid_n_Sm_data_350_sub = (TH1D*) IMnpipi_wSid_n_Sm_data_350->Clone("IMnpipi_wSid_n_Sm_data_350_sub");
  IMnpipi_wSid_n_Sm_data_350_sub->Add(IMnpipi_wSid_n_Sm_mix_350,-1);
  IMnpipi_wSid_n_Sm_data_350_sub->Draw("HE");


  TH2D* IMpippim_IMnpip_vici_woSm_data = (TH2D*)fr->Get("IMpippim_IMnpip_vici_woSm");
  TH2D* IMpippim_IMnpim_vici_woSp_data = (TH2D*)fr->Get("IMpippim_IMnpim_vici_woSp");
  TH2D* IMpippim_IMnpip_vici_woSm_mix = (TH2D*)fmix->Get("IMpippim_IMnpip_vici_woSm");
  TH2D* IMpippim_IMnpim_vici_woSp_mix = (TH2D*)fmix->Get("IMpippim_IMnpim_vici_woSp");
  
  /*
  TCanvas *cIMpippim_IMnpip_vici_woSm_data = new TCanvas("cIMpippim_IMnpip_vici_woSm_data","cIMpippim_IMnpip_vici_woSm_data");
  cIMpippim_IMnpip_vici_woSm_data->cd();
  IMpippim_IMnpip_vici_woSm_data->Draw("colz");

  TCanvas *cIMpippim_IMnpim_vici_woSp_data = new TCanvas("cIMpippim_IMnpim_vici_woSp_data","cIMpippim_IMnpim_vici_woSp_data");
  cIMpippim_IMnpim_vici_woSp_data->cd();
  IMpippim_IMnpim_vici_woSp_data->Draw("colz");

  TCanvas *cIMpippim_IMnpip_vici_woSm_mix = new TCanvas("cIMpippim_IMnpip_vici_woSm_mix","cIMpippim_IMnpip_vici_woSm_mix");
  cIMpippim_IMnpip_vici_woSm_mix->cd();
  IMpippim_IMnpip_vici_woSm_mix->Draw("colz");

  TCanvas *cIMpippim_IMnpim_vici_woSp_mix = new TCanvas("cIMpippim_IMnpim_vici_woSp_mix","cIMpippim_IMnpim_vici_woSp_mix");
  cIMpippim_IMnpim_vici_woSp_mix->cd();
  IMpippim_IMnpim_vici_woSp_mix->Draw("colz");
  
  TCanvas *cIMpippim_vici_woSp_comp = new TCanvas("cIMpippim_vici_woSp_comp","cIMpippim_vici_woSp_comp");
  cIMpippim_vici_woSp_comp->cd();
  TH1D* IMpippim_vici_woSp_data = (TH1D*)IMpippim_IMnpim_vici_woSp_data->ProjectionY("IMpippim_vici_woSp_data");
  TH1D* IMpippim_vici_woSp_mix  = (TH1D*)IMpippim_IMnpim_vici_woSp_mix->ProjectionY("IMpippim_vici_woSp_mix");
  TH1D* IMpippim_vici_woSp_comp = (TH1D*)IMpippim_vici_woSp_data->Clone("IMpippim_vici_woSp_comp");
  IMpippim_vici_woSp_comp->Rebin(8);
  IMpippim_vici_woSp_mix->Rebin(8);
  IMpippim_vici_woSp_comp->Add(IMpippim_vici_woSp_mix,-1);
  IMpippim_vici_woSp_comp->Draw("HE");
  
  TCanvas *cIMpippim_vici_woSm_comp = new TCanvas("cIMpippim_vici_woSm_comp","cIMpippim_vici_woSm_comp");
  cIMpippim_vici_woSm_comp->cd();
  TH1D* IMpippim_vici_woSm_data = (TH1D*)IMpippim_IMnpip_vici_woSm_data->ProjectionY("IMpippim_vici_woSm_data");
  TH1D* IMpippim_vici_woSm_mix  = (TH1D*)IMpippim_IMnpip_vici_woSm_mix->ProjectionY("IMpippim_vici_woSm_mix");
  TH1D* IMpippim_vici_woSm_comp = (TH1D*)IMpippim_vici_woSm_data->Clone("IMpippim_vici_woSm_comp");
  IMpippim_vici_woSm_comp->Rebin(8);
  IMpippim_vici_woSm_mix->Rebin(8);
  IMpippim_vici_woSm_comp->Add(IMpippim_vici_woSm_mix,-1);
  IMpippim_vici_woSm_comp->Draw("HE");
  */
  
  TH2D* IMpippim_IMnpip_wSid_n_woSm_data = (TH2D*)fr->Get("IMpippim_IMnpip_wSid_n_woSm");
  TH2D* IMpippim_IMnpim_wSid_n_woSp_data = (TH2D*)fr->Get("IMpippim_IMnpim_wSid_n_woSp");
  TH2D* IMpippim_IMnpip_wSid_n_woSm_mix = (TH2D*)fmix->Get("IMpippim_IMnpip_wSid_n_woSm");
  TH2D* IMpippim_IMnpim_wSid_n_woSp_mix = (TH2D*)fmix->Get("IMpippim_IMnpim_wSid_n_woSp");

  TCanvas *cIMpippim_IMnpip_wSid_n_woSm_data = new TCanvas("cIMpippim_IMnpip_wSid_n_woSm_data","cIMpippim_IMnpip_wSid_n_woSm_data");
  cIMpippim_IMnpip_wSid_n_woSm_data->cd();
  IMpippim_IMnpip_wSid_n_woSm_data->Draw("colz");

  TCanvas *cIMpippim_IMnpim_wSid_n_woSp_data = new TCanvas("cIMpippim_IMnpim_wSid_n_woSp_data","cIMpippim_IMnpim_wSid_n_woSp_data");
  cIMpippim_IMnpim_wSid_n_woSp_data->cd();
  IMpippim_IMnpim_wSid_n_woSp_data->Draw("colz");

  TCanvas *cIMpippim_IMnpip_wSid_n_woSm_mix = new TCanvas("cIMpippim_IMnpip_wSid_n_woSm_mix","cIMpippim_IMnpip_wSid_n_woSm_mix");
  cIMpippim_IMnpip_wSid_n_woSm_mix->cd();
  IMpippim_IMnpip_wSid_n_woSm_mix->Draw("colz");

  TCanvas *cIMpippim_IMnpim_wSid_n_woSp_mix = new TCanvas("cIMpippim_IMnpim_wSid_n_woSp_mix","cIMpippim_IMnpim_wSid_n_woSp_mix");
  cIMpippim_IMnpim_wSid_n_woSp_mix->cd();
  IMpippim_IMnpim_wSid_n_woSp_mix->Draw("colz");
  
  TCanvas *cIMpippim_wSid_n_woSp_comp = new TCanvas("cIMpippim_wSid_n_woSp_comp","cIMpippim_wSid_n_woSp_comp");
  cIMpippim_wSid_n_woSp_comp->cd();
  TH1D* IMpippim_wSid_n_woSp_data = (TH1D*)IMpippim_IMnpim_wSid_n_woSp_data->ProjectionY("IMpippim_wSid_n_woSp_data");
  TH1D* IMpippim_wSid_n_woSp_mix  = (TH1D*)IMpippim_IMnpim_wSid_n_woSp_mix->ProjectionY("IMpippim_wSid_n_woSp_mix");
  TH1D* IMpippim_wSid_n_woSp_comp = (TH1D*)IMpippim_wSid_n_woSp_data->Clone("IMpippim_wSid_n_woSp_comp");
  IMpippim_wSid_n_woSp_comp->Rebin(8);
  IMpippim_wSid_n_woSp_mix->Rebin(8);
  IMpippim_wSid_n_woSp_comp->Add(IMpippim_wSid_n_woSp_mix,-1);
  IMpippim_wSid_n_woSp_comp->Draw("HE");
  
  TCanvas *cIMpippim_wSid_n_woSm_comp = new TCanvas("cIMpippim_wSid_n_woSm_comp","cIMpippim_wSid_n_woSm_comp");
  cIMpippim_wSid_n_woSm_comp->cd();
  TH1D* IMpippim_wSid_n_woSm_data = (TH1D*)IMpippim_IMnpip_wSid_n_woSm_data->ProjectionY("IMpippim_wSid_n_woSm_data");
  TH1D* IMpippim_wSid_n_woSm_mix  = (TH1D*)IMpippim_IMnpip_wSid_n_woSm_mix->ProjectionY("IMpippim_wSid_n_woSm_mix");
  TH1D* IMpippim_wSid_n_woSm_comp = (TH1D*)IMpippim_wSid_n_woSm_data->Clone("IMpippim_wSid_n_woSm_comp");
  IMpippim_wSid_n_woSm_comp->Rebin(8);
  IMpippim_wSid_n_woSm_mix->Rebin(8);
  IMpippim_wSid_n_woSm_comp->Add(IMpippim_wSid_n_woSm_mix,-1);
  IMpippim_wSid_n_woSm_comp->Draw("HE");
  


  TCanvas *cIMpippim_IMnpipi_n_wSid_woSp = new TCanvas("cIMpippim_IMnpipi_n_wSid_woSp","cIMpippim_IMnpipi_n_wSid_woSp");
  TH2D* IMpippim_IMnpipi_n_wSid_woSp_data = (TH2D*)fr->Get("IMpippim_IMnpipi_n_wSid_woSp");
  TH2D* IMpippim_IMnpipi_n_wSid_woSp_mix  = (TH2D*)fmix->Get("IMpippim_IMnpipi_n_wSid_woSp");
  IMpippim_IMnpipi_n_wSid_woSp_data->RebinX(2);
  IMpippim_IMnpipi_n_wSid_woSp_mix->RebinX(2);
  IMpippim_IMnpipi_n_wSid_woSp_data->RebinY(10);
  IMpippim_IMnpipi_n_wSid_woSp_mix->RebinY(10);
  TH2D* IMpippim_IMnpipi_n_wSid_woSp_sub  = (TH2D*)IMpippim_IMnpipi_n_wSid_woSp_data->Clone("IMpippim_IMnpipi_n_wSid_woSp_sub");
  IMpippim_IMnpipi_n_wSid_woSp_sub->Add(IMpippim_IMnpipi_n_wSid_woSp_mix,-1);
  IMpippim_IMnpipi_n_wSid_woSp_sub->Draw("colz");

  TCanvas *cIMpippim_IMnpipi_n_wSid_woSm = new TCanvas("cIMpippim_IMnpipi_n_wSid_woSm","cIMpippim_IMnpipi_n_wSid_woSm");
  TH2D* IMpippim_IMnpipi_n_wSid_woSm_data = (TH2D*)fr->Get("IMpippim_IMnpipi_n_wSid_woSm");
  TH2D* IMpippim_IMnpipi_n_wSid_woSm_mix  = (TH2D*)fmix->Get("IMpippim_IMnpipi_n_wSid_woSm");
  IMpippim_IMnpipi_n_wSid_woSm_data->RebinX(2);
  IMpippim_IMnpipi_n_wSid_woSm_mix->RebinX(2);
  IMpippim_IMnpipi_n_wSid_woSm_data->RebinY(10);
  IMpippim_IMnpipi_n_wSid_woSm_mix->RebinY(10);
  TH2D* IMpippim_IMnpipi_n_wSid_woSm_sub  = (TH2D*)IMpippim_IMnpipi_n_wSid_woSm_data->Clone("IMpippim_IMnpipi_n_wSid_woSm_sub");
  IMpippim_IMnpipi_n_wSid_woSm_sub->Add(IMpippim_IMnpipi_n_wSid_woSm_mix,-1);
  IMpippim_IMnpipi_n_wSid_woSm_sub->Draw("colz");
  
  TCanvas *cIMpippim_n_wSid_woSm1450 = new TCanvas("cIMpippim_n_wSid_woSm1450","cIMpippim_n_wSid_woSm1450");
  const int bin1450 = IMpippim_IMnpipi_n_wSid_woSm_sub->GetXaxis()->FindBin(1.450);
  TH1D* IMpippim_n_wSid_woSm1450 = (TH1D*)IMpippim_IMnpipi_n_wSid_woSm_sub->ProjectionY("IMpippim_n_wSid_woSm1450",bin1450,bin1450);
  IMpippim_n_wSid_woSm1450->SetTitle("IMpippim_n_wSid_woSm1450");
  IMpippim_n_wSid_woSm1450->Draw("HE");

  TCanvas *cIMpippim_n_wSid_woSm1520 = new TCanvas("cIMpippim_n_wSid_woSm1520","cIMpippim_n_wSid_woSm1520");
  const int bin1520 = IMpippim_IMnpipi_n_wSid_woSm_sub->GetXaxis()->FindBin(1.520);
  TH1D* IMpippim_n_wSid_woSm1520 = (TH1D*)IMpippim_IMnpipi_n_wSid_woSm_sub->ProjectionY("IMpippim_n_wSid_woSm1520",bin1520,bin1520);
  IMpippim_n_wSid_woSm1520->SetTitle("IMpippim_n_wSid_woSm1520");
  IMpippim_n_wSid_woSm1520->Draw("HE");

  
  TH2D* q_IMnpipi_wK0_wSid_n_data = (TH2D*)fr->Get("q_IMnpipi_wK0_wSid_n");
  TH2D* q_IMnpipi_wK0_wSid_n_mix  = (TH2D*)fmix->Get("q_IMnpipi_wK0_wSid_n");
  TCanvas *cq_IMnpipi_wK0_wSid_n_sub = new TCanvas("cq_IMnpipi_wK0_wSid_n_sub","cq_IMnpipi_wK0_wSid_n_sub");
  TH2D* q_IMnpipi_wK0_wSid_n_sub = (TH2D*)q_IMnpipi_wK0_wSid_n_data->Clone("q_IMnpipi_wK0_wSid_n_sub");
  q_IMnpipi_wK0_wSid_n_sub->Add(q_IMnpipi_wK0_wSid_n_mix,-1);
  q_IMnpipi_wK0_wSid_n_sub->Draw("colz");

  TCanvas *cIMnpipi_wK0_wSid_n_sub_0 = new TCanvas("cIMnpipi_wK0_wSid_n_sub_0","cIMnpipi_wK0_wSid_n_sub_0");
  TH1D* IMnpipi_wK0_wSid_n_sub_0 = (TH1D*)q_IMnpipi_wK0_wSid_n_sub->ProjectionX("IMnpipi_wK0_wSid_n_sub_0",1,bin350-1);
  IMnpipi_wK0_wSid_n_sub_0->SetTitle("IMnpipi_wK0_wSid_n_sub_0");
  IMnpipi_wK0_wSid_n_sub_0->Draw("HE");

  TCanvas *cIMnpipi_wK0_wSid_n_sub_350 = new TCanvas("cIMnpipi_wK0_wSid_n_sub_350","cIMnpipi_wK0_wSid_n_sub_350");
  TH1D* IMnpipi_wK0_wSid_n_sub_350 = (TH1D*)q_IMnpipi_wK0_wSid_n_sub->ProjectionX("IMnpipi_wK0_wSid_n_sub_350",bin350,5000);
  IMnpipi_wK0_wSid_n_sub_350->SetTitle("IMnpipi_wK0_wSid_n_sub_350");
  IMnpipi_wK0_wSid_n_sub_350->Draw("HE");
  
  TH2D* q_IMnpipi_wSid_n_data = (TH2D*)fr->Get("q_IMnpipi_wSid_n");
  TH2D* q_IMnpipi_wSid_n_mix  = (TH2D*)fmix->Get("q_IMnpipi_wSid_n");
  //TH2D* q_IMnpipi_wSid_n_sub = 


  
  TCanvas *c = NULL;
  TSeqCollection *SCol = gROOT->GetListOfCanvases();
  int size = SCol->GetSize();
  TIter next(SCol);
  TString pdfname = "disp_mixcomp.pdf";
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
  



  return;
  

}
