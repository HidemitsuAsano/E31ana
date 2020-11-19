const bool showBG = true;
#include "anacuts.h"

void disp_2Dcompmix()
{
  TFile *fr = TFile::Open("evanaIMpisigma_npippim_v202_out_iso.root");
  TFile *fmix = TFile::Open("evanaIMpisigma_npippim_v202_MIX_cut3_out_iso.root");
  std::cout << fr->Print() << std::endl;
  std::cout << fmix->Print() << std::endl;

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TH1::SetDefaultSumw2();
  //const char opt[10]="cont3";
  //const char opt[10]="cont1";
  //const char opt[10]="colz";
  gStyle->SetErrorX(0.);  

  TH2D* MMnmiss_IMnpip_dE_woK0_woSm_data = (TH2D*)fr->Get("MMnmiss_IMnpip_dE_woK0_woSm");
  TCanvas *cMMnmiss_IMnpip_dE_woK0_woSm_data = new TCanvas("cMMnmiss_IMnpip_dE_woK0_woSm_data","cMMnmiss_IMnpip_dE_woK0_woSm_data",800,800);
  MMnmiss_IMnpip_dE_woK0_woSm_data->Draw("colz");


  TH2D* MMnmiss_IMnpim_dE_woK0_woSp_data = (TH2D*)fr->Get("MMnmiss_IMnpim_dE_woK0_woSp");
  TCanvas *cMMnmiss_IMnpim_dE_woK0_woSp_data = new TCanvas("cMMnmiss_IMnpim_dE_woK0_woSp_data","cMMnmiss_IMnpim_dE_woK0_woSp_data",800,800);
  MMnmiss_IMnpim_dE_woK0_woSp_data->Draw("colz");

  TH2D* MMnmiss_IMnpip_dE_woK0_woSm_mix = (TH2D*)fmix->Get("MMnmiss_IMnpip_dE_woK0_woSm");
  TCanvas *cMMnmiss_IMnpip_dE_woK0_woSm_mix = new TCanvas("cMMnmiss_IMnpip_dE_woK0_woSm_mix","cMMnmiss_IMnpip_dE_woK0_woSm_mix",800,800);
  MMnmiss_IMnpip_dE_woK0_woSm_mix->Draw("colz");
  
  TH2D* MMnmiss_IMnpim_dE_woK0_woSp_mix = (TH2D*)fmix->Get("MMnmiss_IMnpim_dE_woK0_woSp");
  TCanvas *cMMnmiss_IMnpim_dE_woK0_woSp_mix = new TCanvas("cMMnmiss_IMnpim_dE_woK0_woSp_mix","cMMnmiss_IMnpim_dE_woK0_woSp_mix",800,800);
  MMnmiss_IMnpim_dE_woK0_woSp_mix->Draw("colz");
  
  TH2D* IMnpim_IMnpip_dE_woK0_woSp_vici_data = (TH2D*)fr->Get("IMnpim_IMnpip_dE_woK0_woSp_vici");
  TCanvas *cIMnpim_IMnpip_dE_woK0_woSp_vici_data = new TCanvas("cIMnpim_IMnpip_dE_woK0_woSp_vici_data","cIMnpim_IMnpip_dE_woK0_woSp_vici_data",800,800);
  IMnpim_IMnpip_dE_woK0_woSp_vici_data->Draw("colz");


  TH2D* IMnpim_IMnpip_dE_woK0_woSm_vici_data = (TH2D*)fr->Get("IMnpim_IMnpip_dE_woK0_woSm_vici");
  TCanvas *cIMnpim_IMnpip_dE_woK0_woSm_vici_data = new TCanvas("cIMnpim_IMnpip_dE_woK0_woSm_vici_data","cIMnpim_IMnpip_dE_woK0_woSm_vici_data",800,800);
  IMnpim_IMnpip_dE_woK0_woSm_vici_data->Draw("colz");



  TH2D* MMnmiss_IMnpip_dE_woK0_woSm_vici_data = (TH2D*)fr->Get("MMnmiss_IMnpip_dE_woK0_woSm_vici");
  TH2D* MMnmiss_IMnpim_dE_woK0_woSp_vici_data = (TH2D*)fr->Get("MMnmiss_IMnpim_dE_woK0_woSp_vici");
  TH2D* MMnmiss_IMnpip_dE_woK0_woSm_viciext_data = (TH2D*)fr->Get("MMnmiss_IMnpip_dE_woK0_woSm_viciext");
  TH2D* MMnmiss_IMnpim_dE_woK0_woSp_viciext_data = (TH2D*)fr->Get("MMnmiss_IMnpim_dE_woK0_woSp_viciext");
  TH2D* MMnmiss_IMnpip_dE_woK0_woSm_vici_mix = (TH2D*)fmix->Get("MMnmiss_IMnpip_dE_woK0_woSm_vici");
  TH2D* MMnmiss_IMnpim_dE_woK0_woSp_vici_mix = (TH2D*)fmix->Get("MMnmiss_IMnpim_dE_woK0_woSp_vici");
  TH2D* MMnmiss_IMnpip_dE_woK0_woSm_viciext_mix = (TH2D*)fmix->Get("MMnmiss_IMnpip_dE_woK0_woSm_viciext");
  TH2D* MMnmiss_IMnpim_dE_woK0_woSp_viciext_mix = (TH2D*)fmix->Get("MMnmiss_IMnpim_dE_woK0_woSp_viciext");
  
  const double mixScale = 0.0496125;

  TCanvas *cMMnmiss_IMnpip_dE_woK0_woSm_data_cut = new TCanvas("cMMnmiss_IMnpip_dE_woK0_woSm_data_cut","cMMnmiss_IMnpip_dE_woK0_woSm_data_cut",800,800);
  MMnmiss_IMnpip_dE_woK0_woSm_data->Draw("colz");
  MMnmiss_IMnpip_dE_woK0_woSm_vici_data->SetMaximum(MMnmiss_IMnpip_dE_woK0_woSm_data->GetMaximum());
  MMnmiss_IMnpip_dE_woK0_woSm_vici_data->Draw("contzsame");
  
  TCanvas *cMMnmiss_woK0_woSm_data_vici = new TCanvas("cMMnmiss_woK0_woSm_data_vici","cMMnmiss_woK0_woSm_data_vici",800,800);
  TH1D* MMnmiss_woK0_woSm_vici_data = (TH1D*)MMnmiss_IMnpip_dE_woK0_woSm_vici_data->ProjectionY("MMnmiss_dE_woK0_woSm_vici_data");
  TH1D* MMnmiss_woK0_woSm_vici_mix  = (TH1D*)MMnmiss_IMnpip_dE_woK0_woSm_vici_mix->ProjectionY("MMnmiss_dE_woK0_woSm_vici_mix");
  MMnmiss_woK0_woSm_vici_data->Draw("HE");
  MMnmiss_woK0_woSm_vici_mix->SetLineColor(2);
  MMnmiss_woK0_woSm_vici_mix->Scale(mixScale);
  MMnmiss_woK0_woSm_vici_mix->Draw("HEsame");
  
  TCanvas *cMMnmiss_woK0_woSm_data_vici_sub = new TCanvas("cMMnmiss_woK0_woSm_data_vici_sub","cMMnmiss_woK0_woSm_data_vici_sub");
  TH1D* MMnmiss_woK0_woSm_vici_data_sub = (TH1D*)MMnmiss_woK0_woSm_vici_data->Clone("MMnmiss_woK0_woSm_vici_data_sub");
  MMnmiss_woK0_woSm_vici_data_sub->Add(MMnmiss_woK0_woSm_vici_mix,-1.0);
  MMnmiss_woK0_woSm_vici_data_sub->Draw("HE");
  
  //
  TCanvas *cIMnpip_woK0_woSm_data_vici = new TCanvas("cIMnpip_woK0_woSm_data_vici","cIMnpip_woK0_woSm_data_vici",800,800);
  TH1D* IMnpip_woK0_woSm_vici_data = (TH1D*)MMnmiss_IMnpip_dE_woK0_woSm_vici_data->ProjectionX("IMnpip_dE_woK0_woSm_vici_data");
  TH1D* IMnpip_woK0_woSm_vici_mix  = (TH1D*)MMnmiss_IMnpip_dE_woK0_woSm_vici_mix->ProjectionX("IMnpip_dE_woK0_woSm_vici_mix");
  IMnpip_woK0_woSm_vici_data->Draw("HE");
  IMnpip_woK0_woSm_vici_mix->SetLineColor(2);
  IMnpip_woK0_woSm_vici_mix->Scale(mixScale);
  IMnpip_woK0_woSm_vici_mix->Draw("HEsame");
  
  TCanvas *cIMnpip_woK0_woSm_data_vici_sub = new TCanvas("cIMnpip_woK0_woSm_data_vici_sub","cIMnpip_woK0_woSm_data_vici_sub");
  TH1D* IMnpip_woK0_woSm_data_vici_sub = (TH1D*)IMnpip_woK0_woSm_vici_data->Clone("IMnpip_woK0_woSm_data_vici_sub");
  IMnpip_woK0_woSm_data_vici_sub->Add(IMnpip_woK0_woSm_vici_mix,-1.0);
  IMnpip_woK0_woSm_data_vici_sub->Draw("HE");




  //MMnmiss vs IMnpim
  TCanvas *cMMnmiss_IMnpim_dE_woK0_woSp_data_cut = new TCanvas("cMMnmiss_IMnpim_dE_woK0_woSp_data_cut","cMMnmiss_IMnpim_dE_woK0_woSp_data_cut",800,800);
  MMnmiss_IMnpim_dE_woK0_woSp_data->Draw("colz");
  MMnmiss_IMnpim_dE_woK0_woSp_vici_data->SetMaximum(MMnmiss_IMnpim_dE_woK0_woSp_data->GetMaximum());
  MMnmiss_IMnpim_dE_woK0_woSp_vici_data->Draw("contzsame");


  TCanvas *cMMnmiss_IMnpip_dE_woK0_woSm_mix_cut = new TCanvas("cMMnmiss_IMnpip_dE_woK0_woSm_mix_cut","cMMnmiss_IMnpip_dE_woK0_woSm_mix_cut",800,800);
  MMnmiss_IMnpip_dE_woK0_woSm_mix->Draw("colz");
  MMnmiss_IMnpip_dE_woK0_woSm_vici_mix->SetMaximum(MMnmiss_IMnpip_dE_woK0_woSm_mix->GetMaximum());
  MMnmiss_IMnpip_dE_woK0_woSm_vici_mix->Draw("contzsame");

  TCanvas *cMMnmiss_IMnpim_dE_woK0_woSp_mix_cut = new TCanvas("cMMnmiss_IMnpim_dE_woK0_woSp_mix_cut","cMMnmiss_IMnpim_dE_woK0_woSp_mix_cut",800,800);
  MMnmiss_IMnpim_dE_woK0_woSp_mix->Draw("colz");
  MMnmiss_IMnpim_dE_woK0_woSp_vici_mix->SetMaximum(MMnmiss_IMnpim_dE_woK0_woSp_mix->GetMaximum());
  MMnmiss_IMnpim_dE_woK0_woSp_vici_mix->Draw("contzsame");






  return;

  //draw BG at first
  TH2D* IMnpim_IMnpip_woSid_won_data = (TH2D*)fr->Get("IMnpim_IMnpip_dE_woSid_won");IMnpim_IMnpip_woSid_won_data->Print();
  TH2D* IMnpim_IMnpip_woSid_won_mix= (TH2D*)fmix->Get("IMnpim_IMnpip_dE_woSid_won");IMnpim_IMnpip_woSid_won_mix->Print();
   
  const double ndata = IMnpim_IMnpip_woSid_won_data->Integral();
  const double nmix = IMnpim_IMnpip_woSid_won_mix->Integral();
  IMnpim_IMnpip_woSid_won_mix->Scale(ndata/nmix);

  TCanvas *cIMnpim_IMnpip_woSid_won = new TCanvas("cIMnpim_IMnpip_woSid_won","cIMnpim_IMnpip_woSid_won",1000,1000);
  cIMnpim_IMnpip_woSid_won->Divide(2,2,0,0);
  cIMnpim_IMnpip_woSid_won->cd(3);
  IMnpim_IMnpip_woSid_won_data->SetTitle("");
  IMnpim_IMnpip_woSid_won_data->RebinX(4);
  IMnpim_IMnpip_woSid_won_data->RebinY(4);
  IMnpim_IMnpip_woSid_won_data->GetXaxis()->SetRangeUser(1,1.7);
  IMnpim_IMnpip_woSid_won_data->GetYaxis()->SetRangeUser(1,1.7);
  IMnpim_IMnpip_woSid_won_mix->RebinX(4);
  IMnpim_IMnpip_woSid_won_mix->RebinY(4);
  IMnpim_IMnpip_woSid_won_mix->GetXaxis()->SetRangeUser(1,1.7);
  IMnpim_IMnpip_woSid_won_mix->GetYaxis()->SetRangeUser(1,1.7);
  IMnpim_IMnpip_woSid_won_data->SetContour(8);
  IMnpim_IMnpip_woSid_won_mix->SetContour(8);
  IMnpim_IMnpip_woSid_won_data->SetMinimum(1);
  IMnpim_IMnpip_woSid_won_data->SetLineWidth(2);
  IMnpim_IMnpip_woSid_won_data->Draw("cont3");
  IMnpim_IMnpip_woSid_won_mix->SetLineColor(6);
  IMnpim_IMnpip_woSid_won_mix->SetLineWidth(2);
  IMnpim_IMnpip_woSid_won_mix->SetMinimum(1);
  IMnpim_IMnpip_woSid_won_mix->SetMaximum(IMnpim_IMnpip_woSid_won_data->GetMaximum());
  IMnpim_IMnpip_woSid_won_mix->Draw("cont3same");
  //TPaletteAxis *palette2 = (TPaletteAxis*)IMnpim_IMnpip_woSid_won_data->GetListOfFunctions()->FindObject("palette");
  // the following lines moe the paletter. Choose the values you need for the position.
  //palette2->SetX1NDC(0.78);
  //palette2->SetX2NDC(0.83);
  //palette2->SetY1NDC(0.1);
  //palette2->SetY2NDC(0.90);
   
  
  cIMnpim_IMnpip_woSid_won->cd(1);
  TH1D* IMnpip_woSid_won_data = (TH1D*)IMnpim_IMnpip_woSid_won_data->ProjectionX("IMnpip_woSid_won_data");
  TH1D* IMnpip_woSid_won_mix = (TH1D*)IMnpim_IMnpip_woSid_won_mix->ProjectionX("IMnpip_woSid_won_mix");
  IMnpip_woSid_won_data->GetXaxis()->SetLabelSize(0);
  IMnpip_woSid_won_data->SetMarkerStyle(20);
  IMnpip_woSid_won_data->SetMarkerColor(1);
  IMnpip_woSid_won_data->Draw("E");
  IMnpip_woSid_won_mix->SetLineColor(6);
  IMnpip_woSid_won_mix->SetMarkerStyle(20);
  IMnpip_woSid_won_mix->SetMarkerColor(6);
  IMnpip_woSid_won_mix->Draw("Esame");
  
  cIMnpim_IMnpip_woSid_won->cd(4);
  TH1D* IMnpim_woSid_won_data = (TH1D*)IMnpim_IMnpip_woSid_won_data->ProjectionY("IMnpim_woSid_won_data");
  TH1D* IMnpim_woSid_won_mix = (TH1D*)IMnpim_IMnpip_woSid_won_mix->ProjectionY("IMnpim_woSid_won_mix");
  IMnpim_woSid_won_data->SetTitle("");
  IMnpim_woSid_won_data->GetXaxis()->SetTitle("");
  
  TGraphErrors *gr_IMnpim_woSid_won_data = new TGraphErrors();
  TGraphErrors *gr_IMnpim_woSid_won_mix = new TGraphErrors();
  
  for(int ibin=0;ibin<IMnpim_woSid_won_data->GetNbinsX();ibin++){
    double cont = IMnpim_woSid_won_data->GetBinContent(ibin);
    double err = IMnpim_woSid_won_data->GetBinError(ibin);
    double bincenter = IMnpim_woSid_won_data->GetBinCenter(ibin);
    gr_IMnpim_woSid_won_data->SetPoint(ibin,cont, bincenter);
    gr_IMnpim_woSid_won_data->SetPointError(ibin,err,0);
  }
  for(int ibin=0;ibin<IMnpim_woSid_won_mix->GetNbinsX();ibin++){
    double cont = IMnpim_woSid_won_mix->GetBinContent(ibin);
    double err = IMnpim_woSid_won_mix->GetBinError(ibin);
    double bincenter = IMnpim_woSid_won_mix->GetBinCenter(ibin);
    gr_IMnpim_woSid_won_mix->SetPoint(ibin,cont, bincenter);
    gr_IMnpim_woSid_won_mix->SetPointError(ibin,err,0);
  }
  gr_IMnpim_woSid_won_data->GetYaxis()->SetRangeUser(1,1.7);
  gr_IMnpim_woSid_won_data->SetMarkerStyle(20);
  gr_IMnpim_woSid_won_data->GetYaxis()->SetLabelSize(0);
  gr_IMnpim_woSid_won_data->Draw("AP");
  gr_IMnpim_woSid_won_mix->SetMarkerStyle(20);
  gr_IMnpim_woSid_won_mix->SetMarkerColor(6);
  gr_IMnpim_woSid_won_mix->Draw("P");
   
  TCanvas *cIMnpim_IMnpip_woSid_won_2 = new TCanvas("cIMnpim_IMnpip_woSid_won_2","cIMnpim_IMnpip_woSid_won_2",1000,1000);
  cIMnpim_IMnpip_woSid_won_2->Divide(2,2);
  cIMnpim_IMnpip_woSid_won_2->cd(1);
  TH2D* IMnpim_IMnpip_woSid_won_data_2 = (TH2D*)IMnpim_IMnpip_woSid_won_data->Clone("IMnpim_IMnpip_woSid_won_data_2");
  IMnpim_IMnpip_woSid_won_data_2->SetContour(10);
  IMnpim_IMnpip_woSid_won_data_2->Draw("cont1");
  gPad->SetRightMargin(0);
  cIMnpim_IMnpip_woSid_won_2->cd(2);
  TH2D* IMnpim_IMnpip_woSid_won_mix_2 = (TH2D*)IMnpim_IMnpip_woSid_won_mix->Clone("IMnpim_IMnpip_woSid_won_mix_2");
  IMnpim_IMnpip_woSid_won_mix_2->SetTitle("");
  IMnpim_IMnpip_woSid_won_mix_2->SetContour(10);
  IMnpim_IMnpip_woSid_won_mix_2->SetMaximum(IMnpim_IMnpip_woSid_won_data->GetMaximum());
  IMnpim_IMnpip_woSid_won_mix_2->Draw("cont1z");
  gPad->SetLeftMargin(0);
  cIMnpim_IMnpip_woSid_won_2->cd(3);
  IMnpim_IMnpip_woSid_won_data_2->Draw("col");
  gPad->SetRightMargin(0);
  cIMnpim_IMnpip_woSid_won_2->cd(4);
  IMnpim_IMnpip_woSid_won_mix_2->Draw("colz");
  gPad->SetLeftMargin(0);


  TCanvas *cIMnpim_IMnpip_woSid_won_sub = new TCanvas("cIMnpim_IMnpip_woSid_won_sub","cIMnpim_IMnpip_woSid_won_sub",800,800);
  TH2D* IMnpim_IMnpip_woSid_won_sub = (TH2D*)IMnpim_IMnpip_woSid_won_data_2->Clone("IMnpim_IMnpip_woSid_won_sub");
  IMnpim_IMnpip_woSid_won_sub->Add(IMnpim_IMnpip_woSid_won_mix,-1.0);
  IMnpim_IMnpip_woSid_won_sub->Divide(IMnpim_IMnpip_woSid_won_mix);
  IMnpim_IMnpip_woSid_won_sub->SetTitle("(reco. data - MC)/MC");
  //IMnpim_IMnpip_woSid_won_sub->RebinX(2);
  //IMnpim_IMnpip_woSid_won_sub->RebinY(2);
  //IMnpim_IMnpip_woSid_won_sub->Scale(1./4.);
  IMnpim_IMnpip_woSid_won_sub->GetXaxis()->SetRangeUser(1.0,1.7);
  IMnpim_IMnpip_woSid_won_sub->GetYaxis()->SetRangeUser(1.0,1.7);
  IMnpim_IMnpip_woSid_won_sub->GetZaxis()->SetRangeUser(-0.5,0.5);
  TBox *box_sigmap = new TBox(anacuts::Sigmap_MIN_wide,1.0,anacuts::Sigmap_MAX_wide,1.7);
  TBox *box_sigmam = new TBox(1.0,anacuts::Sigmam_MIN_wide,1.7,anacuts::Sigmam_MAX_wide);
  box_sigmap->SetFillColor(1);
  box_sigmam->SetFillColor(1);
  IMnpim_IMnpip_woSid_won_sub->Draw("colz");
  box_sigmap->Draw();
  box_sigmam->Draw();
  
  //
  //BG2 MMnmiss vs IMpipim
  //
  TH2D* MMnmiss_IMpippim_woSid_won_data = (TH2D*)fr->Get("MMnmiss_IMpippim_dE_woSid_won");MMnmiss_IMpippim_woSid_won_data->Print();
  TH2D* MMnmiss_IMpippim_woSid_won_mix = (TH2D*)fmix->Get("MMnmiss_IMpippim_dE_woSid_won");MMnmiss_IMpippim_woSid_won_mix->Print();
  MMnmiss_IMpippim_woSid_won_mix->Scale(ndata/nmix);

  TCanvas *cMMnmiss_IMpippim_woSid_won = new TCanvas("cMMnmiss_IMpippim_woSid_won","cMMnmiss_IMpippim_woSid_won",1000,1000);
  cMMnmiss_IMpippim_woSid_won ->Divide(2,2,0,0);
  
 
  cMMnmiss_IMpippim_woSid_won->cd(3);
  MMnmiss_IMpippim_woSid_won_data->SetTitle("");
  MMnmiss_IMpippim_woSid_won_data->RebinX(4);
  MMnmiss_IMpippim_woSid_won_data->RebinY(3);
  MMnmiss_IMpippim_woSid_won_data->GetXaxis()->SetRangeUser(0.2,1.0);
  MMnmiss_IMpippim_woSid_won_mix->RebinX(4);
  MMnmiss_IMpippim_woSid_won_mix->RebinY(3);
  MMnmiss_IMpippim_woSid_won_mix->GetXaxis()->SetRangeUser(0.2,1.0);
  MMnmiss_IMpippim_woSid_won_data->SetContour(10);
  MMnmiss_IMpippim_woSid_won_mix->SetContour(10);

  MMnmiss_IMpippim_woSid_won_data->SetLineWidth(2);
  MMnmiss_IMpippim_woSid_won_data->SetMinimum(1);
  MMnmiss_IMpippim_woSid_won_data->Draw("cont3");
  MMnmiss_IMpippim_woSid_won_mix->SetLineColor(6);
  MMnmiss_IMpippim_woSid_won_mix->SetLineWidth(2);
  MMnmiss_IMpippim_woSid_won_mix->SetMinimum(1);
  MMnmiss_IMpippim_woSid_won_mix->SetMaximum(MMnmiss_IMpippim_woSid_won_data->GetMaximum());
  MMnmiss_IMpippim_woSid_won_mix->Draw("cont3same");
   
  
  cMMnmiss_IMpippim_woSid_won->cd(1);
  TH1D* IMpippim_woSid_won_data = (TH1D*)MMnmiss_IMpippim_woSid_won_data->ProjectionX("IMpippim_woSid_won_data");
  TH1D* IMpippim_woSid_won_mix = (TH1D*)MMnmiss_IMpippim_woSid_won_mix->ProjectionX("IMpippim_woSid_won_mix");
  IMpippim_woSid_won_data->SetTitle("");
  IMpippim_woSid_won_data->GetXaxis()->SetLabelSize(0);
  IMpippim_woSid_won_data->SetMarkerStyle(20);
  IMpippim_woSid_won_data->SetMarkerColor(1);
  IMpippim_woSid_won_data->Draw("E");
  IMpippim_woSid_won_mix->SetLineColor(6);
  IMpippim_woSid_won_mix->SetMarkerStyle(20);
  IMpippim_woSid_won_mix->SetMarkerColor(6);
  IMpippim_woSid_won_mix->Draw("Esame");
  
  cMMnmiss_IMpippim_woSid_won->cd(4);
  TH1D* MMnmiss_woSid_won_data = (TH1D*)MMnmiss_IMpippim_woSid_won_data->ProjectionY("MMnmiss_woSid_won_data");
  TH1D* MMnmiss_woSid_won_mix = (TH1D*)MMnmiss_IMpippim_woSid_won_mix->ProjectionY("MMnmiss_woSid_won_mix");
  MMnmiss_woSid_won_data->SetTitle("");
  MMnmiss_woSid_won_data->GetXaxis()->SetTitle("");
  
  TGraphErrors *gr_MMnmiss_woSid_won_data = new TGraphErrors();
  TGraphErrors *gr_MMnmiss_woSid_won_mix = new TGraphErrors();
  
  for(int ibin=0;ibin<MMnmiss_woSid_won_data->GetNbinsX();ibin++){
    double cont = MMnmiss_woSid_won_data->GetBinContent(ibin);
    double err = MMnmiss_woSid_won_data->GetBinError(ibin);
    double bincenter = MMnmiss_woSid_won_data->GetBinCenter(ibin);
    gr_MMnmiss_woSid_won_data->SetPoint(ibin,cont, bincenter);
    gr_MMnmiss_woSid_won_data->SetPointError(ibin,err,0);
  }
  for(int ibin=0;ibin<MMnmiss_woSid_won_mix->GetNbinsX();ibin++){
    double cont = MMnmiss_woSid_won_mix->GetBinContent(ibin);
    double err = MMnmiss_woSid_won_mix->GetBinError(ibin);
    double bincenter = MMnmiss_woSid_won_mix->GetBinCenter(ibin);
    gr_MMnmiss_woSid_won_mix->SetPoint(ibin,cont, bincenter);
    gr_MMnmiss_woSid_won_mix->SetPointError(ibin,err,0);
  }
  gr_MMnmiss_woSid_won_data->GetYaxis()->SetRangeUser(0,1.5);
  gr_MMnmiss_woSid_won_data->SetMarkerStyle(20);
  gr_MMnmiss_woSid_won_data->GetYaxis()->SetLabelSize(0);
  gr_MMnmiss_woSid_won_data->Draw("AP");
  gr_MMnmiss_woSid_won_mix->SetMarkerStyle(20);
  gr_MMnmiss_woSid_won_mix->SetMarkerColor(6);
  gr_MMnmiss_woSid_won_mix->Draw("P");

  TCanvas *cMMnmiss_IMpippim_woSid_won_2 = new TCanvas("cMMnmiss_IMpippim_woSid_won_2","cMMnmiss_IMpippim_woSid_won_2",1000,1000);
  cMMnmiss_IMpippim_woSid_won_2->Divide(2,2);
  cMMnmiss_IMpippim_woSid_won_2->cd(1);
  TH2D* MMnmiss_IMpippim_woSid_won_data_2 = (TH2D*)MMnmiss_IMpippim_woSid_won_data->Clone("MMnmiss_IMpippim_woSid_won_data_2");
  TH2D* MMnmiss_IMpippim_woSid_won_mix_2 = (TH2D*)MMnmiss_IMpippim_woSid_won_mix->Clone("MMnmiss_IMpippim_woSid_won_mix_2");
  MMnmiss_IMpippim_woSid_won_data_2->SetTitle("reco data");
  MMnmiss_IMpippim_woSid_won_data_2->SetContour(10); 
  MMnmiss_IMpippim_woSid_won_data_2->SetMinimum(1);
  MMnmiss_IMpippim_woSid_won_data_2->Draw("cont1z");
  gPad->SetRightMargin(0);
  cMMnmiss_IMpippim_woSid_won_2->cd(2);
  MMnmiss_IMpippim_woSid_won_mix_2->SetTitle("MC");
  MMnmiss_IMpippim_woSid_won_mix_2->SetContour(10); 
  MMnmiss_IMpippim_woSid_won_mix_2->SetMaximum(MMnmiss_IMpippim_woSid_won_data_2->GetMaximum());
  MMnmiss_IMpippim_woSid_won_mix_2->SetMinimum(1);
  MMnmiss_IMpippim_woSid_won_mix_2->Draw("cont1z");
  gPad->SetLeftMargin(0);
  cMMnmiss_IMpippim_woSid_won_2->cd(3);
  MMnmiss_IMpippim_woSid_won_data_2->Draw("col");
  gPad->SetRightMargin(0);
  cMMnmiss_IMpippim_woSid_won_2->cd(4);
  MMnmiss_IMpippim_woSid_won_mix_2->Draw("colz");
  gPad->SetLeftMargin(0);
  

  TCanvas *cMMnmiss_IMpippim_woSid_won_sub = new TCanvas("cMMnmiss_IMpippim_woSid_won_sub","cMMnmiss_IMpippim_woSid_won_sub",800,800);
  TH2D* MMnmiss_IMpippim_woSid_won_sub = (TH2D*)MMnmiss_IMpippim_woSid_won_data->Clone("MMnmiss_IMpippim_woSid_won_sub");
  MMnmiss_IMpippim_woSid_won_sub->Add(MMnmiss_IMpippim_woSid_won_mix,-1.0);
  MMnmiss_IMpippim_woSid_won_sub->Divide(MMnmiss_IMpippim_woSid_won_mix);
  MMnmiss_IMpippim_woSid_won_sub->SetTitle("(data-MC)/MC");
  MMnmiss_IMpippim_woSid_won_sub->RebinX(4);
  MMnmiss_IMpippim_woSid_won_sub->RebinY(2);
  MMnmiss_IMpippim_woSid_won_sub->Scale(1./8.);
  MMnmiss_IMpippim_woSid_won_sub->GetZaxis()->SetRangeUser(-0.3,0.3);
  MMnmiss_IMpippim_woSid_won_sub->Draw("colz");
  TBox *box_pipi = new TBox(anacuts::pipi_MIN,0,anacuts::pipi_MAX,1.5);
  TBox *box_MMnmiss = new TBox(0,anacuts::neutron_MIN_wide,1.0,anacuts::neutron_MAX_wide);
  box_pipi->SetFillColor(1);
  box_pipi->SetFillStyle(3002);
  box_MMnmiss->SetFillColor(1);
  //box_pipi->Draw();
  box_MMnmiss->Draw();
 
  return;
  //
  //BG3 q vs nmom
  //
  TH2D* q_nmom_woK0_woSid_won_data = (TH2D*)f->Get("q_nmom_woK0_woSid_won_data");q_nmom_woK0_woSid_won_data->Print();
  TH2D* q_nmom_wK0_woSid_won_data = (TH2D*)f->Get("q_nmom_wK0_woSid_won_data");q_nmom_wK0_woSid_won_data->Print();

  TH2D* q_nmom_woK0_woSid_won_mix = (TH2D*)f->Get("q_nmom_woK0_woSid_won_mix");q_nmom_woK0_woSid_won_mix->Print();

  TH2D* q_nmom_wK0_woSid_won_mix = (TH2D*)f->Get("q_nmom_wK0_woSid_won_mix");q_nmom_wK0_woSid_won_mix->Print();
  TH2D* q_nmom_wK0_woSid_won_mixgeta = (TH2D*)f->Get("q_nmom_wK0_woSid_won_mixgeta");q_nmom_wK0_woSid_won_mixgeta->Print();

  TH2D* q_nmom_woSid_won_data = (TH2D*)q_nmom_woK0_woSid_won_data->Clone("q_nmom_woSid_won_data");
  q_nmom_woSid_won_data->Add(q_nmom_wK0_woSid_won_data);
  TH2D* q_nmom_woSid_won_mix = (TH2D*)q_nmom_woK0_woSid_won_mix->Clone("q_nmom_woSid_won_mix");
  q_nmom_woSid_won_mix->Add(q_nmom_wK0_woSid_won_mix);
   
  //divide pipin"X" BG and K0n"X" BG
  TH2D* q_nmom_woSid_won_mixpipinX = (TH2D*)q_nmom_woK0_woSid_won_mix->Clone("q_nmom_woSid_won_mixpipinX");
  q_nmom_woSid_won_mixpipinX->Add(q_nmom_wK0_woSid_won_mixgeta,1.0);
  TH2D* q_nmom_woSid_won_mixK0nX = (TH2D*)q_nmom_wK0_woSid_won_mix->Clone("q_nmom_wSid_won_mixK0nX");
  q_nmom_woSid_won_mixK0nX->Add(q_nmom_wK0_woSid_won_mixgeta,-1.0);
  
  TCanvas *cq_nmom_woSid_won = new TCanvas("cq_nmom_woSid_won","cq_nmom_woSid_won",1000,1000);
  cq_nmom_woSid_won ->Divide(2,2,0,0);
  
 
  cq_nmom_woSid_won->cd(3);
  q_nmom_woSid_won_data->SetTitle("");
  q_nmom_woSid_won_data->RebinX(4);
  q_nmom_woSid_won_data->RebinY(4);
  q_nmom_woSid_won_data->GetXaxis()->SetRangeUser(0.2,1.0);
  //q_nmom_woSid_won_data->GetYaxis()->SetRangeUser(1,1.7);
  q_nmom_woSid_won_mix->RebinX(4);
  q_nmom_woSid_won_mix->RebinY(4);
  q_nmom_woSid_won_mixpipinX->RebinX(4);
  q_nmom_woSid_won_mixpipinX->RebinY(4);
  q_nmom_woSid_won_mixK0nX->RebinX(4);
  q_nmom_woSid_won_mixK0nX->RebinY(4);
  q_nmom_woSid_won_mix->GetXaxis()->SetRangeUser(0.2,1.0);
  q_nmom_woSid_won_mixpipinX->GetXaxis()->SetRangeUser(0.2,1.0);
  q_nmom_woSid_won_mixK0nX->GetXaxis()->SetRangeUser(0.2,1.0);
  q_nmom_woSid_won_data->SetContour(10);
  q_nmom_woSid_won_mix->SetContour(10);
  q_nmom_woSid_won_data->SetMinimum(1);
  q_nmom_woSid_won_data->SetLineWidth(2);
  q_nmom_woSid_won_data->Draw("cont3");
  q_nmom_woSid_won_mix->SetLineColor(6);
  q_nmom_woSid_won_mix->SetLineWidth(2);
  q_nmom_woSid_won_mix->SetMaximum(q_nmom_woSid_won_data->GetMaximum());
  q_nmom_woSid_won_mix->SetMinimum(1);
  q_nmom_woSid_won_mix->Draw("cont3same");
   
  
  cq_nmom_woSid_won->cd(1);
  TH1D* nmom_woSid_won_data = (TH1D*)q_nmom_woSid_won_data->ProjectionX("nmom_woSid_won_data");
  TH1D* nmom_woSid_won_mix = (TH1D*)q_nmom_woSid_won_mix->ProjectionX("nmom_woSid_won_mix");
  TH1D* nmom_woSid_won_mixpipinX = (TH1D*)q_nmom_woSid_won_mixpipinX->ProjectionX("nmom_woSid_won_mixpipinX");
  TH1D* nmom_woSid_won_mixK0nX = (TH1D*)q_nmom_woSid_won_mixK0nX->ProjectionX("nmom_woSid_won_mixK0nX");
  nmom_woSid_won_data->SetTitle("");
  nmom_woSid_won_data->GetXaxis()->SetLabelSize(0);
  nmom_woSid_won_data->SetMarkerStyle(20);
  nmom_woSid_won_data->SetMarkerColor(1);
  nmom_woSid_won_data->Draw("E");
  nmom_woSid_won_mix->SetLineColor(6);
  nmom_woSid_won_mix->SetMarkerStyle(20);
  nmom_woSid_won_mix->SetMarkerColor(6);
  nmom_woSid_won_mix->Draw("Esame");
  nmom_woSid_won_mixpipinX->SetLineColor(2);
  nmom_woSid_won_mixpipinX->SetMarkerStyle(20);
  nmom_woSid_won_mixpipinX->SetMarkerColor(2);
  if(showBG)nmom_woSid_won_mixpipinX->Draw("Esame");
  nmom_woSid_won_mixK0nX->SetLineColor(3);
  nmom_woSid_won_mixK0nX->SetMarkerStyle(20);
  nmom_woSid_won_mixK0nX->SetMarkerColor(3);
  if(showBG)nmom_woSid_won_mixK0nX->Draw("Esame");
  
  cq_nmom_woSid_won->cd(4);
  TH1D* q_woSid_won_data = (TH1D*)q_nmom_woSid_won_data->ProjectionY("q_woSid_won_data");
  TH1D* q_woSid_won_mix = (TH1D*)q_nmom_woSid_won_mix->ProjectionY("q_woSid_won_mix");
  TH1D* q_woSid_won_mixpipinX = (TH1D*)q_nmom_woSid_won_mixpipinX->ProjectionY("q_woSid_won_mixpipinX");
  TH1D* q_woSid_won_mixK0nX = (TH1D*)q_nmom_woSid_won_mixK0nX->ProjectionY("q_woSid_won_mixK0nX");
  q_woSid_won_data->SetTitle("");
  q_woSid_won_data->GetXaxis()->SetTitle("");
  
  TGraphErrors *gr_q_woSid_won_data = new TGraphErrors();
  TGraphErrors *gr_q_woSid_won_mix = new TGraphErrors();
  TGraphErrors *gr_q_woSid_won_mixpipinX = new TGraphErrors();
  TGraphErrors *gr_q_woSid_won_mixK0nX = new TGraphErrors();
  
  for(int ibin=0;ibin<q_woSid_won_data->GetNbinsX();ibin++){
    double cont = q_woSid_won_data->GetBinContent(ibin);
    double err = q_woSid_won_data->GetBinError(ibin);
    double bincenter = q_woSid_won_data->GetBinCenter(ibin);
    gr_q_woSid_won_data->SetPoint(ibin,cont, bincenter);
    gr_q_woSid_won_data->SetPointError(ibin,err,0);
  }
  for(int ibin=0;ibin<q_woSid_won_mix->GetNbinsX();ibin++){
    double cont = q_woSid_won_mix->GetBinContent(ibin);
    double err = q_woSid_won_mix->GetBinError(ibin);
    double bincenter = q_woSid_won_mix->GetBinCenter(ibin);
    gr_q_woSid_won_mix->SetPoint(ibin,cont, bincenter);
    gr_q_woSid_won_mix->SetPointError(ibin,err,0);
  }
  for(int ibin=0;ibin<q_woSid_won_mixpipinX->GetNbinsX();ibin++){
    double cont = q_woSid_won_mixpipinX->GetBinContent(ibin);
    double err = q_woSid_won_mixpipinX->GetBinError(ibin);
    double bincenter = q_woSid_won_mixpipinX->GetBinCenter(ibin);
    gr_q_woSid_won_mixpipinX->SetPoint(ibin,cont, bincenter);
    gr_q_woSid_won_mixpipinX->SetPointError(ibin,err,0);
  }
  for(int ibin=0;ibin<q_woSid_won_mixK0nX->GetNbinsX();ibin++){
    double cont = q_woSid_won_mixK0nX->GetBinContent(ibin);
    double err = q_woSid_won_mixK0nX->GetBinError(ibin);
    double bincenter = q_woSid_won_mixK0nX->GetBinCenter(ibin);
    gr_q_woSid_won_mixK0nX->SetPoint(ibin,cont, bincenter);
    gr_q_woSid_won_mixK0nX->SetPointError(ibin,err,0);
  }

  gr_q_woSid_won_data->GetYaxis()->SetRangeUser(0,1.5);
  gr_q_woSid_won_data->SetMarkerStyle(20);
  gr_q_woSid_won_data->GetYaxis()->SetLabelSize(0);
  gr_q_woSid_won_data->Draw("AP");
  gr_q_woSid_won_mix->SetMarkerStyle(20);
  gr_q_woSid_won_mix->SetMarkerColor(6);
  gr_q_woSid_won_mix->Draw("P");
  gr_q_woSid_won_mixpipinX->SetMarkerStyle(20);
  gr_q_woSid_won_mixpipinX->SetMarkerColor(2);
  if(showBG)gr_q_woSid_won_mixpipinX->Draw("P");
  gr_q_woSid_won_mixK0nX->SetMarkerStyle(20);
  gr_q_woSid_won_mixK0nX->SetMarkerColor(3);
  if(showBG)gr_q_woSid_won_mixK0nX->Draw("P");

  TCanvas *cq_nmom_woSid_won_2 = new TCanvas("cq_nmom_woSid_won_2","cq_nmom_woSid_won_2",1000,1000);
  cq_nmom_woSid_won_2->Divide(2,2);
  cq_nmom_woSid_won_2->cd(1);
  TH2D* q_nmom_woSid_won_data_2 = (TH2D*)q_nmom_woSid_won_data->Clone("q_nmom_woSid_won_data_2");
  TH2D* q_nmom_woSid_won_mix_2 = (TH2D*)q_nmom_woSid_won_mix->Clone("q_nmom_woSid_won_mix_2");
  q_nmom_woSid_won_data_2->SetTitle("real data");
  q_nmom_woSid_won_data_2->SetContour(10); 
  q_nmom_woSid_won_data_2->SetMinimum(1);
  q_nmom_woSid_won_data_2->Draw("cont1");
  gPad->SetRightMargin(0);
  cq_nmom_woSid_won_2->cd(2);
  q_nmom_woSid_won_mix_2->SetTitle("MC");
  q_nmom_woSid_won_mix_2->SetContour(10); 
  q_nmom_woSid_won_mix_2->SetMaximum(q_nmom_woSid_won_data_2->GetMaximum());
  q_nmom_woSid_won_mix_2->SetMinimum(1);
  q_nmom_woSid_won_mix_2->Draw("cont1z");
  gPad->SetLeftMargin(0);
  cq_nmom_woSid_won_2->cd(3);
  q_nmom_woSid_won_data_2->Draw("colz");
  gPad->SetRightMargin(0);
  cq_nmom_woSid_won_2->cd(4);
  q_nmom_woSid_won_mix_2->Draw("colz");
  gPad->SetLeftMargin(0);
  
  TCanvas *cq_nmom_woSid_won_sub = new TCanvas("cq_nmom_woSid_won_sub","cq_nmom_woSid_won_sub",800,800);
  TH2D* q_nmom_woSid_won_sub = (TH2D*)q_nmom_woSid_won_data->Clone("q_nmom_woSid_won_data");
  q_nmom_woSid_won_sub->Add(q_nmom_woSid_won_mix,-1.0);
  q_nmom_woSid_won_sub->Divide(q_nmom_woSid_won_mix);
  q_nmom_woSid_won_sub->SetTitle("(data-MC)/MC");
  q_nmom_woSid_won_sub->RebinX(2);
  q_nmom_woSid_won_sub->RebinY(2);
  q_nmom_woSid_won_sub->Scale(0.25);
  q_nmom_woSid_won_sub->GetZaxis()->SetRangeUser(-0.3,0.3);
  q_nmom_woSid_won_sub->Draw("colz");
  
  //
  //BG check q vs IMnpipi
  //
  TH2D* q_IMnpipi_woK0_woSid_won_data = (TH2D*)f->Get("q_IMnpipi_woK0_woSid_won_data");q_IMnpipi_woK0_woSid_won_data->Print();
  TH2D* q_IMnpipi_wK0_woSid_won_data = (TH2D*)f->Get("q_IMnpipi_wK0_woSid_won_data");q_IMnpipi_wK0_woSid_won_data->Print();

  TH2D* q_IMnpipi_woK0_woSid_won_mix = (TH2D*)f->Get("q_IMnpipi_woK0_woSid_won_mix");q_IMnpipi_woK0_woSid_won_mix->Print();

  TH2D* q_IMnpipi_wK0_woSid_won_mix = (TH2D*)f->Get("q_IMnpipi_wK0_woSid_won_mix");q_IMnpipi_wK0_woSid_won_mix->Print();

  TH2D* q_IMnpipi_woSid_won_data = (TH2D*)q_IMnpipi_woK0_woSid_won_data->Clone("q_IMnpipi_woSid_won_data");
  q_IMnpipi_woSid_won_data->Add(q_IMnpipi_wK0_woSid_won_data);
  TH2D* q_IMnpipi_woSid_won_mix = (TH2D*)q_IMnpipi_woK0_woSid_won_mix->Clone("q_IMnpipi_woSid_won_mix");
  q_IMnpipi_woSid_won_mix->Add(q_IMnpipi_wK0_woSid_won_mix);
  //q_IMnpipi_woSid_won_mix->Add(q_IMnpipi_wK0_woSid_won_mixgeta);
   
  
  TCanvas *cq_IMnpipi_woSid_won = new TCanvas("cq_IMnpipi_woSid_won","cq_IMnpipi_woSid_won",1000,1000);
  cq_IMnpipi_woSid_won ->Divide(2,2,0,0);
  
 
  cq_IMnpipi_woSid_won->cd(3);
  q_IMnpipi_woSid_won_data->SetTitle("");
  q_IMnpipi_woSid_won_data->RebinX(2);
  std::cout << "bin width x q_IMnpipi_woSid_won_data " << q_IMnpipi_woSid_won_data->GetXaxis()->GetBinWidth(1) << std::endl;
  q_IMnpipi_woSid_won_data->RebinY(4);
  std::cout << "bin width y q_IMnpipi_woSid_won_data " << q_IMnpipi_woSid_won_data->GetYaxis()->GetBinWidth(1) << std::endl;
  q_IMnpipi_woSid_won_data->GetXaxis()->SetRangeUser(0.2,1.0);
  //q_IMnpipi_woSid_won_data->GetYaxis()->SetRangeUser(1,1.7);
  q_IMnpipi_woSid_won_mix->RebinX(2);
  q_IMnpipi_woSid_won_mix->RebinY(4);
  q_IMnpipi_woSid_won_mix->GetXaxis()->SetRangeUser(0.2,1.0);
  //q_IMnpipi_woSid_won_mix->GetYaxis()->SetRangeUser(1,1.7);
  q_IMnpipi_woSid_won_data->SetContour(10);
  q_IMnpipi_woSid_won_mix->SetContour(10);
  q_IMnpipi_woSid_won_data->SetLineWidth(2);
  q_IMnpipi_woSid_won_data->Draw("cont3");
  q_IMnpipi_woSid_won_mix->SetLineColor(6);
  q_IMnpipi_woSid_won_mix->SetLineWidth(2);
  q_IMnpipi_woSid_won_mix->Draw("cont3same");
   
  
  cq_IMnpipi_woSid_won->cd(1);
  TH1D* IMnpipi_woSid_won_data = (TH1D*)q_IMnpipi_woSid_won_data->ProjectionX("IMnpipi_woSid_won_data");
  TH1D* IMnpipi_woSid_won_mix = (TH1D*)q_IMnpipi_woSid_won_mix->ProjectionX("IMnpipi_woSid_won_mix");
  IMnpipi_woSid_won_data->SetTitle("");
  IMnpipi_woSid_won_data->GetXaxis()->SetLabelSize(0);
  IMnpipi_woSid_won_data->SetMarkerStyle(20);
  IMnpipi_woSid_won_data->SetMarkerColor(1);
  IMnpipi_woSid_won_data->Draw("E");
  IMnpipi_woSid_won_mix->SetLineColor(6);
  IMnpipi_woSid_won_mix->SetMarkerStyle(20);
  IMnpipi_woSid_won_mix->SetMarkerColor(6);
  IMnpipi_woSid_won_mix->Draw("Esame");
  
  
  cq_IMnpipi_woSid_won->cd(4);
  TH1D* q_woSid_won_data_2 = (TH1D*)q_IMnpipi_woSid_won_data->ProjectionY("q_woSid_won_data");
  TH1D* q_woSid_won_mix_2 = (TH1D*)q_IMnpipi_woSid_won_mix->ProjectionY("q_woSid_won_mix");
  q_woSid_won_data_2->SetTitle("");
  q_woSid_won_data_2->GetXaxis()->SetTitle("");
  
  TGraphErrors *gr_q_woSid_won_data_2 = new TGraphErrors();
  TGraphErrors *gr_q_woSid_won_mix_2 = new TGraphErrors();
  
  for(int ibin=0;ibin<q_woSid_won_data_2->GetNbinsX();ibin++){
    double cont = q_woSid_won_data_2->GetBinContent(ibin);
    double err = q_woSid_won_data_2->GetBinError(ibin);
    double bincenter = q_woSid_won_data_2->GetBinCenter(ibin);
    gr_q_woSid_won_data_2->SetPoint(ibin,cont, bincenter);
    gr_q_woSid_won_data_2->SetPointError(ibin,err,0);
  }
  for(int ibin=0;ibin<q_woSid_won_mix_2->GetNbinsX();ibin++){
    double cont = q_woSid_won_mix_2->GetBinContent(ibin);
    double err = q_woSid_won_mix_2->GetBinError(ibin);
    double bincenter = q_woSid_won_mix_2->GetBinCenter(ibin);
    gr_q_woSid_won_mix_2->SetPoint(ibin,cont, bincenter);
    gr_q_woSid_won_mix_2->SetPointError(ibin,err,0);
  }
  gr_q_woSid_won_data_2->GetYaxis()->SetRangeUser(0,1.5);
  gr_q_woSid_won_data_2->SetMarkerStyle(20);
  gr_q_woSid_won_data_2->GetYaxis()->SetLabelSize(0);
  gr_q_woSid_won_data_2->Draw("AP");
  gr_q_woSid_won_mix_2->SetMarkerStyle(20);
  gr_q_woSid_won_mix_2->SetMarkerColor(6);
  gr_q_woSid_won_mix_2->Draw("P");

  TCanvas *cq_IMnpipi_woSid_won_2 = new TCanvas("cq_IMnpipi_woSid_won_2","cq_IMnpipi_woSid_won_2",1000,1000);
  cq_IMnpipi_woSid_won_2->Divide(2,2);
  cq_IMnpipi_woSid_won_2->cd(1);
  TH2D* q_IMnpipi_woSid_won_data_2 = (TH2D*)q_IMnpipi_woSid_won_data->Clone("q_IMnpipi_woSid_won_data_2");
  TH2D* q_IMnpipi_woSid_won_mix_2 = (TH2D*)q_IMnpipi_woSid_won_mix->Clone("q_IMnpipi_woSid_won_mix_2");
  q_IMnpipi_woSid_won_data_2->SetTitle("real data");
  q_IMnpipi_woSid_won_data_2->SetContour(10);
  q_IMnpipi_woSid_won_data_2->SetMinimum(1);
  q_IMnpipi_woSid_won_data_2->Draw("cont1");
  gPad->SetRightMargin(0);
  cq_IMnpipi_woSid_won_2->cd(2);
  q_IMnpipi_woSid_won_mix_2->SetTitle("MC");
  q_IMnpipi_woSid_won_mix_2->SetContour(10); 
  q_IMnpipi_woSid_won_mix_2->SetMinimum(1);
  q_IMnpipi_woSid_won_mix_2->SetMaximum(q_IMnpipi_woSid_won_data_2->GetMaximum());
  q_IMnpipi_woSid_won_mix_2->Draw("cont1z");
  gPad->SetLeftMargin(0);
  cq_IMnpipi_woSid_won_2->cd(3);
  q_IMnpipi_woSid_won_data_2->Draw("col");
  gPad->SetRightMargin(0);
  cq_IMnpipi_woSid_won_2->cd(4);
  q_IMnpipi_woSid_won_mix_2->Draw("colz");
  gPad->SetLeftMargin(0);


  
  TCanvas *cq_IMnpipi_woSid_won_sub = new TCanvas("cq_IMnpipi_woSid_won_sub","cq_IMnpipi_woSid_won_sub",800,800);
  TH2D* q_IMnpipi_woSid_won_sub = (TH2D*)q_IMnpipi_woSid_won_data->Clone("q_IMnpipi_woSid_won_data");
  q_IMnpipi_woSid_won_sub->Add(q_IMnpipi_woSid_won_mix,-1.0);
  q_IMnpipi_woSid_won_sub->Divide(q_IMnpipi_woSid_won_mix);
  q_IMnpipi_woSid_won_sub->SetTitle("(data-MC)/mc");
  q_IMnpipi_woSid_won_sub->RebinX(2);
  q_IMnpipi_woSid_won_sub->RebinY(2);
  q_IMnpipi_woSid_won_sub->Scale(0.25);
  q_IMnpipi_woSid_won_sub->GetZaxis()->SetRangeUser(-0.3,0.3);
  q_IMnpipi_woSid_won_sub->Draw("colz");

  
  //signal check IMnpim vs IMnpip
  const double Slow = 1.0;
  const double Shigh = 1.7;
  
  TCanvas *cIMnpim_IMnpip_n = new TCanvas("cIMnpim_IMnpip_n", "cIMnpim_IMnpip_n",1000,1000);     
  cIMnpim_IMnpip_n->Divide(2,2,0,0);                                                             
  TH2D* IMnpim_IMnpip_n_rdata = (TH2D*)f->Get("IMnpim_IMnpip_n_rdata");                          
  TH2D* IMnpim_IMnpip_n_mix    = (TH2D*)f->Get("IMnpim_IMnpip_n_mix");                             
                                                                                                 
  cIMnpim_IMnpip_n->cd(3);                                                                       
  //IMnpim_IMnpip_n_rdata->RebinX(4);                                                            
  //IMnpim_IMnpip_n_rdata->RebinY(4);                                                            
  //IMnpim_IMnpip_n_mix->RebinX(4);                                                               
  //IMnpim_IMnpip_n_mix->RebinY(4);                                                               
  IMnpim_IMnpip_n_rdata->SetTitle("");                                                           
  IMnpim_IMnpip_n_rdata->GetXaxis()->SetRangeUser(Slow,Shigh);                                   
  IMnpim_IMnpip_n_rdata->GetYaxis()->SetRangeUser(Slow,Shigh);                                   
  IMnpim_IMnpip_n_rdata->Draw("cont1");                                                          
                                                                                                 
  cIMnpim_IMnpip_n->cd(1);                                                                       
  TH1D* IMnpip_n_rdata = (TH1D*)IMnpim_IMnpip_n_rdata->ProjectionX("IMnpip_n_rdata");            
  TH1D* IMnpip_n_mix = (TH1D*)IMnpim_IMnpip_n_mix->ProjectionX("IMnpip_n_mix");                     
  IMnpip_n_rdata->SetTitle("");                                                                  
  IMnpip_n_rdata->GetXaxis()->SetTitle("");                                                      
  IMnpip_n_rdata->GetXaxis()->SetLabelSize(0);                                                   
  IMnpip_n_rdata->SetMarkerStyle(20);                                                            
  IMnpip_n_rdata->Draw("E1");                                                                    
  IMnpip_n_mix->SetLineColor(6);                                                                  
  IMnpip_n_mix->SetMarkerColor(6);                                                                
  IMnpip_n_mix->SetMarkerStyle(20);                                                               
  IMnpip_n_mix->Draw("E1same");                                                                   
                                                                                                 
  cIMnpim_IMnpip_n->cd(4);                                                                       
  TH1D* IMnpim_n_rdata = (TH1D*)IMnpim_IMnpip_n_rdata->ProjectionY("IMnpim_n_rdata");            
  TH1D* IMnpim_n_mix = (TH1D*)IMnpim_IMnpip_n_mix->ProjectionY("IMnpim_n_mix");                     
  IMnpim_n_rdata->SetTitle("");                                                                  
  IMnpim_n_rdata->GetXaxis()->SetTitle("");                                                      
  TGraphErrors *gr_IMnpim_n_rdata = new TGraphErrors();                                          
  TGraphErrors *gr_IMnpim_n_mix = new TGraphErrors();                                             
  for(int ibin=0;ibin<IMnpim_n_rdata->GetNbinsX();ibin++){                                       
    double cont = IMnpim_n_rdata->GetBinContent(ibin);                                           
    double bincenter = IMnpim_n_rdata->GetBinCenter(ibin);                                       
    double err = IMnpim_n_rdata->GetBinError(ibin);                                              
    gr_IMnpim_n_rdata->SetPoint(ibin,cont,bincenter);                                            
    gr_IMnpim_n_rdata->SetPointError(ibin,err,0);                                                
  }                                                                                              
  for(int ibin=0;ibin<IMnpim_n_mix->GetNbinsX();ibin++){                                          
    double cont = IMnpim_n_mix->GetBinContent(ibin);                                              
    double bincenter = IMnpim_n_mix->GetBinCenter(ibin);                                          
    double err = IMnpim_n_mix->GetBinError(ibin);                                                 
    gr_IMnpim_n_mix->SetPoint(ibin,cont,bincenter);                                               
    gr_IMnpim_n_mix->SetPointError(ibin,err,0);                                                   
  }                                                                                              
  gr_IMnpim_n_rdata->GetYaxis()->SetLabelOffset(999);                                            
  gr_IMnpim_n_rdata->GetYaxis()->SetLabelSize(0);                                                
  gr_IMnpim_n_rdata->GetYaxis()->SetRangeUser(IMnpim_IMnpip_n_rdata->GetYaxis()->GetXmin(),1.7);
  gr_IMnpim_n_rdata->SetLineWidth(2);                                                            
  gr_IMnpim_n_rdata->SetMarkerStyle(20);                                                         
  gr_IMnpim_n_rdata->Draw("AP");                                                                 
  gr_IMnpim_n_mix->SetMarkerStyle(20);                                                            
  gr_IMnpim_n_mix->SetMarkerColor(6);                                                             
  gr_IMnpim_n_mix->SetLineColor(6);                                                               
  gr_IMnpim_n_mix->Draw("P");                                                                     
  
  
  
  TH2D* IMnpim_IMnpip_woK0_wSid_n_data = (TH2D*)f->Get("IMnpim_IMnpip_woK0_wSid_n_data");IMnpim_IMnpip_woK0_wSid_n_data->Print();
  TH2D* IMnpim_IMnpip_wK0_wSid_n_data = (TH2D*)f->Get("IMnpim_IMnpip_wK0_wSid_n_data");IMnpim_IMnpip_wK0_wSid_n_data->Print();
  TH2D* IMnpim_IMnpip_woK0_wSid_n_mix = (TH2D*)f->Get("IMnpim_IMnpip_woK0_wSid_n_mix");IMnpim_IMnpip_woK0_wSid_n_mix->Print();//wo geta
  TH2D* IMnpim_IMnpip_wK0_wSid_n_mix = (TH2D*)f->Get("IMnpim_IMnpip_wK0_wSid_n_mix");IMnpim_IMnpip_wK0_wSid_n_mix->Print();//w geta
  TH2D* IMnpim_IMnpip_wK0_wSid_n_mixgeta = (TH2D*)f->Get("IMnpim_IMnpip_wK0_wSid_n_mixgeta");IMnpim_IMnpip_wK0_wSid_n_mixgeta->Print();
 
  //adding woK0 + wK0
  TH2D* IMnpim_IMnpip_wSid_n_data = (TH2D*)IMnpim_IMnpip_woK0_wSid_n_data->Clone("IMnpim_IMnpip_wSid_n_data");
  IMnpim_IMnpip_wSid_n_data->Add(IMnpim_IMnpip_wK0_wSid_n_data);
  TH2D* IMnpim_IMnpip_wSid_n_mix = (TH2D*)IMnpim_IMnpip_woK0_wSid_n_mix->Clone("IMnpim_IMnpip_wSid_n_mix");
  IMnpim_IMnpip_wSid_n_mix->Add(IMnpim_IMnpip_wK0_wSid_n_mix);
  IMnpim_IMnpip_wSid_n_mix->Add(IMnpim_IMnpip_wK0_wSid_n_mixgeta,1.0);
  IMnpim_IMnpip_wSid_n_mix->Print();
  //divide pipin"X" BG and K0n"X" BG
  TH2D* IMnpim_IMnpip_wSid_n_mixpipinX = (TH2D*)IMnpim_IMnpip_woK0_wSid_n_mix->Clone("IMnpim_IMnpip_wSid_n_mixpipinX");
  IMnpim_IMnpip_wSid_n_mixpipinX->Add(IMnpim_IMnpip_wK0_wSid_n_mixgeta,1.0);
  IMnpim_IMnpip_wSid_n_mixpipinX->Print();
  TH2D* IMnpim_IMnpip_wSid_n_mixK0nX = (TH2D*)IMnpim_IMnpip_wK0_wSid_n_mix->Clone("IMnpim_IMnpip_wSid_n_mixK0nX");
  
  TCanvas *cIMnpim_IMnpip_wSid_n = new TCanvas("cIMnpim_IMnpip_wSid_n","cIMnpim_IMnpip_wSid_n",1000,1000);
  cIMnpim_IMnpip_wSid_n->Divide(2,2,0,0);
  cIMnpim_IMnpip_wSid_n->cd(3);
  IMnpim_IMnpip_wSid_n_data->SetTitle("");
  //IMnpim_IMnpip_wSid_n_data->RebinX(4);
  //IMnpim_IMnpip_wSid_n_data->RebinY(4);
  IMnpim_IMnpip_wSid_n_data->GetXaxis()->SetRangeUser(Slow,Shigh);
  IMnpim_IMnpip_wSid_n_data->GetYaxis()->SetRangeUser(Slow,Shigh);
  //IMnpim_IMnpip_wSid_n_mix->RebinX(4);
  //IMnpim_IMnpip_wSid_n_mix->RebinY(4);
  //IMnpim_IMnpip_wSid_n_mixpipinX->RebinX(4);
  //IMnpim_IMnpip_wSid_n_mixK0nX->RebinX(4);
  //IMnpim_IMnpip_wSid_n_mixpipinX->RebinY(4);
  //IMnpim_IMnpip_wSid_n_mixK0nX->RebinY(4);
  //IMnpim_IMnpip_wSid_n_mix->GetXaxis()->SetRangeUser(1.1,1.3);
  //IMnpim_IMnpip_wSid_n_mix->GetYaxis()->SetRangeUser(1.1,1.3);
  //IMnpim_IMnpip_wSid_n_data->GetZaxis()->SetNdivisions(503);
  //IMnpim_IMnpip_wSid_n_mix->GetZaxis()->SetNdivisions(503);
  IMnpim_IMnpip_wSid_n_data->SetContour(8);
  IMnpim_IMnpip_wSid_n_mix->SetContour(8);
  IMnpim_IMnpip_wSid_n_data->SetMinimum(1);
  IMnpim_IMnpip_wSid_n_data->SetLineWidth(2);
  //IMnpim_IMnpip_wSid_n_data->Draw("cont3");
  IMnpim_IMnpip_wSid_n_data->Draw("colz");
  //IMnpim_IMnpip_wSid_n_data->Draw("col");
  IMnpim_IMnpip_wSid_n_mix->SetLineColor(6);
  IMnpim_IMnpip_wSid_n_mix->SetLineWidth(2);
  IMnpim_IMnpip_wSid_n_mix->SetMinimum(1);
  //IMnpim_IMnpip_wSid_n_mix->SetMaximum(IMnpim_IMnpip_wSid_n_mix->GetMaximum());
  //IMnpim_IMnpip_wSid_n_mix->Draw("cont3same");
  //TPaletteAxis *palette2 = (TPaletteAxis*)IMnpim_IMnpip_wSid_n_data->GetListOfFunctions()->FindObject("palette");
  // the following lines moe the paletter. Choose the values you need for the position.
  //palette2->SetX1NDC(0.78);
  //palette2->SetX2NDC(0.83);
  //palette2->SetY1NDC(0.1);
  //palette2->SetY2NDC(0.90);
   
  
  cIMnpim_IMnpip_wSid_n->cd(1);
  TH1D* IMnpip_wSid_n_data = (TH1D*)IMnpim_IMnpip_wSid_n_data->ProjectionX("IMnpip_wSid_n_data");
  TH1D* IMnpip_wSid_n_mix = (TH1D*)IMnpim_IMnpip_wSid_n_mix->ProjectionX("IMnpip_wSid_n_mix");
  TH1D* IMnpip_wSid_n_mixpipinX = (TH1D*)IMnpim_IMnpip_wSid_n_mixpipinX->ProjectionX("IMnpip_wSid_n_mixpipinX");
  TH1D* IMnpip_wSid_n_mixK0nX = (TH1D*)IMnpim_IMnpip_wSid_n_mixK0nX->ProjectionX("IMnpip_wSid_n_mixK0nX");
  IMnpip_wSid_n_data->GetXaxis()->SetLabelSize(0);
  IMnpip_wSid_n_data->SetMarkerStyle(20);
  IMnpip_wSid_n_data->SetMarkerColor(1);
  IMnpip_wSid_n_data->Draw("E");
  IMnpip_wSid_n_mix->SetLineColor(6);
  IMnpip_wSid_n_mix->SetMarkerStyle(20);
  IMnpip_wSid_n_mix->SetMarkerColor(6);
  IMnpip_wSid_n_mix->Draw("Esame");IMnpip_wSid_n_mix->Print();
  IMnpip_wSid_n_mixpipinX->SetLineColor(2);
  IMnpip_wSid_n_mixpipinX->SetMarkerStyle(20);
  IMnpip_wSid_n_mixpipinX->SetMarkerColor(2);IMnpip_wSid_n_mixpipinX->Print();
  if(showBG)IMnpip_wSid_n_mixpipinX->Draw("Esame");
  IMnpip_wSid_n_mixK0nX->SetLineColor(3);
  IMnpip_wSid_n_mixK0nX->SetMarkerStyle(20);
  IMnpip_wSid_n_mixK0nX->SetMarkerColor(3);
  if(showBG)IMnpip_wSid_n_mixK0nX->Draw("Esame");


  cIMnpim_IMnpip_wSid_n->cd(4);
  TH1D* IMnpim_wSid_n_data = (TH1D*)IMnpim_IMnpip_wSid_n_data->ProjectionY("IMnpim_wSid_n_data");
  TH1D* IMnpim_wSid_n_mix = (TH1D*)IMnpim_IMnpip_wSid_n_mix->ProjectionY("IMnpim_wSid_n_mix");
  TH1D* IMnpim_wSid_n_mixpipinX = (TH1D*)IMnpim_IMnpip_wSid_n_mixpipinX->ProjectionY("IMnpim_wSid_n_mixpipinX");
  TH1D* IMnpim_wSid_n_mixK0nX = (TH1D*)IMnpim_IMnpip_wSid_n_mixK0nX->ProjectionY("IMnpim_wSid_n_mixK0nX");
  IMnpim_wSid_n_data->SetTitle("");
  IMnpim_wSid_n_data->GetXaxis()->SetTitle("");
  
  TGraphErrors *gr_IMnpim_wSid_n_data = new TGraphErrors();
  TGraphErrors *gr_IMnpim_wSid_n_mix = new TGraphErrors();
  TGraphErrors *gr_IMnpim_wSid_n_mixpipinX = new TGraphErrors();
  TGraphErrors *gr_IMnpim_wSid_n_mixK0nX = new TGraphErrors();
  
  for(int ibin=0;ibin<IMnpim_wSid_n_data->GetNbinsX();ibin++){
    double cont = IMnpim_wSid_n_data->GetBinContent(ibin);
    double err = IMnpim_wSid_n_data->GetBinError(ibin);
    double bincenter = IMnpim_wSid_n_data->GetBinCenter(ibin);
    gr_IMnpim_wSid_n_data->SetPoint(ibin,cont, bincenter);
    gr_IMnpim_wSid_n_data->SetPointError(ibin,err,0);
  }
  for(int ibin=0;ibin<IMnpim_wSid_n_mix->GetNbinsX();ibin++){
    double cont = IMnpim_wSid_n_mix->GetBinContent(ibin);
    double err = IMnpim_wSid_n_mix->GetBinError(ibin);
    double bincenter = IMnpim_wSid_n_mix->GetBinCenter(ibin);
    gr_IMnpim_wSid_n_mix->SetPoint(ibin,cont, bincenter);
    gr_IMnpim_wSid_n_mix->SetPointError(ibin,err,0);
  }
  for(int ibin=0;ibin<IMnpim_wSid_n_mixpipinX->GetNbinsX();ibin++){
    double cont = IMnpim_wSid_n_mixpipinX->GetBinContent(ibin);
    double err = IMnpim_wSid_n_mixpipinX->GetBinError(ibin);
    double bincenter = IMnpim_wSid_n_mixpipinX->GetBinCenter(ibin);
    gr_IMnpim_wSid_n_mixpipinX->SetPoint(ibin,cont, bincenter);
    gr_IMnpim_wSid_n_mixpipinX->SetPointError(ibin,err,0);
  }
  for(int ibin=0;ibin<IMnpim_wSid_n_mixK0nX->GetNbinsX();ibin++){
    double cont = IMnpim_wSid_n_mixK0nX->GetBinContent(ibin);
    double err = IMnpim_wSid_n_mixK0nX->GetBinError(ibin);
    double bincenter = IMnpim_wSid_n_mixK0nX->GetBinCenter(ibin);
    gr_IMnpim_wSid_n_mixK0nX->SetPoint(ibin,cont, bincenter);
    gr_IMnpim_wSid_n_mixK0nX->SetPointError(ibin,err,0);
  }
  gr_IMnpim_wSid_n_data->GetYaxis()->SetRangeUser(Slow,Shigh);
  gr_IMnpim_wSid_n_data->SetMarkerStyle(20);
  gr_IMnpim_wSid_n_data->GetYaxis()->SetLabelSize(0);
  gr_IMnpim_wSid_n_data->Draw("AP");
  gr_IMnpim_wSid_n_mix->SetMarkerStyle(20);
  gr_IMnpim_wSid_n_mix->SetMarkerColor(6);
  gr_IMnpim_wSid_n_mix->Draw("P");
  gr_IMnpim_wSid_n_mixpipinX->SetMarkerStyle(20);
  gr_IMnpim_wSid_n_mixpipinX->SetMarkerColor(2);
  if(showBG)gr_IMnpim_wSid_n_mixpipinX->Draw("P");
  gr_IMnpim_wSid_n_mixK0nX->SetMarkerStyle(20);
  gr_IMnpim_wSid_n_mixK0nX->SetMarkerColor(3);
  if(showBG)gr_IMnpim_wSid_n_mixK0nX->Draw("P");
  
  //BG-subtracted 
  TCanvas *cIMnpim_IMnpip_wSid_n_sub = new TCanvas("cIMnpim_IMnpip_wSid_n_sub","cIMnpim_IMnpip_wSid_n_sub",1000,1000);
  cIMnpim_IMnpip_wSid_n_sub->Divide(2,2,0,0);
  cIMnpim_IMnpip_wSid_n_sub->cd(3);
  TH2D* IMnpim_IMnpip_wSid_n_data_sub = (TH2D*)IMnpim_IMnpip_wSid_n_data->Clone("IMnpim_IMnpip_wSid_n_data_sub");
  IMnpim_IMnpip_wSid_n_data_sub->Add(IMnpim_IMnpip_wSid_n_mix,-1.0);
  IMnpim_IMnpip_wSid_n_data_sub->SetMinimum(1);
  IMnpim_IMnpip_wSid_n_data_sub->Draw("col");

  cIMnpim_IMnpip_wSid_n_sub->cd(1);
  TH1D* IMnpip_wSid_n_data_sub = (TH1D*)IMnpim_IMnpip_wSid_n_data_sub->ProjectionX("IMnpip_wSid_n_data_sub");
  IMnpip_wSid_n_data_sub->SetMarkerStyle(20);
  IMnpip_wSid_n_data_sub->Draw();

  cIMnpim_IMnpip_wSid_n_sub->cd(4);
  TH1D* IMnpim_wSid_n_data_sub = (TH1D*)IMnpim_IMnpip_wSid_n_data_sub->ProjectionY("IMnpim_wSid_n_data_sub");
  TGraphErrors *gr_IMnpim_wSid_n_data_sub = new TGraphErrors();
  for(int ibin=0;ibin<IMnpim_wSid_n_data_sub->GetNbinsX();ibin++){
    double cont = IMnpim_wSid_n_data_sub->GetBinContent(ibin);
    double err = IMnpim_wSid_n_data_sub->GetBinError(ibin);
    double bincenter = IMnpim_wSid_n_data_sub->GetBinCenter(ibin);
    gr_IMnpim_wSid_n_data_sub->SetPoint(ibin,cont, bincenter);
    gr_IMnpim_wSid_n_data_sub->SetPointError(ibin,err,0);
  }
  gr_IMnpim_wSid_n_data_sub->SetMarkerStyle(20);
  gr_IMnpim_wSid_n_data_sub->GetYaxis()->SetRangeUser(Slow,Shigh);
  gr_IMnpim_wSid_n_data_sub->Draw("AP");
  
 
  //
  //signal + BG check MMnmiss vs IMpippim wSid 
  //
  TH2D* MMnmiss_IMpippim_woK0_wSid_data = (TH2D*)f->Get("MMnmiss_IMpippim_woK0_wSid_data");MMnmiss_IMpippim_woK0_wSid_data->Print();
  TH2D* MMnmiss_IMpippim_wK0_wSid_data = (TH2D*)f->Get("MMnmiss_IMpippim_wK0_wSid_data");MMnmiss_IMpippim_wK0_wSid_data->Print();
  TH2D* MMnmiss_IMpippim_woK0_wSid_mix = (TH2D*)f->Get("MMnmiss_IMpippim_woK0_wSid_mix");MMnmiss_IMpippim_woK0_wSid_mix->Print();//wo geta
  TH2D* MMnmiss_IMpippim_wK0_wSid_mix = (TH2D*)f->Get("MMnmiss_IMpippim_wK0_wSid_mix");MMnmiss_IMpippim_wK0_wSid_mix->Print();//w geta
  TH2D* MMnmiss_IMpippim_wK0_wSid_mixgeta = (TH2D*)f->Get("MMnmiss_IMpippim_wK0_wSid_mixgeta");MMnmiss_IMpippim_wK0_wSid_mixgeta->Print();
  
  //adding woK0 + wK0
  TH2D* MMnmiss_IMpippim_wSid_data = (TH2D*)MMnmiss_IMpippim_woK0_wSid_data->Clone("MMnmiss_IMpippim_wSid_data");
  MMnmiss_IMpippim_wSid_data->Add(MMnmiss_IMpippim_wK0_wSid_data);
  TH2D* MMnmiss_IMpippim_wSid_mix = (TH2D*)MMnmiss_IMpippim_woK0_wSid_mix->Clone("MMnmiss_IMpippim_wSid_mix");
  MMnmiss_IMpippim_wSid_mix->Add(MMnmiss_IMpippim_wK0_wSid_mix);
  MMnmiss_IMpippim_wSid_mix->Add(MMnmiss_IMpippim_wK0_wSid_mixgeta,1.0);
  MMnmiss_IMpippim_wSid_mix->Print();
  
  TCanvas *cMMnmiss_IMpippim_wSid = new TCanvas("cMMnmiss_IMpippim_wSid","cMMnmiss_IMpippim_wSid",1000,1000);
  cMMnmiss_IMpippim_wSid->Divide(2,2,0,0);
  cMMnmiss_IMpippim_wSid->cd(3);
  MMnmiss_IMpippim_wSid_data->SetTitle("");
  MMnmiss_IMpippim_wSid_data->RebinX(4);
  MMnmiss_IMpippim_wSid_mix->RebinX(4);
  MMnmiss_IMpippim_wSid_data->SetContour(8);
  MMnmiss_IMpippim_wSid_mix->SetContour(8);
  MMnmiss_IMpippim_wSid_data->SetMinimum(1);
  MMnmiss_IMpippim_wSid_data->SetLineWidth(2);
  MMnmiss_IMpippim_wSid_data->Draw("colz");
  MMnmiss_IMpippim_wSid_mix->SetLineColor(6);
  MMnmiss_IMpippim_wSid_mix->SetLineWidth(2);
  MMnmiss_IMpippim_wSid_mix->SetMinimum(1);
 
 
  cMMnmiss_IMpippim_wSid->cd(1);
  TH1D* IMpippim_wSid_data = (TH1D*)MMnmiss_IMpippim_wSid_data->ProjectionX("IMpippim_wSid_data");
  TH1D* IMpippim_wSid_mix = (TH1D*)MMnmiss_IMpippim_wSid_mix->ProjectionX("IMpippim_wSid_mix");
  IMpippim_wSid_data->GetXaxis()->SetLabelSize(0);
  IMpippim_wSid_data->SetMarkerStyle(20);
  IMpippim_wSid_data->SetMarkerColor(1);
  IMpippim_wSid_data->Draw("E");
  IMpippim_wSid_mix->SetLineColor(6);
  IMpippim_wSid_mix->SetMarkerStyle(20);
  IMpippim_wSid_mix->SetMarkerColor(6);
  IMpippim_wSid_mix->Draw("Esame");IMpippim_wSid_mix->Print();
 
 
  cMMnmiss_IMpippim_wSid->cd(4);
  TH1D* MMnmiss_wSid_data = (TH1D*)MMnmiss_IMpippim_wSid_data->ProjectionY("MMnmiss_wSid_data");
  TH1D* MMnmiss_wSid_mix = (TH1D*)MMnmiss_IMpippim_wSid_mix->ProjectionY("MMnmiss_wSid_mix");
  MMnmiss_wSid_data->SetTitle("");
  MMnmiss_wSid_data->GetXaxis()->SetTitle("");
  
  TGraphErrors *gr_MMnmiss_wSid_data = new TGraphErrors();
  TGraphErrors *gr_MMnmiss_wSid_mix = new TGraphErrors();
  
  for(int ibin=0;ibin<MMnmiss_wSid_data->GetNbinsX();ibin++){
    double cont = MMnmiss_wSid_data->GetBinContent(ibin);
    double err = MMnmiss_wSid_data->GetBinError(ibin);
    double bincenter = MMnmiss_wSid_data->GetBinCenter(ibin);
    gr_MMnmiss_wSid_data->SetPoint(ibin,cont, bincenter);
    gr_MMnmiss_wSid_data->SetPointError(ibin,err,0);
  }
  for(int ibin=0;ibin<MMnmiss_wSid_mix->GetNbinsX();ibin++){
    double cont = MMnmiss_wSid_mix->GetBinContent(ibin);
    double err = MMnmiss_wSid_mix->GetBinError(ibin);
    double bincenter = MMnmiss_wSid_mix->GetBinCenter(ibin);
    gr_MMnmiss_wSid_mix->SetPoint(ibin,cont, bincenter);
    gr_MMnmiss_wSid_mix->SetPointError(ibin,err,0);
  }
  gr_MMnmiss_wSid_data->SetMarkerStyle(20);
  gr_MMnmiss_wSid_data->GetYaxis()->SetLabelSize(0);
  gr_MMnmiss_wSid_data->Draw("AP");
  gr_MMnmiss_wSid_mix->SetMarkerStyle(20);
  gr_MMnmiss_wSid_mix->SetMarkerColor(6);
  gr_MMnmiss_wSid_mix->Draw("P");
 


  //signal check MMnmiss vs IMpippim
  TH2D* MMnmiss_IMpippim_woK0_wSid_n_data = (TH2D*)f->Get("MMnmiss_IMpippim_woK0_wSid_n_data");MMnmiss_IMpippim_woK0_wSid_n_data->Print();
  TH2D* MMnmiss_IMpippim_wK0_wSid_n_data = (TH2D*)f->Get("MMnmiss_IMpippim_wK0_wSid_n_data");MMnmiss_IMpippim_wK0_wSid_n_data->Print();
  TH2D* MMnmiss_IMpippim_woK0_wSid_n_mix = (TH2D*)f->Get("MMnmiss_IMpippim_woK0_wSid_n_mix");MMnmiss_IMpippim_woK0_wSid_n_mix->Print();//wo geta
  TH2D* MMnmiss_IMpippim_wK0_wSid_n_mix = (TH2D*)f->Get("MMnmiss_IMpippim_wK0_wSid_n_mix");MMnmiss_IMpippim_wK0_wSid_n_mix->Print();//w geta
  TH2D* MMnmiss_IMpippim_wK0_wSid_n_mixgeta = (TH2D*)f->Get("MMnmiss_IMpippim_wK0_wSid_n_mixgeta");MMnmiss_IMpippim_wK0_wSid_n_mixgeta->Print();
  
  //adding woK0 + wK0
  TH2D* MMnmiss_IMpippim_wSid_n_data = (TH2D*)MMnmiss_IMpippim_woK0_wSid_n_data->Clone("MMnmiss_IMpippim_wSid_n_data");
  MMnmiss_IMpippim_wSid_n_data->Add(MMnmiss_IMpippim_wK0_wSid_n_data);
  TH2D* MMnmiss_IMpippim_wSid_n_mix = (TH2D*)MMnmiss_IMpippim_woK0_wSid_n_mix->Clone("MMnmiss_IMpippim_wSid_n_mix");
  MMnmiss_IMpippim_wSid_n_mix->Add(MMnmiss_IMpippim_wK0_wSid_n_mix);
  MMnmiss_IMpippim_wSid_n_mix->Add(MMnmiss_IMpippim_wK0_wSid_n_mixgeta,1.0);
  MMnmiss_IMpippim_wSid_n_mix->Print();
  //divide pipin"X" BG and K0n"X" BG
  TH2D* MMnmiss_IMpippim_wSid_n_mixpipinX = (TH2D*)MMnmiss_IMpippim_woK0_wSid_n_mix->Clone("MMnmiss_IMpippim_wSid_n_mixpipinX");
  MMnmiss_IMpippim_wSid_n_mixpipinX->Add(MMnmiss_IMpippim_wK0_wSid_n_mixgeta,1.0);
  MMnmiss_IMpippim_wSid_n_mixpipinX->Print();
  TH2D* MMnmiss_IMpippim_wSid_n_mixK0nX = (TH2D*)MMnmiss_IMpippim_wK0_wSid_n_mix->Clone("MMnmiss_IMpippim_wSid_n_mixK0nX");
  
  TCanvas *cMMnmiss_IMpippim_wSid_n = new TCanvas("cMMnmiss_IMpippim_wSid_n","cMMnmiss_IMpippim_wSid_n",1000,1000);
  cMMnmiss_IMpippim_wSid_n->Divide(2,2,0,0);
  cMMnmiss_IMpippim_wSid_n->cd(3);
  MMnmiss_IMpippim_wSid_n_data->SetTitle("");
  MMnmiss_IMpippim_wSid_n_data->RebinX(4);
  //MMnmiss_IMpippim_wSid_n_data->RebinY(4);
  //MMnmiss_IMpippim_wSid_n_data->GetXaxis()->SetRangeUser(1.1,1.3);
  //MMnmiss_IMpippim_wSid_n_data->GetYaxis()->SetRangeUser(1.1,1.3);
  MMnmiss_IMpippim_wSid_n_mix->RebinX(4);
  //MMnmiss_IMpippim_wSid_n_mix->RebinY(4);
  MMnmiss_IMpippim_wSid_n_mixpipinX->RebinX(4);
  MMnmiss_IMpippim_wSid_n_mixK0nX->RebinX(4);
  //MMnmiss_IMpippim_wSid_n_mixpipinX->RebinY(4);
  //MMnmiss_IMpippim_wSid_n_mixK0nX->RebinY(4);
  //MMnmiss_IMpippim_wSid_n_mix->GetXaxis()->SetRangeUser(1.1,1.3);
  //MMnmiss_IMpippim_wSid_n_mix->GetYaxis()->SetRangeUser(1.1,1.3);
  MMnmiss_IMpippim_wSid_n_data->SetContour(8);
  MMnmiss_IMpippim_wSid_n_mix->SetContour(8);
  MMnmiss_IMpippim_wSid_n_data->SetMinimum(1);
  MMnmiss_IMpippim_wSid_n_data->SetLineWidth(2);
  //MMnmiss_IMpippim_wSid_n_data->Draw("cont3");
  MMnmiss_IMpippim_wSid_n_data->Draw("colz");
  //MMnmiss_IMpippim_wSid_n_data->Draw("col");
  MMnmiss_IMpippim_wSid_n_mix->SetLineColor(6);
  MMnmiss_IMpippim_wSid_n_mix->SetLineWidth(2);
  MMnmiss_IMpippim_wSid_n_mix->SetMinimum(1);
   
  
  cMMnmiss_IMpippim_wSid_n->cd(1);
  TH1D* IMpippim_wSid_n_data = (TH1D*)MMnmiss_IMpippim_wSid_n_data->ProjectionX("IMpippim_wSid_n_data");
  TH1D* IMpippim_wSid_n_mix = (TH1D*)MMnmiss_IMpippim_wSid_n_mix->ProjectionX("IMpippim_wSid_n_mix");
  TH1D* IMpippim_wSid_n_mixpipinX = (TH1D*)MMnmiss_IMpippim_wSid_n_mixpipinX->ProjectionX("IMpippim_wSid_n_mixpipinX");
  TH1D* IMpippim_wSid_n_mixK0nX = (TH1D*)MMnmiss_IMpippim_wSid_n_mixK0nX->ProjectionX("IMpippim_wSid_n_mixK0nX");
  IMpippim_wSid_n_data->GetXaxis()->SetLabelSize(0);
  IMpippim_wSid_n_data->SetMarkerStyle(20);
  IMpippim_wSid_n_data->SetMarkerColor(1);
  IMpippim_wSid_n_data->Draw("E");
  IMpippim_wSid_n_mix->SetLineColor(6);
  IMpippim_wSid_n_mix->SetMarkerStyle(20);
  IMpippim_wSid_n_mix->SetMarkerColor(6);
  IMpippim_wSid_n_mix->Draw("Esame");IMpippim_wSid_n_mix->Print();
  IMpippim_wSid_n_mixpipinX->SetLineColor(2);
  IMpippim_wSid_n_mixpipinX->SetMarkerStyle(20);
  IMpippim_wSid_n_mixpipinX->SetMarkerColor(2);IMpippim_wSid_n_mixpipinX->Print();
  if(showBG)IMpippim_wSid_n_mixpipinX->Draw("Esame");
  IMpippim_wSid_n_mixK0nX->SetLineColor(3);
  IMpippim_wSid_n_mixK0nX->SetMarkerStyle(20);
  IMpippim_wSid_n_mixK0nX->SetMarkerColor(3);
  if(showBG)IMpippim_wSid_n_mixK0nX->Draw("Esame");


  cMMnmiss_IMpippim_wSid_n->cd(4);
  TH1D* MMnmiss_wSid_n_data = (TH1D*)MMnmiss_IMpippim_wSid_n_data->ProjectionY("MMnmiss_wSid_n_data");
  TH1D* MMnmiss_wSid_n_mix = (TH1D*)MMnmiss_IMpippim_wSid_n_mix->ProjectionY("MMnmiss_wSid_n_mix");
  TH1D* MMnmiss_wSid_n_mixpipinX = (TH1D*)MMnmiss_IMpippim_wSid_n_mixpipinX->ProjectionY("MMnmiss_wSid_n_mixpipinX");
  TH1D* MMnmiss_wSid_n_mixK0nX = (TH1D*)MMnmiss_IMpippim_wSid_n_mixK0nX->ProjectionY("MMnmiss_wSid_n_mixK0nX");
  MMnmiss_wSid_n_data->SetTitle("");
  MMnmiss_wSid_n_data->GetXaxis()->SetTitle("");
  
  TGraphErrors *gr_MMnmiss_wSid_n_data = new TGraphErrors();
  TGraphErrors *gr_MMnmiss_wSid_n_mix = new TGraphErrors();
  TGraphErrors *gr_MMnmiss_wSid_n_mixpipinX = new TGraphErrors();
  TGraphErrors *gr_MMnmiss_wSid_n_mixK0nX = new TGraphErrors();
  
  for(int ibin=0;ibin<MMnmiss_wSid_n_data->GetNbinsX();ibin++){
    double cont = MMnmiss_wSid_n_data->GetBinContent(ibin);
    double err = MMnmiss_wSid_n_data->GetBinError(ibin);
    double bincenter = MMnmiss_wSid_n_data->GetBinCenter(ibin);
    gr_MMnmiss_wSid_n_data->SetPoint(ibin,cont, bincenter);
    gr_MMnmiss_wSid_n_data->SetPointError(ibin,err,0);
  }
  for(int ibin=0;ibin<MMnmiss_wSid_n_mix->GetNbinsX();ibin++){
    double cont = MMnmiss_wSid_n_mix->GetBinContent(ibin);
    double err = MMnmiss_wSid_n_mix->GetBinError(ibin);
    double bincenter = MMnmiss_wSid_n_mix->GetBinCenter(ibin);
    gr_MMnmiss_wSid_n_mix->SetPoint(ibin,cont, bincenter);
    gr_MMnmiss_wSid_n_mix->SetPointError(ibin,err,0);
  }
  for(int ibin=0;ibin<MMnmiss_wSid_n_mixpipinX->GetNbinsX();ibin++){
    double cont = MMnmiss_wSid_n_mixpipinX->GetBinContent(ibin);
    double err = MMnmiss_wSid_n_mixpipinX->GetBinError(ibin);
    double bincenter = MMnmiss_wSid_n_mixpipinX->GetBinCenter(ibin);
    gr_MMnmiss_wSid_n_mixpipinX->SetPoint(ibin,cont, bincenter);
    gr_MMnmiss_wSid_n_mixpipinX->SetPointError(ibin,err,0);
  }
  for(int ibin=0;ibin<MMnmiss_wSid_n_mixK0nX->GetNbinsX();ibin++){
    double cont = MMnmiss_wSid_n_mixK0nX->GetBinContent(ibin);
    double err = MMnmiss_wSid_n_mixK0nX->GetBinError(ibin);
    double bincenter = MMnmiss_wSid_n_mixK0nX->GetBinCenter(ibin);
    gr_MMnmiss_wSid_n_mixK0nX->SetPoint(ibin,cont, bincenter);
    gr_MMnmiss_wSid_n_mixK0nX->SetPointError(ibin,err,0);
  }
  //gr_MMnmiss_wSid_n_data->GetYaxis()->SetRangeUser(1,1.3);
  gr_MMnmiss_wSid_n_data->SetMarkerStyle(20);
  gr_MMnmiss_wSid_n_data->GetYaxis()->SetLabelSize(0);
  gr_MMnmiss_wSid_n_data->Draw("AP");
  gr_MMnmiss_wSid_n_mix->SetMarkerStyle(20);
  gr_MMnmiss_wSid_n_mix->SetMarkerColor(6);
  gr_MMnmiss_wSid_n_mix->Draw("P");
  gr_MMnmiss_wSid_n_mixpipinX->SetMarkerStyle(20);
  gr_MMnmiss_wSid_n_mixpipinX->SetMarkerColor(2);
  if(showBG)gr_MMnmiss_wSid_n_mixpipinX->Draw("P");
  gr_MMnmiss_wSid_n_mixK0nX->SetMarkerStyle(20);
  gr_MMnmiss_wSid_n_mixK0nX->SetMarkerColor(3);
  if(showBG)gr_MMnmiss_wSid_n_mixK0nX->Draw("P");

  //BG-subtracted 
  TCanvas *cMMnmiss_IMpippim_wSid_n_sub = new TCanvas("cMMnmiss_IMpippim_wSid_n_sub","cMMnmiss_IMpippim_wSid_n_sub",1000,1000);
  cMMnmiss_IMpippim_wSid_n_sub->Divide(2,2,0,0);
  cMMnmiss_IMpippim_wSid_n_sub->cd(3);
  TH2D* MMnmiss_IMpippim_wSid_n_data_sub = (TH2D*)MMnmiss_IMpippim_wSid_n_data->Clone("MMnmiss_IMpippim_wSid_n_data_sub");
  MMnmiss_IMpippim_wSid_n_data_sub->Add(MMnmiss_IMpippim_wSid_n_mix,-1.0);
  MMnmiss_IMpippim_wSid_n_data_sub->SetMinimum(1);
  MMnmiss_IMpippim_wSid_n_data_sub->Draw("col");

  cMMnmiss_IMpippim_wSid_n_sub->cd(1);
  TH1D* IMpippim_wSid_n_data_sub = (TH1D*)MMnmiss_IMpippim_wSid_n_data_sub->ProjectionX("IMpippim_wSid_n_data_sub");
  IMpippim_wSid_n_data_sub->SetMarkerStyle(20);
  IMpippim_wSid_n_data_sub->Draw();

  cMMnmiss_IMpippim_wSid_n_sub->cd(4);
  TH1D* MMnmiss_wSid_n_data_sub = (TH1D*)MMnmiss_IMpippim_wSid_n_data_sub->ProjectionY("MMnmiss_wSid_n_data_sub");
  TGraphErrors *gr_MMnmiss_wSid_n_data_sub = new TGraphErrors();
  for(int ibin=0;ibin<MMnmiss_wSid_n_data_sub->GetNbinsX();ibin++){
    double cont = MMnmiss_wSid_n_data_sub->GetBinContent(ibin);
    double err = MMnmiss_wSid_n_data_sub->GetBinError(ibin);
    double bincenter = MMnmiss_wSid_n_data_sub->GetBinCenter(ibin);
    gr_MMnmiss_wSid_n_data_sub->SetPoint(ibin,cont, bincenter);
    gr_MMnmiss_wSid_n_data_sub->SetPointError(ibin,err,0);
  }
  gr_MMnmiss_wSid_n_data_sub->SetMarkerStyle(20);
  gr_MMnmiss_wSid_n_data_sub->GetYaxis()->SetRangeUser(0,1.5);
  gr_MMnmiss_wSid_n_data_sub->Draw("AP");
  

  //signam check q vs nmom
  TH2D* q_nmom_woK0_wSid_n_data = (TH2D*)f->Get("q_nmom_woK0_wSid_n_data");q_nmom_woK0_wSid_n_data->Print();
  TH2D* q_nmom_wK0_wSid_n_data = (TH2D*)f->Get("q_nmom_wK0_wSid_n_data");q_nmom_wK0_wSid_n_data->Print();
  TH2D* q_nmom_woK0_wSid_n_mix = (TH2D*)f->Get("q_nmom_woK0_wSid_n_mix");q_nmom_woK0_wSid_n_mix->Print();//wo geta
  TH2D* q_nmom_wK0_wSid_n_mix = (TH2D*)f->Get("q_nmom_wK0_wSid_n_mix");q_nmom_wK0_wSid_n_mix->Print();//w geta
  TH2D* q_nmom_wK0_wSid_n_mixgeta = (TH2D*)f->Get("q_nmom_wK0_wSid_n_mixgeta");q_nmom_wK0_wSid_n_mixgeta->Print();
  
  //adding woK0 + wK0
  TH2D* q_nmom_wSid_n_data = (TH2D*)q_nmom_woK0_wSid_n_data->Clone("q_nmom_wSid_n_data");
  q_nmom_wSid_n_data->Add(q_nmom_wK0_wSid_n_data);
  TH2D* q_nmom_wSid_n_mix = (TH2D*)q_nmom_woK0_wSid_n_mix->Clone("q_nmom_wSid_n_mix");
  q_nmom_wSid_n_mix->Add(q_nmom_wK0_wSid_n_mix);
  q_nmom_wSid_n_mix->Add(q_nmom_wK0_wSid_n_mixgeta,1.0);
  q_nmom_wSid_n_mix->Print();
  //divide pipin"X" BG and K0n"X" BG
  TH2D* q_nmom_wSid_n_mixpipinX = (TH2D*)q_nmom_woK0_wSid_n_mix->Clone("q_nmom_wSid_n_mixpipinX");
  q_nmom_wSid_n_mixpipinX->Add(q_nmom_wK0_wSid_n_mixgeta,1.0);
  q_nmom_wSid_n_mixpipinX->Print();
  TH2D* q_nmom_wSid_n_mixK0nX = (TH2D*)q_nmom_wK0_wSid_n_mix->Clone("q_nmom_wSid_n_mixK0nX");
  
  TCanvas *cq_nmom_wSid_n = new TCanvas("cq_nmom_wSid_n","cq_nmom_wSid_n",1000,1000);
  cq_nmom_wSid_n->Divide(2,2,0,0);
  cq_nmom_wSid_n->cd(3);
  q_nmom_wSid_n_data->SetTitle("");
  q_nmom_wSid_n_data->RebinX(4);
  //q_nmom_wSid_n_data->RebinY(4);
  //q_nmom_wSid_n_data->GetXaxis()->SetRangeUser(1.1,1.3);
  //q_nmom_wSid_n_data->GetYaxis()->SetRangeUser(1.1,1.3);
  q_nmom_wSid_n_mix->RebinX(4);
  //q_nmom_wSid_n_mix->RebinY(4);
  q_nmom_wSid_n_mixpipinX->RebinX(4);
  q_nmom_wSid_n_mixK0nX->RebinX(4);
  //q_nmom_wSid_n_mixpipinX->RebinY(4);
  //q_nmom_wSid_n_mixK0nX->RebinY(4);
  //q_nmom_wSid_n_mix->GetXaxis()->SetRangeUser(1.1,1.3);
  //q_nmom_wSid_n_mix->GetYaxis()->SetRangeUser(1.1,1.3);
  q_nmom_wSid_n_data->SetContour(8);
  q_nmom_wSid_n_mix->SetContour(8);
  q_nmom_wSid_n_data->SetMinimum(1);
  q_nmom_wSid_n_data->SetLineWidth(2);
  //q_nmom_wSid_n_data->Draw("cont3");
  q_nmom_wSid_n_data->Draw("colz");
  //q_nmom_wSid_n_data->Draw("col");
  q_nmom_wSid_n_mix->SetLineColor(6);
  q_nmom_wSid_n_mix->SetLineWidth(2);
  q_nmom_wSid_n_mix->SetMinimum(1);
   
  
  cq_nmom_wSid_n->cd(1);
  TH1D* nmom_wSid_n_data = (TH1D*)q_nmom_wSid_n_data->ProjectionX("nmom_wSid_n_data");
  TH1D* nmom_wSid_n_mix = (TH1D*)q_nmom_wSid_n_mix->ProjectionX("nmom_wSid_n_mix");
  TH1D* nmom_wSid_n_mixpipinX = (TH1D*)q_nmom_wSid_n_mixpipinX->ProjectionX("nmom_wSid_n_mixpipinX");
  TH1D* nmom_wSid_n_mixK0nX = (TH1D*)q_nmom_wSid_n_mixK0nX->ProjectionX("nmom_wSid_n_mixK0nX");
  nmom_wSid_n_data->GetXaxis()->SetLabelSize(0);
  nmom_wSid_n_data->SetMarkerStyle(20);
  nmom_wSid_n_data->SetMarkerColor(1);
  nmom_wSid_n_data->Draw("E");
  nmom_wSid_n_mix->SetLineColor(6);
  nmom_wSid_n_mix->SetMarkerStyle(20);
  nmom_wSid_n_mix->SetMarkerColor(6);
  nmom_wSid_n_mix->Draw("Esame");nmom_wSid_n_mix->Print();
  nmom_wSid_n_mixpipinX->SetLineColor(2);
  nmom_wSid_n_mixpipinX->SetMarkerStyle(20);
  nmom_wSid_n_mixpipinX->SetMarkerColor(2);nmom_wSid_n_mixpipinX->Print();
  if(showBG)nmom_wSid_n_mixpipinX->Draw("Esame");
  nmom_wSid_n_mixK0nX->SetLineColor(3);
  nmom_wSid_n_mixK0nX->SetMarkerStyle(20);
  nmom_wSid_n_mixK0nX->SetMarkerColor(3);
  if(showBG)nmom_wSid_n_mixK0nX->Draw("Esame");


  cq_nmom_wSid_n->cd(4);
  TH1D* q_wSid_n_data = (TH1D*)q_nmom_wSid_n_data->ProjectionY("q_wSid_n_data");
  TH1D* q_wSid_n_mix = (TH1D*)q_nmom_wSid_n_mix->ProjectionY("q_wSid_n_mix");
  TH1D* q_wSid_n_mixpipinX = (TH1D*)q_nmom_wSid_n_mixpipinX->ProjectionY("q_wSid_n_mixpipinX");
  TH1D* q_wSid_n_mixK0nX = (TH1D*)q_nmom_wSid_n_mixK0nX->ProjectionY("q_wSid_n_mixK0nX");
  q_wSid_n_data->SetTitle("");
  q_wSid_n_data->GetXaxis()->SetTitle("");
  
  TGraphErrors *gr_q_wSid_n_data = new TGraphErrors();
  TGraphErrors *gr_q_wSid_n_mix = new TGraphErrors();
  TGraphErrors *gr_q_wSid_n_mixpipinX = new TGraphErrors();
  TGraphErrors *gr_q_wSid_n_mixK0nX = new TGraphErrors();
  
  for(int ibin=0;ibin<q_wSid_n_data->GetNbinsX();ibin++){
    double cont = q_wSid_n_data->GetBinContent(ibin);
    double err = q_wSid_n_data->GetBinError(ibin);
    double bincenter = q_wSid_n_data->GetBinCenter(ibin);
    gr_q_wSid_n_data->SetPoint(ibin,cont, bincenter);
    gr_q_wSid_n_data->SetPointError(ibin,err,0);
  }
  for(int ibin=0;ibin<q_wSid_n_mix->GetNbinsX();ibin++){
    double cont = q_wSid_n_mix->GetBinContent(ibin);
    double err = q_wSid_n_mix->GetBinError(ibin);
    double bincenter = q_wSid_n_mix->GetBinCenter(ibin);
    gr_q_wSid_n_mix->SetPoint(ibin,cont, bincenter);
    gr_q_wSid_n_mix->SetPointError(ibin,err,0);
  }
  for(int ibin=0;ibin<q_wSid_n_mixpipinX->GetNbinsX();ibin++){
    double cont = q_wSid_n_mixpipinX->GetBinContent(ibin);
    double err = q_wSid_n_mixpipinX->GetBinError(ibin);
    double bincenter = q_wSid_n_mixpipinX->GetBinCenter(ibin);
    gr_q_wSid_n_mixpipinX->SetPoint(ibin,cont, bincenter);
    gr_q_wSid_n_mixpipinX->SetPointError(ibin,err,0);
  }
  for(int ibin=0;ibin<q_wSid_n_mixK0nX->GetNbinsX();ibin++){
    double cont = q_wSid_n_mixK0nX->GetBinContent(ibin);
    double err = q_wSid_n_mixK0nX->GetBinError(ibin);
    double bincenter = q_wSid_n_mixK0nX->GetBinCenter(ibin);
    gr_q_wSid_n_mixK0nX->SetPoint(ibin,cont, bincenter);
    gr_q_wSid_n_mixK0nX->SetPointError(ibin,err,0);
  }
  //gr_q_wSid_n_data->GetYaxis()->SetRangeUser(1,1.3);
  gr_q_wSid_n_data->SetMarkerStyle(20);
  gr_q_wSid_n_data->GetYaxis()->SetLabelSize(0);
  gr_q_wSid_n_data->Draw("AP");
  gr_q_wSid_n_mix->SetMarkerStyle(20);
  gr_q_wSid_n_mix->SetMarkerColor(6);
  gr_q_wSid_n_mix->Draw("P");
  gr_q_wSid_n_mixpipinX->SetMarkerStyle(20);
  gr_q_wSid_n_mixpipinX->SetMarkerColor(2);
  if(showBG)gr_q_wSid_n_mixpipinX->Draw("P");
  gr_q_wSid_n_mixK0nX->SetMarkerStyle(20);
  gr_q_wSid_n_mixK0nX->SetMarkerColor(3);
  if(showBG)gr_q_wSid_n_mixK0nX->Draw("P");

  //BG-subtracted 
  TCanvas *cq_nmom_wSid_n_sub = new TCanvas("cq_nmom_wSid_n_sub","cq_nmom_wSid_n_sub",1000,1000);
  cq_nmom_wSid_n_sub->Divide(2,2,0,0);
  cq_nmom_wSid_n_sub->cd(3);
  TH2D* q_nmom_wSid_n_data_sub = (TH2D*)q_nmom_wSid_n_data->Clone("q_nmom_wSid_n_data_sub");
  q_nmom_wSid_n_data_sub->Add(q_nmom_wSid_n_mix,-1.0);
  q_nmom_wSid_n_data_sub->SetMinimum(1);
  q_nmom_wSid_n_data_sub->Draw("col");

  cq_nmom_wSid_n_sub->cd(1);
  TH1D* nmom_wSid_n_data_sub = (TH1D*)q_nmom_wSid_n_data_sub->ProjectionX("nmom_wSid_n_data_sub");
  nmom_wSid_n_data_sub->SetMarkerStyle(20);
  nmom_wSid_n_data_sub->Draw();

  cq_nmom_wSid_n_sub->cd(4);
  TH1D* q_wSid_n_data_sub = (TH1D*)q_nmom_wSid_n_data_sub->ProjectionY("q_wSid_n_data_sub");
  TGraphErrors *gr_q_wSid_n_data_sub = new TGraphErrors();
  for(int ibin=0;ibin<q_wSid_n_data_sub->GetNbinsX();ibin++){
    double cont = q_wSid_n_data_sub->GetBinContent(ibin);
    double err = q_wSid_n_data_sub->GetBinError(ibin);
    double bincenter = q_wSid_n_data_sub->GetBinCenter(ibin);
    gr_q_wSid_n_data_sub->SetPoint(ibin,cont, bincenter);
    gr_q_wSid_n_data_sub->SetPointError(ibin,err,0);
  }
  gr_q_wSid_n_data_sub->SetMarkerStyle(20);
  gr_q_wSid_n_data_sub->GetYaxis()->SetRangeUser(0,1.5);
  gr_q_wSid_n_data_sub->Draw("AP");

  //
  //signal check q vs IMnpipi
  //
  TH2D* q_IMnpipi_woK0_wSid_n_data = (TH2D*)f->Get("q_IMnpipi_woK0_wSid_n_data");q_IMnpipi_woK0_wSid_n_data->Print();
  TH2D* q_IMnpipi_wK0_wSid_n_data = (TH2D*)f->Get("q_IMnpipi_wK0_wSid_n_data");q_IMnpipi_wK0_wSid_n_data->Print();
  TH2D* q_IMnpipi_woK0_wSid_n_mix = (TH2D*)f->Get("q_IMnpipi_woK0_wSid_n_mix");q_IMnpipi_woK0_wSid_n_mix->Print();//wo geta
  TH2D* q_IMnpipi_wK0_wSid_n_mix = (TH2D*)f->Get("q_IMnpipi_wK0_wSid_n_mix");q_IMnpipi_wK0_wSid_n_mix->Print();//w geta
  TH2D* q_IMnpipi_wK0_wSid_n_mixgeta = (TH2D*)f->Get("q_IMnpipi_wK0_wSid_n_mixgeta");q_IMnpipi_wK0_wSid_n_mixgeta->Print();
  
  //adding woK0 + wK0
  TH2D* q_IMnpipi_wSid_n_data = (TH2D*)q_IMnpipi_woK0_wSid_n_data->Clone("q_IMnpipi_wSid_n_data");
  q_IMnpipi_wSid_n_data->Add(q_IMnpipi_wK0_wSid_n_data);
  TH2D* q_IMnpipi_wSid_n_mix = (TH2D*)q_IMnpipi_woK0_wSid_n_mix->Clone("q_IMnpipi_wSid_n_mix");
  q_IMnpipi_wSid_n_mix->Add(q_IMnpipi_wK0_wSid_n_mix);
  q_IMnpipi_wSid_n_mix->Add(q_IMnpipi_wK0_wSid_n_mixgeta,1.0);
  q_IMnpipi_wSid_n_mix->Print();
  //divide pipin"X" BG and K0n"X" BG
  TH2D* q_IMnpipi_wSid_n_mixpipinX = (TH2D*)q_IMnpipi_woK0_wSid_n_mix->Clone("q_IMnpipi_wSid_n_mixpipinX");
  q_IMnpipi_wSid_n_mixpipinX->Add(q_IMnpipi_wK0_wSid_n_mixgeta,1.0);
  q_IMnpipi_wSid_n_mixpipinX->Print();
  TH2D* q_IMnpipi_wSid_n_mixK0nX = (TH2D*)q_IMnpipi_wK0_wSid_n_mix->Clone("q_IMnpipi_wSid_n_mixK0nX");
  
  TCanvas *cq_IMnpipi_wSid_n = new TCanvas("cq_IMnpipi_wSid_n","cq_IMnpipi_wSid_n",1000,1000);
  cq_IMnpipi_wSid_n->Divide(2,2,0,0);
  cq_IMnpipi_wSid_n->cd(3);
  q_IMnpipi_wSid_n_data->SetTitle("");
  q_IMnpipi_wSid_n_data->RebinX(2);
  //q_IMnpipi_wSid_n_data->RebinY(4);
  //q_IMnpipi_wSid_n_data->GetXaxis()->SetRangeUser(1.1,1.3);
  //q_IMnpipi_wSid_n_data->GetYaxis()->SetRangeUser(1.1,1.3);
  q_IMnpipi_wSid_n_mix->RebinX(2);
  //q_IMnpipi_wSid_n_mix->RebinY(4);
  q_IMnpipi_wSid_n_mixpipinX->RebinX(2);
  q_IMnpipi_wSid_n_mixK0nX->RebinX(2);
  //q_IMnpipi_wSid_n_mix->GetXaxis()->SetRangeUser(1.1,1.3);
  //q_IMnpipi_wSid_n_mix->GetYaxis()->SetRangeUser(1.1,1.3);
  q_IMnpipi_wSid_n_data->SetContour(8);
  q_IMnpipi_wSid_n_mix->SetContour(8);
  q_IMnpipi_wSid_n_data->SetMinimum(1);
  q_IMnpipi_wSid_n_data->SetLineWidth(2);
  //q_IMnpipi_wSid_n_data->Draw("cont3");
  q_IMnpipi_wSid_n_data->Draw("col");
  //q_IMnpipi_wSid_n_data->Draw("col");
  q_IMnpipi_wSid_n_mix->SetLineColor(6);
  q_IMnpipi_wSid_n_mix->SetLineWidth(2);
  q_IMnpipi_wSid_n_mix->SetMinimum(1);
   
  
  cq_IMnpipi_wSid_n->cd(1);
  TH1D* IMnpipi_wSid_n_data = (TH1D*)q_IMnpipi_wSid_n_data->ProjectionX("IMnpipi_wSid_n_data");
  TH1D* IMnpipi_wSid_n_mix = (TH1D*)q_IMnpipi_wSid_n_mix->ProjectionX("IMnpipi_wSid_n_mix");
  TH1D* IMnpipi_wSid_n_mixpipinX = (TH1D*)q_IMnpipi_wSid_n_mixpipinX->ProjectionX("IMnpipi_wSid_n_mixpipinX");
  TH1D* IMnpipi_wSid_n_mixK0nX = (TH1D*)q_IMnpipi_wSid_n_mixK0nX->ProjectionX("IMnpipi_wSid_n_mixK0nX");
  IMnpipi_wSid_n_data->GetXaxis()->SetLabelSize(0);
  IMnpipi_wSid_n_data->SetMarkerStyle(20);
  IMnpipi_wSid_n_data->SetMarkerColor(1);
  IMnpipi_wSid_n_data->Draw("E");
  IMnpipi_wSid_n_mix->SetLineColor(6);
  IMnpipi_wSid_n_mix->SetMarkerStyle(20);
  IMnpipi_wSid_n_mix->SetMarkerColor(6);
  IMnpipi_wSid_n_mix->Draw("Esame");IMnpipi_wSid_n_mix->Print();
  IMnpipi_wSid_n_mixpipinX->SetLineColor(2);
  IMnpipi_wSid_n_mixpipinX->SetMarkerStyle(20);
  IMnpipi_wSid_n_mixpipinX->SetMarkerColor(2);IMnpipi_wSid_n_mixpipinX->Print();
  if(showBG)IMnpipi_wSid_n_mixpipinX->Draw("Esame");
  IMnpipi_wSid_n_mixK0nX->SetLineColor(3);
  IMnpipi_wSid_n_mixK0nX->SetMarkerStyle(20);
  IMnpipi_wSid_n_mixK0nX->SetMarkerColor(3);
  if(showBG)IMnpipi_wSid_n_mixK0nX->Draw("Esame");

  cq_IMnpipi_wSid_n->cd(4);
  q_wSid_n_data->SetTitle("");
  q_wSid_n_data->GetXaxis()->SetTitle("");
  
  gr_q_wSid_n_data->SetMarkerStyle(20);
  gr_q_wSid_n_data->GetYaxis()->SetLabelSize(0);
  gr_q_wSid_n_data->Draw("AP");
  gr_q_wSid_n_mix->SetMarkerStyle(20);
  gr_q_wSid_n_mix->SetMarkerColor(6);
  gr_q_wSid_n_mix->Draw("P");
  gr_q_wSid_n_mixpipinX->SetMarkerStyle(20);
  gr_q_wSid_n_mixpipinX->SetMarkerColor(2);
  if(showBG)gr_q_wSid_n_mixpipinX->Draw("P");
  gr_q_wSid_n_mixK0nX->SetMarkerStyle(20);
  gr_q_wSid_n_mixK0nX->SetMarkerColor(3);
  if(showBG)gr_q_wSid_n_mixK0nX->Draw("P");


  TCanvas *cq_IMnpipi_wSid_n_2 = new TCanvas("cq_IMnpipi_wSid_n_2","cq_IMnpipi_wSid_n_2",1000,800);
  cq_IMnpipi_wSid_n_2->Divide(2,1);
  cq_IMnpipi_wSid_n_2->cd(1);
  TH2D* q_IMnpipi_wSid_n_data_2 = (TH2D*)q_IMnpipi_wSid_n_data->Clone("q_IMnpipi_wSid_n_data_2");
  TH2D* q_IMnpipi_wSid_n_mix_2 = (TH2D*)q_IMnpipi_wSid_n_mix->Clone("q_IMnpipi_wSid_n_mix_2");
  q_IMnpipi_wSid_n_data_2->SetTitle("real data");
  q_IMnpipi_wSid_n_data_2->SetContour(10);
  q_IMnpipi_wSid_n_data_2->SetMinimum(1);
  q_IMnpipi_wSid_n_data_2->Draw("col");
  gPad->SetRightMargin(0);
  cq_IMnpipi_wSid_n_2->cd(2);
  q_IMnpipi_wSid_n_mix_2->SetTitle("MC");
  q_IMnpipi_wSid_n_mix_2->SetContour(10); 
  q_IMnpipi_wSid_n_mix_2->SetMinimum(1);
  q_IMnpipi_wSid_n_mix_2->SetMaximum(q_IMnpipi_wSid_n_data_2->GetMaximum());
  q_IMnpipi_wSid_n_mix_2->Draw("colz");
  gPad->SetLeftMargin(0);

  TCanvas *cq_IMnpipi_wSid_n_0 = new TCanvas("cq_IMnpipi_wSid_n_0","cq_IMnpipi_wSid_n_0",1000,1000);
  cq_IMnpipi_wSid_n_0->Divide(2,2,0,0);
  cq_IMnpipi_wSid_n_0->cd(3);
  TH2D* q_IMnpipi_wSid_n_data_0 = (TH2D*)q_IMnpipi_wSid_n_data->Clone("q_IMnpipi_wSid_n_data_0");
  TH2D* q_IMnpipi_wSid_n_mix_0 = (TH2D*)q_IMnpipi_wSid_n_mix->Clone("q_IMnpipi_wSid_n_mix_0");
  TH2D* q_IMnpipi_wSid_n_mixpipinX_0 = (TH2D*)q_IMnpipi_wSid_n_mixpipinX->Clone("q_IMnpipi_wSid_n_mixpipinX_0");
  TH2D* q_IMnpipi_wSid_n_mixK0nX_0 = (TH2D*)q_IMnpipi_wSid_n_mixK0nX->Clone("q_IMnpipi_wSid_n_mixK0nX_0");
  for(int ix=0;ix<q_IMnpipi_wSid_n_data_0->GetNbinsX();ix++){
    for(int iy=q_IMnpipi_wSid_n_data_0->GetYaxis()->FindBin(0.35);iy<q_IMnpipi_wSid_n_data_0->GetNbinsY();iy++){
      q_IMnpipi_wSid_n_data_0->SetBinContent(ix,iy,0);
      q_IMnpipi_wSid_n_data_0->SetBinError(ix,iy,0);
      q_IMnpipi_wSid_n_mix_0->SetBinContent(ix,iy,0);
      q_IMnpipi_wSid_n_mixpipinX_0->SetBinContent(ix,iy,0);
      q_IMnpipi_wSid_n_mixK0nX_0->SetBinContent(ix,iy,0);
    }
  }
  q_IMnpipi_wSid_n_data_0->Draw("col");
  
  cq_IMnpipi_wSid_n_0->cd(1);
  TH1D* IMnpipi_wSid_n_data_0 = (TH1D*)q_IMnpipi_wSid_n_data_0->ProjectionX("IMnpipi_wSid_n_data_0");
  TH1D* IMnpipi_wSid_n_mix_0 = (TH1D*)q_IMnpipi_wSid_n_mix_0->ProjectionX("IMnpipi_wSid_n_mix_0");
  TH1D* IMnpipi_wSid_n_mixpipinX_0 = (TH1D*)q_IMnpipi_wSid_n_mixpipinX_0->ProjectionX("IMnpipi_wSid_n_mixpipinX_0");
  TH1D* IMnpipi_wSid_n_mixK0nX_0 = (TH1D*)q_IMnpipi_wSid_n_mixK0nX_0->ProjectionX("IMnpipi_wSid_n_mixK0nX_0");
  IMnpipi_wSid_n_data_0->GetXaxis()->SetLabelSize(0);
  IMnpipi_wSid_n_data_0->SetMarkerStyle(20);
  IMnpipi_wSid_n_data_0->SetMarkerColor(1);
  IMnpipi_wSid_n_data_0->Draw("E");
  IMnpipi_wSid_n_mix_0->SetLineColor(6);
  IMnpipi_wSid_n_mix_0->SetMarkerStyle(20);
  IMnpipi_wSid_n_mix_0->SetMarkerColor(6);
  IMnpipi_wSid_n_mix_0->Draw("Esame");IMnpipi_wSid_n_mix->Print();
  IMnpipi_wSid_n_mixpipinX_0->SetLineColor(2);
  IMnpipi_wSid_n_mixpipinX_0->SetMarkerStyle(20);
  IMnpipi_wSid_n_mixpipinX_0->SetMarkerColor(2);IMnpipi_wSid_n_mixpipinX->Print();
  if(showBG)IMnpipi_wSid_n_mixpipinX_0->Draw("Esame");
  IMnpipi_wSid_n_mixK0nX_0->SetLineColor(3);
  IMnpipi_wSid_n_mixK0nX_0->SetMarkerStyle(20);
  IMnpipi_wSid_n_mixK0nX_0->SetMarkerColor(3);
  if(showBG)IMnpipi_wSid_n_mixK0nX_0->Draw("Esame");
  
  cq_IMnpipi_wSid_n_0->cd(4);
  TH1D* q_wSid_n_data_0 = (TH1D*)q_IMnpipi_wSid_n_data_0->ProjectionY("q_wSid_n_data_0");
  TH1D* q_wSid_n_mix_0 = (TH1D*)q_IMnpipi_wSid_n_mix_0->ProjectionY("q_wSid_n_mix_0");
  TH1D* q_wSid_n_mixpipinX_0 = (TH1D*)q_IMnpipi_wSid_n_mixpipinX_0->ProjectionY("q_wSid_n_mixpipinX_0");
  TH1D* q_wSid_n_mixK0nX_0 = (TH1D*)q_IMnpipi_wSid_n_mixK0nX_0->ProjectionY("q_wSid_n_mixK0nX_0");
  TGraphErrors *gr_q_wSid_n_data_0 = new TGraphErrors();
  TGraphErrors *gr_q_wSid_n_mix_0 = new TGraphErrors();
  TGraphErrors *gr_q_wSid_n_mixpipinX_0 = new TGraphErrors();
  TGraphErrors *gr_q_wSid_n_mixK0nX_0 = new TGraphErrors();
  
  for(int ibin=0;ibin<q_wSid_n_data_0->GetNbinsX();ibin++){
    double cont = q_wSid_n_data_0->GetBinContent(ibin);
    double err = q_wSid_n_data_0->GetBinError(ibin);
    double bincenter = q_wSid_n_data_0->GetBinCenter(ibin);
    gr_q_wSid_n_data_0->SetPoint(ibin,cont, bincenter);
    gr_q_wSid_n_data_0->SetPointError(ibin,err,0);
  }
  for(int ibin=0;ibin<q_wSid_n_mix_0->GetNbinsX();ibin++){
    double cont = q_wSid_n_mix_0->GetBinContent(ibin);
    double err = q_wSid_n_mix_0->GetBinError(ibin);
    double bincenter = q_wSid_n_mix_0->GetBinCenter(ibin);
    gr_q_wSid_n_mix_0->SetPoint(ibin,cont, bincenter);
    gr_q_wSid_n_mix_0->SetPointError(ibin,err,0);
  }
  for(int ibin=0;ibin<q_wSid_n_mixpipinX_0->GetNbinsX();ibin++){
    double cont = q_wSid_n_mixpipinX_0->GetBinContent(ibin);
    double err = q_wSid_n_mixpipinX_0->GetBinError(ibin);
    double bincenter = q_wSid_n_mixpipinX_0->GetBinCenter(ibin);
    gr_q_wSid_n_mixpipinX_0->SetPoint(ibin,cont, bincenter);
    gr_q_wSid_n_mixpipinX_0->SetPointError(ibin,err,0);
  }
  for(int ibin=0;ibin<q_wSid_n_mixK0nX_0->GetNbinsX();ibin++){
    double cont = q_wSid_n_mixK0nX_0->GetBinContent(ibin);
    double err = q_wSid_n_mixK0nX_0->GetBinError(ibin);
    double bincenter = q_wSid_n_mixK0nX_0->GetBinCenter(ibin);
    gr_q_wSid_n_mixK0nX_0->SetPoint(ibin,cont, bincenter);
    gr_q_wSid_n_mixK0nX_0->SetPointError(ibin,err,0);
  }
  
  q_wSid_n_data_0->SetTitle("");
  q_wSid_n_data_0->GetXaxis()->SetTitle("");
  gr_q_wSid_n_data_0->SetMarkerStyle(20);
  gr_q_wSid_n_data_0->GetYaxis()->SetLabelSize(0);
  gr_q_wSid_n_data_0->Draw("AP");
  gr_q_wSid_n_mix_0->SetMarkerStyle(20);
  gr_q_wSid_n_mix_0->SetMarkerColor(6);
  gr_q_wSid_n_mix_0->Draw("P");
  gr_q_wSid_n_mixpipinX_0->SetMarkerStyle(20);
  gr_q_wSid_n_mixpipinX_0->SetMarkerColor(2);
  if(showBG)gr_q_wSid_n_mixpipinX_0->Draw("P");
  gr_q_wSid_n_mixK0nX_0->SetMarkerStyle(20);
  gr_q_wSid_n_mixK0nX_0->SetMarkerColor(3);
  if(showBG)gr_q_wSid_n_mixK0nX_0->Draw("P");
  q_wSid_n_data_0->SetTitle("");
  q_wSid_n_data_0->GetXaxis()->SetTitle("");
  
  gr_q_wSid_n_data_0->SetMarkerStyle(20);
  gr_q_wSid_n_data_0->GetYaxis()->SetLabelSize(0);
  gr_q_wSid_n_data_0->GetYaxis()->SetRangeUser(0,1.5);
  gr_q_wSid_n_data_0->Draw("AP");
  gr_q_wSid_n_mix_0->SetMarkerStyle(20);
  gr_q_wSid_n_mix_0->SetMarkerColor(6);
  gr_q_wSid_n_mix_0->Draw("P");
  gr_q_wSid_n_mixpipinX_0->SetMarkerStyle(20);
  gr_q_wSid_n_mixpipinX_0->SetMarkerColor(2);
  if(showBG)gr_q_wSid_n_mixpipinX_0->Draw("P");
  gr_q_wSid_n_mixK0nX_0->SetMarkerStyle(20);
  gr_q_wSid_n_mixK0nX_0->SetMarkerColor(3);
  if(showBG)gr_q_wSid_n_mixK0nX_0->Draw("P");


  TCanvas *cq_IMnpipi_wSid_n_350 = new TCanvas("cq_IMnpipi_wSid_n_350","cq_IMnpipi_wSid_n_350",1000,1000);
  cq_IMnpipi_wSid_n_350->Divide(2,2,0,0);
  cq_IMnpipi_wSid_n_350->cd(3);
  TH2D* q_IMnpipi_wSid_n_data_350 = (TH2D*)q_IMnpipi_wSid_n_data->Clone("q_IMnpipi_wSid_n_data_350");
  TH2D* q_IMnpipi_wSid_n_mix_350 = (TH2D*)q_IMnpipi_wSid_n_mix->Clone("q_IMnpipi_wSid_n_mix_350");
  TH2D* q_IMnpipi_wSid_n_mixpipinX_350 = (TH2D*)q_IMnpipi_wSid_n_mixpipinX->Clone("q_IMnpipi_wSid_n_mixpipinX_350");
  TH2D* q_IMnpipi_wSid_n_mixK0nX_350 = (TH2D*)q_IMnpipi_wSid_n_mixK0nX->Clone("q_IMnpipi_wSid_n_mixK0nX_350");
  for(int ix=0;ix<q_IMnpipi_wSid_n_data_350->GetNbinsX();ix++){
    for(int iy=0;iy<q_IMnpipi_wSid_n_data_350->GetYaxis()->FindBin(0.35);iy++){
      q_IMnpipi_wSid_n_data_350->SetBinContent(ix,iy,0);
      q_IMnpipi_wSid_n_data_350->SetBinError(ix,iy,0);
      q_IMnpipi_wSid_n_mix_350->SetBinContent(ix,iy,0);
      q_IMnpipi_wSid_n_mixpipinX_350->SetBinContent(ix,iy,0);
      q_IMnpipi_wSid_n_mixK0nX_350->SetBinContent(ix,iy,0);
    }
  }
  q_IMnpipi_wSid_n_data_350->Draw("col");
  
  cq_IMnpipi_wSid_n_350->cd(1);
  TH1D* IMnpipi_wSid_n_data_350 = (TH1D*)q_IMnpipi_wSid_n_data_350->ProjectionX("IMnpipi_wSid_n_data_350");
  TH1D* IMnpipi_wSid_n_mix_350 = (TH1D*)q_IMnpipi_wSid_n_mix_350->ProjectionX("IMnpipi_wSid_n_mix_350");
  TH1D* IMnpipi_wSid_n_mixpipinX_350 = (TH1D*)q_IMnpipi_wSid_n_mixpipinX_350->ProjectionX("IMnpipi_wSid_n_mixpipinX_350");
  TH1D* IMnpipi_wSid_n_mixK0nX_350 = (TH1D*)q_IMnpipi_wSid_n_mixK0nX_350->ProjectionX("IMnpipi_wSid_n_mixK0nX_350");
  IMnpipi_wSid_n_data_350->GetXaxis()->SetLabelSize(0);
  IMnpipi_wSid_n_data_350->SetMarkerStyle(20);
  IMnpipi_wSid_n_data_350->SetMarkerColor(1);
  IMnpipi_wSid_n_data_350->Draw("E");
  IMnpipi_wSid_n_mix_350->SetLineColor(6);
  IMnpipi_wSid_n_mix_350->SetMarkerStyle(20);
  IMnpipi_wSid_n_mix_350->SetMarkerColor(6);
  IMnpipi_wSid_n_mix_350->Draw("Esame");IMnpipi_wSid_n_mix->Print();
  IMnpipi_wSid_n_mixpipinX_350->SetLineColor(2);
  IMnpipi_wSid_n_mixpipinX_350->SetMarkerStyle(20);
  IMnpipi_wSid_n_mixpipinX_350->SetMarkerColor(2);IMnpipi_wSid_n_mixpipinX->Print();
  if(showBG)IMnpipi_wSid_n_mixpipinX_350->Draw("Esame");
  IMnpipi_wSid_n_mixK0nX_350->SetLineColor(3);
  IMnpipi_wSid_n_mixK0nX_350->SetMarkerStyle(20);
  IMnpipi_wSid_n_mixK0nX_350->SetMarkerColor(3);
  if(showBG)IMnpipi_wSid_n_mixK0nX_350->Draw("Esame");
  
  cq_IMnpipi_wSid_n_350->cd(4);
  TH1D* q_wSid_n_data_350 = (TH1D*)q_IMnpipi_wSid_n_data_350->ProjectionY("q_wSid_n_data_350");
  TH1D* q_wSid_n_mix_350 = (TH1D*)q_IMnpipi_wSid_n_mix_350->ProjectionY("q_wSid_n_mix_350");
  TH1D* q_wSid_n_mixpipinX_350 = (TH1D*)q_IMnpipi_wSid_n_mixpipinX_350->ProjectionY("q_wSid_n_mixpipinX_350");
  TH1D* q_wSid_n_mixK0nX_350 = (TH1D*)q_IMnpipi_wSid_n_mixK0nX_350->ProjectionY("q_wSid_n_mixK0nX_350");
  TGraphErrors *gr_q_wSid_n_data_350 = new TGraphErrors();
  TGraphErrors *gr_q_wSid_n_mix_350 = new TGraphErrors();
  TGraphErrors *gr_q_wSid_n_mixpipinX_350 = new TGraphErrors();
  TGraphErrors *gr_q_wSid_n_mixK0nX_350 = new TGraphErrors();
  
  for(int ibin=0;ibin<q_wSid_n_data_350->GetNbinsX();ibin++){
    double cont = q_wSid_n_data_350->GetBinContent(ibin);
    double err = q_wSid_n_data_350->GetBinError(ibin);
    double bincenter = q_wSid_n_data_350->GetBinCenter(ibin);
    gr_q_wSid_n_data_350->SetPoint(ibin,cont, bincenter);
    gr_q_wSid_n_data_350->SetPointError(ibin,err,0);
  }
  for(int ibin=0;ibin<q_wSid_n_mix_350->GetNbinsX();ibin++){
    double cont = q_wSid_n_mix_350->GetBinContent(ibin);
    double err = q_wSid_n_mix_350->GetBinError(ibin);
    double bincenter = q_wSid_n_mix_350->GetBinCenter(ibin);
    gr_q_wSid_n_mix_350->SetPoint(ibin,cont, bincenter);
    gr_q_wSid_n_mix_350->SetPointError(ibin,err,0);
  }
  for(int ibin=0;ibin<q_wSid_n_mixpipinX_350->GetNbinsX();ibin++){
    double cont = q_wSid_n_mixpipinX_350->GetBinContent(ibin);
    double err = q_wSid_n_mixpipinX_350->GetBinError(ibin);
    double bincenter = q_wSid_n_mixpipinX_350->GetBinCenter(ibin);
    gr_q_wSid_n_mixpipinX_350->SetPoint(ibin,cont, bincenter);
    gr_q_wSid_n_mixpipinX_350->SetPointError(ibin,err,0);
  }
  for(int ibin=0;ibin<q_wSid_n_mixK0nX_350->GetNbinsX();ibin++){
    double cont = q_wSid_n_mixK0nX_350->GetBinContent(ibin);
    double err = q_wSid_n_mixK0nX_350->GetBinError(ibin);
    double bincenter = q_wSid_n_mixK0nX_350->GetBinCenter(ibin);
    gr_q_wSid_n_mixK0nX_350->SetPoint(ibin,cont, bincenter);
    gr_q_wSid_n_mixK0nX_350->SetPointError(ibin,err,0);
  }
  
  q_wSid_n_data_350->SetTitle("");
  q_wSid_n_data_350->GetXaxis()->SetTitle("");
  gr_q_wSid_n_data_350->SetMarkerStyle(20);
  gr_q_wSid_n_data_350->GetYaxis()->SetLabelSize(0);
  gr_q_wSid_n_data_350->Draw("AP");
  gr_q_wSid_n_mix_350->SetMarkerStyle(20);
  gr_q_wSid_n_mix_350->SetMarkerColor(6);
  gr_q_wSid_n_mix_350->Draw("P");
  gr_q_wSid_n_mixpipinX_350->SetMarkerStyle(20);
  gr_q_wSid_n_mixpipinX_350->SetMarkerColor(2);
  if(showBG)gr_q_wSid_n_mixpipinX_350->Draw("P");
  gr_q_wSid_n_mixK0nX_350->SetMarkerStyle(20);
  gr_q_wSid_n_mixK0nX_350->SetMarkerColor(3);
  if(showBG)gr_q_wSid_n_mixK0nX_350->Draw("P");
  q_wSid_n_data_350->SetTitle("");
  q_wSid_n_data_350->GetXaxis()->SetTitle("");
  
  gr_q_wSid_n_data_350->SetMarkerStyle(20);
  gr_q_wSid_n_data_350->GetYaxis()->SetLabelSize(0);
  gr_q_wSid_n_data_350->GetYaxis()->SetRangeUser(0,1.5);
  gr_q_wSid_n_data_350->Draw("AP");
  gr_q_wSid_n_mix_350->SetMarkerStyle(20);
  gr_q_wSid_n_mix_350->SetMarkerColor(6);
  gr_q_wSid_n_mix_350->Draw("P");
  gr_q_wSid_n_mixpipinX_350->SetMarkerStyle(20);
  gr_q_wSid_n_mixpipinX_350->SetMarkerColor(2);
  if(showBG)gr_q_wSid_n_mixpipinX_350->Draw("P");
  gr_q_wSid_n_mixK0nX_350->SetMarkerStyle(20);
  gr_q_wSid_n_mixK0nX_350->SetMarkerColor(3);
  if(showBG)gr_q_wSid_n_mixK0nX_350->Draw("P");

  //subtracted 
  TCanvas *cq_IMnpipi_wSid_n_sub = new TCanvas("cq_IMnpipi_wSid_n_sub","cq_IMnpipi_wSid_n_sub",1000,1000);
  TH2D* q_IMnpipi_wSid_n_sub = (TH2D*)q_IMnpipi_wSid_n_data->Clone("q_IMnpipi_wSid_n_data");
  TH2D* q_IMnpipi_wSid_n_mix_no0 = (TH2D*)q_IMnpipi_wSid_n_mix->Clone("q_IMnpipi_wSid_n_mix_no0");
  //clean up 0 count bin
  for(int ix=0;ix<q_IMnpipi_wSid_n_sub->GetNbinsX();ix++){
    for(int iy=0;iy<q_IMnpipi_wSid_n_sub->GetNbinsY();iy++){
      double cont = q_IMnpipi_wSid_n_sub->GetBinContent(ix,iy);
      if(cont <1.0) q_IMnpipi_wSid_n_mix_no0->SetBinContent(ix,iy,0);
    }
  }

  q_IMnpipi_wSid_n_sub->Add(q_IMnpipi_wSid_n_mix_no0,-1.0);
  //q_IMnpipi_wSid_n_sub->Divide(q_IMnpipi_wSid_n_data);
  //q_IMnpipi_wSid_n_sub->SetTitle("(data-MC)/data");
  //q_IMnpipi_wSid_n_sub->RebinX(2);
  //q_IMnpipi_wSid_n_sub->RebinY(2);
  //q_IMnpipi_wSid_n_sub->Scale(0.25);
  //q_IMnpipi_wSid_n_sub->GetZaxis()->SetRangeUser(-0.3,0.3);
  //q_IMnpipi_wSid_n_sub->SetMinimum(0);
  q_IMnpipi_wSid_n_sub->Draw("colz");
   
  TCanvas *cq_IMnpipi_wSid_n_sub_proj0 = new TCanvas("cq_IMnpipi_wSid_n_sub_proj0","cq_IMnpipi_wSid_n_sub_proj0");
  TH1D* IMnpipi_wSid_n_sub_0 = (TH1D*) q_IMnpipi_wSid_n_sub->ProjectionX("IMnpipi_wSid_n_sub_0",0,q_IMnpipi_wSid_n_sub->GetYaxis()->FindBin(0.35)-1);
  IMnpipi_wSid_n_sub_0->SetMarkerStyle(20);
  IMnpipi_wSid_n_sub_0->Draw();

  TCanvas *cq_IMnpipi_wSid_n_sub_proj350 = new TCanvas("cq_IMnpipi_wSid_n_sub_proj350_600","cq_IMnpipi_wSid_n_sub_proj350_600");
  TH1D* IMnpipi_wSid_n_sub_350_600 = (TH1D*) q_IMnpipi_wSid_n_sub->ProjectionX("IMnpipi_wSid_n_sub_350_600",q_IMnpipi_wSid_n_sub->GetYaxis()->FindBin(0.35),q_IMnpipi_wSid_n_sub->GetYaxis()->FindBin(0.60));
  IMnpipi_wSid_n_sub_350_600->SetMarkerStyle(20);
  IMnpipi_wSid_n_sub_350_600->Draw();
  TF1 *f_L1520 = new TF1("f_L1520","gaus(0)+gaus(3)",1,2);
  f_L1520->SetParameter(0,180);
  f_L1520->SetParameter(1,1.517);
  f_L1520->SetParameter(2,0.0156);
  f_L1520->SetParameter(3,90);
  f_L1520->SetParameter(4,1.405);
  f_L1520->SetParameter(5,0.050);
  f_L1520->SetLineColor(2);
  f_L1520->SetNpx(1000);
  //f_L1520->Draw("same");

  TF1* f_bf = new TF1("f_bf","4.1*TMath::BreitWigner(x, 1.5195, 0.0156)+7.2*TMath::BreitWigner(x, 1.4051, 0.050)",1,2);
  f_bf->SetLineColor(3);
  f_bf->SetNpx(5000);
  //f_bf->Draw("same");


  //
  //signal check q vs IMnpipi wK0 selection BG subtracted
  //
  TH2D* q_IMnpipi_wK0_wSid_n_data_sub = (TH2D*)q_IMnpipi_wK0_wSid_n_data->Clone("q_IMnpipi_wK0_wSid_n_data_sub");
  q_IMnpipi_wK0_wSid_n_data_sub->Add(q_IMnpipi_wK0_wSid_n_mix,-1.0);
  q_IMnpipi_wK0_wSid_n_data_sub->Add(q_IMnpipi_wK0_wSid_n_mixgeta,-1.0);
  TCanvas *cq_IMnpipi_wK0_wSid_n_data_sub = new TCanvas("cq_IMnpipi_wK0_wSid_n_data_sub","cq_IMnpipi_wK0_wSid_n_data_sub",1000,1000);
  cq_IMnpipi_wK0_wSid_n_data_sub->Divide(2,2,0,0);
  cq_IMnpipi_wK0_wSid_n_data_sub->cd(3);
  q_IMnpipi_wK0_wSid_n_data_sub->RebinX(2);
  q_IMnpipi_wK0_wSid_n_data_sub->SetTitle("");
  q_IMnpipi_wK0_wSid_n_data_sub->Draw("colz");
  cq_IMnpipi_wK0_wSid_n_data_sub->cd(1);
  TH1D* IMnpipi_wK0_wSid_n_data_sub = (TH1D*)q_IMnpipi_wK0_wSid_n_data_sub->ProjectionX("IMnpipi_wK0_wSid_n_data_sub");
  IMnpipi_wK0_wSid_n_data_sub->SetMarkerStyle(20);
  IMnpipi_wK0_wSid_n_data_sub->Draw("");
  cq_IMnpipi_wK0_wSid_n_data_sub->cd(4);
  TH1D* q_wK0_wSid_n_data_sub = (TH1D*)q_IMnpipi_wK0_wSid_n_data_sub->ProjectionY("q_wK0_wSid_n_data_sub");
  TGraphErrors *gr_q_wK0_wSid_n_data_sub = new TGraphErrors();
  for(int ibin=0;ibin<q_wK0_wSid_n_data_sub->GetNbinsX();ibin++){
    double cont = q_wK0_wSid_n_data_sub->GetBinContent(ibin);
    double err = q_wK0_wSid_n_data_sub->GetBinError(ibin);
    double bincenter = q_wK0_wSid_n_data_sub->GetBinCenter(ibin);
    gr_q_wK0_wSid_n_data_sub->SetPoint(ibin,cont, bincenter);
    gr_q_wK0_wSid_n_data_sub->SetPointError(ibin,err,0);
  }
  gr_q_wK0_wSid_n_data_sub->GetYaxis()->SetRangeUser(0,1.5);
  gr_q_wK0_wSid_n_data_sub->SetMarkerStyle(20);
  gr_q_wK0_wSid_n_data_sub->Draw("AP");
  
  TCanvas *cq_BG_comp = new TCanvas("cq_BG_comp","cq_BG_comp",1600,800);
  cq_BG_comp->Divide(3,1);
  //q_IMnpipi_wSid_n_mix->Draw("colz");
  //gPad->SetRightMargin(0.1);
  //gPad->SetLeftMargin(0);
  //TFile *fGSp = TFile::Open("../simIMpisigma_nSppim_pippimn_v108_out.root","READ");
  //TFile *fGSp = TFile::Open("../simIMpisigma_nSppim_pippimn_v113_out.root","READ");
  //TFile *fGSp_old = TFile::Open("../simIMpisigma_nSppim_pippimn_v117_out_v1.root","READ");
  TFile *fGSp = TFile::Open("../simIMpisigma_nSppim_pippimn_v124_out.root","READ");
  TH2D* q_IMnpipi_wSid_n = (TH2D*)fGSp->Get("q_IMnpipi_wSid_n");
  TH2D* q_IMnpipi_wSid_n_true = (TH2D*)q_IMnpipi_wSid_n->Clone("q_IMnpipi_wSid_n_true");
  TH2D* q_IMnpipi_wSid_n_fake = (TH2D*)fGSp->Get("q_IMnpipi_wSid_n_fake");
  //TH2D* q_IMnpipi_wSid_n_fake_v1 = (TH2D*)fGSp_old->Get("q_IMnpipi_wSid_n_fake");
  q_IMnpipi_wSid_n_true->Add(q_IMnpipi_wSid_n_fake,-1.0);
  q_IMnpipi_wSid_n_true->SetTitle("GEANT true");
  cq_BG_comp->cd(1);
  q_IMnpipi_wSid_n_true->RebinX(2);
  q_IMnpipi_wSid_n_true->Draw("colz");
  cq_BG_comp->cd(2);
  q_IMnpipi_wSid_n_fake->RebinX(2);
  //q_IMnpipi_wSid_n_fake_v1->RebinX(2);
  q_IMnpipi_wSid_n_fake->SetMaximum(q_IMnpipi_wSid_n_true->GetMaximum());
  q_IMnpipi_wSid_n_fake->SetTitle("GEANT fake");
  q_IMnpipi_wSid_n_fake->Draw("colz");
  cq_BG_comp->cd(3);
  q_IMnpipi_wSid_n_mix->SetMaximum(q_IMnpipi_wSid_n_true->GetMaximum());
  q_IMnpipi_wSid_n_mix->SetTitle("#pi^{+}#pi^{-}n_{FAKE}X");
  q_IMnpipi_wSid_n_mix->Draw("colz");

  std::cout << "model BG" << q_IMnpipi_wSid_n_mix->Integral() << std::endl;
  std::cout << "GEANT BG" << q_IMnpipi_wSid_n_fake->Integral() << std::endl;

  TCanvas *cq_BG_comp_X = new TCanvas("cq_BG_comp_X","cq_BG_comp_X");
  TH1D* IMnpipi_true = (TH1D*)q_IMnpipi_wSid_n_true->ProjectionX("IMnpipi_true");
  IMnpipi_true->Draw("HE");
  TH1D* IMnpipi_BGmodel = (TH1D*) q_IMnpipi_wSid_n_mix->ProjectionX("IMnpipi_BGmodel");
  IMnpipi_BGmodel->SetLineColor(3);
  IMnpipi_BGmodel->Draw("HEsame");
  TH1D* IMnpipi_BGGEANT = (TH1D*) q_IMnpipi_wSid_n_fake->ProjectionX("IMnpipi_BGGEANT");
  IMnpipi_BGGEANT->SetLineColor(2);
  IMnpipi_BGGEANT->Draw("HEsame");
  //TH1D* IMnpipi_BGGEANT_v1 = (TH1D*) q_IMnpipi_wSid_n_fake_v1->ProjectionX("IMnpipi_BGGEANT_v1");
  //IMnpipi_BGGEANT_v1->SetLineColor(46);
  //IMnpipi_BGGEANT_v1->Draw("HEsame");
  
  const int bin350 = q_IMnpipi_wSid_n_true->GetYaxis()->FindBin(0.35);
  const int bin800 = q_IMnpipi_wSid_n_true->GetYaxis()->FindBin(0.80);

  TCanvas *cq_BG_comp_X_0 = new TCanvas("cq_BG_comp_X_0","cq_BG_comp_X_0");
  TH1D* IMnpipi_true_0 = (TH1D*)q_IMnpipi_wSid_n_true->ProjectionX("IMnpipi_true_0",0,bin350-1);
  IMnpipi_true_0->Draw("HE");
  TH1D* IMnpipi_BGmodel_0 = (TH1D*) q_IMnpipi_wSid_n_mix->ProjectionX("IMnpipi_BGmodel_0",0,bin350-1);
  IMnpipi_BGmodel_0->SetLineColor(3);
  IMnpipi_BGmodel_0->Draw("HEsame");
  TH1D* IMnpipi_BGGEANT_0 = (TH1D*) q_IMnpipi_wSid_n_fake->ProjectionX("IMnpipi_BGGEANT_0",0,bin350-1);
  IMnpipi_BGGEANT_0->SetLineColor(2);
  IMnpipi_BGGEANT_0->Draw("HEsame");
  //TH1D* IMnpipi_BGGEANT_0_v1 = (TH1D*) q_IMnpipi_wSid_n_fake_v1->ProjectionX("IMnpipi_BGGEANT_0_v1",0,bin350-1);
  //IMnpipi_BGGEANT_0_v1->SetLineColor(46);
  //IMnpipi_BGGEANT_0_v1->Draw("HEsame");

  TCanvas *cq_BG_comp_X_350 = new TCanvas("cq_BG_comp_X_350","cq_BG_comp_X_350");
  TH1D* IMnpipi_true_350 = (TH1D*)q_IMnpipi_wSid_n_true->ProjectionX("IMnpipi_true_350",bin350,bin800);
  IMnpipi_true_350->Draw("HE");
  TH1D* IMnpipi_BGmodel_350 = (TH1D*) q_IMnpipi_wSid_n_mix->ProjectionX("IMnpipi_BGmodel_350",bin350,bin800);
  IMnpipi_BGmodel_350->SetLineColor(3);
  IMnpipi_BGmodel_350->Draw("HEsame");
  TH1D* IMnpipi_BGGEANT_350 = (TH1D*) q_IMnpipi_wSid_n_fake->ProjectionX("IMnpipi_BGGEANT_350",bin350,bin800);
  IMnpipi_BGGEANT_350->SetLineColor(2);
  IMnpipi_BGGEANT_350->Draw("HEsame");
  //TH1D* IMnpipi_BGGEANT_350_v1 = (TH1D*) q_IMnpipi_wSid_n_fake_v1->ProjectionX("IMnpipi_BGGEANT_350_v1",bin350,bin800);
  //IMnpipi_BGGEANT_350_v1->SetLineColor(46);
  //IMnpipi_BGGEANT_350_v1->Draw("HEsame");

  //other signal checks
  //MMnmiss 
  TCanvas *cMMnmiss_BG_comp = new TCanvas("cMMnmiss_BG_comp","cMMnmiss_BG_comp");
  TH2D* MMnmiss_IMpippim_dE_wSid_reco = (TH2D*)fGSp->Get("MMnmiss_IMpippim_dE_wSid");
  TH2D* MMnmiss_IMpippim_dE_wSid_fake = (TH2D*)fGSp->Get("MMnmiss_IMpippim_dE_wSid_fake");
  //MMnmiss_IMpippim_dE_wSid_true->Add(MMnmiss_IMpippim_dE_wSid_fake,-1);
  MMnmiss_IMpippim_dE_wSid_reco->RebinX(4);
  MMnmiss_IMpippim_dE_wSid_fake->RebinX(4);
  TH1D* MMnmiss_reco = (TH1D*)MMnmiss_IMpippim_dE_wSid_reco->ProjectionY("MMnmiss_reco");
  TH1D* MMnmiss_fake = (TH1D*)MMnmiss_IMpippim_dE_wSid_fake->ProjectionY("MMnmiss_fake");
  MMnmiss_reco->Draw("HE");
  MMnmiss_fake->SetLineColor(2);
  MMnmiss_fake->Draw("HEsame");
  MMnmiss_wSid_mix->SetLineColor(3);
  MMnmiss_wSid_mix->Draw("HEsame");

  //IMnpip 
  TCanvas *cIMnpip_BG_comp = new TCanvas("cIMnpip_BG_comp","cIMnpip_BG_comp");
  TH2D* IMnpim_IMnpip_dE_n_reco = (TH2D*)fGSp->Get("IMnpim_IMnpip_dE_n");
  TH2D* IMnpim_IMnpip_dE_n_fake = (TH2D*)fGSp->Get("IMnpim_IMnpip_dE_n_fake");
  //MMnmiss_IMpippim_dE_wSid_true->Add(MMnmiss_IMpippim_dE_wSid_fake,-1);
  //MMnmiss_IMpippim_dE_wSid_true->RebinX(4);
  //MMnmiss_IMpippim_dE_wSid_fake->RebinX(4);
  TH1D* IMnpip_reco = (TH1D*)IMnpim_IMnpip_dE_n_reco->ProjectionX("IMnpip_reco");
  TH1D* IMnpip_fake = (TH1D*)IMnpim_IMnpip_dE_n_fake->ProjectionX("IMnpip_fake");
  IMnpip_reco->Draw("HE");
  IMnpip_fake->SetLineColor(2);
  IMnpip_fake->Draw("HEsame");
  IMnpip_n_mix->SetLineColor(3);
  IMnpip_n_mix->SetMarkerColor(3);
  IMnpip_n_mix->Draw("HEsame");
  
  TCanvas *cMMnmiss_sub = new TCanvas("cMMnmiss_sub","cMMnmiss_sub");
  cMMnmiss_sub->cd();
  TH1D* MMnmiss_sub = (TH1D*)MMnmiss_reco->Clone("MMnmiss_sub");
  MMnmiss_sub->Add(MMnmiss_fake,-1.0);
  MMnmiss_sub->Draw("HE");
  TF1 *fgaus_MM = new TF1("fgaus_MM","gaus",0.5,1.5);
  fgaus_MM->SetParameter(0,8000);
  fgaus_MM->SetParameter(1,anacuts::neutron_center);
  fgaus_MM->SetParameter(2,anacuts::neutron_sigma);
  fgaus_MM->SetLineColor(3);
  fgaus_MM->SetNpx(10000);
  fgaus_MM->Draw("same");

  TCanvas *cIMnip_sub = new TCanvas("cIMnpip_sub","cIMnpip_sub");
  cIMnpip_sub->cd();
  TH1D* IMnpip_sub = (TH1D*)IMnpip_reco->Clone("IMnpip_sub");
  IMnpip_sub->Add(IMnpip_fake,-1.0);
  IMnpip_sub->Draw("HE");
  TF1 *fgaus_S = new TF1("fgaus_S","gaus",0.5,1.5);
  fgaus_S->SetParameter(0,11000);
  fgaus_S->SetParameter(1,anacuts::Sigmap_center);
  fgaus_S->SetParameter(2,anacuts::Sigmap_sigma);
  fgaus_S->SetLineColor(3);
  fgaus_S->SetNpx(10000);
  fgaus_S->Draw("same");


}
