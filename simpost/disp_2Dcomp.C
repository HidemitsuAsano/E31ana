#include "../post/weightfunc.h"

void disp_2Dcomp(const char *filename="comp_fakedata_out.root")
{
 
  TFile *f = TFile::Open(filename,"READ");
  f->cd();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  //const char opt[10]="cont3";
  const char opt[10]="cont1";
  //const char opt[10]="colz";
  gStyle->SetErrorX(0.);  
  TCanvas *cIMnpim_IMnpip_n = new TCanvas("cIMnpim_IMnpip_n", "cIMnpim_IMnpip_n",1200,1000);
  TPad *center_pad = new TPad("center_pad", "center_pad",0.0,0.0,0.6,0.6);
  center_pad->Draw();

  TPad *right_pad = new TPad("right_pad", "right_pad",0.55,0.0,1.0,0.6);
  right_pad->Draw();

  TPad *top_pad = new TPad("top_pad", "top_pad",0.0,0.55,0.6,1.0);
  top_pad->Draw();

  TH2D* IMnpim_IMnpip_n_rdata = f->Get("IMnpim_IMnpip_n_rdata");
  TH2D* IMnpim_IMnpip_n_mc    = f->Get("IMnpim_IMnpip_n_mc");
  
  center_pad->cd();
  //center_pad->SetTopMargin(0.0);
  //center_pad->SetRightMargin(0.1);
  IMnpim_IMnpip_n_rdata->RebinX(4);
  IMnpim_IMnpip_n_rdata->RebinY(4);
  IMnpim_IMnpip_n_rdata->GetXaxis()->SetRangeUser(1,1.7);
  IMnpim_IMnpip_n_rdata->GetYaxis()->SetRangeUser(1,1.7);
  IMnpim_IMnpip_n_rdata->Draw(opt);
  
  gPad->Update();
  //TPaletteAxis *palette = (TPaletteAxis*)IMnpim_IMnpip_n_rdata->GetListOfFunctions()->FindObject("palette");
  // the following lines moe the paletter. Choose the values you need for the position.
  //palette->SetX1NDC(0.78);
  //palette->SetX2NDC(0.83);
  //palette->SetY1NDC(0.1);
  //palette->SetY2NDC(0.90);
  gPad->Modified();
  gPad->Update();
   
  top_pad->cd();
  TH1D* IMnpip_n_rdata = (TH1D*)IMnpim_IMnpip_n_rdata->ProjectionX("IMnpip_n_rdata");
  IMnpip_n_rdata->SetTitle("");
  IMnpip_n_rdata->GetXaxis()->SetTitle("");
  IMnpip_n_rdata->Draw("E1");
  top_pad->SetBottomMargin(0);
  IMnpip_n_rdata->GetXaxis()->SetLabelSize(0);
  
  right_pad->cd();
  TH1D* IMnpim_n_rdata = (TH1D*)IMnpim_IMnpip_n_rdata->ProjectionY("IMnpim_n_rdata");
  IMnpim_n_rdata->SetTitle("");
  IMnpim_n_rdata->GetXaxis()->SetTitle("");
  TGraphErrors *gr_IMnpim_n_rdata = new TGraphErrors();
  for(int ibin=0;ibin<IMnpim_n_rdata->GetNbinsX();ibin++){
    double cont = IMnpim_n_rdata->GetBinContent(ibin);
    double bincenter = IMnpim_n_rdata->GetBinCenter(ibin);
    double err = IMnpim_n_rdata->GetBinError(ibin);
    gr_IMnpim_n_rdata->SetPoint(ibin,cont,bincenter);
    gr_IMnpim_n_rdata->SetPointError(ibin,err,0);
  }
  gr_IMnpim_n_rdata->GetYaxis()->SetLabelOffset(999);
  gr_IMnpim_n_rdata->GetYaxis()->SetLabelSize(0);
  right_pad->SetLeftMargin(0);
  //IMnpim_n_rdata->Draw("HE");
  gr_IMnpim_n_rdata->GetYaxis()->SetRangeUser(IMnpim_IMnpip_n_rdata->GetYaxis()->GetXmin(),1.7);//IMnpim_IMnpip_n_rdata->GetYaxis()->GetXmax());
  gr_IMnpim_n_rdata->SetLineWidth(2);
  gr_IMnpim_n_rdata->Draw("AP");
  



  //draw BG at first
  TH2D* IMnpim_IMnpip_woK0_woSid_won_data = (TH2D*)f->Get("IMnpim_IMnpip_woK0_woSid_won_data");IMnpim_IMnpip_woK0_woSid_won_data->Print();
  TH2D* IMnpim_IMnpip_wK0_woSid_won_data = (TH2D*)f->Get("IMnpim_IMnpip_wK0_woSid_won_data");IMnpim_IMnpip_wK0_woSid_won_data->Print();

  TH2D* IMnpim_IMnpip_woK0_woSid_won_mc = (TH2D*)f->Get("IMnpim_IMnpip_woK0_woSid_won_mc");IMnpim_IMnpip_woK0_woSid_won_mc->Print();

  TH2D* IMnpim_IMnpip_wK0_woSid_won_mc = (TH2D*)f->Get("IMnpim_IMnpip_wK0_woSid_won_mc");IMnpim_IMnpip_wK0_woSid_won_mc->Print();
  //TH2D* IMnpim_IMnpip_wK0_woSid_won_mcgeta = (TH2D*)f->Get("IMnpim_IMnpip_wK0_woSid_won_mcgeta");IMnpim_IMnpip_wK0_woSid_won_mcgeta->Print();

  TH2D* IMnpim_IMnpip_woSid_won_data = (TH2D*)IMnpim_IMnpip_woK0_woSid_won_data->Clone("IMnpim_IMnpip_woSid_won_data");
  IMnpim_IMnpip_woSid_won_data->Add(IMnpim_IMnpip_wK0_woSid_won_data);
  TH2D* IMnpim_IMnpip_woSid_won_mc = (TH2D*)IMnpim_IMnpip_woK0_woSid_won_mc->Clone("IMnpim_IMnpip_woSid_won_mc");
  IMnpim_IMnpip_woSid_won_mc->Add(IMnpim_IMnpip_wK0_woSid_won_mc);
  //IMnpim_IMnpip_woSid_won_mc->Add(IMnpim_IMnpip_wK0_woSid_won_mcgeta);
   
  
  TCanvas *cIMnpim_IMnpip_woSid_won = new TCanvas("cIMnpim_IMnpip_woSid_won","cIMnpim_IMnpip_woSid_won",1200,1000);
  TPad *center_pad2 = new TPad("center_pad2", "center_pad2",0.0,0.0,0.6,0.6);
  center_pad2->Draw();
  TPad *right_pad2 = new TPad("right_pad2", "right_pad2",0.55,0.0,1.0,0.6);
  right_pad2->Draw();

  TPad  *top_pad2 = new TPad("top_pad2", "top_pad2",0.0,0.55,0.6,1.0);
  top_pad2->Draw();
 
  center_pad2->cd();
  IMnpim_IMnpip_woSid_won_data->RebinX(4);
  IMnpim_IMnpip_woSid_won_data->RebinY(4);
  IMnpim_IMnpip_woSid_won_data->GetXaxis()->SetRangeUser(1,1.7);
  IMnpim_IMnpip_woSid_won_data->GetYaxis()->SetRangeUser(1,1.7);
  IMnpim_IMnpip_woSid_won_mc->RebinX(4);
  IMnpim_IMnpip_woSid_won_mc->RebinY(4);
  IMnpim_IMnpip_woSid_won_mc->GetXaxis()->SetRangeUser(1,1.7);
  IMnpim_IMnpip_woSid_won_mc->GetYaxis()->SetRangeUser(1,1.7);
  //IMnpim_IMnpip_woSid_won_data->GetZaxis()->SetNdivisions(503);
  //IMnpim_IMnpip_woSid_won_mc->GetZaxis()->SetNdivisions(503);
  IMnpim_IMnpip_woSid_won_data->SetContour(7);
  IMnpim_IMnpip_woSid_won_mc->SetContour(7);

  IMnpim_IMnpip_woSid_won_data->SetLineWidth(0.5);
  IMnpim_IMnpip_woSid_won_data->Draw("cont3");
  IMnpim_IMnpip_woSid_won_mc->SetLineColor(2);
  IMnpim_IMnpip_woSid_won_mc->SetLineWidth(0.5);
  IMnpim_IMnpip_woSid_won_mc->Draw("cont3same");
  TPaletteAxis *palette2 = (TPaletteAxis*)IMnpim_IMnpip_woSid_won_data->GetListOfFunctions()->FindObject("palette");
  // the following lines moe the paletter. Choose the values you need for the position.
  palette2->SetX1NDC(0.78);
  palette2->SetX2NDC(0.83);
  palette2->SetY1NDC(0.1);
  palette2->SetY2NDC(0.90);
  gPad->Modified();
  gPad->Update();
   
  
  top_pad2->cd();
  TH1D* IMnpip_woSid_won_data = (TH1D*)IMnpim_IMnpip_woSid_won_data->ProjectionX("IMnpip_woSid_won_data");
  TH1D* IMnpip_woSid_won_mc = (TH1D*)IMnpim_IMnpip_woSid_won_mc->ProjectionX("IMnpip_woSid_won_mc");
  top_pad2->SetBottomMargin(0);
  IMnpip_woSid_won_data->GetXaxis()->SetLabelSize(0);
  IMnpip_woSid_won_data->SetMarkerStyle(20);
  IMnpip_woSid_won_data->SetMarkerColor(1);
  IMnpip_woSid_won_data->Draw("E");
  IMnpip_woSid_won_mc->SetLineColor(2);
  IMnpip_woSid_won_mc->SetMarkerStyle(20);
  IMnpip_woSid_won_mc->SetMarkerColor(2);
  IMnpip_woSid_won_mc->Draw("Esame");
  
  
  right_pad2->cd();
  right_pad2->SetLeftMargin(0);
  TH1D* IMnpim_woSid_won_data = (TH1D*)IMnpim_IMnpip_woSid_won_data->ProjectionY("IMnpim_woSid_won_data");
  TH1D* IMnpim_woSid_won_mc = (TH1D*)IMnpim_IMnpip_woSid_won_mc->ProjectionY("IMnpim_woSid_won_mc");
  IMnpim_woSid_won_data->SetTitle("");
  IMnpim_woSid_won_data->GetXaxis()->SetTitle("");
  
  TGraphErrors *gr_IMnpim_woSid_won_data = new TGraphErrors();
  TGraphErrors *gr_IMnpim_woSid_won_mc = new TGraphErrors();
  
  for(int ibin=0;ibin<IMnpim_woSid_won_data->GetNbinsX();ibin++){
    double cont = IMnpim_woSid_won_data->GetBinContent(ibin);
    double err = IMnpim_woSid_won_data->GetBinError(ibin);
    double bincenter = IMnpim_woSid_won_data->GetBinCenter(ibin);
    gr_IMnpim_woSid_won_data->SetPoint(ibin,cont, bincenter);
    gr_IMnpim_woSid_won_data->SetPointError(ibin,err,0);
  }
  for(int ibin=0;ibin<IMnpim_woSid_won_mc->GetNbinsX();ibin++){
    double cont = IMnpim_woSid_won_mc->GetBinContent(ibin);
    double err = IMnpim_woSid_won_mc->GetBinError(ibin);
    double bincenter = IMnpim_woSid_won_mc->GetBinCenter(ibin);
    gr_IMnpim_woSid_won_mc->SetPoint(ibin,cont, bincenter);
    gr_IMnpim_woSid_won_mc->SetPointError(ibin,err,0);
  }
  gr_IMnpim_woSid_won_data->GetYaxis()->SetRangeUser(1,1.7);
  gr_IMnpim_woSid_won_data->SetMarkerStyle(20);
  gr_IMnpim_woSid_won_data->GetYaxis()->SetLabelSize(0);
  gr_IMnpim_woSid_won_data->Draw("AP");
  gr_IMnpim_woSid_won_mc->SetMarkerStyle(20);
  gr_IMnpim_woSid_won_mc->SetMarkerColor(2);
  gr_IMnpim_woSid_won_mc->Draw("P");

  TCanvas *cIMnpim_IMnpip_woSid_won_sub = new TCanvas("cIMnpim_IMnpip_woSid_won_sub","cIMnpim_IMnpip_woSid_won_sub");
  TH2D* IMnpim_IMnpip_woSid_won_sub = (TH2D*)IMnpim_IMnpip_woSid_won_data->Clone("IMnpim_IMnpip_woSid_won_data");
  IMnpim_IMnpip_woSid_won_sub->Add(IMnpim_IMnpip_woSid_won_mc,-1.0);
  IMnpim_IMnpip_woSid_won_sub->Divide(IMnpim_IMnpip_woSid_won_data);
  IMnpim_IMnpip_woSid_won_sub->SetTitle("(data-MC)/data");
  IMnpim_IMnpip_woSid_won_sub->Draw("colz");
  
   gStyle->SetPadBorderMode(0);
  //
  //BG
  //
  TH2D* MMnmiss_IMpippim_woK0_woSid_won_data = (TH2D*)f->Get("MMnmiss_IMpippim_woK0_woSid_won_data");MMnmiss_IMpippim_woK0_woSid_won_data->Print();
  TH2D* MMnmiss_IMpippim_wK0_woSid_won_data = (TH2D*)f->Get("MMnmiss_IMpippim_wK0_woSid_won_data");MMnmiss_IMpippim_wK0_woSid_won_data->Print();

  TH2D* MMnmiss_IMpippim_woK0_woSid_won_mc = (TH2D*)f->Get("MMnmiss_IMpippim_woK0_woSid_won_mc");MMnmiss_IMpippim_woK0_woSid_won_mc->Print();

  TH2D* MMnmiss_IMpippim_wK0_woSid_won_mc = (TH2D*)f->Get("MMnmiss_IMpippim_wK0_woSid_won_mc");MMnmiss_IMpippim_wK0_woSid_won_mc->Print();
  //TH2D* MMnmiss_IMpippim_wK0_woSid_won_mcgeta = (TH2D*)f->Get("MMnmiss_IMpippim_wK0_woSid_won_mcgeta");MMnmiss_IMpippim_wK0_woSid_won_mcgeta->Print();

  TH2D* MMnmiss_IMpippim_woSid_won_data = (TH2D*)MMnmiss_IMpippim_woK0_woSid_won_data->Clone("MMnmiss_IMpippim_woSid_won_data");
  MMnmiss_IMpippim_woSid_won_data->Add(MMnmiss_IMpippim_wK0_woSid_won_data);
  TH2D* MMnmiss_IMpippim_woSid_won_mc = (TH2D*)MMnmiss_IMpippim_woK0_woSid_won_mc->Clone("MMnmiss_IMpippim_woSid_won_mc");
  MMnmiss_IMpippim_woSid_won_mc->Add(MMnmiss_IMpippim_wK0_woSid_won_mc);
  //MMnmiss_IMpippim_woSid_won_mc->Add(MMnmiss_IMpippim_wK0_woSid_won_mcgeta);
   
  
  TCanvas *cMMnmiss_IMpippim_woSid_won = new TCanvas("cMMnmiss_IMpippim_woSid_won","cMMnmiss_IMpippim_woSid_won",1200,1000);
  cMMnmiss_IMpippim_woSid_won ->Divide(2,2,0,0);
  
  //TPad *center_pad3 = new TPad("center_pad2", "center_pad2",0.0,0.0,0.6,0.6);
  //center_pad3->Draw();
  //TPad *right_pad3 = new TPad("right_pad2", "right_pad2",0.55,0.0,1.0,0.6);
  //right_pad3->Draw();

  //TPad  *top_pad3 = new TPad("top_pad2", "top_pad2",0.0,0.55,0.6,1.0);
  //top_pad3->Draw();
 
  //center_pad3->cd();
  //center_pad3->SetTopMargin(0.1);
  //center_pad3->SetRightMargin(0.1);
  cMMnmiss_IMpippim_woSid_won->cd(3);
  MMnmiss_IMpippim_woSid_won_data->SetTitle("");
  MMnmiss_IMpippim_woSid_won_data->RebinX(4);
  MMnmiss_IMpippim_woSid_won_data->RebinY(3);
  MMnmiss_IMpippim_woSid_won_data->GetXaxis()->SetRangeUser(0.2,1.0);
  //MMnmiss_IMpippim_woSid_won_data->GetYaxis()->SetRangeUser(1,1.7);
  MMnmiss_IMpippim_woSid_won_mc->RebinX(4);
  MMnmiss_IMpippim_woSid_won_mc->RebinY(3);
  MMnmiss_IMpippim_woSid_won_mc->GetXaxis()->SetRangeUser(0.2,1.0);
  //MMnmiss_IMpippim_woSid_won_mc->GetYaxis()->SetRangeUser(1,1.7);
  //MMnmiss_IMpippim_woSid_won_data->GetZaxis()->SetNdivisions(503);
  //MMnmiss_IMpippim_woSid_won_mc->GetZaxis()->SetNdivisions(503);
  MMnmiss_IMpippim_woSid_won_data->SetContour(9);
  MMnmiss_IMpippim_woSid_won_mc->SetContour(9);

  MMnmiss_IMpippim_woSid_won_data->SetLineWidth(0.5);
  MMnmiss_IMpippim_woSid_won_data->Draw("cont3");
  MMnmiss_IMpippim_woSid_won_mc->SetLineColor(2);
  MMnmiss_IMpippim_woSid_won_mc->SetLineWidth(0.5);
  MMnmiss_IMpippim_woSid_won_mc->Draw("cont3same");
  gPad->Modified();
  gPad->Update();
   
  
  //top_pad3->cd();
  cMMnmiss_IMpippim_woSid_won->cd(1);
  TH1D* IMpippim_woSid_won_data = (TH1D*)MMnmiss_IMpippim_woSid_won_data->ProjectionX("IMpippim_woSid_won_data");
  TH1D* IMpippim_woSid_won_mc = (TH1D*)MMnmiss_IMpippim_woSid_won_mc->ProjectionX("IMpippim_woSid_won_mc");
  //top_pad3->SetBottomMargin(0.0);
  //top_pad3->SetRightMargin(0);
  IMpippim_woSid_won_data->SetTitle("");
  IMpippim_woSid_won_data->GetXaxis()->SetLabelSize(0);
  IMpippim_woSid_won_data->SetMarkerStyle(20);
  IMpippim_woSid_won_data->SetMarkerColor(1);
  IMpippim_woSid_won_data->Draw("E");
  IMpippim_woSid_won_mc->SetLineColor(2);
  IMpippim_woSid_won_mc->SetMarkerStyle(20);
  IMpippim_woSid_won_mc->SetMarkerColor(2);
  IMpippim_woSid_won_mc->Draw("Esame");
  
  
  //right_pad3->cd();
  //right_pad3->SetTopMargin(0);
  //right_pad3->SetLeftMargin(0);
  cMMnmiss_IMpippim_woSid_won->cd(4);
  TH1D* MMnmiss_woSid_won_data = (TH1D*)MMnmiss_IMpippim_woSid_won_data->ProjectionY("MMnmiss_woSid_won_data");
  TH1D* MMnmiss_woSid_won_mc = (TH1D*)MMnmiss_IMpippim_woSid_won_mc->ProjectionY("MMnmiss_woSid_won_mc");
  MMnmiss_woSid_won_data->SetTitle("");
  MMnmiss_woSid_won_data->GetXaxis()->SetTitle("");
  
  TGraphErrors *gr_MMnmiss_woSid_won_data = new TGraphErrors();
  TGraphErrors *gr_MMnmiss_woSid_won_mc = new TGraphErrors();
  
  for(int ibin=0;ibin<MMnmiss_woSid_won_data->GetNbinsX();ibin++){
    double cont = MMnmiss_woSid_won_data->GetBinContent(ibin);
    double err = MMnmiss_woSid_won_data->GetBinError(ibin);
    double bincenter = MMnmiss_woSid_won_data->GetBinCenter(ibin);
    gr_MMnmiss_woSid_won_data->SetPoint(ibin,cont, bincenter);
    gr_MMnmiss_woSid_won_data->SetPointError(ibin,err,0);
  }
  for(int ibin=0;ibin<MMnmiss_woSid_won_mc->GetNbinsX();ibin++){
    double cont = MMnmiss_woSid_won_mc->GetBinContent(ibin);
    double err = MMnmiss_woSid_won_mc->GetBinError(ibin);
    double bincenter = MMnmiss_woSid_won_mc->GetBinCenter(ibin);
    gr_MMnmiss_woSid_won_mc->SetPoint(ibin,cont, bincenter);
    gr_MMnmiss_woSid_won_mc->SetPointError(ibin,err,0);
  }
  gr_MMnmiss_woSid_won_data->GetYaxis()->SetRangeUser(0,1.5);
  gr_MMnmiss_woSid_won_data->SetMarkerStyle(20);
  gr_MMnmiss_woSid_won_data->GetYaxis()->SetLabelSize(0);
  gr_MMnmiss_woSid_won_data->Draw("AP");
  gr_MMnmiss_woSid_won_mc->SetMarkerStyle(20);
  gr_MMnmiss_woSid_won_mc->SetMarkerColor(2);
  gr_MMnmiss_woSid_won_mc->Draw("P");

  TCanvas *cMMnmiss_IMpippim_woSid_won_sub = new TCanvas("cMMnmiss_IMpippim_woSid_won_sub","cMMnmiss_IMpippim_woSid_won_sub");
  TH2D* MMnmiss_IMpippim_woSid_won_sub = (TH2D*)MMnmiss_IMpippim_woSid_won_data->Clone("MMnmiss_IMpippim_woSid_won_data");
  MMnmiss_IMpippim_woSid_won_sub->Add(MMnmiss_IMpippim_woSid_won_mc,-1.0);
  MMnmiss_IMpippim_woSid_won_sub->Divide(MMnmiss_IMpippim_woSid_won_data);
  MMnmiss_IMpippim_woSid_won_sub->SetTitle("(data-MC)/data");
  MMnmiss_IMpippim_woSid_won_sub->Draw("colz");
   


}
