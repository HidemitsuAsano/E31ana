#include "../../post/weightfuncGSm.h"


void disp_comp(const char *filename="comp_fakedata_out.root")
{
  TFile *f = TFile::Open(filename,"READ");
  f->cd();

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TH1D* q_woK0_woSid_won_data = (TH1D*)f->Get("q_woK0_woSid_won_data");
  TH1D* q_woK0_woSid_won_mc = (TH1D*)f->Get("q_woK0_woSid_won_mc");
  TH1D* q_wK0_woSid_won_data = (TH1D*)f->Get("q_wK0_woSid_won_data");
  TH1D* q_wK0_woSid_won_mc = (TH1D*)f->Get("q_wK0_woSid_won_mc");
  TH1D* q_wK0_woSid_won_mcgeta = (TH1D*)f->Get("q_wK0_woSid_won_mcgeta");
  TH1D* q_woK0_woSid_won_ratio = (TH1D*)f->Get("q_woK0_woSid_won_ratio");
  TH1D* q_wK0_woSid_won_ratio = (TH1D*)f->Get("q_wK0_woSid_won_ratio");
  
  TCanvas *cq = new TCanvas("cq","cq",1000,1000);
  cq->Divide(2,2);
  cq->cd(1);
  q_woK0_woSid_won_data->Draw("HE");
  q_woK0_woSid_won_mc->SetLineColor(2);
  q_woK0_woSid_won_mc->Draw("HEsame");
  cq->cd(2);
  q_wK0_woSid_won_data->Draw("HE");
  q_wK0_woSid_won_mc->SetLineColor(3);
  q_wK0_woSid_won_mc->Draw("HEsame");
  q_wK0_woSid_won_mcgeta->SetLineColor(2);
  q_wK0_woSid_won_mcgeta->Draw("HEsame");
  TH1D* q_wK0_woSid_won_mcsum = (TH1D*)q_wK0_woSid_won_mc->Clone("q_wK0_woSid_won_mcsum");
  q_wK0_woSid_won_mcsum->Add(q_wK0_woSid_won_mcgeta,1.0);
  q_wK0_woSid_won_mcsum->SetLineColor(6);
  q_wK0_woSid_won_mcsum->Draw("HEsame");
  cq->cd(3);
  q_woK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,2);
  q_woK0_woSid_won_ratio->Draw("HE");
  TCanvas *c_q_mod = new TCanvas("c_q_mod","c_q_mod");
  q_woK0_woSid_won_ratio->Draw("HE");
  TF1 *f_q_mod = new TF1("f_q_mod",func_q_mod,0.0,1.5,8);
  //q_woK0_woSid_won_ratio->Fit("f_q_mod","","",0,1.5);
  f_q_mod->Draw("same");
  cq->cd(4);
  q_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,2);
  q_wK0_woSid_won_ratio->Draw("HE");
  cq->Update();
  cq->Modified();
  //TCanvas *c_q_wK0_mod = new TCanvas("c_q_wK0_mod","c_q_wK0_mod");
  //q_wK0_woSid_won_ratio->Draw("HE");
  //TF1 *f_q_wK0_mod = new TF1("f_q_wK0_mod",func_q_mod,0.0,1.5,8);
  //q_wK0_woSid_won_ratio->Fit("f_q_wK0_mod","","",0,1.35);
  //f_q_wK0_mod->Draw("same");
  //TF1 *f_q_wK0_mod_ref = new TF1("f_q_wK0_mod_ref",func_q_mod,0.0,1.5,8);
  //f_q_wK0_mod_ref->SetParameters(param_q_wK0_mod);
  //f_q_wK0_mod_ref->SetLineColor(3);
  //f_q_wK0_mod_ref->Draw("same");
   
  TCanvas *cq_sum = new TCanvas("cq_sum","cq_sum",1000,1000);
  TH1D* q_woSid_won_data = (TH1D*) q_woK0_woSid_won_data->Clone("q_woSid_won_data");
  q_woSid_won_data->Add(q_wK0_woSid_won_data);
  q_woSid_won_data->Draw("HE");
  TH1D* q_woSid_won_pipinxmc = (TH1D*)q_woK0_woSid_won_mc->Clone("q_woSid_won_pipinxmc");
  q_woSid_won_pipinxmc->Add(q_wK0_woSid_won_mcgeta);
  q_woSid_won_pipinxmc->SetLineColor(2);
  q_woSid_won_pipinxmc->Draw("HEsame");
  q_wK0_woSid_won_mc->SetLineColor(3);
  q_wK0_woSid_won_mc->Draw("HEsame");
  TH1D* q_woSid_won_mcsum = (TH1D*)q_woSid_won_pipinxmc->Clone("q_woSid_won_mcsum");
  q_woSid_won_mcsum->Add(q_wK0_woSid_won_mc);
  q_woSid_won_mcsum->SetLineColor(6);
  q_woSid_won_mcsum->Draw("HEsame");


  TH1D* MMnmiss_woK0_woSid_won_data = (TH1D*)f->Get("MMnmiss_woK0_woSid_won_data");
  TH1D* MMnmiss_woK0_woSid_won_mc = (TH1D*)f->Get("MMnmiss_woK0_woSid_won_mc");
  TH1D* MMnmiss_wK0_woSid_won_data = (TH1D*)f->Get("MMnmiss_wK0_woSid_won_data");
  TH1D* MMnmiss_wK0_woSid_won_mc = (TH1D*)f->Get("MMnmiss_wK0_woSid_won_mc");
  TH1D* MMnmiss_wK0_woSid_won_mcgeta = (TH1D*)f->Get("MMnmiss_wK0_woSid_won_mcgeta");
  TH1D* MMnmiss_woK0_woSid_won_ratio = (TH1D*)f->Get("MMnmiss_woK0_woSid_won_ratio");
  TH1D* MMnmiss_wK0_woSid_won_ratio = (TH1D*)f->Get("MMnmiss_wK0_woSid_won_ratio");
  
  TCanvas *cMMnmiss = new TCanvas("cMMnmiss","cMMnmiss",1000,1000);
  cMMnmiss->Divide(2,2);
  cMMnmiss->cd(1);
  MMnmiss_woK0_woSid_won_data->Draw("HE");
  MMnmiss_woK0_woSid_won_mc->SetLineColor(2);
  MMnmiss_woK0_woSid_won_mc->Draw("HEsame");
  cMMnmiss->cd(2);
  MMnmiss_wK0_woSid_won_data->Draw("HE");
  MMnmiss_wK0_woSid_won_mc->SetLineColor(3);
  MMnmiss_wK0_woSid_won_mc->Draw("HEsame");
  MMnmiss_wK0_woSid_won_mcgeta->SetLineColor(2);
  MMnmiss_wK0_woSid_won_mcgeta->Draw("HEsame");
  TH1D* MMnmiss_wK0_woSid_won_mcsum = (TH1D*)MMnmiss_wK0_woSid_won_mc->Clone("MMnmiss_wK0_woSid_won_mcsum");
  MMnmiss_wK0_woSid_won_mcsum->Add(MMnmiss_wK0_woSid_won_mcgeta,1.0);
  MMnmiss_wK0_woSid_won_mcsum->SetLineColor(6);
  MMnmiss_wK0_woSid_won_mcsum->Draw("HEsame");
  cMMnmiss->cd(3);
  //MMnmiss_woK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,2);
  MMnmiss_woK0_woSid_won_ratio->Draw("HE");
  TCanvas *cMMnmiss_mod = new TCanvas("cMMnmiss_mod","cMMnmiss_mod");
  MMnmiss_woK0_woSid_won_ratio->Draw("HE");
  TF1 *f_MMnmiss_mod = new TF1("f_MMnmiss_mod",func_MMnmiss_mod,0.0,1.5,20);
  f_MMnmiss_mod->SetParameters(param_MMnmiss_mod);
  f_MMnmiss_mod->FixParameter(10,0.018);
  //f_MMnmiss_mod->SetParLimits(11,0,28000);
  //f_MMnmiss_mod->SetParLimits(12,0,2800);
  //f_MMnmiss_mod->SetParLimits(13,0,2800);
  f_MMnmiss_mod->FixParameter(14,0.010);
  MMnmiss_woK0_woSid_won_ratio->Fit("f_MMnmiss_mod","","",0.0,1.5);
  //f_MMnmiss_mod->Draw("same");
  cMMnmiss->cd(4);
  MMnmiss_wK0_woSid_won_ratio->Rebin(2);
  MMnmiss_wK0_woSid_won_ratio->Scale(0.5);
  MMnmiss_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,2);
  MMnmiss_wK0_woSid_won_ratio->Draw("HE");
  TCanvas *cMMnmiss_wK0_mod = new TCanvas("cMMnmiss_wK0_mod","cMMnmiss_wK0_mod");
  MMnmiss_wK0_woSid_won_ratio->Draw("HE");
  TF1 *f_MMnmiss_wK0_mod = new TF1("f_MMnmiss_wK0_mod",func_MMnmiss_wK0_mod,0.0,1.28,11);
  //f_MMnmiss_wK0_mod->SetParameter(0,-0.0156614);
  //f_MMnmiss_wK0_mod->SetParameter(1,2.94145);
  //f_MMnmiss_wK0_mod->SetParameter(2,-47.1231);
  //f_MMnmiss_wK0_mod->SetParameter(3,491.928);
  //f_MMnmiss_wK0_mod->SetParameter(4,-2519.55);
  //f_MMnmiss_wK0_mod->SetParameter(5,7184.26);
  //f_MMnmiss_wK0_mod->SetParameter(6,-11965.7);
  //f_MMnmiss_wK0_mod->SetParameter(7,11581.9);
  //f_MMnmiss_wK0_mod->SetParameter(8,-6029.28);
  //f_MMnmiss_wK0_mod->SetParameter(9,1302.35);
  //f_MMnmiss_wK0_mod->FixParameter(6,0.01);
  //MMnmiss_wK0_woSid_won_ratio->Fit(f_MMnmiss_wK0_mod,"","",0.0,1.28);

  TCanvas *cMMnmiss_sum = new TCanvas("cMMnmiss_sum","cMMnmiss_sum",1000,1000);
  TH1D* MMnmiss_woSid_won_data = (TH1D*) MMnmiss_woK0_woSid_won_data->Clone("MMnmiss_woSid_won_data");
  MMnmiss_woSid_won_data->Add(MMnmiss_wK0_woSid_won_data);
  MMnmiss_woSid_won_data->Draw("HE");
  TH1D* MMnmiss_woSid_won_pipinxmc = (TH1D*)MMnmiss_woK0_woSid_won_mc->Clone("MMnmiss_woSid_won_pipinxmc");
  MMnmiss_woSid_won_pipinxmc->Add(MMnmiss_wK0_woSid_won_mcgeta);
  MMnmiss_woSid_won_pipinxmc->SetLineColor(2);
  MMnmiss_woSid_won_pipinxmc->Draw("HEsame");
  MMnmiss_wK0_woSid_won_mc->SetLineColor(3);
  MMnmiss_wK0_woSid_won_mc->Draw("HEsame");
  TH1D* MMnmiss_woSid_won_mcsum = (TH1D*)MMnmiss_woSid_won_pipinxmc->Clone("MMnmiss_woSid_won_mcsum");
  MMnmiss_woSid_won_mcsum->Add(MMnmiss_wK0_woSid_won_mc);
  MMnmiss_woSid_won_mcsum->SetLineColor(6);
  MMnmiss_woSid_won_mcsum->Draw("HEsame");

  TH1D* IMnpip_woK0_woSid_won_data = (TH1D*)f->Get("IMnpip_woK0_woSid_won_data");
  TH1D* IMnpip_woK0_woSid_won_mc = (TH1D*)f->Get("IMnpip_woK0_woSid_won_mc");
  TH1D* IMnpip_wK0_woSid_won_data = (TH1D*)f->Get("IMnpip_wK0_woSid_won_data");
  TH1D* IMnpip_wK0_woSid_won_mc = (TH1D*)f->Get("IMnpip_wK0_woSid_won_mc");
  TH1D* IMnpip_wK0_woSid_won_mcgeta = (TH1D*)f->Get("IMnpip_wK0_woSid_won_mcgeta");
  TH1D* IMnpip_woK0_woSid_won_ratio = (TH1D*)f->Get("IMnpip_woK0_woSid_won_ratio");
  TH1D* IMnpip_wK0_woSid_won_ratio = (TH1D*)f->Get("IMnpip_wK0_woSid_won_ratio");
  
  TCanvas *cIMnpip = new TCanvas("cIMnpip","cIMnpip",1000,1000);
  cIMnpip->Divide(2,2);
  cIMnpip->cd(1);
  IMnpip_woK0_woSid_won_data->Draw("HE");
  IMnpip_woK0_woSid_won_mc->SetLineColor(2);
  IMnpip_woK0_woSid_won_mc->Draw("HEsame");
  cIMnpip->cd(2);
  IMnpip_wK0_woSid_won_data->Draw("HE");
  IMnpip_wK0_woSid_won_mc->SetLineColor(3);
  IMnpip_wK0_woSid_won_mc->Draw("HEsame");
  IMnpip_wK0_woSid_won_mcgeta->SetLineColor(2);
  IMnpip_wK0_woSid_won_mcgeta->Draw("HEsame");
  TH1D* IMnpip_wK0_woSid_won_mcsum = (TH1D*)IMnpip_wK0_woSid_won_mc->Clone("IMnpip_wK0_woSid_won_mcsum");
  IMnpip_wK0_woSid_won_mcsum->Add(IMnpip_wK0_woSid_won_mcgeta,1.0);
  IMnpip_wK0_woSid_won_mcsum->SetLineColor(6);
  IMnpip_wK0_woSid_won_mcsum->Draw("HEsame");
  cIMnpip->cd(3);
  IMnpip_woK0_woSid_won_ratio->Rebin(2);
  IMnpip_woK0_woSid_won_ratio->Scale(0.5);
  IMnpip_woK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,2);
  IMnpip_woK0_woSid_won_ratio->Draw("HE");
  
  TCanvas *cIMnpip_mod = new TCanvas("cIMnpip_mod","cIMnpip_mod");
  cIMnpip_mod->cd();
  IMnpip_woK0_woSid_won_ratio->Draw("HE");
  TF1 *f_IMnpip_mod = new TF1("f_IMnpip_mod",func_IMnpip_mod,1,2,17);
  //TF1 *f_IMnpip_mod_ref = new TF1("f_IMnpip_mod_ref",func_IMnpip_mod,1,2,12);
  //f_IMnpip_mod->SetParameters(param_IMnpip_mod);
  //f_IMnpip_mod_ref->SetParameters(param_IMnpip_mod);
  f_IMnpip_mod->FixParameter(0,0.454092);
  f_IMnpip_mod->FixParameter(1,1.08315);
  f_IMnpip_mod->FixParameter(2,7.26318e-03);
  f_IMnpip_mod->FixParameter(3,0.008);
  //f_IMnpip_mod->SetParameter(3,5.75432e+01);
  //f_IMnpip_mod->SetParameter(4,96.6903);
  //f_IMnpip_mod->SetParameter(5,0.445582);
  //f_IMnpip_mod->SetParameter(6,0.261843);
  f_IMnpip_mod->FixParameter(11,0.02);
  //f_IMnpip_mod->SetParameter(9,0.326883);
  //f_IMnpip_mod->SetParameter(10,1.52682);
  //f_IMnpip_mod->SetParameter(11,0.103997);
  
  //f_IMnpip_mod_ref->SetLineColor(3);
  //f_IMnpip_mod_ref->Draw("same");
  IMnpip_woK0_woSid_won_ratio->Fit("f_IMnpip_mod","","",1.,2.0);
  
  cIMnpip->cd(4);
  //IMnpip_wK0_woSid_won_ratio->RebinX(2);
  //IMnpip_wK0_woSid_won_ratio->Scale(0.5);
  IMnpip_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,4);
  IMnpip_wK0_woSid_won_ratio->Draw("HE");
  //TCanvas *cIMnpip_wK0_mod = new TCanvas("cIMnpip_wK0_mod","cIMnpip_wK0_mod");
  //cIMnpip_wK0_mod->cd();
  //IMnpip_wK0_woSid_won_ratio->Draw("HE");
  //TF1 *f_IMnpip_wK0 = new TF1("f_IMnpip_wK0",func_IMnpip_wK0_mod,1,2,14);
  //f_IMnpip_wK0->SetParameters(param_IMnpip_wK0_mod);
  //f_IMnpip_wK0->FixParameter(4,0.01);
  //f_IMnpip_wK0->FixParameter(8,0.01);

  //f_IMnpip_wK0->SetLineColor(2);
  //f_IMnpip_wK0->Draw("same");
  //IMnpip_wK0_woSid_won_ratio->Fit("f_IMnpip_wK0","","",1.,1.8);
  
  
  TCanvas *cIMnpip_sum = new TCanvas("cIMnpip_sum","cIMnpip_sum",1000,1000);
  TH1D* IMnpip_woSid_won_data = (TH1D*) IMnpip_woK0_woSid_won_data->Clone("IMnpip_woSid_won_data");
  IMnpip_woSid_won_data->Add(IMnpip_wK0_woSid_won_data);
  IMnpip_woSid_won_data->Draw("HE");
  TH1D* IMnpip_woSid_won_pipinxmc = (TH1D*)IMnpip_woK0_woSid_won_mc->Clone("IMnpip_woSid_won_pipinxmc");
  IMnpip_woSid_won_pipinxmc->Add(IMnpip_wK0_woSid_won_mcgeta);
  IMnpip_woSid_won_pipinxmc->SetLineColor(2);
  IMnpip_woSid_won_pipinxmc->Draw("HEsame");
  IMnpip_wK0_woSid_won_mc->SetLineColor(3);
  IMnpip_wK0_woSid_won_mc->Draw("HEsame");
  TH1D* IMnpip_woSid_won_mcsum = (TH1D*)IMnpip_woSid_won_pipinxmc->Clone("IMnpip_woSid_won_mcsum");
  IMnpip_woSid_won_mcsum->Add(IMnpip_wK0_woSid_won_mc);
  IMnpip_woSid_won_mcsum->SetLineColor(6);
  IMnpip_woSid_won_mcsum->Draw("HEsame");
  

  TH1D* IMnpim_woK0_woSid_won_data = (TH1D*)f->Get("IMnpim_woK0_woSid_won_data");
  TH1D* IMnpim_woK0_woSid_won_mc = (TH1D*)f->Get("IMnpim_woK0_woSid_won_mc");
  TH1D* IMnpim_wK0_woSid_won_data = (TH1D*)f->Get("IMnpim_wK0_woSid_won_data");
  TH1D* IMnpim_wK0_woSid_won_mc = (TH1D*)f->Get("IMnpim_wK0_woSid_won_mc");
  TH1D* IMnpim_wK0_woSid_won_mcgeta = (TH1D*)f->Get("IMnpim_wK0_woSid_won_mcgeta");
  TH1D* IMnpim_woK0_woSid_won_ratio = (TH1D*)f->Get("IMnpim_woK0_woSid_won_ratio");
  TH1D* IMnpim_wK0_woSid_won_ratio = (TH1D*)f->Get("IMnpim_wK0_woSid_won_ratio");
  
  TCanvas *cIMnpim = new TCanvas("cIMnpim","cIMnpim",1000,1000);
  cIMnpim->Divide(2,2);
  cIMnpim->cd(1);
  IMnpim_woK0_woSid_won_data->Draw("HE");
  IMnpim_woK0_woSid_won_mc->SetLineColor(2);
  IMnpim_woK0_woSid_won_mc->Draw("HEsame");
  cIMnpim->cd(2);
  IMnpim_wK0_woSid_won_data->Draw("HE");
  IMnpim_wK0_woSid_won_mc->SetLineColor(3);
  IMnpim_wK0_woSid_won_mc->Draw("HEsame");
  IMnpim_wK0_woSid_won_mcgeta->SetLineColor(2);
  IMnpim_wK0_woSid_won_mcgeta->Draw("HEsame");
  std::cout << __LINE__ << std::endl;
  TH1D* IMnpim_wK0_woSid_won_mcsum = IMnpim_wK0_woSid_won_mcsum = (TH1D*)IMnpim_wK0_woSid_won_mc->Clone("IMnpim_wK0_woSid_won_mcsum");
  IMnpim_wK0_woSid_won_mcsum->Add(IMnpim_wK0_woSid_won_mcgeta,1.0);
  IMnpim_wK0_woSid_won_mcsum->SetLineColor(6);
  IMnpim_wK0_woSid_won_mcsum->Draw("HEsame");
  cIMnpim->cd(3);
  //IMnpim_woK0_woSid_won_ratio->RebinX(2);
  //IMnpim_woK0_woSid_won_ratio->Scale(0.5);
  IMnpim_woK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,2.5);
  IMnpim_woK0_woSid_won_ratio->Draw("HE");
  TCanvas *cIMnpim_mod = new TCanvas("cIMnpim_mod","cIMnpim_mod");
  cIMnpim_mod->cd();
  IMnpim_woK0_woSid_won_ratio->Draw("HE");
  TF1 *f_IMnpim_mod = new TF1("f_IMnpim_mod",func_IMnpim_mod,1,2.0,14);
  //TF1 *f_IMnpim_mod = new TF1("f_IMnpim_mod",func_IMnpim_mod,1,2.0,6);
  //TF1 *f_IMnpim_mod_ref = new TF1("f_IMnpim_mod_ref",func_IMnpim_mod,1,2.0,12);
  //f_IMnpim_mod->SetParameters(param_IMnpim_mod);
  //f_IMnpim_mod_ref->SetParameters(param_IMnpim_mod);
  f_IMnpim_mod->FixParameter(0,-6053.91);
  f_IMnpim_mod->FixParameter(1,15687.3);
  f_IMnpim_mod->FixParameter(2,-13543.8);
  f_IMnpim_mod->FixParameter(3,3897.16);
  f_IMnpim_mod->FixParameter(4,0.010); 
  f_IMnpim_mod->FixParameter(5,3380.89); 
  f_IMnpim_mod->FixParameter(6,-7965.97); 
  f_IMnpim_mod->FixParameter(7,4858.63); 
  f_IMnpim_mod->FixParameter(8,1700.92); 
  f_IMnpim_mod->FixParameter(9,-2654.37); 
  f_IMnpim_mod->FixParameter(10,702.541); 
  f_IMnpim_mod->FixParameter(11,0.010); 
  f_IMnpim_mod->FixParameter(12,14.4862); 
  f_IMnpim_mod->FixParameter(13,-11.6483); 


  //f_IMnpim_mod->FixParameter(0,1.98042);
  //f_IMnpim_mod->FixParameter(1,1.11227);
  //f_IMnpim_mod->FixParameter(2,0.0240267);
  //f_IMnpim_mod->FixParameter(3,0.005); 
  IMnpim_woK0_woSid_won_ratio->Fit("f_IMnpim_mod","","",1.07,1.8);
  //f_IMnpim_mod_ref->SetLineColor(3);
  //f_IMnpim_mod_ref->Draw("same");
  cIMnpim->cd(4);
  IMnpim_wK0_woSid_won_ratio->RebinX(4);
  IMnpim_wK0_woSid_won_ratio->Scale(0.25);
  IMnpim_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,2);
  IMnpim_wK0_woSid_won_ratio->Draw("HE");
  TCanvas *cIMnpim_wK0_mod = new TCanvas("cIMnpim_wK0_mod","cIMnpim_wK0_mod");
  cIMnpim_wK0_mod->cd();
  IMnpim_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  IMnpim_wK0_woSid_won_ratio->Draw("HE");
  //TF1 *f_IMnpim_wK0_mod = new TF1("f_IMnpim_wK0_mod",func_IMnpim_wK0_mod,1,2.0,11);
  TF1 *f_IMnpim_wK0_mod = new TF1("f_IMnpim_wK0_mod",func_IMnpim_wK0_mod,1,2.0,10);
  //TF1 *f_IMnpim_mod = new TF1("f_IMnpim_mod",func_IMnpim_mod,1,2.0,10);
  //TF1 *f_IMnpim_mod_ref = new TF1("f_IMnpim_mod_ref",func_IMnpim_mod,1,2.0,10);
  f_IMnpim_wK0_mod->SetParameters(param_IMnpim_wK0_mod);
  f_IMnpim_wK0_mod->FixParameter(3,0.005);
  f_IMnpim_wK0_mod->FixParameter(7,0.005);
  //f_IMnpim_wK0_mod->Print();
  //f_IMnpim_mod_ref->SetParameters(param_IMnpim_mod);
  f_IMnpim_wK0_mod->FixParameter(0,0.0858578);
  f_IMnpim_wK0_mod->FixParameter(1,1.10463);
  f_IMnpim_wK0_mod->FixParameter(2,0.00904288);
  //f_IMnpim_wK0_mod->FixParameter(4,3296.35);
  //f_IMnpim_wK0_mod->FixParameter(5,-8069.16);
  //f_IMnpim_wK0_mod->FixParameter(6,6579.4);
  //f_IMnpim_wK0_mod->FixParameter(7,-1786.53);
  f_IMnpim_wK0_mod->FixParameter(8, 8.48006e+00);
  f_IMnpim_wK0_mod->FixParameter(9,-6.98101e+00);

  //f_IMnpim_mod->FixParameter(3,0.02); 
  //IMnpim_wK0_woSid_won_ratio->Fit("f_IMnpim_wK0_mod","","",1.07,1.8);
  //f_IMnpim_wK0_mod->SetNpx(1000);
  //f_IMnpim_wK0_mod->SetLineColor(3);
  //f_IMnpim_wK0_mod->Draw("same");

  TCanvas *cIMnpim_sum = new TCanvas("cIMnpim_sum","cIMnpim_sum",1000,1000);
  TH1D* IMnpim_woSid_won_data = (TH1D*) IMnpim_woK0_woSid_won_data->Clone("IMnpim_woSid_won_data");
  IMnpim_woSid_won_data->Add(IMnpim_wK0_woSid_won_data);
  IMnpim_woSid_won_data->Draw("HE");
  TH1D* IMnpim_woSid_won_pipinxmc = (TH1D*)IMnpim_woK0_woSid_won_mc->Clone("IMnpim_woSid_won_pipinxmc");
  IMnpim_woSid_won_pipinxmc->Add(IMnpim_wK0_woSid_won_mcgeta);
  IMnpim_woSid_won_pipinxmc->SetLineColor(2);
  IMnpim_woSid_won_pipinxmc->Draw("HEsame");
  IMnpim_wK0_woSid_won_mc->SetLineColor(3);
  IMnpim_wK0_woSid_won_mc->Draw("HEsame");
  TH1D* IMnpim_woSid_won_mcsum = (TH1D*)IMnpim_woSid_won_pipinxmc->Clone("IMnpim_woSid_won_mcsum");
  IMnpim_woSid_won_mcsum->Add(IMnpim_wK0_woSid_won_mc);
  IMnpim_woSid_won_mcsum->SetLineColor(6);
  IMnpim_woSid_won_mcsum->Draw("HEsame");
  
  /*
  TH1D* Momnpip_woK0_woSid_won_data = (TH1D*)f->Get("Momnpip_woK0_woSid_won_data");
  TH1D* Momnpip_woK0_woSid_won_mc = (TH1D*)f->Get("Momnpip_woK0_woSid_won_mc");
  TH1D* Momnpip_wK0_woSid_won_data = (TH1D*)f->Get("Momnpip_wK0_woSid_won_data");
  TH1D* Momnpip_wK0_woSid_won_mc = (TH1D*)f->Get("Momnpip_wK0_woSid_won_mc");
  TH1D* Momnpip_woK0_woSid_won_ratio = (TH1D*)f->Get("Momnpip_woK0_woSid_won_ratio");
  TH1D* Momnpip_wK0_woSid_won_ratio = (TH1D*)f->Get("Momnpip_wK0_woSid_won_ratio");
  
  TCanvas *cMomnpip = new TCanvas("cMomnpip","cMomnpip",1000,1000);
  cMomnpip->Divide(2,2);
  cMomnpip->cd(1);
  Momnpip_woK0_woSid_won_data->Draw("HE");
  Momnpip_woK0_woSid_won_mc->Draw("HEsame");
  cMomnpip->cd(2);
  Momnpip_wK0_woSid_won_data->Draw("HE");
  Momnpip_wK0_woSid_won_mc->Draw("HEsame");
  cMomnpip->cd(3);
  Momnpip_woK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,2);
  Momnpip_woK0_woSid_won_ratio->Draw("HE");
  cMomnpip->cd(4);
  Momnpip_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,2);
  Momnpip_wK0_woSid_won_ratio->Draw("HE");

  TH1D* Momnpim_woK0_woSid_won_data = (TH1D*)f->Get("Momnpim_woK0_woSid_won_data");
  TH1D* Momnpim_woK0_woSid_won_mc = (TH1D*)f->Get("Momnpim_woK0_woSid_won_mc");
  TH1D* Momnpim_wK0_woSid_won_data = (TH1D*)f->Get("Momnpim_wK0_woSid_won_data");
  TH1D* Momnpim_wK0_woSid_won_mc = (TH1D*)f->Get("Momnpim_wK0_woSid_won_mc");
  TH1D* Momnpim_woK0_woSid_won_ratio = (TH1D*)f->Get("Momnpim_woK0_woSid_won_ratio");
  TH1D* Momnpim_wK0_woSid_won_ratio = (TH1D*)f->Get("Momnpim_wK0_woSid_won_ratio");
  
  TCanvas *cMomnpim = new TCanvas("cMomnpim","cMomnpim",1000,1000);
  cMomnpim->Divide(2,2);
  cMomnpim->cd(1);
  Momnpim_woK0_woSid_won_data->Draw("HE");
  Momnpim_woK0_woSid_won_mc->Draw("HEsame");
  cMomnpim->cd(2);
  Momnpim_wK0_woSid_won_data->Draw("HE");
  Momnpim_wK0_woSid_won_mc->Draw("HEsame");
  cMomnpim->cd(3);
  Momnpim_woK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,2);
  Momnpim_woK0_woSid_won_ratio->Draw("HE");
  cMomnpim->cd(4);
  Momnpim_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,2);
  Momnpim_wK0_woSid_won_ratio->Draw("HE");
  */

  TH1D* IMnpipi_woK0_woSid_won_data = (TH1D*)f->Get("IMnpipi_woK0_woSid_won_data");
  TH1D* IMnpipi_woK0_woSid_won_mc = (TH1D*)f->Get("IMnpipi_woK0_woSid_won_mc");
  TH1D* IMnpipi_wK0_woSid_won_data = (TH1D*)f->Get("IMnpipi_wK0_woSid_won_data");
  TH1D* IMnpipi_wK0_woSid_won_mc = (TH1D*)f->Get("IMnpipi_wK0_woSid_won_mc");
  TH1D* IMnpipi_wK0_woSid_won_mcgeta = (TH1D*)f->Get("IMnpipi_wK0_woSid_won_mcgeta");
  TH1D* IMnpipi_woK0_woSid_won_ratio = (TH1D*)f->Get("IMnpipi_woK0_woSid_won_ratio");
  TH1D* IMnpipi_wK0_woSid_won_ratio = (TH1D*)f->Get("IMnpipi_wK0_woSid_won_ratio");

  TCanvas *cIMnpipi = new TCanvas("cIMnpipi","cIMnpipi",1000,1000);
  cIMnpipi->Divide(2,2);
  cIMnpipi->cd(1);
  IMnpipi_woK0_woSid_won_data->Draw("HE");
  IMnpipi_woK0_woSid_won_mc->SetLineColor(2);
  IMnpipi_woK0_woSid_won_mc->Draw("HEsame");
  cIMnpipi->cd(2);
  IMnpipi_wK0_woSid_won_data->Draw("HE");
  IMnpipi_wK0_woSid_won_mc->SetLineColor(3);
  IMnpipi_wK0_woSid_won_mc->Draw("HEsame");
  IMnpipi_wK0_woSid_won_mcgeta->SetLineColor(2);
  IMnpipi_wK0_woSid_won_mcgeta->Draw("HEsame");
  TH1D* IMnpipi_wK0_woSid_won_mcsum = (TH1D*)IMnpipi_wK0_woSid_won_mc->Clone("IMnpipi_wK0_woSid_won_mcsum");
  IMnpipi_wK0_woSid_won_mcsum->Add(IMnpipi_wK0_woSid_won_mcgeta,1.0);
  IMnpipi_wK0_woSid_won_mcsum->SetLineColor(6);
  IMnpipi_wK0_woSid_won_mcsum->Draw("HEsame");
  cIMnpipi->cd(3);
  IMnpipi_woK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,2);
  IMnpipi_woK0_woSid_won_ratio->Draw("HE");
  cIMnpipi->cd(4);
  IMnpipi_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,2);
  IMnpipi_wK0_woSid_won_ratio->Draw("HE");

  /*
  TH1D* Momnpipi_woK0_woSid_won_data = (TH1D*)f->Get("Momnpipi_woK0_woSid_won_data");
  TH1D* Momnpipi_woK0_woSid_won_mc = (TH1D*)f->Get("Momnpipi_woK0_woSid_won_mc");
  TH1D* Momnpipi_wK0_woSid_won_data = (TH1D*)f->Get("Momnpipi_wK0_woSid_won_data");
  TH1D* Momnpipi_wK0_woSid_won_mc = (TH1D*)f->Get("Momnpipi_wK0_woSid_won_mc");
  TH1D* Momnpipi_woK0_woSid_won_ratio = (TH1D*)f->Get("Momnpipi_woK0_woSid_won_ratio");
  TH1D* Momnpipi_wK0_woSid_won_ratio = (TH1D*)f->Get("Momnpipi_wK0_woSid_won_ratio");

  TCanvas *cMomnpipi = new TCanvas("cMomnpipi","cMomnpipi",1000,1000);
  cMomnpipi->Divide(2,2);
  cMomnpipi->cd(1);
  Momnpipi_woK0_woSid_won_data->Draw("HE");
  Momnpipi_woK0_woSid_won_mc->Draw("HEsame");
  cMomnpipi->cd(2);
  Momnpipi_wK0_woSid_won_data->Draw("HE");
  Momnpipi_wK0_woSid_won_mc->Draw("HEsame");
  cMomnpipi->cd(3);
  Momnpipi_woK0_woSid_won_ratio->Draw("HE");
  cMomnpipi->cd(4);
  Momnpipi_wK0_woSid_won_ratio->Draw("HE");
  */

  TH1D* IMpippim_woK0_woSid_won_data = (TH1D*)f->Get("IMpippim_woK0_woSid_won_data");
  TH1D* IMpippim_woK0_woSid_won_mc = (TH1D*)f->Get("IMpippim_woK0_woSid_won_mc");
  TH1D* IMpippim_wK0_woSid_won_data = (TH1D*)f->Get("IMpippim_wK0_woSid_won_data");
  TH1D* IMpippim_wK0_woSid_won_mc = (TH1D*)f->Get("IMpippim_wK0_woSid_won_mc");
  TH1D* IMpippim_wK0_woSid_won_mcgeta = (TH1D*)f->Get("IMpippim_wK0_woSid_won_mcgeta");
  TH1D* IMpippim_woK0_woSid_won_ratio = (TH1D*)f->Get("IMpippim_woK0_woSid_won_ratio");
  TH1D* IMpippim_wK0_woSid_won_ratio = (TH1D*)f->Get("IMpippim_wK0_woSid_won_ratio");

  TCanvas *cIMpippim = new TCanvas("cIMpippim","cIMpippim",1000,1000);
  cIMpippim->Divide(2,2);
  cIMpippim->cd(1);
  //IMpippim_woK0_woSid_won_data->RebinX(4);
  IMpippim_woK0_woSid_won_data->Draw("HE");
  //IMpippim_woK0_woSid_won_mc->RebinX(4);
  IMpippim_woK0_woSid_won_mc->SetLineColor(2);
  IMpippim_woK0_woSid_won_mc->Draw("HEsame");
  cIMpippim->cd(2);
  //IMpippim_wK0_woSid_won_data->RebinX(4);
  //IMpippim_wK0_woSid_won_data->Scale(0.25);
  IMpippim_wK0_woSid_won_data->GetXaxis()->SetRangeUser(0.44,0.55);
  IMpippim_wK0_woSid_won_data->Draw("HE");

  IMpippim_wK0_woSid_won_mc->SetLineColor(3);
  //IMpippim_wK0_woSid_won_mc->RebinX(4);
  //IMpippim_wK0_woSid_won_mc->Scale(0.25);
  IMpippim_wK0_woSid_won_mc->Draw("HEsame");
  IMpippim_wK0_woSid_won_mcgeta->SetLineColor(2);
  //IMpippim_wK0_woSid_won_mcgeta->RebinX(4);
  //IMpippim_wK0_woSid_won_mcgeta->Scale(0.25);
  IMpippim_wK0_woSid_won_mcgeta->Draw("HEsame");
  TH1D* IMpippim_wK0_woSid_won_mcsum = (TH1D*)IMpippim_wK0_woSid_won_mc->Clone("IMpippim_wK0_woSid_won_mcsum");
  IMpippim_wK0_woSid_won_mcsum->Add(IMpippim_wK0_woSid_won_mcgeta,1.0);
  IMpippim_wK0_woSid_won_mcsum->SetLineColor(6);
  IMpippim_wK0_woSid_won_mcsum->Draw("HEsame");
  cIMpippim->cd(3);
  IMpippim_woK0_woSid_won_ratio->RebinX(4);
  IMpippim_woK0_woSid_won_ratio->Scale(0.25);
  IMpippim_woK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,2);
  IMpippim_woK0_woSid_won_ratio->Draw("HE");
  TCanvas *cIMpippim_mod = new TCanvas("cIMpippim_mod","cIMpippim_mod");
  cIMpippim_mod->cd();
  TF1 *f_IMpippim_mod = new TF1("f_IMpippim_mod",func_IMpippim_mod,0,1.0,15);
  f_IMpippim_mod->SetParameters(param_IMpippim_mod);
  //f_IMpippim_mod->SetParameter(0,0.606399);
  //f_IMpippim_mod->SetParameter(1,0.307431);
  //f_IMpippim_mod->SetParameter(2,0.024839);
  f_IMpippim_mod->FixParameter(3,0.015);
  f_IMpippim_mod->FixParameter(11,0.03);
  //f_IMpippim_mod->SetParameter(12,1.03224);
  //f_IMpippim_mod->SetParameter(13,0.763145);
  //f_IMpippim_mod->SetParameter(14,0.0590810);
  //IMpippim_woK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,2);
  IMpippim_woK0_woSid_won_ratio->Draw("HE");
  IMpippim_woK0_woSid_won_ratio->Fit("f_IMpippim_mod","","",0.28,1.0);
  
  cIMpippim->cd(4);
  //IMpippim_wK0_woSid_won_ratio->RebinX(2);
  //IMpippim_wK0_woSid_won_ratio->Scale(0.5);
  IMpippim_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,2);
  IMpippim_wK0_woSid_won_ratio->GetXaxis()->SetRangeUser(0.44,0.55);
  IMpippim_wK0_woSid_won_ratio->Draw("HE");
  
  TCanvas *cIMpippim_sum = new TCanvas("cIMpippim_sum","cIMpippim_sum",1000,1000);
  TH1D* IMpippim_woSid_won_data = (TH1D*) IMpippim_woK0_woSid_won_data->Clone("IMpippim_woSid_won_data");
  IMpippim_woSid_won_data->Add(IMpippim_wK0_woSid_won_data);
  IMpippim_woSid_won_data->GetXaxis()->SetRangeUser(0.2,1.0);
  IMpippim_woSid_won_data->Draw("HE");
  TH1D* IMpippim_woSid_won_pipinxmc = (TH1D*)IMpippim_woK0_woSid_won_mc->Clone("IMpippim_woSid_won_pipinxmc");
  IMpippim_woSid_won_pipinxmc->Add(IMpippim_wK0_woSid_won_mcgeta);
  IMpippim_woSid_won_pipinxmc->SetLineColor(2);
  IMpippim_wK0_woSid_won_mc->SetLineColor(3);
  IMpippim_wK0_woSid_won_mc->Draw("HEsame");
  TH1D* IMpippim_woSid_won_mcsum = (TH1D*)IMpippim_woSid_won_pipinxmc->Clone("IMpippim_woSid_won_mcsum");
  IMpippim_woSid_won_mcsum->Add(IMpippim_wK0_woSid_won_mc);
  IMpippim_woSid_won_mcsum->SetLineColor(6);
  IMpippim_woSid_won_mcsum->Draw("HEsame");
  IMpippim_woSid_won_pipinxmc->Draw("HEsame");

  /*
  TH1D* Mompippim_woK0_woSid_won_data = (TH1D*)f->Get("Mompippim_woK0_woSid_won_data");
  TH1D* Mompippim_woK0_woSid_won_mc = (TH1D*)f->Get("Mompippim_woK0_woSid_won_mc");
  TH1D* Mompippim_wK0_woSid_won_data = (TH1D*)f->Get("Mompippim_wK0_woSid_won_data");
  TH1D* Mompippim_wK0_woSid_won_mc = (TH1D*)f->Get("Mompippim_wK0_woSid_won_mc");
  TH1D* Mompippim_wK0_woSid_won_mcgeta = (TH1D*)f->Get("Mompippim_wK0_woSid_won_mcgeta");
  TH1D* Mompippim_woK0_woSid_won_ratio = (TH1D*)f->Get("Mompippim_woK0_woSid_won_ratio");
  TH1D* Mompippim_wK0_woSid_won_ratio = (TH1D*)f->Get("Mompippim_wK0_woSid_won_ratio");

  TCanvas *cMompippim = new TCanvas("cMompippim","cMompippim",1000,1000);
  cMompippim->Divide(2,2);
  cMompippim->cd(1);
  Mompippim_woK0_woSid_won_data->Draw("HE");
  Mompippim_woK0_woSid_won_mc->Draw("HEsame");
  cMompippim->cd(2);
  Mompippim_wK0_woSid_won_data->Draw("HE");
  Mompippim_wK0_woSid_won_mc->SetLineColor(3);
  Mompippim_wK0_woSid_won_mc->Draw("HEsame");
  Mompippim_wK0_woSid_won_mcgeta->SetLineColor(2);
  Mompippim_wK0_woSid_won_mcgeta->Draw("HEsame");
  TH1D* Mompippim_wK0_woSid_won_mcsum = (TH1D*)Mompippim_wK0_woSid_won_mcgeta->Clone("Mompippim_wK0_woSid_won_mcsum");
  Mompippim_wK0_woSid_won_mcsum->Add(Mompippim_wK0_woSid_won_mc);
  Mompippim_wK0_woSid_won_mcsum->SetLineColor(6);
  Mompippim_wK0_woSid_won_mcsum->Draw("HEsame");
  cMompippim->cd(3);
  Mompippim_woK0_woSid_won_ratio->Draw("HE");
  TCanvas *cMompippim_woK0_woSid_won_ratio = new TCanvas("cMompippim_woK0_woSid_won_ratio","cMompippim_woK0_woSid_won_ratio");
  cMompippim_woK0_woSid_won_ratio->cd();
  Mompippim_woK0_woSid_won_ratio->Draw("HE");
  //TF1 *f_Mompippim = new TF1("f_Mompippim",func_Mompippim,0,1,7);
  //Mompippim_woK0_woSid_won_ratio->Fit("f_Mompippim","","",0,1);
  //cMompippim->cd(4);
  //Mompippim_wK0_woSid_won_ratio->Draw("HE");
  */

  TH1D* nmom_woK0_woSid_won_data = (TH1D*)f->Get("nmom_woK0_woSid_won_data");
  TH1D* nmom_woK0_woSid_won_mc = (TH1D*)f->Get("nmom_woK0_woSid_won_mc");
  TH1D* nmom_wK0_woSid_won_data = (TH1D*)f->Get("nmom_wK0_woSid_won_data");
  TH1D* nmom_wK0_woSid_won_mc = (TH1D*)f->Get("nmom_wK0_woSid_won_mc");
  TH1D* nmom_wK0_woSid_won_mcgeta = (TH1D*)f->Get("nmom_wK0_woSid_won_mcgeta");
  TH1D* nmom_woK0_woSid_won_ratio = (TH1D*)f->Get("nmom_woK0_woSid_won_ratio");
  TH1D* nmom_wK0_woSid_won_ratio = (TH1D*)f->Get("nmom_wK0_woSid_won_ratio");

  TCanvas *cnmom = new TCanvas("cnmom","cnmom",1000,1000);
  cnmom->Divide(2,2);
  cnmom->cd(1);
  nmom_woK0_woSid_won_data->Draw("HE");
  nmom_woK0_woSid_won_mc->SetLineColor(2);
  nmom_woK0_woSid_won_mc->Draw("HEsame");
  cnmom->cd(2);
  nmom_wK0_woSid_won_data->Draw("HE");
  nmom_wK0_woSid_won_mc->SetLineColor(3);
  nmom_wK0_woSid_won_mc->Draw("HEsame");
  nmom_wK0_woSid_won_mcgeta->SetLineColor(2);
  nmom_wK0_woSid_won_mcgeta->Draw("HEsame");
  TH1D* nmom_wK0_woSid_won_mcsum = (TH1D*)nmom_wK0_woSid_won_mc->Clone("nmom_wK0_woSid_won_mc");
  nmom_wK0_woSid_won_mcsum->Add(nmom_wK0_woSid_won_mcgeta,1.0);
  nmom_wK0_woSid_won_mcsum->SetLineColor(6);
  nmom_wK0_woSid_won_mcsum->Draw("HEsame");
  cnmom->cd(3);
  nmom_woK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,2);
  nmom_woK0_woSid_won_ratio->Draw("HE");
  TCanvas *c_nmom_mod = new TCanvas("c_nmom_mod","c_nmom_mod");
  c_nmom_mod->cd();
  nmom_woK0_woSid_won_ratio->Draw("HE");
  TF1 *f_nmom_mod = new TF1("f_nmom_mod",func_nmom_mod,0.13,1.0,12);
  //f_nmom_mod->SetParameter(0,4.30264);
  //f_nmom_mod->SetParameter(1,-16.3496);
  f_nmom_mod->FixParameter(4,0.01);//woods-saxon
  nmom_woK0_woSid_won_ratio->Fit("f_nmom_mod","","",0.140,1.0);

  cnmom->cd(4);
  nmom_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,2);
  nmom_wK0_woSid_won_ratio->Draw("HE");
  //TCanvas *c_nmom_wK0_mod = new TCanvas("c_nmom_wK0_mod","c_nmom_wK0_mod");
  //c_nmom_wK0_mod->cd();
  //nmom_wK0_woSid_won_ratio->Draw("HE");
  //TF1 *f_nmom_wK0_mod = new TF1("f_nmom_wK0_mod",func_nmom_mod,0.13,1.0,12);
  //f_nmom_wK0_mod->SetParameters(param_nmom_wK0_mod);
  //f_nmom_wK0_mod->FixParameter(4,0.01);//woods-saxon
  //nmom_wK0_woSid_won_ratio->Fit("f_nmom_wK0_mod","","",0.14,1.0);

  TCanvas *cnmom_sum = new TCanvas("cnmom_sum","cnmom_sum",1000,1000);
  TH1D* nmom_woSid_won_data = (TH1D*) nmom_woK0_woSid_won_data->Clone("nmom_woSid_won_data");
  nmom_woSid_won_data->Add(nmom_wK0_woSid_won_data);
  nmom_woSid_won_data->Draw("HE");
  TH1D* nmom_woSid_won_pipinxmc = (TH1D*)nmom_woK0_woSid_won_mc->Clone("nmom_woSid_won_pipinxmc");
  nmom_woSid_won_pipinxmc->Add(nmom_wK0_woSid_won_mcgeta);
  nmom_woSid_won_pipinxmc->SetLineColor(2);
  nmom_woSid_won_pipinxmc->Draw("HEsame");
  nmom_wK0_woSid_won_mc->SetLineColor(3);
  nmom_wK0_woSid_won_mc->Draw("HEsame");
  TH1D* nmom_woSid_won_mcsum = (TH1D*)nmom_woSid_won_pipinxmc->Clone("nmom_woSid_won_mcsum");
  nmom_woSid_won_mcsum->Add(nmom_wK0_woSid_won_mc);
  nmom_woSid_won_mcsum->SetLineColor(6);
  nmom_woSid_won_mcsum->Draw("HEsame");
  
  
  
  std::ofstream os;
  os.open("param_corr.txt");
   //os << "IMnpip" << endl;
  //os << "MMnmiss " << endl;
  //os << "q " << endl;
  //for(int i=0;i<f_q_mod->GetNpar();i++){
  //for(int i=0;i<f_nmom_mod->GetNpar();i++){
  //for(int i=0;i<f_q_wK0_mod->GetNpar();i++){
  //for(int i=0;i<f_MMnmiss_wK0_mod->GetNpar();i++){
  //for(int i=0;i<f_MMnmiss_mod->GetNpar();i++){
  //for(int i=0;i<f_IMpippim_mod->GetNpar();i++){
  //for(int i=0;i<f_IMnpip_mod->GetNpar();i++){
  //for(int i=0;i<f_IMnpip_wK0->GetNpar();i++){
  for(int i=0;i<f_IMnpim_mod->GetNpar();i++){
  //for(int i=0;i<f_IMnpim_wK0_mod->GetNpar();i++){
  //for(int i=0;i<f_nmom_wK0_mod->GetNpar();i++){
  //for(int i=0;i<f_Mompippim->GetNpar();i++){
    os << std::setprecision(6);
    //os << f_nmom_mod->GetParameter(i) << ",";
    //os << f_q_wK0_mod->GetParameter(i) << ",";
    //os << f_q_mod->GetParameter(i) << ",";
    //os << f_MMnmiss_wK0_mod->GetParameter(i) << ",";
    //os << f_nmom_wK0_mod->GetParameter(i) << ",";
    //os << f_MMnmiss_mod->GetParameter(i) << ",";
    //os << f_IMpippim_mod->GetParameter(i) << ",";
    //os << f_IMnpip_mod->GetParameter(i) << ",";
    os << f_IMnpim_mod->GetParameter(i) << ",";
    //os << f_IMnpip_wK0->GetParameter(i) << ",";
    //os << f_IMnpim_wK0_mod->GetParameter(i) << ",";
    //os << f_Mompippim->GetParameter(i) << ",";
    //os << f_IMnpip_wK0->GetParameter(i) << ",";
    os << endl;
  }
  // os << f_q_mod->GetParameter(i) << ",";
  //  os << endl;
  //}
  //for(int i=0;i<f_nmom_mod->GetNpar();i++){
  //  os << std::setprecision(6);
  //  os << f_nmom_mod->GetParameter(i) << ",";
  //  os << endl;
  //}

}
