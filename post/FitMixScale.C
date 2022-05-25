//display vicinity of signal region of data and mixed events
//determine scaling factor of mixed events by fitting
#include "anacuts.h"

TH1D* MMnmiss_woK0_woSm_mix;
TH1D* IMnpip_woK0_woSm_mix;
TH1D* MMnmiss_woK0_woSp_mix;
TH1D* IMnpim_woK0_woSp_mix;

double sysScale=0.05;
const int dEcut=2;
const int Version=241;
const bool Scale=false;
const int qcut = 0;

void FitMixScale()
{
  gStyle->SetOptStat(0);
  fr = TFile::Open(Form("evanaIMpisigma_npippim_v%d_out_dE%d_iso_nostop.root",Version,dEcut));
  fmix = TFile::Open(Form("evanaIMpisigma_npippim_v%d_MIX_cut4_out_dE%d_iso_nostop_sys0.root",Version,dEcut));
  fr->Print();
  fmix->Print();

  TH2D* MMnmiss_IMnpip_woK0_woSm_vici_data = (TH2D*)fr->Get("MMnmiss_IMnpip_dE_woK0_woSm_vici");
  TH2D* MMnmiss_IMnpim_woK0_woSp_vici_data = (TH2D*)fr->Get("MMnmiss_IMnpim_dE_woK0_woSp_vici");

  TH2D* MMnmiss_IMnpip_woK0_woSm_vici_mix = (TH2D*)fmix->Get("MMnmiss_IMnpip_dE_woK0_woSm_vici");
  TH2D* MMnmiss_IMnpim_woK0_woSp_vici_mix = (TH2D*)fmix->Get("MMnmiss_IMnpim_dE_woK0_woSp_vici");

  MMnmiss_IMnpip_woK0_woSm_vici_data->GetXaxis()->SetRangeUser(1.1,1.3);
  MMnmiss_IMnpip_woK0_woSm_vici_data->GetYaxis()->SetRangeUser(0.7,1.1);
  MMnmiss_IMnpip_woK0_woSm_vici_mix->GetXaxis()->SetRangeUser(1.1,1.3);
  MMnmiss_IMnpip_woK0_woSm_vici_mix->GetYaxis()->SetRangeUser(0.7,1.1);
  MMnmiss_IMnpim_woK0_woSp_vici_data->GetXaxis()->SetRangeUser(1.1,1.3);
  MMnmiss_IMnpim_woK0_woSp_vici_data->GetYaxis()->SetRangeUser(0.7,1.1);
  MMnmiss_IMnpim_woK0_woSp_vici_mix->GetXaxis()->SetRangeUser(1.1,1.3);
  MMnmiss_IMnpim_woK0_woSp_vici_mix->GetYaxis()->SetRangeUser(0.7,1.1);
  
  TCanvas *c1 = new TCanvas("c1","c1");
  MMnmiss_IMnpip_woK0_woSm_vici_data->Draw("colz");

  TCanvas *c2 = new TCanvas("c2","c2");
  MMnmiss_IMnpim_woK0_woSp_vici_data->Draw("colz");

  TCanvas *c3 = new TCanvas("c3","c3");
  MMnmiss_IMnpip_woK0_woSm_vici_mix->Draw("colz");

  TCanvas *c4 = new TCanvas("c4","c4");
  MMnmiss_IMnpim_woK0_woSp_vici_mix->Draw("colz");
  
  TCanvas *c5 = new TCanvas("c5","c5");
  TH1D* IMnpip_woK0_woSm_data = (TH1D*)MMnmiss_IMnpip_woK0_woSm_vici_data->ProjectionX("IMnpip_woK0_woSm_vici_data");
  IMnpip_woK0_woSm_mix  = (TH1D*)MMnmiss_IMnpip_woK0_woSm_vici_mix->ProjectionX("IMnpip_woK0_woSm_vici_mix");
  
  IMnpip_woK0_woSm_data->GetXaxis()->SetRangeUser(1.14,1.24);
  IMnpip_woK0_woSm_data->Draw();
  TF1 *fIMnpip = new TF1("fIMnpip",fit_IMnpip,1.14,1.24,1);
  TF1 *fIMnpip_up = new TF1("fIMnpip_up",fit_IMnpip,1.14,1.24,1);
  TF1 *fIMnpip_down = new TF1("fIMnpip_down",fit_IMnpip,1.14,1.24,1);
  fIMnpip->SetParameter(1,1.0);
  IMnpip_woK0_woSm_data->Fit("fIMnpip","r","",1.212,1.218);
  Double_t ret = fIMnpip->GetParameter(0);
  Double_t reterr = fIMnpip->GetParError(0);
  fIMnpip->SetParameter(0,ret);
  fIMnpip->SetLineColor(3);
  fIMnpip->Draw("same");
   

  TCanvas *c6 = new TCanvas("c6","c6");
  TH1D* MMnmiss_woK0_woSm_data = (TH1D*)MMnmiss_IMnpip_woK0_woSm_vici_data->ProjectionY("MMnmiss_woK0_woSm_vici_data");
  MMnmiss_woK0_woSm_mix = (TH1D*)MMnmiss_IMnpip_woK0_woSm_vici_mix->ProjectionY("MMnmiss_woK0_woSm_vici_mix");
  
  MMnmiss_woK0_woSm_data->Draw();
  TF1 *fMMnmiss_Sp = new TF1("fMMnmiss_Sp",fit_MMnmiss_Sp,0.7,1.1,1);
  TF1 *fMMnmiss_Sp_up = new TF1("fMMnmiss_Sp_up",fit_MMnmiss_Sp,0.7,1.1,1);
  TF1 *fMMnmiss_Sp_down = new TF1("fMMnmiss_Sp_down",fit_MMnmiss_Sp,0.7,1.1,1);
  fMMnmiss_Sp->SetParameter(1,1.0);
  //MMnmiss_woK0_woSm_data->Fit("fMMnmiss_Sp","r","",1.0,1.05);
  MMnmiss_woK0_woSm_data->Fit("fMMnmiss_Sp","r","",0.75,0.85);
  Double_t ret2 = fMMnmiss_Sp->GetParameter(0);
  Double_t ret2err = fMMnmiss_Sp->GetParError(0);
  fMMnmiss_Sp->SetParameter(0,ret2);
  fMMnmiss_Sp->SetLineColor(3);
  fMMnmiss_Sp->Draw("same");
  
  TCanvas *c7 = new TCanvas("c7","c7");
  TH1D* IMnpim_woK0_woSp_data = (TH1D*)MMnmiss_IMnpim_woK0_woSp_vici_data->ProjectionX("IMnpim_woK0_woSp_vici_data");
  IMnpim_woK0_woSp_mix  = (TH1D*)MMnmiss_IMnpim_woK0_woSp_vici_mix->ProjectionX("IMnpim_woK0_woSp_vici_mix");
  
  IMnpim_woK0_woSp_data->GetXaxis()->SetRangeUser(1.14,1.24);
  IMnpim_woK0_woSp_data->Draw();
  TF1 *fIMnpim = new TF1("fIMnpim",fit_IMnpim,1.14,1.24,1);
  TF1 *fIMnpim_up = new TF1("f1Mnpim_up",fit_IMnpim,1.14,1.24,1);
  TF1 *fIMnpim_down = new TF1("fIMnpim_down",fit_IMnpim,1.14,1.24,1);
  fIMnpim->SetParameter(1,1.0);
  IMnpim_woK0_woSp_data->Fit("fIMnpim","r","",1.217,1.226);
  Double_t ret3 = fIMnpim->GetParameter(0);
  Double_t ret3err = fIMnpim->GetParError(0);
  fIMnpim->SetParameter(0,ret3);
  fIMnpim->SetLineColor(3);
  fIMnpim->Draw("same");

  TCanvas *c8 = new TCanvas("c8","c8");
  TH1D* MMnmiss_woK0_woSp_data = (TH1D*)MMnmiss_IMnpim_woK0_woSp_vici_data->ProjectionY("MMnmiss_woK0_woSp_vici_data");
  MMnmiss_woK0_woSp_mix = (TH1D*)MMnmiss_IMnpim_woK0_woSp_vici_mix->ProjectionY("MMnmiss_woK0_woSp_vici_mix");
  
  MMnmiss_woK0_woSp_data->Draw();
  TF1 *fMMnmiss_Sm = new TF1("fMMnmiss_Sm",fit_MMnmiss_Sm,0.7,1.1,1);
  TF1 *fMMnmiss_Sm_up = new TF1("fMMnmiss_Sm_up",fit_MMnmiss_Sm,0.7,1.1,1);
  TF1 *fMMnmiss_Sm_down = new TF1("fMMnmiss_Sm_down",fit_MMnmiss_Sm,0.7,1.1,1);
  fMMnmiss_Sm->SetParameter(1,1.0);
  //MMnmiss_woK0_woSp_data->Fit("fMMnmiss_Sm","r","",1.0,1.05);
  MMnmiss_woK0_woSp_data->Fit("fMMnmiss_Sm","r","",0.75,0.85);
  Double_t ret4 = fMMnmiss_Sm->GetParameter(0);
  Double_t ret4err = fMMnmiss_Sm->GetParError(0);
  fMMnmiss_Sm->SetParameter(0,ret4);
  fMMnmiss_Sm->SetLineColor(3);
  fMMnmiss_Sm->Draw("same");

  
  double avgSp = (ret*ret2err+ret2*reterr)/(reterr+ret2err);
  double avgSm = (ret3*ret4err+ret4*ret3err)/(ret3err+ret4err);
  double avg = (avgSp+avgSm)/2.0;
  std::cout << "Sp mode scale avgSp: " << avgSp << std::endl;
  fIMnpip->SetParameter(0,avg);
  fIMnpip_up->SetParameter(0,avg*(1.0+sysScale));
  fIMnpip_down->SetParameter(0,avg*(1.0-sysScale));
  fIMnpip->SetLineColor(3);
  fIMnpip_up->SetLineColor(4);
  fIMnpip_down->SetLineColor(4);
  c5->cd();
  fIMnpip->Draw("same");
  fIMnpip_up->Draw("same");
  fIMnpip_down->Draw("same");
  fMMnmiss_Sp->SetParameter(0,avg);
  fMMnmiss_Sp_up->SetParameter(0,avg*(1.0+sysScale));
  fMMnmiss_Sp_down->SetParameter(0,avg*(1.0-sysScale));
  fMMnmiss_Sp->SetLineColor(3);
  fMMnmiss_Sp_up->SetLineColor(4);
  fMMnmiss_Sp_down->SetLineColor(4);
  c6->cd();
  fMMnmiss_Sp->Draw("same");
  fMMnmiss_Sp_up->Draw("same");
  fMMnmiss_Sp_down->Draw("same");
 

  //double avgSm = ret4;
  std::cout << "Sm mode scale avg: " << avgSm << std::endl;
  fIMnpim->SetParameter(0,avg);
  fIMnpim_up->SetParameter(0,avg*(1.0+sysScale));
  fIMnpim_down->SetParameter(0,avg*(1.0-sysScale));
  fIMnpim->SetLineColor(3);
  fIMnpim_up->SetLineColor(4);
  fIMnpim_down->SetLineColor(4);
  c7->cd();
  fIMnpim->Draw("same");
  fIMnpim_up->Draw("same");
  fIMnpim_down->Draw("same");
  fMMnmiss_Sm->SetParameter(0,avg);
  fMMnmiss_Sm_up->SetParameter(0,avg*(1.0+sysScale));
  fMMnmiss_Sm_down->SetParameter(0,avg*(1.0-sysScale));
  fMMnmiss_Sm->SetLineColor(3);
  fMMnmiss_Sm_up->SetLineColor(4);
  fMMnmiss_Sm_down->SetLineColor(4);
  c8->cd();
  fMMnmiss_Sm->Draw("same");
  fMMnmiss_Sm_up->Draw("same");
  fMMnmiss_Sm_down->Draw("same");
 
  std::cout << "average of Sp and Sm " << avg << std::endl;

  TH2D* MMnmiss_IMnpip_woK0_woSm_data = (TH2D*)fr->Get("MMnmiss_IMnpip_dE_woK0_woSm");
  TH2D* MMnmiss_IMnpim_woK0_woSp_data = (TH2D*)fr->Get("MMnmiss_IMnpim_dE_woK0_woSp");
  TH2D* MMnmiss_IMnpip_woK0_woSm_mix = (TH2D*)fmix->Get("MMnmiss_IMnpip_dE_woK0_woSm");
  TH2D* MMnmiss_IMnpim_woK0_woSp_mix = (TH2D*)fmix->Get("MMnmiss_IMnpim_dE_woK0_woSp");
  TH2D* MMnmiss_IMpippim_woSid_data = (TH2D*)fr->Get("MMnmiss_IMpippim_dE_woSid");
  TH2D* MMnmiss_IMpippim_woSid_mix = (TH2D*)fmix->Get("MMnmiss_IMpippim_dE_woSid");
  TH2D* q_IMnpipi_wSid_n_data = (TH2D*)fr->Get("q_IMnpipi_wSid_n");
  TH2D* q_IMnpipi_wSid_n_mix = (TH2D*)fmix->Get("q_IMnpipi_wSid_n");
   
  if(Scale){
    std::cout << "Scaling mixed events ..... " << std::endl;
    MMnmiss_IMnpip_woK0_woSm_mix->Scale(avg);
    MMnmiss_IMnpim_woK0_woSp_mix->Scale(avg);
    MMnmiss_IMpippim_woSid_mix->Scale(avg);
    q_IMnpipi_wSid_n_mix->Scale(avg); 
  }else{
    std::cout << "mixed events are not scaled " << std::endl;
  }

  TCanvas *c8_1 = new TCanvas("c8_1","c8_1");
  TH2D* MMnmiss_IMnpip_woK0_woSm_sub = (TH2D*)MMnmiss_IMnpip_woK0_woSm_data->Clone("MMnmiss_IMnpip_woK0_woSm_sub");
  MMnmiss_IMnpip_woK0_woSm_sub->Add(MMnmiss_IMnpip_woK0_woSm_mix,-1.0);
  MMnmiss_IMnpip_woK0_woSm_sub->Rebin2D(4,2);
  MMnmiss_IMnpip_woK0_woSm_sub->GetXaxis()->SetRangeUser(1.0,1.7);
  MMnmiss_IMnpip_woK0_woSm_sub->Draw("colz");
  
  TCanvas *c8_2 = new TCanvas("c8_2","c8_2");
  TH2D* MMnmiss_IMnpim_woK0_woSp_sub = (TH2D*)MMnmiss_IMnpim_woK0_woSp_data->Clone("MMnmiss_IMnpim_woK0_woSp_sub");
  MMnmiss_IMnpim_woK0_woSp_sub->Add(MMnmiss_IMnpim_woK0_woSp_mix,-1.0);
  MMnmiss_IMnpim_woK0_woSp_sub->Rebin2D(4,2);
  MMnmiss_IMnpim_woK0_woSp_sub->GetXaxis()->SetRangeUser(1.0,1.7);
  MMnmiss_IMnpim_woK0_woSp_sub->Draw("colz");
  
  TCanvas *c8_3 = new TCanvas("c8_3","c8_3");
  TH2D* MMnmiss_IMpippim_woSid_sub = (TH2D*)MMnmiss_IMpippim_woSid_data->Clone("MMnmiss_IMpippim_woSid_sub");
  MMnmiss_IMpippim_woSid_sub->Add(MMnmiss_IMpippim_woSid_mix,-1.0);
  MMnmiss_IMpippim_woSid_sub->Rebin2D(4,2);
  MMnmiss_IMpippim_woSid_sub->Draw("colz");

  int Spbinlow = MMnmiss_IMnpip_woK0_woSm_data->GetXaxis()->FindBin(anacuts::Sigmap_MIN);
  int Spbinhi = MMnmiss_IMnpip_woK0_woSm_data->GetXaxis()->FindBin(anacuts::Sigmap_MAX);
  int Smbinlow = MMnmiss_IMnpim_woK0_woSp_data->GetXaxis()->FindBin(anacuts::Sigmam_MIN);
  int Smbinhi = MMnmiss_IMnpim_woK0_woSp_data->GetXaxis()->FindBin(anacuts::Sigmam_MAX);
  int nlow = MMnmiss_IMnpip_woK0_woSm_data->GetYaxis()->FindBin(anacuts::neutron_MIN);
  int nhi = MMnmiss_IMnpip_woK0_woSm_data->GetYaxis()->FindBin(anacuts::neutron_MAX);
  int q350 = q_IMnpipi_wSid_n_data->GetYaxis()->FindBin(0.35);
  int q650 = q_IMnpipi_wSid_n_data->GetYaxis()->FindBin(0.65);
  TH1D* MMnmiss_Sp_data = (TH1D*)MMnmiss_IMnpip_woK0_woSm_data->ProjectionY("MMnmiss_Sp_data",Spbinlow,Spbinhi);
  TH1D* MMnmiss_Sp_mix[3];
  for(int isys=0;isys<3;isys++){  
    MMnmiss_Sp_mix[isys] = (TH1D*)MMnmiss_IMnpip_woK0_woSm_mix->ProjectionY(Form("MMnmiss_Sp_mix_%d",isys),Spbinlow,Spbinhi);
    MMnmiss_Sp_mix[isys]->Scale(1+(isys-1)*sysScale);
  }
  TH1D* IMnpip_n_data = (TH1D*)MMnmiss_IMnpip_woK0_woSm_data->ProjectionX("IMnpip_n_data",nlow,nhi);
  TH1D* IMnpip_n_mix[3];
  for(int isys=0;isys<3;isys++){  
    IMnpip_n_mix[isys] = (TH1D*)MMnmiss_IMnpip_woK0_woSm_mix->ProjectionX(Form("IMnpip_n_mix_%d",isys),nlow,nhi);
    IMnpip_n_mix[isys]->Scale(1+(isys-1)*sysScale);
  }
  TH1D* MMnmiss_Sm_data = (TH1D*)MMnmiss_IMnpim_woK0_woSp_data->ProjectionY("MMnmiss_Sm_data",Smbinlow,Smbinhi);
  TH1D* MMnmiss_Sm_mix[3];
  for(int isys=0;isys<3;isys++){
    MMnmiss_Sm_mix[isys] = (TH1D*)MMnmiss_IMnpim_woK0_woSp_mix->ProjectionY(Form("MMnmiss_Sm_mix_%d",isys),Smbinlow,Smbinhi);
    MMnmiss_Sm_mix[isys]->Scale(1+(isys-1)*sysScale);
  }
  TH1D* IMnpim_n_data = (TH1D*)MMnmiss_IMnpim_woK0_woSp_data->ProjectionX("IMnpim_n_data",nlow,nhi);
  TH1D* IMnpim_n_mix[3]; 
  for(int isys=0;isys<3;isys++){
    IMnpim_n_mix[isys] = (TH1D*)MMnmiss_IMnpim_woK0_woSp_mix->ProjectionX(Form("IMnpim_n_mix_%d",isys),nlow,nhi);
    IMnpim_n_mix[isys]->Scale(1+(isys-1)*sysScale);
  }
  
  TH1D* IMpippim_woSid_data = (TH1D*)MMnmiss_IMpippim_woSid_data->ProjectionX("IMpippim_data",nlow,nhi);
  TH1D* IMpippim_woSid_mix[3];
  for(int isys=0;isys<3;isys++){
    IMpippim_woSid_mix[isys] = (TH1D*)MMnmiss_IMpippim_woSid_mix->ProjectionX(Form("IMpippim_mix_%d",isys),nlow,nhi);
    IMpippim_woSid_mix[isys]->Scale(1+(isys-1)*sysScale);
  }
  TH1D* IMnpipi_wSid_n_data_qhi = (TH1D*)q_IMnpipi_wSid_n_data->ProjectionX("IMnpipi_wSid_n_data_qhi",q350,q650);
  TH1D* IMnpipi_wSid_n_mix_qhi[3]; 
  for(int isys=0;isys<3;isys++){
    IMnpipi_wSid_n_mix_qhi[isys] = (TH1D*)q_IMnpipi_wSid_n_mix->ProjectionX(Form("IMnpipi_wSid_n_mix_qhi_%d",isys),q350,q650);   
    IMnpipi_wSid_n_mix_qhi[isys]->Scale(1+(isys-1)*sysScale);
  }
  TH1D* IMnpipi_wSid_n_data_qlo = (TH1D*)q_IMnpipi_wSid_n_data->ProjectionX("IMnpipi_wSid_n_data_qlo",0,q350);
  TH1D* IMnpipi_wSid_n_mix_qlo[3];
  for(int isys=0;isys<3;isys++){
    IMnpipi_wSid_n_mix_qlo[isys] = (TH1D*)q_IMnpipi_wSid_n_mix->ProjectionX(Form("IMnpipi_wSid_n_mix_qlo_%d",isys),0,q350);
    IMnpipi_wSid_n_mix_qlo[isys]->Scale(1+(isys-1)*sysScale);
  }
  

  TCanvas *c9 = new TCanvas("c9","c9");
  MMnmiss_Sp_data->Draw();
  MMnmiss_Sp_mix[1]->SetLineColor(2);
  MMnmiss_Sp_mix[0]->SetLineColor(4);
  MMnmiss_Sp_mix[2]->SetLineColor(4);
  MMnmiss_Sp_mix[1]->Draw("same");
  MMnmiss_Sp_mix[0]->Draw("same");
  MMnmiss_Sp_mix[2]->Draw("same");

  TCanvas *c9_1 = new TCanvas("c9_1","c9_1");
  TH1D* MMnmiss_Sp_sub = (TH1D*)MMnmiss_Sp_data->Clone("MMnmiss_Sp_sub");
  MMnmiss_Sp_sub->Add(MMnmiss_Sp_mix[1],-1);
  MMnmiss_Sp_sub->Draw("HE");
  TF1* fMMSp = new TF1("fMMSp","gaus",0,1.5);
  MMnmiss_Sp_sub->Fit("fMMSp","","",0.85,1.0);
 
  TCanvas *c10 = new TCanvas("c10","c10");
  IMnpip_n_data->RebinX(2);
  IMnpip_n_data->Draw();
  IMnpip_n_mix[1]->SetLineColor(2);
  IMnpip_n_mix[0]->SetLineColor(4);
  IMnpip_n_mix[2]->SetLineColor(4);
  for(int isys=0;isys<3;isys++){
    IMnpip_n_mix[isys]->RebinX(2);
    IMnpip_n_mix[isys]->Draw("same");
  }
  
  TCanvas *c10_1 = new TCanvas("c10_1","c10_1");
  TH1D* IMnpip_n_sub = (TH1D*)IMnpip_n_data->Clone("IMnpip_n_sub");
  IMnpip_n_sub->Add(IMnpip_n_mix[1],-1);
  IMnpip_n_sub->GetXaxis()->SetRangeUser(1.1,1.3);
  IMnpip_n_sub->Draw("HE");
  TF1* fIMSp = new TF1("fIMSp","gaus");
  IMnpip_n_sub->Fit("fIMSp","","",1.16,1.22);


  TCanvas *c11 = new TCanvas("c11","c11");
  MMnmiss_Sm_data->Draw();
  MMnmiss_Sm_mix[1]->SetLineColor(2);
  MMnmiss_Sm_mix[0]->SetLineColor(4);
  MMnmiss_Sm_mix[2]->SetLineColor(4);
  for(int isys=0;isys<3;isys++){
    MMnmiss_Sm_mix[isys]->Draw("same");
  }
  TCanvas *c11_1 = new TCanvas("c11_1","c11_1");
  TH1D* MMnmiss_Sm_sub = (TH1D*)MMnmiss_Sm_data->Clone("MMnmiss_Sm_sub");
  MMnmiss_Sm_sub->Add(MMnmiss_Sm_mix[1],-1);
  MMnmiss_Sm_sub->Draw("HE");
  TF1* fMMSm = new TF1("fMMSm","gaus");
  MMnmiss_Sm_sub->Fit("fMMSm","","",0.85,1.0);

  TCanvas *c12 = new TCanvas("c12","c12");
  IMnpim_n_data->RebinX(2);
  IMnpim_n_data->Draw();
  IMnpim_n_mix[1]->SetLineColor(2);
  IMnpim_n_mix[0]->SetLineColor(4);
  IMnpim_n_mix[2]->SetLineColor(4);
  for(int isys=0;isys<3;isys++){
    IMnpim_n_mix[isys]->RebinX(2);
    IMnpim_n_mix[isys]->Draw("same");
  }
  TCanvas *c12_1 = new TCanvas("c12_1","c12_1");
  TH1D* IMnpim_n_sub = (TH1D*)IMnpim_n_data->Clone("IMnpim_n_sub");
  IMnpim_n_sub->Add(IMnpim_n_mix[1],-1);
  IMnpim_n_sub->GetXaxis()->SetRangeUser(1.1,1.3);
  IMnpim_n_sub->Draw("HE");
  TF1* fIMSm = new TF1("fIMSm","gaus");
  IMnpim_n_sub->Fit("fIMSm","","",1.17,1.22);
  
  TCanvas *c13 = new TCanvas("c13","c13");
  IMpippim_woSid_data->RebinX(2);
  IMpippim_woSid_data->Draw();
  IMpippim_woSid_mix[1]->SetLineColor(2);
  IMpippim_woSid_mix[0]->SetLineColor(4);
  IMpippim_woSid_mix[2]->SetLineColor(4);
  for(int isys=0;isys<3;isys++){
    IMpippim_woSid_mix[isys]->RebinX(2);
    IMpippim_woSid_mix[isys]->Draw("same");
  }

  TCanvas *c14 = new TCanvas("c14","c14");
  IMnpipi_wSid_n_data_qhi->RebinX(3);
  IMnpipi_wSid_n_data_qhi->Draw();
  IMnpipi_wSid_n_mix_qhi[1]->SetLineColor(2);
  IMnpipi_wSid_n_mix_qhi[0]->SetLineColor(4);
  IMnpipi_wSid_n_mix_qhi[2]->SetLineColor(4);
  for(int isys=0;isys<3;isys++){
    IMnpipi_wSid_n_mix_qhi[isys]->RebinX(3);
    IMnpipi_wSid_n_mix_qhi[isys]->Draw("same");
  }

  TCanvas *c15 = new TCanvas("c15","c15");
  IMnpipi_wSid_n_data_qlo->RebinX(3);
  IMnpipi_wSid_n_data_qlo->Draw();
  IMnpipi_wSid_n_mix_qlo[1]->SetLineColor(2);
  IMnpipi_wSid_n_mix_qlo[0]->SetLineColor(4);
  IMnpipi_wSid_n_mix_qlo[2]->SetLineColor(4);
  for(int isys=0;isys<3;isys++){
    IMnpipi_wSid_n_mix_qlo[isys]->RebinX(3);
    IMnpipi_wSid_n_mix_qlo[isys]->Draw("same");
  }

  TCanvas *c16 = new TCanvas("c16","c16");
  q_IMnpipi_wSid_n_data->Draw("colz");
  
  TCanvas *c17 = new TCanvas("c17","c17");
  q_IMnpipi_wSid_n_mix->Draw("colz");

  TCanvas *c18 = new TCanvas("c18","c18");
  TH2D* MMnmiss_IMnpipi_wSid_data = (TH2D*)fr->Get("MMnmiss_IMnpipi_wSid");
  TH2D* MMnmiss_IMnpipi_wSid_mix = (TH2D*)fmix->Get("MMnmiss_IMnpipi_wSid");
  TH1D* MMnmiss_wSid_data = (TH1D*)MMnmiss_IMnpipi_wSid_data->ProjectionY("MMnmiss_wSid_data");
  TH1D* MMnmiss_wSid_mix[3];
  for(int isys=0;isys<3;isys++){
    MMnmiss_wSid_mix[isys] = (TH1D*)MMnmiss_IMnpipi_wSid_mix->ProjectionY(Form("MMnmiss_wSid_mix_%d",isys));
  }
  MMnmiss_wSid_data->Draw("E");
  MMnmiss_wSid_mix[1]->SetLineColor(2);
  MMnmiss_wSid_mix[0]->SetLineColor(4);
  MMnmiss_wSid_mix[2]->SetLineColor(4);
  for(int isys=0;isys<3;isys++){
    MMnmiss_wSid_mix[isys]->Scale(1+(isys-1)*sysScale);
    MMnmiss_wSid_mix[isys]->Draw("same");
  }
  TCanvas *c19 = new TCanvas("c19","c19");
  TH1D* MMnmiss_wSid_sub = (TH1D*)MMnmiss_wSid_data->Clone("MMnmiss_wSid_sub");
  MMnmiss_wSid_sub->Add(MMnmiss_wSid_mix[1],-1.0);
  MMnmiss_wSid_sub->Draw();

  TF1* fMM = new TF1("fMM","gaus",0.85,1.05);
  MMnmiss_wSid_sub->Fit("fMM","","",0.85,1.05);
  
  TF1* fMMLam = new TF1("fMMLam","gaus",1.05,1.15);
  MMnmiss_wSid_sub->Fit("fMMLam","","",1.05,1.15);

  MMnmiss_wSid_sub->Draw();
  fMM->Draw("same");
  fMMLam->Draw("same");

  TCanvas *c = NULL;
  TSeqCollection *SCol = gROOT->GetListOfCanvases();
  int size = SCol->GetSize();
  TIter next(SCol);
  TString pdfname = Form("fitmixscale.pdf");
  for(int i=0;i<size;i++){
    //pdf->NewPage();
    c= (TCanvas*)next();
    c->Modified();
    c->Update();
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


  


}

Double_t fit_IMnpip(Double_t *x,Double_t *par)
{
  if(!IMnpip_woK0_woSm_mix){
    std::cout << "cannot find template !! " << std::endl;
    return 0;
  }
  Double_t IMnpip=x[0]; 
  Int_t bin = IMnpip_woK0_woSm_mix->GetXaxis()->FindBin(IMnpip);
  Double_t br= par[0]*(IMnpip_woK0_woSm_mix->GetBinContent(bin));
   
  return br;
}

Double_t fit_MMnmiss_Sp(Double_t *x,Double_t *par)
{
  if(!MMnmiss_woK0_woSm_mix){
    std::cout << "cannot find template !! " << std::endl;
    return 0;
  }
  Double_t MM=x[0]; 
  Int_t bin = MMnmiss_woK0_woSm_mix->GetXaxis()->FindBin(MM);
  Double_t br= par[0]*(MMnmiss_woK0_woSm_mix->GetBinContent(bin));
   
  return br;
}

Double_t fit_IMnpim(Double_t *x,Double_t *par)
{
  if(!IMnpim_woK0_woSp_mix){
    std::cout << "cannot find template !! " << std::endl;
    return 0;
  }
  Double_t IMnpim=x[0]; 
  Int_t bin = IMnpim_woK0_woSp_mix->GetXaxis()->FindBin(IMnpim);
  Double_t br= par[0]*(IMnpim_woK0_woSp_mix->GetBinContent(bin));
   
  return br;
}

Double_t fit_MMnmiss_Sm(Double_t *x,Double_t *par)
{
  if(!MMnmiss_woK0_woSp_mix){
    std::cout << "cannot find template !! " << std::endl;
    return 0;
  }
  Double_t MM=x[0]; 
  Int_t bin = MMnmiss_woK0_woSp_mix->GetXaxis()->FindBin(MM);
  Double_t br= par[0]*(MMnmiss_woK0_woSp_mix->GetBinContent(bin));
   
  return br;
}


