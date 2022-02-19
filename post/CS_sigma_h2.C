#include "anacuts.h"

void CS_sigma_h2()
{
  TH1::SetDefaultSumw2();
  TFile *f = TFile::Open("evanaIMsigma_npi_h2_v9_out_iso_nostop_sub.root");
  //TFile *f = TFile::Open("evanaIMsigma_npi_h2_v7_out_iso_nostop_sub.root");
  
  
  TH2F* Cospicm_IMnpip_pi = (TH2F*)f->Get("Cospicm_IMnpip_pi");
  TH2F* Cospicm_IMnpim_pi = (TH2F*)f->Get("Cospicm_IMnpim_pi");
  
  
  TFile *facc = TFile::Open("accH2.root","READ");
  TH1D* accSp = (TH1D*)facc->Get("accCosSp");
  TH1D* accSm = (TH1D*)facc->Get("accCosSm");
  TH1D* accSp2 = (TH1D*)facc->Get("accCosSp2");
  TH1D* accSm2 = (TH1D*)facc->Get("accCosSm2");
   
  const int Splow = Cospicm_IMnpip_pi->GetXaxis()->FindBin(anacuts::Sigmap_MIN);
  const int Sphigh = Cospicm_IMnpip_pi->GetXaxis()->FindBin(anacuts::Sigmap_MAX);
  TH1D* CS_Sp = (TH1D*)Cospicm_IMnpip_pi->ProjectionY("CS_Sp",Splow,Sphigh);
  const int Smlow = Cospicm_IMnpim_pi->GetXaxis()->FindBin(anacuts::Sigmam_MIN);
  const int Smhigh = Cospicm_IMnpim_pi->GetXaxis()->FindBin(anacuts::Sigmam_MAX);
  TH1D* CS_Sm = (TH1D*)Cospicm_IMnpim_pi->ProjectionY("CS_Sm",Smlow,Smhigh);
  TH1D* CS_Sp2 = (TH1D*)f->Get("Cospicm_pi_Sp");
  TH1D* CS_Sm2 = (TH1D*)f->Get("Cospicm_pi_Sm");
  std::cout << "bin width: " << CS_Sp->GetBinWidth(2) << std::endl;
  CS_Sp->RebinX(5);
  CS_Sm->RebinX(5);
  CS_Sp2->RebinX(5);
  CS_Sm2->RebinX(5);
  CS_Sp2->Print("all");
  CS_Sm2->Print("all");
  CS_Sp->Divide(accSp);
  CS_Sm->Divide(accSm);
  CS_Sp2->Divide(accSp2);
  CS_Sm2->Divide(accSm2);
 
  double Lumi=379.8;
  TCanvas *cCS_Sp = new TCanvas("cCS_Sp","cCS_Sp",1000,800);
  //CS_Sp->RebinX(5);
  double CospiBinW = CS_Sp->GetBinWidth(2);
  double CospiBinW2 = CS_Sp2->GetBinWidth(2);
  std::cout << "Bin width " << CospiBinW2 << std::endl;
  const double cosToStrbin = 2.*3.1415*(CospiBinW); 
  const double cosToStrbin2 = 2.*3.1415*(CospiBinW2); 
  CS_Sp->Scale(1./cosToStrbin/Lumi);
  CS_Sp2->Scale(1./cosToStrbin2/Lumi);
  CS_Sp->GetXaxis()->SetRangeUser(0.3,1);
  //CS_Sp->Draw("E");
  TGraphAsymmErrors* gCS_Sp = new TGraphAsymmErrors(CS_Sp);
  gCS_Sp->GetXaxis()->SetRangeUser(0.3,1);
  //gCS_Sp->Draw("AP");
  TGraphAsymmErrors* gCS_Sp2 = new TGraphAsymmErrors(CS_Sp2);
  gCS_Sp2->GetXaxis()->SetRangeUser(0.3,1);
  gCS_Sp2->Draw("AP*");
  TCanvas *cCS_Sm = new TCanvas("cCS_Sm","cCS_Sm",1000,800);
  //CS_Sm->RebinX(5);
  CS_Sm->Scale(1./cosToStrbin/Lumi);
  CS_Sm2->Scale(1./cosToStrbin2/Lumi);
  CS_Sm->GetXaxis()->SetRangeUser(0.3,1);
  //CS_Sm->Draw("E");
  TGraphAsymmErrors* gCS_Sm = new TGraphAsymmErrors(CS_Sm);
  gCS_Sm->GetXaxis()->SetRangeUser(0.3,1);
  //gCS_Sm->Draw("AP");
  TGraphAsymmErrors* gCS_Sm2 = new TGraphAsymmErrors(CS_Sm2);
  gCS_Sm2->GetXaxis()->SetRangeUser(0.3,1);
  gCS_Sm2->Draw("AP*");


  TFile *fout = new TFile("CSsigma_H2.root","RECREATE");
  fout->cd();
  gCS_Sp->Write();
  gCS_Sm->Write();
  gCS_Sp2->Write();
  gCS_Sm2->Write();

}
