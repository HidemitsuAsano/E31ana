//fit macro to decompose K0/S+/S- final states
#include <TF2.h>
#include <TH2.h>
#include <TMath.h>
#include "anacuts.h"

Double_t K0fit2d(Double_t *x, Double_t *par)
{
  //Double_t r1 = (x[0]-par[1])/par[2];
  //Double_t r2 = (x[1]-par[3])/par[4];
  Double_t r1 = x[0]*par[1];
  Double_t r2 = x[1]*par[2];

  return par[0]*TMath::Exp(r1+r2);
}


void Fit2DK0()
{
  TFile *fr = TFile::Open("evanaIMpisigma_npippim_v202_out_iso_qhi_sub.root","READ");

  auto *IMnpim_IMnpip_dE_wK0_woSid_n = (TH2F*)fr->Get("IMnpim_IMnpip_dE_wK0_woSid_n");
  auto *IMnpim_IMnpip_dE_wK0_woSid_n_45rot = (TH2F*)fr->Get("IMnpim_IMnpip_dE_wK0_woSid_n_45rot");
  auto *cIMnpim_IMnpip_dE_wK0_woSid_n = new TCanvas("cIMnpim_IMnpip_dE_wK0_woSid_n","cIMnpim_IMnpip_dE_wK0_woSid_n",800,800);
  cIMnpim_IMnpip_dE_wK0_woSid_n->Divide(2,2,0.,0.);
  cIMnpim_IMnpip_dE_wK0_woSid_n->cd(3);
  
  IMnpim_IMnpip_dE_wK0_woSid_n->RebinX(4);
  IMnpim_IMnpip_dE_wK0_woSid_n->RebinY(4);
  //IMnpim_IMnpip_dE_wK0_woSid_n->GetXaxis()->SetRangeUser(1.0,1.8);
  //IMnpim_IMnpip_dE_wK0_woSid_n->GetYaxis()->SetRangeUser(1.0,1.8);
  IMnpim_IMnpip_dE_wK0_woSid_n->SetMinimum(0);
  IMnpim_IMnpip_dE_wK0_woSid_n->Draw("colz");
  
  cIMnpim_IMnpip_dE_wK0_woSid_n->cd(1);
  TH1D *IMnpip_wK0_woSid_n = (TH1D*)IMnpim_IMnpip_dE_wK0_woSid_n->ProjectionX("IMnpip_wK0_woSid_n");
  IMnpip_wK0_woSid_n->Draw("HE");
  
  cIMnpim_IMnpip_dE_wK0_woSid_n->cd(4);
  TH1D *IMnpim_wK0_woSid_n = (TH1D*)IMnpim_IMnpip_dE_wK0_woSid_n->ProjectionY("IMnpim_wK0_woSid_n");
  IMnpim_wK0_woSid_n->Draw("HE");
   
  auto *cIMnpim_IMnpip_dE_wK0_woSid_n2 = new TCanvas("cIMnpim_IMnpip_dE_wK0_woSid_n2","cIMnpim_IMnpip_dE_wK0_woSid_n2",800,800);
  cIMnpim_IMnpip_dE_wK0_woSid_n2->Divide(2,2,0.,0.);
  cIMnpim_IMnpip_dE_wK0_woSid_n2->cd(3);
  TH2D* IMnpim_IMnpip_dE_wK0_woSid_n2 = (TH2D*)IMnpim_IMnpip_dE_wK0_woSid_n->Clone("IMnpim_IMnpip_dE_wK0_woSid_n2");
  IMnpim_IMnpip_dE_wK0_woSid_n2->Draw("colz");
  const int nbinx = IMnpim_IMnpip_dE_wK0_woSid_n2->GetXaxis()->FindBin(1.2);
  const int nbiny = IMnpim_IMnpip_dE_wK0_woSid_n2->GetYaxis()->FindBin(1.2);
  cIMnpim_IMnpip_dE_wK0_woSid_n2->cd(1);
  TH1D *IMnpip_wK0_woSid_n2 = (TH1D*)IMnpim_IMnpip_dE_wK0_woSid_n2->ProjectionX("IMnpip_wK0_woSid_n2",nbinx,300);
  IMnpip_wK0_woSid_n2->Draw("E");
  TF1 *fnpip = new TF1("fnpip","expo",1.13,1.26);
  IMnpip_wK0_woSid_n2->Fit(fnpip,"","",1.13,1.26);
  double parnpip[2];
  fnpip->GetParameters(parnpip);
  std::cout << "par0: " << parnpip[0] << std::endl;
  std::cout << "par1: " << parnpip[1] << std::endl;
  
  cIMnpim_IMnpip_dE_wK0_woSid_n2->cd(4);
  TH1D *IMnpim_wK0_woSid_n2 = (TH1D*)IMnpim_IMnpip_dE_wK0_woSid_n2->ProjectionY("IMnpim_wK0_woSid_n2",nbiny,300);
  IMnpim_wK0_woSid_n2->Draw("E");
  TF1 *fnpim = new TF1("fnpim","expo",1.12,1.36);
  IMnpim_wK0_woSid_n2->Fit(fnpim,"","",1.12,1.36);
  double parnpim[2];
  fnpim->GetParameters(parnpim);
  std::cout << "par0: " << parnpim[0] << std::endl;
  std::cout << "par1: " << parnpim[1] << std::endl;
  
  cIMnpim_IMnpip_dE_wK0_woSid_n2->cd(3);
  const Int_t npar = 3;
  //TF2 *f2 = new TF2("f2",K0fit2d,1.13,1.26,1.12,1.36,3);   
  TF2 *f2 = new TF2("f2",K0fit2d,1.,1.,1.5,1.5,3);   
  f2->SetParameters(parnpip[0],parnpip[1],parnpim[1]);
  //IMnpim_IMnpip_dE_wK0_woSid_n2->Fit("f2");
  f2->Draw("cont1 same");
  TH2F* hfit = (TH2F*)f2->GetHistogram();
  
  TCanvas *cfitcomp = new TCanvas("cfitcomp","cfitcomp",800,800);
  cfitcomp->Divide(2,2,0,0);
  cfitcomp->cd(3);
  //IMnpim_IMnpip_dE_wK0_woSid_n2->Draw("colz");
  hfit->Draw("colz");

  cfitcomp->cd(1);
  IMnpip_wK0_woSid_n->Draw("E");
  TH1D* hfitnpip = (TH1D*)hfit->ProjectionX("hfitnpip");
  hfitnpip->SetLineColor(2);
  hfitnpip->Draw("HEsame");

  cfitcomp->cd(4);
  IMnpim_wK0_woSid_n->Draw("E");
  
  //try interpolation
  TCanvas *cinter = new TCanvas("cinter","cinter",800,800);
  cinter->Divide(2,2,0,0);
  cinter->cd(3);
  const int binxx = IMnpim_IMnpip_dE_wK0_woSid_n2->GetXaxis()->FindBin(anacuts::Sigmap_center);
  const int binyy = IMnpim_IMnpip_dE_wK0_woSid_n2->GetYaxis()->FindBin(anacuts::Sigmam_center);
  TH2F *IMnpim_IMnpip_wK0_woSid_n_inter = (TH2F*)IMnpim_IMnpip_dE_wK0_woSid_n2->Clone("IMnpim_IMnpip_wK0_woSid_n_inter");
  for(int ix=0;ix<IMnpim_IMnpip_dE_wK0_woSid_n2->GetNbinsX();ix++){
    double bincx = IMnpim_IMnpip_dE_wK0_woSid_n2->GetXaxis()->GetBinCenter(ix);
    double bincy = IMnpim_IMnpip_dE_wK0_woSid_n2->GetYaxis()->GetBinCenter(binyy);
    double contx = IMnpim_IMnpip_dE_wK0_woSid_n2->Interpolate(bincx,bincy);
    IMnpim_IMnpip_wK0_woSid_n_inter->SetBinContent(ix,binyy,contx);
    //std::cout << ix  << " " << binyy << std::endl;
    //std::cout << bincx << "  " << bincy << std::endl;
    //std::cout << contx << std::endl;
    //std::cout << std::endl;
  }
  //for(int iy=0;iy<IMnpim_IMnpip_dE_wK0_woSid_n2->GetNbinsY();iy++){
  //  double conty = IMnpim_IMnpip_dE_wK0_woSid_n2->Interpolate(binxx,iy);
  //  IMnpim_IMnpip_wK0_woSid_n_inter->SetBinContent(binxx,iy,conty);
  //}
  IMnpim_IMnpip_wK0_woSid_n_inter->Draw("colz");
 

  //TF1Convolution *f_conv = new TF1Convolution("expo","gaus",1.0,1.8,true);
  //f_conv->SetRange(1.0,1.8);
  //f_conv->SetNofPointsFFT(1000);
  //TF1 *f1 = new TF1("f1",*f_conv,1.0,1.8,f_conv->GetNpar());
  //f1->SetParameters(100,0,1.0,-5.0);
  //f1->SetParLimits(3,-10.0,0.0);
  //IMnpim_wK0_woSid_n->Fit("f1","","",1.1,1.3);

  auto *cIMnpim_IMnpip_dE_wK0_woSid_n_45rot = new TCanvas("cIMnpim_IMnpip_dE_wK0_woSid_n_45rot","cIMnpim_IMnpip_dE_wK0_woSid_n_45rot",800,800);
  cIMnpim_IMnpip_dE_wK0_woSid_n_45rot->Divide(2,2,0,0);
  cIMnpim_IMnpip_dE_wK0_woSid_n_45rot->cd(3);
  IMnpim_IMnpip_dE_wK0_woSid_n_45rot->SetMinimum(0);
  IMnpim_IMnpip_dE_wK0_woSid_n_45rot->Draw("colz");
  
  cIMnpim_IMnpip_dE_wK0_woSid_n_45rot->cd(1);
  IMnpim_IMnpip_dE_wK0_woSid_n_45rot->ProjectionX()->Draw();
  
  cIMnpim_IMnpip_dE_wK0_woSid_n_45rot->cd(4);
  IMnpim_IMnpip_dE_wK0_woSid_n_45rot->ProjectionY()->Draw();

  /*
  auto *cIMnpim_IMnpip_dE_woK0_wSid_n_woSm = new TCanvas("cIMnpim_IMnpip_dE_woK0_wSid_n_woSm","cIMnpim_IMnpip_dE_woK0_wSid_n_woSm",800,400);
  cIMnpim_IMnpip_dE_woK0_wSid_n_woSm->Divide(2,1);
  cIMnpim_IMnpip_dE_woK0_wSid_n_woSm->cd(1);
  TH2D* IMnpim_IMnpip_dE_woK0_wSid_n_woSm = (TH2D*)fr->Get("IMnpim_IMnpip_dE_woK0_wSid_n_woSm");
  IMnpim_IMnpip_dE_woK0_wSid_n_woSm->SetMinimum(0);
  IMnpim_IMnpip_dE_woK0_wSid_n_woSm->RebinX(5);
  IMnpim_IMnpip_dE_woK0_wSid_n_woSm->RebinY(5);
  IMnpim_IMnpip_dE_woK0_wSid_n_woSm->GetXaxis()->SetRangeUser(1.0,1.8);
  IMnpim_IMnpip_dE_woK0_wSid_n_woSm->GetYaxis()->SetRangeUser(1.0,1.8);
  IMnpim_IMnpip_dE_woK0_wSid_n_woSm->Draw("colz");
  cIMnpim_IMnpip_dE_woK0_wSid_n_woSm->cd(2);
  TH1D* IMnpim_dE_woK0_wSid_n_woSm = (TH1D*)IMnpim_IMnpip_dE_woK0_wSid_n_woSm->ProjectionY("IMnpim_dE_woK0_wSid_n_woSm");
  IMnpim_dE_woK0_wSid_n_woSm->Draw("HE");

  auto *cIMnpim_IMnpip_dE_woK0_wSid_n_woSp = new TCanvas("cIMnpim_IMnpip_dE_woK0_wSid_n_woSp","cIMnpim_IMnpip_dE_woK0_wSid_n_woSp",800,400);
  cIMnpim_IMnpip_dE_woK0_wSid_n_woSp->Divide(2,1);
  cIMnpim_IMnpip_dE_woK0_wSid_n_woSp->cd(1);
  TH2D* IMnpim_IMnpip_dE_woK0_wSid_n_woSp = (TH2D*)fr->Get("IMnpim_IMnpip_dE_woK0_wSid_n_woSp");
  IMnpim_IMnpip_dE_woK0_wSid_n_woSp->SetMinimum(0);
  IMnpim_IMnpip_dE_woK0_wSid_n_woSp->RebinX(5);
  IMnpim_IMnpip_dE_woK0_wSid_n_woSp->RebinY(5);
  IMnpim_IMnpip_dE_woK0_wSid_n_woSp->GetXaxis()->SetRangeUser(1.0,1.8);
  IMnpim_IMnpip_dE_woK0_wSid_n_woSp->GetYaxis()->SetRangeUser(1.0,1.8);
  IMnpim_IMnpip_dE_woK0_wSid_n_woSp->Draw("colz");
  cIMnpim_IMnpip_dE_woK0_wSid_n_woSp->cd(2);
  TH1D* IMnpip_dE_woK0_wSid_n_woSp = (TH1D*)IMnpim_IMnpip_dE_woK0_wSid_n_woSp->ProjectionX("IMnpip_dE_woK0_wSid_n_woSp");
  IMnpip_dE_woK0_wSid_n_woSp->Draw("HE");
  */

}
