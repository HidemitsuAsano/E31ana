//fit macro to decompose K0/S+/S- final states
#include <TF2.h>
#include <TH2.h>
#include <TMath.h>
#include "anacuts.h"


//2-dimensional fit
//x gaus
//y expo x gaus convoluted fit
Double_t K0fit2d(Double_t *x, Double_t *par)
{
  Double_t r1 = (x[0]-par[1])/par[2];
  //Double_t r1 = x[0]*par[1];
  Double_t r2 = TMath::Exp(x[1]*par[3]);

  return par[0]*TMath::Exp(r1*r1)*TMath::Exp(r2);
}


void Fit2DK0()
{
  TFile *fr = TFile::Open("evanaIMpisigma_npippim_v202_out_iso_qhi_sub.root","READ");

  auto *IMnpim_IMnpip_dE_wK0_woSid_n = (TH2F*)fr->Get("IMnpim_IMnpip_dE_wK0_woSid_n");
  auto *IMnpim_IMnpip_dE_wK0_woSid_n_45rot = (TH2F*)fr->Get("IMnpim_IMnpip_dE_wK0_woSid_n_45rot");
  auto *IMnpim_IMnpip_dE_wK0_woSid_n_45rot3 = (TH2F*)fr->Get("IMnpim_IMnpip_dE_wK0_woSid_n_45rot3");
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
  const int binSp = IMnpim_IMnpip_dE_wK0_woSid_n2->GetXaxis()->FindBin(anacuts::Sigmap_center);
  const int binSm = IMnpim_IMnpip_dE_wK0_woSid_n2->GetYaxis()->FindBin(anacuts::Sigmam_center);
  cIMnpim_IMnpip_dE_wK0_woSid_n2->cd(1);
  TH1D *IMnpip_wK0_woSid_n2 = (TH1D*)IMnpim_IMnpip_dE_wK0_woSid_n2->ProjectionX("IMnpip_wK0_woSid_n2",binSp,300);
  IMnpip_wK0_woSid_n2->Draw("E");
  TF1 *fnpip = new TF1("fnpip","expo",1.13,1.26);
  IMnpip_wK0_woSid_n2->Fit(fnpip,"","",1.13,1.26);
  //double parnpip[2];
  //fnpip->GetParameters(parnpip);
  //std::cout << "par0: " << parnpip[0] << std::endl;
  //std::cout << "par1: " << parnpip[1] << std::endl;
  
  cIMnpim_IMnpip_dE_wK0_woSid_n2->cd(4);
  TH1D *IMnpim_wK0_woSid_n2 = (TH1D*)IMnpim_IMnpip_dE_wK0_woSid_n2->ProjectionY("IMnpim_wK0_woSid_n2",binSm,300);
  IMnpim_wK0_woSid_n2->Draw("E");
  

  //TF1 *fnpim = new TF1("fnpim","expo",1.12,1.36);
  //IMnpim_wK0_woSid_n2->Fit(fnpim,"","",1.12,1.36);
  //double parnpim[2];
  //fnpim->GetParameters(parnpim);
  //std::cout << "par0: " << parnpim[0] << std::endl;
  //std::cout << "par1: " << parnpim[1] << std::endl;
  
  //cIMnpim_IMnpip_dE_wK0_woSid_n2->cd(3);
  //const Int_t npar = 3;
  //TF2 *f2 = new TF2("f2",K0fit2d,1.13,1.26,1.12,1.36,3);   
  //TF2 *f2 = new TF2("f2",K0fit2d,1.,1.,1.5,1.5,3);   
  //f2->SetParameters(parnpip[0],parnpip[1],parnpim[1]);
  //IMnpim_IMnpip_dE_wK0_woSid_n2->Fit("f2");
  //f2->Draw("cont1 same");
  //TH2F* hfit = (TH2F*)f2->GetHistogram();
  
  
  //try interpolation
  TH2F *IMnpim_IMnpip_wK0_woSid_n_inter = (TH2F*)IMnpim_IMnpip_dE_wK0_woSid_n2->Clone("IMnpim_IMnpip_wK0_woSid_n_inter");
  
  const int nbinsX = IMnpim_IMnpip_dE_wK0_woSid_n2->GetNbinsX();
  TCanvas *cIMnpim_wK0_woSid_n[nbinsX];
  TH1D* hIMnpim[nbinsX];
  TGraphErrors* grIMnpim_wK0_woSid_n[nbinsX];
  TSpline3 *spIMnpim3[nbinsX];
  TSpline5 *spIMnpim5[nbinsX];
  
  for(int ix=8;ix<nbinsX;ix++){
    hIMnpim[ix] = (TH1D*)IMnpim_IMnpip_dE_wK0_woSid_n2->ProjectionY(Form("IMnpim_wK0_woSid_n%d",ix),ix,ix);
    grIMnpim_wK0_woSid_n[ix] = new TGraphErrors(hIMnpim[ix]);
    double IMnpim_low = IMnpim_IMnpip_dE_wK0_woSid_n2->GetXaxis()->GetBinLowEdge(ix);
    double IMnpim_hi = IMnpim_IMnpip_dE_wK0_woSid_n2->GetXaxis()->GetBinLowEdge(ix+1);
    grIMnpim_wK0_woSid_n[ix]->RemovePoint(7);
    grIMnpim_wK0_woSid_n[ix]->SetTitle(Form("IMnpim_wK0_woSid_n %0.2f-%0.2f",IMnpim_low,IMnpim_hi));
    cIMnpim_wK0_woSid_n[ix] = new TCanvas(Form("cIMnpim%d",ix),Form("cIMnpim%d",ix));
    grIMnpim_wK0_woSid_n[ix]->Draw("AP");
    TBox *b = new TBox(1.1791920,0,1.2152680,30);
    b->SetFillColor(4);
    b->SetFillStyle(3003);
    b->Draw();
    spIMnpim3[ix] = new TSpline3(Form("spIMnpim3%d",ix),grIMnpim_wK0_woSid_n[ix]); 
    spIMnpim3[ix]->SetLineColor(2);
    spIMnpim3[ix]->Draw("same");
    spIMnpim5[ix] = new TSpline5(Form("spIMnpim5%d",ix),grIMnpim_wK0_woSid_n[ix]); 
    spIMnpim5[ix]->SetLineColor(3);
    spIMnpim5[ix]->Draw("same");
    double bincx = IMnpim_IMnpip_dE_wK0_woSid_n2->GetXaxis()->GetBinCenter(ix);
    double bincy = IMnpim_IMnpip_dE_wK0_woSid_n2->GetYaxis()->GetBinCenter(binSm);
    double contx = spIMnpim3[ix]->Eval(bincx);
    if(contx<0.0) contx=0.0;
    //double contx = IMnpim_IMnpip_dE_wK0_woSid_n2->Interpolate(bincx,bincy);
    IMnpim_IMnpip_wK0_woSid_n_inter->SetBinContent(ix,binSm,contx);
    //std::cout << ix  << " " << binyy << std::endl;
    //std::cout << bincx << "  " << bincy << std::endl;
    //std::cout << contx << std::endl;
    //std::cout << std::endl;
  }
  //for(int iy=0;iy<IMnpim_IMnpip_dE_wK0_woSid_n2->GetNbinsY();iy++){
  //  double conty = IMnpim_IMnpip_dE_wK0_woSid_n2->Interpolate(binxx,iy);
  //  IMnpim_IMnpip_wK0_woSid_n_inter->SetBinContent(binxx,iy,conty);
  //}
  

  const int nbinsY = IMnpim_IMnpip_dE_wK0_woSid_n2->GetNbinsY();
  TCanvas *cIMnpip_wK0_woSid_n[nbinsY];
  TH1D* hIMnpip[nbinsY];
  TGraphErrors* grIMnpip_wK0_woSid_n[nbinsY];
  TSpline3 *spIMnpip3[nbinsY];
  TSpline5 *spIMnpip5[nbinsY];
  
  for(int iy=9;iy<nbinsX;iy++){
    hIMnpip[iy] = (TH1D*)IMnpim_IMnpip_dE_wK0_woSid_n2->ProjectionX(Form("IMnpip_wK0_woSid_n%d",iy),iy,iy);
    grIMnpip_wK0_woSid_n[iy] = new TGraphErrors(hIMnpip[iy]);
    double IMnpip_low = IMnpim_IMnpip_dE_wK0_woSid_n2->GetYaxis()->GetBinLowEdge(iy);
    double IMnpip_hi = IMnpim_IMnpip_dE_wK0_woSid_n2->GetYaxis()->GetBinLowEdge(iy+1);
    grIMnpip_wK0_woSid_n[iy]->RemovePoint(6);
    grIMnpip_wK0_woSid_n[iy]->SetTitle(Form("IMnpip_wK0_woSid_n %0.2f-%0.2f",IMnpip_low,IMnpip_hi));
    cIMnpip_wK0_woSid_n[iy] = new TCanvas(Form("cIMnpip%d",iy),Form("cIMnpip%d",iy));
    grIMnpip_wK0_woSid_n[iy]->Draw("AP");
    TBox *b = new TBox(1.1728847,0,1.2053353,30);
    b->SetFillColor(4);
    b->SetFillStyle(3003);
    b->Draw();
    spIMnpip3[iy] = new TSpline3(Form("spIMnpip3%d",iy),grIMnpip_wK0_woSid_n[iy]); 
    spIMnpip3[iy]->SetLineColor(2);
    spIMnpip3[iy]->Draw("same");
    spIMnpip5[iy] = new TSpline5(Form("spIMnpip5%d",iy),grIMnpip_wK0_woSid_n[iy]); 
    spIMnpip5[iy]->SetLineColor(3);
    spIMnpip5[iy]->Draw("same");
    double bincx = IMnpim_IMnpip_dE_wK0_woSid_n2->GetXaxis()->GetBinCenter(binSp);
    double bincy = IMnpim_IMnpip_dE_wK0_woSid_n2->GetYaxis()->GetBinCenter(iy);
    double conty = spIMnpip3[iy]->Eval(bincx);
    if(conty<0) conty=0.0;
    //double contx = IMnpip_IMnpip_dE_wK0_woSid_n2->Interpolate(bincx,bincy);
    IMnpim_IMnpip_wK0_woSid_n_inter->SetBinContent(binSp,iy,conty);
    //std::cout << iy  << " " << binyy << std::endl;
    //std::cout << bincx << "  " << bincy << std::endl;
    //std::cout << contx << std::endl;
    //std::cout << std::endl;
  }
  TCanvas *cinter = new TCanvas("cinter","cinter",800,800);
  cinter->Divide(2,2,0,0);
  cinter->cd(3);
  IMnpim_IMnpip_wK0_woSid_n_inter->Draw("colz");
  cinter->cd(1);
  auto *IMnpip_inter = (TH1D*)IMnpim_IMnpip_wK0_woSid_n_inter->ProjectionX("IMnpip_inter");
  IMnpip_inter->Draw("HE");
  cinter->cd(4);
  auto *IMnpim_inter = (TH1D*)IMnpim_IMnpip_wK0_woSid_n_inter->ProjectionY("IMnpim_inter");
  IMnpim_inter->Draw("HE");
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

  auto *cIMnpim_IMnpip_dE_wK0_woSid_n_45rot3 = new TCanvas("cIMnpim_IMnpip_dE_wK0_woSid_n_45rot3","cIMnpim_IMnpip_dE_wK0_woSid_n_45rot3",800,800);
  cIMnpim_IMnpip_dE_wK0_woSid_n_45rot3->Divide(2,2,0,0);
  cIMnpim_IMnpip_dE_wK0_woSid_n_45rot3->cd(3);
  IMnpim_IMnpip_dE_wK0_woSid_n_45rot3->SetMinimum(0);
  IMnpim_IMnpip_dE_wK0_woSid_n_45rot3->Draw("colz");
  
  cIMnpim_IMnpip_dE_wK0_woSid_n_45rot3->cd(1);
  IMnpim_IMnpip_dE_wK0_woSid_n_45rot3->ProjectionX()->Draw();
  
  cIMnpim_IMnpip_dE_wK0_woSid_n_45rot3->cd(4);
  IMnpim_IMnpip_dE_wK0_woSid_n_45rot3->ProjectionY()->Draw();
   
  auto *crot3 = new TCanvas("crot3","crot3",800,800);
  crot3->Divide(2,2);
  crot3->cd(3);
  IMnpim_IMnpip_dE_wK0_woSid_n_45rot3->Draw("colz");
  const int biny17 = IMnpim_IMnpip_dE_wK0_woSid_n_45rot3->GetYaxis()->FindBin(1.7);
  const int binx00 = IMnpim_IMnpip_dE_wK0_woSid_n_45rot3->GetXaxis()->FindBin(0.0);
  crot3->cd(4);
  auto *py_rot3_right = (TH1D*)IMnpim_IMnpip_dE_wK0_woSid_n_45rot3->ProjectionY("px_rot3_right",binx00,100);
  auto *py_rot3_left = (TH1D*)IMnpim_IMnpip_dE_wK0_woSid_n_45rot3->ProjectionY("px_rot3_left",0,binx00-1);
  py_rot3_right->Draw("E");
  py_rot3_left->SetLineColor(2);
  py_rot3_left->Draw("Esame");
  auto leg = new TLegend(0.6,0.7,0.9,0.9);
  leg->AddEntry(py_rot3_right,"right");
  leg->AddEntry(py_rot3_left,"left");
  leg->Draw();
  
  crot3->cd(1);
  auto *px_rot3 = (TH1D*)IMnpim_IMnpip_dE_wK0_woSid_n_45rot3->ProjectionX("px_rot3");
  px_rot3->Draw("HE");
  crot3->cd(2);
  auto px_rot3_17 = (TH1D*)IMnpim_IMnpip_dE_wK0_woSid_n_45rot3->ProjectionX("xx_rot3",0,biny17);
  px_rot3_17->Draw("HE");


  //crot3->BuildLegend();

  TCanvas *cfitcomp = new TCanvas("cfitcomp","cfitcomp",800,800);
  cfitcomp->Divide(2,2,0,0);
  cfitcomp->cd(3);
  IMnpim_IMnpip_dE_wK0_woSid_n_45rot3->Draw("colz");
  //hfit->Draw("colz");
  TF2 *f2 = new TF2("f2",K0fit2d,1.,1.,1.5,1.5,3);   
  //par0 : scaling factor
  //par1 : x,gaus mean
  //par2 : x,gaus sigma
  //par3 : y,exp slope

  f2->SetParameters(100,0.05,0.01,0.01);
  IMnpim_IMnpip_dE_wK0_woSid_n_45rot3->Fit("f2");
  f2->Draw("cont1 same");

  cfitcomp->cd(1);
  IMnpip_wK0_woSid_n->Draw("E");
  //TH1D* hfitnpip = (TH1D*)hfit->ProjectionX("hfitnpip");
  //hfitnpip->SetLineColor(2);
  //hfitnpip->Draw("HEsame");

  cfitcomp->cd(4);
  IMnpim_wK0_woSid_n->Draw("E");
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
