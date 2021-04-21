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
  Double_t r2 = (x[1]*par[3]);

  return par[0]*TMath::Exp(-0.5*r1*r1)*TMath::Exp(r2);
}


void Fit2DK0()
{
  TFile *fr = TFile::Open("evanaIMpisigma_npippim_v202_out_iso_qhi_sub.root","READ");

  auto *IMnpim_IMnpip_dE_wK0_woSid_n_1 = (TH2F*)fr->Get("IMnpim_IMnpip_dE_wK0_woSid_n");
  auto *IMnpim_IMnpip_dE_wK0_woSid_n_45rot = (TH2F*)fr->Get("IMnpim_IMnpip_dE_wK0_woSid_n_45rot");
  auto *IMnpim_IMnpip_dE_wK0_woSid_n_45rot3 = (TH2F*)fr->Get("IMnpim_IMnpip_dE_wK0_woSid_n_45rot3");
  auto *cIMnpim_IMnpip_dE_wK0_woSid_n = new TCanvas("cIMnpim_IMnpip_dE_wK0_woSid_n","cIMnpim_IMnpip_dE_wK0_woSid_n",800,800);
  
  cIMnpim_IMnpip_dE_wK0_woSid_n->Divide(2,2,0.,0.);
  cIMnpim_IMnpip_dE_wK0_woSid_n->cd(3);
  //IMnpim_IMnpip_dE_wK0_woSid_n_1->RebinX(4);
  //IMnpim_IMnpip_dE_wK0_woSid_n_1->RebinY(4);
  //IMnpim_IMnpip_dE_wK0_woSid_n->GetXaxis()->SetRangeUser(1.0,1.8);
  //IMnpim_IMnpip_dE_wK0_woSid_n->GetYaxis()->SetRangeUser(1.0,1.8);
  IMnpim_IMnpip_dE_wK0_woSid_n_1->SetMinimum(0);
  IMnpim_IMnpip_dE_wK0_woSid_n_1->Draw("colz");
  
  cIMnpim_IMnpip_dE_wK0_woSid_n->cd(1);
  TH1D *IMnpip_wK0_woSid_n = (TH1D*)IMnpim_IMnpip_dE_wK0_woSid_n_1->ProjectionX("IMnpip_wK0_woSid_n");
  IMnpip_wK0_woSid_n->Draw("HE");
  
  cIMnpim_IMnpip_dE_wK0_woSid_n->cd(4);
  TH1D *IMnpim_wK0_woSid_n = (TH1D*)IMnpim_IMnpip_dE_wK0_woSid_n_1->ProjectionY("IMnpim_wK0_woSid_n");
  IMnpim_wK0_woSid_n->Draw("HE");
   
  auto *cIMnpim_IMnpip_dE_wK0_woSid_n2 = new TCanvas("cIMnpim_IMnpip_dE_wK0_woSid_n2","cIMnpim_IMnpip_dE_wK0_woSid_n2",800,800);
  cIMnpim_IMnpip_dE_wK0_woSid_n2->Divide(2,2,0.,0.);
  cIMnpim_IMnpip_dE_wK0_woSid_n2->cd(3);
  TH2D* IMnpim_IMnpip_dE_wK0_woSid_n2 = (TH2D*)IMnpim_IMnpip_dE_wK0_woSid_n_1->Clone("IMnpim_IMnpip_dE_wK0_woSid_n2");
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


  TCanvas *cfitcomp = new TCanvas("cfitcomp","cfitcomp",800,800);
  cfitcomp->Divide(2,2);
  cfitcomp->cd(3);
  auto *IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_2 = (TH2D*)IMnpim_IMnpip_dE_wK0_woSid_n_45rot3->Clone("IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_2");
  for(int ix=0;ix<IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_2->GetNbinsX();ix++){
    for(int iy=0;iy<IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_2->GetNbinsY();iy++){
      double cont = IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_2->GetBinContent(ix,iy);
      if(cont<0.00001){
        IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_2->SetBinContent(ix,iy,0);
        IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_2->SetBinError(ix,iy,0);
      }
    }
  }
  IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_2->Draw("colz");
  IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_2->GetYaxis()->SetRangeUser(1.660,2.1);
  //hfit->Draw("colz");
  TF2 *f2 = new TF2("f2",K0fit2d,-0.5,0.5,1.6,2.2,4);
  //par0 : scaling factor
  //par1 : x,gaus mean
  //par2 : x,gaus sigma
  //par3 : y,exp slope

  for(int ix=0;ix<IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_2->GetNbinsX();ix++){
    for(int iy=0;iy<IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_2->GetNbinsY();iy++){
      double cont = IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_2->GetBinContent(ix,iy);
      if(fabs(cont)<0.001){ 
        //IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_2->SetBinContent(ix,iy,0);
      }
    }
  }
  //f2->SetRange(-0.5,0.5,1.660,2.1,4); 
  f2->SetRange(-0.4,1.666,0.4,1.85); 
  f2->SetParameters(8.0e9,0.005,0.16,-11.2);
  //f2->SetParLimits(0,3.9e9,4.5e9);
  //f2->SetParLimits(1,0.0,0.1);
  //f2->FixParameter(1,0.005);
  //f2->SetParLimits(2,0.15,0.2);
  //f2->SetParLimits(3,-13,-12);
  //f2->FixParameter(3,-12.5);
  IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_2->Fit("f2","R","");
  f2->Draw("cont1 same");
  TF2 *f2wide = new TF2("f2wide",K0fit2d,-0.5,0.5,1.6,2.2,4);
  Double_t param[4];
  f2->GetParameters(param);
  f2wide->SetParameters(param);
  f2wide->SetNpx(100);
  f2wide->SetNpy(60);
  std::cout<<f2wide->GetExpFormula() << std::endl ;
  TH2D *hf2  =  (TH2D*)f2->GetHistogram();
  hf2->SetName("hf2");
  TH2D *hf2wide  =  (TH2D*)f2wide->GetHistogram();
  hf2wide->SetName("hf2wide");
  TH2D *hf2wide_nosub = (TH2D*)hf2wide->Clone("hf2wide_nosub");
  for(int ix=0;ix<IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_2->GetNbinsX();ix++){
    for(int iy=0;iy<IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_2->GetNbinsY();iy++){
      double cont = IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_2->GetBinContent(ix,iy);
      if((cont)<0.00000001) hf2wide->SetBinContent(ix,iy,0);
      double xcen=  IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_2->GetXaxis()->GetBinCenter(ix);
      double ycen=  IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_2->GetYaxis()->GetBinCenter(iy);

      if( (fabs(xcen)<0.02) && (ycen < 1.68)){
         hf2wide->SetBinContent(ix,iy,0);
         IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_2->SetBinContent(ix,iy,0);
      }
    }
  }

  IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_2->Print("base");
  hf2->Print("base");
  cfitcomp->cd(1);
  TH1D* pxrot3 = (TH1D*)IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_2->ProjectionX("pxrot3");
  pxrot3->Draw("HE");
  hf2->SetFillColor(0);
  hf2->ProjectionX()->Draw("HISTsame");
  
  cfitcomp->cd(4);
  TH1D* pyrot3 = (TH1D*)IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_2->ProjectionY("pyrot3");
  pyrot3->Draw("HE");
  hf2->ProjectionY()->Draw("HISTsame");
  
  //subtract and plot wide fit
  TCanvas *cfitcompsub = new TCanvas("cfitcompsub","cfitcompsub",800,800);
  cfitcompsub->Divide(2,2);
  cfitcompsub->cd(3);
  IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_2->Draw("colz");
  hf2wide->Draw("box same");

  cfitcompsub->cd(1);
  pxrot3->Draw("HE");
  hf2wide->SetFillColor(0);
  hf2wide->ProjectionX()->Draw("HISTsame");

  cfitcompsub->cd(4);
  pyrot3->Draw("HE");
  hf2wide->ProjectionY()->Draw("HISTsame");
  
  auto *IMnpim_IMnpip_dE_wK0_woSid_n_2 = (TH2F*)fr->Get("IMnpim_IMnpip_dE_wK0_woSid_n");
  IMnpim_IMnpip_dE_wK0_woSid_n_2->SetName("IMnpim_IMnpip_dE_wK0_woSid_n_2");
  TCanvas *cinter = new TCanvas("cinter","cinter",800,800);
  cinter->Divide(2,2);
  
  for(int ix=0;ix<IMnpim_IMnpip_dE_wK0_woSid_n_2->GetNbinsX();ix++){
    for(int iy=0;iy<IMnpim_IMnpip_dE_wK0_woSid_n_2->GetNbinsY();iy++){
      double xcent = IMnpim_IMnpip_dE_wK0_woSid_n_2->GetXaxis()->GetBinCenter(ix);
      double ycent = IMnpim_IMnpip_dE_wK0_woSid_n_2->GetYaxis()->GetBinCenter(iy);
      if( (anacuts::Sigmap_MIN_wide < xcent && xcent < anacuts::Sigmap_MAX_wide)
        ||(anacuts::Sigmam_MIN_wide < ycent && ycent < anacuts::Sigmam_MAX_wide)){
        double xx = 1./sqrt(2.0)*(xcent-ycent);
        double yy = 1./sqrt(2.0)*(xcent+ycent);
        double yy2 = yy-(cosh(1.96*xx)-1.0);
        double evalK0 = f2wide->Eval(xx,yy2);
        double scale=0.2;
        evalK0 *= scale;
        if(yy2>1.645 && xcent < 1.7 && ycent<1.7)
        IMnpim_IMnpip_dE_wK0_woSid_n_2->SetBinContent(ix,iy,evalK0);
      }
    }
  }
  IMnpim_IMnpip_dE_wK0_woSid_n_2->Rebin2D(4,4);
  cinter->cd(3);
  IMnpim_IMnpip_dE_wK0_woSid_n_2->Draw("colz");

  cinter->cd(1);
  IMnpim_IMnpip_dE_wK0_woSid_n_2->ProjectionX("pxinter")->Draw("EH");

  cinter->cd(4);
  IMnpim_IMnpip_dE_wK0_woSid_n_2->ProjectionY("pyinter")->Draw("EH");

}
