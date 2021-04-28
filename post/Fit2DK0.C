//fit macro to decompose K0/S+/S- final states
#include <TF2.h>
#include <TH2.h>
#include <TMath.h>
#include "anacuts.h"


//2-dimensional fit
//x gaus
//y expo x gaus convoluted fit
const int nparfit=4;
Double_t K0fit2d(Double_t *x, Double_t *par)
{
  Double_t r1 = (x[0]-par[1])/par[2];
  //Double_t r1 = x[0]*par[1];
  Double_t r2 = (x[1]*par[3]);
  //Double_t r2 = (x[1]-par[3])/par[4];

  return par[0]*TMath::Exp(-0.5*r1*r1)*TMath::Exp(-1.0*r2*r2);
  //return par[0]*TMath::Exp(-0.5*r1*r1)*TMath::Exp(-0.5*r2*r2);
}


void Fit2DK0(const int qcut=2)
{
  
  TFile *fr = NULL;
  if(qcut==1){
    fr = TFile::Open("evanaIMpisigma_npippim_v202_out_iso_qlo_sub.root","READ");
  }else if(qcut==2){
    fr = TFile::Open("evanaIMpisigma_npippim_v202_out_iso_qhi_sub.root","READ");
  }else{
    std::cout << "no file" << std::endl;
    return;
  }
  auto *IMnpim_IMnpip_dE_wK0_woSid_n_1 = (TH2F*)fr->Get("IMnpim_IMnpip_dE_wK0_woSid_n");
  auto *IMnpim_IMnpip_dE_wK0_woSid_n_45rot = (TH2F*)fr->Get("IMnpim_IMnpip_dE_wK0_woSid_n_45rot");
  auto *IMnpim_IMnpip_dE_wK0_woSid_n_45rot3 = (TH2F*)fr->Get("IMnpim_IMnpip_dE_wK0_woSid_n_45rot3");
  auto *cIMnpim_IMnpip_dE_wK0_woSid_n = new TCanvas("cIMnpim_IMnpip_dE_wK0_woSid_n","cIMnpim_IMnpip_dE_wK0_woSid_n",800,800);
  
  cIMnpim_IMnpip_dE_wK0_woSid_n->Divide(2,2,0.,0.);
  cIMnpim_IMnpip_dE_wK0_woSid_n->cd(3);
  //IMnpim_IMnpip_dE_wK0_woSid_n_1->Rebin2D(4,4);
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
      if(cont<0.001){
        IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_2->SetBinContent(ix,iy,0);
        IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_2->SetBinError(ix,iy,0);
      }
    }
  }
  IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_2->Draw("colz");
  IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_2->GetYaxis()->SetRangeUser(1.660,2.1);
  TF2 *f2 = new TF2("f2",K0fit2d,-0.5,0.5,1.6,2.2,nparfit);
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
  //f2->SetRange(-0.4,1.666,0.4,1.85); 
  f2->SetRange(-0.4,1.666,0.4,1.85); 
  //f2->SetParameters(8.0e9,0.005,0.16,-15.2);
  f2->SetParameters(2.0e5,0.005,0.16,1.9);
  f2->SetParLimits(0,0,4.5e10);
  f2->FixParameter(0,2.0e5);
  f2->SetParLimits(1,0.0,0.1);
  //f2->FixParameter(1,0.005);
  f2->SetParLimits(2,0.15,0.2);
  //f2->SetParameter(3,0.5);
  //f2->SetParLimits(3,1.66,3.00);
  f2->FixParameter(3,1.9);
  IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_2->Fit("f2","R","");
  f2->Draw("cont1 same");
  TF2 *f2wide = new TF2("f2wide",K0fit2d,-0.5,0.5,1.6,2.2,nparfit);
  Double_t param[nparfit];
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
  
  auto *IMnpim_IMnpip_dE_wK0_woSid_n_2 = (TH2D*)IMnpim_IMnpip_dE_wK0_woSid_n_1->Clone("IMnpim_IMnpip_dE_wK0_woSid_n_2");
  IMnpim_IMnpip_dE_wK0_woSid_n_2->SetName("IMnpim_IMnpip_dE_wK0_woSid_n_2");
  TCanvas *cinter = new TCanvas("cinter","cinter",800,800);
  cinter->Divide(2,2);
  TH2D* h2K0inter = (TH2D*)IMnpim_IMnpip_dE_wK0_woSid_n_2->Clone("h2K0inter");
  h2K0inter->Reset();
  for(int ix=0;ix<IMnpim_IMnpip_dE_wK0_woSid_n_2->GetNbinsX();ix++){
    for(int iy=0;iy<IMnpim_IMnpip_dE_wK0_woSid_n_2->GetNbinsY();iy++){
      double xcent = IMnpim_IMnpip_dE_wK0_woSid_n_2->GetXaxis()->GetBinCenter(ix);
      double ycent = IMnpim_IMnpip_dE_wK0_woSid_n_2->GetYaxis()->GetBinCenter(iy);
      double cont  = IMnpim_IMnpip_dE_wK0_woSid_n_2->GetBinContent(ix,iy);
      if(cont < 0.0001){
        if( (anacuts::Sigmap_MIN_wide < xcent && xcent < anacuts::Sigmap_MAX_wide)
            ||(anacuts::Sigmam_MIN_wide < ycent && ycent < anacuts::Sigmam_MAX_wide)){
          double xx = 1./sqrt(2.0)*(xcent-ycent);
          double yy = 1./sqrt(2.0)*(xcent+ycent);
          double yy2 = yy-(cosh(1.96*xx)-1.0);
          double evalK0 = f2wide->Eval(xx,yy2);
          double scale=0.15;
          evalK0 *= scale;
          if(yy2>1.645 && xcent < 1.7 && ycent<1.7){
            IMnpim_IMnpip_dE_wK0_woSid_n_2->SetBinContent(ix,iy,evalK0);
            h2K0inter->SetBinContent(ix,iy,evalK0);
          }
        }
      }
    }
  }
  IMnpim_IMnpip_dE_wK0_woSid_n_2->Rebin2D(4,4);
  h2K0inter->Rebin2D(4,4);
  cinter->cd(3);
  IMnpim_IMnpip_dE_wK0_woSid_n_2->Draw("colz");

  cinter->cd(1);
  TH1D* pxinter= (TH1D*)IMnpim_IMnpip_dE_wK0_woSid_n_2->ProjectionX("pxinter");
  pxinter->Draw("EH");
  TH1D* K0interpx = (TH1D*)h2K0inter->ProjectionX("K0interpx");
  K0interpx->SetFillColor(2);
  K0interpx->Draw("HISTsame");
  cinter->cd(4);
  TH1D* pyinter= (TH1D*)IMnpim_IMnpip_dE_wK0_woSid_n_2->ProjectionY("pyinter");
  pyinter->Draw("EH");
  TH1D* K0interpy = (TH1D*)h2K0inter->ProjectionY("K0interpy");
  K0interpy->SetFillColor(2);
  K0interpy->Draw("HISTsame");

  TH2F* IMnpim_IMnpip_dE_wK0orwSid_n = (TH2F*)fr->Get("IMnpim_IMnpip_dE_wK0orwSid_n");
  //TH2F* IMnpim_IMnpip_dE_wSid_n = (TH2F*)fr->Get("IMnpim_IMnpip_dE_wSid_n");
  auto *cwK0orwSid_n = new TCanvas("cwK0orwSid_n","cwK0orwSid_n",1600,800);
  cwK0orwSid_n->Divide(2,1);
  cwK0orwSid_n->cd(1);
  IMnpim_IMnpip_dE_wK0orwSid_n->Rebin2D(4,4);
  IMnpim_IMnpip_dE_wK0orwSid_n->Draw("colz");
  //IMnpim_IMnpip_dE_wSid_n->Rebin2D(4,4);
  //IMnpim_IMnpip_dE_wSid_n->Draw("colz");
  
  cwK0orwSid_n->cd(2);
  TH2F* IMnpim_IMnpip_dE_wK0orwSid_n_K0sub = (TH2F*)IMnpim_IMnpip_dE_wK0orwSid_n->Clone("IMnpim_IMnpip_dE_wK0orwSid_n_K0sub");
  IMnpim_IMnpip_dE_wK0orwSid_n_K0sub->Add(IMnpim_IMnpip_dE_wK0_woSid_n_2,-1.0);
  IMnpim_IMnpip_dE_wK0orwSid_n_K0sub->SetMaximum(IMnpim_IMnpip_dE_wK0orwSid_n->GetMaximum());
  IMnpim_IMnpip_dE_wK0orwSid_n_K0sub->SetMinimum(0);
  IMnpim_IMnpip_dE_wK0orwSid_n_K0sub->Draw("colz");
  

  auto *cwK0orwSid_n_pro = new TCanvas("cwK0orwSid_n_pro","cwK0orwSid_n_pro",1600,800);
  cwK0orwSid_n_pro->Divide(2,1);
  cwK0orwSid_n_pro->cd(1);
  TH1D* IMnpip_wK0orwSid_n = (TH1D*)IMnpim_IMnpip_dE_wK0orwSid_n->ProjectionX("IMnpip_wK0orwSid_n");
  IMnpip_wK0orwSid_n->Draw("HE");
  pxinter->SetLineColor(2);
  pxinter->Draw("HEsame");

  cwK0orwSid_n_pro->cd(2);
  TH1D* IMnpim_wK0orwSid_n = (TH1D*)IMnpim_IMnpip_dE_wK0orwSid_n->ProjectionY("IMnpim_wK0orwSid_n");
  IMnpim_wK0orwSid_n->Draw("HE");
  pyinter->SetLineColor(2);
  pyinter->Draw("HEsame");
   

  //K0 subtracted events
  auto *cwSid_n_K0sub = new TCanvas("cwSid_n_K0sub","cwK0orwSid_n",800,800);
  cwSid_n_K0sub->Divide(2,2);
  cwSid_n_K0sub->cd(3);
  IMnpim_IMnpip_dE_wK0orwSid_n_K0sub->Draw("colz");

  cwSid_n_K0sub->cd(1);
  const int SmbinStart = IMnpim_IMnpip_dE_wK0orwSid_n_K0sub->GetYaxis()->FindBin(anacuts::Sigmam_MIN);
  std::cout << "SmbinStart low Edge " << IMnpim_IMnpip_dE_wK0orwSid_n_K0sub->GetYaxis()->GetBinLowEdge(SmbinStart) << std::endl;
  std::cout << "Sigmam_MIN " << anacuts::Sigmam_MIN << std::endl;
  std::cout << "Smbin width " << IMnpim_IMnpip_dE_wK0orwSid_n_K0sub->GetYaxis()->GetBinWidth(SmbinStart) << std::endl;
  std::cout << "Sigmam_sigma " << anacuts::Sigmam_sigma << std::endl;
  const int SmbinEnd = IMnpim_IMnpip_dE_wK0orwSid_n_K0sub->GetYaxis()->FindBin(anacuts::Sigmam_MAX);
  std::cout << "SmbinEnd low Edge " << IMnpim_IMnpip_dE_wK0orwSid_n_K0sub->GetYaxis()->GetBinLowEdge(SmbinEnd) << std::endl;
  std::cout << "SmbinEnd low Edge+1 " << IMnpim_IMnpip_dE_wK0orwSid_n_K0sub->GetYaxis()->GetBinLowEdge(SmbinEnd+1) << std::endl;
  std::cout << "Sigmam_MAX " << anacuts::Sigmam_MAX << std::endl;
  const int Spbin = IMnpim_IMnpip_dE_wK0orwSid_n_K0sub->GetXaxis()->FindBin(anacuts::Sigmap_center);
  TH1D* IMnpip_K0sub = (TH1D*)IMnpim_IMnpip_dE_wK0orwSid_n_K0sub->ProjectionX("IMnpip_K0sub",SmbinStart,SmbinEnd);
  TH1D* IMnpip_K0sub_Sp = (TH1D*)IMnpim_IMnpip_dE_wK0orwSid_n_K0sub->ProjectionX("IMnpip_K0sub_Sp",SmbinStart,SmbinEnd);
  IMnpip_K0sub->Draw("HE");
  IMnpip_K0sub_Sp->SetFillColor(2);
  IMnpip_K0sub_Sp->GetXaxis()->SetRange(Spbin,Spbin);
  IMnpip_K0sub_Sp->SetFillStyle(3002);
  IMnpip_K0sub_Sp->Draw("HEsame");

  cwSid_n_K0sub->cd(4);
  const int SpbinStart = IMnpim_IMnpip_dE_wK0orwSid_n_K0sub->GetXaxis()->FindBin(anacuts::Sigmap_MIN);
  std::cout << "SpbinStart low Edge " << IMnpim_IMnpip_dE_wK0orwSid_n_K0sub->GetXaxis()->GetBinLowEdge(SpbinStart) << std::endl;
  std::cout << "Sigmap_MIN " << anacuts::Sigmap_MIN << std::endl;
  const int SpbinEnd = IMnpim_IMnpip_dE_wK0orwSid_n_K0sub->GetXaxis()->FindBin(anacuts::Sigmap_MAX);
  std::cout << "SpbinEnd low Edge " << IMnpim_IMnpip_dE_wK0orwSid_n_K0sub->GetXaxis()->GetBinLowEdge(SpbinEnd) << std::endl;
  std::cout << "SpbinEnd low Edge+1 " << IMnpim_IMnpip_dE_wK0orwSid_n_K0sub->GetXaxis()->GetBinLowEdge(SpbinEnd+1) << std::endl;
  std::cout << "Sigmap_MAX " << anacuts::Sigmap_MAX << std::endl;
  const int Smbin = IMnpim_IMnpip_dE_wK0orwSid_n_K0sub->GetYaxis()->FindBin(anacuts::Sigmam_center);
  TH1D* IMnpim_K0sub = (TH1D*)IMnpim_IMnpip_dE_wK0orwSid_n_K0sub->ProjectionY("IMnpim_K0sub",SpbinStart,SpbinEnd);
  IMnpim_K0sub->Draw("HE");
  TH1D* IMnpim_K0sub_Sm = (TH1D*)IMnpim_IMnpip_dE_wK0orwSid_n_K0sub->ProjectionY("IMnpim_K0sub_Sm",SpbinStart,SpbinEnd);
  IMnpim_K0sub_Sm->SetFillColor(2);
  IMnpim_K0sub_Sm->GetXaxis()->SetRange(Smbin,Smbin);
  IMnpim_K0sub_Sm->SetFillStyle(3002);
  IMnpim_K0sub_Sm->Draw("HEsame");
  std::cout << "overlap events"  << IMnpim_K0sub_Sm->GetBinContent(Smbin) << std::endl;
  std::cout << "overlap error"  << IMnpim_K0sub_Sm->GetBinError(Smbin) << std::endl;
  //remove signal each other and apply interpolation in crossing region
  TH1D* IMnpip_K0sub_woSp = (TH1D*)IMnpip_K0sub->Clone("IMnpip_K0sub_woSp");
  TH1D* IMnpim_K0sub_woSm = (TH1D*)IMnpim_K0sub->Clone("IMnpim_K0sub_woSm");
  
  for(int ibin=SpbinStart;ibin<=SpbinEnd;ibin++){
    IMnpip_K0sub_woSp->SetBinContent(ibin,0);
    IMnpip_K0sub_woSp->SetBinError(ibin,0);
  }

  for(int ibin=SmbinStart;ibin<=SmbinEnd;ibin++){
    IMnpim_K0sub_woSm->SetBinContent(ibin,0);
    IMnpim_K0sub_woSm->SetBinError(ibin,0);
  }

  /*
  for(int ibin=0;ibin<IMnpip_K0sub_woSp->GetNbinsX();ibin++){
    double bincen = IMnpip_K0sub_woSp->GetBinCenter(ibin);
    if(anacuts::Sigmap_MIN_wide < bincen && bincen < anacuts::Sigmap_MAX_wide){
       IMnpip_K0sub_woSp->SetBinContent(ibin,0);
       IMnpip_K0sub_woSp->SetBinError(ibin,0);
       std::cout << "Sp" << std::endl;
    }
  }
  

  for(int ibin=0;ibin<IMnpim_K0sub_woSm->GetNbinsX();ibin++){
    double bincen = IMnpim_K0sub_woSm->GetBinCenter(ibin);
    if(anacuts::Sigmam_MIN_wide < bincen && bincen < anacuts::Sigmam_MAX_wide){
       IMnpim_K0sub_woSm->SetBinContent(ibin,0);
       IMnpim_K0sub_woSm->SetBinError(ibin,0);
       std::cout << "Sm" << std::endl;
    }
  }*/


  auto *cwSid_n_K0sub_wo = new TCanvas("cwSid_n_K0sub_wo","cwK0orwSid_n_K0sub_wo",1600,800);
  cwSid_n_K0sub_wo->Divide(2,1);
  cwSid_n_K0sub_wo->cd(1);
  IMnpip_K0sub_woSp->Draw("HE");
  //IMnpip_K0sub_woSp->Rebin(2);
  TGraphErrors *gIMnpip_K0sub_woSp = new TGraphErrors(IMnpip_K0sub_woSp);
  //gIMnpip_K0sub_woSp->Draw("c");
  //gIMnpip_K0sub_woSp->Print("all");
  gIMnpip_K0sub_woSp->RemovePoint(6);
 // gIMnpip_K0sub_woSp->RemovePoint(12);
  TSpline3 *sIMnpip_K0sub_woSp = new TSpline3("snpip",gIMnpip_K0sub_woSp);
  TSpline5 *s5IMnpip_K0sub_woSp = new TSpline5("s5npip",gIMnpip_K0sub_woSp);
  
  sIMnpip_K0sub_woSp->SetLineColor(3);
  sIMnpip_K0sub_woSp->SetLineWidth(3);
  sIMnpip_K0sub_woSp->Draw("same");

  cwSid_n_K0sub_wo->cd(2);
  //IMnpim_K0sub_woSm->Rebin(2);
  IMnpim_K0sub_woSm->Draw("HE");
  TGraphErrors *gIMnpim_K0sub_woSm = new TGraphErrors(IMnpim_K0sub_woSm);
  //gIMnpip_K0sub_woSp->Print("all");
  gIMnpim_K0sub_woSm->RemovePoint(7);
 // gIMnpim_K0sub_woSm->RemovePoint(6);
  //gIMnpim_K0sub_woSm->Draw("c");
  TSpline3 *sIMnpim_K0sub_woSm = new TSpline3("sppim",gIMnpim_K0sub_woSm);
  TSpline5 *s5IMnpim_K0sub_woSm = new TSpline5("s5ppim",gIMnpim_K0sub_woSm);
  
  sIMnpim_K0sub_woSm->SetLineColor(3);
  sIMnpim_K0sub_woSm->SetLineWidth(3);
  sIMnpim_K0sub_woSm->Draw("same");


  cwSid_n_K0sub->cd(1);
  //sIMnpip_K0sub_woSp->Scale(0.5);
  sIMnpip_K0sub_woSp->Draw("same");
  std::cout << "Sp estimated:" << sIMnpip_K0sub_woSp->Eval(anacuts::Sigmap_center) << std::endl;
  cwSid_n_K0sub->cd(4);
  sIMnpim_K0sub_woSm->Draw("same");
  std::cout << "Sm estimated:" << sIMnpim_K0sub_woSm->Eval(anacuts::Sigmam_center) << std::endl;


}
