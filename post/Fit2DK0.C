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
  Double_t r2 = (x[1]*par[3]);
  
  //woods-saxon shape K0 threshold in x[1] axis
  const Double_t th = 1.00;
  const Double_t a  = 0.002;
  
  return par[0]*TMath::Exp(-0.5*r1*r1)*TMath::Exp(-1.0*r2*r2)/(1.0+TMath::Exp((-x[1]+th)/a));
  //return par[0]*TMath::Exp(-0.5*r1*r1)*TMath::Exp(-0.5*r2*r2);
}


//2D fit for K0nn modeling 
//x[0]: IM(npi+)
//x[1]: IM(npi-)
Double_t K0fit2dNoconvert(Double_t *x,Double_t *par)
{
  double xx = 1.0/sqrt(2.0)*(x[0]-x[1]);
  double yy = 1.0/sqrt(2.0)*(x[0]+x[1]);
  //double yy2 = yy-(sqrt(6.76*xx*xx+2.725)-1.0);
  double yy2 = yy-(sqrt(6.76*xx*xx+2.765)-1.0);//slightly tuned from original value by looking final fitting result

  Double_t r1 = (xx-par[1])/par[2];
  Double_t r2 = yy2*par[3];
  const Double_t th = 1.00;
  const Double_t a  = 0.002;
  const Double_t thend = 2.00;
  const Double_t aend = 0.02;
  double ret = par[0]*TMath::Exp(-0.5*r1*r1)*TMath::Exp(-1.0*r2*r2)/(1.0+TMath::Exp((-yy2+th)/a))/(1.0+TMath::Exp((yy-thend)/aend));    
  //double ret = par[0]*TMath::Exp(-0.5*r1*r1)*TMath::Exp(-1.0*r2*r2)/(1.0+TMath::Exp((-yy2+th)/a))/(1.0+TMath::Exp((yy-thend)/aend));    
  return ret;
  /*
  //if(yy2>1.0){
  if(yy2>1.0){
    return  ret;
   }else{
    return 0;
  }*/
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
  //IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_2->GetYaxis()->SetRangeUser(1.0,1.5);
  TF2 *f2 = new TF2("f2",K0fit2d,-0.5,0.5,0.9,1.5,nparfit);
  //par0 : scaling factor
  //par1 : x,gaus mean
  //par2 : x,gaus sigma
  //par3 : y,exp slope

  //f2->SetRange(-0.5,0.5,1.660,2.1,4); 
  //f2->SetRange(-0.4,1.666,0.4,1.85); 
  f2->SetRange(-0.4,0.4,0.9,1.4); 
  //f2->SetParameters(8.0e9,0.005,0.16,-15.2);
  f2->SetParameters(2.0e5,0.005,0.16,1.9);
  f2->SetParLimits(0,0,4.5e10);
  //f2->FixParameter(0,2.0e5);
  f2->SetParLimits(1,0.0,0.1);
  //f2->FixParameter(1,0.005);
  f2->SetParLimits(2,0.15,0.2);
  //f2->SetParameter(3,0.5);
  //f2->SetParLimits(3,1.66,3.00);
  f2->FixParameter(3,1.9);
  IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_2->Fit("f2","R","");
  IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_2->Print("base");
  f2->Draw("cont1 same");
  TF2 *f2wide = new TF2("f2wide",K0fit2d,-0.5,0.5,0.9,1.5,nparfit);
  Double_t param[nparfit];
  f2->GetParameters(param);
  f2wide->SetParameters(param);
  f2wide->SetNpx(100);//=NbinsX of rot3 histogram
  f2wide->SetNpy(120);//=NbinsY of rot3 histogram
  std::cout<<f2wide->GetExpFormula() << std::endl ;
  TH2D *hf2  =  (TH2D*)f2->GetHistogram();
  hf2->SetName("hf2");
  TH2D *hf2wide  =  (TH2D*)f2wide->GetHistogram();
  hf2wide->SetName("hf2wide");
  TH2D *hf2wide_nosub = (TH2D*)hf2wide->Clone("hf2wide_nosub");
  for(int ix=0;ix<=IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_2->GetNbinsX();ix++){
    for(int iy=0;iy<=IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_2->GetNbinsY();iy++){
      double cont = IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_2->GetBinContent(ix,iy);
      if((cont)<0.00000001) hf2wide->SetBinContent(ix,iy,0);
      double xcen=  IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_2->GetXaxis()->GetBinCenter(ix);
      double ycen=  IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_2->GetYaxis()->GetBinCenter(iy);
      //remove edge area
      if( (fabs(xcen)<0.02) && (ycen < 1.02)){
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
   

  TCanvas *ctest = new TCanvas("ctest","ctest",1600,800);
  ctest->Divide(2,1);
  ctest->cd(1);
  hf2wide->Draw("colz");
  ctest->cd(2);
  hf2wide_nosub->SetTitle("hf2wide_nosub");
  hf2wide_nosub->Draw("colz");


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
          double yy2 = yy-(sqrt(6.76*xx*xx+2.725)-1.0);
          double evalK0 = f2wide->Eval(xx,yy2);
          double scale=0.28;//this scaling factor is arbitrary and determined by eye at this moment
          evalK0 *= scale;
          if(yy2>0.98 && xcent < 1.7 && ycent<1.7){
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
  
  //fit to original histogram
  auto *IMnpim_IMnpip_dE_wK0_woSid_n_3 = (TH2D*)IMnpim_IMnpip_dE_wK0_woSid_n_1->Clone("IMnpim_IMnpip_dE_wK0_woSid_n_3");
  auto *cK0fit = new TCanvas("cK0fit","cK0fit",800,800);
  cK0fit->Divide(2,2);
  cK0fit->cd(3);
  for(int ix=0;ix<IMnpim_IMnpip_dE_wK0_woSid_n_3->GetNbinsX();ix++){
    for(int iy=0;iy<IMnpim_IMnpip_dE_wK0_woSid_n_3->GetNbinsY();iy++){
      double cont = IMnpim_IMnpip_dE_wK0_woSid_n_3->GetBinContent(ix,iy);
      double xcent = IMnpim_IMnpip_dE_wK0_woSid_n_3->GetXaxis()->GetBinCenter(ix);
      double ycent = IMnpim_IMnpip_dE_wK0_woSid_n_3->GetYaxis()->GetBinCenter(iy);
      //remove edge region and negative bin
      //if((cont < 0.0001) || (xcent<1.18 && ycent<1.18)) {       
      if( (xcent<1.18 && ycent<1.18)) {       
        IMnpim_IMnpip_dE_wK0_woSid_n_3->SetBinContent(ix,iy,0);
        IMnpim_IMnpip_dE_wK0_woSid_n_3->SetBinError(ix,iy,0);
      }
    }
  }
  
  //default binning +/- 0.5 sigma=1 sigma.
  //woSid +/-3 sigma cut (no rounding) 
  //
  //fit with wider bin (+/- 3sigma) ->fail
  //IMnpim_IMnpip_dE_wK0_woSid_n_3->Rebin2D(3,3);
  IMnpim_IMnpip_dE_wK0_woSid_n_3->Draw("colz");
   
  const float xmin = IMnpim_IMnpip_dE_wK0_woSid_n_3->GetXaxis()->GetXmin();
  const float xmax = IMnpim_IMnpip_dE_wK0_woSid_n_3->GetXaxis()->GetXmax();
  const float ymin = IMnpim_IMnpip_dE_wK0_woSid_n_3->GetYaxis()->GetXmin();
  const float ymax = IMnpim_IMnpip_dE_wK0_woSid_n_3->GetYaxis()->GetXmax();
  
  const int nbinsX =IMnpim_IMnpip_dE_wK0_woSid_n_3->GetNbinsX(); 
  const int nbinsY =IMnpim_IMnpip_dE_wK0_woSid_n_3->GetNbinsY(); 
  //fit to the original 2d histo
  TF2 *f3 = new TF2("f3",K0fit2dNoconvert,xmin,xmax,ymin,ymax,nparfit);
  f3->SetNpx(nbinsX);//use same nbin to compare the projection
  f3->SetNpy(nbinsY);//use same nbin to compare the projection
  f3->SetParameters(param);
  //f3->SetParameter(0,param[0]/2.0);
  //f3->FixParameter(0,param[0]*0.3);
  f3->SetParLimits(0,param[0]*0.2,param[0]*0.5);
  //f2->FixParamter(0,1.25448e+02);
  //f3->FixParameter(0,1.05448e+02*0.70);
  //f3->SetParLimits(0,1.05448e+02*0.52,1.05448e+02*0.55);
  //f3->SetParLimits(0,1.05448e+02*0.4,1.05448e+02*0.6);
  //f3->SetParLimits(0,1.14645e+03*0.05,1.14645e+03*0.06);
  //f3->SetParLimits(0,1.14645e+03*0.07,1.14645e+03*0.1);
  //f3->SetParLimits(1,0.0,0.1);
  //f3->SetParLimits(2,0.15,0.2);
  f3->FixParameter(1,param[1]);
  f3->FixParameter(2,param[2]);
  f3->FixParameter(3,param[3]);
  //f3->FixParameter(3,1.9);
  //f3->SetParLimits(0,param[0]/3.0,param[0]/1.5);
  f3->SetRange(1.1,1.1,1.4,1.4);
  IMnpim_IMnpip_dE_wK0_woSid_n_3->Fit("f3","R","");
  //f3->Draw("cont2 same");
   
  TH2D* f3hist = (TH2D*)f3->GetHistogram();
  IMnpim_IMnpip_dE_wK0_woSid_n_3->Print("base");
  f3hist->Print("base");
  f3hist->SetLineColor(2);
  f3hist->SetFillColor(0);
  cK0fit->cd(1);
  TH1D* IMnpip_3 = (TH1D*)IMnpim_IMnpip_dE_wK0_woSid_n_3->ProjectionX("IMnpip_3");
  IMnpip_3->Draw("HE");
  f3hist->ProjectionX()->Draw("HISTsame");

  cK0fit->cd(4);
  TH1D* IMnpim_3 = (TH1D*)IMnpim_IMnpip_dE_wK0_woSid_n_3->ProjectionY("IMnpim_3");
  IMnpim_3->Draw("HE");
  f3hist->ProjectionY()->Draw("HISTsame");
   
  double paramf3[nparfit];
  f3->GetParameters(paramf3);

  TF2 *f3wide = new TF2("f3wide",K0fit2dNoconvert,xmin,xmax,ymin,ymax,nparfit);
  f3wide->SetNpx(nbinsX);//use same nbin to compare the projection
  f3wide->SetNpy(nbinsY);//use same nbin to compare the projection
  f3wide->SetParameters(paramf3);
  
  TF2 *f3wide_nocut = new TF2("f3wide_nocut",K0fit2dNoconvert,xmin,xmax,ymin,ymax,nparfit);
  f3wide_nocut->SetNpx(nbinsX);//use same nbin to compare the projection
  f3wide_nocut->SetNpy(nbinsY);//use same nbin to compare the projection
  f3wide_nocut->SetParameters(paramf3); 

  TH2D* f3widehist = (TH2D*)f3wide->GetHistogram();
  TH2D* f3widehist_nocut = (TH2D*)f3wide->GetHistogram();
  f3widehist->Print("base");
  
  const int SpbinMIN = IMnpim_IMnpip_dE_wK0_woSid_n_3->GetXaxis()->FindBin(anacuts::Sigmap_MIN_wide);
  const int SpbinMAX = IMnpim_IMnpip_dE_wK0_woSid_n_3->GetXaxis()->FindBin(anacuts::Sigmap_MAX_wide);
  const int SmbinMIN = IMnpim_IMnpip_dE_wK0_woSid_n_3->GetYaxis()->FindBin(anacuts::Sigmam_MIN_wide);
  const int SmbinMAX = IMnpim_IMnpip_dE_wK0_woSid_n_3->GetYaxis()->FindBin(anacuts::Sigmam_MAX_wide);
  
  //removing Sp OR Sm
  //idea:  also remove negative bin of IMnpim_IMnpip histogram ? 
  //-->Conclusion: Non-sense
  for(int ixbin=0;ixbin<=nbinsX;ixbin++){
    for(int iybin=0;iybin<=nbinsY;iybin++){
      
      //double cont =  IMnpim_IMnpip_dE_wK0_woSid_n_3->GetBinContent(ixbin,iybin);
      //if(cont<0.001){
      //  f3widehist->SetBinContent(ixbin,iybin,0);
      //  f3widehist->SetBinError(ixbin,iybin,0);
      //}
      if(SpbinMIN <= ixbin && ixbin<=SpbinMAX){
        f3widehist->SetBinContent(ixbin,iybin,0);
        f3widehist->SetBinError(ixbin,iybin,0);
      }
      if(SmbinMIN <= iybin && iybin<=SmbinMAX){
        f3widehist->SetBinContent(ixbin,iybin,0);
        f3widehist->SetBinError(ixbin,iybin,0);
      }

    }
  }
  auto *ctemp = new TCanvas("ctemp","ctemp",800,800);
  ctemp->Divide(2,2);
  ctemp->cd(3);
  //f3widehist->Draw("colz");
  TH2D* hf3widehist_nocut = (TH2D*)f3wide_nocut->GetHistogram();
  hf3widehist_nocut->SetFillColor(0);
  hf3widehist_nocut->Draw("colz");

  ctemp->cd(1);
  //TH2D* hf3wide = (TH2D*)f3wide->GetHistogram();
  hf3widehist_nocut->ProjectionX("hf3wide_px")->Draw("HIST");

  ctemp->cd(4);
  //hf3wide->ProjectionY("hf3wide_py")->Draw("HIST");
  hf3widehist_nocut->ProjectionY("hf3wide_py")->Draw("HIST");

  TCanvas *c3wide = new TCanvas("c3wide","c3wide",800,800);
  c3wide->Divide(2,2);
  c3wide->cd(3);
  IMnpim_IMnpip_dE_wK0_woSid_n_3->Draw("colz");
  f3widehist->SetFillColor(0);
  f3widehist->Draw("cont2same");

  c3wide->cd(1);
  IMnpip_3->Draw("HE");
  TH1D* f3widehist_px = (TH1D*)f3widehist->ProjectionX("f3widehist_px");
  f3widehist_px->Draw("HISTsame");

  c3wide->cd(4);
  IMnpim_3->Draw("HE");
  TH1D* f3widehist_py = (TH1D*)f3widehist->ProjectionY("f3widehist_py");
  f3widehist_py->Draw("HISTsame");
  
  TH2D* IMnpim_IMnpip_dE_wK0_woSid_n_3_inter = (TH2D*)IMnpim_IMnpip_dE_wK0_woSid_n_1->Clone("IMnpim_IMnpip_dE_wK0_woSid_n_3_inter");
  TH2D* h2K0inter_3 = (TH2D*)IMnpim_IMnpip_dE_wK0_woSid_n_3_inter->Clone("h2K0inter_3");
  h2K0inter_3->Reset();
  auto cinter_3 = new TCanvas("cinter_3","cinter_3",800,800);
  cinter_3->Divide(2,2);
  for(int ixbin=0;ixbin<=nbinsX;ixbin++){
    for(int iybin=0;iybin<=nbinsY;iybin++){
      double xcent = IMnpim_IMnpip_dE_wK0_woSid_n_3_inter->GetXaxis()->GetBinCenter(ixbin);
      double ycent = IMnpim_IMnpip_dE_wK0_woSid_n_3_inter->GetYaxis()->GetBinCenter(iybin);
      double cont = f3wide_nocut->Eval(xcent,ycent);   
      double xx = 1.0/sqrt(2.0)*(xcent-ycent);
      double yy = 1.0/sqrt(2.0)*(xcent+ycent);
      double yy2 = yy-(sqrt(6.76*xx*xx+2.765)-1.0);//slightly tuned from original value by looking final fitting result
      if(yy2<0.98 ) continue;
      if(SpbinMIN <= ixbin && ixbin<=SpbinMAX){
        //std::cout << xcent << " " << ycent << "  " << cont << std::endl;
        IMnpim_IMnpip_dE_wK0_woSid_n_3_inter->SetBinContent(ixbin,iybin,cont);
      }
      if(SmbinMIN <= iybin && iybin<=SmbinMAX){
        IMnpim_IMnpip_dE_wK0_woSid_n_3_inter->SetBinContent(ixbin,iybin,cont);
      }
    }
  }
  
  cinter_3->cd(3);
  IMnpim_IMnpip_dE_wK0_woSid_n_3_inter->Rebin2D(4,4);//+/- 2sigma binning
  IMnpim_IMnpip_dE_wK0_woSid_n_3_inter->Draw("colz");

  cinter_3->cd(1);
  auto *IMnpip_3_inter = (TH1D*)IMnpim_IMnpip_dE_wK0_woSid_n_3_inter->ProjectionX("IMnpip_3_inter");
  IMnpip_3_inter->Draw("HE");
  

  cinter_3->cd(4);
  auto *IMnpim_3_inter = (TH1D*)IMnpim_IMnpip_dE_wK0_woSid_n_3_inter->ProjectionY("IMnpim_3_inter");
  IMnpim_3_inter->Draw("HE");
  //auto IMnpim_3_inter->Draw("HE");
  


  //next step
  //subtract K0 and solve Sp/Sm overlap region

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
  std::cout << "Sigma m mass center " << anacuts::Sigmam_center << std::endl;
  std::cout << "SmbinStart low Edge " << IMnpim_IMnpip_dE_wK0orwSid_n_K0sub->GetYaxis()->GetBinLowEdge(SmbinStart) << std::endl;
  std::cout << "Sigmam_MIN " << anacuts::Sigmam_MIN << std::endl;
  std::cout << "Smbin width " << IMnpim_IMnpip_dE_wK0orwSid_n_K0sub->GetYaxis()->GetBinWidth(SmbinStart) << std::endl;
  std::cout << "Sigmam_sigma " << anacuts::Sigmam_sigma << std::endl;
  const int SmbinEnd = IMnpim_IMnpip_dE_wK0orwSid_n_K0sub->GetYaxis()->FindBin(anacuts::Sigmam_MAX);
  std::cout << "SmbinEnd low Edge " << IMnpim_IMnpip_dE_wK0orwSid_n_K0sub->GetYaxis()->GetBinLowEdge(SmbinEnd) << std::endl;
  std::cout << "SmbinEnd low Edge+1 " << IMnpim_IMnpip_dE_wK0orwSid_n_K0sub->GetYaxis()->GetBinLowEdge(SmbinEnd+1) << std::endl;
  std::cout << "Sigmam_MAX " << anacuts::Sigmam_MAX << std::endl;
  const int Spbin = IMnpim_IMnpip_dE_wK0orwSid_n_K0sub->GetXaxis()->FindBin(anacuts::Sigmap_center);
  TH1D* IMnpip_K0sub = (TH1D*)IMnpim_IMnpip_dE_wK0orwSid_n_K0sub->ProjectionX("IMnpip_K0sub",SmbinEnd,SmbinEnd);
  TH1D* IMnpip_K0sub_Sp = (TH1D*)IMnpim_IMnpip_dE_wK0orwSid_n_K0sub->ProjectionX("IMnpip_K0sub_Sp",SmbinEnd,SmbinEnd);
  IMnpip_K0sub->Draw("HE");
  IMnpip_K0sub_Sp->SetFillColor(2);
  IMnpip_K0sub_Sp->GetXaxis()->SetRange(Spbin,Spbin);
  IMnpip_K0sub_Sp->SetFillStyle(3002);
  IMnpip_K0sub_Sp->Draw("HEsame");

  cwSid_n_K0sub->cd(4);
  const int SpbinStart = IMnpim_IMnpip_dE_wK0orwSid_n_K0sub->GetXaxis()->FindBin(anacuts::Sigmap_MIN);
  std::cout << std::endl;
  std::cout << "Sigma p mass center " << anacuts::Sigmap_center << std::endl;
  std::cout << "SpbinStart low Edge " << IMnpim_IMnpip_dE_wK0orwSid_n_K0sub->GetXaxis()->GetBinLowEdge(SpbinStart) << std::endl;
  std::cout << "Sigmap_MIN " << anacuts::Sigmap_MIN << std::endl;
  const int SpbinEnd = IMnpim_IMnpip_dE_wK0orwSid_n_K0sub->GetXaxis()->FindBin(anacuts::Sigmap_MAX);
  std::cout << "SpbinEnd low Edge " << IMnpim_IMnpip_dE_wK0orwSid_n_K0sub->GetXaxis()->GetBinLowEdge(SpbinEnd) << std::endl;
  std::cout << "SpbinEnd low Edge+1 " << IMnpim_IMnpip_dE_wK0orwSid_n_K0sub->GetXaxis()->GetBinLowEdge(SpbinEnd+1) << std::endl;
  std::cout << "Sigmap_MAX " << anacuts::Sigmap_MAX << std::endl;
  const int Smbin = IMnpim_IMnpip_dE_wK0orwSid_n_K0sub->GetYaxis()->FindBin(anacuts::Sigmam_center);
  TH1D* IMnpim_K0sub = (TH1D*)IMnpim_IMnpip_dE_wK0orwSid_n_K0sub->ProjectionY("IMnpim_K0sub",SpbinEnd,SpbinEnd);
  IMnpim_K0sub->Draw("HE");
  TH1D* IMnpim_K0sub_Sm = (TH1D*)IMnpim_IMnpip_dE_wK0orwSid_n_K0sub->ProjectionY("IMnpim_K0sub_Sm",SpbinEnd,SpbinEnd);
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
