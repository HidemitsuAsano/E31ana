//fit macro to decompose K0/S+/S- final states
//step 1. determine the shape of K0 distribution by deformed IM(npi-) vs IM(npi+) distribution under the condition which negative bin is set to 0.
//step 2. refit on the original IM(npi-) vs IM(npi+) without zero suppression, 
//this is necesarry because the scaling parameter due to the change of the binning is umbiguous.
//
//TODO : decompose K0 & S+ & S- three states overlap 
#include <TF2.h>
#include <TH2.h>
#include <TMath.h>
#include "anacuts.h"

const double sysscale=0.2;//systematic error for K0 inter ratio
const bool NoRebin = false;
const int RebinFactor = 4;
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
  const Double_t a  = 0.004;
  
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
  const Double_t thend = 1.95;
  const Double_t aend = 0.05;
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


const int Version = 245;

void Fit2DK0(const int qcut=1,const int dEcut=2,const int sysud=0)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TFile *fr = NULL;
  //Because the statistics is limited, we just divide data into q<0.35 and q>0.35.
  if(qcut==1){//qlo
    //fr = TFile::Open("evanaIMpisigma_npippim_v206_out_iso_qlo_sub.root","READ");
    fr = TFile::Open(Form("evanaIMpisigma_npippim_v%d_out_dE%d_iso_qlo_nostop_sys%d_sub.root",Version,dEcut,sysud),"READ");
  }else if(qcut==2){//qhi
    //fr = TFile::Open("evanaIMpisigma_npippim_v206_out_iso_qhi_sub.root","READ");
    fr = TFile::Open(Form("evanaIMpisigma_npippim_v%d_out_dE%d_iso_qhi_nostop_sys%d_sub.root",Version,dEcut,sysud),"READ");
  }else{
    std::cout << "no file" << std::endl;
    return;
  }
  fr->Print();


  //
  TH2F* IMnpim_IMnpip_dE_wK0_woSid_n_1 = (TH2F*)fr->Get("IMnpim_IMnpip_dE_wK0_woSid_n");
  //45 degeree rotation of IMnpi- vs IMnpi+
  TH2F* IMnpim_IMnpip_dE_wK0_woSid_n_45rot = (TH2F*)fr->Get("IMnpim_IMnpip_dE_wK0_woSid_n_45rot");
  //deform the distribution so as to align the edge of x-axis
  TH2F* IMnpim_IMnpip_dE_wK0_woSid_n_45rot3 = (TH2F*)fr->Get("IMnpim_IMnpip_dE_wK0_woSid_n_45rot3");
  
  
  TCanvas *cIMnpim_IMnpip_dE_wK0_woSid_n = new TCanvas("cIMnpim_IMnpip_dE_wK0_woSid_n","cIMnpim_IMnpip_dE_wK0_woSid_n",800,800);
  cIMnpim_IMnpip_dE_wK0_woSid_n->Divide(2,2,0.,0.);
  cIMnpim_IMnpip_dE_wK0_woSid_n->cd(3);
  //IMnpim_IMnpip_dE_wK0_woSid_n_1->Rebin2D(4,4);
  //IMnpim_IMnpip_dE_wK0_woSid_n_1->RebinX(4);
  //IMnpim_IMnpip_dE_wK0_woSid_n_1->RebinY(4);
  IMnpim_IMnpip_dE_wK0_woSid_n_1->SetMinimum(0);
  IMnpim_IMnpip_dE_wK0_woSid_n_1->Draw("colz");
  
  cIMnpim_IMnpip_dE_wK0_woSid_n->cd(1);
  TH1D *IMnpip_wK0_woSid_n = (TH1D*)IMnpim_IMnpip_dE_wK0_woSid_n_1->ProjectionX("IMnpip_wK0_woSid_n");
  IMnpip_wK0_woSid_n->Draw("HE");
  
  cIMnpim_IMnpip_dE_wK0_woSid_n->cd(4);
  TH1D *IMnpim_wK0_woSid_n = (TH1D*)IMnpim_IMnpip_dE_wK0_woSid_n_1->ProjectionY("IMnpim_wK0_woSid_n");
  IMnpim_wK0_woSid_n->Draw("HE");
   

  TCanvas *cIMnpim_IMnpip_dE_wK0_woSid_n_45rot = new TCanvas("cIMnpim_IMnpip_dE_wK0_woSid_n_45rot","cIMnpim_IMnpip_dE_wK0_woSid_n_45rot",800,800);
  cIMnpim_IMnpip_dE_wK0_woSid_n_45rot->Divide(2,2,0,0);
  cIMnpim_IMnpip_dE_wK0_woSid_n_45rot->cd(3);
  IMnpim_IMnpip_dE_wK0_woSid_n_45rot->SetMinimum(0);
  IMnpim_IMnpip_dE_wK0_woSid_n_45rot->Draw("colz");
  
  cIMnpim_IMnpip_dE_wK0_woSid_n_45rot->cd(1);
  IMnpim_IMnpip_dE_wK0_woSid_n_45rot->ProjectionX()->Draw();
  
  cIMnpim_IMnpip_dE_wK0_woSid_n_45rot->cd(4);
  IMnpim_IMnpip_dE_wK0_woSid_n_45rot->ProjectionY()->Draw();

  TCanvas *cIMnpim_IMnpip_dE_wK0_woSid_n_45rot3 = new TCanvas("cIMnpim_IMnpip_dE_wK0_woSid_n_45rot3","cIMnpim_IMnpip_dE_wK0_woSid_n_45rot3",800,800);
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
  TH2D *IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_0suppress = (TH2D*)IMnpim_IMnpip_dE_wK0_woSid_n_45rot3->Clone("IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_0suppress");
  for(int ix=0;ix<IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_0suppress->GetNbinsX();ix++){
    for(int iy=0;iy<IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_0suppress->GetNbinsY();iy++){
      double cont = IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_0suppress->GetBinContent(ix,iy);
      //Is this OK ? 
      //-> without this treatment,
      //negative bin with very small statistic error which comes from the mixed events
      //strongly constraints the entire shape of fitting.
      if(cont<0.001){
        IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_0suppress->SetBinContent(ix,iy,0);
        IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_0suppress->SetBinError(ix,iy,0);
      }
    }
  }
  //IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_0suppress->GetYaxis()->SetRangeUser(1.0,1.5);
  IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_0suppress->Draw("colz");
  
  //fitting on deformed K0 distribution
  TF2 *f2 = new TF2("f2",K0fit2d,-0.5,0.5,0.9,1.5,nparfit);
  //par0 : scaling factor
  //par1 : x,gaus mean
  //par2 : x,gaus sigma
  //par3 : y,exp slope
  
  //xmin,ymin,xmax,ymax
  f2->SetRange(-0.3,0.9,0.3,1.4); 
  //f2->SetParameters(8.0e9,0.005,0.16,-15.2);
  f2->SetParameters(3.0e5,0.005,0.20,1.9);
  if(qcut==2)f2->SetParLimits(0,150,5.5e12);
  //else if(qcut==1)f2->SetParLimits(0,3.57224e+04*1.8,3.57224e+04*2.0);
  else if(qcut==1)f2->SetParLimits(0,2.19399e+04*1.85,2.19399e+04*2.55 );
  //f2->FixParameter(0,2.0e5);
  if(qcut==2)f2->SetParLimits(1,-0.005,0.005);
  else if(qcut==1)f2->SetParLimits(1,-0.005,0.005);
  //f2->FixParameter(1,0.005);
  //f2->SetParLimits(2,0.20,0.25);
  if(qcut==2)f2->FixParameter(2,0.20);
  //f2->SetParameter(3,0.5);
  //f2->SetParLimits(3,1.66,3.00);
  //f2->SetNpx(100);//=NbinsX of rot3 histogram
  //f2->SetNpy(120);//=NbinsY of rot3 histogram
  if(qcut==2)f2->FixParameter(3,1.9);
  else if(qcut==1)f2->FixParameter(3,2.9);
  //Log likehoog is to be used when the histogram represents counts,
  //however this data is not the case because the BG is subtracted by event. (?)
  //Ignoring this because statistic is not enough and chi2 fitting is not stable.
  IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_0suppress->Fit("f2","RL","");
  IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_0suppress->Print("base");
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
  for(int ix=0;ix<=IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_0suppress->GetNbinsX();ix++){
    for(int iy=0;iy<=IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_0suppress->GetNbinsY();iy++){
      double cont = IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_0suppress->GetBinContent(ix,iy);
      if((cont)<0.00000001) hf2wide->SetBinContent(ix,iy,0);
      double xcen=  IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_0suppress->GetXaxis()->GetBinCenter(ix);
      double ycen=  IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_0suppress->GetYaxis()->GetBinCenter(iy);
      //remove edge area
      if( (fabs(xcen)<0.02) && (ycen < 1.02)){
        hf2wide->SetBinContent(ix,iy,0);
        hf2wide->SetBinError(ix,iy,0);
        IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_0suppress->SetBinContent(ix,iy,0);
        IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_0suppress->SetBinError(ix,iy,0);
      }
    }
  }

  IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_0suppress->Print("base");
  hf2->Print("base");
  cfitcomp->cd(1);
  TH1D* pxrot3 = (TH1D*)IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_0suppress->ProjectionX("pxrot3");
  pxrot3->Draw("HE");
  hf2->SetFillColor(0);
  hf2->ProjectionX()->Draw("HISTsame");
  
  cfitcomp->cd(4);
  TH1D* pyrot3 = (TH1D*)IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_0suppress->ProjectionY("pyrot3");
  pyrot3->Draw("HE");
  hf2->ProjectionY()->Draw("HISTsame");
  
  //subtract and plot wide fit
  TCanvas *cZeroSuppress = new TCanvas("cZeroSuppress","cZerosuppress",800,800);
  cZeroSuppress->Divide(2,2);
  cZeroSuppress->cd(3);
  IMnpim_IMnpip_dE_wK0_woSid_n_45rot3_0suppress->Draw("colz");
  hf2wide->Draw("box same");

  cZeroSuppress->cd(1);
  pxrot3->Draw("HE");
  hf2wide->SetFillColor(0);
  hf2wide->ProjectionX()->Draw("HISTsame");

  cZeroSuppress->cd(4);
  pyrot3->Draw("HE");
  hf2wide->ProjectionY()->Draw("HISTsame");
   
  TCanvas *ctest = new TCanvas("ctest","ctest",1600,800);
  ctest->Divide(2,1);
  ctest->cd(1);
  hf2wide->Draw("colz");
  ctest->cd(2);
  hf2wide_nosub->SetTitle("hf2wide_nosub");
  hf2wide_nosub->Draw("colz");

  
  //no rotation plots
  TH2D *IMnpim_IMnpip_dE_wK0_woSid_n_2 = (TH2D*)IMnpim_IMnpip_dE_wK0_woSid_n_1->Clone("IMnpim_IMnpip_dE_wK0_woSid_n_2");
  IMnpim_IMnpip_dE_wK0_woSid_n_2->SetName("IMnpim_IMnpip_dE_wK0_woSid_n_2");
  TCanvas *cinter = new TCanvas("cinter","cinter",800,800);
  cinter->Divide(2,2);
  //Interpolation histogram in the Sigma+ and Simga- mass region
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
          double yy2 = yy-(sqrt(6.76*xx*xx+2.765)-1.0);
          double evalK0 = f2wide->Eval(xx,yy2);
          double scale=0.28;//this scaling factor is arbitrary and determined by eye->this is OK, because final parametes are determined by next step
          evalK0 *= scale;
          if(yy2>0.98 && xcent < 1.7 && ycent<1.7){
            IMnpim_IMnpip_dE_wK0_woSid_n_2->SetBinContent(ix,iy,evalK0);
            h2K0inter->SetBinContent(ix,iy,evalK0);
          }
        }
      }//cont
    }
  }
  std::cout  << "hf2wide integral " << hf2wide->Integral() << std::endl;
  std::cout << " K0inter         " << h2K0inter->Integral() << std::endl;

  if(!NoRebin)IMnpim_IMnpip_dE_wK0_woSid_n_2->Rebin2D(RebinFactor,RebinFactor);
  if(!NoRebin)h2K0inter->Rebin2D(RebinFactor,RebinFactor);
  cinter->cd(3);
  IMnpim_IMnpip_dE_wK0_woSid_n_2->Draw("colz");

  cinter->cd(1);
  TH1D* IMnpip_wK0_woSid= (TH1D*)IMnpim_IMnpip_dE_wK0_woSid_n_2->ProjectionX("IMnpip_wK0_woSid");
  IMnpip_wK0_woSid->Draw("EH");
  TH1D* K0interpx = (TH1D*)h2K0inter->ProjectionX("K0interpx");
  K0interpx->SetFillColor(2);
  K0interpx->Draw("HISTsame");
  cinter->cd(4);
  TH1D* IMnpim_wK0_woSid= (TH1D*)IMnpim_IMnpip_dE_wK0_woSid_n_2->ProjectionY("IMnpim_wK0_woSid");
  IMnpim_wK0_woSid->Draw("EH");
  TH1D* K0interpy = (TH1D*)h2K0inter->ProjectionY("K0interpy");
  K0interpy->SetFillColor(2);
  K0interpy->Draw("HISTsame");
  

  //In the steps up to this point, the scaling parameter were determined by eye. 
  //From here, the scaling parameter will be determined using fits.
  TH2D *IMnpim_IMnpip_dE_wK0_woSid_n_3 = (TH2D*)IMnpim_IMnpip_dE_wK0_woSid_n_1->Clone("IMnpim_IMnpip_dE_wK0_woSid_n_3");
  TCanvas *cK0fit = new TCanvas("cK0fit","cK0fit",800,800);
  cK0fit->Divide(2,2);
  cK0fit->cd(3);
  for(int ix=0;ix<IMnpim_IMnpip_dE_wK0_woSid_n_3->GetNbinsX();ix++){
    for(int iy=0;iy<IMnpim_IMnpip_dE_wK0_woSid_n_3->GetNbinsY();iy++){
      double cont = IMnpim_IMnpip_dE_wK0_woSid_n_3->GetBinContent(ix,iy);
      double xcent = IMnpim_IMnpip_dE_wK0_woSid_n_3->GetXaxis()->GetBinCenter(ix);
      double ycent = IMnpim_IMnpip_dE_wK0_woSid_n_3->GetYaxis()->GetBinCenter(iy);
      //remove edge region which K0nn cannot go.
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
  if(qcut==2)f3->SetParLimits(0,param[0]*0.25,param[0]*0.30);
  //else if(qcut==1)f3->SetParLimits(0,param[0]*0.35,param[0]*0.9);
  else if(qcut==1)f3->SetParLimits(0,param[0]*0.3,param[0]*0.6);
  f3->FixParameter(1,param[1]);
  f3->FixParameter(2,param[2]);
  f3->FixParameter(3,param[3]);
  //f3->FixParameter(3,1.9);
  //f3->SetParLimits(0,param[0]/3.0,param[0]/1.5);
  f3->SetRange(1.1,1.1,1.4,1.4);
  //f3->SetRange(1.1,1.1,1.45,1.45);
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
  
  //Sigma region will be removed later
  TF2 *f3wide = new TF2("f3wide",K0fit2dNoconvert,xmin,xmax,ymin,ymax,nparfit);
  f3wide->SetNpx(nbinsX);//use same nbin to compare the projection
  f3wide->SetNpy(nbinsY);//use same nbin to compare the projection
  f3wide->SetParameters(paramf3);
  
  TF2 *f3wide_nocut = new TF2("f3wide_nocut",K0fit2dNoconvert,xmin,xmax,ymin,ymax,nparfit);
  f3wide_nocut->SetNpx(nbinsX);//use same nbin to compare the projection
  f3wide_nocut->SetNpy(nbinsY);//use same nbin to compare the projection
  f3wide_nocut->SetParameters(paramf3); 
  
  TF2 *f3wide_nocut_sysup = new TF2("f3wide_nocut_sysup",K0fit2dNoconvert,xmin,xmax,ymin,ymax,nparfit);
  f3wide_nocut_sysup->SetNpx(nbinsX);//use same nbin to compare the projection
  f3wide_nocut_sysup->SetNpy(nbinsY);//use same nbin to compare the projection
  f3wide_nocut_sysup->SetParameters(paramf3); 
  f3wide_nocut_sysup->SetParameter(0,paramf3[0]*(1.0+sysscale)); 
  
  TF2 *f3wide_nocut_sysdown = new TF2("f3wide_nocut_sysdown",K0fit2dNoconvert,xmin,xmax,ymin,ymax,nparfit);
  f3wide_nocut_sysdown->SetNpx(nbinsX);//use same nbin to compare the projection
  f3wide_nocut_sysdown->SetNpy(nbinsY);//use same nbin to compare the projection
  f3wide_nocut_sysdown->SetParameters(paramf3); 
  f3wide_nocut_sysdown->SetParameter(0,paramf3[0]*(1.0-sysscale)); 
  
  f3wide_nocut->Print();
  f3wide_nocut_sysup->Print();
  f3wide_nocut_sysdown->Print();

  TH2D* f3widehist = (TH2D*)f3wide->GetHistogram();
  TH2D* f3widehist_nocut = (TH2D*)f3wide->GetHistogram();
  f3widehist->Print("base");
  
  //wide (= 3sigma cut)
  const int SpbinMIN_wide = IMnpim_IMnpip_dE_wK0_woSid_n_3->GetXaxis()->FindBin(anacuts::Sigmap_MIN_wide);
  const int SpbinMAX_wide = IMnpim_IMnpip_dE_wK0_woSid_n_3->GetXaxis()->FindBin(anacuts::Sigmap_MAX_wide);
  const int SmbinMIN_wide = IMnpim_IMnpip_dE_wK0_woSid_n_3->GetYaxis()->FindBin(anacuts::Sigmam_MIN_wide);
  const int SmbinMAX_wide = IMnpim_IMnpip_dE_wK0_woSid_n_3->GetYaxis()->FindBin(anacuts::Sigmam_MAX_wide);
  
  //removing Sp OR Sm
  //idea:  also remove negative bin of IMnpim_IMnpip histogram ? 
  //-->Conclusion: Non-sense
  for(int ixbin=0;ixbin<=nbinsX;ixbin++){
    for(int iybin=0;iybin<=nbinsY;iybin++){
      //double cont =  IMnpim_IMnpip_dE_wK0_woSid_n_3->GetBinContent(ixbin,iybin);
      //if(cont<0.001){
        //f3widehist->SetBinContent(ixbin,iybin,0);
        //f3widehist->SetBinError(ixbin,iybin,0);
      //}
      //3 sigma cut
      if(SpbinMIN_wide <= ixbin && ixbin<=SpbinMAX_wide){
        f3widehist->SetBinContent(ixbin,iybin,0);
        f3widehist->SetBinError(ixbin,iybin,0);
      }
      //3 sigma cut
      if(SmbinMIN_wide <= iybin && iybin<=SmbinMAX_wide){
        f3widehist->SetBinContent(ixbin,iybin,0);
        f3widehist->SetBinError(ixbin,iybin,0);
      }
    }
  }

  TCanvas *ctemp = new TCanvas("ctemp","ctemp",800,800);
  ctemp->Divide(2,2);
  ctemp->cd(3);
  //f3widehist->Draw("colz");
  TH2D* hf3widehist_nocut = (TH2D*)f3wide_nocut->GetHistogram();
  hf3widehist_nocut->SetFillColor(0);
  hf3widehist_nocut->Draw("colz");

  TH2D* hf3widehist_nocut_sysup = (TH2D*)f3wide_nocut_sysup->GetHistogram();
  TH2D* hf3widehist_nocut_sysdown = (TH2D*)f3wide_nocut_sysdown->GetHistogram();
  ctemp->cd(1);
  //TH2D* hf3wide = (TH2D*)f3wide->GetHistogram();
  hf3widehist_nocut->ProjectionX("hf3wide_px")->Draw("HIST");
  hf3widehist_nocut_sysup->ProjectionX("hf3wide_sysup_px")->Draw("HISTsame");
  hf3widehist_nocut_sysdown->ProjectionX("hf3wide_sysdown_px")->Draw("HISTsame");

  ctemp->cd(4);
  //hf3wide->ProjectionY("hf3wide_py")->Draw("HIST");
  hf3widehist_nocut->ProjectionY("hf3wide_py")->Draw("HIST");
  hf3widehist_nocut_sysup->ProjectionY("hf3wide_sysup_py")->Draw("HISTsame");
  hf3widehist_nocut_sysdown->ProjectionY("hf3wide_sysdown_py")->Draw("HISTsame");

  TCanvas *cNoZeroSuppress = new TCanvas("cNoZeroSuppress","cNoZeroSuppress",800,800);
  cNoZeroSuppress->Divide(2,2);
  cNoZeroSuppress->cd(3);
  IMnpim_IMnpip_dE_wK0_woSid_n_3->Draw("colz");
  f3widehist->SetFillColor(0);
  f3widehist->Draw("cont2same");
  TH2D* f3widehist_sysup = (TH2D*)f3widehist->Clone("f3widehist_sysup");
  f3widehist_sysup->Scale(1.0+sysscale);
  f3widehist_sysup->SetLineColor(3);
  TH2D* f3widehist_sysdown = (TH2D*)f3widehist->Clone("f3widehist_sysdown");
  f3widehist_sysdown->Scale(1.0-sysscale);
  f3widehist_sysdown->SetLineColor(3);

  cNoZeroSuppress->cd(1);
  IMnpip_3->Draw("HE");
  TH1D* f3widehist_px = (TH1D*)f3widehist->ProjectionX("f3widehist_px");
  TH1D* f3widehist_sysup_px = (TH1D*)f3widehist_sysup->ProjectionX("f3widehist_sysup_px");
  TH1D* f3widehist_sysdown_px = (TH1D*)f3widehist_sysdown->ProjectionX("f3widehist_sysdown_px");
  f3widehist_px->Draw("HISTsame");
  f3widehist_sysup_px->Draw("HISTsame");
  f3widehist_sysdown_px->Draw("HISTsame");

  cNoZeroSuppress->cd(4);
  IMnpim_3->Draw("HE");
  TH1D* f3widehist_py = (TH1D*)f3widehist->ProjectionY("f3widehist_py");
  TH1D* f3widehist_sysup_py = (TH1D*)f3widehist_sysup->ProjectionY("f3widehist_sysup_py");
  TH1D* f3widehist_sysdown_py = (TH1D*)f3widehist_sysdown->ProjectionY("f3widehist_sysdown_py");
  f3widehist_py->Draw("HISTsame");
  f3widehist_sysup_py->Draw("HISTsame");
  f3widehist_sysdown_py->Draw("HISTsame");
   
  cNoZeroSuppress->cd(2);
  TPaveText *pt = new TPaveText(.05,.05,.95,.7);
  pt->AddText(Form("Negative bin is NOT removed."));
  pt->AddText(Form("(x<1.18 and y <1.18) data are removed for comparison "));
  pt->AddText(Form("Final scaling parameter is determined by this fit. "));
  pt->Draw();


  TH2D* IMnpim_IMnpip_dE_wK0_woSid_n_3_inter = (TH2D*)IMnpim_IMnpip_dE_wK0_woSid_n_1->Clone("IMnpim_IMnpip_dE_wK0_woSid_n_3_inter");
  TH2D* IMnpim_IMnpip_dE_wK0_woSid_n_3_inter_sysup = (TH2D*)IMnpim_IMnpip_dE_wK0_woSid_n_1->Clone("IMnpim_IMnpip_dE_wK0_woSid_n_3_inter_sysup");
  TH2D* IMnpim_IMnpip_dE_wK0_woSid_n_3_inter_sysdown = (TH2D*)IMnpim_IMnpip_dE_wK0_woSid_n_1->Clone("IMnpim_IMnpip_dE_wK0_woSid_n_3_inter_sysdown");
  TH2D* h2K0inter_3 = (TH2D*)IMnpim_IMnpip_dE_wK0_woSid_n_3_inter->Clone("h2K0inter_3");
  TH2D* h2K0inter_3_sysup = (TH2D*)IMnpim_IMnpip_dE_wK0_woSid_n_3_inter_sysup->Clone("h2K0inter_3_sysup");
  TH2D* h2K0inter_3_sysdown = (TH2D*)IMnpim_IMnpip_dE_wK0_woSid_n_3_inter_sysdown->Clone("h2K0inter_3_sysdown");
  h2K0inter_3->Reset();
  h2K0inter_3_sysup->Reset();
  h2K0inter_3_sysdown->Reset();
  TCanvas* cinter_3 = new TCanvas("cinter_3","cinter_3",800,800);
  cinter_3->Divide(2,2);
  for(int ixbin=0;ixbin<=nbinsX;ixbin++){
    for(int iybin=0;iybin<=nbinsY;iybin++){
      double xcent = IMnpim_IMnpip_dE_wK0_woSid_n_3_inter->GetXaxis()->GetBinCenter(ixbin);
      double ycent = IMnpim_IMnpip_dE_wK0_woSid_n_3_inter->GetYaxis()->GetBinCenter(iybin);
      double cont = f3wide_nocut->Eval(xcent,ycent);   
      double contup = f3wide_nocut_sysup->Eval(xcent,ycent);   
      double contdown = f3wide_nocut_sysdown->Eval(xcent,ycent);   
      double xx = 1.0/sqrt(2.0)*(xcent-ycent);
      double yy = 1.0/sqrt(2.0)*(xcent+ycent);
      double yy2 = yy-(sqrt(6.76*xx*xx+2.765)-1.0);//
      if(yy2<0.98 ) continue;
      if(SpbinMIN_wide <= ixbin && ixbin<=SpbinMAX_wide){
        //std::cout << xcent << " " << ycent << "  " << cont << std::endl;
        IMnpim_IMnpip_dE_wK0_woSid_n_3_inter->SetBinContent(ixbin,iybin,cont);
        IMnpim_IMnpip_dE_wK0_woSid_n_3_inter_sysup->SetBinContent(ixbin,iybin,contup);
        IMnpim_IMnpip_dE_wK0_woSid_n_3_inter_sysdown->SetBinContent(ixbin,iybin,contdown);
        h2K0inter_3->SetBinContent(ixbin,iybin,cont);
        h2K0inter_3_sysup->SetBinContent(ixbin,iybin,contup);
        h2K0inter_3_sysdown->SetBinContent(ixbin,iybin,contdown);
      }
      if(SmbinMIN_wide <= iybin && iybin<=SmbinMAX_wide){
        IMnpim_IMnpip_dE_wK0_woSid_n_3_inter->SetBinContent(ixbin,iybin,cont);
        IMnpim_IMnpip_dE_wK0_woSid_n_3_inter_sysup->SetBinContent(ixbin,iybin,contup);
        IMnpim_IMnpip_dE_wK0_woSid_n_3_inter_sysdown->SetBinContent(ixbin,iybin,contdown);
        h2K0inter_3->SetBinContent(ixbin,iybin,cont);
        h2K0inter_3_sysup->SetBinContent(ixbin,iybin,contup);
        h2K0inter_3_sysdown->SetBinContent(ixbin,iybin,contdown);
      }
    }
  }
  
  cinter_3->cd(3);
  //result of before rebin
  TH2D* h2K0inter_3fine = (TH2D*)h2K0inter_3->Clone("h2K0inter_3fine");
  TH2D* h2K0inter_3fine_sysup = (TH2D*)h2K0inter_3_sysup->Clone("h2K0inter_3fine_sysup");
  TH2D* h2K0inter_3fine_sysdown = (TH2D*)h2K0inter_3_sysdown->Clone("h2K0inter_3fine_sysdown");
  IMnpim_IMnpip_dE_wK0_woSid_n_3_inter->Rebin2D(RebinFactor,RebinFactor);//+/- 2sigma binning
  IMnpim_IMnpip_dE_wK0_woSid_n_3_inter_sysup->Rebin2D(RebinFactor,RebinFactor);//+/- 2sigma binning
  IMnpim_IMnpip_dE_wK0_woSid_n_3_inter_sysdown->Rebin2D(RebinFactor,RebinFactor);//+/- 2sigma binning
  h2K0inter_3->Rebin2D(RebinFactor,RebinFactor);
  h2K0inter_3_sysup->Rebin2D(RebinFactor,RebinFactor);
  h2K0inter_3_sysdown->Rebin2D(RebinFactor,RebinFactor);
  IMnpim_IMnpip_dE_wK0_woSid_n_3_inter->Draw("colz");
  
  h2K0inter_3->SetFillColor(2);
  cinter_3->cd(1);
  //K0nn final
  TH1D *IMnpip_3_inter = (TH1D*)IMnpim_IMnpip_dE_wK0_woSid_n_3_inter->ProjectionX("IMnpip_3_inter");
  TH1D *IMnpip_3_inter_sysup = (TH1D*)IMnpim_IMnpip_dE_wK0_woSid_n_3_inter_sysup->ProjectionX("IMnpip_3_inter_sysup");
  TH1D *IMnpip_3_inter_sysdown = (TH1D*)IMnpim_IMnpip_dE_wK0_woSid_n_3_inter_sysdown->ProjectionX("IMnpip_3_inter_sysdown");
  IMnpip_3_inter->Draw("HE");
  h2K0inter_3_sysup->SetLineColor(3);
  h2K0inter_3_sysdown->SetLineColor(4);
  h2K0inter_3->ProjectionX()->Draw("HEsame");
  h2K0inter_3_sysup->ProjectionX()->Draw("HEsame");
  h2K0inter_3_sysdown->ProjectionX()->Draw("HEsame");

  cinter_3->cd(4);
  //K0nn final
  TH1D *IMnpim_3_inter = (TH1D*)IMnpim_IMnpip_dE_wK0_woSid_n_3_inter->ProjectionY("IMnpim_3_inter");
  TH1D *IMnpim_3_inter_sysup = (TH1D*)IMnpim_IMnpip_dE_wK0_woSid_n_3_inter_sysup->ProjectionY("IMnpim_3_inter_sysup");
  TH1D *IMnpim_3_inter_sysdown = (TH1D*)IMnpim_IMnpip_dE_wK0_woSid_n_3_inter_sysdown->ProjectionY("IMnpim_3_inter_sysdown");
  IMnpim_3_inter->Draw("HE");//estimated K0nn whole region
  h2K0inter_3->ProjectionY()->Draw("HEsame");//estimated K0nn on Sp or Sm region 
  h2K0inter_3_sysup->ProjectionY()->Draw("HEsame");
  h2K0inter_3_sysdown->ProjectionY()->Draw("HEsame");


  TH2F* IMnpim_IMnpip_dE_wK0orwSid_n = (TH2F*)fr->Get("IMnpim_IMnpip_dE_wK0orwSid_n");
  TH2F* IMnpim_IMnpip_dE_wK0orwSid_n_norebin = (TH2F*)fr->Get("IMnpim_IMnpip_dE_wK0orwSid_n");
  IMnpim_IMnpip_dE_wK0orwSid_n_norebin->SetName("IMnpim_IMnpip_dE_wK0orwSid_n_norebin");
  TH2F* IMnpim_IMnpip_dE_wK0_wSid_n_Sp = (TH2F*)fr->Get("IMnpim_IMnpip_dE_wK0_wSid_n_Sp");
  TH2F* IMnpim_IMnpip_dE_wK0_wSid_n_Sp_wbin1 = (TH2F*)fr->Get("IMnpim_IMnpip_dE_wK0_wSid_n_Sp_wbin1");
  TH2F* IMnpim_IMnpip_dE_wK0_wSid_n_Sp_wbin2 = (TH2F*)fr->Get("IMnpim_IMnpip_dE_wK0_wSid_n_Sp_wbin2");
  TH2F* IMnpim_IMnpip_dE_wK0_wSid_n_Sm = (TH2F*)fr->Get("IMnpim_IMnpip_dE_wK0_wSid_n_Sm");
  TH2F* IMnpim_IMnpip_dE_wK0_wSid_n_Sm_wbin1 = (TH2F*)fr->Get("IMnpim_IMnpip_dE_wK0_wSid_n_Sm_wbin1");
  TH2F* IMnpim_IMnpip_dE_wK0_wSid_n_Sm_wbin2 = (TH2F*)fr->Get("IMnpim_IMnpip_dE_wK0_wSid_n_Sm_wbin2");
  TH2F* IMnpim_IMnpip_dE_wK0orwSid_n_rebin = (TH2F*)IMnpim_IMnpip_dE_wK0orwSid_n->Clone("IMnpim_IMnpip_dE_wK0orwSid_n_rebin");
  IMnpim_IMnpip_dE_wK0orwSid_n_rebin->Rebin2D(RebinFactor,RebinFactor);
  TCanvas *ccomp = new TCanvas("ccomp","ccomp",1600,800);
  ccomp->Divide(2,1);
  ccomp->cd(1);
  h2K0inter_3->Draw("colz");
  ccomp->cd(2);
  //IMnpim_IMnpip_dE_wK0orwSid_n_norebin->SetMinimum(0);
  IMnpim_IMnpip_dE_wK0orwSid_n_norebin->Draw("colz");
  TCanvas *ccomp2 = new TCanvas("ccomp2","ccomp2",1600,800);
  ccomp2->Divide(2,1);
  ccomp2->cd(1);
  IMnpim_IMnpip_dE_wK0_wSid_n_Sp->Rebin2D(RebinFactor,RebinFactor);
  IMnpim_IMnpip_dE_wK0_wSid_n_Sp_wbin1->Rebin2D(RebinFactor,RebinFactor);
  IMnpim_IMnpip_dE_wK0_wSid_n_Sp_wbin2->Rebin2D(RebinFactor,RebinFactor);
  IMnpim_IMnpip_dE_wK0_wSid_n_Sp->Draw("colz");
  ccomp2->cd(2);
  IMnpim_IMnpip_dE_wK0_wSid_n_Sm->Rebin2D(RebinFactor,RebinFactor);
  IMnpim_IMnpip_dE_wK0_wSid_n_Sm_wbin1->Rebin2D(RebinFactor,RebinFactor);
  IMnpim_IMnpip_dE_wK0_wSid_n_Sm_wbin2->Rebin2D(RebinFactor,RebinFactor);
  IMnpim_IMnpip_dE_wK0_wSid_n_Sm->Draw("colz");
 
  //projection in Sigma+ region to IM(npi-)
  TCanvas *ccomp_Spregion = new TCanvas("ccomp_Spregion","ccomp_Spregion");
  IMnpim_IMnpip_dE_wK0_wSid_n_Sp->ProjectionY()->Draw();
  const int Spbinlow = h2K0inter_3->GetXaxis()->FindBin(anacuts::Sigmap_MIN);
  const int Spbinhi = h2K0inter_3->GetXaxis()->FindBin(anacuts::Sigmap_MAX);
  TH1D* IMnpimK0inter_Sp = (TH1D*)h2K0inter_3->ProjectionY("IMnpim_K0inter_Sp",Spbinlow,Spbinhi);
  TH1D* IMnpimK0inter_Sp_sysup = (TH1D*)h2K0inter_3_sysup->ProjectionY("IMnpim_K0inter_Sp_sysup",Spbinlow,Spbinhi);
  TH1D* IMnpimK0inter_Sp_sysdown = (TH1D*)h2K0inter_3_sysdown->ProjectionY("IMnpim_K0inter_Sp_sysdown",Spbinlow,Spbinhi);
  IMnpimK0inter_Sp->SetLineColor(2);
  IMnpimK0inter_Sp->Draw("Hsame");
  TH1D* wbin1_Sp = (TH1D*)IMnpim_IMnpip_dE_wK0_wSid_n_Sp_wbin1->ProjectionY();
  wbin1_Sp->SetLineColor(3);
  wbin1_Sp->Draw("same");
  TH1D* wbin2_Sp = (TH1D*)IMnpim_IMnpip_dE_wK0_wSid_n_Sp_wbin2->ProjectionY();
  wbin2_Sp->SetLineColor(4);
  wbin2_Sp->Draw("same");

  int spbinlow[2];  
  if(qcut==2){
    spbinlow[0]=6;
    spbinlow[1]=11;
  }else if(qcut==1){
    spbinlow[0]=6;
    spbinlow[1]=11;
  }
  int spbinhigh[2];
  if(qcut==2){
    spbinhigh[0]=10;
    spbinhigh[1]=28;
  }else if(qcut==1){
    spbinhigh[0]=10;
    spbinhigh[1]=21;
  }
  
  double nSporK0[2]={0.0,0.0};
  double nK0[2][3];
  
  nSporK0[0]=wbin1_Sp->Integral(spbinlow[0],spbinhigh[0]);
  nSporK0[1]=wbin2_Sp->Integral(spbinlow[1],spbinhigh[1]);
  nK0[0][0]=IMnpimK0inter_Sp_sysdown->Integral(spbinlow[0],spbinhigh[0]);
  nK0[1][0]=IMnpimK0inter_Sp_sysdown->Integral(spbinlow[1],spbinhigh[1]);
  nK0[0][1]=IMnpimK0inter_Sp->Integral(spbinlow[0],spbinhigh[0]);
  nK0[1][1]=IMnpimK0inter_Sp->Integral(spbinlow[1],spbinhigh[1]);
  nK0[0][2]=IMnpimK0inter_Sp_sysup->Integral(spbinlow[0],spbinhigh[0]);
  nK0[1][2]=IMnpimK0inter_Sp_sysup->Integral(spbinlow[1],spbinhigh[1]);


  std::cout << "total_wbin1: " <<  nSporK0[0] << std::endl;
  std::cout << "nK0_wbin1: " <<  nK0[0][1] << std::endl;
  std::cout << "nK0_wbin1 sysdown: " <<  nK0[0][0] << std::endl;
  std::cout << "nK0_wbin1 sysup: " <<  nK0[0][2] << std::endl;
  std::cout << "total_wbin2: " <<  nSporK0[1] << std::endl;
  std::cout << "nK0_wbin2: " <<  nK0[1][1] << std::endl;
  std::cout << "nK0_wbin2 sysdown: " <<  nK0[1][0] << std::endl;
  std::cout << "nK0_wbin2 sysup: " <<  nK0[1][2] << std::endl;
  TGraphErrors *gr_Spratio_SporK0[3]; 
  for(int isys=0;isys<3;isys++){
    gr_Spratio_SporK0[isys] = new TGraphErrors();
    gr_Spratio_SporK0[isys]->SetName(Form("gr_Spratio_SporK0_sys%d",isys-1));
    gr_Spratio_SporK0[isys]->SetTitle(Form("gr_Spratio_SporK0_sys%d",isys-1));
    gr_Spratio_SporK0[isys]->SetPoint(0,1,(nSporK0[0]-nK0[0][isys])/nSporK0[0]);
    gr_Spratio_SporK0[isys]->SetPoint(1,2,(nSporK0[1]-nK0[1][isys])/nSporK0[1]);
  }
   
  TCanvas *cratioSp = new TCanvas("cratioSp","cratioSp");
  for(int isys=0;isys<3;isys++){
    gr_Spratio_SporK0[isys]->SetMarkerStyle(20);
  }
  gr_Spratio_SporK0[1]->GetYaxis()->SetRangeUser(0.1,0.6);
  gr_Spratio_SporK0[1]->Draw("ap");
  gr_Spratio_SporK0[0]->Draw("p");
  gr_Spratio_SporK0[2]->Draw("p");
  gr_Spratio_SporK0[1]->Print();
  gr_Spratio_SporK0[0]->Print();
  gr_Spratio_SporK0[2]->Print();
  


  TCanvas *ccomp_Smregion = new TCanvas("ccomp_Smregion");
  IMnpim_IMnpip_dE_wK0_wSid_n_Sm->ProjectionX()->Draw();
  const int Smbinlow = h2K0inter_3->GetYaxis()->FindBin(anacuts::Sigmam_MIN);
  const int Smbinhi = h2K0inter_3->GetYaxis()->FindBin(anacuts::Sigmam_MAX);
  TH1D* IMnpipK0inter_Sm = (TH1D*)h2K0inter_3->ProjectionX("IMnpip_K0inter_Sm",Smbinlow,Smbinhi);
  TH1D* IMnpipK0inter_Sm_sysdown = (TH1D*)h2K0inter_3_sysdown->ProjectionX("IMnpip_K0inter_Sm_sysdown",Smbinlow,Smbinhi);
  TH1D* IMnpipK0inter_Sm_sysup = (TH1D*)h2K0inter_3_sysup->ProjectionX("IMnpip_K0inter_Sm_sysup",Smbinlow,Smbinhi);
  IMnpipK0inter_Sm->SetLineColor(2);
  IMnpipK0inter_Sm->Draw("Hsame");
  TH1D* wbin1_Sm = (TH1D*)IMnpim_IMnpip_dE_wK0_wSid_n_Sm_wbin1->ProjectionX();
  wbin1_Sm->SetLineColor(3);
  wbin1_Sm->Draw("same");
  TH1D* wbin2_Sm = (TH1D*)IMnpim_IMnpip_dE_wK0_wSid_n_Sm_wbin2->ProjectionX();
  wbin2_Sm->SetLineColor(4);
  wbin2_Sm->Draw("same");

  int smbinlow[2];  
  if(qcut==2){
    smbinlow[0]=5;
    smbinlow[1]=12;
  }else if(qcut==1){
    smbinlow[0]=5;
    smbinlow[1]=10;
  }
  int smbinhigh[2];
  if(qcut==2){
    smbinhigh[0]=11;
    smbinhigh[1]=28;
  }else if(qcut==1){
    smbinhigh[0]=9;
    smbinhigh[1]=25;
  }
  
  
  double nSmorK0[2]={0.0,0.0};
  double nK0_SmorK0[2][3];
  
  nSmorK0[0]=wbin1_Sm->Integral(smbinlow[0],smbinhigh[0]);
  nSmorK0[1]=wbin2_Sm->Integral(smbinlow[1],smbinhigh[1]);
  nK0_SmorK0[0][0]=IMnpipK0inter_Sm_sysdown->Integral(smbinlow[0],smbinhigh[0]);
  nK0_SmorK0[1][0]=IMnpipK0inter_Sm_sysdown->Integral(smbinlow[1],smbinhigh[1]);
  nK0_SmorK0[0][1]=IMnpipK0inter_Sm->Integral(smbinlow[0],smbinhigh[0]);
  nK0_SmorK0[1][1]=IMnpipK0inter_Sm->Integral(smbinlow[1],smbinhigh[1]);
  nK0_SmorK0[0][2]=IMnpipK0inter_Sm_sysup->Integral(smbinlow[0],smbinhigh[0]);
  nK0_SmorK0[1][2]=IMnpipK0inter_Sm_sysup->Integral(smbinlow[1],smbinhigh[1]);
  
  double ratioSm_SmorK0[3][2];
  

  TGraphErrors *gr_Smratio_SmorK0[3]; 
  for(int isys=0;isys<3;isys++){
    gr_Smratio_SmorK0[isys] = new TGraphErrors();
    gr_Smratio_SmorK0[isys]->SetName(Form("gr_Smratio_SmorK0_sys%d",isys-1));
    gr_Smratio_SmorK0[isys]->SetTitle(Form("gr_Smratio_SmorK0_sys%d",isys-1));
    ratioSm_SmorK0[isys][0] = (nSmorK0[0]-nK0_SmorK0[0][isys])/nSmorK0[0];
    if(ratioSm_SmorK0[isys][0]<0) ratioSm_SmorK0[isys][0]=0.0;
    ratioSm_SmorK0[isys][1] = (nSmorK0[1]-nK0_SmorK0[1][isys])/nSmorK0[1];
    if(ratioSm_SmorK0[isys][1]<0) ratioSm_SmorK0[isys][1]=0.0;
    gr_Smratio_SmorK0[isys]->SetPoint(0,1,ratioSm_SmorK0[isys][0]);
    gr_Smratio_SmorK0[isys]->SetPoint(1,2,ratioSm_SmorK0[isys][1]);
  }
   
  TCanvas *cratioSm = new TCanvas("cratioSm","cratioSm");
  for(int isys=0;isys<3;isys++){
    gr_Smratio_SmorK0[isys]->SetMarkerStyle(20);
  }
  gr_Smratio_SmorK0[1]->GetYaxis()->SetRangeUser(0.0,0.6);
  gr_Smratio_SmorK0[1]->Draw("ap");
  gr_Smratio_SmorK0[0]->Draw("p");
  gr_Smratio_SmorK0[2]->Draw("p");
  std::cout << "total_wbin1: " <<  nSmorK0[0] << std::endl;
  std::cout << "nK0_wbin1: " <<  nK0_SmorK0[0][1] << std::endl;
  std::cout << "nK0_wbin1 sysdown: " <<  nK0_SmorK0[0][0] << std::endl;
  std::cout << "nK0_wbin1 sysup: " <<  nK0_SmorK0[0][2] << std::endl;
  std::cout << "total_wbin2: " <<  nSmorK0[1] << std::endl;
  std::cout << "nK0_wbin2: " <<  nK0_SmorK0[1][1] << std::endl;
  std::cout << "nK0_wbin2 sysdown: " <<  nK0_SmorK0[1][0] << std::endl;
  std::cout << "nK0_wbin2 sysup: " <<  nK0_SmorK0[1][2] << std::endl;
  gr_Smratio_SmorK0[1]->Print();
  gr_Smratio_SmorK0[0]->Print();
  gr_Smratio_SmorK0[2]->Print();



  //
  /////next step/////
  //
  //subtract K0 and solve Sp/Sm overlap region
  TCanvas *cwK0orwSid_n = new TCanvas("cwK0orwSid_n","cwK0orwSid_n",1600,800);
  cwK0orwSid_n->Divide(2,1);
  cwK0orwSid_n->cd(1);
  if(!NoRebin)IMnpim_IMnpip_dE_wK0orwSid_n->Rebin2D(RebinFactor,RebinFactor);
  //Sigma is bright in this histogram. K0 is included.
  IMnpim_IMnpip_dE_wK0orwSid_n->Draw("colz");
  //IMnpim_IMnpip_dE_wSid_n->Rebin2D(4,4);
  //IMnpim_IMnpip_dE_wSid_n->Draw("colz");
  
  cwK0orwSid_n->cd(2);
  TH2F* IMnpim_IMnpip_dE_wK0orwSid_n_K0sub = (TH2F*)IMnpim_IMnpip_dE_wK0orwSid_n->Clone("IMnpim_IMnpip_dE_wK0orwSid_n_K0sub");
  IMnpim_IMnpip_dE_wK0orwSid_n_K0sub->SetTitle("IMnpim_IMnpip_dE_wK0orwSid_n_K0sub");
  //subtract K0
  //IMnpim_IMnpip_dE_wK0orwSid_n_K0sub->Add(IMnpim_IMnpip_dE_wK0_woSid_n_2,-1.0);
  TH2F* IMnpim_IMnpip_dE_wK0_woSid_n_3_rebin = (TH2F*)IMnpim_IMnpip_dE_wK0_woSid_n_3->Clone("IMnpim_IMnpip_dE_wK0_woSid_n_3_rebin");
  IMnpim_IMnpip_dE_wK0_woSid_n_3_rebin->Rebin2D(RebinFactor,RebinFactor);
  IMnpim_IMnpip_dE_wK0orwSid_n_K0sub->Add(IMnpim_IMnpip_dE_wK0_woSid_n_3_rebin,-1.0);
  IMnpim_IMnpip_dE_wK0orwSid_n_K0sub->SetMaximum(IMnpim_IMnpip_dE_wK0orwSid_n->GetMaximum());
  //IMnpim_IMnpip_dE_wK0orwSid_n_K0sub->SetMinimum(0);
  IMnpim_IMnpip_dE_wK0orwSid_n_K0sub->Draw("colz");
  

  TCanvas *cwK0orwSid_n_pro = new TCanvas("cwK0orwSid_n_pro","cwK0orwSid_n_pro",1600,800);
  cwK0orwSid_n_pro->Divide(2,1);
  cwK0orwSid_n_pro->cd(1);
  TH1D* IMnpip_wK0orwSid_n = (TH1D*)IMnpim_IMnpip_dE_wK0orwSid_n->ProjectionX("IMnpip_wK0orwSid_n");
  IMnpip_wK0orwSid_n->Draw("E");
  IMnpip_wK0_woSid->SetLineColor(2);
  IMnpip_wK0_woSid->Draw("Esame");
  TLegend *leg = new TLegend(0.6,0.6,0.8,0.8);
  leg->AddEntry(IMnpip_wK0_woSid,"K0");
  leg->AddEntry(IMnpip_wK0orwSid_n,"Sigma+/-/K0 total");
  leg->Draw();

  cwK0orwSid_n_pro->cd(2);
  TH1D* IMnpim_wK0orwSid_n = (TH1D*)IMnpim_IMnpip_dE_wK0orwSid_n->ProjectionY("IMnpim_wK0orwSid_n");
  IMnpim_wK0orwSid_n->Draw("E");
  IMnpim_wK0_woSid->SetLineColor(2);
  IMnpim_wK0_woSid->Draw("Esame");


  //K0 subtracted events
  TCanvas *cwSid_n_K0sub = new TCanvas("cwSid_n_K0sub","cwK0orwSid_n",800,800);
  cwSid_n_K0sub->Divide(2,2);
  cwSid_n_K0sub->cd(3);
  IMnpim_IMnpip_dE_wK0orwSid_n_K0sub->Draw("colz");

  cwSid_n_K0sub->cd(1);
  //2 sigma, narrow cut
  const int SmbinStart = IMnpim_IMnpip_dE_wK0orwSid_n_K0sub->GetYaxis()->FindBin(anacuts::Sigmam_MIN);
//  const int SmbinStart = IMnpim_IMnpip_dE_wK0orwSid_n_K0sub->GetYaxis()->FindBin(anacuts::Sigmam_MIN_wide);
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
  std::cout << "overlap events "  << IMnpim_K0sub_Sm->GetBinContent(Smbin) << std::endl;
  std::cout << "overlap error "  << IMnpim_K0sub_Sm->GetBinError(Smbin) << std::endl;
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


  TCanvas *cwSid_n_K0sub_wo = new TCanvas("cwSid_n_K0sub_wo","cwK0orwSid_n_K0sub_wo",1600,800);
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
  

  TCanvas *cwSid_n_K0sub_wo_fit = new TCanvas("cwSid_n_K0sub_wo_fit","cwSid_n_K0sub_wo_fit",1600,800);
  cwSid_n_K0sub_wo_fit->Divide(2,1);
  cwSid_n_K0sub_wo_fit->cd(1);
  IMnpip_K0sub_woSp->Draw("HE");
  IMnpip_K0sub_woSp->Fit("pol1","","",anacuts::Sigmap_MIN_wide-anacuts::Sigmap_sigma,anacuts::Sigmap_MAX_wide+anacuts::Sigmap_sigma);
  TF1* pol1_npip = (TF1*)IMnpip_K0sub_woSp->GetFunction("pol1");
  pol1_npip->Print();

  cwSid_n_K0sub_wo_fit->cd(2);
  IMnpim_K0sub_woSm->Draw("HE");
  IMnpim_K0sub_woSm->Fit("pol1","","",anacuts::Sigmam_MIN_wide-anacuts::Sigmam_sigma,anacuts::Sigmam_MAX_wide+anacuts::Sigmam_sigma);
  TF1* pol1_npim = (TF1*)IMnpim_K0sub_woSm->GetFunction("pol1");
  double par_npim[2];
  //double *par_npimerr;
  pol1_npim->GetParameters(par_npim);
  //par_npimerr = pol1_npim->GetParErrors();

  TGraphErrors *gr_SmONnpip_fin_pol1 = new TGraphErrors(IMnpip_K0sub_woSp);
  TGraphErrors *gr_SpONnpim_fin_pol1 = new TGraphErrors(IMnpim_K0sub_woSm);
  gr_SmONnpip_fin_pol1->SetName("gr_SmONnpip_fin_pol1");
  gr_SpONnpim_fin_pol1->SetName("gr_SpONnpim_fin_pol1");
  gr_SmONnpip_fin_pol1->RemovePoint(6);
  gr_SpONnpim_fin_pol1->RemovePoint(7);
   
  //confirm removal bins avove are correct or not
  cwSid_n_K0sub_wo_fit->cd(1);
  gr_SmONnpip_fin_pol1->SetLineColor(3);
  gr_SmONnpip_fin_pol1->Draw("c");
  cwSid_n_K0sub_wo_fit->cd(2);
  gr_SpONnpim_fin_pol1->SetLineColor(3);
  gr_SpONnpim_fin_pol1->Draw("c");

  TGraphErrors *gr_SmONnpip_fin_pol1_inter = (TGraphErrors*)gr_SmONnpip_fin_pol1->Clone("gr_SmONnpip_fin_pol1_inter");
  gr_SmONnpip_fin_pol1_inter->SetName("gr_SmONnpip_fin_pol1_inter");
  TGraphErrors *gr_SpONnpim_fin_pol1_inter = (TGraphErrors*)gr_SpONnpim_fin_pol1->Clone("gr_SpONnpim_fin_pol1_inter");
  gr_SpONnpim_fin_pol1_inter->SetName("gr_SpONnpim_fin_pol1_inter");

  TCanvas *cwSid_n_K0sub_wo_fit2 = new TCanvas("cwSid_n_K0sub_wo_fit2","cwSid_n_K0sub_wo_fit2",1600,800);
  cwSid_n_K0sub_wo_fit2->Divide(2,1);
  cwSid_n_K0sub_wo_fit2->cd(1);
  IMnpip_K0sub->Draw("HE");
  //gr_SmONnpip_fin_pol1->Draw("c");
  double SmEstimate_onnpipcross = pol1_npip->Eval(anacuts::Sigmap_center);
  double SmEstimate_err_low = IMnpim_K0sub_woSm->GetBinError(IMnpim_K0sub_woSm->FindBin(anacuts::Sigmam_MIN_wide-anacuts::Sigmam_sigma));
  double SmEstimate_err_high = IMnpim_K0sub_woSm->GetBinError(IMnpim_K0sub_woSm->FindBin(anacuts::Sigmam_MIN_wide+anacuts::Sigmam_sigma));
  double SmEstimate_err = sqrt(SmEstimate_err_low*SmEstimate_err_low+SmEstimate_err_high*SmEstimate_err_high)/2.0;

  gr_SmONnpip_fin_pol1_inter->SetPoint(gr_SmONnpip_fin_pol1->GetN(),anacuts::Sigmap_center,SmEstimate_onnpipcross);
  gr_SmONnpip_fin_pol1_inter->SetPointError(gr_SmONnpip_fin_pol1->GetN(),0,SmEstimate_err);
  gr_SmONnpip_fin_pol1_inter->SetMarkerStyle(20);
  gr_SmONnpip_fin_pol1_inter->SetMarkerColor(3);
  gr_SmONnpip_fin_pol1_inter->Draw("p");
  
  cwSid_n_K0sub_wo_fit2->cd(2);
  IMnpim_K0sub->Draw("HE");
  gr_SpONnpim_fin_pol1->Draw("p");
  double SpEstimate_onnpimcross = pol1_npim->Eval(anacuts::Sigmam_center);
  double SpEstimate_err_low = IMnpip_K0sub_woSp->GetBinError(IMnpip_K0sub_woSp->FindBin(anacuts::Sigmap_MIN_wide-anacuts::Sigmap_sigma));
  double SpEstimate_err_high = IMnpip_K0sub_woSp->GetBinError(IMnpip_K0sub_woSp->FindBin(anacuts::Sigmap_MIN_wide+anacuts::Sigmap_sigma));
  double SpEstimate_err = sqrt(SpEstimate_err_low*SpEstimate_err_low+SpEstimate_err_high*SpEstimate_err_high)/2.0;
  gr_SpONnpim_fin_pol1_inter->SetPoint(gr_SpONnpim_fin_pol1->GetN(),anacuts::Sigmam_center,SpEstimate_onnpimcross);
  gr_SpONnpim_fin_pol1_inter->SetPointError(gr_SpONnpim_fin_pol1->GetN(),0,SpEstimate_err);
  gr_SpONnpim_fin_pol1_inter->SetMarkerStyle(20);
  gr_SpONnpim_fin_pol1_inter->SetMarkerColor(3);
  gr_SpONnpim_fin_pol1_inter->Draw("p");

  TCanvas *cwSid_n_K0sub_wo_fit3 = new TCanvas("cwSid_n_K0sub_wo_fit3","cwSid_n_K0sub_wo_fit3",1600,800);
  cwSid_n_K0sub_wo_fit3->Divide(2,1);
  cwSid_n_K0sub_wo_fit3->cd(1);
  IMnpip_K0sub->Draw("HE");
  cwSid_n_K0sub_wo_fit3->cd(2);
  IMnpim_K0sub->Draw("HE");
  double crossCount = IMnpip_K0sub->GetBinContent(IMnpip_K0sub->FindBin(anacuts::Sigmap_center));
  double SpEstimate_final = (SpEstimate_onnpimcross*SpEstimate_err+(crossCount-SmEstimate_onnpipcross)*SmEstimate_err)/(SpEstimate_err+SmEstimate_err);
  // 
  //double SpEstimate_devi = (SpEstimate_onnpimcross*SpEstimate_err-(crossCount-SmEstimate_onnpipcross)*SmEstimate_err)/(SpEstimate_err+SmEstimate_err);
  double SpEstimate_devi = fabs(SpEstimate_onnpimcross-(crossCount-SmEstimate_onnpipcross));
  //double SpEstimate_err_final = sqrt(SpEstimate_err*SpEstimate_err+SmEstimate_err*SmEstimate_err);
  //this is overesitmate, but To be conservative.
  double SpEstimate_err_final = sqrt(SpEstimate_err*SpEstimate_err+SmEstimate_err*SmEstimate_err+SpEstimate_devi*SpEstimate_devi);

  std::cout << "Sp on pim " << SpEstimate_onnpimcross << std::endl;
  std::cout << "cross - Sm " << crossCount - SmEstimate_onnpipcross << std::endl;
  std::cout << "Sp on pim final " << SpEstimate_final << std::endl;
  std::cout << "SpEstimate_final devi." <<  SpEstimate_devi  <<   std::endl;
  std::cout << "SpEstimate_final err" <<  SpEstimate_err_final  <<   std::endl;

  TGraphErrors *gr_SmONnpip_fin_pol1_final = (TGraphErrors*)gr_SmONnpip_fin_pol1->Clone("gr_SmONnpip_fin_pol1_final");
  gr_SmONnpip_fin_pol1_final->SetName("gr_SmONnpip_fin_pol1_final");
  gr_SmONnpip_fin_pol1_final->SetPoint(gr_SmONnpip_fin_pol1->GetN(),anacuts::Sigmap_center,crossCount-SpEstimate_final);
  gr_SmONnpip_fin_pol1_final->SetPointError(gr_SmONnpip_fin_pol1->GetN(),0,SpEstimate_err_final);
  
  TGraphErrors *gr_SpONnpim_fin_pol1_final = (TGraphErrors*)gr_SpONnpim_fin_pol1->Clone("gr_SpONnpim_fin_pol1_final");
  gr_SpONnpim_fin_pol1_final->SetName("gr_SpONnpim_fin_pol1_final");
  gr_SpONnpim_fin_pol1_final->SetPoint(gr_SpONnpim_fin_pol1->GetN(),anacuts::Sigmam_center,SpEstimate_final);
  gr_SpONnpim_fin_pol1_final->SetPointError(gr_SpONnpim_fin_pol1->GetN(),0,SpEstimate_err_final);
  cwSid_n_K0sub_wo_fit3->cd(1);
  gr_SmONnpip_fin_pol1_final->Draw("p");
  cwSid_n_K0sub_wo_fit3->cd(2);
  gr_SpONnpim_fin_pol1_final->Draw("p");
  TH2F* q_IMnpipi_wK0_wSid_n_SpSm = (TH2F*)fr->Get("q_IMnpipi_wK0_wSid_n_SpSm");
  TH2F* IMnpim_IMnpip_wK0_wSid_n_SpSm = (TH2F*)fr->Get("IMnpim_IMnpip_dE_wK0_wSid_n_SpSm");
  
  TCanvas *c3overlap = new TCanvas("c3overlap","c3overlap",800,800);
  q_IMnpipi_wK0_wSid_n_SpSm->ProjectionX("IMnpipi_wK0_wSid_n")->Draw();
  
  TCanvas *c3overlap2 = new TCanvas("c3overlap2","c3overlap2",800,800);
  IMnpim_IMnpip_wK0_wSid_n_SpSm->SetMinimum(0);
  IMnpim_IMnpip_wK0_wSid_n_SpSm->Draw("colz");
   
  cwK0orwSid_n_pro->cd(1);
  gr_SmONnpip_fin_pol1_final->Draw("p");
  cwK0orwSid_n_pro->cd(2);
  gr_SpONnpim_fin_pol1_final->Draw("p");
  
  //final decomposition summary plots for the paper 
  TCanvas *cFinalDeco = new TCanvas("cFinalDeco","cFinalDeco",1600,800);
  cFinalDeco->Divide(2,1);
  cFinalDeco->cd(1);
  IMnpip_wK0orwSid_n->Draw("E");
  //IMnpip_wK0_woSid->SetLineColor(2);
  //IMnpip_wK0_woSid->Draw("Esame");
  //IMnpip_3_inter->Draw("same");
  //IMnpip_3_inter_sysup->Draw("same");
  //IMnpip_3_inter_sysdown->Draw("same");
  TGraphAsymmErrors *gr_K0inter_onnpip = new TGraphAsymmErrors(IMnpip_3_inter);
  gr_K0inter_onnpip->SetName("gr_K0inter_onnpip");
  TGraphAsymmErrors *gr_K0inter_onnpip_sysup = new TGraphAsymmErrors(IMnpip_3_inter_sysup);
  gr_K0inter_onnpip_sysup->SetName("gr_K0inter_onnpip_sysup");
  TGraphAsymmErrors *gr_K0inter_onnpip_sysdown = new TGraphAsymmErrors(IMnpip_3_inter_sysdown);
  gr_K0inter_onnpip_sysdown->SetName("gr_K0inter_onnpip_sysdown");
  //TGraphAsymmErrors *gr_K0inter_onnpip = new TGraphAsymmErrors(f3widehist_px);
  //TGraphAsymmErrors *gr_K0inter_onnpip_sysup = new TGraphAsymmErrors(f3widehist_sysup_px);
  //TGraphAsymmErrors *gr_K0inter_onnpip_sysdown = new TGraphAsymmErrors(f3widehist_sysdown_px);
  for(int ip=0;ip<gr_K0inter_onnpip->GetN();ip++){
    double cont = gr_K0inter_onnpip->GetPointY(ip);
    double contup = gr_K0inter_onnpip_sysup->GetPointY(ip);
    double contdown = gr_K0inter_onnpip_sysdown->GetPointY(ip);
    double ystat = gr_K0inter_onnpip->GetErrorY(ip);
    double yehi = contup - cont; 
    double yelo = cont - contdown;
    gr_K0inter_onnpip->SetPointEYhigh(ip,sqrt(yehi*yehi+ystat*ystat));
    gr_K0inter_onnpip->SetPointEYlow(ip,sqrt(yelo*yelo+ystat*ystat));
  }
  gr_K0inter_onnpip->SetMarkerColor(2);
  gr_K0inter_onnpip->SetLineColor(2);
  gr_K0inter_onnpip->Draw("p");
 
  TGraphErrors *gr_SmONnpip_fin_pol1_final_sort = (TGraphErrors*)gr_SmONnpip_fin_pol1_final->Clone("gr_SmONnpip_fin_pol1_final_sort");
  gr_SmONnpip_fin_pol1_final_sort->SetName("gr_SmONnpip_fin_pol1_final_sort");
  gr_SmONnpip_fin_pol1_final_sort->Sort();
  gr_SmONnpip_fin_pol1_final_sort->SetFillStyle(0);
  gr_SmONnpip_fin_pol1_final_sort->SetFillColor(3);
  gr_SmONnpip_fin_pol1_final_sort->Draw("p");
  //auto* h_SmONnpip_fin_pol1_final = gr_SmONnpip_fin_pol1_final->GetHistogram();
  //h_SmONnpip_fin_pol1_final->SetLineColor(3);
  //h_SmONnpip_fin_pol1_final->Draw("Esame");
  //h_SmONnpip_fin_pol1_final->Print("");
  TGraphErrors *gr_K0inter_plusSm_onnpip = (TGraphErrors*)gr_K0inter_onnpip->Clone("gr_K0inter_plusSm_onnpip");
  gr_K0inter_plusSm_onnpip->SetName("gr_K0inter_plusSm_onnpip");
  
  //->current procudure does not consider tail components of signal outside of +/- 2 sigma
  //but such an arbitrary factor is not good...
  //so I fix it at 1
  //original strategy is to determine the ratio of K0/S+/S- on overlap region, fixing the counts.
  //so that summary plots which conserve the counts will be made later, converting xaxis to MeV.
  const double TailFactor = 1.0;// 
  for(int ip=0;ip<gr_K0inter_plusSm_onnpip->GetN();ip++){
    double x = gr_K0inter_plusSm_onnpip->GetPointX(ip);
    double y = gr_K0inter_plusSm_onnpip->GetPointY(ip);
    double addval = gr_SmONnpip_fin_pol1_final_sort->Eval(x);
    gr_K0inter_plusSm_onnpip->SetPointY(ip,y+addval*TailFactor);
  }
  gr_K0inter_plusSm_onnpip->SetLineColor(4);
  gr_K0inter_plusSm_onnpip->SetMarkerColor(4);
  gr_K0inter_plusSm_onnpip->Draw("P");

  TLegend *legf = new TLegend(0.6,0.6,0.8,0.8);
  //legf->AddEntry(IMnpip_wK0_woSid,"K0");
  legf->AddEntry(IMnpip_3_inter,"K0");
  legf->AddEntry(IMnpip_wK0orwSid_n,"Sigma+/-/K0 total");
  legf->Draw();
  
  cFinalDeco->cd(2);
  IMnpim_wK0orwSid_n->Draw("E");
  //IMnpim_wK0_woSid->SetLineColor(2);
  //IMnpim_wK0_woSid->Draw("Esame");
  //IMnpim_3_inter->Draw("same");
  //IMnpim_3_inter_sysup->Draw("same");
  //IMnpim_3_inter_sysdown->Draw("same");
  TGraphAsymmErrors *gr_K0inter_onnpim = new TGraphAsymmErrors(IMnpim_3_inter);
  gr_K0inter_onnpim->SetName("gr_K0inter_onnpim");
  TGraphAsymmErrors *gr_K0inter_onnpim_sysup = new TGraphAsymmErrors(IMnpim_3_inter_sysup);
  gr_K0inter_onnpim_sysup->SetName("gr_K0inter_onnpim_sysup");
  TGraphAsymmErrors *gr_K0inter_onnpim_sysdown = new TGraphAsymmErrors(IMnpim_3_inter_sysdown);
  gr_K0inter_onnpim_sysdown->SetName("gr_K0inter_onnpim_sysdown");
  for(int ip=0;ip<gr_K0inter_onnpim->GetN();ip++){
    double cont = gr_K0inter_onnpim->GetPointY(ip);
    double contup = gr_K0inter_onnpim_sysup->GetPointY(ip);
    double contdown = gr_K0inter_onnpim_sysdown->GetPointY(ip);
    double ystat = gr_K0inter_onnpim->GetErrorY(ip);
    double yehi = contup - cont; 
    double yelo = cont - contdown;
    gr_K0inter_onnpim->SetPointEYhigh(ip,sqrt(yehi*yehi+ystat*ystat));
    gr_K0inter_onnpim->SetPointEYlow(ip,sqrt(yelo*yelo+ystat*ystat));
  }
  gr_K0inter_onnpim->SetMarkerColor(2);
  gr_K0inter_onnpim->SetLineColor(2);
  gr_K0inter_onnpim->Draw("p");
  //gr_SpONnpim_fin_pol1_final->Draw("p");
  TGraphErrors *gr_SpONnpim_fin_pol1_final_sort = (TGraphErrors*)gr_SpONnpim_fin_pol1_final->Clone("gr_SpONnpim_fin_pol1_final_sort");
  gr_SpONnpim_fin_pol1_final_sort->SetName("gr_SpONnpim_fin_pol1_final_sort");
  gr_SpONnpim_fin_pol1_final_sort->Sort();
  gr_SpONnpim_fin_pol1_final_sort->SetFillStyle(0);
  gr_SpONnpim_fin_pol1_final_sort->SetFillColor(3);
  gr_SpONnpim_fin_pol1_final_sort->Draw("p");
  TGraphErrors *gr_K0inter_plusSp_onnpim = (TGraphErrors*)gr_K0inter_onnpim->Clone("gr_K0inter_plusSp_onnpim");
  gr_K0inter_plusSp_onnpim->SetName("gr_K0inter_plusSp_onnpim");
  for(int ip=0;ip<gr_K0inter_plusSp_onnpim->GetN();ip++){
    double x = gr_K0inter_plusSp_onnpim->GetPointX(ip);
    double y = gr_K0inter_plusSp_onnpim->GetPointY(ip);
    double addval = gr_SpONnpim_fin_pol1_final_sort->Eval(x);
    gr_K0inter_plusSp_onnpim->SetPointY(ip,y+addval*TailFactor);
  }
  gr_K0inter_plusSp_onnpim->SetLineColor(4);
  gr_K0inter_plusSp_onnpim->SetMarkerColor(4);
  gr_K0inter_plusSp_onnpim->Draw("P");


  //summary plot
  //gev->mev conversion for ptep paper
  const double HISTMAX = 1800;
  const int fillstyle = 1001;
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  TCanvas *cFinalDeco_mev = new TCanvas("cFinalDeco_mev","cFinalDeco_mev",1800,800);
  cFinalDeco_mev->Divide(2,1);
  cFinalDeco_mev->cd(1);
  const int npipbin = IMnpip_wK0orwSid_n->GetNbinsX();
  const double npipmin = (IMnpip_wK0orwSid_n->GetXaxis()->GetXmin())*1000;
  const double npipmax = (IMnpip_wK0orwSid_n->GetXaxis()->GetXmax())*1000;
  TH1D* IMnpip_wK0orwSid_n_mev = new TH1D("IMnpip_wK0orwSid_n_mev","",npipbin,npipmin,npipmax);
  TH1D* IMnpip_wK0orwSid_n_mev_disp = new TH1D("IMnpip_wK0orwSid_n_mev_disp","",npipbin,npipmin,npipmax);
  for(int ibin = 0;ibin<IMnpip_wK0orwSid_n->GetNbinsX();ibin++){
    double x =  IMnpip_wK0orwSid_n->GetBinCenter(ibin+1)*1000.;
    double y =  IMnpip_wK0orwSid_n->GetBinContent(ibin+1);
    double ye = IMnpip_wK0orwSid_n->GetBinError(ibin+1);
    int Spbin = IMnpip_wK0orwSid_n->GetXaxis()->FindBin(anacuts::Sigmap_center);
    IMnpip_wK0orwSid_n_mev->SetBinContent(ibin+1,y);
    IMnpip_wK0orwSid_n_mev->SetBinError(ibin+1,ye);
    if(Spbin == ibin+1){
      IMnpip_wK0orwSid_n_mev_disp->SetBinContent(ibin+1,y);
      IMnpip_wK0orwSid_n_mev_disp->SetBinError(ibin+1,ye);
    }
  }
  gPad->SetRightMargin(0.08);
  gPad->SetLeftMargin(0.14);  

  IMnpip_wK0orwSid_n_mev_disp->GetXaxis()->SetTitle("IM(n#pi^{+}) [MeV/c^{2}]");
  IMnpip_wK0orwSid_n_mev_disp->GetXaxis()->CenterTitle();
  IMnpip_wK0orwSid_n_mev_disp->GetXaxis()->SetRangeUser(1050,1500);
  IMnpip_wK0orwSid_n_mev_disp->GetYaxis()->SetTitle("counts");
  IMnpip_wK0orwSid_n_mev_disp->GetYaxis()->CenterTitle();
  IMnpip_wK0orwSid_n_mev_disp->GetYaxis()->SetTitleOffset(1.6);
//  IMnpip_wK0orwSid_n_mev_disp->SetMaximum(IMnpip_wK0orwSid_n_mev->GetMaximum()*1.11);
  IMnpip_wK0orwSid_n_mev_disp->SetMaximum(HISTMAX);
  IMnpip_wK0orwSid_n_mev_disp->SetMinimum(-100);
  IMnpip_wK0orwSid_n_mev_disp->SetLineWidth(2);
  IMnpip_wK0orwSid_n_mev_disp->SetLineColor(0);
  IMnpip_wK0orwSid_n_mev_disp->SetFillColor(46);
  IMnpip_wK0orwSid_n_mev_disp->SetFillStyle(fillstyle);
  IMnpip_wK0orwSid_n_mev_disp->Draw("HE");
  TGraphErrors *gr_K0inter_onnpip_mev = new TGraphErrors();
  for(int ip=0;ip<gr_K0inter_onnpip->GetN();ip++){
    double x = gr_K0inter_onnpip->GetPointX(ip);
    double y = gr_K0inter_onnpip->GetPointY(ip);
    double xe = gr_K0inter_onnpip->GetErrorX(2);
    double ye = gr_K0inter_onnpip->GetErrorY(ip);
    gr_K0inter_onnpip_mev->SetPoint(ip,x*1000.,y);
    gr_K0inter_onnpip_mev->SetPointError(ip,xe*1000.,ye);
  }
  gr_K0inter_onnpip_mev->SetMarkerColor(3);
  gr_K0inter_onnpip_mev->SetLineColor(3);
  gr_K0inter_onnpip_mev->SetLineWidth(0);
  TGraphErrors *gr_SmOnnpip_fin_mev = new TGraphErrors();
  for(int ip=0;ip<gr_SmONnpip_fin_pol1_final_sort->GetN();ip++){
    double x = gr_SmONnpip_fin_pol1_final_sort->GetPointX(ip);
    double y = gr_SmONnpip_fin_pol1_final_sort->GetPointY(ip);
    double xe = gr_SmONnpip_fin_pol1_final_sort->GetErrorX(2);
    double ye = gr_SmONnpip_fin_pol1_final_sort->GetErrorY(ip);
    gr_SmOnnpip_fin_mev->SetPoint(ip,x*1000.,y);
    gr_SmOnnpip_fin_mev->SetPointError(ip,xe*1000.,ye);
  }
  gr_SmOnnpip_fin_mev->SetMarkerColor(3);
  gr_SmOnnpip_fin_mev->SetLineColor(3);
  gr_SmOnnpip_fin_mev->SetLineWidth(2);
  gr_SmOnnpip_fin_mev->RemovePoint(0);
  gr_SmOnnpip_fin_mev->RemovePoint(0);
  TH1D* h_SmOnnpip_mev = (TH1D*)IMnpip_wK0orwSid_n_mev->Clone("h_SmOnnpip_mev");
  h_SmOnnpip_mev->Reset();
  for(int ip=0;ip<gr_SmOnnpip_fin_mev->GetN();ip++){
    double x,y;
    gr_SmOnnpip_fin_mev->GetPoint(ip,x,y);
    double ye = gr_SmOnnpip_fin_mev->GetErrorY(ip);
    double hbin = h_SmOnnpip_mev->GetXaxis()->FindBin(x);
    h_SmOnnpip_mev->SetBinContent(hbin,y);
    h_SmOnnpip_mev->SetBinError(hbin,ye);
  }
  h_SmOnnpip_mev->SetFillStyle(fillstyle);
  h_SmOnnpip_mev->SetMarkerColor(3);
  h_SmOnnpip_mev->SetLineColor(3);
  h_SmOnnpip_mev->SetFillColor(38);
  h_SmOnnpip_mev->SetLineWidth(2);
  
  TH1D* h_K0inter_onnpip_mev = (TH1D*)IMnpip_wK0orwSid_n_mev->Clone("h_K0inter_onnpip_mev");
  h_K0inter_onnpip_mev->Reset();
  for(int ip=0;ip<gr_K0inter_onnpip_mev->GetN();ip++){
    double x,y;
    gr_K0inter_onnpip_mev->GetPoint(ip,x,y);
    double ye = gr_K0inter_onnpip_mev->GetErrorY(ip);
    double hbin = h_K0inter_onnpip_mev->GetXaxis()->FindBin(x);
    h_K0inter_onnpip_mev->SetBinContent(hbin,y);
    h_K0inter_onnpip_mev->SetBinError(hbin,ye);
  }
  h_K0inter_onnpip_mev->SetFillStyle(fillstyle);
  h_K0inter_onnpip_mev->SetMarkerColor(3);
  h_K0inter_onnpip_mev->SetLineColor(3);
  h_K0inter_onnpip_mev->SetFillColor(30);
  h_K0inter_onnpip_mev->SetLineWidth(2);


  TLine *p = new TLine(1050,0,1500,0);
  p->SetLineColor(1);
  //p->SetLineWidth(2.0);
  p->SetLineStyle(2);
  p->Draw();
  TGraphErrors *gr_K0inter_plusSm_onnpip_mev = new TGraphErrors();
  for(int ip=0;ip<gr_K0inter_plusSm_onnpip->GetN();ip++){
    double x = gr_K0inter_plusSm_onnpip->GetPointX(ip);
    double y = gr_K0inter_plusSm_onnpip->GetPointY(ip);
    double xe = gr_K0inter_plusSm_onnpip->GetErrorX(ip);
    double ye = gr_K0inter_plusSm_onnpip->GetErrorY(ip);
    int ibin = IMnpip_wK0orwSid_n_mev->GetXaxis()->FindBin(x*1000.);
    double val = IMnpip_wK0orwSid_n_mev->GetBinContent(ibin);
    double vale = IMnpip_wK0orwSid_n_mev->GetBinError(ibin);
    int Spbin = IMnpip_wK0orwSid_n_mev->GetXaxis()->FindBin(anacuts::Sigmap_center*1000.);
    if(Spbin == ibin){
       val = y;
       vale = ye;
    }
    gr_K0inter_plusSm_onnpip_mev->SetPoint(ip,x*1000.,val);
    gr_K0inter_plusSm_onnpip_mev->SetPointError(ip,xe*1000.,vale);
  }
  gr_K0inter_plusSm_onnpip_mev->SetMarkerColor(4);
  gr_K0inter_plusSm_onnpip_mev->SetLineColor(4);
  gr_K0inter_plusSm_onnpip_mev->SetFillColor(38);
  gr_K0inter_plusSm_onnpip_mev->SetFillStyle(fillstyle);
  //gr_K0inter_plusSm_onnpip_mev->SetLineWidth(2);
  gr_K0inter_plusSm_onnpip_mev->SetLineWidth(0);
  gr_K0inter_plusSm_onnpip_mev->Draw("Fc");
  gr_SmOnnpip_fin_mev->SetFillStyle(fillstyle);
  gr_SmOnnpip_fin_mev->SetFillColor(30);
  gr_SmOnnpip_fin_mev->SetLineColor(30);
  gr_K0inter_onnpip_mev->SetFillStyle(fillstyle);
  gr_K0inter_onnpip_mev->SetFillColor(30);
  gr_K0inter_onnpip_mev->SetLineColor(30);
  //gr_SmOnnpip_fin_mev->Draw("FB1");
  //h_SmOnnpip_mev->Draw("E1Hsame");
  gr_K0inter_onnpip_mev->Draw("Fc");
  //h_K0inter_onnpip_mev->Draw("HE0same");
  
  TGraphErrors *g1 = new TGraphErrors(IMnpip_wK0orwSid_n_mev);
  g1->SetLineColor(1);
  g1->SetMarkerColor(1);
  g1->Draw("p");  


  TLegend *legp = new TLegend(0.60,0.6,0.8,0.75);
  legp->SetLineWidth(0);   
  legp->SetFillStyle(0);   
  //legp->AddEntry(g1,"data");
  legp->AddEntry(IMnpip_wK0orwSid_n_mev_disp,"#pi^{-}#Sigma^{+}n","f");
  legp->AddEntry(gr_K0inter_plusSm_onnpip_mev,"#pi^{+}#Sigma^{-}n","f");
  legp->AddEntry(gr_K0inter_onnpip_mev,"#bar{K}^{0}nn","f");
  legp->SetTextSize(0.04);
  legp->Draw();  

  TLatex *texnpip = new TLatex();
  double texnpip_ymax = HISTMAX*0.9;
  texnpip->SetTextSize(0.04);
  texnpip->SetTextColor(1);
  if(qcut == 1) texnpip->DrawLatex( 1100, texnpip_ymax, "(c)   0 < q < 300 MeV/c" );
  if(qcut == 2) texnpip->DrawLatex( 1100, texnpip_ymax, "(a)   300 < q < 650 MeV/c" );

  cFinalDeco_mev->cd(2);
  const int npimbin = IMnpim_wK0orwSid_n->GetNbinsX();
  const double npimmin = (IMnpim_wK0orwSid_n->GetXaxis()->GetXmin())*1000;
  const double npimmax = (IMnpim_wK0orwSid_n->GetXaxis()->GetXmax())*1000;
  TH1D* IMnpim_wK0orwSid_n_mev = new TH1D("IMnpim_wK0orwSid_n_mev","",npimbin,npimmin,npimmax);
  TH1D* IMnpim_wK0orwSid_n_mev_disp = new TH1D("IMnpim_wK0orwSid_n_mev_disp","",npimbin,npimmin,npimmax);
  for(int ibin = 0;ibin<IMnpim_wK0orwSid_n->GetNbinsX();ibin++){
    double x =  IMnpim_wK0orwSid_n->GetBinCenter(ibin+1)*1000.;
    double y =  IMnpim_wK0orwSid_n->GetBinContent(ibin+1);
    double ye = IMnpim_wK0orwSid_n->GetBinError(ibin+1);
    int Smbin = IMnpim_wK0orwSid_n->GetXaxis()->FindBin(anacuts::Sigmam_center);
    IMnpim_wK0orwSid_n_mev->SetBinContent(ibin+1,y);
    IMnpim_wK0orwSid_n_mev->SetBinError(ibin+1,ye);
    if(Smbin == ibin+1){
      IMnpim_wK0orwSid_n_mev_disp->SetBinContent(ibin+1,y);
      IMnpim_wK0orwSid_n_mev_disp->SetBinError(ibin+1,ye);
    }
  }
  gPad->SetRightMargin(0.08);
  gPad->SetLeftMargin(0.14);  

  IMnpim_wK0orwSid_n_mev_disp->GetXaxis()->SetTitle("IM(n#pi^{-}) [MeV/c^{2}]");
  IMnpim_wK0orwSid_n_mev_disp->GetXaxis()->CenterTitle();
  IMnpim_wK0orwSid_n_mev_disp->GetXaxis()->SetRangeUser(1050,1500);
  IMnpim_wK0orwSid_n_mev_disp->GetYaxis()->SetTitle("counts");
  IMnpim_wK0orwSid_n_mev_disp->GetYaxis()->CenterTitle();
  IMnpim_wK0orwSid_n_mev_disp->GetYaxis()->SetTitleOffset(1.6);
//  IMnpim_wK0orwSid_n_mev_disp->SetMaximum(IMnpim_wK0orwSid_n_mev->GetMaximum()*1.11);
  IMnpim_wK0orwSid_n_mev_disp->SetMinimum(-100);
  //IMnpim_wK0orwSid_n_mev_disp->SetMaximum(IMnpip_wK0orwSid_n_mev->GetMaximum()*1.11);
  IMnpim_wK0orwSid_n_mev_disp->SetMaximum(HISTMAX);
  IMnpim_wK0orwSid_n_mev_disp->SetLineWidth(2);
  IMnpim_wK0orwSid_n_mev_disp->SetLineColor(0);
  IMnpim_wK0orwSid_n_mev_disp->SetFillColor(38);
  IMnpim_wK0orwSid_n_mev_disp->SetFillStyle(fillstyle);
  IMnpim_wK0orwSid_n_mev_disp->Draw("HE");
  TGraphErrors *gr_K0inter_onnpim_mev = new TGraphErrors();
  for(int ip=0;ip<gr_K0inter_onnpim->GetN();ip++){
    double x = gr_K0inter_onnpim->GetPointX(ip);
    double y = gr_K0inter_onnpim->GetPointY(ip);
    double xe = gr_K0inter_onnpim->GetErrorX(2);
    double ye = gr_K0inter_onnpim->GetErrorY(ip);
    gr_K0inter_onnpim_mev->SetPoint(ip,x*1000.,y);
    gr_K0inter_onnpim_mev->SetPointError(ip,xe*1000.,ye);
  }
  gr_K0inter_onnpim_mev->SetMarkerColor(3);
  gr_K0inter_onnpim_mev->SetLineColor(3);
  gr_K0inter_onnpim_mev->SetLineWidth(0);
  TGraphErrors *gr_SpOnnpim_fin_mev = new TGraphErrors();
  for(int ip=0;ip<gr_SpONnpim_fin_pol1_final_sort->GetN();ip++){
    double x = gr_SpONnpim_fin_pol1_final_sort->GetPointX(ip);
    double y = gr_SpONnpim_fin_pol1_final_sort->GetPointY(ip);
    double xe = gr_SpONnpim_fin_pol1_final_sort->GetErrorX(2);
    double ye = gr_SpONnpim_fin_pol1_final_sort->GetErrorY(ip);
    gr_SpOnnpim_fin_mev->SetPoint(ip,x*1000.,y);
    gr_SpOnnpim_fin_mev->SetPointError(ip,xe*1000.,ye);
  }
  gr_SpOnnpim_fin_mev->SetMarkerColor(3);
  gr_SpOnnpim_fin_mev->SetLineColor(3);
  gr_SpOnnpim_fin_mev->SetLineWidth(2);
  gr_SpOnnpim_fin_mev->RemovePoint(0);
  gr_SpOnnpim_fin_mev->RemovePoint(0);
  TH1D* h_SpOnnpim_mev = (TH1D*)IMnpim_wK0orwSid_n_mev->Clone("h_SpOnnpim_mev");
  h_SpOnnpim_mev->Reset();
  for(int ip=0;ip<gr_SpOnnpim_fin_mev->GetN();ip++){
    double x,y;
    gr_SpOnnpim_fin_mev->GetPoint(ip,x,y);
    double ye = gr_SpOnnpim_fin_mev->GetErrorY(ip);
    double hbin = h_SpOnnpim_mev->GetXaxis()->FindBin(x);
    h_SpOnnpim_mev->SetBinContent(hbin,y);
    h_SpOnnpim_mev->SetBinError(hbin,ye);
  }
  h_SpOnnpim_mev->SetFillStyle(fillstyle);
  h_SpOnnpim_mev->SetMarkerColor(3);
  h_SpOnnpim_mev->SetLineColor(3);
  h_SpOnnpim_mev->SetFillColor(30);
  h_SpOnnpim_mev->SetLineWidth(2);
  
  TH1D* h_K0inter_onnpim_mev = (TH1D*)IMnpim_wK0orwSid_n_mev->Clone("h_K0inter_onnpim_mev");
  h_K0inter_onnpim_mev->Reset();
  for(int ip=0;ip<gr_K0inter_onnpim_mev->GetN();ip++){
    double x,y;
    gr_K0inter_onnpim_mev->GetPoint(ip,x,y);
    double ye = gr_K0inter_onnpim_mev->GetErrorY(ip);
    double hbin = h_K0inter_onnpim_mev->GetXaxis()->FindBin(x);
    h_K0inter_onnpim_mev->SetBinContent(hbin,y);
    h_K0inter_onnpim_mev->SetBinError(hbin,ye);
  }
  h_K0inter_onnpim_mev->SetFillStyle(fillstyle);
  h_K0inter_onnpim_mev->SetMarkerColor(3);
  h_K0inter_onnpim_mev->SetLineColor(30);
  h_K0inter_onnpim_mev->SetFillColor(30);
  h_K0inter_onnpim_mev->SetLineWidth(2);


  TGraphErrors *gr_K0inter_plusSp_onnpim_mev = new TGraphErrors();
  for(int ip=0;ip<gr_K0inter_plusSp_onnpim->GetN();ip++){
    double x = gr_K0inter_plusSp_onnpim->GetPointX(ip);
    double y = gr_K0inter_plusSp_onnpim->GetPointY(ip);
    double xe = gr_K0inter_plusSp_onnpim->GetErrorX(ip);
    double ye = gr_K0inter_plusSp_onnpim->GetErrorY(ip);
    int ibin = IMnpim_wK0orwSid_n_mev->GetXaxis()->FindBin(x*1000.);
    double val = IMnpim_wK0orwSid_n_mev->GetBinContent(ibin);
    double vale = IMnpim_wK0orwSid_n_mev->GetBinError(ibin);
    int Spbin = IMnpim_wK0orwSid_n_mev->GetXaxis()->FindBin(anacuts::Sigmap_center*1000.);
    if(Spbin == ibin){
       val = y;
       vale = ye;
    }
    gr_K0inter_plusSp_onnpim_mev->SetPoint(ip,x*1000.,val);
    gr_K0inter_plusSp_onnpim_mev->SetPointError(ip,xe*1000.,vale);
  }
  gr_K0inter_plusSp_onnpim_mev->SetMarkerColor(2);
  gr_K0inter_plusSp_onnpim_mev->SetLineColor(2);
  gr_K0inter_plusSp_onnpim_mev->SetFillColor(46);
  gr_K0inter_plusSp_onnpim_mev->SetFillStyle(fillstyle);
  //gr_K0inter_plusSp_onnpim_mev->SetLineWidth(2);
  gr_K0inter_plusSp_onnpim_mev->SetLineWidth(0);
  gr_K0inter_plusSp_onnpim_mev->Draw("Fc");
  gr_SpOnnpim_fin_mev->SetFillStyle(fillstyle);
  gr_SpOnnpim_fin_mev->SetFillColor(46);
  gr_SpOnnpim_fin_mev->SetLineColor(46);
  gr_K0inter_onnpim_mev->SetFillStyle(fillstyle);
  gr_K0inter_onnpim_mev->SetFillColor(30);
  gr_K0inter_onnpim_mev->SetLineColor(30);
  //gr_SpOnnpim_fin_mev->Draw("FB1");
  //h_SpOnnpim_mev->Draw("E1Hsame");
  gr_K0inter_onnpim_mev->Draw("Fc");
  //h_K0inter_onnpim_mev->Draw("HE0same");
  
  TGraphErrors *g2 = new TGraphErrors(IMnpim_wK0orwSid_n_mev);
  g2->SetLineColor(1);
  g2->SetMarkerColor(1);
  g2->Draw("p");  
  p->Draw();

  TLatex *texnpim = new TLatex();
  double texnpim_ymax = HISTMAX*0.9;
  texnpim->SetTextSize(0.04);
  texnpim->SetTextColor(1);
  if(qcut == 1) texnpim->DrawLatex( 1100, texnpim_ymax, "(d)    0 < q < 300 MeV/c" );
  if(qcut == 2) texnpim->DrawLatex( 1250, texnpim_ymax, "(b)    300 < q < 650 MeV/c" );



  TFile *fout = NULL;
  if(qcut==1){
     fout = TFile::Open(Form("fout_qlo_v%d_dE%d_sys%d.root",Version,dEcut,sysud),"RECREATE");
  }else if(qcut==2){
     fout = TFile::Open(Form("fout_qhi_v%d_dE%d_sys%d.root",Version,dEcut,sysud),"RECREATE");
  }
  //IMnpip_K0sub->Write();
  //IMnpim_K0sub->Write();
  //IMnpip_K0sub_woSp->Write();
  //IMnpim_K0sub_woSm->Write();
  //IMnpim_IMnpip_dE_wK0orwSid_n_K0sub->Write();
  //IMnpim_IMnpip_dE_wK0_woSid_n_3_inter->Write();

  f3wide->Write();
  f3widehist->Write();
  gr_SmONnpip_fin_pol1_final->Write();
  gr_SpONnpim_fin_pol1_final->Write();
  h2K0inter_3fine->Write(); 
  h2K0inter_3fine_sysup->Write(); 
  h2K0inter_3fine_sysdown->Write(); 
  for(int isys=0;isys<3;isys++){
    gr_Spratio_SporK0[isys]->Write();
    gr_Smratio_SmorK0[isys]->Write();
  }

  TCanvas *c = NULL;
  TSeqCollection *SCol = gROOT->GetListOfCanvases();
  int size = SCol->GetSize();
  TIter next(SCol);
  TString pdfname;
  if(qcut==1){
    pdfname= Form("Fit2DK0_qlo_datav%d_dE%d_sys%d.pdf",Version,dEcut,sysud);
  }else if(qcut==2){
    pdfname= Form("Fit2DK0_qhi_datav%d_dE%d_sys%d.pdf",Version,dEcut,sysud);
  }


  for(int i=0;i<size;i++){
    //pdf->NewPage();
    c= (TCanvas*)next();
    c->Draw();
    c->cd();
    //inside the canvas
    //TPaveText *pt = new TPaveText(.74,.81,0.9,0.90,"NDC");
    c->Modified();
    c->Update();
    std::cout << c->GetName() << std::endl;
    //make 1 pdf file
    if(i==0) c->Print(pdfname+"(",Form("pdf Title:%s",c->GetTitle()));
    else if(i==size-1)c->Print(pdfname+")",Form("pdf Title:%s",c->GetTitle())); 
    else c->Print(pdfname,Form("pdf Title:%s",c->GetTitle())); 
    
    //make separated pdf files
    //c->Print(Form("pdf/%s.pdf",c->GetTitle()));
  }

  if(qcut==1) cFinalDeco_mev->SaveAs("cDeco_qlo.pdf");
  else if(qcut==2) cFinalDeco_mev->SaveAs("cDeco_qhi.pdf");


}
