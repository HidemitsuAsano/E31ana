#include "anacuts.h"

const double CDHR= 55.9;//center cm

Double_t pimleftp(Double_t *x,Double_t *par)
{
  Double_t r1 = (anacuts::Isonpim_zcut)*sqrt(1.0 - pow(x[0]-anacuts::Isonpim_shift,2)/pow(anacuts::Isonpim_phicut_left,2));
 
  return r1;

}

Double_t pimleftm(Double_t *x,Double_t *par)
{
  Double_t r1 = -1.0*(anacuts::Isonpim_zcut)*sqrt(1.0 - pow(x[0]-anacuts::Isonpim_shift,2)/pow(anacuts::Isonpim_phicut_left,2));
 
  return r1;

}

Double_t pimrightp(Double_t *x,Double_t *par)
{
  Double_t r1 = (anacuts::Isonpim_zcut)*sqrt(1.0 - pow(x[0]-anacuts::Isonpim_shift,2)/pow(anacuts::Isonpim_phicut_right,2));
 
  return r1;

}

Double_t pimrightm(Double_t *x,Double_t *par)
{
  Double_t r1 = -1.0*(anacuts::Isonpim_zcut)*sqrt(1.0 - pow(x[0]-anacuts::Isonpim_shift,2)/pow(anacuts::Isonpim_phicut_right,2));
 
  return r1;
}

Double_t pipleftp(Double_t *x,Double_t *par)
{
  Double_t r1 = (anacuts::Isonpip_zcut)*sqrt(1.0 - pow(x[0]-anacuts::Isonpip_shift,2)/pow(anacuts::Isonpip_phicut_left,2));
 
  return r1;

}

Double_t pipleftm(Double_t *x,Double_t *par)
{
  Double_t r1 = -1.0*(anacuts::Isonpip_zcut)*sqrt(1.0 - pow(x[0]-anacuts::Isonpip_shift,2)/pow(anacuts::Isonpip_phicut_left,2));
 
  return r1;

}

Double_t piprightp(Double_t *x,Double_t *par)
{
  Double_t r1 = (anacuts::Isonpip_zcut)*sqrt(1.0 - pow(x[0]-anacuts::Isonpip_shift,2)/pow(anacuts::Isonpip_phicut_right,2));
 
  return r1;

}

Double_t piprightm(Double_t *x,Double_t *par)
{
  Double_t r1 = -1.0*(anacuts::Isonpip_zcut)*sqrt(1.0 - pow(x[0]-anacuts::Isonpip_shift,2)/pow(anacuts::Isonpip_phicut_right,2));
 
  return r1;
}

//convert phi radian to cm
Double_t pimleftp_len(Double_t *x,Double_t *par)
{
  Double_t r1 = (anacuts::Isonpim_zcut)*sqrt(1.0 - pow(x[0]-anacuts::Isonpim_shift*CDHR,2)/pow(anacuts::Isonpim_phicut_left*CDHR,2));
 
  return r1;

}

Double_t pimleftm_len(Double_t *x,Double_t *par)
{
  Double_t r1 = -1.0*(anacuts::Isonpim_zcut)*sqrt(1.0 - pow(x[0]-anacuts::Isonpim_shift*CDHR,2)/pow(anacuts::Isonpim_phicut_left*CDHR,2));
 
  return r1;

}

Double_t pimrightp_len(Double_t *x,Double_t *par)
{
  Double_t r1 = (anacuts::Isonpim_zcut)*sqrt(1.0 - pow(x[0]-anacuts::Isonpim_shift*CDHR,2)/pow(anacuts::Isonpim_phicut_right*CDHR,2));
 
  return r1;

}

Double_t pimrightm_len(Double_t *x,Double_t *par)
{
  Double_t r1 = -1.0*(anacuts::Isonpim_zcut)*sqrt(1.0 - pow(x[0]-anacuts::Isonpim_shift*CDHR,2)/pow(anacuts::Isonpim_phicut_right*CDHR,2));
 
  return r1;
}

Double_t pipleftp_len(Double_t *x,Double_t *par)
{
  Double_t r1 = (anacuts::Isonpip_zcut)*sqrt(1.0 - pow(x[0]-anacuts::Isonpip_shift*CDHR,2)/pow(anacuts::Isonpip_phicut_left*CDHR,2));
 
  return r1;

}

Double_t pipleftm_len(Double_t *x,Double_t *par)
{
  Double_t r1 = -1.0*(anacuts::Isonpip_zcut)*sqrt(1.0 - pow(x[0]-anacuts::Isonpip_shift*CDHR,2)/pow(anacuts::Isonpip_phicut_left*CDHR,2));
 
  return r1;

}

Double_t piprightp_len(Double_t *x,Double_t *par)
{
  Double_t r1 = (anacuts::Isonpip_zcut)*sqrt(1.0 - pow(x[0]-anacuts::Isonpip_shift*CDHR,2)/pow(anacuts::Isonpip_phicut_right*CDHR,2));
 
  return r1;

}

Double_t piprightm_len(Double_t *x,Double_t *par)
{
  Double_t r1 = -1.0*(anacuts::Isonpip_zcut)*sqrt(1.0 - pow(x[0]-anacuts::Isonpip_shift*CDHR,2)/pow(anacuts::Isonpip_phicut_right*CDHR,2));
 
  return r1;
}

void disp_IsolationCut()
{
  //before subtracting mixed events
  //TFile *_file0 = TFile::Open("evanaIMpisigma_npippim_v245_out_dE2_iso_nostop.root");
  gStyle->SetPalette(1);
  gStyle->SetFuncColor(kRed);
  TFile *_file0 = TFile::Open("evanaIMpisigma_v245.root");
  
  TH2D* diff2d_CDC_CDH_pim = (TH2D*)_file0->Get("diff2d_CDC_CDH_pim");
  TH2D* diff2d_CDC_CDH_pip = (TH2D*)_file0->Get("diff2d_CDC_CDH_pip");
  gStyle->SetOptStat(0);
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  TCanvas *c1 = new TCanvas("c1","c1",1600,600);
  c1->Divide(2,1,0.000001,0.0001);
  //gStyle->SetPadRightMargin(1.3);
  gStyle->SetPadLeftMargin(0.12);
  c1->cd(1);
  diff2d_CDC_CDH_pim->SetTitle("");
  diff2d_CDC_CDH_pim->GetYaxis()->SetNdivisions(5,5,0,kTRUE);
  diff2d_CDC_CDH_pip->GetYaxis()->SetNdivisions(5,5,0,kTRUE);
  diff2d_CDC_CDH_pim->GetZaxis()->SetNdivisions(5,5,0,kTRUE);
  diff2d_CDC_CDH_pip->GetZaxis()->SetNdivisions(5,5,0,kTRUE);
  diff2d_CDC_CDH_pim->SetXTitle("CDH hit - #pi^{-} track (phi) [radian]");
  diff2d_CDC_CDH_pim->GetXaxis()->CenterTitle();
  diff2d_CDC_CDH_pim->GetYaxis()->CenterTitle();
  diff2d_CDC_CDH_pip->GetXaxis()->CenterTitle();
  diff2d_CDC_CDH_pip->GetYaxis()->CenterTitle();
  diff2d_CDC_CDH_pip->SetXTitle("CDH hit - #pi^{+} track (phi) [radian]");
  diff2d_CDC_CDH_pim->SetYTitle("CDH hit - #pi^{-} track (z) [cm]");
  diff2d_CDC_CDH_pim->Draw("col");
   
  TF1 *f1 = new TF1("f1",pimleftp,-0.85425,0,1);
  TF1 *f2 = new TF1("f2",pimleftm,-0.85425,0,1);
  TF1 *f3 = new TF1("f3",pimrightp,0,0.85427,1);
  TF1 *f4 = new TF1("f4",pimrightm,0,0.85427,1);
  f1->SetLineWidth(4);
  f2->SetLineWidth(4);
  f3->SetLineWidth(4);
  f4->SetLineWidth(4);
  f1->SetLineStyle(2);
  f2->SetLineStyle(2);
  f3->SetLineStyle(2);
  f4->SetLineStyle(2);
  f1->Draw("same");
  f2->Draw("same");
  f3->Draw("same");
  f4->Draw("same");
  gPad->SetLogz();
  
  c1->cd(2);
  diff2d_CDC_CDH_pip->SetTitle("");
  double maxpim = diff2d_CDC_CDH_pim->GetMaximum();
  diff2d_CDC_CDH_pip->SetMaximum(maxpim);
  diff2d_CDC_CDH_pip->Draw("colz");
  TF1 *f5 = new TF1("f5",pipleftp,-0.5625,0,1);
  TF1 *f6 = new TF1("f6",pipleftm,-0.5625,0,1);
  TF1 *f7 = new TF1("f7",piprightp,0,0.7533,1);
  TF1 *f8 = new TF1("f8",piprightm,0,0.7533,1);
  f5->SetLineWidth(4);
  f6->SetLineWidth(4);
  f7->SetLineWidth(4);
  f8->SetLineWidth(4);
  f5->SetLineStyle(2);
  f6->SetLineStyle(2);
  f7->SetLineStyle(2);
  f8->SetLineStyle(2);
  f5->Draw("same");
  f6->Draw("same");
  f7->Draw("same");
  f8->Draw("same");
  gPad->SetLogz();

  c1->Print("Isolationcut.pdf");
  
  //convert phi radian to cm, and display projection (?, trial)
  TH2D* diff2d_CDC_CDH_pim_philen = new TH2D("diff2d_CDC_CDH_pim_philen","",100,-1.0*TMath::Pi()*CDHR,TMath::Pi()*CDHR,100,-100,100);
  for(int ix=0;ix<diff2d_CDC_CDH_pim->GetNbinsX();ix++){
    for(int iy=0;iy<diff2d_CDC_CDH_pim->GetNbinsY();iy++){
      double cont = diff2d_CDC_CDH_pim->GetBinContent(ix,iy);
      diff2d_CDC_CDH_pim_philen->SetBinContent(ix,iy,cont);
    }
  }
  TH2D* diff2d_CDC_CDH_pip_philen = new TH2D("diff2d_CDC_CDH_pip_philen","",100,-1.0*TMath::Pi()*CDHR,TMath::Pi()*CDHR,100,-100,100);
  for(int ix=0;ix<diff2d_CDC_CDH_pip->GetNbinsX();ix++){
    for(int iy=0;iy<diff2d_CDC_CDH_pip->GetNbinsY();iy++){
      double cont = diff2d_CDC_CDH_pip->GetBinContent(ix,iy);
      diff2d_CDC_CDH_pip_philen->SetBinContent(ix,iy,cont);
    }
  }
  diff2d_CDC_CDH_pim_philen->GetYaxis()->SetNdivisions(5,5,0,kTRUE);
  diff2d_CDC_CDH_pip_philen->GetYaxis()->SetNdivisions(5,5,0,kTRUE);
  diff2d_CDC_CDH_pim_philen->GetZaxis()->SetNdivisions(5,5,0,kTRUE);
  diff2d_CDC_CDH_pip_philen->GetZaxis()->SetNdivisions(5,5,0,kTRUE);
  diff2d_CDC_CDH_pim_philen->SetXTitle("CDH hit - #pi^{-} track (phi) [cm]");
  diff2d_CDC_CDH_pim_philen->GetXaxis()->CenterTitle();
  diff2d_CDC_CDH_pim_philen->GetYaxis()->CenterTitle();
  diff2d_CDC_CDH_pip_philen->GetXaxis()->CenterTitle();
  diff2d_CDC_CDH_pip_philen->GetYaxis()->CenterTitle();
  diff2d_CDC_CDH_pip_philen->SetXTitle("CDH hit - #pi^{+} track (phi) [cm]");
  diff2d_CDC_CDH_pim_philen->SetYTitle("CDH hit - #pi^{-} track (z) [cm]");
  diff2d_CDC_CDH_pip_philen->SetYTitle("CDH hit - #pi^{+} track (z) [cm]");
  
  
  TCanvas *cm = new TCanvas("cm","cm",1600,1700);
  cm->Divide(2,2,0,0);
  cm->cd(3);
  //diff2d_CDC_CDH_pim_philen->GetXaxis()->SetLabelSize(0.04);
  //diff2d_CDC_CDH_pim_philen->GetYaxis()->SetLabelSize(0.04);
  //diff2d_CDC_CDH_pim_philen->GetXaxis()->SetLabelOffset(1.4);
  //diff2d_CDC_CDH_pim_philen->GetYaxis()->SetLabelOffset(1.6);
  diff2d_CDC_CDH_pim_philen->GetXaxis()->SetTitleSize(0.04);
  diff2d_CDC_CDH_pim_philen->GetYaxis()->SetTitleSize(0.04);
  diff2d_CDC_CDH_pim_philen->GetXaxis()->SetTitleOffset(1.2);
  diff2d_CDC_CDH_pim_philen->GetYaxis()->SetTitleOffset(1.2);

  diff2d_CDC_CDH_pim_philen->Draw("col");
  TF1 *f1len = new TF1("f1len",pimleftp_len,-0.85425*CDHR,0,1);
  TF1 *f2len = new TF1("f2len",pimleftm_len,-0.85425*CDHR,0,1);
  TF1 *f3len = new TF1("f3len",pimrightp_len,0,0.85427*CDHR,1);
  TF1 *f4len = new TF1("f4len",pimrightm_len,0,0.85427*CDHR,1);
  f1len->SetLineWidth(4);
  f2len->SetLineWidth(4);
  f3len->SetLineWidth(4);
  f4len->SetLineWidth(4);
  f1len->SetLineStyle(2);
  f2len->SetLineStyle(2);
  f3len->SetLineStyle(2);
  f4len->SetLineStyle(2);
  f1len->Draw("same");
  f2len->Draw("same");
  f3len->Draw("same");
  f4len->Draw("same");
  gPad->SetLogz();
  cm->cd(1);
  TH1D* hpimx = diff2d_CDC_CDH_pim_philen->ProjectionX();
  hpimx->SetMarkerStyle(20);
  hpimx->SetMarkerSize(1.1);
  hpimx->Draw("E");
  cm->cd(4);
  TH1D* hpimy = diff2d_CDC_CDH_pim_philen->ProjectionY();//->Draw();
  TGraphErrors *gy = new TGraphErrors();
  for(int ix=0;ix<hpimy->GetNbinsX();ix++){
    double x= hpimy->GetXaxis()->GetBinCenter(ix);
    double cont = hpimy->GetBinContent(ix);
    gy->AddPoint(cont,x);
  }
  gy->SetLineColor(1);
  gy->SetMarkerStyle(20);
  gy->SetMarkerSize(1.2);
  gy->Draw("Ap");


  TCanvas *cp = new TCanvas("cp","cp",1800,1200);
  cp->Divide(2,2,0.,0.);
  cp->cd(3);
  double maxpimlen = diff2d_CDC_CDH_pim_philen->GetMaximum();
  diff2d_CDC_CDH_pip_philen->SetMaximum(maxpimlen);
  diff2d_CDC_CDH_pip_philen->GetXaxis()->SetTitleSize(0.04);
  diff2d_CDC_CDH_pip_philen->GetYaxis()->SetTitleSize(0.04);
  diff2d_CDC_CDH_pip_philen->GetXaxis()->SetTitleOffset(1.2);
  diff2d_CDC_CDH_pip_philen->GetYaxis()->SetTitleOffset(1.2);
  //diff2d_CDC_CDH_pip_philen->Draw("colz");
  TF1 *f5len = new TF1("f5len",pipleftp_len,-0.5625*CDHR,0,1);
  TF1 *f6len = new TF1("f6len",pipleftm_len,-0.5625*CDHR,0,1);
  TF1 *f7len = new TF1("f7len",piprightp_len,0,0.7536*CDHR,1);
  TF1 *f8len = new TF1("f8len",piprightm_len,0,0.7536*CDHR,1);
  f5len->SetLineWidth(4);
  f6len->SetLineWidth(4);
  f7len->SetLineWidth(4);
  f8len->SetLineWidth(4);
  f5len->SetLineStyle(2);
  f6len->SetLineStyle(2);
  f7len->SetLineStyle(2);
  f8len->SetLineStyle(2);
  f5len->Draw("same");
  f6len->Draw("same");
  f7len->Draw("same");
  f8len->Draw("same");
  gPad->SetLogz();

 
   
  TCanvas *c2 = new TCanvas("c2","c2",1600,600);
  c2->Divide(2,1);
  c2->cd(1);
  diff2d_CDC_CDH_pim_philen->Draw("col");
  f1len->Draw("same");
  f2len->Draw("same");
  f3len->Draw("same");
  f4len->Draw("same");
  gPad->SetLogz();
  c2->cd(2);
  diff2d_CDC_CDH_pip_philen->Draw("colz");
  f5len->Draw("same");
  f6len->Draw("same");
  f7len->Draw("same");
  f8len->Draw("same");
  gPad->SetLogz();
}






