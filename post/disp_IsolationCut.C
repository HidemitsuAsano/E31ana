#include "anacuts.h"

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

void disp_IsolationCut()
{
  //before subtracting mixed events
  //TFile *_file0 = TFile::Open("evanaIMpisigma_npippim_v245_out_dE2_iso_nostop.root");
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
  diff2d_CDC_CDH_pip->SetYTitle("CDH hit - #pi^{+} track (z) [cm]");
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
}






