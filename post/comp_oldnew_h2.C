#include "anacuts.h"

void comp_oldnew_h2()
{
  gStyle->SetOptStat(0);
  TFile *_file0 = TFile::Open("oldversion/evanaIMsigma_npi_h2_v4_out_iso_nostop_sub.root");
  TFile *_file1 = TFile::Open("evanaIMsigma_npi_h2_v14_out_dE2_iso_nostop_sub_sys0.root");

  TH2F* Cospicm_IMnpip_pi_old = (TH2F*)_file0->Get("Cospicm_IMnpip_pi");
  TH2F* Cospicm_IMnpip_pi_new = (TH2F*)_file1->Get("Cospicm_IMnpip_pi");
  TH2F* Cospicm_IMnpim_pi_old = (TH2F*)_file0->Get("Cospicm_IMnpim_pi");
  TH2F* Cospicm_IMnpim_pi_new = (TH2F*)_file1->Get("Cospicm_IMnpim_pi");
  Cospicm_IMnpip_pi_old->RebinY(5);
  Cospicm_IMnpip_pi_new->RebinY(5);
  Cospicm_IMnpim_pi_old->RebinY(5);
  Cospicm_IMnpim_pi_new->RebinY(5);

  const int Splow = Cospicm_IMnpip_pi_old->GetXaxis()->FindBin(anacuts::Sigmap_MIN);
  const int Sphigh = Cospicm_IMnpip_pi_old->GetXaxis()->FindBin(anacuts::Sigmap_MAX);
  const int Smlow = Cospicm_IMnpim_pi_old->GetXaxis()->FindBin(anacuts::Sigmam_MIN);
  const int Smhigh = Cospicm_IMnpim_pi_old->GetXaxis()->FindBin(anacuts::Sigmam_MAX);

  TCanvas *c1 = new TCanvas("c1","c1");
  TH1D* hcosnew = (TH1D*)Cospicm_IMnpip_pi_new->ProjectionY("hcosnew",Splow,Sphigh);
  hcosnew->SetLineColor(2);
  hcosnew->Draw("");

  TH1D* hcosold = (TH1D*)Cospicm_IMnpip_pi_old->ProjectionY("hcosold",Splow,Sphigh);
  hcosold->Draw("same");
  
  TCanvas *c2 = new TCanvas("c2","c2");
  TH1D* hcosnew_Sm = (TH1D*)Cospicm_IMnpim_pi_new->ProjectionY("hcosnew_Sm",Smlow,Smhigh);
  hcosnew_Sm->SetLineColor(2);
  hcosnew_Sm->Draw("");

  TH1D* hcosold_Sm = (TH1D*)Cospicm_IMnpim_pi_old->ProjectionY("hcosold_Sm",Smlow,Smhigh);
  hcosold_Sm->Draw("same");

  const int dEcut=2;
  TFile *fSp = TFile::Open(Form("../simpost/simIMsigma_H2_Sppim_npi_v15_out_dE%d_iso_rej_nostop.root",dEcut));
  TFile *fSm = TFile::Open(Form("../simpost/simIMsigma_H2_Smpip_npi_v15_out_dE%d_iso_rej_nostop.root",dEcut));
  TFile *fGenSp = TFile::Open("../simpost/simIMsigma_H2_Sppim_v15.root");
  TFile *fGenSm = TFile::Open("../simpost/simIMsigma_H2_Smpip_v15.root");
  //TFile *fSp = TFile::Open("../simpost/simIMsigma_H2_Sppim_npi_v3_out_iso_rej_nostop.root");
  //TFile *fSm = TFile::Open("../simpost/simIMsigma_H2_Smpip_npi_v3_out_iso_rej_nostop.root");
  //TFile *fGenSp = TFile::Open("../simpost/simIMsigma_H2_Sppim_v3.root");
  //TFile *fGenSm = TFile::Open("../simpost/simIMsigma_H2_Smpip_v3.root");
  TH2F* Cospicm_IMnpip_pi_Sp = (TH2F*)fSp->Get("Cospicm_IMnpip_pi");//0.02 cos bin 
  TH2F* Cospicm_IMnpim_pi_Sm = (TH2F*)fSm->Get("Cospicm_IMnpim_pi");//0.02 cos bin

  TH1F* CosGenSp = (TH1F*)fGenSp->Get("ReactCosCM_0");//0.02 cos bin
  TH1F* CosGenSm = (TH1F*)fGenSm->Get("ReactCosCM_0");//0.02 cos bin
  TH1F* BLAnaPassedSp = (TH1F*)fGenSp->Get("BLAnaPassed");
  TH1F* BLAnaPassedSm = (TH1F*)fGenSm->Get("BLAnaPassed");
  double SimBeamSurvivalOKSp = BLAnaPassedSp->GetBinContent(2);//passed
  double SimBeamSurvivalFailSp = BLAnaPassedSp->GetBinContent(1);//not passed
  double SimBeamSurvivalRateSp = SimBeamSurvivalOKSp / (SimBeamSurvivalOKSp+SimBeamSurvivalFailSp);
  double SimBeamSurvivalOKSm = BLAnaPassedSm->GetBinContent(2);//passed
  double SimBeamSurvivalFailSm = BLAnaPassedSm->GetBinContent(1);//not passed
  double SimBeamSurvivalRateSm = SimBeamSurvivalOKSm / (SimBeamSurvivalOKSm+SimBeamSurvivalFailSm);
  TH1D* Cospicm_pi_Sp = (TH1D*)Cospicm_IMnpip_pi_Sp->ProjectionY("Cospicm_piSp",Splow,Sphigh);
  TH1D* Cospicm_pi_Sm = (TH1D*)Cospicm_IMnpim_pi_Sm->ProjectionY("Cospicm_piSm",Smlow,Smhigh);

  Cospicm_pi_Sp->RebinX(5);
  Cospicm_pi_Sm->RebinX(5);
  TH1D* accCosSp = (TH1D*)Cospicm_pi_Sp->Clone("accCosSp");
  TH1D* accCosSm = (TH1D*)Cospicm_pi_Sm->Clone("accCosSm");
  
  accCosSp->Scale(1./SimBeamSurvivalRateSp);
  accCosSm->Scale(1./SimBeamSurvivalRateSm);
  
  
  TH1F* CosGenSp2 = new TH1F("CosGenSp2","CosGenSp2",50,0,1);
  TH1F* CosGenSm2 = new TH1F("CosGenSm2","CosGenSm2",50,0,1);
  
  for(int ibin=0;ibin<51;ibin++){
    double cont = CosGenSp->GetBinContent(ibin+50);
    double conterr = CosGenSp->GetBinError(ibin+50);
    CosGenSp2->SetBinContent(ibin,cont);
    CosGenSp2->SetBinError(ibin,conterr);
    double cont2 = CosGenSm->GetBinContent(ibin+50);
    double conterr2 = CosGenSm->GetBinError(ibin+50);
    CosGenSm2->SetBinContent(ibin,cont2);
    CosGenSm2->SetBinError(ibin,conterr2);
  }
  CosGenSp2->RebinX(5);
  CosGenSm2->RebinX(5);

  accCosSp->Divide(accCosSp,CosGenSp2,1.,1.,"b");
  accCosSm->Divide(accCosSm,CosGenSm2,1.,1.,"b");



  TFile *fSpold = TFile::Open("../simpost/simIMsigma_H2_Sppim_npi_v3_out_iso_rej_nostop.root");
  TFile *fSmold = TFile::Open("../simpost/simIMsigma_H2_Smpip_npi_v3_out_iso_rej_nostop.root");
  TFile *fGenSpold = TFile::Open("../simpost/simIMsigma_H2_Sppim_v3.root");
  TFile *fGenSmold = TFile::Open("../simpost/simIMsigma_H2_Smpip_v3.root");
  TH2F* Cospicm_IMnpip_pi_Spold = (TH2F*)fSpold->Get("Cospicm_IMnpip_pi");//0.02 cos bin 
  TH2F* Cospicm_IMnpim_pi_Smold = (TH2F*)fSmold->Get("Cospicm_IMnpim_pi");//0.02 cos bin

  TH1F* CosGenSpold = (TH1F*)fGenSpold->Get("ReactCosCM_0");//0.02 cos bin
  TH1F* CosGenSmold = (TH1F*)fGenSmold->Get("ReactCosCM_0");//0.02 cos bin
  TH1F* BLAnaPassedSpold = (TH1F*)fGenSpold->Get("BLAnaPassed");
  TH1F* BLAnaPassedSmold = (TH1F*)fGenSmold->Get("BLAnaPassed");
  double SimBeamSurvivalOKSpold = BLAnaPassedSpold->GetBinContent(2);//passed
  double SimBeamSurvivalFailSpold = BLAnaPassedSpold->GetBinContent(1);//not passed
  double SimBeamSurvivalRateSpold = SimBeamSurvivalOKSpold / (SimBeamSurvivalOKSpold+SimBeamSurvivalFailSpold);
  double SimBeamSurvivalOKSmold = BLAnaPassedSmold->GetBinContent(2);//passed
  double SimBeamSurvivalFailSmold = BLAnaPassedSmold->GetBinContent(1);//not passed
  double SimBeamSurvivalRateSmold = SimBeamSurvivalOKSmold / (SimBeamSurvivalOKSmold+SimBeamSurvivalFailSmold);
  TH1D* Cospicm_pi_Spold = (TH1D*)Cospicm_IMnpip_pi_Spold->ProjectionY("Cospicm_piSpold",Splow,Sphigh);
  TH1D* Cospicm_pi_Smold = (TH1D*)Cospicm_IMnpim_pi_Smold->ProjectionY("Cospicm_piSmold",Smlow,Smhigh);

  Cospicm_pi_Sp->SetLineColor(2);
  Cospicm_pi_Sm->SetLineColor(2);
  Cospicm_pi_Spold->RebinX(5);
  Cospicm_pi_Smold->RebinX(5);
  
  TH1D* accCosSpold = (TH1D*)Cospicm_pi_Spold->Clone("accCosSpold");
  TH1D* accCosSmold = (TH1D*)Cospicm_pi_Smold->Clone("accCosSmold");
  
  accCosSpold->Scale(1./SimBeamSurvivalRateSpold);
  accCosSmold->Scale(1./SimBeamSurvivalRateSmold);
  
  TH1F* CosGenSp2old = new TH1F("CosGenSp2old","CosGenSp2old",50,0,1);
  TH1F* CosGenSm2old = new TH1F("CosGenSm2old","CosGenSm2old",50,0,1);
  
  for(int ibin=0;ibin<51;ibin++){
    double cont = CosGenSpold->GetBinContent(ibin+50);
    double conterr = CosGenSpold->GetBinError(ibin+50);
    CosGenSp2old->SetBinContent(ibin,cont);
    CosGenSp2old->SetBinError(ibin,conterr);
    double cont2 = CosGenSmold->GetBinContent(ibin+50);
    double conterr2 = CosGenSmold->GetBinError(ibin+50);
    CosGenSm2old->SetBinContent(ibin,cont2);
    CosGenSm2old->SetBinError(ibin,conterr2);
  }
  CosGenSp2old->RebinX(5);
  CosGenSm2old->RebinX(5);

  accCosSpold->Divide(accCosSpold,CosGenSp2old,1.,1.,"b");
  accCosSmold->Divide(accCosSmold,CosGenSm2old,1.,1.,"b");

  
  TCanvas *c3 = new TCanvas("c3","c3");
  accCosSp->SetLineColor(2);
  accCosSp->Draw("");
  accCosSpold->Draw("same");

  TCanvas *c4 = new TCanvas("c4","c4");
  accCosSm->SetLineColor(2);
  accCosSm->Draw("");
  accCosSmold->Draw("same");

}
