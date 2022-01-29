#include "anacuts.h"

void GetAccH2()
{
  TFile *fSp = TFile::Open("../simpost/simIMsigma_H2_Sppim_npi_v1_out_iso_rej_nostop.root");
  TFile *fSm = TFile::Open("../simpost/simIMsigma_H2_Smpip_npi_v1_out_iso_rej_nostop.root");
  TFile *fGetSp = TFile::Open("../simpost/simIMsigma_H2_Sppim_v1.root");
  TFile *fGetSm = TFile::Open("../simpost/simIMsigma_H2_Smpip_v1.root");

  TH2F* Cospicm_IMnpip_pi_Sp = (TH2F*)fSp->Get("Cospicm_IMnpip_pi");//0.02 cos bin 
  TH2F* Cospicm_IMnpim_pi_Sm = (TH2F*)fSm->Get("Cospicm_IMnpim_pi");//0.02 cos bin

  TH1F* CosGenSp = (TH1F*)fGetSp->Get("ReactCosCM_0");//0.02 cos bin
  TH1F* CosGenSm = (TH1F*)fGetSm->Get("ReactCosCM_0");//0.02 cos bin

  TH1F* CosGenSp2 = new TH1F("CosGenSp2","CosGenSp2",50,0,1);
  TH1F* CosGenSm2 = new TH1F("CosGenSm2","CosGenSm2",50,0,1);

  for(int ibin=0;ibin<50;ibin++){
    double cont = CosGenSp->GetBinContent(ibin+50);
    CosGenSp2->SetBinContent(ibin,cont);
    double cont = CosGenSm->GetBinContent(ibin+50);
    CosGenSm2->SetBinContent(ibin,cont);
  }


  TH1D* Cospicm_pi_Sp = (TH1D*)Cospicm_IMnpip_pi_Sp->ProjectionY("Cospicm_pi_Sp");
  TH1D* Cospicm_pi_Sm = (TH1D*)Cospicm_IMnpim_pi_Sm->ProjectionY("Cospicm_pi_Sm");
  
  Cospicm_pi_Sm->Print("base");
  CosGenSm2->Print("base");

  TH1D* accCosSp = (TH1D*)Cospicm_pi_Sp->Clone("accCosSp");
  TH1D* accCosSm = (TH1D*)Cospicm_pi_Sm->Clone("accCosSm");
  accCosSp->Divide(accCosSp,CosGenSp2,1.,1.,"b");
  accCosSm->Divide(accCosSm,CosGenSm2,1.,1.,"b");

  TCanvas *c1 = new TCanvas("c1","c1");
  accCosSp->Draw("E");

  TCanvas *c2 = new TCanvas("c2","c2");
  accCosSm->Draw("E");

  TFile *fout = new TFile("accH2.root","RECREATE");
  fout->cd();
  accCosSp->Write();
  accCosSm->Write();



}
