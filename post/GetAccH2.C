#include "anacuts.h"

void GetAccH2()
{
  TH1::SetDefaultSumw2();
  TFile *fSp = TFile::Open("../simpost/simIMsigma_H2_Sppim_npi_v4_out_iso_rej_nostop.root");
  TFile *fSm = TFile::Open("../simpost/simIMsigma_H2_Smpip_npi_v4_out_iso_rej_nostop.root");
  TFile *fGenSp = TFile::Open("../simpost/simIMsigma_H2_Sppim_v4.root");
  TFile *fGenSm = TFile::Open("../simpost/simIMsigma_H2_Smpip_v4.root");

  TH2F* Cospicm_IMnpip_pi_Sp = (TH2F*)fSp->Get("Cospicm_IMnpip_pi");//0.02 cos bin 
  TH2F* Cospicm_IMnpim_pi_Sm = (TH2F*)fSm->Get("Cospicm_IMnpim_pi");//0.02 cos bin

  TH1F* CosGenSp = (TH1F*)fGenSp->Get("ReactCosCM_0");//0.02 cos bin
  TH1F* CosGenSm = (TH1F*)fGenSm->Get("ReactCosCM_0");//0.02 cos bin

  TH1F* CosGenSp2 = new TH1F("CosGenSp2","CosGenSp2",50,0,1);
  TH1F* CosGenSm2 = new TH1F("CosGenSm2","CosGenSm2",50,0,1);
  
  TH1F* BLAnaPassedSp = (TH1F*)fGenSp->Get("BLAnaPassed");
  TH1F* BLAnaPassedSm = (TH1F*)fGenSm->Get("BLAnaPassed");
  double SimBeamSurvivalOKSp = BLAnaPassedSp->GetBinContent(2);//passed
  double SimBeamSurvivalFailSp = BLAnaPassedSp->GetBinContent(1);//not passed
  double SimBeamSurvivalRateSp = SimBeamSurvivalOKSp / (SimBeamSurvivalOKSp+SimBeamSurvivalFailSp);
  double SimBeamSurvivalOKSm = BLAnaPassedSm->GetBinContent(2);//passed
  double SimBeamSurvivalFailSm = BLAnaPassedSm->GetBinContent(1);//not passed
  double SimBeamSurvivalRateSm = SimBeamSurvivalOKSm / (SimBeamSurvivalOKSm+SimBeamSurvivalFailSm);
  
  for(int ibin=0;ibin<51;ibin++){
    double cont = CosGenSp->GetBinContent(ibin+50);
    CosGenSp2->SetBinContent(ibin,cont);
    double cont2 = CosGenSm->GetBinContent(ibin+50);
    CosGenSm2->SetBinContent(ibin,cont2);
  }
  std::cout << CosGenSp->GetBinCenter(51) << std::endl;
  std::cout << CosGenSp2->GetBinCenter(1) << std::endl;
  std::cout << Cospicm_IMnpip_pi_Sp->GetYaxis()->GetBinCenter(1) << std::endl;

  const int Splow = Cospicm_IMnpip_pi_Sp->GetXaxis()->FindBin(anacuts::Sigmap_MIN);
  const int Sphigh = Cospicm_IMnpip_pi_Sp->GetXaxis()->FindBin(anacuts::Sigmap_MAX);
  const int Smlow = Cospicm_IMnpim_pi_Sm->GetXaxis()->FindBin(anacuts::Sigmam_MIN);
  const int Smhigh = Cospicm_IMnpim_pi_Sm->GetXaxis()->FindBin(anacuts::Sigmam_MAX);
  
  TH1D* Cospicm_pi_Sp = (TH1D*)Cospicm_IMnpip_pi_Sp->ProjectionY("Cospicm_piSp",Splow,Sphigh);
  TH1D* Cospicm_pi_Sm = (TH1D*)Cospicm_IMnpim_pi_Sm->ProjectionY("Cospicm_piSm",Smlow,Smhigh);
    
  
  Cospicm_pi_Sm->Print("base");
  CosGenSm2->Print("base");
  
  TCanvas *cgenSp = new TCanvas("cgenSp","cgenSp");
  CosGenSp2->SetMinimum(0);
  CosGenSp2->Draw();
  TCanvas *cgenSm = new TCanvas("cgenSm","cgenSm");
  CosGenSm2->SetMinimum(0);
  CosGenSm2->Draw();

  TH1D* accCosSp = (TH1D*)Cospicm_pi_Sp->Clone("accCosSp");
  TH1D* accCosSm = (TH1D*)Cospicm_pi_Sm->Clone("accCosSm");
  TH1D* accCosSp2 = (TH1D*)fSp->Get("Cospicm_pi_Sp");
  TH1D* accCosSm2 = (TH1D*)fSm->Get("Cospicm_pi_Sm");
  TCanvas *ctest1 = new TCanvas("ctes1","ctest1");
  Cospicm_pi_Sp->RebinX(5);
  Cospicm_pi_Sp->Draw("E");
  TCanvas *ctest2 = new TCanvas("ctes2","ctest2");
  Cospicm_pi_Sm->RebinX(5);
  Cospicm_pi_Sm->Draw("E");
  std::cout << "Beam Survival Rate Sp " << SimBeamSurvivalRateSp << std::endl;
  std::cout << "Beam Survival Rate Sm " << SimBeamSurvivalRateSm << std::endl;
  accCosSp->Scale(1./SimBeamSurvivalRateSp);
  accCosSm->Scale(1./SimBeamSurvivalRateSm);
  accCosSp2->Scale(1./SimBeamSurvivalRateSp);
  accCosSm2->Scale(1./SimBeamSurvivalRateSm);
  
  accCosSp->RebinX(5);
  accCosSm->RebinX(5);
  accCosSp2->RebinX(5);
  accCosSm2->RebinX(5);
  CosGenSp2->RebinX(5);
  CosGenSm2->RebinX(5);
  
  accCosSp->Divide(accCosSp,CosGenSp2,1.,1.,"b");
  accCosSm->Divide(accCosSm,CosGenSm2,1.,1.,"b");
  accCosSp2->Divide(accCosSp2,CosGenSp2,1.,1.,"b");
  accCosSm2->Divide(accCosSm2,CosGenSm2,1.,1.,"b");

  TCanvas *c1 = new TCanvas("c1","c1");
  accCosSp->Draw("E");
  accCosSp2->SetLineColor(2);
  accCosSp2->Draw("Esame");

  TCanvas *c2 = new TCanvas("c2","c2");
  accCosSm->Draw("E");
  accCosSm2->SetLineColor(2);
  accCosSm2->Draw("Esame");

  TFile *fout = new TFile("accH2.root","RECREATE");
  fout->cd();
  accCosSp->SetName("accCosSp");
  accCosSm->SetName("accCosSm");
  accCosSp->Write();
  accCosSm->Write();
  accCosSp2->SetName("accCosSp2");
  accCosSm2->SetName("accCosSm2");
  accCosSp2->Write();
  accCosSm2->Write();


}
