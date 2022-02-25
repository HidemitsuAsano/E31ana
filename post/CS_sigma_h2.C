#include "anacuts.h"

void CS_sigma_h2()
{
  TH1::SetDefaultSumw2();
  TFile *f[3][3];
   
  const int dEcut[3]={2,4,6};
  const int sysud[3]={0,1,-1};
  for(int iEcut=0;iEcut<3;iEcut++){
    for(int isys=0;isys<3;isys++){
      f[iEcut][isys] = TFile::Open(Form("evanaIMsigma_npi_h2_v9_out_dE%d_iso_nostop_sub_sys%d.root",dEcut[iEcut],sysud[isys]));
    }
  }
  //TFile *f = TFile::Open("evanaIMsigma_npi_h2_v7_out_iso_nostop_sub.root");
  
  
  TH2F* Cospicm_IMnpip_pi[3][3];
  TH2F* Cospicm_IMnpim_pi[3][3];  
   
  for(int iEcut=0;iEcut<3;iEcut++){
    for(int isys=0;isys<3;isys++){
      Cospicm_IMnpip_pi[iEcut][isys]  = (TH2F*)f[iEcut][isys]->Get("Cospicm_IMnpip_pi");
      Cospicm_IMnpim_pi[iEcut][isys]  = (TH2F*)f[iEcut][isys]->Get("Cospicm_IMnpim_pi");
    }
  }
  
  TFile *facc[3];
  for(int iEcut=0;iEcut<3;iEcut++){
    facc[iEcut] = TFile::Open(Form("accH2dE%d.root",dEcut[iEcut]),"READ");
  }
  TH1D* accSp[3];
  TH1D* accSm[3]; 
  TH1D* accSp2[3];
  TH1D* accSm2[3];
  for(int iEcut=0;iEcut<3;iEcut++){
    accSp[iEcut] = (TH1D*)facc[iEcut]->Get("accCosSp");
    accSm[iEcut] = (TH1D*)facc[iEcut]->Get("accCosSm");
    accSp2[iEcut] = (TH1D*)facc[iEcut]->Get("accCosSp2");
    accSm2[iEcut] = (TH1D*)facc[iEcut]->Get("accCosSm2");
  }

  const int Splow = Cospicm_IMnpip_pi[0][0]->GetXaxis()->FindBin(anacuts::Sigmap_MIN);
  const int Sphigh = Cospicm_IMnpip_pi[0][0]->GetXaxis()->FindBin(anacuts::Sigmap_MAX);
  const int Smlow = Cospicm_IMnpim_pi[0][0]->GetXaxis()->FindBin(anacuts::Sigmam_MIN);
  const int Smhigh = Cospicm_IMnpim_pi[0][0]->GetXaxis()->FindBin(anacuts::Sigmam_MAX);
  
  TH1D* CS_Sp[3][3];   
  TH1D* CS_Sm[3][3];
  TH1D* CS_Sp2[3][3]; 
  TH1D* CS_Sm2[3][3]; 
  for(int iEcut=0;iEcut<3;iEcut++){
    for(int isys=0;isys<3;isys++){
      CS_Sp[iEcut][isys] = (TH1D*)Cospicm_IMnpip_pi[iEcut][isys]->ProjectionY(Form("CS_Sp%d%d",iEcut,isys),Splow,Sphigh);    
      CS_Sm[iEcut][isys] = (TH1D*)Cospicm_IMnpim_pi[iEcut][isys]->ProjectionY(Form("CS_Sm%d%d",iEcut,isys),Splow,Sphigh);    
      CS_Sp2[iEcut][isys] = (TH1D*)f[iEcut][isys]->Get("Cospicm_pi_Sp");    
      CS_Sm2[iEcut][isys] = (TH1D*)f[iEcut][isys]->Get("Cospicm_pi_Sm");    
      CS_Sp[iEcut][isys]->RebinX(5);
      CS_Sm[iEcut][isys]->RebinX(5);
      CS_Sp2[iEcut][isys]->RebinX(5);
      CS_Sm2[iEcut][isys]->RebinX(5);
      //CS_Sp2->Print("all");
      //CS_Sm2->Print("all");
      CS_Sp[iEcut][isys]->Divide(accSp[iEcut]);
      CS_Sm[iEcut][isys]->Divide(accSm[iEcut]);
      CS_Sp2[iEcut][isys]->Divide(accSp2[iEcut]);
      CS_Sm2[iEcut][isys]->Divide(accSm2[iEcut]);
    }
  }
  const double Lumi=379.8;
  TCanvas *cCS_Sp = new TCanvas("cCS_Sp","cCS_Sp",1000,800);
  //CS_Sp->RebinX(5);
  double CospiBinW = CS_Sp[0][0]->GetBinWidth(2);
  double CospiBinW2 = CS_Sp2[0][0]->GetBinWidth(2);
  std::cout << "Bin width " << CospiBinW2 << std::endl;
  const double cosToStrbin = 2.*3.1415*(CospiBinW); 
  const double cosToStrbin2 = 2.*3.1415*(CospiBinW2); 
  
  TGraphAsymmErrors* gCS_Sp[3][3];
  TGraphAsymmErrors* gCS_Sp2[3][3];
  for(int iEcut=0;iEcut<3;iEcut++){
    for(int isys=0;isys<3;isys++){
      CS_Sp[iEcut][isys]->Scale(1./cosToStrbin/Lumi);
      CS_Sp2[iEcut][isys]->Scale(1./cosToStrbin2/Lumi);
      CS_Sp[iEcut][isys]->GetXaxis()->SetRangeUser(0.3,1);
      gCS_Sp[iEcut][isys] = new TGraphAsymmErrors(CS_Sp[iEcut][isys]);
      gCS_Sp[iEcut][isys]->GetXaxis()->SetRangeUser(0.3,1);
      gCS_Sp2[iEcut][isys] = new TGraphAsymmErrors(CS_Sp2);
      gCS_Sp2[iEcut][isys]->GetXaxis()->SetRangeUser(0.3,1);
      if(iECut==0 && isys == 0)gCS_Sp2[iEcut][isys]->Draw("AP*");
    }
  }
  TCanvas *cCS_Sm = new TCanvas("cCS_Sm","cCS_Sm",1000,800);
  //CS_Sm->RebinX(5);
  CS_Sm->Scale(1./cosToStrbin/Lumi);
  CS_Sm2->Scale(1./cosToStrbin2/Lumi);
  CS_Sm->GetXaxis()->SetRangeUser(0.3,1);
  //CS_Sm->Draw("E");
  TGraphAsymmErrors* gCS_Sm = new TGraphAsymmErrors(CS_Sm);
  gCS_Sm->GetXaxis()->SetRangeUser(0.3,1);
  //gCS_Sm->Draw("AP");
  TGraphAsymmErrors* gCS_Sm2 = new TGraphAsymmErrors(CS_Sm2);
  gCS_Sm2->GetXaxis()->SetRangeUser(0.3,1);
  gCS_Sm2->Draw("AP*");


  TFile *fout = new TFile("CSsigma_H2.root","RECREATE");
  fout->cd();
  gCS_Sp->Write();
  gCS_Sm->Write();
  gCS_Sp2->Write();
  gCS_Sm2->Write();

}
