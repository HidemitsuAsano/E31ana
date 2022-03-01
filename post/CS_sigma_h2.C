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
  TCanvas *cCS_Sp[3];
  //CS_Sp->RebinX(5);
  double CospiBinW = CS_Sp[0][0]->GetBinWidth(2);
  double CospiBinW2 = CS_Sp2[0][0]->GetBinWidth(2);
  std::cout << "Bin width " << CospiBinW2 << std::endl;
  const double cosToStrbin = 2.*3.1415*(CospiBinW); 
  const double cosToStrbin2 = 2.*3.1415*(CospiBinW2); 
  
  TGraphAsymmErrors* gCS_Sp[3][3];
  TGraphAsymmErrors* gCS_Sp2[3][3];
  for(int iEcut=0;iEcut<3;iEcut++){
    cCS_Sp[iEcut]= new TCanvas(Form("cCS_Sp%d",iEcut),Form("cCS_Sp%d",iEcut),1000,800);   
    cCS_Sp[iEcut]->cd();
    for(int isys=0;isys<3;isys++){
      CS_Sp[iEcut][isys]->Scale(1./cosToStrbin/Lumi);
      CS_Sp2[iEcut][isys]->Scale(1./cosToStrbin2/Lumi);
      CS_Sp[iEcut][isys]->GetXaxis()->SetRangeUser(0.3,1);
      gCS_Sp[iEcut][isys] = new TGraphAsymmErrors(CS_Sp[iEcut][isys]);
      gCS_Sp[iEcut][isys]->GetXaxis()->SetRangeUser(0.3,1);
      gCS_Sp2[iEcut][isys] = new TGraphAsymmErrors(CS_Sp2[iEcut][isys]);
      gCS_Sp2[iEcut][isys]->GetXaxis()->SetRangeUser(0.3,1);
      gCS_Sp2[iEcut][isys]->GetYaxis()->SetRangeUser(0,600);
      gCS_Sp2[iEcut][isys]->GetXaxis()->SetTitle("Cos#theta_{CM} miss-#pi^{-}");
      gCS_Sp2[iEcut][isys]->GetYaxis()->SetTitle("d#sigma/d#Omega [#mub/sr]");
      gCS_Sp2[iEcut][isys]->GetXaxis()->CenterTitle();
      gCS_Sp2[iEcut][isys]->GetYaxis()->CenterTitle();
      if(isys == 0)gCS_Sp2[iEcut][isys]->Draw("AP*");
      else                     gCS_Sp2[iEcut][isys]->Draw("p*");
    }
  }

  TCanvas *cCS_Sm[3]; 
  
  TGraphAsymmErrors* gCS_Sm[3][3];
  TGraphAsymmErrors* gCS_Sm2[3][3];
  for(int iEcut=0;iEcut<3;iEcut++){
    cCS_Sm[iEcut] = new TCanvas(Form("cCS_Sm%d",iEcut),Form("cCS_Sm%d",iEcut),1000,800);  
    cCS_Sm[iEcut]->cd();
    for(int isys=0;isys<3;isys++){
      CS_Sm[iEcut][isys]->Scale(1./cosToStrbin/Lumi);
      CS_Sm2[iEcut][isys]->Scale(1./cosToStrbin2/Lumi);
      CS_Sm[iEcut][isys]->GetXaxis()->SetRangeUser(0.3,1);
      gCS_Sm[iEcut][isys] = new TGraphAsymmErrors(CS_Sm[iEcut][isys]);
      gCS_Sm[iEcut][isys]->GetXaxis()->SetRangeUser(0.3,1);
      gCS_Sm2[iEcut][isys] = new TGraphAsymmErrors(CS_Sm2[iEcut][isys]);
      gCS_Sm2[iEcut][isys]->GetXaxis()->SetRangeUser(0.3,1);
      gCS_Sm2[iEcut][isys]->GetYaxis()->SetRangeUser(0,300);
      gCS_Sm2[iEcut][isys]->GetXaxis()->SetTitle("Cos#theta_{CM} miss-#pi^{+}");
      gCS_Sm2[iEcut][isys]->GetYaxis()->SetTitle("d#sigma/d#Omega [#mub/sr]");
      gCS_Sm2[iEcut][isys]->GetXaxis()->CenterTitle();
      gCS_Sm2[iEcut][isys]->GetYaxis()->CenterTitle();
      if(isys == 0) gCS_Sm2[iEcut][isys]->Draw("AP*");
      else                      gCS_Sm2[iEcut][isys]->Draw("P*");
    }
  }

  
  TGraphAsymmErrors *gCS_SpdE[3] ;//err:stat
  TGraphAsymmErrors *gCS_SmdE[3] ;//err:stat 
  TGraphAsymmErrors *gCS_SpdE_sys[3] ;//err: sys of mix.
  TGraphAsymmErrors *gCS_SmdE_sys[3] ;//err: sys of mix.
  
  TCanvas *cCS_SpdE[3];
  TCanvas *cCS_SmdE[3];
  for(int iEcut=0;iEcut<3;iEcut++){
    gCS_SpdE[iEcut] = (TGraphAsymmErrors*)gCS_Sp[iEcut][0]->Clone(Form("gCS_SpdE%d",iEcut));
    gCS_SmdE[iEcut] = (TGraphAsymmErrors*)gCS_Sm[iEcut][0]->Clone(Form("gCS_SmdE%d",iEcut));
    gCS_SpdE_sys[iEcut] = (TGraphAsymmErrors*) gCS_Sp[iEcut][0]->Clone(Form("gCS_SpdE_sys%d",iEcut));
    gCS_SmdE_sys[iEcut] = (TGraphAsymmErrors*) gCS_Sm[iEcut][0]->Clone(Form("gCS_SmdE_sys%d",iEcut));
    gCS_SpdE[iEcut]->GetXaxis()->SetTitle("Cos#theta_{CM} miss-#pi^{-}");
    gCS_SpdE[iEcut]->GetYaxis()->SetTitle("d#sigma/d#Omega [#mub/sr]");
    gCS_SpdE[iEcut]->GetXaxis()->CenterTitle();
    gCS_SpdE[iEcut]->GetYaxis()->CenterTitle();
    gCS_SmdE[iEcut]->GetXaxis()->SetTitle("Cos#theta_{CM} miss-#pi^{+}");
    gCS_SmdE[iEcut]->GetYaxis()->SetTitle("d#sigma/d#Omega [#mub/sr]");
    gCS_SmdE[iEcut]->GetXaxis()->CenterTitle();
    gCS_SmdE[iEcut]->GetYaxis()->CenterTitle();
    double *Xerr  = gCS_Sp[iEcut][0]->GetEXhigh();
    double *YcentSp = gCS_Sp[iEcut][0]->GetY();
    double *errYhiSp = gCS_Sp[iEcut][2]->GetY();
    double *errYloSp = gCS_Sp[iEcut][1]->GetY();
    double *YcentSm = gCS_Sm[iEcut][0]->GetY();
    double *errYhiSm = gCS_Sm[iEcut][2]->GetY();
    double *errYloSm = gCS_Sm[iEcut][1]->GetY();
    if(iEcut==0){
      //gCS_Sp[iEcut][0]->Print();
      //gCS_Sp[iEcut][1]->Print();
      //gCS_Sp[iEcut][2]->Print();
    }

    for(int ip=0;ip<gCS_SpdE_sys[iEcut]->GetN();ip++){
      gCS_SpdE_sys[iEcut]->SetPointEYhigh(ip,errYhiSp[ip]-YcentSp[ip]);
      gCS_SpdE_sys[iEcut]->SetPointEYlow(ip,YcentSp[ip]-errYloSp[ip]);
      gCS_SpdE_sys[iEcut]->SetPointEXhigh(ip,Xerr[ip]*0.5);
      gCS_SpdE_sys[iEcut]->SetPointEXlow(ip,Xerr[ip]*0.5);
    }
    cCS_SpdE[iEcut] = new TCanvas(Form("cCS_SpdE%d",iEcut),Form("cCS_SpdE%d",iEcut));
    gCS_SpdE[iEcut]->GetXaxis()->SetRangeUser(0.3,1);
    gCS_SpdE[iEcut]->GetYaxis()->SetRangeUser(0,600);
    gCS_SpdE[iEcut]->Draw("AP");
    gCS_SpdE_sys[iEcut]->SetLineColor(3);
    gCS_SpdE_sys[iEcut]->SetLineWidth(3);
    gCS_SpdE_sys[iEcut]->SetMarkerStyle(20);
    gCS_SpdE_sys[iEcut]->SetFillStyle(0);
    gCS_SpdE_sys[iEcut]->Draw("5 ");
    if(iEcut==0){
     // gCS_SpdE_sys[iEcut]->Print();
    }

    for(int ip=0;ip<gCS_SmdE_sys[iEcut]->GetN();ip++){
      gCS_SmdE_sys[iEcut]->SetPointEYhigh(ip,errYhiSm[ip]-YcentSm[ip]);
      gCS_SmdE_sys[iEcut]->SetPointEYlow(ip,YcentSm[ip]-errYloSm[ip]);
      gCS_SmdE_sys[iEcut]->SetPointEXhigh(ip,Xerr[ip]*0.5);
      gCS_SmdE_sys[iEcut]->SetPointEXlow(ip,Xerr[ip]*0.5);
    }
    cCS_SmdE[iEcut] = new TCanvas(Form("cCS_SmdE%d",iEcut),Form("cCS_SmdE%d",iEcut));
    gCS_SmdE[iEcut]->GetXaxis()->SetRangeUser(0.3,1);
    gCS_SmdE[iEcut]->GetYaxis()->SetRangeUser(0,300);
    gCS_SmdE[iEcut]->Draw("AP");
    gCS_SmdE_sys[iEcut]->SetLineColor(3);
    gCS_SmdE_sys[iEcut]->SetLineWidth(3);
    gCS_SmdE_sys[iEcut]->SetMarkerStyle(20);
    gCS_SmdE_sys[iEcut]->SetFillStyle(0);
    gCS_SmdE_sys[iEcut]->Draw("5 ");
  }
  
  
  TCanvas *cCS_SpdEsys = new TCanvas("cCS_SpdEsys","cCS_SpdEsys");
  gCS_SpdE[0]->Draw("AP");
  gCS_SpdE[1]->SetLineColor(2);
  gCS_SpdE[2]->SetLineColor(3);
  gCS_SpdE[1]->Draw("P");
  gCS_SpdE[2]->Draw("P");

  //gCS_SpdE[0]->Print();
  //gCS_SpdE[1]->Print();
  //gCS_SpdE[2]->Print();
  
  TCanvas *cCS_SmdEsys = new TCanvas("cCS_SmdEsys","cCS_SmdEsys");
  gCS_SmdE[0]->Draw("AP");
  gCS_SmdE[1]->SetLineColor(2);
  gCS_SmdE[2]->SetLineColor(3);
  gCS_SmdE[1]->Draw("P");
  gCS_SmdE[2]->Draw("P");
  

  
  TGraphAsymmErrors *gCS_SpdEsys = (TGraphAsymmErrors*)gCS_SpdE[0]->Clone("gCS_SpdEsys");
  double *Xerr  = gCS_SpdE[0]->GetEXhigh();
  double *YvalSp[3];
  for(int iEcut=0;iEcut<3;iEcut++){
    YvalSp[iEcut] = gCS_SpdE[iEcut]->GetY();
  }

  for(int ip=0;ip<gCS_SpdE[0]->GetN();ip++){
    double yeh = 0.00001;
    double yehcan1 = YvalSp[1][ip]-YvalSp[0][ip];
    //std::cout << ip << " "  << YvalSp[0][ip] << "  " << YvalSp[1][ip] << "  " << YvalSp[2][ip]  << std::endl;
    double yehcan2 = YvalSp[2][ip]-YvalSp[0][ip];
    if(yehcan1>yeh) yeh = yehcan1;
    if(yehcan2>yeh) yeh = yehcan2;
    
    double yel =0.00001;
    double yelcan1 = YvalSp[1][ip]-YvalSp[0][ip];
    double yelcan2 = YvalSp[2][ip]-YvalSp[0][ip];
    if(yel>yelcan1) yel=yelcan1;
    if(yel>yelcan2) yel=yelcan2;
    
    gCS_SpdEsys->SetPointEXhigh(ip,Xerr[ip]*0.2);
    gCS_SpdEsys->SetPointEXlow(ip,Xerr[ip]*0.2);
    gCS_SpdEsys->SetPointEYhigh(ip,yeh);
    gCS_SpdEsys->SetPointEYlow(ip,fabs(yel));
  }
  cCS_SpdEsys->cd();
  gCS_SpdEsys->SetLineColor(4);
  gCS_SpdEsys->SetFillStyle(0);
  //gCS_SpdEsys->Print();
  gCS_SpdEsys->Draw("5");


  TGraphAsymmErrors *gCS_SmdEsys = (TGraphAsymmErrors*)gCS_SmdE[0]->Clone("gCS_SmdEsys");
  double *YvalSm[3];
  for(int iEcut=0;iEcut<3;iEcut++){
    YvalSm[iEcut] = gCS_SmdE[iEcut]->GetY();
  }

  for(int ip=0;ip<gCS_SmdE[0]->GetN();ip++){
    double yeh = 0.00001;
    double yehcan1 = YvalSm[1][ip]-YvalSm[0][ip];
    //std::cout << ip << " "  << YvalSm[0][ip] << "  " << YvalSm[1][ip] << "  " << YvalSm[2][ip]  << std::endl;
    double yehcan2 = YvalSm[2][ip]-YvalSm[0][ip];
    if(yehcan1>yeh) yeh = yehcan1;
    if(yehcan2>yeh) yeh = yehcan2;
    
    double yel =0.00001;
    double yelcan1 = YvalSm[1][ip]-YvalSm[0][ip];
    double yelcan2 = YvalSm[2][ip]-YvalSm[0][ip];
    if(yel>yelcan1) yel=yelcan1;
    if(yel>yelcan2) yel=yelcan2;
    
    gCS_SmdEsys->SetPointEXhigh(ip,Xerr[ip]*0.2);
    gCS_SmdEsys->SetPointEXlow(ip,Xerr[ip]*0.2);
    gCS_SmdEsys->SetPointEYhigh(ip,yeh);
    gCS_SmdEsys->SetPointEYlow(ip,fabs(yel));
  }
  cCS_SmdEsys->cd();
  gCS_SmdEsys->SetLineColor(4);
  gCS_SmdEsys->SetFillStyle(0);
  //gCS_SmdEsys->Print();
  gCS_SmdEsys->Draw("5");

  
  TCanvas *cCS_SpToTal = new TCanvas("cCS_SpTotal","cCS_SpTotal");
  gCS_SpdE[0]->Draw("AP");
  gCS_SpdE_sys[0]->Draw("5 ");
  gCS_SpdEsys->Draw("5 ");
  
  TCanvas *cCS_SpToTal_syssum = new TCanvas("cCS_SpTotal_syssum","cCS_SpTotal_syssum");
  gCS_SpdE[0]->Draw("AP");
  TGraphAsymmErrors *gCS_Sp_syssum = (TGraphAsymmErrors*)gCS_SpdE_sys[0]->Clone("gCS_Sp_syssum");
  
  double *yehSp = gCS_Sp_syssum->GetEYhigh();
  double *yelSp = gCS_Sp_syssum->GetEYlow();
  double *yehSpadd = gCS_SpdEsys->GetEYhigh();
  double *yelSpadd = gCS_SpdEsys->GetEYlow();
  for(int ip=0;ip<gCS_Sp_syssum->GetN();ip++){
    double eh = sqrt(yehSp[ip]*yehSp[ip]+yehSpadd[ip]*yehSpadd[ip]);
    double el = sqrt(yelSp[ip]*yelSp[ip]+yelSpadd[ip]*yelSpadd[ip]);
    gCS_Sp_syssum->SetPointEYhigh(ip,eh);
    gCS_Sp_syssum->SetPointEYlow(ip,el);
  }
  gCS_SpdE[0]->Draw("AP");
  gCS_Sp_syssum->Draw("5");


  TCanvas *cCS_SmToTal = new TCanvas("cCS_SmTotal","cCS_SmTotal");
  gCS_SmdE[0]->Draw("AP");
  gCS_SmdE_sys[0]->Draw("5 ");
  gCS_SmdEsys->Draw("5 ");
 
  
  TCanvas *cCS_SmToTal_syssum = new TCanvas("cCS_SmTotal_syssum","cCS_SmTotal_syssum");
  gCS_SmdE[0]->Draw("AP");
  TGraphAsymmErrors *gCS_Sm_syssum = (TGraphAsymmErrors*)gCS_SmdE_sys[0]->Clone("gCS_Sm_syssum");
  
  double *yehSm = gCS_Sm_syssum->GetEYhigh();
  double *yelSm = gCS_Sm_syssum->GetEYlow();
  double *yehSmadd = gCS_SmdEsys->GetEYhigh();
  double *yelSmadd = gCS_SmdEsys->GetEYlow();
  for(int ip=0;ip<gCS_Sm_syssum->GetN();ip++){
    double eh = sqrt(yehSm[ip]*yehSm[ip]+yehSmadd[ip]*yehSmadd[ip]);
    double el = sqrt(yelSm[ip]*yelSm[ip]+yelSmadd[ip]*yelSmadd[ip]);
    gCS_Sm_syssum->SetPointEYhigh(ip,eh);
    gCS_Sm_syssum->SetPointEYlow(ip,el);
  }
  gCS_SmdE[0]->Draw("AP");
  gCS_Sm_syssum->Draw("5");
  
  TFile *fK0 = TFile::Open("CS_K0contami_H2.root","READ");
  TGraphErrors *gCSK0_Sp = (TGraphErrors*)fK0->Get("gCSK0_Sp");
  TGraphErrors *gCSK0_Sm = (TGraphErrors*)fK0->Get("gCSK0_Sm");
  double *xK0Sp = gCSK0_Sp->GetX();
  double *yK0Sp = gCSK0_Sp->GetY();
  TGraphAsymmErrors *gCS_Sp_K0sub = (TGraphAsymmErrors*)gCS_SpdE[0]->Clone("gCS_Sp_K0sub");
  TGraphAsymmErrors *gCS_Sp_syssum_K0sub = (TGraphAsymmErrors*)gCS_Sp_syssum->Clone("gCS_Sp_syssum_K0sub");
  double *ySp = gCS_Sp_K0sub->GetY();
  double *ySpsys = gCS_Sp_syssum_K0sub->GetY();
  for(int ip=0;ip<gCS_Sp_K0sub->GetN();ip++){
    gCS_Sp_K0sub->SetPoint(ip,xK0Sp[ip],ySp[ip]-yK0Sp[ip]);
    gCS_Sp_syssum_K0sub->SetPoint(ip,xK0Sp[ip],ySpsys[ip]-yK0Sp[ip]);
  }

  double *yK0Sm = gCSK0_Sm->GetY();
  TGraphAsymmErrors *gCS_Sm_K0sub = (TGraphAsymmErrors*)gCS_SmdE[0]->Clone("gCS_Sm_K0sub");
  TGraphAsymmErrors *gCS_Sm_syssum_K0sub = (TGraphAsymmErrors*)gCS_Sm_syssum->Clone("gCS_Sm_syssum_K0sub");
  double *ySm = gCS_Sm_K0sub->GetY();
  double *ySmsys = gCS_Sm_syssum_K0sub->GetY();
  for(int ip=0;ip<gCS_Sm_K0sub->GetN();ip++){
    gCS_Sm_K0sub->SetPoint(ip,xK0Sp[ip],ySm[ip]-yK0Sm[ip]);
    gCS_Sm_syssum_K0sub->SetPoint(ip,xK0Sp[ip],ySmsys[ip]-yK0Sm[ip]);
  }


  TFile *fout = new TFile("CSsigma_H2.root","RECREATE");
  fout->cd();
  gCS_Sp_K0sub->Write();
  gCS_Sm_K0sub->Write();
  gCS_Sp_syssum_K0sub->Write();
  gCS_Sm_syssum_K0sub->Write();
  gCS_SpdEsys->Write();
  gCS_SmdEsys->Write();
  gCS_Sp_syssum->Write();
  gCS_Sm_syssum->Write();
  for(int iEcut=0;iEcut<3;iEcut++){
    gCS_SpdE[iEcut]->Write();
    gCS_SmdE[iEcut]->Write();
    gCS_SpdE_sys[iEcut]->Write();
    gCS_SmdE_sys[iEcut]->Write();
    
    for(int isys=0;isys<3;isys++){
      gCS_Sp[iEcut][isys]->Write();
      gCS_Sm[iEcut][isys]->Write();
      gCS_Sp2[iEcut][isys]->Write();
      gCS_Sm2[iEcut][isys]->Write();
    }
  }
}
