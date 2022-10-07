#include "../src/GlobalVariables.h"

void CS_finals()
{
  if(gROOT->GetVersionInt() < 60000){
    std::cout << "Use ROOT6 !!" << std::endl;
    return;
  }

  const int Version = 245;
  const int dEcut[3]={2,4,6};

  //gStyle->SetErrorX(1.);  
  //gStyle->SetEndErrorSize(10);  
  TFile *fpisigma[3][3];//dEcut, sysud
  for(int iEcut=0;iEcut<3;iEcut++){
    fpisigma[iEcut][1] = TFile::Open(Form("cs_pisigma_v%d_dE%d_sys0.root",Version,dEcut[iEcut]),"READ");
  }
  fpisigma[0][0] = TFile::Open(Form("cs_pisigma_v%d_dE2_sys-1.root",Version),"READ");
  fpisigma[0][2] = TFile::Open(Form("cs_pisigma_v%d_dE2_sys1.root",Version),"READ");

  gStyle->SetOptStat(0);
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  //gROOT->ForceStyle();
  TH2D* q_IMnpipi_Sp_cs[4][3][3];//iq,dEcut,sysud
  TH2D* q_IMnpipi_Sm_cs[4][3][3];
  TH2D* q_IMnpipi_K0_cs[4][3][3];
  TH2D* q_IMnpipi_SpSmSum[4][3][3];
  TH1D* IMnpipi_Sp_cs[4][3][3];
  TH1D* IMnpipi_Sm_cs[4][3][3];
  TH1D* IMnpipi_K0_cs[4][3][3];
  TH1D* IMnpipi_SpSmSum[4][3][3];
  TGraphAsymmErrors* gIMnpipi_Sp_cs_Etotal[4];//CS with bin by bin total error
  TGraphAsymmErrors* gIMnpipi_Sm_cs_Etotal[4];//CS with bin by bin total error
  TGraphAsymmErrors* gIMnpipi_K0_cs_Etotal[4];//CS with bin by bin total error
  TGraphAsymmErrors* gIMnpipi_SpSmSum_cs_Etotal[4];//CS with bin by bin total error
  TGraphAsymmErrors  *gDecoErrorSp_CS[4];
  TGraphAsymmErrors  *gDecoErrorSm_CS[4];
  TGraphAsymmErrors  *gDecoErrorK0_CS[4];  
  TGraphAsymmErrors  *gDecoErrorSpSm_CS[4];
  
  for(int iq=0;iq<4;iq++){
    for(int iEcut=0;iEcut<3;iEcut++){
      q_IMnpipi_Sp_cs[iq][iEcut][1] = (TH2D*)fpisigma[iEcut][1]->Get(Form("q_IMnpipi_Sp_cs%d_sys0",iq));
      q_IMnpipi_Sm_cs[iq][iEcut][1] = (TH2D*)fpisigma[iEcut][1]->Get(Form("q_IMnpipi_Sm_cs%d_sys0",iq));
      q_IMnpipi_K0_cs[iq][iEcut][1] = (TH2D*)fpisigma[iEcut][1]->Get(Form("q_IMnpipi_K0_cs%d_sys0",iq));
      q_IMnpipi_SpSmSum[iq][iEcut][1] = (TH2D*)fpisigma[iEcut][1]->Get(Form("q_IMnpipi_SpSmSum%d_sys0",iq));
      IMnpipi_Sp_cs[iq][iEcut][1] = (TH1D*)fpisigma[iEcut][1]->Get(Form("IMnpipi_Sp_cs%d_sys0",iq));
      IMnpipi_Sm_cs[iq][iEcut][1] = (TH1D*)fpisigma[iEcut][1]->Get(Form("IMnpipi_Sm_cs%d_sys0",iq));
      IMnpipi_K0_cs[iq][iEcut][1] = (TH1D*)fpisigma[iEcut][1]->Get(Form("IMnpipi_K0_cs%d_sys0",iq));
      IMnpipi_SpSmSum[iq][iEcut][1] = (TH1D*)fpisigma[iEcut][1]->Get(Form("IMnpipi_SpSmSum%d_sys0",iq));
    }
    q_IMnpipi_Sp_cs[iq][0][0] = (TH2D*)fpisigma[0][0]->Get(Form("q_IMnpipi_Sp_cs%d_sys0",iq));
    q_IMnpipi_Sm_cs[iq][0][0] = (TH2D*)fpisigma[0][0]->Get(Form("q_IMnpipi_Sm_cs%d_sys0",iq));
    q_IMnpipi_K0_cs[iq][0][0] = (TH2D*)fpisigma[0][0]->Get(Form("q_IMnpipi_K0_cs%d_sys0",iq));
    q_IMnpipi_SpSmSum[iq][0][0] = (TH2D*)fpisigma[0][0]->Get(Form("q_IMnpipi_SpSmSum%d_sys0",iq));
    IMnpipi_Sp_cs[iq][0][0] = (TH1D*)fpisigma[0][0]->Get(Form("IMnpipi_Sp_cs%d_sys0",iq));
    IMnpipi_Sm_cs[iq][0][0] = (TH1D*)fpisigma[0][0]->Get(Form("IMnpipi_Sm_cs%d_sys0",iq));
    IMnpipi_K0_cs[iq][0][0] = (TH1D*)fpisigma[0][0]->Get(Form("IMnpipi_K0_cs%d_sys0",iq));
    IMnpipi_SpSmSum[iq][0][0] = (TH1D*)fpisigma[0][0]->Get(Form("IMnpipi_SpSmSum%d_sys0",iq));
    
    q_IMnpipi_Sp_cs[iq][0][2] = (TH2D*)fpisigma[0][2]->Get(Form("q_IMnpipi_Sp_cs%d_sys0",iq));
    q_IMnpipi_Sm_cs[iq][0][2] = (TH2D*)fpisigma[0][2]->Get(Form("q_IMnpipi_Sm_cs%d_sys0",iq));
    q_IMnpipi_K0_cs[iq][0][2] = (TH2D*)fpisigma[0][2]->Get(Form("q_IMnpipi_K0_cs%d_sys0",iq));
    q_IMnpipi_SpSmSum[iq][0][2] = (TH2D*)fpisigma[0][2]->Get(Form("q_IMnpipi_SpSmSum%d_sys0",iq));
    IMnpipi_Sp_cs[iq][0][2] = (TH1D*)fpisigma[0][2]->Get(Form("IMnpipi_Sp_cs%d_sys0",iq));
    IMnpipi_Sm_cs[iq][0][2] = (TH1D*)fpisigma[0][2]->Get(Form("IMnpipi_Sm_cs%d_sys0",iq));
    IMnpipi_K0_cs[iq][0][2] = (TH1D*)fpisigma[0][2]->Get(Form("IMnpipi_K0_cs%d_sys0",iq));
    IMnpipi_SpSmSum[iq][0][2] = (TH1D*)fpisigma[0][2]->Get(Form("IMnpipi_SpSmSum%d_sys0",iq));
    
    gDecoErrorSp_CS[iq] = (TGraphAsymmErrors*)fpisigma[0][1]->Get(Form("Graph_from_IMnpipi_Sp_cs_single%d_sys0",iq));
    gDecoErrorSm_CS[iq] = (TGraphAsymmErrors*)fpisigma[0][1]->Get(Form("Graph_from_IMnpipi_Sm_cs_single%d_sys0",iq));
    gDecoErrorK0_CS[iq] = (TGraphAsymmErrors*)fpisigma[0][1]->Get(Form("Graph_from_IMnpipi_K0_cs_single%d_sys0",iq));
    gDecoErrorSpSm_CS[iq] = (TGraphAsymmErrors*)fpisigma[0][1]->Get(Form("Graph_from_IMnpipi_SpSmSum%d_sys0",iq));
  }

  std::cout << "FILE GET " << std::endl;
  TCanvas *cqM_Sp[4];
  TCanvas *cqM_Sm[4];
  TCanvas *cqM_K0[4];
  TFile *fnuSp = new TFile("../simpost/NumericalRootFinder_fine20_Sp.root");
  TMultiGraph *mgSp = (TMultiGraph*)fnuSp->Get("mg");
  
  TFile *fnuSm = new TFile("../simpost/NumericalRootFinder_fine20_Sm.root");
  TMultiGraph *mgSm = (TMultiGraph*)fnuSm->Get("mg");
  for(int iq=0;iq<4;iq++){
    cqM_Sp[iq] = new TCanvas(Form("cqM_Sp%d",iq),Form("cqM_Sp%d",iq),1000,800);
    q_IMnpipi_Sp_cs[iq][0][1]->SetMinimum(0);
    q_IMnpipi_Sp_cs[iq][0][1]->GetYaxis()->SetRangeUser(0,0.65);
    q_IMnpipi_Sp_cs[iq][0][1]->GetXaxis()->SetRangeUser(1.3,1.9);
    q_IMnpipi_Sp_cs[iq][0][1]->Draw("colz");
    mgSp->Draw("p");
    cqM_Sm[iq] = new TCanvas(Form("cqM_Sm%d",iq),Form("cqM_Sm%d",iq),1000,800);
    q_IMnpipi_Sm_cs[iq][0][1]->SetMinimum(0);
    q_IMnpipi_Sm_cs[iq][0][1]->GetYaxis()->SetRangeUser(0,0.65);
    q_IMnpipi_Sm_cs[iq][0][1]->GetXaxis()->SetRangeUser(1.3,1.9);
    q_IMnpipi_Sm_cs[iq][0][1]->Draw("colz");
    mgSm->Draw("p");

    cqM_K0[iq] = new TCanvas(Form("cqM_K0%d",iq),Form("cqM_K0%d",iq),1000,800);
    q_IMnpipi_K0_cs[iq][0][1]->SetMinimum(0);
    q_IMnpipi_K0_cs[iq][0][1]->GetYaxis()->SetRangeUser(0,0.65);
    q_IMnpipi_K0_cs[iq][0][1]->GetXaxis()->SetRangeUser(1.3,1.9);
    q_IMnpipi_K0_cs[iq][0][1]->Draw("colz");
  }

  //zero suppresion of gDecoError
  for(int iq=0;iq<4;iq++){
    double n = gDecoErrorSp_CS[iq]->GetN();
    double *yval = gDecoErrorSp_CS[iq]->GetEYlow();
    for(int ip=n;ip>=0;ip--){
      if(yval[ip]<=0.00001 ) gDecoErrorSp_CS[iq]->RemovePoint(ip);
    }
  }
  for(int iq=0;iq<4;iq++){
    double n = gDecoErrorSm_CS[iq]->GetN();
    double *yval = gDecoErrorSm_CS[iq]->GetEYlow();
    for(int ip=n;ip>=0;ip--){
      if(yval[ip]<=0.00001) gDecoErrorSm_CS[iq]->RemovePoint(ip);
    }
  }
  for(int iq=0;iq<4;iq++){
    double n = gDecoErrorK0_CS[iq]->GetN();
    double *yval = gDecoErrorK0_CS[iq]->GetEYlow();
    for(int ip=n;ip>=0;ip--){
      if(yval[ip]<=0.00001) gDecoErrorK0_CS[iq]->RemovePoint(ip);
    }
  }
  for(int iq=0;iq<4;iq++){
    double n = gDecoErrorSpSm_CS[iq]->GetN();
    double *yval = gDecoErrorSpSm_CS[iq]->GetEYlow();
    for(int ip=n;ip>=0;ip--){
      if(yval[ip]<=0.00001) gDecoErrorSpSm_CS[iq]->RemovePoint(ip);
    }
  }


  //Sys Error summary
  TGraphAsymmErrors *gMIXErrorSp_CS[4];
  TGraphAsymmErrors *gMIXErrorSm_CS[4];
  TGraphAsymmErrors *gMIXErrorK0_CS[4];
  TGraphAsymmErrors *gMIXErrorSpSm_CS[4];
  TGraphAsymmErrors *gdEErrorSp_CS[4];
  TGraphAsymmErrors *gdEErrorSm_CS[4];
  TGraphAsymmErrors *gdEErrorK0_CS[4];

  
  for(int iq=0;iq<4;iq++){
    gMIXErrorSp_CS[iq] = new TGraphAsymmErrors(IMnpipi_Sp_cs[iq][0][1]);
    gMIXErrorSm_CS[iq] = new TGraphAsymmErrors(IMnpipi_Sm_cs[iq][0][1]);
    gMIXErrorK0_CS[iq] = new TGraphAsymmErrors(IMnpipi_K0_cs[iq][0][1]);
    gMIXErrorSpSm_CS[iq] = new TGraphAsymmErrors(IMnpipi_SpSmSum[iq][0][1]);
    gMIXErrorSp_CS[iq]->SetName(Form("gMIXErrorSp_CS%d",iq));
    gMIXErrorSm_CS[iq]->SetName(Form("gMIXErrorSm_CS%d",iq));
    gMIXErrorK0_CS[iq]->SetName(Form("gMIXErrorK0_CS%d",iq));
    gMIXErrorSpSm_CS[iq]->SetName(Form("gMIXErrorSpSm_CS%d",iq));

    gdEErrorSp_CS[iq] = new TGraphAsymmErrors(IMnpipi_Sp_cs[iq][0][1]);
    gdEErrorSm_CS[iq] = new TGraphAsymmErrors(IMnpipi_Sm_cs[iq][0][1]);
    gdEErrorK0_CS[iq] = new TGraphAsymmErrors(IMnpipi_K0_cs[iq][0][1]);
    gdEErrorSp_CS[iq]->SetName(Form("gdEErrorSp_CS%d",iq)); 
    gdEErrorSm_CS[iq]->SetName(Form("gdEErrorSm_CS%d",iq)); 
    gdEErrorK0_CS[iq]->SetName(Form("gdEErrorK0_CS%d",iq)); 


    //shift MIX scale and evaluate systematics
    for(int ip=0;ip<( gMIXErrorSp_CS[iq]->GetN());ip++){
       
      double xe = IMnpipi_Sp_cs[iq][0][0]->GetBinWidth(ip+1);
      gMIXErrorSp_CS[iq]->SetPointEXhigh(ip,xe/4);
      gMIXErrorSp_CS[iq]->SetPointEXlow(ip,xe/4);
  
      double  valdown = IMnpipi_Sp_cs[iq][0][0]->GetBinContent(ip+1);
      double  valdef = IMnpipi_Sp_cs[iq][0][1]->GetBinContent(ip+1);
      double  valup  = IMnpipi_Sp_cs[iq][0][2]->GetBinContent(ip+1);
      
      double yeh = valup-valdef;
      double yel = valdown-valdef;
         
      gMIXErrorSp_CS[iq]->SetPointEYhigh(ip,fabs(yeh));
      gMIXErrorSp_CS[iq]->SetPointEYlow(ip,fabs(yel));
    }
    for(int ip=0;ip<( gMIXErrorSm_CS[iq]->GetN());ip++){
      
      double xe = IMnpipi_Sm_cs[iq][0][0]->GetBinWidth(ip+1);
      gMIXErrorSm_CS[iq]->SetPointEXhigh(ip,xe/4);
      gMIXErrorSm_CS[iq]->SetPointEXlow(ip,xe/4);
      
      double  valdown = IMnpipi_Sm_cs[iq][0][0]->GetBinContent(ip+1);
      double  valdef = IMnpipi_Sm_cs[iq][0][1]->GetBinContent(ip+1);
      double  valup  = IMnpipi_Sm_cs[iq][0][2]->GetBinContent(ip+1);
      
      double yeh = valup-valdef;
      double yel = valdown-valdef;
         
      gMIXErrorSm_CS[iq]->SetPointEYhigh(ip,fabs(yeh));
      gMIXErrorSm_CS[iq]->SetPointEYlow(ip,fabs(yel));
    }
    for(int ip=0;ip<( gMIXErrorK0_CS[iq]->GetN());ip++){
      double  valdown = IMnpipi_K0_cs[iq][0][0]->GetBinContent(ip+1);
      double  valdef = IMnpipi_K0_cs[iq][0][1]->GetBinContent(ip+1);
      double  valup  = IMnpipi_K0_cs[iq][0][2]->GetBinContent(ip+1);
      
      double yeh = valup-valdef;
      double yel = valdown-valdef;
         
      gMIXErrorK0_CS[iq]->SetPointEYhigh(ip,fabs(yeh));
      gMIXErrorK0_CS[iq]->SetPointEYlow(ip,fabs(yel));
    }
    
    for(int ip=0;ip<( gMIXErrorSpSm_CS[iq]->GetN());ip++){
      double xe = IMnpipi_SpSmSum[iq][0][0]->GetBinWidth(ip+1);
      gMIXErrorSpSm_CS[iq]->SetPointEXhigh(ip,xe/4);
      gMIXErrorSpSm_CS[iq]->SetPointEXlow(ip,xe/4);
      
      double  valdown = IMnpipi_SpSmSum[iq][0][0]->GetBinContent(ip+1);
      double  valdef = IMnpipi_SpSmSum[iq][0][1]->GetBinContent(ip+1);
      double  valup  = IMnpipi_SpSmSum[iq][0][2]->GetBinContent(ip+1);
      
      double yeh = valup-valdef;
      double yel = valdown-valdef;
         
      gMIXErrorSpSm_CS[iq]->SetPointEYhigh(ip,fabs(yeh));
      gMIXErrorSpSm_CS[iq]->SetPointEYlow(ip,fabs(yel));
    }

    //shift dE cut and estimate systematics
    for(int ip=0;ip<( gdEErrorSp_CS[iq]->GetN());ip++){
      double  valdE6 = IMnpipi_Sp_cs[iq][2][1]->GetBinContent(ip+1);
      double  valdE4 = IMnpipi_Sp_cs[iq][1][1]->GetBinContent(ip+1);
      double  valdef = IMnpipi_Sp_cs[iq][0][1]->GetBinContent(ip+1);
      
      double ye6 = valdE6-valdef;
      double ye4 = valdE4-valdef;
      
      if(fabs(ye6)>fabs(ye4) && ye6>0){
        gdEErrorSp_CS[iq]->SetPointEYhigh(ip,fabs(ye6));
        gdEErrorSp_CS[iq]->SetPointEYlow(ip,0);
      }
      if(fabs(ye6)>fabs(ye4) && ye6<0){
        gdEErrorSp_CS[iq]->SetPointEYhigh(ip,0);
        gdEErrorSp_CS[iq]->SetPointEYlow(ip,fabs(ye6));
      }
      
      if(fabs(ye6)<fabs(ye4) && ye4>0){
        gdEErrorSp_CS[iq]->SetPointEYhigh(ip,fabs(ye4));
        gdEErrorSp_CS[iq]->SetPointEYlow(ip,0);
      }

      if(fabs(ye6)<fabs(ye4) && ye4<0){
        gdEErrorSp_CS[iq]->SetPointEYhigh(ip,0);
        gdEErrorSp_CS[iq]->SetPointEYlow(ip,fabs(ye4));
      
      }
    }

    for(int ip=0;ip<( gdEErrorSm_CS[iq]->GetN());ip++){
      double  valdE6 = IMnpipi_Sm_cs[iq][2][1]->GetBinContent(ip+1);
      double  valdE4 = IMnpipi_Sm_cs[iq][1][1]->GetBinContent(ip+1);
      double  valdef = IMnpipi_Sm_cs[iq][0][1]->GetBinContent(ip+1);
      
      double ye6 = valdE6-valdef;
      double ye4 = valdE4-valdef;
      
      if(fabs(ye6)>fabs(ye4) && ye6>0) gdEErrorSm_CS[iq]->SetPointEYhigh(ip,fabs(ye6));
      if(fabs(ye6)>fabs(ye4) && ye6<0) gdEErrorSm_CS[iq]->SetPointEYlow(ip,fabs(ye6));
      if(fabs(ye6)<fabs(ye4) && ye4>0) gdEErrorSm_CS[iq]->SetPointEYhigh(ip,fabs(ye4));
      if(fabs(ye6)<fabs(ye4) && ye4<0) gdEErrorSm_CS[iq]->SetPointEYlow(ip,fabs(ye4));
    }
    for(int ip=0;ip<( gdEErrorK0_CS[iq]->GetN());ip++){
      double  valdE6 = IMnpipi_K0_cs[iq][2][1]->GetBinContent(ip+1);
      double  valdE4 = IMnpipi_K0_cs[iq][1][1]->GetBinContent(ip+1);
      double  valdef = IMnpipi_K0_cs[iq][0][1]->GetBinContent(ip+1);
      
      double ye6 = valdE6-valdef;
      double ye4 = valdE4-valdef;
      
      if(fabs(ye6)>fabs(ye4) && ye6>0) gdEErrorK0_CS[iq]->SetPointEYhigh(ip,fabs(ye6));
      if(fabs(ye6)>fabs(ye4) && ye6<0) gdEErrorK0_CS[iq]->SetPointEYlow(ip,fabs(ye6));
      if(fabs(ye6)<fabs(ye4) && ye4>0) gdEErrorK0_CS[iq]->SetPointEYhigh(ip,fabs(ye4));
      if(fabs(ye6)<fabs(ye4) && ye4<0) gdEErrorK0_CS[iq]->SetPointEYlow(ip,fabs(ye4));
    }
  }

  
  TFile *flpim = TFile::Open("CSLpimFit_calc.root","READ");
  TH1D* CS_S1385_ToSp[3][3];
  TH1D* CS_S1385_ToSm[3][3];
  TH1D* CS_S1385_ToSpSm[3][3];
  TGraphAsymmErrors *gr_S1385_ToSqlow  = (TGraphAsymmErrors*)flpim->Get("gr_S1385_ToSqlow");
  TGraphAsymmErrors *gr_S1385_ToSqhi  = (TGraphAsymmErrors*)flpim->Get("gr_S1385_ToSqhi");
  TGraphAsymmErrors *gr_S1385_ToSqlowSum  = (TGraphAsymmErrors*)flpim->Get("gr_S1385_ToSqlowSum");
  TGraphAsymmErrors *gr_S1385_ToSqhiSum  = (TGraphAsymmErrors*)flpim->Get("gr_S1385_ToSqhiSum");
  

  for(int iq=0;iq<3;iq++){
    for(int isys=0;isys<3;isys++){
      CS_S1385_ToSp[iq][isys] = (TH1D*)flpim->Get(Form("CS_S1385_ToSp%d_sys%d",iq,isys));
      CS_S1385_ToSm[iq][isys] = (TH1D*)flpim->Get(Form("CS_S1385_ToSm%d_sys%d",iq,isys));
      CS_S1385_ToSpSm[iq][isys] = (TH1D*)flpim->Get(Form("CS_S1385_ToSpSm%d_sys%d",iq,isys));
    }
  }

  TGraphAsymmErrors *gS1385ErrorSp[3];//S1385 -> pi-Sigma+  with syserror 
  TGraphAsymmErrors *gS1385ErrorSm[3];//S1385 -> pi+Sigma-
  TGraphAsymmErrors *gS1385ErrorSpSm[3];//S1385 -> Sigma+ Sigma-

  for(int iq=0;iq<3;iq++){
    gS1385ErrorSp[iq] = new TGraphAsymmErrors(CS_S1385_ToSp[iq][1]);//def
    gS1385ErrorSm[iq] = new TGraphAsymmErrors(CS_S1385_ToSm[iq][1]);//def
    gS1385ErrorSpSm[iq] = new TGraphAsymmErrors(CS_S1385_ToSpSm[iq][1]);//def
    
    
    for(int ip=0;ip<(gS1385ErrorSp[iq]->GetN());ip++){
      double valup = CS_S1385_ToSp[iq][2]->GetBinContent(ip+1); 
      double valdown = CS_S1385_ToSp[iq][0]->GetBinContent(ip+1);
      double valdef = CS_S1385_ToSp[iq][1]->GetBinContent(ip+1);

      double yeup = valup - valdef;
      double yedown = valdef - valdown;
   
      gS1385ErrorSp[iq]->SetPointEYhigh(ip,yeup);
      gS1385ErrorSp[iq]->SetPointEYlow(ip,yedown);
    }

    for(int ip=0;ip<(gS1385ErrorSm[iq]->GetN());ip++){
      double valup = CS_S1385_ToSm[iq][2]->GetBinContent(ip+1); 
      double valdown = CS_S1385_ToSm[iq][0]->GetBinContent(ip+1);
      double valdef = CS_S1385_ToSm[iq][1]->GetBinContent(ip+1);

      double yeup = valup - valdef;
      double yedown = valdef - valdown;

      gS1385ErrorSm[iq]->SetPointEYhigh(ip,yeup);
      gS1385ErrorSm[iq]->SetPointEYlow(ip,yedown);
    }

    for(int ip=0;ip<(gS1385ErrorSpSm[iq]->GetN());ip++){
      double valup = CS_S1385_ToSpSm[iq][2]->GetBinContent(ip+1); 
      double valdown = CS_S1385_ToSpSm[iq][0]->GetBinContent(ip+1);
      double valdef = CS_S1385_ToSpSm[iq][1]->GetBinContent(ip+1);

      double yeup = valup - valdef;
      double yedown = valdef - valdown;

      gS1385ErrorSpSm[iq]->SetPointEYhigh(ip,yeup);
      gS1385ErrorSpSm[iq]->SetPointEYlow(ip,yedown);
    }
    
  }//iq
  
  //std::cout << "Nbin " << gS1385ErrorSp[0]->GetN() << std::endl;

  //gS1385ErrorSp[1]->Print();
  //  CS_S1385_ToSp[1][1]->Print("all");
  //  CS_S1385_ToSp[1][0]->Print("all");
  //  CS_S1385_ToSp[1][2]->Print("all");
  
  std::cout << __LINE__ << std::endl;
  
  
  TCanvas *cSysSp[4];
  for(int iq=0;iq<4;iq++){  
    cSysSp[iq] = new TCanvas(Form("cSysSp%d",iq),Form("cSysSp%d",iq),1100,800);
    IMnpipi_Sp_cs[iq][0][1]->SetLineColor(1);
    IMnpipi_Sp_cs[iq][0][1]->SetMarkerColor(1);
    IMnpipi_Sp_cs[iq][0][1]->GetXaxis()->SetRangeUser(1.3,1.6);
    IMnpipi_Sp_cs[iq][0][1]->SetTitle("");
    IMnpipi_Sp_cs[iq][0][1]->Draw("E");
    gDecoErrorSp_CS[iq]->Draw("5");
    //gdEErrorSp_CS[iq]->SetFillStyle(0);
    //gdEErrorSp_CS[iq]->SetMarkerColor(2);
    //gdEErrorSp_CS[iq]->SetLineColor(2);
    //gdEErrorSp_CS[iq]->Draw("5");
    gMIXErrorSp_CS[iq]->SetFillStyle(3002);
    gMIXErrorSp_CS[iq]->SetFillColor(4);
    gMIXErrorSp_CS[iq]->SetMarkerColor(4);
    gMIXErrorSp_CS[iq]->SetLineColor(4);
    gMIXErrorSp_CS[iq]->Draw("3");
    if(iq<3){
      gS1385ErrorSp[iq]->SetFillStyle(3001);
      gS1385ErrorSp[iq]->SetFillColor(0);
      gS1385ErrorSp[iq]->SetMarkerColor(6);
      gS1385ErrorSp[iq]->SetLineColor(6);
      //CS_S1385_ToSp[iq][1]->SetLineColor(6);
      //gS1385ErrorSp[iq]->Draw("3");
      //gS1385ErrorSp[iq]->Draw("5");
    }
    if(iq==1){
      gr_S1385_ToSqlow->SetLineColor(5);
      gr_S1385_ToSqlow->SetFillStyle(3002);
      gr_S1385_ToSqlow->SetFillColor(0);
      gr_S1385_ToSqlow->Draw("5");
    }
    if(iq==2){
      gr_S1385_ToSqhi->SetLineColor(5);
      gr_S1385_ToSqhi->SetFillStyle(3002);
      gr_S1385_ToSqhi->SetFillColor(0);
      gr_S1385_ToSqhi->Draw("5");
    }
    TLine *p = new TLine(1.29,0,1.605,0);
    p->SetLineColor(1);
    //p->SetLineWidth(2.0);
    p->SetLineStyle(2);
    p->Draw();
  }



  TCanvas *cSysSm[4];
  for(int iq=0;iq<4;iq++){  
    cSysSm[iq] = new TCanvas(Form("cSysSm%d",iq),Form("cSysSm%d",iq),1100,800);
    IMnpipi_Sm_cs[iq][0][1]->SetLineColor(1);
    IMnpipi_Sm_cs[iq][0][1]->SetMarkerColor(1);
    IMnpipi_Sm_cs[iq][0][1]->GetXaxis()->SetRangeUser(1.3,1.6);
    //IMnpipi_Sm_cs[iq][0][1]->GetYaxis()->SetRangeUser(0,IMnpipi_Sm_cs[iq][0][1]->GetMaximum()*1.2 );
    //gDecoErrorSm_CS[iq]->GetXaxis()->SetRangeUser(1.3,1.6);
    //gDecoErrorSm_CS[iq]->GetYaxis()->SetRangeUser(0,IMnpipi_Sm_cs[iq][0][1]->GetMaximum()*1.2 );
    IMnpipi_Sm_cs[iq][0][1]->SetTitle("");
    IMnpipi_Sm_cs[iq][0][1]->Draw("E");
    gDecoErrorSm_CS[iq]->Draw("5");
    //gdEErrorSm_CS[iq]->SetFillStyle(0);
    //gdEErrorSm_CS[iq]->SetMarkerColor(2);
    //gdEErrorSm_CS[iq]->SetLineColor(2);
    //gdEErrorSm_CS[iq]->Draw("5");
    gMIXErrorSm_CS[iq]->SetFillStyle(3002);
    gMIXErrorSm_CS[iq]->SetFillColor(4);
    gMIXErrorSm_CS[iq]->SetMarkerColor(4);
    gMIXErrorSm_CS[iq]->SetLineColor(4);
    gMIXErrorSm_CS[iq]->Draw("3");
    if(iq==1){
      gr_S1385_ToSqlow->SetLineColor(5);
      gr_S1385_ToSqlow->SetFillStyle(3002);
      gr_S1385_ToSqlow->SetFillColor(0);
      gr_S1385_ToSqlow->Draw("5");
    }
    if(iq==2){
      gr_S1385_ToSqhi->SetLineColor(5);
      gr_S1385_ToSqhi->SetFillStyle(3002);
      gr_S1385_ToSqhi->SetFillColor(0);
      gr_S1385_ToSqhi->Draw("5");
    }
    /*
    if(iq<3){
      gS1385ErrorSm[iq]->SetFillStyle(3001);
      gS1385ErrorSm[iq]->SetFillColor(0);
      gS1385ErrorSm[iq]->SetMarkerColor(6);
      gS1385ErrorSm[iq]->SetLineColor(6);
      //CS_S1385_ToSm[iq][1]->SetLineColor(6);
      gS1385ErrorSm[iq]->Draw("5");
    }*/

    TLine *p = new TLine(1.29,0,1.605,0);
    p->SetLineColor(1);
    //p->SetLineWidth(2.0);
    p->SetLineStyle(2);
    p->Draw();
  }

  TCanvas *cSysK0[4];
  for(int iq=0;iq<4;iq++){  
    cSysK0[iq] = new TCanvas(Form("cSysK0%d",iq),Form("cSysK0%d",iq),1100,800);
    IMnpipi_K0_cs[iq][0][1]->SetLineColor(1);
    IMnpipi_K0_cs[iq][0][1]->SetMarkerColor(1);
    IMnpipi_K0_cs[iq][0][1]->GetXaxis()->SetRangeUser(1.3,1.6);
    IMnpipi_K0_cs[iq][0][1]->Draw("E");
    gDecoErrorK0_CS[iq]->Draw("5");
    //gdEErrorSm_CS[iq]->SetFillStyle(0);
    //gdEErrorSm_CS[iq]->SetMarkerColor(2);
    //gdEErrorSm_CS[iq]->SetLineColor(2);
    //gdEErrorSm_CS[iq]->Draw("5");
    gMIXErrorK0_CS[iq]->SetFillStyle(3002);
    gMIXErrorK0_CS[iq]->SetFillColor(4);
    gMIXErrorK0_CS[iq]->SetMarkerColor(4);
    gMIXErrorK0_CS[iq]->SetLineColor(4);
    gMIXErrorK0_CS[iq]->Draw("3");
  }

  //gStyle->SetTitleYOffset(0.6);
  TCanvas *cSysSpSmSum[4];
  for(int iq=0;iq<4;iq++){  
    cSysSpSmSum[iq] = new TCanvas(Form("cSysSpSmSum%d",iq),Form("cSysSpSmSum%d",iq),1100,800);
    IMnpipi_SpSmSum[iq][0][1]->SetLineColor(1);
    IMnpipi_SpSmSum[iq][0][1]->SetMarkerColor(1);
    IMnpipi_SpSmSum[iq][0][1]->GetXaxis()->SetRangeUser(1.3,1.6);
    IMnpipi_SpSmSum[iq][0][1]->GetXaxis()->SetTitle("IM(#pi#Sigma) [GeV/c^{2}]");
    IMnpipi_SpSmSum[iq][0][1]->SetTitle("");
    IMnpipi_SpSmSum[iq][0][1]->SetYTitle("d#sigma/dM  [#mub/MeV^{2}]");
    IMnpipi_SpSmSum[iq][0][1]->GetYaxis()->CenterTitle();
    IMnpipi_SpSmSum[iq][0][1]->Draw("E");
    gDecoErrorSpSm_CS[iq]->Draw("5");
    gMIXErrorSpSm_CS[iq]->SetFillStyle(3002);
    gMIXErrorSpSm_CS[iq]->SetFillColor(4);
    gMIXErrorSpSm_CS[iq]->SetMarkerColor(4);
    gMIXErrorSpSm_CS[iq]->SetLineColor(4);
    gMIXErrorSpSm_CS[iq]->Draw("3");
    if(iq==1){
      gr_S1385_ToSqlowSum->SetLineColor(5);
      gr_S1385_ToSqlowSum->SetFillStyle(3002);
      gr_S1385_ToSqlowSum->SetFillColor(0);
      gr_S1385_ToSqlowSum->Draw("5");
    }
    if(iq==2){
      gr_S1385_ToSqhiSum->SetLineColor(5);
      gr_S1385_ToSqhiSum->SetFillStyle(3002);
      gr_S1385_ToSqhiSum->SetFillColor(0);
      gr_S1385_ToSqhiSum->Draw("5");
    }
    /*
    if(iq<3){
      gS1385ErrorSpSm[iq]->SetFillStyle(3001);
      gS1385ErrorSpSm[iq]->SetFillColor(0);
      gS1385ErrorSpSm[iq]->SetMarkerColor(6);
      gS1385ErrorSpSm[iq]->SetLineColor(6);
      //CS_S1385_ToSm[iq][1]->SetLineColor(6);
      gS1385ErrorSpSm[iq]->Draw("5");
    }*/
    TLine *p = new TLine(1.29,0,1.605,0);
    p->SetLineColor(1);
    //p->SetLineWidth(2.0);
    p->SetLineStyle(2);
    p->Draw();
  }

  
  TCanvas *cSysdESp[4];
  for(int iq=0;iq<4;iq++){
    cSysdESp[iq] = new TCanvas(Form("cSysdESp%d",iq),Form("cSysdESp%d",iq),1000,800);
    IMnpipi_Sp_cs[iq][0][1]->SetLineColor(1);
    IMnpipi_Sp_cs[iq][0][1]->SetMarkerColor(1);
    IMnpipi_Sp_cs[iq][0][1]->Draw();
    IMnpipi_Sp_cs[iq][1][1]->SetLineColor(2);
    IMnpipi_Sp_cs[iq][1][1]->SetMarkerColor(2);
    IMnpipi_Sp_cs[iq][1][1]->Draw("same");
    IMnpipi_Sp_cs[iq][2][1]->SetLineColor(3);
    IMnpipi_Sp_cs[iq][2][1]->SetMarkerColor(3);
    IMnpipi_Sp_cs[iq][2][1]->Draw("same");
  }

  TCanvas *cSysdESm[4];
  for(int iq=0;iq<4;iq++){
    cSysdESm[iq] = new TCanvas(Form("cSysdESm%d",iq),Form("cSysdESm%d",iq),1000,800);
    IMnpipi_Sm_cs[iq][0][1]->SetLineColor(1);
    IMnpipi_Sm_cs[iq][0][1]->SetMarkerColor(1);
    IMnpipi_Sm_cs[iq][0][1]->Draw();
    IMnpipi_Sm_cs[iq][1][1]->SetLineColor(2);
    IMnpipi_Sm_cs[iq][1][1]->SetMarkerColor(2);
    IMnpipi_Sm_cs[iq][1][1]->Draw("same");
    IMnpipi_Sm_cs[iq][2][1]->SetLineColor(3);
    IMnpipi_Sm_cs[iq][2][1]->SetMarkerColor(3);
    IMnpipi_Sm_cs[iq][2][1]->Draw("same");
  }          


  TCanvas *cSysdEK0[4];
  for(int iq=0;iq<4;iq++){
    cSysdEK0[iq] = new TCanvas(Form("cSysdEK0%d",iq),Form("cSysdEK0%d",iq),1000,800);
    IMnpipi_K0_cs[iq][0][1]->SetLineColor(1);
    IMnpipi_K0_cs[iq][0][1]->SetMarkerColor(1);
    IMnpipi_K0_cs[iq][0][1]->Draw();
    IMnpipi_K0_cs[iq][1][1]->SetLineColor(2);
    IMnpipi_K0_cs[iq][1][1]->SetMarkerColor(2);
    IMnpipi_K0_cs[iq][1][1]->Draw("same");
    IMnpipi_K0_cs[iq][2][1]->SetLineColor(3);
    IMnpipi_K0_cs[iq][2][1]->SetMarkerColor(3);
    IMnpipi_K0_cs[iq][2][1]->Draw("same");
  }


  //const double solidAngleCoscut = 2.0*3.1415926535*(1.00-0.99657);//theta 0.5 = 
  const double solidAngleCoscut = 2.0*3.1415926535*(0.90-0.60);//theta 0.5 = 
  const double solidAngleReso = 0.0015;//
  std::cout << "solid Angle. Coscut " << solidAngleCoscut << std::endl;
  std::cout << "solid Angle. Err   " << solidAngleReso << std::endl;
  std::cout << "1/Sr          " << 1./solidAngleCoscut << std::endl;
  const double SrErrup = 1./(2.0*3.1415926535*(1.00-0.99657+solidAngleReso));
  const double SrErrdown = 1./(2.0*3.1415926535*(1.00-0.99657-solidAngleReso));
  std::cout << "1/Sr Err up " << SrErrup  << std::endl;
  std::cout << "1/Sr Err down" << SrErrdown << std::endl;
   
  TH1D *CS_Spcomp = (TH1D*)IMnpipi_Sp_cs[3][0][1]->Clone("CS_Spcomp");
  TH1D *CS_Smcomp = (TH1D*)IMnpipi_Sm_cs[3][0][1]->Clone("CS_Smcomp");
  CS_Spcomp->Scale(1.0/solidAngleCoscut);
  CS_Smcomp->Scale(1.0/solidAngleCoscut);
  CS_Spcomp->GetXaxis()->SetRangeUser(1.3,1.6);
  CS_Smcomp->GetXaxis()->SetRangeUser(1.3,1.6);
  TGraphAsymmErrors *CS_SpcompMIX = (TGraphAsymmErrors*)gMIXErrorSp_CS[3]->Clone("CS_SpcompMIX");
  TGraphAsymmErrors *CS_SpcompDeco = (TGraphAsymmErrors*)gDecoErrorSp_CS[3]->Clone("CS_SpcompDeco");
  TGraphAsymmErrors *CS_SpcompSrErr = new TGraphAsymmErrors(IMnpipi_Sp_cs[3][0][1]);
  CS_SpcompSrErr->SetName("CS_SpcompSrErr");
  TGraphAsymmErrors *CS_SmcompSrErr = new TGraphAsymmErrors(IMnpipi_Sm_cs[3][0][1]);
  CS_SmcompSrErr->SetName("CS_SmcompSrErr");
  TH2F *CS_S1385_ToSp_coscut[3][3];
  TH2F *CS_S1385_ToSm_coscut[3][3];
  TH2F *CS_S1385_ToSpSm_coscut[3][3];
  
  for(int iq=0;iq<3;iq++){
    for(int isys=0;isys<3;isys++){
      CS_S1385_ToSp_coscut[iq][isys] = (TH2F*)flpim->Get(Form("CS_S1385_ToSp_coscut%d_%d",iq,isys));
      CS_S1385_ToSm_coscut[iq][isys] = (TH2F*)flpim->Get(Form("CS_S1385_ToSm_coscut%d_%d",iq,isys));
      CS_S1385_ToSpSm_coscut[iq][isys] = (TH2F*)flpim->Get(Form("CS_S1385_ToSpSm_coscut%d_%d",iq,isys));
    }
  }

  TH1D* CS_S1385_ToSp_coscut_px[3][3];
  TH1D* CS_S1385_ToSm_coscut_px[3][3];
  TH1D* CS_S1385_ToSpSm_coscut_px[3][3];

  for(int iq=0;iq<3;iq++){
    for(int isys=0;isys<3;isys++){
      CS_S1385_ToSp_coscut_px[iq][isys] = (TH1D*)CS_S1385_ToSp_coscut[iq][isys]->ProjectionX(Form("CS_S1385ToSpcoscut_px%d_%d",iq,isys));
      CS_S1385_ToSm_coscut_px[iq][isys] = (TH1D*)CS_S1385_ToSm_coscut[iq][isys]->ProjectionX(Form("CS_S1385ToSmcoscut_px%d_%d",iq,isys));
      CS_S1385_ToSpSm_coscut_px[iq][isys] = (TH1D*)CS_S1385_ToSpSm_coscut[iq][isys]->ProjectionX(Form("CS_S1385ToSpSmcoscut_px%d_%d",iq,isys));
      double binqwidth = CS_S1385_ToSp_coscut_px[iq][isys]->GetYaxis()->GetBinWidth(1);
      CS_S1385_ToSp_coscut_px[iq][isys]->Scale(1.0/solidAngleCoscut);
      CS_S1385_ToSm_coscut_px[iq][isys]->Scale(1.0/solidAngleCoscut);
      CS_S1385_ToSpSm_coscut_px[iq][isys]->Scale(1.0/solidAngleCoscut);
    }
  }

  TGraphAsymmErrors *gCS_coscutSp[3];
  TGraphAsymmErrors *gCS_coscutSm[3];
  TGraphAsymmErrors *gCS_coscutSpSm[3];

  for(int iq=0;iq<3;iq++){
    gCS_coscutSp[iq] = new TGraphAsymmErrors(CS_S1385_ToSp_coscut_px[iq][1]);
    gCS_coscutSm[iq] = new TGraphAsymmErrors(CS_S1385_ToSm_coscut_px[iq][1]);
    gCS_coscutSpSm[iq] = new TGraphAsymmErrors(CS_S1385_ToSpSm_coscut_px[iq][1]);
    gCS_coscutSp[iq]->SetName(Form("gCS_coscutSp%d",iq));
    gCS_coscutSm[iq]->SetName(Form("gCS_coscutSm%d",iq));
    gCS_coscutSpSm[iq]->SetName(Form("gCS_coscutSpSm%d",iq));
  }

  
  for(int iq=0;iq<3;iq++){
    for(int ip=0;ip<(gCS_coscutSp[0]->GetN());ip++){
      double valdown = CS_S1385_ToSp_coscut_px[iq][0]->GetBinContent(ip+1);
      double valdef = CS_S1385_ToSp_coscut_px[iq][1]->GetBinContent(ip+1);
      double valup = CS_S1385_ToSp_coscut_px[iq][2]->GetBinContent(ip+1);

      double yeh = valup-valdef;
      double yel = valdown-valdef;
      gCS_coscutSp[iq]->SetPointEYhigh(ip,fabs(yeh));
      gCS_coscutSp[iq]->SetPointEYlow(ip,fabs(yel));
    }
  }

  TCanvas *cgCS_coscutSp = new TCanvas("cgCS_coscutSp","cgCS_coscutSp");
  gCS_coscutSp[0]->Draw("ap");

  for(int iq=0;iq<3;iq++){
    for(int ip=0;ip<(gCS_coscutSm[0]->GetN());ip++){
      double valdown = CS_S1385_ToSm_coscut_px[iq][0]->GetBinContent(ip+1);
      double valdef = CS_S1385_ToSm_coscut_px[iq][1]->GetBinContent(ip+1);
      double valup = CS_S1385_ToSm_coscut_px[iq][2]->GetBinContent(ip+1);

      double yeh = valup-valdef;
      double yel = valdown-valdef;
      gCS_coscutSm[iq]->SetPointEYhigh(ip,fabs(yeh));
      gCS_coscutSm[iq]->SetPointEYlow(ip,fabs(yel));
    }
  }

  TCanvas *cgCS_coscutSm = new TCanvas("cgCS_coscutSm","cgCS_coscutSm");
  gCS_coscutSm[0]->Draw("ap");

  for(int iq=0;iq<3;iq++){
    for(int ip=0;ip<(gCS_coscutSpSm[0]->GetN());ip++){
      double valdown = CS_S1385_ToSpSm_coscut_px[iq][0]->GetBinContent(ip+1);
      double valdef = CS_S1385_ToSpSm_coscut_px[iq][1]->GetBinContent(ip+1);
      double valup = CS_S1385_ToSpSm_coscut_px[iq][2]->GetBinContent(ip+1);

      double yeh = valup-valdef;
      double yel = valdown-valdef;
      gCS_coscutSpSm[iq]->SetPointEYhigh(ip,fabs(yeh));
      gCS_coscutSpSm[iq]->SetPointEYlow(ip,fabs(yel));
    }
  }


  for(int i=0;i<CS_SpcompMIX->GetN();i++){
    double x = CS_SpcompMIX->GetPointX(i);
    double y = CS_SpcompMIX->GetPointY(i);
    y = y/solidAngleCoscut;
    double ex = CS_SpcompMIX->GetErrorX(i);
    double eyh = CS_SpcompMIX->GetErrorYhigh(i);
    double eyl = CS_SpcompMIX->GetErrorYlow(i);
    eyh = eyh/solidAngleCoscut;
    eyl = eyl/solidAngleCoscut;

    CS_SpcompMIX->SetPoint(i,x,y);
    CS_SpcompMIX->SetPointEYhigh(i,eyh);
    CS_SpcompMIX->SetPointEYlow(i,eyl);
  }
  for(int i=0;i<CS_SpcompDeco->GetN();i++){
    double x = CS_SpcompDeco->GetPointX(i);
    double y = CS_SpcompDeco->GetPointY(i);
    y = y/solidAngleCoscut;
    double ex = CS_SpcompDeco->GetErrorX(i);
    double eyh = CS_SpcompDeco->GetErrorYhigh(i);
    double eyl = CS_SpcompDeco->GetErrorYlow(i);
    eyh = eyh/solidAngleCoscut;
    eyl = eyl/solidAngleCoscut;

    CS_SpcompDeco->SetPoint(i,x,y);
    CS_SpcompDeco->SetPointEYhigh(i,eyh);
    CS_SpcompDeco->SetPointEYlow(i,eyl);
  }
  
  for(int i=0;i<CS_SpcompSrErr->GetN();i++){
    double x = CS_SpcompSrErr->GetPointX(i);
    double y = CS_SpcompSrErr->GetPointY(i);
    y = y/solidAngleCoscut;
    double ex = CS_SpcompSrErr->GetErrorX(i);
    double eyh = CS_SpcompSrErr->GetPointY(i);
    double eyl = CS_SpcompSrErr->GetPointY(i);
    eyh = eyh*SrErrdown-y;
    eyl = y-eyl*SrErrup;

    CS_SpcompSrErr->SetPoint(i,x,y);
    CS_SpcompSrErr->SetPointEYhigh(i,eyh);
    CS_SpcompSrErr->SetPointEYlow(i,eyl);
  }
  //CS_Spcomp->Draw();
  //CS_Smcomp->SetLineColor(3);
  //CS_Smcomp->Draw("same");

  TFile *finouepS = TFile::Open("pSinoue_CS.root","READ");
  TGraphErrors *grinoueSp = (TGraphErrors*)finouepS->Get("gr_inoueSp");
  TGraphErrors *grinoueSpcs = (TGraphErrors*)finouepS->Get("gr_inoueSpcs");
  TGraphErrors *grinoueSm = (TGraphErrors*)finouepS->Get("gr_inoueSm");
  TGraphErrors *grinoueSmcs = (TGraphErrors*)finouepS->Get("gr_inoueSmcs");
  
  TCanvas *cthetacompSp = new TCanvas("cthetacompSp","cthetacompSp",1000,800);
  cthetacompSp->cd();
  //CS_IMppipi_p_wL_mc_coscut->Draw("HEsame");
  //CS_Spcomp->SetMaximum(35);
  CS_Spcomp->SetYTitle("d^{2}#rho/dM d#Omega [#mub/MeVsr]");
  //CS_Spcomp->GetYaxis()->SetRangeUser(-1,30);
  CS_Spcomp->Draw();
  //grinoueSpcs->Draw("p");
  CS_SpcompMIX->Draw("3");
  CS_SpcompSrErr->SetFillStyle(3002);
  CS_SpcompSrErr->SetFillColor(5);
  CS_SpcompSrErr->GetXaxis()->SetRangeUser(1.3,1.6);
  //CS_SpcompSrErr->Draw("3");
  CS_SpcompDeco->SetLineColor(3);
  CS_SpcompDeco->Draw("5");
  gCS_coscutSp[0]->SetFillStyle(3002);
  gCS_coscutSp[0]->SetFillColor(6);
  gCS_coscutSp[0]->GetXaxis()->SetRangeUser(1.3,1.6);
  gCS_coscutSp[0]->Draw("3");
  
  


  TGraphAsymmErrors *CS_SmcompMIX = (TGraphAsymmErrors*)gMIXErrorSm_CS[3]->Clone("CS_SmcompMIX");
  TGraphAsymmErrors *CS_SmcompDeco = (TGraphAsymmErrors*)gDecoErrorSm_CS[3]->Clone("CS_SmcompDeco");
  for(int i=0;i<CS_SmcompMIX->GetN();i++){
    double x = CS_SmcompMIX->GetPointX(i);
    double y = CS_SmcompMIX->GetPointY(i);
    y = y/solidAngleCoscut;
    double ex = CS_SmcompMIX->GetErrorX(i);
    double eyh = CS_SmcompMIX->GetErrorYhigh(i);
    double eyl = CS_SmcompMIX->GetErrorYlow(i);
    eyh = eyh/solidAngleCoscut;
    eyl = eyl/solidAngleCoscut;

    CS_SmcompMIX->SetPoint(i,x,y);
    CS_SmcompMIX->SetPointEYhigh(i,eyh);
    CS_SmcompMIX->SetPointEYlow(i,eyl);
  }
  for(int i=0;i<CS_SmcompDeco->GetN();i++){
    double x = CS_SmcompDeco->GetPointX(i);
    double y = CS_SmcompDeco->GetPointY(i);
    y = y/solidAngleCoscut;
    double ex = CS_SmcompDeco->GetErrorX(i);
    double eyh = CS_SmcompDeco->GetErrorYhigh(i);
    double eyl = CS_SmcompDeco->GetErrorYlow(i);
    eyh = eyh/solidAngleCoscut;
    eyl = eyl/solidAngleCoscut;

    CS_SmcompDeco->SetPoint(i,x,y);
    CS_SmcompDeco->SetPointEYhigh(i,eyh);
    CS_SmcompDeco->SetPointEYlow(i,eyl);
  }
  for(int i=0;i<CS_SmcompSrErr->GetN();i++){
    double x = CS_SmcompSrErr->GetPointX(i);
    double y = CS_SmcompSrErr->GetPointY(i);
    y = y/solidAngleCoscut;
    double ex = CS_SmcompSrErr->GetErrorX(i);
    double eyh = CS_SmcompSrErr->GetPointY(i);
    double eyl = CS_SmcompSrErr->GetPointY(i);
    eyh = eyh*SrErrdown-y;
    eyl = y-eyl*SrErrup;

    CS_SmcompSrErr->SetPoint(i,x,y);
    CS_SmcompSrErr->SetPointEYhigh(i,eyh);
    CS_SmcompSrErr->SetPointEYlow(i,eyl);
  }
 
  TCanvas *cthetacompSm = new TCanvas("cthetacompSm","cthetacompSm",1000,800);
  cthetacompSm->cd();
  //CS_IMppipi_p_wL_mc_coscut->Draw("HEsame");
  //CS_Smcomp->SetMaximum(10);
  CS_Smcomp->SetYTitle("d^{2}#rho/dMd#Omega [#mub/MeVsr]");
  //CS_Smcomp->GetYaxis()->SetRangeUser(-1,10);
  CS_Smcomp->Draw();
  CS_SmcompDeco->Draw("5");
  CS_SmcompSrErr->SetFillStyle(3002);
  CS_SmcompSrErr->SetFillColor(5);
  CS_SmcompSrErr->GetXaxis()->SetRangeUser(1.3,1.6);
  //CS_SmcompSrErr->Draw("3");
  CS_SmcompMIX->Draw("3");
  //grinoueSmcs->Draw("p");
  //grinoueSm->Draw("P");
  gCS_coscutSm[0]->SetFillStyle(3002);
  gCS_coscutSm[0]->SetFillColor(6);
  gCS_coscutSm[0]->GetXaxis()->SetRangeUser(1.3,1.6);
  gCS_coscutSm[0]->Draw("3");

  
  cSysSp[0]->SaveAs("csSp0.pdf","PDF");
  cSysSp[1]->SaveAs("csSp1.pdf","PDF");
  cSysSp[2]->SaveAs("csSp2.pdf","PDF");
  cSysSm[0]->SaveAs("csSm0.pdf","PDF");
  cSysSm[1]->SaveAs("csSm1.pdf","PDF");
  cSysSm[2]->SaveAs("csSm2.pdf","PDF");
  cSysSpSmSum[0]->SaveAs("csSpSm0.pdf","PDF");
  cSysSpSmSum[1]->SaveAs("csSpSm1.pdf","PDF");
  cSysSpSmSum[2]->SaveAs("csSpSm2.pdf","PDF");
  
  
  TCanvas *cSysSpEtotal[4];
  for(int iq=1;iq<3;iq++){  
    cSysSpEtotal[iq] = new TCanvas(Form("cSysSpEtotal%d",iq),Form("cSysSpEtotal%d",iq),1100,800);
    gIMnpipi_Sp_cs_Etotal[iq] = new TGraphAsymmErrors(IMnpipi_Sp_cs[iq][0][1]);
    gIMnpipi_Sp_cs_Etotal[iq]->SetName(Form("gIMnpipi_Sp_cs_Etotal%d",iq));
    gIMnpipi_Sp_cs_Etotal[iq]->SetMarkerColor(1);
    gIMnpipi_Sp_cs_Etotal[iq]->GetXaxis()->SetRangeUser(1.3,1.6);
    gIMnpipi_Sp_cs_Etotal[iq]->SetTitle("");
    gIMnpipi_Sp_cs_Etotal[iq]->GetXaxis()->SetTitle("IM(#pi^{-}#Sigma^{+}) [GeV/c^{2}]");
    gIMnpipi_Sp_cs_Etotal[iq]->GetXaxis()->CenterTitle();
    gIMnpipi_Sp_cs_Etotal[iq]->GetYaxis()->SetTitle("d#sigma/dM  [#mub/(MeV/c^{2})]");
    gIMnpipi_Sp_cs_Etotal[iq]->GetYaxis()->CenterTitle();
    if(iq==1){
      for(int ip=0;ip<10;ip++){
        gIMnpipi_Sp_cs_Etotal[iq]->RemovePoint(0);
        gMIXErrorSp_CS[iq]->RemovePoint(0);
      }
    }
    if(iq==2){
      for(int ip=0;ip<8;ip++){
        gIMnpipi_Sp_cs_Etotal[iq]->RemovePoint(0);
        gMIXErrorSp_CS[iq]->RemovePoint(0);
      }
    }
   
    double start1 = gIMnpipi_Sp_cs_Etotal[iq]->GetPointX(0);
    double start2 = gDecoErrorSp_CS[iq]->GetPointX(0);
    int diffstart = (start2-start1)/0.015 + 1 ;
    std::cout << diffstart << std::endl;

    for(int ip=0;ip<gIMnpipi_Sp_cs_Etotal[iq]->GetN()-diffstart;ip++){
      double statElow  = gIMnpipi_Sp_cs_Etotal[iq]->GetErrorYlow(ip+diffstart);
      double statEhigh = gIMnpipi_Sp_cs_Etotal[iq]->GetErrorYhigh(ip+diffstart);
      double decoElow  = gDecoErrorSp_CS[iq]->GetErrorYlow(ip);
      double decoEhigh = gDecoErrorSp_CS[iq]->GetErrorYhigh(ip);
      double Elow = sqrt(pow(statElow,2.0)+pow(decoElow,2.0));  
      double Ehigh = sqrt(pow(statEhigh,2.0)+pow(decoEhigh,2.0));  
      gIMnpipi_Sp_cs_Etotal[iq]->SetPointEYlow(ip+diffstart,Elow);
      gIMnpipi_Sp_cs_Etotal[iq]->SetPointEYhigh(ip+diffstart,Ehigh);
    }    
    gIMnpipi_Sp_cs_Etotal[iq]->Draw("ap");
    gMIXErrorSp_CS[iq]->SetFillStyle(3001);
    gMIXErrorSp_CS[iq]->SetFillColor(12);
    gMIXErrorSp_CS[iq]->SetMarkerColor(12);
    gMIXErrorSp_CS[iq]->SetLineColor(12);
    gMIXErrorSp_CS[iq]->Draw("5");
    if(iq==1){
      gr_S1385_ToSqlow->SetLineColor(3);
      gr_S1385_ToSqlow->SetFillStyle(3002);
      gr_S1385_ToSqlow->SetFillColor(0);
      gr_S1385_ToSqlow->Draw("5");
    }
    if(iq==2){
      gr_S1385_ToSqhi->SetLineColor(3);
      gr_S1385_ToSqhi->SetFillStyle(3002);
      gr_S1385_ToSqhi->SetFillColor(0);
      gr_S1385_ToSqhi->Draw("5");
    }
    TLine *p = new TLine(1.29,0,1.605,0);
    p->SetLineColor(1);
    //p->SetLineWidth(2.0);
    p->SetLineStyle(2);
    p->Draw();
  }
  
  TCanvas *cSysSmEtotal[4];
  for(int iq=1;iq<3;iq++){  
    cSysSmEtotal[iq] = new TCanvas(Form("cSysSmEtotal%d",iq),Form("cSysSmEtotal%d",iq),1100,800);
    gIMnpipi_Sm_cs_Etotal[iq] = new TGraphAsymmErrors(IMnpipi_Sm_cs[iq][0][1]);
    gIMnpipi_Sm_cs_Etotal[iq]->SetName(Form("gIMnpipi_Sm_cs_Etotal%d",iq));
    gIMnpipi_Sm_cs_Etotal[iq]->SetMarkerColor(1);
    gIMnpipi_Sm_cs_Etotal[iq]->GetXaxis()->SetRangeUser(1.3,1.6);
    gIMnpipi_Sm_cs_Etotal[iq]->SetTitle("");
    gIMnpipi_Sm_cs_Etotal[iq]->GetXaxis()->SetTitle("IM(#pi^{+}#Sigma^{-}) [GeV/c^{2}]");
    gIMnpipi_Sm_cs_Etotal[iq]->GetXaxis()->CenterTitle();
    gIMnpipi_Sm_cs_Etotal[iq]->GetYaxis()->SetTitle("d#sigma/dM  [#mub/(MeV/c^{2})]");
    gIMnpipi_Sm_cs_Etotal[iq]->GetYaxis()->CenterTitle();
    if(iq==1){
      for(int ip=0;ip<10;ip++){
        gIMnpipi_Sm_cs_Etotal[iq]->RemovePoint(0);
        gMIXErrorSm_CS[iq]->RemovePoint(0);
      }
    }
    if(iq==2){
      for(int ip=0;ip<8;ip++){
        gIMnpipi_Sm_cs_Etotal[iq]->RemovePoint(0);
        gMIXErrorSm_CS[iq]->RemovePoint(0);
      }
    }
   
    double start1 = gIMnpipi_Sm_cs_Etotal[iq]->GetPointX(0);
    double start2 = gDecoErrorSm_CS[iq]->GetPointX(0);
    int diffstart = (start2-start1)/0.015 + 1 ;
    std::cout << diffstart << std::endl;
    
    const double dx=0.0015;
    for(int ip=0;ip<gIMnpipi_Sm_cs_Etotal[iq]->GetN()-diffstart;ip++){
      double x = gIMnpipi_Sm_cs_Etotal[iq]->GetPointX(ip);
      gIMnpipi_Sm_cs_Etotal[iq]->SetPointX(ip,x+dx);
      double xe = gIMnpipi_Sm_cs_Etotal[iq]->GetErrorXlow(ip);
      gIMnpipi_Sm_cs_Etotal[iq]->SetPointEXlow(ip,xe+dx);
      gIMnpipi_Sm_cs_Etotal[iq]->SetPointEXhigh(ip,xe-dx);
    }
    for(int ip=0;ip<gMIXErrorSm_CS[iq]->GetN();ip++){
      double x = gMIXErrorSm_CS[iq]->GetPointX(ip);
      gMIXErrorSm_CS[iq]->SetPointX(ip,x+dx);
      double xe = gMIXErrorSm_CS[iq]->GetErrorXlow(ip);
      gMIXErrorSm_CS[iq]->SetPointEXlow(ip,xe+dx);
      gMIXErrorSm_CS[iq]->SetPointEXhigh(ip,xe-dx);
    }

    for(int ip=0;ip<gIMnpipi_Sm_cs_Etotal[iq]->GetN()-diffstart;ip++){
      double statElow  = gIMnpipi_Sm_cs_Etotal[iq]->GetErrorYlow(ip+diffstart);
      double statEhigh = gIMnpipi_Sm_cs_Etotal[iq]->GetErrorYhigh(ip+diffstart);
      double decoElow  = gDecoErrorSm_CS[iq]->GetErrorYlow(ip);
      double decoEhigh = gDecoErrorSm_CS[iq]->GetErrorYhigh(ip);
      double Elow = sqrt(pow(statElow,2.0)+pow(decoElow,2.0));  
      double Ehigh = sqrt(pow(statEhigh,2.0)+pow(decoEhigh,2.0));  
      gIMnpipi_Sm_cs_Etotal[iq]->SetPointEYlow(ip+diffstart,Elow);
      gIMnpipi_Sm_cs_Etotal[iq]->SetPointEYhigh(ip+diffstart,Ehigh);
    }    
    gIMnpipi_Sm_cs_Etotal[iq]->Draw("ap");
    gMIXErrorSm_CS[iq]->SetFillStyle(3001);
    gMIXErrorSm_CS[iq]->SetFillColor(12);
    gMIXErrorSm_CS[iq]->SetMarkerColor(12);
    gMIXErrorSm_CS[iq]->SetLineColor(12);
    gMIXErrorSm_CS[iq]->Draw("5");
    if(iq==1){
      gr_S1385_ToSqlow->SetLineColor(3);
      gr_S1385_ToSqlow->SetFillStyle(3002);
      gr_S1385_ToSqlow->SetFillColor(0);
      gr_S1385_ToSqlow->Draw("5");
    }
    if(iq==2){
      gr_S1385_ToSqhi->SetLineColor(3);
      gr_S1385_ToSqhi->SetFillStyle(3002);
      gr_S1385_ToSqhi->SetFillColor(0);
      gr_S1385_ToSqhi->Draw("5");
    }
    TLine *p = new TLine(1.29,0,1.605,0);
    p->SetLineColor(1);
    //p->SetLineWidth(2.0);
    p->SetLineStyle(2);
    p->Draw();
  }

  IMnpipi_Sp_cs[1][0][1]->Print("all");
  gIMnpipi_Sp_cs_Etotal[1]->Print();
  gDecoErrorSp_CS[1]->Print();
  gMIXErrorSp_CS[1]->Print();

  TCanvas *cSysSpSmSumEtotal[4];
  for(int iq=1;iq<3;iq++){  
    cSysSpSmSumEtotal[iq] = new TCanvas(Form("cSysSpSmSumEtotal%d",iq),Form("cSysSpSmSumEtotal%d",iq),1100,1100);
    
    cSysSpSmSumEtotal[iq]->SetBottomMargin(0.12);
    cSysSpSmSumEtotal[iq]->SetLeftMargin(0.14);
    cSysSpSmSumEtotal[iq]->SetRightMargin(0.08);


    gIMnpipi_SpSmSum_cs_Etotal[iq] = new TGraphAsymmErrors(IMnpipi_SpSmSum[iq][0][1]);
    gIMnpipi_SpSmSum_cs_Etotal[iq]->SetName(Form("gIMnpipi_SpSmSum_cs_Etotal%d",iq));
    gIMnpipi_SpSmSum_cs_Etotal[iq]->SetMarkerColor(1);
    gIMnpipi_SpSmSum_cs_Etotal[iq]->GetXaxis()->SetRangeUser(1.3,1.6);
    gIMnpipi_SpSmSum_cs_Etotal[iq]->SetTitle("");
    gIMnpipi_SpSmSum_cs_Etotal[iq]->GetXaxis()->SetTitle("IM(#pi#Sigma) [GeV/c^{2}]");
    gIMnpipi_SpSmSum_cs_Etotal[iq]->GetXaxis()->CenterTitle();
    gIMnpipi_SpSmSum_cs_Etotal[iq]->GetYaxis()->SetTitle("d#sigma/dM  [#mub/(MeV/c^{2})]");
    gIMnpipi_SpSmSum_cs_Etotal[iq]->GetYaxis()->CenterTitle();

    gIMnpipi_SpSmSum_cs_Etotal[iq]->GetXaxis()->SetTitleSize(0.04);
    gIMnpipi_SpSmSum_cs_Etotal[iq]->GetYaxis()->SetTitleSize(0.04);
    gIMnpipi_SpSmSum_cs_Etotal[iq]->GetXaxis()->SetTitleOffset(1.4);
    gIMnpipi_SpSmSum_cs_Etotal[iq]->GetYaxis()->SetTitleOffset(1.6);
   if(iq==1){
      for(int ip=0;ip<10;ip++){
        gIMnpipi_SpSmSum_cs_Etotal[iq]->RemovePoint(0);
        gMIXErrorSpSm_CS[iq]->RemovePoint(0);
      }
    }
    if(iq==2){
      for(int ip=0;ip<8;ip++){
        gIMnpipi_SpSmSum_cs_Etotal[iq]->RemovePoint(0);
        gMIXErrorSpSm_CS[iq]->RemovePoint(0);
      }
    }
   
    double start1 = gIMnpipi_SpSmSum_cs_Etotal[iq]->GetPointX(0);
    double start2 = gDecoErrorSpSm_CS[iq]->GetPointX(0);
    int diffstart = (start2-start1)/0.015 + 1 ;
    std::cout << diffstart << std::endl;

    for(int ip=0;ip<gIMnpipi_SpSmSum_cs_Etotal[iq]->GetN()-diffstart;ip++){
      double statElow  = gIMnpipi_SpSmSum_cs_Etotal[iq]->GetErrorYlow(ip+diffstart);
      double statEhigh = gIMnpipi_SpSmSum_cs_Etotal[iq]->GetErrorYhigh(ip+diffstart);
      double decoElow  = gDecoErrorSpSm_CS[iq]->GetErrorYlow(ip);
      double decoEhigh = gDecoErrorSpSm_CS[iq]->GetErrorYhigh(ip);
      double Elow = sqrt(pow(statElow,2.0)+pow(decoElow,2.0));  
      double Ehigh = sqrt(pow(statEhigh,2.0)+pow(decoEhigh,2.0));  
      gIMnpipi_SpSmSum_cs_Etotal[iq]->SetPointEYlow(ip+diffstart,Elow);
      gIMnpipi_SpSmSum_cs_Etotal[iq]->SetPointEYhigh(ip+diffstart,Ehigh);
    }    
    gIMnpipi_SpSmSum_cs_Etotal[iq]->Draw("ap");
    gMIXErrorSpSm_CS[iq]->SetFillStyle(3001);
    //gMIXErrorSpSm_CS[iq]->SetFillColor(4);
    gMIXErrorSpSm_CS[iq]->SetFillColor(12);
    gMIXErrorSpSm_CS[iq]->SetMarkerColor(4);
    gMIXErrorSpSm_CS[iq]->SetLineColor(12);
    gMIXErrorSpSm_CS[iq]->Draw("5");
    if(iq==1){
      gr_S1385_ToSqlowSum->SetLineColor(5);
      gr_S1385_ToSqlowSum->SetFillStyle(3002);
      gr_S1385_ToSqlowSum->SetFillColor(0);
      gr_S1385_ToSqlowSum->Draw("5");
    }
    if(iq==2){
      gr_S1385_ToSqhiSum->SetLineColor(5);
      gr_S1385_ToSqhiSum->SetFillStyle(3002);
      gr_S1385_ToSqhiSum->SetFillColor(0);
      gr_S1385_ToSqhiSum->Draw("5");
    }
    TLine *p = new TLine(1.29,0,1.605,0);
    p->SetLineColor(1);
    //p->SetLineWidth(2.0);
    p->SetLineStyle(2);
    p->Draw();
  }
  cSysSpSmSumEtotal[1]->SaveAs("Sumqlow.pdf");
  cSysSpSmSumEtotal[2]->SaveAs("Sumqhi.pdf");

  
  TCanvas *cboth[4];
  for(int iq=1;iq<3;iq++){
    cboth[iq] = new TCanvas(Form("cboth%d",iq),Form("cboth%d",iq),1100,1100);
    cboth[iq]->SetBottomMargin(0.12);
    cboth[iq]->SetLeftMargin(0.14);
    cboth[iq]->SetRightMargin(0.08);
    gIMnpipi_Sp_cs_Etotal[iq]->SetLineColor(kGreen+2);
    gIMnpipi_Sp_cs_Etotal[iq]->SetFillColor(0);
    gIMnpipi_Sp_cs_Etotal[iq]->GetXaxis()->SetTitle("IM(#pi^{#pm}#Sigma^{#mp}) [GeV/c^{2}]");
    gIMnpipi_Sp_cs_Etotal[iq]->GetXaxis()->SetTitleSize(0.04);
    gIMnpipi_Sp_cs_Etotal[iq]->GetYaxis()->SetTitleSize(0.04);
    gIMnpipi_Sp_cs_Etotal[iq]->GetXaxis()->SetTitleOffset(1.4);
    gIMnpipi_Sp_cs_Etotal[iq]->GetYaxis()->SetTitleOffset(1.6);
    gIMnpipi_Sp_cs_Etotal[iq]->SetMarkerColor(kGreen+2);
    gMIXErrorSp_CS[iq]->SetLineColor(kGreen+2);
    gMIXErrorSp_CS[iq]->SetFillColor(0);
    gMIXErrorSp_CS[iq]->SetFillStyle(0);
    gMIXErrorSp_CS[iq]->SetMarkerColor(kGreen+2);
    gIMnpipi_Sp_cs_Etotal[iq]->Draw("ap");
    gMIXErrorSp_CS[iq]->Draw("5");
    gIMnpipi_Sm_cs_Etotal[iq]->SetLineColor(4);
    gIMnpipi_Sm_cs_Etotal[iq]->SetFillColor(0);
    gIMnpipi_Sm_cs_Etotal[iq]->SetMarkerColor(4);
    gIMnpipi_Sm_cs_Etotal[iq]->SetMarkerStyle(24);
    gMIXErrorSm_CS[iq]->SetLineColor(4);
    gMIXErrorSm_CS[iq]->SetFillColor(0);
    gMIXErrorSm_CS[iq]->SetFillStyle(0);
    gMIXErrorSm_CS[iq]->SetMarkerColor(4);
    gIMnpipi_Sm_cs_Etotal[iq]->Draw("p");
    gMIXErrorSm_CS[iq]->Draw("5");
    if(iq==1){
      gr_S1385_ToSqlow->SetLineColor(5);
      gr_S1385_ToSqlow->SetFillStyle(3002);
      gr_S1385_ToSqlow->SetFillColor(0);
      gr_S1385_ToSqlow->Draw("5");
    }
    if(iq==2){
      gr_S1385_ToSqhi->SetLineColor(5);
      gr_S1385_ToSqhi->SetFillStyle(3002);
      gr_S1385_ToSqhi->SetFillColor(0);
      gr_S1385_ToSqhi->Draw("5");
    }
    TLine *p = new TLine(1.29,0,1.605,0);
    p->SetLineColor(1);
    //p->SetLineWidth(2.0);
    p->SetLineStyle(2);
    p->Draw();
    TLine *pkp = new TLine(kpMass+pMass,0,kpMass+pMass,3);
    pkp->SetLineColor(1);
    pkp->SetLineStyle(2);
    pkp->Draw();
    
    TLine *pk0n = new TLine(k0Mass+nMass,0,k0Mass+nMass,3);
    pk0n->SetLineColor(1);
    pk0n->SetLineStyle(2);
    pk0n->Draw();

    //gIMnpipi_Sp_cs_Etotal[iq]->Draw("p");
  }
  cboth[1]->SaveAs("cbothqlow.pdf");
  cboth[2]->SaveAs("cbothqhi.pdf");

  std::cout << "k-p mass " <<   kpMass+pMass << std::endl;
  std::cout << "k0n mass " <<   k0Mass+nMass << std::endl;
  TCanvas *csumqboth = new TCanvas("csumqboth","csumqboth",1100,1100);
  csumqboth->SetBottomMargin(0.12);
  csumqboth->SetLeftMargin(0.14);
  csumqboth->SetRightMargin(0.08);
  gIMnpipi_SpSmSum_cs_Etotal[1]->Draw("ap");
  gIMnpipi_SpSmSum_cs_Etotal[2]->Draw("p");

  TLine *p = new TLine(1.29,0,1.605,0);
  p->SetLineColor(1);
  //p->SetLineWidth(2.0);
  p->SetLineStyle(2);
  p->Draw();

  TH2D* Cosn_IMnpipi_Sp_cs[3][3];//deco sys, mix sys
  TH2D* Cosn_IMnpipi_Sm_cs[3][3];//deco sys, mix sys
  TH2D* Cosn_IMnpipi_SpSmSum[3][3];//deco sys, mix sys
  TH1D* CosL1405[3][3];//deco sys, mix sys
   
  
  for(int idecosys=0;idecosys<3;idecosys++){
    for(int imixsys=0;imixsys<3;imixsys++){
      Cosn_IMnpipi_Sp_cs[idecosys][imixsys] = (TH2D*)fpisigma[0][imixsys]->Get(Form("Cosn_IMnpipi_Sp_cs_sys%d",idecosys-1));
      Cosn_IMnpipi_Sm_cs[idecosys][imixsys] = (TH2D*)fpisigma[0][imixsys]->Get(Form("Cosn_IMnpipi_Sm_cs_sys%d",idecosys-1));
      Cosn_IMnpipi_SpSmSum[idecosys][imixsys] = (TH2D*)fpisigma[0][imixsys]->Get(Form("Cosn_IMnpipi_SpSmSum_sys%d",idecosys-1));
      CosL1405[idecosys][imixsys] = (TH1D*)fpisigma[0][imixsys]->Get(Form("CosL1405_sys%d",idecosys-1));
    }
  }
    

  TCanvas *cCosIMpimSp = new TCanvas("cCosIMpimSp","cCosIMpimSp",1200,800);
  TH2D* Cosn_IMnpipi_Sp_cs_disp = (TH2D*)Cosn_IMnpipi_Sp_cs[1][1]->Clone("Cosn_IMnpipi_Sp_cs_disp");
  Cosn_IMnpipi_Sp_cs_disp->SetMinimum(0);
  Cosn_IMnpipi_Sp_cs_disp->RebinY(2);
  Cosn_IMnpipi_Sp_cs_disp->GetXaxis()->SetRangeUser(1.3,1.6);
  Cosn_IMnpipi_Sp_cs_disp->GetYaxis()->SetRangeUser(0.6,1);
  //gPad->SetLogz();
  Cosn_IMnpipi_Sp_cs_disp->Draw("colz");

  TCanvas *cCosIMpipSm = new TCanvas("cCosIMpipSm","cCosIMpipSm",1200,800);
  TH2D* Cosn_IMnpipi_Sm_cs_disp = (TH2D*)Cosn_IMnpipi_Sm_cs[1][1]->Clone("Cosn_IMnpipi_Sm_cs_disp");
  Cosn_IMnpipi_Sm_cs_disp->SetMinimum(0);
  Cosn_IMnpipi_Sm_cs_disp->RebinY(2);
  Cosn_IMnpipi_Sm_cs_disp->GetXaxis()->SetRangeUser(1.3,1.6);
  Cosn_IMnpipi_Sm_cs_disp->GetYaxis()->SetRangeUser(0.6,1);
  //gPad->SetLogz();
  Cosn_IMnpipi_Sm_cs_disp->Draw("colz");
  
  TGraphAsymmErrors *grCosL1405 = new TGraphAsymmErrors(CosL1405[1][1]);
  TGraphAsymmErrors *grCosL1405_decosysup = new TGraphAsymmErrors(CosL1405[2][1]);
  TGraphAsymmErrors *grCosL1405_decosysdown = new TGraphAsymmErrors(CosL1405[0][1]);
  for(int ip=0;ip<grCosL1405->GetN();ip++){
    double valup = grCosL1405_decosysup->GetPointY(ip);
    double valdef = grCosL1405->GetPointY(ip);
    double yeh = grCosL1405->GetErrorYhigh(ip);
    double yel = grCosL1405->GetErrorYlow(ip);
    double valdown = grCosL1405_decosysdown->GetPointY(ip);
    double diffup = valup - valdef;
    double diffdown = valdown - valdef;
     
    if(diffup > diffdown){
      yeh = sqrt(yeh*yeh + diffup*diffup);
      yel = sqrt(yel*yel + diffdown*diffdown);
    }else{
      yeh = sqrt(yeh*yeh + diffdown*diffdown);
      yel = sqrt(yel*yel + diffup*diffup);
    }
    grCosL1405->SetPointEYhigh(ip,yeh);
    grCosL1405->SetPointEYlow(ip,yel);
  }
  TGraphAsymmErrors *gMIXErrorCosL1405 = new TGraphAsymmErrors(CosL1405[1][1]);
  TGraphAsymmErrors *gMIXErrorCosL1405_sysup = new TGraphAsymmErrors(CosL1405[1][2]);
  TGraphAsymmErrors *gMIXErrorCosL1405_sysdown = new TGraphAsymmErrors(CosL1405[1][0]);
  for(int ip=0;ip<gMIXErrorCosL1405->GetN();ip++){
    double valup = gMIXErrorCosL1405_sysup->GetPointY(ip);
    double valdef = gMIXErrorCosL1405->GetPointY(ip);
    double yeh = gMIXErrorCosL1405->GetErrorYhigh(ip);
    double yel = gMIXErrorCosL1405->GetErrorYlow(ip);
    double valdown = gMIXErrorCosL1405_sysdown->GetPointY(ip);
    double diffup = valup - valdef;
    double diffdown = valdown - valdef;
     
    if(diffup > diffdown){
      yeh = fabs(diffup);
      yel = fabs(diffdown);
    }else{
      yeh = fabs(diffdown);
      yel = fabs(diffup);
    }
    gMIXErrorCosL1405->SetPointEYhigh(ip,yeh);
    gMIXErrorCosL1405->SetPointEYlow(ip,yel);
    
    double xe = CosL1405[1][1]->GetBinWidth(ip+1);
    gMIXErrorCosL1405->SetPointEXhigh(ip,xe/4);
    gMIXErrorCosL1405->SetPointEXlow(ip,xe/4);
  }
  
  TCanvas *cCosL1405 = new TCanvas("cCosL1405","cCosL1405",1200,800);
  grCosL1405->GetXaxis()->SetRangeUser(0.6,1);
  grCosL1405->GetXaxis()->SetTitle("cos#theta_{n} (CM)");
  grCosL1405->GetXaxis()->CenterTitle();
  grCosL1405->GetYaxis()->SetTitle("d#sigma / dcos#theta_{n} [#mub]  ");
  grCosL1405->GetYaxis()->CenterTitle();
  grCosL1405->Draw("ap");
  gMIXErrorCosL1405->SetLineColor(12);
  gMIXErrorCosL1405->SetFillStyle(3001);
  gMIXErrorCosL1405->SetFillColor(12);
  gMIXErrorCosL1405->Draw("5");
  TLine *pL1405 = new TLine(0.6,0,1,0);
  pL1405->SetLineColor(1);
  //p->SetLineWidth(2.0);
  pL1405->SetLineStyle(2);
  pL1405->Draw();

  TFile *f = TFile::Open("yamagataL1405.root");
  TGraph *gry = (TGraph*)f->Get("gr_yamagata");
  gry->SetLineColor(3);
  gry->SetLineWidth(3);
  //gry->Draw("c");


  TCanvas *c = NULL;
  TSeqCollection *SCol = gROOT->GetListOfCanvases();
  int size = SCol->GetSize();
  TIter next(SCol);
  TString pdfname = Form("csfinal.pdf");
  for(int i=0;i<size;i++){
    //pdf->NewPage();
    c= (TCanvas*)next();
    c->Modified();
    c->Update();
    c->Draw();
    c->cd();
    //inside the canvas
    //TPaveText *pt = new TPaveText(.74,.81,0.9,0.90,"NDC");
    c->Modified();
    c->Update();
    //std::cout << c->GetName() << std::endl;
    //make 1 pdf file
    if(i==0) c->Print(pdfname+"(",Form("pdf Title:%s",c->GetTitle()));
    else if(i==size-1)c->Print(pdfname+")",Form("pdf Title:%s",c->GetTitle())); 
    else c->Print(pdfname,Form("pdf Title:%s",c->GetTitle())); 
    //make separated pdf files
    //c->Print(Form("pdf/%s.pdf",c->GetTitle()));
  }

  TFile *fout = new TFile("csfinal.root","RECREATE");
  fout->cd();
  for(int iq=0;iq<3;iq++){
    gS1385ErrorSp[iq]->Write();
    gS1385ErrorSm[iq]->Write();
    gS1385ErrorSpSm[iq]->Write();
    gCS_coscutSp[iq]->Write();
    gCS_coscutSm[iq]->Write();
    gCS_coscutSpSm[iq]->Write();
  }

  





  return;
  
  /*
  TCanvas *c1 = new TCanvas("c1","c1",1000,800);
  c1->cd();
  IMnpipi_Sp_cs[1]->SetMarkerStyle(20);
  Mnpipi_Sp_cs[1]->SetMarkerColor(3);
  IMnpipi_Sp_cs[1]->SetXTitle("IM(#pi^{-}#Sigma^{+}) [GeV/c^{2}]");
  IMnpipi_Sm_cs[1]->SetTitle("#pi^{-}#Sigma^{+} mode");
  IMnpipi_Sp_cs[1]->Draw("E");
  TH1D* CS_IMppipi_p_wL_sum_350_ToSp = (TH1D*)CS_IMppipi_p_wL_sum_350->Clone("CS_IMppipi_p_wL_sum_350_ToSp");
  CS_IMppipi_p_wL_sum_350_ToSp->Scale(br_s1385TopiSigma/2.0/br_s1385ToLambdapi); 
  CS_IMppipi_p_wL_sum_350_ToSp->SetMarkerStyle(20); 
  CS_IMppipi_p_wL_sum_350_ToSp->SetMarkerColor(6); 
  CS_IMppipi_p_wL_sum_350_ToSp->SetMarkerColor(6); 
  CS_IMppipi_p_wL_sum_350_ToSp->Draw("Esame"); 


  TCanvas *c2 = new TCanvas("c2","c2",1000,800);
  c2->cd();
  IMnpipi_Sm_cs[1]->SetMarkerStyle(20);
  IMnpipi_Sm_cs[1]->SetMarkerColor(4);
  IMnpipi_Sm_cs[1]->SetXTitle("IM(#pi^{+}#Sigma^{-}) [GeV/c^{2}]");
  IMnpipi_Sm_cs[1]->SetTitle("#pi^{+}#Sigma^{-} mode");
  IMnpipi_Sm_cs[1]->Draw("E");
  TH1D* CS_IMppipi_p_wL_sum_350_ToSm = (TH1D*)CS_IMppipi_p_wL_sum_350->Clone("CS_IMppipi_p_wL_sum_350_ToSm");
  CS_IMppipi_p_wL_sum_350_ToSm->Scale(br_s1385TopiSigma/2.0/br_s1385ToLambdapi); 
  CS_IMppipi_p_wL_sum_350_ToSm->SetMarkerStyle(20); 
  CS_IMppipi_p_wL_sum_350_ToSm->SetMarkerColor(6); 
  CS_IMppipi_p_wL_sum_350_ToSm->Draw("Esame"); 


  TCanvas *c3 = new TCanvas("c3","c3",1000,800);
  c3->cd();
  TH1D* IMnpipi_sum_cs[2];
  IMnpipi_sum_cs[1] = (TH1D*)IMnpipi_Sp_cs[1]->Clone("IMnpipi_sum_cs1");
  IMnpipi_sum_cs[1]->SetTitle("Charge Sum");
  IMnpipi_sum_cs[1]->SetXTitle("IM(#pi#Sigma) [GeV/c^{2}]");
  IMnpipi_sum_cs[1]->Add(IMnpipi_Sm_cs[1]);
  IMnpipi_sum_cs[1]->SetLineColor(1);
  IMnpipi_sum_cs[1]->SetMarkerColor(1);
  IMnpipi_sum_cs[1]->SetYTitle("d#rho/dM [#mu b (MeV/c^{2})]");
  IMnpipi_sum_cs[1]->Draw("E");
  TH1D* IMLpim_sum_350 = (TH1D*)CS_IMppipi_p_wL_sum_350_ToSp->Clone("CS_IMppipi_p_wL_sum");
  IMLpim_sum_350->Add(CS_IMppipi_p_wL_sum_350_ToSm);
  //IMLpim_sum_350->GetXaxis()->SetRangeUser(1.35,1.5);
  IMLpim_sum_350->Scale(10.);
  IMLpim_sum_350->Draw("Esame");
  */
}


