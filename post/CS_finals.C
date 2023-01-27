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
  gStyle->SetNdivisions(505,"x");
  gStyle->SetNdivisions(505,"y");
  gStyle->SetTitleXSize(0.05);
  gStyle->SetTitleYSize(0.05);
  //gROOT->ForceStyle();
  //deco. err. is already added to statstical err in quadruture
  TH2D* q_IMnpipi_Sp_cs[4][3][3];//iq,dEcut,sysud of mix
  TH2D* q_IMnpipi_Sm_cs[4][3][3];//iq,dEcut,sysud of mix
  TH2D* q_IMnpipi_K0_cs[4][3][3];//iq,dEcut,sysud of mix
  TH2D* q_IMnpipi_SpSmSum[4][3][3];//iq,dEcut,sysud of mix
  TH1D* IMnpipi_Sp_cs[4][3][3];
  TH1D* IMnpipi_Sm_cs[4][3][3];
  TH1D* IMnpipi_K0_cs[4][3][3];
  TH1D* IMnpipi_SpSmSum[4][3][3];
  TH1D* IMnpipi_SpSmAvg[4][3][3];
  TGraphAsymmErrors* gIMnpipi_Sp_cs_Etotal[4];//CS with bin by bin total error
  TGraphAsymmErrors* gIMnpipi_Sm_cs_Etotal[4];//CS with bin by bin total error
  TGraphAsymmErrors* gIMnpipi_K0_cs_Etotal[4];//CS with bin by bin total error
  TGraphAsymmErrors* gIMnpipi_SpSmSum_cs_Etotal[4];//CS with bin by bin total error
  TGraphAsymmErrors  *gDecoErrorSp_CS[4];
  TGraphAsymmErrors  *gDecoErrorSm_CS[4];
  TGraphAsymmErrors  *gDecoErrorK0_CS[4];  
  TGraphAsymmErrors  *gDecoErrorSpSm_CS[4];
  
  for(int iq=0;iq<4;iq++){ 
    //mix center val
    for(int iEcut=0;iEcut<3;iEcut++){
      q_IMnpipi_Sp_cs[iq][iEcut][1] = (TH2D*)fpisigma[iEcut][1]->Get(Form("q_IMnpipi_Sp_cs%d_sys0",iq));
      q_IMnpipi_Sm_cs[iq][iEcut][1] = (TH2D*)fpisigma[iEcut][1]->Get(Form("q_IMnpipi_Sm_cs%d_sys0",iq));
      q_IMnpipi_K0_cs[iq][iEcut][1] = (TH2D*)fpisigma[iEcut][1]->Get(Form("q_IMnpipi_K0_cs%d_sys0",iq));
      q_IMnpipi_SpSmSum[iq][iEcut][1] = (TH2D*)fpisigma[iEcut][1]->Get(Form("q_IMnpipi_SpSmSum%d_sys0",iq));
      IMnpipi_Sp_cs[iq][iEcut][1] = (TH1D*)fpisigma[iEcut][1]->Get(Form("IMnpipi_Sp_cs%d_sys0",iq));
      IMnpipi_Sm_cs[iq][iEcut][1] = (TH1D*)fpisigma[iEcut][1]->Get(Form("IMnpipi_Sm_cs%d_sys0",iq));
      IMnpipi_K0_cs[iq][iEcut][1] = (TH1D*)fpisigma[iEcut][1]->Get(Form("IMnpipi_K0_cs%d_sys0",iq));
      IMnpipi_SpSmSum[iq][iEcut][1] = (TH1D*)fpisigma[iEcut][1]->Get(Form("IMnpipi_SpSmSum%d_sys0",iq));
      IMnpipi_SpSmAvg[iq][iEcut][1] = (TH1D*)IMnpipi_SpSmSum[iq][iEcut][1]->Clone(Form("IMnpipi_SpSmAvg%d_dE%d_sys0",iq,iEcut));
      IMnpipi_SpSmAvg[iq][iEcut][1]->Scale(1./2.);
    }
    //mix sys. down
    q_IMnpipi_Sp_cs[iq][0][0] = (TH2D*)fpisigma[0][0]->Get(Form("q_IMnpipi_Sp_cs%d_sys0",iq));
    q_IMnpipi_Sm_cs[iq][0][0] = (TH2D*)fpisigma[0][0]->Get(Form("q_IMnpipi_Sm_cs%d_sys0",iq));
    q_IMnpipi_K0_cs[iq][0][0] = (TH2D*)fpisigma[0][0]->Get(Form("q_IMnpipi_K0_cs%d_sys0",iq));
    q_IMnpipi_SpSmSum[iq][0][0] = (TH2D*)fpisigma[0][0]->Get(Form("q_IMnpipi_SpSmSum%d_sys0",iq));
    IMnpipi_Sp_cs[iq][0][0] = (TH1D*)fpisigma[0][0]->Get(Form("IMnpipi_Sp_cs%d_sys0",iq));
    IMnpipi_Sm_cs[iq][0][0] = (TH1D*)fpisigma[0][0]->Get(Form("IMnpipi_Sm_cs%d_sys0",iq));
    IMnpipi_K0_cs[iq][0][0] = (TH1D*)fpisigma[0][0]->Get(Form("IMnpipi_K0_cs%d_sys0",iq));
    IMnpipi_SpSmSum[iq][0][0] = (TH1D*)fpisigma[0][0]->Get(Form("IMnpipi_SpSmSum%d_sys0",iq));
    IMnpipi_SpSmAvg[iq][0][0] = (TH1D*)IMnpipi_SpSmSum[iq][0][0]->Clone(Form("IMnpipi_SpSmAvg%d_sys-1",iq));
    IMnpipi_SpSmAvg[iq][0][0]->Scale(1./2.);
    //mix sys. up
    q_IMnpipi_Sp_cs[iq][0][2] = (TH2D*)fpisigma[0][2]->Get(Form("q_IMnpipi_Sp_cs%d_sys0",iq));
    q_IMnpipi_Sm_cs[iq][0][2] = (TH2D*)fpisigma[0][2]->Get(Form("q_IMnpipi_Sm_cs%d_sys0",iq));
    q_IMnpipi_K0_cs[iq][0][2] = (TH2D*)fpisigma[0][2]->Get(Form("q_IMnpipi_K0_cs%d_sys0",iq));
    q_IMnpipi_SpSmSum[iq][0][2] = (TH2D*)fpisigma[0][2]->Get(Form("q_IMnpipi_SpSmSum%d_sys0",iq));
    IMnpipi_Sp_cs[iq][0][2] = (TH1D*)fpisigma[0][2]->Get(Form("IMnpipi_Sp_cs%d_sys0",iq));
    IMnpipi_Sm_cs[iq][0][2] = (TH1D*)fpisigma[0][2]->Get(Form("IMnpipi_Sm_cs%d_sys0",iq));
    IMnpipi_K0_cs[iq][0][2] = (TH1D*)fpisigma[0][2]->Get(Form("IMnpipi_K0_cs%d_sys0",iq));
    IMnpipi_SpSmSum[iq][0][2] = (TH1D*)fpisigma[0][2]->Get(Form("IMnpipi_SpSmSum%d_sys0",iq));
    IMnpipi_SpSmAvg[iq][0][2] = (TH1D*)IMnpipi_SpSmSum[iq][0][2]->Clone(Form("IMnpipi_SpSmAvg%d_sys1",iq));
    IMnpipi_SpSmAvg[iq][0][2]->Scale(1./2.);
    
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
  TGraphAsymmErrors *gMIXErrorSpSmAvg_CS[4];
  TGraphAsymmErrors *gdEErrorSp_CS[4];
  TGraphAsymmErrors *gdEErrorSm_CS[4];
  TGraphAsymmErrors *gdEErrorK0_CS[4];

  
  for(int iq=0;iq<4;iq++){
    gMIXErrorSp_CS[iq] = new TGraphAsymmErrors(IMnpipi_Sp_cs[iq][0][1]);
    gMIXErrorSm_CS[iq] = new TGraphAsymmErrors(IMnpipi_Sm_cs[iq][0][1]);
    gMIXErrorK0_CS[iq] = new TGraphAsymmErrors(IMnpipi_K0_cs[iq][0][1]);
    gMIXErrorSpSm_CS[iq] = new TGraphAsymmErrors(IMnpipi_SpSmSum[iq][0][1]);
    gMIXErrorSpSmAvg_CS[iq] = new TGraphAsymmErrors(IMnpipi_SpSmAvg[iq][0][1]);
    gMIXErrorSp_CS[iq]->SetName(Form("gMIXErrorSp_CS%d",iq));
    gMIXErrorSm_CS[iq]->SetName(Form("gMIXErrorSm_CS%d",iq));
    gMIXErrorK0_CS[iq]->SetName(Form("gMIXErrorK0_CS%d",iq));
    gMIXErrorSpSm_CS[iq]->SetName(Form("gMIXErrorSpSm_CS%d",iq));
    gMIXErrorSpSmAvg_CS[iq]->SetName(Form("gMIXErrorSpSmAvg_CS%d",iq));

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
    
    for(int ip=0;ip<( gMIXErrorSpSmAvg_CS[iq]->GetN());ip++){
      double xe = IMnpipi_SpSmSum[iq][0][0]->GetBinWidth(ip+1);
      gMIXErrorSpSmAvg_CS[iq]->SetPointEXhigh(ip,xe/4);
      gMIXErrorSpSmAvg_CS[iq]->SetPointEXlow(ip,xe/4);
      
      double  valdown = IMnpipi_SpSmAvg[iq][0][0]->GetBinContent(ip+1);
      double  valdef = IMnpipi_SpSmAvg[iq][0][1]->GetBinContent(ip+1);
      double  valup  = IMnpipi_SpSmAvg[iq][0][2]->GetBinContent(ip+1);
      
      double yeh = valup-valdef;
      double yel = valdown-valdef;
         
      gMIXErrorSpSmAvg_CS[iq]->SetPointEYhigh(ip,fabs(yeh));
      gMIXErrorSpSmAvg_CS[iq]->SetPointEYlow(ip,fabs(yel));
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
  
  TFile *flpimnofit = TFile::Open("cs_lpim_killcombi.root","READ");
  TH1D* CS_IMppipi_p_wL_sum_0 = (TH1D*)flpimnofit->Get("CS_IMppipi_p_wL_sum_0");
  TH1D* CS_IMppipi_p_wL_sum_300 = (TH1D*)flpimnofit->Get("CS_IMppipi_p_wL_sum_300");
  TCanvas *cSysSpSmAvgEtotal[4];
  TGraphAsymmErrors *gIMnpipi_SpSmAvg_cs_Etotal[4];
  for(int iq=1;iq<3;iq++){  
    cSysSpSmAvgEtotal[iq] = new TCanvas(Form("cSysSpSmAvgEtotal%d",iq),Form("cSysSpSmAvgEtotal%d",iq),1100,1100);
    
    cSysSpSmAvgEtotal[iq]->SetBottomMargin(0.12);
    cSysSpSmAvgEtotal[iq]->SetLeftMargin(0.14);
    cSysSpSmAvgEtotal[iq]->SetRightMargin(0.08);
    
    gIMnpipi_SpSmAvg_cs_Etotal[iq] = new TGraphAsymmErrors(IMnpipi_SpSmAvg[iq][0][1]);
    gIMnpipi_SpSmAvg_cs_Etotal[iq]->SetName(Form("gIMnpipi_SpSmAvg_cs_Etotal%d",iq));
    gIMnpipi_SpSmAvg_cs_Etotal[iq]->SetMarkerColor(1);
    gIMnpipi_SpSmAvg_cs_Etotal[iq]->GetXaxis()->SetRangeUser(1.3,1.6);
    gIMnpipi_SpSmAvg_cs_Etotal[iq]->SetTitle("");
    gIMnpipi_SpSmAvg_cs_Etotal[iq]->GetXaxis()->SetTitle("IM(#pi#Sigma) [GeV/c^{2}]");
    gIMnpipi_SpSmAvg_cs_Etotal[iq]->GetXaxis()->CenterTitle();
    gIMnpipi_SpSmAvg_cs_Etotal[iq]->GetYaxis()->SetTitle("d#sigma/dM  [#mub/(MeV/c^{2})]");
    gIMnpipi_SpSmAvg_cs_Etotal[iq]->GetYaxis()->CenterTitle();

    gIMnpipi_SpSmAvg_cs_Etotal[iq]->GetXaxis()->SetTitleSize(0.04);
    gIMnpipi_SpSmAvg_cs_Etotal[iq]->GetYaxis()->SetTitleSize(0.04);
    gIMnpipi_SpSmAvg_cs_Etotal[iq]->GetXaxis()->SetTitleOffset(1.4);
    gIMnpipi_SpSmAvg_cs_Etotal[iq]->GetYaxis()->SetTitleOffset(1.6);
   if(iq==1){
      for(int ip=0;ip<10;ip++){
        gIMnpipi_SpSmAvg_cs_Etotal[iq]->RemovePoint(0);
      }
    }
    if(iq==2){
      for(int ip=0;ip<8;ip++){
        gIMnpipi_SpSmAvg_cs_Etotal[iq]->RemovePoint(0);
      }
    }
   
    double start1 = gIMnpipi_SpSmAvg_cs_Etotal[iq]->GetPointX(0);
    double start2 = gDecoErrorSpSm_CS[iq]->GetPointX(0);
    int diffstart = (start2-start1)/0.015 + 1 ;

    for(int ip=0;ip<gIMnpipi_SpSmAvg_cs_Etotal[iq]->GetN()-diffstart;ip++){
      double statElow  = gIMnpipi_SpSmAvg_cs_Etotal[iq]->GetErrorYlow(ip+diffstart);
      double statEhigh = gIMnpipi_SpSmAvg_cs_Etotal[iq]->GetErrorYhigh(ip+diffstart);
      double decoElow  = (gDecoErrorSpSm_CS[iq]->GetErrorYlow(ip))/2.0;
      double decoEhigh = (gDecoErrorSpSm_CS[iq]->GetErrorYhigh(ip))/2.0;
      double Elow = sqrt(pow(statElow,2.0)+pow(decoElow,2.0));  
      double Ehigh = sqrt(pow(statEhigh,2.0)+pow(decoEhigh,2.0));  
      gIMnpipi_SpSmAvg_cs_Etotal[iq]->SetPointEYlow(ip+diffstart,Elow);
      gIMnpipi_SpSmAvg_cs_Etotal[iq]->SetPointEYhigh(ip+diffstart,Ehigh);
    }    
    gIMnpipi_SpSmAvg_cs_Etotal[iq]->Draw("ap");
    gMIXErrorSpSmAvg_CS[iq]->SetFillStyle(3001);
    //gMIXErrorSpSmAvg_CS[iq]->SetFillColor(4);
    gMIXErrorSpSmAvg_CS[iq]->SetFillColor(12);
    gMIXErrorSpSmAvg_CS[iq]->SetMarkerColor(4);
    gMIXErrorSpSmAvg_CS[iq]->SetLineColor(12);
    gMIXErrorSpSmAvg_CS[iq]->Draw("5");
    if(iq==1){
      gr_S1385_ToSqlow->SetLineColor(5);
      gr_S1385_ToSqlow->SetFillStyle(3002);
      gr_S1385_ToSqlow->SetFillColor(0);
      //gr_S1385_ToSqlow->Draw("5");
      CS_IMppipi_p_wL_sum_0->SetLineColor(5);
      CS_IMppipi_p_wL_sum_0->SetFillColor(5);
      CS_IMppipi_p_wL_sum_0->SetMarkerColor(5);
      CS_IMppipi_p_wL_sum_0->Draw("Esame");
    }
    if(iq==2){
      gr_S1385_ToSqhi->SetLineColor(5);
      gr_S1385_ToSqhi->SetFillStyle(3002);
      gr_S1385_ToSqhi->SetFillColor(0);
      //gr_S1385_ToSqhi->Draw("5");
      CS_IMppipi_p_wL_sum_300->SetLineColor(5);
      CS_IMppipi_p_wL_sum_300->SetFillColor(5);
      CS_IMppipi_p_wL_sum_300->SetMarkerColor(5);
      CS_IMppipi_p_wL_sum_300->Draw("Esame");
    }
    TLine *p = new TLine(1.29,0,1.605,0);
    p->SetLineColor(1);
    //p->SetLineWidth(2.0);
    p->SetLineStyle(2);
    p->Draw();
  }
  cSysSpSmAvgEtotal[1]->SaveAs("Avgqlow.pdf");
  cSysSpSmAvgEtotal[2]->SaveAs("Avgqhi.pdf");
   
  


  TCanvas *cSysSpSmAvgEtotalSide = new TCanvas("cSysSpSmAvgEtotalSide","cSysSpSmAvgEtotalSide",1700,800);
  cSysSpSmAvgEtotalSide->SetBottomMargin(0.12);
  cSysSpSmAvgEtotalSide->SetLeftMargin(0.14);
  cSysSpSmAvgEtotalSide->SetRightMargin(0.08);
  cSysSpSmAvgEtotalSide->Divide(2,1,0,0);
  cSysSpSmAvgEtotalSide->cd(1);
  gIMnpipi_SpSmAvg_cs_Etotal[1]->GetYaxis()->SetRangeUser(-0.2,1.6);
  gIMnpipi_SpSmAvg_cs_Etotal[2]->GetYaxis()->SetRangeUser(-0.2,1.6);
  gIMnpipi_SpSmAvg_cs_Etotal[1]->Draw("ap");
  gMIXErrorSpSmAvg_CS[1]->Draw("5");
  CS_IMppipi_p_wL_sum_0->SetLineColor(5);
  CS_IMppipi_p_wL_sum_0->SetFillColor(5);
  CS_IMppipi_p_wL_sum_0->SetMarkerColor(5);
  CS_IMppipi_p_wL_sum_0->SetMarkerStyle(20);
  CS_IMppipi_p_wL_sum_0->SetMarkerSize(1.2);
  CS_IMppipi_p_wL_sum_0->Draw("Esame");
  //gr_S1385_ToSqlow->Draw("5");
  TLine *ps = new TLine(1.29,0,1.605,0);
  ps->SetLineColor(1);
  //p->SetLineWidth(2.0);
  ps->SetLineStyle(2);
  ps->Draw();
  TLine *pkp = new TLine(kpMass+pMass,0,kpMass+pMass,3);
  pkp->SetLineColor(1);
  pkp->SetLineStyle(2);
  pkp->Draw();

  TLine *pk0n = new TLine(k0Mass+nMass,0,k0Mass+nMass,3);
  pk0n->SetLineColor(1);
  pk0n->SetLineStyle(2);
  pk0n->Draw();
  cSysSpSmAvgEtotalSide->cd(2);
  gIMnpipi_SpSmAvg_cs_Etotal[2]->Draw("ap");
  gMIXErrorSpSmAvg_CS[2]->Draw("5");
  //gr_S1385_ToSqhi->Draw("5");
  CS_IMppipi_p_wL_sum_300->SetLineColor(5);
  CS_IMppipi_p_wL_sum_300->SetLineWidth(2);
  CS_IMppipi_p_wL_sum_300->SetFillColor(5);
  CS_IMppipi_p_wL_sum_300->SetMarkerColor(5);
  CS_IMppipi_p_wL_sum_300->SetMarkerStyle(20);
  CS_IMppipi_p_wL_sum_300->SetMarkerSize(1.2);
  CS_IMppipi_p_wL_sum_300->Draw("Esame");
  ps->Draw();
  pkp->Draw();
  pk0n->Draw();
  cSysSpSmAvgEtotalSide->SaveAs("cSpSmAvgSide.pdf");
  

  TGraphAsymmErrors* gIMnpipi_SpSmAvg_cs_Etotal_mev[2];
  TGraphAsymmErrors* gMIXErrorSpSmAvg_CS_mev[2];
  
  for(int iqlh=0;iqlh<2;iqlh++){
    gIMnpipi_SpSmAvg_cs_Etotal_mev[iqlh] = new TGraphAsymmErrors();
    for(int ip=0;ip<gIMnpipi_SpSmAvg_cs_Etotal[iqlh+1]->GetN();ip++){
      double x = gIMnpipi_SpSmAvg_cs_Etotal[iqlh+1]->GetPointX(ip);
      double y = gIMnpipi_SpSmAvg_cs_Etotal[iqlh+1]->GetPointY(ip);
      double xe = gIMnpipi_SpSmAvg_cs_Etotal[iqlh+1]->GetErrorXhigh(ip);
      double yeh = gIMnpipi_SpSmAvg_cs_Etotal[iqlh+1]->GetErrorYhigh(ip);
      double yel = gIMnpipi_SpSmAvg_cs_Etotal[iqlh+1]->GetErrorYlow(ip);
      gIMnpipi_SpSmAvg_cs_Etotal_mev[iqlh]->SetPoint(ip,x*1000.,y);
      gIMnpipi_SpSmAvg_cs_Etotal_mev[iqlh]->SetPointEXhigh(ip,xe*1000.);
      gIMnpipi_SpSmAvg_cs_Etotal_mev[iqlh]->SetPointEXlow(ip,xe*1000.);
      gIMnpipi_SpSmAvg_cs_Etotal_mev[iqlh]->SetPointEYlow(ip,yel);
      gIMnpipi_SpSmAvg_cs_Etotal_mev[iqlh]->SetPointEYhigh(ip,yeh);
    }//for ip
    gMIXErrorSpSmAvg_CS_mev[iqlh] = new TGraphAsymmErrors();
    for(int ip=0;ip<gMIXErrorSpSmAvg_CS[iqlh+1]->GetN();ip++){
      double x = gMIXErrorSpSmAvg_CS[iqlh+1]->GetPointX(ip);
      double y = gMIXErrorSpSmAvg_CS[iqlh+1]->GetPointY(ip);
      double xe = gMIXErrorSpSmAvg_CS[iqlh+1]->GetErrorXhigh(ip);
      double yeh = gMIXErrorSpSmAvg_CS[iqlh+1]->GetErrorYhigh(ip);
      double yel = gMIXErrorSpSmAvg_CS[iqlh+1]->GetErrorYlow(ip);
      gMIXErrorSpSmAvg_CS_mev[iqlh]->SetPoint(ip,x*1000.,y);
      gMIXErrorSpSmAvg_CS_mev[iqlh]->SetPoint(ip,x*1000.,y);
      gMIXErrorSpSmAvg_CS_mev[iqlh]->SetPointEXhigh(ip,xe*1000.);
      gMIXErrorSpSmAvg_CS_mev[iqlh]->SetPointEXlow(ip,xe*1000.);
      gMIXErrorSpSmAvg_CS_mev[iqlh]->SetPointEYlow(ip,yel);
      gMIXErrorSpSmAvg_CS_mev[iqlh]->SetPointEYhigh(ip,yeh);
    }
  }




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


  TCanvas *cboth_mev[2];
  TGraphAsymmErrors* gIMnpipi_Sp_cs_Etotal_mev[2];    
  TGraphAsymmErrors* gIMnpipi_Sm_cs_Etotal_mev[2];
  TGraphAsymmErrors* gMIXErrorSp_CS_mev[2];
  TGraphAsymmErrors* gMIXErrorSm_CS_mev[2];
  TGraphAsymmErrors* gr_S1385_ToSqlow_mev;
  TGraphAsymmErrors* gr_S1385_ToSqhi_mev;
  
  for(int iqlh=0;iqlh<2;iqlh++){
    gIMnpipi_Sp_cs_Etotal_mev[iqlh] = new TGraphAsymmErrors();
    for(int ip=0;ip<gIMnpipi_Sp_cs_Etotal[iqlh+1]->GetN();ip++){
      double x = gIMnpipi_Sp_cs_Etotal[iqlh+1]->GetPointX(ip);
      double y = gIMnpipi_Sp_cs_Etotal[iqlh+1]->GetPointY(ip);
      double xe = gIMnpipi_Sp_cs_Etotal[iqlh+1]->GetErrorXhigh(ip);
      double yeh = gIMnpipi_Sp_cs_Etotal[iqlh+1]->GetErrorYhigh(ip);
      double yel = gIMnpipi_Sp_cs_Etotal[iqlh+1]->GetErrorYlow(ip);
      gIMnpipi_Sp_cs_Etotal_mev[iqlh]->SetPoint(ip,x*1000.,y);
      gIMnpipi_Sp_cs_Etotal_mev[iqlh]->SetPoint(ip,x*1000.,y);
      gIMnpipi_Sp_cs_Etotal_mev[iqlh]->SetPointEXhigh(ip,xe*1000.);
      gIMnpipi_Sp_cs_Etotal_mev[iqlh]->SetPointEXlow(ip,xe*1000.);
      gIMnpipi_Sp_cs_Etotal_mev[iqlh]->SetPointEYlow(ip,yel);
      gIMnpipi_Sp_cs_Etotal_mev[iqlh]->SetPointEYhigh(ip,yeh);
    }//for ip
    if(iqlh==0)gIMnpipi_Sp_cs_Etotal_mev[iqlh]->GetYaxis()->SetRangeUser(-0.2,2.5);
    if(iqlh==1)gIMnpipi_Sp_cs_Etotal_mev[iqlh]->GetYaxis()->SetRangeUser(-0.4,1);
    gIMnpipi_Sp_cs_Etotal_mev[iqlh]->SetLineColor(kGreen+2);
    gIMnpipi_Sp_cs_Etotal_mev[iqlh]->SetFillColor(0);
    gIMnpipi_Sp_cs_Etotal_mev[iqlh]->GetXaxis()->SetTitle("IM(#pi^{#pm}#Sigma^{#mp}) [MeV/c^{2}]");
    gIMnpipi_Sp_cs_Etotal_mev[iqlh]->GetXaxis()->CenterTitle();
    gIMnpipi_Sp_cs_Etotal_mev[iqlh]->GetYaxis()->SetTitle("d#sigma/dM  [#mub/MeV^{2}]");
    gIMnpipi_Sp_cs_Etotal_mev[iqlh]->GetYaxis()->CenterTitle();
    gIMnpipi_Sp_cs_Etotal_mev[iqlh]->GetXaxis()->SetTitleSize(0.04);
    gIMnpipi_Sp_cs_Etotal_mev[iqlh]->SetLineWidth(2);
    gIMnpipi_Sp_cs_Etotal_mev[iqlh]->GetYaxis()->SetTitleSize(0.04);
    gIMnpipi_Sp_cs_Etotal_mev[iqlh]->GetXaxis()->SetTitleOffset(1.4);
    gIMnpipi_Sp_cs_Etotal_mev[iqlh]->GetYaxis()->SetTitleOffset(1.6);
    gIMnpipi_Sp_cs_Etotal_mev[iqlh]->SetMarkerColor(kGreen+2);
    if(iqlh==1){
      gIMnpipi_Sp_cs_Etotal_mev[1]->SetMaximum(gIMnpipi_Sp_cs_Etotal_mev[0]->GetHistogram()->GetMaximum());
      gIMnpipi_Sp_cs_Etotal_mev[0]->SetMinimum(gIMnpipi_Sp_cs_Etotal_mev[1]->GetHistogram()->GetMinimum());
      //double maxY = TMath::MaxElement(n,gIMnpipi_Sp_cs_Etotal_mev[0]->GetY());
      //std::cout << "MAX" << gIMnpipi_Sp_cs_Etotal_mev[0]->GetMaximum() << std::endl;
    }
    gMIXErrorSp_CS_mev[iqlh] = new TGraphAsymmErrors();
    for(int ip=0;ip<gMIXErrorSp_CS[iqlh+1]->GetN();ip++){
      double x = gMIXErrorSp_CS[iqlh+1]->GetPointX(ip);
      double y = gMIXErrorSp_CS[iqlh+1]->GetPointY(ip);
      double xe = gMIXErrorSp_CS[iqlh+1]->GetErrorXhigh(ip);
      double yeh = gMIXErrorSp_CS[iqlh+1]->GetErrorYhigh(ip);
      double yel = gMIXErrorSp_CS[iqlh+1]->GetErrorYlow(ip);
      gMIXErrorSp_CS_mev[iqlh]->SetPoint(ip,x*1000.,y);
      gMIXErrorSp_CS_mev[iqlh]->SetPoint(ip,x*1000.,y);
      gMIXErrorSp_CS_mev[iqlh]->SetPointEXhigh(ip,xe*1000.);
      gMIXErrorSp_CS_mev[iqlh]->SetPointEXlow(ip,xe*1000.);
      gMIXErrorSp_CS_mev[iqlh]->SetPointEYlow(ip,yel);
      gMIXErrorSp_CS_mev[iqlh]->SetPointEYhigh(ip,yeh);
    }
    gMIXErrorSp_CS_mev[iqlh]->SetLineColor(kGreen+2);
    gMIXErrorSp_CS_mev[iqlh]->SetLineWidth(2);
    gMIXErrorSp_CS_mev[iqlh]->SetFillColor(0);
    gMIXErrorSp_CS_mev[iqlh]->SetFillStyle(0);
    gMIXErrorSp_CS_mev[iqlh]->SetMarkerColor(kGreen+2);
    gIMnpipi_Sm_cs_Etotal_mev[iqlh] = new TGraphAsymmErrors();
    for(int ip=0;ip<gIMnpipi_Sm_cs_Etotal[iqlh+1]->GetN();ip++){
      double x = gIMnpipi_Sm_cs_Etotal[iqlh+1]->GetPointX(ip);
      double y = gIMnpipi_Sm_cs_Etotal[iqlh+1]->GetPointY(ip);
      double xe = gIMnpipi_Sm_cs_Etotal[iqlh+1]->GetErrorXhigh(ip);
      double yeh = gIMnpipi_Sm_cs_Etotal[iqlh+1]->GetErrorYhigh(ip);
      double yel = gIMnpipi_Sm_cs_Etotal[iqlh+1]->GetErrorYlow(ip);
      gIMnpipi_Sm_cs_Etotal_mev[iqlh]->SetPoint(ip,x*1000.,y);
      gIMnpipi_Sm_cs_Etotal_mev[iqlh]->SetPoint(ip,x*1000.,y);
      gIMnpipi_Sm_cs_Etotal_mev[iqlh]->SetPointEXhigh(ip,xe*1000.);
      gIMnpipi_Sm_cs_Etotal_mev[iqlh]->SetPointEXlow(ip,xe*1000.);
      gIMnpipi_Sm_cs_Etotal_mev[iqlh]->SetPointEYlow(ip,yel);
      gIMnpipi_Sm_cs_Etotal_mev[iqlh]->SetPointEYhigh(ip,yeh);
    }
    gIMnpipi_Sm_cs_Etotal_mev[iqlh]->SetLineColor(kBlue+2);
    gIMnpipi_Sm_cs_Etotal_mev[iqlh]->SetLineWidth(2);
    gIMnpipi_Sm_cs_Etotal_mev[iqlh]->SetFillColor(0);
    gIMnpipi_Sm_cs_Etotal_mev[iqlh]->SetMarkerColor(kBlue+2);
    gMIXErrorSm_CS_mev[iqlh] = new TGraphAsymmErrors();
    for(int ip=0;ip<gMIXErrorSm_CS[iqlh+1]->GetN();ip++){
      double x = gMIXErrorSm_CS[iqlh+1]->GetPointX(ip);
      double y = gMIXErrorSm_CS[iqlh+1]->GetPointY(ip);
      double xe = gMIXErrorSm_CS[iqlh+1]->GetErrorXhigh(ip);
      double yeh = gMIXErrorSm_CS[iqlh+1]->GetErrorYhigh(ip);
      double yel = gMIXErrorSm_CS[iqlh+1]->GetErrorYlow(ip);
      gMIXErrorSm_CS_mev[iqlh]->SetPoint(ip,x*1000.,y);
      gMIXErrorSm_CS_mev[iqlh]->SetPoint(ip,x*1000.,y);
      gMIXErrorSm_CS_mev[iqlh]->SetPointEXhigh(ip,xe*1000.);
      gMIXErrorSm_CS_mev[iqlh]->SetPointEXlow(ip,xe*1000.);
      gMIXErrorSm_CS_mev[iqlh]->SetPointEYlow(ip,yel);
      gMIXErrorSm_CS_mev[iqlh]->SetPointEYhigh(ip,yeh);
    }
    gMIXErrorSm_CS_mev[iqlh]->SetLineColor(kBlue+2);
    gMIXErrorSm_CS_mev[iqlh]->SetLineWidth(2);
    gMIXErrorSm_CS_mev[iqlh]->SetFillColor(0);
    gMIXErrorSm_CS_mev[iqlh]->SetFillStyle(0);
    gMIXErrorSm_CS_mev[iqlh]->SetMarkerColor(kBlue+2);
    if(iqlh==0){
      gr_S1385_ToSqlow_mev = new TGraphAsymmErrors();
      for(int ip=0;ip<gr_S1385_ToSqlow->GetN();ip++){
        double x = gr_S1385_ToSqlow->GetPointX(ip);
        double y = gr_S1385_ToSqlow->GetPointY(ip);
        double xe = gr_S1385_ToSqlow->GetErrorXhigh(ip);
        double yeh = gr_S1385_ToSqlow->GetErrorYhigh(ip);
        double yel = gr_S1385_ToSqlow->GetErrorYlow(ip);
        gr_S1385_ToSqlow_mev->SetPoint(ip,x*1000.,y);
        gr_S1385_ToSqlow_mev->SetPoint(ip,x*1000.,y);
        gr_S1385_ToSqlow_mev->SetPointEXhigh(ip,xe*1000.);
        gr_S1385_ToSqlow_mev->SetPointEXlow(ip,xe*1000.);
        gr_S1385_ToSqlow_mev->SetPointEYlow(ip,yel);
        gr_S1385_ToSqlow_mev->SetPointEYhigh(ip,yeh);
      }
      gr_S1385_ToSqlow_mev->SetLineColor(5);
      gr_S1385_ToSqlow_mev->SetFillStyle(3002);
      gr_S1385_ToSqlow_mev->SetFillColor(0);
    }else if(iqlh==1){
      gr_S1385_ToSqhi_mev = new TGraphAsymmErrors();
      for(int ip=0;ip<gr_S1385_ToSqhi->GetN();ip++){
        double x = gr_S1385_ToSqhi->GetPointX(ip);
        double y = gr_S1385_ToSqhi->GetPointY(ip);
        double xe = gr_S1385_ToSqhi->GetErrorXhigh(ip);
        double yeh = gr_S1385_ToSqhi->GetErrorYhigh(ip);
        double yel = gr_S1385_ToSqhi->GetErrorYlow(ip);
        gr_S1385_ToSqhi_mev->SetPoint(ip,x*1000.,y);
        gr_S1385_ToSqhi_mev->SetPoint(ip,x*1000.,y);
        gr_S1385_ToSqhi_mev->SetPointEXhigh(ip,xe*1000.);
        gr_S1385_ToSqhi_mev->SetPointEXlow(ip,xe*1000.);
        gr_S1385_ToSqhi_mev->SetPointEYlow(ip,yel);
        gr_S1385_ToSqhi_mev->SetPointEYhigh(ip,yeh);
      }
      gr_S1385_ToSqhi_mev->SetLineColor(5);
      gr_S1385_ToSqhi_mev->SetFillStyle(3002);
      gr_S1385_ToSqhi_mev->SetFillColor(0);
    }
  }

  for(int iqlh=0;iqlh<2;iqlh++){
    cboth_mev[iqlh] = new TCanvas(Form("cboth_mev%d",iqlh),Form("cboth_mev%d",iqlh),1100,1100);
    cboth_mev[iqlh]->SetBottomMargin(0.12);
    cboth_mev[iqlh]->SetLeftMargin(0.14);
    cboth_mev[iqlh]->SetRightMargin(0.08);
    cboth_mev[iqlh]->cd(); 
    gIMnpipi_Sp_cs_Etotal_mev[iqlh]->GetXaxis()->SetLimits(1300,1610);
    gIMnpipi_Sp_cs_Etotal_mev[iqlh]->Draw("ap");
    gMIXErrorSp_CS_mev[iqlh]->Draw("5");
    gIMnpipi_Sm_cs_Etotal_mev[iqlh]->Draw("p");
    gMIXErrorSm_CS_mev[iqlh]->Draw("5");
     
    if(iqlh==0){
      //gr_S1385_ToSqlow_mev->Draw("5");
      TLine *p = new TLine(1290,0,1605,0);
      p->SetLineColor(1);
      //p->SetLineWidth(2.0);
      p->SetLineStyle(2);
      p->Draw();
      TLine *pkp = new TLine((kpMass+pMass)*1000,0,(kpMass+pMass)*1000,gIMnpipi_Sp_cs_Etotal_mev[0]->GetHistogram()->GetMaximum());
      pkp->SetLineColor(1);
      pkp->SetLineStyle(2);
      pkp->Draw();

      TLine *pk0n = new TLine((k0Mass+nMass)*1000,0,(k0Mass+nMass)*1000,gIMnpipi_Sp_cs_Etotal_mev[0]->GetHistogram()->GetMaximum());
      pk0n->SetLineColor(1);
      pk0n->SetLineStyle(2);
      pk0n->Draw(); 
      TLatex *tex = new TLatex();
      double tex_ymax = gIMnpipi_Sp_cs_Etotal_mev[0]->GetHistogram()->GetMaximum();
      tex->SetTextSize(0.05);
      tex->SetTextColor(1);
      tex->DrawLatex( 1320,tex_ymax*0.85 , "(a)" );
    }else if(iqlh==1){
      //gr_S1385_ToSqhi_mev->Draw("5");
      TLine *p = new TLine(1290,0,1605,0);
      p->SetLineColor(1);
      //p->SetLineWidth(2.0);
      p->SetLineStyle(2);
      p->Draw();
      TLine *pkp = new TLine((kpMass+pMass)*1000,0,(kpMass+pMass)*1000,gIMnpipi_Sp_cs_Etotal_mev[0]->GetHistogram()->GetMaximum());
      pkp->SetLineColor(1);
      pkp->SetLineStyle(2);
      pkp->Draw();

      TLine *pk0n = new TLine((k0Mass+nMass)*1000,0,(k0Mass+nMass)*1000,gIMnpipi_Sp_cs_Etotal_mev[0]->GetHistogram()->GetMaximum());
      pk0n->SetLineColor(1);
      pk0n->SetLineStyle(2);
      pk0n->Draw(); 

      TLatex *tex = new TLatex();
      double tex_ymax = gIMnpipi_Sp_cs_Etotal_mev[1]->GetHistogram()->GetMaximum();
      tex->SetTextSize(0.05);
      tex->SetTextColor(1);
      tex->DrawLatex( 1320,tex_ymax*0.85 , "(b)" );
    }

    cboth_mev[iqlh]->Update(); 
    cboth_mev[iqlh]->Modified(); 

   
  }
  cboth_mev[0]->SaveAs("cbothqlow_mev.pdf");
  cboth_mev[1]->SaveAs("cbothqhi_mev.pdf");

  TCanvas *cboth_mev_side = new TCanvas("cboth_mev_side","cboth_mev_side",2200,1050);
  cboth_mev_side->SetBottomMargin(0.23);
  cboth_mev_side->SetTopMargin(-0.50);
  cboth_mev_side->SetLeftMargin(0.21);
  cboth_mev_side->Divide(2,1,0,0);
  for(int iqlh=0;iqlh<2;iqlh++){
    cboth_mev_side->cd(iqlh+1);
    if(iqlh==0)gIMnpipi_Sp_cs_Etotal_mev[iqlh]->GetXaxis()->SetLimits(1300,1620);
    if(iqlh==1)gIMnpipi_Sp_cs_Etotal_mev[iqlh]->GetXaxis()->SetLimits(1301,1620);
    gIMnpipi_Sp_cs_Etotal_mev[iqlh]->GetXaxis()->SetLabelSize(0.05);
    gIMnpipi_Sp_cs_Etotal_mev[iqlh]->GetYaxis()->SetLabelSize(0.05);
    gIMnpipi_Sp_cs_Etotal_mev[iqlh]->GetXaxis()->SetLabelOffset(0.05);
    gIMnpipi_Sp_cs_Etotal_mev[iqlh]->GetYaxis()->SetLabelOffset(0.05);
    gIMnpipi_Sp_cs_Etotal_mev[iqlh]->GetXaxis()->SetTickLength(0.02);
    gIMnpipi_Sp_cs_Etotal_mev[iqlh]->GetYaxis()->SetTickLength(0.02);
    gIMnpipi_Sp_cs_Etotal_mev[iqlh]->GetXaxis()->SetTitleSize(0.07);
    gIMnpipi_Sp_cs_Etotal_mev[iqlh]->GetYaxis()->SetTitleSize(0.07);
    gIMnpipi_Sp_cs_Etotal_mev[iqlh]->GetXaxis()->SetTitleOffset(1.4);
    gIMnpipi_Sp_cs_Etotal_mev[iqlh]->GetYaxis()->SetTitleOffset(1.4);
   

    gIMnpipi_Sp_cs_Etotal_mev[iqlh]->SetMaximum(gIMnpipi_Sp_cs_Etotal_mev[iqlh]->GetHistogram()->GetMaximum()*1.03);
    gIMnpipi_Sp_cs_Etotal_mev[iqlh]->Draw("ap");
    gMIXErrorSp_CS_mev[iqlh]->Draw("5");
    gIMnpipi_Sm_cs_Etotal_mev[iqlh]->Draw("p");
    gMIXErrorSm_CS_mev[iqlh]->Draw("5");
    if(iqlh==0){
      //gr_S1385_ToSqlow_mev->Draw("5");
      TLine *p = new TLine(1300,0,1605,0);
      p->SetLineColor(1);
      p->SetLineStyle(2);
      p->Draw();
      TLine *pkp = new TLine((kpMass+pMass)*1000,0,(kpMass+pMass)*1000,gIMnpipi_Sp_cs_Etotal_mev[0]->GetHistogram()->GetMaximum());
      pkp->SetLineColor(1);
      pkp->SetLineStyle(2);
      pkp->Draw();

      TLine *pk0n = new TLine((k0Mass+nMass)*1000,0,(k0Mass+nMass)*1000,gIMnpipi_Sp_cs_Etotal_mev[0]->GetHistogram()->GetMaximum());
      pk0n->SetLineColor(1);
      pk0n->SetLineStyle(2);
      pk0n->Draw(); 
      TLatex *tex = new TLatex();
      double tex_ymax = gIMnpipi_Sp_cs_Etotal_mev[0]->GetHistogram()->GetMaximum();
      tex->SetTextSize(0.06);
      tex->SetTextColor(1);
      tex->DrawLatex( 1320,tex_ymax*0.85 , "(a)" );
      tex->DrawLatex( 1490,tex_ymax*0.85 , "q < 300 MeV/c" );
      TLegend *leg = new TLegend(0.25,0.6,0.55,0.8);
      leg->AddEntry(gMIXErrorSp_CS_mev[iqlh],"#pi^{-}#Sigma^{+}");
      leg->AddEntry(gMIXErrorSm_CS_mev[iqlh],"#pi^{+}#Sigma^{-}");
      leg->SetLineWidth(0);
      leg->SetFillStyle(0);
      leg->Draw();
    }else if(iqlh==1){
      //gr_S1385_ToSqhi_mev->Draw("5");
      TLine *p = new TLine(1300,0,1605,0);
      p->SetLineColor(1);
      //p->SetLineWidth(2.0);
      p->SetLineStyle(2);
      p->Draw();
      TLine *pkp = new TLine((kpMass+pMass)*1000,0,(kpMass+pMass)*1000,gIMnpipi_Sp_cs_Etotal_mev[0]->GetHistogram()->GetMaximum());
      pkp->SetLineColor(1);
      pkp->SetLineStyle(2);
      pkp->Draw();

      TLine *pk0n = new TLine((k0Mass+nMass)*1000,0,(k0Mass+nMass)*1000,gIMnpipi_Sp_cs_Etotal_mev[0]->GetHistogram()->GetMaximum());
      pk0n->SetLineColor(1);
      pk0n->SetLineStyle(2);
      pk0n->Draw(); 

      TLatex *tex = new TLatex();
      double tex_ymax = gIMnpipi_Sp_cs_Etotal_mev[1]->GetHistogram()->GetMaximum();
      tex->SetTextSize(0.06);
      tex->SetTextColor(1);
      tex->DrawLatex( 1320,tex_ymax*0.85 , "(b)" );
      tex->DrawLatex( 1445,tex_ymax*0.85 , "300 < q < 650 MeV/c" );
    }
  }
  cboth_mev_side->SaveAs("cboth_mev.pdf");

  TCanvas *cAvgmev_side = new TCanvas("cAvgmev_side","cAvgmevv_side",2200,1050);
  cAvgmev_side->SetBottomMargin(0.23);
  cAvgmev_side->SetTopMargin(-0.50);
  cAvgmev_side->SetLeftMargin(0.21);
  //cAvgmev_side->SetRightMargin(0.08);
  for(int iqlh=0;iqlh<2;iqlh++){
    gIMnpipi_SpSmAvg_cs_Etotal_mev[iqlh]->GetYaxis()->SetRangeUser(-0.2,1.65);
    gMIXErrorSpSmAvg_CS_mev[iqlh]->SetFillStyle(3001);
    gMIXErrorSpSmAvg_CS_mev[iqlh]->SetFillColor(12);
    gMIXErrorSpSmAvg_CS_mev[iqlh]->SetMarkerColor(4);
    gMIXErrorSpSmAvg_CS_mev[iqlh]->SetLineColor(12);
  }
  gIMnpipi_SpSmAvg_cs_Etotal_mev[0]->GetXaxis()->SetLimits(1300,1620);
  
  gIMnpipi_SpSmAvg_cs_Etotal_mev[0]->GetXaxis()->SetLabelSize(0.05);
  gIMnpipi_SpSmAvg_cs_Etotal_mev[0]->GetYaxis()->SetLabelSize(0.05);
  gIMnpipi_SpSmAvg_cs_Etotal_mev[0]->GetXaxis()->SetLabelOffset(0.03);
  gIMnpipi_SpSmAvg_cs_Etotal_mev[0]->GetYaxis()->SetLabelOffset(0.03);
  gIMnpipi_SpSmAvg_cs_Etotal_mev[0]->GetXaxis()->SetTickLength(0.02);
  gIMnpipi_SpSmAvg_cs_Etotal_mev[0]->GetYaxis()->SetTickLength(0.02);
  gIMnpipi_SpSmAvg_cs_Etotal_mev[1]->GetXaxis()->SetTickLength(0.02);
  gIMnpipi_SpSmAvg_cs_Etotal_mev[1]->GetYaxis()->SetTickLength(0.02);
  gIMnpipi_SpSmAvg_cs_Etotal_mev[1]->GetXaxis()->SetLabelSize(0.05);
  gIMnpipi_SpSmAvg_cs_Etotal_mev[1]->GetYaxis()->SetLabelSize(0.05);
  gIMnpipi_SpSmAvg_cs_Etotal_mev[1]->GetXaxis()->SetLabelOffset(0.03);
  gIMnpipi_SpSmAvg_cs_Etotal_mev[1]->GetYaxis()->SetLabelOffset(0.03);
  gIMnpipi_SpSmAvg_cs_Etotal_mev[1]->GetXaxis()->SetLimits(1301,1620);
  gIMnpipi_SpSmAvg_cs_Etotal_mev[0]->GetXaxis()->SetTitle("IM(#pi#Sigma) [MeV/c^{2}]");
  gIMnpipi_SpSmAvg_cs_Etotal_mev[1]->GetXaxis()->SetTitle("IM(#pi#Sigma) [MeV/c^{2}]");
  
  gIMnpipi_SpSmAvg_cs_Etotal_mev[0]->GetXaxis()->SetTitleSize(0.07);
  gIMnpipi_SpSmAvg_cs_Etotal_mev[0]->GetXaxis()->SetTitleOffset(1.4);
  gIMnpipi_SpSmAvg_cs_Etotal_mev[0]->GetYaxis()->SetTitleSize(0.07);
  gIMnpipi_SpSmAvg_cs_Etotal_mev[0]->GetYaxis()->SetTitleOffset(1.2);
  gIMnpipi_SpSmAvg_cs_Etotal_mev[1]->GetXaxis()->SetTitleSize(0.07);
  gIMnpipi_SpSmAvg_cs_Etotal_mev[1]->GetXaxis()->SetTitleOffset(1.4);
  gIMnpipi_SpSmAvg_cs_Etotal_mev[1]->GetYaxis()->SetTitleSize(0.07);
  gIMnpipi_SpSmAvg_cs_Etotal_mev[1]->GetYaxis()->SetTitleOffset(1.6);
  
  gIMnpipi_SpSmAvg_cs_Etotal_mev[0]->GetXaxis()->CenterTitle();
  gIMnpipi_SpSmAvg_cs_Etotal_mev[1]->GetXaxis()->CenterTitle();
  gIMnpipi_SpSmAvg_cs_Etotal_mev[0]->SetTitle("");
  gIMnpipi_SpSmAvg_cs_Etotal_mev[0]->GetYaxis()->SetTitle("d#sigma/dM  [#mub/MeV^{2}]");
  gIMnpipi_SpSmAvg_cs_Etotal_mev[0]->GetYaxis()->CenterTitle();
  cAvgmev_side->Divide(2,1,0,0);
  //cAvgmev_side->Divide(2,1,0.0001,0.001);
  cAvgmev_side->cd(1);
  gIMnpipi_SpSmAvg_cs_Etotal_mev[0]->SetMaximum(1.8);
  gIMnpipi_SpSmAvg_cs_Etotal_mev[1]->SetMaximum(1.8);
  gIMnpipi_SpSmAvg_cs_Etotal_mev[0]->Draw("ap");
  gMIXErrorSpSmAvg_CS_mev[0]->Draw("5");
  gr_S1385_ToSqlow_mev->Draw("5");
  TLine *pmev = new TLine(1300,0,1620,0);
  pmev->SetLineColor(1);
  //p->SetLineWidth(2.0);
  pmev->SetLineStyle(2);
  pmev->Draw();
  {
    TLatex *tex = new TLatex();
    double tex_ymax = gIMnpipi_SpSmAvg_cs_Etotal_mev[0]->GetHistogram()->GetMaximum();
    tex->SetTextSize(0.06);
    tex->SetTextColor(1);
    tex->DrawLatex( 1320,tex_ymax*0.85 , "(a)" );
    tex->DrawLatex( 1480,tex_ymax*0.85 , "q < 300 MeV/c" );
    TLine *pkp = new TLine((kpMass+pMass)*1000,0,(kpMass+pMass)*1000,gIMnpipi_Sp_cs_Etotal_mev[0]->GetHistogram()->GetMaximum());
    pkp->SetLineColor(1);
    pkp->SetLineStyle(2);
    pkp->Draw();

    TLine *pk0n = new TLine((k0Mass+nMass)*1000,0,(k0Mass+nMass)*1000,gIMnpipi_Sp_cs_Etotal_mev[0]->GetHistogram()->GetMaximum());
    pk0n->SetLineColor(1);
    pk0n->SetLineStyle(2);
    pk0n->Draw(); 
    //TArrow *ar3 = new TArrow(1365,tex_ymax*0.85 ,1425,tex_ymax*0.85 ,0.04,"<|>");
    //ar3->SetAngle(20);
    //ar3->SetLineWidth(2);
    //ar3->Draw();
  }
  cAvgmev_side->cd(2);
  gIMnpipi_SpSmAvg_cs_Etotal_mev[1]->Draw("ap");
  gMIXErrorSpSmAvg_CS_mev[1]->Draw("5");
  gr_S1385_ToSqhi_mev->Draw("5");
  pmev->Draw();
  {
    TLatex *tex = new TLatex();
    double tex_ymax = gIMnpipi_SpSmAvg_cs_Etotal_mev[0]->GetHistogram()->GetMaximum();
    tex->SetTextSize(0.06);
    tex->SetTextColor(1);
    tex->DrawLatex( 1320,tex_ymax*0.85 , "(b)" );
    tex->DrawLatex( 1445,tex_ymax*0.85 , "300 < q < 650 MeV/c" );
    TLine *pkp = new TLine((kpMass+pMass)*1000,0,(kpMass+pMass)*1000,gIMnpipi_Sp_cs_Etotal_mev[0]->GetHistogram()->GetMaximum());
    pkp->SetLineColor(1);
    pkp->SetLineStyle(2);
    pkp->Draw();

    TLine *pk0n = new TLine((k0Mass+nMass)*1000,0,(k0Mass+nMass)*1000,gIMnpipi_Sp_cs_Etotal_mev[0]->GetHistogram()->GetMaximum());
    pk0n->SetLineColor(1);
    pk0n->SetLineStyle(2);
    pk0n->Draw(); 
  }
  cAvgmev_side->SaveAs("cAvgmev_side.pdf","PDF");

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
  TH1D* CosL1520[3][3];//deco sys, mix sys
  TH1D* qL1405[3][3];//deco sys, mix sys
  TH1D* qL1520[3][3];//deco sys, mix sys
  TH1D* IMnpipi_SpSmAvgCosCut[3][3][2];//deco sys, mix sys, cos cut 0.95-1


  for(int idecosys=0;idecosys<3;idecosys++){
    for(int imixsys=0;imixsys<3;imixsys++){
      Cosn_IMnpipi_Sp_cs[idecosys][imixsys] = (TH2D*)fpisigma[0][imixsys]->Get(Form("Cosn_IMnpipi_Sp_cs_sys%d",idecosys-1));
      Cosn_IMnpipi_Sm_cs[idecosys][imixsys] = (TH2D*)fpisigma[0][imixsys]->Get(Form("Cosn_IMnpipi_Sm_cs_sys%d",idecosys-1));
      Cosn_IMnpipi_SpSmSum[idecosys][imixsys] = (TH2D*)fpisigma[0][imixsys]->Get(Form("Cosn_IMnpipi_SpSmSum_sys%d",idecosys-1));
      CosL1405[idecosys][imixsys] = (TH1D*)fpisigma[0][imixsys]->Get(Form("CosL1405_sys%d",idecosys-1));
      CosL1520[idecosys][imixsys] = (TH1D*)fpisigma[0][imixsys]->Get(Form("CosL1520_sys%d",idecosys-1));
      qL1405[idecosys][imixsys] = (TH1D*)fpisigma[0][imixsys]->Get(Form("qL1405_sys%d",idecosys-1));
      qL1520[idecosys][imixsys] = (TH1D*)fpisigma[0][imixsys]->Get(Form("qL1520_sys%d",idecosys-1));
      IMnpipi_SpSmAvgCosCut[idecosys][imixsys][0] = (TH1D*)fpisigma[0][imixsys]->Get(Form("IMnpipi_SpSmAvgCosCut%d_0",idecosys-1));
      IMnpipi_SpSmAvgCosCut[idecosys][imixsys][1] = (TH1D*)fpisigma[0][imixsys]->Get(Form("IMnpipi_SpSmAvgCosCut%d_1",idecosys-1));
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
  grCosL1405->SetName("grCosL1405");
  TGraphAsymmErrors *grCosL1405_decosysup = new TGraphAsymmErrors(CosL1405[2][1]);
  grCosL1405_decosysup->SetName("grCosL1405_decosysup");
  TGraphAsymmErrors *grCosL1405_decosysdown = new TGraphAsymmErrors(CosL1405[0][1]);
  grCosL1405_decosysdown->SetName("grCosL1405_decosysdown");

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
  gMIXErrorCosL1405->SetName("gMIXErrorCosL1405");

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
  gry->SetLineColor(2);
  gry->SetLineWidth(2);
  gry->Draw("c");

  TFile *fele = TFile::Open("feleangle.root");
  fele->cd();
  //TGraphAsymmErrors *grsumn = (TGraphAsymmErrors*)fele->Get("grsumn");
  //grsumn->Draw("4p");
  


  TGraphAsymmErrors *grCosL1520 = new TGraphAsymmErrors(CosL1520[1][1]);
  grCosL1520->SetName("grCosL1520");
  TGraphAsymmErrors *grCosL1520_decosysup = new TGraphAsymmErrors(CosL1520[2][1]);
  TGraphAsymmErrors *grCosL1520_decosysdown = new TGraphAsymmErrors(CosL1520[0][1]);
  for(int ip=0;ip<grCosL1520->GetN();ip++){
    double valup = grCosL1520_decosysup->GetPointY(ip);
    double valdef = grCosL1520->GetPointY(ip);
    double yeh = grCosL1520->GetErrorYhigh(ip);
    double yel = grCosL1520->GetErrorYlow(ip);
    double valdown = grCosL1520_decosysdown->GetPointY(ip);
    double diffup = valup - valdef;
    double diffdown = valdown - valdef;
     
    if(diffup > diffdown){
      yeh = sqrt(yeh*yeh + diffup*diffup);
      yel = sqrt(yel*yel + diffdown*diffdown);
    }else{
      yeh = sqrt(yeh*yeh + diffdown*diffdown);
      yel = sqrt(yel*yel + diffup*diffup);
    }
    grCosL1520->SetPointEYhigh(ip,yeh);
    grCosL1520->SetPointEYlow(ip,yel);
  }
  TGraphAsymmErrors *gMIXErrorCosL1520 = new TGraphAsymmErrors(CosL1520[1][1]);
  gMIXErrorCosL1520->SetName("gMIXErrorCosL1520");
  TGraphAsymmErrors *gMIXErrorCosL1520_sysup = new TGraphAsymmErrors(CosL1520[1][2]);
  TGraphAsymmErrors *gMIXErrorCosL1520_sysdown = new TGraphAsymmErrors(CosL1520[1][0]);
  for(int ip=0;ip<gMIXErrorCosL1520->GetN();ip++){
    double valup = gMIXErrorCosL1520_sysup->GetPointY(ip);
    double valdef = gMIXErrorCosL1520->GetPointY(ip);
    double yeh = gMIXErrorCosL1520->GetErrorYhigh(ip);
    double yel = gMIXErrorCosL1520->GetErrorYlow(ip);
    double valdown = gMIXErrorCosL1520_sysdown->GetPointY(ip);
    double diffup = valup - valdef;
    double diffdown = valdown - valdef;
     
    if(diffup > diffdown){
      yeh = fabs(diffup);
      yel = fabs(diffdown);
    }else{
      yeh = fabs(diffdown);
      yel = fabs(diffup);
    }
    gMIXErrorCosL1520->SetPointEYhigh(ip,yeh);
    gMIXErrorCosL1520->SetPointEYlow(ip,yel);
    
    double xe = CosL1520[1][1]->GetBinWidth(ip+1);
    gMIXErrorCosL1520->SetPointEXhigh(ip,xe/4);
    gMIXErrorCosL1520->SetPointEXlow(ip,xe/4);
  }
  
  TCanvas *cCosL1520 = new TCanvas("cCosL1520","cCosL1520",1200,800);
  grCosL1520->GetXaxis()->SetRangeUser(0.6,1);
  grCosL1520->GetXaxis()->SetTitle("cos#theta_{n} (CM)");
  grCosL1520->GetXaxis()->CenterTitle();
  grCosL1520->GetYaxis()->SetTitle("d#sigma / dcos#theta_{n} [#mub]  ");
  grCosL1520->GetYaxis()->CenterTitle();
  grCosL1520->Draw("ap");
  gMIXErrorCosL1520->SetLineColor(12);
  gMIXErrorCosL1520->SetFillStyle(3001);
  gMIXErrorCosL1520->SetFillColor(12);
  gMIXErrorCosL1520->Draw("5");
  TLine *pL1520 = new TLine(0.6,0,1,0);
  pL1520->SetLineColor(1);
  //p->SetLineWidth(2.0);
  pL1520->SetLineStyle(2);
  pL1520->Draw();

  TGraphAsymmErrors *grqL1405 = new TGraphAsymmErrors(qL1405[1][1]);
  TGraphAsymmErrors *grqL1405_decosysup = new TGraphAsymmErrors(qL1405[2][1]);
  TGraphAsymmErrors *grqL1405_decosysdown = new TGraphAsymmErrors(qL1405[0][1]);
  for(int ip=0;ip<grqL1405->GetN();ip++){
    double valup = grqL1405_decosysup->GetPointY(ip);
    double valdef = grqL1405->GetPointY(ip);
    double yeh = grqL1405->GetErrorYhigh(ip);
    double yel = grqL1405->GetErrorYlow(ip);
    double valdown = grqL1405_decosysdown->GetPointY(ip);
    double diffup = valup - valdef;
    double diffdown = valdown - valdef;
     
    if(diffup > diffdown){
      yeh = sqrt(yeh*yeh + diffup*diffup);
      yel = sqrt(yel*yel + diffdown*diffdown);
    }else{
      yeh = sqrt(yeh*yeh + diffdown*diffdown);
      yel = sqrt(yel*yel + diffup*diffup);
    }
    grqL1405->SetPointEYhigh(ip,yeh);
    grqL1405->SetPointEYlow(ip,yel);
  }
  TGraphAsymmErrors *gMIXErrorqL1405 = new TGraphAsymmErrors(qL1405[1][1]);
  TGraphAsymmErrors *gMIXErrorqL1405_sysup = new TGraphAsymmErrors(qL1405[1][2]);
  TGraphAsymmErrors *gMIXErrorqL1405_sysdown = new TGraphAsymmErrors(qL1405[1][0]);
  for(int ip=0;ip<gMIXErrorqL1405->GetN();ip++){
    double valup = gMIXErrorqL1405_sysup->GetPointY(ip);
    double valdef = gMIXErrorqL1405->GetPointY(ip);
    double yeh = gMIXErrorqL1405->GetErrorYhigh(ip);
    double yel = gMIXErrorqL1405->GetErrorYlow(ip);
    double valdown = gMIXErrorqL1405_sysdown->GetPointY(ip);
    double diffup = valup - valdef;
    double diffdown = valdown - valdef;
     
    if(diffup > diffdown){
      yeh = fabs(diffup);
      yel = fabs(diffdown);
    }else{
      yeh = fabs(diffdown);
      yel = fabs(diffup);
    }
    gMIXErrorqL1405->SetPointEYhigh(ip,yeh);
    gMIXErrorqL1405->SetPointEYlow(ip,yel);
    
    double xe = qL1405[1][1]->GetBinWidth(ip+1);
    gMIXErrorqL1405->SetPointEXhigh(ip,xe/4);
    gMIXErrorqL1405->SetPointEXlow(ip,xe/4);
  }
  
  TCanvas *cqL1405 = new TCanvas("cqL1405","cqL1405",1200,800);
  grqL1405->SetTitle("CS L1405 M 1400-1440");
  grqL1405->GetXaxis()->SetRangeUser(0.,0.65);
  grqL1405->GetXaxis()->SetTitle("Mom. Transfer [GeV/c]");
  grqL1405->GetXaxis()->CenterTitle();
  grqL1405->GetYaxis()->SetTitle("d#sigma / dq [#mub/(MeV/c)]  ");
  grqL1405->GetYaxis()->CenterTitle();
  grqL1405->Draw("ap");
  gMIXErrorqL1405->SetLineColor(12);
  gMIXErrorqL1405->SetFillStyle(3001);
  gMIXErrorqL1405->SetFillColor(12);
  gMIXErrorqL1405->Draw("5");
  TLine *pqL1405 = new TLine(0.0,0,0.65,0);
  pqL1405->SetLineColor(1);
  //p->SetLineWidth(2.0);
  pqL1405->SetLineStyle(2);
  pqL1405->Draw();


  TGraphAsymmErrors *grqL1520 = new TGraphAsymmErrors(qL1520[1][1]);
  TGraphAsymmErrors *grqL1520_decosysup = new TGraphAsymmErrors(qL1520[2][1]);
  TGraphAsymmErrors *grqL1520_decosysdown = new TGraphAsymmErrors(qL1520[0][1]);
  for(int ip=0;ip<grqL1520->GetN();ip++){
    double valup = grqL1520_decosysup->GetPointY(ip);
    double valdef = grqL1520->GetPointY(ip);
    double yeh = grqL1520->GetErrorYhigh(ip);
    double yel = grqL1520->GetErrorYlow(ip);
    double valdown = grqL1520_decosysdown->GetPointY(ip);
    double diffup = valup - valdef;
    double diffdown = valdown - valdef;
     
    if(diffup > diffdown){
      yeh = sqrt(yeh*yeh + diffup*diffup);
      yel = sqrt(yel*yel + diffdown*diffdown);
    }else{
      yeh = sqrt(yeh*yeh + diffdown*diffdown);
      yel = sqrt(yel*yel + diffup*diffup);
    }
    grqL1520->SetPointEYhigh(ip,yeh);
    grqL1520->SetPointEYlow(ip,yel);
  }
  TGraphAsymmErrors *gMIXErrorqL1520 = new TGraphAsymmErrors(qL1520[1][1]);
  TGraphAsymmErrors *gMIXErrorqL1520_sysup = new TGraphAsymmErrors(qL1520[1][2]);
  TGraphAsymmErrors *gMIXErrorqL1520_sysdown = new TGraphAsymmErrors(qL1520[1][0]);
  for(int ip=0;ip<gMIXErrorqL1520->GetN();ip++){
    double valup = gMIXErrorqL1520_sysup->GetPointY(ip);
    double valdef = gMIXErrorqL1520->GetPointY(ip);
    double yeh = gMIXErrorqL1520->GetErrorYhigh(ip);
    double yel = gMIXErrorqL1520->GetErrorYlow(ip);
    double valdown = gMIXErrorqL1520_sysdown->GetPointY(ip);
    double diffup = valup - valdef;
    double diffdown = valdown - valdef;
     
    if(diffup > diffdown){
      yeh = fabs(diffup);
      yel = fabs(diffdown);
    }else{
      yeh = fabs(diffdown);
      yel = fabs(diffup);
    }
    gMIXErrorqL1520->SetPointEYhigh(ip,yeh);
    gMIXErrorqL1520->SetPointEYlow(ip,yel);
    
    double xe = qL1520[1][1]->GetBinWidth(ip+1);
    gMIXErrorqL1520->SetPointEXhigh(ip,xe/4);
    gMIXErrorqL1520->SetPointEXlow(ip,xe/4);
  }
  
  TCanvas *cqL1520 = new TCanvas("cqL1520","cqL1520",1200,800);
  grqL1520->SetTitle("CS L1520 M 1500-1545");
  grqL1520->GetXaxis()->SetRangeUser(0.0,0.65);
  grqL1520->GetXaxis()->SetTitle("Mom. Transfer [GeV/c]");
  grqL1520->GetXaxis()->CenterTitle();
  grqL1520->GetYaxis()->SetTitle("d#sigma / dq [#mub/(MeV/c)]  ");
  grqL1520->GetYaxis()->CenterTitle();
  grqL1520->Draw("ap");
  gMIXErrorqL1520->SetLineColor(12);
  gMIXErrorqL1520->SetFillStyle(3001);
  gMIXErrorqL1520->SetFillColor(12);
  gMIXErrorqL1520->Draw("5");
  TLine *pqL1520 = new TLine(0.0,0,0.65,0);
  pqL1520->SetLineColor(1);
  //p->SetLineWidth(2.0);
  pqL1520->SetLineStyle(2);
  pqL1520->Draw();


  TH1D* CS_CosS1385Lpim = (TH1D*)flpimnofit->Get("CS_CosS1385Lpim");
  CS_CosS1385Lpim->Print("base");
  TH1D* CS_qS1385Lpim = (TH1D*)flpimnofit->Get("CS_qS1385Lpim");
  CS_qS1385Lpim->Print("base");
  
  const double br_s1385ToLambdapi = 0.87;
  const double br_s1385TopiSigma = 0.117;
  const double br_s1385TopiSigma_err = 0.015;
  const double IsospinCGFactor = 2.0;  
  
  CS_CosS1385Lpim->Scale(1./IsospinCGFactor);
  CS_qS1385Lpim->Scale(1./IsospinCGFactor);

  TCanvas *cCosdistSPDwave = new TCanvas("cCosdistSPDwave","cCosdistSPDwave");

  grCosL1405->Draw("ap");
  gMIXErrorCosL1405->Draw("5");
  grCosL1520->SetLineColor(3);
  grCosL1520->SetMarkerColor(3);
  grCosL1520->Draw("p");
  gMIXErrorCosL1520->SetLineColor(3);
  gMIXErrorCosL1520->SetFillColor(3);
  gMIXErrorCosL1520->Draw("5");
  CS_CosS1385Lpim->SetMarkerColor(5);
  CS_CosS1385Lpim->SetLineColor(5);
  CS_CosS1385Lpim->Draw("same");
  pL1405->Draw();
  gry->Draw("c");
  //gPad->SetLogy();
  

  TCanvas *cqdistSPDwave = new TCanvas("cqdistSPDwave","cqdistSPDwave");
  grqL1405->RemovePoint(0);
  grqL1405->RemovePoint(0);
  grqL1405->RemovePoint(0);
  grqL1405->RemovePoint(0);
  grqL1405->Draw("ap");
  gMIXErrorqL1405->RemovePoint(0);
  gMIXErrorqL1405->RemovePoint(0);
  gMIXErrorqL1405->RemovePoint(0);
  gMIXErrorqL1405->RemovePoint(0);
  gMIXErrorqL1405->Draw("5");
  grqL1520->RemovePoint(0);
  grqL1520->RemovePoint(0);
  grqL1520->SetLineColor(3);
  grqL1520->SetMarkerColor(3);
  grqL1520->Draw("p");
  gMIXErrorqL1520->RemovePoint(0);
  gMIXErrorqL1520->RemovePoint(0);
  gMIXErrorqL1520->SetLineColor(3);
  gMIXErrorqL1520->SetFillColor(3);
  gMIXErrorqL1520->Draw("5");
  CS_qS1385Lpim->SetMarkerColor(5);
  CS_qS1385Lpim->SetLineColor(5);
  CS_qS1385Lpim->Draw("same");
  pqL1405->Draw();


  TGraphAsymmErrors *grIMnpipi_SpSmAvgCosCut0 = new TGraphAsymmErrors(IMnpipi_SpSmAvgCosCut[1][1][0]);
  TGraphAsymmErrors *grIMnpipi_SpSmAvgCosCut0_decosysup = new TGraphAsymmErrors(IMnpipi_SpSmAvgCosCut[2][1][0]);
  TGraphAsymmErrors *grIMnpipi_SpSmAvgCosCut0_decosysdown = new TGraphAsymmErrors(IMnpipi_SpSmAvgCosCut[0][1][0]);
  for(int ip=0;ip<grIMnpipi_SpSmAvgCosCut0->GetN();ip++){
    double valup = grIMnpipi_SpSmAvgCosCut0_decosysup->GetPointY(ip);
    double valdef = grIMnpipi_SpSmAvgCosCut0->GetPointY(ip);
    double yeh = grIMnpipi_SpSmAvgCosCut0->GetErrorYhigh(ip);
    double yel = grIMnpipi_SpSmAvgCosCut0->GetErrorYlow(ip);
    double valdown = grIMnpipi_SpSmAvgCosCut0_decosysdown->GetPointY(ip);
    double diffup = valup - valdef;
    double diffdown = valdown - valdef;
     
    if(diffup > diffdown){
      yeh = sqrt(yeh*yeh + diffup*diffup);
      yel = sqrt(yel*yel + diffdown*diffdown);
    }else{
      yeh = sqrt(yeh*yeh + diffdown*diffdown);
      yel = sqrt(yel*yel + diffup*diffup);
    }
    grIMnpipi_SpSmAvgCosCut0->SetPointEYhigh(ip,yeh);
    grIMnpipi_SpSmAvgCosCut0->SetPointEYlow(ip,yel);
  }


  TGraphAsymmErrors *gMIXErrorIMnpipi_SpSmAvgCosCut0 = new TGraphAsymmErrors(IMnpipi_SpSmAvgCosCut[1][1][0]);
  TGraphAsymmErrors *gMIXErrorIMnpipi_SpSmAvgCosCut0_sysup = new TGraphAsymmErrors(IMnpipi_SpSmAvgCosCut[1][2][0]);
  TGraphAsymmErrors *gMIXErrorIMnpipi_SpSmAvgCosCut0_sysdown = new TGraphAsymmErrors(IMnpipi_SpSmAvgCosCut[1][0][0]);
  for(int ip=0;ip<gMIXErrorIMnpipi_SpSmAvgCosCut0->GetN();ip++){
    double valup = gMIXErrorIMnpipi_SpSmAvgCosCut0_sysup->GetPointY(ip);
    double valdef = gMIXErrorIMnpipi_SpSmAvgCosCut0->GetPointY(ip);
    double yeh = gMIXErrorIMnpipi_SpSmAvgCosCut0->GetErrorYhigh(ip);
    double yel = gMIXErrorIMnpipi_SpSmAvgCosCut0->GetErrorYlow(ip);
    double valdown = gMIXErrorIMnpipi_SpSmAvgCosCut0_sysdown->GetPointY(ip);
    double diffup = valup - valdef;
    double diffdown = valdown - valdef;
     
    if(diffup > diffdown){
      yeh = fabs(diffup);
      yel = fabs(diffdown);
    }else{
      yeh = fabs(diffdown);
      yel = fabs(diffup);
    }
    gMIXErrorIMnpipi_SpSmAvgCosCut0->SetPointEYhigh(ip,yeh);
    gMIXErrorIMnpipi_SpSmAvgCosCut0->SetPointEYlow(ip,yel);
    
    double xe = IMnpipi_SpSmAvgCosCut[1][1][0]->GetBinWidth(ip+1);
    gMIXErrorIMnpipi_SpSmAvgCosCut0->SetPointEXhigh(ip,xe/4);
    gMIXErrorIMnpipi_SpSmAvgCosCut0->SetPointEXlow(ip,xe/4);
  }
  
  {
    double n = grIMnpipi_SpSmAvgCosCut0->GetN();
    double *yval = grIMnpipi_SpSmAvgCosCut0->GetEYlow();
    for(int ip=n;ip>=0;ip--){
      if(yval[ip]<=0.00001){
        grIMnpipi_SpSmAvgCosCut0->RemovePoint(ip);
        gMIXErrorIMnpipi_SpSmAvgCosCut0->RemovePoint(ip);
      }
    }
  }
  
  TCanvas *cIMnpipi_SpSmAvgCosCut0 = new TCanvas("cIMnpipi_SpSmAvgCosCut0","cIMnpipi_SpSmAvgCosCut0",1200,800);
  grIMnpipi_SpSmAvgCosCut0->SetTitle("");
  grIMnpipi_SpSmAvgCosCut0->GetXaxis()->SetRangeUser(1.3,1.6);
  grIMnpipi_SpSmAvgCosCut0->GetXaxis()->SetTitle("IM(#pi#Sigma) [GeV/c^{2}]");
  grIMnpipi_SpSmAvgCosCut0->GetXaxis()->CenterTitle();
  grIMnpipi_SpSmAvgCosCut0->GetYaxis()->SetTitle("d#sigma/dM [#mub/MeV^{2}]");
  grIMnpipi_SpSmAvgCosCut0->GetYaxis()->CenterTitle();
  grIMnpipi_SpSmAvgCosCut0->Draw("ap");
  gMIXErrorIMnpipi_SpSmAvgCosCut0->SetLineColor(12);
  gMIXErrorIMnpipi_SpSmAvgCosCut0->SetFillStyle(3001);
  gMIXErrorIMnpipi_SpSmAvgCosCut0->SetFillColor(12);
  gMIXErrorIMnpipi_SpSmAvgCosCut0->Draw("5");
  TLine *pIMnpipi_SpSmAvgCosCut0 = new TLine(1.3,0,1.6,0);
  pIMnpipi_SpSmAvgCosCut0->SetLineColor(1);
  //p->SetLineWidth(2.0);
  pIMnpipi_SpSmAvgCosCut0->SetLineStyle(2);
  pIMnpipi_SpSmAvgCosCut0->Draw();


  TGraphAsymmErrors *grIMnpipi_SpSmAvgCosCut1 = new TGraphAsymmErrors(IMnpipi_SpSmAvgCosCut[1][1][1]);
  TGraphAsymmErrors *grIMnpipi_SpSmAvgCosCut1_decosysup = new TGraphAsymmErrors(IMnpipi_SpSmAvgCosCut[2][1][1]);
  TGraphAsymmErrors *grIMnpipi_SpSmAvgCosCut1_decosysdown = new TGraphAsymmErrors(IMnpipi_SpSmAvgCosCut[0][1][1]);
  for(int ip=0;ip<grIMnpipi_SpSmAvgCosCut1->GetN();ip++){
    double valup = grIMnpipi_SpSmAvgCosCut1_decosysup->GetPointY(ip);
    double valdef = grIMnpipi_SpSmAvgCosCut1->GetPointY(ip);
    double yeh = grIMnpipi_SpSmAvgCosCut1->GetErrorYhigh(ip);
    double yel = grIMnpipi_SpSmAvgCosCut1->GetErrorYlow(ip);
    double valdown = grIMnpipi_SpSmAvgCosCut1_decosysdown->GetPointY(ip);
    double diffup = valup - valdef;
    double diffdown = valdown - valdef;
     
    if(diffup > diffdown){
      yeh = sqrt(yeh*yeh + diffup*diffup);
      yel = sqrt(yel*yel + diffdown*diffdown);
    }else{
      yeh = sqrt(yeh*yeh + diffdown*diffdown);
      yel = sqrt(yel*yel + diffup*diffup);
    }
    grIMnpipi_SpSmAvgCosCut1->SetPointEYhigh(ip,yeh);
    grIMnpipi_SpSmAvgCosCut1->SetPointEYlow(ip,yel);
  }


  TGraphAsymmErrors *gMIXErrorIMnpipi_SpSmAvgCosCut1 = new TGraphAsymmErrors(IMnpipi_SpSmAvgCosCut[1][1][1]);
  TGraphAsymmErrors *gMIXErrorIMnpipi_SpSmAvgCosCut1_sysup = new TGraphAsymmErrors(IMnpipi_SpSmAvgCosCut[1][2][1]);
  TGraphAsymmErrors *gMIXErrorIMnpipi_SpSmAvgCosCut1_sysdown = new TGraphAsymmErrors(IMnpipi_SpSmAvgCosCut[1][0][1]);
  for(int ip=0;ip<gMIXErrorIMnpipi_SpSmAvgCosCut1->GetN();ip++){
    double valup = gMIXErrorIMnpipi_SpSmAvgCosCut1_sysup->GetPointY(ip);
    double valdef = gMIXErrorIMnpipi_SpSmAvgCosCut1->GetPointY(ip);
    double yeh = gMIXErrorIMnpipi_SpSmAvgCosCut1->GetErrorYhigh(ip);
    double yel = gMIXErrorIMnpipi_SpSmAvgCosCut1->GetErrorYlow(ip);
    double valdown = gMIXErrorIMnpipi_SpSmAvgCosCut1_sysdown->GetPointY(ip);
    double diffup = valup - valdef;
    double diffdown = valdown - valdef;
     
    if(diffup > diffdown){
      yeh = fabs(diffup);
      yel = fabs(diffdown);
    }else{
      yeh = fabs(diffdown);
      yel = fabs(diffup);
    }
    gMIXErrorIMnpipi_SpSmAvgCosCut1->SetPointEYhigh(ip,yeh);
    gMIXErrorIMnpipi_SpSmAvgCosCut1->SetPointEYlow(ip,yel);
    
    double xe = IMnpipi_SpSmAvgCosCut[1][1][1]->GetBinWidth(ip+1);
    gMIXErrorIMnpipi_SpSmAvgCosCut1->SetPointEXhigh(ip,xe/4);
    gMIXErrorIMnpipi_SpSmAvgCosCut1->SetPointEXlow(ip,xe/4);
  }
  
  {
    double n = grIMnpipi_SpSmAvgCosCut1->GetN();
    double *yval = grIMnpipi_SpSmAvgCosCut1->GetEYlow();
    for(int ip=n;ip>=0;ip--){
      if(yval[ip]<=0.00001){
        grIMnpipi_SpSmAvgCosCut1->RemovePoint(ip);
        gMIXErrorIMnpipi_SpSmAvgCosCut1->RemovePoint(ip);
      }
    }
  }
  
  TCanvas *cIMnpipi_SpSmAvgCosCut1 = new TCanvas("cIMnpipi_SpSmAvgCosCut1","cIMnpipi_SpSmAvgCosCut1",1200,800);
  grIMnpipi_SpSmAvgCosCut1->SetTitle("");
  grIMnpipi_SpSmAvgCosCut1->GetXaxis()->SetRangeUser(1.3,1.6);
  grIMnpipi_SpSmAvgCosCut1->GetXaxis()->SetTitle("IM(#pi#Sigma) [GeV/c^{2}]");
  grIMnpipi_SpSmAvgCosCut1->GetXaxis()->CenterTitle();
  grIMnpipi_SpSmAvgCosCut1->GetYaxis()->SetTitle("d#sigma/dM [#mub/(MeV/c^{2}])");
  grIMnpipi_SpSmAvgCosCut1->GetYaxis()->CenterTitle();
  grIMnpipi_SpSmAvgCosCut1->Draw("ap");
  gMIXErrorIMnpipi_SpSmAvgCosCut1->SetLineColor(12);
  gMIXErrorIMnpipi_SpSmAvgCosCut1->SetFillStyle(3001);
  gMIXErrorIMnpipi_SpSmAvgCosCut1->SetFillColor(12);
  gMIXErrorIMnpipi_SpSmAvgCosCut1->Draw("5");
  TLine *pIMnpipi_SpSmAvgCosCut1 = new TLine(1.3,0,1.6,0);
  pIMnpipi_SpSmAvgCosCut1->SetLineColor(1);
  //p->SetLineWidth(2.0);
  pIMnpipi_SpSmAvgCosCut1->SetLineStyle(2);
  pIMnpipi_SpSmAvgCosCut1->Draw();

  TCanvas *cIMnpipi_SpSmAvgCosCut = new TCanvas("cIMnpipi_SpSmAvgCosCut","cIMnpipi_SpSmAvgCosCut",1600,800);
  cIMnpipi_SpSmAvgCosCut->SetBottomMargin(0.12);
  cIMnpipi_SpSmAvgCosCut->SetLeftMargin(0.14);
  cIMnpipi_SpSmAvgCosCut->SetRightMargin(0.08);
  cIMnpipi_SpSmAvgCosCut->Divide(2,1,0.,0.0);
  cIMnpipi_SpSmAvgCosCut->cd(1);
  grIMnpipi_SpSmAvgCosCut0->GetYaxis()->SetRangeUser(-0.1,1.4);
  grIMnpipi_SpSmAvgCosCut0->Draw("ap");
  gMIXErrorIMnpipi_SpSmAvgCosCut0->Draw("5");
  pIMnpipi_SpSmAvgCosCut0->Draw();
  pkp->Draw();
  pk0n->Draw();
  cIMnpipi_SpSmAvgCosCut->cd(2);
  grIMnpipi_SpSmAvgCosCut1->GetYaxis()->SetRangeUser(-0.1,1.4);
  grIMnpipi_SpSmAvgCosCut1->Draw("ap");
  gMIXErrorIMnpipi_SpSmAvgCosCut1->Draw("5");
  pIMnpipi_SpSmAvgCosCut1->Draw();
  pkp->Draw();
  pk0n->Draw();
  

  TCanvas *cq_IMnpipi_SpSmAvg = new TCanvas("cq_IMnpipi_SpSmAvg","cq_IMnpipi_SpSmAvg",1400,1000);
  cq_IMnpipi_SpSmAvg->SetBottomMargin(0.12);
  cq_IMnpipi_SpSmAvg->SetLeftMargin(0.10);
  cq_IMnpipi_SpSmAvg->SetRightMargin(0.10);

  q_IMnpipi_SpSmSum[0][0][1]->Scale(1./2);
  q_IMnpipi_SpSmSum[0][0][1]->SetTitle("#pi^{+}#Sigma^{-} + #pi^{-}#Sigma^{+} charge avg.");
  q_IMnpipi_SpSmSum[0][0][1]->GetXaxis()->SetTitle("IM(#pi#Sigma) [GeV/c^{2}]");
  q_IMnpipi_SpSmSum[0][0][1]->SetMinimum(0);
  q_IMnpipi_SpSmSum[0][0][1]->GetZaxis()->SetMaxDigits(3);
  q_IMnpipi_SpSmSum[0][0][1]->GetZaxis()->SetLabelSize(0.03);
  q_IMnpipi_SpSmSum[0][0][1]->GetXaxis()->SetRangeUser(1.3,1.6);
  q_IMnpipi_SpSmSum[0][0][1]->GetYaxis()->SetRangeUser(0,0.65);
  q_IMnpipi_SpSmSum[0][0][1]->Draw("colz");
  mgSp->Draw("c");
  
  const double Kp_mass = pMass + kpMass;
  const double Kn_mass = nMass + k0Mass;
  TF1 *fkp = new TF1("f1", "sqrt(((x*x-[0]*[0]-[1]*[1])/(2*[0]))*((x*x-[0]*[0]-[1]*[1])/(2*[0]))-[1]*[1])",Kp_mass+0.001,1.8);
  fkp->SetParameter(0,nMass);
  fkp->SetParameter(1,kpMass);

  //fkp->SetLineColor(4);
  fkp->SetLineWidth(5);
  fkp->SetLineStyle(4);
  fkp->SetLineColorAlpha(kPink, 0.35);
  fkp->Draw("same");
   
  /*
  TF1 *fk0 = new TF1("f2", "sqrt(((x*x-[0]*[0]-[1]*[1])/(2*[0]))*((x*x-[0]*[0]-[1]*[1])/(2*[0]))-[1]*[1])",Kn_mass-0.0005,1.8);
  fk0->SetParameter(0,pMass);
  fk0->SetParameter(1,k0Mass);

  //fkp->SetLineColor(4);
  fk0->SetLineWidth(5);
  fk0->SetLineStyle(4);
  fk0->SetLineColorAlpha(kBlue, 0.35);
  fk0->Draw("same");
  */

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
  grCosL1405->Write();
  gMIXErrorCosL1405->Write();
  grCosL1520->Write();
  gMIXErrorCosL1520->Write();
  CS_CosS1385Lpim->Write();


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


