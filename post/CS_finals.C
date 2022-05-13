void CS_finals()
{
  if(gROOT->GetVersionInt() < 60000){
    std::cout << "Use ROOT6 !!" << std::endl;
    return;
  }

  const int Version = 241;
  const int dEcut[3]={2,4,6};

  //gStyle->SetErrorX(0.);  
  TFile *fpisigma[3][3];//dEcut, sysud
  for(int iEcut=0;iEcut<3;iEcut++){
    fpisigma[iEcut][1] = TFile::Open(Form("cs_pisigma_v%d_dE%d_sys0.root",Version,dEcut[iEcut]),"READ");
  }
  fpisigma[0][0] = TFile::Open(Form("cs_pisigma_v%d_dE2_sys-1.root",Version),"READ");
  fpisigma[0][2] = TFile::Open(Form("cs_pisigma_v%d_dE2_sys1.root",Version),"READ");

  //TFile *flpim    = TFile::Open("cs_lpim_killcombi.root","READ");
  TFile *flpim = TFile::Open("CSLpimFit.root","READ");

  const double br_s1385ToLambdapi = 0.87;
  const double br_s1385TopiSigma = 0.117;
  const double br_s1385TopiSigma_err = 0.015;
  const double br_SpToNpi = 0.4831;
  const double br_SpToNpi_err = 0.003;
  const double br_SmToNpi = 0.99848;
  const double br_SmToNpi_err = 0.00005;
  
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
  for(int iq=0;iq<4;iq++){
    cqM_Sp[iq] = new TCanvas(Form("cqM_Sp%d",iq),Form("cqM_Sp%d",iq),1000,800);
    q_IMnpipi_Sp_cs[iq][0][0]->SetMinimum(0);
    q_IMnpipi_Sp_cs[iq][0][0]->GetYaxis()->SetRangeUser(0,0.8);
    q_IMnpipi_Sp_cs[iq][0][0]->Draw("colz");
    cqM_Sm[iq] = new TCanvas(Form("cqM_Sm%d",iq),Form("cqM_Sm%d",iq),1000,800);
    q_IMnpipi_Sm_cs[iq][0][0]->SetMinimum(0);
    q_IMnpipi_Sm_cs[iq][0][0]->GetYaxis()->SetRangeUser(0,0.8);
    q_IMnpipi_Sm_cs[iq][0][0]->Draw("colz");
    cqM_K0[iq] = new TCanvas(Form("cqM_K0%d",iq),Form("cqM_K0%d",iq),1000,800);
    q_IMnpipi_K0_cs[iq][0][0]->SetMinimum(0);
    q_IMnpipi_K0_cs[iq][0][0]->GetYaxis()->SetRangeUser(0,0.8);
    q_IMnpipi_K0_cs[iq][0][0]->Draw("colz");
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
  



  //TH2D* CS_q_IMppipi_p_wL_sum;
  //TH1D* CS_IMppipi_p_wL_sum_0;
  //TH1D* CS_IMppipi_p_wL_sum_350;
  //CS_q_IMppipi_p_wL_sum = (TH2D*)flpim->Get("CS_q_IMppipi_p_wL_sum");
  //CS_IMppipi_p_wL_sum_0 = (TH1D*)flpim->Get("CS_IMppipi_p_wL_sum_0");
  //CS_IMppipi_p_wL_sum_350 = (TH1D*)flpim->Get("CS_IMppipi_p_wL_sum_350");
  
  //get Lpim
  TH2F* CS_lpim_sum = (TH2F*)flpim->Get("CS_sum");
  TH2F* CS_lpim_fit = (TH2F*)flpim->Get("Func");

  TCanvas *clpim = new TCanvas("clpim","clpim",1000,800);
  clpim->cd();
  CS_lpim_sum->Draw("colz");
  CS_lpim_fit->Draw("boxsame");
  
  TH2F* CS_lpim_qcut[3];//iq
  for(int iq=0;iq<3;iq++){
    CS_lpim_qcut[iq] = (TH2F*) CS_lpim_fit->Clone(Form("CS_lpim_qcut%d",iq));
    CS_lpim_qcut[iq]->SetMinimum(0);
    CS_lpim_qcut[iq]->SetFillColor(0);
  }
  CS_lpim_qcut[0]->GetYaxis()->SetRangeUser(0,0.6);
  CS_lpim_qcut[1]->GetYaxis()->SetRangeUser(0,0.35);
  CS_lpim_qcut[2]->GetYaxis()->SetRangeUser(0.35,0.6);
  double binwidthq = CS_lpim_qcut[0]->GetYaxis()->GetBinWidth(1)*1000.0; 
  std::cout << "binq width " << binwidthq  << std::endl;
  TH1D* CS_S1385_ToSp[3]; 
  TH1D* CS_S1385_ToSm[3]; 
  TH1D* CS_S1385_ToSpSm[3]; 
   
  //assume C.S. Sigma(1385)- ~ Sigma(1385)0
  for(int iq=0;iq<3;iq++){
    CS_lpim_qcut[iq]->Scale(br_s1385TopiSigma/2.0/br_s1385ToLambdapi*binwidthq);
    CS_S1385_ToSp[iq] = (TH1D*)CS_lpim_qcut[iq]->ProjectionX(Form("CS_S1385_ToSp%d",iq));
    CS_S1385_ToSm[iq] = (TH1D*)CS_lpim_qcut[iq]->ProjectionX(Form("CS_S1385_ToSm%d",iq));
    CS_S1385_ToSpSm[iq] = (TH1D*)CS_lpim_qcut[iq]->ProjectionX(Form("CS_S1385_ToSpSm%d",iq));
    CS_S1385_ToSpSm[iq]->Scale(2.0);
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
      double  valdown = IMnpipi_Sp_cs[iq][0][0]->GetBinContent(ip+1);
      double  valdef = IMnpipi_Sp_cs[iq][0][1]->GetBinContent(ip+1);
      double  valup  = IMnpipi_Sp_cs[iq][0][2]->GetBinContent(ip+1);
      
      double yeh = valup-valdef;
      double yel = valdown-valdef;
         
      gMIXErrorSp_CS[iq]->SetPointEYhigh(ip,fabs(yeh));
      gMIXErrorSp_CS[iq]->SetPointEYlow(ip,fabs(yel));
    }
    for(int ip=0;ip<( gMIXErrorSm_CS[iq]->GetN());ip++){
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
  std::cout << __LINE__ << std::endl;
  TCanvas *cSysSp[4];
  for(int iq=0;iq<4;iq++){  
    cSysSp[iq] = new TCanvas(Form("cSysSp%d",iq),Form("cSysSp%d",iq),1000,800);
    IMnpipi_Sp_cs[iq][0][1]->SetLineColor(1);
    IMnpipi_Sp_cs[iq][0][1]->SetMarkerColor(1);
    IMnpipi_Sp_cs[iq][0][1]->GetXaxis()->SetRangeUser(1.3,1.6);
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
      CS_S1385_ToSp[iq]->SetLineColor(6);
      CS_S1385_ToSp[iq]->Draw("same");
    }
  }

  TCanvas *cSysSm[4];
  for(int iq=0;iq<4;iq++){  
    cSysSm[iq] = new TCanvas(Form("cSysSm%d",iq),Form("cSysSm%d",iq),1000,800);
    IMnpipi_Sm_cs[iq][0][1]->SetLineColor(1);
    IMnpipi_Sm_cs[iq][0][1]->SetMarkerColor(1);
    IMnpipi_Sm_cs[iq][0][1]->GetXaxis()->SetRangeUser(1.3,1.6);
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
    if(iq<3){
      CS_S1385_ToSm[iq]->SetLineColor(6);
      CS_S1385_ToSm[iq]->Draw("same");
    }
  }

  TCanvas *cSysK0[4];
  for(int iq=0;iq<4;iq++){  
    cSysK0[iq] = new TCanvas(Form("cSysK0%d",iq),Form("cSysK0%d",iq),1000,800);
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

  TCanvas *cSysSpSmSum[4];
  for(int iq=0;iq<4;iq++){  
    cSysSpSmSum[iq] = new TCanvas(Form("cSysSpSmSum%d",iq),Form("cSysSpSmSum%d",iq),1000,800);
    IMnpipi_SpSmSum[iq][0][1]->SetLineColor(1);
    IMnpipi_SpSmSum[iq][0][1]->SetMarkerColor(1);
    IMnpipi_SpSmSum[iq][0][1]->GetXaxis()->SetRangeUser(1.3,1.6);
    IMnpipi_SpSmSum[iq][0][1]->Draw("E");
    gDecoErrorSpSm_CS[iq]->Draw("5");
    gMIXErrorSpSm_CS[iq]->SetFillStyle(3002);
    gMIXErrorSpSm_CS[iq]->SetFillColor(4);
    gMIXErrorSpSm_CS[iq]->SetMarkerColor(4);
    gMIXErrorSpSm_CS[iq]->SetLineColor(4);
    gMIXErrorSpSm_CS[iq]->Draw("3");
    if(iq<3){
      CS_S1385_ToSpSm[iq]->SetFillStyle(0);
      CS_S1385_ToSpSm[iq]->SetFillColor(0);
      CS_S1385_ToSpSm[iq]->SetLineColor(6);
      CS_S1385_ToSpSm[iq]->Draw("same");
    }
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


  const double solidAngleCoscut = 2.0*3.1415926535*(1.00-0.99657);//theta 0.5 = 
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
  CS_Spcomp->SetMaximum(35);
  CS_Spcomp->SetYTitle("d^{2}#rho/dM d#Omega [#mu b/MeVsr]");
  CS_Spcomp->Draw();
  CS_SpcompMIX->Draw("3");
  CS_SpcompSrErr->SetFillStyle(3002);
  CS_SpcompSrErr->SetFillColor(5);
  CS_SpcompSrErr->GetXaxis()->SetRangeUser(1.3,1.6);
  CS_SpcompSrErr->Draw("3");
  CS_SpcompDeco->SetLineColor(3);
  CS_SpcompDeco->Draw("5");
  //grinoueSp->Draw("P");
  grinoueSpcs->Draw("p");
  
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
  //CS_Smcomp->Draw();
 
  TCanvas *cthetacompSm = new TCanvas("cthetacompSm","cthetacompSm",1000,800);
  cthetacompSm->cd();
  //CS_IMppipi_p_wL_mc_coscut->Draw("HEsame");
  CS_Smcomp->SetMaximum(10);
  CS_Smcomp->SetYTitle("d^{2}#rho/dMd#Omega [#mu b/MeVsr]");
  CS_Smcomp->Draw();
  CS_SmcompDeco->Draw("5");
  CS_SmcompSrErr->SetFillStyle(3002);
  CS_SmcompSrErr->SetFillColor(5);
  CS_SmcompSrErr->GetXaxis()->SetRangeUser(1.3,1.6);
  CS_SmcompSrErr->Draw("3");
  CS_SmcompMIX->Draw("3");
  //grinoueSm->Draw("P");
  grinoueSmcs->Draw("p");






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
