void CS_finals()
{
  const int Version = 241;
  const int dEcut[3]={2,4,6};

  //gStyle->SetErrorX(0.);  
  TFile *fpisigma[3][3];//dEcut, sysud
  for(int iEcut=0;iEcut<3;iEcut++){
    fpisigma[iEcut][1] = TFile::Open(Form("cs_pisigma_v%d_dE%d_sys0.root",Version,dEcut[iEcut]),"READ");
  }
  fpisigma[0][0] = TFile::Open(Form("cs_pisigma_v%d_dE2_sys-1.root",Version),"READ");
  fpisigma[0][2] = TFile::Open(Form("cs_pisigma_v%d_dE2_sys1.root",Version),"READ");

  TFile *flpim    = TFile::Open("cs_lpim_killcombi.root","READ");

  const double br_s1385ToLambdapi = 0.87;
  const double br_s1385TopiSigma = 0.117;
  const double br_s1385TopiSigma_err = 0.015;
  const double br_SpToNpi = 0.4831;
  const double br_SpToNpi_err = 0.003;
  const double br_SmToNpi = 0.99848;
  const double br_SmToNpi_err = 0.00005;
  

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
  }
  std::cout << "FILE GET " << std::endl;

  TH2D* CS_q_IMppipi_p_wL_sum;
  TH1D* CS_IMppipi_p_wL_sum_0;
  TH1D* CS_IMppipi_p_wL_sum_350;
  CS_q_IMppipi_p_wL_sum = (TH2D*)flpim->Get("CS_q_IMppipi_p_wL_sum");
  CS_IMppipi_p_wL_sum_0 = (TH1D*)flpim->Get("CS_IMppipi_p_wL_sum_0");
  CS_IMppipi_p_wL_sum_350 = (TH1D*)flpim->Get("CS_IMppipi_p_wL_sum_350");

  
  //Syserror calculation
  TGraphAsymmErrors *gMIXErrorSp_CS[4];
  TGraphAsymmErrors *gMIXErrorSm_CS[4];
  TGraphAsymmErrors *gMIXErrorK0_CS[4];
  TGraphAsymmErrors *gdEErrorSp_CS[4];
  TGraphAsymmErrors *gdEErrorSm_CS[4];
  TGraphAsymmErrors *gdEErrorK0_CS[4];
  
  for(int iq=0;iq<4;iq++){
    gMIXErrorSp_CS[iq] = new TGraphAsymmErrors(IMnpipi_Sp_cs[iq][0][1]);
    gMIXErrorSm_CS[iq] = new TGraphAsymmErrors(IMnpipi_Sm_cs[iq][0][1]);
    gMIXErrorK0_CS[iq] = new TGraphAsymmErrors(IMnpipi_K0_cs[iq][0][1]);
    gMIXErrorSp_CS[iq]->SetName(Form("gMIXErrorSp_CS%d",iq));
    gMIXErrorSm_CS[iq]->SetName(Form("gMIXErrorSm_CS%d",iq));
    gMIXErrorK0_CS[iq]->SetName(Form("gMIXErrorK0_CS%d",iq));

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
    IMnpipi_Sp_cs[iq][0][1]->Draw("E");
    //gDecoErrorSp_CS[iq]->Draw("5");
    gdEErrorSp_CS[iq]->SetFillStyle(0);
    gdEErrorSp_CS[iq]->SetMarkerColor(2);
    gdEErrorSp_CS[iq]->SetLineColor(2);
    gdEErrorSp_CS[iq]->Draw("5");
    gMIXErrorSp_CS[iq]->SetFillStyle(3002);
    gMIXErrorSp_CS[iq]->SetFillColor(4);
    gMIXErrorSp_CS[iq]->SetMarkerColor(4);
    gMIXErrorSp_CS[iq]->SetLineColor(4);
    gMIXErrorSp_CS[iq]->Draw("3");
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
