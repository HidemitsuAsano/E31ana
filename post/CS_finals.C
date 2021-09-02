void CS_finals()
{

  TFile *fpisigma = TFile::Open("cs_pisigma.root","READ");
  TFile *flpim    = TFile::Open("cs_lpim.root","READ");

  const double br_s1385TopiSigma = 0.117;
  const double br_s1385TopiSigma_err = 0.015;
  const double br_SpToNpi = 0.4831;
  const double br_SpToNpi_err = 0.003;
  const double br_SmToNpi = 0.99848;
  const double br_SmToNpi_err = 0.00005;
  
  TH2D* q_IMnpipi_Sp_cs[2];
  TH2D* q_IMnpipi_Sm_cs[2];
  TH2D* q_IMnpipi_K0_cs[2];
  TH1D* IMnpipi_Sp_cs[2];
  TH1D* IMnpipi_Sm_cs[2];
  TH1D* IMnpipi_K0_cs[2];
  TH2D* CS_q_IMppipi_p_wL_sum;
  TH1D* CS_IMppipi_p_wL_sum_0;
  TH1D* CS_IMppipi_p_wL_sum_350;

  for(int iq=0;iq<2;iq++){
    q_IMnpipi_Sp_cs[iq] = (TH2D*)fpisigma->Get(Form("q_IMnpipi_Sp_cs%d",iq));
    q_IMnpipi_Sm_cs[iq] = (TH2D*)fpisigma->Get(Form("q_IMnpipi_Sm_cs%d",iq));
    q_IMnpipi_K0_cs[iq] = (TH2D*)fpisigma->Get(Form("q_IMnpipi_K0_cs%d",iq));
    IMnpipi_Sp_cs[iq] = (TH1D*)fpisigma->Get(Form("IMnpipi_Sp_cs%d",iq));
    IMnpipi_Sm_cs[iq] = (TH1D*)fpisigma->Get(Form("IMnpipi_Sm_cs%d",iq));
    IMnpipi_K0_cs[iq] = (TH1D*)fpisigma->Get(Form("IMnpipi_K0_cs%d",iq));
  }
  CS_q_IMppipi_p_wL_sum = (TH2D*)flpim->Get("CS_q_IMppipi_p_wL_sum");
  CS_IMppipi_p_wL_sum_0 = (TH1D*)flpim->Get("CS_IMppipi_p_wL_sum_0");
  CS_IMppipi_p_wL_sum_350 = (TH1D*)flpim->Get("CS_IMppipi_p_wL_sum_350");



  TCanvas *c1 = new TCanvas("c1","c1",1000,800);
  c1->cd();
  IMnpipi_Sp_cs[1]->Draw("HE");
  TH1D* CS_IMppipi_p_wL_sum_350_ToSp = (TH1D*)CS_IMppipi_p_wL_sum_350->Clone("CS_IMppipi_p_wL_sum_350_ToSp");
  CS_IMppipi_p_wL_sum_350_ToSp->Scale(br_s1385TopiSigma*br_SpToNpi/3.0); 
  CS_IMppipi_p_wL_sum_350_ToSp->Draw("HEsame"); 


  TCanvas *c2 = new TCanvas("c2","c2",1000,800);
  c2->cd();
  IMnpipi_Sm_cs[1]->Draw("HE");
  TH1D* CS_IMppipi_p_wL_sum_350_ToSm = (TH1D*)CS_IMppipi_p_wL_sum_350->Clone("CS_IMppipi_p_wL_sum_350_ToSm");
  CS_IMppipi_p_wL_sum_350_ToSm->Scale(br_s1385TopiSigma*br_SmToNpi/3.0); 
  CS_IMppipi_p_wL_sum_350_ToSm->Draw("HEsame"); 


}
