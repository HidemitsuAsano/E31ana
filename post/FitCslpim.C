const int npartfit=4;

Double_t BWandFormS(Double_t *x,Double_t *par)
{
  Double_t 





}


Double_t BWandFormP(Double_t *x,Double_t *par)
{
  Double_t r1 = par[0]*pow(par[1]/2.0,2.0)/((pow(x[0]-par[2]),2.0)+pow(par[1]/2.0,2.0));
 
  Double_t r2 = pow(x[1]/par[4],2.0)*exp(-pow(x[1]/par[4],2.0));

}



void FitCslpim()
{
  TFile *file = TFile::Open("cs_lpim_killcombi.root");
  file->cd();
 
  TH2F* CS_q_IMppipi_p_wL_nop2 = (TH2F*)file->Get(Form("CS_q_IMppipi_p_wL_nop2_acc%d",0));
  TH2F* CS_q_IMppipi_p_wL_wp2 = (TH2F*)file->Get(Form("CS_q_IMppipi_p_wL_wp2_acc%d",0));
  
  TCanvas *cCS_q_IMppipi_p_wL_nop2 = new TCanvas("cCS_q_IMppipi_p_wL_nop2","cCS_q_IMppipi_p_wL_nop2",1000,800);
  CS_q_IMppipi_p_wL_nop2->Draw("colz");

  TCanvas *cCS_q_IMppipi_p_wL_wp2 = new TCanvas("cCS_q_IMppipi_p_wL_wp2","cCS_q_IMppipi_p_wL_wp2",1000,800);
  CS_q_IMppipi_p_wL_wp2->Draw("colz");
 
  TCanvas *cCS_sum = new TCanvas("cCS_sum","cCS_sum",1000,800);
  //CS_q_IMppipi_p_wL_nop2->SetMaximum(0.0045);
  //CS_q_IMppipi_p_wL_nop2->Draw("colz");
  //CS_q_IMppipi_p_wL_wp2->SetMaximum(CS_q_IMppipi_p_wL_nop2->GetMaximum());
  //CS_q_IMppipi_p_wL_wp2->Draw("colzsame");
  TH2F* CS_sum = (TH2F*)CS_q_IMppipi_p_wL_nop2->Clone("CS_sum");
  CS_sum->SetTitle("CS_sum");
  CS_sum->Add(CS_q_IMppipi_p_wL_wp2);
  CS_sum->Draw("colz");

  TCanvas *cCS_acc_sum = new TCanvas("cCS_acc_sum","cCS_acc_sum",1000,800);
  TH2F* acc_nop2 = (TH2F*)file->Get("q_IMppipi_p_wL_nop2_acc_clean0");
  TH2F* acc_wp2 = (TH2F*)file->Get("q_IMppipi_p_wL_wp2_acc_clean0");
  TH2F* acc_sum = (TH2F*)acc_nop2->Clone("acc_sum");
  acc_sum->Add(acc_wp2);
  acc_sum->SetTitle("acc_sum");
  acc_sum->GetZaxis()->SetRangeUser(0,0.1);
  acc_sum->Draw("colz");
  
  TCanvas *cCS_sum_fit = new TCanvas("cCS_sum_fit","cCS_sum_fit",1000,800);
  TH2F* CS_sum_fit = (TH2F*)CS_sum->Clone("CS_sum_fit");
  CS_sum_fit->SetName("CS_sum_fit");
  CS_sum_fit->SetTitle("CS_sum_fit");
  for(int ix=0;ix<CS_sum_fit->GetNbinsX();ix++){
    for(int iy=0;iy<CS_sum_fit->GetNbinsY();iy++){
      double acccont = acc_sum->GetBinContent(ix,iy);
      if(acccont<0.01){
         CS_sum_fit->SetBinContent(ix,iy,0);
         CS_sum_fit->SetBinError(ix,iy,0);
      }
    }
  }
  CS_sum_fit->Draw("colz");
  
  
   


  TCanvas *cCS_q = new TCanvas("cCS_q","cCS_q",1000,800);
  TH1D* CS_q_nop20 = (TH1D*)file->Get("CS_q_nop20");
  CS_q_nop20->Draw("");
  
  TH1D* CS_q_wp20 = (TH1D*)file->Get("CS_q_wp20");
  CS_q_wp20->Draw("same");
  const double qbinwidth = CS_q_nop20->GetXaxis()->GetBinWidth(2);
  CS_q_nop20->Scale(1./qbinwidth);
  CS_q_wp20->Scale(1./qbinwidth);
  TCanvas *cCS_q_gr = new TCanvas("cCS_q_gr","cCS_q_gr",1000,800);
  TGraphErrors *gr_nop2 = new TGraphErrors();
  int ip=0;
  for(int i=0;i<CS_q_nop20->GetNbinsX();i++){
    double cont = CS_q_nop20->GetBinContent(i);
    double err = CS_q_nop20->GetBinError(i);
    double bincent = CS_q_nop20->GetBinCenter(i);
    if(cont>0.0 && 0.4<bincent && bincent<0.75){
      gr_nop2->SetPoint(ip,bincent,cont);
      gr_nop2->SetPointError(ip,0.015,err);
      ip++;
    }
  }
  //inoue
  gr_nop2->SetPoint(ip,0.315,1.00582900000000008e+00);
  gr_nop2->SetPointError(ip,0.015,0.203361);
  gr_nop2->SetMarkerStyle(20);
  //gr_nop2->SetMarkerColor(20);
  gr_nop2->GetXaxis()->SetLimits(0.0,1.5);
  gr_nop2->SetMinimum(0.0);
  gr_nop2->SetMaximum(CS_q_nop20->GetMaximum());
  gr_nop2->GetXaxis()->SetTitle("Mom. Transfer [GeV/c]");
  gr_nop2->GetXaxis()->CenterTitle();
  gr_nop2->Draw("AP");

  TGraphErrors *gr_wp2 = new TGraphErrors();
  ip=0;
  for(int i=0;i<CS_q_wp20->GetNbinsX();i++){
    double cont = CS_q_wp20->GetBinContent(i);
    double err = CS_q_wp20->GetBinError(i);
    if(cont>0.0){
      double bincent = CS_q_wp20->GetBinCenter(i);
      gr_wp2->SetPoint(ip,bincent,cont);
      gr_wp2->SetPointError(ip,0.015,err);
      ip++;
    }
  }
  gr_wp2->SetMarkerStyle(20);
  gr_wp2->SetMarkerColor(2);
  gr_wp2->SetLineColor(2);
  gr_wp2->Draw("P");
  
  


}
