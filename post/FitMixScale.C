//display vicinity of signal region of data and mixed events
//determine scaling factor of mixed events by fitting

TH1D* MMnmiss_woK0_woSm_mix;
TH1D* IMnpip_woK0_woSm_mix;
TH1D* MMnmiss_woK0_woSp_mix;
TH1D* IMnpim_woK0_woSp_mix;

const int dEcut=2;
const int Version=232;

void FitMixScale()
{
  fr = TFile::Open(Form("evanaIMpisigma_npippim_v%d_out_dE%d_iso_nostop.root",Version,dEcut));
  fmix = TFile::Open(Form("evanaIMpisigma_npippim_v%d_MIX_cut4_out_dE%d_iso_nostop_sys0.root",Version,dEcut));

  TH2D* MMnmiss_IMnpip_woK0_woSm_vici_data = (TH2D*)fr->Get("MMnmiss_IMnpip_dE_woK0_woSm_vici");
  TH2D* MMnmiss_IMnpim_woK0_woSp_vici_data = (TH2D*)fr->Get("MMnmiss_IMnpim_dE_woK0_woSp_vici");

  TH2D* MMnmiss_IMnpip_woK0_woSm_vici_mix = (TH2D*)fmix->Get("MMnmiss_IMnpip_dE_woK0_woSm_vici");
  TH2D* MMnmiss_IMnpim_woK0_woSp_vici_mix = (TH2D*)fmix->Get("MMnmiss_IMnpim_dE_woK0_woSp_vici");

  MMnmiss_IMnpip_woK0_woSm_vici_data->GetXaxis()->SetRangeUser(1.1,1.3);
  MMnmiss_IMnpip_woK0_woSm_vici_data->GetYaxis()->SetRangeUser(0.7,1.1);
  MMnmiss_IMnpip_woK0_woSm_vici_mix->GetXaxis()->SetRangeUser(1.1,1.3);
  MMnmiss_IMnpip_woK0_woSm_vici_mix->GetYaxis()->SetRangeUser(0.7,1.1);
  MMnmiss_IMnpim_woK0_woSp_vici_data->GetXaxis()->SetRangeUser(1.1,1.3);
  MMnmiss_IMnpim_woK0_woSp_vici_data->GetYaxis()->SetRangeUser(0.7,1.1);
  MMnmiss_IMnpim_woK0_woSp_vici_mix->GetXaxis()->SetRangeUser(1.1,1.3);
  MMnmiss_IMnpim_woK0_woSp_vici_mix->GetYaxis()->SetRangeUser(0.7,1.1);
  
  TCanvas *c1 = new TCanvas("c1","c1");
  MMnmiss_IMnpip_woK0_woSm_vici_data->Draw("colz");

  TCanvas *c2 = new TCanvas("c2","c2");
  MMnmiss_IMnpim_woK0_woSp_vici_data->Draw("colz");

  TCanvas *c3 = new TCanvas("c3","c3");
  MMnmiss_IMnpip_woK0_woSm_vici_mix->Draw("colz");

  TCanvas *c4 = new TCanvas("c4","c4");
  MMnmiss_IMnpim_woK0_woSp_vici_mix->Draw("colz");
  
  TCanvas *c5 = new TCanvas("c5","c5");
  TH1D* IMnpip_woK0_woSm_data = (TH1D*)MMnmiss_IMnpip_woK0_woSm_vici_data->ProjectionX("IMnpip_woK0_woSm_vici_data");
  IMnpip_woK0_woSm_mix  = (TH1D*)MMnmiss_IMnpip_woK0_woSm_vici_mix->ProjectionX("IMnpip_woK0_woSm_vici_mix");
  
  IMnpip_woK0_woSm_data->Draw();
  TF1 *fIMnpip = new TF1("fIMnpip",fit_IMnpip,-0.2,0.2,1);
  TF1 *fIMnpip_up = new TF1("f1IMnpip_up",fit_IMnpip,-0.2,0.2,1);
  TF1 *fIMnpip_down = new TF1("f1IMnip_down",fit_IMnpip,-0.2,0.2,1);
  fIMnpip->SetParameter(1,1.0);
  IMnpip_woK0_woSm_data->Fit("fIMnpip","r","",-0.1,-0.06);
  Double_t ret = fIMnpip->GetParameter(0);
  Double_t reterr = fIMnpip->GetParError(0);
  fIMnpip->SetParameter(0,ret);
  fIMnpip->SetLineColor(3);
  fIMnpip->Draw("same");


  /*
  double avg = (ret*ret2err+ret2*reterr)/(reterr+ret2err);
  std::cout << "Sp mode scale avg: " << avg << std::endl;
  f1->SetParameter(0,avg);
  f1u->SetParameter(0,avg*1.1);
  f1d->SetParameter(0,avg*0.9);
  f1->SetLineColor(3);
  f1u->SetLineColor(4);
  f1d->SetLineColor(4);
  c5->cd();
  f1->Draw("same");
  f1u->Draw("same");
  f1d->Draw("same");
  f2->SetParameter(0,avg);
  f2u->SetParameter(0,avg*1.1);
  f2d->SetParameter(0,avg*0.9);
  f2->SetLineColor(3);
  f2u->SetLineColor(4);
  f2d->SetLineColor(4);
  c6->cd();
  f2->Draw("same");
  f2u->Draw("same");
  f2d->Draw("same");
  //IMnpip_fmix->SetLineColor(4);
  //IMnpip_fmix->Draw("same");

  TCanvas *c7 = new TCanvas("c7","c7");
  MM2npi_f_pim->GetXaxis()->SetRangeUser(-0.2,0.2);
  MM2npi_f_pim->Draw("HE");

  TF1 *f3 = new TF1("f3",fit_MM2npi_pim,-0.2,0.2,1);
  TF1 *f3u = new TF1("f3u",fit_MM2npi_pim,-0.2,0.2,1);
  TF1 *f3d = new TF1("f3d",fit_MM2npi_pim,-0.2,0.2,1);
  f3->SetParameter(1,1.0);
  MM2npi_f_pim->Fit("f3","r","",-0.1,-0.075);
  Double_t ret3 = f3->GetParameter(0);
  Double_t ret3err = f3->GetParError(0);
  f3->SetParameter(0,ret3);
  f3->SetLineColor(3);
  f3->Draw("same");


  TCanvas *c8 = new TCanvas("c8","c8");
  IMnpim_f->GetXaxis()->SetRangeUser(1.14,1.24);
  IMnpim_f->Draw("HE");

  TF1 *f4 = new TF1("f4",fit_IMnpim,1,2,1);
  TF1 *f4u = new TF1("f4u",fit_IMnpim,1,2,1);
  TF1 *f4d = new TF1("f4d",fit_IMnpim,1,2,1);
  IMnpim_f->Fit("f4","r","",1.21,1.225);
  Double_t ret4 = f4->GetParameter(0);
  Double_t ret4err = f4->GetParError(0);
  f4->SetParameter(0,ret4);
  f4->SetLineColor(3);
  f4->Draw("same");

  double avg2 = (ret3*ret4err+ret4*ret3err)/(ret3err+ret4err);
  std::cout << "Sm mode scale avg: " << avg2 << std::endl;
  f3->SetParameter(0,avg2);
  f3u->SetParameter(0,avg2*1.1);
  f3d->SetParameter(0,avg2*0.9);
  f3->SetLineColor(3);
  f3u->SetLineColor(4);
  f3d->SetLineColor(4);
  c7->cd();
  f3->Draw("same");
  f3u->Draw("same");
  f3d->Draw("same");
  f4->SetParameter(0,avg2);
  f4u->SetParameter(0,avg2*1.1);
  f4d->SetParameter(0,avg2*0.9);
  f4->SetLineColor(3);
  f4u->SetLineColor(4);
  f4d->SetLineColor(4);
  c8->cd();
  f4->Draw("same");
  f4u->Draw("same");
  f4d->Draw("same");
  //IMnpip_fmix->SetLineColor(4);
  //IMnpip_fmix->Draw("same");
  */
}

Double_t fit_IMnpip(Double_t *x,Double_t *par)
{
  if(!IMnpip_woK0_woSm_mix){
    std::cout << "cannot find template !! " << std::endl;
    return 0;
  }
  Double_t IMnpip=x[0]; 
  Int_t bin = IMnpip_woK0_woSm_mix->GetXaxis()->FindBin(IMnpip);
  Double_t br= par[0]*(IMnpip_woK0_woSm_mix->GetBinContent(bin));
   
  return br;
}

