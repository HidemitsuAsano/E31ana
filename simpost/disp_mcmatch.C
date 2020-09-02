void disp_mcmatch( )
{
  TFile *fGSp = TFile::Open("simIMpisigma_nSppim_pippimn_v108_out.root","READ");
  fGSp->cd();
  TCanvas *c1 = new TCanvas();
  TH2F* diff2D_IMnpip_Momnpip_wSid_n_reactmc_Sp 
  = (TH2F*)fGSp->Get("diff2D_IMnpip_Momnpip_wSid_n_reactmc");
  diff2D_IMnpip_Momnpip_wSid_n_reactmc_Sp->Draw("colz");
  c1->SetLogz();
  TCanvas *c3 = new TCanvas();
  diff2D_IMnpip_Momnpip_wSid_n_reactmc_Sp->ProjectionX()->Draw("H");
  c3->SetLogy();
  TCanvas *c5 = new TCanvas();
  diff2D_IMnpip_Momnpip_wSid_n_reactmc_Sp->ProjectionY()->Draw("H");
  c5->SetLogy();
  TCanvas *c7 = new TCanvas();
  TH2F* diffMomnpim_Momnpip_recomc_wSid_n_Sp 
  = (TH2F*)fGSp->Get("diffMomnpim_Momnpip_recomc_wSid_n");
  diffMomnpim_Momnpip_recomc_wSid_n_Sp->ProjectionX()->Draw("H");


  

  TFile *fGSm = TFile::Open("simIMpisigma_nSmpip_pippimn_v108_out.root","READ");
  fGSm->cd();
  TCanvas *c2 = new TCanvas();
  TH2F* diff2D_IMnpim_Momnpim_wSid_n_reactmc_Sm 
  = fGSm->Get("diff2D_IMnpim_Momnpim_wSid_n_reactmc");
  diff2D_IMnpim_Momnpim_wSid_n_reactmc_Sm->Draw("colz");
  c2->SetLogz();
  TCanvas *c4 = new TCanvas();
  diff2D_IMnpim_Momnpim_wSid_n_reactmc_Sm->ProjectionX()->Draw("H");
  c4->SetLogy();
  TCanvas *c6 = new TCanvas();
  diff2D_IMnpim_Momnpim_wSid_n_reactmc_Sm->ProjectionY()->Draw("H");
  c6->SetLogy();
  TCanvas *c8 = new TCanvas();
  TH2F* diffMomnpim_Momnpip_recomc_wSid_n_Sm 
  = fGSm->Get("diffMomnpim_Momnpip_recomc_wSid_n");
  diffMomnpim_Momnpip_recomc_wSid_n_Sm->ProjectionY()->Draw("H");











}
