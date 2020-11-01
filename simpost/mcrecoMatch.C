void mcrecoMatch(){

  TFile *fileSp = TFile::Open("simIMpisigma_nSppim_pippimn_v108_out.root","READ");
  fileSp->cd();
  
  //mcData -reaction Data
  TCanvas *cdiff_nmiss_reactmc = new TCanvas("cdiff_nmiss_reactmc","cdiff_nmiss_reactmc");
  TH1D* diff_nmiss_reactmc = (TH1D*)fileSp->Get("diff_nmiss_reactmc");
  diff_nmiss_reactmc->Draw("HE");
  cdiff_nmiss_reactmc->SetLogy();

  TCanvas *cdiff_cosnmiss_reactmc = new TCanvas("cdiff_cosnmiss_reactmc","cdiff_cosnmiss_reactmc");
  TH1D* diff_cosnmiss_reactmc = (TH1D*)fileSp->Get("diff_cosnmiss_reactmc");
  diff_cosnmiss_reactmc->Draw("HE");
  cdiff_cosnmiss_reactmc->SetLogy();

  TCanvas *cdiff_IMnpip_reactmc = new TCanvas("cdiff_IMnpip_reactmc","cdiff_IMnpip_reactmc");
  TH1D* diff_IMnpip_reactmc = (TH1D*)fileSp->Get("diff_IMnpip_reactmc");
  diff_IMnpip_reactmc->Draw("HE");
  cdiff_IMnpip_reactmc->SetLogy();

  //reco - mcData 
  TCanvas *cSpdiff2D_MMnmiss_IMnpip = new TCanvas("cSpdiff_MMnmiss_IMnpip","cSpdiff_MMnmiss_IMnpip");
  TH2D* diff2D_MMnmiss_IMnpip_recomc_wSid_n = (TH2D*)fileSp->Get("diff2D_MMnmiss_IMnpip_recomc_wSid_n");
  TH2D* diff2D_MMnmiss_IMnpip_recomc_wSid_n_fake1 = (TH2D*)fileSp->Get("diff2D_MMnmiss_IMnpip_recomc_wSid_n_fake1");
  diff2D_MMnmiss_IMnpip_recomc_wSid_n_fake1->Draw("colz");
  cSpdiff2D_MMnmiss_IMnpip->SetLogz();

  TCanvas *cSpdiff2D_MMnmiss = new TCanvas("cSpdiff_MMnmiss","cSpdiff_MMnmiss");
  cSpdiff2D_MMnmiss->cd();
  diff2D_MMnmiss_IMnpip_recomc_wSid_n->ProjectionY()->Draw("HE");
  TH1D* MMnmiss_fake1= (TH1D*)diff2D_MMnmiss_IMnpip_recomc_wSid_n_fake1->ProjectionY();
  MMnmiss_fake1->SetLineColor(3);
  MMnmiss_fake1->Draw("HEsame");
  cSpdiff2D_MMnmiss->SetLogy();
  
  //mass check is applied
  TCanvas *cSpdiff2D_nmom_IMnpip = new TCanvas("cSpdiff_nmom_IMnpip","cSpdiff_nmom_IMnpip");
  TH2D* diff2D_nmom_IMnpip_recomc_wSid_n = (TH2D*)fileSp->Get("diff2D_nmom_IMnpip_recomc_wSid_n");
  TH2D* diff2D_nmom_IMnpip_recomc_wSid_n_fake1 = (TH2D*)fileSp->Get("diff2D_nmom_IMnpip_recomc_wSid_n_fake1");
  diff2D_nmom_IMnpip_recomc_wSid_n_fake1->Draw("colz");
  cSpdiff2D_nmom_IMnpip->SetLogz();
  
  TCanvas *cSpdiff2D_nmom = new TCanvas("cSpdiff_nmom","cSpdiff_nmom");
  cSpdiff2D_nmom->cd();
  diff2D_nmom_IMnpip_recomc_wSid_n->ProjectionY()->Draw("HE");
  TH1D* nmom_fake1= (TH1D*)diff2D_nmom_IMnpip_recomc_wSid_n_fake1->ProjectionY();
  nmom_fake1->SetLineColor(3);
  nmom_fake1->Draw("HEsame");
  nmom_fake1->Fit("gaus","","",-0.1,0.1);
  cSpdiff2D_nmom->SetLogy();

  TCanvas *cSpdiff2D_IMnpip = new TCanvas("cSpdiff_IMnpip","cSpdiff_IMnpip");
  cSpdiff2D_IMnpip->cd();
  TH1D* IMnpip= (TH1D*)diff2D_nmom_IMnpip_recomc_wSid_n->ProjectionX();
  IMnpip->GetXaxis()->SetRangeUser(-0.15,0.15);
  IMnpip->Draw("HE");
  TH1D* IMnpip_fake1= (TH1D*)diff2D_nmom_IMnpip_recomc_wSid_n_fake1->ProjectionX();
  IMnpip_fake1->SetLineColor(3);
  IMnpip_fake1->Draw("HEsame");
  //IMnpip_fake1->Fit("gaus","","",-0.1,0.1);
  cSpdiff2D_IMnpip->SetLogy();


  //Sm mode
  TFile *fileSm = TFile::Open("simIMpisigma_nSmpip_pippimn_v108_out.root","READ");
  
  //mcData -reaction Data
  TCanvas *cSmdiff_nmiss_reactmc = new TCanvas("cSmdiff_nmiss_reactmc","cSmdiff_nmiss_reactmc");
  TH1D* Smdiff_nmiss_reactmc = (TH1D*)fileSm->Get("diff_nmiss_reactmc");
  Smdiff_nmiss_reactmc->Draw("HE");
  cSmdiff_nmiss_reactmc->SetLogy();
  
  TCanvas *cSmdiff_cosnmiss_reactmc = new TCanvas("cSmdiff_cosnmiss_reactmc","cSmdiff_cosnmiss_reactmc");
  TH1D* Smdiff_cosnmiss_reactmc = (TH1D*)fileSm->Get("diff_cosnmiss_reactmc");
  Smdiff_cosnmiss_reactmc->Draw("HE");
  cSmdiff_cosnmiss_reactmc->SetLogy();
 
  TCanvas *cdiff_IMnpim_reactmc = new TCanvas("cdiff_IMnpim_reactmc","cdiff_IMnpim_reactmc");
  TH1D* diff_IMnpim_reactmc = (TH1D*)fileSm->Get("diff_IMnpim_reactmc");
  diff_IMnpim_reactmc->Draw("HE");
  cdiff_IMnpim_reactmc->SetLogy();


  TCanvas *cSmdiff2D_MMnmiss_IMnpim = new TCanvas("cSmdiff_MMnmiss_IMnpim","cSmdiff_MMnmiss_IMnpim");
  TH2D* diff2D_MMnmiss_IMnpim_recomc_wSid_n = (TH2D*)fileSm->Get("diff2D_MMnmiss_IMnpim_recomc_wSid_n");
  TH2D* diff2D_MMnmiss_IMnpim_recomc_wSid_n_fake1 = (TH2D*)fileSm->Get("diff2D_MMnmiss_IMnpim_recomc_wSid_n_fake1");
  diff2D_MMnmiss_IMnpim_recomc_wSid_n_fake1->Draw("colz");
  cSmdiff2D_MMnmiss_IMnpim->SetLogz();

  TCanvas *cSmdiff2D_MMnmiss = new TCanvas("cSmdiff_MMnmiss","cSmdiff_MMnmiss");
  cSmdiff2D_MMnmiss->cd();
  diff2D_MMnmiss_IMnpim_recomc_wSid_n->ProjectionY()->Draw("HE");
  TH1D* SmMMnmiss_fake1= (TH1D*)diff2D_MMnmiss_IMnpim_recomc_wSid_n_fake1->ProjectionY();
  SmMMnmiss_fake1->SetLineColor(3);
  SmMMnmiss_fake1->Draw("HEsame");
  cSmdiff2D_MMnmiss->SetLogy();
  
  TCanvas *cSmdiff2D_nmom_IMnpim = new TCanvas("cSmdiff_nmom_IMnpim","cSmdiff_nmom_IMnpim");
  TH2D* diff2D_nmom_IMnpim_recomc_wSid_n = (TH2D*)fileSm->Get("diff2D_nmom_IMnpim_recomc_wSid_n");
  TH2D* diff2D_nmom_IMnpim_recomc_wSid_n_fake1 = (TH2D*)fileSm->Get("diff2D_nmom_IMnpim_recomc_wSid_n_fake1");
  diff2D_nmom_IMnpim_recomc_wSid_n_fake1->Draw("colz");
  cSmdiff2D_nmom_IMnpim->SetLogz();


  TCanvas *cSmdiff2D_nmom = new TCanvas("cSmdiff_nmom","cSmdiff_nmom");
  cSmdiff2D_nmom->cd();
  diff2D_nmom_IMnpim_recomc_wSid_n->ProjectionY()->Draw("HE");
  TH1D* Smnmom_fake1= (TH1D*)diff2D_nmom_IMnpim_recomc_wSid_n_fake1->ProjectionY();
  Smnmom_fake1->SetLineColor(3);
  Smnmom_fake1->Draw("HEsame");
  Smnmom_fake1->Fit("gaus","","",-0.1,0.1);
  cSmdiff2D_nmom->SetLogy();


  TCanvas *cSmdiff2D_IMnpim = new TCanvas("cSmdiff_IMnpim","cSmdiff_IMnpim");
  cSmdiff2D_IMnpim->cd();
  TH1D* SmIMnpim = (TH1D*)diff2D_nmom_IMnpim_recomc_wSid_n->ProjectionX();
  SmIMnpim->Draw("HE");
  TH1D* SmIMnpim_fake1= (TH1D*)diff2D_nmom_IMnpim_recomc_wSid_n_fake1->ProjectionX();
  SmIMnpim_fake1->SetLineColor(3);
  SmIMnpim->GetXaxis()->SetRangeUser(-0.15,0.15);
  SmIMnpim_fake1->Draw("HEsame");
  //SmIMnpim_fake1->Fit("gaus","","",-0.1,0.1);
  cSmdiff2D_IMnpim->SetLogy();












};
