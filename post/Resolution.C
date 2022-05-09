void Resolution()
{
  //const char* filename="../simpost/simIMpisigma_nSppim_pippimn_v156_out_dE2_iso_rej_nostop.root";
  const int Version = 156;
  const char* filenameSp=Form("../simpost/simIMpisigma_nSppim_pippimn_v%d_out_dE2_iso_rej_nostop.root",Version);
  const char* filenameSm=Form("../simpost/simIMpisigma_nSmpip_pippimn_v%d_out_dE2_iso_rej_nostop.root",Version);
  const char* filenameSptheta15=Form("../simpost/simIMpisigma_nSppim_pippimn_v%d_out_dE2_iso_theta15_rej_nostop.root",Version);
  const char* filenameSmtheta15=Form("../simpost/simIMpisigma_nSmpip_pippimn_v%d_out_dE2_iso_theta15_rej_nostop.root",Version);
  std::cout << "infile " << filenameSp <<std::endl;
  std::cout << "infile " << filenameSm <<std::endl;
  //TString pdfname = std::string(filename);
  TString pdfname = std::string(Form("resolution_v%d",Version));
  
  TFile *fSp = new TFile(filenameSp);
  TFile *fSm = new TFile(filenameSm);
  TFile *fSptheta15 = new TFile(filenameSptheta15);
  TFile *fSmtheta15 = new TFile(filenameSmtheta15);
  
  const int nbinIMnpipi = 180;//1.2-2 GeV/c^2
  const int nbinq = 30;// 25;//0-1.5 GeV/c

  //Sp mode caluculation
  TCanvas *cq_IMnpipi_wSid_n_Sp = new TCanvas("cq_IMnpipi_wSid_n_Sp","cq_IMnpipi_wSid_n_Sp",800,800);
  cq_IMnpipi_wSid_n_Sp->cd();
  TH2D* q_IMnpipi_wSid_n_Sp = (TH2D*)fSp->Get("q_IMnpipi_wSid_n_Sp");
  q_IMnpipi_wSid_n_Sp->Draw("colz");

  TCanvas *cdiff_IMnpipi_wSid_n_Sp = new TCanvas("cdiff_IMnpipi_wSid_n_Sp","diff_IMnpipi_wSid_n_Sp",1000,800);
  cdiff_IMnpipi_wSid_n_Sp->cd();
  TH2D* diff_IMnpipi_wSid_n_Sp = (TH2D*)fSp->Get("diff_IMnpipi_wSid_n_Sp");
  diff_IMnpipi_wSid_n_Sp->GetYaxis()->SetRangeUser(-0.1,0.1);
  diff_IMnpipi_wSid_n_Sp->Draw("colz");
  //s: errors are standard deviation 
  TProfile *pfxSp = (TProfile*)diff_IMnpipi_wSid_n_Sp->ProfileX("pfxSp",1,-1,"s");
  pfxSp->SetLineColor(2);
  pfxSp->SetMarkerStyle(33);
  pfxSp->Draw("same");
  TCanvas *cpfxSp = new TCanvas("cpfxSp","cpfxSp",1000,800);
  pfxSp->Draw();
    
  TCanvas *cgaus = new TCanvas("cgaus","cgaus");
  double recomass[nbinIMnpipi];
  double cent[nbinIMnpipi];
  double cent_err[nbinIMnpipi];
  double sigma[nbinIMnpipi];
  double sigma_err[nbinIMnpipi];
  //cgaus->Divide(10,10);
  TH1D *px[nbinIMnpipi];
  for(int i=0; i<nbinIMnpipi; i++) {
    px[i] = (TH1D*)diff_IMnpipi_wSid_n_Sp->ProjectionY(Form("px%d",i),i+1,i+2,"");
    recomass[i] = diff_IMnpipi_wSid_n_Sp->GetXaxis()->GetBinCenter(i+1);
    cgaus->cd(i+1);
    if(px[i]->GetEntries()>200) {
      //px[i]->Draw("HE");
      px[i]->Fit("gaus","q");
      cent[i] = px[i]->GetFunction("gaus")->GetParameter(1);
      cent_err[i] = px[i]->GetFunction("gaus")->GetParError(1);
      sigma[i]= px[i]->GetFunction("gaus")->GetParameter(2);
      sigma_err[i]= px[i]->GetFunction("gaus")->GetParError(2);
    } else {
      cent[i]=-999.;
      cent_err[i]=-0.;
      sigma[i]=-999.;
      sigma_err[i]=0.;
    }
  }
  TCanvas *cfitmean = new TCanvas("cfitmean","fitmean",1000,800);
  cfitmean->cd();
  TGraphErrors *grcent = new TGraphErrors(nbinIMnpipi,recomass,cent,0,cent_err);
  grcent->SetTitle("gaussian mean #pi^{-}#Sigma^{+}");
  grcent->GetXaxis()->SetTitle("IM(#pi^{-}#Sigma^{+}) [GeV/c^{2}]");
  grcent->GetXaxis()->CenterTitle();
  grcent->GetYaxis()->SetTitle("mass center [GeV/c^{2}]");
  grcent->GetYaxis()->SetRangeUser(-0.01,0.03);
  grcent->GetYaxis()->CenterTitle();
  grcent->SetMarkerStyle(20);
  grcent->GetYaxis()->SetMaxDigits(2);
  //grcent->GetYaxis()->SetTitleOffset(1.5);
  grcent->Draw("APE");
  std::cout << __LINE__ << std::endl;
  TCanvas *cfitsigma = new TCanvas("cfitsigma","fitsigma",1000,800);
  cfitsigma->cd();
  TGraphErrors *grsigma = new TGraphErrors(nbinIMnpipi,recomass,sigma,0,sigma_err);
  grsigma->SetTitle("mass resolution");
  grsigma->GetXaxis()->SetTitle("IM(#pi^{-}#Sigma^{+}) [GeV/c^{2}]");
  grsigma->GetXaxis()->CenterTitle();
  grsigma->GetYaxis()->SetTitle("mass resolution [GeV/c^{2}]");
  grsigma->GetYaxis()->CenterTitle();
  //grsigma->GetYaxis()->SetTitleOffset(1.5);
  grsigma->GetYaxis()->SetMaxDigits(2);
  grsigma->GetYaxis()->SetRangeUser(0.0,0.03);
  grsigma->SetMarkerStyle(20);
  grsigma->Draw("APE");


  TCanvas *cdiff_q_wSid_n_Sp = new TCanvas("cdiff_q_wSid_n_Sp","diff_q_wSid_n_Sp",1000,800);
  cdiff_q_wSid_n_Sp->cd();
  TH2D* diff_q_wSid_n_Sp = (TH2D*)fSp->Get("diff_q_wSid_n_Sp");
  diff_q_wSid_n_Sp->GetYaxis()->SetRangeUser(-0.1,0.1);
  diff_q_wSid_n_Sp->Draw("colz");
  TProfile *pfxSp_q = (TProfile*)diff_q_wSid_n_Sp->ProfileX("pfxSp_q",1,-1,"s");
  pfxSp_q->SetLineColor(2);
  pfxSp_q->SetMarkerStyle(33);
  pfxSp_q->Draw("same");
  TCanvas *cgaus_q = new TCanvas("cgaus_q","cgaus_q",1000,800);
  double recoq[nbinq];
  double cent_q[nbinq];
  double cent_q_err[nbinq];
  double sigma_q[nbinq];
  double sigma_q_err[nbinq];
  //cgaus->Divide(10,10);
  TH1D *px_q[nbinq];
  for(int i=0; i<nbinq; i++) {
    px_q[i] = (TH1D*)diff_q_wSid_n_Sp->ProjectionY(Form("px_q%d",i),i+1,i+2,"");
    recoq[i] = diff_q_wSid_n_Sp->GetXaxis()->GetBinCenter(i+1);
    cgaus_q->cd();
    if(px_q[i]->GetEntries()>200) {
      //px[i]->Draw("HE");
      px_q[i]->Fit("gaus","q");
      cent_q[i] = px_q[i]->GetFunction("gaus")->GetParameter(1);
      cent_q_err[i] = px_q[i]->GetFunction("gaus")->GetParError(1);
      sigma_q[i]= px_q[i]->GetFunction("gaus")->GetParameter(2);
      sigma_q_err[i]= px_q[i]->GetFunction("gaus")->GetParError(2);
    } else {
      cent_q[i]=-999.;
      cent_q_err[i]=-0.;
      sigma_q[i]=-999.;
      sigma_q_err[i]=0.;
    }
  }

  TCanvas *cfitmean_q = new TCanvas("cfitmean_q","fitmean_q",1000,800);
  cfitmean_q->cd();
  TGraphErrors *grcent_q = new TGraphErrors(nbinq,recoq,cent_q,0,cent_q_err);
  grcent_q->SetTitle("gaussian mean #pi^{-}#Sigma^{+}");
  grcent_q->GetXaxis()->SetTitle("Mom. Transfer [GeV/c]");
  grcent_q->GetXaxis()->CenterTitle();
  grcent_q->GetYaxis()->SetTitle("Mom. center [GeV/c]");
  grcent_q->GetYaxis()->SetRangeUser(-0.01,0.03);
  grcent_q->GetYaxis()->CenterTitle();
  grcent_q->SetMarkerStyle(20);
  grcent_q->GetYaxis()->SetMaxDigits(2);
  grcent_q->Draw("APE");
  TCanvas *cfitsigma_q = new TCanvas("cfitsigma_q","fitsigma_q",1000,800);
  cfitsigma_q->cd();
  TGraphErrors *grsigma_q = new TGraphErrors(nbinq,recoq,sigma_q,0,sigma_q_err);
  grsigma_q->SetTitle("Mom. resolution #pi^{-}#Sigma^{+} ");
  grsigma_q->GetXaxis()->SetTitle("Mom. Transfer [GeV/c]");
  grsigma_q->GetXaxis()->CenterTitle();
  grsigma_q->GetYaxis()->SetTitle("Mom. resolution [GeV/c]");
  grsigma_q->GetYaxis()->CenterTitle();
  //grsigma_q->GetYaxis()->SetTitleOffset(1.5);
  grsigma_q->GetYaxis()->SetMaxDigits(2);
  grsigma_q->GetYaxis()->SetRangeUser(0.0,0.05);
  grsigma_q->SetMarkerStyle(20);
  grsigma_q->Draw("APE");

  TCanvas *ccos = new TCanvas("ccos","ccos",1200,800);
  //TH2D* diffncos_ncos = (TH2D*)f->Get("diffncos_ncos");
  //diffncos_ncos->Draw("colz");
  TH2D* diffncos_IMnpipi = (TH2D*)fSptheta15->Get("diffncos_IMnpipi");
  diffncos_IMnpipi->GetYaxis()->SetRangeUser(-0.01,0.02);
  diffncos_IMnpipi->GetYaxis()->SetMaxDigits(2);
  diffncos_IMnpipi->Draw("colz");
  //diffncos_ncos->RebinY(2);
  TProfile *diffcos_IMnpipi_pfx = (TProfile*)diffncos_IMnpipi->ProfileX("pfxSpcos",1,-1,"s");
  diffcos_IMnpipi_pfx->SetLineColor(2);
  diffcos_IMnpipi_pfx->Draw("same");

  const int nbinsx = diffncos_IMnpipi->GetNbinsX();
  std::cout << nbinsx << std::endl;
  TH1D* cospx[1000];
  double ncos[1000];
  double rms[1000];
  double gaussigma[1000];
  double gaussigmaerr[1000];
  double rmserror[1000];
  for(int i=0;i<nbinsx;i++){
    cospx[i] = (TH1D*)diffncos_IMnpipi->ProjectionY(Form("cospx%d",i),i+1,i+2);
    ncos[i] = diffncos_IMnpipi->GetXaxis()->GetBinCenter(i+1);
    if(cospx[i]->GetEntries()>10){
      //cospx[i]->Fit("gaus","q","",-0.005,0.005);
      //gaussigma[i] = cospx[i]->GetFunction("gaus")->GetParameter(2);
      //gaussigmaerr[i] = cospx[i]->GetFunction("gaus")->GetParError(2);
      rms[i] = cospx[i]->GetRMS();
      rmserror[i] = cospx[i]->GetRMSError();
    }else{
      rms[i] = 0.0;
      rmserror[i] = 0.0;
    }
  }

  TCanvas *ccosres = new TCanvas("ccosres","ccosres");
  //TGraphErrors *gr = new TGraphErrors(1000,pcos,gaussigma,0,gaussigmaerr);
  TGraphErrors *gr = new TGraphErrors(180,ncos,rms,0,rmserror);
  gr->SetName("ncosreso");
  gr->SetTitle("nmisscos resolution #pi^{-}#Sigma^{+}");
  gr->SetMarkerColor(2);
  gr->SetMarkerStyle(20);
  gr->GetYaxis()->SetTitle("cos#theta_n Resolution");
  gr->GetXaxis()->SetTitle("miss. n cos#theta_{lab}");
  gr->GetXaxis()->CenterTitle();
  gr->GetYaxis()->CenterTitle();
  gr->GetYaxis()->SetMaxDigits(2);
  //gr->Print();
  gr->Draw("AP");

  TCanvas *cq_IMnpipi_wSid_n_Sm = new TCanvas("cq_IMnpipi_wSid_n_Sm","cq_IMnpipi_wSid_n_Sm",800,800);
  cq_IMnpipi_wSid_n_Sm->cd();
  TH2D* q_IMnpipi_wSid_n_Sm = (TH2D*)fSm->Get("q_IMnpipi_wSid_n_Sm");
  q_IMnpipi_wSid_n_Sm->Draw("colz");
    
  TCanvas *cdiff_IMnpipi_wSid_n_Sm = new TCanvas("cdiff_IMnpipi_wSid_n_Sm","diff_IMnpipi_wSid_n_Sm",1000,800);
  cdiff_IMnpipi_wSid_n_Sm->cd();
  TH2D* diff_IMnpipi_wSid_n_Sm = (TH2D*)fSm->Get("diff_IMnpipi_wSid_n_Sm");
  diff_IMnpipi_wSid_n_Sm->GetYaxis()->SetRangeUser(-0.1,0.1);
  diff_IMnpipi_wSid_n_Sm->Draw("colz");
  TProfile *pfxSm = diff_IMnpipi_wSid_n_Sm->ProfileX("pfxSm",1,-1,"s");
  pfxSm->SetLineColor(2);
  pfxSm->SetMarkerStyle(33);
  pfxSm->Draw("same");
    
  TCanvas *cgausSm = new TCanvas("cgausSm","cgausSm");
  double recomassSm[nbinIMnpipi];
  double centSm[nbinIMnpipi];
  double cent_errSm[nbinIMnpipi];
  double sigmaSm[nbinIMnpipi];
  double sigma_errSm[nbinIMnpipi];
  TH1D *pxSm[nbinIMnpipi];
  for(int i=0; i<nbinIMnpipi; i++) {
    pxSm[i] = (TH1D*)diff_IMnpipi_wSid_n_Sm->ProjectionY(Form("pxSm%d",i),i+1,i+2,"");
    recomassSm[i] = diff_IMnpipi_wSid_n_Sm->GetXaxis()->GetBinCenter(i+1);
    //cgaus->cd(i+1);
    if(pxSm[i]->GetEntries()>200) {
      //px[i]->Draw("HE");
      pxSm[i]->Fit("gaus","q");
      centSm[i] = pxSm[i]->GetFunction("gaus")->GetParameter(1);
      cent_errSm[i] = pxSm[i]->GetFunction("gaus")->GetParError(1);
      sigmaSm[i]= pxSm[i]->GetFunction("gaus")->GetParameter(2);
      sigma_errSm[i]= pxSm[i]->GetFunction("gaus")->GetParError(2);
    } else {
      centSm[i]=-999.;
      cent_errSm[i]=0;
      sigmaSm[i]=-999.;
      sigma_errSm[i]=0;
    }
  }

  TCanvas *cfitmeanSm = new TCanvas("cfitmeanSm","fitmeanSm",1000,800);
  cfitmeanSm->cd();
  TGraphErrors *grcentSm = new TGraphErrors(nbinIMnpipi,recomassSm,centSm,0,cent_errSm);
  grcentSm->SetTitle("gaussian mean #pi^{+}#Sigma^{-}");
  grcentSm->GetXaxis()->SetTitle("IM(#pi^{+}#Sigma^{-}) [GeV/c^{2}]");
  grcentSm->GetXaxis()->CenterTitle();
  grcentSm->GetYaxis()->SetTitle("mass center [GeV/c^{2}]");
  grcentSm->GetYaxis()->SetRangeUser(-0.02,0.03);
  grcentSm->GetYaxis()->CenterTitle();
  grcentSm->SetMarkerStyle(20);
  //grcentSm->GetYaxis()->SetTitleOffset(1.5);
  grcentSm->GetYaxis()->SetMaxDigits(2);
  grcentSm->Draw("APE");

  TCanvas *cfitsigmaSm = new TCanvas("cfitsigmaSm","fitsigmaSm",1000,800);
  cfitsigmaSm->cd();
  TGraphErrors *grsigmaSm = new TGraphErrors(nbinIMnpipi,recomassSm,sigmaSm,0,sigma_errSm);
  grsigmaSm->SetTitle("mass resolution #pi^{+}#Sigma^{-}");
  grsigmaSm->GetXaxis()->SetTitle("IM(#pi^{+}#Sigma^{-}) [GeV/c^{2}]");
  grsigmaSm->GetXaxis()->CenterTitle();
  grsigmaSm->GetYaxis()->SetTitle("mass resolution [GeV/c^{2}]");
  grsigmaSm->GetYaxis()->CenterTitle();
  //grsigmaSm->GetYaxis()->SetTitleOffset(1.5);
  grsigmaSm->GetYaxis()->SetMaxDigits(2);
  grsigmaSm->GetYaxis()->SetRangeUser(0.0,0.03);
  grsigmaSm->SetMarkerStyle(20);
  grsigmaSm->Draw("APE");

  TCanvas *cdiff_q_wSid_n_Sm = new TCanvas("cdiff_q_wSid_n_Sm","diff_q_wSid_n_Sm",1000,800);
  cdiff_q_wSid_n_Sm->cd();
  TH2D* diff_q_wSid_n_Sm = (TH2D*)fSm->Get("diff_q_wSid_n_Sm");
  diff_q_wSid_n_Sm->GetYaxis()->SetRangeUser(-0.1,0.1);
  diff_q_wSid_n_Sm->Draw("colz");
  TProfile *pfxSm_q = (TProfile*)diff_q_wSid_n_Sm->ProfileX("pfxSm_q",1,-1,"s");
  pfxSm_q->SetLineColor(2);
  pfxSm_q->SetMarkerStyle(33);
  pfxSm_q->Draw("same");
  TCanvas *cgaus_qSm = new TCanvas("cgaus_qSm","cgaus_qSm",800,800);
  double recoqSm[nbinq];
  double cent_qSm[nbinq];
  double cent_q_errSm[nbinq];
  double sigma_qSm[nbinq];
  double sigma_q_errSm[nbinq];
  TH1D *px_qSm[nbinq];
  for(int i=0; i<nbinq; i++) {
    px_qSm[i] = (TH1D*)diff_q_wSid_n_Sm->ProjectionY(Form("px_qSm%d",i),i+1,i+2);
    recoqSm[i] = diff_q_wSid_n_Sm->GetXaxis()->GetBinCenter(i+1);
    cgaus_qSm->cd();
    if(px_qSm[i]->GetEntries()>200) {
      px_qSm[i]->Fit("gaus","q");
      cent_qSm[i] = px_qSm[i]->GetFunction("gaus")->GetParameter(1);
      cent_q_errSm[i] = px_qSm[i]->GetFunction("gaus")->GetParError(1);
      sigma_qSm[i]= px_qSm[i]->GetFunction("gaus")->GetParameter(2);
      sigma_q_errSm[i]= px_qSm[i]->GetFunction("gaus")->GetParError(2);
    } else {
      cent_qSm[i]=-999.;
      cent_q_errSm[i]=-0.;
      sigma_qSm[i]=-999.;
      sigma_q_errSm[i]=0.;
    }
  }

  TCanvas *cfitmean_qSm = new TCanvas("cfitmean_qSm","fitmean_qSm",1000,800);
  cfitmean_qSm->cd();
  TGraphErrors *grcent_qSm = new TGraphErrors(nbinq,recoqSm,cent_qSm,0,cent_q_errSm);
  grcent_qSm->SetTitle("gaussian mean #pi^{+}#Sigma^{-}");
  grcent_qSm->GetXaxis()->SetTitle("Mom. Transfer [GeV/c]");
  grcent_qSm->GetXaxis()->CenterTitle();
  grcent_qSm->GetYaxis()->SetTitle("Mom. center [GeV/c]");
  grcent_qSm->GetYaxis()->SetRangeUser(-0.02,0.05);
  grcent_qSm->GetYaxis()->CenterTitle();
  grcent_qSm->SetMarkerStyle(20);
  //grcent_qSm->GetYaxis()->SetTitleOffset(1.5);
  grcent_qSm->GetYaxis()->SetMaxDigits(2);
  grcent_qSm->Draw("APE");
  TCanvas *cfitsigma_qSm = new TCanvas("cfitsigma_qSm","fitsigma_qSm",1000,800);
  cfitsigma_qSm->cd();
  TGraphErrors *grsigma_qSm = new TGraphErrors(nbinq,recoqSm,sigma_qSm,0,sigma_q_errSm);
  grsigma_qSm->SetTitle("Mom. resolution #pi^{+}#Sigma^{-}");
  grsigma_qSm->GetXaxis()->SetTitle("Mom. Transfer [GeV/c]");
  grsigma_qSm->GetXaxis()->CenterTitle();
  grsigma_qSm->GetYaxis()->SetTitle("Mom. resolution [GeV/c]");
  grsigma_qSm->GetYaxis()->CenterTitle();
  //grsigma_qSm->GetYaxis()->SetTitleOffset(1.5);
  grsigma_qSm->GetYaxis()->SetMaxDigits(2);
  grsigma_qSm->GetYaxis()->SetRangeUser(0.0,0.05);
  grsigma_qSm->SetMarkerStyle(20);
  grsigma_qSm->Draw("APE");



}
