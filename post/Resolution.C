void Resolution(const char* filename="../simpost/simIMpisigma_nSppim_pippimn_v156_out_dE2_iso_rej_nostop.root")
{
  
  std::cout << "infile " << filename <<std::endl;
  TString pdfname = std::string(filename);

  bool SimSpmode = (std::string(filename).find("Sp")!= std::string::npos);
  bool SimSmmode = (std::string(filename).find("Sm")!= std::string::npos);
  
  if(SimSpmode) std::cout << "Sim Sigma+ mode" << std::endl;
  if(SimSmmode) std::cout << "Sim Sigma- mode" << std::endl;

  TFile *f = new TFile(filename);
  
  const int nbinIMnpipi = 80;//1.2-2 GeV/c^2
  const int nbinq = 100;// 25;//0-1.5 GeV/c
  if(SimSpmode) {
    TCanvas *cq_IMnpipi_wSid_n_Sp = new TCanvas("cq_IMnpipi_wSid_n_Sp","cq_IMnpipi_wSid_n_Sp",800,800);
    cq_IMnpipi_wSid_n_Sp->cd();
    TH2D* q_IMnpipi_wSid_n_Sp = (TH2D*)f->Get("q_IMnpipi_wSid_n_Sp");
    q_IMnpipi_wSid_n_Sp->Draw("colz");

    TCanvas *cdiff_IMnpipi_wSid_n_Sp = new TCanvas("cdiff_IMnpipi_wSid_n_Sp","diff_IMnpipi_wSid_n_Sp",1000,800);
    cdiff_IMnpipi_wSid_n_Sp->cd();
    TH2D* diff_IMnpipi_wSid_n_Sp = (TH2D*)f->Get("diff_IMnpipi_wSid_n_Sp");
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
    grcent->SetTitle("gaussian mean");
    grcent->GetXaxis()->SetTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
    grcent->GetXaxis()->CenterTitle();
    grcent->GetYaxis()->SetTitle("mass center [GeV/c^{2}]");
    grcent->GetYaxis()->SetRangeUser(-0.01,0.03);
    grcent->GetYaxis()->CenterTitle();
    grcent->SetMarkerStyle(20);
    grcent->GetYaxis()->SetTitleOffset(1.5);
    grcent->Draw("APE");
    TCanvas *cfitsigma = new TCanvas("cfitsigma","fitsigma",1000,800);
    cfitsigma->cd();
    TGraphErrors *grsigma = new TGraphErrors(nbinIMnpipi,recomass,sigma,0,sigma_err);
    grsigma->SetTitle("mass resolution");
    grsigma->GetXaxis()->SetTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
    grsigma->GetXaxis()->CenterTitle();
    grsigma->GetYaxis()->SetTitle("mass resolution [GeV/c^{2}]");
    grsigma->GetYaxis()->CenterTitle();
    grsigma->GetYaxis()->SetTitleOffset(1.5);
    grsigma->GetYaxis()->SetRangeUser(0.0,0.03);
    grsigma->SetMarkerStyle(20);
    grsigma->Draw("APE");


    TCanvas *cdiff_q_wSid_n_Sp = new TCanvas("cdiff_q_wSid_n_Sp","diff_q_wSid_n_Sp",1000,800);
    cdiff_q_wSid_n_Sp->cd();
    TH2D* diff_q_wSid_n_Sp = (TH2D*)f->Get("diff_q_wSid_n_Sp");
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
    grcent_q->SetTitle("gaussian mean");
    grcent_q->GetXaxis()->SetTitle("Mom. Transfer [GeV/c]");
    grcent_q->GetXaxis()->CenterTitle();
    grcent_q->GetYaxis()->SetTitle("Mom. center [GeV/c]");
    grcent_q->GetYaxis()->SetRangeUser(-0.01,0.05);
    grcent_q->GetYaxis()->CenterTitle();
    grcent_q->SetMarkerStyle(20);
    grcent_q->GetYaxis()->SetTitleOffset(1.5);
    grcent_q->Draw("APE");
    TCanvas *cfitsigma_q = new TCanvas("cfitsigma_q","fitsigma_q",1000,800);
    cfitsigma_q->cd();
    TGraphErrors *grsigma_q = new TGraphErrors(nbinq,recoq,sigma_q,0,sigma_q_err);
    grsigma_q->SetTitle("Mom. resolution");
    grsigma_q->GetXaxis()->SetTitle("Mom. Transfer [GeV/c]");
    grsigma_q->GetXaxis()->CenterTitle();
    grsigma_q->GetYaxis()->SetTitle("Mom. resolution [GeV/c]");
    grsigma_q->GetYaxis()->CenterTitle();
    grsigma_q->GetYaxis()->SetTitleOffset(1.5);
    grsigma_q->GetYaxis()->SetRangeUser(0.0,0.05);
    grsigma_q->SetMarkerStyle(20);
    grsigma_q->Draw("APE");
  
    TCanvas *c1 = new TCanvas("c1","c1",1200,800);
    TH2D* diffncos_ncos = (TH2D*)f->Get("diffncos_ncos");
    diffncos_ncos->Draw("colz");
    diffncos_ncos->RebinY(2);
    TProfile *diffcos_ncos_pfx = (TProfile*)diffncos_ncos->ProfileX();
    diffcos_ncos_pfx->SetLineColor(2);
    diffcos_ncos_pfx->Draw("same");
  
    const int nbinsx = diffncos_ncos->GetNbinsX();
    std::cout << nbinsx << std::endl;
    TH1D* cospx[1000];
    double ncos[1000];
    double rms[1000];
    double gaussigma[1000];
    double gaussigmaerr[1000];
    double rmserror[1000];
    for(int i=0;i<nbinsx;i++){
      cospx[i] = (TH1D*)diffncos_ncos->ProjectionY(Form("cospx%d",i),i+1,i+2);
      ncos[i] = diffncos_ncos->GetXaxis()->GetBinCenter(i+1);
      if(cospx[i]->GetEntries()>10){
        //cospx[i]->Fit("gaus","q","",-0.005,0.005);
        //gaussigma[i] = cospx[i]->GetFunction("gaus")->GetParameter(2);
        //gaussigmaerr[i] = cospx[i]->GetFunction("gaus")->GetParError(2);
        rms[i] = cospx[i]->GetRMS();
        rmserror[i] = cospx[i]->GetRMSError();
      }else{
        rms[i] = 0;
        rmserror[i] = 0;
      }
    }

    TCanvas *c2 = new TCanvas("c2","c2");
    //TGraphErrors *gr = new TGraphErrors(1000,pcos,gaussigma,0,gaussigmaerr);
    TGraphErrors *gr = new TGraphErrors(1000,ncos,rms,0,rmserror);
    gr->SetName("ncosreso");
    gr->SetTitle("nmisscos resolution");
    gr->SetMarkerColor(2);
    gr->SetMarkerStyle(20);
    gr->GetYaxis()->SetTitle("cos#theta_n Resolution");
    gr->GetXaxis()->SetTitle("miss. n cos#theta_{lab}");
    gr->GetXaxis()->CenterTitle();
    gr->GetYaxis()->CenterTitle();
    gr->Print();
    gr->Draw("AP");
  
  
  
  
  }//if Spmode

  if(SimSmmode) {
    TCanvas *cq_IMnpipi_wSid_n_Sm = new TCanvas("cq_IMnpipi_wSid_n_Sm","cq_IMnpipi_wSid_n_Sm",800,800);
    cq_IMnpipi_wSid_n_Sm->cd();
    TH2D* q_IMnpipi_wSid_n_Sm = (TH2D*)f->Get("q_IMnpipi_wSid_n_Sm");
    q_IMnpipi_wSid_n_Sm->Draw("colz");
    
    TCanvas *cdiff_IMnpipi_wSid_n_Sm = new TCanvas("cdiff_IMnpipi_wSid_n_Sm","diff_IMnpipi_wSid_n_Sm",1000,800);
    cdiff_IMnpipi_wSid_n_Sm->cd();
    TH2D* diff_IMnpipi_wSid_n_Sm = (TH2D*)f->Get("diff_IMnpipi_wSid_n_Sm");
    diff_IMnpipi_wSid_n_Sm->GetYaxis()->SetRangeUser(-0.1,0.1);
    diff_IMnpipi_wSid_n_Sm->Draw("colz");
    TProfile *pfxSm = diff_IMnpipi_wSid_n_Sm->ProfileX("pfxSm",1,-1,"s");
    pfxSm->SetLineColor(2);
    pfxSm->SetMarkerStyle(33);
    pfxSm->Draw("same");
    
    TCanvas *cgaus = new TCanvas("cgaus","cgaus");
    double recomass[nbinIMnpipi];
    double cent[nbinIMnpipi];
    double cent_err[nbinIMnpipi];
    double sigma[nbinIMnpipi];
    double sigma_err[nbinIMnpipi];
    //cgaus->Divide(10,10);
    TH1D *px[nbinIMnpipi];
    for(int i=0; i<nbinIMnpipi; i++) {
      px[i] = (TH1D*)diff_IMnpipi_wSid_n_Sm->ProjectionY(Form("px%d",i),i+1,i+2,"");
      recomass[i] = diff_IMnpipi_wSid_n_Sm->GetXaxis()->GetBinCenter(i+1);
      //cgaus->cd(i+1);
      if(px[i]->GetEntries()>200) {
        //px[i]->Draw("HE");
        px[i]->Fit("gaus","q");
        cent[i] = px[i]->GetFunction("gaus")->GetParameter(1);
        cent_err[i] = px[i]->GetFunction("gaus")->GetParError(1);
        sigma[i]= px[i]->GetFunction("gaus")->GetParameter(2);
        sigma_err[i]= px[i]->GetFunction("gaus")->GetParError(2);
      } else {
        cent[i]=-999.;
        cent_err[i]=0;
        sigma[i]=-999.;
        sigma_err[i]=0;
      }
    }
    TCanvas *cfitmean = new TCanvas("cfitmean","fitmean",1000,800);
    cfitmean->cd();
    TGraphErrors *grcent = new TGraphErrors(nbinIMnpipi,recomass,cent,0,cent_err);
    grcent->SetTitle("gaussian mean");
    grcent->GetXaxis()->SetTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
    grcent->GetXaxis()->CenterTitle();
    grcent->GetYaxis()->SetTitle("mass center [GeV/c^{2}]");
    grcent->GetYaxis()->SetRangeUser(-0.02,0.03);
    grcent->GetYaxis()->CenterTitle();
    grcent->SetMarkerStyle(20);
    grcent->GetYaxis()->SetTitleOffset(1.5);
    grcent->Draw("APE");

    TCanvas *cfitsigma = new TCanvas("cfitsigma","fitsigma",1000,800);
    cfitsigma->cd();
    TGraphErrors *grsigma = new TGraphErrors(nbinIMnpipi,recomass,sigma,0,sigma_err);
    grsigma->SetTitle("mass resolution");
    grsigma->GetXaxis()->SetTitle("IM(n#pi^{+}#pi^{-}) [GeV/c^{2}]");
    grsigma->GetXaxis()->CenterTitle();
    grsigma->GetYaxis()->SetTitle("mass resolution [GeV/c^{2}]");
    grsigma->GetYaxis()->CenterTitle();
    grsigma->GetYaxis()->SetTitleOffset(1.5);
    grsigma->GetYaxis()->SetRangeUser(0.0,0.03);
    grsigma->SetMarkerStyle(20);
    grsigma->Draw("APE");

    TCanvas *cdiff_q_wSid_n_Sm = new TCanvas("cdiff_q_wSid_n_Sm","diff_q_wSid_n_Sm",1000,800);
    cdiff_q_wSid_n_Sm->cd();
    TH2D* diff_q_wSid_n_Sm = (TH2D*)f->Get("diff_q_wSid_n_Sm");
    diff_q_wSid_n_Sm->GetYaxis()->SetRangeUser(-0.1,0.1);
    diff_q_wSid_n_Sm->Draw("colz");
    TProfile *pfxSm_q = (TProfile*)diff_q_wSid_n_Sm->ProfileX("pfxSm_q",1,-1,"s");
    pfxSm_q->SetLineColor(2);
    pfxSm_q->SetMarkerStyle(33);
    pfxSm_q->Draw("same");
    TCanvas *cgaus_q = new TCanvas("cgaus_q","cgaus_q",800,800);
    double recoq[nbinq];
    double cent_q[nbinq];
    double cent_q_err[nbinq];
    double sigma_q[nbinq];
    double sigma_q_err[nbinq];
    //cgaus->Divide(10,10);
    TH1D *px_q[nbinq];
    for(int i=0; i<nbinq; i++) {
      px_q[i] = (TH1D*)diff_q_wSid_n_Sm->ProjectionY(Form("px_q%d",i),i+1,i+2,"");
      recoq[i] = diff_q_wSid_n_Sm->GetXaxis()->GetBinCenter(i+1);
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
    grcent_q->SetTitle("gaussian mean");
    grcent_q->GetXaxis()->SetTitle("Mom. Transfer [GeV/c]");
    grcent_q->GetXaxis()->CenterTitle();
    grcent_q->GetYaxis()->SetTitle("Mom. center [GeV/c]");
    grcent_q->GetYaxis()->SetRangeUser(-0.02,0.05);
    grcent_q->GetYaxis()->CenterTitle();
    grcent_q->SetMarkerStyle(20);
    grcent_q->GetYaxis()->SetTitleOffset(1.5);
    grcent_q->Draw("APE");
    TCanvas *cfitsigma_q = new TCanvas("cfitsigma_q","fitsigma_q",1000,800);
    cfitsigma_q->cd();
    TGraphErrors *grsigma_q = new TGraphErrors(nbinq,recoq,sigma_q,0,sigma_q_err);
    grsigma_q->SetTitle("Mom. resolution");
    grsigma_q->GetXaxis()->SetTitle("Mom. Transfer [GeV/c]");
    grsigma_q->GetXaxis()->CenterTitle();
    grsigma_q->GetYaxis()->SetTitle("Mom. resolution [GeV/c]");
    grsigma_q->GetYaxis()->CenterTitle();
    grsigma_q->GetYaxis()->SetTitleOffset(1.5);
    grsigma_q->GetYaxis()->SetRangeUser(0.0,0.05);
    grsigma_q->SetMarkerStyle(20);
    grsigma_q->Draw("APE");
  }//if SimSmmode
  




}
