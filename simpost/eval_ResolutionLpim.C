//eval. Resolution of costheta of missp in K-d->Lpimp sim.
void eval_ResolutionLpim()
{

  TFile *f = TFile::Open("../simpost/simIMLpim_ppimpim_v28_out.root");
  TCanvas *c1 = new TCanvas("c1","c1",1200,800);
  TH2D* diffpcos_pcos = (TH2D*)f->Get("diffpcos_pcos");
  diffpcos_pcos->Draw("colz");
  //diffpcos_pcos->RebinY(2);
  TProfile *diffcos_pcos_pfx = (TProfile*)diffpcos_pcos->ProfileX();
  diffcos_pcos_pfx->SetLineColor(2);
  diffcos_pcos_pfx->Draw("same");
  
  //costheta resolution
  const int nbinsx = diffpcos_pcos->GetNbinsX();
  TH1D* cospx[1000];
  double pcos[1000];
  double rms[1000];
  double gaussigma[1000];
  double gaussigmaerr[1000];
  double rmserror[1000];
  for(int i=0;i<nbinsx;i++){
    cospx[i] = (TH1D*)diffpcos_pcos->ProjectionY(Form("cospx%d",i),i+1,i+2);
    if(cospx[i]->GetEntries()>100){
      cospx[i]->Fit("gaus","q","",-0.005,0.005);
      gaussigma[i] = cospx[i]->GetFunction("gaus")->GetParameter(2);
      gaussigmaerr[i] = cospx[i]->GetFunction("gaus")->GetParError(2);
      rms[i] = cospx[i]->GetRMS();
      rmserror[i] = cospx[i]->GetRMSError();
      pcos[i] = diffpcos_pcos->GetXaxis()->GetBinCenter(i+1);
    }
  }
  
  TCanvas *c2 = new TCanvas("c2","c2");
  TGraphErrors *gr = new TGraphErrors(1000,pcos,gaussigma,0,gaussigmaerr);
  gr->SetName("pcosreso");
  gr->SetTitle("pmisscos resolution");
  gr->SetMarkerColor(2);
  gr->SetMarkerStyle(20);
  gr->GetYaxis()->SetTitle("Resolution");
  gr->GetXaxis()->SetTitle("miss. p cos#theta");
  gr->GetXaxis()->CenterTitle();
  gr->GetYaxis()->CenterTitle();
  gr->Draw("AP");

  //IM(Lpim) resolution
  TCanvas *c3 = new TCanvas("c3","c3");
  TH2D* diffIMppipi_IMppipi = (TH2D*)f->Get("diffIMppipi_IMppipi");
  
  TH1D* pxM[1000];
  double M[1000];
  double rmsM[1000];
  double gaussmeanM[1000];
  double gaussmeanerrM[1000];
  double gausssigmaM[1000];
  double gausssigmaerrM[1000];
  double rmserrorM[1000];
  const int nbinsxM = diffIMppipi_IMppipi->GetNbinsX();
  for(int i=0;i<nbinsxM;i++){
    pxM[i] = (TH1D*)diffIMppipi_IMppipi->ProjectionY(Form("pxM%d",i),i+1,i+2);
    if(pxM[i]->Integral()>100){
      pxM[i]->Fit("gaus","q","",-0.1,0.1);
      if(pxM[i]->GetFunction("gaus")){
        gaussmeanM[i] = pxM[i]->GetFunction("gaus")->GetParameter(1);
        gaussmeanerrM[i] = pxM[i]->GetFunction("gaus")->GetParError(1);
        gausssigmaM[i] = pxM[i]->GetFunction("gaus")->GetParameter(2);
        gausssigmaerrM[i] = pxM[i]->GetFunction("gaus")->GetParError(2);
        M[i] = diffIMppipi_IMppipi->GetXaxis()->GetBinCenter(i+1);
      }
    }
  }
  
  diffIMppipi_IMppipi->Draw("colz");
  TProfile *diffIMppipi_IMppipi_pfx = (TProfile*)diffIMppipi_IMppipi->ProfileX();
  diffIMppipi_IMppipi_pfx->SetLineColor(2);
  diffIMppipi_IMppipi_pfx->Draw("same");

  TCanvas *c3shift = new TCanvas("c3shift","c3shift");
  TGraphErrors *grMshift = new TGraphErrors(nbinsxM,M,gaussmeanM,0,gausssigmaerrM);
  grMshift->SetName("Mshift");
  grMshift->SetTitle("IM #Lambda#pi^{-} shift");
  grMshift->SetMarkerColor(2);
  grMshift->SetMarkerStyle(20);
  grMshift->GetYaxis()->SetTitle("mass shift [GeV/c^{2}]");
  grMshift->GetXaxis()->SetTitle("IM (#Lambda#pi^{-}) [GeV/c^{2}]");
  grMshift->GetXaxis()->SetRangeUser(1.2,2.0);
  grMshift->GetXaxis()->CenterTitle();
  grMshift->GetYaxis()->CenterTitle();
  grMshift->Draw("AP");

  TCanvas *c3reso = new TCanvas("c3reso","c3reso");
  TGraphErrors *grM = new TGraphErrors(nbinsxM,M,gausssigmaM,0,gausssigmaerrM);
  grM->SetName("Mreso");
  grM->SetTitle("IM #Lambda#pi^{-} resolution");
  grM->SetMarkerColor(2);
  grM->SetMarkerStyle(20);
  grM->GetYaxis()->SetTitle("Resolution [GeV/c^{2}]");
  grM->GetXaxis()->SetTitle("IM (#Lambda#pi^{-}) [GeV/c^{2}]");
  grM->GetXaxis()->SetRangeUser(1.2,2.0);
  grM->GetXaxis()->CenterTitle();
  grM->GetYaxis()->CenterTitle();
  grM->Draw("AP");
  
  TCanvas *c4 = new TCanvas("c4","c4");
  TH2D* diffq_q = (TH2D*)f->Get("diffq_q");
  diffq_q->Draw("colz");
  TProfile *diff_q_pfx = (TProfile*)diffq_q->ProfileX();
  diff_q_pfx->SetLineColor(2);
  diff_q_pfx->Draw("same");


  TH1D* pxq[1000];
  double q[1000];
  double rmsq[1000];
  double gaussmeanq[1000];
  double gaussmeanerrq[1000];
  double gausssigmaq[1000];
  double gausssigmaerrq[1000];
  double rmserrorq[1000];
  const int nbinsxq = diffq_q->GetNbinsX();
  for(int i=0;i<nbinsxq;i++){
    pxq[i] = (TH1D*)diffq_q->ProjectionY(Form("pxq%d",i),i+1,i+2);
    pxq[i]->RebinX(2);
    if(pxq[i]->Integral()>100){
      pxq[i]->Fit("gaus","q","",-0.1,0.1);
      if(pxq[i]->GetFunction("gaus")){
        gaussmeanq[i] = pxq[i]->GetFunction("gaus")->GetParameter(1);
        gaussmeanerrq[i] = pxq[i]->GetFunction("gaus")->GetParError(1);
        gausssigmaq[i] = pxq[i]->GetFunction("gaus")->GetParameter(2);
        gausssigmaerrq[i] = pxq[i]->GetFunction("gaus")->GetParError(2);
        q[i] = diffq_q->GetXaxis()->GetBinCenter(i+1);
        rmsq[i] = pxq[i]->GetRMS();
        rmserrorq[i] = pxq[i]->GetRMSError();
      }
    }
  }
  
  TCanvas *c4shift = new TCanvas("c4shift","c4shift");
  TGraphErrors *grqshift = new TGraphErrors(nbinsxq,q,gaussmeanq,0,gausssigmaerrq);
  grqshift->SetName("qshift");
  grqshift->SetTitle("Mom. Transfer shift");
  grqshift->SetMarkerColor(2);
  grqshift->SetMarkerStyle(20);
  grqshift->GetYaxis()->SetTitle("Mom. Transfer shift [GeV/c]");
  grqshift->GetXaxis()->SetTitle("Mom. Transfer [GeV/c]");
 // grqshift->GetXaxis()->SetRangeUser(1.2,2.0);
  grqshift->GetXaxis()->CenterTitle();
  grqshift->GetYaxis()->CenterTitle();
  grqshift->Draw("AP");
  grqshift->Print();

  TCanvas *c4reso = new TCanvas("c4reso","c4reso");
  TGraphErrors *grq = new TGraphErrors(nbinsxq,q,gausssigmaq,0,gausssigmaerrq);
  //TGraphErrors *grq = new TGraphErrors(nbinsxq,q,rmsq,0,rmserrorq);
  grq->SetName("qreso");
  grq->SetTitle("Mom. Transfer resolution");
  grq->SetMarkerColor(2);
  grq->SetMarkerStyle(20);
  grq->GetYaxis()->SetTitle("Resolution [GeV/c]");
  grq->GetXaxis()->SetTitle("Mom. Trasfer [GeV/c]");
//  grq->GetXaxis()->SetRangeUser(1.2,2.0);
  grq->GetXaxis()->CenterTitle();
  grq->GetYaxis()->CenterTitle();
  grq->Draw("AP");
  grq->Print();
}
