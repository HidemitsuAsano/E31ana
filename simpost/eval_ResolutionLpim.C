//eval. Resolution of costheta of missp in K-d->Lpimp sim.
void eval_misscos_Lpim()
{

  TFile *f = TFile::Open("../simpost/simIMLpim_ppimpim_v21_out.root");
  TCanvas *c1 = new TCanvas("c1","c1",1200,800);
  TH2D* diffpcos_pcos = (TH2D*)f->Get("diffpcos_pcos");
  diffpcos_pcos->Draw("colz");
  //diffpcos_pcos->RebinY(2);
  TProfile *diffcos_pcos_pfx = (TProfile*)diffpcos_pcos->ProfileX();
  diffcos_pcos_pfx->SetLineColor(2);
  diffcos_pcos_pfx->Draw("same");

  const int nbinsx = diffpcos_pcos->GetNbinsX();
  TH1D* px[1000];
  double pcos[1000];
  double rms[1000];
  double gaussigma[1000];
  double gaussigmaerr[1000];
  double rmserror[1000];
  for(int i=0;i<nbinsx;i++){
    px[i] = (TH1D*)diffpcos_pcos->ProjectionY(Form("px%d",i),i+1,i+2);
    if(px[i]->GetEntries()>100){
      px[i]->Fit("gaus","q","",-0.005,0.005);
      gaussigma[i] = px[i]->GetFunction("gaus")->GetParameter(2);
      gaussigmaerr[i] = px[i]->GetFunction("gaus")->GetParError(2);
      rms[i] = px[i]->GetRMS();
      rmserror[i] = px[i]->GetRMSError();
      pcos[i] = diffpcos_pcos->GetXaxis()->GetBinCenter(i+1);
    }
  }
  
  TCanvas *c2 = new TCanvas("c2","c2");
  TGraphErrors *gr = new TGraphErrors(1000,pcos,gaussigma,0,gaussigmaerr);
  gr->SetMarkerColor(2);
  gr->SetMarkerStyle(20);
  gr->GetYaxis()->SetTitle("Resolution");
  gr->GetXaxis()->SetTitle("miss. p cos#theta");
  gr->GetXaxis()->CenterTitle();
  gr->GetYaxis()->CenterTitle();
  gr->Draw("AP");

}
