void plotCS(const char *filename= "CS_20180622.root"){
  
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetStatBorderSize(1);
  TFile *f = new TFile(filename,"READ");
  f->cd(); 
  TGraphAsymmErrors *grpimSp = f->Get("pimSp_CS_allErr");
  TGraphAsymmErrors *grpipSm = f->Get("pipSm_CS_allErr");
  
  double *x1,*y1,*ey1h,*ey1l;
  double *x2,*y2,*ey2h,*ey2l;
  x1 = grpimSp->GetX();
  y1 = grpimSp->GetY();
  ey1h = grpimSp->GetEYhigh();
  ey1l = grpimSp->GetEYlow();
  x2 = grpipSm->GetX();
  y2 = grpipSm->GetY();
  ey2h = grpipSm->GetEYhigh();
  ey2l = grpipSm->GetEYlow();
  
  int n = grpimSp->GetN();
  TCanvas *cg = new TCanvas("cg","cg",800,800);
  cg->cd();
  TGraphAsymmErrors *grI0I1 = new TGraphAsymmErrors(n);
  for(int i =0;i< n;i++){
    cout << x1[i] << "  " << x2[i] << endl;
    double xav = (x1[i]+x2[i])/2.;
    cout << y1[i] << "  " << y2[i] << endl;
    double yav = (y1[i]+y2[i])/2.;
    //cout << ex1[i] << endl;
    //double exav = ex1[i];
    double eyavh = sqrt(ey1h[i]*ey1h[i]+ey2h[i]*ey2h[i])/2.0;
    double eyavl = sqrt(ey1l[i]*ey1l[i]+ey2l[i]*ey2l[i])/2.0;
    grI0I1->SetPoint(i,xav,yav);
    grI0I1->SetPointError(i,0,0,eyavl,eyavl);
  }
  grI0I1->SetMarkerStyle(20);
  grI0I1->SetMarkerColor(3);
  grI0I1->SetLineWidth(2);
  grI0I1->SetLineColor(3);
  gStyle->SetErrorX(5);
  gStyle->SetEndErrorSize(2);
  gStyle->SetLineWidth(2);
  grI0I1->SetTitle();
  grI0I1->GetXaxis()->SetTitle("#Sigma #pi mass [GeV/c^{2}]");
  grI0I1->GetYaxis()->SetTitle("#frac{d^{2}#sigma}{d#OmegadM} [#mub/sr MeV/c^{2}]");
  grI0I1->GetXaxis()->SetTitleSize(0.04);
  grI0I1->GetXaxis()->SetTitleOffset(1.2);
  grI0I1->GetYaxis()->SetTitleSize(0.04);
  grI0I1->GetYaxis()->SetTitleOffset(1.5);
  grI0I1->GetXaxis()->CenterTitle();
  grI0I1->GetYaxis()->CenterTitle();
  grI0I1->GetXaxis()->SetRangeUser(1.35,1.5);
  grI0I1->GetYaxis()->SetRangeUser(0,25);
  grI0I1->GetXaxis()->SetLabelOffset(0.001);
  grI0I1->GetYaxis()->SetLabelOffset(0.004);
  grI0I1->GetXaxis()->SetLabelSize(0.04);
  grI0I1->GetYaxis()->SetLabelSize(0.04);
  grI0I1->Draw("Ae1P");
  gPad->SetLeftMargin(0.14);
  gPad->SetBottomMargin(0.14);
  //gPad->SetRightMargin(0.03);
  //grI0I1->GetXaxis()->SetLabelSize(0.1,"x");
  //grI0I1->GetYaxis()->SetLabelSize(0.1,"y");
  //gStyle->SetLabelSize(0.1,"y");
  //gStyle->SetLabelOffset(0.03,"x");
  //gStyle->SetLabelOffset(0.03,"y");
  
  gStyle->SetTitleX(0.1f);
  gStyle->SetTitleW(0.8f);
  gStyle->SetTitleBorderSize(0);
  //gPad->SetLeftMargin(0.13);
  //gPad->SetRightMargin(0.03);

  TGraphAsymmErrors *pimS0_CS_allErr = f->Get("pimS0_CS_allErr");
  double *x3,*y3,*ey3l,*ey3h;
  int n2 = pimS0_CS_allErr->GetN();
  TGraphAsymmErrors *grI1 = new TGraphAsymmErrors(n2);
  x3= pimS0_CS_allErr->GetX();
  y3= pimS0_CS_allErr->GetY();
  ey3l = pimS0_CS_allErr->GetEYlow();
  ey3h = pimS0_CS_allErr->GetEYhigh();
  for(int i=0;i<n2;i++){
    grI1->SetPoint(i,x3[i],y3[i]/2.);
    grI1->SetPointError(i,0,0,ey3l[i]/2.,ey3h[i]/2.);
  }
  grI1->SetMarkerStyle(20);
  grI1->SetLineWidth(2);
  grI1->Draw("Pe1");
  TLegend *leg = new TLegend(0.15,0.6,0.5,0.8);
  leg->SetBorderSize(0);
  leg->AddEntry(grI0I1,"1/2(#Sigma^{+}#pi^{-}+#Sigma^{-}#pi^{+})","l");
  leg->AddEntry(grI1,"1/2(#Sigma^{0}#pi^{-})","l");
  leg->SetEntrySeparation(0.1);
  leg->Draw();

}
