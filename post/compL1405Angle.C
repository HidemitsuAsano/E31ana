void compL1405Angle()
{
  gStyle->SetOptStat(0);
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPadTopMargin(0.02);
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetTitleYOffset(1.8);
  gStyle->SetLabelOffset(0.0002,"y");

  TFile *ff = TFile::Open("csfinal.root");

  auto *CS_CosS1385Lpim = (TH1D*)ff->Get("CS_CosS1385Lpim");
  auto *grCosL1405 = (TGraphAsymmErrors*)ff->Get("grCosL1405");
  auto *gMIXErrorCosL1405 = (TGraphAsymmErrors*)ff->Get("gMIXErrorCosL1405");
  auto *grCosL1520 = (TGraphAsymmErrors*)ff->Get("grCosL1520");
  auto *gMIXErrorCosL1520 = (TGraphAsymmErrors*)ff->Get("gMIXErrorCosL1520");
  auto *grCosQF = (TGraphAsymmErrors*)ff->Get("grCosQF");
  auto *gMIXErrorCosQF = (TGraphAsymmErrors*)ff->Get("gMIXErrorCosQF");

  TCanvas *c1 = new TCanvas("c1","c1",1200,1000);
  grCosL1405->SetTitle("");
  grCosQF->SetTitle("");
  grCosL1405->GetYaxis()->SetTitleOffset(1.5);
  grCosQF->SetLineColor(3);
  grCosQF->SetMarkerColor(3);
  grCosQF->Draw("ap");
  gMIXErrorCosQF->SetLineColor(3);
  gMIXErrorCosQF->SetFillColor(3);
  gMIXErrorCosQF->Draw("5same");
  grCosL1405->Draw("p");
  gMIXErrorCosL1405->Draw("5");
  //grCosL1520->Draw("p");
  //gMIXErrorCosL1520->Draw("5");
  //CS_CosS1385Lpim->Draw("same");

  gPad->SetLeftMargin(0.15);
  TLine *pL1405 = new TLine(0.6,0,1,0);
  pL1405->SetLineColor(1);
  //p->SetLineWidth(2.0);
  pL1405->SetLineStyle(2);
  pL1405->Draw();

  
  TFile *f = TFile::Open("yamagataL1405.root");
  TGraph *gry = (TGraph*)f->Get("gr_yamagata");
  gry->SetLineColor(2);
  gry->SetLineWidth(2);
  gry->Scale(1./2.7);
  //gry->Draw("c");

  TFile *fele = TFile::Open("feleangle.root");
  fele->cd();
  TGraphAsymmErrors *grsumn = (TGraphAsymmErrors*)fele->Get("grsumn");
  grsumn->SetLineColor(4);
  grsumn->SetFillColor(4);
  //grsumn->Scale(0.025*2.7*1.05);
  //grsumn->Scale(0.025*2.7*0.92);
  //const double scaleP = 0.068;
  const double scaleP = 0.01;
  grsumn->Scale(scaleP);
  grsumn->Draw("c3");
  //TGraphAsymmErrors *grK0n = (TGraphAsymmErrors*)fele->Get("grK0n");
  //TGraphAsymmErrors *grKmn = (TGraphAsymmErrors*)fele->Get("grKmn");
  //grK0n->Scale(0.05);
  //grKmn->Scale(0.015);
  //grK0n->Draw("c");
  //grKmn->Draw("c");
   
  //TLatex *tex = new TLatex();
  //double texmax = grCosL1405->GetMaximum();
  //tex->SetTextSize(0.05);
  //tex->SetTextColor(1);
  //tex->DrawLatex( 0.63,1000 , "(a)" );
  TLegend *leg = new TLegend(0.18,0.50,0.70,0.8);
  //leg->AddEntry(gMIXErrorCosQF,"#splitline{K^{-}d #rightarrow #pi^{+}#Sigma^{-}n & K^{-}d #rightarrow #pi^{-}#Sigma^{+}n}{average} ","l");
  //leg->AddEntry(gMIXErrorCosL1405,"#splitline{K^{-}d #rightarrow #pi^{+}#Sigma^{-}n & K^{-}d #rightarrow #pi^{-}#Sigma^{+}n}{average} ","l");
  leg->AddEntry(gMIXErrorCosQF,"#splitline{K^{-}d #rightarrow #pi^{+}#Sigma^{-}n, #pi^{-}#Sigma^{+}n  average}{M:1440-1500}","l");
  leg->AddEntry(gMIXErrorCosL1405,"#splitline{K^{-}d #rightarrow #pi^{+}#Sigma^{-}n, #pi^{-}#Sigma^{+}n  average}{M:1365-1425}","l");
  leg->AddEntry(grsumn,Form ("#splitline{(K^{-}p #rightarrow #bar{K}^{0}n) & (K^{-}n #rightarrow K^{-}n)}{average (x %0.3f scaled)}",scaleP/2.),"l");
  leg->SetLineWidth(0);
  leg->Draw();
  
  c1->SaveAs("angle.pdf");
}
