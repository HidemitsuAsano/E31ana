void compL1405Angle()
{

  TFile *ff = TFile::Open("csfinal.root");

  auto *grCosL1405 = (TGraphAsymmErrors*)ff->Get("grCosL1405");
  auto *gMIXErrorCosL1405 = (TGraphAsymmErrors*)ff->Get("gMIXErrorCosL1405");

  TCanvas *c1 = new TCanvas("c1","c1");
  grCosL1405->Draw("ap");
  gMIXErrorCosL1405->Draw("5");

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
  grsumn->Scale(0.069);
  grsumn->Draw("c3");
  TGraphAsymmErrors *grK0n = (TGraphAsymmErrors*)fele->Get("grK0n");
  TGraphAsymmErrors *grKmn = (TGraphAsymmErrors*)fele->Get("grKmn");
  //grK0n->Scale(0.05);
  //grKmn->Scale(0.015);
  //grK0n->Draw("c");
  //grKmn->Draw("c");


}
