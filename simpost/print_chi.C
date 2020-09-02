{
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(111111);

  TFile *f = new TFile("kpp_50_50_Sp_reconst.root");
  //TFile *f = new TFile("Lpn_reconst.root");
  //TFile *f = new TFile("yamaga_reconst.root");
  //TFile *f = new TFile("tmp.root");

  TCanvas *c1 = new TCanvas("c1", "", 1200, 600);
  c1->Divide(2,1);
  c1->cd(1);
  his = (TH1F*)f->Get(Form("Lpn_kfit_chi2"));
  his->Draw();
  his->SetXTitle("#chi^{2}");
  TF1 *f1 = new TF1("f1","[0]*ROOT::Math::chisquared_pdf(x, [1])"); 
  f1->SetParameter(0, 100);
  f1->SetParameter(1, 2);
  his->Fit("f1");

  c1->cd(2);
  his = (TH1F*)f->Get(Form("Lpn_kfit_pvalue"));
  his->Draw();
  his->SetXTitle("p-value");
  cerr<<his->Integral(6,100)<<" "<<his->Integral()<<" "<<(double)his->Integral(6,100)/his->Integral()<<endl;;
  
  c1->Print("tmp.pdf");

  //double v1 = ROOT::Math::chisquared_pdf(3, 2);
  //double v2 = ROOT::Math::chisquared_cdf_c(3, 2);
  //cerr<<v1<<" "<<v2<<endl;

}
