void CalcCDH3trgEff()
{
  TFile *_file0 = TFile::Open("evanaIMpisigma_v225.root");
  TH2D* lasttime_mul_CDH = (TH2D*)_file0->Get("lasttime_mul_CDH");
  TH2D* lasttime_mul_CDH_CDH3trg = (TH2D*)_file0->Get("lasttime_mul_CDH_CDH3trg");
  TCanvas *c1 = new TCanvas("c1","c1");
  TH1D* lasttime = (TH1D*)lasttime_mul_CDH->ProjectionY("lasttime");
  lasttime->SetXTitle("CDH 3rd hit time [ns]");
  lasttime->GetXaxis()->CenterTitle();
  lasttime->Draw("HIST");
  TH1D* lasttime_fire = (TH1D*)lasttime_mul_CDH_CDH3trg->ProjectionY("lasttime_fire");
  lasttime_fire->SetLineColor(2);
  lasttime_fire->Draw("HISTsame");
  c1->SetLogy();

  int bin20ns = lasttime->FindBin(20);
  int bin50ns = lasttime->FindBin(50);
  double deno = lasttime->Integral(bin20ns,bin50ns-1);
  double nume = lasttime_fire->Integral(bin20ns,bin50ns-1);
  std::cout << "fired: " << deno << std::endl;
  std::cout << "req: " << nume << std::endl;
  double up = TEfficiency::Wilson(deno,nume,0.90,1);
  double down = TEfficiency::Wilson(deno,nume,0.90,0);
  std::cout << "DAQ trg. eff: (90% confidence interval)  " << down << " - " << up << std::endl;


}
