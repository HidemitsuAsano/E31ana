void CalcIntegratedLumi()
{
  //4.819 x 10e23 +/- 4.58 x 10e21 (run100-133)
  //4.830 x 10e23 +/- 4.58 x 10e21 (run134-812)
  //
  const double D2density = 4.830e+23*1.25;// #deuteron /10 cm target length * 12.5
  const double D2density_err = 4.58e+21*1.25;
  const double CDH3trgeff = 0.996122;
  const double CDH3trgeff_err = 3.0*0.001235; 
  
  //
  const double scalerKaon = 5.80154e+10;//total num.
  const double scalerKaon_err = sqrt(5.80154e+10);//total num.
  const double Survival = 0.302809;//run average
  const double Survival_err = 0.0002;//run average
  const double DAQeff = 0.772601;//run average
  const double DAQerr_err = 9.1763e-06;
  const double microbarnInv = 1e-30;
  const double InteLumi = scalerKaon*Survival*D2density*DAQeff*microbarnInv;

  
  std::cout << "IF run dep. is ignored "  << std::endl;
  std::cout << "Integrated Lumi " << InteLumi << std::endl;
  const double simBeamOK = 2.47625376000000000e+08;
  const double simBeamNO = 3.81361560000000000e+07;
  std::cout << "Sim. Beam survival factor " << simBeamOK / (simBeamOK+simBeamNO) << std::endl;
  std::cout << "Err " << sqrt(pow(D2density_err*CDH3trgeff*scalerKaon*Survival*DAQeff,2.0)
                             +pow(D2density*CDH3trgeff_err*scalerKaon*Survival*DAQeff,2.0)   
                             +pow(D2density*CDH3trgeff*scalerKaon_err*Survival*DAQeff,2.0)   
                             +pow(D2density*CDH3trgeff*scalerKaon*Survival_err*DAQeff,2.0) 
      )*microbarnInv << std::endl;
  std::cout << std::endl;
  
  TFile *fdaq = TFile::Open("daqeff.root");
  TFile *fscalerK = TFile::Open("lumi.root");
  TFile *fbeam = TFile::Open("beamSurvival.root");
  
  TH1D* heff  = (TH1D*)fdaq->Get("heff");
  TH1D* hnk   = (TH1D*)fscalerK->Get("hnK");
  TH1D* hbeam = (TH1D*)fbeam->Get("hBeamSurvivalR");
  TH1D* hmulti =(TH1D*)heff->Clone("hmulti");
  hmulti->Multiply(hnk);
  hmulti->Multiply(hbeam);

  TH1D* hInteLumi = (TH1D*)heff->Clone("hInteLumi");
  for(int ibin=0;ibin<heff->GetNbinsX();ibin++){
    double cont = hmulti->GetBinContent(ibin);
    if(cont<0.00001)continue;
    double err = hmulti->GetBinError(ibin);
    double Lumi = D2density*CDH3trgeff*cont*microbarnInv;
    double errtotal = sqrt((pow(err*D2density*CDH3trgeff,2.0)+pow(cont*D2density_err*CDH3trgeff,2.0)+pow(cont*D2density*CDH3trgeff_err,2.0)))*microbarnInv;
    hInteLumi->SetBinContent(ibin,Lumi);
    hInteLumi->SetBinError(ibin,errtotal);
  }
  hInteLumi->Draw("E");
  double errorInte=0.0;
  gStyle->SetOptStat("i");
  Double_t A = hInteLumi->IntegralAndError(1,hInteLumi->GetNbinsX(),errorInte, "");
  hInteLumi->SetTitle("Run by Run Integrated Luminosity");
  std::cout << "Integral " << A << std::endl;
  std::cout << "Integral err " << errorInte << std::endl;  
  
  TParameter<double> IntegLumi("IntegLumi",A);
  TParameter<double> IntegErr("Err",errorInte);

  TFile *fout = new TFile("InteLumi.root","RECREATE");
  IntegLumi.Write();
  IntegErr.Write();
  hInteLumi->Write();

}
