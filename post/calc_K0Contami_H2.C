void calc_K0Contami_H2()
{
  const double K0CS = 6454.0 ; //ub
  
  TFile *f = TFile::Open("../simpost/simIMsigma_H2_K0n_npi_v2_out_dE2_iso_rej_nostop.root");
  TFile *fgen = TFile::Open("../simpost/simIMsigma_H2_K0n_v2.root");
  TFile *facc = TFile::Open("accH2dE2.root");


  TCanvas *c1 = new TCanvas("c1","c1");
  TH1D* ReactCosCM_0 = (TH1D*)fgen->Get("ReactCosCM_0");
  ReactCosCM_0->SetXTitle("Cos_{CM} K0");
  ReactCosCM_0->GetXaxis()->CenterTitle();
  //ReactCosCM_0->SetYTitle("d#sigma/d(cos) ");
  ReactCosCM_0->Draw();
  double Ent = ReactCosCM_0->GetEntries();
  std::cout << "gen. entries" << Ent << std::endl;
  // d(sigma)/d(cos) -> d(sigma)/d(Omega) = *(2/NBIN)/(2pi*2/NBIN) = *(1./2pi)

  TCanvas *c2 = new TCanvas("c2","c2");
  TH1D* Cospicm_pi_Sp = (TH1D*)f->Get("Cospicm_pi_Sp");
  Cospicm_pi_Sp->RebinX(5);
  Cospicm_pi_Sp->SetYTitle("counts");
  Cospicm_pi_Sp->Draw();

  TCanvas *c2_2 = new TCanvas("c2_2","c2_2");
  double CospiBinW = Cospicm_pi_Sp->GetBinWidth(2);
  const double cosToStrbin = 2.*3.1415*(CospiBinW); 
  TH1D* COspicm_pi_Sp_CS = (TH1D*)Cospicm_pi_Sp->Clone("Cospicm_pi_Sp_CS");
  Cospicm_pi_Sp_CS->Scale(K0CS/Ent/cosToStrbin);
  Cospicm_pi_Sp_CS->SetYTitle("d#sigma/d#Omega [#mub/sr]");
  Cospicm_pi_Sp_CS->GetYaxis()->CenterTitle();
  TH1D* accCosSp2 = (TH1D*)facc->Get("accCosSp2");
  Cospicm_pi_Sp_CS->Divide(accCosSp2);
  Cospicm_pi_Sp_CS->GetXaxis()->SetRangeUser(0.3,1);
 
  Cospicm_pi_Sp_CS->Draw();
  TGraphErrors *gCSK0_Sp = new TGraphErrors(Cospicm_pi_Sp_CS);
  gCSK0_Sp->SetName("gCSK0_Sp");

  TCanvas *c3 = new TCanvas("c3","c3");
  TH1D* Cospicm_pi_Sm = (TH1D*)f->Get("Cospicm_pi_Sm");
  Cospicm_pi_Sm->RebinX(5);
  Cospicm_pi_Sm->SetYTitle("counts");
  Cospicm_pi_Sm->Draw();

  TCanvas *c3_2 = new TCanvas("c3_2","c3_2");
  double CospiBinW = Cospicm_pi_Sm->GetBinWidth(2);
  TH1D* COspicm_pi_Sm_CS = (TH1D*)Cospicm_pi_Sm->Clone("Cospicm_pi_Sm_CS");
  Cospicm_pi_Sm_CS->Scale(K0CS/Ent/cosToStrbin);
  Cospicm_pi_Sm_CS->SetYTitle("d#sigma/d#Omega [#mub/sr]");
  Cospicm_pi_Sm_CS->GetYaxis()->CenterTitle();
  TH1D* accCosSm2 = (TH1D*)facc->Get("accCosSm2");
  Cospicm_pi_Sm_CS->Divide(accCosSm2);
  Cospicm_pi_Sm_CS->GetXaxis()->SetRangeUser(0.3,1);
  Cospicm_pi_Sm_CS->Draw();
  TGraphErrors *gCSK0_Sm = new TGraphErrors(Cospicm_pi_Sm_CS);
  gCSK0_Sm->SetName("gCSK0_Sm");


  TFile *fout = new TFile("CS_K0contami_H2.root","RECREATE");
  gCSK0_Sp->Write();
  gCSK0_Sm->Write();



}
