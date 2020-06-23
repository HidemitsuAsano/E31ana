#include "../post/weightfunc.h"
#include "../post/anacuts.h"

//woK0
TF1 *fweight_q_v363 = NULL;
TF1 *fweight_MMnmiss_v362 = NULL;
TF1 *fweight_nmom_v353 = NULL;
TF1 *fweight_IMnpip_v346 = NULL;
TF1 *fweight_IMnpim_v356 = NULL;
TF1 *fweight_IMpippim_v357 = NULL;


//wK0

TF1 *fweight_q_wK0_v377 = NULL;
TF1 *fweight_MMnmiss_wK0_v378 = NULL;
TF1 *fweight_nmom_wK0_v379 = NULL;
TF1 *fweight_IMnpip_wK0_v380 = NULL;
TF1 *fweight_IMnpim_wK0_v381 = NULL;


Double_t func_qmul(Double_t *x,Double_t *par)
{
  const Double_t xx=x[0];
  return fweight_q_v363->Eval(xx);
}

Double_t func_MMmul(Double_t *x,Double_t *par)
{
  const Double_t xx=x[0];
  return 
  fweight_MMnmiss_v362->Eval(xx);
}

Double_t func_nmommul(Double_t *x,Double_t *par)
{
  const Double_t xx=x[0];
  return 
  fweight_nmom_v353->Eval(xx);
}

Double_t func_IMnpipmul(Double_t *x,Double_t *par)
{
  const Double_t xx=x[0];
  return 
  fweight_IMnpip_v346->Eval(xx);
}


Double_t func_IMnpimmul(Double_t *x,Double_t *par)
{
  const Double_t xx=x[0];
  return 
  fweight_IMnpim_v356->Eval(xx);
}


Double_t func_IMpippimmul(Double_t *x,Double_t *par)
{
  const Double_t xx=x[0];
  return
  fweight_IMpippim_v357->Eval(xx);
}




void FakeMCWeight()
{
  //TFile *forg = TFile::Open("comp_fakedata_out_org.root","READ");
  //TFile *f = TFile::Open("comp_fakedata_out_v346.root","READ");

  TCanvas *c_woK0_func = new TCanvas("c_woK0_func","c_woK0_func",1800,1000);
  c_woK0_func->Divide(3,2);
  
  //q
  c_woK0_func->cd(1);
  fweight_q_v363 = new TF1("fweight_q_v363",func_q_mod,0,1.5,8);
  fweight_q_v363->SetParameters(param_q_mod);
  //q
  TF1* f_qmul = new TF1("f_qmul",func_qmul,0,1.5,8);
  f_qmul->SetTitle("");
  f_qmul->GetXaxis()->SetTitle("q [GeV/c]");
  f_qmul->GetXaxis()->CenterTitle();
  f_qmul->SetMinimum(0);
  f_qmul->Draw("");
  std::cout << __LINE__ << std::endl;
  //TCanvas *cqmul = new TCanvas("cqmul","cqmul");
  //cqmul->cd();
  //std::cout << __LINE__ << std::endl;
  //TH1D* h_qmul = (TH1D*)f_qmul->GetHistogram();
  //h_qmul->Draw();
  //std::cout << __LINE__ << std::endl;
  //h_qmul->Fit(f_qmul,"","",0,1.5);
 // std::cout << __LINE__ << std::endl;

  //MMnmiss
  c_woK0_func->cd(2);
  fweight_MMnmiss_v362 = new TF1("fweight_MMnmiss_v362",func_MMnmiss_mod,0,1.5,20);
  fweight_MMnmiss_v362->SetParameters(param_MMnmiss_mod);
  
  TF1* f_MMnmissmul = new TF1("MissMass",func_MMmul,0,1.5,20);
  //f_MMnmissmul->SetNpx(1000);
  f_MMnmissmul->SetTitle("");
  f_MMnmissmul->GetXaxis()->SetTitle("Miss. Mass [GeV/c^{2}]");
  f_MMnmissmul->GetXaxis()->CenterTitle();
  //f_MMnmissmul->SetMinimum(0);
  f_MMnmissmul->Draw("");
  
  TBox *box_neutron = new TBox(anacuts::neutron_MIN,0,anacuts::neutron_MAX,3);
  box_neutron->SetFillColor(4);
  box_neutron->SetFillStyle(3002);
  box_neutron->Draw();
  
  
  //nmom
  c_woK0_func->cd(3);

  fweight_nmom_v353 = new TF1("fweight_nmom_v353",func_nmom_mod,0,1.0,12);
  fweight_nmom_v353->SetParameters(param_nmom_mod);

  TF1* f_nmommul = new TF1("f_nmommul",func_nmommul,0,1.0,12);
  //f_nmommul->SetNpx(1000);
  f_nmommul->SetTitle("");
  f_nmommul->GetXaxis()->SetTitle("n_{CDS} mom. [GeV/c^{2}]");
  f_nmommul->GetXaxis()->CenterTitle();
  f_nmommul->Draw();
  TLine *nmomline = new TLine(anacuts::nmomcut,0,anacuts::nmomcut,9);
  nmomline->SetLineColor(3);
  nmomline->SetLineWidth(2.0);
  nmomline->SetLineStyle(10);
  nmomline->Draw();

  //IMnpip
  c_woK0_func->cd(4);
   
  fweight_IMnpip_v346 = new TF1("fweight_IMnpip_v346",func_IMnpipmul_s,0,2.0,12);
  fweight_IMnpip_v346->SetParameters(param_IMnpip_s);
  TF1* f_IMnpipmul = new TF1("f_IMnpipmul",func_IMnpipmul,1.06,2.0,12);
  f_IMnpipmul->SetTitle("");
  f_IMnpipmul->GetXaxis()->SetTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  f_IMnpipmul->GetXaxis()->CenterTitle();
  f_IMnpipmul->SetLineColor(2);
  f_IMnpipmul->Draw("");
  TBox *box_sigmap = new TBox(anacuts::Sigmap_MIN,0,anacuts::Sigmap_MAX,3);
  box_sigmap->SetFillColor(4);
  box_sigmap->SetFillStyle(3002);
  box_sigmap->Draw();

  c_woK0_func->cd(5);
  
  fweight_IMnpim_v356 = new TF1("fweight_IMnpim_v356",func_IMnpim_mod,1,2.0,15);
  fweight_IMnpim_v356->SetParameters(param_IMnpim_mod);
  TF1* f_IMnpimmul = new TF1("f_IMnpimmul",func_IMnpimmul,1.08,2.0,15);
  f_IMnpimmul->SetTitle("");
  f_IMnpimmul->GetXaxis()->SetTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  f_IMnpimmul->GetXaxis()->CenterTitle();
  f_IMnpimmul->Draw("");
  TBox *box_sigmam = new TBox(anacuts::Sigmam_MIN,0,anacuts::Sigmam_MAX,3);
  box_sigmam->SetFillColor(4);
  box_sigmam->SetFillStyle(3002);
  box_sigmam->Draw();

  //IMpippim
  c_woK0_func->cd(6);
  
  fweight_IMpippim_v357 = new TF1("fweight_IMpippim_v357",func_IMpippim_mod,0,1.0,15);
  fweight_IMpippim_v357->SetParameters(param_IMpippim_mod);

  TF1 *f_IMpippimmul = new TF1("f_IMpippimmul",func_IMpippimmul,0,1.0,15);
  f_IMpippimmul->SetTitle("");
  f_IMpippimmul->GetXaxis()->SetTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  f_IMpippimmul->GetXaxis()->CenterTitle();
  f_IMpippimmul->Draw("c");
  //TBox *box_K0 = new TBox(anacuts::pipi_MIN,0,anacuts::pipi_MAX,3);
  //box_K0->SetFillColor(4);
  //box_K0->SetFillStyle(3002);
  //box_K0->Draw();

  
  TCanvas *c_wK0_func = new TCanvas("c_wK0_func","c_wK0_func",1800,1000);
  c_wK0_func->Divide(3,2);
  
  //q
  c_wK0_func->cd(1);
  TF1* fweight_q_wK0_v377 = new TF1("fweight_q_wK0_v377",func_q_wK0_mod,0,1.5,8);
  fweight_q_wK0_v377->SetParameters(param_q_wK0_mod);
  fweight_q_wK0_v377->SetTitle("");
  fweight_q_wK0_v377->GetXaxis()->SetTitle("q [GeV/c]");
  fweight_q_wK0_v377->GetXaxis()->CenterTitle();
  fweight_q_wK0_v377->Draw("");


  //MMnmiss
  c_wK0_func->cd(2);
  TF1* fweight_MMnmiss_wK0_v378 = new TF1("fweight_MMnmiss_wK0_v378",func_MMnmiss_wK0_mod,0,1.5,11);
  fweight_MMnmiss_wK0_v378->SetParameters(param_MMnmiss_wK0_mod);
  fweight_MMnmiss_wK0_v378->SetTitle("");
  fweight_MMnmiss_wK0_v378->GetXaxis()->SetTitle("Miss. Mass [GeV/c^{2}]");
  fweight_MMnmiss_wK0_v378->GetXaxis()->CenterTitle();
  fweight_MMnmiss_wK0_v378->Draw("");
  box_neutron->Draw();
  //nmom
  c_wK0_func->cd(3);
  TF1* fweight_nmom_wK0_v379 = new TF1("fweight_nmom_wK0_v379",func_nmom_mod,0,1.0,12);
  fweight_nmom_wK0_v379->SetParameters(param_nmom_wK0_mod);
  //fweight_nmom_wK0_v379->SetNpx(1000);
  fweight_nmom_wK0_v379->SetTitle("");
  fweight_nmom_wK0_v379->GetXaxis()->SetTitle("n_{CDS} mom. [GeV/c^{2}]");
  fweight_nmom_wK0_v379->GetXaxis()->CenterTitle();
  fweight_nmom_wK0_v379->GetYaxis()->SetRangeUser(0,53);
  fweight_nmom_wK0_v379->Draw();
  
  TLine *nmomline2 = new TLine(anacuts::nmomcut,0,anacuts::nmomcut,53);
  nmomline2->SetLineColor(3);
  nmomline2->SetLineWidth(2.0);
  nmomline2->SetLineStyle(10);
  nmomline2->Draw();
  //IMnpip
  c_wK0_func->cd(4);
  TF1* fweight_IMnpip_wK0_v380 = new TF1("fweight_IMnpip_wK0_v380",func_IMnpip_wK0_mod,1,2,16);
  fweight_IMnpip_wK0_v380->SetParameters(param_IMnpip_wK0_mod);
  fweight_IMnpip_wK0_v380->SetTitle("");
  fweight_IMnpip_wK0_v380->GetXaxis()->SetTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  fweight_IMnpip_wK0_v380->GetXaxis()->CenterTitle();
  fweight_IMnpip_wK0_v380->Draw("");
  box_sigmap->Draw();
  
  std::cout << __LINE__ << std::endl;
  c_wK0_func->cd(5);
  TF1* fweight_IMnpim_wK0_v381 = new TF1("fweight_IMnpim_wK0_v381",func_IMnpim_wK0_mod,1,2,11);
  fweight_IMnpim_wK0_v381->SetParameters(param_IMnpim_wK0_mod);
  fweight_IMnpim_wK0_v381->SetTitle("");
  fweight_IMnpim_wK0_v381->GetXaxis()->SetTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  fweight_IMnpim_wK0_v381->GetXaxis()->CenterTitle();
  fweight_IMnpim_wK0_v381->Draw("");
  box_sigmam->Draw();
  //IMpippim (N/A)
  //c_wK0_func->cd(6);
  
  


};
