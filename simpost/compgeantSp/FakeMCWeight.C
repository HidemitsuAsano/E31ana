#include "../../post/weightfuncGSp.h"
#include "../../post/anacuts.h"

//woK0
TF1 *fweight_q = NULL;
TF1 *fweight_MMnmiss = NULL;
TF1 *fweight_nmom = NULL;
TF1 *fweight_IMnpip= NULL;
TF1 *fweight_IMnpim = NULL;
TF1 *fweight_IMpippim = NULL;
TF1 *fweight_Mompippim = NULL;
TF1 *fweight_Momnpip = NULL;
TF1 *fweight_Momnpim = NULL;

//wK0

TF1 *fweight_q_wK0 = NULL;
TF1 *fweight_MMnmiss_wK0 = NULL;
TF1 *fweight_nmom_wK0 = NULL;
TF1 *fweight_IMnpip_wK0 = NULL;
TF1 *fweight_IMnpim_wK0 = NULL;



void FakeMCWeight()
{
  //TFile *forg = TFile::Open("comp_fakedata_out_org.root","READ");
  //TFile *f = TFile::Open("comp_fakedata_out_v346.root","READ");

  TCanvas *c_woK0_func = new TCanvas("c_woK0_func","c_woK0_func",1800,1000);
  c_woK0_func->Divide(4,2);
  
  //q
  /*
  c_woK0_func->cd(1);
  fweight_q = new TF1("fweight_q",func_q_mod,0,1.5,8);
  fweight_q->SetParameters(param_q_mod);
  fweight_q->SetTitle("");
  fweight_q->SetMinimum(0);
  fweight_q->SetNpx(10000);
  fweight_q->GetXaxis()->SetTitle("q [GeV/c]");
  fweight_q->GetXaxis()->CenterTitle();
  fweight_q->Draw("c");
  std::cout << __LINE__ << std::endl;
  */

  //MMnmiss
  c_woK0_func->cd(1);
  fweight_MMnmiss = new TF1("fweight_MMnmiss",func_MMnmiss_mod,0,1.5,20);
  fweight_MMnmiss->SetParameters(param_MMnmiss_mod);
  fweight_MMnmiss->SetNpx(10000);
  fweight_MMnmiss->SetTitle("");
  fweight_MMnmiss->GetXaxis()->SetTitle("Miss. Mass [GeV/c^{2}]");
  fweight_MMnmiss->GetXaxis()->CenterTitle();
  fweight_MMnmiss->Draw("c");
  TBox *box_neutron = new TBox(anacuts::neutron_MIN,0,anacuts::neutron_MAX,3);
  box_neutron->SetFillColor(4);
  box_neutron->SetFillStyle(3002);
  box_neutron->Draw();
  
  
  //nmom
  c_woK0_func->cd(2);
  fweight_nmom = new TF1("fweight_nmom",func_nmom_mod,0,1.0,12);
  fweight_nmom->SetParameters(param_nmom_mod);
  fweight_nmom->SetNpx(10000);
  fweight_nmom->SetTitle("");
  fweight_nmom->GetXaxis()->SetTitle("n_{CDS} mom. [GeV/c^{2}]");
  fweight_nmom->GetXaxis()->CenterTitle();
  fweight_nmom->Draw("c");
  TLine *nmomline = new TLine(anacuts::nmomcut,0,anacuts::nmomcut,10);
  nmomline->SetLineColor(3);
  nmomline->SetLineWidth(2.0);
  nmomline->SetLineStyle(10);
  nmomline->Draw();

  //IMnpip
  c_woK0_func->cd(3);
  fweight_IMnpip = new TF1("fweight_IMnpip",func_IMnpip_mod,1,2.0,6);
  fweight_IMnpip->SetParameters(param_IMnpip_mod);
  fweight_IMnpip->SetTitle("");
  fweight_IMnpip->SetLineColor(2);
  fweight_IMnpip->SetNpx(10000);
  fweight_IMnpip->GetXaxis()->SetTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  fweight_IMnpip->GetXaxis()->CenterTitle();
  fweight_IMnpip->Draw("c");
  TBox *box_sigmap = new TBox(anacuts::Sigmap_MIN,0,anacuts::Sigmap_MAX,5);
  box_sigmap->SetFillColor(4);
  box_sigmap->SetFillStyle(3002);
  box_sigmap->Draw();
  
  //IMnpim
  c_woK0_func->cd(5);
  fweight_IMnpim = new TF1("fweight_IMnpim",func_IMnpim_mod,1,2.0,10);
  fweight_IMnpim->SetParameters(param_IMnpim_mod);
  fweight_IMnpim->SetTitle("");
  fweight_IMnpim->SetNpx(10000);
  fweight_IMnpim->GetXaxis()->SetTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  fweight_IMnpim->GetXaxis()->CenterTitle();
  fweight_IMnpim->Draw("c");
  TBox *box_sigmam = new TBox(anacuts::Sigmam_MIN,0,anacuts::Sigmam_MAX,3);
  box_sigmam->SetFillColor(4);
  box_sigmam->SetFillStyle(3002);
  box_sigmam->Draw();

  //IMpippim
  c_woK0_func->cd(7);
  fweight_IMpippim = new TF1("fweight_IMpippim",func_IMpippim_mod,0,1.0,15);
  fweight_IMpippim->SetParameters(param_IMpippim_mod);
  fweight_IMpippim->SetTitle("");
  fweight_IMpippim->SetNpx(10000);
  fweight_IMpippim->GetXaxis()->SetTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  fweight_IMpippim->GetXaxis()->CenterTitle();
  fweight_IMpippim->Draw("c");
  //TBox *box_K0 = new TBox(anacuts::pipi_MIN,0,anacuts::pipi_MAX,3);
  //box_K0->SetFillColor(4);
  //box_K0->SetFillStyle(3002);
  //box_K0->Draw();

  c_woK0_func->cd(8);
  fweight_Mompippim = new TF1("fweight_Mompippim",func_Mompippim,0,1.5,7);
  fweight_Mompippim->SetParameters(param_Mompippim);
  fweight_Mompippim->SetTitle("");
  fweight_Mompippim->SetNpx(10000);
  fweight_Mompippim->GetXaxis()->SetTitle("Mom(#pi^{+}#pi^{-}) [GeV/c]");
  fweight_Mompippim->GetXaxis()->CenterTitle();
  fweight_Mompippim->Draw("c");


  c_woK0_func->cd(4);
  fweight_Momnpip = new TF1("fweight_Momnpip",func_Momnpip,0,1.5,6);
  fweight_Momnpip->SetParameters(param_Momnpip);
  fweight_Momnpip->SetTitle("");
  fweight_Momnpip->SetNpx(10000);
  fweight_Momnpip->GetXaxis()->SetTitle("Mom(n#pi^{+}) [GeV/c]");
  fweight_Momnpip->GetXaxis()->CenterTitle();
  fweight_Momnpip->Draw("c");


  c_woK0_func->cd(6);
  fweight_Momnpim = new TF1("fweight_Momnpim",func_Momnpim,0,1.5,8);
  fweight_Momnpim->SetParameters(param_Momnpim);
  fweight_Momnpim->SetTitle("");
  fweight_Momnpim->SetNpx(10000);
  fweight_Momnpim->GetXaxis()->SetTitle("Mom(n#pi^{-}) [GeV/c]");
  fweight_Momnpim->GetXaxis()->CenterTitle();
  fweight_Momnpim->Draw("c");



  TCanvas *c_wK0_func = new TCanvas("c_wK0_func","c_wK0_func",1800,1000);
  c_wK0_func->Divide(3,2);
  
  //q
  c_wK0_func->cd(1);
  TF1* fweight_q_wK0 = new TF1("fweight_q_wK0",func_q_wK0_mod,0,1.5,8);
  fweight_q_wK0->SetParameters(param_q_wK0_mod);
  fweight_q_wK0->SetTitle("");
  fweight_q_wK0->SetNpx(10000);
  fweight_q_wK0->GetXaxis()->SetTitle("q [GeV/c]");
  fweight_q_wK0->GetXaxis()->CenterTitle();
  fweight_q_wK0->Draw("c");


  //MMnmiss
  c_wK0_func->cd(2);
  TF1* fweight_MMnmiss_wK0 = new TF1("fweight_MMnmiss_wK0",func_MMnmiss_wK0_mod,0,1.5,20);
  fweight_MMnmiss_wK0->SetParameters(param_MMnmiss_wK0_mod);
  fweight_MMnmiss_wK0->SetTitle("");
  fweight_MMnmiss_wK0->SetNpx(10000);
  fweight_MMnmiss_wK0->GetXaxis()->SetTitle("Miss. Mass [GeV/c^{2}]");
  fweight_MMnmiss_wK0->GetXaxis()->CenterTitle();
  fweight_MMnmiss_wK0->Draw("c");
  box_neutron->Draw();
  //nmom
  c_wK0_func->cd(3);
  TF1* fweight_nmom_wK0 = new TF1("fweight_nmom_wK0",func_nmom_mod,0,1.0,12);
  fweight_nmom_wK0->SetParameters(param_nmom_wK0_mod);
  fweight_nmom_wK0->SetNpx(10000);
  fweight_nmom_wK0->SetTitle("");
  fweight_nmom_wK0->GetXaxis()->SetTitle("n_{CDS} mom. [GeV/c^{2}]");
  fweight_nmom_wK0->GetXaxis()->CenterTitle();
  //fweight_nmom_wK0->GetYaxis()->SetRangeUser(0,53);
  fweight_nmom_wK0->Draw("c");
  
  TLine *nmomline2 = new TLine(anacuts::nmomcut,0,anacuts::nmomcut,53);
  nmomline2->SetLineColor(3);
  nmomline2->SetLineWidth(2.0);
  nmomline2->SetLineStyle(10);
  nmomline2->Draw();

  //IMnpip
  c_wK0_func->cd(4);
  TF1* fweight_IMnpip_wK0 = new TF1("fweight_IMnpip_wK0",func_IMnpip_wK0_mod,1,2,14);
  fweight_IMnpip_wK0->SetParameters(param_IMnpip_wK0_mod);
  fweight_IMnpip_wK0->SetTitle("");
  fweight_IMnpip_wK0->SetNpx(10000);
  fweight_IMnpip_wK0->GetXaxis()->SetTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  fweight_IMnpip_wK0->GetXaxis()->CenterTitle();
  fweight_IMnpip_wK0->Draw("c");
  box_sigmap->Draw();
  
  std::cout << __LINE__ << std::endl;
  c_wK0_func->cd(5);
  TF1* fweight_IMnpim_wK0 = new TF1("fweight_IMnpim_wK0_v381",func_IMnpim_wK0_mod,1,2,10);
  fweight_IMnpim_wK0->SetParameters(param_IMnpim_wK0_mod);
  fweight_IMnpim_wK0->SetTitle("");
  fweight_IMnpim_wK0->SetNpx(10000);
  fweight_IMnpim_wK0->GetXaxis()->SetTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  fweight_IMnpim_wK0->GetXaxis()->CenterTitle();
  fweight_IMnpim_wK0->Draw("c");
  box_sigmam->Draw();
  //IMpippim (N/A)
  //c_wK0_func->cd(6);
  
  


};
