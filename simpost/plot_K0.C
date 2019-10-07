#include "../src/GlobalVariables.h"

void plot_K0()
{

  //TFile *_file0 = TFile::Open("../post/evanaIMpisigma_npippim_v177_K015_out.root");
  TFile *_file0 = TFile::Open("../post/evanaIMpisigma_npippim_v162_outncutK015.root");
  _file0->cd();
  //TH2F* IMpippim_IMnpipi_n_rdata = (TH2F*)_file0->Get("IMpippim_IMnpipi_n");
  //TH2F* mnmom_IMpippim_n = (TH2F*)_file0->Get("mnmom_IMpippim_n");
  TH2F* q_IMpippim_n_rdata = (TH2F*)_file0->Get("q_IMpippim_n");
  TH2F* nmom_IMnpipi_wK0_n_rdata = (TH2F*)_file0->Get("nmom_IMnpipi_wK0_n");
  TH2F* q_IMnpipi_wK0_n_rdata = (TH2F*)_file0->Get("q_IMnpipi_wK0_n");
  TH2F* Mompippim_IMnpipi_dE_wK0_n_rdata = (TH2F*)_file0->Get("Mompippim_IMnpipi_dE_wK0_n");
  TH2F* MMnmiss_IMpippim_dE_rdata = (TH2F*)_file0->Get("MMnmiss_IMpippim_dE");

  //TCanvas *cnmom_IMnpipi_wK0_n_rdata = new TCanvas("cnmom_IMnpipi_wK0_n_rdata","nmom_IMnpipi_wK0_n_rdata",1000,800);
  //cnmom_IMnpipi_wK0_n_rdata->Divide(2,2);
  //cnmom_IMnpipi_wK0_n_rdata->cd(1);
  //nmom_IMnpipi_wK0_n_rdata->Draw("colz");
  const double Kp_mass = pMass + kpMass;  
  TF1 *fkp = new TF1("f", "sqrt(((x*x-[0]*[0]-[1]*[1])/(2*[0]))*((x*x-[0]*[0]-[1]*[1])/(2*[0]))-[1]*[1])",Kp_mass-0.0001,2);
  fkp->SetParameter(0,nMass);
  fkp->SetParameter(1,kpMass);
  //fkp->SetLineColor(4);
  fkp->SetLineWidth(5);
  fkp->SetLineStyle(4);
  fkp->SetLineColorAlpha(kPink, 0.35);
  // fkp->Draw("same");
  //cnmom_IMnpipi_wK0_n_rdata->cd(3);
  TH1D* IMnpipi_rdata = (TH1D*)nmom_IMnpipi_wK0_n_rdata->ProjectionX("IMnpipi_rdata");
  //IMnpipi_rdata->Draw();
  //cnmom_IMnpipi_wK0_n_rdata->cd(2);
  TH1D* nmom_IMnpipi_wK0_n_rdata_py = (TH1D*) nmom_IMnpipi_wK0_n_rdata->ProjectionY("nmom_IMnpipi_wK0_n_rdata_py");
  //nmom_IMnpipi_wK0_n_rdata_py->Draw();
  
  //TCanvas *cq_IMnpipi_wK0_n_rdata = new TCanvas("cq_IMnpipi_wK0_n_rdata","q_IMnpipi_wK0_n_rdata",1000,800);
  //cq_IMnpipi_wK0_n_rdata->Divide(2,2);
  //cq_IMnpipi_wK0_n_rdata->cd(1);
  //q_IMnpipi_wK0_n_rdata->Draw("colz");
  //fkp->Draw("same");
  //cq_IMnpipi_wK0_n_rdata->cd(3);
  //q_IMnpipi_wK0_n_rdata->ProjectionX()->Draw();
  //cq_IMnpipi_wK0_n_rdata->cd(2);
  TH1D* q_rdata = (TH1D*)q_IMnpipi_wK0_n_rdata->ProjectionY("q_rdata");
  //q_rdata->Draw();
  
  //TCanvas* cMompippim_IMnpipi_dE_wK0_n_rdata = new TCanvas("cMompippim_IMnpipi_dE_wK0_n_rdata","cMompippim_IMnpipi_dE_wK0_n_rdata",1000,800);
  //cMompippim_IMnpipi_dE_wK0_n_rdata->Divide(2,2);
  //cMompippim_IMnpipi_dE_wK0_n_rdata->cd(1);
  //Mompippim_IMnpipi_dE_wK0_n_rdata->Draw("colz");
  //cMompippim_IMnpipi_dE_wK0_n_rdata->cd(3);
  //Mompippim_IMnpipi_dE_wK0_n_rdata->ProjectionX()->Draw();
  //cMompippim_IMnpipi_dE_wK0_n_rdata->cd(2);
  TH1D* Mompippim_rdata = (TH1D*)Mompippim_IMnpipi_dE_wK0_n_rdata->ProjectionY("Mompippim_rdata"); 
  //Mompippim_rdata->Draw();
   
  const double pipi_MIN = 0.487;// (almost 1.5 sigma)
  const double pipi_MAX = 0.510;// (almost 1.5 sigma)
  int pipimin = MMnmiss_IMpippim_dE_rdata->GetXaxis()->FindBin(pipi_MIN);
  int pipimax = MMnmiss_IMpippim_dE_rdata->GetXaxis()->FindBin(pipi_MAX);
  TH1D* MMnmiss_rdata = (TH1D*) MMnmiss_IMpippim_dE_rdata->ProjectionY("MMnmiss_rdata",pipimin,pipimax) ;

  //K-p -> K0n -> n-n elastic
  TFile *_file1 = TFile::Open("simIMpisigma_K0_nnts_pippimn_v10_outncutK015.root");
  //TFile *_file1 = TFile::Open("simIMpisigma_K0_nnts_pippimn_v7_outncutK015.root");
  TH2F* q_IMpippim_n_nnts = (TH2F*)_file1->Get("q_IMpippim_n");
  TH2F* nmom_IMnpipi_wK0_n_nnts = (TH2F*)_file1->Get("nmom_IMnpipi_wK0_n");
  TH2F* q_IMnpipi_wK0_n_nnts = (TH2F*)_file1->Get("q_IMnpipi_wK0_n");
  TH2F* Mompippim_IMnpipi_dE_wK0_n_nnts = (TH2F*)_file1->Get("Mompippim_IMnpipi_dE_wK0_n");
  TH2F* MMnmiss_IMpippim_dE_nnts = (TH2F*)_file1->Get("MMnmiss_IMpippim_dE");

  //TCanvas *cnmom_IMnpipi_wK0_n_nnts = new TCanvas("cnmom_IMnpipi_wK0_n_nnts","nmom_IMnpipi_wK0_n_nnts",1000,800);
  //cnmom_IMnpipi_wK0_n_nnts->Divide(2,2);
  //cnmom_IMnpipi_wK0_n_nnts->cd(1);
  //nmom_IMnpipi_wK0_n_nnts->Draw("colz");
  //cnmom_IMnpipi_wK0_n_nnts->cd(3);
  TH1D* IMnpipi_nnts = (TH1D*) nmom_IMnpipi_wK0_n_nnts->ProjectionX("IMnpipi_nnts");
  //IMnpipi_nnts->Draw();
  //cnmom_IMnpipi_wK0_n_nnts->cd(2);
  TH1D* nmom_IMnpipi_wK0_n_nnts_py = (TH1D*)nmom_IMnpipi_wK0_n_nnts->ProjectionY("nmom_IMnpipi_wK0_n_nnts_py");
  //nmom_IMnpipi_wK0_n_nnts_py->Draw();
  
  //TCanvas *cq_IMnpipi_wK0_n_nnts = new TCanvas("cq_IMnpipi_wK0_n_nnts","q_IMnpipi_wK0_n_nnts",1000,800);
  //cq_IMnpipi_wK0_n_nnts->Divide(2,2);
  //cq_IMnpipi_wK0_n_nnts->cd(1);
  //q_IMnpipi_wK0_n_nnts->Draw("colz");
  //fkp->Draw("same");
  //cq_IMnpipi_wK0_n_nnts->cd(3);
  //q_IMnpipi_wK0_n_nnts->ProjectionX()->Draw();
  //cq_IMnpipi_wK0_n_nnts->cd(2);
  TH1D* q_nnts = (TH1D*)q_IMnpipi_wK0_n_nnts->ProjectionY("q_nnts");
  //q_nnts->Draw();
  
  //TCanvas* cMompippim_IMnpipi_dE_wK0_n_nnts = new TCanvas("cMompippim_IMnpipi_dE_wK0_n_nnts","cMompippim_IMnpipi_dE_wK0_n_nnts",1000,800);
  //cMompippim_IMnpipi_dE_wK0_n_nnts->Divide(2,2);
  //cMompippim_IMnpipi_dE_wK0_n_nnts->cd(1);
  //Mompippim_IMnpipi_dE_wK0_n_nnts->Draw("colz");
  //cMompippim_IMnpipi_dE_wK0_n_nnts->cd(3);
  //Mompippim_IMnpipi_dE_wK0_n_nnts->ProjectionX()->Draw();
  //cMompippim_IMnpipi_dE_wK0_n_nnts->cd(2);
  TH1D* Mompippim_nnts = (TH1D*) Mompippim_IMnpipi_dE_wK0_n_nnts->ProjectionY("Mompippim_nnts");
  //Mompippim_nnts->Draw();
  TH1D* MMnmiss_nnts = (TH1D*) MMnmiss_IMpippim_dE_nnts->ProjectionY("MMnmiss_nnts",pipimin,pipimax);

  //K-p -> K0n (ns)
  TFile *_file2 = TFile::Open("simIMpisigma_K0n_ns_pippimn_v12_outncutK015.root");
  //TFile *_file2 = TFile::Open("simIMpisigma_K0n_ns_pippimn_v7_outncutK015.root");
  TH2F* q_IMpippim_n_ns = (TH2F*)_file2->Get("q_IMpippim_n");
  TH2F* nmom_IMnpipi_wK0_n_ns = (TH2F*)_file2->Get("nmom_IMnpipi_wK0_n");
  TH2F* q_IMnpipi_wK0_n_ns = (TH2F*)_file2->Get("q_IMnpipi_wK0_n");
  TH2F* Mompippim_IMnpipi_dE_wK0_n_ns = (TH2F*)_file2->Get("Mompippim_IMnpipi_dE_wK0_n");
  TH2F* MMnmiss_IMpippim_dE_ns = (TH2F*)_file2->Get("MMnmiss_IMpippim_dE");


  //TCanvas *cnmom_IMnpipi_wK0_n_ns = new TCanvas("cnmom_IMnpipi_wK0_n_ns","nmom_IMnpipi_wK0_n_ns",1000,800);
  //cnmom_IMnpipi_wK0_n_ns->Divide(2,2);
  //cnmom_IMnpipi_wK0_n_ns->cd(1);
  //nmom_IMnpipi_wK0_n_ns->Draw("colz");
  //cnmom_IMnpipi_wK0_n_ns->cd(3);
  TH1D* IMnpipi_ns = (TH1D*) nmom_IMnpipi_wK0_n_ns->ProjectionX();
  //IMnpipi_ns->Draw();
  //cnmom_IMnpipi_wK0_n_ns->cd(2);
  TH1D* nmom_IMnpipi_wK0_n_ns_py = (TH1D*) nmom_IMnpipi_wK0_n_ns->ProjectionY("nmom_IMnpipi_wK0_n_ns_py"); 
  //nmom_IMnpipi_wK0_n_ns_py->Draw();
  
  //TCanvas *cq_IMnpipi_wK0_n_ns = new TCanvas("cq_IMnpipi_wK0_n_ns","q_IMnpipi_wK0_n_ns",1000,800);
  //cq_IMnpipi_wK0_n_ns->Divide(2,2);
  //cq_IMnpipi_wK0_n_ns->cd(1);
  //q_IMnpipi_wK0_n_ns->Draw("colz");
  //fkp->Draw("same");
  //cq_IMnpipi_wK0_n_ns->cd(3);
  //q_IMnpipi_wK0_n_ns->ProjectionX()->Draw();
  //cq_IMnpipi_wK0_n_ns->cd(2);
  TH1D* q_ns = (TH1D*) q_IMnpipi_wK0_n_ns->ProjectionY();
  //q_ns->Draw();
  
  //TCanvas* cMompippim_IMnpipi_dE_wK0_n_ns = new TCanvas("cMompippim_IMnpipi_dE_wK0_n_ns","cMompippim_IMnpipi_dE_wK0_n_ns",1000,800);
  //cMompippim_IMnpipi_dE_wK0_n_ns->Divide(2,2);
  //cMompippim_IMnpipi_dE_wK0_n_ns->cd(1);
  //Mompippim_IMnpipi_dE_wK0_n_ns->Draw("colz");
  //cMompippim_IMnpipi_dE_wK0_n_ns->cd(3);
  //Mompippim_IMnpipi_dE_wK0_n_ns->ProjectionX()->Draw();
  //cMompippim_IMnpipi_dE_wK0_n_ns->cd(2);
  TH1D* Mompippim_ns = (TH1D*)Mompippim_IMnpipi_dE_wK0_n_ns->ProjectionY();
  //Mompippim_ns->Draw();
  TH1D* MMnmiss_ns = (TH1D*) MMnmiss_IMpippim_dE_ns->ProjectionY("MMnmiss_ns",pipimin,pipimax);
  

  //K-d -> K0nn phase space
  TFile *_file3 = TFile::Open("simIMpisigma_K0nn_pippimn_v7_outncutK015.root");
  //TFile *_file3 = TFile::Open("simIMpisigma_K0nn_pippimn_v5_outncutK015.root");
  TH2F* q_IMpippim_n_K0nn = (TH2F*)_file3->Get("q_IMpippim_n");
  TH2F* nmom_IMnpipi_wK0_n_K0nn = (TH2F*)_file3->Get("nmom_IMnpipi_wK0_n");
  TH2F* q_IMnpipi_wK0_n_K0nn = (TH2F*)_file3->Get("q_IMnpipi_wK0_n");
  TH2F* Mompippim_IMnpipi_dE_wK0_n_K0nn = (TH2F*)_file3->Get("Mompippim_IMnpipi_dE_wK0_n");
  TH2F* MMnmiss_IMpippim_dE_K0nn = (TH2F*)_file3->Get("MMnmiss_IMpippim_dE");


  //TCanvas *cnmom_IMnpipi_wK0_n_K0nn = new TCanvas("cnmom_IMnpipi_wK0_n_K0nn","nmom_IMnpipi_wK0_n_K0nn",1000,800);
  //cnmom_IMnpipi_wK0_n_K0nn->Divide(2,2);
  //cnmom_IMnpipi_wK0_n_K0nn->cd(1);
  //nmom_IMnpipi_wK0_n_K0nn->Draw("colz");
  //cnmom_IMnpipi_wK0_n_K0nn->cd(3);
  TH1D* IMnpipi_K0nn = (TH1D*) nmom_IMnpipi_wK0_n_K0nn->ProjectionX();
  //IMnpipi_K0nn->Draw();
  //cnmom_IMnpipi_wK0_n_K0nn->cd(2);
  TH1D* nmom_IMnpipi_wK0_n_K0nn_py = (TH1D*)nmom_IMnpipi_wK0_n_K0nn->ProjectionY("nmom_IMnpipi_wK0_n_K0nn_py");
  //nmom_IMnpipi_wK0_n_K0nn_py->Draw();
  
  //TCanvas *cq_IMnpipi_wK0_n_K0nn = new TCanvas("cq_IMnpipi_wK0_n_K0nn","q_IMnpipi_wK0_n_K0nn",1000,800);
  //cq_IMnpipi_wK0_n_K0nn->Divide(2,2);
  //cq_IMnpipi_wK0_n_K0nn->cd(1);
  //q_IMnpipi_wK0_n_K0nn->Draw("colz");
  //fkp->Draw("same");
  //cq_IMnpipi_wK0_n_K0nn->cd(3);
  //q_IMnpipi_wK0_n_K0nn->ProjectionX()->Draw();
  //cq_IMnpipi_wK0_n_K0nn->cd(2);
  TH1D* q_K0nn = (TH1D*) q_IMnpipi_wK0_n_K0nn->ProjectionY();
  //q_K0nn->Draw();
  
  //TCanvas* cMompippim_IMnpipi_dE_wK0_n_K0nn = new TCanvas("cMompippim_IMnpipi_dE_wK0_n_K0nn","cMompippim_IMnpipi_dE_wK0_n_K0nn",1000,800);
  //cMompippim_IMnpipi_dE_wK0_n_K0nn->Divide(2,2);
  //cMompippim_IMnpipi_dE_wK0_n_K0nn->cd(1);
  //Mompippim_IMnpipi_dE_wK0_n_K0nn->Draw("colz");
  //cMompippim_IMnpipi_dE_wK0_n_K0nn->cd(3);
  //Mompippim_IMnpipi_dE_wK0_n_K0nn->ProjectionX()->Draw();
  //cMompippim_IMnpipi_dE_wK0_n_K0nn->cd(2);
  TH1D* Mompippim_K0nn = (TH1D*) Mompippim_IMnpipi_dE_wK0_n_K0nn->ProjectionY();
  //Mompippim_K0nn->Draw();
  TH1D* MMnmiss_K0nn = (TH1D*) MMnmiss_IMpippim_dE_K0nn->ProjectionY("MMnmiss_K0nn",pipimin,pipimax);

  //K-p -> K0n (n_s) -> K0-n, K-n ->K-n (p_s),  K-p->K0 n (two step reaction)
  TFile *_file4 = TFile::Open("simIMpisigma_K0_Knts_pippimn_v3_outncutK015.root");
  //TFile *_file4 = TFile::Open("simIMpisigma_K0_Knts_pippimn_v1_outncutK015.root");
  TH2F* q_IMpippim_n_Knts = (TH2F*)_file4->Get("q_IMpippim_n");
  TH2F* nmom_IMnpipi_wK0_n_Knts = (TH2F*)_file4->Get("nmom_IMnpipi_wK0_n");
  TH2F* q_IMnpipi_wK0_n_Knts = (TH2F*)_file4->Get("q_IMnpipi_wK0_n");
  TH2F* Mompippim_IMnpipi_dE_wK0_n_Knts = (TH2F*)_file4->Get("Mompippim_IMnpipi_dE_wK0_n");
  TH2F* MMnmiss_IMpippim_dE_Knts = (TH2F*)_file4->Get("MMnmiss_IMpippim_dE");

  //TCanvas *cnmom_IMnpipi_wK0_n_Knts = new TCanvas("cnmom_IMnpipi_wK0_n_Knts","nmom_IMnpipi_wK0_n_Knts",1000,800);
  //cnmom_IMnpipi_wK0_n_Knts->Divide(2,2);
  //cnmom_IMnpipi_wK0_n_Knts->cd(1);
  //nmom_IMnpipi_wK0_n_Knts->Draw("colz");
  //cnmom_IMnpipi_wK0_n_Knts->cd(3);
  TH1D* IMnpipi_Knts = (TH1D*) nmom_IMnpipi_wK0_n_Knts->ProjectionX();
  //IMnpipi_Knts->Draw();
  //cnmom_IMnpipi_wK0_n_Knts->cd(2);
  TH1D* nmom_IMnpipi_wK0_n_Knts_py = (TH1D*)nmom_IMnpipi_wK0_n_Knts->ProjectionY("nmom_IMnpipi_wK0_n_Knts_py");
  //nmom_IMnpipi_wK0_n_Knts_py->Draw();
  
  //TCanvas *cq_IMnpipi_wK0_n_Knts = new TCanvas("cq_IMnpipi_wK0_n_Knts","q_IMnpipi_wK0_n_Knts",1000,800);
  //cq_IMnpipi_wK0_n_Knts->Divide(2,2);
  //cq_IMnpipi_wK0_n_Knts->cd(1);
  //q_IMnpipi_wK0_n_Knts->Draw("colz");
  //fkp->Draw("same");
  //cq_IMnpipi_wK0_n_Knts->cd(3);
  //q_IMnpipi_wK0_n_Knts->ProjectionX()->Draw();
  //cq_IMnpipi_wK0_n_Knts->cd(2);
  TH1D* q_Knts = q_IMnpipi_wK0_n_Knts->ProjectionY();
  //q_Knts->Draw();
  

  //TCanvas* cMompippim_IMnpipi_dE_wK0_n_Knts = new TCanvas("cMompippim_IMnpipi_dE_wK0_n_Knts","cMompippim_IMnpipi_dE_wK0_n_Knts",1000,800);
  //cMompippim_IMnpipi_dE_wK0_n_Knts->Divide(2,2);
  //cMompippim_IMnpipi_dE_wK0_n_Knts->cd(1);
  //Mompippim_IMnpipi_dE_wK0_n_Knts->Draw("colz");
  //cMompippim_IMnpipi_dE_wK0_n_Knts->cd(3);
  //Mompippim_IMnpipi_dE_wK0_n_Knts->ProjectionX()->Draw();
  //cMompippim_IMnpipi_dE_wK0_n_Knts->cd(2);
  TH1D* Mompippim_Knts = (TH1D*)Mompippim_IMnpipi_dE_wK0_n_Knts->ProjectionY();
  //Mompippim_Knts->Draw();
  TH1D* MMnmiss_Knts = (TH1D*) MMnmiss_IMpippim_dE_Knts->ProjectionY("MMnmiss_Knts",pipimin,pipimax);
  
 


  const double scale_ns = 0.32;
  const double scale_nnts = 0.01;
  const double scale_Knts = 0.12;
  const double scale_K0nn = 0.01;

  //neutron momentum dist.
  TCanvas *cMom_ncds = new TCanvas("cMom_ncds","cMom_ncds");
  cMom_ncds->cd();
  nmom_IMnpipi_wK0_n_rdata_py->Draw("HE");
  nmom_IMnpipi_wK0_n_ns_py->SetLineColor(2);
  nmom_IMnpipi_wK0_n_ns_py->Scale(scale_ns);
  nmom_IMnpipi_wK0_n_ns_py->Draw("HEsame");
  nmom_IMnpipi_wK0_n_nnts_py->SetLineColor(3);
  nmom_IMnpipi_wK0_n_nnts_py->Scale(scale_nnts);
  nmom_IMnpipi_wK0_n_nnts_py->Draw("HEsame");
  nmom_IMnpipi_wK0_n_K0nn_py->Scale(scale_K0nn);
  nmom_IMnpipi_wK0_n_K0nn_py->SetLineColor(4);
  nmom_IMnpipi_wK0_n_K0nn_py->Draw("HEsame");
  nmom_IMnpipi_wK0_n_Knts_py->Scale(scale_Knts);
  nmom_IMnpipi_wK0_n_Knts_py->SetLineColor(5);
  nmom_IMnpipi_wK0_n_Knts_py->Draw("HEsame");
  TH1D* nmom_sum = (TH1D*)nmom_IMnpipi_wK0_n_ns_py->Clone();
  nmom_sum->Add(nmom_IMnpipi_wK0_n_nnts_py);
  nmom_sum->Add(nmom_IMnpipi_wK0_n_K0nn_py);
  nmom_sum->Add(nmom_IMnpipi_wK0_n_Knts_py);
  nmom_sum->SetLineColor(6);
  nmom_sum->Draw("same");

  TLegend *tl = new TLegend(0.8,0.68,0.99,0.99);
  tl->AddEntry(nmom_IMnpipi_wK0_n_rdata_py, "real data","l");
  tl->AddEntry(nmom_IMnpipi_wK0_n_ns_py, "1NA","l");
  tl->AddEntry(nmom_IMnpipi_wK0_n_nnts_py, "two step n-n scat.","l");
  tl->AddEntry(nmom_IMnpipi_wK0_n_K0nn_py, "K0nn phase space","l");
  tl->AddEntry(nmom_IMnpipi_wK0_n_Knts_py, "two step K0-n scat.","l");
  tl->Draw();

  //K0 momentum dist.
  TCanvas *cMom_K0 = new TCanvas("cMom_K0","cMom_K0");
  cMom_K0->cd();
  Mompippim_rdata->Draw("HE");
  Mompippim_ns->SetLineColor(2);
  Mompippim_ns->Scale(scale_ns);
  Mompippim_ns->Draw("HEsame");
  Mompippim_nnts->Scale(scale_nnts);
  Mompippim_nnts->SetLineColor(3);
  Mompippim_nnts->Draw("HEsame");
  Mompippim_K0nn->Scale(scale_K0nn);
  Mompippim_K0nn->SetLineColor(4);
  Mompippim_K0nn->Draw("HEsame");
  Mompippim_Knts->Scale(scale_Knts);
  Mompippim_Knts->SetLineColor(5);
  Mompippim_Knts->Draw("HEsame");
  TH1D* Mompippim_sum = (TH1D*)Mompippim_ns->Clone();
  Mompippim_sum->Add(Mompippim_nnts);
  Mompippim_sum->Add(Mompippim_K0nn);
  Mompippim_sum->Add(Mompippim_Knts);
  Mompippim_sum->SetLineColor(6);
  Mompippim_sum->Draw("HEsame");

  //IMnpipi
  TCanvas *cIMnpipi = new TCanvas("cIMnpipi","cIMnpipi");
  cIMnpipi->cd();
  IMnpipi_rdata->Draw("HE");
  IMnpipi_ns->SetLineColor(2);
  IMnpipi_ns->Scale(scale_ns);
  IMnpipi_ns->Draw("HEsame");
  IMnpipi_nnts->Scale(scale_nnts);
  IMnpipi_nnts->SetLineColor(3);
  IMnpipi_nnts->Draw("HEsame");
  IMnpipi_K0nn->Scale(scale_K0nn);
  IMnpipi_K0nn->SetLineColor(4);
  IMnpipi_K0nn->Draw("HEsame");
  IMnpipi_Knts->Scale(scale_Knts);
  IMnpipi_Knts->SetLineColor(5);
  IMnpipi_Knts->Draw("HEsame");
  TH1D* IMnpipi_sum = (TH1D*)IMnpipi_ns->Clone();
  IMnpipi_sum->Add(IMnpipi_nnts);
  IMnpipi_sum->Add(IMnpipi_K0nn);
  IMnpipi_sum->Add(IMnpipi_Knts);
  IMnpipi_sum->SetLineColor(6);
  IMnpipi_sum->Draw("HEsame");

  //q
  TCanvas *cq = new TCanvas("cq","cq");
  cq->cd();
  q_rdata->Draw("HE");
  q_ns->SetLineColor(2);
  q_ns->Scale(scale_ns);
  q_ns->Draw("HEsame");
  q_nnts->Scale(scale_nnts);
  q_nnts->SetLineColor(3);
  q_nnts->Draw("HEsame");
  q_K0nn->Scale(scale_K0nn);
  q_K0nn->SetLineColor(4);
  q_K0nn->Draw("HEsame");
  q_Knts->Scale(scale_Knts);
  q_Knts->SetLineColor(5);
  q_Knts->Draw("HEsame");
  TH1D* q_sum = (TH1D*)q_ns->Clone();
  q_sum->Add(q_Knts);
  q_sum->Add(q_nnts);
  q_sum->Add(q_K0nn);
  q_sum->SetLineColor(6);
  q_sum->Draw("HEsame");
  
  TCanvas *cnmiss = new TCanvas("cnmiss","cnmiss");
  cnmiss->cd();
  MMnmiss_rdata->Draw("HE");
  MMnmiss_ns->SetLineColor(2);
  MMnmiss_ns->Scale(scale_ns);
  MMnmiss_ns->Draw("HEsame");
  MMnmiss_nnts->SetLineColor(3);
  MMnmiss_nnts->Scale(scale_nnts);
  MMnmiss_nnts->Draw("HEsame");
  MMnmiss_K0nn->SetLineColor(4);
  MMnmiss_K0nn->Scale(scale_K0nn);
  MMnmiss_K0nn->Draw("HEsame");
  MMnmiss_Knts->SetLineColor(5);
  MMnmiss_Knts->Scale(scale_Knts);
  MMnmiss_Knts->Draw("HEsame");
  
  TH1D* MMnmiss_sum = (TH1D*)MMnmiss_ns->Clone();
  MMnmiss_sum->Add(MMnmiss_Knts);
  MMnmiss_sum->Add(MMnmiss_nnts);
  MMnmiss_sum->Add(MMnmiss_K0nn);
  MMnmiss_sum->SetLineColor(6);
  MMnmiss_sum->Draw("HEsame");

  return;
}
