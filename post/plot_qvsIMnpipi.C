
#include "../src/GlobalVariables.h"

void plot_qvsIMnpipi()
{
  TFile *file0 = TFile::Open("evanaIMpisigma_npippim_v195_out.root","READ");
  //TFile *file0 = TFile::Open("evanaIMpisigma_npippim_v195_out_iso.root","READ");
  file0->cd();
  TH2F* q_IMnpipi_wSid_n_Sp = (TH2F*)file0->Get("q_IMnpipi_wSid_n_Sp");
  TH2F* q_IMnpipi_wSid_n_Sm = (TH2F*)file0->Get("q_IMnpipi_wSid_n_Sm");
  TH2F* q_IMnpipi_wK0_wSid_n_Sp = (TH2F*)file0->Get("q_IMnpipi_wK0_wSid_n_Sp");
  TH2F* q_IMnpipi_wK0_wSid_n_Sm = (TH2F*)file0->Get("q_IMnpipi_wK0_wSid_n_Sm");

  TFile *fnu = TFile::Open("../simpost/NumericalRootFinder_K0n.root");

  TCanvas *c1 = new TCanvas("c1","c1");
  c1->cd();
  q_IMnpipi_wSid_n_Sp->Draw("colz");

  TCanvas *c2 = new TCanvas("c2","c2");
  c2->cd();
  q_IMnpipi_wSid_n_Sm->Draw("colz");

  TCanvas *c3 = new TCanvas("c3","c3");
  c3->cd();
  q_IMnpipi_wK0_wSid_n_Sp->Draw("colz");
  
  fnu->cd();
  TMultiGraph *mg = (TMultiGraph*)fnu->Get("mg"); 
  mg->Draw("c");

  TCanvas *c4 = new TCanvas("c4","c4");
  c4->cd();
  q_IMnpipi_wK0_wSid_n_Sm->Draw("colz");
  mg->Draw("c");

  TCanvas *c5 = new TCanvas("c5","c5");
  c5->cd();
  TH1D* IMnpipi_Sp = q_IMnpipi_wSid_n_Sp->ProjectionX("IMnpipi_Sp");
  IMnpipi_Sp->SetLineColor(2);
  IMnpipi_Sp->Draw("HE");
  TH1D* IMnpipi_Sm = q_IMnpipi_wSid_n_Sm->ProjectionX("IMnpipi_Sm");
  IMnpipi_Sm->SetLineColor(3);
  IMnpipi_Sm->Draw("HEsame");
  

  TCanvas *c6 = new TCanvas("c6","c6");
  c6->cd();
  const int binq300 = q_IMnpipi_wSid_n_Sp->GetYaxis()->FindBin(0.3);
  TH1D* IMnpipi_Sp_0_300 = q_IMnpipi_wSid_n_Sp->ProjectionX("IMnpipi_Sp_0_300",0,binq300-1);
  IMnpipi_Sp_0_300->SetLineColor(2);
  IMnpipi_Sp_0_300->Draw("HE");
  TH1D* IMnpipi_Sm_0_300 = q_IMnpipi_wSid_n_Sm->ProjectionX("IMnpipi_Sm_0_300",0,binq300-1);
  IMnpipi_Sm_0_300->SetLineColor(3);
  IMnpipi_Sm_0_300->Draw("HEsame");
  
  /*
  TCanvas *c7 = new TCanvas("c7","c7");
  c7->cd();
  q_IMnpipi_wSid_n_Sp->GetYaxis()->SetRange(0,binq300-1);
  q_IMnpipi_wSid_n_Sp->Draw("colz");
  */
  
  TCanvas *c7 = new TCanvas("c7","c7");
  c7->cd();
  TH1D* IMnpipi_Sp_300 = q_IMnpipi_wSid_n_Sp->ProjectionX("IMnpipi_Sp_300",binq300,100);
  IMnpipi_Sp_300->SetLineColor(2);
  IMnpipi_Sp_300->Draw("HE");
  TH1D* IMnpipi_Sm_300 = q_IMnpipi_wSid_n_Sm->ProjectionX("IMnpipi_Sm_300",binq300,100);
  IMnpipi_Sm_300->SetLineColor(3);
  IMnpipi_Sm_300->Draw("HEsame");

  
  TCanvas *c8 = new TCanvas("c8","c8");
  c8->cd();
  TH1D* IMnpipi_wK0_Sp_0_300 = q_IMnpipi_wK0_wSid_n_Sp->ProjectionX("IMnpipi_wK0_Sp_0_300",0,binq300-1);
  IMnpipi_wK0_Sp_0_300->SetLineColor(4);
  IMnpipi_wK0_Sp_0_300->Draw("HE");
  TH1D* IMnpipi_wK0_Sm_0_300 = q_IMnpipi_wK0_wSid_n_Sm->ProjectionX("IMnpipi_wK0_Sm_0_300",0,binq300-1);
  IMnpipi_wK0_Sm_0_300->SetLineColor(5);
  IMnpipi_wK0_Sm_0_300->Draw("HEsame");

  TCanvas *c9 = new TCanvas("c9","c9");
  c9->cd();
  TH1D* IMnpipi_wK0_Sp_300 = q_IMnpipi_wK0_wSid_n_Sp->ProjectionX("IMnpipi_wK0_Sp_300",binq300,100);
  IMnpipi_wK0_Sp_300->SetLineColor(4);
  IMnpipi_wK0_Sp_300->Draw("HE");
  TH1D* IMnpipi_wK0_Sm_300 = q_IMnpipi_wK0_wSid_n_Sm->ProjectionX("IMnpipi_wK0_Sm_300",binq300,100);
  IMnpipi_wK0_Sm_300->SetLineColor(5);
  IMnpipi_wK0_Sm_300->Draw("HEsame");


  TCanvas *c10 = new TCanvas("c10","c10");
  c10->cd();
  IMnpipi_Sp_0_300->Draw("HE");
  IMnpipi_wK0_Sp_0_300->Draw("HEsame");

  TCanvas *c11 = new TCanvas("c11","c11");
  c11->cd();
  IMnpipi_Sm_0_300->Draw("HE");
  IMnpipi_wK0_Sm_0_300->Draw("HEsame");

  TCanvas *c12 = new TCanvas("c12","c12");
  c12->cd();
  IMnpipi_Sp_300->Draw("HE");
  IMnpipi_wK0_Sp_300->Draw("HEsame");

  TCanvas *c13 = new TCanvas("c13","c13");
  c13->cd();
  IMnpipi_Sm_300->Draw("HE");
  IMnpipi_wK0_Sm_300->Draw("HEsame");
  
  TFile *fileK0n_ns = TFile::Open("../simpost/simIMpisigma_K0n_ns_pippimn_v24_out.root");
  fileK0n_ns->cd();
  TH2F* q_IMnpipi_wSid_n_Sp_K0n_ns_sim = (TH2F*)fileK0n_ns->Get("q_IMnpipi_wSid_n_Sp");
  TH2F* q_IMnpipi_wSid_n_Sm_K0n_ns_sim = (TH2F*)fileK0n_ns->Get("q_IMnpipi_wSid_n_Sm");
  TH2F* q_IMnpipi_wK0_wSid_n_Sp_K0n_ns_sim = (TH2F*)fileK0n_ns->Get("q_IMnpipi_wK0_wSid_n_Sp");
  TH2F* q_IMnpipi_wK0_wSid_n_Sm_K0n_ns_sim = (TH2F*)fileK0n_ns->Get("q_IMnpipi_wK0_wSid_n_Sm");
  
  const double cs_ns = 6.454; //elementary CS [mb] 
  const double scale_1NA = 0.06;
  q_IMnpipi_wSid_n_Sp_K0n_ns_sim->Scale(cs_ns*scale_1NA);
  q_IMnpipi_wSid_n_Sm_K0n_ns_sim->Scale(cs_ns*scale_1NA); 
  q_IMnpipi_wK0_wSid_n_Sp_K0n_ns_sim->Scale(cs_ns*scale_1NA);   
  q_IMnpipi_wK0_wSid_n_Sm_K0n_ns_sim->Scale(cs_ns*scale_1NA);   
  
  
  TCanvas *c14 = new TCanvas("c14","c14");
  q_IMnpipi_wK0_wSid_n_Sp_K0n_ns_sim->Draw("colz");
  const double Kp_mass = pMass + kpMass;  
  TF1 *fkp = new TF1("f", "sqrt(((x*x-[0]*[0]-[1]*[1])/(2*[0]))*((x*x-[0]*[0]-[1]*[1])/(2*[0]))-[1]*[1])",Kp_mass-0.0001,2);
  fkp->SetParameter(0,nMass);
  fkp->SetParameter(1,kpMass);
  //fkp->SetLineColor(4);
  fkp->SetLineWidth(5);
  fkp->SetLineStyle(4);
  fkp->SetLineColorAlpha(kPink, 0.35);
  fkp->Draw("same");
  TCanvas *c15 = new TCanvas("c15","c15");
  q_IMnpipi_wK0_wSid_n_Sm_K0n_ns_sim->Draw("colz");
  fkp->Draw("same");
  


  TCanvas *c16 = new TCanvas("c16","c16");
  c16->cd();
  TH1D* IMnpipi_Sp_K0n_ns_sim = q_IMnpipi_wSid_n_Sp_K0n_ns_sim->ProjectionX("IMnpipi_Sp_K0n_ns_sim");
  IMnpipi_Sp_K0n_ns_sim->SetLineColor(2);
  IMnpipi_Sp_K0n_ns_sim->Draw("HE");
  TH1D* IMnpipi_Sm_K0n_ns_sim = q_IMnpipi_wSid_n_Sm_K0n_ns_sim->ProjectionX("IMnpipi_Sm_K0n_ns_sim");
  IMnpipi_Sm_K0n_ns_sim->SetLineColor(3);
  IMnpipi_Sm_K0n_ns_sim->Draw("HEsame");
  
  TCanvas *c17 = new TCanvas("c17","c17");
  c17->cd();
  TH1D* IMnpipi_Sp_K0n_ns_sim_0_300 = q_IMnpipi_wSid_n_Sp_K0n_ns_sim->ProjectionX("IMnpipi_Sp_K0n_ns_sim_0_300",0,binq300-1);
  IMnpipi_Sp_K0n_ns_sim_0_300->SetLineColor(2);
  IMnpipi_Sp_K0n_ns_sim_0_300->Draw("HE");
  TH1D* IMnpipi_Sm_K0n_ns_sim_0_300 = q_IMnpipi_wSid_n_Sm_K0n_ns_sim->ProjectionX("IMnpipi_Sm_K0n_ns_sim_0_300",0,binq300-1);
  IMnpipi_Sm_K0n_ns_sim_0_300->SetLineColor(3);
  IMnpipi_Sm_K0n_ns_sim_0_300->Draw("HEsame");

  TCanvas *c18 = new TCanvas("c18","c18");
  c18->cd();
  TH1D* IMnpipi_Sp_K0n_ns_sim_300 = q_IMnpipi_wSid_n_Sp_K0n_ns_sim->ProjectionX("IMnpipi_Sp_K0n_ns_sim_300",binq300,100);
  IMnpipi_Sp_K0n_ns_sim_300->SetLineColor(2);
  IMnpipi_Sp_K0n_ns_sim_300->Draw("HE");
  TH1D* IMnpipi_Sm_K0n_ns_sim_300 = q_IMnpipi_wSid_n_Sm_K0n_ns_sim->ProjectionX("IMnpipi_Sm_K0n_ns_sim_300",binq300,100);
  IMnpipi_Sm_K0n_ns_sim_300->SetLineColor(3);
  IMnpipi_Sm_K0n_ns_sim_300->Draw("HEsame");

  

}
