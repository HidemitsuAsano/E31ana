#include "anacuts.h"

void disp_mulcut()
{
  TFile *fbefore = new TFile("evanaIMpisigma_npippim_v202_out.root");//no CDC 3 lay cut
  TFile *fafter  = new TFile("evanaIMpisigma_npippim_v202_out.root");//w CDC 3lay cut

  TH2F* MMnmiss_IMpippim_dE_wSid_b = (TH2F*)fbefore->Get("MMnmiss_IMpippim_dE_wSid");
  TH2F* MMnmiss_IMpippim_dE_wSid_a = (TH2F*)fafter->Get("MMnmiss_IMpippim_dE_wSid");


  TH1D* MMnmiss_b = (TH1D*)MMnmiss_IMpippim_dE_wSid_b->ProjectionY("MMnmiss_b");
  TH1D* MMnmiss_a = (TH1D*)MMnmiss_IMpippim_dE_wSid_a->ProjectionY("MMnmiss_a");
  
  TCanvas *c1 = new TCanvas();
  MMnmiss_b->Draw("HE");
  MMnmiss_a->SetLineColor(2);
  MMnmiss_a->Draw("HEsame");
  TBox *box = new TBox(anacuts::neutron_MIN,0,anacuts::neutron_MAX,2000);
  box->SetFillColor(4);
  box->SetFillStyle(3002);
  box->Draw();






}
