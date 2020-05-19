//plot_miss2D
//purpose : plot Missing mass, IM(npip) and so on. BG fit by modifying MC data.
//Mar.19th , version 2 clean up version

#include "TH1.h"
#include "TH2.h"
#include "TH2D.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TFile.h"
#include "TCanvas.h"
#include <fstream>

void comp_fakedata_v2()
{
  TH1::SetDefaultSumw2();
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  //rdata,nSp,nSm,Lambda,S0pi+pi-,fake,mcsum
  const unsigned int colordef[3]= {1,2,3};
  const char name[][10]= {"rdata","woK0","wK0"};
  //const char opt[10]="cont1z";
  const char opt[10]="colz";

  //real data
  TFile *filerdata = TFile::Open("../post/evanaIMpisigma_npippim_v196_out.root","READ");
  //fake nonK0+fake_n 
  TFile *filefake = TFile::Open("/gpfs/group/had/knucl/e15/asano/sim/fakemc/fakepippim_pippimn_out_sum.root","READ");
  //TFile *filefake = TFile::Open("fakepippim_pippimn_out_sum.root","READ");
  //fake K0(mass is smeared by resolution) + fake_n
  TFile *filefakeK0 = TFile::Open("/gpfs/group/had/knucl/e15/asano/sim/fakemc/fakepippimK0_pippimn_out_sum.root","READ");
  
  //superimpose 2d hists of miss mass vs IM(npi+) or IM(npi-)
  //fake
  //Sp
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sp_fake = (TH2D*) filefake->Get("MMnmiss_IMnpipi_woK0_wSid_Sp");
  double fake_Sp = MMnmiss_IMnpipi_woK0_wSid_Sp_fake->GetEntries();
  //Sm
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sm_fake = (TH2D*) filefake->Get("MMnmiss_IMnpipi_woK0_wSid_Sm");
  double fake_Sm = MMnmiss_IMnpipi_woK0_wSid_Sm_fake->GetEntries();
  
  TH2D *MMnmiss_IMnpipi_woK0_wSid_fake = (TH2D*)filefake->Get("MMnmiss_IMnpipi_woK0_wSid");
  TH2D *MMnmiss_IMnpipi_wK0_wSid_fake = (TH2D*)filefakeK0->Get("MMnmiss_IMnpipi_wK0_wSid");
  TH2D *MMnmiss_IMnpip_woK0_fake = (TH2D*)filefake->Get("MMnmiss_IMnpip_dE_woK0");
  TH2D *MMnmiss_IMnpim_woK0_fake = (TH2D*)filefake->Get("MMnmiss_IMnpim_dE_woK0");
  TH2D *MMnmiss_IMnpip_woK0_woSm_fake = (TH2D*)filefake->Get("MMnmiss_IMnpip_dE_woK0_woSm");
  TH2D *MMnmiss_IMnpip_wK0_woSm_fake = (TH2D*)filefakeK0->Get("MMnmiss_IMnpip_dE_wK0_woSm");
  TH2D *MMnmiss_IMnpim_woK0_woSp_fake = (TH2D*)filefake->Get("MMnmiss_IMnpim_dE_woK0_woSp");
  TH2D *MMnmiss_IMnpim_wK0_woSp_fake = (TH2D*)filefakeK0->Get("MMnmiss_IMnpim_dE_wK0_woSp");
  TH2D *MMnmiss_IMnpip_woK0_woSm_n_fake = (TH2D*)filefake->Get("MMnmiss_IMnpip_dE_woK0_woSm_n");
  TH2D *MMnmiss_IMnpip_wK0_woSm_n_fake = (TH2D*)filefakeK0->Get("MMnmiss_IMnpip_dE_wK0_woSm_n");
  TH2D *MMnmiss_IMnpim_woK0_woSp_n_fake = (TH2D*)filefake->Get("MMnmiss_IMnpim_dE_woK0_woSp_n");
  TH2D *MMnmiss_IMnpim_wK0_woSp_n_fake = (TH2D*)filefakeK0->Get("MMnmiss_IMnpim_dE_wK0_woSp_n");
  TH2D *MMnmiss_IMnpip_woK0_woSid_won_fake = (TH2D*)filefake->Get("MMnmiss_IMnpip_dE_woK0_woSid_won");
  TH2D *MMnmiss_IMnpim_woK0_woSid_won_fake = (TH2D*)filefake->Get("MMnmiss_IMnpim_dE_woK0_woSid_won");
  TH2D *MMnmiss_IMnpip_wK0_woSid_won_fake = (TH2D*)filefakeK0->Get("MMnmiss_IMnpip_dE_wK0_woSid_won");
  TH2D *MMnmiss_IMnpim_wK0_woSid_won_fake = (TH2D*)filefakeK0->Get("MMnmiss_IMnpim_dE_wK0_woSid_won");
  TH2D *MMnmiss_Momnpip_woK0_woSid_won_fake = (TH2D*)filefake->Get("MMnmiss_Momnpip_dE_woK0_woSid_won");
  TH2D *MMnmiss_Momnpim_woK0_woSid_won_fake = (TH2D*)filefake->Get("MMnmiss_Momnpim_dE_woK0_woSid_won");
  TH2D *MMnmiss_Momnpip_wK0_woSid_won_fake = (TH2D*)filefakeK0->Get("MMnmiss_Momnpip_dE_wK0_woSid_won");
  TH2D *MMnmiss_Momnpim_wK0_woSid_won_fake = (TH2D*)filefakeK0->Get("MMnmiss_Momnpim_dE_wK0_woSid_won");
  TH2D *MMnmiss_IMpippim_fake = (TH2D*)filefake->Get("MMnmiss_IMpippim_dE");
  TH2D *MMnmiss_IMpippim_wSid_fake = (TH2D*)filefake->Get("MMnmiss_IMpippim_dE_woK0_wSid");
  TH2D *MMnmiss_IMpippim_wSid_n_fake = (TH2D*)filefake->Get("MMnmiss_IMpippim_dE_woK0_wSid_n");
  TH2D *MMnmiss_IMpippim_wSid_fakeK0 = (TH2D*)filefakeK0->Get("MMnmiss_IMpippim_dE_wK0_wSid");
  TH2D *MMnmiss_IMpippim_wSid_n_fakeK0 = (TH2D*)filefakeK0->Get("MMnmiss_IMpippim_dE_wK0_wSid_n");
  TH2D *MMnmiss_IMpippim_woK0_woSid_won_fake = (TH2D*)filefake->Get("MMnmiss_IMpippim_dE_woK0_woSid_won");
  TH2D *MMnmiss_IMpippim_wK0_woSid_won_fake = (TH2D*)filefakeK0->Get("MMnmiss_IMpippim_dE_wK0_woSid_won");
  TH2D *IMnpim_IMnpip_fake = (TH2D*)filefake->Get("IMnpim_IMnpip_dE");
  TH2D *IMnpim_IMnpip_woK0_fake = (TH2D*)filefake->Get("IMnpim_IMnpip_dE_woK0");
  TH2D *IMnpim_IMnpip_wK0_fakeK0 = (TH2D*)filefakeK0->Get("IMnpim_IMnpip_dE_wK0");
  TH2D *IMnpim_IMnpip_woK0_n_fake = (TH2D*)filefake->Get("IMnpim_IMnpip_dE_woK0_n");
  TH2D *IMnpim_IMnpip_wK0_n_fakeK0 = (TH2D*)filefakeK0->Get("IMnpim_IMnpip_dE_wK0_n");
  TH2D *IMnpim_IMnpip_woK0_woSid_won_fake = (TH2D*)filefake->Get("IMnpim_IMnpip_dE_woK0_woSid_won");
  TH2D *IMnpim_IMnpip_wK0_woSid_won_fake = (TH2D*)filefakeK0->Get("IMnpim_IMnpip_dE_wK0_woSid_won");
  TH2D *IMpippim_IMnpip_woK0_woSid_won_fake = (TH2D*)filefake->Get("IMpippim_IMnpip_woK0_woSid_won");
  TH2D *IMpippim_IMnpip_wK0_woSid_won_fake = (TH2D*)filefakeK0->Get("IMpippim_IMnpip_wK0_woSid_won");
  TH2D *IMpippim_IMnpim_woK0_woSid_won_fake = (TH2D*)filefake->Get("IMpippim_IMnpim_woK0_woSid_won");
  TH2D *IMpippim_IMnpim_wK0_woSid_won_fake = (TH2D*)filefakeK0->Get("IMpippim_IMnpim_wK0_woSid_won");
  TH2D *q_IMnpipi_woK0_woSid_won_fake = (TH2D*)filefake->Get("q_IMnpipi_woK0_woSid_won");
  TH2D *q_IMnpipi_wK0_woSid_won_fake = (TH2D*)filefakeK0->Get("q_IMnpipi_wK0_woSid_won");
  TH2D *MMnmiss_Mompippim_woK0_woSid_won_fake  = (TH2D*)filefake->Get("MMnmiss_Mompippim_dE_woK0_woSid_won");
  TH2D *MMnmiss_Mompippim_wK0_woSid_won_fake = (TH2D*)filefakeK0->Get("MMnmiss_Mompippim_dE_wK0_woSid_won");
  TH2D *pipmom_MMnmiss_woK0_woSid_won_fake = (TH2D*)filefake->Get("pipmom_MMnmiss_dE_woK0_woSid_won");
  TH2D *pipmom_MMnmiss_wK0_woSid_won_fake = (TH2D*)filefakeK0->Get("pipmom_MMnmiss_dE_wK0_woSid_won");
  TH2D *pipmom_pimmom_woK0_woSid_won_fake = (TH2D*)filefake->Get("pipmom_pimmom_dE_woK0_woSid_won");
  TH2D *pipmom_pimmom_wK0_woSid_won_fake = (TH2D*)filefakeK0->Get("pipmom_pimmom_dE_wK0_woSid_won");
  TH2D *pimmom_MMnmiss_woK0_woSid_won_fake = (TH2D*)filefake->Get("pimmom_MMnmiss_dE_woK0_woSid_won");
  TH2D *pimmom_MMnmiss_wK0_woSid_won_fake = (TH2D*)filefakeK0->Get("pimmom_MMnmiss_dE_wK0_woSid_won");
  TH2D *MMom_MMass_woK0_woSid_won_fake = (TH2D*)filefake->Get("MMom_MMass_woK0_woSid_won");
  TH2D *MMom_MMass_wK0_woSid_won_fake = (TH2D*)filefakeK0->Get("MMom_MMass_wK0_woSid_won");
  TH2D *nmom_MMnmiss_woK0_woSid_won_fake = (TH2D*)filefake->Get("nmom_MMnmiss_woK0_woSid_won");
  TH2D *nmom_MMnmiss_wK0_woSid_won_fake = (TH2D*)filefakeK0->Get("nmom_MMnmiss_wK0_woSid_won");
  TH2D *nmom_IMpippim_woK0_woSid_won_fake = (TH2D*)filefake->Get("nmom_IMpippim_woK0_woSid_won");
  TH2D *nmom_IMpippim_wK0_woSid_won_fake = (TH2D*)filefakeK0->Get("nmom_IMpippim_wK0_woSid_won");
  TH2D *nmom_IMnpip_woK0_woSid_won_fake = (TH2D*)filefake->Get("nmom_IMnpip_woK0_woSid_won");
  TH2D *nmom_IMnpip_wK0_woSid_won_fake = (TH2D*)filefakeK0->Get("nmom_IMnpip_wK0_woSid_won");
  TH2D *nmom_IMnpim_woK0_woSid_won_fake = (TH2D*)filefake->Get("nmom_IMnpim_woK0_woSid_won");
  TH2D *nmom_IMnpim_wK0_woSid_won_fake = (TH2D*)filefakeK0->Get("nmom_IMnpim_wK0_woSid_won");
  TH2D *q_IMnpipi_woK0_wSid_n_fake = (TH2D*)filefake->Get("q_IMnpipi_woK0_wSid_n");
  TH2D *q_IMnpipi_wK0_wSid_n_fake = (TH2D*)filefakeK0->Get("q_IMnpipi_wK0_wSid_n");
  TH2D *q_IMnpipi_wSid_n_Sp_fake = (TH2D*)filefake->Get("q_IMnpipi_wSid_n_Sp");
  TH2D *q_IMnpipi_woK0_wSid_n_Sp_fake = (TH2D*)filefake->Get("q_IMnpipi_woK0_wSid_n_Sp");
  TH2D *q_IMnpipi_wK0_wSid_n_Sp_fake = (TH2D*)filefakeK0->Get("q_IMnpipi_wK0_wSid_n_Sp");
  TH2D *q_IMnpipi_wSid_n_Sm_fake = (TH2D*)filefake->Get("q_IMnpipi_wSid_n_Sm");
  TH2D *q_IMnpipi_woK0_wSid_n_Sm_fake = (TH2D*)filefake->Get("q_IMnpipi_woK0_wSid_n_Sm");
  TH2D *q_IMnpipi_wK0_wSid_n_Sm_fake = (TH2D*)filefakeK0->Get("q_IMnpipi_wK0_wSid_n_Sm");
  TH2D *q_MMnmiss_woK0_woSid_won_fake = (TH2D*)filefake->Get("q_MMnmiss_woK0_woSid_won");
  TH2D *q_MMnmiss_wK0_woSid_won_fake = (TH2D*)filefakeK0->Get("q_MMnmiss_wK0_woSid_won");
  TH2D *q_IMnpip_woK0_woSid_won_fake = (TH2D*)filefake->Get("q_IMnpip_woK0_woSid_won");
  TH2D *q_IMnpip_wK0_woSid_won_fake = (TH2D*)filefakeK0->Get("q_IMnpip_wK0_woSid_won");
  TH2D *q_IMnpim_woK0_woSid_won_fake = (TH2D*)filefake->Get("q_IMnpim_woK0_woSid_won");
  TH2D *q_IMnpim_wK0_woSid_won_fake = (TH2D*)filefakeK0->Get("q_IMnpim_wK0_woSid_won");
  TH2D *q_IMpippim_woK0_woSid_won_fake = (TH2D*)filefake->Get("q_IMpippim_woK0_woSid_won");
  TH2D *q_IMpippim_wK0_woSid_won_fake = (TH2D*)filefakeK0->Get("q_IMpippim_wK0_woSid_won");
  TH2D *q_nmom_woK0_woSid_won_fake = (TH2D*)filefake->Get("q_nmom_woK0_woSid_won");
  TH2D *q_nmom_wK0_woSid_won_fake = (TH2D*)filefakeK0->Get("q_nmom_wK0_woSid_won");
  TH2D *nmom_cosn_woK0_woSid_won_fake = (TH2D*)filefake->Get("nmom_cosn_woK0_woSid_won");
  TH2D *nmom_cospip_woK0_woSid_won_fake = (TH2D*)filefake->Get("nmom_cospip_woK0_woSid_won");
  TH2D *nmom_cospim_woK0_woSid_won_fake = (TH2D*)filefake->Get("nmom_cospim_woK0_woSid_won");
  TH2D *nmom_cospippim_woK0_woSid_won_fake = (TH2D*)filefake->Get("nmom_cospippim_woK0_woSid_won");
  TH2D *nmom_pipmom_woK0_woSid_won_fake = (TH2D*)filefake->Get("nmom_pipmom_woK0_woSid_won");
  TH2D *nmom_pimmom_woK0_woSid_won_fake = (TH2D*)filefake->Get("nmom_pimmom_woK0_woSid_won");
  TH2D *nmom_cosn_wK0_woSid_won_fake = (TH2D*)filefakeK0->Get("nmom_cosn_wK0_woSid_won");
  TH2D *nmom_cospip_wK0_woSid_won_fake = (TH2D*)filefakeK0->Get("nmom_cospip_wK0_woSid_won");
  TH2D *nmom_cospim_wK0_woSid_won_fake = (TH2D*)filefakeK0->Get("nmom_cospim_wK0_woSid_won");
  TH2D *nmom_cospippim_wK0_woSid_won_fake = (TH2D*)filefakeK0->Get("nmom_cospippim_wK0_woSid_won");
  TH2D *nmom_pipmom_wK0_woSid_won_fake = (TH2D*)filefakeK0->Get("nmom_pipmom_wK0_woSid_won");
  TH2D *nmom_pimmom_wK0_woSid_won_fake = (TH2D*)filefakeK0->Get("nmom_pimmom_wK0_woSid_won");
  TH2D *MMnmiss_Momnpipi_woK0_woSid_won_fake = (TH2D*)filefake->Get("MMnmiss_Momnpipi_woK0_woSid_won");
  TH2D *MMnmiss_Momnpipi_wK0_woSid_won_fake = (TH2D*)filefakeK0->Get("MMnmiss_Momnpipi_wK0_woSid_won");
  double Nnpip_fake = MMnmiss_IMnpip_woK0_fake->GetEntries();
  double Nnpim_fake = MMnmiss_IMnpim_woK0_fake->GetEntries();
  //real data
  
  //Sp
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sp_rdata = (TH2D*) filerdata->Get("MMnmiss_IMnpipi_woK0_wSid_Sp");
  double  nrdata_Sp = MMnmiss_IMnpipi_woK0_wSid_Sp_rdata->GetEntries();
  std::cout << nrdata_Sp << std::endl;
  //Sm
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sm_rdata = (TH2D*) filerdata->Get("MMnmiss_IMnpipi_woK0_wSid_Sm");
  double nrdata_Sm = MMnmiss_IMnpipi_woK0_wSid_Sm_rdata->GetEntries();
  TH2D *MMnmiss_IMnpipi_woK0_wSid_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpipi_woK0_wSid");
  TH2D *MMnmiss_IMnpipi_wK0_wSid_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpipi_wK0_wSid");
  TH2D *MMnmiss_IMnpip_woK0_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpip_dE_woK0");
  TH2D *MMnmiss_IMnpim_woK0_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpim_dE_woK0");
  TH2D *MMnmiss_IMnpip_woK0_woSm_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpip_dE_woK0_woSm");
  TH2D *MMnmiss_IMnpip_wK0_woSm_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpip_dE_wK0_woSm");
  TH2D *MMnmiss_IMnpim_woK0_woSp_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpim_dE_woK0_woSp");
  TH2D *MMnmiss_IMnpim_wK0_woSp_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpim_dE_wK0_woSp");
  TH2D *MMnmiss_IMnpip_woK0_woSm_n_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpip_dE_woK0_woSm_n");
  TH2D *MMnmiss_IMnpim_woK0_woSp_n_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpim_dE_woK0_woSp_n");
  TH2D *MMnmiss_IMnpip_wK0_woSm_n_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpip_dE_wK0_woSm_n");
  TH2D *MMnmiss_IMnpim_wK0_woSp_n_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpim_dE_wK0_woSp_n");
  TH2D *MMnmiss_IMnpip_woK0_woSid_won_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpip_dE_woK0_woSid_won");
  TH2D *MMnmiss_IMnpim_woK0_woSid_won_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpim_dE_woK0_woSid_won");
  TH2D *MMnmiss_IMnpip_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpip_dE_wK0_woSid_won");
  TH2D *MMnmiss_IMnpim_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpim_dE_wK0_woSid_won");
  TH2D *MMnmiss_Momnpip_woK0_woSid_won_rdata = (TH2D*)filerdata->Get("MMnmiss_Momnpip_dE_woK0_woSid_won");
  TH2D *MMnmiss_Momnpim_woK0_woSid_won_rdata = (TH2D*)filerdata->Get("MMnmiss_Momnpim_dE_woK0_woSid_won");
  TH2D *MMnmiss_Momnpip_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("MMnmiss_Momnpip_dE_wK0_woSid_won");
  TH2D *MMnmiss_Momnpim_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("MMnmiss_Momnpim_dE_wK0_woSid_won");
  TH2D *MMnmiss_IMpippim_rdata = (TH2D*)filerdata->Get("MMnmiss_IMpippim_dE");
  TH2D *MMnmiss_IMpippim_wSid_rdata = (TH2D*)filerdata->Get("MMnmiss_IMpippim_dE_wSid");
  TH2D *MMnmiss_IMpippim_wSid_n_rdata = (TH2D*)filerdata->Get("MMnmiss_IMpippim_dE_wSid_n");
  TH2D *MMnmiss_IMpippim_woK0_woSid_won_rdata = (TH2D*)filerdata->Get("MMnmiss_IMpippim_dE_woK0_woSid_won");
  TH2D *MMnmiss_IMpippim_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("MMnmiss_IMpippim_dE_wK0_woSid_won");
  TH2D *IMnpim_IMnpip_rdata = (TH2D*)filerdata->Get("IMnpim_IMnpip_dE");
  TH2D *IMnpim_IMnpip_woK0_rdata = (TH2D*)filerdata->Get("IMnpim_IMnpip_dE_woK0");
  TH2D *IMnpim_IMnpip_wK0_rdata = (TH2D*)filerdata->Get("IMnpim_IMnpip_dE_wK0");
  //TH2D *IMnpim_IMnpip_n_rdata = (TH2D*)filerdata->Get("IMnpim_IMnpip_dE_n");
  TH2D *IMnpim_IMnpip_woK0_n_rdata = (TH2D*)filerdata->Get("IMnpim_IMnpip_dE_woK0_n");
  TH2D *IMnpim_IMnpip_wK0_n_rdata = (TH2D*)filerdata->Get("IMnpim_IMnpip_dE_wK0_n");
  TH2D *IMnpim_IMnpip_wSid_n_rdata = (TH2D*)filerdata->Get("IMnpim_IMnpip_dE_wSid_n");
  TH2D *IMnpim_IMnpip_woK0_woSid_won_rdata = (TH2D*)filerdata->Get("IMnpim_IMnpip_dE_woK0_woSid_won");
  TH2D *IMnpim_IMnpip_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("IMnpim_IMnpip_dE_wK0_woSid_won");
  TH2D *IMpippim_IMnpip_woK0_woSid_won_rdata = (TH2D*)filerdata->Get("IMpippim_IMnpip_woK0_woSid_won");
  TH2D *IMpippim_IMnpip_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("IMpippim_IMnpip_wK0_woSid_won");
  TH2D *IMpippim_IMnpim_woK0_woSid_won_rdata = (TH2D*)filerdata->Get("IMpippim_IMnpim_woK0_woSid_won");
  TH2D *IMpippim_IMnpim_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("IMpippim_IMnpim_wK0_woSid_won");
  TH2D *q_IMnpipi_woK0_woSid_won_rdata = (TH2D*)filerdata->Get("q_IMnpipi_woK0_woSid_won");
  TH2D *q_IMnpipi_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("q_IMnpipi_wK0_woSid_won");
  TH2D *MMnmiss_Mompippim_woK0_woSid_won_rdata  = (TH2D*)filerdata->Get("MMnmiss_Mompippim_dE_woK0_woSid_won");
  TH2D *MMnmiss_Mompippim_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("MMnmiss_Mompippim_dE_wK0_woSid_won");
  TH2D *pipmom_MMnmiss_woK0_woSid_won_rdata = (TH2D*)filerdata->Get("pipmom_MMnmiss_dE_woK0_woSid_won");
  TH2D *pipmom_MMnmiss_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("pipmom_MMnmiss_dE_wK0_woSid_won");
  TH2D *pipmom_pimmom_woK0_woSid_won_rdata = (TH2D*)filerdata->Get("pipmom_pimmom_dE_woK0_woSid_won");
  TH2D *pipmom_pimmom_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("pipmom_pimmom_dE_wK0_woSid_won");
  TH2D *pimmom_MMnmiss_woK0_woSid_won_rdata = (TH2D*)filerdata->Get("pimmom_MMnmiss_dE_woK0_woSid_won");
  TH2D *pimmom_MMnmiss_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("pimmom_MMnmiss_dE_wK0_woSid_won");
  TH2D *MMom_MMass_woK0_woSid_won_rdata = (TH2D*)filerdata->Get("MMom_MMass_woK0_woSid_won");
  TH2D *MMom_MMass_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("MMom_MMass_wK0_woSid_won");
  TH2D *nmom_MMnmiss_woK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_MMnmiss_woK0_woSid_won");
  TH2D *nmom_MMnmiss_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_MMnmiss_wK0_woSid_won");
  TH2D *nmom_IMpippim_woK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_IMpippim_woK0_woSid_won");
  TH2D *nmom_IMpippim_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_IMpippim_wK0_woSid_won");
  TH2D *nmom_IMnpip_woK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_IMnpip_woK0_woSid_won");
  TH2D *nmom_IMnpip_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_IMnpip_wK0_woSid_won");
  TH2D *nmom_IMnpim_woK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_IMnpim_woK0_woSid_won");
  TH2D *nmom_IMnpim_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_IMnpim_wK0_woSid_won");
  TH2D *q_IMnpipi_wSid_n_rdata = (TH2D*)filerdata->Get("q_IMnpipi_wSid_n");
  TH2D *q_IMnpipi_wSid_n_Sp_rdata = (TH2D*)filerdata->Get("q_IMnpipi_wSid_n_Sp");
  TH2D *q_IMnpipi_wSid_n_Sm_rdata = (TH2D*)filerdata->Get("q_IMnpipi_wSid_n_Sm");
  TH2D *q_IMnpip_woK0_woSid_won_rdata = (TH2D*)filerdata->Get("q_IMnpip_woK0_woSid_won");
  TH2D *q_IMnpip_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("q_IMnpip_wK0_woSid_won");
  TH2D *q_IMnpim_woK0_woSid_won_rdata = (TH2D*)filerdata->Get("q_IMnpim_woK0_woSid_won");
  TH2D *q_IMnpim_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("q_IMnpim_wK0_woSid_won");
  TH2D *q_IMpippim_woK0_woSid_won_rdata = (TH2D*)filerdata->Get("q_IMpippim_woK0_woSid_won");
  TH2D *q_IMpippim_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("q_IMpippim_wK0_woSid_won");
  TH2D *q_nmom_woK0_woSid_won_rdata = (TH2D*)filerdata->Get("q_nmom_woK0_woSid_won");
  TH2D *q_nmom_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("q_nmom_wK0_woSid_won");
  TH2D *nmom_cosn_woK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_cosn_woK0_woSid_won");
  TH2D *nmom_cospip_woK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_cospip_woK0_woSid_won");
  TH2D *nmom_cospim_woK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_cospim_woK0_woSid_won");
  TH2D *nmom_cospippim_woK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_cospippim_woK0_woSid_won");
  TH2D *nmom_pipmom_woK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_pipmom_woK0_woSid_won");
  TH2D *nmom_pimmom_woK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_pimmom_woK0_woSid_won");
  TH2D *nmom_cosn_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_cosn_wK0_woSid_won");
  TH2D *nmom_cospip_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_cospip_wK0_woSid_won");
  TH2D *nmom_cospim_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_cospim_wK0_woSid_won");
  TH2D *nmom_cospippim_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_cospippim_wK0_woSid_won");
  TH2D *nmom_pipmom_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_pipmom_wK0_woSid_won");
  TH2D *nmom_pimmom_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_pimmom_wK0_woSid_won");
  TH2D *MMnmiss_Momnpipi_woK0_woSid_won_rdata = (TH2D*)filerdata->Get("MMnmiss_Momnpipi_woK0_woSid_won");
  TH2D *MMnmiss_Momnpipi_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("MMnmiss_Momnpipi_wK0_woSid_won");
  double Nnpip_rdata = MMnmiss_IMnpip_woK0_rdata->GetEntries();
  double Nnpim_rdata = MMnmiss_IMnpim_woK0_rdata->GetEntries();
   
  //v2 2020113
  const double scaleFactor[3]= {1.0,
                                0.3*nrdata_Sp/fake_Sp,
                                0.3*nrdata_Sp/fake_Sp
                               };


  //arrays
  TH2D* MMnmiss_IMnpipi_woK0_wSid[2]= {
    MMnmiss_IMnpipi_woK0_wSid_rdata,
    MMnmiss_IMnpipi_woK0_wSid_fake
  };
  
  TH2D* MMnmiss_IMnpipi_wK0_wSid[2]= {
    MMnmiss_IMnpipi_wK0_wSid_rdata,
    MMnmiss_IMnpipi_wK0_wSid_fake
  };

  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sp[2]= {
    MMnmiss_IMnpipi_woK0_wSid_Sp_rdata,
    MMnmiss_IMnpipi_woK0_wSid_Sp_fake
  };

  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sm[2]= {
    MMnmiss_IMnpipi_woK0_wSid_Sm_rdata,
    MMnmiss_IMnpipi_woK0_wSid_Sm_fake
  };

  TH2D* MMnmiss_IMnpip_woK0[2]= {
    MMnmiss_IMnpip_woK0_rdata,
    MMnmiss_IMnpip_woK0_fake
  };

  TH2D* MMnmiss_IMnpip_woK0_woSm[2]= {
    MMnmiss_IMnpip_woK0_woSm_rdata,
    MMnmiss_IMnpip_woK0_woSm_fake
  };

  TH2D* MMnmiss_IMnpip_woK0_woSm_n[2]= {
    MMnmiss_IMnpip_woK0_woSm_n_rdata,
    MMnmiss_IMnpip_woK0_woSm_n_fake
  };
  
  TH2D* MMnmiss_IMnpip_wK0_woSm_n[2]= {
    MMnmiss_IMnpip_wK0_woSm_n_rdata,
    MMnmiss_IMnpip_wK0_woSm_n_fake
  };

  TH2D* MMnmiss_IMnpip_woK0_woSid_won[2]= {
    MMnmiss_IMnpip_woK0_woSid_won_rdata,
    MMnmiss_IMnpip_woK0_woSid_won_fake
  };

  TH2D* MMnmiss_IMnpip_wK0_woSid_won[2]= {
    MMnmiss_IMnpip_wK0_woSid_won_rdata,
    MMnmiss_IMnpip_wK0_woSid_won_fake
  };
  
  TH2D* MMnmiss_Momnpip_woK0_woSid_won[2]= {
    MMnmiss_Momnpip_woK0_woSid_won_rdata,
    MMnmiss_Momnpip_woK0_woSid_won_fake
  };

  TH2D* MMnmiss_Momnpip_wK0_woSid_won[2]= {
    MMnmiss_Momnpip_wK0_woSid_won_rdata,
    MMnmiss_Momnpip_wK0_woSid_won_fake
  };

  TH2D* MMnmiss_IMnpim_woK0[2]= {
    MMnmiss_IMnpim_woK0_rdata,
    MMnmiss_IMnpim_woK0_fake
  };

  TH2D* MMnmiss_IMnpim_woK0_woSm[2]= {
    MMnmiss_IMnpim_woK0_woSp_rdata,
    MMnmiss_IMnpim_woK0_woSp_fake
  };
  
  TH2D* MMnmiss_IMnpim_woK0_woSp[2]= {
    MMnmiss_IMnpim_woK0_woSp_rdata,
    MMnmiss_IMnpim_woK0_woSp_fake
  };

  TH2D* MMnmiss_IMnpim_woK0_woSp_n[2]= {
    MMnmiss_IMnpim_woK0_woSp_n_rdata,
    MMnmiss_IMnpim_woK0_woSp_n_fake
  };
  
  TH2D* MMnmiss_IMnpim_wK0_woSp_n[2]= {
    MMnmiss_IMnpim_wK0_woSp_n_rdata,
    MMnmiss_IMnpim_wK0_woSp_n_fake
  };

  TH2D* MMnmiss_IMnpim_woK0_woSid_won[2]= {
    MMnmiss_IMnpim_woK0_woSid_won_rdata,
    MMnmiss_IMnpim_woK0_woSid_won_fake
  };

  TH2D* MMnmiss_IMnpim_wK0_woSid_won[2]= {
    MMnmiss_IMnpim_wK0_woSid_won_rdata,
    MMnmiss_IMnpim_wK0_woSid_won_fake
  };
  
  TH2D* MMnmiss_Momnpim_woK0_woSid_won[2]= {
    MMnmiss_Momnpim_woK0_woSid_won_rdata,
    MMnmiss_Momnpim_woK0_woSid_won_fake
  };

  TH2D* MMnmiss_Momnpim_wK0_woSid_won[2]= {
    MMnmiss_Momnpim_wK0_woSid_won_rdata,
    MMnmiss_Momnpim_wK0_woSid_won_fake
  };

  
  TH2D* IMnpim_IMnpip_woK0_woSid_won[2]= {
    IMnpim_IMnpip_woK0_woSid_won_rdata,
    IMnpim_IMnpip_woK0_woSid_won_fake
  };
  
  TH2D* IMnpim_IMnpip_wK0_woSid_won[2]= {
    IMnpim_IMnpip_wK0_woSid_won_rdata,
    IMnpim_IMnpip_wK0_woSid_won_fake
  };


  TH2D* IMpippim_IMnpip_woK0_woSid_won[2]={
    IMpippim_IMnpip_woK0_woSid_won_rdata,
    IMpippim_IMnpip_woK0_woSid_won_fake
  };
  
  TH2D* IMpippim_IMnpip_wK0_woSid_won[2]={
    IMpippim_IMnpip_wK0_woSid_won_rdata,
    IMpippim_IMnpip_wK0_woSid_won_fake
  };
  
  TH2D* IMpippim_IMnpim_woK0_woSid_won[2]={
    IMpippim_IMnpim_woK0_woSid_won_rdata,
    IMpippim_IMnpim_woK0_woSid_won_fake
  };
  
  TH2D* IMpippim_IMnpim_wK0_woSid_won[2]={
    IMpippim_IMnpim_wK0_woSid_won_rdata,
    IMpippim_IMnpim_wK0_woSid_won_fake
  };

  TH2D* MMnmiss_IMpippim_wSid[3]= {
    MMnmiss_IMpippim_wSid_rdata,
    MMnmiss_IMpippim_wSid_fake,
    MMnmiss_IMpippim_wSid_fakeK0
  };
  
  TH2D* MMnmiss_IMpippim_wSid_n[3]= {
    MMnmiss_IMpippim_wSid_n_rdata,
    MMnmiss_IMpippim_wSid_n_fake,
    MMnmiss_IMpippim_wSid_n_fakeK0
  };

  TH2D* MMnmiss_IMpippim_woK0_woSid_won[2]= {
    MMnmiss_IMpippim_woK0_woSid_won_rdata,
    MMnmiss_IMpippim_woK0_woSid_won_fake
  };

  TH2D *MMnmiss_IMpippim_wK0_woSid_won[2]= {
    MMnmiss_IMpippim_wK0_woSid_won_rdata,
    MMnmiss_IMpippim_wK0_woSid_won_fake
  };

  TH2D* IMnpim_IMnpip_woK0_n[2]={
    IMnpim_IMnpip_woK0_n_rdata,
    IMnpim_IMnpip_woK0_n_fake
  };
  
  TH2D* IMnpim_IMnpip_wK0_n[2]={
    IMnpim_IMnpip_wK0_n_rdata,
    IMnpim_IMnpip_wK0_n_fakeK0
  };

  TH2D* q_IMnpipi_woK0_woSid_won[2]= {
    q_IMnpipi_woK0_woSid_won_rdata,
    q_IMnpipi_woK0_woSid_won_fake
  };

  TH2D* q_IMnpipi_wK0_woSid_won[2]= {
    q_IMnpipi_wK0_woSid_won_rdata,
    q_IMnpipi_wK0_woSid_won_fake
  };

  
  TH2D* MMnmiss_Mompippim_woK0_woSid_won[2]= {
    MMnmiss_Mompippim_woK0_woSid_won_rdata,
    MMnmiss_Mompippim_woK0_woSid_won_fake
  };

  TH2D* MMnmiss_Mompippim_wK0_woSid_won[2]= {
    MMnmiss_Mompippim_wK0_woSid_won_rdata,
    MMnmiss_Mompippim_wK0_woSid_won_fake
  };

  
  TH2D* pipmom_MMnmiss_woK0_woSid_won[2]= {
    pipmom_MMnmiss_woK0_woSid_won_rdata,
    pipmom_MMnmiss_woK0_woSid_won_fake
  };

  TH2D* pipmom_MMnmiss_wK0_woSid_won[2]= {
    pipmom_MMnmiss_wK0_woSid_won_rdata,
    pipmom_MMnmiss_wK0_woSid_won_fake
  };


  TH2D *pipmom_pimmom_woK0_woSid_won[2]={
    pipmom_pimmom_woK0_woSid_won_rdata,
    pipmom_pimmom_woK0_woSid_won_fake
  };
  
  
  TH2D *pipmom_pimmom_wK0_woSid_won[2]={
    pipmom_pimmom_wK0_woSid_won_rdata,
    pipmom_pimmom_wK0_woSid_won_fake
  };

  
  TH2D* pimmom_MMnmiss_woK0_woSid_won[2]= {
    pimmom_MMnmiss_woK0_woSid_won_rdata,
    pimmom_MMnmiss_woK0_woSid_won_fake
  };

  TH2D* pimmom_MMnmiss_wK0_woSid_won[2]= {
    pimmom_MMnmiss_wK0_woSid_won_rdata,
    pimmom_MMnmiss_wK0_woSid_won_fake
  };
  
  TH2D* MMom_MMass_woK0_woSid_won[2]= {
    MMom_MMass_woK0_woSid_won_rdata,
    MMom_MMass_woK0_woSid_won_fake
  };

  TH2D* MMom_MMass_wK0_woSid_won[2]= {
    MMom_MMass_wK0_woSid_won_rdata,
    MMom_MMass_wK0_woSid_won_fake
  };

  
  TH2D* nmom_MMnmiss_woK0_woSid_won[2]= {
    nmom_MMnmiss_woK0_woSid_won_rdata,
    nmom_MMnmiss_woK0_woSid_won_fake
  };

  TH2D* nmom_MMnmiss_wK0_woSid_won[2]= {
    nmom_MMnmiss_wK0_woSid_won_rdata,
    nmom_MMnmiss_wK0_woSid_won_fake
  };

  TH2D* nmom_IMpippim_woK0_woSid_won[2]={
    nmom_IMpippim_woK0_woSid_won_rdata,
    nmom_IMpippim_woK0_woSid_won_fake
  };
  
  TH2D* nmom_IMpippim_wK0_woSid_won[2]={
    nmom_IMpippim_woK0_woSid_won_rdata,
    nmom_IMpippim_woK0_woSid_won_fake
  };

  TH2D* nmom_IMnpip_woK0_woSid_won[2]={
    nmom_IMnpip_woK0_woSid_won_rdata,
    nmom_IMnpip_woK0_woSid_won_fake
  };
  
  TH2D* nmom_IMnpip_wK0_woSid_won[2]={
    nmom_IMnpip_wK0_woSid_won_rdata,
    nmom_IMnpip_wK0_woSid_won_fake
  };

  TH2D* nmom_IMnpim_woK0_woSid_won[2]={
    nmom_IMnpim_woK0_woSid_won_rdata,
    nmom_IMnpim_woK0_woSid_won_fake
  };
  
  TH2D* nmom_IMnpim_wK0_woSid_won[2]={
    nmom_IMnpim_wK0_woSid_won_rdata,
    nmom_IMnpim_wK0_woSid_won_fake
  };


  TH2D* q_IMnpipi_wSid_n[3]= {
    q_IMnpipi_wSid_n_rdata,
    q_IMnpipi_woK0_wSid_n_fake,
    q_IMnpipi_wK0_wSid_n_fake
  };

  TH2D* q_IMnpipi_wSid_n_Sp[3]= {
    q_IMnpipi_wSid_n_Sp_rdata,
    q_IMnpipi_woK0_wSid_n_Sp_fake,
    q_IMnpipi_wK0_wSid_n_Sp_fake
  };

  TH2D* q_IMnpipi_wSid_n_Sm[3]= {
    q_IMnpipi_wSid_n_Sm_rdata,
    q_IMnpipi_woK0_wSid_n_Sm_fake,
    q_IMnpipi_wK0_wSid_n_Sm_fake
  };

  
  TH2D* nmom_cosn_woK0_woSid_won[2]= {
    nmom_cosn_woK0_woSid_won_rdata,
    nmom_cosn_woK0_woSid_won_fake
  };

  
  TH2D* nmom_cospip_woK0_woSid_won[2]= {
    nmom_cospip_woK0_woSid_won_rdata,
    nmom_cospip_woK0_woSid_won_fake
  };

  
  TH2D* nmom_cospim_woK0_woSid_won[2]= {
    nmom_cospim_woK0_woSid_won_rdata,
    nmom_cospim_woK0_woSid_won_fake
  };
  
  TH2D* nmom_cospippim_woK0_woSid_won[2]= {
    nmom_cospippim_woK0_woSid_won_rdata,
    nmom_cospippim_woK0_woSid_won_fake
  };

  
  TH2D* nmom_pipmom_woK0_woSid_won[2]= {
    nmom_pipmom_woK0_woSid_won_rdata,
    nmom_pipmom_woK0_woSid_won_fake
  };

  
  TH2D* nmom_pimmom_woK0_woSid_won[2]= {
    nmom_pimmom_woK0_woSid_won_rdata,
    nmom_pimmom_woK0_woSid_won_fake
  };

  TH2D* nmom_cosn_wK0_woSid_won[2]= {
    nmom_cosn_wK0_woSid_won_rdata,
    nmom_cosn_wK0_woSid_won_fake
  };

  TH2D* nmom_cospip_wK0_woSid_won[2]= {
    nmom_cospip_wK0_woSid_won_rdata,
    nmom_cospip_wK0_woSid_won_fake
  };

  TH2D* nmom_cospim_wK0_woSid_won[2]= {
    nmom_cospim_wK0_woSid_won_rdata,
    nmom_cospim_wK0_woSid_won_fake
  };
  
  TH2D* nmom_cospippim_wK0_woSid_won[2]= {
    nmom_cospippim_wK0_woSid_won_rdata,
    nmom_cospippim_wK0_woSid_won_fake
  };


  TH2D* nmom_pipmom_wK0_woSid_won[2]= {
    nmom_pipmom_wK0_woSid_won_rdata,
    nmom_pipmom_wK0_woSid_won_fake
  };

  TH2D* nmom_pimmom_wK0_woSid_won[2]= {
    nmom_pimmom_wK0_woSid_won_rdata,
    nmom_pimmom_wK0_woSid_won_fake
  };

  TH2D* MMnmiss_Momnpipi_woK0_woSid_won[2]={
    MMnmiss_Momnpipi_woK0_woSid_won_rdata,
    MMnmiss_Momnpipi_woK0_woSid_won_fake
  };
  
  TH2D* MMnmiss_Momnpipi_wK0_woSid_won[2]={
    MMnmiss_Momnpipi_wK0_woSid_won_rdata,
    MMnmiss_Momnpipi_wK0_woSid_won_fake
  };


  for(int i=0; i<2; i++) {
    MMnmiss_IMnpipi_woK0_wSid[i]->Scale(scaleFactor[i]);
    MMnmiss_IMnpipi_wK0_wSid[i]->Scale(scaleFactor[i]);
    MMnmiss_IMnpipi_woK0_wSid_Sp[i]->Scale(scaleFactor[i]);
    MMnmiss_IMnpip_woK0[i]->Scale(scaleFactor[i]);
    MMnmiss_IMnpip_woK0_woSm[i]->Scale(scaleFactor[i]);
    MMnmiss_IMnpip_woK0_woSm_n[i]->Scale(scaleFactor[i]);
    MMnmiss_IMnpip_wK0_woSm_n[i]->Scale(scaleFactor[i]);
    MMnmiss_IMnpip_woK0_woSid_won[i]->Scale(scaleFactor[i]);
    MMnmiss_IMnpip_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    MMnmiss_Momnpip_woK0_woSid_won[i]->Scale(scaleFactor[i]);
    MMnmiss_Momnpip_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    MMnmiss_IMnpipi_woK0_wSid_Sm[i]->Scale(scaleFactor[i]);
    MMnmiss_IMnpim_woK0[i]->Scale(scaleFactor[i]);
    MMnmiss_IMnpim_woK0_woSp[i]->Scale(scaleFactor[i]);
    MMnmiss_IMnpim_woK0_woSp_n[i]->Scale(scaleFactor[i]);
    MMnmiss_IMnpim_wK0_woSp_n[i]->Scale(scaleFactor[i]);
    MMnmiss_IMnpim_woK0_woSid_won[i]->Scale(scaleFactor[i]);
    MMnmiss_IMnpim_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    MMnmiss_Momnpim_woK0_woSid_won[i]->Scale(scaleFactor[i]);
    MMnmiss_Momnpim_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    IMnpim_IMnpip_woK0_woSid_won[i]->Scale(scaleFactor[i]);
    IMnpim_IMnpip_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    IMpippim_IMnpip_woK0_woSid_won[i]->Scale(scaleFactor[i]);
    IMpippim_IMnpip_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    IMpippim_IMnpim_woK0_woSid_won[i]->Scale(scaleFactor[i]);
    IMpippim_IMnpim_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    MMnmiss_IMpippim_woK0_woSid_won[i]->Scale(scaleFactor[i]);
    MMnmiss_IMpippim_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    IMnpim_IMnpip_woK0_n[i]->Scale(scaleFactor[i]);
    IMnpim_IMnpip_wK0_n[i]->Scale(scaleFactor[i]);
    q_IMnpipi_woK0_woSid_won[i]->Scale(scaleFactor[i]);
    q_IMnpipi_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    MMnmiss_Mompippim_woK0_woSid_won[i]->Scale(scaleFactor[i]);
    MMnmiss_Mompippim_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    pipmom_MMnmiss_woK0_woSid_won[i]->Scale(scaleFactor[i]);
    pipmom_MMnmiss_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    pipmom_pimmom_woK0_woSid_won[i]->Scale(scaleFactor[i]);
    pipmom_pimmom_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    pimmom_MMnmiss_woK0_woSid_won[i]->Scale(scaleFactor[i]);
    pimmom_MMnmiss_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    MMom_MMass_woK0_woSid_won[i]->Scale(scaleFactor[i]);
    MMom_MMass_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    nmom_MMnmiss_woK0_woSid_won[i]->Scale(scaleFactor[i]);
    nmom_MMnmiss_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    nmom_IMpippim_woK0_woSid_won[i]->Scale(scaleFactor[i]);
    nmom_IMpippim_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    nmom_IMnpip_woK0_woSid_won[i]->Scale(scaleFactor[i]);
    nmom_IMnpip_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    nmom_IMnpim_woK0_woSid_won[i]->Scale(scaleFactor[i]);
    nmom_IMnpim_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    nmom_cosn_woK0_woSid_won[i]->Scale(scaleFactor[i]);
    nmom_cospip_woK0_woSid_won[i]->Scale(scaleFactor[i]);
    nmom_cospim_woK0_woSid_won[i]->Scale(scaleFactor[i]);
    nmom_cospippim_woK0_woSid_won[i]->Scale(scaleFactor[i]);
    nmom_pipmom_woK0_woSid_won[i]->Scale(scaleFactor[i]);
    nmom_pimmom_woK0_woSid_won[i]->Scale(scaleFactor[i]);
    nmom_cosn_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    nmom_cospip_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    nmom_cospim_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    nmom_cospippim_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    nmom_pipmom_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    nmom_pimmom_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    MMnmiss_Momnpipi_woK0_woSid_won[i]->Scale(scaleFactor[i]);
    MMnmiss_Momnpipi_wK0_woSid_won[i]->Scale(scaleFactor[i]);
  }
  for(int i=0;i<3;i++){
    q_IMnpipi_wSid_n[i]->Scale(scaleFactor[i]);
    q_IMnpipi_wSid_n_Sp[i]->Scale(scaleFactor[i]);
    q_IMnpipi_wSid_n_Sm[i]->Scale(scaleFactor[i]);
    MMnmiss_IMpippim_wSid[i]->Scale(scaleFactor[i]);
    MMnmiss_IMpippim_wSid_n[i]->Scale(scaleFactor[i]);
  }
  
  TCanvas *cMMnmiss_woK0_wSid = new TCanvas("cMMnmiss_woK0_wSid","cMMnmiss_woK0_wSid",800,800);
  cMMnmiss_woK0_wSid->cd();
  TH1D* MMnmiss_woK0_wSid[2];
  for(int i=0; i<2; i++) MMnmiss_woK0_wSid[i] = (TH1D*)MMnmiss_IMnpipi_woK0_wSid[i]->ProjectionY(Form("MMnmiss_woK0_wSid_%s",name[i]));
  MMnmiss_woK0_wSid[0]->Draw("HE");
  TH1D* MMnmiss_woK0_wSid_mc = (TH1D*)MMnmiss_IMnpipi_woK0_wSid[1]->ProjectionY("MMnmiss_woK0_wSid_mc");
  MMnmiss_woK0_wSid_mc->SetLineColor(6);
  MMnmiss_woK0_wSid_mc->Draw("HEsame");
  
  

  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sp_mc = (TH2D*)MMnmiss_IMnpipi_woK0_wSid_Sp[1]->Clone("MMnmiss_IMnpipi_woK0_wSid_Sp_mc");

  TCanvas *cMMnmiss_IMnpipi_woK0_wSid_Sp = new TCanvas("cMMnmiss_IMnpipi_woK0_wSid_Sp","cMMnmiss_IMnpipi_woK0_wSid_Sp",1200,800);
  cMMnmiss_IMnpipi_woK0_wSid_Sp->Divide(2,1);
  cMMnmiss_IMnpipi_woK0_wSid_Sp->cd(1);
  MMnmiss_IMnpipi_woK0_wSid_Sp[0]->SetTitle("#splitline{MMnmiss_IMnpipi_woK0_wSid_Sp}{  real Data}");
  MMnmiss_IMnpipi_woK0_wSid_Sp[0]->SetMinimum(1);
  MMnmiss_IMnpipi_woK0_wSid_Sp[0]->Draw(opt);
  cMMnmiss_IMnpipi_woK0_wSid_Sp->cd(2);
  MMnmiss_IMnpipi_woK0_wSid_Sp_mc->SetTitle("#splitline{MMnmiss_IMnpipi_woK0_wSid_Sp}{  MC sum}");
  MMnmiss_IMnpipi_woK0_wSid_Sp_mc->SetMaximum(MMnmiss_IMnpipi_woK0_wSid_Sp[0]->GetMaximum());
  MMnmiss_IMnpipi_woK0_wSid_Sp_mc->SetMinimum(1);
  MMnmiss_IMnpipi_woK0_wSid_Sp_mc->Draw(opt);

  TCanvas *cMMnmiss_woK0_wSid_Sp = new TCanvas("cMMnmiss_woK0_wSid_Sp","cMMnmiss_woK0_wSid_Sp",800,800);
  cMMnmiss_woK0_wSid_Sp->cd();
  TH1D* MMnmiss_woK0_wSid_Sp[2];
  for(int i=0; i<2; i++) MMnmiss_woK0_wSid_Sp[i] = (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sp[i]->ProjectionX(Form("MMnmiss_woK0_wSid_Sp_%s",name[i]));
  MMnmiss_woK0_wSid_Sp[0]->Draw("HE");
  TH1D* MMnmiss_woK0_wSid_Sp_mc = (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sp_mc->ProjectionX("MMnmiss_woK0_wSid_Sp_mc");
  MMnmiss_woK0_wSid_Sp_mc->SetLineColor(6);
  MMnmiss_woK0_wSid_Sp_mc->Draw("HEsame");
  for(int i=0; i<2; i++) {
    MMnmiss_woK0_wSid_Sp[i]->SetLineColor(colordef[i]);
    //MMnmiss_woK0_wSid_Sp[i]->Draw("HEsame");
  }

  TCanvas *cIMnpim_IMnpip_n = new TCanvas("cIMnpim_IMnpip_n","cIMnpim_IMnpip_n");
  cIMnpim_IMnpip_n->Divide(2,1);
  cIMnpim_IMnpip_n->cd(1);
  TH2D *IMnpim_IMnpip_n_rdata = (TH2D*)IMnpim_IMnpip_woK0_n[0]->Clone("IMnpim_IMnpip_n_rdata");
  IMnpim_IMnpip_n_rdata->Add(IMnpim_IMnpip_wK0_n[0]);
  IMnpim_IMnpip_n_rdata->Draw(opt);
  cIMnpim_IMnpip_n->cd(2);
  TH2D *IMnpim_IMnpip_n_mc = (TH2D*)IMnpim_IMnpip_woK0_n[1]->Clone("IMnpim_IMnpip_n_mc");
  IMnpim_IMnpip_n_mc->Add(IMnpim_IMnpip_wK0_n[1]);
  IMnpim_IMnpip_n_mc->SetMaximum(IMnpim_IMnpip_n_rdata->GetMaximum());
  IMnpim_IMnpip_n_mc->SetMinimum(1);
  IMnpim_IMnpip_n_mc->Draw(opt);




  TCanvas *cIMnpipi_woK0_wSid_Sp = new TCanvas("cIMnpipi_woK0_wSid_Sp","cIMnpipi_woK0_wSid_Sp");
  cIMnpipi_woK0_wSid_Sp->cd();
  TH1D* IMnpipi_woK0_wSid_Sp[2];
  for(int i=0; i<2; i++)IMnpipi_woK0_wSid_Sp[i] = (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sp[i]->ProjectionY(Form("IMnpipi_woK0_wSid_Sp_%s",name[i]));
  IMnpipi_woK0_wSid_Sp[0]->Draw("HE");
  TH1D* IMnpipi_woK0_wSid_Sp_mc = (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sp_mc->ProjectionY("IMnpipi_woK0_wSid_Sp_mc");
  IMnpipi_woK0_wSid_Sp_mc->SetLineColor(6);
  IMnpipi_woK0_wSid_Sp_mc->Draw("HEsame");
  for(int i=1; i<2; i++) {
    IMnpipi_woK0_wSid_Sp[i]->SetLineColor(colordef[i]);
    //IMnpipi_woK0_wSid_Sp[i]->Draw("HEsame");
  }


  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sm_mc = (TH2D*)MMnmiss_IMnpipi_woK0_wSid_Sm[1]->Clone("MMnmiss_IMnpipi_woK0_wSid_Sm_mc");

  TCanvas *cMMnmiss_IMnpipi_woK0_wSid_Sm = new TCanvas("cMMnmiss_IMnpipi_woK0_wSid_Sm","cMMnmiss_IMnpipi_woK0_wSid_Sm",1200,800);
  cMMnmiss_IMnpipi_woK0_wSid_Sm->Divide(2,1);
  cMMnmiss_IMnpipi_woK0_wSid_Sm->cd(1);
  MMnmiss_IMnpipi_woK0_wSid_Sm[0]->SetTitle("#splitline{MMnmiss_IMnpipi_woK0_wSid_Sm}{  real data}");
  MMnmiss_IMnpipi_woK0_wSid_Sm[0]->SetMinimum(1);
  MMnmiss_IMnpipi_woK0_wSid_Sm[0]->Draw(opt);
  cMMnmiss_IMnpipi_woK0_wSid_Sm->cd(2);
  MMnmiss_IMnpipi_woK0_wSid_Sm_mc->SetTitle("#splitline{MMnmiss_IMnpipi_woK0_wSid_Sm}{  MC sum}");
  MMnmiss_IMnpipi_woK0_wSid_Sm_mc->SetMaximum(MMnmiss_IMnpipi_woK0_wSid_Sm[0]->GetMaximum());
  MMnmiss_IMnpipi_woK0_wSid_Sm_mc->SetMinimum(1);
  MMnmiss_IMnpipi_woK0_wSid_Sm_mc->Draw(opt);

  TCanvas *cMMnmiss_woK0_wSid_Sm = new TCanvas("cMMnmiss_woK0_wSid_Sm","cMMnmiss_woK0_wSid_Sm",800,800);
  cMMnmiss_woK0_wSid_Sm->cd();
  TH1D* MMnmiss_woK0_wSid_Sm[2];
  for(int i=0; i<2; i++) MMnmiss_woK0_wSid_Sm[i] = (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sm[i]->ProjectionX(Form("MMnmiss_woK0_wSid_Sm_%s",name[i]));
  MMnmiss_woK0_wSid_Sm[0]->Draw("HE");
  TH1D* MMnmiss_woK0_wSid_Sm_mc = (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sm_mc->ProjectionX("MMnmiss_woK0_wSid_Sm_mc");
  MMnmiss_woK0_wSid_Sm_mc->SetLineColor(6);
  MMnmiss_woK0_wSid_Sm_mc->Draw("HEsame");
  for(int i=1; i<2; i++) {
    MMnmiss_woK0_wSid_Sm[i]->SetLineColor(colordef[i]);
    //MMnmiss_woK0_wSid_Sm[i]->Draw("HEsame");
  }

  TCanvas *cIMnpipi_woK0_wSid_Sm  = new TCanvas("cIMnpipi_woK0_wSid_Sm","cIMnpipi_woK0_wSid_Sm");
  cIMnpipi_woK0_wSid_Sm->cd();
  TH1D* IMnpipi_woK0_wSid_Sm[2];
  for(int i=0; i<2; i++)IMnpipi_woK0_wSid_Sm[i] = (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sm[i]->ProjectionY(Form("IMnpipi_woK0_wSid_Sm_%s",name[i]));
  IMnpipi_woK0_wSid_Sm[0]->Draw("HE");
  TH1D* IMnpipi_woK0_wSid_Sm_mc = (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sm_mc->ProjectionY("IMnpipi_woK0_wSid_Sm_mc");
  IMnpipi_woK0_wSid_Sm_mc->SetLineColor(6);
  IMnpipi_woK0_wSid_Sm_mc->Draw("HEsame");
  for(int i=1; i<2; i++) {
    IMnpipi_woK0_wSid_Sm[i]->SetLineColor(colordef[i]);
    //IMnpipi_woK0_wSid_Sm[i]->Draw("HEsame");
  }


  //missing mass and IMnpip/npim w/ Sp or Sm selection w/o K0
  //
  //ploting Sp mode
  //

  //w/o K0, including Sp/Sm mode
  TH2D *MMnmiss_IMnpip_woK0_mc = (TH2D*)MMnmiss_IMnpip_woK0[1]->Clone("MMnmiss_IMnpip_woK0_mc");

  TCanvas *cMMnmiss_IMnpip_woK0 = new TCanvas("cMMnmiss_IMnpip_woK0","cMMnmiss_IMnpip_woK0",1200,800);
  cMMnmiss_IMnpip_woK0->Divide(2,1);
  cMMnmiss_IMnpip_woK0->cd(1);
  MMnmiss_IMnpip_woK0_rdata->SetTitle("#splitline{MMnmiss_IMnpip_woK0}{  Real data}");
  MMnmiss_IMnpip_woK0_rdata->SetMinimum(1);
  MMnmiss_IMnpip_woK0_rdata->Draw(opt);
  cMMnmiss_IMnpip_woK0->cd(2);
  MMnmiss_IMnpip_woK0_mc->SetTitle("#splitline{MMnmiss_IMnpip_woK0}{  MC sum}");
  MMnmiss_IMnpip_woK0_mc->SetMaximum(MMnmiss_IMnpip_woK0_rdata->GetMaximum());
  MMnmiss_IMnpip_woK0_mc->SetMinimum(1);
  MMnmiss_IMnpip_woK0_mc->Draw(opt);
  

  //w/o K0, w/o Sm mode
  TH2D *MMnmiss_IMnpip_woK0_woSm_mc = (TH2D*)MMnmiss_IMnpip_woK0_woSm[1]->Clone("MMnmiss_IMnpip_woK0_woSm_mc");

  //w/o K0, w/o Sm mode, selecting missing neutron
  TH2D *MMnmiss_IMnpip_woK0_woSm_n_mc = (TH2D*)MMnmiss_IMnpip_woK0_woSm_n[1]->Clone("MMnmiss_IMnpip_woK0_woSm_n_mc");

  TCanvas *cMMnmiss_IMnpip_woK0_woSm = new TCanvas("cMMnmiss_IMnpip_woK0_woSm","cMMnmiss_IMnpip_woK0_woSm",1200,800);
  cMMnmiss_IMnpip_woK0_woSm->Divide(2,1);
  cMMnmiss_IMnpip_woK0_woSm->cd(1);
  MMnmiss_IMnpip_woK0_woSm_rdata->SetTitle("#splitline{MMnmiss_IMnpip_woK0_woSm}{  Real data}");
  MMnmiss_IMnpip_woK0_woSm_rdata->SetMinimum(1);
  MMnmiss_IMnpip_woK0_woSm_rdata->Draw(opt);
  cMMnmiss_IMnpip_woK0_woSm->cd(2);
  MMnmiss_IMnpip_woK0_woSm_mc->SetTitle("#splitline{MMnmiss_IMnpip_woK0_woSm}{  MC sum}");
  MMnmiss_IMnpip_woK0_woSm_mc->SetMinimum(1);
  MMnmiss_IMnpip_woK0_woSm_mc->SetMaximum(MMnmiss_IMnpip_woK0_woSm_rdata->GetMaximum());
  MMnmiss_IMnpip_woK0_woSm_mc->Draw(opt);

  //TCanvas *cIMnpip_woK0_woSm_n = new TCanvas("cIMnpip_woK0_woSm_n","cIMnpip_woK0_woSm_n");
  //cIMnpip_woK0_woSm_n->cd();
  //TH1D *IMnpip_woK0_woSm_n[2];
  //for(int i=0; i<2; i++) {
  //  IMnpip_woK0_woSm_n[i] = (TH1D*)MMnmiss_IMnpip_woK0_woSm_n[i]->ProjectionX(Form("IMnpip_woK0_woSm_n_%s",name[i]));
  //  IMnpip_woK0_woSm_n[i]->SetLineColor(colordef[i]);
  //}
  //TH1D* IMnpip_woK0_woSm_n_mc = (TH1D*)MMnmiss_IMnpip_woK0_woSm_n_mc->ProjectionX("IMnpip_woK0_woSm_n_mc");
  //IMnpip_woK0_woSm_n[0]->Draw("HE");
  //IMnpip_woK0_woSm_n_mc->SetLineColor(6);
  //IMnpip_woK0_woSm_n_mc->Draw("same");
  


  //w/o K0, w/o ((Sp or Sm) & missing neutron)
  TH2D *MMnmiss_IMnpip_woK0_woSid_won_mc = (TH2D*)MMnmiss_IMnpip_woK0_woSid_won[1]->Clone("MMnmiss_IMnpip_woK0_woSid_won_mc");


  //projection to Missing mass (miss n & Sigma+/-)
  TCanvas *cMMnmiss_woK0_woSid_won = new TCanvas("cMMnmiss_woK0_woSid_won","cMMnmiss_woK0_woSid_won");
  cMMnmiss_woK0_woSid_won->cd();
  TH1D* MMnmiss_woK0_woSid_won[2];
  for(int i=0; i<2; i++) MMnmiss_woK0_woSid_won[i] = (TH1D*)MMnmiss_IMnpip_woK0_woSid_won[i]->ProjectionY(Form("MMnmiss_woK0_woSid_won_%s",name[i]));
  MMnmiss_woK0_woSid_won[0]->Draw("HE");
  TH1D* MMnmiss_woK0_woSid_won_mc = (TH1D*)MMnmiss_IMnpip_woK0_woSid_won_mc->ProjectionY("MMnmiss_woK0_woSid_won_mc");
  MMnmiss_woK0_woSid_won_mc->SetLineColor(6);
  MMnmiss_woK0_woSid_won_mc->Draw("HEsame");

  for(int i=1; i<2; i++) {
    MMnmiss_woK0_woSid_won[i]->SetLineColor(colordef[i]);
    //MMnmiss_woK0_woSid_won[i]->Draw("HEsame");
  }

  //Data/MC before modifying MC data
  //TCanvas *cMMnmiss_woK0_woSid_won_ratio = new TCanvas("cMMnmiss_woK0_woSid_won_ratio","cMMnmiss_woK0_woSid_ratio");
  TH1D* MMnmiss_woK0_woSid_won_ratio = (TH1D*)MMnmiss_woK0_woSid_won[0]->Clone("MMnmiss_woK0_woSid_won_ratio");
  MMnmiss_woK0_woSid_won_ratio->Divide(MMnmiss_woK0_woSid_won_mc);
  MMnmiss_woK0_woSid_won_ratio->SetTitle("Data/MC");
  MMnmiss_woK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  //MMnmiss_woK0_woSid_won_ratio->Draw("HE");


  TH2D *MMom_MMass_woK0_woSid_won_mc = (TH2D*)MMom_MMass_woK0_woSid_won[1]->Clone("MMom_MMass_woK0_woSid_won_mc");
  TCanvas *cMMom_MMass_woK0_woSid_won = new TCanvas("cMMom_MMass_woK0_woSid_won","cMMom_MMass_woK0_woSid_won");
  cMMom_MMass_woK0_woSid_won->Divide(2,1);
  cMMom_MMass_woK0_woSid_won->cd(1);
  MMom_MMass_woK0_woSid_won[0]->SetMinimum(1);
  MMom_MMass_woK0_woSid_won[0]->Draw(opt);
  cMMom_MMass_woK0_woSid_won->cd(2);
  MMom_MMass_woK0_woSid_won_mc->SetMinimum(1);
  MMom_MMass_woK0_woSid_won_mc->SetMaximum(MMom_MMass_woK0_woSid_won[0]->GetMaximum());
  MMom_MMass_woK0_woSid_won_mc->Draw(opt);

  TH1D* MMom_woK0_woSid_won_mc = MMom_MMass_woK0_woSid_won_mc->ProjectionY("MMom_woK0_woSid_won_mc");
  TH1D* MMom_woK0_woSid_won[2];
  for(int i=0; i<2; i++) {
    MMom_woK0_woSid_won[i] = (TH1D*)MMom_MMass_woK0_woSid_won[i]->ProjectionY(Form("MMom_woK0_woSid_won_%s",name[i]));
    MMom_woK0_woSid_won[i]->SetLineColor(colordef[i]);
  }
  TCanvas *cMMom_woK0_woSid_won = new TCanvas("cMMom_woK0_woSid_won","cMMom_woK0_woSid_won");
  MMom_woK0_woSid_won[0]->Draw("HE");
  MMom_woK0_woSid_won_mc->SetLineColor(6);
  MMom_woK0_woSid_won_mc->Draw("HEsame");
  //for(int i=1;i<6;i++)MMom_woK0_woSid_won[i]->Draw("HEsame");


  //TCanvas *cMMom_woK0_woSid_won_ratio = new TCanvas("cMMom_woK0_woSid_won_ratio","cMMom_woK0_woSid_won_ratio");
  //cMMom_woK0_woSid_won_ratio->cd();
  //TH1D* MMom_woK0_woSid_won_ratio = (TH1D*)MMom_woK0_woSid_won[0]->Clone("MMom_woK0_woSid_won_ratio");
  //MMom_woK0_woSid_won_ratio->Divide(MMom_woK0_woSid_won_mc);
  //MMom_woK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  //MMom_woK0_woSid_won_ratio->SetTitle("Data/MC");
  //MMom_woK0_woSid_won_ratio->Draw("HE");


  //projection to IMnpip (miss n & Sigma+/-)
  TCanvas *cIMnpip_woK0_woSid_won = new TCanvas("cIMnpip_woK0_woSid_won","cIMnpip_woK0_woSid_won");
  cIMnpip_woK0_woSid_won->cd();
  TH1D* IMnpip_woK0_woSid_won[6];
  for(int i=0; i<2; i++)IMnpip_woK0_woSid_won[i] = (TH1D*)MMnmiss_IMnpip_woK0_woSid_won[i]->ProjectionX(Form("IMnpip_woK0_woSid_won_%s",name[i]));
  IMnpip_woK0_woSid_won[0]->Draw("HE");//rdata
  TH1D* IMnpip_woK0_woSid_won_mc = (TH1D*)MMnmiss_IMnpip_woK0_woSid_won_mc->ProjectionX("IMnpip_woK0_woSid_won_mc");
  IMnpip_woK0_woSid_won_mc->SetLineColor(6);
  IMnpip_woK0_woSid_won_mc->Draw("HEsame");

  //TCanvas *cIMnpip_woK0_woSid_won_ratio = new TCanvas("cIMnpip_woK0_woSid_won_ratio","cIMnpip_woK0_woSid_won_ratio");
  //cIMnpip_woK0_woSid_won_ratio->cd();
  TH1D* IMnpip_woK0_woSid_won_ratio = (TH1D*)IMnpip_woK0_woSid_won[0]->Clone("IMnpip_woK0_woSid_won_ratio");
  IMnpip_woK0_woSid_won_ratio->Divide(IMnpip_woK0_woSid_won_mc);
  IMnpip_woK0_woSid_won_ratio->SetTitle("IMnpip_woK0_woSid_won Data/MC");
  IMnpip_woK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  //IMnpip_woK0_woSid_won_ratio->Draw("HE");
  
  //
  TCanvas *cMomnpip_woK0_woSid_won = new TCanvas("cMomnpip_woK0_woSid_won","cMomnpip_woK0_woSid_won");
  cMomnpip_woK0_woSid_won->cd();
  TH1D* Momnpip_woK0_woSid_won[6];
  for(int i=0; i<2; i++)Momnpip_woK0_woSid_won[i] = (TH1D*)MMnmiss_Momnpip_woK0_woSid_won[i]->ProjectionX(Form("Momnpip_woK0_woSid_won_%s",name[i]));
  Momnpip_woK0_woSid_won[0]->Draw("HE");//rdata
  TH1D* Momnpip_woK0_woSid_won_mc = (TH1D*)MMnmiss_Momnpip_woK0_woSid_won[1]->ProjectionX("Momnpip_woK0_woSid_won_mc");
  Momnpip_woK0_woSid_won_mc->SetLineColor(6);
  Momnpip_woK0_woSid_won_mc->Draw("HEsame");

  //TCanvas *cMomnpip_woK0_woSid_won_ratio = new TCanvas("cMomnpip_woK0_woSid_won_ratio","cMomnpip_woK0_woSid_won_ratio");
  //cMomnpip_woK0_woSid_won_ratio->cd();
  TH1D* Momnpip_woK0_woSid_won_ratio = (TH1D*)Momnpip_woK0_woSid_won[0]->Clone("Momnpip_woK0_woSid_won_ratio");
  Momnpip_woK0_woSid_won_ratio->Divide(Momnpip_woK0_woSid_won_mc);
  Momnpip_woK0_woSid_won_ratio->SetTitle("Momnpip_woK0_woSid_won Data/MC");
  Momnpip_woK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  //Momnpip_woK0_woSid_won_ratio->Draw("HE");

  //
  //Sigma-
  TH2D *MMnmiss_IMnpim_woK0_mc = (TH2D*)MMnmiss_IMnpim_woK0[1]->Clone("MMnmiss_IMnpim_woK0_mc");

  TCanvas *cMMnmiss_IMnpim_woK0 = new TCanvas("cMMnmiss_IMnpim_woK0","cMMnmiss_IMnpim_woK0",1200,800);
  cMMnmiss_IMnpim_woK0->Divide(2,1);
  cMMnmiss_IMnpim_woK0->cd(1);
  MMnmiss_IMnpim_woK0_rdata->SetTitle("#splitline{MMnmiss_IMnpim_woK0}{  Real data}");
  MMnmiss_IMnpim_woK0_rdata->SetMinimum(1);
  MMnmiss_IMnpim_woK0_rdata->Draw(opt);
  cMMnmiss_IMnpim_woK0->cd(2);
  MMnmiss_IMnpim_woK0_mc->SetTitle("#splitline{MMnmiss_IMnpim_woK0}{ MC sum}");
  MMnmiss_IMnpim_woK0_mc->SetMinimum(1);
  MMnmiss_IMnpim_woK0_mc->SetMaximum(MMnmiss_IMnpim_woK0_rdata->GetMaximum());
  MMnmiss_IMnpim_woK0_mc->Draw(opt);

  //w/o K0, w/o Sp mode
  TH2D *MMnmiss_IMnpim_woK0_woSp_mc = (TH2D*)MMnmiss_IMnpim_woK0_woSp[1]->Clone("MMnmiss_IMnpim_woK0_woSp_mc");
  //w/o K0, w/o Sp mode, selecting missing neutron
  TH2D *MMnmiss_IMnpim_woK0_woSp_n_mc = (TH2D*)MMnmiss_IMnpim_woK0_woSp_n[1]->Clone("MMnmiss_IMnpim_woK0_woSp_n_mc");

  TCanvas *cMMnmiss_IMnpim_woK0_woSp = new TCanvas("cMMnmiss_IMnpim_woK0_woSp","cMMnmiss_IMnpim_woK0_woSp",1200,800);
  cMMnmiss_IMnpim_woK0_woSp->Divide(2,1);
  cMMnmiss_IMnpim_woK0_woSp->cd(1);
  MMnmiss_IMnpim_woK0_woSp_rdata->SetTitle("#splitline{MMnmiss_IMnpim_woK0_woSp}{  Real data}");
  MMnmiss_IMnpim_woK0_woSp_rdata->SetMinimum(1);
  MMnmiss_IMnpim_woK0_woSp_rdata->Draw(opt);
  cMMnmiss_IMnpim_woK0_woSp->cd(2);
  MMnmiss_IMnpim_woK0_woSp_mc->SetTitle("#splitline{MMnmiss_IMnpim_woK0_woSp}{  MC sum}");
  MMnmiss_IMnpim_woK0_woSp_mc->SetMinimum(1);
  MMnmiss_IMnpim_woK0_woSp_mc->SetMaximum(MMnmiss_IMnpim_woK0_woSp_rdata->GetMaximum());
  MMnmiss_IMnpim_woK0_woSp_mc->Draw(opt);
  
  TCanvas *cIMnpim_woK0_woSp_n = new TCanvas("cIMnpim_woK0_woSp_n","cIMnpim_woK0_woSp_n");
  cIMnpim_woK0_woSp_n->cd();
  TH1D *IMnpim_woK0_woSp_n[2];
  for(int i=0; i<2; i++) {
    IMnpim_woK0_woSp_n[i] = (TH1D*)MMnmiss_IMnpim_woK0_woSp_n[i]->ProjectionX(Form("IMnpim_woK0_woSp_n_%s",name[i]));
  }
  TH1D *IMnpim_wK0_woSp_n[2];
  for(int i=0; i<2; i++) {
    IMnpim_wK0_woSp_n[i] = (TH1D*)MMnmiss_IMnpim_wK0_woSp_n[i]->ProjectionX(Form("IMnpim_wK0_woSp_n_%s",name[i]));
  }
  //TH1D* IMnpim_woK0_woSp_n_mc = (TH1D*) MMnmiss_IMnpim_woK0_woSp_n_mc->ProjectionX("IMnpim_woK0_woSp_n_mc");
  TH1D* IMnpim_woSp_n_rdata = (TH1D*)IMnpim_woK0_woSp_n[0]->Clone("IMnpim_woSp_n_rdata");
  IMnpim_woSp_n_rdata->Add(IMnpim_wK0_woSp_n[0]);
  IMnpim_woSp_n_rdata->Draw("HE");
  IMnpim_woK0_woSp_n[1]->SetLineColor(2);
  IMnpim_woK0_woSp_n[1]->Draw("HEsame");
  IMnpim_wK0_woSp_n[1]->SetLineColor(3);
  IMnpim_wK0_woSp_n[1]->Draw("HEsame");
  TH1D* IMnpim_woSp_n_mc = (TH1D*)IMnpim_woK0_woSp_n[1]->Clone("IMnpim_woSp_n_mc");
  IMnpim_woSp_n_mc->Add(IMnpim_wK0_woSp_n[1]);
  IMnpim_woSp_n_mc->SetLineColor(6);
  IMnpim_woSp_n_mc->Draw("HEsame");


  TCanvas *cIMnpip_woK0_woSm_n = new TCanvas("cIMnpip_woK0_woSm_n","cIMnpip_woK0_woSm_n");
  cIMnpip_woK0_woSm_n->cd();
  TH1D *IMnpip_woK0_woSm_n[2];
  for(int i=0; i<2; i++) {
    IMnpip_woK0_woSm_n[i] = (TH1D*)MMnmiss_IMnpip_woK0_woSm_n[i]->ProjectionX(Form("IMnpip_woK0_woSm_n_%s",name[i]));
  }
  TH1D *IMnpip_wK0_woSm_n[2];
  for(int i=0; i<2; i++) {
    IMnpip_wK0_woSm_n[i] = (TH1D*)MMnmiss_IMnpip_wK0_woSm_n[i]->ProjectionX(Form("IMnpip_wK0_woSm_n_%s",name[i]));
  }
  TH1D* IMnpip_woSm_n_rdata = (TH1D*)IMnpip_woK0_woSm_n[0]->Clone("IMnpip_woSm_n_rdata");
  IMnpip_woSm_n_rdata->Add(IMnpip_wK0_woSm_n[0]);
  IMnpip_woSm_n_rdata->Draw("HE");
  IMnpip_woK0_woSm_n[1]->SetLineColor(2);
  IMnpip_woK0_woSm_n[1]->Draw("HEsame");
  IMnpip_wK0_woSm_n[1]->SetLineColor(3);
  IMnpip_wK0_woSm_n[1]->Draw("HEsame");
  TH1D* IMnpip_woSm_n_mc = (TH1D*)IMnpip_woK0_woSm_n[1]->Clone("IMnpip_woSm_n_mc");
  IMnpip_woSm_n_mc->Add(IMnpip_wK0_woSm_n[1]);
  IMnpip_woSm_n_mc->SetLineColor(6);
  IMnpip_woSm_n_mc->Draw("HEsame");


  //w/o K0, w/o ((Sp or Sm) & missing neutron)
  TH2D *MMnmiss_IMnpim_woK0_woSid_won_mc = (TH2D*)MMnmiss_IMnpim_woK0_woSid_won[1]->Clone("MMnmiss_IMnpim_woK0_woSid_won_mc");

  TCanvas *cMMnmiss_IMnpim_woK0_woSid_won = new TCanvas("cMMnmiss_IMnpim_woK0_woSid_won","cMMnmiss_IMnpim_woK0_woSid_won",1200,800);
  cMMnmiss_IMnpim_woK0_woSid_won->Divide(2,1);
  cMMnmiss_IMnpim_woK0_woSid_won->cd(1);
  MMnmiss_IMnpim_woK0_woSid_won_rdata->SetTitle("#splitline{MMnmiss_IMnpim_woK0_woSid_won}{  Real data}");
  MMnmiss_IMnpim_woK0_woSid_won_rdata->SetMinimum(1);
  MMnmiss_IMnpim_woK0_woSid_won_rdata->Draw(opt);
  cMMnmiss_IMnpim_woK0_woSid_won->cd(2);
  MMnmiss_IMnpim_woK0_woSid_won_mc->SetTitle("#splitline{MMnmiss_IMnpim_woK0_woSid_won}{  MC sum}");
  MMnmiss_IMnpim_woK0_woSid_won_mc->SetMinimum(1);
  MMnmiss_IMnpim_woK0_woSid_won_mc->SetMaximum(MMnmiss_IMnpim_woK0_woSid_won_rdata->GetMaximum());
  MMnmiss_IMnpim_woK0_woSid_won_mc->Draw(opt);

  TCanvas *cIMnpim_woK0_woSid_won = new TCanvas("IMnpim_woK0_woSid_won","IMnpim_woK0_woSid_won");
  cIMnpim_woK0_woSid_won->cd();
  TH1D* IMnpim_woK0_woSid_won[2];
  for(int i=0; i<2; i++)IMnpim_woK0_woSid_won[i] = (TH1D*)MMnmiss_IMnpim_woK0_woSid_won[i]->ProjectionX(Form("IMnpim_woK0_woSid_won_%s",name[i]));
  TH1D* IMnpim_woK0_woSid_won_mc = MMnmiss_IMnpim_woK0_woSid_won_mc->ProjectionX("IMnpim_woK0_woSid_won_mc");
  IMnpim_woK0_woSid_won[0]->Draw("HE");
  IMnpim_woK0_woSid_won_mc->SetLineColor(6);
  IMnpim_woK0_woSid_won_mc->Draw("HEsame");
  for(int i=1; i<2; i++) {
    IMnpim_woK0_woSid_won[i]->SetLineColor(colordef[i]);
    //IMnpim_woK0_woSid_won[i]->Draw("HEsame");
  }

  //TCanvas *cIMnpim_woK0_woSid_won_ratio = new TCanvas("cIMnpim_woK0_woSid_won_ratio","cIMnpim_woK0_woSid_won_ratio");
  //cIMnpim_woK0_woSid_won_ratio->cd();
  TH1D* IMnpim_woK0_woSid_won_ratio = (TH1D*)IMnpim_woK0_woSid_won[0]->Clone("IMnpim_woK0_woSid_won_ratio");
  IMnpim_woK0_woSid_won_ratio->Divide(IMnpim_woK0_woSid_won_mc);
  IMnpim_woK0_woSid_won_ratio->SetTitle("IMnpim_woK0_woSid_won Data/MC");
  IMnpim_woK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  //IMnpim_woK0_woSid_won_ratio->Draw("HE");
  
  
  TCanvas *cMomnpim_woK0_woSid_won = new TCanvas("Momnpim_woK0_woSid_won","Momnpim_woK0_woSid_won");
  cMomnpim_woK0_woSid_won->cd();
  TH1D* Momnpim_woK0_woSid_won[2];
  for(int i=0; i<2; i++)Momnpim_woK0_woSid_won[i] = (TH1D*)MMnmiss_Momnpim_woK0_woSid_won[i]->ProjectionX(Form("Momnpim_woK0_woSid_won_%s",name[i]));
  TH1D* Momnpim_woK0_woSid_won_mc = MMnmiss_Momnpim_woK0_woSid_won[1]->ProjectionX("Momnpim_woK0_woSid_won_mc");
  Momnpim_woK0_woSid_won[0]->Draw("HE");
  Momnpim_woK0_woSid_won_mc->SetLineColor(6);
  Momnpim_woK0_woSid_won_mc->Draw("HEsame");
  for(int i=1; i<2; i++) {
    Momnpim_woK0_woSid_won[i]->SetLineColor(colordef[i]);
  }

  //TCanvas *cMomnpim_woK0_woSid_won_ratio = new TCanvas("cMomnpim_woK0_woSid_won_ratio","cMomnpim_woK0_woSid_won_ratio");
  //cMomnpim_woK0_woSid_won_ratio->cd();
  //TH1D* Momnpim_woK0_woSid_won_ratio = (TH1D*)Momnpim_woK0_woSid_won[0]->Clone("Momnpim_woK0_woSid_won_ratio");
  //Momnpim_woK0_woSid_won_ratio->Divide(Momnpim_woK0_woSid_won_mc);
  //Momnpim_woK0_woSid_won_ratio->SetTitle("Momnpim_woK0_woSid_won Data/MC");
  //Momnpim_woK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  //Momnpim_woK0_woSid_won_ratio->Draw("HE");
  

  //MMnmiss vs IMpippim wSid
  TH2D* MMnmiss_IMpippim_wSid_mc = (TH2D*)MMnmiss_IMpippim_wSid[1]->Clone("MMnmiss_IMpippim_wSid_mc");
  MMnmiss_IMpippim_wSid_mc->Add(MMnmiss_IMpippim_wSid[2]);

  TCanvas *cMMnmiss_IMpippim_wSid = new TCanvas("cMMnmiss_IMpippim_wSid","cMMnmiss_IMpippim_wSid",1200,800);
  cMMnmiss_IMpippim_wSid->Divide(2,1);
  cMMnmiss_IMpippim_wSid->cd(1);
  MMnmiss_IMpippim_wSid[0]->SetTitle("#splitline{MMnmiss_IMpippim_wSid}{  Real data}");
  MMnmiss_IMpippim_wSid[0]->SetMinimum(1);
  MMnmiss_IMpippim_wSid[0]->Draw(opt);
  cMMnmiss_IMpippim_wSid->cd(2);
  MMnmiss_IMpippim_wSid_mc->SetTitle("#splitline{MMnmiss_IMpippim_wSid}{  MC sum}");
  MMnmiss_IMpippim_wSid_mc->SetMinimum(1);
  MMnmiss_IMpippim_wSid_mc->SetMaximum(MMnmiss_IMpippim_wSid[0]->GetMaximum());
  MMnmiss_IMpippim_wSid_mc->Draw(opt);

  TH1D* IMpippim_wSid_mc = (TH1D*)MMnmiss_IMpippim_wSid_mc->ProjectionX("IMpippim_wSid_mc");
  TH1D* MMnmiss_wSid_mc = (TH1D*)MMnmiss_IMpippim_wSid_mc->ProjectionY("MMnmiss_wSid_mc");
  IMpippim_wSid_mc->SetLineColor(6);
  MMnmiss_wSid_mc->SetLineColor(6);
  TCanvas *cIMpippim_wSid = new TCanvas("cIMpippim_wSid","cIMpippim_wSid");
  cIMpippim_wSid->cd();
  TH1D* IMpippim_wSid[3];
  for(int i=0; i<3; i++) IMpippim_wSid[i] = (TH1D*)MMnmiss_IMpippim_wSid[i]->ProjectionX(Form("IMpippim_wSid_%s",name[i]));
  IMpippim_wSid[0]->Draw("HE");
  IMpippim_wSid_mc->Draw("HEsame");
  for(int i=1; i<3; i++) {
    IMpippim_wSid[i]->SetLineColor(colordef[i]);
    IMpippim_wSid[i]->Draw("HEsame");
  }
  
  TCanvas *cMMnmiss_wSid = new TCanvas("cMMnmiss_wSid","cMMnmiss_wSid");
  cMMnmiss_wSid->cd();
  TH1D* MMnmiss_wSid[3];
  for(int i=0; i<3; i++) MMnmiss_wSid[i] = (TH1D*)MMnmiss_IMpippim_wSid[i]->ProjectionY(Form("MMnmiss_wSid_%s",name[i]));
  MMnmiss_wSid[0]->Draw("HE");
  MMnmiss_wSid_mc->Draw("HEsame");
  for(int i=1; i<3; i++){
    MMnmiss_wSid[i]->SetLineColor(colordef[i]);
    MMnmiss_wSid[i]->Draw("HEsame");
  }

  
  //MMnmiss vs IMpippim wSid w miss n.
  TH2D* MMnmiss_IMpippim_wSid_n_mc = (TH2D*)MMnmiss_IMpippim_wSid_n[1]->Clone("MMnmiss_IMpippim_wSid_n_mc");
  MMnmiss_IMpippim_wSid_n_mc->Add(MMnmiss_IMpippim_wSid_n[2]);

  TCanvas *cMMnmiss_IMpippim_wSid_n = new TCanvas("cMMnmiss_IMpippim_wSid_n","cMMnmiss_IMpippim_wSid_n",1200,800);
  cMMnmiss_IMpippim_wSid_n->Divide(2,1);
  cMMnmiss_IMpippim_wSid_n->cd(1);
  MMnmiss_IMpippim_wSid_n[0]->SetTitle("#splitline{MMnmiss_IMpippim_wSid_n}{  Real data}");
  MMnmiss_IMpippim_wSid_n[0]->SetMinimum(1);
  MMnmiss_IMpippim_wSid_n[0]->Draw(opt);
  cMMnmiss_IMpippim_wSid_n->cd(2);
  MMnmiss_IMpippim_wSid_n_mc->SetTitle("#splitline{MMnmiss_IMpippim_wSid_n}{  MC sum}");
  MMnmiss_IMpippim_wSid_n_mc->SetMinimum(1);
  MMnmiss_IMpippim_wSid_n_mc->SetMaximum(MMnmiss_IMpippim_wSid_n[0]->GetMaximum());
  MMnmiss_IMpippim_wSid_n_mc->Draw(opt);

  TH1D* IMpippim_wSid_n_mc = (TH1D*)MMnmiss_IMpippim_wSid_n_mc->ProjectionX("IMpippim_wSid_n_mc");
  IMpippim_wSid_n_mc->SetLineColor(6);
  TCanvas *cIMpippim_wSid_n = new TCanvas("cIMpippim_wSid_n","cIMpippim_wSid_n");
  cIMpippim_wSid_n->cd();
  TH1D* IMpippim_wSid_n[3];
  for(int i=0; i<3; i++) IMpippim_wSid_n[i] = (TH1D*)MMnmiss_IMpippim_wSid_n[i]->ProjectionX(Form("IMpippim_wSid_n_%s",name[i]));
  IMpippim_wSid_n[0]->Draw("HE");
  IMpippim_wSid_n_mc->Draw("HEsame");
  for(int i=1; i<3; i++) {
    IMpippim_wSid_n[i]->SetLineColor(colordef[i]);
    IMpippim_wSid_n[i]->Draw("HEsame");
  }


  //MMnmiss vs IMpippim w/o K0 w/o (Sid & n);
  TH2D* MMnmiss_IMpippim_woK0_woSid_won_mc = (TH2D*)MMnmiss_IMpippim_woK0_woSid_won[1]->Clone("MMnmiss_IMpippim_woK0_woSid_won_mc");

  TH1D* IMpippim_woK0_woSid_won_mc = (TH1D*)MMnmiss_IMpippim_woK0_woSid_won_mc->ProjectionX("IMpippim_woK0_woSid_won_mc");
  IMpippim_woK0_woSid_won_mc->SetLineColor(6);
  //TCanvas *cIMpippim_woK0_woSid_won = new TCanvas("cIMpippim_woK0_woSid_won","cIMpippim_woK0_woSid_won");
  //cIMpippim_woK0_woSid_won->cd();
  TH1D* IMpippim_woK0_woSid_won[2];
  //for(int i=0; i<2; i++) IMpippim_woK0_woSid_won[i] = (TH1D*)MMnmiss_IMpippim_woK0_woSid_won[i]->ProjectionX(Form("IMpippim_woK0_woSid_won_%s",name[i]));
  //IMpippim_woK0_woSid_won[0]->Draw("HE");
  //IMpippim_woK0_woSid_won_mc->Draw("HEsame");
  //for(int i=1; i<2; i++) {
  //  IMpippim_woK0_woSid_won[i]->SetLineColor(colordef[i]);
    //IMpippim_woK0_woSid_won[i]->Draw("HEsame");
  //}

  //TCanvas *cIMpippim_woK0_woSid_won_ratio = new TCanvas("cIMpippim_woK0_woSid_won_ratio","cIMpippim_woK0_woSid_won_ratio");
  //cIMpippim_woK0_woSid_won_ratio->cd();
  //TH1D* IMpippim_woK0_woSid_won_ratio = (TH1D*)IMpippim_woK0_woSid_won[0]->Clone("IMpippim_woK0_woSid_won_ratio");
  //IMpippim_woK0_woSid_won_ratio->Divide(IMpippim_woK0_woSid_won_mc);
  //IMpippim_woK0_woSid_won_ratio->SetTitle("Data/MC");
  //IMpippim_woK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  //IMpippim_woK0_woSid_won_ratio->Draw("HE");
  
  //MMnmiss vs IMpippim w/o K0 w/o (Sid & n);
  TH2D* MMnmiss_IMpippim_wK0_woSid_won_mc = (TH2D*)MMnmiss_IMpippim_wK0_woSid_won[1]->Clone("MMnmiss_IMpippim_wK0_woSid_won_mc");


  TH1D* IMpippim_wK0_woSid_won_mc = (TH1D*)MMnmiss_IMpippim_wK0_woSid_won_mc->ProjectionX("IMpippim_wK0_woSid_won_mc");
  IMpippim_wK0_woSid_won_mc->SetLineColor(6);
  //TCanvas *cIMpippim_wK0_woSid_won = new TCanvas("cIMpippim_wK0_woSid_won","cIMpippim_wK0_woSid_won");
  //cIMpippim_wK0_woSid_won->cd();
  TH1D* IMpippim_wK0_woSid_won[2];
  for(int i=0; i<2; i++) IMpippim_wK0_woSid_won[i] = (TH1D*)MMnmiss_IMpippim_wK0_woSid_won[i]->ProjectionX(Form("IMpippim_wK0_woSid_won_%s",name[i]));
  //IMpippim_wK0_woSid_won[0]->Draw("HE");
  //IMpippim_wK0_woSid_won_mc->Draw("HEsame");

  //TCanvas *cIMpippim_wK0_woSid_won_ratio = new TCanvas("cIMpippim_wK0_woSid_won_ratio","cIMpippim_wK0_woSid_won_ratio");
  //cIMpippim_wK0_woSid_won_ratio->cd();
  TH1D* IMpippim_wK0_woSid_won_ratio = (TH1D*)IMpippim_wK0_woSid_won[0]->Clone("IMpippim_wK0_woSid_won_ratio");
  IMpippim_wK0_woSid_won_ratio->Divide(IMpippim_wK0_woSid_won_mc);
  IMpippim_wK0_woSid_won_ratio->SetTitle("Data/MC");
  IMpippim_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  //IMpippim_wK0_woSid_won_ratio->Draw("HE");
  
  TCanvas *cMMnmiss_IMpippim_woK0_woSid_won = new TCanvas("cMMnmiss_IMpippim_woK0_woSid_won","cMMnmiss_IMpippim_woK0_woSid_won",1200,800);
  cMMnmiss_IMpippim_woK0_woSid_won->Divide(2,1);
  cMMnmiss_IMpippim_woK0_woSid_won->cd(1);
  MMnmiss_IMpippim_woK0_woSid_won[0]->SetTitle("#splitline{MMnmiss_IMpippim_woK0_woSid_won}{  Real data}");
  MMnmiss_IMpippim_woK0_woSid_won[0]->RebinX(2);
  MMnmiss_IMpippim_woK0_woSid_won[0]->RebinY(4);
  MMnmiss_IMpippim_woK0_woSid_won[0]->Draw(opt);
  cMMnmiss_IMpippim_woK0_woSid_won->cd(2);
  MMnmiss_IMpippim_woK0_woSid_won_mc->SetTitle("#splitline{MMnmiss_IMpippim_woK0_woSid_won}{  MC sum}");
  MMnmiss_IMpippim_woK0_woSid_won_mc->RebinX(2);
  MMnmiss_IMpippim_woK0_woSid_won_mc->RebinY(4);
  MMnmiss_IMpippim_woK0_woSid_won_mc->SetMinimum(1);
  MMnmiss_IMpippim_woK0_woSid_won_mc->SetMaximum(MMnmiss_IMpippim_woK0_woSid_won[0]->GetMaximum());
  MMnmiss_IMpippim_woK0_woSid_won_mc->Draw(opt);


  
  TCanvas *cMMnmiss_IMpippim_wK0_woSid_won = new TCanvas("cMMnmiss_IMpippim_wK0_woSid_won","cMMnmiss_IMpippim_wK0_woSid_won",1200,800);
  cMMnmiss_IMpippim_wK0_woSid_won->Divide(2,1);
  cMMnmiss_IMpippim_wK0_woSid_won->cd(1);
  MMnmiss_IMpippim_wK0_woSid_won[0]->SetTitle("#splitline{MMnmiss_IMpippim_wK0_woSid_won}{  Real data}");
  MMnmiss_IMpippim_wK0_woSid_won[0]->RebinX(2);
  MMnmiss_IMpippim_wK0_woSid_won[0]->RebinY(4);
  MMnmiss_IMpippim_wK0_woSid_won[0]->Draw(opt);
  cMMnmiss_IMpippim_wK0_woSid_won->cd(2);
  MMnmiss_IMpippim_wK0_woSid_won_mc->SetTitle("#splitline{MMnmiss_IMpippim_wK0_woSid_won}{  MC sum}");
  MMnmiss_IMpippim_wK0_woSid_won_mc->RebinX(2);
  MMnmiss_IMpippim_wK0_woSid_won_mc->RebinY(4);
  MMnmiss_IMpippim_wK0_woSid_won_mc->SetMinimum(1);
  MMnmiss_IMpippim_wK0_woSid_won_mc->SetMaximum(MMnmiss_IMpippim_wK0_woSid_won[0]->GetMaximum());
  MMnmiss_IMpippim_wK0_woSid_won_mc->Draw(opt);
  
  
  TCanvas *cMMnmiss_IMpippim_woSid_won = new TCanvas("cMMnmiss_IMpippim_woSid_won","cMMnmiss_IMpippim_woSid_won");
  cMMnmiss_IMpippim_woSid_won->Divide(2,1);
  cMMnmiss_IMpippim_woSid_won->cd(1);
  TH2D* MMnmiss_IMpippim_woSid_won_rdata = (TH2D*)MMnmiss_IMpippim_woK0_woSid_won[0]->Clone("MMnmiss_IMpippim_woSid_won_rdata");
  MMnmiss_IMpippim_woSid_won_rdata->Add(MMnmiss_IMpippim_wK0_woSid_won[0]);
  MMnmiss_IMpippim_woSid_won_rdata->Draw(opt);
  cMMnmiss_IMpippim_woSid_won->cd(2);
  TH2D* MMnmiss_IMpippim_woSid_won_mc = (TH2D*)MMnmiss_IMpippim_woK0_woSid_won[1]->Clone("MMnmiss_IMpippim_woSid_won_mc");
  MMnmiss_IMpippim_woSid_won_mc->Add(MMnmiss_IMpippim_wK0_woSid_won[1]);
  MMnmiss_IMpippim_woSid_won_mc->RebinX(2);
  MMnmiss_IMpippim_woSid_won_mc->RebinY(4);
  MMnmiss_IMpippim_woSid_won_mc->SetMinimum(1);
  MMnmiss_IMpippim_woSid_won_mc->SetMaximum(MMnmiss_IMpippim_woSid_won_rdata->GetMaximum());
  MMnmiss_IMpippim_woSid_won_mc->Draw(opt);


  //q vs IMnpipi w/o K0 w/o (Sid & n);
  TH2D* q_IMnpipi_woK0_woSid_won_mc = (TH2D*)q_IMnpipi_woK0_woSid_won[1]->Clone("q_IMnpipi_woK0_woSid_won_mc");

  TCanvas *cq_IMnpipi_woK0_woSid_won = new TCanvas("cq_IMnpipi_woK0_woSid_won","q_IMnpipi_woK0_woSid_won",1200,800);
  cq_IMnpipi_woK0_woSid_won->Divide(2,1);
  cq_IMnpipi_woK0_woSid_won->cd(1);
  q_IMnpipi_woK0_woSid_won[0]->SetTitle("#splitline{q_IMnpipi_woK0_woSid_won}{  Real data}");
  q_IMnpipi_woK0_woSid_won[0]->SetMinimum(1);
  q_IMnpipi_woK0_woSid_won[0]->Draw(opt);
  cq_IMnpipi_woK0_woSid_won->cd(2);
  q_IMnpipi_woK0_woSid_won_mc->SetTitle("#splitline{q_IMnpipi_woK0_woSid_won}{  MC sum}");
  q_IMnpipi_woK0_woSid_won_mc->SetMinimum(1);
  q_IMnpipi_woK0_woSid_won_mc->SetMaximum(q_IMnpipi_woK0_woSid_won[0]->GetMaximum());
  q_IMnpipi_woK0_woSid_won_mc->Draw(opt);

  TH1D* IMnpipi_woK0_woSid_won_mc = (TH1D*)q_IMnpipi_woK0_woSid_won_mc->ProjectionX("IMnpipi_woK0_woSid_won_mc");
  IMnpipi_woK0_woSid_won_mc->SetLineColor(6);
  TCanvas *cIMnpipi_woK0_woSid_won = new TCanvas("cIMnpipi_woK0_woSid_won","cIMnpipi_woK0_woSid_won");
  cIMnpipi_woK0_woSid_won->cd();
  TH1D* IMnpipi_woK0_woSid_won[2];
  for(int i=0; i<2; i++) IMnpipi_woK0_woSid_won[i] = (TH1D*)q_IMnpipi_woK0_woSid_won[i]->ProjectionX(Form("IMnpipi_woK0_woSid_won_%s",name[i]));
  IMnpipi_woK0_woSid_won[0]->Draw("HE");
  IMnpipi_woK0_woSid_won_mc->Draw("HEsame");
  for(int i=1; i<2; i++) {
    IMnpipi_woK0_woSid_won[i]->SetLineColor(colordef[i]);
    //IMnpipi_woK0_woSid_won[i]->Draw("HEsame");
  }

  //TCanvas *cIMnpipi_woK0_woSid_won_ratio = new TCanvas("cIMnpipi_woK0_woSid_won_ratio","cIMnpipi_woK0_woSid_won_ratio");
  //cIMnpipi_woK0_woSid_won_ratio->cd();
  TH1D* IMnpipi_woK0_woSid_won_ratio = (TH1D*)IMnpipi_woK0_woSid_won[0]->Clone("IMnpipi_woK0_woSid_won_ratio");
  IMnpipi_woK0_woSid_won_ratio->Divide(IMnpipi_woK0_woSid_won_mc);
  IMnpipi_woK0_woSid_won_ratio->SetTitle("Data/MC");
  IMnpipi_woK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  //IMnpipi_woK0_woSid_won_ratio->Draw("HE");
  
  TH1D* q_woK0_woSid_won_mc = (TH1D*)q_IMnpipi_woK0_woSid_won_mc->ProjectionY("q_woK0_woSid_won_mc");
  q_woK0_woSid_won_mc->SetLineColor(6);
  TCanvas *cq_woK0_woSid_won = new TCanvas("cq_woK0_woSid_won","cq_woK0_woSid_won");
  cq_woK0_woSid_won->cd();
  TH1D* q_woK0_woSid_won[2];
  for(int i=0; i<2; i++) q_woK0_woSid_won[i] = (TH1D*)q_IMnpipi_woK0_woSid_won[i]->ProjectionY(Form("q_woK0_woSid_won_%s",name[i]));
  q_woK0_woSid_won[0]->Draw("HE");
  q_woK0_woSid_won_mc->SetLineColor(6);
  q_woK0_woSid_won_mc->Draw("HEsame");
   
  //TCanvas *cq_woK0_woSid_won_ratio = new TCanvas("cq_woK0_woSid_won_ratio","cq_woK0_woSid_won_ratio");
  //cq_woK0_woSid_won_ratio->cd();
  TH1D* q_woK0_woSid_won_ratio = (TH1D*)q_woK0_woSid_won[0]->Clone("q_woK0_woSid_won_ratio");
  q_woK0_woSid_won_ratio->Divide(q_woK0_woSid_won_mc);
  q_woK0_woSid_won_ratio->SetTitle("Data/MC");
  q_woK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  //q_woK0_woSid_won_ratio->Draw("HE");

  
  //missing mass vs Mom(pi+pi-) w/o K0 w/o (Sid & n)
  TH2D* MMnmiss_Mompippim_woK0_woSid_won_mc = (TH2D*)MMnmiss_Mompippim_woK0_woSid_won[1]->Clone("MMnmiss_Mompippim_woK0_woSid_won_mc");

  TCanvas *cMMnmiss_Mompippim_woK0_woSid_won = new TCanvas("cMMnmiss_Mompippim_woK0_woSid_won","cMMnmiss_Mompippim_woK0_woSid_won",1200,800);
  cMMnmiss_Mompippim_woK0_woSid_won->Divide(2,1);
  cMMnmiss_Mompippim_woK0_woSid_won->cd(1);
  MMnmiss_Mompippim_woK0_woSid_won[0]->SetTitle("#splitline{MMnmiss_Mompippim_woK0_woSid_won}{  Real data}");
  MMnmiss_Mompippim_woK0_woSid_won[0]->SetMinimum(1);
  MMnmiss_Mompippim_woK0_woSid_won[0]->Draw(opt);
  cMMnmiss_Mompippim_woK0_woSid_won->cd(2);
  MMnmiss_Mompippim_woK0_woSid_won_mc->SetTitle("#splitline{MMnmiss_Mompippim_woK0_woSid_won}{  MC sum}");
  MMnmiss_Mompippim_woK0_woSid_won_mc->SetMinimum(1);
  MMnmiss_Mompippim_woK0_woSid_won_mc->SetMaximum(MMnmiss_Mompippim_woK0_woSid_won[0]->GetMaximum());  
  MMnmiss_Mompippim_woK0_woSid_won_mc->Draw(opt);

  TH1D* Mompippim_woK0_woSid_won_mc = (TH1D*)MMnmiss_Mompippim_woK0_woSid_won_mc->ProjectionX("Mompippim_woK0_woSid_won_mc");
  Mompippim_woK0_woSid_won_mc->SetLineColor(6);
  TCanvas *cMompippim_woK0_woSid_won = new TCanvas("cMompippim_woK0_woSid_won","cMompippim_woK0_woSid_won");
  cMompippim_woK0_woSid_won->cd();
  TH1D* Mompippim_woK0_woSid_won[2];
  for(int i=0; i<2; i++) Mompippim_woK0_woSid_won[i] = (TH1D*)MMnmiss_Mompippim_woK0_woSid_won[i]->ProjectionX(Form("MMnmiss_Mompippim_woK0_woSid_won_%s",name[i]));
  Mompippim_woK0_woSid_won[0]->SetName("Mompippim_woK0_woSid_won_rdata");
  Mompippim_woK0_woSid_won[0]->Draw("HE");
  Mompippim_woK0_woSid_won_mc->SetLineColor(6);
  Mompippim_woK0_woSid_won_mc->Draw("HEsame");

  TH2D* IMnpim_IMnpip_woK0_woSid_won_mc = (TH2D*)IMnpim_IMnpip_woK0_woSid_won[1]->Clone("IMnpim_IMnpip_woK0_woSid_won_mc");

  TCanvas* cIMnpim_IMnpip_woK0_woSid_won = new TCanvas("cIMnpim_IMnpip_woK0_woSid_won","cIMnpim_IMnpip_woK0_woSid_won");
  cIMnpim_IMnpip_woK0_woSid_won->Divide(2,1);
  cIMnpim_IMnpip_woK0_woSid_won->cd(1);
  IMnpim_IMnpip_woK0_woSid_won[0]->Draw(opt);
  cIMnpim_IMnpip_woK0_woSid_won->cd(2);
  IMnpim_IMnpip_woK0_woSid_won_mc->SetMinimum(1);
  IMnpim_IMnpip_woK0_woSid_won_mc->SetMaximum(IMnpim_IMnpip_woK0_woSid_won[0]->GetMaximum());
  IMnpim_IMnpip_woK0_woSid_won_mc->Draw(opt);

  //nCDSmom
  TH2D* nmom_MMnmiss_woK0_woSid_won_mc = (TH2D*)nmom_MMnmiss_woK0_woSid_won[1]->Clone("nmom_MMnmiss_woK0_woSid_won_mc");

  TH1D* nmom_woK0_woSid_won_mc = (TH1D*)nmom_MMnmiss_woK0_woSid_won_mc->ProjectionY("nmom_woK0_woSid_won_mc");
  TH1D* nmom_woK0_woSid_won[2];
  for(int i=0; i<2; i++) {
    nmom_woK0_woSid_won[i] = (TH1D*)nmom_MMnmiss_woK0_woSid_won[i]->ProjectionY(Form("nmom_woK0_woSid_won_%s",name[i]));
    nmom_woK0_woSid_won[i]->SetLineColor(colordef[i]);
  }
  TCanvas *cnmom_woK0_woSid_won = new TCanvas("cnmom_woK0_woSid_won","cnmom_woK0_woSid_won");
  cnmom_woK0_woSid_won->cd();
  nmom_woK0_woSid_won[0]->Draw("HE");
  nmom_woK0_woSid_won_mc->SetLineColor(6);
  nmom_woK0_woSid_won_mc->Draw("HEsame");

  //TCanvas *cnmom_woK0_woSid_won_ratio = new TCanvas("cnmom_woK0_woSid_won_ratio","cnmom_woK0_woSid_won_ratio");
  //cnmom_woK0_woSid_won_ratio->cd();
  TH1D* nmom_woK0_woSid_won_ratio = (TH1D*)nmom_woK0_woSid_won[0]->Clone("nmom_woK0_woSid_won_ratio");
  nmom_woK0_woSid_won_ratio->Divide(nmom_woK0_woSid_won_mc);
  nmom_woK0_woSid_won_ratio->SetTitle("Data/MC");
  nmom_woK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  //nmom_woK0_woSid_won_ratio->Draw("HEsame");

  
  //weighting to momentum vector 
  TH2D* nmom_cosn_woK0_woSid_won_mc = (TH2D*)nmom_cosn_woK0_woSid_won[1]->Clone("nmom_cosn_woK0_woSid_won_mc");
  TCanvas *cnmom_cosn_woK0_woSid_won = new TCanvas("cnmom_cosn_woK0_woSid_won","cnmom_cosn_woK0_woSid_won",1200,800);
  cnmom_cosn_woK0_woSid_won->Divide(2,1);
  cnmom_cosn_woK0_woSid_won->cd(1);
  nmom_cosn_woK0_woSid_won[0]->SetTitle("#splitline{nmom_cosn_woK0_woSid_won}{Real data}");
  nmom_cosn_woK0_woSid_won[0]->Draw(opt);
  cnmom_cosn_woK0_woSid_won->cd(2);
  nmom_cosn_woK0_woSid_won_mc->SetTitle("#splitline{nmom_cosn_woK0_woSid_won}{MC}");
  nmom_cosn_woK0_woSid_won_mc->SetMinimum(1);
  nmom_cosn_woK0_woSid_won_mc->SetMaximum(nmom_cosn_woK0_woSid_won[0]->GetMaximum());
  nmom_cosn_woK0_woSid_won_mc->Draw(opt);

  TCanvas *ccosn_woK0_woSid_won = new TCanvas("ccosn_woK0_woSid_won","ccosn_woK0_woSid_won");
  ccosn_woK0_woSid_won->cd();
  TH1D* cosn_woK0_woSid_won[2];
  for(int i=0;i<2;i++) cosn_woK0_woSid_won[i] = (TH1D*)nmom_cosn_woK0_woSid_won[i]->ProjectionX(Form("cosn_woK0_woSid_won_%s",name[i]));
  TH1D* cosn_woK0_woSid_won_mc = (TH1D*)nmom_cosn_woK0_woSid_won_mc->ProjectionX("cosn_woK0_woSid_won_mc");
  cosn_woK0_woSid_won_mc->SetLineColor(6);
  cosn_woK0_woSid_won[0]->Draw("HE");
  cosn_woK0_woSid_won_mc->Draw("HEsame");
   
  ////////////////////////////
  // w/K0 + fake neutron
  ///////////////////////////

  //w K0, w/o ((Sp or Sm) & missing neutron)
  TH2D *MMnmiss_IMnpip_wK0_woSid_won_mc = (TH2D*)MMnmiss_IMnpip_wK0_woSid_won[1]->Clone("MMnmiss_IMnpip_wK0_woSid_won_mc");

  TCanvas *cMMnmiss_IMnpip_wK0_woSid_won = new TCanvas("cMMnmiss_IMnpip_wK0_woSid_won","cMMnmiss_IMnpip_wK0_woSid_won",1200,800);
  cMMnmiss_IMnpip_wK0_woSid_won->Divide(2,1);
  cMMnmiss_IMnpip_wK0_woSid_won->cd(1);
  MMnmiss_IMnpip_wK0_woSid_won_rdata->SetTitle("#splitline{MMnmiss_IMnpip_wK0_woSid_won}{  Real data}");
  MMnmiss_IMnpip_wK0_woSid_won_rdata->Draw(opt);
  cMMnmiss_IMnpip_wK0_woSid_won->cd(2);
  MMnmiss_IMnpip_wK0_woSid_won_mc->SetTitle("#splitline{MMnmiss_IMnpip_wK0_woSid_won}{  MC sum}");
  MMnmiss_IMnpip_wK0_woSid_won_mc->SetMinimum(1);
  MMnmiss_IMnpip_wK0_woSid_won_mc->SetMaximum(MMnmiss_IMnpip_wK0_woSid_won_rdata->GetMaximum());
  MMnmiss_IMnpip_wK0_woSid_won_mc->Draw(opt);

  //projection to Missing mass (miss n & Sigma+/-)
  TCanvas *cMMnmiss_wK0_woSid_won = new TCanvas("cMMnmiss_wK0_woSid_won","cMMnmiss_wK0_woSid_won");
  cMMnmiss_wK0_woSid_won->cd();
  TH1D* MMnmiss_wK0_woSid_won[2];
  for(int i=0; i<2; i++) MMnmiss_wK0_woSid_won[i] = (TH1D*)MMnmiss_IMnpip_wK0_woSid_won[i]->ProjectionY(Form("MMnmiss_wK0_woSid_won_%s",name[i]));
  MMnmiss_wK0_woSid_won[0]->Draw("HE");
  TH1D* MMnmiss_wK0_woSid_won_mc = (TH1D*)MMnmiss_IMnpip_wK0_woSid_won_mc->ProjectionY("MMnmiss_wK0_woSid_won_mc");
  MMnmiss_wK0_woSid_won_mc->SetLineColor(6);
  MMnmiss_wK0_woSid_won_mc->Draw("HEsame");

  //Data/MC before modifying MC data
  //TCanvas *cMMnmiss_wK0_woSid_won_ratio = new TCanvas("cMMnmiss_wK0_woSid_won_ratio","cMMnmiss_wK0_woSid_won_ratio");
  //cMMnmiss_wK0_woSid_won_ratio->cd();
  TH1D* MMnmiss_wK0_woSid_won_ratio = (TH1D*)MMnmiss_wK0_woSid_won[0]->Clone("MMnmiss_wK0_woSid_won_ratio");
  MMnmiss_wK0_woSid_won_ratio->Divide(MMnmiss_wK0_woSid_won_mc);
  MMnmiss_wK0_woSid_won_ratio->SetTitle("Data/MC");
  MMnmiss_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  //MMnmiss_wK0_woSid_won_ratio->Draw("HE");


  TH2D *MMom_MMass_wK0_woSid_won_mc = (TH2D*)MMom_MMass_wK0_woSid_won[1]->Clone("MMom_MMass_wK0_woSid_won_mc");
  TCanvas *cMMom_MMass_wK0_woSid_won = new TCanvas("cMMom_MMass_wK0_woSid_won","cMMom_MMass_wK0_woSid_won");
  cMMom_MMass_wK0_woSid_won->Divide(2,1);
  cMMom_MMass_wK0_woSid_won->cd(1);
  MMom_MMass_wK0_woSid_won[0]->Draw(opt);
  cMMom_MMass_wK0_woSid_won->cd(2);
  MMom_MMass_wK0_woSid_won_mc->SetMinimum(1);
  MMom_MMass_wK0_woSid_won_mc->SetMaximum(MMom_MMass_wK0_woSid_won[0]->GetMaximum());
  MMom_MMass_wK0_woSid_won_mc->Draw(opt);

  TH1D* MMom_wK0_woSid_won_mc = (TH1D*)MMom_MMass_wK0_woSid_won_mc->ProjectionY("MMom_wK0_woSid_won_mc");
  TH1D* MMom_wK0_woSid_won[2];
  for(int i=0; i<2; i++) {
    MMom_wK0_woSid_won[i] = (TH1D*)MMom_MMass_wK0_woSid_won[i]->ProjectionY(Form("MMom_wK0_woSid_won_%s",name[i]));
    MMom_wK0_woSid_won[i]->SetLineColor(colordef[i]);
  }
  TCanvas *cMMom_wK0_woSid_won = new TCanvas("cMMom_wK0_woSid_won","cMMom_wK0_woSid_won");
  MMom_wK0_woSid_won[0]->Draw("HE");
  MMom_wK0_woSid_won_mc->SetLineColor(6);
  MMom_wK0_woSid_won_mc->Draw("HEsame");


  //TCanvas *cMMom_wK0_woSid_won_ratio = new TCanvas("cMMom_wK0_woSid_won_ratio","cMMom_wK0_woSid_won_ratio");
  //cMMom_wK0_woSid_won_ratio->cd();
  //TH1D* MMom_wK0_woSid_won_ratio = (TH1D*)MMom_wK0_woSid_won[0]->Clone("MMom_wK0_woSid_won_ratio");
  //MMom_wK0_woSid_won_ratio->Divide(MMom_wK0_woSid_won_mc);
  //MMom_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  //MMom_wK0_woSid_won_ratio->SetTitle("Data/MC");
  //MMom_wK0_woSid_won_ratio->Draw("HE");


  //projection to IMnpip (miss n & Sigma+/-)
  TCanvas *cIMnpip_wK0_woSid_won = new TCanvas("cIMnpip_wK0_woSid_won","cIMnpip_wK0_woSid_won");
  cIMnpip_wK0_woSid_won->cd();
  TH1D* IMnpip_wK0_woSid_won[2];
  for(int i=0; i<2; i++)IMnpip_wK0_woSid_won[i] = (TH1D*)MMnmiss_IMnpip_wK0_woSid_won[i]->ProjectionX(Form("IMnpip_wK0_woSid_won_%s",name[i]));
  IMnpip_wK0_woSid_won[0]->Draw("HE");//rdata
  TH1D* IMnpip_wK0_woSid_won_mc = (TH1D*)MMnmiss_IMnpip_wK0_woSid_won_mc->ProjectionX("IMnpip_wK0_woSid_won_mc");
  IMnpip_wK0_woSid_won_mc->SetLineColor(6);
  IMnpip_wK0_woSid_won_mc->Draw("HEsame");

  //TCanvas *cIMnpip_wK0_woSid_won_ratio = new TCanvas("cIMnpip_wK0_woSid_won_ratio","cIMnpip_wK0_woSid_won_ratio");
  //cIMnpip_wK0_woSid_won_ratio->cd();
  TH1D* IMnpip_wK0_woSid_won_ratio = (TH1D*)IMnpip_wK0_woSid_won[0]->Clone("IMnpip_wK0_woSid_won_ratio");
  IMnpip_wK0_woSid_won_ratio->Divide(IMnpip_wK0_woSid_won_mc);
  IMnpip_wK0_woSid_won_ratio->SetTitle("IMnpip_wK0_woSid_won Data/MC");
  IMnpip_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  //IMnpip_wK0_woSid_won_ratio->Draw("HE");
  
  
  TCanvas *cMomnpip_wK0_woSid_won = new TCanvas("cMomnpip_wK0_woSid_won","cMomnpip_wK0_woSid_won");
  cMomnpip_wK0_woSid_won->cd();
  TH1D* Momnpip_wK0_woSid_won[2];
  for(int i=0; i<2; i++)Momnpip_wK0_woSid_won[i] = (TH1D*)MMnmiss_Momnpip_wK0_woSid_won[i]->ProjectionX(Form("Momnpip_wK0_woSid_won_%s",name[i]));
  Momnpip_wK0_woSid_won[0]->Draw("HE");//rdata
  TH1D* Momnpip_wK0_woSid_won_mc = (TH1D*)MMnmiss_Momnpip_wK0_woSid_won[1]->ProjectionX("Momnpip_wK0_woSid_won_mc");
  Momnpip_wK0_woSid_won_mc->SetLineColor(6);
  Momnpip_wK0_woSid_won_mc->Draw("HEsame");

//  TCanvas *cMomnpip_wK0_woSid_won_ratio = new TCanvas("cMomnpip_wK0_woSid_won_ratio","cMomnpip_wK0_woSid_won_ratio");
//  cMomnpip_wK0_woSid_won_ratio->cd();
//  TH1D* Momnpip_wK0_woSid_won_ratio = (TH1D*)Momnpip_wK0_woSid_won[0]->Clone("Momnpip_wK0_woSid_won_ratio");
//  Momnpip_wK0_woSid_won_ratio->Divide(Momnpip_wK0_woSid_won_mc);
//  Momnpip_wK0_woSid_won_ratio->SetTitle("Momnpip_wK0_woSid_won Data/MC");
//  Momnpip_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
//  Momnpip_wK0_woSid_won_ratio->Draw("HE");

  //Sigma-
  //w K0, w/o ((Sp or Sm) & missing neutron)
  TH2D *MMnmiss_IMnpim_wK0_woSid_won_mc = (TH2D*)MMnmiss_IMnpim_wK0_woSid_won[1]->Clone("MMnmiss_IMnpim_wK0_woSid_won_mc");

  TCanvas *cMMnmiss_IMnpim_wK0_woSid_won = new TCanvas("cMMnmiss_IMnpim_wK0_woSid_won","cMMnmiss_IMnpim_wK0_woSid_won",1200,800);
  cMMnmiss_IMnpim_wK0_woSid_won->Divide(2,1);
  cMMnmiss_IMnpim_wK0_woSid_won->cd(1);
  MMnmiss_IMnpim_wK0_woSid_won_rdata->SetTitle("#splitline{MMnmiss_IMnpim_wK0_woSid_won}{  Real data}");
  MMnmiss_IMnpim_wK0_woSid_won_rdata->Draw(opt);
  cMMnmiss_IMnpim_wK0_woSid_won->cd(2);
  MMnmiss_IMnpim_wK0_woSid_won_mc->SetTitle("#splitline{MMnmiss_IMnpim_wK0_woSid_won}{  MC sum}");
  MMnmiss_IMnpim_wK0_woSid_won_mc->SetMinimum(1);
  MMnmiss_IMnpim_wK0_woSid_won_mc->SetMaximum(MMnmiss_IMnpim_wK0_woSid_won_rdata->GetMaximum());
  MMnmiss_IMnpim_wK0_woSid_won_mc->Draw(opt);

  TCanvas *cIMnpim_wK0_woSid_won = new TCanvas("IMnpim_wK0_woSid_won","IMnpim_wK0_woSid_won");
  cIMnpim_wK0_woSid_won->cd();
  TH1D* IMnpim_wK0_woSid_won[2];
  for(int i=0; i<2; i++)IMnpim_wK0_woSid_won[i] = (TH1D*)MMnmiss_IMnpim_wK0_woSid_won[i]->ProjectionX(Form("IMnpim_wK0_woSid_won_%s",name[i]));
  TH1D* IMnpim_wK0_woSid_won_mc = MMnmiss_IMnpim_wK0_woSid_won_mc->ProjectionX("IMnpim_wK0_woSid_won_mc");
  IMnpim_wK0_woSid_won[0]->Draw("HE");
  IMnpim_wK0_woSid_won_mc->SetLineColor(6);
  IMnpim_wK0_woSid_won_mc->Draw("HEsame");

  //TCanvas *cIMnpim_wK0_woSid_won_ratio = new TCanvas("cIMnpim_wK0_woSid_won_ratio","cIMnpim_wK0_woSid_won_ratio");
  //cIMnpim_wK0_woSid_won_ratio->cd();
  TH1D* IMnpim_wK0_woSid_won_ratio = (TH1D*)IMnpim_wK0_woSid_won[0]->Clone("IMnpim_wK0_woSid_won_ratio");
  IMnpim_wK0_woSid_won_ratio->Divide(IMnpim_wK0_woSid_won_mc);
  IMnpim_wK0_woSid_won_ratio->SetTitle("IMnpim_wK0_woSid_won Data/MC");
  IMnpim_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  //IMnpim_wK0_woSid_won_ratio->Draw("HE");
  
  
  TCanvas *cMomnpim_wK0_woSid_won = new TCanvas("Momnpim_wK0_woSid_won","Momnpim_wK0_woSid_won");
  cMomnpim_wK0_woSid_won->cd();
  TH1D* Momnpim_wK0_woSid_won[2];
  for(int i=0; i<2; i++)Momnpim_wK0_woSid_won[i] = (TH1D*)MMnmiss_Momnpim_wK0_woSid_won[i]->ProjectionX(Form("Momnpim_wK0_woSid_won_%s",name[i]));
  TH1D* Momnpim_wK0_woSid_won_mc = MMnmiss_Momnpim_wK0_woSid_won[1]->ProjectionX("Momnpim_wK0_woSid_won_mc");
  Momnpim_wK0_woSid_won[0]->Draw("HE");
  Momnpim_wK0_woSid_won_mc->SetLineColor(6);
  Momnpim_wK0_woSid_won_mc->Draw("HEsame");
  for(int i=1; i<2; i++) {
    Momnpim_wK0_woSid_won[i]->SetLineColor(colordef[i]);
  }

  //TCanvas *cMomnpim_wK0_woSid_won_ratio = new TCanvas("cMomnpim_wK0_woSid_won_ratio","cMomnpim_wK0_woSid_won_ratio");
  //cMomnpim_wK0_woSid_won_ratio->cd();
  //TH1D* Momnpim_wK0_woSid_won_ratio = (TH1D*)Momnpim_wK0_woSid_won[0]->Clone("Momnpim_wK0_woSid_won_ratio");
  //Momnpim_wK0_woSid_won_ratio->Divide(Momnpim_wK0_woSid_won_mc);
  //Momnpim_wK0_woSid_won_ratio->SetTitle("Momnpim_wK0_woSid_won Data/MC");
  //Momnpim_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  //Momnpim_wK0_woSid_won_ratio->Draw("HE");
 


  //q vs IMnpipi w/o K0 w/o (Sid & n);
  TH2D* q_IMnpipi_wK0_woSid_won_mc = (TH2D*)q_IMnpipi_wK0_woSid_won[1]->Clone("q_IMnpipi_wK0_woSid_won_mc");

  TCanvas *cq_IMnpipi_wK0_woSid_won = new TCanvas("cq_IMnpipi_wK0_woSid_won","q_IMnpipi_wK0_woSid_won",1200,800);
  cq_IMnpipi_wK0_woSid_won->Divide(2,1);
  cq_IMnpipi_wK0_woSid_won->cd(1);
  q_IMnpipi_wK0_woSid_won[0]->SetTitle("#splitline{q_IMnpipi_wK0_woSid_won}{  Real data}");
  q_IMnpipi_wK0_woSid_won[0]->Draw(opt);
  cq_IMnpipi_wK0_woSid_won->cd(2);
  q_IMnpipi_wK0_woSid_won_mc->SetTitle("#splitline{q_IMnpipi_wK0_woSid_won}{  MC sum}");
  q_IMnpipi_wK0_woSid_won_mc->SetMinimum(1);
  q_IMnpipi_wK0_woSid_won_mc->SetMaximum(q_IMnpipi_wK0_woSid_won[0]->GetMaximum());
  q_IMnpipi_wK0_woSid_won_mc->Draw(opt);

  TH1D* IMnpipi_wK0_woSid_won_mc = (TH1D*)q_IMnpipi_wK0_woSid_won_mc->ProjectionX("IMnpipi_wK0_woSid_won_mc");
  IMnpipi_wK0_woSid_won_mc->SetLineColor(6);
  TCanvas *cIMnpipi_wK0_woSid_won = new TCanvas("cIMnpipi_wK0_woSid_won","cIMnpipi_wK0_woSid_won");
  cIMnpipi_wK0_woSid_won->cd();
  TH1D* IMnpipi_wK0_woSid_won[2];
  for(int i=0; i<2; i++) IMnpipi_wK0_woSid_won[i] = (TH1D*)q_IMnpipi_wK0_woSid_won[i]->ProjectionX(Form("IMnpipi_wK0_woSid_won_%s",name[i]));
  IMnpipi_wK0_woSid_won[0]->Draw("HE");
  IMnpipi_wK0_woSid_won_mc->Draw("HEsame");

  //TCanvas *cIMnpipi_wK0_woSid_won_ratio = new TCanvas("cIMnpipi_wK0_woSid_won_ratio","cIMnpipi_wK0_woSid_won_ratio");
  //cIMnpipi_wK0_woSid_won_ratio->cd();
  TH1D* IMnpipi_wK0_woSid_won_ratio = (TH1D*) IMnpipi_wK0_woSid_won[0]->Clone("IMnpipi_wK0_woSid_won_ratio");
  IMnpipi_wK0_woSid_won_ratio->Divide(IMnpipi_wK0_woSid_won_mc);
  IMnpipi_wK0_woSid_won_ratio->SetTitle("Data/MC");
  IMnpipi_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  //IMnpipi_wK0_woSid_won_ratio->Draw("HE");


  TH1D* q_wK0_woSid_won_mc = (TH1D*)q_IMnpipi_wK0_woSid_won_mc->ProjectionY("q_wK0_woSid_won_mc");
  q_wK0_woSid_won_mc->SetLineColor(6);
  TCanvas *cq_wK0_woSid_won = new TCanvas("cq_wK0_woSid_won","cq_wK0_woSid_won");
  cq_wK0_woSid_won->cd();
  TH1D* q_wK0_woSid_won[2];
  for(int i=0; i<2; i++) q_wK0_woSid_won[i] = (TH1D*)q_IMnpipi_wK0_woSid_won[i]->ProjectionY(Form("q_wK0_woSid_won_%s",name[i]));
  q_wK0_woSid_won[0]->Draw("HE");
  q_wK0_woSid_won_mc->SetLineColor(6);
  q_wK0_woSid_won_mc->Draw("HEsame");
 
  //TCanvas *cq_wK0_woSid_won_ratio = new TCanvas("cq_wK0_woSid_won_ratio","cq_wK0_woSid_won_ratio");
  TH1D* q_wK0_woSid_won_ratio = (TH1D*)q_wK0_woSid_won[0]->Clone("q_wK0_woSid_won_ratio");
  q_wK0_woSid_won_ratio->Divide(q_wK0_woSid_won_mc);
  q_wK0_woSid_won_ratio->SetTitle("Data/MC");
  q_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  //q_wK0_woSid_won_ratio->Draw("HE");
  

  //missing mass vs Mom(pi+pi-) w/o K0 w/o (Sid & n)
  TH2D* MMnmiss_Mompippim_wK0_woSid_won_mc = (TH2D*)MMnmiss_Mompippim_wK0_woSid_won[1]->Clone("MMnmiss_Mompippim_wK0_woSid_won_mc");

  TCanvas *cMMnmiss_Mompippim_wK0_woSid_won = new TCanvas("cMMnmiss_Mompippim_wK0_woSid_won","cMMnmiss_Mompippim_wK0_woSid_won",1200,800);
  cMMnmiss_Mompippim_wK0_woSid_won->Divide(2,1);
  cMMnmiss_Mompippim_wK0_woSid_won->cd(1);
  MMnmiss_Mompippim_wK0_woSid_won[0]->SetTitle("#splitline{MMnmiss_Mompippim_wK0_woSid_won}{  Real data}");
  MMnmiss_Mompippim_wK0_woSid_won[0]->Draw(opt);
  cMMnmiss_Mompippim_wK0_woSid_won->cd(2);
  MMnmiss_Mompippim_wK0_woSid_won_mc->SetTitle("#splitline{MMnmiss_Mompippim_wK0_woSid_won}{  MC sum}");
  MMnmiss_Mompippim_wK0_woSid_won_mc->SetMinimum(1);
  MMnmiss_Mompippim_wK0_woSid_won_mc->SetMaximum(MMnmiss_Mompippim_wK0_woSid_won[0]->GetMaximum());
  MMnmiss_Mompippim_wK0_woSid_won_mc->Draw(opt);

  TH1D* Mompippim_wK0_woSid_won_mc = (TH1D*)MMnmiss_Mompippim_wK0_woSid_won_mc->ProjectionX("Mompippim_wK0_woSid_won_mc");
  Mompippim_wK0_woSid_won_mc->SetLineColor(6);
  TCanvas *cMompippim_wK0_woSid_won = new TCanvas("cMompippim_wK0_woSid_won","cMompippim_wK0_woSid_won");
  cMompippim_wK0_woSid_won->cd();
  TH1D* Mompippim_wK0_woSid_won[2];
  for(int i=0; i<2; i++) Mompippim_wK0_woSid_won[i] = (TH1D*)MMnmiss_Mompippim_wK0_woSid_won[i]->ProjectionX(Form("MMnmiss_Mompippim_wK0_woSid_won_%s",name[i]));
  Mompippim_wK0_woSid_won[0]->SetName("Mompippim_wK0_woSid_won_rdata");
  Mompippim_wK0_woSid_won[0]->Draw("HE");
  Mompippim_wK0_woSid_won_mc->SetLineColor(6);

  //TCanvas *cMompippim_wK0_woSid_won_ratio = new TCanvas("cMompippim_wK0_woSid_won_ratio","cMompippim_wK0_woSid_won_ratio");
  //TH1D* Mompippim_wK0_woSid_won_ratio = (TH1D*)Mompippim_wK0_woSid_won[0]->Clone("Mompippim_wK0_woSid_won_ratio");
  //Mompippim_wK0_woSid_won_ratio->Divide(Mompippim_wK0_woSid_won_mc);
  //Mompippim_wK0_woSid_won_ratio->SetTitle("Data/MC");
  //Mompippim_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  //Mompippim_wK0_woSid_won_ratio->Draw("HE");

  //nCDSmom
  TH2D* nmom_MMnmiss_wK0_woSid_won_mc = (TH2D*)nmom_MMnmiss_wK0_woSid_won[1]->Clone("nmom_MMnmiss_wK0_woSid_won_mc");

  TH1D* nmom_wK0_woSid_won_mc = (TH1D*)nmom_MMnmiss_wK0_woSid_won_mc->ProjectionY("nmom_wK0_woSid_won_mc");
  TH1D* nmom_wK0_woSid_won[2];
  for(int i=0; i<2; i++) {
    nmom_wK0_woSid_won[i] = (TH1D*)nmom_MMnmiss_wK0_woSid_won[i]->ProjectionY(Form("nmom_wK0_woSid_won_%s",name[i]));
    nmom_wK0_woSid_won[i]->SetLineColor(colordef[i]);
  }
  TCanvas *cnmom_wK0_woSid_won = new TCanvas("cnmom_wK0_woSid_won","cnmom_wK0_woSid_won");
  cnmom_wK0_woSid_won->cd();
  nmom_wK0_woSid_won[0]->Draw("HE");
  nmom_wK0_woSid_won_mc->SetLineColor(6);
  nmom_wK0_woSid_won_mc->Draw("HEsame");
  //for(int i=1;i<6;i++)nmom_wK0_woSid_won[i]->Draw("same");

  //TCanvas *cnmom_wK0_woSid_won_ratio = new TCanvas("cnmom_wK0_woSid_won_ratio","cnmom_wK0_woSid_won_ratio");
  //cnmom_wK0_woSid_won_ratio->cd();
  TH1D* nmom_wK0_woSid_won_ratio = (TH1D*)nmom_wK0_woSid_won[0]->Clone("nmom_wK0_woSid_won_ratio");
  nmom_wK0_woSid_won_ratio->Divide(nmom_wK0_woSid_won_mc);
  nmom_wK0_woSid_won_ratio->SetTitle("Data/MC");
  nmom_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  //nmom_wK0_woSid_won_ratio->Draw("HEsame");

  
  //weighting to momentum vector 
  TH2D* nmom_cosn_wK0_woSid_won_mc = (TH2D*)nmom_cosn_wK0_woSid_won[1]->Clone("nmom_cosn_wK0_woSid_won_mc");
  TCanvas *cnmom_cosn_wK0_woSid_won = new TCanvas("cnmom_cosn_wK0_woSid_won","cnmom_cosn_wK0_woSid_won",1200,800);
  cnmom_cosn_wK0_woSid_won->Divide(2,1);
  cnmom_cosn_wK0_woSid_won->cd(1);
  nmom_cosn_wK0_woSid_won[0]->SetTitle("#splitline{nmom_cosn_wK0_woSid_won}{Real data}");
  nmom_cosn_wK0_woSid_won[0]->Draw(opt);
  cnmom_cosn_wK0_woSid_won->cd(2);
  nmom_cosn_wK0_woSid_won_mc->SetTitle("#splitline{nmom_cosn_wK0_woSid_won}{MC}");
  nmom_cosn_wK0_woSid_won_mc->SetMinimum(1);
  nmom_cosn_wK0_woSid_won_mc->SetMaximum(nmom_cosn_wK0_woSid_won[0]->GetMaximum());
  nmom_cosn_wK0_woSid_won_mc->Draw(opt);

  TCanvas *ccosn_wK0_woSid_won = new TCanvas("ccosn_wK0_woSid_won","ccosn_wK0_woSid_won");
  ccosn_wK0_woSid_won->cd();
  TH1D* cosn_wK0_woSid_won[2];
  for(int i=0;i<2;i++) cosn_wK0_woSid_won[i] = (TH1D*)nmom_cosn_wK0_woSid_won[i]->ProjectionX(Form("cosn_wK0_woSid_won_%s",name[i]));
  TH1D* cosn_wK0_woSid_won_mc = (TH1D*)nmom_cosn_wK0_woSid_won_mc->ProjectionX("cosn_wK0_woSid_won_mc");
  cosn_wK0_woSid_won_mc->SetLineColor(6);
  cosn_wK0_woSid_won[0]->Draw("HE");
  cosn_wK0_woSid_won_mc->Draw("HEsame");
   
  /*
  TCanvas *ccosn_wK0_woSid_won_ratio = new TCanvas("ccosn_wK0_woSid_won_ratio","ccosn_wK0_woSid_won_ratio");
  ccosn_wK0_woSid_won_ratio->cd();
  TH1D* cosn_wK0_woSid_won_ratio = (TH1D*)cosn_wK0_woSid_won[0]->Clone("cosn_wK0_woSid_won_ratio");
  cosn_wK0_woSid_won_ratio->Divide(cosn_wK0_woSid_won_mc);
  cosn_wK0_woSid_won_ratio->SetTitle("Data/MC");
  cosn_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  cosn_wK0_woSid_won_ratio->Draw("HE");
  */
  
  TCanvas *cMMnmiss_Momnpipi_woK0_woSid_won = new TCanvas("cMMnmiss_Momnpipi_woK0_woSid_won","cMMnmiss_Momnpipi_woK0_woSid_won",1000,800);
  cMMnmiss_Momnpipi_woK0_woSid_won->Divide(2,1);
  cMMnmiss_Momnpipi_woK0_woSid_won->cd(1);
  MMnmiss_Momnpipi_woK0_woSid_won_rdata->Draw(opt);
  cMMnmiss_Momnpipi_woK0_woSid_won->cd(2);
  MMnmiss_Momnpipi_woK0_woSid_won[1]->SetMaximum(MMnmiss_Momnpipi_woK0_woSid_won_rdata->GetMaximum());
  MMnmiss_Momnpipi_woK0_woSid_won[1]->Draw(opt);//fake
   
  /*
  TCanvas *cMomnpipi_woK0_woSid_won_ratio = new TCanvas("cMomnpipi_woK0_woSid_won_ratio","cMomnpipi_woK0_woSid_won_ratio");
  cMomnpipi_woK0_woSid_won_ratio->cd();
  TH1D* Momnpipi_woK0_woSid_won_rdata = (TH1D*)MMnmiss_Momnpipi_woK0_woSid_won_rdata->ProjectionX("Momnpipi_woK0_woSid_won_rdata");
  TH1D* Momnpipi_woK0_woSid_won_mc = (TH1D*)MMnmiss_Momnpipi_woK0_woSid_won[1]->ProjectionX("Momnpipi_woK0_woSid_won_mc");
  Momnpipi_woK0_woSid_won_mc->SetLineColor(6);
  TH1D* Momnpipi_woK0_woSid_won_ratio = (TH1D*)Momnpipi_woK0_woSid_won_rdata->Clone("Momnpipi_woK0_woSid_won_ratio");
  Momnpipi_woK0_woSid_won_ratio->Divide(Momnpipi_woK0_woSid_won_mc);
  Momnpipi_woK0_woSid_won_ratio->SetTitle("Data/MC");
  Momnpipi_woK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  Momnpipi_woK0_woSid_won_ratio->Draw("HE");
  */
  
  TCanvas *cMMnmiss_Momnpipi_wK0_woSid_won = new TCanvas("cMMnmiss_Momnpipi_wK0_woSid_won","cMMnmiss_Momnpipi_wK0_woSid_won",1000,800);
  cMMnmiss_Momnpipi_wK0_woSid_won->Divide(2,1);
  cMMnmiss_Momnpipi_wK0_woSid_won->cd(1);
  MMnmiss_Momnpipi_wK0_woSid_won_rdata->Draw(opt);
  cMMnmiss_Momnpipi_wK0_woSid_won->cd(2);
  MMnmiss_Momnpipi_wK0_woSid_won[1]->SetMaximum(MMnmiss_Momnpipi_wK0_woSid_won_rdata->GetMaximum());
  MMnmiss_Momnpipi_wK0_woSid_won[1]->Draw(opt);//fake
  
  /*
  TCanvas *cMomnpipi_wK0_woSid_won_ratio = new TCanvas("cMomnpipi_wK0_woSid_won_ratio","cMomnpipi_wK0_woSid_won_ratio");
  cMomnpipi_wK0_woSid_won_ratio->cd();
  TH1D* Momnpipi_wK0_woSid_won_rdata = (TH1D*)MMnmiss_Momnpipi_wK0_woSid_won_rdata->ProjectionX("Momnpipi_wK0_woSid_won_rdata");
  TH1D* Momnpipi_wK0_woSid_won_mc = (TH1D*)MMnmiss_Momnpipi_wK0_woSid_won[1]->ProjectionX("Momnpipi_wK0_woSid_won_mc");
  Momnpipi_wK0_woSid_won_mc->SetLineColor(6);
  TH1D* Momnpipi_wK0_woSid_won_ratio = (TH1D*)Momnpipi_wK0_woSid_won_rdata->Clone("Momnpipi_wK0_woSid_won_ratio");
  Momnpipi_wK0_woSid_won_ratio->Divide(Momnpipi_wK0_woSid_won_mc);
  Momnpipi_wK0_woSid_won_ratio->SetTitle("Data/MC");
  Momnpipi_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  Momnpipi_wK0_woSid_won_ratio->Draw("HE");
  */

  //Signal check
  TH2D* q_IMnpipi_wSid_n_Sp_mc = (TH2D*)q_IMnpipi_wSid_n_Sp[1]->Clone();//woK0
  q_IMnpipi_wSid_n_Sp_mc->Add(q_IMnpipi_wSid_n_Sp[2]);

  TCanvas *cq_IMnpipi_wSid_n_Sp = new TCanvas("cq_IMnpipi_wSid_n_Sp","cq_IMnpipi_wSid_n_Sp");
  cq_IMnpipi_wSid_n_Sp->Divide(2,1);
  cq_IMnpipi_wSid_n_Sp->cd(1);
  q_IMnpipi_wSid_n_Sp[0]->SetMinimum(1);
  q_IMnpipi_wSid_n_Sp[0]->Draw(opt);
  cq_IMnpipi_wSid_n_Sp->cd(2);
  q_IMnpipi_wSid_n_Sp_mc->SetMinimum(1);
  q_IMnpipi_wSid_n_Sp_mc->SetMaximum(q_IMnpipi_wSid_n_Sp[0]->GetMaximum());
  q_IMnpipi_wSid_n_Sp_mc->Draw(opt);

  TCanvas *cIMnpipi_wSid_n_Sp_0 = new TCanvas("cIMnpipi_wSid_n_Sp_0","cIMnpipi_wSid_n_Sp_0");
  const int q350bin = q_IMnpipi_wSid_n_Sp[0]->GetYaxis()->FindBin(0.35);
  TH1D* IMnpipi_wSid_n_Sp_0[3];
  for(int i=0; i<3; i++)IMnpipi_wSid_n_Sp_0[i] = q_IMnpipi_wSid_n_Sp[i]->ProjectionX(Form("IMnpipi_wSid_n_Sp_0_%s",name[i]),0,q350bin-1);
  IMnpipi_wSid_n_Sp_0[0]->Draw("HE");
  TH1D* IMnpipi_wSid_n_Sp_mc_0 = q_IMnpipi_wSid_n_Sp_mc->ProjectionX("IMnpipi_wSid_n_Sp_mc_0",0,q350bin-1);
  IMnpipi_wSid_n_Sp_mc_0->SetLineColor(6);
  IMnpipi_wSid_n_Sp_mc_0->Draw("HEsame");
  IMnpipi_wSid_n_Sp_0[1]->SetLineColor(2);
  IMnpipi_wSid_n_Sp_0[1]->Draw("HEsame");
  IMnpipi_wSid_n_Sp_0[2]->SetLineColor(3);
  IMnpipi_wSid_n_Sp_0[2]->Draw("HEsame");


  TCanvas *cIMnpipi_wSid_n_Sp_350 = new TCanvas("cIMnpipi_wSid_n_Sp_350","cIMnpipi_wSid_n_Sp_350");
  TH1D* IMnpipi_wSid_n_Sp_350[3];
  const int q600bin = q_IMnpipi_wSid_n_Sp[0]->GetYaxis()->FindBin(2.0);
  for(int i=0; i<3; i++)IMnpipi_wSid_n_Sp_350[i] = q_IMnpipi_wSid_n_Sp[i]->ProjectionX(Form("IMnpipi_wSid_n_Sp_350_%s",name[i]),q350bin,q600bin);
  IMnpipi_wSid_n_Sp_350[0]->Draw("HE");
  TH1D* IMnpipi_wSid_n_Sp_mc_350 = q_IMnpipi_wSid_n_Sp_mc->ProjectionX("IMnpipi_wSid_n_Sp_mc_350",q350bin,q600bin);
  IMnpipi_wSid_n_Sp_mc_350->SetLineColor(6);
  IMnpipi_wSid_n_Sp_mc_350->Draw("HEsame");
  IMnpipi_wSid_n_Sp_350[1]->SetLineColor(2);
  IMnpipi_wSid_n_Sp_350[1]->Draw("HEsame");
  IMnpipi_wSid_n_Sp_350[2]->SetLineColor(3);
  IMnpipi_wSid_n_Sp_350[2]->Draw("HEsame");


  TH2D* q_IMnpipi_wSid_n_Sm_mc = (TH2D*)q_IMnpipi_wSid_n_Sm[1]->Clone();

  TCanvas *cq_IMnpipi_wSid_n_Sm = new TCanvas("cq_IMnpipi_wSid_n_Sm","cq_IMnpipi_wSid_n_Sm");
  cq_IMnpipi_wSid_n_Sm->Divide(2,1);
  cq_IMnpipi_wSid_n_Sm->cd(1);
  q_IMnpipi_wSid_n_Sm[0]->SetMinimum(1);
  q_IMnpipi_wSid_n_Sm[0]->Draw(opt);
  cq_IMnpipi_wSid_n_Sm->cd(2);
  q_IMnpipi_wSid_n_Sm_mc->SetMinimum(1);
  q_IMnpipi_wSid_n_Sm_mc->SetMaximum(q_IMnpipi_wSid_n_Sm[0]->GetMaximum());
  q_IMnpipi_wSid_n_Sm_mc->Draw(opt);

  TCanvas *cIMnpipi_wSid_n_Sm_0 = new TCanvas("cIMnpipi_wSid_n_Sm_0","cIMnpipi_wSid_n_Sm_0");
  TH1D* IMnpipi_wSid_n_Sm_0[3];
  for(int i=0; i<3; i++)IMnpipi_wSid_n_Sm_0[i] = q_IMnpipi_wSid_n_Sm[i]->ProjectionX(Form("IMnpipi_wSid_n_Sm_0_%s",name[i]),0,q350bin-1);
  IMnpipi_wSid_n_Sm_0[0]->Draw("HE");
  TH1D* IMnpipi_wSid_n_Sm_mc_0 = q_IMnpipi_wSid_n_Sm_mc->ProjectionX("IMnpipi_wSid_n_Sm_mc_0",0,q350bin-1);
  IMnpipi_wSid_n_Sm_mc_0->SetLineColor(6);
  IMnpipi_wSid_n_Sm_mc_0->Draw("HEsame");
  IMnpipi_wSid_n_Sm_0[1]->SetLineColor(2);
  IMnpipi_wSid_n_Sm_0[1]->Draw("HEsame");
  IMnpipi_wSid_n_Sm_0[2]->SetLineColor(3);
  IMnpipi_wSid_n_Sm_0[2]->Draw("HEsame");

  TCanvas *cIMnpipi_wSid_n_Sm_350 = new TCanvas("cIMnpipi_wSid_n_Sm_350","cIMnpipi_wSid_n_Sm_350");
  TH1D* IMnpipi_wSid_n_Sm_350[3];
  for(int i=0; i<3; i++)IMnpipi_wSid_n_Sm_350[i] = q_IMnpipi_wSid_n_Sm[i]->ProjectionX(Form("IMnpipi_wSid_n_Sm_350_%s",name[i]),q350bin,q600bin);
  IMnpipi_wSid_n_Sm_350[0]->Draw("HE");
  TH1D* IMnpipi_wSid_n_Sm_mc_350 = q_IMnpipi_wSid_n_Sm_mc->ProjectionX("IMnpipi_wSid_n_Sm_mc_350",q350bin,q600bin);
  IMnpipi_wSid_n_Sm_mc_350->SetLineColor(6);
  IMnpipi_wSid_n_Sm_mc_350->Draw("HEsame");
  IMnpipi_wSid_n_Sm_350[1]->SetLineColor(2);
  IMnpipi_wSid_n_Sm_350[1]->Draw("HEsame");
  IMnpipi_wSid_n_Sm_350[2]->SetLineColor(3);
  IMnpipi_wSid_n_Sm_350[2]->Draw("HEsame");
   
 
  TCanvas *cq_IMnpipi_wSid_n_rdata = new TCanvas("cq_IMnpipi_wSid_n_rdata","cq_IMnpipi_wSid_n_rdata");
  cq_IMnpipi_wSid_n_rdata->Divide(2,1);
  cq_IMnpipi_wSid_n_rdata->cd(1);
  q_IMnpipi_wSid_n[0]->Draw(opt);
  //TCanvas *cq_IMnpipi_wSid_n_mc = new TCanvas("cq_IMnpipi_wSid_n_mc","cq_IMnpipi_wSid_n_mc");
  cq_IMnpipi_wSid_n_rdata->cd(2);
  TH2D* q_IMnpipi_wSid_n_mc = (TH2D*)q_IMnpipi_wSid_n[1]->Clone("q_IMnpipi_wSid_n_mc");
  q_IMnpipi_wSid_n_mc->SetTitle("q_IMnpipi_wSid_n_mc");
  q_IMnpipi_wSid_n_mc->Add(q_IMnpipi_wSid_n[2]);
  q_IMnpipi_wSid_n_mc->SetMaximum(q_IMnpipi_wSid_n[0]->GetMaximum());
  q_IMnpipi_wSid_n_mc->Draw(opt);
   

  TCanvas *cq_IMnpipi_wSid_mc = new TCanvas("cq_IMnpipi_wSid_mc","cq_IMnpipi_wSid_mc");
  cq_IMnpipi_wSid_mc->Divide(2,1);
  cq_IMnpipi_wSid_mc->cd(1);
  q_IMnpipi_wSid_n[1]->Draw(opt);
  cq_IMnpipi_wSid_mc->cd(2);
  q_IMnpipi_wSid_n[2]->SetMaximum(q_IMnpipi_wSid_n[1]->GetMaximum());
  q_IMnpipi_wSid_n[2]->Draw(opt);
 

  TCanvas *cq_IMnpipi_wSid_n_sub = new TCanvas("cq_IMnpipi_wSid_n_sub","cq_IMnpipi_wSid_n_sub");
  TH1D* q_IMnpipi_wSid_n_sub = (TH1D*)q_IMnpipi_wSid_n[0]->Clone("q_IMnpipi_wSid_n_sub");
  q_IMnpipi_wSid_n_sub->Add(q_IMnpipi_wSid_n_mc,-1.0);
  q_IMnpipi_wSid_n_sub->Draw(opt);




  TCanvas *cIMnpipi_wSid_n = new TCanvas("cIMnpipi_wSid_n","cIMnpipi_wSid_n");
  TH1D* IMnpipi_wSid_n_rdata = (TH1D*)q_IMnpipi_wSid_n[0]->ProjectionX("IMnpipi_wSid_n_rdata");
  TH1D* IMnpipi_wSid_n_woK0 = (TH1D*)q_IMnpipi_wSid_n[1]->ProjectionX("IMnpipi_wSid_n_woK0");
  TH1D* IMnpipi_wSid_n_wK0 = (TH1D*)q_IMnpipi_wSid_n[2]->ProjectionX("IMnpipi_wSid_n_wK0");
  
  IMnpipi_wSid_n_rdata->Draw("HE");
  TH1D* IMnpipi_wSid_n_mc = (TH1D*)IMnpipi_wSid_n_woK0->Clone("IMnpipi_wSid_n_mc");
  IMnpipi_wSid_n_mc->Add(IMnpipi_wSid_n_wK0);
  IMnpipi_wSid_n_mc->SetLineColor(6);
  IMnpipi_wSid_n_mc->Draw("HEsame");
  IMnpipi_wSid_n_woK0->SetLineColor(2);
  IMnpipi_wSid_n_woK0->Draw("HEsame");
  IMnpipi_wSid_n_wK0->SetLineColor(3);
  IMnpipi_wSid_n_wK0->Draw("HEsame");

  TCanvas *cq_wSid_n = new TCanvas("cq_wSid_n","cq_wSid_n");
  TH1D* q_wSid_n_rdata = (TH1D*)q_IMnpipi_wSid_n[0]->ProjectionY("q_wSid_n_rdata");
  TH1D* q_wSid_n_woK0 = (TH1D*)q_IMnpipi_wSid_n[1]->ProjectionY("q_wSid_n_woK0");
  TH1D* q_wSid_n_wK0 = (TH1D*)q_IMnpipi_wSid_n[2]->ProjectionY("q_wSid_n_wK0");
  q_wSid_n_rdata->Draw("HE");
  TH1D* q_wSid_n_mc = (TH1D*)q_wSid_n_woK0->Clone("q_wSid_n_mc");
  q_wSid_n_mc->Add(q_wSid_n_wK0);
  q_wSid_n_mc->SetLineColor(6);
  q_wSid_n_mc->Draw("HEsame");
  q_wSid_n_woK0->SetLineColor(2);
  q_wSid_n_woK0->Draw("HEsame");
  q_wSid_n_wK0->SetLineColor(3);
  q_wSid_n_wK0->Draw("HEsame");



  TCanvas *cIMnpipi_wSid_n_0 = new TCanvas("cIMnpipi_wSid_n_0","cIMnpipi_wSid_n_0");
  TH1D* IMnpipi_wSid_n_0[3];
  for(int i=0;i<3;i++)IMnpipi_wSid_n_0[i] = q_IMnpipi_wSid_n[i]->ProjectionX(Form("IMnpipi_wSid_n_0_%s",name[i]),0,q350bin-1);
  IMnpipi_wSid_n_0[0]->Draw("HE");
  TH1D* IMnpipi_wSid_n_mc_0 = q_IMnpipi_wSid_n_mc->ProjectionX("IMnpipi_wSid_n_mc_0",0,q350bin-1);
  IMnpipi_wSid_n_mc_0->SetLineColor(6);
  IMnpipi_wSid_n_mc_0->Draw("HEsame");
  IMnpipi_wSid_n_0[1]->SetLineColor(2);
  IMnpipi_wSid_n_0[1]->Draw("HEsame");
  IMnpipi_wSid_n_0[2]->SetLineColor(3);
  IMnpipi_wSid_n_0[2]->Draw("HEsame");

  TCanvas *cIMnpipi_wSid_n_350 = new TCanvas("cIMnpipi_wSid_n_350","cIMnpipi_wSid_n_350");
  TH1D* IMnpipi_wSid_n_350[3];
  for(int i=0;i<3;i++)IMnpipi_wSid_n_350[i] = q_IMnpipi_wSid_n[i]->ProjectionX(Form("IMnpipi_wSid_n_350_%s",name[i]),q350bin,q600bin);
  IMnpipi_wSid_n_350[0]->Draw("HE");
  TH1D* IMnpipi_wSid_n_mc_350 = q_IMnpipi_wSid_n_mc->ProjectionX("IMnpipi_wSid_n_mc_350",q350bin,q600bin);
  IMnpipi_wSid_n_mc_350->SetLineColor(6);
  IMnpipi_wSid_n_mc_350->Draw("HEsame");
  IMnpipi_wSid_n_350[1]->SetLineColor(2);
  IMnpipi_wSid_n_350[1]->Draw("HEsame");
  IMnpipi_wSid_n_350[2]->SetLineColor(3);
  IMnpipi_wSid_n_350[2]->Draw("HEsame");


  TCanvas *c = NULL;
  TSeqCollection *SCol = gROOT->GetListOfCanvases();
  int size = SCol->GetSize();
  TIter next(SCol);
  TString pdfname="comp_fakedata_out_v2.pdf";
  for(int i=0; i<size; i++) {
    c= (TCanvas*)next();
    c->Draw();
    c->cd();
    if(i==0) c->Print(pdfname+"(",Form("Title:%s",c->GetTitle()));
    else if(i==size-1)c->Print(pdfname+")",Form("Title:%s",c->GetTitle()));
    else c->Print(pdfname,Form("Title:%s",c->GetTitle()));
  }

  TIter nexthist2(gDirectory->GetList());

  TFile *fout = new TFile("comp_fakedata_out.root","RECREATE");
  TObject *obj = nullptr;
  while( (obj = (TObject*)nexthist2())!=NULL ) {
    obj->Write();
  }
  
  fout->Print();
  fout->cd();
  fout->Close();
}

