//plot_miss2D
//purpose : plot Missing mass, IM(npip) and so on. BG fit by modifying MC data.
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TFile.h"
#include "TCanvas.h"

Double_t func_MMnmiss(Double_t *x,Double_t *par)
{
  if(x[0]<1.116) {
    return par[0]*exp(-0.5*pow((x[0]-par[1])/(par[2]+(x[0]<par[1])*par[3]*(x[0]-par[1])),2.0));
  } else if(1.116<=x[0] && x[0]<1.5) {
    return par[4]*exp(-0.5*pow((x[0]-par[5])/par[6],2.0));
  } else {
    return 0.;
  }
}



Double_t func_IMnpip(Double_t *x,Double_t *par)
{
  if(x[0]<1.08) {
    return 0.;
  } else if(1.08 <= x[0] && x[0]<1.10) {
    return par[0]*exp(-0.5*pow((x[0]-par[1])/par[2],2.0));
  } else if(1.10 <= x[0] && x[0]<1.25) {
    return exp(par[3]+par[4]*x[0]);
  } else if(1.25 <= x[0] && x[0]<2.0) {
    return par[5]*exp(-0.5*pow((x[0]-par[6])/par[7],2.0));
  } else {
    return 0.;
  }
}

Double_t func_IMnpim(Double_t *x,Double_t *par)
{
  if(x[0]<1.07) {
    return 0.;
  } else if(1.07 <= x[0] && x[0]<1.10) {
    return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0);
  } else if(1.10 <= x[0] && x[0]<2.00) {
    return par[4]+par[5]*x[0]+par[6]*pow(x[0],2.0)+par[7]*pow(x[0],3.0);
  } else {
    return 0.;
  }
}

Double_t func_cosn(Double_t *x,Double_t *p)
{
  if(x[0]<-0.9){
    return TMath::Exp(p[0]+p[1]*x[0]);
  }else{
    return (p[2]+p[3]*x[0]+p[4]*TMath::Power(x[0],2)+p[5]*TMath::Power(x[0],3)) ; 
  }
}


Double_t func_cospip(Double_t *x,Double_t *p)
{
  if(x[0]<-0.76){
    return TMath::Exp(p[0]+p[1]*x[0]);
  }else{
    return (p[2]+p[3]*x[0]+p[4]*TMath::Power(x[0],2)+p[5]*TMath::Power(x[0],3)) ;
  }
}



void plot_miss2D()
{
  TH1::SetDefaultSumw2();
  //rdata,nSp,nSm,Lambda,S0pi+pi-,K0nn,mcsum
  const unsigned int colordef[7]= {1,2,3,4,5,28,6};
  const char name[][10]= {"rdata","Sigmap","Sigman","Lambda","S0","K0nn"};

  //real data
  TFile *filerdata = TFile::Open("../post/evanaIMpisigma_npippim_v196_out.root","READ");
  //Sp
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sp_rdata = (TH2D*) filerdata->Get("MMnmiss_IMnpipi_woK0_wSid_Sp");
  double  nrdata_Sp = MMnmiss_IMnpipi_woK0_wSid_Sp_rdata->GetEntries();
  //Sm
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sm_rdata = (TH2D*) filerdata->Get("MMnmiss_IMnpipi_woK0_wSid_Sm");
  double nrdata_Sm = MMnmiss_IMnpipi_woK0_wSid_Sm_rdata->GetEntries();

  //Lambda pi+ pi- n sim.
  TFile *fileLambda = TFile::Open("simIMpisigma_npipiL_pippimn_v23_out.root","READ");
  //Sp
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sp_Lambda = (TH2D*) fileLambda->Get("MMnmiss_IMnpipi_woK0_wSid_Sp");
  double nLambda_Sp = MMnmiss_IMnpipi_woK0_wSid_Sp_Lambda->GetEntries();
  //Sm
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sm_Lambda = (TH2D*) fileLambda->Get("MMnmiss_IMnpipi_woK0_wSid_Sm");
  double nLambda_Sm = MMnmiss_IMnpipi_woK0_wSid_Sm_Lambda->GetEntries();


  //Sigma+ mode sim. ,flat dist in q vs IMnpi+pi-
  TFile *filenSp = TFile::Open("simIMpisigma_nSppim_pippimn_v108_out.root","READ");
  //Sp
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sp_Sigmap = (TH2D*) filenSp->Get("MMnmiss_IMnpipi_woK0_wSid_Sp");
  double nSigmap_Sp = MMnmiss_IMnpipi_woK0_wSid_Sp_Sigmap->GetEntries();
  //Sm
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sm_Sigmap = (TH2D*) filenSp->Get("MMnmiss_IMnpipi_woK0_wSid_Sm");
  double nSigmap_Sm = MMnmiss_IMnpipi_woK0_wSid_Sm_Sigmap->GetEntries();


  //Sigma- mode sim. flat. dist in q vs IMnpi+pi-
  TFile *filenSm = TFile::Open("simIMpisigma_nSmpip_pippimn_v108_out.root","READ");
  //Sp
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sp_Sigmam = (TH2D*) filenSm->Get("MMnmiss_IMnpipi_woK0_wSid_Sp");
  double nSigmam_Sp = MMnmiss_IMnpipi_woK0_wSid_Sp_Sigmam->GetEntries();
  //Sm
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sm_Sigmam = (TH2D*) filenSm->Get("MMnmiss_IMnpipi_woK0_wSid_Sm");
  double nSigmam_Sm = MMnmiss_IMnpipi_woK0_wSid_Sm_Sigmam->GetEntries();


  //nS0pi+pi-
  TFile *fileS0pipi = TFile::Open("simIMpisigma_nS0pippim_pippimn_v5_out.root","READ");
  //Sp
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sp_S0 = (TH2D*) fileS0pipi->Get("MMnmiss_IMnpipi_woK0_wSid_Sp");
  double nS0miss_Sp = MMnmiss_IMnpipi_woK0_wSid_Sp_S0->GetEntries();
  //Sm
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sm_S0 = (TH2D*) fileS0pipi->Get("MMnmiss_IMnpipi_woK0_wSid_Sm");
  double nS0miss_Sm = MMnmiss_IMnpipi_woK0_wSid_Sm_S0->GetEntries();

  //K0nn
  TFile *fileK0nn = TFile::Open("simIMpisigma_K0nn_pippimn_v8_out.root","READ");
  //Sp
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sp_K0nn = (TH2D*) fileK0nn->Get("MMnmiss_IMnpipi_woK0_wSid_Sp");
  double K0nn_Sp = MMnmiss_IMnpipi_woK0_wSid_Sp_K0nn->GetEntries();
  //Sm
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sm_K0nn = (TH2D*) fileK0nn->Get("MMnmiss_IMnpipi_woK0_wSid_Sm");
  double K0nn_Sm = MMnmiss_IMnpipi_woK0_wSid_Sm_K0nn->GetEntries();



  //superimpose 2d hists of miss mass vs IM(npi+) or IM(npi-)
  //Lambda
  TH2D *MMnmiss_IMnpip_woK0_Lambda = (TH2D*)fileLambda->Get("MMnmiss_IMnpip_dE_woK0");
  TH2D *MMnmiss_IMnpim_woK0_Lambda = (TH2D*)fileLambda->Get("MMnmiss_IMnpim_dE_woK0");
  TH2D *MMnmiss_IMnpip_woK0_woSm_Lambda = (TH2D*)fileLambda->Get("MMnmiss_IMnpip_dE_woK0_woSm");
  TH2D *MMnmiss_IMnpim_woK0_woSp_Lambda = (TH2D*)fileLambda->Get("MMnmiss_IMnpim_dE_woK0_woSp");
  TH2D *MMnmiss_IMnpip_woK0_woSm_n_Lambda = (TH2D*)fileLambda->Get("MMnmiss_IMnpip_dE_woK0_woSm_n");
  TH2D *MMnmiss_IMnpim_woK0_woSp_n_Lambda = (TH2D*)fileLambda->Get("MMnmiss_IMnpim_dE_woK0_woSp_n");
  TH2D *MMnmiss_IMnpip_woK0_woSidn_Lambda = (TH2D*)fileLambda->Get("MMnmiss_IMnpip_dE_woK0_woSidn");
  TH2D *MMnmiss_IMnpim_woK0_woSidn_Lambda = (TH2D*)fileLambda->Get("MMnmiss_IMnpim_dE_woK0_woSidn");
  TH2D *MMnmiss_IMnpip_wK0_woSid_won_Lambda = (TH2D*)fileLambda->Get("MMnmiss_IMnpip_dE_wK0_woSid_won");
  TH2D *MMnmiss_IMnpim_wK0_woSid_won_Lambda = (TH2D*)fileLambda->Get("MMnmiss_IMnpim_dE_wK0_woSid_won");
  TH2D *MMnmiss_IMpippim_Lambda = (TH2D*)fileLambda->Get("MMnmiss_IMpippim_dE");
  TH2D *MMnmiss_IMpippim_wSid_Lambda = (TH2D*)fileLambda->Get("MMnmiss_IMpippim_dE_wSid");
  TH2D *MMnmiss_IMpippim_woK0_woSidn_Lambda = (TH2D*)fileLambda->Get("MMnmiss_IMpippim_dE_woK0_woSidn");
  TH2D *MMnmiss_IMpippim_wK0_woSid_won_Lambda = (TH2D*)fileLambda->Get("MMnmiss_IMpippim_dE_wK0_woSid_won");
  TH2D *IMnpim_IMnpip_Lambda = (TH2D*)fileLambda->Get("IMnpim_IMnpip_dE");
  TH2D *IMnpim_IMnpip_woK0_woSidn_Lambda = (TH2D*)fileLambda->Get("IMnpim_IMnpip_dE_woK0_woSidn");
  TH2D *q_IMnpipi_woK0_woSidn_Lambda = (TH2D*)fileLambda->Get("q_IMnpipi_woK0_woSidn");
  TH2D *q_IMnpipi_wK0_woSid_won_Lambda = (TH2D*)fileLambda->Get("q_IMnpipi_wK0_woSid_won");
  TH2D *MMnmiss_Mompippim_woK0_woSidn_Lambda  = (TH2D*)fileLambda->Get("MMnmiss_Mompippim_dE_woK0_woSidn");
  TH2D *MMnmiss_Mompippim_wK0_woSid_won_Lambda = (TH2D*)fileLambda->Get("MMnmiss_Mompippim_dE_wK0_woSid_won");
  TH2D *pipmom_MMnmiss_woK0_woSidn_Lambda = (TH2D*)fileLambda->Get("pipmom_MMnmiss_dE_woK0_woSidn");
  TH2D *pipmom_MMnmiss_wK0_woSid_won_Lambda = (TH2D*)fileLambda->Get("pipmom_MMnmiss_dE_wK0_woSid_won");
  TH2D *pipmom_pimmom_woK0_woSidn_Lambda = (TH2D*)fileLambda->Get("pipmom_pimmom_dE_woK0_woSidn");
  TH2D *pipmom_pimmom_wK0_woSid_won_Lambda = (TH2D*)fileLambda->Get("pipmom_pimmom_dE_wK0_woSid_won");
  TH2D *pimmom_MMnmiss_woK0_woSidn_Lambda = (TH2D*)fileLambda->Get("pimmom_MMnmiss_dE_woK0_woSidn");
  TH2D *pimmom_MMnmiss_wK0_woSid_won_Lambda = (TH2D*)fileLambda->Get("pimmom_MMnmiss_dE_wK0_woSid_won");
  TH2D *MMom_MMass_woK0_woSidn_Lambda = (TH2D*)fileLambda->Get("MMom_MMass_woK0_woSidn");
  TH2D *MMom_MMass_wK0_woSid_won_Lambda = (TH2D*)fileLambda->Get("MMom_MMass_wK0_woSid_won");
  TH2D *nmom_MMnmiss_woK0_woSidn_Lambda = (TH2D*)fileLambda->Get("nmom_MMnmiss_woK0_woSidn");
  TH2D *nmom_MMnmiss_wK0_woSid_won_Lambda = (TH2D*)fileLambda->Get("nmom_MMnmiss_wK0_woSid_won");
  TH2D *q_IMnpipi_wSid_n_Sp_Lambda = (TH2D*)fileLambda->Get("q_IMnpipi_wSid_n_Sp");
  TH2D *q_IMnpipi_wSid_n_Sm_Lambda = (TH2D*)fileLambda->Get("q_IMnpipi_wSid_n_Sm");
  TH2D *nmom_cosn_woK0_woSidn_Lambda = (TH2D*)fileLambda->Get("nmom_cosn_woK0_woSidn");
  TH2D *nmom_cospip_woK0_woSidn_Lambda = (TH2D*)fileLambda->Get("nmom_cospip_woK0_woSidn");
  TH2D *nmom_cospim_woK0_woSidn_Lambda = (TH2D*)fileLambda->Get("nmom_cospim_woK0_woSidn");
  TH2D *nmom_phinpip_woK0_woSidn_Lambda = (TH2D*)fileLambda->Get("nmom_phinpip_woK0_woSidn");
  TH2D *nmom_phinpim_woK0_woSidn_Lambda = (TH2D*)fileLambda->Get("nmom_phinpim_woK0_woSidn");
  TH2D *nmom_pipmom_woK0_woSidn_Lambda = (TH2D*)fileLambda->Get("nmom_pipmom_woK0_woSidn");
  TH2D *nmom_pimmom_woK0_woSidn_Lambda = (TH2D*)fileLambda->Get("nmom_pimmom_woK0_woSidn");
  TH2D *nmom_cosn_wK0_woSid_won_Lambda = (TH2D*)fileLambda->Get("nmom_cosn_wK0_woSid_won");
  TH2D *nmom_cospip_wK0_woSid_won_Lambda = (TH2D*)fileLambda->Get("nmom_cospip_wK0_woSid_won");
  TH2D *nmom_cospim_wK0_woSid_won_Lambda = (TH2D*)fileLambda->Get("nmom_cospim_wK0_woSid_won");
  TH2D *nmom_phinpip_wK0_woSid_won_Lambda = (TH2D*)fileLambda->Get("nmom_phinpip_wK0_woSid_won");
  TH2D *nmom_phinpim_wK0_woSid_won_Lambda = (TH2D*)fileLambda->Get("nmom_phinpim_wK0_woSid_won");
  TH2D *nmom_pipmom_wK0_woSid_won_Lambda = (TH2D*)fileLambda->Get("nmom_pipmom_wK0_woSid_won");
  TH2D *nmom_pimmom_wK0_woSid_won_Lambda = (TH2D*)fileLambda->Get("nmom_pimmom_wK0_woSid_won");

  double Nnpip_Lambda = MMnmiss_IMnpip_woK0_Lambda->GetEntries();
  double Nnpim_Lambda = MMnmiss_IMnpim_woK0_Lambda->GetEntries();
  //nSppi-
  TH2D *MMnmiss_IMnpip_woK0_Sigmap = (TH2D*)filenSp->Get("MMnmiss_IMnpip_dE_woK0");
  TH2D *MMnmiss_IMnpim_woK0_Sigmap = (TH2D*)filenSp->Get("MMnmiss_IMnpim_dE_woK0");
  TH2D *MMnmiss_IMnpip_woK0_woSm_Sigmap = (TH2D*)filenSp->Get("MMnmiss_IMnpip_dE_woK0_woSm");
  TH2D *MMnmiss_IMnpim_woK0_woSp_Sigmap = (TH2D*)filenSp->Get("MMnmiss_IMnpim_dE_woK0_woSp");
  TH2D *MMnmiss_IMnpip_woK0_woSm_n_Sigmap = (TH2D*)filenSp->Get("MMnmiss_IMnpip_dE_woK0_woSm_n");
  TH2D *MMnmiss_IMnpim_woK0_woSp_n_Sigmap = (TH2D*)filenSp->Get("MMnmiss_IMnpim_dE_woK0_woSp_n");
  TH2D *MMnmiss_IMnpip_woK0_woSidn_Sigmap = (TH2D*)filenSp->Get("MMnmiss_IMnpip_dE_woK0_woSidn");
  TH2D *MMnmiss_IMnpim_woK0_woSidn_Sigmap = (TH2D*)filenSp->Get("MMnmiss_IMnpim_dE_woK0_woSidn");
  TH2D *MMnmiss_IMnpip_wK0_woSid_won_Sigmap = (TH2D*)filenSp->Get("MMnmiss_IMnpip_dE_wK0_woSid_won");
  TH2D *MMnmiss_IMnpim_wK0_woSid_won_Sigmap = (TH2D*)filenSp->Get("MMnmiss_IMnpim_dE_wK0_woSid_won");
  TH2D *MMnmiss_IMpippim_Sigmap = (TH2D*)filenSp->Get("MMnmiss_IMpippim_dE");
  TH2D *MMnmiss_IMpippim_wSid_Sigmap = (TH2D*)filenSp->Get("MMnmiss_IMpippim_dE_wSid");
  TH2D *MMnmiss_IMpippim_woK0_woSidn_Sigmap = (TH2D*)filenSp->Get("MMnmiss_IMpippim_dE_woK0_woSidn");
  TH2D *MMnmiss_IMpippim_wK0_woSid_won_Sigmap = (TH2D*)filenSp->Get("MMnmiss_IMpippim_dE_wK0_woSid_won");
  TH2D *IMnpim_IMnpip_Sigmap = (TH2D*)filenSp->Get("IMnpim_IMnpip_dE");
  TH2D *IMnpim_IMnpip_woK0_woSidn_Sigmap = (TH2D*)filenSp->Get("IMnpim_IMnpip_dE_woK0_woSidn");
  TH2D *q_IMnpipi_woK0_woSidn_Sigmap = (TH2D*)filenSp->Get("q_IMnpipi_woK0_woSidn");
  TH2D *q_IMnpipi_wK0_woSid_won_Sigmap = (TH2D*)filenSp->Get("q_IMnpipi_wK0_woSid_won");
  TH2D *MMnmiss_Mompippim_woK0_woSidn_Sigmap  = (TH2D*)filenSp->Get("MMnmiss_Mompippim_dE_woK0_woSidn");
  TH2D *MMnmiss_Mompippim_wK0_woSid_won_Sigmap = (TH2D*)filenSp->Get("MMnmiss_Mompippim_dE_wK0_woSid_won");
  TH2D *pipmom_MMnmiss_woK0_woSidn_Sigmap = (TH2D*)filenSp->Get("pipmom_MMnmiss_dE_woK0_woSidn");
  TH2D *pipmom_MMnmiss_wK0_woSid_won_Sigmap = (TH2D*)filenSp->Get("pipmom_MMnmiss_dE_wK0_woSid_won");
  TH2D *pipmom_pimmom_woK0_woSidn_Sigmap = (TH2D*)filenSp->Get("pipmom_pimmom_dE_woK0_woSidn");
  TH2D *pipmom_pimmom_wK0_woSid_won_Sigmap = (TH2D*)filenSp->Get("pipmom_pimmom_dE_wK0_woSid_won");
  TH2D *pimmom_MMnmiss_woK0_woSidn_Sigmap = (TH2D*)filenSp->Get("pimmom_MMnmiss_dE_woK0_woSidn");
  TH2D *pimmom_MMnmiss_wK0_woSid_won_Sigmap = (TH2D*)filenSp->Get("pimmom_MMnmiss_dE_wK0_woSid_won");
  TH2D *MMom_MMass_woK0_woSidn_Sigmap = (TH2D*)filenSp->Get("MMom_MMass_woK0_woSidn");
  TH2D *MMom_MMass_wK0_woSid_won_Sigmap = (TH2D*)filenSp->Get("MMom_MMass_wK0_woSid_won");
  TH2D *nmom_MMnmiss_woK0_woSidn_Sigmap = (TH2D*)filenSp->Get("nmom_MMnmiss_woK0_woSidn");
  TH2D *nmom_MMnmiss_wK0_woSid_won_Sigmap = (TH2D*)filenSp->Get("nmom_MMnmiss_wK0_woSid_won");
  TH2D *q_IMnpipi_wSid_n_Sp_Sigmap = (TH2D*)filenSp->Get("q_IMnpipi_wSid_n_Sp");
  TH2D *q_IMnpipi_wSid_n_Sm_Sigmap = (TH2D*)filenSp->Get("q_IMnpipi_wSid_n_Sm");
  TH2D *nmom_cosn_woK0_woSidn_Sigmap = (TH2D*)filenSp->Get("nmom_cosn_woK0_woSidn");
  TH2D *nmom_cospip_woK0_woSidn_Sigmap = (TH2D*)filenSp->Get("nmom_cospip_woK0_woSidn");
  TH2D *nmom_cospim_woK0_woSidn_Sigmap = (TH2D*)filenSp->Get("nmom_cospim_woK0_woSidn");
  TH2D *nmom_phinpip_woK0_woSidn_Sigmap = (TH2D*)filenSp->Get("nmom_phinpip_woK0_woSidn");
  TH2D *nmom_phinpim_woK0_woSidn_Sigmap = (TH2D*)filenSp->Get("nmom_phinpim_woK0_woSidn");
  TH2D *nmom_pipmom_woK0_woSidn_Sigmap = (TH2D*)filenSp->Get("nmom_pipmom_woK0_woSidn");
  TH2D *nmom_pimmom_woK0_woSidn_Sigmap = (TH2D*)filenSp->Get("nmom_pimmom_woK0_woSidn");
  TH2D *nmom_cosn_wK0_woSid_won_Sigmap = (TH2D*)filenSp->Get("nmom_cosn_wK0_woSid_won");
  TH2D *nmom_cospip_wK0_woSid_won_Sigmap = (TH2D*)filenSp->Get("nmom_cospip_wK0_woSid_won");
  TH2D *nmom_cospim_wK0_woSid_won_Sigmap = (TH2D*)filenSp->Get("nmom_cospim_wK0_woSid_won");
  TH2D *nmom_phinpip_wK0_woSid_won_Sigmap = (TH2D*)filenSp->Get("nmom_phinpip_wK0_woSid_won");
  TH2D *nmom_phinpim_wK0_woSid_won_Sigmap = (TH2D*)filenSp->Get("nmom_phinpim_wK0_woSid_won");
  TH2D *nmom_pipmom_wK0_woSid_won_Sigmap = (TH2D*)filenSp->Get("nmom_pipmom_wK0_woSid_won");
  TH2D *nmom_pimmom_wK0_woSid_won_Sigmap = (TH2D*)filenSp->Get("nmom_pimmom_wK0_woSid_won");
  double Nnpip_Sigmap = MMnmiss_IMnpip_woK0_Sigmap->GetEntries();
  double Nnpim_Sigmap = MMnmiss_IMnpim_woK0_Sigmap->GetEntries();
  //nSmpi+
  TH2D *MMnmiss_IMnpip_woK0_Sigmam = (TH2D*)filenSm->Get("MMnmiss_IMnpip_dE_woK0");
  TH2D *MMnmiss_IMnpim_woK0_Sigmam = (TH2D*)filenSm->Get("MMnmiss_IMnpim_dE_woK0");
  TH2D *MMnmiss_IMnpip_woK0_woSm_Sigmam = (TH2D*)filenSm->Get("MMnmiss_IMnpip_dE_woK0_woSm");
  TH2D *MMnmiss_IMnpim_woK0_woSp_Sigmam = (TH2D*)filenSm->Get("MMnmiss_IMnpim_dE_woK0_woSp");
  TH2D *MMnmiss_IMnpip_woK0_woSm_n_Sigmam = (TH2D*)filenSm->Get("MMnmiss_IMnpip_dE_woK0_woSm_n");
  TH2D *MMnmiss_IMnpim_woK0_woSp_n_Sigmam = (TH2D*)filenSm->Get("MMnmiss_IMnpim_dE_woK0_woSp_n");
  TH2D *MMnmiss_IMnpip_woK0_woSidn_Sigmam = (TH2D*)filenSm->Get("MMnmiss_IMnpip_dE_woK0_woSidn");
  TH2D *MMnmiss_IMnpim_woK0_woSidn_Sigmam = (TH2D*)filenSm->Get("MMnmiss_IMnpim_dE_woK0_woSidn");
  TH2D *MMnmiss_IMnpip_wK0_woSid_won_Sigmam = (TH2D*)filenSm->Get("MMnmiss_IMnpip_dE_wK0_woSid_won");
  TH2D *MMnmiss_IMnpim_wK0_woSid_won_Sigmam = (TH2D*)filenSm->Get("MMnmiss_IMnpim_dE_wK0_woSid_won");
  TH2D *MMnmiss_IMpippim_Sigmam = (TH2D*)filenSm->Get("MMnmiss_IMpippim_dE");
  TH2D *MMnmiss_IMpippim_wSid_Sigmam = (TH2D*)filenSm->Get("MMnmiss_IMpippim_dE_wSid");
  TH2D *MMnmiss_IMpippim_woK0_woSidn_Sigmam = (TH2D*)filenSm->Get("MMnmiss_IMpippim_dE_woK0_woSidn");
  TH2D *MMnmiss_IMpippim_wK0_woSid_won_Sigmam = (TH2D*)filenSm->Get("MMnmiss_IMpippim_dE_wK0_woSid_won");
  TH2D *IMnpim_IMnpip_Sigmam = (TH2D*)filenSm->Get("IMnpim_IMnpip_dE");
  TH2D *IMnpim_IMnpip_woK0_woSidn_Sigmam = (TH2D*)filenSm->Get("IMnpim_IMnpip_dE_woK0_woSidn");
  TH2D *q_IMnpipi_woK0_woSidn_Sigmam = (TH2D*)filenSm->Get("q_IMnpipi_woK0_woSidn");
  TH2D *q_IMnpipi_wK0_woSid_won_Sigmam = (TH2D*)filenSm->Get("q_IMnpipi_wK0_woSid_won");
  TH2D *MMnmiss_Mompippim_woK0_woSidn_Sigmam  = (TH2D*)filenSm->Get("MMnmiss_Mompippim_dE_woK0_woSidn");
  TH2D *MMnmiss_Mompippim_wK0_woSid_won_Sigmam = (TH2D*)filenSm->Get("MMnmiss_Mompippim_dE_wK0_woSid_won");
  TH2D *pipmom_MMnmiss_woK0_woSidn_Sigmam = (TH2D*)filenSm->Get("pipmom_MMnmiss_dE_woK0_woSidn");
  TH2D *pipmom_MMnmiss_wK0_woSid_won_Sigmam = (TH2D*)filenSm->Get("pipmom_MMnmiss_dE_wK0_woSid_won");
  TH2D *pipmom_pimmom_woK0_woSidn_Sigmam = (TH2D*)filenSm->Get("pipmom_pimmom_dE_woK0_woSidn");
  TH2D *pipmom_pimmom_wK0_woSid_won_Sigmam = (TH2D*)filenSm->Get("pipmom_pimmom_dE_wK0_woSid_won");
  TH2D *pimmom_MMnmiss_woK0_woSidn_Sigmam = (TH2D*)filenSm->Get("pimmom_MMnmiss_dE_woK0_woSidn");
  TH2D *pimmom_MMnmiss_wK0_woSid_won_Sigmam = (TH2D*)filenSm->Get("pimmom_MMnmiss_dE_wK0_woSid_won");
  TH2D *MMom_MMass_woK0_woSidn_Sigmam = (TH2D*)filenSm->Get("MMom_MMass_woK0_woSidn");
  TH2D *MMom_MMass_wK0_woSid_won_Sigmam = (TH2D*)filenSm->Get("MMom_MMass_wK0_woSid_won");
  TH2D *nmom_MMnmiss_woK0_woSidn_Sigmam = (TH2D*)filenSm->Get("nmom_MMnmiss_woK0_woSidn");
  TH2D *nmom_MMnmiss_wK0_woSid_won_Sigmam = (TH2D*)filenSm->Get("nmom_MMnmiss_wK0_woSid_won");
  TH2D *q_IMnpipi_wSid_n_Sp_Sigmam = (TH2D*)filenSm->Get("q_IMnpipi_wSid_n_Sp");
  TH2D *q_IMnpipi_wSid_n_Sm_Sigmam = (TH2D*)filenSm->Get("q_IMnpipi_wSid_n_Sm");
  TH2D *nmom_cosn_woK0_woSidn_Sigmam = (TH2D*)filenSm->Get("nmom_cosn_woK0_woSidn");
  TH2D *nmom_cospip_woK0_woSidn_Sigmam = (TH2D*)filenSm->Get("nmom_cospip_woK0_woSidn");
  TH2D *nmom_cospim_woK0_woSidn_Sigmam = (TH2D*)filenSm->Get("nmom_cospim_woK0_woSidn");
  TH2D *nmom_phinpip_woK0_woSidn_Sigmam = (TH2D*)filenSm->Get("nmom_phinpip_woK0_woSidn");
  TH2D *nmom_phinpim_woK0_woSidn_Sigmam = (TH2D*)filenSm->Get("nmom_phinpim_woK0_woSidn");
  TH2D *nmom_pipmom_woK0_woSidn_Sigmam = (TH2D*)filenSm->Get("nmom_pipmom_woK0_woSidn");
  TH2D *nmom_pimmom_woK0_woSidn_Sigmam = (TH2D*)filenSm->Get("nmom_pimmom_woK0_woSidn");
  TH2D *nmom_cosn_wK0_woSid_won_Sigmam = (TH2D*)filenSm->Get("nmom_cosn_wK0_woSid_won");
  TH2D *nmom_cospip_wK0_woSid_won_Sigmam = (TH2D*)filenSm->Get("nmom_cospip_wK0_woSid_won");
  TH2D *nmom_cospim_wK0_woSid_won_Sigmam = (TH2D*)filenSm->Get("nmom_cospim_wK0_woSid_won");
  TH2D *nmom_phinpip_wK0_woSid_won_Sigmam = (TH2D*)filenSm->Get("nmom_phinpip_wK0_woSid_won");
  TH2D *nmom_phinpim_wK0_woSid_won_Sigmam = (TH2D*)filenSm->Get("nmom_phinpim_wK0_woSid_won");
  TH2D *nmom_pipmom_wK0_woSid_won_Sigmam = (TH2D*)filenSm->Get("nmom_pipmom_wK0_woSid_won");
  TH2D *nmom_pimmom_wK0_woSid_won_Sigmam = (TH2D*)filenSm->Get("nmom_pimmom_wK0_woSid_won");
  double Nnpip_Sigmam = MMnmiss_IMnpip_woK0_Sigmam->GetEntries();
  double Nnpim_Sigmam = MMnmiss_IMnpim_woK0_Sigmam->GetEntries();
  //Sigma0
  TH2D *MMnmiss_IMnpip_woK0_S0 = (TH2D*)fileS0pipi->Get("MMnmiss_IMnpip_dE_woK0");
  TH2D *MMnmiss_IMnpim_woK0_S0 = (TH2D*)fileS0pipi->Get("MMnmiss_IMnpim_dE_woK0");
  TH2D *MMnmiss_IMnpip_woK0_woSm_S0 = (TH2D*)fileS0pipi->Get("MMnmiss_IMnpip_dE_woK0_woSm");
  TH2D *MMnmiss_IMnpim_woK0_woSp_S0 = (TH2D*)fileS0pipi->Get("MMnmiss_IMnpim_dE_woK0_woSp");
  TH2D *MMnmiss_IMnpip_woK0_woSm_n_S0 = (TH2D*)fileS0pipi->Get("MMnmiss_IMnpip_dE_woK0_woSm_n");
  TH2D *MMnmiss_IMnpim_woK0_woSp_n_S0 = (TH2D*)fileS0pipi->Get("MMnmiss_IMnpim_dE_woK0_woSp_n");
  TH2D *MMnmiss_IMnpip_woK0_woSidn_S0 = (TH2D*)fileS0pipi->Get("MMnmiss_IMnpip_dE_woK0_woSidn");
  TH2D *MMnmiss_IMnpim_woK0_woSidn_S0 = (TH2D*)fileS0pipi->Get("MMnmiss_IMnpim_dE_woK0_woSidn");
  TH2D *MMnmiss_IMnpip_wK0_woSid_won_S0 = (TH2D*)fileS0pipi->Get("MMnmiss_IMnpip_dE_wK0_woSid_won");
  TH2D *MMnmiss_IMnpim_wK0_woSid_won_S0 = (TH2D*)fileS0pipi->Get("MMnmiss_IMnpim_dE_wK0_woSid_won");
  TH2D *MMnmiss_IMpippim_S0 = (TH2D*)fileS0pipi->Get("MMnmiss_IMpippim_dE");
  TH2D *MMnmiss_IMpippim_wSid_S0 = (TH2D*)fileS0pipi->Get("MMnmiss_IMpippim_dE_wSid");
  TH2D *MMnmiss_IMpippim_woK0_woSidn_S0 = (TH2D*)fileS0pipi->Get("MMnmiss_IMpippim_dE_woK0_woSidn");
  TH2D *MMnmiss_IMpippim_wK0_woSid_won_S0 = (TH2D*)fileS0pipi->Get("MMnmiss_IMpippim_dE_wK0_woSid_won");
  TH2D *IMnpim_IMnpip_S0 = (TH2D*)fileS0pipi->Get("IMnpim_IMnpip_dE");
  TH2D *IMnpim_IMnpip_woK0_woSidn_S0 = (TH2D*)fileS0pipi->Get("IMnpim_IMnpip_dE_woK0_woSidn");
  TH2D *q_IMnpipi_woK0_woSidn_S0 = (TH2D*)fileS0pipi->Get("q_IMnpipi_woK0_woSidn");
  TH2D *q_IMnpipi_wK0_woSid_won_S0 = (TH2D*)fileS0pipi->Get("q_IMnpipi_wK0_woSid_won");
  TH2D *MMnmiss_Mompippim_woK0_woSidn_S0  = (TH2D*)fileS0pipi->Get("MMnmiss_Mompippim_dE_woK0_woSidn");
  TH2D *MMnmiss_Mompippim_wK0_woSid_won_S0 = (TH2D*)fileS0pipi->Get("MMnmiss_Mompippim_dE_wK0_woSid_won");
  TH2D *pipmom_MMnmiss_woK0_woSidn_S0 = (TH2D*)fileS0pipi->Get("pipmom_MMnmiss_dE_woK0_woSidn");
  TH2D *pipmom_MMnmiss_wK0_woSid_won_S0 = (TH2D*)fileS0pipi->Get("pipmom_MMnmiss_dE_wK0_woSid_won");
  TH2D *pipmom_pimmom_woK0_woSidn_S0 = (TH2D*)fileS0pipi->Get("pipmom_pimmom_dE_woK0_woSidn");
  TH2D *pipmom_pimmom_wK0_woSid_won_S0 = (TH2D*)fileS0pipi->Get("pipmom_pimmom_dE_wK0_woSid_won");
  TH2D *pimmom_MMnmiss_woK0_woSidn_S0 = (TH2D*)fileS0pipi->Get("pimmom_MMnmiss_dE_woK0_woSidn");
  TH2D *pimmom_MMnmiss_wK0_woSid_won_S0 = (TH2D*)fileS0pipi->Get("pimmom_MMnmiss_dE_wK0_woSid_won");
  TH2D *MMom_MMass_woK0_woSidn_S0 = (TH2D*)fileS0pipi->Get("MMom_MMass_woK0_woSidn");
  TH2D *MMom_MMass_wK0_woSid_won_S0 = (TH2D*)fileS0pipi->Get("MMom_MMass_wK0_woSid_won");
  TH2D *nmom_MMnmiss_woK0_woSidn_S0 = (TH2D*)fileS0pipi->Get("nmom_MMnmiss_woK0_woSidn");
  TH2D *nmom_MMnmiss_wK0_woSid_won_S0 = (TH2D*)fileS0pipi->Get("nmom_MMnmiss_wK0_woSid_won");
  TH2D *q_IMnpipi_wSid_n_Sp_S0 = (TH2D*)fileS0pipi->Get("q_IMnpipi_wSid_n_Sp");
  TH2D *q_IMnpipi_wSid_n_Sm_S0 = (TH2D*)fileS0pipi->Get("q_IMnpipi_wSid_n_Sm");
  TH2D *nmom_cosn_woK0_woSidn_S0 = (TH2D*)fileS0pipi->Get("nmom_cosn_woK0_woSidn");
  TH2D *nmom_cospip_woK0_woSidn_S0 = (TH2D*)fileS0pipi->Get("nmom_cospip_woK0_woSidn");
  TH2D *nmom_cospim_woK0_woSidn_S0 = (TH2D*)fileS0pipi->Get("nmom_cospim_woK0_woSidn");
  TH2D *nmom_phinpip_woK0_woSidn_S0 = (TH2D*)fileS0pipi->Get("nmom_phinpip_woK0_woSidn");
  TH2D *nmom_phinpim_woK0_woSidn_S0 = (TH2D*)fileS0pipi->Get("nmom_phinpim_woK0_woSidn");
  TH2D *nmom_pipmom_woK0_woSidn_S0 = (TH2D*)fileS0pipi->Get("nmom_pipmom_woK0_woSidn");
  TH2D *nmom_pimmom_woK0_woSidn_S0 = (TH2D*)fileS0pipi->Get("nmom_pimmom_woK0_woSidn");
  TH2D *nmom_cosn_wK0_woSid_won_S0 = (TH2D*)fileS0pipi->Get("nmom_cosn_wK0_woSid_won");
  TH2D *nmom_cospip_wK0_woSid_won_S0 = (TH2D*)fileS0pipi->Get("nmom_cospip_wK0_woSid_won");
  TH2D *nmom_cospim_wK0_woSid_won_S0 = (TH2D*)fileS0pipi->Get("nmom_cospim_wK0_woSid_won");
  TH2D *nmom_phinpip_wK0_woSid_won_S0 = (TH2D*)fileS0pipi->Get("nmom_phinpip_wK0_woSid_won");
  TH2D *nmom_phinpim_wK0_woSid_won_S0 = (TH2D*)fileS0pipi->Get("nmom_phinpim_wK0_woSid_won");
  TH2D *nmom_pipmom_wK0_woSid_won_S0 = (TH2D*)fileS0pipi->Get("nmom_pipmom_wK0_woSid_won");
  TH2D *nmom_pimmom_wK0_woSid_won_S0 = (TH2D*)fileS0pipi->Get("nmom_pimmom_wK0_woSid_won");
  double Nnpip_S0 = MMnmiss_IMnpip_woK0_S0->GetEntries();
  double Nnpim_S0 = MMnmiss_IMnpim_woK0_S0->GetEntries();

  //K0nn
  TH2D *MMnmiss_IMnpip_woK0_K0nn = (TH2D*)fileK0nn->Get("MMnmiss_IMnpip_dE_woK0");
  TH2D *MMnmiss_IMnpim_woK0_K0nn = (TH2D*)fileK0nn->Get("MMnmiss_IMnpim_dE_woK0");
  TH2D *MMnmiss_IMnpip_woK0_woSm_K0nn = (TH2D*)fileK0nn->Get("MMnmiss_IMnpip_dE_woK0_woSm");
  TH2D *MMnmiss_IMnpim_woK0_woSp_K0nn = (TH2D*)fileK0nn->Get("MMnmiss_IMnpim_dE_woK0_woSp");
  TH2D *MMnmiss_IMnpip_woK0_woSm_n_K0nn = (TH2D*)fileK0nn->Get("MMnmiss_IMnpip_dE_woK0_woSm_n");
  TH2D *MMnmiss_IMnpim_woK0_woSp_n_K0nn = (TH2D*)fileK0nn->Get("MMnmiss_IMnpim_dE_woK0_woSp_n");
  TH2D *MMnmiss_IMnpip_woK0_woSidn_K0nn = (TH2D*)fileK0nn->Get("MMnmiss_IMnpip_dE_woK0_woSidn");
  TH2D *MMnmiss_IMnpim_woK0_woSidn_K0nn = (TH2D*)fileK0nn->Get("MMnmiss_IMnpim_dE_woK0_woSidn");
  TH2D *MMnmiss_IMnpip_wK0_woSid_won_K0nn = (TH2D*)fileK0nn->Get("MMnmiss_IMnpip_dE_wK0_woSid_won");
  TH2D *MMnmiss_IMnpim_wK0_woSid_won_K0nn = (TH2D*)fileK0nn->Get("MMnmiss_IMnpim_dE_wK0_woSid_won");
  TH2D *MMnmiss_IMpippim_K0nn = (TH2D*)fileK0nn->Get("MMnmiss_IMpippim_dE");
  TH2D *MMnmiss_IMpippim_wSid_K0nn = (TH2D*)fileK0nn->Get("MMnmiss_IMpippim_dE_wSid");
  TH2D *MMnmiss_IMpippim_woK0_woSidn_K0nn = (TH2D*)fileK0nn->Get("MMnmiss_IMpippim_dE_woK0_woSidn");
  TH2D *MMnmiss_IMpippim_wK0_woSid_won_K0nn = (TH2D*)fileK0nn->Get("MMnmiss_IMpippim_dE_wK0_woSid_won");
  TH2D *IMnpim_IMnpip_K0nn = (TH2D*)fileK0nn->Get("IMnpim_IMnpip_dE");
  TH2D *IMnpim_IMnpip_woK0_woSidn_K0nn = (TH2D*)fileK0nn->Get("IMnpim_IMnpip_dE_woK0_woSidn");
  TH2D *q_IMnpipi_woK0_woSidn_K0nn = (TH2D*)fileK0nn->Get("q_IMnpipi_woK0_woSidn");
  TH2D *q_IMnpipi_wK0_woSid_won_K0nn = (TH2D*)fileK0nn->Get("q_IMnpipi_wK0_woSid_won");
  TH2D *MMnmiss_Mompippim_woK0_woSidn_K0nn  = (TH2D*)fileK0nn->Get("MMnmiss_Mompippim_dE_woK0_woSidn");
  TH2D *MMnmiss_Mompippim_wK0_woSid_won_K0nn = (TH2D*)fileK0nn->Get("MMnmiss_Mompippim_dE_wK0_woSid_won");
  TH2D *pipmom_MMnmiss_woK0_woSidn_K0nn = (TH2D*)fileK0nn->Get("pipmom_MMnmiss_dE_woK0_woSidn");
  TH2D *pipmom_MMnmiss_wK0_woSid_won_K0nn = (TH2D*)fileK0nn->Get("pipmom_MMnmiss_dE_wK0_woSid_won");
  TH2D *pipmom_pimmom_woK0_woSidn_K0nn = (TH2D*)fileK0nn->Get("pipmom_pimmom_dE_woK0_woSidn");
  TH2D *pipmom_pimmom_wK0_woSid_won_K0nn = (TH2D*)fileK0nn->Get("pipmom_pimmom_dE_wK0_woSid_won");
  TH2D *pimmom_MMnmiss_woK0_woSidn_K0nn = (TH2D*)fileK0nn->Get("pimmom_MMnmiss_dE_woK0_woSidn");
  TH2D *pimmom_MMnmiss_wK0_woSid_won_K0nn = (TH2D*)fileK0nn->Get("pimmom_MMnmiss_dE_wK0_woSid_won");
  TH2D *MMom_MMass_woK0_woSidn_K0nn = (TH2D*)fileK0nn->Get("MMom_MMass_woK0_woSidn");
  TH2D *MMom_MMass_wK0_woSid_won_K0nn = (TH2D*)fileK0nn->Get("MMom_MMass_wK0_woSid_won");
  TH2D *nmom_MMnmiss_woK0_woSidn_K0nn = (TH2D*)fileK0nn->Get("nmom_MMnmiss_woK0_woSidn");
  TH2D *nmom_MMnmiss_wK0_woSid_won_K0nn = (TH2D*)fileK0nn->Get("nmom_MMnmiss_wK0_woSid_won");
  TH2D *q_IMnpipi_wSid_n_Sp_K0nn = (TH2D*)fileK0nn->Get("q_IMnpipi_wSid_n_Sp");
  TH2D *q_IMnpipi_wSid_n_Sm_K0nn = (TH2D*)fileK0nn->Get("q_IMnpipi_wSid_n_Sm");
  TH2D *nmom_cosn_woK0_woSidn_K0nn = (TH2D*)fileK0nn->Get("nmom_cosn_woK0_woSidn");
  TH2D *nmom_cospip_woK0_woSidn_K0nn = (TH2D*)fileK0nn->Get("nmom_cospip_woK0_woSidn");
  TH2D *nmom_cospim_woK0_woSidn_K0nn = (TH2D*)fileK0nn->Get("nmom_cospim_woK0_woSidn");
  TH2D *nmom_phinpip_woK0_woSidn_K0nn = (TH2D*)fileK0nn->Get("nmom_phinpip_woK0_woSidn");
  TH2D *nmom_phinpim_woK0_woSidn_K0nn = (TH2D*)fileK0nn->Get("nmom_phinpim_woK0_woSidn");
  TH2D *nmom_pipmom_woK0_woSidn_K0nn = (TH2D*)fileK0nn->Get("nmom_pipmom_woK0_woSidn");
  TH2D *nmom_pimmom_woK0_woSidn_K0nn = (TH2D*)fileK0nn->Get("nmom_pimmom_woK0_woSidn");
  TH2D *nmom_cosn_wK0_woSid_won_K0nn = (TH2D*)fileK0nn->Get("nmom_cosn_wK0_woSid_won");
  TH2D *nmom_cospip_wK0_woSid_won_K0nn = (TH2D*)fileK0nn->Get("nmom_cospip_wK0_woSid_won");
  TH2D *nmom_cospim_wK0_woSid_won_K0nn = (TH2D*)fileK0nn->Get("nmom_cospim_wK0_woSid_won");
  TH2D *nmom_phinpip_wK0_woSid_won_K0nn = (TH2D*)fileK0nn->Get("nmom_phinpip_wK0_woSid_won");
  TH2D *nmom_phinpim_wK0_woSid_won_K0nn = (TH2D*)fileK0nn->Get("nmom_phinpim_wK0_woSid_won");
  TH2D *nmom_pipmom_wK0_woSid_won_K0nn = (TH2D*)fileK0nn->Get("nmom_pipmom_wK0_woSid_won");
  TH2D *nmom_pimmom_wK0_woSid_won_K0nn = (TH2D*)fileK0nn->Get("nmom_pimmom_wK0_woSid_won");
  double Nnpip_K0nn = MMnmiss_IMnpip_woK0_K0nn->GetEntries();
  double Nnpim_K0nn = MMnmiss_IMnpim_woK0_K0nn->GetEntries();
  //real data
  TH2D *MMnmiss_IMnpip_woK0_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpip_dE_woK0");
  TH2D *MMnmiss_IMnpim_woK0_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpim_dE_woK0");
  TH2D *MMnmiss_IMnpip_woK0_woSm_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpip_dE_woK0_woSm");
  TH2D *MMnmiss_IMnpim_woK0_woSp_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpim_dE_woK0_woSp");
  TH2D *MMnmiss_IMnpip_woK0_woSm_n_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpip_dE_woK0_woSm_n");
  TH2D *MMnmiss_IMnpim_woK0_woSp_n_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpim_dE_woK0_woSp_n");
  TH2D *MMnmiss_IMnpip_woK0_woSidn_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpip_dE_woK0_woSidn");
  TH2D *MMnmiss_IMnpim_woK0_woSidn_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpim_dE_woK0_woSidn");
  TH2D *MMnmiss_IMnpip_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpip_dE_wK0_woSid_won");
  TH2D *MMnmiss_IMnpim_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpim_dE_wK0_woSid_won");
  TH2D *MMnmiss_IMpippim_rdata = (TH2D*)filerdata->Get("MMnmiss_IMpippim_dE");
  TH2D *MMnmiss_IMpippim_wSid_rdata = (TH2D*)filerdata->Get("MMnmiss_IMpippim_dE_wSid");
  TH2D *MMnmiss_IMpippim_woK0_woSidn_rdata = (TH2D*)filerdata->Get("MMnmiss_IMpippim_dE_woK0_woSidn");
  TH2D *MMnmiss_IMpippim_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("MMnmiss_IMpippim_dE_wK0_woSid_won");
  TH2D *IMnpim_IMnpip_rdata = (TH2D*)filerdata->Get("IMnpim_IMnpip_dE");
  TH2D *IMnpim_IMnpip_woK0_woSidn_rdata = (TH2D*)filerdata->Get("IMnpim_IMnpip_dE_woK0_woSidn");
  TH2D *q_IMnpipi_woK0_woSidn_rdata = (TH2D*)filerdata->Get("q_IMnpipi_woK0_woSidn");
  TH2D *q_IMnpipi_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("q_IMnpipi_wK0_woSid_won");
  TH2D *MMnmiss_Mompippim_woK0_woSidn_rdata  = (TH2D*)filerdata->Get("MMnmiss_Mompippim_dE_woK0_woSidn");
  TH2D *MMnmiss_Mompippim_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("MMnmiss_Mompippim_dE_wK0_woSid_won");
  TH2D *pipmom_MMnmiss_woK0_woSidn_rdata = (TH2D*)filerdata->Get("pipmom_MMnmiss_dE_woK0_woSidn");
  TH2D *pipmom_MMnmiss_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("pipmom_MMnmiss_dE_wK0_woSid_won");
  TH2D *pipmom_pimmom_woK0_woSidn_rdata = (TH2D*)filerdata->Get("pipmom_pimmom_dE_woK0_woSidn");
  TH2D *pipmom_pimmom_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("pipmom_pimmom_dE_wK0_woSid_won");
  TH2D *pimmom_MMnmiss_woK0_woSidn_rdata = (TH2D*)filerdata->Get("pimmom_MMnmiss_dE_woK0_woSidn");
  TH2D *pimmom_MMnmiss_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("pimmom_MMnmiss_dE_wK0_woSid_won");
  TH2D *MMom_MMass_woK0_woSidn_rdata = (TH2D*)filerdata->Get("MMom_MMass_woK0_woSidn");
  TH2D *MMom_MMass_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("MMom_MMass_wK0_woSid_won");
  TH2D *nmom_MMnmiss_woK0_woSidn_rdata = (TH2D*)filerdata->Get("nmom_MMnmiss_woK0_woSidn");
  TH2D *nmom_MMnmiss_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_MMnmiss_wK0_woSid_won");
  TH2D *q_IMnpipi_wSid_n_Sp_rdata = (TH2D*)filerdata->Get("q_IMnpipi_wSid_n_Sp");
  TH2D *q_IMnpipi_wSid_n_Sm_rdata = (TH2D*)filerdata->Get("q_IMnpipi_wSid_n_Sm");
  TH2D *nmom_cosn_woK0_woSidn_rdata = (TH2D*)filerdata->Get("nmom_cosn_woK0_woSidn");
  TH2D *nmom_cospip_woK0_woSidn_rdata = (TH2D*)filerdata->Get("nmom_cospip_woK0_woSidn");
  TH2D *nmom_cospim_woK0_woSidn_rdata = (TH2D*)filerdata->Get("nmom_cospim_woK0_woSidn");
  TH2D *nmom_phinpip_woK0_woSidn_rdata = (TH2D*)filerdata->Get("nmom_phinpip_woK0_woSidn");
  TH2D *nmom_phinpim_woK0_woSidn_rdata = (TH2D*)filerdata->Get("nmom_phinpim_woK0_woSidn");
  TH2D *nmom_pipmom_woK0_woSidn_rdata = (TH2D*)filerdata->Get("nmom_pipmom_woK0_woSidn");
  TH2D *nmom_pimmom_woK0_woSidn_rdata = (TH2D*)filerdata->Get("nmom_pimmom_woK0_woSidn");
  TH2D *nmom_cosn_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_cosn_wK0_woSid_won");
  TH2D *nmom_cospip_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_cospip_wK0_woSid_won");
  TH2D *nmom_cospim_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_cospim_wK0_woSid_won");
  TH2D *nmom_phinpip_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_phinpip_wK0_woSid_won");
  TH2D *nmom_phinpim_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_phinpim_wK0_woSid_won");
  TH2D *nmom_pipmom_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_pipmom_wK0_woSid_won");
  TH2D *nmom_pimmom_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_pimmom_wK0_woSid_won");
  double Nnpip_rdata = MMnmiss_IMnpip_woK0_rdata->GetEntries();
  double Nnpim_rdata = MMnmiss_IMnpim_woK0_rdata->GetEntries();

  //v2 2020113
  const double scaleFactor[6]= {1.0,
                                0.085*nrdata_Sp/nSigmap_Sp,
                                0.010*nrdata_Sp/nSigmam_Sp,
                                0.40*nrdata_Sp/nLambda_Sp,//
                                0.08*nrdata_Sp/nS0miss_Sp,//
                                0.015*nrdata_Sp/K0nn_Sp
                               };


  //arrays
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sp[6]= {
    MMnmiss_IMnpipi_woK0_wSid_Sp_rdata,
    MMnmiss_IMnpipi_woK0_wSid_Sp_Sigmap,
    MMnmiss_IMnpipi_woK0_wSid_Sp_Sigmam,
    MMnmiss_IMnpipi_woK0_wSid_Sp_Lambda,
    MMnmiss_IMnpipi_woK0_wSid_Sp_S0,
    MMnmiss_IMnpipi_woK0_wSid_Sp_K0nn
  };

  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sm[6]= {
    MMnmiss_IMnpipi_woK0_wSid_Sm_rdata,
    MMnmiss_IMnpipi_woK0_wSid_Sm_Sigmap,
    MMnmiss_IMnpipi_woK0_wSid_Sm_Sigmam,
    MMnmiss_IMnpipi_woK0_wSid_Sm_Lambda,
    MMnmiss_IMnpipi_woK0_wSid_Sm_S0,
    MMnmiss_IMnpipi_woK0_wSid_Sm_K0nn
  };

  TH2D* MMnmiss_IMnpip_woK0[6]= {
    MMnmiss_IMnpip_woK0_rdata,
    MMnmiss_IMnpip_woK0_Sigmap,
    MMnmiss_IMnpip_woK0_Sigmam,
    MMnmiss_IMnpip_woK0_Lambda,
    MMnmiss_IMnpip_woK0_S0,
    MMnmiss_IMnpip_woK0_K0nn
  };

  TH2D* MMnmiss_IMnpip_woK0_woSm[6]= {
    MMnmiss_IMnpip_woK0_woSm_rdata,
    MMnmiss_IMnpip_woK0_woSm_Sigmap,
    MMnmiss_IMnpip_woK0_woSm_Sigmam,
    MMnmiss_IMnpip_woK0_woSm_Lambda,
    MMnmiss_IMnpip_woK0_woSm_S0,
    MMnmiss_IMnpip_woK0_woSm_K0nn
  };

  TH2D* MMnmiss_IMnpip_woK0_woSm_n[6]= {
    MMnmiss_IMnpip_woK0_woSm_n_rdata,
    MMnmiss_IMnpip_woK0_woSm_n_Sigmap,
    MMnmiss_IMnpip_woK0_woSm_n_Sigmam,
    MMnmiss_IMnpip_woK0_woSm_n_Lambda,
    MMnmiss_IMnpip_woK0_woSm_n_S0,
    MMnmiss_IMnpip_woK0_woSm_n_K0nn
  };

  TH2D* MMnmiss_IMnpip_woK0_woSidn[6]= {
    MMnmiss_IMnpip_woK0_woSidn_rdata,
    MMnmiss_IMnpip_woK0_woSidn_Sigmap,
    MMnmiss_IMnpip_woK0_woSidn_Sigmam,
    MMnmiss_IMnpip_woK0_woSidn_Lambda,
    MMnmiss_IMnpip_woK0_woSidn_S0,
    MMnmiss_IMnpip_woK0_woSidn_K0nn
  };

  TH2D* MMnmiss_IMnpip_wK0_woSid_won[6]= {
    MMnmiss_IMnpip_wK0_woSid_won_rdata,
    MMnmiss_IMnpip_wK0_woSid_won_Sigmap,
    MMnmiss_IMnpip_wK0_woSid_won_Sigmam,
    MMnmiss_IMnpip_wK0_woSid_won_Lambda,
    MMnmiss_IMnpip_wK0_woSid_won_S0,
    MMnmiss_IMnpip_wK0_woSid_won_K0nn
  };


  TH2D* MMnmiss_IMnpim_woK0[6]= {
    MMnmiss_IMnpim_woK0_rdata,
    MMnmiss_IMnpim_woK0_Sigmap,
    MMnmiss_IMnpim_woK0_Sigmam,
    MMnmiss_IMnpim_woK0_Lambda,
    MMnmiss_IMnpim_woK0_S0,
    MMnmiss_IMnpim_woK0_K0nn
  };

  TH2D* MMnmiss_IMnpim_woK0_woSm[6]= {
    MMnmiss_IMnpim_woK0_woSp_rdata,
    MMnmiss_IMnpim_woK0_woSp_Sigmap,
    MMnmiss_IMnpim_woK0_woSp_Sigmam,
    MMnmiss_IMnpim_woK0_woSp_Lambda,
    MMnmiss_IMnpim_woK0_woSp_S0,
    MMnmiss_IMnpim_woK0_woSp_K0nn
  };

  TH2D* MMnmiss_IMnpim_woK0_woSp[6]= {
    MMnmiss_IMnpim_woK0_woSp_rdata,
    MMnmiss_IMnpim_woK0_woSp_Sigmap,
    MMnmiss_IMnpim_woK0_woSp_Sigmam,
    MMnmiss_IMnpim_woK0_woSp_Lambda,
    MMnmiss_IMnpim_woK0_woSp_S0,
    MMnmiss_IMnpim_woK0_woSp_K0nn
  };

  TH2D* MMnmiss_IMnpim_woK0_woSp_n[6]= {
    MMnmiss_IMnpim_woK0_woSp_n_rdata,
    MMnmiss_IMnpim_woK0_woSp_n_Sigmap,
    MMnmiss_IMnpim_woK0_woSp_n_Sigmam,
    MMnmiss_IMnpim_woK0_woSp_n_Lambda,
    MMnmiss_IMnpim_woK0_woSp_n_S0,
    MMnmiss_IMnpim_woK0_woSp_n_K0nn
  };

  TH2D* MMnmiss_IMnpim_woK0_woSidn[6]= {
    MMnmiss_IMnpim_woK0_woSidn_rdata,
    MMnmiss_IMnpim_woK0_woSidn_Sigmap,
    MMnmiss_IMnpim_woK0_woSidn_Sigmam,
    MMnmiss_IMnpim_woK0_woSidn_Lambda,
    MMnmiss_IMnpim_woK0_woSidn_S0,
    MMnmiss_IMnpim_woK0_woSidn_K0nn
  };

  TH2D* MMnmiss_IMnpim_wK0_woSid_won[6]= {
    MMnmiss_IMnpim_wK0_woSid_won_rdata,
    MMnmiss_IMnpim_wK0_woSid_won_Sigmap,
    MMnmiss_IMnpim_wK0_woSid_won_Sigmam,
    MMnmiss_IMnpim_wK0_woSid_won_Lambda,
    MMnmiss_IMnpim_wK0_woSid_won_S0,
    MMnmiss_IMnpim_wK0_woSid_won_K0nn
  };

  TH2D* IMnpim_IMnpip_woK0_woSidn[6]= {
    IMnpim_IMnpip_woK0_woSidn_rdata,
    IMnpim_IMnpip_woK0_woSidn_Sigmap,
    IMnpim_IMnpip_woK0_woSidn_Sigmam,
    IMnpim_IMnpip_woK0_woSidn_Lambda,
    IMnpim_IMnpip_woK0_woSidn_S0,
    IMnpim_IMnpip_woK0_woSidn_K0nn
  };

  TH2D* MMnmiss_IMpippim_wSid[6]= {
    MMnmiss_IMpippim_wSid_rdata,
    MMnmiss_IMpippim_wSid_Sigmap,
    MMnmiss_IMpippim_wSid_Sigmam,
    MMnmiss_IMpippim_wSid_Lambda,
    MMnmiss_IMpippim_wSid_S0,
    MMnmiss_IMpippim_wSid_K0nn
  };

  TH2D* MMnmiss_IMpippim_woK0_woSidn[6]= {
    MMnmiss_IMpippim_woK0_woSidn_rdata,
    MMnmiss_IMpippim_woK0_woSidn_Sigmap,
    MMnmiss_IMpippim_woK0_woSidn_Sigmam,
    MMnmiss_IMpippim_woK0_woSidn_Lambda,
    MMnmiss_IMpippim_woK0_woSidn_S0,
    MMnmiss_IMpippim_woK0_woSidn_K0nn
  };

  TH2D* MMnmiss_IMpippim_wK0_woSid_won[6]= {
    MMnmiss_IMpippim_wK0_woSid_won_rdata,
    MMnmiss_IMpippim_wK0_woSid_won_Sigmap,
    MMnmiss_IMpippim_wK0_woSid_won_Sigmam,
    MMnmiss_IMpippim_wK0_woSid_won_Lambda,
    MMnmiss_IMpippim_wK0_woSid_won_S0,
    MMnmiss_IMpippim_wK0_woSid_won_K0nn
  };

  TH2D* q_IMnpipi_woK0_woSidn[6]= {
    q_IMnpipi_woK0_woSidn_rdata,
    q_IMnpipi_woK0_woSidn_Sigmap,
    q_IMnpipi_woK0_woSidn_Sigmam,
    q_IMnpipi_woK0_woSidn_Lambda,
    q_IMnpipi_woK0_woSidn_S0,
    q_IMnpipi_woK0_woSidn_K0nn
  };

  TH2D* q_IMnpipi_wK0_woSid_won[6]= {
    q_IMnpipi_wK0_woSid_won_rdata,
    q_IMnpipi_wK0_woSid_won_Sigmap,
    q_IMnpipi_wK0_woSid_won_Sigmam,
    q_IMnpipi_wK0_woSid_won_Lambda,
    q_IMnpipi_wK0_woSid_won_S0,
    q_IMnpipi_wK0_woSid_won_K0nn
  };

  TH2D* MMnmiss_Mompippim_woK0_woSidn[6]= {
    MMnmiss_Mompippim_woK0_woSidn_rdata,
    MMnmiss_Mompippim_woK0_woSidn_Sigmap,
    MMnmiss_Mompippim_woK0_woSidn_Sigmam,
    MMnmiss_Mompippim_woK0_woSidn_Lambda,
    MMnmiss_Mompippim_woK0_woSidn_S0,
    MMnmiss_Mompippim_woK0_woSidn_K0nn
  };

  TH2D* MMnmiss_Mompippim_wK0_woSid_won[6]= {
    MMnmiss_Mompippim_wK0_woSid_won_rdata,
    MMnmiss_Mompippim_wK0_woSid_won_Sigmap,
    MMnmiss_Mompippim_wK0_woSid_won_Sigmam,
    MMnmiss_Mompippim_wK0_woSid_won_Lambda,
    MMnmiss_Mompippim_wK0_woSid_won_S0,
    MMnmiss_Mompippim_wK0_woSid_won_K0nn
  };

  TH2D* pipmom_MMnmiss_woK0_woSidn[6]= {
    pipmom_MMnmiss_woK0_woSidn_rdata,
    pipmom_MMnmiss_woK0_woSidn_Sigmap,
    pipmom_MMnmiss_woK0_woSidn_Sigmam,
    pipmom_MMnmiss_woK0_woSidn_Lambda,
    pipmom_MMnmiss_woK0_woSidn_S0,
    pipmom_MMnmiss_woK0_woSidn_K0nn
  };

  TH2D* pipmom_MMnmiss_wK0_woSid_won[6]= {
    pipmom_MMnmiss_wK0_woSid_won_rdata,
    pipmom_MMnmiss_wK0_woSid_won_Sigmap,
    pipmom_MMnmiss_wK0_woSid_won_Sigmam,
    pipmom_MMnmiss_wK0_woSid_won_Lambda,
    pipmom_MMnmiss_wK0_woSid_won_S0,
    pipmom_MMnmiss_wK0_woSid_won_K0nn
  };


  TH2D *pipmom_pimmom_woK0_woSidn[6]={
    pipmom_pimmom_woK0_woSidn_rdata,
    pipmom_pimmom_woK0_woSidn_Sigmap,
    pipmom_pimmom_woK0_woSidn_Sigmam,
    pipmom_pimmom_woK0_woSidn_Lambda,
    pipmom_pimmom_woK0_woSidn_S0,
    pipmom_pimmom_woK0_woSidn_K0nn
  };
  
  
  TH2D *pipmom_pimmom_wK0_woSid_won[6]={
    pipmom_pimmom_wK0_woSid_won_rdata,
    pipmom_pimmom_wK0_woSid_won_Sigmap,
    pipmom_pimmom_wK0_woSid_won_Sigmam,
    pipmom_pimmom_wK0_woSid_won_Lambda,
    pipmom_pimmom_wK0_woSid_won_S0,
    pipmom_pimmom_wK0_woSid_won_K0nn
  };



  TH2D* pimmom_MMnmiss_woK0_woSidn[6]= {
    pimmom_MMnmiss_woK0_woSidn_rdata,
    pimmom_MMnmiss_woK0_woSidn_Sigmap,
    pimmom_MMnmiss_woK0_woSidn_Sigmam,
    pimmom_MMnmiss_woK0_woSidn_Lambda,
    pimmom_MMnmiss_woK0_woSidn_S0,
    pimmom_MMnmiss_woK0_woSidn_K0nn
  };

  TH2D* pimmom_MMnmiss_wK0_woSid_won[6]= {
    pimmom_MMnmiss_wK0_woSid_won_rdata,
    pimmom_MMnmiss_wK0_woSid_won_Sigmap,
    pimmom_MMnmiss_wK0_woSid_won_Sigmam,
    pimmom_MMnmiss_wK0_woSid_won_Lambda,
    pimmom_MMnmiss_wK0_woSid_won_S0,
    pimmom_MMnmiss_wK0_woSid_won_K0nn
  };

  TH2D* MMom_MMass_woK0_woSidn[6]= {
    MMom_MMass_woK0_woSidn_rdata,
    MMom_MMass_woK0_woSidn_Sigmap,
    MMom_MMass_woK0_woSidn_Sigmam,
    MMom_MMass_woK0_woSidn_Lambda,
    MMom_MMass_woK0_woSidn_S0,
    MMom_MMass_woK0_woSidn_K0nn
  };

  TH2D* MMom_MMass_wK0_woSid_won[6]= {
    MMom_MMass_wK0_woSid_won_rdata,
    MMom_MMass_wK0_woSid_won_Sigmap,
    MMom_MMass_wK0_woSid_won_Sigmam,
    MMom_MMass_wK0_woSid_won_Lambda,
    MMom_MMass_wK0_woSid_won_S0,
    MMom_MMass_wK0_woSid_won_K0nn
  };


  TH2D* nmom_MMnmiss_woK0_woSidn[6]= {
    nmom_MMnmiss_woK0_woSidn_rdata,
    nmom_MMnmiss_woK0_woSidn_Sigmap,
    nmom_MMnmiss_woK0_woSidn_Sigmam,
    nmom_MMnmiss_woK0_woSidn_Lambda,
    nmom_MMnmiss_woK0_woSidn_S0,
    nmom_MMnmiss_woK0_woSidn_K0nn
  };

  TH2D* nmom_MMnmiss_wK0_woSid_won[6]= {
    nmom_MMnmiss_wK0_woSid_won_rdata,
    nmom_MMnmiss_wK0_woSid_won_Sigmap,
    nmom_MMnmiss_wK0_woSid_won_Sigmam,
    nmom_MMnmiss_wK0_woSid_won_Lambda,
    nmom_MMnmiss_wK0_woSid_won_S0,
    nmom_MMnmiss_wK0_woSid_won_K0nn
  };

  TH2D* q_IMnpipi_wSid_n_Sp[6]= {
    q_IMnpipi_wSid_n_Sp_rdata,
    q_IMnpipi_wSid_n_Sp_Sigmap,
    q_IMnpipi_wSid_n_Sp_Sigmam,
    q_IMnpipi_wSid_n_Sp_Lambda,
    q_IMnpipi_wSid_n_Sp_S0,
    q_IMnpipi_wSid_n_Sp_K0nn
  };

  TH2D* q_IMnpipi_wSid_n_Sm[6]= {
    q_IMnpipi_wSid_n_Sm_rdata,
    q_IMnpipi_wSid_n_Sm_Sigmap,
    q_IMnpipi_wSid_n_Sm_Sigmam,
    q_IMnpipi_wSid_n_Sm_Lambda,
    q_IMnpipi_wSid_n_Sm_S0,
    q_IMnpipi_wSid_n_Sm_K0nn
  };


  TH2D* nmom_cosn_woK0_woSidn[6]= {
    nmom_cosn_woK0_woSidn_rdata,
    nmom_cosn_woK0_woSidn_Sigmap,
    nmom_cosn_woK0_woSidn_Sigmam,
    nmom_cosn_woK0_woSidn_Lambda,
    nmom_cosn_woK0_woSidn_S0,
    nmom_cosn_woK0_woSidn_K0nn
  };

  TH2D* nmom_cospip_woK0_woSidn[6]= {
    nmom_cospip_woK0_woSidn_rdata,
    nmom_cospip_woK0_woSidn_Sigmap,
    nmom_cospip_woK0_woSidn_Sigmam,
    nmom_cospip_woK0_woSidn_Lambda,
    nmom_cospip_woK0_woSidn_S0,
    nmom_cospip_woK0_woSidn_K0nn
  };

  TH2D* nmom_cospim_woK0_woSidn[6]= {
    nmom_cospim_woK0_woSidn_rdata,
    nmom_cospim_woK0_woSidn_Sigmap,
    nmom_cospim_woK0_woSidn_Sigmam,
    nmom_cospim_woK0_woSidn_Lambda,
    nmom_cospim_woK0_woSidn_S0,
    nmom_cospim_woK0_woSidn_K0nn
  };

  TH2D* nmom_phinpip_woK0_woSidn[6]= {
    nmom_phinpip_woK0_woSidn_rdata,
    nmom_phinpip_woK0_woSidn_Sigmap,
    nmom_phinpip_woK0_woSidn_Sigmam,
    nmom_phinpip_woK0_woSidn_Lambda,
    nmom_phinpip_woK0_woSidn_S0,
    nmom_phinpip_woK0_woSidn_K0nn
  };

  TH2D* nmom_phinpim_woK0_woSidn[6]= {
    nmom_phinpim_woK0_woSidn_rdata,
    nmom_phinpim_woK0_woSidn_Sigmap,
    nmom_phinpim_woK0_woSidn_Sigmam,
    nmom_phinpim_woK0_woSidn_Lambda,
    nmom_phinpim_woK0_woSidn_S0,
    nmom_phinpim_woK0_woSidn_K0nn
  };

  TH2D* nmom_pipmom_woK0_woSidn[6]= {
    nmom_pipmom_woK0_woSidn_rdata,
    nmom_pipmom_woK0_woSidn_Sigmap,
    nmom_pipmom_woK0_woSidn_Sigmam,
    nmom_pipmom_woK0_woSidn_Lambda,
    nmom_pipmom_woK0_woSidn_S0,
    nmom_pipmom_woK0_woSidn_K0nn
  };

  TH2D* nmom_pimmom_woK0_woSidn[6]= {
    nmom_pimmom_woK0_woSidn_rdata,
    nmom_pimmom_woK0_woSidn_Sigmap,
    nmom_pimmom_woK0_woSidn_Sigmam,
    nmom_pimmom_woK0_woSidn_Lambda,
    nmom_pimmom_woK0_woSidn_S0,
    nmom_pimmom_woK0_woSidn_K0nn
  };

  TH2D* nmom_cosn_wK0_woSid_won[6]= {
    nmom_cosn_wK0_woSid_won_rdata,
    nmom_cosn_wK0_woSid_won_Sigmap,
    nmom_cosn_wK0_woSid_won_Sigmam,
    nmom_cosn_wK0_woSid_won_Lambda,
    nmom_cosn_wK0_woSid_won_S0,
    nmom_cosn_wK0_woSid_won_K0nn
  };

  TH2D* nmom_cospip_wK0_woSid_won[6]= {
    nmom_cospip_wK0_woSid_won_rdata,
    nmom_cospip_wK0_woSid_won_Sigmap,
    nmom_cospip_wK0_woSid_won_Sigmam,
    nmom_cospip_wK0_woSid_won_Lambda,
    nmom_cospip_wK0_woSid_won_S0,
    nmom_cospip_wK0_woSid_won_K0nn
  };

  TH2D* nmom_cospim_wK0_woSid_won[6]= {
    nmom_cospim_wK0_woSid_won_rdata,
    nmom_cospim_wK0_woSid_won_Sigmap,
    nmom_cospim_wK0_woSid_won_Sigmam,
    nmom_cospim_wK0_woSid_won_Lambda,
    nmom_cospim_wK0_woSid_won_S0,
    nmom_cospim_wK0_woSid_won_K0nn
  };

  TH2D* nmom_phinpip_wK0_woSid_won[6]= {
    nmom_phinpip_wK0_woSid_won_rdata,
    nmom_phinpip_wK0_woSid_won_Sigmap,
    nmom_phinpip_wK0_woSid_won_Sigmam,
    nmom_phinpip_wK0_woSid_won_Lambda,
    nmom_phinpip_wK0_woSid_won_S0,
    nmom_phinpip_wK0_woSid_won_K0nn
  };

  TH2D* nmom_phinpim_wK0_woSid_won[6]= {
    nmom_phinpim_wK0_woSid_won_rdata,
    nmom_phinpim_wK0_woSid_won_Sigmap,
    nmom_phinpim_wK0_woSid_won_Sigmam,
    nmom_phinpim_wK0_woSid_won_Lambda,
    nmom_phinpim_wK0_woSid_won_S0,
    nmom_phinpim_wK0_woSid_won_K0nn
  };

  TH2D* nmom_pipmom_wK0_woSid_won[6]= {
    nmom_pipmom_wK0_woSid_won_rdata,
    nmom_pipmom_wK0_woSid_won_Sigmap,
    nmom_pipmom_wK0_woSid_won_Sigmam,
    nmom_pipmom_wK0_woSid_won_Lambda,
    nmom_pipmom_wK0_woSid_won_S0,
    nmom_pipmom_wK0_woSid_won_K0nn
  };

  TH2D* nmom_pimmom_wK0_woSid_won[6]= {
    nmom_pimmom_wK0_woSid_won_rdata,
    nmom_pimmom_wK0_woSid_won_Sigmap,
    nmom_pimmom_wK0_woSid_won_Sigmam,
    nmom_pimmom_wK0_woSid_won_Lambda,
    nmom_pimmom_wK0_woSid_won_S0,
    nmom_pimmom_wK0_woSid_won_K0nn
  };


  for(int i=0; i<6; i++) {
    MMnmiss_IMnpipi_woK0_wSid_Sp[i]->Scale(scaleFactor[i]);
    MMnmiss_IMnpip_woK0[i]->Scale(scaleFactor[i]);
    MMnmiss_IMnpip_woK0_woSm[i]->Scale(scaleFactor[i]);
    MMnmiss_IMnpip_woK0_woSm_n[i]->Scale(scaleFactor[i]);
    MMnmiss_IMnpip_woK0_woSidn[i]->Scale(scaleFactor[i]);
    MMnmiss_IMnpip_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    MMnmiss_IMnpipi_woK0_wSid_Sm[i]->Scale(scaleFactor[i]);
    MMnmiss_IMnpim_woK0[i]->Scale(scaleFactor[i]);
    MMnmiss_IMnpim_woK0_woSp[i]->Scale(scaleFactor[i]);
    MMnmiss_IMnpim_woK0_woSp_n[i]->Scale(scaleFactor[i]);
    MMnmiss_IMnpim_woK0_woSidn[i]->Scale(scaleFactor[i]);
    MMnmiss_IMnpim_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    IMnpim_IMnpip_woK0_woSidn[i]->Scale(scaleFactor[i]);
    MMnmiss_IMpippim_wSid[i]->Scale(scaleFactor[i]);
    MMnmiss_IMpippim_woK0_woSidn[i]->Scale(scaleFactor[i]);
    MMnmiss_IMpippim_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    q_IMnpipi_woK0_woSidn[i]->Scale(scaleFactor[i]);
    q_IMnpipi_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    MMnmiss_Mompippim_woK0_woSidn[i]->Scale(scaleFactor[i]);
    MMnmiss_Mompippim_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    pipmom_MMnmiss_woK0_woSidn[i]->Scale(scaleFactor[i]);
    pipmom_MMnmiss_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    pipmom_pimmom_woK0_woSidn[i]->Scale(scaleFactor[i]);
    pipmom_pimmom_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    pimmom_MMnmiss_woK0_woSidn[i]->Scale(scaleFactor[i]);
    pimmom_MMnmiss_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    MMom_MMass_woK0_woSidn[i]->Scale(scaleFactor[i]);
    MMom_MMass_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    nmom_MMnmiss_woK0_woSidn[i]->Scale(scaleFactor[i]);
    nmom_MMnmiss_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    q_IMnpipi_wSid_n_Sp[i]->Scale(scaleFactor[i]);
    q_IMnpipi_wSid_n_Sm[i]->Scale(scaleFactor[i]);
    nmom_cosn_woK0_woSidn[i]->Scale(scaleFactor[i]);
    nmom_cospip_woK0_woSidn[i]->Scale(scaleFactor[i]);
    nmom_cospim_woK0_woSidn[i]->Scale(scaleFactor[i]);
    nmom_phinpip_woK0_woSidn[i]->Scale(scaleFactor[i]);
    nmom_phinpim_woK0_woSidn[i]->Scale(scaleFactor[i]);
    nmom_pipmom_woK0_woSidn[i]->Scale(scaleFactor[i]);
    nmom_pimmom_woK0_woSidn[i]->Scale(scaleFactor[i]);
    nmom_cosn_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    nmom_cospip_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    nmom_cospim_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    nmom_phinpip_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    nmom_phinpim_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    nmom_pipmom_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    nmom_pimmom_wK0_woSid_won[i]->Scale(scaleFactor[i]);
  }

  std::cout << "scaling done " << std::endl;

  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sp_mc = (TH2D*)MMnmiss_IMnpipi_woK0_wSid_Sp[1]->Clone("MMnmiss_IMnpipi_woK0_wSid_Sp_mc");
  for(int i=2; i<6; i++) MMnmiss_IMnpipi_woK0_wSid_Sp_mc->Add(MMnmiss_IMnpipi_woK0_wSid_Sp[i]);

  TCanvas *cMMnmiss_IMnpipi_woK0_wSid_Sp = new TCanvas("cMMnmiss_IMnpipi_woK0_wSid_Sp","cMMnmiss_IMnpipi_woK0_wSid_Sp",1200,800);
  cMMnmiss_IMnpipi_woK0_wSid_Sp->Divide(2,1);
  cMMnmiss_IMnpipi_woK0_wSid_Sp->cd(1);
  MMnmiss_IMnpipi_woK0_wSid_Sp[0]->SetTitle("#splitline{MMnmiss_IMnpipi_woK0_wSid_Sp}{  real Data}");
  MMnmiss_IMnpipi_woK0_wSid_Sp[0]->SetMinimum(1);
  MMnmiss_IMnpipi_woK0_wSid_Sp[0]->Draw("colz");
  cMMnmiss_IMnpipi_woK0_wSid_Sp->cd(2);
  MMnmiss_IMnpipi_woK0_wSid_Sp_mc->SetTitle("#splitline{MMnmiss_IMnpipi_woK0_wSid_Sp}{  MC sum}");
  MMnmiss_IMnpipi_woK0_wSid_Sp_mc->SetMinimum(1);
  MMnmiss_IMnpipi_woK0_wSid_Sp_mc->Draw("colz");

  TCanvas *cMMnmiss_woK0_wSid_Sp = new TCanvas("cMMnmiss_woK0_wSid_Sp","cMMnmiss_woK0_wSid_Sp",800,800);
  cMMnmiss_woK0_wSid_Sp->cd();
  TH1D* MMnmiss_woK0_wSid_Sp[6];
  for(int i=0; i<6; i++) MMnmiss_woK0_wSid_Sp[i] = (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sp[i]->ProjectionX(Form("MMnmiss_woK0_wSid_Sp_%s",name[i]));
  MMnmiss_woK0_wSid_Sp[0]->Draw("HE");
  TH1D* MMnmiss_woK0_wSid_Sp_mc = (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sp_mc->ProjectionX("MMnmiss_woK0_wSid_Sp_mc");
  MMnmiss_woK0_wSid_Sp_mc->SetLineColor(6);
  MMnmiss_woK0_wSid_Sp_mc->Draw("HEsame");
  for(int i=0; i<6; i++) {
    MMnmiss_woK0_wSid_Sp[i]->SetLineColor(colordef[i]);
    //MMnmiss_woK0_wSid_Sp[i]->Draw("HEsame");
  }

  TCanvas *cIMnpipi_woK0_wSid_Sp = new TCanvas("cIMnpipi_woK0_wSid_Sp","cIMnpipi_woK0_wSid_Sp");
  cIMnpipi_woK0_wSid_Sp->cd();
  TH1D* IMnpipi_woK0_wSid_Sp[6];
  for(int i=0; i<6; i++)IMnpipi_woK0_wSid_Sp[i] = (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sp[i]->ProjectionY(Form("IMnpipi_woK0_wSid_Sp_%s",name[i]));
  IMnpipi_woK0_wSid_Sp[0]->Draw("HE");
  TH1D* IMnpipi_woK0_wSid_Sp_mc = (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sp_mc->ProjectionY("IMnpipi_woK0_wSid_Sp_mc");
  IMnpipi_woK0_wSid_Sp_mc->SetLineColor(6);
  IMnpipi_woK0_wSid_Sp_mc->Draw("HEsame");
  for(int i=1; i<6; i++) {
    IMnpipi_woK0_wSid_Sp[i]->SetLineColor(colordef[i]);
    //IMnpipi_woK0_wSid_Sp[i]->Draw("HEsame");
  }


  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sm_mc = (TH2D*)MMnmiss_IMnpipi_woK0_wSid_Sm[1]->Clone("MMnmiss_IMnpipi_woK0_wSid_Sm_mc");
  for(int i=2; i<6; i++) MMnmiss_IMnpipi_woK0_wSid_Sm_mc->Add(MMnmiss_IMnpipi_woK0_wSid_Sm[i]);

  TCanvas *cMMnmiss_IMnpipi_woK0_wSid_Sm = new TCanvas("cMMnmiss_IMnpipi_woK0_wSid_Sm","cMMnmiss_IMnpipi_woK0_wSid_Sm",1200,800);
  cMMnmiss_IMnpipi_woK0_wSid_Sm->Divide(2,1);
  cMMnmiss_IMnpipi_woK0_wSid_Sm->cd(1);
  MMnmiss_IMnpipi_woK0_wSid_Sm[0]->SetTitle("#splitline{MMnmiss_IMnpipi_woK0_wSid_Sm}{  real data}");
  MMnmiss_IMnpipi_woK0_wSid_Sm[0]->SetMinimum(1);
  MMnmiss_IMnpipi_woK0_wSid_Sm[0]->Draw("colz");
  cMMnmiss_IMnpipi_woK0_wSid_Sm->cd(2);
  MMnmiss_IMnpipi_woK0_wSid_Sm_mc->SetTitle("#splitline{MMnmiss_IMnpipi_woK0_wSid_Sm}{  MC sum}");
  MMnmiss_IMnpipi_woK0_wSid_Sm_mc->SetMinimum(1);
  MMnmiss_IMnpipi_woK0_wSid_Sm_mc->Draw("colz");

  TCanvas *cMMnmiss_woK0_wSid_Sm = new TCanvas("cMMnmiss_woK0_wSid_Sm","cMMnmiss_woK0_wSid_Sm",800,800);
  cMMnmiss_woK0_wSid_Sm->cd();
  TH1D* MMnmiss_woK0_wSid_Sm[6];
  for(int i=0; i<6; i++) MMnmiss_woK0_wSid_Sm[i] = (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sm[i]->ProjectionX(Form("MMnmiss_woK0_wSid_Sm_%s",name[i]));
  MMnmiss_woK0_wSid_Sm[0]->Draw("HE");
  TH1D* MMnmiss_woK0_wSid_Sm_mc = (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sm_mc->ProjectionX("MMnmiss_woK0_wSid_Sm_mc");
  MMnmiss_woK0_wSid_Sm_mc->SetLineColor(6);
  MMnmiss_woK0_wSid_Sm_mc->Draw("HEsame");
  for(int i=1; i<6; i++) {
    MMnmiss_woK0_wSid_Sm[i]->SetLineColor(colordef[i]);
    //MMnmiss_woK0_wSid_Sm[i]->Draw("HEsame");
  }

  TCanvas *cIMnpipi_woK0_wSid_Sm  = new TCanvas("cIMnpipi_woK0_wSid_Sm","cIMnpipi_woK0_wSid_Sm");
  cIMnpipi_woK0_wSid_Sm->cd();
  TH1D* IMnpipi_woK0_wSid_Sm[6];
  for(int i=0; i<6; i++)IMnpipi_woK0_wSid_Sm[i] = (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sm[i]->ProjectionY(Form("IMnpipi_woK0_wSid_Sm_%s",name[i]));
  IMnpipi_woK0_wSid_Sm[0]->Draw("HE");
  TH1D* IMnpipi_woK0_wSid_Sm_mc = (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sm_mc->ProjectionY("IMnpipi_woK0_wSid_Sm_mc");
  IMnpipi_woK0_wSid_Sm_mc->SetLineColor(6);
  IMnpipi_woK0_wSid_Sm_mc->Draw("HEsame");
  for(int i=1; i<6; i++) {
    IMnpipi_woK0_wSid_Sm[i]->SetLineColor(colordef[i]);
    //IMnpipi_woK0_wSid_Sm[i]->Draw("HEsame");
  }


  //missing mass and IMnpip/npim w/ Sp or Sm selection w/o K0
  //
  //ploting Sp mode
  //

  //w/o K0, including Sp/Sm mode
  TH2D *MMnmiss_IMnpip_woK0_mc = (TH2D*)MMnmiss_IMnpip_woK0[1]->Clone("MMnmiss_IMnpip_woK0_mc");
  //adding all MC data
  for(int i=2; i<6; i++)  MMnmiss_IMnpip_woK0_mc->Add(MMnmiss_IMnpip_woK0[i]);

  TCanvas *cMMnmiss_IMnpip_woK0 = new TCanvas("cMMnmiss_IMnpip_woK0","cMMnmiss_IMnpip_woK0",1200,800);
  cMMnmiss_IMnpip_woK0->Divide(2,1);
  cMMnmiss_IMnpip_woK0->cd(1);
  MMnmiss_IMnpip_woK0_rdata->SetTitle("#splitline{MMnmiss_IMnpip_woK0}{  Real data}");
  MMnmiss_IMnpip_woK0_rdata->SetMinimum(1);
  MMnmiss_IMnpip_woK0_rdata->Draw("colz");
  cMMnmiss_IMnpip_woK0->cd(2);
  MMnmiss_IMnpip_woK0_mc->SetTitle("#splitline{MMnmiss_IMnpip_woK0}{  MC sum}");
  MMnmiss_IMnpip_woK0_mc->SetMinimum(1);
  MMnmiss_IMnpip_woK0_mc->Draw("colz");

  //w/o K0, w/o Sm mode
  TH2D *MMnmiss_IMnpip_woK0_woSm_mc = (TH2D*)MMnmiss_IMnpip_woK0_woSm[1]->Clone("MMnmiss_IMnpip_woK0_woSm_mc");
  //adding all MC data
  for(int i=2; i<6; i++)MMnmiss_IMnpip_woK0_woSm_mc->Add(MMnmiss_IMnpip_woK0_woSm[i]);

  //w/o K0, w/o Sm mode, selecting missing neutron
  TH2D *MMnmiss_IMnpip_woK0_woSm_n_mc = (TH2D*)MMnmiss_IMnpip_woK0_woSm_n[1]->Clone("MMnmiss_IMnpip_woK0_woSm_n_mc");
  for(int i=2; i<6; i++)MMnmiss_IMnpip_woK0_woSm_n_mc->Add(MMnmiss_IMnpip_woK0_woSm_n[i]);

  TCanvas *cMMnmiss_IMnpip_woK0_woSm = new TCanvas("cMMnmiss_IMnpip_woK0_woSm","cMMnmiss_IMnpip_woK0_woSm",1200,800);
  cMMnmiss_IMnpip_woK0_woSm->Divide(2,1);
  cMMnmiss_IMnpip_woK0_woSm->cd(1);
  MMnmiss_IMnpip_woK0_woSm_rdata->SetTitle("#splitline{MMnmiss_IMnpip_woK0_woSm}{  Real data}");
  MMnmiss_IMnpip_woK0_woSm_rdata->SetMinimum(1);
  MMnmiss_IMnpip_woK0_woSm_rdata->Draw("colz");
  cMMnmiss_IMnpip_woK0_woSm->cd(2);
  MMnmiss_IMnpip_woK0_woSm_mc->SetTitle("#splitline{MMnmiss_IMnpip_woK0_woSm}{  MC sum}");
  MMnmiss_IMnpip_woK0_woSm_mc->SetMinimum(1);
  MMnmiss_IMnpip_woK0_woSm_mc->Draw("colz");

  TCanvas *cIMnpip_woK0_woSm_n = new TCanvas("cIMnpip_woK0_woSm_n","cIMnpip_woK0_woSm_n");
  cIMnpip_woK0_woSm_n->cd();
  TH1D *IMnpip_woK0_woSm_n[6];
  for(int i=0; i<6; i++) {
    IMnpip_woK0_woSm_n[i] = (TH1D*)MMnmiss_IMnpip_woK0_woSm_n[i]->ProjectionX(Form("IMnpip_woK0_woSm_n_%s",name[i]));
    IMnpip_woK0_woSm_n[i]->SetLineColor(colordef[i]);
  }
  TH1D* IMnpip_woK0_woSm_n_mc = (TH1D*)MMnmiss_IMnpip_woK0_woSm_n_mc->ProjectionX("IMnpip_woK0_woSm_n_mc");
  IMnpip_woK0_woSm_n[0]->Draw("HE");
  //for(int i=1;i<6;i++)IMnpip_woK0_woSm_n[i]->Draw("HEsame");
  IMnpip_woK0_woSm_n_mc->SetLineColor(6);
  IMnpip_woK0_woSm_n_mc->Draw("same");

  //w/o K0, w/o ((Sp or Sm) & missing neutron)
  TH2D *MMnmiss_IMnpip_woK0_woSidn_mc = (TH2D*)MMnmiss_IMnpip_woK0_woSidn[1]->Clone("MMnmiss_IMnpip_woK0_woSidn_mc");
  for(int i=2; i<6; i++) MMnmiss_IMnpip_woK0_woSidn_mc->Add(MMnmiss_IMnpip_woK0_woSidn[i]);

  TCanvas *cMMnmiss_IMnpip_woSidn = new TCanvas("cMMnmiss_IMnpip_woSidn","cMMnmiss_IMnpip_woSidn",1200,800);
  cMMnmiss_IMnpip_woSidn->Divide(2,1);
  cMMnmiss_IMnpip_woSidn->cd(1);
  MMnmiss_IMnpip_woK0_woSidn_rdata->SetTitle("#splitline{MMnmiss_IMnpip_woK0_woSidn}{  Real data}");
  MMnmiss_IMnpip_woK0_woSidn_rdata->SetMinimum(1);
  MMnmiss_IMnpip_woK0_woSidn_rdata->Draw("colz");
  cMMnmiss_IMnpip_woSidn->cd(2);
  MMnmiss_IMnpip_woK0_woSidn_mc->SetTitle("#splitline{MMnmiss_IMnpip_woK0_woSidn}{  MC sum}");
  MMnmiss_IMnpip_woK0_woSidn_mc->SetMinimum(1);
  MMnmiss_IMnpip_woK0_woSidn_mc->Draw("colz");

  //projection to Missing mass (miss n & Sigma+/-)
  TCanvas *cMMnmiss_woK0_woSidn = new TCanvas("cMMnmiss_woK0_woSidn","cMMnmiss_woK0_woSidn");
  cMMnmiss_woK0_woSidn->cd();
  TH1D* MMnmiss_woK0_woSidn[6];
  for(int i=0; i<6; i++) MMnmiss_woK0_woSidn[i] = (TH1D*)MMnmiss_IMnpip_woK0_woSidn[i]->ProjectionY(Form("MMnmiss_woK0_woSidn_%s",name[i]));
  MMnmiss_woK0_woSidn[0]->Draw("HE");
  TH1D* MMnmiss_woK0_woSidn_mc = (TH1D*)MMnmiss_IMnpip_woK0_woSidn_mc->ProjectionY("MMnmiss_woK0_woSidn_mc");
  MMnmiss_woK0_woSidn_mc->SetLineColor(6);
  MMnmiss_woK0_woSidn_mc->Draw("HEsame");

  for(int i=1; i<6; i++) {
    MMnmiss_woK0_woSidn[i]->SetLineColor(colordef[i]);
    //MMnmiss_woK0_woSidn[i]->Draw("HEsame");
  }

  //Data/MC before modifying MC data
  TCanvas *cMMnmiss_woK0_woSidn_ratio = new TCanvas("cMMnmiss_woK0_woSidn_ratio","cMMnmiss_woK0_woSid_ratio");
  TH1D* MMnmiss_woK0_woSidn_ratio = (TH1D*)MMnmiss_woK0_woSidn[0]->Clone("MMnmiss_woK0_woSidn_ratio");
  MMnmiss_woK0_woSidn_ratio->Divide(MMnmiss_woK0_woSidn_mc);
  MMnmiss_woK0_woSidn_ratio->SetTitle("Data/MC");
  MMnmiss_woK0_woSidn_ratio->GetYaxis()->SetRangeUser(0,3);
  MMnmiss_woK0_woSidn_ratio->Draw("HE");


  TF1 *sgf_MMnmiss = new TF1("sgf_MMnmiss","[0]*exp(-0.5*pow((x-[1])/([2]+(x<[1])*[3]*(x-[1])),2))");
  TF1 *fgaus_MMnmiss_high = new TF1("fgaus_MMnmiss_high","gaus");
  sgf_MMnmiss->SetParameter(0,1.82171e+00);
  sgf_MMnmiss->SetParameter(1,8.56016e-01);
  sgf_MMnmiss->SetParameter(2,6.81677e-01);
  sgf_MMnmiss->SetLineColor(2);

  MMnmiss_woK0_woSidn_ratio->Fit("sgf_MMnmiss","R","",0.,1.116);
  MMnmiss_woK0_woSidn_ratio->Fit("fgaus_MMnmiss_high","R+","",1.116,1.5);
  Double_t param_MMnmiss[7];
  sgf_MMnmiss->GetParameters(&param_MMnmiss[0]);
  fgaus_MMnmiss_high->GetParameters(&param_MMnmiss[4]);
  TF1 *evalf_MMnmiss = new TF1("evalf_MMnmiss",func_MMnmiss,0,1.5,7);
  evalf_MMnmiss->SetParameters(param_MMnmiss);
  evalf_MMnmiss->SetLineColor(4);
  evalf_MMnmiss->Draw("same");
  evalf_MMnmiss->Print();


  TH2D *MMom_MMass_woK0_woSidn_mc = (TH2D*)MMom_MMass_woK0_woSidn[1]->Clone("MMom_MMass_woK0_woSidn_mc");
  for(int i=2; i<6; i++) MMom_MMass_woK0_woSidn_mc->Add(MMom_MMass_woK0_woSidn[i]);
  TCanvas *cMMom_MMass_woK0_woSidn = new TCanvas("cMMom_MMass_woK0_woSidn","cMMom_MMass_woK0_woSidn");
  cMMom_MMass_woK0_woSidn->Divide(2,1);
  cMMom_MMass_woK0_woSidn->cd(1);
  MMom_MMass_woK0_woSidn[0]->SetMinimum(1);
  MMom_MMass_woK0_woSidn[0]->Draw("colz");
  cMMom_MMass_woK0_woSidn->cd(2);
  MMom_MMass_woK0_woSidn_mc->SetMinimum(1);
  MMom_MMass_woK0_woSidn_mc->Draw("colz");

  TH1D* MMom_woK0_woSidn_mc = MMom_MMass_woK0_woSidn_mc->ProjectionY("MMom_woK0_woSidn_mc");
  TH1D* MMom_woK0_woSidn[6];
  for(int i=0; i<6; i++) {
    MMom_woK0_woSidn[i] = (TH1D*)MMom_MMass_woK0_woSidn[i]->ProjectionY(Form("MMom_woK0_woSidn_%s",name[i]));
    MMom_woK0_woSidn[i]->SetLineColor(colordef[i]);
  }
  TCanvas *cMMom_woK0_woSidn = new TCanvas("cMMom_woK0_woSidn","cMMom_woK0_woSidn");
  MMom_woK0_woSidn[0]->Draw("HE");
  MMom_woK0_woSidn_mc->SetLineColor(6);
  MMom_woK0_woSidn_mc->Draw("HEsame");
  //for(int i=1;i<6;i++)MMom_woK0_woSidn[i]->Draw("HEsame");


  TCanvas *cMMom_woK0_woSidn_ratio = new TCanvas("cMMom_woK0_woSidn_ratio","cMMom_woK0_woSidn_ratio");
  cMMom_woK0_woSidn_ratio->cd();
  TH1D* MMom_woK0_woSidn_ratio = (TH1D*)MMom_woK0_woSidn[0]->Clone("MMom_woK0_woSidn_ratio");
  MMom_woK0_woSidn_ratio->Divide(MMom_woK0_woSidn_mc);
  MMom_woK0_woSidn_ratio->GetYaxis()->SetRangeUser(0,5);
  MMom_woK0_woSidn_ratio->SetTitle("Data/MC");
  MMom_woK0_woSidn_ratio->Draw("HE");

  Double_t param_MMom[5];
  TF1 *evalf_MMom = new TF1("evalf_MMom","pol5",0.4,1.5);
  evalf_MMom->SetLineColor(4);
  evalf_MMom->SetTitle("MMom");
  MMom_woK0_woSidn_ratio->Fit("evalf_MMom","","",0.4,1.5);


  //projection to IMnpip (miss n & Sigma+/-)
  TCanvas *cIMnpip_woK0_woSidn = new TCanvas("cIMnpip_woK0_woSidn","cIMnpip_woK0_woSidn");
  cIMnpip_woK0_woSidn->cd();
  TH1D* IMnpip_woK0_woSidn[6];
  for(int i=0; i<6; i++)IMnpip_woK0_woSidn[i] = (TH1D*)MMnmiss_IMnpip_woK0_woSidn[i]->ProjectionX(Form("IMnpip_woK0_woSidn_%s",name[i]));
  IMnpip_woK0_woSidn[0]->Draw("HE");//rdata
  TH1D* IMnpip_woK0_woSidn_mc = (TH1D*)MMnmiss_IMnpip_woK0_woSidn_mc->ProjectionX("IMnpip_woK0_woSidn_mc");
  IMnpip_woK0_woSidn_mc->SetLineColor(6);
  IMnpip_woK0_woSidn_mc->Draw("HEsame");
  for(int i=1; i<6; i++) {
    IMnpip_woK0_woSidn[i]->SetLineColor(colordef[i]);
    //IMnpip_woK0_woSidn[i]->Draw("HEsame");
  }

  TCanvas *cIMnpip_woK0_woSidn_ratio = new TCanvas("cIMnpip_woK0_woSidn_ratio","cIMnpip_woK0_woSidn_ratio");
  cIMnpip_woK0_woSidn_ratio->cd();
  TH1D* IMnpip_woK0_woSidn_ratio = (TH1D*)IMnpip_woK0_woSidn[0]->Clone("IMnpip_woK0_woSidn_ratio");
  IMnpip_woK0_woSidn_ratio->Divide(IMnpip_woK0_woSidn_mc);
  IMnpip_woK0_woSidn_ratio->SetTitle("IMnpip_woK0_woSidn Data/MC");
  IMnpip_woK0_woSidn_ratio->GetYaxis()->SetRangeUser(0,3);
  IMnpip_woK0_woSidn_ratio->Draw("HE");

  TF1* fgaus_IMnpip_1 = new TF1("fgaus_IMnpip_1","gaus",1.06,1.10);
  fgaus_IMnpip_1->SetParameters(0,1.636);
  fgaus_IMnpip_1->SetParameters(1,1.102);
  fgaus_IMnpip_1->SetParameters(2,0.02845);
  IMnpip_woK0_woSidn_ratio->Fit("fgaus_IMnpip_1","R","",1.08,1.10);

  TF1* fexpo_IMnpip_2 = new TF1("fexpo_IMnpip_2","expo",1.10,1.25);
  fexpo_IMnpip_2->SetParameters(0,1.667);
  fexpo_IMnpip_2->SetParameters(1,-1.117);
  IMnpip_woK0_woSidn_ratio->Fit("fexpo_IMnpip_2","R+","",1.10,1.25);

  TF1* fgaus_IMnpip_3 = new TF1("fgaus_IMnpip_3","gaus",1.25,2.0);
  fgaus_IMnpip_3->SetParameters(0,10.83);
  fgaus_IMnpip_3->SetParameters(1,5.094);
  fgaus_IMnpip_3->SetParameters(2,1.87);
  IMnpip_woK0_woSidn_ratio->Fit("fgaus_IMnpip_3","R+","",1.25,2.0);

  Double_t param_IMnpip[8];
  fgaus_IMnpip_1->GetParameters(&param_IMnpip[0]);
  fexpo_IMnpip_2->GetParameters(&param_IMnpip[3]);
  fgaus_IMnpip_3->GetParameters(&param_IMnpip[5]);

  TF1 *evalf_IMnpip = new TF1("evalf_IMnpip",func_IMnpip,1.06,2.0,8);
  evalf_IMnpip->SetParameters(param_IMnpip);
  evalf_IMnpip->SetLineColor(4);
  evalf_IMnpip->Draw("same");

  //
  //Sigma-
  TH2D *MMnmiss_IMnpim_woK0_mc = (TH2D*)MMnmiss_IMnpim_woK0[1]->Clone("MMnmiss_IMnpim_woK0_mc");
  //adding all MC data
  for(int i=2; i<6; i++)MMnmiss_IMnpim_woK0_mc->Add(MMnmiss_IMnpim_woK0[i]);

  TCanvas *cMMnmiss_IMnpim_woK0 = new TCanvas("cMMnmiss_IMnpim_woK0","cMMnmiss_IMnpim_woK0",1200,800);
  cMMnmiss_IMnpim_woK0->Divide(2,1);
  cMMnmiss_IMnpim_woK0->cd(1);
  MMnmiss_IMnpim_woK0_rdata->SetTitle("#splitline{MMnmiss_IMnpim_woK0}{  Real data}");
  MMnmiss_IMnpim_woK0_rdata->SetMinimum(1);
  MMnmiss_IMnpim_woK0_rdata->Draw("colz");
  cMMnmiss_IMnpim_woK0->cd(2);
  MMnmiss_IMnpim_woK0_mc->SetTitle("#splitline{MMnmiss_IMnpim_woK0}{ MC sum}");
  MMnmiss_IMnpim_woK0_mc->SetMinimum(1);
  MMnmiss_IMnpim_woK0_mc->Draw("colz");

  //w/o K0, w/o Sp mode
  TH2D *MMnmiss_IMnpim_woK0_woSp_mc = (TH2D*)MMnmiss_IMnpim_woK0_woSp[1]->Clone("MMnmiss_IMnpim_woK0_woSp_mc");
  //adding all MC data
  for(int i=2; i<6; i++) {
    //MMnmiss_IMnpim_woK0_woSp[i]->Print("base");
    MMnmiss_IMnpim_woK0_woSp_mc->Add(MMnmiss_IMnpim_woK0_woSp[i]);
  }
  //w/o K0, w/o Sp mode, selecting missing neutron
  TH2D *MMnmiss_IMnpim_woK0_woSp_n_mc = (TH2D*)MMnmiss_IMnpim_woK0_woSp_n[1]->Clone("MMnmiss_IMnpim_woK0_woSp_n_mc");
  for(int i=2; i<6; i++)MMnmiss_IMnpim_woK0_woSp_n_mc->Add(MMnmiss_IMnpim_woK0_woSp_n[i]);

  TCanvas *cMMnmiss_IMnpim_woK0_woSp = new TCanvas("cMMnmiss_IMnpim_woK0_woSp","cMMnmiss_IMnpim_woK0_woSp",1200,800);
  cMMnmiss_IMnpim_woK0_woSp->Divide(2,1);
  cMMnmiss_IMnpim_woK0_woSp->cd(1);
  MMnmiss_IMnpim_woK0_woSp_rdata->SetTitle("#splitline{MMnmiss_IMnpim_woK0_woSp}{  Real data}");
  MMnmiss_IMnpim_woK0_woSp_rdata->SetMinimum(1);
  MMnmiss_IMnpim_woK0_woSp_rdata->Draw("colz");
  cMMnmiss_IMnpim_woK0_woSp->cd(2);
  MMnmiss_IMnpim_woK0_woSp_mc->SetTitle("#splitline{MMnmiss_IMnpim_woK0_woSp}{  MC sum}");
  MMnmiss_IMnpim_woK0_woSp_mc->SetMinimum(1);
  MMnmiss_IMnpim_woK0_woSp_mc->Draw("colz");

  TCanvas *cIMnpim_woK0_woSp_n = new TCanvas("cIMnpim_woK0_woSp_n","cIMnpim_woK0_woSp_n");
  cIMnpim_woK0_woSp_n->cd();
  TH1D *IMnpim_woK0_woSp_n[6];
  for(int i=0; i<6; i++) {
    IMnpim_woK0_woSp_n[i] = (TH1D*)MMnmiss_IMnpim_woK0_woSp_n[i]->ProjectionX(Form("IMnpim_woK0_woSp_n_%s",name[i]));
    IMnpim_woK0_woSp_n[i]->SetLineColor(colordef[i]);
  }
  TH1D* IMnpim_woK0_woSp_n_mc = (TH1D*) MMnmiss_IMnpim_woK0_woSp_n_mc->ProjectionX("IMnpim_woK0_woSp_n_mc");
  IMnpim_woK0_woSp_n[0]->Draw("HE");
  //for(int i=1;i<6;i++)IMnpim_woK0_woSp_n[i]->Draw("HEsame");
  IMnpim_woK0_woSp_n_mc->SetLineColor(6);
  IMnpim_woK0_woSp_n_mc->Draw("same");


  //w/o K0, w/o ((Sp or Sm) & missing neutron)
  TH2D *MMnmiss_IMnpim_woK0_woSidn_mc = (TH2D*)MMnmiss_IMnpim_woK0_woSidn[1]->Clone("MMnmiss_IMnpim_woK0_woSidn_mc");
  for(int i=2; i<6; i++)MMnmiss_IMnpim_woK0_woSidn_mc->Add(MMnmiss_IMnpim_woK0_woSidn[i]);

  TCanvas *cMMnmiss_IMnpim_woK0_woSidn = new TCanvas("cMMnmiss_IMnpim_woK0_woSidn","cMMnmiss_IMnpim_woK0_woSidn",1200,800);
  cMMnmiss_IMnpim_woK0_woSidn->Divide(2,1);
  cMMnmiss_IMnpim_woK0_woSidn->cd(1);
  MMnmiss_IMnpim_woK0_woSidn_rdata->SetTitle("#splitline{MMnmiss_IMnpim_woK0_woSidn}{  Real data}");
  MMnmiss_IMnpim_woK0_woSidn_rdata->SetMinimum(1);
  MMnmiss_IMnpim_woK0_woSidn_rdata->Draw("colz");
  cMMnmiss_IMnpim_woK0_woSidn->cd(2);
  MMnmiss_IMnpim_woK0_woSidn_mc->SetTitle("#splitline{MMnmiss_IMnpim_woK0_woSidn}{  MC sum}");
  MMnmiss_IMnpim_woK0_woSidn_mc->SetMinimum(1);
  MMnmiss_IMnpim_woK0_woSidn_mc->Draw("colz");

  TCanvas *cIMnpim_woK0_woSidn = new TCanvas("IMnpim_woK0_woSidn","IMnpim_woK0_woSidn");
  cIMnpim_woK0_woSidn->cd();
  TH1D* IMnpim_woK0_woSidn[6];
  for(int i=0; i<6; i++)IMnpim_woK0_woSidn[i] = (TH1D*)MMnmiss_IMnpim_woK0_woSidn[i]->ProjectionX(Form("IMnpim_woK0_woSidn_%s",name[i]));
  TH1D* IMnpim_woK0_woSidn_mc = MMnmiss_IMnpim_woK0_woSidn_mc->ProjectionX("IMnpim_woK0_woSidn_mc");
  IMnpim_woK0_woSidn[0]->Draw("HE");
  IMnpim_woK0_woSidn_mc->SetLineColor(6);
  IMnpim_woK0_woSidn_mc->Draw("HEsame");
  for(int i=1; i<6; i++) {
    IMnpim_woK0_woSidn[i]->SetLineColor(colordef[i]);
    //IMnpim_woK0_woSidn[i]->Draw("HEsame");
  }

  TCanvas *cIMnpim_woK0_woSidn_ratio = new TCanvas("cIMnpim_woK0_woSidn_ratio","cIMnpim_woK0_woSidn_ratio");
  cIMnpim_woK0_woSidn_ratio->cd();
  TH1D* IMnpim_woK0_woSidn_ratio = (TH1D*)IMnpim_woK0_woSidn[0]->Clone("IMnpim_woK0_woSidn_ratio");
  IMnpim_woK0_woSidn_ratio->Divide(IMnpim_woK0_woSidn_mc);
  IMnpim_woK0_woSidn_ratio->SetTitle("IMnpim_woK0_woSidn Data/MC");
  IMnpim_woK0_woSidn_ratio->GetYaxis()->SetRangeUser(0,3);
  IMnpim_woK0_woSidn_ratio->Draw("HE");

  TF1* pol3_IMnpim_1 = new TF1("pol3_IMnpim_1","pol3",1.00,1.10);
  pol3_IMnpim_1->SetParameter(0,-32074.9);
  pol3_IMnpim_1->SetParameter(1,85205.3);
  pol3_IMnpim_1->SetParameter(2,-75374.);
  pol3_IMnpim_1->SetParameter(3,22204.);
  IMnpim_woK0_woSidn_ratio->Fit("pol3_IMnpim_1","R","",1.07,1.10);

  TF1* pol3_IMnpim_2 = new TF1("pol3_IMnpim_2","pol3",1.10,2.00);
  pol3_IMnpim_2->SetParameter(0,38.13);
  pol3_IMnpim_2->SetParameter(1,-62.3139);
  pol3_IMnpim_2->SetParameter(2,32.172);
  pol3_IMnpim_2->SetParameter(3,-4.81);
  IMnpim_woK0_woSidn_ratio->Fit("pol3_IMnpim_2","R+","",1.10,2.00);

  Double_t param_IMnpim[8];
  pol3_IMnpim_1->GetParameters(&param_IMnpim[0]);
  pol3_IMnpim_2->GetParameters(&param_IMnpim[4]);
  TF1 *evalf_IMnpim = new TF1("evalf_IMnpim",func_IMnpim,1.00,2.00,8);
  evalf_IMnpim->SetParameters(param_IMnpim);
  evalf_IMnpim->SetTitle("IMnpim");
  evalf_IMnpim->SetLineColor(4);
  evalf_IMnpim->Draw("same");


  //MMnmiss vs IMpippim wSid
  TH2D* MMnmiss_IMpippim_wSid_mc = (TH2D*)MMnmiss_IMpippim_wSid[1]->Clone("MMnmiss_IMpippim_wSid_mc");
  for(int i=1; i<6; i++) MMnmiss_IMpippim_wSid_mc->Add(MMnmiss_IMpippim_wSid[i]);

  TCanvas *cMMnmiss_IMpippim_wSid = new TCanvas("cMMnmiss_IMpippim_wSid","cMMnmiss_IMpippim_wSid",1200,800);
  cMMnmiss_IMpippim_wSid->Divide(2,1);
  cMMnmiss_IMpippim_wSid->cd(1);
  MMnmiss_IMpippim_wSid[0]->SetTitle("#splitline{MMnmiss_IMpippim_wSid}{  Real data}");
  MMnmiss_IMpippim_wSid[0]->SetMinimum(1);
  MMnmiss_IMpippim_wSid[0]->Draw("colz");
  cMMnmiss_IMpippim_wSid->cd(2);
  MMnmiss_IMpippim_wSid_mc->SetTitle("#splitline{MMnmiss_IMpippim_wSid}{  MC sum}");
  MMnmiss_IMpippim_wSid_mc->SetMinimum(1);
  MMnmiss_IMpippim_wSid_mc->Draw("colz");

  TH1D* IMpippim_wSid_mc = (TH1D*)MMnmiss_IMpippim_wSid_mc->ProjectionX("IMpippim_wSid_mc");
  IMpippim_wSid_mc->SetLineColor(6);
  TCanvas *cIMpippim_wSid = new TCanvas("cIMpippim_wSid","cIMpippim_wSid");
  cIMpippim_wSid->cd();
  TH1D* IMpippim_wSid[6];
  for(int i=0; i<6; i++) IMpippim_wSid[i] = (TH1D*)MMnmiss_IMpippim_wSid[i]->ProjectionX(Form("IMpippim_wSid_%s",name[i]));
  IMpippim_wSid[0]->Draw("HE");
  IMpippim_wSid_mc->Draw("HEsame");
  for(int i=1; i<6; i++) {
    IMpippim_wSid[i]->SetLineColor(colordef[i]);
    //IMpippim_wSid[i]->Draw("HEsame");
  }


  //MMnmiss vs IMpippim w/o K0 w/o (Sid & n);
  TH2D* MMnmiss_IMpippim_woK0_woSidn_mc = (TH2D*)MMnmiss_IMpippim_woK0_woSidn[1]->Clone("MMnmiss_IMpippim_woK0_woSidn_mc");
  for(int i=1; i<6; i++) MMnmiss_IMpippim_woK0_woSidn_mc->Add(MMnmiss_IMpippim_woK0_woSidn[i]);

  TCanvas *cMMnmiss_IMpippim_woK0_woSidn = new TCanvas("cMMnmiss_IMpippim_woK0_woSidn","cMMnmiss_IMpippim_woK0_woSidn",1200,800);
  cMMnmiss_IMpippim_woK0_woSidn->Divide(2,1);
  cMMnmiss_IMpippim_woK0_woSidn->cd(1);
  MMnmiss_IMpippim_woK0_woSidn[0]->SetTitle("#splitline{MMnmiss_IMpippim_woK0_woSidn}{  Real data}");
  MMnmiss_IMpippim_woK0_woSidn[0]->SetMinimum(1);
  MMnmiss_IMpippim_woK0_woSidn[0]->Draw("colz");
  cMMnmiss_IMpippim_woK0_woSidn->cd(2);
  MMnmiss_IMpippim_woK0_woSidn_mc->SetTitle("#splitline{MMnmiss_IMpippim_woK0_woSidn}{  MC sum}");
  MMnmiss_IMpippim_woK0_woSidn_mc->SetMinimum(1);
  MMnmiss_IMpippim_woK0_woSidn_mc->Draw("colz");

  TH1D* IMpippim_woK0_woSidn_mc = (TH1D*)MMnmiss_IMpippim_woK0_woSidn_mc->ProjectionX("IMpippim_woK0_woSidn_mc");
  IMpippim_woK0_woSidn_mc->SetLineColor(6);
  TCanvas *cIMpippim_woK0_woSidn = new TCanvas("cIMpippim_woK0_woSidn","cIMpippim_woK0_woSidn");
  cIMpippim_woK0_woSidn->cd();
  TH1D* IMpippim_woK0_woSidn[6];
  for(int i=0; i<6; i++) IMpippim_woK0_woSidn[i] = (TH1D*)MMnmiss_IMpippim_woK0_woSidn[i]->ProjectionX(Form("IMpippim_woK0_woSidn_%s",name[i]));
  IMpippim_woK0_woSidn[0]->Draw("HE");
  IMpippim_woK0_woSidn_mc->Draw("HEsame");
  for(int i=1; i<6; i++) {
    IMpippim_woK0_woSidn[i]->SetLineColor(colordef[i]);
    //IMpippim_woK0_woSidn[i]->Draw("HEsame");
  }

  TCanvas *cIMpippim_woK0_woSidn_ratio = new TCanvas("cIMpippim_woK0_woSidn_ratio","cIMpippim_woK0_woSidn_ratio");
  cIMpippim_woK0_woSidn_ratio->cd();
  TH1D* IMpippim_woK0_woSidn_ratio = (TH1D*)IMpippim_woK0_woSidn[0]->Clone("IMpippim_woK0_woSidn_ratio");
  IMpippim_woK0_woSidn_ratio->Divide(IMpippim_woK0_woSidn_mc);
  IMpippim_woK0_woSidn_ratio->SetTitle("Data/MC");
  IMpippim_woK0_woSidn_ratio->GetYaxis()->SetRangeUser(0,5);
  IMpippim_woK0_woSidn_ratio->Draw("HE");

  TF1 *evalf_IMpippim = new TF1("evalf_IMpippim","pol6",0,1);
  evalf_IMpippim->SetTitle("IMpippim");
  IMpippim_woK0_woSidn_ratio->Fit("evalf_IMpippim","","",0.28,0.97);

  //q vs IMnpipi w/o K0 w/o (Sid & n);
  TH2D* q_IMnpipi_woK0_woSidn_mc = (TH2D*)q_IMnpipi_woK0_woSidn[1]->Clone("q_IMnpipi_woK0_woSidn_mc");
  for(int i=2; i<6; i++)q_IMnpipi_woK0_woSidn_mc->Add(q_IMnpipi_woK0_woSidn[i]);

  TCanvas *cq_IMnpipi_woK0_woSidn = new TCanvas("cq_IMnpipi_woK0_woSidn","q_IMnpipi_woK0_woSidn",1200,800);
  cq_IMnpipi_woK0_woSidn->Divide(2,1);
  cq_IMnpipi_woK0_woSidn->cd(1);
  q_IMnpipi_woK0_woSidn[0]->SetTitle("#splitline{q_IMnpipi_woK0_woSidn}{  Real data}");
  q_IMnpipi_woK0_woSidn[0]->SetMinimum(1);
  q_IMnpipi_woK0_woSidn[0]->Draw("colz");
  cq_IMnpipi_woK0_woSidn->cd(2);
  q_IMnpipi_woK0_woSidn_mc->SetTitle("#splitline{q_IMnpipi_woK0_woSidn}{  MC sum}");
  q_IMnpipi_woK0_woSidn_mc->SetMinimum(1);
  q_IMnpipi_woK0_woSidn_mc->Draw("colz");

  TH1D* IMnpipi_woK0_woSidn_mc = (TH1D*)q_IMnpipi_woK0_woSidn_mc->ProjectionX("IMnpipi_woK0_woSidn_mc");
  IMnpipi_woK0_woSidn_mc->SetLineColor(6);
  TCanvas *cIMnpipi_woK0_woSidn = new TCanvas("cIMnpipi_woK0_woSidn","cIMnpipi_woK0_woSidn");
  cIMnpipi_woK0_woSidn->cd();
  TH1D* IMnpipi_woK0_woSidn[6];
  for(int i=0; i<6; i++) IMnpipi_woK0_woSidn[i] = (TH1D*)q_IMnpipi_woK0_woSidn[i]->ProjectionX(Form("IMnpipi_woK0_woSidn_%s",name[i]));
  IMnpipi_woK0_woSidn[0]->Draw("HE");
  IMnpipi_woK0_woSidn_mc->Draw("HEsame");
  for(int i=1; i<6; i++) {
    IMnpipi_woK0_woSidn[i]->SetLineColor(colordef[i]);
    //IMnpipi_woK0_woSidn[i]->Draw("HEsame");
  }

  TCanvas *cIMnpipi_woK0_woSidn_ratio = new TCanvas("cIMnpipi_woK0_woSidn_ratio","cIMnpipi_woK0_woSidn_ratio");
  cIMnpipi_woK0_woSidn_ratio->cd();
  TH1D* IMnpipi_woK0_woSidn_ratio = (TH1D*)IMnpipi_woK0_woSidn[0]->Clone("IMnpipi_woK0_woSidn_ratio");
  IMnpipi_woK0_woSidn_ratio->Divide(IMnpipi_woK0_woSidn_mc);
  IMnpipi_woK0_woSidn_ratio->SetTitle("Data/MC");
  IMnpipi_woK0_woSidn_ratio->GetYaxis()->SetRangeUser(0,3);
  IMnpipi_woK0_woSidn_ratio->Draw("HE");

  TF1 *evalf_IMnpipi = new TF1("evalf_IMnpipi","pol5",1.22,2.00);
  evalf_IMnpipi->SetTitle("IMnpipi");
  IMnpipi_woK0_woSidn_ratio->Fit("evalf_IMnpipi","","",1.22,2.00);


  TH1D* q_woK0_woSidn_mc = (TH1D*)q_IMnpipi_woK0_woSidn_mc->ProjectionY("q_woK0_woSidn_mc");
  q_woK0_woSidn_mc->SetLineColor(6);
  TCanvas *cq_woK0_woSidn = new TCanvas("cq_woK0_woSidn","cq_woK0_woSidn");
  cq_woK0_woSidn->cd();
  TH1D* q_woK0_woSidn[6];
  for(int i=0; i<6; i++) q_woK0_woSidn[i] = (TH1D*)q_IMnpipi_woK0_woSidn[i]->ProjectionY(Form("q_woK0_woSidn_%s",name[i]));
  q_woK0_woSidn[0]->Draw("HE");
  q_woK0_woSidn_mc->Draw("HEsame");
  for(int i=1; i<6; i++) {
    q_woK0_woSidn[i]->SetLineColor(colordef[i]);
    //q_woK0_woSidn[i]->Draw("HEsame");
  }

  //missing mass vs Mom(pi+pi-) w/o K0 w/o (Sid & n)
  TH2D* MMnmiss_Mompippim_woK0_woSidn_mc = (TH2D*)MMnmiss_Mompippim_woK0_woSidn[1]->Clone("MMnmiss_Mompippim_woK0_woSidn_mc");
  for(int i=1; i<6; i++)MMnmiss_Mompippim_woK0_woSidn_mc->Add(MMnmiss_Mompippim_woK0_woSidn[i]);

  TCanvas *cMMnmiss_Mompippim_woK0_woSidn = new TCanvas("cMMnmiss_Mompippim_woK0_woSidn","cMMnmiss_Mompippim_woK0_woSidn",1200,800);
  cMMnmiss_Mompippim_woK0_woSidn->Divide(2,1);
  cMMnmiss_Mompippim_woK0_woSidn->cd(1);
  MMnmiss_Mompippim_woK0_woSidn[0]->SetTitle("#splitline{MMnmiss_Mompippim_woK0_woSidn}{  Real data}");
  MMnmiss_Mompippim_woK0_woSidn[0]->SetMinimum(1);
  MMnmiss_Mompippim_woK0_woSidn[0]->Draw("colz");
  cMMnmiss_Mompippim_woK0_woSidn->cd(2);
  MMnmiss_Mompippim_woK0_woSidn_mc->SetTitle("#splitline{MMnmiss_Mompippim_woK0_woSidn}{  MC sum}");
  MMnmiss_Mompippim_woK0_woSidn_mc->SetMinimum(1);;
  MMnmiss_Mompippim_woK0_woSidn_mc->Draw("colz");

  TH1D* Mompippim_woK0_woSidn_mc = (TH1D*)MMnmiss_Mompippim_woK0_woSidn_mc->ProjectionX("Mompippim_woK0_woSidn_mc");
  Mompippim_woK0_woSidn_mc->SetLineColor(6);
  TCanvas *cMompippim_woK0_woSidn = new TCanvas("cMompippim_woK0_woSidn","cMompippim_woK0_woSidn");
  cMompippim_woK0_woSidn->cd();
  TH1D* Mompippim_woK0_woSidn[6];
  for(int i=0; i<6; i++) Mompippim_woK0_woSidn[i] = (TH1D*)MMnmiss_Mompippim_woK0_woSidn[i]->ProjectionX(Form("MMnmiss_Mompippim_woK0_woSidn_%s",name[i]));
  Mompippim_woK0_woSidn[0]->Draw("HE");
  Mompippim_woK0_woSidn_mc->Draw("HEsame");
  for(int i=1; i<6; i++) {
    Mompippim_woK0_woSidn[i]->SetLineColor(colordef[i]);
    //Mompippim_woK0_woSidn[i]->Draw("HEsame");
  }

  TCanvas *cMompippim_woK0_woSidn_ratio = new TCanvas("cMompippim_woK0_woSidn_ratio","cMompippim_woK0_woSidn_ratio");
  TH1D* Mompippim_woK0_woSidn_ratio = (TH1D*)Mompippim_woK0_woSidn[0]->Clone("Mompippim_woK0_woSidn_ratio");
  Mompippim_woK0_woSidn_ratio->Divide(Mompippim_woK0_woSidn_mc);
  Mompippim_woK0_woSidn_ratio->SetTitle("Data/MC");
  Mompippim_woK0_woSidn_ratio->GetYaxis()->SetRangeUser(0,3);
  Mompippim_woK0_woSidn_ratio->Draw("HE");

  TF1 *evalf_Mompippim = new TF1("evalf_Mompippim","pol3",0,1);
  evalf_Mompippim->SetTitle("Mompippim");
  Mompippim_woK0_woSidn_ratio->Fit("evalf_Mompippim","","",0,0.97);
  

  //
  TH2D* IMnpim_IMnpip_woK0_woSidn_mc = (TH2D*)IMnpim_IMnpip_woK0_woSidn[1]->Clone("IMnpim_IMnpip_woK0_woSidn_mc");
  for(int i=2; i<6; i++)IMnpim_IMnpip_woK0_woSidn_mc->Add(IMnpim_IMnpip_woK0_woSidn[i]);

  TCanvas* cIMnpim_IMnpip_woK0_woSidn = new TCanvas("cIMnpim_IMnpip_woK0_woSidn","cIMnpim_IMnpip_woK0_woSidn");
  cIMnpim_IMnpip_woK0_woSidn->Divide(2,1);
  cIMnpim_IMnpip_woK0_woSidn->cd(1);
  IMnpim_IMnpip_woK0_woSidn[0]->Draw("colz");
  cIMnpim_IMnpip_woK0_woSidn->cd(2);
  IMnpim_IMnpip_woK0_woSidn_mc->Draw("colz");

  //pipmom
  TH2D* pipmom_MMnmiss_woK0_woSidn_mc = (TH2D*)pipmom_MMnmiss_woK0_woSidn[1]->Clone("pipmom_MMnmiss_woK0_woSidn_mc");
  for(int i=2; i<6; i++) pipmom_MMnmiss_woK0_woSidn_mc->Add(pipmom_MMnmiss_woK0_woSidn[i]);

  TH1D* pipmom_woK0_woSidn_mc = (TH1D*)pipmom_MMnmiss_woK0_woSidn_mc->ProjectionY("pipmom_woK0_woSidn_mc");
  TH1D* pipmom_woK0_woSidn[6];
  for(int i=0; i<6; i++) {
    pipmom_woK0_woSidn[i] = (TH1D*)pipmom_MMnmiss_woK0_woSidn[i]->ProjectionY(Form("pipmom_woK0_woSidn_%s",name[i]));
    pipmom_woK0_woSidn[i]->SetLineColor(colordef[i]);
  }
  TCanvas *cpipmom_woK0_woSidn = new TCanvas("cpipmom_woK0_woSidn","cpipmom_woK0_woSidn");
  cpipmom_woK0_woSidn->cd();
  pipmom_woK0_woSidn[0]->Draw("HE");
  pipmom_woK0_woSidn_mc->SetLineColor(6);
  pipmom_woK0_woSidn_mc->Draw("HEsame");
  //for(int i=1;i<6;i++)pipmom_woK0_woSidn[i]->Draw("same");

  TCanvas *cpipmom_woK0_woSidn_ratio = new TCanvas("cpipmom_woK0_woSidn_ratio","cpipmom_woK0_woSidn_ratio");
  cpipmom_woK0_woSidn_ratio->cd();
  TH1D* pipmom_woK0_woSidn_ratio = (TH1D*)pipmom_woK0_woSidn[0]->Clone("pipmom_woK0_woSidn_ratio");
  pipmom_woK0_woSidn_ratio->Divide(pipmom_woK0_woSidn_mc);
  pipmom_woK0_woSidn_ratio->GetYaxis()->SetRangeUser(-1,6);
  pipmom_woK0_woSidn_ratio->SetTitle("Data/MC");
  pipmom_woK0_woSidn_ratio->Draw("HEsame");

  TF1 *evalf_pipmom = new TF1("evalf_pipmom","pol8",0.08,0.7);
  evalf_pipmom->SetTitle("pipmom");
  pipmom_woK0_woSidn_ratio->Fit("evalf_pipmom","","",0.08,0.7);
  

  TH2D* pipmom_pimmom_woK0_woSidn_mc = (TH2D*)pipmom_pimmom_woK0_woSidn[1]->Clone("pipmom_pimmom_woK0_woSidn_mc");
  for(int i=1;i<6;i++)pipmom_pimmom_woK0_woSidn_mc->Add(pipmom_pimmom_woK0_woSidn[i]);
  TCanvas *cpipmom_pimmom_woK0_woSidn = new TCanvas("cpipmom_pimmom_woK0_woSidn","cpipmom_pimmom_woK0_woSidn");
  cpipmom_pimmom_woK0_woSidn->Divide(2,1);
  cpipmom_pimmom_woK0_woSidn->cd(1);
  pipmom_pimmom_woK0_woSidn[0]->Draw("colz");
  cpipmom_pimmom_woK0_woSidn->cd(2);
  pipmom_pimmom_woK0_woSidn_mc->Draw("colz");


  //pimmom
  TH2D* pimmom_MMnmiss_woK0_woSidn_mc = (TH2D*)pimmom_MMnmiss_woK0_woSidn[1]->Clone("pimmom_MMnmiss_woK0_woSidn_mc");
  for(int i=2; i<6; i++) pimmom_MMnmiss_woK0_woSidn_mc->Add(pimmom_MMnmiss_woK0_woSidn[i]);

  TH1D* pimmom_woK0_woSidn_mc = (TH1D*)pimmom_MMnmiss_woK0_woSidn_mc->ProjectionY("pimmom_woK0_woSidn_mc");
  TH1D* pimmom_woK0_woSidn[6];
  for(int i=0; i<6; i++) {
    pimmom_woK0_woSidn[i] = (TH1D*)pimmom_MMnmiss_woK0_woSidn[i]->ProjectionY(Form("pimmom_woK0_woSidn_%s",name[i]));
    pimmom_woK0_woSidn[i]->SetLineColor(colordef[i]);
  }
  TCanvas *cpimmom_woK0_woSidn = new TCanvas("cpimmom_woK0_woSidn","cpimmom_woK0_woSidn");
  cpimmom_woK0_woSidn->cd();
  pimmom_woK0_woSidn[0]->Draw("HE");
  pimmom_woK0_woSidn_mc->SetLineColor(6);
  pimmom_woK0_woSidn_mc->Draw("HEsame");
  //for(int i=1;i<6;i++)pimmom_woK0_woSidn[i]->Draw("same");

  TCanvas *cpimmom_woK0_woSidn_ratio = new TCanvas("cpimmom_woK0_woSidn_ratio","cpimmom_woK0_woSidn_ratio");
  cpimmom_woK0_woSidn_ratio->cd();
  TH1D* pimmom_woK0_woSidn_ratio = (TH1D*)pimmom_woK0_woSidn[0]->Clone("pimmom_woK0_woSidn_ratio");
  pimmom_woK0_woSidn_ratio->Divide(pimmom_woK0_woSidn_mc);
  pimmom_woK0_woSidn_ratio->GetYaxis()->SetRangeUser(-1,6);
  pimmom_woK0_woSidn_ratio->SetTitle("Data/MC");
  pimmom_woK0_woSidn_ratio->Draw("HEsame");

  TF1 *evalf_pimmom = new TF1("evalf_pimmom","pol8",0.06,0.73);
  evalf_pimmom->SetTitle("pimmom");
  pimmom_woK0_woSidn_ratio->Fit("evalf_pimmom","","",0.06,0.73);

  //nCDSmom
  TH2D* nmom_MMnmiss_woK0_woSidn_mc = (TH2D*)nmom_MMnmiss_woK0_woSidn[1]->Clone("nmom_MMnmiss_woK0_woSidn_mc");
  for(int i=2; i<6; i++) nmom_MMnmiss_woK0_woSidn_mc->Add(nmom_MMnmiss_woK0_woSidn[i]);

  TH1D* nmom_woK0_woSidn_mc = (TH1D*)nmom_MMnmiss_woK0_woSidn_mc->ProjectionY("nmom_woK0_woSidn_mc");
  TH1D* nmom_woK0_woSidn[6];
  for(int i=0; i<6; i++) {
    nmom_woK0_woSidn[i] = (TH1D*)nmom_MMnmiss_woK0_woSidn[i]->ProjectionY(Form("nmom_woK0_woSidn_%s",name[i]));
    nmom_woK0_woSidn[i]->SetLineColor(colordef[i]);
  }
  TCanvas *cnmom_woK0_woSidn = new TCanvas("cnmom_woK0_woSidn","cnmom_woK0_woSidn");
  cnmom_woK0_woSidn->cd();
  nmom_woK0_woSidn[0]->Draw("HE");
  nmom_woK0_woSidn_mc->SetLineColor(6);
  nmom_woK0_woSidn_mc->Draw("HEsame");
  //for(int i=1;i<6;i++)nmom_woK0_woSidn[i]->Draw("same");

  TCanvas *cnmom_woK0_woSidn_ratio = new TCanvas("cnmom_woK0_woSidn_ratio","cnmom_woK0_woSidn_ratio");
  cnmom_woK0_woSidn_ratio->cd();
  TH1D* nmom_woK0_woSidn_ratio = (TH1D*)nmom_woK0_woSidn[0]->Clone("nmom_woK0_woSidn_ratio");
  nmom_woK0_woSidn_ratio->Divide(nmom_woK0_woSidn_mc);
  nmom_woK0_woSidn_ratio->SetTitle("Data/MC");
  nmom_woK0_woSidn_ratio->GetYaxis()->SetRangeUser(-1,4);
  nmom_woK0_woSidn_ratio->Draw("HEsame");

  TF1 *evalf_nmom = new TF1("evalf_nmom","pol8",0.14,1);
  nmom_woK0_woSidn_ratio->Fit("evalf_nmom","","",0.14,1);

  
  //weighting to momentum vector 
  TH2D* nmom_cosn_woK0_woSidn_mc = (TH2D*)nmom_cosn_woK0_woSidn[1]->Clone("nmom_cosn_woK0_woSidn_mc");
  for(int i=2;i<6;i++) nmom_cosn_woK0_woSidn_mc->Add(nmom_cosn_woK0_woSidn[i]);
  TCanvas *cnmom_cosn_woK0_woSidn = new TCanvas("cnmom_cosn_woK0_woSidn","cnmom_cosn_woK0_woSidn",1200,800);
  cnmom_cosn_woK0_woSidn->Divide(2,1);
  cnmom_cosn_woK0_woSidn->cd(1);
  nmom_cosn_woK0_woSidn[0]->SetTitle("#splitline{nmom_cosn_woK0_woSidn}{Real data}");
  nmom_cosn_woK0_woSidn[0]->Draw("colz");
  cnmom_cosn_woK0_woSidn->cd(2);
  nmom_cosn_woK0_woSidn_mc->SetTitle("#splitline{nmom_cosn_woK0_woSidn}{MC}");
  nmom_cosn_woK0_woSidn_mc->Draw("colz");

  TCanvas *ccosn_woK0_woSidn = new TCanvas("ccosn_woK0_woSidn","ccosn_woK0_woSidn");
  ccosn_woK0_woSidn->cd();
  TH1D* cosn_woK0_woSidn[6];
  for(int i=0;i<6;i++) cosn_woK0_woSidn[i] = (TH1D*)nmom_cosn_woK0_woSidn[i]->ProjectionX(Form("cosn_woK0_woSidn_%s",name[i]));
  TH1D* cosn_woK0_woSidn_mc = (TH1D*)nmom_cosn_woK0_woSidn_mc->ProjectionX("cosn_woK0_woSidn_mc");
  cosn_woK0_woSidn_mc->SetLineColor(6);
  cosn_woK0_woSidn[0]->Draw("HE");
  cosn_woK0_woSidn_mc->Draw("HEsame");
   

  TCanvas *ccosn_woK0_woSidn_ratio = new TCanvas("ccosn_woK0_woSidn_ratio","ccosn_woK0_woSidn_ratio");
  ccosn_woK0_woSidn_ratio->cd();
  TH1D* cosn_woK0_woSidn_ratio = (TH1D*)cosn_woK0_woSidn[0]->Clone("cosn_woK0_woSidn_ratio");
  cosn_woK0_woSidn_ratio->Divide(cosn_woK0_woSidn_mc);
  cosn_woK0_woSidn_ratio->SetTitle("Data/MC");
  cosn_woK0_woSidn_ratio->GetYaxis()->SetRangeUser(0,4);
  cosn_woK0_woSidn_ratio->Draw("HE");
  
  TF1 *evalf_cosn_1 = new TF1("evalf_cosn_1","expo",-1,-0.9);
  cosn_woK0_woSidn_ratio->Fit(evalf_cosn_1,"R","",-1,-0.9);
   
  TF1 *evalf_cosn_2 = new TF1("evalf_cosn_2","pol3",-0.9,0.3);
  cosn_woK0_woSidn_ratio->Fit(evalf_cosn_2,"R+","",-0.9,0.3);

  Double_t param_cosn[6];
  evalf_cosn_1->GetParameters(&param_cosn[0]);
  evalf_cosn_2->GetParameters(&param_cosn[2]);
  TF1 *evalf_cosn = new TF1("evalf_cosn",func_cosn,-1.00,0.3,6);
  evalf_cosn->SetParameters(param_cosn);
  evalf_cosn->SetLineColor(3);
  evalf_cosn->Draw("same");
  

  //
  TH2D* nmom_cospip_woK0_woSidn_mc = (TH2D*)nmom_cospip_woK0_woSidn[1]->Clone("nmom_cospip_woK0_woSidn_mc");
  for(int i=2;i<6;i++)nmom_cospip_woK0_woSidn_mc->Add(nmom_cospip_woK0_woSidn[i]);
  TCanvas *cnmom_cospip_woK0_woSidn = new TCanvas("cnmom_cospip_woK0_woSidn","cnmom_cospip_woK0_woSidn",1200,800);
  cnmom_cospip_woK0_woSidn->Divide(2,1);
  cnmom_cospip_woK0_woSidn->cd(1);
  nmom_cospip_woK0_woSidn[0]->SetTitle("#splitline{nmom_cospip_woK0_woSidn}{Real data}");
  nmom_cospip_woK0_woSidn[0]->Draw("colz");
  cnmom_cospip_woK0_woSidn->cd(2);
  nmom_cospip_woK0_woSidn_mc->SetTitle("#splitline{nmom_cospip_woK0_woSidn}{MC}");
  nmom_cospip_woK0_woSidn_mc->Draw("colz");

  TCanvas *ccospip_woK0_woSidn = new TCanvas("ccospip_woK0_woSidn","ccospip_woK0_woSidn");
  ccospip_woK0_woSidn->cd();
  TH1D* cospip_woK0_woSidn[6];
  for(int i=0;i<6;i++) cospip_woK0_woSidn[i] = (TH1D*)nmom_cospip_woK0_woSidn[i]->ProjectionX(Form("cospip_woK0_woSidn_%s",name[i]));
  TH1D* cospip_woK0_woSidn_mc = (TH1D*)nmom_cospip_woK0_woSidn_mc->ProjectionX("cospip_woK0_woSidn_mc");
  cospip_woK0_woSidn_mc->SetLineColor(6);
  cospip_woK0_woSidn[0]->Draw("HE");
  cospip_woK0_woSidn_mc->Draw("HEsame");

  TCanvas *ccospip_woK0_woSidn_ratio = new TCanvas("ccospip_woK0_woSidn_ratio","ccospip_woK0_woSidn_ratio");
  ccospip_woK0_woSidn_ratio->cd();
  TH1D* cospip_woK0_woSidn_ratio = (TH1D*)cospip_woK0_woSidn[0]->Clone("cospip_woK0_woSidn_ratio");
  cospip_woK0_woSidn_ratio->Divide(cospip_woK0_woSidn_mc);
  cospip_woK0_woSidn_ratio->SetTitle("Data/MC");
  cospip_woK0_woSidn_ratio->GetYaxis()->SetRangeUser(0,4);
  cospip_woK0_woSidn_ratio->Draw("HE");
  
  TF1 *evalf_cospip = new TF1("evalf_cospip",func_cospip,-1.0,0.6,6); 
  cospip_woK0_woSidn_ratio->Fit("evalf_cospip","","",-1.0,0.6);


  //
  TH2D* nmom_cospim_woK0_woSidn_mc = (TH2D*)nmom_cospim_woK0_woSidn[1]->Clone("nmom_cospim_woK0_woSidn_mc");
  for(int i=2;i<6;i++)nmom_cospim_woK0_woSidn_mc->Add(nmom_cospim_woK0_woSidn[i]);
  TCanvas *cnmom_cospim_woK0_woSidn = new TCanvas("cnmom_cospim_woK0_woSidn","cnmom_cospim_woK0_woSidn",1200,800);
  cnmom_cospim_woK0_woSidn->Divide(2,1);
  cnmom_cospim_woK0_woSidn->cd(1);
  nmom_cospim_woK0_woSidn[0]->SetTitle("#splitline{nmom_cospim_woK0_woSidn}{Real data}");
  nmom_cospim_woK0_woSidn[0]->Draw("colz");
  cnmom_cospim_woK0_woSidn->cd(2);
  nmom_cospim_woK0_woSidn_mc->SetTitle("#splitline{nmom_cospim_woK0_woSidn}{MC}");
  nmom_cospim_woK0_woSidn_mc->Draw("colz");

  TCanvas *ccospim_woK0_woSidn = new TCanvas("ccospim_woK0_woSidn","ccospim_woK0_woSidn");
  ccospim_woK0_woSidn->cd();
  TH1D* cospim_woK0_woSidn[6];
  for(int i=0;i<6;i++) cospim_woK0_woSidn[i] = (TH1D*)nmom_cospim_woK0_woSidn[i]->ProjectionX(Form("cospim_woK0_woSidn_%s",name[i]));
  TH1D* cospim_woK0_woSidn_mc = (TH1D*)nmom_cospim_woK0_woSidn_mc->ProjectionX("cospim_woK0_woSidn_mc");
  cospim_woK0_woSidn_mc->SetLineColor(6);
  cospim_woK0_woSidn[0]->Draw("HE");
  cospim_woK0_woSidn_mc->Draw("HEsame");

  TCanvas *ccospim_woK0_woSidn_ratio = new TCanvas("ccospim_woK0_woSidn_ratio","ccospim_woK0_woSidn_ratio");
  ccospim_woK0_woSidn_ratio->cd();
  TH1D* cospim_woK0_woSidn_ratio = (TH1D*)cospim_woK0_woSidn[0]->Clone("cospim_woK0_woSidn_ratio");
  cospim_woK0_woSidn_ratio->Divide(cospim_woK0_woSidn_mc);
  cospim_woK0_woSidn_ratio->SetTitle("Data/MC");
  cospim_woK0_woSidn_ratio->GetYaxis()->SetRangeUser(0,4);
  cospim_woK0_woSidn_ratio->Draw("HE");
   
  TF1 *evalf_cospim = new TF1("evalf_cospim","pol8",-0.92,0.50);
  cospim_woK0_woSidn_ratio->Fit(evalf_cospim,"","",-0.92,0.50);
 
  //
  TH2D* nmom_phinpip_woK0_woSidn_mc = (TH2D*)nmom_phinpip_woK0_woSidn[1]->Clone("nmom_phinpip_woK0_woSidn_mc");
  for(int i=2;i<6;i++)nmom_phinpip_woK0_woSidn_mc->Add(nmom_phinpip_woK0_woSidn[i]);
  TCanvas *cnmom_phinpip_woK0_woSidn = new TCanvas("cnmom_phinpip_woK0_woSidn","cnmom_phinpip_woK0_woSidn",1200,800);
  cnmom_phinpip_woK0_woSidn->Divide(2,1);
  cnmom_phinpip_woK0_woSidn->cd(1);
  nmom_phinpip_woK0_woSidn[0]->SetTitle("#splitline{nmom_phinpip_woK0_woSidn}{Real data}");
  nmom_phinpip_woK0_woSidn[0]->Draw("colz");
  cnmom_phinpip_woK0_woSidn->cd(2);
  nmom_phinpip_woK0_woSidn_mc->SetTitle("#splitline{nmom_phinpip_woK0_woSidn}{MC}");
  nmom_phinpip_woK0_woSidn_mc->Draw("colz");

  TCanvas *cphinpip_woK0_woSidn = new TCanvas("cphinpip_woK0_woSidn","cphinpip_woK0_woSidn");
  cphinpip_woK0_woSidn->cd();
  TH1D* phinpip_woK0_woSidn[6];
  for(int i=0;i<6;i++) phinpip_woK0_woSidn[i] = (TH1D*)nmom_phinpip_woK0_woSidn[i]->ProjectionX(Form("phinpip_woK0_woSidn_%s",name[i]));
  TH1D* phinpip_woK0_woSidn_mc = (TH1D*)nmom_phinpip_woK0_woSidn_mc->ProjectionX("phinpip_woK0_woSidn_mc");
  phinpip_woK0_woSidn_mc->SetLineColor(6);
  phinpip_woK0_woSidn[0]->Draw("HE");
  phinpip_woK0_woSidn_mc->Draw("HEsame");

  TCanvas *cphinpip_woK0_woSidn_ratio = new TCanvas("cphinpip_woK0_woSidn_ratio","cphinpip_woK0_woSidn_ratio");
  cphinpip_woK0_woSidn_ratio->cd();
  TH1D* phinpip_woK0_woSidn_ratio = (TH1D*)phinpip_woK0_woSidn[0]->Clone("phinpip_woK0_woSidn_ratio");
  phinpip_woK0_woSidn_ratio->Divide(phinpip_woK0_woSidn_mc);
  phinpip_woK0_woSidn_ratio->SetTitle("Data/MC");
  phinpip_woK0_woSidn_ratio->GetYaxis()->SetRangeUser(0,4);
  phinpip_woK0_woSidn_ratio->Draw("HE");
   
  //
  TH2D* nmom_phinpim_woK0_woSidn_mc = (TH2D*)nmom_phinpim_woK0_woSidn[1]->Clone("nmom_phinpim_woK0_woSidn_mc");
  for(int i=2;i<6;i++)nmom_phinpim_woK0_woSidn_mc->Add(nmom_phinpim_woK0_woSidn[i]);
  TCanvas *cnmom_phinpim_woK0_woSidn = new TCanvas("cnmom_phinpim_woK0_woSidn","cnmom_phinpim_woK0_woSidn",1200,800);
  cnmom_phinpim_woK0_woSidn->Divide(2,1);
  cnmom_phinpim_woK0_woSidn->cd(1);
  nmom_phinpim_woK0_woSidn[0]->SetTitle("#splitline{nmom_phinpim_woK0_woSidn}{Real data}");
  nmom_phinpim_woK0_woSidn[0]->Draw("colz");
  cnmom_phinpim_woK0_woSidn->cd(2);
  nmom_phinpim_woK0_woSidn_mc->SetTitle("#splitline{nmom_phinpim_woK0_woSidn}{MC}");
  nmom_phinpim_woK0_woSidn_mc->Draw("colz");

  TCanvas *cphinpim_woK0_woSidn = new TCanvas("cphinpim_woK0_woSidn","cphinpim_woK0_woSidn");
  cphinpim_woK0_woSidn->cd();
  TH1D* phinpim_woK0_woSidn[6];
  for(int i=0;i<6;i++) phinpim_woK0_woSidn[i] = (TH1D*)nmom_phinpim_woK0_woSidn[i]->ProjectionX(Form("phinpim_woK0_woSidn_%s",name[i]));
  TH1D* phinpim_woK0_woSidn_mc = (TH1D*)nmom_phinpim_woK0_woSidn_mc->ProjectionX("phinpim_woK0_woSidn_mc");
  phinpim_woK0_woSidn_mc->SetLineColor(6);
  phinpim_woK0_woSidn[0]->Draw("HE");
  phinpim_woK0_woSidn_mc->Draw("HEsame");

  TCanvas *cphinpim_woK0_woSidn_ratio = new TCanvas("cphinpim_woK0_woSidn_ratio","cphinpim_woK0_woSidn_ratio");
  cphinpim_woK0_woSidn_ratio->cd();
  TH1D* phinpim_woK0_woSidn_ratio = (TH1D*)phinpim_woK0_woSidn[0]->Clone("phinpim_woK0_woSidn_ratio");
  phinpim_woK0_woSidn_ratio->Divide(phinpim_woK0_woSidn_mc);
  phinpim_woK0_woSidn_ratio->SetTitle("Data/MC");
  phinpim_woK0_woSidn_ratio->GetYaxis()->SetRangeUser(0,4);
  phinpim_woK0_woSidn_ratio->Draw("HE");





  ////////////////////////////
  // w/K0 + fake neutron
  ///////////////////////////

  //w K0, w/o ((Sp or Sm) & missing neutron)
  TH2D *MMnmiss_IMnpip_wK0_woSid_won_mc = (TH2D*)MMnmiss_IMnpip_wK0_woSid_won[1]->Clone("MMnmiss_IMnpip_wK0_woSid_won_mc");
  for(int i=2; i<6; i++) MMnmiss_IMnpip_wK0_woSid_won_mc->Add(MMnmiss_IMnpip_wK0_woSid_won[i]);

  TCanvas *cMMnmiss_IMnpip_wK0_woSid_won = new TCanvas("cMMnmiss_IMnpip_wK0_woSid_won","cMMnmiss_IMnpip_wK0_woSid_won",1200,800);
  cMMnmiss_IMnpip_wK0_woSid_won->Divide(2,1);
  cMMnmiss_IMnpip_wK0_woSid_won->cd(1);
  MMnmiss_IMnpip_wK0_woSid_won_rdata->SetTitle("#splitline{MMnmiss_IMnpip_wK0_woSid_won}{  Real data}");
  MMnmiss_IMnpip_wK0_woSid_won_rdata->Draw("colz");
  cMMnmiss_IMnpip_wK0_woSid_won->cd(2);
  MMnmiss_IMnpip_wK0_woSid_won_mc->SetTitle("#splitline{MMnmiss_IMnpip_wK0_woSid_won}{  MC sum}");
  MMnmiss_IMnpip_wK0_woSid_won_mc->Draw("colz");

  //projection to Missing mass (miss n & Sigma+/-)
  TCanvas *cMMnmiss_wK0_woSid_won = new TCanvas("cMMnmiss_wK0_woSid_won","cMMnmiss_wK0_woSid_won");
  cMMnmiss_wK0_woSid_won->cd();
  TH1D* MMnmiss_wK0_woSid_won[6];
  for(int i=0; i<6; i++) MMnmiss_wK0_woSid_won[i] = (TH1D*)MMnmiss_IMnpip_wK0_woSid_won[i]->ProjectionY(Form("MMnmiss_wK0_woSid_won_%s",name[i]));
  MMnmiss_wK0_woSid_won[0]->Draw("HE");
  TH1D* MMnmiss_wK0_woSid_won_mc = (TH1D*)MMnmiss_IMnpip_wK0_woSid_won_mc->ProjectionY("MMnmiss_wK0_woSid_won_mc");
  MMnmiss_wK0_woSid_won_mc->SetLineColor(6);
  MMnmiss_wK0_woSid_won_mc->Draw("HEsame");

  for(int i=1; i<6; i++) {
    MMnmiss_wK0_woSid_won[i]->SetLineColor(colordef[i]);
    //MMnmiss_wK0_woSid_won[i]->Draw("HEsame");
  }

  //Data/MC before modifying MC data
  TCanvas *cMMnmiss_wK0_woSid_won_ratio = new TCanvas("cMMnmiss_wK0_woSid_won_ratio","cMMnmiss_wK0_woSid_won_ratio");
  TH1D* MMnmiss_wK0_woSid_won_ratio = (TH1D*)MMnmiss_wK0_woSid_won[0]->Clone("MMnmiss_wK0_woSid_won_ratio");
  MMnmiss_wK0_woSid_won_ratio->Divide(MMnmiss_wK0_woSid_won_mc);
  MMnmiss_wK0_woSid_won_ratio->SetTitle("Data/MC");
  MMnmiss_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  MMnmiss_wK0_woSid_won_ratio->Draw("HE");


  TF1 *sgf_MMnmiss_wK0 = new TF1("sgf_MMnmiss_wK0","[0]*exp(-0.5*pow((x-[1])/([2]+(x<[1])*[3]*(x-[1])),2))");
  TF1 *fgaus_MMnmiss_wK0_high = new TF1("fgaus_MMnmiss_wK0_high","gaus");
  sgf_MMnmiss_wK0->SetParameter(0,1.82171e+00);
  sgf_MMnmiss_wK0->SetParameter(1,8.56016e-01);
  sgf_MMnmiss_wK0->SetParameter(2,6.81677e-01);
  sgf_MMnmiss_wK0->SetLineColor(2);

  MMnmiss_wK0_woSid_won_ratio->Fit("sgf_MMnmiss_wK0","R","",0.,1.116);
  MMnmiss_wK0_woSid_won_ratio->Fit("fgaus_MMnmiss_wK0_high","R+","",1.116,1.5);
  Double_t param_MMnmiss_wK0[7];
  sgf_MMnmiss_wK0->GetParameters(&param_MMnmiss_wK0[0]);
  fgaus_MMnmiss_wK0_high->GetParameters(&param_MMnmiss_wK0[4]);
  TF1 *evalf_MMnmiss_wK0 = new TF1("evalf_MMnmiss_wK0",func_MMnmiss,0,1.5,7);
  evalf_MMnmiss_wK0->SetParameters(param_MMnmiss_wK0);
  evalf_MMnmiss_wK0->SetLineColor(4);
  evalf_MMnmiss_wK0->Draw("same");
  evalf_MMnmiss_wK0->Print();


  TH2D *MMom_MMass_wK0_woSid_won_mc = (TH2D*)MMom_MMass_wK0_woSid_won[1]->Clone("MMom_MMass_wK0_woSid_won_mc");
  for(int i=2; i<6; i++) MMom_MMass_wK0_woSid_won_mc->Add(MMom_MMass_wK0_woSid_won[i]);
  TCanvas *cMMom_MMass_wK0_woSid_won = new TCanvas("cMMom_MMass_wK0_woSid_won","cMMom_MMass_wK0_woSid_won");
  cMMom_MMass_wK0_woSid_won->Divide(2,1);
  cMMom_MMass_wK0_woSid_won->cd(1);
  MMom_MMass_wK0_woSid_won[0]->Draw("colz");
  cMMom_MMass_wK0_woSid_won->cd(2);
  MMom_MMass_wK0_woSid_won_mc->Draw("colz");

  TH1D* MMom_wK0_woSid_won_mc = (TH1D*)MMom_MMass_wK0_woSid_won_mc->ProjectionY("MMom_wK0_woSid_won_mc");
  TH1D* MMom_wK0_woSid_won[6];
  for(int i=0; i<6; i++) {
    MMom_wK0_woSid_won[i] = (TH1D*)MMom_MMass_wK0_woSid_won[i]->ProjectionY(Form("MMom_wK0_woSid_won_%s",name[i]));
    MMom_wK0_woSid_won[i]->SetLineColor(colordef[i]);
  }
  TCanvas *cMMom_wK0_woSid_won = new TCanvas("cMMom_wK0_woSid_won","cMMom_wK0_woSid_won");
  MMom_wK0_woSid_won[0]->Draw("HE");
  MMom_wK0_woSid_won_mc->SetLineColor(6);
  MMom_wK0_woSid_won_mc->Draw("HEsame");
  //for(int i=1;i<6;i++)MMom_wK0_woSid_won[i]->Draw("HEsame");


  TCanvas *cMMom_wK0_woSid_won_ratio = new TCanvas("cMMom_wK0_woSid_won_ratio","cMMom_wK0_woSid_won_ratio");
  cMMom_wK0_woSid_won_ratio->cd();
  TH1D* MMom_wK0_woSid_won_ratio = (TH1D*)MMom_wK0_woSid_won[0]->Clone("MMom_wK0_woSid_won_ratio");
  MMom_wK0_woSid_won_ratio->Divide(MMom_wK0_woSid_won_mc);
  MMom_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,5);
  MMom_wK0_woSid_won_ratio->SetTitle("Data/MC");
  MMom_wK0_woSid_won_ratio->Draw("HE");

  Double_t param_MMom_wK0[5];
  TF1 *evalf_MMom_wK0 = new TF1("evalf_MMom_wK0","pol5",0.4,1.5);
  evalf_MMom_wK0->SetLineColor(4);
  MMom_wK0_woSid_won_ratio->Fit("evalf_MMom_wK0","","",0.4,1.5);


  //projection to IMnpip (miss n & Sigma+/-)
  TCanvas *cIMnpip_wK0_woSid_won = new TCanvas("cIMnpip_wK0_woSid_won","cIMnpip_wK0_woSid_won");
  cIMnpip_wK0_woSid_won->cd();
  TH1D* IMnpip_wK0_woSid_won[6];
  for(int i=0; i<6; i++)IMnpip_wK0_woSid_won[i] = (TH1D*)MMnmiss_IMnpip_wK0_woSid_won[i]->ProjectionX(Form("IMnpip_wK0_woSid_won_%s",name[i]));
  IMnpip_wK0_woSid_won[0]->Draw("HE");//rdata
  TH1D* IMnpip_wK0_woSid_won_mc = (TH1D*)MMnmiss_IMnpip_wK0_woSid_won_mc->ProjectionX("IMnpip_wK0_woSid_won_mc");
  IMnpip_wK0_woSid_won_mc->SetLineColor(6);
  IMnpip_wK0_woSid_won_mc->Draw("HEsame");
  for(int i=1; i<6; i++) {
    IMnpip_wK0_woSid_won[i]->SetLineColor(colordef[i]);
    //IMnpip_wK0_woSid_won[i]->Draw("HEsame");
  }

  TCanvas *cIMnpip_wK0_woSid_won_ratio = new TCanvas("cIMnpip_wK0_woSid_won_ratio","cIMnpip_wK0_woSid_won_ratio");
  cIMnpip_wK0_woSid_won_ratio->cd();
  TH1D* IMnpip_wK0_woSid_won_ratio = (TH1D*)IMnpip_wK0_woSid_won[0]->Clone("IMnpip_wK0_woSid_won_ratio");
  IMnpip_wK0_woSid_won_ratio->Divide(IMnpip_wK0_woSid_won_mc);
  IMnpip_wK0_woSid_won_ratio->SetTitle("IMnpip_wK0_woSid_won Data/MC");
  IMnpip_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  IMnpip_wK0_woSid_won_ratio->Draw("HE");

  TF1* fgaus_IMnpip_wK0_1 = new TF1("fgaus_IMnpip_wK0_1","gaus",1.06,1.10);
  fgaus_IMnpip_wK0_1->SetParameters(0,1.636);
  fgaus_IMnpip_wK0_1->SetParameters(1,1.102);
  fgaus_IMnpip_wK0_1->SetParameters(2,0.02845);
  IMnpip_wK0_woSid_won_ratio->Fit("fgaus_IMnpip_wK0_1","R","",1.08,1.10);

  TF1* fexpo_IMnpip_wK0_2 = new TF1("fexpo_IMnpip_wK0_2","expo",1.10,1.25);
  fexpo_IMnpip_wK0_2->SetParameters(0,1.667);
  fexpo_IMnpip_wK0_2->SetParameters(1,-1.117);
  IMnpip_wK0_woSid_won_ratio->Fit("fexpo_IMnpip_wK0_2","R+","",1.10,1.25);

  TF1* fgaus_IMnpip_wK0_3 = new TF1("fgaus_IMnpip_wK0_3","gaus",1.25,2.0);
  fgaus_IMnpip_wK0_3->SetParameters(0,10.83);
  fgaus_IMnpip_wK0_3->SetParameters(1,5.094);
  fgaus_IMnpip_wK0_3->SetParameters(2,1.87);
  IMnpip_wK0_woSid_won_ratio->Fit("fgaus_IMnpip_wK0_3","R+","",1.25,2.0);

  Double_t param_IMnpip_wK0[8];
  fgaus_IMnpip_wK0_1->GetParameters(&param_IMnpip_wK0[0]);
  fexpo_IMnpip_wK0_2->GetParameters(&param_IMnpip_wK0[3]);
  fgaus_IMnpip_wK0_3->GetParameters(&param_IMnpip_wK0[5]);

  TF1 *evalf_IMnpip_wK0 = new TF1("evalf_IMnpip_wK0",func_IMnpip,1.06,2.0,8);
  evalf_IMnpip_wK0->SetParameters(param_IMnpip_wK0);
  evalf_IMnpip_wK0->SetLineColor(4);
  evalf_IMnpip_wK0->Draw("same");

  //Sigma-
  //w K0, w/o ((Sp or Sm) & missing neutron)
  TH2D *MMnmiss_IMnpim_wK0_woSid_won_mc = (TH2D*)MMnmiss_IMnpim_wK0_woSid_won[1]->Clone("MMnmiss_IMnpim_wK0_woSid_won_mc");
  for(int i=2; i<6; i++)MMnmiss_IMnpim_wK0_woSid_won_mc->Add(MMnmiss_IMnpim_wK0_woSid_won[i]);

  TCanvas *cMMnmiss_IMnpim_wK0_woSid_won = new TCanvas("cMMnmiss_IMnpim_wK0_woSid_won","cMMnmiss_IMnpim_wK0_woSid_won",1200,800);
  cMMnmiss_IMnpim_wK0_woSid_won->Divide(2,1);
  cMMnmiss_IMnpim_wK0_woSid_won->cd(1);
  MMnmiss_IMnpim_wK0_woSid_won_rdata->SetTitle("#splitline{MMnmiss_IMnpim_wK0_woSid_won}{  Real data}");
  MMnmiss_IMnpim_wK0_woSid_won_rdata->Draw("colz");
  cMMnmiss_IMnpim_wK0_woSid_won->cd(2);
  MMnmiss_IMnpim_wK0_woSid_won_mc->SetTitle("#splitline{MMnmiss_IMnpim_wK0_woSid_won}{  MC sum}");
  MMnmiss_IMnpim_wK0_woSid_won_mc->Draw("colz");

  TCanvas *cIMnpim_wK0_woSid_won = new TCanvas("IMnpim_wK0_woSid_won","IMnpim_wK0_woSid_won");
  cIMnpim_wK0_woSid_won->cd();
  TH1D* IMnpim_wK0_woSid_won[6];
  for(int i=0; i<6; i++)IMnpim_wK0_woSid_won[i] = (TH1D*)MMnmiss_IMnpim_wK0_woSid_won[i]->ProjectionX(Form("IMnpim_wK0_woSid_won_%s",name[i]));
  TH1D* IMnpim_wK0_woSid_won_mc = MMnmiss_IMnpim_wK0_woSid_won_mc->ProjectionX("IMnpim_wK0_woSid_won_mc");
  IMnpim_wK0_woSid_won[0]->Draw("HE");
  IMnpim_wK0_woSid_won_mc->SetLineColor(6);
  IMnpim_wK0_woSid_won_mc->Draw("HEsame");
  for(int i=1; i<6; i++) {
    IMnpim_wK0_woSid_won[i]->SetLineColor(colordef[i]);
    //IMnpim_wK0_woSid_won[i]->Draw("HEsame");
  }

  TCanvas *cIMnpim_wK0_woSid_won_ratio = new TCanvas("cIMnpim_wK0_woSid_won_ratio","cIMnpim_wK0_woSid_won_ratio");
  cIMnpim_wK0_woSid_won_ratio->cd();
  TH1D* IMnpim_wK0_woSid_won_ratio = (TH1D*)IMnpim_wK0_woSid_won[0]->Clone("IMnpim_wK0_woSid_won_ratio");
  IMnpim_wK0_woSid_won_ratio->Divide(IMnpim_wK0_woSid_won_mc);
  IMnpim_wK0_woSid_won_ratio->SetTitle("IMnpim_wK0_woSid_won Data/MC");
  IMnpim_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  IMnpim_wK0_woSid_won_ratio->Draw("HE");

  TF1* pol3_IMnpim_wK0_1 = new TF1("pol3_IMnpim_wK0_1","pol3",1.00,1.10);
  pol3_IMnpim_wK0_1->SetParameter(0,-32074.9);
  pol3_IMnpim_wK0_1->SetParameter(1,85205.3);
  pol3_IMnpim_wK0_1->SetParameter(2,-75374.);
  pol3_IMnpim_wK0_1->SetParameter(3,22204.);
  IMnpim_wK0_woSid_won_ratio->Fit("pol3_IMnpim_wK0_1","R","",1.07,1.10);

  TF1* pol3_IMnpim_wK0_2 = new TF1("pol3_IMnpim_wK0_2","pol3",1.10,2.00);
  pol3_IMnpim_wK0_2->SetParameter(0,38.13);
  pol3_IMnpim_wK0_2->SetParameter(1,-62.3139);
  pol3_IMnpim_wK0_2->SetParameter(2,32.172);
  pol3_IMnpim_wK0_2->SetParameter(3,-4.81);
  IMnpim_wK0_woSid_won_ratio->Fit("pol3_IMnpim_wK0_2","R+","",1.10,2.00);

  Double_t param_IMnpim_wK0[8];
  pol3_IMnpim_wK0_1->GetParameters(&param_IMnpim_wK0[0]);
  pol3_IMnpim_wK0_2->GetParameters(&param_IMnpim_wK0[4]);
  TF1 *evalf_IMnpim_wK0 = new TF1("evalf_IMnpim_wK0",func_IMnpim,1.00,2.00,8);
  evalf_IMnpim_wK0->SetParameters(param_IMnpim_wK0);
  evalf_IMnpim_wK0->SetLineColor(4);
  evalf_IMnpim_wK0->Draw("same");


  //MMnmiss vs IMpippim w/o K0 w/o (Sid & n);
  TH2D* MMnmiss_IMpippim_wK0_woSid_won_mc = (TH2D*)MMnmiss_IMpippim_wK0_woSid_won[1]->Clone("MMnmiss_IMpippim_wK0_woSid_won_mc");
  for(int i=2; i<6; i++) MMnmiss_IMpippim_wK0_woSid_won_mc->Add(MMnmiss_IMpippim_wK0_woSid_won[i]);

  TCanvas *cMMnmiss_IMpippim_wK0_woSid_won = new TCanvas("cMMnmiss_IMpippim_wK0_woSid_won","cMMnmiss_IMpippim_wK0_woSid_won",1200,800);
  cMMnmiss_IMpippim_wK0_woSid_won->Divide(2,1);
  cMMnmiss_IMpippim_wK0_woSid_won->cd(1);
  MMnmiss_IMpippim_wK0_woSid_won[0]->SetTitle("#splitline{MMnmiss_IMpippim_wK0_woSid_won}{  Real data}");
  MMnmiss_IMpippim_wK0_woSid_won[0]->Draw("colz");
  cMMnmiss_IMpippim_wK0_woSid_won->cd(2);
  MMnmiss_IMpippim_wK0_woSid_won_mc->SetTitle("#splitline{MMnmiss_IMpippim_wK0_woSid_won}{  MC sum}");
  MMnmiss_IMpippim_wK0_woSid_won_mc->Draw("colz");

  TH1D* IMpippim_wK0_woSid_won_mc = (TH1D*)MMnmiss_IMpippim_wK0_woSid_won_mc->ProjectionX("IMpippim_wK0_woSid_won_mc");
  IMpippim_wK0_woSid_won_mc->SetLineColor(6);
  TCanvas *cIMpippim_wK0_woSid_won = new TCanvas("cIMpippim_wK0_woSid_won","cIMpippim_wK0_woSid_won");
  cIMpippim_wK0_woSid_won->cd();
  TH1D* IMpippim_wK0_woSid_won[6];
  for(int i=0; i<6; i++) IMpippim_wK0_woSid_won[i] = (TH1D*)MMnmiss_IMpippim_wK0_woSid_won[i]->ProjectionX(Form("IMpippim_wK0_woSid_won_%s",name[i]));
  IMpippim_wK0_woSid_won[0]->Draw("HE");
  IMpippim_wK0_woSid_won_mc->Draw("HEsame");
  for(int i=1; i<6; i++) {
    IMpippim_wK0_woSid_won[i]->SetLineColor(colordef[i]);
    //IMpippim_wK0_woSid_won[i]->Draw("HEsame");
  }

  TCanvas *cIMpippim_wK0_woSid_won_ratio = new TCanvas("cIMpippim_wK0_woSid_won_ratio","cIMpippim_wK0_woSid_won_ratio");
  cIMpippim_wK0_woSid_won_ratio->cd();
  TH1D* IMpippim_wK0_woSid_won_ratio = (TH1D*)IMpippim_wK0_woSid_won[0]->Clone("IMpippim_wK0_woSid_won_ratio");
  IMpippim_wK0_woSid_won_ratio->Divide(IMpippim_wK0_woSid_won_mc);
  IMpippim_wK0_woSid_won_ratio->SetTitle("Data/MC");
  IMpippim_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,5);
  IMpippim_wK0_woSid_won_ratio->Draw("HE");

  TF1 *evalf_IMpippim_wK0 = new TF1("evalf_IMpippim_wK0","pol6",0,1);
  IMpippim_wK0_woSid_won_ratio->Fit("evalf_IMpippim_wK0","","",0.28,0.97);

  //q vs IMnpipi w/o K0 w/o (Sid & n);
  TH2D* q_IMnpipi_wK0_woSid_won_mc = (TH2D*)q_IMnpipi_wK0_woSid_won[1]->Clone("q_IMnpipi_wK0_woSid_won_mc");
  for(int i=2; i<6; i++)q_IMnpipi_wK0_woSid_won_mc->Add(q_IMnpipi_wK0_woSid_won[i]);

  TCanvas *cq_IMnpipi_wK0_woSid_won = new TCanvas("cq_IMnpipi_wK0_woSid_won","q_IMnpipi_wK0_woSid_won",1200,800);
  cq_IMnpipi_wK0_woSid_won->Divide(2,1);
  cq_IMnpipi_wK0_woSid_won->cd(1);
  q_IMnpipi_wK0_woSid_won[0]->SetTitle("#splitline{q_IMnpipi_wK0_woSid_won}{  Real data}");
  q_IMnpipi_wK0_woSid_won[0]->Draw("colz");
  cq_IMnpipi_wK0_woSid_won->cd(2);
  q_IMnpipi_wK0_woSid_won_mc->SetTitle("#splitline{q_IMnpipi_wK0_woSid_won}{  MC sum}");
  q_IMnpipi_wK0_woSid_won_mc->Draw("colz");

  TH1D* IMnpipi_wK0_woSid_won_mc = (TH1D*)q_IMnpipi_wK0_woSid_won_mc->ProjectionX("IMnpipi_wK0_woSid_won_mc");
  IMnpipi_wK0_woSid_won_mc->SetLineColor(6);
  TCanvas *cIMnpipi_wK0_woSid_won = new TCanvas("cIMnpipi_wK0_woSid_won","cIMnpipi_wK0_woSid_won");
  cIMnpipi_wK0_woSid_won->cd();
  TH1D* IMnpipi_wK0_woSid_won[6];
  for(int i=0; i<6; i++) IMnpipi_wK0_woSid_won[i] = (TH1D*)q_IMnpipi_wK0_woSid_won[i]->ProjectionX(Form("IMnpipi_wK0_woSid_won_%s",name[i]));
  IMnpipi_wK0_woSid_won[0]->Draw("HE");
  IMnpipi_wK0_woSid_won_mc->Draw("HEsame");
  for(int i=1; i<6; i++) {
    IMnpipi_wK0_woSid_won[i]->SetLineColor(colordef[i]);
    //IMnpipi_wK0_woSid_won[i]->Draw("HEsame");
  }

  TCanvas *cIMnpipi_wK0_woSid_won_ratio = new TCanvas("cIMnpipi_wK0_woSid_won_ratio","cIMnpipi_wK0_woSid_won_ratio");
  cIMnpipi_wK0_woSid_won_ratio->cd();
  TH1D* IMnpipi_wK0_woSid_won_ratio = (TH1D*) IMnpipi_wK0_woSid_won[0]->Clone("IMnpipi_wK0_woSid_won_ratio");
  IMnpipi_wK0_woSid_won_ratio->Divide(IMnpipi_wK0_woSid_won_mc);
  IMnpipi_wK0_woSid_won_ratio->SetTitle("Data/MC");
  IMnpipi_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  IMnpipi_wK0_woSid_won_ratio->Draw("HE");

  TF1 *evalf_IMnpipi_wK0 = new TF1("evalf_IMnpipi_wK0","pol4",1.0,2.00);
  IMnpipi_wK0_woSid_won_ratio->Fit("evalf_IMnpipi_wK0","","",1.0,2.00);


  TH1D* q_wK0_woSid_won_mc = (TH1D*)q_IMnpipi_wK0_woSid_won_mc->ProjectionY("q_wK0_woSid_won_mc");
  q_wK0_woSid_won_mc->SetLineColor(6);
  TCanvas *cq_wK0_woSid_won = new TCanvas("cq_wK0_woSid_won","cq_wK0_woSid_won");
  cq_wK0_woSid_won->cd();
  TH1D* q_wK0_woSid_won[6];
  for(int i=0; i<6; i++) q_wK0_woSid_won[i] = (TH1D*)q_IMnpipi_wK0_woSid_won[i]->ProjectionY(Form("q_wK0_woSid_won_%s",name[i]));
  q_wK0_woSid_won[0]->Draw("HE");
  q_wK0_woSid_won_mc->Draw("HEsame");
  for(int i=1; i<6; i++) {
    q_wK0_woSid_won[i]->SetLineColor(colordef[i]);
    //q_wK0_woSid_won[i]->Draw("HEsame");
  }

  //missing mass vs Mom(pi+pi-) w/o K0 w/o (Sid & n)
  TH2D* MMnmiss_Mompippim_wK0_woSid_won_mc = (TH2D*)MMnmiss_Mompippim_wK0_woSid_won[1]->Clone("MMnmiss_Mompippim_wK0_woSid_won_mc");
  for(int i=2; i<6; i++)MMnmiss_Mompippim_wK0_woSid_won_mc->Add(MMnmiss_Mompippim_wK0_woSid_won[i]);

  TCanvas *cMMnmiss_Mompippim_wK0_woSid_won = new TCanvas("cMMnmiss_Mompippim_wK0_woSid_won","cMMnmiss_Mompippim_wK0_woSid_won",1200,800);
  cMMnmiss_Mompippim_wK0_woSid_won->Divide(2,1);
  cMMnmiss_Mompippim_wK0_woSid_won->cd(1);
  MMnmiss_Mompippim_wK0_woSid_won[0]->SetTitle("#splitline{MMnmiss_Mompippim_wK0_woSid_won}{  Real data}");
  MMnmiss_Mompippim_wK0_woSid_won[0]->Draw("colz");
  cMMnmiss_Mompippim_wK0_woSid_won->cd(2);
  MMnmiss_Mompippim_wK0_woSid_won_mc->SetTitle("#splitline{MMnmiss_Mompippim_wK0_woSid_won}{  MC sum}");
  MMnmiss_Mompippim_wK0_woSid_won_mc->Draw("colz");

  TH1D* Mompippim_wK0_woSid_won_mc = (TH1D*)MMnmiss_Mompippim_wK0_woSid_won_mc->ProjectionX("Mompippim_wK0_woSid_won_mc");
  Mompippim_wK0_woSid_won_mc->SetLineColor(6);
  TCanvas *cMompippim_wK0_woSid_won = new TCanvas("cMompippim_wK0_woSid_won","cMompippim_wK0_woSid_won");
  cMompippim_wK0_woSid_won->cd();
  TH1D* Mompippim_wK0_woSid_won[6];
  for(int i=0; i<6; i++) Mompippim_wK0_woSid_won[i] = (TH1D*)MMnmiss_Mompippim_wK0_woSid_won[i]->ProjectionX(Form("MMnmiss_Mompippim_wK0_woSid_won_%s",name[i]));
  Mompippim_wK0_woSid_won[0]->Draw("HE");
  //Mompippim_wK0_woSid_won_mc->Draw("HEsame");
  for(int i=1; i<6; i++) {
    Mompippim_wK0_woSid_won[i]->SetLineColor(colordef[i]);
    //Mompippim_wK0_woSid_won[i]->Draw("HEsame");
  }

  TCanvas *cMompippim_wK0_woSid_won_ratio = new TCanvas("cMompippim_wK0_woSid_won_ratio","cMompippim_wK0_woSid_won_ratio");
  TH1D* Mompippim_wK0_woSid_won_ratio = (TH1D*)Mompippim_wK0_woSid_won[0]->Clone("Mompippim_wK0_woSid_won_ratio");
  Mompippim_wK0_woSid_won_ratio->Divide(Mompippim_wK0_woSid_won_mc);
  Mompippim_wK0_woSid_won_ratio->SetTitle("Data/MC");
  Mompippim_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  Mompippim_wK0_woSid_won_ratio->Draw("HE");

  TF1 *evalf_Mompippim_wK0 = new TF1("evalf_Mompippim_wK0","pol3",0,1);
  Mompippim_wK0_woSid_won_ratio->Fit("evalf_Mompippim_wK0","","",0,0.97);


  //pipmom
  TH2D* pipmom_MMnmiss_wK0_woSid_won_mc = (TH2D*)pipmom_MMnmiss_wK0_woSid_won[1]->Clone("pipmom_MMnmiss_wK0_woSid_won_mc");
  for(int i=2; i<6; i++) pipmom_MMnmiss_wK0_woSid_won_mc->Add(pipmom_MMnmiss_wK0_woSid_won[i]);

  TH1D* pipmom_wK0_woSid_won_mc = (TH1D*)pipmom_MMnmiss_wK0_woSid_won_mc->ProjectionY("pipmom_wK0_woSid_won_mc");
  TH1D* pipmom_wK0_woSid_won[6];
  for(int i=0; i<6; i++) {
    pipmom_wK0_woSid_won[i] = (TH1D*)pipmom_MMnmiss_wK0_woSid_won[i]->ProjectionY(Form("pipmom_wK0_woSid_won_%s",name[i]));
    pipmom_wK0_woSid_won[i]->SetLineColor(colordef[i]);
  }
  TCanvas *cpipmom_wK0_woSid_won = new TCanvas("cpipmom_wK0_woSid_won","cpipmom_wK0_woSid_won");
  cpipmom_wK0_woSid_won->cd();
  pipmom_wK0_woSid_won[0]->Draw("HE");
  pipmom_wK0_woSid_won_mc->SetLineColor(6);
  pipmom_wK0_woSid_won_mc->Draw("HEsame");
  //for(int i=1;i<6;i++)pipmom_wK0_woSid_won[i]->Draw("same");

  TCanvas *cpipmom_wK0_woSid_won_ratio = new TCanvas("cpipmom_wK0_woSid_won_ratio","cpipmom_wK0_woSid_won_ratio");
  cpipmom_wK0_woSid_won_ratio->cd();
  TH1D* pipmom_wK0_woSid_won_ratio = (TH1D*)pipmom_wK0_woSid_won[0]->Clone("pipmom_wK0_woSid_won_ratio");
  pipmom_wK0_woSid_won_ratio->Divide(pipmom_wK0_woSid_won_mc);
  pipmom_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(-1,6);
  pipmom_wK0_woSid_won_ratio->SetTitle("Data/MC");
  pipmom_wK0_woSid_won_ratio->Draw("HEsame");

  TF1 *evalf_pipmom_wK0 = new TF1("evalf_pipmom_wK0","pol8",0.08,0.7);
  pipmom_wK0_woSid_won_ratio->Fit("evalf_pipmom_wK0","","",0.08,0.7);
  
  //
  TH2D* pipmom_pimmom_wK0_woSid_won_mc = (TH2D*)pipmom_pimmom_wK0_woSid_won[1]->Clone("pipmom_pimmom_wK0_woSid_won_mc");
  for(int i=1;i<6;i++)pipmom_pimmom_wK0_woSid_won_mc->Add(pipmom_pimmom_wK0_woSid_won[i]);
  TCanvas *cpipmom_pimmom_wK0_woSid_won = new TCanvas("cpipmom_pimmom_wK0_woSid_won","cpipmom_pimmom_wK0_woSid_won");
  cpipmom_pimmom_wK0_woSid_won->Divide(2,1);
  cpipmom_pimmom_wK0_woSid_won->cd(1);
  pipmom_pimmom_wK0_woSid_won[0]->Draw("colz");
  cpipmom_pimmom_wK0_woSid_won->cd(2);
  pipmom_pimmom_wK0_woSid_won_mc->Draw("colz");
  
  
  //pimmom
  TH2D* pimmom_MMnmiss_wK0_woSid_won_mc = (TH2D*)pimmom_MMnmiss_wK0_woSid_won[1]->Clone("pimmom_MMnmiss_wK0_woSid_won_mc");
  for(int i=2; i<6; i++) pimmom_MMnmiss_wK0_woSid_won_mc->Add(pimmom_MMnmiss_wK0_woSid_won[i]);

  TH1D* pimmom_wK0_woSid_won_mc = pimmom_MMnmiss_wK0_woSid_won_mc->ProjectionY("pimmom_wK0_woSid_won_mc");
  TH1D* pimmom_wK0_woSid_won[6];
  for(int i=0; i<6; i++) {
    pimmom_wK0_woSid_won[i] = (TH1D*)pimmom_MMnmiss_wK0_woSid_won[i]->ProjectionY(Form("pimmom_wK0_woSid_won_%s",name[i]));
    pimmom_wK0_woSid_won[i]->SetLineColor(colordef[i]);
  }
  TCanvas *cpimmom_wK0_woSid_won = new TCanvas("cpimmom_wK0_woSid_won","cpimmom_wK0_woSid_won");
  cpimmom_wK0_woSid_won->cd();
  pimmom_wK0_woSid_won[0]->Draw("HE");
  pimmom_wK0_woSid_won_mc->SetLineColor(6);
  pimmom_wK0_woSid_won_mc->Draw("HEsame");
  //for(int i=1;i<6;i++)pimmom_wK0_woSid_won[i]->Draw("same");

  TCanvas *cpimmom_wK0_woSid_won_ratio = new TCanvas("cpimmom_wK0_woSid_won_ratio","cpimmom_wK0_woSid_won_ratio");
  cpimmom_wK0_woSid_won_ratio->cd();
  TH1D* pimmom_wK0_woSid_won_ratio = (TH1D*)pimmom_wK0_woSid_won[0]->Clone("pimmom_wK0_woSid_won_ratio");
  pimmom_wK0_woSid_won_ratio->Divide(pimmom_wK0_woSid_won_mc);
  pimmom_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(-1,6);
  pimmom_wK0_woSid_won_ratio->SetTitle("Data/MC");
  pimmom_wK0_woSid_won_ratio->Draw("HEsame");

  TF1 *evalf_pimmom_wK0 = new TF1("evalf_pimmom_wK0","pol8",0.06,0.73);
  pimmom_wK0_woSid_won_ratio->Fit("evalf_pimmom_wK0","","",0.06,0.73);
  
  //nCDSmom
  TH2D* nmom_MMnmiss_wK0_woSid_won_mc = (TH2D*)nmom_MMnmiss_wK0_woSid_won[1]->Clone("nmom_MMnmiss_wK0_woSid_won_mc");
  for(int i=2; i<6; i++) nmom_MMnmiss_wK0_woSid_won_mc->Add(nmom_MMnmiss_wK0_woSid_won[i]);

  TH1D* nmom_wK0_woSid_won_mc = (TH1D*)nmom_MMnmiss_wK0_woSid_won_mc->ProjectionY("nmom_wK0_woSid_won_mc");
  TH1D* nmom_wK0_woSid_won[6];
  for(int i=0; i<6; i++) {
    nmom_wK0_woSid_won[i] = (TH1D*)nmom_MMnmiss_wK0_woSid_won[i]->ProjectionY(Form("nmom_wK0_woSid_won_%s",name[i]));
    nmom_wK0_woSid_won[i]->SetLineColor(colordef[i]);
  }
  TCanvas *cnmom_wK0_woSid_won = new TCanvas("cnmom_wK0_woSid_won","cnmom_wK0_woSid_won");
  cnmom_wK0_woSid_won->cd();
  nmom_wK0_woSid_won[0]->Draw("HE");
  nmom_wK0_woSid_won_mc->SetLineColor(6);
  nmom_wK0_woSid_won_mc->Draw("HEsame");
  //for(int i=1;i<6;i++)nmom_wK0_woSid_won[i]->Draw("same");

  TCanvas *cnmom_wK0_woSid_won_ratio = new TCanvas("cnmom_wK0_woSid_won_ratio","cnmom_wK0_woSid_won_ratio");
  cnmom_wK0_woSid_won_ratio->cd();
  TH1D* nmom_wK0_woSid_won_ratio = (TH1D*)nmom_wK0_woSid_won[0]->Clone("nmom_wK0_woSid_won_ratio");
  nmom_wK0_woSid_won_ratio->Divide(nmom_wK0_woSid_won_mc);
  nmom_wK0_woSid_won_ratio->SetTitle("Data/MC");
  nmom_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(-1,4);
  nmom_wK0_woSid_won_ratio->Draw("HEsame");

  TF1 *evalf_nmom_wK0 = new TF1("evalf_nmom_wK0","pol8",0.14,1);
  nmom_wK0_woSid_won_ratio->Fit("evalf_nmom_wK0","","",0.14,1);
  


  //weighting to momentum vector 
  TH2D* nmom_cosn_wK0_woSid_won_mc = (TH2D*)nmom_cosn_wK0_woSid_won[1]->Clone("nmom_cosn_wK0_woSid_won_mc");
  for(int i=2;i<6;i++) nmom_cosn_wK0_woSid_won_mc->Add(nmom_cosn_wK0_woSid_won[i]);
  TCanvas *cnmom_cosn_wK0_woSid_won = new TCanvas("cnmom_cosn_wK0_woSid_won","cnmom_cosn_wK0_woSid_won",1200,800);
  cnmom_cosn_wK0_woSid_won->Divide(2,1);
  cnmom_cosn_wK0_woSid_won->cd(1);
  nmom_cosn_wK0_woSid_won[0]->SetTitle("#splitline{nmom_cosn_wK0_woSid_won}{Real data}");
  nmom_cosn_wK0_woSid_won[0]->Draw("colz");
  cnmom_cosn_wK0_woSid_won->cd(2);
  nmom_cosn_wK0_woSid_won_mc->SetTitle("#splitline{nmom_cosn_wK0_woSid_won}{MC}");
  nmom_cosn_wK0_woSid_won_mc->Draw("colz");

  TCanvas *ccosn_wK0_woSid_won = new TCanvas("ccosn_wK0_woSid_won","ccosn_wK0_woSid_won");
  ccosn_wK0_woSid_won->cd();
  TH1D* cosn_wK0_woSid_won[6];
  for(int i=0;i<6;i++) cosn_wK0_woSid_won[i] = (TH1D*)nmom_cosn_wK0_woSid_won[i]->ProjectionX(Form("cosn_wK0_woSid_won_%s",name[i]));
  TH1D* cosn_wK0_woSid_won_mc = (TH1D*)nmom_cosn_wK0_woSid_won_mc->ProjectionX("cosn_wK0_woSid_won_mc");
  cosn_wK0_woSid_won_mc->SetLineColor(6);
  cosn_wK0_woSid_won[0]->Draw("HE");
  cosn_wK0_woSid_won_mc->Draw("HEsame");
   

  TCanvas *ccosn_wK0_woSid_won_ratio = new TCanvas("ccosn_wK0_woSid_won_ratio","ccosn_wK0_woSid_won_ratio");
  ccosn_wK0_woSid_won_ratio->cd();
  TH1D* cosn_wK0_woSid_won_ratio = (TH1D*)cosn_wK0_woSid_won[0]->Clone("cosn_wK0_woSid_won_ratio");
  cosn_wK0_woSid_won_ratio->Divide(cosn_wK0_woSid_won_mc);
  cosn_wK0_woSid_won_ratio->SetTitle("Data/MC");
  cosn_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,4);
  cosn_wK0_woSid_won_ratio->Draw("HE");

  //TF1 *evalf_cosn_wK0 = new TF1("evalf_cosn_wK0","pol8",-1,1);
  //cosn_wK0_woSid_won_ratio->Fit(evalf_cosn_wK0,"","",-1,0.2);
  
  TF1 *evalf_cosn_wK0 = new TF1("evalf_cosn_wK0",func_cosn,-1.00,0.3,6);
  evalf_cosn_wK0->SetParameters(param_cosn);
  evalf_cosn_wK0->SetLineColor(3);
  cosn_wK0_woSid_won_ratio->Fit(evalf_cosn,"","",-1,0.2);
  


  //
  TH2D* nmom_cospip_wK0_woSid_won_mc = (TH2D*)nmom_cospip_wK0_woSid_won[1]->Clone("nmom_cospip_wK0_woSid_won_mc");
  for(int i=2;i<6;i++)nmom_cospip_wK0_woSid_won_mc->Add(nmom_cospip_wK0_woSid_won[i]);
  TCanvas *cnmom_cospip_wK0_woSid_won = new TCanvas("cnmom_cospip_wK0_woSid_won","cnmom_cospip_wK0_woSid_won",1200,800);
  cnmom_cospip_wK0_woSid_won->Divide(2,1);
  cnmom_cospip_wK0_woSid_won->cd(1);
  nmom_cospip_wK0_woSid_won[0]->SetTitle("#splitline{nmom_cospip_wK0_woSid_won}{Real data}");
  nmom_cospip_wK0_woSid_won[0]->Draw("colz");
  cnmom_cospip_wK0_woSid_won->cd(2);
  nmom_cospip_wK0_woSid_won_mc->SetTitle("#splitline{nmom_cospip_wK0_woSid_won}{MC}");
  nmom_cospip_wK0_woSid_won_mc->Draw("colz");

  TCanvas *ccospip_wK0_woSid_won = new TCanvas("ccospip_wK0_woSid_won","ccospip_wK0_woSid_won");
  ccospip_wK0_woSid_won->cd();
  TH1D* cospip_wK0_woSid_won[6];
  for(int i=0;i<6;i++) cospip_wK0_woSid_won[i] = (TH1D*)nmom_cospip_wK0_woSid_won[i]->ProjectionX(Form("cospip_wK0_woSid_won_%s",name[i]));
  TH1D* cospip_wK0_woSid_won_mc = (TH1D*)nmom_cospip_wK0_woSid_won_mc->ProjectionX("cospip_wK0_woSid_won_mc");
  cospip_wK0_woSid_won_mc->SetLineColor(6);
  cospip_wK0_woSid_won[0]->Draw("HE");
  cospip_wK0_woSid_won_mc->Draw("HEsame");

  TCanvas *ccospip_wK0_woSid_won_ratio = new TCanvas("ccospip_wK0_woSid_won_ratio","ccospip_wK0_woSid_won_ratio");
  ccospip_wK0_woSid_won_ratio->cd();
  TH1D* cospip_wK0_woSid_won_ratio = (TH1D*)cospip_wK0_woSid_won[0]->Clone("cospip_wK0_woSid_won_ratio");
  cospip_wK0_woSid_won_ratio->Divide(cospip_wK0_woSid_won_mc);
  cospip_wK0_woSid_won_ratio->SetTitle("Data/MC");
  cospip_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,4);
  cospip_wK0_woSid_won_ratio->Draw("HE");
  
  TF1 *evalf_cospip_wK0 = new TF1("evalf_cospip_wK0","pol5",-1.0,0.6); 
  cospip_woK0_woSidn_ratio->Fit("evalf_cospip_wK0","","",-1.0,0.6);



  //
  TH2D* nmom_cospim_wK0_woSid_won_mc = (TH2D*)nmom_cospim_wK0_woSid_won[1]->Clone("nmom_cospim_wK0_woSid_won_mc");
  for(int i=2;i<6;i++)nmom_cospim_wK0_woSid_won_mc->Add(nmom_cospim_wK0_woSid_won[i]);
  TCanvas *cnmom_cospim_wK0_woSid_won = new TCanvas("cnmom_cospim_wK0_woSid_won","cnmom_cospim_wK0_woSid_won",1200,800);
  cnmom_cospim_wK0_woSid_won->Divide(2,1);
  cnmom_cospim_wK0_woSid_won->cd(1);
  nmom_cospim_wK0_woSid_won[0]->SetTitle("#splitline{nmom_cospim_wK0_woSid_won}{Real data}");
  nmom_cospim_wK0_woSid_won[0]->Draw("colz");
  cnmom_cospim_wK0_woSid_won->cd(2);
  nmom_cospim_wK0_woSid_won_mc->SetTitle("#splitline{nmom_cospim_wK0_woSid_won}{MC}");
  nmom_cospim_wK0_woSid_won_mc->Draw("colz");

  TCanvas *ccospim_wK0_woSid_won = new TCanvas("ccospim_wK0_woSid_won","ccospim_wK0_woSid_won");
  ccospim_wK0_woSid_won->cd();
  TH1D* cospim_wK0_woSid_won[6];
  for(int i=0;i<6;i++) cospim_wK0_woSid_won[i] = (TH1D*)nmom_cospim_wK0_woSid_won[i]->ProjectionX(Form("cospim_wK0_woSid_won_%s",name[i]));
  TH1D* cospim_wK0_woSid_won_mc = (TH1D*)nmom_cospim_wK0_woSid_won_mc->ProjectionX("cospim_wK0_woSid_won_mc");
  cospim_wK0_woSid_won_mc->SetLineColor(6);
  cospim_wK0_woSid_won[0]->Draw("HE");
  cospim_wK0_woSid_won_mc->Draw("HEsame");

  TCanvas *ccospim_wK0_woSid_won_ratio = new TCanvas("ccospim_wK0_woSid_won_ratio","ccospim_wK0_woSid_won_ratio");
  ccospim_wK0_woSid_won_ratio->cd();
  TH1D* cospim_wK0_woSid_won_ratio = (TH1D*)cospim_wK0_woSid_won[0]->Clone("cospim_wK0_woSid_won_ratio");
  cospim_wK0_woSid_won_ratio->Divide(cospim_wK0_woSid_won_mc);
  cospim_wK0_woSid_won_ratio->SetTitle("Data/MC");
  cospim_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,4);
  cospim_wK0_woSid_won_ratio->Draw("HE");
   
  TF1 *evalf_cospim_wK0 = new TF1("evalf_cospim_wK0","pol8",-0.92,0.55);
  cospim_wK0_woSid_won_ratio->Fit(evalf_cospim_wK0,"","",-0.92,0.55);

  //
  TH2D* nmom_phinpip_wK0_woSid_won_mc = (TH2D*)nmom_phinpip_wK0_woSid_won[1]->Clone("nmom_phinpip_wK0_woSid_won_mc");
  for(int i=2;i<6;i++)nmom_phinpip_wK0_woSid_won_mc->Add(nmom_phinpip_wK0_woSid_won[i]);
  TCanvas *cnmom_phinpip_wK0_woSid_won = new TCanvas("cnmom_phinpip_wK0_woSid_won","cnmom_phinpip_wK0_woSid_won",1200,800);
  cnmom_phinpip_wK0_woSid_won->Divide(2,1);
  cnmom_phinpip_wK0_woSid_won->cd(1);
  nmom_phinpip_wK0_woSid_won[0]->SetTitle("#splitline{nmom_phinpip_wK0_woSid_won}{Real data}");
  nmom_phinpip_wK0_woSid_won[0]->Draw("colz");
  cnmom_phinpip_wK0_woSid_won->cd(2);
  nmom_phinpip_wK0_woSid_won_mc->SetTitle("#splitline{nmom_phinpip_wK0_woSid_won}{MC}");
  nmom_phinpip_wK0_woSid_won_mc->Draw("colz");

  TCanvas *cphinpip_wK0_woSid_won = new TCanvas("cphinpip_wK0_woSid_won","cphinpip_wK0_woSid_won");
  cphinpip_wK0_woSid_won->cd();
  TH1D* phinpip_wK0_woSid_won[6];
  for(int i=0;i<6;i++) phinpip_wK0_woSid_won[i] = (TH1D*)nmom_phinpip_wK0_woSid_won[i]->ProjectionX(Form("phinpip_wK0_woSid_won_%s",name[i]));
  TH1D* phinpip_wK0_woSid_won_mc = (TH1D*)nmom_phinpip_wK0_woSid_won_mc->ProjectionX("phinpip_wK0_woSid_won_mc");
  phinpip_wK0_woSid_won_mc->SetLineColor(6);
  phinpip_wK0_woSid_won[0]->Draw("HE");
  phinpip_wK0_woSid_won_mc->Draw("HEsame");

  TCanvas *cphinpip_wK0_woSid_won_ratio = new TCanvas("cphinpip_wK0_woSid_won_ratio","cphinpip_wK0_woSid_won_ratio");
  cphinpip_wK0_woSid_won_ratio->cd();
  TH1D* phinpip_wK0_woSid_won_ratio = (TH1D*)phinpip_wK0_woSid_won[0]->Clone("phinpip_wK0_woSid_won_ratio");
  phinpip_wK0_woSid_won_ratio->Divide(phinpip_wK0_woSid_won_mc);
  phinpip_wK0_woSid_won_ratio->SetTitle("Data/MC");
  phinpip_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,4);
  phinpip_wK0_woSid_won_ratio->Draw("HE");
   
  //
  TH2D* nmom_phinpim_wK0_woSid_won_mc = (TH2D*)nmom_phinpim_wK0_woSid_won[1]->Clone("nmom_phinpim_wK0_woSid_won_mc");
  for(int i=2;i<6;i++)nmom_phinpim_wK0_woSid_won_mc->Add(nmom_phinpim_wK0_woSid_won[i]);
  TCanvas *cnmom_phinpim_wK0_woSid_won = new TCanvas("cnmom_phinpim_wK0_woSid_won","cnmom_phinpim_wK0_woSid_won",1200,800);
  cnmom_phinpim_wK0_woSid_won->Divide(2,1);
  cnmom_phinpim_wK0_woSid_won->cd(1);
  nmom_phinpim_wK0_woSid_won[0]->SetTitle("#splitline{nmom_phinpim_wK0_woSid_won}{Real data}");
  nmom_phinpim_wK0_woSid_won[0]->Draw("colz");
  cnmom_phinpim_wK0_woSid_won->cd(2);
  nmom_phinpim_wK0_woSid_won_mc->SetTitle("#splitline{nmom_phinpim_wK0_woSid_won}{MC}");
  nmom_phinpim_wK0_woSid_won_mc->Draw("colz");

  TCanvas *cphinpim_wK0_woSid_won = new TCanvas("cphinpim_wK0_woSid_won","cphinpim_wK0_woSid_won");
  cphinpim_wK0_woSid_won->cd();
  TH1D* phinpim_wK0_woSid_won[6];
  for(int i=0;i<6;i++) phinpim_wK0_woSid_won[i] = (TH1D*)nmom_phinpim_wK0_woSid_won[i]->ProjectionX(Form("phinpim_wK0_woSid_won_%s",name[i]));
  TH1D* phinpim_wK0_woSid_won_mc = (TH1D*)nmom_phinpim_wK0_woSid_won_mc->ProjectionX("phinpim_wK0_woSid_won_mc");
  phinpim_wK0_woSid_won_mc->SetLineColor(6);
  phinpim_wK0_woSid_won[0]->Draw("HE");
  phinpim_wK0_woSid_won_mc->Draw("HEsame");

  TCanvas *cphinpim_wK0_woSid_won_ratio = new TCanvas("cphinpim_wK0_woSid_won_ratio","cphinpim_wK0_woSid_won_ratio");
  cphinpim_wK0_woSid_won_ratio->cd();
  TH1D* phinpim_wK0_woSid_won_ratio = (TH1D*)phinpim_wK0_woSid_won[0]->Clone("phinpim_wK0_woSid_won_ratio");
  phinpim_wK0_woSid_won_ratio->Divide(phinpim_wK0_woSid_won_mc);
  phinpim_wK0_woSid_won_ratio->SetTitle("Data/MC");
  phinpim_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,4);
  phinpim_wK0_woSid_won_ratio->Draw("HE");



  //Signal check
  TH2D* q_IMnpipi_wSid_n_Sp_mc = (TH2D*)q_IMnpipi_wSid_n_Sp[1]->Clone();
  for(int i=2; i<6; i++)q_IMnpipi_wSid_n_Sp_mc->Add(q_IMnpipi_wSid_n_Sp[i]);

  TCanvas *cq_IMnpipi_wSid_n_Sp = new TCanvas("cq_IMnpipi_wSid_n_Sp","cq_IMnpipi_wSid_n_Sp");
  cq_IMnpipi_wSid_n_Sp->Divide(2,1);
  cq_IMnpipi_wSid_n_Sp->cd(1);
  q_IMnpipi_wSid_n_Sp[0]->SetMinimum(1);
  q_IMnpipi_wSid_n_Sp[0]->Draw("colz");
  cq_IMnpipi_wSid_n_Sp->cd(2);
  q_IMnpipi_wSid_n_Sp_mc->SetMinimum(1);
  q_IMnpipi_wSid_n_Sp_mc->SetMaximum(q_IMnpipi_wSid_n_Sp[0]->GetMaximum());
  q_IMnpipi_wSid_n_Sp_mc->Draw("colz");

  TCanvas *cIMnpipi_wSid_n_Sp_0 = new TCanvas("cIMnpipi_wSid_n_Sp_0","cIMnpipi_wSid_n_Sp_0");
  const int q350bin = q_IMnpipi_wSid_n_Sp[0]->GetYaxis()->FindBin(0.35);
  TH1D* IMnpipi_wSid_n_Sp_0[6];
  for(int i=0; i<6; i++)IMnpipi_wSid_n_Sp_0[i] = q_IMnpipi_wSid_n_Sp[i]->ProjectionX(Form("IMnpipi_wSid_n_Sp_0_%s",name[i]),0,q350bin);
  IMnpipi_wSid_n_Sp_0[0]->Draw("HE");
  TH1D* IMnpipi_wSid_n_Sp_mc_0 = q_IMnpipi_wSid_n_Sp_mc->ProjectionX("IMnpipi_wSid_n_Sp_mc_0",0,q350bin-1);
  IMnpipi_wSid_n_Sp_mc_0->SetLineColor(6);
  IMnpipi_wSid_n_Sp_mc_0->Draw("HEsame");

  TCanvas *cIMnpipi_wSid_n_Sp_350 = new TCanvas("cIMnpipi_wSid_n_Sp_350","cIMnpipi_wSid_n_Sp_350");
  TH1D* IMnpipi_wSid_n_Sp_350[6];
  for(int i=0; i<6; i++)IMnpipi_wSid_n_Sp_350[i] = q_IMnpipi_wSid_n_Sp[i]->ProjectionX(Form("IMnpipi_wSid_n_Sp_350_%s",name[i]),q350bin,100);
  IMnpipi_wSid_n_Sp_350[0]->Draw("HE");
  TH1D* IMnpipi_wSid_n_Sp_mc_350 = q_IMnpipi_wSid_n_Sp_mc->ProjectionX("IMnpipi_wSid_n_Sp_mc_350",q350bin,100);
  IMnpipi_wSid_n_Sp_mc_350->SetLineColor(6);
  IMnpipi_wSid_n_Sp_mc_350->Draw("HEsame");


  TH2D* q_IMnpipi_wSid_n_Sm_mc = (TH2D*)q_IMnpipi_wSid_n_Sm[1]->Clone();
  for(int i=2; i<6; i++)q_IMnpipi_wSid_n_Sm_mc->Add(q_IMnpipi_wSid_n_Sm[i]);

  TCanvas *cq_IMnpipi_wSid_n_Sm = new TCanvas("cq_IMnpipi_wSid_n_Sm","cq_IMnpipi_wSid_n_Sm");
  cq_IMnpipi_wSid_n_Sm->Divide(2,1);
  cq_IMnpipi_wSid_n_Sm->cd(1);
  q_IMnpipi_wSid_n_Sm[0]->SetMinimum(1);
  q_IMnpipi_wSid_n_Sm[0]->Draw("colz");
  cq_IMnpipi_wSid_n_Sm->cd(2);
  q_IMnpipi_wSid_n_Sm_mc->SetMinimum(1);
  q_IMnpipi_wSid_n_Sm_mc->SetMaximum(q_IMnpipi_wSid_n_Sm[0]->GetMaximum());
  q_IMnpipi_wSid_n_Sm_mc->Draw("colz");

  TCanvas *cIMnpipi_wSid_n_Sm_0 = new TCanvas("cIMnpipi_wSid_n_Sm_0","cIMnpipi_wSid_n_Sm_0");
  TH1D* IMnpipi_wSid_n_Sm_0[6];
  for(int i=0; i<6; i++)IMnpipi_wSid_n_Sm_0[i] = q_IMnpipi_wSid_n_Sm[i]->ProjectionX(Form("IMnpipi_wSid_n_Sm_0_%s",name[i]),0,q350bin);
  IMnpipi_wSid_n_Sm_0[0]->Draw("HE");
  TH1D* IMnpipi_wSid_n_Sm_mc_0 = q_IMnpipi_wSid_n_Sm_mc->ProjectionX("IMnpipi_wSid_n_Sm_mc_0",0,q350bin-1);
  IMnpipi_wSid_n_Sm_mc_0->SetLineColor(6);
  IMnpipi_wSid_n_Sm_mc_0->Draw("HEsame");

  TCanvas *cIMnpipi_wSid_n_Sm_350 = new TCanvas("cIMnpipi_wSid_n_Sm_350","cIMnpipi_wSid_n_Sm_350");
  TH1D* IMnpipi_wSid_n_Sm_350[6];
  for(int i=0; i<6; i++)IMnpipi_wSid_n_Sm_350[i] = q_IMnpipi_wSid_n_Sm[i]->ProjectionX(Form("IMnpipi_wSid_n_Sm_350_%s",name[i]),q350bin,100);
  IMnpipi_wSid_n_Sm_350[0]->Draw("HE");
  TH1D* IMnpipi_wSid_n_Sm_mc_350 = q_IMnpipi_wSid_n_Sm_mc->ProjectionX("IMnpipi_wSid_n_Sm_mc_350",q350bin,100);
  IMnpipi_wSid_n_Sm_mc_350->SetLineColor(6);
  IMnpipi_wSid_n_Sm_mc_350->Draw("HEsame");

  TCanvas *cfunc = new TCanvas("cfunc","cfunc",1200,800);
  cfunc->Divide(4,2);
  cfunc->cd(1);
  evalf_MMnmiss->Draw();
  cfunc->cd(2);
  evalf_MMom->Draw();
  cfunc->cd(2);
  evalf_nmom->Draw();
  cfunc->cd(3);
  evalf_pipmom->Draw();
  cfunc->cd(4);
  evalf_pimmom->Draw();
  cfunc->cd(5);
  evalf_IMpippim->Draw();
  cfunc->cd(6);
  evalf_Mompippim->Draw();
  cfunc->cd(7);
  evalf_Mompippim->Draw();
  cfunc->cd(8);
  evalf_IMnpim->Draw();


  TCanvas *c = NULL;
  TSeqCollection *SCol = gROOT->GetListOfCanvases();
  int size = SCol->GetSize();
  TIter next(SCol);
  TString pdfname="plot_miss2D_out.pdf";
  for(int i=0; i<size; i++) {
    c= (TCanvas*)next();
    c->Draw();
    c->cd();
    if(i==0) c->Print(pdfname+"(",Form("Title:%s",c->GetTitle()));
    else if(i==size-1)c->Print(pdfname+")",Form("Title:%s",c->GetTitle()));
    else c->Print(pdfname,Form("Title:%s",c->GetTitle()));
  }

  TIter nexthist2(gDirectory->GetList());

  TFile *fout = new TFile("plot_miss2D_out.root","RECREATE");
  TObject *obj = nullptr;
  while( (obj = (TObject*)nexthist2())!=NULL ) {
    obj->Write();
  }

 // evalf_MMnmiss->Write();
 // evalf_MMom->Write();
  evalf_nmom->Write();
  evalf_pipmom->Write();
  evalf_pimmom->Write();
//  evalf_IMpippim->Write();
//  evalf_Mompippim->Write();
//  evalf_IMnpipi->Write();
//  evalf_IMnpip->Write();
//  evalf_IMnpim->Write();
  evalf_nmom_wK0->Write();
//  evalf_MMom_wK0->Write();
  evalf_pipmom_wK0->Write();
  evalf_pimmom_wK0->Write();
  evalf_cosn->Write();
  evalf_cosn_wK0->Write();
  evalf_cospip->Write();
  evalf_cospim->Write();
  evalf_cospip_wK0->Write();
  evalf_cospim_wK0->Write();
//  evalf_IMnpipi_wK0->Write();
//  evalf_IMnpip_wK0->Write();
//  evalf_IMnpim_wK0->Write();
  fout->Print();
  fout->cd();
  fout->Close();
}

