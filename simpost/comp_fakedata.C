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

/*
Double_t func_MMnmiss(Double_t *x,Double_t *par)
{
  if(x[0]<1.116) {
    return par[0]*exp(-0.5*pow((x[0]-par[1])/(par[2]+(x[0]<par[1])*par[3]*(x[0]-par[1])),2.0));
  } else if(1.116<=x[0] && x[0]<1.5) {
    return par[4]*exp(-0.5*pow((x[0]-par[5])/par[6],2.0));
  } else {
    return 0.;
  }
}*/

/*
Double_t func_MMnmiss(Double_t *x,Double_t *par)
{
  if( x[0]<1.12) {
    return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)
    +par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0);
  } else if(1.12<=x[0] && x[0]<1.5) {
    return par[6]+par[7]*x[0]+par[8]*pow(x[0],2.0)+par[9]*pow(x[0],3.0)
    +par[10]*pow(x[0],4.0)+par[11]*pow(x[0],5.0);
  } else {
    return 1.;
  }
}*/

Double_t func_MMnmiss(Double_t *x,Double_t *par)
{
   if(0.0<x[0] && x[0]<1.06){
     return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)
       +par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0)+par[6]*pow(x[0],6.0)
       +par[7]*pow(x[0],7.0)+par[8]*pow(x[0],8.0)+par[9]*pow(x[0],9.0);
   }else if(1.06<=x[0] && x[0]<1.115 ){
     return par[10]+par[11]*x[0]+par[12]*pow(x[0],2.0);
   }else if(1.115<x[0] && x[0]<1.5){
     return par[13]+par[14]*x[0]+par[15]*pow(x[0],2.0)+par[16]*pow(x[0],3.0)+par[17]*pow(x[0],4.0);
   }else{
     return 0;
   }
}




Double_t func_MMnmiss_corr(Double_t *x,Double_t *par)
{
  if(0.0 <x[0] && x[0]<0.86){
    return 1.0;
  }else if(0.86<=x[0] && x[0]<1.40){
    return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)
    +par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0)+par[6]*pow(x[0],6.0);
  }else{
    return 1.0;
  }
}

Double_t func_MMnmiss_corr2(Double_t *x,Double_t *par)
{
  if(0.0 <x[0] && x[0]<0.78){
    return 1.0;
  }else if(0.78<=x[0] && x[0]<1.50){
    return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)
    +par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0)+par[6]*pow(x[0],6.0);
  }else{
    return 1.0;
  }
}


Double_t func_MMnmiss_mul(Double_t *x,Double_t *par)
{
   return func_MMnmiss(&x[0],&par[0])*func_MMnmiss_corr(&x[0],&par[12]);
}


Double_t func_MMnmiss_wK0(Double_t *x,Double_t *par)
{
  if(x[0]<1.05){
    return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)
    +par[4]*pow(x[0],4.0);
  }else if(1.05<=x[0] && x[0]<1.40) {
    return par[5]+par[6]*x[0]+par[7]*pow(x[0],2.0)+par[8]*pow(x[0],3.0)
    +par[9]*pow(x[0],4.0);
  }else{ 
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


Double_t func_cosn_corr(Double_t *x,Double_t *par)
{
  if(-1.0<=x[0] && x[0]<0.90){
    return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)+par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0);
  }else{
    return 1.0;
  }
}





Double_t func_cospip(Double_t *x,Double_t *p)
{
  if(-0.92<x[0] && x[0]<-0.75){
    return TMath::Exp(p[0]+p[1]*x[0]);
  }else{
    return (p[2]+p[3]*x[0]+p[4]*TMath::Power(x[0],2)+p[5]*TMath::Power(x[0],3)) ;
  }
}

Double_t func_phinpip(Double_t *x,Double_t *par)
{
  if(x[0]<-0.4){
    return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0);
  }else{
    return par[4]+par[5]*x[0]+par[6]*pow(x[0],2.0)+par[7]*pow(x[0],3.0);
  }
}


Double_t func_phinpim(Double_t *x,Double_t *par)
{
  if(x[0]<0.5){
    return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)
    +par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0);
  }else{
    return par[6]+par[7]*x[0]+par[8]*pow(x[0],2.0)+par[9]*pow(x[0],3.0);
  }
}


Double_t func_phinpim_wK0(Double_t *x,Double_t *par)
{
  if(x[0]<0.3){
    return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)
    +par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0);
  }else{
    return par[6]+par[7]*x[0]+par[8]*pow(x[0],2.0)+par[9]*pow(x[0],3.0);
  }
}



void comp_fakedata()
{
  TH1::SetDefaultSumw2();
  //rdata,nSp,nSm,Lambda,S0pi+pi-,fake,mcsum
  const unsigned int colordef[7]= {1,2,3,4,5,28,6};
  const char name[][10]= {"rdata","Sigmap","Sigman","Lambda","S0","K0nn"};

  //real data
  TFile *filerdata = TFile::Open("../post/evanaIMpisigma_npippim_v196_out.root","READ");
  //Sp
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sp_rdata = (TH2D*) filerdata->Get("MMnmiss_IMnpipi_woK0_wSid_Sp");
  double  nrdata_Sp = MMnmiss_IMnpipi_woK0_wSid_Sp_rdata->GetEntries();
  std::cout << nrdata_Sp << std::endl;
  //Sm
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sm_rdata = (TH2D*) filerdata->Get("MMnmiss_IMnpipi_woK0_wSid_Sm");
  double nrdata_Sm = MMnmiss_IMnpipi_woK0_wSid_Sm_rdata->GetEntries();

  //fake
  TFile *filefake = TFile::Open("fakepippim_pippimn_out_sum.root","READ");
  //Sp
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sp_fake = (TH2D*) filefake->Get("MMnmiss_IMnpipi_woK0_wSid_Sp");
  double fake_Sp = MMnmiss_IMnpipi_woK0_wSid_Sp_fake->GetEntries();
  std::cout << fake_Sp << std::endl;
  //Sm
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sm_fake = (TH2D*) filefake->Get("MMnmiss_IMnpipi_woK0_wSid_Sm");
  double fake_Sm = MMnmiss_IMnpipi_woK0_wSid_Sm_fake->GetEntries();



  //superimpose 2d hists of miss mass vs IM(npi+) or IM(npi-)
  //fake
  TH2D *MMnmiss_IMnpip_woK0_fake = (TH2D*)filefake->Get("MMnmiss_IMnpip_dE_woK0");
  TH2D *MMnmiss_IMnpim_woK0_fake = (TH2D*)filefake->Get("MMnmiss_IMnpim_dE_woK0");
  TH2D *MMnmiss_IMnpip_woK0_woSm_fake = (TH2D*)filefake->Get("MMnmiss_IMnpip_dE_woK0_woSm");
  TH2D *MMnmiss_IMnpim_woK0_woSp_fake = (TH2D*)filefake->Get("MMnmiss_IMnpim_dE_woK0_woSp");
  TH2D *MMnmiss_IMnpip_woK0_woSm_n_fake = (TH2D*)filefake->Get("MMnmiss_IMnpip_dE_woK0_woSm_n");
  TH2D *MMnmiss_IMnpim_woK0_woSp_n_fake = (TH2D*)filefake->Get("MMnmiss_IMnpim_dE_woK0_woSp_n");
  TH2D *MMnmiss_IMnpip_woK0_woSidn_fake = (TH2D*)filefake->Get("MMnmiss_IMnpip_dE_woK0_woSidn");
  TH2D *MMnmiss_IMnpim_woK0_woSidn_fake = (TH2D*)filefake->Get("MMnmiss_IMnpim_dE_woK0_woSidn");
  TH2D *MMnmiss_IMnpip_woK0_woSid_won_fake = (TH2D*)filefake->Get("MMnmiss_IMnpip_dE_woK0_woSid_won");
  TH2D *MMnmiss_IMnpim_woK0_woSid_won_fake = (TH2D*)filefake->Get("MMnmiss_IMnpim_dE_woK0_woSid_won");
  TH2D *MMnmiss_IMnpip_wK0_woSid_won_fake = (TH2D*)filefake->Get("MMnmiss_IMnpip_dE_wK0_woSid_won");
  TH2D *MMnmiss_IMnpim_wK0_woSid_won_fake = (TH2D*)filefake->Get("MMnmiss_IMnpim_dE_wK0_woSid_won");
  TH2D *MMnmiss_IMpippim_fake = (TH2D*)filefake->Get("MMnmiss_IMpippim_dE");
  TH2D *MMnmiss_IMpippim_wSid_fake = (TH2D*)filefake->Get("MMnmiss_IMpippim_dE_wSid");
  TH2D *MMnmiss_IMpippim_woK0_woSidn_fake = (TH2D*)filefake->Get("MMnmiss_IMpippim_dE_woK0_woSidn");
  TH2D *MMnmiss_IMpippim_woK0_woSid_won_fake = (TH2D*)filefake->Get("MMnmiss_IMpippim_dE_woK0_woSid_won");
  TH2D *MMnmiss_IMpippim_wK0_woSid_won_fake = (TH2D*)filefake->Get("MMnmiss_IMpippim_dE_wK0_woSid_won");
  TH2D *IMnpim_IMnpip_fake = (TH2D*)filefake->Get("IMnpim_IMnpip_dE");
  TH2D *IMnpim_IMnpip_woK0_woSidn_fake = (TH2D*)filefake->Get("IMnpim_IMnpip_dE_woK0_woSidn");
  TH2D *IMnpim_IMnpip_woK0_woSid_won_fake = (TH2D*)filefake->Get("IMnpim_IMnpip_dE_woK0_woSid_won");
  TH2D *q_IMnpipi_woK0_woSidn_fake = (TH2D*)filefake->Get("q_IMnpipi_woK0_woSidn");
  TH2D *q_IMnpipi_woK0_woSid_won_fake = (TH2D*)filefake->Get("q_IMnpipi_woK0_woSid_won");
  TH2D *q_IMnpipi_wK0_woSid_won_fake = (TH2D*)filefake->Get("q_IMnpipi_wK0_woSid_won");
  TH2D *MMnmiss_Mompippim_woK0_woSidn_fake  = (TH2D*)filefake->Get("MMnmiss_Mompippim_dE_woK0_woSidn");
  TH2D *MMnmiss_Mompippim_woK0_woSid_won_fake  = (TH2D*)filefake->Get("MMnmiss_Mompippim_dE_woK0_woSid_won");
  TH2D *MMnmiss_Mompippim_wK0_woSid_won_fake = (TH2D*)filefake->Get("MMnmiss_Mompippim_dE_wK0_woSid_won");
  TH2D *pipmom_MMnmiss_woK0_woSidn_fake = (TH2D*)filefake->Get("pipmom_MMnmiss_dE_woK0_woSidn");
  TH2D *pipmom_MMnmiss_woK0_woSid_won_fake = (TH2D*)filefake->Get("pipmom_MMnmiss_dE_woK0_woSid_won");
  TH2D *pipmom_MMnmiss_wK0_woSid_won_fake = (TH2D*)filefake->Get("pipmom_MMnmiss_dE_wK0_woSid_won");
  TH2D *pipmom_pimmom_woK0_woSidn_fake = (TH2D*)filefake->Get("pipmom_pimmom_dE_woK0_woSidn");
  TH2D *pipmom_pimmom_woK0_woSid_won_fake = (TH2D*)filefake->Get("pipmom_pimmom_dE_woK0_woSid_won");
  TH2D *pipmom_pimmom_wK0_woSid_won_fake = (TH2D*)filefake->Get("pipmom_pimmom_dE_wK0_woSid_won");
  TH2D *pimmom_MMnmiss_woK0_woSidn_fake = (TH2D*)filefake->Get("pimmom_MMnmiss_dE_woK0_woSidn");
  TH2D *pimmom_MMnmiss_woK0_woSid_won_fake = (TH2D*)filefake->Get("pimmom_MMnmiss_dE_woK0_woSid_won");
  TH2D *pimmom_MMnmiss_wK0_woSid_won_fake = (TH2D*)filefake->Get("pimmom_MMnmiss_dE_wK0_woSid_won");
  TH2D *MMom_MMass_woK0_woSidn_fake = (TH2D*)filefake->Get("MMom_MMass_woK0_woSidn");
  TH2D *MMom_MMass_woK0_woSid_won_fake = (TH2D*)filefake->Get("MMom_MMass_woK0_woSid_won");
  TH2D *MMom_MMass_wK0_woSid_won_fake = (TH2D*)filefake->Get("MMom_MMass_wK0_woSid_won");
  TH2D *nmom_MMnmiss_woK0_woSidn_fake = (TH2D*)filefake->Get("nmom_MMnmiss_woK0_woSidn");
  TH2D *nmom_MMnmiss_woK0_woSid_won_fake = (TH2D*)filefake->Get("nmom_MMnmiss_woK0_woSid_won");
  TH2D *nmom_MMnmiss_wK0_woSid_won_fake = (TH2D*)filefake->Get("nmom_MMnmiss_wK0_woSid_won");
  TH2D *q_IMnpipi_wSid_n_Sp_fake = (TH2D*)filefake->Get("q_IMnpipi_wSid_n_Sp");
  TH2D *q_IMnpipi_wSid_n_Sm_fake = (TH2D*)filefake->Get("q_IMnpipi_wSid_n_Sm");
  TH2D *nmom_cosn_woK0_woSidn_fake = (TH2D*)filefake->Get("nmom_cosn_woK0_woSidn");
  TH2D *nmom_cospip_woK0_woSidn_fake = (TH2D*)filefake->Get("nmom_cospip_woK0_woSidn");
  TH2D *nmom_cospim_woK0_woSidn_fake = (TH2D*)filefake->Get("nmom_cospim_woK0_woSidn");
  TH2D *nmom_phinpip_woK0_woSidn_fake = (TH2D*)filefake->Get("nmom_phinpip_woK0_woSidn");
  TH2D *nmom_phinpim_woK0_woSidn_fake = (TH2D*)filefake->Get("nmom_phinpim_woK0_woSidn");
  TH2D *nmom_phipip_woK0_woSidn_fake = (TH2D*)filefake->Get("nmom_phipip_woK0_woSidn");
  TH2D *nmom_phipim_woK0_woSidn_fake = (TH2D*)filefake->Get("nmom_phipim_woK0_woSidn");
  TH2D *nmom_phin_woK0_woSidn_fake = (TH2D*)filefake->Get("nmom_phin_woK0_woSidn");
  TH2D *nmom_pipmom_woK0_woSidn_fake = (TH2D*)filefake->Get("nmom_pipmom_woK0_woSidn");
  TH2D *nmom_pimmom_woK0_woSidn_fake = (TH2D*)filefake->Get("nmom_pimmom_woK0_woSidn");
  TH2D *nmom_cosn_woK0_woSid_won_fake = (TH2D*)filefake->Get("nmom_cosn_woK0_woSid_won");
  TH2D *nmom_cospip_woK0_woSid_won_fake = (TH2D*)filefake->Get("nmom_cospip_woK0_woSid_won");
  TH2D *nmom_cospim_woK0_woSid_won_fake = (TH2D*)filefake->Get("nmom_cospim_woK0_woSid_won");
  TH2D *nmom_phinpip_woK0_woSid_won_fake = (TH2D*)filefake->Get("nmom_phinpip_woK0_woSid_won");
  TH2D *nmom_phinpim_woK0_woSid_won_fake = (TH2D*)filefake->Get("nmom_phinpim_woK0_woSid_won");
  TH2D *nmom_phipip_woK0_woSid_won_fake = (TH2D*)filefake->Get("nmom_phipip_woK0_woSid_won");
  TH2D *nmom_phipim_woK0_woSid_won_fake = (TH2D*)filefake->Get("nmom_phipim_woK0_woSid_won");
  TH2D *nmom_phin_woK0_woSid_won_fake = (TH2D*)filefake->Get("nmom_phin_woK0_woSid_won");
  TH2D *nmom_pipmom_woK0_woSid_won_fake = (TH2D*)filefake->Get("nmom_pipmom_woK0_woSid_won");
  TH2D *nmom_pimmom_woK0_woSid_won_fake = (TH2D*)filefake->Get("nmom_pimmom_woK0_woSid_won");
  TH2D *nmom_cosn_wK0_woSid_won_fake = (TH2D*)filefake->Get("nmom_cosn_wK0_woSid_won");
  TH2D *nmom_cospip_wK0_woSid_won_fake = (TH2D*)filefake->Get("nmom_cospip_wK0_woSid_won");
  TH2D *nmom_cospim_wK0_woSid_won_fake = (TH2D*)filefake->Get("nmom_cospim_wK0_woSid_won");
  TH2D *nmom_phinpip_wK0_woSid_won_fake = (TH2D*)filefake->Get("nmom_phinpip_wK0_woSid_won");
  TH2D *nmom_phinpim_wK0_woSid_won_fake = (TH2D*)filefake->Get("nmom_phinpim_wK0_woSid_won");
  TH2D *nmom_phipip_wK0_woSid_won_fake = (TH2D*)filefake->Get("nmom_phipip_wK0_woSid_won");
  TH2D *nmom_phipim_wK0_woSid_won_fake = (TH2D*)filefake->Get("nmom_phipim_wK0_woSid_won");
  TH2D *nmom_phin_wK0_woSid_won_fake = (TH2D*)filefake->Get("nmom_phin_wK0_woSid_won");
  TH2D *nmom_pipmom_wK0_woSid_won_fake = (TH2D*)filefake->Get("nmom_pipmom_wK0_woSid_won");
  TH2D *nmom_pimmom_wK0_woSid_won_fake = (TH2D*)filefake->Get("nmom_pimmom_wK0_woSid_won");
  double Nnpip_fake = MMnmiss_IMnpip_woK0_fake->GetEntries();
  double Nnpim_fake = MMnmiss_IMnpim_woK0_fake->GetEntries();
  //real data
  TH2D *MMnmiss_IMnpip_woK0_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpip_dE_woK0");
  TH2D *MMnmiss_IMnpim_woK0_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpim_dE_woK0");
  TH2D *MMnmiss_IMnpip_woK0_woSm_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpip_dE_woK0_woSm");
  TH2D *MMnmiss_IMnpim_woK0_woSp_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpim_dE_woK0_woSp");
  TH2D *MMnmiss_IMnpip_woK0_woSm_n_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpip_dE_woK0_woSm_n");
  TH2D *MMnmiss_IMnpim_woK0_woSp_n_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpim_dE_woK0_woSp_n");
  TH2D *MMnmiss_IMnpip_woK0_woSidn_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpip_dE_woK0_woSidn");
  TH2D *MMnmiss_IMnpim_woK0_woSidn_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpim_dE_woK0_woSidn");
  TH2D *MMnmiss_IMnpip_woK0_woSid_won_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpip_dE_woK0_woSid_won");
  TH2D *MMnmiss_IMnpim_woK0_woSid_won_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpim_dE_woK0_woSid_won");
  TH2D *MMnmiss_IMnpip_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpip_dE_wK0_woSid_won");
  TH2D *MMnmiss_IMnpim_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpim_dE_wK0_woSid_won");
  TH2D *MMnmiss_IMpippim_rdata = (TH2D*)filerdata->Get("MMnmiss_IMpippim_dE");
  TH2D *MMnmiss_IMpippim_wSid_rdata = (TH2D*)filerdata->Get("MMnmiss_IMpippim_dE_wSid");
  TH2D *MMnmiss_IMpippim_woK0_woSidn_rdata = (TH2D*)filerdata->Get("MMnmiss_IMpippim_dE_woK0_woSidn");
  TH2D *MMnmiss_IMpippim_woK0_woSid_won_rdata = (TH2D*)filerdata->Get("MMnmiss_IMpippim_dE_woK0_woSid_won");
  TH2D *MMnmiss_IMpippim_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("MMnmiss_IMpippim_dE_wK0_woSid_won");
  TH2D *IMnpim_IMnpip_rdata = (TH2D*)filerdata->Get("IMnpim_IMnpip_dE");
  TH2D *IMnpim_IMnpip_woK0_woSidn_rdata = (TH2D*)filerdata->Get("IMnpim_IMnpip_dE_woK0_woSidn");
  TH2D *IMnpim_IMnpip_woK0_woSid_won_rdata = (TH2D*)filerdata->Get("IMnpim_IMnpip_dE_woK0_woSid_won");
  TH2D *q_IMnpipi_woK0_woSidn_rdata = (TH2D*)filerdata->Get("q_IMnpipi_woK0_woSidn");
  TH2D *q_IMnpipi_woK0_woSid_won_rdata = (TH2D*)filerdata->Get("q_IMnpipi_woK0_woSid_won");
  TH2D *q_IMnpipi_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("q_IMnpipi_wK0_woSid_won");
  TH2D *MMnmiss_Mompippim_woK0_woSidn_rdata  = (TH2D*)filerdata->Get("MMnmiss_Mompippim_dE_woK0_woSidn");
  TH2D *MMnmiss_Mompippim_woK0_woSid_won_rdata  = (TH2D*)filerdata->Get("MMnmiss_Mompippim_dE_woK0_woSid_won");
  TH2D *MMnmiss_Mompippim_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("MMnmiss_Mompippim_dE_wK0_woSid_won");
  TH2D *pipmom_MMnmiss_woK0_woSidn_rdata = (TH2D*)filerdata->Get("pipmom_MMnmiss_dE_woK0_woSidn");
  TH2D *pipmom_MMnmiss_woK0_woSid_won_rdata = (TH2D*)filerdata->Get("pipmom_MMnmiss_dE_woK0_woSid_won");
  TH2D *pipmom_MMnmiss_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("pipmom_MMnmiss_dE_wK0_woSid_won");
  TH2D *pipmom_pimmom_woK0_woSidn_rdata = (TH2D*)filerdata->Get("pipmom_pimmom_dE_woK0_woSidn");
  TH2D *pipmom_pimmom_woK0_woSid_won_rdata = (TH2D*)filerdata->Get("pipmom_pimmom_dE_woK0_woSid_won");
  TH2D *pipmom_pimmom_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("pipmom_pimmom_dE_wK0_woSid_won");
  TH2D *pimmom_MMnmiss_woK0_woSidn_rdata = (TH2D*)filerdata->Get("pimmom_MMnmiss_dE_woK0_woSidn");
  TH2D *pimmom_MMnmiss_woK0_woSid_won_rdata = (TH2D*)filerdata->Get("pimmom_MMnmiss_dE_woK0_woSid_won");
  TH2D *pimmom_MMnmiss_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("pimmom_MMnmiss_dE_wK0_woSid_won");
  TH2D *MMom_MMass_woK0_woSidn_rdata = (TH2D*)filerdata->Get("MMom_MMass_woK0_woSidn");
  TH2D *MMom_MMass_woK0_woSid_won_rdata = (TH2D*)filerdata->Get("MMom_MMass_woK0_woSid_won");
  TH2D *MMom_MMass_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("MMom_MMass_wK0_woSid_won");
  TH2D *nmom_MMnmiss_woK0_woSidn_rdata = (TH2D*)filerdata->Get("nmom_MMnmiss_woK0_woSidn");
  TH2D *nmom_MMnmiss_woK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_MMnmiss_woK0_woSid_won");
  TH2D *nmom_MMnmiss_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_MMnmiss_wK0_woSid_won");
  TH2D *q_IMnpipi_wSid_n_Sp_rdata = (TH2D*)filerdata->Get("q_IMnpipi_wSid_n_Sp");
  TH2D *q_IMnpipi_wSid_n_Sm_rdata = (TH2D*)filerdata->Get("q_IMnpipi_wSid_n_Sm");
  TH2D *nmom_cosn_woK0_woSidn_rdata = (TH2D*)filerdata->Get("nmom_cosn_woK0_woSidn");
  TH2D *nmom_cospip_woK0_woSidn_rdata = (TH2D*)filerdata->Get("nmom_cospip_woK0_woSidn");
  TH2D *nmom_cospim_woK0_woSidn_rdata = (TH2D*)filerdata->Get("nmom_cospim_woK0_woSidn");
  TH2D *nmom_phinpip_woK0_woSidn_rdata = (TH2D*)filerdata->Get("nmom_phinpip_woK0_woSidn");
  TH2D *nmom_phinpim_woK0_woSidn_rdata = (TH2D*)filerdata->Get("nmom_phinpim_woK0_woSidn");
  TH2D *nmom_phipip_woK0_woSidn_rdata = (TH2D*)filerdata->Get("nmom_phipip_woK0_woSidn");
  TH2D *nmom_phipim_woK0_woSidn_rdata = (TH2D*)filerdata->Get("nmom_phipim_woK0_woSidn");
  TH2D *nmom_phin_woK0_woSidn_rdata = (TH2D*)filerdata->Get("nmom_phin_woK0_woSidn");
  TH2D *nmom_pipmom_woK0_woSidn_rdata = (TH2D*)filerdata->Get("nmom_pipmom_woK0_woSidn");
  TH2D *nmom_pimmom_woK0_woSidn_rdata = (TH2D*)filerdata->Get("nmom_pimmom_woK0_woSidn");
  TH2D *nmom_cosn_woK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_cosn_woK0_woSid_won");
  TH2D *nmom_cospip_woK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_cospip_woK0_woSid_won");
  TH2D *nmom_cospim_woK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_cospim_woK0_woSid_won");
  TH2D *nmom_phinpip_woK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_phinpip_woK0_woSid_won");
  TH2D *nmom_phinpim_woK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_phinpim_woK0_woSid_won");
  TH2D *nmom_phipip_woK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_phipip_woK0_woSid_won");
  TH2D *nmom_phipim_woK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_phipim_woK0_woSid_won");
  TH2D *nmom_phin_woK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_phin_woK0_woSid_won");
  TH2D *nmom_pipmom_woK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_pipmom_woK0_woSid_won");
  TH2D *nmom_pimmom_woK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_pimmom_woK0_woSid_won");
  TH2D *nmom_cosn_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_cosn_wK0_woSid_won");
  TH2D *nmom_cospip_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_cospip_wK0_woSid_won");
  TH2D *nmom_cospim_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_cospim_wK0_woSid_won");
  TH2D *nmom_phinpip_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_phinpip_wK0_woSid_won");
  TH2D *nmom_phinpim_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_phinpim_wK0_woSid_won");
  TH2D *nmom_phipip_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_phipip_wK0_woSid_won");
  TH2D *nmom_phipim_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_phipim_wK0_woSid_won");
  TH2D *nmom_phin_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_phin_wK0_woSid_won");
  TH2D *nmom_pipmom_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_pipmom_wK0_woSid_won");
  TH2D *nmom_pimmom_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("nmom_pimmom_wK0_woSid_won");
  double Nnpip_rdata = MMnmiss_IMnpip_woK0_rdata->GetEntries();
  double Nnpim_rdata = MMnmiss_IMnpim_woK0_rdata->GetEntries();
   
  //v2 2020113
  const double scaleFactor[2]= {1.0,
                                nrdata_Sp/fake_Sp
                               };


  //arrays
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

  TH2D* MMnmiss_IMnpip_woK0_woSidn[2]= {
    MMnmiss_IMnpip_woK0_woSidn_rdata,
    MMnmiss_IMnpip_woK0_woSidn_fake
  };
  
  TH2D* MMnmiss_IMnpip_woK0_woSid_won[2]= {
    MMnmiss_IMnpip_woK0_woSid_won_rdata,
    MMnmiss_IMnpip_woK0_woSid_won_fake
  };

  TH2D* MMnmiss_IMnpip_wK0_woSid_won[2]= {
    MMnmiss_IMnpip_wK0_woSid_won_rdata,
    MMnmiss_IMnpip_wK0_woSid_won_fake
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

  TH2D* MMnmiss_IMnpim_woK0_woSidn[2]= {
    MMnmiss_IMnpim_woK0_woSidn_rdata,
    MMnmiss_IMnpim_woK0_woSidn_fake
  };
  
  TH2D* MMnmiss_IMnpim_woK0_woSid_won[2]= {
    MMnmiss_IMnpim_woK0_woSid_won_rdata,
    MMnmiss_IMnpim_woK0_woSid_won_fake
  };

  TH2D* MMnmiss_IMnpim_wK0_woSid_won[2]= {
    MMnmiss_IMnpim_wK0_woSid_won_rdata,
    MMnmiss_IMnpim_wK0_woSid_won_fake
  };

  TH2D* IMnpim_IMnpip_woK0_woSidn[2]= {
    IMnpim_IMnpip_woK0_woSidn_rdata,
    IMnpim_IMnpip_woK0_woSidn_fake
  };
  
  TH2D* IMnpim_IMnpip_woK0_woSid_won[2]= {
    IMnpim_IMnpip_woK0_woSid_won_rdata,
    IMnpim_IMnpip_woK0_woSid_won_fake
  };

  TH2D* MMnmiss_IMpippim_wSid[2]= {
    MMnmiss_IMpippim_wSid_rdata,
    MMnmiss_IMpippim_wSid_fake
  };

  TH2D* MMnmiss_IMpippim_woK0_woSidn[2]= {
    MMnmiss_IMpippim_woK0_woSidn_rdata,
    MMnmiss_IMpippim_woK0_woSidn_fake
  };
  
  TH2D* MMnmiss_IMpippim_woK0_woSid_won[2]= {
    MMnmiss_IMpippim_woK0_woSid_won_rdata,
    MMnmiss_IMpippim_woK0_woSid_won_fake
  };

  TH2D* MMnmiss_IMpippim_wK0_woSid_won[2]= {
    MMnmiss_IMpippim_wK0_woSid_won_rdata,
    MMnmiss_IMpippim_wK0_woSid_won_fake
  };

  TH2D* q_IMnpipi_woK0_woSidn[2]= {
    q_IMnpipi_woK0_woSidn_rdata,
    q_IMnpipi_woK0_woSidn_fake
  };
  
  TH2D* q_IMnpipi_woK0_woSid_won[2]= {
    q_IMnpipi_woK0_woSid_won_rdata,
    q_IMnpipi_woK0_woSid_won_fake
  };

  TH2D* q_IMnpipi_wK0_woSid_won[2]= {
    q_IMnpipi_wK0_woSid_won_rdata,
    q_IMnpipi_wK0_woSid_won_fake
  };

  TH2D* MMnmiss_Mompippim_woK0_woSidn[2]= {
    MMnmiss_Mompippim_woK0_woSidn_rdata,
    MMnmiss_Mompippim_woK0_woSidn_fake
  };
  
  TH2D* MMnmiss_Mompippim_woK0_woSid_won[2]= {
    MMnmiss_Mompippim_woK0_woSid_won_rdata,
    MMnmiss_Mompippim_woK0_woSid_won_fake
  };

  TH2D* MMnmiss_Mompippim_wK0_woSid_won[2]= {
    MMnmiss_Mompippim_wK0_woSid_won_rdata,
    MMnmiss_Mompippim_wK0_woSid_won_fake
  };

  TH2D* pipmom_MMnmiss_woK0_woSidn[2]= {
    pipmom_MMnmiss_woK0_woSidn_rdata,
    pipmom_MMnmiss_woK0_woSidn_fake
  };
  
  TH2D* pipmom_MMnmiss_woK0_woSid_won[2]= {
    pipmom_MMnmiss_woK0_woSid_won_rdata,
    pipmom_MMnmiss_woK0_woSid_won_fake
  };

  TH2D* pipmom_MMnmiss_wK0_woSid_won[2]= {
    pipmom_MMnmiss_wK0_woSid_won_rdata,
    pipmom_MMnmiss_wK0_woSid_won_fake
  };


  TH2D *pipmom_pimmom_woK0_woSidn[2]={
    pipmom_pimmom_woK0_woSidn_rdata,
    pipmom_pimmom_woK0_woSidn_fake
  };
  
  TH2D *pipmom_pimmom_woK0_woSid_won[2]={
    pipmom_pimmom_woK0_woSid_won_rdata,
    pipmom_pimmom_woK0_woSid_won_fake
  };
  
  
  TH2D *pipmom_pimmom_wK0_woSid_won[2]={
    pipmom_pimmom_wK0_woSid_won_rdata,
    pipmom_pimmom_wK0_woSid_won_fake
  };

  TH2D* pimmom_MMnmiss_woK0_woSidn[2]= {
    pimmom_MMnmiss_woK0_woSidn_rdata,
    pimmom_MMnmiss_woK0_woSidn_fake
  };
  
  TH2D* pimmom_MMnmiss_woK0_woSid_won[2]= {
    pimmom_MMnmiss_woK0_woSid_won_rdata,
    pimmom_MMnmiss_woK0_woSid_won_fake
  };

  TH2D* pimmom_MMnmiss_wK0_woSid_won[2]= {
    pimmom_MMnmiss_wK0_woSid_won_rdata,
    pimmom_MMnmiss_wK0_woSid_won_fake
  };

  TH2D* MMom_MMass_woK0_woSidn[2]= {
    MMom_MMass_woK0_woSidn_rdata,
    MMom_MMass_woK0_woSidn_fake
  };
  
  TH2D* MMom_MMass_woK0_woSid_won[2]= {
    MMom_MMass_woK0_woSid_won_rdata,
    MMom_MMass_woK0_woSid_won_fake
  };

  TH2D* MMom_MMass_wK0_woSid_won[2]= {
    MMom_MMass_wK0_woSid_won_rdata,
    MMom_MMass_wK0_woSid_won_fake
  };

  TH2D* nmom_MMnmiss_woK0_woSidn[2]= {
    nmom_MMnmiss_woK0_woSidn_rdata,
    nmom_MMnmiss_woK0_woSidn_fake
  };
  
  TH2D* nmom_MMnmiss_woK0_woSid_won[2]= {
    nmom_MMnmiss_woK0_woSid_won_rdata,
    nmom_MMnmiss_woK0_woSid_won_fake
  };

  TH2D* nmom_MMnmiss_wK0_woSid_won[2]= {
    nmom_MMnmiss_wK0_woSid_won_rdata,
    nmom_MMnmiss_wK0_woSid_won_fake
  };

  TH2D* q_IMnpipi_wSid_n_Sp[2]= {
    q_IMnpipi_wSid_n_Sp_rdata,
    q_IMnpipi_wSid_n_Sp_fake
  };

  TH2D* q_IMnpipi_wSid_n_Sm[2]= {
    q_IMnpipi_wSid_n_Sm_rdata,
    q_IMnpipi_wSid_n_Sm_fake
  };

  TH2D* nmom_cosn_woK0_woSidn[2]= {
    nmom_cosn_woK0_woSidn_rdata,
    nmom_cosn_woK0_woSidn_fake
  };
  
  TH2D* nmom_cosn_woK0_woSid_won[2]= {
    nmom_cosn_woK0_woSid_won_rdata,
    nmom_cosn_woK0_woSid_won_fake
  };

  TH2D* nmom_cospip_woK0_woSidn[2]= {
    nmom_cospip_woK0_woSidn_rdata,
    nmom_cospip_woK0_woSidn_fake
  };
  
  TH2D* nmom_cospip_woK0_woSid_won[2]= {
    nmom_cospip_woK0_woSid_won_rdata,
    nmom_cospip_woK0_woSid_won_fake
  };

  TH2D* nmom_cospim_woK0_woSidn[2]= {
    nmom_cospim_woK0_woSidn_rdata,
    nmom_cospim_woK0_woSidn_fake
  };
  
  TH2D* nmom_cospim_woK0_woSid_won[2]= {
    nmom_cospim_woK0_woSid_won_rdata,
    nmom_cospim_woK0_woSid_won_fake
  };

  TH2D* nmom_phinpip_woK0_woSidn[2]= {
    nmom_phinpip_woK0_woSidn_rdata,
    nmom_phinpip_woK0_woSidn_fake
  };
  
  TH2D* nmom_phinpip_woK0_woSid_won[2]= {
    nmom_phinpip_woK0_woSid_won_rdata,
    nmom_phinpip_woK0_woSid_won_fake
  };
  
  TH2D* nmom_phipip_woK0_woSidn[2]= {
    nmom_phipip_woK0_woSidn_rdata,
    nmom_phipip_woK0_woSidn_fake
  };
  
  TH2D* nmom_phipip_woK0_woSid_won[2]= {
    nmom_phipip_woK0_woSid_won_rdata,
    nmom_phipip_woK0_woSid_won_fake
  };

  TH2D* nmom_phinpim_woK0_woSidn[2]= {
    nmom_phinpim_woK0_woSidn_rdata,
    nmom_phinpim_woK0_woSidn_fake
  };
  
  TH2D* nmom_phinpim_woK0_woSid_won[2]= {
    nmom_phinpim_woK0_woSid_won_rdata,
    nmom_phinpim_woK0_woSid_won_fake
  };
  
  TH2D* nmom_phipim_woK0_woSidn[2]= {
    nmom_phipim_woK0_woSidn_rdata,
    nmom_phipim_woK0_woSidn_fake
  };
  
  TH2D* nmom_phipim_woK0_woSid_won[2]= {
    nmom_phipim_woK0_woSid_won_rdata,
    nmom_phipim_woK0_woSid_won_fake
  };
  
  TH2D* nmom_phin_woK0_woSidn[2]= {
    nmom_phin_woK0_woSidn_rdata,
    nmom_phin_woK0_woSidn_fake
  };
  
  TH2D* nmom_phin_woK0_woSid_won[2]= {
    nmom_phin_woK0_woSid_won_rdata,
    nmom_phin_woK0_woSid_won_fake
  };

  TH2D* nmom_pipmom_woK0_woSidn[2]= {
    nmom_pipmom_woK0_woSidn_rdata,
    nmom_pipmom_woK0_woSidn_fake
  };
  
  TH2D* nmom_pipmom_woK0_woSid_won[2]= {
    nmom_pipmom_woK0_woSid_won_rdata,
    nmom_pipmom_woK0_woSid_won_fake
  };

  TH2D* nmom_pimmom_woK0_woSidn[2]= {
    nmom_pimmom_woK0_woSidn_rdata,
    nmom_pimmom_woK0_woSidn_fake
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

  TH2D* nmom_phinpip_wK0_woSid_won[2]= {
    nmom_phinpip_wK0_woSid_won_rdata,
    nmom_phinpip_wK0_woSid_won_fake
  };

  TH2D* nmom_phinpim_wK0_woSid_won[2]= {
    nmom_phinpim_wK0_woSid_won_rdata,
    nmom_phinpim_wK0_woSid_won_fake
  };
  
  TH2D* nmom_phipip_wK0_woSid_won[2]= {
    nmom_phipip_wK0_woSid_won_rdata,
    nmom_phipip_wK0_woSid_won_fake
  };

  TH2D* nmom_phipim_wK0_woSid_won[2]= {
    nmom_phipim_wK0_woSid_won_rdata,
    nmom_phipim_wK0_woSid_won_fake
  };
  
  TH2D* nmom_phin_wK0_woSid_won[2]= {
    nmom_phin_wK0_woSid_won_rdata,
    nmom_phin_wK0_woSid_won_fake
  };


  TH2D* nmom_pipmom_wK0_woSid_won[2]= {
    nmom_pipmom_wK0_woSid_won_rdata,
    nmom_pipmom_wK0_woSid_won_fake
  };

  TH2D* nmom_pimmom_wK0_woSid_won[2]= {
    nmom_pimmom_wK0_woSid_won_rdata,
    nmom_pimmom_wK0_woSid_won_fake
  };


  for(int i=0; i<2; i++) {
    MMnmiss_IMnpipi_woK0_wSid_Sp[i]->Scale(scaleFactor[i]);
    MMnmiss_IMnpip_woK0[i]->Scale(scaleFactor[i]);
    MMnmiss_IMnpip_woK0_woSm[i]->Scale(scaleFactor[i]);
    MMnmiss_IMnpip_woK0_woSm_n[i]->Scale(scaleFactor[i]);
    MMnmiss_IMnpip_woK0_woSidn[i]->Scale(scaleFactor[i]);
    MMnmiss_IMnpip_woK0_woSid_won[i]->Scale(scaleFactor[i]);
    MMnmiss_IMnpip_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    MMnmiss_IMnpipi_woK0_wSid_Sm[i]->Scale(scaleFactor[i]);
    MMnmiss_IMnpim_woK0[i]->Scale(scaleFactor[i]);
    MMnmiss_IMnpim_woK0_woSp[i]->Scale(scaleFactor[i]);
    MMnmiss_IMnpim_woK0_woSp_n[i]->Scale(scaleFactor[i]);
    MMnmiss_IMnpim_woK0_woSidn[i]->Scale(scaleFactor[i]);
    MMnmiss_IMnpim_woK0_woSid_won[i]->Scale(scaleFactor[i]);
    MMnmiss_IMnpim_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    IMnpim_IMnpip_woK0_woSidn[i]->Scale(scaleFactor[i]);
    IMnpim_IMnpip_woK0_woSid_won[i]->Scale(scaleFactor[i]);
    MMnmiss_IMpippim_wSid[i]->Scale(scaleFactor[i]);
    MMnmiss_IMpippim_woK0_woSidn[i]->Scale(scaleFactor[i]);
    MMnmiss_IMpippim_woK0_woSid_won[i]->Scale(scaleFactor[i]);
    MMnmiss_IMpippim_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    q_IMnpipi_woK0_woSidn[i]->Scale(scaleFactor[i]);
    q_IMnpipi_woK0_woSid_won[i]->Scale(scaleFactor[i]);
    q_IMnpipi_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    MMnmiss_Mompippim_woK0_woSidn[i]->Scale(scaleFactor[i]);
    MMnmiss_Mompippim_woK0_woSid_won[i]->Scale(scaleFactor[i]);
    MMnmiss_Mompippim_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    pipmom_MMnmiss_woK0_woSidn[i]->Scale(scaleFactor[i]);
    pipmom_MMnmiss_woK0_woSid_won[i]->Scale(scaleFactor[i]);
    pipmom_MMnmiss_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    pipmom_pimmom_woK0_woSidn[i]->Scale(scaleFactor[i]);
    pipmom_pimmom_woK0_woSid_won[i]->Scale(scaleFactor[i]);
    pipmom_pimmom_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    pimmom_MMnmiss_woK0_woSidn[i]->Scale(scaleFactor[i]);
    pimmom_MMnmiss_woK0_woSid_won[i]->Scale(scaleFactor[i]);
    pimmom_MMnmiss_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    MMom_MMass_woK0_woSidn[i]->Scale(scaleFactor[i]);
    MMom_MMass_woK0_woSid_won[i]->Scale(scaleFactor[i]);
    MMom_MMass_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    nmom_MMnmiss_woK0_woSidn[i]->Scale(scaleFactor[i]);
    nmom_MMnmiss_woK0_woSid_won[i]->Scale(scaleFactor[i]);
    nmom_MMnmiss_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    q_IMnpipi_wSid_n_Sp[i]->Scale(scaleFactor[i]);
    q_IMnpipi_wSid_n_Sm[i]->Scale(scaleFactor[i]);
    nmom_cosn_woK0_woSidn[i]->Scale(scaleFactor[i]);
    nmom_cosn_woK0_woSid_won[i]->Scale(scaleFactor[i]);
    nmom_cospip_woK0_woSidn[i]->Scale(scaleFactor[i]);
    nmom_cospip_woK0_woSid_won[i]->Scale(scaleFactor[i]);
    nmom_cospim_woK0_woSidn[i]->Scale(scaleFactor[i]);
    nmom_cospim_woK0_woSid_won[i]->Scale(scaleFactor[i]);
    nmom_phinpip_woK0_woSidn[i]->Scale(scaleFactor[i]);
    nmom_phinpip_woK0_woSid_won[i]->Scale(scaleFactor[i]);
    nmom_phipip_woK0_woSidn[i]->Scale(scaleFactor[i]);
    nmom_phipip_woK0_woSid_won[i]->Scale(scaleFactor[i]);
    nmom_phinpim_woK0_woSidn[i]->Scale(scaleFactor[i]);
    nmom_phinpim_woK0_woSid_won[i]->Scale(scaleFactor[i]);
    nmom_phipim_woK0_woSidn[i]->Scale(scaleFactor[i]);
    nmom_phipim_woK0_woSid_won[i]->Scale(scaleFactor[i]);
    nmom_phin_woK0_woSidn[i]->Scale(scaleFactor[i]);
    nmom_phin_woK0_woSid_won[i]->Scale(scaleFactor[i]);
    nmom_pipmom_woK0_woSidn[i]->Scale(scaleFactor[i]);
    nmom_pipmom_woK0_woSid_won[i]->Scale(scaleFactor[i]);
    nmom_pimmom_woK0_woSidn[i]->Scale(scaleFactor[i]);
    nmom_pimmom_woK0_woSid_won[i]->Scale(scaleFactor[i]);
    nmom_cosn_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    nmom_cospip_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    nmom_cospim_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    nmom_phinpip_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    nmom_phinpim_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    nmom_phipip_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    nmom_phipim_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    nmom_phin_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    nmom_pipmom_wK0_woSid_won[i]->Scale(scaleFactor[i]);
    nmom_pimmom_wK0_woSid_won[i]->Scale(scaleFactor[i]);
  }


  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sp_mc = (TH2D*)MMnmiss_IMnpipi_woK0_wSid_Sp[1]->Clone("MMnmiss_IMnpipi_woK0_wSid_Sp_mc");

  TCanvas *cMMnmiss_IMnpipi_woK0_wSid_Sp = new TCanvas("cMMnmiss_IMnpipi_woK0_wSid_Sp","cMMnmiss_IMnpipi_woK0_wSid_Sp",1200,800);
  cMMnmiss_IMnpipi_woK0_wSid_Sp->Divide(2,1);
  cMMnmiss_IMnpipi_woK0_wSid_Sp->cd(1);
  MMnmiss_IMnpipi_woK0_wSid_Sp[0]->SetTitle("#splitline{MMnmiss_IMnpipi_woK0_wSid_Sp}{  real Data}");
  MMnmiss_IMnpipi_woK0_wSid_Sp[0]->SetMinimum(1);
  MMnmiss_IMnpipi_woK0_wSid_Sp[0]->Draw("colz");
  cMMnmiss_IMnpipi_woK0_wSid_Sp->cd(2);
  MMnmiss_IMnpipi_woK0_wSid_Sp_mc->SetTitle("#splitline{MMnmiss_IMnpipi_woK0_wSid_Sp}{  MC sum}");
  MMnmiss_IMnpipi_woK0_wSid_Sp_mc->SetMaximum(MMnmiss_IMnpipi_woK0_wSid_Sp[0]->GetMaximum());
  MMnmiss_IMnpipi_woK0_wSid_Sp_mc->SetMinimum(1);
  MMnmiss_IMnpipi_woK0_wSid_Sp_mc->Draw("colz");

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
  MMnmiss_IMnpipi_woK0_wSid_Sm[0]->Draw("colz");
  cMMnmiss_IMnpipi_woK0_wSid_Sm->cd(2);
  MMnmiss_IMnpipi_woK0_wSid_Sm_mc->SetTitle("#splitline{MMnmiss_IMnpipi_woK0_wSid_Sm}{  MC sum}");
  MMnmiss_IMnpipi_woK0_wSid_Sm_mc->SetMaximum(MMnmiss_IMnpipi_woK0_wSid_Sm[0]->GetMaximum());
  MMnmiss_IMnpipi_woK0_wSid_Sm_mc->SetMinimum(1);
  MMnmiss_IMnpipi_woK0_wSid_Sm_mc->Draw("colz");

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
  MMnmiss_IMnpip_woK0_rdata->Draw("colz");
  cMMnmiss_IMnpip_woK0->cd(2);
  MMnmiss_IMnpip_woK0_mc->SetTitle("#splitline{MMnmiss_IMnpip_woK0}{  MC sum}");
  MMnmiss_IMnpip_woK0_mc->SetMaximum(MMnmiss_IMnpip_woK0_rdata->GetMaximum());
  MMnmiss_IMnpip_woK0_mc->SetMinimum(1);
  MMnmiss_IMnpip_woK0_mc->Draw("colz");

  //w/o K0, w/o Sm mode
  TH2D *MMnmiss_IMnpip_woK0_woSm_mc = (TH2D*)MMnmiss_IMnpip_woK0_woSm[1]->Clone("MMnmiss_IMnpip_woK0_woSm_mc");

  //w/o K0, w/o Sm mode, selecting missing neutron
  TH2D *MMnmiss_IMnpip_woK0_woSm_n_mc = (TH2D*)MMnmiss_IMnpip_woK0_woSm_n[1]->Clone("MMnmiss_IMnpip_woK0_woSm_n_mc");

  TCanvas *cMMnmiss_IMnpip_woK0_woSm = new TCanvas("cMMnmiss_IMnpip_woK0_woSm","cMMnmiss_IMnpip_woK0_woSm",1200,800);
  cMMnmiss_IMnpip_woK0_woSm->Divide(2,1);
  cMMnmiss_IMnpip_woK0_woSm->cd(1);
  MMnmiss_IMnpip_woK0_woSm_rdata->SetTitle("#splitline{MMnmiss_IMnpip_woK0_woSm}{  Real data}");
  MMnmiss_IMnpip_woK0_woSm_rdata->SetMinimum(1);
  MMnmiss_IMnpip_woK0_woSm_rdata->Draw("colz");
  cMMnmiss_IMnpip_woK0_woSm->cd(2);
  MMnmiss_IMnpip_woK0_woSm_mc->SetTitle("#splitline{MMnmiss_IMnpip_woK0_woSm}{  MC sum}");
  MMnmiss_IMnpip_woK0_woSm_mc->SetMinimum(1);
  MMnmiss_IMnpip_woK0_woSm_mc->SetMaximum(MMnmiss_IMnpip_woK0_woSm_rdata->GetMaximum());
  MMnmiss_IMnpip_woK0_woSm_mc->Draw("colz");

  TCanvas *cIMnpip_woK0_woSm_n = new TCanvas("cIMnpip_woK0_woSm_n","cIMnpip_woK0_woSm_n");
  cIMnpip_woK0_woSm_n->cd();
  TH1D *IMnpip_woK0_woSm_n[2];
  for(int i=0; i<2; i++) {
    IMnpip_woK0_woSm_n[i] = (TH1D*)MMnmiss_IMnpip_woK0_woSm_n[i]->ProjectionX(Form("IMnpip_woK0_woSm_n_%s",name[i]));
    IMnpip_woK0_woSm_n[i]->SetLineColor(colordef[i]);
  }
  TH1D* IMnpip_woK0_woSm_n_mc = (TH1D*)MMnmiss_IMnpip_woK0_woSm_n_mc->ProjectionX("IMnpip_woK0_woSm_n_mc");
  IMnpip_woK0_woSm_n[0]->Draw("HE");
  IMnpip_woK0_woSm_n_mc->SetLineColor(6);
  IMnpip_woK0_woSm_n_mc->Draw("same");

  //w/o K0, w/o ((Sp or Sm) & missing neutron)
  TH2D *MMnmiss_IMnpip_woK0_woSid_won_mc = (TH2D*)MMnmiss_IMnpip_woK0_woSid_won[1]->Clone("MMnmiss_IMnpip_woK0_woSid_won_mc");

  TCanvas *cMMnmiss_IMnpip_woSidn = new TCanvas("cMMnmiss_IMnpip_woSidn","cMMnmiss_IMnpip_woSidn",1200,800);
  cMMnmiss_IMnpip_woSidn->Divide(2,1);
  cMMnmiss_IMnpip_woSidn->cd(1);
  MMnmiss_IMnpip_woK0_woSid_won_rdata->SetTitle("#splitline{MMnmiss_IMnpip_woK0_woSid_won}{  Real data}");
  MMnmiss_IMnpip_woK0_woSid_won_rdata->SetMinimum(1);
  MMnmiss_IMnpip_woK0_woSid_won_rdata->Draw("colz");
  cMMnmiss_IMnpip_woSidn->cd(2);
  MMnmiss_IMnpip_woK0_woSid_won_mc->SetTitle("#splitline{MMnmiss_IMnpip_woK0_woSid_won}{  MC sum}");
  MMnmiss_IMnpip_woK0_woSid_won_mc->SetMinimum(1);
  MMnmiss_IMnpip_woK0_woSid_won_mc->SetMaximum(MMnmiss_IMnpip_woK0_woSid_won_rdata->GetMaximum());
  MMnmiss_IMnpip_woK0_woSid_won_mc->Draw("colz");

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
  TCanvas *cMMnmiss_woK0_woSid_won_ratio = new TCanvas("cMMnmiss_woK0_woSid_won_ratio","cMMnmiss_woK0_woSid_ratio");
  TH1D* MMnmiss_woK0_woSid_won_ratio = (TH1D*)MMnmiss_woK0_woSid_won[0]->Clone("MMnmiss_woK0_woSid_won_ratio");
  MMnmiss_woK0_woSid_won_ratio->Divide(MMnmiss_woK0_woSid_won_mc);
  MMnmiss_woK0_woSid_won_ratio->SetTitle("Data/MC");
  MMnmiss_woK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  MMnmiss_woK0_woSid_won_ratio->Draw("HE");


  TF1 *sgf_MMnmiss = new TF1("sgf_MMnmiss","[0]*exp(-0.5*pow((x-[1])/([2]+(x<[1])*[3]*(x-[1])),2))");
  //TF1 *fgaus_MMnmiss_high = new TF1("fgaus_MMnmiss_high","gaus");
  //sgf_MMnmiss->SetParameter(0,1.82171e+00);
  //sgf_MMnmiss->SetParameter(1,8.56016e-01);
  //sgf_MMnmiss->SetParameter(2,6.81677e-01);
  //sgf_MMnmiss->SetLineColor(2);

  MMnmiss_woK0_woSid_won_ratio->Fit("sgf_MMnmiss","R","",1.06,1.5);
  //MMnmiss_woK0_woSid_won_ratio->Fit("fgaus_MMnmiss_high","R+","",1.116,1.5);
  Double_t param_MMnmiss[7];
  sgf_MMnmiss->GetParameters(&param_MMnmiss[0]);
  //fgaus_MMnmiss_high->GetParameters(&param_MMnmiss[4]);
  


  //first fit
  TF1 *evalf_MMnmiss = new TF1("evalf_MMnmiss",func_MMnmiss,0,1.5,18);
  /*
  evalf_MMnmiss->SetParameter(0,0.002689);
  evalf_MMnmiss->SetParameter(1,1.115);
  evalf_MMnmiss->SetParameter(2,2.171);
  evalf_MMnmiss->SetParameter(3,-5.516);
  evalf_MMnmiss->SetParameter(4,0.002196);
  evalf_MMnmiss->SetParameter(5,5.456);
  evalf_MMnmiss->SetParameter(6,3.889);
  evalf_MMnmiss->SetParameter(7,-2.033);
  evalf_MMnmiss->SetParameter(8,-5.099);
  evalf_MMnmiss->SetParameter(9,1.647);
  evalf_MMnmiss->SetParameter(10,10);
  evalf_MMnmiss->SetParameter(11,1.115);
  */

  MMnmiss_woK0_woSid_won_ratio->Fit("evalf_MMnmiss");
  evalf_MMnmiss->SetLineColor(4);
  evalf_MMnmiss->Draw("same");
  //evalf_MMnmiss->Print();
  
  //second fit 
  TF1 *evalf_MMnmiss_corr = new TF1("evalf_MMnmiss_corr",func_MMnmiss_corr,0.0,1.50,7);
  evalf_MMnmiss_corr->SetParameter(0,-5260.37);
  evalf_MMnmiss_corr->SetParameter(1,28408.5);
  evalf_MMnmiss_corr->SetParameter(2,-63542.6);
  evalf_MMnmiss_corr->SetParameter(3,75328.0);
  evalf_MMnmiss_corr->SetParameter(4,-49941.3);
  evalf_MMnmiss_corr->SetParameter(5,17554.7);
  evalf_MMnmiss_corr->SetParameter(6,-2556.0);
  //MMnmiss_woK0_woSid_won_ratio->Fit("evalf_MMnmiss_corr");
  
  TF1 *evalf_MMnmiss_corr2 = new TF1("evalf_MMnmiss_corr2",func_MMnmiss_corr2,0.0,1.50,7);
  evalf_MMnmiss_corr2->SetParameter(0,-856.526);
  evalf_MMnmiss_corr2->SetParameter(1,4911.96);
  evalf_MMnmiss_corr2->SetParameter(2,-11598.5);
  evalf_MMnmiss_corr2->SetParameter(3,14448.0);
  evalf_MMnmiss_corr2->SetParameter(4,-10014.4);
  evalf_MMnmiss_corr2->SetParameter(5,3663.12);
  evalf_MMnmiss_corr2->SetParameter(6,-552.666);
  //MMnmiss_woK0_woSid_won_ratio->Fit("evalf_MMnmiss_corr2");


  TH2D *MMom_MMass_woK0_woSid_won_mc = (TH2D*)MMom_MMass_woK0_woSid_won[1]->Clone("MMom_MMass_woK0_woSid_won_mc");
  TCanvas *cMMom_MMass_woK0_woSid_won = new TCanvas("cMMom_MMass_woK0_woSid_won","cMMom_MMass_woK0_woSid_won");
  cMMom_MMass_woK0_woSid_won->Divide(2,1);
  cMMom_MMass_woK0_woSid_won->cd(1);
  MMom_MMass_woK0_woSid_won[0]->SetMinimum(1);
  MMom_MMass_woK0_woSid_won[0]->Draw("colz");
  cMMom_MMass_woK0_woSid_won->cd(2);
  MMom_MMass_woK0_woSid_won_mc->SetMinimum(1);
  MMom_MMass_woK0_woSid_won_mc->SetMaximum(MMom_MMass_woK0_woSid_won[0]->GetMaximum());
  MMom_MMass_woK0_woSid_won_mc->Draw("colz");

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


  TCanvas *cMMom_woK0_woSid_won_ratio = new TCanvas("cMMom_woK0_woSid_won_ratio","cMMom_woK0_woSid_won_ratio");
  cMMom_woK0_woSid_won_ratio->cd();
  TH1D* MMom_woK0_woSid_won_ratio = (TH1D*)MMom_woK0_woSid_won[0]->Clone("MMom_woK0_woSid_won_ratio");
  MMom_woK0_woSid_won_ratio->Divide(MMom_woK0_woSid_won_mc);
  MMom_woK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  MMom_woK0_woSid_won_ratio->SetTitle("Data/MC");
  MMom_woK0_woSid_won_ratio->Draw("HE");

  Double_t param_MMom[5];
  TF1 *evalf_MMom = new TF1("evalf_MMom","pol5",0.4,1.5);
  evalf_MMom->SetLineColor(4);
  evalf_MMom->SetTitle("MMom");
  MMom_woK0_woSid_won_ratio->Fit("evalf_MMom","","",0.4,1.5);


  //projection to IMnpip (miss n & Sigma+/-)
  TCanvas *cIMnpip_woK0_woSid_won = new TCanvas("cIMnpip_woK0_woSid_won","cIMnpip_woK0_woSid_won");
  cIMnpip_woK0_woSid_won->cd();
  TH1D* IMnpip_woK0_woSid_won[6];
  for(int i=0; i<2; i++)IMnpip_woK0_woSid_won[i] = (TH1D*)MMnmiss_IMnpip_woK0_woSid_won[i]->ProjectionX(Form("IMnpip_woK0_woSid_won_%s",name[i]));
  IMnpip_woK0_woSid_won[0]->Draw("HE");//rdata
  TH1D* IMnpip_woK0_woSid_won_mc = (TH1D*)MMnmiss_IMnpip_woK0_woSid_won_mc->ProjectionX("IMnpip_woK0_woSid_won_mc");
  IMnpip_woK0_woSid_won_mc->SetLineColor(6);
  IMnpip_woK0_woSid_won_mc->Draw("HEsame");

  TCanvas *cIMnpip_woK0_woSid_won_ratio = new TCanvas("cIMnpip_woK0_woSid_won_ratio","cIMnpip_woK0_woSid_won_ratio");
  cIMnpip_woK0_woSid_won_ratio->cd();
  TH1D* IMnpip_woK0_woSid_won_ratio = (TH1D*)IMnpip_woK0_woSid_won[0]->Clone("IMnpip_woK0_woSid_won_ratio");
  IMnpip_woK0_woSid_won_ratio->Divide(IMnpip_woK0_woSid_won_mc);
  IMnpip_woK0_woSid_won_ratio->SetTitle("IMnpip_woK0_woSid_won Data/MC");
  IMnpip_woK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  IMnpip_woK0_woSid_won_ratio->Draw("HE");

  TF1* fgaus_IMnpip_1 = new TF1("fgaus_IMnpip_1","gaus",1.06,1.10);
  fgaus_IMnpip_1->SetParameters(0,1.636);
  fgaus_IMnpip_1->SetParameters(1,1.102);
  fgaus_IMnpip_1->SetParameters(2,0.02845);
  IMnpip_woK0_woSid_won_ratio->Fit("fgaus_IMnpip_1","R","",1.08,1.10);

  TF1* fexpo_IMnpip_2 = new TF1("fexpo_IMnpip_2","expo",1.10,1.25);
  fexpo_IMnpip_2->SetParameters(0,1.667);
  fexpo_IMnpip_2->SetParameters(1,-1.117);
  IMnpip_woK0_woSid_won_ratio->Fit("fexpo_IMnpip_2","R+","",1.10,1.25);

  TF1* fgaus_IMnpip_3 = new TF1("fgaus_IMnpip_3","gaus",1.25,2.0);
  fgaus_IMnpip_3->SetParameters(0,10.83);
  fgaus_IMnpip_3->SetParameters(1,5.094);
  fgaus_IMnpip_3->SetParameters(2,1.87);
  IMnpip_woK0_woSid_won_ratio->Fit("fgaus_IMnpip_3","R+","",1.25,2.0);

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

  TCanvas *cMMnmiss_IMnpim_woK0 = new TCanvas("cMMnmiss_IMnpim_woK0","cMMnmiss_IMnpim_woK0",1200,800);
  cMMnmiss_IMnpim_woK0->Divide(2,1);
  cMMnmiss_IMnpim_woK0->cd(1);
  MMnmiss_IMnpim_woK0_rdata->SetTitle("#splitline{MMnmiss_IMnpim_woK0}{  Real data}");
  MMnmiss_IMnpim_woK0_rdata->SetMinimum(1);
  MMnmiss_IMnpim_woK0_rdata->Draw("colz");
  cMMnmiss_IMnpim_woK0->cd(2);
  MMnmiss_IMnpim_woK0_mc->SetTitle("#splitline{MMnmiss_IMnpim_woK0}{ MC sum}");
  MMnmiss_IMnpim_woK0_mc->SetMinimum(1);
  MMnmiss_IMnpim_woK0_mc->SetMaximum(MMnmiss_IMnpim_woK0_rdata->GetMaximum());
  MMnmiss_IMnpim_woK0_mc->Draw("colz");

  //w/o K0, w/o Sp mode
  TH2D *MMnmiss_IMnpim_woK0_woSp_mc = (TH2D*)MMnmiss_IMnpim_woK0_woSp[1]->Clone("MMnmiss_IMnpim_woK0_woSp_mc");
  //w/o K0, w/o Sp mode, selecting missing neutron
  TH2D *MMnmiss_IMnpim_woK0_woSp_n_mc = (TH2D*)MMnmiss_IMnpim_woK0_woSp_n[1]->Clone("MMnmiss_IMnpim_woK0_woSp_n_mc");

  TCanvas *cMMnmiss_IMnpim_woK0_woSp = new TCanvas("cMMnmiss_IMnpim_woK0_woSp","cMMnmiss_IMnpim_woK0_woSp",1200,800);
  cMMnmiss_IMnpim_woK0_woSp->Divide(2,1);
  cMMnmiss_IMnpim_woK0_woSp->cd(1);
  MMnmiss_IMnpim_woK0_woSp_rdata->SetTitle("#splitline{MMnmiss_IMnpim_woK0_woSp}{  Real data}");
  MMnmiss_IMnpim_woK0_woSp_rdata->SetMinimum(1);
  MMnmiss_IMnpim_woK0_woSp_rdata->Draw("colz");
  cMMnmiss_IMnpim_woK0_woSp->cd(2);
  MMnmiss_IMnpim_woK0_woSp_mc->SetTitle("#splitline{MMnmiss_IMnpim_woK0_woSp}{  MC sum}");
  MMnmiss_IMnpim_woK0_woSp_mc->SetMinimum(1);
  MMnmiss_IMnpim_woK0_woSp_mc->SetMaximum(MMnmiss_IMnpim_woK0_woSp_rdata->GetMaximum());
  MMnmiss_IMnpim_woK0_woSp_mc->Draw("colz");

  TCanvas *cIMnpim_woK0_woSp_n = new TCanvas("cIMnpim_woK0_woSp_n","cIMnpim_woK0_woSp_n");
  cIMnpim_woK0_woSp_n->cd();
  TH1D *IMnpim_woK0_woSp_n[2];
  for(int i=0; i<2; i++) {
    IMnpim_woK0_woSp_n[i] = (TH1D*)MMnmiss_IMnpim_woK0_woSp_n[i]->ProjectionX(Form("IMnpim_woK0_woSp_n_%s",name[i]));
    IMnpim_woK0_woSp_n[i]->SetLineColor(colordef[i]);
  }
  TH1D* IMnpim_woK0_woSp_n_mc = (TH1D*) MMnmiss_IMnpim_woK0_woSp_n_mc->ProjectionX("IMnpim_woK0_woSp_n_mc");
  IMnpim_woK0_woSp_n[0]->Draw("HE");
  //for(int i=1;i<6;i++)IMnpim_woK0_woSp_n[i]->Draw("HEsame");
  IMnpim_woK0_woSp_n_mc->SetLineColor(6);
  IMnpim_woK0_woSp_n_mc->Draw("same");


  //w/o K0, w/o ((Sp or Sm) & missing neutron)
  TH2D *MMnmiss_IMnpim_woK0_woSid_won_mc = (TH2D*)MMnmiss_IMnpim_woK0_woSid_won[1]->Clone("MMnmiss_IMnpim_woK0_woSid_won_mc");

  TCanvas *cMMnmiss_IMnpim_woK0_woSid_won = new TCanvas("cMMnmiss_IMnpim_woK0_woSid_won","cMMnmiss_IMnpim_woK0_woSid_won",1200,800);
  cMMnmiss_IMnpim_woK0_woSid_won->Divide(2,1);
  cMMnmiss_IMnpim_woK0_woSid_won->cd(1);
  MMnmiss_IMnpim_woK0_woSid_won_rdata->SetTitle("#splitline{MMnmiss_IMnpim_woK0_woSid_won}{  Real data}");
  MMnmiss_IMnpim_woK0_woSid_won_rdata->SetMinimum(1);
  MMnmiss_IMnpim_woK0_woSid_won_rdata->Draw("colz");
  cMMnmiss_IMnpim_woK0_woSid_won->cd(2);
  MMnmiss_IMnpim_woK0_woSid_won_mc->SetTitle("#splitline{MMnmiss_IMnpim_woK0_woSid_won}{  MC sum}");
  MMnmiss_IMnpim_woK0_woSid_won_mc->SetMinimum(1);
  MMnmiss_IMnpim_woK0_woSid_won_mc->SetMaximum(MMnmiss_IMnpim_woK0_woSid_won_rdata->GetMaximum());
  MMnmiss_IMnpim_woK0_woSid_won_mc->Draw("colz");

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

  TCanvas *cIMnpim_woK0_woSid_won_ratio = new TCanvas("cIMnpim_woK0_woSid_won_ratio","cIMnpim_woK0_woSid_won_ratio");
  cIMnpim_woK0_woSid_won_ratio->cd();
  TH1D* IMnpim_woK0_woSid_won_ratio = (TH1D*)IMnpim_woK0_woSid_won[0]->Clone("IMnpim_woK0_woSid_won_ratio");
  IMnpim_woK0_woSid_won_ratio->Divide(IMnpim_woK0_woSid_won_mc);
  IMnpim_woK0_woSid_won_ratio->SetTitle("IMnpim_woK0_woSid_won Data/MC");
  IMnpim_woK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  IMnpim_woK0_woSid_won_ratio->Draw("HE");

  TF1* pol3_IMnpim_1 = new TF1("pol3_IMnpim_1","pol3",1.00,1.10);
  pol3_IMnpim_1->SetParameter(0,-32074.9);
  pol3_IMnpim_1->SetParameter(1,85205.3);
  pol3_IMnpim_1->SetParameter(2,-75374.);
  pol3_IMnpim_1->SetParameter(3,22204.);
  IMnpim_woK0_woSid_won_ratio->Fit("pol3_IMnpim_1","R","",1.07,1.10);

  TF1* pol3_IMnpim_2 = new TF1("pol3_IMnpim_2","pol3",1.10,2.00);
  pol3_IMnpim_2->SetParameter(0,38.13);
  pol3_IMnpim_2->SetParameter(1,-62.3139);
  pol3_IMnpim_2->SetParameter(2,32.172);
  pol3_IMnpim_2->SetParameter(3,-4.81);
  IMnpim_woK0_woSid_won_ratio->Fit("pol3_IMnpim_2","R+","",1.10,2.00);

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
  for(int i=1; i<2; i++) MMnmiss_IMpippim_wSid_mc->Add(MMnmiss_IMpippim_wSid[i]);

  TCanvas *cMMnmiss_IMpippim_wSid = new TCanvas("cMMnmiss_IMpippim_wSid","cMMnmiss_IMpippim_wSid",1200,800);
  cMMnmiss_IMpippim_wSid->Divide(2,1);
  cMMnmiss_IMpippim_wSid->cd(1);
  MMnmiss_IMpippim_wSid[0]->SetTitle("#splitline{MMnmiss_IMpippim_wSid}{  Real data}");
  MMnmiss_IMpippim_wSid[0]->SetMinimum(1);
  MMnmiss_IMpippim_wSid[0]->Draw("colz");
  cMMnmiss_IMpippim_wSid->cd(2);
  MMnmiss_IMpippim_wSid_mc->SetTitle("#splitline{MMnmiss_IMpippim_wSid}{  MC sum}");
  MMnmiss_IMpippim_wSid_mc->SetMinimum(1);
  MMnmiss_IMpippim_wSid_mc->SetMaximum(MMnmiss_IMpippim_wSid[0]->GetMaximum());
  MMnmiss_IMpippim_wSid_mc->Draw("colz");

  TH1D* IMpippim_wSid_mc = (TH1D*)MMnmiss_IMpippim_wSid_mc->ProjectionX("IMpippim_wSid_mc");
  IMpippim_wSid_mc->SetLineColor(6);
  TCanvas *cIMpippim_wSid = new TCanvas("cIMpippim_wSid","cIMpippim_wSid");
  cIMpippim_wSid->cd();
  TH1D* IMpippim_wSid[2];
  for(int i=0; i<2; i++) IMpippim_wSid[i] = (TH1D*)MMnmiss_IMpippim_wSid[i]->ProjectionX(Form("IMpippim_wSid_%s",name[i]));
  IMpippim_wSid[0]->Draw("HE");
  IMpippim_wSid_mc->Draw("HEsame");
  for(int i=1; i<2; i++) {
    IMpippim_wSid[i]->SetLineColor(colordef[i]);
  }


  //MMnmiss vs IMpippim w/o K0 w/o (Sid & n);
  TH2D* MMnmiss_IMpippim_woK0_woSid_won_mc = (TH2D*)MMnmiss_IMpippim_woK0_woSid_won[1]->Clone("MMnmiss_IMpippim_woK0_woSid_won_mc");

  TCanvas *cMMnmiss_IMpippim_woK0_woSid_won = new TCanvas("cMMnmiss_IMpippim_woK0_woSid_won","cMMnmiss_IMpippim_woK0_woSid_won",1200,800);
  cMMnmiss_IMpippim_woK0_woSid_won->Divide(2,1);
  cMMnmiss_IMpippim_woK0_woSid_won->cd(1);
  MMnmiss_IMpippim_woK0_woSid_won[0]->SetTitle("#splitline{MMnmiss_IMpippim_woK0_woSid_won}{  Real data}");
  MMnmiss_IMpippim_woK0_woSid_won[0]->SetMinimum(1);
  MMnmiss_IMpippim_woK0_woSid_won[0]->Draw("colz");
  cMMnmiss_IMpippim_woK0_woSid_won->cd(2);
  MMnmiss_IMpippim_woK0_woSid_won_mc->SetTitle("#splitline{MMnmiss_IMpippim_woK0_woSid_won}{  MC sum}");
  MMnmiss_IMpippim_woK0_woSid_won_mc->SetMinimum(1);
  MMnmiss_IMpippim_woK0_woSid_won_mc->SetMaximum(MMnmiss_IMpippim_woK0_woSid_won[0]->GetMaximum());
  MMnmiss_IMpippim_woK0_woSid_won_mc->Draw("colz");

  TH1D* IMpippim_woK0_woSid_won_mc = (TH1D*)MMnmiss_IMpippim_woK0_woSid_won_mc->ProjectionX("IMpippim_woK0_woSid_won_mc");
  IMpippim_woK0_woSid_won_mc->SetLineColor(6);
  TCanvas *cIMpippim_woK0_woSid_won = new TCanvas("cIMpippim_woK0_woSid_won","cIMpippim_woK0_woSid_won");
  cIMpippim_woK0_woSid_won->cd();
  TH1D* IMpippim_woK0_woSid_won[2];
  for(int i=0; i<2; i++) IMpippim_woK0_woSid_won[i] = (TH1D*)MMnmiss_IMpippim_woK0_woSid_won[i]->ProjectionX(Form("IMpippim_woK0_woSid_won_%s",name[i]));
  IMpippim_woK0_woSid_won[0]->Draw("HE");
  IMpippim_woK0_woSid_won_mc->Draw("HEsame");
  for(int i=1; i<2; i++) {
    IMpippim_woK0_woSid_won[i]->SetLineColor(colordef[i]);
    //IMpippim_woK0_woSid_won[i]->Draw("HEsame");
  }

  TCanvas *cIMpippim_woK0_woSid_won_ratio = new TCanvas("cIMpippim_woK0_woSid_won_ratio","cIMpippim_woK0_woSid_won_ratio");
  cIMpippim_woK0_woSid_won_ratio->cd();
  TH1D* IMpippim_woK0_woSid_won_ratio = (TH1D*)IMpippim_woK0_woSid_won[0]->Clone("IMpippim_woK0_woSid_won_ratio");
  IMpippim_woK0_woSid_won_ratio->Divide(IMpippim_woK0_woSid_won_mc);
  IMpippim_woK0_woSid_won_ratio->SetTitle("Data/MC");
  IMpippim_woK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  IMpippim_woK0_woSid_won_ratio->Draw("HE");

  TF1 *evalf_IMpippim = new TF1("evalf_IMpippim","pol6",0,1);
  evalf_IMpippim->SetTitle("IMpippim");
  IMpippim_woK0_woSid_won_ratio->Fit("evalf_IMpippim","","",0.28,0.97);

  //q vs IMnpipi w/o K0 w/o (Sid & n);
  TH2D* q_IMnpipi_woK0_woSid_won_mc = (TH2D*)q_IMnpipi_woK0_woSid_won[1]->Clone("q_IMnpipi_woK0_woSid_won_mc");

  TCanvas *cq_IMnpipi_woK0_woSid_won = new TCanvas("cq_IMnpipi_woK0_woSid_won","q_IMnpipi_woK0_woSid_won",1200,800);
  cq_IMnpipi_woK0_woSid_won->Divide(2,1);
  cq_IMnpipi_woK0_woSid_won->cd(1);
  q_IMnpipi_woK0_woSid_won[0]->SetTitle("#splitline{q_IMnpipi_woK0_woSid_won}{  Real data}");
  q_IMnpipi_woK0_woSid_won[0]->SetMinimum(1);
  q_IMnpipi_woK0_woSid_won[0]->Draw("colz");
  cq_IMnpipi_woK0_woSid_won->cd(2);
  q_IMnpipi_woK0_woSid_won_mc->SetTitle("#splitline{q_IMnpipi_woK0_woSid_won}{  MC sum}");
  q_IMnpipi_woK0_woSid_won_mc->SetMinimum(1);
  q_IMnpipi_woK0_woSid_won_mc->SetMaximum(q_IMnpipi_woK0_woSid_won[0]->GetMaximum());
  q_IMnpipi_woK0_woSid_won_mc->Draw("colz");

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

  TCanvas *cIMnpipi_woK0_woSid_won_ratio = new TCanvas("cIMnpipi_woK0_woSid_won_ratio","cIMnpipi_woK0_woSid_won_ratio");
  cIMnpipi_woK0_woSid_won_ratio->cd();
  TH1D* IMnpipi_woK0_woSid_won_ratio = (TH1D*)IMnpipi_woK0_woSid_won[0]->Clone("IMnpipi_woK0_woSid_won_ratio");
  IMnpipi_woK0_woSid_won_ratio->Divide(IMnpipi_woK0_woSid_won_mc);
  IMnpipi_woK0_woSid_won_ratio->SetTitle("Data/MC");
  IMnpipi_woK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  IMnpipi_woK0_woSid_won_ratio->Draw("HE");

  TF1 *evalf_IMnpipi = new TF1("evalf_IMnpipi","pol5",1.22,2.00);
  evalf_IMnpipi->SetTitle("IMnpipi");
  IMnpipi_woK0_woSid_won_ratio->Fit("evalf_IMnpipi","","",1.22,2.00);


  TH1D* q_woK0_woSid_won_mc = (TH1D*)q_IMnpipi_woK0_woSid_won_mc->ProjectionY("q_woK0_woSid_won_mc");
  q_woK0_woSid_won_mc->SetLineColor(6);
  TCanvas *cq_woK0_woSid_won = new TCanvas("cq_woK0_woSid_won","cq_woK0_woSid_won");
  cq_woK0_woSid_won->cd();
  TH1D* q_woK0_woSid_won[2];
  for(int i=0; i<2; i++) q_woK0_woSid_won[i] = (TH1D*)q_IMnpipi_woK0_woSid_won[i]->ProjectionY(Form("q_woK0_woSid_won_%s",name[i]));
  q_woK0_woSid_won[0]->Draw("HE");
  q_woK0_woSid_won_mc->Draw("HEsame");

  //missing mass vs Mom(pi+pi-) w/o K0 w/o (Sid & n)
  TH2D* MMnmiss_Mompippim_woK0_woSid_won_mc = (TH2D*)MMnmiss_Mompippim_woK0_woSid_won[1]->Clone("MMnmiss_Mompippim_woK0_woSid_won_mc");

  TCanvas *cMMnmiss_Mompippim_woK0_woSid_won = new TCanvas("cMMnmiss_Mompippim_woK0_woSid_won","cMMnmiss_Mompippim_woK0_woSid_won",1200,800);
  cMMnmiss_Mompippim_woK0_woSid_won->Divide(2,1);
  cMMnmiss_Mompippim_woK0_woSid_won->cd(1);
  MMnmiss_Mompippim_woK0_woSid_won[0]->SetTitle("#splitline{MMnmiss_Mompippim_woK0_woSid_won}{  Real data}");
  MMnmiss_Mompippim_woK0_woSid_won[0]->SetMinimum(1);
  MMnmiss_Mompippim_woK0_woSid_won[0]->Draw("colz");
  cMMnmiss_Mompippim_woK0_woSid_won->cd(2);
  MMnmiss_Mompippim_woK0_woSid_won_mc->SetTitle("#splitline{MMnmiss_Mompippim_woK0_woSid_won}{  MC sum}");
  MMnmiss_Mompippim_woK0_woSid_won_mc->SetMinimum(1);
  MMnmiss_Mompippim_woK0_woSid_won_mc->SetMaximum(MMnmiss_Mompippim_woK0_woSid_won[0]->GetMaximum());  
  MMnmiss_Mompippim_woK0_woSid_won_mc->Draw("colz");

  TH1D* Mompippim_woK0_woSid_won_mc = (TH1D*)MMnmiss_Mompippim_woK0_woSid_won_mc->ProjectionX("Mompippim_woK0_woSid_won_mc");
  Mompippim_woK0_woSid_won_mc->SetLineColor(6);
  TCanvas *cMompippim_woK0_woSid_won = new TCanvas("cMompippim_woK0_woSid_won","cMompippim_woK0_woSid_won");
  cMompippim_woK0_woSid_won->cd();
  TH1D* Mompippim_woK0_woSid_won[2];
  for(int i=0; i<2; i++) Mompippim_woK0_woSid_won[i] = (TH1D*)MMnmiss_Mompippim_woK0_woSid_won[i]->ProjectionX(Form("MMnmiss_Mompippim_woK0_woSid_won_%s",name[i]));
  Mompippim_woK0_woSid_won[0]->Draw("HE");
  Mompippim_woK0_woSid_won_mc->Draw("HEsame");

  TCanvas *cMompippim_woK0_woSid_won_ratio = new TCanvas("cMompippim_woK0_woSid_won_ratio","cMompippim_woK0_woSid_won_ratio");
  TH1D* Mompippim_woK0_woSid_won_ratio = (TH1D*)Mompippim_woK0_woSid_won[0]->Clone("Mompippim_woK0_woSid_won_ratio");
  Mompippim_woK0_woSid_won_ratio->Divide(Mompippim_woK0_woSid_won_mc);
  Mompippim_woK0_woSid_won_ratio->SetTitle("Data/MC");
  Mompippim_woK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  Mompippim_woK0_woSid_won_ratio->Draw("HE");

  TF1 *evalf_Mompippim = new TF1("evalf_Mompippim","pol3",0,1);
  evalf_Mompippim->SetTitle("Mompippim");
  Mompippim_woK0_woSid_won_ratio->Fit("evalf_Mompippim","","",0,0.97);
  

  //
  TH2D* IMnpim_IMnpip_woK0_woSid_won_mc = (TH2D*)IMnpim_IMnpip_woK0_woSid_won[1]->Clone("IMnpim_IMnpip_woK0_woSid_won_mc");

  TCanvas* cIMnpim_IMnpip_woK0_woSid_won = new TCanvas("cIMnpim_IMnpip_woK0_woSid_won","cIMnpim_IMnpip_woK0_woSid_won");
  cIMnpim_IMnpip_woK0_woSid_won->Divide(2,1);
  cIMnpim_IMnpip_woK0_woSid_won->cd(1);
  IMnpim_IMnpip_woK0_woSid_won[0]->Draw("colz");
  cIMnpim_IMnpip_woK0_woSid_won->cd(2);
  IMnpim_IMnpip_woK0_woSid_won_mc->SetMinimum(1);
  IMnpim_IMnpip_woK0_woSid_won_mc->SetMaximum(IMnpim_IMnpip_woK0_woSid_won[0]->GetMaximum());
  IMnpim_IMnpip_woK0_woSid_won_mc->Draw("colz");

  //pipmom
  TH2D* pipmom_MMnmiss_woK0_woSid_won_mc = (TH2D*)pipmom_MMnmiss_woK0_woSid_won[1]->Clone("pipmom_MMnmiss_woK0_woSid_won_mc");

  TH1D* pipmom_woK0_woSid_won_mc = (TH1D*)pipmom_MMnmiss_woK0_woSid_won_mc->ProjectionY("pipmom_woK0_woSid_won_mc");
  TH1D* pipmom_woK0_woSid_won[2];
  for(int i=0; i<2; i++) {
    pipmom_woK0_woSid_won[i] = (TH1D*)pipmom_MMnmiss_woK0_woSid_won[i]->ProjectionY(Form("pipmom_woK0_woSid_won_%s",name[i]));
  }
  TCanvas *cpipmom_woK0_woSid_won = new TCanvas("cpipmom_woK0_woSid_won","cpipmom_woK0_woSid_won");
  cpipmom_woK0_woSid_won->cd();
  pipmom_woK0_woSid_won[0]->Draw("HE");
  pipmom_woK0_woSid_won_mc->SetLineColor(6);
  pipmom_woK0_woSid_won_mc->Draw("HEsame");

  TCanvas *cpipmom_woK0_woSid_won_ratio = new TCanvas("cpipmom_woK0_woSid_won_ratio","cpipmom_woK0_woSid_won_ratio");
  cpipmom_woK0_woSid_won_ratio->cd();
  TH1D* pipmom_woK0_woSid_won_ratio = (TH1D*)pipmom_woK0_woSid_won[0]->Clone("pipmom_woK0_woSid_won_ratio");
  pipmom_woK0_woSid_won_ratio->Divide(pipmom_woK0_woSid_won_mc);
  pipmom_woK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  pipmom_woK0_woSid_won_ratio->SetTitle("Data/MC");
  pipmom_woK0_woSid_won_ratio->Draw("HEsame");

  TF1 *evalf_pipmom = new TF1("evalf_pipmom","pol8",0.06,0.7);
  evalf_pipmom->SetTitle("pipmom");
  pipmom_woK0_woSid_won_ratio->Fit("evalf_pipmom","","",0.06,0.7);
  

  TH2D* pipmom_pimmom_woK0_woSid_won_mc = (TH2D*)pipmom_pimmom_woK0_woSid_won[1]->Clone("pipmom_pimmom_woK0_woSid_won_mc");
  TCanvas *cpipmom_pimmom_woK0_woSid_won = new TCanvas("cpipmom_pimmom_woK0_woSid_won","cpipmom_pimmom_woK0_woSid_won");
  cpipmom_pimmom_woK0_woSid_won->Divide(2,1);
  cpipmom_pimmom_woK0_woSid_won->cd(1);
  pipmom_pimmom_woK0_woSid_won[0]->Draw("colz");
  cpipmom_pimmom_woK0_woSid_won->cd(2);
  pipmom_pimmom_woK0_woSid_won_mc->SetMinimum(1);
  pipmom_pimmom_woK0_woSid_won_mc->SetMaximum(pipmom_pimmom_woK0_woSid_won[0]->GetMaximum());
  pipmom_pimmom_woK0_woSid_won_mc->Draw("colz");


  //pimmom
  TH2D* pimmom_MMnmiss_woK0_woSid_won_mc = (TH2D*)pimmom_MMnmiss_woK0_woSid_won[1]->Clone("pimmom_MMnmiss_woK0_woSid_won_mc");

  TH1D* pimmom_woK0_woSid_won_mc = (TH1D*)pimmom_MMnmiss_woK0_woSid_won_mc->ProjectionY("pimmom_woK0_woSid_won_mc");
  TH1D* pimmom_woK0_woSid_won[2];
  for(int i=0; i<2; i++) {
    pimmom_woK0_woSid_won[i] = (TH1D*)pimmom_MMnmiss_woK0_woSid_won[i]->ProjectionY(Form("pimmom_woK0_woSid_won_%s",name[i]));
    pimmom_woK0_woSid_won[i]->SetLineColor(colordef[i]);
  }
  TCanvas *cpimmom_woK0_woSid_won = new TCanvas("cpimmom_woK0_woSid_won","cpimmom_woK0_woSid_won");
  cpimmom_woK0_woSid_won->cd();
  pimmom_woK0_woSid_won[0]->Draw("HE");
  pimmom_woK0_woSid_won_mc->SetLineColor(6);
  pimmom_woK0_woSid_won_mc->Draw("HEsame");

  TCanvas *cpimmom_woK0_woSid_won_ratio = new TCanvas("cpimmom_woK0_woSid_won_ratio","cpimmom_woK0_woSid_won_ratio");
  cpimmom_woK0_woSid_won_ratio->cd();
  TH1D* pimmom_woK0_woSid_won_ratio = (TH1D*)pimmom_woK0_woSid_won[0]->Clone("pimmom_woK0_woSid_won_ratio");
  pimmom_woK0_woSid_won_ratio->Divide(pimmom_woK0_woSid_won_mc);
  pimmom_woK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  pimmom_woK0_woSid_won_ratio->SetTitle("Data/MC");
  pimmom_woK0_woSid_won_ratio->Draw("HEsame");

  TF1 *evalf_pimmom = new TF1("evalf_pimmom","pol8",0.06,0.73);
  evalf_pimmom->SetTitle("pimmom");
  pimmom_woK0_woSid_won_ratio->Fit("evalf_pimmom","","",0.06,0.73);

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

  TCanvas *cnmom_woK0_woSid_won_ratio = new TCanvas("cnmom_woK0_woSid_won_ratio","cnmom_woK0_woSid_won_ratio");
  cnmom_woK0_woSid_won_ratio->cd();
  TH1D* nmom_woK0_woSid_won_ratio = (TH1D*)nmom_woK0_woSid_won[0]->Clone("nmom_woK0_woSid_won_ratio");
  nmom_woK0_woSid_won_ratio->Divide(nmom_woK0_woSid_won_mc);
  nmom_woK0_woSid_won_ratio->SetTitle("Data/MC");
  nmom_woK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  nmom_woK0_woSid_won_ratio->Draw("HEsame");

  TF1 *evalf_nmom = new TF1("evalf_nmom","pol8",0.14,1);
  nmom_woK0_woSid_won_ratio->Fit("evalf_nmom","","",0.14,1);

  
  //weighting to momentum vector 
  TH2D* nmom_cosn_woK0_woSid_won_mc = (TH2D*)nmom_cosn_woK0_woSid_won[1]->Clone("nmom_cosn_woK0_woSid_won_mc");
  TCanvas *cnmom_cosn_woK0_woSid_won = new TCanvas("cnmom_cosn_woK0_woSid_won","cnmom_cosn_woK0_woSid_won",1200,800);
  cnmom_cosn_woK0_woSid_won->Divide(2,1);
  cnmom_cosn_woK0_woSid_won->cd(1);
  nmom_cosn_woK0_woSid_won[0]->SetTitle("#splitline{nmom_cosn_woK0_woSid_won}{Real data}");
  nmom_cosn_woK0_woSid_won[0]->Draw("colz");
  cnmom_cosn_woK0_woSid_won->cd(2);
  nmom_cosn_woK0_woSid_won_mc->SetTitle("#splitline{nmom_cosn_woK0_woSid_won}{MC}");
  nmom_cosn_woK0_woSid_won_mc->SetMinimum(1);
  nmom_cosn_woK0_woSid_won_mc->SetMaximum(nmom_cosn_woK0_woSid_won[0]->GetMaximum());
  nmom_cosn_woK0_woSid_won_mc->Draw("colz");

  TCanvas *ccosn_woK0_woSid_won = new TCanvas("ccosn_woK0_woSid_won","ccosn_woK0_woSid_won");
  ccosn_woK0_woSid_won->cd();
  TH1D* cosn_woK0_woSid_won[2];
  for(int i=0;i<2;i++) cosn_woK0_woSid_won[i] = (TH1D*)nmom_cosn_woK0_woSid_won[i]->ProjectionX(Form("cosn_woK0_woSid_won_%s",name[i]));
  TH1D* cosn_woK0_woSid_won_mc = (TH1D*)nmom_cosn_woK0_woSid_won_mc->ProjectionX("cosn_woK0_woSid_won_mc");
  cosn_woK0_woSid_won_mc->SetLineColor(6);
  cosn_woK0_woSid_won[0]->Draw("HE");
  cosn_woK0_woSid_won_mc->Draw("HEsame");
   

  TCanvas *ccosn_woK0_woSid_won_ratio = new TCanvas("ccosn_woK0_woSid_won_ratio","ccosn_woK0_woSid_won_ratio");
  ccosn_woK0_woSid_won_ratio->cd();
  TH1D* cosn_woK0_woSid_won_ratio = (TH1D*)cosn_woK0_woSid_won[0]->Clone("cosn_woK0_woSid_won_ratio");
  cosn_woK0_woSid_won_ratio->Divide(cosn_woK0_woSid_won_mc);
  cosn_woK0_woSid_won_ratio->SetTitle("Data/MC");
  cosn_woK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  cosn_woK0_woSid_won_ratio->Draw("HE");
  
  TF1 *evalf_cosn_1 = new TF1("evalf_cosn_1","expo",-1,-0.9);
  cosn_woK0_woSid_won_ratio->Fit(evalf_cosn_1,"R","",-1,-0.9);
   
  TF1 *evalf_cosn_2 = new TF1("evalf_cosn_2","pol3",-0.9,0.3);
  cosn_woK0_woSid_won_ratio->Fit(evalf_cosn_2,"R+","",-0.9,0.3);

  Double_t param_cosn[6];
  evalf_cosn_1->GetParameters(&param_cosn[0]);
  evalf_cosn_2->GetParameters(&param_cosn[2]);
  TF1 *evalf_cosn = new TF1("evalf_cosn",func_cosn,-1.00,0.3,6);
  evalf_cosn->SetParameters(param_cosn);
  evalf_cosn->SetLineColor(3);
  evalf_cosn->Draw("same");
  
  TF1 *evalf_cosn_corr = new TF1("evalf_cosn_corr",func_cosn_corr,-1.00,1.0,6);
  //evalf_cosn->SetParameters(param_cosn);
  cosn_woK0_woSid_won_ratio->Fit("evalf_cosn_corr");
  evalf_cosn->SetLineColor(3);
  evalf_cosn->Draw("same");

  //
  TH2D* nmom_cospip_woK0_woSid_won_mc = (TH2D*)nmom_cospip_woK0_woSid_won[1]->Clone("nmom_cospip_woK0_woSid_won_mc");
  TCanvas *cnmom_cospip_woK0_woSid_won = new TCanvas("cnmom_cospip_woK0_woSid_won","cnmom_cospip_woK0_woSid_won",1200,800);
  cnmom_cospip_woK0_woSid_won->Divide(2,1);
  cnmom_cospip_woK0_woSid_won->cd(1);
  nmom_cospip_woK0_woSid_won[0]->SetTitle("#splitline{nmom_cospip_woK0_woSid_won}{Real data}");
  nmom_cospip_woK0_woSid_won[0]->Draw("colz");
  cnmom_cospip_woK0_woSid_won->cd(2);
  nmom_cospip_woK0_woSid_won_mc->SetTitle("#splitline{nmom_cospip_woK0_woSid_won}{MC}");
  nmom_cospip_woK0_woSid_won_mc->SetMinimum(1);
  nmom_cospip_woK0_woSid_won_mc->SetMaximum(nmom_cospip_woK0_woSid_won[0]->GetMaximum());
  nmom_cospip_woK0_woSid_won_mc->Draw("colz");

  TCanvas *ccospip_woK0_woSid_won = new TCanvas("ccospip_woK0_woSid_won","ccospip_woK0_woSid_won");
  ccospip_woK0_woSid_won->cd();
  TH1D* cospip_woK0_woSid_won[2];
  for(int i=0;i<2;i++) cospip_woK0_woSid_won[i] = (TH1D*)nmom_cospip_woK0_woSid_won[i]->ProjectionX(Form("cospip_woK0_woSid_won_%s",name[i]));
  TH1D* cospip_woK0_woSid_won_mc = (TH1D*)nmom_cospip_woK0_woSid_won_mc->ProjectionX("cospip_woK0_woSid_won_mc");
  cospip_woK0_woSid_won_mc->SetLineColor(6);
  cospip_woK0_woSid_won[0]->Draw("HE");
  cospip_woK0_woSid_won_mc->Draw("HEsame");

  TCanvas *ccospip_woK0_woSid_won_ratio = new TCanvas("ccospip_woK0_woSid_won_ratio","ccospip_woK0_woSid_won_ratio");
  ccospip_woK0_woSid_won_ratio->cd();
  TH1D* cospip_woK0_woSid_won_ratio = (TH1D*)cospip_woK0_woSid_won[0]->Clone("cospip_woK0_woSid_won_ratio");
  cospip_woK0_woSid_won_ratio->Divide(cospip_woK0_woSid_won_mc);
  cospip_woK0_woSid_won_ratio->SetTitle("Data/MC");
  cospip_woK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  cospip_woK0_woSid_won_ratio->Draw("HE");
  
  TF1 *evalf_cospip = new TF1("evalf_cospip",func_cospip,-1.0,0.6,6); 
  cospip_woK0_woSid_won_ratio->Fit("evalf_cospip","","",-1.0,0.6);


  //
  TH2D* nmom_cospim_woK0_woSid_won_mc = (TH2D*)nmom_cospim_woK0_woSid_won[1]->Clone("nmom_cospim_woK0_woSid_won_mc");
  TCanvas *cnmom_cospim_woK0_woSid_won = new TCanvas("cnmom_cospim_woK0_woSid_won","cnmom_cospim_woK0_woSid_won",1200,800);
  cnmom_cospim_woK0_woSid_won->Divide(2,1);
  cnmom_cospim_woK0_woSid_won->cd(1);
  nmom_cospim_woK0_woSid_won[0]->SetTitle("#splitline{nmom_cospim_woK0_woSid_won}{Real data}");
  nmom_cospim_woK0_woSid_won[0]->Draw("colz");
  cnmom_cospim_woK0_woSid_won->cd(2);
  nmom_cospim_woK0_woSid_won_mc->SetTitle("#splitline{nmom_cospim_woK0_woSid_won}{MC}");
  nmom_cospim_woK0_woSid_won_mc->SetMinimum(1);
  nmom_cospim_woK0_woSid_won_mc->SetMaximum(nmom_cospim_woK0_woSid_won[0]->GetMaximum());
  nmom_cospim_woK0_woSid_won_mc->Draw("colz");

  TCanvas *ccospim_woK0_woSid_won = new TCanvas("ccospim_woK0_woSid_won","ccospim_woK0_woSid_won");
  ccospim_woK0_woSid_won->cd();
  TH1D* cospim_woK0_woSid_won[2];
  for(int i=0;i<2;i++) cospim_woK0_woSid_won[i] = (TH1D*)nmom_cospim_woK0_woSid_won[i]->ProjectionX(Form("cospim_woK0_woSid_won_%s",name[i]));
  TH1D* cospim_woK0_woSid_won_mc = (TH1D*)nmom_cospim_woK0_woSid_won_mc->ProjectionX("cospim_woK0_woSid_won_mc");
  cospim_woK0_woSid_won_mc->SetLineColor(6);
  cospim_woK0_woSid_won[0]->Draw("HE");
  cospim_woK0_woSid_won_mc->Draw("HEsame");

  TCanvas *ccospim_woK0_woSid_won_ratio = new TCanvas("ccospim_woK0_woSid_won_ratio","ccospim_woK0_woSid_won_ratio");
  ccospim_woK0_woSid_won_ratio->cd();
  TH1D* cospim_woK0_woSid_won_ratio = (TH1D*)cospim_woK0_woSid_won[0]->Clone("cospim_woK0_woSid_won_ratio");
  cospim_woK0_woSid_won_ratio->Divide(cospim_woK0_woSid_won_mc);
  cospim_woK0_woSid_won_ratio->SetTitle("Data/MC");
  cospim_woK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  cospim_woK0_woSid_won_ratio->Draw("HE");
   
  TF1 *evalf_cospim = new TF1("evalf_cospim","pol8",-0.92,0.50);
  cospim_woK0_woSid_won_ratio->Fit(evalf_cospim,"","",-0.92,0.50);
 
  //
  TH2D* nmom_phinpip_woK0_woSid_won_mc = (TH2D*)nmom_phinpip_woK0_woSid_won[1]->Clone("nmom_phinpip_woK0_woSid_won_mc");
  TCanvas *cnmom_phinpip_woK0_woSid_won = new TCanvas("cnmom_phinpip_woK0_woSid_won","cnmom_phinpip_woK0_woSid_won",1200,800);
  cnmom_phinpip_woK0_woSid_won->Divide(2,1);
  cnmom_phinpip_woK0_woSid_won->cd(1);
  nmom_phinpip_woK0_woSid_won[0]->SetTitle("#splitline{nmom_phinpip_woK0_woSid_won}{Real data}");
  nmom_phinpip_woK0_woSid_won[0]->Draw("colz");
  cnmom_phinpip_woK0_woSid_won->cd(2);
  nmom_phinpip_woK0_woSid_won_mc->SetTitle("#splitline{nmom_phinpip_woK0_woSid_won}{MC}");
  nmom_phinpip_woK0_woSid_won_mc->SetMaximum(nmom_phinpip_woK0_woSid_won[0]->GetMaximum());
  nmom_phinpip_woK0_woSid_won_mc->SetMinimum(1);
  nmom_phinpip_woK0_woSid_won_mc->Draw("colz");

  TCanvas *cphinpip_woK0_woSid_won = new TCanvas("cphinpip_woK0_woSid_won","cphinpip_woK0_woSid_won");
  cphinpip_woK0_woSid_won->cd();
  TH1D* phinpip_woK0_woSid_won[2];
  for(int i=0;i<2;i++) phinpip_woK0_woSid_won[i] = (TH1D*)nmom_phinpip_woK0_woSid_won[i]->ProjectionX(Form("phinpip_woK0_woSid_won_%s",name[i]));
  TH1D* phinpip_woK0_woSid_won_mc = (TH1D*)nmom_phinpip_woK0_woSid_won_mc->ProjectionX("phinpip_woK0_woSid_won_mc");
  phinpip_woK0_woSid_won_mc->SetLineColor(6);
  phinpip_woK0_woSid_won[0]->Draw("HE");
  phinpip_woK0_woSid_won_mc->Draw("HEsame");

  TCanvas *cphinpip_woK0_woSid_won_ratio = new TCanvas("cphinpip_woK0_woSid_won_ratio","cphinpip_woK0_woSid_won_ratio");
  cphinpip_woK0_woSid_won_ratio->cd();
  TH1D* phinpip_woK0_woSid_won_ratio = (TH1D*)phinpip_woK0_woSid_won[0]->Clone("phinpip_woK0_woSid_won_ratio");
  phinpip_woK0_woSid_won_ratio->Divide(phinpip_woK0_woSid_won_mc);
  phinpip_woK0_woSid_won_ratio->SetTitle("Data/MC");
  phinpip_woK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  phinpip_woK0_woSid_won_ratio->Draw("HE");
  
  TH2D* nmom_phipip_woK0_woSid_won_mc = (TH2D*)nmom_phipip_woK0_woSid_won[1]->Clone("nmom_phipip_woK0_woSid_won_mc");
  TCanvas *cnmom_phipip_woK0_woSid_won = new TCanvas("cnmom_phipip_woK0_woSid_won","cnmom_phipip_woK0_woSid_won",1200,800);
  cnmom_phipip_woK0_woSid_won->Divide(2,1);
  cnmom_phipip_woK0_woSid_won->cd(1);
  nmom_phipip_woK0_woSid_won[0]->SetTitle("#splitline{nmom_phipip_woK0_woSid_won}{Real data}");
  nmom_phipip_woK0_woSid_won[0]->Draw("colz");
  cnmom_phipip_woK0_woSid_won->cd(2);
  nmom_phipip_woK0_woSid_won_mc->SetTitle("#splitline{nmom_phipip_woK0_woSid_won}{MC}");
  nmom_phipip_woK0_woSid_won_mc->SetMaximum(nmom_phipip_woK0_woSid_won[0]->GetMaximum());
  nmom_phipip_woK0_woSid_won_mc->SetMinimum(1);
  nmom_phipip_woK0_woSid_won_mc->Draw("colz");

  TCanvas *cphipip_woK0_woSid_won = new TCanvas("cphipip_woK0_woSid_won","cphipip_woK0_woSid_won");
  cphipip_woK0_woSid_won->cd();
  TH1D* phipip_woK0_woSid_won[2];
  for(int i=0;i<2;i++) phipip_woK0_woSid_won[i] = (TH1D*)nmom_phipip_woK0_woSid_won[i]->ProjectionX(Form("phipip_woK0_woSid_won_%s",name[i]));
  TH1D* phipip_woK0_woSid_won_mc = (TH1D*)nmom_phipip_woK0_woSid_won_mc->ProjectionX("phipip_woK0_woSid_won_mc");
  phipip_woK0_woSid_won_mc->SetLineColor(6);
  phipip_woK0_woSid_won[0]->Draw("HE");
  phipip_woK0_woSid_won_mc->Draw("HEsame");

  TCanvas *cphipip_woK0_woSid_won_ratio = new TCanvas("cphipip_woK0_woSid_won_ratio","cphipip_woK0_woSid_won_ratio");
  cphipip_woK0_woSid_won_ratio->cd();
  TH1D* phipip_woK0_woSid_won_ratio = (TH1D*)phipip_woK0_woSid_won[0]->Clone("phipip_woK0_woSid_won_ratio");
  phipip_woK0_woSid_won_ratio->Divide(phipip_woK0_woSid_won_mc);
  phipip_woK0_woSid_won_ratio->SetTitle("Data/MC");
  phipip_woK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  phipip_woK0_woSid_won_ratio->Draw("HE");
   

  //
  TH2D* nmom_phinpim_woK0_woSid_won_mc = (TH2D*)nmom_phinpim_woK0_woSid_won[1]->Clone("nmom_phinpim_woK0_woSid_won_mc");
  TCanvas *cnmom_phinpim_woK0_woSid_won = new TCanvas("cnmom_phinpim_woK0_woSid_won","cnmom_phinpim_woK0_woSid_won",1200,800);
  cnmom_phinpim_woK0_woSid_won->Divide(2,1);
  cnmom_phinpim_woK0_woSid_won->cd(1);
  nmom_phinpim_woK0_woSid_won[0]->SetTitle("#splitline{nmom_phinpim_woK0_woSid_won}{Real data}");
  nmom_phinpim_woK0_woSid_won[0]->Draw("colz");
  cnmom_phinpim_woK0_woSid_won->cd(2);
  nmom_phinpim_woK0_woSid_won_mc->SetTitle("#splitline{nmom_phinpim_woK0_woSid_won}{MC}");
  nmom_phinpim_woK0_woSid_won_mc->SetMinimum(1);
  nmom_phinpim_woK0_woSid_won_mc->SetMaximum(nmom_phinpim_woK0_woSid_won[0]->GetMaximum());
  nmom_phinpim_woK0_woSid_won_mc->Draw("colz");

  TCanvas *cphinpim_woK0_woSid_won = new TCanvas("cphinpim_woK0_woSid_won","cphinpim_woK0_woSid_won");
  cphinpim_woK0_woSid_won->cd();
  TH1D* phinpim_woK0_woSid_won[2];
  for(int i=0;i<2;i++) phinpim_woK0_woSid_won[i] = (TH1D*)nmom_phinpim_woK0_woSid_won[i]->ProjectionX(Form("phinpim_woK0_woSid_won_%s",name[i]));
  TH1D* phinpim_woK0_woSid_won_mc = (TH1D*)nmom_phinpim_woK0_woSid_won_mc->ProjectionX("phinpim_woK0_woSid_won_mc");
  phinpim_woK0_woSid_won_mc->SetLineColor(6);
  phinpim_woK0_woSid_won[0]->Draw("HE");
  phinpim_woK0_woSid_won_mc->Draw("HEsame");

  TCanvas *cphinpim_woK0_woSid_won_ratio = new TCanvas("cphinpim_woK0_woSid_won_ratio","cphinpim_woK0_woSid_won_ratio");
  cphinpim_woK0_woSid_won_ratio->cd();
  TH1D* phinpim_woK0_woSid_won_ratio = (TH1D*)phinpim_woK0_woSid_won[0]->Clone("phinpim_woK0_woSid_won_ratio");
  phinpim_woK0_woSid_won_ratio->Divide(phinpim_woK0_woSid_won_mc);
  phinpim_woK0_woSid_won_ratio->SetTitle("Data/MC");
  phinpim_woK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  phinpim_woK0_woSid_won_ratio->Draw("HE");

  TF1 *evalf_phinpim = new TF1("evalf_phinpim",func_phinpim,-1.0*TMath::Pi(),TMath::Pi(),10);
  phinpim_woK0_woSid_won_ratio->Fit(evalf_phinpim);

  TH2D* nmom_phipim_woK0_woSid_won_mc = (TH2D*)nmom_phipim_woK0_woSid_won[1]->Clone("nmom_phipim_woK0_woSid_won_mc");
  TCanvas *cnmom_phipim_woK0_woSid_won = new TCanvas("cnmom_phipim_woK0_woSid_won","cnmom_phipim_woK0_woSid_won",1200,800);
  cnmom_phipim_woK0_woSid_won->Divide(2,1);
  cnmom_phipim_woK0_woSid_won->cd(1);
  nmom_phipim_woK0_woSid_won[0]->SetTitle("#splitline{nmom_phipim_woK0_woSid_won}{Real data}");
  nmom_phipim_woK0_woSid_won[0]->Draw("colz");
  cnmom_phipim_woK0_woSid_won->cd(2);
  nmom_phipim_woK0_woSid_won_mc->SetTitle("#splitline{nmom_phipim_woK0_woSid_won}{MC}");
  nmom_phipim_woK0_woSid_won_mc->SetMinimum(1);
  nmom_phipim_woK0_woSid_won_mc->SetMaximum(nmom_phipim_woK0_woSid_won[0]->GetMaximum());
  nmom_phipim_woK0_woSid_won_mc->Draw("colz");

  TCanvas *cphipim_woK0_woSid_won = new TCanvas("cphipim_woK0_woSid_won","cphipim_woK0_woSid_won");
  cphipim_woK0_woSid_won->cd();
  TH1D* phipim_woK0_woSid_won[2];
  for(int i=0;i<2;i++) phipim_woK0_woSid_won[i] = (TH1D*)nmom_phipim_woK0_woSid_won[i]->ProjectionX(Form("phipim_woK0_woSid_won_%s",name[i]));
  TH1D* phipim_woK0_woSid_won_mc = (TH1D*)nmom_phipim_woK0_woSid_won_mc->ProjectionX("phipim_woK0_woSid_won_mc");
  phipim_woK0_woSid_won_mc->SetLineColor(6);
  phipim_woK0_woSid_won[0]->Draw("HE");
  phipim_woK0_woSid_won_mc->Draw("HEsame");

  TCanvas *cphipim_woK0_woSid_won_ratio = new TCanvas("cphipim_woK0_woSid_won_ratio","cphipim_woK0_woSid_won_ratio");
  cphipim_woK0_woSid_won_ratio->cd();
  TH1D* phipim_woK0_woSid_won_ratio = (TH1D*)phipim_woK0_woSid_won[0]->Clone("phipim_woK0_woSid_won_ratio");
  phipim_woK0_woSid_won_ratio->Divide(phipim_woK0_woSid_won_mc);
  phipim_woK0_woSid_won_ratio->SetTitle("Data/MC");
  phipim_woK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  phipim_woK0_woSid_won_ratio->Draw("HE");


  TH2D* nmom_phin_woK0_woSid_won_mc = (TH2D*)nmom_phin_woK0_woSid_won[1]->Clone("nmom_phin_woK0_woSid_won_mc");
  TCanvas *cnmom_phin_woK0_woSid_won = new TCanvas("cnmom_phin_woK0_woSid_won","cnmom_phin_woK0_woSid_won",1200,800);
  cnmom_phin_woK0_woSid_won->Divide(2,1);
  cnmom_phin_woK0_woSid_won->cd(1);
  nmom_phin_woK0_woSid_won[0]->SetTitle("#splitline{nmom_phin_woK0_woSid_won}{Real data}");
  nmom_phin_woK0_woSid_won[0]->Draw("colz");
  cnmom_phin_woK0_woSid_won->cd(2);
  nmom_phin_woK0_woSid_won_mc->SetTitle("#splitline{nmom_phin_woK0_woSid_won}{MC}");
  nmom_phin_woK0_woSid_won_mc->SetMinimum(1);
  nmom_phin_woK0_woSid_won_mc->SetMaximum(nmom_phin_woK0_woSid_won[0]->GetMaximum());
  nmom_phin_woK0_woSid_won_mc->Draw("colz");

  TCanvas *cphin_woK0_woSid_won = new TCanvas("cphin_woK0_woSid_won","cphin_woK0_woSid_won");
  cphin_woK0_woSid_won->cd();
  TH1D* phin_woK0_woSid_won[2];
  for(int i=0;i<2;i++) phin_woK0_woSid_won[i] = (TH1D*)nmom_phin_woK0_woSid_won[i]->ProjectionX(Form("phin_woK0_woSid_won_%s",name[i]));
  TH1D* phin_woK0_woSid_won_mc = (TH1D*)nmom_phin_woK0_woSid_won_mc->ProjectionX("phin_woK0_woSid_won_mc");
  phin_woK0_woSid_won_mc->SetLineColor(6);
  phin_woK0_woSid_won[0]->Draw("HE");
  phin_woK0_woSid_won_mc->Draw("HEsame");

  TCanvas *cphin_woK0_woSid_won_ratio = new TCanvas("cphin_woK0_woSid_won_ratio","cphin_woK0_woSid_won_ratio");
  cphin_woK0_woSid_won_ratio->cd();
  TH1D* phin_woK0_woSid_won_ratio = (TH1D*)phin_woK0_woSid_won[0]->Clone("phin_woK0_woSid_won_ratio");
  phin_woK0_woSid_won_ratio->Divide(phin_woK0_woSid_won_mc);
  phin_woK0_woSid_won_ratio->SetTitle("Data/MC");
  phin_woK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  phin_woK0_woSid_won_ratio->Draw("HE");


  ////////////////////////////
  // w/K0 + fake neutron
  ///////////////////////////

  //w K0, w/o ((Sp or Sm) & missing neutron)
  TH2D *MMnmiss_IMnpip_wK0_woSid_won_mc = (TH2D*)MMnmiss_IMnpip_wK0_woSid_won[1]->Clone("MMnmiss_IMnpip_wK0_woSid_won_mc");

  TCanvas *cMMnmiss_IMnpip_wK0_woSid_won = new TCanvas("cMMnmiss_IMnpip_wK0_woSid_won","cMMnmiss_IMnpip_wK0_woSid_won",1200,800);
  cMMnmiss_IMnpip_wK0_woSid_won->Divide(2,1);
  cMMnmiss_IMnpip_wK0_woSid_won->cd(1);
  MMnmiss_IMnpip_wK0_woSid_won_rdata->SetTitle("#splitline{MMnmiss_IMnpip_wK0_woSid_won}{  Real data}");
  MMnmiss_IMnpip_wK0_woSid_won_rdata->Draw("colz");
  cMMnmiss_IMnpip_wK0_woSid_won->cd(2);
  MMnmiss_IMnpip_wK0_woSid_won_mc->SetTitle("#splitline{MMnmiss_IMnpip_wK0_woSid_won}{  MC sum}");
  MMnmiss_IMnpip_wK0_woSid_won_mc->SetMinimum(1);
  MMnmiss_IMnpip_wK0_woSid_won_mc->SetMaximum(MMnmiss_IMnpip_wK0_woSid_won_rdata->GetMaximum());
  MMnmiss_IMnpip_wK0_woSid_won_mc->Draw("colz");

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
  TCanvas *cMMnmiss_wK0_woSid_won_ratio = new TCanvas("cMMnmiss_wK0_woSid_won_ratio","cMMnmiss_wK0_woSid_won_ratio");
  TH1D* MMnmiss_wK0_woSid_won_ratio = (TH1D*)MMnmiss_wK0_woSid_won[0]->Clone("MMnmiss_wK0_woSid_won_ratio");
  MMnmiss_wK0_woSid_won_ratio->Divide(MMnmiss_wK0_woSid_won_mc);
  MMnmiss_wK0_woSid_won_ratio->SetTitle("Data/MC");
  MMnmiss_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  MMnmiss_wK0_woSid_won_ratio->Draw("HE");


  TF1 *evalf_MMnmiss_wK0 = new TF1("evalf_MMnmiss_wK0",func_MMnmiss_wK0,0,1.5,10);
  MMnmiss_wK0_woSid_won_ratio->Fit("evalf_MMnmiss_wK0");
  evalf_MMnmiss_wK0->SetLineColor(4);
  evalf_MMnmiss_wK0->Draw("same");
  evalf_MMnmiss_wK0->Print();


  TH2D *MMom_MMass_wK0_woSid_won_mc = (TH2D*)MMom_MMass_wK0_woSid_won[1]->Clone("MMom_MMass_wK0_woSid_won_mc");
  TCanvas *cMMom_MMass_wK0_woSid_won = new TCanvas("cMMom_MMass_wK0_woSid_won","cMMom_MMass_wK0_woSid_won");
  cMMom_MMass_wK0_woSid_won->Divide(2,1);
  cMMom_MMass_wK0_woSid_won->cd(1);
  MMom_MMass_wK0_woSid_won[0]->Draw("colz");
  cMMom_MMass_wK0_woSid_won->cd(2);
  MMom_MMass_wK0_woSid_won_mc->SetMinimum(1);
  MMom_MMass_wK0_woSid_won_mc->SetMaximum(MMom_MMass_wK0_woSid_won[0]->GetMaximum());
  MMom_MMass_wK0_woSid_won_mc->Draw("colz");

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


  TCanvas *cMMom_wK0_woSid_won_ratio = new TCanvas("cMMom_wK0_woSid_won_ratio","cMMom_wK0_woSid_won_ratio");
  cMMom_wK0_woSid_won_ratio->cd();
  TH1D* MMom_wK0_woSid_won_ratio = (TH1D*)MMom_wK0_woSid_won[0]->Clone("MMom_wK0_woSid_won_ratio");
  MMom_wK0_woSid_won_ratio->Divide(MMom_wK0_woSid_won_mc);
  MMom_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  MMom_wK0_woSid_won_ratio->SetTitle("Data/MC");
  MMom_wK0_woSid_won_ratio->Draw("HE");

  Double_t param_MMom_wK0[5];
  TF1 *evalf_MMom_wK0 = new TF1("evalf_MMom_wK0","pol5",0.4,1.5);
  evalf_MMom_wK0->SetLineColor(4);
  MMom_wK0_woSid_won_ratio->Fit("evalf_MMom_wK0","","",0.4,1.5);


  //projection to IMnpip (miss n & Sigma+/-)
  TCanvas *cIMnpip_wK0_woSid_won = new TCanvas("cIMnpip_wK0_woSid_won","cIMnpip_wK0_woSid_won");
  cIMnpip_wK0_woSid_won->cd();
  TH1D* IMnpip_wK0_woSid_won[2];
  for(int i=0; i<2; i++)IMnpip_wK0_woSid_won[i] = (TH1D*)MMnmiss_IMnpip_wK0_woSid_won[i]->ProjectionX(Form("IMnpip_wK0_woSid_won_%s",name[i]));
  IMnpip_wK0_woSid_won[0]->Draw("HE");//rdata
  TH1D* IMnpip_wK0_woSid_won_mc = (TH1D*)MMnmiss_IMnpip_wK0_woSid_won_mc->ProjectionX("IMnpip_wK0_woSid_won_mc");
  IMnpip_wK0_woSid_won_mc->SetLineColor(6);
  IMnpip_wK0_woSid_won_mc->Draw("HEsame");

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

  TCanvas *cMMnmiss_IMnpim_wK0_woSid_won = new TCanvas("cMMnmiss_IMnpim_wK0_woSid_won","cMMnmiss_IMnpim_wK0_woSid_won",1200,800);
  cMMnmiss_IMnpim_wK0_woSid_won->Divide(2,1);
  cMMnmiss_IMnpim_wK0_woSid_won->cd(1);
  MMnmiss_IMnpim_wK0_woSid_won_rdata->SetTitle("#splitline{MMnmiss_IMnpim_wK0_woSid_won}{  Real data}");
  MMnmiss_IMnpim_wK0_woSid_won_rdata->Draw("colz");
  cMMnmiss_IMnpim_wK0_woSid_won->cd(2);
  MMnmiss_IMnpim_wK0_woSid_won_mc->SetTitle("#splitline{MMnmiss_IMnpim_wK0_woSid_won}{  MC sum}");
  MMnmiss_IMnpim_wK0_woSid_won_mc->SetMinimum(1);
  MMnmiss_IMnpim_wK0_woSid_won_mc->SetMaximum(MMnmiss_IMnpim_wK0_woSid_won_rdata->GetMaximum());
  MMnmiss_IMnpim_wK0_woSid_won_mc->Draw("colz");

  TCanvas *cIMnpim_wK0_woSid_won = new TCanvas("IMnpim_wK0_woSid_won","IMnpim_wK0_woSid_won");
  cIMnpim_wK0_woSid_won->cd();
  TH1D* IMnpim_wK0_woSid_won[2];
  for(int i=0; i<2; i++)IMnpim_wK0_woSid_won[i] = (TH1D*)MMnmiss_IMnpim_wK0_woSid_won[i]->ProjectionX(Form("IMnpim_wK0_woSid_won_%s",name[i]));
  TH1D* IMnpim_wK0_woSid_won_mc = MMnmiss_IMnpim_wK0_woSid_won_mc->ProjectionX("IMnpim_wK0_woSid_won_mc");
  IMnpim_wK0_woSid_won[0]->Draw("HE");
  IMnpim_wK0_woSid_won_mc->SetLineColor(6);
  IMnpim_wK0_woSid_won_mc->Draw("HEsame");

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

  TCanvas *cMMnmiss_IMpippim_wK0_woSid_won = new TCanvas("cMMnmiss_IMpippim_wK0_woSid_won","cMMnmiss_IMpippim_wK0_woSid_won",1200,800);
  cMMnmiss_IMpippim_wK0_woSid_won->Divide(2,1);
  cMMnmiss_IMpippim_wK0_woSid_won->cd(1);
  MMnmiss_IMpippim_wK0_woSid_won[0]->SetTitle("#splitline{MMnmiss_IMpippim_wK0_woSid_won}{  Real data}");
  MMnmiss_IMpippim_wK0_woSid_won[0]->Draw("colz");
  cMMnmiss_IMpippim_wK0_woSid_won->cd(2);
  MMnmiss_IMpippim_wK0_woSid_won_mc->SetTitle("#splitline{MMnmiss_IMpippim_wK0_woSid_won}{  MC sum}");
  MMnmiss_IMpippim_wK0_woSid_won_mc->SetMinimum(1);
  MMnmiss_IMpippim_wK0_woSid_won_mc->SetMaximum(MMnmiss_IMpippim_wK0_woSid_won[0]->GetMaximum());
  MMnmiss_IMpippim_wK0_woSid_won_mc->Draw("colz");

  TH1D* IMpippim_wK0_woSid_won_mc = (TH1D*)MMnmiss_IMpippim_wK0_woSid_won_mc->ProjectionX("IMpippim_wK0_woSid_won_mc");
  IMpippim_wK0_woSid_won_mc->SetLineColor(6);
  TCanvas *cIMpippim_wK0_woSid_won = new TCanvas("cIMpippim_wK0_woSid_won","cIMpippim_wK0_woSid_won");
  cIMpippim_wK0_woSid_won->cd();
  TH1D* IMpippim_wK0_woSid_won[2];
  for(int i=0; i<2; i++) IMpippim_wK0_woSid_won[i] = (TH1D*)MMnmiss_IMpippim_wK0_woSid_won[i]->ProjectionX(Form("IMpippim_wK0_woSid_won_%s",name[i]));
  IMpippim_wK0_woSid_won[0]->Draw("HE");
  IMpippim_wK0_woSid_won_mc->Draw("HEsame");

  TCanvas *cIMpippim_wK0_woSid_won_ratio = new TCanvas("cIMpippim_wK0_woSid_won_ratio","cIMpippim_wK0_woSid_won_ratio");
  cIMpippim_wK0_woSid_won_ratio->cd();
  TH1D* IMpippim_wK0_woSid_won_ratio = (TH1D*)IMpippim_wK0_woSid_won[0]->Clone("IMpippim_wK0_woSid_won_ratio");
  IMpippim_wK0_woSid_won_ratio->Divide(IMpippim_wK0_woSid_won_mc);
  IMpippim_wK0_woSid_won_ratio->SetTitle("Data/MC");
  IMpippim_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  IMpippim_wK0_woSid_won_ratio->Draw("HE");

  TF1 *evalf_IMpippim_wK0 = new TF1("evalf_IMpippim_wK0","pol6",0,1);
  IMpippim_wK0_woSid_won_ratio->Fit("evalf_IMpippim_wK0","","",0.28,0.97);

  //q vs IMnpipi w/o K0 w/o (Sid & n);
  TH2D* q_IMnpipi_wK0_woSid_won_mc = (TH2D*)q_IMnpipi_wK0_woSid_won[1]->Clone("q_IMnpipi_wK0_woSid_won_mc");

  TCanvas *cq_IMnpipi_wK0_woSid_won = new TCanvas("cq_IMnpipi_wK0_woSid_won","q_IMnpipi_wK0_woSid_won",1200,800);
  cq_IMnpipi_wK0_woSid_won->Divide(2,1);
  cq_IMnpipi_wK0_woSid_won->cd(1);
  q_IMnpipi_wK0_woSid_won[0]->SetTitle("#splitline{q_IMnpipi_wK0_woSid_won}{  Real data}");
  q_IMnpipi_wK0_woSid_won[0]->Draw("colz");
  cq_IMnpipi_wK0_woSid_won->cd(2);
  q_IMnpipi_wK0_woSid_won_mc->SetTitle("#splitline{q_IMnpipi_wK0_woSid_won}{  MC sum}");
  q_IMnpipi_wK0_woSid_won_mc->SetMinimum(1);
  q_IMnpipi_wK0_woSid_won_mc->SetMaximum(q_IMnpipi_wK0_woSid_won[0]->GetMaximum());
  q_IMnpipi_wK0_woSid_won_mc->Draw("colz");

  TH1D* IMnpipi_wK0_woSid_won_mc = (TH1D*)q_IMnpipi_wK0_woSid_won_mc->ProjectionX("IMnpipi_wK0_woSid_won_mc");
  IMnpipi_wK0_woSid_won_mc->SetLineColor(6);
  TCanvas *cIMnpipi_wK0_woSid_won = new TCanvas("cIMnpipi_wK0_woSid_won","cIMnpipi_wK0_woSid_won");
  cIMnpipi_wK0_woSid_won->cd();
  TH1D* IMnpipi_wK0_woSid_won[2];
  for(int i=0; i<2; i++) IMnpipi_wK0_woSid_won[i] = (TH1D*)q_IMnpipi_wK0_woSid_won[i]->ProjectionX(Form("IMnpipi_wK0_woSid_won_%s",name[i]));
  IMnpipi_wK0_woSid_won[0]->Draw("HE");
  IMnpipi_wK0_woSid_won_mc->Draw("HEsame");

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
  TH1D* q_wK0_woSid_won[2];
  for(int i=0; i<2; i++) q_wK0_woSid_won[i] = (TH1D*)q_IMnpipi_wK0_woSid_won[i]->ProjectionY(Form("q_wK0_woSid_won_%s",name[i]));
  q_wK0_woSid_won[0]->Draw("HE");
  q_wK0_woSid_won_mc->Draw("HEsame");

  //missing mass vs Mom(pi+pi-) w/o K0 w/o (Sid & n)
  TH2D* MMnmiss_Mompippim_wK0_woSid_won_mc = (TH2D*)MMnmiss_Mompippim_wK0_woSid_won[1]->Clone("MMnmiss_Mompippim_wK0_woSid_won_mc");

  TCanvas *cMMnmiss_Mompippim_wK0_woSid_won = new TCanvas("cMMnmiss_Mompippim_wK0_woSid_won","cMMnmiss_Mompippim_wK0_woSid_won",1200,800);
  cMMnmiss_Mompippim_wK0_woSid_won->Divide(2,1);
  cMMnmiss_Mompippim_wK0_woSid_won->cd(1);
  MMnmiss_Mompippim_wK0_woSid_won[0]->SetTitle("#splitline{MMnmiss_Mompippim_wK0_woSid_won}{  Real data}");
  MMnmiss_Mompippim_wK0_woSid_won[0]->Draw("colz");
  cMMnmiss_Mompippim_wK0_woSid_won->cd(2);
  MMnmiss_Mompippim_wK0_woSid_won_mc->SetTitle("#splitline{MMnmiss_Mompippim_wK0_woSid_won}{  MC sum}");
  MMnmiss_Mompippim_wK0_woSid_won_mc->SetMinimum(1);
  MMnmiss_Mompippim_wK0_woSid_won_mc->SetMaximum(MMnmiss_Mompippim_wK0_woSid_won[0]->GetMaximum());
  MMnmiss_Mompippim_wK0_woSid_won_mc->Draw("colz");

  TH1D* Mompippim_wK0_woSid_won_mc = (TH1D*)MMnmiss_Mompippim_wK0_woSid_won_mc->ProjectionX("Mompippim_wK0_woSid_won_mc");
  Mompippim_wK0_woSid_won_mc->SetLineColor(6);
  TCanvas *cMompippim_wK0_woSid_won = new TCanvas("cMompippim_wK0_woSid_won","cMompippim_wK0_woSid_won");
  cMompippim_wK0_woSid_won->cd();
  TH1D* Mompippim_wK0_woSid_won[2];
  for(int i=0; i<2; i++) Mompippim_wK0_woSid_won[i] = (TH1D*)MMnmiss_Mompippim_wK0_woSid_won[i]->ProjectionX(Form("MMnmiss_Mompippim_wK0_woSid_won_%s",name[i]));
  Mompippim_wK0_woSid_won[0]->Draw("HE");

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

  TH1D* pipmom_wK0_woSid_won_mc = (TH1D*)pipmom_MMnmiss_wK0_woSid_won_mc->ProjectionY("pipmom_wK0_woSid_won_mc");
  TH1D* pipmom_wK0_woSid_won[2];
  for(int i=0; i<2; i++) {
    pipmom_wK0_woSid_won[i] = (TH1D*)pipmom_MMnmiss_wK0_woSid_won[i]->ProjectionY(Form("pipmom_wK0_woSid_won_%s",name[i]));
    pipmom_wK0_woSid_won[i]->SetLineColor(colordef[i]);
  }
  TCanvas *cpipmom_wK0_woSid_won = new TCanvas("cpipmom_wK0_woSid_won","cpipmom_wK0_woSid_won");
  cpipmom_wK0_woSid_won->cd();
  pipmom_wK0_woSid_won[0]->Draw("HE");
  pipmom_wK0_woSid_won_mc->SetLineColor(6);
  pipmom_wK0_woSid_won_mc->Draw("HEsame");

  TCanvas *cpipmom_wK0_woSid_won_ratio = new TCanvas("cpipmom_wK0_woSid_won_ratio","cpipmom_wK0_woSid_won_ratio");
  cpipmom_wK0_woSid_won_ratio->cd();
  TH1D* pipmom_wK0_woSid_won_ratio = (TH1D*)pipmom_wK0_woSid_won[0]->Clone("pipmom_wK0_woSid_won_ratio");
  pipmom_wK0_woSid_won_ratio->Divide(pipmom_wK0_woSid_won_mc);
  pipmom_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  pipmom_wK0_woSid_won_ratio->SetTitle("Data/MC");
  pipmom_wK0_woSid_won_ratio->Draw("HEsame");

  TF1 *evalf_pipmom_wK0 = new TF1("evalf_pipmom_wK0","pol8",0.06,0.75);
  pipmom_wK0_woSid_won_ratio->Fit("evalf_pipmom_wK0","","",0.06,0.75);
  
  //
  TH2D* pipmom_pimmom_wK0_woSid_won_mc = (TH2D*)pipmom_pimmom_wK0_woSid_won[1]->Clone("pipmom_pimmom_wK0_woSid_won_mc");
  TCanvas *cpipmom_pimmom_wK0_woSid_won = new TCanvas("cpipmom_pimmom_wK0_woSid_won","cpipmom_pimmom_wK0_woSid_won");
  cpipmom_pimmom_wK0_woSid_won->Divide(2,1);
  cpipmom_pimmom_wK0_woSid_won->cd(1);
  pipmom_pimmom_wK0_woSid_won[0]->Draw("colz");
  cpipmom_pimmom_wK0_woSid_won->cd(2);
  pipmom_pimmom_wK0_woSid_won_mc->SetMinimum(1);
  pipmom_pimmom_wK0_woSid_won_mc->SetMaximum(pipmom_pimmom_wK0_woSid_won[0]->GetMaximum());
  pipmom_pimmom_wK0_woSid_won_mc->Draw("colz");
  
  
  //pimmom
  TH2D* pimmom_MMnmiss_wK0_woSid_won_mc = (TH2D*)pimmom_MMnmiss_wK0_woSid_won[1]->Clone("pimmom_MMnmiss_wK0_woSid_won_mc");

  TH1D* pimmom_wK0_woSid_won_mc = pimmom_MMnmiss_wK0_woSid_won_mc->ProjectionY("pimmom_wK0_woSid_won_mc");
  TH1D* pimmom_wK0_woSid_won[2];
  for(int i=0; i<2; i++) {
    pimmom_wK0_woSid_won[i] = (TH1D*)pimmom_MMnmiss_wK0_woSid_won[i]->ProjectionY(Form("pimmom_wK0_woSid_won_%s",name[i]));
    pimmom_wK0_woSid_won[i]->SetLineColor(colordef[i]);
  }
  TCanvas *cpimmom_wK0_woSid_won = new TCanvas("cpimmom_wK0_woSid_won","cpimmom_wK0_woSid_won");
  cpimmom_wK0_woSid_won->cd();
  pimmom_wK0_woSid_won[0]->Draw("HE");
  pimmom_wK0_woSid_won_mc->SetLineColor(6);
  pimmom_wK0_woSid_won_mc->Draw("HEsame");

  TCanvas *cpimmom_wK0_woSid_won_ratio = new TCanvas("cpimmom_wK0_woSid_won_ratio","cpimmom_wK0_woSid_won_ratio");
  cpimmom_wK0_woSid_won_ratio->cd();
  TH1D* pimmom_wK0_woSid_won_ratio = (TH1D*)pimmom_wK0_woSid_won[0]->Clone("pimmom_wK0_woSid_won_ratio");
  pimmom_wK0_woSid_won_ratio->Divide(pimmom_wK0_woSid_won_mc);
  pimmom_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  pimmom_wK0_woSid_won_ratio->SetTitle("Data/MC");
  pimmom_wK0_woSid_won_ratio->Draw("HEsame");

  TF1 *evalf_pimmom_wK0 = new TF1("evalf_pimmom_wK0","pol8",0.06,0.73);
  pimmom_wK0_woSid_won_ratio->Fit("evalf_pimmom_wK0","","",0.06,0.73);
  
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

  TCanvas *cnmom_wK0_woSid_won_ratio = new TCanvas("cnmom_wK0_woSid_won_ratio","cnmom_wK0_woSid_won_ratio");
  cnmom_wK0_woSid_won_ratio->cd();
  TH1D* nmom_wK0_woSid_won_ratio = (TH1D*)nmom_wK0_woSid_won[0]->Clone("nmom_wK0_woSid_won_ratio");
  nmom_wK0_woSid_won_ratio->Divide(nmom_wK0_woSid_won_mc);
  nmom_wK0_woSid_won_ratio->SetTitle("Data/MC");
  nmom_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  nmom_wK0_woSid_won_ratio->Draw("HEsame");

  TF1 *evalf_nmom_wK0 = new TF1("evalf_nmom_wK0","pol8",0.14,1);
  nmom_wK0_woSid_won_ratio->Fit("evalf_nmom_wK0","","",0.14,1);
  
  //weighting to momentum vector 
  TH2D* nmom_cosn_wK0_woSid_won_mc = (TH2D*)nmom_cosn_wK0_woSid_won[1]->Clone("nmom_cosn_wK0_woSid_won_mc");
  TCanvas *cnmom_cosn_wK0_woSid_won = new TCanvas("cnmom_cosn_wK0_woSid_won","cnmom_cosn_wK0_woSid_won",1200,800);
  cnmom_cosn_wK0_woSid_won->Divide(2,1);
  cnmom_cosn_wK0_woSid_won->cd(1);
  nmom_cosn_wK0_woSid_won[0]->SetTitle("#splitline{nmom_cosn_wK0_woSid_won}{Real data}");
  nmom_cosn_wK0_woSid_won[0]->Draw("colz");
  cnmom_cosn_wK0_woSid_won->cd(2);
  nmom_cosn_wK0_woSid_won_mc->SetTitle("#splitline{nmom_cosn_wK0_woSid_won}{MC}");
  nmom_cosn_wK0_woSid_won_mc->SetMinimum(1);
  nmom_cosn_wK0_woSid_won_mc->SetMaximum(nmom_cosn_wK0_woSid_won[0]->GetMaximum());
  nmom_cosn_wK0_woSid_won_mc->Draw("colz");

  TCanvas *ccosn_wK0_woSid_won = new TCanvas("ccosn_wK0_woSid_won","ccosn_wK0_woSid_won");
  ccosn_wK0_woSid_won->cd();
  TH1D* cosn_wK0_woSid_won[2];
  for(int i=0;i<2;i++) cosn_wK0_woSid_won[i] = (TH1D*)nmom_cosn_wK0_woSid_won[i]->ProjectionX(Form("cosn_wK0_woSid_won_%s",name[i]));
  TH1D* cosn_wK0_woSid_won_mc = (TH1D*)nmom_cosn_wK0_woSid_won_mc->ProjectionX("cosn_wK0_woSid_won_mc");
  cosn_wK0_woSid_won_mc->SetLineColor(6);
  cosn_wK0_woSid_won[0]->Draw("HE");
  cosn_wK0_woSid_won_mc->Draw("HEsame");
   

  TCanvas *ccosn_wK0_woSid_won_ratio = new TCanvas("ccosn_wK0_woSid_won_ratio","ccosn_wK0_woSid_won_ratio");
  ccosn_wK0_woSid_won_ratio->cd();
  TH1D* cosn_wK0_woSid_won_ratio = (TH1D*)cosn_wK0_woSid_won[0]->Clone("cosn_wK0_woSid_won_ratio");
  cosn_wK0_woSid_won_ratio->Divide(cosn_wK0_woSid_won_mc);
  cosn_wK0_woSid_won_ratio->SetTitle("Data/MC");
  cosn_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  cosn_wK0_woSid_won_ratio->Draw("HE");

  TF1 *evalf_cosn_wK0 = new TF1("evalf_cosn_wK0",func_cosn,-1.00,0.3,6);
  evalf_cosn_wK0->SetParameters(param_cosn);
  evalf_cosn_wK0->SetLineColor(3);
  cosn_wK0_woSid_won_ratio->Fit(evalf_cosn,"","",-1,0.2);
  
  //
  TH2D* nmom_cospip_wK0_woSid_won_mc = (TH2D*)nmom_cospip_wK0_woSid_won[1]->Clone("nmom_cospip_wK0_woSid_won_mc");
  TCanvas *cnmom_cospip_wK0_woSid_won = new TCanvas("cnmom_cospip_wK0_woSid_won","cnmom_cospip_wK0_woSid_won",1200,800);
  cnmom_cospip_wK0_woSid_won->Divide(2,1);
  cnmom_cospip_wK0_woSid_won->cd(1);
  nmom_cospip_wK0_woSid_won[0]->SetTitle("#splitline{nmom_cospip_wK0_woSid_won}{Real data}");
  nmom_cospip_wK0_woSid_won[0]->Draw("colz");
  cnmom_cospip_wK0_woSid_won->cd(2);
  nmom_cospip_wK0_woSid_won_mc->SetTitle("#splitline{nmom_cospip_wK0_woSid_won}{MC}");
  nmom_cospip_wK0_woSid_won_mc->SetMinimum(1);
  nmom_cospip_wK0_woSid_won_mc->SetMaximum(nmom_cospip_wK0_woSid_won[0]->GetMaximum());
  nmom_cospip_wK0_woSid_won_mc->Draw("colz");

  TCanvas *ccospip_wK0_woSid_won = new TCanvas("ccospip_wK0_woSid_won","ccospip_wK0_woSid_won");
  ccospip_wK0_woSid_won->cd();
  TH1D* cospip_wK0_woSid_won[2];
  for(int i=0;i<2;i++) cospip_wK0_woSid_won[i] = (TH1D*)nmom_cospip_wK0_woSid_won[i]->ProjectionX(Form("cospip_wK0_woSid_won_%s",name[i]));
  TH1D* cospip_wK0_woSid_won_mc = (TH1D*)nmom_cospip_wK0_woSid_won_mc->ProjectionX("cospip_wK0_woSid_won_mc");
  cospip_wK0_woSid_won_mc->SetLineColor(6);
  cospip_wK0_woSid_won[0]->Draw("HE");
  cospip_wK0_woSid_won_mc->Draw("HEsame");

  TCanvas *ccospip_wK0_woSid_won_ratio = new TCanvas("ccospip_wK0_woSid_won_ratio","ccospip_wK0_woSid_won_ratio");
  ccospip_wK0_woSid_won_ratio->cd();
  TH1D* cospip_wK0_woSid_won_ratio = (TH1D*)cospip_wK0_woSid_won[0]->Clone("cospip_wK0_woSid_won_ratio");
  cospip_wK0_woSid_won_ratio->Divide(cospip_wK0_woSid_won_mc);
  cospip_wK0_woSid_won_ratio->SetTitle("Data/MC");
  cospip_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  cospip_wK0_woSid_won_ratio->Draw("HE");
  
  TF1 *evalf_cospip_wK0 = new TF1("evalf_cospip_wK0","pol5",-1.0,0.6); 
  cospip_wK0_woSid_won_ratio->Fit("evalf_cospip_wK0","","",-1.0,0.6);

  //
  TH2D* nmom_cospim_wK0_woSid_won_mc = (TH2D*)nmom_cospim_wK0_woSid_won[1]->Clone("nmom_cospim_wK0_woSid_won_mc");
  TCanvas *cnmom_cospim_wK0_woSid_won = new TCanvas("cnmom_cospim_wK0_woSid_won","cnmom_cospim_wK0_woSid_won",1200,800);
  cnmom_cospim_wK0_woSid_won->Divide(2,1);
  cnmom_cospim_wK0_woSid_won->cd(1);
  nmom_cospim_wK0_woSid_won[0]->SetTitle("#splitline{nmom_cospim_wK0_woSid_won}{Real data}");
  nmom_cospim_wK0_woSid_won[0]->Draw("colz");
  cnmom_cospim_wK0_woSid_won->cd(2);
  nmom_cospim_wK0_woSid_won_mc->SetTitle("#splitline{nmom_cospim_wK0_woSid_won}{MC}");
  nmom_cospim_wK0_woSid_won_mc->SetMinimum(1);
  nmom_cospim_wK0_woSid_won_mc->SetMaximum(nmom_cospim_wK0_woSid_won[0]->GetMaximum());
  nmom_cospim_wK0_woSid_won_mc->Draw("colz");

  TCanvas *ccospim_wK0_woSid_won = new TCanvas("ccospim_wK0_woSid_won","ccospim_wK0_woSid_won");
  ccospim_wK0_woSid_won->cd();
  TH1D* cospim_wK0_woSid_won[2];
  for(int i=0;i<2;i++) cospim_wK0_woSid_won[i] = (TH1D*)nmom_cospim_wK0_woSid_won[i]->ProjectionX(Form("cospim_wK0_woSid_won_%s",name[i]));
  TH1D* cospim_wK0_woSid_won_mc = (TH1D*)nmom_cospim_wK0_woSid_won_mc->ProjectionX("cospim_wK0_woSid_won_mc");
  cospim_wK0_woSid_won_mc->SetLineColor(6);
  cospim_wK0_woSid_won[0]->Draw("HE");
  cospim_wK0_woSid_won_mc->Draw("HEsame");

  TCanvas *ccospim_wK0_woSid_won_ratio = new TCanvas("ccospim_wK0_woSid_won_ratio","ccospim_wK0_woSid_won_ratio");
  ccospim_wK0_woSid_won_ratio->cd();
  TH1D* cospim_wK0_woSid_won_ratio = (TH1D*)cospim_wK0_woSid_won[0]->Clone("cospim_wK0_woSid_won_ratio");
  cospim_wK0_woSid_won_ratio->Divide(cospim_wK0_woSid_won_mc);
  cospim_wK0_woSid_won_ratio->SetTitle("Data/MC");
  cospim_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  cospim_wK0_woSid_won_ratio->Draw("HE");
   
  TF1 *evalf_cospim_wK0 = new TF1("evalf_cospim_wK0","pol8",-0.92,0.55);
  cospim_wK0_woSid_won_ratio->Fit(evalf_cospim_wK0,"","",-0.92,0.55);

  //
  TH2D* nmom_phinpip_wK0_woSid_won_mc = (TH2D*)nmom_phinpip_wK0_woSid_won[1]->Clone("nmom_phinpip_wK0_woSid_won_mc");
  TCanvas *cnmom_phinpip_wK0_woSid_won = new TCanvas("cnmom_phinpip_wK0_woSid_won","cnmom_phinpip_wK0_woSid_won",1200,800);
  cnmom_phinpip_wK0_woSid_won->Divide(2,1);
  cnmom_phinpip_wK0_woSid_won->cd(1);
  nmom_phinpip_wK0_woSid_won[0]->SetTitle("#splitline{nmom_phinpip_wK0_woSid_won}{Real data}");
  nmom_phinpip_wK0_woSid_won[0]->Draw("colz");
  cnmom_phinpip_wK0_woSid_won->cd(2);
  nmom_phinpip_wK0_woSid_won_mc->SetTitle("#splitline{nmom_phinpip_wK0_woSid_won}{MC}");
  nmom_phinpip_wK0_woSid_won_mc->SetMinimum(1);
  nmom_phinpip_wK0_woSid_won_mc->SetMaximum(nmom_phinpip_wK0_woSid_won[0]->GetMaximum());
  nmom_phinpip_wK0_woSid_won_mc->Draw("colz");

  TCanvas *cphinpip_wK0_woSid_won = new TCanvas("cphinpip_wK0_woSid_won","cphinpip_wK0_woSid_won");
  cphinpip_wK0_woSid_won->cd();
  TH1D* phinpip_wK0_woSid_won[2];
  for(int i=0;i<2;i++) phinpip_wK0_woSid_won[i] = (TH1D*)nmom_phinpip_wK0_woSid_won[i]->ProjectionX(Form("phinpip_wK0_woSid_won_%s",name[i]));
  TH1D* phinpip_wK0_woSid_won_mc = (TH1D*)nmom_phinpip_wK0_woSid_won_mc->ProjectionX("phinpip_wK0_woSid_won_mc");
  phinpip_wK0_woSid_won_mc->SetLineColor(6);
  phinpip_wK0_woSid_won[0]->Draw("HE");
  phinpip_wK0_woSid_won_mc->Draw("HEsame");

  TCanvas *cphinpip_wK0_woSid_won_ratio = new TCanvas("cphinpip_wK0_woSid_won_ratio","cphinpip_wK0_woSid_won_ratio");
  cphinpip_wK0_woSid_won_ratio->cd();
  TH1D* phinpip_wK0_woSid_won_ratio = (TH1D*)phinpip_wK0_woSid_won[0]->Clone("phinpip_wK0_woSid_won_ratio");
  phinpip_wK0_woSid_won_ratio->Divide(phinpip_wK0_woSid_won_mc);
  phinpip_wK0_woSid_won_ratio->SetTitle("Data/MC");
  phinpip_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  phinpip_wK0_woSid_won_ratio->Draw("HE");
   
  TF1 *evalf_phinpip_wK0 = new TF1("evalf_phinpip_wK0",func_phinpip,-1.0*TMath::Pi(),TMath::Pi(),8);
  phinpip_wK0_woSid_won_ratio->Fit(evalf_phinpip_wK0);


  TH2D* nmom_phipip_wK0_woSid_won_mc = (TH2D*)nmom_phipip_wK0_woSid_won[1]->Clone("nmom_phipip_wK0_woSid_won_mc");
  TCanvas *cnmom_phipip_wK0_woSid_won = new TCanvas("cnmom_phipip_wK0_woSid_won","cnmom_phipip_wK0_woSid_won",1200,800);
  cnmom_phipip_wK0_woSid_won->Divide(2,1);
  cnmom_phipip_wK0_woSid_won->cd(1);
  nmom_phipip_wK0_woSid_won[0]->SetTitle("#splitline{nmom_phipip_wK0_woSid_won}{Real data}");
  nmom_phipip_wK0_woSid_won[0]->Draw("colz");
  cnmom_phipip_wK0_woSid_won->cd(2);
  nmom_phipip_wK0_woSid_won_mc->SetTitle("#splitline{nmom_phipip_wK0_woSid_won}{MC}");
  nmom_phipip_wK0_woSid_won_mc->SetMinimum(1);
  nmom_phipip_wK0_woSid_won_mc->SetMaximum(nmom_phipip_wK0_woSid_won[0]->GetMaximum());
  nmom_phipip_wK0_woSid_won_mc->Draw("colz");

  TCanvas *cphipip_wK0_woSid_won = new TCanvas("cphipip_wK0_woSid_won","cphipip_wK0_woSid_won");
  cphipip_wK0_woSid_won->cd();
  TH1D* phipip_wK0_woSid_won[2];
  for(int i=0;i<2;i++) phipip_wK0_woSid_won[i] = (TH1D*)nmom_phipip_wK0_woSid_won[i]->ProjectionX(Form("phipip_wK0_woSid_won_%s",name[i]));
  TH1D* phipip_wK0_woSid_won_mc = (TH1D*)nmom_phipip_wK0_woSid_won_mc->ProjectionX("phipip_wK0_woSid_won_mc");
  phipip_wK0_woSid_won_mc->SetLineColor(6);
  phipip_wK0_woSid_won[0]->Draw("HE");
  phipip_wK0_woSid_won_mc->Draw("HEsame");

  TCanvas *cphipip_wK0_woSid_won_ratio = new TCanvas("cphipip_wK0_woSid_won_ratio","cphipip_wK0_woSid_won_ratio");
  cphipip_wK0_woSid_won_ratio->cd();
  TH1D* phipip_wK0_woSid_won_ratio = (TH1D*)phipip_wK0_woSid_won[0]->Clone("phipip_wK0_woSid_won_ratio");
  phipip_wK0_woSid_won_ratio->Divide(phipip_wK0_woSid_won_mc);
  phipip_wK0_woSid_won_ratio->SetTitle("Data/MC");
  phipip_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  phipip_wK0_woSid_won_ratio->Draw("HE");

  //
  TH2D* nmom_phinpim_wK0_woSid_won_mc = (TH2D*)nmom_phinpim_wK0_woSid_won[1]->Clone("nmom_phinpim_wK0_woSid_won_mc");
  TCanvas *cnmom_phinpim_wK0_woSid_won = new TCanvas("cnmom_phinpim_wK0_woSid_won","cnmom_phinpim_wK0_woSid_won",1200,800);
  cnmom_phinpim_wK0_woSid_won->Divide(2,1);
  cnmom_phinpim_wK0_woSid_won->cd(1);
  nmom_phinpim_wK0_woSid_won[0]->SetTitle("#splitline{nmom_phinpim_wK0_woSid_won}{Real data}");
  nmom_phinpim_wK0_woSid_won[0]->Draw("colz");
  cnmom_phinpim_wK0_woSid_won->cd(2);
  nmom_phinpim_wK0_woSid_won_mc->SetTitle("#splitline{nmom_phinpim_wK0_woSid_won}{MC}");
  nmom_phinpim_wK0_woSid_won_mc->SetMinimum(1);
  nmom_phinpim_wK0_woSid_won_mc->SetMaximum(nmom_phinpim_wK0_woSid_won[0]->GetMaximum());
  nmom_phinpim_wK0_woSid_won_mc->Draw("colz");

  TCanvas *cphinpim_wK0_woSid_won = new TCanvas("cphinpim_wK0_woSid_won","cphinpim_wK0_woSid_won");
  cphinpim_wK0_woSid_won->cd();
  TH1D* phinpim_wK0_woSid_won[1];
  for(int i=0;i<1;i++) phinpim_wK0_woSid_won[i] = (TH1D*)nmom_phinpim_wK0_woSid_won[i]->ProjectionX(Form("phinpim_wK0_woSid_won_%s",name[i]));
  TH1D* phinpim_wK0_woSid_won_mc = (TH1D*)nmom_phinpim_wK0_woSid_won_mc->ProjectionX("phinpim_wK0_woSid_won_mc");
  phinpim_wK0_woSid_won_mc->SetLineColor(6);
  phinpim_wK0_woSid_won[0]->Draw("HE");
  phinpim_wK0_woSid_won_mc->Draw("HEsame");

  TCanvas *cphinpim_wK0_woSid_won_ratio = new TCanvas("cphinpim_wK0_woSid_won_ratio","cphinpim_wK0_woSid_won_ratio");
  cphinpim_wK0_woSid_won_ratio->cd();
  TH1D* phinpim_wK0_woSid_won_ratio = (TH1D*)phinpim_wK0_woSid_won[0]->Clone("phinpim_wK0_woSid_won_ratio");
  phinpim_wK0_woSid_won_ratio->Divide(phinpim_wK0_woSid_won_mc);
  phinpim_wK0_woSid_won_ratio->SetTitle("Data/MC");
  phinpim_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  phinpim_wK0_woSid_won_ratio->Draw("HE");

  TF1 *evalf_phinpim_wK0 = new TF1("evalf_phinpim_wK0",func_phinpim_wK0,-1.0*TMath::Pi(),TMath::Pi(),10);
  phinpim_wK0_woSid_won_ratio->Fit(evalf_phinpim_wK0);

  TH2D* nmom_phipim_wK0_woSid_won_mc = (TH2D*)nmom_phipim_wK0_woSid_won[1]->Clone("nmom_phipim_wK0_woSid_won_mc");
  TCanvas *cnmom_phipim_wK0_woSid_won = new TCanvas("cnmom_phipim_wK0_woSid_won","cnmom_phipim_wK0_woSid_won",1200,800);
  cnmom_phipim_wK0_woSid_won->Divide(2,1);
  cnmom_phipim_wK0_woSid_won->cd(1);
  nmom_phipim_wK0_woSid_won[0]->SetTitle("#splitline{nmom_phipim_wK0_woSid_won}{Real data}");
  nmom_phipim_wK0_woSid_won[0]->Draw("colz");
  cnmom_phipim_wK0_woSid_won->cd(2);
  nmom_phipim_wK0_woSid_won_mc->SetTitle("#splitline{nmom_phipim_wK0_woSid_won}{MC}");
  nmom_phipim_wK0_woSid_won_mc->SetMinimum(1);
  nmom_phipim_wK0_woSid_won_mc->SetMaximum(nmom_phipim_wK0_woSid_won[0]->GetMaximum());
  nmom_phipim_wK0_woSid_won_mc->Draw("colz");

  TCanvas *cphipim_wK0_woSid_won = new TCanvas("cphipim_wK0_woSid_won","cphipim_wK0_woSid_won");
  cphipim_wK0_woSid_won->cd();
  TH1D* phipim_wK0_woSid_won[1];
  for(int i=0;i<1;i++) phipim_wK0_woSid_won[i] = (TH1D*)nmom_phipim_wK0_woSid_won[i]->ProjectionX(Form("phipim_wK0_woSid_won_%s",name[i]));
  TH1D* phipim_wK0_woSid_won_mc = (TH1D*)nmom_phipim_wK0_woSid_won_mc->ProjectionX("phipim_wK0_woSid_won_mc");
  phipim_wK0_woSid_won_mc->SetLineColor(6);
  phipim_wK0_woSid_won[0]->Draw("HE");
  phipim_wK0_woSid_won_mc->Draw("HEsame");

  TCanvas *cphipim_wK0_woSid_won_ratio = new TCanvas("cphipim_wK0_woSid_won_ratio","cphipim_wK0_woSid_won_ratio");
  cphipim_wK0_woSid_won_ratio->cd();
  TH1D* phipim_wK0_woSid_won_ratio = (TH1D*)phipim_wK0_woSid_won[0]->Clone("phipim_wK0_woSid_won_ratio");
  phipim_wK0_woSid_won_ratio->Divide(phipim_wK0_woSid_won_mc);
  phipim_wK0_woSid_won_ratio->SetTitle("Data/MC");
  phipim_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  phipim_wK0_woSid_won_ratio->Draw("HE");
  
  
  TH2D* nmom_phin_wK0_woSid_won_mc = (TH2D*)nmom_phin_wK0_woSid_won[1]->Clone("nmom_phin_wK0_woSid_won_mc");
  TCanvas *cnmom_phin_wK0_woSid_won = new TCanvas("cnmom_phin_wK0_woSid_won","cnmom_phin_wK0_woSid_won",1200,800);
  cnmom_phin_wK0_woSid_won->Divide(2,1);
  cnmom_phin_wK0_woSid_won->cd(1);
  nmom_phin_wK0_woSid_won[0]->SetTitle("#splitline{nmom_phin_wK0_woSid_won}{Real data}");
  nmom_phin_wK0_woSid_won[0]->Draw("colz");
  cnmom_phin_wK0_woSid_won->cd(2);
  nmom_phin_wK0_woSid_won_mc->SetTitle("#splitline{nmom_phin_wK0_woSid_won}{MC}");
  nmom_phin_wK0_woSid_won_mc->SetMinimum(1);
  nmom_phin_wK0_woSid_won_mc->SetMaximum(nmom_phin_wK0_woSid_won[0]->GetMaximum());
  nmom_phin_wK0_woSid_won_mc->Draw("colz");

  TCanvas *cphin_wK0_woSid_won = new TCanvas("cphin_wK0_woSid_won","cphin_wK0_woSid_won");
  cphin_wK0_woSid_won->cd();
  TH1D* phin_wK0_woSid_won[1];
  for(int i=0;i<1;i++) phin_wK0_woSid_won[i] = (TH1D*)nmom_phin_wK0_woSid_won[i]->ProjectionX(Form("phin_wK0_woSid_won_%s",name[i]));
  TH1D* phin_wK0_woSid_won_mc = (TH1D*)nmom_phin_wK0_woSid_won_mc->ProjectionX("phin_wK0_woSid_won_mc");
  phin_wK0_woSid_won_mc->SetLineColor(6);
  phin_wK0_woSid_won[0]->Draw("HE");
  phin_wK0_woSid_won_mc->Draw("HEsame");

  TCanvas *cphin_wK0_woSid_won_ratio = new TCanvas("cphin_wK0_woSid_won_ratio","cphin_wK0_woSid_won_ratio");
  cphin_wK0_woSid_won_ratio->cd();
  TH1D* phin_wK0_woSid_won_ratio = (TH1D*)phin_wK0_woSid_won[0]->Clone("phin_wK0_woSid_won_ratio");
  phin_wK0_woSid_won_ratio->Divide(phin_wK0_woSid_won_mc);
  phin_wK0_woSid_won_ratio->SetTitle("Data/MC");
  phin_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,3);
  phin_wK0_woSid_won_ratio->Draw("HE");

  //Signal check
  TH2D* q_IMnpipi_wSid_n_Sp_mc = (TH2D*)q_IMnpipi_wSid_n_Sp[1]->Clone();

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
  TH1D* IMnpipi_wSid_n_Sp_0[1];
  for(int i=0; i<1; i++)IMnpipi_wSid_n_Sp_0[i] = q_IMnpipi_wSid_n_Sp[i]->ProjectionX(Form("IMnpipi_wSid_n_Sp_0_%s",name[i]),0,q350bin);
  IMnpipi_wSid_n_Sp_0[0]->Draw("HE");
  TH1D* IMnpipi_wSid_n_Sp_mc_0 = q_IMnpipi_wSid_n_Sp_mc->ProjectionX("IMnpipi_wSid_n_Sp_mc_0",0,q350bin-1);
  IMnpipi_wSid_n_Sp_mc_0->SetLineColor(6);
  IMnpipi_wSid_n_Sp_mc_0->Draw("HEsame");

  TCanvas *cIMnpipi_wSid_n_Sp_350 = new TCanvas("cIMnpipi_wSid_n_Sp_350","cIMnpipi_wSid_n_Sp_350");
  TH1D* IMnpipi_wSid_n_Sp_350[2];
  for(int i=0; i<1; i++)IMnpipi_wSid_n_Sp_350[i] = q_IMnpipi_wSid_n_Sp[i]->ProjectionX(Form("IMnpipi_wSid_n_Sp_350_%s",name[i]),q350bin,100);
  IMnpipi_wSid_n_Sp_350[0]->Draw("HE");
  TH1D* IMnpipi_wSid_n_Sp_mc_350 = q_IMnpipi_wSid_n_Sp_mc->ProjectionX("IMnpipi_wSid_n_Sp_mc_350",q350bin,100);
  IMnpipi_wSid_n_Sp_mc_350->SetLineColor(6);
  IMnpipi_wSid_n_Sp_mc_350->Draw("HEsame");


  TH2D* q_IMnpipi_wSid_n_Sm_mc = (TH2D*)q_IMnpipi_wSid_n_Sm[1]->Clone();

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
  TH1D* IMnpipi_wSid_n_Sm_0[2];
  for(int i=0; i<2; i++)IMnpipi_wSid_n_Sm_0[i] = q_IMnpipi_wSid_n_Sm[i]->ProjectionX(Form("IMnpipi_wSid_n_Sm_0_%s",name[i]),0,q350bin);
  IMnpipi_wSid_n_Sm_0[0]->Draw("HE");
  TH1D* IMnpipi_wSid_n_Sm_mc_0 = q_IMnpipi_wSid_n_Sm_mc->ProjectionX("IMnpipi_wSid_n_Sm_mc_0",0,q350bin-1);
  IMnpipi_wSid_n_Sm_mc_0->SetLineColor(6);
  IMnpipi_wSid_n_Sm_mc_0->Draw("HEsame");

  TCanvas *cIMnpipi_wSid_n_Sm_350 = new TCanvas("cIMnpipi_wSid_n_Sm_350","cIMnpipi_wSid_n_Sm_350");
  TH1D* IMnpipi_wSid_n_Sm_350[2];
  for(int i=0; i<2; i++)IMnpipi_wSid_n_Sm_350[i] = q_IMnpipi_wSid_n_Sm[i]->ProjectionX(Form("IMnpipi_wSid_n_Sm_350_%s",name[i]),q350bin,100);
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
  
  TCanvas *cfunc_MMnmiss = new TCanvas("cfunc_MMnmiss","cfunc_MMnmiss");
  cfunc_MMnmiss->Divide(2,2);
  cfunc_MMnmiss->cd(1);
  evalf_MMnmiss->Draw();
  cfunc_MMnmiss->cd(2);
  evalf_MMnmiss_corr->Draw();
  cfunc_MMnmiss->cd(3);
  TF1 *evalf_MMnmiss_mul = new TF1("evalf_MMnmiss_mul",func_MMnmiss_mul,0,1.5,19);
  double par_MMnmiss_mul[19];
  evalf_MMnmiss->GetParameters(&par_MMnmiss_mul[0]);
  evalf_MMnmiss_corr->GetParameters(&par_MMnmiss_mul[12]);
  evalf_MMnmiss_mul->SetParameters(par_MMnmiss_mul);
  evalf_MMnmiss_mul->Draw();

  TCanvas *c = NULL;
  TSeqCollection *SCol = gROOT->GetListOfCanvases();
  int size = SCol->GetSize();
  TIter next(SCol);
  TString pdfname="comp_fakedata_out.pdf";
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

  evalf_MMnmiss->Write();
//  evalf_MMnmiss_corr->Write();
//  evalf_MMnmiss_corr2->Write();
  evalf_MMnmiss_wK0->Write();
 // evalf_MMom->Write();
//  evalf_nmom->Write();
//  evalf_nmom_wK0->Write();
//  evalf_pipmom->Write();
//  evalf_pipmom_wK0->Write();
//  evalf_pimmom->Write();
//  evalf_IMpippim->Write();
//  evalf_Mompippim->Write();
//  evalf_IMnpipi->Write();
//  evalf_IMnpip->Write();
//  evalf_IMnpim->Write();
//  evalf_MMom_wK0->Write();
//  evalf_pimmom_wK0->Write();
//  evalf_cosn->Write();
//  evalf_cosn_corr->Write();
//  evalf_cosn_wK0->Write();
//  evalf_cospip->Write();
//  evalf_cospip_wK0->Write();
//  evalf_cospim->Write();
//  evalf_cospim_wK0->Write();
//  evalf_phinpip->Write();
//  evalf_phinpip_wK0->Write();
//  evalf_phinpim->Write();
//  evalf_phinpim_wK0->Write();
//  evalf_IMnpipi_wK0->Write();
//  evalf_IMnpip_wK0->Write();
//  evalf_IMnpim_wK0->Write();
  fout->Print();
  fout->cd();
  fout->Close();
}

