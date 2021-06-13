const bool showBG = true;
#include "anacuts.h"

#include <iostream>
#include <vector>
#include <string>
#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TMultiGraph.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TLatex.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>

#include "../src/GlobalVariables.h"

bool RemoveNegative = true;
bool RebinMode = true;
bool Sidefar=false;
bool FitNoWeight=true;


void plot_AfterDecompos()
{
  TFile *fr[4]={NULL};
  
  fr[0] = TFile::Open("evanaIMpisigma_npippim_v202_out_iso_sub.root");
  fr[1] = TFile::Open("evanaIMpisigma_npippim_v202_out_iso_qlo_sub.root");
  fr[2] = TFile::Open("evanaIMpisigma_npippim_v202_out_iso_qhi_sub.root");
  fr[3] = TFile::Open("evanaIMpisigma_npippim_v202_out_iso_theta15_sub.root");
  fr[0]->Print();
  fr[1]->Print();
  fr[2]->Print();
  fr[3]->Print();
  TFile *fdeco[2]={NULL};
  fdeco[0] = TFile::Open("deco_qlo.root","READ");
  fdeco[1] = TFile::Open("deco_qhi.root","READ");
  fdeco[0]->Print();
  fdeco[1]->Print();

  //gROOT->SetBatch(1);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0.);  

  const unsigned int nbinIMnpipi = 80;
  const int nqcut=4;
  const int qstart=1;

  //for the overlap of S+ & S- & K0 counting 
  TH2F* q_IMnpipi_wK0_wSid_n_SpSm[nqcut];
  TH1D* IMnpipi_wK0_wSid_n_SpSm[nqcut];
  TCanvas *cq_IMnpipi_wK0_wSid_n_SpSm[nqcut];
  TCanvas *cIMnpipi_Sp[nqcut];
  TCanvas *cIMnpipi_Sm[nqcut];
  TCanvas *cIMnpipi_K0[nqcut];
  double OverlapCount[nbinIMnpipi][nqcut]; 
  const unsigned int nwbin = 3;
  TH2D* IMnpim_IMnpip_dE_wK0_wSid_n_Sp_bin[nbinIMnpipi][nqcut];
  TH2D* IMnpim_IMnpip_dE_wK0_wSid_n_Sm_bin[nbinIMnpipi][nqcut];
  TH2D* IMnpim_IMnpip_dE_wK0_wSid_n_Sp_wbin[nwbin][nqcut];
  TH2D* IMnpim_IMnpip_dE_wK0_wSid_n_Sm_wbin[nwbin][nqcut];
  TH2D* q_IMnpipi_wSid_n_Sp[nqcut];
  TH1D* IMnpipi_wSid_n_Sp[nqcut];
  TH2D* q_IMnpipi_wK0_wSid_n_Sp[nqcut];
  TH1D* IMnpipi_wK0_wSid_n_Sp[nqcut];
  TH2D* q_IMnpipi_wSid_n_Sm[nqcut];
  TH1D* IMnpipi_wSid_n_Sm[nqcut];
  TH2D* q_IMnpipi_wK0_wSid_n_Sm[nqcut];
  TH1D* IMnpipi_wK0_wSid_n_Sm[nqcut];
  TH2D* q_IMnpipi_wSid_n_SpSm[nqcut];
  TH1D* IMnpipi_wSid_n_SpSm[nqcut];
  TH2D* q_IMnpipi_wK0_n[nqcut];
  TH1D* IMnpipi_wK0_n[nqcut];


  const char cqcut[][10]= {"all","qlo","qhi","theta"};
  std::cout << __LINE__ << std::endl;
  for(int iq=qstart;iq<nqcut;iq++){
    //get
    q_IMnpipi_wK0_wSid_n_SpSm[iq] = (TH2F*)fr[iq]->Get("q_IMnpipi_wK0_wSid_n_SpSm");
    q_IMnpipi_wSid_n_Sp[iq] = (TH2D*)fr[iq]->Get("q_IMnpipi_wSid_n_Sp");
    q_IMnpipi_wK0_wSid_n_Sp[iq] = (TH2D*)fr[iq]->Get("q_IMnpipi_wK0_wSid_n_Sp");
    q_IMnpipi_wSid_n_Sm[iq] = (TH2D*)fr[iq]->Get("q_IMnpipi_wSid_n_Sm");
    q_IMnpipi_wK0_wSid_n_Sm[iq] = (TH2D*)fr[iq]->Get("q_IMnpipi_wK0_wSid_n_Sm");
    q_IMnpipi_wSid_n_SpSm[iq] = (TH2D*)fr[iq]->Get("q_IMnpipi_wSid_n_SpSm");
    q_IMnpipi_wK0_n[iq] = (TH2D*)fr[iq]->Get("q_IMnpipi_wK0_n");
    
    //Draw
    cIMnpipi_Sp[iq] = new TCanvas(Form("IMnpipi_Sp_%s",cqcut[iq]),Form("IMnpipi_Sp_%s",cqcut[iq]),800,800);
    IMnpipi_wSid_n_Sp[iq] = (TH1D*)q_IMnpipi_wSid_n_Sp[iq]->ProjectionX(Form("IMnpipi_wSid_n_Sp_%d",iq));
    IMnpipi_wSid_n_Sp[iq]->Draw("HE");
    IMnpipi_wK0_wSid_n_Sp[iq] = (TH1D*)q_IMnpipi_wK0_wSid_n_Sp[iq]->ProjectionX(Form("IMnpipi_wK0_wSid_n_Sp_%d",iq));
    IMnpipi_wK0_wSid_n_Sp[iq]->SetLineColor(2);
    IMnpipi_wK0_wSid_n_Sp[iq]->Draw("HEsame");
    IMnpipi_wSid_n_SpSm[iq] = (TH1D*)q_IMnpipi_wSid_n_SpSm[iq]->ProjectionX(Form("IMnpipi_wSid_n_SpSm_%d",iq));
    IMnpipi_wSid_n_SpSm[iq]->SetLineColor(3);
    IMnpipi_wSid_n_SpSm[iq]->Draw("HEsame");

    cIMnpipi_Sm[iq] = new TCanvas(Form("IMnpipi_Sm_%s",cqcut[iq]),Form("IMnpipi_Sm_%s",cqcut[iq]),800,800);
    IMnpipi_wSid_n_Sm[iq] = (TH1D*)q_IMnpipi_wSid_n_Sm[iq]->ProjectionX(Form("IMnpipi_wSid_n_Sm_%d",iq));
    IMnpipi_wSid_n_Sm[iq]->Draw("HE");
    IMnpipi_wK0_wSid_n_Sm[iq] = (TH1D*)q_IMnpipi_wK0_wSid_n_Sm[iq]->ProjectionX(Form("IMnpipi_wK0_wSid_n_Sm_%d",iq));
    IMnpipi_wK0_wSid_n_Sm[iq]->SetLineColor(2);
    IMnpipi_wK0_wSid_n_Sm[iq]->Draw("HEsame");
    IMnpipi_wSid_n_SpSm[iq]->Draw("HEsame");
    
    cIMnpipi_K0[iq] = new TCanvas(Form("IMnpipi_K0_%s",cqcut[iq]),Form("IMnpipi_K0_%s",cqcut[iq]),800,800);
    IMnpipi_wK0_n[iq] = (TH1D*)q_IMnpipi_wK0_n[iq]->ProjectionX(Form("IMnpipi_wK0_n_%d",iq));
    IMnpipi_wK0_n[iq]->Draw("HE");
    IMnpipi_wK0_wSid_n_Sp[iq]->Draw("HEsame");
    IMnpipi_wK0_wSid_n_Sm[iq]->Draw("HEsame");
    
    cq_IMnpipi_wK0_wSid_n_SpSm[iq] = new TCanvas(Form("q_IMnpipi_wK0_wSid_n_SpSm_%s",cqcut[iq]),Form("q_IMnpipi_wK0_wSid_n_SpSm_%s",cqcut[iq]));
    IMnpipi_wK0_wSid_n_SpSm[iq] = (TH1D*)q_IMnpipi_wK0_wSid_n_SpSm[iq]->ProjectionX(Form("IMnpipi_wK0_wSid_n_SpSm_%d",iq));
    IMnpipi_wK0_wSid_n_SpSm[iq]->Draw("HISTE");
    for(unsigned int ibin=0;ibin<nbinIMnpipi;ibin++){
      IMnpim_IMnpip_dE_wK0_wSid_n_Sp_bin[ibin][iq] = (TH2D*)fr[iq]->Get(Form("IMnpim_IMnpip_dE_wK0_wSid_n_Sp_bin%d",ibin));
      IMnpim_IMnpip_dE_wK0_wSid_n_Sm_bin[ibin][iq] = (TH2D*)fr[iq]->Get(Form("IMnpim_IMnpip_dE_wK0_wSid_n_Sm_bin%d",ibin));
      OverlapCount[ibin][iq] = IMnpipi_wK0_wSid_n_SpSm[iq]->GetBinContent(ibin+1);
      if(OverlapCount[ibin][iq]<0.0) OverlapCount[ibin][iq]=0.0;
    }
    for(int iwbin=0;iwbin<nwbin;iwbin++){
      IMnpim_IMnpip_dE_wK0_wSid_n_Sp_wbin[iwbin][iq] = (TH2D*)fr[iq]->Get(Form("IMnpim_IMnpip_dE_wK0_wSid_n_Sp_wbin%d",iwbin));
      IMnpim_IMnpip_dE_wK0_wSid_n_Sm_wbin[iwbin][iq] = (TH2D*)fr[iq]->Get(Form("IMnpim_IMnpip_dE_wK0_wSid_n_Sm_wbin%d",iwbin));
    }
  }

  //only qlo(=0) and qhi(1) decomposition results;
  TH2D* IMnpim_IMnpip_K0inter[2];
  TGraphErrors *gr_SpONnpim_fin_pol1[2];
  TGraphErrors *gr_SmONnpip_fin_pol1[2];
  TGraphErrors *gr_SpONnpim_fin_3rd[2];
  TGraphErrors *gr_SmONnpip_fin_3rd[2];

  //only qlo(=0) and qhi(1) decomposition results;
  for(int iq=0;iq<2;iq++){
    IMnpim_IMnpip_K0inter[iq] = (TH2D*)fdeco[iq]->Get("h2K0inter_3fine");
    gr_SpONnpim_fin_pol1[iq] = (TGraphErrors*)fdeco[iq]->Get("gr_SpONnpim_fin_pol1");
    gr_SmONnpip_fin_pol1[iq] = (TGraphErrors*)fdeco[iq]->Get("gr_SmONnpip_fin_pol1");
    gr_SpONnpim_fin_3rd[iq] = (TGraphErrors*)fdeco[iq]->Get("gr_SpONnpim_fin_3rd");
    gr_SmONnpip_fin_3rd[iq] = (TGraphErrors*)fdeco[iq]->Get("gr_SmONnpip_fin_3rd");
  }
  
  //have data only in Sigma+/- region
  TCanvas *cK0inter = new TCanvas("cK0inter","cK0inter",1600,800);
  cK0inter->Divide(2,1);
  cK0inter->cd(1);
  IMnpim_IMnpip_K0inter[0]->Draw("colz");
  cK0inter->cd(2);
  IMnpim_IMnpip_K0inter[1]->Draw("colz");
  
  


  //only qlo(=0) and qhi(1) decomposition results;
  for(int iq=0;iq<2;iq++){
    //no overlap in iwbin=0;
    for(int iwbin=1;iwbin<nwbin;iwbin++){
      for(int ix=0;ix<IMnpim_IMnpip_dE_wK0_wSid_n_Sp_wbin[iwbin][iq+1]->GetNbinsX();ix++){
        for(int iy=0;iy<IMnpim_IMnpip_dE_wK0_wSid_n_Sp_wbin[iwbin][iq+1]->GetNbinsY();iy++){
          double cont = IMnpim_IMnpip_dE_wK0_wSid_n_Sp_wbin[iwbin][iq+1]->GetBinContent(ix,iy);
        }
      }
    }
  }





}
/*
{
  std::cout << __LINE__ << std::endl;
  int ngroup = 3;
  TH1D* IMnpipi_overlapdeco_K0[ngroup][nqcut];
  TH1D* IMnpipi_overlapdeco_Sp[ngroup][nqcut];
  TH1D* IMnpipi_overlapdeco_Sm[ngroup][nqcut];

  for(int ig=0;ig<ngroup;ig++){
    for(int iq=qstart;iq<nqcut;iq++){
      if(ig!=2)IMnpipi_overlapdeco_K0[ig][iq] = new TH1D(Form("IMnpipi_overlapdeco_K0_g%d_%d",ig,iq),Form("IMnpipi_overlapdeco_K0_g%d_%d",ig,iq),100,1,2);
      if(ig!=1)IMnpipi_overlapdeco_Sp[ig][iq] = new TH1D(Form("IMnpipi_overlapdeco_Sp_g%d_%d",ig,iq),Form("IMnpipi_overlapdeco_Sp_g%d_%d",ig,iq),100,1,2);
      if(ig!=0)IMnpipi_overlapdeco_Sm[ig][iq] = new TH1D(Form("IMnpipi_overlapdeco_Sm_g%d_%d",ig,iq),Form("IMnpipi_overlapdeco_Sm_g%d_%d",ig,iq),100,1,2);
    }
  }
   
  std::cout << __LINE__ << std::endl;
  TCanvas *csum_Sp[ngroup][nqcut];
  TCanvas *csum_Sm[ngroup][nqcut];
  TCanvas *csum_K0[ngroup][nqcut];
  for(int ig=0;ig<ngroup;ig++){
    for(int iq=qstart;iq<nqcut;iq++){
      if(ig!=2){
        csum_K0[ig][iq] = new TCanvas(Form("csum_K0_g%d_q%d",ig,iq),Form("csum_K0_g%d_q%d",ig,iq));
        IMnpipi_overlapdeco_K0[ig][iq]->Draw("H");
      }
      if(ig!=1){  
        csum_Sp[ig][iq] = new TCanvas(Form("csum_Sp_g%d_q%d",ig,iq),Form("csum_Sp_g%d_q%d",ig,iq));
        IMnpipi_overlapdeco_Sp[ig][iq]->Draw("H");
      }
      if(ig!=0){
        csum_Sm[ig][iq] = new TCanvas(Form("csum_Sm_g%d_q%d",ig,iq),Form("csum_Sm_g%d_q%d",ig,iq));
        IMnpipi_overlapdeco_Sm[ig][iq]->Draw("H");
      }
    }
  }

  std::cout << __LINE__ << std::endl;
  TCanvas *cratio_SpSm[nqcut];//Sp/Sm  ->group 2
  TCanvas *cratio_SpK0[nqcut];//Sp/K0  ->group 0
  TCanvas *cratio_SmK0[nqcut];//Sm/K0  ->group 1
  TH1D* IMnpipi_SpoverSm_ratio[nqcut];
  TH1D* IMnpipi_SpoverK0_ratio[nqcut];
  TH1D* IMnpipi_SmoverK0_ratio[nqcut];

  for(int iq=qstart;iq<nqcut;iq++){
    IMnpipi_SpoverSm_ratio[iq] = new TH1D(Form("IMnpipi_SpoverSm_ratio_%d",iq),Form("IMnpipi_SpoverSm_ratio_%s",cqcut[iq]),100,1,2);
    IMnpipi_SpoverK0_ratio[iq] = new TH1D(Form("IMnpipi_SpoverK0_ratio_%d",iq),Form("IMnpipi_SpoverK0_ratio_%s",cqcut[iq]),100,1,2);
    IMnpipi_SmoverK0_ratio[iq] = new TH1D(Form("IMnpipi_SmoverK0_ratio_%d",iq),Form("IMnpipi_SmoverK0_ratio_%s",cqcut[iq]),100,1,2);
  
    IMnpipi_SpoverSm_ratio[iq] = (TH1D*)IMnpipi_overlapdeco_Sp[2][iq]->Clone(Form("IMnpipi_SpoverSm_ratio_%d",iq));
    IMnpipi_SpoverSm_ratio[iq]->Divide(IMnpipi_overlapdeco_Sm[2][iq]);
    cratio_SpSm[iq] = new TCanvas(Form("cratrio_SpSm_%s",cqcut[iq]));
    IMnpipi_SpoverSm_ratio[iq]->SetTitle(Form("IMnpipi_Sp/Sm_ratio_%s",cqcut[iq]));
    IMnpipi_SpoverSm_ratio[iq]->Draw("HIST");
    
    IMnpipi_SpoverK0_ratio[iq] = (TH1D*)IMnpipi_overlapdeco_Sp[0][iq]->Clone(Form("IMnpipi_SpoverK0_ratio_%d",iq));
    IMnpipi_SpoverK0_ratio[iq]->Divide(IMnpipi_overlapdeco_K0[0][iq]);
    cratio_SpK0[iq] = new TCanvas(Form("cratrio_SpK0_%s",cqcut[iq]));
    IMnpipi_SpoverK0_ratio[iq]->SetTitle(Form("IMnpipi_Sp/K0_ratio_%s",cqcut[iq]));
    IMnpipi_SpoverK0_ratio[iq]->Draw("HIST");

    IMnpipi_SmoverK0_ratio[iq] = (TH1D*)IMnpipi_overlapdeco_Sm[1][iq]->Clone(Form("IMnpipi_SmoverK0_ratio_%d",iq));
    IMnpipi_SmoverK0_ratio[iq]->Divide(IMnpipi_overlapdeco_K0[1][iq]);
    cratio_SmK0[iq] = new TCanvas(Form("cratrio_SmK0_%s",cqcut[iq]));
    IMnpipi_SmoverK0_ratio[iq]->SetTitle(Form("IMnpipi_Sm/K0_ratio_%s",cqcut[iq]));
    IMnpipi_SmoverK0_ratio[iq]->Draw("HIST");
  }
  
  //tuned bins 2d histograms
  TH2F* IMnpim_IMnpip_wSid_n_Sp_wbin[nwbin][nqcut];
  TH2F* IMnpim_IMnpip_wSid_n_Sm_wbin[nwbin][nqcut];
  TH2F* IMpippim_IMnpip_wSid_n_Sp_wbin[nwbin][nqcut];
  TH2F* IMpippim_IMnpip_wK0_n_wbin[nwbin][nqcut];
  TH2F* IMpippim_IMnpim_wSid_n_Sm_wbin[nwbin][nqcut];
  TH2F* IMpippim_IMnpim_wK0_n_wbin[nwbin][nqcut];
  TH2F* IMnpim_IMnpip_wK0orwSid_n_wbin[nwbin][nqcut];
  TH2F* IMpippim_IMnpip_wK0orwSid_n_wbin[nwbin][nqcut];
  TH2F* IMpippim_IMnpim_wK0orwSid_n_wbin[nwbin][nqcut];
  for(int iq=0;iq<nqcut;iq++){
    for(int iwbin=0;iwbin<nwbin;iwbin++){
      IMnpim_IMnpip_wSid_n_Sp_wbin[iwbin][iq]= (TH2F*)fr[iq]->Get(Form("IMnpim_IMnpip_wSid_n_Sp_wbin%d",iwbin));
      IMnpim_IMnpip_wSid_n_Sm_wbin[iwbin][iq]= (TH2F*)fr[iq]->Get(Form("IMnpim_IMnpip_wSid_n_Sm_wbin%d",iwbin));
      IMpippim_IMnpip_wSid_n_Sp_wbin[iwbin][iq]= (TH2F*)fr[iq]->Get(Form("IMpippim_IMnpip_wSid_n_Sp_wbin%d",iwbin));
      IMpippim_IMnpip_wK0_n_wbin[iwbin][iq]= (TH2F*)fr[iq]->Get(Form("IMpippim_IMnpip_wK0_n_wbin%d",iwbin));
      IMpippim_IMnpim_wSid_n_Sm_wbin[iwbin][iq]= (TH2F*)fr[iq]->Get(Form("IMpippim_IMnpim_wSid_n_Sm_wbin%d",iwbin));
      IMpippim_IMnpim_wK0_n_wbin[iwbin][iq]= (TH2F*)fr[iq]->Get(Form("IMpippim_IMnpim_wK0_n_wbin%d",iwbin));
      IMnpim_IMnpip_wK0orwSid_n_wbin[iwbin][iq] = (TH2F*)fr[iq]->Get(Form("IMnpim_IMnpip_wK0orwSid_n_wbin%d",iwbin));
      IMpippim_IMnpip_wK0orwSid_n_wbin[iwbin][iq] = (TH2F*)fr[iq]->Get(Form("IMpippim_IMnpip_wK0orwSid_n_wbin%d",iwbin));
      IMpippim_IMnpim_wK0orwSid_n_wbin[iwbin][iq] = (TH2F*)fr[iq]->Get(Form("IMpippim_IMnpim_wK0orwSid_n_wbin%d",iwbin));
    }
  }

  //correlation plots excluding overlap region
  TH2F* IMnpim_IMnpip_wK0_woSid_n_wbin[nwbin][nqcut];
  TH2F* IMnpim_IMnpip_woK0_wSid_n_woSm_wbin[nwbin][nqcut];//
  TH2F* IMnpim_IMnpip_woK0_wSid_n_woSp_wbin[nwbin][nqcut];//
  TH2F* MMnmiss_IMpippim_wK0_woSid_n_wbin[nwbin][nqcut];
  TH2F* MMnmiss_IMnpip_woK0_wSid_n_woSm_wbin[nwbin][nqcut];
  TH2F* MMnmiss_IMnpim_woK0_wSid_n_woSp_wbin[nwbin][nqcut];
  TH2F* q_nmom_wK0_woSid_n_wbin[nwbin][nqcut];
  TH2F* q_nmom_woK0_wSid_n_woSm_wbin[nwbin][nqcut];
  TH2F* q_nmom_woK0_wSid_n_woSp_wbin[nwbin][nqcut];
  TH2F* IMpippim_IMnpip_woK0_wSid_n_woSp_wbin[nwbin][nqcut];//for Sm template
  TH2F* IMpippim_IMnpip_woK0_wSid_n_woSm_wbin[nwbin][nqcut];//for Sp template
  TH2F* IMpippim_IMnpim_woK0_wSid_n_woSp_wbin[nwbin][nqcut];
  TH2F* IMpippim_IMnpim_woK0_wSid_n_woSm_wbin[nwbin][nqcut];
  TH2F* IMpippim_IMnpip_wK0_woSid_n_wbin[nwbin][nqcut];
  TH2F* IMpippim_IMnpim_wK0_woSid_n_wbin[nwbin][nqcut];
  TH2F* q_IMnpipi_wK0_woSid_n[nqcut];
  TH2F* q_IMnpipi_woK0_wSid_n_woSp[nqcut];
  TH2F* q_IMnpipi_woK0_wSid_n_woSm[nqcut];
  for(int iq=0;iq<nqcut;iq++){
    for(int iwbin=0;iwbin<nwbin;iwbin++){
      IMnpim_IMnpip_wK0_woSid_n_wbin[iwbin][iq] = (TH2F*)fr[iq]->Get(Form("IMnpim_IMnpip_wK0_woSid_n_wbin%d",iwbin));
      IMnpim_IMnpip_woK0_wSid_n_woSm_wbin[iwbin][iq] = (TH2F*)fr[iq]->Get(Form("IMnpim_IMnpip_woK0_wSid_n_woSm_wbin%d",iwbin));
      IMnpim_IMnpip_woK0_wSid_n_woSp_wbin[iwbin][iq] = (TH2F*)fr[iq]->Get(Form("IMnpim_IMnpip_woK0_wSid_n_woSp_wbin%d",iwbin));
      MMnmiss_IMpippim_wK0_woSid_n_wbin[iwbin][iq] = (TH2F*)fr[iq]->Get(Form("MMnmiss_IMpippim_wK0_woSid_n_wbin%d",iwbin));
      MMnmiss_IMnpip_woK0_wSid_n_woSm_wbin[iwbin][iq] = (TH2F*)fr[iq]->Get(Form("MMnmiss_IMnpip_woK0_wSid_n_woSm_wbin%d",iwbin));
      MMnmiss_IMnpim_woK0_wSid_n_woSp_wbin[iwbin][iq] = (TH2F*)fr[iq]->Get(Form("MMnmiss_IMnpim_woK0_wSid_n_woSp_wbin%d",iwbin));
      q_nmom_wK0_woSid_n_wbin[iwbin][iq] = (TH2F*)fr[iq]->Get(Form("q_nmom_wK0_woSid_n_wbin%d",iwbin));
      q_nmom_woK0_wSid_n_woSm_wbin[iwbin][iq] = (TH2F*)fr[iq]->Get(Form("q_nmom_woK0_wSid_n_woSm_wbin%d",iwbin));
      q_nmom_woK0_wSid_n_woSp_wbin[iwbin][iq] = (TH2F*)fr[iq]->Get(Form("q_nmom_woK0_wSid_n_woSp_wbin%d",iwbin));
      IMpippim_IMnpip_woK0_wSid_n_woSp_wbin[iwbin][iq] = (TH2F*)fr[iq]->Get(Form("IMpippim_IMnpip_woK0_wSid_n_woSp_wbin%d",iwbin));
      IMpippim_IMnpip_woK0_wSid_n_woSm_wbin[iwbin][iq] = (TH2F*)fr[iq]->Get(Form("IMpippim_IMnpip_woK0_wSid_n_woSm_wbin%d",iwbin));
      IMpippim_IMnpim_woK0_wSid_n_woSp_wbin[iwbin][iq] = (TH2F*)fr[iq]->Get(Form("IMpippim_IMnpim_woK0_wSid_n_woSp_wbin%d",iwbin));
      IMpippim_IMnpim_woK0_wSid_n_woSm_wbin[iwbin][iq] = (TH2F*)fr[iq]->Get(Form("IMpippim_IMnpim_woK0_wSid_n_woSm_wbin%d",iwbin));
      IMpippim_IMnpip_wK0_woSid_n_wbin[iwbin][iq] = (TH2F*)fr[iq]->Get(Form("IMpippim_IMnpip_wK0_woSid_n_wbin%d",iwbin));
      IMpippim_IMnpim_wK0_woSid_n_wbin[iwbin][iq] = (TH2F*)fr[iq]->Get(Form("IMpippim_IMnpim_wK0_woSid_n_wbin%d",iwbin));
    }
    q_IMnpipi_wK0_woSid_n[iq] = (TH2F*)fr[iq]->Get("q_IMnpipi_wK0_woSid_n");
    q_IMnpipi_woK0_wSid_n_woSp[iq] = (TH2F*)fr[iq]->Get("q_IMnpipi_woK0_wSid_n_woSp");
    q_IMnpipi_woK0_wSid_n_woSm[iq] = (TH2F*)fr[iq]->Get("q_IMnpipi_woK0_wSid_n_woSm");
  }
  
  //for K0  (missing n is always selected)
  TH1D* IMnpim_wK0_woSid[nwbin][nqcut];
  TH1D* IMnpip_wK0_woSid[nwbin][nqcut];
  TH1D* MMnmiss_wK0_woSid[nwbin][nqcut];
  TH1D* IMpippim_wK0_woSid[nwbin][nqcut];
  TH1D* q_wK0_woSid[nwbin][nqcut];
  TH1D* nmom_wK0_woSid[nwbin][nqcut];
  TH1D* IMnpipi_wK0_woSid[nqcut];
  //for Sp (missing n is always selected),excluding K0 and Sm
  TH1D* IMnpim_woK0_woSm[nwbin][nqcut];
  TH1D* IMnpip_woK0_woSm[nwbin][nqcut];
  TH1D* MMnmiss_woK0_woSm[nwbin][nqcut];
  TH1D* IMpippim_woK0_woSm[nwbin][nqcut];
  TH1D* q_woK0_woSm[nwbin][nqcut];
  TH1D* nmom_woK0_woSm[nwbin][nqcut];
  TH1D* IMnpipi_woK0_woSm[nqcut];
  //for Sm (missing n is always selected),excluding K0 and Sp
  TH1D* IMnpim_woK0_woSp[nwbin][nqcut];
  TH1D* IMnpip_woK0_woSp[nwbin][nqcut];
  TH1D* MMnmiss_woK0_woSp[nwbin][nqcut];
  TH1D* IMpippim_woK0_woSp[nwbin][nqcut];
  TH1D* q_woK0_woSp[nwbin][nqcut];
  TH1D* nmom_woK0_woSp[nwbin][nqcut];
  TH1D* IMnpipi_woK0_woSp[nqcut];
  



  //projection hists
  TH1D* IMnpip_wK0[10][nqcut];
  TH1D* IMnpip_wK0_select[10][nqcut];
  TH1D* IMnpip_wK0_lohi[10][nqcut];
  TH1D* IMpippim_Sp[10][nqcut];
  TH1D* IMpippim_Sp_select[10][nqcut];
  TH1D* IMpippim_Sp_lohi[10][nqcut];
  TH1D* IMnpim_wK0[10][nqcut];
  TH1D* IMnpim_wK0_select[10][nqcut];
  TH1D* IMnpim_wK0_lohi[10][nqcut];
  TH1D* IMpippim_Sm[10][nqcut];
  TH1D* IMpippim_Sm_select[10][nqcut];
  TH1D* IMpippim_Sm_lohi[10][nqcut];
  TH1D* IMnpip_Sm[10][nqcut];
  TH1D* IMnpip_Sm_select[10][nqcut];
  TH1D* IMnpip_Sm_lohi[10][nqcut];
  TH1D* IMnpim_Sp[10][nqcut];
  TH1D* IMnpim_Sp_select[10][nqcut];
  TH1D* IMnpim_Sp_lohi[10][nqcut];
  
  TCanvas *cIMpippim_IMnpip[10][nqcut];
  TCanvas *cIMpippim_IMnpim[10][nqcut];
  TCanvas *cIMnpim_IMnpip[10][nqcut];
  
  TH1D* IMnpipi_Sp_ratio_wK0[nqcut];
  TH1D* IMnpipi_Sp_ratio_wSm[nqcut];
  TH1D* IMnpipi_Sm_ratio_wK0[nqcut];
  TH1D* IMnpipi_Sm_ratio_wSp[nqcut];
  TH1D* IMnpipi_K0_ratio_wSp[nqcut];
  TH1D* IMnpipi_K0_ratio_wSm[nqcut];
  for(int iq=0;iq<nqcut;iq++){
    IMnpipi_Sp_ratio_wK0[iq] = new TH1D(Form("IMnpipi_Sp_ratio_wK0_%d",iq),Form("IMnpipi_Sp_ratio_wK0_%s",cqcut[iq]),100,1,2);
    IMnpipi_Sp_ratio_wSm[iq] = new TH1D(Form("IMnpipi_Sp_ratio_wSm_%d",iq),Form("IMnpipi_Sp_ratio_wSm_%s",cqcut[iq]),100,1,2);
    IMnpipi_Sm_ratio_wK0[iq] = new TH1D(Form("IMnpipi_Sm_ratio_wK0_%d",iq),Form("IMnpipi_Sm_ratio_wK0_%s",cqcut[iq]),100,1,2);
    IMnpipi_Sm_ratio_wSp[iq] = new TH1D(Form("IMnpipi_Sm_ratio_wSp_%d",iq),Form("IMnpipi_Sm_ratio_wSp_%s",cqcut[iq]),100,1,2);
    IMnpipi_K0_ratio_wSp[iq] = new TH1D(Form("IMnpipi_K0_ratio_wSp_%d",iq),Form("IMnpipi_K0_ratio_wSp_%s",cqcut[iq]),100,1,2);
    IMnpipi_K0_ratio_wSm[iq] = new TH1D(Form("IMnpipi_K0_ratio_wSm_%d",iq),Form("IMnpipi_K0_ratio_wSm_%s",cqcut[iq]),100,1,2);
  }

  std::cout << __LINE__ << std::endl;
  
  for(int imerge=1;imerge<nwbin;imerge++){
    for(int iq=0;iq<nqcut;iq++){
      IMnpip_wK0[imerge][iq] = (TH1D*)IMpippim_IMnpip_wK0_n_wbin[imerge][iq]->ProjectionX(Form("IMnpip_wK0%d_%d",imerge,iq));
      IMpippim_Sp[imerge][iq] = (TH1D*)IMpippim_IMnpip_wSid_n_Sp_wbin[imerge][iq]->ProjectionY(Form("IMpippim_Sp%d_%d",imerge,iq));
      IMnpim_wK0[imerge][iq] = (TH1D*)IMpippim_IMnpim_wK0_n_wbin[imerge][iq]->ProjectionX(Form("IMnpim_wK0%d_%d",imerge,iq));
      IMpippim_Sm[imerge][iq] = (TH1D*)IMpippim_IMnpim_wSid_n_Sm_wbin[imerge][iq]->ProjectionY(Form("IMpippim_Sm%d_%d",imerge,iq));
      IMnpip_Sm[imerge][iq] = (TH1D*)IMnpim_IMnpip_wSid_n_Sm_wbin[imerge][iq]->ProjectionX(Form("IMnpip_Sm%d_%d",imerge,iq));
      IMnpim_Sp[imerge][iq] = (TH1D*)IMnpim_IMnpip_wSid_n_Sp_wbin[imerge][iq]->ProjectionY(Form("IMnpim_Sp%d_%d",imerge,iq));

      std::cout << __LINE__ << std::endl;
      //Draw in Canvas 
      cIMpippim_IMnpip[imerge][iq] = new TCanvas(Form("cIMpippim_IMnpip%d_%d",imerge,iq),Form("cIMpippim_IMnpip%d_%d",imerge,iq));
      cIMpippim_IMnpip[imerge][iq]->Divide(2,2);
      cIMpippim_IMnpip[imerge][iq]->cd(3);
      if(RemoveNegative)IMpippim_IMnpip_wK0orwSid_n_wbin[imerge][iq]->SetMinimum(0);
      IMpippim_IMnpip_wK0orwSid_n_wbin[imerge][iq]->Draw("colz");
      cIMpippim_IMnpip[imerge][iq]->cd(1);
      IMnpip_wK0[imerge][iq]->Draw("HE");
      IMnpip_wK0_select[imerge][iq] = (TH1D*)IMnpip_wK0[imerge][iq]->Clone(Form("IMnpip_wK0_select_%d_%d",imerge,iq));
      IMnpip_wK0_select[imerge][iq]->SetLineColor(2);
      IMnpip_wK0_select[imerge][iq]->GetXaxis()->SetRangeUser(
          IMnpip_wK0_select[imerge][iq]->GetBinLowEdge(9),
          IMnpip_wK0_select[imerge][iq]->GetBinLowEdge(10)
          );
      IMnpip_wK0_lohi[imerge][iq] = (TH1D*)IMnpip_wK0[imerge][iq]->Clone(Form("IMnpip_wK0_lohi_%d_%d",imerge,iq));
      IMnpip_wK0_lohi[imerge][iq]->SetBinContent(9,0);
      IMnpip_wK0_lohi[imerge][iq]->SetBinError(9,0);
      IMnpip_wK0_lohi[imerge][iq]->SetLineColor(4);
      IMnpip_wK0_lohi[imerge][iq]->Draw("HEsame");
      IMnpip_wK0_select[imerge][iq]->Draw("HEsame");
      cIMpippim_IMnpip[imerge][iq]->cd(4);
      IMpippim_Sp[imerge][iq]->Draw("HE");
      IMpippim_Sp_select[imerge][iq] = (TH1D*)IMpippim_Sp[imerge][iq]->Clone(Form("IMpippim_Sp_select_%d_%d",imerge,iq));
      IMpippim_Sp_select[imerge][iq]->SetLineColor(2);
      IMpippim_Sp_select[imerge][iq]->GetXaxis()->SetRangeUser(
          IMpippim_Sp_select[imerge][iq]->GetBinLowEdge(10),
          IMpippim_Sp_select[imerge][iq]->GetBinLowEdge(11));
      IMpippim_Sp_select[imerge][iq]->Draw("HEsame");
      IMpippim_Sp_lohi[imerge][iq] = (TH1D*)IMpippim_Sp[imerge][iq]->Clone(Form("IMpippim_Sp_lohi_%d_%d",imerge,iq));
      IMpippim_Sp_lohi[imerge][iq]->SetBinContent(10,0);
      IMpippim_Sp_lohi[imerge][iq]->SetBinError(10,0);
      IMpippim_Sp_lohi[imerge][iq]->SetLineColor(4);
      IMpippim_Sp_lohi[imerge][iq]->Draw("HEsame");

      cIMpippim_IMnpip[imerge][iq]->cd(2);

      std::cout << __LINE__ << std::endl;

      //Sigma- and K0 overlap 
      cIMpippim_IMnpim[imerge][iq] = new TCanvas(Form("cIMpippim_IMnpim%d_%d",imerge,iq),Form("cIMpippim_IMnpim%d_%d",imerge,iq));
      cIMpippim_IMnpim[imerge][iq]->Divide(2,2);
      cIMpippim_IMnpim[imerge][iq]->cd(3);
      if(RemoveNegative)IMpippim_IMnpim_wK0orwSid_n_wbin[imerge][iq]->SetMinimum(0);
      IMpippim_IMnpim_wK0orwSid_n_wbin[imerge][iq]->Draw("colz");
      cIMpippim_IMnpim[imerge][iq]->cd(1);
      IMnpim_wK0[imerge][iq]->Draw("HE");
      IMnpim_wK0_select[imerge][iq] = (TH1D*)IMnpim_wK0[imerge][iq]->Clone(Form("IMnpim_wK0_select_%d_%d",imerge,iq));
      IMnpim_wK0_select[imerge][iq]->SetLineColor(2);
      IMnpim_wK0_select[imerge][iq]->GetXaxis()->SetRangeUser(
          IMnpim_wK0_select[imerge][iq]->GetBinLowEdge(9),
          IMnpim_wK0_select[imerge][iq]->GetBinLowEdge(10)); 
      IMnpim_wK0_select[imerge][iq]->Draw("HEsame");
      IMnpim_wK0_lohi[imerge][iq] = (TH1D*)IMnpim_wK0[imerge][iq]->Clone(Form("IMnpim_wK0_lohi_%d_%d",imerge,iq));
      IMnpim_wK0_lohi[imerge][iq]->SetBinContent(9,0);
      IMnpim_wK0_lohi[imerge][iq]->SetBinError(9,0);
      IMnpim_wK0_lohi[imerge][iq]->SetLineColor(4);
      IMnpim_wK0_lohi[imerge][iq]->Draw("HEsame");
      cIMpippim_IMnpim[imerge][iq]->cd(4);
      IMpippim_Sm[imerge][iq]->Draw("HE");
      IMpippim_Sm_select[imerge][iq] = (TH1D*)IMpippim_Sm[imerge][iq]->Clone(Form("IMpippim_Sm_select_%d_%d",imerge,iq));
      IMpippim_Sm_select[imerge][iq]->SetLineColor(2);
      IMpippim_Sm_select[imerge][iq]->GetXaxis()->SetRangeUser(
            IMpippim_Sm_select[imerge][iq]->GetBinLowEdge(10),
            IMpippim_Sm_select[imerge][iq]->GetBinLowEdge(11));
      IMpippim_Sm_select[imerge][iq]->Draw("HEsame");
      IMpippim_Sm_lohi[imerge][iq] = (TH1D*)IMpippim_Sm[imerge][iq]->Clone(Form("IMpippim_Sm_lohi_%d_%d",imerge,iq));
      IMpippim_Sm_lohi[imerge][iq]->SetLineColor(4);
      IMpippim_Sm_lohi[imerge][iq]->SetBinContent(10,0);
      IMpippim_Sm_lohi[imerge][iq]->SetBinError(10,0);
      IMpippim_Sm_lohi[imerge][iq]->Draw("HEsame");
      

      cIMnpim_IMnpip[imerge][iq] = new TCanvas(Form("cIMnpim_IMnpip%d_%d",imerge,iq),Form("cIMnpim_IMnpip%d_%d",imerge,iq));
      cIMnpim_IMnpip[imerge][iq]->Divide(2,2);
      cIMnpim_IMnpip[imerge][iq]->cd(3);
      if(RemoveNegative)IMnpim_IMnpip_wK0orwSid_n_wbin[imerge][iq]->SetMinimum(0);
      IMnpim_IMnpip_wK0orwSid_n_wbin[imerge][iq]->GetXaxis()->SetRangeUser(1,1.4);
      IMnpim_IMnpip_wK0orwSid_n_wbin[imerge][iq]->GetYaxis()->SetRangeUser(1,1.4);
      IMnpim_IMnpip_wK0orwSid_n_wbin[imerge][iq]->Draw("colz");
      cIMnpim_IMnpip[imerge][iq]->cd(1);
      IMnpip_Sm[imerge][iq]->GetXaxis()->SetRangeUser(1,1.4);
      IMnpip_Sm[imerge][iq]->SetMinimum(0);
      IMnpip_Sm[imerge][iq]->Draw("HE");
      IMnpip_Sm_select[imerge][iq] = (TH1D*)IMnpip_Sm[imerge][iq]->Clone(Form("IMnpip_Sm_select_%d_%d",imerge,iq));
      IMnpip_Sm_select[imerge][iq]->SetLineColor(2);
      IMnpip_Sm_select[imerge][iq]->GetXaxis()->SetRangeUser(
            IMnpip_Sm_select[imerge][iq]->GetBinLowEdge(9),
            IMnpip_Sm_select[imerge][iq]->GetBinLowEdge(10));
      IMnpip_Sm_select[imerge][iq]->Draw("HEsame");
      IMnpip_Sm_lohi[imerge][iq] = (TH1D*)IMnpip_Sm[imerge][iq]->Clone(Form("IMnpip_Sm_lohi_%d_%d",imerge,iq));
      IMnpip_Sm_lohi[imerge][iq]->SetLineColor(4);
      IMnpip_Sm_lohi[imerge][iq]->SetBinContent(9,0);
      IMnpip_Sm_lohi[imerge][iq]->SetBinError(9,0);
      IMnpip_Sm_lohi[imerge][iq]->Draw("HEsame");

      cIMnpim_IMnpip[imerge][iq]->cd(4);
      IMnpim_Sp[imerge][iq]->GetXaxis()->SetRangeUser(1,1.4);
      IMnpim_Sp[imerge][iq]->SetMinimum(0);
      IMnpim_Sp[imerge][iq]->Draw("HE");
      IMnpim_Sp_select[imerge][iq] = (TH1D*)IMnpim_Sp[imerge][iq]->Clone(Form("IMnpim_Sp_select_%d_%d",imerge,iq));
      IMnpim_Sp_select[imerge][iq]->SetLineColor(2);
      IMnpim_Sp_select[imerge][iq]->GetXaxis()->SetRangeUser(
            IMnpim_Sp_select[imerge][iq]->GetBinLowEdge(9),
            IMnpim_Sp_select[imerge][iq]->GetBinLowEdge(10));
      
      IMnpim_Sp_select[imerge][iq]->Draw("HEsame");
      IMnpim_Sp_lohi[imerge][iq] = (TH1D*)IMnpim_Sp[imerge][iq]->Clone(Form("IMnpim_Sp_lohi_%d_%d",imerge,iq));
      IMnpim_Sp_lohi[imerge][iq]->SetLineColor(4);
      IMnpim_Sp_lohi[imerge][iq]->SetBinContent(9,0);
      IMnpim_Sp_lohi[imerge][iq]->SetBinError(9,0);
      IMnpim_Sp_lohi[imerge][iq]->Draw("HEsame");
    }
  }
  
  TCanvas *cIMnpipi_Sp_ratio_wK0[nqcut];
  TCanvas *cIMnpipi_K0_ratio_wSp[nqcut];
  TCanvas *cIMnpipi_Sm_ratio_wK0[nqcut];
  TCanvas *cIMnpipi_K0_ratio_wSm[nqcut];
  TCanvas *cIMnpipi_Sp_ratio_wSm[nqcut];
  TCanvas *cIMnpipi_Sm_ratio_wSp[nqcut];
  
  for(int iq=0;iq<nqcut;iq++){
    cIMnpipi_Sp_ratio_wK0[iq] = new TCanvas(Form("cIMnpipi_Sp_ratio_wK0_%d",iq),Form("cIMnpipi_Sp_ratio_wK0_%d",iq));
    IMnpipi_Sp_ratio_wK0[iq]->Draw("HIST");

    cIMnpipi_K0_ratio_wSp[iq] = new TCanvas(Form("cIMnpipi_K0_ratio_wSp_%d",iq),Form("cIMnpipi_K0_ratio_wSp_%d",iq));
    IMnpipi_K0_ratio_wSp[iq]->Draw("HIST");
    
    cIMnpipi_Sm_ratio_wK0[iq] = new TCanvas(Form("cIMnpipi_Sm_ratio_wK0_%d",iq),Form("cIMnpipi_Sm_ratio_wK0_%d",iq));
    IMnpipi_Sm_ratio_wK0[iq]->Draw("HIST");

    cIMnpipi_K0_ratio_wSm[iq] = new TCanvas(Form("cIMnpipi_K0_ratio_wSm_%d",iq),Form("cIMnpipi_K0_ratio_wSm_%d",iq));
    IMnpipi_K0_ratio_wSm[iq]->Draw("HIST");

    cIMnpipi_Sp_ratio_wSm[iq] = new TCanvas(Form("cIMnpipi_Sp_ratio_wSm_%d",iq),Form("cIMnpipi_Sp_ratio_wSm_%d",iq));
    IMnpipi_Sp_ratio_wSm[iq]->Draw("HIST");

    cIMnpipi_Sm_ratio_wSp[iq] = new TCanvas(Form("cIMnpipi_Sm_ratio_wSp_%d",iq),Form("cIMnpipi_Sm_ratio_wSp_%d",iq));
    IMnpipi_Sm_ratio_wSp[iq]->Draw("HIST");
  }


  TCanvas *c = NULL;
  TSeqCollection *SCol = gROOT->GetListOfCanvases();
  int size = SCol->GetSize();
  TIter next(SCol);
  TString pdfname = "Decomposition.pdf";
  if(RebinMode) pdfname = "Decomposition_Rebin.pdf";
  for(int i=0;i<size;i++){
    //pdf->NewPage();
    c= (TCanvas*)next();
    c->Draw();
    c->cd();
    //inside the canvas
    c->Modified();
    c->Update();
    //std::cout << c->GetName() << std::endl;
    //make 1 pdf file
    if(i==0) c->Print(pdfname+"(",Form("pdf Title:%s",c->GetTitle()));
    else if(i==size-1)c->Print(pdfname+")",Form("pdf Title:%s",c->GetTitle())); 
    else c->Print(pdfname,Form("pdf Title:%s",c->GetTitle())); 
    
    //make separated pdf files
    //c->Print(Form("pdf/%s.pdf",c->GetTitle()));
  }
  
  TIter nexthist2(gDirectory->GetList());
  TFile *fout;
  if(!RebinMode) fout = new TFile("Decomposition.root","RECREATE");
  else           fout = new TFile("Decomposition_Rebin.root","RECREATE");
  fout->Print();
  fout->cd();
  TObject *obj = nullptr;
  while( (obj = (TObject*)nexthist2())!=nullptr) {
    obj->Write();
  }
  fout->cd();
  fout->Close();

  return;
  

}
*/
