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

void HistToRorateGraph(TH1D* h1, TGraphErrors &gr)
{
  for(int ibin=1;ibin<=h1->GetNbinsX();ibin++){
    double cont = h1->GetBinContent(ibin);
    double err  = h1->GetBinError(ibin);
    double bincenter = h1->GetBinCenter(ibin);
    gr.SetPoint(ibin,cont,bincenter);
    gr.SetPointError(ibin,err,0);
    if(cont < -1000) std::cout << ibin << " " << cont << std::endl;
  }
  gr.GetYaxis()->SetRangeUser(h1->GetXaxis()->GetXmin(),h1->GetXaxis()->GetXmax());
  gr.GetYaxis()->SetLabelSize(0);
  gr.SetMarkerStyle(20);
  return;
}

void K0SigmaTemp()
{
  TFile *fr = TFile::Open("evanaIMpisigma_npippim_v202_out_iso.root");
  TFile *fmix = TFile::Open("evanaIMpisigma_npippim_v202_MIX_cut4_out_iso.root");
  TFile *frqhi = TFile::Open("evanaIMpisigma_npippim_v202_out_iso_qhi.root");
  TFile *fmixqhi = TFile::Open("evanaIMpisigma_npippim_v202_MIX_cut4_out_iso_qhi.root");
  TFile *frqlo = TFile::Open("evanaIMpisigma_npippim_v202_out_iso_qlo.root");
  TFile *fmixqlo = TFile::Open("evanaIMpisigma_npippim_v202_MIX_cut4_out_iso_qlo.root");
  fr->Print() ;
  fmix->Print();
   
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0.);  

  const unsigned int nbintemplate = 50;
  TH2F* IMpippim_IMnpip_n_bin_data[nbintemplate];
  TH2F* IMpippim_IMnpim_n_bin_data[nbintemplate];
  TH2F* IMpippim_IMnpip_n_bin_qhi_data[nbintemplate];
  TH2F* IMpippim_IMnpim_n_bin_qhi_data[nbintemplate];
  TH2F* IMpippim_IMnpip_n_bin_qlo_data[nbintemplate];
  TH2F* IMpippim_IMnpim_n_bin_qlo_data[nbintemplate];
  TH2F* IMpippim_IMnpip_n_bin_mix[nbintemplate];
  TH2F* IMpippim_IMnpim_n_bin_mix[nbintemplate];
  TH2F* IMpippim_IMnpip_n_bin_qhi_mix[nbintemplate];
  TH2F* IMpippim_IMnpim_n_bin_qhi_mix[nbintemplate];
  TH2F* IMpippim_IMnpip_n_bin_qlo_mix[nbintemplate];
  TH2F* IMpippim_IMnpim_n_bin_qlo_mix[nbintemplate];
  TH2F* IMpippim_IMnpip_n_woSm_bin_data[nbintemplate];
  TH2F* IMpippim_IMnpim_n_woSp_bin_data[nbintemplate];
  TH2F* IMpippim_IMnpip_n_woSm_bin_mix[nbintemplate];
  TH2F* IMpippim_IMnpim_n_woSp_bin_mix[nbintemplate];
  TH2F* IMpippim_IMnpip_n_woSm_bin_qhi_data[nbintemplate];
  TH2F* IMpippim_IMnpim_n_woSp_bin_qhi_data[nbintemplate];
  TH2F* IMpippim_IMnpip_n_woSm_bin_qhi_mix[nbintemplate];
  TH2F* IMpippim_IMnpim_n_woSp_bin_qhi_mix[nbintemplate];
  TH2F* IMpippim_IMnpip_n_woSm_bin_qlo_data[nbintemplate];
  TH2F* IMpippim_IMnpim_n_woSp_bin_qlo_data[nbintemplate];
  TH2F* IMpippim_IMnpip_n_woSm_bin_qlo_mix[nbintemplate];
  TH2F* IMpippim_IMnpim_n_woSp_bin_qlo_mix[nbintemplate];
  TH2F* IMpippim_IMnpip_n_Sp_bin_data[nbintemplate];
  TH2F* IMpippim_IMnpim_n_Sm_bin_data[nbintemplate];
  TH2F* IMpippim_IMnpip_n_Sp_bin_mix[nbintemplate];
  TH2F* IMpippim_IMnpim_n_Sm_bin_mix[nbintemplate];
  TH2F* IMpippim_IMnpip_n_Sp_bin_qhi_data[nbintemplate];
  TH2F* IMpippim_IMnpim_n_Sm_bin_qhi_data[nbintemplate];
  TH2F* IMpippim_IMnpip_n_Sp_bin_qhi_mix[nbintemplate];
  TH2F* IMpippim_IMnpim_n_Sm_bin_qhi_mix[nbintemplate];
  TH2F* IMpippim_IMnpip_n_Sp_bin_qlo_data[nbintemplate];
  TH2F* IMpippim_IMnpim_n_Sm_bin_qlo_data[nbintemplate];
  TH2F* IMpippim_IMnpip_n_Sp_bin_qlo_mix[nbintemplate];
  TH2F* IMpippim_IMnpim_n_Sm_bin_qlo_mix[nbintemplate];
  
  for(unsigned int ibin=0;ibin<nbintemplate;ibin++){
    IMpippim_IMnpip_n_bin_data[ibin] = (TH2F*)fr->Get(Form("IMpippim_IMnpip_n_bin%d",ibin));
    IMpippim_IMnpip_n_bin_mix[ibin] = (TH2F*)fmix->Get(Form("IMpippim_IMnpip_n_bin%d",ibin));
    IMpippim_IMnpip_n_bin_qhi_data[ibin] = (TH2F*)frqhi->Get(Form("IMpippim_IMnpip_n_bin%d",ibin));
    IMpippim_IMnpip_n_bin_qhi_mix[ibin] = (TH2F*)fmixqhi->Get(Form("IMpippim_IMnpip_n_bin%d",ibin));
    IMpippim_IMnpip_n_bin_qlo_data[ibin] = (TH2F*)frqlo->Get(Form("IMpippim_IMnpip_n_bin%d",ibin));
    IMpippim_IMnpip_n_bin_qlo_mix[ibin] = (TH2F*)fmixqlo->Get(Form("IMpippim_IMnpip_n_bin%d",ibin));
    IMpippim_IMnpip_n_Sp_bin_data[ibin] = (TH2F*)fr->Get(Form("IMpippim_IMnpip_n_Sp_bin%d",ibin));
    IMpippim_IMnpip_n_Sp_bin_mix[ibin] = (TH2F*)fmix->Get(Form("IMpippim_IMnpip_n_Sp_bin%d",ibin));
    IMpippim_IMnpip_n_Sp_bin_qhi_data[ibin] = (TH2F*)frqhi->Get(Form("IMpippim_IMnpip_n_Sp_bin%d",ibin));
    IMpippim_IMnpip_n_Sp_bin_qhi_mix[ibin] = (TH2F*)fmixqhi->Get(Form("IMpippim_IMnpip_n_Sp_bin%d",ibin));
    IMpippim_IMnpip_n_Sp_bin_qlo_data[ibin] = (TH2F*)frqlo->Get(Form("IMpippim_IMnpip_n_Sp_bin%d",ibin));
    IMpippim_IMnpip_n_Sp_bin_qlo_mix[ibin] = (TH2F*)fmixqlo->Get(Form("IMpippim_IMnpip_n_Sp_bin%d",ibin));
    IMpippim_IMnpim_n_Sm_bin_data[ibin] = (TH2F*)fr->Get(Form("IMpippim_IMnpim_n_Sm_bin%d",ibin));
    IMpippim_IMnpim_n_Sm_bin_mix[ibin] = (TH2F*)fmix->Get(Form("IMpippim_IMnpim_n_Sm_bin%d",ibin));
    IMpippim_IMnpim_n_Sm_bin_qhi_data[ibin] = (TH2F*)frqhi->Get(Form("IMpippim_IMnpim_n_Sm_bin%d",ibin));
    IMpippim_IMnpim_n_Sm_bin_qhi_mix[ibin] = (TH2F*)fmixqhi->Get(Form("IMpippim_IMnpim_n_Sm_bin%d",ibin));
    IMpippim_IMnpim_n_Sm_bin_qlo_data[ibin] = (TH2F*)frqlo->Get(Form("IMpippim_IMnpim_n_Sm_bin%d",ibin));
    IMpippim_IMnpim_n_Sm_bin_qlo_mix[ibin] = (TH2F*)fmixqlo->Get(Form("IMpippim_IMnpim_n_Sm_bin%d",ibin));
  }

  //subtract

  for(unsigned int ibin=0;ibin<nbintemplate;ibin++){

  
  }

  TCanvas *c = NULL;
  TSeqCollection *SCol = gROOT->GetListOfCanvases();
  int size = SCol->GetSize();
  TIter next(SCol);
  TString pdfname = "disp_mixcomp.pdf";
  for(int i=0;i<size;i++){
    //pdf->NewPage();
    c= (TCanvas*)next();
    c->Draw();
    c->cd();
    //inside the canvas
    //TPaveText *pt = new TPaveText(.74,.81,0.9,0.90,"NDC");
    c->Modified();
    c->Update();
    std::cout << c->GetName() << std::endl;
    //make 1 pdf file
    if(i==0) c->Print(pdfname+"(",Form("pdf Title:%s",c->GetTitle()));
    else if(i==size-1)c->Print(pdfname+")",Form("pdf Title:%s",c->GetTitle())); 
    else c->Print(pdfname,Form("pdf Title:%s",c->GetTitle())); 
    
    //make separated pdf files
    //c->Print(Form("pdf/%s.pdf",c->GetTitle()));
  }
  



  return;
  

}
