#include "anacuts.h"

void disp_IMsigma_h2(const int dEcut=2,const int sysud=0,const int simmode=0)
{
  gStyle->SetOptStat(0);
  bool SimSpMode=false;
  bool SimSmMode=false;
  bool SimK0Mode=false;//K0n
  bool SimLppMode=false;//Lambda pi+ pi-
  
  if(simmode==0){

  }else if(simmode==1){
    SimSpMode=true;
  }else if(simmode==2){
    SimSmMode=true;
  }else if(simmode==3){
    SimK0Mode=true;
  }else if(simmode==4){
    SimLppMode=true;
  }
  
  const int version = 14;
  TFile *f;
  TFile *fr;
  TFile *fmix;
  if(!SimSpMode && !SimSmMode){
    f = TFile::Open(Form("evanaIMsigma_npi_h2_v%d_out_dE%d_iso_nostop_sub_sys%d.root",version,dEcut,sysud));
    fr = TFile::Open(Form("evanaIMsigma_npi_h2_v%d_out_dE%d_iso_nostop.root",version,dEcut));
    fmix = TFile::Open(Form("evanaIMsigma_npi_h2_v%d_MIX_out_dE%d_iso_nostop_sys%d.root",version,dEcut,sysud));
  }else if(SimSpMode){
    f = TFile::Open(Form("../simpost/simIMsigma_H2_Sppim_npi_v11_out_iso_dE%d_rej_nostop.root",dEcut));
  }else if(SimSmMode){
    f = TFile::Open(Form("../simpost/simIMsigma_H2_Smpip_npi_v11_out_iso_dE%d_rej_nostop.root",dEcut));
  }else if(SimK0Mode){
    f = TFile::Open(Form("../simpost/simIMsigma_H2_K0n_npi_v2_out_iso_dE%d_rej_nostop.root",dEcut));
  }else if(SimLppMode){
    f = TFile::Open(Form("../simpost/simIMsigma_H2_pipiL_npi_v1_out_iso_dE%d_rej_nostop.root",dEcut));
  }
  
  if(!SimSpMode && !SimSmMode){
    TCanvas *cvici_Sp = new TCanvas("cvici_Sp","cvici_Sp",1000,800);
    TH2F* MM2npi_IMnpip_vici = (TH2F*)fr->Get("MM2npi_IMnpip_vici");
    TH2F* MM2npi_IMnpip_vici_mix = (TH2F*)fmix->Get("MM2npi_IMnpip_vici");
    MM2npi_IMnpip_vici->RebinX(2);
    MM2npi_IMnpip_vici_mix->RebinX(2);
    //MM2npi_IMnpip_vici->GetXaxis()->SetRangeUser(1.1,1.3);
    //MM2npi_IMnpip_vici->GetYaxis()->SetRangeUser(-0.4,0.4);
    MM2npi_IMnpip_vici->Draw("colz");

    TCanvas *cvici_Sp_mix = new TCanvas("cvici_Sp_mix","cvici_Sp_mix",1000,800);
    MM2npi_IMnpip_vici_mix->Draw("colz");

    TCanvas *cvici_Sp_px = new TCanvas("cvici_Sp_px","cvici_Sp_px",1000,800);
    TH1D* IMnpip_vici = (TH1D*)MM2npi_IMnpip_vici->ProjectionX("IMnpip_vici");
    TH1D* IMnpip_vici_mix = (TH1D*)MM2npi_IMnpip_vici_mix->ProjectionX("IMnpip_vici_mix");

    IMnpip_vici->Draw("E");
    IMnpip_vici_mix->SetLineColor(2);
    IMnpip_vici_mix->Draw("HEsame");

    TCanvas *cvici_Sp_px_sum = new TCanvas("cvici_Sp_px_sum","cvici_Sp_px_sum",1000,800);
    TH1D* IMnpip_vici_sub = (TH1D*)IMnpip_vici->Clone("IMnpip_vici_sub");
    IMnpip_vici_sub->Add(IMnpip_vici_mix,-1.0);
    IMnpip_vici_sub->Draw("E");

    TCanvas *cvici_Sp_py = new TCanvas("cvici_Sp_py","cvici_Sp_py",1000,800);
    TH1D* MMnpip_vici = (TH1D*)MM2npi_IMnpip_vici->ProjectionY("MMnpip_vici");
    TH1D* MMnpip_vici_mix = (TH1D*)MM2npi_IMnpip_vici_mix->ProjectionY("MMnpip_vici_mix");

    MMnpip_vici->Draw("E");
    MMnpip_vici_mix->SetLineColor(2);
    MMnpip_vici_mix->Draw("HEsame");

    TCanvas *cvici_Sp_py_sum = new TCanvas("cvici_Sp_py_sum","cvici_Sp_py_sum",1000,800);
    TH1D* MMnpip_vici_sub = (TH1D*)MMnpip_vici->Clone("MMnpip_vici_sub");
    MMnpip_vici_sub->Add(MMnpip_vici_mix,-1.0);
    MMnpip_vici_sub->Draw("E");



    TCanvas *cvici_Sm = new TCanvas("cvici_Sm","cvici_Sm",1000,800);
    TH2F* MM2npi_IMnpim_vici = (TH2F*)fr->Get("MM2npi_IMnpim_vici");
    TH2F* MM2npi_IMnpim_vici_mix = (TH2F*)fmix->Get("MM2npi_IMnpim_vici");
    MM2npi_IMnpim_vici->RebinX(2);
    MM2npi_IMnpim_vici_mix->RebinX(2);
    //MM2npi_IMnpim_vici->GetXaxis()->SetRangeUser(1.1,1.3);
    //MM2npi_IMnpim_vici->GetYaxis()->SetRangeUser(-0.4,0.4);
    MM2npi_IMnpim_vici->Draw("colz");

    TCanvas *cvici_Sm_mix = new TCanvas("cvici_Sm_mix","cvici_Sm_mix",1000,800);
    MM2npi_IMnpim_vici_mix->Draw("colz");

    TCanvas *cvici_Sm_px = new TCanvas("cvici_Sm_px","cvici_Sm_px",1000,800);
    TH1D* IMnpim_vici = (TH1D*)MM2npi_IMnpim_vici->ProjectionX("IMnpim_vici");
    TH1D* IMnpim_vici_mix = (TH1D*)MM2npi_IMnpim_vici_mix->ProjectionX("IMnpim_vici_mix");

    IMnpim_vici->Draw("E");
    IMnpim_vici_mix->SetLineColor(2);
    IMnpim_vici_mix->Draw("HEsame");

    TCanvas *cvici_Sm_px_sum = new TCanvas("cvici_Sm_px_sum","cvici_Sm_px_sum",1000,800);
    TH1D* IMnpim_vici_sub = (TH1D*)IMnpim_vici->Clone("IMnpim_vici_sub");
    IMnpim_vici_sub->Add(IMnpim_vici_mix,-1.0);
    IMnpim_vici_sub->Draw("E");

    TCanvas *cvici_Sm_py = new TCanvas("cvici_Sm_py","cvici_Sm_py",1000,800);
    TH1D* MMnpim_vici = (TH1D*)MM2npi_IMnpim_vici->ProjectionY("MMnpim_vici");
    TH1D* MMnpim_vici_mix = (TH1D*)MM2npi_IMnpim_vici_mix->ProjectionY("MMnpim_vici_mix");

    MMnpim_vici->Draw("E");
    MMnpim_vici_mix->SetLineColor(2);
    MMnpim_vici_mix->Draw("HEsame");

    TCanvas *cvici_Sm_py_sum = new TCanvas("cvici_Sm_py_sum","cvici_Sm_py_sum",1000,800);
    TH1D* MMnpim_vici_sub = (TH1D*)MMnpim_vici->Clone("MMnpim_vici_sub");
    MMnpim_vici_sub->Add(MMnpim_vici_mix,-1.0);
    MMnpim_vici_sub->Draw("E");
    
    TH2F* Cospicm_IMnpip_pi_r = (TH2F*)fr->Get("Cospicm_IMnpip_pi");
    Cospicm_IMnpip_pi_r->RebinX(5);
    Cospicm_IMnpip_pi_r->GetXaxis()->SetRangeUser(1,1.5);
    TH2F* Cospicm_IMnpip_pi_mix = (TH2F*)fmix->Get("Cospicm_IMnpip_pi");
    Cospicm_IMnpip_pi_mix->RebinX(5);
    Cospicm_IMnpip_pi_mix->GetXaxis()->SetRangeUser(1,1.5);
    
    TH2F* Cospicm_IMnpim_pi_r = (TH2F*)fr->Get("Cospicm_IMnpim_pi");
    Cospicm_IMnpim_pi_r->RebinX(5);
    Cospicm_IMnpim_pi_r->GetXaxis()->SetRangeUser(1,1.5);
    TH2F* Cospicm_IMnpim_pi_mix = (TH2F*)fmix->Get("Cospicm_IMnpim_pi");
    Cospicm_IMnpim_pi_mix->RebinX(5);
    Cospicm_IMnpim_pi_mix->GetXaxis()->SetRangeUser(1,1.5);
    
    const int ncosbin=6;
    int cosbin[7];
    cosbin[0] = Cospicm_IMnpip_pi_r->GetYaxis()->FindBin(1.0);
    cosbin[1] = Cospicm_IMnpip_pi_r->GetYaxis()->FindBin(0.9); 
    cosbin[2] = Cospicm_IMnpip_pi_r->GetYaxis()->FindBin(0.8); 
    cosbin[3] = Cospicm_IMnpip_pi_r->GetYaxis()->FindBin(0.7); 
    cosbin[4] = Cospicm_IMnpip_pi_r->GetYaxis()->FindBin(0.6); 
    cosbin[5] = Cospicm_IMnpip_pi_r->GetYaxis()->FindBin(0.5); 
    cosbin[6] = Cospicm_IMnpip_pi_r->GetYaxis()->FindBin(0.4); 

    TCanvas *ccos_Sp_rmix = new TCanvas("ccos_Sp_rmix","ccos_Sp_rmix",1000,800);
    ccos_Sp_rmix->Divide(3,2);
    TH1D* IMnpip_coscut_r[6];
    TH1D* IMnpip_coscut_mix[6];
    const int Splow = Cospicm_IMnpip_pi_r->GetXaxis()->FindBin(anacuts::Sigmap_MIN);
    const int Sphigh = Cospicm_IMnpip_pi_r->GetXaxis()->FindBin(anacuts::Sigmap_MAX);
    for(int icosbin=0;icosbin<ncosbin;icosbin++){
      IMnpip_coscut_r[icosbin] = (TH1D*)Cospicm_IMnpip_pi_r->ProjectionX(Form("IMnpip_coscut_r%d",icosbin),cosbin[icosbin+1],cosbin[icosbin]);
      IMnpip_coscut_mix[icosbin] = (TH1D*)Cospicm_IMnpip_pi_mix->ProjectionX(Form("IMnpip_coscut_mix%d",icosbin),cosbin[icosbin+1],cosbin[icosbin]);
      ccos_Sp_rmix->cd(icosbin+1);
      //std::cout << IMnpip_coscut[icosbin]->Integral(Splow,Sphigh) << std::endl;
      IMnpip_coscut_r[icosbin]->GetXaxis()->SetRangeUser(1.1,1.3);
      IMnpip_coscut_r[icosbin]->SetTitle(Form("cos%0.1f-cos%0.1f ", 1-icosbin*0.1, 0.9-icosbin*0.1));
      IMnpip_coscut_r[icosbin]->Draw("HE");
      IMnpip_coscut_mix[icosbin]->SetLineColor(2);
      IMnpip_coscut_mix[icosbin]->Draw("HEsame");
    }

    TCanvas *ccos_Sm_rmix = new TCanvas("ccos_Sm_rmix","ccos_Sm_rmix",1000,800);
    ccos_Sm_rmix->Divide(3,2);
    TH1D* IMnpim_coscut_r[6];
    TH1D* IMnpim_coscut_mix[6];
    const int Smlow = Cospicm_IMnpim_pi_r->GetXaxis()->FindBin(anacuts::Sigmam_MIN);
    const int Smhigh = Cospicm_IMnpim_pi_r->GetXaxis()->FindBin(anacuts::Sigmam_MAX);
    for(int icosbin=0;icosbin<ncosbin;icosbin++){
      IMnpim_coscut_r[icosbin] = (TH1D*)Cospicm_IMnpim_pi_r->ProjectionX(Form("IMnpim_coscut_r%d",icosbin),cosbin[icosbin+1],cosbin[icosbin]);
      IMnpim_coscut_mix[icosbin] = (TH1D*)Cospicm_IMnpim_pi_mix->ProjectionX(Form("IMnpim_coscut_mix%d",icosbin),cosbin[icosbin+1],cosbin[icosbin]);
      ccos_Sm_rmix->cd(icosbin+1);
      //std::cout << IMnpim_coscut[icosbin]->Integral(Smlow,Smhigh) << std::endl;
      IMnpim_coscut_r[icosbin]->GetXaxis()->SetRangeUser(1.1,1.3);
      IMnpim_coscut_r[icosbin]->SetTitle(Form("cos%0.1f-cos%0.1f ", 1-icosbin*0.1, 0.9-icosbin*0.1));
      IMnpim_coscut_r[icosbin]->Draw("HE");
      IMnpim_coscut_mix[icosbin]->SetLineColor(2);
      IMnpim_coscut_mix[icosbin]->Draw("HEsame");
    }

    TH2F* Cospicm_MM2npi_Sp_r = (TH2F*)fr->Get("Cospicm_MM2npi_Sp");
    TH2F* Cospicm_MM2npi_Sp_mix = (TH2F*)fmix->Get("Cospicm_MM2npi_Sp");

    TCanvas *ccos_Sp_MM2_rmix = new TCanvas("ccos_Sp_MM2_rmix","ccos_Sp_MM2_rmix",1000,800);
    ccos_Sp_MM2_rmix->Divide(3,2);
    TH1D* MM2npi_coscut_Sp_r[6];
    TH1D* MM2npi_coscut_Sp_mix[6];
    for(int icosbin=0;icosbin<ncosbin;icosbin++){
      MM2npi_coscut_Sp_r[icosbin] = (TH1D*)Cospicm_MM2npi_Sp_r->ProjectionX(Form("MM2npi_coscut_r%d_Sp",icosbin),cosbin[icosbin+1],cosbin[icosbin]);
      MM2npi_coscut_Sp_mix[icosbin] = (TH1D*)Cospicm_MM2npi_Sp_mix->ProjectionX(Form("MM2npi_coscut_mix%d_Sp",icosbin),cosbin[icosbin+1],cosbin[icosbin]);
      ccos_Sp_MM2_rmix->cd(icosbin+1);
      //std::cout << IMnpim_coscut[icosbin]->Integral(Smlow,Smhigh) << std::endl;
      MM2npi_coscut_Sp_r[icosbin]->Rebin(2);
      MM2npi_coscut_Sp_mix[icosbin]->Rebin(2);
      MM2npi_coscut_Sp_r[icosbin]->GetXaxis()->SetRangeUser(-0.4,0.4);
      MM2npi_coscut_Sp_r[icosbin]->SetTitle(Form("cos%0.1f-cos%0.1f ", 1-icosbin*0.1, 0.9-icosbin*0.1));
      MM2npi_coscut_Sp_r[icosbin]->Draw("HE");
      MM2npi_coscut_Sp_mix[icosbin]->SetLineColor(2);
      MM2npi_coscut_Sp_mix[icosbin]->Draw("HEsame");
    }

    TH2F* Cospicm_MM2npi_Sm_r = (TH2F*)fr->Get("Cospicm_MM2npi_Sm");
    TH2F* Cospicm_MM2npi_Sm_mix = (TH2F*)fmix->Get("Cospicm_MM2npi_Sm");

    TCanvas *ccos_Sm_MM2_rmix = new TCanvas("ccos_Sm_MM2_rmix","ccos_Sm_MM2_rmix",1000,800);
    ccos_Sm_MM2_rmix->Divide(3,2);
    TH1D* MM2npi_coscut_Sm_r[6];
    TH1D* MM2npi_coscut_Sm_mix[6];
    for(int icosbin=0;icosbin<ncosbin;icosbin++){
      MM2npi_coscut_Sm_r[icosbin] = (TH1D*)Cospicm_MM2npi_Sm_r->ProjectionX(Form("MM2npi_coscut_r%d_Sm",icosbin),cosbin[icosbin+1],cosbin[icosbin]);
      MM2npi_coscut_Sm_mix[icosbin] = (TH1D*)Cospicm_MM2npi_Sm_mix->ProjectionX(Form("MM2npi_coscut_mix%d_Sm",icosbin),cosbin[icosbin+1],cosbin[icosbin]);
      ccos_Sm_MM2_rmix->cd(icosbin+1);
      //std::cout << IMnpim_coscut[icosbin]->Integral(Smlow,Smhigh) << std::endl;
      MM2npi_coscut_Sm_r[icosbin]->Rebin(2);
      MM2npi_coscut_Sm_mix[icosbin]->Rebin(2);
      MM2npi_coscut_Sm_r[icosbin]->GetXaxis()->SetRangeUser(-0.4,0.4);
      MM2npi_coscut_Sm_r[icosbin]->SetTitle(Form("cos%0.1f-cos%0.1f ", 1-icosbin*0.1, 0.9-icosbin*0.1));
      MM2npi_coscut_Sm_r[icosbin]->Draw("HE");
      MM2npi_coscut_Sm_mix[icosbin]->SetLineColor(2);
      MM2npi_coscut_Sm_mix[icosbin]->Draw("HEsame");
    }

    TCanvas *ccos_purity_Sp = new TCanvas("ccos_purity_Sp","ccos_purity_Sp");

    TH1F* Cospicm_pi_Sp_S = (TH1F*)f->Get("Cospicm_pi_Sp");
    Cospicm_pi_Sp_S->RebinX(5);

    TH1F* Cospicm_pi_Sp_SN = (TH1F*)fr->Get("Cospicm_pi_Sp");
    Cospicm_pi_Sp_SN->RebinX(5);
    TH1F* Cospicm_purity_Sp = (TH1F*)Cospicm_pi_Sp_S->Clone("Cospicm_putiry");
    Cospicm_purity_Sp->SetTitle("purity S/(S+N) #Sigma^{+}#pi^{-}");
    Cospicm_purity_Sp->Divide(Cospicm_pi_Sp_S,Cospicm_pi_Sp_SN,1.,1.,"b");
    Cospicm_purity_Sp->Draw();
    
    TCanvas *ccos_purity_Sm = new TCanvas("ccos_purity_Sm","ccos_purity_Sm");
    TH1F* Cospicm_pi_Sm_S = (TH1F*)f->Get("Cospicm_pi_Sm");
    Cospicm_pi_Sm_S->RebinX(5);
    TH1F* Cospicm_pi_Sm_SN = (TH1F*)fr->Get("Cospicm_pi_Sm");
    Cospicm_pi_Sm_SN->RebinX(5);
    TH1F* Cospicm_purity_Sm = (TH1F*)Cospicm_pi_Sm_S->Clone("Cospicm_putiry");
    Cospicm_purity_Sm->SetTitle("purity S/(S+N) #Sigma^{-}#pi^{+}");
    Cospicm_purity_Sm->Divide(Cospicm_pi_Sm_S,Cospicm_pi_Sm_SN,1.,1.,"b");
    Cospicm_purity_Sm->Draw();

  }


  TCanvas *c1 = new TCanvas("c1","c1",1000,800);
  TH2F* MMnpi_IMnpip = (TH2F*)f->Get("MMnpi_IMnpip");
  MMnpi_IMnpip->RebinX(5);
  MMnpi_IMnpip->GetXaxis()->SetRangeUser(1,1.5);
  MMnpi_IMnpip->GetYaxis()->SetRangeUser(-0.6,0.6);
  MMnpi_IMnpip->Draw("colz");
  
  TCanvas *c11 = new TCanvas("c11","c11",1000,800);
  TH2F* MM2npi_IMnpip = (TH2F*)f->Get("MM2npi_IMnpip");
  MM2npi_IMnpip->RebinX(5);
  MM2npi_IMnpip->GetXaxis()->SetRangeUser(1,1.5);
  MM2npi_IMnpip->GetYaxis()->SetRangeUser(-0.6,0.4);
  MM2npi_IMnpip->Draw("colz");
  
  TCanvas *c2 = new TCanvas("c2","c2",1000,800);
  TH2F* MMpi_IMnpip_pi = (TH2F*)f->Get("MMpi_IMnpip_pi");
  MMpi_IMnpip_pi->RebinX(5);
  MMpi_IMnpip_pi->GetXaxis()->SetRangeUser(1,1.5);
  MMpi_IMnpip_pi->GetYaxis()->SetRangeUser(1,1.7);
  MMpi_IMnpip_pi->Draw("colz");
  
  TCanvas *c2_2 = new TCanvas("c2_2","c2_2",1000,800);
  const int sigma_low = MMpi_IMnpip_pi->GetXaxis()->FindBin(1.17);
  const int sigma_hi = MMpi_IMnpip_pi->GetXaxis()->FindBin(1.20);
  TH1D* MMpi_sp = MMpi_IMnpip_pi->ProjectionY("MMpi_sp",sigma_low,sigma_hi);
  MMpi_sp->Draw("E");


  TCanvas *c3 = new TCanvas("c3","c3",1000,800);
  TH2F* MMn_IMnpip_pi = (TH2F*)f->Get("MMn_IMnpip_pi");
  MMn_IMnpip_pi->RebinX(5);
  MMn_IMnpip_pi->GetXaxis()->SetRangeUser(1,1.5);
  MMn_IMnpip_pi->GetYaxis()->SetRangeUser(0,0.7);
  MMn_IMnpip_pi->Draw("colz");

  TCanvas *c4 = new TCanvas("c4","c4",1000,800);
  TH2F*  Cospicm_IMnpip_pi = (TH2F*)f->Get("Cospicm_IMnpip_pi");
  Cospicm_IMnpip_pi->RebinX(5);
  Cospicm_IMnpip_pi->GetXaxis()->SetRangeUser(1,1.5);
  //Cospicm_IMnpip->GetYaxis()->SetRangeUser(0,0.7);
  Cospicm_IMnpip_pi->Draw("colz");

  TCanvas *c4_1 = new TCanvas("c4_1","c4_1",1000,800);
  TH1D* IMnpip_pi = (TH1D*)Cospicm_IMnpip_pi->ProjectionX("IMnpip_pi");
  IMnpip_pi->Draw("E");

  TCanvas *c5 = new TCanvas("c5","c5",1000,800);
  TH2F* MMnpi_IMnpim = (TH2F*)f->Get("MMnpi_IMnpim");
  MMnpi_IMnpim->RebinX(5);
  MMnpi_IMnpim->GetXaxis()->SetRangeUser(1,1.5);
  MMnpi_IMnpim->GetYaxis()->SetRangeUser(-0.6,0.6);
  MMnpi_IMnpim->Draw("colz");
  
  TCanvas *c55 = new TCanvas("c55","c55",1000,800);
  TH2F* MM2npi_IMnpim = (TH2F*)f->Get("MM2npi_IMnpim");
  MM2npi_IMnpim->RebinX(5);
  MM2npi_IMnpim->GetXaxis()->SetRangeUser(1,1.5);
  MM2npi_IMnpim->GetYaxis()->SetRangeUser(-0.6,0.4);
  MM2npi_IMnpim->Draw("colz");

  TCanvas *c6 = new TCanvas("c6","c6",1000,800);
  TH2F* MMpi_IMnpim = (TH2F*)f->Get("MMpi_IMnpim");
  MMpi_IMnpim->RebinX(5);
  MMpi_IMnpim->GetXaxis()->SetRangeUser(1,1.5);
  MMpi_IMnpim->GetYaxis()->SetRangeUser(1,1.7);
  MMpi_IMnpim->Draw("colz");

  TCanvas *c7 = new TCanvas("c7","c7",1000,800);
  TH2F* MMn_IMnpim = (TH2F*)f->Get("MMn_IMnpim");
  MMn_IMnpim->RebinX(5);
  MMn_IMnpim->GetXaxis()->SetRangeUser(1,1.5);
  MMn_IMnpim->GetYaxis()->SetRangeUser(0,0.7);
  MMn_IMnpim->Draw("colz");

  TCanvas *c8 = new TCanvas("c8","c8",1000,800);
  TH2F* Cospicm_IMnpim_pi = (TH2F*)f->Get("Cospicm_IMnpim_pi");
  Cospicm_IMnpim_pi->RebinX(5);
  Cospicm_IMnpim_pi->GetXaxis()->SetRangeUser(1,1.5);
  Cospicm_IMnpim_pi->Draw("colz");

  TCanvas *c8_1 = new TCanvas("c8_1","c8_1",1000,800);
  TH1D* IMnpim_pi = (TH1D*)Cospicm_IMnpim_pi->ProjectionX("IMnpim_pi");
  IMnpim_pi->Draw("E");
  
  const int ncosbin=6;
  int cosbin[7];
  cosbin[0] = Cospicm_IMnpip_pi->GetYaxis()->FindBin(1.0);
  cosbin[1] = Cospicm_IMnpip_pi->GetYaxis()->FindBin(0.9); 
  cosbin[2] = Cospicm_IMnpip_pi->GetYaxis()->FindBin(0.8); 
  cosbin[3] = Cospicm_IMnpip_pi->GetYaxis()->FindBin(0.7); 
  cosbin[4] = Cospicm_IMnpip_pi->GetYaxis()->FindBin(0.6); 
  cosbin[5] = Cospicm_IMnpip_pi->GetYaxis()->FindBin(0.5); 
  cosbin[6] = Cospicm_IMnpip_pi->GetYaxis()->FindBin(0.4); 
  
  TCanvas *ccos_Sp = new TCanvas("ccos_Sp","ccos_Sp",1000,800);
  ccos_Sp->Divide(3,2);
  TH1D* IMnpip_coscut[6];
  const int Splow = Cospicm_IMnpip_pi->GetXaxis()->FindBin(anacuts::Sigmap_MIN);
  const int Sphigh = Cospicm_IMnpip_pi->GetXaxis()->FindBin(anacuts::Sigmap_MAX);
  for(int icosbin=0;icosbin<ncosbin;icosbin++){
    IMnpip_coscut[icosbin] = (TH1D*)Cospicm_IMnpip_pi->ProjectionX(Form("IMnpip_coscut%d",icosbin),cosbin[icosbin+1],cosbin[icosbin]);
    ccos_Sp->cd(icosbin+1);
    //std::cout << IMnpip_coscut[icosbin]->Integral(Splow,Sphigh) << std::endl;
    IMnpip_coscut[icosbin]->GetXaxis()->SetRangeUser(1.1,1.3);
    IMnpip_coscut[icosbin]->SetTitle(Form("cos%0.1f-cos%0.1f ", 1-icosbin*0.1, 0.9-icosbin*0.1));
    IMnpip_coscut[icosbin]->Draw("HE");
  }
  
  TCanvas *ccos_Sm = new TCanvas("ccos_Sm","ccos_Sm",1000,800);
  ccos_Sm->Divide(3,2);
  TH1D* IMnpim_coscut[6];
  const int Smlow = Cospicm_IMnpim_pi->GetXaxis()->FindBin(anacuts::Sigmam_MIN);
  const int Smhigh = Cospicm_IMnpim_pi->GetXaxis()->FindBin(anacuts::Sigmam_MAX);
  for(int icosbin=0;icosbin<ncosbin;icosbin++){
    IMnpim_coscut[icosbin] = (TH1D*)Cospicm_IMnpim_pi->ProjectionX(Form("IMnpim_coscut%d",icosbin),cosbin[icosbin+1],cosbin[icosbin]);
    ccos_Sm->cd(icosbin+1);
    //std::cout << IMnpim_coscut[icosbin]->Integral(Smlow,Smhigh) << std::endl;
    IMnpim_coscut[icosbin]->GetXaxis()->SetRangeUser(1.1,1.3);
    IMnpim_coscut[icosbin]->SetTitle(Form("cos%0.1f-cos%0.1f ", 1-icosbin*0.1, 0.9-icosbin*0.1));
    IMnpim_coscut[icosbin]->Draw("HE");
  }
  
  TCanvas *cCospicim_MM2npi_Sp = new TCanvas("cCospicim_MM2npi_Sp","cCospicim_MM2npi_Sp",1000,800);
  TH2F* Cospicm_MM2npi_Sp = (TH2F*)f->Get("Cospicm_MM2npi_Sp");
  Cospicm_MM2npi_Sp->Draw("colz");

  TCanvas *ccos_Sp_MM2 = new TCanvas("ccos_Sp_MM2","ccos_Sp_MM2",1000,800);
  ccos_Sp_MM2->Divide(3,2);
  TH1D* MM2npi_coscut_Sp[6];
  for(int icosbin=0;icosbin<ncosbin;icosbin++){
    MM2npi_coscut_Sp[icosbin] = (TH1D*)Cospicm_MM2npi_Sp->ProjectionX(Form("MM2npi_coscut%d_Sp",icosbin),cosbin[icosbin+1],cosbin[icosbin]);
    ccos_Sp_MM2->cd(icosbin+1);
    //std::cout << IMnpim_coscut[icosbin]->Integral(Smlow,Smhigh) << std::endl;
    MM2npi_coscut_Sp[icosbin]->Rebin(2);
    MM2npi_coscut_Sp[icosbin]->GetXaxis()->SetRangeUser(-0.4,0.4);
    MM2npi_coscut_Sp[icosbin]->SetTitle(Form("cos%0.1f-cos%0.1f ", 1-icosbin*0.1, 0.9-icosbin*0.1));
    MM2npi_coscut_Sp[icosbin]->Draw("HE");
  }
  
  TCanvas *cCospicim_MM2npi_Sm = new TCanvas("cCospicim_MM2npi_Sm","cCospicim_MM2npi_Sm",1000,800);
  TH2F* Cospicm_MM2npi_Sm = (TH2F*)f->Get("Cospicm_MM2npi_Sm");
  Cospicm_MM2npi_Sm->Draw("colz");

  TCanvas *ccos_Sm_MM2 = new TCanvas("ccos_Sm_MM2","ccos_Sm_MM2",1000,800);
  ccos_Sm_MM2->Divide(3,2);
  TH1D* MM2npi_coscut_Sm[6];
  for(int icosbin=0;icosbin<ncosbin;icosbin++){
    MM2npi_coscut_Sm[icosbin] = (TH1D*)Cospicm_MM2npi_Sm->ProjectionX(Form("MM2npi_coscut%d_Sm",icosbin),cosbin[icosbin+1],cosbin[icosbin]);
    ccos_Sm_MM2->cd(icosbin+1);
    //std::cout << IMnpim_coscut[icosbin]->Integral(Smlow,Smhigh) << std::endl;
    MM2npi_coscut_Sm[icosbin]->Rebin(2);
    MM2npi_coscut_Sm[icosbin]->GetXaxis()->SetRangeUser(-0.4,0.4);
    MM2npi_coscut_Sm[icosbin]->SetTitle(Form("cos%0.1f-cos%0.1f ", 1-icosbin*0.1, 0.9-icosbin*0.1));
    MM2npi_coscut_Sm[icosbin]->Draw("HE");
  }


  
  TCanvas *cCospicm_pi_Sp = new TCanvas("cCospicm_pi_Sp","cCospicm_pi_Sp",1000,800);
  TH1F* Cospicm_pi_Sp = (TH1F*)f->Get("Cospicm_pi_Sp");
  Cospicm_pi_Sp->RebinX(5);
  Cospicm_pi_Sp->Draw("E");

  TCanvas *cCospicm_pi_Sm = new TCanvas("cCospicm_pi_Sm","cCospicm_pi_Sm",1000,800);
  TH1F* Cospicm_pi_Sm = (TH1F*)f->Get("Cospicm_pi_Sm");
  Cospicm_pi_Sm->RebinX(5);
  Cospicm_pi_Sm->Draw("E");
  
  TCanvas *cCospicm_nmom_Sp = new TCanvas("cCospicm_nmom_Sp","cCospicm_nmom_Sp");
  TH2F* Cospicm_nmom_Sp = (TH2F*)f->Get("Cospicm_nmom_Sp");
  Cospicm_nmom_Sp->Draw("colz");
  
  TCanvas *cCospicm_nmom_Sm = new TCanvas("cCospicm_nmom_Sm","cCospicm_nmom_Sm");
  TH2F* Cospicm_nmom_Sm = (TH2F*)f->Get("Cospicm_nmom_Sm");
  Cospicm_nmom_Sm->Draw("colz");


  TCanvas *c = NULL;
  TSeqCollection *SCol = gROOT->GetListOfCanvases();
  int size = SCol->GetSize();
  TIter next(SCol);
  TString pdfname;
  pdfname= Form("H2_datav%d_dE%d_sys%d.pdf",version,dEcut,sysud);
  if(SimSpMode) pdfname= "H2data_sim_Sp.pdf";
  if(SimSmMode) pdfname= "H2data_sim_Sm.pdf";
   
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







}
