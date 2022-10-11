const bool RemoveNotEnough = true;
const double UncertCut = 0.20;
const int ncut=6;

void GetAccMapLpim()
{
   
  gStyle->SetPalette(1);
  gStyle->SetOptStat("ei");
  gStyle->SetStatX(0.9);
  gStyle->SetStatY(0.9);
  TH1::SetDefaultSumw2(true);
  TFile *fLpim=NULL;
  TFile *fgen=NULL;
  const int version = 29;

  fLpim = TFile::Open(Form("../simpost/simIMLpim_ppimpim_v%d_out.root",version));
  fgen = TFile::Open(Form("../simpost/simIMLpim_v%d.root",version));
  TFile *fkin = TFile::Open("../simpost/NumericalRootFinderLPim.root","READ");
//  fLpim = TFile::Open("simIMLpim_ppimpim_pS0pim_v1_out.root");
//  fgen = TFile::Open("simIMLpim_pS0pim_v1.root");
  
  TH2F* q_IMLPim_gen=NULL;
  q_IMLPim_gen = (TH2F*) fgen->Get("React_q_IMLPim");
//  q_IMLPim_gen = (TH2F*) fgen->Get("React_q_IMS0Pim");
  q_IMLPim_gen->SetTitle("generated evt. ");
  q_IMLPim_gen->SetXTitle("true IM(#pi^{-}#Lambda) [GeV/c^{2}]");
  q_IMLPim_gen->SetYTitle("true Mom. Transfer [GeV/c]");
  q_IMLPim_gen->GetXaxis()->CenterTitle();
  q_IMLPim_gen->GetYaxis()->CenterTitle();
  q_IMLPim_gen->Print("base");
  q_IMLPim_gen->RebinX(15);
  q_IMLPim_gen->RebinY(10);
  q_IMLPim_gen->Print("base");
   

  TH1F* BLAnaPassed = (TH1F*)fgen->Get("BLAnaPassed");
  const double SimBeamSurvivalOK = BLAnaPassed->GetBinContent(2);//passed
  const double SimBeamSurvivalFail = BLAnaPassed->GetBinContent(1);//not passed
  const double SimBeamSurvivalRate = SimBeamSurvivalOK / (SimBeamSurvivalOK+SimBeamSurvivalFail);
  std::cout << "L " << __LINE__ << " SimBeamSurvivalRate " << SimBeamSurvivalRate << std::endl;
   
  //get sum of p1 and p2
  TH2F* q_IMppipi_p_wL_sum[ncut];
  for(int icut=0;icut<ncut;icut++){
    q_IMppipi_p_wL_sum[icut]= (TH2F*)fLpim->Get(Form("q_IMppipi_p_wL_sum_nocombi%d",icut));
    q_IMppipi_p_wL_sum[icut]->GetYaxis()->SetRangeUser(0,1.5);
    q_IMppipi_p_wL_sum[icut]->SetTitle("reco. evt.");
    q_IMppipi_p_wL_sum[icut]->Scale(1./SimBeamSurvivalRate);
    //q_IMppipi_p_wL_sum[icut]->RebinX(3);
    //q_IMppipi_p_wL_sum[icut]->RebinY(3);
  }

  TH2F* q_IMppipi_p_wL_acc[ncut];
  for(int icut=0;icut<ncut;icut++){
    q_IMppipi_p_wL_acc[icut] = (TH2F*)q_IMppipi_p_wL_sum[icut]->Clone(Form("q_IMppipi_p_wL_acc%d",icut));
    q_IMppipi_p_wL_acc[icut]->SetTitle(Form("q_IMppipi_p_wL acc. %d",icut));
    q_IMppipi_p_wL_acc[icut]->Print("base");
    q_IMppipi_p_wL_acc[icut]->Divide(q_IMppipi_p_wL_acc[icut],q_IMLPim_gen,1.0,1.0,"b");
  }

  TH2F* q_IMppipi_p_wL_accerr[ncut];
  for(int icut=0;icut<ncut;icut++){
    q_IMppipi_p_wL_accerr[icut] = (TH2F*)q_IMppipi_p_wL_acc[icut]->Clone(Form("q_IMppipi_p_wL_accerr%d",icut));
    q_IMppipi_p_wL_accerr[icut]->Reset();
    q_IMppipi_p_wL_accerr[icut]->SetTitle(Form("q_IMppipi_p_wL precision %d",icut));
    for(int ix=0;ix<q_IMppipi_p_wL_accerr[icut]->GetNbinsX();ix++){
      for(int iy=0;iy<q_IMppipi_p_wL_accerr[icut]->GetNbinsY();iy++){
        double cont = q_IMppipi_p_wL_acc[icut]->GetBinContent(ix,iy);
        double err = q_IMppipi_p_wL_acc[icut]->GetBinError(ix,iy);
        if(cont!=0){
          q_IMppipi_p_wL_accerr[icut]->SetBinContent(ix,iy,err/cont);  
        }
      }
    }
  }//icut
  
  //get no p2 hists
  TH2F* q_IMppipi_p_wL_sum_nop2[ncut];//0:default, 1 half low, 2 half high, 3 sigma0 region, 4 wide range
  for(int icut=0;icut<ncut;icut++){
    q_IMppipi_p_wL_sum_nop2[icut]= (TH2F*)fLpim->Get(Form("q_IMppipi_p_wL_sum_nocombi_nop2%d",icut));
    q_IMppipi_p_wL_sum_nop2[icut]->GetYaxis()->SetRangeUser(0,1.5);
    q_IMppipi_p_wL_sum_nop2[icut]->SetTitle("reco. evt.");
    q_IMppipi_p_wL_sum_nop2[icut]->Scale(1./SimBeamSurvivalRate);
  }//icut 

  TH2F* q_IMppipi_p_wL_nop2_acc[ncut];
  for(int icut=0;icut<ncut;icut++){
    q_IMppipi_p_wL_nop2_acc[icut] = (TH2F*)q_IMppipi_p_wL_sum_nop2[icut]->Clone(Form("q_IMppipi_p_wL_nop2_acc%d",icut));
    q_IMppipi_p_wL_nop2_acc[icut]->SetTitle(Form("q_IMppipi_p_wL nop2 acc. %d",icut));
    q_IMppipi_p_wL_nop2_acc[icut]->Print("base");
    q_IMppipi_p_wL_nop2_acc[icut]->Divide(q_IMppipi_p_wL_nop2_acc[icut],q_IMLPim_gen,1.0,1.0,"b");
  }

  TH2F* q_IMppipi_p_wL_nop2_accerr[ncut];
  for(int icut=0;icut<ncut;icut++){
    q_IMppipi_p_wL_nop2_accerr[icut] = (TH2F*)q_IMppipi_p_wL_nop2_acc[icut]->Clone(Form("q_IMppipi_p_wL_nop2_accerr%d",icut));
    q_IMppipi_p_wL_nop2_accerr[icut]->Reset();
    q_IMppipi_p_wL_nop2_accerr[icut]->SetTitle(Form("q_IMppipi_p_wL_nop2 precision %d",icut));
    for(int ix=0;ix<q_IMppipi_p_wL_nop2_accerr[icut]->GetNbinsX();ix++){
      for(int iy=0;iy<q_IMppipi_p_wL_nop2_accerr[icut]->GetNbinsY();iy++){
        double cont = q_IMppipi_p_wL_nop2_acc[icut]->GetBinContent(ix,iy);
        double err = q_IMppipi_p_wL_nop2_acc[icut]->GetBinError(ix,iy);
        if(cont!=0){
          q_IMppipi_p_wL_nop2_accerr[icut]->SetBinContent(ix,iy,err/cont);  
        }
      }
    }
  }//icut


  //get no p2 hists_mc w/ true val.
  TH2F* q_IMppipi_p_wL_sum_nop2_mc[ncut];//0:default, 1 half low, 2 half high, 3 sigma0 region, 4 wide range
  for(int icut=0;icut<ncut;icut++){
    q_IMppipi_p_wL_sum_nop2_mc[icut]= (TH2F*)fLpim->Get(Form("q_IMppipi_p_wL_sum_nocombi_nop2_mc%d",icut));
    q_IMppipi_p_wL_sum_nop2_mc[icut]->GetYaxis()->SetRangeUser(0,1.5);
    q_IMppipi_p_wL_sum_nop2_mc[icut]->SetTitle("reco. evt.");
    q_IMppipi_p_wL_sum_nop2_mc[icut]->Scale(1./SimBeamSurvivalRate);
  }//icut 

  TH2F* q_IMppipi_p_wL_nop2_mc_acc[ncut];
  for(int icut=0;icut<ncut;icut++){
    q_IMppipi_p_wL_nop2_mc_acc[icut] = (TH2F*)q_IMppipi_p_wL_sum_nop2_mc[icut]->Clone(Form("q_IMppipi_p_wL_nop2_mc_acc%d",icut));
    q_IMppipi_p_wL_nop2_mc_acc[icut]->SetTitle(Form("q_IMppipi_p_wL nop2 true val. acc. %d",icut));
    q_IMppipi_p_wL_nop2_mc_acc[icut]->Print("base");
    q_IMppipi_p_wL_nop2_mc_acc[icut]->Divide(q_IMppipi_p_wL_nop2_mc_acc[icut],q_IMLPim_gen,1.0,1.0,"b");
  }

  TH2F* q_IMppipi_p_wL_nop2_mc_accerr[ncut];
  for(int icut=0;icut<ncut;icut++){
    q_IMppipi_p_wL_nop2_mc_accerr[icut] = (TH2F*)q_IMppipi_p_wL_nop2_mc_acc[icut]->Clone(Form("q_IMppipi_p_wL_nop2_mc_accerr%d",icut));
    q_IMppipi_p_wL_nop2_mc_accerr[icut]->Reset();
    q_IMppipi_p_wL_nop2_mc_accerr[icut]->SetTitle(Form("q_IMppipi_p_wL_nop2_mc precision %d",icut));
    for(int ix=0;ix<q_IMppipi_p_wL_nop2_mc_accerr[icut]->GetNbinsX();ix++){
      for(int iy=0;iy<q_IMppipi_p_wL_nop2_mc_accerr[icut]->GetNbinsY();iy++){
        double cont = q_IMppipi_p_wL_nop2_mc_acc[icut]->GetBinContent(ix,iy);
        double err = q_IMppipi_p_wL_nop2_mc_acc[icut]->GetBinError(ix,iy);
        if(cont!=0){
          q_IMppipi_p_wL_nop2_mc_accerr[icut]->SetBinContent(ix,iy,err/cont);  
        }
      }
    }
  }//icut

  //get w p2 hists
  TH2F* q_IMppipi_p_wL_sum_wp2[ncut];//0:default, 1 half low, 2 half high, 3 sigma0 region, 4 wide range
  for(int icut=0;icut<ncut;icut++){
    q_IMppipi_p_wL_sum_wp2[icut]= (TH2F*)fLpim->Get(Form("q_IMppipi_p_wL_sum_nocombi_wp2%d",icut));
    q_IMppipi_p_wL_sum_wp2[icut]->GetYaxis()->SetRangeUser(0,1.5);
    q_IMppipi_p_wL_sum_wp2[icut]->SetTitle("reco. evt.");
    q_IMppipi_p_wL_sum_wp2[icut]->Scale(1./SimBeamSurvivalRate);
  }//icut 

  TH2F* q_IMppipi_p_wL_wp2_acc[ncut];
  for(int icut=0;icut<ncut;icut++){
    q_IMppipi_p_wL_wp2_acc[icut] = (TH2F*)q_IMppipi_p_wL_sum_wp2[icut]->Clone(Form("q_IMppipi_p_wL_wp2_acc%d",icut));
    q_IMppipi_p_wL_wp2_acc[icut]->SetTitle(Form("q_IMppipi_p_wL wp2 acc. %d",icut));
    q_IMppipi_p_wL_wp2_acc[icut]->Print("base");
    q_IMppipi_p_wL_wp2_acc[icut]->Divide(q_IMppipi_p_wL_wp2_acc[icut],q_IMLPim_gen,1.0,1.0,"b");
  }

  TH2F* q_IMppipi_p_wL_wp2_accerr[ncut];
  for(int icut=0;icut<ncut;icut++){
    q_IMppipi_p_wL_wp2_accerr[icut] = (TH2F*)q_IMppipi_p_wL_wp2_acc[icut]->Clone(Form("q_IMppipi_p_wL_wp2_accerr%d",icut));
    q_IMppipi_p_wL_wp2_accerr[icut]->Reset();
    q_IMppipi_p_wL_wp2_accerr[icut]->SetTitle(Form("q_IMppipi_p_wL_wp2 precision %d",icut));
    for(int ix=0;ix<q_IMppipi_p_wL_wp2_accerr[icut]->GetNbinsX();ix++){
      for(int iy=0;iy<q_IMppipi_p_wL_wp2_accerr[icut]->GetNbinsY();iy++){
        double cont = q_IMppipi_p_wL_wp2_acc[icut]->GetBinContent(ix,iy);
        double err = q_IMppipi_p_wL_wp2_acc[icut]->GetBinError(ix,iy);
        if(cont!=0){
          q_IMppipi_p_wL_wp2_accerr[icut]->SetBinContent(ix,iy,err/cont);  
        }
      }
    }
  }//icut
 
  //get w p2 hists w/ true val.
  TH2F* q_IMppipi_p_wL_sum_wp2_mc[ncut];//0:default, 1 half low, 2 half high, 3 sigma0 region, 4 wide range
  for(int icut=0;icut<ncut;icut++){
    q_IMppipi_p_wL_sum_wp2_mc[icut]= (TH2F*)fLpim->Get(Form("q_IMppipi_p_wL_sum_nocombi_wp2_mc%d",icut));
    q_IMppipi_p_wL_sum_wp2_mc[icut]->GetYaxis()->SetRangeUser(0,1.5);
    q_IMppipi_p_wL_sum_wp2_mc[icut]->SetTitle("reco. evt.");
    q_IMppipi_p_wL_sum_wp2_mc[icut]->Scale(1./SimBeamSurvivalRate);
  }//icut 

  TH2F* q_IMppipi_p_wL_wp2_mc_acc[ncut];
  for(int icut=0;icut<ncut;icut++){
    q_IMppipi_p_wL_wp2_mc_acc[icut] = (TH2F*)q_IMppipi_p_wL_sum_wp2_mc[icut]->Clone(Form("q_IMppipi_p_wL_wp2_mc_acc%d",icut));
    q_IMppipi_p_wL_wp2_mc_acc[icut]->SetTitle(Form("q_IMppipi_p_wL wp2 true val. acc. %d",icut));
    q_IMppipi_p_wL_wp2_mc_acc[icut]->Print("base");
    q_IMppipi_p_wL_wp2_mc_acc[icut]->Divide(q_IMppipi_p_wL_wp2_mc_acc[icut],q_IMLPim_gen,1.0,1.0,"b");
  }

  TH2F* q_IMppipi_p_wL_wp2_mc_accerr[ncut];
  for(int icut=0;icut<ncut;icut++){
    q_IMppipi_p_wL_wp2_mc_accerr[icut] = (TH2F*)q_IMppipi_p_wL_wp2_mc_acc[icut]->Clone(Form("q_IMppipi_p_wL_wp2_mc_accerr%d",icut));
    q_IMppipi_p_wL_wp2_mc_accerr[icut]->Reset();
    q_IMppipi_p_wL_wp2_mc_accerr[icut]->SetTitle(Form("q_IMppipi_p_wL_wp2_mc precision %d",icut));
    for(int ix=0;ix<q_IMppipi_p_wL_wp2_mc_accerr[icut]->GetNbinsX();ix++){
      for(int iy=0;iy<q_IMppipi_p_wL_wp2_mc_accerr[icut]->GetNbinsY();iy++){
        double cont = q_IMppipi_p_wL_wp2_mc_acc[icut]->GetBinContent(ix,iy);
        double err = q_IMppipi_p_wL_wp2_mc_acc[icut]->GetBinError(ix,iy);
        if(cont!=0){
          q_IMppipi_p_wL_wp2_mc_accerr[icut]->SetBinContent(ix,iy,err/cont);  
        }
      }
    }
  }//icut


  TCanvas *cLpim[ncut];
  TGraph *gth = (TGraph*)fkin->Get("th");
  TGraph *gr_0 = (TGraph*)fkin->Get("gr_0");
  TGraph *gr_100 = (TGraph*)fkin->Get("gr_100");
  TGraph *gr_65 = (TGraph*)fkin->Get("gr_65");
  
  for(int icut=0;icut<ncut;icut++){
    cLpim[icut] = new TCanvas(Form("cLpim%d",icut),Form("cLpim%d",icut),2000,1200);
    cLpim[icut]->Divide(2,2);
    cLpim[icut]->cd(1);
    q_IMLPim_gen->Draw("colz");
 
 
    //gth->Draw("pc");
    //gr_0->Draw("pc");
    //gr_65->Draw("pc");
    //gr_100->Draw("pc");
    cLpim[icut]->cd(2);
    q_IMppipi_p_wL_sum[icut]->Draw("colz");
    //gth->Draw("pc");
    //gr_0->Draw("pc");
    //gr_100->Draw("pc");
    //gr_65->Draw("pc");
    cLpim[icut]->cd(3);
    //q_IMppipi_p_wL_acc[icut]->SetMaximum(0.055);
    q_IMppipi_p_wL_acc[icut]->SetMaximum(0.075);
    q_IMppipi_p_wL_acc[icut]->Draw("colz");
    //gth->Draw("pc");
    //gr_0->Draw("pc");
    //gr_100->Draw("pc");
    //gr_65->Draw("pc");
    cLpim[icut]->cd(4);
    q_IMppipi_p_wL_accerr[icut]->SetMaximum(0.5);
    q_IMppipi_p_wL_accerr[icut]->Draw("colz");
  }
  
  TCanvas *cLpim_nop2[ncut];
  for(int icut=0;icut<ncut;icut++){
    cLpim_nop2[icut] = new TCanvas(Form("cLpim_nop2%d",icut),Form("cLpim_nop2%d",icut),2000,1200);
    cLpim_nop2[icut]->Divide(2,2);
    cLpim_nop2[icut]->cd(1);
    q_IMLPim_gen->Draw("colz");
 
 
    //gth->Draw("pc");
    //gr_0->Draw("pc");
    //gr_65->Draw("pc");
    //gr_100->Draw("pc");
    cLpim_nop2[icut]->cd(2);
    q_IMppipi_p_wL_sum_nop2[icut]->Draw("colz");
    //gth->Draw("pc");
    //gr_0->Draw("pc");
    //gr_100->Draw("pc");
    //gr_65->Draw("pc");
    cLpim_nop2[icut]->cd(3);
    //q_IMppipi_p_wL_acc[icut]->SetMaximum(0.055);
    q_IMppipi_p_wL_nop2_acc[icut]->SetMaximum(0.075);
    q_IMppipi_p_wL_nop2_acc[icut]->Draw("colz");
    //gth->Draw("pc");
    //gr_0->Draw("pc");
    //gr_100->Draw("pc");
    //gr_65->Draw("pc");
    cLpim_nop2[icut]->cd(4);
    q_IMppipi_p_wL_nop2_accerr[icut]->SetMaximum(0.5);
    q_IMppipi_p_wL_nop2_accerr[icut]->Draw("colz");
  }
  
  TCanvas *cLpim_nop2_mc[ncut];
  for(int icut=0;icut<ncut;icut++){
    cLpim_nop2_mc[icut] = new TCanvas(Form("cLpim_nop2_mc%d",icut),Form("cLpim_nop2_mc%d",icut),2000,1200);
    cLpim_nop2_mc[icut]->Divide(2,2);
    cLpim_nop2_mc[icut]->cd(1);
    q_IMLPim_gen->Draw("colz");
 
 
    //gth->Draw("pc");
    //gr_0->Draw("pc");
    //gr_65->Draw("pc");
    //gr_100->Draw("pc");
    cLpim_nop2_mc[icut]->cd(2);
    q_IMppipi_p_wL_sum_nop2_mc[icut]->Draw("colz");
    //gth->Draw("pc");
    //gr_0->Draw("pc");
    //gr_100->Draw("pc");
    //gr_65->Draw("pc");
    cLpim_nop2_mc[icut]->cd(3);
    //q_IMppipi_p_wL_acc[icut]->SetMaximum(0.055);
    q_IMppipi_p_wL_nop2_mc_acc[icut]->SetMaximum(0.075);
    q_IMppipi_p_wL_nop2_mc_acc[icut]->Draw("colz");
    //gth->Draw("pc");
    //gr_0->Draw("pc");
    //gr_100->Draw("pc");
    //gr_65->Draw("pc");
    cLpim_nop2_mc[icut]->cd(4);
    q_IMppipi_p_wL_nop2_mc_accerr[icut]->SetMaximum(0.5);
    q_IMppipi_p_wL_nop2_mc_accerr[icut]->Draw("colz");
  }
  
  TCanvas *cLpim_wp2[ncut];
  for(int icut=0;icut<ncut;icut++){
    cLpim_wp2[icut] = new TCanvas(Form("cLpim_wp2%d",icut),Form("cLpim_wp2%d",icut),2000,1200);
    cLpim_wp2[icut]->Divide(2,2);
    cLpim_wp2[icut]->cd(1);
    q_IMLPim_gen->Draw("colz");
 
 
    //gth->Draw("pc");
    //gr_0->Draw("pc");
    //gr_65->Draw("pc");
    //gr_100->Draw("pc");
    cLpim_wp2[icut]->cd(2);
    q_IMppipi_p_wL_sum_wp2[icut]->Draw("colz");
    //gth->Draw("pc");
    //gr_0->Draw("pc");
    //gr_100->Draw("pc");
    //gr_65->Draw("pc");
    cLpim_wp2[icut]->cd(3);
    //q_IMppipi_p_wL_acc[icut]->SetMaximum(0.055);
    q_IMppipi_p_wL_wp2_acc[icut]->SetMaximum(0.075);
    q_IMppipi_p_wL_wp2_acc[icut]->Draw("colz");
    //gth->Draw("pc");
    //gr_0->Draw("pc");
    //gr_100->Draw("pc");
    //gr_65->Draw("pc");
    cLpim_wp2[icut]->cd(4);
    q_IMppipi_p_wL_wp2_accerr[icut]->SetMaximum(0.5);
    q_IMppipi_p_wL_wp2_accerr[icut]->Draw("colz");
  }
  

  TCanvas *cLpim_wp2_mc[ncut];
  for(int icut=0;icut<ncut;icut++){
    cLpim_wp2_mc[icut] = new TCanvas(Form("cLpim_wp2_mc%d",icut),Form("cLpim_wp2_mc%d",icut),2000,1200);
    cLpim_wp2_mc[icut]->Divide(2,2);
    cLpim_wp2_mc[icut]->cd(1);
    q_IMLPim_gen->Draw("colz");
 
 
    //gth->Draw("pc");
    //gr_0->Draw("pc");
    //gr_65->Draw("pc");
    //gr_100->Draw("pc");
    cLpim_wp2_mc[icut]->cd(2);
    q_IMppipi_p_wL_sum_wp2_mc[icut]->Draw("colz");
    //gth->Draw("pc");
    //gr_0->Draw("pc");
    //gr_100->Draw("pc");
    //gr_65->Draw("pc");
    cLpim_wp2_mc[icut]->cd(3);
    //q_IMppipi_p_wL_acc[icut]->SetMaximum(0.055);
    q_IMppipi_p_wL_wp2_mc_acc[icut]->SetMaximum(0.075);
    q_IMppipi_p_wL_wp2_mc_acc[icut]->Draw("colz");
    //gth->Draw("pc");
    //gr_0->Draw("pc");
    //gr_100->Draw("pc");
    //gr_65->Draw("pc");
    cLpim_wp2_mc[icut]->cd(4);
    q_IMppipi_p_wL_wp2_mc_accerr[icut]->SetMaximum(0.5);
    q_IMppipi_p_wL_wp2_mc_accerr[icut]->Draw("colz");
  }
  
  /*
  //test for reco.events/acc map 
  TCanvas *cLpim_gentest = new TCanvas("cLpim_gentest","cLpim_gentest",2000,800);
  cLpim_gentest->Divide(2,1);
  cLpim_gentest->cd(1);
  TH2F* q_IMppipi_p_wL_gentest = (TH2F*)q_IMppipi_p_wL_sum->Clone("q_IMppipi_p_wL_gentest");
  q_IMppipi_p_wL_gentest->Divide(q_IMppipi_p_wL_acc);
  q_IMppipi_p_wL_gentest->SetTitle("TEST reco/acc map");
  q_IMppipi_p_wL_gentest->Draw("colz");
  cLpim_gentest->cd(2);
  q_IMLPim_gen->Draw("colz");

  TH1F* htest = new TH1F("htest","htest",8000,-6000,2000);
  for(int ix=0;ix<q_IMppipi_p_wL_gentest->GetNbinsX();ix++){
    for(int iy=0;iy<q_IMppipi_p_wL_gentest->GetNbinsY();iy++){
      double cont = q_IMppipi_p_wL_gentest->GetBinContent(ix,iy);
      double gen = q_IMLPim_gen->GetBinContent(ix,iy);
      htest->Fill(cont-gen);
    }
  }
  TCanvas *ctest  = new TCanvas("ctest","ctest",1000,800);
  htest->Draw("HE");

  TCanvas *ccleanup[5];
  TH2F* q_IMppipi_p_wL_gentest_clean[5];

  for(int icut=0;icut<5;icut++){
    ccleanup[icut] = new TCanvas(Form("ccleanup%d",icut),Form("ccleanup%d",icut),1000,800);
    q_IMppipi_p_wL_gentest_clean[icut] = (TH2F*)q_IMppipi_p_wL_gentest[icut]->Clone(Form("q_IMppipi_p_wL_gentest_clean%d",icut));
    q_IMppipi_p_wL_acc_clean[icut] = (TH2F*)q_IMppipi_p_wL_acc[icut]->Clone(Form("q_IMppipi_p_wL_acc_clean%d",icut));
    q_IMppipi_p_wL_acc_clean[icut]->SetName(Form("q_IMppipi_p_wL_acc_clean%d",icut));
    q_IMppipi_p_wL_acc_clean[icut]->SetTitle(Form("q_IMppipi_p_wL_acc_clean%d",icut));
    for(int ix=0;ix<q_IMppipi_p_wL_gentest->GetNbinsX();ix++){
      for(int iy=0;iy<q_IMppipi_p_wL_gentest->GetNbinsY();iy++){
        double cont = q_IMppipi_p_wL_gentest->GetBinContent(ix,iy);
        double gen = q_IMLPim_gen->GetBinContent(ix,iy);
        if(cont -gen < -100 ){
          q_IMppipi_p_wL_gentest_clean->SetBinContent(ix,iy,0);
          q_IMppipi_p_wL_gentest_clean->SetBinError(ix,iy,0);
          q_IMppipi_p_wL_acc_clean->SetBinContent(ix,iy,0);
          q_IMppipi_p_wL_acc_clean->SetBinError(ix,iy,0);
        }
      }
    }
  }
  
  */
  
  TH2F* q_IMppipi_p_wL_acc_clean[ncut];     
  for(int icut=0;icut<ncut;icut++){
    q_IMppipi_p_wL_acc_clean[icut] = (TH2F*)q_IMppipi_p_wL_acc[icut]->Clone(Form("q_IMppipi_p_wL_acc_clean%d",icut));
    q_IMppipi_p_wL_acc_clean[icut]->SetName(Form("q_IMppipi_p_wL_acc_clean%d",icut));
    q_IMppipi_p_wL_acc_clean[icut]->SetTitle(Form("q_IMppipi_p_wL_acc_clean%d",icut));
    for(int ix=0;ix<q_IMppipi_p_wL_accerr[icut]->GetNbinsX();ix++){
      for(int iy=0;iy<q_IMppipi_p_wL_accerr[icut]->GetNbinsY();iy++){
        double err = q_IMppipi_p_wL_accerr[icut]->GetBinContent(ix,iy);
        if( RemoveNotEnough && (err>UncertCut) 
            //  || q_IMppipi_p_wL_acc->GetBinContent(ix,iy)>0.053
            //  || q_IMppipi_p_wL_acc->GetXaxis()->GetBinCenter(ix)>1.895 
            || q_IMppipi_p_wL_acc[icut]->GetXaxis()->GetBinCenter(ix)<1.260) {
          q_IMppipi_p_wL_acc_clean[icut]->SetBinContent(ix,iy,0);
          q_IMppipi_p_wL_acc_clean[icut]->SetBinError(ix,iy,0);
        }
      }
    }
  }

  TCanvas *cacccleanup[ncut];
  for(int icut=0;icut<ncut;icut++){
    cacccleanup[icut]= new TCanvas(Form("cacccleanup%d",icut),Form("cacccleanup%d",icut),1000,800);
    q_IMppipi_p_wL_acc_clean[icut]->SetMaximum(q_IMppipi_p_wL_acc[icut]->GetMaximum());
    q_IMppipi_p_wL_acc_clean[icut]->Draw("colz");
    gr_0->Draw("pc");
    gr_100->Draw("pc");
    gr_65->Draw("pc");
    gth->Draw("pc");
  }


  TH2F* q_IMppipi_p_wL_nop2_acc_clean[ncut];     
  TH2F* q_IMppipi_p_wL_nop2_mc_acc_clean[ncut];     
  for(int icut=0;icut<ncut;icut++){
    q_IMppipi_p_wL_nop2_acc_clean[icut] = (TH2F*)q_IMppipi_p_wL_nop2_acc[icut]->Clone(Form("q_IMppipi_p_wL_nop2_acc_clean%d",icut));
    q_IMppipi_p_wL_nop2_acc_clean[icut]->SetName(Form("q_IMppipi_p_wL_nop2_acc_clean%d",icut));
    q_IMppipi_p_wL_nop2_acc_clean[icut]->SetTitle(Form("q_IMppipi_p_wL_nop2_acc_clean%d",icut));
    q_IMppipi_p_wL_nop2_mc_acc_clean[icut] = (TH2F*)q_IMppipi_p_wL_nop2_mc_acc[icut]->Clone(Form("q_IMppipi_p_wL_nop2_mc_acc_clean%d",icut));
    q_IMppipi_p_wL_nop2_mc_acc_clean[icut]->SetName(Form("q_IMppipi_p_wL_nop2_mc_acc_clean%d",icut));
    q_IMppipi_p_wL_nop2_mc_acc_clean[icut]->SetTitle(Form("q_IMppipi_p_wL_nop2_mc_acc_clean%d",icut));
    for(int ix=0;ix<q_IMppipi_p_wL_nop2_accerr[icut]->GetNbinsX();ix++){
      for(int iy=0;iy<q_IMppipi_p_wL_nop2_accerr[icut]->GetNbinsY();iy++){
        double cont = q_IMppipi_p_wL_nop2_acc[icut]->GetBinContent(ix,iy);
        double err = q_IMppipi_p_wL_nop2_accerr[icut]->GetBinContent(ix,iy);
        if( RemoveNotEnough && (err>UncertCut) 
            //  || q_IMppipi_p_wL_acc->GetBinContent(ix,iy)>0.053
            //  || q_IMppipi_p_wL_acc->GetXaxis()->GetBinCenter(ix)>1.895 
            || q_IMppipi_p_wL_nop2_acc[icut]->GetXaxis()->GetBinCenter(ix)<1.260
            || cont == 0.0
            ) {
          q_IMppipi_p_wL_nop2_acc_clean[icut]->SetBinContent(ix,iy,0);
          q_IMppipi_p_wL_nop2_acc_clean[icut]->SetBinError(ix,iy,0);
          q_IMppipi_p_wL_nop2_mc_acc_clean[icut]->SetBinContent(ix,iy,0);
          q_IMppipi_p_wL_nop2_mc_acc_clean[icut]->SetBinError(ix,iy,0);
        }
      }
    }
  }

  TCanvas *cacccleanup_nop2[ncut];
  for(int icut=0;icut<ncut;icut++){
    cacccleanup_nop2[icut]= new TCanvas(Form("cacccleanup_nop2%d",icut),Form("cacccleanup_nop2%d",icut),1000,800);
    q_IMppipi_p_wL_nop2_acc_clean[icut]->SetMaximum(q_IMppipi_p_wL_nop2_acc[icut]->GetMaximum());
    q_IMppipi_p_wL_nop2_acc_clean[icut]->Draw("colz");
    gr_0->Draw("pc");
    gr_100->Draw("pc");
    gr_65->Draw("pc");
    gth->Draw("pc");
  }

  TH2F* q_IMppipi_p_wL_wp2_acc_clean[ncut];     
  TH2F* q_IMppipi_p_wL_wp2_mc_acc_clean[ncut];     
  for(int icut=0;icut<ncut;icut++){
    q_IMppipi_p_wL_wp2_acc_clean[icut] = (TH2F*)q_IMppipi_p_wL_wp2_acc[icut]->Clone(Form("q_IMppipi_p_wL_wp2_acc_clean%d",icut));
    q_IMppipi_p_wL_wp2_acc_clean[icut]->SetName(Form("q_IMppipi_p_wL_wp2_acc_clean%d",icut));
    q_IMppipi_p_wL_wp2_acc_clean[icut]->SetTitle(Form("q_IMppipi_p_wL_wp2_acc_clean%d",icut));
    q_IMppipi_p_wL_wp2_mc_acc_clean[icut] = (TH2F*)q_IMppipi_p_wL_wp2_mc_acc[icut]->Clone(Form("q_IMppipi_p_wL_wp2_mc_acc_clean%d",icut));
    q_IMppipi_p_wL_wp2_mc_acc_clean[icut]->SetName(Form("q_IMppipi_p_wL_wp2_mc_acc_clean%d",icut));
    q_IMppipi_p_wL_wp2_mc_acc_clean[icut]->SetTitle(Form("q_IMppipi_p_wL_wp2_mc_acc_clean%d",icut));
    for(int ix=0;ix<q_IMppipi_p_wL_wp2_accerr[icut]->GetNbinsX();ix++){
      for(int iy=0;iy<q_IMppipi_p_wL_wp2_accerr[icut]->GetNbinsY();iy++){
        double cont = q_IMppipi_p_wL_wp2_acc[icut]->GetBinContent(ix,iy);
        double err = q_IMppipi_p_wL_wp2_accerr[icut]->GetBinContent(ix,iy);
        if( RemoveNotEnough && (err>UncertCut) 
            //  || q_IMppipi_p_wL_acc->GetBinContent(ix,iy)>0.053
            //  || q_IMppipi_p_wL_acc->GetXaxis()->GetBinCenter(ix)>1.895 
            || q_IMppipi_p_wL_wp2_acc[icut]->GetXaxis()->GetBinCenter(ix)<1.260
            || cont == 0.0
            ) {
          q_IMppipi_p_wL_wp2_acc_clean[icut]->SetBinContent(ix,iy,0);
          q_IMppipi_p_wL_wp2_acc_clean[icut]->SetBinError(ix,iy,0);
          q_IMppipi_p_wL_wp2_mc_acc_clean[icut]->SetBinContent(ix,iy,0);
          q_IMppipi_p_wL_wp2_mc_acc_clean[icut]->SetBinError(ix,iy,0);
        }
      }
    }
  }

  TCanvas *cacccleanup_wp2[ncut];
  for(int icut=0;icut<ncut;icut++){
    cacccleanup_wp2[icut]= new TCanvas(Form("cacccleanup_wp2%d",icut),Form("cacccleanup_wp2%d",icut),1000,800);
    q_IMppipi_p_wL_wp2_acc_clean[icut]->SetMaximum(q_IMppipi_p_wL_wp2_acc[icut]->GetMaximum());
    q_IMppipi_p_wL_wp2_acc_clean[icut]->Draw("colz");
    gr_0->Draw("pc");
    gr_100->Draw("pc");
    gr_65->Draw("pc");
    gth->Draw("pc");
  }


  //////////////////////////////////////////////////////////
  //acceptance map for costheta_missing p (LAB frame) v.s. IM(pi-Lambda) 
  //////////////////////////////////////////////////////////
  TH2F* costhetap_IMLpim_gen=NULL;
  costhetap_IMLpim_gen = (TH2F*)fgen->Get("React_costhetap_IMLPim");
  costhetap_IMLpim_gen->SetTitle("generated evt. ");
  costhetap_IMLpim_gen->SetXTitle("true IM(#pi^{-}#Lambda) [GeV/c^{2}]");
  costhetap_IMLpim_gen->SetYTitle("true proton CosTheta");
  costhetap_IMLpim_gen->GetXaxis()->CenterTitle();
  costhetap_IMLpim_gen->GetYaxis()->CenterTitle();
  costhetap_IMLpim_gen->Print("base");
  costhetap_IMLpim_gen->RebinX(15);
  costhetap_IMLpim_gen->Print("base");
 
  TH2F* CosTheta_IMppipi_p_wL=NULL;
  CosTheta_IMppipi_p_wL = (TH2F*)fLpim->Get("CosTheta_IMppipi_p_wL_nocombi");
  CosTheta_IMppipi_p_wL ->SetTitle("reco. evt.");
  CosTheta_IMppipi_p_wL->Scale(1./SimBeamSurvivalRate);

  TH2F* CosTheta_IMppipi_p_wL_acc=NULL;
  CosTheta_IMppipi_p_wL_acc = (TH2F*)CosTheta_IMppipi_p_wL->Clone("CosTheta_IMppipi_p_wL_acc");
  CosTheta_IMppipi_p_wL_acc->SetTitle("CosTheta_IMppipi_p_wL_acc");
  CosTheta_IMppipi_p_wL_acc->Print("base");
  CosTheta_IMppipi_p_wL_acc->Divide(CosTheta_IMppipi_p_wL_acc,costhetap_IMLpim_gen,1.0,1.0,"b");


  TH2F* CosTheta_IMppipi_p_wL_accerr=NULL;
  CosTheta_IMppipi_p_wL_accerr = (TH2F*)CosTheta_IMppipi_p_wL_acc->Clone("CosTheta_IMppipi_p_wL_accerr");
  CosTheta_IMppipi_p_wL_accerr->Reset();
  CosTheta_IMppipi_p_wL_accerr->SetTitle("CosTheta_IMppipi_p_wL_acc precision");
  
  for(int ix=0;ix<CosTheta_IMppipi_p_wL_accerr->GetNbinsX();ix++){
    for(int iy=0;iy<CosTheta_IMppipi_p_wL_accerr->GetNbinsY();iy++){
      double cont = CosTheta_IMppipi_p_wL_acc->GetBinContent(ix,iy);
      double err = CosTheta_IMppipi_p_wL_acc->GetBinError(ix,iy);
      if(cont!=0){
        CosTheta_IMppipi_p_wL_accerr->SetBinContent(ix,iy,err/cont);  
      }
    }
  }
  
  TCanvas *cLpimCos;
  cLpimCos = new TCanvas(Form("cLpimCos"),Form("cLpimCos"),2000,1200);
  cLpimCos->Divide(2,2);
  cLpimCos->cd(1);
  costhetap_IMLpim_gen->GetYaxis()->SetRangeUser(0,1);
  costhetap_IMLpim_gen->Draw("colz");
  gPad->SetLogz();
  cLpimCos->cd(2);
  CosTheta_IMppipi_p_wL->GetYaxis()->SetRangeUser(0,1);
  CosTheta_IMppipi_p_wL->Draw("colz");
  gPad->SetLogz();
  cLpimCos->cd(3);
  //CosThetaIMppipi_p_wL_acc->SetMaximum(0.05);
  CosTheta_IMppipi_p_wL_acc->GetYaxis()->SetRangeUser(0,1);
  CosTheta_IMppipi_p_wL_acc->SetMaximum(0.05);
  CosTheta_IMppipi_p_wL_acc->Draw("colz");
 // gPad->SetLogz();
  cLpimCos->cd(4);
  //CosThetaIMppipi_p_wL_accerr->SetMaximum(0.5);
  CosTheta_IMppipi_p_wL_accerr->GetYaxis()->SetRangeUser(0,1);
  CosTheta_IMppipi_p_wL_accerr->Draw("colz");
  

  TH2F* CosTheta_IMppipi_p_wL_mc=NULL;
  CosTheta_IMppipi_p_wL_mc = (TH2F*)fLpim->Get("CosTheta_IMppipi_p_wL_nocombi_mc");
  CosTheta_IMppipi_p_wL_mc ->SetTitle("reco. evt.");
  CosTheta_IMppipi_p_wL_mc->Scale(1./SimBeamSurvivalRate);

  TH2F* CosTheta_IMppipi_p_wL_mc_acc=NULL;
  CosTheta_IMppipi_p_wL_mc_acc = (TH2F*)CosTheta_IMppipi_p_wL_mc->Clone("CosTheta_IMppipi_p_wL_mc_acc");
  CosTheta_IMppipi_p_wL_mc_acc->SetTitle("CosTheta_IMppipi_p_wL_mc_acc");
  CosTheta_IMppipi_p_wL_mc_acc->Print("base");
  CosTheta_IMppipi_p_wL_mc_acc->Divide(CosTheta_IMppipi_p_wL_mc_acc,costhetap_IMLpim_gen,1.0,1.0,"b");


  TH2F* CosTheta_IMppipi_p_wL_mc_accerr=NULL;
  CosTheta_IMppipi_p_wL_mc_accerr = (TH2F*)CosTheta_IMppipi_p_wL_mc_acc->Clone("CosTheta_IMppipi_p_wL_mc_accerr");
  CosTheta_IMppipi_p_wL_mc_accerr->Reset();
  CosTheta_IMppipi_p_wL_mc_accerr->SetTitle("CosTheta_IMppipi_p_wL_mc_acc precision");
  
  for(int ix=0;ix<CosTheta_IMppipi_p_wL_mc_accerr->GetNbinsX();ix++){
    for(int iy=0;iy<CosTheta_IMppipi_p_wL_mc_accerr->GetNbinsY();iy++){
      double cont = CosTheta_IMppipi_p_wL_mc_acc->GetBinContent(ix,iy);
      double err = CosTheta_IMppipi_p_wL_mc_acc->GetBinError(ix,iy);
      if(cont!=0){
        CosTheta_IMppipi_p_wL_mc_accerr->SetBinContent(ix,iy,err/cont);  
      }
    }
  }
  
  TCanvas *cLpimCos_mc;
  cLpimCos_mc = new TCanvas(Form("cLpimCos_mc"),Form("cLpimCos_mc"),2000,1200);
  cLpimCos_mc->Divide(2,2);
  cLpimCos_mc->cd(1);
  costhetap_IMLpim_gen->GetYaxis()->SetRangeUser(0,1);
  costhetap_IMLpim_gen->Draw("colz");
  gPad->SetLogz();
  cLpimCos_mc->cd(2);
  CosTheta_IMppipi_p_wL_mc->GetYaxis()->SetRangeUser(0,1);
  CosTheta_IMppipi_p_wL_mc->Draw("colz");
  gPad->SetLogz();
  cLpimCos_mc->cd(3);
  //CosThetaIMppipi_p_wL_mc_acc->SetMaximum(0.05);
  CosTheta_IMppipi_p_wL_mc_acc->GetYaxis()->SetRangeUser(0,1);
  CosTheta_IMppipi_p_wL_mc_acc->SetMaximum(0.05);
  CosTheta_IMppipi_p_wL_mc_acc->Draw("colz");
 // gPad->SetLogz();
  cLpimCos_mc->cd(4);
  //CosThetaIMppipi_p_wL_mc_accerr->SetMaximum(0.5);
  CosTheta_IMppipi_p_wL_mc_accerr->GetYaxis()->SetRangeUser(0,1);
  CosTheta_IMppipi_p_wL_mc_accerr->Draw("colz");

  //TCanvas *cLpimCosTest = new TCanvas("cLpimCosTest","cLpimCosTest",1000,800);
  //TH2F* CosTheta_IMppipi_p_wL_test = (TH2F*)CosTheta_IMppipi_p_wL->Clone("CosTheta_IMppipi_p_wL_test");
  //CosTheta_IMppipi_p_wL_test->Divide(CosTheta_IMppipi_p_wL_acc);
  //CosTheta_IMppipi_p_wL_test->Draw("colz");
  
  /*
  TCanvas *cLpimCosTestDiff = new TCanvas("cLpimCosTestDiff","cLpimCosTestDiff",1000,800);
  TH1D* hLpimCosTestDiff = new TH1D("hLpimCosTestDiff","hLpimCosTestDiff",10000,-0.0001,0.0001);
  for(int ix=0;ix<CosTheta_IMppipi_p_wL_test->GetNbinsX();ix++){
    for(int iy=0;iy<CosTheta_IMppipi_p_wL_test->GetNbinsY();iy++){
      double conttest = CosTheta_IMppipi_p_wL_test->GetBinContent(ix,iy);
      double contgen = costhetap_IMLpim_gen->GetBinContent(ix,iy);
      if(contgen!=0){
        hLpimCosTestDiff->Fill((conttest-contgen)/contgen*100.0);
        if(fabs((conttest-contgen)/contgen*100.0)>0.02e-03){
          CosTheta_IMppipi_p_wL_acc_clean->SetBinContent(ix,iy,0);
          CosTheta_IMppipi_p_wL_acc_clean->SetBinError(ix,iy,0);
          //std::cout << "ix:iy" << ix << " " << iy << std::endl;
        }
      }
    }
  }
  hLpimCosTestDiff->Draw("HE");
  */
  

  TH2F* CosTheta_IMppipi_p_wL_acc_clean = (TH2F*)CosTheta_IMppipi_p_wL_acc->Clone("CosTheta_IMppipi_p_wL_acc_clean");
  for(int ix=0;ix<CosTheta_IMppipi_p_wL_acc->GetNbinsX();ix++){
    for(int iy=0;iy<CosTheta_IMppipi_p_wL_acc->GetNbinsY();iy++){
      double error = CosTheta_IMppipi_p_wL_accerr->GetBinContent(ix,iy);
      if(error>0.30){
        CosTheta_IMppipi_p_wL_acc_clean->SetBinContent(ix,iy,0);
        CosTheta_IMppipi_p_wL_acc_clean->SetBinError(ix,iy,0);
      }
    }
  }

  TCanvas *cLpimCosClean = new TCanvas("cLpimCosClean","cLpimCosClean",1000,800);
  CosTheta_IMppipi_p_wL_acc_clean->Draw("colz");
        

  //acceptance map for costheta_missing p (C.M. frame) v.s. IM(pi-Lambda) 
  TH2F* costhetapCM_IMLpim_gen=NULL;
  costhetapCM_IMLpim_gen = (TH2F*)fgen->Get("React_costhetapCM_IMLPim");
  costhetapCM_IMLpim_gen->SetTitle("generated evt. ");
  costhetapCM_IMLpim_gen->SetXTitle("true IM(#pi^{-}#Lambda) [GeV/c^{2}]");
  costhetapCM_IMLpim_gen->SetYTitle("true proton CosTheta (CM)");
  costhetapCM_IMLpim_gen->GetXaxis()->CenterTitle();
  costhetapCM_IMLpim_gen->GetYaxis()->CenterTitle();
  costhetapCM_IMLpim_gen->Print("base");
  costhetapCM_IMLpim_gen->RebinX(15);
  costhetapCM_IMLpim_gen->Print("base");
  
  TH2F* CosThetaCM_IMppipi_p_wL_nocombi=NULL;
  CosThetaCM_IMppipi_p_wL_nocombi = (TH2F*)fLpim->Get("CosThetaCM_IMppipi_p_wL_nocombi");
  CosThetaCM_IMppipi_p_wL_nocombi->SetTitle("reco. evt.");
  CosThetaCM_IMppipi_p_wL_nocombi->Scale(1./SimBeamSurvivalRate);

  TH2F* CosThetaCM_IMppipi_p_wL_acc=NULL;
  CosThetaCM_IMppipi_p_wL_acc = (TH2F*)CosThetaCM_IMppipi_p_wL_nocombi->Clone("CosThetaCM_IMppipi_p_wL_acc");
  CosThetaCM_IMppipi_p_wL_acc->SetTitle("CosThetaCM_IMppipi_p_wL_acc");
  CosThetaCM_IMppipi_p_wL_acc->Print("base");
  CosThetaCM_IMppipi_p_wL_acc->Divide(CosThetaCM_IMppipi_p_wL_acc,costhetapCM_IMLpim_gen,1.0,1.0,"b");


  TH2F* CosThetaCM_IMppipi_p_wL_accerr=NULL;
  CosThetaCM_IMppipi_p_wL_accerr = (TH2F*)CosThetaCM_IMppipi_p_wL_acc->Clone("CosThetaCM_IMppipi_p_wL_accerr");
  CosThetaCM_IMppipi_p_wL_accerr->Reset();
  CosThetaCM_IMppipi_p_wL_accerr->SetTitle("CosThetaCM_IMppipi_p_wL_acc precision");
  
  for(int ix=0;ix<CosThetaCM_IMppipi_p_wL_accerr->GetNbinsX();ix++){
    for(int iy=0;iy<CosThetaCM_IMppipi_p_wL_accerr->GetNbinsY();iy++){
      double cont = CosThetaCM_IMppipi_p_wL_acc->GetBinContent(ix,iy);
      double err = CosThetaCM_IMppipi_p_wL_acc->GetBinError(ix,iy);
      if(cont!=0){
        CosThetaCM_IMppipi_p_wL_accerr->SetBinContent(ix,iy,err/cont);  
      }
    }
  } 

  /*
  TH2F* CosTheta_IMppipi_p_wL=NULL;
  CosTheta_IMppipi_p_wL = (TH2F*)fLpim->Get("CosTheta_IMppipi_p_wL_nocombi");
  CosTheta_IMppipi_p_wL ->SetTitle("reco. evt.");
  CosTheta_IMppipi_p_wL->Scale(1./SimBeamSurvivalRate);

  TH2F* CosTheta_IMppipi_p_wL_acc=NULL;
  CosTheta_IMppipi_p_wL_acc = (TH2F*)CosTheta_IMppipi_p_wL->Clone("CosTheta_IMppipi_p_wL_acc");
  CosTheta_IMppipi_p_wL_acc->SetTitle("CosTheta_IMppipi_p_wL_acc");
  CosTheta_IMppipi_p_wL_acc->Print("base");
  CosTheta_IMppipi_p_wL_acc->Divide(CosTheta_IMppipi_p_wL_acc,costhetap_IMLpim_gen,1.0,1.0,"b");


  TH2F* CosTheta_IMppipi_p_wL_accerr=NULL;
  CosTheta_IMppipi_p_wL_accerr = (TH2F*)CosTheta_IMppipi_p_wL_acc->Clone("CosTheta_IMppipi_p_wL_accerr");
  CosTheta_IMppipi_p_wL_accerr->Reset();
  CosTheta_IMppipi_p_wL_accerr->SetTitle("CosTheta_IMppipi_p_wL_acc precision");
  
  for(int ix=0;ix<CosTheta_IMppipi_p_wL_accerr->GetNbinsX();ix++){
    for(int iy=0;iy<CosTheta_IMppipi_p_wL_accerr->GetNbinsY();iy++){
      double cont = CosTheta_IMppipi_p_wL_acc->GetBinContent(ix,iy);
      double err = CosTheta_IMppipi_p_wL_acc->GetBinError(ix,iy);
      if(cont!=0){
        CosTheta_IMppipi_p_wL_accerr->SetBinContent(ix,iy,err/cont);  
      }
    }
  } 
  */

  
  TCanvas *cLpimCosCM;
  cLpimCosCM = new TCanvas("cLpimCosCM","cLpimCosCM",2000,1200);
  cLpimCosCM->Divide(2,2);
  cLpimCosCM->cd(1);
  costhetapCM_IMLpim_gen->GetYaxis()->SetRangeUser(0,1);
  costhetapCM_IMLpim_gen->Draw("colz");
  gPad->SetLogz();
  cLpimCosCM->cd(2);
  CosThetaCM_IMppipi_p_wL_nocombi->GetYaxis()->SetRangeUser(0,1);
  CosThetaCM_IMppipi_p_wL_nocombi->Draw("colz");
  gPad->SetLogz();
  cLpimCosCM->cd(3);
  //CosThetaIMppipi_p_wL_acc->SetMaximum(0.05);
  CosThetaCM_IMppipi_p_wL_acc->GetYaxis()->SetRangeUser(0,1);
  CosThetaCM_IMppipi_p_wL_acc->SetMaximum(0.05);
  CosThetaCM_IMppipi_p_wL_acc->Draw("colz");
 // gPad->SetLogz();
  cLpimCosCM->cd(4);
  //CosThetaIMppipi_p_wL_accerr->SetMaximum(0.5);
  CosThetaCM_IMppipi_p_wL_accerr->GetYaxis()->SetRangeUser(0,1);
  CosThetaCM_IMppipi_p_wL_accerr->Draw("colz");






  //  TFile *fout = new TFile("accmapLpim_pS0pim.root","RECREATE");
  TFile *fout = new TFile(Form("accmapLpimv%d.root",version),"RECREATE");
  for(int icut=0;icut<ncut;icut++){
    q_IMppipi_p_wL_acc[icut]->Write();
    q_IMppipi_p_wL_nop2_acc[icut]->Write();
    q_IMppipi_p_wL_wp2_acc[icut]->Write();
    q_IMppipi_p_wL_nop2_mc_acc[icut]->Write();
    q_IMppipi_p_wL_wp2_mc_acc[icut]->Write();
    q_IMppipi_p_wL_accerr[icut]->Write();
    q_IMppipi_p_wL_nop2_accerr[icut]->Write();
    q_IMppipi_p_wL_wp2_accerr[icut]->Write();
    q_IMppipi_p_wL_nop2_mc_accerr[icut]->Write();
    q_IMppipi_p_wL_wp2_mc_accerr[icut]->Write();
    q_IMppipi_p_wL_acc_clean[icut]->Write();
    q_IMppipi_p_wL_nop2_acc_clean[icut]->Write();
    q_IMppipi_p_wL_wp2_acc_clean[icut]->Write();
    q_IMppipi_p_wL_nop2_mc_acc_clean[icut]->Write();
    q_IMppipi_p_wL_wp2_mc_acc_clean[icut]->Write();
  }
  CosTheta_IMppipi_p_wL_acc->Write();
  CosTheta_IMppipi_p_wL_acc_clean->Write();
  CosTheta_IMppipi_p_wL_mc_acc->Write();
  fout->Close();

  TString pdfname = "accmaplpim.pdf";
  TCanvas *c = NULL;
  //TPDF *pdf = new TPDF(pdfname);
  TSeqCollection *SCol = gROOT->GetListOfCanvases();
  int size = SCol->GetSize();
  //TIter next(gROOT->GetListOfCanvases());
  TIter next(SCol);
  for(int i=0; i<size; i++) {
    //while((c= (TCanvas*)next())){
    //pdf->NewPage();
    c= (TCanvas*)next();
    c->Draw();
    c->cd();
    TPaveText *pt;
    pt = new TPaveText(.80,0.90,0.98,0.99,"NDC");
    pt->AddText("MC");
    pt->SetFillColor(kCyan-9);
    pt->SetBorderSize(1);
    pt->Draw();

    //gPad->SetLeftMargin(0.13);
    //gPad->SetBottomMargin(0.13);
    c->Modified();
    c->Update();
    if(i==0) c->Print(pdfname+"(",Form("Title:%s",c->GetTitle()));
    else if(i==size-1)c->Print(pdfname+")",Form("Title:%s",c->GetTitle()));
    else c->Print(pdfname,Form("Title:%s",c->GetTitle()));
  }

}
