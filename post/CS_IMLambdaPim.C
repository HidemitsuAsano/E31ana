void CS_IMLambdaPim()
{
  const double ForwardAngle = 0.996;
  gStyle->SetOptStat(0);
  TFile *file = TFile::Open("evanaIMLambdaPim_ppimpim_v16_out.root","READ");
  TFile *facc = TFile::Open("../simpost/accmapLpimv23.root","READ");
  TFile *flumi = TFile::Open("InteLumi.root","READ");
  TFile *fkin = TFile::Open("../simpost/NumericalRootFinderLPim.root","READ");
  TParameter<double>*IntegLumi = (TParameter<double>*)flumi->Get("IntegLumi");
  TParameter<double>*Err = (TParameter<double>*)flumi->Get("Err");
  double lumi = IntegLumi->GetVal();
  double lumierr = Err->GetVal();
  double trigScale = 0.5;
  std::cout << "Lumi:  " << lumi << std::endl;
  std::cout << "Err:   " << lumierr << std::endl;
  TCanvas *cq_IMppipi_p_wL_sum  = new TCanvas("cq_IMppipi_p_wL_sum","cq_IMppipi_p_wL_sum",1000,800);
  TH2F* q_IMppipi_p_wL_sum = (TH2F*)file->Get("q_IMppipi_p_wL_sum0");
  q_IMppipi_p_wL_sum->SetXTitle("IM(#pi^{-}#Lambda) [GeV/c^{2}]");
  q_IMppipi_p_wL_sum->Draw("colz");
  TGraph *gth = (TGraph*)fkin->Get("th");
  gth->Draw("pc");
  TGraph *gr_0 = (TGraph*)fkin->Get("gr_0");
  gr_0->Draw("pc");
  TGraph *gr_100 = (TGraph*)fkin->Get("gr_100");
  gr_100->Draw("pc");
  TGraph *gr_65 = (TGraph*)fkin->Get("gr_65");
  gr_65->Draw("pc");

  TCanvas *cq_IMppipi_p_wL_sum_f  = new TCanvas("cq_IMppipi_p_wL_sum_f","cq_IMppipi_p_wL_sum_f",1000,800);
  TH2F* q_IMppipi_p_wL_sum_f = (TH2F*)file->Get("q_IMppipi_p_wL_sum_forward");
  q_IMppipi_p_wL_sum_f->SetXTitle("IM(#pi^{-}#Lambda) [GeV/c^{2}]");
  q_IMppipi_p_wL_sum_f->Draw("colz");
  
  TCanvas *cacc = new TCanvas("cacc","cacc",1000,800);
  TH2F* q_IMppipi_p_wL_acc = (TH2F*)facc->Get("q_IMppipi_p_wL_acc_clean0");
  q_IMppipi_p_wL_acc->SetXTitle("IM(#pi^{-}#Lambda) [GeV/c^{2}]");
  q_IMppipi_p_wL_acc->Draw("colz");  
  //TGraph *gth = (TGraph*)fkin->Get("th");
  gth->Draw("pc");
  //TGraph *gr_0 = (TGraph*)fkin->Get("gr_0");
  gr_0->Draw("pc");
  //TGraph *gr_100 = (TGraph*)fkin->Get("gr_100");
  gr_100->Draw("pc");
  //TGraph *gr_65 = (TGraph*)fkin->Get("gr_65");
  gr_65->Draw("pc");


  TH2F* CS_q_IMppipi_p_wL_sum = (TH2F*)q_IMppipi_p_wL_sum->Clone("CS_q_IMppipi_p_wL_sum");
  CS_q_IMppipi_p_wL_sum->Divide(q_IMppipi_p_wL_acc);
  TCanvas *cCS = new TCanvas("cCS","cCS",1000,800);
  double binwidth = CS_q_IMppipi_p_wL_sum->ProjectionX()->GetBinWidth(1)*1000.0;
  CS_q_IMppipi_p_wL_sum->Scale(1.0/binwidth/trigScale/lumi);
//  CS_q_IMppipi_p_wL_sum->SetMaximum(0.02);
  CS_q_IMppipi_p_wL_sum->SetXTitle("IM(#pi^{-}#Lambda) [GeV/c^{2}]");
  CS_q_IMppipi_p_wL_sum->Draw("colz");

  
  TCanvas *cCS_px = new TCanvas("cCS_px","cCS_px",1000,800);
  TH1D* CS_IMppipi_p_wL_sum = (TH1D*)CS_q_IMppipi_p_wL_sum->ProjectionX("CS_IMppipi_p_wL_sum");
  CS_IMppipi_p_wL_sum->SetYTitle("d#rho/dM [#mu b /(MeV/c^{2})]");
  CS_IMppipi_p_wL_sum->Draw("HE");
  
  const int bin350 = CS_q_IMppipi_p_wL_sum->GetYaxis()->FindBin(0.35);
  TCanvas *cCS_px_0 = new TCanvas("cCS_px_0","cCS_px_0",1000,800);
  TH1D* CS_IMppipi_p_wL_sum_0 = (TH1D*)CS_q_IMppipi_p_wL_sum->ProjectionX("CS_IMppipi_p_wL_sum_0",1,bin350-1);
  //CS_IMppipi_p_wL_sum_0->GetXaxis()->SetRangeUser(1.2,1.6);
  CS_IMppipi_p_wL_sum_0->SetYTitle("d#rho/dM [#mu b /(MeV/c^{2})]");
  CS_IMppipi_p_wL_sum_0->Draw("HE");

  TCanvas *cCS_px_350 = new TCanvas("cCS_px_350","cCS_px_350",1000,800);
  TH1D* CS_IMppipi_p_wL_sum_350 = (TH1D*)CS_q_IMppipi_p_wL_sum->ProjectionX("CS_IMppipi_p_wL_sum_350",bin350,600);
  //CS_IMppipi_p_wL_sum_350->GetXaxis()->SetRangeUser(1.2,1.6);
  CS_IMppipi_p_wL_sum_350->SetYTitle("d#rho/dM [#mu b /(MeV/c^{2})]");
  CS_IMppipi_p_wL_sum_350->Draw("HE");
  
  /*
  TCanvas *cCS_f = new TCanvas("cCS_f","cCS_f",1000,800);
  TH2F* CS_q_IMppipi_p_wL_sum_f = (TH2F*)q_IMppipi_p_wL_sum_f->Clone("CS_q_IMppipi_p_wL_sum_f");
  CS_q_IMppipi_p_wL_sum_f->Divide(q_IMppipi_p_wL_acc);
  CS_q_IMppipi_p_wL_sum_f->Scale(1.0/binwidth/trigScale/lumi);
//  CS_q_IMppipi_p_wL_sum->SetMaximum(0.02);
  CS_q_IMppipi_p_wL_sum_f->SetXTitle("IM(#pi^{-}#Lambda) [GeV/c^{2}]");
  CS_q_IMppipi_p_wL_sum_f->Draw("colz");
  
  TCanvas *cCS_f_px = new TCanvas("cCS_f_px","cCS_f_px",1000,800);
  TH1D* CS_IMppipi_p_wL_sum_f = (TH1D*)CS_q_IMppipi_p_wL_sum_f->ProjectionX("CS_IMppipi_p_wL_sum_f_0");
  //CS_IMppipi_p_wL_sum_0->GetXaxis()->SetRangeUser(1.2,1.6);
  CS_IMppipi_p_wL_sum_f->SetYTitle("d#rho/dM [#mu b (MeV/c^{2})]");
  CS_IMppipi_p_wL_sum_f->Draw("HE");
  */
 
  TCanvas *cCS_q = new TCanvas("cCS_q","cCS_q",1000,800);
  const int bin1360 = CS_q_IMppipi_p_wL_sum->GetXaxis()->FindBin(1.360);
  const int bin1410 = CS_q_IMppipi_p_wL_sum->GetXaxis()->FindBin(1.410);
  std::cout << CS_q_IMppipi_p_wL_sum->GetXaxis()->GetBinLowEdge(bin1360) << std::endl;
  std::cout << CS_q_IMppipi_p_wL_sum->GetXaxis()->GetBinLowEdge(bin1410+1) << std::endl;
  TH1D* CS_q = (TH1D*)CS_q_IMppipi_p_wL_sum->ProjectionY("CS_q",bin1360,bin1410);
  CS_q->SetMarkerStyle(20);
  CS_q->SetYTitle("d#rho/dM [#mu b /(MeV/c^{2})]");
  CS_q->GetYaxis()->CenterTitle();
  CS_q->Draw("E");
  
  //divided acceptance map
  TH2F* q_IMppipi_p_wL_sum_nop2[5];//0:default, 1 half low, 2 half high, 3 sigma0 region, 4 wide range
  TH2F* q_IMppipi_p_wL_sum_wp2[5];//0:default, 1 half low, 2 half high, 3 sigma0 region, 4 wide range
  for(int icut=0;icut<5;icut++){
    q_IMppipi_p_wL_sum_nop2[icut] = (TH2F*)file->Get(Form("q_IMppipi_p_wL_sum_nop2%d",icut));
    q_IMppipi_p_wL_sum_wp2[icut]  = (TH2F*)file->Get(Form("q_IMppipi_p_wL_sum_wp2%d",icut));
  }
  
  TCanvas *cq_IMppipi_p_wL_sum_nop2[5];//0:default, 1 half low, 2 half high, 3 sigma0 region, 4 wide range
  TCanvas *cq_IMppipi_p_wL_sum_wp2[5];//0:default, 1 half low, 2 half high, 3 sigma0 region, 4 wide range
  for(int icut=0;icut<5;icut++){
    cq_IMppipi_p_wL_sum_nop2[icut] = new TCanvas(Form("cq_IMppipi_p_wL_sum_nop2%d",icut),Form("cq_IMppipi_p_wL_sum_nop2%d",icut),1000,800);
    q_IMppipi_p_wL_sum_nop2[icut]->Draw("colz");
    cq_IMppipi_p_wL_sum_wp2[icut] = new TCanvas(Form("cq_IMppipi_p_wL_sum_wp2%d",icut),Form("cq_IMppipi_p_wL_sum_wp2%d",icut),1000,800);
    q_IMppipi_p_wL_sum_wp2[icut]->Draw("colz");
  }

  //get acceptance map from file
  TH2F* q_IMppipi_p_wL_nop2_acc[5];//normal acc.
  TH2F* q_IMppipi_p_wL_wp2_acc[5];//normal acc.
  TH2F* q_IMppipi_p_wL_nop2_mc_acc[5];//with true val.
  TH2F* q_IMppipi_p_wL_wp2_mc_acc[5];//with true val.
  for(int icut=0;icut<5;icut++){
    q_IMppipi_p_wL_nop2_acc[icut] = (TH2F*)facc->Get(Form("q_IMppipi_p_wL_nop2_acc_clean%d",icut));
    q_IMppipi_p_wL_wp2_acc[icut] = (TH2F*)facc->Get(Form("q_IMppipi_p_wL_wp2_acc_clean%d",icut));
    q_IMppipi_p_wL_nop2_mc_acc[icut] = (TH2F*)facc->Get(Form("q_IMppipi_p_wL_nop2_mc_acc_clean%d",icut));
    q_IMppipi_p_wL_wp2_mc_acc[icut] = (TH2F*)facc->Get(Form("q_IMppipi_p_wL_wp2_mc_acc_clean%d",icut));
  }

  TH2F* CS_q_IMppipi_p_wL_nop2_acc[5];//normal acc.
  TH2F* CS_q_IMppipi_p_wL_wp2_acc[5];//normal acc.
  TH2F* CS_q_IMppipi_p_wL_nop2_mc_acc[5];//with true val.
  TH2F* CS_q_IMppipi_p_wL_wp2_mc_acc[5];//with true val.

  for(int icut=0;icut<5;icut++){
    CS_q_IMppipi_p_wL_nop2_acc[icut] = (TH2F*) q_IMppipi_p_wL_sum_nop2[icut]->Clone(Form("CS_q_IMppipi_p_wL_nop2_acc%d",icut));
    CS_q_IMppipi_p_wL_wp2_acc[icut] = (TH2F*)q_IMppipi_p_wL_sum_wp2[icut]->Clone(Form("CS_q_IMppipi_p_wL_wp2_acc%d",icut));
    CS_q_IMppipi_p_wL_nop2_mc_acc[icut] = (TH2F*)q_IMppipi_p_wL_sum_nop2[icut]->Clone(Form("CS_q_IMppipi_p_wL_nop2_mc_acc%d",icut));
    CS_q_IMppipi_p_wL_wp2_mc_acc[icut] = (TH2F*)q_IMppipi_p_wL_sum_wp2[icut]->Clone(Form("CS_q_IMppipi_p_wL_wp2_mc_acc%d",icut));
  
    CS_q_IMppipi_p_wL_nop2_acc[icut]->Divide(q_IMppipi_p_wL_nop2_acc[icut]);
    CS_q_IMppipi_p_wL_wp2_acc[icut]->Divide(q_IMppipi_p_wL_wp2_acc[icut]);
    CS_q_IMppipi_p_wL_nop2_mc_acc[icut]->Divide(q_IMppipi_p_wL_nop2_mc_acc[icut]);
    CS_q_IMppipi_p_wL_wp2_mc_acc[icut]->Divide(q_IMppipi_p_wL_wp2_mc_acc[icut]);
    CS_q_IMppipi_p_wL_nop2_acc[icut]->Scale(1.0/binwidth/trigScale/lumi);
    CS_q_IMppipi_p_wL_nop2_acc[icut]->SetXTitle("IM(#pi^{-}#Lambda) [GeV/c^{2}]");
    CS_q_IMppipi_p_wL_wp2_acc[icut]->Scale(1.0/binwidth/trigScale/lumi);
    CS_q_IMppipi_p_wL_wp2_acc[icut]->SetXTitle("IM(#pi^{-}#Lambda) [GeV/c^{2}]");
    CS_q_IMppipi_p_wL_nop2_mc_acc[icut]->Scale(1.0/binwidth/trigScale/lumi);
    CS_q_IMppipi_p_wL_nop2_mc_acc[icut]->SetXTitle("IM(#pi^{-}#Lambda) [GeV/c^{2}]");
    CS_q_IMppipi_p_wL_wp2_mc_acc[icut]->Scale(1.0/binwidth/trigScale/lumi);
    CS_q_IMppipi_p_wL_wp2_mc_acc[icut]->SetXTitle("IM(#pi^{-}#Lambda) [GeV/c^{2}]");
  }

  TCanvas *cCS_nop2[5];
  TCanvas *cCS_nop2_mc[5];
  TCanvas *cCS_wp2[5];
  TCanvas *cCS_wp2_mc[5];

  for(int icut=0;icut<5;icut++){
    cCS_nop2[icut] = new TCanvas(Form("cCS_nop2%d",icut),Form("cCS_nop2%d",icut),1000,800);
    cCS_nop2[icut]->cd();
    CS_q_IMppipi_p_wL_nop2_acc[icut]->SetMaximum(0.18);
    CS_q_IMppipi_p_wL_nop2_acc[icut]->Draw("colz");
    
    
    cCS_nop2_mc[icut] = new TCanvas(Form("cCS_nop2_mc%d",icut),Form("cCS_nop2_mc%d",icut),1000,800);
    cCS_nop2_mc[icut]->cd();
    CS_q_IMppipi_p_wL_nop2_mc_acc[icut]->SetMaximum(0.18);
    CS_q_IMppipi_p_wL_nop2_mc_acc[icut]->Draw("colz");


    cCS_wp2[icut] = new TCanvas(Form("cCS_wp2%d",icut),Form("cCS_wp2%d",icut),1000,800);
    cCS_wp2[icut]->cd();
    CS_q_IMppipi_p_wL_wp2_acc[icut]->Draw("colz");

    cCS_wp2_mc[icut] = new TCanvas(Form("cCS_wp2_mc%d",icut),Form("cCS_wp2_mc%d",icut),1000,800);
    cCS_wp2_mc[icut]->cd();
    CS_q_IMppipi_p_wL_wp2_mc_acc[icut]->Draw("colz");
  }

  TCanvas *cCS_q_sum[5];
  TH1D* CS_q_nop2[5];
  TH1D* CS_q_wp2[5];
  TH1D* CS_q_nop2_mc[5];
  TH1D* CS_q_wp2_mc[5];
  for(int icut=0;icut<5;icut++){
    cCS_q_sum[icut]  = new TCanvas(Form("cCS_q_sum%d",icut),Form("cCS_q_sum%d",icut),1000,800);
    CS_q_nop2[icut] = (TH1D*)CS_q_IMppipi_p_wL_nop2_acc[icut]->ProjectionY(Form("CS_q_nop2%d",icut),bin1360,bin1410);
    CS_q_nop2_mc[icut] = (TH1D*)CS_q_IMppipi_p_wL_nop2_mc_acc[icut]->ProjectionY(Form("CS_q_nop2_mc%d",icut),bin1360,bin1410);
    CS_q_wp2[icut] = (TH1D*)CS_q_IMppipi_p_wL_wp2_acc[icut]->ProjectionY(Form("CS_q_wp2%d",icut),bin1360,bin1410);
    CS_q_wp2_mc[icut] = (TH1D*)CS_q_IMppipi_p_wL_wp2_mc_acc[icut]->ProjectionY(Form("CS_q_wp2_mc%d",icut),bin1360,bin1410);
    CS_q_nop2[icut]->SetMarkerStyle(21);
    CS_q_nop2[icut]->Draw("E");
    CS_q_wp2[icut]->SetMarkerStyle(20);
    CS_q_wp2[icut]->Draw("Esame");
    CS_q_nop2_mc[icut]->SetMarkerStyle(22);
    CS_q_nop2_mc[icut]->SetMarkerColor(2);
    CS_q_nop2_mc[icut]->SetLineColor(2);
    CS_q_nop2_mc[icut]->Draw("Esame");
    CS_q_wp2_mc[icut]->SetMarkerStyle(23);
    CS_q_wp2_mc[icut]->SetMarkerColor(2);
    CS_q_wp2_mc[icut]->SetLineColor(2);
    CS_q_wp2_mc[icut]->Draw("Esame");
  }


  //reco.
  TCanvas *cIMLpim_coscut = new TCanvas("cIMLpim_coscut","cIMLpim_coscut",1000,800);
  TH2F* CosTheta_IMppipi_p_wL_sum=NULL;
  CosTheta_IMppipi_p_wL_sum = (TH2F*)file->Get("CosTheta_IMppipi_p_wL_sum");
  CosTheta_IMppipi_p_wL_sum->SetTitle("reco. evt.");
  CosTheta_IMppipi_p_wL_sum->GetYaxis()->SetRangeUser(0.95,1);
  CosTheta_IMppipi_p_wL_sum->Draw("colz"); 
  
  TCanvas *cIMLpim_coscut_px = new TCanvas("cIMLpim_coscut_px","cIMLpim_coscut_px",1000,800);
  const int bin0995 = CosTheta_IMppipi_p_wL_sum->GetYaxis()->FindBin(ForwardAngle);
  const int bin1000 = CosTheta_IMppipi_p_wL_sum->GetYaxis()->FindBin(1.0);
  CosTheta_IMppipi_p_wL_sum->ProjectionX("IMppipi_p_cosf",bin0995,bin1000+1)->Draw("HE");
  
  //acc.
  TCanvas *cCosacc = new TCanvas("cCosacc","cCosacc",1000,800);
  TH2F* CosTheta_IMppipi_p_wL_acc = (TH2F*)facc->Get("CosTheta_IMppipi_p_wL_acc");
  CosTheta_IMppipi_p_wL_acc->GetYaxis()->SetRangeUser(0.95,1);
  CosTheta_IMppipi_p_wL_acc->SetMaximum(0.05);
  CosTheta_IMppipi_p_wL_acc->Draw("colz");

  TCanvas *cCosacc_mc = new TCanvas("cCosacc_mc","cCosacc_mc",1000,800);
  TH2F* CosTheta_IMppipi_p_wL_mc_acc = (TH2F*)facc->Get("CosTheta_IMppipi_p_wL_mc_acc");
  CosTheta_IMppipi_p_wL_mc_acc->GetYaxis()->SetRangeUser(0.95,1);
  CosTheta_IMppipi_p_wL_mc_acc->SetMaximum(0.05);
  CosTheta_IMppipi_p_wL_mc_acc->Draw("colz");

  TCanvas *cCS_CosTheta_IMppipi_p_wL = new TCanvas("cCS_CosTheta_IMppipi_p_wL","cCS_CosTheta_IMppipi_p_wL",1000,800);
  TH2F* CS_CosTheta_IMppipi_p_wL = (TH2F*)CosTheta_IMppipi_p_wL_sum->Clone("CS_CosTheta_IMppipi_p_wL");
  CS_CosTheta_IMppipi_p_wL->Divide(CosTheta_IMppipi_p_wL_acc);
//  CS_CosTheta_IMppipi_p_wL->Scale(1.0/binwidth/trigScale/lumi/0.0314);
  const double solidAngleCoscut = 0.02512;
  CS_CosTheta_IMppipi_p_wL->Scale(1.0/binwidth/trigScale/lumi/solidAngleCoscut);
  CS_CosTheta_IMppipi_p_wL->Draw("colz");

  TCanvas *cCS_CosTheta_IMppipi_p_wL_mc = new TCanvas("cCS_CosTheta_IMppipi_p_wL_mc","cCS_CosTheta_IMppipi_p_wL_mc",1000,800);
  TH2F* CS_CosTheta_IMppipi_p_wL_mc = (TH2F*)CosTheta_IMppipi_p_wL_sum->Clone("CS_CosTheta_IMppipi_p_wL_mc");
  CS_CosTheta_IMppipi_p_wL_mc->Divide(CosTheta_IMppipi_p_wL_mc_acc);
//  CS_CosTheta_IMppipi_p_wL->Scale(1.0/binwidth/trigScale/lumi/0.0314);
  CS_CosTheta_IMppipi_p_wL_mc->Scale(1.0/binwidth/trigScale/lumi/solidAngleCoscut);
  CS_CosTheta_IMppipi_p_wL_mc->Draw("colz");

  TCanvas *cCS_IMppipi_p_wL_coscut = new TCanvas("cCS_IMppipi_p_wL_coscut","cCS_IMppipi_p_wL_coscut",1000,800);
  TH1D* CS_IMppipi_p_wL_coscut = (TH1D*)CS_CosTheta_IMppipi_p_wL->ProjectionX("CS_IMppipi_p_wL_coscut",bin0995,bin1000);
  CS_IMppipi_p_wL_coscut->GetXaxis()->SetRangeUser(1.3,1.6);
  CS_IMppipi_p_wL_coscut->Draw("HE");
  
  //TCanvas *cCS_IMppipi_p_wL_mc_coscut = new TCanvas("cCS_IMppipi_p_wL_mc_coscut","cCS_IMppipi_p_wL_mc_coscut",1000,800);
  TH1D* CS_IMppipi_p_wL_mc_coscut = (TH1D*)CS_CosTheta_IMppipi_p_wL_mc->ProjectionX("CS_IMppipi_p_wL_mc_coscut",bin0995,bin1000);
  CS_IMppipi_p_wL_mc_coscut->GetXaxis()->SetRangeUser(1.3,1.6);
  //CS_IMppipi_p_wL_mc_coscut->Draw("HE");
  CS_IMppipi_p_wL_mc_coscut->SetLineColor(2);
  CS_IMppipi_p_wL_mc_coscut->Draw("HEsame");

  TCanvas *cCompInoue= new TCanvas("cCompInoue","cCompInoue",1000,800);  
  TFile *finoue = TFile::Open("pimL_CS.root","READ");
  CS_IMppipi_p_wL_coscut->SetTitle("C.S.");
  CS_IMppipi_p_wL_coscut->SetYTitle("d^{2}#rho/d#Omega dM [#mu b/sr /(MeV/c^{2})]");
  CS_IMppipi_p_wL_coscut->GetYaxis()->CenterTitle();
  CS_IMppipi_p_wL_coscut->Draw("HE");
  CS_IMppipi_p_wL_mc_coscut->Draw("HEsame");
  TGraphErrors *grinoue = (TGraphErrors*)finoue->Get("Graph");
  grinoue->Draw("P");

  TCanvas *c = NULL;
  TSeqCollection *SCol = gROOT->GetListOfCanvases();
  int size = SCol->GetSize();
  TIter next(SCol);
  TString pdfname = "lpim.pdf";
  for(int i=0;i<size;i++){
    //pdf->NewPage();
    c= (TCanvas*)next();
    c->Modified();
    c->Update();
    c->Draw();
    c->cd();
    //inside the canvas
    //TPaveText *pt = new TPaveText(.74,.81,0.9,0.90,"NDC");
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

  TObject *obj = NULL;
  TFile *fout = new TFile("cs_lpim_killcombi.root","RECREATE");
  TIter nexthist2(gDirectory->GetList());
  fout->Print();
  fout->cd();
  while( (obj = (TObject*)nexthist2())!=NULL) {
    obj->Write();
  }
  CS_q_IMppipi_p_wL_sum->Write();
  CS_q->Write();
  //CS_IMppipi_p_wL_sum_0->Write();
  //CS_IMppipi_p_wL_sum_350->Write();
  fout->Close();
}
