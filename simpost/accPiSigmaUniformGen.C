void accPiSigmaUniformGen(){

  //TFile *file1 = new TFile("simIMpisigma_nSppim_pippimn_DoraAir_v57_v58_acc.root","READ");
  //TFile *file1 = new TFile("simIMpisigma_nSmpip_pippimn_DoraAir_v57_v58_acc.root","READ");
  TFile *file1 = new TFile("simIMpisigma_nSmpip_pippimn_DoraAir_v57_v58_acc_v2.root","READ");
  
  bool Spmode = (std::string(file1->GetName()).find("Sp")!= std::string::npos);
  bool Smmode = (std::string(file1->GetName()).find("Sm")!= std::string::npos);
  gStyle->SetStatX(0.9);     
  gStyle->SetStatY(0.9);      
  gStyle->SetPalette(1);
  gStyle->SetStatBorderSize(1);
  gStyle->SetCanvasDefH(800); gStyle->SetCanvasDefW(1000);
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadLeftMargin(0.12);
  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(0);
  
  TH2F* q_IMpiSigma_gen_2 = (TH2F*)file1->Get("q_IMpiSigma_gen");
  q_IMpiSigma_gen_2->SetName("gen2");
  TCanvas *cq_IMpiSigma_gen_2 = new TCanvas("cq_IMpiSigma_gen_2","q_IMpiSigma_gen_2");
  q_IMpiSigma_gen_2->Draw("colz");
  
  //generated info.
  TH2F* q_IMpiSigma_gen = (TH2F*)file1->Get("q_IMpiSigma_gen");
  TCanvas *cq_IMpiSigma_gen = new TCanvas("cq_IMpiSigma_gen","q_IMpiSigma_gen");
  q_IMpiSigma_gen->RebinX(5);
  q_IMpiSigma_gen->RebinY(12);
  q_IMpiSigma_gen->Draw("colz");
  
  
  //true mass & q  including K0, 
  TH2F* q_IMpiSigma_wSid_n_genacc = (TH2F*)file1->Get("q_IMpiSigma_wSid_n_genacc");
  TCanvas *cq_IMpiSigma_wSid_n_genacc = new TCanvas("cq_IMpiSigma_wSid_n_genacc","q_IMpiSigma_wSid_n_genacc");
  q_IMpiSigma_wSid_n_genacc->RebinX(5);
  q_IMpiSigma_wSid_n_genacc->RebinY(12);
  q_IMpiSigma_wSid_n_genacc->Draw("colz");
  
  //true mass & q, removing K0
  TH2F* q_IMpiSigma_woK0_wSid_n_genacc = (TH2F*)file1->Get("q_IMpiSigma_woK0_wSid_n_genacc");
  TCanvas *cq_IMpiSigma_woK0_wSid_n_genacc = new TCanvas("cq_IMpiSigma_woK0_wSid_n_genacc","q_IMpiSigma_woK0_wSid_n_genacc");
  q_IMpiSigma_woK0_wSid_n_genacc->RebinX(5);
  q_IMpiSigma_woK0_wSid_n_genacc->RebinY(12);
  q_IMpiSigma_woK0_wSid_n_genacc->Draw("colz");
  
  //true mass & q, removing K0, Sp/Sm selection
  TH2F* q_IMpiSigma_woK0_wSid_n_SpSm_genacc; 
  if(Spmode) q_IMpiSigma_woK0_wSid_n_SpSm_genacc = (TH2F*)file1->Get("q_IMpiSigma_woK0_wSid_n_Sp_genacc");
  if(Smmode) q_IMpiSigma_woK0_wSid_n_SpSm_genacc = (TH2F*)file1->Get("q_IMpiSigma_woK0_wSid_n_Sm_genacc");
  TCanvas *cq_IMpiSigma_woK0_wSid_n_SpSm_genacc = new TCanvas("cq_IMpiSigma_woK0_wSid_n_SpSm_genacc","q_IMpiSigma_woK0_wSid_n_SpSm_genacc");
  q_IMpiSigma_woK0_wSid_n_SpSm_genacc->RebinX(5);
  q_IMpiSigma_woK0_wSid_n_SpSm_genacc->RebinY(12);
  q_IMpiSigma_woK0_wSid_n_SpSm_genacc->Draw("colz");

  //reco mass & q including K0
  TH2F* q_IMnpipi_wSid_n_acc_reco = (TH2F*)file1->Get("q_IMnpipi_wSid_n_acc_reco");
  TCanvas *cq_IMnpipi_wSid_n_acc_reco = new TCanvas("cq_IMnpipi_wSid_n_acc_reco","q_IMnpipi_wSid_n_acc_reco");
  q_IMnpipi_wSid_n_acc_reco->RebinX(5);
  q_IMnpipi_wSid_n_acc_reco->RebinY(12);
  q_IMnpipi_wSid_n_acc_reco->Draw("colz");
 
  //reco mass & q, removing K0
  TH2F* q_IMnpipi_woK0_wSid_n_acc_reco = (TH2F*)file1->Get("q_IMnpipi_woK0_wSid_n_acc_reco");
  TCanvas *cq_IMnpipi_woK0_wSid_n_acc_reco = new TCanvas("cq_IMnpipi_woK0_wSid_n_acc_reco","q_IMnpipi_woK0_wSid_n_acc_reco");
  q_IMnpipi_woK0_wSid_n_acc_reco->RebinX(5);
  q_IMnpipi_woK0_wSid_n_acc_reco->RebinY(12);
  q_IMnpipi_woK0_wSid_n_acc_reco->Draw("colz");
  
  //reco mass & q, including K0, Sp/Sm selection
  TH2F* q_IMnpipi_wSid_n_SpSm_acc_reco;
  if(Spmode)q_IMnpipi_wSid_n_SpSm_acc_reco = (TH2F*)file1->Get("q_IMnpipi_wSid_n_Sp_acc_reco");
  if(Smmode)q_IMnpipi_wSid_n_SpSm_acc_reco = (TH2F*)file1->Get("q_IMnpipi_wSid_n_Sm_acc_reco");
  TCanvas *cq_IMnpipi_wSid_n_SpSm_acc_reco = new TCanvas("cq_IMnpipi_wSid_n_SpSm_acc_reco","q_IMnpipi_wSid_n_SpSm_acc_reco");
  q_IMnpipi_wSid_n_SpSm_acc_reco->RebinX(5);
  q_IMnpipi_wSid_n_SpSm_acc_reco->RebinY(12);
  q_IMnpipi_wSid_n_SpSm_acc_reco->Draw("colz");

  //reco mass & q, removing K0, Sp/Sm selection
  TH2F* q_IMnpipi_woK0_wSid_n_SpSm_acc_reco;
  if(Spmode)q_IMnpipi_woK0_wSid_n_SpSm_acc_reco = (TH2F*)file1->Get("q_IMnpipi_woK0_wSid_n_Sp_acc_reco");
  if(Smmode)q_IMnpipi_woK0_wSid_n_SpSm_acc_reco = (TH2F*)file1->Get("q_IMnpipi_woK0_wSid_n_Sm_acc_reco");
  TCanvas *cq_IMnpipi_woK0_wSid_n_SpSm_acc_reco = new TCanvas("cq_IMnpipi_woK0_wSid_n_SpSm_acc_reco","q_IMnpipi_woK0_wSid_n_SpSm_acc_reco");
  q_IMnpipi_woK0_wSid_n_SpSm_acc_reco->RebinX(5);
  q_IMnpipi_woK0_wSid_n_SpSm_acc_reco->RebinY(12);
  q_IMnpipi_woK0_wSid_n_SpSm_acc_reco->Draw("colz");
  
  
  /*
  TEfficiency *pEff;
  pEff = new TEfficiency(*q_IMnpipi_woK0_wSid_n_Sp_acc_add,*React_q_IMPiSigma_Sp_add);
  pEff->SetTitle("eff Sp test");
  pEff->SetName("heff_Sp_add");
  TCanvas *ceff = new TCanvas("ceff","ceff");
  ceff->cd();
  pEff->Draw("EY");
  ceff->Update();
  //TCanvas *ceff_hist = new TCanvas("ceff_hist","ceff_hist");
  auto heff = pEff->GetPaintedHistogram();
  heff->SetMinimum(0.0);
  heff->SetMaximum(0.008);
  ceff->Update();
  heff->Draw("colz");
  heff->Print("base");
  ceff->Modified();
  */
  
  
  TH2F* heff;// = new TH2F("heff","heff",125,1,2,50,0,1.5);
  heff = (TH2F*)q_IMpiSigma_wSid_n_genacc->Clone();
  heff->SetName("eff_q_IMpiSigma_wSid_n");
  heff->SetTitle("eff_q_IMpiSigma_wSid_n");
  TCanvas *ceff_hist = new TCanvas("ceff_hist","ceff_hist");
  heff->Divide(q_IMpiSigma_wSid_n_genacc,q_IMpiSigma_gen,1.0,1.0,"b");
  if(Spmode)heff->SetMaximum(0.005);
  if(Smmode)heff->SetMaximum(0.009);
  heff->Draw("colz");
  
  TH2F* heff_reco;// = new TH2F("heff","heff",125,1,2,50,0,1.5);
  heff_reco = (TH2F*)q_IMnpipi_wSid_n_acc_reco->Clone();
  heff_reco->SetName("eff_q_IMpiSigma_wSid_n_reco");
  heff_reco->SetTitle("eff_q_IMpiSigma_wSid_n_reco");
  TCanvas *ceff_hist_reco = new TCanvas("ceff_hist_reco","ceff_hist_reco");
  heff_reco->Divide(q_IMnpipi_wSid_n_acc_reco,q_IMpiSigma_gen,1.0,1.0,"b");
  if(Spmode)heff_reco->SetMaximum(0.005);
  if(Smmode)heff_reco->SetMaximum(0.009);
  heff_reco->Draw("colz");
  
  TH2F* heff_woK0;// = new TH2F("heff","heff",125,1,2,50,0,1.5);
  heff_woK0 = (TH2F*)q_IMpiSigma_wSid_n_genacc->Clone();
  heff_woK0->SetName("eff_q_IMpiSigma_woK0_wSid_n");
  heff_woK0->SetTitle("eff_q_IMpiSigma_woK0_wSid_n");
  TCanvas *ceff_woK0_hist = new TCanvas("ceff_woK0_hist","ceff_woK0_hist");
  heff_woK0->Divide(q_IMpiSigma_woK0_wSid_n_genacc,q_IMpiSigma_gen,1.0,1.0,"b");
  if(Spmode)heff_woK0->SetMaximum(0.005);
  if(Smmode)heff_woK0->SetMaximum(0.009);
  heff_woK0->Draw("colz");
  
  TH2F* heff_woK0_reco;// = new TH2F("heff","heff",125,1,2,50,0,1.5);
  heff_woK0_reco = (TH2F*)q_IMnpipi_woK0_wSid_n_acc_reco->Clone();
  heff_woK0_reco->SetName("eff_q_IMpiSigma_woK0_wSid_n_reco");
  heff_woK0_reco->SetTitle("eff_q_IMpiSigma_woK0_wSid_n_reco");
  TCanvas *ceff_woK0_hist_reco = new TCanvas("ceff_woK0_hist_reco","ceff_woK0_hist_reco");
  heff_woK0_reco->Divide(q_IMnpipi_woK0_wSid_n_acc_reco,q_IMpiSigma_gen,1.0,1.0,"b");
  if(Spmode)heff_woK0_reco->SetMaximum(0.005);
  if(Smmode)heff_woK0_reco->SetMaximum(0.009);
  heff_woK0_reco->Draw("colz");
  
  TH2F* heff_SpSm;// = new TH2F("heff","heff",125,1,2,50,0,1.5);
  heff_SpSm = (TH2F*)q_IMpiSigma_wSid_n_genacc->Clone();
  heff_SpSm->SetName("eff_q_IMpiSigma_wSid_n_SpSm");
  heff_SpSm->SetTitle("eff_q_IMpiSigma_wSid_n_SpSm");
  TCanvas *ceff_SpSm_hist = new TCanvas("ceff_SpSm_hist","ceff_SpSm_hist");
  heff_SpSm->Divide(q_IMpiSigma_wSid_n_SpSm_genacc,q_IMpiSigma_gen,1.0,1.0,"b");
  if(Spmode)heff_SpSm->SetMaximum(0.005);
  if(Smmode)heff_SpSm->SetMaximum(0.009);
  heff_SpSm->Draw("colz");
  
  TH2F* heff_SpSm_reco;// = new TH2F("heff","heff",125,1,2,50,0,1.5);
  heff_SpSm_reco = (TH2F*)q_IMpiSigma_wSid_n_genacc->Clone();
  heff_SpSm_reco->SetName("eff_q_IMpiSigma_wSid_n_SpSm_reco");
  heff_SpSm_reco->SetTitle("eff_q_IMpiSigma_wSid_n_SpSm_reco");
  TCanvas *ceff_SpSm_reco_hist = new TCanvas("ceff_SpSm_reco_hist","ceff_SpSm_reco_hist");
  heff_SpSm_reco->Divide(q_IMpiSigma_wSid_n_SpSm_acc_reco,q_IMpiSigma_gen,1.0,1.0,"b");
  if(Spmode)heff_SpSm_reco->SetMaximum(0.005);
  if(Smmode)heff_SpSm_reco->SetMaximum(0.009);
  heff_SpSm_reco->Draw("colz");


  TH2F* heff_woK0_SpSm;// = new TH2F("heff","heff",125,1,2,50,0,1.5);
  heff_woK0_SpSm = (TH2F*)q_IMpiSigma_wSid_n_genacc->Clone();
  heff_woK0_SpSm->SetName("eff_q_IMpiSigma_woK0_wSid_n_SpSm");
  heff_woK0_SpSm->SetTitle("eff_q_IMpiSigma_woK0_wSid_n_SpSm");
  TCanvas *ceff_woK0_SpSm_hist = new TCanvas("ceff_woK0_SpSm_hist","ceff_woK0_SpSm_hist");
  heff_woK0_SpSm->Divide(q_IMpiSigma_woK0_wSid_n_SpSm_genacc,q_IMpiSigma_gen,1.0,1.0,"b");
  if(Spmode)heff_woK0_SpSm->SetMaximum(0.005);
  if(Smmode)heff_woK0_SpSm->SetMaximum(0.009);
  heff_woK0_SpSm->Draw("colz");

  TH2F* heff_woK0_SpSm_reco;
  heff_woK0_SpSm_reco = (TH2F*)q_IMnpipi_woK0_wSid_n_SpSm_acc_reco->Clone();
  heff_woK0_SpSm_reco->SetName("eff_q_IMpiSigma_woK0_wSid_n_SpSm_reco");
  heff_woK0_SpSm_reco->SetTitle("eff_q_IMpiSigma_woK0_wSid_n_SpSm_reco");
  TCanvas *ceff_woK0_SpSm_reco_hist = new TCanvas("ceff_woK0_SpSm_reco_hist","ceff_woK0_SpSm_reco_hist");
  heff_woK0_SpSm_reco->Divide(q_IMnpipi_woK0_wSid_n_SpSm_acc_reco,q_IMpiSigma_gen,1.0,1.0,"b");
  if(Spmode)heff_woK0_SpSm_reco->SetMaximum(0.005);
  if(Smmode)heff_woK0_SpSm_reco->SetMaximum(0.009);
  heff_woK0_SpSm_reco->Draw("colz");

  
  TH2F* heff_woK0_SpSm_clone = (TH2F*)heff_woK0_SpSm->Clone();
  heff_woK0_SpSm_clone->SetName("eff_q_IMpiSigma_woK0_wSid_n_SpSm_clone");
  heff_woK0_SpSm_clone->SetTitle("eff_q_IMpiSigma_woK0_wSid_n_SpSm_clone");
  TCanvas *ceff_woK0_SpSm_hist2 = new TCanvas("ceff_woK0_SpSm_hist2","ceff_woK0_SpSm_hist2");
  if(Spmode)heff_woK0_SpSm_clone->SetMaximum(0.005);
  if(Smmode)heff_woK0_SpSm_clone->SetMaximum(0.009);
  heff_woK0_SpSm_clone->Draw("colz");

  TFile *fnu = new TFile("NumericalRootFinder_Spmode.root");
  fnu->cd();
  TMultiGraph *mg = (TMultiGraph*)fnu->Get("mg"); 
  mg->Draw("c");

  
  TH2F* heff_err;// = new TH2F("heff_err","heff_err",125,1,2,50,0,1.5);
  //heff_err = (TH2F*)heff_woK0_SpSm->Clone();
  heff_err = (TH2F*)heff_woK0_SpSm_reco->Clone();
  heff_err->SetName("heff_err");
  heff_err->SetTitle("heff_err");
  for(int ix=0;ix<heff_woK0_SpSm_reco->GetNbinsX();ix++){
    for(int iy=0;iy<heff_woK0_SpSm_reco->GetNbinsY();iy++){
      double err = heff_woK0_SpSm_reco->GetBinErrorUp(ix,iy);
      double cont = heff_woK0_SpSm_reco->GetBinContent(ix,iy);
      if(cont>0.00)heff_err->SetBinContent(ix,iy,err/cont);
    }
  }
  heff_err->SetMaximum(1.0);
  TCanvas *cacc_err = new TCanvas("cacc_err","acc_err");
  cacc_err->cd();
  heff_err->Draw("colz");
  
  //cleaning up
  
  for(int ibinx=0;ibinx<heff_err->GetNbinsX();ibinx++){
    for(int ibiny=0;ibiny<heff_err->GetNbinsY();ibiny++){
      double cont =  heff_woK0_SpSm_reco->GetBinContent(ibinx,ibiny);
      double precision =  heff_err->GetBinContent(ibinx,ibiny);
      if(precision>0.2){ 
        heff_reco->SetBinContent(ibinx,ibiny,0.0);
        heff_reco->SetBinError(ibinx,ibiny,0.0);
        heff_woK0_reco->SetBinContent(ibinx,ibiny,0.0);
        heff_woK0_reco->SetBinError(ibinx,ibiny,0.0);
        heff_SpSm_reco->SetBinContent(ibinx,ibiny,0.0);
        heff_SpSm_reco->SetBinError(ibinx,ibiny,0.0);
        heff_woK0_SpSm_reco->SetBinContent(ibinx,ibiny,0.0);
        heff_woK0_SpSm_reco->SetBinError(ibinx,ibiny,0.0);
      }
    }
  }


  
  TCanvas *ceff_woK0_SpSm_hist_px_q0 = new TCanvas("ceff_woK0_SpSm_hist_px_q0","eff_woK0_SpSm_hist_px_q0");
  //int ybin = heff_woK0_SpSm->GetYaxis()->FindBin(0.35);
  int ybin = heff_woK0_SpSm_reco->GetYaxis()->FindBin(0.35);
  heff_reco->ProjectionX("px0",0,ybin-1);
  heff_SpSm_reco->ProjectionX("px0_SpSm",0,ybin-1);
  heff_woK0_reco->ProjectionX("px0_woK0",0,ybin-1);
  heff_woK0_SpSm_reco->ProjectionX("px0_woK0_SpSm",0,ybin-1);
  px0->SetYTitle("acceptance X efficiency");
  px0->GetYaxis()->CenterTitle();
  px0->Draw();
  px0_woK0->SetLineColor(2);
  px0_woK0->Draw("same");
  px0_woK0_SpSm->SetLineColor(3);
  px0_woK0_SpSm->Draw("same");
  px0_SpSm->SetLineColor(4);
  px0_SpSm->Draw("same");
  
  TCanvas *ceff_woK0_SpSm_hist_px_q350 = new TCanvas("ceff_woK0_SpSm_hist_px_q350","eff_woK0_SpSm_hist_px_q350");
  //int ybin = heff_woK0_SpSm->GetYaxis()->FindBin(0.35);
  heff_reco->ProjectionX("px350",ybin,50);
  heff_woK0_reco->ProjectionX("px350_woK0",ybin,50);
  heff_woK0_SpSm_reco->ProjectionX("px350_woK0_SpSm",ybin,50);
  px350->SetYTitle("acceptance X efficiency");
  px350->GetYaxis()->CenterTitle();
  px350->Draw();
  px350_woK0->SetLineColor(2);
  px350_woK0->Draw("same");
  px350_woK0_SpSm->SetLineColor(3);
  px350_woK0_SpSm->Draw("same");
  px350_SpSm->SetLineColor(4);
  px350_SpSm->Draw("same");

  

  TIter nexthist(gDirectory->GetList());
  TH1F *h1 = NULL;
  TH1D *h1d = NULL;
  TH2F *h2 = NULL;
  TObject *obj = NULL;
  while( (obj = (TObject*)nexthist())!=NULL  ){
    if(obj->InheritsFrom("TH1F")){
      h1 = (TH1F*) obj;
      h1->GetXaxis()->CenterTitle();
      //h1->GetXaxis()->SetTitleSize(0.05);
      //h1->GetXaxis()->SetTitleOffset(0.80);
      h1->GetYaxis()->SetTitleOffset(1.5);
    }
    if(obj->InheritsFrom("TH1D")){
      h1d = (TH1D*) obj;
      h1d->GetXaxis()->CenterTitle();
      //h1d->GetXaxis()->SetTitleSize(0.05);
      //h1d->GetXaxis()->SetTitleOffset(0.80);
      h1d->GetYaxis()->SetTitleOffset(1.5);
    }
    if(obj->InheritsFrom("TH2")){
      h2 = (TH2F*) obj;
      h2->GetXaxis()->CenterTitle();
      h2->GetYaxis()->CenterTitle();
      //h2->GetXaxis()->SetTitleSize(0.05);
      //h2->GetXaxis()->SetTitleOffset(0.80);
      //h2->GetYaxis()->SetTitleSize(0.05);
      h2->GetYaxis()->SetTitleOffset(1.3);
    }
  }

  TString pdfname = std::string(file1->GetName());
  pdfname.Replace(std::string(file1->GetName()).size()-4,5,"pdf");
  TCanvas *c = NULL;
  TSeqCollection *SCol = gROOT->GetListOfCanvases();
  int size = SCol->GetSize();
  TIter next(SCol);
  for(int i=0;i<size;i++){
    c= (TCanvas*)next();
    c->Draw();
    c->cd();
    TPaveText *pt;
    if(Spmode || Smmode){
      pt = new TPaveText(.72,0.90,0.90,0.99,"NDC");    
    }else{
      pt = new TPaveText(.72,0.90,0.90,0.99,"NDC");    
    }
    if(Spmode){
      pt->SetFillColor(kAzure-4);
      pt->AddText("MC #Sigma+#pi- mode");
    }
    else if(Smmode){
      pt->SetFillColor(kAzure-4);
      pt->AddText("MC #Sigma-#pi+ mode"); 
    }else{
      pt->AddText("Real Data");
      pt->SetFillColor(kCyan-9);
    }
    pt->SetBorderSize(1);
    pt->Draw();

    c->Modified();
    c->Update();
    if(i==0) c->Print(pdfname+"(",Form("Title:%s",c->GetTitle()));
    else if(i==size-1)c->Print(pdfname+")",Form("Title:%s",c->GetTitle())); 
    else c->Print(pdfname,Form("Title:%s",c->GetTitle())); 
  }
  TFile *facc = NULL;
  if(Spmode) facc = new TFile("acc_Sp.root","RECREATE");
  if(Smmode) facc = new TFile("acc_Sm.root","RECREATE");
  heff->Write();
  heff_reco->Write();
  heff_woK0->Write();
  heff_woK0_reco->Write();
  heff_SpSm->Write();
  heff_woK0_SpSm->Write();
  heff_SpSm_reco->Write();
  heff_woK0_SpSm_reco->Write();
  heff_err->Write();


}
