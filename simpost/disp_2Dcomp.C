
void disp_2Dcomp(const char *filename="comp_fakedata_out.root")
{
 
  TFile *f = TFile::Open(filename,"READ");
  f->cd();
  gStyle->SetPadBorderMode(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  //const char opt[10]="cont3";
  const char opt[10]="cont1";
  //const char opt[10]="colz";
  gStyle->SetErrorX(0.);  
  TCanvas *cIMnpim_IMnpip_n = new TCanvas("cIMnpim_IMnpip_n", "cIMnpim_IMnpip_n",1200,1000);
  cIMnpim_IMnpip_n->Divide(2,2,0,0);
  TH2D* IMnpim_IMnpip_n_rdata = (TH2D*)f->Get("IMnpim_IMnpip_n_rdata");
  TH2D* IMnpim_IMnpip_n_mc    = (TH2D*)f->Get("IMnpim_IMnpip_n_mc");
  
  cIMnpim_IMnpip_n->cd(3);
  IMnpim_IMnpip_n_rdata->RebinX(4);
  IMnpim_IMnpip_n_rdata->RebinY(4);
  IMnpim_IMnpip_n_rdata->GetXaxis()->SetRangeUser(1,1.7);
  IMnpim_IMnpip_n_rdata->GetYaxis()->SetRangeUser(1,1.7);
  IMnpim_IMnpip_n_rdata->Draw(opt);
  
  cIMnpim_IMnpip_n->cd(1);
  TH1D* IMnpip_n_rdata = (TH1D*)IMnpim_IMnpip_n_rdata->ProjectionX("IMnpip_n_rdata");
  IMnpip_n_rdata->SetTitle("");
  IMnpip_n_rdata->GetXaxis()->SetTitle("");
  IMnpip_n_rdata->Draw("E1");
  IMnpip_n_rdata->GetXaxis()->SetLabelSize(0);
  
  cIMnpim_IMnpip_n->cd(4);
  TH1D* IMnpim_n_rdata = (TH1D*)IMnpim_IMnpip_n_rdata->ProjectionY("IMnpim_n_rdata");
  IMnpim_n_rdata->SetTitle("");
  IMnpim_n_rdata->GetXaxis()->SetTitle("");
  TGraphErrors *gr_IMnpim_n_rdata = new TGraphErrors();
  for(int ibin=0;ibin<IMnpim_n_rdata->GetNbinsX();ibin++){
    double cont = IMnpim_n_rdata->GetBinContent(ibin);
    double bincenter = IMnpim_n_rdata->GetBinCenter(ibin);
    double err = IMnpim_n_rdata->GetBinError(ibin);
    gr_IMnpim_n_rdata->SetPoint(ibin,cont,bincenter);
    gr_IMnpim_n_rdata->SetPointError(ibin,err,0);
  }
  gr_IMnpim_n_rdata->GetYaxis()->SetLabelOffset(999);
  gr_IMnpim_n_rdata->GetYaxis()->SetLabelSize(0);
  gr_IMnpim_n_rdata->GetYaxis()->SetRangeUser(IMnpim_IMnpip_n_rdata->GetYaxis()->GetXmin(),1.7);//IMnpim_IMnpip_n_rdata->GetYaxis()->GetXmax());
  gr_IMnpim_n_rdata->SetLineWidth(2);
  gr_IMnpim_n_rdata->Draw("AP");
  



  //draw BG at first
  TH2D* IMnpim_IMnpip_woK0_woSid_won_data = (TH2D*)f->Get("IMnpim_IMnpip_woK0_woSid_won_data");IMnpim_IMnpip_woK0_woSid_won_data->Print();
  TH2D* IMnpim_IMnpip_wK0_woSid_won_data = (TH2D*)f->Get("IMnpim_IMnpip_wK0_woSid_won_data");IMnpim_IMnpip_wK0_woSid_won_data->Print();

  TH2D* IMnpim_IMnpip_woK0_woSid_won_mc = (TH2D*)f->Get("IMnpim_IMnpip_woK0_woSid_won_mc");IMnpim_IMnpip_woK0_woSid_won_mc->Print();

  TH2D* IMnpim_IMnpip_wK0_woSid_won_mc = (TH2D*)f->Get("IMnpim_IMnpip_wK0_woSid_won_mc");IMnpim_IMnpip_wK0_woSid_won_mc->Print();
  //TH2D* IMnpim_IMnpip_wK0_woSid_won_mcgeta = (TH2D*)f->Get("IMnpim_IMnpip_wK0_woSid_won_mcgeta");IMnpim_IMnpip_wK0_woSid_won_mcgeta->Print();

  TH2D* IMnpim_IMnpip_woSid_won_data = (TH2D*)IMnpim_IMnpip_woK0_woSid_won_data->Clone("IMnpim_IMnpip_woSid_won_data");
  IMnpim_IMnpip_woSid_won_data->Add(IMnpim_IMnpip_wK0_woSid_won_data);
  TH2D* IMnpim_IMnpip_woSid_won_mc = (TH2D*)IMnpim_IMnpip_woK0_woSid_won_mc->Clone("IMnpim_IMnpip_woSid_won_mc");
  IMnpim_IMnpip_woSid_won_mc->Add(IMnpim_IMnpip_wK0_woSid_won_mc);
  //IMnpim_IMnpip_woSid_won_mc->Add(IMnpim_IMnpip_wK0_woSid_won_mcgeta);
   
  
  TCanvas *cIMnpim_IMnpip_woSid_won = new TCanvas("cIMnpim_IMnpip_woSid_won","cIMnpim_IMnpip_woSid_won",1200,1000);
  cIMnpim_IMnpip_woSid_won->Divide(2,2,0,0);
  cIMnpim_IMnpip_woSid_won->cd(3);
  IMnpim_IMnpip_woSid_won_data->SetTitle("");
  IMnpim_IMnpip_woSid_won_data->RebinX(4);
  IMnpim_IMnpip_woSid_won_data->RebinY(4);
  IMnpim_IMnpip_woSid_won_data->GetXaxis()->SetRangeUser(1,1.7);
  IMnpim_IMnpip_woSid_won_data->GetYaxis()->SetRangeUser(1,1.7);
  IMnpim_IMnpip_woSid_won_mc->RebinX(4);
  IMnpim_IMnpip_woSid_won_mc->RebinY(4);
  IMnpim_IMnpip_woSid_won_mc->GetXaxis()->SetRangeUser(1,1.7);
  IMnpim_IMnpip_woSid_won_mc->GetYaxis()->SetRangeUser(1,1.7);
  //IMnpim_IMnpip_woSid_won_data->GetZaxis()->SetNdivisions(503);
  //IMnpim_IMnpip_woSid_won_mc->GetZaxis()->SetNdivisions(503);
  IMnpim_IMnpip_woSid_won_data->SetContour(8);
  IMnpim_IMnpip_woSid_won_mc->SetContour(8);

  IMnpim_IMnpip_woSid_won_data->SetLineWidth(1);
  IMnpim_IMnpip_woSid_won_data->Draw("cont3");
  //IMnpim_IMnpip_woSid_won_data->Draw("col");
  IMnpim_IMnpip_woSid_won_mc->SetLineColor(2);
  IMnpim_IMnpip_woSid_won_mc->SetLineWidth(1);
  IMnpim_IMnpip_woSid_won_mc->Draw("cont3same");
  TPaletteAxis *palette2 = (TPaletteAxis*)IMnpim_IMnpip_woSid_won_data->GetListOfFunctions()->FindObject("palette");
  // the following lines moe the paletter. Choose the values you need for the position.
  palette2->SetX1NDC(0.78);
  palette2->SetX2NDC(0.83);
  palette2->SetY1NDC(0.1);
  palette2->SetY2NDC(0.90);
   
  
  cIMnpim_IMnpip_woSid_won->cd(1);
  TH1D* IMnpip_woSid_won_data = (TH1D*)IMnpim_IMnpip_woSid_won_data->ProjectionX("IMnpip_woSid_won_data");
  TH1D* IMnpip_woSid_won_mc = (TH1D*)IMnpim_IMnpip_woSid_won_mc->ProjectionX("IMnpip_woSid_won_mc");
  IMnpip_woSid_won_data->GetXaxis()->SetLabelSize(0);
  IMnpip_woSid_won_data->SetMarkerStyle(20);
  IMnpip_woSid_won_data->SetMarkerColor(1);
  IMnpip_woSid_won_data->Draw("E");
  IMnpip_woSid_won_mc->SetLineColor(2);
  IMnpip_woSid_won_mc->SetMarkerStyle(20);
  IMnpip_woSid_won_mc->SetMarkerColor(2);
  IMnpip_woSid_won_mc->Draw("Esame");
  
  
  cIMnpim_IMnpip_woSid_won->cd(4);
  TH1D* IMnpim_woSid_won_data = (TH1D*)IMnpim_IMnpip_woSid_won_data->ProjectionY("IMnpim_woSid_won_data");
  TH1D* IMnpim_woSid_won_mc = (TH1D*)IMnpim_IMnpip_woSid_won_mc->ProjectionY("IMnpim_woSid_won_mc");
  IMnpim_woSid_won_data->SetTitle("");
  IMnpim_woSid_won_data->GetXaxis()->SetTitle("");
  
  TGraphErrors *gr_IMnpim_woSid_won_data = new TGraphErrors();
  TGraphErrors *gr_IMnpim_woSid_won_mc = new TGraphErrors();
  
  for(int ibin=0;ibin<IMnpim_woSid_won_data->GetNbinsX();ibin++){
    double cont = IMnpim_woSid_won_data->GetBinContent(ibin);
    double err = IMnpim_woSid_won_data->GetBinError(ibin);
    double bincenter = IMnpim_woSid_won_data->GetBinCenter(ibin);
    gr_IMnpim_woSid_won_data->SetPoint(ibin,cont, bincenter);
    gr_IMnpim_woSid_won_data->SetPointError(ibin,err,0);
  }
  for(int ibin=0;ibin<IMnpim_woSid_won_mc->GetNbinsX();ibin++){
    double cont = IMnpim_woSid_won_mc->GetBinContent(ibin);
    double err = IMnpim_woSid_won_mc->GetBinError(ibin);
    double bincenter = IMnpim_woSid_won_mc->GetBinCenter(ibin);
    gr_IMnpim_woSid_won_mc->SetPoint(ibin,cont, bincenter);
    gr_IMnpim_woSid_won_mc->SetPointError(ibin,err,0);
  }
  gr_IMnpim_woSid_won_data->GetYaxis()->SetRangeUser(1,1.7);
  gr_IMnpim_woSid_won_data->SetMarkerStyle(20);
  gr_IMnpim_woSid_won_data->GetYaxis()->SetLabelSize(0);
  gr_IMnpim_woSid_won_data->Draw("AP");
  gr_IMnpim_woSid_won_mc->SetMarkerStyle(20);
  gr_IMnpim_woSid_won_mc->SetMarkerColor(2);
  gr_IMnpim_woSid_won_mc->Draw("P");
  
  TCanvas *cIMnpim_IMnpip_woSid_won_2 = new TCanvas("cIMnpim_IMnpip_woSid_won_2","cIMnpim_IMnpip_woSid_won_2",1200,1000);
  cIMnpim_IMnpip_woSid_won_2->Divide(2,1);
  cIMnpim_IMnpip_woSid_won_2->cd(1);
  TH2D* IMnpim_IMnpip_woSid_won_data_2 = (TH2D*)IMnpim_IMnpip_woSid_won_data->Clone("IMnpim_IMnpip_woSid_won_data_2");
  IMnpim_IMnpip_woSid_won_data_2->SetContour(10);
  IMnpim_IMnpip_woSid_won_data_2->Draw("cont1z");
  gPad->SetRightMargin(0);
  cIMnpim_IMnpip_woSid_won_2->cd(2);
  TH2D* IMnpim_IMnpip_woSid_won_mc_2 = (TH2D*)IMnpim_IMnpip_woSid_won_mc->Clone("IMnpim_IMnpip_woSid_won_mc_2");
  IMnpim_IMnpip_woSid_won_mc->SetTitle("");
  IMnpim_IMnpip_woSid_won_mc->SetContour(10);
  IMnpim_IMnpip_woSid_won_mc->SetMaximum(IMnpim_IMnpip_woSid_won_data->GetMaximum());
  IMnpim_IMnpip_woSid_won_mc->Draw("cont1z");
  gPad->SetLeftMargin(0);
  //TPaletteAxis *palette = (TPaletteAxis*)IMnpim_IMnpip_woSid_won_mc->GetListOfFunctions()->FindObject("palette");
  // the following lines moe the paletter. Choose the values you need for the position.
  //palette->SetX1NDC(0.8);
  //palette->SetX2NDC(0.85);
  //palette->SetY1NDC(0.1);
  //palette->SetY2NDC(0.90);
  //cIMnpim_IMnpip_woSid_won_2->Modified();
  //cIMnpim_IMnpip_woSid_won_2->Update();

  TCanvas *cIMnpim_IMnpip_woSid_won_sub = new TCanvas("cIMnpim_IMnpip_woSid_won_sub","cIMnpim_IMnpip_woSid_won_sub");
  TH2D* IMnpim_IMnpip_woSid_won_sub = (TH2D*)IMnpim_IMnpip_woSid_won_data->Clone("IMnpim_IMnpip_woSid_won_data");
  IMnpim_IMnpip_woSid_won_sub->Add(IMnpim_IMnpip_woSid_won_mc,-1.0);
  IMnpim_IMnpip_woSid_won_sub->Divide(IMnpim_IMnpip_woSid_won_data);
  IMnpim_IMnpip_woSid_won_sub->SetTitle("(data-MC)/data");
  //IMnpim_IMnpip_woSid_won_sub->RebinX(2);
  //IMnpim_IMnpip_woSid_won_sub->RebinY(2);
  //IMnpim_IMnpip_woSid_won_sub->Scale(1./4.);
  IMnpim_IMnpip_woSid_won_sub->GetZaxis()->SetRangeUser(-0.3,0.3);
  IMnpim_IMnpip_woSid_won_sub->Draw("colz");
  
  //
  //BG2 MMnmiss vs IMpipim
  //
  TH2D* MMnmiss_IMpippim_woK0_woSid_won_data = (TH2D*)f->Get("MMnmiss_IMpippim_woK0_woSid_won_data");MMnmiss_IMpippim_woK0_woSid_won_data->Print();
  TH2D* MMnmiss_IMpippim_wK0_woSid_won_data = (TH2D*)f->Get("MMnmiss_IMpippim_wK0_woSid_won_data");MMnmiss_IMpippim_wK0_woSid_won_data->Print();

  TH2D* MMnmiss_IMpippim_woK0_woSid_won_mc = (TH2D*)f->Get("MMnmiss_IMpippim_woK0_woSid_won_mc");MMnmiss_IMpippim_woK0_woSid_won_mc->Print();

  TH2D* MMnmiss_IMpippim_wK0_woSid_won_mc = (TH2D*)f->Get("MMnmiss_IMpippim_wK0_woSid_won_mc");MMnmiss_IMpippim_wK0_woSid_won_mc->Print();
  //TH2D* MMnmiss_IMpippim_wK0_woSid_won_mcgeta = (TH2D*)f->Get("MMnmiss_IMpippim_wK0_woSid_won_mcgeta");MMnmiss_IMpippim_wK0_woSid_won_mcgeta->Print();

  TH2D* MMnmiss_IMpippim_woSid_won_data = (TH2D*)MMnmiss_IMpippim_woK0_woSid_won_data->Clone("MMnmiss_IMpippim_woSid_won_data");
  MMnmiss_IMpippim_woSid_won_data->Add(MMnmiss_IMpippim_wK0_woSid_won_data);
  TH2D* MMnmiss_IMpippim_woSid_won_mc = (TH2D*)MMnmiss_IMpippim_woK0_woSid_won_mc->Clone("MMnmiss_IMpippim_woSid_won_mc");
  MMnmiss_IMpippim_woSid_won_mc->Add(MMnmiss_IMpippim_wK0_woSid_won_mc);
  //MMnmiss_IMpippim_woSid_won_mc->Add(MMnmiss_IMpippim_wK0_woSid_won_mcgeta);
   
  
  TCanvas *cMMnmiss_IMpippim_woSid_won = new TCanvas("cMMnmiss_IMpippim_woSid_won","cMMnmiss_IMpippim_woSid_won",1200,1000);
  cMMnmiss_IMpippim_woSid_won ->Divide(2,2,0,0);
  
 
  cMMnmiss_IMpippim_woSid_won->cd(3);
  MMnmiss_IMpippim_woSid_won_data->SetTitle("");
  MMnmiss_IMpippim_woSid_won_data->RebinX(4);
  MMnmiss_IMpippim_woSid_won_data->RebinY(3);
  MMnmiss_IMpippim_woSid_won_data->GetXaxis()->SetRangeUser(0.2,1.0);
  //MMnmiss_IMpippim_woSid_won_data->GetYaxis()->SetRangeUser(1,1.7);
  MMnmiss_IMpippim_woSid_won_mc->RebinX(4);
  MMnmiss_IMpippim_woSid_won_mc->RebinY(3);
  MMnmiss_IMpippim_woSid_won_mc->GetXaxis()->SetRangeUser(0.2,1.0);
  //MMnmiss_IMpippim_woSid_won_mc->GetYaxis()->SetRangeUser(1,1.7);
  MMnmiss_IMpippim_woSid_won_data->SetContour(10);
  MMnmiss_IMpippim_woSid_won_mc->SetContour(10);

  MMnmiss_IMpippim_woSid_won_data->SetLineWidth(1);
  MMnmiss_IMpippim_woSid_won_data->Draw("cont3");
  MMnmiss_IMpippim_woSid_won_mc->SetLineColor(2);
  MMnmiss_IMpippim_woSid_won_mc->SetLineWidth(1);
  MMnmiss_IMpippim_woSid_won_mc->Draw("cont3same");
   
  
  cMMnmiss_IMpippim_woSid_won->cd(1);
  TH1D* IMpippim_woSid_won_data = (TH1D*)MMnmiss_IMpippim_woSid_won_data->ProjectionX("IMpippim_woSid_won_data");
  TH1D* IMpippim_woSid_won_mc = (TH1D*)MMnmiss_IMpippim_woSid_won_mc->ProjectionX("IMpippim_woSid_won_mc");
  IMpippim_woSid_won_data->SetTitle("");
  IMpippim_woSid_won_data->GetXaxis()->SetLabelSize(0);
  IMpippim_woSid_won_data->SetMarkerStyle(20);
  IMpippim_woSid_won_data->SetMarkerColor(1);
  IMpippim_woSid_won_data->Draw("E");
  IMpippim_woSid_won_mc->SetLineColor(2);
  IMpippim_woSid_won_mc->SetMarkerStyle(20);
  IMpippim_woSid_won_mc->SetMarkerColor(2);
  IMpippim_woSid_won_mc->Draw("Esame");
  
  
  cMMnmiss_IMpippim_woSid_won->cd(4);
  TH1D* MMnmiss_woSid_won_data = (TH1D*)MMnmiss_IMpippim_woSid_won_data->ProjectionY("MMnmiss_woSid_won_data");
  TH1D* MMnmiss_woSid_won_mc = (TH1D*)MMnmiss_IMpippim_woSid_won_mc->ProjectionY("MMnmiss_woSid_won_mc");
  MMnmiss_woSid_won_data->SetTitle("");
  MMnmiss_woSid_won_data->GetXaxis()->SetTitle("");
  
  TGraphErrors *gr_MMnmiss_woSid_won_data = new TGraphErrors();
  TGraphErrors *gr_MMnmiss_woSid_won_mc = new TGraphErrors();
  
  for(int ibin=0;ibin<MMnmiss_woSid_won_data->GetNbinsX();ibin++){
    double cont = MMnmiss_woSid_won_data->GetBinContent(ibin);
    double err = MMnmiss_woSid_won_data->GetBinError(ibin);
    double bincenter = MMnmiss_woSid_won_data->GetBinCenter(ibin);
    gr_MMnmiss_woSid_won_data->SetPoint(ibin,cont, bincenter);
    gr_MMnmiss_woSid_won_data->SetPointError(ibin,err,0);
  }
  for(int ibin=0;ibin<MMnmiss_woSid_won_mc->GetNbinsX();ibin++){
    double cont = MMnmiss_woSid_won_mc->GetBinContent(ibin);
    double err = MMnmiss_woSid_won_mc->GetBinError(ibin);
    double bincenter = MMnmiss_woSid_won_mc->GetBinCenter(ibin);
    gr_MMnmiss_woSid_won_mc->SetPoint(ibin,cont, bincenter);
    gr_MMnmiss_woSid_won_mc->SetPointError(ibin,err,0);
  }
  gr_MMnmiss_woSid_won_data->GetYaxis()->SetRangeUser(0,1.5);
  gr_MMnmiss_woSid_won_data->SetMarkerStyle(20);
  gr_MMnmiss_woSid_won_data->GetYaxis()->SetLabelSize(0);
  gr_MMnmiss_woSid_won_data->Draw("AP");
  gr_MMnmiss_woSid_won_mc->SetMarkerStyle(20);
  gr_MMnmiss_woSid_won_mc->SetMarkerColor(2);
  gr_MMnmiss_woSid_won_mc->Draw("P");

  TCanvas *cMMnmiss_IMpippim_woSid_won_2 = new TCanvas("cMMnmiss_IMpippim_woSid_won_2","cMMnmiss_IMpippim_woSid_won_2",1200,1000);
  cMMnmiss_IMpippim_woSid_won_2->Divide(2,1);
  cMMnmiss_IMpippim_woSid_won_2->cd(1);
  TH2D* MMnmiss_IMpippim_woSid_won_data_2 = (TH2D*)MMnmiss_IMpippim_woSid_won_data->Clone("MMnmiss_IMpippim_woSid_won_data_2");
  TH2D* MMnmiss_IMpippim_woSid_won_mc_2 = (TH2D*)MMnmiss_IMpippim_woSid_won_mc->Clone("MMnmiss_IMpippim_woSid_won_mc_2");
  MMnmiss_IMpippim_woSid_won_data_2->SetTitle("real data");
  MMnmiss_IMpippim_woSid_won_data_2->SetContour(10); 
  MMnmiss_IMpippim_woSid_won_data_2->Draw("cont1z");
  gPad->SetRightMargin(0);
  cMMnmiss_IMpippim_woSid_won_2->cd(2);
  MMnmiss_IMpippim_woSid_won_mc_2->SetTitle("MC");
  MMnmiss_IMpippim_woSid_won_mc_2->SetContour(10); 
  MMnmiss_IMpippim_woSid_won_mc_2->SetMaximum(MMnmiss_IMpippim_woSid_won_data_2->GetMaximum());
  MMnmiss_IMpippim_woSid_won_mc_2->Draw("cont1z");
  gPad->SetLeftMargin(0);
  

  TCanvas *cMMnmiss_IMpippim_woSid_won_sub = new TCanvas("cMMnmiss_IMpippim_woSid_won_sub","cMMnmiss_IMpippim_woSid_won_sub");
  TH2D* MMnmiss_IMpippim_woSid_won_sub = (TH2D*)MMnmiss_IMpippim_woSid_won_data->Clone("MMnmiss_IMpippim_woSid_won_data");
  MMnmiss_IMpippim_woSid_won_sub->Add(MMnmiss_IMpippim_woSid_won_mc,-1.0);
  MMnmiss_IMpippim_woSid_won_sub->Divide(MMnmiss_IMpippim_woSid_won_data);
  MMnmiss_IMpippim_woSid_won_sub->SetTitle("(data-MC)/data");
  MMnmiss_IMpippim_woSid_won_sub->RebinX(2);
  MMnmiss_IMpippim_woSid_won_sub->RebinY(2);
  MMnmiss_IMpippim_woSid_won_sub->Scale(0.25);
  MMnmiss_IMpippim_woSid_won_sub->GetZaxis()->SetRangeUser(-0.3,0.3);
  MMnmiss_IMpippim_woSid_won_sub->Draw("colz");
 

  //
  //BG3 q vs nmom
  //
  TH2D* q_nmom_woK0_woSid_won_data = (TH2D*)f->Get("q_nmom_woK0_woSid_won_data");q_nmom_woK0_woSid_won_data->Print();
  TH2D* q_nmom_wK0_woSid_won_data = (TH2D*)f->Get("q_nmom_wK0_woSid_won_data");q_nmom_wK0_woSid_won_data->Print();

  TH2D* q_nmom_woK0_woSid_won_mc = (TH2D*)f->Get("q_nmom_woK0_woSid_won_mc");q_nmom_woK0_woSid_won_mc->Print();

  TH2D* q_nmom_wK0_woSid_won_mc = (TH2D*)f->Get("q_nmom_wK0_woSid_won_mc");q_nmom_wK0_woSid_won_mc->Print();

  TH2D* q_nmom_woSid_won_data = (TH2D*)q_nmom_woK0_woSid_won_data->Clone("q_nmom_woSid_won_data");
  q_nmom_woSid_won_data->Add(q_nmom_wK0_woSid_won_data);
  TH2D* q_nmom_woSid_won_mc = (TH2D*)q_nmom_woK0_woSid_won_mc->Clone("q_nmom_woSid_won_mc");
  q_nmom_woSid_won_mc->Add(q_nmom_wK0_woSid_won_mc);
  //q_nmom_woSid_won_mc->Add(q_nmom_wK0_woSid_won_mcgeta);
   
  
  TCanvas *cq_nmom_woSid_won = new TCanvas("cq_nmom_woSid_won","cq_nmom_woSid_won",1200,1000);
  cq_nmom_woSid_won ->Divide(2,2,0,0);
  
 
  cq_nmom_woSid_won->cd(3);
  q_nmom_woSid_won_data->SetTitle("");
  q_nmom_woSid_won_data->RebinX(4);
  q_nmom_woSid_won_data->RebinY(4);
  q_nmom_woSid_won_data->GetXaxis()->SetRangeUser(0.2,1.0);
  //q_nmom_woSid_won_data->GetYaxis()->SetRangeUser(1,1.7);
  q_nmom_woSid_won_mc->RebinX(4);
  q_nmom_woSid_won_mc->RebinY(4);
  q_nmom_woSid_won_mc->GetXaxis()->SetRangeUser(0.2,1.0);
  //q_nmom_woSid_won_mc->GetYaxis()->SetRangeUser(1,1.7);
  q_nmom_woSid_won_data->SetContour(10);
  q_nmom_woSid_won_mc->SetContour(10);

  q_nmom_woSid_won_data->SetLineWidth(1);
  q_nmom_woSid_won_data->Draw("cont3");
  q_nmom_woSid_won_mc->SetLineColor(2);
  q_nmom_woSid_won_mc->SetLineWidth(1);
  q_nmom_woSid_won_mc->Draw("cont3same");
   
  
  cq_nmom_woSid_won->cd(1);
  TH1D* nmom_woSid_won_data = (TH1D*)q_nmom_woSid_won_data->ProjectionX("nmom_woSid_won_data");
  TH1D* nmom_woSid_won_mc = (TH1D*)q_nmom_woSid_won_mc->ProjectionX("nmom_woSid_won_mc");
  nmom_woSid_won_data->SetTitle("");
  nmom_woSid_won_data->GetXaxis()->SetLabelSize(0);
  nmom_woSid_won_data->SetMarkerStyle(20);
  nmom_woSid_won_data->SetMarkerColor(1);
  nmom_woSid_won_data->Draw("E");
  nmom_woSid_won_mc->SetLineColor(2);
  nmom_woSid_won_mc->SetMarkerStyle(20);
  nmom_woSid_won_mc->SetMarkerColor(2);
  nmom_woSid_won_mc->Draw("Esame");
  
  
  cq_nmom_woSid_won->cd(4);
  TH1D* q_woSid_won_data = (TH1D*)q_nmom_woSid_won_data->ProjectionY("q_woSid_won_data");
  TH1D* q_woSid_won_mc = (TH1D*)q_nmom_woSid_won_mc->ProjectionY("q_woSid_won_mc");
  q_woSid_won_data->SetTitle("");
  q_woSid_won_data->GetXaxis()->SetTitle("");
  
  TGraphErrors *gr_q_woSid_won_data = new TGraphErrors();
  TGraphErrors *gr_q_woSid_won_mc = new TGraphErrors();
  
  for(int ibin=0;ibin<q_woSid_won_data->GetNbinsX();ibin++){
    double cont = q_woSid_won_data->GetBinContent(ibin);
    double err = q_woSid_won_data->GetBinError(ibin);
    double bincenter = q_woSid_won_data->GetBinCenter(ibin);
    gr_q_woSid_won_data->SetPoint(ibin,cont, bincenter);
    gr_q_woSid_won_data->SetPointError(ibin,err,0);
  }
  for(int ibin=0;ibin<q_woSid_won_mc->GetNbinsX();ibin++){
    double cont = q_woSid_won_mc->GetBinContent(ibin);
    double err = q_woSid_won_mc->GetBinError(ibin);
    double bincenter = q_woSid_won_mc->GetBinCenter(ibin);
    gr_q_woSid_won_mc->SetPoint(ibin,cont, bincenter);
    gr_q_woSid_won_mc->SetPointError(ibin,err,0);
  }
  gr_q_woSid_won_data->GetYaxis()->SetRangeUser(0,1.5);
  gr_q_woSid_won_data->SetMarkerStyle(20);
  gr_q_woSid_won_data->GetYaxis()->SetLabelSize(0);
  gr_q_woSid_won_data->Draw("AP");
  gr_q_woSid_won_mc->SetMarkerStyle(20);
  gr_q_woSid_won_mc->SetMarkerColor(2);
  gr_q_woSid_won_mc->Draw("P");

  TCanvas *cq_nmom_woSid_won_2 = new TCanvas("cq_nmom_woSid_won_2","cq_nmom_woSid_won_2",1200,1000);
  cq_nmom_woSid_won_2->Divide(2,1);
  cq_nmom_woSid_won_2->cd(1);
  TH2D* q_nmom_woSid_won_data_2 = (TH2D*)q_nmom_woSid_won_data->Clone("q_nmom_woSid_won_data_2");
  TH2D* q_nmom_woSid_won_mc_2 = (TH2D*)q_nmom_woSid_won_mc->Clone("q_nmom_woSid_won_mc_2");
  q_nmom_woSid_won_data_2->SetTitle("real data");
  q_nmom_woSid_won_data_2->SetContour(10); 
  q_nmom_woSid_won_data_2->Draw("cont1z");
  gPad->SetRightMargin(0);
  cq_nmom_woSid_won_2->cd(2);
  q_nmom_woSid_won_mc_2->SetTitle("MC");
  q_nmom_woSid_won_mc_2->SetContour(10); 
  q_nmom_woSid_won_mc_2->SetMaximum(q_nmom_woSid_won_data_2->GetMaximum());
  q_nmom_woSid_won_mc_2->Draw("cont1z");
  gPad->SetLeftMargin(0);
  
  TCanvas *cq_nmom_woSid_won_sub = new TCanvas("cq_nmom_woSid_won_sub","cq_nmom_woSid_won_sub");
  TH2D* q_nmom_woSid_won_sub = (TH2D*)q_nmom_woSid_won_data->Clone("q_nmom_woSid_won_data");
  q_nmom_woSid_won_sub->Add(q_nmom_woSid_won_mc,-1.0);
  q_nmom_woSid_won_sub->Divide(q_nmom_woSid_won_data);
  q_nmom_woSid_won_sub->SetTitle("(data-MC)/data");
  q_nmom_woSid_won_sub->RebinX(2);
  q_nmom_woSid_won_sub->RebinY(2);
  q_nmom_woSid_won_sub->Scale(0.25);
  q_nmom_woSid_won_sub->GetZaxis()->SetRangeUser(-0.3,0.3);
  q_nmom_woSid_won_sub->Draw("colz");
  
  //
  //BG check q vs IMnpipi
  //
  TH2D* q_IMnpipi_woK0_woSid_won_data = (TH2D*)f->Get("q_IMnpipi_woK0_woSid_won_data");q_IMnpipi_woK0_woSid_won_data->Print();
  TH2D* q_IMnpipi_wK0_woSid_won_data = (TH2D*)f->Get("q_IMnpipi_wK0_woSid_won_data");q_IMnpipi_wK0_woSid_won_data->Print();

  TH2D* q_IMnpipi_woK0_woSid_won_mc = (TH2D*)f->Get("q_IMnpipi_woK0_woSid_won_mc");q_IMnpipi_woK0_woSid_won_mc->Print();

  TH2D* q_IMnpipi_wK0_woSid_won_mc = (TH2D*)f->Get("q_IMnpipi_wK0_woSid_won_mc");q_IMnpipi_wK0_woSid_won_mc->Print();

  TH2D* q_IMnpipi_woSid_won_data = (TH2D*)q_IMnpipi_woK0_woSid_won_data->Clone("q_IMnpipi_woSid_won_data");
  q_IMnpipi_woSid_won_data->Add(q_IMnpipi_wK0_woSid_won_data);
  TH2D* q_IMnpipi_woSid_won_mc = (TH2D*)q_IMnpipi_woK0_woSid_won_mc->Clone("q_IMnpipi_woSid_won_mc");
  q_IMnpipi_woSid_won_mc->Add(q_IMnpipi_wK0_woSid_won_mc);
  //q_IMnpipi_woSid_won_mc->Add(q_IMnpipi_wK0_woSid_won_mcgeta);
   
  
  TCanvas *cq_IMnpipi_woSid_won = new TCanvas("cq_IMnpipi_woSid_won","cq_IMnpipi_woSid_won",1200,1000);
  cq_IMnpipi_woSid_won ->Divide(2,2,0,0);
  
 
  cq_IMnpipi_woSid_won->cd(3);
  q_IMnpipi_woSid_won_data->SetTitle("");
  q_IMnpipi_woSid_won_data->RebinX(2);
  std::cout << "bin width x q_IMnpipi_woSid_won_data " << q_IMnpipi_woSid_won_data->GetXaxis()->GetBinWidth(1) << std::endl;
  q_IMnpipi_woSid_won_data->RebinY(4);
  std::cout << "bin width y q_IMnpipi_woSid_won_data " << q_IMnpipi_woSid_won_data->GetYaxis()->GetBinWidth(1) << std::endl;
  q_IMnpipi_woSid_won_data->GetXaxis()->SetRangeUser(0.2,1.0);
  //q_IMnpipi_woSid_won_data->GetYaxis()->SetRangeUser(1,1.7);
  q_IMnpipi_woSid_won_mc->RebinX(2);
  q_IMnpipi_woSid_won_mc->RebinY(4);
  q_IMnpipi_woSid_won_mc->GetXaxis()->SetRangeUser(0.2,1.0);
  //q_IMnpipi_woSid_won_mc->GetYaxis()->SetRangeUser(1,1.7);
  q_IMnpipi_woSid_won_data->SetContour(10);
  q_IMnpipi_woSid_won_mc->SetContour(10);
  q_IMnpipi_woSid_won_data->SetLineWidth(1);
  q_IMnpipi_woSid_won_data->Draw("cont3");
  q_IMnpipi_woSid_won_mc->SetLineColor(2);
  q_IMnpipi_woSid_won_mc->SetLineWidth(1);
  q_IMnpipi_woSid_won_mc->Draw("cont3same");
   
  
  cq_IMnpipi_woSid_won->cd(1);
  TH1D* IMnpipi_woSid_won_data = (TH1D*)q_IMnpipi_woSid_won_data->ProjectionX("IMnpipi_woSid_won_data");
  TH1D* IMnpipi_woSid_won_mc = (TH1D*)q_IMnpipi_woSid_won_mc->ProjectionX("IMnpipi_woSid_won_mc");
  IMnpipi_woSid_won_data->SetTitle("");
  IMnpipi_woSid_won_data->GetXaxis()->SetLabelSize(0);
  IMnpipi_woSid_won_data->SetMarkerStyle(20);
  IMnpipi_woSid_won_data->SetMarkerColor(1);
  IMnpipi_woSid_won_data->Draw("E");
  IMnpipi_woSid_won_mc->SetLineColor(2);
  IMnpipi_woSid_won_mc->SetMarkerStyle(20);
  IMnpipi_woSid_won_mc->SetMarkerColor(2);
  IMnpipi_woSid_won_mc->Draw("Esame");
  
  
  cq_IMnpipi_woSid_won->cd(4);
  TH1D* q_woSid_won_data_2 = (TH1D*)q_IMnpipi_woSid_won_data->ProjectionY("q_woSid_won_data");
  TH1D* q_woSid_won_mc_2 = (TH1D*)q_IMnpipi_woSid_won_mc->ProjectionY("q_woSid_won_mc");
  q_woSid_won_data_2->SetTitle("");
  q_woSid_won_data_2->GetXaxis()->SetTitle("");
  
  TGraphErrors *gr_q_woSid_won_data_2 = new TGraphErrors();
  TGraphErrors *gr_q_woSid_won_mc_2 = new TGraphErrors();
  
  for(int ibin=0;ibin<q_woSid_won_data_2->GetNbinsX();ibin++){
    double cont = q_woSid_won_data_2->GetBinContent(ibin);
    double err = q_woSid_won_data_2->GetBinError(ibin);
    double bincenter = q_woSid_won_data_2->GetBinCenter(ibin);
    gr_q_woSid_won_data_2->SetPoint(ibin,cont, bincenter);
    gr_q_woSid_won_data_2->SetPointError(ibin,err,0);
  }
  for(int ibin=0;ibin<q_woSid_won_mc_2->GetNbinsX();ibin++){
    double cont = q_woSid_won_mc_2->GetBinContent(ibin);
    double err = q_woSid_won_mc_2->GetBinError(ibin);
    double bincenter = q_woSid_won_mc_2->GetBinCenter(ibin);
    gr_q_woSid_won_mc_2->SetPoint(ibin,cont, bincenter);
    gr_q_woSid_won_mc_2->SetPointError(ibin,err,0);
  }
  gr_q_woSid_won_data_2->GetYaxis()->SetRangeUser(0,1.5);
  gr_q_woSid_won_data_2->SetMarkerStyle(20);
  gr_q_woSid_won_data_2->GetYaxis()->SetLabelSize(0);
  gr_q_woSid_won_data_2->Draw("AP");
  gr_q_woSid_won_mc_2->SetMarkerStyle(20);
  gr_q_woSid_won_mc_2->SetMarkerColor(2);
  gr_q_woSid_won_mc_2->Draw("P");

  TCanvas *cq_IMnpipi_woSid_won_2 = new TCanvas("cq_IMnpipi_woSid_won_2","cq_IMnpipi_woSid_won_2",1200,1000);
  cq_IMnpipi_woSid_won_2->Divide(2,1);
  cq_IMnpipi_woSid_won_2->cd(1);
  TH2D* q_IMnpipi_woSid_won_data_2 = (TH2D*)q_IMnpipi_woSid_won_data->Clone("q_IMnpipi_woSid_won_data_2");
  TH2D* q_IMnpipi_woSid_won_mc_2 = (TH2D*)q_IMnpipi_woSid_won_mc->Clone("q_IMnpipi_woSid_won_mc_2");
  q_IMnpipi_woSid_won_data_2->SetTitle("real data");
  q_IMnpipi_woSid_won_data_2->SetContour(10);
  q_IMnpipi_woSid_won_data_2->SetMinimum(1);
  q_IMnpipi_woSid_won_data_2->Draw("cont1z");
  gPad->SetRightMargin(0);
  cq_IMnpipi_woSid_won_2->cd(2);
  q_IMnpipi_woSid_won_mc_2->SetTitle("MC");
  q_IMnpipi_woSid_won_mc_2->SetContour(10); 
  q_IMnpipi_woSid_won_mc_2->SetMinimum(1);
  q_IMnpipi_woSid_won_mc_2->SetMaximum(q_IMnpipi_woSid_won_data_2->GetMaximum());
  q_IMnpipi_woSid_won_mc_2->Draw("cont1z");
  gPad->SetLeftMargin(0);
  
  TCanvas *cq_IMnpipi_woSid_won_sub = new TCanvas("cq_IMnpipi_woSid_won_sub","cq_IMnpipi_woSid_won_sub");
  TH2D* q_IMnpipi_woSid_won_sub = (TH2D*)q_IMnpipi_woSid_won_data->Clone("q_IMnpipi_woSid_won_data");
  q_IMnpipi_woSid_won_sub->Add(q_IMnpipi_woSid_won_mc,-1.0);
  q_IMnpipi_woSid_won_sub->Divide(q_IMnpipi_woSid_won_data);
  q_IMnpipi_woSid_won_sub->SetTitle("(data-MC)/data");
  q_IMnpipi_woSid_won_sub->RebinX(2);
  q_IMnpipi_woSid_won_sub->RebinY(2);
  q_IMnpipi_woSid_won_sub->Scale(0.25);
  q_IMnpipi_woSid_won_sub->GetZaxis()->SetRangeUser(-0.3,0.3);
  q_IMnpipi_woSid_won_sub->Draw("colz");

  
  //
  //signal check q vs IMnpipi
  //
  TH2D* q_IMnpipi_woK0_wSid_n_data = (TH2D*)f->Get("q_IMnpipi_woK0_wSid_n_data");q_IMnpipi_woK0_wSid_n_data->Print();
  TH2D* q_IMnpipi_wK0_wSid_n_data = (TH2D*)f->Get("q_IMnpipi_wK0_wSid_n_data");q_IMnpipi_wK0_wSid_n_data->Print();

  TH2D* q_IMnpipi_woK0_wSid_n_mc = (TH2D*)f->Get("q_IMnpipi_woK0_wSid_n_mc");q_IMnpipi_woK0_wSid_n_mc->Print();

  TH2D* q_IMnpipi_wK0_wSid_n_mc = (TH2D*)f->Get("q_IMnpipi_wK0_wSid_n_mc");q_IMnpipi_wK0_wSid_n_mc->Print();

  TH2D* q_IMnpipi_wSid_n_data = (TH2D*)q_IMnpipi_woK0_wSid_n_data->Clone("q_IMnpipi_wSid_n_data");
  q_IMnpipi_wSid_n_data->Add(q_IMnpipi_wK0_wSid_n_data);
  TH2D* q_IMnpipi_wSid_n_mc = (TH2D*)q_IMnpipi_woK0_wSid_n_mc->Clone("q_IMnpipi_wSid_n_mc");
  q_IMnpipi_wSid_n_mc->Add(q_IMnpipi_wK0_wSid_n_mc);
   
  
  TCanvas *cq_IMnpipi_wSid_n = new TCanvas("cq_IMnpipi_wSid_n","cq_IMnpipi_wSid_n",1200,1000);
  cq_IMnpipi_wSid_n ->Divide(2,2,0,0);
  
 
  cq_IMnpipi_wSid_n->cd(3);
  q_IMnpipi_wSid_n_data->SetTitle("");
  q_IMnpipi_wSid_n_data->RebinX(2);
  std::cout << "bin width x q_IMnpipi_wSid_n_data " << q_IMnpipi_wSid_n_data->GetXaxis()->GetBinWidth(1) << std::endl;
  q_IMnpipi_wSid_n_data->RebinY(4);
  std::cout << "bin width y q_IMnpipi_wSid_n_data " << q_IMnpipi_wSid_n_data->GetYaxis()->GetBinWidth(1) << std::endl;
  q_IMnpipi_wSid_n_data->GetXaxis()->SetRangeUser(0.2,1.0);
  //q_IMnpipi_wSid_n_data->GetYaxis()->SetRangeUser(1,1.7);
  q_IMnpipi_wSid_n_mc->RebinX(2);
  q_IMnpipi_wSid_n_mc->RebinY(4);
  q_IMnpipi_wSid_n_mc->GetXaxis()->SetRangeUser(0.2,1.0);
  //q_IMnpipi_wSid_n_mc->GetYaxis()->SetRangeUser(1,1.7);
  q_IMnpipi_wSid_n_data->SetContour(10);
  q_IMnpipi_wSid_n_mc->SetContour(10);
  q_IMnpipi_wSid_n_data->SetLineWidth(1);
  q_IMnpipi_wSid_n_data->Draw("cont3");
  q_IMnpipi_wSid_n_mc->SetLineColor(2);
  q_IMnpipi_wSid_n_mc->SetLineWidth(1);
  q_IMnpipi_wSid_n_mc->Draw("cont3same");
   
  
  cq_IMnpipi_wSid_n->cd(1);
  TH1D* IMnpipi_wSid_n_data = (TH1D*)q_IMnpipi_wSid_n_data->ProjectionX("IMnpipi_wSid_n_data");
  TH1D* IMnpipi_wSid_n_mc = (TH1D*)q_IMnpipi_wSid_n_mc->ProjectionX("IMnpipi_wSid_n_mc");
  IMnpipi_wSid_n_data->SetTitle("");
  IMnpipi_wSid_n_data->GetXaxis()->SetLabelSize(0);
  IMnpipi_wSid_n_data->SetMarkerStyle(20);
  IMnpipi_wSid_n_data->SetMarkerColor(1);
  IMnpipi_wSid_n_data->Draw("E");
  IMnpipi_wSid_n_mc->SetLineColor(2);
  IMnpipi_wSid_n_mc->SetMarkerStyle(20);
  IMnpipi_wSid_n_mc->SetMarkerColor(2);
  IMnpipi_wSid_n_mc->Draw("Esame");
  
  
  cq_IMnpipi_wSid_n->cd(4);
  TH1D* q_wSid_n_data = (TH1D*)q_IMnpipi_wSid_n_data->ProjectionY("q_wSid_n_data");
  TH1D* q_wSid_n_mc = (TH1D*)q_IMnpipi_wSid_n_mc->ProjectionY("q_wSid_n_mc");
  q_wSid_n_data->SetTitle("");
  q_wSid_n_data->GetXaxis()->SetTitle("");
  
  TGraphErrors *gr_q_wSid_n_data = new TGraphErrors();
  TGraphErrors *gr_q_wSid_n_mc = new TGraphErrors();
  
  for(int ibin=0;ibin<q_wSid_n_data->GetNbinsX();ibin++){
    double cont = q_wSid_n_data->GetBinContent(ibin);
    double err = q_wSid_n_data->GetBinError(ibin);
    double bincenter = q_wSid_n_data->GetBinCenter(ibin);
    gr_q_wSid_n_data->SetPoint(ibin,cont, bincenter);
    gr_q_wSid_n_data->SetPointError(ibin,err,0);
  }
  for(int ibin=0;ibin<q_wSid_n_mc->GetNbinsX();ibin++){
    double cont = q_wSid_n_mc->GetBinContent(ibin);
    double err = q_wSid_n_mc->GetBinError(ibin);
    double bincenter = q_wSid_n_mc->GetBinCenter(ibin);
    gr_q_wSid_n_mc->SetPoint(ibin,cont, bincenter);
    gr_q_wSid_n_mc->SetPointError(ibin,err,0);
  }
  gr_q_wSid_n_data->GetYaxis()->SetRangeUser(0,1.5);
  gr_q_wSid_n_data->SetMarkerStyle(20);
  gr_q_wSid_n_data->GetYaxis()->SetLabelSize(0);
  gr_q_wSid_n_data->Draw("AP");
  gr_q_wSid_n_mc->SetMarkerStyle(20);
  gr_q_wSid_n_mc->SetMarkerColor(2);
  gr_q_wSid_n_mc->Draw("P");

  TCanvas *cq_IMnpipi_wSid_n_2 = new TCanvas("cq_IMnpipi_wSid_n_2","cq_IMnpipi_wSid_n_2",1200,1000);
  cq_IMnpipi_wSid_n_2->Divide(2,1);
  cq_IMnpipi_wSid_n_2->cd(1);
  TH2D* q_IMnpipi_wSid_n_data_2 = (TH2D*)q_IMnpipi_wSid_n_data->Clone("q_IMnpipi_wSid_n_data_2");
  TH2D* q_IMnpipi_wSid_n_mc_2 = (TH2D*)q_IMnpipi_wSid_n_mc->Clone("q_IMnpipi_wSid_n_mc_2");
  q_IMnpipi_wSid_n_data_2->SetTitle("real data");
  q_IMnpipi_wSid_n_data_2->SetContour(10);
  q_IMnpipi_wSid_n_data_2->SetMinimum(1);
  q_IMnpipi_wSid_n_data_2->Draw("cont1z");
  gPad->SetRightMargin(0);
  cq_IMnpipi_wSid_n_2->cd(2);
  q_IMnpipi_wSid_n_mc_2->SetTitle("MC");
  q_IMnpipi_wSid_n_mc_2->SetContour(10); 
  q_IMnpipi_wSid_n_mc_2->SetMinimum(1);
  q_IMnpipi_wSid_n_mc_2->SetMaximum(q_IMnpipi_wSid_n_data_2->GetMaximum());
  q_IMnpipi_wSid_n_mc_2->Draw("cont1z");
  gPad->SetLeftMargin(0);
  
  TCanvas *cq_IMnpipi_wSid_n_sub = new TCanvas("cq_IMnpipi_wSid_n_sub","cq_IMnpipi_wSid_n_sub");
  TH2D* q_IMnpipi_wSid_n_sub = (TH2D*)q_IMnpipi_wSid_n_data->Clone("q_IMnpipi_wSid_n_data");
  q_IMnpipi_wSid_n_sub->Add(q_IMnpipi_wSid_n_mc,-1.0);
  q_IMnpipi_wSid_n_sub->Divide(q_IMnpipi_wSid_n_data);
  q_IMnpipi_wSid_n_sub->SetTitle("(data-MC)/data");
  q_IMnpipi_wSid_n_sub->RebinX(2);
  q_IMnpipi_wSid_n_sub->RebinY(2);
  q_IMnpipi_wSid_n_sub->Scale(0.25);
  q_IMnpipi_wSid_n_sub->GetZaxis()->SetRangeUser(-0.3,0.3);
  q_IMnpipi_wSid_n_sub->Draw("colz");



}
