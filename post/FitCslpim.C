
Double_t BWandFormS(Double_t *x,Double_t *par)
{


  return 0;
}

Double_t FormP(Double_t *x,Double_t *par)
{

  return par[1]*pow(x[0]/par[0],2.0)*exp(-1.0*pow(x[0]/par[0],2.0));

}


Double_t BWandFormP(Double_t *x,Double_t *par)
{
  Double_t r1 = par[0]*pow(par[1]/2.0,2.0)/((pow(x[0]-par[2]),2.0)+pow(par[1]/2.0,2.0));
 
  Double_t r2 = pow(x[1]/par[4],2.0)*exp(-1.0*pow(x[1]/par[4],2.0));

  return r1*r2;

}

//Fit with BW and some phenomenological model
Double_t BWandLandau(Double_t *x,Double_t *par)
{
  //BW and scaling factor
  //Double_t BW = par[0]*pow(par[1]/2.0,2.0)/((pow(x[0]-par[2]),2.0)+pow(par[1]/2.0,2.0));
  Double_t BW = TMath::BreitWigner(x[0],par[2],par[1]);
  //Phenomenological model
  Double_t Ph= TMath::Landau(x[1],par[3],par[4]);


  return par[0]*BW*Ph;

}

//Voigt and Landau
Double_t VGandLandau(Double_t *x,Double_t *par)
{
  Double_t VG = TMath::Voigt(x[0]-par[2],par[5],par[1],4);
  //Phenomenological model
  Double_t Ph= TMath::Landau(x[1],par[3],par[4]);
  
  return par[0]*VG*Ph;

}



void FitCslpim()
{
  TFile *file = TFile::Open("cs_lpim_killcombi.root");
  file->cd();
  gStyle->SetOptStat(0);
  TH2F* CS_q_IMppipi_p_wL_nop2 = (TH2F*)file->Get(Form("CS_q_IMppipi_p_wL_nop2_acc%d",0));
  TH2F* CS_q_IMppipi_p_wL_wp2 = (TH2F*)file->Get(Form("CS_q_IMppipi_p_wL_wp2_acc%d",0));
  
  TCanvas *cCS_q_IMppipi_p_wL_nop2 = new TCanvas("cCS_q_IMppipi_p_wL_nop2","cCS_q_IMppipi_p_wL_nop2",1000,800);
  CS_q_IMppipi_p_wL_nop2->Draw("colz");

  TCanvas *cCS_q_IMppipi_p_wL_wp2 = new TCanvas("cCS_q_IMppipi_p_wL_wp2","cCS_q_IMppipi_p_wL_wp2",1000,800);
  CS_q_IMppipi_p_wL_wp2->Draw("colz");
 
  TCanvas *cCS_sum = new TCanvas("cCS_sum","cCS_sum",1000,800);
  //CS_q_IMppipi_p_wL_nop2->SetMaximum(0.0045);
  //CS_q_IMppipi_p_wL_nop2->Draw("colz");
  //CS_q_IMppipi_p_wL_wp2->SetMaximum(CS_q_IMppipi_p_wL_nop2->GetMaximum());
  //CS_q_IMppipi_p_wL_wp2->Draw("colzsame");
  TH2F* CS_sum = (TH2F*)CS_q_IMppipi_p_wL_nop2->Clone("CS_sum");
  CS_sum->SetTitle("CS_sum");
  CS_sum->Add(CS_q_IMppipi_p_wL_wp2);
  CS_sum->Draw("colz");

  TCanvas *cCS_acc_sum = new TCanvas("cCS_acc_sum","cCS_acc_sum",1000,800);
  TH2F* acc_nop2 = (TH2F*)file->Get("q_IMppipi_p_wL_nop2_acc_clean0");
  TH2F* acc_wp2 = (TH2F*)file->Get("q_IMppipi_p_wL_wp2_acc_clean0");
  TH2F* acc_sum = (TH2F*)acc_nop2->Clone("acc_sum");
  acc_sum->Add(acc_wp2);
  acc_sum->SetTitle("acc_sum");
  acc_sum->GetZaxis()->SetRangeUser(0,0.1);
  acc_sum->Draw("colz");
  
  TCanvas *cCS_sum_fit = new TCanvas("cCS_sum_fit","cCS_sum_fit",1000,800);
  TH2F* CS_sum_fit = (TH2F*)CS_sum->Clone("CS_sum_fit");
  CS_sum_fit->SetName("CS_sum_fit");
  CS_sum_fit->SetTitle("CS_sum_fit");
  for(int ix=0;ix<CS_sum_fit->GetNbinsX();ix++){
    for(int iy=0;iy<CS_sum_fit->GetNbinsY();iy++){
      double acccont = acc_sum->GetBinContent(ix,iy);
      if(acccont<0.01){
        CS_sum_fit->SetBinContent(ix,iy,0);
        CS_sum_fit->SetBinError(ix,iy,0);
      }
    }
  }
  CS_sum_fit->Draw("colz");
  
  //TF2 *f2 = new TF2("f2",BWandFormP,1.30,1.40,0.4,0.8,5);
  //f2->SetParameter(1,0.1);
  //f2->FixParameter(2,1.38);
  //f2->SetRange(1.30,1.40,0.4,0.8);
  //CS_sum_fit->Fit("f2","R","");
  //f2->Draw("cont1 same");
  
  TCanvas *cfittest = new TCanvas("cfittest","cfittest");
  cfittest->Divide(2,2);
  cfittest->cd(3);
  
  /*
  TF2 *f2 = new TF2("f2",BWandLandau,1.32,1.45,0.4,0.75,5);
  f2->SetParLimits(0,0,10000);
  f2->SetParameter(0,0.0001);
  f2->SetParameter(1,0.045);
  f2->FixParameter(2,1.3872);//PDG mass
  f2->SetParameter(3,0.47);
  f2->SetParLimits(3,0.2,0.5);
  f2->SetParameter(4,0.11);
  f2->SetParLimits(4,0.1,1000000);
  f2->SetNpx(100);
  f2->SetNpy(100);
  //CS_sum_fit->Fit("f2","R","");
  //f2->Draw("cont1 same");
  f2->Print("base");
  //f2->Draw("cont1 ");
  TH2D* f2hist = (TH2D*)f2->GetHistogram();
  f2hist->SetLineColor(2);
  f2hist->Draw("colz");
  */

  TF2 *f2 = new TF2("f2",VGandLandau,1.32,1.44,0.39,0.75,6);
  f2->SetParLimits(0,0,10000);
  f2->SetParameter(0,0.00025);
  f2->SetParameter(1,0.045);
  f2->FixParameter(2,1.3872);//PDG mass
  f2->SetParameter(3,0.51);//landau mpv
  f2->SetParLimits(3,0.51,0.55);
  f2->SetParameter(4,0.09);
  f2->SetParLimits(4,0.01,0.1);//landau sigma
  f2->FixParameter(5,0.001);
  f2->SetNpx(8);
  f2->SetNpy(12);
  f2->Print("base");
  CS_sum_fit->Fit("f2","R","");
  //f2->Draw("cont1 ");

  TH2D* f2hist = (TH2D*)f2->GetHistogram();
  f2hist->SetLineColor(2);
  f2hist->Draw("colz");
  CS_sum_fit->Print("base");
  f2hist->Print("base");

  cfittest->cd(1);
  const int binx1320 = CS_sum_fit->GetXaxis()->FindBin(1.32);
  const int binx1440 = CS_sum_fit->GetXaxis()->FindBin(1.44);
  const int biny400 = CS_sum_fit->GetYaxis()->FindBin(0.39);
  const int biny750 = CS_sum_fit->GetYaxis()->FindBin(0.75);
  CS_sum_fit->ProjectionX("CS_sum_fit_px",biny400,biny750)->Draw("HE");
  f2hist->SetFillColor(0);
  TH1D* f2hist_px = (TH1D*)f2hist->ProjectionX("f2hist_px");
  std::cout << "width x " << f2hist_px->GetBinWidth(1) << std::endl;
  //f2hist_px->Rebin(5);
  f2hist_px->Draw("same");

  cfittest->cd(4);
  CS_sum_fit->ProjectionY("CS_sum_fit_py",binx1320,binx1440)->Draw("HE");
  TH1D* f2hist_py = (TH1D*)f2hist->ProjectionY();
  f2hist_py->Draw("same");
  std::cout << "width y " << f2hist_py->GetBinWidth(1) << std::endl;
  
  
  
  

  TF2 *f3 = new TF2("f3",VGandLandau,1.32,1.44,0.21,0.75,6);
  f3->SetParameters(f2->GetParameters());
  f3->SetNpx(8);
  f3->SetNpy(18);
  f3->Print();
  

  TH2D* f3hist = (TH2D*)f3->GetHistogram();
  cfittest->cd(1);
  TH1D* f3hist_px = (TH1D*)f3hist->ProjectionX("f3hist_px");
  f3hist_px->SetLineColor(4);
  f3hist_px->SetFillColor(0);
  f3hist_px->Draw("same");

  cfittest->cd(4);
  TH1D* f3hist_py = (TH1D*)f3hist->ProjectionY("f3hist_py");
  f3hist_py->SetLineColor(4);
  f3hist_py->SetFillColor(0);
  f3hist_py->Draw("same");

  TCanvas *cCS_q_fit = new TCanvas("cCS_q_fit","cCS_q_fit",1000,800);
  const int bin1360 = CS_sum_fit->GetXaxis()->FindBin(1.35);
  const int bin1400 = CS_sum_fit->GetXaxis()->FindBin(1.395);
  TH1D* CS_q_fit = (TH1D*)CS_sum_fit->ProjectionY("CS_q_fit",bin1360,bin1400);
  CS_q_fit->SetMarkerStyle(20);
  CS_q_fit->Draw("E");
  const int bin1360fit = f3hist->GetXaxis()->FindBin(1.35);
  const int bin1400fit = f3hist->GetXaxis()->FindBin(1.395);
  TH1D* f3hist_py2 = (TH1D*)f3hist->ProjectionY("f3hist_py2",bin1360fit,bin1400fit);
  f3hist_py2->SetLineColor(4);
  f3hist_py2->SetFillColor(0);
  f3hist_py2->Draw("same");
  //TF1 *f1 = new TF1("f1",FormP,0.4,0.8,2);
  //CS_q_fit->Fit("f1","","",0.4,0.75);



  TCanvas *c = NULL;
  TSeqCollection *SCol = gROOT->GetListOfCanvases();
  int size = SCol->GetSize();
  TIter next(SCol);
  TString pdfname = "fitlpim.pdf";
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
    c->GetListOfPrimitives();//->Print();
    //std::cout << c->GetName() << std::endl;
    //make 1 pdf file
    if(i==0) c->Print(pdfname+"(",Form("pdf Title:%s",c->GetTitle()));
    else if(i==size-1)c->Print(pdfname+")",Form("pdf Title:%s",c->GetTitle())); 
    else c->Print(pdfname,Form("pdf Title:%s",c->GetTitle())); 
    //make separated pdf files
    //c->Print(Form("pdf/%s.pdf",c->GetTitle()));
  }






  return;


  /*
  TCanvas *cCS_q = new TCanvas("cCS_q","cCS_q",1000,800);
  TH1D* CS_q_nop20 = (TH1D*)file->Get("CS_q_nop20");
  CS_q_nop20->Draw("");
  
  TH1D* CS_q_wp20 = (TH1D*)file->Get("CS_q_wp20");
  CS_q_wp20->Draw("same");
  const double qbinwidth = CS_q_nop20->GetXaxis()->GetBinWidth(2);
  CS_q_nop20->Scale(1./qbinwidth);
  CS_q_wp20->Scale(1./qbinwidth);
  TCanvas *cCS_q_gr = new TCanvas("cCS_q_gr","cCS_q_gr",1000,800);
  TGraphErrors *gr_nop2 = new TGraphErrors();
  int ip=0;
  for(int i=0;i<CS_q_nop20->GetNbinsX();i++){
    double cont = CS_q_nop20->GetBinContent(i);
    double err = CS_q_nop20->GetBinError(i);
    double bincent = CS_q_nop20->GetBinCenter(i);
    if(cont>0.0 && 0.4<bincent && bincent<0.75){
      gr_nop2->SetPoint(ip,bincent,cont);
      gr_nop2->SetPointError(ip,0.015,err);
      ip++;
    }
  }
  //inoue
  gr_nop2->SetPoint(ip,0.315,1.00582900000000008e+00);
  gr_nop2->SetPointError(ip,0.015,0.203361);
  gr_nop2->SetMarkerStyle(20);
  //gr_nop2->SetMarkerColor(20);
  gr_nop2->GetXaxis()->SetLimits(0.0,1.5);
  gr_nop2->SetMinimum(0.0);
  gr_nop2->SetMaximum(CS_q_nop20->GetMaximum());
  gr_nop2->GetXaxis()->SetTitle("Mom. Transfer [GeV/c]");
  gr_nop2->GetXaxis()->CenterTitle();
  gr_nop2->Draw("AP");

  TGraphErrors *gr_wp2 = new TGraphErrors();
  ip=0;
  for(int i=0;i<CS_q_wp20->GetNbinsX();i++){
    double cont = CS_q_wp20->GetBinContent(i);
    double err = CS_q_wp20->GetBinError(i);
    if(cont>0.0){
      double bincent = CS_q_wp20->GetBinCenter(i);
      gr_wp2->SetPoint(ip,bincent,cont);
      gr_wp2->SetPointError(ip,0.015,err);
      ip++;
    }
  }
  gr_wp2->SetMarkerStyle(20);
  gr_wp2->SetMarkerColor(2);
  gr_wp2->SetLineColor(2);
  gr_wp2->Draw("P");
  
  */
  


}
