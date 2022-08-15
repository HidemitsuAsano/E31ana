#include <cmath>
#include <iostream>

Double_t FormP(Double_t *x,Double_t *par)
{

  return par[1]*std::pow(x[0]/par[0],2.0)*std::exp(-1.0*std::pow(x[0]/par[0],2.0));

}


Double_t BWandFormP(Double_t *x,Double_t *par)
{
  Double_t r1 = par[0]*std::pow(par[1]/2.0,2.0)/(std::pow((x[0]-par[2]),2.0)+std::pow(par[1]/2.0,2.0));
 
  Double_t r2 = std::pow(x[1]/par[4],2.0)*std::exp(-1.0*std::pow(x[1]/par[4],2.0));

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
  gROOT->ForceStyle();
  //gROOT->SetBatch();
  TFile *file = TFile::Open("cs_lpim_killcombi.root");
  file->cd();
  gStyle->SetOptStat(0);
  TH2F* CS_q_IMppipi_p_wL_nop2 = (TH2F*)file->Get(Form("CS_q_IMppipi_p_wL_nop2_acc%d",0));//divided by q and M bin width
  TH2F* CS_q_IMppipi_p_wL_wp2 = (TH2F*)file->Get(Form("CS_q_IMppipi_p_wL_wp2_acc%d",0));//divided by q and M bin width
  
  TCanvas *cCS_q_IMppipi_p_wL_nop2 = new TCanvas("cCS_q_IMppipi_p_wL_nop2","cCS_q_IMppipi_p_wL_nop2",1000,1000);
  CS_q_IMppipi_p_wL_nop2->Draw("colz");

  TCanvas *cCS_q_IMppipi_p_wL_wp2 = new TCanvas("cCS_q_IMppipi_p_wL_wp2","cCS_q_IMppipi_p_wL_wp2",1000,1000);
  CS_q_IMppipi_p_wL_wp2->Draw("colz");
 
  TCanvas *cCS_sum = new TCanvas("cCS_sum","cCS_sum",1000,800);
  //CS_q_IMppipi_p_wL_nop2->SetMaximum(0.0045);
  //CS_q_IMppipi_p_wL_nop2->Draw("colz");
  //CS_q_IMppipi_p_wL_wp2->SetMaximum(CS_q_IMppipi_p_wL_nop2->GetMaximum());
  //CS_q_IMppipi_p_wL_wp2->Draw("colzsame");
  TH2F* CS_sum = (TH2F*)CS_q_IMppipi_p_wL_nop2->Clone("CS_sum");
  CS_sum->SetTitle("CS_sum");
  CS_sum->Add(CS_q_IMppipi_p_wL_wp2);
  CS_sum->Draw("colz");//divided by q and M bin width

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
      if(acccont<0.005){//remove small acceptance bin
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

  TF2 *f2 = new TF2("f2",VGandLandau,1.32,1.44,0.35,0.65,6);
  f2->SetParLimits(0,0,0.00095);//0: normalization
  f2->SetParameter(0,0.00025);//
  f2->FixParameter(1,0.040);//Lorenz half width (FWHM/2 ) 39.4 is full width of S(1385)- 
  f2->FixParameter(2,1.3872);//PDG mass
  f2->SetParameter(3,0.51);//landau mpv
  f2->SetParLimits(3,0.47,0.55);
  f2->SetParameter(4,0.09);//Landau sigma
  f2->SetParLimits(4,0.02,0.11);//landau sigma
  //f2->FixParameter(5,0.0006);//Voigt sigma = resolution
  f2->SetParLimits(5,0.005,0.01);//Voigt sigma = resolution
  f2->SetNpx(8);//0.12/8 = 0.015 
  f2->SetNpy(6);//0.3/6 = 0.05
  f2->Print("base");
  CS_sum_fit->Fit("f2","R","");
  std::cout << "Chi2/ndf:  " <<   f2->GetChisquare() << " / " <<  f2->GetNDF() << std::endl;
  //f2->Draw("cont1 ");

  TH2D* f2hist = (TH2D*)f2->GetHistogram();
  f2hist->SetLineColor(2);
  f2hist->Draw("colz");
  CS_sum_fit->Print("base");
  f2hist->Print("base");

  cfittest->cd(1);
  const int binx1320 = CS_sum_fit->GetXaxis()->FindBin(1.32);
  const int binx1440 = CS_sum_fit->GetXaxis()->FindBin(1.44);
  const int biny350 = CS_sum_fit->GetYaxis()->FindBin(0.35);
  const int biny650 = CS_sum_fit->GetYaxis()->FindBin(0.65);
  TH1D* CS_sum_fit_px = (TH1D*)CS_sum_fit->ProjectionX("CS_sum_fit_cut_px",biny350,biny650-1);
  CS_sum_fit_px->SetTitle("projection q: 350-650");
  CS_sum_fit_px->GetXaxis()->SetRangeUser(1.2,1.6);
  CS_sum_fit_px->Draw("HE");
  f2hist->SetFillColor(0);
  TH1D* f2hist_px = (TH1D*)f2hist->ProjectionX("f2hist_px");
  std::cout << "width x " << f2hist_px->GetBinWidth(1) << std::endl;
  //f2hist_px->Rebin(5);
  f2hist_px->Draw("same");


  //q distribution
  cfittest->cd(4);
  TH1D* CS_sum_fit_py = (TH1D*)CS_sum_fit->ProjectionY("CS_sum_fit_cut_py",binx1320,binx1440-1);
  CS_sum_fit_py->SetTitle("M: 1320-1440");
  CS_sum_fit_py->GetXaxis()->SetRangeUser(0,0.75);
  CS_sum_fit_py->Draw("HE");
  TH1D* f2hist_py = (TH1D*)f2hist->ProjectionY("f2hist_py");
  f2hist_py->Draw("same");
  std::cout << "width y " << f2hist_py->GetBinWidth(1) << std::endl;
  
   
  //all range of region of intereset
  TF2 *f3 = new TF2("f3",VGandLandau,1.32,1.44,0.20,0.65,6);
  f3->SetParameters(f2->GetParameters());
  f3->SetNpx(8);//12/8 = 15 MeV bin
  f3->SetNpy(9);//0.45/9 = 50 MeV/c bin
  f3->Print();
  

  TH2D* f3hist = (TH2D*)f3->GetHistogram();
  cfittest->cd(1);
  for(int ix=0;ix<=f3hist->GetNbinsX();ix++){
    for(int iy=0;iy<=f3hist->GetNbinsY();iy++){
      f3hist->SetBinError(ix,iy,0);
    }
  }
  TH1D* f3hist_px = (TH1D*)f3hist->ProjectionX("f3hist_px");
  f3hist_px->SetLineColor(4);
  f3hist_px->SetFillColor(0);
  //f3hist_px->Draw("same");

  cfittest->cd(4);
  const int bin1320_f3 = f3hist->GetXaxis()->FindBin(1.32);
  const int bin1440_f3 = f3hist->GetXaxis()->FindBin(1.44);
  TH1D* f3hist_py = (TH1D*)f3hist->ProjectionY("f3hist_py",bin1320_f3,bin1440_f3);
  f3hist_py->SetLineColor(4);
  f3hist_py->SetFillColor(0);
  //f3hist_py->Draw("same");

  
  TCanvas *callrangecheck = new TCanvas("callrangecheck","callrangecheck");
  callrangecheck->Divide(2,2);
  callrangecheck->cd(3);
  CS_sum_fit->Draw("colz");
  
  callrangecheck->cd(1);
  CS_sum_fit->ProjectionX()->Draw("HE");


  callrangecheck->cd(4);
  CS_sum_fit->ProjectionY()->Draw("HE");




  TCanvas *cCS_q_fit = new TCanvas("cCS_q_fit","cCS_q_fit",1000,800);
  const int bin1320 = CS_sum_fit->GetXaxis()->FindBin(1.320);
  const int bin1440 = CS_sum_fit->GetXaxis()->FindBin(1.44);
  TH1D* CS_q_fit = (TH1D*)CS_sum_fit->ProjectionY("CS_q_fit",bin1320,bin1440-1);
  CS_q_fit->SetMarkerStyle(20);
  CS_q_fit->Draw("E");
  const int bin1320fit = f3hist->GetXaxis()->FindBin(1.32);
  const int bin1440fit = f3hist->GetXaxis()->FindBin(1.44);
  TH1D* f3hist_py2 = (TH1D*)f3hist->ProjectionY("f3hist_py2",bin1320fit,bin1440fit-1);
  
  f3hist_py2->SetLineColor(4);
  f3hist_py2->SetFillColor(0);
  f3hist_py2->Draw("same");
  //TF1 *f1 = new TF1("f1",FormP,0.4,0.8,2);
  //CS_q_fit->Fit("f1","","",0.4,0.75);
  std::cout << __LINE__ << std::endl;
  TCanvas *cCS_q_all = new TCanvas("cCS_q_all","cCS_q_all",1000,800);
  TH1D* CS_q_all = (TH1D*)CS_sum->ProjectionY("CS_q_all",bin1320,bin1440-1);
  CS_q_all->SetMarkerStyle(20);
  CS_q_all->GetXaxis()->SetRangeUser(0,0.65);
  double binwidth_q = CS_q_all->GetXaxis()->GetBinWidth(1)*1000.0; 
  CS_q_all->Scale(binwidth_q);
  f3hist_py2->Scale(binwidth_q);
  CS_q_all->Draw("E");
  f3hist_py2->Draw("same");
  const int binq350 = f3hist_py2->FindBin(0.35);
  const int binq650 = f3hist_py2->FindBin(0.65);
  std::cout << f3hist_py2->Integral(binq350,binq650)  << std::endl;
  std::cout << f3hist_py2->Integral(1,binq350) << std::endl;
  
  //q-dep systematic error
  TCanvas *cCS_q_allsys = new TCanvas("cCS_q_allsys","cCS_q_allsys",1000,800);
  TGraphAsymmErrors *gCS_qdep = new TGraphAsymmErrors();
  //new TGraphAsymmErrors(CS_q_all);
   
  f3hist_py2->Print("all");
  CS_q_all->Print("all");
  gCS_qdep->Print();
  int n  = CS_q_all->GetNbinsX();
  for(int ip=1;ip<9;ip++){
    double valfit = f3hist_py2->GetBinContent(ip+1);
    double valmeasured = CS_q_all->GetBinContent(ip+5);
    double bincen = CS_q_all->GetBinCenter(ip+5);
    double err = valmeasured - valfit;
    if(err>0){
      gCS_qdep->SetPoint(ip,bincen,valmeasured);
      gCS_qdep->SetPointEXlow(ip,0.025);
      gCS_qdep->SetPointEXhigh(ip,0.025);
      gCS_qdep->SetPointEYhigh(ip,fabs(valmeasured*0.1));
      gCS_qdep->SetPointEYlow(ip,fabs(err));
    }else{
      gCS_qdep->SetPoint(ip,bincen,valmeasured);
      gCS_qdep->SetPointEXlow(ip,0.025);
      gCS_qdep->SetPointEXhigh(ip,0.025);
      gCS_qdep->SetPointEYhigh(ip,fabs(err));
      gCS_qdep->SetPointEYlow(ip,valmeasured*0.1);
    }
  }
  
  CS_q_all->Draw("E");
  gCS_qdep->GetXaxis()->SetRangeUser(0,0.65);
  gCS_qdep->SetFillStyle(3002);
  gCS_qdep->SetFillColor(3);
  gCS_qdep->SetMarkerColor(3);
  gCS_qdep->SetLineColor(3);
  gCS_qdep->Draw("5");
  gCS_qdep->Print();
  //gCS_qdep->Draw("ap");
  

  //forget about fitting 2d fitting, just 
 
  TH2F* CS_sum_nofit = (TH2F*)CS_sum->Clone("CS_sum_nofit");
  CS_sum_nofit->SetName("CS_sum_nofit");
  const int binq200nofit = CS_sum_nofit->GetYaxis()->FindBin(0.20);
  const int binq350nofit = CS_sum_nofit->GetYaxis()->FindBin(0.35);
  const int binq650nofit = CS_sum_nofit->GetYaxis()->FindBin(0.65);
  const int binM1440  = CS_sum_nofit->GetXaxis()->FindBin(1.440);
  const double br_s1385ToLambdapi = 0.87;
  const double br_s1385TopiSigma = 0.117;
  const double br_s1385TopiSigma_err = 0.015;
  const double IsospinCGFactor = 2.0;  

  double binwidthM = CS_sum_nofit->GetXaxis()->GetBinWidth(1)*1000.0; 
  double binwidthq = CS_sum_nofit->GetYaxis()->GetBinWidth(1)*1000.0; 
  //CS_sum_nofit->Scale(br_s1385TopiSigma/2.0/br_s1385ToLambdapi/IsospinCGFactor);

  TH1D* CS_sum_nofit_qlow = (TH1D*)CS_sum_nofit->ProjectionX("CS_sum_nofit_qlow",binq200nofit,binq350nofit-1);
  TH1D* CS_sum_nofit_qhi = (TH1D*)CS_sum_nofit->ProjectionX("CS_sum_nofit_qhi",binq350nofit,binq650nofit);
  TH1D* CS_sum_nofit_qall = (TH1D*)CS_sum_nofit->ProjectionX("CS_sum_nofit_qall",binq200nofit,binq650nofit);
   
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  TCanvas *cnofitqall = new TCanvas("cnofitall","cnofitall",1000,800);
  CS_sum_nofit_qall->GetXaxis()->SetRangeUser(1.3,1.6);
  CS_sum_nofit_qall->Scale(binwidthq);//multiply by q bin width
  CS_sum_nofit_qall->Draw();
  
  CS_sum_nofit_qlow->SetLineColor(2);
  CS_sum_nofit_qlow->Scale(binwidthq);//multiply by q bin width
  CS_sum_nofit_qlow->Draw("HEsame");
  CS_sum_nofit_qhi->SetLineColor(3);
  CS_sum_nofit_qhi->Scale(binwidthq);//multiply by q bin width
  CS_sum_nofit_qhi->Draw("HEsame");

  TCanvas *cnofit = new TCanvas("cnofit","cnofit",1000,800);
  cnofit->Divide(2,1);
  cnofit->cd(1);
  CS_sum_nofit_qlow->GetXaxis()->SetRangeUser(1.3,1.6);
  CS_sum_nofit_qlow->Draw("HE");

  cnofit->cd(2);
  CS_sum_nofit_qhi->GetXaxis()->SetRangeUser(1.3,1.6);
  CS_sum_nofit_qhi->Draw("HE");
  

  TH1D* CS_M_measured[8];
  TH1D* CS_M_fit[8];
  
  TCanvas *cM = new TCanvas("cM","cM",1200,800);
  cM->Divide(4,2);
   
  //Note CS_sum and f3hist is divided by binwidth of M and q 
  //must be multiply by bin width after projection
  double binwidth_M = f3hist->GetXaxis()->GetBinWidth(1)*1000.0; 
  for(int iqbin=0;iqbin<8;iqbin++){
    CS_M_measured[iqbin] = (TH1D*)CS_sum->ProjectionX(Form("CS_M_measured_bin%d",iqbin+6),iqbin+6,iqbin+6); 
    std::cout << CS_sum->GetYaxis()->GetBinLowEdge(iqbin+6) << "  " << f3hist->GetYaxis()->GetBinLowEdge(iqbin+2) << std::endl;
    CS_M_fit[iqbin] = (TH1D*)f3hist->ProjectionX(Form("CS_M_fit_bin%d",iqbin+1),iqbin+2,iqbin+2); 
    cM->cd(iqbin+1);
    CS_M_measured[iqbin]->GetXaxis()->SetRangeUser(1.2,1.6);
    CS_M_measured[iqbin]->Scale(binwidthq);
    CS_M_fit[iqbin]->Scale(binwidthq);
    CS_M_measured[iqbin]->SetMaximum(CS_M_measured[iqbin]->GetMaximum()*1.5);
    CS_M_measured[iqbin]->SetTitle(Form("q %0.2f-%0.2f",0.25+0.05*iqbin,0.30+0.05*iqbin));
    CS_M_measured[iqbin]->Draw("HE");
    CS_M_fit[iqbin]->SetFillColor(0);
    CS_M_fit[iqbin]->Draw("Esame");
    std::cout << CS_M_fit[iqbin]->Integral() << std::endl;
  }  
  

  TGraphAsymmErrors *gr_M_qlow = new TGraphAsymmErrors(CS_sum_nofit_qlow);//divided by M
  TGraphAsymmErrors *gr_M_qhi = new TGraphAsymmErrors(CS_sum_nofit_qhi);//divided by M
  gr_M_qlow->SetName("gr_M_qlow");  
  gr_M_qhi->SetName("gr_M_qhi");  


  TCanvas *cMerr = new TCanvas("cMerr","cMerr",1200,800); 
  cMerr->Divide(2,1); 
  cMerr->cd(1);
  gr_M_qlow->GetXaxis()->SetRangeUser(1.2,1.6);
  gr_M_qlow->SetLineColor(1);  
  gr_M_qlow->SetMarkerColor(1);  
  gr_M_qlow->SetMarkerStyle(20);  
  for(int ip=0;ip<8;ip++){
    gr_M_qlow->RemovePoint(0); 
  }
  gr_M_qlow->Draw("AP");  

  
  
  TH1D* CS_M_qlowErr = (TH1D*) CS_M_fit[0]->Clone("CS_M_qlowerr");
  CS_M_qlowErr->Add(CS_M_fit[1]);
  //CS_M_qlowErr->Draw("same");
  
  TGraphAsymmErrors *gr_M_qlowErr = new TGraphAsymmErrors(CS_M_qlowErr);//divided by M
  gr_M_qlowErr->SetName("gr_M_qlowErr");

  for(int ip=0;ip<gr_M_qlowErr->GetN();ip++){
    double yfit = CS_M_qlowErr->GetBinContent(ip+1);
    double xval = CS_M_qlowErr->GetBinCenter(ip+1);
    int fitbin = CS_M_qlowErr->FindBin(xval);
    double ymes = CS_sum_nofit_qlow->GetBinContent(fitbin);
     
    if(ymes<yfit){
      gr_M_qlowErr->SetPointEYhigh(ip,yfit-ymes);
    }else{
      gr_M_qlowErr->SetPointEYlow(ip,ymes-yfit);
    }
  }
  gr_M_qlowErr->Draw("5P");
  gr_M_qlow->Draw("P");  
  
  cMerr->cd(2);
  gr_M_qhi->GetXaxis()->SetRangeUser(1.2,1.6);
  gr_M_qhi->SetLineColor(1);  
  gr_M_qhi->SetMarkerColor(1);  
  gr_M_qhi->SetMarkerStyle(20);  
  gr_M_qhi->Draw("AP");  
  
  TH1D* CS_M_qhiErr = (TH1D*) CS_M_fit[2]->Clone("CS_M_qhierr");
  for(int iq=2;iq<8;iq++){
    CS_M_qhiErr->Add(CS_M_fit[iq]);
  }
  CS_M_qhiErr->Draw("same");
  
  TGraphAsymmErrors *gr_M_qhiErr = new TGraphAsymmErrors(CS_sum_nofit_qhi);//divided by M
  gr_M_qhiErr->SetName("gr_M_qhiErr");


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

  
  TFile *fout = new TFile("CSLpimFit.root","RECREATE");
  fout->cd();
  f3hist->Write();
  CS_sum->Write();


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
