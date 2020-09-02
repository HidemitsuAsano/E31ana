{
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(111111);

  TFile *f = new TFile("kpp_50_50_Lp_reconst.root");
  //TFile *f = new TFile("lambda_reconst.root");

  TCanvas *can[3];
  for( int i=0; i<3; i++ ){
    TString str(Form("can%d", i));
    can[i] = new TCanvas(str, "", 1800, 600);
    can[i]->Divide(6, 3);
  }

  cerr<<endl;
  cerr<<"*** difference of momentum: meas-gene, kfit-gene, meas-kfit"<<endl;
  cerr<<"*** 0:beam_K(K+), 1:L, 2:n, 3:p, 4:p from L, 5:pi- from L "<<endl;
  cerr<<endl;
  TH1F *his;
  TString def[6] = {"beam_K-", "#Lambda", "n", "p", "p from #Lambda", "#pi- from #Lambda"}
  char* name[3][3] = {
    {"p_gene_%d", "p_meas_%d", "p_kfit_%d"},
    {"p_meas_gene_%d", "p_kfit_gene_%d", "p_meas_kfit_%d"},
    {"p_per_meas_gene_%d", "p_per_kfit_gene_%d", "p_per_meas_kfit_%d"}
  };
  char* tit[3][3] = {
    {" :: generated p [GeV/c]", " :: measured p [GeV/c]", " :: kfitted p [GeV/c]"},
    {" :: measured-generated #deltap [GeV/c]", " :: kfitted-generated #deltap [GeV/c]", " :: measured-kfitted #deltap [GeV/c]"},
    {" :: measured-generated #deltap/p [%]", " :: kfitted-generated #deltap/p [%]", " :: measured-kfitted #deltap/p [%]"}
  };

  double mea[3][18];
  double rms[3][18];
  TString title;
  for( int i=0; i<3; i++ ){
    for( int j=0; j<18; j++ ){
      can[i]->cd(j+1);
      his = (TH1F*)f->Get(Form(name[i][j/6], j%6));
      his->Draw();
      title = def[j%6]+tit[i][j/6];
      his->SetXTitle(title);
      if( i<2 ){
	cerr<<j%6<<" | "<<his->GetMean()*1000<<", "<<his->GetRMS()*1000<<" MeV/c"<<endl;
	mea[i][j] = his->GetMean()*1000;
	rms[i][j] = his->GetRMS()*1000;
      }
      else{
	cerr<<j%6<<" | "<<his->GetMean()<<", "<<his->GetRMS()<<" %"<<endl;
	mea[i][j] = his->GetMean();
	rms[i][j] = his->GetRMS();
      }
    } // for( int j=0; j<6; j++ ){
    cerr<<endl;
  } // for( int i=0; i<3; i++ ){


  TCanvas *c1;
  c1 = new TCanvas("c1", "", 600, 600);
  vector <double>  xdata, ydata1;
  vector <TString> xlabel;
  for( int i=0; i<6; i++ ){
    xdata.push_back(i+1);
    ydata1.push_back(rms[2][i]);
    xlabel.push_back(def[i]);
  }
  Double_t* xpointer=&(xdata.begin()[0]);
  Double_t* ypointer1=&(ydata1.begin()[0]);
  TGraph *gr1 = new TGraph(xdata.size(), xpointer, ypointer1);
  gr1->SetMarkerStyle(21);
  gr1->SetMarkerColor(kRed);
  gr1->SetMarkerSize(1);
  gr1->GetXaxis()->SetNdivisions(6);
  gr1->GetXaxis()->SetLabelOffset(99);
  gr1->GetYaxis()->SetTitle("[%]");
  gr1->GetYaxis()->SetTitleOffset(1.3);
  gr1->GetYaxis()->SetRangeUser(-2.0, 4.0);
  gr1->SetTitle("");
  gr1->Draw("AP");
  c1->SetGrid();

  vector <double> ydata2;
  for( int i=0; i<6; i++ ){
    ydata2.push_back(rms[2][i+6]);
  }
  Double_t* ypointer2=&(ydata2.begin()[0]);
  TGraph *gr2 = new TGraph(xdata.size(), xpointer, ypointer2);
  gr2->SetMarkerStyle(21);
  gr2->SetMarkerColor(kBlue);
  gr2->SetMarkerSize(1);
  gr2->Draw("same,P");

  vector <double> ydata3;
  for( int i=0; i<6; i++ ){
    ydata3.push_back(mea[2][i]);
  }
  Double_t* ypointer2=&(ydata3.begin()[0]);
  TGraph *gr3 = new TGraph(xdata.size(), xpointer, ypointer2);
  gr3->SetMarkerStyle(24);
  gr3->SetMarkerColor(kRed);
  gr3->SetMarkerSize(1);
  gr3->Draw("same,P");

  vector <double> ydata4;
  for( int i=0; i<6; i++ ){
    ydata4.push_back(mea[2][i+6]);
  }
  Double_t* ypointer2=&(ydata4.begin()[0]);
  TGraph *gr4 = new TGraph(xdata.size(), xpointer, ypointer2);
  gr4->SetMarkerStyle(24);
  gr4->SetMarkerColor(kBlue);
  gr4->SetMarkerSize(1);
  gr4->Draw("same,P");


  TLatex t;
  t.SetTextAngle(-20);
  t.SetTextSize(0.04);
  t.SetTextAlign(-10);
  Float_t tx,ty;
  //ty = gPad->GetUymin() - 0.06*gr1->GetYaxis()->GetBinWidth(1);
  ty = -2.3;
  for(int i=0;i<xlabel.size();i++){
    tx=xdata.at(i);
    t.DrawLatex(tx,ty,(xlabel.at(i)).Data());
  }

  TLegend* leg;
  leg = new TLegend(0.15,0.8,0.55,0.98);
  leg->AddEntry(gr1,"#deltap/p(measured-generated)","P");
  leg->AddEntry(gr2,"#deltap/p(kfitted-generated)","P");
  leg->AddEntry(gr3,"#deltap(measured-generated)","P");
  leg->AddEntry(gr4,"#deltap(kfitted-generated)","P");
  leg->Draw();



  can[0]->Print("tmp.pdf(");
  can[1]->Print("tmp.pdf");
  can[2]->Print("tmp.pdf");
  c1->Print("tmp.pdf)");
}
