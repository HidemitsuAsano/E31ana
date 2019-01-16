//
// for the covariance matrix evaluation
//
void eval_cov(const char*filename="")
{
  gROOT->SetStyle("Plain");
  //gStyle->SetOptStat(111111);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  gStyle->SetTitleYOffset(1.8);

  TFile *f = new TFile(filename);
  if(!f->IsOpen()){
    cout << "file not found !" << filename << endl;
    return;
  }
  
  TString outfilename = string(filename);
  outfilename.Insert(outfilename.Sizeof()-6,"_eval");
  outfilename.Replace(outfilename.Sizeof()-5,5,"pdf");
  std::cout << outfilename << std::endl;
  //TFile *f = new TFile("simdata012_ana00/sim_0000.root");

  const int num = 6;

  TCanvas *can[num];
  TCanvas *can_vs[num];
  for( int i=0; i<num; i++ ){
    TString str(Form("can%d", i));
    can[i] = new TCanvas(str, "", 800, 800);
    can[i]->Divide(2, 2);
    can_vs[i] = new TCanvas(Form("can_vs%d",i), "",1000,800);
    can_vs[i]->Divide(2, 2);
  }

  // beam_K(K+), n, p, n from L, pi- from L
  double val[num][16];
  double valfit[num][16];
  double pole[num][4];
  double sigm[num][4];
  double pole_e[num][4];
  double sigm_e[num][4];
  for( int i=0; i<num; i++ ){
    int n = 0;
    for( int j=0; j<4; j++ ){
      for( int k=0; k<4; k++ ){
	if( j==k ){ // only diagonal components
	  can[i]->cd(j+1);
	  TH1F *his = (TH1F*)f->Get(Form("cov_%d_%d_%d", i, j, k));
	  TH2F *his2 = (TH2F*)f->Get(Form("cov_mom_%d_%d_%d", i, j, k));
    if(i==0) his->GetXaxis()->SetRangeUser(-0.02,0.02);
    if(i==1){
      his->Rebin(2);
      his->GetXaxis()->SetRangeUser(-0.02,0.02);
    }
    if(i==2){
      his->Rebin(2); 
      his->GetXaxis()->SetRangeUser(-0.02,0.02);
    }
    if(i==3){
      his->Rebin(2);
      his->GetXaxis()->SetRangeUser(-0.02,0.02);
    }
    if(i==4){
      if(j<2)his->Rebin(2);
      his->GetXaxis()->SetRangeUser(-0.02,0.02);
    }
    if(i==5) {
      his->GetXaxis()->SetRangeUser(-0.04,0.04);
	  }
    his->Fit("gaus","q","",-0.01,0.01);
	  his->GetFunction("gaus")->SetLineColor(4);
	  double pol = his->GetFunction("gaus")->GetParameter(1);
	  double sigfit = his->GetFunction("gaus")->GetParameter(2);
	  //his->Fit("gaus","","",pol-3*sig,pol+3*sig);
	  gaus->SetLineColor(3);
    his->GetXaxis()->CenterTitle();
    if(i==0) his->SetTitle("K (beam)");
    if(i==1) his->SetTitle("#pi (prompt)");
    if(i==2) his->SetTitle("#Sigma ");
    if(i==3) his->SetTitle("missing n");
    if(i==4) his->SetTitle("n from #Sigma");
    if(i==5) his->SetTitle("#pi from #Sigma");
    if(k==0){ his->SetXTitle("px (reco.-gen.) [GeV/c]");}
    if(k==1){ his->SetXTitle("py (reco.-gen.) [GeV/c]");}
    if(k==2){ his->SetXTitle("pz (reco.-gen.) [GeV/c]");}
    if(k==3){ his->SetXTitle("E  (reco.-gen.) [GeV]");}
    his->GetXaxis()->CenterTitle();
    his->GetYaxis()->CenterTitle();
    
	  if(i==1) his->Fit("gaus","q","",pol-1.5*sigfit,pol+1.5*sigfit); 
	  else if(k==3 && i==4 ) his->Fit("gaus","q","",pol-1.0*sigfit,pol+1.0*sigfit);
    else if( i==5 ) his->Fit("gaus","q","",pol-2.0*sigfit,pol+2.0*sigfit);
	  else     his->Fit("gaus","q","",pol-3*sigfit,pol+3*sigfit);
    his->Draw("HE");
    //gaus->Draw("same");
    pol = his->GetFunction("gaus")->GetParameter(1);
	  sigfit = his->GetFunction("gaus")->GetParameter(2);
    double sig = his->GetRMS();
	  val[i][n] = sig*sig;
    valfit[i][n] = sigfit*sigfit;
	  pole[i][j] = pol*1000;
	  sigm[i][j] = sigfit*1000;
	  pole_e[i][j] = his->GetFunction("gaus")->GetParError(1)*1000;
	  sigm_e[i][j] = his->GetFunction("gaus")->GetParError(2)*1000;
    
    can_vs[i]->cd(j+1);
    if(i==0) his2->SetTitle("K (beam)");
    if(i==1) his2->SetTitle("#pi (prompt)");
    if(i==2) his2->SetTitle("#Sigma ");
    if(i==3) his2->SetTitle("missing n");
    if(i==4) his2->SetTitle("n from #Sigma");
    if(i==5) his2->SetTitle("#pi from #Sigma");
    if(k==0){ 
      his2->SetXTitle("px (reco.-gen.) [GeV/c]");
      his2->SetYTitle("px (gen.) [GeV/c]");
    }
    if(k==1){ 
      his2->SetXTitle("py (reco.-gen.) [GeV/c]");
      his2->SetYTitle("py (gen.) [GeV/c]");
    }
    if(k==2){ 
      his2->SetXTitle("pz (reco.-gen.) [GeV/c]");
      his2->SetYTitle("pz (gen.) [GeV/c]");
    }
    if(k==3){ 
      his2->SetXTitle("E (reco.-gen.)  [GeV]");
      his2->SetYTitle("mom. (gen.) [GeV/c]");
    }
    his2->GetXaxis()->CenterTitle();
    his2->GetYaxis()->CenterTitle();
    
    if(i==0){ 
      his2->GetXaxis()->SetRangeUser(-0.1,0.1);
      his2->GetYaxis()->SetRangeUser(0,1.2);
    }
    if(i==1){ 
      his2->GetXaxis()->SetRangeUser(-0.2,0.2);
      his2->GetYaxis()->SetRangeUser(0,1.0);
    }
    if(i==2){ 
      his2->GetXaxis()->SetRangeUser(-0.3,0.3);
      his2->GetYaxis()->SetRangeUser(0,1.2);
    }
    if(i==3){ 
      his2->GetXaxis()->SetRangeUser(-0.3,0.3);
      his2->GetYaxis()->SetRangeUser(0,1.2);
    }
    if(i==4){ 
      his2->GetXaxis()->SetRangeUser(-0.2,0.2);
      his2->GetYaxis()->SetRangeUser(0,1.2);
    }
    if(i==5){ 
      his2->GetXaxis()->SetRangeUser(-0.2,0.2);
      his2->GetYaxis()->SetRangeUser(0,0.6);
    }
    his2->Draw("colz");
	  n++;
	}else{
	  val[i][n] = 0;
	  n++;
	}
      }
    }
  }

  can[0]->Print(outfilename+"(");
  for( int i=1; i<num-1; i++ ){
    can[i]->Print(outfilename);
  }
  can[num-1]->Print(outfilename+")");
  
  cout << "cov. matrix calculated from standard deviation of histogram" << endl;
  for( int i=0; i<num; i++ ){
    cout<<"   { "<<val[i][0]<<", "<<val[i][1]<<", "<<val[i][2]<<", "<<val[i][3]<<","<<endl;
    cout<<"     "<<val[i][4]<<", "<<val[i][5]<<", "<<val[i][6]<<", "<<val[i][7]<<","<<endl;
    cout<<"     "<<val[i][8]<<", "<<val[i][9]<<", "<<val[i][10]<<", "<<val[i][11]<<","<<endl;
    cout<<"     "<<val[i][12]<<", "<<val[i][13]<<", "<<val[i][14]<<", "<<val[i][15]<<" },"<<endl;
  }
  
  cout << "cov. matrix calculated from fit" << endl;
  for( int i=0; i<num; i++ ){
    cout<<"   { "<<valfit[i][0]<<", "<<valfit[i][1]<<", "<<valfit[i][2]<<", "<<valfit[i][3]<<","<<endl;
    cout<<"     "<<valfit[i][4]<<", "<<valfit[i][5]<<", "<<valfit[i][6]<<", "<<valfit[i][7]<<","<<endl;
    cout<<"     "<<valfit[i][8]<<", "<<valfit[i][9]<<", "<<valfit[i][10]<<", "<<valfit[i][11]<<","<<endl;
    cout<<"     "<<valfit[i][12]<<", "<<valfit[i][13]<<", "<<valfit[i][14]<<", "<<valfit[i][15]<<" },"<<endl;
  }
  
  cout<<" === pole [MeV] === "<<endl;
  for( int i=0; i<num; i++ ){
    cout<<"{ "<<pole[i][0]<<", "<<pole[i][1]<<", "<<pole[i][2]<<", "<<pole[i][3]<<" },"<<endl;
  }
  cout<<" === sigma [MeV] === "<<endl;
  for( int i=0; i<num; i++ ){
    cout<<"{ "<<sigm[i][0]<<", "<<sigm[i][1]<<", "<<sigm[i][2]<<", "<<sigm[i][3]<<" },"<<endl;
  }

  double xx[4] = {1, 2, 3, 4};
  TGraphErrors *gr1[num];
  TGraphErrors *gr2[num];
  for( int i=0; i<num; i++ ){
    gr1[i] = new TGraphErrors( 4, xx, pole[i], 0, pole_e[i] );
    gr2[i] = new TGraphErrors( 4, xx, sigm[i], 0, sigm_e[i] );
  }

  TLegend *leg;

  char *com[num] = {"beam K", "#pi^{#mp}", "#Sigma^{#pm}", "n_{missing}", "n from #Sigma", "#pi^{#pm} from #Sigma"};

  TCanvas *u1;
  u1 = new TCanvas("u1", "", 600, 600);
  TH2F *h1 = new TH2F("h1", "h1", 4, 0.5, 4.5, 100, -20, 20);
  h1->Draw();
  h1->SetStats(0);
  h1->SetXTitle("px,py,pz,E");
  h1->SetYTitle("measured-generated center [MeV]");
  leg = new TLegend(0.6, 0.6, 0.9, 0.9);
  for( int i=0; i<num; i++ ){
    gr1[i]->Draw("same");
    gr1[i]->SetLineColor(i+1);
    leg->AddEntry(gr1[i], com[i], "L");
  }
  leg->Draw();
  gPad->SetLeftMargin(0.13);
  gPad->SetRightMargin(0.03);
  outfilename.Insert(outfilename.Sizeof()-5,"2");
  u1->Print(outfilename+"(");

  TCanvas *u2;
  u2 = new TCanvas("u2", "", 600, 600);
  TH2F *h2 = new TH2F("h2", "h2", 4, 0.5, 4.5, 100, 0, 40);
  h2->Draw();
  h2->SetStats(0);
  h2->SetXTitle("px,py,pz,E");
  h2->SetYTitle("measured-generated sigma [MeV]");
  leg = new TLegend(0.6, 0.6, 0.9, 0.9);
  for( int i=0; i<num; i++ ){
    gr2[i]->Draw("same");
    gr2[i]->SetLineColor(i+1);
    leg->AddEntry(gr2[i], com[i], "L");
  }
  leg->Draw();
  gPad->SetLeftMargin(0.13);
  gPad->SetRightMargin(0.03);
  u2->Print(outfilename+")");

#endif

}
