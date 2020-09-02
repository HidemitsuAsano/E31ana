//
// for the covariance matrix evaluation
//
{
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(111111);

  TFile *f = new TFile("tmp.root");
  //TFile *f = new TFile("simdata012_ana00/sim_0000.root");

  const int num = 7;

  TCanvas *can[num];
  for( int i=0; i<num; i++ ){
    TString str(Form("can%d", i));
    can[i] = new TCanvas(str, "", 800, 800);
    can[i]->Divide(2, 2);
  }

  // beam_K(K+), L, n, p, p from L, pi- from L
  double val[num][16];
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
	  his->Fit("gaus");
	  his->GetFunction("gaus")->SetLineColor(4);
	  double pol = his->GetFunction("gaus")->GetParameter(1);
	  double sig = his->GetFunction("gaus")->GetParameter(2);
	  his->Fit("gaus","","",pol-3*sig,pol+3*sig);
	  pol = his->GetFunction("gaus")->GetParameter(1);
	  sig = his->GetFunction("gaus")->GetParameter(2);
	  val[i][n] = sig*sig;
	  pole[i][j] = pol*1000;
	  sigm[i][j] = sig*1000;
	  pole_e[i][j] = his->GetFunction("gaus")->GetParError(1)*1000;
	  sigm_e[i][j] = his->GetFunction("gaus")->GetParError(2)*1000;
	  n++;
	}else{
	  val[i][n] = 0;
	  n++;
	}
      }
    }
  }

  can[0]->Print("tmp.pdf(");
  for( int i=1; i<num-1; i++ ){
    can[i]->Print("tmp.pdf");
  }
  can[num-1]->Print("tmp.pdf)");

  for( int i=0; i<num; i++ ){
    cerr<<"   { "<<val[i][0]<<", "<<val[i][1]<<", "<<val[i][2]<<", "<<val[i][3]<<","<<endl;
    cerr<<"     "<<val[i][4]<<", "<<val[i][5]<<", "<<val[i][6]<<", "<<val[i][7]<<","<<endl;
    cerr<<"     "<<val[i][8]<<", "<<val[i][9]<<", "<<val[i][10]<<", "<<val[i][11]<<","<<endl;
    cerr<<"     "<<val[i][12]<<", "<<val[i][13]<<", "<<val[i][14]<<", "<<val[i][15]<<" },"<<endl;
  }
  
#if 1
  cerr<<" === pole [MeV] === "<<endl;
  for( int i=0; i<num; i++ ){
    cerr<<"{ "<<pole[i][0]<<", "<<pole[i][1]<<", "<<pole[i][2]<<", "<<pole[i][3]<<" },"<<endl;
  }
  cerr<<" === sigma [MeV] === "<<endl;
  for( int i=0; i<num; i++ ){
    cerr<<"{ "<<sigm[i][0]<<", "<<sigm[i][1]<<", "<<sigm[i][2]<<", "<<sigm[i][3]<<" },"<<endl;
  }

  double xx[4] = {1, 2, 3, 4};
  TGraphErrors *gr1[num];
  TGraphErrors *gr2[num];
  for( int i=0; i<num; i++ ){
    gr1[i] = new TGraphErrors( 4, xx, pole[i], 0, pole_e[i] );
    gr2[i] = new TGraphErrors( 4, xx, sigm[i], 0, sigm_e[i] );
  }

  TLegend *leg;

#define FLAG 7;

#if FLAG == 0
  char *com[num] = {"beam K", "#Lambda", "n_{missing}", "p", "p from #Lambda", "#pi^{-} from #Lambda"};
#elif FLAG == 1
  char *com[num] = {"beam K", "#pi^{#mp}", "#Sigma^{#pm}", "p", "n_{missing}", "n from #Sigma", "#pi^{#pm} from #Sigma"};
#elif FLAG == 2
  char *com[num] = {"beam K", "#pi^{+}", "#Lambda", "n", "n_{missing}", "p from #Lambda", "#pi^{-} from #Lambda"};
#elif FLAG == 3
  char *com[num] = {"beam K", "#pi^{-}", "#Lambda", "p", "p_{missing}", "p from #Lambda", "#pi^{-} from #Lambda"};
#elif FLAG == 4
  char *com[num] = {"beam K", "K^{-}", "p", "p", "n_{missing}"};
#elif FLAG == 5
  char *com[num] = {"beam K", "#Sigma^{#pm}", "#pi^{#mp}_{missing}", "n from #Sigma", "#pi^{#pm} from #Sigma"};
#elif FLAG == 6
  char *com[num] = {"beam K", "K^{0}", "#pi^{+}", "#pi^{-}", "missing n"};
#elif FLAG == 7
  char *com[num] = {"beam K", "#Lambda", "p", "n", "#pi_{missing}", "p from #Lambda", "#pi^{-} from #Lambda"};
#endif

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
  u1->Print("tmp2.pdf(");

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
  u2->Print("tmp2.pdf)");

#endif

}
