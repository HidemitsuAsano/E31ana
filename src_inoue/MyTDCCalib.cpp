#include "MyTDCCalib.h"

using namespace std;

const int omit=1;

int TDCCalib(ConfMan *conf, int cr, int sl, int ch, double period, TFile *of)
{
  if( !conf->GetGainMapManager() ){ return 2; }
  TH1F *h1=(TH1F*)gFile->Get(Form("TKOc%dn%da%d", cr, sl, ch));
  if( !h1 ){ return 2; }

  TCanvas *c1 = new TCanvas();
  c1-> Divide(2, 2);
  c1-> cd(1);
  h1-> SetLineColor(kBlack);
  double offset, gain;
  conf->GetGainMapManager()->GetParam(cr, sl, ch, gain, offset);
  TLine line;
  int npeak=0;
  double h1_max=0;
  double peakpos[100];
  double fitmean[100];
  double fitsigma[100];
  double fiterror[100];
  double tdiff[100];
  double resi[100];

  for( int i=0; i<100; i++ ) tdiff[i]=(i+1)*period;

  for( int bin=10; bin<=h1->GetNbinsX()-10; bin++ ){
    if( h1_max<h1->GetBinContent(bin) ) h1_max=h1->GetBinContent(bin);
  }
  h1-> GetYaxis()-> SetRangeUser(0, 1.05*h1_max);
  h1-> Draw();

  c1-> cd(3);
  h1-> Draw();
  line.SetLineColor(kBlue);
  for( int bin=10; bin<=h1->GetNbinsX()-10; bin++ ){
    if( 0.3*h1_max<h1->GetBinContent(bin) ){
      peakpos[npeak]=h1->GetBinCenter(bin);
      line.DrawLine(peakpos[npeak], 0, peakpos[npeak], 1.05*h1_max);
      npeak++;
      bin+=5;
    }
  }

  line.SetLineColor(kRed);
  for( int index=0; index<npeak-2*omit; index++ ){
    TF1 *gaus=(TF1*)gROOT->GetFunction("gaus");
    gaus-> SetParameter(0, h1_max);
    gaus-> SetParameter(1, peakpos[index+omit]);
    gaus-> SetParameter(2, 0.5);
    h1-> Fit(gaus, "0q", "0q", peakpos[index+omit]-5, peakpos[index+omit]+5);
    fitmean[index]=gaus->GetParameter(1);
    fitsigma[index]=gaus->GetParameter(2);
    fiterror[index]=gaus->GetParError(1);

    line.DrawLine(fitmean[index], 0, fitmean[index], 1.05*h1_max);
  }

  if( npeak<4 ) return 2;

  c1-> cd(2);
  TGraphErrors *gra = new TGraphErrors(npeak-2*omit, tdiff, fitmean, 0, fiterror);
  gra-> SetMarkerStyle(8);
  TF1 *pol1=(TF1*)gROOT->GetFunction("pol1");
  pol1-> SetParameter(0, 0);
  pol1-> SetParameter(1, 1.0/fabs(gain));
  gra-> Fit(pol1, "q", "q");
  gra-> Draw("AP");
  double tgain=1.0/gra->GetFunction("pol1")->GetParameter(1);

  c1-> cd(4);
  for( int index=0; index<npeak-2*omit; index++ ){
    resi[index]=fitmean[index]-pol1->Eval(tdiff[index]);
  }
  TGraph *gra2 = new TGraph(npeak-2*omit, tdiff, resi);
  gra2-> SetMarkerStyle(8);
  gra2-> GetYaxis()-> SetRangeUser(-5.0, 5.0);
  gra2-> Fit("pol2", "q", "q");
  gra2-> Draw("AP");

  of-> cd();
  c1-> Update();
  c1-> Write(Form("can_TKOc%dn%da%d", cr, sl, ch));

  if( gain<0 ) tgain*=-1.0;
  conf->GetGainMapManager()->SetParam(cr, sl, ch, tgain, offset);

  // string in;
  // cin>>in;
  // delete c1;

  return 0;
}
