#include "MyCalibADC.h"

using namespace std;

std::vector<double> CalibADC(ConfMan *conf, TFile *of, int cid, int seg, double dE)
{
  string name;
  if( dE<=0 ) exit(0);
  if( cid==CID_BHD      ) name="BHD";
  else if( cid==CID_T0  ) name="T0";
  else if( cid==CID_BVC ) name="BVC";
  else if( cid==CID_NC  ) name="NC";
  else if( cid==CID_PC  ) name="PC";
  else if( cid==CID_CVC ) name="CVC";
  else{ cout<<"!!!!! Invailed cid="<<cid<<" !!!!!"<<endl; exit(0); }

  TLine line;
  TH1F *h1_0 = (TH1F*)gFile-> Get(Form("A%s_u%d", name.c_str(), seg));
  TH1F *h1_0t = (TH1F*)gFile-> Get(Form("A%s_u%d_wt", name.c_str(), seg));
  TH1F *h1_1 = (TH1F*)gFile-> Get(Form("A%s_d%d", name.c_str(), seg));
  TH1F *h1_1t = (TH1F*)gFile-> Get(Form("A%s_d%d_wt", name.c_str(), seg));
  h1_0t-> SetLineColor(kRed);
  h1_1t-> SetLineColor(kRed);

  TCanvas *c1 =new TCanvas("c1", "c1");
  c1-> Divide(1, 2);
  c1-> cd(1);

  TF1 *gaus0=MyCalibTools::fitGaus_3sigma(h1_0);
  double pede0=gaus0->GetParameter(1);
  double min=gaus0->GetParameter(1)+5*gaus0->GetParameter(2);
  TF1 *landau0=MyCalibTools::fitLandau(h1_0t, min);
  double mpv0=landau0-> GetParameter(1);
  double gain0=dE/(mpv0-pede0);
  double thr0=0;
  double thre0=0;
  for( int bin=1; bin<=h1_0->GetNbinsX(); bin++ ){
    if( h1_0->GetBinCenter(bin)<min || h1_0-> GetBinContent(bin)==0 ) continue;
    if( h1_0t->GetBinContent(bin)/h1_0->GetBinContent(bin)>0.2 ){
      thre0=h1_0->GetBinCenter(bin);
      thr0=(h1_0->GetBinCenter(bin)-pede0)/(mpv0-pede0);
      break;
    }
  }

  h1_0-> GetXaxis()-> SetRangeUser(0, 2000);
  h1_0-> Draw();
  gPad-> SetLogy();
  gaus0-> Draw("same");
  h1_0t-> Draw("same");
  landau0-> SetLineColor(kBlue);
  landau0-> Draw("same");
  line.SetLineColor(kRed);
  line.DrawLine(pede0, 0, pede0, 2.0*h1_0->GetMaximum());
  line.SetLineColor(kBlue);
  line.DrawLine(mpv0, 0, mpv0, 2.0*h1_0->GetMaximum());
  line.SetLineColor(3);
  line.DrawLine(thre0, 0, thre0, 2.0*h1_0->GetMaximum());

  c1-> cd(2);
  TF1 *gaus1=MyCalibTools::fitGaus_3sigma(h1_1);
  double pede1=gaus1->GetParameter(1);
  min=gaus1->GetParameter(1)+5*gaus1->GetParameter(2);
  TF1 *landau1=MyCalibTools::fitLandau(h1_1t, min);
  double mpv1=landau1-> GetParameter(1);
  double gain1=dE/(mpv1-pede1);
  double thr1=0;
  double thre1=0;
  for( int bin=1; bin<=h1_1->GetNbinsX(); bin++ ){
    if( h1_1->GetBinCenter(bin)<min || h1_1-> GetBinContent(bin)==0 ) continue;
    if( h1_1t->GetBinContent(bin)/h1_1->GetBinContent(bin)>0.2 ){
      thre1=h1_0->GetBinCenter(bin);
      thr1=(h1_1->GetBinCenter(bin)-pede1)/(mpv1-pede1);
      break;
    }
  }

  h1_1-> GetXaxis()-> SetRangeUser(0, 2000);
  h1_1-> Draw();
  gPad-> SetLogy();
  gaus1-> Draw("same");
  h1_1t-> Draw("same");
  landau1-> SetLineColor(kBlue);
  landau1-> Draw("same");
  line.SetLineColor(kRed);
  line.DrawLine(pede1, 0, pede1, 2.0*h1_1->GetMaximum());
  line.SetLineColor(kBlue);
  line.DrawLine(mpv1, 0, mpv1, 2.0*h1_1->GetMaximum());
  line.SetLineColor(3);
  line.DrawLine(thre1, 0, thre1, 2.0*h1_1->GetMaximum());

  // c1-> Update();
  // string in;
  // cin>>in;

  of-> cd();
  c1-> Write(Form("%s_%d_ADC", name.c_str(), seg));
  delete c1;

  vector<double> param;
  int cr, sl, ch;
  conf->GetCounterMapManager()-> GetCNA(cid, seg, 0, 0, cr, sl, ch);
  conf->GetGainMapManager()-> SetParam(cr, sl, ch, gain0, pede0);

  conf->GetCounterMapManager()-> GetCNA(cid, seg, 0, 1, cr, sl, ch);
  conf->GetGainMapManager()-> SetParam(cr, sl, ch, gain1, pede1);

  param.push_back(pede0);
  param.push_back(gain0);
  param.push_back(thr0);
  param.push_back(pede1);
  param.push_back(gain1);
  param.push_back(thr1);
  return param;
}
