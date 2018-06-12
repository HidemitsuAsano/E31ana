#include "MyBLDCCalib.h"

using namespace std;

void Chamber_dxdt(ConfMan *conf, ostream &os)
{
  int cid[]={ CID_BLC1a, CID_BLC1b, CID_BLC2a, CID_BLC2b, CID_BPC, CID_FDC1 };
  for( int i=0; i<sizeof(cid)/sizeof(cid[0]); i++ ){
    Chamber_dxdt(conf, cid[i], os);
  }
}

void Chamber_dxdt(ConfMan *conf, int cid, ostream &os)
{
  string name;
  if( cid==CID_BLC1a ) name="BLC1a";
  else if( cid==CID_BLC1b ) name="BLC1b";
  else if( cid==CID_BLC2a ) name="BLC2a";
  else if( cid==CID_BLC2b ) name="BLC2b";
  else if( cid==CID_BPC   ) name="BPC";
  else if( cid==CID_FDC1  ) name="FDC1";
  double dl=fabs(5.0*conf->GetBLDCWireMapManager()->GetWireMap(cid, 1)->GetdXY());

  int nlay=DetectorList::GetInstance()->GetNlayers(cid);
  int nwire=DetectorList::GetInstance()->GetNwires(cid);
  double threshold[]={ 0,  0.3,  2,  5,  10,  17,  25,  32,  40,  50,  60,  68,  75,  83,  90,  95,  98,  99.7,  100 };
  os<<"XTParam: "<<cid<<" "<<Form("%3.2lf", dl)<<" "<<sizeof(threshold)/sizeof(threshold[0])<<endl;
  os<<"Threshold: ";
  for( int i=0; i<sizeof(threshold)/sizeof(threshold[0]); i++ ) os<<Form("%3.2lf", threshold[i])<<" ";
  os<<endl;
  for( int l=1; l<=nlay; l++ ){
    for( int w=1; w<=nwire; w++ ){
      TH1F *h1=(TH1F*)gFile->Get(Form("%s_dt_%d_%d_1hit", name.c_str(), l, w));
      if( !h1 ){ cout<<"!!!!! histogram not found !!!!!"<<endl; exit(0); }
      if( h1->GetEntries()==0 ) cout<<"  !!! "<<name<<" lay"<<l<<" wire"<<w<<" no data !!!"<<endl;

      TH1F *h1_int=(TH1F*)h1->Clone(Form("xt_%s_%d_%d", name.c_str(), l, w));
      double sum=0;
      for( int bin=1; bin<=h1->GetNbinsX(); bin++ ){
	sum+=h1->GetBinContent(bin);
	h1_int-> SetBinContent(bin, sum);
      }
      h1_int-> Scale(1.0/h1_int->GetMaximum());

      double val[sizeof(threshold)/sizeof(threshold[0])];
      for( int i=0; i<sizeof(threshold)/sizeof(threshold[0]); i++ ){
	for( int bin=1; bin<=h1->GetNbinsX(); bin++ ){
	  if( h1_int->GetBinContent(bin)>=0.01*threshold[i] ){
	    val[i]=h1_int->GetBinCenter(bin);
	    break;
	  }
	}
      }
      h1_int-> Scale(dl);
      os<<setw(5)<<cid<<" "<<setw(4)<<l<<" "<<setw(4)<<w;
      for( int i=0; i<sizeof(threshold)/sizeof(threshold[0]); i++ ) os<<" "<<setw(7)<<Form("%3.2lf", val[i]);
      os<<endl;
    }
  }
}

void Chamber_offset0(ConfMan *conf, TFile *of, int cid, int lay, int wire)
{
  string name;
  if( cid==CID_BLC1a ) name="BLC1a";
  else if( cid==CID_BLC1b ) name="BLC1b";
  else if( cid==CID_BLC2a ) name="BLC2a";
  else if( cid==CID_BLC2b ) name="BLC2b";
  else if( cid==CID_BPC   ) name="BPC";
  else if( cid==CID_FDC1  ) name="FDC1";

  TH1F *h1=0;
  if( wire>0 ) h1=(TH1F*)gFile->Get(Form("%s_dt_%d_%d_1hit", name.c_str(), lay, wire));
  else h1=(TH1F*)gFile->Get(Form("%s_dt_%d_1hit", name.c_str(), lay));
  if( !h1 ){ cout<<"!!!!! histogram not found !!!!!"<<endl; exit(0); }
  if( h1->GetEntries()==0 ){
    cout<<"  !!! "<<name<<" lay"<<lay<<" wire"<<wire<<" no data !!!"<<endl;
    return;
  }

  TH1F *h1_int=(TH1F*)h1->Clone(Form("%s_int_%d_%d", name.c_str(), lay, wire));
  double sum=0;
  for( int bin=1; bin<=h1->GetNbinsX(); bin++ ){
    sum+=h1->GetBinContent(bin);
    h1_int-> SetBinContent(bin, sum);
  }

  TH1F *h1_c=(TH1F*)h1->Clone(Form("%s_clone", h1->GetName()));

  // if( cid==CID_BLC1a && lay==1 && wire==11 ) h1_c-> Rebin(5);
  // else if( h1->GetEntries()<2000 ) h1_c-> Rebin(5);
  // else h1_c-> Rebin(2);

  if( h1->GetEntries()<2000 ) h1_c-> Rebin(10);
  else   h1_c-> Rebin(4);

  TH1F *h1_diff=(TH1F*)h1_c->Clone(Form("%s_diff_%d_%d", name.c_str(), lay, wire));
  for( int bin=1; bin<=h1_c->GetNbinsX(); bin++ ){
    h1_diff-> SetBinContent(bin, h1_c->GetBinContent(bin)-h1_c->GetBinContent(bin-1));
  }
  h1_diff-> Fit("gaus", "0q", "0q", h1_diff->GetBinCenter(h1_diff->GetMaximumBin())-14, h1_diff->GetBinCenter(h1_diff->GetMaximumBin())+14);
  TF1 *gaus=h1_diff->GetFunction("gaus");
  //  TF1 *gaus=MyCalibTools::fitGaus_25sigma(h1_diff, 0, 4);

  h1-> SetLineColor(kBlack);
  h1_int-> SetLineColor(kBlack);
  h1_diff-> SetLineColor(kBlack);

  TCanvas *c1 = new TCanvas("c1", "c1");
  c1-> Divide(1, 3);
  c1-> cd(1);
  h1-> Draw();

  c1-> cd(2);
  h1_int-> Scale(1.0/h1_int->GetMaximum());
  h1_int-> Draw();

  c1-> cd(3);

  h1_diff-> Draw();
  gaus-> Draw("same");
  TLine line;
  line.SetLineColor(kRed);
  double offset=gaus->GetParameter(1);
  line.DrawLine(offset, 1.05*h1_diff->GetMinimum(), offset, 1.05*h1_diff->GetMaximum());

  if( fabs(offset)>10 ) cout<<"  !!! "<<name<<" l"<<lay<<" w"<<wire<<"  offset="<<offset<<" !!!"<<endl;

  if( wire==0 ){
    int nwire=DetectorList::GetInstance()->GetNwires(cid);
    for( int w=1; w<=nwire; w++ ){
      MyCalibTools::setTimeOffset(conf, cid, lay, w, offset);
    }
  }

  // c1-> Update();
  // string in;
  // cin>>in;

  of-> cd();
  if( wire==0 ) c1-> Write(Form("%s_%d", name.c_str(), lay));
  else c1-> Write(Form("%s_%d_%d", name.c_str(), lay, wire));
  delete c1;

  if( cid==CID_BPC && lay==1 && wire==1 ) return;
  if( cid==CID_FDC1 && lay==3 && wire==41 ) return;
  if( cid==CID_FDC1 && lay==5 && wire==47 ) return;
  if( cid==CID_FDC1 && lay==5 && wire==52 ) return;
  else MyCalibTools::setTimeOffset(conf, cid, lay, wire, offset);

  return;
}

