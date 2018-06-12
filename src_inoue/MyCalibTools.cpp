#include "MyCalibTools.h"

using namespace std;

bool MyCalibTools::setTimeOffset(ConfMan *conf, int cid, int lay, int wire, double offset)
{
  int cr, sl, ch;
  double gain, old_offset;
  if( !conf->GetCounterMapManager()->GetCNA(cid, lay, wire, 1, 0, cr, sl, ch) ){
    cout<<"  !!! cid="<<cid<<" layer="<<lay<<" wire="<<wire<<" GetCNA fault !!!"<<endl;
    return false;
  }
  conf-> GetGainMapManager()-> GetParam(cr, sl, ch, gain, old_offset);
  conf-> GetGainMapManager()-> SetParam(cr, sl, ch, gain, old_offset+offset);
  return true;
}

bool MyCalibTools::setTimeOffset(ConfMan *conf, int cid, int seg, double offset)
{
  int cr, sl, ch;
  double gain, old_offset;
  if( !conf->GetCounterMapManager()->GetCNA(cid, seg, 1, 0, cr, sl, ch) ){
    cout<<"  !!! cid="<<cid<<" seg="<<seg<<" up   GetCNA fault !!!"<<endl;
    return false;
  }
  conf-> GetGainMapManager()-> GetParam(cr, sl, ch, gain, old_offset);
  conf-> GetGainMapManager()-> SetParam(cr, sl, ch, gain, old_offset+offset);

  if( !conf->GetCounterMapManager()->GetCNA(cid, seg, 1, 1, cr, sl, ch) ){
    cout<<"  !!! cid="<<cid<<" seg="<<seg<<" down GetCNA fault !!!"<<endl;
    return false;
  }
  conf-> GetGainMapManager()-> GetParam(cr, sl, ch, gain, old_offset);
  conf-> GetGainMapManager()-> SetParam(cr, sl, ch, gain, old_offset+offset);
  return true;
}

TF1 *MyCalibTools::fitLandau(TH1F *h1, double min)
{
  double mpv=min;
  double max_val=0;
  for( int bin=1; bin<=h1->GetNbinsX(); bin++ ){
    if( h1->GetBinCenter(bin)<min ) continue;
    if( max_val<h1->GetBinContent(bin) ){
      max_val=h1->GetBinContent(bin);
      mpv=h1->GetBinCenter(bin);
    }
  }

  int count=0;
  double sigma=(mpv-min)/5.0;
  while( count<1000 ){
    h1-> Fit("landau", "0q", "0q", mpv-5*sigma, mpv+20*sigma);

    if( fabs(mpv-h1->GetFunction("landau")->GetParameter(1))<1.0e-2 && fabs(sigma-h1->GetFunction("landau")->GetParameter(2))<1.0e-2 ){
      break;
    }

    mpv=h1->GetFunction("landau")->GetParameter(1);
    sigma=h1->GetFunction("landau")->GetParameter(2);
    if( sigma<0 ) sigma*=-1;

    count++;
    //    cout<<mpv<<"  "<<sigma<<endl;
  }
  if( count==1000 ) cout<<"    ! MyCalibTools::fitLandau count==1000 !"<<endl;

  return h1->GetFunction("landau");
}

TF1 *MyCalibTools::fitGaus_3sigma(TH1F *h1)
{
  double mean=h1->GetBinCenter(h1->GetMaximumBin());
  double sigma=3*h1->GetBinWidth(h1->GetMaximumBin());

  return MyCalibTools::fitGaus_3sigma(h1, mean, sigma);
}


TF1 *MyCalibTools::fitGaus_3sigma(TH1F *h1, double mean, double sigma)
{
  if( h1->GetEntries()<0 ){
    cout<<"  !!! MyCalibTools::fitGaus_3sigma("<<h1->GetName()<<")  empty !!!"<<endl;
    return 0;
  }

  TF1 *gaus=(TF1*)gROOT->GetFunction("gaus");
  gaus-> SetParameter(0, h1->GetMaximum());
  gaus-> SetParameter(1, mean);
  gaus-> SetParameter(2, sigma);

  int count=0;
  while( count<1000 ){
    h1-> Fit("gaus", "0wq", "0wq", mean-3*sigma, mean+3*sigma);
    if( !h1->GetFunction("gaus") ){
      double mean=h1->GetBinCenter(h1->GetMaximumBin());
      h1-> Fit("gaus", "0wq", "0wq", mean-3*sigma, mean+3*sigma);
    }
    if( fabs(mean-h1->GetFunction("gaus")->GetParameter(1))<1.0e-7 && fabs(sigma-h1->GetFunction("gaus")->GetParameter(2))<1.0e-7 ) break;
    mean=h1->GetFunction("gaus")-> GetParameter(1);
    sigma=h1->GetFunction("gaus")-> GetParameter(2);

    count++;
  }
  if( count==1000 ) cout<<"    ! MyCalibTools::fitGaus_3sigma count==1000 !"<<endl;
  return h1->GetFunction("gaus");
}

TF1 *MyCalibTools::fitGaus_2sigma(TH1F *h1)
{
  double mean=h1->GetBinCenter(h1->GetMaximumBin());
  double sigma=2*h1->GetBinWidth(h1->GetMaximumBin());

  return MyCalibTools::fitGaus_2sigma(h1, mean, sigma);
}


TF1 *MyCalibTools::fitGaus_2sigma(TH1F *h1, double mean, double sigma)
{
  if( h1->GetEntries()<0 ){
    cout<<"  !!! MyCalibTools::fitGaus_2sigma("<<h1->GetName()<<")  empty !!!"<<endl;
    return 0;
  }

  TF1 *gaus=(TF1*)gROOT->GetFunction("gaus");
  gaus-> SetParameter(0, h1->GetMaximum());
  gaus-> SetParameter(1, mean);
  gaus-> SetParameter(2, sigma);

  int count=0;
  while( count<1000 ){
    h1-> Fit("gaus", "0wq", "0wq", mean-2*sigma, mean+2*sigma);
    if( !h1->GetFunction("gaus") ){
      double mean=h1->GetBinCenter(h1->GetMaximumBin());
      h1-> Fit("gaus", "0wq", "0wq", mean-2*sigma, mean+2*sigma);
    }

    if( fabs(mean-h1->GetFunction("gaus")->GetParameter(1))<1.0e-7 && fabs(sigma-h1->GetFunction("gaus")->GetParameter(2))<1.0e-7 ) break;
    mean=h1->GetFunction("gaus")-> GetParameter(1);
    sigma=h1->GetFunction("gaus")-> GetParameter(2);

    count++;
  }
  if( count==1000 ) cout<<"    ! MyCalibTools::fitGaus_2sigma count==1000 !"<<endl;
  return h1->GetFunction("gaus");
}

TF1 *MyCalibTools::fitGaus_25sigma(TH1F *h1)
{
  double mean=h1->GetBinCenter(h1->GetMaximumBin());
  double sigma=2*h1->GetBinWidth(h1->GetMaximumBin());

  return MyCalibTools::fitGaus_2sigma(h1, mean, sigma);
}


TF1 *MyCalibTools::fitGaus_25sigma(TH1F *h1, double mean, double sigma)
{
  if( h1->GetEntries()<0 ){
    cout<<"  !!! MyCalibTools::fitGaus_2sigma("<<h1->GetName()<<")  empty !!!"<<endl;
    return 0;
  }

  TF1 *gaus=(TF1*)gROOT->GetFunction("gaus");
  gaus-> SetParameter(0, h1->GetMaximum());
  gaus-> SetParameter(1, mean);
  gaus-> SetParameter(2, sigma);

  int count=0;
  while( count<1000 ){
    h1-> Fit("gaus", "0wq", "0wq", mean-2.5*sigma, mean+2.5*sigma);
    if( !h1->GetFunction("gaus") ){
      double mean=h1->GetBinCenter(h1->GetMaximumBin());
      h1-> Fit("gaus", "0wq", "0wq", mean-2.5*sigma, mean+2.5*sigma);
    }

    if( fabs(mean-h1->GetFunction("gaus")->GetParameter(1))<1.0e-7 && fabs(sigma-h1->GetFunction("gaus")->GetParameter(2))<1.0e-7 ) break;
    mean=h1->GetFunction("gaus")-> GetParameter(1);
    sigma=h1->GetFunction("gaus")-> GetParameter(2);

    count++;
  }
  if( count==1000 ) cout<<"    ! MyCalibTools::fitGaus_2sigma count==1000 !"<<endl;
  return h1->GetFunction("gaus");
}
