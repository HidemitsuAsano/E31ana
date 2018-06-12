#include "GlobalVariables.h"
#include "ConfMan.h"
#include "MyCalibADC.h"
#include <TRint.h>

using namespace std;

//#define APP 1

const double calcTOF=25.9334;

int main(int argc, char** argv)
{
#if APP
  int argc2=0;
  char **argv2;
  TRint *theApp = new TRint("theApp", &argc2, argv2);
#endif
  int runnum=0;
  string confname;
  string infile;
  string outroot;
  string outparam;

  if( argc==5 ){
    confname=argv[1];
    outroot=argv[2];
    outparam=argv[3];
    infile=argv[4];
  }
  if( outparam.empty() ){ cout<<"Please input $(ConfFile) $(OutROOTFile) $(OutParamFile) $(InROOTFile)"<<endl; return 0; }
  ConfMan *conf=new ConfMan(confname.c_str(), runnum);
  conf-> Initialize(false, false);

  cout<<"> ConfFile    : "<<confname<<endl;
  cout<<"> OutROOTFile : "<<outroot<<endl;
  cout<<"> OutParam    : "<<outparam<<endl;
  cout<<"> InROOTFile  : "<<infile<<endl;

  TFile *f = new TFile(infile.c_str());
  TFile *of = new TFile(outroot.c_str(), "recreate");

  TCanvas *c1 = new TCanvas("c1", "c1");
  c1-> Divide(2, 3);

  double xxx[5];
  double offset[5];
  double sigma[5];
  int index=0;

  double xxx2[5];
  double offset2[5];
  double sigma2[5];
  int index2=0;
    
  for( int seg=1; seg<=5; seg++ ){
    c1-> cd(seg);
    if( seg==1 ) continue;
    TH1F *h1 = (TH1F*)f->Get(Form("BHDT0_BHD9_T0%d_offset", seg));
    TF1 *gaus=MyCalibTools::fitGaus_3sigma(h1);
    if( !gaus ) continue;
    h1->GetXaxis()-> SetRangeUser(-5, 5);
    offset[index]=gaus->GetParameter(1);
    sigma[index]=gaus->GetParameter(2);
    xxx[index]=seg;
    MyCalibTools::setTimeOffset(conf, CID_T0, seg, -offset[index]);
    index++;
  }

  for( int seg=1; seg<=5; seg++ ){
    c1-> cd(seg);
    TH1F *h1 = (TH1F*)f->Get(Form("BHDT0_BHD9_T0%d_offset_pi", seg));
    TF1 *gaus=MyCalibTools::fitGaus_3sigma(h1);
    if( !gaus ) continue;
    h1->GetXaxis()-> SetRangeUser(-5, 5);
    offset2[index2]=gaus->GetParameter(1);
    sigma2[index2]=gaus->GetParameter(2);
    xxx2[index2]=seg;
    if( seg==1 ) MyCalibTools::setTimeOffset(conf, CID_T0, seg, -offset2[index2]);
    index2++;
  }
  c1-> Clear();
  c1-> Divide(1, 2);
  c1-> cd(1);

  TGraph *gra_offset = new TGraph(index, xxx, offset);
  gra_offset-> Write("gra_T0_offset_K");
  TGraph *gra_reso = new TGraph(index, xxx, sigma);
  gra_reso-> Write("gra_T0_reso_K");

  TGraph *gra_offset_pi = new TGraph(index2, xxx2, offset2);
  gra_offset_pi-> Write("gra_T0_offset_pi");

  TGraph *gra_reso_pi = new TGraph(index, xxx, sigma2);
  gra_reso_pi-> Write("gra_T0_reso_pi");

  // c1-> Update();
  // string in;
  // cin>>in;
  //  c1-> Write("can_T0_offset0");
  of-> Write();
  of-> Close();

  ofstream ofs(outparam.c_str());
  conf-> GetGainMapManager()-> PrintMapBL(ofs);
#if APP
  theApp-> Run();
  delete theApp;
#endif
  return 0;
}
