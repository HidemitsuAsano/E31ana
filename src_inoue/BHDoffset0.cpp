#include "GlobalVariables.h"
#include "ConfMan.h"
#include "MyCalibADC.h"
#include <TRint.h>

using namespace std;

#define APP 1

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
  c1-> Divide(2, 5);
  int index=0;

  for( int seg=6; seg<=15; seg++ ){
    c1-> cd(seg-5);
    TH1F *h1 = (TH1F*)f->Get(Form("BHDT0_pi_BHD%d_T03", seg));
    TF1 *gaus=MyCalibTools::fitGaus_3sigma(h1);
    h1->GetXaxis()-> SetRangeUser(20, 30);
    h1-> Draw();
    gaus-> Draw("same");
    double offset=calcTOF-gaus->GetParameter(1);
    MyCalibTools::setTimeOffset(conf, CID_BHD, seg, offset);
  }

  c1-> Update();
  string in;
  cin>>in;

  c1-> Write("can_BHD_offset0");
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
