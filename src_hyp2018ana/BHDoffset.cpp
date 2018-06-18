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

  double xxx[20];
  double offset[20];
  double sigma[20];
  int index=0;
  for( int seg=6; seg<=15; seg++ ){
    TH1F *h1 = (TH1F*)f->Get(Form("BHDT0_BHD%d_T03_offset", seg));
    TF1 *gaus=MyCalibTools::fitGaus_3sigma(h1);
    h1->GetXaxis()-> SetRangeUser(20, 30);

    xxx[index]=seg;
    offset[index]=gaus->GetParameter(1);
    sigma[index]=gaus->GetParameter(2);
    MyCalibTools::setTimeOffset(conf, CID_BHD, seg, offset[index]);
    index++;
  }

  TGraph *gra = new TGraph(index, xxx, offset);
  gra-> Write("gra_BHD_offset");

  TGraph *gra2 = new TGraph(index, xxx, sigma);
  gra2-> Write("gra_BHD_reso");
  // c1-> Update();
  // string in;
  // cin>>in;

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
