#include "GlobalVariables.h"
#include "ConfMan.h"
#include "MyCalibADC.h"
#include <TRint.h>

using namespace std;

//#define APP 1

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

  int index=0;
  double xxx[100];
  double pede_u[100];
  double pede_d[100];
  double gain_u[100];
  double gain_d[100];
  double thre_u[100];
  double thre_d[100];
  for( int seg=3; seg<=18; seg++ ){
    f-> cd();
    vector<double> adc_param=CalibADC(conf, of, CID_BHD, seg, 2.0);
    xxx[index]=seg;
    pede_u[index]=adc_param[0];
    gain_u[index]=adc_param[1];
    thre_u[index]=adc_param[2];
    pede_d[index]=adc_param[3];
    gain_d[index]=adc_param[4];
    thre_d[index]=adc_param[5];
    index++;
  }

  of-> cd();
  TGraph *gra_pede_u=new TGraph(index, xxx, pede_u);
  TGraph *gra_pede_d=new TGraph(index, xxx, pede_d);
  TGraph *gra_gain_u=new TGraph(index, xxx, gain_u);
  TGraph *gra_gain_d=new TGraph(index, xxx, gain_d);
  TGraph *gra_thre_u=new TGraph(index, xxx, thre_u);
  TGraph *gra_thre_d=new TGraph(index, xxx, thre_d);

  gra_pede_u-> Write("gra_BHD_pede_up");
  gra_pede_d-> Write("gra_BHD_pede_down");
  gra_gain_u-> Write("gra_BHD_gain_up");
  gra_gain_d-> Write("gra_BHD_gain_down");
  gra_thre_u-> Write("gra_BHD_thre_up");
  gra_thre_d-> Write("gra_BHD_thre_down");

  ofstream ofs(outparam.c_str());
  conf->GetGainMapManager()-> PrintMapBL(ofs);

  of-> Write();
  of-> Close();

#if APP
  delete theApp;
#endif
  return 0;
}

