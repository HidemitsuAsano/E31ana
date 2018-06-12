#include "GlobalVariables.h"
#include "ConfMan.h"
#include "MyBLDCCalib.h"
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

  TString chmName[]={ "FDC1" };
  int cid[]={ CID_FDC1 };
  for( int i=0; i<sizeof(chmName)/sizeof(chmName[0]); i++ ){
    int nlay=DetectorList::GetInstance()->GetNlayers(cid[i]);
    for( int lay=1; lay<=nlay; lay++ ){
      int nwire=DetectorList::GetInstance()->GetNwires(cid[i]);
      for( int wire=1; wire<=nwire; wire++ ){
	f-> cd();
	Chamber_offset0(conf, of, cid[i], lay, wire);
      }
    }
  }

  ofstream ofs(outparam.c_str());
  conf-> GetGainMapManager()-> PrintMapBL(ofs);

  of-> Write();
  of-> Close();

#if APP
  theApp-> Run();
  delete theApp;
#endif
  return 0;
}
