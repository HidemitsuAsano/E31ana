#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <TString.h>
#include <vector>

#include "ConfMan.h"
#include "MyTDCCalib.h"
#include <TRint.h>

//#define APP 1

using namespace std;

int main(int argc, char **argv)
{
  if( argc!=2 ){
    cout<<"Please input \"./bin/tdccalib ${Conffile} \""<<endl;
    return 0;
  }
#if APP
  int argc2=0;
  char **argv2;
  TApplication *theApp=new TApplication("theApp", &argc2, argv2);
#endif

  string confname;
  string outfile;
  string paramdir;
  string logfile;

  ifstream ifs(argv[1]);
  string line;
  while( getline(ifs, line) ){
    if( line.empty() ) continue;
    if( line[0]=='#' ) continue;

    char c_str[1000];
    if( sscanf( line.c_str(), "ConfMan: %s", c_str)==1 ) confname=c_str;
    else if( sscanf( line.c_str(), "OutFile: %s", c_str)==1 ) outfile=c_str;
    else if( sscanf( line.c_str(), "ParamDir: %s", c_str)==1 ) paramdir=c_str;
    else if( sscanf( line.c_str(), "LogFile: %s", c_str)==1 ) logfile=c_str;
  }

  if( confname.empty() ){ cout<<"!!!!! Conf file not found !!!!!"<<endl; return 0; }
  if( outfile.empty() ){ cout<<"!!!!! Out file directory not found !!!!!"<<endl; return 0; }
  if( paramdir.empty() ){ cout<<"!!!!! Parameter file directory not found !!!!!"<<endl; return 0; }
  if( logfile.empty() ){ cout<<"!!!!! Parameter file directory not found !!!!!"<<endl; return 0; }

  ConfMan *conf=new ConfMan(confname.c_str(), 0);
  conf-> Initialize(false, false);
  TFile *of=new TFile(outfile.c_str(), "recreate");

  ifs.clear();
  ifs.seekg(0, ios_base::beg);

  cout<<"> ConfMan : "<<conf->GetConfFileName()<<endl;
  cout<<"> Out File : "<<outfile<<endl;
  cout<<"> Parameter Diretory : "<<paramdir<<endl;
  cout<<"> LogFile : "<<logfile<<endl;

  ofstream log(logfile.c_str());

  cout<<"TDC Calibration START"<<endl;
  while( getline(ifs, line) ){
    if( line.empty() ) continue;
    if( line[0]=='#' ) continue;
    double period;
    int cr; 
    char c_str[1000], slot[512], channel[512];
    vector<int> sl, ch;
    if( sscanf( line.c_str(), "ConfMan: %s", c_str)==1 ) continue;
    else if( sscanf( line.c_str(), "OutFile: %s", c_str)==1 ) continue;
    else if( sscanf( line.c_str(), "ParamDir: %s", c_str)==1 ) continue;
    else if( sscanf( line.c_str(), "LogFile: %s", c_str)==1 ) continue;
    else if( sscanf( line.c_str(), "%s %lf %d %s %s", c_str, &period, &cr, slot, channel)==5 ){
      int num1, num2;
      if( sscanf(channel, "%d-%d", &num1, &num2)==2 ){
	for( int i=num1; i<=num2; i++ ) ch.push_back(i);
      }

      if( sscanf(slot, "%d-%d", &num1, &num2)==2 ){
	for( int i=num1; i<=num2; i++ ) sl.push_back(i);
      }
      else if( sscanf(slot, "%d,%d", &num1, &num2)==2 ){
	sl.push_back(num1); sl.push_back(num2);
      }
      else if( sscanf(slot, "%d,%d", &num1)==1 ){
	sl.push_back(num1);;
      }
    }
    else if( sscanf( line.c_str(), "%s %lf %d %s", c_str, &period, &cr, slot)==4 ){
      for( int i=0; i<16; i++ ) ch.push_back(i);

      int num1, num2;
      if( sscanf(slot, "%d-%d", &num1, &num2)==2 ){
	for( int i=num1; i<=num2; i++ ) sl.push_back(i);
      }
      else if( sscanf(slot, "%d,%d", &num1, &num2)==2 ){
	sl.push_back(num1); sl.push_back(num2);
      }
      else if( sscanf(slot, "%d,%d", &num1)==1 ){
	sl.push_back(num1);;
      }
    }
    else{ cout<<"!!!!! Invailed file format !!!!! "<<line<<endl; return 0; }

    TFile *f = new TFile(c_str);
    for( int i=0; i<sl.size(); i++ ){
      for( int j=0; j<ch.size(); j++ ){
	f-> cd();
	int status=TDCCalib(conf, cr, sl[i], ch[j], period, of);
	if( status==1 ) log<<"  !!! TKOc"<<cr<<"n"<<sl[i]<<"a"<<ch[j]<<" fitting is strange   parameter is updated !!!"<<endl;
	if( status==2 ) log<<"!!!!! TKOc"<<cr<<"n"<<sl[i]<<"a"<<ch[j]<<" fitting fault   parameter isn't updated !!!!!"<<endl;
      }
    }
  }
  ofstream ofsBL(Form("%s/GainMapBL.param", paramdir.c_str()));
  conf-> GetGainMapManager()-> PrintMapBL(ofsBL);
  ofstream ofsCDS(Form("%s/GainMapCDS.param", paramdir.c_str()));
  conf-> GetGainMapManager()-> PrintMapCDS(ofsCDS);

  cout<<"TDC Calibration finish"<<endl;
#if APP
  delete theApp;
#endif
  return 0;
}
