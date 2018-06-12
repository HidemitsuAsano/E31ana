#include <iostream>
#include <string>

#include <TChain.h>
#include <TFile.h>
#include <TString.h>

#include "AnaData.h"
#include "EventHeader.h"
#include "ScalerMan.h"

int main(int argc, char** argv)
{

  if( argc!=2 ){
    std::cout<<"./readAnaData  ($ReadFile)"<<std::endl;
    return 0;
  }

  std::string str;
  std::ifstream readFile(argv[1]);

  std::vector<int> readrun;
  std::string confname;
  int runnum = -1;
  TString outfile;

  while( !readFile.eof() ){
    char c_str[1024];
    int tmprun, tmpstart, tmpend;
    getline(readFile, str);
    if( str[0]=='#' ) continue;

    if( sscanf(str.c_str(), "ConfFile: %s %d", c_str, &tmprun)==2 ){
      confname = c_str;
      runnum = tmprun;
    }
    else if( sscanf(str.c_str(), "ConfFile: %s", c_str)==1 ){
      confname = c_str;
    }

    if( sscanf(str.c_str(), "run: %d %d", &tmpstart, &tmpend)==2 ){
      if( tmpstart>tmpend){
	std::cout<<"  !!! start>end !!!"<<std::endl;
	return 0;
      }
      for( int i=tmpstart; i<=tmpend; i++ ) readrun.push_back(i);
    }
    else if( sscanf(str.c_str(), "run: %d", &tmpstart)==1 ){
      readrun.push_back(tmpstart);
    }

    if( sscanf(str.c_str(), "OutFile: %s", c_str)==1 ){
      outfile = c_str;
    }
  }

  TFile *of = new TFile(outfile, "recreate");

  TChain *scaTree = new TChain("ScalerTree");
  ScalerMan *scaMan = new ScalerMan();
  scaTree-> SetBranchAddress("ScalerMan", &scaMan);

  TChain *evTree = new TChain("EventTree");
  EventHeader *header = new EventHeader();
  AnaData *anaData = new AnaData();
  evTree-> SetBranchAddress("EventHeader", &header);
  evTree-> SetBranchAddress("AnaData", &anaData);

  for( int i=0; i<readrun.size(); i++ ){
    scaTree-> Add(Form("Run62/root/anaData_%d.root", readrun[i]));
    evTree-> Add(Form("Run62/root/anaData_%d.root", readrun[i]));
  }

  std::cout<<"All Event Number : "<<evTree->GetEntries()<<std::endl;
  std::cout<<"      |=========50%=========|"<<std::endl;
  std::cout<<"Start |"<<std::flush;
  for( int ev=0; ev<evTree->GetEntries(); ev++ ){
    if( ev%(evTree->GetEntries()/20)==0 ) std::cout<<"="<<std::flush;
    evTree-> GetEntries();
  }
  std::cout<<"|Finish"<<std::endl;

  of-> Write();
  of-> Close();

  return 0;
}
