// AnalysisSim.cpp
// 2014.12.22
// T. Yamaga


#include <iostream>
#include <string>
#include <new>
#include <cstdio>

#include "ConfMan.h"
#include "MyAnalysisBase.h"
#include "MyAnalysisAlloc.h"
#include "CDSHitMan.h"
#include "BeamLineHitMan.h"
#include "EventHeader.h"
#include "ScalerMan.h"
#include "BeamLineTrackMan.h"
#include "CDSTrackingMan.h"

const int MaxChar = 144;

int main( int argc, char **argv )
{
  std::cout << "Start!" << std::endl;

  std::string confFile = "conf/analyzer.conf";
  std::string inFile;
  std::string rootFile;
  std::string cdstrackFile;
  std::string mtdcFile;
  std::string CommandName = argv[0];
  int         RunNum = -9999;
  

  // Arguments check  
  if(argc!=4){
    std::cout << "Input ConfFile, OutFile, and InFile " << std::endl;
    return 0;
  }
  for(int i=0;i<argc;i++)
    std::cout<<i<<"  "<<argv[i]<<std::endl;
  if( argv[1][0] != '!' ) confFile  = argv[1];
  if( argv[2][0] != '!' ) rootFile  = argv[2];
  if( argv[3][0] != '!' ) inFile    = argv[3];

  ConfMan *confManager=0;
  if(argc==4)    confManager = new ConfMan( confFile, RunNum, rootFile );
  if(!confManager){
    std::cout<<" ConfManager cannot be opened properly." <<std::endl;
    return 0;
  }
  confManager->SetProgramName(argv[0]);
  confManager->Initialize();

  std::cout << "InFile="<<inFile << std::endl;

  TFile* rtFile = TFile::Open(inFile.c_str());
  if(!rtFile->IsOpen()) {
    std::cout << "Input file open error !!!" << std::endl;
    return 0;
  }
  rtFile -> cd();
  TTree* tree = (TTree*)gFile->Get("EventTree");
  CDSHitMan* cdsMan = 0;
  BeamLineHitMan* blMan = 0;
  CDSTrackingMan* cdsTrackMan = 0;
  EventHeader* header = 0;
  
  tree -> SetBranchAddress("CDSHitMan",&cdsMan);
  tree -> SetBranchAddress("BeamLineHitMan",&blMan);
  tree -> SetBranchAddress("CDSTrackingMan",&cdsTrackMan);

  MyAnalysisAlloc* myalloc = new MyAnalysisAlloc();
  MyAnalysisBase* myana = myalloc->MyAnalysisAllocator();
  myana -> Initialize(confManager);
  for(int i=0; i<tree->GetEntries(); i++){
    if(i%5000==0)
      std::cout << "Event# : " << i << " / Total# : " << tree->GetEntries() << std::endl;

    tree -> GetEvent(i);
    myana -> SetTrackMan(cdsTrackMan);
    myana -> DoAnalysis(cdsMan,blMan);
  }

  myana -> Finalize(); 

  delete myalloc;
  delete confManager;

  std::cout << "Finish!" << std::endl;
  return 0;
}

