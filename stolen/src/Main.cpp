// Main.cpp

#include <iostream>
#include <string>
#include <new>
#include <cstdio>

#include "File.h"
#include "BLEvent.h"
#include "ConfMan.h"

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
  if(argc!=5&&argc!=6&&argc!=7){
    std::cout << "Input ConfFile, runnumber, OutFile, and InFile , (and CDSTrackRootFile, MTDCRootFile)" << std::endl;
    return 0;
  }
  for(int i=0;i<argc;i++)
    std::cout<<i<<"  "<<argv[i]<<std::endl;
  if( argv[1][0] != '!' ) confFile  = argv[1];
  if( argv[2][0] != '!' ) RunNum    = atoi(argv[2]);
  if( argv[3][0] != '!' ) rootFile  = argv[3];
  if( argv[4][0] != '!' ) inFile    = argv[4];
  if( argc>=6 && argv[5][0] != '!' ) cdstrackFile = argv[5];
  if( argc==7 && argv[6][0] != '!' ) mtdcFile     = argv[6];

  ConfMan *confManager=0;
  if(argc==5)    confManager = new ConfMan( confFile, RunNum, rootFile );
  if(argc==6)    confManager = new ConfMan( confFile, RunNum, rootFile, cdstrackFile );
  if(argc==7)    confManager = new ConfMan( confFile, RunNum, rootFile, cdstrackFile, mtdcFile );
  if(!confManager){
    std::cout<<" ConfManager cannot be opened properly." <<std::endl;
    return 0;
  }
  confManager->SetProgramName(argv[0]);
  confManager->Initialize();

  std::cout << "InFile="<<inFile << std::endl;

  //
  if( confManager->GetDataType()==Type_CDS1 ){
    File *file = new File( inFile );
    file->Processing(confManager);
    delete file;
  }
  else if( confManager->GetDataType()==Type_BL1 ){
    BLEvent *blev = new BLEvent();
    blev->BMain( argc, argv, confManager );
    delete blev;
  }
  else if( confManager->GetDataType()==Type_SDD1 ){
    BLEvent *blev = new BLEvent();
    argc=4;
    blev->BMain( argc , argv, confManager );
    delete blev;
  }

  delete confManager;

  std::cout << "Finish!" << std::endl;
  return 0;
}
