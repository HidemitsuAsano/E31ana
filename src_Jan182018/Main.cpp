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
  std::string vmeFile;
  std::string CommandName = argv[0];

  // Arguments check  
  if(argc!=4&&argc!=5){
    std::cout << "Input ConfFile, OutFile, and InFile , (and VMErootFile)" << std::endl;
    return 0;
  }
  
  if( argv[1][0] != '!' ) confFile = argv[1];
  if( argv[2][0] != '!' ) rootFile = argv[2];
  if( argv[3][0] != '!' ) inFile = argv[3];
  if( argc==5 && argv[4][0] != '!' ) vmeFile = argv[4];

  ConfMan *confManager;
  if(argc==4)    confManager = new ConfMan( confFile, rootFile );
  if(argc==5)    confManager = new ConfMan( confFile, rootFile, vmeFile );
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
