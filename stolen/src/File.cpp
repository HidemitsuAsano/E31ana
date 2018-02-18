// File.cpp

#include <iomanip>
#include <iostream>
#include <cstdio>
#include <string>
#include <new>

#include "File.h"
#include "FileType.h"
//#include "ConfMan.h"
#include "VEvent.h"
#include "GlobalVariables.h"

File::File( const std::string &infile )
  : InputFileName(infile), RunNum(0), EventNum(0)
{
  static const std::string funcname = "File::File";

  //confManager = ConfMan::GetConfManager();
  //std::cout << " Enter : " << funcname << std::endl;

  EventBuffer = new unsigned int [Event_Max_Size];
  if( EventBuffer == 0 ){
    std::cerr << "[[" << funcname << "]]: new fail" << std::endl;
    exit(-1);
  }
}

File::~File()
{
  delete [] EventBuffer;
}

bool File::Processing( ConfMan *confMan )
{
  static const std::string funcname = "File::Processing";

  int size;
  FILE *fpIn;
  bool popenedIn = false;
  
  FileCompressedType type = FileType( InputFileName.c_str() );
  if( type == NotExist ) return false;
  
  if( type == Bzip2Compressed ){
    std::string com = "bzcat "; com += InputFileName;
    if( (fpIn=popen(com.c_str(),"r"))==0 ){
      std::cerr << " file open fail. File:" << InputFileName << std::endl;
      return false;
    }
    popenedIn = true;
  }
  else if( type == GzipCompressed ){
    std::string com = "zcat "; com += InputFileName;
    if( (fpIn=popen(com.c_str(),"r"))==0 ){
      std::cerr << " file open fail. File:" << InputFileName << std::endl;
      return false;
    }
    popenedIn = true;
  }
  else{
    if( (fpIn=fopen(InputFileName.c_str(),"r"))==0 ){
      std::cerr << " file open fail. File:" << InputFileName << std::endl;
      return false;
    }
  }
      
  std::cout << " === input file : " << InputFileName << " === " << std::endl;
  int stop  = confMan->GetStopEvNum();
  int blstart  = confMan->GetStartBlockEvNum();
  int blstop  = confMan->GetStopBlockEvNum();
  // int stop  = ConfMan::GetConfManager()->GetStopEvNum();
//   int blstart  = ConfMan::GetConfManager()->GetStartBlockEvNum();
//   int blstop  = ConfMan::GetConfManager()->GetStopBlockEvNum();

  int prev_evnum=0;
  bool status=true;
  VEvent *event = new VEvent();
  while( 0<(size=ReadOneBlockEvent( EventBuffer, fpIn )) ){
    if( 0<blstart && (int)EventBuffer[3]<blstart ) continue;
#if 0
    PrintBuffer(EventBuffer);
#endif
    //VEvent *event = confManager->EventAllocator(EventBuffer,prev_evnum);
    //VEvent *event = new VEvent( EventBuffer, prev_evnum );
    event->SetBuffer( EventBuffer, prev_evnum );
    if(event){
      //event->Processing();
      status = event->Processing3(confMan);
      //event->EventCheck();
      prev_evnum = event->GetEventNumber();
#if 0
      std::cout << " PrevEvNum:" << prev_evnum << std::endl;
#endif
      //delete event;
      if( !status ) break;
      if( 0<stop && stop<prev_evnum ) break;
      if( 0<blstop && blstop<(int)EventBuffer[3] ){ std::cout << "stop" << std::endl; break;}
    }
    else{
      std::cerr << " Error in Event Construction : "
		<< funcname << std::endl;
    }
  
  }
  delete event;
      
  if( popenedIn ) pclose(fpIn);
  else            fclose(fpIn);

  return true;
}

int File::ReadOneBlockEvent( unsigned int *buf, FILE *fp )
{
  static const std::string funcname = "File::ReadOneBlockEvent";
  int size,ret;
#if 0
  std::cout << " Enter: " << funcname << std::endl;
#endif

  if( (ret=fread(buf,sizeof(unsigned int),1,fp))!=1 )
    return 0;
  size = buf[0]-1;
  if( (ret=fread(&buf[1],sizeof(unsigned int),size,fp))!=size ){
    std::cerr << " Read Error : " << funcname
	      << " (val,size)=(" << ret << ", " << size << ")" << std::endl;
    return 0;
  }
  return ret+1;
}

int File::ReadOneEvent( unsigned int *buf, FILE *fp )
{
  return 0;
}

void File::PrintBuffer( unsigned int *buf )
{
  static const std::string funcname = "File::PrintBuffer";
  int size;
  std::cout << " Enter: " << funcname << std::endl;

  size = buf[0];

  std::cout << "DAQ HEAD:";
  std::cout << std::hex;
  for( int i=0; i<6; i++ ){
    std::cout << "  " << std::setw(10) << buf[i];
  }
  std::cout << "\n";
  std::cout << std::hex;
  for( int i=6; i<size; i++ ){
    std::cout << "  " << std::setw(10) << buf[i];
    if( (i-5)%4==0 ) std::cout << "\n";
  }
  std::cout << std::dec << "\n";

}
