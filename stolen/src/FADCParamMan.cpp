// FADCParamMan.cpp

#include <string>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <new>

#include "FADCParamMan.h"
#include "GlobalVariables.h"

ClassImp(FADCParamMan);

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
FADCParamMan::FADCParamMan()
{
  Clear();
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
FADCParamMan::FADCParamMan( const std::string & filename )
{
  FileName = filename;
  Clear();
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
FADCParamMan::FADCParamMan( const FADCParamMan &right )
{
  FileName=right.FileName;
  BaseFunc=right.BaseFunc;
  PeakFunc=right.PeakFunc;
  PostFunc=right.PostFunc;
  BASE_L=right.BASE_L;
  BASE_U=right.BASE_U;
  PEAK_L=right.PEAK_L;
  PEAK_U=right.PEAK_U;
  POST_L=right.POST_L;
  POST_U=right.POST_U;
  NPoint=right.NPoint;
  Interval=right.Interval;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
FADCParamMan::~FADCParamMan()
{
  Clear();
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void FADCParamMan::Clear()
{
  BaseFunc="pol0";
  PeakFunc="pol0";
  PostFunc="pol0";
  BASE_L=0;
  BASE_U=0;
  PEAK_L=0;
  PEAK_U=0;
  POST_L=0;
  POST_U=0;
  NPoint=0;
  Interval=0;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void FADCParamMan::SetFileName( const std::string & filename )
{
  FileName = filename;
}
const int MAXCHAR = 144;
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool FADCParamMan::Initialize()
{
  static const std::string funcname = "FADCParamMan::Initialize";
  std::cout << "[" << funcname << "] Initialization start ..." << std::endl;

  int n;
  //  float f;
  int nd;
  char str[MAXCHAR],str1[MAXCHAR];
  //  int ar[8];
  FILE *fp;

  if( (fp=fopen(FileName.c_str(), "r"))==0 ){
    std::cerr << " File open fail. [" << FileName << "]" << std::endl;
    exit(-1);
  }

  Clear();
 
  while( fgets(str,MAXCHAR,fp)!=0 ){
    if( str[0]=='#' ) continue;
    n=0;
    if( (nd=sscanf(str,"UpperLimitBase: %d", &n)) == 1 ) {
      BASE_U = n;
    }
    else if( (nd=sscanf(str,"LowerLimitBase: %d", &n)) == 1 ) {
      BASE_L = n;
    }
    else if( (nd=sscanf(str,"UpperLimitPeak: %d", &n)) == 1 ) {
      PEAK_U = n;
    }
    else if( (nd=sscanf(str,"LowerLimitPeak: %d", &n)) == 1 ) {
      PEAK_L = n;
    }
    else if( (nd=sscanf(str,"UpperLimitPost: %d", &n)) == 1 ) {
      POST_U = n;
    }
    else if( (nd=sscanf(str,"LowerLimitPost: %d", &n)) == 1 ) {
      POST_L = n;
    }
    else if( (nd=sscanf(str,"NPoint: %d", &n)) == 1 ) {
      NPoint = n;
    }
    else if( (nd=sscanf(str,"Interval: %d", &n)) == 1 ) {
      Interval = n;
    }
    else if( (nd=sscanf(str,"BaseFunc: %s", str1)) == 1 ) {
      BaseFunc = str1;
    }
    else if( (nd=sscanf(str,"PeakFunc: %s", str1)) == 1 ) {
      PeakFunc = str1;
    }
    else if( (nd=sscanf(str,"PostFunc: %s", str1)) == 1 ) {
      PostFunc = str1;
    }
    else{
      std::cerr << "[" << funcname << "]: Invalid data format" << std::endl;
      std::cerr << std::string(str) << std::endl;
    }
  }

  fclose(fp);

  std::cout << "[" << funcname << "] Initialization finish." << std::endl;
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
