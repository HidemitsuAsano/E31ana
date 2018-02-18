// CDSFittingParamMan.cpp

#include <string>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <new>

#include "CDSFittingParamMan.h"
#include "GlobalVariables.h"

ClassImp(CDSFittingParamMan);

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
CDSFittingParamMan::CDSFittingParamMan()
{
  MAXCDCHIT=0;
  MAGFIELD=0;
  MAXCHI=0;
  CDCTDC_U=0;
  CDCTDC_L=0;

}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
CDSFittingParamMan::CDSFittingParamMan( const std::string & filename )
{
  FileName = filename;
  MAXCDCHIT=0;
  MAGFIELD=0;
  MAXCHI=0;
  CDCTDC_U=0;
  CDCTDC_L=0;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
CDSFittingParamMan::CDSFittingParamMan( const CDSFittingParamMan &right )
{
  FileName = right.FileName;
  MAXCDCHIT=right.MAXCDCHIT;
  MAGFIELD=right.MAGFIELD;
  MAXCHI=right.MAXCHI;
  CDCTDC_U=right.CDCTDC_U;
  CDCTDC_L=right.CDCTDC_L;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
CDSFittingParamMan::~CDSFittingParamMan()
{
  Clear();
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void CDSFittingParamMan::Clear()
{
  MAXCDCHIT=0;
  MAGFIELD=0;
  MAXCHI=0;
  CDCTDC_U=0;
  CDCTDC_L=0;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void CDSFittingParamMan::SetFileName( const std::string & filename )
{
  FileName = filename;
}
const int MAXCHAR = 144;
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool CDSFittingParamMan::Initialize()
{
  static const std::string funcname = "CDSFittingParamMan::Initialize";
  std::cout << "[" << funcname << "] Initialization start ..."/* << std::endl*/;

  int n;
  float f;
  int nd;
  char str[MAXCHAR];
  FILE *fp;

  if( (fp=fopen(FileName.c_str(), "r"))==0 ){
    std::cerr << " File open fail. [" << FileName << "]" << std::endl;
    exit(-1);
  }

  Clear();
 
  while( fgets(str,MAXCHAR,fp)!=0 ){
    if( str[0]=='#' ) continue;
    n=0;
    if( (nd=sscanf(str,"MaxCDCHit: %d", &n)) == 1 ) {
      MAXCDCHIT = n;
    }
    else if( (nd=sscanf(str,"UpperLimitCDCTDC: %d", &n)) == 1 ) {
      CDCTDC_U = n;
    }
    else if( (nd=sscanf(str,"LowerLimitCDCTDC: %d", &n)) == 1 ) {
      CDCTDC_L = n;
    }
    else if( (nd=sscanf(str,"MagneticField: %f", &f)) == 1 ) {
      //      std::cout<<"mag "<<f<<std::endl;
      MAGFIELD = f;
    }
    else if( (nd=sscanf(str,"MaxChi: %d", &n)) == 1 ) {
      MAXCHI = n;
    }
    else{
      std::cerr << "[" << funcname << "]: Invalid data format" << std::endl;
      std::cerr << std::string(str) << std::endl;
    }
  }

  fclose(fp);

  std::cout << /*"[" << funcname << "] Initialization*/" finish." << std::endl;
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
