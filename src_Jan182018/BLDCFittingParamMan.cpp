// BLDCFittingParamMan.cpp

#include <string>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <new>

#include "BLDCFittingParamMan.h"
#include "GlobalVariables.h"

ClassImp(BLDCFittingParamMan);

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
BLDCFittingParamMan::BLDCFittingParamMan()
{
  Clear();
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
BLDCFittingParamMan::BLDCFittingParamMan( const std::string & filename )
{
  FileName = filename;
  Clear();
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
BLDCFittingParamMan::BLDCFittingParamMan( const BLDCFittingParamMan &right )
{
  FileName=right.FileName;
  MAXBLDCHIT=right.MAXBLDCHIT;
  MAXHITinLAYER=right.MAXHITinLAYER;
  MAXHITinTRACK=right.MAXHITinTRACK;
  MAGFIELD=right.MAGFIELD;
  MAXCHI=right.MAXCHI;
  MAXCHIPRE=right.MAXCHIPRE;
  MAXSLOPE=right.MAXSLOPE;

  MinHitXBLC1 =right.MinHitXBLC1;
  MinHitXBLC1a=right.MinHitXBLC1a;
  MinHitXBLC1b=right.MinHitXBLC1b;
  MinHitXBLC2 =right.MinHitXBLC2;
  MinHitXBLC2a=right.MinHitXBLC2a;
  MinHitXBLC2b=right.MinHitXBLC2b;
  MinHitXBPC =right.MinHitXBPC;


  MinHitYBLC1 =right.MinHitYBLC1;
  MinHitYBLC1a=right.MinHitYBLC1a;
  MinHitYBLC1b=right.MinHitYBLC1b;
  MinHitYBLC2 =right.MinHitYBLC2;
  MinHitYBLC2a=right.MinHitYBLC2a;
  MinHitYBLC2b=right.MinHitYBLC2b;
  MinHitYBPC =right.MinHitYBPC;
  BLDCTDC_U=right.BLDCTDC_U;
  BLDCTDC_L=right.BLDCTDC_L;
  for(int i=0;i<8;i++){
    LayerDeadBLC1a[i]=right.LayerDeadBLC1a[i];
    LayerDeadBLC1b[i]=right.LayerDeadBLC1b[i];
    LayerDeadBLC2a[i]=right.LayerDeadBLC2a[i];
    LayerDeadBLC2b[i]=right.LayerDeadBLC2b[i];
    LayerDeadBPC[i]=right.LayerDeadBPC[i];
    LayerKilledBLC1a[i]=right.LayerKilledBLC1a[i];
    LayerKilledBLC1b[i]=right.LayerKilledBLC1b[i];
    LayerKilledBLC2a[i]=right.LayerKilledBLC2a[i];
    LayerKilledBLC2b[i]=right.LayerKilledBLC2b[i];
    LayerKilledBPC[i]=right.LayerKilledBPC[i];
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
BLDCFittingParamMan::~BLDCFittingParamMan()
{
  Clear();
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void BLDCFittingParamMan::Clear()
{
  MAXBLDCHIT=0;
  MAXHITinLAYER=0;
  MAXHITinTRACK=0;
  MAGFIELD=0;
  MAXCHI=0;
  MAXCHIPRE=0;
  MAXSLOPE=0;

  MinHitXBLC1=0;
  MinHitXBLC1a=0;
  MinHitXBLC1b=0;
  MinHitXBLC2=0;
  MinHitXBLC2a=0;
  MinHitXBLC2b=0;
  MinHitXBPC=0;

  MinHitYBLC1=0;
  MinHitYBLC1a=0;
  MinHitYBLC1b=0;
  MinHitYBLC2=0;
  MinHitYBLC2a=0;
  MinHitYBLC2b=0;
  MinHitYBPC=0;

  BLDCTDC_U=0;
  BLDCTDC_L=0;
  for(int i=0;i<8;i++){
    LayerDeadBLC1a[i]=false;
    LayerDeadBLC1b[i]=false;
    LayerDeadBLC2a[i]=false;
    LayerDeadBLC2b[i]=false;
    LayerDeadBPC[i]=false;
    LayerDeadFDC1[i]=false;
    LayerKilledBLC1a[i]=false;
    LayerKilledBLC1b[i]=false;
    LayerKilledBLC2a[i]=false;
    LayerKilledBLC2b[i]=false;
    LayerKilledBPC[i]=false;
    LayerKilledFDC1[i]=false;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void BLDCFittingParamMan::SetFileName( const std::string & filename )
{
  FileName = filename;
}
const int MAXCHAR = 144;
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool BLDCFittingParamMan::Initialize()
{
  static const std::string funcname = "BLDCFittingParamMan::Initialize";
  std::cout << "[" << funcname << "] Initialization start ..."/* << std::endl*/;

  int n;
  float f;
  int nd;
  char str[MAXCHAR];
  int ar[8];
  FILE *fp;

  if( (fp=fopen(FileName.c_str(), "r"))==0 ){
    std::cerr << " File open fail. [" << FileName << "]" << std::endl;
    exit(-1);
  }

  Clear();
 
  while( fgets(str,MAXCHAR,fp)!=0 ){
    if( str[0]=='#' ) continue;
    n=0;
    if( (nd=sscanf(str,"MaxBLDCHit: %d", &n)) == 1 ) {
      MAXBLDCHIT = n;
    }
    else if( (nd=sscanf(str,"MaxHitInLayer: %d", &n)) == 1 ) {
      MAXHITinLAYER = n;
    }
    else if( (nd=sscanf(str,"MaxHitInTrack: %d", &n)) == 1 ) {
      MAXHITinTRACK = n;
    }
    else if( (nd=sscanf(str,"UpperLimitBLDCTDC: %d", &n)) == 1 ) {
      BLDCTDC_U = n;
    }
    else if( (nd=sscanf(str,"LowerLimitBLDCTDC: %d", &n)) == 1 ) {
      BLDCTDC_L = n;
    }
    else if( (nd=sscanf(str,"MagneticField: %f", &f)) == 1 ) {
      std::cout<<"mag "<<f<<std::endl;
      MAGFIELD = f;
    }
    else if( (nd=sscanf(str,"MaxChi: %d", &n)) == 1 ) {
      MAXCHI = n;
    }
    else if( (nd=sscanf(str,"MaxChiPreTracking: %d", &n)) == 1 ) {
      MAXCHIPRE = n;
    }
    else if( (nd=sscanf(str,"MaxSlope: %f", &f)) == 1 ) {
      MAXSLOPE = f;
    }
    else if( (nd=sscanf(str,"MinHitXBLC1 : %d", &n)) == 1 ) {
      MinHitXBLC1 = n;
    }
    else if( (nd=sscanf(str,"MinHitXBLC1a: %d", &n)) == 1 ) {
      MinHitXBLC1a = n;
    }
    else if( (nd=sscanf(str,"MinHitXBLC1b: %d", &n)) == 1 ) {
      MinHitXBLC1b = n;
    }
    else if( (nd=sscanf(str,"MinHitXBLC2: %d", &n)) == 1 ) {
      MinHitXBLC2 = n;
    }
    else if( (nd=sscanf(str,"MinHitXBLC2a: %d", &n)) == 1 ) {
      MinHitXBLC2a = n;
    }
    else if( (nd=sscanf(str,"MinHitXBLC2b: %d", &n)) == 1 ) {
      MinHitXBLC2b = n;
    }
    else if( (nd=sscanf(str,"MinHitXBPC: %d", &n)) == 1 ) {
      MinHitXBPC = n;
    }
    else if( (nd=sscanf(str,"MinHitXFDC1: %d", &n)) == 1 ) {
      MinHitXFDC1 = n;
    }
    else if( (nd=sscanf(str,"MinHitYBLC1 : %d", &n)) == 1 ) {
      MinHitYBLC1 = n;
    }
    else if( (nd=sscanf(str,"MinHitYBLC1a: %d", &n)) == 1 ) {
      MinHitYBLC1a = n;
    }
    else if( (nd=sscanf(str,"MinHitYBLC1b: %d", &n)) == 1 ) {
      MinHitYBLC1b = n;
    }
    else if( (nd=sscanf(str,"MinHitYBLC2: %d", &n)) == 1 ) {
      MinHitYBLC2 = n;
    }
    else if( (nd=sscanf(str,"MinHitYBLC2a: %d", &n)) == 1 ) {
      MinHitYBLC2a = n;
    }
    else if( (nd=sscanf(str,"MinHitYBLC2b: %d", &n)) == 1 ) {
      MinHitYBLC2b = n;
    }
    else if( (nd=sscanf(str,"MinHitYBPC : %d", &n)) == 1 ) {
      MinHitYBPC = n;
    }
    else if( (nd=sscanf(str,"MinHitYFDC1 : %d", &n)) == 1 ) {
      MinHitYFDC1 = n;
    }
    else if( (nd=sscanf(str,"LayerBLC1a: %d %d %d %d %d %d %d %d ",
			&ar[0],&ar[1],&ar[2],&ar[3],
			&ar[4],&ar[5],&ar[6],&ar[7])) == 8 ) {
      for(int i=0;i<8;i++){
	if(ar[i]==1) continue;
	else if(ar[i]==0) LayerKilledBLC1a[i]=true; 
	else if(ar[i]==-1){
	  LayerKilledBLC1a[i]=true; 
	  LayerDeadBLC1a[i]=true; 
	}
      }
    }
    else if( (nd=sscanf(str,"LayerBLC1b: %d %d %d %d %d %d %d %d ",
			&ar[0],&ar[1],&ar[2],&ar[3],
			&ar[4],&ar[5],&ar[6],&ar[7])) == 8 ) {
      for(int i=0;i<8;i++){
	if(ar[i]==1) continue;
	else if(ar[i]==0) LayerKilledBLC1b[i]=true; 
	else if(ar[i]==-1){
	  LayerKilledBLC1b[i]=true; 
	  LayerDeadBLC1b[i]=true; 
	}
      }
    }
    else if( (nd=sscanf(str,"LayerBLC2a: %d %d %d %d %d %d %d %d ",
			&ar[0],&ar[1],&ar[2],&ar[3],
			&ar[4],&ar[5],&ar[6],&ar[7])) == 8 ) {
      for(int i=0;i<8;i++){
	if(ar[i]==1) continue;
	else if(ar[i]==0) LayerKilledBLC2a[i]=true; 
	else if(ar[i]==-1){
	  LayerKilledBLC2a[i]=true; 
	  LayerDeadBLC2a[i]=true; 
	}
      }
    }
    else if( (nd=sscanf(str,"LayerBLC2b: %d %d %d %d %d %d %d %d ",
			&ar[0],&ar[1],&ar[2],&ar[3],
			&ar[4],&ar[5],&ar[6],&ar[7])) == 8 ) {
      for(int i=0;i<8;i++){
	if(ar[i]==1) continue;
	else if(ar[i]==0) LayerKilledBLC2b[i]=true; 
	else if(ar[i]==-1){
	  LayerKilledBLC2b[i]=true; 
	  LayerDeadBLC2b[i]=true; 
	}
      }
    }
    else if( (nd=sscanf(str,"LayerBPC: %d %d %d %d %d %d %d %d ",
			&ar[0],&ar[1],&ar[2],&ar[3],
			&ar[4],&ar[5],&ar[6],&ar[7])) == 8 ) {
      for(int i=0;i<8;i++){
	if(ar[i]==1) continue;
	else if(ar[i]==0) LayerKilledBPC[i]=true; 
	else if(ar[i]==-1){
	  LayerKilledBPC[i]=true; 
	  LayerDeadBPC[i]=true; 
	}
      }
    }
    else if( (nd=sscanf(str,"LayerFDC1: %d %d %d %d %d %d %d %d ",
			&ar[0],&ar[1],&ar[2],&ar[3],
			&ar[4],&ar[5],&ar[6],&ar[7])) == 8 ) {
      for(int i=0;i<8;i++){
	if(ar[i]==1) continue;
	else if(ar[i]==0) LayerKilledFDC1[i]=true; 
	else if(ar[i]==-1){
	  LayerKilledFDC1[i]=true; 
	  LayerDeadFDC1[i]=true; 
	}
      }
    }
    else{
      std::cerr << "[" << funcname << "]: Invalid data format" << std::endl;
      std::cerr << std::string(str) << std::endl;
    }
  }

  fclose(fp);

  std::cout <</* "[" << funcname << "] Initialization*/" finish." << std::endl;
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool BLDCFittingParamMan::layerdead( const int &id, const int &i )
{
  if( id==CID_BLC1a )      return layerdeadBLC1a(i);
  else if( id==CID_BLC1b ) return layerdeadBLC1b(i);
  else if( id==CID_BLC2a ) return layerdeadBLC2a(i);
  else if( id==CID_BLC2b ) return layerdeadBLC2b(i);
  else if( id==CID_BPC )   return layerdeadBPC(i);
  else if( id==CID_FDC1 )   return layerdeadFDC1(i);
  else{
    std::cout<<"BLDCFittingParamMan: no parameter for id: "<<id<<" and layer: "<<i<<std::endl;
    return false;
  }
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool BLDCFittingParamMan::layerkilled( const int &id, const int &i )
{
  if(      id==CID_BLC1a ) return layerkilledBLC1a(i);
  else if( id==CID_BLC1b ) return layerkilledBLC1b(i);
  else if( id==CID_BLC2a ) return layerkilledBLC2a(i);
  else if( id==CID_BLC2b ) return layerkilledBLC2b(i);
  else if( id==CID_BPC )   return layerkilledBPC(i);
  else if( id==CID_FDC1 )  return layerkilledFDC1(i);
  else{
    std::cout<<"BLDCFittingParamMan: no parameter for id: "<<id<<" and layer: "<<i<<std::endl;
    return false;
  }
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
int BLDCFittingParamMan::GetMinHit( const int &xy, const int &id ){
  if(xy==0)  return GetMinHitX(id);
  else if(xy==1)  return GetMinHitY(id);
  return -1;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
int BLDCFittingParamMan::GetMinHitX( const int &id )
{
  if(      id==CID_BLC1a ) return GetMinHitXBLC1a();
  else if( id==CID_BLC1b ) return GetMinHitXBLC1b();
  else if( id==CID_BLC2a ) return GetMinHitXBLC2a();
  else if( id==CID_BLC2b ) return GetMinHitXBLC2b();
  else if( id==CID_BPC ) return GetMinHitXBPC();
  else if( id==CID_FDC1 ) return GetMinHitXFDC1();
  else{
    std::cout<<"BLDCFittingParamMan: no parameter for id: "<<id<<std::endl;
    return false;
  }
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
int BLDCFittingParamMan::GetMinHitY( const int &id )
{
  if(      id==CID_BLC1a ) return GetMinHitYBLC1a();
  else if( id==CID_BLC1b ) return GetMinHitYBLC1b();
  else if( id==CID_BLC2a ) return GetMinHitYBLC2a();
  else if( id==CID_BLC2b ) return GetMinHitYBLC2b();
  else if( id==CID_BPC ) return GetMinHitYBPC();
  else if( id==CID_FDC1 ) return GetMinHitYFDC1();
  else{
    std::cout<<"BLDCFittingParamMan: no parameter for id: "<<id<<std::endl;
    return false;
  }
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
