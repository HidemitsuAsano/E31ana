// CDSFittingParamMan.cpp

#include <string>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <new>

#include "CDSFittingParamMan.h"
#include "GlobalVariables.h"
#include "TrackTools.h"

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
CDSFittingParamMan::CDSFittingParamMan( const std::string & filename ):FileName(filename)
{
  //FileName = filename; asano
  MAXCDCHIT=0;
  MAGFIELD=0;
  MAXCHI=0;
  CDCTDC_U=0;
  CDCTDC_L=0;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
CDSFittingParamMan::CDSFittingParamMan( const CDSFittingParamMan &right ):FileName(right.FileName)
{
  //FileName = right.FileName; asano
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
  int type;
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
    else if( (nd=sscanf(str,"PID: %d %d", &n, &type)) == 2 ) {
      if( type==1 ){
	double p0, p1;
	fgets(str,MAXCHAR,fp);
	if( (nd=sscanf(str, "%lf %lf", &p0, &p1)) == 2 ){
	  PIDFuncMin[n]= new TF1(Form("pid%d_min",n), "pol0(0)");
	  PIDFuncMin[n]->SetParameter(0, p0);
	  PIDFuncMax[n]= new TF1(Form("pid%d_max",n), "pol0(0)");
	  PIDFuncMax[n]->SetParameter(0, p1);
	}
	else{
	  std::cerr << "[" << funcname << "]: Invalid data format for PID 1D-function" << std::endl;
	  std::cerr << std::string(str) << std::endl;
	}
      }
      else if( type==2 ){
	double p0, p1, p2, p3, p4;
	fgets(str,MAXCHAR,fp);
	if( (nd=sscanf(str, "%lf %lf %lf %lf %lf", &p0, &p1, &p2, &p3, &p4)) == 5 ){
	  //	  std::cout << p0<<" "<<p1<<" "<<p2<<" "<<p3<<" "<<p4<<std::endl;
	  PIDFuncMin[n]= new TF1(Form("pid%d_min",n), "[0]-[1]*sqrt(4*[2]*[5]*[5]*x*x+4*[5]*[5]*[5]*[5]*[3]*(1+([5]*[5])/(x*x))+4*[4]*x*x*([5]*[5]+x*x))");
	  PIDFuncMin[n]->SetParameter(0, p0);
	  PIDFuncMin[n]->SetParameter(1, p1);
	  PIDFuncMin[n]->SetParameter(2, p2);
	  PIDFuncMin[n]->SetParameter(3, p3);
	  PIDFuncMin[n]->SetParameter(4, p4);
	  PIDFuncMin[n]->SetParameter(5, cdsMass[n]);


	  PIDFuncMax[n]= new TF1(Form("pid%d_max",n), "[0]+[1]*sqrt(4*[2]*[5]*[5]*x*x+4*[5]*[5]*[5]*[5]*[3]*(1+([5]*[5])/(x*x))+4*[4]*x*x*([5]*[5]+x*x))");
	  PIDFuncMax[n]->SetParameter(0, p0);
	  PIDFuncMax[n]->SetParameter(1, p1);
	  PIDFuncMax[n]->SetParameter(2, p2);
	  PIDFuncMax[n]->SetParameter(3, p3);
	  PIDFuncMax[n]->SetParameter(4, p4);
	  PIDFuncMin[n]->SetParameter(5, cdsMass[n]);
	}
	else{
	  std::cerr << "[" << funcname << "]: Invalid data format for PID 2D-function" << std::endl;
	  std::cerr << std::string(str) <<nd<< std::endl;
	}
      }
      else{
	std::cerr << "[" << funcname << "]: PID Func not support type=" <<type<< std::endl;
      }
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
bool CDSFittingParamMan::PID(int pid, double mom, double mass2)
{
  if( (pid==CDS_PiPlus || pid==CDS_Proton || pid==CDS_Deuteron || pid==CDS_Triton || pid==CDS_Helium3) && mom<0 ) return false;
  if( (pid==CDS_PiMinus || pid==CDS_Kaon ) && mom>0 ) return false;

  if( PIDFuncMin.find(pid)!=PIDFuncMin.end() ){
    //    std::cout<<" CDS PID use function"<<std::endl;
    double min=PIDFuncMin.at(pid)->Eval(mom);
    double max=PIDFuncMax.at(pid)->Eval(mom);
    //    if( pid==CDS_PiMinus ) std::cout<<"mom  min  max  "<<mom<<"  "<<min<<"  "<<max<<std::endl;
    //    std::cout<<"mass2  mom  "<<mass2<<"  "<<mom<<std::endl;
    if( min<mass2 && mass2<max ) return true;
  }
  else{
    if( TrackTools::PID(mom, mass2)==pid ) return true;
  }

  return false;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
int CDSFittingParamMan::PID(double mom, double mass2)
{
  bool flag[8];
  for( int i=0; i<7; i++ ) flag[i]=PID(i, mom, mass2);

  int nflag=0;
  for( int i=0; i<7; i++ ) if( flag[i] ) nflag++;
  if( nflag==1 ){
    for( int i=0; i<7; i++ ) if( flag[i] ) return i;
  }

  return CDS_Other;
}
