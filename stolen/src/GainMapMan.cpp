// GainMapMan.cpp

#include <string>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <new>

#include "GainMapMan.h"
#include "GlobalVariables.h"
#include "TRandom.h"

ClassImp(GainMap);
ClassImp(GainMapMan);

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
GainMap::GainMap()
{
  Param[0]=0.;
  Param[1]=1.;
  Param[2]=0.;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
GainMapMan::GainMapMan()
{
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
GainMapMan::GainMapMan( const GainMapMan &right )
{
  FileNameCDS = right.FileNameCDS;
  FileNameBL  = right.FileNameBL;
  for( GainMapContainer::const_iterator i=right.gainContainer.begin();
       i!=right.gainContainer.end(); i++ ){
    gainContainer[i->first] = i->second;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
GainMapMan::~GainMapMan()
{
  gainContainer.clear();
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void GainMapMan::SetFileNameCDS( const std::string & filename )
{
  FileNameCDS = filename;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void GainMapMan::SetFileNameBL( const std::string & filename )
{
  FileNameBL = filename;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void GainMapMan::SetFileNameSDD( const std::string & filename )
{
  FileNameSDD = filename;
}

const unsigned int KEYMASK  = 0x000F;
// |0000|0000|1111|0001|1111|0011|1111|0011|
const unsigned int AMASK    = 0x003F;      /* A Mask 6 Bits (0-63) */
const unsigned int NMASK    = 0x001F;      /* N Mask 5 Bits (0-31) */
const unsigned int CMASK    = 0x000F;      /* C Mask 4 Bits (0-15) */
const int          ASHIFT   =  4;
const int          NSHIFT   = 12;
const int          CSHIFT   = 20;
const unsigned int KEYFLAG  = 0x0003;

#define KEY(c,n,a) \
((((c)&CMASK)<<CSHIFT) | (((n)&NMASK)<<NSHIFT) | (((a)&AMASK)<<ASHIFT) | KEYFLAG )

const int MAXCHAR = 144;

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool GainMapMan::Initialize()
{
  static const std::string funcname = "GainMapMan::Initialize";
  std::cout << "[" << funcname << "] Initialization start ... " ;

  gainContainer.clear();

  if(ReadFile(FileNameCDS))  std::cout << " CDS." ;//<< std::endl;
  if(ReadFile(FileNameBL))   std::cout << " BL." ;//<< std::endl;
  if(ReadFile(FileNameSDD))   std::cout << " SDD." ;//<< std::endl;
  std::cout <</* "[" << funcname << "] Initialization*/" finish." << std::endl;
  return true;
}

bool GainMapMan::ReadFile( const std::string &filename)
{
  if( filename==DefaultFileName ) return false;
  int c,n,a;
  double gain, offset;
  double par[3];
  int nd;
  unsigned int key;
  char str[MAXCHAR];
  FILE *fp;

  if( (fp=fopen(filename.c_str(), "r"))==0 ){
    std::cerr << " File open fail. [" << filename << "]" << std::endl;
    exit(-1);
  }
  else{
    while( fgets(str,MAXCHAR,fp)!=0 ){
      if( str[0]=='#' ) continue;	
      if( (nd=sscanf(str,"%d %d %d %lf %lf %lf", &c, &n, &a, &par[0], &par[1], &par[2])) == 6 ) {
	key = KEY(c,n,a);
	GainMap *again = new GainMap();
	again->SetParams(par);
	gainContainer[key] = *again;
	delete again;
      }
      else if( (nd=sscanf(str,"%d %d %d %lf %lf", &c, &n, &a, &gain, &offset)) == 5 ) {
#if 0
	std::cout << c << "  " << n << "  " << a << "  "
		  << gain << "  " << offset << std::endl;
#endif
	key = KEY(c,n,a);
	GainMap *again = new GainMap();
	again->SetGain(gain);
	again->SetOffset(offset);
	gainContainer[key] = *again;
	delete again;
      }
      else{
	std::cerr << "[" << filename << "]: Invalid data format" << std::endl;
	std::cerr << std::string(str) << std::endl;
      }
    } 
    fclose(fp);
  }
  return true;
}


// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool GainMapMan::GetParam( const int &c, const int &n, const int &a, double &gain, double &offset )
{
  static const std::string funcname = "GainMapMan::GetParam";
  unsigned int key;
  GainMap again;

  key = KEY(c,n,a);
  GainMapContainer::iterator ig = gainContainer.find(key);
  if( ig != gainContainer.end() ){
    gain = (ig->second).GetGain();
    offset = (ig->second).GetOffset();
    return true;
  }
  else{
    std::cout << " Invalid value!!! [" << funcname << "]"
	      << " c:" << c << " n:" << n << " a:" << a
	      << std::endl;
    return false;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool GainMapMan::GetParam( const int &c, const int &n, const int &a, double &p0, double &p1, double &p2 )
{
  static const std::string funcname = "GainMapMan::GetParam";
  unsigned int key;
  GainMap again;

  key = KEY(c,n,a);
  GainMapContainer::iterator ig = gainContainer.find(key);
  if( ig != gainContainer.end() ){
    p0 = (ig->second).GetParam(0);
    p1 = (ig->second).GetParam(1);
    p2 = (ig->second).GetParam(2);
    return true;
  }
  else{
    std::cout << " Invalid value!!! [" << funcname << "]"
	      << " c:" << c << " n:" << n << " a:" << a
	      << std::endl;
    return false;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool GainMapMan::SetParam( const int &c, const int &n, const int &a, const double &gain, const double &offset )
{
  static const std::string funcname = "GainMapMan::SetParam";
  unsigned int key;

  key = KEY(c,n,a);
  GainMapContainer::iterator ig = gainContainer.find(key);
  if( ig != gainContainer.end() ){
    (ig->second).SetGain( gain );
    (ig->second).SetOffset( offset );
    return true;
  }
  else{
    std::cout << " Invalid value!!! [" << funcname << "]"
	      << " c:" << c << " n:" << n << " a:" << a
	      << std::endl;
    return false;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool GainMapMan::SetParam( const int &c, const int &n, const int &a, const double &p0, const double &p1, const double &p2 )
{
  static const std::string funcname = "GainMapMan::SetParam";
  unsigned int key;

  key = KEY(c,n,a);
  GainMapContainer::iterator ig = gainContainer.find(key);
  if( ig != gainContainer.end() ){
    (ig->second).SetParam( 0, p0 );
    (ig->second).SetParam( 1, p1 );
    (ig->second).SetParam( 2, p2 );
    return true;
  }
  else{
    std::cout << " Invalid value!!! [" << funcname << "]"
	      << " c:" << c << " n:" << n << " a:" << a
	      << std::endl;
    return false;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
double GainMapMan::CalcCValue( const int &c, const int &n, const int &a,
			       const int &at, const double &val, const bool &smear )
{
  static const std::string funcname = "GainMapMan::CalcCValue";
  unsigned int key;
  double cval,val2;

  key = KEY(c,n,a);

  GainMapContainer::iterator ig = gainContainer.find(key);
  if( ig != gainContainer.end() ){
    /* dE = gain(ADC-pedestal) */
    /* Time = gain*TDC-offset */
    if(smear){
      val2=val+gRandom->Rndm()-0.5;
    }
    else {
      val2=val;
    }
    
//    if( at==0 ) cval = ((ig->second).GetGain())*(val2 - ((ig->second).GetOffset()));
//    else        cval = ((ig->second).GetGain())* val2 - ((ig->second).GetOffset());
    if( at==0 ){ cval = ((ig->second).GetGain())*(val2 - ((ig->second).GetOffset())); }
    else       { cval = ((ig->second).GetGain())* val2 - ((ig->second).GetOffset());  }
    //else       { cval = ((ig->second).GetGain())*val2 +(ig->second).GetParam(2)*val2*val2 - ((ig->second).GetOffset());  }
#if 0
    if(c==1) 
      std::cout << " c:" << c << " n:" << n << " a:" << a << " val:" << val
		<<" gain:" << (ig->second).GetGain() << " offset:" <<(ig->second).GetOffset()
		<<" cval:" << cval
		<< std::endl;
#endif
  }
  else{
    std::cout << " Invalid value!!! [" << funcname << "]"
	      << " c:" << c << " n:" << n << " a:" << a << " val:" << val
	      << std::endl;
    cval = val;
  }

  return cval;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
int GainMapMan::CalcDATValue( const int &c, const int &n, const int &a,
			      const int &at, const double &val )
{
  static const std::string funcname = "GainMapMan::CalcDATValue";
  unsigned int key;
  //double cval;
  int dat;

  key = KEY(c,n,a);

  GainMapContainer::iterator ig = gainContainer.find(key);
  if( ig != gainContainer.end() ){
    /* dE = gain(ADC-pedestal) */
    /* Time = gain*TDC-offset */
#if 0
    std::cout << " c:" << c << " n:" << n << " a:" << a << std::endl;
    std::cout << " gain:" << ((ig->second).GetGain()) << " offset:" << ((ig->second).GetOffset()) << std::endl;
#endif
    if( at==0 ) dat = (int)(val/((ig->second).GetGain()) + ((ig->second).GetOffset()));
    else        dat = (int)((val+((ig->second).GetOffset()))/((ig->second).GetGain()));
  }
  else{
#if 0
    std::cout << " Invalid value!!! [" << funcname << "]"
	      << " c:" << c << " n:" << n << " a:" << a << " val:" << val
	      << std::endl;
#endif
    dat = (int)val;
  }

  return dat;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void GainMapMan::PrintMapHeader( std::ostream &p_out )
{
  static const std::string funcname = "GainMapMan::PrintMapHeader";
  std::cout << " ---- " << funcname << " ---- " << std::endl;
  p_out << "#" << std::endl
	<< "# Gain & Offset ( but roughly ) map for TKO modules. ( Drt, HRTDC, ADC, etc.)" << std::endl
	<< "#" << std::endl
	<< "# for TDC: time[nsec] = par1*TDC + par2*TDC*TDC - par0" << std::endl
	<< "# for ADC: dE[MeV] = par1*(ADC - par0), offset = pedestal" << std::endl
	<< "#       or #photon = par1*(ADC - par0)" << std::endl
	<< "#" << std::endl;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void GainMapMan::PrintMap( const int &crate, std::ostream &p_out )
{
  static const std::string funcname = "GainMapMan::PrintMap";
  unsigned int key;
  GainMap again;
  int cr;
  int sl=0, sl_old=-1;

  std::cout << " ---- " << funcname << " ---- " << std::endl;

  for( GainMapContainer::const_iterator i=gainContainer.begin();
       i!=gainContainer.end(); i++ ){
    key = i->first;
    again = i->second;
    cr = ((key>>CSHIFT)&CMASK);
    sl = ((key>>NSHIFT)&NMASK);
    if( !( cr==crate || crate==-1 ) ) continue;
    if( sl!=sl_old ){
      p_out << "# crate  slot   ch    par0    par1    par2" << std::endl;
      sl_old=sl;
    }
    p_out.setf(std::ios::showpoint);
    p_out << std::setw(5) << cr << " "
	  << std::setw(7) << sl << " "
	  << std::setw(5) << ((key>>ASHIFT)&AMASK) << " "
	  << std::setprecision(5) << std::setw(10) << again.GetParam(0) << " "
	  << std::setprecision(5) << std::setw(10) << again.GetParam(1) << " "
	  << std::setprecision(5) << std::setw(10) << again.GetParam(2)
	  << std::endl;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void GainMapMan::PrintMapBL( std::ostream &p_out )
{
  static const std::string funcname = "GainMapMan::PrintMapBL";
  PrintMapHeader(p_out);
  PrintMap(0,p_out);
  PrintMap(1,p_out);
  PrintMap(2,p_out);
  PrintMap(7,p_out);
  PrintMap(8,p_out);
  PrintMap(9,p_out);
  PrintMap(10,p_out);
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void GainMapMan::PrintMapCDS( std::ostream &p_out )
{
  static const std::string funcname = "GainMapMan::PrintMapCDS";
  PrintMapHeader(p_out);
  PrintMap(3,p_out);
  PrintMap(4,p_out);
  PrintMap(5,p_out);
  PrintMap(6,p_out);
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void GainMapMan::PrintMapSDD( std::ostream &p_out )
{
  static const std::string funcname = "GainMapMan::PrintMapSDD";
  PrintMapHeader(p_out);
  PrintMap(3,p_out);
}
