// GateMapMan.cpp

#include <string>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <new>

#include "GateMapMan.h"
#include "GlobalVariables.h"
#include "TRandom.h"

ClassImp(GateMap);
ClassImp(GateMapMan);

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
GateMap::GateMap():
  Low(0),High(4095),Type(0)
{
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
GateMapMan::GateMapMan()
{
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
GateMapMan::GateMapMan( const GateMapMan &right )
{
  FileNameCDS = right.FileNameCDS;
  FileNameBL  = right.FileNameBL;
  FileNameSDD  = right.FileNameSDD;
  for( GateMapContainer::const_iterator i=right.gateContainer.begin();
       i!=right.gateContainer.end(); i++ ){
    gateContainer[i->first] = i->second;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
GateMapMan::~GateMapMan()
{
  gateContainer.clear();
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void GateMapMan::SetFileNameCDS( const std::string & filename )
{
  FileNameCDS = filename;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void GateMapMan::SetFileNameBL( const std::string & filename )
{
  FileNameBL = filename;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void GateMapMan::SetFileNameSDD( const std::string & filename )
{
  FileNameSDD = filename;
}

const unsigned int KEYMASK  = 0x000F;
// |1111|1110|1111|1011|1111|1110|1111|0011|
const unsigned int IMASK    = 0x000F;      /* I Mask 4 Bits (0-15) */
const unsigned int SMASK    = 0x01FF;      /* S Mask 9 Bits (0-511) */
const unsigned int UMASK    = 0x001F;      /* U Mask 5 Bits (0-31) */
const unsigned int CMASK    = 0x007F;      /* C Mask 7 Bits (0-127) */
const int          ISHIFT   =  4;
const int          SSHIFT   =  9;
const int          USHIFT   = 19;
const int          CSHIFT   = 25;
const unsigned int KEYFLAG  = 0x0003;

#define KEY(cid,seg,ud,i) \
((((cid)&CMASK)<<CSHIFT) | (((seg)&SMASK)<<SSHIFT) | (((ud)&UMASK)<<USHIFT) | (((i)&IMASK)<<ISHIFT) | KEYFLAG )

const int MAXCHAR = 144;

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool GateMapMan::Initialize()
{
  static const std::string funcname = "GateMapMan::Initialize";
  std::cout << "[" << funcname << "] Initialization start ... " ;

  int cid,seg,ud,at,lay,wire;
  int type;
  double low,high;
  int nd;
  unsigned int key;
  char str[MAXCHAR];
  FILE *fp;

  gateContainer.clear();

  if( FileNameCDS!=DefaultFileName ){
    if( (fp=fopen(FileNameCDS.c_str(), "r"))==0 ){
      std::cerr << " File open fail. [" << FileNameCDS << "]" << std::endl;
      exit(-1);
    }
    else{
      while( fgets(str,MAXCHAR,fp)!=0 ){
	if( str[0]=='#' ) continue;
	
	if( (nd=sscanf(str,"%d %d %d %d %d %lf %lf", &cid, &seg, &ud, &at, &type, &low, &high)) == 7 ) {
#if 0
	  std::cout << cid << "  " << seg << "  " << ud << "  " << at << "  "
		    << type << "  "<< low << "  " << high << std::endl;
#endif
	  key = KEY(cid,seg,ud,at);
	  GateMap *gate = new GateMap();
	  gate->SetLow(low);
	  gate->SetHigh(high);
	  gate->SetType(type);
	  gateContainer[key] = *gate;
	  delete gate;
	}
	else if( (nd=sscanf(str,"%d %d %d %d %lf %lf", &cid, &lay, &wire, &type, &low, &high)) == 6 ) {
#if 0
	  std::cout << cid << "  " << lay << "  " << wire << "  " 
		    << type << "  "<< low << "  " << high << std::endl;
#endif
	  key = KEY(cid,wire,lay,0);
	  GateMap *gate = new GateMap();
	  gate->SetLow(low);
	  gate->SetHigh(high);
	  gate->SetType(type);
	  gateContainer[key] = *gate;
	  delete gate;
	}
	else {
	  std::cerr << "[" << funcname << "]: Invalid data format" << std::endl;
	  std::cerr << std::string(str) << std::endl;
	}
      }
      fclose(fp);
      //      std::cout << "[" << funcname << "] CDS is initialized." << std::endl;
      std::cout << " CDS." ;//<< std::endl;
    }
  }

  if( FileNameBL!=DefaultFileName ){
    if( (fp=fopen(FileNameBL.c_str(), "r"))==0 ){
      std::cerr << " File open fail. [" << FileNameBL << "]" << std::endl;
      exit(-1);
    }
    else{
      while( fgets(str,MAXCHAR,fp)!=0 ){
	if( str[0]=='#' ) continue;
	
	if( (nd=sscanf(str,"%d %d %d %d %d %lf %lf", &cid, &seg, &ud, &at, &type, &low, &high)) == 7 ) {
#if 0
	  std::cout << cid << "  " << seg << "  " << ud << "  " << at << "  "
		    << type << "  "<< low << "  " << high << std::endl;
#endif
	  key = KEY(cid,seg,ud,at);
	  GateMap *gate = new GateMap();
	  gate->SetLow(low);
	  gate->SetHigh(high);
	  gate->SetType(type);
	  gateContainer[key] = *gate;
	  delete gate;
	}
	else if( (nd=sscanf(str,"%d %d %d %d %lf %lf", &cid, &lay, &wire, &type, &low, &high)) == 6 ) {
#if 0
	  std::cout << cid << "  " << lay << "  " << wire << "  " 
		    << type << "  "<< low << "  " << high << std::endl;
#endif
	  key = KEY(cid,wire,lay,0);
	  GateMap *gate = new GateMap();
	  gate->SetLow(low);
	  gate->SetHigh(high);
	  gate->SetType(type);
	  gateContainer[key] = *gate;
	  delete gate;
	}
	else {
	  std::cerr << "[" << funcname << "]: Invalid data format" << std::endl;
	  std::cerr << std::string(str) << std::endl;
	}
      }
      fclose(fp);
      //      std::cout << "[" << funcname << "] CDS is initialized." << std::endl;
      std::cout << " BL." ;//<< std::endl;
    }
  }

  std::cout <</* "[" << funcname << "] Initialization*/" finish." << std::endl;
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool GateMapMan::GetParam( const int &cid, const int &seg, const int &ud, const int &at, int &type, double &low, double &high )
{
  static const std::string funcname = "GateMapMan::GetParam";
  unsigned int key;
  GateMap gate;

  key = KEY(cid,seg,ud,at);
  GateMapContainer::iterator ig = gateContainer.find(key);
  if( ig != gateContainer.end() ){
    low = (ig->second).GetLow();
    high = (ig->second).GetHigh();
    type = (ig->second).GetType();
    return true;
  }
  else{
    return false;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool GateMapMan::GetParam( const int &cid, const int &lay, const int &wire, int &type, double &low, double &high )
{
  static const std::string funcname = "GateMapMan::GetParam";
  unsigned int key;
  GateMap agate;

  key = KEY(cid,wire,lay,0);
  GateMapContainer::iterator ig = gateContainer.find(key);
  if( ig != gateContainer.end() ){
    low = (ig->second).GetLow();
    high = (ig->second).GetHigh();
    type = (ig->second).GetType();
    return true;
  }
  else{
    return false;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool GateMapMan::CheckRange( const int &cid, const int &seg, const int &ud, const int &at, const int &type, const double &val )
{
  static const std::string funcname = "GateMapMan::GetParam";
  unsigned int key;
  GateMap agate;

  key = KEY(cid,seg,ud,at);
  GateMapContainer::iterator ig = gateContainer.find(key);
  if( ig != gateContainer.end() ){
    double low = (ig->second).GetLow();
    double high = (ig->second).GetHigh();
    int temptype = (ig->second).GetType();
    if(temptype==type){
#if 0
      if(cid==CID_MISC)
	std::cout << "[" << funcname << "]"
		  << " cid:" << cid << " seg:" << seg << " ud:" << ud << " at:" << at << " type:" << type
		<< " low:" << low << " high:" << high << " val:" << val
		<< std::endl;
#endif
      if(val>low&&val<high) return true;
      else return false;
    }else{
      if(val>0&&val<4095) return true;
      else return false;
    }
  }
  else{
    if(val>0&&val<4095) return true;
    else return false;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool GateMapMan::CheckRange( const int &cid, const int &lay, const int &wire, const int &type, const double &val )
{
  static const std::string funcname = "GateMapMan::GetParam";
  unsigned int key;
  GateMap agate;

  key = KEY(cid,wire,lay,0);
  GateMapContainer::iterator ig = gateContainer.find(key);
  if( ig != gateContainer.end() ){
    double low = (ig->second).GetLow();
    double high = (ig->second).GetHigh();
    int temptype = (ig->second).GetType();
    if(temptype==type){
      if(val>low&&val<high) return true;
      else return false;
    }else{
      if(val>0&&val<4095) return true;
      else return false;
    }
  }
  else{
    if(val>0&&val<4095) return true;
    else return false;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool GateMapMan::SetParam( const int &cid, const int &seg, const int &ud, const int &at, const int &type,
			    const double &low, const double &high )
{
  static const std::string funcname = "GateMapMan::SetParam";
  unsigned int key;

  key = KEY(cid,seg,ud,at);
  GateMapContainer::iterator ig = gateContainer.find(key);
  if( ig != gateContainer.end() ){
    (ig->second).SetLow( low );
    (ig->second).SetHigh( high );
    (ig->second).SetType( type );
    return true;
  }
  else{
    std::cout << " Invalid value!!! [" << funcname << "]"
	      << " cid:" << cid << " seg:" << seg << " ud:" << ud << " at:" << at
	      << std::endl;
    return false;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool GateMapMan::SetParam( const int &cid, const int &lay, const int &wire, const int &type,
			    const double &low, const double &high )
{
  static const std::string funcname = "GateMapMan::SetParam";
  unsigned int key;

  key = KEY(cid,wire,lay,0);
  GateMapContainer::iterator ig = gateContainer.find(key);
  if( ig != gateContainer.end() ){
    (ig->second).SetLow( low );
    (ig->second).SetHigh( high );
    (ig->second).SetType( type );
    return true;
  }
  else{
    std::cout << " Invalid value!!! [" << funcname << "]"
	      << " cid:" << cid << " layer:" << lay << " wire:" << wire
	      << std::endl;
    return false;
  }
}


// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void GateMapMan::PrintMapHeader( std::ostream &p_out )
{
  static const std::string funcname = "GateMapMan::PrintMapHeader";
  std::cout << " ---- " << funcname << " ---- " << std::endl;
  p_out << "#" << std::endl
	<< "# Gain & Offset ( but roughly ) map for TKO modules. ( Drt, HRTDC, ADC, etc.)" << std::endl
	<< "#" << std::endl
	<< "# for TDC: time[nsec] = gain*TDC - offset" << std::endl
	<< "# for ADC: dE[MeV] = gain*(ADC - pedestal), offset = pedestal" << std::endl
	<< "#       or #photon = gain*(ADC - pedestal)" << std::endl
	<< "#" << std::endl;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void GateMapMan::PrintMap( const int &crate, std::ostream &p_out )
{
  //  static const std::string funcname = "GateMapMan::PrintMap";
  //  unsigned int key;
  //  GateMap agate;
  //  int cr;
  //  int sl=0, sl_old=-1;

  //  std::cout << " ---- " << funcname << " ---- " << std::endl;
  /*
  for( GateMapContainer::const_iterator i=gateContainer.begin();
       i!=gateContainer.end(); i++ ){
    key = i->first;
    agate = i->second;
    cr = ((key>>CSHIFT)&CMASK);
    sl = ((key>>NSHIFT)&NMASK);
    if( !( cr==crate || crate==-1 ) ) continue;
    if( sl!=sl_old ){
      p_out << "# crate  slot   ch    gain    offset" << std::endl;
      sl_old=sl;
    }
    p_out.setf(std::ios::showpoint);
    p_out << std::setw(5) << cr << " "
	  << std::setw(7) << sl << " "
	  << std::setw(5) << ((key>>ASHIFT)&AMASK) << " "
	  << std::setprecision(5) << std::setw(9) << agate.GetLow() << " "
	  << std::setprecision(5) << std::setw(10) << again.GetHigh()
	  << std::endl;
  }
  */
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void GateMapMan::PrintMapBL( std::ostream &p_out )
{
  static const std::string funcname = "GateMapMan::PrintMapBL";
  PrintMapHeader(p_out);
  PrintMap(0,p_out);
  PrintMap(1,p_out);
  PrintMap(2,p_out);
  PrintMap(7,p_out);
  PrintMap(8,p_out);
  PrintMap(9,p_out);
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void GateMapMan::PrintMapCDS( std::ostream &p_out )
{
  static const std::string funcname = "GateMapMan::PrintMapCDS";
  PrintMapHeader(p_out);
  PrintMap(3,p_out);
  PrintMap(4,p_out);
  PrintMap(5,p_out);
  PrintMap(6,p_out);
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void GateMapMan::PrintMapSDD( std::ostream &p_out )
{
  static const std::string funcname = "GateMapMan::PrintMapSDD";
  PrintMapHeader(p_out);
  PrintMap(3,p_out);
}
