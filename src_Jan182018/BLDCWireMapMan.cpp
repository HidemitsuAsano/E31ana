// BLDCWireMapMan.cpp

#include <new>
#include <cmath>
#include <iomanip>

#include "BLDCWireMapMan.h"
#include "GlobalVariables.h"

ClassImp(BLDCWireMap);
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
BLDCWireMap::BLDCWireMap()
{
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void BLDCWireMap::SetParam( const int &nw, const double &z, const int &xy, const double &xy0, const double &dxy,
			    const double &wl, const double &tilt, const double &ra )
{
#if 0
      std::cout << nw << "  " << z << "  "
		<< xy << "  " << xy0 << "  " << dxy << "  " << wl << "  " << tilt << "  " << ra << std::endl;
#endif

  nWire = nw;
  Z = z;
  XY = xy; XY0 = xy0; dXY = dxy;
  WireLength = wl; TiltAngle = tilt; RotationAngle = ra;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void BLDCWireMap::SetGParam( const double  &x, const double  &y, const double  &z, 
			     const double &dx, const double &dy, const double &dz )
{
  GX = x; GY = y; GZ = z; dGX = dx; dGY = dy; dGZ = dz;
}

ClassImp(BLDCWireMapMan);
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
BLDCWireMapMan::BLDCWireMapMan()
{
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
BLDCWireMapMan::BLDCWireMapMan( const BLDCWireMapMan &right )
{
  FileName = right.FileName;
  for( BLDCWireMapContainer::const_iterator i = right.bldcContainer.begin();
       i!=right.bldcContainer.end(); i++ ){
    bldcContainer[i->first] = i->second;
  }
}

const unsigned int KEYMASK = 0x000F;
const unsigned int CMASK   = 0x00FF;
const unsigned int LMASK   = 0x00FF;
const int          CSHIFT  = 4;
const int          LSHIFT  = 16;
const unsigned int KEYFLAG = 0x0003;
#define KEY(cid,layer) \
((((cid)&CMASK)<<CSHIFT) | (((layer)&LMASK)<<LSHIFT) | KEYFLAG )

const int MAXCHAR = 144;

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool BLDCWireMapMan::Initialize()
{
  static const std::string funcname = "BLDCWireMapMan::Initialize";
  std::cout << "[" << funcname << "] Initialization start ...";

  int cid, layer;
  int nd;
  int nw, xy;
  double xy0, dxy, wl, tilt, ra;
  double x,y,z,dx,dy,dz;
  double gx,gy,gz,dgx,dgy,dgz;
  unsigned int key;
  
  char str[MAXCHAR];
  FILE *fp;

  bldcContainer.clear();

  if( (fp=fopen(FileName.c_str(), "r"))==0 ){
    std::cerr << " File open fail. [" << FileName << "]" << std::endl;
    exit(-1);
  }
  while( fgets(str,MAXCHAR,fp)!=0 ){
    if( str[0]=='#' ) continue;

    if( (nd=sscanf(str,"GPOS: %lf %lf %lf %lf %lf %lf",&x,&y,&z,&dx,&dy,&dz)) == 6 ){
      gx=x; gy=y; gz=z; dgx=dx; dgy=dy; dgz=dz;
#if 0
      std::cout << gx << "  " << gy << "  " << gz << "  " << dgx << "  "
		<< dgy << "  " << dgz << std::endl;
#endif

    }
    else if( (nd=sscanf(str,"%d %d %d %lf %d %lf %lf %lf %lf %lf", 
		   &cid,&layer,&nw,&z,&xy,&xy0,&dxy,&wl,&tilt,&ra )) == 10 ) {
#if 0
      std::cout << cid << "  " << layer << "  " << nw << "  " << z << "  "
		<< xy << "  " << xy0 << "  " << dxy << "  " << wl << "  " << tilt << "  " << ra << std::endl;
#endif
      BLDCWireMap *amap = new BLDCWireMap();
      amap->SetParam( nw, z, xy, xy0, dxy, wl, tilt, ra );
      amap->SetGParam( gx,gy,gz,dgx,dgy,dgz );
      key = KEY( cid, layer );
      bldcContainer[key] = *amap;
      delete amap;
    }
    else{
      std::cerr << "[" << funcname << "]: Invalid data format" << std::endl;
      std::cerr << std::string(str) << std::endl;
    }
  }
  fclose(fp);

  //PrintMap();

  //  std::cout << "[" << funcname << "] Initialization finish." << std::endl;
  std::cout << " finish." << std::endl;
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool BLDCWireMapMan::GetParam( const int &cid, const int &layer,
			       int &nw, double &z, int &xy, double &xy0, double &dxy, 
			       double &wl, double &tilt, double &ra )
{
  static const std::string funcname = "BLDCWireMapMan::GetParam";
  unsigned int key;

  key = KEY( cid, layer );
  BLDCWireMapContainer::iterator ic = bldcContainer.find(key);
  if( ic != bldcContainer.end() ){
    nw = (ic->second).GetNWire();
    z  = (ic->second).GetZ();
    xy = (ic->second).GetXY();
    xy0= (ic->second).GetXY0();
    dxy= (ic->second).GetdXY();
    wl = (ic->second).GetWireLength();
    tilt=(ic->second).GetTiltAngle();
    ra = (ic->second).GetRotationAngle();
    return true;
  }
  else{
    std::cout << " Invalid value!!! [" << funcname << "]"
	      << " cid:" << cid
	      << " layer:" << layer
	      << std::endl;
    return false;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool BLDCWireMapMan::GetXY0( const int &cid, const int &layer, double &xy0)
{
  static const std::string funcname = "BLDCWireMapMan::GetXY0";
  unsigned int key;

  key = KEY( cid, layer );
  BLDCWireMapContainer::iterator ic = bldcContainer.find(key);
  if( ic != bldcContainer.end() ){
    xy0= (ic->second).GetXY0();
    return true;
  }
  else{
    std::cout << " Invalid value!!! [" << funcname << "]"
	      << " cid:" << cid
	      << " layer:" << layer
	      << std::endl;
    return false;
  }
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool BLDCWireMapMan::SetXY0( const int &cid, const int &layer,const double &xy0)
{
  static const std::string funcname = "BLDCWireMapMan::SetXY0";
  unsigned int key;

  key = KEY( cid, layer );
  BLDCWireMapContainer::iterator ic = bldcContainer.find(key);
  if( ic != bldcContainer.end() ){
    (ic->second).SetXY0(xy0);
    return true;
  }
  else{
    std::cout << " Invalid value!!! [" << funcname << "]"
	      << " cid:" << cid
	      << " layer:" << layer
	      << std::endl;
    return false;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool BLDCWireMapMan::SetParam( const int &cid, const int &layer,
			       const int &nw, const double &z, const int &xy, const double &xy0,const double &dxy, 
			       const double &wl, const double &tilt, const double &ra )
{
  static const std::string funcname = "BLDCWireMapMan::SetParam";
  unsigned int key;

  key = KEY( cid, layer );
  BLDCWireMapContainer::iterator ic = bldcContainer.find(key);
  if( ic != bldcContainer.end() ){
    (ic->second).SetNWire(nw);
    (ic->second).SetZ(z);
    (ic->second).SetXY(xy);
    (ic->second).SetXY0(xy0);
    (ic->second).SetdXY(dxy);
    (ic->second).SetWireLength(wl);
    (ic->second).SetTiltAngle(tilt);
    (ic->second).SetRotationAngle(ra);
    return true;
  }
  else{
    std::cout << " Invalid value!!! [" << funcname << "]"
	      << " cid:" << cid
	      << " layer:" << layer
	      << std::endl;
    return false;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
int BLDCWireMapMan::GetNWire( const int &cid, const int &layer )
{
  static const std::string funcname = "BLDCWireMapMan::GetParam";
  unsigned int key;
  int nw=-1;

  key = KEY( cid, layer );
  BLDCWireMapContainer::iterator ic = bldcContainer.find(key);
  if( ic != bldcContainer.end() ){
    nw = (ic->second).GetNWire();
    //return true;
  }
  else{
    std::cout << " Invalid value!!! [" << funcname << "]"
	      << " cid:" << cid
	      << " layer:" << layer
	      << std::endl;
    //return false;
  }
  return nw;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool BLDCWireMapMan::GetGParam( const int &cid, double &x, double &y, double &z,
				double &dx, double &dy, double &dz )
{
  static const std::string funcname = "BLDCWireMapMan::GetGParam";
  unsigned int key;

  key = KEY( cid, 1 );
  BLDCWireMapContainer::iterator ic = bldcContainer.find(key);
  if( ic != bldcContainer.end() ){
    x = (ic->second).GetGX();
    y = (ic->second).GetGY();
    z = (ic->second).GetGZ();
    dx = (ic->second).GetdGX();
    dy = (ic->second).GetdGY();
    dz = (ic->second).GetdGZ();
    return true;
  }
  else{
    std::cout << " Invalid value!!! [" << funcname << "]"
	      << " cid:" << cid
	      << std::endl;
    return false;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void BLDCWireMapMan::PrintMap( const int &Cid, std::ostream &p_out )
{
  static const std::string funcname = "BLDCWireMapMan::PrintMap";
  unsigned int key;
  BLDCWireMap map;

  std::cout << " ---- " << funcname << " ---- " << std::endl;  
  int cid, cid_old=-1;
  
  for( BLDCWireMapContainer::const_iterator i=bldcContainer.begin();
       i!=bldcContainer.end(); i++ ){
    key = i->first;
    map = i->second;
    cid = ((key>>CSHIFT)&CMASK);
    if( !( cid==Cid || Cid==-1 ) ) continue;
    if( cid!=cid_old ){
      p_out<<"#         x[cm] y[cm] z[cm]    dx    dy    dz"<<std::endl;
      p_out.setf(std::ios::showpoint);
      p_out<<"GPOS: "
	   <<std::setprecision(4)
	   <<std::setw(10)<< map.GetGX()	 
	   <<std::setw(10)<< map.GetGY()
	   <<std::setw(10)<< map.GetGZ()
	   <<std::setw(10)<< map.GetdGX()
	   <<std::setw(10)<< map.GetdGY()
	   <<std::setw(10)<< map.GetdGZ()
	   <<std::endl;
      p_out	<< "# CID  Layer    nwire    z      xy     x0/y0   dx/dy  wirelength   tilt      rotateangle"<<std::endl;
      cid_old = cid;
    }
    p_out << std::setw(5) << cid
	  << std::setw(8) << ((key>>LSHIFT)&LMASK)
	  << std::setw(8) << map.GetNWire()
	  <<std::setprecision(4)
	  << std::setw(8) << map.GetZ()
	  << std::setw(8) << map.GetXY()
	  <<std::setprecision(5)
	  << std::setw(10) << map.GetXY0()
	  <<std::setprecision(2)
	  << std::setw(8) << map.GetdXY()
	  <<std::setprecision(3)
	  << std::setw(8) << map.GetWireLength()
	  <<std::setprecision(5)
	  << std::setw(8) << map.GetTiltAngle()
	  << std::setw(8) << map.GetRotationAngle()
	  << std::endl;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void BLDCWireMapMan::PrintMapBL( std::ostream &p_out )
{
  static const std::string funcname = "BLDCWireMapMan::PrintMapBL";
  //  PrintMapHeader(p_out);
  PrintMap(CID_BLC1a,p_out);
  PrintMap(CID_BLC1b,p_out);
  PrintMap(CID_BLC2a,p_out);
  PrintMap(CID_BLC2b,p_out);
  PrintMap(CID_BPC,p_out);
}
