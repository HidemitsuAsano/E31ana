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
void BLDCWireMap::SetParam(double *param)
{
#if 0
  for(int i=0;i<9;i++)
    std::cout<<param[i]<<"  ";
  std::cout<<std::endl;
#endif
  nWire = (int)param[0];
  Z = param[1];
  XY = (int)param[2]; XY0 = param[3]; dXY = param[4];
  WireLength = param[5]; TiltAngle = param[6];
  RotationAngle = param[7]; WirePhi = param[8];
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
ClassImp(BLDCWireMapMan);
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
BLDCWireMapMan::BLDCWireMapMan()
{
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool BLDCWireMapMan::Initialize()
{
  static const std::string funcname = "BLDCWireMapMan::Initialize";
  std::cout << "[" << funcname << "] Initialization start ...";
  
  std::string filename=FileName; 
  bldcContainer.clear();
  geomContainer.clear();

  const int MAXCHAR = 512;
  const int MAXTOKEN = 20;
  const char* DELIMITER = " ";

  int cid=-1;
  int layer=-1;
  const int nseg=20;
  double par[20];
  unsigned int key;
  char buf[MAXCHAR];
  ifstream fin;
  
  if( filename!=DefaultFileName ){
    fin.open(filename.c_str());
    if (!fin.good()){
      std::cerr << " File open fail. [" << filename << "]" << std::endl;
      exit(-1);
    }
    while (!fin.eof()){
      fin.getline(buf, MAXCHAR);
      int n = 0;
      const char* token[MAXTOKEN] = {};
      token[0] = strtok(buf, DELIMITER);
      if (token[0] != NULL) {
	if( !strcmp(token[0], "#") ) continue;
	for (n = 1; n < MAXTOKEN; n++){
	  token[n] = strtok( NULL, DELIMITER);
	  if (token[n] == NULL ) break;
	}
      }
      for (int i = 0; i < n; i++){
	if(i==0) cid=atoi(token[i]);
	if(i==1) layer=atoi(token[i]);
	if(i >= 2 && i < nseg+2 )
	  par[i-2]=atof(token[i]);
      }
      if(n==11+2){
	GeomMap *ageom = new GeomMap();
	ageom->SetParam(par);
#if 0
	std::cout<<cid<<"  "<<layer;
	ageom->PrintMap();
#endif 
	key = KEY( cid, layer );
	geomContainer[key] = *ageom;
	delete ageom;
      }else if(n==9+2){
	BLDCWireMap *awire = new BLDCWireMap();
	awire->SetParam(par);
#if 0
	std::cout<<cid<<"  "<<layer;
	awire->PrintMap();
#endif 
	key = KEY( cid, layer );
	bldcContainer[key] = *awire;	  
	delete awire;
      }
    }
    fin.close();
  }
  //  std::cout << "[" << funcname << "] Initialization finish." << std::endl;
  std::cout << " finish." << std::endl;
  return true;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
GeomMap * BLDCWireMapMan::GetGMap( const int &cid, const int &seg)
{
  static const std::string funcname = "BLDCWireMapMan::GetGMap";  
  unsigned int key;
  key = KEY( cid, seg );
  GeomMapContainer::iterator ic = geomContainer.find(key);
  if( ic != geomContainer.end() ) return &(ic->second);
  else{
    std::cout << " Invalid value!!! [" << funcname << "]"
	      << " cid:" << cid
	      << " seg:" << seg
	      << std::endl;
    return NULL;
  }
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
BLDCWireMap * BLDCWireMapMan::GetWireMap( const int &cid, const int &layer)
{
  static const std::string funcname = "BLDCWireMapMan::GetWireMap";  
  unsigned int key;
  key = KEY( cid, layer );
  BLDCWireMapContainer::iterator ic = bldcContainer.find(key);
  if( ic != bldcContainer.end() ) return &(ic->second);
  else{
    std::cout << " Invalid value!!! [" << funcname << "]"
	      << " cid:" << cid
	      << " layer:" << layer
	      << std::endl;
    return NULL;
  }
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool BLDCWireMapMan::GetGParam( const int &cid, const int &seg, double *par)
{
  GeomMap *map=GetGMap(cid,seg);
  if(!map) return false;
  map->GetParam(par);
  return true;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool BLDCWireMapMan::SetGParam( const int &cid, const int &seg,const double *par)
{
  GeomMap *map=GetGMap(cid,seg);
  if(!map) return false;
  map->SetParam(par);
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool BLDCWireMapMan::GetParam( const int &cid, const int &layer,
			       int &nw, double &z, int &xy, double &xy0, double &dxy, 
			       double &wl, double &tilt, double &ra )
{
  static const std::string funcname = "BLDCWireMapMan::GetParam";
  BLDCWireMap *map=GetWireMap(cid,layer);
  if(!map) return false;
  nw = map->GetNWire();
  z  = map->GetZ();
  xy = map->GetXY();
  xy0= map->GetXY0();
  dxy= map->GetdXY();
  wl = map->GetWireLength();
  tilt=map->GetTiltAngle();
  ra = map->GetRotationAngle();
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool BLDCWireMapMan::GetXY0( const int &cid, const int &layer, double &xy0)
{
  static const std::string funcname = "BLDCWireMapMan::GetXY0";
  BLDCWireMap *map=GetWireMap(cid,layer);
  if(!map) return false;
  xy0 = map->GetXY0();
  return true;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool BLDCWireMapMan::SetXY0( const int &cid, const int &layer,const double &xy0)
{
  static const std::string funcname = "BLDCWireMapMan::SetXY0";
  BLDCWireMap *map=GetWireMap(cid,layer);
  if(!map) return false;
  map->SetXY0(xy0);
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool BLDCWireMapMan::SetParam( const int &cid, const int &layer,
			       const int &nw, const double &z, const int &xy, const double &xy0,const double &dxy, 
			       const double &wl, const double &tilt, const double &ra )
{
  static const std::string funcname = "BLDCWireMapMan::SetParam";
  BLDCWireMap *map=GetWireMap(cid,layer);
  if(!map) return false;
  map->SetParam(nw,z,xy,xy0,dxy,wl,tilt,ra);
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
int BLDCWireMapMan::GetNWire( const int &cid, const int &layer )
{
  static const std::string funcname = "BLDCWireMapMan::GetParam";
  int nw=-1;
  BLDCWireMap *map=GetWireMap(cid,layer);
  if(map)  nw = map->GetNWire();
  return nw;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool BLDCWireMapMan::GetGParam( const int &cid, TVector3 &pos, TVector3 &rot)
{
  static const std::string funcname = "BLDCWireMapMan::GetGParam";
  GeomMap *map=GetGMap(cid,0);
  if(!map) return false;
  pos=map->GetPos();
  rot=map->GetRot();
  return true;
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
       i!=bldcContainer.end(); ++i ){
    key = i->first;
    map = i->second;
    cid = ((key>>CSHIFT)&CMASK);
    if( !( cid==Cid || Cid==-1 ) ) continue;
    if( cid!=cid_old ){
      p_out<< "# CID  Seg   X         Y         Z        rotX       rotY       rotZ        "
	   << "SizeX(W)  SizeY(L)  SizeZ(T)  Size4  LightV"<<std::endl;
      p_out<< "#            [cm]      [cm]      [cm]     [deg]      [deg]      [deg]       "
	   << "[cm]      [cm]      [cm]             [cm/ns]"<<std::endl;
      p_out.setf(std::ios::showpoint);
      p_out << std::setw(6) << cid
	    << std::setw(5) << 0;
      GetGMap(cid,0)->PrintMap(p_out);
      p_out	<< "# CID  Layer    nwire    z      xy     x0/y0   dx/dy  wirelength   tilt      rotateangle   wirephi"<<std::endl;
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
	  << std::setw(8) << map.GetWirePhi()
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
