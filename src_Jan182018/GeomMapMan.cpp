// GeomMapMan.cpp

#include <string>
#include <cstdio>
#include <iostream>
#include <new>
#include <cmath>

#include "GeomMapMan.h"
#include "GlobalVariables.h"

ClassImp(GeomMap);
ClassImp(GeomMapMan);

const double Deg2Rad = 3.141592/180.;
const double Rad2Deg = 180./3.141592;

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
GeomMap::GeomMap()
{
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
GeomMapMan::GeomMapMan()
{
  FileNameCDS = DefaultFileName;
  FileNameBL  = DefaultFileName;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
GeomMapMan::GeomMapMan( const GeomMapMan &right )
{
  FileNameCDS = right.FileNameCDS;
  FileNameBL  = right.FileNameBL;
  for( GeomMapContainer::const_iterator i=right.geomContainer.begin();
       i!=right.geomContainer.end(); i++ ){
    geomContainer[i->first] = i->second;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
GeomMapMan::~GeomMapMan()
{
  geomContainer.clear();
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void GeomMapMan::SetFileNameCDS( const std::string & filename )
{
  FileNameCDS = filename;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void GeomMapMan::SetFileNameBL( const std::string & filename )
{
  FileNameBL = filename;
}

const unsigned int KEYMASK  = 0x000F;
const unsigned int SMASK    = 0x00FF;
const unsigned int CMASK    = 0x00FF;
const int          SSHIFT   = 4;
const int          CSHIFT   = 12;
const unsigned int KEYFLAG  = 0x0003;
#define KEY(cid,seg) \
((((cid)&CMASK)<<CSHIFT) | (((seg)&SMASK)<<SSHIFT) | KEYFLAG )

const int MAXCHAR = 255;

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool GeomMapMan::Initialize()
{
  static const std::string funcname = "GeomMapMan::Initialize";
  std::cout << "[" << funcname << "] Initialization start ...";// << std::endl;

  int cid, seg;
  double x, y, z, dx, dy, dz, len, wid, thick;
  double lv;
  int gcid;
  double gx=-999.,gy=-999.,gz=-999.,gdx=-999.,gdy=-999.,gdz=-999.;
  unsigned int key;
  int nd;
  char str[MAXCHAR];
  FILE *fp;

  geomContainer.clear();

  if( FileNameCDS!=DefaultFileName ){
    if( (fp=fopen(FileNameCDS.c_str(), "r"))==0 ){
      std::cerr << " File open fail. [" << FileNameCDS << "]" << std::endl;
      exit(-1);
    }
    gcid=-1;
    while( fgets(str,MAXCHAR,fp)!=0 ){
      if( str[0]=='#' ) continue;
      
      if( (nd=sscanf(str,"GPOS: %d %lf %lf %lf %lf %lf %lf",&cid,&x,&y,&z,&dx,&dy,&dz)) == 7 ){
	gcid=cid; gx=x; gy=y; gz=z; gdx=dx; gdy=dy; gdz=dz;
      }
      else if( (nd=sscanf(str,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", 
			  &cid,&seg,&x,&y,&z,&dx,&dy,&dz,&len,&wid,&thick,&lv )) == 12 ) {
#if 0
	std::cout << cid << "  " << seg << "  " << x << "  " << y << "  " << z << "  "
		  << dx << "  " << dy << "  " << dz << "  "
		  << len << "  " << wid << "  " << thick << "  " << lv
		  << std::endl;
#endif
	GeomMap *ageom = new GeomMap();
	ageom->SetParam(x,y,z,dx,dy,dz,len,wid,thick,lv);
	if( cid==gcid ) ageom->SetGParam( gx,gy,gz,gdx,gdy,gdz );
	key = KEY( cid, seg );
	geomContainer[key] = *ageom;
	delete ageom;
      }
      else{
	std::cerr << "[" << funcname << "]: Invalid data format" << std::endl;
	std::cerr << FileNameCDS << std::endl;
	std::cerr << std::string(str) << std::endl;
      }
    }
    fclose(fp);
    //    std::cout << "[" << funcname << "] CDS is initialized." << std::endl;
    std::cout <<" CDS. " ;
  }

  if( FileNameBL!=DefaultFileName ){
    if( (fp=fopen(FileNameBL.c_str(), "r"))==0 ){
      std::cerr << " File open fail. [" << FileNameBL << "]" << std::endl;
      exit(-1);
    }
    gcid=-1;
    while( fgets(str,MAXCHAR,fp)!=0 ){
      if( str[0]=='#' ) continue;
      
      if( (nd=sscanf(str,"GPOS: %d %lf %lf %lf %lf %lf %lf",&cid,&x,&y,&z,&dx,&dy,&dz))==7 ){
	gcid=cid; gx=x; gy=y;gz=z; gdx=dx; gdy=dy; gdz=dz;
      }
      else if( (nd=sscanf(str,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", 
			  &cid,&seg,&x,&y,&z,&dx,&dy,&dz,&len,&wid,&thick,&lv )) == 12 ) {
#if 0
	if(cid==92){
	  std::cout << cid << "  " << seg << "  " << x << "  " << y << "  " << z << "  "
		    << dx << "  " << dy << "  " << dz << "  "
		    << len << "  " << wid << "  " << thick << "  " << lv
		    << std::endl;
	}
#endif
	GeomMap *ageom = new GeomMap();
	ageom->SetParam(x,y,z,dx,dy,dz,len,wid,thick,lv);
	if( cid==gcid ) ageom->SetGParam( gx,gy,gz,gdx,gdy,gdz );
	key = KEY( cid, seg );
	geomContainer[key] = *ageom;
	delete ageom;
      }
      else{
	std::cerr << "[" << funcname << "]: Invalid data format" << std::endl;
	std::cerr << FileNameBL << std::endl;
	std::cerr << std::string(str) << std::endl;
      }
    }
    fclose(fp);
    //    std::cout << "[" << funcname << "] BeamLine is initialized." << std::endl;
    std::cout << " BeamLine. ";
  }

  std::cout <</* "[" << funcname << "] Initialization*/" finish." << std::endl;
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool GeomMapMan::GetParam( const int &cid, const int &seg, double &x, double &y, double &z,
			   double &dx, double &dy, double &dz,
			   double &len, double &wid, double &th )
{
  static const std::string funcname = "GeomMapMan::GetParam";
  unsigned int key;

  key = KEY( cid, seg );
  if(cid==CID_T0pre||cid==CID_T0post)
    key = KEY( CID_T0, seg );
  if(cid==CID_BHDpost)
    key = KEY( CID_BHD, seg );
  GeomMapContainer::iterator ic = geomContainer.find(key);
  if( ic != geomContainer.end() ){
    x = (ic->second).GetX();
    y = (ic->second).GetY();
    z = (ic->second).GetZ();
    dx = (ic->second).GetdX();
    dy = (ic->second).GetdY();
    dz = (ic->second).GetdZ();
    len  = (ic->second).GetLength();
    wid  = (ic->second).GetWidth();
    th   = (ic->second).GetThick();
    return true;
  }
  else{
    std::cout << " Invalid value!!! [" << funcname << "]"
	      << " cid:" << cid
	      << " seg:" << seg
	      << std::endl;
    return false;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool GeomMapMan::GetLightVelocity( const int &cid, const int &seg, double &lv )
{
  static const std::string funcname = "GeomMapMan::GetParam";
  unsigned int key;
  
  key = KEY( cid, seg );
  GeomMapContainer::iterator ic = geomContainer.find(key);
  if( ic != geomContainer.end() ){
    lv   = (ic->second).GetLightVelocity();
    return true;
  }
  else{
    std::cout << " Invalid value!!! [" << funcname << "]"
	      << " cid:" << cid
	      << " seg:" << seg
	      << std::endl;
    return false;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool GeomMapMan::SetLightVelocity( const int &cid, const int &seg, const double &lv )
{
  static const std::string funcname = "GeomMapMan::GetParam";
  unsigned int key;
  
  key = KEY( cid, seg );
  GeomMapContainer::iterator ic = geomContainer.find(key);
  if( ic != geomContainer.end() ){
    (ic->second).SetLightVelocity(lv);
    return true;
  }
  else{
    std::cout << " Invalid value!!! [" << funcname << "]"
	      << " cid:" << cid
	      << " seg:" << seg
	      << std::endl;
    return false;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool GeomMapMan::GetParam( const int &cid, const int &seg, double &x, double &y, double &z,
			   double &dx, double &dy, double &dz,
			   double &len, double &wid, double &th, double &lv )
{
  static const std::string funcname = "GeomMapMan::GetParam";
  unsigned int key;

  key = KEY( cid, seg );
  if(cid==51||cid==52)
    key = KEY( CID_T0, seg );
  if(cid==56)
    key = KEY( CID_BHD, seg );
  GeomMapContainer::iterator ic = geomContainer.find(key);
  if( ic != geomContainer.end() ){
    x = (ic->second).GetX();
    y = (ic->second).GetY();
    z = (ic->second).GetZ();
    dx = (ic->second).GetdX();
    dy = (ic->second).GetdY();
    dz = (ic->second).GetdZ();
    len  = (ic->second).GetLength();
    wid  = (ic->second).GetWidth();
    th   = (ic->second).GetThick();
    lv   = (ic->second).GetLightVelocity();
    return true;
  }
  else{
    std::cout << " Invalid value!!! [" << funcname << "]"
	      << " cid:" << cid
	      << " seg:" << seg
	      << std::endl;
    return false;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool GeomMapMan::SetParam( const int &cid, const int &seg, const double &x, const double &y, const double &z,
			   const double &dx, const double &dy, const double &dz,
			   const double &len, const double &wid, const double &th, const double &lv )
{
  static const std::string funcname = "GeomMapMan::SetParam";
  unsigned int key;

  key = KEY( cid, seg );
  GeomMapContainer::iterator ic = geomContainer.find(key);
  if( ic != geomContainer.end() ){
    (ic->second).SetX(x);
    (ic->second).SetY(y);
    (ic->second).SetZ(z);
    (ic->second).SetdX(dx);
    (ic->second).SetdY(dy);
    (ic->second).SetdZ(dz);
    (ic->second).SetLength(len);
    (ic->second).SetWidth(wid);
    (ic->second).SetThick(th);
    (ic->second).SetLightVelocity(lv);
    return true;
  }
  else{
    std::cout << " Invalid value!!! [" << funcname << "]"
	      << " cid:" << cid
	      << " seg:" << seg
	      << std::endl;
    return false;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool GeomMapMan::GetGParam( const int &cid, double &x, double &y, double &z,
			    double &dx, double &dy, double &dz )
{
  static const std::string funcname = "GeomMapMan::GetGParam";
  unsigned int key;

  key = KEY( cid, 1 );
  if(cid==51||cid==52)
    key = KEY( CID_T0, 1 );
  if(cid==56)
    key = KEY( CID_BHD, 1 );
  GeomMapContainer::iterator ic = geomContainer.find(key);
  if( ic != geomContainer.end() ){
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
void GeomMapMan::PrintMap( const int &Cid, std::ostream &p_out )
{
  static const std::string funcname = "GeomMapMan::PrintMap";
  unsigned int key;
  GeomMap ageom;
  int cid=0, cid_old=-1;

  std::cout << " ---- " << funcname << " ---- " << std::endl;
  for( GeomMapContainer::const_iterator i=geomContainer.begin();
       i!=geomContainer.end(); i++ ){
    key = i->first;
    ageom = i->second;
    cid = ((key>>CSHIFT)&CMASK);
    if( !( cid==Cid || Cid==-1 ) ) continue;
    if( cid!=cid_old ){
      p_out<<"#         CID    x[cm]  y[cm]  z[cm]    dx    dy    dz"<<std::endl;
      p_out.setf(std::ios_base::fixed,std::ios_base::floatfield);
      p_out.setf(std::ios::showpoint);
      p_out<<"GPOS: "
	   <<std::setprecision(3)
	   <<std::setw(10)<< cid	 
	   <<std::setw(10)<< ageom.GetGX()	 
	   <<std::setw(10)<< ageom.GetGY()
	   <<std::setw(10)<< ageom.GetGZ()
	   <<std::setw(10)<< ageom.GetdGX()
	   <<std::setw(10)<< ageom.GetdGY()
	   <<std::setw(10)<< ageom.GetdGZ()
	   <<std::endl;
      p_out<< "# CID  Seg   X         Y         Z         dX        dY        dZ        "
	   << "Length    Width     Thick     LightV"<<std::endl;
      p_out<< "#            [cm]      [cm]      [cm]                                    "
	   << "[cm]      [cm]      [cm]      [cm/ns]"<<std::endl;
      cid_old = cid;
    }
    ageom = i->second;
    p_out << std::setw(6) << ((key>>CSHIFT)&CMASK)
	  << std::setw(5) << ((key>>SSHIFT)&SMASK)
	  << std::setw(10) << ageom.GetX()
	  << std::setw(10) << ageom.GetY()
	  << std::setw(10) << ageom.GetZ()
	  << std::setw(10) << ageom.GetdX()
	  << std::setw(10) << ageom.GetdY()
	  << std::setw(10) << ageom.GetdZ()
	  << std::setw(10) << ageom.GetLength()
	  << std::setw(10) << ageom.GetWidth()
	  << std::setw(10) << ageom.GetThick()
	  << std::setw(10) << ageom.GetLightVelocity()
	  << std::endl;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void GeomMapMan::PrintMapBL( std::ostream &p_out )
{
  static const std::string funcname = "GeomMapMan::PrintMapBL";
  PrintMap(CID_BHD,p_out);
  PrintMap(CID_PA,p_out);
  PrintMap(CID_T0,p_out);
  PrintMap(CID_E0,p_out);
  PrintMap(CID_B1,p_out);
  PrintMap(CID_B2,p_out);
  PrintMap(CID_LC1,p_out);
  PrintMap(CID_LC2,p_out);
  PrintMap(CID_AC,p_out);
  PrintMap(CID_WC,p_out);
  PrintMap(CID_GC,p_out);
  PrintMap(CID_Range,p_out);
  PrintMap(CID_BVC,p_out);
  PrintMap(CID_NC,p_out);
  PrintMap(CID_CVC,p_out);
  PrintMap(CID_LB,p_out);
  PrintMap(CID_PC,p_out);
  PrintMap(CID_BPD,p_out);
  PrintMap(CID_BD,p_out);
  PrintMap(CID_VCC,p_out);
}
