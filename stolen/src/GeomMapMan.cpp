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

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
GeomMap::GeomMap()
{
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void GeomMap::PrintMap( std::ostream &p_out ){
  double tmppar[npar];
  GetParam(tmppar);
  for(int i=0;i<npar;i++)  
    p_out<<std::setw(10)<<"   "<<tmppar[i];
  p_out<<std::endl;
}
  
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
GeomMapMan::GeomMapMan()
{
  FileNameCDS = DefaultFileName;
  FileNameBL  = DefaultFileName;
  FileNameHall  = DefaultFileName;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
GeomMapMan::~GeomMapMan()
{
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
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void GeomMapMan::SetFileNameHall( const std::string & filename )
{
  FileNameHall = filename;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool GeomMapMan::Initialize()
{
  static const std::string funcname = "GeomMapMan::Initialize";
  std::cout << "[" << funcname << "] Initialization start ...";// << std::endl;
  geomContainer.clear();
  ReadFile(FileNameCDS);
  ReadFile(FileNameBL);
  ReadFile(FileNameHall);
  //  std::cout << "[" << funcname << "] Initialization finish." << std::endl;
  std::cout << " finish." << std::endl;
  return true;
}

bool GeomMapMan::ReadFile( const std::string filename )
{
  const int MAXCHAR = 512;
  const int MAXTOKEN = 20;
  const char* DELIMITER = " ";

  int cid=-1;
  int seg=-1;
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
      if(n!=npar+2) continue;
      for (int i = 0; i < n; i++){
	if(i==0) cid=atoi(token[i]);
	if(i==1) seg=atoi(token[i]);
	if(i >= 2 && i < nseg+2 )
	  par[i-2]=atof(token[i]);
      }
      GeomMap *ageom = new GeomMap();
      ageom->SetParam(par);
#if 0
      std::cout<<cid<<"  "<<seg;
      ageom->PrintMap();
#endif 
      key = KEY( cid, seg );
      geomContainer[key] = *ageom;
      delete ageom;
    }
    fin.close();
  }
  return true;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
GeomMap *GeomMapMan::GetMap( const int &cid, const int &seg)
{
  static const std::string funcname = "GeomMapMan::GetMap";
  unsigned int key;
  key = KEY( cid, seg );
  if(cid==CID_T0pre||cid==CID_T0post)
    key = KEY( CID_T0, seg );
  if(cid==CID_BHDpost)
    key = KEY( CID_BHD, seg );
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
bool GeomMapMan::GetParam( const int &cid, const int &seg, double *par)
{
  GeomMap *map=GetMap(cid,seg);
  if(!map) return false;
  map->GetParam(par);
  return true;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool GeomMapMan::GetPos( const int &cid, const int &seg, TVector3 &pos)
{
  GeomMap *map=GetMap(cid,seg);
  if(!map) return false;
  pos=map->GetPos();
  return true;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool GeomMapMan::GetPos( const int &cid, const int &seg, const double &lv, TVector3 &pos)
{
  GeomMap *map=GetMap(cid,seg);
  if(!map) return false;
  pos=map->GetPos(lv);
  return true;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool GeomMapMan::GetRot( const int &cid, const int &seg, TVector3 &rot)
{
  GeomMap *map=GetMap(cid,seg);
  if(!map) return false;
  rot=map->GetRot();
  return true;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool GeomMapMan::GetSize( const int &cid, const int &seg, double *size)
{
  GeomMap *map=GetMap(cid,seg);
  if(!map) return false;
  size=map->GetSize();
  return true;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool GeomMapMan::GetX( const int &cid, const int &seg, double &x)
{
  GeomMap *map=GetMap(cid,seg);
  if(!map) return false;
  double* tmp;
  tmp=map->GetSize();
  x = tmp[0];
  return true;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool GeomMapMan::GetY( const int &cid, const int &seg, double &y)
{
  GeomMap *map=GetMap(cid,seg);
  if(!map) return false;
  double* tmp;
  tmp=map->GetSize();
  y = tmp[2];
  return true;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool GeomMapMan::GetZ( const int &cid, const int &seg, double &z)
{
  GeomMap *map=GetMap(cid,seg);
  if(!map) return false;
  double* tmp;
  tmp=map->GetSize();
  z = tmp[2];
  return true;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool GeomMapMan::GetLightVelocity( const int &cid, const int &seg, double &lv )
{
  GeomMap *map=GetMap(cid,seg);
  if(!map) return false;
  lv=map->GetLightVelocity();
  return true;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool GeomMapMan::GetGPos( const int &cid, const int &seg, TVector3 &pos)
{
  GeomMap *gmap=GetMap(cid,0);
  GeomMap *map=GetMap(cid,seg);
  if(!map) return false;
  if(!gmap) return false;
  TVector3 gpos = gmap->GetPos();
  TVector3 grot = gmap->GetRot();
  pos = map->GetPos();
  pos.RotateX(grot.X()*Deg2Rad);
  pos.RotateY(grot.Y()*Deg2Rad);
  pos.RotateZ(grot.Z()*Deg2Rad);
  pos += gpos; 
  return true;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool GeomMapMan::GetGPos( const int &cid, const int &seg, const double &lv, TVector3 &pos)
{
  GeomMap *gmap=GetMap(cid,0);
  GeomMap *map=GetMap(cid,seg);
  if(!map) return false;
  if(!gmap) return false;
  TVector3 gpos = gmap->GetPos();
  TVector3 grot = gmap->GetRot();
  pos = map->GetPos(lv);
  pos.RotateX(grot.X()*Deg2Rad);
  pos.RotateY(grot.Y()*Deg2Rad);
  pos.RotateZ(grot.Z()*Deg2Rad);
  pos += gpos; 
  return true;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool GeomMapMan::GetXYCDS( const int &cid, const int &seg,TPolyLine &pline)
{
  if(cid!=CID_CDH&&cid!=CID_IH) return false;
  GeomMap *map=GetMap(cid,seg);
  if(!map) return false;
  TVector3 pos=map->GetPos();
  TVector3 rot=map->GetRot();
  double dy=map->GetThick();
  double dx=map->GetWidth();
  double x[5],y[5];
  TVector3 tmppos[5];
  tmppos[0]=TVector3(dx/2,dy/2,0);
  tmppos[1]=TVector3(-dx/2,dy/2,0);
  tmppos[2]=TVector3(-dx/2,-dy/2,0);
  tmppos[3]=TVector3(dx/2,-dy/2,0);
  tmppos[4]=TVector3(dx/2,dy/2,0);
  for(int i=0;i<5;i++){
    tmppos[i].RotateZ(rot.Z()*Deg2Rad);
    tmppos[i]+=pos;
    x[i]=tmppos[i].X();
    y[i]=tmppos[i].Y();
  }
  pline.DrawPolyLine(5,x,y);
  return true;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool GeomMapMan::GetZX( const int &cid, const int &seg,TPolyLine &pline, const bool &GLOBAL)
{
  if(cid==CID_CDH||cid==CID_IH) return false;
  GeomMap *map=GetMap(cid,seg);
  if(!map) return false;
  TVector3 pos=map->GetPos();
  TVector3 rot=map->GetRot();
  double dy=map->GetThick();
  double dx=map->GetWidth();
  double z[5],x[5];
  TVector3 tmppos[5];
  tmppos[0]=TVector3(dx/2,dy/2,0);
  tmppos[1]=TVector3(-dx/2,dy/2,0);
  tmppos[2]=TVector3(-dx/2,-dy/2,0);
  tmppos[3]=TVector3(dx/2,-dy/2,0);
  tmppos[4]=TVector3(dx/2,dy/2,0);
  for(int i=0;i<5;i++){
    tmppos[i].RotateY(rot.Y()*Deg2Rad);
    tmppos[i]+=pos;
    if(GLOBAL){
      GeomMap *gmap=GetMap(cid,0);
      if(gmap) return false;
      TVector3 gpos=gmap->GetPos();
      TVector3 grot=gmap->GetRot();
      tmppos[i].RotateY(grot.Y()*Deg2Rad);
      tmppos[i]+=gpos;
    }
    z[i]=tmppos[i].Z();
    x[i]=tmppos[i].X();
  }
  pline.DrawPolyLine(5,z,x);
  return false;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool GeomMapMan::SetParam( const int &cid, const int &seg,const double *par)
{
  GeomMap *map=GetMap(cid,seg);
  if(!map) return false;
  map->SetParam(par);
  return true;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool GeomMapMan::SetPos( const int &cid, const int &seg,const TVector3 &pos)
{
  GeomMap *map=GetMap(cid,seg);
  if(!map) return false;
  map->SetPos(pos);
  return true;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool GeomMapMan::SetRot( const int &cid, const int &seg,const TVector3 &rot)
{
  GeomMap *map=GetMap(cid,seg);
  if(!map) return false;
  map->SetRot(rot);
  return true;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool GeomMapMan::SetSize( const int &cid, const int &seg,const  double  *size)
{
  GeomMap *map=GetMap(cid,seg);
  if(!map) return false;
  map->SetSize(size);
  return true;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool GeomMapMan::SetLightVelocity( const int &cid, const int &seg, const double &lv )
{
  GeomMap *map=GetMap(cid,seg);
  if(!map) return false;
  map->SetLightVelocity(lv);
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void GeomMapMan::PrintMap( const int &Cid, std::ostream &p_out )
{
  static const std::string funcname = "GeomMapMan::PrintMap";
  unsigned int key;
  GeomMap ageom;
  int cid=0, cid_old=-1;
  int seg=0;
  std::cout << " ---- " << funcname << " ---- " << std::endl;
  for( GeomMapContainer::const_iterator i=geomContainer.begin();
       i!=geomContainer.end(); i++ ){
    key = i->first;
    ageom = i->second;
    cid = ((key>>CSHIFT)&CMASK);
    if( !( cid==Cid || Cid==-1 ) ) continue;
    if( cid!=cid_old ){
      p_out.setf(std::ios_base::fixed,std::ios_base::floatfield);
      p_out.setf(std::ios::showpoint);
      p_out<< "# CID  Seg   X         Y         Z        rotX       rotY       rotZ        "
	   << "SizeX(W)  SizeY(L)  SizeZ(T)  Size4  LightV"<<std::endl;
      p_out<< "#            [cm]      [cm]      [cm]     [deg]      [deg]      [deg]       "
	   << "[cm]      [cm]      [cm]             [cm/ns]"<<std::endl;

      p_out<< "# CID  Seg   R         Phi       Z        rotX       rotY       rotZ        "
	   << "SizeRmin  SizeRmax  SizedPhi  SizeL  LightV"<<std::endl;
      p_out<< "#            [cm]      [deg]     [cm]     [deg]      [deg]      [deg]       "
	   << "[cm]      [cm]      [deg]     [cm]   [cm/ns]"<<std::endl;

      cid_old = cid;
    }
    seg = ((key>>SSHIFT)&SMASK);;
    ageom = i->second;
    p_out << std::setw(6) << cid
	  << std::setw(5) << seg;
    ageom.PrintMap(p_out);
    if(seg==0) p_out<<"###"<<std::endl;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void GeomMapMan::PrintMapBL( std::ostream &p_out )
{
  static const std::string funcname = "GeomMapMan::PrintMapBL";
  PrintMap(CID_BHD,p_out);
  PrintMap(CID_T0,p_out);
  PrintMap(CID_E0,p_out);
  PrintMap(CID_LC1,p_out);
  PrintMap(CID_LC2,p_out);
  PrintMap(CID_AC,p_out);
  PrintMap(CID_WC,p_out);
  PrintMap(CID_GC,p_out);
  PrintMap(CID_BVC,p_out);
  PrintMap(CID_NC,p_out);
  PrintMap(CID_CVC,p_out);
  PrintMap(CID_LB,p_out);
  PrintMap(CID_PC,p_out);
  PrintMap(CID_BPD,p_out);
  PrintMap(CID_BD,p_out);
}
void GeomMapMan::PrintMapCDS( std::ostream &p_out )
{
  static const std::string funcname = "GeomMapMan::PrintMapCDS";
  PrintMap(CID_CDH,p_out);
  PrintMap(CID_IH,p_out);
}
