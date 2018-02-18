// GeomMapMan.h

#ifndef GeomMapMan_h
#define GeomMapMan_h 1

#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>

#include "TObject.h"
#include "TVector3.h"
#include "TMath.h"
#include "TPolyLine.h"
#include "GlobalVariables.h"

class GeomMap : public TObject
{
 public:
  GeomMap();
  ~GeomMap() {};
 private:
  static const int npar=11;
  double R,Theta;
  double Pos[3];
  double Rot[3];
  double Size[4];
  double LightVelocity;
  
 public:
  inline void SetParam(const double *par);
  inline void GetParam(double *par);

  void SetPos( const TVector3 &pos ) { Pos[0]=pos.X(),Pos[1]=pos.Y(),Pos[2]=pos.Z(); }          
  void SetX( const double &x ) { Pos[0]=x; }
  void SetY( const double &y ) { Pos[1]=y; }
  void SetZ( const double &z ) { Pos[2]=z; }
  void SetRTheta( const double &r, const double &theta)
  { R=r; Theta=theta; Pos[0]=r*TMath::Cos(theta*Deg2Rad);  Pos[1]=r*TMath::Sin(theta*Deg2Rad); }
  void SetRot( const TVector3 &rot ) { Rot[0]=rot.X(),Rot[1]=rot.Y(),Rot[2]=rot.Z(); }          
  void SetRotX( const double &x ) { Rot[0]=x; }
  void SetRotY( const double &y ) { Rot[1]=y; }
  void SetRotZ( const double &z ) { Rot[2]=z; }
  void SetSize( const double *size ) { for( int i =0; i<4;i++ ) Size[i]=size[i]; }
  void SetLength( const double &len ) { Size[1] = len; }
  void SetWidth( const double &w ) { Size[0] = w; }
  void SetThick( const double &t ) { Size[2] = t; }
  void SetLightVelocity( const double &lv ) { LightVelocity = lv; }

  TVector3 GetPos() { return TVector3(Pos); }
  inline TVector3 GetPos(const double &ctsub);
  TVector3 GetRot() { return TVector3(Rot); }
  double* const GetSize() { return Size; }

  double GetX() { return Pos[0]; }
  double GetY() { return Pos[1]; }
  double GetZ() { return Pos[2]; }
  double GetR() { return R; }
  double GetTheta() { return Theta; }

  double GetRotX() { return Rot[0]; }
  double GetRotY() { return Rot[1]; }
  double GetRotZ() { return Rot[2]; }

  double GetLength() { return IsBox() ? Size[1] : Size[3]; }
  double GetWidth() { return IsBox() ? Size[0] : Size[0]*Size[2]*Deg2Rad; }
  double GetThick() { return IsBox() ? Size[2] : Size[1]-Size[0]; }

  double GetRmin() { return Size[0]; }
  double GetRmax() { return Size[1]; }
  double GetRmean() { return (Size[0]+Size[0])/2.; }
  double GetdPhi() { return Size[2]; }

  bool IsBox() { return Size[3]>0 ? false : true; }
  bool IsTube() { return Size[3]>0 ? true : false; }
  bool IsCartesian() { return Size[3]<0 ? false : true; }
  bool IsCyl() { return Size[3]<0 ? true : false; }

  double GetLightVelocity() { return LightVelocity; }
  void PrintMap( std::ostream &p_out = std::cout );

  ClassDef( GeomMap, 1 );
};

inline void GeomMap::SetParam( const double *param ){
  for(int i=0;i<3;i++)
    Pos[i]=param[i];
  for(int i=0;i<3;i++)
    Rot[i]=param[i+3];
  for(int i=0;i<4;i++)
    Size[i]=param[i+6];
  if(IsCyl()) SetRTheta(param[0],param[1]);
  if(IsTube()&&Size[2]<90){
    double tmpr=(GetRmin()+GetRmax())/2.; double tmptheta=GetRotZ();
    SetRTheta(tmpr,tmptheta);
    Pos[0]+=param[0];
    Pos[1]+=param[1];
    R=param[0];
    Theta=param[1];
  }
  LightVelocity=param[10];
}

inline void GeomMap::GetParam( double *param ){
  for(int i=0;i<3;i++)
    param[i]=Pos[i];
  for(int i=0;i<3;i++)
    param[i+3]=Rot[i];
  for(int i=0;i<4;i++)
    param[i+6]=Size[i];
  if(IsCyl()||(IsTube()&&Size[2]<90)){ param[0]=R; param[1]=Theta; }
  param[10]=LightVelocity;
}

inline TVector3 GeomMap::GetPos(const double &ctsub){
  TVector3 dis(0,ctsub*LightVelocity,0);
  if(IsTube()) dis.SetXYZ(0,0,ctsub*LightVelocity);
  dis.RotateX(Rot[0]);
  dis.RotateY(Rot[1]);
  dis.RotateZ(Rot[2]);
  return TVector3(Pos)+dis;
}

class GeomMapMan : public TObject
{
 public:
  GeomMapMan();
  ~GeomMapMan();

  void SetFileNameCDS( const std::string & filenameCDS );
  void SetFileNameBL( const std::string & filenameBL );
  void SetFileNameHall( const std::string & filenameHall );
  bool Initialize();

 private:
  static const int npar=11;
  std::string FileNameCDS;
  std::string FileNameBL;
  std::string FileNameHall;
  
  typedef std::map <unsigned int, GeomMap> GeomMapContainer;
  GeomMapContainer geomContainer;

  static const unsigned int KEYMASK  = 0x000F;
  static const unsigned int SMASK    = 0x00FF;
  static const unsigned int CMASK    = 0x00FF;
  static const int          SSHIFT   = 4;
  static const int          CSHIFT   = 12;
  static const unsigned int KEYFLAG  = 0x0003;
  inline int KEY(const int &cid,const int &seg)	
  { return ((((cid)&CMASK)<<CSHIFT) | (((seg)&SMASK)<<SSHIFT) | KEYFLAG ); }
  bool ReadFile(const std::string filename);
  GeomMap *GetMap(const int &cid, const int &seg);  

 public:
  bool GetParam( const int &cid, const int &seg, double *par);
  bool GetPos( const int &cid, const int &seg, TVector3 &pos);
  bool GetPos( const int &cid, const int &seg, const double &lv, TVector3 &pos);
  bool GetRot( const int &cid, const int &seg, TVector3 &rot);
  bool GetSize( const int &cid, const int &seg, double *par );
  bool GetLightVelocity( const int &cid, const int &seg, double &lv);  
  bool GetGPos( const int &cid, const int &seg, TVector3 &pos);
  bool GetGPos( const int &cid, const int &seg, const double &lv, TVector3 &pos);
  //  bool GetGRot( const int &cid, const int &seg, TVector3 &rot);

  bool GetXYCDS( const int &cid, const int &seg, TPolyLine &pline);
  //  bool GetYZSlice( const int &cid, const int &seg, TPolyLine &pline);
  bool GetZX( const int &cid, const int &seg, TPolyLine &pline, const bool &GLOBAL);

  bool SetParam( const int &cid, const int &seg, const double *par);
  bool SetPos( const int &cid, const int &seg, const TVector3 &pos);
  bool SetRot( const int &cid, const int &seg, const TVector3 &rot);
  bool SetSize( const int &cid, const int &seg, const double *par);
  bool SetLightVelocity( const int &cid, const int &seg,const double &lv);

  std::string GetFileNameCDS() { return FileNameCDS; }
  std::string GetFileNameBL() { return FileNameBL; }

  void PrintMap( const int &id,std::ostream &p_out = std::cout );
  void PrintMapBL(std::ostream &p_out = std::cout );
  void PrintMapCDS(std::ostream &p_out = std::cout );

  ClassDef( GeomMapMan, 1 );
};

#endif
