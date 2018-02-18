// GeomMapMan.h

#ifndef GeomMapMan_h
#define GeomMapMan_h 1

#include <map>
#include <string>
#include <iostream>
#include <stdlib.h>

#include "TObject.h"

class GeomMap : public TObject
{
 public:
  GeomMap();
  ~GeomMap() {};
 private:
  double GX,GY,GZ,dGX,dGY,dGZ; // position of origin in global world. but now, not decided. ('09/11/26)
  double X, Y, Z, dX, dY, dZ;
  double Length, Width, Thick;
  double LightVelocity;
 public:
  void SetParam( const double   &x, const double   &y, const double  &z,
		 const double  &dx, const double  &dy, const double &dz,
		 const double &len, const double &wid, const double &th )
  {X=x;Y=y;Z=z;dX=dx;dY=dy;dZ=dz;Length=len;Width=wid;Thick=th;}

  void SetParam( const double   &x, const double   &y, const double  &z,
		 const double  &dx, const double  &dy, const double &dz,
		 const double &len, const double &wid, const double &th,
		 const double &lv )
  {X=x;Y=y;Z=z;dX=dx;dY=dy;dZ=dz;Length=len;Width=wid;Thick=th;LightVelocity=lv;}
  
  void SetGParam( const double  &x, const double  &y, const double  &z, 
		  const double &dx, const double &dy, const double &dz )
  {GX=x;GY=y;GZ=z;dGX=dx;dGY=dy;dGZ=dz;}
  
  void SetX( const double &x ) { X = x; }
  void SetY( const double &y ) { Y = y; }
  void SetZ( const double &z ) { Z = z; }
  void SetdX( const double &x ) { dX = x; }
  void SetdY( const double &y ) { dY = y; }
  void SetdZ( const double &z ) { dZ = z; }
  void SetLength( const double &len ) { Length = len; }
  void SetWidth( const double &w ) { Width = w; }
  void SetThick( const double &t ) { Thick = t; }
  void SetLightVelocity( const double &lv ) { LightVelocity = lv; }

  double GetX() { return X; }
  double GetY() { return Y; }
  double GetZ() { return Z; }
  double GetdX() { return dX; }
  double GetdY() { return dY; }
  double GetdZ() { return dZ; }
  double GetLength() { return Length; }
  double GetWidth() { return Width; }
  double GetThick() { return Thick; }
  double GetLightVelocity() { return LightVelocity; }

  double GetGX() { return GX; }
  double GetGY() { return GY; }
  double GetGZ() { return GZ; }
  double GetdGX() { return dGX; }
  double GetdGY() { return dGY; }
  double GetdGZ() { return dGZ; }

  ClassDef( GeomMap, 1 );
};

class GeomMapMan : public TObject
{
 public:
  GeomMapMan();
  ~GeomMapMan();

  void SetFileNameCDS( const std::string & filenameCDS );
  void SetFileNameBL( const std::string & filenameBL );
  bool Initialize();

  GeomMapMan( const GeomMapMan &right );
 private:

  std::string FileNameCDS;
  std::string FileNameBL;

  typedef std::map <unsigned int, GeomMap> GeomMapContainer;
  GeomMapContainer geomContainer;

 public:
  bool GetLightVelocity( const int &cid, const int &seg,double &lv);
  bool GetParam( const int &cid, const int &seg,
		 double &x, double &y, double &z, double &dx, double &dy, double &dz,
		 double &len, double &wid, double &th );
  bool GetParam( const int &cid, const int &seg,
		 double &x, double &y, double &z, double &dx, double &dy, double &dz,
		 double &len, double &wid, double &th, double &lv );

  bool SetParam( const int &cid, const int &seg,
		 const double &x, const double &y, const double &z, const double &dx, const double &dy,const double &dz,
		 const double &len, const double &wid, const double &th, const double &lv );

  bool GetGParam( const int &cid, double &x, double &y, double &z, double &dx, double &dy, double &dz );

  bool SetLightVelocity( const int &cid, const int &seg,const double &lv);
  std::string GetFileNameCDS() { return FileNameCDS; }
  std::string GetFileNameBL() { return FileNameBL; }

  void PrintMap( const int &id,std::ostream &p_out = std::cout );
  void PrintMapBL(std::ostream &p_out = std::cout );

  ClassDef( GeomMapMan, 1 );
};

#endif
