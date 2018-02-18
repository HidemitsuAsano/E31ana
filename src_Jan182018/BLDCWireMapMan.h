// BLDCWireMapMan.h

#ifndef BLDCWireMapMan_h
#define BLDCWireMapMan_h 1

#include <map>
#include <string>
#include <iostream>
#include <stdlib.h>

#ifndef ROOT_TObject
#include "TObject.h"
#endif

class BLDCWireMap : public TObject
{
 public:
  BLDCWireMap();
  ~BLDCWireMap() {};

 private:
  double GX,GY,GZ,dGX,dGY,dGZ;
  int nWire;
  double Z;
  int XY;
  double XY0, dXY, WireLength, TiltAngle, RotationAngle;

 public:
  void SetParam( const int &nw, const double &z, const int &xy, const double &xy0, const double &dxy,
		 const double &wl, const double &tilt, const double &ra );
  void SetGParam( const double  &x, const double  &y, const double  &z, 
		  const double &dx, const double &dy, const double &dz );
  void SetNWire( const int &nw ) { nWire = nw; }
  void SetZ(     const double &z ) { Z = z; }
  void SetXY( const int &xy ) { XY = xy; }
  void SetXY0( const double &xy0 ) { XY0 = xy0; }
  void SetdXY( const double &dxy ) { dXY = dxy; }
  void SetWireLength( const double &wl ) { WireLength = wl; }
  void SetTiltAngle( const double &tilt ) { TiltAngle = tilt; }
  void SetRotationAngle( const double &ra ) { RotationAngle = ra; }

  int GetNWire() { return nWire; }
  double GetZ() { return Z; }
  int GetXY() { return XY; }
  double GetXY0() { return XY0; }
  double GetdXY() { return dXY; }
  double GetWireLength() { return WireLength; }
  double GetTiltAngle() { return TiltAngle; }
  double GetRotationAngle() { return RotationAngle; }

  double GetGX() { return GX; }
  double GetGY() { return GY; }
  double GetGZ() { return GZ; }
  double GetdGX() { return dGX; }
  double GetdGY() { return dGY; }
  double GetdGZ() { return dGZ; }

  ClassDef( BLDCWireMap, 1 );
};

class BLDCWireMapMan : public TObject
{
 public:
  BLDCWireMapMan();
  ~BLDCWireMapMan() {};

  void SetFileName( const std::string & filename ) { FileName = filename; }
  bool Initialize();

  BLDCWireMapMan( const BLDCWireMapMan &right );

 private:
  std::string FileName;

  typedef std::map < unsigned int, BLDCWireMap > BLDCWireMapContainer;
  BLDCWireMapContainer bldcContainer;

 public:
  bool GetParam( const int &cid, const int &layer,
		 int &nw, double &z, int &xy, double &xy0, double &dxy, 
		 double &wl, double &tilt, double &ra );
  bool GetXY0( const int &cid, const int &layer,double &xy0);
  bool SetXY0( const int &cid, const int &layer,const double &xy0);
  bool SetParam( const int &cid, const int &layer,
		 const int &nw,const double &z,const int &xy,const double &xy0,const double &dxy, 
		 const double &wl,const double &tilt,const double &ra );
  int GetNWire( const int &cid, const int &layer );
  bool GetGParam( const int &cid, double &x, double &y, double &z, double &dx, double &dy, double &dz );

  std::string GetFileName() { return FileName; }

  void PrintMap( const int &id,std::ostream &p_out = std::cout );
  void PrintMapBL(std::ostream &p_out = std::cout );
  ClassDef( BLDCWireMapMan, 1 );
};

#endif
