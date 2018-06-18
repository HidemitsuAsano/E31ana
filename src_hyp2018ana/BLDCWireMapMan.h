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
#include "TVector3.h"
#include "GeomMapMan.h"

class BLDCWireMap : public TObject
{
 public:
  BLDCWireMap();
  ~BLDCWireMap() {};

 private:
  int nWire;
  double Z;
  int XY;
  double XY0, dXY, WireLength, TiltAngle, RotationAngle;
  double WirePhi;

 public:
  void SetParam( const int &nw, const double &z, const int &xy, const double &xy0, const double &dxy,
		 const double &wl, const double &tilt, const double &ra );
  void SetParam(double *param);

  void SetNWire( const int &nw ) { nWire = nw; }
  void SetZ(     const double &z ) { Z = z; }
  void SetXY( const int &xy ) { XY = xy; }
  void SetXY0( const double &xy0 ) { XY0 = xy0; }
  void SetdXY( const double &dxy ) { dXY = dxy; }
  void SetWireLength( const double &wl ) { WireLength = wl; }
  void SetTiltAngle( const double &tilt ) { TiltAngle = tilt; }
  void SetRotationAngle( const double &ra ) { RotationAngle = ra; }
  void SetWirePhi( const double &phi ) { WirePhi = phi; }

  int GetNWire() { return nWire; }
  double GetZ() { return Z; }
  int GetXY() { return XY; }
  double GetXY0() { return XY0; }
  double GetdXY() { return dXY; }
  double GetWireLength() { return WireLength; }
  double GetTiltAngle() { return TiltAngle; }
  double GetRotationAngle() { return RotationAngle; }
  double GetWirePhi() { return WirePhi; }

  ClassDef( BLDCWireMap, 1 );
};

class BLDCWireMapMan : public TObject
{
 public:
  BLDCWireMapMan();
  ~BLDCWireMapMan() {};

  void SetFileName( const std::string & filename ) { FileName = filename; }
  bool Initialize();

 private:
  std::string FileName;

  typedef std::map < unsigned int, BLDCWireMap > BLDCWireMapContainer;
  BLDCWireMapContainer bldcContainer;
  typedef std::map < unsigned int, GeomMap > GeomMapContainer;
  GeomMapContainer geomContainer;

  static const unsigned int KEYMASK = 0x000F;
  static const unsigned int CMASK   = 0x00FF;
  static const unsigned int LMASK   = 0x00FF;
  static const int          CSHIFT  = 4;
  static const int          LSHIFT  = 16;
  static const unsigned int KEYFLAG = 0x0003; 
  inline int KEY(const int &cid,const int &layer)
  {  return ((((cid)&CMASK)<<CSHIFT) | (((layer)&LMASK)<<LSHIFT) | KEYFLAG ); }
  static const int MAXCHAR = 256;
  
 public:
  GeomMap * GetGMap( const int &cid, const int &layer);
  BLDCWireMap * GetWireMap( const int &cid, const int &layer);

  bool GetParam( const int &cid, const int &layer,
		 int &nw, double &z, int &xy, double &xy0, double &dxy, 
		 double &wl, double &tilt, double &ra );
  bool GetXY0( const int &cid, const int &layer,double &xy0);
  bool SetXY0( const int &cid, const int &layer,const double &xy0);
  bool SetParam( const int &cid, const int &layer,
		 const int &nw,const double &z,const int &xy,const double &xy0,const double &dxy, 
		 const double &wl,const double &tilt,const double &ra );
  int  GetNWire( const int &cid, const int &layer );
  bool GetGParam( const int &cid, TVector3 &pos, TVector3 &rot );

  bool GetGParam( const int &cid, const int &seg, double *par);
  bool SetGParam( const int &cid, const int &seg, const double *par);

  std::string GetFileName() { return FileName; }

  void PrintMap( const int &id,std::ostream &p_out = std::cout );
  void PrintMapBL(std::ostream &p_out = std::cout );
  ClassDef( BLDCWireMapMan, 1 );
};
#endif
