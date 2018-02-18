#ifndef ChamberLikeHit_h
#define ChamberLikeHit_h 1

#include <string>
#include <iostream>
#include <vector>

#ifndef ROOT_TObject
#include "TObject.h"
#endif
#include "TMath.h"

#include "ConfMan.h"
#include "TDCHit.h"
#include "GlobalVariables.h"

class ChamberLikeHit : public TObject
{
 protected:
  TDCHit tdchit;
  int CounterID;
  int HitID;
  int Layer,Wire;
  int XY; 			/* x: 0, y: 1 */
  int STATUS;

 public:
  ChamberLikeHit();
  //  ChamberLikeHit( const ChamberLikeHit &right );
  virtual ~ChamberLikeHit() {};

  int status() const { return STATUS; }
  void SetStatus(const int &sta) { STATUS = sta; }

 public:
  int cid()   const { return CounterID; }
  int hid()   const { return HitID; }
  int layer() const { return Layer; }
  int wire()  const { return Wire; }
  int tdc()   const { return tdchit.data(); }
  int cr()    const { return tdchit.cr(); }
  int sl()    const { return tdchit.sl(); }
  int ch()    const { return tdchit.ch(); }
  int xy()    const { return XY; }

  void SetCounterID( const int & cid ) { CounterID = cid; }
  void SetHitID(     const int & hid ) { HitID = hid; }
  void SetLayer(     const int & lay ) { Layer = lay; }
  void SetWire(      const int &wire ) { Wire = wire; }
  void SetTDC(       const int & tdc ) { tdchit.SetData(tdc); }
  void SetCrate(     const int & cr  ) { tdchit.SetCrate(cr); }
  void SetSlot(      const int & sl  ) { tdchit.SetSlot(sl); }
  void SetChannel(   const int & ch  ) { tdchit.SetChannel(ch); }
  void SetXY(        const int & xy  ) { XY = xy; }

  bool SetData( const int &c, const int &n, const int &a,
		const int &cid, const int &layer, const int &wire, const int &data );
  
 private:
  double DriftLength;
  double TimeOffset;
  double CorrDriftLength;

  TVector3 WirePos; // wire position at readout
  double WireLength;
  double TiltAngle;
  double Resolution;

  TVector3 HitPos;
  double dXY;
  int LR; 			/* left: 0    right: 1 */

  double RotationAngle;

  TVector3 GPos,GDir;

 public:
  
  double dt()    const { return tdchit.time(); }
  double cdt()   const { return dt()-toffs(); }
  double dl()    const { return DriftLength; }
  double toffs() const { return TimeOffset; }
  double cdl()   const { return CorrDriftLength; }

  TVector3 wpos() const { return WirePos; }
  TVector3 hitpos() const { return HitPos; }
  TVector3 gpos() const { return GPos; }
  TVector3 gdir() const { return GDir; }

  double wpos(const int &xy) const { return xy ? WirePos.Y() : WirePos.X(); }
  double wx()    const { return WirePos.X(); }
  double wy()    const { return WirePos.Y(); }
  double wz()    const { return WirePos.Z(); }
  double length()const { return WireLength; }
  double tilt()  const { return TiltAngle; }
  double resl()  const { return Resolution; }
  double hitpos(const int &xy) const { return xy ? HitPos.Y() : HitPos.X(); }
  double x()     const { return HitPos.X(); }
  double y()     const { return HitPos.Y(); }
  double z()     const { return HitPos.Z(); }
  double rot()   const { return RotationAngle; }
  double gx()    const { return GPos.X(); }
  double gy()    const { return GPos.Y(); }
  double gz()    const { return GPos.Z(); }
  double dgx()   const { return GDir.X(); }
  double dgy()   const { return GDir.Y(); }
  double dgz()   const { return GDir.Z(); }
  int leftright() const { return LR; }

  double dxy() const { return dXY; }
  void gwpos(double &x,double &theta, double &z,const bool &TILT=false,const bool &GPOS=false);

  void SetDriftTime(   const double &t ) { tdchit.SetTime(t); }
  void SetDriftLength( const double &l ) { DriftLength = l; }
  void SetTimeOffset(  const double &t ) { TimeOffset = t; }
  void SetCorrDriftLength( const double &cl ) { CorrDriftLength = cl; }

  void SetWirePosition( const double &x, const double &y, const double &z )
  { WirePos.SetXYZ(x,y,z); }
  void SetWirePosX(   const double &x ) { WirePos.SetX(x); }
  void SetWirePosY(   const double &y ) { WirePos.SetY(y); }
  void SetWirePosZ(   const double &z ) { WirePos.SetZ(z); }
  void SetWireLength( const double &l ) { WireLength = l; }
  void SetGPos( const double &x, const double &y, const double &z ) { GPos.SetXYZ(x,y,z); }
  void SetGDir( const double &x, const double &y, const double &z ) { GDir.SetXYZ(x,y,z); }
  void SetTiltAngle(     const double &a ) { TiltAngle = a; }
  void SetRotationAngle( const double &a ) { RotationAngle = a; }
  void SetResolution(    const double &r ) { Resolution = r; }
  void SetHitPosition(   const double &x, const double &y, const double &z )
  { HitPos.SetXYZ(x,y,z); }
  void SetHitPosition(   TVector3 &pos )
  { HitPos=pos; }
  void SetHitPos( const int &xy, const double &x ) { 
    if(xy) HitPos.SetY(x); else HitPos.SetX(x); }
  void SetHitPosX( const double &x ) { HitPos.SetX(x); }
  void SetHitPosY( const double &y ) { HitPos.SetY(y); }
  void SetHitPosZ( const double &z ) { HitPos.SetZ(z); }

  void SetLeftRight( const int &lr ) { LR = lr; }
  
  virtual bool Calc( ConfMan *conf, const double &toffs=0 );
  bool SetTimeOffset( ConfMan *conf, const double &toffs=0 );
  bool CheckRange(const double &low=0.001, const double &high=0.999);

  virtual void Clear();

  // simulation
 private:
  double ResolutionTrue; // Resolution for simulation;
  double TimeOffsetTrue; // Time offset, simulation knows true value.
 public:
  double resltrue() { return ResolutionTrue; }
  double toffstrue() const { return TimeOffsetTrue; }
  bool SetSimulatedResolution( ConfMan *conf ); //Resolution for simulation;
  void SetTimeOffsetTrue(  const double &t ) { TimeOffsetTrue = t; }
  void Reverse( ConfMan *conf );

  ClassDef( ChamberLikeHit, 1 );
};

#endif
