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

class ChamberLikeHit : public TObject
{
 private:
  int CounterID;
  int HitID;
  int XY; 			/* x: 0, y: 1 */

 public:
  ChamberLikeHit();
  ChamberLikeHit( const ChamberLikeHit &right );
  virtual ~ChamberLikeHit() {};

 private:
  int Layer, Wire;
  int TDC0;

  int Crate, Slot, Channel;

 public:
  int cid()   const { return CounterID; }
  int hid()   const { return HitID; }
  int layer() const { return Layer; }
  int wire()  const { return Wire; }
  int tdc()   const { return TDC0; }
  int cr()    const { return Crate; }
  int sl()    const { return Slot; }
  int ch()    const { return Channel; }
  int xy()    const { return XY; }

  void SetCounterID( const int & cid ) { CounterID = cid; }
  void SetHitID(     const int & hid ) { HitID = hid; }
  void SetLayer(     const int & lay ) { Layer = lay; }
  void SetWire(      const int &wire ) { Wire = wire; }
  void SetTDC(       const int & tdc ) { TDC0 = tdc; }
  void SetCrate(     const int & cr  ) { Crate = cr; }
  void SetSlot(      const int & sl  ) { Slot = sl; }
  void SetChannel(   const int & ch  ) { Channel = ch; }
  void SetXY(        const int & xy  ) { XY = xy; }

  bool SetData( const int &c, const int &n, const int &a,
		const int &cid, const int &layer, const int &wire, const int &data );
  
 private:
  double DriftTime;
  double DriftLength;
  double TimeOffset;
  double CorrDriftLength;

  double WirePosX, WirePosY, WirePosZ; // wire position at readout
  double WireLength;
  double TiltAngle;
  double Resolution;

  double HitPosX, HitPosY, HitPosZ;
  double dXY;
  int LR; 			/* left: 0    right: 1 */

  double RotationAngle;

  double  GX,  GY,  GZ; // origin position
  double dGX, dGY, dGZ; // origin direction

 public:
  
  double dt()    const { return DriftTime; }
  double dl()    const { return DriftLength; }
  double toffs() const { return TimeOffset; }
  double cdl()   const { return CorrDriftLength; }

  double wpos(const int &xy) const { return xy ? WirePosY : WirePosX; }
  double wx()    const { return WirePosX; }
  double wy()    const { return WirePosY; }
  double wz()    const { return WirePosZ; }
  double length()const { return WireLength; }
  double tilt()  const { return TiltAngle; }
  double resl()  const { return Resolution; }
  double hitpos(const int &xy) const { return xy ? HitPosY : HitPosX; }
  double x()     const { return HitPosX; }
  double y()     const { return HitPosY; }
  double z()     const { return HitPosZ; }
  double rot()   const { return RotationAngle; }
  double gx()    const { return GX; }
  double gy()    const { return GY; }
  double gz()    const { return GZ; }
  double dgx()   const { return dGX; }
  double dgy()   const { return dGY; }
  double dgz()   const { return dGZ; }
  int leftright() const { return LR; }

  double dxy() const { return dXY; }
  void gwpos(double &x,double &theta, double &z,const bool &TILT=false,const bool &GPOS=false);
  void gwpos2(double &x,double &z,const double &slope, const double &intercept,const bool &TILT=false,const bool &YPOS=false);

  void SetDriftTime(   const double &t ) { DriftTime = t; }
  void SetDriftLength( const double &l ) { DriftLength = l; }
  void SetTimeOffset(  const double &t ) { TimeOffset = t; }
  void SetCorrDriftLength( const double &cl ) { CorrDriftLength = cl; }

  void SetWirePosition( const double &x, const double &y, const double &z )
    { WirePosX = x; WirePosY = y; WirePosZ = z; }
  void SetWirePosX(   const double &x ) { WirePosX = x; }
  void SetWirePosY(   const double &y ) { WirePosY = y; }
  void SetWirePosZ(   const double &z ) { WirePosZ = z; }
  void SetWireLength( const double &l ) { WireLength = l; }
  void SetGPos( const double &x, const double &y, const double &z ) {  GX=x;  GY=y;  GZ=z; }
  void SetGDir( const double &x, const double &y, const double &z ) { dGX=x; dGY=y; dGZ=z; }
  void SetTiltAngle(     const double &a ) { TiltAngle = a; }
  void SetRotationAngle( const double &a ) { RotationAngle = a; }
  void SetResolution(    const double &r ) { Resolution = r; }
  void SetHitPosition(   const double &x, const double &y, const double &z )
    { HitPosX = x; HitPosY = y; HitPosZ = z; }
  void SetHitPos( const int &xy, const double &x ) { 
    if(xy) HitPosY = x; else HitPosX = x; }
  void SetHitPosX( const double &x ) { HitPosX = x; }
  void SetHitPosY( const double &y ) { HitPosY = y; }
  void SetHitPosZ( const double &z ) { HitPosZ = z; }

  void SetLeftRight( const int &lr ) { LR = lr; }
  
  //bool Calc( ConfMan *conf );
  bool Calc( ConfMan *conf, const double &toffs=0 );
  bool SetTimeOffset( ConfMan *conf, const double &toffs=0 );
  bool CheckRange(const double &low=0.001, const double &high=0.999);

  void Clear();

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
