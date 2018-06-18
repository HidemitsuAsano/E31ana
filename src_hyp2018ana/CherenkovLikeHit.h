#ifndef CherenkovLikeHit_h
#define CherenkovLikeHit_h 1

#include <iostream>

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#include "ConfMan.h"
#include "ADCHit.h"
#include "TDCHit.h"

class CherenkovLikeHit : public TObject
{
 private:
  int CounterID;
  int HitID;
  int NSensor; 
  int Seg;
  
 public:
  CherenkovLikeHit();
  //  CherenkovLikeHit( const CherenkovLikeHit &hit );
  virtual ~CherenkovLikeHit() {};

 private:
  // raw data
  ADCHit adchit[4];
  TDCHit tdchit[4];

 public:
  int cid()  const { return CounterID; }
  int hid()  const { return HitID; }
  int seg()  const { return Seg; }
  int adc(const int &i) const { return (0<i&&i<=4 ) ? adchit[i-1].data():0; }
  int tdc(const int &i) const { return (0<i&&i<=4 ) ? tdchit[i-1].data():0; }
  int cra(const int &i) const { return (0<i&&i<=4 ) ? adchit[i-1].cr():0; }
  int sla(const int &i) const { return (0<i&&i<=4 ) ? adchit[i-1].sl():0; }
  int cha(const int &i) const { return (0<i&&i<=4 ) ? adchit[i-1].ch():0; }
  int crt(const int &i) const { return (0<i&&i<=4 ) ? tdchit[i-1].cr():0; }
  int slt(const int &i) const { return (0<i&&i<=4 ) ? tdchit[i-1].sl():0; }
  int cht(const int &i) const { return (0<i&&i<=4 ) ? tdchit[i-1].ch():0; }

  int val(const int &at, const int &i) { return at ? tdc(i) : adc(i); }
  int cr(const int &at, const int &i) { return at ? crt(i) : cra(i); }
  int sl(const int &at, const int &i) { return at ? slt(i) : sla(i); }
  int ch(const int &at, const int &i) { return at ? cht(i) : cha(i); }

  void SetCounterID( const int & cid ) { CounterID = cid; }
  void SetHitID( const int & hid ) { HitID = hid; }
  void SetSegment( const int & seg ) { Seg = seg; }
  void SetADC(const int &i, const int & a ) { adchit[i].SetData(a); }
  void SetTDC(const int &i, const int & t ) { tdchit[i].SetData(t); }
  void SetCrateA(const int &i, const int & cr )   { adchit[i].SetCrate(cr); }
  void SetSlotA(const int &i, const int & sl )    { adchit[i].SetSlot(sl); }
  void SetChannelA(const int &i, const int & ch ) { adchit[i].SetChannel(ch); }
  void SetCrateT(const int &i, const int & cr )   { tdchit[i].SetCrate(cr); }
  void SetSlotT(const int &i, const int & sl )    { tdchit[i].SetSlot(sl); }
  void SetChannelT(const int &i, const int & ch ) { tdchit[i].SetChannel(ch); }
  
  bool SetData( const int &c, const int &n, const int &a,
		const int &cid, const int &seg, const int &at, const int &ud,
		const int &data );

  bool CheckRange( int num_required = 2 );

 private:
  /* double  X,  Y,  Z; // position */
  /* double dX, dY, dZ; // direction */
  /* double  GX,  GY,  GZ; // origin position */
  /* double dGX, dGY, dGZ; // origin direction */
  double Length, Width, Thick;
  //  TVector3 HitPos;

 public:
  double np(const int &i) const { return (0<i&&i<=4 ) ? adchit[i-1].energy():0; }
  double time(const int &i) const { return (0<i&&i<=4 ) ? tdchit[i-1].time():0; }
  /* double x()   const { return X; } */
  /* double y()   const { return Y; } */
  /* double z()   const { return Z; } */
  /* double dx()  const { return dX; } */
  /* double dy()  const { return dY; } */
  /* double dz()  const { return dZ; } */
  /* double gx()  const { return GX; } */
  /* double gy()  const { return GY; } */
  /* double gz()  const { return GZ; } */
  /* double dgx() const { return dGX; } */
  /* double dgy() const { return dGY; } */
  /* double dgz() const { return dGZ; } */
  double len() const { return Length; }
  double wid() const { return Width; }
  double thick() const { return Thick; }
  /* void   pos( double &x, double &y, double &z ) { x=X;  y=Y;  z=Z; } */
  /* void   dir( double &x, double &y, double &z ) { x=dX; y=dY; z=dZ; } */
  /* void   gpos( double &x, double &y, double &z ) { x=GX;  y=GY;  z=GZ; } */
  /* void   gdir( double &x, double &y, double &z ) { x=dGX; y=dGY; z=dGZ; } */

  void SetNumPhoton( const int &i, const double &np ) { adchit[i].SetEnergy(np); }
  /* void SetTime(      const int &i, const double &t  ) { tdchit[i].SetTime(t); } */
  /* void SetPos( const double &x, const double &y, const double &z ) { X=x;  Y=y;  Z=z; } */
  /* void SetDir( const double &x, const double &y, const double &z ) { dX=x; dY=y; dZ=z; } */
  /* void SetGPos( const double &x, const double &y, const double &z ) { GX=x;  GY=y;  GZ=z; } */
  /* void SetGDir( const double &x, const double &y, const double &z ) { dGX=x; dGY=y; dGZ=z; } */
  void SetLength( const double &l ) { Length = l; }
  void SetWidth(  const double &w ) { Width  = w; }
  void SetThick(  const double &t ) { Thick  = t; }

  bool Calc( ConfMan *conf );
  void Clear();

  ClassDef(CherenkovLikeHit, 1 );
};

#endif
