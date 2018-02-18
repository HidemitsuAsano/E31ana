#ifndef CherenkovLikeHit_h
#define CherenkovLikeHit_h 1

#include <iostream>

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#include "ConfMan.h"

class CherenkovLikeHit : public TObject
{
 private:
  int CounterID;
  int HitID; 
  
 public:
  CherenkovLikeHit();
  CherenkovLikeHit( const CherenkovLikeHit &hit );
  virtual ~CherenkovLikeHit() {};

 private:
  // raw data
  int Seg;
  int ADC[4];
  int TDC[4];
  int CrateA[4],SlotA[4],ChannelA[4];  
  int CrateT[4],SlotT[4],ChannelT[4];  

 public:
  int cid()  const { return CounterID; }
  int hid()  const { return HitID; }
  int seg()  const { return Seg; }
  int adc(const int &i) const { return (0<i&&i<=4 ) ? ADC[i-1]:0; }
  int tdc(const int &i) const { return (0<i&&i<=4 ) ? TDC[i-1]:0; }
  int cra(const int &i) const { return (0<i&&i<=4 ) ? CrateA[i-1]:0; }
  int sla(const int &i) const { return (0<i&&i<=4 ) ? SlotA[i-1]:0; }
  int cha(const int &i) const { return (0<i&&i<=4 ) ? ChannelA[i-1]:0; }
  int crt(const int &i) const { return (0<i&&i<=4 ) ? CrateT[i-1]:0; }
  int slt(const int &i) const { return (0<i&&i<=4 ) ? SlotT[i-1]:0; }
  int cht(const int &i) const { return (0<i&&i<=4 ) ? ChannelT[i-1]:0; }

  void SetCounterID( const int & cid ) { CounterID = cid; }
  void SetHitID( const int & hid ) { HitID = hid; }
  void SetSegment( const int & seg ) { Seg = seg; }
  void SetADC(const int &i, const int & a ) { ADC[i] = a; }
  void SetTDC(const int &i, const int & t ) { TDC[i] = t; }
  void SetCrateA(const int &i, const int & cr )   { CrateA[i]   = cr; }
  void SetSlotA(const int &i, const int & sl )    { SlotA[i]    = sl; }
  void SetChannelA(const int &i, const int & ch ) { ChannelA[i] = ch; }
  void SetCrateT(const int &i, const int & cr )   { CrateT[i]   = cr; }
  void SetSlotT(const int &i, const int & sl )    { SlotT[i]    = sl; }
  void SetChannelT(const int &i, const int & ch ) { ChannelT[i] = ch; }
  
  bool SetData( const int &c, const int &n, const int &a,
		const int &cid, const int &seg, const int &at, const int &ud,
		const int &data );

  bool CheckRange( int num_required = 2 );

 private:
  // conversion data
  double NumPhoton[4];
  double Time[4];

  double  X,  Y,  Z; // position
  double dX, dY, dZ; // direction
  double  GX,  GY,  GZ; // origin position
  double dGX, dGY, dGZ; // origin direction
  double Length, Width, Thick;

 public:
  double np(const int &i) const { return (0<i&&i<=4 ) ? NumPhoton[i-1]:0; }
  double time(const int &i) const { return (0<i&&i<=4 ) ? Time[i-1]:0; }
  double x()   const { return X; }
  double y()   const { return Y; }
  double z()   const { return Z; }
  double dx()  const { return dX; }
  double dy()  const { return dY; }
  double dz()  const { return dZ; }
  double gx()  const { return GX; }
  double gy()  const { return GY; }
  double gz()  const { return GZ; }
  double dgx() const { return dGX; }
  double dgy() const { return dGY; }
  double dgz() const { return dGZ; }
  double len() const { return Length; }
  double wid() const { return Width; }
  double thick() const { return Thick; }
  void   pos( double &x, double &y, double &z ) { x=X;  y=Y;  z=Z; }
  void   dir( double &x, double &y, double &z ) { x=dX; y=dY; z=dZ; }
  void   gpos( double &x, double &y, double &z ) { x=GX;  y=GY;  z=GZ; }
  void   gdir( double &x, double &y, double &z ) { x=dGX; y=dGY; z=dGZ; }

  void SetNumPhoton( const int &i, const double &np ) { NumPhoton[i] = np; }
  void SetTime(      const int &i, const double &t  ) { Time[i] = t; }
  void SetPos( const double &x, const double &y, const double &z ) { X=x;  Y=y;  Z=z; }
  void SetDir( const double &x, const double &y, const double &z ) { dX=x; dY=y; dZ=z; }
  void SetGPos( const double &x, const double &y, const double &z ) { GX=x;  GY=y;  GZ=z; }
  void SetGDir( const double &x, const double &y, const double &z ) { dGX=x; dGY=y; dGZ=z; }
  void SetLength( const double &l ) { Length = l; }
  void SetWidth(  const double &w ) { Width  = w; }
  void SetThick(  const double &t ) { Thick  = t; }

  bool Calc( ConfMan *conf );
  void Clear();

  ClassDef(CherenkovLikeHit, 1 );
};

#endif
