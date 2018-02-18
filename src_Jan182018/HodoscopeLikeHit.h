#ifndef HodoscopeLikeHit_h
#define HodoscopeLikeHit_h 1

#include <string>
#include <iostream>
#include <cmath>

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#include "ConfMan.h"

class HodoscopeLikeHit : public TObject
{
 private:
  int CounterID;
  int HitID; 
  
 public:
  HodoscopeLikeHit();
  HodoscopeLikeHit( const HodoscopeLikeHit &hit );
  virtual ~HodoscopeLikeHit() {};

 private:
  // raw data
  int Seg;
  int ADCu, ADCd, TDCu, TDCd; // UpDown ( or LeftRight )
  int CrateAu, SlotAu, ChannelAu;
  int CrateAd, SlotAd, ChannelAd;
  int CrateTu, SlotTu, ChannelTu;
  int CrateTd, SlotTd, ChannelTd;

 public:
  int cid()  const { return CounterID; }
  int hid()  const { return HitID; }
  int seg()  const { return Seg; }
  int adcu() const { return ADCu; }
  int adcd() const { return ADCd; }
  int tdcu() const { return TDCu; }
  int tdcd() const { return TDCd; }
  double tdcmean() const { if(TDCd>0&&TDCu>0) return (TDCd+TDCu)/2.; else return -1; }
  int crtu() const { return CrateTu; }
  int sltu() const { return SlotTu; }
  int chtu() const { return ChannelTu; }
  int crtd() const { return CrateTd; }
  int sltd() const { return SlotTd; }
  int chtd() const { return ChannelTd; }
  int crau() const { return CrateAu; }
  int slau() const { return SlotAu; }
  int chau() const { return ChannelAu; }
  int crad() const { return CrateAd; }
  int slad() const { return SlotAd; }
  int chad() const { return ChannelAd; }

  void SetCounterID( const int & cid ) { CounterID = cid; }
  void SetHitID(     const int & hid ) { HitID = hid; }
  void SetSegment(   const int & seg ) { Seg = seg; }
  void SetADCu(      const int & a   ) { ADCu = a; }
  void SetADCd(      const int & a   ) { ADCd = a; }
  void SetTDCu(      const int & t   ) { TDCu = t; }
  void SetTDCd(      const int & t   ) { TDCd = t; }
  void SetCrateTu(   const int & cr  ) { CrateTu = cr; }
  void SetSlotTu(    const int & sl  ) { SlotTu = sl; }
  void SetChannelTu( const int & ch  ) { ChannelTu = ch; }
  void SetCrateTd(   const int & cr  ) { CrateTd = cr; }
  void SetSlotTd(    const int & sl  ) { SlotTd = sl; }
  void SetChannelTd( const int & ch  ) { ChannelTd = ch; }
  void SetCrateAu(   const int & cr  ) { CrateAu = cr; }
  void SetSlotAu(    const int & sl  ) { SlotAu = sl; }
  void SetChannelAu( const int & ch  ) { ChannelAu = ch; }
  void SetCrateAd(   const int & cr  ) { CrateAd = cr; }
  void SetSlotAd(    const int & sl  ) { SlotAd = sl; }
  void SetChannelAd( const int & ch  ) { ChannelAd = ch; }
  
  bool SetData( const int &c, const int &n, const int &a,
		const int &cid, const int &seg, const int &at, const int &ud,
		const int &data );

  bool CheckRange(const int &type=0);
  bool CheckRange(const int &ll,const int &ul, const bool &type=false);
  bool CheckRange2(ConfMan *conf=0);

 private:
  // conversion data
  double Eneu, Ened, Timeu, Timed;
  double CTimeu, CTimed; // corrected time

  double TMean, CTMean;
  double TSub, CTSub;
  double EMean;

  double  X,  Y,  Z; // position
  double dX, dY, dZ; // direction
  double  GX,  GY,  GZ; // origin position
  double dGX, dGY, dGZ; // origin direction
  double Length, Width, Thick;

  double LightVelocity;
  double HitPosition; // hit position from the center along scinti.

  bool CHECKRANGE;

 public:
  double eu()  const { return   Eneu; }
  double ed()  const { return   Ened; }
  double emean()const { return EMean; }
  double tu()  const { return  Timeu; }
  double td()  const { return  Timed; }
  double ctu() const { return CTimeu; }
  double ctd() const { return CTimed; }
  double tmean() const { return  TMean; }
  double ctmean()const { return CTMean; }
  double tsub()  const { return TSub; }
  double ctsub() const { return CTSub; }
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
  double lv() const { return LightVelocity; }
  double hitpos() const { return HitPosition; }
  void   pos( double &x, double &y, double &z ) { x=X;  y=Y;  z=Z; }
  void   dir( double &x, double &y, double &z ) { x=dX; y=dY; z=dZ; }
  void   gpos( double &x, double &y, double &z ) { x=GX;  y=GY;  z=GZ; }
  void   gdir( double &x, double &y, double &z ) { x=dGX; y=dGY; z=dGZ; }

  void SetEneu(   const double &e ) { Eneu = e; }
  void SetEned(   const double &e ) { Ened = e; }
  void SetEMean(  const double &e ) { EMean = e; }
  void SetTimeu(  const double &t ) { Timeu = t; }
  void SetTimed(  const double &t ) { Timed = t; }
  void SetCTimeu( const double &t ) { CTimeu = t; }
  void SetCTimed( const double &t ) { CTimed = t; }
  void SetTMean(  const double &t ) { TMean = t; }
  void SetCTMean( const double &t ) { CTMean = t; }
  void SetTSub(   const double &t ) { TSub = t; }
  void SetCTSub(  const double &t ) { CTSub = t; }
  void SetPos( const double &x, const double &y, const double &z ) { X=x;  Y=y;  Z=z; }
  void SetDir( const double &x, const double &y, const double &z ) { dX=x; dY=y; dZ=z; }
  void SetGPos( const double &x, const double &y, const double &z ) { GX=x;  GY=y;  GZ=z; }
  void SetGDir( const double &x, const double &y, const double &z ) { dGX=x; dGY=y; dGZ=z; }
  void SetLength( const double &l ) { Length = l; }
  void SetWidth(  const double &w ) { Width  = w; }
  void SetThick(  const double &t ) { Thick  = t; }
  void SetLightVelocity( const double &lv ) { LightVelocity = lv; }
  void SetHitPosition( const double &l ) { HitPosition = l; }
  
  bool Calc( ConfMan *conf );
  void Clear();

  // simulation
 private:
  double dT; // resolution for simulation
  double dE; // resolution for simulation

 public:
  double dt() const { return dT; }
  double de() const { return dE; }
  bool SetSimulatedResolution( ConfMan *conf );
  void Reverse( ConfMan *conf );

  ClassDef(HodoscopeLikeHit, 1 );
};

#endif
