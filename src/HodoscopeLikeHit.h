#ifndef HodoscopeLikeHit_h
#define HodoscopeLikeHit_h 1

#include <string>
#include <iostream>
#include <cmath>

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#include "ConfMan.h"
#include "ADCHit.h"
#include "TDCHit.h"

class HodoscopeLikeHit : public TObject
{
 private:
  int CounterID;
  int HitID; 
  int Seg;
  int NSensor;
  
 public:
  HodoscopeLikeHit();
  //  HodoscopeLikeHit( const HodoscopeLikeHit &hit );
  virtual ~HodoscopeLikeHit(){};

 private:
  // raw data
  ADCHit adchit[2];
  TDCHit tdchit[2];

 public:
  int cid()  const { return CounterID; }
  int hid()  const { return HitID; }
  int seg()  const { return Seg; }
  inline int val(const int &at, const int &i);
  int adc(const int &i) const { return (i>=0 && i<2 ) ? adchit[i].data() : -1; }
  int tdc(const int &i) const { return (i>=0 && i<2 ) ? tdchit[i].data() : -1; }

  int adcu() const { return adchit[0].data(); }
  int adcd() const { return adchit[1].data(); }
  int tdcu() const { return tdchit[0].data(); }
  int tdcd() const { return tdchit[1].data(); }

  inline int cr(const int &at, const int &i);
  inline int sl(const int &at, const int &i);
  inline int ch(const int &at, const int &i);
  int nsensor() const { return NSensor; }
  void SetCounterID(const int & cid ) { CounterID = cid; }
  void SetHitID(    const int & hid ) { HitID = hid; }
  void SetSegment(  const int & seg ) { Seg = seg; }
  void SetNSensor(  const int & ns ) { NSensor = ns; }
  inline void SetVal(      const int &at, const int &i, const int & val );
  inline void SetCrate(    const int &at, const int &i, const int & cr  );
  inline void SetSlot(     const int &at, const int &i, const int & sl  );
  inline void SetChannel(  const int &at, const int &i, const int & ch  );
  
  bool SetData( const int &c, const int &n, const int &a,
		const int &cid, const int &seg, const int &at, const int &ud,
		const int &data );

 private:
  // conversion data
  double TDCMean, TMean, CTMean;
  double TSub, CTSub;
  double EMean;

 public:
  double ene(const int &i)   { return (i>=0 && i<2 ) ? adchit[i].energy() : 0; }
  double time(const int &i)  { return (i>=0 && i<2 ) ? tdchit[i].time() : 0; }
  double ctime(const int &i) { return (i>=0 && i<2 ) ? tdchit[i].ctime() : 0; }
  double eu() { return adchit[0].energy(); }
  double ed() { return adchit[1].energy(); }
  double tu() { return tdchit[0].time(); }
  double td() { return tdchit[1].time(); }
  double ctu() { return tdchit[0].ctime(); }
  double ctd() { return tdchit[1].ctime(); }

  double emean()const { return EMean; }
  double tdcmean() const { return  TDCMean; }
  double tmean() const { return  TMean; }
  double ctmean()const { return CTMean; }
  double tsub()  const { return TSub; }
  double ctsub() const { return CTSub; }

  void SetEne(   const int &i,  const double &e ) { adchit[i].SetEnergy(e); }
  void SetEMean( const double &e ) { EMean = e; }
  void SetTime(  const int &i,  const double &t ) { tdchit[i].SetTime(t); }
  void SetCTime( const int &i, const double &t ) { tdchit[i].SetCTime(t); }
  void SetTMean( const double &t ) { TMean = t; }
  void SetCTMean(const double &t ) { CTMean = t; }
  void SetTSub(  const double &t ) { TSub = t; }
  void SetCTSub( const double &t ) { CTSub = t; }

 private:
  TVector3 HitPos;
  double HitPosition; // hit position from the center along scinti.
  /* double  X,  Y,  Z; // position */
  /* double dX, dY, dZ; // direction */
  /* double  GX,  GY,  GZ; // origin position */
  /* double dGX, dGY, dGZ; // origin direction */
  /* double Length, Width, Thick; */  
  //  double LightVelocity;
  
 public:
  TVector3 pos() const { return HitPos; }
  double x()   const { return HitPos.X(); }
  double y()   const { return HitPos.Y(); }
  double z()   const { return HitPos.Z(); }
  
  /* inline double pos( int &i ); */
  /* inline double dpos( int &i ); */
  /* double dx()  const { return dX; } */
  /* double dy()  const { return dY; } */
  /* double dz()  const { return dZ; } */
  /* inline double gpos( int &i ); */
  /* double gx()  const { return GX; } */
  /* double gy()  const { return GY; } */
  /* double gz()  const { return GZ; } */
  /* inline double dgpos( int &i ); */
  /* double dgx() const { return dGX; } */
  /* double dgy() const { return dGY; } */
  /* double dgz() const { return dGZ; } */
  /* double len() const { return Length; } */
  /* double wid() const { return Width; } */
  /* double thick() const { return Thick; } */
  /* void pos( double &x, double &y, double &z ) { x=X;  y=Y;  z=Z; } */
  /* void dir( double &x, double &y, double &z ) { x=dX; y=dY; z=dZ; } */
  /* void gpos( double &x, double &y, double &z ) { x=GX;  y=GY;  z=GZ; } */
  /* void gdir( double &x, double &y, double &z ) { x=dGX; y=dGY; z=dGZ; } */
  
  void SetPos( const TVector3 &pos ){ HitPos=pos; }
  void SetPos( const double &x, const double &y, const double &z ) { HitPos.SetXYZ(x,y,z); } 
  double hitpos() const { return HitPosition; }
  void SetHitPosition( const double &l ) { HitPosition = l; }

  /* void SetDir( const double &x, const double &y, const double &z ) { dX=x; dY=y; dZ=z; } */
  /* void SetGPos( const double &x, const double &y, const double &z ) { GX=x;  GY=y;  GZ=z; } */
  /* void SetGDir( const double &x, const double &y, const double &z ) { dGX=x; dGY=y; dGZ=z; } */
  /* void SetLength( const double &l ) { Length = l; } */
  /* void SetWidth(  const double &w ) { Width  = w; } */
  /* void SetThick(  const double &t ) { Thick  = t; } */


 public:  
  bool Calc( ConfMan *conf );
  bool CheckRange() const;
  bool CheckRange(const int &ll,const int &ul, const bool &type=false) const;
  bool CheckRange2(ConfMan *conf=0) const;
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

inline int HodoscopeLikeHit::val(const int &at, const int &i){
  if(at==0) 
    return adc(i);
  else if(at==1)
    return tdc(i);
  else
    return -1;
}

inline int HodoscopeLikeHit::cr(const int &at, const int &i){
  if(at==0)
    return (i>=0 && i<2 ) ? adchit[i].cr() : -1;
  else if(at==1) 
    return (i>=0 && i<2 ) ? tdchit[i].cr() : -1;
  else 
    return -1;
}

inline int HodoscopeLikeHit::sl(const int &at, const int &i){
  if(at==0)
    return (i>=0 && i<2 ) ? adchit[i].sl() : -1;
  else if(at==1)
    return (i>=0 && i<2 ) ? tdchit[i].sl() : -1;
  else
    return -1;
}

inline int HodoscopeLikeHit::ch(const int &at, const int &i){
  if(at==0)
    return (i>=0 && i<2 ) ? adchit[i].ch() : -1;
  else if(at==1)
    return (i>=0 && i<2 ) ? tdchit[i].ch() : -1;
  else
    return -1;
}

inline void HodoscopeLikeHit::SetVal(const int &at, const int &i, const int &val){
  if (i>=0 && i<2 ){
    if(at==0)      adchit[i].SetData(val); 
    else if(at==1) tdchit[i].SetData(val);
  }
}

inline void HodoscopeLikeHit::SetCrate(const int &at, const int &i, const int &val){
  if (i>=0 && i<2 ){
    if(at==0)      adchit[i].SetCrate(val); 
    else if(at==1) tdchit[i].SetCrate(val);
  }
}

inline void HodoscopeLikeHit::SetSlot(const int &at, const int &i, const int &val){
  if (i>=0 && i<2 ){
    if(at==0)      adchit[i].SetSlot(val); 
    else if(at==1) tdchit[i].SetSlot(val);
  }
}

inline void HodoscopeLikeHit::SetChannel(const int &at, const int &i, const int &val){
  if (i>=0 && i<2 ){
    if(at==0)      adchit[i].SetChannel(val); 
    else if(at==1) tdchit[i].SetChannel(val);
  }
}

#endif
