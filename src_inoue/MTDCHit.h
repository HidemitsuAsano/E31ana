// MTDCHit.h

#ifndef MTDCHit_h
#define MTDCHit_h 1

#include <cstddef>
#include <vector>
#include <new>
#include <algorithm>
#include <iostream>

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#include "ConfMan.h"

class aMTDCHit : public TObject
{
 private:
  int TDC;
  double Time,CTime;

 public:
  aMTDCHit();
  aMTDCHit( const int &data );
  aMTDCHit( const int &data, const double &time, const double &ctime);
  virtual ~aMTDCHit() {}

 public:
  int tdc() const { return TDC; }
  double time() const { return Time; }
  double ctime() const { return CTime; }
  void SetTDC( const int &t ) { TDC=t; }
  void SetTime( const double &t ) { Time=t; }
  void SetCTime( const double &t ) { CTime=t; }

  bool CheckRange(const int &ll=0,const int &ul=20000);
  ClassDef( aMTDCHit, 1 );
};

class MTDCHit : public TObject
{
 private:
  int Crate, Slot, Channel;
  int CounterID;
  int HitID;
  int Seg;
  int NSensor;
  std::vector< aMTDCHit > Data[2];

 public:
  MTDCHit();
  MTDCHit( const int &cr, const int &sl, const int &ch, 
	   const int &cid, const int &seg, const int &ud,
	   const int &data );
  virtual ~MTDCHit() {}

 public:
  int cid()  const { return CounterID; }
  int hid()  const { return HitID; }
  int seg()  const { return Seg; }

  int cr() const { return Crate; }
  int sl() const { return Slot; }
  int ch() const { return Channel; }
  int tdc(const int &ud,const int &i) const { return ( i>=0 && i<(int)Data[ud].size() )? Data[ud][i].tdc() : -1; }
  double time(const int &ud,const int &i) const { return ( i>=0 && i<(int)Data[ud].size() )? Data[ud][i].time() : -999; }
  double ctime(const int &ud,const int &i) const { return ( i>=0 && i<(int)Data[ud].size() )? Data[ud][i].ctime() : -999; }
  int ndata(const int &ud) const { return Data[ud].size(); }

  void SetCounterID(const int & cid ) { CounterID = cid; }
  void SetHitID(const int & hid ) { HitID = hid; }
  void SetSegment(  const int & seg ) { Seg = seg; }
  void SetCrate( const int &cr ) { Crate = cr; }
  void SetSlot(  const int &sl ) { Slot = sl; }
  void SetChannel( const int &ch ) { Channel = ch; }
  void SetData( const int &ud, const int &data ) { Data[ud].push_back( data ); }

  void SetData( const int &cr, const int &sl, const int &ch,
	       const int &cid, const int &seg, const int &ud,const int &data )
  { Crate = cr; Slot = sl; Channel = ch; 
    CounterID=cid; Seg=seg;
    Data[ud].push_back( aMTDCHit( data ) ); }

  void Clear();
  void Calc( ConfMan *conf, const int &ref );
  int CheckRange(const int &ll=0,const int &ul=20000);

  ClassDef( MTDCHit, 1 );
};

#endif
