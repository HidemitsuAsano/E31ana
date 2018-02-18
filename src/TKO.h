// TKO.h

#ifndef TKO_h
#define TKO_h 1

#include <cstddef>
#include <vector>
#include <new>
#include <algorithm>
#include <iostream>

#ifndef ROOT_TObject
#include "TObject.h"
#endif

class TKOHit : public TObject
{
 protected:
  int Crate, Slot, Channel, Data;

 public:
  TKOHit();
  //  TKOHit( TKOHit *hit );
  TKOHit( const int &cr, const int &sl, const int &ch, const int &data );
  virtual ~TKOHit() {}

 public:
  int cr() const { return Crate; }
  int sl() const { return Slot; }
  int ch() const { return Channel; }
  int data() const { return Data; }

  void SetCrate( const int &cr ) { Crate = cr; }
  void SetSlot(  const int &sl ) { Slot = sl; }
  void SetChannel( const int &ch ) { Channel = ch; }
  void SetData( const int &data ) { Data = data; }

  void SetHit( const int &cr, const int &sl, const int &ch, const int &data )
  { Crate = cr; Slot = sl; Channel = ch; Data = data; }

  void Clear();
  bool CheckRange(const int &ll=0,const int &ul=4095) const;

  ClassDef( TKOHit, 1 );
};

class TKOHitCollection : public TObject
{
 private:
  std::vector <TKOHit > theCollection;

 public:
  TKOHitCollection();
  virtual ~TKOHitCollection() {}
  
 public:
  int entries() const { return theCollection.size(); }
  TKOHit *hit( const int &i );
  void AddHit( const TKOHit &hit ) { theCollection.push_back(hit); }

  void Clear();

  ClassDef( TKOHitCollection, 1 );
};

#endif
