#ifndef TDCHit_h
#define TDCHit_h 1

#include <string>
#include <iostream>
#include <cmath>

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#include "ConfMan.h"
#include "TKO.h"

class TDCHit : public TKOHit
{
 public:
  TDCHit();
  //  TDCHit( const TDCHit &hit );
  virtual ~TDCHit() {};

 private:
  // conversion data
  double Time,CTime; // corrected time

 public:
  double time()  const { return  Time; }
  double ctime() const { return  CTime; }

  void SetTime(  const double &t ) { Time = t; }
  void SetCTime( const double &t ) { CTime = t; }

  bool Calc( ConfMan *conf );
  bool Reverse( ConfMan *conf );
  void Clear();

  ClassDef(TDCHit, 1 );
};

#endif
