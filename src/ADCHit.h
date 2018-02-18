#ifndef ADCHit_h
#define ADCHit_h 1

#include <string>
#include <iostream>
#include <cmath>

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#include "ConfMan.h"
#include "TKO.h"

class ADCHit : public TKOHit
{
 public:
  ADCHit();
  //  ADCHit( const ADCHit &hit );

 private:
  // conversion data
  double Energy;

 public:
  double energy()   const { return  Energy; }

  void SetEnergy( const double &e ) { Energy = e; }

  bool Calc( ConfMan *conf );
  void Clear();

  ClassDef(ADCHit, 1 );
};

#endif
