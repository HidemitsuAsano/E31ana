  // ADCHit.cpp

#include <iostream>
#include <iomanip>
#include <string>
#include <new>
#include <cstring>

#include "ADCHit.h"

ClassImp(ADCHit);

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
ADCHit::ADCHit():  Energy(-999)
{
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void ADCHit::Clear()
{
  Energy = -999;
}


// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool ADCHit::Calc( ConfMan *conf )
{
#if 0
  std::string name = conf->GetCounterMapManager()->GetName(Crate,Slot,Channel);
  std::cout << name << "  " << Crate << "  " << Slot << "  " << Channel << "  " <<  Data << std::endl
#endif    
  if( conf->GetGainMapManager() ){
    if( 0<Data ) Energy  = conf->GetGainMapManager()->CalcCValue( Crate, Slot, Channel, 0, Data, 1);
  }

  return true;

}
