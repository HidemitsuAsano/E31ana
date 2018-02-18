// TDCHit.cpp
#include <iostream>
#include <iomanip>
#include <string>
#include <new>
#include <cstring>

#include "TDCHit.h"

ClassImp(TDCHit);

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
TDCHit::TDCHit(): Time(-999.),CTime(-999.)
{
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void TDCHit::Clear()
{
  Time = CTime =-999;
}


// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool TDCHit::Calc( ConfMan *conf )
{
#if 0
  std::string name = conf->GetCounterMapManager()->GetName(Crate,Slot,Channel);
  std::cout << name << "  " << Crate << "  " << Slot << "  " << Channel << "  " <<  Data << std::endl;
#endif    
  if( conf->GetGainMapManager() ){
    if( 0<Data ) Time = conf->GetGainMapManager()->CalcCValue( Crate, Slot, Channel, 1, Data, 1);
  }
  return true;

}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool TDCHit::Reverse( ConfMan *conf )
{
#if 0
  std::string name = conf->GetCounterMapManager()->GetName(Crate,Slot,Channel);
  std::cout << name << "  " << Crate << "  " << Slot << "  " << Channel << "  " <<  Data << std::endl;
#endif    
  if( conf->GetGainMapManager() ){
    Data = conf->GetGainMapManager()->CalcDATValue( Crate, Slot, Channel, 1, Time ) ;
  }
  return true;

}
