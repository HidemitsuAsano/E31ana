// SDDHitCom.cpp

#include <iostream>
#include <iomanip>
#include <string>
#include <new>
#include <cstring>

#include "SDDHitCom.h"

ClassImp(SDDHitCom);

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
SDDHitCom::SDDHitCom()
{
  Clear();
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
SDDHitCom::SDDHitCom( const SDDHitCom &hit )
{
  *this = hit; // I cannot understand why I need this copy constructor.
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void SDDHitCom::Clear()
{
  CounterID=-1;

  CrateT0t = SlotT0t = ChannelT0t = -1;
  CrateE0t1 = SlotE0t1 = ChannelE0t1 = -1;
  CrateE0t2 = SlotE0t2 = ChannelE0t2 = -1;
  CrateE0t3 = SlotE0t3 = ChannelE0t3 = -1;
  CrateFirst1 = SlotFirst1 = ChannelFirst1 = -1;
  CrateFirst2 = SlotFirst2 = ChannelFirst2 = -1;

  T0t = E0t1 = E0t2 = E0t3 = First1 = First2 = -1;


  T0time = E0time1 = E0time2 = E0time3 = Firstime1 = Firstime2 = -999;
  CT0time = CE0time1 = CE0time2 = CE0time3 = CFirstime1 = CFirstime2 = -999;

}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
/*
bool SDDHitCom::CheckRange()
{
  if( 0<TDCsdd && TDCsdd<4095
      )
    return true;
  else
    return false;
}
*/
bool SDDHitCom::SetData( const int &c, const int &n, const int &a,
			 const int &cid, const int &seg, const int &at, const int &ud,
			 const int &data )
{
#if 0
  std::cout << c << " " << n << " " << a << " " << cid << " " << seg << " " << at << " " << ud << " " << data 
	    << std::endl;
#endif
  CounterID = cid;
  if( seg == 0 ){
    if(ud==0){  T0t=data; CrateT0t=c; SlotT0t=n; ChannelT0t=a; } 
    else if(ud==1){  E0t1=data; CrateE0t1=c; SlotE0t1=n; ChannelE0t1=a; } 
    else if(ud==2){  E0t2=data; CrateE0t2=c; SlotE0t2=n; ChannelE0t2=a; } 
    else if(ud==3){  E0t3=data; CrateE0t3=c; SlotE0t3=n; ChannelE0t3=a; } 
    else if(ud==4){  First1=data; CrateFirst1=c; SlotFirst1=n; ChannelFirst1=a; } 
    else{  First2=data; CrateFirst2=c; SlotFirst2=n; ChannelFirst2=a; } 
  }
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool SDDHitCom::Calc( ConfMan *conf )
{
  if( conf->GetGainMapManager() ){
    if( 0<T0t ) T0time = conf->GetGainMapManager()->CalcCValue( CrateT0t, SlotT0t, ChannelT0t, 1, T0t );
    if( 0<E0t1 ) E0time1 = conf->GetGainMapManager()->CalcCValue( CrateE0t1, SlotE0t1, ChannelE0t1, 1, E0t1 );
    if( 0<E0t2 ) E0time2 = conf->GetGainMapManager()->CalcCValue( CrateE0t2, SlotE0t2, ChannelE0t2, 1, E0t2 );
    if( 0<E0t3 ) E0time3 = conf->GetGainMapManager()->CalcCValue( CrateE0t3, SlotE0t3, ChannelE0t3, 1, E0t3 );
    if( 0<First1 ) Firstime1 = conf->GetGainMapManager()->CalcCValue( CrateFirst1, SlotFirst1, ChannelFirst1, 1, First1 );
    if( 0<First2 ) Firstime2 = conf->GetGainMapManager()->CalcCValue( CrateFirst2, SlotFirst2, ChannelFirst2, 1, First2 );
  }
  /*
  if( conf->GetSlewingMapManager() ){
    if( 0<T0t ) CT0time = conf->GetSlewingMapManager()->CalcCValue( CounterID, 0, 0, T0time, Out );
    if( 0<E0t1 ) CE0time1 = conf->GetSlewingMapManager()->CalcCValue( CounterID, 0, 1, E0time1, Out );
    if( 0<E0t2 ) CE0time2 = conf->GetSlewingMapManager()->CalcCValue( CounterID, 0, 2, E0time2, Out );
    if( 0<E0t3 ) CE0time3 = conf->GetSlewingMapManager()->CalcCValue( CounterID, 0, 3, E0time3, Out );
    if( 0<First1 ) CFirstime1 = conf->GetSlewingMapManager()->CalcCValue( CounterID, 0, 4, Firstime1, Out );
    if( 0<First2 ) CFirstime2 = conf->GetSlewingMapManager()->CalcCValue( CounterID, 0, 5, Firstime2, Out );
  }
  */
  return true;
}
