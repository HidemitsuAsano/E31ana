// HodoscopeLikeHit.cpp

#include <iostream>
#include <iomanip>
#include <string>
#include <new>
#include <cstring>

#include "HodoscopeLikeHit.h"

ClassImp(HodoscopeLikeHit);

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
HodoscopeLikeHit::HodoscopeLikeHit()
{
  Clear();
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
HodoscopeLikeHit::HodoscopeLikeHit( const HodoscopeLikeHit &hit )
{
  *this = hit; // I cannot understand why I need this copy constructor.
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void HodoscopeLikeHit::Clear()
{
  CrateTu = SlotTu = ChannelTu = -1;
  CrateTd = SlotTd = ChannelTd = -1;
  CrateAu = SlotAu = ChannelAu = -1;
  CrateAd = SlotAd = ChannelAd = -1;
  Seg = -1;
  ADCu = ADCd = TDCu = TDCd = -1;

  Eneu = Ened = EMean = Timeu = Timed = CTimeu = CTimed = -999;
  X = Y = Z = dX = dY = dZ = -999;
  Length = Width = Thick = -999;
  TMean = CTMean = TSub = CTSub = -999.;
  LightVelocity = 15.;
  HitPosition = 0;

  CHECKRANGE=false;

  dT = dE = 0;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool HodoscopeLikeHit::CheckRange(const int &type)
{
  if(type==0){
    if( 0<TDCu && TDCu<4095 &&
	0<TDCd && TDCd<4095 )
      return true;
    else
      return false;
  }else if(type==1){
    if( 0<TDCu && TDCu<4095 )
      return true;
    else
      return false;
  }else
    return false;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool HodoscopeLikeHit::CheckRange(const int &ll, const int &ul,const bool &type)
{
  if(type==false){
    if( ll<TDCu && TDCu<ul &&
	ll<TDCd && TDCd<ul )
      return true;
    else
      return false;
  }
  else{
    if( ll<TDCu && TDCu<ul )
      return true;
    else
      return false;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool HodoscopeLikeHit::CheckRange2( ConfMan* conf)
{
  if( conf ){
    CHECKRANGE=true;
    if(conf->GetGateMapManager()){
      if(!conf->GetGateMapManager()->CheckRange( CounterID, Seg, 0, 0, 0, (double)ADCu )) CHECKRANGE=false; 
      else if(!conf->GetGateMapManager()->CheckRange( CounterID, Seg, 0, 1, 0, (double)TDCu )) CHECKRANGE=false; 
      else if(!conf->GetGateMapManager()->CheckRange( CounterID, Seg, 1, 0, 0, (double)ADCd )) CHECKRANGE=false; 
      else if(!conf->GetGateMapManager()->CheckRange( CounterID, Seg, 1, 1, 0, (double)TDCd )) CHECKRANGE=false; 
    }
    return CHECKRANGE;
  }
  else return CheckRange();
}

bool HodoscopeLikeHit::SetData( const int &c, const int &n, const int &a,
				const int &cid, const int &seg, const int &at, const int &ud,
				const int &data )
{
#if 0
  std::cout << c << " " << n << " " << a << " " << cid << " " << seg << " " << at << " " << ud << " " << data 
	    << std::endl;
#endif
  CounterID = cid;
  Seg = seg;
  if( at==0 ){
    if( ud==0 ){ ADCu=data; CrateAu=c; SlotAu=n; ChannelAu=a; }
    else{        ADCd=data; CrateAd=c; SlotAd=n; ChannelAd=a; }
  }
  else if( at==1 ){
    if( ud==0 ){ TDCu=data; CrateTu=c; SlotTu=n; ChannelTu=a; }
    else{        TDCd=data; CrateTd=c; SlotTd=n; ChannelTd=a; }
  }  
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool HodoscopeLikeHit::SetSimulatedResolution( ConfMan *conf )
{
  if( conf->GetReslMapManager() ){
    CTSub = HitPosition/LightVelocity;
    conf->GetReslMapManager()->GetResolution( CounterID, 0, Seg, dT, dE ); // first call for CTSub. just for a calculation. maybe not so important.
    CTSub  += dT;
    conf->GetReslMapManager()->GetResolution( CounterID, 0, Seg, dT, dE ); // this value is most important and saved.
    CTMean += dT;
    EMean  += dE; // this treatment is not so good. but maybe, everyone 
    Reverse(conf);
    return true;
  }
  return false;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void HodoscopeLikeHit::Reverse( ConfMan *conf )
{
  CTimeu = CTMean + CTSub;
  CTimed = CTMean - CTSub;
  Eneu = EMean*exp( -HitPosition/200. ); // assuming the attenuation lenght of 200 cm. maybe, geom map man should deal with this value.
  Ened = EMean*exp( +HitPosition/200. ); // assuming the attenuation lenght of 200 cm. maybe, geom map man should deal with this value.

  if( conf->GetSlewingMapManager() ){
    Timeu = conf->GetSlewingMapManager()->CalcDATValue( CID_CDH, Seg, 0, CTimeu, Eneu );
    Timed = conf->GetSlewingMapManager()->CalcDATValue( CID_CDH, Seg, 1, CTimed, Ened );
  }
  if( conf->GetGainMapManager() ){
    ADCu = conf->GetGainMapManager()->CalcDATValue( CrateAu, SlotAu, ChannelAu, 0, Eneu );
    ADCd = conf->GetGainMapManager()->CalcDATValue( CrateAd, SlotAd, ChannelAd, 0, Ened );
    TDCu = conf->GetGainMapManager()->CalcDATValue( CrateTu, SlotTu, ChannelTu, 1, Timeu );
    TDCd = conf->GetGainMapManager()->CalcDATValue( CrateTd, SlotTd, ChannelTd, 1, Timed );
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool HodoscopeLikeHit::Calc( ConfMan *conf )
{
#if 0
  std::string nameau = conf->GetCounterMapManager()->GetName(CrateAu,SlotAu,ChannelAu);
  std::string namead = conf->GetCounterMapManager()->GetName(CrateAd,SlotAd,ChannelAd);
  std::string nametu = conf->GetCounterMapManager()->GetName(CrateTu,SlotTu,ChannelTu);
  std::string nametd = conf->GetCounterMapManager()->GetName(CrateTd,SlotTd,ChannelTd);
  std::cout << nameau << "  " << CrateAu << "  " << SlotAu << "  " << ChannelAu << "  " <<  ADCu << std::endl
	    << namead << "  " << CrateAd << "  " << SlotAd << "  " << ChannelAd << "  " <<  ADCd << std::endl
	    << nametu << "  " << CrateTu << "  " << SlotTu << "  " << ChannelTu << "  " <<  TDCu << std::endl
	    << nametd << "  " << CrateTd << "  " << SlotTd << "  " << ChannelTd << "  " <<  TDCd << std::endl;
#endif    

  if( conf->GetGeomMapManager() ){
    double x,y,z,dx,dy,dz,len,wid,thi,lv;
    if( conf->GetGeomMapManager()->GetParam( CounterID, Seg, x, y, z, dx, dy, dz, len, wid, thi, lv ) ){
      X = x; Y = y; Z = z; dX = dx; dY = dy; dZ = dz;
      Length = len; Width = wid; Thick = thi;
      LightVelocity = lv;
    }


    if( conf->GetGeomMapManager()->GetGParam( CounterID, x, y, z, dx, dy, dz ) ){
      GX = x; GY = y; GZ = z; dGX = dx; dGY = dy; dGZ = dz;
    }
  }

  if( conf->GetGainMapManager() ){
    if( 0<ADCu ) Eneu  = conf->GetGainMapManager()->CalcCValue( CrateAu, SlotAu, ChannelAu, 0, ADCu, 1);
    if( 0<ADCd ) Ened  = conf->GetGainMapManager()->CalcCValue( CrateAd, SlotAd, ChannelAd, 0, ADCd, 1);
    if( 0<TDCu ) Timeu = conf->GetGainMapManager()->CalcCValue( CrateTu, SlotTu, ChannelTu, 1, TDCu, 1);
    if( 0<TDCd ) Timed = conf->GetGainMapManager()->CalcCValue( CrateTd, SlotTd, ChannelTd, 1, TDCd, 1);
  }

  if( conf->GetSlewingMapManager() ){
    if( 0<TDCu ) CTimeu = conf->GetSlewingMapManager()->CalcCValue( CounterID, Seg, 0, Timeu, Eneu );
    if( 0<TDCd ) CTimed = conf->GetSlewingMapManager()->CalcCValue( CounterID, Seg, 1, Timed, Ened );
  }

  EMean  = sqrt(Eneu*Ened);
  TMean  = ( Timeu+ Timed)/2.;
  CTMean = (CTimeu+CTimed)/2.;
  TSub   = ( Timeu- Timed)/2.;
  CTSub  = (CTimeu-CTimed)/2.;
  HitPosition = -CTSub*LightVelocity;

  CheckRange2( conf );

  return true;
}
