// HodoscopeLikeHit.cpp

#include <iostream>
#include <iomanip>
#include <string>
#include <new>
#include <cstring>

#include "HodoscopeLikeHit.h"

ClassImp(HodoscopeLikeHit);

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
HodoscopeLikeHit::HodoscopeLikeHit():
  CounterID(-1),HitID(-1),  Seg(-1),NSensor(2),
  TDCMean(-999),TMean(-999),CTMean(-999),TSub(-999),CTSub(-999),
  EMean(-999),
  HitPos(DEFVECT), HitPosition(-999),
  dT(0.),dE(0.)
{
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void HodoscopeLikeHit::Clear()
{
  Seg = -1;
  NSensor = 2;

  EMean = -999;
  TDCMean = TMean = CTMean = TSub = CTSub = -999.;
  // X = Y = Z = dX = dY = dZ = -999;
  // Length = Width = Thick = -999;
  HitPosition = 0;

  dT = dE = 0;
  for(int i=0;i<2;i++){
    adchit[i].Clear();
    tdchit[i].Clear();
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool HodoscopeLikeHit::CheckRange() const
{
  for(int i=0;i<NSensor;i++)
    if(!tdchit[i].CheckRange())
      return false;
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool HodoscopeLikeHit::CheckRange(const int &ll, const int &ul,const bool &type) const
{
  for(int i=0;i<NSensor;i++)
    if(!tdchit[i].CheckRange(ll,ul))
      return false;
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool HodoscopeLikeHit::CheckRange2( ConfMan* conf) const
{
  if( conf ){
    // for(int i=0;i<NSensor;i++)
    //   if(!tdchit[i].CheckRange(conf))
    //  	return false;
    // for(int i=0;i<NSensor;i++)
    //   if(!adchit[i].CheckRange(conf))
    // 	return false;
    return true;
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
  if( at==0 ) adchit[ud].SetHit(c,n,a,data);
  else if( at==1 ) tdchit[ud].SetHit(c,n,a,data);
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool HodoscopeLikeHit::SetSimulatedResolution( ConfMan *conf )
{
  if( conf->GetReslMapManager() ){
    //    CTSub = HitPosition/LightVelocity; //  temporaly commeted out
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
  tdchit[0].SetCTime( CTMean + CTSub );
  tdchit[1].SetCTime( CTMean - CTSub );
 // assuming the attenuation lenght of 200 cm. maybe, geom map man should deal with this value.
  adchit[0].SetEnergy( EMean*exp( -HitPosition/200. ) );
  adchit[1].SetEnergy( EMean*exp( +HitPosition/200. ) );

  SlewingMapMan* slewman=conf->GetSlewingMapManager();
  if( slewman ){
    for(int i=0;i<2;i++)
      tdchit[i].SetTime(slewman->CalcDATValue( CounterID, Seg, i,ctime(i), ene(i) ));
  }
  if( conf->GetGainMapManager() )
    for(int at=0;at<2;at++)
      for(int i=0;i<2;i++){
	double val=ene(i);
	if(at==1) val=time(i); 
	//	SetVal(at,i,conf->GetGainMapManager()->CalcDATValue( cr(at,i), sl(at,i), ch(at,i), at, val ));
	SetVal(at,i,1);
      }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool HodoscopeLikeHit::Calc( ConfMan *conf )
{
  if( conf->GetGainMapManager() )
    for(int i=0;i<NSensor;i++){
      adchit[i].Calc(conf);
      tdchit[i].Calc(conf);
    }  

  if( conf->GetSlewingMapManager() )
    for(int i=0;i<2;i++)
      if( tdchit[i].CheckRange() && adchit[i].data()>0 ){
	double ctime=conf->GetSlewingMapManager()->CalcCValue( CounterID, Seg, i, tdchit[i].time() , adchit[i].energy() );
	tdchit[i].SetCTime(ctime);
      }
  
  EMean  = sqrt(adchit[0].energy()*adchit[1].energy());
  TDCMean  = (tdchit[0].data() + tdchit[1].data())/2.;
  TMean  = (tdchit[0].time() + tdchit[1].time())/2.;
  CTMean = (tdchit[0].ctime() + tdchit[1].ctime())/2.;
  TSub   = (tdchit[0].time() - tdchit[1].time())/2.;
  CTSub   = (tdchit[0].ctime() - tdchit[1].ctime())/2.;
  double lv;
  if( conf->GetGeomMapManager() ){
    //    conf->GetGeomMapManager()->GetGPos( CounterID, Seg, -CTSub , HitPos);
    conf->GetGeomMapManager()->GetGPos( CounterID, Seg , HitPos);
    if(conf->GetGeomMapManager()->GetLightVelocity( CounterID, Seg, lv))
      HitPosition = -CTSub*lv;
  }

  return true;
}
