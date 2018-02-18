// CherenkovLikeHit.cpp

#include <iostream>
#include <iomanip>
#include <string>
#include <new>
#include <cstring>

#include "CherenkovLikeHit.h"

ClassImp(CherenkovLikeHit);

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
CherenkovLikeHit::CherenkovLikeHit()
{
  Clear();
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
CherenkovLikeHit::CherenkovLikeHit( const CherenkovLikeHit &hit )
{
  *this = hit; // I cannot understand why I need this copy constructor.
}

bool CherenkovLikeHit::SetData( const int &c, const int &n, const int &a,
				const int &cid, const int &seg, const int &at, const int &ud,
				const int &data )
{

#if 0
  std::cout << c << " " << n << " " << a << " " << cid << " " << seg << " " << at << " " << ud << " " << data 
	    << std::endl;
#endif
  if(seg==0){
    std::cout << c << " " << n << " " << a << " " << cid << " " << seg << " " << at << " " << ud << " " << data 
	      << std::endl;
  }
  CounterID = cid;
  Seg = seg;
  if( at==0 ){
    ADC[ud]=data; CrateA[ud]=c; SlotA[ud]=n; ChannelA[ud]=a; 
  }
  else if(at==1){
    TDC[ud]=data; CrateT[ud]=c; SlotT[ud]=n; ChannelT[ud]=a; 
  }
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void CherenkovLikeHit::Clear()
{
  for(int i=0;i<4;i++){
    CrateT[i] = SlotT[i] = ChannelT[i] = -1;
    CrateA[i] = SlotA[i] = ChannelA[i] = -1;
    ADC[i] = 0;
    TDC[i] = 0;
    NumPhoton[i] = 0;
    Time[i] = -999;
  }
  Seg = -1;
  X = Y = Z = dX = dY = dZ = -999;
  Length = Width = Thick = -999;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool CherenkovLikeHit::CheckRange( int num_required )
{
  bool status = true;
  
  for(int i=0;i<4;i++){
    if( 0<TDC[i] && TDC[i]<4095 ){}
    else{ status = false; }
    if( num_required==i+1 ) return status;
  }
  std::cout << " invalid value of #required TDC " << std::endl;
  return false;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool CherenkovLikeHit::Calc( ConfMan *conf )
{
  //  std::cout<<"CherenkovLikeHit::Calc "<<CounterID<<std::endl;
  if( conf->GetGeomMapManager() ){
    double x,y,z,dx,dy,dz,len,wid,thi,lv;
    if( conf->GetGeomMapManager()->GetParam( CounterID, Seg, x, y, z, dx, dy, dz, len, wid, thi, lv ) ){
      X = x; Y = y; Z = z; dX = dx; dY = dy; dZ = dz;
      Length = len; Width = wid; Thick = thi;
    }
    if( conf->GetGeomMapManager()->GetGParam( CounterID, x, y, z, dx, dy, dz ) ){
      GX = x; GY = y; GZ = z; dGX = dx; dGY = dy; dGZ = dz;
    }
  }
  if( conf->GetGainMapManager() ){
    for(int i=0;i<4;i++){
      if( 0<ADC[i] ) NumPhoton[i] = conf->GetGainMapManager()->CalcCValue( CrateA[i], SlotA[i], ChannelA[i], 0, ADC[i], 1 );
      
      //      std::cout<<"ADC="<<ADC[i]<<"\tNumPhoton="<<NumPhoton[i]<<std::endl;
      
      if( 0<TDC[i] ) Time[i]      = conf->GetGainMapManager()->CalcCValue( CrateT[i], SlotT[i], ChannelT[i], 1, TDC[i] );
    }
  }
  
  return true;
}
