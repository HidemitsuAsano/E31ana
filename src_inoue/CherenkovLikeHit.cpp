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
  if( at==0 ) adchit[ud].SetHit(c,n,a,data);
  else if( at==1 ) tdchit[ud].SetHit(c,n,a,data);
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void CherenkovLikeHit::Clear()
{
  for(int i=0;i<4;i++){
    adchit[i].Clear();
    tdchit[i].Clear();
  }
  Seg = -1;
  NSensor=2;
  //  X = Y = Z = dX = dY = dZ = -999;
  Length = Width = Thick = -999;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool CherenkovLikeHit::CheckRange( int num_required )
{
  bool status = true;
  
  for(int i=0;i<4;i++){
    if( tdchit[i].CheckRange() ){}
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
    double size[4];    
    if( conf->GetGeomMapManager()->GetSize( CounterID, Seg, size) ){
      //      X = x; Y = y; Z = z; dX = dx; dY = dy; dZ = dz;
      Length = size[1]; Width = size[0]; Thick = size[2];
    }
  }

  for(int i=0;i<4;i++){
    adchit[i].Calc(conf); 
    tdchit[i].Calc(conf); 
  }  
  return true;
}
