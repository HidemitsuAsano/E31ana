// TKO.cpp

// #include <string>
// #include <cstring>
// #include <iostream>
// #include <iomanip>

#include "TKO.h"

ClassImp( TKOHit );
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
TKOHit::TKOHit() : TObject(),
		   Crate(-1),Slot(-1),Channel(-1),Data(-1)
{
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
TKOHit::TKOHit( const int &cr, const int &sl, const int &ch, const int &data ) : TObject(),
										 Crate(cr),Slot(sl),Channel(ch),Data(data)
{
}

// // + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
// TKOHit::TKOHit( TKOHit *hit )
// {
//   Crate   = hit->cr();
//   Slot    = hit->sl();
//   Channel = hit->ch();
//   Data    = hit->data();
// }
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
void TKOHit::Clear()
{
  Crate = Slot = Channel = Data = -1;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool TKOHit::CheckRange(const int &ll, const int &ul) const
{
  if( ll<Data && Data<ul )
    return true;
  else
    return false;  
}


ClassImp( TKOHitCollection );

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
TKOHitCollection::TKOHitCollection() : TObject()
{
  Clear();
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
TKOHit * TKOHitCollection::hit( const int &i )
{
  unsigned int j = i;
  if( 0<=j && j<theCollection.size() ){
    return &theCollection[j];
  }
  else{
    return 0;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
void TKOHitCollection::Clear()
{
  theCollection.clear();
}

