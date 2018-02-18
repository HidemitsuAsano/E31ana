// MTDCHit.cpp

#include "MTDCHit.h"

ClassImp( aMTDCHit );
ClassImp( MTDCHit );
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
aMTDCHit::aMTDCHit() : TObject(),
		       TDC(-1),Time(-999.),CTime(-999)
{
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
aMTDCHit::aMTDCHit( const int &data ):
  TDC(data),Time(-999.),CTime(-999)
{
}
aMTDCHit::aMTDCHit( const int &data, const double &time , const double &ctime ) : TObject(),
										  TDC(data),Time(time),CTime(ctime)
{
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
bool aMTDCHit::CheckRange(const int &ll, const int &ul)
{
  if(ll<TDC && TDC<ul ) return true;
  else return false;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
MTDCHit::MTDCHit() : TObject(),
		     Crate(-1),Slot(-1),Channel(-1),
		     CounterID(-1),HitID(-1),Seg(-1)		    
{
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
MTDCHit::MTDCHit( const int &cr, const int &sl, const int &ch, const int &cid, const int &seg, const int &ud, const int &data ) : TObject()
{
  CounterID= cid; Seg = seg;
  Crate = cr; Slot = sl; Channel = ch;
  Data[ud].push_back( aMTDCHit( data ) );
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
void MTDCHit::Clear()
{
  Crate = Slot = Channel  = -1;
  CounterID=HitID=Seg=-1;
  Data[0].clear();
  Data[1].clear();
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
int MTDCHit::CheckRange(const int &ll, const int &ul)
{
  int n=0;
  int ud=0;
  for(int i=0;i<(int)Data[ud].size();i++)
    if( Data[ud][i].CheckRange(ll,ul))
      n++;
  return n;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void MTDCHit::Calc( ConfMan *conf, const int &ref )
{
  //  for(int i;i<Data.size();i++)
  //    Data  
}
