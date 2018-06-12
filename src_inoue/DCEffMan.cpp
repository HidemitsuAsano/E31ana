#include "DCEffMan.h"

DCEffMan::DCEffMan(const int &cid, ConfMan *confMan, BeamLineHitMan *blMan, BeamLineTrackMan *bltrackMan) 
  : fConfMan(confMan), fCID(cid), fBLMan(blMan), fBLTrackMan(bltrackMan)
{
  if( fCID==CID_BLC1a || fCID==CID_BLC1b || fCID==CID_BLC2a || fCID==CID_BLC2b || fCID==CID_BPC ){
    fHits.resize(8);
  }
  else if( fCID==CID_BLC1 || fCID==CID_BLC2 ){
    fHits.resize(16);
  }
  else if( fCID==CID_FDC1 ){
    fHits.resize(6);
  }
}

void DCEffMan::get()
{
  clear();
  for( int i=0; i<fHits.size(); i++ ){
    for( int j=0; j<fBLMan->nBLDC(fCID, i+1); j++ ){
      fHits[i].push_back(fBLMan->BLDC(fCID, i+1, j));
    }
  }
}

bool DCEffMan::trig18()
{
  bool flag1=false, flag2=false;
  int max_wire = maxWire();
  for( int j=0; j<fHits[0].size(); j++ ){
    if( fHits[0][j]->wire()>1 && fHits[0][j]->wire()<max_wire ) flag1=true;
  }
  
  int last = fHits.size()-1;
  for( int j=0; j<fHits[last].size(); j++ ){
    if( fHits[last][j]->wire()>1 && fHits[last][j]->wire()<max_wire ) flag2=true;
  }
  if( flag1 && flag2 ) return true;
  else return false;
}

bool DCEffMan::trig18_2()
{
  bool flag1=false, flag2=false;
  int max_wire = maxWire();
  for( int j=0; j<fHits[0].size(); j++ ){
    if( fHits[0][j]->wire()>2 && fHits[0][j]->wire()<max_wire-1 ) flag1=true;
  }

  int last = fHits.size()-1;
  for( int j=0; j<fHits[last].size(); j++ ){
    if( fHits[last][j]->wire()>2 && fHits[last][j]->wire()<max_wire-1 ) flag2=true;
  }
  if( flag1 && flag2 ) return true;
  else return false;
}

bool DCEffMan::trig(const int &lay, const int &max) const
{
  for( int i=0; i<fHits.size(); i++ ){
    if( i==lay-1 ) continue;
    if( fHits[i].size()<1 || fHits[i].size()>max ) return false;
  }
  return true;
}

bool DCEffMan::eff(const int &lay, const int &max) const
{
  if( fHits[lay-1].size()>max )    return false;
  else if( fHits[lay-1].size()>0 ) return true;
  else                             return false;
}

LocalTrack* DCEffMan::getTrack()
{
  if( fBLTrackMan->ntrackBLDC(fCID)==1 ){
    return fBLTrackMan->trackBLDC(fCID, 0);
  }
  else return 0;
}

bool DCEffMan::trigTrack(const int &lay)
{



}

bool DCEffMan::all1hit()
{
  for( int i=0; i<fHits.size(); i++ ){
    if( fHits[i].size()!=1 ) return false;
  }
    return true;
}

void DCEffMan::dump()
{
  std::cout<<"DCEffMan::dump CID="<<fCID<<std::endl;
  int max_wire = maxWire();
  bool flag[fHits.size()][max_wire];

  for( int i=0; i<fHits.size(); i++ ){
    for( int j=0; j<max_wire; j++ ){
      flag[i][j] = false;
    }
  }

  for( int i=0; i<fHits.size(); i++ ){
    for( int j=0; j<fHits[i].size(); j++ ){
      flag[i][fHits[i][j]->wire()-1] = true;
    }
  }

  for( int i=0; i<fHits.size(); i++ ){
    std::string str;
    for( int j=0; j<max_wire; j++ ){
      if( flag[i][j] ) str+= "*";
      else str +="-";
    }
    std::cout<<str<<std::endl;
  }
  char in;
  std::cin>>in;
  if( in=='q' ) exit(-1);
}

int DCEffMan::maxWire() const
{
  int max_wire=0;
  if( fCID==CID_BLC1 || fCID==CID_BLC2 || fCID==CID_BLC1a || fCID==CID_BLC1b ||  fCID==CID_BLC2a || fCID==CID_BLC2b ) max_wire=32;
  else if( fCID==CID_BPC ) max_wire=15;
  else if( fCID==CID_FDC1 ) max_wire=64;

  return max_wire;
}

void DCEffMan::clear()
{
  for( int i=0; i<fHits.size(); i++ ){
    fHits[i].clear();
  }
}
