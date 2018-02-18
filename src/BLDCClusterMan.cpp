// BeamLineTrackMan.cpp
#include <map>

#include "GlobalVariables.h"
#include "BLDCClusterMan.h"
#include "TVector2.h"
#include "TMath.h"

static const double SpatialResolutionOfBLDC=0.02; // [cm]

ClassImp(BLDCClusterMan);
ClassImp(BLDCCluster);

#define DEBUG 0
// ----------------------------- //
// class BLDCCluster             //
// ------------------------------//

BLDCCluster::BLDCCluster():
ClusterID(-1),
TimeMean(-999.),
TimeSub(-999.),
CTime(-999.)
{

}


BLDCCluster::~BLDCCluster()
{
  //  this->Clear();
}
void BLDCCluster::Calc( ConfMan* conf )
{
  if(nhit()==2){
    int cid=hit(0)->cid();
    int lay=hit(0)->layer();
    int wire=hit(0)->wire();
    int lay2=hit(1)->layer();
    if(lay!=lay2){
      TimeMean=(hit(0)->dt()+hit(1)->dt())/2.;
      TimeSub=(hit(0)->dt()-hit(1)->dt())/2.;
      if( lay>lay2 ){
	TimeSub *= -1;
	lay=lay2;
	wire=hit(1)->wire();
      }
      if(conf->GetDCTimeCorrManager()){
	CTime=conf->GetDCTimeCorrManager()->CalcCValue(cid,lay,0,TimeMean,TimeSub);
      }
    }else{
      TimeMean=-999;
      TimeSub=-999;
      CTime=-999;
    }
  }
  else{
    TimeMean=-999;
    TimeSub=-999;
    CTime=-999;
  }
   
}

void BLDCCluster::Clear() { 
  TimeMean=-999;
  TimeSub=-999;
  CTime=-999;
  ClusterID=-1;
  bldcHitCluster.clear(); 
}

// ----------------------------- //
// class BLDCClusterMan          //
// ------------------------------//

BLDCClusterMan::~BLDCClusterMan()
{
  //  this->Clear();
}

void BLDCClusterMan::Calc( ConfMan* conf )
{
  for(int icl=0;icl<4;icl++)
    for( unsigned int i=0;i<bldcClusterContainer[icl].size(); i++ ){
      bldcClusterContainer[icl][i].Calc(conf);
    }
}

void BLDCClusterMan::DeleteCluster( const int &ud, const int &xy, const int &i )
{
  BLDCClusterContainer::iterator it=bldcClusterContainer[ud+2*xy].begin();
  for( int j=0; j<i; j++ ) it++;
  bldcClusterContainer[ud+2*xy].erase( it );
}


void BLDCClusterMan::Clear()
{
  for(int i=0;i<4;i++)
    bldcClusterContainer[i].clear();
}

