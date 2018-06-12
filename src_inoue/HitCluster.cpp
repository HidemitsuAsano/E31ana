#include "HitCluster.h"

ClassImp(HitCluster);


//const double Inclradius=20.40,Outclradius=47.78,cellsize=0.845*2.;
//const double ResolutionOfCDC=0.02;

HitCluster::HitCluster(): TObject()
{
  Clear();
}

// HitCluster::HitCluster(CDSHitMan *cdsman): TObject()
// {
//   Clear();
//   cdsMan=cdsman;
// }

void HitCluster::Clear()
{
  HitInClusterContainer.clear();
}

void HitCluster::Calc(CDSHitMan* cdsMan)
{
  MeanPhi=0; MeanRadius=0;
  MeanX=0; MeanY=0;
  MeanXp=0; MeanYp=0;
  for(int n=0;n<(int)nHit();n++)
    {
      CDCHit *cdc=CDC(cdsMan,n);
      MeanRadius+=cdc->radius();
      MeanX+=cdc->wx();
      MeanY+=cdc->wy();
      MeanXp+=cdc->wxp();
      MeanYp+=cdc->wyp();
    }
  MeanX=MeanX/(double)nHit();
  MeanY=MeanY/(double)nHit();
  MeanXp=MeanXp/(double)nHit();
  MeanYp=MeanYp/(double)nHit();
  MeanRadius=MeanRadius/(double)nHit();
  MeanPhi=MathTools::CalcDeg(MeanX,MeanY);
}

bool HitCluster::DeleteHit( const int &i) 
{
  ClusterHitIDContainer::iterator it = HitInClusterContainer.begin();
  for( int j=0; j<i; j++ ) ++it;
  HitInClusterContainer.erase( it ) ;
  return true;
}

