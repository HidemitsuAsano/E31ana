#include "HitCluster.h"

ClassImp(HitCluster);


//const double Inclradius=20.40,Outclradius=47.78,cellsize=0.845*2.;
//const double ResolutionOfCDC=0.02;

HitCluster::HitCluster(): TObject()
{
  Clear();
}

void HitCluster::Clear()
{
  HitInClusterContainer.clear();
}

void HitCluster::Calc()
{
  MeanPhi=0; MeanRadius=0;
  MeanX=0; MeanY=0;
  MeanXp=0; MeanYp=0;
  for(int n=0;n<(int)HitInClusterContainer.size();n++)
    {
      MeanRadius+=HitInClusterContainer[n].radius();
      MeanX+=HitInClusterContainer[n].wx();
      MeanY+=HitInClusterContainer[n].wy();
      MeanXp+=HitInClusterContainer[n].wxp();
      MeanYp+=HitInClusterContainer[n].wyp();

    }
  
  MeanX=MeanX/(double)HitInClusterContainer.size();
  MeanY=MeanY/(double)HitInClusterContainer.size();
  MeanXp=MeanXp/(double)HitInClusterContainer.size();
  MeanYp=MeanYp/(double)HitInClusterContainer.size();
  MeanRadius=MeanRadius/(double)HitInClusterContainer.size();
  MathTools *tool=new MathTools();
  MeanPhi=tool->CalcDeg(MeanX,MeanY);
  delete tool;
}

