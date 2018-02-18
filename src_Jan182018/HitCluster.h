#ifndef HitCluster_h
#define HitCluster_h 1

#include <vector>
#include <iostream>

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#include "TVector3.h"
#include "ChamberLikeHit.h"
#include "HodoscopeLikeHit.h"
#include "CDCHit.h"
#include "CDSHitMan.h"
#include "GlobalVariables.h"
#include "ConfMan.h"
#include "MathTools.h"


class HitCluster : public TObject
{
 private:
  typedef std::vector <CDCHit> ClusterHitContainer;
  ClusterHitContainer HitInClusterContainer;

  double MeanPhi,MeanRadius;
  double MeanX,MeanY;
  double MeanXp,MeanYp;

 public:
  HitCluster();
  virtual ~HitCluster() {};


 public:
  void Clear();
  double Meanphi()  const { return MeanPhi; }
  double Meanradius()  const { return MeanRadius; }
  double Meanx()  const { return MeanX; }
  double Meany()  const { return MeanY;} 
  double Meanxp()  const { return MeanXp; }
  double Meanyp()  const { return MeanYp;} 
  CDCHit *CDC( const int  &i ) { return &HitInClusterContainer[i]; }
  void Calc();
  void SetMeanphi( const double &meanphi ) { MeanPhi = meanphi; }
  void SetMeanradius( const double &meanradius ) { MeanRadius = meanradius; }
  int nHit() { return  (int)HitInClusterContainer.size(); }
  void SetHit( const CDCHit &hit ) { HitInClusterContainer.push_back(hit); }


  ClassDef(HitCluster, 1 );
};

#endif
