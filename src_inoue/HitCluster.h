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
  //  CDSHitMan* cdsMan;
  typedef std::vector <int> ClusterHitIDContainer;
  ClusterHitIDContainer HitInClusterContainer;

  double MeanPhi,MeanRadius;
  double MeanX,MeanY;
  double MeanXp,MeanYp;

  int key(const int &layer, const int &hid){ return (hid<<5) | layer ; } 

 public:
  HitCluster();
  //  HitCluster( CDSHitMan* cdsman );
  virtual ~HitCluster() {};

 public:
  //  void SetCDSHitMan(CDSHitMan *cdsman){ cdsMan=cdsman; }

  void Clear();
  double Meanphi()  const { return MeanPhi; }
  double Meanradius()  const { return MeanRadius; }
  double Meanx()  const { return MeanX; }
  double Meany()  const { return MeanY;} 
  double Meanxp()  const { return MeanXp; }
  double Meanyp()  const { return MeanYp;} 
  inline CDCHit *CDC( CDSHitMan *cdsMan, const int  &i );

  void Calc(CDSHitMan *cdsMan);
  void SetMeanphi( const double &meanphi ) { MeanPhi = meanphi; }
  void SetMeanradius( const double &meanradius ) { MeanRadius = meanradius; }
  int nHit() { return  (int)HitInClusterContainer.size(); }
  void SetHit( const CDCHit &hit ) { HitInClusterContainer.push_back(key(hit.layer(),hit.hid())); }
  bool DeleteHit( const int &i );  
  
  ClassDef(HitCluster, 1 );
};

inline CDCHit* HitCluster::CDC(CDSHitMan *cdsMan, const int &i){
  if(!cdsMan) return 0;
  int layer= HitInClusterContainer[i] & 0x1F;
  int wire = (HitInClusterContainer[i] >> 5) & 0xFFF;
  return cdsMan->CDC(layer,wire); 
}
#endif
