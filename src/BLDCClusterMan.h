/* BLDCClusterMan.h */

#ifndef BLDCClusterMan_h
#define BLDCClusterMan_h 1

#include <iostream>
#include <vector>

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#include "ChamberLikeHit.h"
#include "ConfMan.h"

class BLDCCluster : public TObject
{ 
 private:
  typedef std::vector<ChamberLikeHit> BLDCHitCluster;
  BLDCHitCluster bldcHitCluster;

 public:
  BLDCCluster();
  //  BLDCCluster( const BLDCCluster &right );
  ~BLDCCluster();

 private:
  int ClusterID;
  double TimeMean;
  double TimeSub;
  double CTime;

 public:
  void Calc( ConfMan *conf );
 
  int GetClusterID() const { return ClusterID; }
  double GetTimeMean() const { return TimeMean; }
  double GetTimeSub() const { return TimeSub; }
  double GetCTime() const { return CTime; }
 
  void SetClusterID( const int &i ){ ClusterID = i; }
  void SetHit( ChamberLikeHit hit ){ bldcHitCluster.push_back( hit ); }
  int nhit() const { return bldcHitCluster.size(); }
  ChamberLikeHit *hit( const int &i ){ return &(bldcHitCluster[i]); }
  void Clear();
  
  ClassDef(BLDCCluster,1);
};


class BLDCClusterMan : public TObject
{
 private:
  typedef std::vector<BLDCCluster> BLDCClusterContainer;
  BLDCClusterContainer bldcClusterContainer[4]; // UD, XY
  
 public:
  BLDCClusterMan(){}
  //  BLDCClusterMan( const BLDCClusterMan &right );
  ~BLDCClusterMan();
 
 public:
  void SetCluster( const int &ud, const int &xy, BLDCCluster cluster ){ bldcClusterContainer[ud+2*xy].push_back( cluster ); }
  void SetCluster( const int &i, BLDCCluster cluster ){ bldcClusterContainer[i].push_back( cluster ); }

  int ncluster( const int &icl) const { return bldcClusterContainer[icl].size(); }
  int ncluster( const int &ud, const int &xy) const { return bldcClusterContainer[ud+2*xy].size(); }

  BLDCCluster *cluster( const int &icl, const int &i ){ return &(bldcClusterContainer[icl][i]); }
  BLDCCluster *cluster( const int &ud, const int &xy, const int &i ){ return &(bldcClusterContainer[ud+2*xy][i]); }

  void DeleteCluster( const int &ud, const int &xy, const int &i );
  
  void Calc( ConfMan *conf );
  void Clear();  
  ClassDef(BLDCClusterMan,1);
};


#endif
