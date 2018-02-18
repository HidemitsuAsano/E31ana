/* NCClusterMan.h */

#ifndef NCClusterMan_h
#define NCClusterMan_h 1

#include <iostream>
#include <vector>

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#include "HodoscopeLikeHit.h"
#include "BeamLineHitMan.h"
#include "CDSHitMan.h"
#include "HitMan.h"
#include "ConfMan.h"
#include "GlobalVariables.h"

#include "TVector3.h"

class HodoHit : public TObject
{
 private:
  int Seg;
  double Time;
  double Edep;
  TVector3 Pos;

 public:
  HodoHit(){ }
  HodoHit( int seg, double time, double edep, const TVector3 &tmppos);
  ~HodoHit(){}

  bool operator<(const HodoHit &rhs) const { return this->seg() < rhs.seg(); }
  int seg() const { return Seg; }
  double time() const { return Time; }
  double edep() const { return Edep; }
  TVector3 pos() const { return Pos; }

  ClassDef(HodoHit,1);
};

class HodoCluster : public TObject
{ 
 private:
  std::vector<HodoHit> hitContainer;

 public:
  HodoCluster(){ Clear(); }
  HodoCluster( const int &cid,const int &cluid );
  ~HodoCluster(){}

 private:
  int CounterID;
  int ClusterID;
  double CTimeMean;
  double EdepTotal;
  int aTrackID;
  bool FLAG;
  
 public:
  void Calc( );
  
  int GetCounterID() const { return CounterID; }
  int GetClusterID() const { return ClusterID; }
  double GetCTimeMean() const { return CTimeMean; }
  double GetEdepTotal() const { return EdepTotal; }
  
  void SetCounterID( const int &i ){ CounterID = i; }
  void SetClusterID( const int &i ){ ClusterID = i; }
  void SetHit( const HodoscopeLikeHit &hit );
  
  int nhit() const { return hitContainer.size(); }
  int nhit( const double &threshold );
  
  bool first(  double threshold , int &seg, double &time, TVector3 &pos);
  bool firsttime(  double threshold , int &seg, double &time, TVector3 &pos);
  bool firsttime(  double threshold , int &seg, double &time, TVector3 &pos,double &edep);
  bool firsttimeinlayer(  double threshold , int &seg, double &time, TVector3 &pos);
  bool second( double threshold , int &seg, double &time, TVector3 &pos);
  
  //  HodoscopeLikeHit *first( HitMan *hitMan,const double &threshold );
  //  HodoscopeLikeHit *second( HitMan *hitMan,const double &threshold );
  
  HodoHit *hit( int i ){ return (i>=0&&i<nhit()) ?  &hitContainer[i]: 0;}
  HodoscopeLikeHit *hit( HitMan *hitMan,const int &i )
  { return (i>=0&&i<nhit()) ? hitMan->Hodo(CounterID,hitContainer[i].seg()) : 0; }
  
  bool IsFlag() { return FLAG; }
  bool IsNeutral() { return FLAG; }
  void Clear();
  
  ClassDef(HodoCluster,1);
};

class HodoClusterMan : public TObject
{
 private:
  typedef std::vector<HodoCluster> ClusterContainer;
  typedef std::map<int,ClusterContainer> HodoClusterContainer;
  HodoClusterContainer Container;
  
  typedef std::vector<int> NCCVCIDContainer;
  NCCVCIDContainer nccvcID;

 public:
  HodoClusterMan(){};
  HodoClusterMan( BeamLineHitMan *blMan, CDSHitMan *cdsMan, ConfMan *conf );
  ~HodoClusterMan(){};
 
 public:
  int SetCluster(const int &cid, const HodoCluster &cluster ){ Container[cid].push_back( cluster ); return ncluster(cid)-1; }

  int ncluster( int cid) { return Container[cid].size(); }
  int ncluster( int cid, HitMan *hitMan, double threshold);
  HodoCluster *ncfirst( double threshold , int &seg, double &time, TVector3 &pos);
  HodoCluster *ncfirsttime( double threshold , int &seg, double &time, TVector3 &pos);
  HodoCluster *ncfirsttime( double threshold , int &seg, double &time, TVector3 &pos, double &edep);
  bool Clustering( BeamLineHitMan *blMan , ConfMan *conf );
  bool Clustering( CDSHitMan *blMan , ConfMan *conf );
  bool Clustering( HitMan *hitMan , ConfMan *conf, const int &cid );

  HodoCluster *cluster(const int &cid, const int &i) { return (i>=0&&i<ncluster(cid)) ? &(Container[cid][i]) : 0; }  
  void DeleteCluster(const int &cid, const int &i);

  int GetCVCID(const int &ncid){ return (ncid>=0&&ncid<(int)nccvcID.size()) ? nccvcID[ncid] : -1; };
  HodoCluster *GetCVC(const int &ncid){ return (ncid>=0&&ncid<(int)nccvcID.size()) ? &(Container[CID_CVC][GetCVCID(ncid)]) : 0; }  
  
  //  void Calc( HitMan *hitMan );
  void Calc();
  void Clear();  
  ClassDef(HodoClusterMan,1);
};


#endif
