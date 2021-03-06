#ifndef CDSTrackingMan_h
#define CDSTrackingMan_h 1

#include <vector>
#include <map>
#include <iostream>
#include <time.h>

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
#include "HitCluster.h"
#include "CDSTrack.h"
#include "BeamLineTrackMan.h"
#include "ConfMan.h"
#include "CDSFittingParamMan.h"
#include "MathTools.h"


class TrackVertex : public TObject
{
 private:
  int trackid[2];
  TVector3 Vertex_mean;
  TVector3 Vertex1, Vertex2;
  TVector3 Momentum1, Momentum2;

 public:
  TrackVertex();
  virtual ~TrackVertex() {};
  void Clear();  
 public:
  //Get Parameters
  int GetTrackID1(){return trackid[0];}
  int GetTrackID2(){return trackid[1];}
  TVector3 GetVertexPos(const int &aid);
  TVector3 GetVertexPos1(){return Vertex1;}
  TVector3 GetVertexPos2(){return Vertex2;}
  TVector3 GetVertexPos_mean(){return Vertex_mean;}
  TVector3 GetMomentum(const int &aid);
  TVector3 GetMomentum1(){return Momentum1;}
  TVector3 GetMomentum2(){return Momentum2;}

  //Set Parameters
  void SetTrackID1(const int &atrackid){trackid[0]=atrackid;}
  void SetTrackID2(const int &atrackid){trackid[1]=atrackid;}
  void SetVertexPos1(const TVector3 &avertex){Vertex1=avertex;}
  void SetVertexPos2(const TVector3 &avertex){Vertex2=avertex;}
  void SetVertexPos_mean(const TVector3 &avertex){Vertex_mean=avertex;}
  void SetMomentum1(const TVector3 &amomentum){Momentum1=amomentum;}
  void SetMomentum2(const TVector3 &amomentum){Momentum2=amomentum;}
  //Calclations

  ClassDef(TrackVertex, 1 );
};



class CDSTrackingMan : public TObject
{
 private:
  typedef std::vector <CDSTrack> TrackContainer;
  TrackContainer trackContainer;
  typedef std::vector <HitCluster> HitClusterContainer;
  HitClusterContainer ClusterContainer[7];
  int KilledLayer;
  double selected_chi;
  std::vector <int> GoodTrack;  
  typedef std::vector <TrackVertex> VertexContainer;
  VertexContainer Vertex;
  VertexContainer Vertex_beam;

 public:
  CDSTrackingMan();
  virtual ~CDSTrackingMan() {};
  void Clear();  
 public:
  //Get Parameters
  int GetKilledLayer(){return KilledLayer;}
  double SelectedChi(){return selected_chi;}//##Default 50
  bool GetVertex(const int &trackid1,const int &trackid2,TrackVertex &vertex);
  bool GetVertex_beam(const int &trackid,const int &beamid,TrackVertex &vertex);
  int nTrack(){return (int)trackContainer.size();}
  int nGoodTrack(){return (int)GoodTrack.size();}
  int GoodTrackID(const int &i){return GoodTrack[i];}
  int nCluster( const int &slayer ) { return (0<slayer&&slayer<=7) ?  (int)ClusterContainer[slayer-1].size() : 0; }
  HitCluster *Cluster(const int &slayer,const int &i){return (0<slayer&&slayer<=7) ? &ClusterContainer[slayer-1][i] : 0 ;}
  CDSTrack *Track(const int &i){return &trackContainer[i];}

  //Set Parameters
  void SetKilledLayer(const int &akilledlayer){KilledLayer=akilledlayer;}
  void SetSelectedChi(const double &achi){selected_chi=achi;}
  void AddTrack(const CDSTrack &atrack );

  //Calclations
  void Calc(ConfMan *conf);
  void CalcVertex(const int &id1, const int &id2);
  void CalcVertex_beam(const int &trackid,BeamLineTrackMan *bltrackman,ConfMan *conf);
  TrackVertex CalcVertex_beam(const int &trackid,LocalTrack *bltrack,ConfMan *conf, const int &bi);
  void CalcBetaMass(TVector3 vertex,LocalTrack *beam,CDSTrack* cdctrack,
		    ConfMan *conf, const int &beam_pid, const double &tof,
		    double &beta, double &mass);
  bool Execute(CDSHitMan *cdsMan,ConfMan *conf);
  void FullCircleFit();
  void FullHelixFit();
  void FullAxalLineFit();
  void FullStereoLineFit();
  bool SearchCluster(CDSHitMan *cdsMan,ConfMan *conf);
  bool FindingTrackCandidate(CDSHitMan *cdsMan,ConfMan *conf);
  bool FindingStereoHit(CDSHitMan *cdsMan,ConfMan *conf);
  void SelectSharedHit();
  bool SearchCDHHit(CDSHitMan *cdsMan);
  void PointToCircle(const double &x,const double &y,const double &radius,const double &x_cen,const double &y_cen,double &dis,double &xest,double &yest);

  ClassDef(CDSTrackingMan, 1 );
};



#endif
