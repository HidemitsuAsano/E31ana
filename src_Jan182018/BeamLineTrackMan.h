/* BeamLineTrackMan.h */

#ifndef BeamLineTrackMan_h
#define BeamLineTrackMan_h 1

#include <iostream>
#include <vector>

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#include "BeamLineHitMan.h"
#include "ChamberLikeHit.h"
#include "LinearTrack.h"
#include "ConfMan.h"
#include "BLDCFittingParamMan.h"

class LocalTrack : public TObject
{
 private:
  double A, B, C;		/* a track in XZ plane : Ax + Bz + C = 0 */
  double D, E, F;		/* a track in YZ plane : Dy + Ez + F = 0 */
  double GA, GB, GC;		/* Parameters for a global track */
  double GD, GE, GF;		/* Parameters for a global track */

  int xzDof, yzDof;
  double xzChi, yzChi;
  double Vtx1, Vty1, Vtz1; 	/* 1 : upstream vertex point */
  double Vtx2, Vty2, Vtz2; 	/* 2 : downstream vertex point */

  bool LAYER[2][8];
  int CID[2];

 private:
  typedef std::vector<ChamberLikeHit> LocalTrackHitContainer;
  LocalTrackHitContainer xzLocalTrackHitContainer;
  LocalTrackHitContainer yzLocalTrackHitContainer;

 public:
  LocalTrack();
  LocalTrack( const LocalTrack &right );
  ~LocalTrack();

 public:
  void SetCID(const int &cid1,const int &cid2){ CID[0]=cid1,CID[1]=cid2; }
  void SetLayerStatus(const int &cid, const bool lay[8])
    { 
      if(cid==CID[0]) 
	for(int i=0;i<8;i++) 
	  LAYER[0][i]=lay[i];
      if(cid==CID[1])
	for(int i=0;i<8;i++) 
	  LAYER[1][i]=lay[i];
    }
  bool layerstatus(const int &cid,const int &layer){
    if(cid==CID[0]) return LAYER[0][layer-1];
    else if(cid==CID[1]) return LAYER[1][layer-1];
    else return true;
  }

  int GetCID( const int &i ) { return CID[i]; }
  void GetLayerStatus( const int &ic, bool *out ) { for(int i=0;i<8;i++) out[i]=LAYER[ic][i]; }
  void SetHitXZ( ChamberLikeHit hit ){ xzLocalTrackHitContainer.push_back(hit); }
  int nhitxz() const { return xzLocalTrackHitContainer.size(); }
  ChamberLikeHit *hitxz( const int &i ){ return &(xzLocalTrackHitContainer[i]); }
  void DeleteHitXZ( const int &i );

  void SetHitYZ( ChamberLikeHit hit ){ yzLocalTrackHitContainer.push_back(hit); }
  int nhityz() const { return yzLocalTrackHitContainer.size(); }
  ChamberLikeHit *hityz( const int &i ){ return &(yzLocalTrackHitContainer[i]); }
  void DeleteHitYZ( const int &i );

  int nhit() const { return (xzLocalTrackHitContainer.size() + yzLocalTrackHitContainer.size()); }
  ChamberLikeHit *hit( const int &i );

  void SetHit( const int &xy, ChamberLikeHit hit ){
    if(xy==0)  xzLocalTrackHitContainer.push_back(hit);
    else if(xy==1)  yzLocalTrackHitContainer.push_back(hit);
  }
  int nhit( const int &xy ) { return xy ? nhityz() : nhitxz(); }
  ChamberLikeHit *hit( const int &xy, const int &i ){ return xy ? &(yzLocalTrackHitContainer[i]) : &(xzLocalTrackHitContainer[i]); }
  void DeleteHit( const int &xy, const int &i );

  void Clear();

  void SetABC( const double &a, const double &b, const double &c ){ A=a; B=b; C=c; }
  void SetDEF( const double &d, const double &e, const double &f ){ D=d; E=e; F=f; }

  void SetChisqrXZ( const double &val ) { xzChi=val; }
  void SetChisqrYZ( const double &val ) { yzChi=val; }
  void SetDofXZ( const int &val ) { xzDof=val; }
  void SetDofYZ( const int &val ) { yzDof=val; }
  void SetVertex1( const double &vx, const double &vy, const double &vz ){ Vtx1=vx; Vty1=vy; Vtz1=vz; }
  void SetVertex2( const double &vx, const double &vy, const double &vz ){ Vtx2=vx; Vty2=vy; Vtz2=vz; }

  double a() const { return A; }
  double b() const { return B; }
  double c() const { return C; }
  double d() const { return D; }
  double e() const { return E; }
  double f() const { return F; }
  void abc( double &a, double &b, double &c ) const { a=A; b=B; c=C; }
  void def( double &d, double &e, double &f ) const { d=D; e=E; f=F; }
  void semiabcdef( double &a, double &b, double &c, double &d, double &e, double &f );
  void gabc( double &a, double &b, double &c ) const { a=GA; b=GB; c=GC; }
  void gdef( double &d, double &e, double &f ) const { d=GD; e=GE; f=GF; }
  double chi2xz() const { return xzChi; }
  double chi2yz() const { return yzChi; }
  double chi2all() const { return (xzChi*xzDof + yzChi*yzDof)/(double)(xzDof+yzDof); }
  int dofxz() const { return xzDof; }
  int dofyz() const { return yzDof; }

  void SetChisqr( const int &xy, const double &val ) { 
    if(xy==0) xzChi=val; if(xy==1) yzChi=val; }
  void SetDof( const int &xy, const int &val ) {
    if(xy==0) xzDof=val; if(xy==1) yzDof=val; }
  double chi2(const int &xy) const { return xy ? yzChi : xzChi; }
  int dof(const int &xy) const { return xy ? yzDof : xzDof ; }

  void vertex1( double &vx, double &vy, double &vz ) const { vx=Vtx1; vy=Vty1; vz=Vtz1; }
  void vertex2( double &vx, double &vy, double &vz ) const { vx=Vtx2; vy=Vty2; vz=Vtz2; }

 public:
  bool XYLocalPosatZ( const double &z, double &x, double &y );
  bool XYSemiLocalPosatZ( const double &z, double &x, double &y );
  bool XYPosatZ( const double &z, double &x, double &y );
  bool Calc( ConfMan *conf );
  void CalcHitPosition( const bool &TILT=false);
  void CalcResidual( const bool &ROT=false);
  bool LeastSquareFit( ConfMan *conf, const int &xy ); /* XZ plane: 0, YZ plane: 1 */
  bool LinearFit( ConfMan *conf, const int &xy ); /* XZ plane: 0, YZ plane: 1 */
  bool PreTracking( ConfMan* conf, const int &xy);  
 public:
  void ConvLocalToGlobal();

  ClassDef(LocalTrack,1);
};

class BLDCCluster : public TObject
{ 
 private:
  typedef std::vector<ChamberLikeHit> BLDCHitCluster;
  BLDCHitCluster bldcHitCluster;

 public:
  BLDCCluster(){}
  BLDCCluster( const BLDCCluster &right );
  ~BLDCCluster();

 public:
  void SetHit( ChamberLikeHit hit ){ bldcHitCluster.push_back( hit ); }
  int nhit() const { return bldcHitCluster.size(); }
  ChamberLikeHit *hit( const int &i ){ return &(bldcHitCluster[i]); }
  void Clear() { bldcHitCluster.clear(); }

 ClassDef(BLDCCluster,1);
};


class BLDCClusterMan : public TObject
{
 private:
  typedef std::vector<BLDCCluster> BLDCClusterContainer;
  BLDCClusterContainer bldcClusterContainer[4]; // UD, XY
  
 public:
  BLDCClusterMan(){}
  BLDCClusterMan( const BLDCClusterMan &right );
  ~BLDCClusterMan();
 
 public:
  void SetCluster( const int &ud, const int &xy, BLDCCluster cluster ){ bldcClusterContainer[ud+2*xy].push_back( cluster ); }
  void SetClusterUpXZ( BLDCCluster cluster ){ bldcClusterContainer[0].push_back( cluster ); }
  void SetClusterDownXZ( BLDCCluster cluster ){ bldcClusterContainer[1].push_back( cluster ); }
  void SetClusterUpYZ( BLDCCluster cluster ){ bldcClusterContainer[2].push_back( cluster ); }
  void SetClusterDownYZ( BLDCCluster cluster ){ bldcClusterContainer[3].push_back( cluster ); }
  int ncluster( const int &ud, const int &xy) const { return bldcClusterContainer[ud+2*xy].size(); }
  int nclusterUpXZ() const { return bldcClusterContainer[0].size(); }
  int nclusterDownXZ() const { return bldcClusterContainer[1].size(); }
  int nclusterUpYZ() const { return bldcClusterContainer[2].size(); }
  int nclusterDownYZ() const { return bldcClusterContainer[3].size(); }
  BLDCCluster *cluster( const int &ud, const int &xy, const int &i ){ return &(bldcClusterContainer[ud+2*xy][i]); }
  BLDCCluster *clusterUpXZ( const int &i ){ return &(bldcClusterContainer[0][i]); }
  BLDCCluster *clusterDownXZ( const int &i ){ return &(bldcClusterContainer[1][i]); }
  BLDCCluster *clusterUpYZ( const int &i ){ return &(bldcClusterContainer[2][i]); }
  BLDCCluster *clusterDownYZ( const int &i ){ return &(bldcClusterContainer[3][i]); }
  void DeleteCluster( const int &ud, const int &xy, const int &i );
  void DeleteClusterUpXZ( const int &i );
  void DeleteClusterDownXZ( const int &i );
  void DeleteClusterUpYZ( const int &i );
  void DeleteClusterDownYZ( const int &i );
  
  void Clear();  
  ClassDef(BLDCClusterMan,1);
};


class BeamLineTrackMan : public TObject
{
 private:
  typedef std::vector<LocalTrack> BeamLineTrackContainer;
  BeamLineTrackContainer BLC1TrackContainer;
  BeamLineTrackContainer BLC1aTrackContainer;
  BeamLineTrackContainer BLC1bTrackContainer;
  BeamLineTrackContainer BLC2TrackContainer;
  BeamLineTrackContainer BLC2aTrackContainer;
  BeamLineTrackContainer BLC2bTrackContainer;
  BeamLineTrackContainer BPCTrackContainer;
  typedef std::vector<LinearTrack> LinearTrackContainer;
  LinearTrackContainer BLC1LinearTrackContainer;
  LinearTrackContainer BLC1aLinearTrackContainer;
  LinearTrackContainer BLC1bLinearTrackContainer;
  LinearTrackContainer BLC2LinearTrackContainer;
  LinearTrackContainer BLC2aLinearTrackContainer;
  LinearTrackContainer BLC2bLinearTrackContainer;
  LinearTrackContainer BPCLinearTrackContainer;
  LinearTrackContainer FDCLinearTrackContainer;

  int STATUS[8];
  bool LAYER[2][8];

 public:
  BeamLineTrackMan();
  BeamLineTrackMan( const BeamLineTrackMan &right );
  ~BeamLineTrackMan();

 public:
  int status(const int &cid);
  void SetStatus(const int &cid, const int &sta);
  int ntrackBLDC(const int &cid);
  LocalTrack *trackBLDC(const int &cid, const unsigned int &i);
  int nltrackBLDC(const int &cid);
  LinearTrack *ltrackBLDC(const int &cid, const unsigned int &i);
  void SetLinearTrack( const int  &cid, LinearTrack track );

  void SetTrackBLC1(  LocalTrack track ){ BLC1TrackContainer.push_back(track); }
  void SetLinearTrackBLC1(  LinearTrack track ){ BLC1LinearTrackContainer.push_back(track); }
  void SetTrackBLC1a( LocalTrack track ){ BLC1aTrackContainer.push_back(track); }
  void SetTrackBLC1b( LocalTrack track ){ BLC1bTrackContainer.push_back(track); }
  void SetLinearTrackBLC1a( LinearTrack track ){ BLC1aLinearTrackContainer.push_back(track); }
  void SetLinearTrackBLC1b( LinearTrack track ){ BLC1bLinearTrackContainer.push_back(track); }

  void SetTrackBLC2(  LocalTrack track ){ BLC2TrackContainer.push_back(track); }
  void SetLinearTrackBLC2(  LinearTrack track ){ BLC2LinearTrackContainer.push_back(track); }
  void SetTrackBLC2a( LocalTrack track ){ BLC2aTrackContainer.push_back(track); }
  void SetTrackBLC2b( LocalTrack track ){ BLC2bTrackContainer.push_back(track); }
  void SetLinearTrackBLC2a( LinearTrack track ){ BLC2aLinearTrackContainer.push_back(track); }
  void SetLinearTrackBLC2b( LinearTrack track ){ BLC2bLinearTrackContainer.push_back(track); }

  void SetTrackBPC(  LocalTrack track ){ BPCTrackContainer.push_back(track); }
  void SetLinearTrackBPC(  LinearTrack track ){ BPCLinearTrackContainer.push_back(track); }

  int ntrackBLC1()  const { return BLC1TrackContainer.size(); }
  int nltrackBLC1() const { return BLC1LinearTrackContainer.size(); }
  int ntrackBLC1a() const { return BLC1aTrackContainer.size(); }
  int ntrackBLC1b() const { return BLC1bTrackContainer.size(); }
  int nltrackBLC1a() const { return BLC1aLinearTrackContainer.size(); }
  int nltrackBLC1b() const { return BLC1bLinearTrackContainer.size(); }

  int ntrackBLC2()  const { return BLC2TrackContainer.size(); }
  int nltrackBLC2() const { return BLC2LinearTrackContainer.size(); }
  int ntrackBLC2a() const { return BLC2aTrackContainer.size(); }
  int ntrackBLC2b() const { return BLC2bTrackContainer.size(); }
  int nltrackBLC2a() const { return BLC2aLinearTrackContainer.size(); }
  int nltrackBLC2b() const { return BLC2bLinearTrackContainer.size(); }

  int ntrackBPC()  const { return BPCTrackContainer.size(); }
  int nltrackBPC()  const { return BPCLinearTrackContainer.size(); }

  LocalTrack  *trackBLC1( const int &i ){ return &(BLC1TrackContainer[i]); }
  LinearTrack *ltrackBLC1( const int &i ){ return &(BLC1LinearTrackContainer[i]); }
  LocalTrack  *trackBLC1a(const int &i ){ return &(BLC1aTrackContainer[i]); }
  LocalTrack  *trackBLC1b(const int &i ){ return &(BLC1bTrackContainer[i]); }
  LinearTrack *ltrackBLC1a(const int &i ){ return &(BLC1aLinearTrackContainer[i]); }
  LinearTrack *ltrackBLC1b(const int &i ){ return &(BLC1bLinearTrackContainer[i]); }

  LocalTrack  *trackBLC2( const int &i ){ return &(BLC2TrackContainer[i]); }
  LinearTrack *ltrackBLC2( const int &i ){ return &(BLC2LinearTrackContainer[i]); }
  LocalTrack  *trackBLC2a(const int &i ){ return &(BLC2aTrackContainer[i]); }
  LocalTrack  *trackBLC2b(const int &i ){ return &(BLC2bTrackContainer[i]); }
  LinearTrack *ltrackBLC2a(const int &i ){ return &(BLC2aLinearTrackContainer[i]); }
  LinearTrack *ltrackBLC2b(const int &i ){ return &(BLC2bLinearTrackContainer[i]); }

  LocalTrack *trackBPC( const int &i ){ return &(BPCTrackContainer[i]); }
  LinearTrack *ltrackBPC( const int &i ){ return &(BPCLinearTrackContainer[i]); }

  void DeleteTrackBLC1(  const int &i );
  void DeleteTrackBLC1a( const int &i );
  void DeleteTrackBLC1b( const int &i );
  void DeleteTrackBLC2(  const int &i );
  void DeleteTrackBLC2a( const int &i );
  void DeleteTrackBLC2b( const int &i );
  void DeleteTrackBPC(  const int &i );

  void Clear();
  void ClearBLC1(){  BLC1TrackContainer.clear(); }
  void ClearBLC1a(){ BLC1aTrackContainer.clear(); }
  void ClearBLC1b(){ BLC1bTrackContainer.clear(); }
  void ClearBLC2(){  BLC1TrackContainer.clear(); }
  void ClearBLC2a(){ BLC1aTrackContainer.clear(); }
  void ClearBLC2b(){ BLC1bTrackContainer.clear(); }
  void ClearBPC(){  BPCTrackContainer.clear(); }

 public:
  bool DoTracking( BeamLineHitMan *blMan, ConfMan *conf );
  
  // private:
  int LocalTracking( BeamLineHitMan *blMan, ConfMan *conf, const int &id );
  int LinearTracking( ConfMan *conf, const int &id );
  int LinearTracking( BeamLineHitMan *blMan, ConfMan *conf, const int &id );
  int LocalLinearTracking( ConfMan *conf, const int &id );
  int SemiLocalTracking( BeamLineHitMan *blMan, ConfMan *conf, const int &id );
  bool ConvertLocalToLinear(const int &cid);

 private:
  bool Clustering( BeamLineHitMan *blMan, BLDCClusterMan *clMan, 
		   ConfMan *conf, const int &id, const int &id2=0 );
  bool Clustering2(BeamLineHitMan *blMan, BLDCClusterMan *clMan, 
		   ConfMan *conf, const int &id );
  bool LeastSquareFit( LocalTrack *track, ConfMan *conf,
		       const int &xy ); /* XZ plane: 0, YZ plane: 1 */
  void ConvLocalToGlobal(){};
  void ConvLocalToGlobal(const int &id);

  ClassDef(BeamLineTrackMan,1);
};

#endif
