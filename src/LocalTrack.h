/* LocalTrack.h */

#ifndef LocalTrack_h
#define LocalTrack_h 1

#include <iostream>
#include <vector>

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#include "ChamberLikeHit.h"
#include "ConfMan.h"
#include "BLDCFittingParamMan.h"
#include "BLDCClusterMan.h"

class LocalTrack : public TObject
{
 private:
  typedef std::vector<BLDCCluster> BLDCClusterContainer;
  BLDCClusterContainer clusterContainer[2];
  //  int AssociateHodoClusterID;

 public:
  LocalTrack();
  //  LocalTrack( const LocalTrack &right );
  ~LocalTrack();

  int SetCluster(const int &xy, BLDCCluster *clu );
  int ncluster( const int &xy ) const { return clusterContainer[xy].size(); }
  BLDCCluster* cluster( const int &xy, const int &i ){ return &(clusterContainer[xy][i]); }
  
  int nhit();
  int nhit( const int &xy );
  ChamberLikeHit *hit( const int &i );
  ChamberLikeHit *hit( const int &xy, const int &i );
  void DeleteHit( const int &xy, const int &i );


 private:
  double TrackTime;
  double TrackTimeRMS;
  double TrackTimeX,TrackTimeY;

 public:
  void CalcTrackTime();
  void DeleteOffTimingCluster();
  double GetTrackTime() const { return TrackTime; }
  double GetTrackTimeRMS() const { return TrackTimeRMS; }
  double GetTrackTime(const int &xy) const { return xy ? TrackTimeY: TrackTimeX; }
  
 private:
  double A, B, C;		/* a track in XZ plane : Ax + Bz + C = 0 */
  double D, E, F;		/* a track in YZ plane : Dy + Ez + F = 0 */
  double GA, GB, GC;		/* Parameters for a global track */
  double GD, GE, GF;		/* Parameters for a global track */

  int xzDof, yzDof;
  double xzChi, yzChi;

 public:  
  void SetABC( const double &a, const double &b, const double &c ){ A=a; B=b; C=c; }
  void SetDEF( const double &d, const double &e, const double &f ){ D=d; E=e; F=f; }
  void SetGParam( const double &a, const double &b, const double &c, const double &d );
  
  void SetChisqrXZ( const double &val ) { xzChi=val; }
  void SetChisqrYZ( const double &val ) { yzChi=val; }
  void SetDofXZ( const int &val ) { xzDof=val; }
  void SetDofYZ( const int &val ) { yzDof=val; }

  double a() const { return A; }
  double b() const { return B; }
  double c() const { return C; }
  double d() const { return D; }
  double e() const { return E; }
  double f() const { return F; }
  double x() const  { return -C/A; }
  double dx() const { return -B/A; }
  double y() const  { return -F/D; }
  double dy() const { return -E/D; }
  double ga() const { return GA; }
  double gb() const { return GB; }
  double gc() const { return GC; }
  double gd() const { return GD; }
  double ge() const { return GE; }
  double gf() const { return GF; }
  double gx() const  { return -GC/GA; }
  double gdx() const { return -GB/GA; }
  double gy() const  { return -GF/GD; }
  double gdy() const { return -GE/GD; }
  void labc( double &a, double &b, double &c ) { a=A; b=B-hit(0,0)->tilt()*Deg2Rad; c=C; }
  void ldef( double &d, double &e, double &f ) { d=D; e=E-hit(1,0)->tilt()*Deg2Rad; f=F; }
  void abc( double &a, double &b, double &c ) const { a=A; b=B; c=C; }
  void def( double &d, double &e, double &f ) const { d=D; e=E; f=F; }
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

 public:
  bool XYLocalPosatZ( const double &z, double &x, double &y, const bool &tilt=false, const bool &rot =false);
  bool XYPosatZ( const double &z, double &x, double &y );
  TVector3 GetPosatZ(const double &z);
  TVector3 GetMomDir();
  bool LeastSquareFit( ConfMan *conf, const int &xy, TString option="" ); /* XZ plane: 0, YZ plane: 1 */
  bool LinearFit( ConfMan *conf, const bool &CHECKLR=false );
  bool CompareCluster( LocalTrack *tr,const bool &HIT=true );
  bool CompareTrackHit( LocalTrack *tr );

 public:
  void ConvLocalToGlobal();
  void ConvLocalToGlobal2();
  bool Calc( ConfMan *conf , TString option="");
  void CalcHitPosition( const bool &TILT=false);
  void CalcResidual( const bool &ROT=false, const TString &option="");
  double GetCalcChisquare();
  bool CheckRange(const double &ll,const double &ul){
    return (TrackTime>ll&&TrackTime<ul) ;
  }
  void Clear();
  void Print();

  ClassDef(LocalTrack,1);
};
#endif
