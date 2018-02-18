#ifndef CDSTrack_h
#define CDSTrack_h 1

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
#include "CDSFittingParamMan.h"
#include "CircleFit.h"
#include "HelixFit.h"
#include "LineFit.h"
#include "MathTools.h"
#include "ELossTools.h"

#define TMPPARAM 0
class CDSTrack : public TObject
{
 private:
  typedef std::vector <short> TrackHitIDContainer;
  TrackHitIDContainer hitidContainer[15];
  TrackHitIDContainer CDHhitID;
  TrackHitIDContainer IHhitID;
  float ChiSqr;
  float dof;
  float Al,Bl,Cl;//ax+by+c=0
  float Ap,Bp,Cp,Dp;//ax+by+cz+d=0
  float theta,jitter;
  float Clusterx[7],Clustery[7];
  double param[5];
  double pt;
  double CirCenterX,CirCenterY,CirRho;  
  double magneticfield;
  bool shortFlag;
  bool goodflag;
  bool CDHflag;
  bool IHflag;
  short fittinglevel;
  short pid;
  TVector3 CDHvertex; 
  TVector3 IHvertex; 
  TVector3 Line_o; // line origin ??
  TVector3 Line_d; // line direction ?
  typedef std::vector<double> parContainer;
  std::map<int,parContainer> parCont;
  std::map<int,TVector3> vtxContainer;

 public:
  CDSTrack();
  virtual ~CDSTrack() {};
  void Clear();  

#if TMPPARAM
 private:
  double ChiSqrtmp[4];
  double paramtmp[5][5];
 public:
  double TmpChi(const int &num) const {return (0<=num&&num<4) ? ChiSqrtmp[num] : 0; }
  void GetTmpParameters(const int &num, double *aparam);
#endif

 private:

  void CheckCharge();
  void CheckCharge2();

 public:
  void CalcELoss(double mass=-1);
  int FittingLevel() const { return fittinglevel; }
  int PID() const { return pid; }
  void SetPID(const int &apid) { pid=apid; }
  double Chi() const { return ChiSqr; }

  double Dof() const { return dof; }
  double A() const { return Al; }
  double B() const { return Bl; }
  double C() const { return Cl; }
  double Jitter() const {return jitter;}
  double Theta() const {return theta;}
  double CenterX() const {return CirCenterX;}
  double CenterY() const {return CirCenterY;}
  double Rho() const {return CirRho;}
  double Pt() const {return pt;}
  double Momentum() const {return pt*sqrt(1+param[4]*param[4]);}
  double Momentum(const int &id);
  double Mass(){ return cdsMass[PID()]; }
  bool GoodFlag() { return goodflag; }
  bool IsShort() { return shortFlag; }
  double ClusterX( const int &slayer)  {return (0<slayer&&slayer<=7) ? Clusterx[slayer-1] : 0;}
  double ClusterY( const int &slayer)  {return (0<slayer&&slayer<=7) ? Clustery[slayer-1] : 0;}
  int charge(){ return Pt() >0 ? 1 :-1; }
  int nTrackHit( const int &layer ) { return (0<layer&&layer<=NumOfCDCLayers) ? hitidContainer[layer-1].size() : 0; }
  CDCHit *TrackHit( CDSHitMan* cdsMan, const int &layer, const int &i )
  { return (0<layer&&layer<=NumOfCDCLayers) ? cdsMan->CDC(layer,(int)hitidContainer[layer-1][i]) : 0; }
  CDCHit *hit( CDSHitMan* cdsMan,const int &layer, const int &i ){ return TrackHit(cdsMan,layer,i); }
  bool CDHFlag() { return CDHflag; }
  bool IHFlag() { return IHflag; }

  int nCDHHit() const { return CDHhitID.size(); }
  HodoscopeLikeHit *CDHHit( CDSHitMan *cdsMan, const int &i){ return (i>=0&&i<(int)CDHhitID.size()) ? cdsMan->CDH(CDHhitID[i]) : 0 ;}
  TVector3 CDHVertex(){return CDHvertex;}

  bool GetCDHHit( CDSHitMan *cdsMan,int &seg, double &time);
  bool GetIHHit( CDSHitMan *cdsMan,int &seg, double &time);

  int nIHHit() const { return IHhitID.size(); }
  HodoscopeLikeHit *IHHit( CDSHitMan *cdsMan, const int &i){ return (i>=0&&i<(int)IHhitID.size()) ? cdsMan->IH(IHhitID[i]) : 0;}
  TVector3 IHVertex(){return IHvertex;}

  TVector3 LineDirection(){return Line_d;}
  TVector3 LineOrigin(){return Line_o;}

  int nParamSets(){ return (int)parCont.size(); }
  void AddParameters(const int &id, const double *par,const TVector3 &vtx);
  void GetParameters(double *aparam);
  bool GetGParameters(double *aparam);
  bool GetParameters(const int &id, double *aparam, TVector3 &vtx);
  bool GetNthParameters(int n, int &id, double *aparam, TVector3 &vtx);

  bool GetVertex(const TVector3 &pos, const TVector3 &dir, TVector3 &lpos, TVector3 &hpos);
  bool CalcVertexTimeLength(const TVector3 &pos,const TVector3 &dir,const double &mass,TVector3 &lpos, TVector3 &hpos,double &time, double &length, bool ADDPAR=false);
  //Set Parameters
  void SetParameters(const double *aparam);
  void SetChiSqr( const double &chi ) { ChiSqr = chi; }
  void SetDof( const double &adof ) { dof = adof; }
  void SetABC( const double &a,const double &b,const double &c ) { Al = a; Bl = b; Cl = c; }
  void SetTheta( const double &atheta ) { theta = atheta; }
  void SetJitter( const double &ajitter ) { jitter = ajitter; }
  void SetClusterX( const int &slayer,const double &clx ) {if(0<slayer&&slayer<=7) Clusterx[slayer-1]=clx;  return;  }
  void SetClusterY( const int &slayer,const double &cly ) {if(0<slayer&&slayer<=7) Clustery[slayer-1]=cly;  return;  }
  void SetGoodFlag( const bool &flag ) { goodflag = flag; }
  void SetShort( const bool &flag ) { goodflag = flag; }

  void AddHit(const CDCHit &hit );
  void SetHodoHit(const HodoscopeLikeHit &hit, const TVector3 &vertex );
  void SetCDHVertex(const TVector3 &avertex ){CDHvertex=avertex; }
  void SetIHVertex(const TVector3 &avertex ){IHvertex=avertex; }
  bool RemoveAllHitInLayer( const int &layer );
  bool DeleteHit( const int &layer,const int &i );

  //Calclations
  void Calc(ConfMan *conf);
  void XYatZ( double &x, double &y, const double &z );
  void XYtoZ( const  double &x, const double &y, double &z );
  bool GetMomentum(const TVector3 &pos,TVector3 &p, bool GLOBAL=true,bool ELOSS=false); 
  void ReconstructHit(CDSHitMan *cdsMan,ConfMan *Conf);
  bool SearchHodoHit2(CDSHitMan *cdsMan, ConfMan *conf);
  bool SearchHodoHit(CDSHitMan *cdsMan, ConfMan *conf,const double &cdhmaxdist=7.,const double &ihmaxdist=10);
  bool SearchHodoHit(CDSHitMan *cdsMan, ConfMan *conf,const gCounterID &cid,double maxdist=0, bool HITONLY=true);

  //Fitting routinu
  void SetHitPos(CDSHitMan *cdsMan,bool PRINT=false);

  void Refit(CDSHitMan *cdsMan,ConfMan *conf);
  bool Retiming(CDSHitMan *cdsMan,ConfMan *conf,double mass,bool SLEW=false);

  bool FirstCircleFitting(CDSHitMan *cdsMan);
  bool FirstHelixFitting(CDSHitMan *cdsMan);
  bool CircleFitting(CDSHitMan *cdsMan);
  bool CircleFitting2(CDSHitMan *cdsMan);
  bool HelixFitting(CDSHitMan *cdsMan);
  bool FirstLineFitting(CDSHitMan *cdsMan);
  bool LineFitting(CDSHitMan *cdsMan);
  bool FirstStereoLineFitting(CDSHitMan *cdsMan);
  bool StereoLineFitting(CDSHitMan *cdsMan);
  bool StereoLineFitting2(CDSHitMan *cdsMan);
  void CalcLineChiSqr(CDSHitMan *cdsMan);
  bool CalcInitialParameters();
  bool DoLineFit( const int &n, const double *x, const double *y,
		  const double *w,	
		  double &a, double &b, double &c, float &chi );  
  //Calc for short track
  bool CalcShortInitialParameters();
  void Print();

  ClassDef(CDSTrack, 1 );
};

#endif
