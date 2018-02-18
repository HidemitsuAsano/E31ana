#ifndef CDSTrack_h
#define CDSTrackh 1

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
#include "MathTools.h"


class CDSTrack : public TObject
{
 private:
  typedef std::vector <CDCHit> TrackHitContainer;
  TrackHitContainer trackContainer[15];
  HodoscopeLikeHit CDHhit;
  double ChiSqr;
  double ChiSqrtmp[4];
  double dof;
  double Al,Bl,Cl;//ax+by+c=0
  double Ap,Bp,Cp,Dp;//ax+by+cz+d=0
  double theta,jitter;
  double Clusterx[7],Clustery[7];
  double param[5];
  double paramtmp[5][5];
  double pt;
  double mom;
  double CirCenterX,CirCenterY,CirRho;  
  double magneticfield;
  bool goodflag;
  bool CDHflag;
  int fittinglevel;
  int pid;
  TVector3 CDHvertex; 
  TVector3 Line_o,Line_d;
  double dfunc_PTH(const TVector3 &pos,const double &helixphi,const double *par);
  double dfunc_LTH(const TVector3 &lpos,const TVector3 &dline,const double &helixphi,const double *par);

 public:
  CDSTrack();
  //  CDSTrack(const CDSTracl &atrack);
  virtual ~CDSTrack() {};
  void Clear();  
 public:

  //Get parameters
  int FittingLevel() const { return fittinglevel; }
  int PID() const { return pid; }
  void SetPID(const int &apid) { pid=apid; }
  double Chi() const { return ChiSqr; }
  double TmpChi(const int &num) const {return (0<=num&&num<4) ? ChiSqrtmp[num] : 0; }
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
  double Momentum() const {return mom;}
  bool GoodFlag() { return goodflag; }
  double ClusterX( const int &slayer)  {return (0<slayer&&slayer<=7) ? Clusterx[slayer-1] : 0;}
  double ClusterY( const int &slayer)  {return (0<slayer&&slayer<=7) ? Clustery[slayer-1] : 0;}

  int nTrackHit( const int &layer ) { return (0<layer&&layer<=NumOfCDCLayers) ? trackContainer[layer-1].size() : 0; }
  CDCHit *TrackHit( const int &layer, const int &i ) { return (0<layer&&layer<=NumOfCDCLayers) ? &trackContainer[layer-1][i] : 0; }
  bool CDHFlag() { return CDHflag; }
  HodoscopeLikeHit *CDHHit(){ return &CDHhit;}
  TVector3 CDHVertex(){return CDHvertex;}
  TVector3 LineDirection(){return Line_d;}
  TVector3 LineOrigin(){return Line_o;}
  void GetParameters(double *aparam);

  void GetGParameters(double *aparam);

  void GetTmpParameters(const int &num, double *aparam);
  TVector3 GetPosition(const double &helixphi); 
  TVector3 GetPosition(const double &helixphi,const double *par); 
  bool GetMomentum(const TVector3 &pos,TVector3 &p); 

  //Set Parameters
  void SetParameters(const double *aparam);
  void SetChiSqr( const double &chi ) { ChiSqr = chi; }
  void SetDof( const double &adof ) { dof = adof; }
  void SetABC( const double &a,const double &b,const double &c ) { Al = a; Bl = b; Cl = c; }
  void SetTheta( const double &atheta ) { theta = atheta; }
  void SetJitter( const double &ajitter ) { jitter = ajitter; }
  void SetClusterX( const int &slayer,const double &clx ) {if(0<slayer&&slayer<=7)Clusterx[slayer-1]=clx;  return;  }
  void SetClusterY( const int &slayer,const double &cly ) {if(0<slayer&&slayer<=7) Clustery[slayer-1]=cly;  return;  }
  void SetGoodFlag( const bool &flag ) { goodflag = flag; }

  void AddHit(const CDCHit &hit );
  void SetCDHHit(const HodoscopeLikeHit &cdhhit );
  void SetCDHVertex(const TVector3 &avertex ){CDHvertex=avertex; }
  bool RemoveAllHitInLayer( const int &layer );
  bool DeleteHit( const int &layer,const int &i );

  //Calclations
  void Calc(ConfMan *conf);
  void  XYatZ( double &x, double &y, const double &z );
  void  XYtoZ( const  double &x, const double &y, double &z );
  void PointToCircle(const double &x,const double &y,const double &radius,
		     const double &x_cen,const double &y_cen,double &dis,
		     double &xest,double &yest);
  bool LineToCircle(const double &a,const double &b, const double &c, 
		    const double &rho ,const double &xc, const double &yc,
		    double &x_p,double &y_p,double &x_n,double &y_n );
  double CalcHelixPhi(const double &x,const double &y); 
  double CalcHelixPhi(const double &x,const double &y,const double *par); 
  bool LineToHelix(const TVector3 &a, const TVector3 &dline, 
		   const double *par, TVector3 &lnest,
		   TVector3 &hnest, double &dis);
  bool PointToHelix(const TVector3 &hitpos, const double *par,
		    TVector3 &fitpos,double &dis); 
  void PointToLine( const TVector3 &p,
		    const TVector3 &x, const TVector3 &a,
		    double &dist,TVector3 &xest );
  void LineToLine( const TVector3 &x1, const TVector3 &a1,
		   const TVector3 &x2, const TVector3 &a2,
		   const double &dl,
		   double &dist,
		   TVector3 &xest, TVector3 &next );

  void ReconstructHit(ConfMan *Conf);

  //Fitting routinue
  void SetHitPos();
  void Refit(ConfMan *conf);
  bool FirstCircleFitting();
  bool FirstHelixFitting();
  bool CircleFitting();
  bool HelixFitting();
  bool FirstLineFitting();
  bool LineFitting();
  bool FirstStereoLineFitting();
  bool StereoLineFitting();
  bool TestLineFitting();
  bool TestStereoLineFitting();
  void CalcLineChiSqr();
  bool CalcInitialParameters();
  bool DoLineFit( const int &n, const double *x, const double *y,
		  const double *w,	
		  double &a, double &b, double &c, double &chi );  

  ClassDef(CDSTrack, 1 );
};




#endif
