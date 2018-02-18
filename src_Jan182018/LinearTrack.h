/* LinearTrack.h */

#ifndef LinearTrack_h
#define LinearTrack_h 1

#include <iostream>
#include <vector>

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#include "BeamLineHitMan.h"
#include "ChamberLikeHit.h"
#include "ConfMan.h"
#include "BLDCFittingParamMan.h"

class LinearTrack : public TObject
{
 private:
  double A, B;		/* a track in XZ plane : x = Az + B */
  double C, D;		/* a track in YZ plane : y = Cz + D */
  double GA, GB;		/* Parameters for a global track */
  double GC, GD;		/* Parameters for a global track */

  int Dof;
  double Chi;

  int CID[2];
  bool LAYER[2][8];

 private:
  typedef std::vector<ChamberLikeHit> vLinearTrackHitContainer;
  vLinearTrackHitContainer LinearTrackHitContainer;

 public:
  LinearTrack();
  LinearTrack( const LinearTrack &right );
  ~LinearTrack();

 public:
  void SetCID(const int &cid1,const int &cid2=0){ CID[0]=cid1,CID[1]=cid2; }
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
  void SetHit( ChamberLikeHit hit ){ LinearTrackHitContainer.push_back(hit); }
  int nhit() const { return LinearTrackHitContainer.size(); }
  ChamberLikeHit *hit( const int &i ){ return &(LinearTrackHitContainer[i]); }
  void DeleteHit( const int &i );

  void Clear();

  void SetAB( const double &a, const double &b){ A=a; B=b; }
  void SetCD( const double &c, const double &d){ C=c; D=d; }
  void SetABCD( const double &a, const double &b, const double &c, const double &d ){ A=a; B=b; C=c; D=d; }
  void SetGABCD( const double &a, const double &b, const double &c, const double &d );
  void SetChisqr( const double &val ) { Chi=val; }
  void SetDof( const int &val ) { Dof=val; }

  double a() const { return A; }
  void abcd( double &a, double &b, double &c, double &d) const { a=A; b=B; c=C; d=D; }
  void gabcd( double &a, double &b, double &c, double &d ) const { a=GA; b=GB; c=GC; d=GD; }
  double chi2all() const { return Chi; }
  int dof() const { return Dof; }

 public:
  bool XYLocalPosatZ( const double &z, double &x, double &y, const bool &ROT=false );
  bool XYPosatZ( const double &z, double &x, double &y );
  bool Calc( ConfMan *conf );
  void CalcHitPosition();
  void CalcResidual(bool PRINT=false);
  double GetCalcChisquare();
  bool LinearFit( ConfMan *conf, const bool &CHECKLR=true ); /* XZ plane: 0, YZ plane: 1 */
  
 public:
  void ConvLocalToGlobal();

  ClassDef(LinearTrack,1);
};


#endif
