/* HelixFit.h */

// ======================================================== //
// Reference for circle(helix) fit :                        //   
//   http://www-jlc.kek.jp/subg/offl/lib/docs/helix_manip/  //
// ======================================================== //

#ifndef HelixFit_h
#define HelixFit_h 1

#include "TROOT.h"
#include "TVector3.h"
#include "TMinuit.h"
#include "MathTools.h"
const int MAX_NUM_OF_HITS_H = 60;

class HelixFit : public TObject
{
  // ------------------------------------------- //
  // Parameters for circle fit :                 //
  //   par[3] = { d_rho, phi_0, alpha/k(=rho) }  //
  // ------------------------------------------- //
 private:
  Double_t Par[5];
  Double_t Err[5];
  
  Int_t NumOfHits;

  TVector3 HitPos[MAX_NUM_OF_HITS_H];
  Double_t Weight[MAX_NUM_OF_HITS_H];
  
  Double_t FitChi2;
  Int_t FitDof;
  Int_t FitStat;

 private:
  TMinuit *minuit;

 public:
  HelixFit();
  HelixFit(const double *initPar, const double *hitx,
	   const double *hity,const double *hitz,
	   const double *weight, const int &numofhit);
  ~HelixFit();

 public:
  void SetParameters( const double *param );
  void SetParameter( const int &i, const double &param){ Par[i]=(Double_t)param; }
  void SetHitPos( const TVector3 hitpos, const int &numofhit );
  void SetWeight( const double *weight, const int &numofhit );
  void SetNumOfHit( const int &i ){ NumOfHits=(Int_t)i; }

  void GetParameters( double *param );
  double param( const int &i ){ return (double)Par[i]; }
  double chisquare() { return (double)FitChi2; }
  int dof() { return (int)FitDof; }
  int stat() { return (int)FitStat; }

 public:
  void SetGlobalVariables();
  void fit();
  void CalcNearestPos( const Double_t &hitx, const Double_t &hity,
		       Double_t &fitx, Double_t &fity);
  void CalcNearestPosForStereo( const TVector3 &hitpos,TVector3 &fitpos);
  void CalcChi2();

  ClassDef(HelixFit,1);
};

#endif 
