/* CircleFit.h */

// ======================================================== //
// Reference for circle(helix) fit :                        //   
//   http://www-jlc.kek.jp/subg/offl/lib/docs/helix_manip/  //
// ======================================================== //

#ifndef CircleFit_h
#define CircleFit_h 1

#include "TROOT.h"
#include "TVector3.h"
#include "TMinuit.h"

const int MAX_NUM_OF_HITS = 60;

class CircleFit : public TObject
{
  // ------------------------------------------- //
  // Parameters for circle fit :                 //
  //   par[3] = { d_rho, phi_0, alpha/k(=rho) }  //
  // ------------------------------------------- //
 private:
  Double_t Par[3];
  Double_t Err[3];
  
  Int_t NumOfHits;

  TVector3 HitPos[MAX_NUM_OF_HITS];
  Double_t Weight[MAX_NUM_OF_HITS];
  
  Double_t FitChi2;
  Int_t FitDof;
  Int_t FitStat;

 private:
  TMinuit *minuit;

 public:
  CircleFit();
  CircleFit(const double *initPar, const double *hitx, const double *hity,
	      const double *weight, const int &numofhit);
  ~CircleFit();

 public:
  void SetParameters( const double *param );
  void SetParameter( const int &i, const double &param){ Par[i]=(Double_t)param; }
  void SetHitPos( const double *hitx, const double *hity, 
		  const int &numofhit );
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
  void CalcChi2();

  ClassDef(CircleFit,1);
};

#endif 
