// CircleFit.cpp
#include <iostream>
#include "CircleFit.h"
#include "TMath.h"
#include "MathTools.h"
#include <cmath>

ClassImp(CircleFit);


// parameters for TMinuit
static const Double_t  FitStep[3] = { 0.001, 0.001, 1e-8 };
static const Double_t LowLimit[3] = { -50, -2*TMath::Pi(), -0.5 };
static const Double_t  UpLimit[3] = { 50, 2*TMath::Pi(), 0.5 };
// static const Double_t  FitStep[3] = { 0.001, 0.0001, 0.001 };
// static const Double_t LowLimit[3] = { -4000, -TMath::Pi(), -2000 };
// static const Double_t  UpLimit[3] = { 4000, TMath::Pi(), 2000 };

// global variables for TMinuit
static Int_t gNumOfHits;
static TVector3 gHitPos[MAX_NUM_OF_HITS];
static Double_t gWeight[MAX_NUM_OF_HITS];


// --------------------------------------------------------------//
// functions for TMinuit


static void fcn( Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag )
{
  Double_t chisq=0.;
  Int_t dof = 0;
  Double_t hitx, hity, fitx, fity, dis;
  for( Int_t i=0; i<gNumOfHits; i++ ){    
    hitx = gHitPos[i].x();
    hity = gHitPos[i].y();
    MathTools::PointToCircle(par, hitx, hity, fitx, fity, dis);
    chisq +=pow( dis/gWeight[i] , 2 );
    dof++;
  }

  f = chisq/(dof-3);
}

// --------------------------------------------------------------//
CircleFit::CircleFit()
{
  for( int i=0; i<3; i++ ){
    Par[i] = Err[i] = -999.;
  }

  NumOfHits = FitDof = FitStat = -999;
  FitChi2 = -999.;

  for( int i=0; i<MAX_NUM_OF_HITS; i++ ){
    HitPos[i].SetXYZ(-999,-999,-999);
    Weight[i] = -999.;
  }
  
  minuit = new TMinuit(3);
  TROOT minexam("CircleFit","Circle fit using TMinuit");
}

CircleFit::CircleFit( const double *initPar, const double *hitx,
			  const double *hity, const double *weight,
			  const int &numofhit)
{
  NumOfHits = (Int_t)numofhit;

  for( int i=0; i<3; i++ ) {
    Par[i] = (Double_t)initPar[i];
    Err[i] = -999.;
  }

  for( int i=0; i<MAX_NUM_OF_HITS; i++ ) {
    if( i<NumOfHits ){
      HitPos[i].SetXYZ( (Double_t)hitx[i], (Double_t)hity[i], -999. );
      Weight[i] = (Double_t)weight[i];
    }
    else {
      HitPos[i].SetXYZ(-999,-999,-999);
      Weight[i] = -999.;
    }
  }

  FitDof = FitStat = -999;
  FitChi2 = -999.;

  minuit = new TMinuit(3);
  TROOT minexam("CircleFit","Circle fit using TMinuit");

  fit();
}

CircleFit::~CircleFit()
{
  delete minuit;
}

void CircleFit::SetParameters( const double *param )
{
  for( int i=0; i<3; i++ ) Par[i] = (Double_t)param[i];
}

void CircleFit::SetHitPos( const double *hitx, const double *hity,
			     const int &numofhit )
{
  for( int i=0; i<numofhit; i++ )
          HitPos[i].SetXYZ( (Double_t)hitx[i], (Double_t)hity[i], 0. );
}

void CircleFit::SetWeight( const double *weight, const int &numofhit )
{
  for( int i=0; i<numofhit; i++ )
    Weight[i] = (Double_t)weight[i];
}

void CircleFit::GetParameters( double *param )
{
  for( int i=0; i<3; i++ ) param[i] = (double)Par[i];
}

void CircleFit::SetGlobalVariables()
{
  gNumOfHits = NumOfHits;
  for( int i=0; i<MAX_NUM_OF_HITS; i++ ){
    gHitPos[i] = HitPos[i];
    gWeight[i] = Weight[i];
  }
}

void CircleFit::fit()
{
#if 0
  std::cout << "!!! CircleFit::fit() !!!" << std::endl;
#endif

#if DEBUG
  Int_t plevel=1;
#else
  Int_t plevel=-1;
#endif

  SetGlobalVariables();

  minuit->SetPrintLevel( plevel );
  minuit->SetFCN( fcn );
  
  Double_t arglist[10];
  Int_t ierflg = 0;

  arglist[0] = 1;
  //  minuit->mnexcm("SET ERR", arglist,1,ierflg);
  minuit->mnexcm("SET NOW", arglist,1,ierflg); // No warnings

  // Set starting values and step sizes for parameters
  TString name[3] = {"rho", "xc", "yc"};
  for( Int_t i=0; i<3; i++){
    minuit->mnparm(i, name[i],Par[i],FitStep[i],LowLimit[i],UpLimit[i],ierflg);
  }

  // Now ready for minimization step
  arglist[0] = 1000;
  arglist[1] = 1.;
  minuit->mnexcm("MIGRAD", arglist, 2, ierflg);

  //minuit->Command("TMProve 100");

  // Print results
  Double_t amin, edm, errdef;
  Int_t nvpar, nparx, icstat;
  minuit->mnstat(amin, edm,  errdef, nvpar, nparx, icstat);
#if DEBUG
  minuit->mnprin(3,amin);
#endif
  //  std::cout<<minuit->fNfcn<<std::endl;  
  Int_t err;
  Double_t bnd1, bnd2;
  for( Int_t i=0; i<3; i++ ){
    minuit->mnpout(i, name[i], Par[i], Err[i], bnd1, bnd2, err);
  }
  
  FitStat = icstat;
  CalcChi2();
}

void CircleFit::CalcChi2()
{  
  Double_t chisq=0.;
  Int_t dof = 0;
  Double_t hitx, hity, fitx, fity, dis;
  for( Int_t i=0; i<NumOfHits; i++ ){
    hitx = HitPos[i].x();
    hity = HitPos[i].y();
    MathTools::PointToCircle( Par, hitx, hity, fitx, fity ,dis);
    chisq +=pow( dis/ Weight[i] , 2 );
    dof++;
  }

  FitDof = dof - 3;
  FitChi2 = chisq;
}


