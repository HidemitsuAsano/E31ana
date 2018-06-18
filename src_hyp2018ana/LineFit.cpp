// Linefit2.cpp
#include <iostream>
#include "LineFit.h"
#include "TMath.h"

#define DEBUG 0
ClassImp(LineFit);

// parameters for TMinuit
static const int npar = 6;
static const Double_t  FitStep[npar] = { 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10 };
static const Double_t LowLimit[npar] = { -100, -100, -100, -100 , -100, -100 };
static const Double_t  UpLimit[npar] = {  100,  100,  100,  100 ,  100,  100 };

// global variables for TMinuit
static Int_t gNumOfHits;
static TVector3 gWirePos[MAX_NUM_OF_HITS_L];
static TVector3 gWireDir[MAX_NUM_OF_HITS_L];
static Double_t gWeight[MAX_NUM_OF_HITS_L];
static Double_t gDrift[MAX_NUM_OF_HITS_L];

// --------------------------------------------------------------//
// functions for TMinuit
static void fcn2( Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag )
{
  Double_t chisq=0.;
  Int_t dof = 0;
  Double_t dist=0.;
  TVector3 fittmp1,fittmp2;
  TVector3 pos(par[0],par[1],par[2]);  
  TVector3 dir(par[3],par[4],par[5]);  
  for( Int_t i=0; i<gNumOfHits; i++ ){
    {
      MathTools::LineToLine(gWirePos[i], gWireDir[i],pos,dir,0,dist,fittmp1,fittmp2);
      chisq +=pow( ( dist - gDrift[i] )/gWeight[i] , 2 );
    }
    dof++;
  }  
  f = chisq/(dof-4);
}

// --------------------------------------------------------------//

LineFit::LineFit()
{
  for( int i=0; i<4; i++ ){
    Par[i] = Err[i] = -999.;
  }

  NumOfHits = FitDof = FitStat = -999;
  FitChi2 = -999.;

  for( int i=0; i<MAX_NUM_OF_HITS_L; i++ ){
    WirePos[i].SetXYZ(-999,-999,-999);
    WireDir[i].SetXYZ(-999,-999,-999);
    Weight[i] = -999.;
    Drift[i] = -999.;
  }
  
  minuit = new TMinuit(5);
  TROOT minexam("LineFit","line fit using TMinuit");
}

LineFit::LineFit( const double *initPar, const TVector3 *wirepos,
		    const TVector3 *wiredir,const double *weight,
		    const double *drift, const int &numofhit)
{
  NumOfHits = (Int_t)numofhit;

  for( int i=0; i<npar; i++ ) {
    Par[i] = (Double_t)initPar[i];
    Err[i] = -999.;
  }

  for( int i=0; i<MAX_NUM_OF_HITS_L; i++ ) {
    if( i<NumOfHits ){
      WirePos[i]= wirepos[i];
      WireDir[i]= wiredir[i];
      Weight[i] = (Double_t)weight[i];      
      Drift[i]  = (Double_t)drift[i];      
    }
    else {
      WirePos[i].SetXYZ(-999,-999,-999);
      WireDir[i].SetXYZ(-999,-999,-999);
      Weight[i] = -999.;
      Drift[i]  = -999.;
    }
  }

  FitDof = FitStat = -999;
  FitChi2 = -999.;

  minuit = new TMinuit(npar);
  TROOT minexam("LineFit","Line fit using TMinuit");

  fit();
}

LineFit::~LineFit()
{
  delete minuit;
}
void LineFit::SetParameters( const double *param )
{
  for( int i=0; i<npar; i++ ) Par[i] = (Double_t)param[i];
}
void LineFit::SetWirePos( const TVector3 hitpos,  const int &numofhit )
{
  for( int i=0; i<numofhit; i++ )    WirePos[i]=hitpos;
}
void LineFit::SetWireDir( const TVector3 hitpos,  const int &numofhit )
{
  for( int i=0; i<numofhit; i++ )    WireDir[i]=hitpos;
}
void LineFit::SetWeight( const double *weight, const int &numofhit )
{
  for( int i=0; i<numofhit; i++ )    Weight[i] = (Double_t)weight[i];
}
void LineFit::SetDrift( const double *weight, const int &numofhit )
{
  for( int i=0; i<numofhit; i++ )    Drift[i] = (Double_t)weight[i];
}
void LineFit::GetParameters( double *param )
{
  for( int i=0; i<npar; i++ ) param[i] = (double)Par[i];
}

void LineFit::SetGlobalVariables()
{
  gNumOfHits = NumOfHits;
  for( int i=0; i<gNumOfHits; i++ ){
    gWirePos[i] = WirePos[i];
    gWireDir[i] = WireDir[i];
    gWeight[i] = Weight[i];
    gDrift[i] = Drift[i];
  }
}

void LineFit::fit()
{
#if 0
  std::cout << "!!! LineFit::fit() !!!" << std::endl;
#endif

#if DEBUG
  Int_t plevel=1;
#else
  Int_t plevel=-1;
#endif

  SetGlobalVariables();

  minuit->SetPrintLevel( plevel );
  minuit->SetFCN( fcn2 );

  Double_t arglist[10];
  Int_t ierflg = 0;

  arglist[0] = 1;
  //  minuit->mnexcm("SET ERR", arglist,1,ierflg);
  minuit->mnexcm("SET NOW", arglist,1,ierflg); // No warnings

  // Set starting values and step sizes for parameters
  TString name[npar] = {"x0", "y0", "z0","xdir","ydir","zdir"};
  for( Int_t i=0; i<npar; i++){
    minuit->mnparm(i, name[i],Par[i],FitStep[i],LowLimit[i],UpLimit[i],ierflg);
  }
  minuit->Command("SET STRategy 0");
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
  minuit->mnprin(npar,amin);
#endif
  
  Int_t err;
  Double_t bnd1, bnd2;
  for( Int_t i=0; i<npar; i++ ){
    minuit->mnpout(i, name[i], Par[i], Err[i], bnd1, bnd2, err);
  }
#if DEBUG
  for( Int_t i=0; i<npar; i++ )
    std::cout<<Par[i]<<"  ";
  std::cout<<std::endl;
#endif
  FitStat = icstat;
  CalcChi2();

}

void LineFit::CalcChi2()
{  
  Double_t chisq=0.;
  Int_t dof = 0;
  Double_t dist=0.;
  TVector3 fittmp1,fittmp2;
  TVector3 pos(Par[0],Par[1],Par[2]);  
  TVector3 dir(Par[3],Par[4],Par[5]);  
  for( Int_t i=0; i<NumOfHits; i++ ){
    {
      MathTools::LineToLine(WirePos[i], WireDir[i],pos,dir,0,dist,fittmp1,fittmp2);
      chisq +=pow( ( dist - Drift[i] )/Weight[i] , 2 );
    }
    dof++;
  }  

  FitDof = dof - 4;
  FitChi2 = chisq;
}

