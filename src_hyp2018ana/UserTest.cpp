#include <iostream>
#include <fstream>
#include <TApplication.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <TSystem.h>
#include <TBranch.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>

#include "ConfMan.h"
#include "TKO.h"
#include "BeamLineHitMan.h"
#include "ScalerMan.h"
#include "HodoscopeLikeHit.h"
#include "CherenkovLikeHit.h"
#include "ChamberLikeHit.h"
#include "EventHeader.h"

double fun1( double *x, double *par )
{
  
  return (*x)*par[0];
}

int main( int argc, char **argv )
{
  TROOT root( "GUI", "GUI" );
  TApplication theApp( "App", &argc, argv );
  gROOT->SetStyle( "Plain" );
  gROOT->cd();



  TF1 *func = new TF1( "fun1", fun1, 0, 1);

  TCanvas *c1 = new TCanvas( "c1", "c1", 600, 600 );

  //gSystem->Exit(1);
  theApp.Run();

  return 0;
}
