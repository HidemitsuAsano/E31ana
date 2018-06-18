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

#include "ConfMan.h"
#include "TKO.h"
#include "BeamLineHitMan.h"
#include "ScalerMan.h"
#include "HodoscopeLikeHit.h"
#include "CherenkovLikeHit.h"
#include "ChamberLikeHit.h"
#include "EventHeader.h"

int main( int argc, char **argv )
{
  TROOT root( "GUI", "GUI" );
  TApplication theApp( "App", &argc, argv );
  gROOT->SetStyle( "Plain" );
  gROOT->cd();

  /*** load library ***/
  gSystem->Load("lib/libAll.so");

  /*** assign input file & call tree ***/
  TFile *f = new TFile("root/730.root");
  TTree *evtree = (TTree*)f->Get("EventTree");
  TTree *scatree = (TTree*)f->Get("ScalerTree");

  /*** conf file for new parameters ***/
  ConfMan *conf = new ConfMan("conf/730.conf");
  conf->Initialize();

  /*** assign output file ***/
  //TFile *of = new TFile("root/out.root","recreate");
  TFile *of = new TFile("root/out730.root","recreate");


  /* --------------- */
  /* ---- start ---- */
  /* --------------- */

  /*** declaration of classes ***/
  //TKOHitCollection *tko = 0;
  BeamLineHitMan *blMan = 0;
  EventHeader *head = 0;
  ScalerMan *scaMan = 0;
  //evtree->SetBranchAddress( "TKOHitCol", &tko );
  evtree->SetBranchAddress( "BeamLineHitMan", &blMan );
  evtree->SetBranchAddress( "EventHeader",  &head );
  scatree->SetBranchAddress( "ScalerMan", &scaMan );
  
  of->cd();
  TH1F *h1;
  TH2F *h2;

  /*                 */
  /* scaler analysis */
  /*                 */
  new TH1F( "sca", "sca", 50, 0, 50 );
  int nev = scatree->GetEntries();
  std::cout << " # Scaler read : " << nev << std::endl;
  for( int iev=0; iev<nev; iev++ ){
    scatree->GetEvent(iev);
    for( int i=0; i<scaMan->nsca(); i++ ){
      int val = scaMan->sca(i)->val();
      TString name = scaMan->sca(i)->name();
      h1=(TH1F*)gFile->Get("sca"); h1->Fill(i,val);
    }
  }

  /*                */
  /* event analysis */
  /*                */
  new TH1F( "Pattern", "Pattern", 10, 0, 10 );
  for( int i=1; i<=32; i++ ){
    new TH1F( Form( "TOF%ddediv", i ), Form( "TOF%ddediv", i ), 200, 0, 10 );

    new TH1F( Form( "TOF%dadcu", i ), Form( "TOF%dadcu", i ), 4000, 0, 4000 );
    new TH1F( Form( "TOF%dadcd", i ), Form( "TOF%dadcd", i ), 4000, 0, 4000 );
    new TH1F( Form( "TOF%dtdcu", i ), Form( "TOF%dtdcu", i ), 4000, 0, 4000 );
    new TH1F( Form( "TOF%dtdcd", i ), Form( "TOF%dtdcd", i ), 4000, 0, 4000 );
    new TH1F( Form( "TOF%ddeu", i ), Form( "TOF%ddeu", i ), 300, -5, 25 );
    new TH1F( Form( "TOF%dded", i ), Form( "TOF%dded", i ), 300, -5, 25 );
    new TH1F( Form( "TOF%dtimeu", i ), Form( "TOF%dtimeu", i ), 400, -20, 20 );
    new TH1F( Form( "TOF%dtimed", i ), Form( "TOF%dtimed", i ), 400, -20, 20 );
    new TH2F( Form( "TOF%ddetimeu", i ), Form( "TOF%ddetimeu", i ), 100, -20, 20, 100, -5, 25 );
    new TH2F( Form( "TOF%ddetimed", i ), Form( "TOF%ddetimed", i ), 100, -20, 20, 100, -5, 25 );
  }
  for( int i=1; i<=20; i++ ){
    new TH1F( Form( "BHD%dtdcu", i ), Form( "BHD%dtdcu", i ), 4000, 0, 4000 );
    new TH1F( Form( "BHD%dtdcd", i ), Form( "BHD%dtdcd", i ), 4000, 0, 4000 );
    new TH1F( Form( "BHD%dtimeu", i ), Form( "BHD%dtimeu", i ), 400, -20, 20 );
    new TH1F( Form( "BHD%dtimed", i ), Form( "BHD%dtimed", i ), 400, -20, 20 );
    new TH1F( Form( "BHD%dtmean", i ), Form( "BHD%dtmean", i ), 400, -20, 20 );
    new TH1F( Form( "BHD%dtsub",  i ), Form( "BHD%dtsub",  i ), 400, -20, 20 );
  }
  for( int i=1; i<=5; i++ ){
    new TH1F( Form( "T0%dtdcu", i ), Form( "T0%dtdcu", i ), 4000, 0, 4000 );
    new TH1F( Form( "T0%dtdcd", i ), Form( "T0%dtdcd", i ), 4000, 0, 4000 );
    new TH1F( Form( "T0%dtimeu", i ), Form( "T0%dtimeu", i ), 400, -20, 20 );
    new TH1F( Form( "T0%dtimed", i ), Form( "T0%dtimed", i ), 400, -20, 20 );
    new TH1F( Form( "T0%dtmean", i ), Form( "T0%dtmean", i ), 400, -20, 20 );
    new TH1F( Form( "T0%dtsub",  i ), Form( "T0%dtsub",  i ), 400, -20, 20 );
  }
  new TH1F( "B1ADCL",  "B1ADCL", 1500, 0, 1500 );
  new TH1F( "B1ADCR",  "B1ADCR", 1500, 0, 1500 );

  new TH1F( "BHD_T0k",  "BHD-T0 Kaon Trig", 300, -10, 20 );
  new TH1F( "BHD_T0pi", "BHD-T0 Pion Trig", 300, -10, 20 );
  new TH1F( "BHD_T0p",  "BHD-T0 Proton Trig", 300, -10, 20 );
  new TH1F( "BHD_T0e",  "BHD-T0 Electron Trig", 300, -10, 20 );

  for( int lay=1; lay<=8; lay++ ){
    for( int wire=1; wire<=32; wire++ ){
      new TH1F( Form("PDC1_Time_L%d_W%d",lay,wire), Form("PDC1_Time_L%d_W%d",lay,wire), 300, -300, 1200 );
      new TH1F( Form("PDC2_Time_L%d_W%d",lay,wire), Form("PDC2_Time_L%d_W%d",lay,wire), 300, -300, 1200 );
      new TH1F( Form("BLC1_Time_L%d_W%d",lay,wire), Form("BLC1_Time_L%d_W%d",lay,wire), 300, -300, 1200 );
      new TH1F( Form("BLC2_Time_L%d_W%d",lay,wire), Form("BLC2_Time_L%d_W%d",lay,wire), 300, -300, 1200 );

      new TH1F( Form("PDC1_DL_L%d_W%d",lay,wire), Form("PDC1_DL_L%d_W%d",lay,wire), 200, -0.5, 1.5 );
      new TH1F( Form("PDC2_DL_L%d_W%d",lay,wire), Form("PDC2_DL_L%d_W%d",lay,wire), 200, -0.5, 1.5 );
      new TH1F( Form("BLC1_DL_L%d_W%d",lay,wire), Form("BLC1_DL_L%d_W%d",lay,wire), 200, -0.5, 1.5 );
      new TH1F( Form("BLC2_DL_L%d_W%d",lay,wire), Form("BLC2_DL_L%d_W%d",lay,wire), 200, -0.5, 1.5 );
    }
    new TH1F( Form("PDC1_HitPat_L%d",lay), Form("PDC1_HitPat_L%d",lay), 35, 0, 35 );
    new TH1F( Form("PDC2_HitPat_L%d",lay), Form("PDC2_HitPat_L%d",lay), 35, 0, 35 );
    new TH1F( Form("BLC1_HitPat_L%d",lay), Form("BLC1_HitPat_L%d",lay), 35, 0, 35 );
    new TH1F( Form("BLC2_HitPat_L%d",lay), Form("BLC2_HitPat_L%d",lay), 35, 0, 35 );
  }

  nev = evtree->GetEntries();
  std::cout << " #Event: " << nev << std::endl;
  //for( int iev=0; iev<nev; iev++ ){
  for( int iev=0; iev<500; iev++ ){
    evtree->GetEvent(iev);

    if( iev%5000==0 )
      std::cout << " Event : " << iev << std::endl;

    /*** if some parameters are newed, you can re-calc all quantities as below, ***/
    blMan->Calc(conf);

    {
      HodoscopeLikeHit *hit = blMan->Hodo(CID_TOFstop,10);
      if( hit!=0 )
	std::cout << hit << "  " << hit->seg() << std::endl;
    }


    /*** after here, we only fill histograms. ***/
    for( int i=0; i<20; i++ ){
      int val = head->pattern(i);
      if( 0<val ){ h1=(TH1F*)gFile->Get("Pattern"); h1->Fill(i); }
    }
    for( int i=0; i<blMan->nTOF(); i++ ){
      int seg  = blMan->TOF(i)->seg();
      //std::cout << " nTOF:" << blMan->nTOF() << " seg:" << seg << std::endl;
      int adcu = blMan->TOF(i)->adcu();
      int adcd = blMan->TOF(i)->adcd();
      h1=(TH1F*)gFile->Get( Form("TOF%dadcu",seg ) ); h1->Fill(adcu);
      h1=(TH1F*)gFile->Get( Form("TOF%dadcd",seg ) ); h1->Fill(adcd);
      int tdcu = blMan->TOF(i)->tdcu();
      int tdcd = blMan->TOF(i)->tdcd();
      if( 0<tdcu ){ h1=(TH1F*)gFile->Get( Form("TOF%dtdcu",seg ) ); h1->Fill(tdcu); }
      if( 0<tdcd ){ h1=(TH1F*)gFile->Get( Form("TOF%dtdcd",seg ) ); h1->Fill(tdcd); }

      double tu = blMan->TOF(i)->tu();
      double td = blMan->TOF(i)->td();
      double eu = blMan->TOF(i)->eu();
      double ed = blMan->TOF(i)->ed();
      if( 20<ed ){
	h1=(TH1F*)gFile->Get( Form("TOF%ddediv",seg ) ); h1->Fill(eu/ed);
      }
      if( tdcu==0 || tdcd==0 ) continue;
      h1=(TH1F*)gFile->Get( Form("TOF%dtimeu",seg ) ); h1->Fill(tu);
      h1=(TH1F*)gFile->Get( Form("TOF%dtimed",seg ) ); h1->Fill(td);
      h1=(TH1F*)gFile->Get( Form("TOF%ddeu",seg ) ); h1->Fill(eu);
      h1=(TH1F*)gFile->Get( Form("TOF%dded",seg ) ); h1->Fill(ed);
      h2=(TH2F*)gFile->Get( Form("TOF%ddetimeu",seg ) ); h1->Fill(tu,eu);
      h2=(TH2F*)gFile->Get( Form("TOF%ddetimed",seg ) ); h1->Fill(td,ed);
    }

    for( int i=0; i<blMan->nBHD(); i++ ){
      int seg  = blMan->BHD(i)->seg();
      int tdcu = blMan->BHD(i)->tdcu();
      int tdcd = blMan->BHD(i)->tdcd();
      if( 0<tdcu ){ h1=(TH1F*)gFile->Get( Form("BHD%dtdcu",seg ) ); h1->Fill(tdcu); }
      if( 0<tdcd ){ h1=(TH1F*)gFile->Get( Form("BHD%dtdcd",seg ) ); h1->Fill(tdcd); }

      double tu = blMan->BHD(i)->tu();
      double td = blMan->BHD(i)->td();
      if( tdcu==0 || tdcd==0 ) continue;
      double tmean = blMan->BHD(i)->tmean(); 
      double tsub  = blMan->BHD(i)->tsub(); 
      h1=(TH1F*)gFile->Get( Form("BHD%dtimeu",seg ) ); h1->Fill(tu);
      h1=(TH1F*)gFile->Get( Form("BHD%dtimed",seg ) ); h1->Fill(td);
      h1=(TH1F*)gFile->Get( Form("BHD%dtmean",seg ) ); h1->Fill(tmean);
      h1=(TH1F*)gFile->Get( Form("BHD%dtsub", seg ) ); h1->Fill(tsub);
    }

    for( int i=0; i<blMan->nT0(); i++ ){
      int seg  = blMan->T0(i)->seg();
      int tdcu = blMan->T0(i)->tdcu();
      int tdcd = blMan->T0(i)->tdcd();
      if( 0<tdcu ){ h1=(TH1F*)gFile->Get( Form("T0%dtdcu",seg ) ); h1->Fill(tdcu); }
      if( 0<tdcd ){ h1=(TH1F*)gFile->Get( Form("T0%dtdcd",seg ) ); h1->Fill(tdcd); }

      double tu = blMan->T0(i)->tu();
      double td = blMan->T0(i)->td(); 
      if( tdcu==0 || tdcd==0 ) continue;
      double tmean = blMan->T0(i)->tmean(); 
      double tsub  = blMan->T0(i)->tsub(); 
      h1=(TH1F*)gFile->Get( Form("T0%dtimeu",seg ) ); h1->Fill(tu);
      h1=(TH1F*)gFile->Get( Form("T0%dtimed",seg ) ); h1->Fill(td);
      h1=(TH1F*)gFile->Get( Form("T0%dtmean",seg ) ); h1->Fill(tmean);
      h1=(TH1F*)gFile->Get( Form("T0%dtsub", seg ) ); h1->Fill(tsub);
    }

    for( int it=0; it<blMan->nT0(); it++ ){
      HodoscopeLikeHit *t0 = blMan->T0(it);
      if( t0->tdcu()==0 || t0->tdcd()==0 ) continue;
      for( int ib=0; ib<blMan->nBHD(); ib++ ){
	HodoscopeLikeHit *bhd = blMan->BHD(ib);
	if( bhd->tdcu()==0 || bhd->tdcd()==0 ) continue;
	double mt0 = t0->tmean();
	double mbhd = bhd->tmean();
	if( head->kaon() ){     h1=(TH1F*)gFile->Get("BHD_T0k");  h1->Fill(mt0-mbhd); }
	if( head->pion() ){     h1=(TH1F*)gFile->Get("BHD_T0pi"); h1->Fill(mt0-mbhd); }
	if( head->proton() ){   h1=(TH1F*)gFile->Get("BHD_T0p");  h1->Fill(mt0-mbhd); }
	if( head->electron() ){ h1=(TH1F*)gFile->Get("BHD_T0e");  h1->Fill(mt0-mbhd); }
      }
    }

    for( int it=0; it<blMan->nB1(); it++ ){
      HodoscopeLikeHit *b1 = blMan->B1(it);
      int adcu = blMan->B1(it)->adcu();
      int adcd = blMan->B1(it)->adcd();
      h1=(TH1F*)gFile->Get( ("B1ADCL" ) ); h1->Fill(adcu);
      h1=(TH1F*)gFile->Get( ("B1ADCR" ) ); h1->Fill(adcd);
    }

    for( int lay=1; lay<=8; lay++ ){
      // PDC1
      int nhit = blMan->nPDC1(lay);
      for( int i=0; i<nhit; i++ ){
	ChamberLikeHit *hit = blMan->PDC1( lay, i );
	int wire = hit->wire();
	int tdc  = hit->tdc();
	double dt = hit->dt();
	double dl = hit->dl();
	h1=(TH1F*)gFile->Get( Form("PDC1_Time_L%d_W%d",lay,wire) ); h1->Fill(dt);
	h1=(TH1F*)gFile->Get( Form("PDC1_DL_L%d_W%d",lay,wire) ); h1->Fill(dl);
	h1=(TH1F*)gFile->Get( Form("PDC1_HitPat_L%d",lay) ); h1->Fill(wire);
      }
      // PDC2
      nhit = blMan->nPDC2(lay);
      for( int i=0; i<nhit; i++ ){
	ChamberLikeHit *hit = blMan->PDC2( lay, i );
	int wire = hit->wire();
	int tdc  = hit->tdc();
	double dt = hit->dt();
	double dl = hit->dl();
	h1=(TH1F*)gFile->Get( Form("PDC2_Time_L%d_W%d",lay,wire) ); h1->Fill(dt);
	h1=(TH1F*)gFile->Get( Form("PDC2_DL_L%d_W%d",lay,wire) ); h1->Fill(dl);
	h1=(TH1F*)gFile->Get( Form("PDC2_HitPat_L%d",lay) ); h1->Fill(wire);
      }
      // BLC1
      nhit = blMan->nBLC1(lay);
      for( int i=0; i<nhit; i++ ){
	ChamberLikeHit *hit = blMan->BLC1( lay, i );
	int wire = hit->wire();
	int tdc  = hit->tdc();
	double dt = hit->dt();
	double dl = hit->dl();
	h1=(TH1F*)gFile->Get( Form("BLC1_Time_L%d_W%d",lay,wire) ); h1->Fill(dt);
	h1=(TH1F*)gFile->Get( Form("BLC1_HitPat_L%d",lay) ); h1->Fill(wire);
	h1=(TH1F*)gFile->Get( Form("BLC1_DL_L%d_W%d",lay,wire) ); h1->Fill(dl);
      }
      // BLC2
      nhit = blMan->nBLC2(lay);
      for( int i=0; i<nhit; i++ ){
	ChamberLikeHit *hit = blMan->BLC2( lay, i );
	int wire = hit->wire();
	int tdc  = hit->tdc();
	double dt = hit->dt();
	double dl = hit->dl();
	h1=(TH1F*)gFile->Get( Form("BLC2_Time_L%d_W%d",lay,wire) ); h1->Fill(dt);
	h1=(TH1F*)gFile->Get( Form("BLC2_HitPat_L%d",lay) ); h1->Fill(wire);
	h1=(TH1F*)gFile->Get( Form("BLC2_DL_L%d_W%d",lay,wire) ); h1->Fill(dl);
      }
    }
  }

  gFile->Write();
  gFile->Close();

  gSystem->Exit(1);
  theApp.Run();

  return 0;
}
