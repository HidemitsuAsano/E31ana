#include <fstream.h>
#include <iostream.h>
#include <vector.h>
#include <math.h>
#include "TMinuit.h"
#include "TMath.h"

// class k18ana
class ConfMan;
class CDSHitMan;
class BeamLineHitMan;
class CDCHit;
class CDSTrackingMan;
class BeamLineTrackMan;
class Display;

// class knucl
class MCData;

#define DEBUG 1
#define DISPLAY 1

void mcread()
{

  /*** load library ***/
  gSystem->Load("./lib/libAll.so");
  gSystem->Load("../../geant/knucl3/libKnuclRootData.so");

  /*** assign input file & call tree ***/
  TFile *f = new TFile("simout.root");

  //#if SIM
  TTree *evtree = (TTree*)f->Get("SimTree");

  /*** conf file for new parameters ***/
  ConfMan *confMan = new ConfMan("./conf/Oct2010/analyzer3057-3060.conf");
  confMan->Initialize();

  /*** declaration of classes ***/
  CDSHitMan *cdsMan = 0;
  BeamLineHitMan* blMan = 0;
  MCData* mcData = 0;
  evtree->SetBranchAddress( "CDSHitMan", &cdsMan );
  evtree->SetBranchAddress( "BeamLineHitMan", &blMan );
  evtree->SetBranchAddress( "MCData", &mcData );

#if DISPLAY
  Display *disp = new Display();
  TCanvas *c_disp = new TCanvas( "c_disp", "c_disp", 0, 0, 1200, 1200 );
  c_disp->Divide(2,2);
  disp->SetCDSFrameXY(-70,70,-70,70);
  disp->SetCDSFrameYZ(-70,70,-70,70);
  const double ZLIM = 30, XYLIM = 10;
  disp->SetBLDCFrameXZ(-ZLIM,ZLIM,-XYLIM,XYLIM);
  disp->SetBLDCFrameYZ(-ZLIM,ZLIM,-XYLIM,XYLIM);
#endif

  /*                */
  /* event analysis */
  /*                */
  int nev = evtree->GetEntries();
  std::cout << " AllEvent : " << nev << std::endl;
  for( int iev=0; iev<nev; iev++ ){
    evtree->GetEvent(iev);

    if( iev%5000==0 )
      std::cout << " Event : " << iev << std::endl;

    CDSTrackingMan *trackMan = new CDSTrackingMan();
    BeamLineTrackMan *blTrackMan = new BeamLineTrackMan();

#if DEBUG
    cerr<<"# of T0 hits   = "<<blMan->nT0()<<endl;
    cerr<<"# of CDH hits  = "<<cdsMan->nCDH()<<endl;
    cerr<<"# of BLC1 hits = "<<blMan->nBLC1()<<endl;
    cerr<<"# of BLC2 hits = "<<blMan->nBLC2()<<endl;
    cerr<<"# of CDC hits  = "<<cdsMan->nCDC()<<endl;

    //=== dump ===//
    for(int i=0; i<blMan->nT0(); i++){
      HodoscopeLikeHit *hit = blMan->T0(i);
      cerr<<"T0 : "<<hit->hid()<<" "<<hit->seg()<<" "<<hit->ctmean()<<" "<<hit->emean()<<endl;
    }
    for(int i=0; i<cdsMan->nCDH(); i++){
      HodoscopeLikeHit *hit = cdsMan->CDH(i);
      cerr<<"CDH : "<<hit->hid()<<" "<<hit->seg()<<" "<<hit->ctmean()<<" "<<hit->emean()<<endl;
    }
    for(int i=1; i<=NumOfBLCLayers; i++){
      for(int j=0; j<blMan->nBLC1(i); j++){
	ChamberLikeHit *chm = blMan->BLC1(i, j);
	cerr<<"BLC1 : "<<chm->hid()<<" "<<chm->layer()<<" "<<chm->wire()<<" "<<chm->dl()<<endl;
      }
    }
    for(int i=1; i<=NumOfBLCLayers; i++){
      for(int j=0; j<blMan->nBLC2(i); j++){
	ChamberLikeHit *chm = blMan->BLC2(i, j);
	cerr<<"BLC2 : "<<chm->hid()<<" "<<chm->layer()<<" "<<chm->wire()<<" "<<chm->dl()<<endl;
      }
    }
    for(int i=1; i<=NumOfCDCLayers; i++){
      for(int j=0; j<cdsMan->nCDC(i); j++){
	CDCHit *cdc = cdsMan->CDC(i, j);
	cerr<<"CDC : "<<cdc->hid()<<" "<<cdc->layer()<<" "<<cdc->wire()<<" "<<cdc->dl()<<endl;
      }
    }
#endif

#if DISPLAY
    TVirtualPad *pad;
    pad = c_disp->GetPad(1);
    disp->DrawCDSFrameXY( pad );
    disp->DrawSegmentsXY( pad, confMan, CID_CDH );
    disp->DrawCDCLayersXY( pad, confMan );
    disp->DrawCDSHitXY( pad, confMan, cdsMan, CID_CDH );
    disp->DrawCDSHitXY( pad, confMan, cdsMan, CID_CDC );
    
    pad = c_disp->GetPad(2);
    disp->DrawCDSFrameYZ( pad );
    disp->DrawCDCLayersYZ( pad, confMan );
    disp->DrawCDSHitYZ( pad, confMan, cdsMan, CID_CDH );
    disp->DrawCDSHitYZ( pad, confMan, cdsMan, CID_CDC );

    pad = c_disp->GetPad(3);
    disp->DrawBLDCFrameXZ( pad );
    disp->DrawBLDCLayersXZ( pad, confMan, CID_BLC1 );
    disp->DrawBLDCLayersXZ( pad, confMan, CID_BLC2 );
    disp->DrawBLDCHit( pad, confMan, blMan, CID_BLC1, 0 );
    disp->DrawBLDCHit( pad, confMan, blMan, CID_BLC2, 0 );

    pad = c_disp->GetPad(4);
    disp->DrawBLDCFrameYZ( pad );
    disp->DrawBLDCLayersYZ( pad, confMan, CID_BLC1 );
    disp->DrawBLDCLayersYZ( pad, confMan, CID_BLC2 );
    disp->DrawBLDCHit( pad, confMan, blMan, CID_BLC1, 1 );
    disp->DrawBLDCHit( pad, confMan, blMan, CID_BLC2, 1 );

    c_disp->Update();
#endif    

    // --- tracking ---//
    trackMan->Execute(cdsMan,confMan);
    blTrackMan->DoTracking(blMan);

    cerr<<"# of CDC tracks  = "<<trackMan->nTrack()
	<<" (good track :"<<trackMan->nGoodTrack()<<")"<<endl;
    cerr<<"# of BLC tracks = "<<blTrackMan->ntrackBLC2()<<endl;

    for(int i=0; i<trackMan->nTrack(); i++){
      cerr<<trackMan->Track(i)->Chi()<<endl;
    }
    for(int i=0; i<blTrackMan->ntrackBLC2(); i++){
      cerr<<blTrackMan->trackBLC2(i)->chi2xz()<<" "<<blTrackMan->trackBLC2(i)->chi2yz()<<endl;
    }


#if DISPLAY
    //--- CDC circle tracks ---//
    double arho[10],ax_c[10],ay_c[10],aPt[10];
    double aparam[10][5];
    for(int i=0; i<trackMan->nTrack(); i++){
      CDSTrack *track=trackMan->Track(i);
      track->GetParameters(aparam[i]);
      arho[i]=fabs(1./aparam[i][2]);
      ax_c[i]=(aparam[i][0]+1./aparam[i][2])*cos(aparam[i][1]);
      ay_c[i]=(aparam[i][0]+1./aparam[i][2])*sin(aparam[i][1]);
      aPt[i]=0.3*0.5*arho[i]/100.;
    }
    pad = c_disp->GetPad(1);
    pad->cd();
    TArc arc[10];
    for(int i=0; i<trackMan->nTrack(); i++){
      arc[i].SetFillStyle(0);
      arc[i].SetLineColor(1+i);
      arc[i].DrawArc(ax_c[i],ay_c[i],arho[i],0,360);
    }

    //--- BLC line tracks ---//
    double la[5],lb[5],lc[5],ld[5],le[5],lf[5];
    for(int i=0; i<blTrackMan->ntrackBLC2(); i++){
      blTrackMan->trackBLC2(i)->abc(la[i],lb[i],lc[i]);
      blTrackMan->trackBLC2(i)->def(ld[i],le[i],lf[i]);
    }
    TLine lineZX[5], lineZY[5];
    pad = c_disp->GetPad(3);
    pad->cd();
    for(int i=0; i<blTrackMan->ntrackBLC2(); i++){
      lineZX[i].SetLineStyle(1);
      lineZX[i].SetLineColor(2+i);
      double xmin = -(lb[i]*(-ZLIM)+lc[i])/la[i];
      double xmax = -(lb[i]*(ZLIM)+lc[i])/la[i];
      lineZX[i].DrawLine(-ZLIM, xmin, ZLIM, xmax);
    }
    pad = c_disp->GetPad(4);
    pad->cd();
    for(int i=0; i<blTrackMan->ntrackBLC2(); i++){
      lineZY[i].SetLineStyle(1);
      lineZY[i].SetLineColor(2+i);
      double ymin = -(le[i]*(-ZLIM)+lf[i])/ld[i];
      double ymax = -(le[i]*(ZLIM)+lf[i])/ld[i];
      lineZY[i].DrawLine(-ZLIM, ymin, ZLIM, ymax);
    }

    c_disp->Update();
    bool status = disp->Wait();
    if( !status ) return;
#endif 

    delete trackMan;
    delete blTrackMan;

  } //for( int iev=0; iev<nev; iev++ ){

}
