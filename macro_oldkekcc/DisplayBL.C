#include <fstream.h>
#include <iostream.h>
#include <vector.h>

#define DISPLAY 1

class ConfMan;
class Display;
class CDSHitMan;

void DisplayBL()
{
  /*** load library ***/
  gSystem->Load("lib/libAll.so");

  /*** assign input file & call tree ***/
  TFile *f = new TFile("root/1164.root");
  TTree *evtree = (TTree*)f->Get("EventTree");

  /*** conf file for new parameters ***/
  ConfMan *conf = new ConfMan("conf/730.conf");
  conf->Initialize();

  /* --------------- */
  /* ---- start ---- */
  /* --------------- */

#if DISPLAY
  Display *disp = new Display();
  TCanvas *c_disp = new TCanvas( "c_disp", "c_disp", 1000, 0, 400, 800 );
  disp->SetBLFrameXZ(-90,90,-650,1550);
#endif

  /*** declaration of classes ***/
  TKOHitCollection *tko = 0;
  BeamLineHitMan *blMan = 0;
  EventHeader *head = 0;
  evtree->SetBranchAddress( "TKOHitCol", &tko );
  evtree->SetBranchAddress( "BeamLineHitMan", &blMan );
  evtree->SetBranchAddress( "EventHeader", &head );
  
  TH1F *h1;
  char dispflag;
  char dispin[100]="";

  /*                */
  /* event analysis */
  /*                */
  int nev = evtree->GetEntries();
  for( int iev=0; iev<nev; iev++ ){
    evtree->GetEvent(iev);

    if( iev%5000==0 )
      std::cout << " Event : " << iev << std::endl;

    /*** if some parameters are newed, you can re-calc all quantities as below, ***/
    blMan->Calc(conf);

    if( !head->proton() ) continue; // select proton trigger

#if DISPLAY
    disp->DrawBLFrameXZ( c_disp );
    disp->DrawSegmentsXZ( c_disp, conf, CID_TOFstop );
    disp->DrawSegmentsXZ( c_disp, conf, CID_T0 );
    disp->DrawSegmentsXZ( c_disp, conf, CID_B1 );
    disp->DrawSegmentsXZ( c_disp, conf, CID_B2 );
    disp->DrawSegmentsXZ( c_disp, conf, CID_PA );
    disp->DrawSegmentsXZ( c_disp, conf, CID_BHD );
    disp->DrawBLHitXZ( c_disp, conf, blMan, CID_TOFstop );
    disp->DrawBLHitXZ( c_disp, conf, blMan, CID_T0 );
    disp->DrawBLHitXZ( c_disp, conf, blMan, CID_B1 );
    disp->DrawBLHitXZ( c_disp, conf, blMan, CID_B2 );
    disp->DrawBLHitXZ( c_disp, conf, blMan, CID_PA );
    disp->DrawBLHitXZ( c_disp, conf, blMan, CID_BHD );
    
    c_disp->Update();

    bool status = disp->Wait();
    if( !status ) return;
    //c_disp->Clear();
#endif    
  }
}

