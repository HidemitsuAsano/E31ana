#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TVirtualPad.h>
#include <TRint.h>

#include "File.h"
#include "ConfMan.h"
#include "TKO.h"
#include "CDSHitMan.h"
#include "BeamLineHitMan.h"
#include "BeamLineTrackMan.h"
#include "EventHeader.h"
#include "ScalerMan.h"

#include "CDSTrackingMan.h"
#include "EventAlloc.h"
#include "EventTemp.h"
#include "Display.h"

#define TRACK 1
#define LINEAR 0

class BLDCDisplay: public EventTemp
{
public:
  BLDCDisplay();
  ~BLDCDisplay();
private:
  Display *disp;
  TCanvas *canvas1;
  TCanvas *canvas2;

  CDSHitMan *cdsMan;
  BeamLineHitMan *blMan;
  //  BeamLineTrackMan *blTrackMan;
  EventHeader *header;
  ScalerMan *scaMan;
public:
  void Initialize( ConfMan *conf );
  void USca( int nsca, unsigned int *sca );
  bool UAna( TKOHitCollection *tko );
  void Finalize();
};

BLDCDisplay::BLDCDisplay()
  : EventTemp()
{
}

BLDCDisplay::~BLDCDisplay()
{
}

const int MaxTreeSize = 1000000000;
void BLDCDisplay::Initialize( ConfMan *conf )
{
#if 1
  std::cout << " Enter BLDCDisplay::Initialize " << std::endl;
#endif
  
  int argc2=1;
  char *temp="bldcdisp";
  //  strcpy(temp,conf->GetProgramName().data());
  char **argv2=&temp;
  new TRint( "theApp", &argc2, argv2 );
  
  confMan = conf;
  std::cout<<"make header"<<std::endl;  
  header = new EventHeader();
  //  cdsMan = new CDSHitMan();
  std::cout<<"make blman"<<std::endl;  
  blMan = new BeamLineHitMan();

  //scaMan = new ScalerMan();
  //  std::cout<<"make bltrackman"<<std::endl;  
  std::cout<<"make disp"<<std::endl;  
  disp = new Display();
  canvas1 = new TCanvas( "canvas1", "canvas1", 800, 800 );
  canvas2 = new TCanvas( "canvas2", "canvas2", 600, 600 );
  std::cout << " disp:" << disp << " canvas1:" << canvas1 << std::endl;
  canvas1->Divide(3,3);
  canvas2->Divide(1,2);
}

void BLDCDisplay::USca( int nsca, unsigned int *sca )
{
#if 0
  std::cout << " Enter BLDCDisplay::USca " << std::endl;
#endif
  Block_Event_Number++;
  header->SetBlockEventNumber( Block_Event_Number );
}

bool BLDCDisplay::UAna( TKOHitCollection *tko )
{
#if 1
  std::cout << " Enter BLDCDisplay::UAna " << std::endl;
#endif
  BeamLineTrackMan* blTrackMan = new BeamLineTrackMan();
  Event_Number++;
  { int status = confMan->CheckEvNum( Event_Number, Block_Event_Number );
  if( status==1 ) return true;
  if( status==2 ) return false; }
  
  if( Event_Number%1==0 )
    std::cout << " Event# : " << Event_Number << std::endl;
  
  header->SetRunNumber(0);
  header->SetEventNumber(Event_Number);
  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );

#if TRACK
  std::cout<<"start tracking "<<std::endl;
  blTrackMan->DoTracking(blMan,confMan,true,true);
#endif
  //  if(1){
  //  if(blTrackMan->ntrackBLC1()==1&&blTrackMan->ntrackBLC2()<1&&blTrackMan->ntrackBPC()==1){
  if(blTrackMan->ntrackBLC2()==1 && blTrackMan->trackBLC2(0)->chi2all()>60){
    //    if(blTrackMan->status(CID_FDC1)>1)
    TVirtualPad *pad;    
    //BLC1
    pad = canvas1->GetPad(1);
    disp->SetBLDCFrameXZ(-25,25,-15,15);
    disp->DrawBLDCFrameXZ( pad );
    disp->DrawBLDCXZ( pad, confMan,blMan, blTrackMan, CID_BLC1a );
    disp->DrawBLDCXZ( pad, confMan,blMan, blTrackMan, CID_BLC1b );
    disp->DrawBLDCTrackXZ( pad, confMan, blTrackMan, CID_BLC1, 4 );

    pad = canvas1->GetPad(4);
    disp->SetBLDCFrameYZ(-25,25,-15,15);
    disp->DrawBLDCFrameYZ( pad );
    disp->DrawBLDCYZ( pad, confMan,blMan, blTrackMan, CID_BLC1a );
    disp->DrawBLDCYZ( pad, confMan,blMan, blTrackMan, CID_BLC1b );
    disp->DrawBLDCTrackYZ( pad, confMan, blTrackMan, CID_BLC1, 4 );


    pad = canvas1->GetPad(7);
    disp->SetBLDCFrameYZ(-25,25,-100,300);
    disp->DrawBLDCFrameYZ( pad );
    disp->DrawClusterTime( pad, confMan,blTrackMan, CID_BLC1 );

    //BLC2
    pad = canvas1->GetPad(2);
    disp->SetBLDCFrameXZ(-25,25,-10,10);
    disp->DrawBLDCFrameXZ( pad );
    disp->DrawSegmentsXZ( pad, confMan,CID_T0 );
    disp->DrawBLHitYZ( pad, confMan, blMan, CID_T0 );
    disp->DrawBLDCXZ( pad, confMan,blMan, blTrackMan, CID_BLC2a );
    disp->DrawBLDCXZ( pad, confMan,blMan, blTrackMan, CID_BLC2b );
    disp->DrawBLDCTrackXZ( pad, confMan, blTrackMan, CID_BLC2, 4 );

    pad = canvas1->GetPad(5);
    disp->SetBLDCFrameYZ(-25,25,-10,10);
    disp->DrawBLDCFrameYZ( pad );
    disp->DrawSegmentsXZ( pad, confMan,CID_T0 );
    disp->DrawBLHitXZ( pad, confMan, blMan, CID_T0 );
    disp->DrawBLDCYZ( pad, confMan,blMan, blTrackMan, CID_BLC2a );
    disp->DrawBLDCYZ( pad, confMan,blMan, blTrackMan, CID_BLC2b );
    disp->DrawBLDCTrackYZ( pad, confMan, blTrackMan, CID_BLC2, 4 );

    pad = canvas1->GetPad(8);
    disp->SetBLDCFrameYZ(-25,25,-100,300);
    disp->DrawBLDCFrameYZ( pad );
    disp->DrawClusterTime( pad, confMan, blTrackMan, CID_BLC2 );

    //BPC
    pad = canvas1->GetPad(3);
    disp->SetBLDCFrameXZ(-10,10,-8,8);
    disp->DrawBLDCFrameXZ( pad );
    disp->DrawBLDCXZ( pad, confMan,blMan, blTrackMan, CID_BPC );           
    pad = canvas1->GetPad(6);
    disp->SetBLDCFrameYZ(-10,10,-8,8);
    disp->DrawBLDCFrameYZ( pad );
    disp->DrawBLDCYZ( pad, confMan,blMan, blTrackMan, CID_BPC );
    pad = canvas1->GetPad(9);
    disp->SetBLDCFrameYZ(-10,10,-100,300);
    disp->DrawBLDCFrameYZ( pad );
    disp->DrawClusterTime( pad, confMan, blTrackMan, CID_BPC );

#if FDC
    pad = canvas2->GetPad(1);
    disp->SetBLDCFrameYZ(-30,30,-20,20);
    disp->DrawBLDCFrameYZ( pad );
    disp->DrawFDC( pad, confMan,blMan );

    pad = canvas2->GetPad(2);
    disp->SetBLDCFrameYZ(-30,30,-10,10);
    disp->DrawBLDCFrameYZ( pad );
    disp->DrawBLDCXZ( pad, confMan,blMan, blTrackMan, CID_FDC1 );           
#endif
    std::cout<<"tof= "<<ctof[0][1]<<std::endl;
    std::cout<<"BPC tracking status: "<<blTrackMan->status(CID_BPC)<<std::endl;
    canvas1->Update();
    canvas2->Update();
    canvas1->SaveAs("tmpdis.ps");
    bool status = disp->Wait();
    if( !status ) return status;
  }
  header->Clear();
  blMan->Clear();
  //  blTrackMan->Clear();
  std::cout<<"clear bltrackman"<<std::endl;
  delete blTrackMan;
  return true;
}

void BLDCDisplay::Finalize()
{
  std::cout << " Enter BLDCDisplay::Finalize " << std::endl;  
  delete blMan;
  delete header;  
  gSystem->Exit(1);
  gROOT->GetApplication()->Run();
}

EventTemp *EventAlloc::EventAllocator()
{
  BLDCDisplay *event = new BLDCDisplay();
  return (EventTemp*)event;
}
