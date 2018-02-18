#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TVirtualPad.h>
#include <TLine.h>
#include <TRint.h>

#include "File.h"
#include "ConfMan.h"
#include "TKO.h"
#include "CDSHitMan.h"
#include "BeamLineHitMan.h"
#include "BeamLineTrackMan.h"
#include "EventHeader.h"
#include "ScalerMan.h"

#include "EventAlloc.h"
#include "EventTemp.h"
#include "Display.h"

class EventDisplay: public EventTemp
{
public:
  EventDisplay();
  ~EventDisplay();
private:
  Display *disp;
  TCanvas *canvas1;

  CDSHitMan *cdsMan;
  BeamLineHitMan *blMan;
  EventHeader *header;
  ScalerMan *scaMan;
public:
  void Initialize( ConfMan *conf );
  void USca( int nsca, unsigned int *sca );
  bool UAna( TKOHitCollection *tko );
  void Finalize();
};

EventDisplay::EventDisplay()
  : EventTemp()
{
}

EventDisplay::~EventDisplay()
{
}

const int MaxTreeSize = 1900000000000;
void EventDisplay::Initialize( ConfMan *conf )
{
#if 1
  std::cout << " Enter EventDisplay::Initialize " << std::endl;
#endif

  int argc2;
  char **argv2;
  TRint *theApp = new TRint( "theApp", &argc2, argv2 );

  confMan = conf;
  header = new EventHeader();
  cdsMan = new CDSHitMan();
  blMan = new BeamLineHitMan();
  //scaMan = new ScalerMan();

  disp = new Display();
  canvas1 = new TCanvas( "canvas1", "canvas1", 800,400 );
  std::cout << " disp:" << disp << " canvas1:" << canvas1 << std::endl;
  canvas1->Divide(2,1);
  //disp->SetCDSFrameXY(-70,70,-70,70);
  //disp->SetCDSFrameYZ(-70,70,-70,70);
  //  disp->SetBLFrameXZ(-90,90,-650,1550);
  disp->SetBLDCFrameXZ(-25,25,-15,15);
  disp->SetBLDCFrameYZ(-25,25,-15,15);
}

void EventDisplay::USca( int nsca, unsigned int *sca )
{
#if 0
  std::cout << " Enter EventDisplay::USca " << std::endl;
#endif
  Block_Event_Number++;
  header->SetBlockEventNumber( Block_Event_Number );
  //scaMan->SetBlockEventNumber( Block_Event_Number );
  //for( int i=0; i<nsca; i++ ){
  //  scaMan->AddHit( sca[i], confMan->GetScalerMapManager()->GetName(i) );
  //}
#if 0
  //std::cout << nsca<< std::endl;
  //for( int i=0; i<nsca; i++ ){
  //  std::cout << "  " << sca[i];
  //}
  //std::cout << std::endl;
#endif
  //scaMan->Clear();
}

bool EventDisplay::UAna( TKOHitCollection *tko )
{
#if 1
  std::cout << " Enter EventDisplay::UAna " << std::endl;
#endif

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
  cdsMan->Convert( tko, confMan );
  /*  
  disp->DrawCDSFrameXY( pad );
  disp->DrawSegmentsXY( pad, confMan, CID_CDH );
  disp->DrawCDCLayersXY( pad, confMan );
  disp->DrawCDSHitXY( pad, confMan, cdsMan, CID_CDH );
  disp->DrawCDSHitXY( pad, confMan, cdsMan, CID_CDC );

  pad = canvas1->GetPad(2);
  disp->DrawCDSFrameYZ( pad );
  disp->DrawCDCLayersYZ( pad, confMan );
  disp->DrawCDSHitYZ( pad, confMan, cdsMan, CID_CDH );
  disp->DrawCDSHitYZ( pad, confMan, cdsMan, CID_CDC );
  */
  /*  
  disp->DrawBLFrameXZ( pad );
  disp->DrawSegmentsXZ( pad, confMan, CID_TOFstop );
  disp->DrawSegmentsXZ( pad, confMan, CID_T0 );
  //  disp->DrawSegmentsXZ( pad, confMan, CID_B1 );
  // disp->DrawSegmentsXZ( pad, confMan, CID_B2 );
  disp->DrawSegmentsXZ( pad, confMan, CID_PA );
  disp->DrawSegmentsXZ( pad, confMan, CID_BHD );
  disp->DrawBLHitXZ( pad, confMan, blMan, CID_TOFstop );
  disp->DrawBLHitXZ( pad, confMan, blMan, CID_T0 );
  // disp->DrawBLHitXZ( pad, confMan, blMan, CID_B1 );
  //disp->DrawBLHitXZ( pad, confMan, blMan, CID_B2 );
  disp->DrawBLHitXZ( pad, confMan, blMan, CID_PA );
  disp->DrawBLHitXZ( pad, confMan, blMan, CID_BHD );
  */ 

  BeamLineTrackMan *blTrackMan=new BeamLineTrackMan();
  if(!blTrackMan->DoTracking(blMan) )   std::cout<<"miss traking !"<<std::endl;
;

  TVirtualPad *pad;
  pad = canvas1->GetPad(1);
  disp->DrawBLDCFrameXZ( pad );
  disp->DrawBLDCLayersXZ( pad, confMan, CID_BLC1 );
  disp->DrawBLDCLayersXZ( pad, confMan, CID_BLC2 );
  disp->DrawBLDCHit( pad, confMan,blMan, CID_BLC1,0 );
  disp->DrawBLDCHit( pad, confMan,blMan, CID_BLC2,0 );
  std::cout<<"nTrack BLC2 "<<blTrackMan->ntrackBLC2()<<std::endl;
  for(int itr=0;itr<blTrackMan->ntrackBLC2();itr++)
    {
      LocalTrack *blc2=blTrackMan->trackBLC2(itr);
      double z1,z2,x1,x2,y1,y2;
      z1=-25;z2=25;
      blc2->XYPosatZ(z1,x1,y1);
      blc2->XYPosatZ(z2,x2,y2);
      std::cout<<"chi zx : zy : all "<<blc2->chi2xz()<<" : "<<blc2->chi2yz()<<" : "<<blc2->chi2all()<<std::endl;
      TLine line;
      line.SetLineColor(itr+2);
      line.DrawLine(z1,x1,z2,x2);
      //      std::cout<<"line "<<itr<<" z,x="<<z1<<" "<<x1<<std::endl;
    }
  
  pad = canvas1->GetPad(2);
  disp->DrawBLDCFrameYZ( pad );
  disp->DrawBLDCLayersYZ( pad, confMan, CID_BLC1 );
  disp->DrawBLDCLayersYZ( pad, confMan, CID_BLC2 );
  disp->DrawBLDCHit( pad, confMan,blMan, CID_BLC1,1 );
  disp->DrawBLDCHit( pad, confMan,blMan, CID_BLC2,1 );
  for(int itr=0;itr<blTrackMan->ntrackBLC2();itr++)
    {
      LocalTrack *blc2=blTrackMan->trackBLC2(itr);
      double z1,z2,x1,x2,y1,y2;
      z1=-25;z2=25;
      blc2->XYPosatZ(z1,x1,y1);
      blc2->XYPosatZ(z2,x2,y2);
      TLine line;
      line.SetLineColor(itr+2);
      line.DrawLine(z1,y1,z2,y2);
    }
  

  canvas1->Update();
  canvas1->SaveAs("tmp.ps");
  bool status = disp->Wait();
  if( !status ) return status;

  header->Clear();
  blMan->Clear();
  cdsMan->Clear();
  return true;
}

void EventDisplay::Finalize()
{
  std::cout << " Enter EventDisplay::Finalize " << std::endl;

  delete blMan;
  delete cdsMan;
  delete header;

  gSystem->Exit(1);
  gROOT->GetApplication()->Run();

}

EventTemp *EventAlloc::EventAllocator()
{
  EventDisplay *event = new EventDisplay();
  return (EventTemp*)event;
}
