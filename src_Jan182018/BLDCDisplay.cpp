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

class BLDCDisplay: public EventTemp
{
public:
  BLDCDisplay();
  ~BLDCDisplay();
private:
  Display *disp;
  TCanvas *canvas1;

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
  
  int argc2;
  char **argv2;
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
  canvas1 = new TCanvas( "canvas1", "canvas1", 1200, 700 );
  std::cout << " disp:" << disp << " canvas1:" << canvas1 << std::endl;
  canvas1->Divide(3,2);
  //  disp->SetCDSFrameXY(-70,70,-70,70);
  //  disp->SetCDSFrameYZ(-70,70,-70,70);
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
  std::cout<<"convert finished"<<std::endl;

  /// TOF analysis start
  int cid[3]={CID_BHD,CID_T0,CID_TOFstop};
  TString name[3]={"BHD","T0","TOF"};
  double ctmHodo[3];
  int nHodo[3]={0,0,0};
  int Hodoseg[3]={0,0,0};
  for(int ic=0;ic<3;ic++)
    for( int i = 0; i< blMan->nHodo(cid[ic]); i++ ){
      HodoscopeLikeHit *hit = blMan->Hodoi(cid[ic],i);
      if(hit->CheckRange()){
	nHodo[ic]++;
	Hodoseg[ic]=hit->seg();
	if(nHodo[ic]==1)
	  ctmHodo[ic]=hit->ctmean();
      }
    }    
  
  double ctof[3][3];    
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      ctof[i][j]= TMath::Abs(ctmHodo[i]-ctmHodo[j]);
  std::cout<<"tof derived"<<std::endl;

  //TOF analysis end
  //  if(header->proton()&&ctof[0][1]>30&&nHodo[1]==1&&nHodo[2]==1){  
  // Mometnum analysis using TOF
  if(0){
    double t0x,t0y,t0z,tofsx,tofsy,tofsz;
    blMan->Hodo(CID_T0,Hodoseg[1])->pos(t0x,t0y,t0z);
    blMan->Hodo(CID_TOFstop,Hodoseg[2])->pos(tofsx,tofsy,tofsz);
    double dist=TMath::Sqrt(1710.*1710.+tofsx*tofsx);
    //  std::cout<<dist<<std::endl;
    double c=29.979245; // cm/ns
    double beta=dist/c/(dist/c+ctof[2][1]-57.0395);
    double mom=beta/TMath::Sqrt(1-beta*beta)*938.3+8.;
    std::cout<<"Momentum: "<<mom<<" MeV/c"<<std::endl;
  }
  // Mometnum analysis using TOF end
#if TRACK
  std::cout<<"start tracking "<<std::endl;
  blTrackMan->LocalTracking(blMan,confMan,CID_BPC);
  blTrackMan->LinearTracking(blMan,confMan,CID_BLC1);
  blTrackMan->LinearTracking(blMan,confMan,CID_BLC2);
  std::cout<<"tracking finished"<<std::endl;
  //  int nGoodTrack=blTrackMan->nGoodTrack();
#endif
  //  if(blTrackMan->ntrackBLC()==1){
  //  if(blTrackMan->ntrackPDC1()==1&&blTrackMan->trackPDC1(0)->chi2all()>100){
  if(1){
  //  if(blTrackMan->status(CID_BLC1)==4){
    TVirtualPad *pad;
    
    //BLC1
    pad = canvas1->GetPad(1);
    disp->SetBLDCFrameXZ(-25,25,-15,15);
    disp->DrawBLDCFrameXZ( pad );
    disp->DrawBLDCXZ( pad, confMan,blMan, blTrackMan, CID_BLC1a );
    disp->DrawBLDCXZ( pad, confMan,blMan, blTrackMan, CID_BLC1b );
    disp->DrawBLDCTrackXZ( pad, confMan, blTrackMan, CID_BLC, 1 );
    pad = canvas1->GetPad(4);
    disp->SetBLDCFrameYZ(-25,25,-15,15);
    disp->DrawBLDCFrameYZ( pad );
    disp->DrawBLDCYZ( pad, confMan,blMan, blTrackMan, CID_BLC1a );
    disp->DrawBLDCYZ( pad, confMan,blMan, blTrackMan, CID_BLC1b );
    disp->DrawBLDCTrackYZ( pad, confMan, blTrackMan, CID_BLC1, 2 );
    
    //BLC2
    pad = canvas1->GetPad(2);
    disp->SetBLDCFrameXZ(-25,25,-10,10);
    disp->DrawBLDCFrameXZ( pad );
    disp->DrawSegmentsXZ( pad, confMan,CID_T0 );
    disp->DrawBLHitYZ( pad, confMan, blMan, CID_T0 );
    disp->DrawBLDCXZ( pad, confMan,blMan, blTrackMan, CID_BLC2a );
    disp->DrawBLDCXZ( pad, confMan,blMan, blTrackMan, CID_BLC2b );
    disp->DrawBLDCTrackXZ( pad, confMan, blTrackMan, CID_BLC2, 2 );
    pad = canvas1->GetPad(5);
    disp->SetBLDCFrameYZ(-25,25,-10,10);
    disp->DrawBLDCFrameYZ( pad );
    disp->DrawSegmentsXZ( pad, confMan,CID_T0 );
    disp->DrawBLHitXZ( pad, confMan, blMan, CID_T0 );
    disp->DrawBLDCYZ( pad, confMan,blMan, blTrackMan, CID_BLC2a );
    disp->DrawBLDCYZ( pad, confMan,blMan, blTrackMan, CID_BLC2b );
    disp->DrawBLDCTrackYZ( pad, confMan, blTrackMan, CID_BLC2, 2 );
    
    //BPC
    pad = canvas1->GetPad(3);
    disp->SetBLDCFrameXZ(-10,10,-8,8);
    disp->DrawBLDCFrameXZ( pad );
    disp->DrawBLDCXZ( pad, confMan,blMan, blTrackMan, CID_BPC );           
    pad = canvas1->GetPad(6);
    disp->SetBLDCFrameYZ(-10,10,-8,8);
    disp->DrawBLDCFrameYZ( pad );
    disp->DrawBLDCYZ( pad, confMan,blMan, blTrackMan, CID_BPC );
    std::cout<<"tof= "<<ctof[0][1]<<std::endl;

    canvas1->Update();
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
  delete cdsMan;
  delete header;  
  gSystem->Exit(1);
  gROOT->GetApplication()->Run();
}

EventTemp *EventAlloc::EventAllocator()
{
  BLDCDisplay *event = new BLDCDisplay();
  return (EventTemp*)event;
}
