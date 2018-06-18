#include "EventTemp.h"
#include "EventAlloc.h"

#include "EventHeader.h"
#include "ScalerMan.h"
#include "BeamSpectrometer.h"
#include "TrackTools.h"
#include "MyHistTools.h"

class EventAnalysisBT : public EventTemp
{
public:
  EventAnalysisBT();
  ~EventAnalysisBT(){};
private:
  int CDC_Event_Number;

  TFile *rtFile;
  TTree *evTree;
  TTree *scaTree;
  EventHeader *header;
  ScalerMan *scaMan;

  BeamLineHitMan *blMan;
  CDSHitMan *cdsMan;
  BeamLineTrackMan *bltrackMan;
  BeamSpectrometer *beamSpec;

public:
  void Initialize(ConfMan *conf);
  void USca(int usca, unsigned int *sca);
  bool UAna(TKOHitCollection *tko);
  bool UAnaBeam();
  void Finalize();

  bool Clear(bool flag=true);
};

EventAnalysisBT::EventAnalysisBT() : EventTemp()
{
}

void EventAnalysisBT::Initialize(ConfMan *conf)
{
  confMan = conf;
  
  rtFile = new TFile(confMan->GetOutFileName().c_str(), "recreate");
  evTree = new TTree("EventTree", "EventTree");
  header = new EventHeader();
  blMan = new BeamLineHitMan();
  bltrackMan = new BeamLineTrackMan();
  beamSpec = new BeamSpectrometer(confMan);
  cdsMan =  new CDSHitMan();

  scaTree = new TTree("ScalerTree", "ScalerTree");
  scaMan = new ScalerMan();
  scaTree-> Branch("ScalerMan", &scaMan);
  rtFile-> cd();
  MyHistTools::initT0();
}

bool EventAnalysisBT::UAna(TKOHitCollection *tko)
{
  Event_Number++;
  if( Event_Number%5000==0 )
    std::cout << " Event# : " << Event_Number << "  BlockEvent# : " << Block_Event_Number << std::endl;


  int status = confMan->CheckEvNum( Event_Number, Block_Event_Number );
  if( status==1 ){
    return Clear();
  }
  if( status==2 ){
    return Clear(false);
  }

  header-> SetRunNumber(confMan->GetRunNumber());
  header-> SetEventNumber(Event_Number);
  header->Convert( tko, confMan );
  if( header->IsTrig(Trig_Cosmic) || header->IsTrig(Trig_Reject) ) return Clear();
  blMan-> Convert(tko, confMan);
  cdsMan-> Convert(tko, confMan);
  UAnaBeam();

  return Clear();
}


bool EventAnalysisBT::Clear(bool flag)
{
  header-> Clear();
  blMan-> Clear();
  bltrackMan-> Clear();
  beamSpec-> Clear();
  cdsMan-> Clear();
  return flag;
}

void EventAnalysisBT::USca(int nsca, unsigned int *sca)
{
  Block_Event_Number++;
  header->SetBlockEventNumber( Block_Event_Number );
  scaMan->SetBlockEventNumber( Block_Event_Number );
  for( int i=0; i<nsca; i++ ){
    scaMan->AddHit( sca[i], confMan->GetScalerMapManager()->GetName(i) );
  }
  scaTree->Fill();
  scaMan->Clear();
}

void EventAnalysisBT::Finalize()
{
  gFile-> Write();
  gFile-> Close();

  delete blMan;
  delete cdsMan;
  delete bltrackMan;
  delete beamSpec;
}
