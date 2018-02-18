#include "EventTemp.h"
#include "EventAlloc.h"

#include "EventHeader.h"
#include "ScalerMan.h"
#include "BeamSpectrometer.h"
#include "TrackTools.h"

class EventAnalysis : public EventTemp
{
public:
  EventAnalysis();
  ~EventAnalysis(){};
private:
  int t0, t1;
  int CDC_Event_Number;
  TFile *cdsFile;
  TTree *cdsTree;
  EventHeader *cdsHeader;
  CDSTrackingMan *cdstrackMan;

  TFile *rtFile;
  TTree *scaTree;
  EventHeader *header;
  ScalerMan *scaMan;

  BeamLineHitMan *blMan;
  BeamLineTrackMan *bltrackMan;
  BeamSpectrometer *beamSpec;
  CDSHitMan *cdsMan;

public:
  void Initialize(ConfMan *conf);
  void USca(int usca, unsigned int *sca);
  bool UAna(TKOHitCollection *tko);
  bool UAna();
  void Finalize();

private:
  void InitializeHistogram();
  bool cdsFileMatching();
  bool Clear(bool flag=true);
};

EventAnalysis::EventAnalysis() : EventTemp()
{
}

bool EventAnalysis::UAna(TKOHitCollection *tko)
{
  Event_Number++;
  if( Event_Number%1000==0 ){
    t1 = time(0);
    std::cout << " Event# : " << Event_Number << "  BlockEvent# : " << Block_Event_Number << " Time (s): "<<(t1-t0)<<std::endl;
  }

  int status = confMan->CheckEvNum( Event_Number, Block_Event_Number );
  if( status==1 ){
    return Clear();
  }
  if( status==2 ){
    return Clear(false);
  }
  header-> SetRunNumber(confMan->GetRunNumber());
  header-> SetEventNumber(Event_Number);
  header-> Convert( tko, confMan );
  if( header->IsTrig(Trig_Cosmic) || header->IsTrig(Trig_Reject) ) return Clear();

  blMan-> Convert(tko, confMan);
  cdsMan-> Convert(tko, confMan);
  bltrackMan-> DoTracking(blMan, confMan, true, true);

  return Clear(UAna());
}

void EventAnalysis::Initialize(ConfMan *conf)
{
  confMan = conf;

  cdsFile = new TFile(confMan->GetCDSTrackFileName().c_str());
  cdstrackMan = new CDSTrackingMan();
  cdsHeader = new EventHeader();
  CDC_Event_Number=0;
  cdsTree = (TTree*)cdsFile-> Get("EventTree");
  cdsTree-> SetBranchAddress("EventHeader", &cdsHeader);  
  cdsTree-> SetBranchAddress("CDSTrackingMan", &cdstrackMan);

  rtFile = new TFile(confMan->GetOutFileName().c_str(), "recreate");
  header = new EventHeader();
  blMan = new BeamLineHitMan();
  bltrackMan = new BeamLineTrackMan();
  beamSpec = new BeamSpectrometer(confMan);
  cdsMan =  new CDSHitMan();

  scaTree = new TTree("ScalerTree", "ScalerTree");
  scaMan = new ScalerMan();
  scaTree-> Branch("ScalerMan", &scaMan);

  t0 = time(0);

  InitializeHistogram();
}

bool EventAnalysis::cdsFileMatching()
{
  cdsTree -> GetEntry(CDC_Event_Number);
  if(CDC_Event_Number>=cdsTree->GetEntries()){
    cdsHeader->Clear();
    cdstrackMan->Clear();
    return false;
  }
  while( cdsHeader->ev()<Event_Number ){
    CDC_Event_Number++;
    if(CDC_Event_Number>cdsTree->GetEntries()){
      cdsHeader->Clear();
      cdstrackMan->Clear();
      return false;
    }
    cdsTree -> GetEntry(CDC_Event_Number);
  }
  if(cdsHeader->ev()>Event_Number){
    cdsHeader->Clear();
    cdstrackMan->Clear();
    return false;
  }
  cdstrackMan->Calc(cdsMan, confMan);

  return true;
}

bool EventAnalysis::Clear(bool flag)
{
  header-> Clear();
  blMan-> Clear();
  bltrackMan-> Clear();
  beamSpec-> Clear();
  cdsMan-> Clear();

  cdsHeader-> Clear();
  cdstrackMan-> Clear();

  return flag;
}

void EventAnalysis::USca(int nsca, unsigned int *sca)
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

void EventAnalysis::Finalize()
{
  rtFile-> cd();
  rtFile-> Write();
  rtFile-> Close();

  delete blMan;
  delete bltrackMan;
  delete beamSpec;
  delete cdsMan;
  delete cdsHeader;
  delete cdstrackMan;
  delete cdsFile;
}
