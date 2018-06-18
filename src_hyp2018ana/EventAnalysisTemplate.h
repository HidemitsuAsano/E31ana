#include "EventTemp.h"
#include "EventAlloc.h"

#include "EventHeader.h"
#include "ScalerMan.h"
#include "BeamLineTrackMan.h"
#include "BeamSpectrometer.h"
#include "TrackTools.h"
#include "GeomTools.h"

#include "UserTools.h"

class EventAnalysisTemplate : public EventTemp
{
public:
  EventAnalysisTemplate();
  ~EventAnalysisTemplate(){};
private:
  int t0, t1;

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
  bool anaFileMatching();
  bool Clear(bool flag=true);
};

EventAnalysisTemplate::EventAnalysisTemplate() : EventTemp()
{
}

bool EventAnalysisTemplate::UAna(TKOHitCollection *tko)
{
  Event_Number++;
  if( Event_Number%1000==0 ){
    t1 = time(0);
    std::cout << " Event# : " << Event_Number << "  BlockEvent# : " << Block_Event_Number << " Time (s): "<<(t1-t0)<<std::endl;
  }

  int status = confMan->CheckEvNum( Event_Number, Block_Event_Number );
  if( status==1 ){ return Clear(); }
  if( status==2 ){ return Clear(false); }

  header-> SetRunNumber(confMan->GetRunNumber());
  header-> SetEventNumber(Event_Number);
  header->Convert( tko, confMan );
  if( header->IsTrig(Trig_Cosmic) || header->IsTrig(Trig_Reject) ) return Clear();

  blMan-> Convert(tko, confMan);
  cdsMan-> Convert(tko, confMan);

  bltrackMan-> DoTracking(blMan, confMan, true, true); //(blMan, confMan, global, timing);

  return Clear(UAna());
}

void EventAnalysisTemplate::Initialize(ConfMan *conf)
{
  confMan = conf;
  
  rtFile     = new TFile(confMan->GetOutFileName().c_str(), "recreate");
  header     = new EventHeader();
  blMan      = new BeamLineHitMan();
  bltrackMan = new BeamLineTrackMan();
  beamSpec   = new BeamSpectrometer(); 
  cdsMan     = new CDSHitMan();

  scaTree = new TTree("ScalerTree", "ScalerTree");
  scaMan = new ScalerMan();
  scaTree-> Branch("ScalerMan", &scaMan);

  t0 = time(0);
  InitializeHistogram();
}

bool EventAnalysisTemplate::Clear(bool flag)
{
  header-> Clear();
  blMan-> Clear();
  cdsMan-> Clear();
  bltrackMan-> Clear();
  beamSpec-> Clear();

  return flag;
}

void EventAnalysisTemplate::USca(int nsca, unsigned int *sca)
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

void EventAnalysisTemplate::Finalize()
{
  delete blMan;
  delete bltrackMan;
  delete beamSpec;
  delete cdsMan;
  delete scaMan;
  delete header;

  rtFile-> cd();
  rtFile-> Write();
  rtFile-> Close();
}
