#include "File.h"
#include "ConfMan.h"
#include "TKO.h"
#include "CDSHitMan.h"
#include "BeamLineHitMan.h"
#include "EventHeader.h"
#include "ScalerMan.h"
#include "BeamLineTrackMan.h"
#include "HistManBeamAna.h"
#include "AnaData.h"

#include "EventAlloc.h"
#include "EventTemp.h"
#include "HistMan.h"


class EventAnalysisAnaData : public EventTemp
{
public:
  EventAnalysisAnaData();
  ~EventAnalysisAnaData();

private:
  int t0, t1;

  int CDC_Event_Num;
  TFile *cdcFile;
  TTree *cdcTree;
  CDSTrackingMan *cdstrackMan;
  EventHeader *cdcHeader;

  TFile *rtFile;
  TTree *scaTree;
  ScalerMan *scaMan;

  TTree *evTree;
  CDSHitMan *cdsMan;
  BeamLineHitMan *blMan;
  EventHeader *header;
  AnaData *anaData;
  HistMan *histMan;

  BeamLineTrackMan *bltrackMan; 
  BeamSpectrometer *D5;

public:
  void Initialize(ConfMan *conf);
  void USca(int nsca, unsigned int *sca);
  bool UAna(TKOHitCollection *tko);
  void Clear();
  void Finalize();
};

EventAnalysisAnaData::EventAnalysisAnaData() : EventTemp()
{
}

EventAnalysisAnaData::~EventAnalysisAnaData()
{
}

void EventAnalysisAnaData::Initialize(ConfMan *conf)
{
  confMan = conf;
  CDC_Event_Num=0;
  cdcFile = new TFile((TString)confMan->GetCDSTrackFileName());
  cdstrackMan = new CDSTrackingMan();
  cdcHeader = new EventHeader();
  cdcTree = (TTree*)cdcFile-> Get("EventTree");
  cdcTree-> SetBranchAddress("EventHeader", &cdcHeader);
  cdcTree-> SetBranchAddress("CDSTrackingMan", &cdstrackMan);

  rtFile = new TFile(confMan-> GetOutFileName().c_str(), "recreate");
  scaMan = new ScalerMan();
  scaTree = new TTree("ScalerTree", "ScalerTree");
  scaTree-> Branch("ScalerMan", &scaMan);

  header = new EventHeader();
  cdsMan = new CDSHitMan();
  blMan = new BeamLineHitMan();
  anaData = new AnaData();
  evTree = new TTree("EventTree", "EventTree");
  evTree-> Branch("EventHeader", &header);
  evTree-> Branch("AnaData", &anaData);

  bltrackMan = new BeamLineTrackMan();
  D5 = new BeamSpectrometer(confMan);
  histMan = new HistMan(rtFile);

  t0 = time(0);
}

void EventAnalysisAnaData::USca(int nsca, unsigned int *sca)
{
  Block_Event_Number++;
  header-> SetBlockEventNumber(Block_Event_Number);
  scaMan-> SetBlockEventNumber(Block_Event_Number);
  for( int i=0; i<nsca; i++ ){
    scaMan-> AddHit( sca[i], confMan->GetScalerMapManager()->GetName(i) );
  }
  scaTree-> Fill();
  scaMan-> Clear();
}

bool EventAnalysisAnaData::UAna(TKOHitCollection *tko)
{
  Event_Number++;
  header-> SetRunNumber(confMan->GetRunNumber());
  header-> SetEventNumber(Event_Number);

  if( Event_Number%5000==0 ){
    t1 = time(0);
    std::cout<<" Event# : "<<Event_Number<<" BlockEvent# : "<<Block_Event_Number<<" Time(s) : "<<(t1-t0)<<std::endl;
  }

  histMan-> fill(header, confMan);

  header-> Convert(tko, confMan);
  if( header->IsTrig(Trig_Cosmic) ){
    Clear();
    return true;
  }

  blMan-> Convert(tko, confMan);
  cdsMan-> Convert(tko, confMan);
  bltrackMan-> DoTracking(blMan, confMan, true, false);
  double D5mom = DBL_MIN;
  if( bltrackMan->ntrackBLC1()==1 && bltrackMan->ntrackBLC2()==1 ){
    D5-> TMinuitFit(bltrackMan->trackBLC1(0), bltrackMan->trackBLC2(0), confMan);
    D5mom = D5->mom();
  }

  histMan->set(blMan);

  if( CDC_Event_Num>cdcTree->GetEntries() ){
    Clear();
    return true;
  }
  cdcTree-> GetEntry(CDC_Event_Num);
  while( cdcHeader->ev()<Event_Number ){
    CDC_Event_Num++;
    if( CDC_Event_Num>cdcTree->GetEntries() ){
      Clear();
      return true;
    }
    cdcTree-> GetEntry(CDC_Event_Num);
  }
  if( cdcHeader->ev()>Event_Number ){
    Clear();
    return true;
  }
  if( cdcHeader->ev()!=Event_Number ){
    std::cout<<"  !!! CDC File Event matching miss !!!"<<std::endl;
    Clear();
    return false;
  }

  if( D5mom>0.0 ){
    anaData-> set(header, confMan, D5mom, blMan, bltrackMan, cdsMan, cdstrackMan);
    evTree-> Fill();
  }

  Clear();
  return true;
}

void EventAnalysisAnaData::Finalize()
{
  cdcFile-> Close();

  rtFile-> cd();
  rtFile-> Write();
  rtFile-> Close();
  confMan-> SaveCode();
  confMan-> SaveParams();

  delete blMan;
  delete cdsMan;
  delete header;
  delete bltrackMan;
  delete cdstrackMan;
  delete scaMan;
  delete D5;

  delete cdcHeader;
}

void EventAnalysisAnaData::Clear()
{
  cdcHeader-> Clear();
  cdstrackMan-> Clear();

  header-> Clear();
  cdsMan-> Clear();
  blMan-> Clear();
  bltrackMan-> Clear();
  anaData-> Clear();
  D5-> Clear();

  histMan-> clear();
}
				

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisAnaData *event = new EventAnalysisAnaData();
  return (EventTemp*)event;
}
