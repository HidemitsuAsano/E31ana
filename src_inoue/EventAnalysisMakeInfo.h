#include "EventTemp.h"
#include "EventAlloc.h"

#include "EventHeader.h"
#include "ScalerMan.h"
#include "BeamSpectrometer.h"
#include "TrackTools.h"

#include "MyTools.h"
#include "AnaInfo.h"

class EventAnalysisMakeInfo : public EventTemp
{
public:
  EventAnalysisMakeInfo();
  ~EventAnalysisMakeInfo(){};
private:
  int t0, t1;
  int Num_Cosmic;
  int Num_Reject;
  int CDC_Event_Number;
  TFile *cdsFile;
  TTree *cdsTree;
  EventHeader *cdsHeader;
  CDSTrackingMan *cdstrackMan;

  TFile *rtFile;
  TTree *evTree;
  TTree *scaTree;
  EventHeader *header;
  ScalerMan *scaMan;

  BeamLineHitMan *blMan;
  CDSHitMan *cdsMan;
  BeamLineTrackMan *bltrackMan;
  BeamSpectrometer *beamSpec;
  AnaInfo *anaInfo;

public:
  void Initialize(ConfMan *conf);
  void USca(int usca, unsigned int *sca);
  bool UAna(TKOHitCollection *tko);
  void Finalize();

private:
  bool cdsFileMatching();
  bool Clear(bool flag=true);
};

EventAnalysisMakeInfo::EventAnalysisMakeInfo() : EventTemp()
{
}

void EventAnalysisMakeInfo::Initialize(ConfMan *conf)
{
  confMan = conf;
  ProtonArm::Initialize(conf);

  cdsFile = new TFile(confMan->GetCDSTrackFileName().c_str());
  cdstrackMan = new CDSTrackingMan();
  cdsHeader = new EventHeader();
  CDC_Event_Number=0;
  cdsTree = (TTree*)cdsFile-> Get("EventTree");
  cdsTree-> SetBranchAddress("EventHeader", &cdsHeader);  
  cdsTree-> SetBranchAddress("CDSTrackingMan", &cdstrackMan);
  
  rtFile = new TFile(confMan->GetOutFileName().c_str(), "recreate");
  evTree = new TTree("EventTree", "EventTree");
  header = new EventHeader();
  blMan = new BeamLineHitMan();
  bltrackMan = new BeamLineTrackMan();
  beamSpec = new BeamSpectrometer(confMan);
  cdsMan =  new CDSHitMan();
  anaInfo = new AnaInfo();
  evTree-> Branch("EventHeader", &header);
  evTree-> Branch("AnaInfo", &anaInfo);

  scaTree = new TTree("ScalerTree", "ScalerTree");
  scaMan = new ScalerMan();
  scaTree-> Branch("ScalerMan", &scaMan);

  Num_Cosmic=0;
  Num_Reject=0;

  t0 = time(0);
}

bool EventAnalysisMakeInfo::cdsFileMatching()
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
  return true;
}

bool EventAnalysisMakeInfo::Clear(bool flag)
{
  header-> Clear();
  blMan-> Clear();
  bltrackMan-> Clear();
  beamSpec-> Clear();
  cdsMan-> Clear();

  cdsHeader-> Clear();
  cdstrackMan-> Clear();

  anaInfo-> Clear();
  return flag;
}

void EventAnalysisMakeInfo::USca(int nsca, unsigned int *sca)
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

void EventAnalysisMakeInfo::Finalize()
{
  std::cout<<"===== EventAnalysisMakeInfo::Finalize ====="<<std::endl;
  std::cout<<"> Event_Number   : "<<Event_Number<<std::endl;
  std::cout<<"> evTree Entries : "<<evTree->GetEntries()<<std::endl;
  std::cout<<">       Cosmic   : "<<Num_Cosmic<<std::endl;
  std::cout<<">       Reject   : "<<Num_Reject<<std::endl;

  rtFile-> cd();
  gFile-> Write();
  gFile-> Close();

  delete blMan;
  delete cdsMan;
  delete bltrackMan;
  delete beamSpec;
  delete cdstrackMan;
  delete anaInfo;
}
