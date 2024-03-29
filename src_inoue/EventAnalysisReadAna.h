#include "EventTemp.h"
#include "EventAlloc.h"

#include "EventHeader.h"
#include "ScalerMan.h"
#include "BeamSpectrometer.h"
#include "TrackTools.h"

#include "Tools.h"
#include "MyTools.h"
#include "AnaInfo.h"
#include "MyHistReadCDS.h"
#include "MyHistReadNC.h"
#include "MyHistReadFC.h"
#include "MyHistReadKNpim.h"
#include "MyHistReadKNpip.h"
#include "MyHistReadKNkm.h"
#include "MyHistReadKNpipi.h"
#include "MyHistReadKNpipi_D2.h"

#include "MyHistReadKPpimpim_D2.h"
#include "MyHistReadCDS_Lpim_D2.h"

#include "MyHistReadKNpipi_H2.h"

#include "MyHistReadVtxDEF.h"

class EventAnalysisReadAna : public EventTemp
{
public:
  EventAnalysisReadAna();
  ~EventAnalysisReadAna(){};
private:
  int t0, t1;
  int CDC_Event_Number;
  TFile *cdsFile;
  TTree *cdsTree;
  EventHeader *cdsHeader;
  CDSTrackingMan *cdstrackMan;

  int Ana_Event_Number;
  TFile *anaFile;
  TTree *anaTree;
  EventHeader *anaHeader;
  AnaInfo *anaInfo;
  bool IsMTDCFile;

  TFile *rtFile;
  //  TTree *scaTree;
  EventHeader *header;
  ScalerMan *scaMan;

  BeamLineHitMan *blMan;
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

EventAnalysisReadAna::EventAnalysisReadAna() : EventTemp()
{
}

bool EventAnalysisReadAna::UAna(TKOHitCollection *tko)
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
  header->Convert( tko, confMan );
  if( header->IsTrig(Trig_Cosmic) || header->IsTrig(Trig_Reject) ) return Clear();

  blMan-> Convert(tko, confMan);
  cdsMan-> Convert(tko, confMan);

  if( !anaFileMatching() ){
    std::cout<<" !!! AnaInfo file matching fault !!!"<<std::endl;
    return Clear(false);
  }

  return Clear(UAna());
}

void EventAnalysisReadAna::Initialize(ConfMan *conf)
{
  confMan = conf;
  ProtonArm::Initialize(conf, 0.9725);

  cdsFile = new TFile(confMan->GetCDSTrackFileName().c_str());
  cdstrackMan = new CDSTrackingMan();
  cdsHeader = new EventHeader();
  CDC_Event_Number=0;
  cdsTree = (TTree*)cdsFile-> Get("EventTree");
  cdsTree-> SetBranchAddress("EventHeader", &cdsHeader);  
  cdsTree-> SetBranchAddress("CDSTrackingMan", &cdstrackMan);
  
  anaFile = new TFile(confMan->GetMTDCTrackFileName().c_str());
  anaInfo = new AnaInfo();
  anaHeader = new EventHeader();
  Ana_Event_Number=0;
  anaTree = (TTree*)anaFile-> Get("EventTree");
  anaTree-> SetBranchAddress("EventHeader", &anaHeader);  
  anaTree-> SetBranchAddress("AnaInfo", &anaInfo);
  
  rtFile = new TFile(confMan->GetOutFileName().c_str(), "recreate");
  header = new EventHeader();
  blMan = new BeamLineHitMan();
  cdsMan =  new CDSHitMan();

  //  scaTree = new TTree("ScalerTree", "ScalerTree");
  scaMan = new ScalerMan();
  //  scaTree-> Branch("ScalerMan", &scaMan);

  t0 = time(0);

  InitializeHistogram();
}

bool EventAnalysisReadAna::cdsFileMatching()
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

bool EventAnalysisReadAna::anaFileMatching()
{
  anaTree -> GetEntry(Ana_Event_Number);
  if(Ana_Event_Number>=anaTree->GetEntries()){
    anaHeader->Clear();
    anaInfo->Clear();
    return false;
  }
  while( anaHeader->ev()<Event_Number ){
    Ana_Event_Number++;
    if(Ana_Event_Number>anaTree->GetEntries()){
      anaHeader->Clear();
      anaInfo->Clear();
      return false;
    }
    anaTree -> GetEntry(Ana_Event_Number);
  }
  if(anaHeader->ev()>Event_Number){
    anaHeader->Clear();
    anaInfo->Clear();
    return false;
  }
  return true;
}

bool EventAnalysisReadAna::Clear(bool flag)
{
  header-> Clear();
  blMan-> Clear();
  cdsMan-> Clear();

  cdsHeader-> Clear();
  cdstrackMan-> Clear();

  anaHeader-> Clear();
  anaInfo-> Clear();
  return flag;
}

void EventAnalysisReadAna::USca(int nsca, unsigned int *sca)
{
  Block_Event_Number++;
  header->SetBlockEventNumber( Block_Event_Number );
  scaMan->SetBlockEventNumber( Block_Event_Number );
  for( int i=0; i<nsca; i++ ){
    scaMan->AddHit( sca[i], confMan->GetScalerMapManager()->GetName(i) );
  }
  //  scaTree->Fill();
  scaMan->Clear();
}

void EventAnalysisReadAna::Finalize()
{
  delete blMan;
  delete cdsMan;
  delete cdsHeader;
  delete cdstrackMan;
  delete anaHeader;
  delete anaInfo;
  delete cdsFile;
  delete anaFile;

  rtFile-> cd();
  rtFile-> Write();
  rtFile-> Close();
}
