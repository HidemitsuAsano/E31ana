#include "EventTemp.h"
#include "EventAlloc.h"

#include "EventHeader.h"
#include "ScalerMan.h"
#include "BeamSpectrometer.h"
#include "TrackTools.h"

#include "MyTools.h"
#include "AnaInfo.h"
#include "MyHistReadBeam.h"
#include "MyHistReadCDS.h"
#include "MyHistCDS_pipi.h"
#include "MyHistCDS_ppim.h"
#include "MyHistReadNC.h" 
#include "MyHistKN_pipi.h"

class EventAnalysisMakeInfo : public EventTemp
{
public:
  EventAnalysisMakeInfo();
  ~EventAnalysisMakeInfo(){};
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

  TFile *rtFile;
  TTree *scaTree;
  EventHeader *header;
  ScalerMan *scaMan;

  BeamLineHitMan *blMan;
  CDSHitMan *cdsMan;

public:
  void Initialize(ConfMan *conf);
  void USca(int usca, unsigned int *sca);
  bool UAna(TKOHitCollection *tko);
  void Finalize();

private:
  void InitializeHistogram();
  bool cdsFileMatching();
  bool anaFileMatching();
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

  anaFile = new TFile(confMan->GetMTDCTrackFileName().c_str());
  anaInfo = new AnaInfo();
  anaHeader = new EventHeader();
  CDC_Event_Number=0;
  anaTree = (TTree*)anaFile-> Get("EventTree");
  anaTree-> SetBranchAddress("EventHeader", &anaHeader);  
  anaTree-> SetBranchAddress("AnaInfo", &anaInfo);
  
  rtFile = new TFile(confMan->GetOutFileName().c_str(), "recreate");
  header = new EventHeader();
  blMan = new BeamLineHitMan();
  cdsMan =  new CDSHitMan();

  scaTree = new TTree("ScalerTree", "ScalerTree");
  scaMan = new ScalerMan();
  scaTree-> Branch("ScalerMan", &scaMan);

  t0 = time(0);

  InitializeHistogram();
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

bool EventAnalysisMakeInfo::anaFileMatching()
{
  anaTree -> GetEntry(Ana_Event_Number);
  std::cout<<"aaa"<<Event_Number<<"  "<<Ana_Event_Number<<"  "<<anaHeader->ev()<<std::endl;
  if(Ana_Event_Number>=anaTree->GetEntries()){
    anaHeader->Clear();
    anaInfo->Clear();
    std::cout<<"bbb"<<std::endl;
    return false;
  }
  while( anaHeader->ev()<Event_Number ){
    Ana_Event_Number++;
    if(Ana_Event_Number>anaTree->GetEntries()){
      anaHeader->Clear();
      anaInfo->Clear();
      std::cout<<"ccc"<<std::endl;
      return false;
    }
    anaTree -> GetEntry(Ana_Event_Number);
  }
  if(anaHeader->ev()>Event_Number){
    std::cout<<"ddd"<<anaHeader->ev()<<std::endl;
    anaHeader->Clear();
    anaInfo->Clear();
    return false;
  }
  return true;
}

/* bool EventAnalysisMakeInfo::anaFileMatching() */
/* { */
/*   anaTree -> GetEntry(Ana_Event_Number); */
/*   if(Ana_Event_Number>=anaTree->GetEntries()){ */
/*     anaHeader->Clear(); */
/*     anaInfo->Clear(); */
/*     return false; */
/*   } */
/*   while( anaHeader->ev()<Event_Number ){ */
/*     Ana_Event_Number++; */
/*     if(Ana_Event_Number>anaTree->GetEntries()){ */
/*       anaHeader->Clear(); */
/*       anaInfo->Clear(); */
/*       return false; */
/*     } */
/*     anaTree -> GetEntry(Ana_Event_Number); */
/*   } */
/*   if(anaHeader->ev()>Event_Number){ */
/*     anaHeader->Clear(); */
/*     anaInfo->Clear(); */
/*     return false; */
/*   } */
/*   return true; */
/* } */

bool EventAnalysisMakeInfo::Clear(bool flag)
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
  rtFile-> cd();
  gFile-> Write();
  gFile-> Close();

  delete blMan;
  delete cdsMan;
  delete cdstrackMan;
  delete anaInfo;
}
