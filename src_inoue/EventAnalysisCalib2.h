#include "File.h"
#include "ConfMan.h"
#include "TKO.h"
#include "CDSHitMan.h"
#include "BeamLineHitMan.h"
#include "EventHeader.h"
#include "ScalerMan.h"

#include "EventAlloc.h"
#include "EventTemp.h"
#include "Display.h"

#include "BeamLineTrackMan.h"
#include "BeamSpectrometer.h"

#include "MyHistT0.h"
#include "MyHistBHD.h"
#include "MyHistPC.h"
#include "MyHistCVC.h"
#include "MyHistNC.h"
#include "MyHistBPD.h"
#include "MyHistDEF.h"
#include "MyHistBVC.h"

#include "MyHistBLC1a.h"
#include "MyHistBLC1b.h"
#include "MyHistBLC2a.h"
#include "MyHistBLC2b.h"
#include "MyHistBPC.h"
#include "MyHistFDC1.h"
#include "MyHistCDC.h"
#include "MyHistBLC1.h"
#include "MyHistBLC2.h"
#include "MyHistCDH.h"
#include "MyHistReduction.h"

#include "MyHistBHDT0.h"
#include "MyHistBeamAna.h"
#include "MyHistCDS.h"

class EventAnalysisCalib2 : public EventTemp
{
 public:
  EventAnalysisCalib2();
  ~EventAnalysisCalib2(){};

 private:
  int t0, t1;

  TFile *rtFile;
  TTree *evTree;
  TTree *scaTree;
  EventHeader *header;
  ScalerMan *scaMan;

  TFile *cdsFile;
  int cdsTreeNum;
  TTree *cdsTree;
  EventHeader *cdsHeader;
  CDSTrackingMan *cdstrackMan;

  BeamLineHitMan *blMan;
  BeamLineTrackMan *bltrackMan;
  CDSHitMan *cdsMan;
  BeamSpectrometer *beamSpec;

 public:
  void Initialize(ConfMan *conf);
  void USca(int nsca, unsigned int *sca);
  bool UAna(TKOHitCollection *tko);
  void Finalize();

 private:
  void InitializeHistogram();
  void anaBL0();
  void anaBL1();
  void anaCDS();

  bool Clear(bool flag=true);
};

EventAnalysisCalib2::EventAnalysisCalib2() : EventTemp()
{
}

void EventAnalysisCalib2::Initialize(ConfMan *conf)
{
  confMan = conf;

  rtFile = new TFile(confMan->GetOutFileName().c_str(), "recreate");
  scaTree = new TTree("ScalerTree", "ScalerTree");
  scaMan = new ScalerMan();
  scaTree->Branch("ScalerMan", &scaMan);

  cdsTreeNum=0;
  cdsFile = new TFile((TString)confMan->GetCDSTrackFileName());
  cdstrackMan = new CDSTrackingMan();
  cdsHeader = new EventHeader();
  cdsTree = (TTree*)cdsFile-> Get("EventTree");
  cdsTree-> SetBranchAddress("EventHeader", &cdsHeader);
  cdsTree-> SetBranchAddress("CDSTrackingMan", &cdstrackMan);

  header = new EventHeader();
  blMan = new BeamLineHitMan();
  cdsMan = new CDSHitMan();
  bltrackMan = new BeamLineTrackMan();
  beamSpec = new BeamSpectrometer(conf);
  evTree = new TTree("EventTree", "EventTree");
  evTree->Branch("EventHeader", &header);

  InitializeHistogram();

  t0 = time(0);
}

void EventAnalysisCalib2::USca(int nsca, unsigned int *sca)
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

bool EventAnalysisCalib2::UAna(TKOHitCollection *tko)
{
  Event_Number++;
  if( Event_Number%1000==0 ){
    t1 = time(0);
    std::cout << " Event# : " << Event_Number << "  BlockEvent# : " << Block_Event_Number << " Time (s): "<<(t1-t0)<<std::endl;
  }

  int status = confMan->CheckEvNum( Event_Number, Block_Event_Number );
  if( status==1 ) return Clear();
  if( status==2 ) return Clear(false);

  header-> SetRunNumber(confMan->GetRunNumber());
  header-> SetEventNumber(Event_Number);
  header->Convert( tko, confMan );

  if( header->IsTrig(Trig_Cosmic) || header->IsTrig(Trig_Reject) ) return Clear();
  blMan-> Convert(tko, confMan);
  cdsMan-> Convert(tko, confMan);

  anaBL0();
  //  std::cout<<" DoTracking call"<<std::endl;
  bltrackMan-> DoTracking(blMan, confMan, true, true);

  LocalTrack *trackBLC1=MyTools::trackBLC1(bltrackMan);
  LocalTrack *trackBLC2=MyTools::trackBLC2(bltrackMan);
  if( trackBLC1 && trackBLC2 ) beamSpec-> TMinuitFit(trackBLC1, trackBLC2, confMan);

  anaBL1();

  if( cdsTreeNum<cdsTree->GetEntries() ){
    cdsTree->GetEntry(cdsTreeNum);
    while( cdsHeader->ev()<Event_Number ){
      cdsTreeNum++;
      if( cdsTreeNum>cdsTree->GetEntries() ){
	return Clear();
      }
      cdsTree->GetEntry(cdsTreeNum);
    }
  }

  /* std::cout<<Event_Number<<std::endl; */
  /* std::cout<<cdsHeader->ev()<<std::endl; */
  if( cdsHeader->ev()==Event_Number ){
    cdstrackMan->Calc(cdsMan, confMan, false);
    anaCDS();
  }


  Clear();
}


void EventAnalysisCalib2::Finalize()
{
  gFile-> cd();
  gFile-> Write();
  gFile-> Close();

  delete header;
  delete blMan;
  delete bltrackMan;
  delete beamSpec;
  delete cdsMan;
  delete cdstrackMan;
  delete scaMan;
}

bool EventAnalysisCalib2::Clear(bool flag)
{
  //  cdsTreeNum=-1;

  header-> Clear();
  blMan-> Clear();
  bltrackMan-> Clear();
  beamSpec-> Clear();
  cdsMan-> Clear();
  cdstrackMan-> Clear();

  return flag;
}
