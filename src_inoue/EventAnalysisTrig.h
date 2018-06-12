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

#include "MyHistTrig.h"

#include "MyHistT0.h"
#include "MyHistBHD.h"
#include "MyHistPC.h"
#include "MyHistCVC.h"
#include "MyHistNC.h"
#include "MyHistBPD.h"

#include "MyHistBLC1a.h"
#include "MyHistBLC1b.h"
#include "MyHistBLC2a.h"
#include "MyHistBLC2b.h"
#include "MyHistBPC.h"
#include "MyHistFDC1.h"
#include "MyHistCDC.h"
#include "MyHistBLC1.h"
#include "MyHistBLC2.h"
#include "MyHistReduction.h"

#include "MyHistBHDT0.h"

class EventAnalysisSca : public EventTemp
{
 public:
  EventAnalysisSca();
  ~EventAnalysisSca(){};

 private:
  int t0, t1;

  TFile *rtFile;
  TTree *evTree;
  TTree *scaTree;
  EventHeader *header;
  ScalerMan *scaMan;

 public:
  void Initialize(ConfMan *conf);
  void USca(int nsca, unsigned int *sca);
  bool UAna(TKOHitCollection *tko);
  void Finalize();

  bool Clear(bool flag=true);
};

EventAnalysisSca::EventAnalysisSca() : EventTemp()
{
}

void EventAnalysisSca::Initialize(ConfMan *conf)
{
  confMan = conf;

  rtFile = new TFile(confMan->GetOutFileName().c_str(), "recreate");
  scaTree = new TTree("ScalerTree", "ScalerTree");
  scaMan = new ScalerMan();
  scaTree->Branch("ScalerMan", &scaMan);

  header = new EventHeader();
  evTree = new TTree("EventTree", "EventTree");
  evTree->Branch("EventHeader", &header);

  initHistTrig();

  t0 = time(0);
}

void EventAnalysisSca::USca(int nsca, unsigned int *sca)
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

bool EventAnalysisSca::UAna(TKOHitCollection *tko)
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
  header-> Convert( tko, confMan );

  fillTrig(header);

  evTree-> Fill();
  Clear();
}


void EventAnalysisSca::Finalize()
{
  gFile-> cd();
  gFile-> Write();
  gFile-> Close();

  delete header;
  delete scaMan;
}

bool EventAnalysisSca::Clear(bool flag)
{
  header-> Clear();

  return flag;
}
