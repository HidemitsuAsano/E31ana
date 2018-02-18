#include "File.h"
#include "ConfMan.h"
#include "TKO.h"
#include "CDSHitMan.h"
#include "BeamLineHitMan.h"
#include "BeamLineTrackMan.h"
#include "EventHeader.h"
#include "ScalerMan.h"
#include "Particle.h"

#include "EventAlloc.h"
#include "EventTemp.h"
#include "Display.h" 
#include "ELossTools.h"
#include "TrackTools.h"
#include "MyAnalysisBL.h"
#include "MyAnalysisCDS.h"
#include "MyAnalysisFWDNeutral.h"

#define Debug 0

class EventAnalysisMyAll: public EventTemp
{
  public:
    EventAnalysisMyAll();
    ~EventAnalysisMyAll();
  private:
    TFile *rtFile;
    TTree *evTree;
    TTree *scaTree;
    TFile *cdcFile;
    TTree *cdcTree;
    TFile *mtdcFile;
    TTree *mtdcTree;
    CDSHitMan *cdsMan;
    CDSTrackingMan *cdstrackMan;
    EventHeader *cdsheader;
    BeamLineHitMan *blMan;
    BeamLineTrackMan *bltrackMan;
    EventHeader *header;
    EventHeader *mtdcheader;
    Particle* particle;
    ScalerMan *scaMan;
    MyAnalysisBL* blAna;
    MyAnalysisCDS* cdsAna;
    MyAnalysisFWDNeutral* fwdnAna;

    int AllGoodTrack;
    int nAllTrack;
    int CDC_Event_Number;
    int MTDC_Event_Number;
    int t0, t1;

  public:
    void Initialize( ConfMan *conf );
    void USca( int nsca, unsigned int *sca );
    bool UAna( TKOHitCollection *tko );
    void Finalize();
    void Clear();

};

  EventAnalysisMyAll::EventAnalysisMyAll()
: EventTemp()
{
}

EventAnalysisMyAll::~EventAnalysisMyAll()
{
}

const int MaxTreeSize = 1900000000000;
void EventAnalysisMyAll::Initialize( ConfMan *conf )
{
#if 1
  std::cout << " Enter EventAnalysisMyAll::Initialize " << std::endl;
#endif
  confMan = conf;


  TString cdcfname=confMan->GetCDSTrackFileName();
  cdcFile = new TFile(cdcfname);
  if(!cdcFile->IsOpen()){
    std::cout<<" failed to open " <<cdcfname<< "  !!!"<<std::endl;
    exit(false);
  }
  std::cout<<" CDC Track File Successfully opend !!! " <<cdcfname<<std::endl;
  cdstrackMan   = new CDSTrackingMan();
  cdsheader    = new EventHeader();
  cdcTree=(TTree*)cdcFile->Get("EventTree");
  cdcTree->SetBranchAddress( "CDSTrackingMan", &cdstrackMan );
  cdcTree->SetBranchAddress( "EventHeader" ,&cdsheader);

  TString mtdcfname=confMan->GetMTDCTrackFileName();
  mtdcFile = new TFile(mtdcfname);
  if(!mtdcFile->IsOpen()){
    std::cout<<" failed to open " <<mtdcfname<< "  !!!"<<std::endl;
    exit(false);
  }
  std::cout<<" MTDC Track File Successfully opend !!! " <<mtdcfname<<std::endl;
  particle   = new Particle();
  mtdcheader    = new EventHeader();
  mtdcTree=(TTree*)mtdcFile->Get("EventTree");
  mtdcTree->SetBranchAddress( "Particle", &particle );
  mtdcTree->SetBranchAddress( "EventHeader" ,&mtdcheader);

  rtFile = new TFile( confMan->GetOutFileName().c_str(), "recreate" );
  header = new EventHeader();
  if( header==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  cdsMan = new CDSHitMan();
  if( cdsMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  blMan = new BeamLineHitMan();
  if( blMan ==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  bltrackMan = new BeamLineTrackMan();
  if( bltrackMan ==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  
  scaMan = new ScalerMan();
  if( scaMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }

  //blAna = new MyAnalysisBL(rtFile, confMan);
  cdsAna  = new MyAnalysisCDS(rtFile, confMan);
  fwdnAna = new MyAnalysisFWDNeutral(rtFile, confMan);
  //particle = new Particle();

  //rtFile->cd();
  //evTree = new TTree("EventTree","EventTree");
  //evTree -> Branch("EventHeader",&header);
  //evTree -> Branch("Particle",&particle);

  AllGoodTrack = 0;
  nAllTrack = 0;
  CDC_Event_Number = 0;
  MTDC_Event_Number = 0;
  t0=time(0);
}

void EventAnalysisMyAll::USca( int nsca, unsigned int *sca )
{
  Block_Event_Number++;
  header->SetBlockEventNumber( Block_Event_Number );
  scaMan->SetBlockEventNumber( Block_Event_Number );
  for( int i=0; i<nsca; i++ ){
    scaMan->AddHit( sca[i], confMan->GetScalerMapManager()->GetName(i) );
  }

  scaMan->Clear();
}

bool EventAnalysisMyAll::UAna( TKOHitCollection *tko )
{
  Event_Number++;
  { int status = confMan->CheckEvNum( Event_Number, Block_Event_Number );
    if( status==1 ){
      Clear();
      return true;
    }
    if( status==2 ){
      Clear();
      return false;
    }
  }

  if( Event_Number%1000==0 )
  {
    t1=time(0);
    std::cout << " Event# : " << Event_Number << "  BlockEvent# : " << Block_Event_Number << " AllTrack# : " << nAllTrack <<" GoodTrack# : " << AllGoodTrack  <<"  "<<CDC_Event_Number<<" / "<<cdcTree->GetEntries() << " Time (s):" << (t1-t0) << std::endl;
  }

  if(CDC_Event_Number>=cdcTree->GetEntries()){
    Clear();
    return false;  
  }
  cdcTree->GetEntry(CDC_Event_Number);
  if(cdsheader->ev()!=Event_Number){
    int tmpnum=CDC_Event_Number;
    while(cdsheader->ev()<Event_Number&&tmpnum<=cdcTree->GetEntries()){
      tmpnum++;    
      cdcTree->GetEntry(tmpnum);
    }
    if(cdsheader->ev()!=Event_Number){
      Clear();
      return true;
    }
    else CDC_Event_Number=tmpnum;
  }
  CDC_Event_Number++;
  int nGoodTrack = cdstrackMan -> nGoodTrack();
  int nAllTrack = cdstrackMan -> nTrack();
  AllGoodTrack += nGoodTrack;
  nAllTrack += nAllTrack;

  if(MTDC_Event_Number>=mtdcTree->GetEntries()){
    Clear();
    return false;  
  }
  mtdcTree->GetEntry(MTDC_Event_Number);
  if(mtdcheader->ev()!=Event_Number){
    int tmpnum=MTDC_Event_Number;
    while(mtdcheader->ev()<Event_Number&&tmpnum<=mtdcTree->GetEntries()){
      tmpnum++;    
      mtdcTree->GetEntry(tmpnum);
    }
    if(mtdcheader->ev()!=Event_Number){
      Clear();
      return true;
    }
    else MTDC_Event_Number=tmpnum;
  }
  MTDC_Event_Number++;

  header->SetRunNumber(0);
  header->SetEventNumber(Event_Number);

  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );
  cdstrackMan -> Calc( cdsMan, confMan );

  rtFile->cd();

  //if(!blAna->DoAnalysis(confMan, header, blMan, bltrackMan, particle)){
  //  Clear();
  //  return true;
  //}
  if(!cdsAna->DoAnalysis(confMan, header, blMan, bltrackMan, cdsMan, cdstrackMan, particle)){
    Clear();
    return true;
  }
  if(!fwdnAna->DoAnalysis(confMan, header, blMan, bltrackMan, cdsMan, cdstrackMan, particle)){
    Clear();
    return true;
  }

  //evTree->Fill();

  Clear();
  return true;
}

void EventAnalysisMyAll::Clear()
{
  //particle->Clear();
  //blAna->Clear();
  cdsAna->Clear();
  fwdnAna->Clear();
  blMan->Clear();
  bltrackMan->Clear();
  cdsMan->Clear();
  header->Clear();
}

void EventAnalysisMyAll::Finalize()
{
  std::cout << " Enter EventAnalysisMyAll::Finalize " << std::endl;

  rtFile->cd();
  gFile->Write();
  gFile->Close();

  //delete particle;
  //delete blAna;
  delete cdsAna;
  delete fwdnAna;
  delete blMan;
  delete bltrackMan;
  delete cdsMan;
  delete header;
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisMyAll *event = new EventAnalysisMyAll();
  return (EventTemp*)event;
}
