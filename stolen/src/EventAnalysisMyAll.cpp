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
#include "MyAnalysisFWDCharged.h"
#include "MyAnalysisK0.h"

#define Debug 0

class EventAnalysisMyAll: public EventTemp
{
  public:
    EventAnalysisMyAll();
    ~EventAnalysisMyAll();
  private:
    TFile *rtFile;
    TFile *rtFile2;
    TTree *evTree;
    TTree *scaTree;
    TFile *cdcFile;
    TTree *cdcTree;
    CDSHitMan *cdsMan;
    CDSTrackingMan *cdstrackMan;
    EventHeader *cdsheader;
    BeamLineHitMan *blMan;
    BeamLineTrackMan *bltrackMan;
    EventHeader *header;
    Particle* particle;
    ScalerMan *scaMan;
    MyAnalysisBL* blAna;
    MyAnalysisCDS* cdsAna;
    MyAnalysisFWDNeutral* fwdnAna;
    MyAnalysisFWDCharged* fwdcAna;
    //MyAnalysisK0* k0Ana;

    int AllGoodTrack;
    int nAllTrack;
    int nTrack;
    int CDC_Event_Number;
    bool trackflag;
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

  particle = new Particle();

  rtFile->cd();
  evTree = new TTree("EventTree","EventTree");
  evTree -> Branch("EventHeader",&header);
  evTree -> Branch("Particle",&particle);

  std::string ofname = confMan->GetOutFileName();
  ofname.insert(ofname.find(".root"),"_anaAll");
  rtFile2 =  new TFile( Form("%s",ofname.c_str()), "RECREATE");
  rtFile2->cd();
  blAna = new MyAnalysisBL(rtFile2, confMan);
  cdsAna  = new MyAnalysisCDS(rtFile2, confMan);
  fwdnAna = new MyAnalysisFWDNeutral(rtFile2, confMan);
  fwdcAna = new MyAnalysisFWDCharged(rtFile2, confMan);

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
  trackflag = false;
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

  if( Event_Number%10000==0 )
  {
    if(cdcFile!=0){
      t1=time(0);
      std::cout << " Event# : " << Event_Number << "  BlockEvent# : " << Block_Event_Number << " AllTrack# : " << nAllTrack <<" GoodTrack# : " << AllGoodTrack  <<"  "<<CDC_Event_Number<<" / "<<cdcTree->GetEntries() << " Time (s):" << (t1-t0) << std::endl;
    }
  }

  if(cdcFile!=0){
    if(CDC_Event_Number>=cdcTree->GetEntries()){
      cdsheader=0;
      cdstrackMan=0;
      cdcFile->Close();
      cdcFile=0;
    }
    while( cdsheader->ev()<Event_Number ){
      CDC_Event_Number++;
      if(CDC_Event_Number>cdcTree->GetEntries()){
        cdsheader=0;
        cdstrackMan=0;
        cdcFile->Close();
        cdcFile=0;
        break;
      }
      cdcTree -> GetEntry(CDC_Event_Number);
    }
  }

  if(cdsheader){
    if(cdsheader->ev()>Event_Number){
      trackflag = false;
    }
    else{
      trackflag = true;
    }
  }
  if(trackflag){
    int nGoodTrack = cdstrackMan -> nGoodTrack();
    int nAllTrack = cdstrackMan -> nTrack();
    AllGoodTrack += nGoodTrack;
    nTrack += nAllTrack;
  }

  header->SetRunNumber(0);
  header->SetEventNumber(Event_Number);

  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );
  if(trackflag){
    cdstrackMan -> Calc(cdsMan,confMan);
  }

  rtFile2->cd();
  if(!blAna->DoAnalysis(confMan, header, blMan, bltrackMan, particle)){
    Clear();
    return true;
  }
  if(trackflag){
    if(!cdsAna->DoAnalysis(confMan, header, blMan, bltrackMan, cdsMan, cdstrackMan, particle)){
      Clear();
      return true;
    }
    if(!fwdnAna->DoAnalysis(confMan, header, blMan, bltrackMan, cdsMan, cdstrackMan, particle)){
      Clear();
      return true;
    }
    if(!fwdcAna->DoAnalysis(confMan, header, blMan, bltrackMan, cdsMan, cdstrackMan, particle)){
      Clear();
      return true;
    }
  }

  rtFile->cd();
  evTree->Fill();

  Clear();
  return true;
}

void EventAnalysisMyAll::Clear()
{
  particle->Clear();
  blAna->Clear();
  cdsAna->Clear();
  fwdnAna->Clear();
  fwdcAna->Clear();
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
  rtFile2->cd();
  gFile->Write();
  gFile->Close();

  delete particle;
  delete blAna;
  delete cdsAna;
  delete fwdnAna;
  delete fwdcAna;
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
