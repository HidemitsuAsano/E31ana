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

#define Debug 0

class EventAnalysisMyReAll: public EventTemp
{
  public:
    EventAnalysisMyReAll();
    ~EventAnalysisMyReAll();
  private:
    TFile *rtFile;
    TFile *rtFile2;
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
    Particle* particle2;
    ScalerMan *scaMan;
    MyAnalysisBL* blAna;
    MyAnalysisCDS* cdsAna;
    MyAnalysisFWDNeutral* fwdnAna;
    MyAnalysisFWDCharged* fwdcAna;

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

  EventAnalysisMyReAll::EventAnalysisMyReAll()
: EventTemp()
{
}

EventAnalysisMyReAll::~EventAnalysisMyReAll()
{
}

const int MaxTreeSize = 1900000000000;
void EventAnalysisMyReAll::Initialize( ConfMan *conf )
{
#if 1
  std::cout << " Enter EventAnalysisMyReAll::Initialize " << std::endl;
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

  particle2 = new Particle();

  rtFile->cd();
  evTree = new TTree("EventTree","EventTree");
  evTree -> Branch("EventHeader",&header);
  evTree -> Branch("Particle",&particle2);

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

void EventAnalysisMyReAll::USca( int nsca, unsigned int *sca )
{
  Block_Event_Number++;
  header->SetBlockEventNumber( Block_Event_Number );
  scaMan->SetBlockEventNumber( Block_Event_Number );
  for( int i=0; i<nsca; i++ ){
    scaMan->AddHit( sca[i], confMan->GetScalerMapManager()->GetName(i) );
  }

  scaMan->Clear();
}

bool EventAnalysisMyReAll::UAna( TKOHitCollection *tko )
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


  if(MTDC_Event_Number>=mtdcTree->GetEntries()){
    Clear();
    return true;  
  }
  while( mtdcheader->ev()<Event_Number ){
    MTDC_Event_Number++;
    if(MTDC_Event_Number>mtdcTree->GetEntries()){
      Clear();  
      return true;
    }
    mtdcTree -> GetEntry(MTDC_Event_Number);
  }
  if(mtdcheader->ev()>Event_Number){
    Clear();
    return true;
  }

  header->SetRunNumber(0);
  header->SetEventNumber(Event_Number);

  //if(Event_Number%100==0){
  //  std::cout << "=======================================" << std::endl;
  //  std::cout << "Event : " << Event_Number << std::endl;
  //  std::cout << "CDC   : " << CDC_Event_Number << std::endl;
  //  std::cout << "MTDC  : " << MTDC_Event_Number << std::endl;
  //  std::cout << "Event : " << Event_Number << std::endl;
  //  std::cout << "CDC   : " << cdsheader->ev() << std::endl;
  //  std::cout << "MTDC  : " << mtdcheader->ev() << std::endl;
  //  std::cout << "=======================================" << std::endl;
  //}

  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );
  if(trackflag){
    cdstrackMan -> Calc(cdsMan,confMan);
  }

  /* pBeam */
  if(particle==0){
    Clear();
    return true;
  }
  if(particle->nBeam()==0){
    Clear();
    return true;
  }
  for(int it=0; it<particle->nBeam(); it++){
    particle2->AddBeam(*particle->beam(it));
  }
  rtFile2->cd();
  if(trackflag){
    if(!cdsAna->DoAnalysis(confMan, header, blMan, bltrackMan, cdsMan, cdstrackMan, particle2)){
      Clear();
      return true;
    }
    if(!fwdnAna->DoAnalysis(confMan, header, blMan, bltrackMan, cdsMan, cdstrackMan, particle2)){
      Clear();
      return true;
    }
    if(!fwdcAna->DoAnalysis(confMan, header, blMan, bltrackMan, cdsMan, cdstrackMan, particle2)){
      Clear();
      return true;
    }
  }

  rtFile->cd();
  evTree->Fill();

  Clear();
  return true;
}

void EventAnalysisMyReAll::Clear()
{
  particle->Clear();
  particle2->Clear();
  blAna->Clear();
  cdsAna->Clear();
  fwdnAna->Clear();
  fwdcAna->Clear();
  blMan->Clear();
  bltrackMan->Clear();
  cdsMan->Clear();
  header->Clear();
}

void EventAnalysisMyReAll::Finalize()
{
  std::cout << " Enter EventAnalysisMyReAll::Finalize " << std::endl;

  rtFile->cd();
  gFile->Write();
  gFile->Close();
  rtFile2->cd();
  gFile->Write();
  gFile->Close();

  delete particle;
  delete particle2;
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
  EventAnalysisMyReAll *event = new EventAnalysisMyReAll();
  return (EventTemp*)event;
}
