#include "File.h"
#include "ConfMan.h"
#include "TKO.h"
#include "CDSHitMan.h"
#include "BeamLineHitMan.h"
#include "BeamLineTrackMan.h"
#include "EventHeader.h"
#include "ScalerMan.h"

#include "EventAlloc.h"
#include "EventTemp.h"
#include "Display.h" 
#include "ELossTools.h"
#include "TrackTools.h"
#include "MyAnalysisBLDCCheckSource.h"
#include "MyAnalysisDEFTest.h"
#include "MyAnalysisHVCTest.h"
#include "MyAnalysisHodoCheck.h"
#include "MyAnalysisDetectorCheckSource.h"
#include "MyAnalysisFWDHodoCheck.h"
#include "MyAnalysisCDSCheckCosmic.h"

#define Debug 0

//#define BLDCCheckSource
//#define DEFTest
//#define HVCTest
//#define HodoCheck
//#define FWDHodoCheck
//#define DetectorCheckSource
#define CDSCheckCosmic

class EventAnalysisMyPreAnaCosmic: public EventTemp
{
  public:
    EventAnalysisMyPreAnaCosmic();
    ~EventAnalysisMyPreAnaCosmic();
  private:
    TFile *rtFile;
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
    ScalerMan *scaMan;

#ifdef BLDCCheckSource
    MyAnalysisBLDCCheckSource* bldccheckAna;
#endif
#ifdef DEFTest
    MyAnalysisDEFTest* deftestAna;
#endif
#ifdef HVCTest
    MyAnalysisHVCTest* hvctestAna;
#endif
#ifdef HodoCheck
    MyAnalysisHodoCheck* hodocheckAna;
#endif
#ifdef DetectorCheckSource
    MyAnalysisDetectorCheckSource* detcheckAna;
#endif
#ifdef FWDHodoCheck
    MyAnalysisFWDHodoCheck* fwdhodocheckAna;
#endif
#ifdef CDSCheckCosmic
    MyAnalysisCDSCheckCosmic* cdscheckAna;
#endif
    int t0, t1;
    int AllGoodTrack;
    int nTrack;
    int CDC_Event_Number;
    bool trackflag;

  public:
    void Initialize( ConfMan *conf );
    void USca( int nsca, unsigned int *sca );
    bool UAna( TKOHitCollection *tko );
    void Finalize();
    void Clear();

};

  EventAnalysisMyPreAnaCosmic::EventAnalysisMyPreAnaCosmic()
: EventTemp()
{
}

EventAnalysisMyPreAnaCosmic::~EventAnalysisMyPreAnaCosmic()
{
}

const int MaxTreeSize = 1900000000000;
void EventAnalysisMyPreAnaCosmic::Initialize( ConfMan *conf )
{
#if 1
  std::cout << " Enter EventAnalysisMyPreAnaCosmic::Initialize " << std::endl;
#endif
  confMan = conf;

  rtFile = new TFile( confMan->GetOutFileName().c_str(), "recreate" );
  header = new EventHeader();
  if( header==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  cdsheader = new EventHeader();
  if( cdsheader==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  cdsMan = new CDSHitMan();
  if( cdsMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  blMan = new BeamLineHitMan();
  if( blMan ==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  bltrackMan = new BeamLineTrackMan();
  if( bltrackMan ==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  cdstrackMan = new CDSTrackingMan();
  if( cdstrackMan ==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  scaMan = new ScalerMan();
  if( scaMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }

  /* CDS Tracking file */
  TString cdcfname=confMan->GetCDSTrackFileName();
  if(cdcfname!="None!!"){
    cdcFile = new TFile(cdcfname);
    if(!cdcFile->IsOpen()){
      std::cout<<" failed to open " <<cdcfname<< "  !!!"<<std::endl;
      exit(false);
    }
    std::cout<<" CDC Track File Successfully opend !!! " <<cdcfname<<std::endl;
    cdcTree=(TTree*)cdcFile->Get("EventTree");
    cdcTree->SetBranchAddress( "CDSTrackingMan", &cdstrackMan );
    cdcTree->SetBranchAddress( "EventHeader" ,&cdsheader);
  }

#ifdef BLDCCheckSource
  bldccheckAna  = new MyAnalysisBLDCCheckSource(rtFile, confMan);
#endif
#ifdef DEFTest
  deftestAna  = new MyAnalysisDEFTest(rtFile, confMan);
#endif
#ifdef HVCTest
  hvctestAna  = new MyAnalysisHVCTest(rtFile, confMan);
#endif
#ifdef HodoCheck
  hodocheckAna  = new MyAnalysisHodoCheck(rtFile, confMan);
#endif
#ifdef DetectorCheckSource
  detcheckAna  = new MyAnalysisDetectorCheckSource(rtFile, confMan);
#endif
#ifdef FWDHodoCheck
  fwdhodocheckAna  = new MyAnalysisFWDHodoCheck(rtFile, confMan);
#endif
#ifdef CDSCheckCosmic
  cdscheckAna  = new MyAnalysisCDSCheckCosmic(rtFile, confMan);
#endif

  AllGoodTrack = 0;
  nTrack = 0;
  CDC_Event_Number = 0;
  t0=time(0);
}

void EventAnalysisMyPreAnaCosmic::USca( int nsca, unsigned int *sca )
{
  Block_Event_Number++;
  header->SetBlockEventNumber( Block_Event_Number );
  scaMan->SetBlockEventNumber( Block_Event_Number );
  for( int i=0; i<nsca; i++ ){
    scaMan->AddHit( sca[i], confMan->GetScalerMapManager()->GetName(i) );
  }

  scaMan->Clear();
}

bool EventAnalysisMyPreAnaCosmic::UAna( TKOHitCollection *tko )
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

  if( Event_Number%10000==0 )
  {
    t1=time(0);
    std::cout << " Event# : " << Event_Number << "  BlockEvent# : " << Block_Event_Number <<  " Time (s):" << (t1-t0);
    std::cout << " All Tracks : " << nTrack << ", All GoodTrack : "
      << AllGoodTrack << std::endl;
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
      //Clear();
      //return true;
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

  rtFile->cd();

#ifdef BLDCCheckSource
  if(!bldccheckAna->DoAnalysis(confMan, header, blMan, bltrackMan)){
    //std::cout << " k0Ana::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef DEFTest
  if(!deftestAna->DoAnalysis(confMan, header, blMan, bltrackMan)){
    //std::cout << " k0Ana::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef HVCTest
  if(!hvctestAna->DoAnalysis(confMan, header, blMan, bltrackMan)){
    //std::cout << " k0Ana::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef HodoCheck
  if(!hodocheckAna->DoAnalysis(confMan, header, blMan, bltrackMan)){
    //std::cout << " k0Ana::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef DetectorCheckSource
  if(!detcheckAna->DoAnalysis(confMan, header, blMan, bltrackMan)){
    //std::cout << " k0Ana::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef FWDHodoCheck
  if(!fwdhodocheckAna->DoAnalysis(confMan, header, blMan, bltrackMan)){
    //std::cout << " k0Ana::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef CDSCheckCosmic
  if(trackflag){
    if(!cdscheckAna->DoAnalysis(confMan, header, blMan, cdsMan, cdstrackMan)){
      //std::cout << " k0Ana::false " << std::endl;
      //Clear();
      //return true;
    }
  }
  else{
    if(!cdscheckAna->DoAnalysis(confMan, header, blMan, cdsMan, 0)){
      //std::cout << " k0Ana::false " << std::endl;
      //Clear();
      //return true;
    }
  }
#endif

  Clear();
  return true;
}

void EventAnalysisMyPreAnaCosmic::Clear()
{
  trackflag = false;
#ifdef BLDCCheckSource
  bldccheckAna->Clear();
#endif
#ifdef DEFTest
  deftestAna->Clear();
#endif
#ifdef HVCTest
  hvctestAna->Clear();
#endif
#ifdef HodoCheck
  hodocheckAna->Clear();
#endif
#ifdef FWDHodoCheck
  fwdhodocheckAna->Clear();
#endif
#ifdef CDSCheckCosmic
  cdscheckAna->Clear();
#endif
  blMan->Clear();
  bltrackMan->Clear();
  cdsMan->Clear();
  header->Clear();
}

void EventAnalysisMyPreAnaCosmic::Finalize()
{
  std::cout << " Enter EventAnalysisMyPreAnaCosmic::Finalize " << std::endl;

  rtFile->cd();
  gFile->Write();
  gFile->Close();

#ifdef BLDCCheckSource
  delete bldccheckAna;
#endif
#ifdef DEFTest
  delete deftestAna;
#endif
#ifdef HVCTest
  delete hvctestAna;
#endif
#ifdef HodoCheck
  delete hodocheckAna;
#endif
#ifdef DetectorCheckSource
  delete detcheckAna;
#endif
#ifdef FWDHodoCheck
  delete fwdhodocheckAna;
#endif
#ifdef CDSCheckCosmic
  delete cdscheckAna;
#endif
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisMyPreAnaCosmic *event = new EventAnalysisMyPreAnaCosmic();
  return (EventTemp*)event;
}
