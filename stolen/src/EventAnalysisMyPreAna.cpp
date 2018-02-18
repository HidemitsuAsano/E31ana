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
#include "MyAnalysisBLDCCheck.h"
#include "MyAnalysisSDCCheck.h"
#include "MyAnalysisUSWKCheck.h"
#include "MyAnalysisCDHIHCheck.h"
#include "MyAnalysisCDCBPCCheck.h"
#include "MyAnalysisNCCVCCheck.h"
#include "MyAnalysisDEFTest.h"
#include "MyAnalysisHVCTest.h"
#include "MyAnalysisBLHodoCheck.h"
#include "MyAnalysisTrigCheck.h"
#include "MyAnalysisDetectorCheck.h"
#include "MyAnalysisFWDHodoCheck.h"
#include "MyAnalysisCVCPCCheck.h"
#include "MyAnalysisCDSCheck.h"
#include "MyAnalysisCDCTrackEffic.h"
#include "MyAnalysisTest.h"

#define Debug 0

//#define DEFTest
//#define HVCTest
//#define BLHodoCheck
#define TrigCheck
//#define BLDCCheck
//#define SDCCheck
//#define USWKCheck
//#define CDHIHCheck
//#define CDCBPCCheck
//#define NCCVCCheck
//#define FWDHodoCheck
//#define CVCPCCheck
//#define DetectorCheck
//#define CDSCheck
//#define CDCTrackEffic
//#define Test

class EventAnalysisMyPreAna: public EventTemp
{
  public:
    EventAnalysisMyPreAna();
    ~EventAnalysisMyPreAna();
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

#ifdef BLDCCheck
    MyAnalysisBLDCCheck* bldccheckAna;
#endif
#ifdef SDCCheck
    MyAnalysisSDCCheck* bldccheckAna;
#endif
#ifdef USWKCheck
    MyAnalysisUSWKCheck* uswkcheckAna;
#endif
#ifdef CDHIHCheck
    MyAnalysisCDHIHCheck* cdhihcheckAna;
#endif
#ifdef CDCBPCCheck
    MyAnalysisCDCBPCCheck* cdcbpccheckAna;
#endif
#ifdef NCCVCCheck
    MyAnalysisNCCVCCheck* nccvccheckAna;
#endif
#ifdef DEFTest
    MyAnalysisDEFTest* deftestAna;
#endif
#ifdef HVCTest
    MyAnalysisHVCTest* hvctestAna;
#endif
#ifdef BLHodoCheck
    MyAnalysisBLHodoCheck* hodocheckAna;
#endif
#ifdef TrigCheck
    MyAnalysisTrigCheck* trigcheckAna;
#endif
#ifdef DetectorCheck
    MyAnalysisDetectorCheck* detcheckAna;
#endif
#ifdef FWDHodoCheck
    MyAnalysisFWDHodoCheck* fwdhodocheckAna;
#endif
#ifdef CVCPCCheck
    MyAnalysisCVCPCCheck* cvcpccheckAna;
#endif
#ifdef CDSCheck
    MyAnalysisCDSCheck* cdscheckAna;
#endif
#ifdef CDCTrackEffic
    MyAnalysisCDCTrackEffic* cdctrackeffcAna;
#endif
#ifdef Test
    MyAnalysisTest* testAna;
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

  EventAnalysisMyPreAna::EventAnalysisMyPreAna()
: EventTemp()
{
}

EventAnalysisMyPreAna::~EventAnalysisMyPreAna()
{
}

const int MaxTreeSize = 1900000000000;
void EventAnalysisMyPreAna::Initialize( ConfMan *conf )
{
#if 1
  std::cout << " Enter EventAnalysisMyPreAna::Initialize " << std::endl;
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

#ifdef BLDCCheck
  bldccheckAna  = new MyAnalysisBLDCCheck(rtFile, confMan);
#endif
#ifdef SDCCheck
  bldccheckAna  = new MyAnalysisSDCCheck(rtFile, confMan);
#endif
#ifdef USWKCheck
  uswkcheckAna  = new MyAnalysisUSWKCheck(rtFile, confMan);
#endif
#ifdef CDHIHCheck
  cdhihcheckAna  = new MyAnalysisCDHIHCheck(rtFile, confMan);
#endif
#ifdef CDCBPCCheck
  cdcbpccheckAna  = new MyAnalysisCDCBPCCheck(rtFile, confMan);
#endif
#ifdef NCCVCCheck
  nccvccheckAna  = new MyAnalysisNCCVCCheck(rtFile, confMan);
#endif
#ifdef DEFTest
  deftestAna  = new MyAnalysisDEFTest(rtFile, confMan);
#endif
#ifdef HVCTest
  hvctestAna  = new MyAnalysisHVCTest(rtFile, confMan);
#endif
#ifdef BLHodoCheck
  hodocheckAna  = new MyAnalysisBLHodoCheck(rtFile, confMan);
#endif
#ifdef TrigCheck
  trigcheckAna  = new MyAnalysisTrigCheck(rtFile, confMan);
#endif
#ifdef DetectorCheck
  detcheckAna  = new MyAnalysisDetectorCheck(rtFile, confMan);
#endif
#ifdef FWDHodoCheck
  fwdhodocheckAna  = new MyAnalysisFWDHodoCheck(rtFile, confMan);
#endif
#ifdef CVCPCCheck
  cvcpccheckAna  = new MyAnalysisCVCPCCheck(rtFile, confMan);
#endif
#ifdef CDSCheck
  cdscheckAna  = new MyAnalysisCDSCheck(rtFile, confMan);
#endif
#ifdef CDCTrackEffic
  cdctrackefficAna  = new MyAnalysisCDCTrackEffic(rtFile, confMan);
#endif
#ifdef Test
  testAna  = new MyAnalysisTest(rtFile, confMan);
#endif

  AllGoodTrack = 0;
  nTrack = 0;
  CDC_Event_Number = 0;
  t0=time(0);
}

void EventAnalysisMyPreAna::USca( int nsca, unsigned int *sca )
{
  Block_Event_Number++;
  header->SetBlockEventNumber( Block_Event_Number );
  scaMan->SetBlockEventNumber( Block_Event_Number );
  for( int i=0; i<nsca; i++ ){
    scaMan->AddHit( sca[i], confMan->GetScalerMapManager()->GetName(i) );
  }

  scaMan->Clear();
}

bool EventAnalysisMyPreAna::UAna( TKOHitCollection *tko )
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
      std::cout << " Event# : " << Event_Number << "  BlockEvent# : " << Block_Event_Number <<  " GoodTrack# : " << AllGoodTrack  <<"  "<<CDC_Event_Number<<" / "<<cdcTree->GetEntries() << " Time (s):" << (t1-t0) << std::endl;
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

#ifdef BLDCCheck
  if(!bldccheckAna->DoAnalysis(confMan, header, blMan, bltrackMan)){
    //std::cout << " k0Ana::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef SDCCheck
  if(!bldccheckAna->DoAnalysis(confMan, header, blMan, bltrackMan)){
    //std::cout << " k0Ana::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef USWKCheck
  if(!uswkcheckAna->DoAnalysis(confMan, header, blMan, bltrackMan)){
    //std::cout << " k0Ana::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef CDHIHCheck
  if(trackflag){
    if(!cdhihcheckAna->DoAnalysis(confMan, header, blMan, bltrackMan, cdsMan, cdstrackMan)){
      //std::cout << " k0Ana::false " << std::endl;
      //Clear();
      //return true;
    }
  }
#endif
#ifdef CDCBPCCheck
  if(trackflag){
    if(!cdcbpccheckAna->DoAnalysis(confMan, header, blMan, bltrackMan, cdsMan, cdstrackMan)){
      //std::cout << " k0Ana::false " << std::endl;
      //Clear();
      //return true;
    }
  }
#endif
#ifdef NCCVCCheck
  if(!nccvccheckAna->DoAnalysis(confMan, header, blMan, bltrackMan, cdsMan, cdstrackMan)){
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
#ifdef BLHodoCheck
  if(!hodocheckAna->DoAnalysis(confMan, header, blMan, bltrackMan)){
    //std::cout << " k0Ana::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef TrigCheck
  if(!trigcheckAna->DoAnalysis(confMan, header, blMan, cdsMan)){
    //std::cout << " k0Ana::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef DetectorCheck
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
#ifdef CVCPCCheck
  if(!cvcpccheckAna->DoAnalysis(confMan, header, blMan, bltrackMan)){
    //std::cout << " k0Ana::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef CDSCheck
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
#ifdef CDCTrackEffic
  if(trackflag){
    if(!cdctrackefficAna->DoAnalysis(confMan, header, blMan, cdsMan, cdstrackMan)){
      //std::cout << " k0Ana::false " << std::endl;
      //Clear();
      //return true;
    }
  }
  else{
    if(!cdctrackefficAna->DoAnalysis(confMan, header, blMan, cdsMan, 0)){
      //std::cout << " k0Ana::false " << std::endl;
      //Clear();
      //return true;
    }
  }
#endif
#ifdef Test
  if(!testAna->DoAnalysis(confMan, header, blMan, bltrackMan)){
    //std::cout << " k0Ana::false " << std::endl;
    //Clear();
    //return true;
  }
#endif

  Clear();
  return true;
}

void EventAnalysisMyPreAna::Clear()
{
  trackflag = false;
#ifdef BLDCCheck
  bldccheckAna->Clear();
#endif
#ifdef SDCCheck
  bldccheckAna->Clear();
#endif
#ifdef USWKCheck
  uswkcheckAna->Clear();
#endif
#ifdef CDHIHCheck
  cdhihcheckAna->Clear();
#endif
#ifdef CDCBPCCheck
  cdcbpccheckAna->Clear();
#endif
#ifdef NCCVCCheck
  nccvccheckAna->Clear();
#endif
#ifdef DEFTest
  deftestAna->Clear();
#endif
#ifdef HVCTest
  hvctestAna->Clear();
#endif
#ifdef BLHodoCheck
  hodocheckAna->Clear();
#endif
#ifdef TrigCheck
  trigcheckAna->Clear();
#endif
#ifdef FWDHodoCheck
  fwdhodocheckAna->Clear();
#endif
#ifdef CVCPCCheck
  cvcpccheckAna->Clear();
#endif
#ifdef CDSCheck
  cdscheckAna->Clear();
#endif
#ifdef CDCTrackEffic
  cdctrackefficAna->Clear();
#endif
#ifdef Test
  testAna->Clear();
#endif
  blMan->Clear();
  bltrackMan->Clear();
  cdsMan->Clear();
  header->Clear();
}

void EventAnalysisMyPreAna::Finalize()
{
  std::cout << " Enter EventAnalysisMyPreAna::Finalize " << std::endl;

  rtFile->cd();
  gFile->Write();
  gFile->Close();

#ifdef BLDCCheck
  delete bldccheckAna;
#endif
#ifdef SDCCheck
  delete bldccheckAna;
#endif
#ifdef USWKCheck
  delete uswkcheckAna;
#endif
#ifdef CDHIHCheck
  delete cdhihcheckAna;
#endif
#ifdef CDCBPCCheck
  delete cdcbpccheckAna;
#endif
#ifdef NCCVCCheck
  delete nccvccheckAna;
#endif
#ifdef DEFTest
  delete deftestAna;
#endif
#ifdef HVCTest
  delete hvctestAna;
#endif
#ifdef BLHodoCheck
  delete hodocheckAna;
#endif
#ifdef TrigCheck
  delete trigcheckAna;
#endif
#ifdef DetectorCheck
  delete detcheckAna;
#endif
#ifdef FWDHodoCheck
  delete fwdhodocheckAna;
#endif
#ifdef CVCPCCheck
  delete cvcpccheckAna;
#endif
#ifdef CDSCheck
  delete cdscheckAna;
#endif
#ifdef CDCTrackEffic
  delete cdctrackefficAna;
#endif
#ifdef Test
  delete testAna;
#endif
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisMyPreAna *event = new EventAnalysisMyPreAna();
  return (EventTemp*)event;
}
