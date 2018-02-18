#include "File.h"
#include "ConfMan.h"
#include "TKO.h"
#include "CDSHitMan.h"
#include "BeamLineHitMan.h"
#include "BeamLineTrackMan.h"
#include "EventHeader.h"
#include "ScalerMan.h"
#include "Particle.h"
#include "CDHCluster.cpp"

#include "EventAlloc.h"
#include "EventTemp.h"
#include "Display.h" 
#include "ELossTools.h"
#include "TrackTools.h"
#include "MyAnalysisBasic.h"
#include "MyAnalysisK0.h"
#include "MyAnalysisK0pn.h"
#include "MyAnalysisKpn.h"
#include "MyAnalysisK0n.h"
#include "MyAnalysispKn.h"
#include "MyAnalysisdKn.h"
#include "MyAnalysisdKnpipi.h"
#include "MyAnalysisdKnppim.h"
#include "MyAnalysis3HeKn.h"
#include "MyAnalysisHeKFWD.h"
#include "MyAnalysisHeKpipip.h"
#include "MyAnalysisHeKpipipp.h"
#include "MyAnalysisHeKFWDNeutral.h"
#include "MyAnalysisHeK2CDS.h"
#include "MyAnalysisCDCSlew.h"
#include "MyAnalysisdKFWD.h"
#include "MyAnalysisdKCDS.h"
#include "MyAnalysisCDCTrackingwSlew.h"
#include "MyAnalysisdKCDH3.h"
#include "MyAnalysisdKpipiFWD.h"
#include "MyAnalysisCheck.h"
#include "MyAnalysisHeKpipp.h"
#include "MyAnalysisHeKpippOneDummy.h"
#include "MyAnalysisHeKpippOneDummyLpair.h"
#include "MyAnalysisHeKpippMixEvent.h"
#include "MyAnalysisHeKLpnFWD.h"
#include "MyAnalysisCheckChi2.h"
#include "MyAnalysisHeKpFWD.h"
#include "MyAnalysisdKppipi.h"
#include "MyAnalysisHeKpn.h"
#include "MyAnalysisHeKppin.h"
#include "MyAnalysisOverVeto.h"
#include "MyAnalysisCDSCalib.h"
#include "MyAnalysisNCCalib.h"
#include "MyAnalysis2Track.h"
#include "MyAnalysisLambda.h"

#define Debug 0

//#define Basic
//#define pKn
//#define dKn
//#define dKnpipi
//#define dKnppim
//#define HeKn
//#define HeKFWD
//#define HeKpipip
//#define HeKpipipp
//#define HeKFWDNeutral
//#define HeK2CDS
//#define CDCSlew
//#define CDCTrackingwSlew
//#define dKFWD
//#define dKCDS
//#define dKCDH3
//#define dKpipiFWD
#define Check
//#define HeKpipp
//#define HeKpippOneDummy
//#define HeKpippOneDummyLpair
//#define HeKpippMixEvent
//#define HeKLpnFWD
//#define CheckChi2
//#define HeKpFWD
//#define dKppipi
//#define HeKpn
//#define HeKppin
//#define OverVeto
//#define CDSCalib
//#define NCCalib
//#define TwoTrack
//#define K0
//#define K0pn
//#define Kpn
//#define K0n
//#define Lambda

class EventAnalysisMyAna: public EventTemp
{
  public:
    EventAnalysisMyAna();
    ~EventAnalysisMyAna();
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

#ifdef Basic
    MyAnalysisBasic* basicAna;
#endif
#ifdef K0
    MyAnalysisK0* k0Ana;
#endif
#ifdef K0pn
    MyAnalysisK0pn* k0pnAna;
#endif
#ifdef Kpn
    MyAnalysisKpn* kpnAna;
#endif
#ifdef K0n
    MyAnalysisK0n* k0nAna;
#endif
#ifdef pKn
    MyAnalysispKn* pKnAna;
#endif
#ifdef dKn
    MyAnalysisdKn* dKnAna;
#endif
#ifdef dKnpipi
    MyAnalysisdKnpipi* dKnpipiAna;
#endif
#ifdef dKnppim
    MyAnalysisdKnppim* dKnppimAna;
#endif
#ifdef HeKn
    MyAnalysis3HeKn* HeKnAna;
#endif
#ifdef HeKFWD
    MyAnalysisHeKFWD* HeKFWDAna;
#endif
#ifdef HeKpipip
    MyAnalysisHeKpipip* HeKpipipAna;
#endif
#ifdef HeKpipipp
    MyAnalysisHeKpipipp* HeKpipippAna;
#endif
#ifdef HeKFWDNeutral
    MyAnalysisHeKFWDNeutral* HeKFWDNeutralAna;
#endif
#ifdef HeK2CDS
    MyAnalysisHeK2CDS* HeK2CDSAna;
#endif
#ifdef CDCSlew
    MyAnalysisCDCSlew* CDCSlewAna;
#endif
#ifdef CDCTrackingwSlew
    MyAnalysisCDCTrackingwSlew* CDCTrackingwSlewAna;
#endif
#ifdef dKFWD
    MyAnalysisdKFWD* dKFWDAna;
#endif
#ifdef dKCDS
    MyAnalysisdKCDS* dKCDSAna;
#endif
#ifdef dKCDH3
    MyAnalysisdKCDH3* dKCDH3Ana;
#endif
#ifdef dKpipiFWD
    MyAnalysisdKpipiFWD* dKpipiFWDAna;
#endif
#ifdef Check
    MyAnalysisCheck* CheckAna;
#endif
#ifdef HeKpipp
    MyAnalysisHeKpipp* HeKpippAna;
#endif
#ifdef HeKpippOneDummy
    MyAnalysisHeKpippOneDummy* HeKpippOneDummyAna;
#endif
#ifdef HeKpippOneDummyLpair
    MyAnalysisHeKpippOneDummyLpair* HeKpippOneDummyLpairAna;
#endif
#ifdef HeKpippMixEvent
    MyAnalysisHeKpippMixEvent* HeKpippMixEventAna;
#endif
#ifdef HeKLpnFWD
    MyAnalysisHeKLpnFWD* HeKLpnFWDAna;
#endif
#ifdef CheckChi2
    MyAnalysisCheckChi2* CheckChi2Ana;
#endif
#ifdef HeKpFWD
    MyAnalysisHeKpFWD* HeKpFWDAna;
#endif
#ifdef dKppipi
    MyAnalysisdKppipi* dKppipiAna;
#endif
#ifdef HeKpn
    MyAnalysisHeKpn* HeKpnAna;
#endif
#ifdef HeKppin
    MyAnalysisHeKppin* HeKppinAna;
#endif
#ifdef OverVeto
    MyAnalysisOverVeto* OverVetoAna;
#endif
#ifdef CDSCalib
    MyAnalysisCDSCalib* CDSCalibAna;
#endif
#ifdef NCCalib
    MyAnalysisNCCalib* NCCalibAna;
#endif
#ifdef TwoTrack
    MyAnalysis2Track* TwoTrackAna;
#endif
#ifdef Lambda
    MyAnalysisLambda* LambdaAna;
#endif

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

  EventAnalysisMyAna::EventAnalysisMyAna()
: EventTemp()
{
}

EventAnalysisMyAna::~EventAnalysisMyAna()
{
}

const int MaxTreeSize = 1900000000000;
void EventAnalysisMyAna::Initialize( ConfMan *conf )
{
#if 1
  std::cout << " Enter EventAnalysisMyAna::Initialize " << std::endl;
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

#ifdef Basic
  basicAna  = new MyAnalysisBasic(rtFile, confMan);
#endif
#ifdef K0
  k0Ana  = new MyAnalysisK0(rtFile, confMan);
#endif
#ifdef K0pn
  k0pnAna  = new MyAnalysisK0pn(rtFile, confMan);
#endif
#ifdef Kpn
  kpnAna  = new MyAnalysisKpn(rtFile, confMan);
#endif
#ifdef K0n
  k0nAna  = new MyAnalysisK0n(rtFile, confMan);
#endif
#ifdef pKn
  pKnAna  = new MyAnalysispKn(rtFile, confMan);
#endif
#ifdef dKn
  dKnAna  = new MyAnalysisdKn(rtFile, confMan);
#endif
#ifdef dKnpipi
  dKnpipiAna  = new MyAnalysisdKnpipi(rtFile, confMan);
#endif
#ifdef dKnppim
  dKnppimAna  = new MyAnalysisdKnppim(rtFile, confMan);
#endif
#ifdef HeKn
  HeKnAna  = new MyAnalysis3HeKn(rtFile, confMan);
#endif
#ifdef HeKFWD
  HeKFWDAna  = new MyAnalysisHeKFWD(rtFile, confMan);
#endif
#ifdef HeKpipip
  HeKpipipAna  = new MyAnalysisHeKpipip(rtFile, confMan);
#endif
#ifdef HeKpipipp
  HeKpipippAna  = new MyAnalysisHeKpipipp(rtFile, confMan);
#endif
#ifdef HeKFWDNeutral
  HeKFWDNeutralAna  = new MyAnalysisHeKFWDNeutral(rtFile, confMan);
#endif
#ifdef HeK2CDS
  HeK2CDSAna  = new MyAnalysisHeK2CDS(rtFile, confMan);
#endif
#ifdef CDCSlew
  CDCSlewAna  = new MyAnalysisCDCSlew(rtFile, confMan);
#endif
#ifdef CDCTrackingwSlew
  CDCTrackingwSlewAna  = new MyAnalysisCDCTrackingwSlew(rtFile, confMan);
#endif
#ifdef dKFWD
  dKFWDAna  = new MyAnalysisdKFWD(rtFile, confMan);
#endif
#ifdef dKCDS
  dKCDSAna  = new MyAnalysisdKCDS(rtFile, confMan);
#endif
#ifdef dKCDH3
  dKCDH3Ana  = new MyAnalysisdKCDH3(rtFile, confMan);
#endif
#ifdef dKpipiFWD
  dKpipiFWDAna  = new MyAnalysisdKpipiFWD(rtFile, confMan);
#endif
#ifdef Check
  CheckAna  = new MyAnalysisCheck(rtFile, confMan);
#endif
#ifdef HeKpipp
  HeKpippAna  = new MyAnalysisHeKpipp(rtFile, confMan);
#endif
#ifdef HeKpippOneDummy
  HeKpippOneDummyAna  = new MyAnalysisHeKpippOneDummy(rtFile, confMan);
#endif
#ifdef HeKpippOneDummyLpair
  HeKpippOneDummyLpairAna  = new MyAnalysisHeKpippOneDummyLpair(rtFile, confMan);
#endif
#ifdef HeKpippMixEvent
  HeKpippMixEventAna  = new MyAnalysisHeKpippMixEvent(rtFile, confMan);
#endif
#ifdef HeKLpnFWD
  HeKLpnFWDAna  = new MyAnalysisHeKLpnFWD(rtFile, confMan);
#endif
#ifdef CheckChi2
  CheckChi2Ana  = new MyAnalysisCheckChi2(rtFile, confMan);
#endif
#ifdef HeKpFWD
  HeKpFWDAna  = new MyAnalysisHeKpFWD(rtFile, confMan);
#endif
#ifdef dKppipi
  dKppipiAna  = new MyAnalysisdKppipi(rtFile, confMan);
#endif
#ifdef HeKpn
  HeKpnAna  = new MyAnalysisHeKpn(rtFile, confMan);
#endif
#ifdef HeKppin
  HeKppinAna  = new MyAnalysisHeKppin(rtFile, confMan);
#endif
#ifdef OverVeto
  OverVetoAna  = new MyAnalysisOverVeto(rtFile, confMan);
#endif
#ifdef CDSCalib
  CDSCalibAna  = new MyAnalysisCDSCalib(rtFile, confMan);
#endif
#ifdef NCCalib
  NCCalibAna  = new MyAnalysisNCCalib(rtFile, confMan);
#endif
#ifdef TwoTrack
  TwoTrackAna  = new MyAnalysis2Track(rtFile, confMan);
#endif
#ifdef Lambda
  LambdaAna  = new MyAnalysisLambda(rtFile, confMan);
#endif

  AllGoodTrack = 0;
  nAllTrack = 0;
  CDC_Event_Number = 0;
  MTDC_Event_Number = 0;
  t0=time(0);
}

void EventAnalysisMyAna::USca( int nsca, unsigned int *sca )
{
  Block_Event_Number++;
  header->SetBlockEventNumber( Block_Event_Number );
  scaMan->SetBlockEventNumber( Block_Event_Number );
  for( int i=0; i<nsca; i++ ){
    scaMan->AddHit( sca[i], confMan->GetScalerMapManager()->GetName(i) );
  }

  scaMan->Clear();
}

bool EventAnalysisMyAna::UAna( TKOHitCollection *tko )
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
    std::cout << " Event# : " << Event_Number << "  BlockEvent# : " << Block_Event_Number << " AllTrack# : " << nAllTrack <<" GoodTrack# : " << AllGoodTrack  <<"  "<<CDC_Event_Number<<" / "<<cdcTree->GetEntries() << " Time (s):" << (t1-t0) << std::endl;
  }

  if(CDC_Event_Number>=cdcTree->GetEntries()){
    Clear();
    return true;  
  }
  while( cdsheader->ev()<Event_Number ){
    CDC_Event_Number++;
    if(CDC_Event_Number>cdcTree->GetEntries()){
      Clear();  
      return true;
    }
    cdcTree -> GetEntry(CDC_Event_Number);
  }
  if(cdsheader->ev()>Event_Number){
    Clear();
    return true;
  }
  int nGoodTrack = cdstrackMan -> nGoodTrack();
  int nAllTrack = cdstrackMan -> nTrack();
  AllGoodTrack += nGoodTrack;
  nAllTrack += nAllTrack;

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

#ifdef NCCalib
  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );
  cdstrackMan -> Calc( cdsMan, confMan );
#endif
#ifdef CDSCalib
  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );
  cdstrackMan -> Calc( cdsMan, confMan );
#endif
#ifdef HeKFWD
  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );
  cdstrackMan -> Calc( cdsMan, confMan );
#endif
#ifdef HeKpipip
  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );
  cdstrackMan -> Calc( cdsMan, confMan );
#endif
#ifdef HeKpipipp
  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );
  cdstrackMan -> Calc( cdsMan, confMan );
#endif
#ifdef HeKFWDNeutral
  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );
  cdstrackMan -> Calc( cdsMan, confMan );
#endif
#ifdef HeK2CDS
  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );
  cdstrackMan -> Calc( cdsMan, confMan );
#endif
#ifdef CDCSlew
  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );
  cdstrackMan -> Calc( cdsMan, confMan );
#endif
#ifdef CDCTrackingwSlew
  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );
  cdstrackMan -> Calc( cdsMan, confMan );
#endif
#ifdef dKFWD
  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );
  cdstrackMan -> Calc( cdsMan, confMan );
#endif
#ifdef dKCDS
  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );
  cdstrackMan -> Calc( cdsMan, confMan );
#endif
#ifdef dKpipiFWD
  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );
  cdstrackMan -> Calc( cdsMan, confMan );
#endif
#ifdef Check
  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );
  cdstrackMan -> Calc( cdsMan, confMan );
#endif
#ifdef HeKpipp
  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );
  cdstrackMan -> Calc( cdsMan, confMan );
#endif
#ifdef HeKpippOneDummy
  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );
  cdstrackMan -> Calc( cdsMan, confMan );
#endif
#ifdef HeKpippOneDummyLpair
  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );
  cdstrackMan -> Calc( cdsMan, confMan );
#endif
#ifdef HeKpippMixEvent
  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );
  cdstrackMan -> Calc( cdsMan, confMan );
#endif
#ifdef HeKLpnFWD
  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );
  cdstrackMan -> Calc( cdsMan, confMan );
#endif
#ifdef CheckChi2
  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );
  cdstrackMan -> Calc( cdsMan, confMan );
#endif
#ifdef HeKpFWD
  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );
  cdstrackMan -> Calc( cdsMan, confMan );
#endif
#ifdef dKppipi
  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );
  cdstrackMan -> Calc( cdsMan, confMan );
#endif
#ifdef HeKpn
  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );
  cdstrackMan -> Calc( cdsMan, confMan );
#endif
#ifdef HeKppin
  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );
  cdstrackMan -> Calc( cdsMan, confMan );
#endif
#ifdef TwoTrack
  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );
  cdstrackMan -> Calc( cdsMan, confMan );
#endif
#ifdef Lambda
  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );
  cdstrackMan -> Calc( cdsMan, confMan );
#endif
#ifdef K0
  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );
  cdstrackMan -> Calc( cdsMan, confMan );
#endif
#ifdef K0pn
  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );
  cdstrackMan -> Calc( cdsMan, confMan );
#endif
#ifdef Kpn
  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );
  cdstrackMan -> Calc( cdsMan, confMan );
#endif
#ifdef K0n
  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );
  cdstrackMan -> Calc( cdsMan, confMan );
#endif

  rtFile->cd();

  //double beamtof = particle->beam(0)->bhdt0tof();
  //if(beamtof<27.7854||29.6031<beamtof){
  //  Clear();
  //  return true;
  //}

  //if(!blAna->DoAnalysis(confMan, header, blMan, bltrackMan, particle)){
  //  Clear();
  //  return true;
  //}
#ifdef Basic
  if(!basicAna->DoAnalysis(confMan, mtdcheader, blMan, bltrackMan, cdsMan, cdstrackMan, particle)){
    //std::cout << " k0Ana::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef K0
  if(!k0Ana->DoAnalysis(confMan, mtdcheader, blMan, bltrackMan, cdsMan, cdstrackMan, particle)){
    //std::cout << " k0Ana::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef K0pn
  if(!k0pnAna->DoAnalysis(confMan, mtdcheader, blMan, bltrackMan, cdsMan, cdstrackMan, particle)){
    //std::cout << " k0Ana::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef Kpn
  if(!kpnAna->DoAnalysis(confMan, mtdcheader, blMan, bltrackMan, cdsMan, cdstrackMan, particle)){
    //std::cout << " k0Ana::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef K0n
  if(!k0nAna->DoAnalysis(confMan, mtdcheader, blMan, bltrackMan, cdsMan, cdstrackMan, particle)){
    //std::cout << " k0Ana::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef pKn
  if(!pKnAna->DoAnalysis(confMan, mtdcheader, blMan, bltrackMan, cdsMan, cdstrackMan, particle)){
    //std::cout << " pKnAna::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef dKn
  if(!dKnAna->DoAnalysis(confMan, mtdcheader, blMan, bltrackMan, cdsMan, cdstrackMan, particle)){
    //std::cout << " dKnAna::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef dKnpipi
  if(!dKnpipiAna->DoAnalysis(confMan, mtdcheader, blMan, bltrackMan, cdsMan, cdstrackMan, particle)){
    //std::cout << " dKnpipiAna::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef dKnppim
  if(!dKnppimAna->DoAnalysis(confMan, mtdcheader, blMan, bltrackMan, cdsMan, cdstrackMan, particle)){
    //std::cout << " dKnppimAna::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef HeKn
  if(!HeKnAna->DoAnalysis(confMan, mtdcheader, blMan, bltrackMan, cdsMan, cdstrackMan, particle)){
    //std::cout << " HeKnAna::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef HeKFWD
  if(!HeKFWDAna->DoAnalysis(confMan, mtdcheader, blMan, bltrackMan, cdsMan, cdstrackMan, particle)){
    //std::cout << " HeKFWDAna::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef HeKpipip
  if(!HeKpipipAna->DoAnalysis(confMan, mtdcheader, blMan, bltrackMan, cdsMan, cdstrackMan, particle)){
    //std::cout << " HeKpipipAna::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef HeKpipipp
  if(!HeKpipippAna->DoAnalysis(confMan, mtdcheader, blMan, bltrackMan, cdsMan, cdstrackMan, particle)){
    //std::cout << " HeKpipippAna::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef HeKFWDNeutral
  if(!HeKFWDNeutralAna->DoAnalysis(confMan, mtdcheader, blMan, bltrackMan, cdsMan, cdstrackMan, particle)){
    //std::cout << " HeKFWDNeutralAna::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef HeK2CDS
  if(!HeK2CDSAna->DoAnalysis(confMan, mtdcheader, blMan, bltrackMan, cdsMan, cdstrackMan, particle)){
    //std::cout << " HeK2CDSAna::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef CDCSlew
  if(!CDCSlewAna->DoAnalysis(confMan, mtdcheader, blMan, bltrackMan, cdsMan, cdstrackMan, particle)){
    //std::cout << " CDCSlewAna::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef CDCTrackingwSlew
  if(!CDCTrackingwSlewAna->DoAnalysis(confMan, mtdcheader, blMan, bltrackMan, cdsMan, cdstrackMan, particle)){
    //std::cout << " CDCTrackingwSlewAna::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef dKFWD
  if(!dKFWDAna->DoAnalysis(confMan, mtdcheader, blMan, bltrackMan, cdsMan, cdstrackMan, particle)){
    //std::cout << " dKFWDAna::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef dKCDS
  if(!dKCDSAna->DoAnalysis(confMan, mtdcheader, blMan, bltrackMan, cdsMan, cdstrackMan, particle)){
    //std::cout << " dKCDSAna::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef dKCDH3
  if(!dKCDH3Ana->DoAnalysis(confMan, mtdcheader, blMan, bltrackMan, cdsMan, cdstrackMan, particle)){
    //std::cout << " dKCDH3Ana::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef dKpipiFWD
  if(!dKpipiFWDAna->DoAnalysis(confMan, mtdcheader, blMan, bltrackMan, cdsMan, cdstrackMan, particle)){
    //std::cout << " dKpipiFWDAna::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef Check
  if(!CheckAna->DoAnalysis(confMan, mtdcheader, blMan, bltrackMan, cdsMan, cdstrackMan, particle)){
    //std::cout << " CheckAna::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef HeKpipp
  if(!HeKpippAna->DoAnalysis(confMan, mtdcheader, blMan, bltrackMan, cdsMan, cdstrackMan, particle)){
    //std::cout << " HeKpippAna::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef HeKpippOneDummy
  if(!HeKpippOneDummyAna->DoAnalysis(confMan, mtdcheader, blMan, bltrackMan, cdsMan, cdstrackMan, particle)){
    //std::cout << " HeKpippOneDummyAna::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef HeKpippOneDummyLpair
  if(!HeKpippOneDummyLpairAna->DoAnalysis(confMan, mtdcheader, blMan, bltrackMan, cdsMan, cdstrackMan, particle)){
    //std::cout << " HeKpippOneDummyLpairAna::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef HeKpippMixEvent
  if(!HeKpippMixEventAna->DoAnalysis(confMan, mtdcheader, blMan, bltrackMan, cdsMan, cdstrackMan, particle)){
    //std::cout << " HeKpippMixEventAna::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef HeKLpnFWD
  if(!HeKLpnFWDAna->DoAnalysis(confMan, mtdcheader, blMan, bltrackMan, cdsMan, cdstrackMan, particle)){
    //std::cout << " HeKLpnFWDAna::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef CheckChi2
  if(!CheckChi2Ana->DoAnalysis(confMan, mtdcheader, blMan, bltrackMan, cdsMan, cdstrackMan, particle)){
    //std::cout << " CheckChi2Ana::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef HeKpFWD
  if(!HeKpFWDAna->DoAnalysis(confMan, mtdcheader, blMan, bltrackMan, cdsMan, cdstrackMan, particle)){
    //std::cout << " HeKpFWDAna::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef dKppipi
  if(!dKppipiAna->DoAnalysis(confMan, mtdcheader, blMan, bltrackMan, cdsMan, cdstrackMan, particle)){
    //std::cout << " dKppipiAna::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef HeKpn
  if(!HeKpnAna->DoAnalysis(confMan, mtdcheader, blMan, bltrackMan, cdsMan, cdstrackMan, particle)){
    //std::cout << " HeKpnAna::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef HeKppin
  if(!HeKppinAna->DoAnalysis(confMan, mtdcheader, blMan, bltrackMan, cdsMan, cdstrackMan, particle)){
    //std::cout << " HeKppinAna::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef OverVeto
  if(!OverVetoAna->DoAnalysis(confMan, mtdcheader, blMan, bltrackMan, cdsMan, cdstrackMan, particle)){
    //std::cout << " OverVetoAna::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef CDSCalib
  if(!CDSCalibAna->DoAnalysis(confMan, mtdcheader, blMan, bltrackMan, cdsMan, cdstrackMan, particle)){
    //std::cout << " CDSCalibAna::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef NCCalib
  if(!NCCalibAna->DoAnalysis(confMan, mtdcheader, blMan, bltrackMan, cdsMan, cdstrackMan, particle)){
    //std::cout << " NCCalibAna::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef TwoTrack
  if(!TwoTrackAna->DoAnalysis(confMan, mtdcheader, blMan, bltrackMan, cdsMan, cdstrackMan, particle)){
    //std::cout << " TwoTrackAna::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef Lambda
  if(!LambdaAna->DoAnalysis(confMan, mtdcheader, blMan, bltrackMan, cdsMan, cdstrackMan, particle)){
    //std::cout << " LambdaAna::false " << std::endl;
    //Clear();
    //return true;
  }
#endif

  //evTree->Fill();

  Clear();
  return true;
}

void EventAnalysisMyAna::Clear()
{
#ifdef Basic
  basicAna->Clear();
#endif
#ifdef K0
  k0Ana->Clear();
#endif
#ifdef K0pn
  k0pnAna->Clear();
#endif
#ifdef Kpn
  kpnAna->Clear();
#endif
#ifdef K0n
  k0nAna->Clear();
#endif
#ifdef pKn
  pKnAna->Clear();
#endif
#ifdef dKn
  dKnAna->Clear();
#endif
#ifdef dKnpipi
  dKnpipiAna->Clear();
#endif
#ifdef dKnppim
  dKnppimAna->Clear();
#endif
#ifdef HeKn
  HeKnAna->Clear();
#endif
#ifdef HeKFWD
  HeKFWDAna->Clear();
#endif
#ifdef HeKpipip
  HeKpipipAna->Clear();
#endif
#ifdef HeKpipipp
  HeKpipippAna->Clear();
#endif
#ifdef HeKFWDNeutral
  HeKFWDNeutralAna->Clear();
#endif
#ifdef HeK2CDS
  HeK2CDSAna->Clear();
#endif
#ifdef CDCSlew
  CDCSlewAna->Clear();
#endif
#ifdef CDCTrackingwSlew
  CDCTrackingwSlewAna->Clear();
#endif
#ifdef dKFWD
  dKFWDAna->Clear();
#endif
#ifdef dKCDS
  dKCDSAna->Clear();
#endif
#ifdef dKCDH3
  dKCDH3Ana->Clear();
#endif
#ifdef dKpipiFWD
  dKpipiFWDAna->Clear();
#endif
#ifdef Check
  CheckAna->Clear();
#endif
#ifdef HeKpipp
  HeKpippAna->Clear();
#endif
#ifdef HeKpippOneDummy
  HeKpippOneDummyAna->Clear();
#endif
#ifdef HeKpippOneDummyLpair
  HeKpippOneDummyLpairAna->Clear();
#endif
#ifdef HeKpippMixEvent
  HeKpippMixEventAna->Clear();
#endif
#ifdef HeKLpnFWD
  HeKLpnFWDAna->Clear();
#endif
#ifdef CheckChi2
  CheckChi2Ana->Clear();
#endif
#ifdef HeKpFWD
  HeKpFWDAna->Clear();
#endif
#ifdef dKppipi
  dKppipiAna->Clear();
#endif
#ifdef HeKpn
  HeKpnAna->Clear();
#endif
#ifdef HeKppin
  HeKppinAna->Clear();
#endif
#ifdef OverVeto
  OverVetoAna->Clear();
#endif
#ifdef CDSCalib
  CDSCalibAna->Clear();
#endif
#ifdef NCCalib
  NCCalibAna->Clear();
#endif
#ifdef TwoTrack
  TwoTrackAna->Clear();
#endif
#ifdef Lambda
  LambdaAna->Clear();
#endif
  blMan->Clear();
  bltrackMan->Clear();
  cdsMan->Clear();
  header->Clear();
}

void EventAnalysisMyAna::Finalize()
{
  std::cout << " Enter EventAnalysisMyAna::Finalize " << std::endl;

  rtFile->cd();
  gFile->Write();
  gFile->Close();

#ifdef Basic
  delete basicAna;
#endif
#ifdef K0
  delete k0Ana;
#endif
#ifdef K0pn
  delete k0pnAna;
#endif
#ifdef Kpn
  delete kpnAna;
#endif
#ifdef K0n
  delete k0nAna;
#endif
#ifdef pKn
  delete pKnAna;
#endif
#ifdef dKn
  delete dKnAna;
#endif
#ifdef dKnpipi
  delete dKnpipiAna;
#endif
#ifdef dKnppim
  delete dKnppimAna;
#endif
#ifdef HeKn
  delete HeKnAna;
#endif
#ifdef HeKFWD
  delete HeKFWDAna;
#endif
#ifdef HeKpipip
  delete HeKpipipAna;
#endif
#ifdef HeKpipipp
  delete HeKpipippAna;
#endif
#ifdef HeKFWDNeutral
  delete HeKFWDNeutralAna;
#endif
#ifdef HeK2CDS
  delete HeK2CDSAna;
#endif
#ifdef CDCSlew
  delete CDCSlewAna;
#endif
#ifdef CDCTrackingwSlew
  delete CDCTrackingwSlewAna;
#endif
#ifdef dKFWD
  delete dKFWDAna;
#endif
#ifdef dKCDS
  delete dKCDSAna;
#endif
#ifdef dKCDH3
  delete dKCDH3Ana;
#endif
#ifdef dKpipiFWD
  delete dKpipiFWDAna;
#endif
#ifdef Check
  delete CheckAna;
#endif
#ifdef HeKpipp
  delete HeKpippAna;
#endif
#ifdef HeKpippOneDummy
  delete HeKpippOneDummyAna;
#endif
#ifdef HeKpippOneDummyLpair
  delete HeKpippOneDummyLpairAna;
#endif
#ifdef HeKpippMixEvent
  delete HeKpippMixEventAna;
#endif
#ifdef HeKLpnFWD
  delete HeKLpnFWDAna;
#endif
#ifdef CheckChi2
  delete CheckChi2Ana;
#endif
#ifdef HeKpFWD
  delete HeKpFWDAna;
#endif
#ifdef dKppipi
  delete dKppipiAna;
#endif
#ifdef HeKpn
  delete HeKpnAna;
#endif
#ifdef HeKppin
  delete HeKppinAna;
#endif
#ifdef OverVeto
  delete OverVetoAna;
#endif
#ifdef CDSCalib
  delete CDSCalibAna;
#endif
#ifdef NCCalib
  delete NCCalibAna;
#endif
#ifdef TwoTrack
  delete TwoTrackAna;
#endif
#ifdef Lambda
  delete LambdaAna;
#endif
  delete blMan;
  delete bltrackMan;
  delete cdsMan;
  delete header;

}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisMyAna *event = new EventAnalysisMyAna();
  return (EventTemp*)event;
}
