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

#define Debug 0


class EventAnalysisMyCDCTrackingwSlew: public EventTemp
{
  public:
    EventAnalysisMyCDCTrackingwSlew();
    ~EventAnalysisMyCDCTrackingwSlew();
  private:
    TFile *rtFile;
    TTree *evTree;
    TTree *scaTree;
    TFile *cdcFile;
    TTree *cdcTree;
    TFile *mtdcFile;
    TTree *mtdcTree;
    CDSHitMan *cdsMan;
    CDSTrackingMan *trackMan;
    CDSTrackingMan *cdstrackMan;
    EventHeader *cdsheader;
    BeamLineHitMan *blMan;
    BeamLineTrackMan *bltrackMan;
    EventHeader *header;
    EventHeader *mtdcheader;
    Particle* particle;
    ScalerMan *scaMan;

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

  EventAnalysisMyCDCTrackingwSlew::EventAnalysisMyCDCTrackingwSlew()
: EventTemp()
{
}

EventAnalysisMyCDCTrackingwSlew::~EventAnalysisMyCDCTrackingwSlew()
{
}

const int MaxTreeSize = 1900000000000;
void EventAnalysisMyCDCTrackingwSlew::Initialize( ConfMan *conf )
{
#if 1
  std::cout << " Enter EventAnalysisMyCDCTrackingwSlew::Initialize " << std::endl;
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
  evTree = new TTree( "EventTree", "EventTree" );
  header = new EventHeader();
  if( header==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  evTree->Branch( "EventHeader", &header );
  cdsMan = new CDSHitMan();
  if( cdsMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  blMan = new BeamLineHitMan();
  if( blMan ==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  bltrackMan = new BeamLineTrackMan();
  if( bltrackMan ==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  trackMan = new CDSTrackingMan();  
  if( trackMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  evTree->Branch( "CDSTrackingMan", &trackMan );

  scaMan = new ScalerMan();
  if( scaMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }

  AllGoodTrack = 0;
  nAllTrack = 0;
  CDC_Event_Number = 0;
  MTDC_Event_Number = 0;
  t0=time(0);
}

void EventAnalysisMyCDCTrackingwSlew::USca( int nsca, unsigned int *sca )
{
  Block_Event_Number++;
  header->SetBlockEventNumber( Block_Event_Number );
  scaMan->SetBlockEventNumber( Block_Event_Number );
  for( int i=0; i<nsca; i++ ){
    scaMan->AddHit( sca[i], confMan->GetScalerMapManager()->GetName(i) );
  }

  scaMan->Clear();
}

bool EventAnalysisMyCDCTrackingwSlew::UAna( TKOHitCollection *tko )
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

  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );
  cdstrackMan -> Calc( cdsMan, confMan );

  rtFile->cd();

  if(particle->nBeam()!=1) return false;
  pBeam* beam = particle->beam(0);

  // ###################### //
  // CDS all track analysis //
  // ###################### //
  for( int id=0; id<cdstrackMan->nGoodTrack(); id++ ){
    CDSTrack* cdc=cdstrackMan->GoodTrack(id);

    if(cdc->nCDHHit()!=1) continue;
    double cdhtime=0,dis=0;
    int cdhseg=0;
    if(!cdc->GetCDHHit(cdsMan,cdhseg,cdhtime)) continue;
    HodoscopeLikeHit* cdh=0;
    for(int i=0; i<cdsMan->nCDH(); i++){
      HodoscopeLikeHit* hit = cdsMan->CDH(i);
      if(hit->CheckRange()){
        if(hit->seg()==cdhseg) cdh=hit;
      }
    }
    double de = cdh->emean();

    TVector3 cdsvtx,beamvtx;
    double par[5];
    cdc->GetParameters(CID_CDC,par,cdsvtx);
    double mom=cdc->Momentum();
    if(!cdc->GetVertex(beam->t0pos(),beam->bpcdir(),beamvtx,cdsvtx)) continue;
    dis = (beamvtx-cdsvtx).Mag();

    double beam_out=0, beam_tof=0;
    ELossTools::CalcElossBeamTGeo(beam->t0pos(),beamvtx,beam->mom(),kpMass,beam_out,beam_tof);
    TVector3 cdhvtx=cdc->CDHVertex();
    double cdc_dis=MathTools::CalcHelixArc(par,cdhvtx,cdsvtx);
    double beta=cdc_dis/(cdhtime-beam->t0time()-beam_tof)/(Const*100);
    double mass2=mom*mom*(1/(beta*beta)-1);

    double calc_beta,tofvtxcdc;
    if(!TrackTools::FindMass2(cdc,beam,cdhtime-beam->t0time(),calc_beta,mass2,tofvtxcdc)) continue;
    beta = calc_beta;

    cdc->Retiming(cdsMan,confMan,beta,true);
  }

  int nTrack=cdstrackMan->nTrack();

  trackMan = cdstrackMan;

  if( nTrack<1 ){
    return true;
  }
  else{
    evTree->Fill();
  }

  Clear();
  return true;
}

void EventAnalysisMyCDCTrackingwSlew::Clear()
{
  blMan->Clear();
  bltrackMan->Clear();
  cdsMan->Clear();
  header->Clear();
  cdstrackMan->Clear();
  trackMan->Clear();
}

void EventAnalysisMyCDCTrackingwSlew::Finalize()
{
  std::cout << " Enter EventAnalysisMyCDCTrackingwSlew::Finalize " << std::endl;

  rtFile->cd();
  gFile->Write();
  gFile->Close();

  delete blMan;
  delete bltrackMan;
  delete cdsMan;
  delete cdstrackMan;
  //delete trackMan;
  delete header;
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisMyCDCTrackingwSlew *event = new EventAnalysisMyCDCTrackingwSlew();
  return (EventTemp*)event;
}
