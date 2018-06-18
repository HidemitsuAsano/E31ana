#include "File.h"
#include "ConfMan.h"
#include "TKO.h"
#include "CDSHitMan.h"
#include "CDSTrackingMan.h"
#include "CircleFit.h"
#include "HelixFit.h"
#include "BeamLineHitMan.h"
#include "BeamLineTrackMan.h"
#include "BeamSpectrometer.h"
#include "EventHeader.h"
#include "ScalerMan.h"

#include "EventAlloc.h"
#include "EventTemp.h"
#include "Display.h"

#include "TrackTools.h"
#include "Tools.h"

using namespace std;

#define HISTOGRAM 1
static const int STOP_FILL_NUM=20000;
//static const int STOP_FILL_NUM=3000;

static const double BLDC_TIME_WINDOW_MIN=-30;
static const double BLDC_TIME_WINDOW_MAX=100;

static const double TOF_K_MIN=27.8588;
static const double TOF_K_MAX=39.5663;

static const double BHD_MATCH_MIN[10]={ 0.970692, 0.9754, 0.981274, 0.986733, 0.992045,
                                        0.99713, 1.00219, 1.00708, 1.01191, 1.01588 };

static const double BHD_MATCH_MAX[10]={ 1.01324, 1.01846, 1.02366, 1.0282, 1.0321,
                                        1.03563, 1.03903, 1.0422, 1.04523, 1.04518 };

class EventAnalysis: public EventTemp
{
public:
  EventAnalysis();
  ~EventAnalysis();
private:
  TFile *rtFile;
  TTree *evTree;
  TTree *scaTree;
  BeamLineHitMan *blMan;
  BeamLineTrackMan *bltrackMan;
  BeamSpectrometer *beamSpec;
  CDSHitMan *cdsMan;
  CDSTrackingMan *trackMan;
  EventHeader *header;
  int t0,t1;
  TString tmpname;

  int AllGoodTrack;
  int nTrack;
  int nFill;

public:
  void Initialize( ConfMan *conf );
  void USca( int nsca, unsigned int *sca );
  bool UAna( TKOHitCollection *tko );
  bool Clear();
  void Finalize();
  void InitializeHistogram();
};

EventAnalysis::EventAnalysis()
  : EventTemp()
{
}

EventAnalysis::~EventAnalysis()
{
}

void EventAnalysis::Initialize( ConfMan *conf )
{
#if 1
  std::cout << " Enter EventAnalysis::Initialize " << std::endl;
#endif
  confMan = conf;
  rtFile = new TFile( confMan->GetOutFileName().c_str() , "recreate" );
  rtFile->cd();

  evTree = new TTree( "EventTree", "EventTree" );
  header = new EventHeader();
  cdsMan = new CDSHitMan();
  trackMan = new CDSTrackingMan();  
  blMan = new BeamLineHitMan();
  bltrackMan = new BeamLineTrackMan();
  evTree->Branch( "BeamLineTrackMan", &bltrackMan );
  beamSpec = new BeamSpectrometer(conf);

  t0=time(0);
  AllGoodTrack=0;
  nTrack=0;
  nFill=0;

#if HISTOGRAM
  std::cout<<" > Histogram mode"<<std::endl;
  InitializeHistogram();
#else
  std::cout<<" > Make tree mode"<<std::endl;
  evTree->Branch( "EventHeader", &header );
  evTree->Branch( "CDSTrackingMan", &trackMan );
#endif
}

void EventAnalysis::USca( int nsca, unsigned int *sca )
{
#if 0
  std::cout << " Enter EventAnalysis::USca " << std::endl;
#endif
  Block_Event_Number++;
  header->SetBlockEventNumber( Block_Event_Number );
}

bool EventAnalysis::UAna( TKOHitCollection *tko )
{
#if 0
  std::cout << " Enter EventAnalysis::UAna " << std::endl;
#endif

  Event_Number++;
  { int status = confMan->CheckEvNum( Event_Number, Block_Event_Number );
    if( status==1 ) return true;
    if( status==2 ) return false; }

  if( Event_Number%100==0)
    {
      t1=time(0);
      std::cout << " Event# : " << Event_Number << "  BlockEvent# : " << Block_Event_Number << " AllTrack# : " << nTrack <<" GoodTrack# : " << AllGoodTrack << " Time (s): " << (t1-t0) << "  nFill : "<<nFill<<std::endl;
    }
  
  header->SetRunNumber(0);
  header->SetEventNumber(Event_Number);
  
  header->Convert( tko, confMan );

  if( header->IsTrig(Trig_Cosmic) ) return Clear();

  blMan-> Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );
  HodoscopeLikeHit *T0_hit=0;
  int nT0=0, nBHD=0;
  for( int i=0; i<blMan->nT0(); i++ ){
    if( blMan->T0(i)->CheckRange() ){
      nT0++;
      T0_hit = blMan->T0(i);
    }
  }
  for( int i=0; i<blMan->nBHD(); i++ ){
    if( blMan->BHD(i)->CheckRange() ) nBHD++;
  }
  if( nT0!=1 || nBHD==0 ) return Clear();

  int nCDH=0;
  for( int i=0; i<cdsMan->nCDH(); i++ ){
    if( cdsMan->CDH(i)->CheckRange() ) nCDH++;
  }
  if( nCDH==0 ) return Clear();

  bltrackMan-> DoTracking(blMan, confMan, true, true);
  int ntrackBPC=0, ntrackBLC1=0, ntrackBLC2=0;
  LocalTrack *BPC=0, *BLC1=0, *BLC2=0;
  for( int i=0; i<bltrackMan->ntrackBLC1(); i++ ){
    LocalTrack *track = bltrackMan->trackBLC1(i);
    if( BLDC_TIME_WINDOW_MIN<track->GetTrackTime() && track->GetTrackTime()<BLDC_TIME_WINDOW_MAX ){
      ntrackBLC1++;
      BLC1=track;
    }
  }
  for( int i=0; i<bltrackMan->ntrackBLC2(); i++ ){
    LocalTrack *track = bltrackMan->trackBLC2(i);
    if( BLDC_TIME_WINDOW_MIN<track->GetTrackTime() && track->GetTrackTime()<BLDC_TIME_WINDOW_MAX ){
      ntrackBLC2++;
      BLC2=track;
    }
  }
  for( int i=0; i<bltrackMan->ntrackBPC(); i++ ){
    LocalTrack *track = bltrackMan->trackBPC(i);
    if( BLDC_TIME_WINDOW_MIN<track->GetTrackTime() && track->GetTrackTime()<BLDC_TIME_WINDOW_MAX ){
      ntrackBPC++;
      BPC=track;
    }
  }
  if( ntrackBPC!=1 || ntrackBLC1!=1 || ntrackBLC2!=1 ) return Clear();
  if( nCDH==0 ) return Clear();

  beamSpec-> TMinuitFit(BLC1, BLC2, confMan);
  double beam_mom = beamSpec->mom();
  int beam_pid=Beam_Kaon;
  double beam_mass = kpMass;
  if( header->IsTrig(Trig_Pion) ){
    beam_pid=Beam_Pion;
    beam_mass=piMass;

  }
  bool goodK=false;
  for( int i=0; i<blMan->nBHD(); i++ ){
    if( blMan->BHD(i)->CheckRange() ){
      double tof = T0_hit->ctmean()-blMan->BHD(i)->ctmean();
      if( TOF_K_MIN<tof && tof<TOF_K_MAX ) goodK=true;
    }
  }

  bool goodBLC1=false, goodBLC2=false, goodBPC=false, goodD5=false;
  if( -10<BPC->GetTrackTime()  && BPC->GetTrackTime()<10  && BPC->chi2all()<10  ) goodBPC=true;
  if( -10<BLC1->GetTrackTime() && BLC1->GetTrackTime()<10 && BLC1->chi2all()<10 ) goodBLC1=true;
  if( -10<BLC2->GetTrackTime() && BLC2->GetTrackTime()<10 && BLC2->chi2all()<10 ) goodBLC2=true;
  if( beamSpec->chisquare()<20 ){
    for( int i=0; i<blMan->nBHD(); i++ ){
      if( blMan->BHD(i)->CheckRange() ){
	int index = blMan-> BHD(i)->seg()-6;
	if( BHD_MATCH_MIN[index]<beam_mom && beam_mom<BHD_MATCH_MAX[index] ) goodD5=true;
      }
    }
  }
  TVector3 T0pos=BPC->GetPosatZ(-110.5);

  trackMan->Execute(cdsMan,confMan);

  for( int i=0; i<trackMan->nGoodTrack(); i++ ){
    CDSTrack *track = trackMan-> GoodTrack(i);
    int CDHseg;
    double CDHtime;
    if( !track-> GetCDHHit(cdsMan, CDHseg, CDHtime ) ) continue;
    TVector3 vtxBeam, vtxCDS;
    if( !track-> GetVertex(T0pos, BPC-> GetMomDir(), vtxBeam, vtxCDS) ) continue;
    double beam_mom_vtx, beam_tof;
    ELossTools::CalcElossBeamTGeo(T0pos, vtxBeam, beam_mom, beam_mass, beam_mom_vtx, beam_tof);

    double mass2, beta, tof;
    if( !TrackTools::FindMass2(track, BPC, CDHtime-T0_hit->ctmean(), beam_mom, beam_pid, beta, mass2, tof) ) continue;

    bool single=true;
    for( int layer=1; layer<=NumOfCDCLayers; layer++ ){
      if( track->nTrackHit(layer)!=1 ) single=false;
    }
    if( single ){
      for( int layer=1; layer<=NumOfCDCLayers; layer++ ){
	CDCHit *cdc = track->TrackHit(cdsMan, layer, 0);
	int wire=cdc->wire();
	double res=cdc->resl();
	double dt=cdc->dt();
	double dl=cdc->dl();
	double dlr=dl-res;
	double dxdt=confMan-> GetXTMapManager()-> CalcDxDt(CID_CDC, layer, 0, dt);

	Tools::Fill1D(Form("CDC_res_%d", layer), res);
	Tools::Fill1D(Form("CDC_res_%d_%d", layer, wire), res);
	Tools::Fill2D(Form("CDC_dt_dl_%d", layer), dt, dlr);
	Tools::Fill2D(Form("CDC_dt_res_%d", layer), dt, res);
	Tools::Fill2D(Form("CDC_dt_dl_%d_%d", layer, wire), dt, dlr);
	Tools::Fill2D(Form("CDC_dt_res_%d_%d", layer, wire), dt, res);
      }
    }

    track-> Retiming(cdsMan, confMan, beta, true);
    for( int i=0; i<3; i++ ) track-> HelixFitting(cdsMan);
  
#if HISTOGRAM
    if( goodBLC1 && goodBLC2 && goodBPC && goodK && goodD5 && track->Chi()<30 ){
      bool single=true;
      for( int layer=1; layer<=NumOfCDCLayers; layer++ ){
	if( track->nTrackHit(layer)!=1 ) single=false;
      }
      if( single ){
	for( int layer=1; layer<=NumOfCDCLayers; layer++ ){
	  CDCHit *cdc = track->TrackHit(cdsMan, layer, 0);
	  int wire=cdc->wire();
	  double res=cdc->resl();
	  double dt=cdc->dt();
	  double dl=cdc->dl();
	  double dlr=dl-res;
	  double dxdt=confMan-> GetXTMapManager()-> CalcDxDt(CID_CDC, layer, wire, dt);

	  Tools::Fill1D(Form("CDC_res2_%d", layer), res);
	  Tools::Fill1D(Form("CDC_res2_%d_%d", layer, wire), res);
	  Tools::Fill2D("CDC_ob2_res", 1/(beta*beta), res/dxdt);
	  Tools::Fill2D(Form("CDC_ob2_res_%d", layer), 1/(beta*beta), res/dxdt);
	  Tools::Fill2D(Form("CDC_ob2_res_%d_%d", layer, wire), 1/(beta*beta), res/dxdt);
	}
	nFill++;
      }
      TVector3 vertexCDH=track->CDHVertex();
      HodoscopeLikeHit *CDH=track-> CDHHit(cdsMan, 0);

      Tools::Fill1D("CDC_Z", vertexCDH.Z());
      Tools::Fill2D(Form("CDC_Z_CDH%d_tsub", CDH->seg()), vertexCDH.Z(), CDH->ctsub());
    }
#endif
  }

  int nGoodTrack=trackMan->nGoodTrack();
  int nallTrack=  trackMan->nTrack();
  AllGoodTrack+=nGoodTrack;
  nTrack+=nallTrack;

#ifndef HISTOGRAM
  if( nallTrack<1 ) return Clear();
  evTree->Fill();
#endif
#ifdef HISTOGRAM
  if( STOP_FILL_NUM<nFill ){
    Clear();
    return false;
  }
#endif
  return Clear();
}

void EventAnalysis::Finalize()
{
  std::cout << " Enter EventAnalysis::Finalize    Event_Number : "<<Event_Number<< std::endl;

  rtFile->cd();
  confMan->SaveParams();
  //  confMan-> Write("ConfMan");
  gFile->Write();
  gFile->Close();

  delete cdsMan;
  delete trackMan;
  delete header;
}

bool EventAnalysis::Clear()
{
  header->Clear();
  blMan-> Clear();
  bltrackMan->Clear();
  cdsMan->Clear();
  trackMan->Clear();
  beamSpec->Clear();

  return true;
}

void EventAnalysis::InitializeHistogram()
{
  rtFile-> cd();
  Tools::newTH1F("CDC_Z", 1000, -100, 100);
  for( int seg=1; seg<=36; seg++ ){
    Tools::newTH2F(Form("CDC_Z_CDH%d_tsub", seg), 1000, -100, 100, 500, -5, 5);
  }


  Tools::newTH2F("CDC_ob2_res", 500, 0, 100, 500, -25, 25);
  for( int layer=1; layer<=NumOfCDCLayers; layer++ ){
    Tools::newTH1F(Form("CDC_res_%d", layer), 1000, -1.0, 1.0);
    Tools::newTH1F(Form("CDC_res2_%d", layer), 1000, -1.0, 1.0);
    Tools::newTH2F(Form("CDC_dt_dl_%d", layer), 200, -20, 380, 200, -0.1, 1.1);
    Tools::newTH2F(Form("CDC_dt_res_%d", layer), 200, -20, 380, 200, -0.25, 0.25);
    Tools::newTH2F(Form("CDC_ob2_res_%d", layer), 200, 0, 100, 200, -25, 25);

    for( int wire=1; wire<=NumOfCDCWiresInLayer[layer-1]; wire++ ){
      Tools::newTH1F(Form("CDC_res_%d_%d", layer, wire), 1000, -1.0, 1.0);
      Tools::newTH1F(Form("CDC_res2_%d_%d", layer, wire), 1000, -1.0, 1.0);
      Tools::newTH2F(Form("CDC_dt_dl_%d_%d", layer, wire), 200, -20, 380, 200, -0.1, 1.1);
      Tools::newTH2F(Form("CDC_dt_res_%d_%d", layer, wire), 200, -20, 380, 200, -0.25, 0.25);
      Tools::newTH2F(Form("CDC_ob2_res_%d_%d", layer, wire), 200, 0, 100, 200, -25, 25);
    }
  }
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysis *event = new EventAnalysis();
  return (EventTemp*)event;
}
