#include "File.h"
#include "ConfMan.h"
#include "TKO.h"
#include "CDSHitMan.h"
#include "BeamLineHitMan.h"
#include "EventHeader.h"
#include "ScalerMan.h"
#include "BeamLineTrackMan.h"
//#include "HistManBeamAna.h"
#include "CDSTrackingMan.h"

#include "EventAlloc.h"
#include "EventTemp.h"
#include "BeamSpectrometer.h"

#include "HistManwMC.h"
#include "HistTools.h"

#define W_TRIG 1
#define BLDC_CHI_CUT 1
#define BLDC_WIDE_TIMING 1

static const double BLDC_CHI2_MAX=30;
static const double TOF_K_MIN=27.8588;
static const double TOF_K_MAX=39.5663;
static const double TOF_PI_MIN=24.5606;
static const double TOF_PI_MAX=26.9866;
static const double BLC2BPC_X_MIN=-0.72645;
static const double BLC2BPC_X_MAX=0.770232;
static const double BLC2BPC_Y_MIN=-0.778481;
static const double BLC2BPC_Y_MAX=0.755978;
static const double BLC2BPC_dX_MIN=-0.0202092;
static const double BLC2BPC_dX_MAX=0.0201483;
static const double BLC2BPC_dY_MIN=-0.0200048;
static const double BLC2BPC_dY_MAX=0.0204583;

static const double BHD_MATCH_MIN[10]={ 0.970692, 0.9754, 0.981274, 0.986733, 0.992045, 
					0.99713, 1.00219, 1.00708, 1.01191, 1.01588 };

static const double BHD_MATCH_MAX[10]={ 1.01324, 1.01846, 1.02366, 1.0282, 1.0321,
					1.03563, 1.03903, 1.0422, 1.04523, 1.04518 };

class EventAnalysisNpipi: public EventTemp
{
public:
  EventAnalysisNpipi();
  ~EventAnalysisNpipi();

private:
  int Event_Number_onSpill;
  int CDC_Event_Num;
  TFile *cdcFile;
  TTree *cdcTree;
  CDSTrackingMan *cdstrackMan;
  EventHeader *cdcHeader;

  TTree *scaTree;

  TFile *rtFile;
  CDSHitMan *cdsMan;
  BeamLineHitMan *blMan;
  EventHeader *header;
  ScalerMan *scaMan;

  BeamLineTrackMan *bltrackMan;
  BeamSpectrometer *beamSpec;

  HistManwMC *histMan;

public:
  void Initialize( ConfMan *conf );
  void InitializeHistogram();
  void USca( int nsca, unsigned int *sca );
  bool UAna( TKOHitCollection *tko );
  void Clear();
  void Finalize();
};

EventAnalysisNpipi::EventAnalysisNpipi()
  : EventTemp()
{
  std::cout<<"EventAnalysisNpipi::Constractor Call"<<std::endl;
}

EventAnalysisNpipi::~EventAnalysisNpipi()
{
}

const int MaxTreeSize = 1900000000000;
void EventAnalysisNpipi::Initialize( ConfMan *conf )
{
  std::cout<<"EventAnalysisNpipi::Initialization START"<<std::endl;
  Event_Number_onSpill = 0;
  CDC_Event_Num = 0;
  confMan = conf;
  cdcFile = new TFile((TString)confMan->GetCDSTrackFileName());
  cdstrackMan = new CDSTrackingMan();
  cdcHeader = new EventHeader();
  cdcTree = (TTree*)cdcFile-> Get("EventTree");
  cdcTree-> SetBranchAddress("EventHeader", &cdcHeader);
  cdcTree-> SetBranchAddress("CDSTrackingMan", &cdstrackMan);
  
  rtFile = new TFile( confMan->GetOutFileName().c_str(), "recreate" );
  scaTree = new TTree("ScalerTree", "ScalerTree");
  scaMan = new ScalerMan();
  scaTree-> Branch("ScalerMan", &scaMan);

  header = new EventHeader();
  cdsMan = new CDSHitMan();
  blMan = new BeamLineHitMan();
  bltrackMan = new BeamLineTrackMan();
  beamSpec = new BeamSpectrometer(confMan);

  histMan = new HistManwMC(rtFile, confMan);
  histMan-> setK18ana(blMan, cdsMan, bltrackMan, cdstrackMan);
  histMan-> initHist();

  InitializeHistogram();
  std::cout<<"EventAnalysisNpipi::Initialization FINISH"<<std::endl;
}

void EventAnalysisNpipi::USca( int nsca, unsigned int *sca )
{
  Block_Event_Number++;
  header->SetBlockEventNumber( Block_Event_Number );
  scaMan->SetBlockEventNumber( Block_Event_Number );
  for( int i=0; i<nsca; i++ ){
   scaMan->AddHit( sca[i], confMan->GetScalerMapManager()->GetName(i) );
  }

  scaTree-> Fill();
  scaMan->Clear();
}

bool EventAnalysisNpipi::UAna( TKOHitCollection *tko )
{
  TH1F *h1;
  TH2F *h2;

  Event_Number++;
  { int status = confMan->CheckEvNum( Event_Number, Block_Event_Number );
    if( status==1 ) return true;
    if( status==2 ) return false; }

  if( Event_Number%5000==0 )
    std::cout << " Event# : " << Event_Number << "  BlockEvent# : " << Block_Event_Number << std::endl;

  TH1F *h1_ER = (TH1F*)gFile-> Get("EventReduction");
  h1_ER-> Fill(-1);
    
  header->SetRunNumber(0);
  header->SetEventNumber(Event_Number);
  
  header->Convert( tko, confMan );
  if( header->IsTrig(Trig_Cosmic) ){
    //    std::cout<<" Error Trig"<<std::endl;
    Clear();
    return true;
  }
  Event_Number_onSpill++;

  TH1F *h1_Kf = (TH1F*)gFile-> Get("Kf_Reduction");
  TH1F *h1_N  = (TH1F*)gFile-> Get("N_Reduction");

  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );

  h1_ER-> Fill(0);
  if( header->IsTrig(Trig_Kf)      ) h1_Kf-> Fill(0);
  if( header->IsTrig(Trig_Neutral) ) h1_N-> Fill(0);

  //***** Event Analysis Tmep Check **************************************************************************//
  bltrackMan-> DoTracking(blMan, confMan, true, true);

  int nlayBLC1a=0;
  int nlayBLC1b=0;
  int nlayBLC2a=0;
  int nlayBLC2b=0;
  int nlayBPC=0;
  for( int l=1; l<=8; l++ ){
    if( blMan->nBLC1a(l)>0 ) nlayBLC1a++;
    if( blMan->nBLC1b(l)>0 ) nlayBLC1b++;
    if( blMan->nBLC2a(l)>0 ) nlayBLC2a++;
    if( blMan->nBLC2b(l)>0 ) nlayBLC2b++;
    if( blMan->nBPC(l)>0 ) nlayBPC++;
  }
#if 0
  std::cout<<"BLC1a nlay:"<<nlayBLC1a<<" ntrack:"<<bltrackMan->ntrackBLC1a()<<std::endl;
  std::cout<<"BLC1b nlay:"<<nlayBLC1b<<" ntrack:"<<bltrackMan->ntrackBLC1b()<<std::endl;
  std::cout<<"BLC1            ntrack:"<<bltrackMan->ntrackBLC1()<<std::endl;
  std::cout<<"BLC2a nlay:"<<nlayBLC1a<<" ntrack:"<<bltrackMan->ntrackBLC2a()<<std::endl;
  std::cout<<"BLC2b nlay:"<<nlayBLC1b<<" ntrack:"<<bltrackMan->ntrackBLC2b()<<std::endl;
  std::cout<<"BLC2            ntrack:"<<bltrackMan->ntrackBLC1()<<std::endl;
  std::cout<<"BPC   nlay:"<<nlayBPC<<" ntrack:"<<bltrackMan->ntrackBPC()<<std::endl;
#endif

  //*** Check T0 1hit ***//
  int nT0=0;
  double T0time = DBL_MIN;
  for( int i=0; i<blMan->nT0(); i++ ){
    if( blMan->T0(i)->CheckRange() ){
      nT0++;
      T0time = blMan->T0(i)->ctmean();
    }
  }
  if( header->IsTrig(Trig_Kf) ){
    h1 = (TH1F*)gFile-> Get("nT0"), h1->Fill(nT0);
  }

  if( nT0!=1 ){
    Clear();
    return true;
  }

  h1_ER-> Fill(1);
  if( header->IsTrig(Trig_Kf)      ) h1_Kf-> Fill(1);
  if( header->IsTrig(Trig_Neutral) ) h1_N-> Fill(1);

  //*** Check BHD-T0 TOF ***//
  bool kaon_flag = false;
  bool pion_flag = false;
  int nBHD=0;
  int BHDseg=-1;
  for( int i=0; i<blMan->nBHD(); i++ ){
    if( blMan->BHD(i)->CheckRange() ){
      if( 5<blMan->BHD(i)->seg() && blMan-> BHD(i)->seg()<16 ){
	nBHD++;
	BHDseg = blMan->BHD(i)->seg();
	if( nT0==1 ){
	  double tof = T0time-blMan->BHD(i)->ctmean();
	  h1 = (TH1F*)gFile-> Get("BHDT0"), h1-> Fill(tof);
	  if( header->IsTrig(Trig_Kf) ) h1 = (TH1F*)gFile-> Get("BHDT0_Kf"), h1-> Fill(tof);
	  if( header->IsTrig(Trig_Kaon) ) h1 = (TH1F*)gFile-> Get("BHDT0_K"), h1-> Fill(tof);
	  if( header->IsTrig(Trig_Pion) ) h1 = (TH1F*)gFile-> Get("BHDT0_pi"), h1-> Fill(tof);
	  if( TOF_K_MIN<tof && tof<TOF_K_MAX ) kaon_flag = true;
	  if( TOF_PI_MIN<tof && tof<TOF_PI_MAX ) pion_flag = true;
	}
      }
    }
  }
  if( header->IsTrig(Trig_Kf) ){
    h1 = (TH1F*)gFile-> Get("nBHD"), h1->Fill(nBHD);
  }

  histMan->setBeamKaon(kaon_flag);
  histMan->setBeamPion(pion_flag);

  int beam_pid=Beam_Other;
#if W_TRIG
  if( header->IsTrig(Trig_Kaon) ){
    if( kaon_flag ) beam_pid=Beam_Kaon;
  } 
  else if( header->IsTrig(Trig_Pion) ){
    if( pion_flag ) beam_pid=Beam_Pion;
  }
  histMan->setBeamPID(beam_pid);
  if( beam_pid!=Beam_Kaon && beam_pid!=Beam_Pion ){
    Clear();
    return true;
  }
  if( beam_pid==Beam_Kaon ){
    h1_ER-> Fill(2);
    if( header->IsTrig(Trig_Kf)      ) h1_Kf-> Fill(2);
    if( header->IsTrig(Trig_Neutral) ) h1_N-> Fill(2);
  }
#else
  if( pion_flag ) beam_pid=Beam_Pion;
  if( kaon_flag ) beam_pid=Beam_Kaon;
  histMan->setBeamPID(beam_pid);
  if( beam_pid!=Beam_Kaon ){
    Clear();
    return true;
  }
  if( beam_pid==Beam_Kaon ){
    h1_ER-> Fill(2);
    if( header->IsTrig(Trig_Kf)      ) h1_Kf-> Fill(2);
    if( header->IsTrig(Trig_Neutral) ) h1_N-> Fill(2);
  }
#endif

  //  bltrackMan-> DoTracking2(blMan, confMan);

  int ntrackBLC1=0;
  int ntrackBLC1_1=0;
  int ntrackBLC1_2=0;
  LocalTrack *trackBLC1=0;
  for( int i=0; i<bltrackMan->ntrackBLC1(); i++ ){
    double tracktime = bltrackMan->trackBLC1(i)->GetTrackTime();
    if( beam_pid==Beam_Kaon ) h1 = (TH1F*)gFile-> Get("BLC1_time"), h1-> Fill(tracktime);
    else if( beam_pid==Beam_Pion ) h1 = (TH1F*)gFile-> Get("BLC1_time_pi"), h1-> Fill(tracktime);
    if( -50<tracktime && tracktime<200 ) ntrackBLC1_1++;
    if( -30<tracktime && tracktime<100 ) ntrackBLC1_2++;
    if( -10<tracktime && tracktime<10 ){
      ntrackBLC1++;
      trackBLC1 = bltrackMan->trackBLC1(i);
    }
  }

  if( beam_pid==Beam_Kaon ){
    h1 = (TH1F*)gFile-> Get("ntrackBLC1"), h1-> Fill(bltrackMan->ntrackBLC1());
    h1 = (TH1F*)gFile-> Get("ntrackBLC1_1"), h1-> Fill(ntrackBLC1_1);
    h1 = (TH1F*)gFile-> Get("ntrackBLC1_2"), h1-> Fill(ntrackBLC1_2);
  }
  else if( beam_pid==Beam_Pion ){
    h1 = (TH1F*)gFile-> Get("ntrackBLC1_pi"), h1-> Fill(bltrackMan->ntrackBLC1());
    h1 = (TH1F*)gFile-> Get("ntrackBLC1_1_pi"), h1-> Fill(ntrackBLC1_1);
    h1 = (TH1F*)gFile-> Get("ntrackBLC1_2_pi"), h1-> Fill(ntrackBLC1_2);
  }

  if( ntrackBLC1!=1 ){
    Clear();
    return true;
  }
#if BLDC_WIDE_TIMING
  if( ntrackBLC1_2!=1 ){
    Clear();
    h1 = (TH1F*)rtFile-> Get("check_hist"), h1-> Fill(2);
    //    std::cout<<" !!! BLC1 timeing cut !!!"<<std::endl;
    return true;
  }
#endif
  if( beam_pid==Beam_Kaon ) h1 = (TH1F*)gFile-> Get("BLC1_chi2"), h1-> Fill(trackBLC1->chi2all());
  else if( beam_pid==Beam_Pion ) h1 = (TH1F*)gFile-> Get("BLC1_chi2_pi"), h1-> Fill(trackBLC1->chi2all());
#if BLDC_CHI_CUT
  if( trackBLC1->chi2all()>BLDC_CHI2_MAX ){
    Clear();
    //    std::cout<<" !!! BLC1 chi2 cut !!!"<<std::endl;
    return true;
  }
#endif

  if( beam_pid==Beam_Kaon){
    h1_ER-> Fill(3);
    if( header->IsTrig(Trig_Kf)      ) h1_Kf-> Fill(3);
    if( header->IsTrig(Trig_Neutral) ) h1_N-> Fill(3);
  }

  int ntrackBLC2=0;
  int ntrackBLC2_1=0;
  int ntrackBLC2_2=0;
  LocalTrack *trackBLC2=0;
  for( int i=0; i<bltrackMan->ntrackBLC2(); i++ ){
    double tracktime = bltrackMan->trackBLC2(i)->GetTrackTime();
    if( beam_pid==Beam_Kaon ) h1 = (TH1F*)gFile-> Get("BLC2_time"), h1-> Fill(tracktime);
    else if( beam_pid==Beam_Pion ) h1 = (TH1F*)gFile-> Get("BLC2_time_pi"), h1-> Fill(tracktime);
    if( -50<tracktime && tracktime<200 ) ntrackBLC2_1++;
    if( -30<tracktime && tracktime<100 ) ntrackBLC2_2++;
    if( -10<tracktime && tracktime<10 ){
      ntrackBLC2++;
      trackBLC2 = bltrackMan->trackBLC2(i);
    }
  }
  if( beam_pid==Beam_Kaon ){
    h1 = (TH1F*)gFile-> Get("ntrackBLC2"), h1-> Fill(bltrackMan->ntrackBLC2());
    h1 = (TH1F*)gFile-> Get("ntrackBLC2_1"), h1-> Fill(ntrackBLC2_1);
    h1 = (TH1F*)gFile-> Get("ntrackBLC2_2"), h1-> Fill(ntrackBLC2_2);
  }
  else if( beam_pid==Beam_Kaon ){
    h1 = (TH1F*)gFile-> Get("ntrackBLC2_pi"), h1-> Fill(bltrackMan->ntrackBLC2());
    h1 = (TH1F*)gFile-> Get("ntrackBLC2_1_pi"), h1-> Fill(ntrackBLC2_1);
    h1 = (TH1F*)gFile-> Get("ntrackBLC2_2_pi"), h1-> Fill(ntrackBLC2_2);
  }

  if( ntrackBLC2!=1 ){
    Clear();
    return true;
  }
#if BLDC_WIDE_TIMING
  if( ntrackBLC2_2!=1 ){
    //    std::cout<<" !!! BLC2 timeingcut !!!"<<std::endl;
    Clear();
    h1 = (TH1F*)rtFile-> Get("check_hist"), h1-> Fill(3);
    return true;
  }
#endif
  if( beam_pid==Beam_Kaon ) h1 = (TH1F*)gFile-> Get("BLC2_chi2"), h1-> Fill(trackBLC2->chi2all());
  else if( beam_pid==Beam_Pion ) h1 = (TH1F*)gFile-> Get("BLC2_chi2_pi"), h1-> Fill(trackBLC2->chi2all());

#if BLDC_CHI_CUT
  if( trackBLC2->chi2all()>BLDC_CHI2_MAX ){
    Clear();    return true;
  }
#endif

  if( beam_pid==Beam_Kaon){
    h1_ER-> Fill(4);
    if( header->IsTrig(Trig_Kf)      ) h1_Kf-> Fill(4);
    if( header->IsTrig(Trig_Neutral) ) h1_N-> Fill(4);
  }

  beamSpec->TMinuitFit(trackBLC1, trackBLC2, confMan);
  if( beam_pid==Beam_Kaon ) h1 = (TH1F*)gFile-> Get("D5_chi2"), h1-> Fill(beamSpec->chisquare());
  else if( beam_pid==Beam_Pion ) h1 = (TH1F*)gFile-> Get("D5_chi2"), h1-> Fill(beamSpec->chisquare());
#if BLDC_CHI_CUT
  if( beamSpec->chisquare()>20 ){
    //    std::cout<<" !!! BLC1 chi2 cut !!!"<<std::endl;
    Clear();
    return true;
  }
#endif
  double D5mom = beamSpec-> mom();
  if( nBHD==1 ){
    h1 = (TH1F*)gFile-> Get(Form("D5_mom_BHD%d", BHDseg)), h1-> Fill(D5mom);
  }
  if( beam_pid==Beam_Kaon )  h1 = (TH1F*)gFile-> Get("D5_mom"), h1->Fill(beamSpec->mom());
  else if( beam_pid==Beam_Pion ) h1= (TH1F*)gFile-> Get("D5_mom_pi"), h1-> Fill(beamSpec->mom());

  bool BHD_flag = false;
  for( int i=0; i<blMan->nBHD(); i++ ){
    if( blMan->BHD(i)->CheckRange() ){
      if( 5<blMan->BHD(i)->seg() && blMan-> BHD(i)->seg()<16 ){
	int seg = blMan->BHD(i)->seg();
	if( BHD_MATCH_MIN[seg-6]<D5mom && D5mom<BHD_MATCH_MAX[seg-6] ) BHD_flag=true;
      }
    }
  }  
  if( !BHD_flag ){
    //    std::cout<<" !!! not match beam mom !!!"<<std::endl;
    return true;
  }
  if( beam_pid==Beam_Kaon){
    h1_ER-> Fill(5);
    if( header->IsTrig(Trig_Kf)      ) h1_Kf-> Fill(5);
    if( header->IsTrig(Trig_Neutral) ) h1_N-> Fill(5);  
  }

  int ntrackBPC=0;
  int ntrackBPC_1=0;
  int ntrackBPC_2=0;
  LocalTrack *trackBPC=0;
  for( int i=0; i<bltrackMan->ntrackBPC(); i++ ){
    double tracktime = bltrackMan->trackBPC(i)->GetTrackTime();
    if( beam_pid==Beam_Kaon ) h1 = (TH1F*)gFile-> Get("BPC_time"), h1-> Fill(tracktime);

    if( -50<tracktime && tracktime<200 ) ntrackBPC_1++;
    if( -30<tracktime && tracktime<100 ) ntrackBPC_2++;
    if( -10<tracktime && tracktime<10 ){
      ntrackBPC++;
      trackBPC = bltrackMan->trackBPC(i);
    }
  }
  if( beam_pid==Beam_Kaon ){
    h1 = (TH1F*)gFile-> Get("ntrackBPC"), h1-> Fill(bltrackMan->ntrackBPC());
    h1 = (TH1F*)gFile-> Get("ntrackBPC_1"), h1-> Fill(ntrackBPC_1);
    h1 = (TH1F*)gFile-> Get("ntrackBPC_2"), h1-> Fill(ntrackBPC_2);
  }
  else if( beam_pid==Beam_Pion ){
    h1 = (TH1F*)gFile-> Get("ntrackBPC_pi"), h1-> Fill(bltrackMan->ntrackBPC());
    h1 = (TH1F*)gFile-> Get("ntrackBPC_1_pi"), h1-> Fill(ntrackBPC_1);
    h1 = (TH1F*)gFile-> Get("ntrackBPC_2_pi"), h1-> Fill(ntrackBPC_2);
  }

  if( ntrackBPC!=1 ){
    Clear();
    return true;
  }
#if BLDC_WIDE_TIMING
  if( ntrackBPC_2!=1 ){
    //    std::cout<<" !!! PC timeing cut !!!"<<std::endl;
    Clear();
    h1 = (TH1F*)rtFile-> Get("check_hist"), h1-> Fill(4);
    return true;
  }
#endif
  if( beam_pid==Beam_Kaon ) h1 = (TH1F*)gFile-> Get("BPC_chi2"), h1-> Fill(trackBPC->chi2all());
  else if( beam_pid==Beam_Pion ) h1 = (TH1F*)gFile-> Get("BPC_chi2_pi"), h1-> Fill(trackBPC->chi2all());

#if BLDC_CHI_CUT
  if( trackBPC->chi2all()>BLDC_CHI2_MAX ){
    std::cout<<" !!! BPC chi2 cut !!!"<<std::endl;
    Clear();
    return true;
  }
#endif
  double center_z = 0.5*(-130-20.3);
  TVector3 BLC2_track_pos = trackBLC2-> GetPosatZ(center_z);
  TVector3 BPC_track_pos  = trackBPC-> GetPosatZ(center_z);
  TVector3 BLC2_mom_dir = trackBLC2-> GetMomDir();
  TVector3 BPC_mom_dir = trackBPC-> GetMomDir();
  if( beam_pid==Beam_Kaon ){
    h2 = (TH2F*)rtFile-> Get("BLC2BPC_diff"), h2-> Fill((BLC2_track_pos-BPC_track_pos).X(), (BLC2_track_pos-BPC_track_pos).Y());
    h2 = (TH2F*)rtFile-> Get("BLC2BPC_dir_diff"), h2-> Fill((BLC2_mom_dir-BPC_mom_dir).X(), (BLC2_mom_dir-BPC_mom_dir).Y());
  } else if( beam_pid==Beam_Pion ){
    h2 = (TH2F*)rtFile-> Get("BLC2BPC_diff_pi"), h2-> Fill((BLC2_track_pos-BPC_track_pos).X(), (BLC2_track_pos-BPC_track_pos).Y());
    h2 = (TH2F*)rtFile-> Get("BLC2BPC_dir_diff_pi"), h2-> Fill((BLC2_mom_dir-BPC_mom_dir).X(), (BLC2_mom_dir-BPC_mom_dir).Y());
  }

  TVector3 BLC2BPC_diff = BLC2_track_pos-BPC_track_pos;
  TVector3 BLC2BPC_dir_diff = BLC2_mom_dir-BPC_mom_dir;
  if( BLC2BPC_diff.X()<BLC2BPC_X_MIN || BLC2BPC_X_MAX<BLC2BPC_diff.X() || BLC2BPC_diff.Y()<BLC2BPC_Y_MIN || BLC2BPC_Y_MAX<BLC2BPC_diff.Y() ||
      BLC2BPC_dir_diff.X()<BLC2BPC_dX_MIN || BLC2BPC_dX_MAX<BLC2BPC_dir_diff.X() || BLC2BPC_dir_diff.Y()<BLC2BPC_dY_MIN || BLC2BPC_dY_MAX<BLC2BPC_dir_diff.Y() ){
    std::cout<<" !!! BLC2 BPC not match !!!"<<std::endl;
    Clear();
    return true;
  }

  if( beam_pid==Beam_Kaon){
    h1_ER-> Fill(6);
    if( header->IsTrig(Trig_Kf)      ) h1_Kf-> Fill(6);
    if( header->IsTrig(Trig_Neutral) ) h1_N-> Fill(6);
  }

  for( int i=0; i<bltrackMan->ntrackBLC1a(); i++ ){
    //    std::cout<<" Fill BLC1a : "<<i<<std::endl;
    fillBLC1a(bltrackMan->trackBLC1a(i));
  }
  for( int i=0; i<bltrackMan->ntrackBLC1b(); i++ ){
    //    std::cout<<" Fill BLC1b : "<<i<<std::endl;
    fillBLC1b(bltrackMan->trackBLC1b(i));
  }
  for( int i=0; i<bltrackMan->ntrackBLC2a(); i++ ){
    //    std::cout<<" Fill BLC2a : "<<i<<std::endl;
    fillBLC2a(bltrackMan->trackBLC2a(i));
  }
  for( int i=0; i<bltrackMan->ntrackBLC2b(); i++ ){
    //    std::cout<<" Fill BLC2b : "<<i<<std::endl;
    fillBLC2b(bltrackMan->trackBLC2b(i));
  }

  //  std::cout<<" Fill BPC : "<<std::endl;
  fillBPC(trackBPC);

  TVector3 BPC_FF = trackBPC-> GetPosatZ(-4);
  if( beam_pid==Beam_Kaon ) h1= (TH1F*)gFile->Get("BeamProf"), h1->Fill(BPC_FF.X(), BPC_FF.Y());
  else if( beam_pid==Beam_Pion ) h1= (TH1F*)gFile->Get("BeamProf_pi"), h1->Fill(BPC_FF.X(), BPC_FF.Y());

  if( header->IsTrig(Trig_Kf) ){
    h1= (TH1F*)gFile->Get("BeamProf_Kf"), h1->Fill(BPC_FF.X(), BPC_FF.Y());
    if( beam_pid==Beam_Kaon){
      if( GeomTools::GetID(BPC_FF)==CID_Fiducial ) h1_Kf-> Fill(7);
    }
  }

  //****************************************//
  //*** CDC Tracking File Event Matching ***//
  //****************************************//
  //  std::cout<<"CDC File Matching START"<<std::endl;
  cdcTree-> GetEntry(CDC_Event_Num);
  if( CDC_Event_Num>cdcTree->GetEntries() ){
    Clear();
    std::cout<<"  CDC File Reach End"<<std::endl;
    return true;
  }
  while( cdcHeader->ev()<Event_Number ){
    CDC_Event_Num++;
    if( CDC_Event_Num>cdcTree->GetEntries() ){
      //      std::cout<<" Event Number > CDCTree's Entries"<<std::endl;
      Clear();
      return true;
    }
    cdcTree-> GetEntry(CDC_Event_Num);
  }
  if( cdcHeader->ev()>Event_Number ){
    //    std::cout<<" CDCHeader > Event Number"<<std::endl;
    Clear();
    return true;
  }
  if( cdcHeader->ev()!=Event_Number ){
    //    std::cout<<"  !!! CDC File Event matching miss !!!"<<std::endl;
    Clear();
    return false;
  }
  std::cout<<" CDC > 1track "<<std::endl;

  histMan-> setT0time(T0time);
  histMan-> setD5mom(D5mom);
  histMan-> setTrackBPC(trackBPC);
  histMan-> ana(header);
  histMan-> fill(header);

  Clear();
  return true;
}

void EventAnalysisNpipi::Clear()
{
  header->Clear();
  histMan-> Clear();
  beamSpec-> Clear();
}

void EventAnalysisNpipi::InitializeHistogram()
{
  rtFile-> cd();
  new TH1F("EventReduction", "Event Reduction", 21, -1.5, 19.5);
  new TH1F("Kf_Reduction", "K/f Event Reduction", 10, -0.5, 9.5);
  new TH1F("N_Reduction", "Neutral Event Reduction", 15, -0.5, 14.5);

  new TH1F("check_hist", "histogram for check", 100, -0.5, 99.5);

  new TH1F("nT0", "nT0", 6, -0.5, 5.5);

  new TH1F("nBHD", "nBHD", 21, -0.5, 21.5);
  new TH1F("BHDT0",    "BHD-T0", 5000, 0.0, 50);
  new TH1F("BHDT0_Kf", "BHD-T0", 5000, 0.0, 50);
  new TH1F("BHDT0_K",  "BHD-T0", 5000, 0.0, 50);
  new TH1F("BHDT0_pi", "BHD-T0", 5000, 0.0, 50);

  new TH1F("ntrackBLC1", "BLC1 n track", 10, -0.5, 9.5);
  new TH1F("ntrackBLC1_1", "BLC1 n track", 10, -0.5, 9.5);
  new TH1F("ntrackBLC1_2", "BLC1 n track", 10, -0.5, 9.5);
  new TH1F("BLC1_time", "BLC1 track time", 10000, -200, 800);
  new TH1F("BLC1_chi2", "BLC1 chi-square", 1000, 0.0, 100);
  new TH1F("ntrackBLC1_pi", "BLC1 n track", 10, -0.5, 9.5);
  new TH1F("ntrackBLC1_1_pi", "BLC1 n track", 10, -0.5, 9.5);
  new TH1F("ntrackBLC1_2_pi", "BLC1 n track", 10, -0.5, 9.5);
  new TH1F("BLC1_time_pi", "BLC1 track time", 10000, -200, 800);
  new TH1F("BLC1_chi2_pi", "BLC1 chi-square", 1000, 0.0, 100);

  new TH1F("ntrackBLC2", "BLC2 n track", 10, -0.5, 9.5);
  new TH1F("ntrackBLC2_1", "BLC2 n track", 10, -0.5, 9.5);
  new TH1F("ntrackBLC2_2", "BLC2 n track", 10, -0.5, 9.5);
  new TH1F("BLC2_time", "BLC2 track time", 10000, -200, 800);
  new TH1F("BLC2_chi2", "BLC2 chi-square", 1000, 0.0, 100);
  new TH1F("ntrackBLC2_pi", "BLC2 n track", 10, -0.5, 9.5);
  new TH1F("ntrackBLC2_1_pi", "BLC2 n track", 10, -0.5, 9.5);
  new TH1F("ntrackBLC2_2_pi", "BLC2 n track", 10, -0.5, 9.5);
  new TH1F("BLC2_time_pi", "BLC2 track time", 10000, -200, 800);
  new TH1F("BLC2_chi2_pi", "BLC2 chi-square", 1000, 0.0, 100);

  new TH1F("ntrackBPC", "BPC n track", 10, -0.5, 9.5);
  new TH1F("ntrackBPC_1", "BPC n track", 10, -0.5, 9.5);
  new TH1F("ntrackBPC_2", "BPC n track", 10, -0.5, 9.5);
  new TH1F("BPC_time", "BPC track time", 10000, -200, 800);
  new TH1F("BPC_chi2", "BPC chi-square", 1000, 0.0, 100);
  new TH1F("ntrackBPC_pi", "BPC n track", 10, -0.5, 9.5);
  new TH1F("ntrackBPC_1_pi", "BPC n track", 10, -0.5, 9.5);
  new TH1F("ntrackBPC_2_pi", "BPC n track", 10, -0.5, 9.5);
  new TH1F("BPC_time_pi", "BPC track time", 10000, -200, 800);
  new TH1F("BPC_chi2_pi", "BPC chi-square", 1000, 0.0, 100);

  new TH2F("BLC2BPC_diff", "BLC2 BPC pos diff", 1000, -10, 10, 1000, -10, 10);
  new TH2F("BLC2BPC_dir_diff", "BLC2 BPC dir diff", 1000, -0.1, 0.1, 1000, -0.1, 0.1);
  new TH2F("BLC2BPC_diff_pi", "BLC2 BPC pos diff", 1000, -10, 10, 1000, -10, 10);
  new TH2F("BLC2BPC_dir_diff_pi", "BLC2 BPC dir diff", 1000, -0.1, 0.1, 1000, -0.1, 0.1);

  new TH1F("D5_mom",  "Beam Momentum by D5", 1000,  0.5, 1.5);
  for( int seg=1; seg<=20; seg++ )  new TH1F(Form("D5_mom_BHD%d", seg),  "Beam Momentum by D5", 1000,  0.5, 1.5);
  new TH1F("D5_chi2", "D5 chisquare", 1000,  0.0, 100);
  new TH1F("D5_mom_pi",  "Beam Momentum by D5", 1000,  0.5, 1.5);
  new TH1F("D5_chi2_pi", "D5 chisquare", 1000,  0.0, 100);

  new TH2F("BeamProf",    "Beam Profile at FF",        1000, -30, 30, 1000, -30, 30);
  new TH2F("BeamProf_pi",    "Beam Profile at FF",        1000, -30, 30, 1000, -30, 30);
  new TH2F("BeamProf_Kf", "Beam Profile at FF w/ K/f", 1000, -30, 30, 1000, -30, 30);

  initHistBLDC();
  initHistT0NC();
}

void EventAnalysisNpipi::Finalize()
{
  std::cout<<"EventAnalysisNpipi Finish   Event Number : "<<Event_Number<<"  Event Number on Spill : "<<Event_Number_onSpill<<std::endl;
  delete cdcFile;
  delete blMan;
  delete cdsMan;
  delete header;

  delete bltrackMan;
  delete beamSpec;

  rtFile-> cd();
  histMan->finit();
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisNpipi *event = new EventAnalysisNpipi();
  return (EventTemp*)event;
}
